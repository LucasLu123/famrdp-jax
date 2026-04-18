# 多块 3D 方法 vs 全局 1D 方法：设计对比分析

> 日期：2026-04-18  
> 作者：Jiang Yi  
> 对比对象：`famrdp_jax/`（3D 多块方法）vs `1d_famrdp/`（全局 1D 方法）

---

## 1. 背景与动机

FAMRDP 是一个 3D 结构多块有限体积 CFD 求解器。将其移植到 JAX 时，核心挑战是：
**如何在保留多块拓扑的同时，让数值核高效地接入 `jax.jit` 和 GPU 向量化？**

目前存在两种实现路线：

| 路线 | 代码位置 | 核心思想 |
|---|---|---|
| **3D 多块法** | `famrdp_jax/` | 每块保持独立 3D 数组，Python 循环遍历块 |
| **全局 1D 法** | `1d_famrdp/` | 所有块拼接为单一全局 1D 数组，前处理生成索引 |

---

## 2. 数据结构对比

### 2.1 3D 多块法（`famrdp_jax/`）

```python
# 每块独立存储
@dataclass
class Block:
    q:   jnp.ndarray  # (5, ni, nj, nk)  守恒变量，含 ghost 层
    xyz: jnp.ndarray  # (3, ni, nj, nk)  坐标

@dataclass
class Metrics:
    jac:  jnp.ndarray  # (ni, nj, nk)
    kxyz: jnp.ndarray  # (3, 3, ni, nj, nk)
    vol:  jnp.ndarray  # (ni, nj, nk)

@dataclass
class State:
    blocks:  tuple[Block, ...]    # 逐块列表
    metrics: tuple[Metrics, ...]  # 逐块列表
    t: float
    step: int
```

每块形状各异（块间 shape 可异构），块间通过 Python 循环访问。

### 2.2 全局 1D 法（`1d_famrdp/`）

```python
class simulationBuffer(NamedTuple):
    prime:               jnp.ndarray  # (N_total, 5)  原始变量，全局拼接
    qc:                  jnp.ndarray  # (N_total, 5)  守恒变量
    grid_derivative:     jnp.ndarray  # (N_total, 9)  度量系数
    grid_face_derivative: tuple       # 3 × (N_face_i/j/k, 3)
    jacobi:              jnp.ndarray  # (N_total, 1)
    cell_coordinates:    jnp.ndarray  # (N_total, 3)
```

所有块的真实单元在前处理时展平拼接，N_total = Σ(nci × ncj × nck)。

---

## 3. 核心机制对比

### 3.1 Ghost 层管理

**3D 多块法**

Ghost 层嵌入在每块的 q 数组中：
```
q.shape = (5, ni+2g, nj+2g, nk+2g)
```
Halo exchange 在每步计算前通过 Python 循环完成，将相邻块的内部单元直接写入当前块的 ghost 层：

```python
def halo_exchange(blocks, topos, ghost):
    for i, topo in enumerate(topos):
        for face, info in topo.neighbors.items():
            qs[i] = _copy_face(qs[i], blocks[src_i].q, face, src_face, ghost)
```

- 时机：每个 RK3 子步前执行，共 3 次/步
- 代价：Python 循环，每次 halo exchange 产生 JAX 追踪

**全局 1D 法**

Ghost 层在前处理（`domain_manager.py`）中由索引映射决定：
```python
pad_width = ((g, g), (g, g), (g, g))
cells_pad = np.pad(block_cells_index, pad_width, mode='edge')
```
块间接口（CUT1TO1）由 `_map_connection_face_to_ghost` 将邻块真实单元的全局编号直接写入目标块 ghost 区域的编号表。在计算时，ghost 层访问通过索引自动指向正确的邻块单元，**无需显式的 halo exchange 函数**。

### 3.2 索引生成与模板操作

**3D 多块法**

MUSCL-2 重构直接操作 3D 数组切片：
```python
# famrdp_jax/physics/reconstruct/muscl2.py
def muscl2_reconstruct(pv, axis, limiter):
    sl = lambda s, e: slice_along(pv, axis, s, e)
    delta_right = sl(2, -1) - sl(1, -2)   # roll-like slice
    delta_left  = sl(1, -2) - sl(0, -3)
    ...
```
结构规整，可被 XLA 自动识别为连续内存访问，不需要预处理。

**全局 1D 法**

前处理（`scheme_factory.py`，约 930 行）生成所有方向的模板索引表：
```python
inter_index[axis]  # shape: (stencil_len, N_faces_in_direction)
```
运行时插值为：
```python
derive = prime[inter_index[1:]] - prime[inter_index[:-1]]  # (4, N_faces, 5)
prime_L = prime[inter_index[2]] + k1 * derive_right + ...
```
用高级索引（gather）替代切片，天然支持非规整访问模式（多块拼接后面索引不连续）。

### 3.3 JIT 接入方式

**3D 多块法**

JIT 以整个 `step` 函数为边界：
```python
# 当前未显式标注 @jit，可以外部包装
step_jit = jax.jit(step, static_argnums=...)
```
内部 Python 循环（`for b in blocks`）在第一次 JIT 时被展开（unroll）。块数决定展开深度：2 块 → 展开 2 次 RHS 计算。每块为独立子图。

**全局 1D 法**

```python
@partial(jax.jit, static_argnums=(0,))
def do_integration_step(self, prime, grid_derivative, grid_face_derivative,
                        volume, qc, dtau):
    ...
```
所有索引表（inter_index 等）作为 `self` 的属性被 JIT 作为常量闭包捕获。流场变量（prime, qc）是唯一的动态参数，计算图固定，**只需编译一次**。

---

## 4. 优劣势详细分析

### 4.1 三维多块法（`famrdp_jax/`）

**优势**

1. **代码可读性高**  
   数组维度直接对应物理方向：`q[rho, i, j, k]`，与 Fortran 代码一一对应，便于验证和调试。slice 操作语义清晰，无需理解索引映射关系。

2. **前处理开销为零**  
   无需生成索引表，初始化时间 O(1)，适合快速实验和参数扫描场景。

3. **对拍 Fortran 直接**  
   3D 布局与原 Fortran 完全对应，逐点 diff 无需转换，验证容差（1e-13）直接可测。

4. **数组访问连续**  
   同一块内，切片访问内存连续，XLA 可识别为规整 gather，不产生分散访问。

5. **异构块支持自然**  
   每块 shape 可不同，无需任何 padding 或 masking，不浪费显存。

6. **可微性干净**  
   `jax.grad` 可直接穿透，无 gather/scatter 引起的梯度聚合问题；测试已验证（梯度 max ≈ 9.8，非零有限）。

**劣势**

1. **多块时 Python 循环无法并行**  
   `step` 中的 `for b in blocks` 在 JIT 时被展开而非向量化，N 块 = N 个独立子图，编译时间随块数线性增长。

2. **块间通信每步显式执行**  
   每个 RK3 子步前需 `halo_exchange`，产生额外的 Python→JAX 边界。块数多时性能下降明显。

3. **块内规模受限**  
   单块若很小（典型多块算例每块仅数千格）则 GPU 利用率低，无法充分流水线。

4. **无法做跨块向量化**  
   多块 Roe 通量计算无法合并为单次 `vmap`，因为各块形状不同。

5. **JIT 重编译风险**  
   若块数或 shape 在不同算例间变化，JAX 会重新追踪和编译。

---

### 4.2 全局一维法（`1d_famrdp/`）

**优势**

1. **单次全局 JIT，无重编译**  
   所有块拼接为同一数组，JIT 生成单一计算图，规模固定后不再重编译。GPU 利用率可达最高。

2. **真正向量化**  
   `prime[inter_index]` 将多块所有面的通量计算合并为单次 gather + 向量运算，XLA 可高效调度。

3. **Halo 隐式完成**  
   Ghost 层通过索引映射"自动"访问邻块数据，每步无需显式 halo exchange 调用，减少 Python 调度开销。

4. **高阶格式扩展性强**  
   模板宽度（stencil_len）只影响 inter_index 的行数，切换 MUSCL-5 或 WCNS-6 只需更改前处理，运行时代码不变。

5. **内存布局可控**  
   拼接顺序可按访问模式优化（例如按方向排列面索引），利于 L2 缓存命中。

**劣势**

1. **前处理复杂度高，调试难**  
   `scheme_factory.py` 约 930 行、`domain_manager.py` 约 588 行，索引正确性难以直观验证。索引 bug 在运行时通常表现为精度异常而非 crash，定位困难。

2. **gather 操作的梯度聚合**  
   高级索引（gather）反向时产生 scatter-add，若同一位置被多次索引则梯度正确但性能下降（`jnp.ndarray.at[].add` 无法并行写同一位置）。这对可微优化场景不友好。

3. **显存利用不精确**  
   若各块大小差异大，1D 拼接后内存连续但各块在 1D 数组中的位置需偏移量管理，调试输出时需反向 reshape。

4. **非规整访问模式（scattered gather）**  
   多块拼接后，沿计算方向的索引在物理内存中通常不连续（块尾→下一块头有跳跃），GPU 的 memory coalescing 效果受影响，取决于块的排列顺序。

5. **与 Fortran 对拍困难**  
   1D 布局与 Fortran 的块-坐标布局不同，对拍时必须先反 reshape，增加验证复杂度。

6. **初始化时间长**  
   对于大规模算例（>1000 块），前处理的索引生成（Python 循环 + numpy 操作）可能需要数分钟。

---

## 5. 关键技术指标对比

| 指标 | 3D 多块法 | 全局 1D 法 |
|---|---|---|
| **JIT 编译次数** | 1 次（展开 N 块） | 1 次（单图） |
| **JIT 编译时间** | ∝ 块数 × 块规模 | ∝ 总单元数 |
| **GPU 利用率** | 受单块规模限制 | 全局优化，理论上限高 |
| **Halo exchange 开销** | 每步 3 次 Python 调用 | 隐式（0 次调用） |
| **前处理时间** | O(1)（无索引生成） | O(总面数)（数秒至数分钟） |
| **代码行数（核心）** | ~350 行（rhs+step+bc） | ~2500 行（含前处理） |
| **对拍 Fortran 容易度** | 直接（同维度布局） | 需 reshape 转换 |
| **可微性（jax.grad）** | 已验证，干净 | gather 反向有聚合代价 |
| **异构块（不同 shape）** | 原生支持 | 需 padding 或分段处理 |
| **高阶格式扩展** | 需修改 slice 逻辑 | 只修改前处理 stencil 参数 |
| **调试难度** | 低 | 高（索引正确性难直观验证） |

---

## 6. 适用场景建议

### 推荐使用 3D 多块法的场景

- **精度验证阶段**：需要与 Fortran 逐点对拍（容差 1e-13），3D 布局使比较直接。
- **小块数算例**（≤ 10 块）：Python 展开开销可接受，代码可读性优势明显。
- **可微优化研究**：`jax.grad` 穿透干净，适合气动反设计、形状优化场景。
- **快速原型 / 新物理模型开发**：无需理解索引映射，直接修改 slice 操作。
- **教学与文档目的**：代码结构清晰，易于理解 CFD 算法逻辑。

### 推荐使用全局 1D 法的场景

- **大规模生产计算**（> 100 块）：消除多块 Python 循环，GPU 吞吐量最大化。
- **高阶格式（WCNS-5/6、dcsh5pi）**：宽模板在 1D 索引下扩展成本低，3D 切片逻辑复杂度高。
- **纯前向推进（无 grad）**：避免 gather 反向的梯度聚合性能问题。
- **形状固定的批量算例**：前处理一次，重复使用同一套索引表。

---

## 7. 融合路径建议

两种方法并非互斥，可以分阶段融合：

**近期（阶段 1 完成后）**  
维持 3D 多块法作为主线，完成精度验证（M1.0–M1.9）和阶段 2 高阶格式（M2.x）的参考实现。

**中期**  
借鉴 1D 法的索引预生成思路，在 3D 多块框架内优化 halo exchange：将每次运行时的 `_copy_face` 替换为前处理生成的（src_indices, dst_indices）对，用 `q.at[dst].set(q_other[src])` 执行，消除 Python 字典查找开销。

**远期**  
若生产计算需要 > 100 块的规模，以 1D 法为基础构建高性能计算后端，3D 法保留作为验证参考和可微优化前端，两者共享 `validation/references/` 中的参考数据。

---

## 8. 结论

| 维度 | 3D 多块法优先 | 1D 法优先 |
|---|---|---|
| 代码可维护性 | ✅ | ❌ |
| 可微性与优化研究 | ✅ | ⚠️ |
| 精度验证 | ✅ | ⚠️ |
| GPU 峰值吞吐 | ⚠️（小块）| ✅ |
| 高阶格式扩展 | ⚠️ | ✅ |
| 大规模（>100 块）| ❌ | ✅ |

**当前阶段（M1.x 脚手架，2 块测试算例）：3D 多块法是正确选择。**  
全局 1D 法提供了宝贵的工程参考，其索引预生成模式值得在后续性能优化中借鉴。两者的核心物理（MUSCL-2 + Roe + TVD-RK3）完全相同，差异仅在数据布局与访问策略。
