# FAMRDP-JAX 全局一维化求解器项目路线图

> 日期：2026-04-18  
> 目标：面向复杂算例计算与神经网络辅助学习的全局 1D JAX CFD 求解器

---

## 项目总体目标

```
可微分 CFD 求解器（全局 1D，单 GPU JIT）
         │
         ├── 复杂算例计算（>100 块，生产级精度）
         └── 神经网络辅助
               ├── 数据驱动湍流/转捩模型
               ├── 神经网络加速代理模型
               └── 气动反设计 / 拓扑优化
```

**现状盘点**

| 模块 | `famrdp_jax/`（3D 多块） | `1d_famrdp/`（1D 原型） |
|---|---|---|
| 网格读取 + 拓扑 | ✅ 完成，已验证 | ✅ 完成（2D） |
| 索引前处理 | ✗ 无需 | ✅ MUSCL-2，2D |
| 无粘 RHS（MUSCL+Roe） | ✅ 3D，已验证 | ✅ 2D，无 Harten 修复 |
| 粘性 RHS（Sutherland） | ✅ 3D，已验证 | ✗ 未完成 |
| TVD-RK3 | ✅ 完成 | ✅ 完成 |
| 边界条件 | ✅ 4 种，已验证 | ⚠️ 2D 部分 |
| 高阶格式（dcsh5pi） | ✗ 未开始 | ✗ 未开始 |
| 可微性（jax.grad） | ✅ 已验证 | ✗ 未测试 |
| 神经网络接口 | ✗ 未开始 | ✗ 未开始 |

---

## 阶段划分

```
Phase 0  ──► Phase 1  ──► Phase 2  ──► Phase 3  ──► Phase 4
基础迁移     3D 完整化    高阶格式     复杂算例     NN 集成
（已有→1D）  （验证对拍）  （dcsh5pi）  （生产级）   （可微优化）
 ~2 周        ~3 周        ~4 周        ~2 周        持续
```

---

## Phase 0：基础迁移（将 famrdp_jax 的验证结果固化为 1D 框架的基础）

**目标**：在 1D 框架内复现当前 3D 多块法已通过的所有验证测试。  
**策略**：直接复用 `1d_famrdp/preset/domain_manager.py` 和 `scheme_factory.py` 的核心索引逻辑，补全 3D + 粘性 + 正确边界条件。

### M0.1 — 3D 索引前处理完整化

`1d_famrdp` 的 `scheme_factory.py` 使用 `conv_general_dilated` 生成索引，当前支持 2D（`dim=2`）。需要验证其 3D 路径（`nck > 1`）的正确性。

- [ ] 用 `test3.grd`（2 块，51×51×11）运行 `DomainManager`，检查 `block_cells_index` shape 和 ghost 映射
- [ ] 用 `SchemeFactory` 生成 MUSCL-2 的 `inter_index`，与 `famrdp_jax` 的切片等价性验证（取同一块内部面，数值应相同）
- [ ] 验证 CUT1TO1 接口的 ghost 索引映射与 `halo.py` 结果一致

### M0.2 — 1D RHS 核（无粘，对拍 Fortran）

- [ ] 将 `famrdp_jax/flux/roe.py`（含 Harten 熵修复）迁移为 1D 接口：`roe_flux_1d(UL, UR, n_face)` 其中 shape 为 `(N_faces, 5)` 和 `(N_faces, 3)`
- [ ] 将 `muscl2_reconstruct` 改为基于 `inter_index` 的 gather 版本
- [ ] 组装 `compute_invscd_rhs_1d(prime, inter_index, face_metric, vol)`
- [ ] 验收：`test3.grd` 单步 RHS，与 `famrdp_jax` 结果 diff ≤ 1e-10

### M0.3 — 1D 粘性 RHS

- [ ] Sutherland 粘性律迁移（已在 `famrdp_jax/physics/eos.py`，接口不变）
- [ ] 中心差粘性通量改写为 `difference_index` gather 形式
- [ ] 验收：静止等温场零粘性 RHS

### M0.4 — 1D 边界条件

- [ ] 无滑移壁面、滑移壁面、对称面、远场：将 `bc_type` 的全局面索引（`start_index`/`end_index`）直接用于 `prime.at[boundary_idx].set(...)` 操作
- [ ] 验收：与 `famrdp_jax/physics/bc/` 单步对拍 ≤ 1e-13

### M0.5 — 端到端 10 步对拍

- [ ] 将 `step_1d(buffer, precomp, dt)` 接入 TVD-RK3
- [ ] 使用 `@jax.jit` 修饰，确认首步编译通过
- [ ] 验收：10 步整场 diff ≤ 1e-10（复用 `validation/references/test3/`）

### M0.6 — 可微性验证

- [ ] `jax.grad(loss)(prime)` 返回非零有限值，`loss = sum(step_1d(...))`
- [ ] 确认无 gather 反向的 NaN 或无穷梯度（用 `jnp.where` 守护 Roe 奇点）

**Phase 0 里程碑验收**：所有 56 个测试（目前 3D 框架的）在 1D 框架下等价通过，JIT 编译成功，`jax.grad` 可用。

---

## Phase 1：3D 生产级完整化

**目标**：1D 框架支持真实复杂算例所需的完整物理和数值能力。

### M1.1 — 完整 BC 体系

- [ ] 特征线入流-出流（Riemann 不变量）精确实现，替代当前 `farfield` 的简化版（ghost = 自由流）
- [ ] 周期边界（旋转对称算例需要）
- [ ] 固壁绝热 / 等温统一接口

### M1.2 — CFL 自适应时间步

- [ ] 逐单元计算局部 CFL 限制的 dt：`dt = CFL * min(vol / (|u|+c) * J)`
- [ ] 全局 `jnp.min` 归约，保证 JIT 内完成（无 Python 条件）

### M1.3 — 残差监控与 I/O

- [ ] 每步 L∞ 残差归约（`jax.debug.callback` 输出，不破坏 JIT 图）
- [ ] Tecplot ASCII 输出（调用已有 `io_output.py` 逻辑，在 JIT 外执行）
- [ ] 检查点 pickle 存储/恢复（`jnp.save` / `jnp.load`）

### M1.4 — 基准测试

- [ ] `test3.grd`（2 块，小算例）：GPU 100 步中位时间
- [ ] 中等算例（~50 块）：与单核 Fortran 吞吐比
- [ ] 验收：GPU vs CPU-Fortran 比值记录在 `benchmarks/` 下

**Phase 1 里程碑验收**：完整运行圆柱绕流算例（test3）收敛至残差 < 1e-6，I/O 完整，benchmark 数据完备。

---

## Phase 2：高阶格式（dcsh5pi + Roe-scmp）

**目标**：实现阶段 2 的核心数值格式，这是高精度复杂算例的基础。

### M2.1 — dcsh5pi 索引前处理

`dcsh5pi` 是 5 阶紧致差分格式，模板宽度为 5。在 1D 框架下，前处理只需将 `cell_range` 扩展至 `[-3,-2,-1,0,1,2]`，生成更宽的 `inter_index`。

- [ ] `SchemeFactory` 中增加 `scheme_type='dcsh5pi'` 分支，生成 6 行 `inter_index`
- [ ] 实现 dcsh5pi 插值系数（参考 Fortran `rhs_invscd.f90` 的 `dcsh5pi` 部分）
- [ ] 块界面退化：近界面（距界面 < 2 格）降阶至 MUSCL-2，索引前处理自动标记退化面
- [ ] 验收：光滑流场上 dcsh5pi 比 MUSCL-2 截断误差低 2 阶

### M2.2 — Roe-scmp（低马数预处理 Roe）

- [ ] 实现 Roe-scmp 通量（参考 Fortran `rhs_invscd.f90` 的 scmp 分支），1D 接口 `roe_scmp_flux_1d`
- [ ] 验收：低马数（Ma < 0.1）算例与 MUSCL+Roe 比较，压力残差下降更快

### M2.3 — 高阶格式端到端验证

- [ ] 生成 dcsh5pi + Roe-scmp 的 Fortran reference（修改 `gen_refs.f90` 调用高阶路径）
- [ ] 10 步整场 diff ≤ 1e-10

**Phase 2 里程碑验收**：dcsh5pi + Roe-scmp 配置跑通算例，精度优于 MUSCL-2，diff 测试通过。

---

## Phase 3：复杂算例生产计算

**目标**：验证求解器在工程级复杂算例（翼型、机翼、旋翼等）上的可靠性。

### M3.1 — 大规模算例适配

- [ ] 测试 > 50 块算例：确认前处理时间可接受（< 5 分钟）
- [ ] 显存估算工具：根据块数 + ghost 层自动估算 GPU 显存需求
- [ ] 若单 GPU 显存不足：提供按块分批计算的 fallback（牺牲 JIT 效率）

### M3.2 — 湍流模型接口（可选占位）

> 范围约束 CR-01 说明不实现湍流模型，但为后续 NN 替代湍流模型留接口

- [ ] 定义 `turbulence_closure(prime, metrics) -> mu_t`：接受流场，返回湍流粘性
- [ ] 粘性 RHS 中加入 `mu + mu_t` 路径（mu_t=0 时退化为层流）
- [ ] 这个接口将来可直接插入神经网络模型

### M3.3 — 收敛性测试

- [ ] NACA0012 翼型（亚声速层流），与文献 Cl/Cd 对比
- [ ] 圆柱绕流 Re=200（层流），与实验 St 数对比

**Phase 3 里程碑验收**：2 个工程算例结果与文献/实验一致，求解器可用于生产。

---

## Phase 4：神经网络集成

**目标**：利用可微分求解器实现三类 NN 应用。

### 应用路线 A — 神经网络湍流/转捩模型

**原理**：用 NN 替代或修正 RANS 湍流模型（SA/Menter）中的代数关系。

```
训练数据：DNS/LES 高精度流场数据
输入特征：局部流场特征（q, grad_q, wall_dist）
输出：修正后的 mu_t 或湍流源项
损失函数：L = ||CFD(q; NN) - q_ref||² 
可微路径：jax.grad(L)(NN_params)，穿透 CFD 的 step 函数
```

- [ ] M4.A1：设计特征提取器（局部 stencil → 特征向量），接口与 Phase 3 的 `turbulence_closure` 对接
- [ ] M4.A2：小型 MLP（Equinox 或 Flax NNX）作为闭合关系，参数量 < 10K
- [ ] M4.A3：从高精度参考数据中监督预训练，再做端到端微调
- [ ] M4.A4：验收：圆柱绕流 Re=3900（湍流），Cl 与 LES 误差 < 5%

### 应用路线 B — 代理模型加速

**原理**：用 NN 拟合 CFD 的输入-输出映射，代替耗时的完整求解。

```
训练数据：多组（参数, 收敛流场）对
输入：几何参数 / 来流条件
输出：表面压力分布 / 气动力系数
损失函数：L = ||NN(params) - CFD(params)||²
用途：设计空间探索、不确定性量化
```

- [ ] M4.B1：设计数据采集流程（拉丁超立方采样 → 批量 CFD 运算）
- [ ] M4.B2：几何/来流参数编码器 + 流场解码器（Graph NN 或 FNO）
- [ ] M4.B3：验收：NACA 系列翼型 Cl 预测，测试集误差 < 2%

### 应用路线 C — 气动反设计（最核心的可微分应用）

**原理**：将几何参数化后，通过 `jax.grad` 直接对 CFD 求解器反向传播求梯度，梯度下降优化几何。

```python
@jax.jit
def loss(shape_params):
    mesh  = parameterize_geometry(shape_params)    # 参数化网格变形
    state = run_cfd(mesh, n_steps=1000)            # 可微 CFD 前向推进
    Cl    = compute_lift(state)                     # 升力系数
    Cd    = compute_drag(state)                     # 阻力系数
    return -Cl / Cd                                 # 最大化升阻比

# 梯度由 JAX 自动计算，穿透整个 CFD 计算图
grad_params = jax.grad(loss)(initial_params)
optimized   = optax.adam(lr=1e-3).update(grad_params, ...)
```

- [ ] M4.C1：Hicks-Henne 或 CST 翼型参数化 → 网格变形（可微实现）
- [ ] M4.C2：升阻力系数计算（表面积分，可微）
- [ ] M4.C3：端到端梯度验证（有限差分 vs `jax.grad`，相对误差 < 1%）
- [ ] M4.C4：单翼型反设计算例（目标 Cl=0.8，优化 10 个 CST 参数）

**Phase 4 里程碑验收**：至少完成路线 A 或路线 C 的一个完整算例，展示 NN 与 CFD 可微耦合。

---

## 关键技术决策

### 决策 1：`famrdp_jax/` 与 `1d_famrdp/` 的关系

**建议**：以 `famrdp_jax/` 的验证结果（physics 模块、参考数据）为基准，以 `1d_famrdp/` 的索引框架（DomainManager、SchemeFactory）为骨架，**在 `famrdp_jax/` 中新建 `famrdp_jax/core1d/` 子包**，而非直接在 `1d_famrdp/` 上扩展。

理由：`1d_famrdp` 是 2D 原型，其物理模块（Roe、BC）未经严格验证；`famrdp_jax` 的 physics 已通过 Fortran 对拍。

### 决策 2：与 Flax/Equinox 的集成

Phase 4 的 NN 组件建议使用 **Equinox**（`eqx.nn`），原因：
- 纯 pytree 结构，与 JAX `jit`/`grad`/`vmap` 无缝集成
- 不依赖 `flax.struct`（当前环境已知无法使用）
- 模型参数和 CFD 状态统一在同一 `jax.tree_util` 框架下可微

### 决策 3：梯度稳定性

长时间积分（>100 步）的梯度会指数爆炸。Phase 4 气动反设计建议：
- 使用 **截断时间反向传播（TBPTT）**：每 K 步截断梯度，分段优化
- 或使用 **伴随方法（adjoint）**：解伴随方程代替完整反向传播，显存 O(1)

JAX 支持通过 `jax.lax.stop_gradient` 和自定义 `jvp/vjp` 实现伴随。

---

## 时间估算与优先级

| 阶段 | 预计周期 | 依赖 | 优先级 |
|---|---|---|---|
| Phase 0：基础迁移 | 2–3 周 | 当前 3D 验证结果 | P0 必须 |
| Phase 1：完整化 | 2–3 周 | Phase 0 | P0 必须 |
| Phase 2：高阶格式 | 3–4 周 | Phase 1 | P1 重要 |
| Phase 3：复杂算例 | 2 周 | Phase 1（高阶可选） | P1 重要 |
| Phase 4A：NN 湍流 | 4–6 周 | Phase 3 + 高精度数据 | P2 |
| Phase 4C：反设计 | 3–4 周 | Phase 1 + 参数化网格 | P1 |
| Phase 4B：代理模型 | 4–6 周 | Phase 3 + 批量数据 | P2 |

**最短关键路径**（直达反设计）：Phase 0 → Phase 1 → Phase 4C  
**最完整路径**（生产 + NN）：Phase 0 → Phase 1 → Phase 2 → Phase 3 → Phase 4

---

## 下一步行动

1. **确认 Phase 0 起点**：运行 `1d_famrdp/main.py` 的 `test3.grd` 算例，记录当前 2D 路径的实际输出，确认 DomainManager 的 3D 分支是否已可用
2. **制定 Phase 0 详细实施计划**（参照 `docs/superpowers/plans/` 的格式）
3. **确定 NN 应用优先路线**（A/B/C 中哪个最先）

---

*文档版本 v1.0 — 初稿，待用户确认后转化为各阶段实施计划*
