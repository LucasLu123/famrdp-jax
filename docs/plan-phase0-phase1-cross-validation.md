# Phase 0 & Phase 1 细化计划：双路并行 + 对比验证

> 日期：2026-04-18  
> 作者：Jiang Yi  
> 依据：`docs/roadmap-1d-jax-solver.md` 阶段 0–1 扩写  
> 测试算例：`test3.grd`（圆柱绕流，2 块，51×51×11，Ma=0.1，Re=200）

---

## 0. 总体策略

```
famrdp_jax（3D 多块法）          1d_famrdp（全局 1D 法）
        │                                  │
        │  已验证：Roe+MUSCL-2+粘性         │  已验证：Steger-Warming+MUSCL-2
        │  56 个 pytest 全通过              │  2D/3D 均可运行，无粘性
        │                                  │
        ▼                                  ▼
 Phase 0: 巩固基准               Phase 0: 对齐 3D + 补全物理
 （补 Steger-Warming 路径）       （对拍 Roe、补粘性）
        │                                  │
        └──────── 对比验证层 ──────────────┘
                       │
          两套代码共用同一算例、同一初场
          三级验证准则（见第 1 节）
```

**两套代码各自的定位：**

| 方面 | `famrdp_jax` (3D 多块) | `1d_famrdp` (全局 1D) |
|---|---|---|
| 角色 | **精度基准 + 梯度参考** | **生产效率目标** |
| 当前通量 | Roe + Harten 熵修复 | Steger-Warming |
| 当前粘性 | 完整 Sutherland | 骨架（未实现） |
| 当前 BC | 4 种已验证 | 3 种（远场简化） |
| 测试数量 | 56 个 pytest | 0 个正式测试 |
| 梯度 | jax.grad 已验证 | 未测试 |

---

## 1. 三级验证准则

所有里程碑按以下层级定义验收标准：

| 级别 | 含义 | 典型容差 | 何时使用 |
|---|---|---|---|
| **L1：Fortran 对拍** | 与 gen_refs.f90 生成的 Fortran 参考数据逐点比较 | ≤ 1e-10（无粘 RHS）、≤ 1e-8（粘性 RHS） | 每个物理量独立验证 |
| **L2：双路等价** | 两套 JAX 代码在同一算例上结果一致 | ≤ 1e-8（同通量方案）、稳态量 ≤ 1%（不同通量） | 端到端步数对比 |
| **L3：物理合理性** | 与文献/解析解对比（St 数、Cl、残差阶数） | 圆柱 St ∈ [0.18, 0.22]，Cl 对称收敛 | 长时运行验证 |

**关于通量方案差异的说明：**  
`famrdp_jax` 用 Roe，`1d_famrdp` 用 Steger-Warming，两者是不同的数值格式，单步结果不相等。因此：
- L1 验证各自独立（每套代码对拍自己的 Fortran 参考）
- L2 验证分两种情况：
  - **同通量方案**：在 `1d_famrdp` 中加入 Roe 选项后，与 `famrdp_jax` 逐步对拍（≤ 1e-10）
  - **不同通量方案**：仅比较稳态收敛解（误差 ≤ 1%）和收敛速率

---

## 2. 当前状态清单（起点快照，2026-04-18）

### 2.1 `famrdp_jax` 现状

```
✅ M1.0  xyz 读取（PLOT3D binary，order='F'）
✅ M1.1  Metrics（4 阶中心差，Jacobian，vol=|jac|）
✅ M1.2  Halo exchange（CUT1TO1，ghost=1，orient=0）
✅ M1.3a 无滑移壁面 BC
✅ M1.3b 滑移壁面 + 对称面 BC
✅ M1.3c 远场 Riemann BC（特征线入射/出射）
✅ M1.4a EOS 理想气体（p, T, a, H 从守恒量提取）
✅ M1.4b MUSCL-2 重构（k=1/3，minmod 限制器）
✅ M1.4c Roe 通量 + Harten 熵修复
✅ M1.4d 无粘 RHS 装配（3 方向法向，坐标变换）
✅ M1.5  粘性 RHS（Sutherland，中心差，绝热壁 BC）
✅ M1.6  TVD-RK3（Shu-Osher 3 步）
✅ M1.7a Solver step 组装（halo → BC → RHS → RK3）
✅ M1.7b test3 10 步端到端对拍（Fortran ref，≤ 1e-10）
✅ M1.8  GPU Benchmark（CPU 有记录）
✅ M1.9  jax.grad smoke test（梯度非零有限）
❌  Steger-Warming 通量（未实现，需 Phase 0 补）
❌  全局 JIT（当前每步单独 jit，未 lax.scan 化）
```

### 2.2 `1d_famrdp` 现状

```
✅ DomainManager（多块索引，ghost=4，CUT1TO1 映射）
✅ SchemeFactory（MUSCL-2 inter_index，2D+3D）
✅ ConvectiveTerm（Steger-Warming，i/j/k 三方向）
✅ TVD-RK3（3 步 Shu-Osher，true_cell_index）
✅ BoundaryManager（Wall/Symmetry/Farfield，2D+3D 索引）
✅ main.py（CPU/GPU 自动选择，dim 驱动 is_3d）
✅ 2D test3：4000 步，avg=0.055s，残差收敛
✅ 3D test3：4000 步，avg=0.105s，残差等价（挤压算例）
⚠️  远场 BC：简化版（ghost=自由流），非特征线
⚠️  粘性 RHS：骨架，函数体为 pass
⚠️  Roe 通量：未实现（当前仅 Steger-Warming）
❌  正式测试套件（无 pytest）
❌  jax.grad 验证
❌  Harten 熵修复
```

---

## 3. Phase 0：基础迁移与验证（目标周期：2–3 周）

**总目标**：两套代码在 `test3.grd` 上物理行为等价，各自通过 Fortran L1 验证，互相 L2 验证。

---

### M0.1 — 建立双路测试基础设施

**背景**：在正式开始物理模块开发前，需建立统一的测试框架，使两套代码共用同一算例、同一参考数据。

#### `famrdp_jax` 侧

- [ ] 整理 `validation/references/test3/` 目录结构：
  ```
  validation/references/test3/
    xyz/block_0001.bin, block_0002.bin   ← 坐标参考
    metric/block_0001_jac.bin            ← Jacobian 参考
    invscd/step_0001_rhs.bin             ← 单步无粘 RHS（Roe）
    invscd_sw/step_0001_rhs.bin          ← 单步无粘 RHS（Steger-Warming，待生成）
    viscous/step_0001_rhs.bin            ← 单步粘性 RHS
    steps/step_0010.bin                  ← 10 步整场（Roe）
    steps_sw/step_0010.bin               ← 10 步整场（Steger-Warming，待生成）
  ```
- [ ] 确认现有 56 个 pytest 全部通过：
  ```bash
  cd famrdp_jax && pytest tests/ -v --tb=short 2>&1 | tail -10
  ```

#### `1d_famrdp` 侧

- [ ] 新建 `1d_famrdp/tests/` 目录，配置 `conftest.py`：
  ```python
  # tests/conftest.py
  import sys, os
  sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
  import pytest
  import jax; jax.config.update("jax_enable_x64", True)
  try:
      jax.config.update("jax_default_device", jax.devices("gpu")[0])
  except RuntimeError:
      jax.config.update("jax_default_device", jax.devices("cpu")[0])
  ```
- [ ] 新建 `tests/test_smoke.py`，验证 2D/3D 各跑 5 步不报错：
  ```python
  def test_2d_smoke():
      # dim=2, 5 steps
      ...
  def test_3d_smoke():
      # dim=3, 5 steps
      ...
  ```
- [ ] 运行通过：`cd 1d_famrdp && pytest tests/ -v`

**验收**：两套代码各自的测试套件均可执行，无报错。

---

### M0.2 — Steger-Warming 参考数据生成（famrdp_jax 补充）

**背景**：`1d_famrdp` 使用 Steger-Warming 通量，为进行 L2 同通量对比，需在 `famrdp_jax` 中也实现 Steger-Warming，并生成对应的 Fortran 参考数据。

#### `famrdp_jax` 侧

- [ ] 在 `famrdp_jax/flux/` 新建 `steger_warming.py`，实现法向 Steger-Warming 通量：
  ```python
  def steger_warming_flux(qL, qR, n, gamma=1.4):
      """
      n: (3,) 面法向（已归一化）
      qL, qR: (5,) 原始变量 [rho, u, v, w, p]
      returns: (5,) 数值通量
      """
      def split_flux(q, sign):
          rho, u, v, w, p = q
          a = jnp.sqrt(gamma * p / rho)
          vn = u*n[0] + v*n[1] + w*n[2]
          lam = jnp.array([vn-a, vn, vn, vn, vn+a])
          lam_pm = 0.5 * (lam + sign * jnp.abs(lam))
          # 特征通量展开（见 Steger-Warming 1979）
          ...
      return split_flux(qL, +1) + split_flux(qR, -1)
  ```
- [ ] 测试：`pytest tests/test_steger_warming.py -v`，光滑初场单步通量非 NaN
- [ ] 在 `gen_refs.f90` 中新增 `--flux SW` 参数分支，生成 `invscd_sw/step_0001_rhs.bin`（或直接读取 `1d_famrdp` Fortran 的参考代码）

#### `1d_famrdp` 侧

- [ ] 新建 `tests/test_steger_warming.py`：
  ```python
  def test_sw_flux_nonneg_pressure():
      """纯右行激波：通量计算后密度和压力为正"""
      ...
  def test_sw_symmetric_uniform():
      """均匀流：左右通量之和等于解析通量"""
      ...
  ```

**验收**：`famrdp_jax` 中 Steger-Warming 实现通过单元测试；`gen_refs.f90` 可生成 SW 参考数据。

---

### M0.3 — 3D 无粘 RHS 对拍（L1 + L2）

**背景**：核心验证步骤。分两个子任务：各自对拍 Fortran（L1），然后互相对比（L2）。

#### `famrdp_jax` 侧（Roe 路径）

- [ ] 验证现有 `test_ten_steps_match` 通过（已完成 M1.7b）：
  ```bash
  pytest tests/validation/test_ten_steps_match.py -v
  ```
- [ ] 新增 `tests/validation/test_rhs_match_fortran.py`：
  ```python
  def test_invscd_rhs_roe_match():
      """单步 Roe RHS 与 Fortran gen_refs 结果 L∞ ≤ 1e-10"""
      ...
  ```

#### `famrdp_jax` 侧（Steger-Warming 路径，M0.2 完成后）

- [ ] 新增 `tests/validation/test_rhs_sw_match_fortran.py`：
  ```python
  def test_invscd_rhs_sw_match():
      """单步 SW RHS 与 Fortran SW 参考 L∞ ≤ 1e-10"""
      ...
  ```

#### `1d_famrdp` 侧

- [ ] 新建 `tests/test_rhs_l1.py`，对拍 `famrdp_jax` 侧已生成的 SW Fortran 参考：
  ```python
  REF_DIR = "../../validation/references/test3/invscd_sw"

  def test_invscd_rhs_step0_block0():
      """1D 路径 Steger-Warming RHS，单步，block 0，L∞ ≤ 1e-10"""
      buf, precomp = load_test3()
      rhs = compute_rhs_1d(buf, precomp, axis_all=True)
      ref = load_binary_ref(REF_DIR + "/step_0001_rhs_block0.bin")
      assert jnp.max(jnp.abs(rhs[:block0_size] - ref)) < 1e-10
  ```

#### L2 双路对比

- [ ] 新建 `tests/test_cross_validation.py`（位于 `famrdp_jax/tests/cross/`）：
  ```python
  def test_sw_rhs_cross():
      """
      famrdp_jax（SW）单步 RHS vs 1d_famrdp 单步 RHS
      对同一初场（均匀来流），全场 L∞ ≤ 1e-8
      """
      q0 = uniform_init(test3_shape)
      rhs_3d = famrdp_step_rhs_sw(q0)          # famrdp_jax SW 路径
      rhs_1d = famrdp_1d_step_rhs(q0_flat)     # 1d_famrdp SW 路径
      diff = compare_rhs_3d_to_1d(rhs_3d, rhs_1d, block_info)
      assert diff < 1e-8
  ```

**验收**：
- `famrdp_jax` Roe 路径：L1 ≤ 1e-10 ✅（已完成）
- `famrdp_jax` SW 路径：L1 ≤ 1e-10
- `1d_famrdp` SW 路径：L1 ≤ 1e-10
- 双路 L2 对比：同初场、同通量，L∞ ≤ 1e-8

---

### M0.4 — 粘性 RHS 实现（`1d_famrdp`）

**背景**：`1d_famrdp/terms/viscous.py` 当前为空骨架，参照 `famrdp_jax/physics/step/viscous_rhs.py` 移植。

#### `famrdp_jax` 侧（参考）

- [ ] 确认 `famrdp_jax` 粘性 RHS 对拍已通过（`tests/physics/test_viscous.py`）：
  ```bash
  pytest tests/physics/test_viscous.py -v
  ```
- [ ] 将 `famrdp_jax/physics/step/viscous_rhs.py` 的核心公式注释整理为接口文档，供 `1d_famrdp` 参考

#### `1d_famrdp` 侧

- [ ] 实现 `1d_famrdp/terms/viscous.py` 中 `_compute_duvwt`：
  ```python
  def _compute_duvwt(self, prime, temperature, grid_derivative, volume):
      """
      计算速度梯度 du/dx, dv/dx, dw/dx, dT/dx（共 12 个分量 per cell）
      使用 inter_index 做中心差插值，与对流项相同的 gather 机制
      """
      # 1. 用 center_index[axis] 插值 prime 到面
      # 2. 差分得到梯度
      # 3. 壁面 BC：无滑移速度梯度不变，温度梯度按绝热壁设为 0
      ...
  ```
- [ ] 实现 `_compute_viscous`：Sutherland 粘性律 + 粘性应力张量 + 热流
  ```python
  def _compute_viscous(self, prime, grid_derivative, vsl, duvwt, ref_data):
      """
      vsl: (N_total,) 动力粘性（Sutherland）
      duvwt: (N_total, 12) 速度梯度 + 温度梯度
      返回粘性 RHS，shape (N_total, 5)
      """
      ...
  ```
- [ ] 新建 `tests/test_viscous_l1.py`：
  ```python
  def test_viscous_zero_uniform():
      """均匀静止等温场：粘性 RHS 应为零（对称性）"""
      ...
  def test_viscous_vs_famrdp_jax():
      """test3 初场单步粘性 RHS，与 famrdp_jax 参考 L∞ ≤ 1e-8"""
      ...
  ```

**验收**：
- 均匀场粘性 RHS = 0（L∞ < 1e-14）
- test3 初场粘性 RHS 与 `famrdp_jax` 参考 L∞ ≤ 1e-8

---

### M0.5 — 远场 BC 精确化（`1d_famrdp`）

**背景**：`1d_famrdp` 的 `Farfield.py` 当前将 ghost 简单设为自由流，而 `famrdp_jax` 使用特征线 Riemann 不变量入射/出射。需对齐。

#### `famrdp_jax` 侧（参考）

- [ ] 梳理 `famrdp_jax/physics/bc/farfield.py` 的 Riemann 不变量实现：
  - 超声速入射：ghost = 自由流
  - 亚声速入射：`R+ = q_int - 2a/(γ-1)` 保持，`p_ghost` 由自由流给定
  - 亚声速出射：`R- = q_int + 2a/(γ-1)` 保持，`p_int` 由外部给定
  - 超声速出射：ghost = 内部值（零梯度外推）

#### `1d_famrdp` 侧

- [ ] 修改 `1d_famrdp/boundary/Farfield.py`，加入特征线判断：
  ```python
  def _riemann_farfield(prime_int, prime_face, n_face, farfield_var, gamma):
      """
      n_face: 面法向（3,）
      根据 Vn 和 a 判断入射/出射和超/亚声速，分 4 种情况设 ghost
      """
      vn = jnp.dot(prime_int[1:4], n_face)
      a  = jnp.sqrt(gamma * prime_int[4] / prime_int[0])
      Rp = vn + 2*a/(gamma-1)   # R+（出流不变量）
      Rm = farfield_var_vn - 2*farfield_a/(gamma-1)  # R-（入流不变量）
      ...
  ```
- [ ] 新建 `tests/test_bc_farfield.py`：
  ```python
  def test_farfield_subsonic_inflow():
      """亚声速入射：ghost 速度与 farfield 一致，压力由内插值决定"""
      ...
  def test_farfield_cross_vs_famrdp():
      """test3 初场远场面，与 famrdp_jax BC 结果 L∞ ≤ 1e-10"""
      ...
  ```

**验收**：
- 特征线远场 BC 4 种情况各有单元测试通过
- test3 远场面 prime 与 `famrdp_jax` L∞ ≤ 1e-10

---

### M0.6 — 端到端 10 步对拍（L2 双路）

**背景**：M0.3–M0.5 各物理量通过后，组装完整的端到端步进器进行全场比较。

#### 测试设计（共用）

- [ ] 新建 `cross_validation/test_end_to_end_cross.py`：

  ```python
  """
  双路端到端对比测试
  初场：test3 均匀来流
  步数：10 步
  通量：均使用 Steger-Warming
  容差：整场 L∞ ≤ 1e-8
  """
  
  def test_ten_steps_3d_vs_1d():
      q0_3d, q0_1d = load_uniform_init_both()
      
      # famrdp_jax：切换到 SW 路径
      q_3d = famrdp_step_n(q0_3d, flux='steger_warming', n=10)
      
      # 1d_famrdp：正常路径
      q_1d = famrdp_1d_step_n(q0_1d, n=10)
      
      # 转换布局：1D → 3D 对应位置
      q_1d_3d = reshape_1d_to_blocks(q_1d, block_info)
      
      diff = max_abs_diff(q_3d, q_1d_3d)
      assert diff < 1e-8, f"max diff = {diff:.3e}"
  ```

#### `famrdp_jax` 侧专项

- [ ] 确认 10 步 Roe 参考测试通过（已完成）
- [ ] 新增 10 步 SW 测试，对拍 SW 参考数据

#### `1d_famrdp` 侧专项

- [ ] 整场 10 步结果写入 `tests/fixtures/1d_step10_sw.npz`（基准快照）
- [ ] 新建 `tests/test_regression.py`：后续改动后快照不变，防止退步

**验收**（Phase 0 里程碑）：
- `famrdp_jax` 56 个原有测试仍全通过
- `1d_famrdp` 新增不少于 15 个测试全通过
- L2 双路对比：同通量 10 步 L∞ ≤ 1e-8
- 不同通量（Roe vs SW）稳态残差收敛曲线形状相似

---

### M0.7 — 可微性验证（`1d_famrdp`）

**背景**：`famrdp_jax` 的 `jax.grad` 已验证（M1.9）。`1d_famrdp` 需做同等测试。

#### `1d_famrdp` 侧

- [ ] 新建 `tests/test_grad.py`：
  ```python
  def test_grad_finite_nonzero():
      """jax.grad 穿透 1D step，梯度有限且非零"""
      buf, precomp = load_test3()
      
      def loss(prime):
          q_next = step_1d(buf._replace(prime=prime), precomp, dt=0.001)
          return jnp.sum(q_next.prime ** 2)
      
      g = jax.grad(loss)(buf.prime)
      assert jnp.isfinite(g).all()
      assert jnp.max(jnp.abs(g)) > 1e-30
  
  def test_grad_no_nan_from_gather():
      """gather 操作的反向传播不产生 NaN（inter_index 可能重复索引）"""
      ...
  ```
- [ ] 若出现 NaN：在 Roe / SW 奇点处加 `jnp.where` 守护

#### `famrdp_jax` 侧（对比）

- [ ] 记录 `jax.grad` 在 Roe 路径下的梯度量级（`max|g|`），与 1D 路径比较量级是否相近

**验收**：两套代码 `jax.grad` 均返回有限非零梯度；梯度量级在同一数量级内（1–2 个数量级差异可接受）。

---

## 4. Phase 1：生产级完整化（目标周期：2–3 周）

**总目标**：两套代码均能在 `test3.grd` 上完整运行至收敛，I/O 完整，性能有基准数据。

---

### M1.1 — CFL 自适应时间步（双路统一接口）

**背景**：当前两套代码都使用固定 `dtau`（全局 `param.json` 设定）。生产级计算需局部 CFL 限制。

#### 接口设计（两套共用）

```python
def compute_dt_global(prime, metrics, cfl, gamma=1.4):
    """
    逐单元计算 dt = cfl * vol / (|u_n|+c) * |J|
    返回全局 min（标量），在 JIT 内完成
    """
    rho = prime[:, 0]
    uvw = prime[:, 1:4]
    p   = prime[:, 4]
    a   = jnp.sqrt(gamma * p / rho)
    # 三方向法向速度
    vn_i = jnp.einsum('ni,ni->n', uvw, metrics.n_i)  # (N,)
    vn_j = jnp.einsum('ni,ni->n', uvw, metrics.n_j)
    vn_k = jnp.einsum('ni,ni->n', uvw, metrics.n_k)
    sr_i = jnp.abs(vn_i) + a * jnp.sqrt(jnp.sum(metrics.n_i**2, axis=1))
    sr_j = jnp.abs(vn_j) + a * jnp.sqrt(jnp.sum(metrics.n_j**2, axis=1))
    sr_k = jnp.abs(vn_k) + a * jnp.sqrt(jnp.sum(metrics.n_k**2, axis=1))
    dt_local = cfl * metrics.vol / (sr_i + sr_j + sr_k)
    return jnp.min(dt_local)
```

#### `famrdp_jax` 侧

- [ ] 在 `famrdp_jax/solver/step.py` 中集成 `compute_dt_global`
- [ ] 新建 `tests/test_cfl.py`：
  ```python
  def test_dt_positive():
      """局部 dt 严格正"""
      ...
  def test_dt_decreases_with_higher_cfl():
      """CFL=0.5 的 dt 是 CFL=1.0 的一半"""
      ...
  ```
- [ ] 验证：CFL=0.9 运行 100 步，残差收敛，无发散

#### `1d_famrdp` 侧

- [ ] 将 `simulation_manager.py` 中的 `spectral_radius` 函数改为 `compute_dt_global` 接口，接受 `cfl` 参数
- [ ] `param.json` 中 `"time_step_method": "CFL"` 模式生效
- [ ] 同等测试：`pytest tests/test_cfl.py -v`

#### L2 对比

- [ ] 同初场、同 CFL，两套代码的 `dt[step=0]` 差异 ≤ 1e-12

**验收**：两套代码 CFL 自适应步长计算一致，100 步测试无发散。

---

### M1.2 — 周期边界条件（双路）

**背景**：旋转对称算例（旋翼、级联翼型）需要周期 BC，这是 Phase 3 复杂算例的必要条件。

#### 接口设计

```
周期对 (block_A, face_A) ↔ (block_B, face_B)
实现：ghost 层直接映射到对侧真实单元（1D 法中与 CUT1TO1 机制相同）
     3D 法中在 halo_exchange 加 periodic 分支
```

#### `famrdp_jax` 侧

- [ ] 在 `famrdp_jax/mesh/halo.py` 中新增 `_periodic_exchange` 分支：
  ```python
  if conn_type == 'PERIODIC':
      qs[i] = _copy_face_periodic(qs[i], qs[j], face_i, face_j, ghost)
  ```
- [ ] 新建 `tests/mesh/test_periodic_bc.py`：正弦波平流，经过一个周期后波形不变（误差 ≤ 1e-6）

#### `1d_famrdp` 侧

- [ ] 在 `DomainManager._map_connection_face_to_ghost` 中处理 PERIODIC 类型：
  ```python
  if connection_type == 'PERIODIC':
      # 目标 ghost = 对侧真实单元全局编号（与 CUT1TO1 相同机制）
      ...
  ```
- [ ] 同等平流波测试：`tests/test_periodic.py`

**验收**：周期 BC 平流波 L2 误差 ≤ 1e-6（两套代码一致）。

---

### M1.3 — 残差监控与 I/O 统一

**背景**：两套代码当前残差输出格式不一致，需标准化以便对比分析。

#### 统一残差格式

```
STEP  res_rho    res_mom_u  res_mom_v  res_mom_w  res_energy  res_max
   0  x.xxxe+xx  ...        ...        ...        ...         x.xxxe+xx
```

#### `famrdp_jax` 侧

- [ ] 新建 `famrdp_jax/io/residual.py`：
  ```python
  def compute_residual(rhs, vol):
      """返回 5 分量 L2 残差和 L∞ 残差"""
      ...
  def print_residual(step, res5, res_max):
      print(f"{step:6d}  " + "  ".join(f"{r:.5e}" for r in res5) + f"  {res_max:.5e}")
  ```
- [ ] 将残差写入 CSV（`results/residual_3d.csv`），便于后续绘图对比

#### `1d_famrdp` 侧

- [ ] 统一 `simulation_manager.py` 中的残差打印为相同格式
- [ ] 写入 `results/residual_1d.csv`

#### 对比分析工具

- [ ] 新建 `tools/compare_residuals.py`：
  ```python
  """
  读取 residual_3d.csv 和 residual_1d.csv
  绘制双路残差收敛曲线，叠加在同一图上
  输出 results/residual_comparison.png
  """
  import matplotlib.pyplot as plt
  ...
  ```

**验收**：两套代码残差 CSV 格式相同；对比图可生成；两路残差曲线收敛趋势一致（步数相差 ≤ 10%）。

---

### M1.4 — Tecplot 输出（双路）

**背景**：可视化结果是用户确认物理正确性的最直观方式。

#### 输出规格

```
格式：Tecplot ASCII (.plt 或 .dat)
坐标：X, Y, Z（物理坐标）
变量：RHO, U, V, W, P（或 Ma, Cp）
每块一个 ZONE
```

#### `famrdp_jax` 侧

- [ ] 新建 `famrdp_jax/io/tecplot.py`：
  ```python
  def write_tecplot(filename, blocks, state, ref_data):
      """将多块解写入 Tecplot ASCII 格式"""
      with open(filename, 'w') as f:
          f.write('TITLE = "FAMRDP-JAX 3D Result"\n')
          f.write('VARIABLES = "X" "Y" "Z" "RHO" "U" "V" "W" "P"\n')
          for nb, (blk, met) in enumerate(zip(blocks, metrics)):
              ni, nj, nk = blk.q.shape[1:4]
              f.write(f'ZONE T="Block {nb+1}", I={ni}, J={nj}, K={nk}\n')
              ...
  ```
- [ ] 测试：写出后文件存在、行数与网格尺寸匹配

#### `1d_famrdp` 侧（现有 `output.py` 扩展）

- [ ] 扩展 `1d_famrdp/output.py` 的 `output_cell_data` 为标准 Tecplot 格式
- [ ] 在 `param.json` 中增加 `"output_interval": 500`，每 500 步自动输出

**验收**：两套代码各自能在 4000 步结束后输出 Tecplot 文件，圆柱绕流流场视觉上一致（速度分布、压力分布形态相同）。

---

### M1.5 — 性能基准测试（双路对比）

**背景**：量化两套方法的计算效率差异，与路线图中的预期吻合。

#### 测试方案

```
环境：CPU-only JAX（conda jax 环境）
算例：test3.grd（2 块，51×51×11 = 28602 真实单元）
步数：100 步（前 10 步 JIT 预热，后 90 步计时）
度量：avg step time (s)，total throughput (cells/s)
```

#### `famrdp_jax` 侧

```python
# benchmarks/bench_3d.py
import time
from famrdp_jax.solver import step_jit

def benchmark(n_warmup=10, n_bench=90):
    state = load_test3_state()
    for _ in range(n_warmup):
        state = step_jit(state, dt=0.001)
        state.q[0].block_until_ready()
    
    t0 = time.perf_counter()
    for _ in range(n_bench):
        state = step_jit(state, dt=0.001)
    state.q[0].block_until_ready()
    t1 = time.perf_counter()
    
    avg = (t1 - t0) / n_bench
    print(f"3D multi-block: {avg*1000:.2f} ms/step  ({28602/avg/1e6:.2f} M cells/s)")
```

#### `1d_famrdp` 侧

```python
# benchmarks/bench_1d.py（参照 simulation_manager.py 的计时逻辑）
# 输出格式与 3D 侧相同
```

#### 预期结果与记录

| 指标 | `famrdp_jax` (3D) | `1d_famrdp` (1D) |
|---|---|---|
| avg/step (2D) | 待测 | 0.055 s（已测） |
| avg/step (3D) | 待测 | 0.105 s（已测） |
| cells/s（3D） | 待测 | ~272K/s（估算） |

- [ ] 结果写入 `benchmarks/results/bench_YYYYMMDD.md`

**验收**：两套代码均有性能基准数据；记录两者的 cell throughput 比值；若 `famrdp_jax` 3D 比 `1d_famrdp` 3D 慢，记录原因（Python 循环、halo 开销等）。

---

### M1.6 — 长时收敛验证（L3 物理合理性）

**背景**：短时间步验证了数值精度，需要长时间运行验证物理正确性。

#### 测试规格

```
算例：test3.grd，Re=200，Ma=0.1，层流圆柱绕流
来流：无量纲化（ρ=1, u=1, v=0, p=1/γ）
运行步数：至少 10000 步（或残差 < 1e-6）
预期物理：周期性涡脱落（Karman 涡街），St ≈ 0.2
```

#### `famrdp_jax` 侧

- [ ] 运行 10000 步，每 100 步记录圆柱后缘升力系数 Cl
- [ ] 绘制 Cl(t) 曲线，确认周期性振荡
- [ ] 计算 Strouhal 数：`St = f*D/U∞`，期望 St ∈ [0.18, 0.22]
- [ ] 保存：`results/3d_cl_history.csv`

#### `1d_famrdp` 侧

- [ ] 同等步数运行，记录 Cl 历史
- [ ] 保存：`results/1d_cl_history.csv`

#### 对比

- [ ] 两路 Cl 振荡频率差异 ≤ 5%（St 数接近）
- [ ] 两路 Cl 振幅差异 ≤ 10%（Roe vs SW 导致的格式差异可接受）

**验收**（Phase 1 里程碑）：两套代码各自能正确捕捉 Karman 涡街，St 数在物理范围内，两路对比差异在容差内。

---

## 5. 测试目录结构

Phase 0 + Phase 1 完成后，测试文件组织如下：

```
famrdp-jax/
├── famrdp_jax/
│   └── tests/
│       ├── mesh/
│       │   ├── test_xyz_match_fortran.py          ✅ 已有
│       │   ├── test_metric_match_fortran.py        ✅ 已有
│       │   └── test_periodic_bc.py                 M1.2 新增
│       ├── physics/
│       │   ├── test_steger_warming.py              M0.2 新增
│       │   ├── test_rhs_match_fortran.py           M0.3 新增
│       │   ├── test_rhs_sw_match_fortran.py        M0.3 新增
│       │   ├── test_cfl.py                         M1.1 新增
│       │   └── test_viscous.py                     ✅ 已有
│       └── validation/
│           ├── test_ten_steps_match.py             ✅ 已有
│           └── test_ten_steps_sw_match.py          M0.6 新增
│
├── 1d_famrdp/
│   └── tests/
│       ├── conftest.py                             M0.1 新增
│       ├── test_smoke.py                           M0.1 新增
│       ├── test_steger_warming.py                  M0.2 新增
│       ├── test_rhs_l1.py                          M0.3 新增
│       ├── test_viscous_l1.py                      M0.4 新增
│       ├── test_bc_farfield.py                     M0.5 新增
│       ├── test_cross_validation.py                M0.6 新增
│       ├── test_grad.py                            M0.7 新增
│       ├── test_cfl.py                             M1.1 新增
│       ├── test_periodic.py                        M1.2 新增
│       └── test_regression.py                      M0.6 新增
│
├── cross_validation/
│   ├── test_rhs_cross.py                           M0.3 新增
│   ├── test_end_to_end_cross.py                    M0.6 新增
│   └── compare_residuals.py                        M1.3 新增
│
├── benchmarks/
│   ├── bench_3d.py                                 M1.5 新增
│   ├── bench_1d.py                                 M1.5 新增
│   └── results/
│
└── validation/
    └── references/test3/
        ├── xyz/                                    ✅ 已有
        ├── metric/                                 ✅ 已有
        ├── invscd/                                 ✅ 已有（Roe）
        ├── invscd_sw/                              M0.2 生成（SW）
        └── steps/                                  ✅ 已有（10步，Roe）
            └── steps_sw/                           M0.6 生成（10步，SW）
```

---

## 6. 里程碑完成标准汇总

| 里程碑 | `famrdp_jax` 任务 | `1d_famrdp` 任务 | L2 双路验证 | 预计周期 |
|---|---|---|---|---|
| M0.1 测试基础设施 | 确认 56 测试通过 | 新建 tests/ + smoke test | — | 0.5 天 |
| M0.2 Steger-Warming | 实现 SW 通量 + 生成 SW 参考 | SW 单元测试 | — | 2 天 |
| M0.3 无粘 RHS 对拍 | SW 路径 L1 测试 | SW L1 测试 | L∞ ≤ 1e-8 | 3 天 |
| M0.4 粘性 RHS | 现有测试确认 | 实现 viscous.py + L1 测试 | L∞ ≤ 1e-8 | 4 天 |
| M0.5 远场 BC | 现有测试确认 | 特征线 Farfield + L1 测试 | L∞ ≤ 1e-10 | 2 天 |
| M0.6 端到端对拍 | SW 10 步测试 | 回归快照 + 双路 10 步 | L∞ ≤ 1e-8 | 2 天 |
| M0.7 可微性 | 确认现有 grad 测试 | grad 测试 + NaN 守护 | 量级相近 | 1 天 |
| M1.1 自适应 CFL | 实现 + 测试 | 实现 + 测试 | dt 一致 ≤ 1e-12 | 2 天 |
| M1.2 周期 BC | 实现 + 平流测试 | 实现 + 平流测试 | 误差 ≤ 1e-6 | 3 天 |
| M1.3 残差 I/O | 统一格式 + CSV | 统一格式 + CSV | 同格式可对比 | 1 天 |
| M1.4 Tecplot 输出 | tecplot.py | 扩展 output.py | 视觉一致 | 2 天 |
| M1.5 性能基准 | bench_3d.py | bench_1d.py | 吞吐比有记录 | 1 天 |
| M1.6 长时收敛 | 10000步 + Cl 历史 | 10000步 + Cl 历史 | St 差异 ≤ 5% | 2 天 |

**Phase 0 总验收**：`1d_famrdp` 测试套件 ≥ 15 个通过；双路 SW 10 步 L∞ ≤ 1e-8；`jax.grad` 两路均可用。  
**Phase 1 总验收**：两套代码圆柱绕流 St 数均在 [0.18, 0.22]；benchmark 数据完备；I/O 完整。

---

*文档版本 v1.0 — 2026-04-18*  
*下一步：按里程碑顺序逐步实施，每个里程碑完成后更新本文档对应条目*
