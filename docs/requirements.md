# FAMRDP → JAX 移植需求文档

- **日期**:2026-04-17
- **关联 spec**:[`docs/superpowers/specs/2026-04-17-famrdp-jax-differentiable-design.md`](superpowers/specs/2026-04-17-famrdp-jax-differentiable-design.md)
- **状态**:等待审阅

本文档把 spec 的设计决策收敛为一份**带唯一 ID、可核查**的需求清单,供实施计划与验收使用。

- **FR**(Functional Requirement):功能需求
- **NFR**(Non-Functional Requirement):非功能需求(性能、精度、可维护性)
- **CR**(Constraint):约束(不做什么 / 外部依赖)
- 每条都标注**优先级**(P0 必须 / P1 重要 / P2 可延后)、**阶段**(S1 / S2)与**验收里程碑**(M1.x / M2.x)

---

## 1. 功能需求(FR)

### 1.1 网格与拓扑

| ID | 需求 | 优先级 | 阶段 | 验收 |
|---|---|---|---|---|
| FR-01 | 读取 Fortran 二进制 `.grd` 网格文件,构造 `State.blocks` 的 `xyz` 坐标 | P0 | S1 | M1.0 |
| FR-02 | 读取 `.inp` topology 文件,还原 `BlockTopology` 中的块间邻接关系与面 BC 类型 | P0 | S1 | M1.0 |
| FR-03 | 支持 3D 结构多块,块数静态、块间 shape 可异构 | P0 | S1 | M1.0 |
| FR-04 | Ghost 层数可配,阶段 1 至少 2,阶段 2 至少 3 | P0 | S1/S2 | M1.2 / M2.2 |
| FR-05 | 块间 halo exchange,支持 6 面每面指向邻块某一面(含 orient) | P0 | S1 | M1.2 |

### 1.2 度量系数

| ID | 需求 | 优先级 | 阶段 | 验收 |
|---|---|---|---|---|
| FR-10 | 由 `xyz` 预计算雅可比 `jac`、度量 `kxyz`、体积 `vol` | P0 | S1 | M1.1 |
| FR-11 | 度量计算与 Fortran `metric.f90` 逐点 diff ≤ 1e-13 | P0 | S1 | M1.1 |

### 1.3 边界条件

| ID | 需求 | 优先级 | 阶段 | 验收 |
|---|---|---|---|---|
| FR-20 | 无滑移壁面(指定壁温) | P0 | S1 | M1.3 |
| FR-21 | 滑移壁面 | P0 | S1 | M1.3 |
| FR-22 | 对称面 | P0 | S1 | M1.3 |
| FR-23 | 远场 / 特征线入流-出流 | P0 | S1 | M1.3 |
| FR-24 | 块界面(ghost 填充来自邻块) | P0 | S1 | M1.2 / M1.3 |
| FR-25 | 每种 BC 与 Fortran 对应 kernel 单步对拍 ≤ 1e-13 | P0 | S1 | M1.3 |

### 1.4 物理与状态方程

| ID | 需求 | 优先级 | 阶段 | 验收 |
|---|---|---|---|---|
| FR-30 | 理想气体状态方程(γ、气体常数由 config 注入) | P0 | S1 | M1.4 |
| FR-31 | Sutherland 粘性律(μ₀、T₀、Tₛ) | P0 | S1 | M1.5 |
| FR-32 | 守恒 ↔ 原始变量相互转换 | P0 | S1 | M1.4 |

### 1.5 空间重构 / 通量

| ID | 需求 | 优先级 | 阶段 | 验收 |
|---|---|---|---|---|
| FR-40 | MUSCL-2 重构(阶段 1 脚手架) | P0 | S1 | M1.4 |
| FR-41 | Roe 通量(阶段 1 脚手架) | P0 | S1 | M1.4 |
| FR-42 | 中心差粘性通量 | P0 | S1 | M1.5 |
| FR-43 | `dcsh5pi` 5 阶紧致差分重构(阶段 2 最终目标) | P0 | S2 | M2.1 |
| FR-44 | Roe-scmp(低马数简化预处理 Roe)通量(阶段 2 最终目标) | P0 | S2 | M2.3 |
| FR-45 | 紧致格式在块界面近处退化处理,与 Fortran 行为一致 | P0 | S2 | M2.2 |

### 1.6 时间推进

| ID | 需求 | 优先级 | 阶段 | 验收 |
|---|---|---|---|---|
| FR-50 | 显式 3 阶 TVD Runge-Kutta | P0 | S1 | M1.6 |
| FR-51 | 全局恒定时间步(CFL 由 config 指定) | P0 | S1 | M1.6 |
| FR-52 | 步级 Python 外循环(不 `lax.scan` 主循环),保留逐步诊断能力 | P0 | S1 | M1.7 |

### 1.7 输入输出

| ID | 需求 | 优先级 | 阶段 | 验收 |
|---|---|---|---|---|
| FR-60 | 解析 `param.inp`(Fortran namelist)到 `Config` dataclass | P0 | S1 | M1.0 |
| FR-61 | 每步残差日志对齐 Fortran `res.plt` 口径(分量顺序、格式) | P1 | S1 | M1.7 |
| FR-62 | 输出 Tecplot `.plt`(ASCII 或 binary,与 Fortran `pltsty` 选项对齐其一即可) | P1 | S1 | M1.7 |
| FR-63 | NaN 检测失败时 dump state pickle 到 `./dumps/`(`Config.dump_on_nan`) | P2 | S1 | M1.7 |

### 1.8 可微性

| ID | 需求 | 优先级 | 阶段 | 验收 |
|---|---|---|---|---|
| FR-70 | `step(State, cfg, dt) -> State'` 整体对 `State` 和 `cfg` 中连续参数可微 | P0 | S1/S2 | M1.9 / M2.6 |
| FR-71 | `jax.grad` 对某标量目标(如某分量积分)返回**非零有限值** | P0 | S1/S2 | M1.9 / M2.6 |
| FR-72 | Roe 熵修正、limiter、开根号、除密度等潜在 NaN 点用 `jnp.where` 守护 | P0 | S1/S2 | M1.4 / M2.3 |

---

## 2. 非功能需求(NFR)

### 2.1 精度

| ID | 需求 | 优先级 |
|---|---|---|
| NFR-01 | 所有计算使用 `float64`(对齐 Fortran `-r8`) | P0 |
| NFR-02 | 程序入口必须 `jax.config.update("jax_enable_x64", True)` | P0 |
| NFR-03 | 构造 `Block` 时运行期校验 dtype,不通过即 raise | P1 |

### 2.2 性能

| ID | 需求 | 优先级 | 验收 |
|---|---|---|---|
| NFR-10 | 主数值核 (`compute_rhs`、`rk_update`、`apply_bc`) 在单 GPU 上 JIT 编译 | P0 | M1.8 |
| NFR-11 | 提交 GPU benchmark 报告:warmup 5 步 + 100 步中位数,含 JAX/jaxlib/CUDA 版本与 GPU 型号 | P0 | M1.8 / M2.5 |
| NFR-12 | 记录每阶段与单核 Fortran 的吞吐比值 ratio | P1 | M1.8 / M2.5 |
| NFR-13 | `Config.nan_check=False` 时关闭 isfinite 归约,benchmark 默认关 | P1 | M1.8 |

### 2.3 对拍容差

| ID | 对象 | 容差 | 验收 |
|---|---|---|---|
| NFR-20 | 网格 `xyz` vs Fortran | 绝对 1e-14 | M1.0 |
| NFR-21 | Metrics 全字段 vs Fortran | 绝对 1e-13 | M1.1 |
| NFR-22 | 单 BC kernel 单步输出 | 绝对 1e-13 | M1.3 |
| NFR-23 | RHS 单步(冻结流场) | 绝对 1e-12 | M1.4 / M2.3 |
| NFR-24 | 端到端 `test3.grd` 10 步整场 | 绝对 1e-10 | M1.7 / M2.4 |
| NFR-25 | 跨语言**不**追求 bit-for-bit;长时间积分允许指数发散 | — | — |

### 2.4 错误处理

| ID | 需求 |
|---|---|
| NFR-30 | JIT 内部**禁止** `raise` / `print` |
| NFR-31 | RK 子步末尾做 `jnp.isfinite(q)` + 正压 / 正密度检查,失败走 `jax.debug.callback` 停机 |
| NFR-32 | IO / 不变量错误走 Python 标准异常 |

### 2.5 可复现性

| ID | 需求 |
|---|---|
| NFR-40 | 依赖版本由 `uv.lock` / `requirements.txt` 锁死 |
| NFR-41 | Reference 二进制全部入 git,版本与生成时的 Fortran commit 对齐 |
| NFR-42 | 可选 `XLA_FLAGS=--xla_gpu_deterministic_ops=true`,默认不开 |

### 2.6 可维护性

| ID | 需求 |
|---|---|
| NFR-50 | Python 模块文件与 Fortran 源文件 1:1 对应,函数名尽量可侧向对照 |
| NFR-51 | `reconstruct/*` 与 `flux/*` 是可替换的函数指针,不动上层 step |
| NFR-52 | 所有数据结构(`Block` / `Metrics` / `State`)为 pytree,`jit`/`grad` 自然穿过 |

---

## 3. 约束(CR)

### 3.1 范围约束

| ID | 约束 |
|---|---|
| CR-01 | **不**实现湍流模型(SA / Menter) |
| CR-02 | **不**实现隐式时间推进(LUSGS) |
| CR-03 | **不**支持 MPI / 多 GPU / `pmap` / `shard_map` |
| CR-04 | **不**支持非结构或重叠网格 |
| CR-05 | **不**实现 `sol_prmat` / `sol_prsca` 独立预处理器 |
| CR-06 | **不**实现声学 / 平均场 / monitor 输出 |

### 3.2 技术约束

| ID | 约束 |
|---|---|
| CR-10 | Python ≥ 3.11 |
| CR-11 | JAX + jaxlib-CUDA12;`flax.struct.dataclass`、`chex`、`scipy.io.FortranFile` |
| CR-12 | 包 `famrdp_jax/` 与现有 Fortran `src/` 并存,不破坏现有 build |
| CR-13 | Fortran 端加 probe(`#ifdef PROBE_xxx`,CMake `-DPROBE=ON`),release build 不含 |

### 3.3 协作约束

| ID | 约束 |
|---|---|
| CR-20 | 每阶段里程碑按顺序完成,前置未通过不得开工后续 |
| CR-21 | 新增 Fortran reference 二进制时需在 PR 中说明触发原因 |

---

## 4. 可追溯矩阵(需求 ↔ 里程碑)

| 里程碑 | 覆盖的需求 |
|---|---|
| M1.0 | FR-01, FR-02, FR-03, FR-60, NFR-20 |
| M1.1 | FR-10, FR-11, NFR-21 |
| M1.2 | FR-04, FR-05, FR-24 |
| M1.3 | FR-20 ~ FR-25, NFR-22 |
| M1.4 | FR-30, FR-32, FR-40, FR-41, FR-72, NFR-23 |
| M1.5 | FR-31, FR-42 |
| M1.6 | FR-50, FR-51 |
| M1.7 | FR-52, FR-61, FR-62, FR-63, NFR-24 |
| M1.8 | NFR-10 ~ NFR-13 |
| M1.9 | FR-70, FR-71 |
| M2.1 | FR-43 |
| M2.2 | FR-04, FR-45 |
| M2.3 | FR-44, FR-72, NFR-23 |
| M2.4 | NFR-24 |
| M2.5 | NFR-10 ~ NFR-13 |
| M2.6 | FR-70, FR-71 |

---

## 5. 变更流程

本文档与 spec 文档同源。修订流程:

1. 发起变更提议(新增/修改/删除需求条目)
2. 更新 spec 对应节
3. 更新本文档并保持 ID 稳定(作废条目标 `~~FR-xx~~ (DEPRECATED <日期>)`,不重用 ID)
4. 更新 §4 可追溯矩阵
5. commit 消息明确写出 ID 列表
