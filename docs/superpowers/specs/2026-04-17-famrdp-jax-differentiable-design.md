# FAMRDP → JAX 可微分 CFD 求解器 · 设计规范

- **日期**:2026-04-17
- **分支**:`famrdp-jax`(worktree)
- **状态**:设计已通过 brainstorming,等待实施计划

---

## 1. 背景与目标

### 1.1 原项目

`FAMRDP` / `HOSTA`(High-Order Simulator for Aerodynamics)是一个 Fortran 90 + MPI 的高阶有限体积气动求解器。现有代码特性:

- 3D 结构多块网格
- 插值/重构:MUSCL-2、WCNS5/7、HDCS5/7、多种紧致格式(SCS/DCS)
- 通量:Steger-Warming、Van Leer、Roe、Roe+preconditioning、SLAU、Roe-scmp
- 时间推进:显式 Runge-Kutta、隐式 LUSGS
- 湍流模型:Spalart-Allmaras、Menter SST
- MPI 块间通信

### 1.2 本次目标

用 JAX 重新实现**一个受控子集**的求解流程,使其:

1. **可微分**:`jax.grad` 可贯穿整个时间步,为后续伴随/优化/反演/闭包学习留出接口
2. **可在 GPU 上 benchmark**
3. **与 Fortran 原代码逐步对拍**,在对应格式/配置下数值结果一致(容差 1e-10~1e-12)

**本次不承诺具体科学应用**(不解决特定优化/反演任务);目标是交付一个**可微 CFD 框架基线**。

### 1.3 非目标(明确不做)

- 🚫 湍流模型(SA / Menter)
- 🚫 LUSGS / 隐式时间推进
- 🚫 MPI / 多 GPU
- 🚫 非结构或重叠网格
- 🚫 `sol_prmat` / `sol_prsca` 以外的预处理
- 🚫 声学 / 平均场 / monitor 输出

后续如需,走新一轮 brainstorming + spec。

---

## 2. 范围锁定

| 维度 | 范围 |
|---|---|
| 物理 | 3D 可压缩 Navier-Stokes,层流 |
| 网格 | 结构多块(Python `list[Block]`,异构 shape) |
| 重构(最终) | `dcsh5pi`(5 阶紧致差分,primitive 插值) |
| 通量(最终) | `Roe-scmp`(低马数简化预处理 Roe) |
| 重构(阶段 1 脚手架) | MUSCL-2 |
| 通量(阶段 1 脚手架) | Roe |
| 时间推进 | 显式 3 阶 TVD Runge-Kutta |
| 并行 | 单 device(CPU 或 1 GPU) |
| 精度 | `float64`(对齐 Fortran `-r8`) |
| 可微性 | 整个 `step(State, dt) -> State'` 可微 |

### 2.1 跨语言 bit-for-bit 风险(已接受)

gfortran + JAX/XLA 在浮点表达式的求和顺序、FMA 融合、GPU 并行归约上**不可能**产生同 bit 结果。可达目标:

- 前 10-100 步,同格式同初值下,残差 diff ≤ 1e-12(float64)
- 长时间积分因混沌放大差异指数增长(物理事实,非 bug)

---

## 3. 总体架构

### 3.1 目录结构

与原 Fortran `src/` 并存,不破坏原 build。

```
famrdp_jax/
├── core/
│   ├── types.py              # Block / State / Metrics NamedTuple(pytree 注册)
│   ├── config.py             # param.inp → dataclass
│   └── constants.py          # nflux_*, nintnon_* 枚举
├── mesh/
│   ├── io_grd.py             # 读 test3.grd / test3.inp
│   ├── metric.py             # 雅可比、度量系数(镜像 metric.f90)
│   └── halo.py               # 块间 ghost 拷贝
├── physics/
│   ├── eos.py                # 理想气体 / Sutherland
│   ├── bc.py                 # BC 分发 + 各 BC 实现
│   └── reconstruct/
│       ├── muscl2.py         # 阶段 1
│       └── dcsh5pi.py        # 阶段 2
├── flux/
│   ├── roe.py                # 阶段 1
│   └── roe_scmp.py           # 阶段 2
├── rhs/
│   ├── invscd.py             # 无粘通量散度
│   └── viscous.py            # 粘性通量
├── time/
│   └── rk.py                 # 显式 RK
├── solver/
│   └── step.py               # 单步编排:halo → BC → RHS → RK
├── io/
│   └── tecplot.py
├── validation/
│   ├── fortran_ref.py        # 调 Fortran + 解析二进制
│   ├── diff.py
│   └── bench.py
└── tests/
```

### 3.2 设计原则

1. **Fortran 模块 1:1 对应**:`rhs_invscd_flux_roe` ↔ `flux.roe.compute`,便于侧向对拍
2. **pytree 化**:`Block`、`Metrics`、`State` 用 `flax.struct.dataclass`,`jit`/`grad` 自然穿过
3. **Stage-1 / Stage-2 格式可替换**:上层不感知,reconstruct/flux 作为函数指针传入

---

## 4. 核心数据结构

```python
# core/types.py

@dataclass(frozen=True)
class Block:
    """单个结构块。ghost 嵌入在 q 内,ni = ni_interior + 2*ghost。"""
    q:       Float[Array, "nvar ni nj nk"]   # 守恒量 (ρ, ρu, ρv, ρw, ρE)
    xyz:     Float[Array, "3 ni nj nk"]

@dataclass(frozen=True)
class Metrics:
    jac:     Float[Array, "ni nj nk"]
    kxyz:    Float[Array, "3 3 ni nj nk"]
    vol:     Float[Array, "ni nj nk"]

@dataclass(frozen=True)
class BlockTopology:
    block_id:     int
    neighbors:    dict[Face, tuple[int, Face, Orient]]
    bc_type:      dict[Face, BCType]

@dataclass(frozen=True)
class State:
    blocks:   tuple[Block, ...]
    metrics:  tuple[Metrics, ...]
    t:        Float
    step:     Int

@dataclass(frozen=True)
class Config:
    gas:         GasModel
    scheme:      SchemeChoice
    bc_params:   dict
    ghost:       int = 3         # 阶段 2 dcsh5pi 需要;阶段 1 也按 3 分配,兼容
```

**约定**:

1. **Ghost 嵌入在 `q` 内**:halo exchange = 切片 `.at[].set()`,不动形状
2. **`State.blocks` 是 `tuple`**:pytree 结构固定,块数静态,块间 shape 可异构
3. **`Config` 不进 pytree**:作 `static_argnames`,避免常数变动触发重编译
4. **`BlockTopology`** 是静态元数据,驱动 halo 索引

---

## 5. 单步数据流与 JIT 边界

### 5.1 单步流水线

```
State_n
  │
  ├─► halo_exchange(blocks)                  ← Python 外层
  ├─► apply_bc(blocks, config)               ← Python 循环,每 BC jit
  ├─► for each block (Python):
  │     rhs_block = compute_rhs(q, metrics, config)   ← JIT ①
  └─► rk_update(q_tuple, rhs_tuple, dt, cfg.rk)       ← JIT ②
```

### 5.2 JIT 边界表

| 层级 | 是否 JIT | 理由 |
|---|---|---|
| `compute_rhs(block, metrics, config)` | ✅ `@jit`,config static | 最重数值核,GPU 主要受益点 |
| `apply_bc_single(q, face_info, config)` | ✅ 每 BC 一 jit | 形状规整,缓存友好 |
| `rk_update` | ✅ tuple 一次 jit | 纯仿射,编译轻 |
| `halo_exchange` | ❌ Python + `.at[].set()` | 控制流灵活,topology 静态但稀疏 |
| 时间主循环 | ❌ Python | 保留逐步对拍能力 |

### 5.3 伪代码

```python
def step(state: State, cfg: Config, dt: float) -> State:
    blocks = halo_exchange(state.blocks, cfg.topology)
    blocks = tuple(apply_bc(b, cfg) for b in blocks)
    rhs    = tuple(compute_rhs(b.q, m, cfg) for b, m in
                   zip(blocks, state.metrics))
    new_q  = rk_update(tuple(b.q for b in blocks), rhs, dt, cfg.rk)
    return state.replace(blocks=replace_q(blocks, new_q),
                         t=state.t+dt, step=state.step+1)
```

**可微性**:所有 JIT 函数均为纯函数;`.at[].set()` 是 scatter,VJP 存在。整步 `step` 可微。

---

## 6. 精度 / 错误处理 / 可复现性

### 6.1 精度

```python
jax.config.update("jax_enable_x64", True)
```

- 默认 `float64`;整型 `int32`
- 禁止隐式降精度
- `core/types.py` 构造 `Block` 时校验 dtype

### 6.2 错误处理

| 类别 | 策略 |
|---|---|
| 不变量 | Python `assert` / `raise`,不在 JIT 内 |
| 数值异常(NaN / 负密度) | RK 子步末尾 `jnp.isfinite` 检查,失败走 `jax.debug.callback` |
| JIT 内部 | 不写 `raise` / `print` |
| IO 错误 | Python 标准异常 |
| 可微 NaN 风险点 | Roe 熵修正 / limiter / sqrt / 除密度:`jnp.where` 守护 |

**开关**:`Config.nan_check: bool`,默认开,benchmark 关。

### 6.3 可复现性

- 无随机数;未来预留 `PRNGKey` 接口
- 可选 `XLA_FLAGS=--xla_gpu_deterministic_ops=true`,默认不开
- 对拍 harness 固定 Fortran commit + param.inp
- 版本锁:`uv.lock` / `requirements.txt`

### 6.4 日志

- Python `logging`,残差输出对齐 Fortran `res.plt`
- `jax.debug.print` 受 `Config.debug` 控制
- `Config.dump_on_nan = True` 时 pickle state 到 `./dumps/`

### 6.5 设备

- `solver.step.make_step(config, device="gpu")` 包装 `jax.jit(backend=...)`
- 测试默认 CPU,benchmark 强制 GPU
- 多 GPU 非目标

---

## 7. 交付里程碑

### 7.1 阶段 1:低阶脚手架(MUSCL + Roe)

| ID | 里程碑 | 验收 |
|---|---|---|
| M1.0 | 读 `.grd` / `.inp` → State | xyz/topology 对拍 1e-14 |
| M1.1 | `metric.py` | Metrics diff < 1e-13 |
| M1.2 | `halo_exchange` | 两块 cube unit test |
| M1.3 | BC(滑移、无滑移、远场、对称、块界面) | 每 BC 单独 unit test + 对拍 1e-13 |
| M1.4 | EOS + MUSCL-2 + Roe | 冻结流场 rhs diff < 1e-12 |
| M1.5 | `rhs_viscous` | 单步对拍 1e-12 |
| M1.6 | 显式 3 阶 TVD-RK | 纯时间算子单元测试 |
| M1.7 | 端到端 `test3.grd` 10 步对拍 | 整场 diff < 1e-10 |
| M1.8 | GPU benchmark | 输出 Markdown 表 |
| M1.9 | 可微性 smoke test(极简) | `jax.grad` 返回非零有限值 |

### 7.2 阶段 2:替换为 `dcsh5pi + Roe-scmp`

| ID | 里程碑 | 验收 |
|---|---|---|
| M2.1 | `reconstruct/dcsh5pi.py` | `tridiagonal_solve` 系数对拍 1e-13 |
| M2.2 | 紧致格式的 halo 扩展 / 边界退化 | 与 Fortran 近界处理一致 |
| M2.3 | `flux/roe_scmp.py` | 单 flux diff 1e-13 |
| M2.4 | 端到端 `test3.grd` 10 步对拍(对齐 param.inp) | 整场 diff < 1e-10 |
| M2.5 | GPU benchmark(阶段 2) | 与阶段 1 对比报告 |
| M2.6 | 可微性 smoke test(阶段 2) | 同 M1.9 |

---

## 8. 测试与对拍 Harness

### 8.1 测试金字塔

| 层级 | 数量级 | 运行时机 |
|---|---|---|
| 单元测试 | ~50+ | 每次提交 |
| 模块对拍 | ~20 | 每次提交(慢的 mark skip) |
| 集成对拍(端到端 N 步) | ~3 | PR 前手动 |
| GPU benchmark | 1 脚本 | M1.8 / M2.5 |
| 可微 smoke test | 1 | M1.9 / M2.6 |

### 8.2 Fortran 对拍 harness

**策略**:修改 Fortran 源码加 probe,CMake flag 控制编译。

```
validation/
├── fortran_probes/
│   ├── probe_metric.f90
│   ├── probe_rhs.f90
│   ├── probe_step.f90
│   └── probes.h              # #ifdef PROBE_xxx
├── references/               # 全部二进制入 git
│   └── test3_10steps/
│       ├── metric.bin
│       ├── step_0001.bin
│       └── ...
├── fortran_ref.py
├── diff.py
└── bench.py
```

**实现要点**:

- Probe 通过 `cmake -DPROBE=ON` 启用;release build 不含
- Fortran `unformatted write`(带 4-byte record marker),Python 用 `scipy.io.FortranFile` 读
- `diff.py` 按字段独立容差(密度/压强 1e-12,残差 1e-10);失败时打印 top-10 diff 坐标

### 8.3 Reference 管理

- **所有 reference 二进制入 git**(用户明确要求)
- `scripts/regenerate_references.sh` 从当前 Fortran commit 重新生成
- PR 中若 reference 有更新,说明触发条件

### 8.4 Benchmark 设计

```
validation/bench.py → docs/benchmarks/M1.8.md

| case       | ncells  | JAX GPU (ms/step) | Fortran CPU (ms/step) | ratio |
|------------|---------|-------------------|-----------------------|-------|
| cube_32^3  | 32K     | ...               | ...                   | ...   |
```

- warmup 5 步 + 100 步中位数
- 表头记录 GPU 型号、JAX/jaxlib/CUDA 版本
- 不做 multi-run 统计

---

## 9. 依赖

- Python ≥ 3.11
- `jax` + `jaxlib-cuda12`
- `flax`(用 `flax.struct.dataclass`)
- `chex`(测试断言)
- `numpy`、`scipy`(读 Fortran unformatted)
- `pytest`、`hypothesis`
- `uv` 做 lock

---

## 10. 后续(超出本 spec)

- 湍流模型(SA / Menter)
- 隐式 LUSGS + 自定义 VJP
- MPI / `shard_map` 多 GPU
- 具体应用:伴随形状优化、数据同化、闭包学习

每一项独立重启 brainstorming。
