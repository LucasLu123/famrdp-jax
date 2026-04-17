# HOSTA / FAMRDP 功能说明文档

- **对象**:原 Fortran 90 + MPI 求解器(`src/*.f90` 共约 73 000 行)
- **目的**:系统描述现有代码的功能、数据流、模块职责,作为 JAX 移植前的基线参考
- **日期**:2026-04-17

---

## 0. 一页总览

**HOSTA**(High-Order Simulator for Aerodynamics)是一套**结构多块、高阶精度、MPI 并行**的 3D 可压缩 Navier-Stokes 求解器,面向外气动/飞行器设计。主要特性:

| 类别 | 支持的选项 |
|---|---|
| 维度 | 1D(轴对称)、2D、3D |
| 物理 | Euler / NS 层流 / RANS(SA、Menter SST、HST) |
| 重构/插值 | 24 种:MUSCL-2、WCNS-5/7、HDCS-5/7、SCSL-3/5、SCSH-3/5、SCSN-2/3/4、DCSH-5/7 |
| 通量 | 7 种:Steger-Warming、SW-mod、Van Leer、Roe、Roe+Precond、SLAU、Roe-SCMP |
| 时间推进 | RK3 显式、LUSGS(标量/矩阵/预处理/SCMP)、PRSGS(标量/矩阵/split/unsplit/SCMP) |
| 网格 | 结构多块(PLOT3D `.grd`)、支持对偶网格(`ncelcet=1`)、周期/极点奇异处理 |
| 并行 | MPI(`PARALLEL`)、OpenMP(`OMP_IMP`)、可混合 |

---

## 1. 程序执行流程

### 1.1 顶层调用链(`main.f90`)

```
hosta program
  ├─ env_initialize       ! MPI_Init, 设置 myid, numprocs
  ├─ preset               ! 预处理:参数、拓扑、网格、内存、初值
  ├─ solve                ! 主时间推进循环
  └─ env_finalize         ! MPI_Finalize
```

3 段耗时通过 `system_clock` 打点并 `msg_seq_and_master` 输出。

### 1.2 Preset 阶段(`preset.f90`,443 行)

按以下顺序执行:

1. **参数读取与广播**:`input_par`(rank 0 读 `param.inp` 12 个 namelist)→ `broadcast_par`(MPI_PACK + MPI_BCAST)→ `reset_par`(基于 `ndim`、`nvis`、`ntursch` 等派生 `npvt`、`nqvst`、`nd1der_con` 等内部维数)
2. **来流初始化**:`init_inflow`(根据 Ma、Re、α、β 以及 `rinf < 0 / vinf < 0` 触发无量纲化或大气表补全)
3. **拓扑分析**(分支):
   - `ncelcet ≠ 1`:`input_top` → `analyze_top` → `set_bc_seq` / `set_inter_bc` / `set_blk_coms`
   - `ncelcet = 1`(对偶):对应的 `_cc` 变体
4. **网格读取**:`input_grd` / `input_grd_cc`(PLOT3D unformatted stream binary)
5. **网格导数与度量**:`grid_derivative`(按 `ncutpol` 选择精度:3 = 标准、5 = 高阶外推)
6. **奇异搜索**:`search_singulars`(极点、周期 BC 的退化点)
7. **内存分配**:`allocate_fld_variables`(为 `mb_qc`、`mb_pv`、`mb_rhs`、`mb_dt`、…分配每块 3D 数组)
8. **辅助量预计算**:`set_wall_dist`(壁距)、`set_sponge_dist`(吸收层距离,按 `nsponge`)
9. **流场初值**:`init_fld_variables`(`nrestrt=0` 用 `initconst` 常数、`=1` 读 `solfile` 重启、`=2` 含湍流)
10. **内存优化(可选)**:若编译 `MemOPT1`,把度量 `sxyz` 存到磁盘以缩减运行期内存

### 1.3 Solve 阶段(`solve.f90`,6018 行)

```
do nstep = nstepst, nstepmx
    update_cfl()                     ! CFL 自适应(cflst → cfled,cflfac 控制)
    call lhs_dispatch(nlhs)          ! 按 nlhs 选择时间推进器
        nlhs=1  → rkutta_3step
        nlhs=2  → lusgs_std_sca
        nlhs=3  → prsgs_ust_sca
        nlhs=4  → prsgs_std_sca
        nlhs=5  → prsgs_ust_mat
        nlhs=6  → prsgs_std_mat
        nlhs=11 → lusgs_std_prec     (+Preconditioner)
        nlhs=12 → lusgs_std_scmp     (+SCMP)
        nlhs=13 → prsgs_ust_sca_scmp
    内部子迭代 ≤ nsubmax,直至 ||Δq|| < tolsub
    条件 IO:
        mod(nstep, nressav) == 0 → res.plt
        mod(nstep, nfcesav) == 0 → force.plt
        mod(nstep, nsolsav) == 0 → flowfield.sol
        mod(nstep, npltsav) == 0 → tecflow.plt
    acoustic/mean/rms 采样(若 nacous / nprms > 0)
end do
```

每个时间推进器内部 **RHS 评估链** 大致相同:
`halo_exchange → apply_bc → rhs_invscd + rhs_viscous (+ rhs_source) → lhs_update`

---

## 2. 模块架构

### 2.1 按职责分组

**数据类型与全局状态**:

| 文件 | 职责 |
|---|---|
| `mod_kndconsts.f90` | `kind_real`(real64)、`kind_int`(int32)、MPI 类型常量 |
| `mod_constants.f90` | 全部枚举:BC 类型、格式编号、LHS 编号、buffer 类型 |
| `mod_datatypes.f90` | 核心类型:`fld_array_t`、`var_block_t`、`top_block_t`、`bc_region_t`、`bc_inter_t` |
| `mod_variables.f90` | `param.inp` 映射变量 + 运行期派生量(`nstep`、`cfl`、`turconsts[32]` 等) |
| `mod_fieldvars.f90` | 全局流场指针:`mb_qc`、`mb_pv`、`mb_rhs`、`mb_dt`、`mb_qt`、`mb_fmean` 等 |
| `mod_parallels.f90` | `myid`、`numprocs`、`master = 0`、MPI 通信器 |
| `mod_runtimers.f90` | CPU / wall 计时器 |
| `mod_opt_vars.f90` | `MemOPT1` 下的磁盘化度量 buffer |
| `mod_singulars.f90` | 奇异点搜索、平均、标记 |
| `mod_interface.f90` | 变量内存管理的对外接口(create / delete / attach pointer) |

**数值核心**:

| 文件 | 职责 |
|---|---|
| `util_math.f90`(~21k 行) | 所有插值/重构/差分算子(24 种 nintnon + node2-8 + edge2-6),包含 TVD limiter、tridiagonal/pentadiagonal 解,是 `dcsh5pi` 等紧致格式的实际 kernel |
| `metric.f90`(~4.9k 行) | 由 `xyz` 计算雅可比 `sxyz`、面积 `area`、体积 `vol`;按 `ncutpol` 选 edge2-6 / ehen4-8 / ehcs6 / scsl4-6 |
| `eos.f90` | 理想气体:`p = ρRT`、音速 `a = sqrt(γp/ρ)`、Sutherland 粘性 `μ(T)` |
| `rhside.f90` | RHS 总装配器:调度对流 + 粘性 + 源项;处理 BC 外法线 |
| `rhs_invscd.f90` | 对流通量:按 `nflux` 分派到 `flux_steger` / `flux_roe` / `flux_roe_scmp` 等;按 `nintnon` 做左右态重构 |
| `rhs_viscous.f90` | 粘性通量:应力张量 + 热通量散度,支持 `AFLUX_VIS_OPT` |
| `rhs_source.f90` | 源项(重力、科氏力等),当前 `param.inp` 未激活 |

**时间推进**:

| 文件 | 职责 |
|---|---|
| `sol_rkutta.f90` | 3 阶显式 TVD Runge-Kutta(`nlhs=1`) |
| `sol_lusgs.f90` | 标量 LUSGS + 前后向扫描(`nlhs=2`) |
| `sol_lusgs_prec.f90` | LUSGS + 低马数预处理器(`nlhs=11`) |
| `sol_lusgs_scmp.f90` | LUSGS + SCMP(`nlhs=12`) |
| `sol_prmat.f90` | PRSGS 矩阵格式(`nlhs=5/6`) |
| `sol_prsca.f90` | PRSGS 标量格式(`nlhs=3/4`) |
| `sol_prsca_scmp.f90` | PRSGS + SCMP(`nlhs=13`) |

**湍流**:

| 文件 | 职责 |
|---|---|
| `turbulent.f90` | 湍流主驱动 + 指针接驳 |
| `tur_spalart.f90` | Spalart-Allmaras RHS |
| `tur_menter.f90` | Menter SST / HST RHS |
| `tur_lusgs.f90` | 湍流量的 LUSGS 时间推进 |
| `tur_bc.f90` | 湍流量的壁面/远场 BC |

**边界条件与并行通信**:

| 文件 | 职责 |
|---|---|
| `bc.f90` | 全部 BC 实现:wall(3 种子类型)、symmetry、farfield(2 种)、inflow、outflow、pole、patched |
| `communicate.f90`(~4.6k 行) | MPI 块间通信:buffer 管理、`halo_exchange`、`exchange_singulars`、全局 reduce |

**IO**:

| 文件 | 职责 |
|---|---|
| `io_input.f90` | 读 `param.inp`、topology `.inp`、网格 `.grd`、重启 `.sol` |
| `io_output.f90` | 写 `res.plt`(ASCII 残差)、`force.plt`(ASCII 力系数)、`flowfield.sol`(binary 重启) |
| `mod_tecplotio.f90` | 写 Tecplot binary `.plt`(多 zone,每块一 zone) |

**辅助**:

| 文件 | 职责 |
|---|---|
| `util_subs.f90` | 文件操作(open/close/rename/delete) |
| `utilities.f90` | 极值、分量转换、向量操作 |
| `misc.f90` | 消息输出、MPI 同步、参数校验 |
| `initialize.f90` | 常数初值、重启读取、流场构造 |
| `preset.f90` | 见 §1.2 |

---

## 3. 核心数据结构

### 3.1 类型定义(`mod_datatypes.f90`)

```fortran
type :: fld_array_t              ! 3D 字段单元
    integer :: fldtype           ! fldtype_i3d=3 / r3d=-3
    integer, pointer :: i3d(:,:,:)
    real,    pointer :: r3d(:,:,:)
end type

type :: var_block_t              ! 可跨块变量
    character(len=*) :: varname
    type(fld_array_t), pointer :: fld(:)   ! fld(1:nblocks)
end type

type :: top_block_t              ! 拓扑块描述
    integer :: pid, pnb                      ! owner rank, rank 内块号
    integer :: nijk(3)                       ! (ni, nj, nk)
    integer :: pst(3), ped(3)                ! 全局索引区间
    integer :: nregions
    type(bc_region_t), pointer :: bcs(:)     ! BC 列表
end type

type :: bc_region_t              ! 一个 BC 区域
    integer :: bctype                        ! bc_wall / bc_farfield / ...
    integer :: s_st(3), s_ed(3)              ! 源面索引
    integer :: t_st(3), t_ed(3)              ! 目标面索引(块接口)
    integer :: s_nd, s_lr                    ! 面法向、左右
    integer, pointer :: mapijk(:,:,:,:)      ! (i,j,k) -> (ii,jj,kk) 1-to-1 映射
    integer, pointer :: bcts(:,:,:)          ! 子类型(如 wall 的 adiabatic/isothermal)
end type
```

### 3.2 流场数组命名约定(`mod_fieldvars.f90`)

守恒与原始量(`neqn = 5`,含 1D/2D 下相应压缩):
- `mb_qc(nb)%fld(m)%r3d(i,j,k)`,`m ∈ 1..5` 对应 `[ρ, ρu, ρv, ρw, ρE]`
- `mb_pv(nb)%fld(m)%r3d(i,j,k)`,`m ∈ 1..5` 对应 `[ρ, u, v, w, p]`

粘性相关:
- `mb_dpv(nb)%fld(m)%r3d`,`m ∈ 1..12` 对应 `∂{u,v,w,T}/∂{x,y,z}`
- `mb_vsl`、`mb_vst`:层流 / 湍流粘性系数

时间推进相关:
- `mb_rhs(nb)%fld(m)`:RHS 向量(`neqn`)
- `mb_dt(nb)%fld(1)`:每单元时间步
- `mb_dq(nb)%fld(m)`:Δq
- `mb_q0(nb)%fld(m)`:RK/PRSGS 的初始副本
- `mb_aiv(nb)%fld(m)`:PRSGS 矩阵版的 `A⁻¹`(`neqn × neqn`)

湍流相关(`neqt = 1`(SA)或 `2`(SST/HST)):
- `mb_qt`、`mb_dqt`、`mb_rhst`、`mb_dtur`
- `mb_pvt`:混合指针 `[ρ, u, v, w, ν̃]` 等(`npvt = 4 + neqt`)
- `mb_qvst`:混合指针 `[ν̃/k, ω, ν_t]`(`nqvst = neqt + 1`)
- `mb_dst`:壁距(2 分量,到最近壁与到最近锐角)
- `mb_bld`:SST 混合函数

统计:
- `mb_fmean`、`mb_frms`(仅 `nprms > 0`)

### 3.3 内存布局

- 每块是 `(ni, nj, nk)` 的 **Fortran 列主序** 3D 数组
- 每变量分量是一个独立的 `fld_array_t`(Fortran 风格的 "struct of arrays")
- 跨块通过 `mb_xxx(nb)` 索引,`nb` 是 **全局块号**(不是 rank 内编号)

---

## 4. 数值格式菜单

### 4.1 重构/插值(`nintnon`,1–24)

| 编号 | 名称 | 阶 | 类型 | 说明 |
|---|---|---|---|---|
| 1 | `muscl2pv` | 2 | TVD | MUSCL 2 阶,原始变量 |
| 2 | `wcns5pv` | 5 | WCNS | 加权紧致非线性,原始变量 |
| 3 | `wcns5cv` | 5 | WCNS | 加权紧致非线性,守恒变量 |
| 4 | `wcns7pv` | 7 | WCNS | 7 阶 WCNS |
| 5 | `wcns7cv` | 7 | WCNS | 7 阶 WCNS(守恒) |
| 6 | `hdcs5ei` | 5 | Hybrid DCNS | 显式插值 |
| 7 | `hdcs7ci` | 7 | Hybrid DCNS | 紧致插值 |
| 8 | `hdcs5ci` | 5 | Hybrid DCNS | 紧致插值 |
| 9–12 | `scsl3ci` … `scsh5ci` | 3/5 | 子单元光滑/混合 | 紧致 |
| 13–15 | `scsn2ci` … `scsn4ci` | 2/3/4 | 子单元非线性 | 紧致 |
| 16–20 | `scsh3pi` … `scsn4pi` | 3/2/3/4 | 提高阶变体 | 紧致 p-i |
| 21–24 | `dcsh5ci` … `dcsh7pi` | 5/7 | 双紧致 | 含 primitive-interp 与 conservative-interp |

`param.inp` 默认:`nintnon = 23 = dcsh5pi`。

### 4.2 通量(`nflux`,1–7)

| 编号 | 名称 | 说明 |
|---|---|---|
| 1 | `steger` | Steger-Warming 分裂 |
| 2 | `sw_mod` | 修正 Steger-Warming |
| 3 | `vanleer` | Van Leer 分裂 |
| 4 | `roe` | 标准 Roe Riemann |
| 5 | `roe_prec` | Roe + 低马数预处理 |
| 6 | `slau` | SLAU(Shear Layer Adopting Upwind) |
| 7 | `scmp` | Roe + Split Compact Matrix Preconditioner |

`param.inp` 默认:`nflux = 7 = scmp`。

### 4.3 紧致格式策略表(`nscheme`,1–14)

`mod_constants.f90` 的 `nscheme_policys(10, 14)` 把一个 `nscheme` 值映射到 10 个子决策(d1 导数/插值、d2 导数/插值、d3 导数/插值、grid d2/d3 导数/插值,分别针对对流与粘性)。常见取值:

| `nscheme` | 对流端 | 粘性端 |
|---|---|---|
| 1 | edge2 | edge2 |
| 3 | edge4 | edge4 |
| 5 | edge6 | edge6 |
| 7 | edge6 + node8 | ehen8 + node8 |
| 11 | edge6 + node8 | ehcs6 |
| 13 | scsl4 | scsl4 |

### 4.4 限制器(`nlimit`)

`0 = 无`、`1 = MinMod`、`2 = Van Leer`、`3 = Van Albada`。

### 4.5 时间推进(`nlhs`)

见 §1.3 表;`1 = RK3`,`2 = 标量 LUSGS`,`11/12/13 = 各种预处理/SCMP`。

### 4.6 非定常时间步(`nunst`,0–5)

`0 = local`(稳态加速)、`1 = global`(统一)、`2 = local weighted`、`3 = local modified`、`4 = global dual-time`、`5 = local quality`。

---

## 5. 边界条件

### 5.1 主类型(`mod_constants.f90:373-414`)

| `bctype` | 名称 | 实现入口(`bc.f90`) |
|---|---|---|
| -1 | `bc_cut1to1` | 块间 1-to-1 接口(ghost 来自邻块) |
| 2 | `bc_wall` | `set_bc_vis_wall` / `set_bc_inv_wall` |
| 3 | `bc_symmetry` | `set_bc_symmetry_plane` |
| 4 | `bc_farfield` | `set_bc_farfield`(Riemann)或 `_charact` |
| 5 | `bc_inflow` | `set_bc_inflow` |
| 6 | `bc_outflow` | `set_bc_outflow` |
| 7 | `bc_pole` | `set_bc_pole`(轴对称退化) |
| 8 | `bc_patched` | `set_bc_patched`(非结构性块接口) |

### 5.2 子类型

- **Wall**:`adiabatic=1`(∂T/∂n = 0)、`isothermal=2`(T = `twall`)、`slip=3`(无粘)
- **Symmetry**:`point=1`、`plane=2`
- **Farfield**:`riemann=1`、`charact=2`
- **Pole**:`pole_i=71`、`pole_j=72`、`pole_k=73`

---

## 6. 输入/输出

### 6.1 输入

**`param.inp`**(Fortran namelist):见 §7 详列 13 组、~80 字段。

**网格 `.grd`**(PLOT3D unformatted stream):
```
int32       nblocks
repeat nblocks 次:
  int32     ni, nj, nk
  real64    x(ni*nj*nk)
  real64    y(ni*nj*nk)
  real64    z(ni*nj*nk)
```

**拓扑 `.inp`**(GridGen 风格 ASCII):
```
NBLOCKS = N
BLOCK ID = id
BNAME = "name"
IMIN IMAX JMIN JMAX KMIN KMAX = ni_min ni_max nj_min nj_max nk_min nk_max
<bctype> <face_spec>                    ! 每条 BC 区域一行
...
```

**重启 `.sol`**:结构与 `.grd` 类似,字段换为守恒量(可选含湍流量)。

### 6.2 输出

**`res.plt`** — 残差(ASCII):
```
variables = NSTEP, CFL, DT, RESAVE, RESMAX, NB, I, J, K, NV, CPU, WALL
```

**`force.plt`** — 力系数(ASCII),每 `nfcesav` 步:
```
variables = NSTEP, Cfx, Cfy, Cfz, Cmx, Cmy, Cmz, Cd, Cl, Xcp, CPU, WALL
```

`Cd`、`Cl` 由攻角/偏角从 `Cfx/Cfy/Cfz` 转换;动压 `q∞ = (1/2)ρ∞ V∞²`,面积 `frefsc`。

**`flowfield.sol`** — 流场备份(binary),与重启格式相同。

**`tecflow.plt`** — Tecplot binary(`mod_tecplotio.f90`):
```
magic "#!TDV102\0" + version
变量列表:X, Y, Z, Density, U, V, W, Pressure, ...
每块一个 zone:zone_mark(299.0)、name、(ni,nj,nk)、real32 数据块
```

---

## 7. `param.inp` 全字段参考

见 [原程序使用手册 §4](user-manual.md#4-paraminp-字段参考) 的表格形式,本文档不重复。

---

## 8. MPI 并行策略

### 8.1 架构

- rank 0 为 `master`,负责 IO 和全局汇总
- 每块 `top_block_t` 自带 `pid`(owner rank)与 `pnb`(rank 内块号)
- 单 rank 可拥有**多块**(`blkcoms(1..nblkcoms)%nb`)

### 8.2 halo exchange

实现于 `communicate.f90`,`MPI_Sendrecv` / `MPI_Isend`+`MPI_Irecv`:

1. **Pre-exchange**:打包发送 buffer(按 `bc_inter_t` 消息表)
2. **Exchange**:块对块点对点通信
3. **Post-exchange**:解包到接收块的 ghost 区

### 8.3 Buffer 类型(`mod_constants.f90:364-371`)

| 编号 | 名称 | 用途 | 分量数 |
|---|---|---|---|
| 0 | `dyn` | 每步动态分配 | 可变 |
| 1 | `pvs` | 原始量 | 5 |
| 2 | `dqc` | RHS | 5 |
| 3 | `dpv` | 粘性梯度 | 12 |
| 4 | `qts` | 湍流量 | `neqt` |
| 5 | `dqt` | 湍流梯度 | `neqt` |
| 6 | `dtur` | 湍流 3 方向梯度 | `3*neqt` |

### 8.4 全局约简

残差、力、CFL 通过 `MPI_REDUCE` + `MPI_BCAST`。奇异点通过 `exchange_singulars` + `average_singulars` 做值平均(避免退化点多块贡献不一致)。

---

## 9. 构建与运行时特性

### 9.1 CMake 构建目标

| `MAKECMDGOALS` | Compiler | Flags | 产物 |
|---|---|---|---|
| `seq` | gfortran | — | `HOSTA.seq` |
| `mpi` | mpif90 | `-DPARALLEL` | `HOSTA.mpi` |
| `omp` | mpif90 | `-DOMP_IMP -openmp` | `HOSTA.omp` |
| `mpiomp` | mpif90 | `-DPARALLEL -DOMP_IMP -openmp` | `HOSTA.mpiomp` |

### 9.2 编译开关

| 宏 | 作用 |
|---|---|
| `PARALLEL` | 启用 MPI 通信 / broadcast / reduce 分支 |
| `OMP_IMP` | OpenMP 并行 do 循环 |
| `O3` | `-O3` 优化 |
| `logo` | 启动时打印编译器 / 版本 |
| `PARAM_OPT` | 编译期固定 `neqn = 5`(省掉一些动态分配) |
| `MemOPT1` | 把度量系数存磁盘,减少常驻内存 |
| `MemOPT2` | 更激进的内存优化(中间量磁盘化) |
| `AFLUX_VIS_OPT` | 粘性通量的特殊路径优化 |

### 9.3 依赖

- Fortran 编译器:gfortran 或 ifort
- MPI 实现:Linux MPICH / Windows MPICH2
- CMake ≥ 3.2

---

## 10. 与 JAX 移植的对应关系(参考)

作为 JAX 端的设计参考,本文档附上 Fortran↔JAX 模块的 1:1 对照(详见 spec §3):

| Fortran | JAX | 说明 |
|---|---|---|
| `main.f90` | `solver/step.py`(主循环) | JAX 版本无 MPI,主循环是 Python for |
| `preset.f90` | `core/config.py` + `mesh/io_grd.py` + `metric.py` | 分拆为更细粒度模块 |
| `metric.f90` | `mesh/metric.py` | 直接对应 |
| `rhs_invscd.f90` | `rhs/invscd.py` + `flux/*` + `physics/reconstruct/*` | 通量与重构分离 |
| `rhs_viscous.f90` | `rhs/viscous.py` | 直接对应 |
| `bc.f90` | `physics/bc.py` | 只实现 wall/farfield/symmetry/interface |
| `sol_rkutta.f90` | `time/rk.py` | 直接对应 |
| `sol_lusgs*.f90` | **不实现** | JAX 版非目标 |
| `tur_*.f90` | **不实现** | JAX 版非目标 |
| `communicate.f90` | `mesh/halo.py` | 单设备,无 MPI |
| `io_input.f90` | `mesh/io_grd.py` | 保留 PLOT3D 格式读取 |
| `io_output.f90` | `io/tecplot.py` | 仅 Tecplot 输出 |
| `util_math.f90` | `physics/reconstruct/*`(仅需要的) | 阶段 1 只需 MUSCL,阶段 2 需 `dcsh5pi` |
| `mod_datatypes.f90` | `core/types.py` | pytree 化 |
| `mod_variables.f90` | `core/config.py` | dataclass |
| `mod_constants.f90` | `core/constants.py` | Python enum/IntEnum |
