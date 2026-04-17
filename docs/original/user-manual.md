# HOSTA / FAMRDP 使用手册

- **对象**:原 Fortran 程序(`HOSTA.mpi` / `.seq` / `.omp` / `.mpiomp`)
- **范围**:环境搭建、构建、输入配置、运行、输出解读、常见问题
- **日期**:2026-04-17

---

## 1. 系统要求

### 1.1 目标平台

- **主支持**:Ubuntu 22.04(CMake 脚本 Linux 分支最完整)
- **备选**:其他 Linux(需手动调整 MPI 路径)、Windows(MPICH2 路径硬编码在 `CMakeLists.txt`,需自行确认)

### 1.2 编译器与依赖

| 依赖 | 版本 | 备注 |
|---|---|---|
| gfortran | ≥ 9 | 默认编译器;可替换为 ifort(修改 `CMAKE_Fortran_COMPILER`) |
| MPICH / MPICH2 | ≥ 3 | Linux 装 mpich;Windows 装 MPICH2 到 `C:/Program Files/MPICH2` |
| CMake | ≥ 3.2 | |
| OpenMP | — | 仅 `omp` / `mpiomp` 目标需要 |

### 1.3 一键安装工具链(Ubuntu)

```bash
sudo apt update
sudo apt install cmake gfortran mpich
```

---

## 2. 获取与构建

### 2.1 克隆

```bash
git clone http://10.136.126.255/xxx/famrdp.git
cd famrdp
```

### 2.2 构建目标选择

`CMakeLists.txt` 通过 `MAKECMDGOALS` 变量切换构建目标。当前脚本默认 `mpi`。要显式指定:

```bash
# 默认 mpi 目标(最常用)
cmake . -B build
cmake --build build -j 16

# 其他目标:修改 CMakeLists.txt 的 SET(MAKECMDGOALS xxx)
#   seq     — 单进程,无 MPI
#   mpi     — MPI 并行 (推荐)
#   omp     — OpenMP
#   mpiomp  — MPI + OpenMP
```

产物位于 `build/src/`,可执行文件名按目标后缀:

| 目标 | 可执行 |
|---|---|
| seq | `HOSTA.seq` |
| mpi | `HOSTA.mpi` |
| omp | `HOSTA.omp` |
| mpiomp | `HOSTA.mpiomp` |

### 2.3 可选编译开关

在 `CMakeLists.txt` 开启额外宏(`SET_OPT_FLAGS*`)时生效:

| 宏 | 何时启用 | 代价/收益 |
|---|---|---|
| `PARAM_OPT` | 固定 5 方程模型 | 稍快,但失去灵活性 |
| `MemOPT1` | 大算例、内存紧张 | 度量系数存磁盘,IO 略增 |
| `MemOPT2` | 极大算例 | 更激进的磁盘化 |
| `AFLUX_VIS_OPT` | 对流 + 粘性通量合并路径 | 某些格式下提速 |

---

## 3. 运行

### 3.1 工作目录准备

可执行读取**当前目录**下的 `param.inp` 与 `param.inp` 中指定的网格/拓扑文件。所以必须拷贝输入文件到可执行同目录:

```bash
cp param.inp test3.grd test3.inp build/src/
cd build/src
```

### 3.2 启动(MPI 目标)

```bash
# 单进程(调试)
./HOSTA.mpi

# 多进程
mpirun -np 4 ./HOSTA.mpi
```

### 3.3 启动(seq 目标)

```bash
./HOSTA.seq
```

### 3.4 启动(OpenMP 目标)

```bash
export OMP_NUM_THREADS=8
./HOSTA.omp
```

### 3.5 进程数建议

- 进程数 ≤ 块数(每块只能属一个 rank)
- 若 `numprocs > nblocks`,多余 rank 会空转;反之 rank 持多块
- 推荐从 `np = nblocks` 开始调优

---

## 4. `param.inp` 字段参考

Fortran namelist 格式,`&组名` 开头,`/` 结尾。13 个组。

### 4.1 `&general` — 案例基础

| 字段 | 含义 | 示例 |
|---|---|---|
| `casname` | 案例名(会进入输出文件名) | `"g01"` |
| `indir` | 输入目录(空串即当前目录) | `" "` |
| `outdir` | 输出目录 | `" "` |

### 4.2 `&inflow` — 来流条件

| 字段 | 含义 | 单位 | 说明 |
|---|---|---|---|
| `moo` | 马赫数 | — | 主控来流速度 |
| `reno` | 雷诺数 | — | 基于 `reflen` |
| `alt` | 高度 | m | 负值 = 海平面或手动 |
| `attack` | 攻角 | 度 | |
| `ayaw` | 偏角 | 度 | |
| `rinf` | 来流密度 | kg/m³ | **负值** = 按 `alt` 自动 |
| `vinf` | 来流速度 | m/s | **负值** = 按 `moo` 自动 |
| `pinf` | 来流静压 | Pa | |
| `tinf` | 来流温度 | K | |

### 4.3 `&reference` — 无量纲/力参考

| 字段 | 含义 | 单位 |
|---|---|---|
| `reflen` | 无量纲化参考长度 | m |
| `reflgrd` | 网格参考长度(计算用) | m |
| `freflen` | 力矩参考长度 | m |
| `frefsc` | 力参考面积 | m² |
| `frefxc/yc/zc` | 力矩参考点坐标 | m |

### 4.4 `&filename` — 文件路径与格式

| 字段 | 含义 |
|---|---|
| `topfile` / `topsty` | topology 文件路径 / 格式(0 = GridGen) |
| `grdfile` / `grdsty` | 网格文件 / 格式(0 = PLOT3D binary) |
| `resfile` | 残差输出(默认 `res.plt`) |
| `fcefile` | 力输出(默认 `force.plt`) |
| `solfile` | 流场重启/备份(默认 `flowfield.sol`) |
| `pltfile` / `pltsty` | Tecplot 输出 / 格式(0 = standard、1 = function、5 = 高阶导数) |

### 4.5 `&initconst` — 初值

| 字段 | 含义 |
|---|---|
| `nincst` | 模式:0 = uniform、1 = Ma 分布、2 = Ma + BC 混合 |
| `mincst` | 初始 Ma |
| `alpcst` / `betcst` | 初始攻 / 偏角(度) |
| `pincst` / `tincst` | 初始压力(Pa)/ 温度(K) |

### 4.6 `&control` — 迭代控制

| 字段 | 含义 |
|---|---|
| `nrestrt` | 0 = 新开始、1 = 读流场重启、2 = 重启 + 湍流 |
| `nstepmx` | 最大迭代步 |
| `nressav` / `nfcesav` / `nsolsav` / `npltsav` | 输出频率(步) |

### 4.7 `&physical` — 物理模型

| 字段 | 含义 |
|---|---|
| `ndim` | 1 = 轴对称、2 = 2D、3 = 3D |
| `nvis` | 0 = Euler、1 = NS 层流、2 = SA、3 = SST、4 = HST |
| `twall` | 壁温(K);**负值** = 绝热 |

### 4.8 `&turbulent` — 湍流求解器

| 字段 | 含义 |
|---|---|
| `nturlhs` | 1 = LUSGS、2 = PRSGS |
| `ntursub` | 湍流子迭代次数 |
| `ntursch` | 湍流格式(1–10) |
| `nturint` | 湍流积分方法 |
| `nrddst` | 0 = 计算壁距、1 = 读取 |
| `ntke2s` | k-ε → SST 转换开关 |
| `nearsm` | 近壁修正开关 |
| `kfstur` / `kmaxtur` | k 初值 / 上限 |
| `vfstur` / `vmaxtur` | ν_t 初值 / 上限 |

### 4.9 `&timestep` — 时间推进

| 字段 | 含义 | 关键取值 |
|---|---|---|
| `nprec` | 预处理开关 | 0 = 无、> 0 = 启用 |
| `nscmp` | SCMP 开关 | 0 = 无、> 0 = 启用 |
| `nlhs` | 时间推进类型 | 1 = RK3,2 = LUSGS 标量,3/4 = PRSGS scalar,5/6 = PRSGS matrix,11 = LUSGS+Prec,12 = LUSGS+SCMP,13 = PRSGS+SCMP |
| `nsubmax` | 最大子迭代数 | |
| `tolsub` | 子迭代收敛容差(相对) | 如 0.01 |
| `nunst` | 非定常时间步模式 | 0 = local(稳态)、1 = global、4 = dual-time |
| `cflst` / `cfled` | CFL 初值 / 终值 | |
| `cflfac` | CFL 渐进系数 | 如 0.001 |
| `cdtmax` | 最大物理时间步 | 负值 = 由 CFL 决定 |
| `relaxs` / `relaxp` | 密度 / 压力松弛(PRSGS 用) | |
| `dtau` / `dtime` | 虚时间步 / 物理时间步 | |

### 4.10 `&method` — 数值方法

| 字段 | 含义 | 关键取值 |
|---|---|---|
| `nscheme` | 紧致格式策略 | 1–14,见功能说明 §4.3 |
| `nflux` | 通量 | 1 = Steger、4 = Roe、6 = SLAU、7 = SCMP(默认) |
| `nintnon` | 非守恒插值 | 1 = MUSCL2、23 = DCSH5pi(默认)、其余见 §4.1 |
| `nlimit` | 限制器 | 0 / 1 = MinMod / 2 = Van Leer / 3 = Van Albada |
| `enfix` | 熵修复系数 | |
| `ckmuscl` / `cbmuscl` | MUSCL κ / β 参数 | |
| `csrvis` | 粘性源修正 | |
| `csafeup` | 上风安全系数 | |
| `rlimit` | 密度 clip 区间 | `[min, max]` |
| `plimit` | 压力 clip 区间 | `[min, max]` |
| `tlimit` | 温度 clip 区间 | `[min, max]` |

### 4.11 `&technic` — 技术选项

| 字段 | 含义 |
|---|---|
| `ncutpol` | 网格导数阶数(3 标准 / 5 高阶) |
| `ncelcet` | 0 = 标准结构块、1 = 对偶网格 |
| `nghnode` | ghost 节点层数(含计算区) |
| `nghedge` | 边 ghost 层数 |
| `nsponge` | 吸收层开关(负值 = 关) |
| `sigspg` | 吸收强度 |
| `nrddsp` | 吸收距离读取开关 |
| `nacous` | 声学采样开关(> 0 启用) |
| `nprms` | 均值 / RMS 采样开关(> 0 启用) |
| `ns_mean` / `ns_prms` | 采样起始 / 结束步 |

### 4.12 `&gasmodel` — 气体属性

| 字段 | 含义 | 空气默认 |
|---|---|---|
| `wgas` | 摩尔质量(kg/kmol) | 28.97 |
| `gamma` | 比热比 | 1.4 |
| `mu0sth` | Sutherland μ₀(Pa·s) | 1.71608e-5 |
| `t0sth` | Sutherland T₀(K) | 273.15 |
| `tssth` | Sutherland Tₛ(K) | 110.4 |
| `prlam` / `prtur` | 层流 / 湍流 Prandtl 数 | 0.72 / 0.90 |
| `sclam` / `sctur` | 层流 / 湍流 Schmidt 数 | 0.50 / 0.50 |

### 4.13 `&monitor` — 监测点

| 字段 | 含义 |
|---|---|
| `nmoni` | 监测点数(负值 = 关) |
| `nbmoni[30]` | 每点所在块号 |
| `nimoni[30]` / `njmoni[30]` / `nkmoni[30]` | 每点 (i, j, k) 索引 |

---

## 5. 网格与拓扑文件

### 5.1 网格 `.grd`(PLOT3D binary stream)

**二进制结构**(小端 / 本地):
```
int32   nblocks
<block 1>
<block 2>
...

每 block:
int32   ni  nj  nk
real64  x(ni*nj*nk)
real64  y(ni*nj*nk)
real64  z(ni*nj*nk)
```

坐标数组存储顺序:`(i, j, k)` Fortran 列主序。

### 5.2 拓扑 `.inp`(GridGen ASCII)

```
NBLOCKS = 2
BLOCK ID = 1
BNAME = "block1"
IMIN IMAX JMIN JMAX KMIN KMAX = 1 65 1 33 1 1
4 IMAX JMIN JMAX KMIN KMAX = 65 65 1 33 1 1
  ! farfield at IMAX 面,范围 i=65 / j=1..33 / k=1
3 IMIN IMAX JMIN JMAX KMIN KMAX = 1 65 1 1 1 1
  ! symmetry at JMIN 面
...

BLOCK ID = 2
BNAME = "block2"
IMIN IMAX JMIN JMAX KMIN KMAX = 1 129 1 97 1 1
-1 IMIN JMIN JMAX KMIN KMAX = 1 1 1 97 1 1
  ! cut1to1 to block 1 at some face
...
```

BC 类型编码:见 功能说明 §5.1。

---

## 6. 输出文件解读

### 6.1 `res.plt`(残差,ASCII)

每 `nressav` 步一行:

```
variables = "NSTEP","CFL","DT","RESAVE","RESMAX","NB","I","J","K","NV","CPU Time","Wall Time"
       1   1.00e+02  1.00e-02  5.00e-03  1.00e-02   1  32  50   1   2  0.50  0.60
       2   ...
```

- `RESAVE`:全域 L2 归一化残差
- `RESMAX` + `NB/I/J/K/NV`:最大残差位置(块号、索引、哪个方程)
- 用于收敛判断(`RESAVE` 下降若干数量级 = 收敛)

### 6.2 `force.plt`(力系数,ASCII)

每 `nfcesav` 步:

```
variables = "NSTEP","Cfx","Cfy","Cfz","Cmx","Cmy","Cmz","Cd","Cl","Xcp","CPU Time","Wall Time"
```

- `Cfx/y/z` — 笛卡尔分量力系数
- `Cmx/y/z` — 力矩系数,参考点 `(frefxc, frefyc, frefzc)`
- `Cd` / `Cl` — 由 attack/yaw 旋转得到的阻力 / 升力系数
- `Xcp` — 压心坐标

### 6.3 `flowfield.sol`(流场备份,binary)

格式与 `.grd` 平行,变量换为守恒量:`(ρ, ρu, ρv, ρw, ρE)`。若启用湍流,附加 `(ν̃)` 或 `(k, ω)`。

用于 `nrestrt ≥ 1` 重启。

### 6.4 `tecflow.plt`(Tecplot binary)

Tecplot 二进制,可用 Tecplot / ParaView(含 Tecplot 插件)/ VisIt 打开。每块作为独立 zone。变量数由 `pltsty` 决定:

- `pltsty = 0`:X, Y, Z, ρ, u, v, w, p
- `pltsty = 1`:function 格式(简化)
- `pltsty = 5`:额外含速度梯度等后处理量(30 分量)

---

## 7. 典型工作流

### 7.1 稳态计算(本文档 `param.inp` 默认配置)

```
&timestep  nunst=0  nlhs=13  cflst=100.0  cfled=100.0  /   ! local time step + PRSGS+SCMP
&method    nscheme=7  nflux=7  nintnon=23  /                ! dcsh5pi + Roe-scmp
&control   nressav=1  nfcesav=10  nsolsav=1e9  npltsav=100 /
```

收敛观察:`tail -f res.plt`,等 `RESAVE` 下降 4–6 数量级。

### 7.2 重启

1. 上轮终止 → 得到 `flowfield.sol`
2. 修改 `&control` 的 `nrestrt = 1`
3. 重新启动

### 7.3 改变来流重新计算

1. 修改 `&inflow` 的 `moo / attack` 等
2. 可选:`nrestrt = 1`(用上次结果作为初值加速收敛)
3. 重启,观察力系数收敛

### 7.4 可视化流场

1. 等 `tecflow.plt` 生成(每 `npltsav` 步)
2. 用 Tecplot 打开:`tecplot tecflow.plt`
3. 或转 VTU:用第三方工具(Plot3D2VTK 等)

---

## 8. 性能调优

### 8.1 进程与块

- 保证 `numprocs ≤ nblocks`;最优 `numprocs = nblocks`
- 若块间 load 严重不平衡,考虑用多 rank 共享块的重新切分(当前代码不自动做)

### 8.2 CFL 策略

- 稳态 PRSGS:`cflst = cfled = 100`,通常能稳定
- 显式 RK3(`nlhs = 1`):CFL 上限取决于格式,一般 ≤ 1.0
- 加速收敛:用 `cflst < cfled`,让 CFL 渐涨

### 8.3 内存

- 大算例使用 `MemOPT1` 或 `MemOPT2` 宏
- 减少 `pltsty` 到 0(默认 30 分量会占较多内存)

### 8.4 编译器优化

- gfortran:`-O3 -march=native -funroll-loops`
- ifort(若可用):`-O3 -xHOST -no-prec-div`

---

## 9. 常见问题

**Q1:启动即 segfault / 内存越界。**
A:通常是 `nghnode / nghedge` 与网格不匹配,或格式 `nintnon` 所需 stencil 超过 ghost 层。检查 `param.inp` 的 `ncutpol`、`nghnode`、`nintnon` 组合。

**Q2:残差不下降 / 发散。**
A:尝试 `cflst` 降到 1.0、`nlimit` 开 MinMod(`1`)、`enfix` 非零(如 0.1 为 Roe 加熵修复)。检查 BC 设定是否正确(尤其 farfield 朝向)。

**Q3:力系数振荡不收敛。**
A:通常是不稳定边界层 / 非定常现象。若是定常问题,检查:
- 来流方向 / 攻角是否合理
- `twall` 设置(绝热 = 负值;等温 = 正 K)
- 网格 y⁺(壁面第一层距离)

**Q4:MPI 启动后 hang。**
A:检查 `topfile` 拓扑文件的 `cut1to1` 接口是否对称配对(A→B 和 B→A 必须一致),以及总块数与 rank 数的关系。

**Q5:输出的 Tecplot 文件打不开。**
A:`pltsty` 版本不匹配。尝试 `pltsty = 0` 并重新计算。或用 `preplot` 工具转换格式。

**Q6:Windows 下编译失败找不到 `fmpich2.lib`。**
A:确认 MPICH2 装在 `C:/Program Files/MPICH2/`;否则修改 `CMakeLists.txt` 的 `MPI_PATH`。

---

## 10. 命令速查

```bash
# 构建
cmake . -B build && cmake --build build -j

# 运行
cp param.inp test3.grd test3.inp build/src/
cd build/src
mpirun -np 4 ./HOSTA.mpi

# 监控
tail -f res.plt        # 残差
tail -f force.plt      # 力系数

# 可视化
tecplot tecflow.plt    # 或 ParaView
```

---

## 11. 参考

- 功能说明:[`docs/original/functional-description.md`](functional-description.md)
- JAX 移植需求:[`docs/requirements.md`](../requirements.md)
- JAX 移植设计 spec:[`docs/superpowers/specs/2026-04-17-famrdp-jax-differentiable-design.md`](../superpowers/specs/2026-04-17-famrdp-jax-differentiable-design.md)
