"""将 .npz 结果导出为 Tecplot 多块结构网格 ASCII .dat 文件。

用法（在 euler1d_jax/ 目录下运行）：
    conda run -n jax python export_tecplot.py [结果文件.npz] [--out 输出.dat]

Tecplot 格式说明
----------------
- ASCII POINT 格式，每行含一个点的全部变量值
- 多个 ZONE 对应多个网格块
- 变量：X  Y  Z  Density  U  V  W  Pressure  Mach  Vmag
- ZONE 尺寸 I/J/K = 块内单元数（与求解器 cell-center 布局一致）
- Tecplot 遍历顺序：I 最快（innermost），K 最慢
"""
from __future__ import annotations
import sys
import json
import pathlib
import argparse

import numpy as np

# ──────────────────────────────────────────────────────────────────────────
# 路径解析：保证以脚本形式或 import 形式均能找到 euler1d_jax 包
# ──────────────────────────────────────────────────────────────────────────
_HERE = pathlib.Path(__file__).parent
_ROOT = _HERE.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from euler1d_jax.visualize import _read_grd_np, _block_true_indices_ordered, _derived


# ──────────────────────────────────────────────────────────────────────────
# 导出核心函数
# ──────────────────────────────────────────────────────────────────────────

def export_tecplot(
    npz_path: pathlib.Path,
    grd_path: pathlib.Path,
    out_dat: pathlib.Path,
    ghost: int = 2,
    gamma: float = 1.4,
) -> pathlib.Path:
    """导出流场为 Tecplot 多块结构网格 ASCII .dat 文件。

    Parameters
    ----------
    npz_path : .npz 结果文件（含 prime: N_true×5，true_idx: N_true）
    grd_path : .grd 网格文件（PLOT3D 二进制，存储 cell-center 坐标）
    out_dat  : 输出 .dat 路径
    ghost    : ghost 层数（须与求解器配置一致）
    gamma    : 比热比（用于计算马赫数）
    """
    # ── 读取结果 ──────────────────────────────────────────────────────────
    data       = np.load(npz_path)
    prime_true = data["prime"]     # (N_true, 5): [rho, u, v, w, p]
    true_idx   = data["true_idx"]  # (N_true,)
    dvars = _derived(prime_true, gamma)

    # ── 读取网格 ──────────────────────────────────────────────────────────
    xyz_blocks = _read_grd_np(grd_path)
    nblocks    = len(xyz_blocks)

    # ── 重建域信息（与 visualize.py 保持一致）────────────────────────────
    g = ghost
    shapes_g: list[tuple[int, int, int]] = []
    offsets:  list[int] = []
    off = 0
    for nb in range(nblocks):
        ni, nj, nk = xyz_blocks[nb].shape[1:]   # grd 存 cell-center，尺寸 = 单元数
        ni_g = ni + 2 * g
        nj_g = nj + 2 * g
        nk_g = nk + 2 * g
        shapes_g.append((ni_g, nj_g, nk_g))
        offsets.append(off)
        off += ni_g * nj_g * nk_g

    block_info = _block_true_indices_ordered(offsets, shapes_g, g)

    # ── 建立 global_idx → prime_true 行号的反查表 ────────────────────────
    idx_to_row: dict[int, int] = {int(true_idx[r]): r for r in range(len(true_idx))}

    # ── 变量定义 ──────────────────────────────────────────────────────────
    var_names  = ["rho", "u", "v", "w", "p", "mach", "vmag"]
    var_labels = ["Density", "U", "V", "W", "Pressure", "Mach", "Vmag"]
    nvars = len(var_names)

    # 预组装 (N_true, nvars) 矩阵，便于批量索引
    all_true_vars = np.stack([dvars[vn] for vn in var_names], axis=1)  # (N_true, nvars)

    # ── 写文件 ────────────────────────────────────────────────────────────
    out_dat.parent.mkdir(parents=True, exist_ok=True)
    with open(out_dat, "w") as f:
        # 文件头
        var_str = " ".join(f'"{lbl}"' for lbl in var_labels)
        f.write('TITLE = "euler1d_jax Flow Field"\n')
        f.write(f'VARIABLES = "X" "Y" "Z" {var_str}\n')

        for nb, (true_idxs_b, ni, nj, nk) in enumerate(block_info):
            n_cells = ni * nj * nk
            xyz_b   = xyz_blocks[nb]                # (3, ni, nj, nk)

            # 坐标：grd 直接存 cell-center 坐标，直接使用
            # 转置为 (ni, nj, nk, 3) 便于后续 concat
            xyz3 = xyz_b.transpose(1, 2, 3, 0)      # (ni, nj, nk, 3)

            # 流场变量：批量查表
            rows_flat = np.fromiter(
                (idx_to_row.get(int(gidx), -1) for gidx in true_idxs_b),
                dtype=np.int64, count=n_cells,
            )
            valid     = rows_flat >= 0
            sol_flat  = np.zeros((n_cells, nvars), dtype=np.float64)
            if valid.any():
                sol_flat[valid] = all_true_vars[rows_flat[valid]]
            sol = sol_flat.reshape(ni, nj, nk, nvars)   # (ni, nj, nk, nvars)

            # 合并坐标 + 流场 → (ni, nj, nk, 3+nvars)
            all_data = np.concatenate([xyz3, sol], axis=-1)

            # Tecplot POINT 顺序：I 最快（i-innermost）
            # data[i,j,k] 以 C-order 存储，转置为 (nk, nj, ni, ...) 后 flatten
            all_tp = all_data.transpose(2, 1, 0, 3).reshape(n_cells, 3 + nvars)

            # ZONE 头
            f.write(f'ZONE T="Block {nb}", I={ni}, J={nj}, K={nk},'
                    f' DATAPACKING=POINT\n')

            # 数据行
            np.savetxt(f, all_tp, fmt="%.6e", delimiter=" ")

    print(f"已导出: {out_dat}  ({nblocks} blocks, "
          f"vars: X Y Z {' '.join(var_labels)})")
    return out_dat


# ──────────────────────────────────────────────────────────────────────────
# CLI 入口
# ──────────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="euler1d_jax 结果 → Tecplot 多块结构网格 .dat 导出",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""示例:
  python export_tecplot.py results/sol_1d_euler.npz
  python export_tecplot.py results/sol_1d_euler.npz --out results/flow.dat
  python export_tecplot.py results/sol_1d_euler.npz --param param.json --out flow.dat
""",
    )
    parser.add_argument(
        "npz", nargs="?",
        default=str(_HERE / "results/sol_1d_euler.npz"),
        help=".npz 结果文件路径（默认: results/sol_1d_euler.npz）",
    )
    parser.add_argument(
        "--param", default=str(_HERE / "param.json"),
        help="param.json 路径（用于自动定位网格文件，默认: param.json）",
    )
    parser.add_argument(
        "--out", default=None,
        help="输出 .dat 路径（默认与 .npz 同名，后缀改为 .dat）",
    )
    args = parser.parse_args()

    npz_path = pathlib.Path(args.npz)
    if not npz_path.exists():
        print(f"错误：找不到结果文件 {npz_path}", file=sys.stderr)
        sys.exit(1)

    param_file = pathlib.Path(args.param)
    if not param_file.exists():
        print(f"错误：找不到 {param_file}", file=sys.stderr)
        sys.exit(1)

    with open(param_file) as pf:
        param = json.load(pf)
    base     = param_file.parent
    grd_path = (base / param["grid_file"]).resolve()
    ghost    = int(param.get("ghost", 2))
    gamma    = float(param.get("gamma", 1.4))

    out_dat = pathlib.Path(args.out) if args.out else npz_path.with_suffix(".dat")

    export_tecplot(npz_path, grd_path, out_dat, ghost=ghost, gamma=gamma)


if __name__ == "__main__":
    main()
