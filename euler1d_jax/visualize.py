"""流场可视化：读取 results/sol_1d_euler.npz，绘制柱绕流 2D 云图。

用法（在 euler1d_jax/ 目录下运行）：
    conda run -n jax python visualize.py [结果文件.npz] [--var rho|p|vmag|mach]

默认绘制 4 个子图：密度 / 压力 / 速度幅值 / 马赫数。
网格坐标从 param.json 所指定的 .grd 文件中读取。
"""
from __future__ import annotations
import sys
import json
import pathlib
import argparse

import numpy as np
import matplotlib
matplotlib.use("Agg")          # 无 GUI 时先用 Agg，有显示再改
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Circle

# ──────────────────────────────────────────────────────────────────────────
# 路径解析
# ──────────────────────────────────────────────────────────────────────────
_HERE = pathlib.Path(__file__).parent
_ROOT = _HERE.parent


def _load_param(param_file: pathlib.Path) -> dict:
    with open(param_file) as f:
        return json.load(f)


# ──────────────────────────────────────────────────────────────────────────
# 域结构重建（不依赖 JAX）
# ──────────────────────────────────────────────────────────────────────────

def _read_grd_np(path) -> list[np.ndarray]:
    """PLOT3D 二进制网格读取（纯 numpy，不需要 JAX）。"""
    with open(path, "rb") as f:
        nblocks = np.frombuffer(f.read(4), dtype="<i4")[0]
        dims = []
        for _ in range(nblocks):
            ni, nj, nk = np.frombuffer(f.read(12), dtype="<i4")
            dims.append((int(ni), int(nj), int(nk)))
        blocks = []
        for ni, nj, nk in dims:
            n = ni * nj * nk
            xyz = np.frombuffer(f.read(n * 3 * 8), dtype="<f8").reshape(3, ni, nj, nk, order="F")
            blocks.append(xyz)
    return blocks


def _block_true_indices_ordered(offsets, shapes_g, ghost):
    """返回每个块真实单元在全局 1D 数组中的索引（行优先，与 true_cell_indices 一致）。"""
    g = ghost
    result = []
    for nb, (ni_g, nj_g, nk_g) in enumerate(shapes_g):
        off = offsets[nb]
        ni = ni_g - 2 * g
        nj = nj_g - 2 * g
        nk = nk_g - 2 * g
        idxs = []
        for i in range(g, g + ni):
            for j in range(g, g + nj):
                for k in range(g, g + nk):
                    idxs.append(off + i * nj_g * nk_g + j * nk_g + k)
        result.append((np.array(idxs, dtype=np.int64), ni, nj, nk))
    return result


# ──────────────────────────────────────────────────────────────────────────
# 流场变量计算
# ──────────────────────────────────────────────────────────────────────────

def _derived(prime_true: np.ndarray, gamma: float = 1.4):
    """计算衍生量。prime_true: (N_true, 5) = [rho, u, v, w, p]"""
    rho = prime_true[:, 0]
    u   = prime_true[:, 1]
    v   = prime_true[:, 2]
    w   = prime_true[:, 3]
    p   = prime_true[:, 4]
    vmag = np.sqrt(u**2 + v**2 + w**2)
    c    = np.sqrt(np.maximum(gamma * p / np.maximum(rho, 1e-30), 1e-30))
    mach = vmag / c
    return dict(rho=rho, u=u, v=v, w=w, p=p, vmag=vmag, mach=mach)


# ──────────────────────────────────────────────────────────────────────────
# 绘图
# ──────────────────────────────────────────────────────────────────────────

_VAR_META = {
    "rho":  dict(label=r"Density $\rho$",          cmap="RdBu_r",  symlog=False),
    "p":    dict(label=r"Pressure $p$",             cmap="plasma",  symlog=False),
    "vmag": dict(label=r"Velocity magnitude $|v|$", cmap="viridis", symlog=False),
    "mach": dict(label=r"Mach number",              cmap="hot_r",   symlog=False),
}


def plot_field(npz_path: pathlib.Path, grd_path: pathlib.Path,
               ghost: int = 2, gamma: float = 1.4,
               vars_to_plot: list | None = None,
               out_png: pathlib.Path | None = None,
               xlim: tuple | None = None,
               ylim: tuple | None = None):
    """绘制流场云图。

    Parameters
    ----------
    npz_path   : .npz 结果文件
    grd_path   : .grd 网格坐标文件
    ghost      : ghost 层数
    gamma      : 比热比
    vars_to_plot : 要绘制的变量列表，默认 ['rho', 'p', 'vmag', 'mach']
    out_png    : 输出 PNG 路径（None → 弹窗显示）
    """
    if vars_to_plot is None:
        vars_to_plot = ["rho", "p", "vmag", "mach"]

    # ── 读取结果 ──
    data = np.load(npz_path)
    prime_true = data["prime"]          # (N_true, 5)
    true_idx   = data["true_idx"]       # (N_true,)
    dvars = _derived(prime_true, gamma)

    # ── 读取网格 ──
    xyz_blocks = _read_grd_np(grd_path)
    nblocks = len(xyz_blocks)

    # ── 重建域信息 ──
    g = ghost
    shapes_g = []
    offsets  = []
    off = 0
    for nb in range(nblocks):
        ni, nj, nk = xyz_blocks[nb].shape[1:]
        ni_g, nj_g, nk_g = ni + 2*g, nj + 2*g, nk + 2*g
        shapes_g.append((ni_g, nj_g, nk_g))
        offsets.append(off)
        off += ni_g * nj_g * nk_g

    block_info = _block_true_indices_ordered(offsets, shapes_g, g)

    # ── 建立真实单元 global_idx → local 位置的反查表 ──
    # prime_true[r] 对应 true_idx[r]；需要知道每个 true_idx 属于哪个块、哪个位置
    idx_to_row = {int(true_idx[r]): r for r in range(len(true_idx))}

    # ── 选取 k 中间切面 ──
    n_vars = len(vars_to_plot)
    fig, axes = plt.subplots(1, n_vars, figsize=(5 * n_vars, 5),
                              constrained_layout=True)
    if n_vars == 1:
        axes = [axes]

    for ax, vname in zip(axes, vars_to_plot):
        var_data = dvars[vname]   # (N_true,)
        meta = _VAR_META.get(vname, dict(label=vname, cmap="viridis", symlog=False))

        # 遍历每个块，取 k 中间切面，先收集全局 vmin/vmax
        all_x, all_y, all_v = [], [], []
        for nb, (true_idxs_b, ni, nj, nk) in enumerate(block_info):
            k_mid = nk // 2
            xyz_b = xyz_blocks[nb]   # (3, ni, nj, nk)
            x2d = xyz_b[0, :, :, k_mid]
            y2d = xyz_b[1, :, :, k_mid]
            val2d = np.zeros((ni, nj), dtype=np.float64)
            for i in range(ni):
                for j in range(nj):
                    row_in_block = i * nj * nk + j * nk + k_mid
                    global_idx = int(true_idxs_b[row_in_block])
                    local_row  = idx_to_row.get(global_idx, -1)
                    if local_row >= 0:
                        val2d[i, j] = var_data[local_row]
            all_x.append(x2d); all_y.append(y2d); all_v.append(val2d)

        vmin = min(v.min() for v in all_v)
        vmax = max(v.max() for v in all_v)
        # 给 colormap 留一点余量，避免完全均匀时失效
        if abs(vmax - vmin) < 1e-10:
            vmin -= 0.01 * abs(vmin + 1)
            vmax += 0.01 * abs(vmax + 1)
        norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

        pc_last = None
        for x2d, y2d, val2d in zip(all_x, all_y, all_v):
            pc_last = ax.pcolormesh(x2d, y2d, val2d,
                                    cmap=meta["cmap"],
                                    norm=norm,
                                    shading="gouraud")

        if pc_last is not None:
            plt.colorbar(pc_last, ax=ax, shrink=0.85, format="%.3g")

        if xlim is not None:
            ax.set_xlim(xlim)
        if ylim is not None:
            ax.set_ylim(ylim)
        ax.set_aspect("equal")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_title(meta["label"])

    fig.suptitle(f"euler1d_jax Flow Field  —  {npz_path.name}", fontsize=13)

    if out_png is None:
        # 尝试弹窗；若无显示则改存文件
        try:
            matplotlib.use("TkAgg")
            import importlib
            importlib.reload(matplotlib.pyplot)
            plt.show()
        except Exception:
            out_png = npz_path.with_suffix(".png")
            print(f"无 GUI，保存至 {out_png}")
            fig.savefig(out_png, dpi=150, bbox_inches="tight")
    else:
        fig.savefig(out_png, dpi=150, bbox_inches="tight")
        print(f"已保存: {out_png}")

    plt.close(fig)
    return out_png


# ──────────────────────────────────────────────────────────────────────────
# CLI 入口
# ──────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="euler1d_jax 流场可视化")
    parser.add_argument("npz", nargs="?",
                        default=str(_HERE / "results/sol_1d_euler.npz"),
                        help=".npz 结果文件路径")
    parser.add_argument("--param", default=str(_HERE / "param.json"),
                        help="param.json 路径（用于定位网格文件）")
    parser.add_argument("--var", nargs="+",
                        default=["rho", "p", "vmag", "mach"],
                        choices=list(_VAR_META.keys()),
                        help="要绘制的变量")
    parser.add_argument("--out", default=None,
                        help="输出 PNG 路径（默认与 .npz 同名）")
    parser.add_argument("--xlim", nargs=2, type=float, default=None,
                        metavar=("XMIN", "XMAX"),
                        help="x 轴范围，例如 --xlim -5 5（默认全域）")
    parser.add_argument("--ylim", nargs=2, type=float, default=None,
                        metavar=("YMIN", "YMAX"),
                        help="y 轴范围，例如 --ylim -5 5")
    parser.add_argument("--zoom", action="store_true",
                        help="自动缩放到圆柱附近 (±5 倍圆柱半径范围)")
    args = parser.parse_args()

    npz_path   = pathlib.Path(args.npz)
    param_file = pathlib.Path(args.param)

    if not npz_path.exists():
        print(f"错误：找不到结果文件 {npz_path}", file=sys.stderr)
        sys.exit(1)

    # 从 param.json 找网格文件
    param = _load_param(param_file)
    base  = param_file.parent
    grd_path = (base / param["grid_file"]).resolve()
    ghost    = int(param.get("ghost", 2))
    gamma    = float(param.get("gamma", 1.4))

    # 估算圆柱半径（网格最小半径）
    xlim = tuple(args.xlim) if args.xlim else None
    ylim = tuple(args.ylim) if args.ylim else None
    if args.zoom and xlim is None:
        xyz_blocks = _read_grd_np(grd_path)
        r_min = min(np.sqrt(xyz_blocks[0][0]**2 + xyz_blocks[0][1]**2).min(),
                    np.sqrt(xyz_blocks[1][0]**2 + xyz_blocks[1][1]**2).min())
        r_zoom = max(5.0 * r_min, 5.0)   # 至少 ±5
        xlim = (-r_zoom, r_zoom)
        ylim = (-r_zoom, r_zoom)
        print(f"自动缩放: r_cylinder≈{r_min:.2f}, 显示范围 ±{r_zoom:.1f}")

    out_png = pathlib.Path(args.out) if args.out else None
    if out_png is None:
        out_png = npz_path.with_suffix(".png")

    plot_field(npz_path, grd_path,
               ghost=ghost, gamma=gamma,
               vars_to_plot=args.var,
               out_png=out_png,
               xlim=xlim, ylim=ylim)


if __name__ == "__main__":
    main()
