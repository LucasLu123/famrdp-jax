"""网格域管理：读取 PLOT3D 网格和 GridGen 拓扑，生成全局 1D 索引。

数据布局：
    所有块（含 ghost 层）展平为一维全局数组。
    设块 nb 的全局偏移为 offset[nb]，则
        global_idx(nb, i, j, k) = offset[nb] + i * nj_g * nk_g + j * nk_g + k
    其中 ni_g = ni + 2*ghost，nj_g / nk_g 类似。

    "真实单元"：ghost 层以内的单元，即
        i in [ghost, ghost+ni)，j/k 类似。
"""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple, Dict, Optional

import numpy as np
import jax.numpy as jnp

from famrdp_jax.core.constants import BCType, Face
from famrdp_jax.mesh.io_grd import read_grd
from famrdp_jax.mesh.io_top import parse_topology


@dataclass
class BlockInfo:
    """单个块的元数据（纯 Python，不含 JAX 数组）。"""
    block_id: int
    ni: int           # 真实格数（不含 ghost）
    nj: int
    nk: int
    ghost: int
    offset: int       # 在全局 1D 数组中的起始位置

    @property
    def ni_g(self): return self.ni + 2 * self.ghost
    @property
    def nj_g(self): return self.nj + 2 * self.ghost
    @property
    def nk_g(self): return self.nk + 2 * self.ghost
    @property
    def n_total(self): return self.ni_g * self.nj_g * self.nk_g

    def local_idx(self, i: int, j: int, k: int) -> int:
        """块内局部线性索引（含 ghost 层偏移）。"""
        return i * self.nj_g * self.nk_g + j * self.nk_g + k

    def global_idx(self, i: int, j: int, k: int) -> int:
        """全局 1D 数组索引。"""
        return self.offset + self.local_idx(i, j, k)

    def true_cell_indices(self) -> np.ndarray:
        """所有真实单元的全局索引，shape (ni*nj*nk,)。"""
        g = self.ghost
        idxs = []
        for i in range(g, g + self.ni):
            for j in range(g, g + self.nj):
                for k in range(g, g + self.nk):
                    idxs.append(self.global_idx(i, j, k))
        return np.array(idxs, dtype=np.int32)

    def cell_index_3d(self) -> np.ndarray:
        """全局索引的 3D 视图（含 ghost），shape (ni_g, nj_g, nk_g)。"""
        arr = np.arange(self.n_total, dtype=np.int32).reshape(
            self.ni_g, self.nj_g, self.nk_g
        ) + self.offset
        return arr


@dataclass
class DomainData:
    """域管理器输出：全局索引表 + 块元数据 + 拓扑 + 坐标 + 边界信息。"""
    blocks: List[BlockInfo]
    xyz_flat: np.ndarray           # (N_total, 3) 所有单元坐标（含 ghost，ghost 用边界值填充）
    topos: list                    # list[BlockTopology]（来自 famrdp_jax）
    true_idx: np.ndarray           # (N_true,) 所有真实单元的全局索引
    n_total: int                   # 全局数组总长度（含所有块的 ghost 层）
    ghost: int


def build_domain(grd_path, inp_path, *, ghost: int = 2) -> DomainData:
    """读取网格 + 拓扑，构建全局 1D 索引结构。

    Parameters
    ----------
    grd_path : PLOT3D .grd 文件路径
    inp_path : GridGen .inp 拓扑文件路径
    ghost    : ghost 层数（MUSCL-2 需要至少 2）

    Returns
    -------
    DomainData
    """
    xyz_blocks = read_grd(grd_path)    # list of (3, ni, nj, nk) float64
    topos      = parse_topology(inp_path)

    nblocks = len(xyz_blocks)
    assert len(topos) == nblocks

    # 各块网格尺寸（真实格数）
    block_infos: List[BlockInfo] = []
    offset = 0
    for nb in range(nblocks):
        _, ni, nj, nk = xyz_blocks[nb].shape   # (3, ni, nj, nk)
        bi = BlockInfo(block_id=nb, ni=ni, nj=nj, nk=nk, ghost=ghost, offset=offset)
        block_infos.append(bi)
        offset += bi.n_total

    n_total = offset

    # 构建全局坐标数组（(N_total, 3)），ghost 层用边界单元坐标填充
    xyz_flat = np.zeros((n_total, 3), dtype=np.float64)
    for nb, bi in enumerate(block_infos):
        xyz_b = xyz_blocks[nb]   # (3, ni, nj, nk)
        g = ghost
        ni, nj, nk = bi.ni, bi.nj, bi.nk

        # 将真实格坐标写入带 ghost 的缓冲区，再展平
        buf = np.zeros((3, bi.ni_g, bi.nj_g, bi.nk_g), dtype=np.float64)
        buf[:, g:g+ni, g:g+nj, g:g+nk] = xyz_b

        # ghost 层用最近的真实边界层填充（edge padding）
        for ax in [1, 2, 3]:
            sl_lo = [slice(None)] * 4
            sl_hi = [slice(None)] * 4
            sl_lo[ax] = slice(0, g)
            sl_hi[ax] = slice(-g, None)
            ref_lo = [slice(None)] * 4
            ref_hi = [slice(None)] * 4
            ref_lo[ax] = g
            ref_hi[ax] = -(g + 1)
            buf[tuple(sl_lo)] = np.expand_dims(buf[tuple(ref_lo)], axis=ax)
            buf[tuple(sl_hi)] = np.expand_dims(buf[tuple(ref_hi)], axis=ax)

        # (3, ni_g, nj_g, nk_g) → (ni_g*nj_g*nk_g, 3)
        flat = buf.reshape(3, -1).T   # (N_block, 3)
        xyz_flat[bi.offset : bi.offset + bi.n_total] = flat

    # 收集所有真实单元全局索引
    true_idx = np.concatenate([bi.true_cell_indices() for bi in block_infos])

    return DomainData(
        blocks=block_infos,
        xyz_flat=xyz_flat,
        topos=topos,
        true_idx=true_idx,
        n_total=n_total,
        ghost=ghost,
    )


# ---------------------------------------------------------------------------
# CUT1TO1 ghost 索引映射
# ---------------------------------------------------------------------------

def build_cut1to1_maps(domain: DomainData):
    """为所有 CUT1TO1 接口预计算 (src_global_idx, dst_global_idx) 对。

    Returns
    -------
    list of (dst_idx, src_idx): np.ndarray pairs, one per interface pair.
        dst_idx[i] <- prime[src_idx[i]] 执行 ghost 层数据传递。
    """
    ghost = domain.ghost
    maps = []
    seen = set()  # 避免重复处理同一接口

    for bi in domain.blocks:
        nb = bi.block_id
        topo = domain.topos[nb]
        for face, bct in topo.bc_type.items():
            if bct != BCType.CUT1TO1:
                continue
            edge_key = (nb, face)
            if edge_key in seen:
                continue

            nbr_info = topo.neighbors[face]
            nbr_nb, nbr_face, _ = nbr_info
            seen.add(edge_key)
            seen.add((nbr_nb, nbr_face))

            nbr_bi = domain.blocks[nbr_nb]

            # 当前块的 ghost 层索引（需被填充）
            dst_idx = _ghost_layer_indices(bi, face, ghost)
            # 邻块对应的真实单元索引（数据源）
            src_idx = _real_layer_indices(nbr_bi, nbr_face, ghost)

            assert dst_idx.shape == src_idx.shape, (
                f"CUT1TO1 形状不匹配: block {nb} face {face} "
                f"dst {dst_idx.shape} vs src {src_idx.shape}"
            )
            maps.append((dst_idx, src_idx))

    return maps


def _face_axis_range(bi: BlockInfo, face: Face):
    """返回 (axis, is_min_face, ghost_layers_i, true_layers_i)。"""
    g = bi.ghost
    if face in (Face.IMIN, Face.IMAX):
        axis = 0
        size = bi.ni
    elif face in (Face.JMIN, Face.JMAX):
        axis = 1
        size = bi.nj
    else:
        axis = 2
        size = bi.nk
    is_min = face in (Face.IMIN, Face.JMIN, Face.KMIN)
    if is_min:
        ghost_range = list(range(0, g))
        true_range  = list(range(g, 2 * g))     # 最近 ghost 层对应的真实单元
    else:
        total = size + 2 * g
        ghost_range = list(range(total - g, total))
        true_range  = list(range(total - 2 * g, total - g))
    return axis, is_min, ghost_range, true_range


def _ghost_layer_indices(bi: BlockInfo, face: Face, ghost: int) -> np.ndarray:
    """当前块指定面的所有 ghost 层单元全局索引，按 ghost 距离从外到内排列。"""
    axis, is_min, ghost_range, _ = _face_axis_range(bi, face)
    idxs = []
    if axis == 0:
        for ig in ghost_range:
            for j in range(bi.nj_g):
                for k in range(bi.nk_g):
                    idxs.append(bi.global_idx(ig, j, k))
    elif axis == 1:
        for j in ghost_range:
            for i in range(bi.ni_g):
                for k in range(bi.nk_g):
                    idxs.append(bi.global_idx(i, j, k))
    else:
        for k in ghost_range:
            for i in range(bi.ni_g):
                for j in range(bi.nj_g):
                    idxs.append(bi.global_idx(i, j, k))
    return np.array(idxs, dtype=np.int32)


def _real_layer_indices(bi: BlockInfo, face: Face, ghost: int) -> np.ndarray:
    """邻块指定面的对应真实单元全局索引（提供数据给 ghost 层）。

    对于 IMIN 面的 ghost 层（从外到内）：源为邻块 IMAX 面从里到外的真实层。
    """
    axis, is_min, _, true_range = _face_axis_range(bi, face)
    # 镜像：若 dst 是当前块的 ghost，则 src 是邻块对侧面的紧邻真实层
    # true_range 已经是紧邻 ghost 的真实层
    idxs = []
    if axis == 0:
        for ir in true_range:
            for j in range(bi.nj_g):
                for k in range(bi.nk_g):
                    idxs.append(bi.global_idx(ir, j, k))
    elif axis == 1:
        for j in true_range:
            for i in range(bi.ni_g):
                for k in range(bi.nk_g):
                    idxs.append(bi.global_idx(i, j, k))
    else:
        for k in true_range:
            for i in range(bi.ni_g):
                for j in range(bi.nj_g):
                    idxs.append(bi.global_idx(i, j, k))
    return np.array(idxs, dtype=np.int32)
