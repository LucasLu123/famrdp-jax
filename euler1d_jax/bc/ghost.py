"""边界条件：对 ghost 层单元赋值。

支持的 BC 类型（与 test3 配置匹配）：
    WALL     → 滑移壁（Euler 方程，无粘）：法向动量反号，其余不变
    SYMMETRY → 对称面：与 WALL 滑移等价
    FARFIELD → 远场：ghost 层直接设为自由流值

所有操作均在原始变量 prime (N_total, 5) = [rho, u, v, w, p] 上进行。

ghost 层的赋值策略（镜像法）：
    对于 IMIN 面（i 轴最小端），ghost 层 g=1..G 对应真实层 2G, 2G-1, ..., G+1
    例如 ghost=2：
        prime[ghost_0] ← f(prime[true_1])
        prime[ghost_1] ← f(prime[true_0])
    其中 f() = 镜像函数（壁面/对称：法向速度取反）
"""
from __future__ import annotations
import numpy as np
import jax.numpy as jnp

from famrdp_jax.core.constants import BCType, Face
from euler1d_jax.mesh.domain import DomainData, BlockInfo


# ---------------------------------------------------------------------------
# 预计算 BC 索引（Python 层，一次性）
# ---------------------------------------------------------------------------

def build_bc_ops(domain: DomainData, q_inf: np.ndarray) -> list:
    """预计算每个 BC 面的操作描述。

    Returns
    -------
    bc_ops : list of dicts，每个 dict 含：
        'type'      : BCType
        'ghost_idx' : (G*M,) int32  需要被写入的 ghost 单元全局索引
        'src_idx'   : (G*M,) int32  镜像来源的真实单元全局索引（WALL/SYM）
                      或 None（FARFIELD）
        'sign'      : (5,) float64  各分量的符号（WALL/SYM 用）
        'q_inf'     : (5,) float64  自由流值（FARFIELD 用）
        'face'      : Face          调试用
    """
    ghost = domain.ghost
    bc_ops = []

    for bi in domain.blocks:
        nb  = bi.block_id
        topo = domain.topos[nb]
        for face, bct in topo.bc_type.items():
            if bct == BCType.CUT1TO1:
                continue  # 由 halo exchange 处理

            ghost_idx, src_idx = _face_ghost_src_indices(bi, face, ghost)

            if bct == BCType.WALL:
                sign = _wall_slip_sign(face)
                bc_ops.append(dict(type=bct, ghost_idx=ghost_idx,
                                   src_idx=src_idx, sign=sign,
                                   q_inf=None, face=face))
            elif bct == BCType.SYMMETRY:
                sign = _wall_slip_sign(face)
                bc_ops.append(dict(type=bct, ghost_idx=ghost_idx,
                                   src_idx=src_idx, sign=sign,
                                   q_inf=None, face=face))
            elif bct == BCType.FARFIELD:
                bc_ops.append(dict(type=bct, ghost_idx=ghost_idx,
                                   src_idx=None, sign=None,
                                   q_inf=np.array(q_inf, dtype=np.float64),
                                   face=face))

    return bc_ops


def _face_ghost_src_indices(bi: BlockInfo, face: Face, ghost: int):
    """返回 ghost 层全局索引 和 对应镜像真实层全局索引。

    ghost 层从外到内：ghost_idx[0] 是最外层。
    src   层从里到外：src_idx[0]   是最近真实层（镜像关系）。
    """
    g = ghost
    if face == Face.IMIN:
        axis = 0
        # ghost: i = 0, 1, ..., g-1（从外到内）
        # src:   i = 2g-1, 2g-2, ..., g（从里到外，镜像）
        ghost_i = list(range(0, g))
        src_i   = list(range(2*g - 1, g - 1, -1))
        def make_idx(i_vals, other_axes_size1, other_axes_size2):
            idxs = []
            for iv in i_vals:
                for j in range(bi.nj_g):
                    for k in range(bi.nk_g):
                        idxs.append(bi.global_idx(iv, j, k))
            return np.array(idxs, dtype=np.int32)
        return make_idx(ghost_i, bi.nj_g, bi.nk_g), make_idx(src_i, bi.nj_g, bi.nk_g)

    elif face == Face.IMAX:
        total = bi.ni_g
        ghost_i = list(range(total - 1, total - 1 - g, -1))
        src_i   = list(range(total - 2*g, total - g))
        def make_idx_i(i_vals):
            idxs = []
            for iv in i_vals:
                for j in range(bi.nj_g):
                    for k in range(bi.nk_g):
                        idxs.append(bi.global_idx(iv, j, k))
            return np.array(idxs, dtype=np.int32)
        return make_idx_i(ghost_i), make_idx_i(src_i)

    elif face == Face.JMIN:
        ghost_j = list(range(0, g))
        src_j   = list(range(2*g - 1, g - 1, -1))
        def make_idx_j(j_vals):
            idxs = []
            for jv in j_vals:
                for i in range(bi.ni_g):
                    for k in range(bi.nk_g):
                        idxs.append(bi.global_idx(i, jv, k))
            return np.array(idxs, dtype=np.int32)
        return make_idx_j(ghost_j), make_idx_j(src_j)

    elif face == Face.JMAX:
        total = bi.nj_g
        ghost_j = list(range(total - 1, total - 1 - g, -1))
        src_j   = list(range(total - 2*g, total - g))
        def make_idx_j2(j_vals):
            idxs = []
            for jv in j_vals:
                for i in range(bi.ni_g):
                    for k in range(bi.nk_g):
                        idxs.append(bi.global_idx(i, jv, k))
            return np.array(idxs, dtype=np.int32)
        return make_idx_j2(ghost_j), make_idx_j2(src_j)

    elif face == Face.KMIN:
        ghost_k = list(range(0, g))
        src_k   = list(range(2*g - 1, g - 1, -1))
        def make_idx_k(k_vals):
            idxs = []
            for kv in k_vals:
                for i in range(bi.ni_g):
                    for j in range(bi.nj_g):
                        idxs.append(bi.global_idx(i, j, kv))
            return np.array(idxs, dtype=np.int32)
        return make_idx_k(ghost_k), make_idx_k(src_k)

    else:  # KMAX
        total = bi.nk_g
        ghost_k = list(range(total - 1, total - 1 - g, -1))
        src_k   = list(range(total - 2*g, total - g))
        def make_idx_k2(k_vals):
            idxs = []
            for kv in k_vals:
                for i in range(bi.ni_g):
                    for j in range(bi.nj_g):
                        idxs.append(bi.global_idx(i, j, kv))
            return np.array(idxs, dtype=np.int32)
        return make_idx_k2(ghost_k), make_idx_k2(src_k)


def _wall_slip_sign(face: Face) -> np.ndarray:
    """滑移壁/对称面：法向速度分量取反，其余保持。"""
    sign = np.ones(5, dtype=np.float64)
    # 原始变量 [rho, u, v, w, p]，速度分量索引 1/2/3
    axis_vel = {Face.IMIN: 1, Face.IMAX: 1,
                Face.JMIN: 2, Face.JMAX: 2,
                Face.KMIN: 3, Face.KMAX: 3}
    sign[axis_vel[face]] = -1.0
    return sign


# ---------------------------------------------------------------------------
# 运行时 BC 施加（纯 JAX）
# ---------------------------------------------------------------------------

def apply_bc_all(prime: jnp.ndarray, bc_ops: list) -> jnp.ndarray:
    """将所有 BC 操作依次施加到 prime。

    Parameters
    ----------
    prime  : (N_total, 5) 当前原始变量（含 ghost 层）
    bc_ops : build_bc_ops() 返回的操作列表

    Returns
    -------
    prime  : 施加 BC 后的数组（ghost 层已更新）
    """
    for op in bc_ops:
        bct       = op['type']
        ghost_idx = jnp.array(op['ghost_idx'])

        if bct in (BCType.WALL, BCType.SYMMETRY):
            src_idx = jnp.array(op['src_idx'])
            sign    = jnp.array(op['sign'])        # (5,)
            mirrored = prime[src_idx] * sign[None, :]  # (G*M, 5)
            prime = prime.at[ghost_idx].set(mirrored)

        elif bct == BCType.FARFIELD:
            q_inf = jnp.array(op['q_inf'])         # (5,)
            n     = ghost_idx.shape[0]
            prime = prime.at[ghost_idx].set(
                jnp.broadcast_to(q_inf[None, :], (n, 5))
            )

    return prime


# ---------------------------------------------------------------------------
# CUT1TO1 Halo exchange
# ---------------------------------------------------------------------------

def apply_halo(prime: jnp.ndarray, cut_maps: list) -> jnp.ndarray:
    """CUT1TO1 块接口 ghost 层数据传递。

    Parameters
    ----------
    prime    : (N_total, 5)
    cut_maps : build_cut1to1_maps() 返回的列表

    Returns
    -------
    prime 更新后的结果
    """
    for dst_idx, src_idx in cut_maps:
        dst = jnp.array(dst_idx)
        src = jnp.array(src_idx)
        prime = prime.at[dst].set(prime[src])
    return prime
