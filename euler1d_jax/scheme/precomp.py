"""预计算阶段：生成 MUSCL-2 所需的所有索引表和面法向量。

核心输出（每个方向 d ∈ {0,1,2}）：
    stencil_idx[d]  : ndarray (4, N_faces_d)  MUSCL-2 四点模板的全局单元索引
    left_cell[d]    : ndarray (N_faces_d,)    面左侧真实单元全局索引
    right_cell[d]   : ndarray (N_faces_d,)    面右侧真实单元全局索引
    face_normal[d]  : ndarray (N_faces_d, 3)  面法向量（含面积权重，即 kxyz * |J|）
    face_vol_L[d]   : ndarray (N_faces_d,)    左侧单元体积
    face_vol_R[d]   : ndarray (N_faces_d,)    右侧单元体积

约定：
    对于方向 d，面 f 介于格点 (d_idx=i) 和 (d_idx=i+1) 之间。
    左单元 = (i)，右单元 = (i+1)。
    只统计内部面（i ∈ [ghost-1, ghost+n-1)），即含一侧 ghost 的所有面，
    这样 ghost 层的 BC 值可以参与通量计算。
    最终 RHS 只累加到真实单元（通过 true_idx 掩码）。
"""
from __future__ import annotations
import numpy as np
import jax.numpy as jnp

from euler1d_jax.mesh.domain import DomainData, BlockInfo


def build_precomp(domain: DomainData, kxyz_flat: jnp.ndarray, vol_flat: jnp.ndarray):
    """生成三个方向的索引表和法向量。

    Parameters
    ----------
    domain    : DomainData，含 blocks 列表
    kxyz_flat : (N_total, 3, 3) jnp 数组，度量导数
    vol_flat  : (N_total,)      jnp 数组，单元体积

    Returns
    -------
    precomp : dict with keys 0, 1, 2 (directions), each dict containing:
        'stencil' : (4, N_faces)  int32 全局索引
        'left'    : (N_faces,)    int32 全局索引
        'right'   : (N_faces,)    int32 全局索引
        'normal'  : (N_faces, 3)  float64 面法向量（面积加权）
        'vol_L'   : (N_faces,)    float64 左单元体积
        'vol_R'   : (N_faces,)    float64 右单元体积
    """
    kxyz_np = np.array(kxyz_flat)   # 在 numpy 中建立索引更快
    vol_np  = np.array(vol_flat)

    precomp = {}
    for d in range(3):
        stencil_list = []
        left_list    = []
        right_list   = []
        normal_list  = []
        vol_L_list   = []
        vol_R_list   = []

        for bi in domain.blocks:
            _collect_faces_one_block(
                bi, d, kxyz_np, vol_np,
                stencil_list, left_list, right_list,
                normal_list, vol_L_list, vol_R_list
            )

        precomp[d] = dict(
            stencil = np.stack(stencil_list, axis=1).astype(np.int32),   # (4, N_faces)
            left    = np.array(left_list,  dtype=np.int32),               # (N_faces,)
            right   = np.array(right_list, dtype=np.int32),
            normal  = np.stack(normal_list, axis=0).astype(np.float64),  # (N_faces, 3)
            vol_L   = np.array(vol_L_list, dtype=np.float64),
            vol_R   = np.array(vol_R_list, dtype=np.float64),
        )

    # 转换为 jnp
    for d in range(3):
        p = precomp[d]
        precomp[d] = dict(
            stencil = jnp.array(p['stencil']),
            left    = jnp.array(p['left']),
            right   = jnp.array(p['right']),
            normal  = jnp.array(p['normal']),
            vol_L   = jnp.array(p['vol_L']),
            vol_R   = jnp.array(p['vol_R']),
        )
    return precomp


def _collect_faces_one_block(
    bi: BlockInfo, d: int,
    kxyz_np: np.ndarray, vol_np: np.ndarray,
    stencil_list, left_list, right_list,
    normal_list, vol_L_list, vol_R_list,
):
    """为块 bi 的方向 d 收集所有内部面的索引和法向量。

    方向映射：d=0 → i 轴 (ni), d=1 → j 轴 (nj), d=2 → k 轴 (nk)
    面的定义：介于 (idx=i) 和 (idx=i+1) 之间，i ∈ [g-1, g+n-1)
    其中 g=ghost，n=真实格数。
    包含紧靠 ghost 层的面，使得边界 ghost 值能参与通量计算。

    MUSCL-2 四点模板（左侧重构）：i-1, i, i+1, i+2
    面 f 在 i 和 i+1 之间：
        stencil[0] = i-1, [1] = i, [2] = i+1, [3] = i+2
        left  = i   (左单元，接收 -flux）
        right = i+1 (右单元，接收 +flux）

    面法向量：kxyz[:, d, :] 在左右单元的平均，再乘以面积 |J|。
    简化做法：用面两侧单元的度量平均值作为面度量。
    """
    g  = bi.ghost
    ni_g, nj_g, nk_g = bi.ni_g, bi.nj_g, bi.nk_g

    # 方向 d 对应的真实格数和总格数
    if d == 0:
        n = bi.ni; total = ni_g
    elif d == 1:
        n = bi.nj; total = nj_g
    else:
        n = bi.nk; total = nk_g

    # 面索引范围：i ∈ [g-1, g+n-1) 使得 i 和 i+1 都在有效范围内
    face_start = g - 1      # 最左面（左侧为 ghost，右侧为真实或 ghost）
    face_end   = g + n - 1  # 最右面（左侧为真实，右侧为 ghost 或真实）

    # 遍历所有面
    for fi in range(face_start, face_end):
        # 四点模板在当前方向的索引
        im1 = max(fi - 1, 0)          # i-1，clamp 到边界
        i0  = fi                       # i
        i1  = fi + 1                   # i+1
        i2  = min(fi + 2, total - 1)  # i+2，clamp

        # 对另外两个方向遍历所有格点（含 ghost）
        if d == 0:
            for j in range(nj_g):
                for k in range(nk_g):
                    _add_face(bi, d, im1, i0, i1, i2, j, k,
                              kxyz_np, vol_np,
                              stencil_list, left_list, right_list,
                              normal_list, vol_L_list, vol_R_list)
        elif d == 1:
            for i in range(ni_g):
                for k in range(nk_g):
                    _add_face(bi, d, im1, i0, i1, i2, i, k,
                              kxyz_np, vol_np,
                              stencil_list, left_list, right_list,
                              normal_list, vol_L_list, vol_R_list)
        else:
            for i in range(ni_g):
                for j in range(nj_g):
                    _add_face(bi, d, im1, i0, i1, i2, i, j,
                              kxyz_np, vol_np,
                              stencil_list, left_list, right_list,
                              normal_list, vol_L_list, vol_R_list)


def _add_face(bi, d, im1, i0, i1, i2, p, q_idx,
              kxyz_np, vol_np,
              stencil_list, left_list, right_list,
              normal_list, vol_L_list, vol_R_list):
    """添加单个面的索引和法向。"""
    def gidx(d_val, p_val, q_val):
        if d == 0:
            return bi.global_idx(d_val, p_val, q_val)
        elif d == 1:
            return bi.global_idx(p_val, d_val, q_val)
        else:
            return bi.global_idx(p_val, q_val, d_val)

    g_im1 = gidx(im1, p, q_idx)
    g_i0  = gidx(i0,  p, q_idx)
    g_i1  = gidx(i1,  p, q_idx)
    g_i2  = gidx(i2,  p, q_idx)

    stencil_list.append(np.array([g_im1, g_i0, g_i1, g_i2], dtype=np.int32))
    left_list.append(g_i0)
    right_list.append(g_i1)

    # 面法向量：取左右单元度量的平均
    # kxyz 的第 d 列 kxyz[:, d, :] 是计算空间方向 d 的物理梯度
    # 法向量 = kxyz[d, :] * |J|（面积权重）
    k_L = kxyz_np[g_i0, d, :]    # (3,) = (∂ξ_d/∂x, ∂ξ_d/∂y, ∂ξ_d/∂z) 左单元
    k_R = kxyz_np[g_i1, d, :]    # (3,) 右单元
    k_avg = 0.5 * (k_L + k_R)
    normal_list.append(k_avg)     # 已含面积权重（kxyz 本身包含 Jacobian 信息）

    vol_L_list.append(vol_np[g_i0])
    vol_R_list.append(vol_np[g_i1])
