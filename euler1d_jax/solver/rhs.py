"""Euler 方程无粘 RHS 装配。

计算策略：
    1. 对三个计算方向 d = 0, 1, 2 依次处理
    2. 用 precomp[d]['stencil'] gather 四点模板 → MUSCL-2 重构 → Roe 通量
    3. 通量散度：rhs[left] -= flux / vol_L, rhs[right] += flux / vol_R
    4. 仅 RHS 累加（不更新 ghost 层）；ghost 层由 BC 在步进函数中处理

注意法向量方向约定：
    precomp[d]['normal'] = kxyz[d, :] * （含Jacobian）
    即 ξ_d 方向的协变基向量，方向为从左到右（i → i+1）。
    Roe 通量 F = flux(UL, UR, n) 定义为从左流向右的通量。
    因此：
        rhs[left]  -= flux   (左单元流出)
        rhs[right] += flux   (右单元流入)
    最终 rhs = -dF/dξ / J，符号已含在除以 vol 中。
"""
from __future__ import annotations
import jax.numpy as jnp

from euler1d_jax.scheme.muscl2 import muscl2_reconstruct_1d
from euler1d_jax.scheme.roe import roe_flux_1d


def compute_euler_rhs(
    prime: jnp.ndarray,
    precomp: dict,
    true_idx: jnp.ndarray,
    *,
    gamma: float = 1.4,
    limiter: str = "minmod",
) -> jnp.ndarray:
    """计算 Euler 方程 RHS（三方向通量散度之和）。

    Parameters
    ----------
    prime    : (N_total, 5) 原始变量 [rho, u, v, w, p]
    precomp  : build_precomp() 返回的预计算数据
    true_idx : (N_true,) 真实单元全局索引（用于初始化 rhs 掩码）
    gamma    : 比热比
    limiter  : MUSCL-2 限制器

    Returns
    -------
    rhs : (N_total, 5)  只有 true_idx 处有非零值
    """
    N = prime.shape[0]
    rhs = jnp.zeros((N, 5), dtype=jnp.float64)

    for d in range(3):
        p = precomp[d]
        stencil  = p['stencil']    # (4, N_faces)
        left_idx = p['left']       # (N_faces,)
        right_idx= p['right']      # (N_faces,)
        normals  = p['normal']     # (N_faces, 3)
        vol_L    = p['vol_L']      # (N_faces,)
        vol_R    = p['vol_R']      # (N_faces,)

        # MUSCL-2 重构
        pL, pR = muscl2_reconstruct_1d(prime, stencil, limiter=limiter)
        # pL, pR : (N_faces, 5)

        # Roe 通量
        flux = roe_flux_1d(pL, pR, normals, gamma=gamma)
        # flux : (N_faces, 5)

        # 散度：flux / vol 累加到左右单元
        flux_over_volL = flux / vol_L[:, None]   # (N_faces, 5)
        flux_over_volR = flux / vol_R[:, None]

        rhs = rhs.at[left_idx].add(-flux_over_volL)
        rhs = rhs.at[right_idx].add(flux_over_volR)

    return rhs
