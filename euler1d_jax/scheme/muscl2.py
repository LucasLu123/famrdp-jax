"""MUSCL-2 重构：基于 gather 索引的 1D 全局数组版本。

接口：
    muscl2_reconstruct_1d(prime, stencil_idx, limiter) -> (pL, pR)

其中：
    prime      : (N_total, 5)  全局原始变量数组
    stencil_idx: (4, N_faces)  四点模板全局索引
    limiter    : 'minmod' | 'vanleer'
    pL, pR     : (N_faces, 5)  左/右重构值
"""
from __future__ import annotations
import jax.numpy as jnp


def _minmod(a: jnp.ndarray, b: jnp.ndarray) -> jnp.ndarray:
    s = jnp.sign(a) * jnp.maximum(jnp.sign(a) * jnp.sign(b), 0.0)
    return s * jnp.minimum(jnp.abs(a), jnp.abs(b))


def _vanleer(a: jnp.ndarray, b: jnp.ndarray) -> jnp.ndarray:
    denom = jnp.where(jnp.abs(a + b) > 1e-30, a + b, 1.0)
    return jnp.where(a * b > 0.0, 2.0 * a * b / denom, 0.0)


_LIMITERS = {"minmod": _minmod, "vanleer": _vanleer}


def muscl2_reconstruct_1d(
    prime: jnp.ndarray,
    stencil_idx: jnp.ndarray,
    *,
    limiter: str = "minmod",
) -> tuple:
    """MUSCL-2 左/右界面重构。

    Parameters
    ----------
    prime       : (N_total, 5) 原始变量 [rho, u, v, w, p]
    stencil_idx : (4, N_faces) int32，四点模板全局索引
                  行 0=i-1, 1=i, 2=i+1, 3=i+2
    limiter     : 限制器名称

    Returns
    -------
    pL : (N_faces, 5) 左界面值（由单元 i 重构）
    pR : (N_faces, 5) 右界面值（由单元 i+1 重构）
    """
    lim = _LIMITERS[limiter]

    # gather: prime[stencil_idx[r]] -> (N_faces, 5) for each row r
    pmm = prime[stencil_idx[0]]   # (N_faces, 5)  i-1
    pc  = prime[stencil_idx[1]]   # i
    pp1 = prime[stencil_idx[2]]   # i+1
    pp2 = prime[stencil_idx[3]]   # i+2

    slope_L = lim(pc - pmm, pp1 - pc)    # 左斜率
    slope_R = lim(pp1 - pc, pp2 - pp1)   # 右斜率

    pL = pc  + 0.5 * slope_L   # 左单元右界面值
    pR = pp1 - 0.5 * slope_R   # 右单元左界面值
    return pL, pR
