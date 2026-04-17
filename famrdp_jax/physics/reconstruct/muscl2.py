from __future__ import annotations
import jax.numpy as jnp


def minmod(a: jnp.ndarray, b: jnp.ndarray) -> jnp.ndarray:
    s = jnp.sign(a) * jnp.maximum(jnp.sign(a) * jnp.sign(b), 0.0)
    return s * jnp.minimum(jnp.abs(a), jnp.abs(b))


def vanleer(a: jnp.ndarray, b: jnp.ndarray) -> jnp.ndarray:
    safe_denom = jnp.where(jnp.abs(a + b) > 1e-30, a + b, 1.0)
    return jnp.where(a * b > 0.0, 2.0 * a * b / safe_denom, 0.0)


_LIMITERS = {"minmod": minmod, "vanleer": vanleer}


def muscl2_reconstruct(pv: jnp.ndarray, *, axis: int, limiter: str = "minmod"):
    lim = _LIMITERS[limiter]
    n = pv.shape[axis]

    def slc(start, end):
        idx = [slice(None)] * pv.ndim
        idx[axis] = slice(start, end)
        return pv[tuple(idx)]

    # 需要 4 点: pv[i-1], pv[i], pv[i+1], pv[i+2]
    # face j 对应左 cell=j+1, 范围 j: 0..n-4
    pmm = slc(0, n-3)   # pv[i-1], i=1..n-3
    pc  = slc(1, n-2)   # pv[i]
    pp1 = slc(2, n-1)   # pv[i+1]
    pp2 = slc(3, n)     # pv[i+2]
    slope_L = lim(pc - pmm, pp1 - pc)
    slope_R = lim(pp1 - pc, pp2 - pp1)
    pv_L = pc  + 0.5 * slope_L
    pv_R = pp1 - 0.5 * slope_R
    return pv_L, pv_R
