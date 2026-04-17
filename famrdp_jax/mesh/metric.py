"""Grid metric computation using 4th-order central differences.

Aligns with Fortran ncutpol=3 (4-order central):
    ∂x/∂ξ = (-f[i+2] + 8 f[i+1] - 8 f[i-1] + f[i-2]) / 12

J_mat[m, l] = ∂x_m/∂ξ_l   (m = x/y/z, l = ξ/η/ζ)
jac  = det(J_mat)
kxyz[l, m] = inv(J_mat)[l, m] = ∂ξ_l/∂x_m
vol  = jac
"""
from __future__ import annotations

import jax.numpy as jnp

from famrdp_jax.core.types import Metrics


def _central4(f: jnp.ndarray, axis: int) -> jnp.ndarray:
    """4th-order central difference along `axis`.

    Note: jnp.roll wraps at boundaries; ghost cells near boundaries will
    have incorrect values.  Interior region (beyond ghost=2 layers) is accurate.
    """
    fp2 = jnp.roll(f, -2, axis=axis)
    fp1 = jnp.roll(f, -1, axis=axis)
    fm1 = jnp.roll(f,  1, axis=axis)
    fm2 = jnp.roll(f,  2, axis=axis)
    return (-fp2 + 8.0 * fp1 - 8.0 * fm1 + fm2) / 12.0


def compute_metrics(xyz: jnp.ndarray, *, ghost: int) -> Metrics:
    """Compute grid metrics from coordinate array.

    Parameters
    ----------
    xyz : jnp.ndarray
        Shape (3, ni, nj, nk).  Axes 1/2/3 correspond to ξ/η/ζ.
    ghost : int
        Number of ghost layers (informational; roll-based difference taints
        the outermost ``ghost`` layers on each side, but no clipping is done
        here — callers should restrict to interior ``[2*ghost:-2*ghost]``).

    Returns
    -------
    Metrics
        jac  : (ni, nj, nk)
        kxyz : (3, 3, ni, nj, nk)  — kxyz[l, m] = ∂ξ_l/∂x_m
        vol  : (ni, nj, nk)  (same as jac for structured grids)
    """
    # J_mat_ml[m, l] = ∂x_m/∂ξ_l,  shape (3, 3, ni, nj, nk)
    # xyz[m] has shape (ni, nj, nk); differentiation axes are 0/1/2
    J_mat_ml = jnp.stack(
        [
            jnp.stack(
                [_central4(xyz[m], axis=l) for l in range(3)],
                axis=0,
            )
            for m in range(3)
        ],
        axis=0,
    )  # (3, 3, ni, nj, nk)

    # Move spatial axes first for linalg operations: (ni, nj, nk, 3, 3)
    J_per_cell = jnp.transpose(J_mat_ml, (2, 3, 4, 0, 1))

    jac = jnp.linalg.det(J_per_cell)          # (ni, nj, nk)
    kxyz_per = jnp.linalg.inv(J_per_cell)     # (ni, nj, nk, 3, 3)  [l, m]
    kxyz = jnp.transpose(kxyz_per, (3, 4, 0, 1, 2))  # (3, 3, ni, nj, nk)

    return Metrics(jac=jac, kxyz=kxyz, vol=jac)
