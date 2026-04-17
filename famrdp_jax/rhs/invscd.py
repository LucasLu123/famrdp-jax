"""Inviscid (Euler) RHS assembler for structured curvilinear grids.

For each computational direction xi (axis=1,2,3):
  1. MUSCL-2 reconstruct prim vars -> pv_L, pv_R on faces
  2. Convert to conserved -> UL, UR
  3. Compute face normal from kxyz metrics
  4. Roe numerical flux F_face
  5. Accumulate flux divergence into RHS
"""
from __future__ import annotations
import jax.numpy as jnp
from famrdp_jax.core.types import Metrics
from famrdp_jax.physics.eos import prim_to_cons, cons_to_prim
from famrdp_jax.physics.reconstruct.muscl2 import muscl2_reconstruct
from famrdp_jax.flux.roe import roe_flux


def compute_invscd_rhs(
    q,
    metrics: Metrics,
    *,
    ghost: int,
    gamma: float,
    reconstruct: str,
    flux: str,
    limiter: str,
):
    """Compute inviscid (Euler) RHS.

    Parameters
    ----------
    q          : conservative vars, shape (5, ni, nj, nk)
    metrics    : grid metrics (jac, kxyz, vol)
    ghost      : number of ghost cells
    gamma      : specific heat ratio
    reconstruct: 'muscl2'
    flux       : 'roe'
    limiter    : 'minmod' or 'vanleer'

    Returns
    -------
    rhs : shape (5, ni, nj, nk)
    """
    pv = cons_to_prim(q, gamma)  # (5, ni, nj, nk)
    rhs = jnp.zeros_like(q)

    # axis in [1, 2, 3] maps to array axes [0, 1, 2] for i/j/k directions
    # Note: q has shape (5, ni, nj, nk), so physical dims are at axes 1,2,3
    # but pv also has shape (5, ni, nj, nk).
    # muscl2_reconstruct works on axis parameter directly on the array.
    # Physical i,j,k axes are at indices 1,2,3 of the 4-D array.

    for axis in [1, 2, 3]:
        n = pv.shape[axis]  # length along this direction (4D pv array, axis in 1,2,3)
        dir_idx = axis - 1  # 0,1,2 for i,j,k

        # MUSCL-2 reconstruction.
        # pv_L, pv_R have axis-dim length n-3
        # They correspond to faces i+1/2 for i = 1 .. n-3
        # i.e. face index f in [0, n-4] corresponds to the face between
        # cells (f+1) and (f+2).
        pv_L, pv_R = muscl2_reconstruct(pv, axis=axis, limiter=limiter)
        UL = prim_to_cons(pv_L, gamma)
        UR = prim_to_cons(pv_R, gamma)

        # Face metric: kxyz[dir_idx, :, ...] gives the metric vector for direction dir_idx
        # k_dir shape: (3, ni, nj, nk) – a 4D array with leading dim=3 for xyz components
        k_dir = metrics.kxyz[dir_idx]  # (3, ni, nj, nk)

        def slc4(arr, start, end, ax=axis):
            """Slice a 4D array (first dim is component dim) along ax in [1,2,3]."""
            idx = [slice(None)] * arr.ndim
            idx[ax] = slice(start, end if end != 0 else None)
            return arr[tuple(idx)]

        def slc3(arr, start, end, ax=dir_idx):
            """Slice a 3D array along ax in [0,1,2]."""
            idx = [slice(None)] * arr.ndim
            idx[ax] = slice(start, end if end != 0 else None)
            return arr[tuple(idx)]

        # Faces between cells i=1..n-3 and i+1=2..n-2.
        # Average the metric at adjacent cell centers to get face metric.
        # muscl2: face f corresponds to cells f+1 (left) and f+2 (right)
        # So k_dir at face f = avg(k_dir[f+1], k_dir[f+2])
        # face f in [0, n-4], so cells [1, n-3] and [2, n-2]
        k_L = slc4(k_dir, 1, n - 2)   # (3, ..., n-3, ...)
        k_R = slc4(k_dir, 2, n - 1)   # (3, ..., n-3, ...)
        n_face = 0.5 * (k_L + k_R)   # (3, ..., n-3, ...)

        # Face normal magnitude (used as area weight)
        # k_dir has shape (3, ni, nj, nk), so sum over axis=0 (the 3 components)
        n_mag = jnp.sqrt(jnp.sum(n_face ** 2, axis=0)) + 1e-300  # (ni_face, nj_face, nk_face)
        n_unit = n_face / n_mag[None]  # unit normal

        # Roe flux: F_face shape (5, ni_face, nj_face, nk_face) where n_face = n-3
        F_face = roe_flux(UL, UR, n_unit, gamma=gamma)

        # Scale by face area (n_mag acts as Jacobian * area factor)
        F_area = F_face * n_mag[None]  # (5, ni_face, nj_face, nk_face)

        # Flux divergence:
        # For cell c, the contribution from direction 'axis' is:
        #   -(F_area[right_face] - F_area[left_face]) / vol[c]
        #
        # F_area has n-3 faces (f = 0..n-4).
        # Face f corresponds to cells (f+1) and (f+2), so:
        #   LEFT face of cell c  = face (c-2)  [if c >= 2]
        #   RIGHT face of cell c = face (c-1)  [if c <= n-2]
        #
        # For cells c = 2..n-3:
        #   right face index in F_area = c-1   (valid: 1..n-4)
        #   left  face index in F_area = c-2   (valid: 0..n-5)
        #
        # div = F_area[c-1] - F_area[c-2]  for c = 2..n-3
        #     = F_area[1:n-3] - F_area[0:n-4]
        # That gives length n-4, covering cells c=2..n-3 (same length n-4).

        F_right = slc4(F_area, 1, None)   # faces 1..n-4, shape (5, ..., n-4, ...)
        F_left  = slc4(F_area, 0, -1)     # faces 0..n-5, shape (5, ..., n-4, ...)
        div = F_right - F_left             # (5, ..., n-4, ...)

        # Volume of cells c = 2..n-3 (3D array, use dir_idx = axis-1)
        vol_cell = slc3(metrics.vol, 2, n - 2)  # (ni-4|nj-4|nk-4, ...) 3D

        rhs_dir = -div / (vol_cell[None] + 1e-300)  # (5, ..., n-4, ...)

        # Accumulate into rhs at cells c = 2..n-3 along axis
        idx_rhs = [slice(None)] * rhs.ndim
        idx_rhs[axis] = slice(2, n - 2)
        rhs = rhs.at[tuple(idx_rhs)].add(rhs_dir)

    return rhs
