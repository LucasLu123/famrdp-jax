from __future__ import annotations
import jax.numpy as jnp
from famrdp_jax.core.constants import Face
from famrdp_jax.physics.bc.base import face_axis, ghost_indices, vel_component


def _mirror_ghost(q, face, ghost, sign_mask):
    ax = face_axis(face)
    ni = q.shape[ax]
    outer, inner = ghost_indices(face, ghost, ni)
    sign_bc = jnp.array(sign_mask, dtype=jnp.float64)
    # After indexing q at a single integer index along ax, the result has ndim-1 dimensions.
    # sign shape: (5,) -> reshape to (5, 1, ..., 1) with (q.ndim - 1) total dims
    # axis 0 is the variable axis; the remaining (q.ndim-2) dims are spatial (excluding ax)
    out_ndim = q.ndim - 1  # shape of qi after integer indexing
    shape = [1] * out_ndim
    shape[0] = 5
    sign_bc = sign_bc.reshape(shape)
    for o, n in zip(outer, inner):
        idx_o = [slice(None)] * q.ndim
        idx_n = [slice(None)] * q.ndim
        idx_o[ax] = o
        idx_n[ax] = n
        qi = q[tuple(idx_n)]
        patched = qi * sign_bc
        q = q.at[tuple(idx_o)].set(patched)
    return q


def apply_wall_noslip_adiabatic(q, face, ghost):
    sign = jnp.array([1.0, -1.0, -1.0, -1.0, 1.0], dtype=jnp.float64)
    return _mirror_ghost(q, face, ghost, sign)


def apply_wall_noslip_isothermal(q, face, ghost, T_wall, gamma, R):
    ax = face_axis(face)
    ni = q.shape[ax]
    outer, inner = ghost_indices(face, ghost, ni)
    for o, n in zip(outer, inner):
        idx_o = [slice(None)] * q.ndim
        idx_n = [slice(None)] * q.ndim
        idx_o[ax] = o
        idx_n[ax] = n
        qi = q[tuple(idx_n)]
        rho_i = qi[0]
        ke_i = 0.5 * (qi[1]**2 + qi[2]**2 + qi[3]**2) / rho_i
        p_i = (gamma - 1.0) * (qi[4] - ke_i)
        p_g = p_i
        rho_g = p_g / (R * T_wall)
        scale = rho_g / rho_i
        rhou_g = -qi[1] * scale
        rhov_g = -qi[2] * scale
        rhow_g = -qi[3] * scale
        ke_g = 0.5 * (rhou_g**2 + rhov_g**2 + rhow_g**2) / rho_g
        rhoE_g = p_g / (gamma - 1.0) + ke_g
        patched = jnp.stack([rho_g, rhou_g, rhov_g, rhow_g, rhoE_g], axis=0)
        q = q.at[tuple(idx_o)].set(patched)
    return q


def apply_wall_slip(q, face, ghost):
    vc = vel_component(face)
    sign = jnp.ones(5, dtype=jnp.float64).at[vc].set(-1.0)
    return _mirror_ghost(q, face, ghost, sign)
