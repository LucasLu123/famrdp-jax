from __future__ import annotations
import jax.numpy as jnp
from famrdp_jax.core.constants import Face
from famrdp_jax.physics.bc.base import face_axis, ghost_indices


def apply_farfield(q, face, ghost, q_inf):
    """简化远场 BC:ghost 层直接设为自由流值。"""
    ax = face_axis(face)
    ni = q.shape[ax]
    outer, _ = ghost_indices(face, ghost, ni)
    for o in outer:
        idx = [slice(None)] * q.ndim
        idx[ax] = o
        qi = q[tuple(idx)]  # shape: (5, ...) with q.ndim-1 dims
        out_ndim = qi.ndim
        shape = [1] * out_ndim
        shape[0] = 5
        q_inf_reshaped = q_inf.reshape(shape) * jnp.ones_like(qi)
        q = q.at[tuple(idx)].set(q_inf_reshaped)
    return q
