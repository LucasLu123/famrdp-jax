from __future__ import annotations
import jax.numpy as jnp
from famrdp_jax.core.constants import Face
from famrdp_jax.physics.bc.base import face_axis, ghost_indices, vel_component
from famrdp_jax.physics.bc.wall import _mirror_ghost


def apply_symmetry_plane(q, face, ghost):
    vc = vel_component(face)
    sign = jnp.ones(5, dtype=jnp.float64).at[vc].set(-1.0)
    return _mirror_ghost(q, face, ghost, sign)
