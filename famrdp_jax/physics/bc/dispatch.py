from __future__ import annotations
import jax.numpy as jnp

from famrdp_jax.core.constants import BCType, Face
from famrdp_jax.core.types import Block, BlockTopology, Config
from famrdp_jax.physics.bc.wall import apply_wall_noslip_adiabatic, apply_wall_noslip_isothermal, apply_wall_slip
from famrdp_jax.physics.bc.symmetry import apply_symmetry_plane
from famrdp_jax.physics.bc.farfield import apply_farfield


def apply_bc(block: Block, topo: BlockTopology, cfg: Config) -> Block:
    q = block.q
    ghost = cfg.ghost
    gas = cfg.gas
    for face in Face:
        bct = topo.bc_type.get(face)
        if bct is None or bct == BCType.CUT1TO1:
            continue
        bc_p = cfg.bc_params.get(face, {})
        if bct == BCType.WALL:
            wall_type = bc_p.get("wall_type", "adiabatic")
            if wall_type == "isothermal":
                T_wall = bc_p.get("T_wall", 300.0)
                q = apply_wall_noslip_isothermal(q, face, ghost, T_wall, gas.gamma, gas.R)
            else:
                q = apply_wall_noslip_adiabatic(q, face, ghost)
        elif bct == BCType.SYMMETRY:
            q = apply_symmetry_plane(q, face, ghost)
        elif bct == BCType.FARFIELD:
            q_inf = jnp.array(bc_p.get("q_inf", [1.0, 0.0, 0.0, 0.0, 2.5]), dtype=jnp.float64)
            q = apply_farfield(q, face, ghost, q_inf)
    return Block(q=q, xyz=block.xyz)
