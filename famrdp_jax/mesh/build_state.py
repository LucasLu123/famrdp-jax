from __future__ import annotations
from pathlib import Path
from typing import Callable

import jax.numpy as jnp
import numpy as np

from famrdp_jax.core.types import Block, State, Metrics
from famrdp_jax.mesh.io_grd import read_grd
from famrdp_jax.mesh.io_top import parse_topology


def _pad_xyz(xyz_interior: np.ndarray, ghost: int) -> np.ndarray:
    return np.pad(
        xyz_interior,
        ((0, 0), (ghost, ghost), (ghost, ghost), (ghost, ghost)),
        mode="edge",
    )


def build_initial_state(
    grd_path,
    inp_path,
    ghost: int,
    q_init_fn: Callable,
) -> tuple[State, tuple]:
    xyz_list = read_grd(grd_path)
    topos = tuple(parse_topology(inp_path))
    blocks = []
    metrics_list = []
    for xyz in xyz_list:
        xyz_pad = _pad_xyz(xyz, ghost)
        _, ni, nj, nk = xyz_pad.shape
        q = q_init_fn((5, ni, nj, nk))
        blocks.append(Block(
            q=jnp.asarray(q, dtype=jnp.float64),
            xyz=jnp.asarray(xyz_pad, dtype=jnp.float64),
        ))
        metrics_list.append(Metrics(
            jac=jnp.ones((ni, nj, nk), dtype=jnp.float64),
            kxyz=jnp.zeros((3, 3, ni, nj, nk), dtype=jnp.float64),
            vol=jnp.ones((ni, nj, nk), dtype=jnp.float64),
        ))
    state = State(blocks=tuple(blocks), metrics=tuple(metrics_list), t=0.0, step=0)
    return state, topos
