import jax
import jax.numpy as jnp
import pytest

from famrdp_jax.core.constants import BCType, Face
from famrdp_jax.core.types import (
    Block, Metrics, State, BlockTopology, Config, GasModel, SchemeChoice,
)

def test_block_is_pytree():
    ni, nj, nk = 8, 6, 4
    b = Block(
        q=jnp.zeros((5, ni, nj, nk), dtype=jnp.float64),
        xyz=jnp.zeros((3, ni, nj, nk), dtype=jnp.float64),
    )
    leaves = jax.tree_util.tree_leaves(b)
    assert len(leaves) == 2

def test_block_rejects_float32():
    with pytest.raises(TypeError):
        Block(
            q=jnp.zeros((5, 2, 2, 2), dtype=jnp.float32),
            xyz=jnp.zeros((3, 2, 2, 2), dtype=jnp.float64),
        )

def test_state_tree_map_preserves_structure():
    ni, nj, nk = 4, 4, 4
    b = Block(
        q=jnp.zeros((5, ni, nj, nk)),
        xyz=jnp.zeros((3, ni, nj, nk)),
    )
    m = Metrics(
        jac=jnp.ones((ni, nj, nk)),
        kxyz=jnp.zeros((3, 3, ni, nj, nk)),
        vol=jnp.ones((ni, nj, nk)),
    )
    s = State(blocks=(b,), metrics=(m,), t=0.0, step=0)
    doubled = jax.tree_util.tree_map(lambda x: 2 * x, s)
    assert doubled.blocks[0].q.shape == (5, ni, nj, nk)

def test_topology_holds_neighbors():
    topo = BlockTopology(
        block_id=0,
        neighbors={Face.IMIN: None, Face.IMAX: (1, Face.IMIN, 0)},
        bc_type={Face.IMIN: BCType.WALL, Face.IMAX: BCType.CUT1TO1},
    )
    assert topo.bc_type[Face.IMIN] == BCType.WALL
