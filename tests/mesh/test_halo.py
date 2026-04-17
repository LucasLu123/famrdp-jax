import jax.numpy as jnp
import numpy as np
import pytest

from famrdp_jax.core.constants import Face, BCType
from famrdp_jax.core.types import Block, BlockTopology
from famrdp_jax.mesh.halo import halo_exchange


def test_two_blocks_imax_to_imin_cut1to1():
    ghost = 2
    ni_in = 4
    ni = ni_in + 2 * ghost  # 8
    shape = (5, ni, 6, 4)
    q0 = jnp.zeros(shape, dtype=jnp.float64).at[:, ghost:ghost+ni_in, :, :].set(
        jnp.arange(1.0, ni_in+1.0, dtype=jnp.float64)[None, :, None, None]
    )
    q1 = jnp.zeros(shape, dtype=jnp.float64).at[:, ghost:ghost+ni_in, :, :].set(
        jnp.arange(10.0, 10.0+ni_in, dtype=jnp.float64)[None, :, None, None]
    )
    xyz = jnp.zeros((3,) + shape[1:], dtype=jnp.float64)
    b0 = Block(q=q0, xyz=xyz)
    b1 = Block(q=q1, xyz=xyz)

    topo0 = BlockTopology(
        block_id=0,
        neighbors={face: None for face in Face},
        bc_type={},
    )
    topo0 = BlockTopology(
        block_id=0,
        neighbors={**{face: None for face in Face}, Face.IMAX: (1, Face.IMIN, 0)},
        bc_type={Face.IMAX: BCType.CUT1TO1},
    )
    topo1 = BlockTopology(
        block_id=1,
        neighbors={**{face: None for face in Face}, Face.IMIN: (0, Face.IMAX, 0)},
        bc_type={Face.IMIN: BCType.CUT1TO1},
    )
    new = halo_exchange((b0, b1), (topo0, topo1), ghost=ghost)
    # 块 0 的 IMAX ghost(最后 2 层) = 块 1 IMIN 内部前 2 层 = 10.0, 11.0
    np.testing.assert_allclose(float(new[0].q[0, -ghost, 0, 0]), 10.0)
    np.testing.assert_allclose(float(new[0].q[0, -1, 0, 0]),     11.0)
    # 块 1 的 IMIN ghost(前 2 层) = 块 0 IMAX 内部最后 2 层 = 3.0, 4.0
    np.testing.assert_allclose(float(new[1].q[0, 0, 0, 0]), 3.0)
    np.testing.assert_allclose(float(new[1].q[0, 1, 0, 0]), 4.0)
