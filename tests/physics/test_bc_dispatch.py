import jax.numpy as jnp
import numpy as np
from famrdp_jax.core.constants import BCType, Face
from famrdp_jax.core.types import Block, BlockTopology, Config, GasModel, SchemeChoice
from famrdp_jax.physics.bc.dispatch import apply_bc


def _make_block(ni=8, nj=4, nk=4, rho=1.0, u=2.0, v=3.0, w=4.0, p=1.0, gamma=1.4):
    rhoE = p / (gamma - 1.0) + 0.5 * rho * (u**2 + v**2 + w**2)
    q = jnp.zeros((5, ni, nj, nk), dtype=jnp.float64)
    q = q.at[0].set(rho).at[1].set(rho*u).at[2].set(rho*v).at[3].set(rho*w).at[4].set(rhoE)
    xyz = jnp.zeros((3, ni, nj, nk), dtype=jnp.float64)
    return Block(q=q, xyz=xyz)


def _make_cfg(bc_type_dict, bc_params_dict=None, ghost=2):
    gas = GasModel()
    scheme = SchemeChoice(reconstruct="MUSCL2PV", flux="ROE", rk_order=3, limiter="minmod")
    topo = BlockTopology(block_id=0, neighbors={}, bc_type=bc_type_dict)
    cfg = Config(
        gas=gas,
        scheme=scheme,
        bc_params=bc_params_dict or {},
        topology=(topo,),
        ghost=ghost,
    )
    return topo, cfg


def test_dispatch_wall_adiabatic_imin():
    block = _make_block()
    topo, cfg = _make_cfg({Face.IMIN: BCType.WALL})
    block_bc = apply_bc(block, topo, cfg)
    # Ghost layers should have velocity negated
    for g in range(cfg.ghost):
        np.testing.assert_allclose(float(block_bc.q[1, g, 0, 0]), -2.0, atol=1e-14)
        np.testing.assert_allclose(float(block_bc.q[0, g, 0, 0]),  1.0, atol=1e-14)


def test_dispatch_symmetry_jmin():
    block = _make_block()
    topo, cfg = _make_cfg({Face.JMIN: BCType.SYMMETRY})
    block_bc = apply_bc(block, topo, cfg)
    # JMIN: v (rho*v index=2) should be negated
    np.testing.assert_allclose(float(block_bc.q[2, 0, 0, 0]), -3.0, atol=1e-14)
    np.testing.assert_allclose(float(block_bc.q[1, 0, 0, 0]),  2.0, atol=1e-14)


def test_dispatch_farfield_imax():
    block = _make_block()
    q_inf_vals = [1.2, 0.5, 0.1, 0.0, 3.0]
    topo, cfg = _make_cfg(
        {Face.IMAX: BCType.FARFIELD},
        bc_params_dict={Face.IMAX: {"q_inf": q_inf_vals}}
    )
    block_bc = apply_bc(block, topo, cfg)
    ni = block.q.shape[1]
    ghost = cfg.ghost
    for g in range(ghost):
        i = ni - 1 - g
        np.testing.assert_allclose(float(block_bc.q[0, i, 0, 0]), q_inf_vals[0], atol=1e-14)
        np.testing.assert_allclose(float(block_bc.q[1, i, 0, 0]), q_inf_vals[1], atol=1e-14)


def test_dispatch_cut1to1_skipped():
    block = _make_block()
    topo, cfg = _make_cfg({Face.IMIN: BCType.CUT1TO1})
    block_bc = apply_bc(block, topo, cfg)
    # CUT1TO1 should be skipped; q should be unchanged
    np.testing.assert_allclose(np.asarray(block_bc.q), np.asarray(block.q), atol=1e-14)
