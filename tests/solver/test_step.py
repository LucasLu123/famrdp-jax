import jax.numpy as jnp
import numpy as np
from famrdp_jax.core.constants import BCType, Face
from famrdp_jax.core.types import (Block, Metrics, State, BlockTopology, Config, GasModel, SchemeChoice)
from famrdp_jax.solver.step import step

def test_single_block_uniform_stays_uniform():
    gamma = 1.4
    ni, nj, nk = 16, 12, 10
    ghost = 2
    rho, u, v, w, p = 1.225, 50.0, 0.0, 0.0, 101325.0
    rhoE = p/(gamma-1.0) + 0.5*rho*(u**2+v**2+w**2)
    q = jnp.zeros((5, ni, nj, nk), dtype=jnp.float64)
    q = q.at[0].set(rho).at[1].set(rho*u).at[2].set(rho*v).at[3].set(rho*w).at[4].set(rhoE)
    xyz = jnp.zeros((3, ni, nj, nk), dtype=jnp.float64)
    jac = jnp.ones((ni, nj, nk), dtype=jnp.float64)
    kxyz = jnp.zeros((3, 3, ni, nj, nk), dtype=jnp.float64)
    for d in range(3): kxyz = kxyz.at[d, d].set(1.0)
    m = Metrics(jac=jac, kxyz=kxyz, vol=jac)
    q_inf = jnp.array([rho, rho*u, rho*v, rho*w, rhoE], dtype=jnp.float64)
    topo = BlockTopology(
        block_id=0,
        neighbors={f: None for f in Face},
        bc_type={f: BCType.FARFIELD for f in Face},
    )
    cfg = Config(
        gas=GasModel(),
        scheme=SchemeChoice(reconstruct="muscl2", flux="roe", rk_order=3, limiter="minmod"),
        bc_params={face: {"q_inf": q_inf} for face in Face},
        topology=(topo,),
        ghost=ghost,
        nan_check=False,
    )
    state = State(blocks=(Block(q=q, xyz=xyz),), metrics=(m,), t=0.0, step=0)
    dt = 1e-6
    for _ in range(3):
        state = step(state, cfg, dt)
    np.testing.assert_allclose(
        np.asarray(state.blocks[0].q[:, 4:-4, 4:-4, 4:-4]),
        np.asarray(q[:, 4:-4, 4:-4, 4:-4]),
        rtol=1e-5, atol=1e-5,
    )
