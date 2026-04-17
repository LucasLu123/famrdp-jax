# tests/solver/test_grad.py
import jax
import jax.numpy as jnp
import numpy as np
import pytest
from famrdp_jax.core.constants import BCType, Face
from famrdp_jax.core.types import Block, Metrics, State, BlockTopology, Config, GasModel, SchemeChoice
from famrdp_jax.solver.step import step

def test_grad_wrt_q_is_nonzero_finite():
    gamma = 1.4
    ni, nj, nk = 12, 10, 8
    ghost = 2
    rho, u, v, w, p = 1.225, 50.0, 0.0, 0.0, 101325.0
    rhoE = p/(gamma-1.0) + 0.5*rho*(u**2+v**2+w**2)
    q_inf = jnp.array([rho, rho*u, rho*v, rho*w, rhoE], dtype=jnp.float64)
    q = jnp.zeros((5, ni, nj, nk), dtype=jnp.float64)
    q = q.at[0].set(rho).at[1].set(rho*u).at[2].set(rho*v).at[3].set(rho*w).at[4].set(rhoE)
    xyz = jnp.zeros((3, ni, nj, nk), dtype=jnp.float64)
    jac = jnp.ones((ni, nj, nk), dtype=jnp.float64)
    kxyz = jnp.zeros((3, 3, ni, nj, nk), dtype=jnp.float64)
    for d in range(3): kxyz = kxyz.at[d, d].set(1.0)
    m = Metrics(jac=jac, kxyz=kxyz, vol=jac)
    topo = BlockTopology(block_id=0, neighbors={f: None for f in Face},
                         bc_type={f: BCType.FARFIELD for f in Face})
    cfg = Config(gas=GasModel(), scheme=SchemeChoice(reconstruct="muscl2", flux="roe", rk_order=3, limiter="minmod"),
                 bc_params={face: {"q_inf": q_inf} for face in Face},
                 topology=(topo,), ghost=ghost, nan_check=False)

    def loss(q0):
        state = State(blocks=(Block(q=q0, xyz=xyz),), metrics=(m,))
        state1 = step(state, cfg, dt=1e-6)
        return jnp.sum(state1.blocks[0].q)

    grad_q = jax.grad(loss)(q)
    assert jnp.all(jnp.isfinite(grad_q)), "gradient contains inf/nan"
    assert jnp.any(grad_q != 0.0), "gradient is all zeros"
    print(f"grad max: {float(jnp.max(jnp.abs(grad_q))):.3e}")
