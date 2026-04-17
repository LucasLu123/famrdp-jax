import jax.numpy as jnp
import numpy as np
from famrdp_jax.core.types import Metrics, GasModel
from famrdp_jax.rhs.viscous import compute_viscous_rhs


def test_static_isothermal_zero_viscous_rhs():
    gas = GasModel()
    ni, nj, nk = 16, 12, 10
    rho, p = 1.225, 101325.0
    rhoE = p / (gas.gamma - 1.0)
    q = jnp.zeros((5, ni, nj, nk), dtype=jnp.float64)
    q = q.at[0].set(rho).at[4].set(rhoE)
    jac = jnp.ones((ni, nj, nk), dtype=jnp.float64)
    kxyz = jnp.zeros((3, 3, ni, nj, nk), dtype=jnp.float64)
    for d in range(3):
        kxyz = kxyz.at[d, d].set(1.0)
    m = Metrics(jac=jac, kxyz=kxyz, vol=jac)
    rhs = compute_viscous_rhs(q, m, gas=gas, ghost=2)
    np.testing.assert_allclose(np.asarray(rhs[:, 4:-4, 4:-4, 4:-4]), 0.0, atol=1e-10)
