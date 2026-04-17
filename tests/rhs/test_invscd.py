import jax.numpy as jnp
import numpy as np
from famrdp_jax.core.types import Metrics
from famrdp_jax.rhs.invscd import compute_invscd_rhs


def test_uniform_flow_zero_rhs():
    ni, nj, nk = 16, 12, 10
    gamma = 1.4
    rho, u, v, w, p = 1.225, 50.0, 0.0, 0.0, 101325.0
    rhoE = p / (gamma - 1.0) + 0.5 * rho * (u**2 + v**2 + w**2)
    q = jnp.zeros((5, ni, nj, nk), dtype=jnp.float64)
    q = q.at[0].set(rho).at[1].set(rho * u).at[2].set(rho * v).at[3].set(rho * w).at[4].set(rhoE)
    jac = jnp.ones((ni, nj, nk), dtype=jnp.float64)
    kxyz = jnp.zeros((3, 3, ni, nj, nk), dtype=jnp.float64)
    for d in range(3):
        kxyz = kxyz.at[d, d].set(1.0)
    m = Metrics(jac=jac, kxyz=kxyz, vol=jac)
    rhs = compute_invscd_rhs(
        q, m, ghost=2, gamma=gamma, reconstruct="muscl2", flux="roe", limiter="minmod"
    )
    # Interior region (excluding ghost + boundary truncation) should be ~0
    np.testing.assert_allclose(np.asarray(rhs[:, 4:-4, 4:-4, 4:-4]), 0.0, atol=1e-8)
