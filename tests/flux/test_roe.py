import jax.numpy as jnp
import numpy as np
from famrdp_jax.flux.roe import roe_flux

def test_roe_symmetric_states_zero_flux_contribution():
    gamma = 1.4
    rho, u, v, w, p = 1.0, 50.0, 0.0, 0.0, 101325.0
    rhoE = p/(gamma-1.0) + 0.5*rho*(u**2+v**2+w**2)
    U = jnp.array([rho, rho*u, rho*v, rho*w, rhoE], dtype=jnp.float64)
    nxyz = jnp.array([1.0, 0.0, 0.0], dtype=jnp.float64)
    F = roe_flux(U, U, nxyz, gamma=gamma)
    un = u
    expected = jnp.array([rho*un, rho*un*un+p, rho*un*v, rho*un*w, u*(rhoE+p)], dtype=jnp.float64)
    np.testing.assert_allclose(np.asarray(F), np.asarray(expected), rtol=1e-12)

def test_roe_arbitrary_normal_identical():
    gamma = 1.4
    rho, u, v, w, p = 1.0, 30.0, 20.0, 10.0, 50000.0
    rhoE = p/(gamma-1.0) + 0.5*rho*(u**2+v**2+w**2)
    U = jnp.array([rho, rho*u, rho*v, rho*w, rhoE], dtype=jnp.float64)
    n = jnp.array([0.6, 0.8, 0.0], dtype=jnp.float64)
    F = roe_flux(U, U, n, gamma=gamma)
    un = u*n[0] + v*n[1] + w*n[2]
    H = (rhoE + p) / rho
    expected = jnp.array([rho*un, rho*u*un+p*n[0], rho*v*un+p*n[1], rho*w*un+p*n[2], rho*H*un], dtype=jnp.float64)
    np.testing.assert_allclose(np.asarray(F), np.asarray(expected), rtol=1e-12)
