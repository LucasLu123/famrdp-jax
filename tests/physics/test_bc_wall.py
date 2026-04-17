import jax.numpy as jnp
import numpy as np
from famrdp_jax.core.constants import Face
from famrdp_jax.physics.bc.wall import apply_wall_noslip_adiabatic, apply_wall_noslip_isothermal

def test_noslip_adiabatic_zeros_velocity_mirrors_p_rho():
    ghost = 2
    ni = 8
    q = jnp.zeros((5, ni, 4, 4), dtype=jnp.float64)
    q = q.at[0].set(1.0).at[1].set(2.0).at[2].set(3.0).at[3].set(4.0).at[4].set(10.0)
    q_bc = apply_wall_noslip_adiabatic(q, face=Face.IMIN, ghost=ghost)
    for g in range(ghost):
        np.testing.assert_allclose(np.asarray(q_bc[0, g]),  1.0, atol=1e-14)
        np.testing.assert_allclose(np.asarray(q_bc[1, g]), -2.0, atol=1e-14)
        np.testing.assert_allclose(np.asarray(q_bc[2, g]), -3.0, atol=1e-14)
        np.testing.assert_allclose(np.asarray(q_bc[3, g]), -4.0, atol=1e-14)

def test_noslip_adiabatic_imax():
    ghost = 2
    ni = 8
    q = jnp.zeros((5, ni, 4, 4), dtype=jnp.float64)
    q = q.at[0].set(1.0).at[1].set(2.0).at[2].set(3.0).at[3].set(4.0).at[4].set(10.0)
    q_bc = apply_wall_noslip_adiabatic(q, face=Face.IMAX, ghost=ghost)
    # outer=[7,6], inner=[4,5]
    for g in [6, 7]:
        np.testing.assert_allclose(np.asarray(q_bc[0, g]),  1.0, atol=1e-14)
        np.testing.assert_allclose(np.asarray(q_bc[1, g]), -2.0, atol=1e-14)
        np.testing.assert_allclose(np.asarray(q_bc[2, g]), -3.0, atol=1e-14)
        np.testing.assert_allclose(np.asarray(q_bc[3, g]), -4.0, atol=1e-14)

def test_noslip_isothermal_pressure_set():
    ghost = 2
    ni = 8
    gamma = 1.4
    R = 287.058
    T_wall = 300.0
    # Set up a simple state: rho=1, u=2, v=3, w=4, E computed from p=1
    rho = 1.0
    u, v, w = 2.0, 3.0, 4.0
    p = 1.0
    rhoE = p / (gamma - 1.0) + 0.5 * rho * (u**2 + v**2 + w**2)
    q = jnp.zeros((5, ni, 4, 4), dtype=jnp.float64)
    q = q.at[0].set(rho).at[1].set(rho*u).at[2].set(rho*v).at[3].set(rho*w).at[4].set(rhoE)
    q_bc = apply_wall_noslip_isothermal(q, face=Face.IMIN, ghost=ghost, T_wall=T_wall, gamma=gamma, R=R)
    # Ghost density should be p/(R*T_wall)
    expected_rho_g = p / (R * T_wall)
    np.testing.assert_allclose(float(q_bc[0, 0, 0, 0]), expected_rho_g, rtol=1e-12)
