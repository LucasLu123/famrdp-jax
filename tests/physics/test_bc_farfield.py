import jax.numpy as jnp
import numpy as np
from famrdp_jax.core.constants import Face
from famrdp_jax.physics.bc.farfield import apply_farfield

def test_farfield_sets_ghost_to_qinf():
    ni = 8
    q = jnp.ones((5, ni, 4, 4), dtype=jnp.float64) * 99.0
    q_inf = jnp.array([1.0, 0.5, 0.0, 0.0, 2.5], dtype=jnp.float64)
    ghost = 2
    q_bc = apply_farfield(q, face=Face.IMAX, ghost=ghost, q_inf=q_inf)
    for g in range(ghost):
        i = ni - 1 - g
        np.testing.assert_allclose(float(q_bc[0, i, 0, 0]), 1.0, atol=1e-14)
        np.testing.assert_allclose(float(q_bc[1, i, 0, 0]), 0.5, atol=1e-14)

def test_farfield_imin_sets_ghost():
    ni = 8
    q = jnp.ones((5, ni, 4, 4), dtype=jnp.float64) * 99.0
    q_inf = jnp.array([1.2, 0.3, 0.1, 0.0, 3.0], dtype=jnp.float64)
    ghost = 2
    q_bc = apply_farfield(q, face=Face.IMIN, ghost=ghost, q_inf=q_inf)
    for g in range(ghost):
        np.testing.assert_allclose(float(q_bc[0, g, 0, 0]), 1.2, atol=1e-14)
        np.testing.assert_allclose(float(q_bc[1, g, 0, 0]), 0.3, atol=1e-14)
    # Interior should be untouched
    np.testing.assert_allclose(float(q_bc[0, ghost, 0, 0]), 99.0, atol=1e-14)

def test_farfield_does_not_change_interior():
    ni = 8
    q = jnp.ones((5, ni, 4, 4), dtype=jnp.float64) * 5.0
    q_inf = jnp.array([1.0, 0.0, 0.0, 0.0, 2.5], dtype=jnp.float64)
    ghost = 2
    q_bc = apply_farfield(q, face=Face.JMIN, ghost=ghost, q_inf=q_inf)
    # Interior j-indices [ghost..nj-1] should be unchanged
    for j in range(ghost, 4):
        np.testing.assert_allclose(float(q_bc[0, 0, j, 0]), 5.0, atol=1e-14)
