import jax.numpy as jnp
import numpy as np
from famrdp_jax.physics.reconstruct.muscl2 import muscl2_reconstruct, minmod

def test_minmod_zero_on_opposite_signs():
    a = jnp.array([1.0, -1.0], dtype=jnp.float64)
    b = jnp.array([-1.0, 1.0], dtype=jnp.float64)
    np.testing.assert_allclose(np.asarray(minmod(a, b)), 0.0, atol=1e-15)

def test_minmod_min_magnitude():
    np.testing.assert_allclose(float(minmod(jnp.array(2.0), jnp.array(5.0))), 2.0)
    np.testing.assert_allclose(float(minmod(jnp.array(-3.0), jnp.array(-7.0))), -3.0)

def test_muscl2_linear_profile_lr_equal():
    ni = 10
    pv = jnp.stack([jnp.linspace(1.0, 10.0, ni, dtype=jnp.float64)] * 5, axis=0)[:, :, None, None]
    pv_L, pv_R = muscl2_reconstruct(pv, axis=1, limiter="minmod")
    np.testing.assert_allclose(
        np.asarray(pv_L[0, :, 0, 0]),
        np.asarray(pv_R[0, :, 0, 0]),
        rtol=1e-12,
    )
