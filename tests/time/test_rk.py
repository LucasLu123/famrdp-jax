import jax.numpy as jnp
import numpy as np
from famrdp_jax.time.rk import rk3_step


def test_rk3_integrates_linear_ode():
    def L(q): return -q
    q0 = jnp.array([1.0], dtype=jnp.float64)
    dt = 0.1
    q1 = rk3_step(q0, L, dt)
    expected = jnp.exp(jnp.array(-dt, dtype=jnp.float64))
    err = float(jnp.abs(q1[0] - expected))
    assert err < 1e-5
