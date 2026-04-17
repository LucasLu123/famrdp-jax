import jax.numpy as jnp
import numpy as np
from famrdp_jax.core.types import GasModel
from famrdp_jax.physics.eos import cons_to_prim, prim_to_cons, speed_of_sound

def test_cons_to_prim_roundtrip():
    gas = GasModel()
    rho, u, v, w, p = 1.225, 50.0, 10.0, -5.0, 101325.0
    prim = jnp.array([rho, u, v, w, p], dtype=jnp.float64)[:, None, None, None]
    q = prim_to_cons(prim, gas.gamma)
    prim_back = cons_to_prim(q, gas.gamma)
    np.testing.assert_allclose(np.asarray(prim_back), np.asarray(prim), rtol=1e-14)

def test_speed_of_sound():
    gas = GasModel()
    a = speed_of_sound(rho=1.225, p=101325.0, gamma=gas.gamma)
    np.testing.assert_allclose(a, (1.4 * 101325.0 / 1.225)**0.5, rtol=1e-12)
