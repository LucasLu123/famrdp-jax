import jax.numpy as jnp
import pytest

from famrdp_jax.core.dtypes import assert_float64, check_array_f64

def test_assert_float64_accepts_f64():
    x = jnp.ones(3, dtype=jnp.float64)
    assert_float64(x, name="x")

def test_assert_float64_rejects_f32():
    x = jnp.ones(3, dtype=jnp.float32)
    with pytest.raises(TypeError, match="float64"):
        assert_float64(x, name="x")

def test_check_array_f64_walks_pytree():
    good = {"a": jnp.zeros(2), "b": [jnp.ones(3)]}
    check_array_f64(good)

    bad = {"a": jnp.zeros(2, dtype=jnp.float32)}
    with pytest.raises(TypeError):
        check_array_f64(bad)
