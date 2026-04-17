# tests/conftest.py
import jax
import pytest

jax.config.update("jax_enable_x64", True)

@pytest.fixture(scope="session")
def rtol():
    return 1e-12

@pytest.fixture(scope="session")
def atol_strict():
    return 1e-13
