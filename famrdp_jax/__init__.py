"""FAMRDP → JAX stage-1 port. Float64-only. GPU-optional."""
import jax

jax.config.update("jax_enable_x64", True)

__version__ = "0.1.0"
