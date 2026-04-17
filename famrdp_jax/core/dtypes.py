import jax
import jax.numpy as jnp


def assert_float64(arr, *, name: str) -> None:
    if arr.dtype != jnp.float64:
        raise TypeError(
            f"{name!r} must be float64, got {arr.dtype}. "
            "检查 jax.config.update('jax_enable_x64', True) 是否在导入前调用。"
        )


def check_array_f64(tree) -> None:
    leaves = jax.tree_util.tree_leaves(tree)
    for i, leaf in enumerate(leaves):
        if hasattr(leaf, "dtype"):
            assert_float64(leaf, name=f"leaf[{i}]")
