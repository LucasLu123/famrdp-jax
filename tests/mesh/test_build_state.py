from pathlib import Path
import jax.numpy as jnp

from famrdp_jax.mesh.build_state import build_initial_state

FG = Path(__file__).parent.parent / "fixtures" / "mini.grd"
FT = Path(__file__).parent.parent / "fixtures" / "mini.inp"

def test_build_state_with_ghost():
    state, topos = build_initial_state(FG, FT, ghost=2, q_init_fn=lambda shape: jnp.zeros(shape, dtype=jnp.float64))
    assert len(state.blocks) == 2
    # block 0: mini.grd 第一块 (4,3,2) + ghost=2 每维各 +4 → (8,7,6)
    assert state.blocks[0].xyz.shape == (3, 8, 7, 6)
    assert state.blocks[0].q.shape == (5, 8, 7, 6)

def test_build_state_types_are_correct():
    state, topos = build_initial_state(FG, FT, ghost=2, q_init_fn=lambda shape: jnp.zeros(shape, dtype=jnp.float64))
    assert state.blocks[0].q.dtype == jnp.float64
    assert state.blocks[0].xyz.dtype == jnp.float64
    assert len(topos) == 2
