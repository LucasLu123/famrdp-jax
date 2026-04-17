from __future__ import annotations
import jax.numpy as jnp
from famrdp_jax.core.types import State, Block, Config
from famrdp_jax.mesh.halo import halo_exchange
from famrdp_jax.physics.bc.dispatch import apply_bc
from famrdp_jax.rhs.invscd import compute_invscd_rhs
from famrdp_jax.rhs.viscous import compute_viscous_rhs


def _rhs_block(q, metrics, cfg):
    inv = compute_invscd_rhs(q, metrics, ghost=cfg.ghost, gamma=cfg.gas.gamma,
                              reconstruct=cfg.scheme.reconstruct, flux=cfg.scheme.flux,
                              limiter=cfg.scheme.limiter)
    vis = compute_viscous_rhs(q, metrics, gas=cfg.gas, ghost=cfg.ghost)
    return inv + vis


def step(state: State, cfg: Config, dt: float) -> State:
    """Advance the simulation state by one time step using TVD-RK3 (Shu-Osher).

    Parameters
    ----------
    state : current simulation state
    cfg   : simulation configuration
    dt    : time step size (seconds)

    Returns
    -------
    State with blocks advanced by dt, t incremented by dt, step incremented by 1.
    """
    # Halo exchange for CUT1TO1 block interfaces
    blocks = halo_exchange(state.blocks, cfg.topology, cfg.ghost)
    # Apply boundary conditions
    blocks = tuple(apply_bc(b, t, cfg) for b, t in zip(blocks, cfg.topology))

    # RK3 Shu-Osher (TVD) time integration:
    #   U1 = U0 + dt * L(U0)
    #   U2 = 3/4 * U0 + 1/4 * (U1 + dt * L(U1))
    #   U_new = 1/3 * U0 + 2/3 * (U2 + dt * L(U2))
    qs0 = tuple(b.q for b in blocks)

    def L(qs_tuple):
        return tuple(_rhs_block(q, m, cfg) for q, m in zip(qs_tuple, state.metrics))

    L0 = L(qs0)
    qs1 = tuple(q + dt * r for q, r in zip(qs0, L0))

    L1 = L(qs1)
    qs2 = tuple(0.75 * q + 0.25 * (q1 + dt * r) for q, q1, r in zip(qs0, qs1, L1))

    L2 = L(qs2)
    qs_new = tuple((1.0 / 3.0) * q + (2.0 / 3.0) * (q2 + dt * r)
                   for q, q2, r in zip(qs0, qs2, L2))

    new_blocks = tuple(Block(q=qn, xyz=b.xyz) for qn, b in zip(qs_new, blocks))
    return State(blocks=new_blocks, metrics=state.metrics, t=state.t + dt, step=state.step + 1)
