"""TVD-RK3 time integration (Shu-Osher 1988)."""
from __future__ import annotations


def rk3_step(q, L, dt):
    """One step of the 3rd-order TVD Runge-Kutta scheme.

    Parameters
    ----------
    q  : current state (JAX array or pytree)
    L  : callable q -> dq/dt  (the RHS operator)
    dt : time-step size (scalar)

    Returns
    -------
    q_new : updated state after one dt
    """
    q1 = q + dt * L(q)
    q2 = 0.75 * q + 0.25 * (q1 + dt * L(q1))
    return (1.0 / 3.0) * q + (2.0 / 3.0) * (q2 + dt * L(q2))
