from __future__ import annotations
import jax.numpy as jnp


def prim_to_cons(prim: jnp.ndarray, gamma: float) -> jnp.ndarray:
    rho = prim[0]; u = prim[1]; v = prim[2]; w = prim[3]; p = prim[4]
    ke = 0.5 * (u*u + v*v + w*w)
    rhoE = p / (gamma - 1.0) + rho * ke
    return jnp.stack([rho, rho*u, rho*v, rho*w, rhoE], axis=0)


def cons_to_prim(q: jnp.ndarray, gamma: float) -> jnp.ndarray:
    rho = q[0]
    u = q[1] / rho; v = q[2] / rho; w = q[3] / rho
    ke = 0.5 * (u*u + v*v + w*w)
    p = (gamma - 1.0) * (q[4] - rho * ke)
    return jnp.stack([rho, u, v, w, p], axis=0)


def speed_of_sound(rho, p, gamma):
    return jnp.sqrt(gamma * p / rho)


def sutherland_viscosity(T, mu0, T0, Ts):
    return mu0 * (T / T0) ** 1.5 * (T0 + Ts) / (T + Ts)
