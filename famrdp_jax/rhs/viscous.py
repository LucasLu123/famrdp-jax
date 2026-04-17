"""Viscous RHS for the compressible Navier-Stokes equations.

Computes the divergence of the viscous flux tensor using central differences
in computational space, mapped back to physical space via the grid metrics.
"""
from __future__ import annotations
import jax.numpy as jnp
from famrdp_jax.core.types import Metrics, GasModel
from famrdp_jax.physics.eos import cons_to_prim, sutherland_viscosity


def _central2(f, axis):
    """2nd-order central difference d f / d xi (axis) using periodic-style roll."""
    fp = jnp.roll(f, -1, axis=axis)
    fm = jnp.roll(f,  1, axis=axis)
    return 0.5 * (fp - fm)


def compute_viscous_rhs(q, metrics: Metrics, *, gas: GasModel, ghost: int):
    """Compute viscous contribution to RHS.

    Parameters
    ----------
    q       : conservative variables, shape (5, ni, nj, nk)
    metrics : grid metrics (jac, kxyz, vol)
    gas     : gas model parameters (gamma, R, mu0, T0, Ts, Pr_lam)
    ghost   : number of ghost cells (currently unused but kept for API parity)

    Returns
    -------
    rhs : shape (5, ni, nj, nk) – viscous RHS  (energy and momentum components)
    """
    pv = cons_to_prim(q, gas.gamma)
    rho = pv[0]; u = pv[1]; v = pv[2]; w = pv[3]; p = pv[4]

    T = p / jnp.maximum(rho * gas.R, 1e-300)
    mu = sutherland_viscosity(T, gas.mu0, gas.T0, gas.Ts)
    cp = gas.gamma * gas.R / (gas.gamma - 1.0)
    kappa = cp * mu / gas.Pr_lam

    # -------------------------------------------------------------------------
    # Physical-space gradient of scalar phi:
    #   d phi / d x_m  =  sum_l  (d phi / d xi_l) * kxyz[l, m]
    # kxyz shape: (3, 3, ni, nj, nk)  where kxyz[l, m] = d xi_l / d x_m
    # -------------------------------------------------------------------------
    def grad_xyz(phi):
        """Return physical gradient, shape (3, ni, nj, nk)."""
        # d phi / d xi_l  for l=0,1,2  (axes 0,1,2 in 3-D array)
        dphi_dxi = jnp.stack([_central2(phi, axis=a) for a in (0, 1, 2)], axis=0)
        # d phi / d x_m = sum_l  dphi_dxi[l] * kxyz[l, m]
        return jnp.einsum("lijk,lmijk->mijk", dphi_dxi, metrics.kxyz)

    du = grad_xyz(u)   # (3, ni, nj, nk)
    dv = grad_xyz(v)
    dw = grad_xyz(w)
    dT = grad_xyz(T)

    div_u = du[0] + dv[1] + dw[2]

    # -------------------------------------------------------------------------
    # Viscous stress tensor  tau[i, j]  (shape (3, 3, ni, nj, nk))
    # tau_ij = mu * (du_i/dx_j + du_j/dx_i) - (2/3) mu (div u) delta_ij
    # -------------------------------------------------------------------------
    tau = jnp.stack([
        jnp.stack([2*mu*(du[0] - div_u/3.0),
                   mu*(du[1] + dv[0]),
                   mu*(du[2] + dw[0])]),
        jnp.stack([mu*(du[1] + dv[0]),
                   2*mu*(dv[1] - div_u/3.0),
                   mu*(dv[2] + dw[1])]),
        jnp.stack([mu*(du[2] + dw[0]),
                   mu*(dv[2] + dw[1]),
                   2*mu*(dw[2] - div_u/3.0)]),
    ])  # (3, 3, ni, nj, nk)

    heat = -kappa * dT  # (3, ni, nj, nk) – heat-flux vector q_j = -kappa dT/dx_j

    vel = jnp.stack([u, v, w])  # (3, ni, nj, nk)

    # -------------------------------------------------------------------------
    # Divergence of viscous flux:
    #   G_j = [0, tau_0j, tau_1j, tau_2j, sum_i u_i tau_ij - heat_j]
    #   RHS_var = sum_j  d G_j / d x_j
    #           = sum_j  sum_l  kxyz[l,j] * d G_j / d xi_l
    # -------------------------------------------------------------------------
    def div_physical_scalar(f):
        """Compute divergence of a vector field f (shape (3, ni, nj, nk))."""
        result = jnp.zeros_like(f[0])
        for j in range(3):
            df_j_dxi = jnp.stack([_central2(f[j], axis=a) for a in (0, 1, 2)], axis=0)
            # d f_j / d x_j = sum_l  kxyz[l, j] * d f_j / d xi_l
            result = result + jnp.einsum("lijk,lijk->ijk", df_j_dxi, metrics.kxyz[:, j])
        return result

    # Momentum viscous divergence: d tau_ij / d x_j  for each i
    rhs_mom = jnp.stack([
        div_physical_scalar(tau[i])   # shape (ni, nj, nk) for each row i
        for i in range(3)
    ])  # (3, ni, nj, nk)

    # Energy viscous divergence: d (u_i tau_ij - heat_j) / d x_j
    energy_flux = jnp.einsum("ijk,mijk->mijk", jnp.zeros_like(u), vel)  # placeholder zero
    # build vector G_energy[j] = sum_i u_i tau_ij - heat_j
    G_energy = jnp.stack([
        jnp.sum(vel * tau[:, j], axis=0) - heat[j]
        for j in range(3)
    ])  # (3, ni, nj, nk)
    rhs_energy = div_physical_scalar(G_energy)  # (ni, nj, nk)

    # Assemble: rhs[0]=mass (0), rhs[1:4]=momentum, rhs[4]=energy
    rhs = jnp.stack([
        jnp.zeros_like(u),    # continuity – no viscous contribution
        rhs_mom[0],
        rhs_mom[1],
        rhs_mom[2],
        rhs_energy,
    ], axis=0)  # (5, ni, nj, nk)

    return rhs
