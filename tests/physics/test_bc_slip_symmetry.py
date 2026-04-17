import jax.numpy as jnp
import numpy as np
from famrdp_jax.core.constants import Face
from famrdp_jax.physics.bc.wall import apply_wall_slip
from famrdp_jax.physics.bc.symmetry import apply_symmetry_plane

def _const_q(ni=8, rho=1.0, u=2.0, v=3.0, w=4.0, p=5.0, gamma=1.4):
    q = jnp.zeros((5, ni, 4, 4), dtype=jnp.float64)
    rhoE = p/(gamma-1.0) + 0.5*rho*(u*u+v*v+w*w)
    q = q.at[0].set(rho).at[1].set(rho*u).at[2].set(rho*v).at[3].set(rho*w).at[4].set(rhoE)
    return q

def test_slip_wall_imin_reflects_u_only():
    q = _const_q()
    q_bc = apply_wall_slip(q, face=Face.IMIN, ghost=2)
    np.testing.assert_allclose(float(q_bc[1, 0, 0, 0]), -2.0, atol=1e-14)
    np.testing.assert_allclose(float(q_bc[2, 0, 0, 0]),  3.0, atol=1e-14)

def test_symmetry_jmin_flips_v():
    q = _const_q()
    q_bc = apply_symmetry_plane(q, face=Face.JMIN, ghost=2)
    np.testing.assert_allclose(float(q_bc[2, 0, 0, 0]), -3.0, atol=1e-14)
    np.testing.assert_allclose(float(q_bc[1, 0, 0, 0]),  2.0, atol=1e-14)

def test_slip_wall_kmin_reflects_w_only():
    q = _const_q()
    q_bc = apply_wall_slip(q, face=Face.KMIN, ghost=2)
    # KMIN: vel_component=3 (rhow), so w flips, u and v unchanged
    np.testing.assert_allclose(float(q_bc[3, 0, 0, 0]), -4.0, atol=1e-14)
    np.testing.assert_allclose(float(q_bc[1, 0, 0, 0]),  2.0, atol=1e-14)
    np.testing.assert_allclose(float(q_bc[2, 0, 0, 0]),  3.0, atol=1e-14)

def test_symmetry_kmax_flips_w():
    q = _const_q()
    nk = q.shape[3]
    q_bc = apply_symmetry_plane(q, face=Face.KMAX, ghost=2)
    # KMAX outer=[nk-1, nk-2]=[3,2]
    np.testing.assert_allclose(float(q_bc[3, 0, 0, nk-1]), -4.0, atol=1e-14)
    np.testing.assert_allclose(float(q_bc[1, 0, 0, nk-1]),  2.0, atol=1e-14)
