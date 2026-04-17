import jax.numpy as jnp
import numpy as np

from famrdp_jax.mesh.metric import compute_metrics


def _uniform_cube_xyz(ni, nj, nk):
    xs = np.linspace(0.0, 1.0, ni)
    ys = np.linspace(0.0, 2.0, nj)
    zs = np.linspace(0.0, 3.0, nk)
    X, Y, Z = np.meshgrid(xs, ys, zs, indexing="ij")
    return jnp.asarray(np.stack([X, Y, Z], axis=0), dtype=jnp.float64)


def test_uniform_cube_jac_equals_volume():
    ni, nj, nk = 20, 20, 20
    xyz = _uniform_cube_xyz(ni, nj, nk)
    m = compute_metrics(xyz, ghost=2)
    expected_j = (1.0/19) * (2.0/19) * (3.0/19)
    inside = m.jac[4:-4, 4:-4, 4:-4]
    np.testing.assert_allclose(np.asarray(inside), expected_j, rtol=1e-10)


def test_uniform_cube_kxyz_diagonal():
    ni, nj, nk = 20, 20, 20
    xyz = _uniform_cube_xyz(ni, nj, nk)
    m = compute_metrics(xyz, ghost=2)
    k = np.asarray(m.kxyz)
    np.testing.assert_allclose(k[0, 0, 4:-4, 4:-4, 4:-4], 19.0, rtol=1e-9)
    np.testing.assert_allclose(k[1, 0, 4:-4, 4:-4, 4:-4], 0.0, atol=1e-9)
