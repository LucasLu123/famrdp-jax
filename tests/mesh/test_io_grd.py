from pathlib import Path
import numpy as np
from famrdp_jax.mesh.io_grd import read_grd

FIXTURE = Path(__file__).parent.parent / "fixtures" / "mini.grd"

def test_read_grd_returns_per_block_xyz():
    blocks = read_grd(FIXTURE)
    assert len(blocks) == 2
    xyz0 = blocks[0]
    assert xyz0.shape == (3, 4, 3, 2)
    assert xyz0.dtype == np.float64

def test_read_grd_coords_match_meshgrid():
    blocks = read_grd(FIXTURE)
    xyz0 = blocks[0]
    np.testing.assert_allclose(xyz0[0, 0, :, :], 0.0, atol=1e-14)
    np.testing.assert_allclose(xyz0[0, -1, :, :], 1.0, atol=1e-14)
