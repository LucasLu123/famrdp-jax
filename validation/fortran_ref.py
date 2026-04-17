from __future__ import annotations
from pathlib import Path
import numpy as np


def read_probe(path, expect_ndim: int) -> np.ndarray:
    """Read a Fortran unformatted stream binary written by probe_utils.f90.

    Format: int32[ndim] shape, float64[prod(shape)] data (Fortran column-major order)
    """
    with open(path, "rb") as f:
        shape = tuple(int(x) for x in np.fromfile(f, dtype="<i4", count=expect_ndim))
        total = int(np.prod(shape))
        arr = np.fromfile(f, dtype="<f8", count=total).reshape(shape, order="F")
    return arr


def read_xyz_probe(probe_dir: Path, block_id: int) -> np.ndarray:
    """Read a block xyz probe written by probe_write_xyz.

    Returns array of shape (3, ni, nj, nk) where axis 0 is [x, y, z].
    File format: int32[3] (ni, nj, nk), float64[ni*nj*nk] x, float64[ni*nj*nk] y,
    float64[ni*nj*nk] z (all Fortran column-major order).
    """
    path = probe_dir / f"block_{block_id:04d}.bin"
    with open(path, "rb") as f:
        ni, nj, nk = np.fromfile(f, dtype="<i4", count=3)
        x = np.fromfile(f, dtype="<f8", count=ni*nj*nk).reshape((ni, nj, nk), order="F")
        y = np.fromfile(f, dtype="<f8", count=ni*nj*nk).reshape((ni, nj, nk), order="F")
        z = np.fromfile(f, dtype="<f8", count=ni*nj*nk).reshape((ni, nj, nk), order="F")
    return np.stack([x, y, z], axis=0)
