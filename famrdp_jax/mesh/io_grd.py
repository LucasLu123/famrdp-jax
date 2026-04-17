from pathlib import Path
import numpy as np


def read_grd(path) -> list[np.ndarray]:
    """Read a PLOT3D .grd file (Fortran unformatted binary stream, little-endian).

    Format:
        int32 nblocks
        repeat nblocks: int32 ni nj nk
        repeat nblocks: float64 x(ni*nj*nk), y(ni*nj*nk), z(ni*nj*nk)  [Fortran column-major]

    Returns a list of nblocks arrays, each shaped (3, ni, nj, nk), dtype float64.
    Axis 0 is the coordinate index: 0=x, 1=y, 2=z.
    """
    with open(path, "rb") as f:
        nblocks = int(np.fromfile(f, dtype="<i4", count=1)[0])
        shapes = []
        for _ in range(nblocks):
            shapes.append(tuple(np.fromfile(f, dtype="<i4", count=3).astype(int)))
        blocks = []
        for (ni, nj, nk) in shapes:
            ncell = ni * nj * nk
            xs = np.fromfile(f, dtype="<f8", count=ncell).reshape((ni, nj, nk))
            ys = np.fromfile(f, dtype="<f8", count=ncell).reshape((ni, nj, nk))
            zs = np.fromfile(f, dtype="<f8", count=ncell).reshape((ni, nj, nk))
            xyz = np.stack([xs, ys, zs], axis=0)
            blocks.append(xyz.astype(np.float64))
    return blocks
