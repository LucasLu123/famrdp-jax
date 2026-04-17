"""M1.0 xyz 对拍: JAX read_grd vs Fortran probe_write_xyz output."""
from pathlib import Path
import numpy as np
import pytest

from famrdp_jax.mesh.io_grd import read_grd
from validation.fortran_ref import read_xyz_probe
from validation.diff import compare, report

REF_DIR = Path(__file__).parent.parent.parent / "validation/references/test3/xyz"
GRD = Path(__file__).parent.parent.parent / "test3.grd"

_ref_exists = REF_DIR.exists() and any(REF_DIR.glob("block_*.bin"))


@pytest.mark.skipif(not GRD.exists(), reason="test3.grd 未就位")
@pytest.mark.skipif(not _ref_exists, reason="Fortran reference 未生成 (needs_reference: 运行带 -DPROBE_XYZ 的 Fortran 可执行并将输出放入 validation/references/test3/xyz/)")
def test_xyz_match_fortran():
    """Verify JAX read_grd produces identical coordinates to Fortran mb_xyz."""
    blocks = read_grd(GRD)
    failures = []
    for nb, xyz in enumerate(blocks, start=1):
        ref_xyz = read_xyz_probe(REF_DIR, nb)
        for ci, coord in enumerate(["x", "y", "z"]):
            res = compare(
                np.asarray(xyz[ci]),
                ref_xyz[ci],
                name=f"block{nb}.{coord}",
                atol=1e-14,
            )
            if not res["ok"]:
                failures.append(report(res))
    assert not failures, "\n".join(failures)
