"""Validation test: compare computed metrics against Fortran reference data.

Skipped automatically when test3.grd or the metric reference directory is absent.
"""
from pathlib import Path

import numpy as np
import pytest

from famrdp_jax.mesh.io_grd import read_grd
from famrdp_jax.mesh.metric import compute_metrics
from validation.fortran_ref import read_probe
from validation.diff import compare, report

REF_DIR = Path(__file__).parent.parent.parent / "validation/references/test3/metric"
GRD = Path(__file__).parent.parent.parent / "test3.grd"


def _ref_exists():
    return REF_DIR.exists() and any(REF_DIR.iterdir())


@pytest.mark.skipif(not GRD.exists(), reason="test3.grd 未就位")
@pytest.mark.skipif(not _ref_exists(), reason="metric reference 未生成")
def test_metric_match_fortran():
    import jax.numpy as jnp

    blocks = read_grd(GRD)
    fails = []
    for nb, xyz in enumerate(blocks, start=1):
        m = compute_metrics(jnp.asarray(xyz), ghost=0)
        ref_jac = read_probe(REF_DIR / f"block_{nb:04d}_jac.bin", expect_ndim=3)
        res = compare(np.asarray(m.jac), ref_jac, name=f"b{nb}.jac", atol=1e-13)
        if not res["ok"]:
            fails.append(report(res))
    assert not fails, "\n".join(fails)
