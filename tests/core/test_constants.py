from famrdp_jax.core.constants import BCType, FluxScheme, ReconstructScheme

def test_bctype_values_match_fortran():
    assert BCType.CUT1TO1.value == -1
    assert BCType.WALL.value == 2
    assert BCType.SYMMETRY.value == 3
    assert BCType.FARFIELD.value == 4
    assert BCType.INFLOW.value == 5
    assert BCType.OUTFLOW.value == 6
    assert BCType.POLE.value == 7
    assert BCType.PATCHED.value == 8

def test_flux_scheme_values():
    assert FluxScheme.ROE.value == 4
    assert FluxScheme.ROE_SCMP.value == 7

def test_reconstruct_scheme_values():
    assert ReconstructScheme.MUSCL2PV.value == 1
    assert ReconstructScheme.DCSH5PI.value == 23
