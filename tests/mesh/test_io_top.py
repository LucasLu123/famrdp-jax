from pathlib import Path
from famrdp_jax.core.constants import BCType, Face
from famrdp_jax.mesh.io_top import parse_topology

FIXTURE = Path(__file__).parent.parent / "fixtures" / "mini.inp"


def test_parse_topology_returns_n_blocks():
    tops = parse_topology(FIXTURE)
    assert len(tops) == 2


def test_block_ids():
    tops = parse_topology(FIXTURE)
    assert tops[0].block_id == 0
    assert tops[1].block_id == 1


def test_all_faces_have_bc_type():
    tops = parse_topology(FIXTURE)
    for top in tops:
        # Each block has 6 faces; all should be assigned
        assert len(top.bc_type) == 6


def test_bc_types_are_valid():
    tops = parse_topology(FIXTURE)
    for top in tops:
        for face, bct in top.bc_type.items():
            assert isinstance(bct, BCType)


def test_interface_faces_are_cut1to1():
    tops = parse_topology(FIXTURE)
    # Block 0 and 1 share a j-face interface
    for top in tops:
        cut_faces = [f for f, bt in top.bc_type.items() if bt == BCType.CUT1TO1]
        assert len(cut_faces) >= 1, f"Block {top.block_id} has no CUT1TO1 faces"


def test_interface_neighbors_are_set():
    tops = parse_topology(FIXTURE)
    for top in tops:
        for face, bct in top.bc_type.items():
            if bct == BCType.CUT1TO1:
                nbr = top.neighbors[face]
                assert nbr is not None, f"Block {top.block_id} face {face} has CUT1TO1 but no neighbor"
                nbr_block_id, nbr_face, orient = nbr
                assert 0 <= nbr_block_id < 2
                assert isinstance(nbr_face, Face)


def test_physical_bc_neighbors_are_none():
    tops = parse_topology(FIXTURE)
    for top in tops:
        for face, bct in top.bc_type.items():
            if bct != BCType.CUT1TO1:
                assert top.neighbors[face] is None


def test_parse_test3_inp():
    """Verify the real test3.inp file is parseable."""
    real_inp = Path(__file__).parent.parent.parent / "test3.inp"
    if not real_inp.exists():
        import pytest
        pytest.skip("test3.inp not found")
    tops = parse_topology(real_inp)
    assert len(tops) == 2
    for top in tops:
        assert len(top.bc_type) == 6
        for face, bct in top.bc_type.items():
            assert isinstance(bct, BCType)
