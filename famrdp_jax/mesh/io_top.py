"""Parse GridGen .inp topology file.

Format (text, Fortran-style column-aligned integers):
    Line 1 : ignored header int (number of zones / version marker)
    Line 2 : nblocks
    --- repeated nblocks times ---
    ni nj nk          (grid dimensions for this block)
    block_label        (1-char or longer name, ignored)
    nbc                (number of BC *face* entries declared for this block;
                        interface entries occupy 2 lines each in the file
                        but count as one logical entry in nbc)
    --- repeated (total face lines, which can exceed nbc for interfaces) ---
    is ie js je ks ke  value

Face-line semantics
-------------------
For physical BCs (WALL, SYMMETRY, FARFIELD …):
    value > 1  →  BCType enum value

For CUT1TO1 (inter-block) interfaces, lines come in *pairs*:
    First line  : current block's face, value == -1 (CUT1TO1 marker)
    Second line : neighbor's face in the neighbor block's coordinate system,
                  value == neighbor_block_id  (1-based)

A degenerate index range (is==ie, js==je, or ks==ke) identifies which
face of the 6 structured faces the entry belongs to:
    is == ie → IMIN (is==1) or IMAX (is==ni)
    js == je → JMIN (js==1) or JMAX (js==nj)
    ks == ke → KMIN (ks==1) or KMAX (ks==nk)

When k values are -nk and -1 (interface placeholder), the face is
determined from the i/j range.
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional, Tuple, List

from famrdp_jax.core.constants import BCType, Face
from famrdp_jax.core.types import BlockTopology


def _face_from_range(
    is_: int, ie: int,
    js: int, je: int,
    ks: int, ke: int,
    ni: int, nj: int, nk: int,
) -> Face:
    """Identify which of the 6 faces this range corresponds to.

    Handles the case where interface lines use placeholder values
    (negative) for the degenerate axis.
    """
    # Degenerate in i?
    if is_ == ie:
        if is_ == 1:
            return Face.IMIN
        else:
            return Face.IMAX
    # Degenerate in j?
    if js == je:
        if js == 1 or (js < 0):
            # js negative is used in interface placeholder lines (abs(js) ~ 1)
            return Face.JMIN
        else:
            return Face.JMAX
    # Degenerate in k?
    if ks == ke:
        if ks == 1 or ks <= 0:
            return Face.KMIN
        else:
            return Face.KMAX
    # Fallback: pick the axis with smallest absolute range
    ranges = [abs(ie - is_), abs(je - js), abs(ke - ks)]
    axis = ranges.index(min(ranges))
    if axis == 0:
        return Face.IMIN if is_ <= ni // 2 else Face.IMAX
    elif axis == 1:
        return Face.JMIN if js <= nj // 2 else Face.JMAX
    else:
        return Face.KMIN if ks <= nk // 2 else Face.KMAX


def _detect_face(vals: List[int], ni: int, nj: int, nk: int) -> Face:
    """Determine face from the 6-value index range, interpreting abs values."""
    is_, ie, js, je, ks, ke = vals
    # For interface placeholder lines, k values may be negative.
    # Use absolute values but keep sign for min/max detection.
    return _face_from_range(abs(is_), abs(ie), abs(js), abs(je), abs(ks), abs(ke),
                            ni, nj, nk)


def _read_tokens(path: Path) -> List[List[int | str]]:
    """Read all lines of the file, stripping blank lines, returning token lists."""
    lines = []
    with open(path, "r") as fh:
        for raw in fh:
            stripped = raw.strip()
            if stripped:
                lines.append(stripped.split())
    return lines


def parse_topology(path) -> list[BlockTopology]:
    """Parse a GridGen .inp file and return one BlockTopology per block.

    Parameters
    ----------
    path : str or Path
        Path to the .inp topology file.

    Returns
    -------
    list[BlockTopology]
        One entry per block, 0-indexed block_id.
    """
    lines = _read_tokens(Path(path))
    idx = 0

    def next_ints(n: int) -> List[int]:
        nonlocal idx
        row = lines[idx]; idx += 1
        return [int(x) for x in row[:n]]

    def next_row_ints() -> List[int]:
        nonlocal idx
        row = lines[idx]; idx += 1
        return [int(x) for x in row]

    def next_label() -> str:
        nonlocal idx
        row = lines[idx]; idx += 1
        return " ".join(row)

    # Line 1: ignored header / zone count
    idx += 1  # skip
    # Line 2: nblocks
    nblocks = int(lines[idx][0]); idx += 1

    results: list[BlockTopology] = []

    for block_id in range(nblocks):
        # Grid dimensions
        dims = next_ints(3)
        ni, nj, nk = dims

        # Block label (may contain letters + spaces)
        label = next_label()

        # Number of BC logical entries
        nbc = int(lines[idx][0]); idx += 1

        # Read face lines: for interfaces, they come in pairs.
        # The total face lines = nbc + (number_of_interface_pairs)
        # We read until we've consumed nbc logical entries.
        bc_type: dict[Face, BCType] = {}
        neighbors: dict[Face, Optional[Tuple]] = {}

        logical_count = 0
        while logical_count < nbc:
            # Read first line of potential pair
            row = next_row_ints()
            is_, ie, js, je, ks, ke = row[0], row[1], row[2], row[3], row[4], row[5]
            value = row[6]

            if value == BCType.CUT1TO1:
                # This is the first line of an interface pair.
                # Identify the current block's face.
                my_face = _detect_face([is_, ie, js, je, ks, ke], ni, nj, nk)

                # Read second line: neighbor's corresponding face
                row2 = next_row_ints()
                is2, ie2, js2, je2, ks2, ke2 = row2[0], row2[1], row2[2], row2[3], row2[4], row2[5]
                neighbor_block_id = row2[6]  # 1-based
                nbr_face = _detect_face([is2, ie2, js2, je2, ks2, ke2], ni, nj, nk)

                bc_type[my_face] = BCType.CUT1TO1
                # Store neighbor as (0-based block_id, neighbor_face, orientation=1)
                neighbors[my_face] = (neighbor_block_id - 1, nbr_face, 1)

                logical_count += 1
            else:
                # Physical BC
                my_face = _detect_face([is_, ie, js, je, ks, ke], ni, nj, nk)
                bc_type[my_face] = BCType(value)
                neighbors[my_face] = None

                logical_count += 1

        results.append(BlockTopology(
            block_id=block_id,
            neighbors=neighbors,
            bc_type=bc_type,
        ))

    return results
