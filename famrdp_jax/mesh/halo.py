"""Halo exchange: copy ghost-cell data between adjacent blocks (1-to-1 cuts).

Each block has q.shape = (5, ni, nj, nk) where the first/last `ghost` layers
along each axis are ghost cells populated from the neighbouring block.
"""
from __future__ import annotations
import jax.numpy as jnp

from famrdp_jax.core.constants import BCType, Face
from famrdp_jax.core.types import Block, BlockTopology

_FACE_AXIS = {
    Face.IMIN: 1, Face.IMAX: 1,
    Face.JMIN: 2, Face.JMAX: 2,
    Face.KMIN: 3, Face.KMAX: 3,
}


def _ghost_slice(face: Face, ghost: int) -> slice:
    """Slice selecting the ghost layer on the given face of the target block."""
    if face in (Face.IMIN, Face.JMIN, Face.KMIN):
        return slice(0, ghost)
    return slice(-ghost, None)


def _interior_adj_slice(face: Face, ghost: int) -> slice:
    """Slice selecting the interior cells adjacent to `face` in the source block.

    These are the cells that will be copied into the neighbouring block's ghosts.
    """
    if face in (Face.IMIN, Face.JMIN, Face.KMIN):
        return slice(ghost, 2 * ghost)
    return slice(-2 * ghost, -ghost)


def _copy_face(
    tgt_q: jnp.ndarray,
    src_q: jnp.ndarray,
    tgt_face: Face,
    src_face: Face,
    ghost: int,
) -> jnp.ndarray:
    """Return `tgt_q` with its ghost layer on `tgt_face` filled from `src_q`."""
    ax = _FACE_AXIS[tgt_face]
    sl_t = _ghost_slice(tgt_face, ghost)
    sl_s = _interior_adj_slice(src_face, ghost)

    idx_t = [slice(None)] * tgt_q.ndim
    idx_s = [slice(None)] * src_q.ndim
    idx_t[ax] = sl_t
    idx_s[_FACE_AXIS[src_face]] = sl_s

    return tgt_q.at[tuple(idx_t)].set(src_q[tuple(idx_s)])


def halo_exchange(
    blocks: tuple[Block, ...],
    topos: tuple[BlockTopology, ...],
    ghost: int,
) -> tuple[Block, ...]:
    """Perform one round of halo exchange for all CUT1TO1 interfaces.

    Parameters
    ----------
    blocks:
        Tuple of Block objects (one per block). Arrays are JAX arrays.
    topos:
        Tuple of BlockTopology objects with matching order/indices.
    ghost:
        Number of ghost layers on each face.

    Returns
    -------
    Tuple of Block objects with ghost layers updated from neighbouring blocks.
    """
    id_to_idx = {t.block_id: i for i, t in enumerate(topos)}
    qs = [b.q for b in blocks]

    for i, topo in enumerate(topos):
        for face, info in topo.neighbors.items():
            if info is None:
                continue
            if topo.bc_type.get(face) != BCType.CUT1TO1:
                continue
            nb_id, src_face, _orient = info
            if _orient != 0:
                raise NotImplementedError(f"orient={_orient} not supported in stage 1")
            src_i = id_to_idx[nb_id]
            # Use original (pre-exchange) source arrays to avoid order-dependency
            qs[i] = _copy_face(qs[i], blocks[src_i].q, face, src_face, ghost)

    return tuple(Block(q=qs[i], xyz=blocks[i].xyz) for i in range(len(blocks)))
