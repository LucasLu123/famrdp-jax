from famrdp_jax.core.constants import Face

_AXIS = {
    Face.IMIN: 1, Face.IMAX: 1,
    Face.JMIN: 2, Face.JMAX: 2,
    Face.KMIN: 3, Face.KMAX: 3,
}

_VEL_IDX = {
    Face.IMIN: 1, Face.IMAX: 1,
    Face.JMIN: 2, Face.JMAX: 2,
    Face.KMIN: 3, Face.KMAX: 3,
}

def face_axis(face: Face) -> int:
    return _AXIS[face]

def vel_component(face: Face) -> int:
    return _VEL_IDX[face]

def ghost_indices(face: Face, ghost: int, ni_total: int):
    """返回 (outer_list, inner_list),逐一配对(mirror ghost layers)。
    outer[g] <- inner[g] 的镜像。
    IMIN: outer=[0,1,...,ghost-1], inner=[ghost-1, ghost-2,...,0] reversed
    实际:outer=[0,1], inner=[2*ghost-1, 2*ghost-2]==[3,2] 时 ghost=2
    """
    if face in (Face.IMIN, Face.JMIN, Face.KMIN):
        outer = list(range(ghost))
        inner = [2*ghost - 1 - g for g in outer]
    else:
        outer = list(range(ni_total - 1, ni_total - 1 - ghost, -1))
        inner = [ni_total - 2*ghost + (ni_total - 1 - o) for o in outer]
    return outer, inner
