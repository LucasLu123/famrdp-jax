from enum import IntEnum


class BCType(IntEnum):
    CUT1TO1   = -1
    WALL      = 2
    SYMMETRY  = 3
    FARFIELD  = 4
    INFLOW    = 5
    OUTFLOW   = 6
    POLE      = 7
    PATCHED   = 8


class WallSubType(IntEnum):
    ADIABATIC  = 1
    ISOTHERMAL = 2
    SLIP       = 3


class SymmetrySubType(IntEnum):
    POINT = 1
    PLANE = 2


class FarfieldSubType(IntEnum):
    RIEMANN = 1
    CHARACT = 2


class FluxScheme(IntEnum):
    STEGER    = 1
    SW_MOD    = 2
    VANLEER   = 3
    ROE       = 4
    ROE_PREC  = 5
    SLAU      = 6
    ROE_SCMP  = 7


class ReconstructScheme(IntEnum):
    MUSCL2PV = 1
    WCNS5PV  = 2
    DCSH5PI  = 23


class TimeScheme(IntEnum):
    RK3            = 1
    LUSGS_STD_SCA  = 2
    LUSGS_STD_PREC = 11
    LUSGS_STD_SCMP = 12


class Face(IntEnum):
    IMIN = 0
    IMAX = 1
    JMIN = 2
    JMAX = 3
    KMIN = 4
    KMAX = 5


N_EQN = 5  # rho, rho*u, rho*v, rho*w, rho*E
N_DIM = 3
