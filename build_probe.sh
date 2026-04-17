#!/bin/bash
# Build FAMRDP Fortran solver with probe support (sequential, no MPI)
set -e

WDIR="/e/famrdp/.claude/worktrees/famrdp-jax"
BDIR="$WDIR/build-probe"
FC="gfortran"
FFLAGS="-cpp -O0 -fno-second-underscore"
PROBE_FLAGS="-DPROBE_XYZ -DPROBE_METRIC"
INC="-I$WDIR -I$WDIR/validation/fortran_probes"

mkdir -p "$BDIR"
cd "$BDIR"

SRC="$WDIR/src"
PROBE_SRC="$WDIR/validation/fortran_probes"

echo "=== Compiling probe_utils module ==="
$FC $FFLAGS $PROBE_FLAGS $INC -c "$PROBE_SRC/probe_utils.f90" -o probe_utils.o

echo "=== Compiling modules (order matters) ==="
$FC $FFLAGS $INC -c "$SRC/mod_kndconsts.f90"         -o mod_kndconsts.o
$FC $FFLAGS $INC -c "$SRC/mod_datatypes.f90"         -o mod_datatypes.o
$FC $FFLAGS $INC -c "$SRC/mod_constants.f90"         -o mod_constants.o
$FC $FFLAGS $INC -c "$SRC/mod_parallels.f90"         -o mod_parallels.o
$FC $FFLAGS $INC -c "$SRC/mod_fieldvars.f90"         -o mod_fieldvars.o
$FC $FFLAGS $INC -c "$SRC/mod_variables.f90"         -o mod_variables.o
$FC $FFLAGS $INC -c "$SRC/mod_opt_vars.f90"          -o mod_opt_vars.o
$FC $FFLAGS $INC -c "$SRC/mod_singulars.f90"         -o mod_singulars.o
$FC $FFLAGS $INC -c "$SRC/mod_runtimers.f90"         -o mod_runtimers.o
$FC $FFLAGS $INC -c "$SRC/mod_tecplotio.f90"         -o mod_tecplotio.o
$FC $FFLAGS $INC -c "$SRC/mod_interface.f90"         -o mod_interface.o

echo "=== Compiling implementation files ==="
for f in \
    util_math.f90 util_subs.f90 utilities.f90 misc.f90 \
    communicate.f90 bc.f90 eos.f90 \
    initialize.f90 preset.f90 \
    io_input.f90 io_output.f90 \
    metric.f90 rhs_invscd.f90 rhs_viscous.f90 rhs_source.f90 rhside.f90 \
    sol_rkutta.f90 sol_lusgs.f90 sol_lusgs_prec.f90 sol_lusgs_scmp.f90 \
    sol_prmat.f90 sol_prsca.f90 sol_prsca_scmp.f90 \
    tur_bc.f90 turbulent.f90 tur_menter.f90 tur_spalart.f90 tur_lusgs.f90 \
    solve.f90 main.f90; do
    echo "  $f"
    $FC $FFLAGS $PROBE_FLAGS $INC -c "$SRC/$f" -o "${f%.f90}.o"
done

echo "=== Linking ==="
$FC $FFLAGS -o "$BDIR/HOSTA.probe" \
    probe_utils.o \
    mod_kndconsts.o mod_datatypes.o mod_constants.o mod_parallels.o \
    mod_fieldvars.o mod_variables.o mod_opt_vars.o mod_singulars.o \
    mod_runtimers.o mod_tecplotio.o mod_interface.o \
    util_math.o util_subs.o utilities.o misc.o \
    communicate.o bc.o eos.o \
    initialize.o preset.o \
    io_input.o io_output.o \
    metric.o rhs_invscd.o rhs_viscous.o rhs_source.o rhside.o \
    sol_rkutta.o sol_lusgs.o sol_lusgs_prec.o sol_lusgs_scmp.o \
    sol_prmat.o sol_prsca.o sol_prsca_scmp.o \
    tur_bc.o turbulent.o tur_menter.o tur_spalart.o tur_lusgs.o \
    solve.o main.o

echo "=== Build complete: $BDIR/HOSTA.probe ==="
