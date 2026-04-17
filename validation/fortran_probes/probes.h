// validation/fortran_probes/probes.h
#ifndef PROBES_H
#define PROBES_H

#ifdef PROBE_METRIC
#define PROBE_DUMP_METRIC(nb, jac, kxyz, vol) \
    call probe_write_metric(nb, jac, kxyz, vol)
#else
#define PROBE_DUMP_METRIC(nb, jac, kxyz, vol)
#endif

#ifdef PROBE_RHS
#define PROBE_DUMP_RHS(nb, rhs) call probe_write_rhs(nb, rhs)
#else
#define PROBE_DUMP_RHS(nb, rhs)
#endif

#ifdef PROBE_STEP
#define PROBE_DUMP_STEP(nstep) call probe_write_step(nstep)
#else
#define PROBE_DUMP_STEP(nstep)
#endif

#endif
