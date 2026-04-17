#include "probes.h"

module probe_utils
    implicit none
    private
    public :: probe_write_metric, probe_write_rhs, probe_write_step

contains

    subroutine probe_write_metric(nb, jac, kxyz, vol)
        integer, intent(in) :: nb
        real(8), intent(in) :: jac(:,:,:)
        real(8), intent(in) :: kxyz(:,:,:,:,:)
        real(8), intent(in) :: vol(:,:,:)
        character(len=256) :: fname
        integer :: unit
        write(fname,'(A,I4.4,A)') "probes/metric/block_", nb, ".bin"
        open(newunit=unit, file=trim(fname), form="unformatted", access="stream")
        write(unit) shape(jac)
        write(unit) jac
        write(unit) shape(kxyz)
        write(unit) kxyz
        write(unit) shape(vol)
        write(unit) vol
        close(unit)
    end subroutine

    subroutine probe_write_rhs(nb, rhs)
        integer, intent(in) :: nb
        real(8), intent(in) :: rhs(:,:,:,:)
        character(len=256) :: fname
        integer :: unit
        write(fname,'(A,I4.4,A)') "probes/rhs/block_", nb, ".bin"
        open(newunit=unit, file=trim(fname), form="unformatted", access="stream")
        write(unit) shape(rhs)
        write(unit) rhs
        close(unit)
    end subroutine

    subroutine probe_write_step(nstep_in)
        integer, intent(in) :: nstep_in
        character(len=256) :: fname
        integer :: unit
        write(fname,'(A,I6.6,A)') "probes/step/step_", nstep_in, ".bin"
        open(newunit=unit, file=trim(fname), form="unformatted", access="stream")
        write(unit) nstep_in
        close(unit)
    end subroutine

end module probe_utils
