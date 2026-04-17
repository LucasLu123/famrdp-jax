#include "probes.h"

module probe_utils
    implicit none
    private
    public :: probe_write_metric, probe_write_rhs, probe_write_step, probe_write_xyz
    public :: probe_write_jac

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

    subroutine probe_write_xyz(nb, x, y, z)
        integer, intent(in) :: nb
        real(8), intent(in) :: x(:,:,:)
        real(8), intent(in) :: y(:,:,:)
        real(8), intent(in) :: z(:,:,:)
        character(len=256) :: fname
        integer :: unit, ios
        integer :: ni, nj, nk
        ni = size(x,1); nj = size(x,2); nk = size(x,3)
        write(fname,'(A,I4.4,A)') "probes/xyz/block_", nb, ".bin"
        write(0,'(A,A)') "[probe] writing ", trim(fname)
        open(newunit=unit, file=trim(fname), form="unformatted", access="stream", iostat=ios)
        if (ios /= 0) then
            write(0,'(A,I0)') "[probe] OPEN ERROR iostat=", ios
            return
        end if
        write(unit) ni, nj, nk
        write(unit) x
        write(unit) y
        write(unit) z
        close(unit)
        write(0,'(A)') "[probe] xyz written"
    end subroutine


    subroutine probe_write_jac(nb, jac)
        integer, intent(in) :: nb
        real(8), intent(in) :: jac(:,:,:)
        character(len=256) :: fname
        integer :: unit
        integer :: ishape(3)
        ishape = shape(jac)
        write(fname,'(A,I4.4,A)') "probes/metric/block_", nb, "_jac.bin"
        open(newunit=unit, file=trim(fname), form="unformatted", access="stream")
        write(unit) ishape
        write(unit) jac
        close(unit)
    end subroutine

end module probe_utils
