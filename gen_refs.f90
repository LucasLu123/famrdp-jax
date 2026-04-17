! gen_refs.f90 — standalone program to generate validation reference binaries
! Reads test3.grd (PLOT3D binary stream) and writes probe-format reference files:
!   validation/references/test3/xyz/block_NNNN.bin
!   validation/references/test3/metric/block_NNNN_jac.bin
!
! Binary format (little-endian, matches read_probe in fortran_ref.py):
!   xyz: int32[3] (ni,nj,nk), float64[ni*nj*nk] x, float64 y, float64 z
!   jac: int32[3] (ni,nj,nk), float64[ni*nj*nk] jac (Fortran order)

program gen_refs
    implicit none
    integer, parameter :: r8 = 8
    integer, parameter :: i4 = 4

    integer :: nb, nblocks, ni, nj, nk, ios
    integer, allocatable :: dims(:,:)
    real(r8), allocatable :: x(:,:,:), y(:,:,:), z(:,:,:)
    real(r8), allocatable :: jac(:,:,:), kxyz(:,:,:,:,:)
    character(len=256) :: fname
    integer :: uin, uout

    ! ---- Read PLOT3D grid ----
    call openfile_stream(uin, "test3.grd")
    read(uin) nblocks
    write(*,'(A,I0,A)') "[gen_refs] nblocks = ", nblocks

    allocate(dims(3, nblocks))
    do nb = 1, nblocks
        read(uin) dims(1,nb), dims(2,nb), dims(3,nb)
    end do

    do nb = 1, nblocks
        ni = dims(1,nb); nj = dims(2,nb); nk = dims(3,nb)
        allocate(x(ni,nj,nk), y(ni,nj,nk), z(ni,nj,nk))
        allocate(jac(ni,nj,nk), kxyz(3,3,ni,nj,nk))

        read(uin) x, y, z

        ! ---- Write xyz reference ----
        write(fname,'(A,I4.4,A)') "validation/references/test3/xyz/block_", nb, ".bin"
        open(newunit=uout, file=trim(fname), form="unformatted", access="stream", iostat=ios)
        if (ios /= 0) then
            write(*,'(A,A)') "[gen_refs] ERROR opening ", trim(fname)
            stop 1
        end if
        write(uout) ni, nj, nk
        write(uout) x
        write(uout) y
        write(uout) z
        close(uout)
        write(*,'(A,A)') "[gen_refs] wrote ", trim(fname)

        ! ---- Compute metrics (4th-order central, matches Fortran ncutpol=3) ----
        call compute_jac(x, y, z, ni, nj, nk, jac, kxyz)

        ! ---- Write jac reference ----
        write(fname,'(A,I4.4,A)') "validation/references/test3/metric/block_", nb, "_jac.bin"
        open(newunit=uout, file=trim(fname), form="unformatted", access="stream", iostat=ios)
        if (ios /= 0) then
            write(*,'(A,A)') "[gen_refs] ERROR opening ", trim(fname)
            stop 1
        end if
        write(uout) ni, nj, nk
        write(uout) jac
        close(uout)
        write(*,'(A,A)') "[gen_refs] wrote ", trim(fname)

        deallocate(x, y, z, jac, kxyz)
    end do

    close(uin)
    write(*,*) "[gen_refs] done."

contains

    subroutine openfile_stream(u, fname)
        integer, intent(out) :: u
        character(len=*), intent(in) :: fname
        integer :: ios
        open(newunit=u, file=fname, form="unformatted", access="stream", status="old", iostat=ios)
        if (ios /= 0) then
            write(*,'(A,A)') "ERROR: cannot open ", fname
            stop 1
        end if
    end subroutine

    ! 4th-order central difference: df/dxi along axis dir (1=i,2=j,3=k)
    ! Uses periodic-wrap at boundaries (matches JAX jnp.roll approach)
    ! 4th-order central difference matching JAX jnp.roll(-2)/roll(-1)/roll(+1)/roll(+2)
    ! JAX roll(-k) at 0-indexed pos p gives element at (p+k) mod n
    ! Equivalent Fortran 1-indexed: for pos i, fp2 = f((i+1) mod ni + 1), etc.
    subroutine central4(f, ni, nj, nk, dir, df)
        real(r8), intent(in)  :: f(ni,nj,nk)
        integer,  intent(in)  :: ni, nj, nk, dir
        real(r8), intent(out) :: df(ni,nj,nk)
        real(r8) :: fp2(ni,nj,nk), fp1(ni,nj,nk), fm1(ni,nj,nk), fm2(ni,nj,nk)
        integer :: i

        select case(dir)
        case(1)  ! i direction
            do i = 1, ni
                fp2(i,:,:) = f(mod(i+1,  ni)+1, :, :)
                fp1(i,:,:) = f(mod(i,    ni)+1, :, :)
                fm1(i,:,:) = f(mod(i+ni-2,ni)+1, :, :)
                fm2(i,:,:) = f(mod(i+ni-3,ni)+1, :, :)
            end do
        case(2)  ! j direction
            do i = 1, nj
                fp2(:,i,:) = f(:, mod(i+1,  nj)+1, :)
                fp1(:,i,:) = f(:, mod(i,    nj)+1, :)
                fm1(:,i,:) = f(:, mod(i+nj-2,nj)+1, :)
                fm2(:,i,:) = f(:, mod(i+nj-3,nj)+1, :)
            end do
        case(3)  ! k direction
            do i = 1, nk
                fp2(:,:,i) = f(:, :, mod(i+1,  nk)+1)
                fp1(:,:,i) = f(:, :, mod(i,    nk)+1)
                fm1(:,:,i) = f(:, :, mod(i+nk-2,nk)+1)
                fm2(:,:,i) = f(:, :, mod(i+nk-3,nk)+1)
            end do
        end select
        df = (-fp2 + 8.0_r8*fp1 - 8.0_r8*fm1 + fm2) / 12.0_r8
    end subroutine

    subroutine compute_jac(x, y, z, ni, nj, nk, jac, kxyz)
        integer, intent(in) :: ni, nj, nk
        real(r8), intent(in)  :: x(ni,nj,nk), y(ni,nj,nk), z(ni,nj,nk)
        real(r8), intent(out) :: jac(ni,nj,nk)
        real(r8), intent(out) :: kxyz(3,3,ni,nj,nk)
        real(r8) :: dxdxi(ni,nj,nk), dxdet(ni,nj,nk), dxdct(ni,nj,nk)
        real(r8) :: dydxi(ni,nj,nk), dydet(ni,nj,nk), dydct(ni,nj,nk)
        real(r8) :: dzdxi(ni,nj,nk), dzdet(ni,nj,nk), dzdct(ni,nj,nk)
        ! J_mat[m,l] = dx_m/dxi_l; det(J)
        call central4(x, ni, nj, nk, 1, dxdxi)
        call central4(x, ni, nj, nk, 2, dxdet)
        call central4(x, ni, nj, nk, 3, dxdct)
        call central4(y, ni, nj, nk, 1, dydxi)
        call central4(y, ni, nj, nk, 2, dydet)
        call central4(y, ni, nj, nk, 3, dydct)
        call central4(z, ni, nj, nk, 1, dzdxi)
        call central4(z, ni, nj, nk, 2, dzdet)
        call central4(z, ni, nj, nk, 3, dzdct)
        ! det(J_mat) where J_mat = [[dxdxi,dxdet,dxdct],[dydxi,dydet,dydct],[dzdxi,dzdet,dzdct]]
        jac = dxdxi*(dydet*dzdct - dydct*dzdet) &
            - dxdet*(dydxi*dzdct - dydct*dzdxi) &
            + dxdct*(dydxi*dzdet - dydet*dzdxi)
    end subroutine

end program gen_refs
