
subroutine rhs_source
    use mod_kndconsts, only : kind_int
    use mod_constants, only : bc_mark_sponge
    use mod_datatypes, only : top_block_t
    use mod_variables, only : nsponge
    use mod_fieldvars, only : nblkcoms,blkcoms
    implicit none
    integer(kind_int)          :: nc,nb
    type(top_block_t), pointer :: top

    if (nsponge > 0) then
        do nc=1,nblkcoms
            nb  =  blkcoms(nc)%nb
            top => blkcoms(nc)%top
            if ( top%mark == bc_mark_sponge ) then
                call sponge_source(nb)
            end if
            !!call acoustic_source(nb)
        end do
    end if

end subroutine rhs_source

subroutine rhs_source_sp
    use mod_kndconsts, only : kind_int
    use mod_constants, only : bc_mark_sponge
    use mod_datatypes, only : top_block_t
    use mod_variables, only : nsponge
    use mod_fieldvars, only : nblkcoms,blkcomssp
    implicit none
    integer(kind_int)          :: nc,nb
    type(top_block_t), pointer :: top

    if (nsponge > 0) then
        do nc=1,nblkcoms
            nb  =  blkcomssp(nc)%nb
            top => blkcomssp(nc)%top
            if ( top%mark == bc_mark_sponge ) then
                call sponge_source_sp(nb)
            end if
            !!call acoustic_source(nb)
        end do
    end if

end subroutine rhs_source_sp

subroutine sponge_source(nb)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nsw_dir_close
    use mod_datatypes, only : fld_array_t
    use mod_variables, only : sigspg,nsw_kdir
    use mod_variables, only : roo,uoo,voo,woo,poo
    use mod_fieldvars, only : mb_top,mb_vol,mb_qc
    use mod_fieldvars, only : mb_dsp,npvs,neqn,mb_rhs
    implicit none
    integer(kind_int), intent(in) :: nb
    integer(kind_int)          :: i,j,k,ni,nj,nk,m,st(3),ed(3),nkst
    real(kind_real)            :: vol0,disp,alpha,sigma,rhs0
    real(kind_real)            :: pvoo(npvs),qcoo(neqn)
    type(fld_array_t), pointer :: vol(:),dsp(:)
    type(fld_array_t), pointer :: qc(:),rhs(:)

    alpha = 2.0

    pvoo(:) = (/roo,uoo,voo,woo,poo/)
    call prim2con(1,npvs,pvoo,1,neqn,qcoo)

    vol => mb_vol(nb)%fld
    dsp => mb_dsp(nb)%fld
    qc  => mb_qc(nb)%fld
    rhs => mb_rhs(nb)%fld

    ni = mb_top(nb)%nijk(1)
    nj = mb_top(nb)%nijk(2)
    nk = mb_top(nb)%nijk(3)

    st(:) = mb_top(nb)%ndst(:)
    ed(:) = mb_top(nb)%nded(:)

    do k=st(3),ed(3)
    do j=st(2),ed(2)
    do i=st(1),ed(1)
        vol0 = vol(1)%r3d(i,j,k)
        disp = dsp(1)%r3d(i,j,k)
        sigma = sigspg*disp**alpha
        do m=1,neqn
            rhs0 = sigma*(qcoo(m) - qc(m)%r3d(i,j,k))
            rhs(m)%r3d(i,j,k) = rhs(m)%r3d(i,j,k) + rhs0*vol0
        end do
    end do
    end do
    end do

    if (nsw_kdir == nsw_dir_close) then
        nkst = mb_top(nb)%ndst(3)

        do m=1,neqn
            do k=1,nk
            do j=1,nj
            do i=1,ni
                rhs(m)%r3d(i,j,k) = rhs(m)%r3d(i,j,nkst)
            end do
            end do
            end do
        end do
    end if

end subroutine sponge_source

subroutine sponge_source_sp(nb)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nsw_dir_close
    use mod_datatypes, only : fld_array_t
    use mod_variables, only : sigspg,nsw_kdir
    use mod_variables, only : roo,uoo,voo,woo,poo
    use mod_fieldvars, only : mb_topsp,mb_volsp,mb_qc
    use mod_fieldvars, only : mb_dsp,npvs,neqn,mb_rhs
    implicit none
    integer(kind_int), intent(in) :: nb
    integer(kind_int)          :: i,j,k,ni,nj,nk,m,st(3),ed(3),nkst
    real(kind_real)            :: vol0,disp,alpha,sigma,rhs0
    real(kind_real)            :: pvoo(npvs),qcoo(neqn)
    type(fld_array_t), pointer :: vol(:),dsp(:)
    type(fld_array_t), pointer :: qc(:),rhs(:)

    alpha = 2.0

    pvoo(:) = (/roo,uoo,voo,woo,poo/)
    call prim2con(1,npvs,pvoo,1,neqn,qcoo)

    vol => mb_volsp(nb)%fld
    dsp => mb_dsp(nb)%fld
    qc  => mb_qc(nb)%fld
    rhs => mb_rhs(nb)%fld

    ni = mb_topsp(nb)%nijk(1)
    nj = mb_topsp(nb)%nijk(2)
    nk = mb_topsp(nb)%nijk(3)

    st(:) = mb_topsp(nb)%ndst(:)
    ed(:) = mb_topsp(nb)%nded(:)

    do k=st(3),ed(3)
    do j=st(2),ed(2)
    do i=st(1),ed(1)
        vol0 = vol(1)%r3d(i,j,k)
        disp = dsp(1)%r3d(i,j,k)
        sigma = sigspg*disp**alpha
        do m=1,neqn
            rhs0 = sigma*(qcoo(m) - qc(m)%r3d(i,j,k))
            rhs(m)%r3d(i,j,k) = rhs(m)%r3d(i,j,k) + rhs0*vol0
        end do
    end do
    end do
    end do

    if (nsw_kdir == nsw_dir_close) then
        nkst = mb_topsp(nb)%ndst(3)

        do m=1,neqn
            do k=1,nk
            do j=1,nj
            do i=1,ni
                rhs(m)%r3d(i,j,k) = rhs(m)%r3d(i,j,nkst)
            end do
            end do
            end do
        end do
    end if

end subroutine sponge_source_sp

subroutine acoustic_source(nb)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nsw_dir_close,eight
    use mod_constants, only : pai,one,two,four,six
    use mod_datatypes, only : fld_array_t
    use mod_variables, only : dtau,nstep,nsw_kdir
    use mod_fieldvars, only : mb_xyz,mb_vol,mb_top
    use mod_fieldvars, only : mb_rhs,neqn
    implicit none
    integer(kind_int), intent(in) :: nb
    integer(kind_int)          :: i,j,k,ni,nj,nk,m,nkst
    real(kind_real)             :: vol0,rhs0,ds
    type(fld_array_t), pointer :: vol(:),rhs(:)
    type(fld_array_t), pointer :: xyz(:)


    vol => mb_vol(nb)%fld
    rhs => mb_rhs(nb)%fld
    xyz => mb_xyz(nb)%fld

    ni = mb_top(nb)%nijk(1)
    nj = mb_top(nb)%nijk(2)
    nk = mb_top(nb)%nijk(3)


    do k=1,nk
    do j=1,nj
    do i=1,ni
        vol0 = vol(1)%r3d(i,j,k)
        ds=xyz(1)%r3d(i,j,k)**two+xyz(2)%r3d(i,j,k)**two
        if(dtau*(nstep*one)<four) then
            rhs0 = 0.4*sin(pai/two*dtau*(nstep*one)/six)**two*exp(-25.0*log(two)*ds) &
                                  *sin(eight*pai*dtau*(nstep*one))
            rhs(5)%r3d(i,j,k) = rhs(5)%r3d(i,j,k) + rhs0*vol0
        else
            rhs0 = 0.4*exp(-25.0*log(two)*ds)*sin(eight*pai*dtau*(nstep*one))
            rhs(5)%r3d(i,j,k) = rhs(5)%r3d(i,j,k) + rhs0*vol0
        end if
    end do
    end do
    end do

    if (nsw_kdir == nsw_dir_close) then
        nkst = mb_top(nb)%ndst(3)

        do m=1,neqn
            do k=1,nk
            do j=1,nj
            do i=1,ni
                rhs(m)%r3d(i,j,k) = rhs(m)%r3d(i,j,nkst)
            end do
            end do
            end do
        end do
    end if

end subroutine acoustic_source

!!subroutine channel_source
!!    use mod_kndconsts, only : kind_int,kind_real
!!    use mod_constants, only : nsw_dir_close,eight
!!    use mod_constants, only : pai,one,two,four,six
!!    use mod_datatypes, only : fld_array_t
!!    use mod_variables, only : dtau,nstep,nsw_kdir,reue
!!    use mod_fieldvars, only : nblkcoms,blkcoms
!!    use mod_fieldvars, only : mb_sxyz,mb_vol,mb_top,neqn
!!    use mod_fieldvars, only : mb_rhs,mb_dpv,mb_pv,mb_vsl
!!    use mod_parallels
!!    implicit none
!!    integer(kind_int)          :: i,j,k,ni,nj,nk,m,nkst,nb,nc,ierr
!!    real(kind_real)             :: rhsc,rhsc1,rhsc2,ds,vis0,visn
!!    real(kind_real)             :: dpv0,dpvn,sxyz0,sxyzn,volo,voln
!!    type(fld_array_t), pointer :: vol(:),rhs(:),vsl(:)
!!    type(fld_array_t), pointer :: sxyz(:),dpv(:),pv(:)
!!
!!
!!!!        do nc=1,nblkcoms
!!!!            nb  =  blkcoms(nc)%nb
!!!!
!!!!            if (nb == 1) then
!!!!                dpv => mb_dpv(1)%fld
!!!!                vol => mb_vol(1)%fld
!!!!                vsl => mb_vsl(1)%fld
!!!!                sxyz => mb_sxyz(1)%fld
!!!!
!!!!                nj0  = mb_top(1)%nijk(2)
!!!!                volo = vol(1)%r3d(1,1,1)
!!!!                voln = vol(1)%r3d(1,nj0,1)
!!!!                vis0 = vsl(1)%r3d(1,1,1)
!!!!                visn = vsl(1)%r3d(1,nj0,1)
!!!!                dpv0 = dpv(2)%r3d(1,1,1)
!!!!                dpvn = dpv(2)%r3d(1,nj0,1)
!!!!                sxyz0 = sxyz(5)%r3d(1,1,1)
!!!!                sxyzn = sxyz(5)%r3d(1,nj0,1)
!!!!
!!!!                do j=1,nj0
!!!!                    rhsc = rhsc + vol(1)%r3d(1,j,1)
!!!!                end do
!!!!                rhsc = 1.0/(reue*rhsc)*(visn*sxyzn*dpvn-vis0*sxyz0*dpv0)
!!!!            end if
!!!!        end do
!!
!!        do nc=1,nblkcoms
!!
!!            nb   =  blkcoms(nc)%nb
!!
!!            rhs  => mb_rhs(nb)%fld
!!            pv   => mb_pv(nb)%fld
!!            vol  => mb_vol(nb)%fld
!!            dpv  => mb_dpv(nb)%fld
!!            vsl  => mb_vsl(nb)%fld
!!            sxyz => mb_sxyz(nb)%fld
!!
!!            ni = mb_top(nb)%nijk(1)
!!            nj = mb_top(nb)%nijk(2)
!!            nk = mb_top(nb)%nijk(3)
!!
!!            volo  = vol(1)%r3d(1,1,1)
!!            voln  = vol(1)%r3d(1,nj,1)
!!            vis0  = vsl(1)%r3d(1,1,1)
!!            visn  = vsl(1)%r3d(1,nj,1)
!!            dpv0  = dpv(2)%r3d(1,1,1)
!!            dpvn  = dpv(2)%r3d(1,nj,1)
!!            sxyz0 = sxyz(5)%r3d(1,1,1)
!!            sxyzn = sxyz(5)%r3d(1,nj,1)
!!
!!            rhsc  = 0.0
!!            rhsc1 = 0.0
!!
!!            do j=1,nj
!!                rhsc  = rhsc + vol(1)%r3d(1,j,1)
!!            end do
!!
!!
!!            do k=1,nk
!!            do i=1,ni
!!                volo  = vol(1)%r3d(i,1,k)
!!                voln  = vol(1)%r3d(i,nj,k)
!!                vis0  = vsl(1)%r3d(i,1,k)
!!                visn  = vsl(1)%r3d(i,nj,k)
!!                dpv0  = dpv(2)%r3d(i,1,k)
!!                dpvn  = dpv(2)%r3d(i,nj,k)
!!                sxyz0 = sxyz(5)%r3d(i,1,k)
!!                sxyzn = sxyz(5)%r3d(i,nj,k)
!!                rhsc1 = rhsc1 + visn*sxyzn*dpvn-vis0*sxyz0*dpv0
!!            end do
!!            end do
!!
!!            rhsc1 = 1.0/(reue*rhsc)*rhsc1/(ni*nk*1.0)
!!        end do
!!
!!#ifdef PARALLEL
!!            call MPI_REDUCE(rhsc1,rhsc2,1,kind_real_mpi,MPI_SUM,master,MPI_COMM_WORLD,ierr)
!!            call MPI_BCAST(rhsc2,1,kind_real_mpi,master,MPI_COMM_WORLD,ierr)
!!#endif
!!
!!#ifdef PARALLEL
!!           if (myid == master) then
!!#endif
!!            write(*,*) rhsc2/4.0
!!#ifdef PARALLEL
!!           end if
!!#endif
!!
!!        do nc=1,nblkcoms
!!            do k=1,nk
!!            do j=1,nj
!!            do i=1,ni
!!                rhs(2)%r3d(i,j,k) = rhs(2)%r3d(i,j,k) - rhsc2*vol(1)%r3d(i,j,k)/4.0
!!                rhs(5)%r3d(i,j,k) = rhs(5)%r3d(i,j,k) - rhsc2*vol(1)%r3d(i,j,k)*pv(2)%r3d(i,j,k)/4.0
!!            end do
!!            end do
!!            end do
!!
!!!!            do k=1,nk
!!!!            do j=1,nj
!!!!            do i=1,ni
!!!!                rhsc  = 0.0
!!!!                do jj=1,nj
!!!!                    rhsc = rhsc + vol(1)%r3d(i,jj,k)
!!!!                end do
!!!!                rhsc = 1.0/(reue*rhsc)*(visn*sxyzn*dpvn-vis0*sxyz0*dpv0)
!!!!
!!!!                rhs(2)%r3d(i,j,k) = rhs(2)%r3d(i,j,k) - rhsc*vol(1)%r3d(i,j,k)
!!!!                rhs(5)%r3d(i,j,k) = rhs(5)%r3d(i,j,k) - rhsc*vol(1)%r3d(i,j,k)*pv(2)%r3d(i,j,k)
!!!!            end do
!!!!            end do
!!!!            end do
!!
!!!!            if (nb == 4) then
!!!!                do k=1,nk
!!!!                do j=1,nj
!!!!                do i=ni,ni
!!!!                    rhsc  = 0.0
!!!!                    do jj=1,nj
!!!!                        rhsc = rhsc + vol(1)%r3d(i,jj,k)
!!!!                    end do
!!!!                    rhsc = 1.0/(reue*rhsc*0.04)*(visn*dpvn-vis0*dpv0)
!!!!
!!!!                    rhs(2)%r3d(i,j,k) = rhs(2)%r3d(i,j,k) - rhsc*vol(1)%r3d(i,j,k)
!!!!                    rhs(5)%r3d(i,j,k) = rhs(5)%r3d(i,j,k) - rhsc*vol(1)%r3d(i,j,k)*pv(2)%r3d(i,j,k)
!!!!                end do
!!!!                end do
!!!!                end do
!!!!            end if
!!
!!        end do
!!
!!end subroutine channel_source


subroutine channel_source
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : bc_wall,nsw_dir_close,eight
    use mod_constants, only : pai,one,two,four,six
    use mod_datatypes, only : var_block_t,fld_array_t
    use mod_datatypes, only : top_block_t,bc_region_t
    use mod_variables, only : dtau,nstep,nsw_kdir,reue
    use mod_fieldvars, only : nblkcoms,blkcoms
    use mod_fieldvars, only : mb_sxyz,mb_vol,mb_top,neqn
    use mod_fieldvars, only : mb_rhs,mb_dpv,mb_pv,mb_vsl
    use mod_parallels
    implicit none
    integer(kind_int)          :: ns,nr,num2,num3
    integer(kind_int)          :: s_st(3),s_ed(3),st(3),ed(3)
    integer(kind_int)          :: nregs,s_nd,s_lr,bctype
    integer(kind_int)          :: i,j,k,ni,nj,nk,m,nkst,nb,nc,ierr
    real(kind_real)             :: rhsc1,rhsc2,ds
    real(kind_real)             :: dpv0,sxyz0,vis0
    type(top_block_t), pointer :: top
    type(fld_array_t), pointer :: vol(:),rhs(:),vsl(:)
    type(fld_array_t), pointer :: sxyz(:),dpv(:),pv(:)
    type(bc_region_t), pointer :: reg

    rhsc2 = 8.22467e-3 !9.1385225936e-4 !2.05616758e-3 !3.65540903744e-3
    rhsc1 = 0.0
    num2  = 0
    do nc=1,nblkcoms

        nb   =  blkcoms(nc)%nb
        top => blkcoms(nc)%top

        rhs  => mb_rhs(nb)%fld
        pv   => mb_pv(nb)%fld
        vol  => mb_vol(nb)%fld
        dpv  => mb_dpv(nb)%fld
        vsl  => mb_vsl(nb)%fld
        sxyz => mb_sxyz(nb)%fld

        nregs = top%nregions
        do ns=1,nregs
            nr  =  top%bcnrs(ns)
            reg => top%bcs(nr)
            bctype  = reg%bctype
            if (bctype == bc_wall) then
                s_st(:) = reg%s_st(:)
                s_ed(:) = reg%s_ed(:)
                s_nd    = reg%s_nd
                s_lr    = reg%s_lr

                st(:) = s_st(:)
                ed(:) = s_ed(:)

                do i=st(1),ed(1)
                    do j=st(2),ed(2)
                       do k=st(3),ed(3)
                           vis0  = vsl(1)%r3d(i,j,k)
                           dpv0  = dpv(2)%r3d(i,j,k)
                           sxyz0 = sxyz(5)%r3d(i,j,k)
                           rhsc1 = rhsc1 - abs(vis0*sxyz0*dpv0)
                           num2  = num2 + 1
                       end do
                   end do
                end do
            end if
        end do
    end do
    rhsc1 = 1.0/(reue*rhsc2)*rhsc1

    num3  = num2
    rhsc2 = rhsc1

#ifdef PARALLEL
    call MPI_REDUCE(num2,num3,1,kind_int_mpi,MPI_SUM,master,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(num3,1,kind_int_mpi,master,MPI_COMM_WORLD,ierr)
    call MPI_REDUCE(rhsc1,rhsc2,1,kind_real_mpi,MPI_SUM,master,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(rhsc2,1,kind_real_mpi,master,MPI_COMM_WORLD,ierr)
#endif

#ifdef PARALLEL
           if (myid == master) then
#endif
            write(*,*) rhsc2/( num3*0.5 )
#ifdef PARALLEL
           end if
#endif

        do nc=1,nblkcoms
            nb  =  blkcoms(nc)%nb
            top => blkcoms(nc)%top

            ni = top%nijk(1)
            nj = top%nijk(2)
            nk = top%nijk(3)

            rhs => mb_rhs(nb)%fld
            pv  => mb_pv(nb)%fld
            vol => mb_vol(nb)%fld

            do k=1,nk
            do j=1,nj
            do i=1,ni
                rhs(2)%r3d(i,j,k) = rhs(2)%r3d(i,j,k) - rhsc2*vol(1)%r3d(i,j,k)/( num3*0.5 )
                rhs(5)%r3d(i,j,k) = rhs(5)%r3d(i,j,k) - rhsc2*vol(1)%r3d(i,j,k)*pv(2)%r3d(i,j,k)/( num3*0.5 )
            end do
            end do
            end do
        end do

end subroutine channel_source

