
subroutine sol_lusgs_std_sca
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero
    use mod_constants, only : ntimeadv_steady
    use mod_variables, only : nghnode,nstep
    use mod_fieldvars, only : npvs,mb_pv
    use mod_fieldvars, only : nblkcoms,blkcoms
    use mod_fieldvars, only : neqn,mb_dq
    use mod_interface, only : assign_mb_var_uniform,assign_bc_var_nb
    use mod_interface, only : assign_var_via_com_nb
    implicit none
    integer(kind_int) :: nc,nb
    real(kind_real)   :: dq0(1:neqn)
    
    call prepare_linear_system

    dq0(:) = zero
    call assign_mb_var_uniform(mb_dq,1,neqn,nghnode,dq0)

    do nc=1,nblkcoms
        nb = blkcoms(nc)%nb

        call lusgs_std_sca_foreward(nb)
        call assign_bc_var_nb(nb,mb_dq,1,neqn,0,dq0)

        call lusgs_std_sca_backward(nb)
        call assign_bc_var_nb(nb,mb_dq,1,neqn,0,dq0)

        call assign_var_via_com_nb(nb,mb_dq,1,neqn)
    end do
    
    call residual(ntimeadv_steady)

    call update

end subroutine sol_lusgs_std_sca

subroutine sol_lusgs_std_sca_sp
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero
    use mod_constants, only : ntimeadv_steady
    use mod_variables, only : nghnode,nstep
    use mod_fieldvars, only : npvs,mb_pv
    use mod_fieldvars, only : nblkcoms,blkcomssp
    use mod_fieldvars, only : neqn,mb_dq
    use mod_interface, only : assign_mb_var_uniform_sp
    use mod_interface, only : assign_var_via_com_nb_sp
    implicit none
    integer(kind_int) :: nc,nb
    real(kind_real)   :: dq0(1:neqn)
    
    call prepare_linear_system_sp

    dq0(:) = zero
    call assign_mb_var_uniform_sp(mb_dq,1,neqn,nghnode,dq0)

    do nc=1,nblkcoms
        nb = blkcomssp(nc)%nb

        call lusgs_std_sca_foreward_sp(nb)
        !call assign_bc_var_nb(nb,mb_dq,1,neqn,0,dq0)

        call lusgs_std_sca_backward_sp(nb)
        !call assign_bc_var_nb(nb,mb_dq,1,neqn,0,dq0)

        call assign_var_via_com_nb_sp(nb,mb_dq,1,neqn)
    end do
    
    call residual_sp(ntimeadv_steady)

    call update_sp

end subroutine sol_lusgs_std_sca_sp

subroutine lusgs_std_sca_foreward(nb)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,one,half
    use mod_constants, only : nvis_euler,nsw_dir_close
    use mod_datatypes, only : fld_array_t
    use mod_variables, only : nvis,nsw_kdir
    use mod_fieldvars, only : mb_top,mb_sxyz,mb_pv
    use mod_fieldvars, only : mb_dq,npvs,neqn,mb_rhs
    use mod_fieldvars, only : mb_src,mb_srv,mb_rdt
    use mod_openmp
    implicit none
    integer(kind_int), intent(in) :: nb
    integer(kind_int)          :: i,j,k,i1,j1,k1,m,ierr
    integer(kind_int)          :: st(3),ed(3),dijk
    real(kind_real)            :: srci,srcj,srck,srvi,srvj,srvk
    real(kind_real)            :: kt,kx,ky,kz,et,ex,ey,ez,ct,cx,cy,cz
    real(kind_real)            :: odia,rhs0(1:neqn)
    real(kind_real)            :: pvi(1:npvs),pvj(1:npvs),pvk(1:npvs)
    real(kind_real)            :: dqi(1:neqn),dqj(1:neqn),dqk(1:neqn)
    real(kind_real)            :: dfi(1:neqn),dfj(1:neqn),dfk(1:neqn)
    type(fld_array_t), pointer :: sxyz(:),pv(:),rhs(:)
    type(fld_array_t), pointer :: dq(:),src(:),srv(:),rdt(:)
#ifdef OMP_IMP
    integer(kind_int)          :: idt,numt,nflag(0:omp_max_num_threads)
#endif

    sxyz => mb_sxyz(nb)%fld
    pv   => mb_pv(nb)%fld
    rhs  => mb_rhs(nb)%fld
    dq   => mb_dq(nb)%fld
    src  => mb_src(nb)%fld
    rdt  => mb_rdt(nb)%fld

    if (nvis > nvis_euler) then
        srv  => mb_srv(nb)%fld
    end if

    st(:) = mb_top(nb)%ndst(:)
    ed(:) = mb_top(nb)%nded(:)
    dijk = 1

!$OMP parallel private(k,k1,j,j1,i,i1,idt,pvi,pvj,pvk,dqi,dqj,dqk,rhs0,dfi,dfj,dfk,&
!$OMP                  kx,ky,kz,ex,ey,ez,cx,cy,cz,srvi,srvj,srvk,srci,srcj,srck,odia)

#ifdef OMP_IMP
!$OMP master
    numt = omp_get_num_threads()
    if (numt > omp_max_num_threads) then
        write(*,*) "OpenMP num_threads > omp_max_num_threads, make sure you want do it, check size of nflag"
    end if
!$OMP end master
    idt = omp_get_thread_num()
    nflag(idt) = 0
#endif

!$OMP barrier

    do k=st(3),ed(3),dijk

#ifdef OMP_IMP
    if (idt > 0 .and. idt < numt) then
        do while (nflag(idt-1) == 0)
            !$OMP flush(nflag)
        end do
        nflag(idt-1) = 0
        !$OMP flush(nflag)
    end if
#endif

!$OMP do schedule(static)
    do j=st(2),ed(2),dijk
    do i=st(1),ed(1),dijk
        i1 = i - 1
        j1 = j - 1
        k1 = k - 1

        srci = src(1)%r3d(i1,j,k)!ksa方向对流项谱半径
        srcj = src(2)%r3d(i,j1,k)
        srck = src(3)%r3d(i,j,k1)

        do m=1,npvs
            pvi(m) = pv(m)%r3d(i1,j,k)!原始变量
            pvj(m) = pv(m)%r3d(i,j1,k)
            pvk(m) = pv(m)%r3d(i,j,k1)
        end do

        do m=1,neqn
            dqi(m) = dq(m)%r3d(i1,j,k)
            dqj(m) = dq(m)%r3d(i,j1,k)
            dqk(m) = dq(m)%r3d(i,j,k1)
        end do

        kx = sxyz(1)%r3d(i1,j,k)!网格导数
        ky = sxyz(2)%r3d(i1,j,k)
        kz = sxyz(3)%r3d(i1,j,k)
        call mxdq_std(1,npvs,pvi,kt,kx,ky,kz,1,neqn,dqi,dfi,srci,one)

        ex = sxyz(4)%r3d(i,j1,k)
        ey = sxyz(5)%r3d(i,j1,k)
        ez = sxyz(6)%r3d(i,j1,k)
        call mxdq_std(1,npvs,pvj,et,ex,ey,ez,1,neqn,dqj,dfj,srcj,one)

        if (nsw_kdir == nsw_dir_close) then
            dfk(:) = zero
        else
            cx = sxyz(7)%r3d(i,j,k1)
            cy = sxyz(8)%r3d(i,j,k1)
            cz = sxyz(9)%r3d(i,j,k1)
            call mxdq_std(1,npvs,pvk,ct,cx,cy,cz,1,neqn,dqk,dfk,srck,one)
        end if

        rhs0(:) = dfi(:) + dfj(:) + dfk(:)


        if (nvis > nvis_euler) then
            srvi = srv(1)%r3d(i1,j,k)
            srvj = srv(2)%r3d(i,j1,k)
            srvk = srv(3)%r3d(i,j,k1)
            do m=1,neqn
                rhs0(m) = rhs0(m) + &
                          half*(srvi*dqi(m) + srvj*dqj(m) + srvk*dqk(m))
            end do
        end if

        odia = one/rdt(1)%r3d(i,j,k)
        do m=1,neqn
            dq(m)%r3d(i,j,k) = (rhs(m)%r3d(i,j,k) + rhs0(m))*odia
        end do
    end do
    end do
!$OMP end do nowait

#ifdef OMP_IMP
    if (idt < min(numt-1,ed(2)-st(2)+1)) then
        !$OMP flush(nflag)
        do while (nflag(idt) == 1)
            !$OMP flush(nflag)
        end do
        nflag(idt) = 1
        !$OMP flush(nflag)
    end if
#endif

    end do  ! K loop

!$OMP end parallel
end subroutine lusgs_std_sca_foreward

subroutine lusgs_std_sca_foreward_sp(nb)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,one,half
    use mod_constants, only : nvis_euler,nsw_dir_close
    use mod_datatypes, only : fld_array_t
    use mod_variables, only : nvis,nsw_kdir
    use mod_fieldvars, only : mb_topsp,mb_sxyzsp,mb_pv
    use mod_fieldvars, only : mb_dq,npvs,neqn,mb_rhs
    use mod_fieldvars, only : mb_src,mb_srv,mb_rdt
    use mod_openmp
    implicit none
    integer(kind_int), intent(in) :: nb
    integer(kind_int)          :: i,j,k,i1,j1,k1,m,ierr
    integer(kind_int)          :: st(3),ed(3),dijk
    real(kind_real)            :: srci,srcj,srck,srvi,srvj,srvk
    real(kind_real)            :: kt,kx,ky,kz,et,ex,ey,ez,ct,cx,cy,cz
    real(kind_real)            :: odia,rhs0(1:neqn)
    real(kind_real)            :: pvi(1:npvs),pvj(1:npvs),pvk(1:npvs)
    real(kind_real)            :: dqi(1:neqn),dqj(1:neqn),dqk(1:neqn)
    real(kind_real)            :: dfi(1:neqn),dfj(1:neqn),dfk(1:neqn)
    type(fld_array_t), pointer :: sxyz(:),pv(:),rhs(:)
    type(fld_array_t), pointer :: dq(:),src(:),srv(:),rdt(:)
#ifdef OMP_IMP
    integer(kind_int)          :: idt,numt,nflag(0:omp_max_num_threads)
#endif

    sxyz => mb_sxyzsp(nb)%fld
    pv   => mb_pv(nb)%fld
    rhs  => mb_rhs(nb)%fld
    dq   => mb_dq(nb)%fld
    src  => mb_src(nb)%fld
    rdt  => mb_rdt(nb)%fld

    if (nvis > nvis_euler) then
        srv  => mb_srv(nb)%fld
    end if

    st(:) = mb_topsp(nb)%ndst(:)
    ed(:) = mb_topsp(nb)%nded(:)
    dijk = 1

!$OMP parallel private(k,k1,j,j1,i,i1,idt,pvi,pvj,pvk,dqi,dqj,dqk,rhs0,dfi,dfj,dfk,&
!$OMP                  kx,ky,kz,ex,ey,ez,cx,cy,cz,srvi,srvj,srvk,srci,srcj,srck,odia)

#ifdef OMP_IMP
!$OMP master
    numt = omp_get_num_threads()
    if (numt > omp_max_num_threads) then
        write(*,*) "OpenMP num_threads > omp_max_num_threads, make sure you want do it, check size of nflag"
    end if
!$OMP end master
    idt = omp_get_thread_num()
    nflag(idt) = 0
#endif

!$OMP barrier

    do k=st(3),ed(3),dijk

#ifdef OMP_IMP
    if (idt > 0 .and. idt < numt) then
        do while (nflag(idt-1) == 0)
            !$OMP flush(nflag)
        end do
        nflag(idt-1) = 0
        !$OMP flush(nflag)
    end if
#endif

!$OMP do schedule(static)
    do j=st(2),ed(2),dijk
    do i=st(1),ed(1),dijk
        i1 = i - 1
        j1 = j - 1
        k1 = k - 1

        srci = src(1)%r3d(i1,j,k)!ksa方向对流项谱半径
        srcj = src(2)%r3d(i,j1,k)
        srck = src(3)%r3d(i,j,k1)

        do m=1,npvs
            pvi(m) = pv(m)%r3d(i1,j,k)!原始变量
            pvj(m) = pv(m)%r3d(i,j1,k)
            pvk(m) = pv(m)%r3d(i,j,k1)
        end do

        do m=1,neqn
            dqi(m) = dq(m)%r3d(i1,j,k)
            dqj(m) = dq(m)%r3d(i,j1,k)
            dqk(m) = dq(m)%r3d(i,j,k1)
        end do

        kx = sxyz(1)%r3d(i1,j,k)!网格导数
        ky = sxyz(2)%r3d(i1,j,k)
        kz = sxyz(3)%r3d(i1,j,k)
        call mxdq_std(1,npvs,pvi,kt,kx,ky,kz,1,neqn,dqi,dfi,srci,one)

        ex = sxyz(4)%r3d(i,j1,k)
        ey = sxyz(5)%r3d(i,j1,k)
        ez = sxyz(6)%r3d(i,j1,k)
        call mxdq_std(1,npvs,pvj,et,ex,ey,ez,1,neqn,dqj,dfj,srcj,one)

        if (nsw_kdir == nsw_dir_close) then
            dfk(:) = zero
        else
            cx = sxyz(7)%r3d(i,j,k1)
            cy = sxyz(8)%r3d(i,j,k1)
            cz = sxyz(9)%r3d(i,j,k1)
            call mxdq_std(1,npvs,pvk,ct,cx,cy,cz,1,neqn,dqk,dfk,srck,one)
        end if

        rhs0(:) = dfi(:) + dfj(:) + dfk(:)


        if (nvis > nvis_euler) then
            srvi = srv(1)%r3d(i1,j,k)
            srvj = srv(2)%r3d(i,j1,k)
            srvk = srv(3)%r3d(i,j,k1)
            do m=1,neqn
                rhs0(m) = rhs0(m) + &
                          half*(srvi*dqi(m) + srvj*dqj(m) + srvk*dqk(m))
            end do
        end if

        odia = one/rdt(1)%r3d(i,j,k)
        do m=1,neqn
            dq(m)%r3d(i,j,k) = (rhs(m)%r3d(i,j,k) + rhs0(m))*odia
        end do
    end do
    end do
!$OMP end do nowait

#ifdef OMP_IMP
    if (idt < min(numt-1,ed(2)-st(2)+1)) then
        !$OMP flush(nflag)
        do while (nflag(idt) == 1)
            !$OMP flush(nflag)
        end do
        nflag(idt) = 1
        !$OMP flush(nflag)
    end if
#endif

    end do  ! K loop

!$OMP end parallel
end subroutine lusgs_std_sca_foreward_sp

subroutine lusgs_std_sca_backward(nb)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,one,half
    use mod_constants, only : nvis_euler,nsw_dir_close
    use mod_datatypes, only : fld_array_t
    use mod_variables, only : nvis,nsw_kdir
    use mod_fieldvars, only : mb_top,mb_sxyz,mb_pv
    use mod_fieldvars, only : mb_dq,npvs,neqn
    use mod_fieldvars, only : mb_src,mb_srv,mb_rdt
    use mod_openmp
    implicit none
    integer(kind_int), intent(in) :: nb
    integer(kind_int)          :: i,j,k,i1,j1,k1,m,ierr
    integer(kind_int)          :: st(3),ed(3),dijk
    real(kind_real)            :: srci,srcj,srck,srvi,srvj,srvk
    real(kind_real)            :: kt,kx,ky,kz,et,ex,ey,ez,ct,cx,cy,cz
    real(kind_real)            :: odia,rhs0(1:neqn)
    real(kind_real)            :: pvi(1:npvs),pvj(1:npvs),pvk(1:npvs)
    real(kind_real)            :: dqi(1:neqn),dqj(1:neqn),dqk(1:neqn)
    real(kind_real)            :: dfi(1:neqn),dfj(1:neqn),dfk(1:neqn)
    type(fld_array_t), pointer :: sxyz(:),pv(:)
    type(fld_array_t), pointer :: dq(:),src(:),srv(:),rdt(:)
#ifdef OMP_IMP
    integer(kind_int)          :: idt,numt,nflag(0:omp_max_num_threads)
#endif

    sxyz => mb_sxyz(nb)%fld
    pv   => mb_pv(nb)%fld
    dq   => mb_dq(nb)%fld
    src  => mb_src(nb)%fld
    rdt  => mb_rdt(nb)%fld

    if (nvis > nvis_euler) then
        srv  => mb_srv(nb)%fld
    end if


    st(:) = mb_top(nb)%nded(:)
    ed(:) = mb_top(nb)%ndst(:)
    dijk = -1

!$OMP parallel private(k,k1,j,j1,i,i1,idt,pvi,pvj,pvk,dqi,dqj,dqk,rhs0,dfi,dfj,dfk,&
!$OMP                  kx,ky,kz,ex,ey,ez,cx,cy,cz,srvi,srvj,srvk,srci,srcj,srck,odia)

#ifdef OMP_IMP
!$OMP master
    numt = omp_get_num_threads()
    if (numt > omp_max_num_threads) then
        write(*,*) "OpenMP num_threads > omp_max_num_threads, make sure you want do it, check size of nflag"
    end if
!$OMP end master
    idt  = omp_get_thread_num()
    nflag(idt) = 0
#endif

!$OMP barrier

    do k=st(3),ed(3),dijk

#ifdef OMP_IMP
    if (idt > 0 .and. idt < numt) then
        do while (nflag(idt-1) == 0)
            !$OMP flush(nflag)
        end do
        nflag(idt-1) = 0
        !$OMP flush(nflag)
    end if
#endif

!$OMP do schedule(static)
    do j=st(2),ed(2),dijk
    do i=st(1),ed(1),dijk
        i1 = i + 1
        j1 = j + 1
        k1 = k + 1

        srci = src(1)%r3d(i1,j,k)
        srcj = src(2)%r3d(i,j1,k)
        srck = src(3)%r3d(i,j,k1)

        do m=1,npvs
            pvi(m) = pv(m)%r3d(i1,j,k)
            pvj(m) = pv(m)%r3d(i,j1,k)
            pvk(m) = pv(m)%r3d(i,j,k1)
        end do

        do m=1,neqn
            dqi(m) = dq(m)%r3d(i1,j,k)
            dqj(m) = dq(m)%r3d(i,j1,k)
            dqk(m) = dq(m)%r3d(i,j,k1)
        end do

        kx = sxyz(1)%r3d(i1,j,k)
        ky = sxyz(2)%r3d(i1,j,k)
        kz = sxyz(3)%r3d(i1,j,k)
        call mxdq_std(1,npvs,pvi,kt,kx,ky,kz,1,neqn,dqi,dfi,srci,-one)

        ex = sxyz(4)%r3d(i,j1,k)
        ey = sxyz(5)%r3d(i,j1,k)
        ez = sxyz(6)%r3d(i,j1,k)
        call mxdq_std(1,npvs,pvj,et,ex,ey,ez,1,neqn,dqj,dfj,srcj,-one)

        if (nsw_kdir == nsw_dir_close) then
            dfk(:) = zero
        else
            cx = sxyz(7)%r3d(i,j,k1)
            cy = sxyz(8)%r3d(i,j,k1)
            cz = sxyz(9)%r3d(i,j,k1)
            call mxdq_std(1,npvs,pvk,ct,cx,cy,cz,1,neqn,dqk,dfk,srck,-one)
        end if

        rhs0(:) = dfi(:) + dfj(:) + dfk(:)


        if (nvis > nvis_euler) then
            srvi = srv(1)%r3d(i1,j,k)
            srvj = srv(2)%r3d(i,j1,k)
            srvk = srv(3)%r3d(i,j,k1)
            do m=1,neqn
                rhs0(m) = rhs0(m) - &
                          half*(srvi*dqi(m) + srvj*dqj(m) + srvk*dqk(m))
            end do
        end if

        odia = one/rdt(1)%r3d(i,j,k)
        do m=1,neqn
            dq(m)%r3d(i,j,k) = dq(m)%r3d(i,j,k) - rhs0(m)*odia
        end do
    end do
    end do
!$OMP end do nowait

#ifdef OMP_IMP
    if (idt < min(numt-1,st(2)-ed(2)+1)) then
        !$OMP flush(nflag)
        do while (nflag(idt) == 1)
            !$OMP flush(nflag)
        end do
        nflag(idt) = 1
        !$OMP flush(nflag)
    end if
#endif

    end do  ! K loop

!$OMP end parallel
end subroutine lusgs_std_sca_backward

subroutine lusgs_std_sca_backward_sp(nb)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,one,half
    use mod_constants, only : nvis_euler,nsw_dir_close
    use mod_datatypes, only : fld_array_t
    use mod_variables, only : nvis,nsw_kdir
    use mod_fieldvars, only : mb_topsp,mb_sxyzsp,mb_pv
    use mod_fieldvars, only : mb_dq,npvs,neqn
    use mod_fieldvars, only : mb_src,mb_srv,mb_rdt
    use mod_openmp
    implicit none
    integer(kind_int), intent(in) :: nb
    integer(kind_int)          :: i,j,k,i1,j1,k1,m,ierr
    integer(kind_int)          :: st(3),ed(3),dijk
    real(kind_real)            :: srci,srcj,srck,srvi,srvj,srvk
    real(kind_real)            :: kt,kx,ky,kz,et,ex,ey,ez,ct,cx,cy,cz
    real(kind_real)            :: odia,rhs0(1:neqn)
    real(kind_real)            :: pvi(1:npvs),pvj(1:npvs),pvk(1:npvs)
    real(kind_real)            :: dqi(1:neqn),dqj(1:neqn),dqk(1:neqn)
    real(kind_real)            :: dfi(1:neqn),dfj(1:neqn),dfk(1:neqn)
    type(fld_array_t), pointer :: sxyz(:),pv(:)
    type(fld_array_t), pointer :: dq(:),src(:),srv(:),rdt(:)
#ifdef OMP_IMP
    integer(kind_int)          :: idt,numt,nflag(0:omp_max_num_threads)
#endif

    sxyz => mb_sxyzsp(nb)%fld
    pv   => mb_pv(nb)%fld
    dq   => mb_dq(nb)%fld
    src  => mb_src(nb)%fld
    rdt  => mb_rdt(nb)%fld

    if (nvis > nvis_euler) then
        srv  => mb_srv(nb)%fld
    end if


    st(:) = mb_topsp(nb)%nded(:)
    ed(:) = mb_topsp(nb)%ndst(:)
    dijk = -1

!$OMP parallel private(k,k1,j,j1,i,i1,idt,pvi,pvj,pvk,dqi,dqj,dqk,rhs0,dfi,dfj,dfk,&
!$OMP                  kx,ky,kz,ex,ey,ez,cx,cy,cz,srvi,srvj,srvk,srci,srcj,srck,odia)

#ifdef OMP_IMP
!$OMP master
    numt = omp_get_num_threads()
    if (numt > omp_max_num_threads) then
        write(*,*) "OpenMP num_threads > omp_max_num_threads, make sure you want do it, check size of nflag"
    end if
!$OMP end master
    idt  = omp_get_thread_num()
    nflag(idt) = 0
#endif

!$OMP barrier

    do k=st(3),ed(3),dijk

#ifdef OMP_IMP
    if (idt > 0 .and. idt < numt) then
        do while (nflag(idt-1) == 0)
            !$OMP flush(nflag)
        end do
        nflag(idt-1) = 0
        !$OMP flush(nflag)
    end if
#endif

!$OMP do schedule(static)
    do j=st(2),ed(2),dijk
    do i=st(1),ed(1),dijk
        i1 = i + 1
        j1 = j + 1
        k1 = k + 1

        srci = src(1)%r3d(i1,j,k)
        srcj = src(2)%r3d(i,j1,k)
        srck = src(3)%r3d(i,j,k1)

        do m=1,npvs
            pvi(m) = pv(m)%r3d(i1,j,k)
            pvj(m) = pv(m)%r3d(i,j1,k)
            pvk(m) = pv(m)%r3d(i,j,k1)
        end do

        do m=1,neqn
            dqi(m) = dq(m)%r3d(i1,j,k)
            dqj(m) = dq(m)%r3d(i,j1,k)
            dqk(m) = dq(m)%r3d(i,j,k1)
        end do

        kx = sxyz(1)%r3d(i1,j,k)
        ky = sxyz(2)%r3d(i1,j,k)
        kz = sxyz(3)%r3d(i1,j,k)
        call mxdq_std(1,npvs,pvi,kt,kx,ky,kz,1,neqn,dqi,dfi,srci,-one)

        ex = sxyz(4)%r3d(i,j1,k)
        ey = sxyz(5)%r3d(i,j1,k)
        ez = sxyz(6)%r3d(i,j1,k)
        call mxdq_std(1,npvs,pvj,et,ex,ey,ez,1,neqn,dqj,dfj,srcj,-one)

        if (nsw_kdir == nsw_dir_close) then
            dfk(:) = zero
        else
            cx = sxyz(7)%r3d(i,j,k1)
            cy = sxyz(8)%r3d(i,j,k1)
            cz = sxyz(9)%r3d(i,j,k1)
            call mxdq_std(1,npvs,pvk,ct,cx,cy,cz,1,neqn,dqk,dfk,srck,-one)
        end if

        rhs0(:) = dfi(:) + dfj(:) + dfk(:)


        if (nvis > nvis_euler) then
            srvi = srv(1)%r3d(i1,j,k)
            srvj = srv(2)%r3d(i,j1,k)
            srvk = srv(3)%r3d(i,j,k1)
            do m=1,neqn
                rhs0(m) = rhs0(m) - &
                          half*(srvi*dqi(m) + srvj*dqj(m) + srvk*dqk(m))
            end do
        end if

        odia = one/rdt(1)%r3d(i,j,k)
        do m=1,neqn
            dq(m)%r3d(i,j,k) = dq(m)%r3d(i,j,k) - rhs0(m)*odia
        end do
    end do
    end do
!$OMP end do nowait

#ifdef OMP_IMP
    if (idt < min(numt-1,st(2)-ed(2)+1)) then
        !$OMP flush(nflag)
        do while (nflag(idt) == 1)
            !$OMP flush(nflag)
        end do
        nflag(idt) = 1
        !$OMP flush(nflag)
    end if
#endif

    end do  ! K loop

!$OMP end parallel
end subroutine lusgs_std_sca_backward_sp


