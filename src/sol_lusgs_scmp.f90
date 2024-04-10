
subroutine sol_lusgs_std_scmp
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

        call lusgs_std_scmp_foreward(nb)
        call assign_bc_var_nb(nb,mb_dq,1,neqn,0,dq0)

        call lusgs_std_scmp_backward(nb)
        call assign_bc_var_nb(nb,mb_dq,1,neqn,0,dq0)

        call assign_var_via_com_nb(nb,mb_dq,1,neqn)
    end do
    
    call residual(ntimeadv_steady)

    call update_scmp

end subroutine sol_lusgs_std_scmp

subroutine sol_lusgs_std_scmp_sp
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

        call lusgs_std_scmp_foreward_sp(nb)
        !call assign_bc_var_nb(nb,mb_dq,1,neqn,0,dq0)

        call lusgs_std_scmp_backward_sp(nb)
        !call assign_bc_var_nb(nb,mb_dq,1,neqn,0,dq0)

        call assign_var_via_com_nb_sp(nb,mb_dq,1,neqn)
    end do
    
    call residual_sp(ntimeadv_steady)

    call update_scmp_sp

end subroutine sol_lusgs_std_scmp_sp

subroutine lusgs_std_scmp_foreward(nb)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,one,half,nprec_non
    use mod_constants, only : nvis_euler,nsw_dir_close
    use mod_datatypes, only : fld_array_t
    use mod_variables, only : nvis,nsw_kdir,nprec,fsw_kdir,poo
    use mod_fieldvars, only : mb_top,mb_sxyz,mb_pv
    use mod_fieldvars, only : mb_dq,npvs,neqn,mb_rhs
    use mod_fieldvars, only : mb_src,mb_srv,mb_rdt,mb_vsl,mb_vst,mb_vol,mb_dt
    !use mod_fieldvars, only : mb_nonprec,mb_prec_diag_inv
    use mod_interface, only : assign_bc_var_nb
    use mod_openmp
    implicit none
    integer(kind_int), intent(in) :: nb
    integer(kind_int)          :: i,j,k,l,i1,j1,k1,m,ierr,oned_index,nsp,nep
    integer(kind_int)          :: st(3),ed(3),dijk
    real(kind_real)            :: srci,srcj,srck,srvi,srvj,srvk
    real(kind_real)            :: kt,kx,ky,kz,et,ex,ey,ez,ct,cx,cy,cz
    real(kind_real)            :: odia,rhs0(1:neqn)
    real(kind_real)            :: pvi(1:npvs),pvj(1:npvs),pvk(1:npvs)
    real(kind_real)            :: dqi(1:neqn),dqj(1:neqn),dqk(1:neqn)
    real(kind_real)            :: dfi(1:neqn),dfj(1:neqn),dfk(1:neqn)
    type(fld_array_t), pointer :: sxyz(:),pv(:),rhs(:)
    type(fld_array_t), pointer :: dq(:),src(:),srv(:),rdt(:),vsl(:),vst(:),vol(:),dt(:)
    !type(fld_array_t), pointer :: nonprec(:),diagInv(:)
    real(kind_real)            :: point_nonprec_i(1:4),point_nonprec_j(1:4),point_nonprec_k(1:4)
    real(kind_real)            :: dUi(1:5),dUj(1:5),dUk(1:5),dfvi(1:5),dfvj(1:5),dfvk(1:5)
    real(kind_real)            :: gamaInvDfi(1:5),gamaInvDfj(1:5),gamaInvDfk(1:5),gamaInvDf(1:5)
    real(kind_real)            :: gamaInvDfvi(1:5),gamaInvDfvj(1:5),gamaInvDfvk(1:5),gamaInvDfv(1:5)
    real(kind_real)            :: gamaInvRhs(1:5)
    real(kind_real)            :: diagInvMatrix(1:25),prim(1:5)
    real(kind_real)            :: visl,vist,vis,length,rca,rcb,rcc,dia_first,odia_first
#ifdef OMP_IMP
    integer(kind_int)          :: idt,numt,nflag(0:omp_max_num_threads)
#endif

    sxyz => mb_sxyz(nb)%fld
    pv   => mb_pv(nb)%fld
    rhs  => mb_rhs(nb)%fld
    dq   => mb_dq(nb)%fld
    src  => mb_src(nb)%fld
    rdt  => mb_rdt(nb)%fld
    dt   => mb_dt(nb)%fld

    if (nvis > nvis_euler) then
        srv  => mb_srv(nb)%fld
        vol  => mb_vol(nb)%fld
        vsl  => mb_vsl(nb)%fld
        vst  => mb_vst(nb)%fld
    end if
    
    nsp = 1
    nep = 5
    vis = 0.
    length = 1.0

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
        
        !> SCM-P 将原始变量转化成新的SCM-P的变量 (ρ,u,v,w,p) ---> (ρ,u,v,w,p')
        pvi(npvs) = pvi(npvs) - poo     !> Done
        pvj(npvs) = pvj(npvs) - poo     !> Done
        pvk(npvs) = pvk(npvs) - poo     !> Done
        
        !> SCM-P守恒变量增量ΔT = Δ(p',ρu,ρv,ρw)
        do m=1,neqn-1
            dqi(m) = dq(m)%r3d(i1,j,k)  !> Done
            dqj(m) = dq(m)%r3d(i,j1,k)  !> Done
            dqk(m) = dq(m)%r3d(i,j,k1)  !> Done
        end do

        !todo 采用新的特征系统分解计算L*ΔT(i-1,j,k)
        kx = sxyz(1)%r3d(i1,j,k)!网格导数
        ky = sxyz(2)%r3d(i1,j,k)
        kz = sxyz(3)%r3d(i1,j,k)
        call mxdq_scmp(1,npvs,pvi,kt,kx,ky,kz,1,neqn,dqi,dfi,srci,one)

        !todo 采用新的特征系统分解计算L*ΔT(i,j-1,k)
        ex = sxyz(4)%r3d(i,j1,k)
        ey = sxyz(5)%r3d(i,j1,k)
        ez = sxyz(6)%r3d(i,j1,k)
        call mxdq_scmp(1,npvs,pvj,et,ex,ey,ez,1,neqn,dqj,dfj,srcj,one)

        !todo 采用新的特征系统分解计算L*ΔT(i,j,k-1)
        if (nsw_kdir == nsw_dir_close) then
            dfk(:) = zero
        else
            cx = sxyz(7)%r3d(i,j,k1)
            cy = sxyz(8)%r3d(i,j,k1)
            cz = sxyz(9)%r3d(i,j,k1)
            call mxdq_scmp(1,npvs,pvk,ct,cx,cy,cz,1,neqn,dqk,dfk,srck,one)
        end if
        
        rhs0(:) = dfi(:) + dfj(:) + dfk(:)
        
        if (nvis > nvis_euler) then
            srvi = srv(1)%r3d(i1,j,k)
            srvj = srv(2)%r3d(i,j1,k)
            srvk = srv(3)%r3d(i,j,k1)
            !计算黏性项
            dfvi(:) = 0.5*srvi*dqi(:)
            dfvj(:) = 0.5*srvj*dqj(:)
            dfvk(:) = 0.5*srvk*dqk(:)
            
            !< SCMP 只有四个方程p',ρu,ρv,ρw，其中第一个方程p'对应的粘性通量为0
            do m=2,4
                rhs0(m) = rhs0(m) + dfvi(m) + dfvj(m) + dfvk(m)     !> Done
            end do
        end if        
        
        !< SCMP 只有四个方程
        do m=1,4
            rhs0(m) = rhs0(m) + rhs(m)%r3d(i,j,k)
        end do
        
        !< SCMP D矩阵第1行对角元素需去掉粘性谱半径
        odia = one/rdt(1)%r3d(i,j,k)
        if (nvis > nvis_euler) then
            dia_first = rdt(1)%r3d(i,j,k) - (srv(1)%r3d(i,j,k) + srv(2)%r3d(i,j,k) + srv(3)%r3d(i,j,k)) !> Done
        else
            dia_first = rdt(1)%r3d(i,j,k)
        end if
        odia_first = 1.0/dia_first                                                                  !> Done
        do m=2,4
            dq(m)%r3d(i,j,k) = rhs0(m)*odia
        end do
        dq(1)%r3d(i,j,k) = rhs0(1)*odia_first     !> Done
        
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
end subroutine lusgs_std_scmp_foreward

subroutine lusgs_std_scmp_foreward_sp(nb)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,one,half,nprec_non
    use mod_constants, only : nvis_euler,nsw_dir_close
    use mod_datatypes, only : fld_array_t
    use mod_variables, only : nvis,nsw_kdir,nprec,fsw_kdir,poo
    use mod_fieldvars, only : mb_topsp,mb_sxyzsp,mb_pv
    use mod_fieldvars, only : mb_dq,npvs,neqn,mb_rhs
    use mod_fieldvars, only : mb_src,mb_srv,mb_rdt,mb_vsl,mb_vst,mb_volsp,mb_dt
    !use mod_fieldvars, only : mb_nonprec,mb_prec_diag_inv
    use mod_openmp
    implicit none
    integer(kind_int), intent(in) :: nb
    integer(kind_int)          :: i,j,k,l,i1,j1,k1,m,ierr,oned_index,nsp,nep
    integer(kind_int)          :: st(3),ed(3),dijk
    real(kind_real)            :: srci,srcj,srck,srvi,srvj,srvk
    real(kind_real)            :: kt,kx,ky,kz,et,ex,ey,ez,ct,cx,cy,cz
    real(kind_real)            :: odia,rhs0(1:neqn)
    real(kind_real)            :: pvi(1:npvs),pvj(1:npvs),pvk(1:npvs)
    real(kind_real)            :: dqi(1:neqn),dqj(1:neqn),dqk(1:neqn)
    real(kind_real)            :: dfi(1:neqn),dfj(1:neqn),dfk(1:neqn)
    type(fld_array_t), pointer :: sxyz(:),pv(:),rhs(:)
    type(fld_array_t), pointer :: dq(:),src(:),srv(:),rdt(:),vsl(:),vst(:),vol(:),dt(:)
    !type(fld_array_t), pointer :: nonprec(:),diagInv(:)
    real(kind_real)            :: point_nonprec_i(1:4),point_nonprec_j(1:4),point_nonprec_k(1:4)
    real(kind_real)            :: dUi(1:5),dUj(1:5),dUk(1:5),dfvi(1:5),dfvj(1:5),dfvk(1:5)
    real(kind_real)            :: gamaInvDfi(1:5),gamaInvDfj(1:5),gamaInvDfk(1:5),gamaInvDf(1:5)
    real(kind_real)            :: gamaInvDfvi(1:5),gamaInvDfvj(1:5),gamaInvDfvk(1:5),gamaInvDfv(1:5)
    real(kind_real)            :: gamaInvRhs(1:5)
    real(kind_real)            :: diagInvMatrix(1:25),prim(1:5)
    real(kind_real)            :: visl,vist,vis,length,rca,rcb,rcc,dia_first,odia_first
#ifdef OMP_IMP
    integer(kind_int)          :: idt,numt,nflag(0:omp_max_num_threads)
#endif

    sxyz => mb_sxyzsp(nb)%fld
    pv   => mb_pv(nb)%fld
    rhs  => mb_rhs(nb)%fld
    dq   => mb_dq(nb)%fld
    src  => mb_src(nb)%fld
    rdt  => mb_rdt(nb)%fld
    dt   => mb_dt(nb)%fld

    if (nvis > nvis_euler) then
        srv  => mb_srv(nb)%fld
        vol  => mb_volsp(nb)%fld
        vsl  => mb_vsl(nb)%fld
        vst  => mb_vst(nb)%fld
    end if
    
    nsp = 1
    nep = 5
    vis = 0.
    length = 1.0

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
        
        !> SCM-P 将原始变量转化成新的SCM-P的变量 (ρ,u,v,w,p) ---> (ρ,u,v,w,p')
        pvi(npvs) = pvi(npvs) - poo     !> Done
        pvj(npvs) = pvj(npvs) - poo     !> Done
        pvk(npvs) = pvk(npvs) - poo     !> Done
        
        !> SCM-P守恒变量增量ΔT = Δ(p',ρu,ρv,ρw)
        do m=1,neqn-1
            dqi(m) = dq(m)%r3d(i1,j,k)  !> Done
            dqj(m) = dq(m)%r3d(i,j1,k)  !> Done
            dqk(m) = dq(m)%r3d(i,j,k1)  !> Done
        end do

        !todo 采用新的特征系统分解计算L*ΔT(i-1,j,k)
        kx = sxyz(1)%r3d(i1,j,k)!网格导数
        ky = sxyz(2)%r3d(i1,j,k)
        kz = sxyz(3)%r3d(i1,j,k)
        call mxdq_scmp(1,npvs,pvi,kt,kx,ky,kz,1,neqn,dqi,dfi,srci,one)

        !todo 采用新的特征系统分解计算L*ΔT(i,j-1,k)
        ex = sxyz(4)%r3d(i,j1,k)
        ey = sxyz(5)%r3d(i,j1,k)
        ez = sxyz(6)%r3d(i,j1,k)
        call mxdq_scmp(1,npvs,pvj,et,ex,ey,ez,1,neqn,dqj,dfj,srcj,one)

        !todo 采用新的特征系统分解计算L*ΔT(i,j,k-1)
        if (nsw_kdir == nsw_dir_close) then
            dfk(:) = zero
        else
            cx = sxyz(7)%r3d(i,j,k1)
            cy = sxyz(8)%r3d(i,j,k1)
            cz = sxyz(9)%r3d(i,j,k1)
            call mxdq_scmp(1,npvs,pvk,ct,cx,cy,cz,1,neqn,dqk,dfk,srck,one)
        end if
        
        rhs0(:) = dfi(:) + dfj(:) + dfk(:)
        
        if (nvis > nvis_euler) then
            srvi = srv(1)%r3d(i1,j,k)
            srvj = srv(2)%r3d(i,j1,k)
            srvk = srv(3)%r3d(i,j,k1)
            !计算黏性项
            dfvi(:) = 0.5*srvi*dqi(:)
            dfvj(:) = 0.5*srvj*dqj(:)
            dfvk(:) = 0.5*srvk*dqk(:)
            
            !< SCMP 只有四个方程p',ρu,ρv,ρw，其中第一个方程p'对应的粘性通量为0
            do m=2,4
                rhs0(m) = rhs0(m) + dfvi(m) + dfvj(m) + dfvk(m)     !> Done
            end do
        end if        
        
        !< SCMP 只有四个方程
        do m=1,4
            rhs0(m) = rhs0(m) + rhs(m)%r3d(i,j,k)
        end do
        
        !< SCMP D矩阵第1行对角元素需去掉粘性谱半径
        odia = one/rdt(1)%r3d(i,j,k)
        if (nvis > nvis_euler) then
            dia_first = rdt(1)%r3d(i,j,k) - (srv(1)%r3d(i,j,k) + srv(2)%r3d(i,j,k) + srv(3)%r3d(i,j,k)) !> Done
        else
            dia_first = rdt(1)%r3d(i,j,k)
        end if
        odia_first = 1.0/dia_first                                                                  !> Done
        do m=2,4
            dq(m)%r3d(i,j,k) = rhs0(m)*odia
        end do
        dq(1)%r3d(i,j,k) = rhs0(1)*odia_first     !> Done
        
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
end subroutine lusgs_std_scmp_foreward_sp

subroutine lusgs_std_scmp_backward(nb)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,one,half
    use mod_constants, only : nvis_euler,nsw_dir_close,nprec_non
    use mod_datatypes, only : fld_array_t
    use mod_variables, only : nvis,nsw_kdir,nprec,fsw_kdir
    use mod_variables, only : refbeta,gamma,poo
    use mod_fieldvars, only : mb_top,mb_sxyz,mb_pv
    use mod_fieldvars, only : mb_dq,npvs,neqn
    use mod_fieldvars, only : mb_src,mb_srv,mb_rdt,mb_vsl,mb_vst,mb_dt,mb_vol
    use mod_interface, only : assign_bc_var_nb
    use mod_openmp
    implicit none
    integer(kind_int), intent(in) :: nb
    integer(kind_int)          :: i,j,k,l,i1,j1,k1,m,ierr,oned_index,nsp,nep
    integer(kind_int)          :: st(3),ed(3),dijk
    real(kind_real)            :: srci,srcj,srck,srvi,srvj,srvk
    real(kind_real)            :: kt,kx,ky,kz,et,ex,ey,ez,ct,cx,cy,cz
    real(kind_real)            :: odia,rhs0(1:neqn)
    real(kind_real)            :: pvi(1:npvs),pvj(1:npvs),pvk(1:npvs),pvijk(1:npvs)
    real(kind_real)            :: dqi(1:neqn),dqj(1:neqn),dqk(1:neqn)
    real(kind_real)            :: dfi(1:neqn),dfj(1:neqn),dfk(1:neqn)
    type(fld_array_t), pointer :: sxyz(:),pv(:)
    type(fld_array_t), pointer :: dq(:),src(:),srv(:),rdt(:),vsl(:),vst(:),dt(:),vol(:)
    real(kind_real)            :: point_nonprec_i(1:4),point_nonprec_j(1:4),point_nonprec_k(1:4)
    real(kind_real)            :: dUi(1:5),dUj(1:5),dUk(1:5),dfvi(1:5),dfvj(1:5),dfvk(1:5)
    real(kind_real)            :: gamaInvDfi(1:5),gamaInvDfj(1:5),gamaInvDfk(1:5),gamaInvDf(1:5)
    real(kind_real)            :: gamaInvDfvi(1:5),gamaInvDfvj(1:5),gamaInvDfvk(1:5),gamaInvDfv(1:5)
    real(kind_real)            :: gamaInvRhs(1:5)
    real(kind_real)            :: diagInvMatrix(1:25),prim(1:5)
    real(kind_real)            :: rca,rcb,rcc,length,visl,vist,vis,dia_first,odia_first
    real(kind_real)            :: ro,u,v,w,T,pressure,uu,dro,du,dv,dw,dp,dTemp
#ifdef OMP_IMP
    integer(kind_int)          :: idt,numt,nflag(0:omp_max_num_threads)
#endif

    sxyz => mb_sxyz(nb)%fld
    pv   => mb_pv(nb)%fld
    dq   => mb_dq(nb)%fld
    src  => mb_src(nb)%fld
    rdt  => mb_rdt(nb)%fld
    dt   => mb_dt(nb)%fld

    if (nvis > nvis_euler) then
        srv  => mb_srv(nb)%fld
        vol  => mb_vol(nb)%fld
        vsl  => mb_vsl(nb)%fld
        vst  => mb_vst(nb)%fld
    end if
    
    nsp = 1
    nep = 5
    vis = 0.
    length = 1.

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
        
        !> SCM-P 将原始变量转化成新的SCM-P的变量 (ρ,u,v,w,p) ---> (ρ,u,v,w,p')
        pvi(npvs) = pvi(npvs) - poo     !> Done
        pvj(npvs) = pvj(npvs) - poo     !> Done
        pvk(npvs) = pvk(npvs) - poo     !> Done
        
        !> SCM-P守恒变量增量ΔT = Δ(p',ρu,ρv,ρw)
        do m=1,neqn-1
            dqi(m) = dq(m)%r3d(i1,j,k)  !> Done
            dqj(m) = dq(m)%r3d(i,j1,k)  !> Done
            dqk(m) = dq(m)%r3d(i,j,k1)  !> Done
        end do

        kx = sxyz(1)%r3d(i1,j,k)
        ky = sxyz(2)%r3d(i1,j,k)
        kz = sxyz(3)%r3d(i1,j,k)
        !todo 采用新的特征系统分解计算U*ΔT(i+1,j,k)
        call mxdq_scmp(1,npvs,pvi,kt,kx,ky,kz,1,neqn,dqi,dfi,srci,-one)

        ex = sxyz(4)%r3d(i,j1,k)
        ey = sxyz(5)%r3d(i,j1,k)
        ez = sxyz(6)%r3d(i,j1,k)
        !todo 采用新的特征系统分解计算U*ΔT(i,j+1,k)
        call mxdq_scmp(1,npvs,pvj,et,ex,ey,ez,1,neqn,dqj,dfj,srcj,-one)

        if (nsw_kdir == nsw_dir_close) then
            dfk(:) = zero
        else
            cx = sxyz(7)%r3d(i,j,k1)
            cy = sxyz(8)%r3d(i,j,k1)
            cz = sxyz(9)%r3d(i,j,k1)
            !todo 采用新的特征系统分解计算U*ΔT(i,j,k+1)
            call mxdq_scmp(1,npvs,pvk,ct,cx,cy,cz,1,neqn,dqk,dfk,srck,-one)
        end if

        rhs0(:) = dfi(:) + dfj(:) + dfk(:)

        if (nvis > nvis_euler) then
            srvi = srv(1)%r3d(i1,j,k)
            srvj = srv(2)%r3d(i,j1,k)
            srvk = srv(3)%r3d(i,j,k1)
            
            !计算黏性项
            dfvi(:) = 0.5*srvi*dqi(:)
            dfvj(:) = 0.5*srvj*dqj(:)
            dfvk(:) = 0.5*srvk*dqk(:)
            
            !< SCMP 只有四个方程p',ρu,ρv,ρw，其中第一个方程p'对应的粘性通量为0
            do m=2,4
                rhs0(m) = rhs0(m) - dfvi(m) -dfvj(m) -dfvk(m)     !> Done
            end do
        end if   
        
        !< SCMP 只有四个方程
        !do m=1,4
        !    rhs0(m) = rhs0(m) + rhs(m)%r3d(i,j,k)
        !end do
        
        !< SCMP D矩阵第1行对角元素需去掉粘性谱半径
        odia = one/rdt(1)%r3d(i,j,k)
        if (nvis > nvis_euler) then
            dia_first = rdt(1)%r3d(i,j,k) - (srv(1)%r3d(i,j,k) + srv(2)%r3d(i,j,k) + srv(3)%r3d(i,j,k)) !> Done
        else
            dia_first = rdt(1)%r3d(i,j,k)
        end if
        odia_first = 1.0/dia_first                                                                  !> Done
        do l = 2,4
            dq(l)%r3d(i,j,k) = dq(l)%r3d(i,j,k) - rhs0(l)*odia
        end do
        dq(1)%r3d(i,j,k) = dq(1)%r3d(i,j,k) - rhs0(1)*odia_first     !> Done

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
    
    !do k=st(3),ed(3),dijk
    !do j=st(2),ed(2),dijk
    !do i=st(1),ed(1),dijk
    !    do m = 1,npvs
    !        prim(m) = pv(m)%r3d(i,j,k)
    !    end do
    !    dp = dq(1)%r3d(i,j,k)
    !    du = dq(2)%r3d(i,j,k)
    !    dv = dq(3)%r3d(i,j,k)
    !    dw = dq(4)%r3d(i,j,k)
    !    dTemp = dq(5)%r3d(i,j,k)
    !    
    !    !增量dp,du,dv,dw,dT增量
    !    !转换成增量dro,du,dv,dw,dp
    !    pressure = prim(5)
    !    ro = prim(1)
    !    u = prim(2)
    !    v = prim(3)
    !    w = prim(4)
    !    T = pressure/(refbeta*ro)
    !    dro = ro/pressure*dp - ro/T*dTemp
    !    uu = u*u + v*v + w*w
    !    
    !    dq(1)%r3d(i,j,k) = dro
    !    dq(2)%r3d(i,j,k) = u*dro + ro*du
    !    dq(3)%r3d(i,j,k) = v*dro + ro*dv
    !    dq(4)%r3d(i,j,k) = w*dro + ro*dw
    !    dq(5)%r3d(i,j,k) = 1/(gamma - 1.)*dp + 0.5*uu*dro + ro*(u*du + v*dv + w*dw)
    !end do
    !end do
    !end do

!$OMP end parallel
end subroutine lusgs_std_scmp_backward

subroutine lusgs_std_scmp_backward_sp(nb)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,one,half
    use mod_constants, only : nvis_euler,nsw_dir_close,nprec_non
    use mod_datatypes, only : fld_array_t
    use mod_variables, only : nvis,nsw_kdir,nprec,fsw_kdir
    use mod_variables, only : refbeta,gamma,poo
    use mod_fieldvars, only : mb_topsp,mb_sxyzsp,mb_pv
    use mod_fieldvars, only : mb_dq,npvs,neqn
    use mod_fieldvars, only : mb_src,mb_srv,mb_rdt,mb_vsl,mb_vst,mb_dt,mb_volsp
    use mod_openmp
    implicit none
    integer(kind_int), intent(in) :: nb
    integer(kind_int)          :: i,j,k,l,i1,j1,k1,m,ierr,oned_index,nsp,nep
    integer(kind_int)          :: st(3),ed(3),dijk
    real(kind_real)            :: srci,srcj,srck,srvi,srvj,srvk
    real(kind_real)            :: kt,kx,ky,kz,et,ex,ey,ez,ct,cx,cy,cz
    real(kind_real)            :: odia,rhs0(1:neqn)
    real(kind_real)            :: pvi(1:npvs),pvj(1:npvs),pvk(1:npvs),pvijk(1:npvs)
    real(kind_real)            :: dqi(1:neqn),dqj(1:neqn),dqk(1:neqn)
    real(kind_real)            :: dfi(1:neqn),dfj(1:neqn),dfk(1:neqn)
    type(fld_array_t), pointer :: sxyz(:),pv(:)
    type(fld_array_t), pointer :: dq(:),src(:),srv(:),rdt(:),vsl(:),vst(:),dt(:),vol(:)
    real(kind_real)            :: point_nonprec_i(1:4),point_nonprec_j(1:4),point_nonprec_k(1:4)
    real(kind_real)            :: dUi(1:5),dUj(1:5),dUk(1:5),dfvi(1:5),dfvj(1:5),dfvk(1:5)
    real(kind_real)            :: gamaInvDfi(1:5),gamaInvDfj(1:5),gamaInvDfk(1:5),gamaInvDf(1:5)
    real(kind_real)            :: gamaInvDfvi(1:5),gamaInvDfvj(1:5),gamaInvDfvk(1:5),gamaInvDfv(1:5)
    real(kind_real)            :: gamaInvRhs(1:5)
    real(kind_real)            :: diagInvMatrix(1:25),prim(1:5)
    real(kind_real)            :: rca,rcb,rcc,length,visl,vist,vis,dia_first,odia_first
    real(kind_real)            :: ro,u,v,w,T,pressure,uu,dro,du,dv,dw,dp,dTemp
#ifdef OMP_IMP
    integer(kind_int)          :: idt,numt,nflag(0:omp_max_num_threads)
#endif

    sxyz => mb_sxyzsp(nb)%fld
    pv   => mb_pv(nb)%fld
    dq   => mb_dq(nb)%fld
    src  => mb_src(nb)%fld
    rdt  => mb_rdt(nb)%fld
    dt   => mb_dt(nb)%fld

    if (nvis > nvis_euler) then
        srv  => mb_srv(nb)%fld
        vol  => mb_volsp(nb)%fld
        vsl  => mb_vsl(nb)%fld
        vst  => mb_vst(nb)%fld
    end if
    
    nsp = 1
    nep = 5
    vis = 0.
    length = 1.

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
        
        !> SCM-P 将原始变量转化成新的SCM-P的变量 (ρ,u,v,w,p) ---> (ρ,u,v,w,p')
        pvi(npvs) = pvi(npvs) - poo     !> Done
        pvj(npvs) = pvj(npvs) - poo     !> Done
        pvk(npvs) = pvk(npvs) - poo     !> Done
        
        !> SCM-P守恒变量增量ΔT = Δ(p',ρu,ρv,ρw)
        do m=1,neqn-1
            dqi(m) = dq(m)%r3d(i1,j,k)  !> Done
            dqj(m) = dq(m)%r3d(i,j1,k)  !> Done
            dqk(m) = dq(m)%r3d(i,j,k1)  !> Done
        end do

        kx = sxyz(1)%r3d(i1,j,k)
        ky = sxyz(2)%r3d(i1,j,k)
        kz = sxyz(3)%r3d(i1,j,k)
        !todo 采用新的特征系统分解计算U*ΔT(i+1,j,k)
        call mxdq_scmp(1,npvs,pvi,kt,kx,ky,kz,1,neqn,dqi,dfi,srci,-one)

        ex = sxyz(4)%r3d(i,j1,k)
        ey = sxyz(5)%r3d(i,j1,k)
        ez = sxyz(6)%r3d(i,j1,k)
        !todo 采用新的特征系统分解计算U*ΔT(i,j+1,k)
        call mxdq_scmp(1,npvs,pvj,et,ex,ey,ez,1,neqn,dqj,dfj,srcj,-one)

        if (nsw_kdir == nsw_dir_close) then
            dfk(:) = zero
        else
            cx = sxyz(7)%r3d(i,j,k1)
            cy = sxyz(8)%r3d(i,j,k1)
            cz = sxyz(9)%r3d(i,j,k1)
            !todo 采用新的特征系统分解计算U*ΔT(i,j,k+1)
            call mxdq_scmp(1,npvs,pvk,ct,cx,cy,cz,1,neqn,dqk,dfk,srck,-one)
        end if

        rhs0(:) = dfi(:) + dfj(:) + dfk(:)

        if (nvis > nvis_euler) then
            srvi = srv(1)%r3d(i1,j,k)
            srvj = srv(2)%r3d(i,j1,k)
            srvk = srv(3)%r3d(i,j,k1)
            
            !计算黏性项
            dfvi(:) = 0.5*srvi*dqi(:)
            dfvj(:) = 0.5*srvj*dqj(:)
            dfvk(:) = 0.5*srvk*dqk(:)
            
            !< SCMP 只有四个方程p',ρu,ρv,ρw，其中第一个方程p'对应的粘性通量为0
            do m=2,4
                rhs0(m) = rhs0(m) - dfvi(m) -dfvj(m) -dfvk(m)     !> Done
            end do
        end if   
        
        !< SCMP 只有四个方程
        !do m=1,4
        !    rhs0(m) = rhs0(m) + rhs(m)%r3d(i,j,k)
        !end do
        
        !< SCMP D矩阵第1行对角元素需去掉粘性谱半径
        odia = one/rdt(1)%r3d(i,j,k)
        if (nvis > nvis_euler) then
            dia_first = rdt(1)%r3d(i,j,k) - (srv(1)%r3d(i,j,k) + srv(2)%r3d(i,j,k) + srv(3)%r3d(i,j,k)) !> Done
        else
            dia_first = rdt(1)%r3d(i,j,k)
        end if
        odia_first = 1.0/dia_first                                                                  !> Done
        do l = 2,4
            dq(l)%r3d(i,j,k) = dq(l)%r3d(i,j,k) - rhs0(l)*odia
        end do
        dq(1)%r3d(i,j,k) = dq(1)%r3d(i,j,k) - rhs0(1)*odia_first     !> Done

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
    
    !do k=st(3),ed(3),dijk
    !do j=st(2),ed(2),dijk
    !do i=st(1),ed(1),dijk
    !    do m = 1,npvs
    !        prim(m) = pv(m)%r3d(i,j,k)
    !    end do
    !    dp = dq(1)%r3d(i,j,k)
    !    du = dq(2)%r3d(i,j,k)
    !    dv = dq(3)%r3d(i,j,k)
    !    dw = dq(4)%r3d(i,j,k)
    !    dTemp = dq(5)%r3d(i,j,k)
    !    
    !    !增量dp,du,dv,dw,dT增量
    !    !转换成增量dro,du,dv,dw,dp
    !    pressure = prim(5)
    !    ro = prim(1)
    !    u = prim(2)
    !    v = prim(3)
    !    w = prim(4)
    !    T = pressure/(refbeta*ro)
    !    dro = ro/pressure*dp - ro/T*dTemp
    !    uu = u*u + v*v + w*w
    !    
    !    dq(1)%r3d(i,j,k) = dro
    !    dq(2)%r3d(i,j,k) = u*dro + ro*du
    !    dq(3)%r3d(i,j,k) = v*dro + ro*dv
    !    dq(4)%r3d(i,j,k) = w*dro + ro*dw
    !    dq(5)%r3d(i,j,k) = 1/(gamma - 1.)*dp + 0.5*uu*dro + ro*(u*du + v*dv + w*dw)
    !end do
    !end do
    !end do

!$OMP end parallel
end subroutine lusgs_std_scmp_backward_sp


