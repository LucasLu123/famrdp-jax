
subroutine sol_prsgs_ust_mat
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,nbc_inter_buf_dqc
    use mod_constants, only : nsgl_aver_art,nsgl_buffer_pvs
    use mod_constants, only : ntimeadv_steady,ntimeadv_unstdy
    use mod_constants, only : nsweep_foreward,nsweep_backward
    use mod_variables, only : nsubmax,nghnode,enfix,niter
    use mod_fieldvars, only : nblkcoms,blkcoms
    use mod_fieldvars, only : neqn,mb_dq,mb_qc,mb_q0,mb_qm
    use mod_interface, only : pre_exchange_bc_var,post_exchange_bc_var,average_bc_var
    use mod_interface, only : exchange_singulars,average_singulars
    use mod_interface, only : assign_mb_var_uniform,assign_bc_var_nb
    use mod_interface, only : calc_mb_var_via_sub
    use mod_interface, only : assign_var_via_com_nb
    implicit none
    integer(kind_int)          :: iter,nc,nb,nstop
    real(kind_real), parameter :: cpr(3)=(/1.5,0.5,1.5/)
    real(kind_real)            :: dq0(1:neqn),cdd
    external                   :: v1_eq_v2

    cdd = enfix

    dq0(:) = zero
    call assign_mb_var_uniform(mb_dq,1,neqn,nghnode,dq0)

    call calc_mb_var_via_sub(mb_qc,1,neqn,v1_eq_v2,mb_qm,1,neqn,nghnode)

    do iter=1,nsubmax

        niter = iter
        
        call prepare_linear_system

        do nc=1,nblkcoms
            nb = blkcoms(nc)%nb

            call prsgs_ust_mat_matrix(nb,cpr,ntimeadv_unstdy,cdd)

            call prsgs_ust_mat_sweep(nb,nsweep_foreward,cpr,ntimeadv_unstdy,cdd)
            call assign_bc_var_nb(nb,mb_dq,1,neqn,0,dq0)

            call prsgs_ust_mat_sweep(nb,nsweep_backward,cpr,ntimeadv_unstdy,cdd)
            call assign_bc_var_nb(nb,mb_dq,1,neqn,0,dq0)

            call assign_var_via_com_nb(nb,mb_dq,1,neqn)
        end do

        call pre_exchange_bc_var(mb_dq,1,neqn,1,nbc_inter_buf_dqc,nsgl_aver_art)
        call exchange_singulars(mb_dq,1,neqn,nsgl_buffer_pvs,nsgl_aver_art)
        call post_exchange_bc_var(mb_dq,1,neqn,1,nbc_inter_buf_dqc,nsgl_aver_art)
        call average_bc_var(mb_dq,1,neqn,nbc_inter_buf_dqc,nsgl_aver_art)
        call average_singulars(mb_dq,1,neqn,nsgl_buffer_pvs,nsgl_aver_art)

        call residual(ntimeadv_unstdy)
        call stop_subiter(iter,nstop)

        call update

        if (nstop /= 0) exit
    end do

    call calc_mb_var_via_sub(mb_qm,1,neqn,v1_eq_v2,mb_q0,1,neqn,nghnode)

end subroutine sol_prsgs_ust_mat

subroutine sol_prsgs_ust_mat_sp
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero
    use mod_constants, only : ntimeadv_steady,ntimeadv_unstdy
    use mod_constants, only : nsweep_foreward,nsweep_backward
    use mod_variables, only : nsubmax,nghnode,enfix,niter
    use mod_fieldvars, only : nblkcoms,blkcomssp
    use mod_fieldvars, only : neqn,mb_dq,mb_qc,mb_q0,mb_qm
    use mod_interface, only : assign_mb_var_uniform_sp
    use mod_interface, only : calc_mb_var_via_sub_sp
    use mod_interface, only : assign_var_via_com_nb_sp
    implicit none
    integer(kind_int)          :: iter,nc,nb,nstop
    real(kind_real), parameter :: cpr(3)=(/1.5,0.5,1.5/)
    real(kind_real)            :: dq0(1:neqn),cdd
    external                   :: v1_eq_v2

    cdd = enfix

    dq0(:) = zero
    call assign_mb_var_uniform_sp(mb_dq,1,neqn,nghnode,dq0)

    call calc_mb_var_via_sub_sp(mb_qc,1,neqn,v1_eq_v2,mb_qm,1,neqn,nghnode)

    do iter=1,nsubmax

        niter = iter
        
        call prepare_linear_system_sp

        do nc=1,nblkcoms
            nb = blkcomssp(nc)%nb

            call prsgs_ust_mat_matrix_sp(nb,cpr,ntimeadv_unstdy,cdd)

            call prsgs_ust_mat_sweep_sp(nb,nsweep_foreward,cpr,ntimeadv_unstdy,cdd)
            !call assign_bc_var_nb(nb,mb_dq,1,neqn,0,dq0)

            call prsgs_ust_mat_sweep_sp(nb,nsweep_backward,cpr,ntimeadv_unstdy,cdd)
            !call assign_bc_var_nb(nb,mb_dq,1,neqn,0,dq0)

            call assign_var_via_com_nb_sp(nb,mb_dq,1,neqn)
        end do

        call residual_sp(ntimeadv_unstdy)
        call stop_subiter_sp(iter,nstop)

        call update_sp

        if (nstop /= 0) exit
    end do

    call calc_mb_var_via_sub_sp(mb_qm,1,neqn,v1_eq_v2,mb_q0,1,neqn,nghnode)

end subroutine sol_prsgs_ust_mat_sp

subroutine sol_prsgs_std_mat
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,nbc_inter_buf_dqc
    use mod_constants, only : nsgl_aver_art,nsgl_buffer_pvs
    use mod_constants, only : ntimeadv_steady,ntimeadv_unstdy
    use mod_constants, only : nsweep_foreward,nsweep_backward
    use mod_variables, only : nsubmax,nghnode,enfix
    use mod_fieldvars, only : nblkcoms,blkcoms,neqn,mb_dq
    use mod_interface, only : pre_exchange_bc_var,post_exchange_bc_var,average_bc_var
    use mod_interface, only : exchange_singulars,average_singulars
    use mod_interface, only : assign_mb_var_uniform,assign_bc_var_nb
    use mod_interface, only : assign_var_via_com_nb
    implicit none
    integer(kind_int)          :: iter,nc,nb,nstop
    real(kind_real), parameter :: cpr(3)=(/0.0,0.0,1.0/)
    real(kind_real)            :: dq0(1:neqn),cdd

    cdd = enfix

    call prepare_linear_system

    dq0(:) = zero
    call assign_mb_var_uniform(mb_dq,1,neqn,nghnode,dq0)

    do nc=1,nblkcoms
        nb = blkcoms(nc)%nb
        call prsgs_ust_mat_matrix(nb,cpr,ntimeadv_steady,cdd)
    end do

    do iter=1,nsubmax
        do nc=1,nblkcoms
            nb = blkcoms(nc)%nb

            call prsgs_ust_mat_sweep(nb,nsweep_foreward,cpr,ntimeadv_steady,cdd)
            call assign_bc_var_nb(nb,mb_dq,1,neqn,0,dq0)

            call prsgs_ust_mat_sweep(nb,nsweep_backward,cpr,ntimeadv_steady,cdd)
            call assign_bc_var_nb(nb,mb_dq,1,neqn,0,dq0)

            call assign_var_via_com_nb(nb,mb_dq,1,neqn)
        end do

        call pre_exchange_bc_var(mb_dq,1,neqn,1,nbc_inter_buf_dqc,nsgl_aver_art)
        call exchange_singulars(mb_dq,1,neqn,nsgl_buffer_pvs,nsgl_aver_art)
        call post_exchange_bc_var(mb_dq,1,neqn,1,nbc_inter_buf_dqc,nsgl_aver_art)
        call average_bc_var(mb_dq,1,neqn,nbc_inter_buf_dqc,nsgl_aver_art)
        call average_singulars(mb_dq,1,neqn,nsgl_buffer_pvs,nsgl_aver_art)

        call residual(ntimeadv_steady)
        call stop_subiter(iter,nstop)

        if (nstop /= 0) exit
    end do

    call update

end subroutine sol_prsgs_std_mat

subroutine sol_prsgs_std_mat_sp
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero
    use mod_constants, only : ntimeadv_steady,ntimeadv_unstdy
    use mod_constants, only : nsweep_foreward,nsweep_backward
    use mod_variables, only : nsubmax,nghnode,enfix
    use mod_fieldvars, only : nblkcoms,blkcomssp,neqn,mb_dq
    use mod_interface, only : assign_mb_var_uniform_sp
    use mod_interface, only : assign_var_via_com_nb_sp
    implicit none
    integer(kind_int)          :: iter,nc,nb,nstop
    real(kind_real), parameter :: cpr(3)=(/0.0,0.0,1.0/)
    real(kind_real)            :: dq0(1:neqn),cdd

    cdd = enfix

    call prepare_linear_system_sp

    dq0(:) = zero
    call assign_mb_var_uniform_sp(mb_dq,1,neqn,nghnode,dq0)

    do nc=1,nblkcoms
        nb = blkcomssp(nc)%nb
        call prsgs_ust_mat_matrix_sp(nb,cpr,ntimeadv_steady,cdd)
    end do

    do iter=1,nsubmax
        do nc=1,nblkcoms
            nb = blkcomssp(nc)%nb

            call prsgs_ust_mat_sweep_sp(nb,nsweep_foreward,cpr,ntimeadv_steady,cdd)
            !call assign_bc_var_nb(nb,mb_dq,1,neqn,0,dq0)

            call prsgs_ust_mat_sweep_sp(nb,nsweep_backward,cpr,ntimeadv_steady,cdd)
            !call assign_bc_var_nb(nb,mb_dq,1,neqn,0,dq0)

            call assign_var_via_com_nb_sp(nb,mb_dq,1,neqn)
        end do

        call residual_sp(ntimeadv_steady)
        call stop_subiter_sp(iter,nstop)

        if (nstop /= 0) exit
    end do

    call update_sp

end subroutine sol_prsgs_std_mat_sp

subroutine prsgs_ust_mat_matrix(nb,cpr,nsw,cdd)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,one,half,thr2nd
    use mod_constants, only : nvis_euler,nsw_dir_close
    use mod_datatypes, only : fld_array_t
    use mod_variables, only : nvis,dtime,nsw_kdir
    use mod_fieldvars, only : mb_top,mb_sxyz,mb_pv
    use mod_fieldvars, only : npvs,neqn,mb_dt
    use mod_fieldvars, only : mb_srv,mb_aiv
    use mod_fieldvars, only : mb_vol,mb_vsl,mb_vst
    implicit none
    integer(kind_int), intent(in) :: nb
    real(kind_real),   intent(in) :: cpr(3)
    integer(kind_int), intent(in) :: nsw
    real(kind_real),   intent(in) :: cdd
    integer(kind_int)          :: i,j,k,m,n,st(3),ed(3),nkst,nked,ierr
    real(kind_real)            :: srvi,srvj,srvk,odt,visl,vist,vol0
    real(kind_real)            :: kt,kx,ky,kz,et,ex,ey,ez,ct,cx,cy,cz
    real(kind_real)            :: pvi(1:npvs)
    real(kind_real)            :: af0(neqn,neqn),afi(neqn,neqn)
    real(kind_real)            :: afj(neqn,neqn),afk(neqn,neqn)
    type(fld_array_t), pointer :: sxyz(:),pv(:),dt(:),vol(:)
    type(fld_array_t), pointer :: srv(:),aiv(:),vsl(:),vst(:)

    st(:) = mb_top(nb)%ndst(:)
    ed(:) = mb_top(nb)%nded(:)

    vol  => mb_vol(nb)%fld
    sxyz => mb_sxyz(nb)%fld
    pv   => mb_pv(nb)%fld
    dt   => mb_dt(nb)%fld
    aiv  => mb_aiv(nb)%fld
    if (nvis > nvis_euler) then
        srv  => mb_srv(nb)%fld
        vsl  => mb_vsl(nb)%fld
        vst  => mb_vst(nb)%fld
    end if

    do k=st(3),ed(3)
    do j=st(2),ed(2)
    do i=st(1),ed(1)
        do m=1,npvs
            pvi(m) = pv(m)%r3d(i,j,k)
        end do

        kx = sxyz(1)%r3d(i,j,k)
        ky = sxyz(2)%r3d(i,j,k)
        kz = sxyz(3)%r3d(i,j,k)
        call aflux_inv(1,npvs,pvi,kt,kx,ky,kz,1,neqn,afi,cdd)

        ex = sxyz(4)%r3d(i,j,k)
        ey = sxyz(5)%r3d(i,j,k)
        ez = sxyz(6)%r3d(i,j,k)
        call aflux_inv(1,npvs,pvi,et,ex,ey,ez,1,neqn,afj,cdd)

        if (nsw_kdir == nsw_dir_close) then
            afk(:,:) = zero
        else
            cx = sxyz(7)%r3d(i,j,k)
            cy = sxyz(8)%r3d(i,j,k)
            cz = sxyz(9)%r3d(i,j,k)
            call aflux_inv(1,npvs,pvi,ct,cx,cy,cz,1,neqn,afk,cdd)
        end if
        af0(:,:) = afi(:,:) + afj(:,:) + afk(:,:)

        vol0 = vol(1)%r3d(i,j,k)
        if (nvis > nvis_euler) then
#define PR_VIS_JACOBIAN
#ifdef PR_VIS_JACOBIAN
            visl = vsl(1)%r3d(i,j,k)
            vist = vst(1)%r3d(i,j,k)
            call aflux_vis(1,npvs,pvi,kt,kx,ky,kz,vol0,visl,vist,1,neqn,afi)
            call aflux_vis(1,npvs,pvi,et,ex,ey,ez,vol0,visl,vist,1,neqn,afj)
            if (nsw_kdir == nsw_dir_close) then
                afk(:,:) = zero
            else
                call aflux_vis(1,npvs,pvi,ct,cx,cy,cz,vol0,visl,vist,1,neqn,afk)
            end if

            do m=1,neqn
            do n=1,neqn
                af0(n,m) = af0(n,m) + (afi(n,m) + afj(n,m) + afk(n,m))
            end do
            end do

#else
!! 最大特征值代替粘性通量Jacbi矩阵。
            srvi = srv(1)%r3d(i,j,k)
            srvj = srv(2)%r3d(i,j,k)
            srvk = srv(3)%r3d(i,j,k)
            do m=1,neqn
                af0(m,m) = af0(m,m) + srvi + srvj + srvk
            end do
#endif
        end if


        odt = one/dt(1)%r3d(i,j,k)
        if (nsw > 0) then
            odt = odt + cpr(3)*vol0/dtime
        end if
        do m=1,neqn
            af0(m,m) = af0(m,m) + odt
        end do

        call brinv(neqn,af0,ierr)

        do m=1,neqn
        do n=1,neqn
            aiv((m-1)*neqn + n)%r3d(i,j,k) = af0(m,n)
        end do
        end do
    end do
    end do
    end do

end subroutine prsgs_ust_mat_matrix

subroutine prsgs_ust_mat_matrix_sp(nb,cpr,nsw,cdd)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,one,half,thr2nd
    use mod_constants, only : nvis_euler,nsw_dir_close
    use mod_datatypes, only : fld_array_t
    use mod_variables, only : nvis,dtime,nsw_kdir
    use mod_fieldvars, only : mb_topsp,mb_sxyzsp,mb_pv
    use mod_fieldvars, only : npvs,neqn,mb_dt
    use mod_fieldvars, only : mb_srv,mb_aiv
    use mod_fieldvars, only : mb_volsp,mb_vsl,mb_vst
    implicit none
    integer(kind_int), intent(in) :: nb
    real(kind_real),   intent(in) :: cpr(3)
    integer(kind_int), intent(in) :: nsw
    real(kind_real),   intent(in) :: cdd
    integer(kind_int)          :: i,j,k,m,n,st(3),ed(3),nkst,nked,ierr
    real(kind_real)            :: srvi,srvj,srvk,odt,visl,vist,vol0
    real(kind_real)            :: kt,kx,ky,kz,et,ex,ey,ez,ct,cx,cy,cz
    real(kind_real)            :: pvi(1:npvs)
    real(kind_real)            :: af0(neqn,neqn),afi(neqn,neqn)
    real(kind_real)            :: afj(neqn,neqn),afk(neqn,neqn)
    type(fld_array_t), pointer :: sxyz(:),pv(:),dt(:),vol(:)
    type(fld_array_t), pointer :: srv(:),aiv(:),vsl(:),vst(:)

    st(:) = mb_topsp(nb)%ndst(:)
    ed(:) = mb_topsp(nb)%nded(:)

    vol  => mb_volsp(nb)%fld
    sxyz => mb_sxyzsp(nb)%fld
    pv   => mb_pv(nb)%fld
    dt   => mb_dt(nb)%fld
    aiv  => mb_aiv(nb)%fld
    if (nvis > nvis_euler) then
        srv  => mb_srv(nb)%fld
        vsl  => mb_vsl(nb)%fld
        vst  => mb_vst(nb)%fld
    end if

    do k=st(3),ed(3)
    do j=st(2),ed(2)
    do i=st(1),ed(1)
        do m=1,npvs
            pvi(m) = pv(m)%r3d(i,j,k)
        end do

        kx = sxyz(1)%r3d(i,j,k)
        ky = sxyz(2)%r3d(i,j,k)
        kz = sxyz(3)%r3d(i,j,k)
        call aflux_inv(1,npvs,pvi,kt,kx,ky,kz,1,neqn,afi,cdd)

        ex = sxyz(4)%r3d(i,j,k)
        ey = sxyz(5)%r3d(i,j,k)
        ez = sxyz(6)%r3d(i,j,k)
        call aflux_inv(1,npvs,pvi,et,ex,ey,ez,1,neqn,afj,cdd)

        if (nsw_kdir == nsw_dir_close) then
            afk(:,:) = zero
        else
            cx = sxyz(7)%r3d(i,j,k)
            cy = sxyz(8)%r3d(i,j,k)
            cz = sxyz(9)%r3d(i,j,k)
            call aflux_inv(1,npvs,pvi,ct,cx,cy,cz,1,neqn,afk,cdd)
        end if
        af0(:,:) = afi(:,:) + afj(:,:) + afk(:,:)

        vol0 = vol(1)%r3d(i,j,k)
        if (nvis > nvis_euler) then
#define PR_VIS_JACOBIAN
#ifdef PR_VIS_JACOBIAN
            visl = vsl(1)%r3d(i,j,k)
            vist = vst(1)%r3d(i,j,k)
            call aflux_vis(1,npvs,pvi,kt,kx,ky,kz,vol0,visl,vist,1,neqn,afi)
            call aflux_vis(1,npvs,pvi,et,ex,ey,ez,vol0,visl,vist,1,neqn,afj)
            if (nsw_kdir == nsw_dir_close) then
                afk(:,:) = zero
            else
                call aflux_vis(1,npvs,pvi,ct,cx,cy,cz,vol0,visl,vist,1,neqn,afk)
            end if

            do m=1,neqn
            do n=1,neqn
                af0(n,m) = af0(n,m) + (afi(n,m) + afj(n,m) + afk(n,m))
            end do
            end do

#else
!! 最大特征值代替粘性通量Jacbi矩阵。
            srvi = srv(1)%r3d(i,j,k)
            srvj = srv(2)%r3d(i,j,k)
            srvk = srv(3)%r3d(i,j,k)
            do m=1,neqn
                af0(m,m) = af0(m,m) + srvi + srvj + srvk
            end do
#endif
        end if


        odt = one/dt(1)%r3d(i,j,k)
        if (nsw > 0) then
            odt = odt + cpr(3)*vol0/dtime
        end if
        do m=1,neqn
            af0(m,m) = af0(m,m) + odt
        end do

        call brinv(neqn,af0,ierr)

        do m=1,neqn
        do n=1,neqn
            aiv((m-1)*neqn + n)%r3d(i,j,k) = af0(m,n)
        end do
        end do
    end do
    end do
    end do

end subroutine prsgs_ust_mat_matrix_sp

subroutine prsgs_ust_mat_sweep(nb,swp,cpr,nsw,cdd)
    ! swp: 控制扫描方向（swp=0向前，否则向后）
    ! cpr: 时间离散系数
    ! nsw: 控制时间离散（nsw=0定常，否则非定常）
    ! cdd: Steger分裂修正系数
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,one,half,thr2nd
    use mod_constants, only : nvis_euler,nsw_dir_close
    use mod_datatypes, only : fld_array_t
    use mod_variables, only : nvis,dtime,nsw_kdir
    use mod_fieldvars, only : mb_top,mb_sxyz,mb_pv
    use mod_fieldvars, only : mb_dq,npvs,neqn,mb_dt
    use mod_fieldvars, only : mb_qc,mb_q0,mb_rhs,mb_qm
    use mod_fieldvars, only : mb_src,mb_srv,mb_aiv
    use mod_fieldvars, only : mb_vol,mb_vsl,mb_vst
    implicit none
    integer(kind_int), intent(in) :: nb,swp
    real(kind_real),   intent(in) :: cpr(3)
    integer(kind_int), intent(in) :: nsw
    real(kind_real),   intent(in) :: cdd
    integer(kind_int)          :: i,j,k,i1,j1,k1,m,n,ierr
    integer(kind_int)          :: st(3),ed(3),dijk
    real(kind_real)            :: srci,srcj,srck,srvi,srvj,srvk
    real(kind_real)            :: kt,kx,ky,kz,et,ex,ey,ez,ct,cx,cy,cz
    real(kind_real)            :: odt,dq1,dq2,rhs0(1:neqn)
    real(kind_real)            :: pvi(1:npvs),pvj(1:npvs),pvk(1:npvs)
    real(kind_real)            :: dqi(1:neqn),dqj(1:neqn),dqk(1:neqn)
    real(kind_real)            :: dfi(1:neqn),dfj(1:neqn),dfk(1:neqn)
    real(kind_real)            :: af0(neqn,neqn),afi(neqn,neqn)
    real(kind_real)            :: afj(neqn,neqn),afk(neqn,neqn)
    real(kind_real)            :: vol0,visl,vist
    type(fld_array_t), pointer :: vol(:),sxyz(:),pv(:),dt(:)
    type(fld_array_t), pointer :: qc(:),q0(:),rhs(:),qm(:)
    type(fld_array_t), pointer :: dq(:),src(:),srv(:)
    type(fld_array_t), pointer :: aiv(:),vsl(:),vst(:)

    vol  => mb_vol(nb)%fld
    sxyz => mb_sxyz(nb)%fld
    pv   => mb_pv(nb)%fld
    qc   => mb_qc(nb)%fld
    dt   => mb_dt(nb)%fld
    dq   => mb_dq(nb)%fld
    rhs  => mb_rhs(nb)%fld
    src  => mb_src(nb)%fld
    aiv  => mb_aiv(nb)%fld

    if (nvis > nvis_euler) then
        srv  => mb_srv(nb)%fld
        vsl  => mb_vsl(nb)%fld
        vst  => mb_vst(nb)%fld
    end if

    if (nsw > 0) then
        q0 => mb_q0(nb)%fld
        qm => mb_qm(nb)%fld
    end if

    if (swp > 0) then
        st(:) = mb_top(nb)%ndst(:)
        ed(:) = mb_top(nb)%nded(:)
        dijk  = 1
    else
        st(:) = mb_top(nb)%nded(:)
        ed(:) = mb_top(nb)%ndst(:)
        dijk  = -1
    end if

    do k=st(3),ed(3),dijk
    do j=st(2),ed(2),dijk
    do i=st(1),ed(1),dijk
        do m=1,neqn
            rhs0(m) = rhs(m)%r3d(i,j,k)
        end do

        if (nsw > 0) then
            odt = vol(1)%r3d(i,j,k)/dtime
            do m=1,neqn
                dq1 = qm(m)%r3d(i,j,k) - q0(m)%r3d(i,j,k)
                dq2 = qc(m)%r3d(i,j,k) - qm(m)%r3d(i,j,k)
                rhs0(m) = rhs0(m) - odt*(cpr(1)*dq2 - cpr(2)*dq1)
            end do
        end if

        i1 = i - 1
        j1 = j - 1
        k1 = k - 1
        !!srci = src(1)%r3d(i1,j,k)
        !!srcj = src(2)%r3d(i,j1,k)
        !!srck = src(3)%r3d(i,j,k1)

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
        call mxdq_sw(1,npvs,pvi,kt,kx,ky,kz,1,neqn,dqi,dfi,cdd,one)

        ex = sxyz(4)%r3d(i,j1,k)
        ey = sxyz(5)%r3d(i,j1,k)
        ez = sxyz(6)%r3d(i,j1,k)
        call mxdq_sw(1,npvs,pvj,et,ex,ey,ez,1,neqn,dqj,dfj,cdd,one)

        if (nsw_kdir == nsw_dir_close) then
            dfk(:) = zero
        else
            cx = sxyz(7)%r3d(i,j,k1)
            cy = sxyz(8)%r3d(i,j,k1)
            cz = sxyz(9)%r3d(i,j,k1)
            call mxdq_sw(1,npvs,pvk,ct,cx,cy,cz,1,neqn,dqk,dfk,cdd,one)
        end if

        rhs0(:) = rhs0(:) + (dfi(:) + dfj(:) + dfk(:))


        if (nvis > nvis_euler) then
#ifdef PR_VIS_JACOBIAN
            vol0 = vol(1)%r3d(i1,j,k)
            visl = vsl(1)%r3d(i1,j,k)
            vist = vst(1)%r3d(i1,j,k)
            call aflux_vis(1,npvs,pvi,kt,kx,ky,kz,vol0,visl,vist,1,neqn,afi)

            vol0 = vol(1)%r3d(i,j1,k)
            visl = vsl(1)%r3d(i,j1,k)
            vist = vst(1)%r3d(i,j1,k)
            call aflux_vis(1,npvs,pvj,et,ex,ey,ez,vol0,visl,vist,1,neqn,afj)

            if (nsw_kdir == nsw_dir_close) then
                afk(:,:) = zero
            else
                vol0 = vol(1)%r3d(i,j,k1)
                visl = vsl(1)%r3d(i,j,k1)
                vist = vst(1)%r3d(i,j,k1)
                call aflux_vis(1,npvs,pvk,ct,cx,cy,cz,vol0,visl,vist,1,neqn,afk)
            end if

            do m=1,neqn
                rhs0(m) = rhs0(m) + &
                          half*(sum(afi(m,:)*dqi(:)) + &
                                sum(afj(m,:)*dqj(:)) + &
                                sum(afk(m,:)*dqk(:)))
            end do
#else
            srvi = srv(1)%r3d(i1,j,k)
            srvj = srv(2)%r3d(i,j1,k)
            srvk = srv(3)%r3d(i,j,k1)
            do m=1,neqn
                rhs0(m) = rhs0(m) + &
                            half*(srvi*dqi(m) + srvj*dqj(m) + srvk*dqk(m))
            end do
#endif
        end if


        i1 = i + 1
        j1 = j + 1
        k1 = k + 1
        !!srci = src(1)%r3d(i1,j,k)
        !!srcj = src(2)%r3d(i,j1,k)
        !!srck = src(3)%r3d(i,j,k1)

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

        kx = half*(sxyz(1)%r3d(i,j,k) + sxyz(1)%r3d(i1,j,k))
        ky = half*(sxyz(2)%r3d(i,j,k) + sxyz(2)%r3d(i1,j,k))
        kz = half*(sxyz(3)%r3d(i,j,k) + sxyz(3)%r3d(i1,j,k))
        call mxdq_sw(1,npvs,pvi,kt,kx,ky,kz,1,neqn,dqi,dfi,cdd,-one)

        ex = half*(sxyz(4)%r3d(i,j,k) + sxyz(4)%r3d(i,j1,k))
        ey = half*(sxyz(5)%r3d(i,j,k) + sxyz(5)%r3d(i,j1,k))
        ez = half*(sxyz(6)%r3d(i,j,k) + sxyz(6)%r3d(i,j1,k))
        call mxdq_sw(1,npvs,pvj,et,ex,ey,ez,1,neqn,dqj,dfj,cdd,-one)

        if (nsw_kdir == nsw_dir_close) then
            dfk(:) = zero
        else
            cx = half*(sxyz(7)%r3d(i,j,k) + sxyz(7)%r3d(i,j,k1))
            cy = half*(sxyz(8)%r3d(i,j,k) + sxyz(8)%r3d(i,j,k1))
            cz = half*(sxyz(9)%r3d(i,j,k) + sxyz(9)%r3d(i,j,k1))
            call mxdq_sw(1,npvs,pvk,ct,cx,cy,cz,1,neqn,dqk,dfk,cdd,-one)
        end if

        rhs0(:) = rhs0(:) - (dfi(:) + dfj(:) + dfk(:))

        if (nvis > nvis_euler) then
#ifdef PR_VIS_JACOBIAN
            vol0 = vol(1)%r3d(i1,j,k)
            visl = vsl(1)%r3d(i1,j,k)
            vist = vst(1)%r3d(i1,j,k)
            call aflux_vis(1,npvs,pvi,kt,kx,ky,kz,vol0,visl,vist,1,neqn,afi)

            vol0 = vol(1)%r3d(i,j1,k)
            visl = vsl(1)%r3d(i,j1,k)
            vist = vst(1)%r3d(i,j1,k)
            call aflux_vis(1,npvs,pvj,et,ex,ey,ez,vol0,visl,vist,1,neqn,afj)

            if (nsw_kdir == nsw_dir_close) then
                afk(:,:) = zero
            else
                vol0 = vol(1)%r3d(i,j,k1)
                visl = vsl(1)%r3d(i,j,k1)
                vist = vst(1)%r3d(i,j,k1)
                call aflux_vis(1,npvs,pvk,ct,cx,cy,cz,vol0,visl,vist,1,neqn,afk)
            end if

            do m=1,neqn
                rhs0(m) = rhs0(m) + &
                          half*(sum(afi(m,:)*dqi(:)) + &
                                sum(afj(m,:)*dqj(:)) + &
                                sum(afk(m,:)*dqk(:)))
            end do
#else
            srvi = srv(1)%r3d(i1,j,k)
            srvj = srv(2)%r3d(i,j1,k)
            srvk = srv(3)%r3d(i,j,k1)
            do m=1,neqn
                rhs0(m) = rhs0(m) + &
                            half*(srvi*dqi(m) + srvj*dqj(m) + srvk*dqk(m))
            end do
#endif
        end if

        do m=1,neqn
        do n=1,neqn
            af0(m,n) = aiv((m-1)*neqn + n)%r3d(i,j,k)
        end do
        end do

        do m=1,neqn
            dq(m)%r3d(i,j,k) = sum(af0(m,:)*rhs0(:))
        end do
    end do
    end do
    end do

end subroutine prsgs_ust_mat_sweep

subroutine prsgs_ust_mat_sweep_sp(nb,swp,cpr,nsw,cdd)
    ! swp: 控制扫描方向（swp=0向前，否则向后）
    ! cpr: 时间离散系数
    ! nsw: 控制时间离散（nsw=0定常，否则非定常）
    ! cdd: Steger分裂修正系数
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,one,half,thr2nd
    use mod_constants, only : nvis_euler,nsw_dir_close
    use mod_datatypes, only : fld_array_t
    use mod_variables, only : nvis,dtime,nsw_kdir
    use mod_fieldvars, only : mb_topsp,mb_sxyzsp,mb_pv
    use mod_fieldvars, only : mb_dq,npvs,neqn,mb_dt
    use mod_fieldvars, only : mb_qc,mb_q0,mb_rhs,mb_qm
    use mod_fieldvars, only : mb_src,mb_srv,mb_aiv
    use mod_fieldvars, only : mb_volsp,mb_vsl,mb_vst
    implicit none
    integer(kind_int), intent(in) :: nb,swp
    real(kind_real),   intent(in) :: cpr(3)
    integer(kind_int), intent(in) :: nsw
    real(kind_real),   intent(in) :: cdd
    integer(kind_int)          :: i,j,k,i1,j1,k1,m,n,ierr
    integer(kind_int)          :: st(3),ed(3),dijk
    real(kind_real)            :: srci,srcj,srck,srvi,srvj,srvk
    real(kind_real)            :: kt,kx,ky,kz,et,ex,ey,ez,ct,cx,cy,cz
    real(kind_real)            :: odt,dq1,dq2,rhs0(1:neqn)
    real(kind_real)            :: pvi(1:npvs),pvj(1:npvs),pvk(1:npvs)
    real(kind_real)            :: dqi(1:neqn),dqj(1:neqn),dqk(1:neqn)
    real(kind_real)            :: dfi(1:neqn),dfj(1:neqn),dfk(1:neqn)
    real(kind_real)            :: af0(neqn,neqn),afi(neqn,neqn)
    real(kind_real)            :: afj(neqn,neqn),afk(neqn,neqn)
    real(kind_real)            :: vol0,visl,vist
    type(fld_array_t), pointer :: vol(:),sxyz(:),pv(:),dt(:)
    type(fld_array_t), pointer :: qc(:),q0(:),rhs(:),qm(:)
    type(fld_array_t), pointer :: dq(:),src(:),srv(:)
    type(fld_array_t), pointer :: aiv(:),vsl(:),vst(:)

    vol  => mb_volsp(nb)%fld
    sxyz => mb_sxyzsp(nb)%fld
    pv   => mb_pv(nb)%fld
    qc   => mb_qc(nb)%fld
    dt   => mb_dt(nb)%fld
    dq   => mb_dq(nb)%fld
    rhs  => mb_rhs(nb)%fld
    src  => mb_src(nb)%fld
    aiv  => mb_aiv(nb)%fld

    if (nvis > nvis_euler) then
        srv  => mb_srv(nb)%fld
        vsl  => mb_vsl(nb)%fld
        vst  => mb_vst(nb)%fld
    end if

    if (nsw > 0) then
        q0 => mb_q0(nb)%fld
        qm => mb_qm(nb)%fld
    end if

    if (swp > 0) then
        st(:) = mb_topsp(nb)%ndst(:)
        ed(:) = mb_topsp(nb)%nded(:)
        dijk  = 1
    else
        st(:) = mb_topsp(nb)%nded(:)
        ed(:) = mb_topsp(nb)%ndst(:)
        dijk  = -1
    end if

    do k=st(3),ed(3),dijk
    do j=st(2),ed(2),dijk
    do i=st(1),ed(1),dijk
        do m=1,neqn
            rhs0(m) = rhs(m)%r3d(i,j,k)
        end do

        if (nsw > 0) then
            odt = vol(1)%r3d(i,j,k)/dtime
            do m=1,neqn
                dq1 = qm(m)%r3d(i,j,k) - q0(m)%r3d(i,j,k)
                dq2 = qc(m)%r3d(i,j,k) - qm(m)%r3d(i,j,k)
                rhs0(m) = rhs0(m) - odt*(cpr(1)*dq2 - cpr(2)*dq1)
            end do
        end if

        i1 = i - 1
        j1 = j - 1
        k1 = k - 1
        !!srci = src(1)%r3d(i1,j,k)
        !!srcj = src(2)%r3d(i,j1,k)
        !!srck = src(3)%r3d(i,j,k1)

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
        call mxdq_sw(1,npvs,pvi,kt,kx,ky,kz,1,neqn,dqi,dfi,cdd,one)

        ex = sxyz(4)%r3d(i,j1,k)
        ey = sxyz(5)%r3d(i,j1,k)
        ez = sxyz(6)%r3d(i,j1,k)
        call mxdq_sw(1,npvs,pvj,et,ex,ey,ez,1,neqn,dqj,dfj,cdd,one)

        if (nsw_kdir == nsw_dir_close) then
            dfk(:) = zero
        else
            cx = sxyz(7)%r3d(i,j,k1)
            cy = sxyz(8)%r3d(i,j,k1)
            cz = sxyz(9)%r3d(i,j,k1)
            call mxdq_sw(1,npvs,pvk,ct,cx,cy,cz,1,neqn,dqk,dfk,cdd,one)
        end if

        rhs0(:) = rhs0(:) + (dfi(:) + dfj(:) + dfk(:))


        if (nvis > nvis_euler) then
#ifdef PR_VIS_JACOBIAN
            vol0 = vol(1)%r3d(i1,j,k)
            visl = vsl(1)%r3d(i1,j,k)
            vist = vst(1)%r3d(i1,j,k)
            call aflux_vis(1,npvs,pvi,kt,kx,ky,kz,vol0,visl,vist,1,neqn,afi)

            vol0 = vol(1)%r3d(i,j1,k)
            visl = vsl(1)%r3d(i,j1,k)
            vist = vst(1)%r3d(i,j1,k)
            call aflux_vis(1,npvs,pvj,et,ex,ey,ez,vol0,visl,vist,1,neqn,afj)

            if (nsw_kdir == nsw_dir_close) then
                afk(:,:) = zero
            else
                vol0 = vol(1)%r3d(i,j,k1)
                visl = vsl(1)%r3d(i,j,k1)
                vist = vst(1)%r3d(i,j,k1)
                call aflux_vis(1,npvs,pvk,ct,cx,cy,cz,vol0,visl,vist,1,neqn,afk)
            end if

            do m=1,neqn
                rhs0(m) = rhs0(m) + &
                          half*(sum(afi(m,:)*dqi(:)) + &
                                sum(afj(m,:)*dqj(:)) + &
                                sum(afk(m,:)*dqk(:)))
            end do
#else
            srvi = srv(1)%r3d(i1,j,k)
            srvj = srv(2)%r3d(i,j1,k)
            srvk = srv(3)%r3d(i,j,k1)
            do m=1,neqn
                rhs0(m) = rhs0(m) + &
                            half*(srvi*dqi(m) + srvj*dqj(m) + srvk*dqk(m))
            end do
#endif
        end if


        i1 = i + 1
        j1 = j + 1
        k1 = k + 1
        !!srci = src(1)%r3d(i1,j,k)
        !!srcj = src(2)%r3d(i,j1,k)
        !!srck = src(3)%r3d(i,j,k1)

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

        kx = half*(sxyz(1)%r3d(i,j,k) + sxyz(1)%r3d(i1,j,k))
        ky = half*(sxyz(2)%r3d(i,j,k) + sxyz(2)%r3d(i1,j,k))
        kz = half*(sxyz(3)%r3d(i,j,k) + sxyz(3)%r3d(i1,j,k))
        call mxdq_sw(1,npvs,pvi,kt,kx,ky,kz,1,neqn,dqi,dfi,cdd,-one)

        ex = half*(sxyz(4)%r3d(i,j,k) + sxyz(4)%r3d(i,j1,k))
        ey = half*(sxyz(5)%r3d(i,j,k) + sxyz(5)%r3d(i,j1,k))
        ez = half*(sxyz(6)%r3d(i,j,k) + sxyz(6)%r3d(i,j1,k))
        call mxdq_sw(1,npvs,pvj,et,ex,ey,ez,1,neqn,dqj,dfj,cdd,-one)

        if (nsw_kdir == nsw_dir_close) then
            dfk(:) = zero
        else
            cx = half*(sxyz(7)%r3d(i,j,k) + sxyz(7)%r3d(i,j,k1))
            cy = half*(sxyz(8)%r3d(i,j,k) + sxyz(8)%r3d(i,j,k1))
            cz = half*(sxyz(9)%r3d(i,j,k) + sxyz(9)%r3d(i,j,k1))
            call mxdq_sw(1,npvs,pvk,ct,cx,cy,cz,1,neqn,dqk,dfk,cdd,-one)
        end if

        rhs0(:) = rhs0(:) - (dfi(:) + dfj(:) + dfk(:))

        if (nvis > nvis_euler) then
#ifdef PR_VIS_JACOBIAN
            vol0 = vol(1)%r3d(i1,j,k)
            visl = vsl(1)%r3d(i1,j,k)
            vist = vst(1)%r3d(i1,j,k)
            call aflux_vis(1,npvs,pvi,kt,kx,ky,kz,vol0,visl,vist,1,neqn,afi)

            vol0 = vol(1)%r3d(i,j1,k)
            visl = vsl(1)%r3d(i,j1,k)
            vist = vst(1)%r3d(i,j1,k)
            call aflux_vis(1,npvs,pvj,et,ex,ey,ez,vol0,visl,vist,1,neqn,afj)

            if (nsw_kdir == nsw_dir_close) then
                afk(:,:) = zero
            else
                vol0 = vol(1)%r3d(i,j,k1)
                visl = vsl(1)%r3d(i,j,k1)
                vist = vst(1)%r3d(i,j,k1)
                call aflux_vis(1,npvs,pvk,ct,cx,cy,cz,vol0,visl,vist,1,neqn,afk)
            end if

            do m=1,neqn
                rhs0(m) = rhs0(m) + &
                          half*(sum(afi(m,:)*dqi(:)) + &
                                sum(afj(m,:)*dqj(:)) + &
                                sum(afk(m,:)*dqk(:)))
            end do
#else
            srvi = srv(1)%r3d(i1,j,k)
            srvj = srv(2)%r3d(i,j1,k)
            srvk = srv(3)%r3d(i,j,k1)
            do m=1,neqn
                rhs0(m) = rhs0(m) + &
                            half*(srvi*dqi(m) + srvj*dqj(m) + srvk*dqk(m))
            end do
#endif
        end if

        do m=1,neqn
        do n=1,neqn
            af0(m,n) = aiv((m-1)*neqn + n)%r3d(i,j,k)
        end do
        end do

        do m=1,neqn
            dq(m)%r3d(i,j,k) = sum(af0(m,:)*rhs0(:))
        end do
    end do
    end do
    end do

end subroutine prsgs_ust_mat_sweep_sp

subroutine brinv(n,a,l)  !> 矩阵求逆
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,one
    implicit none
    integer(kind_int),  intent(in) :: n
    real(kind_real), intent(inout) :: a(n,n)
    integer(kind_int), intent(out) :: l
    real(kind_real)   :: t,d
    integer(kind_int) :: i,j,k
    integer(kind_int) :: is(n),js(n)

    l = 1
    do k=1,n
        d = zero
        do i=k,n
        do j=k,n
            if (abs(a(i,j)) > d) then
              d = abs(a(i,j))
              is(k) = i
              js(k) = j
            end if
        end do
        end do

        if (d+one == one) then
            l = 0
            write(*,*)' err**not inv'
            exit
        end if

        do j=1,n
            t = a(k,j)
            a(k,j) = a(is(k),j)
            a(is(k),j) = t
        end do

        do i=1,n
          t = a(i,k)
          a(i,k) = a(i,js(k))
          a(i,js(k)) = t
        end do

        a(k,k) = one/a(k,k)
        do j=1,n
            if (j /= k) then
                a(k,j) = a(k,j)*a(k,k)
            end if
        end do

        do i=1,n
            if (i /= k) then
                do j=1,n
                    if (j /= k) then
                        a(i,j)=a(i,j)-a(i,k)*a(k,j)
                    end if
                end do
            end if
        end do

        do i=1,n
            if (i /= k) then
                a(i,k) = -a(i,k)*a(k,k)
            end if
        end do
    end do

    if (l == 1) then
        do k=n,1,-1
            do j=1,n
                t = a(k,j)
                a(k,j) = a(js(k),j)
                a(js(k),j) = t
            end do
            do i=1,n
                t = a(i,k)
                a(i,k) = a(i,is(k))
                a(i,is(k)) = t
            end do
        end do
    end if

end subroutine brinv

subroutine stop_subiter(iter,nstop)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,one,small
    use mod_constants, only : nlhs_prsgs_ust_sca
    use mod_constants, only : nlhs_prsgs_std_sca
    use mod_constants, only : nlhs_prsgs_ust_mat
    use mod_constants, only : nlhs_prsgs_std_mat
    use mod_constants, only : nlhs_prsgs_ust_sca_scmp
    use mod_variables, only : nlhs,nsubmax,tolsub
    use mod_variables, only : rcflmax,dtime,nstep,nressav
    use mod_variables, only : restot,resdts,nsub,nsubp,nsubp1
    use mod_variables, only : restp0,restpn,restpm
    use mod_fieldvars, only : ntotpts
    use mod_parallels
    implicit none
    integer(kind_int), intent(in)  :: iter
    integer(kind_int), intent(out) :: nstop
    integer(kind_int) :: m,ierr
    real(kind_real)   :: resave,tol,rcflmax0
    real(kind_real)   :: resave_dts

    resave = sqrt(restot/ntotpts)

    if (iter == 1) then
        restp0 = resave
        restpn = resave
        restpm = resave
        tol = one

        call calc_cfl_realtime
#ifdef PARALLEL
        call MPI_REDUCE(rcflmax,rcflmax0,1,kind_real_mpi,MPI_MAX,master,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(rcflmax0,1,kind_real_mpi,master,MPI_COMM_WORLD,ierr)
        rcflmax = rcflmax0
#endif
    else
        select case(nlhs)
        case(nlhs_prsgs_std_mat,nlhs_prsgs_std_sca)
            restpn = restpm
            restpm = resave
            tol = abs(restpm-restpn)/(restpm + small)
        case(nlhs_prsgs_ust_mat,nlhs_prsgs_ust_sca,nlhs_prsgs_ust_sca_scmp)
            resave_dts = sqrt(resdts/ntotpts)

            restpn = resave_dts
            restpm = resave
            tol = restpn/(restpm + small)
        end select
    end if
    
    if (tol < tolsub .or. iter == nsubmax) then
        nstop = 1
        nsubp1 = nsubp
        nsubp  = nsub
        nsub   = iter
    else
        nstop = 0
    end if

#ifdef PARALLEL
    if (myid == master) then
#endif
        if (nstop /= 0 .and. mod(nstep,nressav) == 0) then
            write(*,20)iter,rcflmax,dtime,restpn,restpm,tol
        end if
#ifdef PARALLEL
    end if
#endif

20  format(7x,"->",i3,2(1x,e9.2),3(1x,e12.5))

end subroutine stop_subiter

subroutine stop_subiter_sp(iter,nstop)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,one,small
    use mod_constants, only : nlhs_prsgs_ust_sca
    use mod_constants, only : nlhs_prsgs_std_sca
    use mod_constants, only : nlhs_prsgs_ust_mat
    use mod_constants, only : nlhs_prsgs_std_mat
    use mod_constants, only : nlhs_prsgs_ust_sca_scmp
    use mod_variables, only : nlhs,nsubmax,tolsub
    use mod_variables, only : rcflmax,dtime,nstep,nressav
    use mod_variables, only : restot,resdts,nsub,nsubp,nsubp1
    use mod_variables, only : restp0,restpn,restpm
    use mod_fieldvars, only : ntotpts
    use mod_parallels
    implicit none
    integer(kind_int), intent(in)  :: iter
    integer(kind_int), intent(out) :: nstop
    integer(kind_int) :: m,ierr
    real(kind_real)   :: resave,tol,rcflmax0
    real(kind_real)   :: resave_dts

    resave = sqrt(restot/ntotpts)

    if (iter == 1) then
        restp0 = resave
        restpn = resave
        restpm = resave
        tol = one

        call calc_cfl_realtime_sp
#ifdef PARALLEL
        call MPI_REDUCE(rcflmax,rcflmax0,1,kind_real_mpi,MPI_MAX,master,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(rcflmax0,1,kind_real_mpi,master,MPI_COMM_WORLD,ierr)
        rcflmax = rcflmax0
#endif
    else
        select case(nlhs)
        case(nlhs_prsgs_std_mat,nlhs_prsgs_std_sca)
            restpn = restpm
            restpm = resave
            tol = abs(restpm-restpn)/(restpm + small)
        case(nlhs_prsgs_ust_mat,nlhs_prsgs_ust_sca,nlhs_prsgs_ust_sca_scmp)
            resave_dts = sqrt(resdts/ntotpts)

            restpn = resave_dts
            restpm = resave
            tol = restpn/(restpm + small)
        end select
    end if
    
    if (tol < tolsub .or. iter == nsubmax) then
        nstop = 1
        nsubp1 = nsubp
        nsubp  = nsub
        nsub   = iter
    else
        nstop = 0
    end if

#ifdef PARALLEL
    if (myid == master) then
#endif
        if (nstop /= 0 .and. mod(nstep,nressav) == 0) then
            write(*,20)iter,rcflmax,dtime,restpn,restpm,tol
        end if
#ifdef PARALLEL
    end if
#endif

20  format(7x,"->",i3,2(1x,e9.2),3(1x,e12.5))

end subroutine stop_subiter_sp

subroutine calc_cfl_realtime
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : small,large,one
    use mod_datatypes, only : fld_array_t,top_block_t
    use mod_variables, only : rcflmin,rcflmax,dtime
    use mod_fieldvars, only : nblkcoms,blkcoms
    use mod_fieldvars, only : mb_vol,mb_rdt
    implicit none
    integer(kind_int)          :: nc,nb,i,j,k,ierr
    integer(kind_int)          :: ni,nj,nk
    type(top_block_t), pointer :: top
    real(kind_real)            :: srad,vol0,dt0,cfl0
    type(fld_array_t), pointer :: vol(:),rdt(:)

    rcflmin = large
    rcflmax = small
    do nc=1,nblkcoms
        nb  =  blkcoms(nc)%nb
        top => blkcoms(nc)%top

        ni = top%nijk(1)
        nj = top%nijk(2)
        nk = top%nijk(3)

        vol => mb_vol(nb)%fld
        rdt => mb_rdt(nb)%fld

        do k=1,nk
        do j=1,nj
        do i=1,ni
            srad = rdt(1)%r3d(i,j,k)
            vol0 = vol(1)%r3d(i,j,k)
            dt0  = dtime/vol0

            cfl0 = dt0*srad
            rcflmin = min(rcflmin,cfl0)
            rcflmax = max(rcflmax,cfl0)
        end do
        end do
        end do
    end do

end subroutine calc_cfl_realtime

subroutine calc_cfl_realtime_sp
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : small,large,one
    use mod_datatypes, only : fld_array_t,top_block_t
    use mod_variables, only : rcflmin,rcflmax,dtime
    use mod_fieldvars, only : nblkcoms,blkcomssp
    use mod_fieldvars, only : mb_volsp,mb_rdt
    implicit none
    integer(kind_int)          :: nc,nb,i,j,k,ierr
    integer(kind_int)          :: ni,nj,nk
    type(top_block_t), pointer :: top
    real(kind_real)            :: srad,vol0,dt0,cfl0
    type(fld_array_t), pointer :: vol(:),rdt(:)

    rcflmin = large
    rcflmax = small
    do nc=1,nblkcoms
        nb  =  blkcomssp(nc)%nb
        top => blkcomssp(nc)%top

        ni = top%nijk(1)
        nj = top%nijk(2)
        nk = top%nijk(3)

        vol => mb_volsp(nb)%fld
        rdt => mb_rdt(nb)%fld

        do k=1,nk
        do j=1,nj
        do i=1,ni
            srad = rdt(1)%r3d(i,j,k)
            vol0 = vol(1)%r3d(i,j,k)
            dt0  = dtime/vol0

            cfl0 = dt0*srad
            rcflmin = min(rcflmin,cfl0)
            rcflmax = max(rcflmax,cfl0)
        end do
        end do
        end do
    end do

end subroutine calc_cfl_realtime_sp
