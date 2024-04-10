
subroutine solve
    use mod_constants, only : zero,one
    use mod_constants, only : nrestrt_restart
    use mod_constants, only : nlhs_rkutta_3step
    use mod_constants, only : nlhs_lusgs_std_sca,nlhs_lusgs_std_prec
    use mod_constants, only : nlhs_prsgs_ust_sca,nlhs_prsgs_std_sca
    use mod_constants, only : nlhs_prsgs_ust_mat,nlhs_prsgs_std_mat
    use mod_constants, only : nlhs_lusgs_std_scmp,nlhs_prsgs_ust_sca_scmp
    use mod_variables, only : nlhs,nstep,nstepsav,nacous
    use mod_variables, only : nstepst,nsteped,nprms,ns_mean,ns_prms
    use mod_variables, only : nrestrt,nstepmx,nsolsav,nmoni
    use mod_variables, only : nressav,nfcesav,npltsav  
    use mod_variables, only : cflst,cfled,cflfac,cfl
    use mod_variables, only : ncelcet
    use mod_fieldvars, only : nstepmean,nsteprms,ntime
    use mod_runtimers
    implicit none

    select case(nrestrt)
    case(nrestrt_restart)
        time_cpu_sav  = zero
        time_wall_sav = zero

        nstepst = 1
        nsteped = nstepmx

        cfl = cflst
    case default
        nstepst = 1 + nstepsav
        nsteped = nstepmx

        cfl = min(max(cfl,cflst),cfled)
    end select

    call start_timer


    do nstep=nstepst,nsteped
        if (cflfac > one) then
           if (cfl < cfled) then
              cfl = min(cfl*cflfac,cfled)
           else
              cfl = max(cfl/cflfac,cfled)
           end if
        else
           cfl = cfled
        end if
        
        !todo SCM ����һ���µ�ʱ���ʽ
        if(ncelcet == 1) then
            select case(nlhs)
            case(nlhs_rkutta_3step)
                call sol_rkutta_3step_sp
            case(nlhs_lusgs_std_sca)
                call sol_lusgs_std_sca_sp
            case(nlhs_prsgs_ust_sca)
                call sol_prsgs_ust_sca_sp
            case(nlhs_prsgs_std_sca)
                call sol_prsgs_std_sca_sp
            case(nlhs_prsgs_ust_mat)
                call sol_prsgs_ust_mat_sp
            case(nlhs_prsgs_std_mat)
                call sol_prsgs_std_mat_sp
            case(nlhs_lusgs_std_prec)
                call sol_lusgs_std_prec_sp
            case(nlhs_lusgs_std_scmp)
                call sol_lusgs_std_scmp_sp
            case(nlhs_prsgs_ust_sca_scmp)
                call sol_prsgs_ust_sca_scmp_sp    
            end select
            
            if (nacous > 0 .and. nstep > ns_mean) then
                ntime = ntime+1
                call to_noise_write
                call to_noise_output_sp
            end if
            
            if (nmoni>0) then
                if (mod(nstep,nressav) == 0) call monitor_output_sp
            end if
            
            !todo call aeroforce5_sp
            
            if (nprms > 0) then
                if (nstep > ns_mean .and. nstep <= ns_prms) then
                    call flow_mean_sp
                    call solve_fmean_sp
                end if

                if (nstep > ns_prms) then
                    call flow_rms_sp
                    call solve_frms_sp
                end if
            end if
            
            call get_run_timer
            
            if (mod(nstep,nressav) == 0) call output_res
            if (mod(nstep,nfcesav) == 0) call output_fce
            if (mod(nstep+nsolsav,nsolsav) == 0) call output_sol_sp
            if (mod(nstep+npltsav,npltsav) == 0) call output_plt_sp
        else
            select case(nlhs)
            case(nlhs_rkutta_3step)
                call sol_rkutta_3step
            case(nlhs_lusgs_std_sca)
                call sol_lusgs_std_sca
            case(nlhs_prsgs_ust_sca)
                call sol_prsgs_ust_sca
            case(nlhs_prsgs_std_sca)
                call sol_prsgs_std_sca
            case(nlhs_prsgs_ust_mat)
                call sol_prsgs_ust_mat
            case(nlhs_prsgs_std_mat)
                call sol_prsgs_std_mat
            case(nlhs_lusgs_std_prec)
                call sol_lusgs_std_prec
            case(nlhs_lusgs_std_scmp)
                call sol_lusgs_std_scmp
            case(nlhs_prsgs_ust_sca_scmp)
                call sol_prsgs_ust_sca_scmp
            end select  
            
            if (nacous > 0 .and. nstep > ns_mean) then
                ntime = ntime+1
                call to_noise_write
                call to_noise_output
            end if
        
            if (nmoni>0) then
                if (mod(nstep,nressav) == 0) call monitor_output
            end if

            call aeroforce5

            if (nprms > 0) then
                if (nstep > ns_mean .and. nstep <= ns_prms) then
                    call flow_mean
                    call solve_fmean
                end if

                if (nstep > ns_prms) then
                    call flow_rms
                    call solve_frms
                end if
            end if

            call get_run_timer

            if (mod(nstep,nressav) == 0) call output_res
            if (mod(nstep,nfcesav) == 0) call output_fce
            if (mod(nstep+nsolsav,nsolsav) == 0) call output_sol
            if (mod(nstep+npltsav,npltsav) == 0) call output_plt
        
            if (nprms > 0) then
                if (nstep > ns_mean .and. nstep <= ns_prms) then
                    if (mod(nstep+nsolsav,nsolsav) == 0 ) call output_mean
                end if            
                if (nstep == ns_prms) call output_fmean

                if (nstep > ns_prms) then
                    if (mod(nstep+nsolsav,nsolsav) == 0 ) call output_rms
                end if
                if (nstep == nstepmx ) call output_frms
            end if
        end if
        
    end do

end subroutine solve

subroutine prepare_linear_system
    use mod_constants, only : nsgl_buffer_pvs
    use mod_constants, only : nsgl_aver_art,nbc_inter_buf_pvs
    use mod_constants, only : nprec_non,nscmp_non
    use mod_fieldvars, only : npvs,mb_pv
    use mod_variables, only : nghnode
    use mod_variables, only : nprec,nscmp
    use mod_variables, only : dtime
    use mod_interface, only : pre_exchange_bc_var,post_exchange_bc_var,average_bc_var
    use mod_interface, only : patched_ghost_points,exchange_singulars,average_singulars
    implicit none

    call boundary_conditions

    call pre_exchange_bc_var(mb_pv,1,npvs,nghnode,nbc_inter_buf_pvs,nsgl_aver_art)
    call exchange_singulars(mb_pv,1,npvs,nsgl_buffer_pvs,nsgl_aver_art)
    call post_exchange_bc_var(mb_pv,1,npvs,nghnode,nbc_inter_buf_pvs,nsgl_aver_art)
    call patched_ghost_points(mb_pv,1,npvs,nghnode)
    call average_bc_var(mb_pv,1,npvs,nbc_inter_buf_pvs,nsgl_aver_art)
    call average_singulars(mb_pv,1,npvs,nsgl_buffer_pvs,nsgl_aver_art)

    call fill_corner_mb_pv

    call calc_other_variables

    !Ԥ��������£���Ҫ����ԭʼ�����ſɱȾ����װ뾶
    if (nprec > nprec_non) then
        call calc_spectral_radius
        call calc_spectral_radius_prec
    else if (nscmp > nscmp_non) then !< todo SCM-P ��Ҫ�����µĶ����ſɱȾ����װ뾶
        call calc_spectral_radius_scmp
    else
        call calc_spectral_radius
    end if

    call calc_time_step
    
    !if (nprec > nprec_non) call calc_diag_inv_matrix

    call rhside

    !call interface_C    

    !call boundary_conditions_C

end subroutine prepare_linear_system

subroutine prepare_linear_system_sp
    use mod_constants, only : nsgl_aver_art,nbc_inter_buf_pvs
    use mod_constants, only : nprec_non,nscmp_non
    use mod_fieldvars, only : npvs,mb_pv
    use mod_variables, only : nghnode
    use mod_variables, only : nprec,nscmp
    use mod_interface, only : pre_exchange_bc_var_sp,post_exchange_bc_var_sp
    implicit none

    !! to do call boundary_conditions

    call pre_exchange_bc_var_sp(mb_pv,1,npvs,nghnode,nbc_inter_buf_pvs,nsgl_aver_art)
    call post_exchange_bc_var_sp(mb_pv,1,npvs,nghnode,nbc_inter_buf_pvs,nsgl_aver_art)

    !call fill_corner_mb_pv_sp

    call calc_other_variables_sp

    !Ԥ��������£���Ҫ����ԭʼ�����ſɱȾ����װ뾶
    if (nprec > nprec_non) then
        call calc_spectral_radius_sp
        call calc_spectral_radius_prec_sp
    else if (nscmp > nscmp_non) then !< todo SCM-P ��Ҫ�����µĶ����ſɱȾ����װ뾶
        call calc_spectral_radius_scmp_sp
    else
        call calc_spectral_radius_sp
    end if

    call calc_time_step_sp
    
    !if (nprec > nprec_non) call calc_diag_inv_matrix

    call rhside_sp

    !call interface_C    

    !call boundary_conditions_C

end subroutine prepare_linear_system_sp

subroutine calc_other_variables
    use mod_datatypes, only : var_block_t,kind_int
    use mod_constants, only : nvis_euler,nsgl_buffer_pvs
    use mod_constants, only : nsgl_aver_art,nbc_inter_buf_pvs
    use mod_variables, only : nvis,nghnode
    use mod_fieldvars, only : npvs,mb_pv,neqn,mb_qc
    use mod_fieldvars, only : mb_t,mb_c,mb_vsl,mb_vst
    use mod_interface, only : calc_bc_var_via_sub,calc_mb_var_via_sub
    use mod_interface, only : mb_var_pointer_create,mb_var_pointer_assign
    use mod_interface, only : mb_var_pointer_delete
    
    use mod_constants, only : nscmp_non     !> Done
    use mod_variables, only : nscmp         !> Done
    
    implicit none
    external                   :: prim2con,prim2tc,t2visl
    external                   :: prim2tc_scmp,t2visl_scmp
    type(var_block_t), pointer :: mb_rp(:),mb_tc(:)
    integer(kind_int)          :: m

    call calc_bc_var_via_sub(mb_pv,1,npvs,prim2con,mb_qc,1,neqn,0,nghnode)

    call mb_var_pointer_create(mb_rp,1,2)
    call mb_var_pointer_assign(mb_rp,1,1,mb_pv,1,1)
    call mb_var_pointer_assign(mb_rp,2,2,mb_pv,5,5)

    call mb_var_pointer_create(mb_tc,1,2)
    call mb_var_pointer_assign(mb_tc,1,1,mb_t,1,1)
    call mb_var_pointer_assign(mb_tc,2,2,mb_c,1,1)

    if (nscmp > nscmp_non) then
        call calc_mb_var_via_sub(mb_rp,1,2,prim2tc_SCMP,mb_tc,1,2,nghnode)  !> Done
    else
        call calc_mb_var_via_sub(mb_rp,1,2,prim2tc     ,mb_tc,1,2,nghnode)  !> Done
    end if

    call mb_var_pointer_delete(mb_rp)
    call mb_var_pointer_delete(mb_tc)

    if (nvis > nvis_euler) then
        if (nscmp > nscmp_non) then
            call calc_mb_var_via_sub(mb_t,1,1,t2visl_SCMP,mb_vsl,1,1,nghnode)  !> Done ճ��ϵ��ȡ1.0
        else
            call calc_mb_var_via_sub(mb_t,1,1,t2visl     ,mb_vsl,1,1,nghnode)
        end if
    end if


end subroutine calc_other_variables

subroutine calc_other_variables_sp
    use mod_datatypes, only : var_block_t,kind_int
    use mod_constants, only : nvis_euler,nsgl_buffer_pvs
    use mod_constants, only : nsgl_aver_art,nbc_inter_buf_pvs
    use mod_variables, only : nvis,nghnode
    use mod_fieldvars, only : npvs,mb_pv,neqn,mb_qc
    use mod_fieldvars, only : mb_t,mb_c,mb_vsl,mb_vst
    use mod_interface, only : calc_bc_var_via_sub_sp,calc_mb_var_via_sub_sp
    use mod_interface, only : mb_var_pointer_create,mb_var_pointer_assign
    use mod_interface, only : mb_var_pointer_delete
    
    use mod_constants, only : nscmp_non     !> Done
    use mod_variables, only : nscmp         !> Done
    
    implicit none
    external                   :: prim2con,prim2tc,t2visl
    external                   :: prim2tc_scmp,t2visl_scmp
    type(var_block_t), pointer :: mb_rp(:),mb_tc(:)
    integer(kind_int)          :: m

    call calc_bc_var_via_sub_sp(mb_pv,1,npvs,prim2con,mb_qc,1,neqn,0,nghnode)

    call mb_var_pointer_create(mb_rp,1,2)
    call mb_var_pointer_assign(mb_rp,1,1,mb_pv,1,1)
    call mb_var_pointer_assign(mb_rp,2,2,mb_pv,5,5)

    call mb_var_pointer_create(mb_tc,1,2)
    call mb_var_pointer_assign(mb_tc,1,1,mb_t,1,1)
    call mb_var_pointer_assign(mb_tc,2,2,mb_c,1,1)

    if (nscmp > nscmp_non) then
        call calc_mb_var_via_sub_sp(mb_rp,1,2,prim2tc_SCMP,mb_tc,1,2,nghnode)  !> Done
    else
        call calc_mb_var_via_sub_sp(mb_rp,1,2,prim2tc     ,mb_tc,1,2,nghnode)  !> Done
    end if

    call mb_var_pointer_delete(mb_rp)
    call mb_var_pointer_delete(mb_tc)

    if (nvis > nvis_euler) then
        if (nscmp > nscmp_non) then
            call calc_mb_var_via_sub_sp(mb_t,1,1,t2visl_SCMP,mb_vsl,1,1,nghnode)  !> Done ճ��ϵ��ȡ1.0
        else
            call calc_mb_var_via_sub_sp(mb_t,1,1,t2visl     ,mb_vsl,1,1,nghnode)
        end if
    end if


end subroutine calc_other_variables_sp

!todoԤ������� ����Ԥ��������Ͷ�Ӧ���װ뾶    
subroutine calc_spectral_radius_prec
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : one,two,four3rd,nvis_euler
    use mod_datatypes, only : fld_array_t,top_block_t
    use mod_variables, only : nvis,reue,csrvis,fsw_kdir
    use mod_variables, only : gamma,prlam,prtur,refbeta
    use mod_fieldvars, only : nblkcoms,blkcoms
    use mod_fieldvars, only : mb_pv,mb_sxyz,mb_vol
    use mod_fieldvars, only : mb_vsl,mb_vst
    !use mod_fieldvars, only : mb_nonprec,mb_prec_diag_inv
    use mod_fieldvars, only : mb_c,mb_src,mb_src_prec,mb_srv,mb_rdt
    use mod_openmp
    implicit none
    integer(kind_int)          :: nc,nb,i,j,k,ierr
    integer(kind_int)          :: ni,nj,nk
    type(top_block_t), pointer :: top
    real(kind_real)            :: vx,vy,vz,kx,ky,kz
    real(kind_real)            :: ex,ey,ez,cx,cy,cz,sw2d
    real(kind_real)            :: vnk,vne,vnc,snk,sne,snc
    real(kind_real)            :: rca,rcb,rcc,rca_prec,rcb_prec,rcc_prec
    real(kind_real)            :: a,kcp,rov,coe,gama,ae,cv
    real(kind_real)            :: cp,cp_prl,cp_prt,c5
    type(fld_array_t), pointer :: pv(:),sxyz(:),vol(:),vsl(:),vst(:)
    type(fld_array_t), pointer :: c(:)
    type(fld_array_t), pointer :: src(:),src_prec(:),srv(:),rdt(:)
    !type(fld_array_t), pointer :: nonprec(:)
    
    !todo Ԥ�������
    integer(kind_int)          :: oned_index,l,m,s,nsp,nep
    real(kind_real)            :: ro,p,sn(1:3,1:3),prim(1:5)
    real(kind_real)            :: nonpreMatrix(1:4)
    real(kind_real)            :: eigenValue(1:3,1:5)
    real(kind_real)            :: sumSpectralRadius,length,visl,vist,vis

    nsp = 1
    nep = 5
    gama = gamma
    ae = gama - one

    cv = refbeta/ae
    cp = gama*cv
    cp_prl = cp/prlam
    cp_prt = cp/prtur
    
    vis = 0.
    length = 1.

    do nc=1,nblkcoms
        nb  =  blkcoms(nc)%nb
        top => blkcoms(nc)%top

        ni = top%nijk(1)
        nj = top%nijk(2)
        nk = top%nijk(3)

        pv   => mb_pv(nb)%fld
        sxyz => mb_sxyz(nb)%fld
        vol  => mb_vol(nb)%fld
        c    => mb_c(nb)%fld
        src  => mb_src(nb)%fld                !<Ac�����װ뾶
        src_prec  => mb_src_prec(nb)%fld      !<����Ԥ������Ķ����ſɱȾ����װ뾶
        rdt  => mb_rdt(nb)%fld
        !nonprec => mb_nonprec(nb)%fld

        if (nvis > nvis_euler) then
            vsl => mb_vsl(nb)%fld
            vst => mb_vst(nb)%fld
            srv => mb_srv(nb)%fld
        end if
!$OMP parallel  private(k,j,i,kx,ky,kz,ex,ey,ez,cx,cy,cz,vx,vy,vz,vnk,vne,vnc,snk,sne,snc,a, &
!$OMP                   rca,rcb,rcc,rov,visl,vist,vis,kcp,c5,coe,rva,rvb,rvc)
!$OMP do
        do k=0,nk+1
        do j=0,nj+1
        do i=0,ni+1
            
            kx = sxyz(1)%r3d(i,j,k)
            ky = sxyz(2)%r3d(i,j,k)
            kz = sxyz(3)%r3d(i,j,k)
            
            do m = 1,3
                sn(1,m) = sxyz(m)%r3d(i,j,k)
            end do

            ex = sxyz(4)%r3d(i,j,k)
            ey = sxyz(5)%r3d(i,j,k)
            ez = sxyz(6)%r3d(i,j,k)
            
            do m = 1,3
                sn(2,m) = sxyz(m+3)%r3d(i,j,k)
            end do

            cx = sxyz(7)%r3d(i,j,k)
            cy = sxyz(8)%r3d(i,j,k)
            cz = sxyz(9)%r3d(i,j,k)
            
            do m = 1,3
                sn(3,m) = sxyz(m+6)%r3d(i,j,k)
            end do

            ro = pv(1)%r3d(i,j,k)
            vx = pv(2)%r3d(i,j,k)
            vy = pv(3)%r3d(i,j,k)
            vz = pv(4)%r3d(i,j,k)
            p  = pv(5)%r3d(i,j,k)
            
            do m = 1,5
                prim(m) = pv(m)%r3d(i,j,k)
            end do
            
            if (nvis > nvis_euler) then
                visl = vsl(1)%r3d(i,j,k)
                vist = vst(1)%r3d(i,j,k)
                vis = visl + vist
                length = vol(1)%r3d(i,j,k)**(1.0/3.0)
            end if
            
            call calculate_precondition_eigenvalue_point(nsp,nep,sn,prim,vis,length,eigenvalue)
            
            !nonprec(:)%r3d(i,j,k) = nonpreMatrix(:)

            rca_prec = abs(eigenvalue(1,4))
            rcb_prec = abs(eigenvalue(2,4))
            rcc_prec = abs(eigenvalue(3,4))
            src_prec(1)%r3d(i,j,k) = rca_prec
            src_prec(2)%r3d(i,j,k) = rcb_prec
            src_prec(3)%r3d(i,j,k) = rcc_prec*fsw_kdir
            
            rdt(1)%r3d(i,j,k) = rca_prec + rcb_prec + rcc_prec*fsw_kdir
            
            !���ڼ���ʱ�䲽�����װ뾶
            if (nvis > nvis_euler) then
                rdt(1)%r3d(i,j,k) = rdt(1)%r3d(i,j,k) +srv(1)%r3d(i,j,k) + srv(2)%r3d(i,j,k) + srv(3)%r3d(i,j,k)
            end if
            
        end do
        end do
        end do
!$OMP end  do
!$OMP end parallel
    end do

end subroutine calc_spectral_radius_prec

!todoԤ������� ����Ԥ��������Ͷ�Ӧ���װ뾶    
subroutine calc_spectral_radius_prec_sp
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : one,two,four3rd,nvis_euler
    use mod_datatypes, only : fld_array_t,top_block_t
    use mod_variables, only : nvis,reue,csrvis,fsw_kdir
    use mod_variables, only : gamma,prlam,prtur,refbeta
    use mod_fieldvars, only : nblkcoms,blkcomssp
    use mod_fieldvars, only : mb_pv,mb_sxyzsp,mb_volsp
    use mod_fieldvars, only : mb_vsl,mb_vst
    !use mod_fieldvars, only : mb_nonprec,mb_prec_diag_inv
    use mod_fieldvars, only : mb_c,mb_src,mb_src_prec,mb_srv,mb_rdt
    use mod_openmp
    implicit none
    integer(kind_int)          :: nc,nb,i,j,k,ierr
    integer(kind_int)          :: ni,nj,nk
    type(top_block_t), pointer :: top
    real(kind_real)            :: vx,vy,vz,kx,ky,kz
    real(kind_real)            :: ex,ey,ez,cx,cy,cz,sw2d
    real(kind_real)            :: vnk,vne,vnc,snk,sne,snc
    real(kind_real)            :: rca,rcb,rcc,rca_prec,rcb_prec,rcc_prec
    real(kind_real)            :: a,kcp,rov,coe,gama,ae,cv
    real(kind_real)            :: cp,cp_prl,cp_prt,c5
    type(fld_array_t), pointer :: pv(:),sxyz(:),vol(:),vsl(:),vst(:)
    type(fld_array_t), pointer :: c(:)
    type(fld_array_t), pointer :: src(:),src_prec(:),srv(:),rdt(:)
    !type(fld_array_t), pointer :: nonprec(:)
    
    !todo Ԥ�������
    integer(kind_int)          :: oned_index,l,m,s,nsp,nep
    real(kind_real)            :: ro,p,sn(1:3,1:3),prim(1:5)
    real(kind_real)            :: nonpreMatrix(1:4)
    real(kind_real)            :: eigenValue(1:3,1:5)
    real(kind_real)            :: sumSpectralRadius,length,visl,vist,vis

    nsp = 1
    nep = 5
    gama = gamma
    ae = gama - one

    cv = refbeta/ae
    cp = gama*cv
    cp_prl = cp/prlam
    cp_prt = cp/prtur
    
    vis = 0.
    length = 1.

    do nc=1,nblkcoms
        nb  =  blkcomssp(nc)%nb
        top => blkcomssp(nc)%top

        ni = top%nijk(1)
        nj = top%nijk(2)
        nk = top%nijk(3)

        pv   => mb_pv(nb)%fld
        sxyz => mb_sxyzsp(nb)%fld
        vol  => mb_volsp(nb)%fld
        c    => mb_c(nb)%fld
        src  => mb_src(nb)%fld                !<Ac�����װ뾶
        src_prec  => mb_src_prec(nb)%fld      !<����Ԥ������Ķ����ſɱȾ����װ뾶
        rdt  => mb_rdt(nb)%fld
        !nonprec => mb_nonprec(nb)%fld

        if (nvis > nvis_euler) then
            vsl => mb_vsl(nb)%fld
            vst => mb_vst(nb)%fld
            srv => mb_srv(nb)%fld
        end if
!$OMP parallel  private(k,j,i,kx,ky,kz,ex,ey,ez,cx,cy,cz,vx,vy,vz,vnk,vne,vnc,snk,sne,snc,a, &
!$OMP                   rca,rcb,rcc,rov,visl,vist,vis,kcp,c5,coe,rva,rvb,rvc)
!$OMP do
        do k=0,nk+1
        do j=0,nj+1
        do i=0,ni+1
            
            kx = sxyz(1)%r3d(i,j,k)
            ky = sxyz(2)%r3d(i,j,k)
            kz = sxyz(3)%r3d(i,j,k)
            
            do m = 1,3
                sn(1,m) = sxyz(m)%r3d(i,j,k)
            end do

            ex = sxyz(4)%r3d(i,j,k)
            ey = sxyz(5)%r3d(i,j,k)
            ez = sxyz(6)%r3d(i,j,k)
            
            do m = 1,3
                sn(2,m) = sxyz(m+3)%r3d(i,j,k)
            end do

            cx = sxyz(7)%r3d(i,j,k)
            cy = sxyz(8)%r3d(i,j,k)
            cz = sxyz(9)%r3d(i,j,k)
            
            do m = 1,3
                sn(3,m) = sxyz(m+6)%r3d(i,j,k)
            end do

            ro = pv(1)%r3d(i,j,k)
            vx = pv(2)%r3d(i,j,k)
            vy = pv(3)%r3d(i,j,k)
            vz = pv(4)%r3d(i,j,k)
            p  = pv(5)%r3d(i,j,k)
            
            do m = 1,5
                prim(m) = pv(m)%r3d(i,j,k)
            end do
            
            if (nvis > nvis_euler) then
                visl = vsl(1)%r3d(i,j,k)
                vist = vst(1)%r3d(i,j,k)
                vis = visl + vist
                length = vol(1)%r3d(i,j,k)**(1.0/3.0)
            end if
            
            call calculate_precondition_eigenvalue_point(nsp,nep,sn,prim,vis,length,eigenvalue)
            
            !nonprec(:)%r3d(i,j,k) = nonpreMatrix(:)

            rca_prec = abs(eigenvalue(1,4))
            rcb_prec = abs(eigenvalue(2,4))
            rcc_prec = abs(eigenvalue(3,4))
            src_prec(1)%r3d(i,j,k) = rca_prec
            src_prec(2)%r3d(i,j,k) = rcb_prec
            src_prec(3)%r3d(i,j,k) = rcc_prec*fsw_kdir
            
            rdt(1)%r3d(i,j,k) = rca_prec + rcb_prec + rcc_prec*fsw_kdir
            
            !���ڼ���ʱ�䲽�����װ뾶
            if (nvis > nvis_euler) then
                rdt(1)%r3d(i,j,k) = rdt(1)%r3d(i,j,k) +srv(1)%r3d(i,j,k) + srv(2)%r3d(i,j,k) + srv(3)%r3d(i,j,k)
            end if
            
        end do
        end do
        end do
!$OMP end  do
!$OMP end parallel
    end do

end subroutine calc_spectral_radius_prec_sp    
    
subroutine calc_diag_inv_matrix
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : one,two,four3rd,nvis_euler,nlhs_prsgs_ust_sca
    use mod_datatypes, only : fld_array_t,top_block_t
    use mod_variables, only : nvis,reue,csrvis,fsw_kdir
    use mod_variables, only : gamma,prlam,prtur,refbeta
    use mod_variables, only : nlhs,dtime
    use mod_fieldvars, only : nblkcoms,blkcoms
    use mod_fieldvars, only : mb_pv,mb_dt,mb_vol
    use mod_fieldvars, only : mb_vsl,mb_vst
    !use mod_fieldvars, only : mb_prec_diag_inv
    use mod_fieldvars, only : mb_c,mb_src,mb_src_prec,mb_srv,mb_rdt
    use mod_openmp
    implicit none
    integer(kind_int)          :: nc,nb,i,j,k,ierr
    integer(kind_int)          :: ni,nj,nk
    type(top_block_t), pointer :: top
    real(kind_real)            :: rca,rcb,rcc
    real(kind_real)            :: gama
    type(fld_array_t), pointer :: pv(:)
    type(fld_array_t), pointer :: src(:),srv(:),dt(:),vol(:),vsl(:),vst(:)
    type(fld_array_t), pointer :: diagInv(:)
    
    !todo Ԥ�������
    integer(kind_int)          :: oned_index,l,m,s,nsp,nep
    real(kind_real)            :: ro,p,prim(1:5),visl,vist,vis,length
    real(kind_real)            :: diagInvMatrix(1:25)
    real(kind_real)            :: sumSpectralRadius
    real(kind_real)            :: diagElement

    nsp = 1
    nep = 5
    
    vis = 0.
    length = 1.

    do nc=1,nblkcoms
        nb  =  blkcoms(nc)%nb
        top => blkcoms(nc)%top

        ni = top%nijk(1)
        nj = top%nijk(2)
        nk = top%nijk(3)

        pv   => mb_pv(nb)%fld
        src  => mb_src(nb)%fld                !<Ac�����װ뾶
        dt   => mb_dt(nb)%fld
        !diagInv => mb_prec_diag_inv(nb)%fld
        vol  => mb_vol(nb)%fld

        if (nvis > nvis_euler) then
            vsl => mb_vsl(nb)%fld
            vst => mb_vst(nb)%fld
            srv => mb_srv(nb)%fld
        end if
!$OMP parallel  private(k,j,i,kx,ky,kz,ex,ey,ez,cx,cy,cz,vx,vy,vz,vnk,vne,vnc,snk,sne,snc,a, &
!$OMP                   rca,rcb,rcc,rov,visl,vist,vis,kcp,c5,coe,rva,rvb,rvc)
!$OMP do
        do k=0,nk+1
        do j=0,nj+1
        do i=0,ni+1
            
            do m=1,5
                prim(m) = pv(m)%r3d(i,j,k)
            end do
            
            !����Ԥ������ĶԽǾ���D
            rca = src(1)%r3d(i,j,k)
            rcb = src(2)%r3d(i,j,k)
            rcc = src(3)%r3d(i,j,k)
            
            sumSpectralRadius = rca + rcb + rcc*fsw_kdir
            
            if (nvis > nvis_euler) then
                visl = vsl(1)%r3d(i,j,k)
                vist = vst(1)%r3d(i,j,k)
                vis = visl + vist
                length = vol(1)%r3d(i,j,k)**(1.0/3.0)
                sumSpectralRadius = sumSpectralRadius +srv(1)%r3d(i,j,k) + srv(2)%r3d(i,j,k) + srv(3)%r3d(i,j,k)
            end if
            
            if (nlhs == nlhs_prsgs_ust_sca) then
                sumSpectralRadius = sumSpectralRadius + 1.5*vol(1)%r3d(i,j,k)/dtime
            end if
            
            !����D����
            diagElement = dt(1)%r3d(i,j,k)
            
            call calculateDiagMatrixInv(prim,diagElement,sumSpectralRadius,vis,length,diagInvMatrix)
            
            do m = 1,25
                diagInv(m)%r3d(i,j,k) = diagInvMatrix(m)
            end do

        end do
        end do
        end do
!$OMP end  do
!$OMP end parallel
    end do

end subroutine calc_diag_inv_matrix

subroutine calc_diag_inv_matrix_sp
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : one,two,four3rd,nvis_euler,nlhs_prsgs_ust_sca
    use mod_datatypes, only : fld_array_t,top_block_t
    use mod_variables, only : nvis,reue,csrvis,fsw_kdir
    use mod_variables, only : gamma,prlam,prtur,refbeta
    use mod_variables, only : nlhs,dtime
    use mod_fieldvars, only : nblkcoms,blkcomssp
    use mod_fieldvars, only : mb_pv,mb_dt,mb_vol
    use mod_fieldvars, only : mb_vsl,mb_vst
    !use mod_fieldvars, only : mb_prec_diag_inv
    use mod_fieldvars, only : mb_c,mb_src,mb_src_prec,mb_srv,mb_rdt
    use mod_openmp
    implicit none
    integer(kind_int)          :: nc,nb,i,j,k,ierr
    integer(kind_int)          :: ni,nj,nk
    type(top_block_t), pointer :: top
    real(kind_real)            :: rca,rcb,rcc
    real(kind_real)            :: gama
    type(fld_array_t), pointer :: pv(:)
    type(fld_array_t), pointer :: src(:),srv(:),dt(:),vol(:),vsl(:),vst(:)
    type(fld_array_t), pointer :: diagInv(:)
    
    !todo Ԥ�������
    integer(kind_int)          :: oned_index,l,m,s,nsp,nep
    real(kind_real)            :: ro,p,prim(1:5),visl,vist,vis,length
    real(kind_real)            :: diagInvMatrix(1:25)
    real(kind_real)            :: sumSpectralRadius
    real(kind_real)            :: diagElement

    nsp = 1
    nep = 5
    
    vis = 0.
    length = 1.

    do nc=1,nblkcoms
        nb  =  blkcomssp(nc)%nb
        top => blkcomssp(nc)%top

        ni = top%nijk(1)
        nj = top%nijk(2)
        nk = top%nijk(3)

        pv   => mb_pv(nb)%fld
        src  => mb_src(nb)%fld                !<Ac�����װ뾶
        dt   => mb_dt(nb)%fld
        !diagInv => mb_prec_diag_inv(nb)%fld
        vol  => mb_vol(nb)%fld

        if (nvis > nvis_euler) then
            vsl => mb_vsl(nb)%fld
            vst => mb_vst(nb)%fld
            srv => mb_srv(nb)%fld
        end if
!$OMP parallel  private(k,j,i,kx,ky,kz,ex,ey,ez,cx,cy,cz,vx,vy,vz,vnk,vne,vnc,snk,sne,snc,a, &
!$OMP                   rca,rcb,rcc,rov,visl,vist,vis,kcp,c5,coe,rva,rvb,rvc)
!$OMP do
        do k=0,nk+1
        do j=0,nj+1
        do i=0,ni+1
            
            do m=1,5
                prim(m) = pv(m)%r3d(i,j,k)
            end do
            
            !����Ԥ������ĶԽǾ���D
            rca = src(1)%r3d(i,j,k)
            rcb = src(2)%r3d(i,j,k)
            rcc = src(3)%r3d(i,j,k)
            
            sumSpectralRadius = rca + rcb + rcc*fsw_kdir
            
            if (nvis > nvis_euler) then
                visl = vsl(1)%r3d(i,j,k)
                vist = vst(1)%r3d(i,j,k)
                vis = visl + vist
                length = vol(1)%r3d(i*2,j*2,k*2)**(1.0/3.0)
                sumSpectralRadius = sumSpectralRadius +srv(1)%r3d(i,j,k) + srv(2)%r3d(i,j,k) + srv(3)%r3d(i,j,k)
            end if
            
            if (nlhs == nlhs_prsgs_ust_sca) then
                sumSpectralRadius = sumSpectralRadius + 1.5*vol(1)%r3d(i*2,j*2,k*2)/dtime
            end if
            
            !����D����
            diagElement = dt(1)%r3d(i,j,k)
            
            call calculateDiagMatrixInv(prim,diagElement,sumSpectralRadius,vis,length,diagInvMatrix)
            
            do m = 1,25
                diagInv(m)%r3d(i,j,k) = diagInvMatrix(m)
            end do

        end do
        end do
        end do
!$OMP end  do
!$OMP end parallel
    end do

end subroutine calc_diag_inv_matrix_sp

subroutine calc_spectral_radius
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : one,two,four3rd,nvis_euler
    use mod_datatypes, only : fld_array_t,top_block_t
    use mod_variables, only : nvis,reue,csrvis,fsw_kdir
    use mod_variables, only : gamma,prlam,prtur,refbeta
    use mod_fieldvars, only : nblkcoms,blkcoms
    use mod_fieldvars, only : mb_pv,mb_sxyz,mb_vol
    use mod_fieldvars, only : mb_vsl,mb_vst
    use mod_fieldvars, only : mb_c,mb_src,mb_srv,mb_rdt
    use mod_openmp
    implicit none
    integer(kind_int)          :: nc,nb,i,j,k,ierr
    integer(kind_int)          :: ni,nj,nk
    type(top_block_t), pointer :: top
    real(kind_real)            :: vx,vy,vz,kx,ky,kz
    real(kind_real)            :: ex,ey,ez,cx,cy,cz,sw2d
    real(kind_real)            :: vnk,vne,vnc,snk,sne,snc
    real(kind_real)            :: rca,rcb,rcc,rva,rvb,rvc
    real(kind_real)            :: a,vis,kcp,rov,coe,gama,ae,cv
    real(kind_real)            :: cp,cp_prl,cp_prt,visl,vist,c5
    type(fld_array_t), pointer :: pv(:),sxyz(:),vol(:)
    type(fld_array_t), pointer :: vsl(:),vst(:),c(:)
    type(fld_array_t), pointer :: src(:),srv(:),rdt(:)

    gama = gamma
    ae = gama - one

    cv = refbeta/ae
    cp = gama*cv
    cp_prl = cp/prlam
    cp_prt = cp/prtur

    do nc=1,nblkcoms
        nb  =  blkcoms(nc)%nb
        top => blkcoms(nc)%top

        ni = top%nijk(1)
        nj = top%nijk(2)
        nk = top%nijk(3)

        pv   => mb_pv(nb)%fld
        sxyz => mb_sxyz(nb)%fld
        vol  => mb_vol(nb)%fld
        c    => mb_c(nb)%fld
        src  => mb_src(nb)%fld
        rdt  => mb_rdt(nb)%fld

        if (nvis > nvis_euler) then
            vsl => mb_vsl(nb)%fld
            vst => mb_vst(nb)%fld
            srv => mb_srv(nb)%fld
        end if
!$OMP parallel  private(k,j,i,kx,ky,kz,ex,ey,ez,cx,cy,cz,vx,vy,vz,vnk,vne,vnc,snk,sne,snc,a, &
!$OMP                   rca,rcb,rcc,rov,visl,vist,vis,kcp,c5,coe,rva,rvb,rvc)
!$OMP do
        do k=0,nk+1
        do j=0,nj+1
        do i=0,ni+1
            kx = sxyz(1)%r3d(i,j,k)
            ky = sxyz(2)%r3d(i,j,k)
            kz = sxyz(3)%r3d(i,j,k)

            ex = sxyz(4)%r3d(i,j,k)
            ey = sxyz(5)%r3d(i,j,k)
            ez = sxyz(6)%r3d(i,j,k)

            cx = sxyz(7)%r3d(i,j,k)
            cy = sxyz(8)%r3d(i,j,k)
            cz = sxyz(9)%r3d(i,j,k)

            vx = pv(2)%r3d(i,j,k)
            vy = pv(3)%r3d(i,j,k)
            vz = pv(4)%r3d(i,j,k)

            vnk = kx*vx + ky*vy + kz*vz
            vne = ex*vx + ey*vy + ez*vz
            vnc = cx*vx + cy*vy + cz*vz

            snk = kx*kx + ky*ky + kz*kz
            sne = ex*ex + ey*ey + ez*ez
            snc = cx*cx + cy*cy + cz*cz

            a = c(1)%r3d(i,j,k)

            rca = abs(vnk) + a*sqrt(snk)
            rcb = abs(vne) + a*sqrt(sne)
            rcc = abs(vnc) + a*sqrt(snc)
            src(1)%r3d(i,j,k) = rca
            src(2)%r3d(i,j,k) = rcb
            src(3)%r3d(i,j,k) = rcc*fsw_kdir

            rdt(1)%r3d(i,j,k) = rca + rcb + rcc

            if (nvis > nvis_euler) then
                rov = pv(1)%r3d(i,j,k)*vol(1)%r3d(i,j,k)

                visl = vsl(1)%r3d(i,j,k)
                vist = vst(1)%r3d(i,j,k)
                vis = visl + vist
                kcp = visl*cp_prl + vist*cp_prt

                c5 = max(kcp/(cv*vis),four3rd)

                coe = csrvis*two*c5*vis/(reue*rov)

                rva = snk*coe
                rvb = sne*coe
                rvc = snc*coe
                srv(1)%r3d(i,j,k) = rva
                srv(2)%r3d(i,j,k) = rvb
                srv(3)%r3d(i,j,k) = rvc*fsw_kdir

                rdt(1)%r3d(i,j,k) = rdt(1)%r3d(i,j,k) + (rva + rvb + rvc)
            end if
        end do
        end do
        end do
!$OMP end  do
!$OMP end parallel
    end do

end subroutine calc_spectral_radius

subroutine calc_spectral_radius_sp
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : one,two,four3rd,nvis_euler
    use mod_datatypes, only : fld_array_t,top_block_t
    use mod_variables, only : nvis,reue,csrvis,fsw_kdir
    use mod_variables, only : gamma,prlam,prtur,refbeta
    use mod_fieldvars, only : nblkcoms,blkcomssp
    use mod_fieldvars, only : mb_pv,mb_sxyzsp,mb_volsp
    use mod_fieldvars, only : mb_vsl,mb_vst
    use mod_fieldvars, only : mb_c,mb_src,mb_srv,mb_rdt
    use mod_openmp
    implicit none
    integer(kind_int)          :: nc,nb,i,j,k,ierr
    integer(kind_int)          :: ni,nj,nk
    type(top_block_t), pointer :: top
    real(kind_real)            :: vx,vy,vz,kx,ky,kz
    real(kind_real)            :: ex,ey,ez,cx,cy,cz,sw2d
    real(kind_real)            :: vnk,vne,vnc,snk,sne,snc
    real(kind_real)            :: rca,rcb,rcc,rva,rvb,rvc
    real(kind_real)            :: a,vis,kcp,rov,coe,gama,ae,cv
    real(kind_real)            :: cp,cp_prl,cp_prt,visl,vist,c5
    type(fld_array_t), pointer :: pv(:),sxyz(:),vol(:)
    type(fld_array_t), pointer :: vsl(:),vst(:),c(:)
    type(fld_array_t), pointer :: src(:),srv(:),rdt(:)

    gama = gamma
    ae = gama - one

    cv = refbeta/ae
    cp = gama*cv
    cp_prl = cp/prlam
    cp_prt = cp/prtur

    do nc=1,nblkcoms
        nb  =  blkcomssp(nc)%nb
        top => blkcomssp(nc)%top

        ni = top%nijk(1)
        nj = top%nijk(2)
        nk = top%nijk(3)

        pv   => mb_pv(nb)%fld
        sxyz => mb_sxyzsp(nb)%fld
        vol  => mb_volsp(nb)%fld
        c    => mb_c(nb)%fld
        src  => mb_src(nb)%fld
        rdt  => mb_rdt(nb)%fld

        if (nvis > nvis_euler) then
            vsl => mb_vsl(nb)%fld
            vst => mb_vst(nb)%fld
            srv => mb_srv(nb)%fld
        end if
!$OMP parallel  private(k,j,i,kx,ky,kz,ex,ey,ez,cx,cy,cz,vx,vy,vz,vnk,vne,vnc,snk,sne,snc,a, &
!$OMP                   rca,rcb,rcc,rov,visl,vist,vis,kcp,c5,coe,rva,rvb,rvc)
!$OMP do
        do k=0,nk+1
        do j=0,nj+1
        do i=0,ni+1
            kx = sxyz(1)%r3d(i,j,k)
            ky = sxyz(2)%r3d(i,j,k)
            kz = sxyz(3)%r3d(i,j,k)

            ex = sxyz(4)%r3d(i,j,k)
            ey = sxyz(5)%r3d(i,j,k)
            ez = sxyz(6)%r3d(i,j,k)

            cx = sxyz(7)%r3d(i,j,k)
            cy = sxyz(8)%r3d(i,j,k)
            cz = sxyz(9)%r3d(i,j,k)

            vx = pv(2)%r3d(i,j,k)
            vy = pv(3)%r3d(i,j,k)
            vz = pv(4)%r3d(i,j,k)

            vnk = kx*vx + ky*vy + kz*vz
            vne = ex*vx + ey*vy + ez*vz
            vnc = cx*vx + cy*vy + cz*vz

            snk = kx*kx + ky*ky + kz*kz
            sne = ex*ex + ey*ey + ez*ez
            snc = cx*cx + cy*cy + cz*cz

            a = c(1)%r3d(i,j,k)

            rca = abs(vnk) + a*sqrt(snk)
            rcb = abs(vne) + a*sqrt(sne)
            rcc = abs(vnc) + a*sqrt(snc)
            src(1)%r3d(i,j,k) = rca
            src(2)%r3d(i,j,k) = rcb
            src(3)%r3d(i,j,k) = rcc*fsw_kdir

            rdt(1)%r3d(i,j,k) = rca + rcb + rcc

            if (nvis > nvis_euler) then
                rov = pv(1)%r3d(i,j,k)*vol(1)%r3d(i,j,k)

                visl = vsl(1)%r3d(i,j,k)
                vist = vst(1)%r3d(i,j,k)
                vis = visl + vist
                kcp = visl*cp_prl + vist*cp_prt

                c5 = max(kcp/(cv*vis),four3rd)

                coe = csrvis*two*c5*vis/(reue*rov)

                rva = snk*coe
                rvb = sne*coe
                rvc = snc*coe
                srv(1)%r3d(i,j,k) = rva
                srv(2)%r3d(i,j,k) = rvb
                srv(3)%r3d(i,j,k) = rvc*fsw_kdir

                rdt(1)%r3d(i,j,k) = rdt(1)%r3d(i,j,k) + (rva + rvb + rvc)
            end if
        end do
        end do
        end do
!$OMP end  do
!$OMP end parallel
    end do

end subroutine calc_spectral_radius_sp

!< todo SCM-P �����װ뾶
subroutine calc_spectral_radius_scmp
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : one,two,four3rd,nvis_euler
    use mod_datatypes, only : fld_array_t,top_block_t
    use mod_variables, only : nvis,reue,csrvis,fsw_kdir
    use mod_variables, only : gamma,prlam,prtur,refbeta
    use mod_variables, only : moo   !> Done
    use mod_fieldvars, only : nblkcoms,blkcoms
    use mod_fieldvars, only : mb_pv,mb_sxyz,mb_vol
    use mod_fieldvars, only : mb_vsl,mb_vst
    use mod_fieldvars, only : mb_c,mb_src,mb_srv,mb_rdt
    use mod_openmp
    implicit none
    integer(kind_int)          :: nc,nb,i,j,k,ierr
    integer(kind_int)          :: ni,nj,nk
    type(top_block_t), pointer :: top
    real(kind_real)            :: vx,vy,vz,kx,ky,kz
    real(kind_real)            :: ex,ey,ez,cx,cy,cz,sw2d
    real(kind_real)            :: vnk,vne,vnc,snk,sne,snc
    real(kind_real)            :: rca,rcb,rcc,rva,rvb,rvc
    real(kind_real)            :: a,vis,kcp,rov,coe,gama,ae,cv
    real(kind_real)            :: cp,cp_prl,cp_prt,visl,vist,c5
    type(fld_array_t), pointer :: pv(:),sxyz(:),vol(:)
    type(fld_array_t), pointer :: vsl(:),vst(:),c(:)
    type(fld_array_t), pointer :: src(:),srv(:),rdt(:)

    real(kind_real)            :: theta, a2        !> Done
    real(kind_real)            :: dk2, Delta_k     !> Done
    real(kind_real)            :: de2, Delta_e     !> Done
    real(kind_real)            :: dc2, Delta_c     !> Done
    
    gama = gamma
    ae = gama - one

    cv = refbeta/ae
    cp = gama*cv
    cp_prl = cp/prlam
    cp_prt = cp/prtur

    do nc=1,nblkcoms
        nb  =  blkcoms(nc)%nb
        top => blkcoms(nc)%top

        ni = top%nijk(1)
        nj = top%nijk(2)
        nk = top%nijk(3)

        pv   => mb_pv(nb)%fld
        sxyz => mb_sxyz(nb)%fld
        vol  => mb_vol(nb)%fld
        c    => mb_c(nb)%fld
        src  => mb_src(nb)%fld
        rdt  => mb_rdt(nb)%fld

        if (nvis > nvis_euler) then
            vsl => mb_vsl(nb)%fld
            vst => mb_vst(nb)%fld
            srv => mb_srv(nb)%fld
        end if
!$OMP parallel  private(k,j,i,kx,ky,kz,ex,ey,ez,cx,cy,cz,vx,vy,vz,vnk,vne,vnc,snk,sne,snc,a, &
!$OMP                   rca,rcb,rcc,rov,visl,vist,vis,kcp,c5,coe,rva,rvb,rvc)
!$OMP do
        do k=0,nk+1
        do j=0,nj+1
        do i=0,ni+1
            kx = sxyz(1)%r3d(i,j,k)
            ky = sxyz(2)%r3d(i,j,k)
            kz = sxyz(3)%r3d(i,j,k)

            ex = sxyz(4)%r3d(i,j,k)
            ey = sxyz(5)%r3d(i,j,k)
            ez = sxyz(6)%r3d(i,j,k)

            cx = sxyz(7)%r3d(i,j,k)
            cy = sxyz(8)%r3d(i,j,k)
            cz = sxyz(9)%r3d(i,j,k)

            vx = pv(2)%r3d(i,j,k)
            vy = pv(3)%r3d(i,j,k)
            vz = pv(4)%r3d(i,j,k)

            vnk = kx*vx + ky*vy + kz*vz
            vne = ex*vx + ey*vy + ez*vz
            vnc = cx*vx + cy*vy + cz*vz

            snk = kx*kx + ky*ky + kz*kz
            sne = ex*ex + ey*ey + ez*ez
            snc = cx*cx + cy*cy + cz*cz

            a = c(1)%r3d(i,j,k)
            
            dk2 = kx*kx + ky*ky + kz*kz     !> Done
            de2 = ex*ex + ey*ey + ez*ez     !> Done
            dc2 = cx*cx + cy*cy + cz*cz     !> Done
            theta = 0.0                     !> Done
            a2    = a*a                     !> Done
            Delta_k = sqrt( vnk*vnk*(1 - 1.0/a2) + dk2 )     !> Done
            Delta_e = sqrt( vne*vne*(1 - 1.0/a2) + de2 )     !> Done
            Delta_c = sqrt( vnc*vnc*(1 - 1.0/a2) + dc2 )     !> Done

            !rca = abs(vnk) + a*sqrt(snk)
            !rcb = abs(vne) + a*sqrt(sne)
            !rcc = abs(vnc) + a*sqrt(snc)
            rca = abs(vnk) + Delta_k     !> Done
            rcb = abs(vne) + Delta_e     !> Done
            rcc = abs(vnc) + Delta_c     !> Done
            src(1)%r3d(i,j,k) = rca
            src(2)%r3d(i,j,k) = rcb
            src(3)%r3d(i,j,k) = rcc*fsw_kdir

            rdt(1)%r3d(i,j,k) = rca + rcb + rcc

            if (nvis > nvis_euler) then
                rov = pv(1)%r3d(i,j,k)*vol(1)%r3d(i,j,k)
                coe = csrvis*two/(reue*rov)

                rva = snk*coe
                rvb = sne*coe
                rvc = snc*coe
                srv(1)%r3d(i,j,k) = rva
                srv(2)%r3d(i,j,k) = rvb
                srv(3)%r3d(i,j,k) = rvc*fsw_kdir

                rdt(1)%r3d(i,j,k) = rdt(1)%r3d(i,j,k) + (rva + rvb + rvc)
            end if
        end do
        end do
        end do
!$OMP end  do
!$OMP end parallel
    end do

end subroutine calc_spectral_radius_scmp

!< todo SCM-P �����װ뾶
subroutine calc_spectral_radius_scmp_sp
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : one,two,four3rd,nvis_euler
    use mod_datatypes, only : fld_array_t,top_block_t
    use mod_variables, only : nvis,reue,csrvis,fsw_kdir
    use mod_variables, only : gamma,prlam,prtur,refbeta
    use mod_variables, only : moo   !> Done
    use mod_fieldvars, only : nblkcoms,blkcomssp
    use mod_fieldvars, only : mb_pv,mb_sxyzsp,mb_volsp
    use mod_fieldvars, only : mb_vsl,mb_vst
    use mod_fieldvars, only : mb_c,mb_src,mb_srv,mb_rdt
    use mod_openmp
    implicit none
    integer(kind_int)          :: nc,nb,i,j,k,ierr
    integer(kind_int)          :: ni,nj,nk
    type(top_block_t), pointer :: top
    real(kind_real)            :: vx,vy,vz,kx,ky,kz
    real(kind_real)            :: ex,ey,ez,cx,cy,cz,sw2d
    real(kind_real)            :: vnk,vne,vnc,snk,sne,snc
    real(kind_real)            :: rca,rcb,rcc,rva,rvb,rvc
    real(kind_real)            :: a,vis,kcp,rov,coe,gama,ae,cv
    real(kind_real)            :: cp,cp_prl,cp_prt,visl,vist,c5
    type(fld_array_t), pointer :: pv(:),sxyz(:),vol(:)
    type(fld_array_t), pointer :: vsl(:),vst(:),c(:)
    type(fld_array_t), pointer :: src(:),srv(:),rdt(:)

    real(kind_real)            :: theta, a2        !> Done
    real(kind_real)            :: dk2, Delta_k     !> Done
    real(kind_real)            :: de2, Delta_e     !> Done
    real(kind_real)            :: dc2, Delta_c     !> Done
    
    gama = gamma
    ae = gama - one

    cv = refbeta/ae
    cp = gama*cv
    cp_prl = cp/prlam
    cp_prt = cp/prtur

    do nc=1,nblkcoms
        nb  =  blkcomssp(nc)%nb
        top => blkcomssp(nc)%top

        ni = top%nijk(1)
        nj = top%nijk(2)
        nk = top%nijk(3)

        pv   => mb_pv(nb)%fld
        sxyz => mb_sxyzsp(nb)%fld
        vol  => mb_volsp(nb)%fld
        c    => mb_c(nb)%fld
        src  => mb_src(nb)%fld
        rdt  => mb_rdt(nb)%fld

        if (nvis > nvis_euler) then
            vsl => mb_vsl(nb)%fld
            vst => mb_vst(nb)%fld
            srv => mb_srv(nb)%fld
        end if
!$OMP parallel  private(k,j,i,kx,ky,kz,ex,ey,ez,cx,cy,cz,vx,vy,vz,vnk,vne,vnc,snk,sne,snc,a, &
!$OMP                   rca,rcb,rcc,rov,visl,vist,vis,kcp,c5,coe,rva,rvb,rvc)
!$OMP do
        do k=0,nk+1
        do j=0,nj+1
        do i=0,ni+1
            kx = sxyz(1)%r3d(i,j,k)
            ky = sxyz(2)%r3d(i,j,k)
            kz = sxyz(3)%r3d(i,j,k)

            ex = sxyz(4)%r3d(i,j,k)
            ey = sxyz(5)%r3d(i,j,k)
            ez = sxyz(6)%r3d(i,j,k)

            cx = sxyz(7)%r3d(i,j,k)
            cy = sxyz(8)%r3d(i,j,k)
            cz = sxyz(9)%r3d(i,j,k)

            vx = pv(2)%r3d(i,j,k)
            vy = pv(3)%r3d(i,j,k)
            vz = pv(4)%r3d(i,j,k)

            vnk = kx*vx + ky*vy + kz*vz
            vne = ex*vx + ey*vy + ez*vz
            vnc = cx*vx + cy*vy + cz*vz

            snk = kx*kx + ky*ky + kz*kz
            sne = ex*ex + ey*ey + ez*ez
            snc = cx*cx + cy*cy + cz*cz

            a = c(1)%r3d(i,j,k)
            
            dk2 = kx*kx + ky*ky + kz*kz     !> Done
            de2 = ex*ex + ey*ey + ez*ez     !> Done
            dc2 = cx*cx + cy*cy + cz*cz     !> Done
            theta = 0.0                     !> Done
            a2    = a*a                     !> Done
            Delta_k = sqrt( vnk*vnk*(1 - 1.0/a2) + dk2 )     !> Done
            Delta_e = sqrt( vne*vne*(1 - 1.0/a2) + de2 )     !> Done
            Delta_c = sqrt( vnc*vnc*(1 - 1.0/a2) + dc2 )     !> Done

            !rca = abs(vnk) + a*sqrt(snk)
            !rcb = abs(vne) + a*sqrt(sne)
            !rcc = abs(vnc) + a*sqrt(snc)
            rca = abs(vnk) + Delta_k     !> Done
            rcb = abs(vne) + Delta_e     !> Done
            rcc = abs(vnc) + Delta_c     !> Done
            src(1)%r3d(i,j,k) = rca
            src(2)%r3d(i,j,k) = rcb
            src(3)%r3d(i,j,k) = rcc*fsw_kdir

            rdt(1)%r3d(i,j,k) = rca + rcb + rcc

            if (nvis > nvis_euler) then
                rov = pv(1)%r3d(i,j,k)*vol(1)%r3d(i,j,k)
                coe = csrvis*two/(reue*rov)

                rva = snk*coe
                rvb = sne*coe
                rvc = snc*coe
                srv(1)%r3d(i,j,k) = rva
                srv(2)%r3d(i,j,k) = rvb
                srv(3)%r3d(i,j,k) = rvc*fsw_kdir

                rdt(1)%r3d(i,j,k) = rdt(1)%r3d(i,j,k) + (rva + rvb + rvc)
            end if
        end do
        end do
        end do
!$OMP end  do
!$OMP end parallel
    end do

end subroutine calc_spectral_radius_scmp_sp

subroutine calc_time_step
    use mod_constants, only : nunst_tst_local,nunst_tst_global
    use mod_constants, only : nunst_tst_localw,nunst_tst_localm
    use mod_constants, only : nunst_tst_globald,nunst_tst_localq
    use mod_variables, only : nunst
    implicit none

    select case(nunst)
    case(nunst_tst_local)
        call calc_timestep_local
    case(nunst_tst_global)
        call calc_timestep_global
    case(nunst_tst_localw)
        call calc_timestep_localw
    case(nunst_tst_localm)
        call calc_timestep_localm
    case(nunst_tst_globald)
        call calc_timestep_globald
    case(nunst_tst_localq)
        call calc_timestep_localq
    case default

    end select

end subroutine calc_time_step

subroutine calc_time_step_sp
    use mod_constants, only : nunst_tst_local,nunst_tst_global
    use mod_constants, only : nunst_tst_localw,nunst_tst_localm
    use mod_constants, only : nunst_tst_globald,nunst_tst_localq
    use mod_variables, only : nunst
    implicit none

    select case(nunst)
    case(nunst_tst_local)
        call calc_timestep_local_sp
    case(nunst_tst_global)
        call calc_timestep_global_sp
    case(nunst_tst_localw)
        call calc_timestep_localw_sp
    case(nunst_tst_localm)
        call calc_timestep_localm_sp
    case(nunst_tst_globald)
        call calc_timestep_globald_sp
    case(nunst_tst_localq)
        call calc_timestep_localq_sp
    case default

    end select

end subroutine calc_time_step_sp

subroutine calc_timestep_localq
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : one,small,large
    use mod_datatypes, only : fld_array_t,top_block_t
    use mod_variables, only : nvis,fsw_kdir
    use mod_variables, only : cfl,cflmax,cflfac
    use mod_variables, only : relaxs,relaxp,cdtmax,dtaumax
    use mod_variables, only : niter,nsub,nsubp,nsubp1,nlhs    
    use mod_fieldvars, only : nblkcoms,blkcoms
    use mod_fieldvars, only : mb_vol,mb_rdt,mb_dt,mb_rhs
    use mod_openmp
    implicit none
    integer(kind_int)          :: nc,nb,i,j,k,ierr
    integer(kind_int)          :: ni,nj,nk
    type(top_block_t), pointer :: top
    real(kind_real)            :: srad,vol0,dt0,rhs0,rhsmin,rhsmax
    real(kind_real)            :: fac,fac1,cflfacr,relax,geom,geomf
    type(fld_array_t), pointer :: vol(:),rdt(:),dt(:)
    type(fld_array_t), pointer :: src(:),srv(:),c(:),rhs(:)

    if(relaxs>0) then
        geomf   = (nsub + nsubp + nsubp1)/3.0
        relax   = geomf*relaxs
        geom    = (cflfac-abs(cdtmax))/abs(relax)
        
        if(niter == 1) then
            cflfacr = abs(cdtmax)
        else
            cflfacr = max(abs(cdtmax)+geom*(niter*1.0),cflfac)
        end if
    else
        cflfacr = cflfac
    end if
    
    do nc=1,nblkcoms
        nb  =  blkcoms(nc)%nb
        top => blkcoms(nc)%top

        ni = top%nijk(1)
        nj = top%nijk(2)
        nk = top%nijk(3)

        vol => mb_vol(nb)%fld
        rdt => mb_rdt(nb)%fld
        dt  => mb_dt(nb)%fld
        rhs => mb_rhs(nb)%fld
               
        rhsmin  = large
        rhsmax  = small
        
        do k=1,nk
        do j=1,nj
        do i=1,ni
            rhsmin = min((abs(rhs(1)%r3d(i,j,k))+abs(rhs(2)%r3d(i,j,k))+ &
                          abs(rhs(3)%r3d(i,j,k))+abs(rhs(4)%r3d(i,j,k))+ &
                          abs(rhs(5)%r3d(i,j,k))),rhsmin)
            rhsmax = max((abs(rhs(1)%r3d(i,j,k))+abs(rhs(2)%r3d(i,j,k))+ &
                          abs(rhs(3)%r3d(i,j,k))+abs(rhs(4)%r3d(i,j,k))+ &
                          abs(rhs(5)%r3d(i,j,k))),rhsmax)
        end do
        end do
        end do
        
        if(rhsmin/rhsmax<cflfac**(1.0/relaxp)) then
            rhsmin=rhsmax*cflfac**(1.0/relaxp)
        end if      

        do k=1,nk
        do j=1,nj
        do i=1,ni
            vol0 = vol(1)%r3d(i,j,k)
            rhs0 = (abs(rhs(1)%r3d(i,j,k))+abs(rhs(2)%r3d(i,j,k))+ &
                    abs(rhs(3)%r3d(i,j,k))+abs(rhs(4)%r3d(i,j,k))+ &
                    abs(rhs(5)%r3d(i,j,k)))

            srad = rdt(1)%r3d(i,j,k)
            
            fac1 = 1.0
            
            if (rhs0 > rhsmin) then
                fac1 = (rhsmin/rhs0)**relaxp
            end if
            
            fac  = max(fac1,cflfacr)

            dt0  = cfl*fac*vol0/srad
            
            dt0 = dt0/vol0
            dt(1)%r3d(i,j,k) = dt0
            rdt(1)%r3d(i,j,k) = rdt(1)%r3d(i,j,k) + one/dt0

        end do
        end do
        end do
        
    end do
    
    cflmax  = cfl
    dtaumax = dt0
    

end subroutine calc_timestep_localq

subroutine calc_timestep_localq_sp
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : one,small,large
    use mod_datatypes, only : fld_array_t,top_block_t
    use mod_variables, only : nvis,fsw_kdir
    use mod_variables, only : cfl,cflmax,cflfac
    use mod_variables, only : relaxs,relaxp,cdtmax,dtaumax
    use mod_variables, only : niter,nsub,nsubp,nsubp1,nlhs    
    use mod_fieldvars, only : nblkcoms,blkcomssp
    use mod_fieldvars, only : mb_vol,mb_rdt,mb_dt,mb_rhs
    use mod_openmp
    implicit none
    integer(kind_int)          :: nc,nb,i,j,k,ierr
    integer(kind_int)          :: ni,nj,nk
    type(top_block_t), pointer :: top
    real(kind_real)            :: srad,vol0,dt0,rhs0,rhsmin,rhsmax
    real(kind_real)            :: fac,fac1,cflfacr,relax,geom,geomf
    type(fld_array_t), pointer :: vol(:),rdt(:),dt(:)
    type(fld_array_t), pointer :: src(:),srv(:),c(:),rhs(:)

    if(relaxs>0) then
        geomf   = (nsub + nsubp + nsubp1)/3.0
        relax   = geomf*relaxs
        geom    = (cflfac-abs(cdtmax))/abs(relax)
        
        if(niter == 1) then
            cflfacr = abs(cdtmax)
        else
            cflfacr = max(abs(cdtmax)+geom*(niter*1.0),cflfac)
        end if
    else
        cflfacr = cflfac
    end if
    
    do nc=1,nblkcoms
        nb  =  blkcomssp(nc)%nb
        top => blkcomssp(nc)%top

        ni = top%nijk(1)
        nj = top%nijk(2)
        nk = top%nijk(3)

        vol => mb_vol(nb)%fld
        rdt => mb_rdt(nb)%fld
        dt  => mb_dt(nb)%fld
        rhs => mb_rhs(nb)%fld
               
        rhsmin  = large
        rhsmax  = small
        
        do k=1,nk
        do j=1,nj
        do i=1,ni
            rhsmin = min((abs(rhs(1)%r3d(i,j,k))+abs(rhs(2)%r3d(i,j,k))+ &
                          abs(rhs(3)%r3d(i,j,k))+abs(rhs(4)%r3d(i,j,k))+ &
                          abs(rhs(5)%r3d(i,j,k))),rhsmin)
            rhsmax = max((abs(rhs(1)%r3d(i,j,k))+abs(rhs(2)%r3d(i,j,k))+ &
                          abs(rhs(3)%r3d(i,j,k))+abs(rhs(4)%r3d(i,j,k))+ &
                          abs(rhs(5)%r3d(i,j,k))),rhsmax)
        end do
        end do
        end do
        
        if(rhsmin/rhsmax<cflfac**(1.0/relaxp)) then
            rhsmin=rhsmax*cflfac**(1.0/relaxp)
        end if      

        do k=1,nk
        do j=1,nj
        do i=1,ni
            vol0 = vol(1)%r3d(i*2,j*2,k*2)
            rhs0 = (abs(rhs(1)%r3d(i,j,k))+abs(rhs(2)%r3d(i,j,k))+ &
                    abs(rhs(3)%r3d(i,j,k))+abs(rhs(4)%r3d(i,j,k))+ &
                    abs(rhs(5)%r3d(i,j,k)))

            srad = rdt(1)%r3d(i,j,k)
            
            fac1 = 1.0
            
            if (rhs0 > rhsmin) then
                fac1 = (rhsmin/rhs0)**relaxp
            end if
            
            fac  = max(fac1,cflfacr)

            dt0  = cfl*fac*vol0/srad
            
            dt0 = dt0/vol0
            dt(1)%r3d(i,j,k) = dt0
            rdt(1)%r3d(i,j,k) = rdt(1)%r3d(i,j,k) + one/dt0

        end do
        end do
        end do
        
    end do
    
    cflmax  = cfl
    dtaumax = dt0
    

end subroutine calc_timestep_localq_sp

subroutine calc_timestep_local
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,one,small,large
    use mod_constants, only : nvis_euler
    use mod_datatypes, only : fld_array_t,top_block_t
    use mod_variables, only : nvis,fsw_kdir
    use mod_variables, only : cfl,cflmin,cflmax
    use mod_variables, only : cdtmax,dtaumin,dtaumax
    use mod_fieldvars, only : nblkcoms,blkcoms
    use mod_fieldvars, only : mb_vol,mb_rdt,mb_dt
    use mod_openmp
    implicit none
    integer(kind_int)          :: nc,nb,i,j,k,ierr
    integer(kind_int)          :: ni,nj,nk
    type(top_block_t), pointer :: top
    real(kind_real)            :: srad,vol0,dt0,cfl0
    type(fld_array_t), pointer :: vol(:),rdt(:),dt(:)
    type(fld_array_t), pointer :: src(:),srv(:)
#ifdef OMP_IMP
    integer(kind_int)          :: idt,numt
    real(kind_real), dimension(0:omp_max_num_threads) :: dmin, dmax, cmin, cmax
#endif

#ifdef OMP_IMP
    dmin(:) = large
    dmax(:) = small
#endif

    dtaumin = large
    dtaumax = small
    do nc=1,nblkcoms
        nb  =  blkcoms(nc)%nb
        top => blkcoms(nc)%top

        ni = top%nijk(1)
        nj = top%nijk(2)
        nk = top%nijk(3)

        vol => mb_vol(nb)%fld
        rdt => mb_rdt(nb)%fld
        dt  => mb_dt(nb)%fld

#ifdef OMP_IMP
!$OMP parallel private(vol0,srad,dt0,idt)
!$OMP master
        numt = omp_get_num_threads()
!$OMP end master
        idt  = omp_get_thread_num()
!$OMP do
#endif

        do k=1,nk
        do j=1,nj
        do i=1,ni
            vol0 = vol(1)%r3d(i,j,k)

            srad = rdt(1)%r3d(i,j,k)

            dt0  = cfl*vol0/srad

#ifdef OMP_IMP
            dmin(idt) = min(dmin(idt),dt0)
            dmax(idt) = max(dmax(idt),dt0)
#else
            dtaumin = min(dtaumin,dt0)
            dtaumax = max(dtaumax,dt0)
#endif
        end do
        end do
        end do
#ifdef OMP_IMP
!$OMP end do
!$OMP end parallel

        do i=0,numt
          dtaumin=min(dmin(i),dtaumin)
          dtaumax=max(dmax(i),dtaumax)
        end do
#endif
    end do

    call broadcast_timestep

    if (cdtmax > zero) then
        dtaumax = min(dtaumax, cdtmax*dtaumin)
    end if

    cflmin = large
    cflmax = small
#ifdef OMP_IMP
    cmin(:)=large
    cmax(:)=small
#endif
    do nc=1,nblkcoms
        nb  =  blkcoms(nc)%nb
        top => blkcoms(nc)%top

        ni = top%nijk(1)
        nj = top%nijk(2)
        nk = top%nijk(3)

        vol => mb_vol(nb)%fld
        rdt => mb_rdt(nb)%fld
        dt  => mb_dt(nb)%fld
#ifdef OMP_IMP
!$OMP parallel  private(vol0,srad,dt0,idt)
!$OMP master
        numt = omp_get_num_threads()
!$OMP end master
        idt  = omp_get_thread_num()
!$OMP do
#endif

        do k=1,nk
        do j=1,nj
        do i=1,ni
            vol0 = vol(1)%r3d(i,j,k)

            srad = rdt(1)%r3d(i,j,k)

            dt0  = cfl*vol0/srad

            dt0 = min(dt0,dtaumax)
            dt0 = dt0/vol0
            dt(1)%r3d(i,j,k) = dt0
            rdt(1)%r3d(i,j,k) = rdt(1)%r3d(i,j,k) + one/dt0

            cfl0 = dt0*srad
#ifdef OMP_IMP
            cmin(idt) = min(cmin(idt),cfl0)
            cmax(idt) = max(cmax(idt),cfl0)
#else
            cflmin = min(cflmin,cfl0)
            cflmax = max(cflmax,cfl0)
#endif
        end do
        end do
        end do
#ifdef OMP_IMP
!$OMP end do
!$OMP end parallel
        do i=0,numt
          cflmin=min(cmin(i),cflmin)
          cflmax=max(cmax(i),cflmax)
        end do
#endif
    end do

end subroutine calc_timestep_local

subroutine calc_timestep_local_sp
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,one,small,large
    use mod_constants, only : nvis_euler
    use mod_datatypes, only : fld_array_t,top_block_t
    use mod_variables, only : nvis,fsw_kdir
    use mod_variables, only : cfl,cflmin,cflmax
    use mod_variables, only : cdtmax,dtaumin,dtaumax
    use mod_fieldvars, only : nblkcoms,blkcomssp
    use mod_fieldvars, only : mb_volsp,mb_rdt,mb_dt
    use mod_openmp
    implicit none
    integer(kind_int)          :: nc,nb,i,j,k,ierr
    integer(kind_int)          :: ni,nj,nk
    type(top_block_t), pointer :: top
    real(kind_real)            :: srad,vol0,dt0,cfl0
    type(fld_array_t), pointer :: vol(:),rdt(:),dt(:)
    type(fld_array_t), pointer :: src(:),srv(:)
#ifdef OMP_IMP
    integer(kind_int)          :: idt,numt
    real(kind_real), dimension(0:omp_max_num_threads) :: dmin, dmax, cmin, cmax
#endif

#ifdef OMP_IMP
    dmin(:) = large
    dmax(:) = small
#endif

    dtaumin = large
    dtaumax = small
    do nc=1,nblkcoms
        nb  =  blkcomssp(nc)%nb
        top => blkcomssp(nc)%top

        ni = top%nijk(1)
        nj = top%nijk(2)
        nk = top%nijk(3)

        vol => mb_volsp(nb)%fld
        rdt => mb_rdt(nb)%fld
        dt  => mb_dt(nb)%fld

#ifdef OMP_IMP
!$OMP parallel private(vol0,srad,dt0,idt)
!$OMP master
        numt = omp_get_num_threads()
!$OMP end master
        idt  = omp_get_thread_num()
!$OMP do
#endif

        do k=1,nk
        do j=1,nj
        do i=1,ni
            vol0 = vol(1)%r3d(i,j,k)

            srad = rdt(1)%r3d(i,j,k)

            dt0  = cfl*vol0/srad

#ifdef OMP_IMP
            dmin(idt) = min(dmin(idt),dt0)
            dmax(idt) = max(dmax(idt),dt0)
#else
            dtaumin = min(dtaumin,dt0)
            dtaumax = max(dtaumax,dt0)
#endif
        end do
        end do
        end do
#ifdef OMP_IMP
!$OMP end do
!$OMP end parallel

        do i=0,numt
          dtaumin=min(dmin(i),dtaumin)
          dtaumax=max(dmax(i),dtaumax)
        end do
#endif
    end do

    call broadcast_timestep

    if (cdtmax > zero) then
        dtaumax = min(dtaumax, cdtmax*dtaumin)
    end if

    cflmin = large
    cflmax = small
#ifdef OMP_IMP
    cmin(:)=large
    cmax(:)=small
#endif
    do nc=1,nblkcoms
        nb  =  blkcomssp(nc)%nb
        top => blkcomssp(nc)%top

        ni = top%nijk(1)
        nj = top%nijk(2)
        nk = top%nijk(3)

        vol => mb_volsp(nb)%fld
        rdt => mb_rdt(nb)%fld
        dt  => mb_dt(nb)%fld
#ifdef OMP_IMP
!$OMP parallel  private(vol0,srad,dt0,idt)
!$OMP master
        numt = omp_get_num_threads()
!$OMP end master
        idt  = omp_get_thread_num()
!$OMP do
#endif

        do k=1,nk
        do j=1,nj
        do i=1,ni
            vol0 = vol(1)%r3d(i,j,k)

            srad = rdt(1)%r3d(i,j,k)

            dt0  = cfl*vol0/srad

            dt0 = min(dt0,dtaumax)
            dt0 = dt0/vol0
            dt(1)%r3d(i,j,k) = dt0
            rdt(1)%r3d(i,j,k) = rdt(1)%r3d(i,j,k) + one/dt0

            cfl0 = dt0*srad
#ifdef OMP_IMP
            cmin(idt) = min(cmin(idt),cfl0)
            cmax(idt) = max(cmax(idt),cfl0)
#else
            cflmin = min(cflmin,cfl0)
            cflmax = max(cflmax,cfl0)
#endif
        end do
        end do
        end do
#ifdef OMP_IMP
!$OMP end do
!$OMP end parallel
        do i=0,numt
          cflmin=min(cmin(i),cflmin)
          cflmax=max(cmax(i),cflmax)
        end do
#endif
    end do

end subroutine calc_timestep_local_sp

subroutine calc_timestep_global
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : small,large,one
    use mod_datatypes, only : fld_array_t,top_block_t
    use mod_variables, only : cflmin,cflmax,dtau
    use mod_variables, only : dtaumin,dtaumax
    use mod_fieldvars, only : nblkcoms,blkcoms
    use mod_fieldvars, only : mb_vol,mb_rdt,mb_dt
    implicit none
    integer(kind_int)          :: nc,nb,i,j,k,ierr
    integer(kind_int)          :: ni,nj,nk
    type(top_block_t), pointer :: top
    real(kind_real)            :: srad,vol0,dt0,cfl0
    type(fld_array_t), pointer :: vol(:),rdt(:),dt(:)

    cflmin = large
    cflmax = small
    do nc=1,nblkcoms
        nb  =  blkcoms(nc)%nb
        top => blkcoms(nc)%top

        ni = top%nijk(1)
        nj = top%nijk(2)
        nk = top%nijk(3)

        vol => mb_vol(nb)%fld
        rdt => mb_rdt(nb)%fld
        dt  => mb_dt(nb)%fld

        do k=1,nk
        do j=1,nj
        do i=1,ni
            srad = rdt(1)%r3d(i,j,k)
            vol0 = vol(1)%r3d(i,j,k)
            dt0  = dtau/vol0
            dt(1)%r3d(i,j,k) = dt0
            rdt(1)%r3d(i,j,k) = rdt(1)%r3d(i,j,k) + one/dt0

            cfl0 = dt0*srad
            cflmin = min(cflmin,cfl0)
            cflmax = max(cflmax,cfl0)
        end do
        end do
        end do
    end do
    dtaumin = dtau
    dtaumax = dtau

end subroutine calc_timestep_global

subroutine calc_timestep_global_sp
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : small,large,one
    use mod_datatypes, only : fld_array_t,top_block_t
    use mod_variables, only : cflmin,cflmax,dtau
    use mod_variables, only : dtaumin,dtaumax
    use mod_fieldvars, only : nblkcoms,blkcomssp
    use mod_fieldvars, only : mb_volsp,mb_rdt,mb_dt
    implicit none
    integer(kind_int)          :: nc,nb,i,j,k,ierr
    integer(kind_int)          :: ni,nj,nk
    type(top_block_t), pointer :: top
    real(kind_real)            :: srad,vol0,dt0,cfl0
    type(fld_array_t), pointer :: vol(:),rdt(:),dt(:)

    cflmin = large
    cflmax = small
    do nc=1,nblkcoms
        nb  =  blkcomssp(nc)%nb
        top => blkcomssp(nc)%top

        ni = top%nijk(1)
        nj = top%nijk(2)
        nk = top%nijk(3)

        vol => mb_volsp(nb)%fld
        rdt => mb_rdt(nb)%fld
        dt  => mb_dt(nb)%fld

        do k=1,nk
        do j=1,nj
        do i=1,ni
            srad = rdt(1)%r3d(i,j,k)
            vol0 = vol(1)%r3d(i,j,k)
            dt0  = dtau/vol0
            dt(1)%r3d(i,j,k) = dt0
            rdt(1)%r3d(i,j,k) = rdt(1)%r3d(i,j,k) + one/dt0

            cfl0 = dt0*srad
            cflmin = min(cflmin,cfl0)
            cflmax = max(cflmax,cfl0)
        end do
        end do
        end do
    end do
    dtaumin = dtau
    dtaumax = dtau

end subroutine calc_timestep_global_sp

subroutine calc_timestep_localw
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,one,small,large
    use mod_constants, only : nvis_euler
    use mod_datatypes, only : fld_array_t,top_block_t
    use mod_variables, only : nvis,fsw_kdir
    use mod_variables, only : cfl,cflmin,cflmax
    use mod_variables, only : cdtmax,dtaumin,dtaumax
    use mod_fieldvars, only : nblkcoms,blkcoms
    use mod_fieldvars, only : mb_vol,mb_rdt,mb_dt
    use mod_fieldvars, only : mb_src,mb_srv
    implicit none
    integer(kind_int)          :: nc,nb,i,j,k,ierr
    integer(kind_int)          :: ni,nj,nk
    type(top_block_t), pointer :: top
    real(kind_real)            :: ra,rb,rc,osr
    real(kind_real)            :: srad,vol0,dt0,cfl0
    type(fld_array_t), pointer :: vol(:),rdt(:),dt(:)
    type(fld_array_t), pointer :: src(:),srv(:)

    dtaumin = large
    dtaumax = small
    do nc=1,nblkcoms
        nb  =  blkcoms(nc)%nb
        top => blkcoms(nc)%top

        ni = top%nijk(1)
        nj = top%nijk(2)
        nk = top%nijk(3)

        vol => mb_vol(nb)%fld
        rdt => mb_rdt(nb)%fld
        dt  => mb_dt(nb)%fld
        src => mb_src(nb)%fld
        if (nvis > nvis_euler) then
            srv => mb_srv(nb)%fld
        end if

        do k=1,nk
        do j=1,nj
        do i=1,ni
            vol0 = vol(1)%r3d(i,j,k)

            ra = src(1)%r3d(i,j,k)
            rb = src(2)%r3d(i,j,k)
            rc = src(3)%r3d(i,j,k) + small
            if (nvis > nvis_euler) then
                ra = ra + srv(1)%r3d(i,j,k)
                rb = rb + srv(2)%r3d(i,j,k)
                rc = rc + srv(3)%r3d(i,j,k)
            end if
            osr = one/ra + one/rb + fsw_kdir/rc
            srad = one/osr

            dt0  = cfl*vol0/srad

            dtaumin = min(dtaumin,dt0)
            dtaumax = max(dtaumax,dt0)
        end do
        end do
        end do
    end do

    call broadcast_timestep

    if (cdtmax > zero) then
        dtaumax = min(dtaumax, cdtmax*dtaumin)
    end if

    cflmin = large
    cflmax = small
    do nc=1,nblkcoms
        nb  =  blkcoms(nc)%nb
        top => blkcoms(nc)%top

        ni = top%nijk(1)
        nj = top%nijk(2)
        nk = top%nijk(3)

        vol => mb_vol(nb)%fld
        rdt => mb_rdt(nb)%fld
        dt  => mb_dt(nb)%fld
        src => mb_src(nb)%fld
        if (nvis > nvis_euler) then
            srv => mb_srv(nb)%fld
        end if

        do k=1,nk
        do j=1,nj
        do i=1,ni
            vol0 = vol(1)%r3d(i,j,k)

            ra = src(1)%r3d(i,j,k)
            rb = src(2)%r3d(i,j,k)
            rc = src(3)%r3d(i,j,k) + small
            if (nvis > nvis_euler) then
                ra = ra + srv(1)%r3d(i,j,k)
                rb = rb + srv(2)%r3d(i,j,k)
                rc = rc + srv(3)%r3d(i,j,k)
            end if
            osr = one/ra + one/rb + fsw_kdir/rc
            srad = one/osr

            dt0  = cfl*vol0/srad

            dt0 = min(dt0,dtaumax)
            dt0 = dt0/vol0
            dt(1)%r3d(i,j,k) = dt0
            rdt(1)%r3d(i,j,k) = rdt(1)%r3d(i,j,k) + one/dt0

            cfl0 = dt0*srad
            cflmin = min(cflmin,cfl0)
            cflmax = max(cflmax,cfl0)
        end do
        end do
        end do
    end do

end subroutine calc_timestep_localw

subroutine calc_timestep_localw_sp
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,one,small,large
    use mod_constants, only : nvis_euler
    use mod_datatypes, only : fld_array_t,top_block_t
    use mod_variables, only : nvis,fsw_kdir
    use mod_variables, only : cfl,cflmin,cflmax
    use mod_variables, only : cdtmax,dtaumin,dtaumax
    use mod_fieldvars, only : nblkcoms,blkcomssp
    use mod_fieldvars, only : mb_volsp,mb_rdt,mb_dt
    use mod_fieldvars, only : mb_src,mb_srv
    implicit none
    integer(kind_int)          :: nc,nb,i,j,k,ierr
    integer(kind_int)          :: ni,nj,nk
    type(top_block_t), pointer :: top
    real(kind_real)            :: ra,rb,rc,osr
    real(kind_real)            :: srad,vol0,dt0,cfl0
    type(fld_array_t), pointer :: vol(:),rdt(:),dt(:)
    type(fld_array_t), pointer :: src(:),srv(:)

    dtaumin = large
    dtaumax = small
    do nc=1,nblkcoms
        nb  =  blkcomssp(nc)%nb
        top => blkcomssp(nc)%top

        ni = top%nijk(1)
        nj = top%nijk(2)
        nk = top%nijk(3)

        vol => mb_volsp(nb)%fld
        rdt => mb_rdt(nb)%fld
        dt  => mb_dt(nb)%fld
        src => mb_src(nb)%fld
        if (nvis > nvis_euler) then
            srv => mb_srv(nb)%fld
        end if

        do k=1,nk
        do j=1,nj
        do i=1,ni
            vol0 = vol(1)%r3d(i,j,k)

            ra = src(1)%r3d(i,j,k)
            rb = src(2)%r3d(i,j,k)
            rc = src(3)%r3d(i,j,k) + small
            if (nvis > nvis_euler) then
                ra = ra + srv(1)%r3d(i,j,k)
                rb = rb + srv(2)%r3d(i,j,k)
                rc = rc + srv(3)%r3d(i,j,k)
            end if
            osr = one/ra + one/rb + fsw_kdir/rc
            srad = one/osr

            dt0  = cfl*vol0/srad

            dtaumin = min(dtaumin,dt0)
            dtaumax = max(dtaumax,dt0)
        end do
        end do
        end do
    end do

    call broadcast_timestep

    if (cdtmax > zero) then
        dtaumax = min(dtaumax, cdtmax*dtaumin)
    end if

    cflmin = large
    cflmax = small
    do nc=1,nblkcoms
        nb  =  blkcomssp(nc)%nb
        top => blkcomssp(nc)%top

        ni = top%nijk(1)
        nj = top%nijk(2)
        nk = top%nijk(3)

        vol => mb_volsp(nb)%fld
        rdt => mb_rdt(nb)%fld
        dt  => mb_dt(nb)%fld
        src => mb_src(nb)%fld
        if (nvis > nvis_euler) then
            srv => mb_srv(nb)%fld
        end if

        do k=1,nk
        do j=1,nj
        do i=1,ni
            vol0 = vol(1)%r3d(i,j,k)

            ra = src(1)%r3d(i,j,k)
            rb = src(2)%r3d(i,j,k)
            rc = src(3)%r3d(i,j,k) + small
            if (nvis > nvis_euler) then
                ra = ra + srv(1)%r3d(i,j,k)
                rb = rb + srv(2)%r3d(i,j,k)
                rc = rc + srv(3)%r3d(i,j,k)
            end if
            osr = one/ra + one/rb + fsw_kdir/rc
            srad = one/osr

            dt0  = cfl*vol0/srad

            dt0 = min(dt0,dtaumax)
            dt0 = dt0/vol0
            dt(1)%r3d(i,j,k) = dt0
            rdt(1)%r3d(i,j,k) = rdt(1)%r3d(i,j,k) + one/dt0

            cfl0 = dt0*srad
            cflmin = min(cflmin,cfl0)
            cflmax = max(cflmax,cfl0)
        end do
        end do
        end do
    end do

end subroutine calc_timestep_localw_sp

subroutine calc_timestep_localm
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,one,small,large
    use mod_constants, only : nvis_euler
    use mod_datatypes, only : fld_array_t,top_block_t
    use mod_variables, only : nvis,fsw_kdir
    use mod_variables, only : cfl,cflmin,cflmax
    use mod_variables, only : cdtmax,dtaumin,dtaumax
    use mod_fieldvars, only : nblkcoms,blkcoms
    use mod_fieldvars, only : mb_vol,mb_rdt,mb_dt
    use mod_fieldvars, only : mb_src,mb_srv
    implicit none
    integer(kind_int)          :: nc,nb,i,j,k,ierr
    integer(kind_int)          :: ni,nj,nk
    type(top_block_t), pointer :: top
    real(kind_real)            :: ra,rb,rc,rsq
    real(kind_real)            :: srad,vol0,dt0,cfl0
    type(fld_array_t), pointer :: vol(:),rdt(:),dt(:)
    type(fld_array_t), pointer :: src(:),srv(:)

    dtaumin = large
    dtaumax = small
    do nc=1,nblkcoms
        nb  =  blkcoms(nc)%nb
        top => blkcoms(nc)%top

        ni = top%nijk(1)
        nj = top%nijk(2)
        nk = top%nijk(3)

        vol => mb_vol(nb)%fld
        rdt => mb_rdt(nb)%fld
        dt  => mb_dt(nb)%fld
        src => mb_src(nb)%fld
        if (nvis > nvis_euler) then
            srv => mb_srv(nb)%fld
        end if

        do k=1,nk
        do j=1,nj
        do i=1,ni
            vol0 = vol(1)%r3d(i,j,k)

            ra = src(1)%r3d(i,j,k)
            rb = src(2)%r3d(i,j,k)
            rc = src(3)%r3d(i,j,k)
            rsq = sqrt(ra*ra + rb*rb + rc*rc)
            srad = rsq
            if (nvis > nvis_euler) then
                ra = srv(1)%r3d(i,j,k)
                rb = srv(2)%r3d(i,j,k)
                rc = srv(3)%r3d(i,j,k)
                rsq = sqrt(ra*ra + rb*rb + rc*rc)
                srad = max(srad,rsq)
            end if

            dt0  = cfl*vol0/srad

            dtaumin = min(dtaumin,dt0)
            dtaumax = max(dtaumax,dt0)
        end do
        end do
        end do
    end do

    call broadcast_timestep

    if (cdtmax > zero) then
        dtaumax = min(dtaumax, cdtmax*dtaumin)
    end if

    cflmin = large
    cflmax = small
    do nc=1,nblkcoms
        nb  =  blkcoms(nc)%nb
        top => blkcoms(nc)%top

        ni = top%nijk(1)
        nj = top%nijk(2)
        nk = top%nijk(3)

        vol => mb_vol(nb)%fld
        rdt => mb_rdt(nb)%fld
        dt  => mb_dt(nb)%fld
        src => mb_src(nb)%fld
        if (nvis > nvis_euler) then
            srv => mb_srv(nb)%fld
        end if

        do k=1,nk
        do j=1,nj
        do i=1,ni
            vol0 = vol(1)%r3d(i,j,k)

            ra = src(1)%r3d(i,j,k)
            rb = src(2)%r3d(i,j,k)
            rc = src(3)%r3d(i,j,k)
            rsq = sqrt(ra*ra + rb*rb + rc*rc)
            srad = rsq
            if (nvis > nvis_euler) then
                ra = srv(1)%r3d(i,j,k)
                rb = srv(2)%r3d(i,j,k)
                rc = srv(3)%r3d(i,j,k)
                rsq = sqrt(ra*ra + rb*rb + rc*rc)
                srad = max(srad,rsq)
            end if

            dt0  = cfl*vol0/srad

            dt0 = min(dt0,dtaumax)
            dt0 = dt0/vol0
            dt(1)%r3d(i,j,k) = dt0
            rdt(1)%r3d(i,j,k) = rdt(1)%r3d(i,j,k) + one/dt0

            cfl0 = dt0*srad
            cflmin = min(cflmin,cfl0)
            cflmax = max(cflmax,cfl0)
        end do
        end do
        end do
    end do

end subroutine calc_timestep_localm

subroutine calc_timestep_localm_sp
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,one,small,large
    use mod_constants, only : nvis_euler
    use mod_datatypes, only : fld_array_t,top_block_t
    use mod_variables, only : nvis,fsw_kdir
    use mod_variables, only : cfl,cflmin,cflmax
    use mod_variables, only : cdtmax,dtaumin,dtaumax
    use mod_fieldvars, only : nblkcoms,blkcomssp
    use mod_fieldvars, only : mb_volsp,mb_rdt,mb_dt
    use mod_fieldvars, only : mb_src,mb_srv
    implicit none
    integer(kind_int)          :: nc,nb,i,j,k,ierr
    integer(kind_int)          :: ni,nj,nk
    type(top_block_t), pointer :: top
    real(kind_real)            :: ra,rb,rc,rsq
    real(kind_real)            :: srad,vol0,dt0,cfl0
    type(fld_array_t), pointer :: vol(:),rdt(:),dt(:)
    type(fld_array_t), pointer :: src(:),srv(:)

    dtaumin = large
    dtaumax = small
    do nc=1,nblkcoms
        nb  =  blkcomssp(nc)%nb
        top => blkcomssp(nc)%top

        ni = top%nijk(1)
        nj = top%nijk(2)
        nk = top%nijk(3)

        vol => mb_volsp(nb)%fld
        rdt => mb_rdt(nb)%fld
        dt  => mb_dt(nb)%fld
        src => mb_src(nb)%fld
        if (nvis > nvis_euler) then
            srv => mb_srv(nb)%fld
        end if

        do k=1,nk
        do j=1,nj
        do i=1,ni
            vol0 = vol(1)%r3d(i,j,k)

            ra = src(1)%r3d(i,j,k)
            rb = src(2)%r3d(i,j,k)
            rc = src(3)%r3d(i,j,k)
            rsq = sqrt(ra*ra + rb*rb + rc*rc)
            srad = rsq
            if (nvis > nvis_euler) then
                ra = srv(1)%r3d(i,j,k)
                rb = srv(2)%r3d(i,j,k)
                rc = srv(3)%r3d(i,j,k)
                rsq = sqrt(ra*ra + rb*rb + rc*rc)
                srad = max(srad,rsq)
            end if

            dt0  = cfl*vol0/srad

            dtaumin = min(dtaumin,dt0)
            dtaumax = max(dtaumax,dt0)
        end do
        end do
        end do
    end do

    call broadcast_timestep

    if (cdtmax > zero) then
        dtaumax = min(dtaumax, cdtmax*dtaumin)
    end if

    cflmin = large
    cflmax = small
    do nc=1,nblkcoms
        nb  =  blkcomssp(nc)%nb
        top => blkcomssp(nc)%top

        ni = top%nijk(1)
        nj = top%nijk(2)
        nk = top%nijk(3)

        vol => mb_volsp(nb)%fld
        rdt => mb_rdt(nb)%fld
        dt  => mb_dt(nb)%fld
        src => mb_src(nb)%fld
        if (nvis > nvis_euler) then
            srv => mb_srv(nb)%fld
        end if

        do k=1,nk
        do j=1,nj
        do i=1,ni
            vol0 = vol(1)%r3d(i,j,k)

            ra = src(1)%r3d(i,j,k)
            rb = src(2)%r3d(i,j,k)
            rc = src(3)%r3d(i,j,k)
            rsq = sqrt(ra*ra + rb*rb + rc*rc)
            srad = rsq
            if (nvis > nvis_euler) then
                ra = srv(1)%r3d(i,j,k)
                rb = srv(2)%r3d(i,j,k)
                rc = srv(3)%r3d(i,j,k)
                rsq = sqrt(ra*ra + rb*rb + rc*rc)
                srad = max(srad,rsq)
            end if

            dt0  = cfl*vol0/srad

            dt0 = min(dt0,dtaumax)
            dt0 = dt0/vol0
            dt(1)%r3d(i,j,k) = dt0
            rdt(1)%r3d(i,j,k) = rdt(1)%r3d(i,j,k) + one/dt0

            cfl0 = dt0*srad
            cflmin = min(cflmin,cfl0)
            cflmax = max(cflmax,cfl0)
        end do
        end do
        end do
    end do

end subroutine calc_timestep_localm_sp

subroutine calc_timestep_globald
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,one,small,large
    use mod_constants, only : nvis_euler
    use mod_datatypes, only : fld_array_t,top_block_t
    use mod_variables, only : nvis,fsw_kdir
    use mod_variables, only : cfl,cflmin,cflmax,cflfac
    use mod_variables, only : cdtmax,dtaumin,dtaumax
    use mod_fieldvars, only : nblkcoms,blkcoms
    use mod_fieldvars, only : mb_vol,mb_rdt,mb_dt,mb_dst
    use mod_openmp
    implicit none
    integer(kind_int)          :: nc,nb,i,j,k,ierr
    integer(kind_int)          :: ni,nj,nk
    type(top_block_t), pointer :: top
    real(kind_real)            :: srad,vol0,dt0,cfl0,dist
    type(fld_array_t), pointer :: vol(:),rdt(:),dt(:)
    type(fld_array_t), pointer :: src(:),srv(:),dst(:)
#ifdef OMP_IMP
    integer(kind_int)          :: idt,numt
    real(kind_real), dimension(0:omp_max_num_threads) :: dmin, dmax, cmin, cmax
#endif

#ifdef OMP_IMP
    dmin(:) = large
    dmax(:) = small
#endif

    dtaumin = large
    dtaumax = small
    do nc=1,nblkcoms
        nb  =  blkcoms(nc)%nb
        top => blkcoms(nc)%top

        ni = top%nijk(1)
        nj = top%nijk(2)
        nk = top%nijk(3)

        vol => mb_vol(nb)%fld
        rdt => mb_rdt(nb)%fld
        dt  => mb_dt(nb)%fld

#ifdef OMP_IMP
!$OMP parallel private(vol0,srad,dt0,idt)
!$OMP master
        numt = omp_get_num_threads()
!$OMP end master
        idt  = omp_get_thread_num()
!$OMP do
#endif

        do k=1,nk
        do j=1,nj
        do i=1,ni
            vol0 = vol(1)%r3d(i,j,k)

            srad = rdt(1)%r3d(i,j,k)

            dt0  = cfl*vol0/srad

#ifdef OMP_IMP
            dmin(idt) = min(dmin(idt),dt0)
            dmax(idt) = max(dmax(idt),dt0)
#else
            dtaumin = min(dtaumin,dt0)
            dtaumax = max(dtaumax,dt0)
#endif
        end do
        end do
        end do
#ifdef OMP_IMP
!$OMP end do
!$OMP end parallel

        do i=0,numt
          dtaumin=min(dmin(i),dtaumin)
          dtaumax=max(dmax(i),dtaumax)
        end do
#endif
    end do

    call broadcast_timestep

    if (cdtmax > zero) then
        dtaumax = min(dtaumax, cdtmax*dtaumin)
    end if

    cflmin = large
    cflmax = small
#ifdef OMP_IMP
    cmin(:)=large
    cmax(:)=small
#endif
    do nc=1,nblkcoms
        nb  =  blkcoms(nc)%nb
        top => blkcoms(nc)%top

        ni = top%nijk(1)
        nj = top%nijk(2)
        nk = top%nijk(3)

        vol => mb_vol(nb)%fld
        rdt => mb_rdt(nb)%fld
        dt  => mb_dt(nb)%fld
        dst => mb_dst(nb)%fld
#ifdef OMP_IMP
!$OMP parallel  private(vol0,srad,dt0,idt)
!$OMP master
        numt = omp_get_num_threads()
!$OMP end master
        idt  = omp_get_thread_num()
!$OMP do
#endif

        do k=1,nk
        do j=1,nj
        do i=1,ni
            vol0 = vol(1)%r3d(i,j,k)

            srad = rdt(1)%r3d(i,j,k)
            
            dist = dst(1)%r3d(i,j,k)

            dt0  = cfl*vol0/srad

            dt0 = min(dt0,dtaumax)
            dt0 = dt0/vol0 + cfl/srad*(dist/cflfac)**2.0
            dt(1)%r3d(i,j,k) = dt0
            rdt(1)%r3d(i,j,k) = rdt(1)%r3d(i,j,k) + one/dt0

            cfl0 = dt0*srad
#ifdef OMP_IMP
            cmin(idt) = min(cmin(idt),cfl0)
            cmax(idt) = max(cmax(idt),cfl0)
#else
            cflmin = min(cflmin,cfl0)
            cflmax = max(cflmax,cfl0)
#endif
        end do
        end do
        end do
#ifdef OMP_IMP
!$OMP end do
!$OMP end parallel
        do i=0,numt
          cflmin=min(cmin(i),cflmin)
          cflmax=max(cmax(i),cflmax)
        end do
#endif
    end do

end subroutine calc_timestep_globald

subroutine calc_timestep_globald_sp
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,one,small,large
    use mod_constants, only : nvis_euler
    use mod_datatypes, only : fld_array_t,top_block_t
    use mod_variables, only : nvis,fsw_kdir
    use mod_variables, only : cfl,cflmin,cflmax,cflfac
    use mod_variables, only : cdtmax,dtaumin,dtaumax
    use mod_fieldvars, only : nblkcoms,blkcomssp
    use mod_fieldvars, only : mb_volsp,mb_rdt,mb_dt,mb_dst
    use mod_openmp
    implicit none
    integer(kind_int)          :: nc,nb,i,j,k,ierr
    integer(kind_int)          :: ni,nj,nk
    type(top_block_t), pointer :: top
    real(kind_real)            :: srad,vol0,dt0,cfl0,dist
    type(fld_array_t), pointer :: vol(:),rdt(:),dt(:)
    type(fld_array_t), pointer :: src(:),srv(:),dst(:)
#ifdef OMP_IMP
    integer(kind_int)          :: idt,numt
    real(kind_real), dimension(0:omp_max_num_threads) :: dmin, dmax, cmin, cmax
#endif

#ifdef OMP_IMP
    dmin(:) = large
    dmax(:) = small
#endif

    dtaumin = large
    dtaumax = small
    do nc=1,nblkcoms
        nb  =  blkcomssp(nc)%nb
        top => blkcomssp(nc)%top

        ni = top%nijk(1)
        nj = top%nijk(2)
        nk = top%nijk(3)

        vol => mb_volsp(nb)%fld
        rdt => mb_rdt(nb)%fld
        dt  => mb_dt(nb)%fld

#ifdef OMP_IMP
!$OMP parallel private(vol0,srad,dt0,idt)
!$OMP master
        numt = omp_get_num_threads()
!$OMP end master
        idt  = omp_get_thread_num()
!$OMP do
#endif

        do k=1,nk
        do j=1,nj
        do i=1,ni
            vol0 = vol(1)%r3d(i,j,k)

            srad = rdt(1)%r3d(i,j,k)

            dt0  = cfl*vol0/srad

#ifdef OMP_IMP
            dmin(idt) = min(dmin(idt),dt0)
            dmax(idt) = max(dmax(idt),dt0)
#else
            dtaumin = min(dtaumin,dt0)
            dtaumax = max(dtaumax,dt0)
#endif
        end do
        end do
        end do
#ifdef OMP_IMP
!$OMP end do
!$OMP end parallel

        do i=0,numt
          dtaumin=min(dmin(i),dtaumin)
          dtaumax=max(dmax(i),dtaumax)
        end do
#endif
    end do

    call broadcast_timestep

    if (cdtmax > zero) then
        dtaumax = min(dtaumax, cdtmax*dtaumin)
    end if

    cflmin = large
    cflmax = small
#ifdef OMP_IMP
    cmin(:)=large
    cmax(:)=small
#endif
    do nc=1,nblkcoms
        nb  =  blkcomssp(nc)%nb
        top => blkcomssp(nc)%top

        ni = top%nijk(1)
        nj = top%nijk(2)
        nk = top%nijk(3)

        vol => mb_volsp(nb)%fld
        rdt => mb_rdt(nb)%fld
        dt  => mb_dt(nb)%fld
        dst => mb_dst(nb)%fld
#ifdef OMP_IMP
!$OMP parallel  private(vol0,srad,dt0,idt)
!$OMP master
        numt = omp_get_num_threads()
!$OMP end master
        idt  = omp_get_thread_num()
!$OMP do
#endif

        do k=1,nk
        do j=1,nj
        do i=1,ni
            vol0 = vol(1)%r3d(i,j,k)

            srad = rdt(1)%r3d(i,j,k)
            
            dist = dst(1)%r3d(i,j,k)

            dt0  = cfl*vol0/srad

            dt0 = min(dt0,dtaumax)
            dt0 = dt0/vol0 + cfl/srad*(dist/cflfac)**2.0
            dt(1)%r3d(i,j,k) = dt0
            rdt(1)%r3d(i,j,k) = rdt(1)%r3d(i,j,k) + one/dt0

            cfl0 = dt0*srad
#ifdef OMP_IMP
            cmin(idt) = min(cmin(idt),cfl0)
            cmax(idt) = max(cmax(idt),cfl0)
#else
            cflmin = min(cflmin,cfl0)
            cflmax = max(cflmax,cfl0)
#endif
        end do
        end do
        end do
#ifdef OMP_IMP
!$OMP end do
!$OMP end parallel
        do i=0,numt
          cflmin=min(cmin(i),cflmin)
          cflmax=max(cmax(i),cflmax)
        end do
#endif
    end do

end subroutine calc_timestep_globald_sp

!������������
subroutine update
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : small,zero,one,sml_ssf
    use mod_constants, only : id_w,nsw_dir_close
    use mod_datatypes, only : fld_array_t,top_block_t
    use mod_variables, only : csafeup,nsw_kdir,nghnode
    use mod_variables, only : rlimit,plimit
    use mod_fieldvars, only : nblkcoms,blkcoms
    use mod_fieldvars, only : mb_pv,mb_qc,mb_dq
    use mod_fieldvars, only : npvs,neqn,mb_sxyz
    use mod_interface, only : assign_mb_var_uniform
    use mod_openmp
    implicit none
    integer(kind_int)          :: nc,nb,i,j,k,m,ierr
    integer(kind_int)          :: st(3),ed(3),i1,j1,k1,n
    integer(kind_int)          :: ni,nj,nk,nkst,nupd
    integer(kind_int)          :: nmius,nbijk(4)
    type(top_block_t), pointer :: top
    real(kind_real)            :: kscal,ro,ps,w0(id_w:id_w)
    real(kind_real)            :: prim(1:npvs),con(1:neqn)
    type(fld_array_t), pointer :: pv(:),qc(:),dq(:),sxyz(:)
#ifdef OMP_IMP
    integer(kind_int)          :: nnew(0:omp_max_num_threads)
    integer(kind_int)          :: numt,idt,icount
#endif

    do nc=1,nblkcoms
        nb  =  blkcoms(nc)%nb
        top => blkcoms(nc)%top

        ni = top%nijk(1)
        nj = top%nijk(2)
        nk = top%nijk(3)

        st(:) = top%ndst(:)
        ed(:) = top%nded(:)

        pv => mb_pv(nb)%fld
        qc => mb_qc(nb)%fld
        dq => mb_dq(nb)%fld

        sxyz => mb_sxyz(nb)%fld


        nmius = 0
        nbijk(:) = 0
#ifdef OMP_IMP
!$OMP parallel  private(kscal,con,ro,ps,prim,i1,j1,k1),reduction(+:n)
!$OMP master
        numt = omp_get_num_threads()
!$OMP end master
        idt  = omp_get_thread_num()
!$OMP  do
#endif
        do k=st(3),ed(3)
        do j=st(2),ed(2)
        do i=st(1),ed(1)
            if (csafeup > zero) then
                kscal = min(one,csafeup*pv(1)%r3d(i,j,k)/(abs(dq(1)%r3d(i,j,k)) + small))
            else
                kscal = one
            end if

            do m=1,neqn
                con(m) = qc(m)%r3d(i,j,k) + kscal*dq(m)%r3d(i,j,k)
            end do

            call con2prim(1,neqn,con,1,npvs,prim)

            ro = prim(1)
            ps = prim(5)

            if (isnan(ps) .or. isnan(ro)) then
                write(*,*) nb,i,j,k,ro,ps
                call error_check(1, &
                     "The data of flow field is NAN in subroutine update")
            end if

            if (ro < rlimit(1) .or. ps < plimit(1)) then
#ifdef OMP_IMP
                nnew(:) = 0
#else
                nmius = nmius + 1
                if (nmius == 1) nbijk(:) = (/nb,i,j,k/)
#endif

                n = 0
                prim(:) = zero

                i1 = i - 1
                if (max(1,min(i1,ni)) == i1) then
#ifdef OMP_IMP
                    nnew(idt) = nnew(idt) + 1
#else
                    n = n + 1
#endif
                    do m=1,npvs
                        prim(m) = prim(m) + pv(m)%r3d(i1,j,k)
                    end do
                end if

                i1 = i + 1
                if (max(1,min(i1,ni)) == i1) then
#ifdef OMP_IMP
                    nnew(idt) = nnew(idt) + 1
#else
                    n = n + 1
#endif
                    do m=1,npvs
                        prim(m) = prim(m) + pv(m)%r3d(i1,j,k)
                    end do
                end if

                j1 = j - 1
                if (max(1,min(j1,nj)) == j1) then
#ifdef OMP_IMP
                    nnew(idt) = nnew(idt) + 1
#else
                    n = n + 1
#endif
                    do m=1,npvs
                        prim(m) = prim(m) + pv(m)%r3d(i,j1,k)
                    end do
                end if

                j1 = j + 1
                if (max(1,min(j1,nj)) == j1) then
#ifdef OMP_IMP
                    nnew(idt) = nnew(idt) + 1
#else
                    n = n + 1
#endif
                    do m=1,npvs
                        prim(m) = prim(m) + pv(m)%r3d(i,j1,k)
                    end do
                end if

                k1 = k - 1
                if (max(1,min(k1,nk)) == k1) then
#ifdef OMP_IMP
                    nnew(idt) = nnew(idt) + 1
#else
                    n = n + 1
#endif
                    do m=1,npvs
                        prim(m) = prim(m) + pv(m)%r3d(i,j,k1)
                    end do
                end if

                k1 = k + 1
                if (max(1,min(k1,nk)) == k1) then
#ifdef OMP_IMP
                    nnew(idt) = nnew(idt) + 1
#else
                    n = n + 1
#endif
                    do m=1,npvs
                        prim(m) = prim(m) + pv(m)%r3d(i,j,k1)
                    end do
                end if

#ifdef OMP_IMP
                do icount=0,numt
                  n=n+nnew(icount)
                end do
#endif

                prim(:) = prim(:)/n

                call prim2con(1,npvs,prim,1,neqn,con)
            end if

            do m=1,npvs
                pv(m)%r3d(i,j,k) = prim(m)
            end do

            do m=1,neqn
                dq(m)%r3d(i,j,k) = con(m) - qc(m)%r3d(i,j,k)
                qc(m)%r3d(i,j,k) = con(m)
            end do
        end do
        end do
        end do
#ifdef OMP_IMP
!$OMP end  do
!$OMP end parallel
#endif

        if (nmius > 0) then
            write(*,'(a,4i5)') " Warning P<0 or Ro<0:",nbijk(:)
        end if

        if (nsw_kdir == nsw_dir_close) then
            nkst = top%ndst(3)
            do j=1,nj
            do i=1,ni
                pv(id_w)%r3d(i,j,nkst) = zero
                qc(id_w)%r3d(i,j,nkst) = zero
                dq(id_w)%r3d(i,j,nkst) = zero
            end do
            end do

            do m=1,neqn
                do k=1,nk
                do j=1,nj
                do i=1,ni
                    pv(m)%r3d(i,j,k) = pv(m)%r3d(i,j,nkst)
                    qc(m)%r3d(i,j,k) = qc(m)%r3d(i,j,nkst)
                    dq(m)%r3d(i,j,k) = dq(m)%r3d(i,j,nkst)
                end do
                end do
                end do
            end do
        end if
    end do

end subroutine update

!������������
subroutine update_sp
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : small,zero,one,sml_ssf
    use mod_constants, only : id_w,nsw_dir_close
    use mod_datatypes, only : fld_array_t,top_block_t
    use mod_variables, only : csafeup,nsw_kdir,nghnode
    use mod_variables, only : rlimit,plimit
    use mod_fieldvars, only : nblkcoms,blkcomssp
    use mod_fieldvars, only : mb_pv,mb_qc,mb_dq
    use mod_fieldvars, only : npvs,neqn
    use mod_openmp
    implicit none
    integer(kind_int)          :: nc,nb,i,j,k,m,ierr
    integer(kind_int)          :: st(3),ed(3),i1,j1,k1,n
    integer(kind_int)          :: ni,nj,nk,nkst,nupd
    integer(kind_int)          :: nmius,nbijk(4)
    type(top_block_t), pointer :: top
    real(kind_real)            :: kscal,ro,ps,w0(id_w:id_w)
    real(kind_real)            :: prim(1:npvs),con(1:neqn)
    type(fld_array_t), pointer :: pv(:),qc(:),dq(:)
#ifdef OMP_IMP
    integer(kind_int)          :: nnew(0:omp_max_num_threads)
    integer(kind_int)          :: numt,idt,icount
#endif

    do nc=1,nblkcoms
        nb  =  blkcomssp(nc)%nb
        top => blkcomssp(nc)%top

        ni = top%nijk(1)
        nj = top%nijk(2)
        nk = top%nijk(3)

        st(:) = top%ndst(:)
        ed(:) = top%nded(:)

        pv => mb_pv(nb)%fld
        qc => mb_qc(nb)%fld
        dq => mb_dq(nb)%fld


        nmius = 0
        nbijk(:) = 0
#ifdef OMP_IMP
!$OMP parallel  private(kscal,con,ro,ps,prim,i1,j1,k1),reduction(+:n)
!$OMP master
        numt = omp_get_num_threads()
!$OMP end master
        idt  = omp_get_thread_num()
!$OMP  do
#endif
        do k=st(3),ed(3)
        do j=st(2),ed(2)
        do i=st(1),ed(1)
            if (csafeup > zero) then
                kscal = min(one,csafeup*pv(1)%r3d(i,j,k)/(abs(dq(1)%r3d(i,j,k)) + small))
            else
                kscal = one
            end if

            do m=1,neqn
                con(m) = qc(m)%r3d(i,j,k) + kscal*dq(m)%r3d(i,j,k)
            end do

            call con2prim(1,neqn,con,1,npvs,prim)

            ro = prim(1)
            ps = prim(5)

            if (isnan(ps) .or. isnan(ro)) then
                write(*,*) nb,i,j,k,ro,ps
                call error_check(1, &
                     "The data of flow field is NAN in subroutine update")
            end if

            if (ro < rlimit(1) .or. ps < plimit(1)) then
#ifdef OMP_IMP
                nnew(:) = 0
#else
                nmius = nmius + 1
                if (nmius == 1) nbijk(:) = (/nb,i,j,k/)
#endif

                n = 0
                prim(:) = zero

                i1 = i - 1
                if (max(1,min(i1,ni)) == i1) then
#ifdef OMP_IMP
                    nnew(idt) = nnew(idt) + 1
#else
                    n = n + 1
#endif
                    do m=1,npvs
                        prim(m) = prim(m) + pv(m)%r3d(i1,j,k)
                    end do
                end if

                i1 = i + 1
                if (max(1,min(i1,ni)) == i1) then
#ifdef OMP_IMP
                    nnew(idt) = nnew(idt) + 1
#else
                    n = n + 1
#endif
                    do m=1,npvs
                        prim(m) = prim(m) + pv(m)%r3d(i1,j,k)
                    end do
                end if

                j1 = j - 1
                if (max(1,min(j1,nj)) == j1) then
#ifdef OMP_IMP
                    nnew(idt) = nnew(idt) + 1
#else
                    n = n + 1
#endif
                    do m=1,npvs
                        prim(m) = prim(m) + pv(m)%r3d(i,j1,k)
                    end do
                end if

                j1 = j + 1
                if (max(1,min(j1,nj)) == j1) then
#ifdef OMP_IMP
                    nnew(idt) = nnew(idt) + 1
#else
                    n = n + 1
#endif
                    do m=1,npvs
                        prim(m) = prim(m) + pv(m)%r3d(i,j1,k)
                    end do
                end if

                k1 = k - 1
                if (max(1,min(k1,nk)) == k1) then
#ifdef OMP_IMP
                    nnew(idt) = nnew(idt) + 1
#else
                    n = n + 1
#endif
                    do m=1,npvs
                        prim(m) = prim(m) + pv(m)%r3d(i,j,k1)
                    end do
                end if

                k1 = k + 1
                if (max(1,min(k1,nk)) == k1) then
#ifdef OMP_IMP
                    nnew(idt) = nnew(idt) + 1
#else
                    n = n + 1
#endif
                    do m=1,npvs
                        prim(m) = prim(m) + pv(m)%r3d(i,j,k1)
                    end do
                end if

#ifdef OMP_IMP
                do icount=0,numt
                  n=n+nnew(icount)
                end do
#endif

                prim(:) = prim(:)/n

                call prim2con(1,npvs,prim,1,neqn,con)
            end if

            do m=1,npvs
                pv(m)%r3d(i,j,k) = prim(m)
            end do

            do m=1,neqn
                dq(m)%r3d(i,j,k) = con(m) - qc(m)%r3d(i,j,k)
                qc(m)%r3d(i,j,k) = con(m)
            end do
        end do
        end do
        end do
#ifdef OMP_IMP
!$OMP end  do
!$OMP end parallel
#endif

        if (nmius > 0) then
            write(*,'(a,4i5)') " Warning P<0 or Ro<0:",nbijk(:)
        end if

        if (nsw_kdir == nsw_dir_close) then
            nkst = top%ndst(3)
            do j=1,nj
            do i=1,ni
                pv(id_w)%r3d(i,j,nkst) = zero
                qc(id_w)%r3d(i,j,nkst) = zero
                dq(id_w)%r3d(i,j,nkst) = zero
            end do
            end do

            do m=1,neqn
                do k=1,nk
                do j=1,nj
                do i=1,ni
                    pv(m)%r3d(i,j,k) = pv(m)%r3d(i,j,nkst)
                    qc(m)%r3d(i,j,k) = qc(m)%r3d(i,j,nkst)
                    dq(m)%r3d(i,j,k) = dq(m)%r3d(i,j,nkst)
                end do
                end do
                end do
            end do
        end if
    end do

end subroutine update_sp
    
    
!������������������Ԥ������������������Ϊdro,du,dv,dw,dT
subroutine update_prec
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : small,zero,one,sml_ssf
    use mod_constants, only : id_w,nsw_dir_close
    use mod_datatypes, only : fld_array_t,top_block_t
    use mod_variables, only : csafeup,nsw_kdir,nghnode,refbeta
    use mod_variables, only : rlimit,plimit
    use mod_fieldvars, only : nblkcoms,blkcoms
    use mod_fieldvars, only : mb_pv,mb_qc,mb_dq
    use mod_fieldvars, only : npvs,neqn,mb_sxyz
    use mod_interface, only : assign_mb_var_uniform
    use mod_openmp
    implicit none
    integer(kind_int)          :: nc,nb,i,j,k,m,ierr
    integer(kind_int)          :: st(3),ed(3),i1,j1,k1,n
    integer(kind_int)          :: ni,nj,nk,nkst,nupd
    integer(kind_int)          :: nmius,nbijk(4)
    type(top_block_t), pointer :: top
    real(kind_real)            :: kscal,ro,ps,w0(id_w:id_w)
    real(kind_real)            :: prim(1:npvs),con(1:neqn)
    real(kind_real)            :: dro,du,dv,dw,dT,dp,pressure,T
    type(fld_array_t), pointer :: pv(:),qc(:),dq(:),sxyz(:)
#ifdef OMP_IMP
    integer(kind_int)          :: nnew(0:omp_max_num_threads)
    integer(kind_int)          :: numt,idt,icount
#endif

    do nc=1,nblkcoms
        nb  =  blkcoms(nc)%nb
        top => blkcoms(nc)%top

        ni = top%nijk(1)
        nj = top%nijk(2)
        nk = top%nijk(3)

        st(:) = top%ndst(:)
        ed(:) = top%nded(:)

        pv => mb_pv(nb)%fld
        qc => mb_qc(nb)%fld
        dq => mb_dq(nb)%fld

        sxyz => mb_sxyz(nb)%fld


        nmius = 0
        nbijk(:) = 0
#ifdef OMP_IMP
!$OMP parallel  private(kscal,con,ro,ps,prim,i1,j1,k1),reduction(+:n)
!$OMP master
        numt = omp_get_num_threads()
!$OMP end master
        idt  = omp_get_thread_num()
!$OMP  do
#endif
        do k=st(3),ed(3)
        do j=st(2),ed(2)
        do i=st(1),ed(1)
            if (csafeup > zero) then
                kscal = min(one,csafeup*pv(1)%r3d(i,j,k)/(abs(dq(1)%r3d(i,j,k)) + small))
            else
                kscal = one
            end if
            
            dp = dq(1)%r3d(i,j,k)
            du = dq(2)%r3d(i,j,k)
            dv = dq(3)%r3d(i,j,k)
            dw = dq(4)%r3d(i,j,k)
            dT = dq(5)%r3d(i,j,k)
            
            !����dp,du,dv,dw,dT����
            !ת��������dro,du,dv,dw,dp
            pressure = pv(5)%r3d(i,j,k)
            ro = pv(1)%r3d(i,j,k)
            T = pressure/(refbeta*ro)
            dro = ro/pressure*dp - ro/T*dT
            
            
            prim(1) = pv(1)%r3d(i,j,k) + kscal*dro
            prim(2) = pv(2)%r3d(i,j,k) + kscal*du
            prim(3) = pv(3)%r3d(i,j,k) + kscal*dv
            prim(4) = pv(4)%r3d(i,j,k) + kscal*dw
            prim(5) = pv(5)%r3d(i,j,k) + kscal*dp
            
            
            !����ԭʼ����
            !do m=1,neqn
            !    prim(m) = pv(m)%r3d(i,j,k) + kscal*dq(m)%r3d(i,j,k)
            !end do

            !call prim2con(1,neqn,con,1,npvs,prim)

            ro = prim(1)
            ps = prim(5)

            if (isnan(ps) .or. isnan(ro)) then
                write(*,*) nb,i,j,k,ro,ps
                call error_check(1, &
                     "The data of flow field is NAN in subroutine update_prec")
            end if

            if (ro < rlimit(1) .or. ps < plimit(1)) then
#ifdef OMP_IMP
                nnew(:) = 0
#else
                nmius = nmius + 1
                if (nmius == 1) nbijk(:) = (/nb,i,j,k/)
#endif

                n = 0
                prim(:) = zero

                i1 = i - 1
                if (max(1,min(i1,ni)) == i1) then
#ifdef OMP_IMP
                    nnew(idt) = nnew(idt) + 1
#else
                    n = n + 1
#endif
                    do m=1,npvs
                        prim(m) = prim(m) + pv(m)%r3d(i1,j,k)
                    end do
                end if

                i1 = i + 1
                if (max(1,min(i1,ni)) == i1) then
#ifdef OMP_IMP
                    nnew(idt) = nnew(idt) + 1
#else
                    n = n + 1
#endif
                    do m=1,npvs
                        prim(m) = prim(m) + pv(m)%r3d(i1,j,k)
                    end do
                end if

                j1 = j - 1
                if (max(1,min(j1,nj)) == j1) then
#ifdef OMP_IMP
                    nnew(idt) = nnew(idt) + 1
#else
                    n = n + 1
#endif
                    do m=1,npvs
                        prim(m) = prim(m) + pv(m)%r3d(i,j1,k)
                    end do
                end if

                j1 = j + 1
                if (max(1,min(j1,nj)) == j1) then
#ifdef OMP_IMP
                    nnew(idt) = nnew(idt) + 1
#else
                    n = n + 1
#endif
                    do m=1,npvs
                        prim(m) = prim(m) + pv(m)%r3d(i,j1,k)
                    end do
                end if

                k1 = k - 1
                if (max(1,min(k1,nk)) == k1) then
#ifdef OMP_IMP
                    nnew(idt) = nnew(idt) + 1
#else
                    n = n + 1
#endif
                    do m=1,npvs
                        prim(m) = prim(m) + pv(m)%r3d(i,j,k1)
                    end do
                end if

                k1 = k + 1
                if (max(1,min(k1,nk)) == k1) then
#ifdef OMP_IMP
                    nnew(idt) = nnew(idt) + 1
#else
                    n = n + 1
#endif
                    do m=1,npvs
                        prim(m) = prim(m) + pv(m)%r3d(i,j,k1)
                    end do
                end if

#ifdef OMP_IMP
                do icount=0,numt
                  n=n+nnew(icount)
                end do
#endif

                prim(:) = prim(:)/n

                !call prim2con(1,npvs,prim,1,neqn,con)
            end if
            
            call prim2con(1,npvs,prim,1,neqn,con)

            !�غ������ֵ
            do m=1,neqn
                qc(m)%r3d(i,j,k) = con(m)
            end do
            
            do m=1,npvs
                !dq(m)%r3d(i,j,k) = prim(m) - pv(m)%r3d(i,j,k)
                pv(m)%r3d(i,j,k) = prim(m)
            end do

            !do m=1,neqn
            !    dq(m)%r3d(i,j,k) = con(m) - qc(m)%r3d(i,j,k)
            !    qc(m)%r3d(i,j,k) = con(m)
            !end do
        end do
        end do
        end do
#ifdef OMP_IMP
!$OMP end  do
!$OMP end parallel
#endif

        if (nmius > 0) then
            write(*,'(a,4i5)') " Warning P<0 or Ro<0:",nbijk(:)
        end if

        if (nsw_kdir == nsw_dir_close) then
            nkst = top%ndst(3)
            do j=1,nj
            do i=1,ni
                pv(id_w)%r3d(i,j,nkst) = zero
                qc(id_w)%r3d(i,j,nkst) = zero
                dq(id_w)%r3d(i,j,nkst) = zero
            end do
            end do

            do m=1,neqn
                do k=1,nk
                do j=1,nj
                do i=1,ni
                    pv(m)%r3d(i,j,k) = pv(m)%r3d(i,j,nkst)
                    qc(m)%r3d(i,j,k) = qc(m)%r3d(i,j,nkst)
                    dq(m)%r3d(i,j,k) = dq(m)%r3d(i,j,nkst)
                end do
                end do
                end do
            end do
        end if
    end do

end subroutine update_prec

subroutine update_prec_sp
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : small,zero,one,sml_ssf
    use mod_constants, only : id_w,nsw_dir_close
    use mod_datatypes, only : fld_array_t,top_block_t
    use mod_variables, only : csafeup,nsw_kdir,nghnode,refbeta
    use mod_variables, only : rlimit,plimit
    use mod_fieldvars, only : nblkcoms,blkcomssp
    use mod_fieldvars, only : mb_pv,mb_qc,mb_dq
    use mod_fieldvars, only : npvs,neqn,mb_sxyzsp
    use mod_interface, only : assign_mb_var_uniform
    use mod_openmp
    implicit none
    integer(kind_int)          :: nc,nb,i,j,k,m,ierr
    integer(kind_int)          :: st(3),ed(3),i1,j1,k1,n
    integer(kind_int)          :: ni,nj,nk,nkst,nupd
    integer(kind_int)          :: nmius,nbijk(4)
    type(top_block_t), pointer :: top
    real(kind_real)            :: kscal,ro,ps,w0(id_w:id_w)
    real(kind_real)            :: prim(1:npvs),con(1:neqn)
    real(kind_real)            :: dro,du,dv,dw,dT,dp,pressure,T
    type(fld_array_t), pointer :: pv(:),qc(:),dq(:),sxyz(:)
#ifdef OMP_IMP
    integer(kind_int)          :: nnew(0:omp_max_num_threads)
    integer(kind_int)          :: numt,idt,icount
#endif

    do nc=1,nblkcoms
        nb  =  blkcomssp(nc)%nb
        top => blkcomssp(nc)%top

        ni = top%nijk(1)
        nj = top%nijk(2)
        nk = top%nijk(3)

        st(:) = top%ndst(:)
        ed(:) = top%nded(:)

        pv => mb_pv(nb)%fld
        qc => mb_qc(nb)%fld
        dq => mb_dq(nb)%fld

        sxyz => mb_sxyzsp(nb)%fld


        nmius = 0
        nbijk(:) = 0
#ifdef OMP_IMP
!$OMP parallel  private(kscal,con,ro,ps,prim,i1,j1,k1),reduction(+:n)
!$OMP master
        numt = omp_get_num_threads()
!$OMP end master
        idt  = omp_get_thread_num()
!$OMP  do
#endif
        do k=st(3),ed(3)
        do j=st(2),ed(2)
        do i=st(1),ed(1)
            if (csafeup > zero) then
                kscal = min(one,csafeup*pv(1)%r3d(i,j,k)/(abs(dq(1)%r3d(i,j,k)) + small))
            else
                kscal = one
            end if
            
            dp = dq(1)%r3d(i,j,k)
            du = dq(2)%r3d(i,j,k)
            dv = dq(3)%r3d(i,j,k)
            dw = dq(4)%r3d(i,j,k)
            dT = dq(5)%r3d(i,j,k)
            
            !����dp,du,dv,dw,dT����
            !ת��������dro,du,dv,dw,dp
            pressure = pv(5)%r3d(i,j,k)
            ro = pv(1)%r3d(i,j,k)
            T = pressure/(refbeta*ro)
            dro = ro/pressure*dp - ro/T*dT
            
            
            prim(1) = pv(1)%r3d(i,j,k) + kscal*dro
            prim(2) = pv(2)%r3d(i,j,k) + kscal*du
            prim(3) = pv(3)%r3d(i,j,k) + kscal*dv
            prim(4) = pv(4)%r3d(i,j,k) + kscal*dw
            prim(5) = pv(5)%r3d(i,j,k) + kscal*dp
            
            
            !����ԭʼ����
            !do m=1,neqn
            !    prim(m) = pv(m)%r3d(i,j,k) + kscal*dq(m)%r3d(i,j,k)
            !end do

            !call prim2con(1,neqn,con,1,npvs,prim)

            ro = prim(1)
            ps = prim(5)

            if (isnan(ps) .or. isnan(ro)) then
                write(*,*) nb,i,j,k,ro,ps
                call error_check(1, &
                     "The data of flow field is NAN in subroutine update_prec")
            end if

            if (ro < rlimit(1) .or. ps < plimit(1)) then
#ifdef OMP_IMP
                nnew(:) = 0
#else
                nmius = nmius + 1
                if (nmius == 1) nbijk(:) = (/nb,i,j,k/)
#endif

                n = 0
                prim(:) = zero

                i1 = i - 1
                if (max(1,min(i1,ni)) == i1) then
#ifdef OMP_IMP
                    nnew(idt) = nnew(idt) + 1
#else
                    n = n + 1
#endif
                    do m=1,npvs
                        prim(m) = prim(m) + pv(m)%r3d(i1,j,k)
                    end do
                end if

                i1 = i + 1
                if (max(1,min(i1,ni)) == i1) then
#ifdef OMP_IMP
                    nnew(idt) = nnew(idt) + 1
#else
                    n = n + 1
#endif
                    do m=1,npvs
                        prim(m) = prim(m) + pv(m)%r3d(i1,j,k)
                    end do
                end if

                j1 = j - 1
                if (max(1,min(j1,nj)) == j1) then
#ifdef OMP_IMP
                    nnew(idt) = nnew(idt) + 1
#else
                    n = n + 1
#endif
                    do m=1,npvs
                        prim(m) = prim(m) + pv(m)%r3d(i,j1,k)
                    end do
                end if

                j1 = j + 1
                if (max(1,min(j1,nj)) == j1) then
#ifdef OMP_IMP
                    nnew(idt) = nnew(idt) + 1
#else
                    n = n + 1
#endif
                    do m=1,npvs
                        prim(m) = prim(m) + pv(m)%r3d(i,j1,k)
                    end do
                end if

                k1 = k - 1
                if (max(1,min(k1,nk)) == k1) then
#ifdef OMP_IMP
                    nnew(idt) = nnew(idt) + 1
#else
                    n = n + 1
#endif
                    do m=1,npvs
                        prim(m) = prim(m) + pv(m)%r3d(i,j,k1)
                    end do
                end if

                k1 = k + 1
                if (max(1,min(k1,nk)) == k1) then
#ifdef OMP_IMP
                    nnew(idt) = nnew(idt) + 1
#else
                    n = n + 1
#endif
                    do m=1,npvs
                        prim(m) = prim(m) + pv(m)%r3d(i,j,k1)
                    end do
                end if

#ifdef OMP_IMP
                do icount=0,numt
                  n=n+nnew(icount)
                end do
#endif

                prim(:) = prim(:)/n

                !call prim2con(1,npvs,prim,1,neqn,con)
            end if
            
            call prim2con(1,npvs,prim,1,neqn,con)

            !�غ������ֵ
            do m=1,neqn
                qc(m)%r3d(i,j,k) = con(m)
            end do
            
            do m=1,npvs
                !dq(m)%r3d(i,j,k) = prim(m) - pv(m)%r3d(i,j,k)
                pv(m)%r3d(i,j,k) = prim(m)
            end do

            !do m=1,neqn
            !    dq(m)%r3d(i,j,k) = con(m) - qc(m)%r3d(i,j,k)
            !    qc(m)%r3d(i,j,k) = con(m)
            !end do
        end do
        end do
        end do
#ifdef OMP_IMP
!$OMP end  do
!$OMP end parallel
#endif

        if (nmius > 0) then
            write(*,'(a,4i5)') " Warning P<0 or Ro<0:",nbijk(:)
        end if

        if (nsw_kdir == nsw_dir_close) then
            nkst = top%ndst(3)
            do j=1,nj
            do i=1,ni
                pv(id_w)%r3d(i,j,nkst) = zero
                qc(id_w)%r3d(i,j,nkst) = zero
                dq(id_w)%r3d(i,j,nkst) = zero
            end do
            end do

            do m=1,neqn
                do k=1,nk
                do j=1,nj
                do i=1,ni
                    pv(m)%r3d(i,j,k) = pv(m)%r3d(i,j,nkst)
                    qc(m)%r3d(i,j,k) = qc(m)%r3d(i,j,nkst)
                    dq(m)%r3d(i,j,k) = dq(m)%r3d(i,j,nkst)
                end do
                end do
                end do
            end do
        end if
    end do

end subroutine update_prec_sp
    
subroutine update_scmp
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : small,zero,one,sml_ssf
    use mod_constants, only : id_w,nsw_dir_close
    use mod_datatypes, only : fld_array_t,top_block_t
    use mod_variables, only : csafeup,nsw_kdir,nghnode
    use mod_variables, only : rlimit,plimit
    use mod_fieldvars, only : nblkcoms,blkcoms
    use mod_fieldvars, only : mb_pv,mb_qc,mb_dq
    use mod_fieldvars, only : npvs,neqn,mb_sxyz
    use mod_interface, only : assign_mb_var_uniform
    use mod_openmp
    use mod_variables, only : gamma,moo,poo     !> Done
    use mod_constants, only : SCMP_sigma        !> Done
    implicit none
    integer(kind_int)          :: nc,nb,i,j,k,m,ierr
    integer(kind_int)          :: st(3),ed(3),i1,j1,k1,n
    integer(kind_int)          :: ni,nj,nk,nkst,nupd
    integer(kind_int)          :: nmius,nbijk(4)
    type(top_block_t), pointer :: top
    real(kind_real)            :: kscal,ro,ps,w0(id_w:id_w), pPrime,pres
    real(kind_real)            :: prim(1:npvs),con(1:neqn)
    type(fld_array_t), pointer :: pv(:),qc(:),dq(:),sxyz(:)
#ifdef OMP_IMP
    integer(kind_int)          :: nnew(0:omp_max_num_threads)
    integer(kind_int)          :: numt,idt,icount
#endif

    do nc=1,nblkcoms
        nb  =  blkcoms(nc)%nb
        top => blkcoms(nc)%top

        ni = top%nijk(1)
        nj = top%nijk(2)
        nk = top%nijk(3)

        st(:) = top%ndst(:)
        ed(:) = top%nded(:)

        pv => mb_pv(nb)%fld
        qc => mb_qc(nb)%fld
        dq => mb_dq(nb)%fld

        sxyz => mb_sxyz(nb)%fld


        nmius = 0
        nbijk(:) = 0
#ifdef OMP_IMP
!$OMP parallel  private(kscal,con,ro,ps,prim,i1,j1,k1),reduction(+:n)
!$OMP master
        numt = omp_get_num_threads()
!$OMP end master
        idt  = omp_get_thread_num()
!$OMP  do
#endif
        do k=st(3),ed(3)
        do j=st(2),ed(2)
        do i=st(1),ed(1)
            if (csafeup > zero) then
                kscal = min(one,csafeup*pv(1)%r3d(i,j,k)/(abs(dq(1)%r3d(i,j,k)) + small))
            else
                kscal = one
            end if

            
            pPrime = pv(5)%r3d(i,j,k) - poo + kscal*dq(1)%r3d(i,j,k)    !> Done: p' = p' + ��p'
            ro = (1 + gamma*moo*moo*pPrime)**(1.0/SCMP_sigma)           !> Done: p' --> ��
            pres = pPrime + poo
            
            con(1) = ro
            do m=2,4
                con(m) = qc(m)%r3d(i,j,k) + kscal*dq(m)%r3d(i,j,k)    !> Done: ��u = ��u + ����u
            end do
            con(5) = pres/(gamma - 1.) + 0.5*(con(2)*con(2) + con(3)*con(3) + con(4)*con(4))/ro

            call con2prim(1,neqn,con,1,npvs,prim)

            ro = prim(1)
            ps = prim(5)

            if (isnan(ps) .or. isnan(ro)) then
                write(*,*) nb,i,j,k,ro,ps
                call error_check(1, &
                     "The data of flow field is NAN in subroutine update")
            end if

            if (ro < rlimit(1) .or. ps < plimit(1)) then
#ifdef OMP_IMP
                nnew(:) = 0
#else
                nmius = nmius + 1
                if (nmius == 1) nbijk(:) = (/nb,i,j,k/)
#endif

                n = 0
                prim(:) = zero

                i1 = i - 1
                if (max(1,min(i1,ni)) == i1) then
#ifdef OMP_IMP
                    nnew(idt) = nnew(idt) + 1
#else
                    n = n + 1
#endif
                    do m=1,npvs
                        prim(m) = prim(m) + pv(m)%r3d(i1,j,k)
                    end do
                end if

                i1 = i + 1
                if (max(1,min(i1,ni)) == i1) then
#ifdef OMP_IMP
                    nnew(idt) = nnew(idt) + 1
#else
                    n = n + 1
#endif
                    do m=1,npvs
                        prim(m) = prim(m) + pv(m)%r3d(i1,j,k)
                    end do
                end if

                j1 = j - 1
                if (max(1,min(j1,nj)) == j1) then
#ifdef OMP_IMP
                    nnew(idt) = nnew(idt) + 1
#else
                    n = n + 1
#endif
                    do m=1,npvs
                        prim(m) = prim(m) + pv(m)%r3d(i,j1,k)
                    end do
                end if

                j1 = j + 1
                if (max(1,min(j1,nj)) == j1) then
#ifdef OMP_IMP
                    nnew(idt) = nnew(idt) + 1
#else
                    n = n + 1
#endif
                    do m=1,npvs
                        prim(m) = prim(m) + pv(m)%r3d(i,j1,k)
                    end do
                end if

                k1 = k - 1
                if (max(1,min(k1,nk)) == k1) then
#ifdef OMP_IMP
                    nnew(idt) = nnew(idt) + 1
#else
                    n = n + 1
#endif
                    do m=1,npvs
                        prim(m) = prim(m) + pv(m)%r3d(i,j,k1)
                    end do
                end if

                k1 = k + 1
                if (max(1,min(k1,nk)) == k1) then
#ifdef OMP_IMP
                    nnew(idt) = nnew(idt) + 1
#else
                    n = n + 1
#endif
                    do m=1,npvs
                        prim(m) = prim(m) + pv(m)%r3d(i,j,k1)
                    end do
                end if

#ifdef OMP_IMP
                do icount=0,numt
                  n=n+nnew(icount)
                end do
#endif

                prim(:) = prim(:)/n

                call prim2con(1,npvs,prim,1,neqn,con)
            end if

            do m=1,npvs
                pv(m)%r3d(i,j,k) = prim(m)
            end do

            do m=1,neqn
                dq(m)%r3d(i,j,k) = con(m) - qc(m)%r3d(i,j,k)
                qc(m)%r3d(i,j,k) = con(m)
            end do
        end do
        end do
        end do
#ifdef OMP_IMP
!$OMP end  do
!$OMP end parallel
#endif

        if (nmius > 0) then
            write(*,'(a,4i5)') " Warning P<0 or Ro<0:",nbijk(:)
        end if

        if (nsw_kdir == nsw_dir_close) then
            nkst = top%ndst(3)
            do j=1,nj
            do i=1,ni
                pv(id_w)%r3d(i,j,nkst) = zero
                qc(id_w)%r3d(i,j,nkst) = zero
                dq(id_w)%r3d(i,j,nkst) = zero
            end do
            end do

            do m=1,neqn
                do k=1,nk
                do j=1,nj
                do i=1,ni
                    pv(m)%r3d(i,j,k) = pv(m)%r3d(i,j,nkst)
                    qc(m)%r3d(i,j,k) = qc(m)%r3d(i,j,nkst)
                    dq(m)%r3d(i,j,k) = dq(m)%r3d(i,j,nkst)
                end do
                end do
                end do
            end do
        end if
    end do

end subroutine update_scmp

subroutine update_scmp_sp
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : small,zero,one,sml_ssf
    use mod_constants, only : id_w,nsw_dir_close
    use mod_datatypes, only : fld_array_t,top_block_t
    use mod_variables, only : csafeup,nsw_kdir,nghnode
    use mod_variables, only : rlimit,plimit
    use mod_fieldvars, only : nblkcoms,blkcomssp
    use mod_fieldvars, only : mb_pv,mb_qc,mb_dq
    use mod_fieldvars, only : npvs,neqn,mb_sxyzsp
    use mod_interface, only : assign_mb_var_uniform
    use mod_openmp
    use mod_variables, only : gamma,moo,poo     !> Done
    use mod_constants, only : SCMP_sigma        !> Done
    implicit none
    integer(kind_int)          :: nc,nb,i,j,k,m,ierr
    integer(kind_int)          :: st(3),ed(3),i1,j1,k1,n
    integer(kind_int)          :: ni,nj,nk,nkst,nupd
    integer(kind_int)          :: nmius,nbijk(4)
    type(top_block_t), pointer :: top
    real(kind_real)            :: kscal,ro,ps,w0(id_w:id_w), pPrime,pres
    real(kind_real)            :: prim(1:npvs),con(1:neqn)
    type(fld_array_t), pointer :: pv(:),qc(:),dq(:),sxyz(:)
#ifdef OMP_IMP
    integer(kind_int)          :: nnew(0:omp_max_num_threads)
    integer(kind_int)          :: numt,idt,icount
#endif

    do nc=1,nblkcoms
        nb  =  blkcomssp(nc)%nb
        top => blkcomssp(nc)%top

        ni = top%nijk(1)
        nj = top%nijk(2)
        nk = top%nijk(3)

        st(:) = top%ndst(:)
        ed(:) = top%nded(:)

        pv => mb_pv(nb)%fld
        qc => mb_qc(nb)%fld
        dq => mb_dq(nb)%fld

        sxyz => mb_sxyzsp(nb)%fld


        nmius = 0
        nbijk(:) = 0
#ifdef OMP_IMP
!$OMP parallel  private(kscal,con,ro,ps,prim,i1,j1,k1),reduction(+:n)
!$OMP master
        numt = omp_get_num_threads()
!$OMP end master
        idt  = omp_get_thread_num()
!$OMP  do
#endif
        do k=st(3),ed(3)
        do j=st(2),ed(2)
        do i=st(1),ed(1)
            if (csafeup > zero) then
                kscal = min(one,csafeup*pv(1)%r3d(i,j,k)/(abs(dq(1)%r3d(i,j,k)) + small))
            else
                kscal = one
            end if

            
            pPrime = pv(5)%r3d(i,j,k) - poo + kscal*dq(1)%r3d(i,j,k)    !> Done: p' = p' + ��p'
            ro = (1 + gamma*moo*moo*pPrime)**(1.0/SCMP_sigma)           !> Done: p' --> ��
            pres = pPrime + poo
            
            con(1) = ro
            do m=2,4
                con(m) = qc(m)%r3d(i,j,k) + kscal*dq(m)%r3d(i,j,k)    !> Done: ��u = ��u + ����u
            end do
            con(5) = pres/(gamma - 1.) + 0.5*(con(2)*con(2) + con(3)*con(3) + con(4)*con(4))/ro

            call con2prim(1,neqn,con,1,npvs,prim)

            ro = prim(1)
            ps = prim(5)

            if (isnan(ps) .or. isnan(ro)) then
                write(*,*) nb,i,j,k,ro,ps
                call error_check(1, &
                     "The data of flow field is NAN in subroutine update")
            end if

            if (ro < rlimit(1) .or. ps < plimit(1)) then
#ifdef OMP_IMP
                nnew(:) = 0
#else
                nmius = nmius + 1
                if (nmius == 1) nbijk(:) = (/nb,i,j,k/)
#endif

                n = 0
                prim(:) = zero

                i1 = i - 1
                if (max(1,min(i1,ni)) == i1) then
#ifdef OMP_IMP
                    nnew(idt) = nnew(idt) + 1
#else
                    n = n + 1
#endif
                    do m=1,npvs
                        prim(m) = prim(m) + pv(m)%r3d(i1,j,k)
                    end do
                end if

                i1 = i + 1
                if (max(1,min(i1,ni)) == i1) then
#ifdef OMP_IMP
                    nnew(idt) = nnew(idt) + 1
#else
                    n = n + 1
#endif
                    do m=1,npvs
                        prim(m) = prim(m) + pv(m)%r3d(i1,j,k)
                    end do
                end if

                j1 = j - 1
                if (max(1,min(j1,nj)) == j1) then
#ifdef OMP_IMP
                    nnew(idt) = nnew(idt) + 1
#else
                    n = n + 1
#endif
                    do m=1,npvs
                        prim(m) = prim(m) + pv(m)%r3d(i,j1,k)
                    end do
                end if

                j1 = j + 1
                if (max(1,min(j1,nj)) == j1) then
#ifdef OMP_IMP
                    nnew(idt) = nnew(idt) + 1
#else
                    n = n + 1
#endif
                    do m=1,npvs
                        prim(m) = prim(m) + pv(m)%r3d(i,j1,k)
                    end do
                end if

                k1 = k - 1
                if (max(1,min(k1,nk)) == k1) then
#ifdef OMP_IMP
                    nnew(idt) = nnew(idt) + 1
#else
                    n = n + 1
#endif
                    do m=1,npvs
                        prim(m) = prim(m) + pv(m)%r3d(i,j,k1)
                    end do
                end if

                k1 = k + 1
                if (max(1,min(k1,nk)) == k1) then
#ifdef OMP_IMP
                    nnew(idt) = nnew(idt) + 1
#else
                    n = n + 1
#endif
                    do m=1,npvs
                        prim(m) = prim(m) + pv(m)%r3d(i,j,k1)
                    end do
                end if

#ifdef OMP_IMP
                do icount=0,numt
                  n=n+nnew(icount)
                end do
#endif

                prim(:) = prim(:)/n

                call prim2con(1,npvs,prim,1,neqn,con)
            end if

            do m=1,npvs
                pv(m)%r3d(i,j,k) = prim(m)
            end do

            do m=1,neqn
                dq(m)%r3d(i,j,k) = con(m) - qc(m)%r3d(i,j,k)
                qc(m)%r3d(i,j,k) = con(m)
            end do
        end do
        end do
        end do
#ifdef OMP_IMP
!$OMP end  do
!$OMP end parallel
#endif

        if (nmius > 0) then
            write(*,'(a,4i5)') " Warning P<0 or Ro<0:",nbijk(:)
        end if

        if (nsw_kdir == nsw_dir_close) then
            nkst = top%ndst(3)
            do j=1,nj
            do i=1,ni
                pv(id_w)%r3d(i,j,nkst) = zero
                qc(id_w)%r3d(i,j,nkst) = zero
                dq(id_w)%r3d(i,j,nkst) = zero
            end do
            end do

            do m=1,neqn
                do k=1,nk
                do j=1,nj
                do i=1,ni
                    pv(m)%r3d(i,j,k) = pv(m)%r3d(i,j,nkst)
                    qc(m)%r3d(i,j,k) = qc(m)%r3d(i,j,nkst)
                    dq(m)%r3d(i,j,k) = dq(m)%r3d(i,j,nkst)
                end do
                end do
                end do
            end do
        end if
    end do

end subroutine update_scmp_sp

subroutine residual(nsw)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : small,zero,one
    use mod_constants, only : ntimeadv_steady
    use mod_datatypes, only : top_block_t,fld_array_t
    use mod_variables, only : dtime,respos
    use mod_variables, only : resmax,restot,resdts
    use mod_fieldvars, only : nblkcoms,blkcoms,neqn
    use mod_fieldvars, only : mb_dt,mb_vol,mb_dq
    use mod_fieldvars, only : mb_qc,mb_qm
    use mod_constants, only : nscmp_non
    use mod_variables, only : nscmp
    use mod_parallels
    use mod_openmp
    implicit none
    integer(kind_int), intent(in) :: nsw
    integer(kind_int)          :: nc,nb,i,j,k,m,ierr
    integer(kind_int)          :: ni,nj,nk
    integer(kind_int)          :: nRes
    type(top_block_t), pointer :: top
    real(kind_real)            :: dt0,vol0,dq0
    real(kind_real)            :: odtvol,rdq,rdq2
    real(kind_real)            :: dq_dts,rdq_dts
    type(fld_array_t), pointer :: dt(:),vol(:),dq(:)
    type(fld_array_t), pointer :: qc(:),qm(:)
#ifdef PARALLEL
    integer(kind_int) :: id_max
    real(kind_real)   :: resproc_max(numprocs)
    real(kind_real)   :: resproc_sum(numprocs)
    real(kind_real)   :: restot_sum,resdts_sum
#endif

#ifdef OMP_IMP
    real(kind_real)   :: rmx(0:omp_max_num_threads)
    integer(kind_int) :: mypos(0:omp_max_num_threads,5)
    integer(kind_int) :: idt,numt

    rmx(:) = zero

!$OMP parallel
!$OMP master
    numt = omp_get_num_threads()
    if (numt > omp_max_num_threads) then
        write(*,*) "OpenMP num_threads > omp_max_num_threads, make sure you want do it, check size of nflag"
    end if
!$OMP end master
!$OMP end parallel
#endif

    respos(:) = 1
    resmax = zero
    restot = zero
    resdts = zero
    do nc=1,nblkcoms
        nb  =  blkcoms(nc)%nb
        top => blkcoms(nc)%top

        ni = top%nijk(1)
        nj = top%nijk(2)
        nk = top%nijk(3)

        dt  => mb_dt(nb)%fld
        vol => mb_vol(nb)%fld
        dq  => mb_dq(nb)%fld

        if (nsw /= ntimeadv_steady) then
            qc => mb_qc(nb)%fld
            qm => mb_qm(nb)%fld
        end if

        if ( nscmp > nscmp_non ) then
            nRes = neqn - 1
        else
            nRes = neqn
        end if
        
        do m=1,nRes

!$OMP parallel private(k,j,i,vol0,odtvol,dq0,dq_dts,rdq_dts,dt0,rdq,rdq2,idt),reduction(+:resdts,restot)

#ifdef OMP_IMP
        idt  = omp_get_thread_num()
#endif

!$OMP do
        do k=1,nk
        do j=1,nj
        do i=1,ni
            vol0 = vol(1)%r3d(i,j,k)

            if (nsw /= ntimeadv_steady) then
                odtvol = one/(dtime + small)

                dq0 = abs(qc(m)%r3d(i,j,k) - qm(m)%r3d(i,j,k))

                dq_dts  = abs(dq(m)%r3d(i,j,k))
                rdq_dts = dq_dts*odtvol
                resdts  = resdts + rdq_dts*rdq_dts
            else
                dt0  = dt(1)%r3d(i,j,k)
                dq0  = abs(dq(m)%r3d(i,j,k))
                odtvol = one/(dt0*vol0 + small)
            end if

            rdq = dq0*odtvol
            rdq2 = rdq*rdq

#ifdef OMP_IMP
            if (rdq > rmx(idt)) then
                mypos(idt,:) = (/nb,i,j,k,m/)
                rmx(idt) = rdq
            end if
#else
            if (rdq > resmax) then
                respos(:) = (/nb,i,j,k,m/)
                resmax = rdq
            end if
#endif
            restot = restot + rdq2
        end do
        end do
        end do
!$OMP end do

!$OMP end parallel

        end do
    end do

#ifdef OMP_IMP
    do idt = 0, numt-1
        if (rmx(idt) > resmax) then
            resmax = rmx(idt)
            respos(:) = mypos(idt,:)
        end if
    end do
#endif

#ifdef PARALLEL
    resproc_max(:) = zero
    resproc_max(myid+1) = resmax
    call MPI_REDUCE(resproc_max,resproc_sum,numprocs,kind_real_mpi,MPI_SUM,master,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(resproc_sum,numprocs,kind_real_mpi,master,MPI_COMM_WORLD,ierr)
    id_max = maxloc(resproc_sum,dim=1)
    resmax = resproc_sum(id_max)
    call MPI_BCAST(respos,5,kind_int_mpi,id_max-1,MPI_COMM_WORLD,ierr)

    call MPI_REDUCE(restot,restot_sum,1,kind_real_mpi,MPI_SUM,master,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(restot_sum,1,kind_real_mpi,master,MPI_COMM_WORLD,ierr)
    restot = restot_sum

    if (nsw /= ntimeadv_steady) then
        call MPI_REDUCE(resdts,resdts_sum,1,kind_real_mpi,MPI_SUM,master,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(resdts_sum,1,kind_real_mpi,master,MPI_COMM_WORLD,ierr)
        resdts = resdts_sum
    end if
#endif

end subroutine residual

subroutine residual_sp(nsw)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : small,zero,one
    use mod_constants, only : ntimeadv_steady
    use mod_datatypes, only : top_block_t,fld_array_t
    use mod_variables, only : dtime,respos
    use mod_variables, only : resmax,restot,resdts
    use mod_fieldvars, only : nblkcoms,blkcomssp,neqn
    use mod_fieldvars, only : mb_dt,mb_volsp,mb_dq
    use mod_fieldvars, only : mb_qc,mb_qm
    use mod_constants, only : nscmp_non
    use mod_variables, only : nscmp
    use mod_parallels
    use mod_openmp
    implicit none
    integer(kind_int), intent(in) :: nsw
    integer(kind_int)          :: nc,nb,i,j,k,m,ierr
    integer(kind_int)          :: ni,nj,nk
    integer(kind_int)          :: nRes
    type(top_block_t), pointer :: top
    real(kind_real)            :: dt0,vol0,dq0
    real(kind_real)            :: odtvol,rdq,rdq2
    real(kind_real)            :: dq_dts,rdq_dts
    type(fld_array_t), pointer :: dt(:),vol(:),dq(:)
    type(fld_array_t), pointer :: qc(:),qm(:)
#ifdef PARALLEL
    integer(kind_int) :: id_max
    real(kind_real)   :: resproc_max(numprocs)
    real(kind_real)   :: resproc_sum(numprocs)
    real(kind_real)   :: restot_sum,resdts_sum
#endif

#ifdef OMP_IMP
    real(kind_real)   :: rmx(0:omp_max_num_threads)
    integer(kind_int) :: mypos(0:omp_max_num_threads,5)
    integer(kind_int) :: idt,numt

    rmx(:) = zero

!$OMP parallel
!$OMP master
    numt = omp_get_num_threads()
    if (numt > omp_max_num_threads) then
        write(*,*) "OpenMP num_threads > omp_max_num_threads, make sure you want do it, check size of nflag"
    end if
!$OMP end master
!$OMP end parallel
#endif

    respos(:) = 1
    resmax = zero
    restot = zero
    resdts = zero
    do nc=1,nblkcoms
        nb  =  blkcomssp(nc)%nb
        top => blkcomssp(nc)%top

        ni = top%nijk(1)
        nj = top%nijk(2)
        nk = top%nijk(3)

        dt  => mb_dt(nb)%fld
        vol => mb_volsp(nb)%fld
        dq  => mb_dq(nb)%fld

        if (nsw /= ntimeadv_steady) then
            qc => mb_qc(nb)%fld
            qm => mb_qm(nb)%fld
        end if

        if ( nscmp > nscmp_non ) then
            nRes = neqn - 1
        else
            nRes = neqn
        end if
        
        do m=1,nRes

!$OMP parallel private(k,j,i,vol0,odtvol,dq0,dq_dts,rdq_dts,dt0,rdq,rdq2,idt),reduction(+:resdts,restot)

#ifdef OMP_IMP
        idt  = omp_get_thread_num()
#endif

!$OMP do
        do k=1,nk
        do j=1,nj
        do i=1,ni
            vol0 = vol(1)%r3d(i,j,k)

            if (nsw /= ntimeadv_steady) then
                odtvol = one/(dtime + small)

                dq0 = abs(qc(m)%r3d(i,j,k) - qm(m)%r3d(i,j,k))

                dq_dts  = abs(dq(m)%r3d(i,j,k))
                rdq_dts = dq_dts*odtvol
                resdts  = resdts + rdq_dts*rdq_dts
            else
                dt0  = dt(1)%r3d(i,j,k)
                dq0  = abs(dq(m)%r3d(i,j,k))
                odtvol = one/(dt0*vol0 + small)
            end if

            rdq = dq0*odtvol
            rdq2 = rdq*rdq

#ifdef OMP_IMP
            if (rdq > rmx(idt)) then
                mypos(idt,:) = (/nb,i,j,k,m/)
                rmx(idt) = rdq
            end if
#else
            if (rdq > resmax) then
                respos(:) = (/nb,i,j,k,m/)
                resmax = rdq
            end if
#endif
            restot = restot + rdq2
        end do
        end do
        end do
!$OMP end do

!$OMP end parallel

        end do
    end do

#ifdef OMP_IMP
    do idt = 0, numt-1
        if (rmx(idt) > resmax) then
            resmax = rmx(idt)
            respos(:) = mypos(idt,:)
        end if
    end do
#endif

#ifdef PARALLEL
    resproc_max(:) = zero
    resproc_max(myid+1) = resmax
    call MPI_REDUCE(resproc_max,resproc_sum,numprocs,kind_real_mpi,MPI_SUM,master,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(resproc_sum,numprocs,kind_real_mpi,master,MPI_COMM_WORLD,ierr)
    id_max = maxloc(resproc_sum,dim=1)
    resmax = resproc_sum(id_max)
    call MPI_BCAST(respos,5,kind_int_mpi,id_max-1,MPI_COMM_WORLD,ierr)

    call MPI_REDUCE(restot,restot_sum,1,kind_real_mpi,MPI_SUM,master,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(restot_sum,1,kind_real_mpi,master,MPI_COMM_WORLD,ierr)
    restot = restot_sum

    if (nsw /= ntimeadv_steady) then
        call MPI_REDUCE(resdts,resdts_sum,1,kind_real_mpi,MPI_SUM,master,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(resdts_sum,1,kind_real_mpi,master,MPI_COMM_WORLD,ierr)
        resdts = resdts_sum
    end if
#endif

end subroutine residual_sp

subroutine aeroforce
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,fourth,one
    use mod_constants, only : two,two3rd,mcyc,m3x3
    use mod_constants, only : sml_ssf,bc_wall,nvis_euler
    use mod_datatypes, only : var_block_t,fld_array_t
    use mod_datatypes, only : top_block_t,bc_region_t
    use mod_variables, only : gamma,reue,refbeta,poo,prlam,prtur
    use mod_variables, only : nvis,cfx,cfy,cfz,cmx,cmy,cmz
    use mod_variables, only : frefxc,frefyc,frefzc,reflgrd
    use mod_fieldvars, only : nblkcoms,blkcoms,mb_xyz,mb_sxyz
    use mod_fieldvars, only : mb_pv,mb_dpv,mb_vsl,mb_vst
    use mod_parallels
    use mod_tecplotio
    implicit none
    integer(kind_int)          :: nc,nb,nr,i,j,k,m,n,m1,m2,m3,ierr
    integer(kind_int)          :: s_st(3),s_ed(3),st(3),ed(3),ijk(3,4)
    integer(kind_int)          :: nregs,s_nd,s_lr,bctype,nd2,nd3,is,js,ks
    real(kind_real)            :: vis,kcp,clr,clrre,nx,ny,nz,sn,osn
    real(kind_real)            :: dux,duy,duz,dvx,dvy,dvz,dwx,dwy,dwz
    real(kind_real)            :: tauxx,tauyy,tauzz,tauxy,tauxz,tauyz
    real(kind_real)            :: re,vis2p3,dfx,dfy,dfz,xc,yc,zc
    real(kind_real)            :: fp,fvx,fvy,fvz,dmx,dmy,dmz,xw,yw,zw
    real(kind_real)            :: rxyz(3,4),r13(3),r24(3)
    real(kind_real),   pointer :: dfxyz(:,:,:,:)
    real(kind_real),  external :: area_cross
    type(top_block_t), pointer :: top
    type(bc_region_t), pointer :: reg
    type(fld_array_t), pointer :: xyz(:),sxyz(:),pv(:)
    type(fld_array_t), pointer :: dpv(:),vsl(:),vst(:)

    re = one/reue

    xc = frefxc/reflgrd
    yc = frefyc/reflgrd
    zc = frefzc/reflgrd

    cfx = zero
    cfy = zero
    cfz = zero
    cmx = zero
    cmy = zero
    cmz = zero
    do nc=1,nblkcoms
        nb  =  blkcoms(nc)%nb
        top => blkcoms(nc)%top

        xyz  => mb_xyz(nb)%fld
        sxyz => mb_sxyz(nb)%fld
        pv   => mb_pv(nb)%fld
        if (nvis > nvis_euler) then
            vsl  => mb_vsl(nb)%fld
            vst  => mb_vst(nb)%fld
            dpv  => mb_dpv(nb)%fld
        end if

        nregs = top%nregions
        do nr=1,nregs
            reg => top%bcs(nr)

            bctype  = reg%bctype
            if (bctype == bc_wall) then
                s_st(:) = reg%s_st(:)
                s_ed(:) = reg%s_ed(:)
                s_nd    = reg%s_nd
                s_lr    = reg%s_lr

                st(:) = s_st(:)
                ed(:) = s_ed(:)

                m1 = m3x3(1,s_nd)
                m2 = m3x3(2,s_nd)
                m3 = m3x3(3,s_nd)

                clr   = two*s_lr
                clrre = -clr*re

                allocate(dfxyz(1:3,st(1):ed(1),st(2):ed(2),st(3):ed(3)), stat=ierr)

                do k=st(3),ed(3)
                do j=st(2),ed(2)
                do i=st(1),ed(1)
                    nx = sxyz(m1)%r3d(i,j,k)
                    ny = sxyz(m2)%r3d(i,j,k)
                    nz = sxyz(m3)%r3d(i,j,k)
                    sn = max(sqrt(nx*nx + ny*ny + nz*nz),sml_ssf)
                    osn = one/sn
                    nx = nx*osn
                    ny = ny*osn
                    nz = nz*osn

                    fp  = clr*(pv(5)%r3d(i,j,k) - poo)
                    dfx = fp*nx
                    dfy = fp*ny
                    dfz = fp*nz

                    if (nvis > nvis_euler) then
                        dux = dpv(1)%r3d(i,j,k)
                        duy = dpv(2)%r3d(i,j,k)
                        duz = dpv(3)%r3d(i,j,k)

                        dvx = dpv(4)%r3d(i,j,k)
                        dvy = dpv(5)%r3d(i,j,k)
                        dvz = dpv(6)%r3d(i,j,k)

                        dwx = dpv(7)%r3d(i,j,k)
                        dwy = dpv(8)%r3d(i,j,k)
                        dwz = dpv(9)%r3d(i,j,k)

                        vis = vsl(1)%r3d(i,j,k)
                        vis = vis + vst(1)%r3d(i,j,k)

                        vis2p3 = two3rd*vis
                        tauxx = vis2p3 * ( two*dux - dvy - dwz )
                        tauyy = vis2p3 * ( two*dvy - dwz - dux )
                        tauzz = vis2p3 * ( two*dwz - dux - dvy )
                        tauxy = vis * ( duy + dvx )
                        tauxz = vis * ( duz + dwx )
                        tauyz = vis * ( dvz + dwy )

                        fvx = clrre*(nx*tauxx + ny*tauxy + nz*tauxz)
                        fvy = clrre*(nx*tauxy + ny*tauyy + nz*tauyz)
                        fvz = clrre*(nx*tauxz + ny*tauyz + nz*tauzz)
                        dfx = dfx + fvx
                        dfy = dfy + fvy
                        dfz = dfz + fvz
                    end if
                    dfxyz(1,i,j,k) = dfx
                    dfxyz(2,i,j,k) = dfy
                    dfxyz(3,i,j,k) = dfz
                end do
                end do
                end do

                nd2 = mcyc(2,s_nd)
                nd3 = mcyc(3,s_nd)
                ed(nd2) = ed(nd2) - 1
                ed(nd3) = ed(nd3) - 1
                do k=st(3),ed(3)
                do j=st(2),ed(2)
                do i=st(1),ed(1)
                    ijk(:,1) = (/i,j,k/)
                    ijk(:,2) = ijk(:,1)
                    ijk(nd2,2) = ijk(nd2,2) + 1

                    ijk(:,3) = ijk(:,2)
                    ijk(nd3,3) = ijk(nd3,3) + 1

                    ijk(:,4) = ijk(:,1)
                    ijk(nd3,4) = ijk(nd3,4) + 1

                    do n=1,4
                        is = ijk(1,n)
                        js = ijk(2,n)
                        ks = ijk(3,n)
                        do m=1,3
                            rxyz(m,n) = xyz(m)%r3d(is,js,ks)
                        end do
                    end do

                    r13(:) = rxyz(:,1)-rxyz(:,3)
                    r24(:) = rxyz(:,4)-rxyz(:,2)
                    sn = area_cross(r13, r24)
                    sn = fourth*sn

                    do n=1,4
                        is = ijk(1,n)
                        js = ijk(2,n)
                        ks = ijk(3,n)
                        dfx = dfxyz(1,is,js,ks)*sn
                        dfy = dfxyz(2,is,js,ks)*sn
                        dfz = dfxyz(3,is,js,ks)*sn
                        cfx = cfx + dfx
                        cfy = cfy + dfy
                        cfz = cfz + dfz

                        xw = rxyz(1,n)
                        yw = rxyz(2,n)
                        zw = rxyz(3,n)
                        dmx = (yw - yc)*dfz - (zw - zc)*dfy
                        dmy = (zw - zc)*dfx - (xw - xc)*dfz
                        dmz = (xw - xc)*dfy - (yw - yc)*dfx
                        cmx = cmx + dmx
                        cmy = cmy + dmy
                        cmz = cmz + dmz
                    end do
                end do
                end do
                end do

                deallocate(dfxyz, stat=ierr)
            end if
        end do
    end do

end subroutine aeroforce

function area_cross(a,b) result(s)
    use mod_kndconsts, only : kind_real
    use mod_constants, only : half
    implicit none
    real(kind_real) :: a(3),b(3)
    real(kind_real) :: s
    real(kind_real) :: c(3)

    c(1) = a(2) * b(3) - a(3) * b(2)
    c(2) = a(3) * b(1) - a(1) * b(3)
    c(3) = a(1) * b(2) - a(2) * b(1)

    s = half*sqrt(sum(c(:)*c(:)))

end function area_cross


subroutine aeroforce5
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,half,one,two,two3rd,m3x3
    use mod_constants, only : sml_ssf,bc_wall,nvis_euler
    use mod_datatypes, only : var_block_t,fld_array_t
    use mod_datatypes, only : top_block_t,bc_region_t
    use mod_variables, only : gamma,reue,refbeta,poo,prlam,prtur
    use mod_variables, only : nvis,cfx,cfy,cfz,cmx,cmy,cmz
    use mod_variables, only : frefxc,frefyc,frefzc,reflgrd
    use mod_fieldvars, only : nblkcoms,blkcoms,mb_xyz,mb_sxyz
    use mod_fieldvars, only : mb_pv,mb_dpv,mb_vsl,mb_vst
    use mod_parallels
    use mod_tecplotio
    implicit none
    integer(kind_int)          :: nc,nb,nr,i,j,k,m,m1,m2,m3,ierr
    integer(kind_int)          :: s_st(3),s_ed(3),st(3),ed(3)
    integer(kind_int)          :: nregs,s_nd,s_lr,bctype,ijk(3)
    real(kind_real)            :: vis,kcp,clr,clrre,nx,ny,nz,sn,osn
    real(kind_real)            :: dux,duy,duz,dvx,dvy,dvz,dwx,dwy,dwz
    real(kind_real)            :: tauxx,tauyy,tauzz,tauxy,tauxz,tauyz
    real(kind_real)            :: re,vis2p3,dfx,dfy,dfz,xc,yc,zc,acsn
    real(kind_real)            :: fp,fvx,fvy,fvz,dmx,dmy,dmz,xw,yw,zw
    real(kind_real),  external :: area_coeff
    type(top_block_t), pointer :: top
    type(bc_region_t), pointer :: reg
    type(fld_array_t), pointer :: xyz(:),sxyz(:),pv(:)
    type(fld_array_t), pointer :: dpv(:),vsl(:),vst(:)

    re = one/reue

    xc = frefxc/reflgrd
    yc = frefyc/reflgrd
    zc = frefzc/reflgrd

    cfx = zero
    cfy = zero
    cfz = zero
    cmx = zero
    cmy = zero
    cmz = zero
    do nc=1,nblkcoms
        nb  =  blkcoms(nc)%nb
        top => blkcoms(nc)%top

        xyz  => mb_xyz(nb)%fld
        sxyz => mb_sxyz(nb)%fld
        pv   => mb_pv(nb)%fld
        if (nvis > nvis_euler) then
            vsl => mb_vsl(nb)%fld
            vst => mb_vst(nb)%fld
            dpv => mb_dpv(nb)%fld
        end if

        nregs = top%nregions
        do nr=1,nregs
            reg => top%bcs(nr)

            bctype  = reg%bctype
            if (bctype == bc_wall) then
                s_st(:) = reg%s_st(:)
                s_ed(:) = reg%s_ed(:)
                s_nd    = reg%s_nd
                s_lr    = reg%s_lr

                st(:) = s_st(:)
                ed(:) = s_ed(:)

                m1 = m3x3(1,s_nd)
                m2 = m3x3(2,s_nd)
                m3 = m3x3(3,s_nd)

                clr   = two*s_lr
                clrre = -clr*re

                do k=st(3),ed(3)
                do j=st(2),ed(2)
                do i=st(1),ed(1)
                    nx = sxyz(m1)%r3d(i,j,k)
                    ny = sxyz(m2)%r3d(i,j,k)
                    nz = sxyz(m3)%r3d(i,j,k)
                    sn = max(sqrt(nx*nx + ny*ny + nz*nz),sml_ssf)
                    osn = one/sn
                    nx = nx*osn
                    ny = ny*osn
                    nz = nz*osn

                    fp  = clr*(pv(5)%r3d(i,j,k) - poo)
                    dfx = fp*nx
                    dfy = fp*ny
                    dfz = fp*nz

                    if (nvis > nvis_euler) then
                        dux = dpv(1)%r3d(i,j,k)
                        duy = dpv(2)%r3d(i,j,k)
                        duz = dpv(3)%r3d(i,j,k)

                        dvx = dpv(4)%r3d(i,j,k)
                        dvy = dpv(5)%r3d(i,j,k)
                        dvz = dpv(6)%r3d(i,j,k)

                        dwx = dpv(7)%r3d(i,j,k)
                        dwy = dpv(8)%r3d(i,j,k)
                        dwz = dpv(9)%r3d(i,j,k)

                        vis = vsl(1)%r3d(i,j,k)
                        vis = vis + vst(1)%r3d(i,j,k)

                        vis2p3 = two3rd*vis
                        tauxx = vis2p3 * ( two*dux - dvy - dwz )
                        tauyy = vis2p3 * ( two*dvy - dwz - dux )
                        tauzz = vis2p3 * ( two*dwz - dux - dvy )
                        tauxy = vis * ( duy + dvx )
                        tauxz = vis * ( duz + dwx )
                        tauyz = vis * ( dvz + dwy )

                        fvx = clrre*(nx*tauxx + ny*tauxy + nz*tauxz)
                        fvy = clrre*(nx*tauxy + ny*tauyy + nz*tauyz)
                        fvz = clrre*(nx*tauxz + ny*tauyz + nz*tauzz)
                        dfx = dfx + fvx
                        dfy = dfy + fvy
                        dfz = dfz + fvz
                    end if
                    ijk(:) = (/i,j,k/)
                    acsn = sn*area_coeff(ijk,st,ed,s_nd)
                    dfx = dfx*acsn
                    dfy = dfy*acsn
                    dfz = dfz*acsn

                    cfx = cfx + dfx
                    cfy = cfy + dfy
                    cfz = cfz + dfz

                    xw = xyz(1)%r3d(i,j,k)
                    yw = xyz(2)%r3d(i,j,k)
                    zw = xyz(3)%r3d(i,j,k)
                    dmx = (yw - yc)*dfz - (zw - zc)*dfy
                    dmy = (zw - zc)*dfx - (xw - xc)*dfz
                    dmz = (xw - xc)*dfy - (yw - yc)*dfx
                    cmx = cmx + dmx
                    cmy = cmy + dmy
                    cmz = cmz + dmz
                end do
                end do
                end do
            end if
        end do
    end do

end subroutine aeroforce5


function area_coeff(ijk,st,ed,nd) result(ac)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : half,fourth,one,mcyc
    implicit none
    integer(kind_int) :: ijk(3),st(3),ed(3),nd
    real(kind_real)   :: ac
    integer(kind_int) :: m,n,d,n0

    n0 = 0
    do m=2,3
       n = mcyc(m,nd)
       d = min(abs(ijk(n) - st(n)), abs(ijk(n) - ed(n)) )
       if (d == 0) n0 = n0 + 1
    end do

    if (n0 == 1) then
       ac = half
    else if (n0 == 2) then
       ac = fourth
    else
       ac = one
    end if

end function area_coeff

!!subroutine flow_mean
!!    use mod_kndconsts, only : kind_int,kind_real
!!    use mod_datatypes, only : fld_array_t,top_block_t
!!    use mod_fieldvars, only : npvs
!!    use mod_fieldvars, only : nblkcoms,blkcoms
!!    use mod_fieldvars, only : mb_pv,mb_fmean,nstepmean
!!    implicit none
!!    integer(kind_int)          :: nc,nb,i,j,k,ierr
!!    integer(kind_int)          :: m,ni,nj,nk
!!    type(top_block_t), pointer :: top
!!    real(kind_real)            :: temp
!!    type(fld_array_t), pointer :: pv(:),fmean(:)
!!
!!    nstepmean = nstepmean + 1
!!
!!    do nc=1,nblkcoms
!!        nb  =  blkcoms(nc)%nb
!!        top => blkcoms(nc)%top
!!
!!        ni = top%nijk(1)
!!        nj = top%nijk(2)
!!        nk = top%nijk(3)
!!
!!        pv    => mb_pv(nb)%fld
!!        fmean => mb_fmean(nb)%fld
!!
!!        do k=1,nk
!!        do j=1,nj
!!        do i=1,ni
!!        do m=1,npvs
!!            temp = pv(m)%r3d(i,j,k)
!!            fmean(m)%r3d(i,j,k) = (nstepmean-1)*fmean(m)%r3d(i,j,k) + temp
!!        end do
!!        end do
!!        end do
!!        end do
!!
!!        do k=2,nk
!!        do j=1,nj
!!        do i=1,ni
!!        do m=1,npvs
!!            fmean(m)%r3d(i,j,1) = fmean(m)%r3d(i,j,1) + fmean(m)%r3d(i,j,k)
!!        end do
!!        end do
!!        end do
!!        end do
!!        
!!        do k=2,nk
!!        do j=1,nj
!!        do i=1,ni
!!        do m=1,npvs
!!            fmean(m)%r3d(i,j,k) = fmean(m)%r3d(i,j,1)
!!        end do
!!        end do
!!        end do
!!        end do        
!!
!!    end do
!!
!!end subroutine flow_mean
!!
!!
!!subroutine flow_rms
!!    use mod_kndconsts, only : kind_int,kind_real
!!    use mod_datatypes, only : fld_array_t,top_block_t
!!    use mod_fieldvars, only : npvs
!!    use mod_fieldvars, only : nblkcoms,blkcoms
!!    use mod_fieldvars, only : mb_pv,nsteprms
!!    use mod_fieldvars, only : mb_frms,mb_fmean
!!    implicit none
!!    integer(kind_int)          :: nc,nb,i,j,k,ierr
!!    integer(kind_int)          :: ni,nj,nk,m
!!    type(top_block_t), pointer :: top
!!    real(kind_real)            :: temp,mean,flow
!!    type(fld_array_t), pointer :: pv(:),frms(:),fmean(:)
!!
!!    nsteprms = nsteprms + 1
!!
!!    do nc=1,nblkcoms
!!        nb  =  blkcoms(nc)%nb
!!        top => blkcoms(nc)%top
!!
!!        ni = top%nijk(1)
!!        nj = top%nijk(2)
!!        nk = top%nijk(3)
!!
!!        pv    => mb_pv(nb)%fld
!!        fmean => mb_fmean(nb)%fld
!!        frms  => mb_frms(nb)%fld
!!
!!        do k=1,nk
!!        do j=1,nj
!!        do i=1,ni
!!        do m=1,npvs
!!            flow = pv(m)%r3d(i,j,k)
!!            mean = fmean(m)%r3d(i,j,k)
!!            temp = flow - mean
!!            frms(m)%r3d(i,j,k) = (nsteprms-1)*frms(m)%r3d(i,j,k)**2.0 + temp*temp
!!        end do
!!        end do
!!        end do
!!        end do
!!
!!        do k=2,nk
!!        do j=1,nj
!!        do i=1,ni
!!        do m=1,npvs
!!            frms(m)%r3d(i,j,1) = frms(m)%r3d(i,j,1) + frms(m)%r3d(i,j,k)
!!        end do
!!        end do
!!        end do
!!        end do
!!        
!!        do k=2,nk
!!        do j=1,nj
!!        do i=1,ni
!!        do m=1,npvs
!!            frms(m)%r3d(i,j,k) = frms(m)%r3d(i,j,1)
!!        end do
!!        end do
!!        end do
!!        end do        
!!
!!    end do
!!
!!end subroutine flow_rms
!!
!!subroutine solve_fmean
!!    use mod_kndconsts, only : kind_int,kind_real
!!    use mod_datatypes, only : fld_array_t,top_block_t
!!    use mod_fieldvars, only : npvs
!!    use mod_fieldvars, only : nblkcoms,blkcoms
!!    use mod_fieldvars, only : mb_fmean,nstepmean
!!    implicit none
!!    integer(kind_int)          :: nc,nb,i,j,k,ierr
!!    integer(kind_int)          :: m,ni,nj,nk
!!    type(top_block_t), pointer :: top
!!    type(fld_array_t), pointer :: fmean(:)
!!
!!
!!    do nc=1,nblkcoms
!!        nb  =  blkcoms(nc)%nb
!!        top => blkcoms(nc)%top
!!
!!        ni = top%nijk(1)
!!        nj = top%nijk(2)
!!        nk = top%nijk(3)
!!
!!        fmean => mb_fmean(nb)%fld
!!
!!        do k=1,nk
!!        do j=1,nj
!!        do i=1,ni
!!        do m=1,npvs
!!            fmean(m)%r3d(i,j,k) = fmean(m)%r3d(i,j,k)/(nstepmean*1.0)/(nk*1.0)
!!        end do
!!        end do
!!        end do
!!        end do
!!    end do
!!
!!end subroutine solve_fmean
!!
!!subroutine solve_frms
!!    use mod_kndconsts, only : kind_int,kind_real
!!    use mod_datatypes, only : fld_array_t,top_block_t
!!    use mod_fieldvars, only : npvs
!!    use mod_fieldvars, only : nblkcoms,blkcoms
!!    use mod_fieldvars, only : nsteprms
!!    use mod_fieldvars, only : mb_frms
!!    implicit none
!!    integer(kind_int)          :: nc,nb,i,j,k,ierr
!!    integer(kind_int)          :: ni,nj,nk,m
!!    type(top_block_t), pointer :: top
!!    type(fld_array_t), pointer :: frms(:)
!!
!!    do nc=1,nblkcoms
!!        nb  =  blkcoms(nc)%nb
!!        top => blkcoms(nc)%top
!!
!!        ni = top%nijk(1)
!!        nj = top%nijk(2)
!!        nk = top%nijk(3)
!!
!!        frms  => mb_frms(nb)%fld
!!
!!        do k=1,nk
!!        do j=1,nj
!!        do i=1,ni
!!        do m=1,npvs
!!            frms(m)%r3d(i,j,k) = sqrt( frms(m)%r3d(i,j,k)/(nsteprms*1.0)/(nk*1.0) )
!!        end do
!!        end do
!!        end do
!!        end do
!!    end do
!!
!!end subroutine solve_frms

subroutine flow_mean
    use mod_kndconsts, only : kind_int,kind_real
    use mod_datatypes, only : fld_array_t,top_block_t
    use mod_fieldvars, only : npvs
    use mod_fieldvars, only : nblkcoms,blkcoms
    use mod_fieldvars, only : mb_pv,mb_fmean,nstepmean
    implicit none
    integer(kind_int)          :: nc,nb,i,j,k,ierr
    integer(kind_int)          :: m,ni,nj,nk
    type(top_block_t), pointer :: top
    real(kind_real)            :: temp
    type(fld_array_t), pointer :: pv(:),fmean(:)

    nstepmean = nstepmean + 1

    do nc=1,nblkcoms
        nb  =  blkcoms(nc)%nb
        top => blkcoms(nc)%top

        ni = top%nijk(1)
        nj = top%nijk(2)
        nk = top%nijk(3)

        pv    => mb_pv(nb)%fld
        fmean => mb_fmean(nb)%fld

        do k=1,nk
        do j=1,nj
        do i=1,ni
        do m=1,npvs
            temp = pv(m)%r3d(i,j,k)
            fmean(m)%r3d(i,j,k) = (nstepmean-1)*fmean(m)%r3d(i,j,k) + temp
        end do
        end do
        end do
        end do

    end do

end subroutine flow_mean

subroutine flow_mean_sp
    use mod_kndconsts, only : kind_int,kind_real
    use mod_datatypes, only : fld_array_t,top_block_t
    use mod_fieldvars, only : npvs
    use mod_fieldvars, only : nblkcoms,blkcomssp
    use mod_fieldvars, only : mb_pv,mb_fmean,nstepmean
    implicit none
    integer(kind_int)          :: nc,nb,i,j,k,ierr
    integer(kind_int)          :: m,ni,nj,nk
    type(top_block_t), pointer :: top
    real(kind_real)            :: temp
    type(fld_array_t), pointer :: pv(:),fmean(:)

    nstepmean = nstepmean + 1

    do nc=1,nblkcoms
        nb  =  blkcomssp(nc)%nb
        top => blkcomssp(nc)%top

        ni = top%nijk(1)
        nj = top%nijk(2)
        nk = top%nijk(3)

        pv    => mb_pv(nb)%fld
        fmean => mb_fmean(nb)%fld

        do k=1,nk
        do j=1,nj
        do i=1,ni
        do m=1,npvs
            temp = pv(m)%r3d(i,j,k)
            fmean(m)%r3d(i,j,k) = (nstepmean-1)*fmean(m)%r3d(i,j,k) + temp
        end do
        end do
        end do
        end do

    end do

end subroutine flow_mean_sp

subroutine flow_rms
    use mod_kndconsts, only : kind_int,kind_real
    use mod_datatypes, only : fld_array_t,top_block_t
    use mod_fieldvars, only : npvs
    use mod_fieldvars, only : nblkcoms,blkcoms
    use mod_fieldvars, only : mb_pv,nsteprms
    use mod_fieldvars, only : mb_frms,mb_fmean
    implicit none
    integer(kind_int)          :: nc,nb,i,j,k,ierr
    integer(kind_int)          :: ni,nj,nk,m
    type(top_block_t), pointer :: top
    real(kind_real)            :: temp,mean,flow
    type(fld_array_t), pointer :: pv(:),frms(:),fmean(:)

    nsteprms = nsteprms + 1

    do nc=1,nblkcoms
        nb  =  blkcoms(nc)%nb
        top => blkcoms(nc)%top

        ni = top%nijk(1)
        nj = top%nijk(2)
        nk = top%nijk(3)

        pv    => mb_pv(nb)%fld
        fmean => mb_fmean(nb)%fld
        frms  => mb_frms(nb)%fld

        do k=1,nk
        do j=1,nj
        do i=1,ni
        do m=1,npvs
            flow = pv(m)%r3d(i,j,k)
            mean = fmean(m)%r3d(i,j,k)
            temp = flow - mean
            frms(m)%r3d(i,j,k) = (nsteprms-1)*frms(m)%r3d(i,j,k)**2.0 + temp*temp
        end do
        end do
        end do
        end do

    end do

end subroutine flow_rms

subroutine flow_rms_sp
    use mod_kndconsts, only : kind_int,kind_real
    use mod_datatypes, only : fld_array_t,top_block_t
    use mod_fieldvars, only : npvs
    use mod_fieldvars, only : nblkcoms,blkcomssp
    use mod_fieldvars, only : mb_pv,nsteprms
    use mod_fieldvars, only : mb_frms,mb_fmean
    implicit none
    integer(kind_int)          :: nc,nb,i,j,k,ierr
    integer(kind_int)          :: ni,nj,nk,m
    type(top_block_t), pointer :: top
    real(kind_real)            :: temp,mean,flow
    type(fld_array_t), pointer :: pv(:),frms(:),fmean(:)

    nsteprms = nsteprms + 1

    do nc=1,nblkcoms
        nb  =  blkcomssp(nc)%nb
        top => blkcomssp(nc)%top

        ni = top%nijk(1)
        nj = top%nijk(2)
        nk = top%nijk(3)

        pv    => mb_pv(nb)%fld
        fmean => mb_fmean(nb)%fld
        frms  => mb_frms(nb)%fld

        do k=1,nk
        do j=1,nj
        do i=1,ni
        do m=1,npvs
            flow = pv(m)%r3d(i,j,k)
            mean = fmean(m)%r3d(i,j,k)
            temp = flow - mean
            frms(m)%r3d(i,j,k) = (nsteprms-1)*frms(m)%r3d(i,j,k)**2.0 + temp*temp
        end do
        end do
        end do
        end do

    end do

end subroutine flow_rms_sp

subroutine solve_fmean
    use mod_kndconsts, only : kind_int,kind_real
    use mod_datatypes, only : fld_array_t,top_block_t
    use mod_fieldvars, only : npvs
    use mod_fieldvars, only : nblkcoms,blkcoms
    use mod_fieldvars, only : mb_fmean,nstepmean
    implicit none
    integer(kind_int)          :: nc,nb,i,j,k,ierr
    integer(kind_int)          :: m,ni,nj,nk
    type(top_block_t), pointer :: top
    type(fld_array_t), pointer :: fmean(:)


    do nc=1,nblkcoms
        nb  =  blkcoms(nc)%nb
        top => blkcoms(nc)%top

        ni = top%nijk(1)
        nj = top%nijk(2)
        nk = top%nijk(3)

        fmean => mb_fmean(nb)%fld

        do k=1,nk
        do j=1,nj
        do i=1,ni
        do m=1,npvs
            fmean(m)%r3d(i,j,k) = fmean(m)%r3d(i,j,k)/(nstepmean*1.0)
        end do
        end do
        end do
        end do
    end do

end subroutine solve_fmean

subroutine solve_fmean_sp
    use mod_kndconsts, only : kind_int,kind_real
    use mod_datatypes, only : fld_array_t,top_block_t
    use mod_fieldvars, only : npvs
    use mod_fieldvars, only : nblkcoms,blkcomssp
    use mod_fieldvars, only : mb_fmean,nstepmean
    implicit none
    integer(kind_int)          :: nc,nb,i,j,k,ierr
    integer(kind_int)          :: m,ni,nj,nk
    type(top_block_t), pointer :: top
    type(fld_array_t), pointer :: fmean(:)


    do nc=1,nblkcoms
        nb  =  blkcomssp(nc)%nb
        top => blkcomssp(nc)%top

        ni = top%nijk(1)
        nj = top%nijk(2)
        nk = top%nijk(3)

        fmean => mb_fmean(nb)%fld

        do k=1,nk
        do j=1,nj
        do i=1,ni
        do m=1,npvs
            fmean(m)%r3d(i,j,k) = fmean(m)%r3d(i,j,k)/(nstepmean*1.0)
        end do
        end do
        end do
        end do
    end do

end subroutine solve_fmean_sp

subroutine solve_frms
    use mod_kndconsts, only : kind_int,kind_real
    use mod_datatypes, only : fld_array_t,top_block_t
    use mod_fieldvars, only : npvs
    use mod_fieldvars, only : nblkcoms,blkcoms
    use mod_fieldvars, only : nsteprms
    use mod_fieldvars, only : mb_frms
    implicit none
    integer(kind_int)          :: nc,nb,i,j,k,ierr
    integer(kind_int)          :: ni,nj,nk,m
    type(top_block_t), pointer :: top
    type(fld_array_t), pointer :: frms(:)

    do nc=1,nblkcoms
        nb  =  blkcoms(nc)%nb
        top => blkcoms(nc)%top

        ni = top%nijk(1)
        nj = top%nijk(2)
        nk = top%nijk(3)

        frms  => mb_frms(nb)%fld

        do k=1,nk
        do j=1,nj
        do i=1,ni
        do m=1,npvs
            frms(m)%r3d(i,j,k) = sqrt( frms(m)%r3d(i,j,k)/(nsteprms*1.0) )
        end do
        end do
        end do
        end do
    end do

end subroutine solve_frms

subroutine solve_frms_sp
    use mod_kndconsts, only : kind_int,kind_real
    use mod_datatypes, only : fld_array_t,top_block_t
    use mod_fieldvars, only : npvs
    use mod_fieldvars, only : nblkcoms,blkcomssp
    use mod_fieldvars, only : nsteprms
    use mod_fieldvars, only : mb_frms
    implicit none
    integer(kind_int)          :: nc,nb,i,j,k,ierr
    integer(kind_int)          :: ni,nj,nk,m
    type(top_block_t), pointer :: top
    type(fld_array_t), pointer :: frms(:)

    do nc=1,nblkcoms
        nb  =  blkcomssp(nc)%nb
        top => blkcomssp(nc)%top

        ni = top%nijk(1)
        nj = top%nijk(2)
        nk = top%nijk(3)

        frms  => mb_frms(nb)%fld

        do k=1,nk
        do j=1,nj
        do i=1,ni
        do m=1,npvs
            frms(m)%r3d(i,j,k) = sqrt( frms(m)%r3d(i,j,k)/(nsteprms*1.0) )
        end do
        end do
        end do
        end do
    end do

end subroutine solve_frms_sp

subroutine to_noise_write()
use mod_constants,only : io_unit_acou
use mod_fieldvars,only : npoints,ntime
use mod_variables,only : casname,nsw_kdir,nstep
use mod_variables,only : vinf,nstep,reflgrd,dtime,moo
use mod_parallels
implicit none

#ifdef PARALLEL
    if (myid == master) then
#endif

   open (io_unit_acou,file='noise/'//trim(casname)//'_size.dat',status='unknown')
   write(io_unit_acou,*) npoints,ntime,abs(1-nsw_kdir),dtime*(nstep*1.0)/moo
   close(io_unit_acou)

#ifdef PARALLEL
    end if
#endif
end subroutine to_noise_write

subroutine calculate_points_acoustic
    use mod_kndconsts, only : kind_int
    use mod_constants, only : nsw_dir_close
    use mod_datatypes, only : top_block_t,bc_region_t
    use mod_variables, only : nsw_kdir
    use mod_fieldvars, only : nblkcoms,blkcoms,npoints
    use mod_parallels
    implicit none
    integer(kind_int)          :: nc,nb,ns,nr,npoints0
    integer(kind_int)          :: nregs,s_nd,subtype_A,ierr
    type(top_block_t), pointer :: top
    type(bc_region_t), pointer :: reg
    logical :: doit

    do nc=1,nblkcoms
        nb  =  blkcoms(nc)%nb
        top => blkcoms(nc)%top

        nregs = top%nregions
        do ns=1,nregs
            nr  =  top%bcnrs(ns)
            reg => top%bcs(nr)

            s_nd      = reg%s_nd
            subtype_A = reg%subtype_A

            doit = .true.
            if (nsw_kdir == nsw_dir_close .and. s_nd == 3) then
               doit = .false.
            end if

            if (doit) then
                select case(subtype_A)
                case(1)
                    call set_points_acoustic(nb,nr)
                case default

                end select
            end if
        end do
    end do
    npoints0 = npoints
#ifdef PARALLEL
    call MPI_REDUCE(npoints,npoints0,1,kind_int_mpi,MPI_SUM,master,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(npoints0,1,kind_int_mpi,master,MPI_COMM_WORLD,ierr)
    npoints = npoints0
#endif
end subroutine calculate_points_acoustic

subroutine set_points_acoustic(nb,nr)
    use mod_kndconsts, only : kind_int
    use mod_datatypes, only : bc_region_t
    use mod_fieldvars, only : mb_top,npoints
    implicit none
    integer(kind_int), intent(in) :: nb,nr
    integer(kind_int)            :: s_st(3),s_ed(3)
    type(bc_region_t), pointer   :: reg

    reg     => mb_top(nb)%bcs(nr)
    s_st(:) = reg%s_st(:)
    s_ed(:) = reg%s_ed(:)

    npoints = npoints + ( s_ed(3) - s_st(3) + 1 )* &
                        ( s_ed(2) - s_st(2) + 1 )* &
                        ( s_ed(1) - s_st(1) + 1 )   

end subroutine set_points_acoustic

subroutine calculate_points_acoustic_sp
    use mod_kndconsts, only : kind_int
    use mod_constants, only : nsw_dir_close
    use mod_datatypes, only : top_block_t,bc_region_t
    use mod_variables, only : nsw_kdir
    use mod_fieldvars, only : nblkcoms,blkcomssp,npoints
    use mod_parallels
    implicit none
    integer(kind_int)          :: nc,nb,ns,nr,npoints0
    integer(kind_int)          :: nregs,s_nd,subtype_A,ierr
    type(top_block_t), pointer :: top
    type(bc_region_t), pointer :: reg
    logical :: doit

    do nc=1,nblkcoms
        nb  =  blkcomssp(nc)%nb
        top => blkcomssp(nc)%top

        nregs = top%nregions
        do ns=1,nregs
            nr  =  top%bcnrs(ns)
            reg => top%bcs(nr)

            s_nd      = reg%s_nd
            subtype_A = reg%subtype_A

            doit = .true.
            if (nsw_kdir == nsw_dir_close .and. s_nd == 3) then
               doit = .false.
            end if

            if (doit) then
                select case(subtype_A)
                case(1)
                    call set_points_acoustic_sp(nb,nr)
                case default

                end select
            end if
        end do
    end do
    npoints0 = npoints
#ifdef PARALLEL
    call MPI_REDUCE(npoints,npoints0,1,kind_int_mpi,MPI_SUM,master,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(npoints0,1,kind_int_mpi,master,MPI_COMM_WORLD,ierr)
    npoints = npoints0
#endif
end subroutine calculate_points_acoustic_sp

subroutine set_points_acoustic_sp(nb,nr)
    use mod_kndconsts, only : kind_int
    use mod_datatypes, only : bc_region_t
    use mod_fieldvars, only : mb_topsp,npoints
    implicit none
    integer(kind_int), intent(in) :: nb,nr
    integer(kind_int)            :: s_st(3),s_ed(3)
    type(bc_region_t), pointer   :: reg

    reg     => mb_topsp(nb)%bcs(nr)
    s_st(:) = reg%s_st(:)
    s_ed(:) = reg%s_ed(:)

    npoints = npoints + ( s_ed(3) - s_st(3) + 1 )* &
                        ( s_ed(2) - s_st(2) + 1 )* &
                        ( s_ed(1) - s_st(1) + 1 )    

end subroutine set_points_acoustic_sp

subroutine to_noise_output()
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : m3x3
    use mod_constants, only : io_unit_nois
    use mod_datatypes, only : var_block_t,fld_array_t
    use mod_datatypes, only : top_block_t,bc_region_t
    use mod_variables, only : casname,vinf,nstep,reflgrd,dtime,poo,moo
    use mod_fieldvars, only : nblocks,mb_top,mb_xyz,mb_sxyz
    use mod_fieldvars, only : mb_pv,ntime,npoints
    use mod_parallels
    use mod_tecplotio
    implicit none
    integer(kind_int)          :: nb,nr,i,j,k,m,m1,m2,m3,ierr,id
    integer(kind_int)          :: s_st(3),s_ed(3),st(3),ed(3)
    integer(kind_int)          :: nregs,s_nd,s_lr,subtype_A
    real(kind_real)            :: nx,ny,nz
    type(top_block_t), pointer :: top
    type(bc_region_t), pointer :: reg
    type(fld_array_t), pointer :: xyz(:),sxyz(:)
    type(fld_array_t), pointer :: pv(:)
    real(kind_real),   pointer :: vars(:,:,:,:)
    real(kind_real),   pointer :: tvar(:,:)
    character(len=32 )         :: str1
#ifdef PARALLEL
    integer(kind_int) :: pid,packsize,tag
    integer(kind_int) :: status(MPI_STATUS_SIZE)
#endif

#ifdef PARALLEL
    if (myid == master) then
#endif
    write(str1, *) ntime + 10000000
    str1 = adjustl(str1)
    open (io_unit_nois,file='noise/'//trim(casname)//'_'//trim(str1(2:8))//'.sfl',form='unformatted',status='unknown',iostat=ierr)
    call error_check(ierr,trim(adjustl('noise/'//trim(casname)//'_'//trim(str1(2:8))//'.sfl'))//" file cannot be opened")
    write(io_unit_nois) dtime*(nstep*1.0)/moo
    id = 0
    allocate(tvar(1:11,1:npoints), stat=ierr)
#ifdef PARALLEL
    end if
#endif
    do nb=1,nblocks
        top => mb_top(nb)

        xyz  => mb_xyz(nb)%fld
        sxyz => mb_sxyz(nb)%fld
        pv   => mb_pv(nb)%fld

        nregs = top%nregions
        do nr=1,nregs
            reg => top%bcs(nr)

            subtype_A  = reg%subtype_A
            if (subtype_A == 1) then
                s_st(:) = reg%s_st(:)
                s_ed(:) = reg%s_ed(:)
                s_nd    = reg%s_nd
                s_lr    = reg%s_lr

                st(:) = s_st(:)
                ed(:) = s_ed(:)

                m1 = m3x3(1,s_nd)
                m2 = m3x3(2,s_nd)
                m3 = m3x3(3,s_nd)

#ifdef PARALLEL
                packsize = product(ed(:)-st(:)+1)*11

                tag = nb*1000+nr

                pid = mb_top(nb)%pid
                if ((pid /= master .and. myid == pid) .or. myid == master) then
#endif
                    allocate(vars(st(1):ed(1),st(2):ed(2),st(3):ed(3),1:11), stat=ierr)
#ifdef PARALLEL
                end if

                if (myid == pid) then
#endif
                    do k=st(3),ed(3)
                    do j=st(2),ed(2)
                    do i=st(1),ed(1)
                        nx = sxyz(m1)%r3d(i,j,k)
                        ny = sxyz(m2)%r3d(i,j,k)
                        nz = sxyz(m3)%r3d(i,j,k)

                        vars(i,j,k,1)  = xyz(1)%r3d(i,j,k)
                        vars(i,j,k,2)  = xyz(2)%r3d(i,j,k)
                        vars(i,j,k,3)  = xyz(3)%r3d(i,j,k)
                        vars(i,j,k,4)  = nx
                        vars(i,j,k,5)  = ny
                        vars(i,j,k,6)  = nz
                        vars(i,j,k,7)  = pv(2)%r3d(i,j,k)*moo
                        vars(i,j,k,8)  = pv(3)%r3d(i,j,k)*moo
                        vars(i,j,k,9)  = pv(4)%r3d(i,j,k)*moo
                        vars(i,j,k,10) = pv(1)%r3d(i,j,k)
                        vars(i,j,k,11) = ( pv(5)%r3d(i,j,k) - poo )*moo*moo
                    end do
                    end do
                    end do
#ifdef PARALLEL
                end if

                if (pid /= master) then
                    if (myid == pid) then
                        call MPI_SEND(vars,packsize,kind_real_mpi, &
                                      master,tag,MPI_COMM_WORLD,ierr)
                    end if

                    if (myid == master) then
                        call MPI_RECV(vars,packsize,kind_real_mpi, &
                                      pid,tag,MPI_COMM_WORLD,status,ierr)
                    end if
                end if

                if (myid == master) then
#endif
                    do k=st(3),ed(3)
                    do j=st(2),ed(2)
                    do i=st(1),ed(1)
                        id=id+1
                    do m=1,11
                        tvar(m,id) = vars(i,j,k,m)
                    end do
                    end do
                    end do
                    end do
#ifdef PARALLEL
                end if

                if ((pid /= master .and. myid == pid) .or. myid == master) then
#endif
                    deallocate(vars, stat=ierr)
#ifdef PARALLEL
                end if
#endif
            end if
        end do
    end do
#ifdef PARALLEL
    if (myid == master) then
#endif
        write(io_unit_nois)((tvar(m,i),m=1,11),i=1,npoints)
        deallocate(tvar, stat=ierr)
        close(io_unit_nois,iostat=ierr)
        call error_check(ierr,trim(adjustl('noise/'//trim(casname)//'_'//trim(str1(2:8))//'.sfl'))//" file cannot be closed")
#ifdef PARALLEL
    end if
#endif
end subroutine to_noise_output

subroutine to_noise_output_sp()
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : m3x3
    use mod_constants, only : io_unit_nois
    use mod_datatypes, only : var_block_t,fld_array_t
    use mod_datatypes, only : top_block_t,bc_region_t
    use mod_variables, only : casname,vinf,nstep,reflgrd,dtime,poo,moo
    use mod_fieldvars, only : nblocks,mb_topsp,mb_xyzsp,mb_sxyzsp
    use mod_fieldvars, only : mb_pv,ntime,npoints
    use mod_parallels
    use mod_tecplotio
    implicit none
    integer(kind_int)          :: nb,nr,i,j,k,m,m1,m2,m3,ierr,id
    integer(kind_int)          :: s_st(3),s_ed(3),st(3),ed(3)
    integer(kind_int)          :: nregs,s_nd,s_lr,subtype_A
    real(kind_real)            :: nx,ny,nz
    type(top_block_t), pointer :: top
    type(bc_region_t), pointer :: reg
    type(fld_array_t), pointer :: xyz(:),sxyz(:)
    type(fld_array_t), pointer :: pv(:)
    real(kind_real),   pointer :: vars(:,:,:,:)
    real(kind_real),   pointer :: tvar(:,:)
    character(len=32 )         :: str1
#ifdef PARALLEL
    integer(kind_int) :: pid,packsize,tag
    integer(kind_int) :: status(MPI_STATUS_SIZE)
#endif

#ifdef PARALLEL
    if (myid == master) then
#endif
    write(str1, *) ntime + 10000000
    str1 = adjustl(str1)
    open (io_unit_nois,file='noise/'//trim(casname)//'_'//trim(str1(2:8))//'.sfl',form='unformatted',status='unknown',iostat=ierr)
    call error_check(ierr,trim(adjustl('noise/'//trim(casname)//'_'//trim(str1(2:8))//'.sfl'))//" file cannot be opened")
    write(io_unit_nois) dtime*(nstep*1.0)/moo
    id = 0
    allocate(tvar(1:11,1:npoints), stat=ierr)
#ifdef PARALLEL
    end if
#endif
    do nb=1,nblocks
        top => mb_topsp(nb)

        xyz  => mb_xyzsp(nb)%fld
        sxyz => mb_sxyzsp(nb)%fld
        pv   => mb_pv(nb)%fld

        nregs = top%nregions
        do nr=1,nregs
            reg => top%bcs(nr)

            subtype_A  = reg%subtype_A
            if (subtype_A == 1) then
                s_st(:) = reg%s_st(:)
                s_ed(:) = reg%s_ed(:)
                s_nd    = reg%s_nd
                s_lr    = reg%s_lr
                
                st(:) = s_st(:)
                ed(:) = s_ed(:)

                m1 = m3x3(1,s_nd)
                m2 = m3x3(2,s_nd)
                m3 = m3x3(3,s_nd)

#ifdef PARALLEL
                packsize = product(ed(:)-st(:)+1)*11

                tag = nb*1000+nr

                pid = mb_topsp(nb)%pid
                if ((pid /= master .and. myid == pid) .or. myid == master) then
#endif
                    allocate(vars(st(1):ed(1),st(2):ed(2),st(3):ed(3),1:11), stat=ierr)
#ifdef PARALLEL
                end if

                if (myid == pid) then
#endif
                    do k=st(3),ed(3)
                    do j=st(2),ed(2)
                    do i=st(1),ed(1)
                        nx = sxyz(m1)%r3d(i,j,k)
                        ny = sxyz(m2)%r3d(i,j,k)
                        nz = sxyz(m3)%r3d(i,j,k)

                        vars(i,j,k,1)  = xyz(1)%r3d(i,j,k)
                        vars(i,j,k,2)  = xyz(2)%r3d(i,j,k)
                        vars(i,j,k,3)  = xyz(3)%r3d(i,j,k)
                        vars(i,j,k,4)  = nx
                        vars(i,j,k,5)  = ny
                        vars(i,j,k,6)  = nz
                        vars(i,j,k,7)  = pv(2)%r3d(i,j,k)*moo
                        vars(i,j,k,8)  = pv(3)%r3d(i,j,k)*moo
                        vars(i,j,k,9)  = pv(4)%r3d(i,j,k)*moo
                        vars(i,j,k,10) = pv(1)%r3d(i,j,k)
                        vars(i,j,k,11) = ( pv(5)%r3d(i,j,k) - poo )*moo*moo
                    end do
                    end do
                    end do
#ifdef PARALLEL
                end if

                if (pid /= master) then
                    if (myid == pid) then
                        call MPI_SEND(vars,packsize,kind_real_mpi, &
                                      master,tag,MPI_COMM_WORLD,ierr)
                    end if

                    if (myid == master) then
                        call MPI_RECV(vars,packsize,kind_real_mpi, &
                                      pid,tag,MPI_COMM_WORLD,status,ierr)
                    end if
                end if

                if (myid == master) then
#endif
                    do k=st(3),ed(3)
                    do j=st(2),ed(2)
                    do i=st(1),ed(1)
                        id=id+1
                    do m=1,11
                        tvar(m,id) = vars(i,j,k,m)
                    end do
                    end do
                    end do
                    end do
#ifdef PARALLEL
                end if

                if ((pid /= master .and. myid == pid) .or. myid == master) then
#endif
                    deallocate(vars, stat=ierr)
#ifdef PARALLEL
                end if
#endif
            end if
        end do
    end do
#ifdef PARALLEL
    if (myid == master) then
#endif
        write(io_unit_nois)((tvar(m,i),m=1,11),i=1,npoints)
        deallocate(tvar, stat=ierr)
        close(io_unit_nois,iostat=ierr)
        call error_check(ierr,trim(adjustl('noise/'//trim(casname)//'_'//trim(str1(2:8))//'.sfl'))//" file cannot be closed")
#ifdef PARALLEL
    end if
#endif
end subroutine to_noise_output_sp

subroutine monitor_output()
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : io_unit_moni
    use mod_datatypes, only : var_block_t,fld_array_t
    use mod_datatypes, only : top_block_t
    use mod_variables, only : nstep,nmoni,nbmoni,nimoni,njmoni,nkmoni
    use mod_fieldvars, only : nblocks,mb_top,mb_xyz,mb_fmean
    use mod_fieldvars, only : mb_pv
    use mod_parallels
    use mod_tecplotio
    implicit none
    integer(kind_int)          :: nb,nr,i,j,k,m,n,ierr,id
    type(top_block_t), pointer :: top
    type(fld_array_t), pointer :: xyz(:)
    type(fld_array_t), pointer :: pv(:),fmean(:)
    real(kind_real),   pointer :: vars(:)
    real(kind_real),   pointer :: tvar(:,:)
    real(kind_real)            :: x,y,z,r,u,v,w,p,fr,fu,fv,fw,fp
#ifdef PARALLEL
    integer(kind_int) :: pid,packsize,tag
    integer(kind_int) :: status(MPI_STATUS_SIZE)
#endif

#ifdef PARALLEL
    if (myid == master) then
#endif
    call openfile(io_unit_moni,"monitor.dat","unknown","formatted","append")
    
    id = 0
    allocate(tvar(1:nmoni,1:8), stat=ierr)
#ifdef PARALLEL
    end if
#endif
    do m=1,nmoni
        nb  =  nbmoni(m)
        top => mb_top(nb)

        xyz  => mb_xyz(nb)%fld
        pv   => mb_pv(nb)%fld
        fmean=> mb_fmean(nb)%fld

#ifdef PARALLEL
        
        packsize = 8
        
        tag = nb*1000

        pid = mb_top(nb)%pid
        if ((pid /= master .and. myid == pid) .or. myid == master) then
#endif            
        allocate(vars(1:8), stat=ierr)
        
#ifdef PARALLEL
        end if

        if (myid == pid) then
#endif
            i        = nimoni(m)
            j        = njmoni(m)
            k        = nkmoni(m)
            vars(1)  = xyz(1)%r3d(i,j,k)
            vars(2)  = xyz(2)%r3d(i,j,k)
            vars(3)  = xyz(3)%r3d(i,j,k)
            vars(4)  = pv(1)%r3d(i,j,k)
            vars(5)  = pv(2)%r3d(i,j,k)
            vars(6)  = pv(3)%r3d(i,j,k)
            vars(7)  = pv(4)%r3d(i,j,k)
            vars(8)  = pv(5)%r3d(i,j,k)
#ifdef PARALLEL
        end if

        if (pid /= master) then
            if (myid == pid) then
                call MPI_SEND(vars,packsize,kind_real_mpi, &
                              master,tag,MPI_COMM_WORLD,ierr)
            end if

            if (myid == master) then
                call MPI_RECV(vars,packsize,kind_real_mpi, &
                              pid,tag,MPI_COMM_WORLD,status,ierr)
            end if
        end if

        if (myid == master) then
#endif
            id=id+1
            do n=1,8
                tvar(id,n) = vars(n)
            end do                
#ifdef PARALLEL
        end if

        if ((pid /= master .and. myid == pid) .or. myid == master) then
#endif
            deallocate(vars, stat=ierr)
#ifdef PARALLEL
        end if
#endif
    end do
#ifdef PARALLEL
    if (myid == master) then
#endif
        write(io_unit_moni,10) nstep
        write(io_unit_moni,20)((tvar(i,m),m=1,8),i=1,nmoni)
        deallocate(tvar, stat=ierr)
        close(io_unit_moni,iostat=ierr)
        call error_check(ierr,"monitor.dat file cannot be closed")
#ifdef PARALLEL
    end if
#endif


10 format(i7)
20 format(8(1x,e12.5))
    
end subroutine monitor_output

subroutine monitor_output_sp()
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : io_unit_moni
    use mod_datatypes, only : var_block_t,fld_array_t
    use mod_datatypes, only : top_block_t
    use mod_variables, only : nstep,nmoni,nbmoni,nimoni,njmoni,nkmoni
    use mod_fieldvars, only : nblocks,mb_topsp,mb_xyzsp,mb_fmean
    use mod_fieldvars, only : mb_pv
    use mod_parallels
    use mod_tecplotio
    implicit none
    integer(kind_int)          :: nb,nr,i,j,k,m,n,ierr,id
    type(top_block_t), pointer :: top
    type(fld_array_t), pointer :: xyz(:)
    type(fld_array_t), pointer :: pv(:),fmean(:)
    real(kind_real),   pointer :: vars(:)
    real(kind_real),   pointer :: tvar(:,:)
    real(kind_real)            :: x,y,z,r,u,v,w,p,fr,fu,fv,fw,fp
#ifdef PARALLEL
    integer(kind_int) :: pid,packsize,tag
    integer(kind_int) :: status(MPI_STATUS_SIZE)
#endif

#ifdef PARALLEL
    if (myid == master) then
#endif
    call openfile(io_unit_moni,"monitor.dat","unknown","formatted","append")
    
    id = 0
    allocate(tvar(1:nmoni,1:8), stat=ierr)
#ifdef PARALLEL
    end if
#endif
    do m=1,nmoni
        nb  =  nbmoni(m)
        top => mb_topsp(nb)

        xyz  => mb_xyzsp(nb)%fld
        pv   => mb_pv(nb)%fld
        fmean=> mb_fmean(nb)%fld

#ifdef PARALLEL
        
        packsize = 8
        
        tag = nb*1000

        pid = mb_topsp(nb)%pid
        if ((pid /= master .and. myid == pid) .or. myid == master) then
#endif            
        allocate(vars(1:8), stat=ierr)
        
#ifdef PARALLEL
        end if

        if (myid == pid) then
#endif
            i        = nimoni(m)
            j        = njmoni(m)
            k        = nkmoni(m)
            vars(1)  = xyz(1)%r3d(i,j,k)
            vars(2)  = xyz(2)%r3d(i,j,k)
            vars(3)  = xyz(3)%r3d(i,j,k)
            vars(4)  = pv(1)%r3d(i,j,k)
            vars(5)  = pv(2)%r3d(i,j,k)
            vars(6)  = pv(3)%r3d(i,j,k)
            vars(7)  = pv(4)%r3d(i,j,k)
            vars(8)  = pv(5)%r3d(i,j,k)
#ifdef PARALLEL
        end if

        if (pid /= master) then
            if (myid == pid) then
                call MPI_SEND(vars,packsize,kind_real_mpi, &
                              master,tag,MPI_COMM_WORLD,ierr)
            end if

            if (myid == master) then
                call MPI_RECV(vars,packsize,kind_real_mpi, &
                              pid,tag,MPI_COMM_WORLD,status,ierr)
            end if
        end if

        if (myid == master) then
#endif
            id=id+1
            do n=1,8
                tvar(id,n) = vars(n)
            end do                
#ifdef PARALLEL
        end if

        if ((pid /= master .and. myid == pid) .or. myid == master) then
#endif
            deallocate(vars, stat=ierr)
#ifdef PARALLEL
        end if
#endif
    end do
#ifdef PARALLEL
    if (myid == master) then
#endif
        write(io_unit_moni,10) nstep
        write(io_unit_moni,20)((tvar(i,m),m=1,8),i=1,nmoni)
        deallocate(tvar, stat=ierr)
        close(io_unit_moni,iostat=ierr)
        call error_check(ierr,"monitor.dat file cannot be closed")
#ifdef PARALLEL
    end if
#endif


10 format(i7)
20 format(8(1x,e12.5))
    
end subroutine monitor_output_sp

!!subroutine monitor_output()
!!    use mod_kndconsts, only : kind_int,kind_real
!!    use mod_variables, only : nmoni
!!    use mod_variables, only : nstep,ns_mean,nprms
!!    implicit none
!!    integer(kind_int)          :: m
!!
!!    if(nprms>0 .and. nstep>ns_mean) then
!!        do m=1,nmoni
!!            call output_monitor_unsteady(m)
!!        end do
!!    end if
!!    
!!end subroutine monitor_output


