
subroutine preset
    use mod_variables, only : ncutpol,ncelcet,nvis,nsponge
    use mod_fieldvars, only : ninters,ninterf
    implicit none
    external :: input_par,input_top,input_top_cc
    external :: reset_file_status

    call run_seq_and_master(input_par)
    call msg_seq_and_master("Input the parameters successfully")
    call broadcast_par
    call msg_on_master("Broadcast the parameters successfully")

    call reset_par
    call reset_turconsts

    call init_inflow
    
    if(ncelcet /= 1) then
        call run_seq_and_master(input_top)
        call msg_seq_and_master("Input the topology successfully")
        call broadcast_top
        call msg_on_master("Broadcast the topology successfully")  
        
        call analyze_top
        call set_bc_sequence
        call set_bc_inters
        call set_block_coms
        call msg_seq_and_master("Analyze the topology successfully")    
        
        call set_flag_surface
        call msg_seq_and_master("Set the block surface flags successfully")
        
        call input_grd
        call msg_seq_and_master("Input the grid successfully")        
    else  
        call run_seq_and_master(input_top_cc)
        call msg_seq_and_master("Input the topology successfully")
        call broadcast_top_cc
        call msg_on_master("Broadcast the dual-mesh topology successfully") 
        
        call analyze_top_cc
        call set_bc_sequence_cc
        call set_bc_inters_cc
        call set_block_coms_cc
        call msg_seq_and_master("Analyze the dual-mesh topology successfully")    
        
        call set_flag_surface_cc
        call msg_seq_and_master("Set the dual-mesh block surface flags successfully")
        
        call input_grd_cc
        call msg_seq_and_master("Input the grid successfully")
    end if

    if(ninters > ninterf) then
        call set_bc_patched
        call analyze_top_patched
        call msg_seq_and_master("Analyze the patched topology successfully")
        call set_bc_inters_patched
        call set_patched_offsets
        call msg_seq_and_master("Initialize the patched boundaries successfully")
    end if
    
    if(ncutpol == 5) then
        call grid_derivative_exp
    else
        if(ncelcet == 1) then
            call grid_derivative_cc
            call reset_grid_derivative_cc
            call msg_seq_and_master("Compute the dual-mesh grid derivatives successfully")
        else
            call grid_derivative
            call reset_grid_derivative
            call msg_seq_and_master("Compute the grid derivatives successfully")
        end if        
    end if
    
    if(ncelcet /= 1) then
        call search_singulars
        call msg_seq_and_master("Search the singulars successfully")
    
        call allocate_fld_variables

        call set_wall_dist

        call set_sponge_dist

        call init_fld_variables

        call run_seq_and_master(reset_file_status)
    else
        call search_singulars
        call msg_seq_and_master("Search the singulars successfully") 
        
        call allocate_fld_variables_sp
        
        !todo call set_wall_dist

        !todo call set_sponge_dist
        if(nvis > 1 .or. nsponge > 0) then
            call msg_seq_and_master("Turbulence models or sponge boundaries are not ready for dual-mesh grid")
            call env_finalize
        end if
        
        call init_fld_variables_sp
        
        call run_seq_and_master(reset_file_status)
    end if


!  CHEyg added
! #define MemOPT1	! 优化sxyz访问,在makefile中定义

#ifdef  MemOPT1	
!   将本进程各块的 mb_sxyz 缓存在 sxyz1 中   
    !write(*,*) "in preset before calling save_mod_opt_vars"
    call save_mod_opt_vars() 
    !write(*,*) "in preset after calling save_mod_opt_vars"
#endif	
		

end subroutine preset

subroutine reset_par
    use mod_constants, only : zero,one
    use mod_constants, only : nscheme_policys
    use mod_constants, only : nscheme_con_d1der,nscheme_con_d1int
    use mod_constants, only : nscheme_vis_d2der,nscheme_vis_d2int
    use mod_constants, only : nscheme_vis_d3der,nscheme_vis_d3int
    use mod_constants, only : nscheme_grd_d2der,nscheme_grd_d2int
    use mod_constants, only : nscheme_grd_d3der,nscheme_grd_d3int
    use mod_constants, only : nlhs_prsgs_ust_mat,nlhs_prsgs_std_mat
    use mod_constants, only : nvis_tur_sa,nvis_tur_sst,nvis_tur_hst
    use mod_constants, only : ndim_2d,nsw_dir_close,nsw_dir_open
    use mod_variables, only : nvis,ndim,nsw_kdir,fsw_kdir
    use mod_variables, only : nscheme,nd1der_con,nd1int_con
    use mod_variables, only : nd2der_vis,nd2int_vis,nd3der_vis,nd3int_vis
    use mod_variables, only : nd2der_grd,nd2int_grd,nd3der_grd,nd3int_grd
    use mod_variables, only : ntursch,nsw_dst,nsw_bld,nunst
    use mod_variables, only : nd1der_tur_con,nd1int_tur_con
    use mod_variables, only : nd2der_tur_vis,nd2int_tur_vis
    use mod_variables, only : nd3der_tur_vis,nd3int_tur_vis
    use mod_fieldvars, only : nblocks,nprocs,neqn,npvs
    use mod_fieldvars, only : neqt,npvt,nqvst
    implicit none

    select case(ndim)
    case(ndim_2d)
        nsw_kdir = nsw_dir_close
    case default
        nsw_kdir = nsw_dir_open
    end select

    if (nsw_kdir == nsw_dir_close) then
        fsw_kdir = zero
    else
        fsw_kdir = one
    end if

    nd1der_con = nscheme_policys(nscheme_con_d1der, nscheme)
    nd1int_con = nscheme_policys(nscheme_con_d1int, nscheme)
    nd2der_vis = nscheme_policys(nscheme_vis_d2der, nscheme)
    nd2int_vis = nscheme_policys(nscheme_vis_d2int, nscheme)
    nd3der_vis = nscheme_policys(nscheme_vis_d3der, nscheme)
    nd3int_vis = nscheme_policys(nscheme_vis_d3int, nscheme)
    nd2der_grd = nscheme_policys(nscheme_grd_d2der, nscheme)
    nd2int_grd = nscheme_policys(nscheme_grd_d2int, nscheme)
    nd3der_grd = nscheme_policys(nscheme_grd_d3der, nscheme)
    nd3int_grd = nscheme_policys(nscheme_grd_d3int, nscheme)

    nblocks = 0
    nprocs  = 0

!#define PARAM_OPT     !CHEyg  参数化优化  在makefile中定义

#ifndef PARAM_OPT
     neqn = 5
     npvs = neqn  ! 要求（npvs=neqn）。
#endif	 

    select case(nvis)
    case(nvis_tur_sa)
        neqt    = 1
        nsw_dst = 1
        nsw_bld = 0
    case(nvis_tur_sst,nvis_tur_hst)
        neqt    = 2
        nsw_dst = 1
        nsw_bld = 1
    case default
        neqt    = 0
        nsw_dst = 0
        nsw_bld = 0
    end select
    
    if(nunst==4) nsw_dst=1

    npvt  = 4 + neqt
    nqvst = 1 + neqt

    nd1der_tur_con = nscheme_policys(nscheme_con_d1der, ntursch)
    nd1int_tur_con = nscheme_policys(nscheme_con_d1int, ntursch)
    nd2der_tur_vis = nscheme_policys(nscheme_vis_d2der, ntursch)
    nd2int_tur_vis = nscheme_policys(nscheme_vis_d2int, ntursch)
    nd3der_tur_vis = nscheme_policys(nscheme_vis_d3der, ntursch)
    nd3int_tur_vis = nscheme_policys(nscheme_vis_d3int, ntursch)

end subroutine reset_par

subroutine allocate_fld_variables
    use mod_kndconsts, only : kind_int
    use mod_constants, only : nvis_euler,nprec_non
    use mod_constants, only : nlhs_rkutta_3step,nlhs_prsgs_std_mat
    use mod_constants, only : nlhs_prsgs_ust_mat,nlhs_prsgs_ust_sca,nlhs_prsgs_ust_sca_scmp
    use mod_variables, only : nvis,nlhs,nghnode,nprms,pltsty,nprec
    use mod_variables, only : nsw_dst,nsw_bld,nsponge
    use mod_fieldvars, only : neqn,npvs,neqt,npvt,nqvst
    use mod_fieldvars, only : mb_qc,mb_pv,mb_ws,mb_t,mb_c
    use mod_fieldvars, only : mb_rhs,mb_dpv,mb_dws,mb_vsl,mb_vst
    use mod_fieldvars, only : mb_src,mb_src_prec,mb_srv,mb_rdt,mb_qm
    !use mod_fieldvars, only : mb_nonprec,mb_prec_diag_inv
    use mod_fieldvars, only : mb_dt,mb_dq,mb_q0,mb_aiv
    use mod_fieldvars, only : mb_qt,mb_pvt,mb_dqt,mb_rhst
    use mod_fieldvars, only : mb_srt,mb_rtur,mb_dtur
    use mod_fieldvars, only : mb_frms,mb_fmean
    use mod_fieldvars, only : mb_dst,mb_bld,mb_qvst,mb_dsp
    use mod_interface, only : mb_var_create
    use mod_interface, only : mb_var_pointer_create
    use mod_interface, only : mb_var_pointer_assign
    implicit none

    call mb_var_create(mb_qc,1,neqn,nghnode)
    call mb_var_create(mb_pv,1,npvs,nghnode)
    
    if (pltsty == 5) then
        call mb_var_create(mb_ws,1,10,nghnode)
    end if

    call mb_var_create(mb_t  ,1,1,nghnode)
    call mb_var_create(mb_c  ,1,1,nghnode)

    call mb_var_create(mb_rhs,1,neqn,nghnode)

    call mb_var_create(mb_src,1,3,nghnode)
    
    !todo 预处理相关
    if (nprec > nprec_non) then
        !call mb_var_create(mb_prec,1,25,nghnode)!创建预处理矩阵
        !call mb_var_create(mb_prec_inv,1,25,nghnode)!创建预处理矩阵的逆矩阵
        !call mb_var_create(mb_nonprec,1,4,nghnode)!创建守恒原始变量转换矩阵
        call mb_var_create(mb_src_prec,1,3,nghnode)!创建预处理后的对流项谱半径
        !call mb_var_create(mb_prec_diag_inv,1,25,nghnode)!创建预处理后的D矩阵的逆矩阵
    end if

    if (nvis > nvis_euler) then
        call mb_var_create(mb_dpv,1,12,nghnode)
        
        if (pltsty == 5) then
            call mb_var_create(mb_dws,1,30,nghnode)
        end if
        
        call mb_var_create(mb_vsl,1,1,nghnode)
        call mb_var_create(mb_vst,1,1,nghnode)
        call mb_var_create(mb_srv ,1,3,nghnode)
    end if

    call mb_var_create(mb_rdt,1,1,nghnode)

    call mb_var_create(mb_dt ,1,1,nghnode)
    call mb_var_create(mb_dq ,1,neqn,nghnode)

    select case(nlhs)
    case(nlhs_rkutta_3step)
        call mb_var_create(mb_q0,1,neqn,nghnode)
    case(nlhs_prsgs_ust_sca,nlhs_prsgs_ust_sca_scmp)
        call mb_var_create(mb_q0,1,neqn,nghnode)
        call mb_var_create(mb_qm,1,neqn,nghnode)        
    case(nlhs_prsgs_ust_mat)
        call mb_var_create(mb_q0,1,neqn,nghnode)
        call mb_var_create(mb_qm,1,neqn,nghnode) 
        call mb_var_create(mb_aiv,1,neqn*neqn,nghnode) 
    case(nlhs_prsgs_std_mat)
        call mb_var_create(mb_aiv,1,neqn*neqn,nghnode) 
    end select

    if (neqt > 0) then
        call mb_var_create(mb_qt ,1,neqt,nghnode)
        call mb_var_create(mb_dqt,1,neqt,nghnode)

        call mb_var_create(mb_rhst,1,neqt,nghnode)

        call mb_var_create(mb_srt,1,3,nghnode)

        call mb_var_create(mb_rtur,1,neqt,nghnode)

        call mb_var_create(mb_dtur,1,neqt*3,nghnode)

        call mb_var_pointer_create(mb_pvt,1,npvt)
        call mb_var_pointer_assign(mb_pvt,1,4,mb_pv,1,4)
        call mb_var_pointer_assign(mb_pvt,5,npvt,mb_qt,1,neqt)

        call mb_var_pointer_create(mb_qvst,1,nqvst)
        call mb_var_pointer_assign(mb_qvst,1,neqt,mb_qt,1,neqt)
        call mb_var_pointer_assign(mb_qvst,neqt+1,neqt+1,mb_vst,1,1)
    end if

    if (nsw_dst > 0) then
        call mb_var_create(mb_dst,1,2,nghnode)
    end if

    if (nsw_bld > 0) then
        call mb_var_create(mb_bld,1,1,nghnode)
    end if

    if (nsponge > 0) then
        call mb_var_create(mb_dsp,1,2,nghnode)
    end if
    
    if (nprms > 0) then
        call mb_var_create(mb_fmean,1,npvs,nghnode)
        call mb_var_create(mb_frms ,1,npvs,nghnode)
    end if

end subroutine allocate_fld_variables

subroutine allocate_fld_variables_sp
    use mod_kndconsts, only : kind_int
    use mod_constants, only : nvis_euler,nprec_non
    use mod_constants, only : nlhs_rkutta_3step,nlhs_prsgs_std_mat
    use mod_constants, only : nlhs_prsgs_ust_mat,nlhs_prsgs_ust_sca,nlhs_prsgs_ust_sca_scmp
    use mod_variables, only : nvis,nlhs,nghnode,nprms,pltsty,nprec
    use mod_variables, only : nsw_dst,nsw_bld,nsponge
    use mod_fieldvars, only : neqn,npvs,neqt,npvt,nqvst
    use mod_fieldvars, only : mb_qc,mb_pv,mb_pvfp,mb_t,mb_c !,mb_ws
    use mod_fieldvars, only : mb_rhs,mb_dpv,mb_vsl,mb_vst !,mb_dws
    use mod_fieldvars, only : mb_src,mb_src_prec,mb_srv,mb_rdt,mb_qm
    !use mod_fieldvars, only : mb_nonprec,mb_prec_diag_inv
    use mod_fieldvars, only : mb_dt,mb_dq,mb_q0,mb_aiv
    use mod_fieldvars, only : mb_qt,mb_pvt,mb_dqt,mb_rhst
    use mod_fieldvars, only : mb_srt,mb_rtur,mb_dtur
    use mod_fieldvars, only : mb_frms,mb_fmean
    use mod_fieldvars, only : mb_dst,mb_bld,mb_qvst,mb_dsp
    use mod_interface, only : mb_var_create_sp,mb_var_create
    use mod_interface, only : mb_var_pointer_create
    use mod_interface, only : mb_var_pointer_assign
    implicit none

    call mb_var_create_sp(mb_qc,1,neqn,nghnode)
    call mb_var_create_sp(mb_pv,1,npvs,nghnode)
    call mb_var_create(mb_pvfp,1,npvs,nghnode)
    
    if (pltsty == 5) then
        !call mb_var_create_sp(mb_ws,1,10,nghnode)
        call msg_seq_and_master("The Tecplot output style is not ready for dual-mesh grid")
    end if

    call mb_var_create_sp(mb_t  ,1,1,nghnode)
    call mb_var_create_sp(mb_c  ,1,1,nghnode)

    call mb_var_create_sp(mb_rhs,1,neqn,nghnode)

    call mb_var_create_sp(mb_src,1,3,nghnode)
    
    !todo 预处理相关
    if (nprec > nprec_non) then
        !call mb_var_create(mb_prec,1,25,nghnode)!创建预处理矩阵
        !call mb_var_create(mb_prec_inv,1,25,nghnode)!创建预处理矩阵的逆矩阵
        !call mb_var_create(mb_nonprec,1,4,nghnode)!创建守恒原始变量转换矩阵
        call mb_var_create_sp(mb_src_prec,1,3,nghnode)!创建预处理后的对流项谱半径
        !call mb_var_create(mb_prec_diag_inv,1,25,nghnode)!创建预处理后的D矩阵的逆矩阵
    end if

    if (nvis > nvis_euler) then
        call mb_var_create_sp(mb_dpv,1,12,nghnode)
        
        if (pltsty == 5) then
            ! call mb_var_create_sp(mb_dws,1,30,nghnode)
            call msg_seq_and_master("The Tecplot output style is not ready for dual-mesh grid")
        end if
        
        call mb_var_create_sp(mb_vsl,1,1,nghnode)
        call mb_var_create_sp(mb_vst,1,1,nghnode)
        call mb_var_create_sp(mb_srv,1,3,nghnode)
    end if

    call mb_var_create_sp(mb_rdt,1,1,nghnode)

    call mb_var_create_sp(mb_dt ,1,1,nghnode)
    call mb_var_create_sp(mb_dq ,1,neqn,nghnode)

    select case(nlhs)
    case(nlhs_rkutta_3step)
        call mb_var_create_sp(mb_q0,1,neqn,nghnode)
    case(nlhs_prsgs_ust_sca,nlhs_prsgs_ust_sca_scmp)
        call mb_var_create_sp(mb_q0,1,neqn,nghnode)
        call mb_var_create_sp(mb_qm,1,neqn,nghnode)        
    case(nlhs_prsgs_ust_mat)
        call mb_var_create_sp(mb_q0,1,neqn,nghnode)
        call mb_var_create_sp(mb_qm,1,neqn,nghnode) 
        call mb_var_create_sp(mb_aiv,1,neqn*neqn,nghnode) 
    case(nlhs_prsgs_std_mat)
        call mb_var_create_sp(mb_aiv,1,neqn*neqn,nghnode) 
    end select

    if (neqt > 0) then
        call mb_var_create_sp(mb_qt ,1,neqt,nghnode)
        call mb_var_create_sp(mb_dqt,1,neqt,nghnode)

        call mb_var_create_sp(mb_rhst,1,neqt,nghnode)

        call mb_var_create_sp(mb_srt,1,3,nghnode)

        call mb_var_create_sp(mb_rtur,1,neqt,nghnode)

        call mb_var_create_sp(mb_dtur,1,neqt*3,nghnode)

        call mb_var_pointer_create(mb_pvt,1,npvt)
        call mb_var_pointer_assign(mb_pvt,1,4,mb_pv,1,4)
        call mb_var_pointer_assign(mb_pvt,5,npvt,mb_qt,1,neqt)

        call mb_var_pointer_create(mb_qvst,1,nqvst)
        call mb_var_pointer_assign(mb_qvst,1,neqt,mb_qt,1,neqt)
        call mb_var_pointer_assign(mb_qvst,neqt+1,neqt+1,mb_vst,1,1)
    end if

    if (nsw_dst > 0) then
        call mb_var_create_sp(mb_dst,1,2,nghnode)
    end if

    if (nsw_bld > 0) then
        call mb_var_create_sp(mb_bld,1,1,nghnode)
    end if

    if (nsponge > 0) then
        call mb_var_create_sp(mb_dsp,1,2,nghnode)
    end if
    
    if (nprms > 0) then
        call mb_var_create_sp(mb_fmean,1,npvs,nghnode)
        call mb_var_create_sp(mb_frms ,1,npvs,nghnode)
    end if

end subroutine allocate_fld_variables_sp

