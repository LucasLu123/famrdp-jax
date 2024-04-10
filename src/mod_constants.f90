
module mod_constants
    use mod_kndconsts,only : kind_double,kind_int
    implicit none

    real(kind_double), parameter :: small = 1.0e-30
    real(kind_double), parameter :: large = 1.0e+30
    real(kind_double), parameter :: epsil = 1.0e-6
    real(kind_double), parameter :: SCMP_sigma  = 1.56

    real(kind_double), parameter :: zero   = 0.0
    real(kind_double), parameter :: tenth  = 0.1
    real(kind_double), parameter :: eighth = 0.125
    real(kind_double), parameter :: fourth = 0.25
    real(kind_double), parameter :: half   = 0.5
    real(kind_double), parameter :: one    = 1.0
    real(kind_double), parameter :: thr2nd = 1.5
    real(kind_double), parameter :: two    = 2.0
    real(kind_double), parameter :: three  = 3.0
    real(kind_double), parameter :: four   = 4.0
    real(kind_double), parameter :: five   = 5.0
    real(kind_double), parameter :: six    = 6.0
    real(kind_double), parameter :: seven  = 7.0
    real(kind_double), parameter :: eight  = 8.0
    real(kind_double), parameter :: nine   = 9.0
    real(kind_double), parameter :: ten    = 10.0
    real(kind_double), parameter :: twenty = 20.0
    real(kind_double), parameter :: sixth  = one/six
    real(kind_double), parameter :: third  = one/three
    real(kind_double), parameter :: two3rd = two/three
    real(kind_double), parameter :: four3rd= four/three
    real(kind_double), parameter :: two9th = two/nine

    real(kind_double), parameter :: pai     = 3.1415926535897932384626433832795_kind_double
    real(kind_double), parameter :: deg2rad = pai/180.0
    real(kind_double), parameter :: rad2deg = 180.0/pai

    real(kind_double), parameter :: runiv = 8.31434

    real(kind_double), parameter :: sml_ssf  = 1.0e-30
    real(kind_double), parameter :: sml_vol  = 1.0e-30
    real(kind_double), parameter :: pole_ssf = 1.0e-16
    real(kind_double), parameter :: pole_vol = 1.0e-16
    real(kind_double), parameter :: zero_vol = 2.0*pole_vol

    integer(kind_int), parameter :: mide(3,3) = reshape((/1,0,0,0,1,0,0,0,1/),(/3,3/))
    integer(kind_int), parameter :: mcyc(3,3) = reshape((/1,2,3,2,3,1,3,1,2/),(/3,3/))
    integer(kind_int), parameter :: m3x3(3,3) = reshape((/1,2,3,4,5,6,7,8,9/),(/3,3/))
    !integer(kind_int), parameter :: bcsp(3,3) = reshape((/0,1,1,1,0,1,1,1,0/),(/3,3/))

    integer(kind_int), parameter :: io_unit_par = 11
    integer(kind_int), parameter :: io_unit_top = 12
    integer(kind_int), parameter :: io_unit_grd = 13
    integer(kind_int), parameter :: io_unit_fld = 14
    integer(kind_int), parameter :: io_unit_res = 15
    integer(kind_int), parameter :: io_unit_fce = 16
    integer(kind_int), parameter :: io_unit_tec = 17
    integer(kind_int), parameter :: io_unit_wfl = 18
    integer(kind_int), parameter :: io_unit_aux = 19
    integer(kind_int), parameter :: io_unit_tur = 20
    integer(kind_int), parameter :: io_unit_dst = 21
    integer(kind_int), parameter :: io_unit_dsp = 22
    integer(kind_int), parameter :: io_unit_men = 23
    integer(kind_int), parameter :: io_unit_rms = 24
    integer(kind_int), parameter :: io_unit_fmen = 25
    integer(kind_int), parameter :: io_unit_frms = 26
    integer(kind_int), parameter :: io_unit_sin = 27
    integer(kind_int), parameter :: io_unit_acou = 28
    integer(kind_int), parameter :: io_unit_nois = 29
    integer(kind_int), parameter :: io_unit_moni = 30
    integer(kind_int), parameter :: io_unit_tfld = 31

    integer(kind_int), parameter :: fldtype_i3d  = 3
    integer(kind_int), parameter :: fldtype_none = 0
    integer(kind_int), parameter :: fldtype_r3d  = -3

    integer(kind_int), parameter :: nplot3d_std = 0
    integer(kind_int), parameter :: nplot3d_fun = 1

    integer(kind_int), parameter :: ndim_axi = 1
    integer(kind_int), parameter :: ndim_2d  = 2
    integer(kind_int), parameter :: ndim_3d  = 3

    integer(kind_int), parameter :: nsw_dir_close = 0
    integer(kind_int), parameter :: nsw_dir_open  = 1

    integer(kind_int), parameter :: id_t = 0    ! direction time
    integer(kind_int), parameter :: id_x = 1    ! direction x
    integer(kind_int), parameter :: id_y = 2    ! direction y
    integer(kind_int), parameter :: id_z = 3    ! direction z

    integer(kind_int), parameter :: id_i = 1    ! direction i
    integer(kind_int), parameter :: id_j = 2    ! direction j
    integer(kind_int), parameter :: id_k = 3    ! direction k

    integer(kind_int), parameter :: id_ms = 1    ! density (ro)
    integer(kind_int), parameter :: id_mx = 2    ! ro * u
    integer(kind_int), parameter :: id_my = 3    ! ro * v
    integer(kind_int), parameter :: id_mz = 4    ! ro * w
    integer(kind_int), parameter :: id_en = 5    ! ro * e
    integer(kind_int), parameter :: id_rs = 6    ! ro * fs

    integer(kind_int), parameter :: id_ro = 1    ! density (ro)
    integer(kind_int), parameter :: id_u  = 2    ! u : velocity components
    integer(kind_int), parameter :: id_v  = 3    ! v
    integer(kind_int), parameter :: id_w  = 4    ! w
    integer(kind_int), parameter :: id_ps = 5    ! p : static pressure
    integer(kind_int), parameter :: id_fs = 6    ! fs

    integer(kind_int), parameter :: id_ux = 1    ! gradients: du/dx
    integer(kind_int), parameter :: id_vx = 2    ! dv/dx
    integer(kind_int), parameter :: id_wx = 3    ! dw/dx
    integer(kind_int), parameter :: id_tx = 4    ! dt/dx
    integer(kind_int), parameter :: id_uy = 5    ! du/dy
    integer(kind_int), parameter :: id_vy = 6    ! dv/dy
    integer(kind_int), parameter :: id_wy = 7    ! dw/dy
    integer(kind_int), parameter :: id_ty = 8    ! dt/dy
    integer(kind_int), parameter :: id_uz = 9    ! du/dz
    integer(kind_int), parameter :: id_vz = 10   ! dv/dz
    integer(kind_int), parameter :: id_wz = 11   ! dw/dz
    integer(kind_int), parameter :: id_tz = 12   ! dt/dz

    integer(kind_int), parameter :: nvis_euler   = 0
    integer(kind_int), parameter :: nvis_ns_lam  = 1
    integer(kind_int), parameter :: nvis_tur_sa  = 2
    integer(kind_int), parameter :: nvis_tur_sst = 3
    integer(kind_int), parameter :: nvis_tur_hst = 4

    integer(kind_int), parameter :: nrestrt_restart  = 0
    integer(kind_int), parameter :: nrestrt_cont_lam = 1
    integer(kind_int), parameter :: nrestrt_cont_tur = 2
    
    !todo 预处理相关
    integer(kind_int), parameter :: nprec_non = 0!预处理标识
    
    !todo SCM-P相关
    integer(kind_int), parameter :: nscmp_non = 0!SCM-P标识

    integer(kind_int), parameter :: nincst_close = 0
    integer(kind_int), parameter :: nincst_inter = 1
    integer(kind_int), parameter :: nincst_intbc = 2

    integer(kind_int), parameter :: nfsf_id_max = 11
    integer(kind_int), parameter :: nfsf_con_d1der  = 1   ! 无粘项通量导数
    integer(kind_int), parameter :: nfsf_con_d1int  = 2   ! 无粘项线性插值(网格导数)
    integer(kind_int), parameter :: nfsf_con_d1nint = 3   ! 无粘项非线性插值
    integer(kind_int), parameter :: nfsf_vis_d2der  = 4   ! 粘性项通量导数
    integer(kind_int), parameter :: nfsf_vis_d2int  = 5   ! 粘性项线性插值(变量、速度/温度导数)
    integer(kind_int), parameter :: nfsf_vis_d3der  = 6   ! 速度/温度导数之求导
    integer(kind_int), parameter :: nfsf_vis_d3int  = 7   ! 速度/温度导数之插值
    integer(kind_int), parameter :: nfsf_grd_d2der  = 8   ! 网格导数d2求导
    integer(kind_int), parameter :: nfsf_grd_d2int  = 9   ! 网格导数d2插值
    integer(kind_int), parameter :: nfsf_grd_d3der  = 10  ! 网格导数d3求导
    integer(kind_int), parameter :: nfsf_grd_d3int  = 11  ! 网格导数d3插值


    integer(kind_int), parameter :: nbc_intbc_scheme = -1   ! 内边界格式标记（新增，JUN 10,2011）
    integer(kind_int), parameter :: nbc_inter_scheme = 0    ! 内点格式标记
    integer(kind_int), parameter :: nbc_bound_scheme = 1    ! 外边界格式标记

    integer(kind_int), parameter :: nbc_policy_max = 5
    integer(kind_int), parameter :: nbc_policy_fullext = 1    ! 全边界格式
    integer(kind_int), parameter :: nbc_policy_partext = 2    ! 部分边界格式
    integer(kind_int), parameter :: nbc_policy_partint = 3    ! 部分内点格式
    integer(kind_int), parameter :: nbc_policy_fullint = 4    ! 全内点格式
    integer(kind_int), parameter :: nbc_policy_fullexp = 5    ! 网格外推全内点格式
    integer(kind_int), parameter :: nfsf_bc_policys(nfsf_id_max,nbc_policy_max) = &
        reshape((/ nbc_bound_scheme, nbc_bound_scheme, nbc_bound_scheme, &   ! nbc_policy = nbc_policy_fullext
                   nbc_bound_scheme, nbc_bound_scheme, nbc_bound_scheme, nbc_bound_scheme, &
                   nbc_bound_scheme, nbc_bound_scheme, nbc_bound_scheme, nbc_bound_scheme, &
                   nbc_bound_scheme, nbc_bound_scheme, nbc_inter_scheme, &   ! nbc_policy = nbc_policy_partext
                   nbc_bound_scheme, nbc_bound_scheme, nbc_bound_scheme, nbc_inter_scheme, &
                   nbc_bound_scheme, nbc_bound_scheme, nbc_bound_scheme, nbc_bound_scheme, &
                   nbc_intbc_scheme, nbc_intbc_scheme, nbc_inter_scheme, &   ! nbc_policy = nbc_policy_partint
                   nbc_intbc_scheme, nbc_inter_scheme, nbc_intbc_scheme, nbc_inter_scheme, &
                   nbc_intbc_scheme, nbc_intbc_scheme, nbc_intbc_scheme, nbc_intbc_scheme, &
                   nbc_inter_scheme, nbc_inter_scheme, nbc_inter_scheme, &   ! nbc_policy = nbc_policy_fullint
                   nbc_inter_scheme, nbc_inter_scheme, nbc_inter_scheme, nbc_inter_scheme, &
                   nbc_inter_scheme, nbc_inter_scheme, nbc_inter_scheme, nbc_inter_scheme, &
                   nbc_inter_scheme, nbc_inter_scheme, nbc_inter_scheme, &   ! nbc_policy = nbc_policy_fullexp
                   nbc_inter_scheme, nbc_inter_scheme, nbc_inter_scheme, nbc_inter_scheme, &
                   nbc_inter_scheme, nbc_inter_scheme, nbc_inter_scheme, nbc_inter_scheme /),(/nfsf_id_max,nbc_policy_max/))


    integer(kind_int), parameter :: nderive_edge2  = 1
    integer(kind_int), parameter :: nderive_edge4  = 2
    integer(kind_int), parameter :: nderive_edge6  = 3
    integer(kind_int), parameter :: nderive_ehen4  = 4
    integer(kind_int), parameter :: nderive_ehen6  = 5
    integer(kind_int), parameter :: nderive_ehen6e = 6
    integer(kind_int), parameter :: nderive_ehcs6  = 7
    integer(kind_int), parameter :: nderive_ehcs6e = 8
    integer(kind_int), parameter :: nderive_ehen8  = 9
    integer(kind_int), parameter :: nderive_ehen8e = 10
    integer(kind_int), parameter :: nderive_scsl4  = 11
    integer(kind_int), parameter :: nderive_scsl4e = 12
    integer(kind_int), parameter :: nderive_scsl6  = 13
    integer(kind_int), parameter :: nderive_scsl6e = 14

    integer(kind_int), parameter :: nintplt_node2   = 1
    integer(kind_int), parameter :: nintplt_node4   = 2
    integer(kind_int), parameter :: nintplt_node6   = 3
    integer(kind_int), parameter :: nintplt_node6w  = 4
    integer(kind_int), parameter :: nintplt_node6e  = 5
    integer(kind_int), parameter :: nintplt_node8   = 6
    integer(kind_int), parameter :: nintplt_node8e  = 7
    integer(kind_int), parameter :: nintplt_scsl4   = 8
    integer(kind_int), parameter :: nintplt_scsl4e  = 9
    integer(kind_int), parameter :: nintplt_scsl6   = 10
    integer(kind_int), parameter :: nintplt_scsl6e  = 11

    integer(kind_int), parameter :: nintnon_muscl2pv = 1
    integer(kind_int), parameter :: nintnon_wcns5pv  = 2
    integer(kind_int), parameter :: nintnon_wcns5cv  = 3
    integer(kind_int), parameter :: nintnon_wcns7pv  = 4
    integer(kind_int), parameter :: nintnon_wcns7cv  = 5
    integer(kind_int), parameter :: nintnon_hdcs5ei  = 6
    integer(kind_int), parameter :: nintnon_hdcs7ci  = 7
    integer(kind_int), parameter :: nintnon_hdcs5ci  = 8
    integer(kind_int), parameter :: nintnon_scsl3ci  = 9
    integer(kind_int), parameter :: nintnon_scsl5ci  = 10
    integer(kind_int), parameter :: nintnon_scsh3ci  = 11
    integer(kind_int), parameter :: nintnon_scsh5ci  = 12
    integer(kind_int), parameter :: nintnon_scsn2ci  = 13
    integer(kind_int), parameter :: nintnon_scsn3ci  = 14
    integer(kind_int), parameter :: nintnon_scsn4ci  = 15
    integer(kind_int), parameter :: nintnon_scsh3pi  = 16
    integer(kind_int), parameter :: nintnon_scsh5pi  = 17
    integer(kind_int), parameter :: nintnon_scsn2pi  = 18
    integer(kind_int), parameter :: nintnon_scsn3pi  = 19
    integer(kind_int), parameter :: nintnon_scsn4pi  = 20
    integer(kind_int), parameter :: nintnon_dcsh5ci  = 21
    integer(kind_int), parameter :: nintnon_dcsh7ci  = 22
    integer(kind_int), parameter :: nintnon_dcsh5pi  = 23
    integer(kind_int), parameter :: nintnon_dcsh7pi  = 24

    integer(kind_int), parameter :: nflux_steger   = 1
    integer(kind_int), parameter :: nflux_sw_mod   = 2
    integer(kind_int), parameter :: nflux_vanleer  = 3
    integer(kind_int), parameter :: nflux_roe      = 4
    integer(kind_int), parameter :: nflux_roe_prec = 5
    integer(kind_int), parameter :: nflux_slau     = 6
    integer(kind_int), parameter :: nflux_scmp     = 7

    integer(kind_int), parameter :: nlimit_none      = 0
    integer(kind_int), parameter :: nlimit_minmod    = 1
    integer(kind_int), parameter :: nlimit_vanleer   = 2
    integer(kind_int), parameter :: nlimit_vanalbada = 3

    integer(kind_int), parameter :: nscheme_id_max = 10
    integer(kind_int), parameter :: nscheme_con_d1der = 1
    integer(kind_int), parameter :: nscheme_con_d1int = 2
    integer(kind_int), parameter :: nscheme_vis_d2der = 3
    integer(kind_int), parameter :: nscheme_vis_d2int = 4
    integer(kind_int), parameter :: nscheme_vis_d3der = 5
    integer(kind_int), parameter :: nscheme_vis_d3int = 6
    integer(kind_int), parameter :: nscheme_grd_d2der = 7
    integer(kind_int), parameter :: nscheme_grd_d2int = 8
    integer(kind_int), parameter :: nscheme_grd_d3der = 9
    integer(kind_int), parameter :: nscheme_grd_d3int = 10

    integer(kind_int), parameter :: nscheme_policy_max = 14
    integer(kind_int), parameter :: nscheme_policy_edge2  = 1
    integer(kind_int), parameter :: nscheme_policy_edge4  = 2
    integer(kind_int), parameter :: nscheme_policy_edge6  = 3
    integer(kind_int), parameter :: nscheme_policy_ehen4  = 4
    integer(kind_int), parameter :: nscheme_policy_ehen6  = 5    
    integer(kind_int), parameter :: nscheme_policy_ehcs6  = 6
    integer(kind_int), parameter :: nscheme_policy_ehen8  = 7
    integer(kind_int), parameter :: nscheme_policy_scsl4  = 8
    integer(kind_int), parameter :: nscheme_policy_scsl6  = 9
    integer(kind_int), parameter :: nscheme_policy_ehen6e = 10
    integer(kind_int), parameter :: nscheme_policy_ehcs6e = 11    
    integer(kind_int), parameter :: nscheme_policy_ehen8e = 12    
    integer(kind_int), parameter :: nscheme_policy_scsl4e = 13    
    integer(kind_int), parameter :: nscheme_policy_scsl6e = 14
    
    integer(kind_int), parameter :: nscheme_policys(nscheme_id_max,nscheme_policy_max) = &
        reshape((/ nderive_edge2, nintplt_node2, &                   ! nscheme_policy = nscheme_policy_edge2
                   nderive_edge2, nintplt_node2, nderive_edge2, nintplt_node2, &
                   nderive_edge2, nintplt_node2, nderive_edge2, nintplt_node2, &
                   nderive_edge4, nintplt_node4, &                   ! nscheme_policy = nscheme_policy_edge4
                   nderive_edge4, nintplt_node4, nderive_edge4, nintplt_node4, &
                   nderive_edge4, nintplt_node4, nderive_edge4, nintplt_node4, &
                   nderive_edge6, nintplt_node6, &                   ! nscheme_policy = nscheme_policy_edge6
                   nderive_edge6, nintplt_node6, nderive_edge6, nintplt_node6, &
                   nderive_edge6, nintplt_node6, nderive_edge6, nintplt_node6, &
                   nderive_ehen4, nintplt_node4, &                   ! nscheme_policy = nscheme_policy_ehen4
                   nderive_ehen4, nintplt_node4, nderive_ehen4, nintplt_node4, &
                   nderive_ehen4, nintplt_node4, nderive_ehen4, nintplt_node4, &
                   nderive_ehen6, nintplt_node6, &                   ! nscheme_policy = nscheme_policy_ehen6
                   nderive_ehen6, nintplt_node6, nderive_ehen6, nintplt_node6, &
                   nderive_ehen6, nintplt_node6, nderive_ehen6, nintplt_node6, &
                   nderive_ehcs6, nintplt_node6, &                   ! nscheme_policy = nscheme_policy_ehcs6
                   nderive_ehcs6, nintplt_node6, nderive_ehcs6, nintplt_node6, &
                   nderive_ehcs6, nintplt_node6, nderive_ehcs6, nintplt_node6, &
                   nderive_ehen8, nintplt_node8, &                   ! nscheme_policy = nscheme_policy_ehen8
                   nderive_ehen8, nintplt_node8, nderive_ehen8, nintplt_node8, &
                   nderive_ehen8, nintplt_node8, nderive_ehen8, nintplt_node8, &
                   nderive_scsl4, nintplt_scsl4, &                   ! nscheme_policy = nscheme_policy_scsl4
                   nderive_scsl4, nintplt_scsl4, nderive_scsl4, nintplt_scsl4, &
                   nderive_scsl4, nintplt_scsl4, nderive_scsl4, nintplt_scsl4, &
                   nderive_scsl6, nintplt_scsl6, &                   ! nscheme_policy = nscheme_policy_scsl6
                   nderive_scsl6, nintplt_scsl6, nderive_scsl6, nintplt_scsl6, &
                   nderive_scsl6, nintplt_scsl6, nderive_scsl6, nintplt_scsl6, &
                   nderive_ehen6e, nintplt_node6e, &                 ! nscheme_policy = nscheme_policy_ehen6e
                   nderive_ehen6e, nintplt_node6e, nderive_ehen6e, nintplt_node6e, &
                   nderive_ehen6e, nintplt_node6e, nderive_ehen6e, nintplt_node6e, &
                   nderive_ehcs6e, nintplt_node6e, &                 ! nscheme_policy = nscheme_policy_ehcs6e
                   nderive_ehcs6e, nintplt_node6e, nderive_ehcs6e, nintplt_node6e, &
                   nderive_ehcs6e, nintplt_node6e, nderive_ehcs6e, nintplt_node6e, &
                   nderive_ehen8e, nintplt_node8e, &                 ! nscheme_policy = nscheme_policy_ehen8e
                   nderive_ehen8e, nintplt_node8e, nderive_ehen8e, nintplt_node8e, &
                   nderive_ehen8e, nintplt_node8e, nderive_ehen8e, nintplt_node8e, &
                   nderive_scsl4e, nintplt_scsl4e, &                 ! nscheme_policy = nscheme_policy_scsl4e
                   nderive_scsl4e, nintplt_scsl4e, nderive_scsl4e, nintplt_scsl4e, &
                   nderive_scsl4e, nintplt_scsl4e, nderive_scsl4e, nintplt_scsl4e, &
                   nderive_scsl6e, nintplt_scsl6e, &                 ! nscheme_policy = nscheme_policy_scsl6e
                   nderive_scsl6e, nintplt_scsl6e, nderive_scsl6e, nintplt_scsl6e, &
                   nderive_scsl6e, nintplt_scsl6e, nderive_scsl6e, nintplt_scsl6e /),(/nscheme_id_max,nscheme_policy_max/))

    integer(kind_int), parameter :: nunst_tst_local   = 0
    integer(kind_int), parameter :: nunst_tst_global  = 1
    integer(kind_int), parameter :: nunst_tst_localw  = 2
    integer(kind_int), parameter :: nunst_tst_localm  = 3
    integer(kind_int), parameter :: nunst_tst_globald = 4
    integer(kind_int), parameter :: nunst_tst_localq  = 5    

    integer(kind_int), parameter :: nlhs_rkutta_3step  = 1
    integer(kind_int), parameter :: nlhs_lusgs_std_sca = 2
    integer(kind_int), parameter :: nlhs_prsgs_ust_sca = 3
    integer(kind_int), parameter :: nlhs_prsgs_std_sca = 4
    integer(kind_int), parameter :: nlhs_prsgs_ust_mat = 5
    integer(kind_int), parameter :: nlhs_prsgs_std_mat = 6
    !todo 预处理相关
    integer(kind_int), parameter :: nlhs_lusgs_std_prec = 11     !< 预处理LUSGS标识
    !todo SCMP
    integer(kind_int), parameter :: nlhs_lusgs_std_scmp = 12     !< 预处理LUSGS-SCMP标识
    !todo SCMP
    integer(kind_int), parameter :: nlhs_prsgs_ust_sca_scmp = 13 !< 预处理非定常双时间步SCMP
    
    integer(kind_int), parameter :: nturlhs_lusgs_std_sca = 1
    integer(kind_int), parameter :: nturlhs_prsgs_std_sca = 2

    integer(kind_int), parameter :: ntimeadv_steady = 0
    integer(kind_int), parameter :: ntimeadv_unstdy = 1

    integer(kind_int), parameter :: nsweep_backward = 0
    integer(kind_int), parameter :: nsweep_foreward = 1

    integer(kind_int), parameter :: nupd_field_no  = 0
    integer(kind_int), parameter :: nupd_field_yes = 1

    integer(kind_int), parameter :: nsgl_aver_art = 0   ! 算术奇点平均
    integer(kind_int), parameter :: nsgl_aver_vol = 1   ! 体积奇点平均

    integer(kind_int), parameter :: nsgl_buffer_max  = 5
    integer(kind_int), parameter :: nsgl_buffer_pvs  = 1   ! 使用静态内存1（用于平均流场pv通讯）
    integer(kind_int), parameter :: nsgl_buffer_dpv  = 2   ! 使用静态内存1（用于平均流场pv通讯）
    integer(kind_int), parameter :: nsgl_buffer_qts  = 3   ! 使用静态内存2（用于湍流流场qvst通讯）
    integer(kind_int), parameter :: nsgl_buffer_dtur = 4   ! 使用静态内存2（用于湍流流场qvst通讯）
    integer(kind_int), parameter :: nsgl_buffer_dqt  = 5   ! 使用静态内存3（用于湍流流场dqt通讯）

    integer(kind_int), parameter :: nbc_inter_buf_max  = 6
    integer(kind_int), parameter :: nbc_inter_buf_dyn  = 0   ! 使用动态内存 （对接边界信息通讯）
    integer(kind_int), parameter :: nbc_inter_buf_pvs  = 1   ! 使用静态内存1（用于pv通讯）
    integer(kind_int), parameter :: nbc_inter_buf_dqc  = 2   ! 使用静态内存2（用于rhs通讯）
    integer(kind_int), parameter :: nbc_inter_buf_dpv  = 3   ! 使用静态内存3（用于uvwt梯度通讯）
    integer(kind_int), parameter :: nbc_inter_buf_qts  = 4   ! 使用静态内存4（用于湍流变量通讯）
    integer(kind_int), parameter :: nbc_inter_buf_dqt  = 5   ! 使用静态内存5（用于湍流右短项通讯）
    integer(kind_int), parameter :: nbc_inter_buf_dtur = 6   ! 使用静态内存6（用于湍流梯度通讯）

    integer(kind_int), parameter :: nbclistmax = 7
    integer(kind_int), parameter :: bc_cut1to1  = -1    
    integer(kind_int), parameter :: bc_wall     = 2
    integer(kind_int), parameter :: bc_symmetry = 3
    integer(kind_int), parameter :: bc_farfield = 4
    integer(kind_int), parameter :: bc_inflow   = 5
    integer(kind_int), parameter :: bc_outflow  = 6
    integer(kind_int), parameter :: bc_pole     = 7
    integer(kind_int), parameter :: bclists(nbclistmax) = &
        (/ bc_cut1to1, bc_farfield, bc_inflow  , &
           bc_outflow, bc_wall    , bc_symmetry, &
           bc_pole /)

    integer(kind_int), parameter :: subc_none = 0

    integer(kind_int), parameter :: bc_cut1to1A       = -11
    integer(kind_int), parameter :: bc_cut1to1_split  = -2
    integer(kind_int), parameter :: bc_cut1to1_splitA = -21
    integer(kind_int), parameter :: bc_patched        = 8
    integer(kind_int), parameter :: bc_patchedA       = 81
    integer(kind_int), parameter :: bc_wallA          = 21
    integer(kind_int), parameter :: subc_cut_original = 0
    integer(kind_int), parameter :: subc_cut_spliting = 1
    integer(kind_int), parameter :: subc_cut_patched  = 2

    integer(kind_int), parameter :: bc_wall_adiabatic  = 1
    integer(kind_int), parameter :: bc_wall_isothermal = 2
    integer(kind_int), parameter :: bc_wall_slip       = 3

    integer(kind_int), parameter :: bc_symmetry_point = 1
    integer(kind_int), parameter :: bc_symmetry_plane = 2

    integer(kind_int), parameter :: bc_farfield_riemann = 1
    integer(kind_int), parameter :: bc_farfield_charact = 2

    integer(kind_int), parameter :: bc_pole_i = 71
    integer(kind_int), parameter :: bc_pole_j = 72
    integer(kind_int), parameter :: bc_pole_k = 73

    integer(kind_int), parameter :: subc_pole_i = 1
    integer(kind_int), parameter :: subc_pole_j = 2
    integer(kind_int), parameter :: subc_pole_k = 3

    integer(kind_int), parameter :: blk_mark_dstmin = 1
    integer(kind_int), parameter :: blk_mark_sponge = 2

    integer(kind_int), parameter :: bc_mark_dstmin = 1
    integer(kind_int), parameter :: bc_mark_sponge = 2

end module mod_constants
