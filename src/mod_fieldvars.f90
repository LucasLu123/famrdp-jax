
module mod_fieldvars
    use mod_kndconsts, only : kind_int,kind_real
    use mod_datatypes, only : top_block_t,var_block_t
    use mod_datatypes, only : bc_inter_t,block_compute_t
    use mod_datatypes, only : wll_point_t,patched_t
    implicit none

    !------------------网格拓扑信息变量声明--------------------
    integer(kind_int)          :: nblocks          !< 网格块数
    integer(kind_int)          :: nprocs           !< 网格剖分数
    type(top_block_t), pointer :: mb_top(:)        !< 网格拓扑信息
    type(top_block_t), pointer :: mb_topc(:)       !< 双网格拓扑信息
    type(top_block_t), pointer :: mb_topsp(:)      !< 双网格求解点拓扑信息

    integer(kind_int)          :: ninters          !< 对接+拼接边界的个数
    integer(kind_int)          :: ninterf          !< 对接边界的个数
    integer(kind_int)          :: nmsg_intf        !< 对接边界的通讯个数
    type(bc_inter_t ), pointer :: inters(:)        !< 对接+拼接边界的信息
    type(bc_inter_t ), pointer :: interscc(:)      !< 双网格对接边界的信息
    type(bc_inter_t ), pointer :: interssp(:)      !< 双网格求解点对接边界的信息
    type(patched_t),   pointer :: patched(:)       !< 拼接边界的插值坐标
    real(kind_real),    pointer :: scal(:,:)          !< 拼接边界的缩放比例
    real(kind_real),    pointer :: scalt(:,:)         !< 目标拼接边界的缩放比例
    integer(kind_int), pointer :: bc_seq(:,:)
    real(kind_real),    pointer :: surf_pxyz(:,:,:)
    real(kind_real),    pointer :: surf_check(:)
    integer(kind_int), pointer :: s_dir_all(:,:)
    real(kind_real),    pointer :: L_check_all(:,:)
    real(kind_real),    pointer :: P_check_all(:,:)

    integer(kind_int)              :: nblkcoms     !< 计算块个数
    type(block_compute_t), pointer :: blkcoms(:)   !< 计算块信息
    type(block_compute_t), pointer :: blkcomscc(:) !< 双网格边界点计算块信息
    type(block_compute_t), pointer :: blkcomssp(:) !< 双网格求解点计算块信息
    integer(kind_int)              :: ntotpts      !< 计算块总点数
    integer(kind_int)              :: ntotptscc    !< 双网格边界点计算块总点数
    integer(kind_int)              :: ntotptssp    !< 双网格求解点计算块总点数

    type(var_block_t), pointer :: mb_xyz(:)        !< 网格坐标
    type(var_block_t), pointer :: mb_xyzcc(:)      !< 双网格边界点网格坐标
    type(var_block_t), pointer :: mb_xyzsp(:)      !< 双网格求解点网格坐标

    type(var_block_t), pointer :: mb_fsf(:,:)      !< 块边界标记
    type(var_block_t), pointer :: mb_fsfcc(:,:)    !< 双网格边界点块边界标记
    type(var_block_t), pointer :: mb_fsfsp(:,:)    !< 双网格求解点块边界标记
    type(var_block_t), pointer :: mb_fsffp(:,:)    !< 双网格通量点块边界标记

    type(var_block_t), pointer :: mb_vol(:)        !< 网格雅可比
    type(var_block_t), pointer :: mb_sxyz(:)       !< 网格导数
    type(var_block_t), pointer :: mb_volsp(:)      !< 双网格求解点网格雅可比
    type(var_block_t), pointer :: mb_sxyzsp(:)     !< 双网格求解点网格导数   
    type(var_block_t), pointer :: mb_volcc(:)      !< 双网格边界点网格雅可比
    type(var_block_t), pointer :: mb_sxyzcc(:)     !< 双网格边界点网格导数     


    !--------------------NS方程变量声明----------------------
!#define PARAM_OPT     !CHEyg  参数化优化  在makefile中定义

#ifndef PARAM_OPT
    integer(kind_int)  :: neqn                    !< 守恒变量个数
    integer(kind_int)  :: npvs                     !< 原始变量个数	
#else	
	integer(kind_int),parameter :: neqn = 5          !< 守恒变量个数
    integer(kind_int),parameter :: npvs = neqn     !< 原始变量个数	
#endif

    type(var_block_t), pointer :: mb_qc(:)         !< 守恒变量
    type(var_block_t), pointer :: mb_pv(:)         !< 原始变量
    type(var_block_t), pointer :: mb_pvfp(:)       !< 双网格中通量点上原始变量
    type(var_block_t), pointer :: mb_ws(:)         !< 涡熵变量

    type(var_block_t), pointer :: mb_t(:)          !< 静温
    type(var_block_t), pointer :: mb_c(:)          !< 声速

    type(var_block_t), pointer :: mb_rhs(:)        !< 右端项

    type(var_block_t), pointer :: mb_dpv(:)        !< 速度、温度梯度
    !type(var_block_t), pointer :: mb_dpvfp(:)      !< 双网格中通量点上速度、温度梯度
    type(var_block_t), pointer :: mb_dpvfp(:,:)    !< 双网格通量点上速度、温度梯度
    type(var_block_t), pointer :: mb_dws(:)        !< 涡量、熵梯度

    type(var_block_t), pointer :: mb_vsl(:)        !< 分子运动粘性系数
    type(var_block_t), pointer :: mb_vst(:)        !< 湍流粘性系数
    
    !todo 预处理相关
    !type(var_block_t), pointer :: mb_prec(:)       !< 预处理矩阵5*5(转成一维数组 先行后列)
    !type(var_block_t), pointer :: mb_prec_inv(:)   !< 预处理矩阵的逆矩阵5*5(转成一维数组 先行后列)
    !type(var_block_t), pointer :: mb_nonprec(:)    !< 守恒原始变量转换矩阵5*5(转成一维数组 先行后列)
    type(var_block_t), pointer :: mb_src_prec(:)   !< 考虑预处理后的对流项谱半径
    !type(var_block_t), pointer :: mb_prec_diag_inv(:) !<考虑预处理后的左端项D矩阵的逆矩阵
    
    type(var_block_t), pointer :: mb_src(:)        !< 对流项谱半径
    type(var_block_t), pointer :: mb_srv(:)        !< 粘性项谱半径

    type(var_block_t), pointer :: mb_rdt(:)        !< 谱半径之和

    type(var_block_t), pointer :: mb_dt(:)         !< 时间步

    type(var_block_t), pointer :: mb_dq(:)         !< 守恒变量增量

    type(var_block_t), pointer :: mb_q0(:)         !< 守恒变量（rkutta,prsgs_std）

    type(var_block_t), pointer :: mb_qm(:)         !< 守恒变量（prsgs_std）

    type(var_block_t), pointer :: mb_aiv(:)        !< 对角逆矩阵（prsgs_std，prsgs_blu）


    !----------------------湍流方程变量声明-----------------------
    integer(kind_int)          :: neqt             !< 湍流方程个数
    type(var_block_t), pointer :: mb_qt(:)         !< 湍流变量

    integer(kind_int)          :: npvt             !< 湍流复合变量个数（4[ro,u,v,w]+neqt[k,w,...]）
    type(var_block_t), pointer :: mb_pvt(:)        !< 湍流复合变量(以指针形式存在)

    integer(kind_int)          :: nqvst            !< 湍流复合变量个数（neqt[k,w,...]+1[vst]）
    type(var_block_t), pointer :: mb_qvst(:)       !< 湍流复合变量(以指针形式存在)

    type(var_block_t), pointer :: mb_rhst(:)       !< 湍流右端项
    type(var_block_t), pointer :: mb_dqt(:)        !< 湍流变量增量
    type(var_block_t), pointer :: mb_dtur(:)       !< 湍流变量导数

    type(var_block_t), pointer :: mb_srt(:)        !< 湍流粘性项谱半径

    type(var_block_t), pointer :: mb_rtur(:)       !< 湍流谱半径之和

    type(var_block_t), pointer :: mb_dst(:)        !< 壁面最小距离
    type(var_block_t), pointer :: mb_bld(:)        !< SST模型混合函数

    integer(kind_int)          :: nwallpnts        !< 壁面点个数
    type(wll_point_t), pointer :: wallpnts(:)      !< 壁面点信息

    !----------------------气动声学和统计变量声明-----------------------
    type(var_block_t), pointer :: mb_dsp(:)        !< 吸波层界面距离

    integer(kind_int)          :: nstepmean        !< 平均流场统计步数
    integer(kind_int)          :: nsteprms         !< 脉动量统计步数
    integer(kind_int)          :: ntime            !< 声学量统计步数
    integer(kind_int)          :: npoints          !< 声学积分面的点数
    type(var_block_t), pointer :: mb_fmean(:)      !< 流场统计平均量
    type(var_block_t), pointer :: mb_frms(:)       !< 脉动量统计平均

end module mod_fieldvars
