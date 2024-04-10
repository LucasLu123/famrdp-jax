
subroutine sol_prsgs_ust_sca
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,nbc_inter_buf_dqc
    use mod_constants, only : nsgl_aver_art,nsgl_buffer_pvs
    use mod_constants, only : ntimeadv_steady,ntimeadv_unstdy
    use mod_constants, only : nsweep_foreward,nsweep_backward,nprec_non
    use mod_variables, only : nsubmax,nghnode,niter,nprec
    use mod_fieldvars, only : nblkcoms,blkcoms
    use mod_fieldvars, only : neqn,mb_dq,mb_qc,mb_q0,mb_qm
    use mod_interface, only : pre_exchange_bc_var,post_exchange_bc_var,average_bc_var
    use mod_interface, only : exchange_singulars,average_singulars
    use mod_interface, only : assign_mb_var_uniform,assign_bc_var_nb
    use mod_interface, only : calc_mb_var_via_sub
    use mod_interface, only : assign_var_via_com_nb
    implicit none
    integer(kind_int)          :: iter,nc,nb,nstop,itersgs
    real(kind_real), parameter :: cpr(3)=(/1.5,0.5,1.5/)
    real(kind_real)            :: dq0(1:neqn)
    external                   :: v1_eq_v2

    dq0(:) = zero
    call assign_mb_var_uniform(mb_dq,1,neqn,nghnode,dq0)

    call calc_mb_var_via_sub(mb_qc,1,neqn,v1_eq_v2,mb_qm,1,neqn,nghnode)

    do iter=1,nsubmax

        niter = iter
        
        call prepare_linear_system

        !do itersgs=1,2
        do nc=1,nblkcoms
            nb = blkcoms(nc)%nb

            call prsgs_ust_sca_sweep(nb,nsweep_foreward,cpr,ntimeadv_unstdy)
            call assign_bc_var_nb(nb,mb_dq,1,neqn,0,dq0)

            call prsgs_ust_sca_sweep(nb,nsweep_backward,cpr,ntimeadv_unstdy)
            call assign_bc_var_nb(nb,mb_dq,1,neqn,0,dq0)

            call assign_var_via_com_nb(nb,mb_dq,1,neqn)
        end do
        !end do

! 0611
        call pre_exchange_bc_var(mb_dq,1,neqn,1,nbc_inter_buf_dqc,nsgl_aver_art)
        call exchange_singulars(mb_dq,1,neqn,nsgl_buffer_pvs,nsgl_aver_art)
        call post_exchange_bc_var(mb_dq,1,neqn,1,nbc_inter_buf_dqc,nsgl_aver_art)
        call average_bc_var(mb_dq,1,neqn,nbc_inter_buf_dqc,nsgl_aver_art)
        call average_singulars(mb_dq,1,neqn,nsgl_buffer_pvs,nsgl_aver_art)

        call residual(ntimeadv_unstdy)
        call stop_subiter(iter,nstop)

        if (nprec > nprec_non) then
            call update_prec
        else
            call update
        end if

        if (nstop /= 0) exit
    end do

    call calc_mb_var_via_sub(mb_qm,1,neqn,v1_eq_v2,mb_q0,1,neqn,nghnode)

end subroutine sol_prsgs_ust_sca

subroutine sol_prsgs_ust_sca_sp
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero
    use mod_constants, only : ntimeadv_steady,ntimeadv_unstdy
    use mod_constants, only : nsweep_foreward,nsweep_backward,nprec_non
    use mod_variables, only : nsubmax,nghnode,niter,nprec
    use mod_fieldvars, only : nblkcoms,blkcomssp
    use mod_fieldvars, only : neqn,mb_dq,mb_qc,mb_q0,mb_qm
    use mod_interface, only : assign_mb_var_uniform_sp
    use mod_interface, only : calc_mb_var_via_sub_sp
    use mod_interface, only : assign_var_via_com_nb_sp
    implicit none
    integer(kind_int)          :: iter,nc,nb,nstop,itersgs
    real(kind_real), parameter :: cpr(3)=(/1.5,0.5,1.5/)
    real(kind_real)            :: dq0(1:neqn)
    external                   :: v1_eq_v2

    dq0(:) = zero
    call assign_mb_var_uniform_sp(mb_dq,1,neqn,nghnode,dq0)

    call calc_mb_var_via_sub_sp(mb_qc,1,neqn,v1_eq_v2,mb_qm,1,neqn,nghnode)

    do iter=1,nsubmax

        niter = iter
        
        call prepare_linear_system_sp

        !do itersgs=1,2
        do nc=1,nblkcoms
            nb = blkcomssp(nc)%nb

            call prsgs_ust_sca_sweep_sp(nb,nsweep_foreward,cpr,ntimeadv_unstdy)
            !call assign_bc_var_nb(nb,mb_dq,1,neqn,0,dq0)

            call prsgs_ust_sca_sweep_sp(nb,nsweep_backward,cpr,ntimeadv_unstdy)
            !call assign_bc_var_nb(nb,mb_dq,1,neqn,0,dq0)

            call assign_var_via_com_nb_sp(nb,mb_dq,1,neqn)
        end do
        !end do

        call residual_sp(ntimeadv_unstdy)
        call stop_subiter_sp(iter,nstop)

        if (nprec > nprec_non) then
            call update_prec_sp
        else
            call update_sp
        end if

        if (nstop /= 0) exit
    end do

    call calc_mb_var_via_sub_sp(mb_qm,1,neqn,v1_eq_v2,mb_q0,1,neqn,nghnode)

end subroutine sol_prsgs_ust_sca_sp

subroutine sol_prsgs_std_sca
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero
    use mod_constants, only : nbc_inter_buf_dqc
    use mod_constants, only : nsgl_aver_art,nsgl_buffer_pvs
    use mod_constants, only : ntimeadv_steady,ntimeadv_unstdy
    use mod_constants, only : nsweep_foreward,nsweep_backward
    use mod_variables, only : nghnode,nsubmax
    use mod_fieldvars, only : nblkcoms,blkcoms,neqn,mb_dq
    use mod_interface, only : pre_exchange_bc_var,post_exchange_bc_var,average_bc_var
    use mod_interface, only : exchange_singulars,average_singulars
    use mod_interface, only : assign_mb_var_uniform,assign_bc_var_nb
    use mod_interface, only : assign_var_via_com_nb
    implicit none
    integer(kind_int) :: iter,nc,nb,nstop
    real(kind_real), parameter :: cpr(3)=(/0.0,0.0,1.0/)
    real(kind_real)   :: dq0(1:neqn)

    call prepare_linear_system

    dq0(:) = zero
    call assign_mb_var_uniform(mb_dq,1,neqn,nghnode,dq0)

    do iter=1,nsubmax
        do nc=1,nblkcoms
            nb = blkcoms(nc)%nb

            call prsgs_ust_sca_sweep(nb,nsweep_foreward,cpr,ntimeadv_steady)
            call assign_bc_var_nb(nb,mb_dq,1,neqn,0,dq0)

            call prsgs_ust_sca_sweep(nb,nsweep_backward,cpr,ntimeadv_steady)
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

end subroutine sol_prsgs_std_sca

subroutine sol_prsgs_std_sca_sp
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero
    use mod_constants, only : ntimeadv_steady,ntimeadv_unstdy
    use mod_constants, only : nsweep_foreward,nsweep_backward
    use mod_variables, only : nghnode,nsubmax
    use mod_fieldvars, only : nblkcoms,blkcomssp,neqn,mb_dq
    use mod_interface, only : assign_mb_var_uniform_sp
    use mod_interface, only : assign_var_via_com_nb_sp
    implicit none
    integer(kind_int) :: iter,nc,nb,nstop
    real(kind_real), parameter :: cpr(3)=(/0.0,0.0,1.0/)
    real(kind_real)   :: dq0(1:neqn)

    call prepare_linear_system_sp

    dq0(:) = zero
    call assign_mb_var_uniform_sp(mb_dq,1,neqn,nghnode,dq0)

    do iter=1,nsubmax
        do nc=1,nblkcoms
            nb = blkcomssp(nc)%nb

            call prsgs_ust_sca_sweep_sp(nb,nsweep_foreward,cpr,ntimeadv_steady)
            !call assign_bc_var_nb(nb,mb_dq,1,neqn,0,dq0)

            call prsgs_ust_sca_sweep_sp(nb,nsweep_backward,cpr,ntimeadv_steady)
            !call assign_bc_var_nb(nb,mb_dq,1,neqn,0,dq0)

            call assign_var_via_com_nb_sp(nb,mb_dq,1,neqn)
        end do

        call residual_sp(ntimeadv_steady)
        call stop_subiter_sp(iter,nstop)

        if (nstop /= 0) exit
    end do

    call update_sp

end subroutine sol_prsgs_std_sca_sp

subroutine prsgs_ust_sca_sweep(nb,swp,cpr,nsw)
!   CHEyg
!#define MemOPT1	! 优化sxyz访问   在makefile中定义

!#define MemOPT2	! 优化rhs,pv,qm,q0,qc访问

! #define MemOPT3	! 优化src和srv访问

#ifdef MemOPT1	
    use mod_fieldvars, only : nblkcoms, blkcoms    ! CHEyg 本地块数和块指针
    use mod_opt_vars,  only : dim4_fld_array_t, sxyz1
#endif
	
	
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,one,half
    use mod_constants, only : nvis_euler,nsw_dir_close,nprec_non,nlhs_prsgs_ust_sca
    use mod_datatypes, only : fld_array_t
    use mod_variables, only : nvis,dtime,nsw_kdir,nprec,nlhs,fsw_kdir
    use mod_fieldvars, only : mb_top,mb_vol,mb_sxyz
    !use mod_fieldvars, only : mb_nonprec,mb_prec_diag_inv
    use mod_fieldvars, only : mb_pv,npvs,mb_qc,neqn
    use mod_fieldvars, only : mb_dq,mb_rhs,mb_dt
    use mod_fieldvars, only : mb_src,mb_srv,mb_rdt
    use mod_fieldvars, only : mb_vsl,mb_vst
    use mod_fieldvars, only : mb_q0,mb_qm
    use mod_openmp
    implicit none
    integer(kind_int), intent(in) :: nb,swp
    real(kind_real),   intent(in) :: cpr(3)
    integer(kind_int), intent(in) :: nsw
    integer(kind_int)          :: i,j,k,i1,j1,k1,m,ierr,l,oned_index,nsp,nep
    integer(kind_int)          :: st(3),ed(3),dijk
    real(kind_real)            :: srci,srcj,srck,srvi,srvj,srvk
    real(kind_real)            :: kt,kx,ky,kz,et,ex,ey,ez,ct,cx,cy,cz
    real(kind_real)            :: dia,odia,odt,dq1,dq2,rhs0(1:neqn)
    real(kind_real)            :: pvi(1:npvs),pvj(1:npvs),pvk(1:npvs)
    real(kind_real)            :: dqi(1:neqn),dqj(1:neqn),dqk(1:neqn)
    real(kind_real)            :: dUi(1:neqn),dUj(1:neqn),dUk(1:neqn)
    real(kind_real)            :: dfi(1:neqn),dfj(1:neqn),dfk(1:neqn)
    real(kind_real)            :: point_nonprec_i(1:4),point_nonprec_j(1:4),point_nonprec_k(1:4)
    type(fld_array_t), pointer :: vol(:),sxyz(:),pv(:),nonprec(:),diagInv(:)
    type(fld_array_t), pointer :: dt(:),dq(:),rhs(:)
    type(fld_array_t), pointer :: src(:),srv(:),rdt(:),vsl(:),vst(:)
    type(fld_array_t), pointer :: qc(:),q0(:),qm(:)
    real(kind_real)            :: diagInvMatrix(1:25),prim(1:5),sn(1:3,1:3),eigenValue(1:3,1:5)
    real(kind_real)            :: tmp,rca,rcb,rcc,sumSpectralRadius,visl,vist,vis,length,diagElement
#ifdef OMP_IMP
    integer(kind_int)          :: idt, numt, nflag(0:omp_max_num_threads)
#endif

! CHEyg
#ifdef MemOPT1	
    integer(kind_int)               :: nc,nc1    !局部块编号
    type(dim4_fld_array_t), pointer :: sxyz2
#endif

#ifdef MemOPT2    
    real(kind_real), allocatable,dimension(:,:,:,:) :: qc_r4d, q0_r4d, rhs_r4d, qm_r4d ,pv_r4d
    integer(kind_int)                               :: nerror(5), i2,j2,k2

#endif

#ifdef MemOPT3    
    real(kind_real), allocatable,dimension(:,:,:,:) :: srv_r4d, src_r4d
#endif

    vol  => mb_vol(nb)%fld
    sxyz => mb_sxyz(nb)%fld
    pv   => mb_pv(nb)%fld
    qc   => mb_qc(nb)%fld
    dt   => mb_dt(nb)%fld
    dq   => mb_dq(nb)%fld
    rhs  => mb_rhs(nb)%fld
    src  => mb_src(nb)%fld
    rdt  => mb_rdt(nb)%fld

    if (nvis > nvis_euler) then
        srv  => mb_srv(nb)%fld
        vsl => mb_vsl(nb)%fld
        vst => mb_vst(nb)%fld
        srv => mb_srv(nb)%fld
    end if
    
    !if (nprec > nprec_non) then
    !    nonprec => mb_nonprec(nb)%fld     !<守恒变量对原始变量的转换矩阵P
    !    diagInv => mb_prec_diag_inv(nb)%fld   !<对角矩阵的逆矩阵
    !end if

    if (nsw > 0) then
        q0 => mb_q0(nb)%fld !< Q(n-1)
        qm => mb_qm(nb)%fld !< Q(n)
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

! CHEyg
#ifdef MemOPT1	
    ! prsgs_ust_mat_sweep等被调用时给出的 nb 是全局的块编号，需要根据给出的全局 nb 编号算出在本进程上的块号
      do nc1 = 1,nblkcoms
         if(nb .EQ. blkcoms(nc1)%nb) then
            nc = nc1
         end if
      enddo

      sxyz2 => sxyz1(nc)   ! 指向当前要用的块
      
#endif     
	  
#ifdef MemOPT2 
    !给临时数组分配存储空间
      if(dijk .EQ. 1) then
         i1 = st(1) 
         i2 = ed(1) 
         
         j1 = st(2)  
         j2 = ed(2)   
                 
         k1 = st(3)  
         k2 = ed(3)  
      else 
         i2 = st(1) 
         i1 = ed(1) 
         
         j2 = st(2)  
         j1 = ed(2)   
                 
         k2 = st(3)  
         k1 = ed(3)        
      endif

       if (nsw > 0) then
          allocate( q0_r4d(neqn,  i1:i2,  j1:j2, k1:k2), stat= nerror(1) )                        
          allocate( qm_r4d(neqn,  i1:i2,  j1:j2, k1:k2), stat= nerror(2) ) 
       
          do i = 1, 2
             if(nerror(i) .NE. 0) then
                write(*,*) "In prsgs_ust_mat_sweep. Allocating memory for X_r4d error. Stop! "
                stop
             endif
          enddo                 
       endif
    
       allocate( qc_r4d(  neqn,  i1:i2, j1:j2, k1:k2), stat= nerror(3) )
       allocate( rhs_r4d( neqn,  i1:i2, j1:j2, k1:k2), stat= nerror(4) )
       allocate( pv_r4d(  neqn,  (i1-1):(i2+1),  (j1-1):(j2+1), (k1-1):(k2+1) ), stat= nerror(5) )    
       ! pv的维数上下界均多一个
       
       do i = 3, 5
          if(nerror(i) .NE. 0) then
             write(*,*) "In prsgs_ust_mat_sweep. Allocating memory for X_r4d error. Stop! "
             stop
          endif
       enddo
       
    !将数据复制到临时数组中。逐个数组的赋值以避免cache冲突
         
       do k = k1,k2
       do j = j1,j2
       do i = i1,i2 
          do m = 1,neqn
             qc_r4d(m,i,j,k) = qc(m)%r3d(i,j,k)
          enddo
       enddo
       enddo
       enddo

   
    if (nsw > 0) then   
       do k = k1,k2
       do j = j1,j2
       do i = i1,i2 
          do m = 1,neqn
             qm_r4d(m,i,j,k) = qm(m)%r3d(i,j,k)
          enddo
       enddo
       enddo
       enddo   
                     
       do k = k1,k2
       do j = j1,j2
       do i = i1,i2 
          do m = 1,neqn
             q0_r4d(m,i,j,k) = q0(m)%r3d(i,j,k)
          enddo
       enddo
       enddo
       enddo

     endif

       do k = k1,k2
       do j = j1,j2
       do i = i1,i2 
          do m = 1,neqn
             rhs_r4d(m,i,j,k) = rhs(m)%r3d(i,j,k)
          enddo
       enddo
       enddo
       enddo
            
       do k = (k1-1),(k2+1)
       do j = (j1-1),(j2+1)
       do i = (i1-1),(i2+1) 
          do m = 1,neqn
             pv_r4d(m,i,j,k) = pv(m)%r3d(i,j,k)
          enddo
       enddo
       enddo
       enddo                           
        
#endif    

#ifdef MemOPT3    
    !给临时数组分配存储空间
      if(dijk .EQ. 1) then
         i1 = st(1) 
         i2 = ed(1) 
         
         j1 = st(2)  
         j2 = ed(2)   
                 
         k1 = st(3)  
         k2 = ed(3)  
      else 
         i2 = st(1) 
         i1 = ed(1) 
         
         j2 = st(2)  
         j1 = ed(2)   
                 
         k2 = st(3)  
         k1 = ed(3)        
      endif

       allocate( src_r4d(  neqn,  (i1-1):(i2+1),  (j1-1):(j2+1), (k1-1):(k2+1) ), stat= nerror(1) )
       allocate( srv_r4d(  neqn,  (i1-1):(i2+1),  (j1-1):(j2+1), (k1-1):(k2+1) ), stat= nerror(2) )

       do k = (k1-1),(k2+1)
       do j = (j1-1),(j2+1)
       do i = (i1-1),(i2+1) 
          do m = 1,3
             src_r4d(m,i,j,k) = src(m)%r3d(i,j,k)
          enddo
       enddo
       enddo
       enddo 

       do k = (k1-1),(k2+1)
       do j = (j1-1),(j2+1)
       do i = (i1-1),(i2+1) 
          do m = 1,3
             srv_r4d(m,i,j,k) = src(m)%r3d(i,j,k)
          enddo
       enddo
       enddo
       enddo 
	   
#endif
	
	

!$OMP parallel private(k,k1,j,j1,i,i1,odt,idt,pvi,pvj,pvk,dqi,dqj,dqk,rhs0,dfi,dfj,dfk,kx,ky,kz,&
!$OMP                  ex,ey,ez,cx,cy,cz,srvi,srvj,srvk,srci,srcj,srck,odia,dq1,dq2,dia)

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
        do m=1,neqn
            
#ifndef MemOPT2          
            rhs0(m) = rhs(m)%r3d(i,j,k)  
#else
            rhs0(m) = rhs_r4d(m, i,j,k)
            
#endif    			
        end do

        dia = rdt(1)%r3d(i,j,k)

        if (nsw > 0) then
            odt = vol(1)%r3d(i,j,k)/dtime
            do m=1,neqn

#ifndef MemOPT2                       
                dq1 = qm(m)%r3d(i,j,k) - q0(m)%r3d(i,j,k)
                dq2 = qc(m)%r3d(i,j,k) - qm(m)%r3d(i,j,k)
#else
                dq1 = qm_r4d(m, i,j,k) - q0_r4d(m, i,j,k)
                dq2 = qc_r4d(m, i,j,k) - qm_r4d(m, i,j,k)

#endif      

                rhs0(m) = rhs0(m) - odt*(cpr(1)*dq2 - cpr(2)*dq1)
            end do

            dia = dia + cpr(3)*odt
        end if

        i1 = i - 1
        j1 = j - 1
        k1 = k - 1

#ifndef MemOPT3    
        srci = src(1)%r3d(i1,j,k)
        srcj = src(2)%r3d(i,j1,k)
        srck = src(3)%r3d(i,j,k1)
#else 
        srci = src_r4d(1, i1,j,k)
        srcj = src_r4d(2, i,j1,k)
        srck = src_r4d(3, i,j,k1)
#endif		
		

        do m=1,npvs

#ifndef MemOPT2        
            pvi(m) = pv(m)%r3d(i1,j,k)
            pvj(m) = pv(m)%r3d(i,j1,k)
            pvk(m) = pv(m)%r3d(i,j,k1)
#else
            pvi(m) = pv_r4d(m, i1,j,k)
            pvj(m) = pv_r4d(m, i,j1,k)
            pvk(m) = pv_r4d(m, i,j,k1)   

#endif      

        end do

        do m=1,neqn
            dqi(m) = dq(m)%r3d(i1,j,k)
            dqj(m) = dq(m)%r3d(i,j1,k)
            dqk(m) = dq(m)%r3d(i,j,k1)
			dUi(m) = dqi(m)
            dUj(m) = dqj(m)
            dUk(m) = dqk(m)
        end do
 
        nsp = 1
        nep = 5
        sn(:,:) = 1.0
        vis = 1.0
        length = 1.0
        
        !todo 原始变量增量转换成守恒变量增量ΔQ=P*ΔW
        if (nprec > nprec_non) then
            call calculate_precondition_matrix_point(nsp,nep,pvi,point_nonprec_i)
            call calculate_precondition_matrix_point(nsp,nep,pvj,point_nonprec_j)
            call calculate_precondition_matrix_point(nsp,nep,pvk,point_nonprec_k)

            !将原始变量增量转换成守恒变量增量ΔQ=P*ΔW
            call primitiveDeltaWToConservativeDeltaQ(point_nonprec_i,dqi,pvi,dUi)
            call primitiveDeltaWToConservativeDeltaQ(point_nonprec_j,dqj,pvj,dUj)
            if (nsw_kdir == nsw_dir_close) then
                dUk(:) = zero
            else
                call primitiveDeltaWToConservativeDeltaQ(point_nonprec_k,dqk,pvk,dUk)
            end if
        end if

#ifndef MemOPT1
        kx = sxyz(1)%r3d(i1,j,k)
        ky = sxyz(2)%r3d(i1,j,k)
        kz = sxyz(3)%r3d(i1,j,k)
#else
        kx = sxyz2%r4d(1,i1,j,k)
        ky = sxyz2%r4d(2,i1,j,k)
        kz = sxyz2%r4d(3,i1,j,k)
#endif
		
        call mxdq_std(1,npvs,pvi,kt,kx,ky,kz,1,neqn,dUi,dfi,srci,one)
        !call mxdq_std(1,npvs,pvi,kt,kx,ky,kz,1,neqn,dqi,dfi,srci,one)


#ifndef MemOPT1        
        ex = sxyz(4)%r3d(i,j1,k)
        ey = sxyz(5)%r3d(i,j1,k)
        ez = sxyz(6)%r3d(i,j1,k)
#else
        ex = sxyz2%r4d(4,i,j1,k)
        ey = sxyz2%r4d(5,i,j1,k)
        ez = sxyz2%r4d(6,i,j1,k)
#endif
		
        call mxdq_std(1,npvs,pvj,et,ex,ey,ez,1,neqn,dUj,dfj,srcj,one)
        !call mxdq_std(1,npvs,pvj,et,ex,ey,ez,1,neqn,dqj,dfj,srcj,one)

        if (nsw_kdir == nsw_dir_close) then
            dfk(:) = zero
        else

#ifndef MemOPT1        
          cx = sxyz(7)%r3d(i,j,k1)
          cy = sxyz(8)%r3d(i,j,k1)
          cz = sxyz(9)%r3d(i,j,k1)
#else
          cx = sxyz2%r4d(7,i,j,k1)
          cy = sxyz2%r4d(8,i,j,k1)
          cz = sxyz2%r4d(9,i,j,k1)
#endif  
		  
            call mxdq_std(1,npvs,pvk,ct,cx,cy,cz,1,neqn,dUk,dfk,srck,one)
            !call mxdq_std(1,npvs,pvk,ct,cx,cy,cz,1,neqn,dqk,dfk,srck,one)
        end if

        rhs0(:) = rhs0(:) + (dfi(:) + dfj(:) + dfk(:))

        if (nvis > nvis_euler) then
		
#ifndef MemOPT3    
            srvi = srv(1)%r3d(i1,j,k)
            srvj = srv(2)%r3d(i,j1,k)
            srvk = srv(3)%r3d(i,j,k1)
#else
			srvi = srv_r4d(1, i1,j,k)
            srvj = srv_r4d(2, i,j1,k)
            srvk = srv_r4d(3, i,j,k1)
#endif
			do m=1,neqn
                rhs0(m) = rhs0(m) + &
                          half*(srvi*dUi(m) + srvj*dUj(m) + srvk*dUk(m))
            end do
            !do m=1,neqn
            !    rhs0(m) = rhs0(m) + &
            !              half*(srvi*dqi(m) + srvj*dqj(m) + srvk*dqk(m))
            !end do
            
        end if


        i1 = i + 1
        j1 = j + 1
        k1 = k + 1

#ifndef MemOPT3 
        srci = src(1)%r3d(i1,j,k)
        srcj = src(2)%r3d(i,j1,k)
        srck = src(3)%r3d(i,j,k1)
#else
        srci = src_r4d(1, i1,j,k)
        srcj = src_r4d(2, i,j1,k)
        srck = src_r4d(3, i,j,k1)
		
#endif 	

        do m=1,npvs

#ifndef MemOPT2
            pvi(m)  = pv(m)%r3d(i1,j,k)
            pvj(m)  = pv(m)%r3d(i,j1,k)
            pvk(m)  = pv(m)%r3d(i,j,k1)
#else
            pvi(m)  = pv_r4d(m, i1,j,k)
            pvj(m)  = pv_r4d(m, i,j1,k)
            pvk(m)  = pv_r4d(m, i,j,k1)
             
#endif
            prim(m) = pv(m)%r3d(i,j,k)
			 
        end do

        do m=1,neqn
            dqi(m) = dq(m)%r3d(i1,j,k)
            dqj(m) = dq(m)%r3d(i,j1,k)
            dqk(m) = dq(m)%r3d(i,j,k1)
            dUi(m) = dqi(m)
            dUj(m) = dqj(m)
            dUk(m) = dqk(m)
        end do

        

        !todo 原始变量增量转换成守恒变量增量
        if (nprec > nprec_non) then
            call calculate_precondition_matrix_point(nsp,nep,pvi,point_nonprec_i)
            call calculate_precondition_matrix_point(nsp,nep,pvj,point_nonprec_j)
            call calculate_precondition_matrix_point(nsp,nep,pvk,point_nonprec_k)
            
            !将原始变量增量转换成守恒变量增量ΔQ=P*ΔW
            call primitiveDeltaWToConservativeDeltaQ(point_nonprec_i,dqi,pvi,dUi)
            call primitiveDeltaWToConservativeDeltaQ(point_nonprec_j,dqj,pvj,dUj)
            if (nsw_kdir == nsw_dir_close) then
                dUk(:) = zero
            else
                call primitiveDeltaWToConservativeDeltaQ(point_nonprec_k,dqk,pvk,dUk)
            end if
        end if

#ifndef MemOPT1
        kx = sxyz(1)%r3d(i1,j,k)
        ky = sxyz(2)%r3d(i1,j,k)
        kz = sxyz(3)%r3d(i1,j,k)
#else
        kx = sxyz2%r4d(1,i1,j,k)
        ky = sxyz2%r4d(2,i1,j,k)
        kz = sxyz2%r4d(3,i1,j,k)
#endif
		
        call mxdq_std(1,npvs,pvi,kt,kx,ky,kz,1,neqn,dUi,dfi,srci,-one)
        !call mxdq_std(1,npvs,pvi,kt,kx,ky,kz,1,neqn,dqi,dfi,srci,-one)


#ifndef MemOPT1        
        ex = sxyz(4)%r3d(i,j1,k)
        ey = sxyz(5)%r3d(i,j1,k)
        ez = sxyz(6)%r3d(i,j1,k)
#else
        ex = sxyz2%r4d(4,i,j1,k)
        ey = sxyz2%r4d(5,i,j1,k)
        ez = sxyz2%r4d(6,i,j1,k)
#endif
		
        call mxdq_std(1,npvs,pvj,et,ex,ey,ez,1,neqn,dUj,dfj,srcj,-one)
        !call mxdq_std(1,npvs,pvj,et,ex,ey,ez,1,neqn,dqj,dfj,srcj,-one)

        if (nsw_kdir == nsw_dir_close) then
            dfk(:) = zero
        else


#ifndef MemOPT1        
          cx = sxyz(7)%r3d(i,j,k1)
          cy = sxyz(8)%r3d(i,j,k1)
          cz = sxyz(9)%r3d(i,j,k1)
#else
          cx = sxyz2%r4d(7,i,j,k1)
          cy = sxyz2%r4d(8,i,j,k1)
          cz = sxyz2%r4d(9,i,j,k1)
#endif  
	        call mxdq_std(1,npvs,pvk,ct,cx,cy,cz,1,neqn,dUk,dfk,srck,-one)	  
            !call mxdq_std(1,npvs,pvk,ct,cx,cy,cz,1,neqn,dqk,dfk,srck,-one)
        end if

        rhs0(:) = rhs0(:) - (dfi(:) + dfj(:) + dfk(:))


        if (nvis > nvis_euler) then
#ifndef MemOPT3 
            srvi = srv(1)%r3d(i1,j,k)
            srvj = srv(2)%r3d(i,j1,k)
            srvk = srv(3)%r3d(i,j,k1)
#else
            srvi = srv_r4d(1, i1,j,k)
            srvj = srv_r4d(2, i,j1,k)
            srvk = srv_r4d(3, i,j,k1)
			
#endif	
            do m=1,neqn
                rhs0(m) = rhs0(m) + &
                          half*(srvi*dUi(m) + srvj*dUj(m) + srvk*dUk(m))
            end do
            !do m=1,neqn
            !    rhs0(m) = rhs0(m) + &
            !              half*(srvi*dqi(m) + srvj*dqj(m) + srvk*dqk(m))
            !end do
        end if

        !todo odia需要替换成D矩阵的逆矩阵
        odia = one/dia
        !do m = 1,25
        !    diagInvMatrix(m) = diagInv(m)%r3d(i,j,k)
        !end do
        
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
        
        !计算D矩阵
        diagElement = dt(1)%r3d(i,j,k)
        
        call calculateDiagMatrixInv(prim,diagElement,sumSpectralRadius,vis,length,diagInvMatrix)
        
        if (nprec > nprec_non) then
            do l = 1,5
                tmp = 0.
                do m = 1,5
                    oned_index = 5*(l-1)+m
                    tmp = tmp + diagInvMatrix(oned_index)*rhs0(m)
                end do
                dq(l)%r3d(i,j,k) = tmp
            end do
        else
            do m=1,neqn
                dq(m)%r3d(i,j,k) = odia*rhs0(m)
            end do
        end if
        
    end do
    end do
!$OMP end do nowait

#ifdef OMP_IMP
    if (idt < min(numt-1, max(ed(2)-st(2)+1, st(2)-ed(2)+1))) then
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

#ifdef MemOPT2 
    !释放临时数组分配的存储空间

    if (nsw > 0) then
       deallocate( q0_r4d )                        
       deallocate( qm_r4d )         
    endif
    
    deallocate( qc_r4d  )
    deallocate( rhs_r4d )
    deallocate( pv_r4d  )    
#endif	

#ifdef MemOPT3 
    !释放临时数组分配的存储空间
   
    deallocate( src_r4d  )
    deallocate( srv_r4d )
 
#endif		
 
end subroutine prsgs_ust_sca_sweep

subroutine prsgs_ust_sca_sweep_sp(nb,swp,cpr,nsw)
!   CHEyg
!#define MemOPT1	! 优化sxyz访问   在makefile中定义

!#define MemOPT2	! 优化rhs,pv,qm,q0,qc访问

! #define MemOPT3	! 优化src和srv访问

#ifdef MemOPT1	
    use mod_fieldvars, only : nblkcoms, blkcomssp    ! CHEyg 本地块数和块指针
    use mod_opt_vars,  only : dim4_fld_array_t, sxyz1
#endif
	
	
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,one,half
    use mod_constants, only : nvis_euler,nsw_dir_close,nprec_non,nlhs_prsgs_ust_sca
    use mod_datatypes, only : fld_array_t
    use mod_variables, only : nvis,dtime,nsw_kdir,nprec,nlhs,fsw_kdir
    use mod_fieldvars, only : mb_topsp,mb_volsp,mb_sxyzsp
    !use mod_fieldvars, only : mb_nonprec,mb_prec_diag_inv
    use mod_fieldvars, only : mb_pv,npvs,mb_qc,neqn
    use mod_fieldvars, only : mb_dq,mb_rhs,mb_dt
    use mod_fieldvars, only : mb_src,mb_srv,mb_rdt
    use mod_fieldvars, only : mb_vsl,mb_vst
    use mod_fieldvars, only : mb_q0,mb_qm
    use mod_openmp
    implicit none
    integer(kind_int), intent(in) :: nb,swp
    real(kind_real),   intent(in) :: cpr(3)
    integer(kind_int), intent(in) :: nsw
    integer(kind_int)          :: i,j,k,i1,j1,k1,m,ierr,l,oned_index,nsp,nep
    integer(kind_int)          :: st(3),ed(3),dijk
    real(kind_real)            :: srci,srcj,srck,srvi,srvj,srvk
    real(kind_real)            :: kt,kx,ky,kz,et,ex,ey,ez,ct,cx,cy,cz
    real(kind_real)            :: dia,odia,odt,dq1,dq2,rhs0(1:neqn)
    real(kind_real)            :: pvi(1:npvs),pvj(1:npvs),pvk(1:npvs)
    real(kind_real)            :: dqi(1:neqn),dqj(1:neqn),dqk(1:neqn)
    real(kind_real)            :: dUi(1:neqn),dUj(1:neqn),dUk(1:neqn)
    real(kind_real)            :: dfi(1:neqn),dfj(1:neqn),dfk(1:neqn)
    real(kind_real)            :: point_nonprec_i(1:4),point_nonprec_j(1:4),point_nonprec_k(1:4)
    type(fld_array_t), pointer :: vol(:),sxyz(:),pv(:),nonprec(:),diagInv(:)
    type(fld_array_t), pointer :: dt(:),dq(:),rhs(:)
    type(fld_array_t), pointer :: src(:),srv(:),rdt(:),vsl(:),vst(:)
    type(fld_array_t), pointer :: qc(:),q0(:),qm(:)
    real(kind_real)            :: diagInvMatrix(1:25),prim(1:5),sn(1:3,1:3),eigenValue(1:3,1:5)
    real(kind_real)            :: tmp,rca,rcb,rcc,sumSpectralRadius,visl,vist,vis,length,diagElement
#ifdef OMP_IMP
    integer(kind_int)          :: idt, numt, nflag(0:omp_max_num_threads)
#endif

! CHEyg
#ifdef MemOPT1	
    integer(kind_int)               :: nc,nc1    !局部块编号
    type(dim4_fld_array_t), pointer :: sxyz2
#endif

#ifdef MemOPT2    
    real(kind_real), allocatable,dimension(:,:,:,:) :: qc_r4d, q0_r4d, rhs_r4d, qm_r4d ,pv_r4d
    integer(kind_int)                               :: nerror(5), i2,j2,k2

#endif

#ifdef MemOPT3    
    real(kind_real), allocatable,dimension(:,:,:,:) :: srv_r4d, src_r4d
#endif

    vol  => mb_volsp(nb)%fld
    sxyz => mb_sxyzsp(nb)%fld
    pv   => mb_pv(nb)%fld
    qc   => mb_qc(nb)%fld
    dt   => mb_dt(nb)%fld
    dq   => mb_dq(nb)%fld
    rhs  => mb_rhs(nb)%fld
    src  => mb_src(nb)%fld
    rdt  => mb_rdt(nb)%fld

    if (nvis > nvis_euler) then
        srv  => mb_srv(nb)%fld
        vsl => mb_vsl(nb)%fld
        vst => mb_vst(nb)%fld
        srv => mb_srv(nb)%fld
    end if
    
    !if (nprec > nprec_non) then
    !    nonprec => mb_nonprec(nb)%fld     !<守恒变量对原始变量的转换矩阵P
    !    diagInv => mb_prec_diag_inv(nb)%fld   !<对角矩阵的逆矩阵
    !end if

    if (nsw > 0) then
        q0 => mb_q0(nb)%fld !< Q(n-1)
        qm => mb_qm(nb)%fld !< Q(n)
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

! CHEyg
#ifdef MemOPT1	
    ! prsgs_ust_mat_sweep等被调用时给出的 nb 是全局的块编号，需要根据给出的全局 nb 编号算出在本进程上的块号
      do nc1 = 1,nblkcoms
         if(nb .EQ. blkcomssp(nc1)%nb) then
            nc = nc1
         end if
      enddo

      sxyz2 => sxyz1(nc)   ! 指向当前要用的块
      
#endif     
	  
#ifdef MemOPT2 
    !给临时数组分配存储空间
      if(dijk .EQ. 1) then
         i1 = st(1) 
         i2 = ed(1) 
         
         j1 = st(2)  
         j2 = ed(2)   
                 
         k1 = st(3)  
         k2 = ed(3)  
      else 
         i2 = st(1) 
         i1 = ed(1) 
         
         j2 = st(2)  
         j1 = ed(2)   
                 
         k2 = st(3)  
         k1 = ed(3)        
      endif

       if (nsw > 0) then
          allocate( q0_r4d(neqn,  i1:i2,  j1:j2, k1:k2), stat= nerror(1) )                        
          allocate( qm_r4d(neqn,  i1:i2,  j1:j2, k1:k2), stat= nerror(2) ) 
       
          do i = 1, 2
             if(nerror(i) .NE. 0) then
                write(*,*) "In prsgs_ust_mat_sweep. Allocating memory for X_r4d error. Stop! "
                stop
             endif
          enddo                 
       endif
    
       allocate( qc_r4d(  neqn,  i1:i2, j1:j2, k1:k2), stat= nerror(3) )
       allocate( rhs_r4d( neqn,  i1:i2, j1:j2, k1:k2), stat= nerror(4) )
       allocate( pv_r4d(  neqn,  (i1-1):(i2+1),  (j1-1):(j2+1), (k1-1):(k2+1) ), stat= nerror(5) )    
       ! pv的维数上下界均多一个
       
       do i = 3, 5
          if(nerror(i) .NE. 0) then
             write(*,*) "In prsgs_ust_mat_sweep. Allocating memory for X_r4d error. Stop! "
             stop
          endif
       enddo
       
    !将数据复制到临时数组中。逐个数组的赋值以避免cache冲突
         
       do k = k1,k2
       do j = j1,j2
       do i = i1,i2 
          do m = 1,neqn
             qc_r4d(m,i,j,k) = qc(m)%r3d(i,j,k)
          enddo
       enddo
       enddo
       enddo

   
    if (nsw > 0) then   
       do k = k1,k2
       do j = j1,j2
       do i = i1,i2 
          do m = 1,neqn
             qm_r4d(m,i,j,k) = qm(m)%r3d(i,j,k)
          enddo
       enddo
       enddo
       enddo   
                     
       do k = k1,k2
       do j = j1,j2
       do i = i1,i2 
          do m = 1,neqn
             q0_r4d(m,i,j,k) = q0(m)%r3d(i,j,k)
          enddo
       enddo
       enddo
       enddo

     endif

       do k = k1,k2
       do j = j1,j2
       do i = i1,i2 
          do m = 1,neqn
             rhs_r4d(m,i,j,k) = rhs(m)%r3d(i,j,k)
          enddo
       enddo
       enddo
       enddo
            
       do k = (k1-1),(k2+1)
       do j = (j1-1),(j2+1)
       do i = (i1-1),(i2+1) 
          do m = 1,neqn
             pv_r4d(m,i,j,k) = pv(m)%r3d(i,j,k)
          enddo
       enddo
       enddo
       enddo                           
        
#endif    

#ifdef MemOPT3    
    !给临时数组分配存储空间
      if(dijk .EQ. 1) then
         i1 = st(1) 
         i2 = ed(1) 
         
         j1 = st(2)  
         j2 = ed(2)   
                 
         k1 = st(3)  
         k2 = ed(3)  
      else 
         i2 = st(1) 
         i1 = ed(1) 
         
         j2 = st(2)  
         j1 = ed(2)   
                 
         k2 = st(3)  
         k1 = ed(3)        
      endif

       allocate( src_r4d(  neqn,  (i1-1):(i2+1),  (j1-1):(j2+1), (k1-1):(k2+1) ), stat= nerror(1) )
       allocate( srv_r4d(  neqn,  (i1-1):(i2+1),  (j1-1):(j2+1), (k1-1):(k2+1) ), stat= nerror(2) )

       do k = (k1-1),(k2+1)
       do j = (j1-1),(j2+1)
       do i = (i1-1),(i2+1) 
          do m = 1,3
             src_r4d(m,i,j,k) = src(m)%r3d(i,j,k)
          enddo
       enddo
       enddo
       enddo 

       do k = (k1-1),(k2+1)
       do j = (j1-1),(j2+1)
       do i = (i1-1),(i2+1) 
          do m = 1,3
             srv_r4d(m,i,j,k) = src(m)%r3d(i,j,k)
          enddo
       enddo
       enddo
       enddo 
	   
#endif
	
	

!$OMP parallel private(k,k1,j,j1,i,i1,odt,idt,pvi,pvj,pvk,dqi,dqj,dqk,rhs0,dfi,dfj,dfk,kx,ky,kz,&
!$OMP                  ex,ey,ez,cx,cy,cz,srvi,srvj,srvk,srci,srcj,srck,odia,dq1,dq2,dia)

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
        do m=1,neqn
            
#ifndef MemOPT2          
            rhs0(m) = rhs(m)%r3d(i,j,k)  
#else
            rhs0(m) = rhs_r4d(m, i,j,k)
            
#endif    			
        end do

        dia = rdt(1)%r3d(i,j,k)

        if (nsw > 0) then
            odt = vol(1)%r3d(i,j,k)/dtime
            do m=1,neqn

#ifndef MemOPT2                       
                dq1 = qm(m)%r3d(i,j,k) - q0(m)%r3d(i,j,k)
                dq2 = qc(m)%r3d(i,j,k) - qm(m)%r3d(i,j,k)
#else
                dq1 = qm_r4d(m, i,j,k) - q0_r4d(m, i,j,k)
                dq2 = qc_r4d(m, i,j,k) - qm_r4d(m, i,j,k)

#endif      

                rhs0(m) = rhs0(m) - odt*(cpr(1)*dq2 - cpr(2)*dq1)
            end do

            dia = dia + cpr(3)*odt
        end if

        i1 = i - 1
        j1 = j - 1
        k1 = k - 1

#ifndef MemOPT3    
        srci = src(1)%r3d(i1,j,k)
        srcj = src(2)%r3d(i,j1,k)
        srck = src(3)%r3d(i,j,k1)
#else 
        srci = src_r4d(1, i1,j,k)
        srcj = src_r4d(2, i,j1,k)
        srck = src_r4d(3, i,j,k1)
#endif		
		

        do m=1,npvs

#ifndef MemOPT2        
            pvi(m) = pv(m)%r3d(i1,j,k)
            pvj(m) = pv(m)%r3d(i,j1,k)
            pvk(m) = pv(m)%r3d(i,j,k1)
#else
            pvi(m) = pv_r4d(m, i1,j,k)
            pvj(m) = pv_r4d(m, i,j1,k)
            pvk(m) = pv_r4d(m, i,j,k1)   

#endif      

        end do

        do m=1,neqn
            dqi(m) = dq(m)%r3d(i1,j,k)
            dqj(m) = dq(m)%r3d(i,j1,k)
            dqk(m) = dq(m)%r3d(i,j,k1)
			dUi(m) = dqi(m)
            dUj(m) = dqj(m)
            dUk(m) = dqk(m)
        end do
 
        nsp = 1
        nep = 5
        sn(:,:) = 1.0
        vis = 1.0
        length = 1.0
        
        !todo 原始变量增量转换成守恒变量增量ΔQ=P*ΔW
        if (nprec > nprec_non) then
            call calculate_precondition_matrix_point(nsp,nep,pvi,point_nonprec_i)
            call calculate_precondition_matrix_point(nsp,nep,pvj,point_nonprec_j)
            call calculate_precondition_matrix_point(nsp,nep,pvk,point_nonprec_k)

            !将原始变量增量转换成守恒变量增量ΔQ=P*ΔW
            call primitiveDeltaWToConservativeDeltaQ(point_nonprec_i,dqi,pvi,dUi)
            call primitiveDeltaWToConservativeDeltaQ(point_nonprec_j,dqj,pvj,dUj)
            if (nsw_kdir == nsw_dir_close) then
                dUk(:) = zero
            else
                call primitiveDeltaWToConservativeDeltaQ(point_nonprec_k,dqk,pvk,dUk)
            end if
        end if

#ifndef MemOPT1
        kx = sxyz(1)%r3d(i1,j,k)
        ky = sxyz(2)%r3d(i1,j,k)
        kz = sxyz(3)%r3d(i1,j,k)
#else
        kx = sxyz2%r4d(1,i1,j,k)
        ky = sxyz2%r4d(2,i1,j,k)
        kz = sxyz2%r4d(3,i1,j,k)
#endif
		
        call mxdq_std(1,npvs,pvi,kt,kx,ky,kz,1,neqn,dUi,dfi,srci,one)
        !call mxdq_std(1,npvs,pvi,kt,kx,ky,kz,1,neqn,dqi,dfi,srci,one)


#ifndef MemOPT1        
        ex = sxyz(4)%r3d(i,j1,k)
        ey = sxyz(5)%r3d(i,j1,k)
        ez = sxyz(6)%r3d(i,j1,k)
#else
        ex = sxyz2%r4d(4,i,j1,k)
        ey = sxyz2%r4d(5,i,j1,k)
        ez = sxyz2%r4d(6,i,j1,k)
#endif
		
        call mxdq_std(1,npvs,pvj,et,ex,ey,ez,1,neqn,dUj,dfj,srcj,one)
        !call mxdq_std(1,npvs,pvj,et,ex,ey,ez,1,neqn,dqj,dfj,srcj,one)

        if (nsw_kdir == nsw_dir_close) then
            dfk(:) = zero
        else

#ifndef MemOPT1        
          cx = sxyz(7)%r3d(i,j,k1)
          cy = sxyz(8)%r3d(i,j,k1)
          cz = sxyz(9)%r3d(i,j,k1)
#else
          cx = sxyz2%r4d(7,i,j,k1)
          cy = sxyz2%r4d(8,i,j,k1)
          cz = sxyz2%r4d(9,i,j,k1)
#endif  
		  
            call mxdq_std(1,npvs,pvk,ct,cx,cy,cz,1,neqn,dUk,dfk,srck,one)
            !call mxdq_std(1,npvs,pvk,ct,cx,cy,cz,1,neqn,dqk,dfk,srck,one)
        end if

        rhs0(:) = rhs0(:) + (dfi(:) + dfj(:) + dfk(:))

        if (nvis > nvis_euler) then
		
#ifndef MemOPT3    
            srvi = srv(1)%r3d(i1,j,k)
            srvj = srv(2)%r3d(i,j1,k)
            srvk = srv(3)%r3d(i,j,k1)
#else
			srvi = srv_r4d(1, i1,j,k)
            srvj = srv_r4d(2, i,j1,k)
            srvk = srv_r4d(3, i,j,k1)
#endif
			do m=1,neqn
                rhs0(m) = rhs0(m) + &
                          half*(srvi*dUi(m) + srvj*dUj(m) + srvk*dUk(m))
            end do
            !do m=1,neqn
            !    rhs0(m) = rhs0(m) + &
            !              half*(srvi*dqi(m) + srvj*dqj(m) + srvk*dqk(m))
            !end do
            
        end if


        i1 = i + 1
        j1 = j + 1
        k1 = k + 1

#ifndef MemOPT3 
        srci = src(1)%r3d(i1,j,k)
        srcj = src(2)%r3d(i,j1,k)
        srck = src(3)%r3d(i,j,k1)
#else
        srci = src_r4d(1, i1,j,k)
        srcj = src_r4d(2, i,j1,k)
        srck = src_r4d(3, i,j,k1)
		
#endif 	

        do m=1,npvs

#ifndef MemOPT2
            pvi(m)  = pv(m)%r3d(i1,j,k)
            pvj(m)  = pv(m)%r3d(i,j1,k)
            pvk(m)  = pv(m)%r3d(i,j,k1)
#else
            pvi(m)  = pv_r4d(m, i1,j,k)
            pvj(m)  = pv_r4d(m, i,j1,k)
            pvk(m)  = pv_r4d(m, i,j,k1)
             
#endif
            prim(m) = pv(m)%r3d(i,j,k)
			 
        end do

        do m=1,neqn
            dqi(m) = dq(m)%r3d(i1,j,k)
            dqj(m) = dq(m)%r3d(i,j1,k)
            dqk(m) = dq(m)%r3d(i,j,k1)
            dUi(m) = dqi(m)
            dUj(m) = dqj(m)
            dUk(m) = dqk(m)
        end do

        

        !todo 原始变量增量转换成守恒变量增量
        if (nprec > nprec_non) then
            call calculate_precondition_matrix_point(nsp,nep,pvi,point_nonprec_i)
            call calculate_precondition_matrix_point(nsp,nep,pvj,point_nonprec_j)
            call calculate_precondition_matrix_point(nsp,nep,pvk,point_nonprec_k)
            
            !将原始变量增量转换成守恒变量增量ΔQ=P*ΔW
            call primitiveDeltaWToConservativeDeltaQ(point_nonprec_i,dqi,pvi,dUi)
            call primitiveDeltaWToConservativeDeltaQ(point_nonprec_j,dqj,pvj,dUj)
            if (nsw_kdir == nsw_dir_close) then
                dUk(:) = zero
            else
                call primitiveDeltaWToConservativeDeltaQ(point_nonprec_k,dqk,pvk,dUk)
            end if
        end if

#ifndef MemOPT1
        kx = sxyz(1)%r3d(i1,j,k)
        ky = sxyz(2)%r3d(i1,j,k)
        kz = sxyz(3)%r3d(i1,j,k)
#else
        kx = sxyz2%r4d(1,i1,j,k)
        ky = sxyz2%r4d(2,i1,j,k)
        kz = sxyz2%r4d(3,i1,j,k)
#endif
		
        call mxdq_std(1,npvs,pvi,kt,kx,ky,kz,1,neqn,dUi,dfi,srci,-one)
        !call mxdq_std(1,npvs,pvi,kt,kx,ky,kz,1,neqn,dqi,dfi,srci,-one)


#ifndef MemOPT1        
        ex = sxyz(4)%r3d(i,j1,k)
        ey = sxyz(5)%r3d(i,j1,k)
        ez = sxyz(6)%r3d(i,j1,k)
#else
        ex = sxyz2%r4d(4,i,j1,k)
        ey = sxyz2%r4d(5,i,j1,k)
        ez = sxyz2%r4d(6,i,j1,k)
#endif
		
        call mxdq_std(1,npvs,pvj,et,ex,ey,ez,1,neqn,dUj,dfj,srcj,-one)
        !call mxdq_std(1,npvs,pvj,et,ex,ey,ez,1,neqn,dqj,dfj,srcj,-one)

        if (nsw_kdir == nsw_dir_close) then
            dfk(:) = zero
        else


#ifndef MemOPT1        
          cx = sxyz(7)%r3d(i,j,k1)
          cy = sxyz(8)%r3d(i,j,k1)
          cz = sxyz(9)%r3d(i,j,k1)
#else
          cx = sxyz2%r4d(7,i,j,k1)
          cy = sxyz2%r4d(8,i,j,k1)
          cz = sxyz2%r4d(9,i,j,k1)
#endif  
	        call mxdq_std(1,npvs,pvk,ct,cx,cy,cz,1,neqn,dUk,dfk,srck,-one)	  
            !call mxdq_std(1,npvs,pvk,ct,cx,cy,cz,1,neqn,dqk,dfk,srck,-one)
        end if

        rhs0(:) = rhs0(:) - (dfi(:) + dfj(:) + dfk(:))


        if (nvis > nvis_euler) then
#ifndef MemOPT3 
            srvi = srv(1)%r3d(i1,j,k)
            srvj = srv(2)%r3d(i,j1,k)
            srvk = srv(3)%r3d(i,j,k1)
#else
            srvi = srv_r4d(1, i1,j,k)
            srvj = srv_r4d(2, i,j1,k)
            srvk = srv_r4d(3, i,j,k1)
			
#endif	
            do m=1,neqn
                rhs0(m) = rhs0(m) + &
                          half*(srvi*dUi(m) + srvj*dUj(m) + srvk*dUk(m))
            end do
            !do m=1,neqn
            !    rhs0(m) = rhs0(m) + &
            !              half*(srvi*dqi(m) + srvj*dqj(m) + srvk*dqk(m))
            !end do
        end if

        !todo odia需要替换成D矩阵的逆矩阵
        odia = one/dia
        !do m = 1,25
        !    diagInvMatrix(m) = diagInv(m)%r3d(i,j,k)
        !end do
        
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
        
        !计算D矩阵
        diagElement = dt(1)%r3d(i,j,k)
        
        call calculateDiagMatrixInv(prim,diagElement,sumSpectralRadius,vis,length,diagInvMatrix)
        
        if (nprec > nprec_non) then
            do l = 1,5
                tmp = 0.
                do m = 1,5
                    oned_index = 5*(l-1)+m
                    tmp = tmp + diagInvMatrix(oned_index)*rhs0(m)
                end do
                dq(l)%r3d(i,j,k) = tmp
            end do
        else
            do m=1,neqn
                dq(m)%r3d(i,j,k) = odia*rhs0(m)
            end do
        end if
        
    end do
    end do
!$OMP end do nowait

#ifdef OMP_IMP
    if (idt < min(numt-1, max(ed(2)-st(2)+1, st(2)-ed(2)+1))) then
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

#ifdef MemOPT2 
    !释放临时数组分配的存储空间

    if (nsw > 0) then
       deallocate( q0_r4d )                        
       deallocate( qm_r4d )         
    endif
    
    deallocate( qc_r4d  )
    deallocate( rhs_r4d )
    deallocate( pv_r4d  )    
#endif	

#ifdef MemOPT3 
    !释放临时数组分配的存储空间
   
    deallocate( src_r4d  )
    deallocate( srv_r4d )
 
#endif		
 
end subroutine prsgs_ust_sca_sweep_sp


