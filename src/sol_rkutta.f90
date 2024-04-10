
subroutine sol_rkutta_3step
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : ntimeadv_steady
    use mod_fieldvars, only : mb_qc,mb_q0,neqn
    use mod_interface, only : calc_mb_var_via_sub
    implicit none
    external                   :: v1_eq_v2
    real(kind_real), parameter :: crk(3)=(/0.d0,0.75d0,1.d0/3.d0/)
    integer(kind_int)          :: n

    call calc_mb_var_via_sub(mb_qc,1,neqn,v1_eq_v2,mb_q0,1,neqn,0)

    do n=1,3
        call rkutta3s_one_step(crk(n))
    end do

    call residual(ntimeadv_steady)

end subroutine sol_rkutta_3step

subroutine sol_rkutta_3step_sp
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : ntimeadv_steady
    use mod_fieldvars, only : mb_qc,mb_q0,neqn
    use mod_interface, only : calc_mb_var_via_sub_sp
    implicit none
    external                   :: v1_eq_v2
    real(kind_real), parameter :: crk(3)=(/0.d0,0.75d0,1.d0/3.d0/)
    integer(kind_int)          :: n

    call calc_mb_var_via_sub_sp(mb_qc,1,neqn,v1_eq_v2,mb_q0,1,neqn,0)

    do n=1,3
        call rkutta3s_one_step_sp(crk(n))
    end do

    call residual_sp(ntimeadv_steady)

end subroutine sol_rkutta_3step_sp

subroutine rkutta3s_one_step(crk)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,one,nsw_dir_close
    use mod_datatypes, only : top_block_t,fld_array_t
    use mod_variables, only : nsw_kdir
    use mod_fieldvars, only : nblkcoms,blkcoms
    use mod_fieldvars, only : mb_qc,mb_q0,neqn
    use mod_fieldvars, only : mb_dq,mb_dt,mb_rhs
    use mod_interface, only : assign_bc_var_nb
    use mod_interface, only : assign_var_via_com_nb
    implicit none
    real(kind_real), intent(in) :: crk
    integer(kind_int)          :: nc,nb,i,j,k,m,ierr
    integer(kind_int)          :: st(3),ed(3)
    real(kind_real)            :: dq0c,dtdq,dqn
    real(kind_real)            :: crk1,crk2,dq0(neqn)
    type(top_block_t), pointer :: top
    type(fld_array_t), pointer :: q0(:),qc(:),dt(:)
    type(fld_array_t), pointer :: dq(:),rhs(:)

    crk1 = crk
    crk2 = one - crk

    dq0(:) = zero

    call prepare_linear_system

    do nc=1,nblkcoms
        nb  =  blkcoms(nc)%nb
        top => blkcoms(nc)%top

        st(:) = top%ndst(:)
        ed(:) = top%nded(:)

        q0  => mb_q0(nb)%fld
        qc  => mb_qc(nb)%fld
        dt  => mb_dt(nb)%fld
        rhs => mb_rhs(nb)%fld
        dq  => mb_dq(nb)%fld

        do m=1,neqn
            do k=st(3),ed(3)
            do j=st(2),ed(2)
            do i=st(1),ed(1)
                dq0c = q0(m)%r3d(i,j,k) - qc(m)%r3d(i,j,k)
                dtdq = dt(1)%r3d(i,j,k)*rhs(m)%r3d(i,j,k)
                dqn  = crk1*dq0c + crk2*dtdq
                dq(m)%r3d(i,j,k) = dqn
            end do
            end do
            end do
        end do

        call assign_bc_var_nb(nb,mb_dq,1,neqn,0,dq0)

        call assign_var_via_com_nb(nb,mb_dq,1,neqn)
    end do

    call update
    
    !!call modify_boundary

end subroutine rkutta3s_one_step

subroutine rkutta3s_one_step_sp(crk)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,one,nsw_dir_close,nscmp_non
    use mod_datatypes, only : top_block_t,fld_array_t
    use mod_variables, only : nsw_kdir,nscmp
    use mod_fieldvars, only : nblkcoms,blkcomssp
    use mod_fieldvars, only : mb_qc,mb_q0,neqn
    use mod_fieldvars, only : mb_dq,mb_dt,mb_rhs
    use mod_interface, only : assign_var_via_com_nb_sp
    implicit none
    real(kind_real), intent(in) :: crk
    integer(kind_int)          :: nc,nb,i,j,k,m,ierr
    integer(kind_int)          :: st(3),ed(3)
    real(kind_real)            :: dq0c,dtdq,dqn
    real(kind_real)            :: crk1,crk2,dq0(neqn)
    type(top_block_t), pointer :: top
    type(fld_array_t), pointer :: q0(:),qc(:),dt(:)
    type(fld_array_t), pointer :: dq(:),rhs(:)

    crk1 = crk
    crk2 = one - crk

    dq0(:) = zero

    call prepare_linear_system_sp

    do nc=1,nblkcoms
        nb  =  blkcomssp(nc)%nb
        top => blkcomssp(nc)%top

        st(:) = top%ndst(:)
        ed(:) = top%nded(:)

        q0  => mb_q0(nb)%fld
        qc  => mb_qc(nb)%fld
        dt  => mb_dt(nb)%fld
        rhs => mb_rhs(nb)%fld
        dq  => mb_dq(nb)%fld

        do m=1,neqn
            do k=st(3),ed(3)
            do j=st(2),ed(2)
            do i=st(1),ed(1)
                dq0c = q0(m)%r3d(i,j,k) - qc(m)%r3d(i,j,k)
                dtdq = dt(1)%r3d(i,j,k)*rhs(m)%r3d(i,j,k)
                dqn  = crk1*dq0c + crk2*dtdq
                dq(m)%r3d(i,j,k) = dqn
            end do
            end do
            end do
        end do

        !call assign_bc_var_nb_sp(nb,mb_dq,1,neqn,0,dq0)

        call assign_var_via_com_nb_sp(nb,mb_dq,1,neqn)
    end do
    
    if (nscmp > nscmp_non) then
        call update_scmp_sp
    else
        call update_sp
    end if
    !!call modify_boundary
    !!call output_plt_sp

end subroutine rkutta3s_one_step_sp

