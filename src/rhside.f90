
subroutine rhside
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,nvis_euler,nsgl_buffer_pvs
    use mod_constants, only : nsgl_aver_vol,nbc_inter_buf_dqc
    use mod_variables, only : nvis,nghnode
    use mod_fieldvars, only : neqn,mb_rhs
    use mod_interface, only : exchange_bc_var,average_bc_var,pre_exchange_bc_var,post_exchange_bc_var
    use mod_interface, only : exchange_singulars,average_singulars
    use mod_interface, only : assign_mb_var_uniform,assign_bc_var_uniform
    implicit none
    real(kind_real) :: rhs0(1:neqn)
    integer :: m

    rhs0(:) = zero
    call assign_mb_var_uniform(mb_rhs,1,neqn,nghnode,rhs0)

    if (nvis > nvis_euler) then
        call pre_rhs_viscous
    end if

    call rhs_invscd

    if (nvis > nvis_euler) then
        call post_rhs_viscous
    end if

    call rhs_source
    !!call channel_source

!    call pre_exchange_bc_var(mb_rhs,1,neqn,1,nbc_inter_buf_dqc,nsgl_aver_vol)
!    call exchange_singulars(mb_rhs,1,neqn,nsgl_buffer_pvs,nsgl_aver_vol)
!    call post_exchange_bc_var(mb_rhs,1,neqn,1,nbc_inter_buf_dqc,nsgl_aver_vol)
!    call average_bc_var(mb_rhs,1,neqn,nbc_inter_buf_dqc,nsgl_aver_vol)
!    call average_singulars(mb_rhs,1,neqn,nsgl_buffer_pvs,nsgl_aver_vol)

!   0611
    call pre_exchange_bc_var(mb_rhs,1,neqn,1,nbc_inter_buf_dqc,nsgl_aver_vol)
    call exchange_singulars(mb_rhs,1,neqn,nsgl_buffer_pvs,nsgl_aver_vol)
    call post_exchange_bc_var(mb_rhs,1,neqn,1,nbc_inter_buf_dqc,nsgl_aver_vol)
    call average_bc_var(mb_rhs,1,neqn,nbc_inter_buf_dqc,nsgl_aver_vol)
    call average_singulars(mb_rhs,1,neqn,nsgl_buffer_pvs,nsgl_aver_vol)
    
    call assign_bc_var_uniform(mb_rhs,1,neqn,nghnode,rhs0)

end subroutine rhside

subroutine rhside_sp
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,nvis_euler,nsgl_buffer_pvs
    use mod_constants, only : nsgl_aver_vol,nbc_inter_buf_dqc
    use mod_variables, only : nvis,nghnode
    use mod_fieldvars, only : neqn,mb_rhs
    use mod_interface, only : pre_exchange_bc_var_sp,post_exchange_bc_var_sp
    use mod_interface, only : assign_mb_var_uniform_sp,assign_bc_var_uniform_sp
    implicit none
    real(kind_real) :: rhs0(1:neqn)
    integer :: m

    rhs0(:) = zero
    call assign_mb_var_uniform_sp(mb_rhs,1,neqn,nghnode,rhs0)

    if (nvis > nvis_euler) then
        call pre_rhs_viscous_sp
    end if

    call rhs_invscd_sp

    if (nvis > nvis_euler) then
        call post_rhs_viscous_sp
    end if

    call rhs_source_sp

!   0611
    !call pre_exchange_bc_var_sp(mb_rhs,1,neqn,1,nbc_inter_buf_dqc,nsgl_aver_vol)
    !call post_exchange_bc_var_sp(mb_rhs,1,neqn,1,nbc_inter_buf_dqc,nsgl_aver_vol)
    
    !call assign_bc_var_uniform_sp(mb_rhs,1,neqn,nghnode,rhs0)

end subroutine rhside_sp




