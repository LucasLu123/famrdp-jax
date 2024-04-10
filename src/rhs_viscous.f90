
subroutine pre_rhs_viscous
    use mod_kndconsts, only : kind_int
    use mod_fieldvars, only : mb_dpv
    use mod_variables, only : nghnode
    use mod_constants, only : nbc_inter_buf_dpv
    use mod_constants, only : nsgl_buffer_dpv,nsgl_aver_art
    use mod_interface, only : pre_exchange_bc_var,exchange_singulars
    implicit none

    call pv_gradient

    call pre_exchange_bc_var(mb_dpv,1,12,nghnode,nbc_inter_buf_dpv,nsgl_aver_art)
    call exchange_singulars(mb_dpv,1,12,nsgl_buffer_dpv,nsgl_aver_art)
end subroutine pre_rhs_viscous

subroutine pre_rhs_viscous_sp
    use mod_kndconsts, only : kind_int
    use mod_fieldvars, only : mb_dpv
    use mod_variables, only : nghnode
    use mod_constants, only : nbc_inter_buf_dpv
    use mod_constants, only : nsgl_buffer_dpv,nsgl_aver_art
    use mod_interface, only : pre_exchange_bc_var_sp
    implicit none

    call pv_gradient_sp

    call pre_exchange_bc_var_sp(mb_dpv,1,12,nghnode,nbc_inter_buf_dpv,nsgl_aver_art)
end subroutine pre_rhs_viscous_sp

subroutine post_rhs_viscous
    use mod_kndconsts, only : kind_int
    use mod_constants, only : nvis_ns_lam
    use mod_constants, only : nintplt_node2
    use mod_constants, only : nintplt_node4,nintplt_node6
    use mod_constants, only : nintplt_node6w,nintplt_node8
    use mod_constants, only : nintplt_scsl4,nintplt_scsl6
    use mod_constants, only : nintplt_node6e,nintplt_node8e
    use mod_constants, only : nintplt_scsl4e,nintplt_scsl6e
    use mod_variables, only : nvis,nd2int_vis
    use mod_fieldvars, only : nblkcoms,blkcoms
    use mod_fieldvars, only : mb_dpv
    use mod_variables, only : nghnode
    use mod_constants, only : nbc_inter_buf_dpv
    use mod_constants, only : nsgl_buffer_dpv,nsgl_aver_art
    use mod_interface, only : patched_ghost_points,post_exchange_bc_var
    use mod_interface, only : average_bc_var,average_singulars
    implicit none
    integer(kind_int) :: nc,nb
    external          :: ve_via_node2,ve_via_node4
    external          :: ve_via_node6,ve_via_node6w
    external          :: ve_via_node8
    external          :: ve_via_scsl4,ve_via_scsl6
    external          :: ve_via_node6e,ve_via_node8e
    external          :: ve_via_scsl4e,ve_via_scsl6e
    

    call post_exchange_bc_var(mb_dpv,1,12,nghnode,nbc_inter_buf_dpv,nsgl_aver_art)
    call patched_ghost_points(mb_dpv,1,12,nghnode)
    call average_bc_var(mb_dpv,1,12,nbc_inter_buf_dpv,nsgl_aver_art)
    call average_singulars(mb_dpv,1,12,nsgl_buffer_dpv,nsgl_aver_art)

    if (nvis > nvis_ns_lam) then
        call turbulent
    end if

    do nc=1,nblkcoms
        nb = blkcoms(nc)%nb
        select case(nd2int_vis)
        case(nintplt_node2)
            call viscous_scheme(nb,ve_via_node2)
        case(nintplt_node4)
            call viscous_scheme(nb,ve_via_node4)
        case(nintplt_node6,nintplt_node6w)
            call viscous_scheme(nb,ve_via_node6)
        case(nintplt_node6e)
            call viscous_scheme(nb,ve_via_node6e)
        case(nintplt_node8)
            call viscous_scheme(nb,ve_via_node8)
        case(nintplt_node8e)
            call viscous_scheme(nb,ve_via_node8e)
        case(nintplt_scsl4)
            call viscous_scheme(nb,ve_via_scsl4)
        case(nintplt_scsl4e)
            call viscous_scheme(nb,ve_via_scsl4e)
        case(nintplt_scsl6)
            call viscous_scheme(nb,ve_via_scsl6)
        case(nintplt_scsl6e)
            call viscous_scheme(nb,ve_via_scsl6e)
        case default

        end select
    end do

end subroutine post_rhs_viscous

subroutine post_rhs_viscous_sp
    use mod_kndconsts, only : kind_int
    use mod_constants, only : nvis_ns_lam
    use mod_constants, only : nintplt_node2
    use mod_constants, only : nintplt_node4,nintplt_node6
    use mod_constants, only : nintplt_node6w,nintplt_node8
    use mod_constants, only : nintplt_scsl4,nintplt_scsl6
    use mod_constants, only : nintplt_node6e,nintplt_node8e
    use mod_constants, only : nintplt_scsl4e,nintplt_scsl6e
    use mod_variables, only : nvis,nd2int_vis
    use mod_fieldvars, only : nblkcoms,blkcomssp
    use mod_fieldvars, only : mb_dpv
    use mod_variables, only : nghnode
    use mod_constants, only : nbc_inter_buf_dpv
    use mod_constants, only : nsgl_buffer_dpv,nsgl_aver_art
    use mod_interface, only : post_exchange_bc_var_sp
    implicit none
    integer(kind_int) :: nc,nb
    external          :: ve_via_node2cc,ve_via_node4cc
    external          :: ve_via_node6cc,ve_via_node6wcc
    external          :: ve_via_node8cc
    external          :: ve_via_scsl4cc,ve_via_scsl6cc
    external          :: ve_via_node6ecc,ve_via_node8ecc
    external          :: ve_via_scsl4ecc,ve_via_scsl6ecc
    

    call post_exchange_bc_var_sp(mb_dpv,1,12,nghnode,nbc_inter_buf_dpv,nsgl_aver_art)

    if (nvis > nvis_ns_lam) then
        !! to do call turbulent
        call msg_seq_and_master("Turbulence models are not ready for dual-mesh grid")
    end if

    do nc=1,nblkcoms
        nb = blkcomssp(nc)%nb
        select case(nd2int_vis)
        case(nintplt_node2)
            call viscous_scheme_sp(nb,ve_via_node2cc)
        case(nintplt_node4)
            call viscous_scheme_sp(nb,ve_via_node4cc)
        case(nintplt_node6,nintplt_node6w)
            call viscous_scheme_sp(nb,ve_via_node6cc)
        case(nintplt_node6e)
            call viscous_scheme_sp(nb,ve_via_node6ecc)
        case(nintplt_node8)
            call viscous_scheme_sp(nb,ve_via_node8cc)
        case(nintplt_node8e)
            call viscous_scheme_sp(nb,ve_via_node8ecc)
        case(nintplt_scsl4)
            call viscous_scheme_sp(nb,ve_via_scsl4cc)
        case(nintplt_scsl4e)
            call viscous_scheme_sp(nb,ve_via_scsl4ecc)
        case(nintplt_scsl6)
            call viscous_scheme_sp(nb,ve_via_scsl6cc)
        case(nintplt_scsl6e)
            call viscous_scheme_sp(nb,ve_via_scsl6ecc)
        case default

        end select
    end do

end subroutine post_rhs_viscous_sp

subroutine pv_gradient
    use mod_kndconsts, only : kind_int
    use mod_constants, only : nintplt_node2,nintplt_node4
    use mod_constants, only : nintplt_node6,nintplt_node6w,nintplt_node8
    use mod_constants, only : nintplt_scsl4,nintplt_scsl6
    use mod_constants, only : nintplt_node6e,nintplt_node8e
    use mod_constants, only : nintplt_scsl4e,nintplt_scsl6e
    use mod_variables, only : nd3int_vis
    implicit none
    external :: ve_via_node2,ve_via_node4
    external :: ve_via_node6,ve_via_node6w,ve_via_node8
    external :: ve_via_scsl4,ve_via_scsl6
    external :: ve_via_node6e,ve_via_node8e
    external :: ve_via_scsl4e,ve_via_scsl6e
    
    select case(nd3int_vis)
    case(nintplt_node2)
        call grad_pv_dnvis(ve_via_node2)
    case(nintplt_node4)
        call grad_pv_dnvis(ve_via_node4)
    case(nintplt_node6)
        call grad_pv_dnvis(ve_via_node6)
    case(nintplt_node6w)
        call grad_pv_dnvis(ve_via_node6w)
    case(nintplt_node6e)
        call grad_pv_dnvis(ve_via_node6e)
    case(nintplt_node8)
        call grad_pv_dnvis(ve_via_node8)
    case(nintplt_node8e)
        call grad_pv_dnvis(ve_via_node8e)
    case(nintplt_scsl4)
        call grad_pv_dnvis(ve_via_scsl4)
    case(nintplt_scsl4e)
        call grad_pv_dnvis(ve_via_scsl4e)
    case(nintplt_scsl6)
        call grad_pv_dnvis(ve_via_scsl6)
    case(nintplt_scsl6e)
        call grad_pv_dnvis(ve_via_scsl6e)                
    case default

    end select

end subroutine pv_gradient

subroutine pv_gradient_sp
    use mod_kndconsts, only : kind_int
    use mod_constants, only : nintplt_node2,nintplt_node4
    use mod_constants, only : nintplt_node6,nintplt_node6w,nintplt_node8
    use mod_constants, only : nintplt_scsl4,nintplt_scsl6
    use mod_constants, only : nintplt_node6e,nintplt_node8e
    use mod_constants, only : nintplt_scsl4e,nintplt_scsl6e
    use mod_variables, only : nd3int_vis
    implicit none
    external :: ve_via_node2cc,ve_via_node4cc
    external :: ve_via_node6cc,ve_via_node6wcc,ve_via_node8cc
    external :: ve_via_scsl4cc,ve_via_scsl6cc
    external :: ve_via_node6ecc,ve_via_node8ecc
    external :: ve_via_scsl4ecc,ve_via_scsl6ecc
    
    select case(nd3int_vis)
    case(nintplt_node2)
        call grad_pv_dnvis_sp(ve_via_node2cc)
    case(nintplt_node4)
        call grad_pv_dnvis_sp(ve_via_node4cc)
    case(nintplt_node6)
        call grad_pv_dnvis_sp(ve_via_node6cc)
    case(nintplt_node6w)
        call grad_pv_dnvis_sp(ve_via_node6wcc)
    case(nintplt_node6e)
        call grad_pv_dnvis_sp(ve_via_node6ecc)
    case(nintplt_node8)
        call grad_pv_dnvis_sp(ve_via_node8cc)
    case(nintplt_node8e)
        call grad_pv_dnvis_sp(ve_via_node8ecc)
    case(nintplt_scsl4)
        call grad_pv_dnvis_sp(ve_via_scsl4cc)
    case(nintplt_scsl4e)
        call grad_pv_dnvis_sp(ve_via_scsl4ecc)
    case(nintplt_scsl6)
        call grad_pv_dnvis_sp(ve_via_scsl6cc)
    case(nintplt_scsl6e)
        call grad_pv_dnvis_sp(ve_via_scsl6ecc)                
    case default

    end select

end subroutine pv_gradient_sp


subroutine grad_pv_dnvis(sub_intvis)
    use mod_constants, only : nderive_edge6,nderive_edge4,nderive_edge2
    use mod_constants, only : nderive_ehen6,nderive_ehen4,nderive_ehen8,nderive_ehcs6
    use mod_constants, only : nderive_scsl4,nderive_scsl6
    use mod_constants, only : nderive_ehen6e,nderive_ehen8e,nderive_ehcs6e
    use mod_constants, only : nderive_scsl4e,nderive_scsl6e
    use mod_variables, only : nd3der_vis,nghnode,nghedge
    implicit none
    external :: sub_intvis
    external :: dn_via_edge6,dn_via_edge4,dn_via_edge2
    external :: dn_via_ehen6,dn_via_ehen4,dn_via_ehen8,dn_via_ehcs6
    external :: dn_via_scsl4,dn_via_scsl6
    external :: dn_via_ehen6e,dn_via_ehen8e,dn_via_ehcs6e
    external :: dn_via_scsl4e,dn_via_scsl6e

    select case(nd3der_vis)
    case(nderive_edge2)
        call calc_grad_pv(nghnode,nghedge,sub_intvis,dn_via_edge2)
    case(nderive_edge4)
        call calc_grad_pv(nghnode,nghedge,sub_intvis,dn_via_edge4)
    case(nderive_edge6)
        call calc_grad_pv(nghnode,nghedge,sub_intvis,dn_via_edge6)
    case(nderive_ehen4)
        call calc_grad_pv(nghnode,nghedge,sub_intvis,dn_via_ehen4)
    case(nderive_ehen6)
        call calc_grad_pv(nghnode,nghedge,sub_intvis,dn_via_ehen6)
    case(nderive_ehen6e)
        call calc_grad_pv(nghnode,nghedge,sub_intvis,dn_via_ehen6e)
    case(nderive_ehcs6)
        call calc_grad_pv(nghnode,nghedge,sub_intvis,dn_via_ehcs6)
    case(nderive_ehcs6e)
        call calc_grad_pv(nghnode,nghedge,sub_intvis,dn_via_ehcs6e)
    case(nderive_ehen8)
        call calc_grad_pv(nghnode,nghedge,sub_intvis,dn_via_ehen8)
    case(nderive_ehen8e)
        call calc_grad_pv(nghnode,nghedge,sub_intvis,dn_via_ehen8e)
    case(nderive_scsl4)
        call calc_grad_pv(nghnode,nghedge,sub_intvis,dn_via_scsl4)
    case(nderive_scsl4e)
        call calc_grad_pv(nghnode,nghedge,sub_intvis,dn_via_scsl4e)
    case(nderive_scsl6)
        call calc_grad_pv(nghnode,nghedge,sub_intvis,dn_via_scsl6)
    case(nderive_scsl6e)
        call calc_grad_pv(nghnode,nghedge,sub_intvis,dn_via_scsl6e)                
    case default

    end select

end subroutine grad_pv_dnvis

subroutine grad_pv_dnvis_sp(sub_intvis)
    use mod_constants, only : nderive_edge6,nderive_edge4,nderive_edge2
    use mod_constants, only : nderive_ehen6,nderive_ehen4,nderive_ehen8,nderive_ehcs6
    use mod_constants, only : nderive_scsl4,nderive_scsl6
    use mod_constants, only : nderive_ehen6e,nderive_ehen8e,nderive_ehcs6e
    use mod_constants, only : nderive_scsl4e,nderive_scsl6e
    use mod_variables, only : nd3der_vis,nghnode,nghedge
    implicit none
    external :: sub_intvis
    external :: dn_via_edge6cc,dn_via_edge4cc,dn_via_edge2cc
    external :: dn_via_ehen6cc,dn_via_ehen4cc,dn_via_ehen8cc,dn_via_ehcs6cc
    external :: dn_via_scsl4cc,dn_via_scsl6cc
    external :: dn_via_ehen6ecc,dn_via_ehen8ecc,dn_via_ehcs6ecc
    external :: dn_via_scsl4ecc,dn_via_scsl6ecc

    select case(nd3der_vis)
    case(nderive_edge2)
        call calc_grad_pv_sp(nghnode,nghedge,sub_intvis,dn_via_edge2cc)
    case(nderive_edge4)
        call calc_grad_pv_sp(nghnode,nghedge,sub_intvis,dn_via_edge4cc)
    case(nderive_edge6)
        call calc_grad_pv_sp(nghnode,nghedge,sub_intvis,dn_via_edge6cc)
    case(nderive_ehen4)
        call calc_grad_pv_sp(nghnode,nghedge,sub_intvis,dn_via_ehen4cc)
    case(nderive_ehen6)
        call calc_grad_pv_sp(nghnode,nghedge,sub_intvis,dn_via_ehen6cc)
    case(nderive_ehen6e)
        call calc_grad_pv_sp(nghnode,nghedge,sub_intvis,dn_via_ehen6ecc)
    case(nderive_ehcs6)
        call calc_grad_pv_sp(nghnode,nghedge,sub_intvis,dn_via_ehcs6cc)
    case(nderive_ehcs6e)
        call calc_grad_pv_sp(nghnode,nghedge,sub_intvis,dn_via_ehcs6ecc)
    case(nderive_ehen8)
        call calc_grad_pv_sp(nghnode,nghedge,sub_intvis,dn_via_ehen8cc)
    case(nderive_ehen8e)
        call calc_grad_pv_sp(nghnode,nghedge,sub_intvis,dn_via_ehen8ecc)
    case(nderive_scsl4)
        call calc_grad_pv_sp(nghnode,nghedge,sub_intvis,dn_via_scsl4cc)
    case(nderive_scsl4e)
        call calc_grad_pv_sp(nghnode,nghedge,sub_intvis,dn_via_scsl4ecc)
    case(nderive_scsl6)
        call calc_grad_pv_sp(nghnode,nghedge,sub_intvis,dn_via_scsl6cc)
    case(nderive_scsl6e)
        call calc_grad_pv_sp(nghnode,nghedge,sub_intvis,dn_via_scsl6ecc)                
    case default

    end select

end subroutine grad_pv_dnvis_sp


subroutine calc_grad_pv(ngn,nge,sub_intvis,sub_dnvis)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : one,bc_cut1to1
    use mod_constants, only : nfsf_vis_d3int,nfsf_vis_d3der
    use mod_constants, only : nbc_inter_buf_dpv
    use mod_constants, only : nsgl_buffer_dpv,nsgl_aver_art
    use mod_datatypes, only : fld_array_t,var_block_t,top_block_t
    use mod_fieldvars, only : nblkcoms,blkcoms
    use mod_fieldvars, only : mb_pv,mb_t,mb_sxyz,mb_vol,mb_dpv
    use mod_interface, only : mb_var_pointer_create,mb_var_pointer_assign
    use mod_interface, only : mb_var_pointer_delete,calc_mb_dn_via_node3
    use mod_openmp
    implicit none
    integer(kind_int), intent(in) :: ngn,nge
    external                      :: sub_intvis,sub_dnvis
    integer(kind_int)          :: nc,nb,i,j,k,m,m1,m2,m3
    integer(kind_int)          :: st(3),ed(3)
    real(kind_real)            :: kx,ky,kz,ex,ey,ez,cx,cy,cz
    real(kind_real)            :: ovol,der(12)
    type(top_block_t), pointer :: top
    type(fld_array_t), pointer :: sxyz(:),vol(:),dpv(:)
    type(var_block_t), pointer :: mb_uvwt(:)

    call mb_var_pointer_create(mb_uvwt,1,4)
    call mb_var_pointer_assign(mb_uvwt,1,3,mb_pv,2,4)
    call mb_var_pointer_assign(mb_uvwt,4,4,mb_t ,1,1)
    do m=1,4
        m1 = 3*m - 2
        m2 = m1 + 1
        m3 = m2 + 1

        call calc_mb_dn_via_node3(mb_uvwt,m,ngn,sub_intvis,nge,sub_dnvis,nfsf_vis_d3int,nfsf_vis_d3der,mb_dpv,m1,1)
        call calc_mb_dn_via_node3(mb_uvwt,m,ngn,sub_intvis,nge,sub_dnvis,nfsf_vis_d3int,nfsf_vis_d3der,mb_dpv,m2,2)
        call calc_mb_dn_via_node3(mb_uvwt,m,ngn,sub_intvis,nge,sub_dnvis,nfsf_vis_d3int,nfsf_vis_d3der,mb_dpv,m3,3)
    end do
    call mb_var_pointer_delete(mb_uvwt)

    do nc=1,nblkcoms
        nb  =  blkcoms(nc)%nb
        top => blkcoms(nc)%top

        sxyz => mb_sxyz(nb)%fld
        vol  => mb_vol(nb)%fld
        dpv  => mb_dpv(nb)%fld

        st(:) = 1
        ed(:) = top%nijk(:)
!$OMP parallel do private(k,j,i,kx,ky,kz,ex,ey,ez,cx,cy,cz,ovol,m1,m2,m3,der)
        do k=st(3),ed(3)
        do j=st(2),ed(2)
        do i=st(1),ed(1)
            kx = sxyz(1)%r3d(i,j,k)
            ky = sxyz(2)%r3d(i,j,k)
            kz = sxyz(3)%r3d(i,j,k)

            ex = sxyz(4)%r3d(i,j,k)
            ey = sxyz(5)%r3d(i,j,k)
            ez = sxyz(6)%r3d(i,j,k)

            cx = sxyz(7)%r3d(i,j,k)
            cy = sxyz(8)%r3d(i,j,k)
            cz = sxyz(9)%r3d(i,j,k)

            ovol = one/vol(1)%r3d(i,j,k)

            do m=1,12
               der(m) = dpv(m)%r3d(i,j,k)
            end do

            do m=1,4
                m1 = 3*m - 2
                m2 = m1 + 1
                m3 = m2 + 1

                dpv(m1)%r3d(i,j,k) = ( kx*der(m1) + ex*der(m2) + cx*der(m3) ) * ovol
                dpv(m2)%r3d(i,j,k) = ( ky*der(m1) + ey*der(m2) + cy*der(m3) ) * ovol
                dpv(m3)%r3d(i,j,k) = ( kz*der(m1) + ez*der(m2) + cz*der(m3) ) * ovol
            end do
        end do
        end do
        end do
!$OMP end parallel do
    end do


end subroutine calc_grad_pv

subroutine calc_grad_pv_sp(ngn,nge,sub_intvis,sub_dnvis)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : one,bc_cut1to1
    use mod_constants, only : nfsf_vis_d3int,nfsf_vis_d3der
    use mod_constants, only : nbc_inter_buf_dpv
    use mod_constants, only : nsgl_buffer_dpv,nsgl_aver_art
    use mod_datatypes, only : fld_array_t,var_block_t,top_block_t
    use mod_fieldvars, only : nblkcoms,blkcomssp
    use mod_fieldvars, only : mb_pv,mb_t,mb_sxyzsp,mb_volsp,mb_dpv
    use mod_interface, only : mb_var_pointer_create,mb_var_pointer_assign
    use mod_interface, only : mb_var_pointer_delete,calc_mb_duvwt_via_node3_sp
    use mod_openmp
    implicit none
    integer(kind_int), intent(in) :: ngn,nge
    external                      :: sub_intvis,sub_dnvis
    integer(kind_int)          :: nc,nb,i,j,k,m,m1,m2,m3
    integer(kind_int)          :: st(3),ed(3)
    real(kind_real)            :: kx,ky,kz,ex,ey,ez,cx,cy,cz
    real(kind_real)            :: ovol,der(12)
    type(top_block_t), pointer :: top
    type(fld_array_t), pointer :: sxyz(:),vol(:),dpv(:)
    type(var_block_t), pointer :: mb_uvwt(:)

    call mb_var_pointer_create(mb_uvwt,1,4)
    call mb_var_pointer_assign(mb_uvwt,1,3,mb_pv,2,4)
    call mb_var_pointer_assign(mb_uvwt,4,4,mb_t ,1,1)
    do m=1,4
        m1 = 3*m - 2
        m2 = m1 + 1
        m3 = m2 + 1

        call calc_mb_duvwt_via_node3_sp(mb_uvwt,m,ngn,sub_intvis,nge,sub_dnvis,nfsf_vis_d3int,nfsf_vis_d3der,mb_dpv,m1,1)
        call calc_mb_duvwt_via_node3_sp(mb_uvwt,m,ngn,sub_intvis,nge,sub_dnvis,nfsf_vis_d3int,nfsf_vis_d3der,mb_dpv,m2,2)
        call calc_mb_duvwt_via_node3_sp(mb_uvwt,m,ngn,sub_intvis,nge,sub_dnvis,nfsf_vis_d3int,nfsf_vis_d3der,mb_dpv,m3,3)
    end do
    call mb_var_pointer_delete(mb_uvwt)

    do nc=1,nblkcoms
        nb  =  blkcomssp(nc)%nb
        top => blkcomssp(nc)%top

        sxyz => mb_sxyzsp(nb)%fld
        vol  => mb_volsp(nb)%fld
        dpv  => mb_dpv(nb)%fld

        st(:) = 1
        ed(:) = top%nijk(:)
!$OMP parallel do private(k,j,i,kx,ky,kz,ex,ey,ez,cx,cy,cz,ovol,m1,m2,m3,der)
        do k=st(3),ed(3)
        do j=st(2),ed(2)
        do i=st(1),ed(1)
            kx = sxyz(1)%r3d(i,j,k)
            ky = sxyz(2)%r3d(i,j,k)
            kz = sxyz(3)%r3d(i,j,k)

            ex = sxyz(4)%r3d(i,j,k)
            ey = sxyz(5)%r3d(i,j,k)
            ez = sxyz(6)%r3d(i,j,k)

            cx = sxyz(7)%r3d(i,j,k)
            cy = sxyz(8)%r3d(i,j,k)
            cz = sxyz(9)%r3d(i,j,k)

            ovol = one/vol(1)%r3d(i,j,k)

            do m=1,12
               der(m) = dpv(m)%r3d(i,j,k)
            end do

            do m=1,4
                m1 = 3*m - 2
                m2 = m1 + 1
                m3 = m2 + 1

                dpv(m1)%r3d(i,j,k) = ( kx*der(m1) + ex*der(m2) + cx*der(m3) ) * ovol
                dpv(m2)%r3d(i,j,k) = ( ky*der(m1) + ey*der(m2) + cy*der(m3) ) * ovol
                dpv(m3)%r3d(i,j,k) = ( kz*der(m1) + ez*der(m2) + cz*der(m3) ) * ovol
            end do
        end do
        end do
        end do
!$OMP end parallel do
    end do


end subroutine calc_grad_pv_sp

subroutine viscous_scheme(nb,sub_intplt)
    use mod_kndconsts, only : kind_int
    use mod_constants, only : nderive_edge6,nderive_edge4,nderive_edge2
    use mod_constants, only : nderive_ehen6,nderive_ehen4,nderive_ehen8,nderive_ehcs6
    use mod_constants, only : nderive_scsl4,nderive_scsl6
    use mod_constants, only : nderive_ehen6e,nderive_ehen8e,nderive_ehcs6e
    use mod_constants, only : nderive_scsl4e,nderive_scsl6e
    use mod_variables, only : nd2der_vis,nghnode,nghedge
    implicit none
    integer(kind_int), intent(in) :: nb
    external                      :: sub_intplt
    external :: dn_via_edge6,dn_via_edge4,dn_via_edge2
    external :: dn_via_ehen6,dn_via_ehen4,dn_via_ehen8,dn_via_ehcs6
    external :: dn_via_scsl4,dn_via_scsl6
    external :: dn_via_ehen6e,dn_via_ehen8e,dn_via_ehcs6e
    external :: dn_via_scsl4e,dn_via_scsl6e

    select case(nd2der_vis)
    case(nderive_edge2)
       call calc_viscous(nb,nghnode,nghedge,sub_intplt,dn_via_edge2)
    case(nderive_edge4)
       call calc_viscous(nb,nghnode,nghedge,sub_intplt,dn_via_edge4)
    case(nderive_edge6)
       call calc_viscous(nb,nghnode,nghedge,sub_intplt,dn_via_edge6)
    case(nderive_ehen4)
       call calc_viscous(nb,nghnode,nghedge,sub_intplt,dn_via_ehen4)
    case(nderive_ehen6)
       call calc_viscous(nb,nghnode,nghedge,sub_intplt,dn_via_ehen6)
    case(nderive_ehen6e)
       call calc_viscous(nb,nghnode,nghedge,sub_intplt,dn_via_ehen6e)
    case(nderive_ehcs6)
       call calc_viscous(nb,nghnode,nghedge,sub_intplt,dn_via_ehcs6)
    case(nderive_ehcs6e)
       call calc_viscous(nb,nghnode,nghedge,sub_intplt,dn_via_ehcs6e)
    case(nderive_ehen8)
       call calc_viscous(nb,nghnode,nghedge,sub_intplt,dn_via_ehen8)
    case(nderive_ehen8e)
       call calc_viscous(nb,nghnode,nghedge,sub_intplt,dn_via_ehen8e)
    case(nderive_scsl4)
       call calc_viscous(nb,nghnode,nghedge,sub_intplt,dn_via_scsl4)
    case(nderive_scsl4e)
       call calc_viscous(nb,nghnode,nghedge,sub_intplt,dn_via_scsl4e)
    case(nderive_scsl6)
       call calc_viscous(nb,nghnode,nghedge,sub_intplt,dn_via_scsl6)
    case(nderive_scsl6e)
       call calc_viscous(nb,nghnode,nghedge,sub_intplt,dn_via_scsl6e)
    case default

    end select

end subroutine viscous_scheme

subroutine viscous_scheme_sp(nb,sub_intplt)
    use mod_kndconsts, only : kind_int
    use mod_constants, only : nderive_edge6,nderive_edge4,nderive_edge2
    use mod_constants, only : nderive_ehen6,nderive_ehen4,nderive_ehen8,nderive_ehcs6
    use mod_constants, only : nderive_scsl4,nderive_scsl6
    use mod_constants, only : nderive_ehen6e,nderive_ehen8e,nderive_ehcs6e
    use mod_constants, only : nderive_scsl4e,nderive_scsl6e
    use mod_variables, only : nd2der_vis,nghnode,nghedge
    implicit none
    integer(kind_int), intent(in) :: nb
    external                      :: sub_intplt
    external :: dn_via_edge6cc,dn_via_edge4cc,dn_via_edge2cc
    external :: dn_via_ehen6cc,dn_via_ehen4cc,dn_via_ehen8cc,dn_via_ehcs6cc
    external :: dn_via_scsl4cc,dn_via_scsl6cc
    external :: dn_via_ehen6ecc,dn_via_ehen8ecc,dn_via_ehcs6ecc
    external :: dn_via_scsl4ecc,dn_via_scsl6ecc

    select case(nd2der_vis)
    case(nderive_edge2)
       call calc_viscous_sp(nb,nghnode,nghedge,sub_intplt,dn_via_edge2cc)
    case(nderive_edge4)
       call calc_viscous_sp(nb,nghnode,nghedge,sub_intplt,dn_via_edge4cc)
    case(nderive_edge6)
       call calc_viscous_sp(nb,nghnode,nghedge,sub_intplt,dn_via_edge6cc)
    case(nderive_ehen4)
       call calc_viscous_sp(nb,nghnode,nghedge,sub_intplt,dn_via_ehen4cc)
    case(nderive_ehen6)
       call calc_viscous_sp(nb,nghnode,nghedge,sub_intplt,dn_via_ehen6cc)
    case(nderive_ehen6e)
       call calc_viscous_sp(nb,nghnode,nghedge,sub_intplt,dn_via_ehen6ecc)
    case(nderive_ehcs6)
       call calc_viscous_sp(nb,nghnode,nghedge,sub_intplt,dn_via_ehcs6cc)
    case(nderive_ehcs6e)
       call calc_viscous_sp(nb,nghnode,nghedge,sub_intplt,dn_via_ehcs6ecc)
    case(nderive_ehen8)
       call calc_viscous_sp(nb,nghnode,nghedge,sub_intplt,dn_via_ehen8cc)
    case(nderive_ehen8e)
       call calc_viscous_sp(nb,nghnode,nghedge,sub_intplt,dn_via_ehen8ecc)
    case(nderive_scsl4)
       call calc_viscous_sp(nb,nghnode,nghedge,sub_intplt,dn_via_scsl4cc)
    case(nderive_scsl4e)
       call calc_viscous_sp(nb,nghnode,nghedge,sub_intplt,dn_via_scsl4ecc)
    case(nderive_scsl6)
       call calc_viscous_sp(nb,nghnode,nghedge,sub_intplt,dn_via_scsl6cc)
    case(nderive_scsl6e)
       call calc_viscous_sp(nb,nghnode,nghedge,sub_intplt,dn_via_scsl6ecc)
    case default

    end select

end subroutine viscous_scheme_sp

subroutine calc_viscous(nb,ngn,nge,sub_intplt,sub_scheme)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,one
!    use mod_constants, only : nbc_bound_scheme,nsw_dir_close
    use mod_constants, only : nbc_inter_scheme,nsw_dir_close
    use mod_constants, only : nbc_intbc_scheme
    
    use mod_constants, only : nfsf_vis_d2int,nfsf_vis_d2der
    use mod_datatypes, only : fld_array_t
    use mod_variables, only : nsw_kdir,gamma,prlam,prtur,refbeta,reue
    use mod_fieldvars, only : mb_top,npvs,mb_pv,neqn,mb_t,mb_sxyz,mb_vol
    use mod_fieldvars, only : mb_vsl,mb_vst,mb_rhs,mb_fsf,mb_dpv
    use mod_openmp
    implicit none
    integer(kind_int), intent(in) :: nb,ngn,nge
    external                      :: sub_intplt,sub_scheme
    integer(kind_int)          :: i,j,k,m,ierr
    integer(kind_int)          :: ni,nj,nk,stn,edn,ste,ede
    integer(kind_int)          :: stn0,edn0,ste0,ede0,nfs,nfe
    integer(kind_int)          :: nkst,nked,njed
    real(kind_real)            :: re,cp,gama,ae,vx,vy,vz,nt
    real(kind_real)            :: kx,ky,kz,ex,ey,ez,cx,cy,cz
    real(kind_real)            :: vis,kcp,cp_prl,cp_prt
    real(kind_real)            :: f(1:neqn),der(12)
    real(kind_real), pointer   :: vn(:,:),ve(:,:)
    real(kind_real), pointer   :: fn(:,:),fc(:,:),dn(:,:)
    type(fld_array_t), pointer :: dpv(:),pv(:),rhs(:)
    type(fld_array_t), pointer :: sxyz(:),vsl(:),vst(:)
    type(fld_array_t), pointer :: fsfs(:),fsfe(:)

    gama = gamma
    ae = gama - one

    cp = gama*refbeta/ae
    cp_prl = cp/prlam
    cp_prt = cp/prtur

    re = one/reue

    ni = mb_top(nb)%nijk(1)
    nj = mb_top(nb)%nijk(2)
    nk = mb_top(nb)%nijk(3)

    nkst = mb_top(nb)%ndst(3)
    nked = mb_top(nb)%nded(3)
    if (nsw_kdir == nsw_dir_close) then
       njed = -1
    else
       njed = nj
    end if

    dpv  => mb_dpv(nb)%fld
    pv   => mb_pv(nb)%fld
    rhs  => mb_rhs(nb)%fld
    sxyz => mb_sxyz(nb)%fld
    vsl  => mb_vsl(nb)%fld
    vst  => mb_vst(nb)%fld

    ! I-direction
    fsfs => mb_fsf(nb,1)%fld
    fsfe => mb_fsf(nb,2)%fld

    stn = 1 - ngn
    edn = ni + ngn
    ste = -nge
    ede = ni + nge
!$OMP parallel private(vn,ve,fn,fc,dn,k,j,i,nfs,nfe,stn0,ste0,edn0,ede0,vis,kcp,kx,ky,kz,vx,vy,vz,f,m,der)
    allocate(vn(stn:edn,1:20), stat=ierr)
    allocate(ve(stn:edn,1:20), stat=ierr)
    allocate(fn(stn:edn,1:neqn), stat=ierr)
    allocate(fc(stn:edn,1:neqn), stat=ierr)
    allocate(dn(stn:edn,1:neqn), stat=ierr)
!$OMP do
    do k=nkst,nked
    do j=1,nj
        nfs = fsfs(nfsf_vis_d2int)%i3d( 1,j,k)
        nfe = fsfe(nfsf_vis_d2int)%i3d(ni,j,k)
        if (nfs == nbc_inter_scheme .or. nfs == nbc_intbc_scheme) then
            stn0 = stn
            ste0 = ste
        else
            stn0 = 1
            ste0 = 0
        end if

        if (nfe == nbc_inter_scheme .or. nfe == nbc_intbc_scheme) then
            edn0 = edn
            ede0 = ede
        else
            edn0 = ni
            ede0 = ni
        end if

        do i=stn0,edn0
            vis = vsl(1)%r3d(i,j,k)
            kcp = vis*cp_prl

            vis = vis + vst(1)%r3d(i,j,k)
            kcp = kcp + vst(1)%r3d(i,j,k)*cp_prt

            vn(i,1) = sxyz(1)%r3d(i,j,k)
            vn(i,2) = sxyz(2)%r3d(i,j,k)
            vn(i,3) = sxyz(3)%r3d(i,j,k)

            do m=4,15
                vn(i,m) = dpv(m-3)%r3d(i,j,k)
            end do

            vn(i,16) = vis
            vn(i,17) = kcp

            vn(i,18) = pv(2)%r3d(i,j,k)
            vn(i,19) = pv(3)%r3d(i,j,k)
            vn(i,20) = pv(4)%r3d(i,j,k)
        end do

        call sub_intplt(ni,1,20,ngn,vn,nfs,nfe,nge,ve)

        do i=ste0,ede0
            kx = ve(i,1)
            ky = ve(i,2)
            kz = ve(i,3)

            do m=1,12
                der(m) = ve(i,m+3)
            end do

            vis = ve(i,16)
            kcp = ve(i,17)

            vx = ve(i,18)
            vy = ve(i,19)
            vz = ve(i,20)

            call flux_vis(vx,vy,vz,nt,kx,ky,kz,kcp,vis,1,12,der,1,neqn,f)

            do m=1,neqn
                fc(i,m) = f(m)
            end do
        end do

        do i=stn0,edn0
            kx = vn(i,1)
            ky = vn(i,2)
            kz = vn(i,3)

            do m=1,12
                der(m) = vn(i,m+3)
            end do

            vis = vn(i,16)
            kcp = vn(i,17)

            vx = vn(i,18)
            vy = vn(i,19)
            vz = vn(i,20)

            call flux_vis(vx,vy,vz,nt,kx,ky,kz,kcp,vis,1,12,der,1,neqn,f)

            do m=1,neqn
                fn(i,m) = f(m)
            end do
        end do

        nfs = fsfs(nfsf_vis_d2der)%i3d( 1,j,k)
        nfe = fsfe(nfsf_vis_d2der)%i3d(ni,j,k)
        call sub_scheme(ni,1,neqn,ngn,fn,nge,fc,nfs,nfe,dn)

        do m=1,neqn
            do i=1,ni
                rhs(m)%r3d(i,j,k) = rhs(m)%r3d(i,j,k) + re*dn(i,m)
            end do
        end do
    end do
    end do
!$OMP end do nowait
    deallocate(dn, stat=ierr)
    deallocate(fc, stat=ierr)
    deallocate(fn, stat=ierr)
    deallocate(ve, stat=ierr)
    deallocate(vn, stat=ierr)
!$OMP end parallel



    ! J-direction
    fsfs => mb_fsf(nb,3)%fld
    fsfe => mb_fsf(nb,4)%fld

    stn = 1 - ngn
    edn = nj + ngn
    ste = -nge
    ede = nj + nge
!$OMP parallel private(vn,ve,fn,fc,dn,k,j,i,nfs,nfe,stn0,ste0,edn0,ede0,vis,kcp,ex,ey,ez,vx,vy,vz,f,m,der)
    allocate(vn(stn:edn,1:20), stat=ierr)
    allocate(ve(stn:edn,1:20), stat=ierr)
    allocate(fn(stn:edn,1:neqn), stat=ierr)
    allocate(fc(stn:edn,1:neqn), stat=ierr)
    allocate(dn(stn:edn,1:neqn), stat=ierr)
!$OMP do
    do k=nkst,nked
    do i=1,ni
        nfs = fsfs(nfsf_vis_d2int)%i3d(i, 1,k)
        nfe = fsfe(nfsf_vis_d2int)%i3d(i,nj,k)
        if (nfs == nbc_inter_scheme .or. nfs == nbc_intbc_scheme) then
            stn0 = stn
            ste0 = ste
        else
            stn0 = 1
            ste0 = 0
        end if

        if (nfe == nbc_inter_scheme .or. nfe == nbc_intbc_scheme) then
            edn0 = edn
            ede0 = ede
        else
            edn0 = nj
            ede0 = nj
        end if

        do j=stn0,edn0
            vis = vsl(1)%r3d(i,j,k)
            kcp = vis*cp_prl

            vis = vis + vst(1)%r3d(i,j,k)
            kcp = kcp + vst(1)%r3d(i,j,k)*cp_prt

            vn(j,1) = sxyz(4)%r3d(i,j,k)
            vn(j,2) = sxyz(5)%r3d(i,j,k)
            vn(j,3) = sxyz(6)%r3d(i,j,k)

            do m=4,15
                vn(j,m) = dpv(m-3)%r3d(i,j,k)
            end do

            vn(j,16) = vis
            vn(j,17) = kcp

            vn(j,18) = pv(2)%r3d(i,j,k)
            vn(j,19) = pv(3)%r3d(i,j,k)
            vn(j,20) = pv(4)%r3d(i,j,k)
        end do

        call sub_intplt(nj,1,20,ngn,vn,nfs,nfe,nge,ve)

        do j=ste0,ede0
            ex = ve(j,1)
            ey = ve(j,2)
            ez = ve(j,3)

            do m=1,12
                der(m) = ve(j,m+3)
            end do

            vis = ve(j,16)
            kcp = ve(j,17)

            vx = ve(j,18)
            vy = ve(j,19)
            vz = ve(j,20)

            call flux_vis(vx,vy,vz,nt,ex,ey,ez,kcp,vis,1,12,der,1,neqn,f)

            do m=1,neqn
                fc(j,m) = f(m)
            end do
        end do

        do j=stn0,edn0
            ex = vn(j,1)
            ey = vn(j,2)
            ez = vn(j,3)

            do m=1,12
                der(m) = vn(j,m+3)
            end do

            vis = vn(j,16)
            kcp = vn(j,17)

            vx = vn(j,18)
            vy = vn(j,19)
            vz = vn(j,20)

            call flux_vis(vx,vy,vz,nt,ex,ey,ez,kcp,vis,1,12,der,1,neqn,f)

            do m=1,neqn
                fn(j,m) = f(m)
            end do
        end do

        nfs = fsfs(nfsf_vis_d2der)%i3d(i, 1,k)
        nfe = fsfe(nfsf_vis_d2der)%i3d(i,nj,k)
        call sub_scheme(nj,1,neqn,ngn,fn,nge,fc,nfs,nfe,dn)

        do m=1,neqn
            do j=1,nj
                rhs(m)%r3d(i,j,k) = rhs(m)%r3d(i,j,k) + re*dn(j,m)
            end do
        end do
    end do
    end do
!$OMP end do
    deallocate(dn, stat=ierr)
    deallocate(fc, stat=ierr)
    deallocate(fn, stat=ierr)
    deallocate(ve, stat=ierr)
    deallocate(vn, stat=ierr)
!$OMP end parallel


    ! K-direction
    fsfs => mb_fsf(nb,5)%fld
    fsfe => mb_fsf(nb,6)%fld

    stn = 1 - ngn
    edn = nk + ngn
    ste = -nge
    ede = nk + nge
!$OMP parallel private(vn,ve,fn,fc,dn,k,j,i,nfs,nfe,stn0,ste0,edn0,ede0,vis,kcp,cx,cy,cz,vx,vy,vz,f,der,m)
    allocate(vn(stn:edn,1:20), stat=ierr)
    allocate(ve(stn:edn,1:20), stat=ierr)
    allocate(fn(stn:edn,1:neqn), stat=ierr)
    allocate(fc(stn:edn,1:neqn), stat=ierr)
    allocate(dn(stn:edn,1:neqn), stat=ierr)
!$OMP do
    do j=1,njed
    do i=1,ni
        nfs = fsfs(nfsf_vis_d2int)%i3d(i,j, 1)
        nfe = fsfe(nfsf_vis_d2int)%i3d(i,j,nk)
        if (nfs == nbc_inter_scheme .or. nfs == nbc_intbc_scheme) then
            stn0 = stn
            ste0 = ste
        else
            stn0 = 1
            ste0 = 0
        end if

        if (nfe == nbc_inter_scheme .or. nfe == nbc_intbc_scheme) then
            edn0 = edn
            ede0 = ede
        else
            edn0 = nk
            ede0 = nk
        end if

        do k=stn0,edn0
            vis = vsl(1)%r3d(i,j,k)
            kcp = vis*cp_prl

            vis = vis + vst(1)%r3d(i,j,k)
            kcp = kcp + vst(1)%r3d(i,j,k)*cp_prt

            vn(k,1) = sxyz(7)%r3d(i,j,k)
            vn(k,2) = sxyz(8)%r3d(i,j,k)
            vn(k,3) = sxyz(9)%r3d(i,j,k)

            do m=4,15
                vn(k,m) = dpv(m-3)%r3d(i,j,k)
            end do

            vn(k,16) = vis
            vn(k,17) = kcp

            vn(k,18) = pv(2)%r3d(i,j,k)
            vn(k,19) = pv(3)%r3d(i,j,k)
            vn(k,20) = pv(4)%r3d(i,j,k)
        end do

        call sub_intplt(nk,1,20,ngn,vn,nfs,nfe,nge,ve)

        do k=ste0,ede0
            cx = ve(k,1)
            cy = ve(k,2)
            cz = ve(k,3)

            do m=1,12
                der(m) = ve(k,m+3)
            end do

            vis = ve(k,16)
            kcp = ve(k,17)

            vx = ve(k,18)
            vy = ve(k,19)
            vz = ve(k,20)

            call flux_vis(vx,vy,vz,nt,cx,cy,cz,kcp,vis,1,12,der,1,neqn,f)

            do m=1,neqn
                fc(k,m) = f(m)
            end do
        end do

        do k=stn0,edn0
            cx = vn(k,1)
            cy = vn(k,2)
            cz = vn(k,3)

            do m=1,12
                der(m) = vn(k,m+3)
            end do

            vis = vn(k,16)
            kcp = vn(k,17)

            vx = vn(k,18)
            vy = vn(k,19)
            vz = vn(k,20)

            call flux_vis(vx,vy,vz,nt,cx,cy,cz,kcp,vis,1,12,der,1,neqn,f)

            do m=1,neqn
                fn(k,m) = f(m)
            end do
        end do

        nfs = fsfs(nfsf_vis_d2der)%i3d(i,j, 1)
        nfe = fsfe(nfsf_vis_d2der)%i3d(i,j,nk)
        call sub_scheme(nk,1,neqn,ngn,fn,nge,fc,nfs,nfe,dn)

        do m=1,neqn
            do k=1,nk
                rhs(m)%r3d(i,j,k) = rhs(m)%r3d(i,j,k) + re*dn(k,m)
            end do
        end do
    end do
    end do
!$OMP end do
    deallocate(dn, stat=ierr)
    deallocate(fc, stat=ierr)
    deallocate(fn, stat=ierr)
    deallocate(ve, stat=ierr)
    deallocate(vn, stat=ierr)
!$OMP end parallel

    if (nsw_kdir == nsw_dir_close) then
        do m=1,neqn
            do k=1,nk
            do j=1,nj
            do i=1,ni
                rhs(m)%r3d(i,j,k) = rhs(m)%r3d(i,j,nkst)
            end do
            end do
            end do
        end do
    end if

end subroutine calc_viscous

subroutine calc_viscous_sp(nb,ngn,nge,sub_intplt,sub_scheme)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,one
!    use mod_constants, only : nbc_bound_scheme,nsw_dir_close
    use mod_constants, only : nbc_inter_scheme,nsw_dir_close
    use mod_constants, only : nbc_intbc_scheme    
    use mod_constants, only : nfsf_vis_d2int,nfsf_vis_d2der
    use mod_datatypes, only : fld_array_t
    use mod_variables, only : nsw_kdir,gamma,prlam,prtur,refbeta,reue
    use mod_fieldvars, only : mb_topsp,npvs,mb_pv,neqn,mb_t,mb_sxyzsp
    use mod_fieldvars, only : mb_vsl,mb_vst,mb_rhs,mb_fsfsp,mb_dpv,mb_fsffp,mb_dpvfp
    use mod_openmp
    implicit none
    integer(kind_int), intent(in) :: nb,ngn,nge
    external                      :: sub_intplt,sub_scheme
    integer(kind_int)          :: i,j,k,m,ierr
    integer(kind_int)          :: ni,nj,nk,stn,edn,ste,ede
    integer(kind_int)          :: stn0,edn0,ste0,ede0,nfs,nfe
    integer(kind_int)          :: nkst,nked,njed
    integer(kind_int)          :: nr,bctype
    real(kind_real)            :: re,cp,gama,ae,vx,vy,vz,nt
    real(kind_real)            :: kx,ky,kz,ex,ey,ez,cx,cy,cz
    real(kind_real)            :: vis,kcp,cp_prl,cp_prt
    real(kind_real)            :: f(1:neqn),der(12),gradvars(12)
    real(kind_real), pointer   :: vn(:,:),ve(:,:)
    real(kind_real), pointer   :: fn(:,:),fc(:,:),dn(:,:)
    type(fld_array_t), pointer :: dpv(:),pv(:),rhs(:)
    type(fld_array_t), pointer :: sxyz(:),vsl(:),vst(:)
    type(fld_array_t), pointer :: fsfs(:),fsfe(:),fsfps(:),fsfpe(:),dpvfps(:),dpvfpe(:)

    gama = gamma
    ae = gama - one

    cp = gama*refbeta/ae
    cp_prl = cp/prlam
    cp_prt = cp/prtur

    re = one/reue

    ni = mb_topsp(nb)%nijk(1)
    nj = mb_topsp(nb)%nijk(2)
    nk = mb_topsp(nb)%nijk(3)

    nkst = mb_topsp(nb)%ndst(3)
    nked = mb_topsp(nb)%nded(3)
    if (nsw_kdir == nsw_dir_close) then
       njed = -1
    else
       njed = nj
    end if

    dpv  => mb_dpv(nb)%fld
    pv   => mb_pv(nb)%fld
    rhs  => mb_rhs(nb)%fld
    sxyz => mb_sxyzsp(nb)%fld
    vsl  => mb_vsl(nb)%fld
    vst  => mb_vst(nb)%fld

    ! I-direction
    fsfs => mb_fsfsp(nb,1)%fld
    fsfe => mb_fsfsp(nb,2)%fld
    fsfps=> mb_fsffp(nb,1)%fld    
    fsfpe=> mb_fsffp(nb,2)%fld 
    dpvfps => mb_dpvfp(nb,1)%fld
    dpvfpe => mb_dpvfp(nb,2)%fld    

    stn = 1 - ngn
    edn = ni + ngn
    ste = 1 - nge
    ede = ni + 1 + nge
!$OMP parallel private(vn,ve,fn,fc,dn,k,j,i,nfs,nfe,stn0,ste0,edn0,ede0,vis,kcp,kx,ky,kz,vx,vy,vz,f,m,der)
    allocate(vn(stn:edn,1:20), stat=ierr)
    allocate(ve(stn:edn,1:20), stat=ierr)
    allocate(fn(stn:edn,1:neqn), stat=ierr)
    allocate(fc(stn:edn,1:neqn), stat=ierr)
    allocate(dn(stn:edn,1:neqn), stat=ierr)
!$OMP do
    do k=nkst,nked
    do j=1,nj
        nfs = fsfs(nfsf_vis_d2int)%i3d( 1,j,k)
        nfe = fsfe(nfsf_vis_d2int)%i3d(ni,j,k)
        if (nfs == nbc_inter_scheme .or. nfs == nbc_intbc_scheme) then
            stn0 = stn
            ste0 = ste
        else
            stn0 = 1
            ste0 = 1
        end if

        if (nfe == nbc_inter_scheme .or. nfe == nbc_intbc_scheme) then
            edn0 = edn
            ede0 = ede
        else
            edn0 = ni
            ede0 = ni+1
        end if

        do i=stn0,edn0
            vis = vsl(1)%r3d(i,j,k)
            kcp = vis*cp_prl

            vis = vis + vst(1)%r3d(i,j,k)
            kcp = kcp + vst(1)%r3d(i,j,k)*cp_prt

            vn(i,1) = sxyz(1)%r3d(i,j,k)
            vn(i,2) = sxyz(2)%r3d(i,j,k)
            vn(i,3) = sxyz(3)%r3d(i,j,k)

            do m=4,15
                vn(i,m) = dpv(m-3)%r3d(i,j,k)
            end do

            vn(i,16) = vis
            vn(i,17) = kcp

            vn(i,18) = pv(2)%r3d(i,j,k)
            vn(i,19) = pv(3)%r3d(i,j,k)
            vn(i,20) = pv(4)%r3d(i,j,k)
        end do

        call sub_intplt(ni,1,20,ngn,vn,nfs,nfe,nge,ve)

        do i=ste0,ede0
            kx = ve(i,1)
            ky = ve(i,2)
            kz = ve(i,3)

            do m=1,12
                der(m) = ve(i,m+3)
            end do

            vis = ve(i,16)
            kcp = ve(i,17)

            vx = ve(i,18)
            vy = ve(i,19)
            vz = ve(i,20)
            
            if(i==1) then
                nr     = fsfps(1)%i3d(1,j,k)
                bctype = fsfps(2)%i3d(1,j,k)
                do m=1,12
                    dpvfps(m)%r3d(i,j,k) = ve(i,m+3)
                end do            
                if(bctype == 2) then
                    call set_boundary_gradvars(1,12,nb,nr,bctype,i,j,k,kx,ky,kz,gradvars)
                    do m=1,12
                        dpvfps(m)%r3d(i,j,k) = gradvars(m)
                        ve(i,m+3) = gradvars(m)
                    end do
                    ve(i,18) = 0.0
                    ve(i,19) = 0.0
                    ve(i,20) = 0.0
                end if                        
            else if(i==ni+1) then
                nr     = fsfpe(1)%i3d(ni+1,j,k)
                bctype = fsfpe(2)%i3d(ni+1,j,k)
                do m=1,12
                    dpvfpe(m)%r3d(i,j,k) = ve(i,m+3)
                end do            
                if(bctype == 2) then
                    call set_boundary_gradvars(1,12,nb,nr,bctype,i,j,k,kx,ky,kz,gradvars)
                    do m=1,12
                        dpvfpe(m)%r3d(i,j,k) = gradvars(m)
                        ve(i,m+3) = gradvars(m)
                    end do
                    ve(i,18) = 0.0
                    ve(i,19) = 0.0
                    ve(i,20) = 0.0
                end if             
            end if            

            call flux_vis(vx,vy,vz,nt,kx,ky,kz,kcp,vis,1,12,der,1,neqn,f)

            do m=1,neqn
                fc(i,m) = f(m)
            end do
        end do

        do i=stn0,edn0
            kx = vn(i,1)
            ky = vn(i,2)
            kz = vn(i,3)

            do m=1,12
                der(m) = vn(i,m+3)
            end do

            vis = vn(i,16)
            kcp = vn(i,17)

            vx = vn(i,18)
            vy = vn(i,19)
            vz = vn(i,20)

            call flux_vis(vx,vy,vz,nt,kx,ky,kz,kcp,vis,1,12,der,1,neqn,f)

            do m=1,neqn
                fn(i,m) = f(m)
            end do
        end do

        nfs = fsfs(nfsf_vis_d2der)%i3d( 1,j,k)
        nfe = fsfe(nfsf_vis_d2der)%i3d(ni,j,k)
        call sub_scheme(ni,1,neqn,ngn,fn,nge,fc,nfs,nfe,dn)

        do m=1,neqn
            do i=1,ni
                rhs(m)%r3d(i,j,k) = rhs(m)%r3d(i,j,k) + re*dn(i,m)
            end do
        end do
    end do
    end do
!$OMP end do nowait
    deallocate(dn, stat=ierr)
    deallocate(fc, stat=ierr)
    deallocate(fn, stat=ierr)
    deallocate(ve, stat=ierr)
    deallocate(vn, stat=ierr)
!$OMP end parallel



    ! J-direction
    fsfs => mb_fsfsp(nb,3)%fld
    fsfe => mb_fsfsp(nb,4)%fld
    fsfps=> mb_fsffp(nb,3)%fld    
    fsfpe=> mb_fsffp(nb,4)%fld 
    dpvfps => mb_dpvfp(nb,3)%fld
    dpvfpe => mb_dpvfp(nb,4)%fld     

    stn = 1 - ngn
    edn = nj + ngn
    ste = 1 - nge
    ede = nj + 1 + nge
!$OMP parallel private(vn,ve,fn,fc,dn,k,j,i,nfs,nfe,stn0,ste0,edn0,ede0,vis,kcp,ex,ey,ez,vx,vy,vz,f,m,der)
    allocate(vn(stn:edn,1:20), stat=ierr)
    allocate(ve(stn:edn,1:20), stat=ierr)
    allocate(fn(stn:edn,1:neqn), stat=ierr)
    allocate(fc(stn:edn,1:neqn), stat=ierr)
    allocate(dn(stn:edn,1:neqn), stat=ierr)
!$OMP do
    do k=nkst,nked
    do i=1,ni
        nfs = fsfs(nfsf_vis_d2int)%i3d(i, 1,k)
        nfe = fsfe(nfsf_vis_d2int)%i3d(i,nj,k)
        if (nfs == nbc_inter_scheme .or. nfs == nbc_intbc_scheme) then
            stn0 = stn
            ste0 = ste
        else
            stn0 = 1
            ste0 = 1
        end if

        if (nfe == nbc_inter_scheme .or. nfe == nbc_intbc_scheme) then
            edn0 = edn
            ede0 = ede
        else
            edn0 = nj
            ede0 = nj+1
        end if

        do j=stn0,edn0
            vis = vsl(1)%r3d(i,j,k)
            kcp = vis*cp_prl

            vis = vis + vst(1)%r3d(i,j,k)
            kcp = kcp + vst(1)%r3d(i,j,k)*cp_prt

            vn(j,1) = sxyz(4)%r3d(i,j,k)
            vn(j,2) = sxyz(5)%r3d(i,j,k)
            vn(j,3) = sxyz(6)%r3d(i,j,k)

            do m=4,15
                vn(j,m) = dpv(m-3)%r3d(i,j,k)
            end do

            vn(j,16) = vis
            vn(j,17) = kcp

            vn(j,18) = pv(2)%r3d(i,j,k)
            vn(j,19) = pv(3)%r3d(i,j,k)
            vn(j,20) = pv(4)%r3d(i,j,k)
        end do

        call sub_intplt(nj,1,20,ngn,vn,nfs,nfe,nge,ve)       

        do j=ste0,ede0
            ex = ve(j,1)
            ey = ve(j,2)
            ez = ve(j,3)

            do m=1,12
                der(m) = ve(j,m+3)
            end do

            vis = ve(j,16)
            kcp = ve(j,17)

            vx = ve(j,18)
            vy = ve(j,19)
            vz = ve(j,20)
            
            if(j==1) then
                nr     = fsfps(1)%i3d(i,1,k)
                bctype = fsfps(2)%i3d(i,1,k)
                do m=1,12
                    dpvfps(m)%r3d(i,j,k) = ve(j,m+3)
                end do            
                if(bctype == 2) then
                    call set_boundary_gradvars(1,12,nb,nr,bctype,i,j,k,ex,ey,ez,gradvars)
                    do m=1,12
                        dpvfps(m)%r3d(i,j,k) = gradvars(m)
                        ve(j,m+3) = gradvars(m)
                    end do
                    ve(j,18) = 0.0
                    ve(j,19) = 0.0
                    ve(j,20) = 0.0                    
                end if                        
            else if(j==nj+1) then
                nr     = fsfpe(1)%i3d(i,nj+1,k)
                bctype = fsfpe(2)%i3d(i,nj+1,k)
                do m=1,12
                    dpvfpe(m)%r3d(i,j,k) = ve(j,m+3)
                end do            
                if(bctype == 2) then
                    call set_boundary_gradvars(1,12,nb,nr,bctype,i,j,k,ex,ey,ez,gradvars)
                    do m=1,12
                        dpvfpe(m)%r3d(i,j,k) = gradvars(m)
                        ve(j,m+3) = gradvars(m)
                    end do
                    ve(j,18) = 0.0
                    ve(j,19) = 0.0
                    ve(j,20) = 0.0                    
                end if             
            end if             

            call flux_vis(vx,vy,vz,nt,ex,ey,ez,kcp,vis,1,12,der,1,neqn,f)

            do m=1,neqn
                fc(j,m) = f(m)
            end do
        end do

        do j=stn0,edn0
            ex = vn(j,1)
            ey = vn(j,2)
            ez = vn(j,3)

            do m=1,12
                der(m) = vn(j,m+3)
            end do

            vis = vn(j,16)
            kcp = vn(j,17)

            vx = vn(j,18)
            vy = vn(j,19)
            vz = vn(j,20)

            call flux_vis(vx,vy,vz,nt,ex,ey,ez,kcp,vis,1,12,der,1,neqn,f)

            do m=1,neqn
                fn(j,m) = f(m)
            end do
        end do

        nfs = fsfs(nfsf_vis_d2der)%i3d(i, 1,k)
        nfe = fsfe(nfsf_vis_d2der)%i3d(i,nj,k)
        call sub_scheme(nj,1,neqn,ngn,fn,nge,fc,nfs,nfe,dn)

        do m=1,neqn
            do j=1,nj
                rhs(m)%r3d(i,j,k) = rhs(m)%r3d(i,j,k) + re*dn(j,m)
            end do
        end do
    end do
    end do
!$OMP end do
    deallocate(dn, stat=ierr)
    deallocate(fc, stat=ierr)
    deallocate(fn, stat=ierr)
    deallocate(ve, stat=ierr)
    deallocate(vn, stat=ierr)
!$OMP end parallel


    ! K-direction
    fsfs => mb_fsfsp(nb,5)%fld
    fsfe => mb_fsfsp(nb,6)%fld
    fsfps=> mb_fsffp(nb,5)%fld    
    fsfpe=> mb_fsffp(nb,6)%fld 
    dpvfps => mb_dpvfp(nb,5)%fld
    dpvfpe => mb_dpvfp(nb,6)%fld     

    stn = 1 - ngn
    edn = nk + ngn
    ste = 1 - nge
    ede = nk + 1 + nge
!$OMP parallel private(vn,ve,fn,fc,dn,k,j,i,nfs,nfe,stn0,ste0,edn0,ede0,vis,kcp,cx,cy,cz,vx,vy,vz,f,der,m)
    allocate(vn(stn:edn,1:20), stat=ierr)
    allocate(ve(stn:edn,1:20), stat=ierr)
    allocate(fn(stn:edn,1:neqn), stat=ierr)
    allocate(fc(stn:edn,1:neqn), stat=ierr)
    allocate(dn(stn:edn,1:neqn), stat=ierr)
!$OMP do
    do j=1,njed
    do i=1,ni
        nfs = fsfs(nfsf_vis_d2int)%i3d(i,j, 1)
        nfe = fsfe(nfsf_vis_d2int)%i3d(i,j,nk)
        if (nfs == nbc_inter_scheme .or. nfs == nbc_intbc_scheme) then
            stn0 = stn
            ste0 = ste
        else
            stn0 = 1
            ste0 = 1
        end if

        if (nfe == nbc_inter_scheme .or. nfe == nbc_intbc_scheme) then
            edn0 = edn
            ede0 = ede
        else
            edn0 = nk
            ede0 = nk+1
        end if

        do k=stn0,edn0
            vis = vsl(1)%r3d(i,j,k)
            kcp = vis*cp_prl

            vis = vis + vst(1)%r3d(i,j,k)
            kcp = kcp + vst(1)%r3d(i,j,k)*cp_prt

            vn(k,1) = sxyz(7)%r3d(i,j,k)
            vn(k,2) = sxyz(8)%r3d(i,j,k)
            vn(k,3) = sxyz(9)%r3d(i,j,k)

            do m=4,15
                vn(k,m) = dpv(m-3)%r3d(i,j,k)
            end do

            vn(k,16) = vis
            vn(k,17) = kcp

            vn(k,18) = pv(2)%r3d(i,j,k)
            vn(k,19) = pv(3)%r3d(i,j,k)
            vn(k,20) = pv(4)%r3d(i,j,k)
        end do

        call sub_intplt(nk,1,20,ngn,vn,nfs,nfe,nge,ve)

        do k=ste0,ede0
            cx = ve(k,1)
            cy = ve(k,2)
            cz = ve(k,3)

            do m=1,12
                der(m) = ve(k,m+3)
            end do

            vis = ve(k,16)
            kcp = ve(k,17)

            vx = ve(k,18)
            vy = ve(k,19)
            vz = ve(k,20)
            
            if(k==1) then
                nr     = fsfps(1)%i3d(i,j,1)
                bctype = fsfps(2)%i3d(i,j,1)
                do m=1,12
                    dpvfps(m)%r3d(i,j,k) = ve(k,m+3)
                end do            
                if(bctype == 2) then
                    call set_boundary_gradvars(1,12,nb,nr,bctype,i,j,k,cx,cy,cz,gradvars)
                    do m=1,12
                        dpvfps(m)%r3d(i,j,k) = gradvars(m)
                        ve(k,m+3) = gradvars(m)
                    end do
                    ve(k,18) = 0.0
                    ve(k,19) = 0.0
                    ve(k,20) = 0.0                    
                end if                        
            else if(k==nk+1) then
                nr     = fsfpe(1)%i3d(i,j,nk+1)
                bctype = fsfpe(2)%i3d(i,j,nk+1)
                do m=1,12
                    dpvfpe(m)%r3d(i,j,k) = ve(k,m+3)
                end do            
                if(bctype == 2) then
                    call set_boundary_gradvars(1,12,nb,nr,bctype,i,j,k,cx,cy,cz,gradvars)
                    do m=1,12
                        dpvfpe(m)%r3d(i,j,k) = gradvars(m)
                        ve(k,m+3) = gradvars(m)
                    end do
                    ve(k,18) = 0.0
                    ve(k,19) = 0.0
                    ve(k,20) = 0.0                    
                end if             
            end if            

            call flux_vis(vx,vy,vz,nt,cx,cy,cz,kcp,vis,1,12,der,1,neqn,f)

            do m=1,neqn
                fc(k,m) = f(m)
            end do
        end do

        do k=stn0,edn0
            cx = vn(k,1)
            cy = vn(k,2)
            cz = vn(k,3)

            do m=1,12
                der(m) = vn(k,m+3)
            end do

            vis = vn(k,16)
            kcp = vn(k,17)

            vx = vn(k,18)
            vy = vn(k,19)
            vz = vn(k,20)

            call flux_vis(vx,vy,vz,nt,cx,cy,cz,kcp,vis,1,12,der,1,neqn,f)

            do m=1,neqn
                fn(k,m) = f(m)
            end do
        end do

        nfs = fsfs(nfsf_vis_d2der)%i3d(i,j, 1)
        nfe = fsfe(nfsf_vis_d2der)%i3d(i,j,nk)
        call sub_scheme(nk,1,neqn,ngn,fn,nge,fc,nfs,nfe,dn)

        do m=1,neqn
            do k=1,nk
                rhs(m)%r3d(i,j,k) = rhs(m)%r3d(i,j,k) + re*dn(k,m)
            end do
        end do
    end do
    end do
!$OMP end do
    deallocate(dn, stat=ierr)
    deallocate(fc, stat=ierr)
    deallocate(fn, stat=ierr)
    deallocate(ve, stat=ierr)
    deallocate(vn, stat=ierr)
!$OMP end parallel

    if (nsw_kdir == nsw_dir_close) then
        do m=1,neqn
            do k=1,nk
            do j=1,nj
            do i=1,ni
                rhs(m)%r3d(i,j,k) = rhs(m)%r3d(i,j,nkst)
            end do
            end do
            end do
        end do
    end if

end subroutine calc_viscous_sp

