
module mod_interface

    implicit none

    interface mb_var_create
        subroutine mb_var_create(mb_var,nst,ned,ngh)
            use mod_kndconsts, only : kind_int
            use mod_datatypes, only : var_block_t
            implicit none
            type(var_block_t), pointer, intent(inout) :: mb_var(:)
            integer(kind_int)         , intent(in)    :: nst,ned,ngh
        end subroutine mb_var_create
    end interface mb_var_create
    
    interface mb_var_dg_create
        subroutine mb_var_dg_create(mb_var,nst,ned,ngh)
            use mod_kndconsts, only : kind_int
            use mod_datatypes, only : var_block_t
            implicit none
            type(var_block_t), pointer, intent(inout) :: mb_var(:)
            integer(kind_int)         , intent(in)    :: nst,ned,ngh
        end subroutine mb_var_dg_create
    end interface mb_var_dg_create      
    
    interface mb_var_cc_create
        subroutine mb_var_cc_create(mb_var,nst,ned,ngh)
            use mod_kndconsts, only : kind_int
            use mod_datatypes, only : var_block_t
            implicit none
            type(var_block_t), pointer, intent(inout) :: mb_var(:)
            integer(kind_int)         , intent(in)    :: nst,ned,ngh
        end subroutine mb_var_cc_create
    end interface mb_var_cc_create  
    
    interface mb_var_create_sp
        subroutine mb_var_create_sp(mb_var,nst,ned,ngh)
            use mod_kndconsts, only : kind_int
            use mod_datatypes, only : var_block_t
            implicit none
            type(var_block_t), pointer, intent(inout) :: mb_var(:)
            integer(kind_int)         , intent(in)    :: nst,ned,ngh
        end subroutine mb_var_create_sp
    end interface mb_var_create_sp

    interface mb_var_delete
        subroutine mb_var_delete(mb_var)
            use mod_datatypes, only : var_block_t
            implicit none
            type(var_block_t), pointer, intent(inout) :: mb_var(:)
        end subroutine mb_var_delete
    end interface mb_var_delete

    interface exchange_bc_var
        subroutine exchange_bc_var(mb_var,nst,ned,ngh,nswmem,nswave)
            use mod_kndconsts, only : kind_int
            use mod_datatypes, only : var_block_t
            implicit none
            type(var_block_t), pointer, intent(in) :: mb_var(:)
            integer(kind_int)         , intent(in) :: nst,ned,ngh
            integer(kind_int)         , intent(in) :: nswmem
            integer(kind_int)         , intent(in) :: nswave
        end subroutine exchange_bc_var
    end interface exchange_bc_var
    interface pre_exchange_bc_var
        subroutine pre_exchange_bc_var(mb_var,nst,ned,ngh,nswmem,nswave)
            use mod_kndconsts, only : kind_int
            use mod_datatypes, only : var_block_t
            implicit none
            type(var_block_t), pointer             :: mb_var(:)
            integer(kind_int)         , intent(in) :: nst,ned,ngh
            integer(kind_int)         , intent(in) :: nswmem
            integer(kind_int)         , intent(in) :: nswave
        end subroutine pre_exchange_bc_var
    end interface pre_exchange_bc_var
    
    interface pre_exchange_bc_var_sp
        subroutine pre_exchange_bc_var_sp(mb_var,nst,ned,ngh,nswmem,nswave)
            use mod_kndconsts, only : kind_int
            use mod_datatypes, only : var_block_t
            implicit none
            type(var_block_t), pointer             :: mb_var(:)
            integer(kind_int)         , intent(in) :: nst,ned,ngh
            integer(kind_int)         , intent(in) :: nswmem
            integer(kind_int)         , intent(in) :: nswave
        end subroutine pre_exchange_bc_var_sp
    end interface pre_exchange_bc_var_sp    

    interface post_exchange_bc_var
        subroutine post_exchange_bc_var(mb_var,nst,ned,ngh,nswmem,nswave)
            use mod_kndconsts, only : kind_int
            use mod_datatypes, only : var_block_t
            implicit none
            type(var_block_t), pointer, intent(in) :: mb_var(:)
            integer(kind_int)         , intent(in) :: nst,ned,ngh
            integer(kind_int)         , intent(in) :: nswmem
            integer(kind_int)         , intent(in) :: nswave
        end subroutine post_exchange_bc_var
    end interface post_exchange_bc_var
    
    interface post_exchange_bc_var_sp
        subroutine post_exchange_bc_var_sp(mb_var,nst,ned,ngh,nswmem,nswave)
            use mod_kndconsts, only : kind_int
            use mod_datatypes, only : var_block_t
            implicit none
            type(var_block_t), pointer, intent(in) :: mb_var(:)
            integer(kind_int)         , intent(in) :: nst,ned,ngh
            integer(kind_int)         , intent(in) :: nswmem
            integer(kind_int)         , intent(in) :: nswave
        end subroutine post_exchange_bc_var_sp
    end interface post_exchange_bc_var_sp    
    
    interface pre_exchange_cc_bc_var
        subroutine pre_exchange_cc_bc_var(mb_var,nst,ned,ngh,nswmem,nswave)
            use mod_kndconsts, only : kind_int
            use mod_datatypes, only : var_block_t
            implicit none
            type(var_block_t), pointer             :: mb_var(:)
            integer(kind_int)         , intent(in) :: nst,ned,ngh
            integer(kind_int)         , intent(in) :: nswmem
            integer(kind_int)         , intent(in) :: nswave
        end subroutine pre_exchange_cc_bc_var
    end interface pre_exchange_cc_bc_var

    interface post_exchange_cc_bc_var
        subroutine post_exchange_cc_bc_var(mb_var,nst,ned,ngh,nswmem,nswave)
            use mod_kndconsts, only : kind_int
            use mod_datatypes, only : var_block_t
            implicit none
            type(var_block_t), pointer, intent(in) :: mb_var(:)
            integer(kind_int)         , intent(in) :: nst,ned,ngh
            integer(kind_int)         , intent(in) :: nswmem
            integer(kind_int)         , intent(in) :: nswave
        end subroutine post_exchange_cc_bc_var
    end interface post_exchange_cc_bc_var    

    interface average_bc_var
        subroutine average_bc_var(mb_var,nst,ned,nswmem,nswave)
            use mod_kndconsts, only : kind_int
            use mod_datatypes, only : var_block_t
            implicit none
            type(var_block_t), pointer, intent(in) :: mb_var(:)
            integer(kind_int)         , intent(in) :: nst,ned
            integer(kind_int)         , intent(in) :: nswmem
            integer(kind_int)         , intent(in) :: nswave
        end subroutine average_bc_var
    end interface average_bc_var
    
    interface average_bc_cc_var
        subroutine average_bc_cc_var(mb_var,nst,ned,nswmem,nswave)
            use mod_kndconsts, only : kind_int
            use mod_datatypes, only : var_block_t
            implicit none
            type(var_block_t), pointer, intent(in) :: mb_var(:)
            integer(kind_int)         , intent(in) :: nst,ned
            integer(kind_int)         , intent(in) :: nswmem
            integer(kind_int)         , intent(in) :: nswave
        end subroutine average_bc_cc_var
    end interface average_bc_cc_var    

    interface extending_ghost_points
        subroutine extending_ghost_points(mb_var,nst,ned,ngh)
            use mod_kndconsts, only : kind_int
            use mod_datatypes, only : var_block_t
            implicit none
            type(var_block_t), pointer, intent(in) :: mb_var(:)
            integer(kind_int)         , intent(in) :: nst,ned,ngh
        end subroutine extending_ghost_points
    end interface extending_ghost_points

    interface patched_ghost_points
        subroutine patched_ghost_points(mb_var,nst,ned,ngh)
            use mod_kndconsts, only : kind_int
            use mod_datatypes, only : var_block_t
            implicit none
            type(var_block_t), pointer, intent(in) :: mb_var(:)
            integer(kind_int)         , intent(in) :: nst,ned,ngh
        end subroutine patched_ghost_points
    end interface patched_ghost_points

    interface exchange_bc_der
        subroutine exchange_bc_der(mb_der,nst,ned,ngh,nswmem)
            use mod_kndconsts, only : kind_int
            use mod_datatypes, only : var_block_t
            implicit none
            type(var_block_t), pointer, intent(in) :: mb_der(:)
            integer(kind_int)         , intent(in) :: nst,ned
            integer(kind_int)         , intent(in) :: ngh
            integer(kind_int)         , intent(in) :: nswmem
        end subroutine exchange_bc_der
    end interface exchange_bc_der
    
    interface exchange_cc_bc_der
        subroutine exchange_cc_bc_der(mb_der,nst,ned,ngh,nswmem)
            use mod_kndconsts, only : kind_int
            use mod_datatypes, only : var_block_t
            implicit none
            type(var_block_t), pointer, intent(in) :: mb_der(:)
            integer(kind_int)         , intent(in) :: nst,ned
            integer(kind_int)         , intent(in) :: ngh
            integer(kind_int)         , intent(in) :: nswmem
        end subroutine exchange_cc_bc_der
    end interface exchange_cc_bc_der  
    
    interface exchange_sp_bc_der
        subroutine exchange_sp_bc_der(mb_der,nst,ned,ngh,nswmem)
            use mod_kndconsts, only : kind_int
            use mod_datatypes, only : var_block_t
            implicit none
            type(var_block_t), pointer, intent(in) :: mb_der(:)
            integer(kind_int)         , intent(in) :: nst,ned
            integer(kind_int)         , intent(in) :: ngh
            integer(kind_int)         , intent(in) :: nswmem
        end subroutine exchange_sp_bc_der
    end interface exchange_sp_bc_der     

    interface exchange_bc_sxyz
        subroutine exchange_bc_sxyz(mb_sxyz,nst,ned,ngh)
            use mod_kndconsts, only : kind_int
            use mod_datatypes, only : var_block_t
            implicit none
            type(var_block_t), pointer, intent(in) :: mb_sxyz(:)
            integer(kind_int)         , intent(in) :: nst,ned,ngh
        end subroutine exchange_bc_sxyz
    end interface exchange_bc_sxyz
    
    interface exchange_cc_bc_sxyz
        subroutine exchange_cc_bc_sxyz(mb_sxyzcc,nst,ned,ngh)
            use mod_kndconsts, only : kind_int
            use mod_datatypes, only : var_block_t
            implicit none
            type(var_block_t), pointer, intent(in) :: mb_sxyzcc(:)
            integer(kind_int)         , intent(in) :: nst,ned,ngh
        end subroutine exchange_cc_bc_sxyz
    end interface exchange_cc_bc_sxyz
    
    interface exchange_sp_bc_sxyz
        subroutine exchange_sp_bc_sxyz(mb_sxyzsp,nst,ned,ngh)
            use mod_kndconsts, only : kind_int
            use mod_datatypes, only : var_block_t
            implicit none
            type(var_block_t), pointer, intent(in) :: mb_sxyzsp(:)
            integer(kind_int)         , intent(in) :: nst,ned,ngh
        end subroutine exchange_sp_bc_sxyz
    end interface exchange_sp_bc_sxyz     

    interface scatter_input_mb_var
        subroutine scatter_input_mb_var(io_unit,mb_var,nst,ned,nfun)
            use mod_kndconsts, only : kind_int
            use mod_datatypes, only : var_block_t
            implicit none
            integer(kind_int),          intent(in) :: io_unit
            type(var_block_t), pointer, intent(in) :: mb_var(:)
            integer(kind_int),          intent(in) :: nst,ned,nfun
        end subroutine scatter_input_mb_var
    end interface scatter_input_mb_var
    
    interface scatter_input_mb_var_cc
        subroutine scatter_input_mb_var_cc(io_unit,mb_var,nst,ned,nfun)
            use mod_kndconsts, only : kind_int
            use mod_datatypes, only : var_block_t
            implicit none
            integer(kind_int),          intent(in) :: io_unit
            type(var_block_t), pointer, intent(in) :: mb_var(:)
            integer(kind_int),          intent(in) :: nst,ned,nfun
        end subroutine scatter_input_mb_var_cc
    end interface scatter_input_mb_var_cc    
    
    interface scatter_input_mb_var_sp
        subroutine scatter_input_mb_var_sp(io_unit,mb_var,nst,ned,nfun)
            use mod_kndconsts, only : kind_int
            use mod_datatypes, only : var_block_t
            implicit none
            integer(kind_int),          intent(in) :: io_unit
            type(var_block_t), pointer, intent(in) :: mb_var(:)
            integer(kind_int),          intent(in) :: nst,ned,nfun
        end subroutine scatter_input_mb_var_sp
    end interface scatter_input_mb_var_sp    

    interface gather_output_mb_var
        subroutine gather_output_mb_var(io_unit,mb_var,nst,ned,nfun)
            use mod_kndconsts, only : kind_int
            use mod_datatypes, only : var_block_t
            implicit none
            integer(kind_int),          intent(in) :: io_unit
            type(var_block_t), pointer, intent(in) :: mb_var(:)
            integer(kind_int),          intent(in) :: nst,ned,nfun
        end subroutine gather_output_mb_var
    end interface gather_output_mb_var
    
    interface gather_output_mb_var_sp
        subroutine gather_output_mb_var_sp(io_unit,mb_var,nst,ned,nfun)
            use mod_kndconsts, only : kind_int
            use mod_datatypes, only : var_block_t
            implicit none
            integer(kind_int),          intent(in) :: io_unit
            type(var_block_t), pointer, intent(in) :: mb_var(:)
            integer(kind_int),          intent(in) :: nst,ned,nfun
        end subroutine gather_output_mb_var_sp
    end interface gather_output_mb_var_sp

    interface calc_mb_dn_via_node3
        subroutine calc_mb_dn_via_node3(mb_vn,nvn, &
                                        ngn,sub_ve,nge,sub_dn, &
                                        nfsf_ve,nfsf_dn, &
                                        mb_dn,ndn,ndir)
            use mod_kndconsts, only : kind_int
            use mod_datatypes, only : var_block_t
            implicit none
            type(var_block_t), pointer, intent(in) :: mb_vn(:)
            integer(kind_int),          intent(in) :: nvn
            integer(kind_int),          intent(in) :: ngn
            external                               :: sub_ve
            integer(kind_int),          intent(in) :: nge
            external                               :: sub_dn
            integer(kind_int),          intent(in) :: nfsf_ve,nfsf_dn
            type(var_block_t), pointer, intent(in) :: mb_dn(:)
            integer(kind_int),          intent(in) :: ndn
            integer(kind_int),          intent(in) :: ndir
        end subroutine calc_mb_dn_via_node3
    end interface calc_mb_dn_via_node3

    interface calc_cc_mb_dn_via_node3
        subroutine calc_cc_mb_dn_via_node3(mb_vn,nvn, &
                                           ngn,sub_ve,nge,sub_dn, &
                                           nfsf_ve,nfsf_dn, &
                                           mb_dn,ndn,ndir)
            use mod_kndconsts, only : kind_int
            use mod_datatypes, only : var_block_t
            implicit none
            type(var_block_t), pointer, intent(in) :: mb_vn(:)
            integer(kind_int),          intent(in) :: nvn
            integer(kind_int),          intent(in) :: ngn
            external                               :: sub_ve
            integer(kind_int),          intent(in) :: nge
            external                               :: sub_dn
            integer(kind_int),          intent(in) :: nfsf_ve,nfsf_dn
            type(var_block_t), pointer, intent(in) :: mb_dn(:)
            integer(kind_int),          intent(in) :: ndn
            integer(kind_int),          intent(in) :: ndir
        end subroutine calc_cc_mb_dn_via_node3
    end interface calc_cc_mb_dn_via_node3
    
    interface calc_sp_mb_dn_via_node3
        subroutine calc_sp_mb_dn_via_node3(mb_vn,nvn, &
                                           ngn,sub_ve,nge,sub_dn, &
                                           nfsf_ve,nfsf_dn, &
                                           mb_dn,ndn,ndir)
            use mod_kndconsts, only : kind_int
            use mod_datatypes, only : var_block_t
            implicit none
            type(var_block_t), pointer, intent(in) :: mb_vn(:)
            integer(kind_int),          intent(in) :: nvn
            integer(kind_int),          intent(in) :: ngn
            external                               :: sub_ve
            integer(kind_int),          intent(in) :: nge
            external                               :: sub_dn
            integer(kind_int),          intent(in) :: nfsf_ve,nfsf_dn
            type(var_block_t), pointer, intent(in) :: mb_dn(:)
            integer(kind_int),          intent(in) :: ndn
            integer(kind_int),          intent(in) :: ndir
        end subroutine calc_sp_mb_dn_via_node3
    end interface calc_sp_mb_dn_via_node3    
    
    interface calc_mb_dn_via_node3_exp
        subroutine calc_mb_dn_via_node3_exp(mb_vn,nvn, &
                                             ngn,sub_ve,nge,sub_dn, &
                                             nfsf_ve,nfsf_dn, &
                                             mb_dn,ndn,ndir)
            use mod_kndconsts, only : kind_int
            use mod_datatypes, only : var_block_t
            implicit none
            type(var_block_t), pointer, intent(in) :: mb_vn(:)
            integer(kind_int),          intent(in) :: nvn
            integer(kind_int),          intent(in) :: ngn
            external                               :: sub_ve
            integer(kind_int),          intent(in) :: nge
            external                               :: sub_dn
            integer(kind_int),          intent(in) :: nfsf_ve,nfsf_dn
            type(var_block_t), pointer, intent(in) :: mb_dn(:)
            integer(kind_int),          intent(in) :: ndn
            integer(kind_int),          intent(in) :: ndir
        end subroutine calc_mb_dn_via_node3_exp
    end interface calc_mb_dn_via_node3_exp
    
    interface calc_mb_dn_via_node3_sp
        subroutine calc_mb_dn_via_node3_sp(mb_vn,nvn, &
                                           ngn,sub_ve,nge,sub_dn, &
                                           nfsf_ve,nfsf_dn, &
                                           mb_dn,ndn,ndir)
            use mod_kndconsts, only : kind_int
            use mod_datatypes, only : var_block_t
            implicit none
            type(var_block_t), pointer, intent(in) :: mb_vn(:)
            integer(kind_int),          intent(in) :: nvn
            integer(kind_int),          intent(in) :: ngn
            external                               :: sub_ve
            integer(kind_int),          intent(in) :: nge
            external                               :: sub_dn
            integer(kind_int),          intent(in) :: nfsf_ve,nfsf_dn
            type(var_block_t), pointer, intent(in) :: mb_dn(:)
            integer(kind_int),          intent(in) :: ndn
            integer(kind_int),          intent(in) :: ndir
        end subroutine calc_mb_dn_via_node3_sp
    end interface calc_mb_dn_via_node3_sp
    
    interface calc_mb_duvwt_via_node3_sp
        subroutine calc_mb_duvwt_via_node3_sp(mb_vn,nvn, &
                                              ngn,sub_ve,nge,sub_dn, &
                                              nfsf_ve,nfsf_dn, &
                                              mb_dn,ndn,ndir)
            use mod_kndconsts, only : kind_int
            use mod_datatypes, only : var_block_t
            implicit none
            type(var_block_t), pointer, intent(in) :: mb_vn(:)
            integer(kind_int),          intent(in) :: nvn
            integer(kind_int),          intent(in) :: ngn
            external                               :: sub_ve
            integer(kind_int),          intent(in) :: nge
            external                               :: sub_dn
            integer(kind_int),          intent(in) :: nfsf_ve,nfsf_dn
            type(var_block_t), pointer, intent(in) :: mb_dn(:)
            integer(kind_int),          intent(in) :: ndn
            integer(kind_int),          intent(in) :: ndir
        end subroutine calc_mb_duvwt_via_node3_sp
    end interface calc_mb_duvwt_via_node3_sp    

    interface calc_mb_dn_via_node2
        subroutine calc_mb_dn_via_node2(mb_vn1,nvn1, &
                                        mb_vn2,nvn2, &
                                        ngn,sub_ve,nge,sub_dn, &
                                        nfsf_ve,nfsf_dn, &
                                        mb_dn,ndn,ndir)
            use mod_kndconsts, only : kind_int
            use mod_datatypes, only : var_block_t
            implicit none
            type(var_block_t), pointer, intent(in) :: mb_vn1(:)
            integer(kind_int),          intent(in) :: nvn1
            type(var_block_t), pointer, intent(in) :: mb_vn2(:)
            integer(kind_int),          intent(in) :: nvn2
            integer(kind_int),          intent(in) :: ngn
            external                               :: sub_ve
            integer(kind_int),          intent(in) :: nge
            external                               :: sub_dn
            integer(kind_int),          intent(in) :: nfsf_ve,nfsf_dn
            type(var_block_t), pointer, intent(in) :: mb_dn(:)
            integer(kind_int),          intent(in) :: ndn
            integer(kind_int),          intent(in) :: ndir
        end subroutine calc_mb_dn_via_node2
    end interface calc_mb_dn_via_node2
    
    interface calc_cc_mb_dn_via_node2
        subroutine calc_cc_mb_dn_via_node2(mb_vn1,nvn1, &
                                           mb_vn2,nvn2, &
                                           ngn,sub_ve,nge,sub_dn, &
                                           nfsf_ve,nfsf_dn, &
                                           mb_dn,ndn,ndir)
            use mod_kndconsts, only : kind_int
            use mod_datatypes, only : var_block_t
            implicit none
            type(var_block_t), pointer, intent(in) :: mb_vn1(:)
            integer(kind_int),          intent(in) :: nvn1
            type(var_block_t), pointer, intent(in) :: mb_vn2(:)
            integer(kind_int),          intent(in) :: nvn2
            integer(kind_int),          intent(in) :: ngn
            external                               :: sub_ve
            integer(kind_int),          intent(in) :: nge
            external                               :: sub_dn
            integer(kind_int),          intent(in) :: nfsf_ve,nfsf_dn
            type(var_block_t), pointer, intent(in) :: mb_dn(:)
            integer(kind_int),          intent(in) :: ndn
            integer(kind_int),          intent(in) :: ndir
        end subroutine calc_cc_mb_dn_via_node2
    end interface calc_cc_mb_dn_via_node2
    
    interface calc_sp_mb_dn_via_node2
        subroutine calc_sp_mb_dn_via_node2(mb_vn1,nvn1, &
                                           mb_vn2,nvn2, &
                                           ngn,sub_ve,nge,sub_dn, &
                                           nfsf_ve,nfsf_dn, &
                                           mb_dn,ndn,ndir)
            use mod_kndconsts, only : kind_int
            use mod_datatypes, only : var_block_t
            implicit none
            type(var_block_t), pointer, intent(in) :: mb_vn1(:)
            integer(kind_int),          intent(in) :: nvn1
            type(var_block_t), pointer, intent(in) :: mb_vn2(:)
            integer(kind_int),          intent(in) :: nvn2
            integer(kind_int),          intent(in) :: ngn
            external                               :: sub_ve
            integer(kind_int),          intent(in) :: nge
            external                               :: sub_dn
            integer(kind_int),          intent(in) :: nfsf_ve,nfsf_dn
            type(var_block_t), pointer, intent(in) :: mb_dn(:)
            integer(kind_int),          intent(in) :: ndn
            integer(kind_int),          intent(in) :: ndir
        end subroutine calc_sp_mb_dn_via_node2
    end interface calc_sp_mb_dn_via_node2

    interface calc_mb_dn_via_node2_exp
        subroutine calc_mb_dn_via_node2_exp(mb_vn1,nvn1, &
                                             mb_vn2,nvn2, &
                                             ngn,sub_ve,nge,sub_dn, &
                                             nfsf_ve,nfsf_dn, &
                                             mb_dn,ndn,ndir)
            use mod_kndconsts, only : kind_int
            use mod_datatypes, only : var_block_t
            implicit none
            type(var_block_t), pointer, intent(in) :: mb_vn1(:)
            integer(kind_int),          intent(in) :: nvn1
            type(var_block_t), pointer, intent(in) :: mb_vn2(:)
            integer(kind_int),          intent(in) :: nvn2
            integer(kind_int),          intent(in) :: ngn
            external                               :: sub_ve
            integer(kind_int),          intent(in) :: nge
            external                               :: sub_dn
            integer(kind_int),          intent(in) :: nfsf_ve,nfsf_dn
            type(var_block_t), pointer, intent(in) :: mb_dn(:)
            integer(kind_int),          intent(in) :: ndn
            integer(kind_int),          intent(in) :: ndir
        end subroutine calc_mb_dn_via_node2_exp
    end interface calc_mb_dn_via_node2_exp

    interface calc_mb_dn_via_vec_node2
        subroutine calc_mb_dn_via_vec_node2(mb_vn1,nst1,ned1, &
                                            mb_vn2,nst2,ned2, &
                                            ngn,sub_ve,nge,sub_dn, &
                                            nfsf_ve,nfsf_dn, &
                                            mb_dn,ndn,ndir)
            use mod_kndconsts, only : kind_int
            use mod_datatypes, only : var_block_t
            implicit none
            type(var_block_t), pointer, intent(in) :: mb_vn1(:)
            integer(kind_int),          intent(in) :: nst1,ned1
            type(var_block_t), pointer, intent(in) :: mb_vn2(:)
            integer(kind_int),          intent(in) :: nst2,ned2
            integer(kind_int),          intent(in) :: ngn
            external                               :: sub_ve
            integer(kind_int),          intent(in) :: nge
            external                               :: sub_dn
            integer(kind_int),          intent(in) :: nfsf_ve,nfsf_dn
            type(var_block_t), pointer, intent(in) :: mb_dn(:)
            integer(kind_int),          intent(in) :: ndn
            integer(kind_int),          intent(in) :: ndir
        end subroutine calc_mb_dn_via_vec_node2
    end interface calc_mb_dn_via_vec_node2
    
    interface calc_cc_mb_dn_via_vec_node2
        subroutine calc_cc_mb_dn_via_vec_node2(mb_vn1,nst1,ned1, &
                                               mb_vn2,nst2,ned2, &
                                               ngn,sub_ve,nge,sub_dn, &
                                               nfsf_ve,nfsf_dn, &
                                               mb_dn,ndn,ndir)
            use mod_kndconsts, only : kind_int
            use mod_datatypes, only : var_block_t
            implicit none
            type(var_block_t), pointer, intent(in) :: mb_vn1(:)
            integer(kind_int),          intent(in) :: nst1,ned1
            type(var_block_t), pointer, intent(in) :: mb_vn2(:)
            integer(kind_int),          intent(in) :: nst2,ned2
            integer(kind_int),          intent(in) :: ngn
            external                               :: sub_ve
            integer(kind_int),          intent(in) :: nge
            external                               :: sub_dn
            integer(kind_int),          intent(in) :: nfsf_ve,nfsf_dn
            type(var_block_t), pointer, intent(in) :: mb_dn(:)
            integer(kind_int),          intent(in) :: ndn
            integer(kind_int),          intent(in) :: ndir
        end subroutine calc_cc_mb_dn_via_vec_node2
    end interface calc_cc_mb_dn_via_vec_node2
    
    interface calc_sp_mb_dn_via_vec_node2
        subroutine calc_sp_mb_dn_via_vec_node2(mb_vn1,nst1,ned1, &
                                               mb_vn2,nst2,ned2, &
                                               ngn,sub_ve,nge,sub_dn, &
                                               nfsf_ve,nfsf_dn, &
                                               mb_dn,ndn,ndir)
            use mod_kndconsts, only : kind_int
            use mod_datatypes, only : var_block_t
            implicit none
            type(var_block_t), pointer, intent(in) :: mb_vn1(:)
            integer(kind_int),          intent(in) :: nst1,ned1
            type(var_block_t), pointer, intent(in) :: mb_vn2(:)
            integer(kind_int),          intent(in) :: nst2,ned2
            integer(kind_int),          intent(in) :: ngn
            external                               :: sub_ve
            integer(kind_int),          intent(in) :: nge
            external                               :: sub_dn
            integer(kind_int),          intent(in) :: nfsf_ve,nfsf_dn
            type(var_block_t), pointer, intent(in) :: mb_dn(:)
            integer(kind_int),          intent(in) :: ndn
            integer(kind_int),          intent(in) :: ndir
        end subroutine calc_sp_mb_dn_via_vec_node2
    end interface calc_sp_mb_dn_via_vec_node2    

    interface calc_mb_sxyz_vol
        subroutine calc_mb_sxyz_vol(mb_der3,nst3,ned3, &
                                    mb_der2,nst2,ned2, &
                                    mb_sxyz,nst,ned, &
                                    mb_vol,nvol)
            use mod_kndconsts, only : kind_int
            use mod_datatypes, only : var_block_t
            implicit none
            type(var_block_t), pointer, intent(in) :: mb_der3(:)
            integer(kind_int),          intent(in) :: nst3,ned3
            type(var_block_t), pointer, intent(in) :: mb_der2(:)
            integer(kind_int),          intent(in) :: nst2,ned2
            type(var_block_t), pointer, intent(in) :: mb_sxyz(:)
            integer(kind_int),          intent(in) :: nst,ned
            type(var_block_t), pointer, intent(in) :: mb_vol(:)
            integer(kind_int),          intent(in) :: nvol
        end subroutine calc_mb_sxyz_vol
    end interface calc_mb_sxyz_vol
    
    interface calc_cc_mb_sxyz_vol
        subroutine calc_cc_mb_sxyz_vol(mb_der3,nst3,ned3, &
                                    mb_der2,nst2,ned2, &
                                    mb_sxyzsp,nst,ned, &
                                    mb_volsp,nvol)
            use mod_kndconsts, only : kind_int
            use mod_datatypes, only : var_block_t
            implicit none
            type(var_block_t), pointer, intent(in) :: mb_der3(:)
            integer(kind_int),          intent(in) :: nst3,ned3
            type(var_block_t), pointer, intent(in) :: mb_der2(:)
            integer(kind_int),          intent(in) :: nst2,ned2
            type(var_block_t), pointer, intent(in) :: mb_sxyzsp(:)
            integer(kind_int),          intent(in) :: nst,ned
            type(var_block_t), pointer, intent(in) :: mb_volsp(:)
            integer(kind_int),          intent(in) :: nvol
        end subroutine calc_cc_mb_sxyz_vol
    end interface calc_cc_mb_sxyz_vol
    
    interface calc_sp_mb_sxyz_vol
        subroutine calc_sp_mb_sxyz_vol(mb_der3,nst3,ned3, &
                                    mb_der2,nst2,ned2, &
                                    mb_sxyzcc,nst,ned, &
                                    mb_volcc,nvol)
            use mod_kndconsts, only : kind_int
            use mod_datatypes, only : var_block_t
            implicit none
            type(var_block_t), pointer, intent(in) :: mb_der3(:)
            integer(kind_int),          intent(in) :: nst3,ned3
            type(var_block_t), pointer, intent(in) :: mb_der2(:)
            integer(kind_int),          intent(in) :: nst2,ned2
            type(var_block_t), pointer, intent(in) :: mb_sxyzcc(:)
            integer(kind_int),          intent(in) :: nst,ned
            type(var_block_t), pointer, intent(in) :: mb_volcc(:)
            integer(kind_int),          intent(in) :: nvol
        end subroutine calc_sp_mb_sxyz_vol
    end interface calc_sp_mb_sxyz_vol    

    interface calc_mb_sxyz_vol_exp
        subroutine calc_mb_sxyz_vol_exp(mb_der3,nst3,ned3, &
                                         mb_der2,nst2,ned2, &
                                         mb_sxyz,nst,ned, &
                                         mb_vol,nvol)
            use mod_kndconsts, only : kind_int
            use mod_datatypes, only : var_block_t
            implicit none
            type(var_block_t), pointer, intent(in) :: mb_der3(:)
            integer(kind_int),          intent(in) :: nst3,ned3
            type(var_block_t), pointer, intent(in) :: mb_der2(:)
            integer(kind_int),          intent(in) :: nst2,ned2
            type(var_block_t), pointer, intent(in) :: mb_sxyz(:)
            integer(kind_int),          intent(in) :: nst,ned
            type(var_block_t), pointer, intent(in) :: mb_vol(:)
            integer(kind_int),          intent(in) :: nvol
        end subroutine calc_mb_sxyz_vol_exp
    end interface calc_mb_sxyz_vol_exp

    interface calc_mb_vol
        subroutine calc_mb_vol(mb_der1,nst1,ned1,mb_vol,nvol)
            use mod_kndconsts, only : kind_int
            use mod_datatypes, only : var_block_t
            implicit none
            type(var_block_t), pointer, intent(in) :: mb_der1(:)
            integer(kind_int),          intent(in) :: nst1,ned1
            type(var_block_t), pointer, intent(in) :: mb_vol(:)
            integer(kind_int),          intent(in) :: nvol
        end subroutine calc_mb_vol
    end interface calc_mb_vol
    
    interface calc_cc_mb_vol
        subroutine calc_cc_mb_vol(mb_der1,nst1,ned1,mb_volcc,nvol)
            use mod_kndconsts, only : kind_int
            use mod_datatypes, only : var_block_t
            implicit none
            type(var_block_t), pointer, intent(in) :: mb_der1(:)
            integer(kind_int),          intent(in) :: nst1,ned1
            type(var_block_t), pointer, intent(in) :: mb_volcc(:)
            integer(kind_int),          intent(in) :: nvol
        end subroutine calc_cc_mb_vol
    end interface calc_cc_mb_vol
    
    interface calc_sp_mb_vol
        subroutine calc_sp_mb_vol(mb_der1,nst1,ned1,mb_volsp,nvol)
            use mod_kndconsts, only : kind_int
            use mod_datatypes, only : var_block_t
            implicit none
            type(var_block_t), pointer, intent(in) :: mb_der1(:)
            integer(kind_int),          intent(in) :: nst1,ned1
            type(var_block_t), pointer, intent(in) :: mb_volsp(:)
            integer(kind_int),          intent(in) :: nvol
        end subroutine calc_sp_mb_vol
    end interface calc_sp_mb_vol    

    interface calc_mb_var_via_sub
        subroutine calc_mb_var_via_sub(mb_vin,nsin,nein,sub_via, &
                                       mb_vout,nsout,neout,ngh)
            use mod_kndconsts, only : kind_int
            use mod_datatypes, only : var_block_t
            implicit none
            type(var_block_t), pointer, intent(in) :: mb_vin(:)
            integer(kind_int)         , intent(in) :: nsin,nein
            external                               :: sub_via
            type(var_block_t), pointer, intent(in) :: mb_vout(:)
            integer(kind_int)         , intent(in) :: nsout,neout
            integer(kind_int)         , intent(in) :: ngh
        end subroutine calc_mb_var_via_sub
    end interface calc_mb_var_via_sub
    
    interface calc_mb_var_via_sub_sp
        subroutine calc_mb_var_via_sub_sp(mb_vin,nsin,nein,sub_via, &
                                       mb_vout,nsout,neout,ngh)
            use mod_kndconsts, only : kind_int
            use mod_datatypes, only : var_block_t
            implicit none
            type(var_block_t), pointer, intent(in) :: mb_vin(:)
            integer(kind_int)         , intent(in) :: nsin,nein
            external                               :: sub_via
            type(var_block_t), pointer, intent(in) :: mb_vout(:)
            integer(kind_int)         , intent(in) :: nsout,neout
            integer(kind_int)         , intent(in) :: ngh
        end subroutine calc_mb_var_via_sub_sp
    end interface calc_mb_var_via_sub_sp    

    interface calc_bc_var_via_sub
        subroutine calc_bc_var_via_sub(mb_vin,nsin,nein,sub_via, &
                                       mb_vout,nsout,neout,ngst,nged)
            use mod_kndconsts, only : kind_int
            use mod_datatypes, only : var_block_t
            implicit none
            type(var_block_t), pointer, intent(in) :: mb_vin(:)
            integer(kind_int)         , intent(in) :: nsin,nein
            external                               :: sub_via
            type(var_block_t), pointer, intent(in) :: mb_vout(:)
            integer(kind_int)         , intent(in) :: nsout,neout
            integer(kind_int)         , intent(in) :: ngst,nged
        end subroutine calc_bc_var_via_sub
    end interface calc_bc_var_via_sub
    
    interface calc_bc_var_via_sub_sp
        subroutine calc_bc_var_via_sub_sp(mb_vin,nsin,nein,sub_via, &
                                       mb_vout,nsout,neout,ngst,nged)
            use mod_kndconsts, only : kind_int
            use mod_datatypes, only : var_block_t
            implicit none
            type(var_block_t), pointer, intent(in) :: mb_vin(:)
            integer(kind_int)         , intent(in) :: nsin,nein
            external                               :: sub_via
            type(var_block_t), pointer, intent(in) :: mb_vout(:)
            integer(kind_int)         , intent(in) :: nsout,neout
            integer(kind_int)         , intent(in) :: ngst,nged
        end subroutine calc_bc_var_via_sub_sp
    end interface calc_bc_var_via_sub_sp    

    interface mb_var_pointer_create
        subroutine mb_var_pointer_create(mb_var,nst,ned)
            use mod_kndconsts, only : kind_int
            use mod_datatypes, only : var_block_t
            implicit none
            type(var_block_t), pointer, intent(inout) :: mb_var(:)
            integer(kind_int)         , intent(in)    :: nst,ned
        end subroutine mb_var_pointer_create
    end interface mb_var_pointer_create

    interface mb_var_pointer_delete
        subroutine mb_var_pointer_delete(mb_var)
            use mod_datatypes, only : var_block_t
            implicit none
            type(var_block_t), pointer, intent(inout) :: mb_var(:)
        end subroutine mb_var_pointer_delete
    end interface mb_var_pointer_delete

    interface mb_var_pointer_assign
        subroutine mb_var_pointer_assign(mb_vars,nvs_st,nvs_ed, &
                                         mb_vard,nvd_st,nvd_ed)
            use mod_kndconsts, only : kind_int
            use mod_datatypes, only : var_block_t
            implicit none
            type(var_block_t), pointer, intent(inout) :: mb_vars(:)
            integer(kind_int)         , intent(in)    :: nvs_st,nvs_ed
            type(var_block_t), pointer, intent(inout) :: mb_vard(:)
            integer(kind_int)         , intent(in)    :: nvd_st,nvd_ed
        end subroutine mb_var_pointer_assign
    end interface mb_var_pointer_assign  

    interface assign_mb_var_uniform
        subroutine assign_mb_var_uniform(mb_var,nst,ned,ngh,var)
            use mod_kndconsts, only : kind_int,kind_real
            use mod_datatypes, only : var_block_t
            implicit none
            type(var_block_t), pointer, intent(in) :: mb_var(:)
            integer(kind_int)         , intent(in) :: nst,ned,ngh
            real(kind_real)           , intent(in) :: var(nst:ned)
        end subroutine assign_mb_var_uniform
    end interface assign_mb_var_uniform
    
    interface assign_mb_var_uniform_sp
        subroutine assign_mb_var_uniform_sp(mb_var,nst,ned,ngh,var)
            use mod_kndconsts, only : kind_int,kind_real
            use mod_datatypes, only : var_block_t
            implicit none
            type(var_block_t), pointer, intent(in) :: mb_var(:)
            integer(kind_int)         , intent(in) :: nst,ned,ngh
            real(kind_real)           , intent(in) :: var(nst:ned)
        end subroutine assign_mb_var_uniform_sp
    end interface assign_mb_var_uniform_sp    

    interface assign_bc_var_nb
        subroutine assign_bc_var_nb(nb,mb_var,nst,ned,ngh,var)
            use mod_kndconsts, only : kind_int,kind_real
            use mod_datatypes, only : var_block_t
            implicit none
            integer(kind_int)         , intent(in) :: nb
            type(var_block_t), pointer, intent(in) :: mb_var(:)
            integer(kind_int)         , intent(in) :: nst,ned,ngh
            real(kind_real)           , intent(in) :: var(nst:ned)
        end subroutine assign_bc_var_nb
    end interface assign_bc_var_nb
    
    interface assign_bc_var_nb_sp
        subroutine assign_bc_var_nb_sp(nb,mb_var,nst,ned,ngh,var)
            use mod_kndconsts, only : kind_int,kind_real
            use mod_datatypes, only : var_block_t
            implicit none
            integer(kind_int)         , intent(in) :: nb
            type(var_block_t), pointer, intent(in) :: mb_var(:)
            integer(kind_int)         , intent(in) :: nst,ned,ngh
            real(kind_real)           , intent(in) :: var(nst:ned)
        end subroutine assign_bc_var_nb_sp
    end interface assign_bc_var_nb_sp    

    interface assign_bc_var_uniform
        subroutine assign_bc_var_uniform(mb_var,nst,ned,ngh,var)
            use mod_kndconsts, only : kind_int,kind_real
            use mod_datatypes, only : var_block_t
            implicit none
            type(var_block_t), pointer, intent(in) :: mb_var(:)
            integer(kind_int)         , intent(in) :: nst,ned,ngh
            real(kind_real)           , intent(in) :: var(nst:ned)
        end subroutine assign_bc_var_uniform
    end interface assign_bc_var_uniform
    
    interface assign_bc_var_uniform_sp
        subroutine assign_bc_var_uniform_sp(mb_var,nst,ned,ngh,var)
            use mod_kndconsts, only : kind_int,kind_real
            use mod_datatypes, only : var_block_t
            implicit none
            type(var_block_t), pointer, intent(in) :: mb_var(:)
            integer(kind_int)         , intent(in) :: nst,ned,ngh
            real(kind_real)           , intent(in) :: var(nst:ned)
        end subroutine assign_bc_var_uniform_sp
    end interface assign_bc_var_uniform_sp    

    interface ghost_bc_var_nb
        subroutine ghost_bc_var_nb(nb,mb_var,nst,ned)
            use mod_kndconsts, only : kind_int
            use mod_datatypes, only : var_block_t
            implicit none
            integer(kind_int)         , intent(in) :: nb
            type(var_block_t), pointer, intent(in) :: mb_var(:)
            integer(kind_int)         , intent(in) :: nst,ned
        end subroutine ghost_bc_var_nb
    end interface ghost_bc_var_nb
    
    interface ghost_bc_var_nb_cc
        subroutine ghost_bc_var_nb_cc(nb,mb_var,nst,ned)
            use mod_kndconsts, only : kind_int
            use mod_datatypes, only : var_block_t
            implicit none
            integer(kind_int)         , intent(in) :: nb
            type(var_block_t), pointer, intent(in) :: mb_var(:)
            integer(kind_int)         , intent(in) :: nst,ned
        end subroutine ghost_bc_var_nb_cc
    end interface ghost_bc_var_nb_cc
    
    interface ghost_bc_var_nb_sp
        subroutine ghost_bc_var_nb_sp(nb,mb_var,nst,ned)
            use mod_kndconsts, only : kind_int
            use mod_datatypes, only : var_block_t
            implicit none
            integer(kind_int)         , intent(in) :: nb
            type(var_block_t), pointer, intent(in) :: mb_var(:)
            integer(kind_int)         , intent(in) :: nst,ned
        end subroutine ghost_bc_var_nb_sp
    end interface ghost_bc_var_nb_sp    

    interface fill_corner_var_nb
        subroutine fill_corner_var_nb(nb,mb_var,nst,ned)
            use mod_kndconsts, only : kind_int
            use mod_datatypes, only : var_block_t
            implicit none
            integer(kind_int)         , intent(in) :: nb
            type(var_block_t), pointer, intent(in) :: mb_var(:)
            integer(kind_int)         , intent(in) :: nst,ned
        end subroutine fill_corner_var_nb
    end interface fill_corner_var_nb
    
    interface fill_corner_var_nb_cc
        subroutine fill_corner_var_nb_cc(nb,mb_var,nst,ned)
            use mod_kndconsts, only : kind_int
            use mod_datatypes, only : var_block_t
            implicit none
            integer(kind_int)         , intent(in) :: nb
            type(var_block_t), pointer, intent(in) :: mb_var(:)
            integer(kind_int)         , intent(in) :: nst,ned
        end subroutine fill_corner_var_nb_cc
    end interface fill_corner_var_nb_cc  
    
    interface fill_corner_var_nb_sp
        subroutine fill_corner_var_nb_sp(nb,mb_var,nst,ned)
            use mod_kndconsts, only : kind_int
            use mod_datatypes, only : var_block_t
            implicit none
            integer(kind_int)         , intent(in) :: nb
            type(var_block_t), pointer, intent(in) :: mb_var(:)
            integer(kind_int)         , intent(in) :: nst,ned
        end subroutine fill_corner_var_nb_sp
    end interface fill_corner_var_nb_sp    

    interface ghost_bc_var_all
        subroutine ghost_bc_var_all(mb_var,nst,ned)
            use mod_kndconsts, only : kind_int
            use mod_datatypes, only : var_block_t
            implicit none
            type(var_block_t), pointer, intent(in) :: mb_var(:)
            integer(kind_int)         , intent(in) :: nst,ned
        end subroutine ghost_bc_var_all
    end interface ghost_bc_var_all
    
    interface ghost_bc_var_all_cc
        subroutine ghost_bc_var_all_cc(mb_var,nst,ned)
            use mod_kndconsts, only : kind_int
            use mod_datatypes, only : var_block_t
            implicit none
            type(var_block_t), pointer, intent(in) :: mb_var(:)
            integer(kind_int)         , intent(in) :: nst,ned
        end subroutine ghost_bc_var_all_cc
    end interface ghost_bc_var_all_cc
    
    interface ghost_bc_var_all_sp
        subroutine ghost_bc_var_all_sp(mb_var,nst,ned)
            use mod_kndconsts, only : kind_int
            use mod_datatypes, only : var_block_t
            implicit none
            type(var_block_t), pointer, intent(in) :: mb_var(:)
            integer(kind_int)         , intent(in) :: nst,ned
        end subroutine ghost_bc_var_all_sp
    end interface ghost_bc_var_all_sp    

    interface assign_com_var_nb
        subroutine assign_com_var_nb(nb,mb_var,nst,ned,var)
            use mod_kndconsts, only : kind_int,kind_real
            use mod_datatypes, only : var_block_t
            use mod_fieldvars, only : mb_top
            implicit none
            integer(kind_int)         , intent(in) :: nb
            type(var_block_t), pointer, intent(in) :: mb_var(:)
            integer(kind_int)         , intent(in) :: nst,ned
            real(kind_real)           , intent(in) :: var(nst:ned)
        end subroutine assign_com_var_nb
    end interface assign_com_var_nb

    interface assign_var_via_com_nb
        subroutine assign_var_via_com_nb(nb,mb_var,nst,ned)
            use mod_kndconsts, only : kind_int
            use mod_datatypes, only : var_block_t
            use mod_fieldvars, only : mb_top
            implicit none
            integer(kind_int)         , intent(in) :: nb
            type(var_block_t), pointer, intent(in) :: mb_var(:)
            integer(kind_int)         , intent(in) :: nst,ned
        end subroutine assign_var_via_com_nb
    end interface assign_var_via_com_nb
    
    interface assign_var_via_com_nb_sp
        subroutine assign_var_via_com_nb_sp(nb,mb_var,nst,ned)
            use mod_kndconsts, only : kind_int
            use mod_datatypes, only : var_block_t
            use mod_fieldvars, only : mb_top
            implicit none
            integer(kind_int)         , intent(in) :: nb
            type(var_block_t), pointer, intent(in) :: mb_var(:)
            integer(kind_int)         , intent(in) :: nst,ned
        end subroutine assign_var_via_com_nb_sp
    end interface assign_var_via_com_nb_sp    

    interface exchange_singulars
        subroutine exchange_singulars(mb_var,nst,ned,nswmem,nswave)
            use mod_kndconsts, only : kind_int
            use mod_datatypes, only : var_block_t
            implicit none
            type(var_block_t), pointer, intent(in) :: mb_var(:)
            integer(kind_int),          intent(in) :: nst,ned
            integer(kind_int),          intent(in) :: nswmem
            integer(kind_int),          intent(in) :: nswave
        end subroutine exchange_singulars
    end interface exchange_singulars
    
    interface exchange_singulars_cc
        subroutine exchange_singulars_cc(mb_var,nst,ned,nswmem,nswave)
            use mod_kndconsts, only : kind_int
            use mod_datatypes, only : var_block_t
            implicit none
            type(var_block_t), pointer, intent(in) :: mb_var(:)
            integer(kind_int),          intent(in) :: nst,ned
            integer(kind_int),          intent(in) :: nswmem
            integer(kind_int),          intent(in) :: nswave
        end subroutine exchange_singulars_cc
    end interface exchange_singulars_cc    

    interface average_singulars
        subroutine average_singulars(mb_var,nst,ned,nswmem,nswave)
            use mod_kndconsts, only : kind_int
            use mod_datatypes, only : var_block_t
            implicit none
            type(var_block_t), pointer, intent(in) :: mb_var(:)
            integer(kind_int),          intent(in) :: nst,ned
            integer(kind_int),          intent(in) :: nswmem
            integer(kind_int),          intent(in) :: nswave
        end subroutine average_singulars
    end interface average_singulars
    
    interface average_singulars_cc
        subroutine average_singulars_cc(mb_var,nst,ned,nswmem,nswave)
            use mod_kndconsts, only : kind_int
            use mod_datatypes, only : var_block_t
            implicit none
            type(var_block_t), pointer, intent(in) :: mb_var(:)
            integer(kind_int),          intent(in) :: nst,ned
            integer(kind_int),          intent(in) :: nswmem
            integer(kind_int),          intent(in) :: nswave
        end subroutine average_singulars_cc
    end interface average_singulars_cc

    interface calc_inviscd
        subroutine calc_inviscd(nb,ngn,nge, &
                                npvs,mb_pv,neqn,mb_rhs, &
                                sub_limit,sub_intnon,sub_intplt, &
                                sub_fvs,sub_flux,sub_scheme, &
                                chkid,chklim)
            use mod_kndconsts, only : kind_int,kind_real
            use mod_datatypes, only : var_block_t
            implicit none
            integer(kind_int), intent(in) :: nb
            integer(kind_int), intent(in) :: ngn,nge
            integer(kind_int), intent(in) :: npvs
            type(var_block_t),    pointer :: mb_pv(:)
            integer(kind_int), intent(in) :: neqn
            type(var_block_t),    pointer :: mb_rhs(:)
            real(kind_real),     external :: sub_limit
            external                      :: sub_intnon,sub_intplt
            external                      :: sub_fvs,sub_flux
            external                      :: sub_scheme
            integer(kind_int), intent(in) :: chkid(npvs)
            real(kind_real),   intent(in) :: chklim(npvs)
        end subroutine calc_inviscd
    end interface calc_inviscd
    
    interface calc_inviscd_sp
        subroutine calc_inviscd_sp(nb,ngn,nge, &
                                   npvs,mb_pv,neqn,mb_rhs, &
                                   sub_limit,sub_intnon,sub_intplt, &
                                   sub_fvs,sub_flux,sub_scheme, &
                                   chkid,chklim)
            use mod_kndconsts, only : kind_int,kind_real
            use mod_datatypes, only : var_block_t
            implicit none
            integer(kind_int), intent(in) :: nb
            integer(kind_int), intent(in) :: ngn,nge
            integer(kind_int), intent(in) :: npvs
            type(var_block_t),    pointer :: mb_pv(:)
            integer(kind_int), intent(in) :: neqn
            type(var_block_t),    pointer :: mb_rhs(:)
            real(kind_real),     external :: sub_limit
            external                      :: sub_intnon,sub_intplt
            external                      :: sub_fvs,sub_flux
            external                      :: sub_scheme
            integer(kind_int), intent(in) :: chkid(npvs)
            real(kind_real),   intent(in) :: chklim(npvs)
        end subroutine calc_inviscd_sp
    end interface calc_inviscd_sp    

    interface cal_sptocp
        subroutine cal_sptocp(nb,ngn,nge,npvs,mb_pv,neqn,sub_intnon)
            use mod_kndconsts, only : kind_int,kind_real
            use mod_constants, only : nfsf_con_d1nint,nsw_dir_close
            use mod_constants, only : nbc_inter_scheme,nbc_intbc_scheme
            use mod_datatypes, only : var_block_t,fld_array_t
            use mod_variables, only : nsw_kdir
            use mod_fieldvars, only : mb_topsp,mb_fsfsp,mb_pvfp,mb_sxyzsp,mb_vsl,mb_vst
            implicit none
            integer(kind_int), intent(in) :: nb
            integer(kind_int), intent(in) :: ngn,nge
            integer(kind_int), intent(in) :: npvs
            type(var_block_t),    pointer :: mb_pv(:)
            integer(kind_int), intent(in) :: neqn
            external                      :: sub_intnon
            real(kind_real), external  :: nolimiter
            integer(kind_int)          :: i,j,k,m,nfs,nfe,ierr
            integer(kind_int)          :: ni,nj,nk,stn,edn,ste,ede
            integer(kind_int)          :: stn0,edn0,ste0,ede0
            integer(kind_int)          :: nkst,nked,njed,ndec
            integer(kind_int)          :: chkid(npvs)
            real(kind_real)            :: chklim(npvs)
            real(kind_real)            :: vl(1:npvs),vr(1:npvs)
            type(fld_array_t), pointer :: sxyzsp(:),pv(:),pvfp(:),vsl(:),vst(:)
            type(fld_array_t), pointer :: fsfs(:),fsfe(:)
            real(kind_real),   pointer :: pvn(:,:),pvl(:,:),pvr(:,:),sn(:,:)
        end subroutine cal_sptocp
    end interface cal_sptocp
    
    interface calc_wall_dist
        subroutine calc_wall_dist(mb_dst,bc_target,mark)
            use mod_kndconsts, only : kind_int
            use mod_datatypes, only : var_block_t
            implicit none
            type(var_block_t), pointer, intent(in) :: mb_dst(:)
            integer(kind_int),          intent(in) :: bc_target,mark
        end subroutine calc_wall_dist
    end interface calc_wall_dist

end module mod_interface
