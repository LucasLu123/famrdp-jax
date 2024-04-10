
subroutine rhs_invscd
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nlimit_none,nlimit_minmod
    use mod_constants, only : nlimit_vanleer,nlimit_vanalbada
    use mod_variables, only : nlimit
    use mod_fieldvars, only : nblkcoms,blkcoms
    implicit none
    integer(kind_int)         :: nc,nb,ierr
    real(kind_real), external :: nolimiter,minmod
    real(kind_real), external :: vanleer,vanalbada

    do nc=1,nblkcoms
        nb = blkcoms(nc)%nb
        select case(nlimit)
        case(nlimit_none)
            call invscd_intnon(nb,nolimiter)
        case(nlimit_minmod)
            call invscd_intnon(nb,minmod   )
        case(nlimit_vanleer)
            call invscd_intnon(nb,vanleer  )
        case(nlimit_vanalbada)
            call invscd_intnon(nb,vanalbada)
        end select
    end do

end subroutine rhs_invscd

subroutine rhs_invscd_sp
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nlimit_none,nlimit_minmod
    use mod_constants, only : nlimit_vanleer,nlimit_vanalbada
    use mod_variables, only : nlimit
    use mod_fieldvars, only : nblkcoms,blkcomssp
    implicit none
    integer(kind_int)         :: nc,nb,ierr
    real(kind_real), external :: nolimiter,minmod
    real(kind_real), external :: vanleer,vanalbada

    do nc=1,nblkcoms
        nb = blkcomssp(nc)%nb
        select case(nlimit)
        case(nlimit_none)
            call invscd_intnon_sp(nb,nolimiter)
        case(nlimit_minmod)
            call invscd_intnon_sp(nb,minmod   )
        case(nlimit_vanleer)
            call invscd_intnon_sp(nb,vanleer  )
        case(nlimit_vanalbada)
            call invscd_intnon_sp(nb,vanalbada)
        end select
    end do

end subroutine rhs_invscd_sp

subroutine invscd_intnon(nb,sub_limit)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nintnon_muscl2pv
    use mod_constants, only : nintnon_wcns5pv,nintnon_wcns5cv
    use mod_constants, only : nintnon_wcns7pv,nintnon_wcns7cv
    use mod_constants, only : nintnon_hdcs5ei,nintnon_hdcs7ci
    use mod_constants, only : nintnon_hdcs5ci,nintnon_scsl3ci
    use mod_constants, only : nintnon_scsl5ci
    use mod_constants, only : nintnon_scsh3ci,nintnon_scsh5ci
    use mod_constants, only : nintnon_scsn2ci,nintnon_scsn3ci
    use mod_constants, only : nintnon_scsn4ci,nintnon_scsh3pi
    use mod_constants, only : nintnon_scsh5pi,nintnon_scsn2pi
    use mod_constants, only : nintnon_scsn3pi,nintnon_scsn4pi
    use mod_constants, only : nintnon_dcsh5ci,nintnon_dcsh7ci
    use mod_constants, only : nintnon_dcsh5pi,nintnon_dcsh7pi
    use mod_variables, only : nintnon
    implicit none
    integer(kind_int), intent(in) :: nb
    real(kind_real), external     :: sub_limit
    external :: muscl2pv,wcns5pv,wcns5cv
    external :: wcns7pv,wcns7cv
    external :: hdcs5ei,hdcs7ci
    external :: hdcs5ci,scsl3ci
    external :: scsl5ci
    external :: scsh3ci,scsh5ci
    external :: scsn2ci,scsn3ci
    external :: scsn4ci,scsh3pi
    external :: scsh5pi,scsn2pi
    external :: scsn3pi,scsn4pi
    external :: dcsh5ci,dcsh7ci
    external :: dcsh5pi,dcsh7pi


    select case(nintnon)
    case(nintnon_muscl2pv)
       call invscd_intplt(nb,sub_limit,muscl2pv)
    case(nintnon_wcns5pv)
       call invscd_intplt(nb,sub_limit,wcns5pv)
    case(nintnon_wcns5cv)
       call invscd_intplt(nb,sub_limit,wcns5cv)
    case(nintnon_wcns7pv)
       call invscd_intplt(nb,sub_limit,wcns7pv)
    case(nintnon_wcns7cv)
       call invscd_intplt(nb,sub_limit,wcns7cv)
    case(nintnon_hdcs5ei)
       call invscd_intplt(nb,sub_limit,hdcs5ei)
    case(nintnon_hdcs7ci)
       call invscd_intplt(nb,sub_limit,hdcs7ci)
    case(nintnon_hdcs5ci)
       call invscd_intplt(nb,sub_limit,hdcs5ci)
    case(nintnon_scsl3ci)
       call invscd_intplt(nb,sub_limit,scsl3ci)
    case(nintnon_scsl5ci)
       call invscd_intplt(nb,sub_limit,scsl5ci)
    case(nintnon_scsh3ci)
       call invscd_intplt(nb,sub_limit,scsh3ci)
    case(nintnon_scsh5ci)
       call invscd_intplt(nb,sub_limit,scsh5ci)
    case(nintnon_scsn2ci)
       call invscd_intplt(nb,sub_limit,scsn2ci)       
    case(nintnon_scsn3ci)
       call invscd_intplt(nb,sub_limit,scsn3ci)
    case(nintnon_scsn4ci)
       call invscd_intplt(nb,sub_limit,scsn4ci) 
    case(nintnon_scsh3pi)
       call invscd_intplt(nb,sub_limit,scsh3pi)
    case(nintnon_scsh5pi)
       call invscd_intplt(nb,sub_limit,scsh5pi)
    case(nintnon_scsn2pi)
       call invscd_intplt(nb,sub_limit,scsn2pi)       
    case(nintnon_scsn3pi)
       call invscd_intplt(nb,sub_limit,scsn3pi)
    case(nintnon_scsn4pi)
       call invscd_intplt(nb,sub_limit,scsn4pi)
    case(nintnon_dcsh5ci)
       call invscd_intplt(nb,sub_limit,dcsh5ci)
    case(nintnon_dcsh7ci)
       call invscd_intplt(nb,sub_limit,dcsh7ci)
    case(nintnon_dcsh5pi)
       call invscd_intplt(nb,sub_limit,dcsh5pi)
    case(nintnon_dcsh7pi)
       call invscd_intplt(nb,sub_limit,dcsh7pi)
    case default

    end select

end subroutine invscd_intnon

subroutine invscd_intnon_sp(nb,sub_limit)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nintnon_muscl2pv
    use mod_constants, only : nintnon_wcns5pv,nintnon_wcns5cv
    use mod_constants, only : nintnon_wcns7pv,nintnon_wcns7cv
    use mod_constants, only : nintnon_hdcs5ei,nintnon_hdcs7ci
    use mod_constants, only : nintnon_hdcs5ci,nintnon_scsl3ci
    use mod_constants, only : nintnon_scsl5ci
    use mod_constants, only : nintnon_scsh3ci,nintnon_scsh5ci
    use mod_constants, only : nintnon_scsn2ci,nintnon_scsn3ci
    use mod_constants, only : nintnon_scsn4ci,nintnon_scsh3pi
    use mod_constants, only : nintnon_scsh5pi,nintnon_scsn2pi
    use mod_constants, only : nintnon_scsn3pi,nintnon_scsn4pi
    use mod_constants, only : nintnon_dcsh5ci,nintnon_dcsh7ci
    use mod_constants, only : nintnon_dcsh5pi,nintnon_dcsh7pi
    use mod_variables, only : nintnon
    implicit none
    integer(kind_int), intent(in) :: nb
    real(kind_real), external     :: sub_limit
    external :: muscl2pvcc,wcns5pvcc,wcns5cvcc
    external :: wcns7pvcc,wcns7cvcc
    external :: hdcs5eicc,hdcs7cicc
    external :: hdcs5cicc,scsl3cicc
    external :: scsl5cicc
    external :: scsh3cicc,scsh5cicc
    external :: scsn2cicc,scsn3cicc
    external :: scsn4cicc,scsh3picc
    external :: scsh5picc,scsn2picc
    external :: scsn3picc,scsn4picc
    external :: dcsh5cicc,dcsh7cicc
    external :: dcsh5picc,dcsh7picc


    select case(nintnon)
    case(nintnon_muscl2pv)
       call invscd_intplt_sp(nb,sub_limit,muscl2pvcc)
    case(nintnon_wcns5pv)
       call invscd_intplt_sp(nb,sub_limit,wcns5pvcc)
    case(nintnon_wcns5cv)
       call invscd_intplt_sp(nb,sub_limit,wcns5cvcc)
    case(nintnon_wcns7pv)
       call invscd_intplt_sp(nb,sub_limit,wcns7pvcc)
    case(nintnon_wcns7cv)
       call invscd_intplt_sp(nb,sub_limit,wcns7cvcc)
    case(nintnon_hdcs5ei)
       call invscd_intplt_sp(nb,sub_limit,hdcs5eicc)
    case(nintnon_hdcs7ci)
       call invscd_intplt_sp(nb,sub_limit,hdcs7cicc)
    case(nintnon_hdcs5ci)
       call invscd_intplt_sp(nb,sub_limit,hdcs5cicc)
    case(nintnon_scsl3ci)
       call invscd_intplt_sp(nb,sub_limit,scsl3cicc)
    case(nintnon_scsl5ci)
       call invscd_intplt_sp(nb,sub_limit,scsl5cicc)
    case(nintnon_scsh3ci)
       call invscd_intplt_sp(nb,sub_limit,scsh3cicc)
    case(nintnon_scsh5ci)
       call invscd_intplt_sp(nb,sub_limit,scsh5cicc)
    case(nintnon_scsn2ci)
       call invscd_intplt_sp(nb,sub_limit,scsn2cicc)       
    case(nintnon_scsn3ci)
       call invscd_intplt_sp(nb,sub_limit,scsn3cicc)
    case(nintnon_scsn4ci)
       call invscd_intplt_sp(nb,sub_limit,scsn4cicc) 
    case(nintnon_scsh3pi)
       call invscd_intplt_sp(nb,sub_limit,scsh3picc)
    case(nintnon_scsh5pi)
       call invscd_intplt_sp(nb,sub_limit,scsh5picc)
    case(nintnon_scsn2pi)
       call invscd_intplt_sp(nb,sub_limit,scsn2picc)       
    case(nintnon_scsn3pi)
       call invscd_intplt_sp(nb,sub_limit,scsn3picc)
    case(nintnon_scsn4pi)
       call invscd_intplt_sp(nb,sub_limit,scsn4picc)
    case(nintnon_dcsh5ci)
       call invscd_intplt_sp(nb,sub_limit,dcsh5cicc)
    case(nintnon_dcsh7ci)
       call invscd_intplt_sp(nb,sub_limit,dcsh7cicc)
    case(nintnon_dcsh5pi)
       call invscd_intplt_sp(nb,sub_limit,dcsh5picc)
    case(nintnon_dcsh7pi)
       call invscd_intplt_sp(nb,sub_limit,dcsh7picc)
    case default

    end select

end subroutine invscd_intnon_sp

subroutine invscd_intplt(nb,sub_limit,sub_intnon)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nintplt_node2
    use mod_constants, only : nintplt_node4,nintplt_node6
    use mod_constants, only : nintplt_node6w,nintplt_node8
    use mod_constants, only : nintplt_scsl4,nintplt_scsl6
    use mod_constants, only : nintplt_node6e,nintplt_node8e
    use mod_constants, only : nintplt_scsl4e,nintplt_scsl6e
    use mod_variables, only : nd1int_con
    implicit none
    integer(kind_int), intent(in) :: nb
    real(kind_real), external     :: sub_limit
    external                      :: sub_intnon
    external :: ve_via_node2,ve_via_node4
    external :: ve_via_node6,ve_via_node8
    external :: ve_via_scsl4,ve_via_scsl6
    external :: ve_via_node6e,ve_via_node8e
    external :: ve_via_scsl4e,ve_via_scsl6e

    select case(nd1int_con)
    case(nintplt_node2)
        call invscd_flux(nb,sub_limit,sub_intnon,ve_via_node2)
    case(nintplt_node4)
        call invscd_flux(nb,sub_limit,sub_intnon,ve_via_node4)
    case(nintplt_node6,nintplt_node6w)
        call invscd_flux(nb,sub_limit,sub_intnon,ve_via_node6)
    case(nintplt_node6e)
        call invscd_flux(nb,sub_limit,sub_intnon,ve_via_node6e)
    case(nintplt_node8)
        call invscd_flux(nb,sub_limit,sub_intnon,ve_via_node8)
    case(nintplt_node8e)
        call invscd_flux(nb,sub_limit,sub_intnon,ve_via_node8e)
    case(nintplt_scsl4)
        call invscd_flux(nb,sub_limit,sub_intnon,ve_via_scsl4)
    case(nintplt_scsl4e)
        call invscd_flux(nb,sub_limit,sub_intnon,ve_via_scsl4e)
    case(nintplt_scsl6)
        call invscd_flux(nb,sub_limit,sub_intnon,ve_via_scsl6)
    case(nintplt_scsl6e)
        call invscd_flux(nb,sub_limit,sub_intnon,ve_via_scsl6e)
    case default

    end select

end subroutine invscd_intplt

subroutine invscd_intplt_sp(nb,sub_limit,sub_intnon)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nintplt_node2
    use mod_constants, only : nintplt_node4,nintplt_node6
    use mod_constants, only : nintplt_node6w,nintplt_node8
    use mod_constants, only : nintplt_scsl4,nintplt_scsl6
    use mod_constants, only : nintplt_node6e,nintplt_node8e
    use mod_constants, only : nintplt_scsl4e,nintplt_scsl6e
    use mod_variables, only : nd1int_con
    implicit none
    integer(kind_int), intent(in) :: nb
    real(kind_real), external     :: sub_limit
    external                      :: sub_intnon
    external :: ve_via_node2cc,ve_via_node4cc
    external :: ve_via_node6cc,ve_via_node8cc
    external :: ve_via_scsl4cc,ve_via_scsl6cc
    external :: ve_via_node6ecc,ve_via_node8ecc
    external :: ve_via_scsl4ecc,ve_via_scsl6ecc

    select case(nd1int_con)
    case(nintplt_node2)
        call invscd_flux_sp(nb,sub_limit,sub_intnon,ve_via_node2cc)
    case(nintplt_node4)
        call invscd_flux_sp(nb,sub_limit,sub_intnon,ve_via_node4cc)
    case(nintplt_node6,nintplt_node6w)
        call invscd_flux_sp(nb,sub_limit,sub_intnon,ve_via_node6cc)
    case(nintplt_node6e)
        call invscd_flux_sp(nb,sub_limit,sub_intnon,ve_via_node6ecc)
    case(nintplt_node8)
        call invscd_flux_sp(nb,sub_limit,sub_intnon,ve_via_node8cc)
    case(nintplt_node8e)
        call invscd_flux_sp(nb,sub_limit,sub_intnon,ve_via_node8ecc)
    case(nintplt_scsl4)
        call invscd_flux_sp(nb,sub_limit,sub_intnon,ve_via_scsl4cc)
    case(nintplt_scsl4e)
        call invscd_flux_sp(nb,sub_limit,sub_intnon,ve_via_scsl4ecc)
    case(nintplt_scsl6)
        call invscd_flux_sp(nb,sub_limit,sub_intnon,ve_via_scsl6cc)
    case(nintplt_scsl6e)
        call invscd_flux_sp(nb,sub_limit,sub_intnon,ve_via_scsl6ecc)
    case default

    end select

end subroutine invscd_intplt_sp

subroutine invscd_flux(nb,sub_limit,sub_intnon,sub_intplt)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nflux_steger,nflux_sw_mod
    use mod_constants, only : nflux_vanleer,nflux_roe,nflux_roe_prec,nflux_slau,nflux_scmp
    use mod_variables, only : nflux
    implicit none
    integer(kind_int), intent(in) :: nb
    real(kind_real), external     :: sub_limit
    external                      :: sub_intnon,sub_intplt
    external :: flux_steger,flux_sw_mod,flux_vanleer,flux_roe,flux_roe_prec,flux_slau,flux_roe_scmp

    !todo SCM 新增一个通量
    select case(nflux)
    case(nflux_steger)
        call invscd_scheme(nb,sub_limit,sub_intnon,sub_intplt,flux_steger)
    case(nflux_sw_mod)
        call invscd_scheme(nb,sub_limit,sub_intnon,sub_intplt,flux_sw_mod)
    case(nflux_vanleer)
        call invscd_scheme(nb,sub_limit,sub_intnon,sub_intplt,flux_vanleer)
    case(nflux_roe)
        call invscd_scheme(nb,sub_limit,sub_intnon,sub_intplt,flux_roe)
    case(nflux_roe_prec)
        call invscd_scheme(nb,sub_limit,sub_intnon,sub_intplt,flux_roe_prec)
    case(nflux_slau)
        call invscd_scheme(nb,sub_limit,sub_intnon,sub_intplt,flux_slau)
    case(nflux_scmp)
        call invscd_scheme(nb,sub_limit,sub_intnon,sub_intplt,flux_roe_scmp)
    case default
    end select

end subroutine invscd_flux

subroutine invscd_flux_sp(nb,sub_limit,sub_intnon,sub_intplt)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nflux_steger,nflux_sw_mod
    use mod_constants, only : nflux_vanleer,nflux_roe,nflux_roe_prec,nflux_slau,nflux_scmp
    use mod_variables, only : nflux
    implicit none
    integer(kind_int), intent(in) :: nb
    real(kind_real), external     :: sub_limit
    external                      :: sub_intnon,sub_intplt
    external :: flux_steger,flux_sw_mod,flux_vanleer,flux_roe,flux_roe_prec,flux_slau,flux_roe_scmp

    !todo SCM 新增一个通量
    select case(nflux)
    case(nflux_steger)
        call invscd_scheme_sp(nb,sub_limit,sub_intnon,sub_intplt,flux_steger)
    case(nflux_sw_mod)
        call invscd_scheme_sp(nb,sub_limit,sub_intnon,sub_intplt,flux_sw_mod)
    case(nflux_vanleer)
        call invscd_scheme_sp(nb,sub_limit,sub_intnon,sub_intplt,flux_vanleer)
    case(nflux_roe)
        call invscd_scheme_sp(nb,sub_limit,sub_intnon,sub_intplt,flux_roe)
    case(nflux_roe_prec)
        call invscd_scheme_sp(nb,sub_limit,sub_intnon,sub_intplt,flux_roe_prec)
    case(nflux_slau)
        call invscd_scheme_sp(nb,sub_limit,sub_intnon,sub_intplt,flux_slau)
    case(nflux_scmp)
        call invscd_scheme_sp(nb,sub_limit,sub_intnon,sub_intplt,flux_roe_scmp)
    case default
    end select

end subroutine invscd_flux_sp

subroutine invscd_scheme(nb,sub_limit,sub_intnon,sub_intplt,sub_fvs)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : id_ro,id_ps,nderive_edge2
    use mod_constants, only : nderive_edge4,nderive_edge6
    use mod_constants, only : nderive_ehen4,nderive_ehen6
    use mod_constants, only : nderive_ehen8,nderive_ehcs6
    use mod_constants, only : nderive_scsl4,nderive_scsl6
    use mod_constants, only : nderive_ehen6e,nderive_ehcs6e
    use mod_constants, only : nderive_ehen8e,nscmp_non
    use mod_constants, only : nderive_scsl4e,nderive_scsl6e
    use mod_variables, only : rlimit,plimit,nscmp
    use mod_variables, only : nd1der_con,nghnode,nghedge
    use mod_fieldvars, only : npvs,mb_pv,neqn,mb_rhs
    use mod_interface, only : calc_inviscd
    implicit none
    integer(kind_int), intent(in) :: nb
    real(kind_real), external     :: sub_limit
    external                      :: sub_intnon,sub_intplt
    external                      :: sub_fvs
    integer(kind_int) :: chkid(npvs)
    real(kind_real)   :: chklim(npvs)
    external          :: flux_euler
    external          :: dn_via_edge6,dn_via_edge4,dn_via_edge2
    external          :: dn_via_ehen4,dn_via_ehen6,dn_via_ehen8,dn_via_ehcs6
    external          :: dn_via_scsl4,dn_via_scsl6
    external          :: dn_via_ehen6e,dn_via_ehen8e,dn_via_ehcs6e
    external          :: dn_via_scsl4e,dn_via_scsl6e

    if (nscmp > nscmp_non) then
        chkid(:) = 0
        chkid(id_ro) = 1

        chklim(id_ro) = rlimit(1)
    else
        chkid(:) = 0
        chkid(id_ro) = 1
        chkid(id_ps) = 1

        chklim(id_ro) = rlimit(1)
        chklim(id_ps) = plimit(1)        
    end if

    select case(nd1der_con)
    case(nderive_edge2)
        call calc_inviscd(nb,nghnode,nghedge, &
                          npvs,mb_pv,neqn,mb_rhs, &
                          sub_limit,sub_intnon,sub_intplt, &
                          sub_fvs,flux_euler,dn_via_edge2, &
                          chkid,chklim)
    case(nderive_edge4)
        call calc_inviscd(nb,nghnode,nghedge, &
                          npvs,mb_pv,neqn,mb_rhs, &
                          sub_limit,sub_intnon,sub_intplt, &
                          sub_fvs,flux_euler,dn_via_edge4, &
                          chkid,chklim)
    case(nderive_edge6)
        call calc_inviscd(nb,nghnode,nghedge, &
                          npvs,mb_pv,neqn,mb_rhs, &
                          sub_limit,sub_intnon,sub_intplt, &
                          sub_fvs,flux_euler,dn_via_edge6, &
                          chkid,chklim)
    case(nderive_ehen4)
        call calc_inviscd(nb,nghnode,nghedge, &
                          npvs,mb_pv,neqn,mb_rhs, &
                          sub_limit,sub_intnon,sub_intplt, &
                          sub_fvs,flux_euler,dn_via_ehen4, &
                          chkid,chklim)
    case(nderive_ehen6)
        call calc_inviscd(nb,nghnode,nghedge, &
                          npvs,mb_pv,neqn,mb_rhs, &
                          sub_limit,sub_intnon,sub_intplt, &
                          sub_fvs,flux_euler,dn_via_ehen6, &
                          chkid,chklim)
    case(nderive_ehen6e)
        call calc_inviscd(nb,nghnode,nghedge, &
                          npvs,mb_pv,neqn,mb_rhs, &
                          sub_limit,sub_intnon,sub_intplt, &
                          sub_fvs,flux_euler,dn_via_ehen6e, &
                          chkid,chklim)
    case(nderive_ehcs6)
        call calc_inviscd(nb,nghnode,nghedge, &
                          npvs,mb_pv,neqn,mb_rhs, &
                          sub_limit,sub_intnon,sub_intplt, &
                          sub_fvs,flux_euler,dn_via_ehcs6, &
                          chkid,chklim)
    case(nderive_ehcs6e)
        call calc_inviscd(nb,nghnode,nghedge, &
                          npvs,mb_pv,neqn,mb_rhs, &
                          sub_limit,sub_intnon,sub_intplt, &
                          sub_fvs,flux_euler,dn_via_ehcs6e, &
                          chkid,chklim)
    case(nderive_ehen8)
        call calc_inviscd(nb,nghnode,nghedge, &
                          npvs,mb_pv,neqn,mb_rhs, &
                          sub_limit,sub_intnon,sub_intplt, &
                          sub_fvs,flux_euler,dn_via_ehen8, &
                          chkid,chklim)
    case(nderive_ehen8e)
        call calc_inviscd(nb,nghnode,nghedge, &
                          npvs,mb_pv,neqn,mb_rhs, &
                          sub_limit,sub_intnon,sub_intplt, &
                          sub_fvs,flux_euler,dn_via_ehen8e, &
                          chkid,chklim)
    case(nderive_scsl4)
        call calc_inviscd(nb,nghnode,nghedge, &
                          npvs,mb_pv,neqn,mb_rhs, &
                          sub_limit,sub_intnon,sub_intplt, &
                          sub_fvs,flux_euler,dn_via_scsl4, &
                          chkid,chklim)
    case(nderive_scsl4e)
        call calc_inviscd(nb,nghnode,nghedge, &
                          npvs,mb_pv,neqn,mb_rhs, &
                          sub_limit,sub_intnon,sub_intplt, &
                          sub_fvs,flux_euler,dn_via_scsl4e, &
                          chkid,chklim)
    case(nderive_scsl6)
        call calc_inviscd(nb,nghnode,nghedge, &
                          npvs,mb_pv,neqn,mb_rhs, &
                          sub_limit,sub_intnon,sub_intplt, &
                          sub_fvs,flux_euler,dn_via_scsl6, &
                          chkid,chklim)
    case(nderive_scsl6e)
        call calc_inviscd(nb,nghnode,nghedge, &
                          npvs,mb_pv,neqn,mb_rhs, &
                          sub_limit,sub_intnon,sub_intplt,  &
                          sub_fvs,flux_euler,dn_via_scsl6e, &
                          chkid,chklim)
    case default
    end select

end subroutine invscd_scheme

subroutine invscd_scheme_sp(nb,sub_limit,sub_intnon,sub_intplt,sub_fvs)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : id_ro,id_ps,nderive_edge2
    use mod_constants, only : nderive_edge4,nderive_edge6
    use mod_constants, only : nderive_ehen4,nderive_ehen6
    use mod_constants, only : nderive_ehen8,nderive_ehcs6
    use mod_constants, only : nderive_scsl4,nderive_scsl6
    use mod_constants, only : nderive_ehen6e,nderive_ehcs6e
    use mod_constants, only : nderive_ehen8e,nscmp_non
    use mod_constants, only : nderive_scsl4e,nderive_scsl6e
    use mod_variables, only : rlimit,plimit,nscmp
    use mod_variables, only : nd1der_con,nghnode,nghedge
    use mod_fieldvars, only : npvs,mb_pv,neqn,mb_rhs
    use mod_interface, only : calc_inviscd_sp
    implicit none
    integer(kind_int), intent(in) :: nb
    real(kind_real), external     :: sub_limit
    external                      :: sub_intnon,sub_intplt
    external                      :: sub_fvs
    integer(kind_int) :: chkid(npvs)
    real(kind_real)   :: chklim(npvs)
    external          :: flux_euler
    external          :: dn_via_edge6cc,dn_via_edge4cc,dn_via_edge2cc
    external          :: dn_via_ehen4cc,dn_via_ehen6cc,dn_via_ehen8cc,dn_via_ehcs6cc
    external          :: dn_via_scsl4cc,dn_via_scsl6cc
    external          :: dn_via_ehen6ecc,dn_via_ehen8ecc,dn_via_ehcs6ecc
    external          :: dn_via_scsl4ecc,dn_via_scsl6ecc

    if (nscmp > nscmp_non) then
        chkid(:) = 0
        chkid(id_ro) = 1

        chklim(id_ro) = rlimit(1)
    else
        chkid(:) = 0
        chkid(id_ro) = 1
        chkid(id_ps) = 1

        chklim(id_ro) = rlimit(1)
        chklim(id_ps) = plimit(1)        
    end if

    select case(nd1der_con)
    case(nderive_edge2)
        call calc_inviscd_sp(nb,nghnode,nghedge, &
                             npvs,mb_pv,neqn,mb_rhs, &
                             sub_limit,sub_intnon,sub_intplt, &
                             sub_fvs,flux_euler,dn_via_edge2cc, &
                             chkid,chklim)
    case(nderive_edge4)
        call calc_inviscd_sp(nb,nghnode,nghedge, &
                             npvs,mb_pv,neqn,mb_rhs, &
                             sub_limit,sub_intnon,sub_intplt, &
                             sub_fvs,flux_euler,dn_via_edge4cc, &
                             chkid,chklim)
    case(nderive_edge6)
        call calc_inviscd_sp(nb,nghnode,nghedge, &
                             npvs,mb_pv,neqn,mb_rhs, &
                             sub_limit,sub_intnon,sub_intplt, &
                             sub_fvs,flux_euler,dn_via_edge6cc, &
                             chkid,chklim)
    case(nderive_ehen4)
        call calc_inviscd_sp(nb,nghnode,nghedge, &
                             npvs,mb_pv,neqn,mb_rhs, &
                             sub_limit,sub_intnon,sub_intplt, &
                             sub_fvs,flux_euler,dn_via_ehen4cc, &
                             chkid,chklim)
    case(nderive_ehen6)
        call calc_inviscd_sp(nb,nghnode,nghedge, &
                             npvs,mb_pv,neqn,mb_rhs, &
                             sub_limit,sub_intnon,sub_intplt, &
                             sub_fvs,flux_euler,dn_via_ehen6cc, &
                             chkid,chklim)
    case(nderive_ehen6e)
        call calc_inviscd_sp(nb,nghnode,nghedge, &
                             npvs,mb_pv,neqn,mb_rhs, &
                             sub_limit,sub_intnon,sub_intplt, &
                             sub_fvs,flux_euler,dn_via_ehen6ecc, &
                             chkid,chklim)
    case(nderive_ehcs6)
        call calc_inviscd_sp(nb,nghnode,nghedge, &
                             npvs,mb_pv,neqn,mb_rhs, &
                             sub_limit,sub_intnon,sub_intplt, &
                             sub_fvs,flux_euler,dn_via_ehcs6cc, &
                             chkid,chklim)
    case(nderive_ehcs6e)
        call calc_inviscd_sp(nb,nghnode,nghedge, &
                             npvs,mb_pv,neqn,mb_rhs, &
                             sub_limit,sub_intnon,sub_intplt, &
                             sub_fvs,flux_euler,dn_via_ehcs6ecc, &
                             chkid,chklim)
    case(nderive_ehen8)
        call calc_inviscd_sp(nb,nghnode,nghedge, &
                             npvs,mb_pv,neqn,mb_rhs, &
                             sub_limit,sub_intnon,sub_intplt, &
                             sub_fvs,flux_euler,dn_via_ehen8cc, &
                             chkid,chklim)
    case(nderive_ehen8e)
        call calc_inviscd_sp(nb,nghnode,nghedge, &
                             npvs,mb_pv,neqn,mb_rhs, &
                             sub_limit,sub_intnon,sub_intplt, &
                             sub_fvs,flux_euler,dn_via_ehen8ecc, &
                             chkid,chklim)
    case(nderive_scsl4)
        call calc_inviscd_sp(nb,nghnode,nghedge, &
                             npvs,mb_pv,neqn,mb_rhs, &
                             sub_limit,sub_intnon,sub_intplt, &
                             sub_fvs,flux_euler,dn_via_scsl4cc, &
                             chkid,chklim)
    case(nderive_scsl4e)
        call calc_inviscd_sp(nb,nghnode,nghedge, &
                             npvs,mb_pv,neqn,mb_rhs, &
                             sub_limit,sub_intnon,sub_intplt, &
                             sub_fvs,flux_euler,dn_via_scsl4ecc, &
                             chkid,chklim)
    case(nderive_scsl6)
        call calc_inviscd_sp(nb,nghnode,nghedge, &
                             npvs,mb_pv,neqn,mb_rhs, &
                             sub_limit,sub_intnon,sub_intplt, &
                             sub_fvs,flux_euler,dn_via_scsl6cc, &
                             chkid,chklim)
    case(nderive_scsl6e)
        call calc_inviscd_sp(nb,nghnode,nghedge, &
                             npvs,mb_pv,neqn,mb_rhs, &
                             sub_limit,sub_intnon,sub_intplt,  &
                             sub_fvs,flux_euler,dn_via_scsl6ecc, &
                             chkid,chklim)
    case default
    end select

end subroutine invscd_scheme_sp

!todo SCM 修改插值变量
subroutine calc_inviscd(nb,ngn,nge, &
                        npvs,mb_pv,neqn,mb_rhs, &
                        sub_limit,sub_intnon,sub_intplt, &
                        sub_fvs,sub_flux,sub_scheme, &
                        chkid,chklim)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nfsf_con_d1nint,nfsf_con_d1int
    use mod_constants, only : nfsf_con_d1der,nsw_dir_close
    use mod_constants, only : nbc_inter_scheme,nbc_intbc_scheme
    use mod_constants, only : nflux_roe_prec,nvis_euler,nscmp_non
    use mod_datatypes, only : var_block_t,fld_array_t
    use mod_variables, only : enfix,nsw_kdir,nflux,nvis,nscmp
    use mod_variables, only : poo
    use mod_fieldvars, only : mb_top,mb_fsf,mb_sxyz,mb_vol,mb_vsl,mb_vst
    use mod_openmp
    implicit none
    integer(kind_int), intent(in) :: nb
    integer(kind_int), intent(in) :: ngn,nge
    integer(kind_int), intent(in) :: npvs
    type(var_block_t),    pointer :: mb_pv(:)
    integer(kind_int), intent(in) :: neqn
    type(var_block_t),    pointer :: mb_rhs(:)
    real(kind_real),     external :: sub_limit
    integer(kind_int), intent(in) :: chkid(npvs)
    real(kind_real),   intent(in) :: chklim(npvs)
    external                      :: sub_intnon,sub_intplt
    external                      :: sub_fvs,sub_flux
    external                      :: sub_scheme
    integer(kind_int)          :: i,j,k,m,nfs,nfe,ierr
    integer(kind_int)          :: ni,nj,nk,stn,edn,ste,ede
    integer(kind_int)          :: stn0,edn0,ste0,ede0
    integer(kind_int)          :: nkst,nked,njed
    integer(kind_int)          :: ndec,ndecmax,ijkdec(3,3)
    integer(kind_int)          :: idec,jdec,kdec
    real(kind_real)            :: nx,ny,nz,nt,efix
    real(kind_real)            :: vl(1:npvs),vr(1:npvs)
    real(kind_real)            :: v(1:npvs),f(1:neqn)
    type(fld_array_t), pointer :: sxyz(:),pv(:),rhs(:),vol(:),vsl(:),vst(:)
    type(fld_array_t), pointer :: fsfs(:),fsfe(:)
    real(kind_real),   pointer :: pvn(:,:),pvl(:,:),pvr(:,:)
    real(kind_real),   pointer :: sn(:,:),se(:,:)
    real(kind_real),   pointer :: fn(:,:),fc(:,:),dn(:,:)
    real(kind_real)            :: length,visl,vist
    
    length = 0.

    efix = enfix

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

    ndecmax     = 0
    ijkdec(:,:) = 0

    sxyz => mb_sxyz(nb)%fld
    pv   => mb_pv(nb)%fld
    rhs  => mb_rhs(nb)%fld
    vol  => mb_vol(nb)%fld

    ! I-direction
    idec = 0

    fsfs => mb_fsf(nb,1)%fld
    fsfe => mb_fsf(nb,2)%fld
    stn = 1 - ngn
    edn = ni + ngn
    ste = -nge
    ede = ni + nge
!$OMP parallel private(k,j,i,m,pvn,pvl,pvr,sn,se,fn,fc,dn,nfs,nfe,stn0,ste0,edn0,ede0,v,vr,vl,nx,ny,nz,f)
    allocate(pvn(stn:edn,1:npvs), stat=ierr)
    allocate(pvl(stn:edn,1:npvs), stat=ierr)
    allocate(pvr(stn:edn,1:npvs), stat=ierr)
    allocate(sn(stn:edn,1:3), stat=ierr)
    allocate(se(stn:edn,1:3), stat=ierr)
    allocate(fn(stn:edn,1:neqn), stat=ierr)
    allocate(fc(stn:edn,1:neqn), stat=ierr)
    allocate(dn(stn:edn,1:neqn), stat=ierr)
!$OMP do
    do k=nkst,nked
    do j=1,nj
        do m=1,npvs
            do i=stn,edn
                pvn(i,m) = pv(m)%r3d(i,j,k)
            end do
        end do
        
        !> SCM-P 将原始变量转化成新的SCM-P的变量 (ρ,u,v,w,p) ---> (ρ,u,v,w,p')
        if (nscmp > nscmp_non) then
            do i=stn,edn
                pvn(i,npvs) = pvn(i,npvs) - poo     !> Done
            end do
        end if

        do m=1,3
            do i=stn,edn
               sn(i,m) = sxyz(m)%r3d(i,j,k)
            end do
        end do
        nfs = fsfs(nfsf_con_d1int)%i3d(1 ,j,k)
        nfe = fsfe(nfsf_con_d1int)%i3d(ni,j,k)
        call sub_intplt(ni,1,3,ngn,sn,nfs,nfe,nge,se)

        nfs = fsfs(nfsf_con_d1nint)%i3d(1 ,j,k)
        nfe = fsfe(nfsf_con_d1nint)%i3d(ni,j,k)
        call sub_intnon(ni,1,npvs,ngn,pvn,sn, &
                        nfs,nfe, &
                        nge,pvl,pvr,sub_limit, &
                        chkid,chklim,ndec)
#ifdef OMP_IMP
#else
        if (ndec > 0) then
            idec    = idec + 1
            ndecmax = ndecmax + ndec
            if (idec == 1) then
                ijkdec(:,1) = (/0,j,k/)
            end if
        end if
#endif

        nfs = fsfs(nfsf_con_d1der)%i3d(1 ,j,k)
        nfe = fsfe(nfsf_con_d1der)%i3d(ni,j,k)
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
            do m=1,npvs
                v(m) = pvn(i,m)
            end do
            !!nt = sn(i,0)
            nx = sn(i,1)
            ny = sn(i,2)
            nz = sn(i,3)

            call sub_flux(1,npvs,v,nt,nx,ny,nz,1,neqn,f)

            do m=1,neqn
                fn(i,m) = f(m)
            end do
        end do

        do i=ste0,ede0
            do m=1,npvs
                vl(m) = pvl(i,m)
                vr(m) = pvr(i,m)
            end do
            !!nt = se(i,0)
            nx = se(i,1)
            ny = se(i,2)
            nz = se(i,3)
            
            length = vol(1)%r3d(i,j,k)**(1.0/3.0)

            if (nflux == nflux_roe_prec) then
                call sub_fvs(1,npvs,vl,vr,length,nt,nx,ny,nz,1,neqn,f,efix)
            else
                call sub_fvs(1,npvs,vl,vr,nt,nx,ny,nz,1,neqn,f,efix)
            end if

            do m=1,neqn
                fc(i,m) = f(m)
            end do
        end do

        call sub_scheme(ni,1,neqn,ngn,fn,nge,fc,nfs,nfe,dn)

        do m=1,neqn
            do i=1,ni
                rhs(m)%r3d(i,j,k) = rhs(m)%r3d(i,j,k) - dn(i,m)
            end do
        end do
    end do
    end do
!$OMP end do nowait
    deallocate(dn, stat=ierr)
    deallocate(fc, stat=ierr)
    deallocate(fn, stat=ierr)
    deallocate(sn, stat=ierr)
    deallocate(se, stat=ierr)
    deallocate(pvr, stat=ierr)
    deallocate(pvl, stat=ierr)
    deallocate(pvn, stat=ierr)
!$OMP end parallel

    ! J-direction
    jdec = 0

    fsfs => mb_fsf(nb,3)%fld
    fsfe => mb_fsf(nb,4)%fld
    stn = 1 - ngn
    edn = nj + ngn
    ste = -nge
    ede = nj + nge
!$OMP parallel private(k,j,i,m,pvn,pvl,pvr,sn,se,fn,fc,dn,nfs,nfe,stn0,ste0,edn0,ede0,v,vr,vl,nx,ny,nz,f)
    allocate(pvn(stn:edn,1:npvs), stat=ierr)
    allocate(pvl(stn:edn,1:npvs), stat=ierr)
    allocate(pvr(stn:edn,1:npvs), stat=ierr)
    allocate(sn(stn:edn,1:3), stat=ierr)
    allocate(se(stn:edn,1:3), stat=ierr)
    allocate(fn(stn:edn,1:neqn), stat=ierr)
    allocate(fc(stn:edn,1:neqn), stat=ierr)
    allocate(dn(stn:edn,1:neqn), stat=ierr)
!$OMP do
    do k=nkst,nked
    do i=1,ni
        do m=1,npvs
            do j=stn,edn
                pvn(j,m) = pv(m)%r3d(i,j,k)
            end do
        end do
        
        !> SCM-P 将原始变量转化成新的SCM-P的变量 (ρ,u,v,w,p) ---> (ρ,u,v,w,p')
        if (nscmp > nscmp_non) then
            do j=stn,edn
                pvn(j,npvs) = pvn(j,npvs) - poo     !> Done
            end do
        end if

        do m=1,3
            do j=stn,edn
               sn(j,m) = sxyz(m+3)%r3d(i,j,k)
            end do
        end do
        nfs = fsfs(nfsf_con_d1int)%i3d(i, 1,k)
        nfe = fsfe(nfsf_con_d1int)%i3d(i,nj,k)
        call sub_intplt(nj,1,3,ngn,sn,nfs,nfe,nge,se)

        nfs = fsfs(nfsf_con_d1nint)%i3d(i, 1,k)
        nfe = fsfe(nfsf_con_d1nint)%i3d(i,nj,k)
        call sub_intnon(nj,1,npvs,ngn,pvn,sn, &
                        nfs,nfe, &
                        nge,pvl,pvr,sub_limit, &
                        chkid,chklim,ndec)
#ifdef OMP_IMP
#else
        if (ndec > 0) then
            jdec    = jdec + 1
            ndecmax = ndecmax + ndec
            if (jdec == 1) then
                ijkdec(:,2) = (/i,0,k/)
            end if
        end if
#endif

        nfs = fsfs(nfsf_con_d1der)%i3d(i, 1,k)
        nfe = fsfe(nfsf_con_d1der)%i3d(i,nj,k)
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
            do m=1,npvs
                v(m) = pvn(j,m)
            end do
            !!nt = sn(j,0)
            nx = sn(j,1)
            ny = sn(j,2)
            nz = sn(j,3)
            call sub_flux(1,npvs,v,nt,nx,ny,nz,1,neqn,f)

            do m=1,neqn
                fn(j,m) = f(m)
            end do
        end do

        do j=ste0,ede0
            do m=1,npvs
                vl(m) = pvl(j,m)
                vr(m) = pvr(j,m)
            end do
            !!nt = se(i,0)
            nx = se(j,1)
            ny = se(j,2)
            nz = se(j,3)
            
            length = vol(1)%r3d(i,j,k)**(1.0/3.0)

            if (nflux == nflux_roe_prec) then
                call sub_fvs(1,npvs,vl,vr,length,nt,nx,ny,nz,1,neqn,f,efix)
            else
                call sub_fvs(1,npvs,vl,vr,nt,nx,ny,nz,1,neqn,f,efix)
            end if

            do m=1,neqn
                fc(j,m) = f(m)
            end do
        end do

        call sub_scheme(nj,1,neqn,ngn,fn,nge,fc,nfs,nfe,dn)

        do m=1,neqn
            do j=1,nj
                rhs(m)%r3d(i,j,k) = rhs(m)%r3d(i,j,k) - dn(j,m)
            end do
        end do
    end do
    end do
!$OMP end do nowait
    deallocate(dn, stat=ierr)
    deallocate(fc, stat=ierr)
    deallocate(fn, stat=ierr)
    deallocate(sn, stat=ierr)
    deallocate(se, stat=ierr)
    deallocate(pvr, stat=ierr)
    deallocate(pvl, stat=ierr)
    deallocate(pvn, stat=ierr)
!$OMP end parallel

    ! K-direction
    kdec = 0

    fsfs => mb_fsf(nb,5)%fld
    fsfe => mb_fsf(nb,6)%fld
    stn = 1 - ngn
    edn = nk + ngn
    ste = -nge
    ede = nk + nge
!$OMP parallel private(k,j,i,m,pvn,pvr,pvl,sn,nfs,nfe,stn0,ste0,edn0,ede0,v,vr,vl,nx,ny,nz,f,se,fc,fn,dn)
    allocate(pvn(stn:edn,1:npvs), stat=ierr)
    allocate(pvl(stn:edn,1:npvs), stat=ierr)
    allocate(pvr(stn:edn,1:npvs), stat=ierr)
    allocate(sn(stn:edn,1:3), stat=ierr)
    allocate(se(stn:edn,1:3), stat=ierr)
    allocate(fn(stn:edn,1:neqn), stat=ierr)
    allocate(fc(stn:edn,1:neqn), stat=ierr)
    allocate(dn(stn:edn,1:neqn), stat=ierr)
!$OMP do
    do j=1,njed
    do i=1,ni
        do m=1,npvs
            do k=stn,edn
                pvn(k,m) = pv(m)%r3d(i,j,k)
            end do
        end do
        
        !> SCM-P 将原始变量转化成新的SCM-P的变量 (ρ,u,v,w,p) ---> (ρ,u,v,w,p')
        if (nscmp > nscmp_non) then
            do k=stn,edn
                pvn(k,npvs) = pvn(k,npvs) - poo     !> Done
            end do
        end if

        do m=1,3
            do k=stn,edn
               sn(k,m) = sxyz(m+6)%r3d(i,j,k)
            end do
        end do
        nfs = fsfs(nfsf_con_d1int)%i3d(i,j, 1)
        nfe = fsfe(nfsf_con_d1int)%i3d(i,j,nk)
        call sub_intplt(nk,1,3,ngn,sn,nfs,nfe,nge,se)

        nfs = fsfs(nfsf_con_d1nint)%i3d(i,j, 1)
        nfe = fsfe(nfsf_con_d1nint)%i3d(i,j,nk)
        call sub_intnon(nk,1,npvs,ngn,pvn,sn, &
                        nfs,nfe, &
                        nge,pvl,pvr,sub_limit, &
                        chkid,chklim,ndec)
#ifdef OMP_IMP
#else
        if (ndec > 0) then
            kdec    = kdec + 1
            ndecmax = ndecmax + ndec
            if (kdec == 1) then
                ijkdec(:,3) = (/i,j,0/)
            end if
        end if
#endif

        nfs = fsfs(nfsf_con_d1der)%i3d(i,j, 1)
        nfe = fsfe(nfsf_con_d1der)%i3d(i,j,nk)
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
            do m=1,npvs
                v(m) = pvn(k,m)
            end do
            !!nt = sn(k,0)
            nx = sn(k,1)
            ny = sn(k,2)
            nz = sn(k,3)
            call sub_flux(1,npvs,v,nt,nx,ny,nz,1,neqn,f)

            do m=1,neqn
                fn(k,m) = f(m)
            end do
        end do

        do k=ste0,ede0
            do m=1,npvs
                vl(m) = pvl(k,m)
                vr(m) = pvr(k,m)
            end do
            !!nt = se(k,0)
            nx = se(k,1)
            ny = se(k,2)
            nz = se(k,3)

            length = vol(1)%r3d(i,j,k)**(1.0/3.0)

            if (nflux == nflux_roe_prec) then
                call sub_fvs(1,npvs,vl,vr,length,nt,nx,ny,nz,1,neqn,f,efix)
            else
                call sub_fvs(1,npvs,vl,vr,nt,nx,ny,nz,1,neqn,f,efix)
            end if

            do m=1,neqn
                fc(k,m) = f(m)
            end do
        end do

        call sub_scheme(nk,1,neqn,ngn,fn,nge,fc,nfs,nfe,dn)

        do m=1,neqn
            do k=1,nk
                rhs(m)%r3d(i,j,k) = rhs(m)%r3d(i,j,k) - dn(k,m)
            end do
        end do
    end do
    end do
!$OMP end do nowait
    deallocate(dn, stat=ierr)
    deallocate(fc, stat=ierr)
    deallocate(fn, stat=ierr)
    deallocate(sn, stat=ierr)
    deallocate(se, stat=ierr)
    deallocate(pvr, stat=ierr)
    deallocate(pvl, stat=ierr)
    deallocate(pvn, stat=ierr)
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

    if (ndecmax > 0) then
        write(*,'(a,i5)') " Warning: variables interpolation failure:",ndecmax
        write(*,'(3x,i5,3(1x,"(",3(1x,i5),")"))') nb,ijkdec
    end if

end subroutine calc_inviscd

subroutine calc_inviscd_sp(nb,ngn,nge, &
                           npvs,mb_pv,neqn,mb_rhs, &
                           sub_limit,sub_intnon,sub_intplt, &
                           sub_fvs,sub_flux,sub_scheme, &
                           chkid,chklim)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nfsf_con_d1nint,nfsf_con_d1int
    use mod_constants, only : nfsf_con_d1der,nsw_dir_close
    use mod_constants, only : nbc_inter_scheme,nbc_intbc_scheme
    use mod_constants, only : nflux_roe_prec,nvis_euler,nscmp_non
    use mod_datatypes, only : var_block_t,fld_array_t
    use mod_variables, only : enfix,nsw_kdir,nflux,nvis,nscmp
    use mod_variables, only : poo
    use mod_fieldvars, only : mb_topsp,mb_fsfsp,mb_fsffp,mb_pvfp,mb_sxyzsp,mb_volsp,mb_vsl,mb_vst
    use mod_openmp
    implicit none
    integer(kind_int), intent(in) :: nb
    integer(kind_int), intent(in) :: ngn,nge
    integer(kind_int), intent(in) :: npvs
    type(var_block_t),    pointer :: mb_pv(:)
    integer(kind_int), intent(in) :: neqn
    type(var_block_t),    pointer :: mb_rhs(:)
    real(kind_real),     external :: sub_limit
    integer(kind_int), intent(in) :: chkid(npvs)
    real(kind_real),   intent(in) :: chklim(npvs)
    external                      :: sub_intnon,sub_intplt
    external                      :: sub_fvs,sub_flux
    external                      :: sub_scheme
    integer(kind_int)          :: i,j,k,m,nfs,nfe,ierr
    integer(kind_int)          :: ni,nj,nk,stn,edn,ste,ede
    integer(kind_int)          :: stn0,edn0,ste0,ede0
    integer(kind_int)          :: nkst,nked,njed
    integer(kind_int)          :: ndec,ndecmax,ijkdec(3,3)
    integer(kind_int)          :: idec,jdec,kdec
    integer(kind_int)          :: nr,bctype
    real(kind_real)            :: nx,ny,nz,nt,efix
    real(kind_real)            :: vl(1:npvs),vr(1:npvs)
    real(kind_real)            :: v(1:npvs),f(1:neqn),varL(1:npvs),varR(1:npvs)
    type(fld_array_t), pointer :: sxyzcc(:),sxyzsp(:),pv(:),pvfp(:),rhs(:),volsp(:),vsl(:),vst(:)
    type(fld_array_t), pointer :: fsfs(:),fsfe(:),fsfps(:),fsfpe(:)
    real(kind_real),   pointer :: pvn(:,:),pvl(:,:),pvr(:,:)
    real(kind_real),   pointer :: sn(:,:),se(:,:)
    real(kind_real),   pointer :: fn(:,:),fc(:,:),dn(:,:)
    real(kind_real)            :: length,visl,vist
    
    length = 0.

    efix = enfix

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

    ndecmax     = 0
    ijkdec(:,:) = 0

    sxyzsp => mb_sxyzsp(nb)%fld
    pv   => mb_pv(nb)%fld
    pvfp => mb_pvfp(nb)%fld
    rhs  => mb_rhs(nb)%fld
    volsp  => mb_volsp(nb)%fld

    ! I-direction
    idec = 0

    fsfs => mb_fsfsp(nb,1)%fld
    fsfe => mb_fsfsp(nb,2)%fld
    fsfps=> mb_fsffp(nb,1)%fld    
    fsfpe=> mb_fsffp(nb,2)%fld    
    stn = 1 - ngn
    edn = ni + ngn
    ste = 1 - nge
    ede = ni + 1 + nge
!$OMP parallel private(k,j,i,m,pvn,pvl,pvr,sn,se,fn,fc,dn,nfs,nfe,stn0,ste0,edn0,ede0,v,vr,vl,nx,ny,nz,f)
    allocate(pvn(stn:edn,1:npvs), stat=ierr)
    allocate(pvl(stn:edn,1:npvs), stat=ierr)
    allocate(pvr(stn:edn,1:npvs), stat=ierr)
    allocate(sn(stn:edn,1:3), stat=ierr)
    allocate(se(stn:edn,1:3), stat=ierr)
    allocate(fn(stn:edn,1:neqn), stat=ierr)
    allocate(fc(stn:edn,1:neqn), stat=ierr)
    allocate(dn(stn:edn,1:neqn), stat=ierr)
!$OMP do
    do k=nkst,nked
    do j=1,nj
        do m=1,npvs
            do i=stn,edn
                pvn(i,m) = pv(m)%r3d(i,j,k)
            end do
        end do
        
        !> SCM-P 将原始变量转化成新的SCM-P的变量 (ρ,u,v,w,p) ---> (ρ,u,v,w,p')
        if (nscmp > nscmp_non) then
            do i=stn,edn
                pvn(i,npvs) = pvn(i,npvs) - poo     !> Done
            end do
        end if

        do m=1,3
            do i=stn,edn
               sn(i,m) = sxyzsp(m)%r3d(i,j,k)
            end do
        end do
        nfs = fsfs(nfsf_con_d1int)%i3d(1 ,j,k)
        nfe = fsfe(nfsf_con_d1int)%i3d(ni,j,k)
        call sub_intplt(ni,1,3,ngn,sn,nfs,nfe,nge,se)
        !do m=1,3
        !    do i=ste,ede
        !       se(i,m) = sxyzcc(m)%r3d(i,j,k)
        !    end do
        !end do        

        nfs = fsfs(nfsf_con_d1nint)%i3d(1 ,j,k)
        nfe = fsfe(nfsf_con_d1nint)%i3d(ni,j,k)
        call sub_intnon(ni,1,npvs,ngn,pvn,sn, &
                        nfs,nfe, &
                        nge,pvl,pvr,sub_limit, &
                        chkid,chklim,ndec)
        do m=1,npvs
            do i=1,ni+1
                pvfp(m)%r3d(i,j,k) = 0.5*(pvl(i,m) + pvr(i,m))
            end do
        end do
#ifdef OMP_IMP
#else
        if (ndec > 0) then
            idec    = idec + 1
            ndecmax = ndecmax + ndec
            if (idec == 1) then
                ijkdec(:,1) = (/0,j,k/)
            end if
        end if
#endif

        nfs = fsfs(nfsf_con_d1der)%i3d(1 ,j,k)
        nfe = fsfe(nfsf_con_d1der)%i3d(ni,j,k)
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
            do m=1,npvs
                v(m) = pvn(i,m)
            end do
            !!nt = sn(i,0)
            nx = sn(i,1)
            ny = sn(i,2)
            nz = sn(i,3)

            call sub_flux(1,npvs,v,nt,nx,ny,nz,1,neqn,f)

            do m=1,neqn
                fn(i,m) = f(m)
            end do
        end do

        do i=ste0,ede0
            do m=1,npvs
                vl(m) = pvl(i,m)
                vr(m) = pvr(i,m)
            end do
            !!nt = se(i,0)
            nx = se(i,1)
            ny = se(i,2)
            nz = se(i,3)
            
            if(i==1) then              
                nr     = fsfps(1)%i3d(1,j,k)
                bctype = fsfps(2)%i3d(1,j,k)
                if(bctype /= -1) then
                    call set_boundary_vars(1,npvs,nb,nr,bctype,i,j,k,nx,ny,nz,varL,varR)
                    do m=1,neqn
                        vl(m) = varL(m)
                        vr(m) = varR(m)
                        pvfp(m)%r3d(i,j,k) = 0.5*(vl(m) + vr(m))
                    end do
                end if
            end if
            if(i==ni+1) then
                nr     = fsfpe(1)%i3d(ni+1,j,k)
                bctype = fsfpe(2)%i3d(ni+1,j,k)
                if(bctype /= -1) then
                    call set_boundary_vars(1,npvs,nb,nr,bctype,i,j,k,nx,ny,nz,varL,varR)
                    do m=1,neqn
                        vl(m) = varL(m)
                        vr(m) = varR(m)
                        pvfp(m)%r3d(i,j,k) = 0.5*(vl(m) + vr(m))
                    end do
                end if
            end if   
            
            length = volsp(1)%r3d(i,j,k)**(1.0/3.0)

            if (nflux == nflux_roe_prec) then
                call sub_fvs(1,npvs,vl,vr,length,nt,nx,ny,nz,1,neqn,f,efix)
            else
                call sub_fvs(1,npvs,vl,vr,nt,nx,ny,nz,1,neqn,f,efix)
            end if

            do m=1,neqn
                fc(i,m) = f(m)
            end do
        end do

        call sub_scheme(ni,1,neqn,ngn,fn,nge,fc,nfs,nfe,dn)

        do m=1,neqn
            do i=1,ni
                rhs(m)%r3d(i,j,k) = rhs(m)%r3d(i,j,k) - dn(i,m)
            end do
        end do
    end do
    end do
!$OMP end do nowait
    deallocate(dn, stat=ierr)
    deallocate(fc, stat=ierr)
    deallocate(fn, stat=ierr)
    deallocate(sn, stat=ierr)
    deallocate(se, stat=ierr)
    deallocate(pvr, stat=ierr)
    deallocate(pvl, stat=ierr)
    deallocate(pvn, stat=ierr)
!$OMP end parallel

    ! J-direction
    jdec = 0

    fsfs => mb_fsfsp(nb,3)%fld
    fsfe => mb_fsfsp(nb,4)%fld
    fsfps=> mb_fsffp(nb,3)%fld
    fsfpe=> mb_fsffp(nb,4)%fld    
    stn = 1 - ngn
    edn = nj + ngn
    ste = 1 - nge
    ede = nj + 1 + nge
!$OMP parallel private(k,j,i,m,pvn,pvl,pvr,sn,se,fn,fc,dn,nfs,nfe,stn0,ste0,edn0,ede0,v,vr,vl,nx,ny,nz,f)
    allocate(pvn(stn:edn,1:npvs), stat=ierr)
    allocate(pvl(stn:edn,1:npvs), stat=ierr)
    allocate(pvr(stn:edn,1:npvs), stat=ierr)
    allocate(sn(stn:edn,1:3), stat=ierr)
    allocate(se(stn:edn,1:3), stat=ierr)
    allocate(fn(stn:edn,1:neqn), stat=ierr)
    allocate(fc(stn:edn,1:neqn), stat=ierr)
    allocate(dn(stn:edn,1:neqn), stat=ierr)
!$OMP do
    do k=nkst,nked
    do i=1,ni
        do m=1,npvs
            do j=stn,edn
                pvn(j,m) = pv(m)%r3d(i,j,k)
            end do
        end do
        
        !> SCM-P 将原始变量转化成新的SCM-P的变量 (ρ,u,v,w,p) ---> (ρ,u,v,w,p')
        if (nscmp > nscmp_non) then
            do j=stn,edn
                pvn(j,npvs) = pvn(j,npvs) - poo     !> Done
            end do
        end if

        do m=1,3
            do j=stn,edn
               sn(j,m) = sxyzsp(m+3)%r3d(i,j,k)
            end do
        end do
        nfs = fsfs(nfsf_con_d1int)%i3d(i, 1,k)
        nfe = fsfe(nfsf_con_d1int)%i3d(i,nj,k)
        call sub_intplt(nj,1,3,ngn,sn,nfs,nfe,nge,se)
        !do m=1,3
        !    do j=ste,ede
        !       se(j,m) = sxyzcc(m+3)%r3d(i,j,k)
        !    end do
        !end do        

        nfs = fsfs(nfsf_con_d1nint)%i3d(i, 1,k)
        nfe = fsfe(nfsf_con_d1nint)%i3d(i,nj,k)
        call sub_intnon(nj,1,npvs,ngn,pvn,sn, &
                        nfs,nfe, &
                        nge,pvl,pvr,sub_limit, &
                        chkid,chklim,ndec)
       do m=1,npvs
            do j=1,nj+1
                pvfp(m)%r3d(i,j,k) = 0.5*(pvl(j,m) + pvr(j,m))
            end do
       end do
#ifdef OMP_IMP
#else
        if (ndec > 0) then
            jdec    = jdec + 1
            ndecmax = ndecmax + ndec
            if (jdec == 1) then
                ijkdec(:,2) = (/i,0,k/)
            end if
        end if
#endif

        nfs = fsfs(nfsf_con_d1der)%i3d(i, 1,k)
        nfe = fsfe(nfsf_con_d1der)%i3d(i,nj,k)
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
            do m=1,npvs
                v(m) = pvn(j,m)
            end do
            !!nt = sn(j,0)
            nx = sn(j,1)
            ny = sn(j,2)
            nz = sn(j,3)
            call sub_flux(1,npvs,v,nt,nx,ny,nz,1,neqn,f)

            do m=1,neqn
                fn(j,m) = f(m)
            end do
        end do

        do j=ste0,ede0
            do m=1,npvs
                vl(m) = pvl(j,m)
                vr(m) = pvr(j,m)
            end do
            !!nt = se(i,0)
            nx = se(j,1)
            ny = se(j,2)
            nz = se(j,3)
            
            if(j==1) then
                nr     = fsfps(1)%i3d(i,1,k)
                bctype = fsfps(2)%i3d(i,1,k)
                if(bctype /= -1) then
                    call set_boundary_vars(1,npvs,nb,nr,bctype,i,j,k,nx,ny,nz,varL,varR)
                    do m=1,neqn
                        vl(m) = varL(m)
                        vr(m) = varR(m)
                        pvfp(m)%r3d(i,j,k) = 0.5*(vl(m) + vr(m))
                    end do
                end if
            end if
            if(j==nj+1) then                
                nr     = fsfpe(1)%i3d(i,nj+1,k)
                bctype = fsfpe(2)%i3d(i,nj+1,k)
                if(bctype /= -1) then
                    call set_boundary_vars(1,npvs,nb,nr,bctype,i,j,k,nx,ny,nz,varL,varR)
                    do m=1,neqn
                        vl(m) = varL(m)
                        vr(m) = varR(m)
                        pvfp(m)%r3d(i,j,k) = 0.5*(vl(m) + vr(m))
                    end do
                end if
            end if
            
            length = volsp(1)%r3d(i,j,k)**(1.0/3.0)

            if (nflux == nflux_roe_prec) then
                call sub_fvs(1,npvs,vl,vr,length,nt,nx,ny,nz,1,neqn,f,efix)
            else
                call sub_fvs(1,npvs,vl,vr,nt,nx,ny,nz,1,neqn,f,efix)
            end if

            do m=1,neqn
                fc(j,m) = f(m)
            end do
        end do

        call sub_scheme(nj,1,neqn,ngn,fn,nge,fc,nfs,nfe,dn)

        do m=1,neqn
            do j=1,nj
                rhs(m)%r3d(i,j,k) = rhs(m)%r3d(i,j,k) - dn(j,m)
            end do
        end do
    end do
    end do
!$OMP end do nowait
    deallocate(dn, stat=ierr)
    deallocate(fc, stat=ierr)
    deallocate(fn, stat=ierr)
    deallocate(sn, stat=ierr)
    deallocate(se, stat=ierr)
    deallocate(pvr, stat=ierr)
    deallocate(pvl, stat=ierr)
    deallocate(pvn, stat=ierr)
!$OMP end parallel

    ! K-direction
    kdec = 0

    fsfs => mb_fsfsp(nb,5)%fld
    fsfe => mb_fsfsp(nb,6)%fld
    fsfps=> mb_fsffp(nb,5)%fld
    fsfpe=> mb_fsffp(nb,6)%fld    
    stn = 1 - ngn
    edn = nk + ngn
    ste = 1 - nge
    ede = nk + 1 + nge
!$OMP parallel private(k,j,i,m,pvn,pvr,pvl,sn,nfs,nfe,stn0,ste0,edn0,ede0,v,vr,vl,nx,ny,nz,f,se,fc,fn,dn)
    allocate(pvn(stn:edn,1:npvs), stat=ierr)
    allocate(pvl(stn:edn,1:npvs), stat=ierr)
    allocate(pvr(stn:edn,1:npvs), stat=ierr)
    allocate(sn(stn:edn,1:3), stat=ierr)
    allocate(se(stn:edn,1:3), stat=ierr)
    allocate(fn(stn:edn,1:neqn), stat=ierr)
    allocate(fc(stn:edn,1:neqn), stat=ierr)
    allocate(dn(stn:edn,1:neqn), stat=ierr)
!$OMP do
    do j=1,njed
    do i=1,ni
        do m=1,npvs
            do k=stn,edn
                pvn(k,m) = pv(m)%r3d(i,j,k)
            end do
        end do
        
        !> SCM-P 将原始变量转化成新的SCM-P的变量 (ρ,u,v,w,p) ---> (ρ,u,v,w,p')
        if (nscmp > nscmp_non) then
            do k=stn,edn
                pvn(k,npvs) = pvn(k,npvs) - poo     !> Done
            end do
        end if

        do m=1,3
            do k=stn,edn
               sn(k,m) = sxyzsp(m+6)%r3d(i,j,k)
            end do
        end do
        nfs = fsfs(nfsf_con_d1int)%i3d(i,j, 1)
        nfe = fsfe(nfsf_con_d1int)%i3d(i,j,nk)
        call sub_intplt(nk,1,3,ngn,sn,nfs,nfe,nge,se)
        !do m=1,3
        !    do k=ste,ede
        !       se(k,m) = sxyzcc(m+6)%r3d(i,j,k)
        !    end do
        !end do         

        nfs = fsfs(nfsf_con_d1nint)%i3d(i,j, 1)
        nfe = fsfe(nfsf_con_d1nint)%i3d(i,j,nk)
        call sub_intnon(nk,1,npvs,ngn,pvn,sn, &
                        nfs,nfe, &
                        nge,pvl,pvr,sub_limit, &
                        chkid,chklim,ndec)
       do m=1,npvs
            do k=1,nk+1
                pvfp(m)%r3d(i,j,k) = 0.5*(pvl(k,m) + pvr(k,m))
            end do
       end do        
#ifdef OMP_IMP
#else
        if (ndec > 0) then
            kdec    = kdec + 1
            ndecmax = ndecmax + ndec
            if (kdec == 1) then
                ijkdec(:,3) = (/i,j,0/)
            end if
        end if
#endif

        nfs = fsfs(nfsf_con_d1der)%i3d(i,j, 1)
        nfe = fsfe(nfsf_con_d1der)%i3d(i,j,nk)
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
            do m=1,npvs
                v(m) = pvn(k,m)
            end do
            !!nt = sn(k,0)
            nx = sn(k,1)
            ny = sn(k,2)
            nz = sn(k,3)
            call sub_flux(1,npvs,v,nt,nx,ny,nz,1,neqn,f)

            do m=1,neqn
                fn(k,m) = f(m)
            end do
        end do

        do k=ste0,ede0
            do m=1,npvs
                vl(m) = pvl(k,m)
                vr(m) = pvr(k,m)
            end do
            !!nt = se(k,0)
            nx = se(k,1)
            ny = se(k,2)
            nz = se(k,3)
            
            if(k==1) then
                nr     = fsfps(1)%i3d(i,j,1)
                bctype = fsfps(2)%i3d(i,j,1)
                if(bctype /= -1) then
                    call set_boundary_vars(1,npvs,nb,nr,bctype,i,j,k,nx,ny,nz,varL,varR)
                    do m=1,neqn
                        vl(m) = varL(m)
                        vr(m) = varR(m)
                        pvfp(m)%r3d(i,j,k) = 0.5*(vl(m) + vr(m))
                    end do
                end if
            end if
            if(k==nk+1) then
                nr     = fsfpe(1)%i3d(i,j,nk+1)
                bctype = fsfpe(2)%i3d(i,j,nk+1)
                if(bctype /= -1) then
                    call set_boundary_vars(1,npvs,nb,nr,bctype,i,j,k,nx,ny,nz,varL,varR)
                    do m=1,neqn
                        vl(m) = varL(m)
                        vr(m) = varR(m)
                        pvfp(m)%r3d(i,j,k) = 0.5*(vl(m) + vr(m))
                    end do
                end if
            end if             

            length = volsp(1)%r3d(i,j,k)**(1.0/3.0)

            if (nflux == nflux_roe_prec) then
                call sub_fvs(1,npvs,vl,vr,length,nt,nx,ny,nz,1,neqn,f,efix)
            else
                call sub_fvs(1,npvs,vl,vr,nt,nx,ny,nz,1,neqn,f,efix)
            end if

            do m=1,neqn
                fc(k,m) = f(m)
            end do
        end do

        call sub_scheme(nk,1,neqn,ngn,fn,nge,fc,nfs,nfe,dn)

        do m=1,neqn
            do k=1,nk
                rhs(m)%r3d(i,j,k) = rhs(m)%r3d(i,j,k) - dn(k,m)
            end do
        end do
    end do
    end do
!$OMP end do nowait
    deallocate(dn, stat=ierr)
    deallocate(fc, stat=ierr)
    deallocate(fn, stat=ierr)
    deallocate(sn, stat=ierr)
    deallocate(se, stat=ierr)
    deallocate(pvr, stat=ierr)
    deallocate(pvl, stat=ierr)
    deallocate(pvn, stat=ierr)
!$OMP end parallel

    if (nsw_kdir == nsw_dir_close) then
        do m=1,neqn
            do k=1,nk
            do j=1,nj
            do i=1,ni
                rhs(m)%r3d(i,j,k) = rhs(m)%r3d(i,j,nkst)
                pvfp(m)%r3d(i,j,k) = pvfp(m)%r3d(i,j,nkst)
            end do
            end do
            end do
        end do
    end if

    if (ndecmax > 0) then
        write(*,'(a,i5)') " Warning: variables interpolation failure:",ndecmax
        write(*,'(3x,i5,3(1x,"(",3(1x,i5),")"))') nb,ijkdec
    end if

end subroutine calc_inviscd_sp

