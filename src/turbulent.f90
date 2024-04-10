

subroutine turbulent
    use mod_constants, only : nturlhs_lusgs_std_sca
    use mod_constants, only : nturlhs_prsgs_std_sca
    use mod_variables, only : nturlhs
    implicit none

    select case(nturlhs)
    case(nturlhs_lusgs_std_sca)
        call tur_lusgs_std_sca
    case(nturlhs_prsgs_std_sca)
        call tur_prsgs_std_sca
    end select

end subroutine turbulent

subroutine reset_turconsts
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : one,nvis_tur_sa
    use mod_constants, only : nvis_tur_sst,nvis_tur_hst
    use mod_variables, only : nvis,turconsts
    implicit none
    real(kind_real) :: kappa,sigma,cb1,cb2,cv1,cv2
    real(kind_real) :: cw1,cw2,cw3,cv5,sdilim,nulim
    real(kind_real) :: betas,a1,pklim,klim,owlim
    real(kind_real) :: sigk1,sigw1,beta1,alfa1,sigd1
    real(kind_real) :: sigk2,sigw2,beta2,alfa2,sigd2

    select case(nvis)
    case(nvis_tur_sa)
        kappa = 0.41
        sigma = 2.0/3.0
        turconsts(1) = kappa
        turconsts(2) = sigma
       
        cb1 = 0.1355
        cb2 = 0.622
        turconsts(3) = cb1
        turconsts(4) = cb2

        cv1 = 7.1
        cv2 = 5.0
        cv5 = 3.5
        turconsts(5) = cv1
        turconsts(6) = cv2
        turconsts(7) = cv5

        cw1 = cb1/(kappa*kappa) + (one+cb2)/sigma
        cw2 = 0.3
        cw3 = 2.0
        turconsts(8) = cw1
        turconsts(9) = cw2
        turconsts(10) = cw3
       
        sdilim= 20.0
        nulim = 1.0e-8
        turconsts(11) = sdilim
        turconsts(12) = nulim
    case(nvis_tur_sst)
        kappa = 0.41
        betas = 0.09
        turconsts(1) = kappa
        turconsts(2) = betas
        
        sigk1 = 0.85
        sigw1 = 0.5
        beta1 = 0.075
        alfa1 = 0.5532
        turconsts(3) = sigk1
        turconsts(4) = sigw1
        turconsts(5) = beta1
        turconsts(6) = alfa1

        sigk2 = 1.0
        sigw2 = 0.856
        beta2 = 0.0828
        alfa2 = 0.4404
        turconsts(7) = sigk2
        turconsts(8) = sigw2
        turconsts(9) = beta2
        turconsts(10) = alfa2
        
        a1 = 0.31
        turconsts(11) = a1

        pklim = 20.0
        klim  = 1.0e-16
        owlim = 1.0e-8
        turconsts(12) = pklim
        turconsts(13) = klim
        turconsts(14) = owlim
    case(nvis_tur_hst)
        kappa = 0.41
        betas = 0.09
        turconsts(1) = kappa
        turconsts(2) = betas
        
        sigk1 = 1.1
        sigw1 = 0.53
        beta1 = 0.0747
        alfa1 = 0.518
        turconsts(3) = sigk1
        turconsts(4) = sigw1
        turconsts(5) = beta1
        turconsts(6) = alfa1

        sigk2 = 1.1
        sigw2 = 1.0
        beta2 = 0.0828
        alfa2 = 0.44
        turconsts(7) = sigk2
        turconsts(8) = sigw2
        turconsts(9) = beta2
        turconsts(10) = alfa2
        
        a1 = 0.31
        turconsts(11) = a1

        pklim = 20.0
        klim  = 1.0e-16
        owlim = 1.0e-8
        turconsts(12) = pklim
        turconsts(13) = klim
        turconsts(14) = owlim

        sigd1 = 1.0
        sigd2 = 0.4
        turconsts(15) = sigd1
        turconsts(16) = sigd2
    end select

end subroutine reset_turconsts

subroutine tur_invscd(nb)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nlimit_none,nlimit_minmod
    use mod_constants, only : nlimit_vanleer,nlimit_vanalbada
    use mod_variables, only : nlimit
    implicit none
    integer(kind_int), intent(in) :: nb
    real(kind_real), external :: nolimiter,minmod
    real(kind_real), external :: vanleer,vanalbada

    select case(nlimit)
    case(nlimit_none)
        call tur_invscd_intnon(nb,nolimiter)
    case(nlimit_minmod)
        call tur_invscd_intnon(nb,minmod   )
    case(nlimit_vanleer)
        call tur_invscd_intnon(nb,vanleer  )
    case(nlimit_vanalbada)
        call tur_invscd_intnon(nb,vanalbada)
    end select

end subroutine tur_invscd

subroutine tur_invscd_intnon(nb,sub_limit)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nintnon_muscl2pv,nintnon_wcns5pv
    use mod_constants, only : nintnon_wcns5cv
    use mod_variables, only : nturint
    implicit none
    integer(kind_int), intent(in) :: nb
    real(kind_real), external     :: sub_limit
    external :: muscl2pv,wcns5pv

    select case(nturint)
    case(nintnon_muscl2pv)
       call tur_invscd_intplt(nb,sub_limit,muscl2pv)
    case(nintnon_wcns5pv)
       call tur_invscd_intplt(nb,sub_limit,wcns5pv)
    case default

    end select

end subroutine tur_invscd_intnon

subroutine tur_invscd_intplt(nb,sub_limit,sub_intnon)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nintplt_node2,nintplt_node4
    use mod_constants, only : nintplt_node6,nintplt_node6w,nintplt_node8
    use mod_constants, only : nintplt_scsl4,nintplt_scsl6
    use mod_constants, only : nintplt_node6e,nintplt_node8e
    use mod_constants, only : nintplt_scsl4e,nintplt_scsl6e
    use mod_variables, only : nd1int_tur_con
    implicit none
    integer(kind_int), intent(in) :: nb
    real(kind_real), external     :: sub_limit
    external                      :: sub_intnon
    external :: ve_via_node2,ve_via_node4,ve_via_node6,ve_via_node8
    external :: ve_via_scsl4,ve_via_scsl6
    external :: ve_via_node6e,ve_via_node8e
    external :: ve_via_scsl4e,ve_via_scsl6e

    select case(nd1int_tur_con)
    case(nintplt_node2)
        call tur_invscd_scheme(nb,sub_limit,sub_intnon,ve_via_node2)
    case(nintplt_node4)
        call tur_invscd_scheme(nb,sub_limit,sub_intnon,ve_via_node4)
    case(nintplt_node6,nintplt_node6w)
        call tur_invscd_scheme(nb,sub_limit,sub_intnon,ve_via_node6)
    case(nintplt_node6e)
        call tur_invscd_scheme(nb,sub_limit,sub_intnon,ve_via_node6e)
    case(nintplt_node8)
        call tur_invscd_scheme(nb,sub_limit,sub_intnon,ve_via_node8)
    case(nintplt_node8e)
        call tur_invscd_scheme(nb,sub_limit,sub_intnon,ve_via_node8e)
    case(nintplt_scsl4)
        call tur_invscd_scheme(nb,sub_limit,sub_intnon,ve_via_scsl4)
    case(nintplt_scsl4e)
        call tur_invscd_scheme(nb,sub_limit,sub_intnon,ve_via_scsl4e)
    case(nintplt_scsl6)
        call tur_invscd_scheme(nb,sub_limit,sub_intnon,ve_via_scsl6)
    case(nintplt_scsl6e)
        call tur_invscd_scheme(nb,sub_limit,sub_intnon,ve_via_scsl6e)
    case default

    end select

end subroutine tur_invscd_intplt

subroutine tur_invscd_scheme(nb,sub_limit,sub_intnon,sub_intplt)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nvis_tur_sa,nvis_tur_sst,nvis_tur_hst
    use mod_constants, only : nderive_edge2,nderive_edge4
    use mod_constants, only : nderive_edge6,nderive_ehen4
    use mod_constants, only : nderive_ehen6,nderive_ehen8,nderive_ehcs6
    use mod_constants, only : nderive_scsl4,nderive_scsl6
    use mod_constants, only : nderive_ehen6e,nderive_ehen8e,nderive_ehcs6e
    use mod_constants, only : nderive_scsl4e,nderive_scsl6e
    use mod_datatypes, only : var_block_t
    use mod_variables, only : nd1der_tur_con,nghnode,nghedge
    use mod_variables, only : nvis,turconsts
    use mod_fieldvars, only : npvt,mb_pvt,neqt,mb_rhst
    use mod_interface, only : calc_inviscd
    implicit none
    integer(kind_int), intent(in) :: nb
    real(kind_real), external     :: sub_limit
    external                      :: sub_intnon,sub_intplt
    integer(kind_int) :: m
    integer(kind_int) :: chkid(npvt)
    real(kind_real)   :: chklim(npvt)
    external          :: dn_via_edge6,dn_via_edge4,dn_via_edge2
    external          :: dn_via_ehen6,dn_via_ehen4,dn_via_ehen8,dn_via_ehcs6
    external          :: dn_via_scsl4,dn_via_scsl6
    external          :: dn_via_ehen6e,dn_via_ehen8e,dn_via_ehcs6e
    external          :: dn_via_scsl4e,dn_via_scsl6e
    external          :: tur_flux_steger,tur_flux_euler

    chkid(:) = 0
    chkid(1) = 1
    chkid(5:npvt) = 0

    select case(nvis)
    case(nvis_tur_sa)
        chklim(5) = turconsts(12)
    case(nvis_tur_sst,nvis_tur_hst)
        chklim(5) = turconsts(13)
        chklim(6) = turconsts(14)
    end select

    select case(nd1der_tur_con)
    case(nderive_edge2)
        call calc_inviscd(nb,nghnode,nghedge, &
                          npvt,mb_pvt,neqt,mb_rhst, &
                          sub_limit,sub_intnon,sub_intplt, &
                          tur_flux_steger,tur_flux_euler,dn_via_edge2, &
                          chkid,chklim)
    case(nderive_edge4)
        call calc_inviscd(nb,nghnode,nghedge, &
                          npvt,mb_pvt,neqt,mb_rhst, &
                          sub_limit,sub_intnon,sub_intplt, &
                          tur_flux_steger,tur_flux_euler,dn_via_edge4, &
                          chkid,chklim)
    case(nderive_edge6)
        call calc_inviscd(nb,nghnode,nghedge, &
                          npvt,mb_pvt,neqt,mb_rhst, &
                          sub_limit,sub_intnon,sub_intplt, &
                          tur_flux_steger,tur_flux_euler,dn_via_edge6, &
                          chkid,chklim)
    case(nderive_ehen4)
        call calc_inviscd(nb,nghnode,nghedge, &
                          npvt,mb_pvt,neqt,mb_rhst, &
                          sub_limit,sub_intnon,sub_intplt, &
                          tur_flux_steger,tur_flux_euler,dn_via_ehen4, &
                          chkid,chklim)
    case(nderive_ehen6)
        call calc_inviscd(nb,nghnode,nghedge, &
                          npvt,mb_pvt,neqt,mb_rhst, &
                          sub_limit,sub_intnon,sub_intplt, &
                          tur_flux_steger,tur_flux_euler,dn_via_ehen6, &
                          chkid,chklim)
    case(nderive_ehen6e)
        call calc_inviscd(nb,nghnode,nghedge, &
                          npvt,mb_pvt,neqt,mb_rhst, &
                          sub_limit,sub_intnon,sub_intplt, &
                          tur_flux_steger,tur_flux_euler,dn_via_ehen6e, &
                          chkid,chklim)
    case(nderive_ehcs6)
        call calc_inviscd(nb,nghnode,nghedge, &
                          npvt,mb_pvt,neqt,mb_rhst, &
                          sub_limit,sub_intnon,sub_intplt, &
                          tur_flux_steger,tur_flux_euler,dn_via_ehcs6, &
                          chkid,chklim)
    case(nderive_ehcs6e)
        call calc_inviscd(nb,nghnode,nghedge, &
                          npvt,mb_pvt,neqt,mb_rhst, &
                          sub_limit,sub_intnon,sub_intplt, &
                          tur_flux_steger,tur_flux_euler,dn_via_ehcs6e, &
                          chkid,chklim)
    case(nderive_ehen8)
        call calc_inviscd(nb,nghnode,nghedge, &
                          npvt,mb_pvt,neqt,mb_rhst, &
                          sub_limit,sub_intnon,sub_intplt, &
                          tur_flux_steger,tur_flux_euler,dn_via_ehen8, &
                          chkid,chklim)
    case(nderive_ehen8e)
        call calc_inviscd(nb,nghnode,nghedge, &
                          npvt,mb_pvt,neqt,mb_rhst, &
                          sub_limit,sub_intnon,sub_intplt, &
                          tur_flux_steger,tur_flux_euler,dn_via_ehen8e, &
                          chkid,chklim)
    case(nderive_scsl4)
        call calc_inviscd(nb,nghnode,nghedge, &
                          npvt,mb_pvt,neqt,mb_rhst, &
                          sub_limit,sub_intnon,sub_intplt, &
                          tur_flux_steger,tur_flux_euler,dn_via_scsl4, &
                          chkid,chklim)
    case(nderive_scsl4e)
        call calc_inviscd(nb,nghnode,nghedge, &
                          npvt,mb_pvt,neqt,mb_rhst, &
                          sub_limit,sub_intnon,sub_intplt, &
                          tur_flux_steger,tur_flux_euler,dn_via_scsl4e, &
                          chkid,chklim)
    case(nderive_scsl6)
        call calc_inviscd(nb,nghnode,nghedge, &
                          npvt,mb_pvt,neqt,mb_rhst, &
                          sub_limit,sub_intnon,sub_intplt, &
                          tur_flux_steger,tur_flux_euler,dn_via_scsl6, &
                          chkid,chklim)
    case(nderive_scsl6e)
        call calc_inviscd(nb,nghnode,nghedge, &
                          npvt,mb_pvt,neqt,mb_rhst, &
                          sub_limit,sub_intnon,sub_intplt, &
                          tur_flux_steger,tur_flux_euler,dn_via_scsl6e, &
                          chkid,chklim)
    case default
    end select

end subroutine tur_invscd_scheme

subroutine tur_flux_euler(nsp,nep,prim,nt,nx,ny,nz,nsf,nef,f)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : sml_ssf,one,two,half
    use mod_variables, only : gamma
    implicit none
    integer(kind_int), intent(in)  :: nsp,nep
    real(kind_real)  , intent(in)  :: prim(nsp:nep)
    real(kind_real)  , intent(in)  :: nt,nx,ny,nz
    integer(kind_int), intent(in)  :: nsf,nef
    real(kind_real)  , intent(out) :: f(nsf:nef)
    integer(kind_int) :: m
    real(kind_real)   :: ro,vx,vy,vz,qt(nsf:nef)
    real(kind_real)   :: gama,ae,v2,h0,vn,rvn

    ro = prim(nsp)
    vx = prim(nsp+1)
    vy = prim(nsp+2)
    vz = prim(nsp+3)

    do m=nsf,nef
        qt(m) = prim(nsp+4+m-nsf)
    end do

    vn = nx*vx + ny*vy + nz*vz

    rvn = ro*vn
    do m=nsf,nef
        f(m) = rvn*qt(m)
    end do

end subroutine tur_flux_euler

subroutine tur_flux_steger(nsp,nep,priml,primr,nt,nx,ny,nz,nsf,nef,flr,efix)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : half,sml_ssf
    implicit none
    integer(kind_int), intent(in)  :: nsp,nep
    real(kind_real)  , intent(in)  :: priml(nsp:nep),primr(nsp:nep)
    real(kind_real)  , intent(in)  :: nt,nx,ny,nz
    integer(kind_int), intent(in)  :: nsf,nef
    real(kind_real)  , intent(out) :: flr(nsf:nef)
    real(kind_real)  , intent(in)  :: efix
    integer(kind_int)         :: m,n
    real(kind_real)           :: rl,ul,vl,wl,ql(nsf:nef),vnl,rvnl
    real(kind_real)           :: rr,ur,vr,wr,qr(nsf:nef),vnr,rvnr
    real(kind_real)           :: sn,eps,avnl,avnr
    real(kind_real), external :: enfix_harten

    rl = priml(nsp)
    ul = priml(nsp+1)
    vl = priml(nsp+2)
    wl = priml(nsp+3)

    rr = primr(nsp)
    ur = primr(nsp+1)
    vr = primr(nsp+2)
    wr = primr(nsp+3)

    do m=nsf,nef
        n = nsp+4+m-nsf
        ql(m) = priml(n)
        qr(m) = primr(n)
    end do

    vnl = nx*ul + ny*vl + nz*wl
    vnr = nx*ur + ny*vr + nz*wr

    sn  = max(sqrt(nx*nx+ny*ny+nz*nz),sml_ssf)
    eps = efix*sn

    avnl = enfix_harten(vnl,eps)
    avnr = enfix_harten(vnr,eps)

    rvnl = half*rl*(vnl + avnl)
    rvnr = half*rr*(vnr - avnr)
    do m=nsf,nef
        flr(m) = rvnl*ql(m) + rvnr*qr(m)
    end do

end subroutine tur_flux_steger



subroutine tur_gradient
    use mod_kndconsts, only : kind_int
    use mod_constants, only : nintplt_node2,nintplt_node4
    use mod_constants, only : nintplt_node6,nintplt_node6w
    use mod_constants, only : nintplt_node8
    use mod_constants, only : nintplt_scsl4,nintplt_scsl6
    use mod_constants, only : nintplt_node6e,nintplt_node8e
    use mod_constants, only : nintplt_scsl4e,nintplt_scsl6e
    use mod_variables, only : nd3int_tur_vis
    implicit none
    external :: ve_via_node2,ve_via_node4
    external :: ve_via_node6,ve_via_node6w,ve_via_node8
    external :: ve_via_scsl4,ve_via_scsl6
    external :: ve_via_node6e,ve_via_node8e
    external :: ve_via_scsl4e,ve_via_scsl6e

    select case(nd3int_tur_vis)
    case(nintplt_node2)
        call grad_tur_dnvis(ve_via_node2)
    case(nintplt_node4)
        call grad_tur_dnvis(ve_via_node4)
    case(nintplt_node6)
        call grad_tur_dnvis(ve_via_node6)
    case(nintplt_node6w)
        call grad_tur_dnvis(ve_via_node6w)
    case(nintplt_node6e)
        call grad_tur_dnvis(ve_via_node6e)
    case(nintplt_node8)
        call grad_tur_dnvis(ve_via_node8)
    case(nintplt_node8e)
        call grad_tur_dnvis(ve_via_node8e)
    case(nintplt_scsl4)
        call grad_tur_dnvis(ve_via_scsl4)
    case(nintplt_scsl4e)
        call grad_tur_dnvis(ve_via_scsl4e)
    case(nintplt_scsl6)
        call grad_tur_dnvis(ve_via_scsl6)
    case(nintplt_scsl6e)
        call grad_tur_dnvis(ve_via_scsl6e)
    case default

    end select

end subroutine tur_gradient

subroutine grad_tur_dnvis(sub_intvis)
    use mod_constants, only : nderive_edge6,nderive_edge4,nderive_edge2
    use mod_constants, only : nderive_ehen6,nderive_ehen4,nderive_ehen8,nderive_ehcs6
    use mod_constants, only : nderive_scsl4,nderive_scsl6
    use mod_constants, only : nderive_ehen6e,nderive_ehen8e,nderive_ehcs6e
    use mod_constants, only : nderive_scsl4e,nderive_scsl6e
    use mod_variables, only : nd3der_tur_vis,nghnode,nghedge
    implicit none
    external :: sub_intvis
    external :: dn_via_edge6,dn_via_edge4,dn_via_edge2
    external :: dn_via_ehen6,dn_via_ehen4,dn_via_ehen8,dn_via_ehcs6
    external :: dn_via_scsl4,dn_via_scsl6
    external :: dn_via_ehen6e,dn_via_ehen8e,dn_via_ehcs6e
    external :: dn_via_scsl4e,dn_via_scsl6e

    select case(nd3der_tur_vis)
    case(nderive_edge2)
        call calc_grad_tur(nghnode,nghedge,sub_intvis,dn_via_edge2)
    case(nderive_edge4)
        call calc_grad_tur(nghnode,nghedge,sub_intvis,dn_via_edge4)
    case(nderive_edge6)
        call calc_grad_tur(nghnode,nghedge,sub_intvis,dn_via_edge6)
    case(nderive_ehen4)
        call calc_grad_tur(nghnode,nghedge,sub_intvis,dn_via_ehen4)
    case(nderive_ehen6)
        call calc_grad_tur(nghnode,nghedge,sub_intvis,dn_via_ehen6)
    case(nderive_ehen6e)
        call calc_grad_tur(nghnode,nghedge,sub_intvis,dn_via_ehen6e)
    case(nderive_ehcs6)
        call calc_grad_tur(nghnode,nghedge,sub_intvis,dn_via_ehcs6)
    case(nderive_ehcs6e)
        call calc_grad_tur(nghnode,nghedge,sub_intvis,dn_via_ehcs6e)
    case(nderive_ehen8)
        call calc_grad_tur(nghnode,nghedge,sub_intvis,dn_via_ehen8)
    case(nderive_ehen8e)
        call calc_grad_tur(nghnode,nghedge,sub_intvis,dn_via_ehen8e)
    case(nderive_scsl4)
        call calc_grad_tur(nghnode,nghedge,sub_intvis,dn_via_scsl4)
    case(nderive_scsl4e)
        call calc_grad_tur(nghnode,nghedge,sub_intvis,dn_via_scsl4e)
    case(nderive_scsl6)
        call calc_grad_tur(nghnode,nghedge,sub_intvis,dn_via_scsl6)
    case(nderive_scsl6e)
        call calc_grad_tur(nghnode,nghedge,sub_intvis,dn_via_scsl6e)            
    case default

    end select

end subroutine grad_tur_dnvis


subroutine calc_grad_tur(ngn,nge,sub_intvis,sub_dnvis)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : one,bc_cut1to1
    use mod_constants, only : nfsf_vis_d3int,nfsf_vis_d3der
    use mod_constants, only : nbc_inter_buf_dtur
    use mod_constants, only : nsgl_buffer_dtur,nsgl_aver_art
    use mod_datatypes, only : fld_array_t,top_block_t
    use mod_fieldvars, only : nblkcoms,blkcoms
    use mod_fieldvars, only : mb_sxyz,mb_vol,neqt,mb_qt,mb_dtur
    use mod_interface, only : calc_mb_dn_via_node3
    use mod_interface, only : pre_exchange_bc_var,post_exchange_bc_var,exchange_singulars
    use mod_interface, only : average_bc_var,average_singulars
    implicit none
    integer(kind_int), intent(in) :: ngn,nge
    external                      :: sub_intvis,sub_dnvis
    integer(kind_int)          :: nc,nb,i,j,k,m,m1,m2,m3
    integer(kind_int)          :: st(3),ed(3)
    real(kind_real)            :: kx,ky,kz,ex,ey,ez,cx,cy,cz
    real(kind_real)            :: ovol,der(3*neqt)
    type(top_block_t), pointer :: top
    type(fld_array_t), pointer :: dtur(:),sxyz(:),vol(:)

    do m=1,neqt
        m1 = 3*m - 2
        m2 = m1 + 1
        m3 = m2 + 1
        call calc_mb_dn_via_node3(mb_qt,m,ngn,sub_intvis,nge,sub_dnvis,nfsf_vis_d3int,nfsf_vis_d3der,mb_dtur,m1,1)
        call calc_mb_dn_via_node3(mb_qt,m,ngn,sub_intvis,nge,sub_dnvis,nfsf_vis_d3int,nfsf_vis_d3der,mb_dtur,m2,2)
        call calc_mb_dn_via_node3(mb_qt,m,ngn,sub_intvis,nge,sub_dnvis,nfsf_vis_d3int,nfsf_vis_d3der,mb_dtur,m3,3)
    end do

    do nc=1,nblkcoms
        nb  =  blkcoms(nc)%nb
        top => blkcoms(nc)%top

        sxyz => mb_sxyz(nb)%fld
        vol  => mb_vol(nb)%fld
        dtur => mb_dtur(nb)%fld

        st(:) = 1
        ed(:) = top%nijk(:)

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

            do m=1,3*neqt
               der(m) = dtur(m)%r3d(i,j,k)
            end do

            do m=1,neqt
                m1 = 3*m - 2
                m2 = m1 + 1
                m3 = m2 + 1
          
                dtur(m1)%r3d(i,j,k) = ( kx*der(m1) + ex*der(m2) + cx*der(m3) ) * ovol
                dtur(m2)%r3d(i,j,k) = ( ky*der(m1) + ey*der(m2) + cy*der(m3) ) * ovol
                dtur(m3)%r3d(i,j,k) = ( kz*der(m1) + ez*der(m2) + cz*der(m3) ) * ovol
            end do
        end do
        end do
        end do
    end do

    call pre_exchange_bc_var(mb_dtur,1,3*neqt,ngn,nbc_inter_buf_dtur,nsgl_aver_art)
    call exchange_singulars(mb_dtur,1,3*neqt,nsgl_buffer_dtur,nsgl_aver_art)
    call post_exchange_bc_var(mb_dtur,1,3*neqt,ngn,nbc_inter_buf_dtur,nsgl_aver_art) 
    call average_bc_var(mb_dtur,1,3*neqt,nbc_inter_buf_dtur,nsgl_aver_art)
    call average_singulars(mb_dtur,1,3*neqt,nsgl_buffer_dtur,nsgl_aver_art)

end subroutine calc_grad_tur


subroutine tur_rhside
    use mod_constants, only : nvis_tur_sa
    use mod_constants, only : nvis_tur_sst,nvis_tur_hst
    use mod_variables, only : nvis
    implicit none

    select case(nvis)
    case(nvis_tur_sa)
       call spalart_rhside
    case(nvis_tur_sst,nvis_tur_hst)
       call menter_rhside
    end select

end subroutine tur_rhside


subroutine tur_update
    use mod_constants, only : nvis_tur_sa
    use mod_constants, only : nvis_tur_sst,nvis_tur_hst
    use mod_variables, only : nvis
    implicit none
    external :: spalart_update,menter_update

    select case(nvis)
    case(nvis_tur_sa)
       call run_on_blkcoms(spalart_update)
    case(nvis_tur_sst,nvis_tur_hst)
       call run_on_blkcoms(menter_update)
    end select

end subroutine tur_update
