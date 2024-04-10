
subroutine spalart_rhside
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,nvis_euler,nsgl_buffer_dqt
    use mod_constants, only : nsgl_aver_vol,nbc_inter_buf_dqt
    use mod_variables, only : nghnode
    use mod_fieldvars, only : nblkcoms,blkcoms,neqt,mb_rhst,mb_rtur
    use mod_interface, only : pre_exchange_bc_var,post_exchange_bc_var,average_bc_var
    use mod_interface, only : exchange_singulars,average_singulars
    use mod_interface, only : assign_mb_var_uniform,assign_bc_var_uniform
    implicit none
    integer(kind_int) :: nc,nb
    real(kind_real)   :: rhs0(1:neqt)

    rhs0(:) = zero
    call assign_mb_var_uniform(mb_rhst,1,neqt,nghnode,rhs0)
    call assign_mb_var_uniform(mb_rtur,1,neqt,nghnode,rhs0)

    do nc=1,nblkcoms
        nb = blkcoms(nc)%nb

        call tur_invscd(nb)
    end do


    call tur_gradient


    do nc=1,nblkcoms
        nb = blkcoms(nc)%nb

        call spalart_viscous(nb)

        call spalart_source(nb)

        call spalart_spectrad(nb)
    end do

    call pre_exchange_bc_var(mb_rhst,1,neqt,1,nbc_inter_buf_dqt,nsgl_aver_vol)
    call exchange_singulars(mb_rhst,1,neqt,nsgl_buffer_dqt,nsgl_aver_vol)
    call post_exchange_bc_var(mb_rhst,1,neqt,1,nbc_inter_buf_dqt,nsgl_aver_vol) 
    call average_bc_var(mb_rhst,1,neqt,nbc_inter_buf_dqt,nsgl_aver_vol)
    call average_singulars(mb_rhst,1,neqt,nsgl_buffer_dqt,nsgl_aver_vol)

    call assign_bc_var_uniform(mb_rhst,1,neqt,nghnode,rhs0)

end subroutine spalart_rhside

subroutine spalart_viscous(nb)
    use mod_kndconsts, only : kind_int
    use mod_constants, only : nintplt_node2,nintplt_node4
    use mod_constants, only : nintplt_node6,nintplt_node6w
    use mod_constants, only : nintplt_node8
    use mod_constants, only : nintplt_scsl4,nintplt_scsl6
    use mod_constants, only : nintplt_node6e,nintplt_node8e
    use mod_constants, only : nintplt_scsl4e,nintplt_scsl6e
    use mod_variables, only : nd2int_tur_vis
    implicit none
    integer(kind_int), intent(in) :: nb
    external :: ve_via_node2,ve_via_node4
    external :: ve_via_node6,ve_via_node6w,ve_via_node8
    external :: ve_via_scsl4,ve_via_scsl6
    external :: ve_via_node6e,ve_via_node8e
    external :: ve_via_scsl4e,ve_via_scsl6e

    select case(nd2int_tur_vis)
    case(nintplt_node2)
        call spalart_viscous_scheme(nb,ve_via_node2)
    case(nintplt_node4)
        call spalart_viscous_scheme(nb,ve_via_node4)
    case(nintplt_node6,nintplt_node6w)
        call spalart_viscous_scheme(nb,ve_via_node6)
    case(nintplt_node6e)
        call spalart_viscous_scheme(nb,ve_via_node6e)
    case(nintplt_node8)
        call spalart_viscous_scheme(nb,ve_via_node8)
    case(nintplt_node8e)
        call spalart_viscous_scheme(nb,ve_via_node8e)
    case(nintplt_scsl4)
        call spalart_viscous_scheme(nb,ve_via_scsl4)
    case(nintplt_scsl4e)
        call spalart_viscous_scheme(nb,ve_via_scsl4e)
    case(nintplt_scsl6)
        call spalart_viscous_scheme(nb,ve_via_scsl6)
    case(nintplt_scsl6e)
        call spalart_viscous_scheme(nb,ve_via_scsl6e)
    case default

    end select

end subroutine spalart_viscous


subroutine spalart_viscous_scheme(nb,sub_intplt)
    use mod_kndconsts, only : kind_int
    use mod_constants, only : nderive_edge6,nderive_edge4,nderive_edge2
    use mod_constants, only : nderive_ehen6,nderive_ehen4,nderive_ehen8,nderive_ehcs6
    use mod_constants, only : nderive_scsl4,nderive_scsl6
    use mod_constants, only : nderive_ehen6e,nderive_ehen8e,nderive_ehcs6e
    use mod_constants, only : nderive_scsl4e,nderive_scsl6e
    use mod_variables, only : nd2der_tur_vis,nghnode,nghedge
    implicit none
    integer(kind_int), intent(in) :: nb
    external                      :: sub_intplt
    external :: dn_via_edge6,dn_via_edge4,dn_via_edge2
    external :: dn_via_ehen6,dn_via_ehen4,dn_via_ehen8,dn_via_ehcs6
    external :: dn_via_scsl4,dn_via_scsl6
    external :: dn_via_ehen6e,dn_via_ehen8e,dn_via_ehcs6e
    external :: dn_via_scsl4e,dn_via_scsl6e

    select case(nd2der_tur_vis)
    case(nderive_edge2)
       call calc_spalart_viscous(nb,nghnode,nghedge,sub_intplt,dn_via_edge2)
    case(nderive_edge4)
       call calc_spalart_viscous(nb,nghnode,nghedge,sub_intplt,dn_via_edge4)
    case(nderive_edge6)
       call calc_spalart_viscous(nb,nghnode,nghedge,sub_intplt,dn_via_edge6)
    case(nderive_ehen4)
       call calc_spalart_viscous(nb,nghnode,nghedge,sub_intplt,dn_via_ehen4)
    case(nderive_ehen6)
       call calc_spalart_viscous(nb,nghnode,nghedge,sub_intplt,dn_via_ehen6)
    case(nderive_ehen6e)
       call calc_spalart_viscous(nb,nghnode,nghedge,sub_intplt,dn_via_ehen6e)
    case(nderive_ehcs6)
       call calc_spalart_viscous(nb,nghnode,nghedge,sub_intplt,dn_via_ehcs6)
    case(nderive_ehcs6e)
       call calc_spalart_viscous(nb,nghnode,nghedge,sub_intplt,dn_via_ehcs6e)
    case(nderive_ehen8)
       call calc_spalart_viscous(nb,nghnode,nghedge,sub_intplt,dn_via_ehen8)
    case(nderive_ehen8e)
       call calc_spalart_viscous(nb,nghnode,nghedge,sub_intplt,dn_via_ehen8e)
    case(nderive_scsl4)
       call calc_spalart_viscous(nb,nghnode,nghedge,sub_intplt,dn_via_scsl4)
    case(nderive_scsl4e)
       call calc_spalart_viscous(nb,nghnode,nghedge,sub_intplt,dn_via_scsl4e)
    case(nderive_scsl6)
       call calc_spalart_viscous(nb,nghnode,nghedge,sub_intplt,dn_via_scsl6)
    case(nderive_scsl6e)
       call calc_spalart_viscous(nb,nghnode,nghedge,sub_intplt,dn_via_scsl6e)
    case default
    end select

end subroutine spalart_viscous_scheme

subroutine calc_spalart_viscous(nb,ngn,nge,sub_intplt,sub_scheme)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,one
    use mod_constants, only : nbc_inter_scheme,nsw_dir_close
    use mod_constants, only : nfsf_vis_d2int,nfsf_vis_d2der
    use mod_datatypes, only : fld_array_t
    use mod_variables, only : nsw_kdir,reue,turconsts
    use mod_fieldvars, only : mb_top,mb_sxyz,neqt,mb_rhst
    use mod_fieldvars, only : mb_vsl,mb_pvt,mb_fsf,mb_dtur
    implicit none
    integer(kind_int), intent(in) :: nb,ngn,nge
    external                      :: sub_intplt,sub_scheme
    integer(kind_int)          :: i,j,k,m,ierr
    integer(kind_int)          :: ni,nj,nk,stn,edn,ste,ede
    integer(kind_int)          :: stn0,edn0,ste0,ede0,nfs,nfe
    integer(kind_int)          :: nkst,nked,njed
    real(kind_real)            :: nt,kx,ky,kz,ex,ey,ez,cx,cy,cz
    real(kind_real)            :: sig,cb2,osig,re,ro,nu,visl
    real(kind_real)            :: visb,visu,f(1:neqt),der(3)
    real(kind_real), pointer   :: vn(:,:),ve(:,:)
    real(kind_real), pointer   :: fn(:,:),fc(:,:),dn(:,:)
    type(fld_array_t), pointer :: sxyz(:),dtur(:),rhst(:)
    type(fld_array_t), pointer :: vsl(:),pvt(:)
    type(fld_array_t), pointer :: fsfs(:),fsfe(:)

    sig = turconsts(2)
    cb2 = turconsts(4)

    osig = one/sig
    re   = one/reue

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

    sxyz => mb_sxyz(nb)%fld
    dtur => mb_dtur(nb)%fld
    rhst => mb_rhst(nb)%fld
    vsl  => mb_vsl(nb)%fld
    pvt  => mb_pvt(nb)%fld
           
    ! I-direction      
    fsfs => mb_fsf(nb,1)%fld
    fsfe => mb_fsf(nb,2)%fld

    stn = 1 - ngn
    edn = ni + ngn    
    ste = -nge
    ede = ni + nge
    allocate(vn(stn:edn,1:7), stat=ierr)
    allocate(ve(stn:edn,1:7), stat=ierr)
    allocate(fn(stn:edn,1:neqt), stat=ierr)
    allocate(fc(stn:edn,1:neqt), stat=ierr)
    allocate(dn(stn:edn,1:neqt), stat=ierr)
    do k=nkst,nked
    do j=1,nj
        nfs = fsfs(nfsf_vis_d2int)%i3d( 1,j,k)
        nfe = fsfe(nfsf_vis_d2int)%i3d(ni,j,k)
        if (nfs /= nbc_inter_scheme) then
            stn0 = 1
            ste0 = 0
        else
            stn0 = stn
            ste0 = ste
        end if 

        if (nfe /= nbc_inter_scheme) then
            edn0 = ni
            ede0 = ni
        else
            edn0 = edn
            ede0 = ede
        end if 
            
        do i=stn0,edn0
            ro = pvt(1)%r3d(i,j,k)
            nu = pvt(5)%r3d(i,j,k)
            visl = vsl(1)%r3d(i,j,k)
            visb = ro*nu

            visu = visl + visb*osig

            vn(i,1) = sxyz(1)%r3d(i,j,k)
            vn(i,2) = sxyz(2)%r3d(i,j,k)
            vn(i,3) = sxyz(3)%r3d(i,j,k)

            do m=4,6
                vn(i,m) = dtur(m-3)%r3d(i,j,k)
            end do

            vn(i,7) = visu
        end do

        call sub_intplt(ni,1,7,ngn,vn,nfs,nfe,nge,ve)

        do i=ste0,ede0
            kx = ve(i,1)
            ky = ve(i,2)
            kz = ve(i,3)

            do m=1,3
                der(m) = ve(i,m+3)
            end do

            visu = ve(i,7)

            call spalart_flux_vis(nt,kx,ky,kz,visu,1,3,der,1,neqt,f)

            do m=1,neqt
                fc(i,m) = f(m)
            end do
        end do

        do i=stn0,edn0
            kx = vn(i,1)
            ky = vn(i,2)
            kz = vn(i,3)

            do m=1,3
                der(m) = vn(i,m+3)
            end do

            visu = vn(i,7)

            call spalart_flux_vis(nt,kx,ky,kz,visu,1,3,der,1,neqt,f)
                
            do m=1,neqt
                fn(i,m) = f(m)
            end do
        end do

        nfs = fsfs(nfsf_vis_d2der)%i3d( 1,j,k)
        nfe = fsfe(nfsf_vis_d2der)%i3d(ni,j,k)
        call sub_scheme(ni,1,neqt,ngn,fn,nge,fc,nfs,nfe,dn)

        do m=1,neqt
            do i=1,ni
                rhst(m)%r3d(i,j,k) = rhst(m)%r3d(i,j,k) + re*dn(i,m)
            end do
        end do
    end do
    end do
    deallocate(dn, stat=ierr)
    deallocate(fc, stat=ierr)
    deallocate(fn, stat=ierr)
    deallocate(ve, stat=ierr)
    deallocate(vn, stat=ierr)


    ! J-direction      
    fsfs => mb_fsf(nb,3)%fld
    fsfe => mb_fsf(nb,4)%fld

    stn = 1 - ngn
    edn = nj + ngn    
    ste = -nge
    ede = nj + nge
    allocate(vn(stn:edn,1:7), stat=ierr)
    allocate(ve(stn:edn,1:7), stat=ierr)
    allocate(fn(stn:edn,1:neqt), stat=ierr)
    allocate(fc(stn:edn,1:neqt), stat=ierr)
    allocate(dn(stn:edn,1:neqt), stat=ierr)
    do k=nkst,nked
    do i=1,ni
        nfs = fsfs(nfsf_vis_d2int)%i3d(i, 1,k)
        nfe = fsfe(nfsf_vis_d2int)%i3d(i,nj,k)
        if (nfs /= nbc_inter_scheme) then
            stn0 = 1
            ste0 = 0
        else
            stn0 = stn
            ste0 = ste
        end if 

        if (nfe /= nbc_inter_scheme) then
            edn0 = nj
            ede0 = nj
        else
            edn0 = edn
            ede0 = ede
        end if 
            
        do j=stn0,edn0
            ro = pvt(1)%r3d(i,j,k)
            nu = pvt(5)%r3d(i,j,k)
            visl = vsl(1)%r3d(i,j,k)
            visb = ro*nu

            visu = visl + visb*osig

            vn(j,1) = sxyz(4)%r3d(i,j,k)
            vn(j,2) = sxyz(5)%r3d(i,j,k)
            vn(j,3) = sxyz(6)%r3d(i,j,k)

            do m=4,6
                vn(j,m) = dtur(m-3)%r3d(i,j,k)
            end do

            vn(j,7) = visu
        end do

        call sub_intplt(nj,1,7,ngn,vn,nfs,nfe,nge,ve)

        do j=ste0,ede0
            ex = ve(j,1)
            ey = ve(j,2)
            ez = ve(j,3)

            do m=1,3
                der(m) = ve(j,m+3)
            end do

            visu = ve(j,7)

            call spalart_flux_vis(nt,ex,ey,ez,visu,1,3,der,1,neqt,f)
                
            do m=1,neqt
                fc(j,m) = f(m)
            end do
        end do

        do j=stn0,edn0
            ex = vn(j,1)
            ey = vn(j,2)
            ez = vn(j,3)

            do m=1,3
                der(m) = vn(j,m+3)
            end do

            visu = vn(j,7)

            call spalart_flux_vis(nt,ex,ey,ez,visu,1,3,der,1,neqt,f)
                
            do m=1,neqt
                fn(j,m) = f(m)
            end do
        end do

        nfs = fsfs(nfsf_vis_d2der)%i3d(i, 1,k)
        nfe = fsfe(nfsf_vis_d2der)%i3d(i,nj,k)
        call sub_scheme(nj,1,neqt,ngn,fn,nge,fc,nfs,nfe,dn)
            
        do m=1,neqt
            do j=1,nj
                rhst(m)%r3d(i,j,k) = rhst(m)%r3d(i,j,k) + re*dn(j,m)
            end do
        end do
    end do
    end do
    deallocate(dn, stat=ierr)
    deallocate(fc, stat=ierr)
    deallocate(fn, stat=ierr)
    deallocate(ve, stat=ierr)
    deallocate(vn, stat=ierr)


    ! K-direction      
    fsfs => mb_fsf(nb,5)%fld
    fsfe => mb_fsf(nb,6)%fld

    stn = 1 - ngn
    edn = nk + ngn    
    ste = -nge
    ede = nk + nge
    allocate(vn(stn:edn,1:7), stat=ierr)
    allocate(ve(stn:edn,1:7), stat=ierr)
    allocate(fn(stn:edn,1:neqt), stat=ierr)
    allocate(fc(stn:edn,1:neqt), stat=ierr)
    allocate(dn(stn:edn,1:neqt), stat=ierr)
    do j=1,njed
    do i=1,ni
        nfs = fsfs(nfsf_vis_d2int)%i3d(i,j, 1)
        nfe = fsfe(nfsf_vis_d2int)%i3d(i,j,nk)
        if (nfs /= nbc_inter_scheme) then
            stn0 = 1
            ste0 = 0
        else
            stn0 = stn
            ste0 = ste
        end if 

        if (nfe /= nbc_inter_scheme) then
            edn0 = nk
            ede0 = nk
        else
            edn0 = edn
            ede0 = ede
        end if 
            
        do k=stn0,edn0
            ro = pvt(1)%r3d(i,j,k)
            nu = pvt(5)%r3d(i,j,k)
            visl = vsl(1)%r3d(i,j,k)
            visb = ro*nu

            visu = visl + visb*osig

            vn(k,1) = sxyz(7)%r3d(i,j,k)
            vn(k,2) = sxyz(8)%r3d(i,j,k)
            vn(k,3) = sxyz(9)%r3d(i,j,k)

            do m=4,6
                vn(k,m) = dtur(m-3)%r3d(i,j,k)
            end do

            vn(k,7) = visu
        end do

        call sub_intplt(nk,1,7,ngn,vn,nfs,nfe,nge,ve)

        do k=ste0,ede0
            cx = ve(k,1)
            cy = ve(k,2)
            cz = ve(k,3)

            do m=1,3
                der(m) = ve(k,m+3)
            end do

            visu = ve(k,7)

            call spalart_flux_vis(nt,cx,cy,cz,visu,1,3,der,1,neqt,f)
                
            do m=1,neqt
                fc(k,m) = f(m)
            end do
        end do

        do k=stn0,edn0
            cx = vn(k,1)
            cy = vn(k,2)
            cz = vn(k,3)

            do m=1,3
                der(m) = vn(k,m+3)
            end do

            visu = vn(k,7)

            call spalart_flux_vis(nt,cx,cy,cz,visu,1,3,der,1,neqt,f)
                
            do m=1,neqt
                fn(k,m) = f(m)
            end do
        end do

        nfs = fsfs(nfsf_vis_d2der)%i3d(i,j, 1)
        nfe = fsfe(nfsf_vis_d2der)%i3d(i,j,nk)
        call sub_scheme(nk,1,neqt,ngn,fn,nge,fc,nfs,nfe,dn)
            
        do m=1,neqt
            do k=1,nk
                rhst(m)%r3d(i,j,k) = rhst(m)%r3d(i,j,k) + re*dn(k,m)
            end do
        end do
    end do
    end do
    deallocate(dn, stat=ierr)
    deallocate(fc, stat=ierr)
    deallocate(fn, stat=ierr)
    deallocate(ve, stat=ierr)
    deallocate(vn, stat=ierr)

    if (nsw_kdir == nsw_dir_close) then
        do m=1,neqt
            do k=1,nk
            do j=1,nj
            do i=1,ni
                rhst(m)%r3d(i,j,k) = rhst(m)%r3d(i,j,nkst)
            end do
            end do
            end do
        end do
    end if

end subroutine calc_spalart_viscous

subroutine spalart_flux_vis(nt,nx,ny,nz,visu,nsd,ned,der,nsf,nef,f)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,one,two,two3rd
    use mod_variables, only : gamma,turconsts
    implicit none
    real(kind_real)  , intent(in)  :: nt,nx,ny,nz
    real(kind_real)  , intent(in)  :: visu
    integer(kind_int), intent(in)  :: nsd,ned
    real(kind_real)  , intent(in)  :: der(nsd:ned)
    integer(kind_int), intent(in)  :: nsf,nef
    real(kind_real)  , intent(out) :: f(nsf:nef)
    real(kind_real) :: dnux,dnuy,dnuz,dnun
    
    dnux = der(nsd  )
    dnuy = der(nsd+1)
    dnuz = der(nsd+2)

    dnun = nx*dnux + ny*dnuy + nz*dnuz
                
    f(nsf) = visu*dnun

end subroutine spalart_flux_vis

subroutine spalart_source(nb)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,small,half,third,two3rd,sixth
    use mod_constants, only : one,two,three,six,ten,nsw_dir_close
    use mod_datatypes, only : fld_array_t
    use mod_variables, only : reue,turconsts,nsw_kdir
    use mod_fieldvars, only : mb_top,mb_pv,mb_dpv,mb_vsl
    use mod_fieldvars, only : neqt,mb_qt,mb_rhst,mb_dtur
    use mod_fieldvars, only : mb_dst,mb_rtur,mb_vol
    implicit none
    integer(kind_int), intent(in) :: nb
    integer(kind_int)          :: i,j,k,m,st(3),ed(3)
    integer(kind_int)          :: ni,nj,nk,nkst
    real(kind_real)            :: ro,nu,visl,visb,dist,dist2,vol0,re
    real(kind_real)            :: kap,cb1,cb2,cv1,cv2,cv5,cw2,cw3,sig,ct3,ct4
    real(kind_real)            :: osig,sdilim,kap2,cw1,cv13,cw36,od2,oka2,oky2
    real(kind_real)            :: dux,duy,duz,dvx,dvy,dvz,dwx,dwy,dwz
    real(kind_real)            :: s11,s22,s33,s12,s13,s23,w12,w13,w23
    real(kind_real)            :: divv,sij2,wij2,divvp3,sij2s,str
    real(kind_real)            :: sp,sd,ft2,ftrans,dnux,dnuy,dnuz,grd2
    real(kind_real)            :: ld,ld2,ld3,rr,rr5,gg,gg6,fv1,fv2,fv3,std,ost
    real(kind_real)            :: fw0,fw,dfv1,dfv2,dfv3,dsdu,dfwdg,dgdr,drdu,dfwdu
    real(kind_real)            :: dft2du,srct,sdi,spprim,sdprim,srctpri,fsnu
    type(fld_array_t), pointer :: pv(:),qt(:),vsl(:),dpv(:),vol(:)
    type(fld_array_t), pointer :: dtur(:),rtur(:),dst(:),rhst(:)

    kap = turconsts(1)
    sig = turconsts(2)

    cb1 = turconsts(3)
    cb2 = turconsts(4)

    cv1 = turconsts(5)
    cv2 = turconsts(6)
    cv5 = turconsts(7)

    cw1 = turconsts(8) 
    cw2 = turconsts(9) 
    cw3 = turconsts(10)
       
    sdilim= turconsts(11)

    kap2 = kap*kap
    osig = one/sig
    cv13  = cv1**3
    cw36  = cw3**6
    oka2  = one/kap2

    ct3 = 1.1
    ct4 = 2.0

    re = one/reue

    pv  => mb_pv(nb)%fld
    qt  => mb_qt(nb)%fld
    vsl => mb_vsl(nb)%fld
    dpv => mb_dpv(nb)%fld
    vol => mb_vol(nb)%fld
    dst => mb_dst(nb)%fld
    dtur => mb_dtur(nb)%fld
    rtur => mb_rtur(nb)%fld
    rhst => mb_rhst(nb)%fld

    ni = mb_top(nb)%nijk(1)
    nj = mb_top(nb)%nijk(2)
    nk = mb_top(nb)%nijk(3)

    st(:) = mb_top(nb)%ndst(:)
    ed(:) = mb_top(nb)%nded(:)

    do k=st(3),ed(3)
    do j=st(2),ed(2)
    do i=st(1),ed(1)
        ro = pv(1)%r3d(i,j,k)
        nu = qt(1)%r3d(i,j,k)
        vol0 = vol(1)%r3d(i,j,k)
        dist = dst(1)%r3d(i,j,k)
        visl = vsl(1)%r3d(i,j,k)

        visb = ro*nu

        dist2 = dist*dist
        od2  = one/(dist2+small)
        oky2 = oka2*od2

        dux = dpv(1)%r3d(i,j,k)
        duy = dpv(2)%r3d(i,j,k)
        duz = dpv(3)%r3d(i,j,k)
        dvx = dpv(4)%r3d(i,j,k)
        dvy = dpv(5)%r3d(i,j,k)
        dvz = dpv(6)%r3d(i,j,k)
        dwx = dpv(7)%r3d(i,j,k)
        dwy = dpv(8)%r3d(i,j,k)
        dwz = dpv(9)%r3d(i,j,k)

        s11 = dux
        s22 = dvy
        s33 = dwz
        s12 = half*(duy + dvx)
        s13 = half*(duz + dwx)
        s23 = half*(dvz + dwy)
        w12 = half*(duy - dvx)
        w13 = half*(duz - dwx)
        w23 = half*(dvz - dwy)
        divv = dux + dvy + dwz
        sij2 = s11*s11 + s22*s22 + s33*s33 + &
               two*(s12*s12 + s13*s13 + s23*s23)
        wij2 = two*(w12*w12 + w13*w13 + w23*w23)

        divvp3 = third*divv
        sij2s  = sij2 - divv*divvp3

        str = sqrt(two*wij2)

        ld  = visb/visl
        ld  = max(ld,1.0e-6)
        ld2 = ld*ld
        ld3 = ld*ld2
        fv1 = ld3/(ld3 + cv13)
        !Spalart's modifications
        !new fv2 which avoids negative production near the wall
        fv2 = (cv2/(cv2 + ld))**3
        fv3 = (one/ld+fv1)*(one-fv2)
        std = fv3*str + nu*fv2*oky2*re
        std = max(std,small)  !!max(std,str)
        ost = sign(one,std)/(abs(std) + small)

        rr  = nu*ost*oky2*re
        rr  = min(ten,rr)
        rr5 = rr**5
        gg  = rr*(one + cw2*(rr5 - one))
        gg  = max(gg,small)
        gg6 = gg**6
        fw0 = ((one + cw36)/(gg6 + cw36))**sixth
        fw  = gg*fw0

        !Spalart's modifications
        dfv1 = three*ld2*(one-fv1)/(cv13 + ld3)
        dfv2 = -three*fv2/(cv2+ld)
        dfv3 = -(one/ld+fv1)*dfv2 + (dfv1-one/ld2)*(one-fv2)
        dsdu = dfv3*str/visl + (fv2+ld*dfv2)*oky2*re/ro

        dfwdg = fw0*(one - gg6/(gg6+cw36))
        dgdr  = one + cw2*(six*rr5 - one)
        drdu  = rr/visb - rr*ost*dsdu
        dfwdu = dfwdg*dgdr*drdu
        
        ftrans = zero
        ft2    = ftrans*ct3*exp(-ct4*ld2)
        dft2du = -ftrans*two*ct3*ct4*ld*exp(-ct4*ld2)/visl

        sp = cb1*(one - ft2)*std
        sd = (cw1*fw-cb1*oka2*ft2)*nu*od2*re

        dnux = dtur(1)%r3d(i,j,k)
        dnuy = dtur(2)%r3d(i,j,k)
        dnuz = dtur(3)%r3d(i,j,k)
        grd2 = dnux*dnux + dnuy*dnuy + dnuz*dnuz
        sdi  = cb2*osig*re*ro*grd2

        sdi = min(sdilim*sp*visb, sdi)
                    
        srct = sp - sd + sdi/visb

        rhst(1)%r3d(i,j,k) = rhst(1)%r3d(i,j,k) + srct*visb*vol0

        spprim  = cb1*((one-ft2)*dsdu - dft2du*std)
        sdprim  = sd/visb + (cw1*dfwdu-cb1*oka2*dft2du)*nu*od2*re
        srctpri = spprim - sdprim
        fsnu = half*(srct+abs(srct) + (srctpri+abs(srctpri))*visb)

        rtur(1)%r3d(i,j,k) = rtur(1)%r3d(i,j,k) + fsnu*vol0
    end do
    end do
    end do

    if (nsw_kdir == nsw_dir_close) then
        nkst = mb_top(nb)%ndst(3)

        do m=1,neqt
            do k=1,nk
            do j=1,nj
            do i=1,ni
                rhst(m)%r3d(i,j,k) = rhst(m)%r3d(i,j,nkst)
                rtur(m)%r3d(i,j,k) = rtur(m)%r3d(i,j,nkst)
            end do
            end do
            end do
        end do
    end if

end subroutine spalart_source

subroutine spalart_spectrad(nb)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : one,two
    use mod_constants, only : bc_cut1to1
    use mod_datatypes, only : fld_array_t,bc_region_t
    use mod_variables, only : reue,csrvis,nghnode
    use mod_variables, only : fsw_kdir,turconsts
    use mod_fieldvars, only : nblkcoms,mb_top
    use mod_fieldvars, only : mb_sxyz,mb_vol
    use mod_fieldvars, only : mb_pvt,mb_vsl
    use mod_fieldvars, only : mb_rtur,mb_srt,mb_dt
    use mod_interface, only : ghost_bc_var_nb
    implicit none
    integer(kind_int), intent(in) :: nb
    integer(kind_int)          :: i,j,k,nr
    integer(kind_int)          :: s_st(3),s_ed(3),st(3),ed(3)
    integer(kind_int)          :: nregs,s_nd,s_lr,bctype
    real(kind_real)            :: sig,osig,ro,nu,visl,visb
    real(kind_real)            :: visu,rov,coe,rsum,rva,rvb,rvc
    real(kind_real)            :: vx,vy,vz,kx,ky,kz
    real(kind_real)            :: ex,ey,ez,cx,cy,cz
    real(kind_real)            :: vnk,vne,vnc,snk,sne,snc
    type(bc_region_t), pointer :: reg
    type(fld_array_t), pointer :: pvt(:),sxyz(:),vol(:)
    type(fld_array_t), pointer :: vsl(:),srt(:),rtur(:),dt(:)

    sig  = turconsts(2)
    osig = one/sig

    pvt  => mb_pvt(nb)%fld
    sxyz => mb_sxyz(nb)%fld
    vol  => mb_vol(nb)%fld
    rtur => mb_rtur(nb)%fld
    vsl  => mb_vsl(nb)%fld
    srt  => mb_srt(nb)%fld
    dt   => mb_dt(nb)%fld

    st(:) = 1
    ed(:) = mb_top(nb)%nijk(:)

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

        vx = pvt(2)%r3d(i,j,k)
        vy = pvt(3)%r3d(i,j,k)
        vz = pvt(4)%r3d(i,j,k)
            
        vnk = kx*vx + ky*vy + kz*vz
        vne = ex*vx + ey*vy + ez*vz
        vnc = cx*vx + cy*vy + cz*vz

        snk = kx*kx + ky*ky + kz*kz
        sne = ex*ex + ey*ey + ez*ez
        snc = cx*cx + cy*cy + cz*cz

        rsum = abs(vnk) + abs(vne) + abs(vnc)

        ro = pvt(1)%r3d(i,j,k)
        nu = pvt(5)%r3d(i,j,k)
        visl = vsl(1)%r3d(i,j,k)
        visb = ro*nu

        visu = visl + visb*osig

        rov = ro*vol(1)%r3d(i,j,k)
        coe = csrvis*two*visu/(reue*rov)
        
        rva = snk*coe
        rvb = sne*coe
        rvc = snc*coe
        srt(1)%r3d(i,j,k) = rva
        srt(2)%r3d(i,j,k) = rvb
        srt(3)%r3d(i,j,k) = rvc*fsw_kdir

        rsum = rsum + rva + rvb + rvc

        rsum = rsum + one/dt(1)%r3d(i,j,k)

        rtur(1)%r3d(i,j,k) = rtur(1)%r3d(i,j,k) + rsum
    end do
    end do
    end do

    nregs = mb_top(nb)%nregions
    do nr=1,nregs
        reg => mb_top(nb)%bcs(nr)
                        
        bctype  = reg%bctype
        if (bctype == bc_cut1to1) then
            s_st(:) = reg%s_st(:)
            s_ed(:) = reg%s_ed(:)
            s_nd    = reg%s_nd
            s_lr    = reg%s_lr       
             
            call bc_extend_ghost(s_st,s_ed,s_nd,s_lr,1,nghnode,st,ed)

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

                snk = kx*kx + ky*ky + kz*kz
                sne = ex*ex + ey*ey + ez*ez
                snc = cx*cx + cy*cy + cz*cz
           
                ro = pvt(1)%r3d(i,j,k)
                nu = pvt(5)%r3d(i,j,k)
                visl = vsl(1)%r3d(i,j,k)
                visb = ro*nu

                visu = visl + visb*osig

                rov = ro*vol(1)%r3d(i,j,k)
                coe = csrvis*two*visu/(reue*rov)
        
                rva = snk*coe
                rvb = sne*coe
                rvc = snc*coe
                srt(1)%r3d(i,j,k) = rva
                srt(2)%r3d(i,j,k) = rvb
                srt(3)%r3d(i,j,k) = rvc*fsw_kdir
            end do
            end do
            end do
        end if
    end do

    call ghost_bc_var_nb(nb,mb_srt,1,3)

end subroutine spalart_spectrad

subroutine spalart_update(nb)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : bc_cut1to1,nsw_dir_close
    use mod_datatypes, only : fld_array_t
    use mod_variables, only : reue,turconsts,ndim,nghnode
    use mod_variables, only : kmaxtur,vmaxtur,nsw_kdir
    use mod_fieldvars, only : mb_top,mb_pv,mb_vsl,mb_vst
    use mod_fieldvars, only : neqt,mb_qt,mb_dqt,nqvst,mb_qvst
    implicit none
    integer(kind_int), intent(in) :: nb
    integer(kind_int)          :: i,j,k,ni,nj,nk,m,st(3),ed(3),nkst
    real(kind_real)            :: ro,nu,visb,visl,dnu,nulim
    real(kind_real)            :: cv1,cv13,ld,ld3,fv1,vist
    type(fld_array_t), pointer :: pv(:),qt(:),dqt(:)
    type(fld_array_t), pointer :: vsl(:),vst(:),qvst(:)

    cv1   = turconsts(5)
    nulim = turconsts(12)

    cv13 = cv1**3

    pv  => mb_pv(nb)%fld
    qt  => mb_qt(nb)%fld
    dqt => mb_dqt(nb)%fld
    vsl => mb_vsl(nb)%fld
    vst => mb_vst(nb)%fld
    qvst => mb_qvst(nb)%fld

    ni = mb_top(nb)%nijk(1)
    nj = mb_top(nb)%nijk(2)
    nk = mb_top(nb)%nijk(3)

    st(:) = mb_top(nb)%ndst(:)
    ed(:) = mb_top(nb)%nded(:)

    do k=st(3),ed(3)
    do j=st(2),ed(2)
    do i=st(1),ed(1)
        ro   = pv(1)%r3d(i,j,k)
        visl = vsl(1)%r3d(i,j,k)
       
        dnu = dqt(1)%r3d(i,j,k)/ro
        nu  = qt(1)%r3d(i,j,k) + dnu
        nu  = max(nu,nulim)
        qt(1)%r3d(i,j,k) = nu
        
        visb = ro*nu
        ld   = visb/visl
        ld3  = ld**3
        fv1  = ld3/(ld3 + cv13)
        vist = fv1*visb
          
        vist = min(vist,vmaxtur*visl)
        vst(1)%r3d(i,j,k) = vist
    end do
    end do
    end do

    if (nsw_kdir == nsw_dir_close) then
        nkst = mb_top(nb)%ndst(3)

        do m=1,nqvst
            do k=1,nk
            do j=1,nj
            do i=1,ni
                qvst(m)%r3d(i,j,k) = qvst(m)%r3d(i,j,nkst)
            end do
            end do
            end do
        end do
    end if

end subroutine spalart_update