
subroutine menter_rhside
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

        call menter_blend(nb)

        call menter_viscous(nb)

        call menter_source(nb)

        call menter_spectrad(nb)
    end do

    call pre_exchange_bc_var(mb_rhst,1,neqt,1,nbc_inter_buf_dqt,nsgl_aver_vol)
    call exchange_singulars(mb_rhst,1,neqt,nsgl_buffer_dqt,nsgl_aver_vol)
    call post_exchange_bc_var(mb_rhst,1,neqt,1,nbc_inter_buf_dqt,nsgl_aver_vol)  
    call average_bc_var(mb_rhst,1,neqt,nbc_inter_buf_dqt,nsgl_aver_vol)
    call average_singulars(mb_rhst,1,neqt,nsgl_buffer_dqt,nsgl_aver_vol)

    call assign_bc_var_uniform(mb_rhst,1,neqt,nghnode,rhs0)

end subroutine menter_rhside


subroutine menter_blend(nb)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : half,third,two3rd,one,two
    use mod_constants, only : nvis_tur_sst,nvis_tur_hst
    use mod_constants, only : large,bc_cut1to1
    use mod_datatypes, only : fld_array_t,bc_region_t
    use mod_variables, only : reue,turconsts,nvis,nghnode
    use mod_fieldvars, only : mb_top,mb_pv,mb_vsl
    use mod_fieldvars, only : neqt,mb_qt,mb_dtur
    use mod_fieldvars, only : mb_dst,mb_bld
    implicit none
    integer(kind_int), intent(in) :: nb
    integer(kind_int)          :: i,j,k,nr
    integer(kind_int)          :: s_st(3),s_ed(3),st(3),ed(3)
    integer(kind_int)          :: nregs,s_nd,s_lr,bctype
    real(kind_real)            :: ro,tke,omeg,crss,cdkwmax,cdkwmin
    real(kind_real)            :: cdkw,tmp1,tmp2,tmp3,arg1,fbsl
    real(kind_real)            :: dist,dist2,visl,betas,sigw2,ct3,cmx
    real(kind_real)            :: dkx,dky,dkz,dwx,dwy,dwz
    type(fld_array_t), pointer :: pv(:),qt(:),vsl(:)
    type(fld_array_t), pointer :: dtur(:),dst(:),bld(:)
    type(bc_region_t), pointer :: reg

    betas = turconsts(2) 
    sigw2 = turconsts(8)

    select case(nvis)
    case(nvis_tur_sst)
        ct3 = 4.0
        cmx = 1.0
    case(nvis_tur_hst)
        ct3 = 40.0
        cmx = 1.5
    end select
    
    pv  => mb_pv(nb)%fld
    qt  => mb_qt(nb)%fld
    vsl => mb_vsl(nb)%fld
    dst => mb_dst(nb)%fld
    bld => mb_bld(nb)%fld
    dtur => mb_dtur(nb)%fld

    st(:) = 1
    ed(:) = mb_top(nb)%nijk(:)

    cdkwmax = -large
    do k=st(3),ed(3)
    do j=st(2),ed(2)
    do i=st(1),ed(1)
        ro   = pv(1)%r3d(i,j,k)
        omeg = qt(2)%r3d(i,j,k)

        dkx = dtur(1)%r3d(i,j,k)
        dky = dtur(2)%r3d(i,j,k)
        dkz = dtur(3)%r3d(i,j,k)
        dwx = dtur(4)%r3d(i,j,k)
        dwy = dtur(5)%r3d(i,j,k)
        dwz = dtur(6)%r3d(i,j,k)
        crss = dkx*dwx + dky*dwy + dkz*dwz

        ! calculate cd_kw
        cdkw = two*ro*sigw2*crss/omeg
        cdkwmax = max(cdkwmax,crss)
    end do
    end do
    end do

    !!cdkwmin = 1.0e-20
    cdkwmin = 1.0e-8*cdkwmax

    do k=st(3),ed(3)
    do j=st(2),ed(2)
    do i=st(1),ed(1)
        ro   = pv(1)%r3d(i,j,k)
        tke  = qt(1)%r3d(i,j,k)
        omeg = qt(2)%r3d(i,j,k)
        visl = vsl(1)%r3d(i,j,k)
        dist = dst(1)%r3d(i,j,k)

        dist2 = dist*dist

        dkx = dtur(1)%r3d(i,j,k)
        dky = dtur(2)%r3d(i,j,k)
        dkz = dtur(3)%r3d(i,j,k)
        dwx = dtur(4)%r3d(i,j,k)
        dwy = dtur(5)%r3d(i,j,k)
        dwz = dtur(6)%r3d(i,j,k)
        crss = dkx*dwx + dky*dwy + dkz*dwz

        ! calculate cd_kw
        cdkw = two*ro*sigw2*crss/omeg
        cdkw = max(cdkw,cdkwmin)

        ! calculate arg1
        tmp1 = sqrt(tke)/(betas*dist*omeg)
        tmp2 = 500.0*visl/(ro*omeg*dist2*reue)
        tmp3 = ct3*ro*sigw2*tke/(cdkw*dist2)
        arg1 = min(max(tmp1,tmp2),tmp3)
        arg1 = min(arg1,10.0)

        !Calculate the blending function fbsl
        fbsl = tanh(cmx*arg1**4)
        bld(1)%r3d(i,j,k) = fbsl
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
                ro   = pv(1)%r3d(i,j,k)
                tke  = qt(1)%r3d(i,j,k)
                omeg = qt(2)%r3d(i,j,k)
                visl = vsl(1)%r3d(i,j,k)
                dist = dst(1)%r3d(i,j,k)

                dist2 = dist*dist

                dkx = dtur(1)%r3d(i,j,k)
                dky = dtur(2)%r3d(i,j,k)
                dkz = dtur(3)%r3d(i,j,k)
                dwx = dtur(4)%r3d(i,j,k)
                dwy = dtur(5)%r3d(i,j,k)
                dwz = dtur(6)%r3d(i,j,k)
                crss = dkx*dwx + dky*dwy + dkz*dwz

                ! calculate cd_kw
                cdkw = two*ro*sigw2*crss/omeg
                cdkw = max(cdkw,cdkwmin)

                ! calculate arg1
                tmp1 = sqrt(tke)/(betas*dist*omeg)
                tmp2 = 500.0*visl/(ro*omeg*dist2*reue)
                tmp3 = ct3*ro*sigw2*tke/(cdkw*dist2)
                arg1 = min(max(tmp1,tmp2),tmp3)
                arg1 = min(arg1,10.0)

                !Calculate the blending function fbsl
                fbsl = tanh(cmx*arg1**4)
                bld(1)%r3d(i,j,k) = fbsl
            end do
            end do
            end do
        end if
    end do

end subroutine menter_blend

subroutine menter_viscous(nb)
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
        call menter_viscous_scheme(nb,ve_via_node2)
    case(nintplt_node4)
        call menter_viscous_scheme(nb,ve_via_node4)
    case(nintplt_node6,nintplt_node6w)
        call menter_viscous_scheme(nb,ve_via_node6)
    case(nintplt_node6e)
        call menter_viscous_scheme(nb,ve_via_node6e)
    case(nintplt_node8)
        call menter_viscous_scheme(nb,ve_via_node8)
    case(nintplt_node8e)
        call menter_viscous_scheme(nb,ve_via_node8e)
    case(nintplt_scsl4)
        call menter_viscous_scheme(nb,ve_via_scsl4)
    case(nintplt_scsl4e)
        call menter_viscous_scheme(nb,ve_via_scsl4e)
    case(nintplt_scsl6)
        call menter_viscous_scheme(nb,ve_via_scsl6)
    case(nintplt_scsl6e)
        call menter_viscous_scheme(nb,ve_via_scsl6e)
    case default

    end select

end subroutine menter_viscous


subroutine menter_viscous_scheme(nb,sub_intplt)
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
       call calc_menter_viscous(nb,nghnode,nghedge,sub_intplt,dn_via_edge2)
    case(nderive_edge4)
       call calc_menter_viscous(nb,nghnode,nghedge,sub_intplt,dn_via_edge4)
    case(nderive_edge6)
       call calc_menter_viscous(nb,nghnode,nghedge,sub_intplt,dn_via_edge6)
    case(nderive_ehen4)
       call calc_menter_viscous(nb,nghnode,nghedge,sub_intplt,dn_via_ehen4)
    case(nderive_ehen6)
       call calc_menter_viscous(nb,nghnode,nghedge,sub_intplt,dn_via_ehen6)
    case(nderive_ehen6e)
       call calc_menter_viscous(nb,nghnode,nghedge,sub_intplt,dn_via_ehen6e)
    case(nderive_ehcs6)
       call calc_menter_viscous(nb,nghnode,nghedge,sub_intplt,dn_via_ehcs6)
    case(nderive_ehcs6e)
       call calc_menter_viscous(nb,nghnode,nghedge,sub_intplt,dn_via_ehcs6e)
    case(nderive_ehen8)
       call calc_menter_viscous(nb,nghnode,nghedge,sub_intplt,dn_via_ehen8)
    case(nderive_ehen8e)
       call calc_menter_viscous(nb,nghnode,nghedge,sub_intplt,dn_via_ehen8e)
    case(nderive_scsl4)
       call calc_menter_viscous(nb,nghnode,nghedge,sub_intplt,dn_via_scsl4)
    case(nderive_scsl4e)
       call calc_menter_viscous(nb,nghnode,nghedge,sub_intplt,dn_via_scsl4e)
    case(nderive_scsl6)
       call calc_menter_viscous(nb,nghnode,nghedge,sub_intplt,dn_via_scsl6)
    case(nderive_scsl6e)
       call calc_menter_viscous(nb,nghnode,nghedge,sub_intplt,dn_via_scsl6e)              
    case default
    end select

end subroutine menter_viscous_scheme

subroutine calc_menter_viscous(nb,ngn,nge,sub_intplt,sub_scheme)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,one
    use mod_constants, only : nbc_inter_scheme,nsw_dir_close
    use mod_constants, only : nfsf_vis_d2int,nfsf_vis_d2der
    use mod_datatypes, only : fld_array_t
    use mod_variables, only : nsw_kdir,reue,turconsts
    use mod_fieldvars, only : mb_top,mb_sxyz,neqt,mb_rhst,mb_bld
    use mod_fieldvars, only : mb_vsl,mb_vst,mb_fsf,mb_dtur
    implicit none
    integer(kind_int), intent(in) :: nb,ngn,nge
    external                      :: sub_intplt,sub_scheme
    integer(kind_int)          :: i,j,k,m,ierr
    integer(kind_int)          :: ni,nj,nk,stn,edn,ste,ede
    integer(kind_int)          :: stn0,edn0,ste0,ede0,nfs,nfe
    integer(kind_int)          :: nkst,nked,njed
    real(kind_real)            :: nt,kx,ky,kz,ex,ey,ez,cx,cy,cz
    real(kind_real)            :: re,visl,vist,fbsl,visk,visw
    real(kind_real)            :: sigk1,sigk2,sigk,sigw1,sigw2,sigw
    real(kind_real)            :: f(1:neqt),der(6)
    real(kind_real), pointer   :: vn(:,:),ve(:,:)
    real(kind_real), pointer   :: fn(:,:),fc(:,:),dn(:,:)
    type(fld_array_t), pointer :: sxyz(:),dtur(:),rhst(:)
    type(fld_array_t), pointer :: bld(:),vsl(:),vst(:)
    type(fld_array_t), pointer :: fsfs(:),fsfe(:)

    sigk1 = turconsts(3) 
    sigw1 = turconsts(4) 
    sigk2 = turconsts(7) 
    sigw2 = turconsts(8) 

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

    sxyz => mb_sxyz(nb)%fld
    dtur => mb_dtur(nb)%fld
    rhst => mb_rhst(nb)%fld
    vsl  => mb_vsl(nb)%fld
    vst  => mb_vst(nb)%fld
    bld  => mb_bld(nb)%fld

           
    ! I-direction      
    fsfs => mb_fsf(nb,1)%fld
    fsfe => mb_fsf(nb,2)%fld

    stn = 1 - ngn
    edn = ni + ngn    
    ste = -nge
    ede = ni + nge
    allocate(vn(stn:edn,1:11), stat=ierr)
    allocate(ve(stn:edn,1:11), stat=ierr)
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
            fbsl = bld(1)%r3d(i,j,k)
            sigk = fbsl*sigk1 + (one-fbsl)*sigk2
            sigw = fbsl*sigw1 + (one-fbsl)*sigw2

            visl = vsl(1)%r3d(i,j,k)
            vist = vst(1)%r3d(i,j,k)
            visk = visl + vist*sigk
            visw = visl + vist*sigw

            vn(i,1) = sxyz(1)%r3d(i,j,k)
            vn(i,2) = sxyz(2)%r3d(i,j,k)
            vn(i,3) = sxyz(3)%r3d(i,j,k)

            do m=4,9
                vn(i,m) = dtur(m-3)%r3d(i,j,k)
            end do

            vn(i,10) = visk
            vn(i,11) = visw
        end do

        call sub_intplt(ni,1,11,ngn,vn,nfs,nfe,nge,ve)

        do i=ste0,ede0
            kx = ve(i,1)
            ky = ve(i,2)
            kz = ve(i,3)

            do m=1,6
                der(m) = ve(i,m+3)
            end do

            visk = ve(i,10)
            visw = ve(i,11)

            call menter_flux_vis(nt,kx,ky,kz,visk,visw,1,6,der,1,neqt,f)

            do m=1,neqt
                fc(i,m) = f(m)
            end do
        end do

        do i=stn0,edn0
            kx = vn(i,1)
            ky = vn(i,2)
            kz = vn(i,3)

            do m=1,6
                der(m) = vn(i,m+3)
            end do

            visk = vn(i,10)
            visw = vn(i,11)

            call menter_flux_vis(nt,kx,ky,kz,visk,visw,1,6,der,1,neqt,f)
                
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
    allocate(vn(stn:edn,1:11), stat=ierr)
    allocate(ve(stn:edn,1:11), stat=ierr)
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
            fbsl = bld(1)%r3d(i,j,k)
            sigk = fbsl*sigk1 + (one-fbsl)*sigk2
            sigw = fbsl*sigw1 + (one-fbsl)*sigw2

            visl = vsl(1)%r3d(i,j,k)
            vist = vst(1)%r3d(i,j,k)
            visk = visl + vist*sigk
            visw = visl + vist*sigw

            vn(j,1) = sxyz(4)%r3d(i,j,k)
            vn(j,2) = sxyz(5)%r3d(i,j,k)
            vn(j,3) = sxyz(6)%r3d(i,j,k)

            do m=4,9
                vn(j,m) = dtur(m-3)%r3d(i,j,k)
            end do

            vn(j,10) = visk
            vn(j,11) = visw
        end do

        call sub_intplt(nj,1,11,ngn,vn,nfs,nfe,nge,ve)

        do j=ste0,ede0
            ex = ve(j,1)
            ey = ve(j,2)
            ez = ve(j,3)

            do m=1,6
                der(m) = ve(j,m+3)
            end do

            visk = ve(j,10)
            visw = ve(j,11)

            call menter_flux_vis(nt,ex,ey,ez,visk,visw,1,6,der,1,neqt,f)
                
            do m=1,neqt
                fc(j,m) = f(m)
            end do
        end do

        do j=stn0,edn0
            ex = vn(j,1)
            ey = vn(j,2)
            ez = vn(j,3)

            do m=1,6
                der(m) = vn(j,m+3)
            end do

            visk = vn(j,10)
            visw = vn(j,11)

            call menter_flux_vis(nt,ex,ey,ez,visk,visw,1,6,der,1,neqt,f)
                
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
    allocate(vn(stn:edn,1:11), stat=ierr)
    allocate(ve(stn:edn,1:11), stat=ierr)
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
            fbsl = bld(1)%r3d(i,j,k)
            sigk = fbsl*sigk1 + (one-fbsl)*sigk2
            sigw = fbsl*sigw1 + (one-fbsl)*sigw2

            visl = vsl(1)%r3d(i,j,k)
            vist = vst(1)%r3d(i,j,k)
            visk = visl + vist*sigk
            visw = visl + vist*sigw

            vn(k,1) = sxyz(7)%r3d(i,j,k)
            vn(k,2) = sxyz(8)%r3d(i,j,k)
            vn(k,3) = sxyz(9)%r3d(i,j,k)

            do m=4,9
                vn(k,m) = dtur(m-3)%r3d(i,j,k)
            end do

            vn(k,10) = visk
            vn(k,11) = visw
        end do

        call sub_intplt(nk,1,11,ngn,vn,nfs,nfe,nge,ve)

        do k=ste0,ede0
            cx = ve(k,1)
            cy = ve(k,2)
            cz = ve(k,3)

            do m=1,6
                der(m) = ve(k,m+3)
            end do

            visk = ve(k,10)
            visw = ve(k,11)

            call menter_flux_vis(nt,cx,cy,cz,visk,visw,1,6,der,1,neqt,f)
                
            do m=1,neqt
                fc(k,m) = f(m)
            end do
        end do

        do k=stn0,edn0
            cx = vn(k,1)
            cy = vn(k,2)
            cz = vn(k,3)

            do m=1,6
                der(m) = vn(k,m+3)
            end do

            visk = vn(k,10)
            visw = vn(k,11)

            call menter_flux_vis(nt,cx,cy,cz,visk,visw,1,6,der,1,neqt,f)
                
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

end subroutine calc_menter_viscous

subroutine menter_flux_vis(nt,nx,ny,nz,visk,visw,nsd,ned,der,nsf,nef,f)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,one,two,two3rd
    use mod_variables, only : gamma,turconsts
    implicit none
    real(kind_real)  , intent(in)  :: nt,nx,ny,nz
    real(kind_real)  , intent(in)  :: visk,visw
    integer(kind_int), intent(in)  :: nsd,ned
    real(kind_real)  , intent(in)  :: der(nsd:ned)
    integer(kind_int), intent(in)  :: nsf,nef
    real(kind_real)  , intent(out) :: f(nsf:nef)
    real(kind_real) :: dkx,dky,dkz,dkn
    real(kind_real) :: dwx,dwy,dwz,dwn
    
    dkx = der(nsd  )
    dky = der(nsd+1)
    dkz = der(nsd+2)

    dwx = der(nsd+3)
    dwy = der(nsd+4)
    dwz = der(nsd+5)
                
    dkn = nx*dkx + ny*dky + nz*dkz
    dwn = nx*dwx + ny*dwy + nz*dwz
                
    f(nsf  ) = visk*dkn
    f(nsf+1) = visw*dwn

end subroutine menter_flux_vis

subroutine menter_source(nb)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,half,third,two3rd
    use mod_constants, only : one,two,nsw_dir_close
    use mod_constants, only : nvis_tur_sst,nvis_tur_hst
    use mod_datatypes, only : fld_array_t
    use mod_variables, only : reue,turconsts,nvis,nsw_kdir
    use mod_fieldvars, only : mb_top,mb_pv,mb_dpv,mb_vst
    use mod_fieldvars, only : neqt,mb_qt,mb_rhst,mb_dtur
    use mod_fieldvars, only : mb_dst,mb_rtur,mb_bld,mb_vol
    implicit none
    integer(kind_int), intent(in) :: nb
    integer(kind_int)          :: i,j,k,m,st(3),ed(3)
    integer(kind_int)          :: ni,nj,nk,nkst
    real(kind_real)            :: ro,tke,omeg,vist,vol0,re
    real(kind_real)            :: alfa1,alfa2,beta1,beta2,sigw2
    real(kind_real)            :: dist,crss,fbsl,rok,row,alfa,beta
    real(kind_real)            :: dux,duy,duz,dvx,dvy,dvz,dwx,dwy,dwz
    real(kind_real)            :: s11,s22,s33,s12,s13,s23,w12,w13,w23
    real(kind_real)            :: divv,sij2,wij2,divvp3,sij2s,pklim
    real(kind_real)            :: betas,betascr,beta_cr,prodk,dissk
    real(kind_real)            :: prodw,dissw,srck,srcw,diak,diaw,fskn,fswn
    real(kind_real)            :: sigd,sigd1,sigd2,cdkw,dkx,dky,dkz
    type(fld_array_t), pointer :: pv(:),qt(:),vst(:),dpv(:),vol(:)
    type(fld_array_t), pointer :: dtur(:),rtur(:),dst(:),bld(:),rhst(:)

    betas = turconsts(2)

    beta1 = turconsts(5)
    alfa1 = turconsts(6)

    sigw2 = turconsts(8)

    beta2 = turconsts(9)  
    alfa2 = turconsts(10) 

    sigd1 = turconsts(15) 
    sigd2 = turconsts(16) 

    pklim = turconsts(12)

    re = one/reue

    pv  => mb_pv(nb)%fld
    qt  => mb_qt(nb)%fld
    vst => mb_vst(nb)%fld
    dpv => mb_dpv(nb)%fld
    vol => mb_vol(nb)%fld
    dst => mb_dst(nb)%fld
    bld => mb_bld(nb)%fld
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
        ro   = pv(1)%r3d(i,j,k)
        tke  = qt(1)%r3d(i,j,k)
        omeg = qt(2)%r3d(i,j,k)
        vist = vst(1)%r3d(i,j,k)
        vol0 = vol(1)%r3d(i,j,k)
        dist = dst(1)%r3d(i,j,k)
        fbsl = bld(1)%r3d(i,j,k)

        rok = ro*tke
        row = ro*omeg
          
        alfa = alfa1*fbsl + alfa2*(one - fbsl)
        beta = beta1*fbsl + beta2*(one - fbsl)

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

        betascr = betas
        beta_cr = beta

        prodk = two*vist*sij2s*re - two3rd*divv*rok
        dissk = betascr*rok*omeg
         
        prodw = alfa*reue*ro*prodk/vist
        dissw = beta_cr*row*omeg

        dkx = dtur(1)%r3d(i,j,k)
        dky = dtur(2)%r3d(i,j,k)
        dkz = dtur(3)%r3d(i,j,k)
        dwx = dtur(4)%r3d(i,j,k)
        dwy = dtur(5)%r3d(i,j,k)
        dwz = dtur(6)%r3d(i,j,k)
        crss = dkx*dwx + dky*dwy + dkz*dwz

        select case(nvis)
        case(nvis_tur_sst)
            sigd = (one-fbsl)*two*sigw2
            cdkw = sigd*ro*crss/omeg
        case(nvis_tur_hst)
            sigd = sigd1*fbsl + sigd2*(one - fbsl)
            cdkw = sigd*ro*max(crss,zero)/omeg
        end select

        prodk = min(prodk, pklim*dissk)
          
        srck = prodk - dissk
        srcw = prodw - dissw + cdkw

        rhst(1)%r3d(i,j,k) = rhst(1)%r3d(i,j,k) + srck*vol0
        rhst(2)%r3d(i,j,k) = rhst(2)%r3d(i,j,k) + srcw*vol0

        !The linearization proposed by Menter is used.
        diak = two*betascr*omeg
        diaw = two*beta_cr*omeg + abs(cdkw)/row
        !The linearization proposed by Moryossef.
        !!diak = betascr*rok*reue/vist
        !!diaw = two*alfa*sij2s/omeg + beta_cr*omeg + cdkw/row
          
        fskn = half*( (srck + abs(srck))/rok + (diak + abs(diak)) )
        fswn = half*( (srcw + abs(srcw))/row + (diaw + abs(diaw)) )

        rtur(1)%r3d(i,j,k) = rtur(1)%r3d(i,j,k) + fskn*vol0
        rtur(2)%r3d(i,j,k) = rtur(2)%r3d(i,j,k) + fswn*vol0
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

end subroutine menter_source

subroutine menter_spectrad(nb)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : one,two
    use mod_constants, only : bc_cut1to1
    use mod_datatypes, only : fld_array_t,bc_region_t
    use mod_variables, only : reue,csrvis,nghnode
    use mod_variables, only : fsw_kdir,turconsts
    use mod_fieldvars, only : nblkcoms,mb_top
    use mod_fieldvars, only : mb_pv,mb_sxyz,mb_vol
    use mod_fieldvars, only : mb_vsl,mb_vst,mb_bld
    use mod_fieldvars, only : mb_rtur,mb_srt,mb_dt
    use mod_interface, only : ghost_bc_var_nb
    implicit none
    integer(kind_int), intent(in) :: nb
    integer(kind_int)          :: i,j,k,nr
    integer(kind_int)          :: s_st(3),s_ed(3),st(3),ed(3)
    integer(kind_int)          :: nregs,s_nd,s_lr,bctype
    real(kind_real)            :: sigk1,sigk2,sigw1,sigw2
    real(kind_real)            :: sigk,sigw,fbsl
    real(kind_real)            :: vx,vy,vz,kx,ky,kz
    real(kind_real)            :: ex,ey,ez,cx,cy,cz
    real(kind_real)            :: vnk,vne,vnc,snk,sne,snc
    real(kind_real)            :: visk,visw,visl,vist
    real(kind_real)            :: rov,coe,rsum,rva,rvb,rvc
    type(bc_region_t), pointer :: reg
    type(fld_array_t), pointer :: pv(:),sxyz(:),vol(:)
    type(fld_array_t), pointer :: vsl(:),vst(:),srt(:)
    type(fld_array_t), pointer :: rtur(:),bld(:),dt(:)

    sigk1 = turconsts(3) 
    sigw1 = turconsts(4) 
    sigk2 = turconsts(7) 
    sigw2 = turconsts(8) 

    pv   => mb_pv(nb)%fld
    sxyz => mb_sxyz(nb)%fld
    vol  => mb_vol(nb)%fld
    rtur => mb_rtur(nb)%fld
    vsl  => mb_vsl(nb)%fld
    vst  => mb_vst(nb)%fld
    bld  => mb_bld(nb)%fld
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

        vx = pv(2)%r3d(i,j,k)
        vy = pv(3)%r3d(i,j,k)
        vz = pv(4)%r3d(i,j,k)
            
        vnk = kx*vx + ky*vy + kz*vz
        vne = ex*vx + ey*vy + ez*vz
        vnc = cx*vx + cy*vy + cz*vz

        snk = kx*kx + ky*ky + kz*kz
        sne = ex*ex + ey*ey + ez*ez
        snc = cx*cx + cy*cy + cz*cz

        rsum = abs(vnk) + abs(vne) + abs(vnc)

        fbsl = bld(1)%r3d(i,j,k)
        sigk = fbsl*sigk1 + (one-fbsl)*sigk2
        sigw = fbsl*sigw1 + (one-fbsl)*sigw2

        visl = vsl(1)%r3d(i,j,k)
        vist = vst(1)%r3d(i,j,k)
        visk = visl + vist*sigk
        visw = visl + vist*sigw

        rov = pv(1)%r3d(i,j,k)*vol(1)%r3d(i,j,k)
        coe = csrvis*two*max(visk,visw)/(reue*rov)
        
        rva = snk*coe
        rvb = sne*coe
        rvc = snc*coe
        srt(1)%r3d(i,j,k) = rva
        srt(2)%r3d(i,j,k) = rvb
        srt(3)%r3d(i,j,k) = rvc*fsw_kdir

        rsum = rsum + rva + rvb + rvc

        rsum = rsum + one/dt(1)%r3d(i,j,k)

        rtur(1)%r3d(i,j,k) = rtur(1)%r3d(i,j,k) + rsum
        rtur(2)%r3d(i,j,k) = rtur(2)%r3d(i,j,k) + rsum
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
           
                fbsl = bld(1)%r3d(i,j,k)
                sigk = fbsl*sigk1 + (one-fbsl)*sigk2
                sigw = fbsl*sigw1 + (one-fbsl)*sigw2

                visl = vsl(1)%r3d(i,j,k)
                vist = vst(1)%r3d(i,j,k)
                visk = visl + vist*sigk
                visw = visl + vist*sigw

                rov = pv(1)%r3d(i,j,k)*vol(1)%r3d(i,j,k)
                coe = csrvis*two*max(visk,visw)/(reue*rov)
        
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

end subroutine menter_spectrad

subroutine menter_update(nb)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : half,third,two3rd,one,two
    use mod_constants, only : nvis_tur_sst,nvis_tur_hst
    use mod_constants, only : bc_cut1to1,nsw_dir_close
    use mod_datatypes, only : fld_array_t
    use mod_variables, only : reue,turconsts,nghnode
    use mod_variables, only : kmaxtur,vmaxtur,nvis,nsw_kdir
    use mod_fieldvars, only : mb_top,mb_pv,mb_dpv,mb_vsl,mb_vst
    use mod_fieldvars, only : neqt,mb_qt,mb_dqt,mb_dst,nqvst,mb_qvst
    implicit none
    integer(kind_int), intent(in) :: nb
    integer(kind_int)          :: i,j,k,ni,nj,nk,m,st(3),ed(3),nkst
    real(kind_real)            :: ro,tke,omeg,dk,dw,klim,owlim
    real(kind_real)            :: dux,duy,duz,dvx,dvy,dvz,dwx,dwy,dwz
    real(kind_real)            :: s11,s22,s33,s12,s13,s23,w12,w13,w23
    real(kind_real)            :: divv,sij2,wij2,divvp3,sij2s
    real(kind_real)            :: a1,tmp1,tmp2,arg1,fsst,cmx,vort
    real(kind_real)            :: dist,dist2,visl,vist,betas
    type(fld_array_t), pointer :: qt(:),vsl(:),vst(:),qvst(:)
    type(fld_array_t), pointer :: pv(:),dpv(:),dqt(:),dst(:)

    betas = turconsts(2) 
    a1    = turconsts(11) 
    klim  = turconsts(13)
    owlim = turconsts(14)

    select case(nvis)
    case(nvis_tur_sst)
        cmx = 1.0
    case(nvis_tur_hst)
        cmx = 1.5
    end select
    
    pv  => mb_pv(nb)%fld
    dpv => mb_dpv(nb)%fld
    qt  => mb_qt(nb)%fld
    dqt => mb_dqt(nb)%fld
    vsl => mb_vsl(nb)%fld
    vst => mb_vst(nb)%fld
    dst => mb_dst(nb)%fld
    qvst => mb_qvst(nb)%fld

    ni = mb_top(nb)%nijk(1)
    nj = mb_top(nb)%nijk(2)
    nk = mb_top(nb)%nijk(3)

    st(:) = mb_top(nb)%ndst(:)
    ed(:) = mb_top(nb)%nded(:)

    do k=st(3),ed(3)
    do j=st(2),ed(2)
    do i=st(1),ed(1)
        ro = pv(1)%r3d(i,j,k)
        dk = dqt(1)%r3d(i,j,k)/ro
        dw = dqt(2)%r3d(i,j,k)/ro

        tke  = qt(1)%r3d(i,j,k) + dk
        omeg = qt(2)%r3d(i,j,k) + dw

        tke  = min(max(tke,klim),kmaxtur)
        omeg = max(omeg,owlim)
        qt(1)%r3d(i,j,k) = tke 
        qt(2)%r3d(i,j,k) = omeg

        visl = vsl(1)%r3d(i,j,k)
        dist = dst(1)%r3d(i,j,k)

        dist2 = dist*dist

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

        ! Menter's limiter,|Wij|
        !!vort = sqrt(two*wij2)
        ! Hellsten's limiter,|Sij|
        vort = sqrt(two*sij2s)
         
        tmp1 = sqrt(tke)/(betas*omeg*dist)
        tmp2 = 500.D0*visl/(ro*omeg*dist2*reue)
        arg1 = max(2.D0*tmp1,tmp2)
        arg1 = min(arg1,100.D0)
        fsst = tanh(cmx*arg1*arg1)


        vist = reue*ro*tke/max(omeg,fsst*vort/a1)
          
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

end subroutine menter_update