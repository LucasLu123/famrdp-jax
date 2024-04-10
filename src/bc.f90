
subroutine boundary_conditions
    use mod_kndconsts, only : kind_int
    use mod_constants, only : bc_cut1to1
    use mod_constants, only : bc_wall,bc_symmetry
    use mod_constants, only : bc_farfield,bc_inflow
    use mod_constants, only : bc_outflow,bc_pole
    use mod_constants, only : nsw_dir_close,nvis_euler,nscmp_non
    use mod_datatypes, only : top_block_t,bc_region_t
    use mod_variables, only : nsw_kdir,nvis,nscmp
    use mod_fieldvars, only : nblkcoms,blkcoms
    implicit none
    integer(kind_int)          :: nc,nb,ns,nr,ierr
    integer(kind_int)          :: nregs,s_nd,bctype,subtype
    type(top_block_t), pointer :: top
    type(bc_region_t), pointer :: reg
    logical :: doit

    do nc=1,nblkcoms
        nb  =  blkcoms(nc)%nb
        top => blkcoms(nc)%top

        nregs = top%nregions
        do ns=1,nregs
            nr  =  top%bcnrs(ns)
            reg => top%bcs(nr)

            s_nd    = reg%s_nd
            bctype  = reg%bctype
            subtype = reg%subtype

            doit = .true.
            if (nsw_kdir == nsw_dir_close .and. s_nd == 3) then
               doit = .false.
            end if

            if (doit) then
                select case(bctype)
                case(bc_symmetry)
                    call set_bc_symmetry_plane(nb,nr)
                case(bc_farfield)
                    if (nscmp > nscmp_non) then
                        call set_bc_farfield_SCMP(nb,nr)
                    else
                        call set_bc_farfield(nb,nr)
                    end if
                case(bc_inflow)
                    call set_bc_inflow(nb,nr)
                case(bc_outflow)
                    call set_bc_outflow(nb,nr)
                case(bc_pole)
                    call set_bc_pole(nb,nr)
                case default

                end select
            end if
        end do
    end do
    
    do nc=1,nblkcoms
        nb  =  blkcoms(nc)%nb
        top => blkcoms(nc)%top

        nregs = top%nregions
        do ns=1,nregs
            nr  =  top%bcnrs(ns)
            reg => top%bcs(nr)

            s_nd    = reg%s_nd
            bctype  = reg%bctype
            subtype = reg%subtype

            doit = .true.
            if (nsw_kdir == nsw_dir_close .and. s_nd == 3) then
               doit = .false.
            end if

            if (doit) then
                select case(bctype)
                case(bc_wall)
                    if(nvis > nvis_euler) then
                        if (nscmp > nscmp_non) then
                            call set_bc_vis_wall_scmp(nb,nr)
                        else
                            call set_bc_vis_wall(nb,nr)
                        end if
                    else
                        if (nscmp > nscmp_non) then
                            !!call set_bc_inv_wall_scmp(nb,nr)
                            call set_bc_symmetry(nb,nr)
                        else
                            call set_bc_inv_wall(nb,nr)
                        end if
                    end if
                case default

                end select
            end if
        end do
    end do    

end subroutine boundary_conditions

subroutine fill_corner_mb_pv
    use mod_kndconsts, only : kind_int
    use mod_fieldvars, only : nblkcoms,blkcoms,npvs,mb_pv
    use mod_interface, only : fill_corner_var_nb
    implicit none
    integer(kind_int) :: nc,nb

    do nc=1,nblkcoms
        nb  =  blkcoms(nc)%nb

        call fill_corner_var_nb(nb,mb_pv,1,npvs)
    end do

end subroutine fill_corner_mb_pv

subroutine fill_corner_mb_pv_sp
    use mod_kndconsts, only : kind_int
    use mod_fieldvars, only : nblkcoms,blkcomssp,npvs,mb_pvfp
    use mod_interface, only : fill_corner_var_nb_sp
    implicit none
    integer(kind_int) :: nc,nb

    do nc=1,nblkcoms
        nb  =  blkcomssp(nc)%nb

        call fill_corner_var_nb_sp(nb,mb_pvfp,1,npvs)
    end do

end subroutine fill_corner_mb_pv_sp

subroutine set_bc_vis_wall(nb,nr)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,two,mide
    use mod_constants, only : id_ro,id_u,id_v,id_w,id_ps
    use mod_constants, only : nscheme_policy_edge2
    use mod_datatypes, only : bc_region_t,fld_array_t
    use mod_variables, only : tlimit,plimit
    use mod_variables, only : refbeta,twlnd,nscheme
    use mod_fieldvars, only : mb_top,mb_pv,npvs
    implicit none
    integer(kind_int), intent(in) :: nb,nr
    integer(kind_int)          :: i,j,k,m,n,ierr
    integer(kind_int)          :: s_st(3),s_ed(3)
    integer(kind_int)          :: s_nd,s_lr,s_lr3d(3)
    integer(kind_int)          :: ijkg(3),ijkb(3),ijk2(3),ijk3(3),ijk4(3)
    real(kind_real)            :: pv2(npvs),pv3(npvs),pv4(npvs)
    real(kind_real)            :: rb,ub,vb,wb,pb,tb,t2,t3,t4
    type(bc_region_t), pointer :: reg
    type(fld_array_t), pointer :: pv(:)

    reg     => mb_top(nb)%bcs(nr)
    s_st(:) = reg%s_st(:)
    s_ed(:) = reg%s_ed(:)
    s_nd    = reg%s_nd
    s_lr    = reg%s_lr

    s_lr3d(:) = s_lr*mide(:,s_nd)

    pv => mb_pv(nb)%fld

    do k=s_st(3),s_ed(3)
    do j=s_st(2),s_ed(2)
    do i=s_st(1),s_ed(1)
        ijkb(:) = (/i,j,k/)
        ijk2(:) = ijkb(:) - s_lr3d(:)
        ijk3(:) = ijk2(:) - s_lr3d(:)
        ijk4(:) = ijk3(:) - s_lr3d(:)

        do m=1,npvs
            pv2(m) = pv(m)%r3d(ijk2(1),ijk2(2),ijk2(3))
            pv3(m) = pv(m)%r3d(ijk3(1),ijk3(2),ijk3(3))
            pv4(m) = pv(m)%r3d(ijk4(1),ijk4(2),ijk4(3))
        end do

        if (nscheme == nscheme_policy_edge2) then
            pb = ( 4.0*pv2(id_ps) - pv3(id_ps) ) / 3.0
            !!pb = pv2(id_ps)
        else
            !!pb = ( 4.0*pv2(id_ps) - pv3(id_ps) ) / 3.0
            pb = ( 18.0*pv2(id_ps) - 9.0*pv3(id_ps) + 2.0*pv4(id_ps)) / 11.0
        end if

        if (pb < plimit(1)) then
            pb = pv2(id_ps)
        end if

        if (twlnd > zero) then
            tb = twlnd
        else
            t2 = pv2(id_ps)/(refbeta*pv2(id_ro))
            t3 = pv3(id_ps)/(refbeta*pv3(id_ro))
            t4 = pv4(id_ps)/(refbeta*pv4(id_ro))
            if (nscheme == nscheme_policy_edge2) then
                tb = ( 4.0*t2 - t3 ) / 3.0
                !!tb = t2
            else
                !!tb = ( 4.0*t2 - t3 ) / 3.0
                tb = ( 18.0*t2 - 9.0*t3 + 2.0*t4) / 11.0
            end if

            if (tb < tlimit(1)) then
                tb = t2
            end if
        end if

        rb = pb / (refbeta*tb)
        ub = zero
        vb = zero
        wb = zero

        pv(id_ro)%r3d(i,j,k) = rb
        pv(id_u )%r3d(i,j,k) = ub
        pv(id_v )%r3d(i,j,k) = vb
        pv(id_w )%r3d(i,j,k) = wb
        pv(id_ps)%r3d(i,j,k) = pb

        do n=1,3
            ijk2(:) = ijkb(:) - n*s_lr3d(:)
            ijkg(:) = ijkb(:) + n*s_lr3d(:)

            pb = pv(id_ps)%r3d(ijk2(1),ijk2(2),ijk2(3))
            if (twlnd > zero) then
                ijk3(:) = ijkg(:) - s_lr3d(:)
                ijk4(:) = ijk3(:) - s_lr3d(:)

                do m=id_ro,id_ps,id_ps-id_ro
                    pv3(m) = pv(m)%r3d(ijk3(1),ijk3(2),ijk3(3))
                    pv4(m) = pv(m)%r3d(ijk4(1),ijk4(2),ijk4(3))
                end do
                t3 = pv3(id_ps)/(refbeta*pv3(id_ro))
                t4 = pv4(id_ps)/(refbeta*pv4(id_ro))

                tb = 2.0*t3 - t4
                if (tb < tlimit(1)) then
                    tb = t3
                end if
                rb = pb / (refbeta*tb)
            else
                rb = pv(id_ro)%r3d(ijk2(1),ijk2(2),ijk2(3))
            end if

            pv(id_ro)%r3d(ijkg(1),ijkg(2),ijkg(3)) = rb
            pv(id_u )%r3d(ijkg(1),ijkg(2),ijkg(3)) = -pv(id_u)%r3d(ijk2(1),ijk2(2),ijk2(3))
            pv(id_v )%r3d(ijkg(1),ijkg(2),ijkg(3)) = -pv(id_v)%r3d(ijk2(1),ijk2(2),ijk2(3))
            pv(id_w )%r3d(ijkg(1),ijkg(2),ijkg(3)) = -pv(id_w)%r3d(ijk2(1),ijk2(2),ijk2(3))
            pv(id_ps)%r3d(ijkg(1),ijkg(2),ijkg(3)) = pb
        end do
    end do
    end do
    end do

end subroutine set_bc_vis_wall
    
subroutine set_bc_vis_wall_scmp(nb,nr)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,two,mide,scmp_sigma
    use mod_constants, only : id_ro,id_u,id_v,id_w,id_ps
    use mod_constants, only : nscheme_policy_edge2
    use mod_datatypes, only : bc_region_t,fld_array_t
    use mod_variables, only : tlimit,plimit,gamma,moo,poo
    use mod_variables, only : refbeta,twlnd,nscheme
    use mod_fieldvars, only : mb_top,mb_pv,npvs
    implicit none
    integer(kind_int), intent(in) :: nb,nr
    integer(kind_int)          :: i,j,k,m,n,ierr
    integer(kind_int)          :: s_st(3),s_ed(3)
    integer(kind_int)          :: s_nd,s_lr,s_lr3d(3)
    integer(kind_int)          :: ijkg(3),ijkb(3),ijk2(3),ijk3(3),ijk4(3)
    real(kind_real)            :: pv2(npvs),pv3(npvs),pv4(npvs)
    real(kind_real)            :: rb,ub,vb,wb,pb,tb,t2,t3,t4,pPrime,gama
    type(bc_region_t), pointer :: reg
    type(fld_array_t), pointer :: pv(:)

    reg     => mb_top(nb)%bcs(nr)
    s_st(:) = reg%s_st(:)
    s_ed(:) = reg%s_ed(:)
    s_nd    = reg%s_nd
    s_lr    = reg%s_lr

    s_lr3d(:) = s_lr*mide(:,s_nd)

    pv => mb_pv(nb)%fld
    gama = gamma

    do k=s_st(3),s_ed(3)
    do j=s_st(2),s_ed(2)
    do i=s_st(1),s_ed(1)
        ijkb(:) = (/i,j,k/)
        ijk2(:) = ijkb(:) - s_lr3d(:)
        ijk3(:) = ijk2(:) - s_lr3d(:)
        ijk4(:) = ijk3(:) - s_lr3d(:)

        do m=1,npvs
            pv2(m) = pv(m)%r3d(ijk2(1),ijk2(2),ijk2(3))
            pv3(m) = pv(m)%r3d(ijk3(1),ijk3(2),ijk3(3))
            pv4(m) = pv(m)%r3d(ijk4(1),ijk4(2),ijk4(3))
        end do

        if (nscheme == nscheme_policy_edge2) then
            pb = ( 4.0*pv2(id_ps) - pv3(id_ps) ) / 3.0
            !!pb = pv2(id_ps)
        else
            !!pb = ( 4.0*pv2(id_ps) - pv3(id_ps) ) / 3.0
            pb = ( 18.0*pv2(id_ps) - 9.0*pv3(id_ps) + 2.0*pv4(id_ps)) / 11.0
        end if

        if (pb < plimit(1)) then
            pb = pv2(id_ps)
        end if

        pPrime = pb - poo

        rb = (1. + gama*moo*moo*pPrime)**(1./scmp_sigma)
        ub = zero
        vb = zero
        wb = zero

        pv(id_ro)%r3d(i,j,k) = rb
        pv(id_u )%r3d(i,j,k) = ub
        pv(id_v )%r3d(i,j,k) = vb
        pv(id_w )%r3d(i,j,k) = wb
        pv(id_ps)%r3d(i,j,k) = pb

        do n=1,3
            ijk2(:) = ijkb(:) - n*s_lr3d(:)
            ijkg(:) = ijkb(:) + n*s_lr3d(:)

            pb = pv(id_ps)%r3d(ijk2(1),ijk2(2),ijk2(3))
            pPrime = pb - poo
            rb = (1. + gama*moo*moo*pPrime)**(1./scmp_sigma)

            pv(id_ro)%r3d(ijkg(1),ijkg(2),ijkg(3)) = rb
            pv(id_u )%r3d(ijkg(1),ijkg(2),ijkg(3)) = -pv(id_u)%r3d(ijk2(1),ijk2(2),ijk2(3))
            pv(id_v )%r3d(ijkg(1),ijkg(2),ijkg(3)) = -pv(id_v)%r3d(ijk2(1),ijk2(2),ijk2(3))
            pv(id_w )%r3d(ijkg(1),ijkg(2),ijkg(3)) = -pv(id_w)%r3d(ijk2(1),ijk2(2),ijk2(3))
            pv(id_ps)%r3d(ijkg(1),ijkg(2),ijkg(3)) = pb
        end do
    end do
    end do
    end do

end subroutine set_bc_vis_wall_scmp

subroutine set_bc_inv_wall(nb,nr)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : half,one,two,mide,m3x3
    use mod_constants, only : sml_ssf,nscheme_policy_edge2
    use mod_constants, only : id_ro,id_u,id_v,id_w,id_ps
    use mod_datatypes, only : bc_region_t,fld_array_t
    use mod_variables, only : nscheme,gamma,rlimit,plimit
    use mod_fieldvars, only : mb_top,npvs,mb_pv,mb_sxyz
    implicit none
    integer(kind_int), intent(in) :: nb,nr
    integer(kind_int)          :: i,j,k,m,n,ierr
    integer(kind_int)          :: s_st(3),s_ed(3),m1,m2,m3
    integer(kind_int)          :: s_nd,s_lr,s_lr3d(3)
    integer(kind_int)          :: ijkg(3),ijkb(3),ijk2(3),ijk3(3),ijk4(3)
    real(kind_real)            :: pve(npvs),pv2(npvs),pv3(npvs),pv4(npvs)
    real(kind_real)            :: nx,ny,nz,on,rb,ub,vb,wb,pb
    real(kind_real)            :: gama,ae,oae,vn,ce,se
    type(bc_region_t), pointer   :: reg
    type(fld_array_t), pointer   :: pv(:),sxyz(:)

    gama = gamma
    ae   = gama - one
    oae  = one/ae

    reg     => mb_top(nb)%bcs(nr)
    s_st(:) = reg%s_st(:)
    s_ed(:) = reg%s_ed(:)
    s_nd    = reg%s_nd
    s_lr    = reg%s_lr

    s_lr3d(:) = s_lr*mide(:,s_nd)

    m1 = m3x3(1,s_nd)
    m2 = m3x3(2,s_nd)
    m3 = m3x3(3,s_nd)

    pv   => mb_pv(nb)%fld
    sxyz => mb_sxyz(nb)%fld

    do k=s_st(3),s_ed(3)
    do j=s_st(2),s_ed(2)
    do i=s_st(1),s_ed(1)
        ijkb(:) = (/i,j,k/)
        ijk2(:) = ijkb(:) - s_lr3d(:)
        ijk3(:) = ijk2(:) - s_lr3d(:)
        ijk4(:) = ijk3(:) - s_lr3d(:)

        do m=1,npvs
            pv2(m) = pv(m)%r3d(ijk2(1),ijk2(2),ijk2(3))
            pv3(m) = pv(m)%r3d(ijk3(1),ijk3(2),ijk3(3))
            pv4(m) = pv(m)%r3d(ijk4(1),ijk4(2),ijk4(3))
        end do

        if (nscheme == nscheme_policy_edge2) then
            do m=1,npvs
                pve(m) = 2.0*pv2(m) - pv3(m)
            end do
        else
            do m=1,npvs
                pve(m) = 3.0*pv2(m) - 3.0*pv3(m) + pv4(m)
            end do
        end if

        if (pve(id_ro) < rlimit(1)) then
            pve(id_ro) = pv2(id_ro)
        end if
        if (pve(id_ps) < plimit(1)) then
            pve(id_ps) = pv2(id_ps)
        end if

        nx = sxyz(m1)%r3d(i,j,k)
        ny = sxyz(m2)%r3d(i,j,k)
        nz = sxyz(m3)%r3d(i,j,k)
        on = s_lr/max(sqrt(nx*nx + ny*ny + nz*nz),sml_ssf)
        nx = nx*on
        ny = ny*on
        nz = nz*on

        vn = nx*pve(id_u) + ny*pve(id_v) + nz*pve(id_w)
        ce = sqrt(gama*pve(id_ps)/pve(id_ro))
        se = pve(id_ps)/(pve(id_ro)**gama)

        rb = ((ce + half*ae*vn)**2/(gama*se))**oae
        ub = pve(id_u ) - nx*vn
        vb = pve(id_v ) - ny*vn
        wb = pve(id_w ) - nz*vn
        pb = se*(rb**gama)

        pv(id_ro)%r3d(i,j,k) = rb
        pv(id_u )%r3d(i,j,k) = ub
        pv(id_v )%r3d(i,j,k) = vb
        pv(id_w )%r3d(i,j,k) = wb
        pv(id_ps)%r3d(i,j,k) = pb

        do n=1,3
            ijk2(:) = ijkb(:) - n*s_lr3d(:)
            ijkg(:) = ijkb(:) + n*s_lr3d(:)

            ub = pv(id_u)%r3d(ijk2(1),ijk2(2),ijk2(3))
            vb = pv(id_v)%r3d(ijk2(1),ijk2(2),ijk2(3))
            wb = pv(id_w)%r3d(ijk2(1),ijk2(2),ijk2(3))
            vn = two*(nx*ub + ny*vb + nz*wb)

            pv(id_ro)%r3d(ijkg(1),ijkg(2),ijkg(3)) = pv(id_ro)%r3d(ijk2(1),ijk2(2),ijk2(3))
            pv(id_u )%r3d(ijkg(1),ijkg(2),ijkg(3)) = ub - nx*vn
            pv(id_v )%r3d(ijkg(1),ijkg(2),ijkg(3)) = vb - ny*vn
            pv(id_w )%r3d(ijkg(1),ijkg(2),ijkg(3)) = wb - nz*vn
            pv(id_ps)%r3d(ijkg(1),ijkg(2),ijkg(3)) = pv(id_ps)%r3d(ijk2(1),ijk2(2),ijk2(3))
        end do
    end do
    end do
    end do

end subroutine set_bc_inv_wall

subroutine set_bc_inv_wall_scmp(nb,nr)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,one,two,mide,m3x3,sml_ssf,scmp_sigma
    use mod_constants, only : id_ro,id_u,id_v,id_w,id_ps
    use mod_constants, only : nscheme_policy_edge2
    use mod_datatypes, only : bc_region_t,fld_array_t
    use mod_variables, only : nscheme,rlimit,plimit,moo,poo
    use mod_fieldvars, only : mb_top,npvs,mb_pv,mb_sxyz
    implicit none
    integer(kind_int), intent(in) :: nb,nr
    integer(kind_int)            :: i,j,k,m,n,ierr
    integer(kind_int)            :: s_st(3),s_ed(3),m1,m2,m3
    integer(kind_int)            :: s_nd,s_lr,s_lr3d(3)
    integer(kind_int)            :: ijkg(3),ijkb(3),ijk2(3)
    integer(kind_int), parameter :: nmax=3
    real(kind_real)              :: pve(npvs),pv2(npvs,nmax)
    real(kind_real)              :: vn,nx,ny,nz,on,rb,ub,vb,wb,pb,pPrime,gama
    type(bc_region_t), pointer   :: reg
    type(fld_array_t), pointer   :: pv(:),sxyz(:)

    reg     => mb_top(nb)%bcs(nr)
    s_st(:) = reg%s_st(:)
    s_ed(:) = reg%s_ed(:)
    s_nd    = reg%s_nd
    s_lr    = reg%s_lr

    s_lr3d(:) = s_lr*mide(:,s_nd)

    m1 = m3x3(1,s_nd)
    m2 = m3x3(2,s_nd)
    m3 = m3x3(3,s_nd)

    pv   => mb_pv(nb)%fld
    sxyz => mb_sxyz(nb)%fld

    nx = zero
    ny = zero
    nz = zero
    do k=s_st(3),s_ed(3)
    do j=s_st(2),s_ed(2)
    do i=s_st(1),s_ed(1)
        nx = nx + sxyz(m1)%r3d(i,j,k)
        ny = ny + sxyz(m2)%r3d(i,j,k)
        nz = nz + sxyz(m3)%r3d(i,j,k)
    end do
    end do
    end do
    on = one/max(sqrt(nx*nx + ny*ny + nz*nz),sml_ssf)
    nx = nx*on
    ny = ny*on
    nz = nz*on

    do k=s_st(3),s_ed(3)
    do j=s_st(2),s_ed(2)
    do i=s_st(1),s_ed(1)
        ijkb(:) = (/i,j,k/)

        do n=1,nmax
            ijk2(:) = ijkb(:) - n*s_lr3d(:)

            do m=1,npvs
                pv2(m,n) = pv(m)%r3d(ijk2(1),ijk2(2),ijk2(3))
            end do

            vn = nx*pv2(id_u,n) + ny*pv2(id_v,n) + nz*pv2(id_w,n)
            pv2(id_u,n) = pv2(id_u,n) - nx*vn
            pv2(id_v,n) = pv2(id_v,n) - ny*vn
            pv2(id_w,n) = pv2(id_w,n) - nz*vn
        end do

        if (nscheme == nscheme_policy_edge2) then
            do m=1,npvs
                pve(m) = (4.0*pv2(m,1) - pv2(m,2))/3.0
            end do
        else
            do m=1,npvs
                pve(m) = (4.0*pv2(m,1) - pv2(m,2))/3.0
                !!pve(m) = ( 18.0*pv2(m,1) -  9.0*pv2(m,2) +  2.0*pv2(m,3)) / 11.0
                !!pve(m) = ( 48.0*pv2(m,1) - 36.0*pv2(m,2) + 16.0*pv2(m,3) - 3.0*pv2(m,4)) / 25.0
            end do
        end if

        if (pve(id_ro) < rlimit(1)) then
            pve(id_ro) = pv2(id_ro,1)
        end if

        if (pve(id_ps) < plimit(1)) then
            pve(id_ps) = pv2(id_ps,1)
        end if

        rb = pve(id_ro)
        ub = pve(id_u )
        vb = pve(id_v )
        wb = pve(id_w )
        pb = pve(id_ps)

        pPrime = pb - poo
        rb = (1. + gama*moo*moo*pPrime)**(1./scmp_sigma)
        
        pv(id_ro)%r3d(i,j,k) = rb
        pv(id_u )%r3d(i,j,k) = ub
        pv(id_v )%r3d(i,j,k) = vb
        pv(id_w )%r3d(i,j,k) = wb
        pv(id_ps)%r3d(i,j,k) = pb

        do n=1,3
            ijk2(:) = ijkb(:) - n*s_lr3d(:)
            ijkg(:) = ijkb(:) + n*s_lr3d(:)

            ub = pv(id_u)%r3d(ijk2(1),ijk2(2),ijk2(3))
            vb = pv(id_v)%r3d(ijk2(1),ijk2(2),ijk2(3))
            wb = pv(id_w)%r3d(ijk2(1),ijk2(2),ijk2(3))
            vn = two*(nx*ub + ny*vb + nz*wb)
            
            pb = pv(id_ps)%r3d(ijk2(1),ijk2(2),ijk2(3))
            pPrime = pb - poo
            rb = (1. + gama*moo*moo*pPrime)**(1./scmp_sigma)

            pv(id_ro)%r3d(ijkg(1),ijkg(2),ijkg(3)) = rb
            pv(id_u )%r3d(ijkg(1),ijkg(2),ijkg(3)) = ub - nx*vn
            pv(id_v )%r3d(ijkg(1),ijkg(2),ijkg(3)) = vb - ny*vn
            pv(id_w )%r3d(ijkg(1),ijkg(2),ijkg(3)) = wb - nz*vn
            pv(id_ps)%r3d(ijkg(1),ijkg(2),ijkg(3)) = pb
        end do
    end do
    end do
    end do

end subroutine set_bc_inv_wall_scmp

subroutine set_bc_symmetry(nb,nr)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,one,two,mide,m3x3,sml_ssf
    use mod_constants, only : id_ro,id_u,id_v,id_w,id_ps
    use mod_constants, only : nscheme_policy_edge2
    use mod_datatypes, only : bc_region_t,fld_array_t
    use mod_variables, only : nscheme,rlimit,plimit
    use mod_fieldvars, only : mb_top,npvs,mb_pv,mb_sxyz
    implicit none
    integer(kind_int), intent(in) :: nb,nr
    integer(kind_int)            :: i,j,k,m,n,ierr
    integer(kind_int)            :: s_st(3),s_ed(3),m1,m2,m3
    integer(kind_int)            :: s_nd,s_lr,s_lr3d(3)
    integer(kind_int)            :: ijkg(3),ijkb(3),ijk2(3)
    integer(kind_int), parameter :: nmax=3
    real(kind_real)              :: pve(npvs),pv2(npvs,nmax)
    real(kind_real)              :: vn,nx,ny,nz,on,rb,ub,vb,wb,pb
    type(bc_region_t), pointer   :: reg
    type(fld_array_t), pointer   :: pv(:),sxyz(:)

    reg     => mb_top(nb)%bcs(nr)
    s_st(:) = reg%s_st(:)
    s_ed(:) = reg%s_ed(:)
    s_nd    = reg%s_nd
    s_lr    = reg%s_lr

    s_lr3d(:) = s_lr*mide(:,s_nd)

    m1 = m3x3(1,s_nd)
    m2 = m3x3(2,s_nd)
    m3 = m3x3(3,s_nd)

    pv   => mb_pv(nb)%fld
    sxyz => mb_sxyz(nb)%fld

    do k=s_st(3),s_ed(3)
    do j=s_st(2),s_ed(2)
    do i=s_st(1),s_ed(1)
        ijkb(:) = (/i,j,k/)

        nx = sxyz(m1)%r3d(i,j,k)
        ny = sxyz(m2)%r3d(i,j,k)
        nz = sxyz(m3)%r3d(i,j,k)
        on = one/max(sqrt(nx*nx + ny*ny + nz*nz),sml_ssf)
        nx = nx*on
        ny = ny*on
        nz = nz*on

        do n=1,nmax
            ijk2(:) = ijkb(:) - n*s_lr3d(:)

            do m=1,npvs
                pv2(m,n) = pv(m)%r3d(ijk2(1),ijk2(2),ijk2(3))
            end do

            vn = nx*pv2(id_u,n) + ny*pv2(id_v,n) + nz*pv2(id_w,n)
            pv2(id_u,n) = pv2(id_u,n) - nx*vn
            pv2(id_v,n) = pv2(id_v,n) - ny*vn
            pv2(id_w,n) = pv2(id_w,n) - nz*vn
        end do

        if (nscheme == nscheme_policy_edge2) then
            do m=1,npvs
                pve(m) = (4.0*pv2(m,1) - pv2(m,2))/3.0
            end do
        else
            do m=1,npvs
                !!pve(m) = (4.0*pv2(m,1) - pv2(m,2))/3.0
                pve(m) = ( 18.0*pv2(m,1) -  9.0*pv2(m,2) +  2.0*pv2(m,3)) / 11.0
                !!pve(m) = ( 48.0*pv2(m,1) - 36.0*pv2(m,2) + 16.0*pv2(m,3) - 3.0*pv2(m,4)) / 25.0
            end do
        end if

        if (pve(id_ro) < rlimit(1)) then
            pve(id_ro) = pv2(id_ro,1)
        end if

        if (pve(id_ps) < plimit(1)) then
            pve(id_ps) = pv2(id_ps,1)
        end if

        rb = pve(id_ro)
        ub = pve(id_u )
        vb = pve(id_v )
        wb = pve(id_w )
        pb = pve(id_ps)

        pv(id_ro)%r3d(i,j,k) = rb
        pv(id_u )%r3d(i,j,k) = ub
        pv(id_v )%r3d(i,j,k) = vb
        pv(id_w )%r3d(i,j,k) = wb
        pv(id_ps)%r3d(i,j,k) = pb

        do n=1,3
            ijk2(:) = ijkb(:) - n*s_lr3d(:)
            ijkg(:) = ijkb(:) + n*s_lr3d(:)

            ub = pv(id_u)%r3d(ijk2(1),ijk2(2),ijk2(3))
            vb = pv(id_v)%r3d(ijk2(1),ijk2(2),ijk2(3))
            wb = pv(id_w)%r3d(ijk2(1),ijk2(2),ijk2(3))
            vn = two*(nx*ub + ny*vb + nz*wb)

            pv(id_ro)%r3d(ijkg(1),ijkg(2),ijkg(3)) = pv(id_ro)%r3d(ijk2(1),ijk2(2),ijk2(3))
            pv(id_u )%r3d(ijkg(1),ijkg(2),ijkg(3)) = ub - nx*vn
            pv(id_v )%r3d(ijkg(1),ijkg(2),ijkg(3)) = vb - ny*vn
            pv(id_w )%r3d(ijkg(1),ijkg(2),ijkg(3)) = wb - nz*vn
            pv(id_ps)%r3d(ijkg(1),ijkg(2),ijkg(3)) = pv(id_ps)%r3d(ijk2(1),ijk2(2),ijk2(3))
        end do
    end do
    end do
    end do

end subroutine set_bc_symmetry

subroutine set_bc_symmetry_plane(nb,nr)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,one,two,mide,m3x3,sml_ssf
    use mod_constants, only : id_ro,id_u,id_v,id_w,id_ps
    use mod_constants, only : nscheme_policy_edge2
    use mod_datatypes, only : bc_region_t,fld_array_t
    use mod_variables, only : nscheme,rlimit,plimit
    use mod_fieldvars, only : mb_top,npvs,mb_pv,mb_sxyz
    implicit none
    integer(kind_int), intent(in) :: nb,nr
    integer(kind_int)            :: i,j,k,m,n,ierr
    integer(kind_int)            :: s_st(3),s_ed(3),m1,m2,m3
    integer(kind_int)            :: s_nd,s_lr,s_lr3d(3)
    integer(kind_int)            :: ijkg(3),ijkb(3),ijk2(3)
    integer(kind_int), parameter :: nmax=3
    real(kind_real)              :: pve(npvs),pv2(npvs,nmax)
    real(kind_real)              :: vn,nx,ny,nz,on,rb,ub,vb,wb,pb
    type(bc_region_t), pointer   :: reg
    type(fld_array_t), pointer   :: pv(:),sxyz(:)

    reg     => mb_top(nb)%bcs(nr)
    s_st(:) = reg%s_st(:)
    s_ed(:) = reg%s_ed(:)
    s_nd    = reg%s_nd
    s_lr    = reg%s_lr

    s_lr3d(:) = s_lr*mide(:,s_nd)

    m1 = m3x3(1,s_nd)
    m2 = m3x3(2,s_nd)
    m3 = m3x3(3,s_nd)

    pv   => mb_pv(nb)%fld
    sxyz => mb_sxyz(nb)%fld

    nx = zero
    ny = zero
    nz = zero
    do k=s_st(3),s_ed(3)
    do j=s_st(2),s_ed(2)
    do i=s_st(1),s_ed(1)
        nx = nx + sxyz(m1)%r3d(i,j,k)
        ny = ny + sxyz(m2)%r3d(i,j,k)
        nz = nz + sxyz(m3)%r3d(i,j,k)
    end do
    end do
    end do
    on = one/max(sqrt(nx*nx + ny*ny + nz*nz),sml_ssf)
    nx = nx*on
    ny = ny*on
    nz = nz*on

    do k=s_st(3),s_ed(3)
    do j=s_st(2),s_ed(2)
    do i=s_st(1),s_ed(1)
        ijkb(:) = (/i,j,k/)

        do n=1,nmax
            ijk2(:) = ijkb(:) - n*s_lr3d(:)

            do m=1,npvs
                pv2(m,n) = pv(m)%r3d(ijk2(1),ijk2(2),ijk2(3))
            end do

            vn = nx*pv2(id_u,n) + ny*pv2(id_v,n) + nz*pv2(id_w,n)
            pv2(id_u,n) = pv2(id_u,n) - nx*vn
            pv2(id_v,n) = pv2(id_v,n) - ny*vn
            pv2(id_w,n) = pv2(id_w,n) - nz*vn
        end do

        if (nscheme == nscheme_policy_edge2) then
            do m=1,npvs
                pve(m) = (4.0*pv2(m,1) - pv2(m,2))/3.0
            end do
        else
            do m=1,npvs
                pve(m) = (4.0*pv2(m,1) - pv2(m,2))/3.0
                !!pve(m) = ( 18.0*pv2(m,1) -  9.0*pv2(m,2) +  2.0*pv2(m,3)) / 11.0
                !!pve(m) = ( 48.0*pv2(m,1) - 36.0*pv2(m,2) + 16.0*pv2(m,3) - 3.0*pv2(m,4)) / 25.0
            end do
        end if

        if (pve(id_ro) < rlimit(1)) then
            pve(id_ro) = pv2(id_ro,1)
        end if

        if (pve(id_ps) < plimit(1)) then
            pve(id_ps) = pv2(id_ps,1)
        end if

        rb = pve(id_ro)
        ub = pve(id_u )
        vb = pve(id_v )
        wb = pve(id_w )
        pb = pve(id_ps)

        pv(id_ro)%r3d(i,j,k) = rb
        pv(id_u )%r3d(i,j,k) = ub
        pv(id_v )%r3d(i,j,k) = vb
        pv(id_w )%r3d(i,j,k) = wb
        pv(id_ps)%r3d(i,j,k) = pb

        do n=1,3
            ijk2(:) = ijkb(:) - n*s_lr3d(:)
            ijkg(:) = ijkb(:) + n*s_lr3d(:)

            ub = pv(id_u)%r3d(ijk2(1),ijk2(2),ijk2(3))
            vb = pv(id_v)%r3d(ijk2(1),ijk2(2),ijk2(3))
            wb = pv(id_w)%r3d(ijk2(1),ijk2(2),ijk2(3))
            vn = two*(nx*ub + ny*vb + nz*wb)

            pv(id_ro)%r3d(ijkg(1),ijkg(2),ijkg(3)) = pv(id_ro)%r3d(ijk2(1),ijk2(2),ijk2(3))
            pv(id_u )%r3d(ijkg(1),ijkg(2),ijkg(3)) = ub - nx*vn
            pv(id_v )%r3d(ijkg(1),ijkg(2),ijkg(3)) = vb - ny*vn
            pv(id_w )%r3d(ijkg(1),ijkg(2),ijkg(3)) = wb - nz*vn
            pv(id_ps)%r3d(ijkg(1),ijkg(2),ijkg(3)) = pv(id_ps)%r3d(ijk2(1),ijk2(2),ijk2(3))
        end do
    end do
    end do
    end do

end subroutine set_bc_symmetry_plane

subroutine set_bc_symmetry_plane1(nb,nr)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : mide    
    use mod_constants, only : id_ro,id_u,id_v,id_w,id_ps
    use mod_constants, only : nscheme_policy_edge2    
    use mod_datatypes, only : bc_region_t,fld_array_t
    use mod_variables, only : nscheme,rlimit,plimit    
    use mod_variables, only : rinf,vinf,pinf,tinf,nstep
    use mod_variables, only : gamma,roo,uoo,voo,woo,poo,moo
    use mod_fieldvars, only : mb_top,mb_pv,mb_qc,npvs
    implicit none
    integer(kind_int), intent(in) :: nb,nr
    integer(kind_int)          :: i,j,k,m,n,ierr
    integer(kind_int)          :: s_st(3),s_ed(3)
    integer(kind_int)          :: s_nd,s_lr,s_lr3d(3)
    integer(kind_int)          :: ijkg(3),ijkb(3),ijk2(3)
    integer(kind_int), parameter :: nmax=5
    real(kind_real)              :: pve(npvs),pv2(npvs,nmax)
    real(kind_real)            :: rb,ub,vb,wb,pb,ab,tb
    real(kind_real)            :: p0,t0,a0,Rmf,temp,cos,cp
    real(kind_real)            :: prim(1:npvs),con(1:npvs)
    type(bc_region_t), pointer :: reg
    type(fld_array_t), pointer :: pv(:),qc(:)

    reg     => mb_top(nb)%bcs(nr)
    s_st(:) = reg%s_st(:)
    s_ed(:) = reg%s_ed(:)
    s_nd    = reg%s_nd
    s_lr    = reg%s_lr 
    
    s_lr3d(:) = s_lr*mide(:,s_nd) 

    pv   => mb_pv(nb)%fld
    qc   => mb_qc(nb)%fld

    p0 = 1.3494300518134715025906735751295/gamma/moo/moo
    t0 = 845.6/tinf
    a0 = sqrt(t0/moo/moo)
    
    cp = 1.0/moo/moo/(gamma-1.0)
    
    do k=s_st(3),s_ed(3)
    do j=s_st(2),s_ed(2)
    do i=s_st(1),s_ed(1)
        
        ijkb(:) = (/i,j,k/)
        do n=1,nmax
            ijk2(:) = ijkb(:) - n*s_lr3d(:)

            do m=1,npvs
                pv2(m,n) = pv(m)%r3d(ijk2(1),ijk2(2),ijk2(3))
            end do
            
        end do

        if (nscheme == nscheme_policy_edge2) then
            do m=1,npvs
                pve(m) = pv2(m,2) !2.0*pv2(m,1) - pv2(m,2)
            end do
        else
            do m=1,npvs
                !!pve(m) = 2.0*pv2(m,1) - pv2(m,2)
                pve(m) = 3.0*pv2(m,1) -  3.0*pv2(m,2) +  pv2(m,3)
            end do
        end if
        
        !!if (pve(id_ro) < rlimit(1)) then
        !!    pve(id_ro) = pv2(id_ro,1)
        !!end if
        !!
        !!if (pve(id_ps) < plimit(1)) then
        !!    pve(id_ps) = pv2(id_ps,1)
        !!end if        
        rb   = pve(id_ro)
        ub   = pve(id_u )
        vb   = pve(id_v )
        wb   = pve(id_w )
        pb   = pve(id_ps)
        if(gamma*pb/rb>0.0) then
            ab   = sqrt(gamma*pb/rb)
            Rmf  = ub-2.0*ab/(gamma-1.0)
            cos  = -ub/(sqrt(ub*ub+vb*vb+wb*wb)+1.0e-30)
            temp = cos*cos + 2.0/(gamma-1.0)
            if(temp*a0*a0-(gamma-1.0)/2.0*Rmf*Rmf>0.0) then
                ab   = 1.0/temp*( -Rmf + cos*sqrt(temp*a0*a0-(gamma-1.0)/2.0*Rmf*Rmf) )
                tb   = t0*(ab/a0)*(ab/a0)
                if(t0-tb>0.0) then
                    pb   = p0*(tb/t0)**(gamma/(gamma-1.0))
                    rb   = pb*gamma*moo*moo/tb
                    ub   = sqrt(2.0*cp*(t0-tb))
                    vb   = 0.0
                    wb   = 0.0
                else
                    pb  = p0
                    rb  = pb*gamma*moo*moo/tb
                    ub  = 0.0
                    vb  = 0.0
                    wb  = 0.0                    
                end if
            else
                pb  = p0
                rb  = pb*gamma*moo*moo/tb
                ub  = 0.0
                vb  = 0.0
                wb  = 0.0                
            end if
        else
            pb  = p0
            rb  = pb*gamma*moo*moo/tb
            ub  = 0.0
            vb  = 0.0
            wb  = 0.0
        endif

        pv(id_ro)%r3d(i,j,k) = rb
        pv(id_u )%r3d(i,j,k) = ub
        pv(id_v )%r3d(i,j,k) = vb
        pv(id_w )%r3d(i,j,k) = wb
        pv(id_ps)%r3d(i,j,k) = pb
        do n=1,3
            ijkg(:) = ijkb(:) + n*s_lr3d(:)

            pv(id_ro)%r3d(ijkg(1),ijkg(2),ijkg(3)) = rb
            pv(id_u )%r3d(ijkg(1),ijkg(2),ijkg(3)) = ub
            pv(id_v )%r3d(ijkg(1),ijkg(2),ijkg(3)) = vb
            pv(id_w )%r3d(ijkg(1),ijkg(2),ijkg(3)) = wb
            pv(id_ps)%r3d(ijkg(1),ijkg(2),ijkg(3)) = pb
        end do
        
        do n=-1,0
            ijkg(:) = ijkb(:) + n*s_lr3d(:)

            pv(id_ro)%r3d(ijkg(1),ijkg(2),ijkg(3)) = rb
            pv(id_u )%r3d(ijkg(1),ijkg(2),ijkg(3)) = ub
            pv(id_v )%r3d(ijkg(1),ijkg(2),ijkg(3)) = vb
            pv(id_w )%r3d(ijkg(1),ijkg(2),ijkg(3)) = wb
            pv(id_ps)%r3d(ijkg(1),ijkg(2),ijkg(3)) = pb
            
            prim(id_ro) = rb
            prim(id_u ) = ub
            prim(id_v ) = vb
            prim(id_w ) = wb
            prim(id_ps) = pb
            
            call prim2con(1,npvs,prim,1,npvs,con)
            do m=1,npvs
                qc(m)%r3d(ijkg(1),ijkg(2),ijkg(3)) = con(m)
            end do
        end do
        
    end do
    end do
    end do

end subroutine set_bc_symmetry_plane1

!!subroutine set_bc_symmetry_plane(nb,nr)
!!    use mod_kndconsts, only : kind_int,kind_real
!!    use mod_constants, only : mide,zero,half,one,two,sml_ssf,m3x3
!!    use mod_constants, only : id_ro,id_u,id_v,id_w,id_ps
!!    use mod_datatypes, only : bc_region_t,fld_array_t
!!    use mod_variables, only : rinf,vinf,pinf,tinf
!!    use mod_variables, only : sigspg,alpcst,betcst
!!    use mod_variables, only : gamma,roo,uoo,voo,woo,poo,moo
!!    use mod_fieldvars, only : mb_top,mb_pv,mb_sxyz
!!    implicit none
!!    integer(kind_int), intent(in) :: nb,nr
!!    integer(kind_int)          :: i,j,k,m,n,ierr
!!    integer(kind_int)          :: s_st(3),s_ed(3)
!!    integer(kind_int)          :: s_nd,s_lr,s_lr3d(3)
!!    integer(kind_int)          :: ijkg(3),ijkb(3)
!!    integer(kind_int)          :: m1,m2,m3    
!!    real(kind_real)            :: rb,ub,vb,wb,pb,tb
!!    real(kind_real)            :: p0,t0,a0,gam1,gam2,ogam1
!!    real(kind_real)            :: v2in,vnin,c02mg,nx,ny,nz,on,vabs
!!    real(kind_real)            :: eta,eta2,rie,disc,cb,cb2,c02,vnb
!!    type(bc_region_t), pointer :: reg
!!    type(fld_array_t), pointer :: pv(:),sxyz(:)
!!
!!    reg     => mb_top(nb)%bcs(nr)
!!    s_st(:) = reg%s_st(:)
!!    s_ed(:) = reg%s_ed(:)
!!    s_nd    = reg%s_nd
!!    s_lr    = reg%s_lr 
!!    
!!    s_lr3d(:) = s_lr*mide(:,s_nd) 
!!    
!!    m1 = m3x3(1,s_nd)    
!!    m2 = m3x3(2,s_nd)    
!!    m3 = m3x3(3,s_nd)    
!!
!!    pv   => mb_pv(nb)%fld
!!    sxyz => mb_sxyz(nb)%fld
!!
!!    p0  = 1.3494300518134715025906735751295/gamma/moo/moo
!!    t0  = 845.6/tinf
!!    a0  = sqrt(t0/moo/moo)
!!    c02 = a0*a0
!!    
!!    gam1  = gamma-1.0
!!    gam2  = gamma/gam1
!!    ogam1 = one/gam1
!!    
!!    do k=s_st(3),s_ed(3)
!!    do j=s_st(2),s_ed(2)
!!    do i=s_st(1),s_ed(1)
!!        
!!        ijkb(:) = (/i,j,k/)
!!        ijkg(:) = ijkb(:) - s_lr3d(:)         
!!        rb   = pv(id_ro)%r3d(ijkg(1),ijkg(2),ijkg(3))
!!        ub   = pv(id_u )%r3d(ijkg(1),ijkg(2),ijkg(3))
!!        vb   = pv(id_v )%r3d(ijkg(1),ijkg(2),ijkg(3))
!!        wb   = pv(id_w )%r3d(ijkg(1),ijkg(2),ijkg(3))
!!        pb   = pv(id_ps)%r3d(ijkg(1),ijkg(2),ijkg(3))
!!        
!!        nx = sxyz(m1)%r3d(i,j,k)
!!        ny = sxyz(m2)%r3d(i,j,k)
!!        nz = sxyz(m3)%r3d(i,j,k)
!!        on = s_lr/max(sqrt(nx*nx + ny*ny + nz*nz),sml_ssf)
!!        nx = nx*on
!!        ny = ny*on
!!        nz = nz*on        
!!        
!!        v2in = ub*ub + vb*vb + wb*wb 
!!        vnin = nx*ub + ny*vb + nz*wb 
!!        c02mg = sqrt(c02 - half*gam1*v2in) 
!!        if ( .true. ) then 
!!            vabs = sqrt(v2in) 
!!            if ( vabs > 1.0e-6 ) then ! avoid ill-defined angle computation 
!!                eta = vnin/vabs 
!!            else 
!!                eta = -one 
!!            end if ! sl 
!!        else 
!!            eta = -one 
!!        end if 
!!        eta2 = eta*eta 
!!        
!!        rie = c02mg - half*gam1*vnin 
!!        disc = half*gam1*eta2*(c02/(rie*rie)*(one + half*gam1*eta2) - one) 
!!        
!!        if ( disc < zero ) then ! discriminant cannot be negative 
!!            tb = t0 
!!            pb = p0 
!!            vnb = zero 
!!        else 
!!            cb = rie/(one + half*gam1*eta2)*(one+sqrt(disc)) 
!!            cb2 = cb*cb 
!!            tb = moo*moo*cb2 
!!            pb = p0*(tb/t0)**gam2 
!!            if ( c02 > cb2 ) then 
!!                vnb = sqrt(two*ogam1*(c02 - cb2)) 
!!            else 
!!                vnb = zero 
!!            end if 
!!        end if 
!!        
!!        rb = gamma*moo*moo*pb/tb 
!!        ! ub = -vnb*nx 
!!        ! vb = -vnb*ny 
!!        ! wb = -vnb*nz 
!!        ub = vnb*cos(alpcst)*cos(betcst) 
!!        vb = vnb*sin(alpcst)*cos(betcst) 
!!        wb = vnb*sin(betcst)
!!
!!        pv(id_ro)%r3d(i,j,k) = rb
!!        pv(id_u )%r3d(i,j,k) = ub
!!        pv(id_v )%r3d(i,j,k) = vb
!!        pv(id_w )%r3d(i,j,k) = wb
!!        pv(id_ps)%r3d(i,j,k) = pb
!!        do n=1,3
!!            ijkg(:) = ijkb(:) + n*s_lr3d(:)
!!
!!            pv(id_ro)%r3d(ijkg(1),ijkg(2),ijkg(3)) = rb
!!            pv(id_u )%r3d(ijkg(1),ijkg(2),ijkg(3)) = ub
!!            pv(id_v )%r3d(ijkg(1),ijkg(2),ijkg(3)) = vb
!!            pv(id_w )%r3d(ijkg(1),ijkg(2),ijkg(3)) = wb
!!            pv(id_ps)%r3d(ijkg(1),ijkg(2),ijkg(3)) = pb
!!        end do
!!
!!    end do
!!    end do
!!    end do
!!
!!end subroutine set_bc_symmetry_plane

!!subroutine set_bc_symmetry_plane(nb,nr)
!!    use mod_kndconsts, only : kind_int,kind_real
!!    use mod_constants, only : mide    
!!    use mod_constants, only : id_ro,id_u,id_v,id_w,id_ps
!!    use mod_datatypes, only : bc_region_t,fld_array_t
!!    use mod_variables, only : rinf,vinf,pinf,tinf,nstep
!!    use mod_variables, only : gamma,roo,uoo,voo,woo,poo,moo
!!    use mod_fieldvars, only : mb_top,mb_pv
!!    implicit none
!!    integer(kind_int), intent(in) :: nb,nr
!!    integer(kind_int)          :: i,j,k,m,n,ierr
!!    integer(kind_int)          :: s_st(3),s_ed(3)
!!    integer(kind_int)          :: s_nd,s_lr,s_lr3d(3)
!!    integer(kind_int)          :: ijkg(3),ijkb(3)
!!    real(kind_real)            :: rb,ub,vb,wb,pb,ab,tb
!!    real(kind_real)            :: p0,t0,a0,Rmf,temp,cos,cp
!!    type(bc_region_t), pointer :: reg
!!    type(fld_array_t), pointer :: pv(:)
!!
!!    reg     => mb_top(nb)%bcs(nr)
!!    s_st(:) = reg%s_st(:)
!!    s_ed(:) = reg%s_ed(:)
!!    s_nd    = reg%s_nd
!!    s_lr    = reg%s_lr 
!!    
!!    s_lr3d(:) = s_lr*mide(:,s_nd) 
!!
!!    pv   => mb_pv(nb)%fld
!!
!!    p0 = 1.3494300518134715025906735751295/gamma/moo/moo
!!    t0 = 845.6/tinf
!!    a0 = sqrt(t0/moo/moo)
!!    
!!    cp = 1.0/moo/moo/(gamma-1.0)
!!    
!!    do k=s_st(3),s_ed(3)
!!    do j=s_st(2),s_ed(2)
!!    do i=s_st(1),s_ed(1)
!!        
!!        ijkb(:) = (/i,j,k/)
!!        ijkg(:) = ijkb(:) - s_lr3d(:)
!!        pb   = pv(id_ps)%r3d(ijkg(1),ijkg(2),ijkg(3))
!!        
!!        tb   = (pb/p0)**((gamma-1.0)/gamma)*t0
!!        rb   = pb*gamma*moo*moo/tb
!!        ub   = sqrt(2.0*cp*(t0-tb))
!!        vb   = 0.0
!!        wb   = 0.0
!!
!!        pv(id_ro)%r3d(i,j,k) = rb
!!        pv(id_u )%r3d(i,j,k) = ub
!!        pv(id_v )%r3d(i,j,k) = vb
!!        pv(id_w )%r3d(i,j,k) = wb
!!        pv(id_ps)%r3d(i,j,k) = pb
!!        do n=1,3
!!            ijkg(:) = ijkb(:) + n*s_lr3d(:)
!!
!!            pv(id_ro)%r3d(ijkg(1),ijkg(2),ijkg(3)) = rb
!!            pv(id_u )%r3d(ijkg(1),ijkg(2),ijkg(3)) = ub
!!            pv(id_v )%r3d(ijkg(1),ijkg(2),ijkg(3)) = vb
!!            pv(id_w )%r3d(ijkg(1),ijkg(2),ijkg(3)) = wb
!!            pv(id_ps)%r3d(ijkg(1),ijkg(2),ijkg(3)) = pb
!!        end do       
!!
!!    end do
!!    end do
!!    end do
!!
!!end subroutine set_bc_symmetry_plane


subroutine set_bc_farfield(nb,nr)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,fourth,half,one,two,mide,m3x3
    use mod_constants, only : id_ro,id_u,id_v,id_w,id_ps
    use mod_constants, only : nscheme_policy_edge2,sml_ssf
    use mod_datatypes, only : bc_region_t,fld_array_t
    use mod_variables, only : gamma,nscheme,rlimit,plimit
    use mod_variables, only : roo,uoo,voo,woo,poo
    use mod_fieldvars, only : mb_top,npvs,mb_pv,mb_sxyz
    implicit none
    integer(kind_int), intent(in) :: nb,nr
    integer(kind_int)          :: i,j,k,m,n,ierr
    integer(kind_int)          :: s_st(3),s_ed(3)
    integer(kind_int)          :: s_nd,s_lr,s_lr3d(3)
    integer(kind_int)          :: ijkg(3),ijkb(3),ijk2(3),ijk3(3),ijk4(3)
    integer(kind_int)          :: m1,m2,m3
    real(kind_real)            :: pve(npvs),pv2(npvs),pv3(npvs),pv4(npvs)
    real(kind_real)            :: gama,ogam,ae,oae,vne,ce,ma,vno,nx,ny,nz,on
    real(kind_real)            :: rb,ub,vb,wb,pb,vnb,ra,ua,va,wa,pa,vna,dvn
    real(kind_real)            :: r0,p0,rc0
    type(bc_region_t), pointer :: reg
    type(fld_array_t), pointer :: pv(:),sxyz(:)

    gama = gamma
    ogam = one/gama
    ae   = gamma - one
    oae  = one/ae

    reg     => mb_top(nb)%bcs(nr)
    s_st(:) = reg%s_st(:)
    s_ed(:) = reg%s_ed(:)
    s_nd    = reg%s_nd
    s_lr    = reg%s_lr

    s_lr3d(:) = s_lr*mide(:,s_nd)

    m1 = m3x3(1,s_nd)
    m2 = m3x3(2,s_nd)
    m3 = m3x3(3,s_nd)

    pv   => mb_pv(nb)%fld !< todo SCM-P ����ת��
    sxyz => mb_sxyz(nb)%fld

    do k=s_st(3),s_ed(3)
    do j=s_st(2),s_ed(2)
    do i=s_st(1),s_ed(1)
        ijkb(:) = (/i,j,k/)
        ijk2(:) = ijkb(:) - s_lr3d(:)
        ijk3(:) = ijk2(:) - s_lr3d(:)
        ijk4(:) = ijk3(:) - s_lr3d(:)

        do m=1,npvs
            pv2(m) = pv(m)%r3d(ijk2(1),ijk2(2),ijk2(3))
            pv3(m) = pv(m)%r3d(ijk3(1),ijk3(2),ijk3(3))
            pv4(m) = pv(m)%r3d(ijk4(1),ijk4(2),ijk4(3))
        end do

        if (nscheme == nscheme_policy_edge2) then
            do m=1,npvs
                pve(m) = 2.0*pv2(m) - pv3(m)
            end do
        else
            do m=1,npvs
                pve(m) = 3.0*pv2(m) - 3.0*pv3(m) + pv4(m)
            end do
        end if

        if (pve(id_ro) < rlimit(1)) then
            pve(id_ro) = pv2(id_ro)
        end if
        if (pve(id_ps) < plimit(1)) then
            pve(id_ps) = pv2(id_ps)
        end if

        nx = sxyz(m1)%r3d(i,j,k)
        ny = sxyz(m2)%r3d(i,j,k)
        nz = sxyz(m3)%r3d(i,j,k)
        on = s_lr/max(sqrt(nx*nx + ny*ny + nz*nz),sml_ssf)
        nx = nx*on
        ny = ny*on
        nz = nz*on

        vne = nx*pve(id_u) + ny*pve(id_v) + nz*pve(id_w)
        ce  = sqrt(gama*pve(id_ps)/pve(id_ro))
        ma  = abs(vne)/ce

        if (ma < one) then
            vno = nx*uoo + ny*voo + nz*woo
            if (vne > zero) then
                ra  = pve(id_ro)
                ua  = pve(id_u )
                va  = pve(id_v )
                wa  = pve(id_w )
                pa  = pve(id_ps)

                vna = vne
            else
                ra  = roo
                ua  = uoo
                va  = voo
                wa  = woo
                pa  = poo

                vna = vno
            end if
            r0 = half*(pve(id_ro) + roo)
            p0 = half*(pve(id_ps) + poo)

            rc0 = sqrt(gama*p0*r0)

            pb  = half*(  pve(id_ps) + poo      + rc0*(vne - vno) )
            vnb = half*( (pve(id_ps) - poo)/rc0 +     (vne + vno) )

            rb = ra*(pb/pa)**ogam

            dvn = (vnb - vna)
            ub = ua + nx*dvn
            vb = va + ny*dvn
            wb = wa + nz*dvn
        else
           if (vne > zero) then
              rb = pve(id_ro)
              ub = pve(id_u )
              vb = pve(id_v )
              wb = pve(id_w )
              pb = pve(id_ps)
           else
              rb = roo
              ub = uoo
              vb = voo
              wb = woo
              pb = poo
           end if
        end if

        pv(id_ro)%r3d(i,j,k) = rb
        pv(id_u )%r3d(i,j,k) = ub
        pv(id_v )%r3d(i,j,k) = vb
        pv(id_w )%r3d(i,j,k) = wb
        pv(id_ps)%r3d(i,j,k) = pb

        do n=1,3
            ijkg(:) = ijkb(:) + n*s_lr3d(:)

            pv(id_ro)%r3d(ijkg(1),ijkg(2),ijkg(3)) = rb
            pv(id_u )%r3d(ijkg(1),ijkg(2),ijkg(3)) = ub
            pv(id_v )%r3d(ijkg(1),ijkg(2),ijkg(3)) = vb
            pv(id_w )%r3d(ijkg(1),ijkg(2),ijkg(3)) = wb
            pv(id_ps)%r3d(ijkg(1),ijkg(2),ijkg(3)) = pb
        end do
    end do
    end do
    end do

end subroutine set_bc_farfield
    
    
subroutine set_bc_farfield_SCMP(nb,nr)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,fourth,half,one,two,mide,m3x3
    use mod_constants, only : id_ro,id_u,id_v,id_w,id_ps
    use mod_constants, only : nscheme_policy_edge2,sml_ssf
    use mod_constants, only : SCMP_sigma
    use mod_datatypes, only : bc_region_t,fld_array_t
    use mod_variables, only : gamma,nscheme,rlimit,plimit
    use mod_variables, only : roo,uoo,voo,woo,poo, moo
    use mod_fieldvars, only : mb_top,npvs,mb_pv,mb_sxyz
    implicit none
    integer(kind_int), intent(in) :: nb,nr
    integer(kind_int)          :: i,j,k,m,n,ierr
    integer(kind_int)          :: s_st(3),s_ed(3)
    integer(kind_int)          :: s_nd,s_lr,s_lr3d(3)
    integer(kind_int)          :: ijkg(3),ijkb(3),ijk2(3),ijk3(3),ijk4(3)
    integer(kind_int)          :: m1,m2,m3
    real(kind_real)            :: pvb(npvs),pvo(npvs),pve(npvs),pv2(npvs),pv3(npvs),pv4(npvs)
    real(kind_real)            :: cvb( 4  ),cvo( 4  ),cve( 4  ),cv2( 4  ),cv3( 4  ),cv4( 4  )
    real(kind_real)            :: gama,ogam,ae,oae,vne,ce,ma,vno,nx,ny,nz,on
    real(kind_real)            :: rb,ub,vb,wb,pb,vnb,ra,ua,va,wa,pa,vna,dvn
    real(kind_real)            :: r0,p0,rc0
    real(kind_real)            :: sn2,pPrime,oa2,Delta,ro, l1,l3,l4
    type(bc_region_t), pointer :: reg
    type(fld_array_t), pointer :: pv(:),sxyz(:)

    gama = gamma
    ogam = one/gama
    ae   = gamma - one
    oae  = one/ae

    reg     => mb_top(nb)%bcs(nr)
    s_st(:) = reg%s_st(:)
    s_ed(:) = reg%s_ed(:)
    s_nd    = reg%s_nd
    s_lr    = reg%s_lr

    s_lr3d(:) = s_lr*mide(:,s_nd)

    m1 = m3x3(1,s_nd)
    m2 = m3x3(2,s_nd)
    m3 = m3x3(3,s_nd)

    pv   => mb_pv(nb)%fld !< todo SCM-P ����ת��
    sxyz => mb_sxyz(nb)%fld

    do k=s_st(3),s_ed(3)
    do j=s_st(2),s_ed(2)
    do i=s_st(1),s_ed(1)
        ijkb(:) = (/i,j,k/)
        ijk2(:) = ijkb(:) - s_lr3d(:)
        ijk3(:) = ijk2(:) - s_lr3d(:)
        ijk4(:) = ijk3(:) - s_lr3d(:)

        do m=1,npvs
            pv2(m) = pv(m)%r3d(ijk2(1),ijk2(2),ijk2(3))
            pv3(m) = pv(m)%r3d(ijk3(1),ijk3(2),ijk3(3))
            pv4(m) = pv(m)%r3d(ijk4(1),ijk4(2),ijk4(3))
        end do
        
        
        nx = s_lr * sxyz(m1)%r3d(i,j,k)
        ny = s_lr * sxyz(m2)%r3d(i,j,k)
        nz = s_lr * sxyz(m3)%r3d(i,j,k)
        
        pvo(1) = roo
        pvo(2) = uoo
        pvo(3) = voo
        pvo(4) = woo
        pvo(5) = poo
        call pv2cv_SCMP ( npvs, nx,ny,nz, pv2, pvo, cvo )
        call pv2cv_SCMP ( npvs, nx,ny,nz, pv2, pv2, cv2 )
        call pv2cv_SCMP ( npvs, nx,ny,nz, pv2, pv3, cv3 )
        call pv2cv_SCMP ( npvs, nx,ny,nz, pv2, pv4, cv4 )
        
        if (nscheme == nscheme_policy_edge2) then
            do m=1,npvs
                pve(m) = 2.0*pv2(m) - pv3(m)
            end do
            do m=1,4
                cve(m) = 2.0*cv2(m) - cv3(m)
            end do
        else
            do m=1,npvs
                pve(m) = 3.0*pv2(m) - 3.0*pv3(m) + pv4(m)
            end do
            do m=1,4
                cve(m) = 3.0*cv2(m) - 3.0*cv3(m) + cv4(m)
            end do
        end if
        
        !call cv2pv_SCMP ( npvs, nx,ny,nz, pv2, cve, pve )
        sn2 = nx*nx + ny*ny + nz*nz
        vne = nx*pve(id_u) + ny*pve(id_v) + nz*pve(id_w)
        
        pPrime = pve(id_ps) - poo
        ro  = (1.0 + gama*moo*moo*pPrime)**(1.0/SCMP_sigma)
        oa2 = (gama*moo*moo*ro**(1-SCMP_sigma)) / SCMP_sigma
        Delta = sqrt( vne*vne*(1 - oa2) + sn2 )
        
        l1 = vne
        l3 = vne + Delta
        l4 = vne - Delta
        
        if (l1 > zero) then
            cvb(1) = cve(1)
            cvb(2) = cve(2)
        else
            cvb(1) = cvo(1)
            cvb(2) = cvo(2)
        end if
        
        if (l3 > zero) then
            cvb(3) = cve(3)
        else
            cvb(3) = cvo(3)
        end if
        
        if (l4 > zero) then
            cvb(4) = cve(4)
        else
            cvb(4) = cvo(4)
        end if
        
        call cv2pv_SCMP ( npvs, nx,ny,nz, pv2, cvb, pvb )
        pv(id_ro)%r3d(i,j,k) = pvb(1)
        pv(id_u )%r3d(i,j,k) = pvb(2)
        pv(id_v )%r3d(i,j,k) = pvb(3)
        pv(id_w )%r3d(i,j,k) = pvb(4)
        pv(id_ps)%r3d(i,j,k) = pvb(5)

        do n=1,3
            ijkg(:) = ijkb(:) + n*s_lr3d(:)

            pv(id_ro)%r3d(ijkg(1),ijkg(2),ijkg(3)) = pvb(1)
            pv(id_u )%r3d(ijkg(1),ijkg(2),ijkg(3)) = pvb(2)
            pv(id_v )%r3d(ijkg(1),ijkg(2),ijkg(3)) = pvb(3)
            pv(id_w )%r3d(ijkg(1),ijkg(2),ijkg(3)) = pvb(4)
            pv(id_ps)%r3d(ijkg(1),ijkg(2),ijkg(3)) = pvb(5)
        end do
    end do
    end do
    end do

end subroutine set_bc_farfield_SCMP

subroutine set_bc_inflow_O(nb,nr)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : mide
    use mod_constants, only : id_ro,id_u,id_v,id_w,id_ps
    use mod_datatypes, only : bc_region_t,fld_array_t
    use mod_variables, only : roo,uoo,voo,woo,poo
    use mod_fieldvars, only : mb_top,mb_pv
    implicit none
    integer(kind_int), intent(in) :: nb,nr
    integer(kind_int)          :: i,j,k,m,n,ierr
    integer(kind_int)          :: s_st(3),s_ed(3)
    integer(kind_int)          :: s_nd,s_lr,s_lr3d(3)
    integer(kind_int)          :: ijkg(3),ijkb(3)
    real(kind_real)            :: rb,ub,vb,wb,pb
    type(bc_region_t), pointer :: reg
    type(fld_array_t), pointer :: pv(:)

    reg     => mb_top(nb)%bcs(nr)
    s_st(:) = reg%s_st(:)
    s_ed(:) = reg%s_ed(:)
    s_nd    = reg%s_nd
    s_lr    = reg%s_lr 
    
    s_lr3d(:) = s_lr*mide(:,s_nd) 

    pv   => mb_pv(nb)%fld

    rb = roo
    ub = uoo
    vb = voo
    wb = woo
    pb = poo
    do k=s_st(3),s_ed(3)
    do j=s_st(2),s_ed(2)
    do i=s_st(1),s_ed(1)
        ijkb(:) = (/i,j,k/)

        pv(id_ro)%r3d(i,j,k) = rb
        pv(id_u )%r3d(i,j,k) = ub
        pv(id_v )%r3d(i,j,k) = vb
        pv(id_w )%r3d(i,j,k) = wb
        pv(id_ps)%r3d(i,j,k) = pb

        do n=1,3
            ijkg(:) = ijkb(:) + n*s_lr3d(:)

            pv(id_ro)%r3d(ijkg(1),ijkg(2),ijkg(3)) = rb
            pv(id_u )%r3d(ijkg(1),ijkg(2),ijkg(3)) = ub
            pv(id_v )%r3d(ijkg(1),ijkg(2),ijkg(3)) = vb
            pv(id_w )%r3d(ijkg(1),ijkg(2),ijkg(3)) = wb
            pv(id_ps)%r3d(ijkg(1),ijkg(2),ijkg(3)) = pb
        end do
    end do
    end do
    end do

end subroutine set_bc_inflow_O

subroutine set_bc_inflow2(nb,nr)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : mide    
    use mod_constants, only : id_ro,id_u,id_v,id_w,id_ps
    use mod_constants, only : nscheme_policy_edge2    
    use mod_datatypes, only : bc_region_t,fld_array_t
    use mod_variables, only : nscheme,rlimit,plimit
    use mod_variables, only : rinf,vinf,pinf,tinf
    use mod_variables, only : gamma,roo,uoo,voo,woo,poo,moo
    use mod_fieldvars, only : mb_top,mb_pv,mb_qc,npvs
    implicit none
    integer(kind_int), intent(in) :: nb,nr
    integer(kind_int)          :: i,j,k,m,n,ierr
    integer(kind_int)          :: s_st(3),s_ed(3)
    integer(kind_int)          :: s_nd,s_lr,s_lr3d(3)
    integer(kind_int)          :: ijkg(3),ijkb(3),ijk2(3)
    integer(kind_int), parameter :: nmax=5
    real(kind_real)              :: pve(npvs),pv2(npvs,nmax)
    real(kind_real)            :: rb,ub,vb,wb,pb,ab,tb
    real(kind_real)            :: p0,t0,a0,Rmf,temp,cos,cp
    real(kind_real)            :: prim(1:npvs),con(1:npvs)
    type(bc_region_t), pointer :: reg
    type(fld_array_t), pointer :: pv(:),qc(:)

    reg     => mb_top(nb)%bcs(nr)
    s_st(:) = reg%s_st(:)
    s_ed(:) = reg%s_ed(:)
    s_nd    = reg%s_nd
    s_lr    = reg%s_lr 
    
    s_lr3d(:) = s_lr*mide(:,s_nd) 

    pv   => mb_pv(nb)%fld
    qc   => mb_qc(nb)%fld

    p0 = 1.6160177646188008882309400444115/gamma/moo/moo
    t0 = 335.73/tinf
    a0 = sqrt(t0/moo/moo)
    
    cp = 1.0/moo/moo/(gamma-1.0)
    
    do k=s_st(3),s_ed(3)
    do j=s_st(2),s_ed(2)
    do i=s_st(1),s_ed(1)
        
        ijkb(:) = (/i,j,k/)
        do n=1,nmax
            ijk2(:) = ijkb(:) - n*s_lr3d(:)

            do m=1,npvs
                pv2(m,n) = pv(m)%r3d(ijk2(1),ijk2(2),ijk2(3))
            end do
            
        end do

        if (nscheme == nscheme_policy_edge2) then
            do m=1,npvs
                pve(m) = pv2(m,2) !2.0*pv2(m,1) - pv2(m,2)
            end do
        else
            do m=1,npvs
                !!pve(m) = 2.0*pv2(m,1) - pv2(m,2)
                pve(m) = 3.0*pv2(m,1) -  3.0*pv2(m,2) +  pv2(m,3)
            end do
        end if
        
        !!if (pve(id_ro) < rlimit(1)) then
        !!    pve(id_ro) = pv2(id_ro,1)
        !!end if
        !!if (pve(id_ps) < plimit(1)) then
        !!    pve(id_ps) = pv2(id_ps,1)
        !!end if
        
        rb   = pve(id_ro)
        ub   = pve(id_u )
        vb   = pve(id_v )
        wb   = pve(id_w )
        pb   = pve(id_ps)
        if(gamma*pb/rb>0.0) then
            ab   = sqrt(gamma*pb/rb)
            Rmf  = ub-2.0*ab/(gamma-1.0)
            cos  = -ub/(sqrt(ub*ub+vb*vb+wb*wb)+1.0e-30)
            temp = cos*cos + 2.0/(gamma-1.0)
            if(temp*a0*a0-(gamma-1.0)/2.0*Rmf*Rmf>0.0) then
                ab   = 1.0/temp*( -Rmf + cos*sqrt(temp*a0*a0-(gamma-1.0)/2.0*Rmf*Rmf) )
                tb   = t0*(ab/a0)*(ab/a0)
                if(t0-tb>0.0) then
                    pb   = p0*(tb/t0)**(gamma/(gamma-1.0))
                    rb   = pb*gamma*moo*moo/tb
                    ub   = sqrt(2.0*cp*(t0-tb))
                    vb   = 0.0
                    wb   = 0.0
                else
                    pb  = p0
                    rb  = pb*gamma*moo*moo/tb
                    ub  = 0.0
                    vb  = 0.0
                    wb  = 0.0                    
                end if
            else
                pb  = p0
                rb  = pb*gamma*moo*moo/tb
                ub  = 0.0
                vb  = 0.0
                wb  = 0.0                
            end if
        else
            pb  = p0
            rb  = pb*gamma*moo*moo/tb
            ub  = 0.0
            vb  = 0.0
            wb  = 0.0
        endif

        pv(id_ro)%r3d(i,j,k) = rb
        pv(id_u )%r3d(i,j,k) = ub
        pv(id_v )%r3d(i,j,k) = vb
        pv(id_w )%r3d(i,j,k) = wb
        pv(id_ps)%r3d(i,j,k) = pb
        do n=1,3
            ijkg(:) = ijkb(:) + n*s_lr3d(:)

            pv(id_ro)%r3d(ijkg(1),ijkg(2),ijkg(3)) = rb
            pv(id_u )%r3d(ijkg(1),ijkg(2),ijkg(3)) = ub
            pv(id_v )%r3d(ijkg(1),ijkg(2),ijkg(3)) = vb
            pv(id_w )%r3d(ijkg(1),ijkg(2),ijkg(3)) = wb
            pv(id_ps)%r3d(ijkg(1),ijkg(2),ijkg(3)) = pb
        end do
        
        do n=-1,0
            ijkg(:) = ijkb(:) + n*s_lr3d(:)

            pv(id_ro)%r3d(ijkg(1),ijkg(2),ijkg(3)) = rb
            pv(id_u )%r3d(ijkg(1),ijkg(2),ijkg(3)) = ub
            pv(id_v )%r3d(ijkg(1),ijkg(2),ijkg(3)) = vb
            pv(id_w )%r3d(ijkg(1),ijkg(2),ijkg(3)) = wb
            pv(id_ps)%r3d(ijkg(1),ijkg(2),ijkg(3)) = pb
            
            prim(id_ro) = rb
            prim(id_u ) = ub
            prim(id_v ) = vb
            prim(id_w ) = wb
            prim(id_ps) = pb
            
            call prim2con(1,npvs,prim,1,npvs,con)
            do m=1,npvs
                qc(m)%r3d(ijkg(1),ijkg(2),ijkg(3)) = con(m)
            end do            
        end do        

    end do
    end do
    end do

end subroutine set_bc_inflow2

subroutine set_bc_inflow3(nb,nr)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : mide
    use mod_constants, only : id_ro,id_u,id_v,id_w,id_ps
    use mod_datatypes, only : bc_region_t,fld_array_t
    use mod_variables, only : gamma,roo,uoo,voo,woo,poo,moo,nstep
    use mod_fieldvars, only : mb_top,mb_pv
    implicit none
    integer(kind_int), intent(in) :: nb,nr
    integer(kind_int)          :: i,j,k,m,n,ierr
    integer(kind_int)          :: s_st(3),s_ed(3)
    integer(kind_int)          :: s_nd,s_lr,s_lr3d(3)
    integer(kind_int)          :: ijkg(3),ijkb(3)
    real(kind_real)            :: rb,ub,vb,wb,pb,ab,tb,rf,qt
    real(kind_real)            :: dm,t0,a0,cp,h0,hb
    type(bc_region_t), pointer :: reg
    type(fld_array_t), pointer :: pv(:)

    reg     => mb_top(nb)%bcs(nr)
    s_st(:) = reg%s_st(:)
    s_ed(:) = reg%s_ed(:)
    s_nd    = reg%s_nd
    s_lr    = reg%s_lr 
    
    s_lr3d(:) = s_lr*mide(:,s_nd) 

    pv   => mb_pv(nb)%fld

    rb = roo
    ub = uoo
    vb = voo
    wb = woo
    pb = poo

!!    if(nstep<=5000) then
!!        dm = (1.0+0.2277818064620095076947868826041*sin(3.14159265*nstep/10000.0))*poo
!!        t0 = 1.0+1.7881733062764783367154402079106*sin(3.14159265*nstep/10000.0)
!!    else
        dm = 0.10575722760824735840144890044608 !0.12 !
        t0 = 1.1125 !1.0 !
!!    endif
    cp = 1.0/moo/moo/(gamma-1.0)
    h0 = cp*t0
    do k=s_st(3),s_ed(3)
    do j=s_st(2),s_ed(2)
    do i=s_st(1),s_ed(1)
        ijkb(:) = (/i,j,k/)
        ijkg(:) = ijkb(:) - s_lr3d(:)
        pb = pv(id_ps)%r3d(ijkg(1),ijkg(2),ijkg(3))
        tb = gamma*moo*moo*pb*(-gamma*pb+sqrt(gamma*gamma*pb*pb+2.0*(gamma-1.0)**2.0*dm*dm*h0))/dm/dm/(gamma-1.0)
        hb = cp*tb
        ub = sqrt(2.0*h0-2.0*hb)
        rb = dm/ub
        vb = 0.0
        wb = 0.0

        pv(id_ro)%r3d(i,j,k) = rb
        pv(id_u )%r3d(i,j,k) = ub
        pv(id_v )%r3d(i,j,k) = vb
        pv(id_w )%r3d(i,j,k) = wb
        pv(id_ps)%r3d(i,j,k) = pb
        do n=1,3
            ijkg(:) = ijkb(:) + n*s_lr3d(:)

            pv(id_ro)%r3d(ijkg(1),ijkg(2),ijkg(3)) = rb
            pv(id_u )%r3d(ijkg(1),ijkg(2),ijkg(3)) = ub
            pv(id_v )%r3d(ijkg(1),ijkg(2),ijkg(3)) = vb
            pv(id_w )%r3d(ijkg(1),ijkg(2),ijkg(3)) = wb
            pv(id_ps)%r3d(ijkg(1),ijkg(2),ijkg(3)) = pb
        end do

    end do
    end do
    end do

end subroutine set_bc_inflow3

subroutine set_bc_inflow(nb,nr)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : mide,zero,half,one,two,sml_ssf,m3x3
    use mod_constants, only : id_ro,id_u,id_v,id_w,id_ps
    use mod_datatypes, only : bc_region_t,fld_array_t
    use mod_variables, only : rinf,vinf,pinf,tinf
    use mod_variables, only : sigspg,alpcst,betcst
    use mod_variables, only : gamma,roo,uoo,voo,woo,poo,moo
    use mod_fieldvars, only : mb_top,mb_pv,mb_sxyz
    implicit none
    integer(kind_int), intent(in) :: nb,nr
    integer(kind_int)          :: i,j,k,m,n,ierr
    integer(kind_int)          :: s_st(3),s_ed(3)
    integer(kind_int)          :: s_nd,s_lr,s_lr3d(3)
    integer(kind_int)          :: ijkg(3),ijkb(3)
    integer(kind_int)          :: m1,m2,m3    
    real(kind_real)            :: rb,ub,vb,wb,pb,tb
    real(kind_real)            :: p0,t0,a0,gam1,gam2,ogam1
    real(kind_real)            :: v2in,vnin,c02mg,nx,ny,nz,on,vabs
    real(kind_real)            :: eta,eta2,rie,disc,cb,cb2,c02,vnb
    type(bc_region_t), pointer :: reg
    type(fld_array_t), pointer :: pv(:),sxyz(:)

    reg     => mb_top(nb)%bcs(nr)
    s_st(:) = reg%s_st(:)
    s_ed(:) = reg%s_ed(:)
    s_nd    = reg%s_nd
    s_lr    = reg%s_lr 
    
    s_lr3d(:) = s_lr*mide(:,s_nd) 
    
    m1 = m3x3(1,s_nd)    
    m2 = m3x3(2,s_nd)    
    m3 = m3x3(3,s_nd)    

    pv   => mb_pv(nb)%fld
    sxyz => mb_sxyz(nb)%fld

    p0  = 0.95/gamma/moo/moo
    t0  = 216.65/tinf
    a0  = sqrt(t0/moo/moo)
    c02 = a0*a0 
    
    gam1  = gamma-1.0
    gam2  = gamma/gam1
    ogam1 = one/gam1
    
    do k=s_st(3),s_ed(3)
    do j=s_st(2),s_ed(2)
    do i=s_st(1),s_ed(1)
        
        ijkb(:) = (/i,j,k/)
        ijkg(:) = ijkb(:) - s_lr3d(:)         
        rb   = pv(id_ro)%r3d(ijkg(1),ijkg(2),ijkg(3))
        ub   = pv(id_u )%r3d(ijkg(1),ijkg(2),ijkg(3))
        vb   = pv(id_v )%r3d(ijkg(1),ijkg(2),ijkg(3))
        wb   = pv(id_w )%r3d(ijkg(1),ijkg(2),ijkg(3))
        pb   = pv(id_ps)%r3d(ijkg(1),ijkg(2),ijkg(3))
        
        nx = sxyz(m1)%r3d(i,j,k)
        ny = sxyz(m2)%r3d(i,j,k)
        nz = sxyz(m3)%r3d(i,j,k)
        on = s_lr/max(sqrt(nx*nx + ny*ny + nz*nz),sml_ssf)
        nx = nx*on
        ny = ny*on
        nz = nz*on        
        
        v2in = ub*ub + vb*vb + wb*wb 
        vnin = nx*ub + ny*vb + nz*wb 
        c02mg = sqrt(c02 - half*gam1*v2in) 
        if ( .true. ) then 
            vabs = sqrt(v2in) 
            if ( vabs > 1.0e-6 ) then ! avoid ill-defined angle computation 
                eta = vnin/vabs 
            else 
                eta = -one 
            end if ! sl 
        else 
            eta = -one 
        end if 
        eta2 = eta*eta 
        
        rie = c02mg - half*gam1*vnin 
        disc = half*gam1*eta2*(c02/(rie*rie)*(one + half*gam1*eta2) - one) 
        
        if ( disc < zero ) then ! discriminant cannot be negative 
            tb = t0 
            pb = p0 
            vnb = zero 
        else 
            cb = rie/(one + half*gam1*eta2)*(one+sqrt(disc)) 
            cb2 = cb*cb 
            tb = moo*moo*cb2 
            pb = p0*(tb/t0)**gam2 
            if ( c02 > cb2 ) then 
                vnb = sqrt(two*ogam1*(c02 - cb2)) 
            else 
                vnb = zero 
            end if 
        end if 
        
        rb = gamma*moo*moo*pb/tb 
        ! ub = -vnb*nx 
        ! vb = -vnb*ny 
        ! wb = -vnb*nz 
        ub = vnb*cos(alpcst)*cos(betcst) 
        vb = vnb*sin(alpcst)*cos(betcst) 
        wb = vnb*sin(betcst)

        pv(id_ro)%r3d(i,j,k) = rb
        pv(id_u )%r3d(i,j,k) = ub
        pv(id_v )%r3d(i,j,k) = vb
        pv(id_w )%r3d(i,j,k) = wb
        pv(id_ps)%r3d(i,j,k) = pb
        do n=1,3
            ijkg(:) = ijkb(:) + n*s_lr3d(:)

            pv(id_ro)%r3d(ijkg(1),ijkg(2),ijkg(3)) = rb
            pv(id_u )%r3d(ijkg(1),ijkg(2),ijkg(3)) = ub
            pv(id_v )%r3d(ijkg(1),ijkg(2),ijkg(3)) = vb
            pv(id_w )%r3d(ijkg(1),ijkg(2),ijkg(3)) = wb
            pv(id_ps)%r3d(ijkg(1),ijkg(2),ijkg(3)) = pb
        end do

    end do
    end do
    end do

end subroutine set_bc_inflow

!!subroutine set_bc_inflow(nb,nr)
!!    use mod_kndconsts, only : kind_int,kind_real
!!    use mod_constants, only : mide    
!!    use mod_constants, only : id_ro,id_u,id_v,id_w,id_ps
!!    use mod_datatypes, only : bc_region_t,fld_array_t
!!    use mod_variables, only : rinf,vinf,pinf,tinf
!!    use mod_variables, only : gamma,roo,uoo,voo,woo,poo,moo
!!    use mod_fieldvars, only : mb_top,mb_pv
!!    implicit none
!!    integer(kind_int), intent(in) :: nb,nr
!!    integer(kind_int)          :: i,j,k,m,n,ierr
!!    integer(kind_int)          :: s_st(3),s_ed(3)
!!    integer(kind_int)          :: s_nd,s_lr,s_lr3d(3)
!!    integer(kind_int)          :: ijkg(3),ijkb(3)
!!    real(kind_real)            :: rb,ub,vb,wb,pb,ab,tb
!!    real(kind_real)            :: p0,t0,a0,Rmf,temp,cos,cp
!!    type(bc_region_t), pointer :: reg
!!    type(fld_array_t), pointer :: pv(:)
!!
!!    reg     => mb_top(nb)%bcs(nr)
!!    s_st(:) = reg%s_st(:)
!!    s_ed(:) = reg%s_ed(:)
!!    s_nd    = reg%s_nd
!!    s_lr    = reg%s_lr 
!!    
!!    s_lr3d(:) = s_lr*mide(:,s_nd) 
!!
!!    pv   => mb_pv(nb)%fld
!!
!!    p0 = 1.6160177646188008882309400444115/gamma/moo/moo
!!    t0 = 335.73/tinf
!!    a0 = sqrt(t0/moo/moo)
!!    
!!    cp = 1.0/moo/moo/(gamma-1.0)
!!    
!!    do k=s_st(3),s_ed(3)
!!    do j=s_st(2),s_ed(2)
!!    do i=s_st(1),s_ed(1)
!!        
!!        ijkb(:) = (/i,j,k/)
!!        ijkg(:) = ijkb(:) - s_lr3d(:)
!!        pb   = pv(id_ps)%r3d(ijkg(1),ijkg(2),ijkg(3))
!!        
!!        tb   = (pb/p0)**((gamma-1.0)/gamma)*t0
!!        rb   = pb*gamma*moo*moo/tb
!!        ub   = sqrt(2.0*cp*(t0-tb))
!!        vb   = 0.0
!!        wb   = 0.0        
!!
!!        pv(id_ro)%r3d(i,j,k) = rb
!!        pv(id_u )%r3d(i,j,k) = ub
!!        pv(id_v )%r3d(i,j,k) = vb
!!        pv(id_w )%r3d(i,j,k) = wb
!!        pv(id_ps)%r3d(i,j,k) = pb
!!        do n=1,3
!!            ijkg(:) = ijkb(:) + n*s_lr3d(:)
!!
!!            pv(id_ro)%r3d(ijkg(1),ijkg(2),ijkg(3)) = rb
!!            pv(id_u )%r3d(ijkg(1),ijkg(2),ijkg(3)) = ub
!!            pv(id_v )%r3d(ijkg(1),ijkg(2),ijkg(3)) = vb
!!            pv(id_w )%r3d(ijkg(1),ijkg(2),ijkg(3)) = wb
!!            pv(id_ps)%r3d(ijkg(1),ijkg(2),ijkg(3)) = pb
!!        end do        
!!
!!    end do
!!    end do
!!    end do
!!
!!end subroutine set_bc_inflow

subroutine set_bc_outflow_O(nb,nr)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,fourth,half,one,two,mide,m3x3
    use mod_constants, only : id_ro,id_u,id_v,id_w,id_ps
    use mod_constants, only : nscheme_policy_edge2,sml_ssf
    use mod_datatypes, only : bc_region_t,fld_array_t
    use mod_variables, only : gamma,nscheme,rlimit,plimit
    use mod_variables, only : roo,uoo,voo,woo,poo
    use mod_fieldvars, only : mb_top,npvs,mb_pv,mb_sxyz
    implicit none
    integer(kind_int), intent(in) :: nb,nr
    integer(kind_int)          :: i,j,k,m,n,ierr
    integer(kind_int)          :: s_st(3),s_ed(3)
    integer(kind_int)          :: s_nd,s_lr,s_lr3d(3)
    integer(kind_int)          :: ijkg(3),ijkb(3),ijk2(3)
    integer(kind_int)          :: m1,m2,m3
    real(kind_real)            :: pve(npvs),pv2(npvs),pv3(npvs),pv4(npvs)
    real(kind_real)            :: gama,ogam,ae,oae,vne,ce,ma,vno,nx,ny,nz,on
    real(kind_real)            :: rb,ub,vb,wb,pb,vnb,ra,ua,va,wa,pa,vna,dvn
    real(kind_real)            :: r0,p0,rc0
    type(bc_region_t), pointer :: reg
    type(fld_array_t), pointer :: pv(:),sxyz(:)

    gama = gamma
    ogam = one/gama
    ae   = gamma - one
    oae  = one/ae

    reg     => mb_top(nb)%bcs(nr)
    s_st(:) = reg%s_st(:)
    s_ed(:) = reg%s_ed(:)
    s_nd    = reg%s_nd
    s_lr    = reg%s_lr 
    
    s_lr3d(:) = s_lr*mide(:,s_nd) 
    
    m1 = m3x3(1,s_nd)    
    m2 = m3x3(2,s_nd)    
    m3 = m3x3(3,s_nd)  
      
    pv   => mb_pv(nb)%fld
    sxyz => mb_sxyz(nb)%fld
        
    do k=s_st(3),s_ed(3)
    do j=s_st(2),s_ed(2)
    do i=s_st(1),s_ed(1)
        ijkb(:) = (/i,j,k/)
        ijk2(:) = ijkb(:) - s_lr3d(:)

        do m=1,npvs
            pv2(m) = pv(m)%r3d(ijk2(1),ijk2(2),ijk2(3))
        end do

        do m=1,npvs
            pve(m) = pv2(m)
        end do

        if (pve(id_ro) < rlimit(1)) then
            pve(id_ro) = pv2(id_ro)
        end if
        if (pve(id_ps) < plimit(1)) then
            pve(id_ps) = pv2(id_ps)
        end if

        nx = sxyz(m1)%r3d(i,j,k)
        ny = sxyz(m2)%r3d(i,j,k)
        nz = sxyz(m3)%r3d(i,j,k)
        on = s_lr/max(sqrt(nx*nx + ny*ny + nz*nz),sml_ssf)
        nx = nx*on
        ny = ny*on
        nz = nz*on

        vne = nx*pve(id_u) + ny*pve(id_v) + nz*pve(id_w)
        ce  = sqrt(gama*pve(id_ps)/pve(id_ro))
        ma  = abs(vne)/ce

        if (ma < one) then
            rb = pve(id_ro)
            ub = pve(id_u )
            vb = pve(id_v )
            wb = pve(id_w )
            pb = pve(id_ps)

            pv(id_ro)%r3d(i,j,k) = rb
            pv(id_u )%r3d(i,j,k) = ub
            pv(id_v )%r3d(i,j,k) = vb
            pv(id_w )%r3d(i,j,k) = wb
            pv(id_ps)%r3d(i,j,k) = poo

            do n=1,3
                ijkg(:) = ijkb(:) + n*s_lr3d(:)

                pv(id_ro)%r3d(ijkg(1),ijkg(2),ijkg(3)) = rb
                pv(id_u )%r3d(ijkg(1),ijkg(2),ijkg(3)) = ub
                pv(id_v )%r3d(ijkg(1),ijkg(2),ijkg(3)) = vb
                pv(id_w )%r3d(ijkg(1),ijkg(2),ijkg(3)) = wb
                pv(id_ps)%r3d(ijkg(1),ijkg(2),ijkg(3)) = poo
            end do
        else
            rb = pve(id_ro)
            ub = pve(id_u )
            vb = pve(id_v )
            wb = pve(id_w )
            pb = pve(id_ps)

            pv(id_ro)%r3d(i,j,k) = rb
            pv(id_u )%r3d(i,j,k) = ub
            pv(id_v )%r3d(i,j,k) = vb
            pv(id_w )%r3d(i,j,k) = wb
            pv(id_ps)%r3d(i,j,k) = pb

            do n=1,3
                ijkg(:) = ijkb(:) + n*s_lr3d(:)

                pv(id_ro)%r3d(ijkg(1),ijkg(2),ijkg(3)) = rb
                pv(id_u )%r3d(ijkg(1),ijkg(2),ijkg(3)) = ub
                pv(id_v )%r3d(ijkg(1),ijkg(2),ijkg(3)) = vb
                pv(id_w )%r3d(ijkg(1),ijkg(2),ijkg(3)) = wb
                pv(id_ps)%r3d(ijkg(1),ijkg(2),ijkg(3)) = pb
            end do
        end if
    end do
    end do
    end do

end subroutine set_bc_outflow_O

subroutine set_bc_outflow(nb,nr)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : one,mide
    use mod_constants, only : id_ro,id_u,id_v,id_w,id_ps
    use mod_datatypes, only : bc_region_t,fld_array_t
    use mod_variables, only : rlimit,plimit,poo
    use mod_fieldvars, only : mb_top,npvs,mb_pv
    implicit none
    integer(kind_int), intent(in) :: nb,nr
    integer(kind_int)          :: i,j,k,m,n,ierr
    integer(kind_int)          :: s_st(3),s_ed(3)
    integer(kind_int)          :: s_nd,s_lr,s_lr3d(3)
    integer(kind_int)          :: ijkg(3),ijkb(3),ijk2(3),ijk3(3)
    real(kind_real)            :: pve(npvs),pv2(npvs),pv3(npvs)
    real(kind_real)            :: rb,ub,vb,wb,pb
    type(bc_region_t), pointer :: reg
    type(fld_array_t), pointer :: pv(:)

    reg     => mb_top(nb)%bcs(nr)
    s_st(:) = reg%s_st(:)
    s_ed(:) = reg%s_ed(:)
    s_nd    = reg%s_nd
    s_lr    = reg%s_lr 
    
    s_lr3d(:) = s_lr*mide(:,s_nd) 
    
    pv   => mb_pv(nb)%fld

    do k=s_st(3),s_ed(3)
    do j=s_st(2),s_ed(2)
    do i=s_st(1),s_ed(1)
        ijkb(:) = (/i,j,k/)
        ijk2(:) = ijkb(:) - s_lr3d(:)
        ijk3(:) = ijk2(:) - s_lr3d(:)

        do m=1,npvs
            pv2(m) = pv(m)%r3d(ijk2(1),ijk2(2),ijk2(3))
            pv3(m) = pv(m)%r3d(ijk3(1),ijk3(2),ijk3(3))
        end do

        do m=1,npvs
            pve(m) = pv2(m)
            !!pve(m) = 2.0*pv2(m) - pv3(m)
        end do
        if (pve(id_ro) < rlimit(1)) then
            pve(id_ro) = pv2(id_ro)
        end if
        if (pve(id_ps) < plimit(1)) then
            pve(id_ps) = pv2(id_ps)
        end if

        rb = pve(id_ro)
        ub = pve(id_u )
        vb = pve(id_v )
        wb = pve(id_w )
        pb = pve(id_ps)

        pv(id_ro)%r3d(i,j,k) = rb
        pv(id_u )%r3d(i,j,k) = ub
        pv(id_v )%r3d(i,j,k) = vb
        pv(id_w )%r3d(i,j,k) = wb
        pv(id_ps)%r3d(i,j,k) = poo

        do n=1,3
            ijkg(:) = ijkb(:) + n*s_lr3d(:)

            pv(id_ro)%r3d(ijkg(1),ijkg(2),ijkg(3)) = rb
            pv(id_u )%r3d(ijkg(1),ijkg(2),ijkg(3)) = ub
            pv(id_v )%r3d(ijkg(1),ijkg(2),ijkg(3)) = vb
            pv(id_w )%r3d(ijkg(1),ijkg(2),ijkg(3)) = wb
            pv(id_ps)%r3d(ijkg(1),ijkg(2),ijkg(3)) = poo
        end do
    end do
    end do
    end do

end subroutine set_bc_outflow

subroutine set_bc_pole(nb,nr)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,one,two,mide,m3x3,sml_ssf
    use mod_constants, only : bc_cut1to1,nscheme_policy_edge2
    use mod_constants, only : id_ro,id_u,id_v,id_w,id_ps
    use mod_datatypes, only : bc_region_t,fld_array_t
    use mod_variables, only : nscheme,rlimit,plimit
    use mod_fieldvars, only : mb_top,npvs,mb_pv,mb_sxyz
    implicit none
    integer(kind_int), intent(in) :: nb,nr
    integer(kind_int)          :: i,j,k,m,n,ierr
    integer(kind_int)          :: ns,m1,m2,m3,npts
    integer(kind_int)          :: nregs,bctype,subtype
    integer(kind_int)          :: npole,naxi,nrot
    integer(kind_int)          :: ist,ied,jst,jed,kst,ked
    integer(kind_int)          :: s_st(3),s_ed(3)
    integer(kind_int)          :: s_nd,s_lr,s_lr3d(3)
    integer(kind_int)          :: ijkg(3),ijkb(3),ijk2(3)
    integer(kind_int)          :: ijk3(3),ijk4(3)
    real(kind_real)            :: pve(npvs),pv2(npvs)
    real(kind_real)            :: pv3(npvs),pv4(npvs)
    real(kind_real)            :: nx,ny,nz,on,vn,sumpv(npvs)
    type(bc_region_t), pointer :: reg
    type(fld_array_t), pointer :: pv(:),sxyz(:)

    reg     => mb_top(nb)%bcs(nr)
    s_st(:) = reg%s_st(:)
    s_ed(:) = reg%s_ed(:)
    s_nd    = reg%s_nd
    s_lr    = reg%s_lr
    subtype = reg%subtype

    s_lr3d(:) = s_lr*mide(:,s_nd)

    pv   => mb_pv(nb)%fld
    sxyz => mb_sxyz(nb)%fld

    nrot = subtype
    naxi = 6 - nrot - s_nd

    npole = 0
    nregs = mb_top(nb)%nregions
    do ns=1,nregs
        reg => mb_top(nb)%bcs(ns)

        bctype = reg%bctype
        if (bctype == bc_cut1to1) then
            if ( reg%nbt == nb .and. reg%t_nd == nrot) then
                npole = 1
                exit
            end if
        end if
    end do

    if (npole == 0) then
        ijkb(:) = (mb_top(nb)%nijk(:)+1)/2
        ijkb(nrot) = 1
        i = ijkb(1)
        j = ijkb(2)
        k = ijkb(3)
        m1 = m3x3(1,nrot)
        m2 = m3x3(2,nrot)
        m3 = m3x3(3,nrot)
        nx = sxyz(m1)%r3d(i,j,k)
        ny = sxyz(m2)%r3d(i,j,k)
        nz = sxyz(m3)%r3d(i,j,k)
        on = one/max(sqrt(nx*nx + ny*ny + nz*nz),sml_ssf)
        nx = nx*on
        ny = ny*on
        nz = nz*on
    else
        nx = zero
        ny = zero
        nz = zero
    end if

    ist = s_st(naxi)
    ied = s_ed(naxi)
    jst = s_st(s_nd)
    jed = s_ed(s_nd)
    kst = s_st(nrot)
    ked = s_ed(nrot)

    npts = ked - kst + 1
    do i=ist,ied
    do j=jst,jed
        sumpv(:) = zero
        do k=kst,ked
            ijkb(1) = i*mide(naxi,1) + j*mide(s_nd,1) + k*mide(nrot,1)
            ijkb(2) = i*mide(naxi,2) + j*mide(s_nd,2) + k*mide(nrot,2)
            ijkb(3) = i*mide(naxi,3) + j*mide(s_nd,3) + k*mide(nrot,3)
            ijk2(:) = ijkb(:) - s_lr3d(:)
            ijk3(:) = ijk2(:) - s_lr3d(:)
            ijk4(:) = ijk3(:) - s_lr3d(:)

            do m=1,npvs
                pv2(m) = pv(m)%r3d(ijk2(1),ijk2(2),ijk2(3))
                pv3(m) = pv(m)%r3d(ijk3(1),ijk3(2),ijk3(3))
                pv4(m) = pv(m)%r3d(ijk4(1),ijk4(2),ijk4(3))
            end do

            if (nscheme == nscheme_policy_edge2) then
                do m=1,npvs
                    !!pve(m) = pv2(m)
                    pve(m) = 2.0*pv2(m) - pv3(m)
                end do
            else
                do m=1,npvs
                    pve(m) = 2.0*pv2(m) - pv3(m)
                    !!pve(m) = 3.0*pv2(m) - 3.0*pv3(m) + pv4(m)
                end do
            end if

            if (pve(id_ro) < rlimit(1)) then
                pve(id_ro) = pv2(id_ro)
            end if
            if (pve(id_ps) < plimit(1)) then
                pve(id_ps) = pv2(id_ps)
            end if

            sumpv(:) = sumpv(:) + pve(:)
        end do

        pve(:) = sumpv(:)/npts

        vn = nx*pve(id_u) + ny*pve(id_v) + nz*pve(id_w)
        pve(id_u) = pve(id_u) - nx*vn
        pve(id_v) = pve(id_v) - ny*vn
        pve(id_w) = pve(id_w) - nz*vn

       do k=kst,ked
            ijkb(1) = i*mide(naxi,1) + j*mide(s_nd,1) + k*mide(nrot,1)
            ijkb(2) = i*mide(naxi,2) + j*mide(s_nd,2) + k*mide(nrot,2)
            ijkb(3) = i*mide(naxi,3) + j*mide(s_nd,3) + k*mide(nrot,3)

            do m=1,npvs
                pv(m)%r3d(ijkb(1),ijkb(2),ijkb(3)) = pve(m)
            end do
        end do
    end do
    end do

    npts = ked + kst
    do k=s_st(3),s_ed(3)
    do j=s_st(2),s_ed(2)
    do i=s_st(1),s_ed(1)
        ijkb(:) = (/i,j,k/)

        do n=1,3
            ijkg(:) = ijkb(:) + n*s_lr3d(:)
            ijk2(:) = ijkb(:) - n*s_lr3d(:)
            if (npole == 0) then
                ijk2(nrot) = npts - ijk2(nrot)
            else
                m1 = ijk2(nrot)
                if (m1 < npts/2) then
                    ijk2(nrot) = npts/2 + m1 - 1
                else

                    ijk2(nrot) = m1 - npts/2 + 1
                end if
            end if

            do m=1,npvs
                pv2(m) = pv(m)%r3d(ijk2(1),ijk2(2),ijk2(3))
            end do

            vn = two*(nx*pv2(id_u) + ny*pv2(id_v) + nz*pv2(id_w))
            pv2(id_u) = pv2(id_u) - nx*vn
            pv2(id_v) = pv2(id_v) - ny*vn
            pv2(id_w) = pv2(id_w) - nz*vn

            do m=1,npvs
                pv(m)%r3d(ijkg(1),ijkg(2),ijkg(3)) = pv2(m)
            end do
        end do
    end do
    end do
    end do

end subroutine set_bc_pole

subroutine boundary_conditions_C
    use mod_kndconsts, only : kind_int
    use mod_constants, only : bc_cut1to1
    use mod_constants, only : bc_wall,bc_symmetry
    use mod_constants, only : bc_farfield,bc_inflow
    use mod_constants, only : bc_outflow,bc_pole
    use mod_constants, only : nsw_dir_close,nvis_euler
    use mod_datatypes, only : top_block_t,bc_region_t
    use mod_variables, only : nsw_kdir,nvis
    use mod_fieldvars, only : nblkcoms,blkcoms
    implicit none
    integer(kind_int)          :: nc,nb,ns,nr,ierr
    integer(kind_int)          :: nregs,s_nd,bctype,subtype
    type(top_block_t), pointer :: top
    type(bc_region_t), pointer :: reg
    logical :: doit

    do nc=1,nblkcoms
        nb  =  blkcoms(nc)%nb
        top => blkcoms(nc)%top

        nregs = top%nregions
        do ns=1,nregs
            nr  =  top%bcnrs(ns)
            reg => top%bcs(nr)

            s_nd    = reg%s_nd
            bctype  = reg%bctype
            subtype = reg%subtype

            doit = .true.
            if (nsw_kdir == nsw_dir_close .and. s_nd == 3) then
               doit = .false.
            end if

            if (doit) then
                select case(bctype)
                case(bc_wall)
                    if(nvis > nvis_euler) then
                        
                    else
                        call inviscid_wall_C(nb,nr)
                    end if
                case default

                end select
            end if
        end do
    end do

end subroutine boundary_conditions_C


subroutine modify_boundary
    use mod_kndconsts, only : kind_int
    use mod_constants, only : bc_cut1to1
    use mod_constants, only : bc_wall,bc_symmetry
    use mod_constants, only : bc_farfield,bc_inflow
    use mod_constants, only : bc_outflow,bc_pole
    use mod_constants, only : nsw_dir_close,nvis_euler
    use mod_datatypes, only : top_block_t,bc_region_t
    use mod_variables, only : nsw_kdir,nvis
    use mod_fieldvars, only : nblkcoms,blkcoms
    implicit none
    integer(kind_int)          :: nc,nb,ns,nr,ierr
    integer(kind_int)          :: nregs,s_nd,bctype,subtype
    type(top_block_t), pointer :: top
    type(bc_region_t), pointer :: reg
    logical :: doit 
    
    do nc=1,nblkcoms
        nb  =  blkcoms(nc)%nb
        top => blkcoms(nc)%top

        nregs = top%nregions
        do ns=1,nregs
            nr  =  top%bcnrs(ns)
            reg => top%bcs(nr)

            s_nd    = reg%s_nd
            bctype  = reg%bctype
            subtype = reg%subtype

            doit = .true.
            if (nsw_kdir == nsw_dir_close .and. s_nd == 3) then
               doit = .false.
            end if

            if (doit) then
                select case(bctype)
                case(bc_farfield)
                    call modify_farfield(nb,nr)
                case(bc_wall)
                    if(nvis == nvis_euler) then
                        call modify_inviscid_wall(nb,nr)
                    end if
                case default

                end select
            end if
        end do
    end do          

end subroutine modify_boundary

subroutine inviscid_wall_C(nb,nr)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : half,one,two,mide,m3x3
    use mod_constants, only : sml_ssf,nscheme_policy_edge2
    use mod_constants, only : id_ro,id_u,id_v,id_w,id_ps
    use mod_datatypes, only : bc_region_t,fld_array_t
    use mod_variables, only : nscheme,gamma,rlimit,plimit
    use mod_fieldvars, only : mb_top,npvs,mb_pv,mb_rhs,mb_sxyz
    implicit none
    integer(kind_int), intent(in) :: nb,nr
    integer(kind_int)          :: i,j,k,m,n,ierr
    integer(kind_int)          :: s_st(3),s_ed(3),m1,m2,m3
    integer(kind_int)          :: s_nd,s_lr,s_lr3d(3)
    integer(kind_int)          :: ijkb(3),ijk2(3)
    integer(kind_int), parameter :: nmax=4
    real(kind_real)            :: nx,ny,nz,on
    real(kind_real)            :: l1,l2,l3,gama,vn
    real(kind_real)            :: vector_e(1:5),vector_f(1:5),vector_g(1:5)
    real(kind_real)            :: rm,um,vm,wm,pm,cm,em,gamam1
	real(kind_real)            :: matrix_p(1:5,1:5),matrix_p1(1:5,1:5),vv,temp
    real(kind_real)            :: temp1(1:5),temp2(1:5),v2,aa,bb,nn(3)
    real(kind_real)            :: pve(npvs),pv2(npvs,nmax)
    type(bc_region_t), pointer   :: reg
    type(fld_array_t), pointer   :: pv(:),sxyz(:),rhs(:)

    gama = gamma

    reg     => mb_top(nb)%bcs(nr)
    s_st(:) = reg%s_st(:)
    s_ed(:) = reg%s_ed(:)
    s_nd    = reg%s_nd
    s_lr    = reg%s_lr

    s_lr3d(:) = s_lr*mide(:,s_nd)

    m1 = m3x3(1,s_nd)
    m2 = m3x3(2,s_nd)
    m3 = m3x3(3,s_nd)

    pv   => mb_pv(nb)%fld
    rhs  => mb_rhs(nb)%fld
    sxyz => mb_sxyz(nb)%fld

    do k=s_st(3),s_ed(3)
    do j=s_st(2),s_ed(2)
    do i=s_st(1),s_ed(1)
        nx = sxyz(m1)%r3d(i,j,k)
        ny = sxyz(m2)%r3d(i,j,k)
        nz = sxyz(m3)%r3d(i,j,k)
        on = one/max(sqrt(nx*nx + ny*ny + nz*nz),sml_ssf)
        nx = nx*on
        ny = ny*on
        nz = nz*on

        ijkb(:) = (/i,j,k/)

        do n=1,nmax
            ijk2(:) = ijkb(:) - n*s_lr3d(:)
            do m=1,npvs
                pv2(m,n) = pv(m)%r3d(ijk2(1),ijk2(2),ijk2(3))
            end do
        end do

        do m=1,npvs
            pve(m) = pv2(m,1) !( 48.0*pv2(m,1) - 36.0*pv2(m,2) + 16.0*pv2(m,3) - 3.0*pv2(m,4)) / 25.0
        end do

        rm=pv(id_ro)%r3d(i,j,k)!pve(1) !
	    um=pv(id_u)%r3d(i,j,k)!pve(2) !
	    vm=pv(id_v)%r3d(i,j,k)!pve(3) !
	    wm=pv(id_w)%r3d(i,j,k)!pve(4) !
	    pm=pv(id_ps)%r3d(i,j,k)!pve(5) !
		cm=sqrt(gama*pm/rm)

	    gamam1 = gama - 1.0
		aa=rm/2.0/cm
		bb=1.0/rm/cm
		v2=0.5*gamam1*(um*um+vm*vm+wm*wm)
	    em=gama*pm/gamam1+0.5*rm*(um*um+vm*vm+wm*wm)


	    vector_e(1)=rm*um
	    vector_e(2)=vector_e(1)*um+pm 
	    vector_e(3)=vector_e(1)*vm    
	    vector_e(4)=vector_e(1)*wm    
	    vector_e(5)=em*um             

		vector_f(1)=rm*vm
		vector_f(2)=vector_e(3)
		vector_f(3)=vector_f(1)*vm+pm 
		vector_f(4)=vector_f(1)*wm    
		vector_f(5)=em*vm             


		vector_g(1)=rm*wm
		vector_g(2)=vector_e(4)
		vector_g(3)=vector_f(4)
		vector_g(4)=vector_g(1)*wm+pm 
		vector_g(5)=em*wm

		l1=nx 
		l2=ny 
		l3=nz 

		vv=l1*um+l2*vm+l3*wm

        call matrix(rm,um,vm,wm,pm,em,cm,l1,l2,l3,vv,v2,aa,bb,gamam1,matrix_p,matrix_p1)
                      
        if(s_lr==-1) then		 
		    do m=1,npvs
                temp1(m)=0.0
		        do n=1,npvs
			        temp1(m)=temp1(m)+matrix_p(m,4)*(matrix_p1(5,n)-matrix_p1(4,n))*rhs(n)%r3d(i,j,k)
			    enddo
		    enddo

		    do m=1,npvs
		        rhs(m)%r3d(i,j,k)=rhs(m)%r3d(i,j,k)+temp1(m)
		    enddo
		elseif(s_lr==1) then
		 
		    do m=1,npvs
                temp1(m)=0.0
		        do n=1,npvs
			        temp1(m)=temp1(m)-matrix_p(m,5)*(matrix_p1(5,n)-matrix_p1(4,n))*rhs(n)%r3d(i,j,k)
			    enddo
		    enddo

		    do m=1,npvs
		        rhs(m)%r3d(i,j,k)=rhs(m)%r3d(i,j,k)+temp1(m)
		    enddo
		end if
    end do
    end do
    end do

end subroutine inviscid_wall_C

subroutine interface_C
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : one,m3x3,sml_ssf
    use mod_constants, only : id_ro,id_u,id_v,id_w,id_ps
    use mod_constants, only : nbc_inter_buf_dqc
    use mod_datatypes, only : fld_array_t
    use mod_variables, only : nscheme,gamma,rlimit,plimit
    use mod_fieldvars, only : inters,ninterf
    use mod_fieldvars, only : mb_top,mb_vol,npvs
    use mod_fieldvars, only : mb_pv,mb_rhs,mb_sxyz
    use mod_parallels    
    implicit none
    integer(kind_int)          :: nint
    integer(kind_int)          :: nbs,nrs,s_st(3),s_ed(3),s_nd,s_lr
    integer(kind_int)          :: nbt,nrt,t_st(3),t_ed(3),t_nd,t_lr
    integer(kind_int)          :: i,j,k,m,n,l,is,js,ks,it,jt,kt,ierr
    integer(kind_int)          :: ijks(3),ijkt(3),id_src,id_des
    integer(kind_int)          :: m1,m2,m3
    real(kind_real)            :: nx,ny,nz,on,gama,vn,eps
    real(kind_real)            :: rm,um,vm,wm,pm,cm,em,gamam1
	real(kind_real)            :: matrix_p(1:5,1:5),matrix_p1(1:5,1:5)
    real(kind_real)            :: diag_p(1:5),diag_n(1:5)
    real(kind_real)            :: a_p(1:5,1:5),a_n(1:5,1:5)
    real(kind_real)            :: temp1(1:5),temp2(1:5),temp,vv,v2,aa,bb
    real(kind_real)            :: chart(1:5),rhss(1:5),vol
    type(fld_array_t), pointer   :: pv(:),sxyz(:),rhs(:)

    gama = gamma
    eps  = 1.0e-4

    do nint=1,ninterf
        nbs     = inters(nint)%bc%nbs
        nrs     = inters(nint)%bc%nrs
        s_st(:) = inters(nint)%bc%s_st(:)
        s_ed(:) = inters(nint)%bc%s_ed(:)
        s_nd    = inters(nint)%bc%s_nd
        s_lr    = inters(nint)%bc%s_lr

        nbt     = inters(nint)%bc%nbt
        nrt     = inters(nint)%bc%nrt
        t_st(:) = inters(nint)%bc%t_st(:)
        t_ed(:) = inters(nint)%bc%t_ed(:)
        t_nd    = inters(nint)%bc%t_nd
        t_lr    = inters(nint)%bc%t_lr
        
#ifdef PARALLEL
        id_src = mb_top(nbs)%pid
        id_des = mb_top(nbt)%pid

        if (myid == id_des) then
#endif        

        m1 = m3x3(1,t_nd)
        m2 = m3x3(2,t_nd)
        m3 = m3x3(3,t_nd)

        pv   => mb_pv(nbt)%fld
        rhs  => mb_rhs(nbt)%fld
        sxyz => mb_sxyz(nbt)%fld
        
        inters(nint)%dat => inters(nint)%buf(nbc_inter_buf_dqc)%dat

        do k=t_st(3),t_ed(3)
        do j=t_st(2),t_ed(2)
        do i=t_st(1),t_ed(1)
            
            ijkt(:) = (/i,j,k/)
            ijks(:) = mb_top(nbt)%bcs(nrt)%mapijk(i,j,k,:)
            
            is = ijks(1)
            js = ijks(2)
            ks = ijks(3)
            it = ijkt(1)
            jt = ijkt(2)
            kt = ijkt(3)            
            
            nx = sxyz(m1)%r3d(it,jt,kt)
            ny = sxyz(m2)%r3d(it,jt,kt)
            nz = sxyz(m3)%r3d(it,jt,kt)
            on = one/max(sqrt(nx*nx + ny*ny + nz*nz),sml_ssf)
            nx = nx*on
            ny = ny*on
            nz = nz*on

            rm=pv(id_ro)%r3d(it,jt,kt)
	        um=pv(id_u)%r3d(it,jt,kt)
	        vm=pv(id_v)%r3d(it,jt,kt)
	        wm=pv(id_w)%r3d(it,jt,kt)
	        pm=pv(id_ps)%r3d(it,jt,kt)
	    	cm=sqrt(gama*pm/rm)

	        gamam1 = gama - 1.0
	    	aa=rm/2.0/cm
	    	bb=1.0/rm/cm
	    	v2=0.5*gamam1*(um*um+vm*vm+wm*wm)
	        em=gama*pm/gamam1+0.5*rm*(um*um+vm*vm+wm*wm)

	    	vv=nx*um+ny*vm+nz*wm
            
            chart(1)=vv
            chart(2)=chart(1)
            chart(3)=chart(1)
            chart(4)=chart(1)+cm
            chart(5)=chart(1)-cm
            
            do m=1,5
                diag_p(m)=(abs(chart(m))+chart(m)+eps)/2.0/(abs(chart(m))+eps)
                diag_n(m)=(abs(chart(m))-chart(m)+eps)/2.0/(abs(chart(m))+eps)
            enddo            
            
            call matrix(rm,um,vm,wm,pm,em,cm,nx,ny,nz,vv,v2,aa,bb,gamam1,matrix_p,matrix_p1)
            
            do m=1,5				     
                do n=1,5
                    a_p(m,n)=0.0
                    a_n(m,n)=0.0
                    do l=1,5		     
                        a_p(m,n)=a_p(m,n)+matrix_p(m,l)*diag_p(l)*matrix_p1(l,n)
                        a_n(m,n)=a_n(m,n)+matrix_p(m,l)*diag_n(l)*matrix_p1(l,n)
                    enddo
                enddo
            enddo
            
            do m=1,5
                rhss(m)=inters(nint)%dat(is,js,ks,m)
            enddo
            
            vol = mb_vol(nbt)%fld(1)%r3d(it,jt,kt)
            
            do m=1,5
                temp=0.0
                do n=1,5
                    if(t_lr==1) then
                        temp=temp+a_p(m,n)*rhs(n)%r3d(it,jt,kt)
                        temp=temp+a_n(m,n)*rhss(n)*vol
                    endif
                    if(t_lr==-1) then
                        temp=temp+a_p(m,n)*rhss(n)*vol
                        temp=temp+a_n(m,n)*rhs(n)%r3d(it,jt,kt)
                    endif
                enddo
                rhs(m)%r3d(it,jt,kt)=temp
            enddo
            
        end do
        end do
        end do

        nullify(inters(nint)%dat)
#ifdef PARALLEL
        end if
#endif        
    end do

end subroutine interface_C

subroutine matrix(rm,um,vm,wm,pm,em,cm,l1,l2,l3,vv,v2,aa,bb,gamam1,matrix_p,matrix_p1)
    use mod_kndconsts, only : kind_real
	implicit none
    real(kind_real) :: matrix_p(5,5),matrix_p1(5,5),l1,l2,l3,vv,v2,aa,bb,gamam1
	real(kind_real) :: rm,um,vm,wm,pm,em,cm

    matrix_p(1,1)=l1
	matrix_p(1,2)=l2
	matrix_p(1,3)=l3
	matrix_p(1,4)=aa
	matrix_p(1,5)=aa
	matrix_p(2,1)=um*l1
	matrix_p(2,2)=um*l2-rm*l3
	matrix_p(2,3)=um*l3+rm*l2
	matrix_p(2,4)=aa*(um+l1*cm)
	matrix_p(2,5)=aa*(um-l1*cm)
	matrix_p(3,1)=vm*l1+rm*l3
	matrix_p(3,2)=vm*l2
	matrix_p(3,3)=vm*l3-rm*l1
	matrix_p(3,4)=aa*(vm+l2*cm)
	matrix_p(3,5)=aa*(vm-l2*cm)
	matrix_p(4,1)=wm*l1-rm*l2
	matrix_p(4,2)=wm*l2+rm*l1
	matrix_p(4,3)=wm*l3
	matrix_p(4,4)=aa*(wm+l3*cm)
	matrix_p(4,5)=aa*(wm-l3*cm)
	matrix_p(5,1)=l1*v2/gamam1+rm*(vm*l3-wm*l2)
	matrix_p(5,2)=l2*v2/gamam1+rm*(wm*l1-um*l3)
	matrix_p(5,3)=l3*v2/gamam1+rm*(um*l2-vm*l1)				 
	matrix_p(5,4)=aa*((v2+cm*cm)/gamam1+cm*vv)
	matrix_p(5,5)=aa*((v2+cm*cm)/gamam1-cm*vv)

    matrix_p1(1,1)=l1*(1.0-v2/cm/cm)-(l3*vm-l2*wm)/rm
    matrix_p1(1,2)=l1*gamam1*um/cm/cm
    matrix_p1(1,3)=l3/rm+l1*gamam1*vm/cm/cm
	matrix_p1(1,4)=-l2/rm+l1*gamam1*wm/cm/cm
	matrix_p1(1,5)=-l1*gamam1/cm/cm
	matrix_p1(2,1)=l2*(1.0-v2/cm/cm)-(l1*wm-l3*um)/rm
	matrix_p1(2,2)=-l3/rm+l2*gamam1*um/cm/cm
	matrix_p1(2,3)=l2*gamam1*vm/cm/cm
	matrix_p1(2,4)=l1/rm+l2*gamam1*wm/cm/cm
	matrix_p1(2,5)=-l2*gamam1/cm/cm
	matrix_p1(3,1)=l3*(1.0-v2/cm/cm)-(l2*um-l1*vm)/rm
	matrix_p1(3,2)=l2/rm+l3*gamam1*um/cm/cm
	matrix_p1(3,3)=-l1/rm+l3*gamam1*vm/cm/cm
	matrix_p1(3,4)=l3*gamam1*wm/cm/cm
	matrix_p1(3,5)=-l3*gamam1/cm/cm
	matrix_p1(4,1)=bb*(v2-cm*vv)
	matrix_p1(4,2)=bb*(l1*cm-gamam1*um)
	matrix_p1(4,3)=bb*(l2*cm-gamam1*vm)
	matrix_p1(4,4)=bb*(l3*cm-gamam1*wm)
	matrix_p1(4,5)=bb*gamam1
	matrix_p1(5,1)=bb*(v2+cm*vv)
	matrix_p1(5,2)=-bb*(l1*cm+gamam1*um)
	matrix_p1(5,3)=-bb*(l2*cm+gamam1*vm)
	matrix_p1(5,4)=-bb*(l3*cm+gamam1*wm)
	matrix_p1(5,5)=bb*gamam1
    return
end subroutine matrix

subroutine modify_inviscid_wall(nb,nr)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : half,one,two,mide,m3x3
    use mod_constants, only : sml_ssf,nscheme_policy_edge2
    use mod_constants, only : id_ro,id_u,id_v,id_w,id_ps
    use mod_variables, only : gamma,poo,roo
    use mod_datatypes, only : bc_region_t,fld_array_t
    use mod_fieldvars, only : mb_top,npvs,mb_pv,mb_sxyz,mb_qc,mb_dq
    implicit none
    integer(kind_int), intent(in) :: nb,nr
    integer(kind_int)          :: i,j,k,m,n,ierr
    integer(kind_int)          :: s_st(3),s_ed(3),m1,m2,m3
    integer(kind_int)          :: s_nd,s_lr,s_lr3d(3)
    real(kind_real)            :: nx,ny,nz,on,vvv(2),vvm1,vvm2
    real(kind_real)            :: gama,vn,con(1:npvs)
    real(kind_real)            :: rm,um,vm,wm,pm,em
    type(bc_region_t), pointer   :: reg
    type(fld_array_t), pointer   :: pv(:),qc(:),dq(:),sxyz(:),rhs(:)

    gama = gamma

    reg     => mb_top(nb)%bcs(nr)
    s_st(:) = reg%s_st(:)
    s_ed(:) = reg%s_ed(:)
    s_nd    = reg%s_nd
    s_lr    = reg%s_lr

    s_lr3d(:) = s_lr*mide(:,s_nd)

    m1 = m3x3(1,s_nd)
    m2 = m3x3(2,s_nd)
    m3 = m3x3(3,s_nd)

    pv   => mb_pv(nb)%fld
    qc => mb_qc(nb)%fld
    dq => mb_dq(nb)%fld
    sxyz => mb_sxyz(nb)%fld

    do k=s_st(3),s_ed(3)
    do j=s_st(2),s_ed(2)
    do i=s_st(1),s_ed(1)
        nx = sxyz(m1)%r3d(i,j,k)
        ny = sxyz(m2)%r3d(i,j,k)
        nz = sxyz(m3)%r3d(i,j,k)
        on = one/max(sqrt(nx*nx + ny*ny + nz*nz),sml_ssf)
        nx = nx*on
        ny = ny*on
        nz = nz*on

        rm=pv(id_ro)%r3d(i,j,k)
        um=pv(id_u)%r3d(i,j,k)
        vm=pv(id_v)%r3d(i,j,k)
        wm=pv(id_w)%r3d(i,j,k)
        pm=pv(id_ps)%r3d(i,j,k)
        vvm1=sqrt(um*um+vm*vm+wm*wm)
        vvv(1)=0.0 !!disv(1,1)*nx+disv(1,2)*ny+disv(1,3)*nz
        vvv(2)=um*nx+vm*ny+wm*nz
!!        vvv(1)=disv(1,1)*nx+disv(1,2)*ny+disv(1,3)*nz
!!        vvv(2)=-nx*vm+ny*um
!!        vvv(3)=nx*wm-ny*um
!!        um=l1*vvv(1)+l2*vvv(2)-l3*vvv(3)
!!        vm=-(-l1*l2*vvv(1)+(l1*l1+l3*l3)*vvv(2)+l2*l3*vvv(3))/l1
!!        wm=(l1*l3*vvv(1)+l2*l3*vvv(2)+(l1*l1+l2*l2)*vvv(3))/l1
        pv(id_ro)%r3d(i,j,k)=rm        
        pv(id_u)%r3d(i,j,k)=um +nx*( vvv(1) -vvv(2) )
        pv(id_v)%r3d(i,j,k)=vm +ny*( vvv(1) -vvv(2) )
        pv(id_w)%r3d(i,j,k)=wm +nz*( vvv(1) -vvv(2) )
        pv(id_ps)%r3d(i,j,k)=pm !poo*rm**1.4 !
        vvm2=sqrt(pv(id_u)%r3d(i,j,k)*pv(id_u)%r3d(i,j,k)+pv(id_v)%r3d(i,j,k)*pv(id_v)%r3d(i,j,k)+ &
                  pv(id_w)%r3d(i,j,k)*pv(id_w)%r3d(i,j,k))
        pv(id_u)%r3d(i,j,k)=pv(id_u)%r3d(i,j,k) !/vvm2*vvm1
        pv(id_v)%r3d(i,j,k)=pv(id_v)%r3d(i,j,k) !/vvm2*vvm1
        pv(id_w)%r3d(i,j,k)=pv(id_w)%r3d(i,j,k) !/vvm2*vvm1
        rm=pv(id_ro)%r3d(i,j,k)
        um=pv(id_u)%r3d(i,j,k)
        vm=pv(id_v)%r3d(i,j,k)
        wm=pv(id_w)%r3d(i,j,k)
        pm=pv(id_ps)%r3d(i,j,k)
        em = pm/(gama-1.0) + 0.5*rm*( um*um + vm*vm + wm*wm )
        con(id_ro) = rm
        con(id_u) = rm * um
        con(id_v) = rm * vm
        con(id_w) = rm * wm
        con(id_ps) = em

        dq(id_ro)%r3d(i,j,k) = con(id_ro) - qc(id_ro)%r3d(i,j,k)
        dq(id_u)%r3d(i,j,k) = con(id_u) - qc(id_u)%r3d(i,j,k)
        dq(id_v)%r3d(i,j,k) = con(id_v) - qc(id_v)%r3d(i,j,k)
        dq(id_w)%r3d(i,j,k) = con(id_w) - qc(id_w)%r3d(i,j,k)
        dq(id_ps)%r3d(i,j,k) = con(id_ps) - qc(id_ps)%r3d(i,j,k)

        qc(id_ro)%r3d(i,j,k) = con(id_ro)
        qc(id_u)%r3d(i,j,k) = con(id_u)
        qc(id_v)%r3d(i,j,k) = con(id_v)
        qc(id_w)%r3d(i,j,k) = con(id_w)
        qc(id_ps)%r3d(i,j,k) = con(id_ps)

    enddo
    enddo
    enddo

end subroutine modify_inviscid_wall

subroutine modify_farfield(nb,nr)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : half,one,two,mide,m3x3
    use mod_constants, only : sml_ssf,nscheme_policy_edge2
    use mod_constants, only : id_ro,id_u,id_v,id_w,id_ps
    use mod_variables, only : gamma,roo,uoo,voo,woo,poo
    use mod_datatypes, only : bc_region_t,fld_array_t
    use mod_fieldvars, only : mb_top,npvs,mb_pv,mb_sxyz,mb_qc,mb_dq
    implicit none
    integer(kind_int), intent(in) :: nb,nr
    integer(kind_int)          :: i,j,k,m,n,ierr
    integer(kind_int)          :: s_st(3),s_ed(3),m1,m2,m3
    integer(kind_int)          :: s_nd,s_lr,s_lr3d(3)
    real(kind_real)            :: nx,ny,nz,on,vvv(2),matrix_p(1:5,1:5),matrix_p1(1:5,1:5)
    real(kind_real)            :: gama,vn,con(1:npvs),gamam1,w(2,5),matrix1_p1(5,5)
    real(kind_real)            :: rm,um,vm,wm,pm,cm,em,aa,bb,v2,vv,chart(5),rm1
    type(bc_region_t), pointer   :: reg
    type(fld_array_t), pointer   :: pv(:),qc(:),dq(:),sxyz(:),rhs(:)

    gama = gamma

    reg     => mb_top(nb)%bcs(nr)
    s_st(:) = reg%s_st(:)
    s_ed(:) = reg%s_ed(:)
    s_nd    = reg%s_nd
    s_lr    = reg%s_lr

    s_lr3d(:) = s_lr*mide(:,s_nd)

    m1 = m3x3(1,s_nd)
    m2 = m3x3(2,s_nd)
    m3 = m3x3(3,s_nd)

    pv   => mb_pv(nb)%fld
    qc => mb_qc(nb)%fld
    dq => mb_dq(nb)%fld
    sxyz => mb_sxyz(nb)%fld

    do k=s_st(3),s_ed(3)
    do j=s_st(2),s_ed(2)
    do i=s_st(1),s_ed(1)
        nx = sxyz(m1)%r3d(i,j,k)
        ny = sxyz(m2)%r3d(i,j,k)
        nz = sxyz(m3)%r3d(i,j,k)
        on = one/max(sqrt(nx*nx + ny*ny + nz*nz),sml_ssf)
        nx = nx*on
        ny = ny*on
        nz = nz*on

        rm=pv(id_ro)%r3d(i,j,k)
        um=pv(id_u)%r3d(i,j,k)
        vm=pv(id_v)%r3d(i,j,k)
        wm=pv(id_w)%r3d(i,j,k)
        pm=pv(id_ps)%r3d(i,j,k)

        rm=roo
		um=uoo
	    vm=voo
	    wm=woo
		pm=poo
		cm=sqrt(gama*poo/roo)
		gamam1 = gama - 1.0
		aa=rm/2.0/cm
		bb=1.0/rm/cm
		v2=0.5*gamam1*(um*um+vm*vm+wm*wm)
		vv=nx*um+ny*vm+nz*wm
	    em=pm/(rm*gamam1)+0.5*(um*um+vm*vm+wm*wm)

		con(id_ro)=rm
		con(id_u)=rm*um
		con(id_v)=rm*vm
		con(id_w)=rm*wm
		con(id_ps)=rm*em

		call matrix(rm,um,vm,wm,pm,em,cm,nx,ny,nz,vv,v2,aa,bb,gamam1,matrix_p,matrix_p1)			     

		do m=1,5
		    w(2,m)=0.0
		    do n=1,5
                w(2,m)=w(2,m)+matrix_p1(m,n)*con(n)
			enddo
		enddo

        rm=pv(id_ro)%r3d(i,j,k)
        um=pv(id_u)%r3d(i,j,k)
        vm=pv(id_v)%r3d(i,j,k)
        wm=pv(id_w)%r3d(i,j,k)
        pm=pv(id_ps)%r3d(i,j,k)
		cm=sqrt(gama*pm/rm)
		aa=rm/2.0/cm
		bb=1.0/rm/cm
		v2=0.5*gamam1*(um*um+vm*vm+wm*wm)
		vv=nx*um+ny*vm+nz*wm
	    em=pm/(rm*gamam1)+0.5*(um*um+vm*vm+wm*wm)
		con(id_ro)=rm
		con(id_u)=rm*um
		con(id_v)=rm*vm
		con(id_w)=rm*wm
		con(id_ps)=rm*em

		chart(1)=sxyz(m1)%r3d(i,j,k)*um+sxyz(m2)%r3d(i,j,k)*vm+sxyz(m3)%r3d(i,j,k)*wm+10e-30
		chart(2)=chart(1)
		chart(3)=chart(1)
		chart(4)=chart(1)+cm*sqrt(sxyz(m1)%r3d(i,j,k)*sxyz(m1)%r3d(i,j,k)+sxyz(m2)%r3d(i,j,k)*sxyz(m2)%r3d(i,j,k)+ &
		         sxyz(m3)%r3d(i,j,k)*sxyz(m3)%r3d(i,j,k))
		chart(5)=chart(1)-cm*sqrt(sxyz(m1)%r3d(i,j,k)*sxyz(m1)%r3d(i,j,k)+sxyz(m2)%r3d(i,j,k)*sxyz(m2)%r3d(i,j,k)+ &
		         sxyz(m3)%r3d(i,j,k)*sxyz(m3)%r3d(i,j,k))				 

		
        call matrix(rm,um,vm,wm,pm,em,cm,nx,ny,nz,vv,v2,aa,bb,gamam1,matrix_p,matrix1_p1)


		do m=1,5
		    w(1,m)=0.0
		    do n=1,5
                w(1,m)=w(1,m)+matrix1_p1(m,n)*con(n)
			 enddo
		enddo

		
		
        do m=1,5
!            if( (1.0_prec*s_lr)*chart(m)/abs(chart(m))<0.0_prec) then	
             if( (1.0*s_lr)*sign(one,chart(m))<0.0) then					     
			     w(1,m)=w(2,m)
			     do n=1,5
			         matrix1_p1(m,n)=matrix_p1(m,n)
			     enddo
			 endif
		enddo
        
        
        call INV(matrix1_p1,matrix_p,npvs)

		do m=1,5
		    con(m)=0.0
			 do n=1,5
			     con(m)=con(m)+matrix_p(m,n)*w(1,n)
			 enddo
		enddo

        dq(id_ro)%r3d(i,j,k) = con(id_ro) - qc(id_ro)%r3d(i,j,k)
        dq(id_u)%r3d(i,j,k) = con(id_u) - qc(id_u)%r3d(i,j,k)
        dq(id_v)%r3d(i,j,k) = con(id_v) - qc(id_v)%r3d(i,j,k)
        dq(id_w)%r3d(i,j,k) = con(id_w) - qc(id_w)%r3d(i,j,k)
        dq(id_ps)%r3d(i,j,k) = con(id_ps) - qc(id_ps)%r3d(i,j,k)

        qc(id_ro)%r3d(i,j,k) = con(id_ro)
        qc(id_u)%r3d(i,j,k) = con(id_u)
        qc(id_v)%r3d(i,j,k) = con(id_v)
        qc(id_w)%r3d(i,j,k) = con(id_w)
        qc(id_ps)%r3d(i,j,k) = con(id_ps)

		rm = con(id_ro)
        rm1 = 1.0/rm
        um = con(id_u) * rm1
        vm = con(id_v) * rm1
        wm = con(id_w) * rm1
        em = con(id_ps)
        pm = gamam1 * ( em - 0.5*rm*( um*um + vm*vm + wm*wm ) )

        pv(id_ro)%r3d(i,j,k)=rm
        pv(id_u)%r3d(i,j,k) =um
        pv(id_v)%r3d(i,j,k)=vm
        pv(id_w)%r3d(i,j,k)=wm
        pv(id_ps)%r3d(i,j,k)=pm

    enddo
    enddo
    enddo
end subroutine modify_farfield

!*********************************************************************************************************************************************************
subroutine INV(A,A1,ND)
!*********************************************************************************************************************************************************
!�������档
!     Input:
!     a                inputing matrix
!     nd               dimension of matrix
!     Output:
!     a1               inverted matrix
!*********************************************************************************************************************************************************
use mod_kndconsts, only : kind_int,kind_real
implicit none
integer(kind_int):: ND,I,J,K,N,IP
real(kind_real):: P,W,W1,W2
real(kind_real):: A(ND,ND),A1(ND,ND)

DO J=1,ND
	DO K=1,ND
		A1(K,J)=0.0
        A1(J,J)=1.0
	ENDDO
ENDDO

DO 10 N=1,ND
	P=0.0
	DO 2 I=N,ND
		IF(P-ABS(A(I,N))) 1,2,2
1			P=ABS(A(I,N))
            IP=I
2   CONTINUE
	W1=A(IP,N)
	if(n.ne.ip) then
		do 21 j=n,ND
			w=a(n,j)
			a(n,j)=a(ip,j)
21          a(ip,j)=w
		do 22 j=1,ND
			w=a1(n,j)
			a1(n,j)=a1(ip,j)
22          a1(ip,j)=w
	endif
	DO 3 J=N,ND
		A(N,J)=A(n,J)/W1
3   continue
	DO 4 J=1,ND
		A1(N,J)=A1(n,J)/W1
4   continue
	DO 8 I=1,ND
		W2=A(I,N)
		IF(I-N) 5,8,5
5   DO 6 J=N,ND
6       A(I,J)=A(I,J)-W2*A(N,J)
    DO 7 J=1,ND
7       A1(I,J)=A1(I,J)-W2*A1(N,J)
8   CONTINUE
10 CONTINUE

RETURN
END SUBROUTINE

subroutine set_boundary_vars(nst,ned,nb,nr,bctype,iindex,jindex,kindex,nx,ny,nz,varl,varr)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : bc_cut1to1
    use mod_constants, only : bc_wall,bc_symmetry
    use mod_constants, only : bc_farfield,bc_inflow
    use mod_constants, only : bc_outflow,bc_pole    
    use mod_constants, only : sml_ssf,m3x3,scmp_sigma
    use mod_constants, only : nvis_euler,nscmp_non
    use mod_constants, only : id_ro,id_u,id_v,id_w,id_ps
    use mod_datatypes, only : bc_region_t,fld_array_t
    use mod_variables, only : nscmp,gamma,moo,nvis,refbeta,twlnd
    use mod_variables, only : roo,uoo,voo,woo,poo
    use mod_fieldvars, only : mb_topc,mb_pvfp,mb_pv
    implicit none    
    integer(kind_int), intent(in)  :: nst,ned
    integer(kind_int), intent(in)  :: nb,nr,bctype,iindex,jindex,kindex
    real(kind_real), intent(in)  :: nx,ny,nz
    real(kind_real), intent(out) :: varl(nst:ned),varr(nst:ned)
    integer(kind_int)          :: s_lr,s_nd,m
    real(kind_real)            :: nxx,nyy,nzz,on,vn
    type(bc_region_t), pointer   :: reg
    type(fld_array_t), pointer   :: pvfp(:),pv(:)

    pvfp => mb_pvfp(nb)%fld
    pv => mb_pv(nb)%fld
    
    select case(bctype) 
    case(bc_wall)
        if(nvis > nvis_euler) then
            reg     => mb_topc(nb)%bcs(nr)
            s_lr    = reg%s_lr
            s_nd    = reg%s_nd
            
            if(s_lr<0) then
                do m=1,5
                    varr(m) = pvfp(m)%r3d(iindex,jindex,kindex)
                end do
                varl(id_ro) = varr(id_ro)
                varl(id_ps) = varr(id_ps)
                varl(id_u) = -varr(id_u)
                varl(id_v) = -varr(id_v)
                varl(id_w) = -varr(id_w)                
                if(twlnd>0.0) then
                    varl(id_ps) = varl(id_ro)*refbeta*twlnd
                end if
            else
                do m=1,5
                    varl(m) = pvfp(m)%r3d(iindex,jindex,kindex)
                end do
                varr(id_ro) = varl(id_ro)
                varr(id_ps) = varl(id_ps)                
                varr(id_u) = -varl(id_u)
                varr(id_v) = -varl(id_v)
                varr(id_w) = -varl(id_w)              
                if(twlnd>0.0) then
                    varr(id_ps) = varr(id_ro)*refbeta*twlnd
                end if                
            end if
        else
            reg     => mb_topc(nb)%bcs(nr)
            s_lr    = reg%s_lr
            on = 1.0/max(sqrt(nx*nx + ny*ny + nz*nz),sml_ssf)
            nxx = nx*on
            nyy = ny*on
            nzz = nz*on
            
            if(s_lr<0) then
                do m=1,5
                    varr(m) = pvfp(m)%r3d(iindex,jindex,kindex)
                end do
                vn = nxx*varr(id_u) + nyy*varr(id_v) + nzz*varr(id_w)
                varl(id_ro) = varr(id_ro)
                varl(id_ps) = varr(id_ps)
                !if(nscmp > nscmp_non) then
                !    varl(id_ps) = (varr(id_ro)**scmp_sigma -1.0)/gamma/moo/moo
                !end if
                varl(id_u) = varr(id_u) - 2.0*vn*nxx
                varl(id_v) = varr(id_v) - 2.0*vn*nyy
                varl(id_w) = varr(id_w) - 2.0*vn*nzz
            else
                do m=1,5
                    varl(m) = pvfp(m)%r3d(iindex,jindex,kindex)
                end do
                vn = nxx*varl(id_u) + nyy*varl(id_v) + nzz*varl(id_w)
                varr(id_ro) = varl(id_ro)
                varr(id_ps) = varl(id_ps)
                !if(nscmp > nscmp_non) then
                !    varr(id_ps) = (varl(id_ro)**scmp_sigma -1.0)/gamma/moo/moo
                !end if                 
                varr(id_u) = varl(id_u) - 2.0*vn*nxx
                varr(id_v) = varl(id_v) - 2.0*vn*nyy
                varr(id_w) = varl(id_w) - 2.0*vn*nzz
            end if
        end if        
    case(bc_symmetry)
        reg     => mb_topc(nb)%bcs(nr)
        s_lr    = reg%s_lr
        on = 1.0/max(sqrt(nx*nx + ny*ny + nz*nz),sml_ssf)
        nxx = nx*on
        nyy = ny*on
        nzz = nz*on
        
        if(s_lr<0) then
            do m=1,5
                varr(m) = pvfp(m)%r3d(iindex,jindex,kindex)
            end do
            vn = nxx*varr(id_u) + nyy*varr(id_v) + nzz*varr(id_w)
            if (nscmp > nscmp_non) then
                varl(id_ro) = varr(id_ro)
                varl(id_ps) = varr(id_ps)
            else
                varl(id_ro) = varr(id_ro)
                varl(id_ps) = varr(id_ps)
            end if
            varl(id_u) = varr(id_u) - 2.0*vn*nxx
            varl(id_v) = varr(id_v) - 2.0*vn*nyy
            varl(id_w) = varr(id_w) - 2.0*vn*nzz
        else
            do m=1,5
                varl(m) = pvfp(m)%r3d(iindex,jindex,kindex)
            end do
            vn = nxx*varl(id_u) + nyy*varl(id_v) + nzz*varl(id_w)
            if (nscmp > nscmp_non) then
                varr(id_ro) = varl(id_ro)
                varr(id_ps) = varl(id_ps)
            else
                varr(id_ro) = varl(id_ro)
                varr(id_ps) = varl(id_ps)
            end if            
            varr(id_u) = varl(id_u) - 2.0*vn*nxx
            varr(id_v) = varl(id_v) - 2.0*vn*nyy
            varr(id_w) = varl(id_w) - 2.0*vn*nzz
        end if
    case(bc_farfield)
        reg     => mb_topc(nb)%bcs(nr)
        s_lr    = reg%s_lr
        on = 1.0/max(sqrt(nx*nx + ny*ny + nz*nz),sml_ssf)
        nxx = nx*on
        nyy = ny*on
        nzz = nz*on
        
        if(s_lr<0) then
            do m=1,5
                varr(m) = pvfp(m)%r3d(iindex,jindex,kindex)
            end do
            varl(id_ro) = roo
            if (nscmp > nscmp_non) then
                varl(id_ps) = 0.0
            else
                varl(id_ps) = poo
            end if
            varl(id_u) = uoo
            varl(id_v) = voo
            varl(id_w) = woo
        else
            do m=1,5
                varl(m) = pvfp(m)%r3d(iindex,jindex,kindex)
            end do
            varr(id_ro) = roo
            if (nscmp > nscmp_non) then
                varr(id_ps) = 0.0
            else
                varr(id_ps) = poo
            end if
            varr(id_u) = uoo
            varr(id_v) = voo
            varr(id_w) = woo
        end if
    case(bc_inflow)
        !todocall set_bcvar_inflow(nb,nr)
    case(bc_outflow)
        !todocall set_bcvar_outflow(nb,nr)
    case(bc_pole)
        !todocall set_bcvar_pole(nb,nr)
    case default
    
    end select
    
end subroutine set_boundary_vars

subroutine set_boundary_gradvars(nst,ned,nb,nr,bctype,iindex,jindex,kindex,nx,ny,nz,var)
    use mod_kndconsts, only : kind_int,kind_real  
    use mod_constants, only : sml_ssf,m3x3
    use mod_datatypes, only : bc_region_t,fld_array_t
    use mod_variables, only : twlnd
    use mod_fieldvars, only : mb_topc,mb_dpvfp
    implicit none    
    integer(kind_int), intent(in)  :: nst,ned
    integer(kind_int), intent(in)  :: nb,nr,bctype,iindex,jindex,kindex
    real(kind_real), intent(in)  :: nx,ny,nz
    real(kind_real), intent(out) :: var(nst:ned)
    integer(kind_int)          :: s_lr,s_nd,m
    real(kind_real)            :: nxx,nyy,nzz,on,vn
    type(bc_region_t), pointer   :: reg
    type(fld_array_t), pointer   :: dpvfps(:),dpvfpe(:)
    
    if(twlnd<0.0) then
        reg     => mb_topc(nb)%bcs(nr)
        s_lr    = reg%s_lr
        s_nd    = reg%s_nd
        if(s_nd ==1) then
            dpvfps => mb_dpvfp(nb,1)%fld
            dpvfpe => mb_dpvfp(nb,2)%fld 
        else if(s_nd ==2) then
            dpvfps => mb_dpvfp(nb,3)%fld
            dpvfpe => mb_dpvfp(nb,4)%fld        
        else
            dpvfps => mb_dpvfp(nb,5)%fld
            dpvfpe => mb_dpvfp(nb,6)%fld         
        end if
        
        on = 1.0/max(sqrt(nx*nx + ny*ny + nz*nz),sml_ssf)
        nxx = nx*on
        nyy = ny*on
        nzz = nz*on
        
        if(s_lr<0) then
            do m=1,5
                var(m) = dpvfps(m)%r3d(iindex,jindex,kindex)
            end do
            vn = nxx*var(10) + nyy*var(11) + nzz*var(12)
            var(10) = var(10) - vn*nxx
            var(11) = var(11) - vn*nyy
            var(12) = var(12) - vn*nzz
        else
            do m=1,5
                var(m) = dpvfpe(m)%r3d(iindex,jindex,kindex)
            end do
            vn = nxx*var(10) + nyy*var(11) + nzz*var(12)
            var(10) = var(10) - vn*nxx
            var(11) = var(11) - vn*nyy
            var(12) = var(12) - vn*nzz
        end if
    end if
    
end subroutine set_boundary_gradvars

