
subroutine tur_boundary_conditions
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
                    call tur_set_bc_wall(nb,nr)
                case(bc_symmetry)
                    call tur_set_bc_symmetry(nb,nr)
                case(bc_farfield)
                    call tur_set_bc_farfield(nb,nr)
                case(bc_inflow)
                    call tur_set_bc_inflow(nb,nr)
                case(bc_outflow)
                    call tur_set_bc_outflow(nb,nr)
                case(bc_pole)
                    call tur_set_bc_pole(nb,nr)
                case default

                end select
            end if
        end do
    end do

end subroutine tur_boundary_conditions

subroutine fill_corner_mb_qvst
    use mod_kndconsts, only : kind_int
    use mod_fieldvars, only : nblkcoms,blkcoms
    use mod_fieldvars, only : nqvst,mb_qvst
    use mod_interface, only : fill_corner_var_nb
    implicit none
    integer(kind_int) :: nc,nb

    do nc=1,nblkcoms
        nb  =  blkcoms(nc)%nb

        call fill_corner_var_nb(nb,mb_qvst,1,nqvst)
    end do

end subroutine fill_corner_mb_qvst

subroutine tur_set_bc_wall(nb,nr)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : fourth,half,two3rd,one,two
    use mod_constants, only : zero,small,sml_ssf
    use mod_constants, only : mide,m3x3,id_u,id_v,id_w
    use mod_constants, only : nvis_tur_sa
    use mod_constants, only : nvis_tur_sst,nvis_tur_hst
    use mod_datatypes, only : bc_region_t,fld_array_t
    use mod_variables, only : ntursch,nvis,reue
    use mod_fieldvars, only : mb_top,mb_xyz,mb_sxyz,mb_vsl
    use mod_fieldvars, only : npvt,mb_pvt,nqvst,mb_qvst
    implicit none
    integer(kind_int), intent(in) :: nb,nr
    integer(kind_int)          :: i,j,k,m,n,ierr
    integer(kind_int)          :: s_st(3),s_ed(3)
    integer(kind_int)          :: s_nd,s_lr,s_lr3d(3)
    integer(kind_int)          :: m1,m2,m3,ijkg(3),ijkb(3)
    integer(kind_int)          :: ijk2(3),ijk3(3)
    real(kind_real)            :: nx,ny,nz,on,vx,vy,vz,vn,vt
    real(kind_real)            :: ro,visl,dist,mul,muoy2
    real(kind_real)            :: dx1,dy1,dz1,tauw,utau,ypl1,kpl1
    real(kind_real)            :: omegw,omegg,omeg2,omeg3
    type(bc_region_t), pointer :: reg
    type(fld_array_t), pointer :: xyz(:),sxyz(:),vsl(:)
    type(fld_array_t), pointer :: pvt(:),qvst(:)

    reg     => mb_top(nb)%bcs(nr)
    s_st(:) = reg%s_st(:)
    s_ed(:) = reg%s_ed(:)
    s_nd    = reg%s_nd
    s_lr    = reg%s_lr

    s_lr3d(:) = s_lr*mide(:,s_nd)

    m1 = m3x3(1,s_nd)
    m2 = m3x3(2,s_nd)
    m3 = m3x3(3,s_nd)

    xyz  => mb_xyz(nb)%fld
    sxyz => mb_sxyz(nb)%fld
    vsl  => mb_vsl(nb)%fld
    pvt  => mb_pvt(nb)%fld
    qvst => mb_qvst(nb)%fld

    do k=s_st(3),s_ed(3)
    do j=s_st(2),s_ed(2)
    do i=s_st(1),s_ed(1)
        ijkb(:) = (/i,j,k/)

        do m=1,nqvst
            qvst(m)%r3d(i,j,k) = small
        end do

        do n=1,3
            ijkg(:) = ijkb(:) + n*s_lr3d(:)
            ijk2(:) = ijkb(:) - n*s_lr3d(:)

            do m=1,nqvst
                qvst(m)%r3d(ijkg(1),ijkg(2),ijkg(3)) = -qvst(m)%r3d(ijk2(1),ijk2(2),ijk2(3))
            end do
        end do

        select case(nvis)
        case(nvis_tur_sst,nvis_tur_hst)
            ijk2(:) = ijkb(:) - s_lr3d(:)

            ro   = pvt(1)%r3d(i,j,k)
            visl = vsl(1)%r3d(i,j,k)
            dx1  = xyz(1)%r3d(ijk2(1),ijk2(2),ijk2(3)) - xyz(1)%r3d(i,j,k)
            dy1  = xyz(2)%r3d(ijk2(1),ijk2(2),ijk2(3)) - xyz(2)%r3d(i,j,k)
            dz1  = xyz(3)%r3d(ijk2(1),ijk2(2),ijk2(3)) - xyz(3)%r3d(i,j,k)

            mul  = visl/ro
            dist = sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1)
            muoy2 = mul/(dist*dist*reue + small)

            ! Menter's wall approach
            omegw = 3200.0*muoy2

            ! Hellsten's rough wall approach
            nx = sxyz(m1)%r3d(i,j,k)
            ny = sxyz(m2)%r3d(i,j,k)
            nz = sxyz(m3)%r3d(i,j,k)
            on = one/max(sqrt(nx*nx + ny*ny + nz*nz),sml_ssf)
            nx = nx*on
            ny = ny*on
            nz = nz*on

            vx = pvt(id_u)%r3d(ijk2(1),ijk2(2),ijk2(3))
            vy = pvt(id_v)%r3d(ijk2(1),ijk2(2),ijk2(3))
            vz = pvt(id_w)%r3d(ijk2(1),ijk2(2),ijk2(3))

            vn = nx*vx + ny*vy + nz*vz
            vx = vx - nx*vn
            vy = vy - ny*vn
            vz = vz - nz*vn

            vt = sqrt(vx*vx + vy*vy + vz*vz)

            utau = sqrt(mul*vt/(dist*reue))
            ypl1 = reue*dist*utau/mul

            !!kpl1 = 0.88*ypl1  !!correspond to the Menter's wall approach.
            !!kpl1 = 2.5*ypl1
            kpl1 = min(2.4*(ypl1**0.85),8.0)
            omegw = 2500.0*muoy2*(ypl1/kpl1)**2

            qvst(2)%r3d(i,j,k) = omegw

            do n=1,3
                ijkg(:) = ijkb(:) + n*s_lr3d(:)
                ijk2(:) = ijkg(:) - s_lr3d(:)
                ijk3(:) = ijk2(:) - s_lr3d(:)

                omeg2 = qvst(2)%r3d(ijk2(1),ijk2(2),ijk2(3))
                omeg3 = qvst(2)%r3d(ijk3(1),ijk3(2),ijk3(3))

                omegg = 2.0*omeg2 - omeg3
                if (omegg < small) then
                    omegg = omeg2
                end if

                qvst(2)%r3d(ijkg(1),ijkg(2),ijkg(3)) = omegg
            end do
        end select
    end do
    end do
    end do

end subroutine tur_set_bc_wall

subroutine tur_set_bc_symmetry(nb,nr)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,small,sml_ssf
    use mod_constants, only : one,two,mide
    use mod_constants, only : nvis_tur_sa
    use mod_constants, only : nvis_tur_sst,nvis_tur_hst
    use mod_constants, only : nscheme_policy_edge2
    use mod_datatypes, only : bc_region_t,fld_array_t
    use mod_variables, only : nvis,ntursch,turconsts
    use mod_fieldvars, only : mb_top,nqvst,mb_qvst
    implicit none
    integer(kind_int), intent(in) :: nb,nr
    integer(kind_int)            :: i,j,k,m,n,ierr
    integer(kind_int)            :: s_st(3),s_ed(3),m1,m2,m3
    integer(kind_int)            :: s_nd,s_lr,s_lr3d(3)
    integer(kind_int)            :: ijkg(3),ijkb(3)
    integer(kind_int)            :: ijk2(3),ijk3(3),ijk4(3)
    real(kind_real)              :: pv2(nqvst),pv3(nqvst),pv4(nqvst)
    real(kind_real)              :: pve(nqvst),turlim(nqvst)
    type(bc_region_t), pointer   :: reg
    type(fld_array_t), pointer   :: qvst(:)

    select case(nvis)
    case(nvis_tur_sa)
        turlim(1) = turconsts(12)
        turlim(2) = small
    case(nvis_tur_sst,nvis_tur_hst)
        turlim(1) = turconsts(13)
        turlim(2) = turconsts(14)
        turlim(3) = small
    end select

    reg     => mb_top(nb)%bcs(nr)
    s_st(:) = reg%s_st(:)
    s_ed(:) = reg%s_ed(:)
    s_nd    = reg%s_nd
    s_lr    = reg%s_lr

    s_lr3d(:) = s_lr*mide(:,s_nd)

    qvst => mb_qvst(nb)%fld

    do k=s_st(3),s_ed(3)
    do j=s_st(2),s_ed(2)
    do i=s_st(1),s_ed(1)
        ijkb(:) = (/i,j,k/)
        ijk2(:) = ijkb(:) - s_lr3d(:)
        ijk3(:) = ijk2(:) - s_lr3d(:)
        ijk4(:) = ijk3(:) - s_lr3d(:)

        do m=1,nqvst
            pv2(m) = qvst(m)%r3d(ijk2(1),ijk2(2),ijk2(3))
            pv3(m) = qvst(m)%r3d(ijk3(1),ijk3(2),ijk3(3))
            pv4(m) = qvst(m)%r3d(ijk4(1),ijk4(2),ijk4(3))
        end do

        if (ntursch == nscheme_policy_edge2) then
            do m=1,nqvst
                pve(m) = (4.0*pv2(m) - pv3(m))/3.0
            end do
        else
            do m=1,nqvst
                pve(m) = (18.0*pv2(m) -  9.0*pv3(m) +  2.0*pv4(m)) / 11.0
            end do
        end if

        do m=1,nqvst
            if ( pve(m) < turlim(m) ) then
                pve(m) = pv2(m)
            end if

            qvst(m)%r3d(i,j,k) = pve(m)
        end do

        do n=1,3
            ijk2(:) = ijkb(:) - n*s_lr3d(:)
            ijkg(:) = ijkb(:) + n*s_lr3d(:)

            do m=1,nqvst
                qvst(m)%r3d(ijkg(1),ijkg(2),ijkg(3)) = qvst(m)%r3d(ijk2(1),ijk2(2),ijk2(3))
            end do
        end do
    end do
    end do
    end do

end subroutine tur_set_bc_symmetry

subroutine tur_set_bc_farfield(nb,nr)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : fourth,half,one,two
    use mod_constants, only : zero,small,sml_ssf
    use mod_constants, only : mide,m3x3,id_u,id_v,id_w
    use mod_constants, only : nvis_tur_sa
    use mod_constants, only : nvis_tur_sst,nvis_tur_hst
    use mod_constants, only : nscheme_policy_edge2
    use mod_datatypes, only : bc_region_t,fld_array_t
    use mod_variables, only : ntursch,nvis,turconsts
    use mod_variables, only : edvisoo,nuoo,tkeoo,omeoo
    use mod_fieldvars, only : mb_top,mb_sxyz
    use mod_fieldvars, only : npvt,mb_pvt,nqvst,mb_qvst
    implicit none
    integer(kind_int), intent(in) :: nb,nr
    integer(kind_int)          :: i,j,k,m,n,ierr
    integer(kind_int)          :: s_st(3),s_ed(3)
    integer(kind_int)          :: s_nd,s_lr,s_lr3d(3)
    integer(kind_int)          :: m1,m2,m3,ijkg(3),ijkb(3)
    integer(kind_int)          :: ijk2(3),ijk3(3),ijk4(3)
    real(kind_real)            :: vve(id_u:id_w),vv2(id_u:id_w)
    real(kind_real)            :: vv3(id_u:id_w),vv4(id_u:id_w)
    real(kind_real)            :: pve(nqvst),pv2(nqvst)
    real(kind_real)            :: pv3(nqvst),pv4(nqvst)
    real(kind_real)            :: vne,nx,ny,nz,on
    real(kind_real)            :: turoo(nqvst),turlim(nqvst)
    type(bc_region_t), pointer :: reg
    type(fld_array_t), pointer :: sxyz(:),pvt(:),qvst(:)

    select case(nvis)
    case(nvis_tur_sa)
        turoo(1) = nuoo
        turoo(2) = edvisoo

        turlim(1) = turconsts(12)
        turlim(2) = small
    case(nvis_tur_sst,nvis_tur_hst)
        turoo(1) = tkeoo
        turoo(2) = omeoo
        turoo(3) = edvisoo

        turlim(1) = turconsts(13)
        turlim(2) = turconsts(14)
        turlim(3) = small
    end select

    reg     => mb_top(nb)%bcs(nr)
    s_st(:) = reg%s_st(:)
    s_ed(:) = reg%s_ed(:)
    s_nd    = reg%s_nd
    s_lr    = reg%s_lr

    s_lr3d(:) = s_lr*mide(:,s_nd)

    m1 = m3x3(1,s_nd)
    m2 = m3x3(2,s_nd)
    m3 = m3x3(3,s_nd)

    sxyz => mb_sxyz(nb)%fld
    pvt  => mb_pvt(nb)%fld
    qvst => mb_qvst(nb)%fld

    do k=s_st(3),s_ed(3)
    do j=s_st(2),s_ed(2)
    do i=s_st(1),s_ed(1)
        ijkb(:) = (/i,j,k/)
        ijk2(:) = ijkb(:) - s_lr3d(:)
        ijk3(:) = ijk2(:) - s_lr3d(:)
        ijk4(:) = ijk3(:) - s_lr3d(:)

        do m=id_u,id_w
            vv2(m) = pvt(m)%r3d(ijk2(1),ijk2(2),ijk2(3))
            vv3(m) = pvt(m)%r3d(ijk3(1),ijk3(2),ijk3(3))
            vv4(m) = pvt(m)%r3d(ijk4(1),ijk4(2),ijk4(3))
        end do

        if (ntursch == nscheme_policy_edge2) then
            do m=id_u,id_w
                vve(m) = 2.0*vv2(m) - vv3(m)
            end do
        else
            do m=id_u,id_w
                vve(m) = 3.0*vv2(m) - 3.0*vv3(m) + vv4(m)
            end do
        end if

        nx = sxyz(m1)%r3d(i,j,k)
        ny = sxyz(m2)%r3d(i,j,k)
        nz = sxyz(m3)%r3d(i,j,k)
        on = s_lr/max(sqrt(nx*nx + ny*ny + nz*nz),sml_ssf)
        nx = nx*on
        ny = ny*on
        nz = nz*on

        vne = nx*vve(id_u) + ny*vve(id_v) + nz*vve(id_w)

        if (vne > zero) then
            do m=1,nqvst
                pv2(m) = qvst(m)%r3d(ijk2(1),ijk2(2),ijk2(3))
                pv3(m) = qvst(m)%r3d(ijk3(1),ijk3(2),ijk3(3))
                pv4(m) = qvst(m)%r3d(ijk4(1),ijk4(2),ijk4(3))
            end do

            if (ntursch == nscheme_policy_edge2) then
                do m=1,nqvst
                    pve(m) = 2.0*pv2(m) - pv3(m)
                end do
            else
                do m=1,nqvst
                    pve(m) = 3.0*pv2(m) - 3.0*pv3(m) + pv4(m)
                end do
            end if

            do m=1,nqvst
                if ( pve(m) < turlim(m) ) then
                    pve(m) = pv2(m)
                end if
            end do
        else
            do m=1,nqvst
                pve(m) = turoo(m)
            end do
        end if

        do m=1,nqvst
            qvst(m)%r3d(i,j,k) = pve(m)
        end do

        do n=1,3
            ijkg(:) = ijkb(:) + n*s_lr3d(:)

            do m=1,nqvst
                qvst(m)%r3d(ijkg(1),ijkg(2),ijkg(3)) = pve(m)
            end do
        end do
    end do
    end do
    end do

end subroutine tur_set_bc_farfield


subroutine tur_set_bc_inflow(nb,nr)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : mide,nvis_tur_sa
    use mod_constants, only : nvis_tur_sst,nvis_tur_hst
    use mod_datatypes, only : bc_region_t,fld_array_t
    use mod_variables, only : nvis,edvisoo,nuoo,tkeoo,omeoo
    use mod_fieldvars, only : mb_top,nqvst,mb_qvst
    implicit none
    integer(kind_int), intent(in) :: nb,nr
    integer(kind_int)          :: i,j,k,m,n,ierr
    integer(kind_int)          :: s_st(3),s_ed(3)
    integer(kind_int)          :: s_nd,s_lr,s_lr3d(3)
    integer(kind_int)          :: ijkg(3),ijkb(3)
    real(kind_real)            :: turoo(nqvst)
    type(bc_region_t), pointer :: reg
    type(fld_array_t), pointer :: qvst(:)

    select case(nvis)
    case(nvis_tur_sa)
        turoo(1) = nuoo
        turoo(2) = edvisoo
    case(nvis_tur_sst,nvis_tur_hst)
        turoo(1) = tkeoo
        turoo(2) = omeoo
        turoo(3) = edvisoo
    end select

    reg     => mb_top(nb)%bcs(nr)
    s_st(:) = reg%s_st(:)
    s_ed(:) = reg%s_ed(:)
    s_nd    = reg%s_nd
    s_lr    = reg%s_lr

    s_lr3d(:) = s_lr*mide(:,s_nd)

    qvst => mb_qvst(nb)%fld

    do k=s_st(3),s_ed(3)
    do j=s_st(2),s_ed(2)
    do i=s_st(1),s_ed(1)
        ijkb(:) = (/i,j,k/)

        do m=1,nqvst
            qvst(m)%r3d(i,j,k) = turoo(m)
        end do

        do n=1,3
            ijkg(:) = ijkb(:) + n*s_lr3d(:)

            do m=1,nqvst
                qvst(m)%r3d(ijkg(1),ijkg(2),ijkg(3)) = turoo(m)
            end do
        end do
    end do
    end do
    end do

end subroutine tur_set_bc_inflow

subroutine tur_set_bc_outflow(nb,nr)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : one,mide,small
    use mod_constants, only : nvis_tur_sa
    use mod_constants, only : nvis_tur_sst,nvis_tur_hst
    use mod_datatypes, only : bc_region_t,fld_array_t
    use mod_variables, only : nvis,turconsts
    use mod_fieldvars, only : mb_top,nqvst,mb_qvst
    implicit none
    integer(kind_int), intent(in) :: nb,nr
    integer(kind_int)          :: i,j,k,m,n,ierr
    integer(kind_int)          :: s_st(3),s_ed(3)
    integer(kind_int)          :: s_nd,s_lr,s_lr3d(3)
    integer(kind_int)          :: ijkg(3),ijkb(3),ijk2(3),ijk3(3)
    real(kind_real)            :: pve(nqvst),pv2(nqvst),pv3(nqvst)
    real(kind_real)            :: turlim(nqvst)
    type(bc_region_t), pointer :: reg
    type(fld_array_t), pointer :: qvst(:)

    select case(nvis)
    case(nvis_tur_sa)
        turlim(1) = turconsts(12)
        turlim(2) = small
    case(nvis_tur_sst,nvis_tur_hst)
        turlim(1) = turconsts(13)
        turlim(2) = turconsts(14)
        turlim(3) = small
    end select

    reg     => mb_top(nb)%bcs(nr)
    s_st(:) = reg%s_st(:)
    s_ed(:) = reg%s_ed(:)
    s_nd    = reg%s_nd
    s_lr    = reg%s_lr

    s_lr3d(:) = s_lr*mide(:,s_nd)

    qvst => mb_qvst(nb)%fld

    do k=s_st(3),s_ed(3)
    do j=s_st(2),s_ed(2)
    do i=s_st(1),s_ed(1)
        ijkb(:) = (/i,j,k/)
        ijk2(:) = ijkb(:) - s_lr3d(:)
        ijk3(:) = ijk2(:) - s_lr3d(:)

        do m=1,nqvst
            pv2(m) = qvst(m)%r3d(ijk2(1),ijk2(2),ijk2(3))
            pv3(m) = qvst(m)%r3d(ijk3(1),ijk3(2),ijk3(3))
        end do

        do m=1,nqvst
            !!pve(m) = pv2(m)
            pve(m) = 2.0*pv2(m) - pv3(m)

            if ( pve(m) < turlim(m) ) then
                pve(m) = pv2(m)
            end if

            qvst(m)%r3d(i,j,k) = pve(m)
        end do

        do n=1,3
            ijkg(:) = ijkb(:) + n*s_lr3d(:)

            do m=1,nqvst
                qvst(m)%r3d(ijkg(1),ijkg(2),ijkg(3)) = pve(m)
            end do
        end do
    end do
    end do
    end do

end subroutine tur_set_bc_outflow

subroutine tur_set_bc_pole(nb,nr)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,one,two,mide,m3x3
    use mod_constants, only : small,sml_ssf,bc_cut1to1
    use mod_constants, only : nvis_tur_sa
    use mod_constants, only : nvis_tur_sst,nvis_tur_hst
    use mod_datatypes, only : bc_region_t,fld_array_t
    use mod_variables, only : nvis,turconsts
    use mod_fieldvars, only : mb_top,nqvst,mb_qvst,mb_sxyz
    implicit none
    integer(kind_int), intent(in) :: nb,nr
    integer(kind_int)          :: i,j,k,m,n,ierr
    integer(kind_int)          :: ns,m1,m2,m3,npts
    integer(kind_int)          :: nregs,bctype,subtype
    integer(kind_int)          :: npole,naxi,nrot
    integer(kind_int)          :: ist,ied,jst,jed,kst,ked
    integer(kind_int)          :: s_st(3),s_ed(3)
    integer(kind_int)          :: s_nd,s_lr,s_lr3d(3)
    integer(kind_int)          :: ijkg(3),ijkb(3),ijk2(3),ijk3(3)
    real(kind_real)            :: pve(nqvst),pv2(nqvst),pv3(nqvst)
    real(kind_real)            :: turlim(nqvst),sumpv(nqvst)
    type(bc_region_t), pointer :: reg
    type(fld_array_t), pointer :: qvst(:),sxyz(:)

    select case(nvis)
    case(nvis_tur_sa)
        turlim(1) = turconsts(12)
        turlim(2) = small
    case(nvis_tur_sst,nvis_tur_hst)
        turlim(1) = turconsts(13)
        turlim(2) = turconsts(14)
        turlim(3) = small
    end select

    reg     => mb_top(nb)%bcs(nr)
    s_st(:) = reg%s_st(:)
    s_ed(:) = reg%s_ed(:)
    s_nd    = reg%s_nd
    s_lr    = reg%s_lr
    subtype = reg%subtype

    s_lr3d(:) = s_lr*mide(:,s_nd)

    qvst => mb_qvst(nb)%fld
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

            do m=1,nqvst
                pv2(m) = qvst(m)%r3d(ijk2(1),ijk2(2),ijk2(3))
                pv3(m) = qvst(m)%r3d(ijk3(1),ijk3(2),ijk3(3))
            end do

            do m=1,nqvst
                !!pve(m) = pv2(m)
                pve(m) = 2.0*pv2(m) - pv3(m)
                if ( pve(m) < turlim(m) ) then
                    pve(m) = pv2(m)
                end if
            end do

            sumpv(:) = sumpv(:) + pve(:)
        end do

        pve(:) = sumpv(:)/npts

       do k=kst,ked
            ijkb(1) = i*mide(naxi,1) + j*mide(s_nd,1) + k*mide(nrot,1)
            ijkb(2) = i*mide(naxi,2) + j*mide(s_nd,2) + k*mide(nrot,2)
            ijkb(3) = i*mide(naxi,3) + j*mide(s_nd,3) + k*mide(nrot,3)

            do m=1,nqvst
                qvst(m)%r3d(ijkb(1),ijkb(2),ijkb(3)) = pve(m)
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

            do m=1,nqvst
                pv2(m) = qvst(m)%r3d(ijk2(1),ijk2(2),ijk2(3))
                qvst(m)%r3d(ijkg(1),ijkg(2),ijkg(3)) = pv2(m)
            end do
        end do
    end do
    end do
    end do

end subroutine tur_set_bc_pole
