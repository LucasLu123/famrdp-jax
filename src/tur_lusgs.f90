
subroutine tur_lusgs_std_sca
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,nbc_inter_buf_qts
    use mod_constants, only : nsgl_aver_art,nsgl_buffer_qts
    use mod_variables, only : nghnode
    use mod_fieldvars, only : nblkcoms,blkcoms
    use mod_fieldvars, only : neqt,mb_dqt,nqvst,mb_qvst
    use mod_interface, only : pre_exchange_bc_var,post_exchange_bc_var,average_bc_var
    use mod_interface, only : exchange_singulars,average_singulars
    use mod_interface, only : assign_mb_var_uniform,assign_bc_var_nb
    use mod_interface, only : assign_var_via_com_nb
    implicit none
    integer(kind_int) :: nc,nb
    real(kind_real)   :: dq0(1:neqt)

    call tur_boundary_conditions

    call pre_exchange_bc_var(mb_qvst,1,nqvst,nghnode,nbc_inter_buf_qts,nsgl_aver_art)
    call exchange_singulars(mb_qvst,1,nqvst,nsgl_buffer_qts,nsgl_aver_art)
    call post_exchange_bc_var(mb_qvst,1,nqvst,nghnode,nbc_inter_buf_qts,nsgl_aver_art)
    call average_bc_var(mb_qvst,1,nqvst,nbc_inter_buf_qts,nsgl_aver_art)
    call average_singulars(mb_qvst,1,nqvst,nsgl_buffer_qts,nsgl_aver_art)

    call fill_corner_mb_qvst

    call tur_rhside

    dq0(:) = zero
    call assign_mb_var_uniform(mb_dqt,1,neqt,nghnode,dq0)

    do nc=1,nblkcoms
        nb = blkcoms(nc)%nb

        call tur_lusgs_std_sca_foreward(nb)
        call assign_bc_var_nb(nb,mb_dqt,1,neqt,0,dq0)

        call tur_lusgs_std_sca_backward(nb)
        call assign_bc_var_nb(nb,mb_dqt,1,neqt,0,dq0)

        call assign_var_via_com_nb(nb,mb_dqt,1,neqt)
    end do

    call tur_update

end subroutine tur_lusgs_std_sca

subroutine tur_lusgs_std_sca_foreward(nb)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,one,half
    use mod_constants, only : nvis_euler,nsw_dir_close
    use mod_datatypes, only : fld_array_t
    use mod_variables, only : nsw_kdir
    use mod_fieldvars, only : mb_top,mb_sxyz,mb_pvt
    use mod_fieldvars, only : mb_dqt,npvt,neqt,mb_rhst
    use mod_fieldvars, only : mb_rtur,mb_srt
    use mod_interface, only : assign_bc_var_nb
    implicit none
    integer(kind_int), intent(in) :: nb
    integer(kind_int)          :: i,j,k,i1,j1,k1,m,ierr
    integer(kind_int)          :: st(3),ed(3),dijk
    real(kind_real)            :: srvi,srvj,srvk
    real(kind_real)            :: kt,kx,ky,kz,et,ex,ey,ez,ct,cx,cy,cz
    real(kind_real)            :: odia,rhs0(1:neqt)
    real(kind_real)            :: pvi(1:4),pvj(1:4),pvk(1:4)
    real(kind_real)            :: dqi(1:neqt),dqj(1:neqt),dqk(1:neqt)
    real(kind_real)            :: dfi(1:neqt),dfj(1:neqt),dfk(1:neqt)
    type(fld_array_t), pointer :: sxyz(:),pvt(:),rhst(:)
    type(fld_array_t), pointer :: dqt(:),rtur(:),srt(:)

    sxyz => mb_sxyz(nb)%fld
    pvt  => mb_pvt(nb)%fld
    rhst => mb_rhst(nb)%fld
    dqt  => mb_dqt(nb)%fld
    rtur => mb_rtur(nb)%fld
    srt  => mb_srt(nb)%fld

    st(:) = mb_top(nb)%ndst(:)
    ed(:) = mb_top(nb)%nded(:)
    dijk = 1
    do k=st(3),ed(3),dijk
    do j=st(2),ed(2),dijk
    do i=st(1),ed(1),dijk
        i1 = i - 1
        j1 = j - 1
        k1 = k - 1

        do m=1,4
            pvi(m) = pvt(m)%r3d(i1,j,k)
            pvj(m) = pvt(m)%r3d(i,j1,k)
            pvk(m) = pvt(m)%r3d(i,j,k1)
        end do

        do m=1,neqt
            dqi(m) = dqt(m)%r3d(i1,j,k)
            dqj(m) = dqt(m)%r3d(i,j1,k)
            dqk(m) = dqt(m)%r3d(i,j,k1)
        end do

        kx = sxyz(1)%r3d(i1,j,k)
        ky = sxyz(2)%r3d(i1,j,k)
        kz = sxyz(3)%r3d(i1,j,k)
        call tur_mxdq_std(1,4,pvi,kt,kx,ky,kz,1,neqt,dqi,dfi,one)

        ex = sxyz(4)%r3d(i,j1,k)
        ey = sxyz(5)%r3d(i,j1,k)
        ez = sxyz(6)%r3d(i,j1,k)
        call tur_mxdq_std(1,4,pvj,et,ex,ey,ez,1,neqt,dqj,dfj,one)

        if (nsw_kdir == nsw_dir_close) then
            dfk(:) = zero
        else
            cx = sxyz(7)%r3d(i,j,k1)
            cy = sxyz(8)%r3d(i,j,k1)
            cz = sxyz(9)%r3d(i,j,k1)
            call tur_mxdq_std(1,4,pvk,ct,cx,cy,cz,1,neqt,dqk,dfk,one)
        end if

        rhs0(:) = dfi(:) + dfj(:) + dfk(:)

        srvi = srt(1)%r3d(i1,j,k)
        srvj = srt(2)%r3d(i,j1,k)
        srvk = srt(3)%r3d(i,j,k1)
        do m=1,neqt
            rhs0(m) = rhs0(m) + &
                      half*(srvi*dqi(m) + srvj*dqj(m) + srvk*dqk(m))
        end do

        do m=1,neqt
            odia = one/rtur(m)%r3d(i,j,k)
            dqt(m)%r3d(i,j,k) = (rhst(m)%r3d(i,j,k) + rhs0(m))*odia
        end do
    end do
    end do
    end do

end subroutine tur_lusgs_std_sca_foreward


subroutine tur_lusgs_std_sca_backward(nb)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,one,half
    use mod_constants, only : nvis_euler,nsw_dir_close
    use mod_datatypes, only : fld_array_t
    use mod_variables, only : nsw_kdir
    use mod_fieldvars, only : mb_top,mb_sxyz,mb_pvt
    use mod_fieldvars, only : mb_dqt,npvt,neqt,mb_rhst
    use mod_fieldvars, only : mb_rtur,mb_srt
    use mod_interface, only : assign_bc_var_nb
    implicit none
    integer(kind_int), intent(in) :: nb
    integer(kind_int)          :: i,j,k,i1,j1,k1,m,ierr
    integer(kind_int)          :: st(3),ed(3),dijk
    real(kind_real)            :: srvi,srvj,srvk
    real(kind_real)            :: kt,kx,ky,kz,et,ex,ey,ez,ct,cx,cy,cz
    real(kind_real)            :: odia,rhs0(1:neqt)
    real(kind_real)            :: pvi(1:4),pvj(1:4),pvk(1:4)
    real(kind_real)            :: dqi(1:neqt),dqj(1:neqt),dqk(1:neqt)
    real(kind_real)            :: dfi(1:neqt),dfj(1:neqt),dfk(1:neqt)
    type(fld_array_t), pointer :: sxyz(:),pvt(:),rhst(:)
    type(fld_array_t), pointer :: dqt(:),rtur(:),srt(:)

    sxyz => mb_sxyz(nb)%fld
    pvt  => mb_pvt(nb)%fld
    rhst => mb_rhst(nb)%fld
    dqt  => mb_dqt(nb)%fld
    rtur => mb_rtur(nb)%fld
    srt  => mb_srt(nb)%fld

    st(:) = mb_top(nb)%nded(:)
    ed(:) = mb_top(nb)%ndst(:)
    dijk = -1
    do k=st(3),ed(3),dijk
    do j=st(2),ed(2),dijk
    do i=st(1),ed(1),dijk
        i1 = i + 1
        j1 = j + 1
        k1 = k + 1

        do m=1,4
            pvi(m) = pvt(m)%r3d(i1,j,k)
            pvj(m) = pvt(m)%r3d(i,j1,k)
            pvk(m) = pvt(m)%r3d(i,j,k1)
        end do

        do m=1,neqt
            dqi(m) = dqt(m)%r3d(i1,j,k)
            dqj(m) = dqt(m)%r3d(i,j1,k)
            dqk(m) = dqt(m)%r3d(i,j,k1)
        end do

        kx = sxyz(1)%r3d(i1,j,k)
        ky = sxyz(2)%r3d(i1,j,k)
        kz = sxyz(3)%r3d(i1,j,k)
        call tur_mxdq_std(1,4,pvi,kt,kx,ky,kz,1,neqt,dqi,dfi,-one)

        ex = sxyz(4)%r3d(i,j1,k)
        ey = sxyz(5)%r3d(i,j1,k)
        ez = sxyz(6)%r3d(i,j1,k)
        call tur_mxdq_std(1,4,pvj,et,ex,ey,ez,1,neqt,dqj,dfj,-one)

        if (nsw_kdir == nsw_dir_close) then
            dfk(:) = zero
        else
            cx = sxyz(7)%r3d(i,j,k1)
            cy = sxyz(8)%r3d(i,j,k1)
            cz = sxyz(9)%r3d(i,j,k1)
            call tur_mxdq_std(1,4,pvk,ct,cx,cy,cz,1,neqt,dqk,dfk,-one)
        end if

        rhs0(:) = dfi(:) + dfj(:) + dfk(:)

        srvi = srt(1)%r3d(i1,j,k)
        srvj = srt(2)%r3d(i,j1,k)
        srvk = srt(3)%r3d(i,j,k1)
        do m=1,neqt
            rhs0(m) = rhs0(m) - &
                      half*(srvi*dqi(m) + srvj*dqj(m) + srvk*dqk(m))
        end do

        do m=1,neqt
            odia = one/rtur(m)%r3d(i,j,k)
            dqt(m)%r3d(i,j,k) = dqt(m)%r3d(i,j,k) - rhs0(m)*odia
        end do
    end do
    end do
    end do

end subroutine tur_lusgs_std_sca_backward


subroutine tur_mxdq_std(nsp,nep,prim,nt,nx,ny,nz,nsf,nef,dq,df,fsw)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : one,half,two,sml_ssf
    use mod_variables, only : gamma
    implicit none
    integer(kind_int), intent(in)  :: nsp,nep
    real(kind_real),   intent(in)  :: prim(nsp:nep)
    real(kind_real),   intent(in)  :: nt,nx,ny,nz
    integer(kind_int), intent(in)  :: nsf,nef
    real(kind_real),   intent(in)  :: dq(nsf:nef)
    real(kind_real),   intent(out) :: df(nsf:nef)
    real(kind_real),   intent(in)  :: fsw
    integer(kind_int) :: m
    real(kind_real)   :: ro,vx,vy,vz,vn,vna

    ro = prim(nsp)
    vx = prim(nsp+1)
    vy = prim(nsp+2)
    vz = prim(nsp+3)

    vn = nx*vx + ny*vy + nz*vz

    vna = half*(vn + fsw*abs(vn))
    do m=nsf,nef
        df(m) = vna*dq(m)
    end do

end subroutine tur_mxdq_std

subroutine tur_prsgs_std_sca
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,nsgl_aver_art
    use mod_constants, only : nsgl_buffer_qts,nsgl_buffer_dqt
    use mod_constants, only : nbc_inter_buf_qts,nbc_inter_buf_dqt
    use mod_constants, only : nsweep_foreward,nsweep_backward
    use mod_variables, only : nghnode,ntursub
    use mod_fieldvars, only : nblkcoms,blkcoms
    use mod_fieldvars, only : neqt,mb_dqt,nqvst,mb_qvst
    use mod_interface, only : pre_exchange_bc_var,post_exchange_bc_var,average_bc_var
    use mod_interface, only : exchange_singulars,average_singulars
    use mod_interface, only : assign_mb_var_uniform,assign_bc_var_nb
    use mod_interface, only : assign_var_via_com_nb
    implicit none
    integer(kind_int) :: iter,nc,nb,m
    real(kind_real)   :: dq0(1:neqt)
    external          :: tur_prsgs_std_sca_pre

    call tur_boundary_conditions

    call pre_exchange_bc_var(mb_qvst,1,nqvst,nghnode,nbc_inter_buf_qts,nsgl_aver_art)
    call exchange_singulars(mb_qvst,1,nqvst,nsgl_buffer_qts,nsgl_aver_art)
    call post_exchange_bc_var(mb_qvst,1,nqvst,nghnode,nbc_inter_buf_qts,nsgl_aver_art)
    call average_bc_var(mb_qvst,1,nqvst,nbc_inter_buf_qts,nsgl_aver_art)
    call average_singulars(mb_qvst,1,nqvst,nsgl_buffer_qts,nsgl_aver_art)

    call fill_corner_mb_qvst


    call tur_rhside


    dq0(:) = zero
    call assign_mb_var_uniform(mb_dqt,1,neqt,nghnode,dq0)

    call run_on_blkcoms(tur_prsgs_std_sca_pre)

    call pre_exchange_bc_var(mb_dqt,1,neqt,1,nbc_inter_buf_dqt,nsgl_aver_art)
    call exchange_singulars(mb_dqt,1,neqt,nsgl_buffer_dqt,nsgl_aver_art)
    call post_exchange_bc_var(mb_dqt,1,neqt,1,nbc_inter_buf_dqt,nsgl_aver_art)
    call average_bc_var(mb_dqt,1,neqt,nbc_inter_buf_dqt,nsgl_aver_art)
    call average_singulars(mb_dqt,1,neqt,nsgl_buffer_dqt,nsgl_aver_art)

    do iter=1,ntursub
        do nc=1,nblkcoms
            nb = blkcoms(nc)%nb

            call tur_prsgs_std_sca_sweep(nb,nsweep_foreward)
            call assign_bc_var_nb(nb,mb_dqt,1,neqt,0,dq0)

            call tur_prsgs_std_sca_sweep(nb,nsweep_backward)
            call assign_bc_var_nb(nb,mb_dqt,1,neqt,0,dq0)

            call assign_var_via_com_nb(nb,mb_dqt,1,neqt)
        end do

        call pre_exchange_bc_var(mb_dqt,1,neqt,1,nbc_inter_buf_dqt,nsgl_aver_art)
        call exchange_singulars(mb_dqt,1,neqt,nsgl_buffer_dqt,nsgl_aver_art)
        call post_exchange_bc_var(mb_dqt,1,neqt,1,nbc_inter_buf_dqt,nsgl_aver_art)
        call average_bc_var(mb_dqt,1,neqt,nbc_inter_buf_dqt,nsgl_aver_art)
        call average_singulars(mb_dqt,1,neqt,nsgl_buffer_dqt,nsgl_aver_art)
    end do

    call tur_update

end subroutine tur_prsgs_std_sca

subroutine tur_prsgs_std_sca_pre(nb)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : one
    use mod_datatypes, only : fld_array_t
    use mod_fieldvars, only : mb_top,mb_sxyz
    use mod_fieldvars, only : mb_dqt,neqt,mb_rhst
    use mod_fieldvars, only : mb_rtur
    implicit none
    integer(kind_int), intent(in) :: nb
    integer(kind_int)          :: i,j,k,m
    integer(kind_int)          :: st(3),ed(3)
    real(kind_real)            :: odia
    type(fld_array_t), pointer :: rhst(:),dqt(:),rtur(:)

    rhst => mb_rhst(nb)%fld
    dqt  => mb_dqt(nb)%fld
    rtur => mb_rtur(nb)%fld

    st(:) = mb_top(nb)%ndst(:)
    ed(:) = mb_top(nb)%nded(:)

    do m=1,neqt
        do k=st(3),ed(3)
        do j=st(2),ed(2)
        do i=st(1),ed(1)
            odia = one/rtur(m)%r3d(i,j,k)
            dqt(m)%r3d(i,j,k) = rhst(m)%r3d(i,j,k)*odia
        end do
        end do
        end do
    end do

end subroutine tur_prsgs_std_sca_pre

subroutine tur_prsgs_std_sca_sweep(nb,swp)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,one,half
    use mod_constants, only : nvis_euler,nsw_dir_close
    use mod_datatypes, only : fld_array_t
    use mod_variables, only : nvis,nsw_kdir
    use mod_fieldvars, only : mb_top,mb_sxyz,mb_pvt
    use mod_fieldvars, only : mb_dqt,npvt,neqt,mb_rhst
    use mod_fieldvars, only : mb_srt,mb_rtur
    implicit none
    integer(kind_int), intent(in) :: nb,swp
    integer(kind_int)          :: i,j,k,i1,j1,k1,m,ierr
    integer(kind_int)          :: st(3),ed(3),dijk
    real(kind_real)            :: srvi,srvj,srvk
    real(kind_real)            :: kt,kx,ky,kz,et,ex,ey,ez,ct,cx,cy,cz
    real(kind_real)            :: odia,rhs0(1:neqt)
    real(kind_real)            :: pvi(1:4),pvj(1:4),pvk(1:4)
    real(kind_real)            :: dqi(1:neqt),dqj(1:neqt),dqk(1:neqt)
    real(kind_real)            :: dfi(1:neqt),dfj(1:neqt),dfk(1:neqt)
    type(fld_array_t), pointer :: sxyz(:),pvt(:),rhst(:)
    type(fld_array_t), pointer :: dqt(:),srt(:),rtur(:)

    sxyz => mb_sxyz(nb)%fld
    pvt  => mb_pvt(nb)%fld
    rhst => mb_rhst(nb)%fld
    dqt  => mb_dqt(nb)%fld
    srt  => mb_srt(nb)%fld
    rtur => mb_rtur(nb)%fld

    if (swp > 0) then
        st(:) = mb_top(nb)%ndst(:)
        ed(:) = mb_top(nb)%nded(:)
        dijk  = 1
    else
        st(:) = mb_top(nb)%nded(:)
        ed(:) = mb_top(nb)%ndst(:)
        dijk  = -1
    end if

    do k=st(3),ed(3),dijk
    do j=st(2),ed(2),dijk
    do i=st(1),ed(1),dijk
        i1 = i - 1
        j1 = j - 1
        k1 = k - 1

        do m=1,4
            pvi(m) = pvt(m)%r3d(i1,j,k)
            pvj(m) = pvt(m)%r3d(i,j1,k)
            pvk(m) = pvt(m)%r3d(i,j,k1)
        end do

        do m=1,neqt
            dqi(m) = dqt(m)%r3d(i1,j,k)
            dqj(m) = dqt(m)%r3d(i,j1,k)
            dqk(m) = dqt(m)%r3d(i,j,k1)
        end do

        kx = sxyz(1)%r3d(i1,j,k)
        ky = sxyz(2)%r3d(i1,j,k)
        kz = sxyz(3)%r3d(i1,j,k)
        call tur_mxdq_std(1,4,pvi,kt,kx,ky,kz,1,neqt,dqi,dfi,one)

        ex = sxyz(4)%r3d(i,j1,k)
        ey = sxyz(5)%r3d(i,j1,k)
        ez = sxyz(6)%r3d(i,j1,k)
        call tur_mxdq_std(1,4,pvj,et,ex,ey,ez,1,neqt,dqj,dfj,one)

        if (nsw_kdir == nsw_dir_close) then
            dfk(:) = zero
        else
            cx = sxyz(7)%r3d(i,j,k1)
            cy = sxyz(8)%r3d(i,j,k1)
            cz = sxyz(9)%r3d(i,j,k1)
            call tur_mxdq_std(1,4,pvk,ct,cx,cy,cz,1,neqt,dqk,dfk,one)
        end if

        rhs0(:) = dfi(:) + dfj(:) + dfk(:)

        srvi = srt(1)%r3d(i1,j,k)
        srvj = srt(2)%r3d(i,j1,k)
        srvk = srt(3)%r3d(i,j,k1)
        do m=1,neqt
            rhs0(m) = rhs0(m) + &
                      half*(srvi*dqi(m) + srvj*dqj(m) + srvk*dqk(m))
        end do


        i1 = i + 1
        j1 = j + 1
        k1 = k + 1

        do m=1,4
            pvi(m) = pvt(m)%r3d(i1,j,k)
            pvj(m) = pvt(m)%r3d(i,j1,k)
            pvk(m) = pvt(m)%r3d(i,j,k1)
        end do

        do m=1,neqt
            dqi(m) = dqt(m)%r3d(i1,j,k)
            dqj(m) = dqt(m)%r3d(i,j1,k)
            dqk(m) = dqt(m)%r3d(i,j,k1)
        end do

        kx = sxyz(1)%r3d(i1,j,k)
        ky = sxyz(2)%r3d(i1,j,k)
        kz = sxyz(3)%r3d(i1,j,k)
        call tur_mxdq_std(1,4,pvi,kt,kx,ky,kz,1,neqt,dqi,dfi,-one)

        ex = sxyz(4)%r3d(i,j1,k)
        ey = sxyz(5)%r3d(i,j1,k)
        ez = sxyz(6)%r3d(i,j1,k)
        call tur_mxdq_std(1,4,pvj,et,ex,ey,ez,1,neqt,dqj,dfj,-one)

        if (nsw_kdir == nsw_dir_close) then
            dfk(:) = zero
        else
            cx = sxyz(7)%r3d(i,j,k1)
            cy = sxyz(8)%r3d(i,j,k1)
            cz = sxyz(9)%r3d(i,j,k1)
            call tur_mxdq_std(1,4,pvk,ct,cx,cy,cz,1,neqt,dqk,dfk,-one)
        end if

        rhs0(:) = rhs0(:) - (dfi(:) + dfj(:) + dfk(:))

        srvi = srt(1)%r3d(i1,j,k)
        srvj = srt(2)%r3d(i,j1,k)
        srvk = srt(3)%r3d(i,j,k1)
        do m=1,neqt
            rhs0(m) = rhs0(m) + &
                      half*(srvi*dqi(m) + srvj*dqj(m) + srvk*dqk(m))
        end do

        do m=1,neqt
            odia = one/rtur(m)%r3d(i,j,k)
            dqt(m)%r3d(i,j,k) = (rhst(m)%r3d(i,j,k) + rhs0(m))*odia
        end do
    end do
    end do
    end do

end subroutine tur_prsgs_std_sca_sweep

