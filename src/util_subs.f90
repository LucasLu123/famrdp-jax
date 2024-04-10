
    
subroutine assign_mb_var_uniform(mb_var,nst,ned,ngh,var)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_datatypes, only : var_block_t
    use mod_fieldvars, only : nblkcoms,blkcoms
    use mod_openmp
    implicit none
    type(var_block_t), pointer, intent(in) :: mb_var(:)
    integer(kind_int)         , intent(in) :: nst,ned,ngh
    real(kind_real)           , intent(in) :: var(nst:ned)
    integer(kind_int) :: nc,nb,i,j,k,m,ierr
    integer(kind_int) :: st(3),ed(3)

    do nc=1,nblkcoms
        nb = blkcoms(nc)%nb

        st(:) = 1 - ngh
        ed(:) = blkcoms(nc)%top%nijk(:) + ngh
        do m=nst,ned
!$OMP parallel do
            do k=st(3),ed(3)
            do j=st(2),ed(2)
            do i=st(1),ed(1)
                mb_var(nb)%fld(m)%r3d(i,j,k) = var(m)
            end do
            end do
            end do
!$OMP end parallel do
        end do
    end do

end subroutine assign_mb_var_uniform

subroutine assign_mb_var_uniform_sp(mb_var,nst,ned,ngh,var)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_datatypes, only : var_block_t
    use mod_fieldvars, only : nblkcoms,blkcomssp
    use mod_openmp
    implicit none
    type(var_block_t), pointer, intent(in) :: mb_var(:)
    integer(kind_int)         , intent(in) :: nst,ned,ngh
    real(kind_real)           , intent(in) :: var(nst:ned)
    integer(kind_int) :: nc,nb,i,j,k,m,ierr
    integer(kind_int) :: st(3),ed(3)

    do nc=1,nblkcoms
        nb = blkcomssp(nc)%nb

        st(:) = 1 - ngh
        ed(:) = blkcomssp(nc)%top%nijk(:) + ngh
        do m=nst,ned
!$OMP parallel do
            do k=st(3),ed(3)
            do j=st(2),ed(2)
            do i=st(1),ed(1)
                mb_var(nb)%fld(m)%r3d(i,j,k) = var(m)
            end do
            end do
            end do
!$OMP end parallel do
        end do
    end do

end subroutine assign_mb_var_uniform_sp

subroutine assign_bc_var_nb(nb,mb_var,nst,ned,ngh,var)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : bc_cut1to1
    use mod_datatypes, only : var_block_t,bc_region_t
    use mod_fieldvars, only : mb_top
    use mod_openmp
    implicit none
    integer(kind_int)         , intent(in) :: nb
    type(var_block_t), pointer, intent(in) :: mb_var(:)
    integer(kind_int)         , intent(in) :: nst,ned,ngh
    real(kind_real)           , intent(in) :: var(nst:ned)
    integer(kind_int)          :: nr,i,j,k,m,ierr
    integer(kind_int)          :: s_st(3),s_ed(3),st(3),ed(3)
    integer(kind_int)          :: nregs,s_nd,s_lr,bctype
    type(bc_region_t), pointer :: reg

    nregs = mb_top(nb)%nregions
    do nr=1,nregs
        reg => mb_top(nb)%bcs(nr)

        bctype  = reg%bctype

        if (bctype /= bc_cut1to1) then
            s_st(:) = reg%s_st(:)
            s_ed(:) = reg%s_ed(:)
            s_nd    = reg%s_nd
            s_lr    = reg%s_lr

            call bc_extend_outward(s_st,s_ed,s_nd,s_lr,ngh,st,ed)

            do m=nst,ned
!OMP parallel do
                do k=st(3),ed(3)
                do j=st(2),ed(2)
                do i=st(1),ed(1)
                    mb_var(nb)%fld(m)%r3d(i,j,k) = var(m)
                end do
                end do
                end do
!OMP end parallel do
            end do
        end if
    end do

end subroutine assign_bc_var_nb

subroutine assign_bc_var_nb_sp(nb,mb_var,nst,ned,ngh,var)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : bc_cut1to1
    use mod_datatypes, only : var_block_t,bc_region_t
    use mod_fieldvars, only : mb_topsp
    use mod_openmp
    implicit none
    integer(kind_int)         , intent(in) :: nb
    type(var_block_t), pointer, intent(in) :: mb_var(:)
    integer(kind_int)         , intent(in) :: nst,ned,ngh
    real(kind_real)           , intent(in) :: var(nst:ned)
    integer(kind_int)          :: nr,i,j,k,m,ierr
    integer(kind_int)          :: s_st(3),s_ed(3),st(3),ed(3)
    integer(kind_int)          :: nregs,s_nd,s_lr,bctype
    type(bc_region_t), pointer :: reg

    nregs = mb_topsp(nb)%nregions
    do nr=1,nregs
        reg => mb_topsp(nb)%bcs(nr)

        bctype  = reg%bctype

        if (bctype /= bc_cut1to1) then
            s_st(:) = reg%s_st(:)
            s_ed(:) = reg%s_ed(:)
            s_nd    = reg%s_nd
            s_lr    = reg%s_lr

            call bc_extend_outward(s_st,s_ed,s_nd,s_lr,ngh,st,ed)

            do m=nst,ned
!OMP parallel do
                do k=st(3),ed(3)
                do j=st(2),ed(2)
                do i=st(1),ed(1)
                    mb_var(nb)%fld(m)%r3d(i,j,k) = var(m)
                end do
                end do
                end do
!OMP end parallel do
            end do
        end if
    end do

end subroutine assign_bc_var_nb_sp

subroutine assign_bc_var_uniform(mb_var,nst,ned,ngh,var)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_datatypes, only : var_block_t
    use mod_fieldvars, only : nblkcoms,blkcoms
    use mod_interface, only : assign_bc_var_nb
    implicit none
    type(var_block_t), pointer, intent(in) :: mb_var(:)
    integer(kind_int)         , intent(in) :: nst,ned,ngh
    real(kind_real)           , intent(in) :: var(nst:ned)
    integer(kind_int) :: nc,nb

    do nc=1,nblkcoms
        nb  =  blkcoms(nc)%nb

        call assign_bc_var_nb(nb,mb_var,nst,ned,ngh,var)
    end do

end subroutine assign_bc_var_uniform

subroutine assign_bc_var_uniform_sp(mb_var,nst,ned,ngh,var)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_datatypes, only : var_block_t
    use mod_fieldvars, only : nblkcoms,blkcomssp
    use mod_interface, only : assign_bc_var_nb_sp
    implicit none
    type(var_block_t), pointer, intent(in) :: mb_var(:)
    integer(kind_int)         , intent(in) :: nst,ned,ngh
    real(kind_real)           , intent(in) :: var(nst:ned)
    integer(kind_int) :: nc,nb

    do nc=1,nblkcoms
        nb  =  blkcomssp(nc)%nb

        call assign_bc_var_nb_sp(nb,mb_var,nst,ned,ngh,var)
    end do

end subroutine assign_bc_var_uniform_sp

subroutine ghost_bc_var_nb(nb,mb_var,nst,ned)
    use mod_kndconsts, only : kind_int
    use mod_constants, only : mide,bc_cut1to1
    use mod_datatypes, only : var_block_t
    use mod_datatypes, only : bc_region_t,fld_array_t
    use mod_fieldvars, only : mb_top
    use mod_interface, only : fill_corner_var_nb
    use mod_openmp
    implicit none
    integer(kind_int)         , intent(in) :: nb
    type(var_block_t), pointer, intent(in) :: mb_var(:)
    integer(kind_int)         , intent(in) :: nst,ned
    integer(kind_int)          :: nr,i,j,k,m,ierr
    integer(kind_int)          :: s_st(3),s_ed(3),s_lr3d(3)
    integer(kind_int)          :: ijkb(3),ijkg(3),ig,jg,kg
    integer(kind_int)          :: nregs,s_nd,s_lr,bctype
    type(bc_region_t), pointer :: reg
    type(fld_array_t), pointer :: var(:)

    var => mb_var(nb)%fld

    nregs = mb_top(nb)%nregions
    do nr=1,nregs
        reg => mb_top(nb)%bcs(nr)

        bctype  = reg%bctype
        if (bctype /= bc_cut1to1) then
            s_st(:) = reg%s_st(:)
            s_ed(:) = reg%s_ed(:)
            s_nd    = reg%s_nd
            s_lr    = reg%s_lr

            s_lr3d(:) = s_lr*mide(:,s_nd)
!$OMP parallel do private(ijkg,ijkb,ig,jg,kg,m)
            do k=s_st(3),s_ed(3)
            do j=s_st(2),s_ed(2)
            do i=s_st(1),s_ed(1)
                ijkb(:) = (/i,j,k/)
                ijkg(:) = ijkb(:) + s_lr3d(:)

                ig = ijkg(1)
                jg = ijkg(2)
                kg = ijkg(3)

                do m=nst,ned
                   var(m)%r3d(ig,jg,kg) = var(m)%r3d(i,j,k)
                end do
            end do
            end do
            end do
!$OMP end parallel do
        end if
    end do

    call fill_corner_var_nb(nb,mb_var,nst,ned)

end subroutine ghost_bc_var_nb

subroutine ghost_bc_var_nb_cc(nb,mb_var,nst,ned)
    use mod_kndconsts, only : kind_int
    use mod_constants, only : mide,bc_cut1to1
    use mod_datatypes, only : var_block_t
    use mod_datatypes, only : bc_region_t,fld_array_t
    use mod_fieldvars, only : mb_topc
    use mod_interface, only : fill_corner_var_nb_cc
    use mod_openmp
    implicit none
    integer(kind_int)         , intent(in) :: nb
    type(var_block_t), pointer, intent(in) :: mb_var(:)
    integer(kind_int)         , intent(in) :: nst,ned
    integer(kind_int)          :: nr,i,j,k,m,ierr
    integer(kind_int)          :: s_st(3),s_ed(3),s_lr3d(3)
    integer(kind_int)          :: ijkb(3),ijkg(3),ig,jg,kg
    integer(kind_int)          :: nregs,s_nd,s_lr,bctype
    type(bc_region_t), pointer :: reg
    type(fld_array_t), pointer :: var(:)

    var => mb_var(nb)%fld

    nregs = mb_topc(nb)%nregions
    do nr=1,nregs
        reg => mb_topc(nb)%bcs(nr)

        bctype  = reg%bctype
        if (bctype /= bc_cut1to1) then
            s_st(:) = reg%s_st(:)
            s_ed(:) = reg%s_ed(:)
            s_nd    = reg%s_nd
            s_lr    = reg%s_lr

            s_lr3d(:) = s_lr*mide(:,s_nd)
!$OMP parallel do private(ijkg,ijkb,ig,jg,kg,m)
            do k=s_st(3),s_ed(3)
            do j=s_st(2),s_ed(2)
            do i=s_st(1),s_ed(1)
                ijkb(:) = (/i,j,k/)
                ijkg(:) = ijkb(:) + s_lr3d(:)

                ig = ijkg(1)
                jg = ijkg(2)
                kg = ijkg(3)

                do m=nst,ned
                   var(m)%r3d(ig,jg,kg) = var(m)%r3d(i,j,k)
                end do
            end do
            end do
            end do
!$OMP end parallel do
        end if
    end do

    call fill_corner_var_nb_cc(nb,mb_var,nst,ned)

end subroutine ghost_bc_var_nb_cc

subroutine ghost_bc_var_nb_sp(nb,mb_var,nst,ned)
    use mod_kndconsts, only : kind_int
    use mod_constants, only : mide,bc_cut1to1
    use mod_datatypes, only : var_block_t
    use mod_datatypes, only : bc_region_t,fld_array_t
    use mod_fieldvars, only : mb_topsp
    use mod_interface, only : fill_corner_var_nb_sp
    use mod_openmp
    implicit none
    integer(kind_int)         , intent(in) :: nb
    type(var_block_t), pointer, intent(in) :: mb_var(:)
    integer(kind_int)         , intent(in) :: nst,ned
    integer(kind_int)          :: nr,i,j,k,m,ierr
    integer(kind_int)          :: s_st(3),s_ed(3),s_lr3d(3)
    integer(kind_int)          :: ijkb(3),ijkg(3),ig,jg,kg
    integer(kind_int)          :: nregs,s_nd,s_lr,bctype
    type(bc_region_t), pointer :: reg
    type(fld_array_t), pointer :: var(:)

    var => mb_var(nb)%fld

    nregs = mb_topsp(nb)%nregions
    do nr=1,nregs
        reg => mb_topsp(nb)%bcs(nr)

        bctype  = reg%bctype
        if (bctype /= bc_cut1to1) then
            s_st(:) = reg%s_st(:)
            s_ed(:) = reg%s_ed(:)
            s_nd    = reg%s_nd
            s_lr    = reg%s_lr

            s_lr3d(:) = s_lr*mide(:,s_nd)
!$OMP parallel do private(ijkg,ijkb,ig,jg,kg,m)
            do k=s_st(3),s_ed(3)
            do j=s_st(2),s_ed(2)
            do i=s_st(1),s_ed(1)
                ijkb(:) = (/i,j,k/)
                ijkg(:) = ijkb(:) + s_lr3d(:)

                ig = ijkg(1)
                jg = ijkg(2)
                kg = ijkg(3)

                do m=nst,ned
                   var(m)%r3d(ig,jg,kg) = var(m)%r3d(i,j,k)
                end do
            end do
            end do
            end do
!$OMP end parallel do
        end if
    end do

    call fill_corner_var_nb_sp(nb,mb_var,nst,ned)

end subroutine ghost_bc_var_nb_sp

subroutine fill_corner_var_nb(nb,mb_var,nst,ned)
    use mod_kndconsts, only : kind_int
    use mod_constants, only : third
    use mod_datatypes, only : var_block_t,fld_array_t
    use mod_fieldvars, only : mb_top
    use mod_openmp
    implicit none
    integer(kind_int)         , intent(in) :: nb
    type(var_block_t), pointer, intent(in) :: mb_var(:)
    integer(kind_int)         , intent(in) :: nst,ned
    integer(kind_int)          :: i,j,k,m,ierr
    integer(kind_int)          :: ig,jg,kg,ni,nj,nk
    type(fld_array_t), pointer :: var(:)

    var => mb_var(nb)%fld

    ni = mb_top(nb)%nijk(1)
    nj = mb_top(nb)%nijk(2)
    nk = mb_top(nb)%nijk(3)
!$OMP parallel do private(ig,jg,kg)
    do m=nst,ned
        do j=0,nj+1,nj+1
        do i=0,ni+1,ni+1
            jg = minabs(1,nj,j)
            ig = minabs(1,ni,i)
            do k=1,nk
                var(m)%r3d(i,j,k) = -var(m)%r3d(ig,jg,k) + &
                                     var(m)%r3d(i ,jg,k) + &
                                     var(m)%r3d(ig,j ,k)
            end do
        end do
        end do

        do k=0,nk+1,nk+1
        do i=0,ni+1,ni+1
            kg = minabs(1,nk,k)
            ig = minabs(1,ni,i)
            do j=1,nj
                var(m)%r3d(i,j,k) = -var(m)%r3d(ig,j,kg) + &
                                     var(m)%r3d(i ,j,kg) + &
                                     var(m)%r3d(ig,j,k )
            end do
        end do
        end do

        do k=0,nk+1,nk+1
        do j=0,nj+1,nj+1
            kg = minabs(1,nk,k)
            jg = minabs(1,nj,j)
            do i=1,ni
                var(m)%r3d(i,j,k) = -var(m)%r3d(i,jg,kg) + &
                                     var(m)%r3d(i,j ,kg) + &
                                     var(m)%r3d(i,jg,k )
            end do
        end do
        end do

        do k=0,nk+1,nk+1
        do j=0,nj+1,nj+1
        do i=0,ni+1,ni+1
            kg = minabs(1,nk,k)
            jg = minabs(1,nj,j)
            ig = minabs(1,ni,i)
            var(m)%r3d(i,j,k) = third*( var(m)%r3d(ig,j,k) + &
                                        var(m)%r3d(i,jg,k) + &
                                        var(m)%r3d(i,j,kg) )
        end do
        end do
        end do
    end do
!$OMP end parallel do

    contains

    function minabs(a,b,c) result(d)
        use mod_kndconsts, only : kind_int
        implicit none
        integer(kind_int) :: a,b,c
        integer(kind_int) :: d
        integer(kind_int) :: d1,d2

        d1 = abs(a-c)
        d2 = abs(b-c)
        if (d1 == d2) then
            call error_check(1, &
                             "The distance of C-A and C-B is equal in function minabs")
        else if (d1 > d2) then
            d = b
        else
            d = a
        end if

    end function minabs

end subroutine fill_corner_var_nb

subroutine fill_corner_var_nb_cc(nb,mb_var,nst,ned)
    use mod_kndconsts, only : kind_int
    use mod_constants, only : third
    use mod_datatypes, only : var_block_t,fld_array_t
    use mod_fieldvars, only : mb_topc
    use mod_openmp
    implicit none
    integer(kind_int)         , intent(in) :: nb
    type(var_block_t), pointer, intent(in) :: mb_var(:)
    integer(kind_int)         , intent(in) :: nst,ned
    integer(kind_int)          :: i,j,k,m,ierr
    integer(kind_int)          :: ig,jg,kg,ni,nj,nk
    type(fld_array_t), pointer :: var(:)

    var => mb_var(nb)%fld

    ni = mb_topc(nb)%nijk(1)
    nj = mb_topc(nb)%nijk(2)
    nk = mb_topc(nb)%nijk(3)
!$OMP parallel do private(ig,jg,kg)
    do m=nst,ned
        do j=0,nj+1,nj+1
        do i=0,ni+1,ni+1
            jg = minabscc(1,nj,j)
            ig = minabscc(1,ni,i)
            do k=1,nk
                var(m)%r3d(i,j,k) = -var(m)%r3d(ig,jg,k) + &
                                     var(m)%r3d(i ,jg,k) + &
                                     var(m)%r3d(ig,j ,k)
            end do
        end do
        end do

        do k=0,nk+1,nk+1
        do i=0,ni+1,ni+1
            kg = minabscc(1,nk,k)
            ig = minabscc(1,ni,i)
            do j=1,nj
                var(m)%r3d(i,j,k) = -var(m)%r3d(ig,j,kg) + &
                                     var(m)%r3d(i ,j,kg) + &
                                     var(m)%r3d(ig,j,k )
            end do
        end do
        end do

        do k=0,nk+1,nk+1
        do j=0,nj+1,nj+1
            kg = minabscc(1,nk,k)
            jg = minabscc(1,nj,j)
            do i=1,ni
                var(m)%r3d(i,j,k) = -var(m)%r3d(i,jg,kg) + &
                                     var(m)%r3d(i,j ,kg) + &
                                     var(m)%r3d(i,jg,k )
            end do
        end do
        end do

        do k=0,nk+1,nk+1
        do j=0,nj+1,nj+1
        do i=0,ni+1,ni+1
            kg = minabscc(1,nk,k)
            jg = minabscc(1,nj,j)
            ig = minabscc(1,ni,i)
            var(m)%r3d(i,j,k) = third*( var(m)%r3d(ig,j,k) + &
                                        var(m)%r3d(i,jg,k) + &
                                        var(m)%r3d(i,j,kg) )
        end do
        end do
        end do
    end do
!$OMP end parallel do

    contains

    function minabscc(a,b,c) result(d)
        use mod_kndconsts, only : kind_int
        implicit none
        integer(kind_int) :: a,b,c
        integer(kind_int) :: d
        integer(kind_int) :: d1,d2

        d1 = abs(a-c)
        d2 = abs(b-c)
        if (d1 == d2) then
            call error_check(1, &
                             "The distance of C-A and C-B is equal in function minabs")
        else if (d1 > d2) then
            d = b
        else
            d = a
        end if

    end function minabscc

end subroutine fill_corner_var_nb_cc

subroutine fill_corner_var_nb_sp(nb,mb_var,nst,ned)
    use mod_kndconsts, only : kind_int
    use mod_constants, only : third
    use mod_datatypes, only : var_block_t,fld_array_t
    use mod_fieldvars, only : mb_topsp
    use mod_openmp
    implicit none
    integer(kind_int)         , intent(in) :: nb
    type(var_block_t), pointer, intent(in) :: mb_var(:)
    integer(kind_int)         , intent(in) :: nst,ned
    integer(kind_int)          :: i,j,k,m,ierr
    integer(kind_int)          :: ig,jg,kg,ni,nj,nk
    type(fld_array_t), pointer :: var(:)

    var => mb_var(nb)%fld

    ni = mb_topsp(nb)%nijk(1)
    nj = mb_topsp(nb)%nijk(2)
    nk = mb_topsp(nb)%nijk(3)
!$OMP parallel do private(ig,jg,kg)
    do m=nst,ned
        do j=0,nj+1,nj+1
        do i=0,ni+1,ni+1
            jg = minabssp(1,nj,j)
            ig = minabssp(1,ni,i)
            do k=1,nk
                var(m)%r3d(i,j,k) = -var(m)%r3d(ig,jg,k) + &
                                     var(m)%r3d(i ,jg,k) + &
                                     var(m)%r3d(ig,j ,k)
            end do
        end do
        end do

        do k=0,nk+1,nk+1
        do i=0,ni+1,ni+1
            kg = minabssp(1,nk,k)
            ig = minabssp(1,ni,i)
            do j=1,nj
                var(m)%r3d(i,j,k) = -var(m)%r3d(ig,j,kg) + &
                                     var(m)%r3d(i ,j,kg) + &
                                     var(m)%r3d(ig,j,k )
            end do
        end do
        end do

        do k=0,nk+1,nk+1
        do j=0,nj+1,nj+1
            kg = minabssp(1,nk,k)
            jg = minabssp(1,nj,j)
            do i=1,ni
                var(m)%r3d(i,j,k) = -var(m)%r3d(i,jg,kg) + &
                                     var(m)%r3d(i,j ,kg) + &
                                     var(m)%r3d(i,jg,k )
            end do
        end do
        end do

        do k=0,nk+1,nk+1
        do j=0,nj+1,nj+1
        do i=0,ni+1,ni+1
            kg = minabssp(1,nk,k)
            jg = minabssp(1,nj,j)
            ig = minabssp(1,ni,i)
            var(m)%r3d(i,j,k) = third*( var(m)%r3d(ig,j,k) + &
                                        var(m)%r3d(i,jg,k) + &
                                        var(m)%r3d(i,j,kg) )
        end do
        end do
        end do
    end do
!$OMP end parallel do

    contains

    function minabssp(a,b,c) result(d)
        use mod_kndconsts, only : kind_int
        implicit none
        integer(kind_int) :: a,b,c
        integer(kind_int) :: d
        integer(kind_int) :: d1,d2

        d1 = abs(a-c)
        d2 = abs(b-c)
        if (d1 == d2) then
            call error_check(1, &
                             "The distance of C-A and C-B is equal in function minabssp")
        else if (d1 > d2) then
            d = b
        else
            d = a
        end if

    end function minabssp

end subroutine fill_corner_var_nb_sp

subroutine ghost_bc_var_all(mb_var,nst,ned)
    use mod_kndconsts, only : kind_int
    use mod_datatypes, only : var_block_t
    use mod_fieldvars, only : nblkcoms,blkcoms
    use mod_interface, only : ghost_bc_var_nb
    implicit none
    type(var_block_t), pointer, intent(in) :: mb_var(:)
    integer(kind_int)         , intent(in) :: nst,ned
    integer(kind_int) :: nc,nb

    do nc=1,nblkcoms
        nb = blkcoms(nc)%nb

        call ghost_bc_var_nb(nb,mb_var,nst,ned)
    end do

end subroutine ghost_bc_var_all

subroutine ghost_bc_var_all_cc(mb_var,nst,ned)
    use mod_kndconsts, only : kind_int
    use mod_datatypes, only : var_block_t
    use mod_fieldvars, only : nblkcoms,blkcomscc
    use mod_interface, only : ghost_bc_var_nb_cc
    implicit none
    type(var_block_t), pointer, intent(in) :: mb_var(:)
    integer(kind_int)         , intent(in) :: nst,ned
    integer(kind_int) :: nc,nb

    do nc=1,nblkcoms
        nb = blkcomscc(nc)%nb

        call ghost_bc_var_nb_cc(nb,mb_var,nst,ned)
    end do

end subroutine ghost_bc_var_all_cc

subroutine ghost_bc_var_all_sp(mb_var,nst,ned)
    use mod_kndconsts, only : kind_int
    use mod_datatypes, only : var_block_t
    use mod_fieldvars, only : nblkcoms,blkcomssp
    use mod_interface, only : ghost_bc_var_nb_sp
    implicit none
    type(var_block_t), pointer, intent(in) :: mb_var(:)
    integer(kind_int)         , intent(in) :: nst,ned
    integer(kind_int) :: nc,nb

    do nc=1,nblkcoms
        nb = blkcomssp(nc)%nb

        call ghost_bc_var_nb_sp(nb,mb_var,nst,ned)
    end do

end subroutine ghost_bc_var_all_sp

subroutine assign_com_var_nb(nb,mb_var,nst,ned,var)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nsw_dir_close
    use mod_datatypes, only : var_block_t,fld_array_t
    use mod_variables, only : nsw_kdir
    use mod_fieldvars, only : mb_top
    implicit none
    integer(kind_int)         , intent(in) :: nb
    type(var_block_t), pointer, intent(in) :: mb_var(:)
    integer(kind_int)         , intent(in) :: nst,ned
    real(kind_real)           , intent(in) :: var(nst:ned)
    integer(kind_int)          :: i,j,k,m,ierr
    integer(kind_int)          :: st(3),ed(3)
    type(fld_array_t), pointer :: fld(:)

    if (nsw_kdir == nsw_dir_close) then
        fld => mb_var(nb)%fld

        st(:) = mb_top(nb)%ndst(:)
        ed(:) = mb_top(nb)%nded(:)

        do m=nst,ned
            do k=st(3),ed(3)
            do j=st(2),ed(2)
            do i=st(1),ed(1)
                fld(m)%r3d(i,j,k) = var(m)
            end do
            end do
            end do
        end do
    end if

end subroutine assign_com_var_nb

subroutine assign_var_via_com_nb(nb,mb_var,nst,ned)  !! 该函数只能作用于dq增量
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nsw_dir_close
    use mod_constants, only : zero,nupd_field_no
    use mod_datatypes, only : var_block_t,fld_array_t
    use mod_variables, only : nsw_kdir
    use mod_fieldvars, only : mb_top
    use mod_openmp
    implicit none
    integer(kind_int)         , intent(in) :: nb
    type(var_block_t), pointer, intent(in) :: mb_var(:)
    integer(kind_int)         , intent(in) :: nst,ned
    integer(kind_int) :: i,j,k,m,nkst,ierr
    integer(kind_int) :: st(3),ed(3),nupd
    type(fld_array_t), pointer :: fld(:)


    if (nsw_kdir == nsw_dir_close) then
        fld => mb_var(nb)%fld

        st(:) = 1
        ed(:) = mb_top(nb)%nijk(:)

        nkst = mb_top(nb)%ndst(3)
!OMP parallel do
        do m=nst,ned
            do k=st(3),ed(3)
            do j=st(2),ed(2)
            do i=st(1),ed(1)
                fld(m)%r3d(i,j,k) = fld(m)%r3d(i,j,nkst)
            end do
            end do
            end do
        end do
!OMP end parallel do
    end if

    nupd = mb_top(nb)%nupd
    if (nupd == nupd_field_no) then
        st(:) = 1
        ed(:) = mb_top(nb)%nijk(:)

        nkst = mb_top(nb)%ndst(3)
!OMP parallel do
        do m=nst,ned
            do k=st(3),ed(3)
            do j=st(2),ed(2)
            do i=st(1),ed(1)
                fld(m)%r3d(i,j,k) = zero
            end do
            end do
            end do
        end do
!OMP end parallel do
    end if

end subroutine assign_var_via_com_nb

subroutine assign_var_via_com_nb_sp(nb,mb_var,nst,ned)  !! 该函数只能作用于dq增量
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nsw_dir_close
    use mod_constants, only : zero,nupd_field_no
    use mod_datatypes, only : var_block_t,fld_array_t
    use mod_variables, only : nsw_kdir
    use mod_fieldvars, only : mb_topsp
    use mod_openmp
    implicit none
    integer(kind_int)         , intent(in) :: nb
    type(var_block_t), pointer, intent(in) :: mb_var(:)
    integer(kind_int)         , intent(in) :: nst,ned
    integer(kind_int) :: i,j,k,m,nkst,ierr
    integer(kind_int) :: st(3),ed(3),nupd
    type(fld_array_t), pointer :: fld(:)


    if (nsw_kdir == nsw_dir_close) then
        fld => mb_var(nb)%fld

        st(:) = 1
        ed(:) = mb_topsp(nb)%nijk(:)

        nkst = mb_topsp(nb)%ndst(3)
!OMP parallel do
        do m=nst,ned
            do k=st(3),ed(3)
            do j=st(2),ed(2)
            do i=st(1),ed(1)
                fld(m)%r3d(i,j,k) = fld(m)%r3d(i,j,nkst)
            end do
            end do
            end do
        end do
!OMP end parallel do
    end if

    nupd = mb_topsp(nb)%nupd
    if (nupd == nupd_field_no) then
        st(:) = 1
        ed(:) = mb_topsp(nb)%nijk(:)

        nkst = mb_topsp(nb)%ndst(3)
!OMP parallel do
        do m=nst,ned
            do k=st(3),ed(3)
            do j=st(2),ed(2)
            do i=st(1),ed(1)
                fld(m)%r3d(i,j,k) = zero
            end do
            end do
            end do
        end do
!OMP end parallel do
    end if

end subroutine assign_var_via_com_nb_sp

subroutine calc_mb_dn_via_node3(mb_vn,nvn, &
                                ngn,sub_ve,nge,sub_dn, &
                                nfsf_ve,nfsf_dn, &
                                mb_dn,ndn,ndir)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_datatypes, only : var_block_t,fld_array_t
    use mod_fieldvars, only : nblkcoms,blkcoms,mb_fsf
    use mod_openmp
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
    integer(kind_int) :: nc,nb,i,j,k,nfs,nfe,ierr
    integer(kind_int) :: ni,nj,nk,st,ed
    real(kind_real), pointer :: vn3d(:,:,:)
    real(kind_real), pointer :: dn3d(:,:,:)
    real(kind_real), pointer :: vn(:,:),ve(:,:),dn(:,:)
    type(fld_array_t), pointer :: fsfs(:),fsfe(:)

    do nc=1,nblkcoms
        nb = blkcoms(nc)%nb

        ni = blkcoms(nc)%top%nijk(1)
        nj = blkcoms(nc)%top%nijk(2)
        nk = blkcoms(nc)%top%nijk(3)

        vn3d => mb_vn(nb)%fld(nvn)%r3d
        dn3d => mb_dn(nb)%fld(ndn)%r3d

        if (ndir == 1) then
            fsfs => mb_fsf(nb,1)%fld
            fsfe => mb_fsf(nb,2)%fld

            st = 1 - ngn
            ed = ni + ngn
!$OMP parallel private(vn,ve,dn,k,j,nfs,nfe)
            allocate(vn(st:ed,1:1), stat=ierr)
            allocate(ve(st:ed,1:1), stat=ierr)
            allocate(dn(st:ed,1:1), stat=ierr)
!$OMP do
            do k=1,nk
            do j=1,nj
                do i=st,ed
                    vn(i,1) = vn3d(i,j,k)
                end do

                nfs = fsfs(nfsf_ve)%i3d(1 ,j,k)
                nfe = fsfe(nfsf_ve)%i3d(ni,j,k)
                call sub_ve(ni,1,1,ngn,vn,nfs,nfe,nge,ve)

                nfs = fsfs(nfsf_dn)%i3d(1 ,j,k)
                nfe = fsfe(nfsf_dn)%i3d(ni,j,k)
                call sub_dn(ni,1,1,ngn,vn,nge,ve,nfs,nfe,dn)

                do i=1,ni
                    dn3d(i,j,k) = dn(i,1)
                end do
            end do
            end do
!$OMP end do
            deallocate(dn, stat=ierr)
            deallocate(ve, stat=ierr)
            deallocate(vn, stat=ierr)
!$OMP end parallel

        else if (ndir == 2) then
            fsfs => mb_fsf(nb,3)%fld
            fsfe => mb_fsf(nb,4)%fld

            st = 1 - ngn
            ed = nj + ngn
!$OMP parallel private(vn,ve,dn,k,j,nfs,nfe)
            allocate(vn(st:ed,1:1), stat=ierr)
            allocate(ve(st:ed,1:1), stat=ierr)
            allocate(dn(st:ed,1:1), stat=ierr)
!$OMP do
            do k=1,nk
            do i=1,ni
                do j=st,ed
                    vn(j,1) = vn3d(i,j,k)
                end do

                nfs = fsfs(nfsf_ve)%i3d(i,1 ,k)
                nfe = fsfe(nfsf_ve)%i3d(i,nj,k)
                call sub_ve(nj,1,1,ngn,vn,nfs,nfe,nge,ve)

                nfs = fsfs(nfsf_dn)%i3d(i,1 ,k)
                nfe = fsfe(nfsf_dn)%i3d(i,nj,k)
                call sub_dn(nj,1,1,ngn,vn,nge,ve,nfs,nfe,dn)

                do j=1,nj
                    dn3d(i,j,k) = dn(j,1)
                end do
            end do
            end do
!$OMP end do
            deallocate(dn, stat=ierr)
            deallocate(ve, stat=ierr)
            deallocate(vn, stat=ierr)
!$OMP end parallel

        else if (ndir == 3) then
            fsfs => mb_fsf(nb,5)%fld
            fsfe => mb_fsf(nb,6)%fld

            st = 1 - ngn
            ed = nk + ngn
!$OMP parallel private(vn,ve,dn,k,j,nfs,nfe)
            allocate(vn(st:ed,1:1), stat=ierr)
            allocate(ve(st:ed,1:1), stat=ierr)
            allocate(dn(st:ed,1:1), stat=ierr)
!$OMP do
            do j=1,nj
            do i=1,ni
                do k=st,ed
                    vn(k,1) = vn3d(i,j,k)
                end do

                nfs = fsfs(nfsf_ve)%i3d(i,j,1 )
                nfe = fsfe(nfsf_ve)%i3d(i,j,nk)
                call sub_ve(nk,1,1,ngn,vn,nfs,nfe,nge,ve)

                nfs = fsfs(nfsf_dn)%i3d(i,j,1 )
                nfe = fsfe(nfsf_dn)%i3d(i,j,nk)
                call sub_dn(nk,1,1,ngn,vn,nge,ve,nfs,nfe,dn)

                do k=1,nk
                    dn3d(i,j,k) = dn(k,1)
                end do
            end do
            end do
!$OMP end do
            deallocate(dn, stat=ierr)
            deallocate(ve, stat=ierr)
            deallocate(vn, stat=ierr)
!$OMP end parallel
        end if
    end do

end subroutine calc_mb_dn_via_node3

subroutine calc_cc_mb_dn_via_node3(mb_vn,nvn, &
                                   ngn,sub_ve,nge,sub_dn, &
                                   nfsf_ve,nfsf_dn, &
                                   mb_dn,ndn,ndir)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_datatypes, only : var_block_t,fld_array_t
    use mod_fieldvars, only : nblkcoms,blkcomscc,mb_fsfcc
    use mod_openmp
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
    integer(kind_int) :: nc,nb,i,j,k,nfs,nfe,ierr
    integer(kind_int) :: ni,nj,nk,st,ed
    real(kind_real), pointer :: vn3d(:,:,:)
    real(kind_real), pointer :: dn3d(:,:,:)
    real(kind_real), pointer :: vn(:,:),ve(:,:),dn(:,:)
    type(fld_array_t), pointer :: fsfs(:),fsfe(:)

    do nc=1,nblkcoms
        nb = blkcomscc(nc)%nb

        ni = blkcomscc(nc)%top%nijk(1)
        nj = blkcomscc(nc)%top%nijk(2)
        nk = blkcomscc(nc)%top%nijk(3)

        vn3d => mb_vn(nb)%fld(nvn)%r3d
        dn3d => mb_dn(nb)%fld(ndn)%r3d

        if (ndir == 1) then
            fsfs => mb_fsfcc(nb,1)%fld
            fsfe => mb_fsfcc(nb,2)%fld

            st = 1 - ngn
            ed = ni + ngn
!$OMP parallel private(vn,ve,dn,k,j,nfs,nfe)
            allocate(vn(st:ed,1:1), stat=ierr)
            allocate(ve(st:ed,1:1), stat=ierr)
            allocate(dn(st:ed,1:1), stat=ierr)
!$OMP do
            do k=1,nk
            do j=1,nj
                do i=st,ed
                    vn(i,1) = vn3d(i,j,k)
                end do

                nfs = fsfs(nfsf_ve)%i3d(1 ,j,k)
                nfe = fsfe(nfsf_ve)%i3d(ni,j,k)
                call sub_ve(ni,1,1,ngn,vn,nfs,nfe,nge,ve)

                nfs = fsfs(nfsf_dn)%i3d(1 ,j,k)
                nfe = fsfe(nfsf_dn)%i3d(ni,j,k)
                call sub_dn(ni,1,1,ngn,vn,nge,ve,nfs,nfe,dn)

                do i=1,ni
                    dn3d(i,j,k) = dn(i,1)
                end do
            end do
            end do
!$OMP end do
            deallocate(dn, stat=ierr)
            deallocate(ve, stat=ierr)
            deallocate(vn, stat=ierr)
!$OMP end parallel

        else if (ndir == 2) then
            fsfs => mb_fsfcc(nb,3)%fld
            fsfe => mb_fsfcc(nb,4)%fld

            st = 1 - ngn
            ed = nj + ngn
!$OMP parallel private(vn,ve,dn,k,j,nfs,nfe)
            allocate(vn(st:ed,1:1), stat=ierr)
            allocate(ve(st:ed,1:1), stat=ierr)
            allocate(dn(st:ed,1:1), stat=ierr)
!$OMP do
            do k=1,nk
            do i=1,ni
                do j=st,ed
                    vn(j,1) = vn3d(i,j,k)
                end do

                nfs = fsfs(nfsf_ve)%i3d(i,1 ,k)
                nfe = fsfe(nfsf_ve)%i3d(i,nj,k)
                call sub_ve(nj,1,1,ngn,vn,nfs,nfe,nge,ve)

                nfs = fsfs(nfsf_dn)%i3d(i,1 ,k)
                nfe = fsfe(nfsf_dn)%i3d(i,nj,k)
                call sub_dn(nj,1,1,ngn,vn,nge,ve,nfs,nfe,dn)

                do j=1,nj
                    dn3d(i,j,k) = dn(j,1)
                end do
            end do
            end do
!$OMP end do
            deallocate(dn, stat=ierr)
            deallocate(ve, stat=ierr)
            deallocate(vn, stat=ierr)
!$OMP end parallel

        else if (ndir == 3) then
            fsfs => mb_fsfcc(nb,5)%fld
            fsfe => mb_fsfcc(nb,6)%fld

            st = 1 - ngn
            ed = nk + ngn
!$OMP parallel private(vn,ve,dn,k,j,nfs,nfe)
            allocate(vn(st:ed,1:1), stat=ierr)
            allocate(ve(st:ed,1:1), stat=ierr)
            allocate(dn(st:ed,1:1), stat=ierr)
!$OMP do
            do j=1,nj
            do i=1,ni
                do k=st,ed
                    vn(k,1) = vn3d(i,j,k)
                end do

                nfs = fsfs(nfsf_ve)%i3d(i,j,1 )
                nfe = fsfe(nfsf_ve)%i3d(i,j,nk)
                call sub_ve(nk,1,1,ngn,vn,nfs,nfe,nge,ve)

                nfs = fsfs(nfsf_dn)%i3d(i,j,1 )
                nfe = fsfe(nfsf_dn)%i3d(i,j,nk)
                call sub_dn(nk,1,1,ngn,vn,nge,ve,nfs,nfe,dn)

                do k=1,nk
                    dn3d(i,j,k) = dn(k,1)
                end do
            end do
            end do
!$OMP end do
            deallocate(dn, stat=ierr)
            deallocate(ve, stat=ierr)
            deallocate(vn, stat=ierr)
!$OMP end parallel
        end if
    end do

end subroutine calc_cc_mb_dn_via_node3

subroutine calc_sp_mb_dn_via_node3(mb_vn,nvn, &
                                   ngn,sub_ve,nge,sub_dn, &
                                   nfsf_ve,nfsf_dn, &
                                   mb_dn,ndn,ndir)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_datatypes, only : var_block_t,fld_array_t
    use mod_fieldvars, only : nblkcoms,blkcomssp,mb_fsfsp
    use mod_openmp
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
    integer(kind_int) :: nc,nb,i,j,k,nfs,nfe,ierr
    integer(kind_int) :: ni,nj,nk,st,ed
    real(kind_real), pointer :: vn3d(:,:,:)
    real(kind_real), pointer :: dn3d(:,:,:)
    real(kind_real), pointer :: vn(:,:),ve(:,:),dn(:,:)
    type(fld_array_t), pointer :: fsfs(:),fsfe(:)

    do nc=1,nblkcoms
        nb = blkcomssp(nc)%nb

        ni = blkcomssp(nc)%top%nijk(1)
        nj = blkcomssp(nc)%top%nijk(2)
        nk = blkcomssp(nc)%top%nijk(3)

        vn3d => mb_vn(nb)%fld(nvn)%r3d
        dn3d => mb_dn(nb)%fld(ndn)%r3d

        if (ndir == 1) then
            fsfs => mb_fsfsp(nb,1)%fld
            fsfe => mb_fsfsp(nb,2)%fld

            st = 1 - ngn
            ed = ni + ngn
!$OMP parallel private(vn,ve,dn,k,j,nfs,nfe)
            allocate(vn(st:ed,1:1), stat=ierr)
            allocate(ve(st:ed,1:1), stat=ierr)
            allocate(dn(st:ed,1:1), stat=ierr)
!$OMP do
            do k=1,nk
            do j=1,nj
                do i=st,ed
                    vn(i,1) = vn3d(i,j,k)
                end do

                nfs = fsfs(nfsf_ve)%i3d(1 ,j,k)
                nfe = fsfe(nfsf_ve)%i3d(ni,j,k)
                call sub_ve(ni,1,1,ngn,vn,nfs,nfe,nge,ve)

                nfs = fsfs(nfsf_dn)%i3d(1 ,j,k)
                nfe = fsfe(nfsf_dn)%i3d(ni,j,k)
                call sub_dn(ni,1,1,ngn,vn,nge,ve,nfs,nfe,dn)

                do i=1,ni
                    dn3d(i,j,k) = dn(i,1)
                end do
            end do
            end do
!$OMP end do
            deallocate(dn, stat=ierr)
            deallocate(ve, stat=ierr)
            deallocate(vn, stat=ierr)
!$OMP end parallel

        else if (ndir == 2) then
            fsfs => mb_fsfsp(nb,3)%fld
            fsfe => mb_fsfsp(nb,4)%fld

            st = 1 - ngn
            ed = nj + ngn
!$OMP parallel private(vn,ve,dn,k,j,nfs,nfe)
            allocate(vn(st:ed,1:1), stat=ierr)
            allocate(ve(st:ed,1:1), stat=ierr)
            allocate(dn(st:ed,1:1), stat=ierr)
!$OMP do
            do k=1,nk
            do i=1,ni
                do j=st,ed
                    vn(j,1) = vn3d(i,j,k)
                end do

                nfs = fsfs(nfsf_ve)%i3d(i,1 ,k)
                nfe = fsfe(nfsf_ve)%i3d(i,nj,k)
                call sub_ve(nj,1,1,ngn,vn,nfs,nfe,nge,ve)

                nfs = fsfs(nfsf_dn)%i3d(i,1 ,k)
                nfe = fsfe(nfsf_dn)%i3d(i,nj,k)
                call sub_dn(nj,1,1,ngn,vn,nge,ve,nfs,nfe,dn)

                do j=1,nj
                    dn3d(i,j,k) = dn(j,1)
                end do
            end do
            end do
!$OMP end do
            deallocate(dn, stat=ierr)
            deallocate(ve, stat=ierr)
            deallocate(vn, stat=ierr)
!$OMP end parallel

        else if (ndir == 3) then
            fsfs => mb_fsfsp(nb,5)%fld
            fsfe => mb_fsfsp(nb,6)%fld

            st = 1 - ngn
            ed = nk + ngn
!$OMP parallel private(vn,ve,dn,k,j,nfs,nfe)
            allocate(vn(st:ed,1:1), stat=ierr)
            allocate(ve(st:ed,1:1), stat=ierr)
            allocate(dn(st:ed,1:1), stat=ierr)
!$OMP do
            do j=1,nj
            do i=1,ni
                do k=st,ed
                    vn(k,1) = vn3d(i,j,k)
                end do

                nfs = fsfs(nfsf_ve)%i3d(i,j,1 )
                nfe = fsfe(nfsf_ve)%i3d(i,j,nk)
                call sub_ve(nk,1,1,ngn,vn,nfs,nfe,nge,ve)

                nfs = fsfs(nfsf_dn)%i3d(i,j,1 )
                nfe = fsfe(nfsf_dn)%i3d(i,j,nk)
                call sub_dn(nk,1,1,ngn,vn,nge,ve,nfs,nfe,dn)

                do k=1,nk
                    dn3d(i,j,k) = dn(k,1)
                end do
            end do
            end do
!$OMP end do
            deallocate(dn, stat=ierr)
            deallocate(ve, stat=ierr)
            deallocate(vn, stat=ierr)
!$OMP end parallel
        end if
    end do

end subroutine calc_sp_mb_dn_via_node3

subroutine calc_mb_dn_via_node3_exp(mb_vn,nvn, &
                                     ngn,sub_ve,nge,sub_dn, &
                                     nfsf_ve,nfsf_dn, &
                                     mb_dn,ndn,ndir)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_inter_scheme,nbc_intbc_scheme
    use mod_datatypes, only : var_block_t,fld_array_t
    use mod_variables, only : nghnode
    use mod_fieldvars, only : nblkcoms,blkcoms,mb_fsf
    use mod_openmp
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
    integer(kind_int) :: nc,nb,i,j,k,nfs,nfe,ierr
    integer(kind_int) :: ni,nj,nk,st,ed,st1,ed1,st2,ed2,st0,ed0,iflg,jflg,kflg
    real(kind_real), pointer :: vn3d(:,:,:)
    real(kind_real), pointer :: dn3d(:,:,:)
    real(kind_real), pointer :: vn(:,:),ve(:,:),dn(:,:)
    type(fld_array_t), pointer :: fsfs(:),fsfe(:)

    do nc=1,nblkcoms
        nb = blkcoms(nc)%nb

        ni = blkcoms(nc)%top%nijk(1)
        nj = blkcoms(nc)%top%nijk(2)
        nk = blkcoms(nc)%top%nijk(3)

        vn3d => mb_vn(nb)%fld(nvn)%r3d
        dn3d => mb_dn(nb)%fld(ndn)%r3d

        if (ndir == 1) then
            fsfs => mb_fsf(nb,1)%fld
            fsfe => mb_fsf(nb,2)%fld

            st = 1 - ngn
            ed = ni + ngn
            st1 = 1 - ngn + nghnode
            ed1 = nk + ngn - nghnode
            st2 = 1 - ngn + nghnode
            ed2 = nj + ngn - nghnode
!$OMP parallel private(vn,ve,dn,k,j,nfs,nfe)
            allocate(vn(st:ed,1:1), stat=ierr)
            allocate(ve(st:ed,1:1), stat=ierr)
            allocate(dn(st:ed,1:1), stat=ierr)
!$OMP do
            do k=st1,ed1
            kflg=k
            if(kflg<1) then
                kflg=1
            else if(kflg>nk) then
                kflg=nk
            end if
            do j=st2,ed2
                jflg=j
                if(jflg<1) then
                    jflg=1
                else if(jflg>nj) then
                    jflg=nj
                end if
                do i=st,ed
                    vn(i,1) = vn3d(i,j,k)
                end do

                nfs = fsfs(nfsf_ve)%i3d(1 ,jflg,kflg)
                nfe = fsfe(nfsf_ve)%i3d(ni,jflg,kflg)
                call sub_ve(ni,1,1,ngn,vn,nfs,nfe,nge,ve)

                nfs = fsfs(nfsf_dn)%i3d(1 ,jflg,kflg)
                nfe = fsfe(nfsf_dn)%i3d(ni,jflg,kflg)
                call sub_dn(ni,1,1,ngn,vn,nge,ve,nfs,nfe,dn)

                if (nfs == nbc_inter_scheme .or. nfs == nbc_intbc_scheme) then
                    st0 = st+nghnode
                else
                    st0 = 1
                end if

                if (nfe == nbc_inter_scheme .or. nfe == nbc_intbc_scheme) then
                    ed0 = ed-nghnode
                else
                    ed0 = ni
                end if

                do i=st0,ed0
                    dn3d(i,j,k) = dn(i,1)
                end do
            end do
            end do
!$OMP end do
            deallocate(dn, stat=ierr)
            deallocate(ve, stat=ierr)
            deallocate(vn, stat=ierr)
!$OMP end parallel

        else if (ndir == 2) then
            fsfs => mb_fsf(nb,3)%fld
            fsfe => mb_fsf(nb,4)%fld

            st = 1 - ngn
            ed = nj + ngn
            st1 = 1 - ngn + nghnode
            ed1 = nk + ngn - nghnode
            st2 = 1 - ngn + nghnode
            ed2 = ni + ngn - nghnode
!$OMP parallel private(vn,ve,dn,k,j,nfs,nfe)
            allocate(vn(st:ed,1:1), stat=ierr)
            allocate(ve(st:ed,1:1), stat=ierr)
            allocate(dn(st:ed,1:1), stat=ierr)
!$OMP do
            do k=st1,ed1
            kflg=k
            if(kflg<1) then
                kflg=1
            else if(kflg>nk) then
                kflg=nk
            end if
            do i=st2,ed2
                iflg=i
                if(iflg<1) then
                    iflg=1
                else if(iflg>ni) then
                    iflg=ni
                end if
                do j=st,ed
                    vn(j,1) = vn3d(i,j,k)
                end do

                nfs = fsfs(nfsf_ve)%i3d(iflg,1 ,kflg)
                nfe = fsfe(nfsf_ve)%i3d(iflg,nj,kflg)
                call sub_ve(nj,1,1,ngn,vn,nfs,nfe,nge,ve)

                nfs = fsfs(nfsf_dn)%i3d(iflg,1 ,kflg)
                nfe = fsfe(nfsf_dn)%i3d(iflg,nj,kflg)
                call sub_dn(nj,1,1,ngn,vn,nge,ve,nfs,nfe,dn)

                if (nfs == nbc_inter_scheme .or. nfs == nbc_intbc_scheme) then
                    st0 = st+nghnode
                else
                    st0 = 1
                end if

                if (nfe == nbc_inter_scheme .or. nfe == nbc_intbc_scheme) then
                    ed0 = ed-nghnode
                else
                    ed0 = nj
                end if

                do j=st0,ed0
                    dn3d(i,j,k) = dn(j,1)
                end do
            end do
            end do
!$OMP end do
            deallocate(dn, stat=ierr)
            deallocate(ve, stat=ierr)
            deallocate(vn, stat=ierr)
!$OMP end parallel

        else if (ndir == 3) then
            fsfs => mb_fsf(nb,5)%fld
            fsfe => mb_fsf(nb,6)%fld

            st = 1 - ngn
            ed = nk + ngn
            st1 = 1 - ngn + nghnode
            ed1 = nj + ngn - nghnode
            st2 = 1 - ngn + nghnode
            ed2 = ni + ngn - nghnode
!$OMP parallel private(vn,ve,dn,k,j,nfs,nfe)
            allocate(vn(st:ed,1:1), stat=ierr)
            allocate(ve(st:ed,1:1), stat=ierr)
            allocate(dn(st:ed,1:1), stat=ierr)
!$OMP do
            do j=st1,ed1
            jflg=j
            if(jflg<1) then
                jflg=1
            else if(jflg>nj) then
                jflg=nj
            end if
            do i=st2,ed2
                iflg=i
                if(iflg<1) then
                    iflg=1
                else if(iflg>ni) then
                    iflg=ni
                end if
                do k=st,ed
                    vn(k,1) = vn3d(i,j,k)
                end do

                nfs = fsfs(nfsf_ve)%i3d(iflg,jflg,1 )
                nfe = fsfe(nfsf_ve)%i3d(iflg,jflg,nk)
                call sub_ve(nk,1,1,ngn,vn,nfs,nfe,nge,ve)

                nfs = fsfs(nfsf_dn)%i3d(iflg,jflg,1 )
                nfe = fsfe(nfsf_dn)%i3d(iflg,jflg,nk)
                call sub_dn(nk,1,1,ngn,vn,nge,ve,nfs,nfe,dn)

                if (nfs == nbc_inter_scheme .or. nfs == nbc_intbc_scheme) then
                    st0 = st+nghnode
                else
                    st0 = 1
                end if

                if (nfe == nbc_inter_scheme .or. nfe == nbc_intbc_scheme) then
                    ed0 = ed-nghnode
                else
                    ed0 = nk
                end if

                do k=st0,ed0
                    dn3d(i,j,k) = dn(k,1)
                end do
            end do
            end do
!$OMP end do
            deallocate(dn, stat=ierr)
            deallocate(ve, stat=ierr)
            deallocate(vn, stat=ierr)
!$OMP end parallel
        end if
    end do

end subroutine calc_mb_dn_via_node3_exp

subroutine calc_mb_dn_via_node3_sp(mb_vn,nvn, &
                                   ngn,sub_ve,nge,sub_dn, &
                                   nfsf_ve,nfsf_dn, &
                                   mb_dn,ndn,ndir)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_datatypes, only : var_block_t,fld_array_t
    use mod_fieldvars, only : nblkcoms,blkcomssp,mb_fsfsp
    use mod_openmp
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
    integer(kind_int) :: nc,nb,i,j,k,nfs,nfe,ierr
    integer(kind_int) :: ni,nj,nk,st,ed
    real(kind_real), pointer :: vn3d(:,:,:)
    real(kind_real), pointer :: dn3d(:,:,:)
    real(kind_real), pointer :: vn(:,:),ve(:,:),dn(:,:)
    type(fld_array_t), pointer :: fsfs(:),fsfe(:)

    do nc=1,nblkcoms
        nb = blkcomssp(nc)%nb

        ni = blkcomssp(nc)%top%nijk(1)
        nj = blkcomssp(nc)%top%nijk(2)
        nk = blkcomssp(nc)%top%nijk(3)

        vn3d => mb_vn(nb)%fld(nvn)%r3d
        dn3d => mb_dn(nb)%fld(ndn)%r3d

        if (ndir == 1) then
            fsfs => mb_fsfsp(nb,1)%fld
            fsfe => mb_fsfsp(nb,2)%fld

            st = 1 - ngn
            ed = ni + ngn
!$OMP parallel private(vn,ve,dn,k,j,nfs,nfe)
            allocate(vn(st:ed,1:1), stat=ierr)
            allocate(ve(st:ed,1:1), stat=ierr)
            allocate(dn(st:ed,1:1), stat=ierr)
!$OMP do
            do k=1,nk
            do j=1,nj
                do i=st,ed
                    vn(i,1) = vn3d(i,j,k)
                end do

                nfs = fsfs(nfsf_ve)%i3d(1 ,j,k)
                nfe = fsfe(nfsf_ve)%i3d(ni,j,k)
                call sub_ve(ni,1,1,ngn,vn,nfs,nfe,nge,ve)

                nfs = fsfs(nfsf_dn)%i3d(1 ,j,k)
                nfe = fsfe(nfsf_dn)%i3d(ni,j,k)
                call sub_dn(ni,1,1,ngn,vn,nge,ve,nfs,nfe,dn)

                do i=1,ni
                    dn3d(i,j,k) = dn(i,1)
                end do
            end do
            end do
!$OMP end do
            deallocate(dn, stat=ierr)
            deallocate(ve, stat=ierr)
            deallocate(vn, stat=ierr)
!$OMP end parallel

        else if (ndir == 2) then
            fsfs => mb_fsfsp(nb,3)%fld
            fsfe => mb_fsfsp(nb,4)%fld

            st = 1 - ngn
            ed = nj + ngn
!$OMP parallel private(vn,ve,dn,k,j,nfs,nfe)
            allocate(vn(st:ed,1:1), stat=ierr)
            allocate(ve(st:ed,1:1), stat=ierr)
            allocate(dn(st:ed,1:1), stat=ierr)
!$OMP do
            do k=1,nk
            do i=1,ni
                do j=st,ed
                    vn(j,1) = vn3d(i,j,k)
                end do

                nfs = fsfs(nfsf_ve)%i3d(i,1 ,k)
                nfe = fsfe(nfsf_ve)%i3d(i,nj,k)
                call sub_ve(nj,1,1,ngn,vn,nfs,nfe,nge,ve)

                nfs = fsfs(nfsf_dn)%i3d(i,1 ,k)
                nfe = fsfe(nfsf_dn)%i3d(i,nj,k)
                call sub_dn(nj,1,1,ngn,vn,nge,ve,nfs,nfe,dn)

                do j=1,nj
                    dn3d(i,j,k) = dn(j,1)
                end do
            end do
            end do
!$OMP end do
            deallocate(dn, stat=ierr)
            deallocate(ve, stat=ierr)
            deallocate(vn, stat=ierr)
!$OMP end parallel

        else if (ndir == 3) then
            fsfs => mb_fsfsp(nb,5)%fld
            fsfe => mb_fsfsp(nb,6)%fld

            st = 1 - ngn
            ed = nk + ngn
!$OMP parallel private(vn,ve,dn,k,j,nfs,nfe)
            allocate(vn(st:ed,1:1), stat=ierr)
            allocate(ve(st:ed,1:1), stat=ierr)
            allocate(dn(st:ed,1:1), stat=ierr)
!$OMP do
            do j=1,nj
            do i=1,ni
                do k=st,ed
                    vn(k,1) = vn3d(i,j,k)
                end do

                nfs = fsfs(nfsf_ve)%i3d(i,j,1 )
                nfe = fsfe(nfsf_ve)%i3d(i,j,nk)
                call sub_ve(nk,1,1,ngn,vn,nfs,nfe,nge,ve)

                nfs = fsfs(nfsf_dn)%i3d(i,j,1 )
                nfe = fsfe(nfsf_dn)%i3d(i,j,nk)
                call sub_dn(nk,1,1,ngn,vn,nge,ve,nfs,nfe,dn)

                do k=1,nk
                    dn3d(i,j,k) = dn(k,1)
                end do
            end do
            end do
!$OMP end do
            deallocate(dn, stat=ierr)
            deallocate(ve, stat=ierr)
            deallocate(vn, stat=ierr)
!$OMP end parallel
        end if
    end do

end subroutine calc_mb_dn_via_node3_sp

subroutine calc_mb_duvwt_via_node3_sp(mb_vn,nvn, &
                                      ngn,sub_ve,nge,sub_dn, &
                                      nfsf_ve,nfsf_dn, &
                                      mb_dn,ndn,ndir)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_datatypes, only : var_block_t,fld_array_t
    use mod_fieldvars, only : nblkcoms,blkcomssp,mb_fsfsp,mb_fsffp
    use mod_openmp
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
    integer(kind_int) :: nc,nb,i,j,k,nfs,nfe,ierr
    integer(kind_int) :: ni,nj,nk,st,ed,nr,bctype
    real(kind_real), pointer :: vn3d(:,:,:)
    real(kind_real), pointer :: dn3d(:,:,:)
    real(kind_real), pointer :: vn(:,:),ve(:,:),dn(:,:)
    type(fld_array_t), pointer :: fsfs(:),fsfe(:),fsfps(:),fsfpe(:)

    do nc=1,nblkcoms
        nb = blkcomssp(nc)%nb

        ni = blkcomssp(nc)%top%nijk(1)
        nj = blkcomssp(nc)%top%nijk(2)
        nk = blkcomssp(nc)%top%nijk(3)

        vn3d => mb_vn(nb)%fld(nvn)%r3d
        dn3d => mb_dn(nb)%fld(ndn)%r3d

        if (ndir == 1) then
            fsfs => mb_fsfsp(nb,1)%fld
            fsfe => mb_fsfsp(nb,2)%fld
            fsfps=> mb_fsffp(nb,1)%fld    
            fsfpe=> mb_fsffp(nb,2)%fld             

            st = 1 - ngn
            ed = ni + ngn
!$OMP parallel private(vn,ve,dn,k,j,nfs,nfe)
            allocate(vn(st:ed,1:1), stat=ierr)
            allocate(ve(st:ed,1:1), stat=ierr)
            allocate(dn(st:ed,1:1), stat=ierr)
!$OMP do
            do k=1,nk
            do j=1,nj
                do i=st,ed
                    vn(i,1) = vn3d(i,j,k)
                end do

                nfs = fsfs(nfsf_ve)%i3d(1 ,j,k)
                nfe = fsfe(nfsf_ve)%i3d(ni,j,k)
                call sub_ve(ni,1,1,ngn,vn,nfs,nfe,nge,ve)

                do i=1,ed,ni
                    if(i==1 .and. nvn/=4) then              
                        nr     = fsfps(1)%i3d(1,j,k)
                        bctype = fsfps(2)%i3d(1,j,k)
                        if(bctype == 2) then
                            ve(1,1) = 0.0
                        end if
                    end if
                    if(i==ni+1 .and. nvn/=4) then
                        nr     = fsfpe(1)%i3d(ni+1,j,k)
                        bctype = fsfpe(2)%i3d(ni+1,j,k)
                        if(bctype == 2) then
                            ve(ni+1,1) = 0.0
                        end if
                    end if                    
                end do
                nfs = fsfs(nfsf_dn)%i3d(1 ,j,k)
                nfe = fsfe(nfsf_dn)%i3d(ni,j,k)
                call sub_dn(ni,1,1,ngn,vn,nge,ve,nfs,nfe,dn)

                do i=1,ni
                    dn3d(i,j,k) = dn(i,1)
                end do
            end do
            end do
!$OMP end do
            deallocate(dn, stat=ierr)
            deallocate(ve, stat=ierr)
            deallocate(vn, stat=ierr)
!$OMP end parallel

        else if (ndir == 2) then
            fsfs => mb_fsfsp(nb,3)%fld
            fsfe => mb_fsfsp(nb,4)%fld
            fsfps=> mb_fsffp(nb,3)%fld
            fsfpe=> mb_fsffp(nb,4)%fld            

            st = 1 - ngn
            ed = nj + ngn
!$OMP parallel private(vn,ve,dn,k,j,nfs,nfe)
            allocate(vn(st:ed,1:1), stat=ierr)
            allocate(ve(st:ed,1:1), stat=ierr)
            allocate(dn(st:ed,1:1), stat=ierr)
!$OMP do
            do k=1,nk
            do i=1,ni
                do j=st,ed
                    vn(j,1) = vn3d(i,j,k)
                end do

                nfs = fsfs(nfsf_ve)%i3d(i,1 ,k)
                nfe = fsfe(nfsf_ve)%i3d(i,nj,k)
                call sub_ve(nj,1,1,ngn,vn,nfs,nfe,nge,ve)

                do j=1,ed,nj
                    if(j==1 .and. nvn/=4) then              
                        nr     = fsfps(1)%i3d(i,1,k)
                        bctype = fsfps(2)%i3d(i,1,k)
                        if(bctype == 2) then
                            ve(1,1) = 0.0
                        end if
                    end if
                    if(j==nj+1 .and. nvn/=4) then
                        nr     = fsfpe(1)%i3d(i,nj+1,k)
                        bctype = fsfpe(2)%i3d(i,nj+1,k)
                        if(bctype == 2) then
                            ve(nj+1,1) = 0.0
                        end if
                    end if                    
                end do                
                nfs = fsfs(nfsf_dn)%i3d(i,1 ,k)
                nfe = fsfe(nfsf_dn)%i3d(i,nj,k)
                call sub_dn(nj,1,1,ngn,vn,nge,ve,nfs,nfe,dn)

                do j=1,nj
                    dn3d(i,j,k) = dn(j,1)
                end do
            end do
            end do
!$OMP end do
            deallocate(dn, stat=ierr)
            deallocate(ve, stat=ierr)
            deallocate(vn, stat=ierr)
!$OMP end parallel

        else if (ndir == 3) then
            fsfs => mb_fsfsp(nb,5)%fld
            fsfe => mb_fsfsp(nb,6)%fld
            fsfps=> mb_fsffp(nb,5)%fld
            fsfpe=> mb_fsffp(nb,6)%fld             

            st = 1 - ngn
            ed = nk + ngn
!$OMP parallel private(vn,ve,dn,k,j,nfs,nfe)
            allocate(vn(st:ed,1:1), stat=ierr)
            allocate(ve(st:ed,1:1), stat=ierr)
            allocate(dn(st:ed,1:1), stat=ierr)
!$OMP do
            do j=1,nj
            do i=1,ni
                do k=st,ed
                    vn(k,1) = vn3d(i,j,k)
                end do

                nfs = fsfs(nfsf_ve)%i3d(i,j,1 )
                nfe = fsfe(nfsf_ve)%i3d(i,j,nk)
                call sub_ve(nk,1,1,ngn,vn,nfs,nfe,nge,ve)

                do k=1,ed,nk
                    if(k==1 .and. nvn/=4) then              
                        nr     = fsfps(1)%i3d(i,j,1)
                        bctype = fsfps(2)%i3d(i,j,1)
                        if(bctype == 2) then
                            ve(1,1) = 0.0
                        end if
                    end if
                    if(k==nk+1 .and. nvn/=4) then
                        nr     = fsfpe(1)%i3d(i,j,nk+1)
                        bctype = fsfpe(2)%i3d(i,j,nk+1)
                        if(bctype == 2) then
                            ve(nk+1,1) = 0.0
                        end if
                    end if                    
                end do                
                nfs = fsfs(nfsf_dn)%i3d(i,j,1 )
                nfe = fsfe(nfsf_dn)%i3d(i,j,nk)
                call sub_dn(nk,1,1,ngn,vn,nge,ve,nfs,nfe,dn)

                do k=1,nk
                    dn3d(i,j,k) = dn(k,1)
                end do
            end do
            end do
!$OMP end do
            deallocate(dn, stat=ierr)
            deallocate(ve, stat=ierr)
            deallocate(vn, stat=ierr)
!$OMP end parallel
        end if
    end do

end subroutine calc_mb_duvwt_via_node3_sp

subroutine calc_mb_dn_via_node2(mb_vn1,nvn1, &
                                mb_vn2,nvn2, &
                                ngn,sub_ve,nge,sub_dn, &
                                nfsf_ve,nfsf_dn, &
                                mb_dn,ndn,ndir)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_datatypes, only : var_block_t,fld_array_t
    use mod_fieldvars, only : nblkcoms,blkcoms,mb_fsf
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
    integer(kind_int) :: nc,nb,i,j,k,nfs,nfe,ierr
    integer(kind_int) :: ni,nj,nk,st,ed
    real(kind_real), pointer :: vn3d1(:,:,:),vn3d2(:,:,:)
    real(kind_real), pointer :: dn3d(:,:,:)
    real(kind_real), pointer :: vn(:,:),ve(:,:),dn(:,:)
    type(fld_array_t), pointer :: fsfs(:),fsfe(:)

    do nc=1,nblkcoms
        nb = blkcoms(nc)%nb

        ni = blkcoms(nc)%top%nijk(1)
        nj = blkcoms(nc)%top%nijk(2)
        nk = blkcoms(nc)%top%nijk(3)

        vn3d1 => mb_vn1(nb)%fld(nvn1)%r3d
        vn3d2 => mb_vn2(nb)%fld(nvn2)%r3d
        dn3d  => mb_dn(nb)%fld(ndn)%r3d

        if (ndir == 1) then
            fsfs => mb_fsf(nb,1)%fld
            fsfe => mb_fsf(nb,2)%fld

            st = 1 - ngn
            ed = ni + ngn
            allocate(vn(st:ed,1:1), stat=ierr)
            allocate(ve(st:ed,1:1), stat=ierr)
            allocate(dn(st:ed,1:1), stat=ierr)
            do k=1,nk
            do j=1,nj
                do i=st,ed
                    vn(i,1) = vn3d1(i,j,k)*vn3d2(i,j,k)
                end do

                nfs = fsfs(nfsf_ve)%i3d(1 ,j,k)
                nfe = fsfe(nfsf_ve)%i3d(ni,j,k)
                call sub_ve(ni,1,1,ngn,vn,nfs,nfe,nge,ve)

                nfs = fsfs(nfsf_dn)%i3d(1 ,j,k)
                nfe = fsfe(nfsf_dn)%i3d(ni,j,k)
                call sub_dn(ni,1,1,ngn,vn,nge,ve,nfs,nfe,dn)

                do i=1,ni
                    dn3d(i,j,k) = dn(i,1)
                end do
            end do
            end do
            deallocate(dn, stat=ierr)
            deallocate(ve, stat=ierr)
            deallocate(vn, stat=ierr)

        else if (ndir == 2) then
            fsfs => mb_fsf(nb,3)%fld
            fsfe => mb_fsf(nb,4)%fld

            st = 1 - ngn
            ed = nj + ngn
            allocate(vn(st:ed,1:1), stat=ierr)
            allocate(ve(st:ed,1:1), stat=ierr)
            allocate(dn(st:ed,1:1), stat=ierr)
            do k=1,nk
            do i=1,ni
                do j=st,ed
                    vn(j,1) = vn3d1(i,j,k)*vn3d2(i,j,k)
                end do

                nfs = fsfs(nfsf_ve)%i3d(i,1 ,k)
                nfe = fsfe(nfsf_ve)%i3d(i,nj,k)
                call sub_ve(nj,1,1,ngn,vn,nfs,nfe,nge,ve)

                nfs = fsfs(nfsf_dn)%i3d(i,1 ,k)
                nfe = fsfe(nfsf_dn)%i3d(i,nj,k)
                call sub_dn(nj,1,1,ngn,vn,nge,ve,nfs,nfe,dn)

                do j=1,nj
                    dn3d(i,j,k) = dn(j,1)
                end do
            end do
            end do
            deallocate(dn, stat=ierr)
            deallocate(ve, stat=ierr)
            deallocate(vn, stat=ierr)

        else if (ndir == 3) then
            fsfs => mb_fsf(nb,5)%fld
            fsfe => mb_fsf(nb,6)%fld

            st = 1 - ngn
            ed = nk + ngn
            allocate(vn(st:ed,1:1), stat=ierr)
            allocate(ve(st:ed,1:1), stat=ierr)
            allocate(dn(st:ed,1:1), stat=ierr)
            do j=1,nj
            do i=1,ni
                do k=st,ed
                    vn(k,1) = vn3d1(i,j,k)*vn3d2(i,j,k)
                end do

                nfs = fsfs(nfsf_ve)%i3d(i,j,1 )
                nfe = fsfe(nfsf_ve)%i3d(i,j,nk)
                call sub_ve(nk,1,1,ngn,vn,nfs,nfe,nge,ve)

                nfs = fsfs(nfsf_dn)%i3d(i,j,1 )
                nfe = fsfe(nfsf_dn)%i3d(i,j,nk)
                call sub_dn(nk,1,1,ngn,vn,nge,ve,nfs,nfe,dn)

                do k=1,nk
                    dn3d(i,j,k) = dn(k,1)
                end do
            end do
            end do
            deallocate(dn, stat=ierr)
            deallocate(ve, stat=ierr)
            deallocate(vn, stat=ierr)
        end if
    end do

end subroutine calc_mb_dn_via_node2

subroutine calc_cc_mb_dn_via_node2(mb_vn1,nvn1, &
                                   mb_vn2,nvn2, &
                                   ngn,sub_ve,nge,sub_dn, &
                                   nfsf_ve,nfsf_dn, &
                                   mb_dn,ndn,ndir)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_datatypes, only : var_block_t,fld_array_t
    use mod_fieldvars, only : nblkcoms,blkcomscc,mb_fsfcc
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
    integer(kind_int) :: nc,nb,i,j,k,nfs,nfe,ierr
    integer(kind_int) :: ni,nj,nk,st,ed
    real(kind_real), pointer :: vn3d1(:,:,:),vn3d2(:,:,:)
    real(kind_real), pointer :: dn3d(:,:,:)
    real(kind_real), pointer :: vn(:,:),ve(:,:),dn(:,:)
    type(fld_array_t), pointer :: fsfs(:),fsfe(:)

    do nc=1,nblkcoms
        nb = blkcomscc(nc)%nb

        ni = blkcomscc(nc)%top%nijk(1)
        nj = blkcomscc(nc)%top%nijk(2)
        nk = blkcomscc(nc)%top%nijk(3)

        vn3d1 => mb_vn1(nb)%fld(nvn1)%r3d
        vn3d2 => mb_vn2(nb)%fld(nvn2)%r3d
        dn3d  => mb_dn(nb)%fld(ndn)%r3d

        if (ndir == 1) then
            fsfs => mb_fsfcc(nb,1)%fld
            fsfe => mb_fsfcc(nb,2)%fld

            st = 1 - ngn
            ed = ni + ngn
            allocate(vn(st:ed,1:1), stat=ierr)
            allocate(ve(st:ed,1:1), stat=ierr)
            allocate(dn(st:ed,1:1), stat=ierr)
            do k=1,nk
            do j=1,nj
                do i=st,ed
                    vn(i,1) = vn3d1(i,j,k)*vn3d2(i,j,k)
                end do

                nfs = fsfs(nfsf_ve)%i3d(1 ,j,k)
                nfe = fsfe(nfsf_ve)%i3d(ni,j,k)
                call sub_ve(ni,1,1,ngn,vn,nfs,nfe,nge,ve)

                nfs = fsfs(nfsf_dn)%i3d(1 ,j,k)
                nfe = fsfe(nfsf_dn)%i3d(ni,j,k)
                call sub_dn(ni,1,1,ngn,vn,nge,ve,nfs,nfe,dn)

                do i=1,ni
                    dn3d(i,j,k) = dn(i,1)
                end do
            end do
            end do
            deallocate(dn, stat=ierr)
            deallocate(ve, stat=ierr)
            deallocate(vn, stat=ierr)

        else if (ndir == 2) then
            fsfs => mb_fsfcc(nb,3)%fld
            fsfe => mb_fsfcc(nb,4)%fld

            st = 1 - ngn
            ed = nj + ngn
            allocate(vn(st:ed,1:1), stat=ierr)
            allocate(ve(st:ed,1:1), stat=ierr)
            allocate(dn(st:ed,1:1), stat=ierr)
            do k=1,nk
            do i=1,ni
                do j=st,ed
                    vn(j,1) = vn3d1(i,j,k)*vn3d2(i,j,k)
                end do

                nfs = fsfs(nfsf_ve)%i3d(i,1 ,k)
                nfe = fsfe(nfsf_ve)%i3d(i,nj,k)
                call sub_ve(nj,1,1,ngn,vn,nfs,nfe,nge,ve)

                nfs = fsfs(nfsf_dn)%i3d(i,1 ,k)
                nfe = fsfe(nfsf_dn)%i3d(i,nj,k)
                call sub_dn(nj,1,1,ngn,vn,nge,ve,nfs,nfe,dn)

                do j=1,nj
                    dn3d(i,j,k) = dn(j,1)
                end do
            end do
            end do
            deallocate(dn, stat=ierr)
            deallocate(ve, stat=ierr)
            deallocate(vn, stat=ierr)

        else if (ndir == 3) then
            fsfs => mb_fsfcc(nb,5)%fld
            fsfe => mb_fsfcc(nb,6)%fld

            st = 1 - ngn
            ed = nk + ngn
            allocate(vn(st:ed,1:1), stat=ierr)
            allocate(ve(st:ed,1:1), stat=ierr)
            allocate(dn(st:ed,1:1), stat=ierr)
            do j=1,nj
            do i=1,ni
                do k=st,ed
                    vn(k,1) = vn3d1(i,j,k)*vn3d2(i,j,k)
                end do

                nfs = fsfs(nfsf_ve)%i3d(i,j,1 )
                nfe = fsfe(nfsf_ve)%i3d(i,j,nk)
                call sub_ve(nk,1,1,ngn,vn,nfs,nfe,nge,ve)

                nfs = fsfs(nfsf_dn)%i3d(i,j,1 )
                nfe = fsfe(nfsf_dn)%i3d(i,j,nk)
                call sub_dn(nk,1,1,ngn,vn,nge,ve,nfs,nfe,dn)

                do k=1,nk
                    dn3d(i,j,k) = dn(k,1)
                end do
            end do
            end do
            deallocate(dn, stat=ierr)
            deallocate(ve, stat=ierr)
            deallocate(vn, stat=ierr)
        end if
    end do

end subroutine calc_cc_mb_dn_via_node2

subroutine calc_sp_mb_dn_via_node2(mb_vn1,nvn1, &
                                   mb_vn2,nvn2, &
                                   ngn,sub_ve,nge,sub_dn, &
                                   nfsf_ve,nfsf_dn, &
                                   mb_dn,ndn,ndir)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_datatypes, only : var_block_t,fld_array_t
    use mod_fieldvars, only : nblkcoms,blkcomssp,mb_fsfsp
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
    integer(kind_int) :: nc,nb,i,j,k,nfs,nfe,ierr
    integer(kind_int) :: ni,nj,nk,st,ed
    real(kind_real), pointer :: vn3d1(:,:,:),vn3d2(:,:,:)
    real(kind_real), pointer :: dn3d(:,:,:)
    real(kind_real), pointer :: vn(:,:),ve(:,:),dn(:,:)
    type(fld_array_t), pointer :: fsfs(:),fsfe(:)

    do nc=1,nblkcoms
        nb = blkcomssp(nc)%nb

        ni = blkcomssp(nc)%top%nijk(1)
        nj = blkcomssp(nc)%top%nijk(2)
        nk = blkcomssp(nc)%top%nijk(3)

        vn3d1 => mb_vn1(nb)%fld(nvn1)%r3d
        vn3d2 => mb_vn2(nb)%fld(nvn2)%r3d
        dn3d  => mb_dn(nb)%fld(ndn)%r3d

        if (ndir == 1) then
            fsfs => mb_fsfsp(nb,1)%fld
            fsfe => mb_fsfsp(nb,2)%fld

            st = 1 - ngn
            ed = ni + ngn
            allocate(vn(st:ed,1:1), stat=ierr)
            allocate(ve(st:ed,1:1), stat=ierr)
            allocate(dn(st:ed,1:1), stat=ierr)
            do k=1,nk
            do j=1,nj
                do i=st,ed
                    vn(i,1) = vn3d1(i,j,k)*vn3d2(i,j,k)
                end do

                nfs = fsfs(nfsf_ve)%i3d(1 ,j,k)
                nfe = fsfe(nfsf_ve)%i3d(ni,j,k)
                call sub_ve(ni,1,1,ngn,vn,nfs,nfe,nge,ve)

                nfs = fsfs(nfsf_dn)%i3d(1 ,j,k)
                nfe = fsfe(nfsf_dn)%i3d(ni,j,k)
                call sub_dn(ni,1,1,ngn,vn,nge,ve,nfs,nfe,dn)

                do i=1,ni
                    dn3d(i,j,k) = dn(i,1)
                end do
            end do
            end do
            deallocate(dn, stat=ierr)
            deallocate(ve, stat=ierr)
            deallocate(vn, stat=ierr)

        else if (ndir == 2) then
            fsfs => mb_fsfsp(nb,3)%fld
            fsfe => mb_fsfsp(nb,4)%fld

            st = 1 - ngn
            ed = nj + ngn
            allocate(vn(st:ed,1:1), stat=ierr)
            allocate(ve(st:ed,1:1), stat=ierr)
            allocate(dn(st:ed,1:1), stat=ierr)
            do k=1,nk
            do i=1,ni
                do j=st,ed
                    vn(j,1) = vn3d1(i,j,k)*vn3d2(i,j,k)
                end do

                nfs = fsfs(nfsf_ve)%i3d(i,1 ,k)
                nfe = fsfe(nfsf_ve)%i3d(i,nj,k)
                call sub_ve(nj,1,1,ngn,vn,nfs,nfe,nge,ve)

                nfs = fsfs(nfsf_dn)%i3d(i,1 ,k)
                nfe = fsfe(nfsf_dn)%i3d(i,nj,k)
                call sub_dn(nj,1,1,ngn,vn,nge,ve,nfs,nfe,dn)

                do j=1,nj
                    dn3d(i,j,k) = dn(j,1)
                end do
            end do
            end do
            deallocate(dn, stat=ierr)
            deallocate(ve, stat=ierr)
            deallocate(vn, stat=ierr)

        else if (ndir == 3) then
            fsfs => mb_fsfsp(nb,5)%fld
            fsfe => mb_fsfsp(nb,6)%fld

            st = 1 - ngn
            ed = nk + ngn
            allocate(vn(st:ed,1:1), stat=ierr)
            allocate(ve(st:ed,1:1), stat=ierr)
            allocate(dn(st:ed,1:1), stat=ierr)
            do j=1,nj
            do i=1,ni
                do k=st,ed
                    vn(k,1) = vn3d1(i,j,k)*vn3d2(i,j,k)
                end do

                nfs = fsfs(nfsf_ve)%i3d(i,j,1 )
                nfe = fsfe(nfsf_ve)%i3d(i,j,nk)
                call sub_ve(nk,1,1,ngn,vn,nfs,nfe,nge,ve)

                nfs = fsfs(nfsf_dn)%i3d(i,j,1 )
                nfe = fsfe(nfsf_dn)%i3d(i,j,nk)
                call sub_dn(nk,1,1,ngn,vn,nge,ve,nfs,nfe,dn)

                do k=1,nk
                    dn3d(i,j,k) = dn(k,1)
                end do
            end do
            end do
            deallocate(dn, stat=ierr)
            deallocate(ve, stat=ierr)
            deallocate(vn, stat=ierr)
        end if
    end do

end subroutine calc_sp_mb_dn_via_node2

subroutine calc_mb_dn_via_node2_exp(mb_vn1,nvn1, &
                                     mb_vn2,nvn2, &
                                     ngn,sub_ve,nge,sub_dn, &
                                     nfsf_ve,nfsf_dn, &
                                     mb_dn,ndn,ndir)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_inter_scheme,nbc_intbc_scheme
    use mod_datatypes, only : var_block_t,fld_array_t
    use mod_variables, only : nghnode
    use mod_fieldvars, only : nblkcoms,blkcoms,mb_fsf
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
    integer(kind_int) :: nc,nb,i,j,k,nfs,nfe,ierr
    integer(kind_int) :: ni,nj,nk,st,ed,st1,ed1,st2,ed2,st0,ed0,iflg,jflg,kflg
    real(kind_real), pointer :: vn3d1(:,:,:),vn3d2(:,:,:)
    real(kind_real), pointer :: dn3d(:,:,:)
    real(kind_real), pointer :: vn(:,:),ve(:,:),dn(:,:)
    type(fld_array_t), pointer :: fsfs(:),fsfe(:)

    do nc=1,nblkcoms
        nb = blkcoms(nc)%nb

        ni = blkcoms(nc)%top%nijk(1)
        nj = blkcoms(nc)%top%nijk(2)
        nk = blkcoms(nc)%top%nijk(3)

        vn3d1 => mb_vn1(nb)%fld(nvn1)%r3d
        vn3d2 => mb_vn2(nb)%fld(nvn2)%r3d
        dn3d  => mb_dn(nb)%fld(ndn)%r3d

        if (ndir == 1) then
            fsfs => mb_fsf(nb,1)%fld
            fsfe => mb_fsf(nb,2)%fld

            st = 1 - ngn
            ed = ni + ngn
            st1 = 1 - ngn + nghnode
            ed1 = nk + ngn - nghnode
            st2 = 1 - ngn + nghnode
            ed2 = nj + ngn - nghnode
            allocate(vn(st:ed,1:1), stat=ierr)
            allocate(ve(st:ed,1:1), stat=ierr)
            allocate(dn(st:ed,1:1), stat=ierr)
            do k=st1,ed1
            kflg=k
            if(kflg<1) then
                kflg=1
            else if(kflg>nk) then
                kflg=nk
            end if
            do j=st2,ed2
                jflg=j
                if(jflg<1) then
                    jflg=1
                else if(jflg>nj) then
                    jflg=nj
                end if
                do i=st,ed
                    vn(i,1) = vn3d1(i,j,k)*vn3d2(i,j,k)
                end do

                nfs = fsfs(nfsf_ve)%i3d(1 ,jflg,kflg)
                nfe = fsfe(nfsf_ve)%i3d(ni,jflg,kflg)
                call sub_ve(ni,1,1,ngn,vn,nfs,nfe,nge,ve)

                nfs = fsfs(nfsf_dn)%i3d(1 ,jflg,kflg)
                nfe = fsfe(nfsf_dn)%i3d(ni,jflg,kflg)
                call sub_dn(ni,1,1,ngn,vn,nge,ve,nfs,nfe,dn)

                if (nfs == nbc_inter_scheme .or. nfs == nbc_intbc_scheme) then
                    st0 = st+nghnode
                else
                    st0 = 1
                end if

                if (nfe == nbc_inter_scheme .or. nfe == nbc_intbc_scheme) then
                    ed0 = ed-nghnode
                else
                    ed0 = ni
                end if

                do i=st0,ed0
                    dn3d(i,j,k) = dn(i,1)
                end do
            end do
            end do
            deallocate(dn, stat=ierr)
            deallocate(ve, stat=ierr)
            deallocate(vn, stat=ierr)

        else if (ndir == 2) then
            fsfs => mb_fsf(nb,3)%fld
            fsfe => mb_fsf(nb,4)%fld

            st = 1 - ngn
            ed = nj + ngn
            st1 = 1 - ngn + nghnode
            ed1 = nk + ngn - nghnode
            st2 = 1 - ngn + nghnode
            ed2 = ni + ngn - nghnode
            allocate(vn(st:ed,1:1), stat=ierr)
            allocate(ve(st:ed,1:1), stat=ierr)
            allocate(dn(st:ed,1:1), stat=ierr)
            do k=st1,ed1
            kflg=k
            if(kflg<1) then
                kflg=1
            else if(kflg>nk) then
                kflg=nk
            end if
            do i=st2,ed2
                iflg=i
                if(iflg<1) then
                    iflg=1
                else if(iflg>ni) then
                    iflg=ni
                end if
                do j=st,ed
                    vn(j,1) = vn3d1(i,j,k)*vn3d2(i,j,k)
                end do

                nfs = fsfs(nfsf_ve)%i3d(iflg,1 ,kflg)
                nfe = fsfe(nfsf_ve)%i3d(iflg,nj,kflg)
                call sub_ve(nj,1,1,ngn,vn,nfs,nfe,nge,ve)

                nfs = fsfs(nfsf_dn)%i3d(iflg,1 ,kflg)
                nfe = fsfe(nfsf_dn)%i3d(iflg,nj,kflg)
                call sub_dn(nj,1,1,ngn,vn,nge,ve,nfs,nfe,dn)

                if (nfs == nbc_inter_scheme .or. nfs == nbc_intbc_scheme) then
                    st0 = st+nghnode
                else
                    st0 = 1
                end if

                if (nfe == nbc_inter_scheme .or. nfe == nbc_intbc_scheme) then
                    ed0 = ed-nghnode
                else
                    ed0 = nj
                end if

                do j=st0,ed0
                    dn3d(i,j,k) = dn(j,1)
                end do
            end do
            end do
            deallocate(dn, stat=ierr)
            deallocate(ve, stat=ierr)
            deallocate(vn, stat=ierr)

        else if (ndir == 3) then
            fsfs => mb_fsf(nb,5)%fld
            fsfe => mb_fsf(nb,6)%fld

            st = 1 - ngn
            ed = nk + ngn
            st1 = 1 - ngn + nghnode
            ed1 = nj + ngn - nghnode
            st2 = 1 - ngn + nghnode
            ed2 = ni + ngn - nghnode
            allocate(vn(st:ed,1:1), stat=ierr)
            allocate(ve(st:ed,1:1), stat=ierr)
            allocate(dn(st:ed,1:1), stat=ierr)
            do j=st1,ed1
            jflg=j
            if(jflg<1) then
                jflg=1
            else if(jflg>nj) then
                jflg=nj
            end if
            do i=st2,ed2
                iflg=i
                if(iflg<1) then
                    iflg=1
                else if(iflg>ni) then
                    iflg=ni
                end if
                do k=st,ed
                    vn(k,1) = vn3d1(i,j,k)*vn3d2(i,j,k)
                end do

                nfs = fsfs(nfsf_ve)%i3d(iflg,jflg,1 )
                nfe = fsfe(nfsf_ve)%i3d(iflg,jflg,nk)
                call sub_ve(nk,1,1,ngn,vn,nfs,nfe,nge,ve)

                nfs = fsfs(nfsf_dn)%i3d(iflg,jflg,1 )
                nfe = fsfe(nfsf_dn)%i3d(iflg,jflg,nk)
                call sub_dn(nk,1,1,ngn,vn,nge,ve,nfs,nfe,dn)

                if (nfs == nbc_inter_scheme .or. nfs == nbc_intbc_scheme) then
                    st0 = st+nghnode
                else
                    st0 = 1
                end if

                if (nfe == nbc_inter_scheme .or. nfe == nbc_intbc_scheme) then
                    ed0 = ed-nghnode
                else
                    ed0 = nk
                end if

                do k=st0,ed0
                    dn3d(i,j,k) = dn(k,1)
                end do
            end do
            end do
            deallocate(dn, stat=ierr)
            deallocate(ve, stat=ierr)
            deallocate(vn, stat=ierr)
        end if
    end do

end subroutine calc_mb_dn_via_node2_exp

subroutine calc_mb_dn_via_vec_node2(mb_vn1,nst1,ned1, &
                                    mb_vn2,nst2,ned2, &
                                    ngn,sub_ve,nge,sub_dn, &
                                    nfsf_ve,nfsf_dn ,&
                                    mb_dn,ndn,ndir)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_datatypes, only : var_block_t,fld_array_t
    use mod_fieldvars, only : nblkcoms,blkcoms,mb_fsf
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
    integer(kind_int) :: nc,nb,i,j,k,m,m2,nfs,nfe,ierr
    integer(kind_int) :: ni,nj,nk,st,ed
    type(fld_array_t), pointer :: vn1(:),vn2(:)
    real(kind_real), pointer :: dn3d(:,:,:)
    real(kind_real), pointer :: vn(:,:),ve(:,:),dn(:,:)
    type(fld_array_t), pointer :: fsfs(:),fsfe(:)

    ierr = (ned1-nst1) - (ned2-nst2)
    call error_check(ierr, &
                     "The size of array isn't 9 in subroutine calc_mb_dn_via_vec_node2")

    do nc=1,nblkcoms
        nb = blkcoms(nc)%nb

        ni = blkcoms(nc)%top%nijk(1)
        nj = blkcoms(nc)%top%nijk(2)
        nk = blkcoms(nc)%top%nijk(3)

        vn1 => mb_vn1(nb)%fld
        vn2 => mb_vn2(nb)%fld
        dn3d => mb_dn(nb)%fld(ndn)%r3d

        if (ndir == 1) then
            fsfs => mb_fsf(nb,1)%fld
            fsfe => mb_fsf(nb,2)%fld

            st = 1 - ngn
            ed = ni + ngn
            allocate(vn(st:ed,1:1), stat=ierr)
            allocate(ve(st:ed,1:1), stat=ierr)
            allocate(dn(st:ed,1:1), stat=ierr)
            do k=1,nk
            do j=1,nj
                do i=st,ed
                    vn(i,1) = vn1(nst1)%r3d(i,j,k)* &
                              vn2(nst2)%r3d(i,j,k)
                    do m=nst1+1,ned1
                       m2 = nst2+m-nst1
                       vn(i,1) = vn(i,1) + &
                                 vn1(m )%r3d(i,j,k)* &
                                 vn2(m2)%r3d(i,j,k)
                    end do
                end do

                nfs = fsfs(nfsf_ve)%i3d(1 ,j,k)
                nfe = fsfe(nfsf_ve)%i3d(ni,j,k)
                call sub_ve(ni,1,1,ngn,vn,nfs,nfe,nge,ve)

                nfs = fsfs(nfsf_dn)%i3d(1 ,j,k)
                nfe = fsfe(nfsf_dn)%i3d(ni,j,k)
                call sub_dn(ni,1,1,ngn,vn,nge,ve,nfs,nfe,dn)

                do i=1,ni
                    dn3d(i,j,k) = dn(i,1)
                end do
            end do
            end do
            deallocate(dn, stat=ierr)
            deallocate(ve, stat=ierr)
            deallocate(vn, stat=ierr)

        else if (ndir == 2) then
            fsfs => mb_fsf(nb,3)%fld
            fsfe => mb_fsf(nb,4)%fld

            st = 1 - ngn
            ed = nj + ngn
            allocate(vn(st:ed,1:1), stat=ierr)
            allocate(ve(st:ed,1:1), stat=ierr)
            allocate(dn(st:ed,1:1), stat=ierr)
            do k=1,nk
            do i=1,ni
                do j=st,ed
                    vn(j,1) = vn1(nst1)%r3d(i,j,k)* &
                              vn2(nst2)%r3d(i,j,k)
                    do m=nst1+1,ned1
                       m2 = nst2+m-nst1
                       vn(j,1) = vn(j,1) + &
                                 vn1(m )%r3d(i,j,k)* &
                                 vn2(m2)%r3d(i,j,k)
                    end do
                end do

                nfs = fsfs(nfsf_ve)%i3d(i,1 ,k)
                nfe = fsfe(nfsf_ve)%i3d(i,nj,k)
                call sub_ve(nj,1,1,ngn,vn,nfs,nfe,nge,ve)

                nfs = fsfs(nfsf_dn)%i3d(i,1 ,k)
                nfe = fsfe(nfsf_dn)%i3d(i,nj,k)
                call sub_dn(nj,1,1,ngn,vn,nge,ve,nfs,nfe,dn)

                do j=1,nj
                    dn3d(i,j,k) = dn(j,1)
                end do
            end do
            end do
            deallocate(dn, stat=ierr)
            deallocate(ve, stat=ierr)
            deallocate(vn, stat=ierr)

        else if (ndir == 3) then
            fsfs => mb_fsf(nb,5)%fld
            fsfe => mb_fsf(nb,6)%fld

            st = 1 - ngn
            ed = nk + ngn
            allocate(vn(st:ed,1:1), stat=ierr)
            allocate(ve(st:ed,1:1), stat=ierr)
            allocate(dn(st:ed,1:1), stat=ierr)
            do j=1,nj
            do i=1,ni
                do k=st,ed
                    vn(k,1) = vn1(nst1)%r3d(i,j,k)* &
                              vn2(nst2)%r3d(i,j,k)
                    do m=nst1+1,ned1
                       m2 = nst2+m-nst1
                       vn(k,1) = vn(k,1) + &
                                 vn1(m )%r3d(i,j,k)* &
                                 vn2(m2)%r3d(i,j,k)
                    end do
                end do

                nfs = fsfs(nfsf_ve)%i3d(i,j,1 )
                nfe = fsfe(nfsf_ve)%i3d(i,j,nk)
                call sub_ve(nk,1,1,ngn,vn,nfs,nfe,nge,ve)

                nfs = fsfs(nfsf_dn)%i3d(i,j,1 )
                nfe = fsfe(nfsf_dn)%i3d(i,j,nk)
                call sub_dn(nk,1,1,ngn,vn,nge,ve,nfs,nfe,dn)

                do k=1,nk
                    dn3d(i,j,k) = dn(k,1)
                end do
            end do
            end do
            deallocate(dn, stat=ierr)
            deallocate(ve, stat=ierr)
            deallocate(vn, stat=ierr)
        end if
    end do

end subroutine calc_mb_dn_via_vec_node2

subroutine calc_cc_mb_dn_via_vec_node2(mb_vn1,nst1,ned1, &
                                       mb_vn2,nst2,ned2, &
                                       ngn,sub_ve,nge,sub_dn, &
                                       nfsf_ve,nfsf_dn ,&
                                       mb_dn,ndn,ndir)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_datatypes, only : var_block_t,fld_array_t
    use mod_fieldvars, only : nblkcoms,blkcomscc,mb_fsfcc
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
    integer(kind_int) :: nc,nb,i,j,k,m,m2,nfs,nfe,ierr
    integer(kind_int) :: ni,nj,nk,st,ed
    type(fld_array_t), pointer :: vn1(:),vn2(:)
    real(kind_real), pointer :: dn3d(:,:,:)
    real(kind_real), pointer :: vn(:,:),ve(:,:),dn(:,:)
    type(fld_array_t), pointer :: fsfs(:),fsfe(:)

    ierr = (ned1-nst1) - (ned2-nst2)
    call error_check(ierr, &
                     "The size of array isn't 9 in subroutine calc_mb_dn_via_vec_node2")

    do nc=1,nblkcoms
        nb = blkcomscc(nc)%nb

        ni = blkcomscc(nc)%top%nijk(1)
        nj = blkcomscc(nc)%top%nijk(2)
        nk = blkcomscc(nc)%top%nijk(3)

        vn1 => mb_vn1(nb)%fld
        vn2 => mb_vn2(nb)%fld
        dn3d => mb_dn(nb)%fld(ndn)%r3d

        if (ndir == 1) then
            fsfs => mb_fsfcc(nb,1)%fld
            fsfe => mb_fsfcc(nb,2)%fld

            st = 1 - ngn
            ed = ni + ngn
            allocate(vn(st:ed,1:1), stat=ierr)
            allocate(ve(st:ed,1:1), stat=ierr)
            allocate(dn(st:ed,1:1), stat=ierr)
            do k=1,nk
            do j=1,nj
                do i=st,ed
                    vn(i,1) = vn1(nst1)%r3d(i,j,k)* &
                              vn2(nst2)%r3d(i,j,k)
                    do m=nst1+1,ned1
                       m2 = nst2+m-nst1
                       vn(i,1) = vn(i,1) + &
                                 vn1(m )%r3d(i,j,k)* &
                                 vn2(m2)%r3d(i,j,k)
                    end do
                end do

                nfs = fsfs(nfsf_ve)%i3d(1 ,j,k)
                nfe = fsfe(nfsf_ve)%i3d(ni,j,k)
                call sub_ve(ni,1,1,ngn,vn,nfs,nfe,nge,ve)

                nfs = fsfs(nfsf_dn)%i3d(1 ,j,k)
                nfe = fsfe(nfsf_dn)%i3d(ni,j,k)
                call sub_dn(ni,1,1,ngn,vn,nge,ve,nfs,nfe,dn)

                do i=1,ni
                    dn3d(i,j,k) = dn(i,1)
                end do
            end do
            end do
            deallocate(dn, stat=ierr)
            deallocate(ve, stat=ierr)
            deallocate(vn, stat=ierr)

        else if (ndir == 2) then
            fsfs => mb_fsfcc(nb,3)%fld
            fsfe => mb_fsfcc(nb,4)%fld

            st = 1 - ngn
            ed = nj + ngn
            allocate(vn(st:ed,1:1), stat=ierr)
            allocate(ve(st:ed,1:1), stat=ierr)
            allocate(dn(st:ed,1:1), stat=ierr)
            do k=1,nk
            do i=1,ni
                do j=st,ed
                    vn(j,1) = vn1(nst1)%r3d(i,j,k)* &
                              vn2(nst2)%r3d(i,j,k)
                    do m=nst1+1,ned1
                       m2 = nst2+m-nst1
                       vn(j,1) = vn(j,1) + &
                                 vn1(m )%r3d(i,j,k)* &
                                 vn2(m2)%r3d(i,j,k)
                    end do
                end do

                nfs = fsfs(nfsf_ve)%i3d(i,1 ,k)
                nfe = fsfe(nfsf_ve)%i3d(i,nj,k)
                call sub_ve(nj,1,1,ngn,vn,nfs,nfe,nge,ve)

                nfs = fsfs(nfsf_dn)%i3d(i,1 ,k)
                nfe = fsfe(nfsf_dn)%i3d(i,nj,k)
                call sub_dn(nj,1,1,ngn,vn,nge,ve,nfs,nfe,dn)

                do j=1,nj
                    dn3d(i,j,k) = dn(j,1)
                end do
            end do
            end do
            deallocate(dn, stat=ierr)
            deallocate(ve, stat=ierr)
            deallocate(vn, stat=ierr)

        else if (ndir == 3) then
            fsfs => mb_fsfcc(nb,5)%fld
            fsfe => mb_fsfcc(nb,6)%fld

            st = 1 - ngn
            ed = nk + ngn
            allocate(vn(st:ed,1:1), stat=ierr)
            allocate(ve(st:ed,1:1), stat=ierr)
            allocate(dn(st:ed,1:1), stat=ierr)
            do j=1,nj
            do i=1,ni
                do k=st,ed
                    vn(k,1) = vn1(nst1)%r3d(i,j,k)* &
                              vn2(nst2)%r3d(i,j,k)
                    do m=nst1+1,ned1
                       m2 = nst2+m-nst1
                       vn(k,1) = vn(k,1) + &
                                 vn1(m )%r3d(i,j,k)* &
                                 vn2(m2)%r3d(i,j,k)
                    end do
                end do

                nfs = fsfs(nfsf_ve)%i3d(i,j,1 )
                nfe = fsfe(nfsf_ve)%i3d(i,j,nk)
                call sub_ve(nk,1,1,ngn,vn,nfs,nfe,nge,ve)

                nfs = fsfs(nfsf_dn)%i3d(i,j,1 )
                nfe = fsfe(nfsf_dn)%i3d(i,j,nk)
                call sub_dn(nk,1,1,ngn,vn,nge,ve,nfs,nfe,dn)

                do k=1,nk
                    dn3d(i,j,k) = dn(k,1)
                end do
            end do
            end do
            deallocate(dn, stat=ierr)
            deallocate(ve, stat=ierr)
            deallocate(vn, stat=ierr)
        end if
    end do

end subroutine calc_cc_mb_dn_via_vec_node2

subroutine calc_sp_mb_dn_via_vec_node2(mb_vn1,nst1,ned1, &
                                       mb_vn2,nst2,ned2, &
                                       ngn,sub_ve,nge,sub_dn, &
                                       nfsf_ve,nfsf_dn ,&
                                       mb_dn,ndn,ndir)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_datatypes, only : var_block_t,fld_array_t
    use mod_fieldvars, only : nblkcoms,blkcomssp,mb_fsfsp
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
    integer(kind_int) :: nc,nb,i,j,k,m,m2,nfs,nfe,ierr
    integer(kind_int) :: ni,nj,nk,st,ed
    type(fld_array_t), pointer :: vn1(:),vn2(:)
    real(kind_real), pointer :: dn3d(:,:,:)
    real(kind_real), pointer :: vn(:,:),ve(:,:),dn(:,:)
    type(fld_array_t), pointer :: fsfs(:),fsfe(:)

    ierr = (ned1-nst1) - (ned2-nst2)
    call error_check(ierr, &
                     "The size of array isn't 9 in subroutine calc_mb_dn_via_vec_node2")

    do nc=1,nblkcoms
        nb = blkcomssp(nc)%nb

        ni = blkcomssp(nc)%top%nijk(1)
        nj = blkcomssp(nc)%top%nijk(2)
        nk = blkcomssp(nc)%top%nijk(3)

        vn1 => mb_vn1(nb)%fld
        vn2 => mb_vn2(nb)%fld
        dn3d => mb_dn(nb)%fld(ndn)%r3d

        if (ndir == 1) then
            fsfs => mb_fsfsp(nb,1)%fld
            fsfe => mb_fsfsp(nb,2)%fld

            st = 1 - ngn
            ed = ni + ngn
            allocate(vn(st:ed,1:1), stat=ierr)
            allocate(ve(st:ed,1:1), stat=ierr)
            allocate(dn(st:ed,1:1), stat=ierr)
            do k=1,nk
            do j=1,nj
                do i=st,ed
                    vn(i,1) = vn1(nst1)%r3d(i,j,k)* &
                              vn2(nst2)%r3d(i,j,k)
                    do m=nst1+1,ned1
                       m2 = nst2+m-nst1
                       vn(i,1) = vn(i,1) + &
                                 vn1(m )%r3d(i,j,k)* &
                                 vn2(m2)%r3d(i,j,k)
                    end do
                end do

                nfs = fsfs(nfsf_ve)%i3d(1 ,j,k)
                nfe = fsfe(nfsf_ve)%i3d(ni,j,k)
                call sub_ve(ni,1,1,ngn,vn,nfs,nfe,nge,ve)

                nfs = fsfs(nfsf_dn)%i3d(1 ,j,k)
                nfe = fsfe(nfsf_dn)%i3d(ni,j,k)
                call sub_dn(ni,1,1,ngn,vn,nge,ve,nfs,nfe,dn)

                do i=1,ni
                    dn3d(i,j,k) = dn(i,1)
                end do
            end do
            end do
            deallocate(dn, stat=ierr)
            deallocate(ve, stat=ierr)
            deallocate(vn, stat=ierr)

        else if (ndir == 2) then
            fsfs => mb_fsfsp(nb,3)%fld
            fsfe => mb_fsfsp(nb,4)%fld

            st = 1 - ngn
            ed = nj + ngn
            allocate(vn(st:ed,1:1), stat=ierr)
            allocate(ve(st:ed,1:1), stat=ierr)
            allocate(dn(st:ed,1:1), stat=ierr)
            do k=1,nk
            do i=1,ni
                do j=st,ed
                    vn(j,1) = vn1(nst1)%r3d(i,j,k)* &
                              vn2(nst2)%r3d(i,j,k)
                    do m=nst1+1,ned1
                       m2 = nst2+m-nst1
                       vn(j,1) = vn(j,1) + &
                                 vn1(m )%r3d(i,j,k)* &
                                 vn2(m2)%r3d(i,j,k)
                    end do
                end do

                nfs = fsfs(nfsf_ve)%i3d(i,1 ,k)
                nfe = fsfe(nfsf_ve)%i3d(i,nj,k)
                call sub_ve(nj,1,1,ngn,vn,nfs,nfe,nge,ve)

                nfs = fsfs(nfsf_dn)%i3d(i,1 ,k)
                nfe = fsfe(nfsf_dn)%i3d(i,nj,k)
                call sub_dn(nj,1,1,ngn,vn,nge,ve,nfs,nfe,dn)

                do j=1,nj
                    dn3d(i,j,k) = dn(j,1)
                end do
            end do
            end do
            deallocate(dn, stat=ierr)
            deallocate(ve, stat=ierr)
            deallocate(vn, stat=ierr)

        else if (ndir == 3) then
            fsfs => mb_fsfsp(nb,5)%fld
            fsfe => mb_fsfsp(nb,6)%fld

            st = 1 - ngn
            ed = nk + ngn
            allocate(vn(st:ed,1:1), stat=ierr)
            allocate(ve(st:ed,1:1), stat=ierr)
            allocate(dn(st:ed,1:1), stat=ierr)
            do j=1,nj
            do i=1,ni
                do k=st,ed
                    vn(k,1) = vn1(nst1)%r3d(i,j,k)* &
                              vn2(nst2)%r3d(i,j,k)
                    do m=nst1+1,ned1
                       m2 = nst2+m-nst1
                       vn(k,1) = vn(k,1) + &
                                 vn1(m )%r3d(i,j,k)* &
                                 vn2(m2)%r3d(i,j,k)
                    end do
                end do

                nfs = fsfs(nfsf_ve)%i3d(i,j,1 )
                nfe = fsfe(nfsf_ve)%i3d(i,j,nk)
                call sub_ve(nk,1,1,ngn,vn,nfs,nfe,nge,ve)

                nfs = fsfs(nfsf_dn)%i3d(i,j,1 )
                nfe = fsfe(nfsf_dn)%i3d(i,j,nk)
                call sub_dn(nk,1,1,ngn,vn,nge,ve,nfs,nfe,dn)

                do k=1,nk
                    dn3d(i,j,k) = dn(k,1)
                end do
            end do
            end do
            deallocate(dn, stat=ierr)
            deallocate(ve, stat=ierr)
            deallocate(vn, stat=ierr)
        end if
    end do

end subroutine calc_sp_mb_dn_via_vec_node2

subroutine calc_mb_var_via_sub(mb_vin,nsin,nein,sub_via, &
                               mb_vout,nsout,neout,ngh)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_datatypes, only : var_block_t,fld_array_t
    use mod_fieldvars, only : nblkcoms,blkcoms
    use mod_openmp
    implicit none
    type(var_block_t), pointer, intent(in) :: mb_vin(:)
    integer(kind_int)         , intent(in) :: nsin,nein
    external                               :: sub_via
    type(var_block_t), pointer, intent(in) :: mb_vout(:)
    integer(kind_int)         , intent(in) :: nsout,neout
    integer(kind_int)         , intent(in) :: ngh
    integer(kind_int)          :: nc,nb,i,j,k,m,ierr
    integer(kind_int)          :: st(3),ed(3)
    real(kind_real)            :: vin(nsin:nein)
    real(kind_real)            :: vout(nsout:neout)
    type(fld_array_t), pointer :: fin(:),fout(:)

    do nc=1,nblkcoms
        nb = blkcoms(nc)%nb

        st(:) = 1 - ngh
        ed(:) = blkcoms(nc)%top%nijk(:) + ngh

        fin  => mb_vin(nb)%fld
        fout => mb_vout(nb)%fld
!$OMP parallel do private(vin,vout)
        do k=st(3),ed(3)
        do j=st(2),ed(2)
        do i=st(1),ed(1)
            do m=nsin,nein
               vin(m) = fin(m)%r3d(i,j,k)
            end do

            call sub_via(nsin,nein,vin,nsout,neout,vout)

            do m=nsout,neout
               fout(m)%r3d(i,j,k) = vout(m)
            end do
        end do
        end do
        end do
!$OMP end parallel do
    end do

end subroutine calc_mb_var_via_sub

subroutine calc_mb_var_via_sub_sp(mb_vin,nsin,nein,sub_via, &
                                  mb_vout,nsout,neout,ngh)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_datatypes, only : var_block_t,fld_array_t
    use mod_fieldvars, only : nblkcoms,blkcomssp
    use mod_openmp
    implicit none
    type(var_block_t), pointer, intent(in) :: mb_vin(:)
    integer(kind_int)         , intent(in) :: nsin,nein
    external                               :: sub_via
    type(var_block_t), pointer, intent(in) :: mb_vout(:)
    integer(kind_int)         , intent(in) :: nsout,neout
    integer(kind_int)         , intent(in) :: ngh
    integer(kind_int)          :: nc,nb,i,j,k,m,ierr
    integer(kind_int)          :: st(3),ed(3)
    real(kind_real)            :: vin(nsin:nein)
    real(kind_real)            :: vout(nsout:neout)
    type(fld_array_t), pointer :: fin(:),fout(:)

    do nc=1,nblkcoms
        nb = blkcomssp(nc)%nb

        st(:) = 1 - ngh
        ed(:) = blkcomssp(nc)%top%nijk(:) + ngh

        fin  => mb_vin(nb)%fld
        fout => mb_vout(nb)%fld
!$OMP parallel do private(vin,vout)
        do k=st(3),ed(3)
        do j=st(2),ed(2)
        do i=st(1),ed(1)
            do m=nsin,nein
               vin(m) = fin(m)%r3d(i,j,k)
            end do

            call sub_via(nsin,nein,vin,nsout,neout,vout)

            do m=nsout,neout
               fout(m)%r3d(i,j,k) = vout(m)
            end do
        end do
        end do
        end do
!$OMP end parallel do
    end do

end subroutine calc_mb_var_via_sub_sp

subroutine calc_bc_var_via_sub(mb_vin,nsin,nein,sub_via, &
                               mb_vout,nsout,neout,ngst,nged)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_datatypes, only : var_block_t,fld_array_t
    use mod_datatypes, only : top_block_t,bc_region_t
    use mod_fieldvars, only : nblkcoms,blkcoms
    use mod_openmp
    implicit none
    type(var_block_t), pointer, intent(in) :: mb_vin(:)
    integer(kind_int)         , intent(in) :: nsin,nein
    external                               :: sub_via
    type(var_block_t), pointer, intent(in) :: mb_vout(:)
    integer(kind_int)         , intent(in) :: nsout,neout
    integer(kind_int)         , intent(in) :: ngst,nged
    integer(kind_int)          :: nc,nb,nr,i,j,k,m,ierr
    integer(kind_int)          :: s_st(3),s_ed(3),st(3),ed(3)
    integer(kind_int)          :: nregs,s_nd,s_lr,bctype
    real(kind_real)            :: vin(nsin:nein)
    real(kind_real)            :: vout(nsout:neout)
    type(top_block_t), pointer :: top
    type(bc_region_t), pointer :: reg
    type(fld_array_t), pointer :: fin(:),fout(:)

    do nc=1,nblkcoms
        nb  =  blkcoms(nc)%nb
        top => blkcoms(nc)%top

        fin  => mb_vin(nb)%fld
        fout => mb_vout(nb)%fld

        nregs = top%nregions
        do nr=1,nregs
            reg => top%bcs(nr)

            bctype  = reg%bctype
            s_st(:) = reg%s_st(:)
            s_ed(:) = reg%s_ed(:)
            s_nd    = reg%s_nd
            s_lr    = reg%s_lr

            call bc_extend_ghost(s_st,s_ed,s_nd,s_lr,ngst,nged,st,ed)

!$OMP parallel do private(vin,vout)
            do k=st(3),ed(3)
            do j=st(2),ed(2)
            do i=st(1),ed(1)
                do m=nsin,nein
                   vin(m) = fin(m)%r3d(i,j,k)
                end do

                call sub_via(nsin,nein,vin,nsout,neout,vout)

                do m=nsout,neout
                   fout(m)%r3d(i,j,k) = vout(m)
                end do
            end do
            end do
            end do
!$OMP end parallel do
        end do
    end do

end subroutine calc_bc_var_via_sub

subroutine calc_bc_var_via_sub_sp(mb_vin,nsin,nein,sub_via, &
                                  mb_vout,nsout,neout,ngst,nged)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_datatypes, only : var_block_t,fld_array_t
    use mod_datatypes, only : top_block_t,bc_region_t
    use mod_fieldvars, only : nblkcoms,blkcomssp
    use mod_openmp
    implicit none
    type(var_block_t), pointer, intent(in) :: mb_vin(:)
    integer(kind_int)         , intent(in) :: nsin,nein
    external                               :: sub_via
    type(var_block_t), pointer, intent(in) :: mb_vout(:)
    integer(kind_int)         , intent(in) :: nsout,neout
    integer(kind_int)         , intent(in) :: ngst,nged
    integer(kind_int)          :: nc,nb,nr,i,j,k,m,ierr
    integer(kind_int)          :: s_st(3),s_ed(3),st(3),ed(3)
    integer(kind_int)          :: nregs,s_nd,s_lr,bctype
    real(kind_real)            :: vin(nsin:nein)
    real(kind_real)            :: vout(nsout:neout)
    type(top_block_t), pointer :: top
    type(bc_region_t), pointer :: reg
    type(fld_array_t), pointer :: fin(:),fout(:)

    do nc=1,nblkcoms
        nb  =  blkcomssp(nc)%nb
        top => blkcomssp(nc)%top

        fin  => mb_vin(nb)%fld
        fout => mb_vout(nb)%fld

        nregs = top%nregions
        do nr=1,nregs
            reg => top%bcs(nr)

            bctype  = reg%bctype
            s_st(:) = reg%s_st(:)
            s_ed(:) = reg%s_ed(:)
            s_nd    = reg%s_nd
            s_lr    = reg%s_lr

            call bc_extend_ghost(s_st,s_ed,s_nd,s_lr,ngst,nged,st,ed)

!$OMP parallel do private(vin,vout)
            do k=st(3),ed(3)
            do j=st(2),ed(2)
            do i=st(1),ed(1)
                do m=nsin,nein
                   vin(m) = fin(m)%r3d(i,j,k)
                end do

                call sub_via(nsin,nein,vin,nsout,neout,vout)

                do m=nsout,neout
                   fout(m)%r3d(i,j,k) = vout(m)
                end do
            end do
            end do
            end do
!$OMP end parallel do
        end do
    end do

end subroutine calc_bc_var_via_sub_sp


