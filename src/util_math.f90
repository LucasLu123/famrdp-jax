subroutine ve_via_node2(ni,nst,ned,ngn,vn,nfs,nfe,nge,ve)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_inter_scheme,nbc_intbc_scheme
    implicit none
    integer(kind_int),intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  ,intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)  :: nfs,nfe
    integer(kind_int),intent(in)  :: nge
    real(kind_real)  ,intent(out) :: ve(1-ngn:ni+ngn,nst:ned)
    integer(kind_int) :: i,m,st,ed
    
    if (nfs == nbc_inter_scheme) then
        st = -nge
    else if (nfs == nbc_intbc_scheme) then
        st = 1
    
        do m=nst,ned
            ve(0,m) = ( 3.0*vn(1,m) - vn(2,m) ) / 2.0
            ve(-1,m) = ( 5.0*vn(1,m) - 3.0*vn(2,m) ) / 2.0
            ve(-2,m) = ( 7.0*vn(1,m) - 5.0*vn(2,m) ) / 2.0
        end do
    else
        st = 1

        do m=nst,ned
            ve(0,m) = ( 3.0*vn(1,m) - vn(2,m) ) / 2.0
        end do
    end if

    if (nfe == nbc_inter_scheme) then
        ed = ni + nge
    else if (nfe == nbc_intbc_scheme) then
        ed = ni - 1
    
        do m=nst,ned
            ve(ni,m) = ( 3.0*vn(ni,m) - vn(ni-1,m) ) / 2.0
            ve(ni+1,m) = ( 5.0*vn(ni,m) - 3.0*vn(ni-1,m) ) / 2.0
            ve(ni+2,m) = ( 7.0*vn(ni,m) - 5.0*vn(ni-1,m) ) / 2.0
        end do
    else
        ed = ni - 1

        do m=nst,ned
            ve(ni,m) = ( 3.0*vn(ni,m) - vn(ni-1,m) ) / 2.0
        end do
    end if
    
    do m=nst,ned
        do i=st,ed
            ve(i,m) = ( vn(i+1,m) + vn(i,m) ) / 2.0
        end do
    end do
    
end subroutine ve_via_node2

subroutine ve_via_node2cc(ni,nst,ned,ngn,vn,nfs,nfe,nge,ve)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_inter_scheme,nbc_intbc_scheme
    implicit none
    integer(kind_int),intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  ,intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)  :: nfs,nfe
    integer(kind_int),intent(in)  :: nge
    real(kind_real)  ,intent(out) :: ve(1-ngn:ni+ngn,nst:ned)
    integer(kind_int) :: i,m,st,ed
    
    if (nfs == nbc_inter_scheme) then
        st = 1 - nge
    else if (nfs == nbc_intbc_scheme) then
        st = 2
    
        do m=nst,ned
            ve(1,m) = ( 3.0*vn(1,m) - vn(2,m) ) / 2.0
            ve(0,m) = ( 5.0*vn(1,m) - 3.0*vn(2,m) ) / 2.0
            ve(-1,m) = ( 7.0*vn(1,m) - 5.0*vn(2,m) ) / 2.0
        end do
    else
        st = 2

        do m=nst,ned
            ve(1,m) = ( 3.0*vn(1,m) - vn(2,m) ) / 2.0
        end do
    end if

    if (nfe == nbc_inter_scheme) then
        ed = ni + 1 + nge
    else if (nfe == nbc_intbc_scheme) then
        ed = ni
    
        do m=nst,ned
            ve(ni+1,m) = ( 3.0*vn(ni,m) - vn(ni-1,m) ) / 2.0
            ve(ni+2,m) = ( 5.0*vn(ni,m) - 3.0*vn(ni-1,m) ) / 2.0
            ve(ni+3,m) = ( 7.0*vn(ni,m) - 5.0*vn(ni-1,m) ) / 2.0
        end do
    else
        ed = ni

        do m=nst,ned
            ve(ni+1,m) = ( 3.0*vn(ni,m) - vn(ni-1,m) ) / 2.0
        end do
    end if
    
    do m=nst,ned
        do i=st-1,ed-1
            ve(i+1,m) = ( vn(i+1,m) + vn(i,m) ) / 2.0
        end do
    end do
    
end subroutine ve_via_node2cc

subroutine ve_via_node4(ni,nst,ned,ngn,vn,nfs,nfe,nge,ve)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_inter_scheme
    implicit none
    integer(kind_int),intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  ,intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)  :: nfs,nfe
    integer(kind_int),intent(in)  :: nge
    real(kind_real)  ,intent(out) :: ve(1-ngn:ni+ngn,nst:ned)
    integer(kind_int) :: i,m,st,ed
    
    if (nfs == nbc_inter_scheme) then
        st = -nge
    else
        st = 2

        do m=nst,ned
            !!ve(1,m) = (  5.0*vn(1,m) + 15.0*vn(2,m) -  5.0*vn(3,m) +     vn(4,m) ) / 16.0   ! fourth-order
            !!ve(0,m) = ( 35.0*vn(1,m) - 35.0*vn(2,m) + 21.0*vn(3,m) - 5.0*vn(4,m) ) / 16.0   ! fourth-order
            ve(1,m) = (  3.0*vn(1,m) +  6.0*vn(2,m) -      vn(3,m) ) / 8.0                  ! third-order
            ve(0,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                  ! third-order
        end do
    end if

    if (nfe == nbc_inter_scheme) then
        ed = ni + nge
    else
        ed = ni - 2

        do m=nst,ned
            !!ve(ni-1,m) = (  5.0*vn(ni,m) + 15.0*vn(ni-1,m) -  5.0*vn(ni-2,m) +     vn(ni-3,m) ) / 16.0  ! fourth-order
            !!ve(ni  ,m) = ( 35.0*vn(ni,m) - 35.0*vn(ni-1,m) + 21.0*vn(ni-2,m) - 5.0*vn(ni-3,m) ) / 16.0  ! fourth-order
            ve(ni-1,m) = (  3.0*vn(ni,m) +  6.0*vn(ni-1,m) -      vn(ni-2,m) ) / 8.0                ! third-order
            ve(ni  ,m) = ( 15.0*vn(ni,m) - 10.0*vn(ni-1,m) +  3.0*vn(ni-2,m) ) / 8.0                ! third-order
        end do
    end if
    
    do m=nst,ned
        do i=st,ed
            ve(i,m) = ( 9.0*(vn(i+1,m) + vn(i  ,m)) - &
                            (vn(i+2,m) + vn(i-1,m)) )/16.0
        end do
    end do
    
end subroutine ve_via_node4

subroutine ve_via_node4cc(ni,nst,ned,ngn,vn,nfs,nfe,nge,ve)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_inter_scheme
    implicit none
    integer(kind_int),intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  ,intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)  :: nfs,nfe
    integer(kind_int),intent(in)  :: nge
    real(kind_real)  ,intent(out) :: ve(1-ngn:ni+ngn,nst:ned)
    integer(kind_int) :: i,m,st,ed
    
    if (nfs == nbc_inter_scheme) then
        st = 1 - nge
    else
        st = 3

        do m=nst,ned
            !!ve(2,m) = (  5.0*vn(1,m) + 15.0*vn(2,m) -  5.0*vn(3,m) +     vn(4,m) ) / 16.0   ! fourth-order
            !!ve(1,m) = ( 35.0*vn(1,m) - 35.0*vn(2,m) + 21.0*vn(3,m) - 5.0*vn(4,m) ) / 16.0   ! fourth-order
            ve(2,m) = (  3.0*vn(1,m) +  6.0*vn(2,m) -      vn(3,m) ) / 8.0                  ! third-order
            ve(1,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                  ! third-order
        end do
    end if

    if (nfe == nbc_inter_scheme) then
        ed = ni + 1 + nge
    else
        ed = ni - 1

        do m=nst,ned
            !!ve(ni  ,m) = (  5.0*vn(ni,m) + 15.0*vn(ni-1,m) -  5.0*vn(ni-2,m) +     vn(ni-3,m) ) / 16.0  ! fourth-order
            !!ve(ni+1,m) = ( 35.0*vn(ni,m) - 35.0*vn(ni-1,m) + 21.0*vn(ni-2,m) - 5.0*vn(ni-3,m) ) / 16.0  ! fourth-order
            ve(ni  ,m) = (  3.0*vn(ni,m) +  6.0*vn(ni-1,m) -      vn(ni-2,m) ) / 8.0                ! third-order
            ve(ni+1,m) = ( 15.0*vn(ni,m) - 10.0*vn(ni-1,m) +  3.0*vn(ni-2,m) ) / 8.0                ! third-order
        end do
    end if
    
    do m=nst,ned
        do i=st-1,ed-1
            ve(i+1,m) = ( 9.0*(vn(i+1,m) + vn(i  ,m)) - &
                              (vn(i+2,m) + vn(i-1,m)) )/16.0
        end do
    end do
    
end subroutine ve_via_node4cc

subroutine ve_via_node6(ni,nst,ned,ngn,vn,nfs,nfe,nge,ve)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_inter_scheme,nbc_intbc_scheme
    implicit none
    integer(kind_int),intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  ,intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)  :: nfs,nfe
    integer(kind_int),intent(in)  :: nge
    real(kind_real)  ,intent(out) :: ve(1-ngn:ni+ngn,nst:ned)
    integer(kind_int) :: i,m,st,ed
    
    if (ni < 4) then
        call ve_via_node2(ni,nst,ned,ngn,vn,nfs,nfe,nge,ve)
        return
    end if

    if (nfs == nbc_inter_scheme) then
        st = -nge
    
    else if (nfs == nbc_intbc_scheme) then
    
        st = 3
        
        do m=nst,ned
        
            ve(-2,m) = (3003.0*vn(1,m) - 8580.0*vn(2,m) + 10010.0*vn(3,m) - 5460.0*vn(4,m) + 1155.0*vn(5,m) ) / 128.0   ! fifth-order
            ve(-1,m) = (1155.0*vn(1,m) - 2772.0*vn(2,m) +  2970.0*vn(3,m) - 1540.0*vn(4,m) +  315.0*vn(5,m) ) / 128.0   ! fifth-order
            ve(0 ,m) = ( 315.0*vn(1,m) -  420.0*vn(2,m) +   378.0*vn(3,m) -  180.0*vn(4,m) +   35.0*vn(5,m) ) / 128.0   ! fifth-order
            ve(1 ,m) = (  35.0*vn(1,m) +  140.0*vn(2,m) -    70.0*vn(3,m) +   28.0*vn(4,m) -    5.0*vn(5,m) ) / 128.0   ! fifth-order
            ve(2 ,m) = (      -vn(1,m) +    9.0*vn(2,m) +     9.0*vn(3,m) -        vn(4,m) ) / 16.0                     ! fourth-order
!!            ve(-2,m) = ( 63.0*vn(1,m) - 90.0*vn(2,m) +  35.0*vn(3,m) ) / 8.0                                       ! third-order
!!            ve(-1,m) = ( 35.0*vn(1,m) - 42.0*vn(2,m) +  15.0*vn(3,m) ) / 8.0                                       ! third-order
!!            ve(0 ,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +   3.0*vn(3,m) ) / 8.0                                       ! third-order
!!            ve(1 ,m) = (  3.0*vn(1,m) +  6.0*vn(2,m) -       vn(3,m) ) / 8.0                                       ! third-order
!!            ve(2 ,m) = (     -vn(1,m) +  9.0*vn(2,m) +   9.0*vn(3,m) - vn(4,m) ) / 16.0                            ! fourth-order
            
        end do
    
    else
        st = 3

        do m=nst,ned
            !!ve(2,m) = (      -vn(1,m) +   9.0*vn(2,m) +   9.0*vn(3,m) -       vn(4,m) ) / 16.0                   ! fourth-order
            !!ve(1,m) = (   5.0*vn(1,m) +  15.0*vn(2,m) -   5.0*vn(3,m) +       vn(4,m) ) / 16.0                   ! fourth-order
            !!ve(0,m) = (  35.0*vn(1,m) -  35.0*vn(2,m) +  21.0*vn(3,m) -   5.0*vn(4,m) ) / 16.0                   ! fourth-order
            ve(0,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                                       ! third-order
            ve(2,m) = (  -5.0*vn(1,m) +  60.0*vn(2,m) +  90.0*vn(3,m) -  20.0*vn(4,m) +  3.0*vn(5,m) ) / 128.0   ! fifth-order
            ve(1,m) = (  35.0*vn(1,m) + 140.0*vn(2,m) -  70.0*vn(3,m) +  28.0*vn(4,m) -  5.0*vn(5,m) ) / 128.0   ! fifth-order
            !!ve(0,m) = ( 315.0*vn(1,m) - 420.0*vn(2,m) + 378.0*vn(3,m) - 180.0*vn(4,m) + 35.0*vn(5,m) ) / 128.0   ! fifth-order
            !!ve(1,m) = ( 3.0*vn(1,m) + 6.0*vn(2,m) -     vn(3,m) ) / 8.0                                          ! third-order
            !!ve(2,m) = (    -vn(1,m) + 9.0*vn(2,m) + 9.0*vn(3,m) - vn(4,m) ) / 16.0                               ! fourth-order
        end do
    end if
    
    if (nfe == nbc_inter_scheme) then
        ed = ni + nge
    
    else if (nfe == nbc_intbc_scheme) then
    
        ed = ni - 3

        do m=nst,ned
            ve(ni+2,m) = &
            ( 3003.0*vn(ni,m) - 8580.0*vn(ni-1,m) + 10010.0*vn(ni-2,m) - 5460.0*vn(ni-3,m) + 1155.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            ve(ni+1,m) = &
            ( 1155.0*vn(ni,m) - 2772.0*vn(ni-1,m) +  2970.0*vn(ni-2,m) - 1540.0*vn(ni-3,m) +  315.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            ve(ni  ,m) = &
            (  315.0*vn(ni,m) -  420.0*vn(ni-1,m) +   378.0*vn(ni-2,m) -  180.0*vn(ni-3,m) +   35.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            ve(ni-1,m) = &
            (   35.0*vn(ni,m) +  140.0*vn(ni-1,m) -    70.0*vn(ni-2,m) +   28.0*vn(ni-3,m) -    5.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            ve(ni-2,m) = &
            (       -vn(ni,m) +    9.0*vn(ni-1,m) +     9.0*vn(ni-2,m) -        vn(ni-3,m) ) / 16.0                        ! fourth-order
!!            ve(ni+2,m) = (  63.0*vn(ni,m) -  90.0*vn(ni-1,m) +  35.0*vn(ni-2,m) ) / 8.0                                          ! third-order
!!            ve(ni+1,m) = (  35.0*vn(ni,m) -  42.0*vn(ni-1,m) +  15.0*vn(ni-2,m) ) / 8.0                                          ! third-order
!!            ve(ni  ,m) = (  15.0*vn(ni,m) -  10.0*vn(ni-1,m) +   3.0*vn(ni-2,m) ) / 8.0                                          ! third-order
!!            ve(ni-1,m) = (   3.0*vn(ni,m) +   6.0*vn(ni-1,m) -       vn(ni-2,m) ) / 8.0                                          ! third-order
!!            ve(ni-2,m) = (      -vn(ni,m) +   9.0*vn(ni-1,m) +   9.0*vn(ni-2,m) - vn(ni-3,m) ) / 16.0                            ! fourth-order
        end do    
    
    else
        ed = ni - 3

        do m=nst,ned
            !!ve(ni-2,m) = (      -vn(ni,m) +   9.0*vn(ni-1,m) +   9.0*vn(ni-2,m) -       vn(ni-3,m) ) / 16.0                      ! fourth-order
            !!ve(ni-1,m) = (   5.0*vn(ni,m) +  15.0*vn(ni-1,m) -   5.0*vn(ni-2,m) +       vn(ni-3,m) ) / 16.0                      ! fourth-order
            !!ve(ni  ,m) = (  35.0*vn(ni,m) -  35.0*vn(ni-1,m) +  21.0*vn(ni-2,m) -   5.0*vn(ni-3,m) ) / 16.0                      ! fourth-order
            ve(ni  ,m) = (  15.0*vn(ni,m) -  10.0*vn(ni-1,m) +   3.0*vn(ni-2,m) ) / 8.0                                          ! third-order
            ve(ni-2,m) = (  -5.0*vn(ni,m) +  60.0*vn(ni-1,m) +  90.0*vn(ni-2,m) -  20.0*vn(ni-3,m) +  3.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            ve(ni-1,m) = (  35.0*vn(ni,m) + 140.0*vn(ni-1,m) -  70.0*vn(ni-2,m) +  28.0*vn(ni-3,m) -  5.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            !!ve(ni  ,m) = ( 315.0*vn(ni,m) - 420.0*vn(ni-1,m) + 378.0*vn(ni-2,m) - 180.0*vn(ni-3,m) + 35.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            !!ve(ni-1,m) = (  3.0*vn(ni,m) +   6.0*vn(ni-1,m) -       vn(ni-2,m) ) / 8.0                                          ! third-order
            !!ve(ni-2,m) = (     -vn(ni,m) +   9.0*vn(ni-1,m) +   9.0*vn(ni-2,m) -       vn(ni-3,m) ) / 16.0                      ! fourth-order
        end do
    end if

    do m=nst,ned
        do i=st,ed
            ve(i,m) = ( 150.0*(vn(i+1,m) + vn(i  ,m)) - &
                         25.0*(vn(i+2,m) + vn(i-1,m)) + &
                          3.0*(vn(i+3,m) + vn(i-2,m)) )/256.0
        end do
    end do
    
end subroutine ve_via_node6

subroutine ve_via_node6cc(ni,nst,ned,ngn,vn,nfs,nfe,nge,ve)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_inter_scheme,nbc_intbc_scheme
    implicit none
    integer(kind_int),intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  ,intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)  :: nfs,nfe
    integer(kind_int),intent(in)  :: nge
    real(kind_real)  ,intent(out) :: ve(1-ngn:ni+ngn,nst:ned)
    integer(kind_int) :: i,m,st,ed
    
    if (ni < 4) then
        call ve_via_node2cc(ni,nst,ned,ngn,vn,nfs,nfe,nge,ve)
        return
    end if

    if (nfs == nbc_inter_scheme) then
        st = 1 - nge
    
    else if (nfs == nbc_intbc_scheme) then
    
        st = 4
        
        do m=nst,ned
        
            ve(-1,m) = (3003.0*vn(1,m) - 8580.0*vn(2,m) + 10010.0*vn(3,m) - 5460.0*vn(4,m) + 1155.0*vn(5,m) ) / 128.0   ! fifth-order
            ve(0,m) = (1155.0*vn(1,m) - 2772.0*vn(2,m) +  2970.0*vn(3,m) - 1540.0*vn(4,m) +  315.0*vn(5,m) ) / 128.0   ! fifth-order
            ve(1 ,m) = ( 315.0*vn(1,m) -  420.0*vn(2,m) +   378.0*vn(3,m) -  180.0*vn(4,m) +   35.0*vn(5,m) ) / 128.0   ! fifth-order
            ve(2 ,m) = (  35.0*vn(1,m) +  140.0*vn(2,m) -    70.0*vn(3,m) +   28.0*vn(4,m) -    5.0*vn(5,m) ) / 128.0   ! fifth-order
            ve(3 ,m) = (      -vn(1,m) +    9.0*vn(2,m) +     9.0*vn(3,m) -        vn(4,m) ) / 16.0                     ! fourth-order
!!            ve(-1,m) = ( 63.0*vn(1,m) - 90.0*vn(2,m) +  35.0*vn(3,m) ) / 8.0                                       ! third-order
!!            ve(0,m) = ( 35.0*vn(1,m) - 42.0*vn(2,m) +  15.0*vn(3,m) ) / 8.0                                       ! third-order
!!            ve(1 ,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +   3.0*vn(3,m) ) / 8.0                                       ! third-order
!!            ve(2 ,m) = (  3.0*vn(1,m) +  6.0*vn(2,m) -       vn(3,m) ) / 8.0                                       ! third-order
!!            ve(3 ,m) = (     -vn(1,m) +  9.0*vn(2,m) +   9.0*vn(3,m) - vn(4,m) ) / 16.0                            ! fourth-order
            
        end do
    
    else
        st = 4

        do m=nst,ned
            !!ve(3,m) = (      -vn(1,m) +   9.0*vn(2,m) +   9.0*vn(3,m) -       vn(4,m) ) / 16.0                   ! fourth-order
            !!ve(2,m) = (   5.0*vn(1,m) +  15.0*vn(2,m) -   5.0*vn(3,m) +       vn(4,m) ) / 16.0                   ! fourth-order
            !!ve(1,m) = (  35.0*vn(1,m) -  35.0*vn(2,m) +  21.0*vn(3,m) -   5.0*vn(4,m) ) / 16.0                   ! fourth-order
            ve(1,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                                       ! third-order
            ve(3,m) = (  -5.0*vn(1,m) +  60.0*vn(2,m) +  90.0*vn(3,m) -  20.0*vn(4,m) +  3.0*vn(5,m) ) / 128.0   ! fifth-order
            ve(2,m) = (  35.0*vn(1,m) + 140.0*vn(2,m) -  70.0*vn(3,m) +  28.0*vn(4,m) -  5.0*vn(5,m) ) / 128.0   ! fifth-order
            !!ve(1,m) = ( 315.0*vn(1,m) - 420.0*vn(2,m) + 378.0*vn(3,m) - 180.0*vn(4,m) + 35.0*vn(5,m) ) / 128.0   ! fifth-order
            !!ve(2,m) = ( 3.0*vn(1,m) + 6.0*vn(2,m) -     vn(3,m) ) / 8.0                                          ! third-order
            !!ve(3,m) = (    -vn(1,m) + 9.0*vn(2,m) + 9.0*vn(3,m) - vn(4,m) ) / 16.0                               ! fourth-order
        end do
    end if
    
    if (nfe == nbc_inter_scheme) then
        ed = ni + 1 + nge
    
    else if (nfe == nbc_intbc_scheme) then
    
        ed = ni - 2

        do m=nst,ned
            ve(ni+3,m) = &
            ( 3003.0*vn(ni,m) - 8580.0*vn(ni-1,m) + 10010.0*vn(ni-2,m) - 5460.0*vn(ni-3,m) + 1155.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            ve(ni+2,m) = &
            ( 1155.0*vn(ni,m) - 2772.0*vn(ni-1,m) +  2970.0*vn(ni-2,m) - 1540.0*vn(ni-3,m) +  315.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            ve(ni+1,m) = &
            (  315.0*vn(ni,m) -  420.0*vn(ni-1,m) +   378.0*vn(ni-2,m) -  180.0*vn(ni-3,m) +   35.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            ve(ni  ,m) = &
            (   35.0*vn(ni,m) +  140.0*vn(ni-1,m) -    70.0*vn(ni-2,m) +   28.0*vn(ni-3,m) -    5.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            ve(ni-1,m) = &
            (       -vn(ni,m) +    9.0*vn(ni-1,m) +     9.0*vn(ni-2,m) -        vn(ni-3,m) ) / 16.0                        ! fourth-order
!!            ve(ni+3,m) = (  63.0*vn(ni,m) -  90.0*vn(ni-1,m) +  35.0*vn(ni-2,m) ) / 8.0                                          ! third-order
!!            ve(ni+2,m) = (  35.0*vn(ni,m) -  42.0*vn(ni-1,m) +  15.0*vn(ni-2,m) ) / 8.0                                          ! third-order
!!            ve(ni+1,m) = (  15.0*vn(ni,m) -  10.0*vn(ni-1,m) +   3.0*vn(ni-2,m) ) / 8.0                                          ! third-order
!!            ve(ni  ,m) = (   3.0*vn(ni,m) +   6.0*vn(ni-1,m) -       vn(ni-2,m) ) / 8.0                                          ! third-order
!!            ve(ni-1,m) = (      -vn(ni,m) +   9.0*vn(ni-1,m) +   9.0*vn(ni-2,m) - vn(ni-3,m) ) / 16.0                            ! fourth-order
        end do    
    
    else
        ed = ni - 2

        do m=nst,ned
            !!ve(ni-1,m) = (      -vn(ni,m) +   9.0*vn(ni-1,m) +   9.0*vn(ni-2,m) -       vn(ni-3,m) ) / 16.0                      ! fourth-order
            !!ve(ni  ,m) = (   5.0*vn(ni,m) +  15.0*vn(ni-1,m) -   5.0*vn(ni-2,m) +       vn(ni-3,m) ) / 16.0                      ! fourth-order
            !!ve(ni+1,m) = (  35.0*vn(ni,m) -  35.0*vn(ni-1,m) +  21.0*vn(ni-2,m) -   5.0*vn(ni-3,m) ) / 16.0                      ! fourth-order
            ve(ni+1,m) = (  15.0*vn(ni,m) -  10.0*vn(ni-1,m) +   3.0*vn(ni-2,m) ) / 8.0                                          ! third-order
            ve(ni-1,m) = (  -5.0*vn(ni,m) +  60.0*vn(ni-1,m) +  90.0*vn(ni-2,m) -  20.0*vn(ni-3,m) +  3.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            ve(ni  ,m) = (  35.0*vn(ni,m) + 140.0*vn(ni-1,m) -  70.0*vn(ni-2,m) +  28.0*vn(ni-3,m) -  5.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            !!ve(ni+1,m) = ( 315.0*vn(ni,m) - 420.0*vn(ni-1,m) + 378.0*vn(ni-2,m) - 180.0*vn(ni-3,m) + 35.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            !!ve(ni  ,m) = (  3.0*vn(ni,m) +   6.0*vn(ni-1,m) -       vn(ni-2,m) ) / 8.0                                          ! third-order
            !!ve(ni-1,m) = (     -vn(ni,m) +   9.0*vn(ni-1,m) +   9.0*vn(ni-2,m) -       vn(ni-3,m) ) / 16.0                      ! fourth-order
        end do
    end if

    do m=nst,ned
        do i=st-1,ed-1
            ve(i+1,m) = ( 150.0*(vn(i+1,m) + vn(i  ,m)) - &
                           25.0*(vn(i+2,m) + vn(i-1,m)) + &
                            3.0*(vn(i+3,m) + vn(i-2,m)) )/256.0
        end do
    end do
    
end subroutine ve_via_node6cc

subroutine ve_via_node6e(ni,nst,ned,ngn,vn,nfs,nfe,nge,ve)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_inter_scheme
    implicit none
    integer(kind_int),intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  ,intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)  :: nfs,nfe
    integer(kind_int),intent(in)  :: nge
    real(kind_real)  ,intent(out) :: ve(1-ngn:ni+ngn,nst:ned)
    integer(kind_int) :: i,m,st,ed
    
    if (ni < 4) then
        call ve_via_node2(ni,nst,ned,ngn,vn,nfs,nfe,nge,ve)
        return
    end if

    if (nfs == nbc_inter_scheme) then
        
        st = -nge
    
    else
        st = 3

        do m=nst,ned
            !!ve(2,m) = (      -vn(1,m) +   9.0*vn(2,m) +   9.0*vn(3,m) -       vn(4,m) ) / 16.0                   ! fourth-order
            !!ve(1,m) = (   5.0*vn(1,m) +  15.0*vn(2,m) -   5.0*vn(3,m) +       vn(4,m) ) / 16.0                   ! fourth-order
            !!ve(0,m) = (  35.0*vn(1,m) -  35.0*vn(2,m) +  21.0*vn(3,m) -   5.0*vn(4,m) ) / 16.0                   ! fourth-order
            ve(0,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                                       ! third-order
            ve(2,m) = (  -5.0*vn(1,m) +  60.0*vn(2,m) +  90.0*vn(3,m) -  20.0*vn(4,m) +  3.0*vn(5,m) ) / 128.0   ! fifth-order
            ve(1,m) = (  35.0*vn(1,m) + 140.0*vn(2,m) -  70.0*vn(3,m) +  28.0*vn(4,m) -  5.0*vn(5,m) ) / 128.0   ! fifth-order
            !!ve(0,m) = ( 315.0*vn(1,m) - 420.0*vn(2,m) + 378.0*vn(3,m) - 180.0*vn(4,m) + 35.0*vn(5,m) ) / 128.0   ! fifth-order
            !!ve(1,m) = ( 3.0*vn(1,m) + 6.0*vn(2,m) -     vn(3,m) ) / 8.0                                          ! third-order
            !!ve(2,m) = (    -vn(1,m) + 9.0*vn(2,m) + 9.0*vn(3,m) - vn(4,m) ) / 16.0                               ! fourth-order
        end do
    end if
    
    if (nfe == nbc_inter_scheme) then

        ed = ni + nge  
    
    else
        ed = ni - 3

        do m=nst,ned
            !!ve(ni-2,m) = (      -vn(ni,m) +   9.0*vn(ni-1,m) +   9.0*vn(ni-2,m) -       vn(ni-3,m) ) / 16.0                      ! fourth-order
            !!ve(ni-1,m) = (   5.0*vn(ni,m) +  15.0*vn(ni-1,m) -   5.0*vn(ni-2,m) +       vn(ni-3,m) ) / 16.0                      ! fourth-order
            !!ve(ni  ,m) = (  35.0*vn(ni,m) -  35.0*vn(ni-1,m) +  21.0*vn(ni-2,m) -   5.0*vn(ni-3,m) ) / 16.0                      ! fourth-order
            ve(ni  ,m) = (  15.0*vn(ni,m) -  10.0*vn(ni-1,m) +   3.0*vn(ni-2,m) ) / 8.0                                          ! third-order
            ve(ni-2,m) = (  -5.0*vn(ni,m) +  60.0*vn(ni-1,m) +  90.0*vn(ni-2,m) -  20.0*vn(ni-3,m) +  3.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            ve(ni-1,m) = (  35.0*vn(ni,m) + 140.0*vn(ni-1,m) -  70.0*vn(ni-2,m) +  28.0*vn(ni-3,m) -  5.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            !!ve(ni  ,m) = ( 315.0*vn(ni,m) - 420.0*vn(ni-1,m) + 378.0*vn(ni-2,m) - 180.0*vn(ni-3,m) + 35.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            !!ve(ni-1,m) = (  3.0*vn(ni,m) +   6.0*vn(ni-1,m) -       vn(ni-2,m) ) / 8.0                                          ! third-order
            !!ve(ni-2,m) = (     -vn(ni,m) +   9.0*vn(ni-1,m) +   9.0*vn(ni-2,m) -       vn(ni-3,m) ) / 16.0                      ! fourth-order
        end do
    end if

    do m=nst,ned
        do i=st,ed
            ve(i,m) = ( 150.0*(vn(i+1,m) + vn(i  ,m)) - &
                         25.0*(vn(i+2,m) + vn(i-1,m)) + &
                          3.0*(vn(i+3,m) + vn(i-2,m)) )/256.0
        end do
    end do
    
end subroutine ve_via_node6e

subroutine ve_via_node6ecc(ni,nst,ned,ngn,vn,nfs,nfe,nge,ve)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_inter_scheme
    implicit none
    integer(kind_int),intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  ,intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)  :: nfs,nfe
    integer(kind_int),intent(in)  :: nge
    real(kind_real)  ,intent(out) :: ve(1-ngn:ni+ngn,nst:ned)
    integer(kind_int) :: i,m,st,ed
    
    if (ni < 4) then
        call ve_via_node2cc(ni,nst,ned,ngn,vn,nfs,nfe,nge,ve)
        return
    end if

    if (nfs == nbc_inter_scheme) then
        
        st = 1 - nge
    
    else
        st = 4

        do m=nst,ned
            !!ve(3,m) = (      -vn(1,m) +   9.0*vn(2,m) +   9.0*vn(3,m) -       vn(4,m) ) / 16.0                   ! fourth-order
            !!ve(2,m) = (   5.0*vn(1,m) +  15.0*vn(2,m) -   5.0*vn(3,m) +       vn(4,m) ) / 16.0                   ! fourth-order
            !!ve(1,m) = (  35.0*vn(1,m) -  35.0*vn(2,m) +  21.0*vn(3,m) -   5.0*vn(4,m) ) / 16.0                   ! fourth-order
            ve(1,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                                       ! third-order
            ve(3,m) = (  -5.0*vn(1,m) +  60.0*vn(2,m) +  90.0*vn(3,m) -  20.0*vn(4,m) +  3.0*vn(5,m) ) / 128.0   ! fifth-order
            ve(2,m) = (  35.0*vn(1,m) + 140.0*vn(2,m) -  70.0*vn(3,m) +  28.0*vn(4,m) -  5.0*vn(5,m) ) / 128.0   ! fifth-order
            !!ve(1,m) = ( 315.0*vn(1,m) - 420.0*vn(2,m) + 378.0*vn(3,m) - 180.0*vn(4,m) + 35.0*vn(5,m) ) / 128.0   ! fifth-order
            !!ve(2,m) = ( 3.0*vn(1,m) + 6.0*vn(2,m) -     vn(3,m) ) / 8.0                                          ! third-order
            !!ve(3,m) = (    -vn(1,m) + 9.0*vn(2,m) + 9.0*vn(3,m) - vn(4,m) ) / 16.0                               ! fourth-order
        end do
    end if
    
    if (nfe == nbc_inter_scheme) then

        ed = ni + 1 + nge  
    
    else
        ed = ni - 2

        do m=nst,ned
            !!ve(ni-1,m) = (      -vn(ni,m) +   9.0*vn(ni-1,m) +   9.0*vn(ni-2,m) -       vn(ni-3,m) ) / 16.0                      ! fourth-order
            !!ve(ni  ,m) = (   5.0*vn(ni,m) +  15.0*vn(ni-1,m) -   5.0*vn(ni-2,m) +       vn(ni-3,m) ) / 16.0                      ! fourth-order
            !!ve(ni+1,m) = (  35.0*vn(ni,m) -  35.0*vn(ni-1,m) +  21.0*vn(ni-2,m) -   5.0*vn(ni-3,m) ) / 16.0                      ! fourth-order
            ve(ni+1,m) = (  15.0*vn(ni,m) -  10.0*vn(ni-1,m) +   3.0*vn(ni-2,m) ) / 8.0                                          ! third-order
            ve(ni-1,m) = (  -5.0*vn(ni,m) +  60.0*vn(ni-1,m) +  90.0*vn(ni-2,m) -  20.0*vn(ni-3,m) +  3.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            ve(ni  ,m) = (  35.0*vn(ni,m) + 140.0*vn(ni-1,m) -  70.0*vn(ni-2,m) +  28.0*vn(ni-3,m) -  5.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            !!ve(ni+1,m) = ( 315.0*vn(ni,m) - 420.0*vn(ni-1,m) + 378.0*vn(ni-2,m) - 180.0*vn(ni-3,m) + 35.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            !!ve(ni  ,m) = (  3.0*vn(ni,m) +   6.0*vn(ni-1,m) -       vn(ni-2,m) ) / 8.0                                          ! third-order
            !!ve(ni-1,m) = (     -vn(ni,m) +   9.0*vn(ni-1,m) +   9.0*vn(ni-2,m) -       vn(ni-3,m) ) / 16.0                      ! fourth-order
        end do
    end if

    do m=nst,ned
        do i=st-1,ed-1
            ve(i+1,m) = ( 150.0*(vn(i+1,m) + vn(i  ,m)) - &
                           25.0*(vn(i+2,m) + vn(i-1,m)) + &
                            3.0*(vn(i+3,m) + vn(i-2,m)) )/256.0
        end do
    end do
    
end subroutine ve_via_node6ecc

subroutine ve_via_node6w(ni,nst,ned,ngn,vn,nfs,nfe,nge,ve)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_inter_scheme
    implicit none
    integer(kind_int),intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  ,intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)  :: nfs,nfe
    integer(kind_int),intent(in)  :: nge
    real(kind_real)  ,intent(out) :: ve(1-ngn:ni+ngn,nst:ned)
    real(kind_real), parameter :: ck(3)=(/3.D0/16.D0,10.D0/16.D0,3.D0/16.D0/)
    real(kind_real), parameter :: eps=1.0e-6
    real(kind_real)            :: f(3),s(3),t(3)
    real(kind_real)            :: is,bk(3),bksum,wk(3),rvk(3)
    integer(kind_int)          :: i,m,n,st,ed
    
    if (nfs == nbc_inter_scheme) then
        st = -nge
    else
        st = 3

        do m=nst,ned
            ! ve(2,m) = (     -vn(1,m) +  9.0*vn(2,m) +  9.0*vn(3,m) -      vn(4,m) ) / 16.0                       ! fourth-order
            ! ve(1,m) = (   5.0*vn(1,m) +  15.0*vn(2,m) -   5.0*vn(3,m) +       vn(4,m) ) / 16.0                   ! fourth-order
            ! ve(0,m) = (  35.0*vn(1,m) -  35.0*vn(2,m) +  21.0*vn(3,m) -   5.0*vn(4,m) ) / 16.0                   ! fourth-order
            ! ve(0,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                                       ! third-order
            ve(2,m) = (  -5.0*vn(1,m) +  60.0*vn(2,m) +  90.0*vn(3,m) -  20.0*vn(4,m) +  3.0*vn(5,m) ) / 128.0   ! fifth-order
            ve(1,m) = (  35.0*vn(1,m) + 140.0*vn(2,m) -  70.0*vn(3,m) +  28.0*vn(4,m) -  5.0*vn(5,m) ) / 128.0   ! fifth-order
            ve(0,m) = ( 315.0*vn(1,m) - 420.0*vn(2,m) + 378.0*vn(3,m) - 180.0*vn(4,m) + 35.0*vn(5,m) ) / 128.0   ! fifth-order
        end do
    end if
    
    if (nfe == nbc_inter_scheme) then
        ed = ni + nge
    else
        ed = ni - 3

        do m=nst,ned
            ! ve(ni-2,m) = (      -vn(ni,m) +   9.0*vn(ni-1,m) +   9.0*vn(ni-2,m) -       vn(ni-3,m) ) / 16.0                      ! fourth-order
            ! ve(ni-1,m) = (   5.0*vn(ni,m) +  15.0*vn(ni-1,m) -   5.0*vn(ni-2,m) +       vn(ni-3,m) ) / 16.0                      ! fourth-order
            ! ve(ni,m)   = (  35.0*vn(ni,m) -  35.0*vn(ni-1,m) +  21.0*vn(ni-2,m) -   5.0*vn(ni-3,m) ) / 16.0                      ! fourth-order
            ! ve(ni,m)   = (  15.0*vn(ni,m) -  10.0*vn(ni-1,m) +   3.0*vn(ni-2,m) ) / 8.0                                          ! third-order
            ve(ni-2,m) = (  -5.0*vn(ni,m) +  60.0*vn(ni-1,m) +  90.0*vn(ni-2,m) -  20.0*vn(ni-3,m) +  3.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            ve(ni-1,m) = (  35.0*vn(ni,m) + 140.0*vn(ni-1,m) -  70.0*vn(ni-2,m) +  28.0*vn(ni-3,m) -  5.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            ve(ni,m)   = ( 315.0*vn(ni,m) - 420.0*vn(ni-1,m) + 378.0*vn(ni-2,m) - 180.0*vn(ni-3,m) + 35.0*vn(ni-4,m) ) / 128.0   ! fifth-order
        end do
    end if

    do m=nst,ned
        do i=st,ed
            f(1) = (       vn(i-2,m) -  3.0*vn(i-1,m) - 21.0*vn(i  ,m) + 23.0*vn(i+1,m) ) / 24.0        ! third-order
            f(2) = (       vn(i-1,m) - 27.0*vn(i  ,m) + 27.0*vn(i+1,m) -      vn(i+2,m) ) / 24.0
            f(3) = ( -23.0*vn(i  ,m) + 21.0*vn(i+1,m) +  3.0*vn(i+2,m) -      vn(i+3,m) ) / 24.0
       
            s(1) = (    -vn(i-2,m) + 5.0*vn(i-1,m) - 7.0*vn(i  ,m) + 3.0*vn(i+1,m) ) / 2.0              ! second-order
            s(2) = (     vn(i-1,m) -     vn(i  ,m) -     vn(i+1,m) +     vn(i+2,m) ) / 2.0
            s(3) = ( 3.0*vn(i  ,m) - 7.0*vn(i+1,m) + 5.0*vn(i+2,m) -     vn(i+3,m) ) / 2.0
       
            t(1) = -vn(i-2,m) + 3.0*vn(i-1,m) - 3.0*vn(i  ,m) + vn(i+1,m)                               ! first-order
            t(2) = -vn(i-1,m) + 3.0*vn(i  ,m) - 3.0*vn(i+1,m) + vn(i+2,m)
            t(3) = -vn(i  ,m) + 3.0*vn(i+1,m) - 3.0*vn(i+2,m) + vn(i+3,m)
            do n=1,3
                is = f(n)*f(n) + t(n)*t(n)
                is = is + eps
                bk(n) = ck(n)/(is*is)
            end do
    
            bksum = bk(1) + bk(2) + bk(3)
            do n=1,3
                !!wk(n) = ck(n)
                wk(n) = bk(n)/bksum
            end do
            rvk(1) = (     vn(i-2,m) -  5.0*vn(i-1,m) + 15.0*vn(i  ,m) + 5.0*vn(i+1,m) ) / 16.0
            rvk(2) = (    -vn(i-1,m) +  9.0*vn(i  ,m) +  9.0*vn(i+1,m) -     vn(i+2,m) ) / 16.0
            rvk(3) = ( 5.0*vn(i  ,m) + 15.0*vn(i+1,m) -  5.0*vn(i+2,m) +     vn(i+3,m) ) / 16.0
       
            ve(i,m) = wk(1)*rvk(1) + wk(2)*rvk(2) + wk(3)*rvk(3)                                       ! sixth-order
        end do
    end do

end subroutine ve_via_node6w

subroutine ve_via_node6wcc(ni,nst,ned,ngn,vn,nfs,nfe,nge,ve)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_inter_scheme
    implicit none
    integer(kind_int),intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  ,intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)  :: nfs,nfe
    integer(kind_int),intent(in)  :: nge
    real(kind_real)  ,intent(out) :: ve(1-ngn:ni+ngn,nst:ned)
    real(kind_real), parameter :: ck(3)=(/3.D0/16.D0,10.D0/16.D0,3.D0/16.D0/)
    real(kind_real), parameter :: eps=1.0e-6
    real(kind_real)            :: f(3),s(3),t(3)
    real(kind_real)            :: is,bk(3),bksum,wk(3),rvk(3)
    integer(kind_int)          :: i,m,n,st,ed
    
    if (nfs == nbc_inter_scheme) then
        st = 1 - nge
    else
        st = 4

        do m=nst,ned
            ! ve(3,m) = (     -vn(1,m) +  9.0*vn(2,m) +  9.0*vn(3,m) -      vn(4,m) ) / 16.0                       ! fourth-order
            ! ve(2,m) = (   5.0*vn(1,m) +  15.0*vn(2,m) -   5.0*vn(3,m) +       vn(4,m) ) / 16.0                   ! fourth-order
            ! ve(1,m) = (  35.0*vn(1,m) -  35.0*vn(2,m) +  21.0*vn(3,m) -   5.0*vn(4,m) ) / 16.0                   ! fourth-order
            ! ve(1,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                                       ! third-order
            ve(3,m) = (  -5.0*vn(1,m) +  60.0*vn(2,m) +  90.0*vn(3,m) -  20.0*vn(4,m) +  3.0*vn(5,m) ) / 128.0   ! fifth-order
            ve(2,m) = (  35.0*vn(1,m) + 140.0*vn(2,m) -  70.0*vn(3,m) +  28.0*vn(4,m) -  5.0*vn(5,m) ) / 128.0   ! fifth-order
            ve(1,m) = ( 315.0*vn(1,m) - 420.0*vn(2,m) + 378.0*vn(3,m) - 180.0*vn(4,m) + 35.0*vn(5,m) ) / 128.0   ! fifth-order
        end do
    end if
    
    if (nfe == nbc_inter_scheme) then
        ed = ni + 1 + nge
    else
        ed = ni - 2

        do m=nst,ned
            ! ve(ni-1,m) = (      -vn(ni,m) +   9.0*vn(ni-1,m) +   9.0*vn(ni-2,m) -       vn(ni-3,m) ) / 16.0                      ! fourth-order
            ! ve(ni  ,m) = (   5.0*vn(ni,m) +  15.0*vn(ni-1,m) -   5.0*vn(ni-2,m) +       vn(ni-3,m) ) / 16.0                      ! fourth-order
            ! ve(ni+1,m)   = (  35.0*vn(ni,m) -  35.0*vn(ni-1,m) +  21.0*vn(ni-2,m) -   5.0*vn(ni-3,m) ) / 16.0                      ! fourth-order
            ! ve(ni+1,m)   = (  15.0*vn(ni,m) -  10.0*vn(ni-1,m) +   3.0*vn(ni-2,m) ) / 8.0                                          ! third-order
            ve(ni-1,m) = (  -5.0*vn(ni,m) +  60.0*vn(ni-1,m) +  90.0*vn(ni-2,m) -  20.0*vn(ni-3,m) +  3.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            ve(ni  ,m) = (  35.0*vn(ni,m) + 140.0*vn(ni-1,m) -  70.0*vn(ni-2,m) +  28.0*vn(ni-3,m) -  5.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            ve(ni+1,m)   = ( 315.0*vn(ni,m) - 420.0*vn(ni-1,m) + 378.0*vn(ni-2,m) - 180.0*vn(ni-3,m) + 35.0*vn(ni-4,m) ) / 128.0   ! fifth-order
        end do
    end if

    do m=nst,ned
        do i=st-1,ed-1
            f(1) = (       vn(i-2,m) -  3.0*vn(i-1,m) - 21.0*vn(i  ,m) + 23.0*vn(i+1,m) ) / 24.0        ! third-order
            f(2) = (       vn(i-1,m) - 27.0*vn(i  ,m) + 27.0*vn(i+1,m) -      vn(i+2,m) ) / 24.0
            f(3) = ( -23.0*vn(i  ,m) + 21.0*vn(i+1,m) +  3.0*vn(i+2,m) -      vn(i+3,m) ) / 24.0
       
            s(1) = (    -vn(i-2,m) + 5.0*vn(i-1,m) - 7.0*vn(i  ,m) + 3.0*vn(i+1,m) ) / 2.0              ! second-order
            s(2) = (     vn(i-1,m) -     vn(i  ,m) -     vn(i+1,m) +     vn(i+2,m) ) / 2.0
            s(3) = ( 3.0*vn(i  ,m) - 7.0*vn(i+1,m) + 5.0*vn(i+2,m) -     vn(i+3,m) ) / 2.0
       
            t(1) = -vn(i-2,m) + 3.0*vn(i-1,m) - 3.0*vn(i  ,m) + vn(i+1,m)                               ! first-order
            t(2) = -vn(i-1,m) + 3.0*vn(i  ,m) - 3.0*vn(i+1,m) + vn(i+2,m)
            t(3) = -vn(i  ,m) + 3.0*vn(i+1,m) - 3.0*vn(i+2,m) + vn(i+3,m)
            do n=1,3
                is = f(n)*f(n) + t(n)*t(n)
                is = is + eps
                bk(n) = ck(n)/(is*is)
            end do
    
            bksum = bk(1) + bk(2) + bk(3)
            do n=1,3
                !!wk(n) = ck(n)
                wk(n) = bk(n)/bksum
            end do
            rvk(1) = (     vn(i-2,m) -  5.0*vn(i-1,m) + 15.0*vn(i  ,m) + 5.0*vn(i+1,m) ) / 16.0
            rvk(2) = (    -vn(i-1,m) +  9.0*vn(i  ,m) +  9.0*vn(i+1,m) -     vn(i+2,m) ) / 16.0
            rvk(3) = ( 5.0*vn(i  ,m) + 15.0*vn(i+1,m) -  5.0*vn(i+2,m) +     vn(i+3,m) ) / 16.0
       
            ve(i+1,m) = wk(1)*rvk(1) + wk(2)*rvk(2) + wk(3)*rvk(3)                                       ! sixth-order
        end do
    end do

end subroutine ve_via_node6wcc

subroutine ve_via_node8i(ni,nst,ned,ngn,vn,nfs,nfe,nge,ve)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_intbc_scheme,nbc_bound_scheme
    implicit none
    integer(kind_int),intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  ,intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)  :: nfs,nfe
    integer(kind_int),intent(in)  :: nge
    real(kind_real)  ,intent(out) :: ve(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: ma(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵ����Aϵ��
    real(kind_real)   :: mb(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵ����Bϵ��
    real(kind_real)   :: mc(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵ����Cϵ��
    real(kind_real)   :: md(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵ��Ҳ�Dϵ��
    real(kind_real)   :: ms(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵļ�����     
    integer(kind_int) :: i,m,n,st,ed,stc,edc

    if (ni < 4) then
        call ve_via_node2(ni,nst,ned,ngn,vn,nfs,nfe,nge,ve)
        return
    end if
    
    
    !>
    !!�����߽�HDCS��߽�������Ҳ��ֵ
    if (nfs == nbc_bound_scheme) then
        st = 3       !< ����nfs=1ʱ���ڵ��ֵ��ʽ�������ʼ����
        stc = 1      !< ����nfs=1ʱ����ֵ�������ʼ����
        ma(st-4)=0.0
        ma(st-3)=0.0
        ma(st-2)=0.0                           !< ���Խ�׷�ϵ����Aϵ������߽�����ֵ��
        ma(st-1)=0.0                           !< ���Խ�׷�ϵ����Aϵ�����ٽ���߽�����ֵ��
        mc(st-4)=0.0
        mc(st-3)=0.0
        mc(st-2)=0.0                           !< ���Խ�׷�ϵ����Cϵ������߽�����ֵ��
        mc(st-1)=0.0                           !< ���Խ�׷�ϵ����Cϵ�����ٽ���߽�����ֵ��

        do m=nst,ned
            ve(1,m) = ( 3.0*vn(1,m) +   6.0*vn(2,m) -      vn(3,m) ) / 8.0                        !< ��������ֵ��߽�     
            ve(2,m) = (    -vn(1,m) +   9.0*vn(2,m) +  9.0*vn(3,m) -      vn(4,m) )/16.0          !< �Ľ�����ֵ+���׺�ɢ�����ٽ��߽磨i=2��
        end do                                                                                   
    else if (nfs == nbc_intbc_scheme) then
    
        st = 3
        stc = -2

        ma(st-5)=0.0
        ma(st-4)=0.0
        ma(st-3)=0.0
        ma(st-2)=0.0                           !< ���Խ�׷�ϵ����Aϵ������߽�����ֵ��
        ma(st-1)=0.0                           !< ���Խ�׷�ϵ����Aϵ�����ٽ���߽�����ֵ��
        mc(st-5)=0.0
        mc(st-4)=0.0
        mc(st-3)=0.0
        mc(st-2)=0.0                           !< ���Խ�׷�ϵ����Cϵ������߽�����ֵ��
        mc(st-1)=0.0                           !< ���Խ�׷�ϵ����Cϵ�����ٽ���߽�����ֵ��
                
        do m=nst,ned
        
            ve(-2,m) = (3003.0*vn(1,m) - 8580.0*vn(2,m) + 10010.0*vn(3,m) - 5460.0*vn(4,m) + 1155.0*vn(5,m) ) / 128.0   ! fifth-order
            ve(-1,m) = (1155.0*vn(1,m) - 2772.0*vn(2,m) +  2970.0*vn(3,m) - 1540.0*vn(4,m) +  315.0*vn(5,m) ) / 128.0   ! fifth-order
            ve(0 ,m) = ( 315.0*vn(1,m) -  420.0*vn(2,m) +   378.0*vn(3,m) -  180.0*vn(4,m) +   35.0*vn(5,m) ) / 128.0   ! fifth-order
            ve(1 ,m) = (  35.0*vn(1,m) +  140.0*vn(2,m) -    70.0*vn(3,m) +   28.0*vn(4,m) -    5.0*vn(5,m) ) / 128.0   ! fifth-order
            ve(2 ,m) = (  -5.0*vn(1,m) +   60.0*vn(2,m) +    90.0*vn(3,m) -   20.0*vn(4,m) +    3.0*vn(5,m) ) / 128.0
            !!ve(2 ,m) = (      -vn(1,m) +    9.0*vn(2,m) +     9.0*vn(3,m) -        vn(4,m) ) / 16.0                     ! fourth-order            
!!            ve(3 ,m) = (150.0*(vn(4,m) + vn(3,m)) - 25.0*(vn(5,m) + vn(2,m)) + 3.0*(vn(6,m) + vn(1,m)) )/256.0          ! sixth-order
!!            ve(-2,m) = ( 63.0*vn(1,m) - 90.0*vn(2,m) +  35.0*vn(3,m) ) / 8.0                                       ! third-order
!!            ve(-1,m) = ( 35.0*vn(1,m) - 42.0*vn(2,m) +  15.0*vn(3,m) ) / 8.0                                       ! third-order
!!            ve(0 ,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +   3.0*vn(3,m) ) / 8.0                                       ! third-order
!!            ve(1 ,m) = (  3.0*vn(1,m) +  6.0*vn(2,m) -       vn(3,m) ) / 8.0                                       ! third-order
!!            ve(2 ,m) = (     -vn(1,m) +  9.0*vn(2,m) +   9.0*vn(3,m) - vn(4,m) ) / 16.0                            ! fourth-order
!!            ve(3 ,m) = (150.0*(vn(4,m) + vn(3,m)) - 25.0*(vn(5,m) + vn(2,m)) + 3.0*(vn(6,m) + vn(1,m)) )/256.0     ! sixth-order
            
        end do
    else
        st = -3    !< ����nfs=0������ȫ�ڵ�߽��ʽ����-1ʱ���ڵ��ֵ��ʽ�������ʼ����
        stc = -4   !< ����nfs=0������ȫ�ڵ�߽��ʽ����-1ʱ����ֵ�������ʼ����
        ma(st-1)=0.0
        mc(st-1)=0.0
        do m=nst,ned
            ve(-4,m) = (  3.0*vn(-4,m) +  6.0*vn(-3,m) -      vn(-2,m) ) / 8.0   !< ��������ֵ��߽� @attention ��ֵ��δʹ�ã�ֻ��Ϊ���������Խ�׷��
        end do

    end if
    !>

    !>
    !!�����߽�HDCS�ұ߽�������Ҳ��ֵ
    if (nfe == nbc_bound_scheme) then
        ed = ni - 3      !< ����nfe=1ʱ���ڵ��ֵ��ʽ�������ֹ����
        edc = ni - 1     !< ����nfe=1ʱ����ֵ�������ֹ����

        ma(ed+1)=0.0                           !< ���Խ�׷�ϵ����Aϵ�����ٽ��ұ߽�����ֵ��
        ma(ed+2)=0.0                           !< ���Խ�׷�ϵ����Aϵ�����ұ߽�����ֵ��
        ma(ed+3)=0.0
        ma(ed+4)=0.0
        mc(ed+1)=0.0                           !< ���Խ�׷�ϵ����Cϵ�����ٽ��ұ߽�����ֵ��
        mc(ed+2)=0.0                           !< ���Խ�׷�ϵ����Cϵ�����ұ߽�����ֵ��
        mc(ed+3)=0.0
        mc(ed+4)=0.0

        do m=nst,ned
            ve(ni-1,m) = (  3.0*vn(ni,m) +   6.0*vn(ni-1,m) -      vn(ni-2,m) ) / 8.0   ! third-order         !< �����Ҳ��ֵ��߽�
            ve(ni-2,m) = (     -vn(ni,m) +   9.0*vn(ni-1,m) +  9.0*vn(ni-2,m) -      vn(ni-3,m) )/16.0                              
        end do 
    else if (nfe == nbc_intbc_scheme) then
        ed = ni - 3
        edc = ni + 2

        ma(ed+1)=0.0                           !< ���Խ�׷�ϵ����Aϵ�����ٽ��ұ߽�����ֵ��
        ma(ed+2)=0.0                           !< ���Խ�׷�ϵ����Aϵ�����ұ߽�����ֵ��
        ma(ed+3)=0.0
        ma(ed+4)=0.0
        ma(ed+5)=0.0
        mc(ed+1)=0.0                           !< ���Խ�׷�ϵ����Cϵ�����ٽ��ұ߽�����ֵ��
        mc(ed+2)=0.0                           !< ���Խ�׷�ϵ����Cϵ�����ұ߽�����ֵ��
        mc(ed+3)=0.0
        mc(ed+4)=0.0
        mc(ed+5)=0.0

        do m=nst,ned
            ve(ni+2,m) = &
            ( 3003.0*vn(ni,m) - 8580.0*vn(ni-1,m) + 10010.0*vn(ni-2,m) - 5460.0*vn(ni-3,m) + 1155.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            ve(ni+1,m) = &
            ( 1155.0*vn(ni,m) - 2772.0*vn(ni-1,m) +  2970.0*vn(ni-2,m) - 1540.0*vn(ni-3,m) +  315.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            ve(ni  ,m) = &
            (  315.0*vn(ni,m) -  420.0*vn(ni-1,m) +   378.0*vn(ni-2,m) -  180.0*vn(ni-3,m) +   35.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            ve(ni-1,m) = &
            (   35.0*vn(ni,m) +  140.0*vn(ni-1,m) -    70.0*vn(ni-2,m) +   28.0*vn(ni-3,m) -    5.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            ve(ni-2,m) = &
            (   -5.0*vn(ni,m) +   60.0*vn(ni-1,m) +    90.0*vn(ni-2,m) -   20.0*vn(ni-3,m) +    3.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            !!ve(ni-2,m) = (       -vn(ni,m) +    9.0*vn(ni-1,m) +     9.0*vn(ni-2,m) -        vn(ni-3,m) ) / 16.0                        ! fourth-order
!!            ve(ni-3,m) = ( 150.0*(vn(ni-2,m) + vn(ni-3,m)) - 25.0*(vn(ni-1,m) + vn(ni-4,m)) + 3.0*(vn(ni ,m) + vn(ni-5,m)) )/256.0      ! sixth-order
!!            ve(ni+2,m) = (  63.0*vn(ni,m) -  90.0*vn(ni-1,m) +  35.0*vn(ni-2,m) ) / 8.0                                            ! third-order
!!            ve(ni+1,m) = (  35.0*vn(ni,m) -  42.0*vn(ni-1,m) +  15.0*vn(ni-2,m) ) / 8.0                                            ! third-order
!!            ve(ni  ,m) = (  15.0*vn(ni,m) -  10.0*vn(ni-1,m) +   3.0*vn(ni-2,m) ) / 8.0                                            ! third-order
!!            ve(ni-1,m) = (   3.0*vn(ni,m) +   6.0*vn(ni-1,m) -       vn(ni-2,m) ) / 8.0                                            ! third-order
!!            ve(ni-2,m) = (      -vn(ni,m) +   9.0*vn(ni-1,m) +   9.0*vn(ni-2,m) - vn(ni-3,m) ) / 16.0                              ! fourth-order
!!            ve(ni-3,m) = ( 150.0*(vn(ni-2,m) + vn(ni-3,m)) - 25.0*(vn(ni-1,m) + vn(ni-4,m)) + 3.0*(vn(ni ,m) + vn(ni-5,m)) )/256.0 ! sixth-order
        end do                                                                                                       
    else
        ed = ni + 3    !< ����nfe=0������ȫ�ڵ�߽��ʽ����-1ʱ���ڵ��ֵ��ʽ�������ֹ����
        edc = ni + 4   !< ����nfe=0������ȫ�ڵ�߽��ʽ����-1ʱ����ֵ�������ֹ����
        ma(ed+1)=0.0
        mc(ed+1)=0.0

        do m=nst,ned
            ve(ni+4,m) = (  3.0*vn(ni+5,m) +  6.0*vn(ni+4,m) -      vn(ni+3,m) ) / 8.0   !< ��������ֵ�ұ߽� @attention ��ֵ��δʹ�ã�ֻ��Ϊ���������Խ�׷��
        end do
    end if
    !>

    !>
    !!�����߽�HDCS�����Ҳ��ڵ��ֵ
    do m=nst,ned
        do i=st,ed
            ve(i,m) =    ( 350.0 *( vn(i+1,m) + vn(i,m  ) )          +  &
                            35.0 *( vn(i+2,m) + vn(i-1,m) )          -  &
                                  ( vn(i+3,m) + vn(i-2,m) ) ) /448.0       
        end do
    end do
    
    do i=stc,edc
        mb(i)=1.0                          !< ���Խ�׷�ϵ����Bϵ��
    enddo
  
    do i=st,ed
        ma(i)=5.0/14.0                     !< ���Խ�׷�ϵ����Aϵ�����ڵ�����ֵ��
    end do

    do i=st,ed                                                
        mc(i)=5.0/14.0                     !< ���Խ�׷�ϵ����Cϵ�����ڵ�����ֵ��
    end do                                               


    
    do m=nst,ned
        do i=stc,edc
            md(i) = ve(i,m)               !< ���Խ�׷�ϵ��Ҳ�Dϵ��������ֵ��
        end do
        call tri_pursue(ma,mb,mc,md,ms,stc,edc,1-ngn,ni+ngn) !< ���Խ�׷��
        do i=stc,edc
            ve(i,m) = ms(i)               !< �������Խ�׷�Ͻ⣨����ֵ��
        end do
    end do
              
end subroutine ve_via_node8i

subroutine ve_via_node8(ni,nst,ned,ngn,vn,nfs,nfe,nge,ve)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_inter_scheme,nbc_intbc_scheme
    implicit none
    integer(kind_int),intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  ,intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)  :: nfs,nfe
    integer(kind_int),intent(in)  :: nge
    real(kind_real)  ,intent(out) :: ve(1-ngn:ni+ngn,nst:ned)
    integer(kind_int) :: i,m,st,ed
    
    if (ni < 4) then
        call ve_via_node2(ni,nst,ned,ngn,vn,nfs,nfe,nge,ve)
        return
    end if

    if (nfs == nbc_inter_scheme) then
        st = -nge
    
    else if (nfs == nbc_intbc_scheme) then
    
        st = 4
        
        do m=nst,ned
        
            ve(-2,m) = (3003.0*vn(1,m) - 8580.0*vn(2,m) + 10010.0*vn(3,m) - 5460.0*vn(4,m) + 1155.0*vn(5,m) ) / 128.0   ! fifth-order
            ve(-1,m) = (1155.0*vn(1,m) - 2772.0*vn(2,m) +  2970.0*vn(3,m) - 1540.0*vn(4,m) +  315.0*vn(5,m) ) / 128.0   ! fifth-order
            ve(0 ,m) = ( 315.0*vn(1,m) -  420.0*vn(2,m) +   378.0*vn(3,m) -  180.0*vn(4,m) +   35.0*vn(5,m) ) / 128.0   ! fifth-order
            ve(1 ,m) = (  35.0*vn(1,m) +  140.0*vn(2,m) -    70.0*vn(3,m) +   28.0*vn(4,m) -    5.0*vn(5,m) ) / 128.0   ! fifth-order
            ve(2 ,m) = (  -5.0*vn(1,m) +   60.0*vn(2,m) +    90.0*vn(3,m) -   20.0*vn(4,m) +    3.0*vn(5,m) ) / 128.0
            !!ve(2 ,m) = (      -vn(1,m) +    9.0*vn(2,m) +     9.0*vn(3,m) -        vn(4,m) ) / 16.0                     ! fourth-order            
            ve(3 ,m) = (150.0*(vn(4,m) + vn(3,m)) - 25.0*(vn(5,m) + vn(2,m)) + 3.0*(vn(6,m) + vn(1,m)) )/256.0          ! sixth-order
!!            ve(-2,m) = ( 63.0*vn(1,m) - 90.0*vn(2,m) +  35.0*vn(3,m) ) / 8.0                                       ! third-order
!!            ve(-1,m) = ( 35.0*vn(1,m) - 42.0*vn(2,m) +  15.0*vn(3,m) ) / 8.0                                       ! third-order
!!            ve(0 ,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +   3.0*vn(3,m) ) / 8.0                                       ! third-order
!!            ve(1 ,m) = (  3.0*vn(1,m) +  6.0*vn(2,m) -       vn(3,m) ) / 8.0                                       ! third-order
!!            ve(2 ,m) = (     -vn(1,m) +  9.0*vn(2,m) +   9.0*vn(3,m) - vn(4,m) ) / 16.0                            ! fourth-order
!!            ve(3 ,m) = (150.0*(vn(4,m) + vn(3,m)) - 25.0*(vn(5,m) + vn(2,m)) + 3.0*(vn(6,m) + vn(1,m)) )/256.0     ! sixth-order
            
        end do
    
    else
        st = 4

        do m=nst,ned
            !!ve(2,m) = (      -vn(1,m) +   9.0*vn(2,m) +   9.0*vn(3,m) -       vn(4,m) ) / 16.0                   ! fourth-order
            !!ve(1,m) = (   5.0*vn(1,m) +  15.0*vn(2,m) -   5.0*vn(3,m) +       vn(4,m) ) / 16.0                   ! fourth-order
            !!ve(0,m) = (  35.0*vn(1,m) -  35.0*vn(2,m) +  21.0*vn(3,m) -   5.0*vn(4,m) ) / 16.0                   ! fourth-order
            !!ve(0,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                                       ! third-order
            ve(2,m) = (  -5.0*vn(1,m) +  60.0*vn(2,m) +  90.0*vn(3,m) -  20.0*vn(4,m) +  3.0*vn(5,m) ) / 128.0   ! fifth-order
            ve(1,m) = (  35.0*vn(1,m) + 140.0*vn(2,m) -  70.0*vn(3,m) +  28.0*vn(4,m) -  5.0*vn(5,m) ) / 128.0   ! fifth-order
            ve(0,m) = ( 315.0*vn(1,m) - 420.0*vn(2,m) + 378.0*vn(3,m) - 180.0*vn(4,m) + 35.0*vn(5,m) ) / 128.0   ! fifth-order
            !!ve(1,m) = ( 3.0*vn(1,m) + 6.0*vn(2,m) -     vn(3,m) ) / 8.0                                          ! third-order
            !!ve(2,m) = (    -vn(1,m) + 9.0*vn(2,m) + 9.0*vn(3,m) - vn(4,m) ) / 16.0                               ! fourth-order
            ve(3 ,m) = (150.0*(vn(4,m) + vn(3,m)) - 25.0*(vn(5,m) + vn(2,m)) + 3.0*(vn(6,m) + vn(1,m)) )/256.0   ! sixth-order
        end do
    end if
    
    if (nfe == nbc_inter_scheme) then
        ed = ni + nge
    
    else if (nfe == nbc_intbc_scheme) then
    
        ed = ni - 4

        do m=nst,ned
            ve(ni+2,m) = &
            ( 3003.0*vn(ni,m) - 8580.0*vn(ni-1,m) + 10010.0*vn(ni-2,m) - 5460.0*vn(ni-3,m) + 1155.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            ve(ni+1,m) = &
            ( 1155.0*vn(ni,m) - 2772.0*vn(ni-1,m) +  2970.0*vn(ni-2,m) - 1540.0*vn(ni-3,m) +  315.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            ve(ni  ,m) = &
            (  315.0*vn(ni,m) -  420.0*vn(ni-1,m) +   378.0*vn(ni-2,m) -  180.0*vn(ni-3,m) +   35.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            ve(ni-1,m) = &
            (   35.0*vn(ni,m) +  140.0*vn(ni-1,m) -    70.0*vn(ni-2,m) +   28.0*vn(ni-3,m) -    5.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            ve(ni-2,m) = &
            (   -5.0*vn(ni,m) +   60.0*vn(ni-1,m) +    90.0*vn(ni-2,m) -   20.0*vn(ni-3,m) +    3.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            !!ve(ni-2,m) = (       -vn(ni,m) +    9.0*vn(ni-1,m) +     9.0*vn(ni-2,m) -        vn(ni-3,m) ) / 16.0                        ! fourth-order
            ve(ni-3,m) = &
            ( 150.0*(vn(ni-2,m) + vn(ni-3,m)) - 25.0*(vn(ni-1,m) + vn(ni-4,m)) + 3.0*(vn(ni ,m) + vn(ni-5,m)) )/256.0      ! sixth-order
!!            ve(ni+2,m) = (  63.0*vn(ni,m) -  90.0*vn(ni-1,m) +  35.0*vn(ni-2,m) ) / 8.0                                            ! third-order
!!            ve(ni+1,m) = (  35.0*vn(ni,m) -  42.0*vn(ni-1,m) +  15.0*vn(ni-2,m) ) / 8.0                                            ! third-order
!!            ve(ni  ,m) = (  15.0*vn(ni,m) -  10.0*vn(ni-1,m) +   3.0*vn(ni-2,m) ) / 8.0                                            ! third-order
!!            ve(ni-1,m) = (   3.0*vn(ni,m) +   6.0*vn(ni-1,m) -       vn(ni-2,m) ) / 8.0                                            ! third-order
!!            ve(ni-2,m) = (      -vn(ni,m) +   9.0*vn(ni-1,m) +   9.0*vn(ni-2,m) - vn(ni-3,m) ) / 16.0                              ! fourth-order
!!            ve(ni-3,m) = ( 150.0*(vn(ni-2,m) + vn(ni-3,m)) - 25.0*(vn(ni-1,m) + vn(ni-4,m)) + 3.0*(vn(ni ,m) + vn(ni-5,m)) )/256.0 ! sixth-order
        end do    
    
    else
        ed = ni - 4

        do m=nst,ned
            !!ve(ni-2,m) = (      -vn(ni,m) +   9.0*vn(ni-1,m) +   9.0*vn(ni-2,m) -       vn(ni-3,m) ) / 16.0                      ! fourth-order
            !!ve(ni-1,m) = (   5.0*vn(ni,m) +  15.0*vn(ni-1,m) -   5.0*vn(ni-2,m) +       vn(ni-3,m) ) / 16.0                      ! fourth-order
            !!ve(ni  ,m) = (  35.0*vn(ni,m) -  35.0*vn(ni-1,m) +  21.0*vn(ni-2,m) -   5.0*vn(ni-3,m) ) / 16.0                      ! fourth-order
            !!ve(ni  ,m) = (  15.0*vn(ni,m) -  10.0*vn(ni-1,m) +   3.0*vn(ni-2,m) ) / 8.0                                          ! third-order
            ve(ni-2,m) = &
            (  -5.0*vn(ni,m) +  60.0*vn(ni-1,m) +  90.0*vn(ni-2,m) -  20.0*vn(ni-3,m) +  3.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            ve(ni-1,m) = &
            (  35.0*vn(ni,m) + 140.0*vn(ni-1,m) -  70.0*vn(ni-2,m) +  28.0*vn(ni-3,m) -  5.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            ve(ni  ,m) = &
            ( 315.0*vn(ni,m) - 420.0*vn(ni-1,m) + 378.0*vn(ni-2,m) - 180.0*vn(ni-3,m) + 35.0*vn(ni-4,m) ) / 128.0     ! fifth-order
            !!ve(ni-1,m) = (  3.0*vn(ni,m) +   6.0*vn(ni-1,m) -       vn(ni-2,m) ) / 8.0                                             ! third-order
            !!ve(ni-2,m) = (     -vn(ni,m) +   9.0*vn(ni-1,m) +   9.0*vn(ni-2,m) -       vn(ni-3,m) ) / 16.0                         ! fourth-order
            ve(ni-3,m) = &
            ( 150.0*(vn(ni-2,m) + vn(ni-3,m)) - 25.0*(vn(ni-1,m) + vn(ni-4,m)) + 3.0*(vn(ni ,m) + vn(ni-5,m)) )/256.0 ! sixth-order            
        end do
    end if

    do m=nst,ned
        do i=st,ed
            ve(i,m) =    ( 1225.0*( vn(i+1,m) + vn(i  ,m) )            -  &
                            245.0*( vn(i+2,m) + vn(i-1,m) )            +  &
                             49.0*( vn(i+3,m) + vn(i-2,m) )            -  &
                              5.0*( vn(i+4,m) + vn(i-3,m) ) ) / 2048.0 
        end do
    end do
    
end subroutine ve_via_node8

subroutine ve_via_node8cc(ni,nst,ned,ngn,vn,nfs,nfe,nge,ve)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_inter_scheme,nbc_intbc_scheme
    implicit none
    integer(kind_int),intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  ,intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)  :: nfs,nfe
    integer(kind_int),intent(in)  :: nge
    real(kind_real)  ,intent(out) :: ve(1-ngn:ni+ngn,nst:ned)
    integer(kind_int) :: i,m,st,ed
    
    if (ni < 4) then
        call ve_via_node2cc(ni,nst,ned,ngn,vn,nfs,nfe,nge,ve)
        return
    end if

    if (nfs == nbc_inter_scheme) then
        st = 1 - nge
    
    else if (nfs == nbc_intbc_scheme) then
    
        st = 5
        
        do m=nst,ned
        
            ve(-1,m) = (3003.0*vn(1,m) - 8580.0*vn(2,m) + 10010.0*vn(3,m) - 5460.0*vn(4,m) + 1155.0*vn(5,m) ) / 128.0   ! fifth-order
            ve(0,m) = (1155.0*vn(1,m) - 2772.0*vn(2,m) +  2970.0*vn(3,m) - 1540.0*vn(4,m) +  315.0*vn(5,m) ) / 128.0   ! fifth-order
            ve(1 ,m) = ( 315.0*vn(1,m) -  420.0*vn(2,m) +   378.0*vn(3,m) -  180.0*vn(4,m) +   35.0*vn(5,m) ) / 128.0   ! fifth-order
            ve(2 ,m) = (  35.0*vn(1,m) +  140.0*vn(2,m) -    70.0*vn(3,m) +   28.0*vn(4,m) -    5.0*vn(5,m) ) / 128.0   ! fifth-order
            ve(3 ,m) = (  -5.0*vn(1,m) +   60.0*vn(2,m) +    90.0*vn(3,m) -   20.0*vn(4,m) +    3.0*vn(5,m) ) / 128.0
            !!ve(3 ,m) = (      -vn(1,m) +    9.0*vn(2,m) +     9.0*vn(3,m) -        vn(4,m) ) / 16.0                     ! fourth-order            
            ve(4 ,m) = (150.0*(vn(4,m) + vn(3,m)) - 25.0*(vn(5,m) + vn(2,m)) + 3.0*(vn(6,m) + vn(1,m)) )/256.0          ! sixth-order
!!            ve(-1,m) = ( 63.0*vn(1,m) - 90.0*vn(2,m) +  35.0*vn(3,m) ) / 8.0                                       ! third-order
!!            ve(0,m) = ( 35.0*vn(1,m) - 42.0*vn(2,m) +  15.0*vn(3,m) ) / 8.0                                       ! third-order
!!            ve(1 ,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +   3.0*vn(3,m) ) / 8.0                                       ! third-order
!!            ve(2 ,m) = (  3.0*vn(1,m) +  6.0*vn(2,m) -       vn(3,m) ) / 8.0                                       ! third-order
!!            ve(3 ,m) = (     -vn(1,m) +  9.0*vn(2,m) +   9.0*vn(3,m) - vn(4,m) ) / 16.0                            ! fourth-order
!!            ve(4 ,m) = (150.0*(vn(4,m) + vn(3,m)) - 25.0*(vn(5,m) + vn(2,m)) + 3.0*(vn(6,m) + vn(1,m)) )/256.0     ! sixth-order
            
        end do
    
    else
        st = 5

        do m=nst,ned
            !!ve(3,m) = (      -vn(1,m) +   9.0*vn(2,m) +   9.0*vn(3,m) -       vn(4,m) ) / 16.0                   ! fourth-order
            !!ve(2,m) = (   5.0*vn(1,m) +  15.0*vn(2,m) -   5.0*vn(3,m) +       vn(4,m) ) / 16.0                   ! fourth-order
            !!ve(1,m) = (  35.0*vn(1,m) -  35.0*vn(2,m) +  21.0*vn(3,m) -   5.0*vn(4,m) ) / 16.0                   ! fourth-order
            !!ve(1,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                                       ! third-order
            ve(3,m) = (  -5.0*vn(1,m) +  60.0*vn(2,m) +  90.0*vn(3,m) -  20.0*vn(4,m) +  3.0*vn(5,m) ) / 128.0   ! fifth-order
            ve(2,m) = (  35.0*vn(1,m) + 140.0*vn(2,m) -  70.0*vn(3,m) +  28.0*vn(4,m) -  5.0*vn(5,m) ) / 128.0   ! fifth-order
            ve(1,m) = ( 315.0*vn(1,m) - 420.0*vn(2,m) + 378.0*vn(3,m) - 180.0*vn(4,m) + 35.0*vn(5,m) ) / 128.0   ! fifth-order
            !!ve(2,m) = ( 3.0*vn(1,m) + 6.0*vn(2,m) -     vn(3,m) ) / 8.0                                          ! third-order
            !!ve(3,m) = (    -vn(1,m) + 9.0*vn(2,m) + 9.0*vn(3,m) - vn(4,m) ) / 16.0                               ! fourth-order
            ve(4 ,m) = (150.0*(vn(4,m) + vn(3,m)) - 25.0*(vn(5,m) + vn(2,m)) + 3.0*(vn(6,m) + vn(1,m)) )/256.0   ! sixth-order
        end do
    end if
    
    if (nfe == nbc_inter_scheme) then
        ed = ni + 1 + nge
    
    else if (nfe == nbc_intbc_scheme) then
    
        ed = ni - 3

        do m=nst,ned
            ve(ni+3,m) = &
            ( 3003.0*vn(ni,m) - 8580.0*vn(ni-1,m) + 10010.0*vn(ni-2,m) - 5460.0*vn(ni-3,m) + 1155.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            ve(ni+2,m) = &
            ( 1155.0*vn(ni,m) - 2772.0*vn(ni-1,m) +  2970.0*vn(ni-2,m) - 1540.0*vn(ni-3,m) +  315.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            ve(ni+1,m) = &
            (  315.0*vn(ni,m) -  420.0*vn(ni-1,m) +   378.0*vn(ni-2,m) -  180.0*vn(ni-3,m) +   35.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            ve(ni  ,m) = &
            (   35.0*vn(ni,m) +  140.0*vn(ni-1,m) -    70.0*vn(ni-2,m) +   28.0*vn(ni-3,m) -    5.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            ve(ni-1,m) = &
            (   -5.0*vn(ni,m) +   60.0*vn(ni-1,m) +    90.0*vn(ni-2,m) -   20.0*vn(ni-3,m) +    3.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            !!ve(ni-2,m) = (       -vn(ni,m) +    9.0*vn(ni-1,m) +     9.0*vn(ni-2,m) -        vn(ni-3,m) ) / 16.0                        ! fourth-order
            ve(ni-2,m) = &
            ( 150.0*(vn(ni-2,m) + vn(ni-3,m)) - 25.0*(vn(ni-1,m) + vn(ni-4,m)) + 3.0*(vn(ni ,m) + vn(ni-5,m)) )/256.0      ! sixth-order
!!            ve(ni+3,m) = (  63.0*vn(ni,m) -  90.0*vn(ni-1,m) +  35.0*vn(ni-2,m) ) / 8.0                                            ! third-order
!!            ve(ni+2,m) = (  35.0*vn(ni,m) -  42.0*vn(ni-1,m) +  15.0*vn(ni-2,m) ) / 8.0                                            ! third-order
!!            ve(ni+1,m) = (  15.0*vn(ni,m) -  10.0*vn(ni-1,m) +   3.0*vn(ni-2,m) ) / 8.0                                            ! third-order
!!            ve(ni  ,m) = (   3.0*vn(ni,m) +   6.0*vn(ni-1,m) -       vn(ni-2,m) ) / 8.0                                            ! third-order
!!            ve(ni-1,m) = (      -vn(ni,m) +   9.0*vn(ni-1,m) +   9.0*vn(ni-2,m) - vn(ni-3,m) ) / 16.0                              ! fourth-order
!!            ve(ni-2,m) = ( 150.0*(vn(ni-2,m) + vn(ni-3,m)) - 25.0*(vn(ni-1,m) + vn(ni-4,m)) + 3.0*(vn(ni ,m) + vn(ni-5,m)) )/256.0 ! sixth-order
        end do    
    
    else
        ed = ni - 3

        do m=nst,ned
            !!ve(ni-1,m) = (      -vn(ni,m) +   9.0*vn(ni-1,m) +   9.0*vn(ni-2,m) -       vn(ni-3,m) ) / 16.0                      ! fourth-order
            !!ve(ni  ,m) = (   5.0*vn(ni,m) +  15.0*vn(ni-1,m) -   5.0*vn(ni-2,m) +       vn(ni-3,m) ) / 16.0                      ! fourth-order
            !!ve(ni+1,m) = (  35.0*vn(ni,m) -  35.0*vn(ni-1,m) +  21.0*vn(ni-2,m) -   5.0*vn(ni-3,m) ) / 16.0                      ! fourth-order
            !!ve(ni+1,m) = (  15.0*vn(ni,m) -  10.0*vn(ni-1,m) +   3.0*vn(ni-2,m) ) / 8.0                                          ! third-order
            ve(ni-1,m) = &
            (  -5.0*vn(ni,m) +  60.0*vn(ni-1,m) +  90.0*vn(ni-2,m) -  20.0*vn(ni-3,m) +  3.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            ve(ni  ,m) = &
            (  35.0*vn(ni,m) + 140.0*vn(ni-1,m) -  70.0*vn(ni-2,m) +  28.0*vn(ni-3,m) -  5.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            ve(ni+1,m) = &
            ( 315.0*vn(ni,m) - 420.0*vn(ni-1,m) + 378.0*vn(ni-2,m) - 180.0*vn(ni-3,m) + 35.0*vn(ni-4,m) ) / 128.0     ! fifth-order
            !!ve(ni  ,m) = (  3.0*vn(ni,m) +   6.0*vn(ni-1,m) -       vn(ni-2,m) ) / 8.0                                             ! third-order
            !!ve(ni-1,m) = (     -vn(ni,m) +   9.0*vn(ni-1,m) +   9.0*vn(ni-2,m) -       vn(ni-3,m) ) / 16.0                         ! fourth-order
            ve(ni-2,m) = &
            ( 150.0*(vn(ni-2,m) + vn(ni-3,m)) - 25.0*(vn(ni-1,m) + vn(ni-4,m)) + 3.0*(vn(ni ,m) + vn(ni-5,m)) )/256.0 ! sixth-order            
        end do
    end if

    do m=nst,ned
        do i=st-1,ed-1
            ve(i+1,m) =    ( 1225.0*( vn(i+1,m) + vn(i  ,m) )            -  &
                              245.0*( vn(i+2,m) + vn(i-1,m) )            +  &
                               49.0*( vn(i+3,m) + vn(i-2,m) )            -  &
                                5.0*( vn(i+4,m) + vn(i-3,m) ) ) / 2048.0 
        end do
    end do
    
end subroutine ve_via_node8cc

subroutine ve_via_node8e(ni,nst,ned,ngn,vn,nfs,nfe,nge,ve)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_inter_scheme
    implicit none
    integer(kind_int),intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  ,intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)  :: nfs,nfe
    integer(kind_int),intent(in)  :: nge
    real(kind_real)  ,intent(out) :: ve(1-ngn:ni+ngn,nst:ned)
    integer(kind_int) :: i,m,st,ed
    
    if (ni < 4) then
        call ve_via_node2(ni,nst,ned,ngn,vn,nfs,nfe,nge,ve)
        return
    end if

    if (nfs == nbc_inter_scheme) then

        st = -nge
    
    else
        st = 4

        do m=nst,ned
            !!ve(2,m) = (      -vn(1,m) +   9.0*vn(2,m) +   9.0*vn(3,m) -       vn(4,m) ) / 16.0                   ! fourth-order
            !!ve(1,m) = (   5.0*vn(1,m) +  15.0*vn(2,m) -   5.0*vn(3,m) +       vn(4,m) ) / 16.0                   ! fourth-order
            !!ve(0,m) = (  35.0*vn(1,m) -  35.0*vn(2,m) +  21.0*vn(3,m) -   5.0*vn(4,m) ) / 16.0                   ! fourth-order
            !!ve(0,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                                       ! third-order
            ve(2,m) = (  -5.0*vn(1,m) +  60.0*vn(2,m) +  90.0*vn(3,m) -  20.0*vn(4,m) +  3.0*vn(5,m) ) / 128.0   ! fifth-order
            ve(1,m) = (  35.0*vn(1,m) + 140.0*vn(2,m) -  70.0*vn(3,m) +  28.0*vn(4,m) -  5.0*vn(5,m) ) / 128.0   ! fifth-order
            ve(0,m) = ( 315.0*vn(1,m) - 420.0*vn(2,m) + 378.0*vn(3,m) - 180.0*vn(4,m) + 35.0*vn(5,m) ) / 128.0   ! fifth-order
            !!ve(1,m) = ( 3.0*vn(1,m) + 6.0*vn(2,m) -     vn(3,m) ) / 8.0                                          ! third-order
            !!ve(2,m) = (    -vn(1,m) + 9.0*vn(2,m) + 9.0*vn(3,m) - vn(4,m) ) / 16.0                               ! fourth-order
            ve(3 ,m) = (150.0*(vn(4,m) + vn(3,m)) - 25.0*(vn(5,m) + vn(2,m)) + 3.0*(vn(6,m) + vn(1,m)) )/256.0   ! sixth-order
        end do
    end if
    
    if (nfe == nbc_inter_scheme) then

        ed = ni + nge   
    
    else
        ed = ni - 4

        do m=nst,ned
            !!ve(ni-2,m) = (      -vn(ni,m) +   9.0*vn(ni-1,m) +   9.0*vn(ni-2,m) -       vn(ni-3,m) ) / 16.0                      ! fourth-order
            !!ve(ni-1,m) = (   5.0*vn(ni,m) +  15.0*vn(ni-1,m) -   5.0*vn(ni-2,m) +       vn(ni-3,m) ) / 16.0                      ! fourth-order
            !!ve(ni  ,m) = (  35.0*vn(ni,m) -  35.0*vn(ni-1,m) +  21.0*vn(ni-2,m) -   5.0*vn(ni-3,m) ) / 16.0                      ! fourth-order
            !!ve(ni  ,m) = (  15.0*vn(ni,m) -  10.0*vn(ni-1,m) +   3.0*vn(ni-2,m) ) / 8.0                                          ! third-order
            ve(ni-2,m) = (  -5.0*vn(ni,m) +  60.0*vn(ni-1,m) +  90.0*vn(ni-2,m) -  20.0*vn(ni-3,m) +  3.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            ve(ni-1,m) = (  35.0*vn(ni,m) + 140.0*vn(ni-1,m) -  70.0*vn(ni-2,m) +  28.0*vn(ni-3,m) -  5.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            ve(ni  ,m) = ( 315.0*vn(ni,m) - 420.0*vn(ni-1,m) + 378.0*vn(ni-2,m) - 180.0*vn(ni-3,m) + 35.0*vn(ni-4,m) ) / 128.0     ! fifth-order
            !!ve(ni-1,m) = (  3.0*vn(ni,m) +   6.0*vn(ni-1,m) -       vn(ni-2,m) ) / 8.0                                             ! third-order
            !!ve(ni-2,m) = (     -vn(ni,m) +   9.0*vn(ni-1,m) +   9.0*vn(ni-2,m) -       vn(ni-3,m) ) / 16.0                         ! fourth-order
            ve(ni-3,m) = ( 150.0*(vn(ni-2,m) + vn(ni-3,m)) - 25.0*(vn(ni-1,m) + vn(ni-4,m)) + 3.0*(vn(ni ,m) + vn(ni-5,m)) )/256.0 ! sixth-order            
        end do
    end if

    do m=nst,ned
        do i=st,ed
            ve(i,m) =    ( 1225.0*( vn(i+1,m) + vn(i  ,m) )            -  &
                            245.0*( vn(i+2,m) + vn(i-1,m) )            +  &
                             49.0*( vn(i+3,m) + vn(i-2,m) )            -  &
                              5.0*( vn(i+4,m) + vn(i-3,m) ) ) / 2048.0 
        end do
    end do
    
end subroutine ve_via_node8e

subroutine ve_via_node8ecc(ni,nst,ned,ngn,vn,nfs,nfe,nge,ve)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_inter_scheme
    implicit none
    integer(kind_int),intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  ,intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)  :: nfs,nfe
    integer(kind_int),intent(in)  :: nge
    real(kind_real)  ,intent(out) :: ve(1-ngn:ni+ngn,nst:ned)
    integer(kind_int) :: i,m,st,ed
    
    if (ni < 4) then
        call ve_via_node2cc(ni,nst,ned,ngn,vn,nfs,nfe,nge,ve)
        return
    end if

    if (nfs == nbc_inter_scheme) then

        st = 1 - nge
    
    else
        st = 5

        do m=nst,ned
            !!ve(3,m) = (      -vn(1,m) +   9.0*vn(2,m) +   9.0*vn(3,m) -       vn(4,m) ) / 16.0                   ! fourth-order
            !!ve(2,m) = (   5.0*vn(1,m) +  15.0*vn(2,m) -   5.0*vn(3,m) +       vn(4,m) ) / 16.0                   ! fourth-order
            !!ve(1,m) = (  35.0*vn(1,m) -  35.0*vn(2,m) +  21.0*vn(3,m) -   5.0*vn(4,m) ) / 16.0                   ! fourth-order
            !!ve(1,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                                       ! third-order
            ve(3,m) = (  -5.0*vn(1,m) +  60.0*vn(2,m) +  90.0*vn(3,m) -  20.0*vn(4,m) +  3.0*vn(5,m) ) / 128.0   ! fifth-order
            ve(2,m) = (  35.0*vn(1,m) + 140.0*vn(2,m) -  70.0*vn(3,m) +  28.0*vn(4,m) -  5.0*vn(5,m) ) / 128.0   ! fifth-order
            ve(1,m) = ( 315.0*vn(1,m) - 420.0*vn(2,m) + 378.0*vn(3,m) - 180.0*vn(4,m) + 35.0*vn(5,m) ) / 128.0   ! fifth-order
            !!ve(2,m) = ( 3.0*vn(1,m) + 6.0*vn(2,m) -     vn(3,m) ) / 8.0                                          ! third-order
            !!ve(3,m) = (    -vn(1,m) + 9.0*vn(2,m) + 9.0*vn(3,m) - vn(4,m) ) / 16.0                               ! fourth-order
            ve(4 ,m) = (150.0*(vn(4,m) + vn(3,m)) - 25.0*(vn(5,m) + vn(2,m)) + 3.0*(vn(6,m) + vn(1,m)) )/256.0   ! sixth-order
        end do
    end if
    
    if (nfe == nbc_inter_scheme) then

        ed = ni + 1 + nge   
    
    else
        ed = ni - 3

        do m=nst,ned
            !!ve(ni-1,m) = (      -vn(ni,m) +   9.0*vn(ni-1,m) +   9.0*vn(ni-2,m) -       vn(ni-3,m) ) / 16.0                      ! fourth-order
            !!ve(ni  ,m) = (   5.0*vn(ni,m) +  15.0*vn(ni-1,m) -   5.0*vn(ni-2,m) +       vn(ni-3,m) ) / 16.0                      ! fourth-order
            !!ve(ni+1,m) = (  35.0*vn(ni,m) -  35.0*vn(ni-1,m) +  21.0*vn(ni-2,m) -   5.0*vn(ni-3,m) ) / 16.0                      ! fourth-order
            !!ve(ni+1,m) = (  15.0*vn(ni,m) -  10.0*vn(ni-1,m) +   3.0*vn(ni-2,m) ) / 8.0                                          ! third-order
            ve(ni-1,m) = (  -5.0*vn(ni,m) +  60.0*vn(ni-1,m) +  90.0*vn(ni-2,m) -  20.0*vn(ni-3,m) +  3.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            ve(ni  ,m) = (  35.0*vn(ni,m) + 140.0*vn(ni-1,m) -  70.0*vn(ni-2,m) +  28.0*vn(ni-3,m) -  5.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            ve(ni+1,m) = ( 315.0*vn(ni,m) - 420.0*vn(ni-1,m) + 378.0*vn(ni-2,m) - 180.0*vn(ni-3,m) + 35.0*vn(ni-4,m) ) / 128.0     ! fifth-order
            !!ve(ni  ,m) = (  3.0*vn(ni,m) +   6.0*vn(ni-1,m) -       vn(ni-2,m) ) / 8.0                                             ! third-order
            !!ve(ni-1,m) = (     -vn(ni,m) +   9.0*vn(ni-1,m) +   9.0*vn(ni-2,m) -       vn(ni-3,m) ) / 16.0                         ! fourth-order
            ve(ni-2,m) = ( 150.0*(vn(ni-2,m) + vn(ni-3,m)) - 25.0*(vn(ni-1,m) + vn(ni-4,m)) + 3.0*(vn(ni ,m) + vn(ni-5,m)) )/256.0 ! sixth-order            
        end do
    end if

    do m=nst,ned
        do i=st-1,ed-1
            ve(i+1,m) =    ( 1225.0*( vn(i+1,m) + vn(i  ,m) )            -  &
                              245.0*( vn(i+2,m) + vn(i-1,m) )            +  &
                               49.0*( vn(i+3,m) + vn(i-2,m) )            -  &
                                5.0*( vn(i+4,m) + vn(i-3,m) ) ) / 2048.0 
        end do
    end do
    
end subroutine ve_via_node8ecc

subroutine ve_via_scsl4(ni,nst,ned,ngn,vn,nfs,nfe,nge,ve)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_inter_scheme,nbc_intbc_scheme
    implicit none
    integer(kind_int),intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  ,intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)  :: nfs,nfe
    integer(kind_int),intent(in)  :: nge
    real(kind_real)  ,intent(out) :: ve(1-ngn:ni+ngn,nst:ned)
!    real(kind_real)               :: ve2(1-ngn:ni+ngn,nst:ned)   
    integer(kind_int) :: i,m,st,ed,ierr
    
    if (ni < 6) then
        call ve_via_node2(ni,nst,ned,ngn,vn,nfs,nfe,nge,ve)
        return
    end if    
    
!    call ve2_via_T2(ni,nst,ned,ngn,vn,nfs,nfe,nge,ve2)
    
    if (nfs == nbc_inter_scheme) then
        st = -1
    
    else if (nfs == nbc_intbc_scheme) then
    
        st = 3
        
        do m=nst,ned
        
            ve(-2,m) = (3003.0*vn(1,m) - 8580.0*vn(2,m) + 10010.0*vn(3,m) - 5460.0*vn(4,m) + 1155.0*vn(5,m) ) / 128.0   ! fifth-order
            ve(-1,m) = (1155.0*vn(1,m) - 2772.0*vn(2,m) +  2970.0*vn(3,m) - 1540.0*vn(4,m) +  315.0*vn(5,m) ) / 128.0   ! fifth-order
            ve(0 ,m) = ( 315.0*vn(1,m) -  420.0*vn(2,m) +   378.0*vn(3,m) -  180.0*vn(4,m) +   35.0*vn(5,m) ) / 128.0   ! fifth-order
            ve(1 ,m) = (  35.0*vn(1,m) +  140.0*vn(2,m) -    70.0*vn(3,m) +   28.0*vn(4,m) -    5.0*vn(5,m) ) / 128.0   ! fifth-order
            ve(2 ,m) = (      -vn(1,m) +    9.0*vn(2,m) +     9.0*vn(3,m) -        vn(4,m) ) / 16.0                     ! fourth-order
!!            ve(-2,m) = ( 63.0*vn(1,m) - 90.0*vn(2,m) +  35.0*vn(3,m) ) / 8.0                                       ! third-order
!!            ve(-1,m) = ( 35.0*vn(1,m) - 42.0*vn(2,m) +  15.0*vn(3,m) ) / 8.0                                       ! third-order
!!            ve(0 ,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +   3.0*vn(3,m) ) / 8.0                                       ! third-order
!!            ve(1 ,m) = (  3.0*vn(1,m) +  6.0*vn(2,m) -       vn(3,m) ) / 8.0                                       ! third-order
!!            ve(2 ,m) = (     -vn(1,m) +  9.0*vn(2,m) +   9.0*vn(3,m) - vn(4,m) ) / 16.0                            ! fourth-order
            
        end do
    
    else
        
        st = 3

        do m=nst,ned
            !!ve(2,m) = (      -vn(1,m) +   9.0*vn(2,m) +   9.0*vn(3,m) -       vn(4,m) ) / 16.0                   ! fourth-order
            !!ve(1,m) = (   5.0*vn(1,m) +  15.0*vn(2,m) -   5.0*vn(3,m) +       vn(4,m) ) / 16.0                   ! fourth-order
            !!ve(0,m) = (  35.0*vn(1,m) -  35.0*vn(2,m) +  21.0*vn(3,m) -   5.0*vn(4,m) ) / 16.0                   ! fourth-order
            ve(0,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                                       ! third-order
            ve(2,m) = (  -5.0*vn(1,m) +  60.0*vn(2,m) +  90.0*vn(3,m) -  20.0*vn(4,m) +  3.0*vn(5,m) ) / 128.0   ! fifth-order
            ve(1,m) = (  35.0*vn(1,m) + 140.0*vn(2,m) -  70.0*vn(3,m) +  28.0*vn(4,m) -  5.0*vn(5,m) ) / 128.0   ! fifth-order
            !!ve(0,m) = ( 315.0*vn(1,m) - 420.0*vn(2,m) + 378.0*vn(3,m) - 180.0*vn(4,m) + 35.0*vn(5,m) ) / 128.0   ! fifth-order
            !!ve(1,m) = ( 3.0*vn(1,m) + 6.0*vn(2,m) -     vn(3,m) ) / 8.0                                          ! third-order
            !!ve(2,m) = (    -vn(1,m) + 9.0*vn(2,m) + 9.0*vn(3,m) - vn(4,m) ) / 16.0                               ! fourth-order
        end do
    end if
    
    if (nfe == nbc_inter_scheme) then
        
        ed = ni + 1
    
    else if (nfe == nbc_intbc_scheme) then
    
        ed = ni - 3

        do m=nst,ned
            ve(ni+2,m) = &
            ( 3003.0*vn(ni,m) - 8580.0*vn(ni-1,m) + 10010.0*vn(ni-2,m) - 5460.0*vn(ni-3,m) + 1155.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            ve(ni+1,m) = &
            ( 1155.0*vn(ni,m) - 2772.0*vn(ni-1,m) +  2970.0*vn(ni-2,m) - 1540.0*vn(ni-3,m) +  315.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            ve(ni  ,m) = &
            (  315.0*vn(ni,m) -  420.0*vn(ni-1,m) +   378.0*vn(ni-2,m) -  180.0*vn(ni-3,m) +   35.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            ve(ni-1,m) = &
            (   35.0*vn(ni,m) +  140.0*vn(ni-1,m) -    70.0*vn(ni-2,m) +   28.0*vn(ni-3,m) -    5.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            ve(ni-2,m) = &
            (       -vn(ni,m) +    9.0*vn(ni-1,m) +     9.0*vn(ni-2,m) -        vn(ni-3,m) ) / 16.0                        ! fourth-order
!!            ve(ni+2,m) = (  63.0*vn(ni,m) -  90.0*vn(ni-1,m) +  35.0*vn(ni-2,m) ) / 8.0                                          ! third-order
!!            ve(ni+1,m) = (  35.0*vn(ni,m) -  42.0*vn(ni-1,m) +  15.0*vn(ni-2,m) ) / 8.0                                          ! third-order
!!            ve(ni  ,m) = (  15.0*vn(ni,m) -  10.0*vn(ni-1,m) +   3.0*vn(ni-2,m) ) / 8.0                                          ! third-order
!!            ve(ni-1,m) = (   3.0*vn(ni,m) +   6.0*vn(ni-1,m) -       vn(ni-2,m) ) / 8.0                                          ! third-order
!!            ve(ni-2,m) = (      -vn(ni,m) +   9.0*vn(ni-1,m) +   9.0*vn(ni-2,m) - vn(ni-3,m) ) / 16.0                            ! fourth-order
        end do    
    
    else
        ed = ni - 3

        do m=nst,ned
            !!ve(ni-2,m) = (      -vn(ni,m) +   9.0*vn(ni-1,m) +   9.0*vn(ni-2,m) -       vn(ni-3,m) ) / 16.0                      ! fourth-order
            !!ve(ni-1,m) = (   5.0*vn(ni,m) +  15.0*vn(ni-1,m) -   5.0*vn(ni-2,m) +       vn(ni-3,m) ) / 16.0                      ! fourth-order
            !!ve(ni  ,m) = (  35.0*vn(ni,m) -  35.0*vn(ni-1,m) +  21.0*vn(ni-2,m) -   5.0*vn(ni-3,m) ) / 16.0                      ! fourth-order
            ve(ni  ,m) = (  15.0*vn(ni,m) -  10.0*vn(ni-1,m) +   3.0*vn(ni-2,m) ) / 8.0                                          ! third-order
            ve(ni-2,m) = (  -5.0*vn(ni,m) +  60.0*vn(ni-1,m) +  90.0*vn(ni-2,m) -  20.0*vn(ni-3,m) +  3.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            ve(ni-1,m) = (  35.0*vn(ni,m) + 140.0*vn(ni-1,m) -  70.0*vn(ni-2,m) +  28.0*vn(ni-3,m) -  5.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            !!ve(ni  ,m) = ( 315.0*vn(ni,m) - 420.0*vn(ni-1,m) + 378.0*vn(ni-2,m) - 180.0*vn(ni-3,m) + 35.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            !!ve(ni-1,m) = (  3.0*vn(ni,m) +   6.0*vn(ni-1,m) -       vn(ni-2,m) ) / 8.0                                          ! third-order
            !!ve(ni-2,m) = (     -vn(ni,m) +   9.0*vn(ni-1,m) +   9.0*vn(ni-2,m) -       vn(ni-3,m) ) / 16.0                      ! fourth-order
        end do
    end if

    do m=nst,ned
        do i=st,ed
!            ve(i,m) = 0.5*(vn(i+1,m) + vn(i  ,m)) - &
!                      1.0/8.0*ve2(i,m)
            ve(i,m) = ( -vn(i-1,m) + 9.0*vn(i,m) + 9.0*vn(i+1,m) - vn(i+2,m) ) / 16.0
        end do
    end do

end subroutine ve_via_scsl4

subroutine ve_via_scsl4cc(ni,nst,ned,ngn,vn,nfs,nfe,nge,ve)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_inter_scheme,nbc_intbc_scheme
    implicit none
    integer(kind_int),intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  ,intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)  :: nfs,nfe
    integer(kind_int),intent(in)  :: nge
    real(kind_real)  ,intent(out) :: ve(1-ngn:ni+ngn,nst:ned)
!    real(kind_real)               :: ve2(1-ngn:ni+ngn,nst:ned)   
    integer(kind_int) :: i,m,st,ed,ierr
    
    if (ni < 6) then
        call ve_via_node2cc(ni,nst,ned,ngn,vn,nfs,nfe,nge,ve)
        return
    end if    
    
!    call ve2_via_T2(ni,nst,ned,ngn,vn,nfs,nfe,nge,ve2)
    
    if (nfs == nbc_inter_scheme) then
        st = 0
    
    else if (nfs == nbc_intbc_scheme) then
    
        st = 4
        
        do m=nst,ned
        
            ve(-1,m) = (3003.0*vn(1,m) - 8580.0*vn(2,m) + 10010.0*vn(3,m) - 5460.0*vn(4,m) + 1155.0*vn(5,m) ) / 128.0   ! fifth-order
            ve(0 ,m) = (1155.0*vn(1,m) - 2772.0*vn(2,m) +  2970.0*vn(3,m) - 1540.0*vn(4,m) +  315.0*vn(5,m) ) / 128.0   ! fifth-order
            ve(1 ,m) = ( 315.0*vn(1,m) -  420.0*vn(2,m) +   378.0*vn(3,m) -  180.0*vn(4,m) +   35.0*vn(5,m) ) / 128.0   ! fifth-order
            ve(2 ,m) = (  35.0*vn(1,m) +  140.0*vn(2,m) -    70.0*vn(3,m) +   28.0*vn(4,m) -    5.0*vn(5,m) ) / 128.0   ! fifth-order
            ve(3 ,m) = (      -vn(1,m) +    9.0*vn(2,m) +     9.0*vn(3,m) -        vn(4,m) ) / 16.0                     ! fourth-order
!!            ve(-1,m) = ( 63.0*vn(1,m) - 90.0*vn(2,m) +  35.0*vn(3,m) ) / 8.0                                       ! third-order
!!            ve(0,m) = ( 35.0*vn(1,m) - 42.0*vn(2,m) +  15.0*vn(3,m) ) / 8.0                                       ! third-order
!!            ve(1 ,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +   3.0*vn(3,m) ) / 8.0                                       ! third-order
!!            ve(2 ,m) = (  3.0*vn(1,m) +  6.0*vn(2,m) -       vn(3,m) ) / 8.0                                       ! third-order
!!            ve(3 ,m) = (     -vn(1,m) +  9.0*vn(2,m) +   9.0*vn(3,m) - vn(4,m) ) / 16.0                            ! fourth-order
            
        end do
    
    else
        
        st = 4

        do m=nst,ned
            !!ve(3,m) = (      -vn(1,m) +   9.0*vn(2,m) +   9.0*vn(3,m) -       vn(4,m) ) / 16.0                   ! fourth-order
            !!ve(2,m) = (   5.0*vn(1,m) +  15.0*vn(2,m) -   5.0*vn(3,m) +       vn(4,m) ) / 16.0                   ! fourth-order
            !!ve(1,m) = (  35.0*vn(1,m) -  35.0*vn(2,m) +  21.0*vn(3,m) -   5.0*vn(4,m) ) / 16.0                   ! fourth-order
            ve(1,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                                       ! third-order
            ve(3,m) = (  -5.0*vn(1,m) +  60.0*vn(2,m) +  90.0*vn(3,m) -  20.0*vn(4,m) +  3.0*vn(5,m) ) / 128.0   ! fifth-order
            ve(2,m) = (  35.0*vn(1,m) + 140.0*vn(2,m) -  70.0*vn(3,m) +  28.0*vn(4,m) -  5.0*vn(5,m) ) / 128.0   ! fifth-order
            !!ve(1,m) = ( 315.0*vn(1,m) - 420.0*vn(2,m) + 378.0*vn(3,m) - 180.0*vn(4,m) + 35.0*vn(5,m) ) / 128.0   ! fifth-order
            !!ve(2,m) = ( 3.0*vn(1,m) + 6.0*vn(2,m) -     vn(3,m) ) / 8.0                                          ! third-order
            !!ve(3,m) = (    -vn(1,m) + 9.0*vn(2,m) + 9.0*vn(3,m) - vn(4,m) ) / 16.0                               ! fourth-order
        end do
    end if
    
    if (nfe == nbc_inter_scheme) then
        
        ed = ni + 2
    
    else if (nfe == nbc_intbc_scheme) then
    
        ed = ni - 2

        do m=nst,ned
            ve(ni+3,m) = &
            ( 3003.0*vn(ni,m) - 8580.0*vn(ni-1,m) + 10010.0*vn(ni-2,m) - 5460.0*vn(ni-3,m) + 1155.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            ve(ni+2,m) = &
            ( 1155.0*vn(ni,m) - 2772.0*vn(ni-1,m) +  2970.0*vn(ni-2,m) - 1540.0*vn(ni-3,m) +  315.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            ve(ni+1,m) = &
            (  315.0*vn(ni,m) -  420.0*vn(ni-1,m) +   378.0*vn(ni-2,m) -  180.0*vn(ni-3,m) +   35.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            ve(ni  ,m) = &
            (   35.0*vn(ni,m) +  140.0*vn(ni-1,m) -    70.0*vn(ni-2,m) +   28.0*vn(ni-3,m) -    5.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            ve(ni-1,m) = &
            (       -vn(ni,m) +    9.0*vn(ni-1,m) +     9.0*vn(ni-2,m) -        vn(ni-3,m) ) / 16.0                        ! fourth-order
!!            ve(ni+3,m) = (  63.0*vn(ni,m) -  90.0*vn(ni-1,m) +  35.0*vn(ni-2,m) ) / 8.0                                          ! third-order
!!            ve(ni+2,m) = (  35.0*vn(ni,m) -  42.0*vn(ni-1,m) +  15.0*vn(ni-2,m) ) / 8.0                                          ! third-order
!!            ve(ni+1,m) = (  15.0*vn(ni,m) -  10.0*vn(ni-1,m) +   3.0*vn(ni-2,m) ) / 8.0                                          ! third-order
!!            ve(ni  ,m) = (   3.0*vn(ni,m) +   6.0*vn(ni-1,m) -       vn(ni-2,m) ) / 8.0                                          ! third-order
!!            ve(ni-1,m) = (      -vn(ni,m) +   9.0*vn(ni-1,m) +   9.0*vn(ni-2,m) - vn(ni-3,m) ) / 16.0                            ! fourth-order
        end do    
    
    else
        ed = ni - 2

        do m=nst,ned
            !!ve(ni-1,m) = (      -vn(ni,m) +   9.0*vn(ni-1,m) +   9.0*vn(ni-2,m) -       vn(ni-3,m) ) / 16.0                      ! fourth-order
            !!ve(ni  ,m) = (   5.0*vn(ni,m) +  15.0*vn(ni-1,m) -   5.0*vn(ni-2,m) +       vn(ni-3,m) ) / 16.0                      ! fourth-order
            !!ve(ni+1,m) = (  35.0*vn(ni,m) -  35.0*vn(ni-1,m) +  21.0*vn(ni-2,m) -   5.0*vn(ni-3,m) ) / 16.0                      ! fourth-order
            ve(ni+1,m) = (  15.0*vn(ni,m) -  10.0*vn(ni-1,m) +   3.0*vn(ni-2,m) ) / 8.0                                          ! third-order
            ve(ni-1,m) = (  -5.0*vn(ni,m) +  60.0*vn(ni-1,m) +  90.0*vn(ni-2,m) -  20.0*vn(ni-3,m) +  3.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            ve(ni  ,m) = (  35.0*vn(ni,m) + 140.0*vn(ni-1,m) -  70.0*vn(ni-2,m) +  28.0*vn(ni-3,m) -  5.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            !!ve(ni+1,m) = ( 315.0*vn(ni,m) - 420.0*vn(ni-1,m) + 378.0*vn(ni-2,m) - 180.0*vn(ni-3,m) + 35.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            !!ve(ni  ,m) = (  3.0*vn(ni,m) +   6.0*vn(ni-1,m) -       vn(ni-2,m) ) / 8.0                                          ! third-order
            !!ve(ni-1,m) = (     -vn(ni,m) +   9.0*vn(ni-1,m) +   9.0*vn(ni-2,m) -       vn(ni-3,m) ) / 16.0                      ! fourth-order
        end do
    end if

    do m=nst,ned
        do i=st-1,ed-1
!            ve(i+1,m) = 0.5*(vn(i+1,m) + vn(i  ,m)) - &
!                        1.0/8.0*ve2(i,m)
            ve(i+1,m) = ( -vn(i-1,m) + 9.0*vn(i,m) + 9.0*vn(i+1,m) - vn(i+2,m) ) / 16.0
        end do
    end do

end subroutine ve_via_scsl4cc

subroutine ve_via_scsl4e(ni,nst,ned,ngn,vn,nfs,nfe,nge,ve)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_inter_scheme
    implicit none
    integer(kind_int),intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  ,intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)  :: nfs,nfe
    integer(kind_int),intent(in)  :: nge
    real(kind_real)  ,intent(out) :: ve(1-ngn:ni+ngn,nst:ned)  
    integer(kind_int) :: i,m,st,ed,ierr
    
    if (ni < 6) then
        call ve_via_node2(ni,nst,ned,ngn,vn,nfs,nfe,nge,ve)
        return
    end if
    
    if (nfs == nbc_inter_scheme) then

        st = -nge
    
    else
        
        st = 3

        do m=nst,ned
            !!ve(2,m) = (      -vn(1,m) +   9.0*vn(2,m) +   9.0*vn(3,m) -       vn(4,m) ) / 16.0                   ! fourth-order
            !!ve(1,m) = (   5.0*vn(1,m) +  15.0*vn(2,m) -   5.0*vn(3,m) +       vn(4,m) ) / 16.0                   ! fourth-order
            !!ve(0,m) = (  35.0*vn(1,m) -  35.0*vn(2,m) +  21.0*vn(3,m) -   5.0*vn(4,m) ) / 16.0                   ! fourth-order
            ve(0,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                                       ! third-order
            ve(2,m) = (  -5.0*vn(1,m) +  60.0*vn(2,m) +  90.0*vn(3,m) -  20.0*vn(4,m) +  3.0*vn(5,m) ) / 128.0   ! fifth-order
            ve(1,m) = (  35.0*vn(1,m) + 140.0*vn(2,m) -  70.0*vn(3,m) +  28.0*vn(4,m) -  5.0*vn(5,m) ) / 128.0   ! fifth-order
            !!ve(0,m) = ( 315.0*vn(1,m) - 420.0*vn(2,m) + 378.0*vn(3,m) - 180.0*vn(4,m) + 35.0*vn(5,m) ) / 128.0   ! fifth-order
            !!ve(1,m) = ( 3.0*vn(1,m) + 6.0*vn(2,m) -     vn(3,m) ) / 8.0                                          ! third-order
            !!ve(2,m) = (    -vn(1,m) + 9.0*vn(2,m) + 9.0*vn(3,m) - vn(4,m) ) / 16.0                               ! fourth-order
        end do
    end if
    
    if (nfe == nbc_inter_scheme) then
        
        ed = ni + nge
    
    else
        ed = ni - 3

        do m=nst,ned
            !!ve(ni-2,m) = (      -vn(ni,m) +   9.0*vn(ni-1,m) +   9.0*vn(ni-2,m) -       vn(ni-3,m) ) / 16.0                      ! fourth-order
            !!ve(ni-1,m) = (   5.0*vn(ni,m) +  15.0*vn(ni-1,m) -   5.0*vn(ni-2,m) +       vn(ni-3,m) ) / 16.0                      ! fourth-order
            !!ve(ni  ,m) = (  35.0*vn(ni,m) -  35.0*vn(ni-1,m) +  21.0*vn(ni-2,m) -   5.0*vn(ni-3,m) ) / 16.0                      ! fourth-order
            ve(ni  ,m) = (  15.0*vn(ni,m) -  10.0*vn(ni-1,m) +   3.0*vn(ni-2,m) ) / 8.0                                          ! third-order
            ve(ni-2,m) = (  -5.0*vn(ni,m) +  60.0*vn(ni-1,m) +  90.0*vn(ni-2,m) -  20.0*vn(ni-3,m) +  3.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            ve(ni-1,m) = (  35.0*vn(ni,m) + 140.0*vn(ni-1,m) -  70.0*vn(ni-2,m) +  28.0*vn(ni-3,m) -  5.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            !!ve(ni  ,m) = ( 315.0*vn(ni,m) - 420.0*vn(ni-1,m) + 378.0*vn(ni-2,m) - 180.0*vn(ni-3,m) + 35.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            !!ve(ni-1,m) = (  3.0*vn(ni,m) +   6.0*vn(ni-1,m) -       vn(ni-2,m) ) / 8.0                                          ! third-order
            !!ve(ni-2,m) = (     -vn(ni,m) +   9.0*vn(ni-1,m) +   9.0*vn(ni-2,m) -       vn(ni-3,m) ) / 16.0                      ! fourth-order
        end do
    end if

    do m=nst,ned
        do i=st,ed
            ve(i,m) = ( -vn(i-1,m) + 9.0*vn(i,m) + 9.0*vn(i+1,m) - vn(i+2,m) ) / 16.0
        end do
    end do

end subroutine ve_via_scsl4e

subroutine ve_via_scsl4ecc(ni,nst,ned,ngn,vn,nfs,nfe,nge,ve)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_inter_scheme
    implicit none
    integer(kind_int),intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  ,intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)  :: nfs,nfe
    integer(kind_int),intent(in)  :: nge
    real(kind_real)  ,intent(out) :: ve(1-ngn:ni+ngn,nst:ned)  
    integer(kind_int) :: i,m,st,ed,ierr
    
    if (ni < 6) then
        call ve_via_node2cc(ni,nst,ned,ngn,vn,nfs,nfe,nge,ve)
        return
    end if
    
    if (nfs == nbc_inter_scheme) then

        st = 1 - nge
    
    else
        
        st = 4

        do m=nst,ned
            !!ve(3,m) = (      -vn(1,m) +   9.0*vn(2,m) +   9.0*vn(3,m) -       vn(4,m) ) / 16.0                   ! fourth-order
            !!ve(2,m) = (   5.0*vn(1,m) +  15.0*vn(2,m) -   5.0*vn(3,m) +       vn(4,m) ) / 16.0                   ! fourth-order
            !!ve(1,m) = (  35.0*vn(1,m) -  35.0*vn(2,m) +  21.0*vn(3,m) -   5.0*vn(4,m) ) / 16.0                   ! fourth-order
            ve(1,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                                       ! third-order
            ve(3,m) = (  -5.0*vn(1,m) +  60.0*vn(2,m) +  90.0*vn(3,m) -  20.0*vn(4,m) +  3.0*vn(5,m) ) / 128.0   ! fifth-order
            ve(2,m) = (  35.0*vn(1,m) + 140.0*vn(2,m) -  70.0*vn(3,m) +  28.0*vn(4,m) -  5.0*vn(5,m) ) / 128.0   ! fifth-order
            !!ve(1,m) = ( 315.0*vn(1,m) - 420.0*vn(2,m) + 378.0*vn(3,m) - 180.0*vn(4,m) + 35.0*vn(5,m) ) / 128.0   ! fifth-order
            !!ve(2,m) = ( 3.0*vn(1,m) + 6.0*vn(2,m) -     vn(3,m) ) / 8.0                                          ! third-order
            !!ve(3,m) = (    -vn(1,m) + 9.0*vn(2,m) + 9.0*vn(3,m) - vn(4,m) ) / 16.0                               ! fourth-order
        end do
    end if
    
    if (nfe == nbc_inter_scheme) then
        
        ed = ni + 1 + nge
    
    else
        ed = ni - 2

        do m=nst,ned
            !!ve(ni-1,m) = (      -vn(ni,m) +   9.0*vn(ni-1,m) +   9.0*vn(ni-2,m) -       vn(ni-3,m) ) / 16.0                      ! fourth-order
            !!ve(ni  ,m) = (   5.0*vn(ni,m) +  15.0*vn(ni-1,m) -   5.0*vn(ni-2,m) +       vn(ni-3,m) ) / 16.0                      ! fourth-order
            !!ve(ni+1,m) = (  35.0*vn(ni,m) -  35.0*vn(ni-1,m) +  21.0*vn(ni-2,m) -   5.0*vn(ni-3,m) ) / 16.0                      ! fourth-order
            ve(ni+1,m) = (  15.0*vn(ni,m) -  10.0*vn(ni-1,m) +   3.0*vn(ni-2,m) ) / 8.0                                          ! third-order
            ve(ni-1,m) = (  -5.0*vn(ni,m) +  60.0*vn(ni-1,m) +  90.0*vn(ni-2,m) -  20.0*vn(ni-3,m) +  3.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            ve(ni  ,m) = (  35.0*vn(ni,m) + 140.0*vn(ni-1,m) -  70.0*vn(ni-2,m) +  28.0*vn(ni-3,m) -  5.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            !!ve(ni+1,m) = ( 315.0*vn(ni,m) - 420.0*vn(ni-1,m) + 378.0*vn(ni-2,m) - 180.0*vn(ni-3,m) + 35.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            !!ve(ni  ,m) = (  3.0*vn(ni,m) +   6.0*vn(ni-1,m) -       vn(ni-2,m) ) / 8.0                                          ! third-order
            !!ve(ni-1,m) = (     -vn(ni,m) +   9.0*vn(ni-1,m) +   9.0*vn(ni-2,m) -       vn(ni-3,m) ) / 16.0                      ! fourth-order
        end do
    end if

    do m=nst,ned
        do i=st-1,ed-1
            ve(i+1,m) = ( -vn(i-1,m) + 9.0*vn(i,m) + 9.0*vn(i+1,m) - vn(i+2,m) ) / 16.0
        end do
    end do

end subroutine ve_via_scsl4ecc

subroutine ve_via_scsl6(ni,nst,ned,ngn,vn,nfs,nfe,nge,ve)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_inter_scheme,nbc_intbc_scheme
    implicit none
    integer(kind_int),intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  ,intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)  :: nfs,nfe
    integer(kind_int),intent(in)  :: nge
    real(kind_real)  ,intent(out) :: ve(1-ngn:ni+ngn,nst:ned) 
    integer(kind_int) :: i,m,st,ed,ierr
    
    if (ni < 8) then
        call ve_via_node2(ni,nst,ned,ngn,vn,nfs,nfe,nge,ve)
        return
    end if
    
    if (nfs == nbc_inter_scheme) then
        st = -2
    
    else if (nfs == nbc_intbc_scheme) then
    
        st = 4
        
        do m=nst,ned
        
            ve(-2,m) = (3003.0*vn(1,m) - 8580.0*vn(2,m) + 10010.0*vn(3,m) - 5460.0*vn(4,m) + 1155.0*vn(5,m) ) / 128.0   ! fifth-order
            ve(-1,m) = (1155.0*vn(1,m) - 2772.0*vn(2,m) +  2970.0*vn(3,m) - 1540.0*vn(4,m) +  315.0*vn(5,m) ) / 128.0   ! fifth-order
            ve(0 ,m) = ( 315.0*vn(1,m) -  420.0*vn(2,m) +   378.0*vn(3,m) -  180.0*vn(4,m) +   35.0*vn(5,m) ) / 128.0   ! fifth-order
            ve(1 ,m) = (  35.0*vn(1,m) +  140.0*vn(2,m) -    70.0*vn(3,m) +   28.0*vn(4,m) -    5.0*vn(5,m) ) / 128.0   ! fifth-order
            ve(2 ,m) = (      -vn(1,m) +    9.0*vn(2,m) +     9.0*vn(3,m) -        vn(4,m) ) / 16.0                     ! fourth-order
            ve(3 ,m) = (150.0*(vn(4,m) + vn(3,m)) - 25.0*(vn(5,m) + vn(2,m)) + 3.0*(vn(6,m) + vn(1,m)) )/256.0          ! sixth-order
!!            ve(-2,m) = ( 63.0*vn(1,m) - 90.0*vn(2,m) +  35.0*vn(3,m) ) / 8.0                                       ! third-order
!!            ve(-1,m) = ( 35.0*vn(1,m) - 42.0*vn(2,m) +  15.0*vn(3,m) ) / 8.0                                       ! third-order
!!            ve(0 ,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +   3.0*vn(3,m) ) / 8.0                                       ! third-order
!!            ve(1 ,m) = (  3.0*vn(1,m) +  6.0*vn(2,m) -       vn(3,m) ) / 8.0                                       ! third-order
!!            ve(2 ,m) = (     -vn(1,m) +  9.0*vn(2,m) +   9.0*vn(3,m) - vn(4,m) ) / 16.0                            ! fourth-order
            
        end do
    
    else
        
        st = 4

        do m=nst,ned
            !!ve(2,m) = (      -vn(1,m) +   9.0*vn(2,m) +   9.0*vn(3,m) -       vn(4,m) ) / 16.0                   ! fourth-order
            !!ve(1,m) = (   5.0*vn(1,m) +  15.0*vn(2,m) -   5.0*vn(3,m) +       vn(4,m) ) / 16.0                   ! fourth-order
            !!ve(0,m) = (  35.0*vn(1,m) -  35.0*vn(2,m) +  21.0*vn(3,m) -   5.0*vn(4,m) ) / 16.0                   ! fourth-order
            ve(0,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                                       ! third-order
            ve(2,m) = (  -5.0*vn(1,m) +  60.0*vn(2,m) +  90.0*vn(3,m) -  20.0*vn(4,m) +  3.0*vn(5,m) ) / 128.0   ! fifth-order
            ve(1,m) = (  35.0*vn(1,m) + 140.0*vn(2,m) -  70.0*vn(3,m) +  28.0*vn(4,m) -  5.0*vn(5,m) ) / 128.0   ! fifth-order
            ve(3,m) = (150.0*(vn(4,m) + vn(3,m)) - 25.0*(vn(5,m) + vn(2,m)) + 3.0*(vn(6,m) + vn(1,m)) )/256.0    ! sixth-order
            !!ve(0,m) = ( 315.0*vn(1,m) - 420.0*vn(2,m) + 378.0*vn(3,m) - 180.0*vn(4,m) + 35.0*vn(5,m) ) / 128.0   ! fifth-order
            !!ve(1,m) = ( 3.0*vn(1,m) + 6.0*vn(2,m) -     vn(3,m) ) / 8.0                                          ! third-order
            !!ve(2,m) = (    -vn(1,m) + 9.0*vn(2,m) + 9.0*vn(3,m) - vn(4,m) ) / 16.0                               ! fourth-order
        end do
    end if
    
    if (nfe == nbc_inter_scheme) then
        
        ed = ni + 2
    
    else if (nfe == nbc_intbc_scheme) then
    
        ed = ni - 4

        do m=nst,ned
            ve(ni+2,m) = &
            ( 3003.0*vn(ni,m) - 8580.0*vn(ni-1,m) + 10010.0*vn(ni-2,m) - 5460.0*vn(ni-3,m) + 1155.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            ve(ni+1,m) = &
            ( 1155.0*vn(ni,m) - 2772.0*vn(ni-1,m) +  2970.0*vn(ni-2,m) - 1540.0*vn(ni-3,m) +  315.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            ve(ni  ,m) = &
            (  315.0*vn(ni,m) -  420.0*vn(ni-1,m) +   378.0*vn(ni-2,m) -  180.0*vn(ni-3,m) +   35.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            ve(ni-1,m) = &
            (   35.0*vn(ni,m) +  140.0*vn(ni-1,m) -    70.0*vn(ni-2,m) +   28.0*vn(ni-3,m) -    5.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            ve(ni-2,m) = &
            (       -vn(ni,m) +    9.0*vn(ni-1,m) +     9.0*vn(ni-2,m) -        vn(ni-3,m) ) / 16.0                        ! fourth-order
            ve(ni-3,m) = &
            ( 150.0*(vn(ni-2,m) + vn(ni-3,m)) - 25.0*(vn(ni-1,m) + vn(ni-4,m)) + 3.0*(vn(ni ,m) + vn(ni-5,m)) )/256.0      ! sixth-order
!!            ve(ni+2,m) = (  63.0*vn(ni,m) -  90.0*vn(ni-1,m) +  35.0*vn(ni-2,m) ) / 8.0                                          ! third-order
!!            ve(ni+1,m) = (  35.0*vn(ni,m) -  42.0*vn(ni-1,m) +  15.0*vn(ni-2,m) ) / 8.0                                          ! third-order
!!            ve(ni  ,m) = (  15.0*vn(ni,m) -  10.0*vn(ni-1,m) +   3.0*vn(ni-2,m) ) / 8.0                                          ! third-order
!!            ve(ni-1,m) = (   3.0*vn(ni,m) +   6.0*vn(ni-1,m) -       vn(ni-2,m) ) / 8.0                                          ! third-order
!!            ve(ni-2,m) = (      -vn(ni,m) +   9.0*vn(ni-1,m) +   9.0*vn(ni-2,m) - vn(ni-3,m) ) / 16.0                            ! fourth-order
        end do    
    
    else
        ed = ni - 4

        do m=nst,ned
            !!ve(ni-2,m) = (      -vn(ni,m) +   9.0*vn(ni-1,m) +   9.0*vn(ni-2,m) -       vn(ni-3,m) ) / 16.0                      ! fourth-order
            !!ve(ni-1,m) = (   5.0*vn(ni,m) +  15.0*vn(ni-1,m) -   5.0*vn(ni-2,m) +       vn(ni-3,m) ) / 16.0                      ! fourth-order
            !!ve(ni  ,m) = (  35.0*vn(ni,m) -  35.0*vn(ni-1,m) +  21.0*vn(ni-2,m) -   5.0*vn(ni-3,m) ) / 16.0                      ! fourth-order
            ve(ni  ,m) = (  15.0*vn(ni,m) -  10.0*vn(ni-1,m) +   3.0*vn(ni-2,m) ) / 8.0                                            ! third-order
            ve(ni-2,m) = (  -5.0*vn(ni,m) +  60.0*vn(ni-1,m) +  90.0*vn(ni-2,m) -  20.0*vn(ni-3,m) +  3.0*vn(ni-4,m) ) / 128.0     ! fifth-order
            ve(ni-1,m) = (  35.0*vn(ni,m) + 140.0*vn(ni-1,m) -  70.0*vn(ni-2,m) +  28.0*vn(ni-3,m) -  5.0*vn(ni-4,m) ) / 128.0     ! fifth-order
            ve(ni-3,m) = ( 150.0*(vn(ni-2,m) + vn(ni-3,m)) - 25.0*(vn(ni-1,m) + vn(ni-4,m)) + 3.0*(vn(ni ,m) + vn(ni-5,m)) )/256.0 ! sixth-order
            !!ve(ni  ,m) = ( 315.0*vn(ni,m) - 420.0*vn(ni-1,m) + 378.0*vn(ni-2,m) - 180.0*vn(ni-3,m) + 35.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            !!ve(ni-1,m) = (  3.0*vn(ni,m) +   6.0*vn(ni-1,m) -       vn(ni-2,m) ) / 8.0                                          ! third-order
            !!ve(ni-2,m) = (     -vn(ni,m) +   9.0*vn(ni-1,m) +   9.0*vn(ni-2,m) -       vn(ni-3,m) ) / 16.0                      ! fourth-order
        end do
    end if

    do m=nst,ned
        do i=st,ed
            
            ve(i,m) = ( 150.0*(vn(i+1,m) + vn(i  ,m)) - &
                         25.0*(vn(i+2,m) + vn(i-1,m)) + &
                          3.0*(vn(i+3,m) + vn(i-2,m)) )/256.0
            
        end do
    end do

end subroutine ve_via_scsl6

subroutine ve_via_scsl6cc(ni,nst,ned,ngn,vn,nfs,nfe,nge,ve)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_inter_scheme,nbc_intbc_scheme
    implicit none
    integer(kind_int),intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  ,intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)  :: nfs,nfe
    integer(kind_int),intent(in)  :: nge
    real(kind_real)  ,intent(out) :: ve(1-ngn:ni+ngn,nst:ned) 
    integer(kind_int) :: i,m,st,ed,ierr
    
    if (ni < 8) then
        call ve_via_node2cc(ni,nst,ned,ngn,vn,nfs,nfe,nge,ve)
        return
    end if
    
    if (nfs == nbc_inter_scheme) then
        st = -1
    
    else if (nfs == nbc_intbc_scheme) then
    
        st = 5
        
        do m=nst,ned
        
            ve(-1,m) = (3003.0*vn(1,m) - 8580.0*vn(2,m) + 10010.0*vn(3,m) - 5460.0*vn(4,m) + 1155.0*vn(5,m) ) / 128.0   ! fifth-order
            ve(0,m) = (1155.0*vn(1,m) - 2772.0*vn(2,m) +  2970.0*vn(3,m) - 1540.0*vn(4,m) +  315.0*vn(5,m) ) / 128.0   ! fifth-order
            ve(1 ,m) = ( 315.0*vn(1,m) -  420.0*vn(2,m) +   378.0*vn(3,m) -  180.0*vn(4,m) +   35.0*vn(5,m) ) / 128.0   ! fifth-order
            ve(2 ,m) = (  35.0*vn(1,m) +  140.0*vn(2,m) -    70.0*vn(3,m) +   28.0*vn(4,m) -    5.0*vn(5,m) ) / 128.0   ! fifth-order
            ve(3 ,m) = (      -vn(1,m) +    9.0*vn(2,m) +     9.0*vn(3,m) -        vn(4,m) ) / 16.0                     ! fourth-order
            ve(4 ,m) = (150.0*(vn(4,m) + vn(3,m)) - 25.0*(vn(5,m) + vn(2,m)) + 3.0*(vn(6,m) + vn(1,m)) )/256.0          ! sixth-order
!!            ve(-1,m) = ( 63.0*vn(1,m) - 90.0*vn(2,m) +  35.0*vn(3,m) ) / 8.0                                       ! third-order
!!            ve(0,m) = ( 35.0*vn(1,m) - 42.0*vn(2,m) +  15.0*vn(3,m) ) / 8.0                                       ! third-order
!!            ve(1 ,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +   3.0*vn(3,m) ) / 8.0                                       ! third-order
!!            ve(2 ,m) = (  3.0*vn(1,m) +  6.0*vn(2,m) -       vn(3,m) ) / 8.0                                       ! third-order
!!            ve(3 ,m) = (     -vn(1,m) +  9.0*vn(2,m) +   9.0*vn(3,m) - vn(4,m) ) / 16.0                            ! fourth-order
            
        end do
    
    else
        
        st = 5

        do m=nst,ned
            !!ve(3,m) = (      -vn(1,m) +   9.0*vn(2,m) +   9.0*vn(3,m) -       vn(4,m) ) / 16.0                   ! fourth-order
            !!ve(2,m) = (   5.0*vn(1,m) +  15.0*vn(2,m) -   5.0*vn(3,m) +       vn(4,m) ) / 16.0                   ! fourth-order
            !!ve(1,m) = (  35.0*vn(1,m) -  35.0*vn(2,m) +  21.0*vn(3,m) -   5.0*vn(4,m) ) / 16.0                   ! fourth-order
            ve(1,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                                       ! third-order
            ve(3,m) = (  -5.0*vn(1,m) +  60.0*vn(2,m) +  90.0*vn(3,m) -  20.0*vn(4,m) +  3.0*vn(5,m) ) / 128.0   ! fifth-order
            ve(2,m) = (  35.0*vn(1,m) + 140.0*vn(2,m) -  70.0*vn(3,m) +  28.0*vn(4,m) -  5.0*vn(5,m) ) / 128.0   ! fifth-order
            ve(4,m) = (150.0*(vn(4,m) + vn(3,m)) - 25.0*(vn(5,m) + vn(2,m)) + 3.0*(vn(6,m) + vn(1,m)) )/256.0    ! sixth-order
            !!ve(1,m) = ( 315.0*vn(1,m) - 420.0*vn(2,m) + 378.0*vn(3,m) - 180.0*vn(4,m) + 35.0*vn(5,m) ) / 128.0   ! fifth-order
            !!ve(2,m) = ( 3.0*vn(1,m) + 6.0*vn(2,m) -     vn(3,m) ) / 8.0                                          ! third-order
            !!ve(3,m) = (    -vn(1,m) + 9.0*vn(2,m) + 9.0*vn(3,m) - vn(4,m) ) / 16.0                               ! fourth-order
        end do
    end if
    
    if (nfe == nbc_inter_scheme) then
        
        ed = ni + 3
    
    else if (nfe == nbc_intbc_scheme) then
    
        ed = ni - 3

        do m=nst,ned
            ve(ni+3,m) = &
            ( 3003.0*vn(ni,m) - 8580.0*vn(ni-1,m) + 10010.0*vn(ni-2,m) - 5460.0*vn(ni-3,m) + 1155.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            ve(ni+2,m) = &
            ( 1155.0*vn(ni,m) - 2772.0*vn(ni-1,m) +  2970.0*vn(ni-2,m) - 1540.0*vn(ni-3,m) +  315.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            ve(ni+1,m) = &
            (  315.0*vn(ni,m) -  420.0*vn(ni-1,m) +   378.0*vn(ni-2,m) -  180.0*vn(ni-3,m) +   35.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            ve(ni  ,m) = &
            (   35.0*vn(ni,m) +  140.0*vn(ni-1,m) -    70.0*vn(ni-2,m) +   28.0*vn(ni-3,m) -    5.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            ve(ni-1,m) = &
            (       -vn(ni,m) +    9.0*vn(ni-1,m) +     9.0*vn(ni-2,m) -        vn(ni-3,m) ) / 16.0                        ! fourth-order
            ve(ni-2,m) = &
            ( 150.0*(vn(ni-2,m) + vn(ni-3,m)) - 25.0*(vn(ni-1,m) + vn(ni-4,m)) + 3.0*(vn(ni ,m) + vn(ni-5,m)) )/256.0      ! sixth-order
!!            ve(ni+3,m) = (  63.0*vn(ni,m) -  90.0*vn(ni-1,m) +  35.0*vn(ni-2,m) ) / 8.0                                          ! third-order
!!            ve(ni+2,m) = (  35.0*vn(ni,m) -  42.0*vn(ni-1,m) +  15.0*vn(ni-2,m) ) / 8.0                                          ! third-order
!!            ve(ni+1,m) = (  15.0*vn(ni,m) -  10.0*vn(ni-1,m) +   3.0*vn(ni-2,m) ) / 8.0                                          ! third-order
!!            ve(ni  ,m) = (   3.0*vn(ni,m) +   6.0*vn(ni-1,m) -       vn(ni-2,m) ) / 8.0                                          ! third-order
!!            ve(ni-1,m) = (      -vn(ni,m) +   9.0*vn(ni-1,m) +   9.0*vn(ni-2,m) - vn(ni-3,m) ) / 16.0                            ! fourth-order
        end do    
    
    else
        ed = ni - 3

        do m=nst,ned
            !!ve(ni-1,m) = (      -vn(ni,m) +   9.0*vn(ni-1,m) +   9.0*vn(ni-2,m) -       vn(ni-3,m) ) / 16.0                      ! fourth-order
            !!ve(ni  ,m) = (   5.0*vn(ni,m) +  15.0*vn(ni-1,m) -   5.0*vn(ni-2,m) +       vn(ni-3,m) ) / 16.0                      ! fourth-order
            !!ve(ni+1,m) = (  35.0*vn(ni,m) -  35.0*vn(ni-1,m) +  21.0*vn(ni-2,m) -   5.0*vn(ni-3,m) ) / 16.0                      ! fourth-order
            ve(ni+1,m) = (  15.0*vn(ni,m) -  10.0*vn(ni-1,m) +   3.0*vn(ni-2,m) ) / 8.0                                            ! third-order
            ve(ni-1,m) = (  -5.0*vn(ni,m) +  60.0*vn(ni-1,m) +  90.0*vn(ni-2,m) -  20.0*vn(ni-3,m) +  3.0*vn(ni-4,m) ) / 128.0     ! fifth-order
            ve(ni  ,m) = (  35.0*vn(ni,m) + 140.0*vn(ni-1,m) -  70.0*vn(ni-2,m) +  28.0*vn(ni-3,m) -  5.0*vn(ni-4,m) ) / 128.0     ! fifth-order
            ve(ni-2,m) = ( 150.0*(vn(ni-2,m) + vn(ni-3,m)) - 25.0*(vn(ni-1,m) + vn(ni-4,m)) + 3.0*(vn(ni ,m) + vn(ni-5,m)) )/256.0 ! sixth-order
            !!ve(ni+1,m) = ( 315.0*vn(ni,m) - 420.0*vn(ni-1,m) + 378.0*vn(ni-2,m) - 180.0*vn(ni-3,m) + 35.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            !!ve(ni  ,m) = (  3.0*vn(ni,m) +   6.0*vn(ni-1,m) -       vn(ni-2,m) ) / 8.0                                          ! third-order
            !!ve(ni-1,m) = (     -vn(ni,m) +   9.0*vn(ni-1,m) +   9.0*vn(ni-2,m) -       vn(ni-3,m) ) / 16.0                      ! fourth-order
        end do
    end if

    do m=nst,ned
        do i=st-1,ed-1
            
            ve(i+1,m) = ( 150.0*(vn(i+1,m) + vn(i  ,m)) - &
                           25.0*(vn(i+2,m) + vn(i-1,m)) + &
                            3.0*(vn(i+3,m) + vn(i-2,m)) )/256.0
            
        end do
    end do

end subroutine ve_via_scsl6cc

subroutine ve_via_scsl6e(ni,nst,ned,ngn,vn,nfs,nfe,nge,ve)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_inter_scheme
    implicit none
    integer(kind_int),intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  ,intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)  :: nfs,nfe
    integer(kind_int),intent(in)  :: nge
    real(kind_real)  ,intent(out) :: ve(1-ngn:ni+ngn,nst:ned) 
    integer(kind_int) :: i,m,st,ed,ierr
    
    if (ni < 8) then
        call ve_via_node2(ni,nst,ned,ngn,vn,nfs,nfe,nge,ve)
        return
    end if
    
    if (nfs == nbc_inter_scheme) then
        
        st = -nge
    
    else
        
        st = 4

        do m=nst,ned
            !!ve(2,m) = (      -vn(1,m) +   9.0*vn(2,m) +   9.0*vn(3,m) -       vn(4,m) ) / 16.0                   ! fourth-order
            !!ve(1,m) = (   5.0*vn(1,m) +  15.0*vn(2,m) -   5.0*vn(3,m) +       vn(4,m) ) / 16.0                   ! fourth-order
            !!ve(0,m) = (  35.0*vn(1,m) -  35.0*vn(2,m) +  21.0*vn(3,m) -   5.0*vn(4,m) ) / 16.0                   ! fourth-order
            ve(0,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                                       ! third-order
            ve(2,m) = (  -5.0*vn(1,m) +  60.0*vn(2,m) +  90.0*vn(3,m) -  20.0*vn(4,m) +  3.0*vn(5,m) ) / 128.0   ! fifth-order
            ve(1,m) = (  35.0*vn(1,m) + 140.0*vn(2,m) -  70.0*vn(3,m) +  28.0*vn(4,m) -  5.0*vn(5,m) ) / 128.0   ! fifth-order
            ve(3,m) = (150.0*(vn(4,m) + vn(3,m)) - 25.0*(vn(5,m) + vn(2,m)) + 3.0*(vn(6,m) + vn(1,m)) )/256.0    ! sixth-order
            !!ve(0,m) = ( 315.0*vn(1,m) - 420.0*vn(2,m) + 378.0*vn(3,m) - 180.0*vn(4,m) + 35.0*vn(5,m) ) / 128.0   ! fifth-order
            !!ve(1,m) = ( 3.0*vn(1,m) + 6.0*vn(2,m) -     vn(3,m) ) / 8.0                                          ! third-order
            !!ve(2,m) = (    -vn(1,m) + 9.0*vn(2,m) + 9.0*vn(3,m) - vn(4,m) ) / 16.0                               ! fourth-order
        end do
    end if
    
    if (nfe == nbc_inter_scheme) then
        
        ed = ni + nge 
    
    else
        ed = ni - 4

        do m=nst,ned
            !!ve(ni-2,m) = (      -vn(ni,m) +   9.0*vn(ni-1,m) +   9.0*vn(ni-2,m) -       vn(ni-3,m) ) / 16.0                      ! fourth-order
            !!ve(ni-1,m) = (   5.0*vn(ni,m) +  15.0*vn(ni-1,m) -   5.0*vn(ni-2,m) +       vn(ni-3,m) ) / 16.0                      ! fourth-order
            !!ve(ni  ,m) = (  35.0*vn(ni,m) -  35.0*vn(ni-1,m) +  21.0*vn(ni-2,m) -   5.0*vn(ni-3,m) ) / 16.0                      ! fourth-order
            ve(ni  ,m) = (  15.0*vn(ni,m) -  10.0*vn(ni-1,m) +   3.0*vn(ni-2,m) ) / 8.0                                            ! third-order
            ve(ni-2,m) = (  -5.0*vn(ni,m) +  60.0*vn(ni-1,m) +  90.0*vn(ni-2,m) -  20.0*vn(ni-3,m) +  3.0*vn(ni-4,m) ) / 128.0     ! fifth-order
            ve(ni-1,m) = (  35.0*vn(ni,m) + 140.0*vn(ni-1,m) -  70.0*vn(ni-2,m) +  28.0*vn(ni-3,m) -  5.0*vn(ni-4,m) ) / 128.0     ! fifth-order
            ve(ni-3,m) = ( 150.0*(vn(ni-2,m) + vn(ni-3,m)) - 25.0*(vn(ni-1,m) + vn(ni-4,m)) + 3.0*(vn(ni ,m) + vn(ni-5,m)) )/256.0 ! sixth-order
            !!ve(ni  ,m) = ( 315.0*vn(ni,m) - 420.0*vn(ni-1,m) + 378.0*vn(ni-2,m) - 180.0*vn(ni-3,m) + 35.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            !!ve(ni-1,m) = (  3.0*vn(ni,m) +   6.0*vn(ni-1,m) -       vn(ni-2,m) ) / 8.0                                          ! third-order
            !!ve(ni-2,m) = (     -vn(ni,m) +   9.0*vn(ni-1,m) +   9.0*vn(ni-2,m) -       vn(ni-3,m) ) / 16.0                      ! fourth-order
        end do
    end if

    do m=nst,ned
        do i=st,ed
            
            ve(i,m) = ( 150.0*(vn(i+1,m) + vn(i  ,m)) - &
                         25.0*(vn(i+2,m) + vn(i-1,m)) + &
                          3.0*(vn(i+3,m) + vn(i-2,m)) )/256.0
            
        end do
    end do

end subroutine ve_via_scsl6e

subroutine ve_via_scsl6ecc(ni,nst,ned,ngn,vn,nfs,nfe,nge,ve)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_inter_scheme
    implicit none
    integer(kind_int),intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  ,intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)  :: nfs,nfe
    integer(kind_int),intent(in)  :: nge
    real(kind_real)  ,intent(out) :: ve(1-ngn:ni+ngn,nst:ned) 
    integer(kind_int) :: i,m,st,ed,ierr
    
    if (ni < 8) then
        call ve_via_node2cc(ni,nst,ned,ngn,vn,nfs,nfe,nge,ve)
        return
    end if
    
    if (nfs == nbc_inter_scheme) then
        
        st = 1 - nge
    
    else
        
        st = 5

        do m=nst,ned
            !!ve(3,m) = (      -vn(1,m) +   9.0*vn(2,m) +   9.0*vn(3,m) -       vn(4,m) ) / 16.0                   ! fourth-order
            !!ve(2,m) = (   5.0*vn(1,m) +  15.0*vn(2,m) -   5.0*vn(3,m) +       vn(4,m) ) / 16.0                   ! fourth-order
            !!ve(1,m) = (  35.0*vn(1,m) -  35.0*vn(2,m) +  21.0*vn(3,m) -   5.0*vn(4,m) ) / 16.0                   ! fourth-order
            ve(1,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                                       ! third-order
            ve(3,m) = (  -5.0*vn(1,m) +  60.0*vn(2,m) +  90.0*vn(3,m) -  20.0*vn(4,m) +  3.0*vn(5,m) ) / 128.0   ! fifth-order
            ve(2,m) = (  35.0*vn(1,m) + 140.0*vn(2,m) -  70.0*vn(3,m) +  28.0*vn(4,m) -  5.0*vn(5,m) ) / 128.0   ! fifth-order
            ve(4,m) = (150.0*(vn(4,m) + vn(3,m)) - 25.0*(vn(5,m) + vn(2,m)) + 3.0*(vn(6,m) + vn(1,m)) )/256.0    ! sixth-order
            !!ve(1,m) = ( 315.0*vn(1,m) - 420.0*vn(2,m) + 378.0*vn(3,m) - 180.0*vn(4,m) + 35.0*vn(5,m) ) / 128.0   ! fifth-order
            !!ve(2,m) = ( 3.0*vn(1,m) + 6.0*vn(2,m) -     vn(3,m) ) / 8.0                                          ! third-order
            !!ve(3,m) = (    -vn(1,m) + 9.0*vn(2,m) + 9.0*vn(3,m) - vn(4,m) ) / 16.0                               ! fourth-order
        end do
    end if
    
    if (nfe == nbc_inter_scheme) then
        
        ed = ni + 1 + nge 
    
    else
        ed = ni - 3

        do m=nst,ned
            !!ve(ni-1,m) = (      -vn(ni,m) +   9.0*vn(ni-1,m) +   9.0*vn(ni-2,m) -       vn(ni-3,m) ) / 16.0                      ! fourth-order
            !!ve(ni  ,m) = (   5.0*vn(ni,m) +  15.0*vn(ni-1,m) -   5.0*vn(ni-2,m) +       vn(ni-3,m) ) / 16.0                      ! fourth-order
            !!ve(ni+1,m) = (  35.0*vn(ni,m) -  35.0*vn(ni-1,m) +  21.0*vn(ni-2,m) -   5.0*vn(ni-3,m) ) / 16.0                      ! fourth-order
            ve(ni+1,m) = (  15.0*vn(ni,m) -  10.0*vn(ni-1,m) +   3.0*vn(ni-2,m) ) / 8.0                                            ! third-order
            ve(ni-1,m) = (  -5.0*vn(ni,m) +  60.0*vn(ni-1,m) +  90.0*vn(ni-2,m) -  20.0*vn(ni-3,m) +  3.0*vn(ni-4,m) ) / 128.0     ! fifth-order
            ve(ni  ,m) = (  35.0*vn(ni,m) + 140.0*vn(ni-1,m) -  70.0*vn(ni-2,m) +  28.0*vn(ni-3,m) -  5.0*vn(ni-4,m) ) / 128.0     ! fifth-order
            ve(ni-2,m) = ( 150.0*(vn(ni-2,m) + vn(ni-3,m)) - 25.0*(vn(ni-1,m) + vn(ni-4,m)) + 3.0*(vn(ni ,m) + vn(ni-5,m)) )/256.0 ! sixth-order
            !!ve(ni+1,m) = ( 315.0*vn(ni,m) - 420.0*vn(ni-1,m) + 378.0*vn(ni-2,m) - 180.0*vn(ni-3,m) + 35.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            !!ve(ni  ,m) = (  3.0*vn(ni,m) +   6.0*vn(ni-1,m) -       vn(ni-2,m) ) / 8.0                                          ! third-order
            !!ve(ni-1,m) = (     -vn(ni,m) +   9.0*vn(ni-1,m) +   9.0*vn(ni-2,m) -       vn(ni-3,m) ) / 16.0                      ! fourth-order
        end do
    end if

    do m=nst,ned
        do i=st-1,ed-1
            
            ve(i+1,m) = ( 150.0*(vn(i+1,m) + vn(i  ,m)) - &
                           25.0*(vn(i+2,m) + vn(i-1,m)) + &
                            3.0*(vn(i+3,m) + vn(i-2,m)) )/256.0
            
        end do
    end do

end subroutine ve_via_scsl6ecc

subroutine dn_via_edge2(ni,nst,ned,ngn,vn,nge,ve,nfs,nfe,dn)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_inter_scheme
    implicit none
    integer(kind_int),intent(in)    :: ni,nst,ned,ngn
    real(kind_real)  ,intent(in)    :: vn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)    :: nge
    real(kind_real)  ,intent(in)    :: ve(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)    :: nfs,nfe
    real(kind_real)  ,intent(inout) :: dn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int) :: i,m,st,ed
    
    st = 1
    ed = ni

    do m=nst,ned
        do i=st,ed
            dn(i,m) = ve(i,m) - ve(i-1,m)
        end do
    end do
       
end subroutine dn_via_edge2

subroutine dn_via_edge2cc(ni,nst,ned,ngn,vn,nge,ve,nfs,nfe,dn)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_inter_scheme
    implicit none
    integer(kind_int),intent(in)    :: ni,nst,ned,ngn
    real(kind_real)  ,intent(in)    :: vn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)    :: nge
    real(kind_real)  ,intent(in)    :: ve(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)    :: nfs,nfe
    real(kind_real)  ,intent(inout) :: dn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int) :: i,m,st,ed
    
    st = 1
    ed = ni

    do m=nst,ned
        do i=st,ed
            dn(i,m) = ve(i+1,m) - ve(i,m)
        end do
    end do
       
end subroutine dn_via_edge2cc

subroutine dn_via_edge4(ni,nst,ned,ngn,vn,nge,ve,nfs,nfe,dn)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_inter_scheme
    implicit none
    integer(kind_int),intent(in)    :: ni,nst,ned,ngn
    real(kind_real)  ,intent(in)    :: vn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)    :: nge
    real(kind_real)  ,intent(in)    :: ve(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)    :: nfs,nfe
    real(kind_real)  ,intent(inout) :: dn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int) :: i,m,st,ed
    
    if (nfs == nbc_inter_scheme) then
        st = 1
    else
        st = 2

        do m=nst,ned
            dn(1,m) = ( -23.0*ve(0,m) + 21.0*ve(1,m) + 3.0*ve(2,m) - ve(3,m) ) / 24.0   ! third-order
        end do
    end if

    if (nfe == nbc_inter_scheme) then
        ed = ni
    else
        ed = ni - 1

        do m=nst,ned
            dn(ni,m) = -( -23.0*ve(ni,m) + 21.0*ve(ni-1,m) + 3.0*ve(ni-2,m) - ve(ni-3,m) ) / 24.0   ! third-order
        end do
    end if

    do m=nst,ned
        do i=st,ed
            dn(i,m) = ( 27.0*(ve(i  ,m) - ve(i-1,m)) - &
                             (ve(i+1,m) - ve(i-2,m)) ) / 24.0
        end do
    end do
      
end subroutine dn_via_edge4

subroutine dn_via_edge4cc(ni,nst,ned,ngn,vn,nge,ve,nfs,nfe,dn)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_inter_scheme
    implicit none
    integer(kind_int),intent(in)    :: ni,nst,ned,ngn
    real(kind_real)  ,intent(in)    :: vn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)    :: nge
    real(kind_real)  ,intent(in)    :: ve(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)    :: nfs,nfe
    real(kind_real)  ,intent(inout) :: dn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int) :: i,m,st,ed
    
    if (nfs == nbc_inter_scheme) then
        st = 1
    else
        st = 2

        do m=nst,ned
            dn(1,m) = ( -23.0*ve(1,m) + 21.0*ve(2,m) + 3.0*ve(3,m) - ve(4,m) ) / 24.0   ! third-order
        end do
    end if

    if (nfe == nbc_inter_scheme) then
        ed = ni
    else
        ed = ni - 1

        do m=nst,ned
            dn(ni,m) = -( -23.0*ve(ni+1,m) + 21.0*ve(ni,m) + 3.0*ve(ni-1,m) - ve(ni-2,m) ) / 24.0   ! third-order
        end do
    end if

    do m=nst,ned
        do i=st,ed
            dn(i,m) = ( 27.0*(ve(i+1,m) - ve(i  ,m)) - &
                             (ve(i+2,m) - ve(i-1,m)) ) / 24.0
        end do
    end do
      
end subroutine dn_via_edge4cc

subroutine dn_via_edge6(ni,nst,ned,ngn,vn,nge,ve,nfs,nfe,dn)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_inter_scheme
    implicit none
    integer(kind_int),intent(in)    :: ni,nst,ned,ngn
    real(kind_real)  ,intent(in)    :: vn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)    :: nge
    real(kind_real)  ,intent(in)    :: ve(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)    :: nfs,nfe
    real(kind_real)  ,intent(inout) :: dn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int) :: i,m,st,ed
    
    if (nfs == nbc_inter_scheme) then
        st = 1
    else
        st = 3

        do m=nst,ned
            !!dn(2,m) = (  71.0*ve(0,m) - 2115.0*ve(1,m) + 2070.0*ve(2,m) + 10.0*ve(3,m) - 45.0*ve(4,m) + 9.0*ve(5,m) ) / 1920.0      ! fifth-order
            !!dn(1,m) = ( -22.0*ve(0,m) +   17.0*ve(1,m) +    9.0*ve(2,m) -  5.0*ve(3,m) +      ve(4,m) ) / 24.0                    ! fourth-order
            dn(2,m) = (       ve(0,m) -   27.0*ve(1,m) +   27.0*ve(2,m) -      ve(3,m) ) / 24.0      ! fourth-order
            dn(1,m) = ( -23.0*ve(0,m) +   21.0*ve(1,m) +    3.0*ve(2,m) -      ve(3,m) ) / 24.0              ! third-order
        end do
    end if

    if (nfe == nbc_inter_scheme) then
        ed = ni
    else
        ed = ni - 2

        do m=nst,ned
            !!dn(ni-1,m) = -(  71.0*ve(ni,m) - 2115.0*ve(ni-1,m) + 2070.0*ve(ni-2,m) + 10.0*ve(ni-3,m) - 45.0*ve(ni-4,m) + 9.0*ve(ni-5,m) ) / 1920.0      ! fifth-order
            !!dn(ni  ,m) = -( -22.0*ve(ni,m) +   17.0*ve(ni-1,m) +    9.0*ve(ni-2,m) -  5.0*ve(ni-3,m) +      ve(ni-4,m) ) / 24.0                       ! fourth-order
            dn(ni-1,m) = -(       ve(ni,m) -   27.0*ve(ni-1,m) +   27.0*ve(ni-2,m) -      ve(ni-3,m) ) / 24.0      ! fourth-order
            dn(ni,m)   = -( -23.0*ve(ni,m) +   21.0*ve(ni-1,m) +    3.0*ve(ni-2,m) -      ve(ni-3,m) ) / 24.0      ! third-order
        end do
    end if

    do m=nst,ned
        do i=st,ed
            dn(i,m) = ( 2250.0*(ve(i  ,m) - ve(i-1,m)) - &
                         125.0*(ve(i+1,m) - ve(i-2,m)) + &
                           9.0*(ve(i+2,m) - ve(i-3,m)) ) / 1920.0
        end do
    end do
       
end subroutine dn_via_edge6

subroutine dn_via_edge6cc(ni,nst,ned,ngn,vn,nge,ve,nfs,nfe,dn)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_inter_scheme
    implicit none
    integer(kind_int),intent(in)    :: ni,nst,ned,ngn
    real(kind_real)  ,intent(in)    :: vn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)    :: nge
    real(kind_real)  ,intent(in)    :: ve(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)    :: nfs,nfe
    real(kind_real)  ,intent(inout) :: dn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int) :: i,m,st,ed
    
    if (nfs == nbc_inter_scheme) then
        st = 1
    else
        st = 3

        do m=nst,ned
            !!dn(2,m) = (  71.0*ve(1,m) - 2115.0*ve(2,m) + 2070.0*ve(3,m) + 10.0*ve(4,m) - 45.0*ve(5,m) + 9.0*ve(6,m) ) / 1920.0      ! fifth-order
            !!dn(1,m) = ( -22.0*ve(1,m) +   17.0*ve(2,m) +    9.0*ve(3,m) -  5.0*ve(4,m) +      ve(5,m) ) / 24.0                    ! fourth-order
            dn(2,m) = (       ve(1,m) -   27.0*ve(2,m) +   27.0*ve(3,m) -      ve(4,m) ) / 24.0      ! fourth-order
            dn(1,m) = ( -23.0*ve(1,m) +   21.0*ve(2,m) +    3.0*ve(3,m) -      ve(4,m) ) / 24.0              ! third-order
        end do
    end if

    if (nfe == nbc_inter_scheme) then
        ed = ni
    else
        ed = ni - 2

        do m=nst,ned
            !!dn(ni-1,m) = -(  71.0*ve(ni+1,m) - 2115.0*ve(ni,m) + 2070.0*ve(ni-1,m) + 10.0*ve(ni-2,m) - 45.0*ve(ni-3,m) + 9.0*ve(ni-4,m) ) / 1920.0      ! fifth-order
            !!dn(ni  ,m) = -( -22.0*ve(ni+1,m) +   17.0*ve(ni,m) +    9.0*ve(ni-1,m) -  5.0*ve(ni-2,m) +      ve(ni-3,m) ) / 24.0                       ! fourth-order
            dn(ni-1,m) = -(       ve(ni+1,m) -   27.0*ve(ni,m) +   27.0*ve(ni-1,m) -      ve(ni-2,m) ) / 24.0      ! fourth-order
            dn(ni,m)   = -( -23.0*ve(ni+1,m) +   21.0*ve(ni,m) +    3.0*ve(ni-1,m) -      ve(ni-2,m) ) / 24.0      ! third-order
        end do
    end if

    do m=nst,ned
        do i=st,ed
            dn(i,m) = ( 2250.0*(ve(i+1,m) - ve(i  ,m)) - &
                         125.0*(ve(i+2,m) - ve(i-1,m)) + &
                           9.0*(ve(i+3,m) - ve(i-2,m)) ) / 1920.0
        end do
    end do
       
end subroutine dn_via_edge6cc

subroutine dn_via_ehen4(ni,nst,ned,ngn,vn,nge,ve,nfs,nfe,dn)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_inter_scheme
    implicit none
    integer(kind_int),intent(in)    :: ni,nst,ned,ngn
    real(kind_real)  ,intent(in)    :: vn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)    :: nge
    real(kind_real)  ,intent(in)    :: ve(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)    :: nfs,nfe
    real(kind_real)  ,intent(inout) :: dn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int) :: i,m,st,ed
    real(kind_real)   :: a1,b1,ap,a2,b2

    a1 = 4.D0/3.D0
    b1 = -1.D0/6.D0 

    ap = 1.0D0
    !!ap = 4.D0/3.D0
    a2 = (16.D0 - 15.D0*ap)/24.D0
    b2 = (3.D0*ap - 4.D0)/48.D0
       
    if (nfs == nbc_inter_scheme) then
        st = 1
    else
        st = 3

        do m=nst,ned
            dn(2,m) = a1*(ve(2,m) - ve(1,m)) + &
                      b1*(vn(3,m) - vn(1,m))
            ! dn(1,m) = ( -11.0*vn(1,m) + 18.0*ve(1,m) - 9.0*vn(2,m) + 2.0*ve(2,m) ) / 3.0
            ! dn(1,m) = ( -2.0*ve(0,m) - 3.0*vn(1,m) + 6.0*ve(1,m) - vn(2,m) )/3.0
            dn(1,m) = (8.0*ve(1,m) - 15.0*vn(1,m) + 10.0*vn(2,m) - 3.0*vn(3,m) ) / 8.0
        end do
    end if

    if (nfe == nbc_inter_scheme) then
        ed = ni
    else
        ed = ni - 2

        do m=nst,ned
            dn(ni-1,m) = a1*(ve(ni-1,m) - ve(ni-2,m)) + &
                         b1*(vn(ni  ,m) - vn(ni-2,m))
            ! dn(ni,m) = -( -11.0*vn(ni,m) + 18.0*ve(ni-1,m) - 9.0*vn(ni-1,m) + 2.0*ve(ni-2,m) ) / 3.0
            ! dn(ni,m) = -( -2.0*ve(ni,m) - 3.0*vn(ni,m) + 6.0*ve(ni-1,m) - vn(ni-1,m) )/3.0
            dn(ni  ,m) = (8.0*ve(ni-1,m) + 9.0*vn(ni,m) - 22.0*vn(ni-1,m) + 5.0*vn(ni-2,m) ) / 8.0
        end do
    end if

    do m=nst,ned
        do i=st,ed
            dn(i,m) = ap*(ve(i  ,m) - ve(i-1,m)) + &
                      a2*(vn(i+1,m) - vn(i-1,m)) + &
                      b2*(vn(i+2,m) - vn(i-2,m))
        end do
    end do
       
end subroutine dn_via_ehen4

subroutine dn_via_ehen4cc(ni,nst,ned,ngn,vn,nge,ve,nfs,nfe,dn)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_inter_scheme
    implicit none
    integer(kind_int),intent(in)    :: ni,nst,ned,ngn
    real(kind_real)  ,intent(in)    :: vn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)    :: nge
    real(kind_real)  ,intent(in)    :: ve(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)    :: nfs,nfe
    real(kind_real)  ,intent(inout) :: dn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int) :: i,m,st,ed
    real(kind_real)   :: a1,b1,ap,a2,b2

    a1 = 4.D0/3.D0
    b1 = -1.D0/6.D0 

    ap = 1.0D0
    !!ap = 4.D0/3.D0
    a2 = (16.D0 - 15.D0*ap)/24.D0
    b2 = (3.D0*ap - 4.D0)/48.D0
       
    if (nfs == nbc_inter_scheme) then
        st = 1
    else
        st = 3

        do m=nst,ned
            dn(2,m) = a1*(ve(3,m) - ve(2,m)) + &
                      b1*(vn(3,m) - vn(1,m))
            ! dn(1,m) = ( -11.0*vn(1,m) + 18.0*ve(2,m) - 9.0*vn(2,m) + 2.0*ve(3,m) ) / 3.0
            dn(1,m) = ( -2.0*ve(1,m) - 3.0*vn(1,m) + 6.0*ve(2,m) - vn(2,m) )/3.0
            ! dn(1,m) = (8.0*ve(2,m) - 15.0*vn(1,m) + 10.0*vn(2,m) - 3.0*vn(3,m) ) / 8.0
        end do
    end if

    if (nfe == nbc_inter_scheme) then
        ed = ni
    else
        ed = ni - 2

        do m=nst,ned
            dn(ni-1,m) = a1*(ve(ni  ,m) - ve(ni-1,m)) + &
                         b1*(vn(ni  ,m) - vn(ni-2,m))
            ! dn(ni,m) = -( -11.0*vn(ni,m) + 18.0*ve(ni,m) - 9.0*vn(ni-1,m) + 2.0*ve(ni-1,m) ) / 3.0
            dn(ni,m) = -( -2.0*ve(ni+1,m) - 3.0*vn(ni,m) + 6.0*ve(ni,m) - vn(ni-1,m) )/3.0
            ! dn(ni  ,m) = (8.0*ve(ni,m) + 9.0*vn(ni,m) - 22.0*vn(ni-1,m) + 5.0*vn(ni-2,m) ) / 8.0
        end do
    end if

    do m=nst,ned
        do i=st,ed
            dn(i,m) = ap*(ve(i+1,m) - ve(i  ,m)) + &
                      a2*(vn(i+1,m) - vn(i-1,m)) + &
                      b2*(vn(i+2,m) - vn(i-2,m))
        end do
    end do
       
end subroutine dn_via_ehen4cc

subroutine dn_via_ehen6(ni,nst,ned,ngn,vn,nge,ve,nfs,nfe,dn)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_inter_scheme,nbc_intbc_scheme
    implicit none
    integer(kind_int),intent(in)    :: ni,nst,ned,ngn
    real(kind_real)  ,intent(in)    :: vn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)    :: nge
    real(kind_real)  ,intent(in)    :: ve(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)    :: nfs,nfe
    real(kind_real)  ,intent(inout) :: dn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int) :: i,m,st,ed
    real(kind_real)   :: a1,b1,a2,b2,c2,ap,a3,b3,c3

    if (ni < 4) then
        call dn_via_edge2(ni,nst,ned,ngn,vn,nge,ve,nfs,nfe,dn)
        return
    end if

    a1 = 4.D0/3.D0
    b1 = -1.D0/6.D0 

    a2 = 1.0D0 !!4.D0/3.D0
    b2 = (16.D0 - 15.D0*a2)/24.D0
    c2 = (3.D0*a2 - 4.D0)/48.D0

    ap = 1.D0 !!256.D0/175.D0
    a3 = (192.D0 - 175.D0*ap)/256.D0
    b3 = (35.D0*ap - 48.D0)/320.D0
    c3 = (64.D0 - 45.D0*ap)/3840.D0
    
    if (nfs == nbc_inter_scheme) then
        st = 1
    else if (nfs == nbc_intbc_scheme) then
    
        st = 4

        do m=nst,ned
            do i=1,3
                dn(i,m) = ( 2250.0*(ve(i  ,m) - ve(i-1,m)) - &
                             125.0*(ve(i+1,m) - ve(i-2,m)) + &
                               9.0*(ve(i+2,m) - ve(i-3,m)) ) / 1920.0
!!                dn(i,m) = ( 27.0*(ve(i  ,m) - ve(i-1,m)) - &
!!                                 (ve(i+1,m) - ve(i-2,m)) ) / 24.0
            end do
        end do
    
    else
        st = 4

        do m=nst,ned
            dn(3,m) = a2*(ve(3,m) - ve(2,m)) + &
                      b2*(vn(4,m) - vn(2,m)) + &
                      c2*(vn(5,m) - vn(1,m))
            dn(2,m) = a1*(ve(2,m) - ve(1,m)) + &
                      b1*(vn(3,m) - vn(1,m))
            !!dn(1,m) = (-25.0*vn(1,m) + 48.0*ve(1,m) - 36.0*vn(2,m) + 16.0*ve(2,m) - 3.0*vn(3,m) ) / 6.0
            dn(1,m) = (8.0*ve(1,m) - 15.0*vn(1,m) + 10.0*vn(2,m) - 3.0*vn(3,m) ) / 8.0
        end do
    end if

    if (nfe == nbc_inter_scheme) then
        ed = ni
    else if (nfe == nbc_intbc_scheme) then
    
        ed = ni - 3

        do m=nst,ned
            do i=ni-2,ni
                dn(i,m) = ( 2250.0*(ve(i  ,m) - ve(i-1,m)) - &
                             125.0*(ve(i+1,m) - ve(i-2,m)) + &
                               9.0*(ve(i+2,m) - ve(i-3,m)) ) / 1920.0
!!                dn(i,m) = ( 27.0*(ve(i  ,m) - ve(i-1,m)) - &
!!                             (ve(i+1,m) - ve(i-2,m)) ) / 24.0
            end do
        end do
    
    else
        ed = ni - 3

        do m=nst,ned
            dn(ni-2,m) = a2*(ve(ni-2,m) - ve(ni-3,m)) + &
                         b2*(vn(ni-1,m) - vn(ni-3,m)) + &
                         c2*(vn(ni  ,m) - vn(ni-4,m))
            dn(ni-1,m) = a1*(ve(ni-1,m) - ve(ni-2,m)) + &
                         b1*(vn(ni  ,m) - vn(ni-2,m))
            !!dn(ni  ,m) = -(-25.0*vn(ni,m) + 48.0*ve(ni-1,m) - 36.0*vn(ni-1,m) + 16.0*ve(ni-2,m) - 3.0*vn(ni-2,m) ) / 6.0
            dn(ni  ,m) = (8.0*ve(ni-1,m) + 9.0*vn(ni,m) - 22.0*vn(ni-1,m) + 5.0*vn(ni-2,m) ) / 8.0
        end do
    end if

    do m=nst,ned
        do i=st,ed
            dn(i,m) = ap*(ve(i  ,m) - ve(i-1,m)) + &
                      a3*(vn(i+1,m) - vn(i-1,m)) + &
                      b3*(vn(i+2,m) - vn(i-2,m)) + &
                      c3*(vn(i+3,m) - vn(i-3,m))
        end do
    end do
       
end subroutine dn_via_ehen6

subroutine dn_via_ehen6cc(ni,nst,ned,ngn,vn,nge,ve,nfs,nfe,dn)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_inter_scheme,nbc_intbc_scheme
    implicit none
    integer(kind_int),intent(in)    :: ni,nst,ned,ngn
    real(kind_real)  ,intent(in)    :: vn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)    :: nge
    real(kind_real)  ,intent(in)    :: ve(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)    :: nfs,nfe
    real(kind_real)  ,intent(inout) :: dn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int) :: i,m,st,ed
    real(kind_real)   :: a1,b1,a2,b2,c2,ap,a3,b3,c3

    if (ni < 4) then
        call dn_via_edge2cc(ni,nst,ned,ngn,vn,nge,ve,nfs,nfe,dn)
        return
    end if

    a1 = 4.D0/3.D0
    b1 = -1.D0/6.D0 

    a2 = 1.0D0 !!4.D0/3.D0
    b2 = (16.D0 - 15.D0*a2)/24.D0
    c2 = (3.D0*a2 - 4.D0)/48.D0

    ap = 1.D0 !!256.D0/175.D0
    a3 = (192.D0 - 175.D0*ap)/256.D0
    b3 = (35.D0*ap - 48.D0)/320.D0
    c3 = (64.D0 - 45.D0*ap)/3840.D0
    
    if (nfs == nbc_inter_scheme) then
        st = 1
    else if (nfs == nbc_intbc_scheme) then
    
        st = 4

        do m=nst,ned
            do i=1,3
                dn(i,m) = ( 2250.0*(ve(i+1,m) - ve(i  ,m)) - &
                             125.0*(ve(i+2,m) - ve(i-1,m)) + &
                               9.0*(ve(i+3,m) - ve(i-2,m)) ) / 1920.0
!!                dn(i,m) = ( 27.0*(ve(i+1,m) - ve(i  ,m)) - &
!!                                 (ve(i+2,m) - ve(i-1,m)) ) / 24.0
            end do
        end do
    
    else
        st = 4

        do m=nst,ned
            dn(3,m) = a2*(ve(4,m) - ve(3,m)) + &
                      b2*(vn(4,m) - vn(2,m)) + &
                      c2*(vn(5,m) - vn(1,m))
            dn(2,m) = a1*(ve(3,m) - ve(2,m)) + &
                      b1*(vn(3,m) - vn(1,m))
            !!dn(1,m) = (-25.0*vn(1,m) + 48.0*ve(2,m) - 36.0*vn(2,m) + 16.0*ve(3,m) - 3.0*vn(3,m) ) / 6.0
            dn(1,m) = ( -2.0*ve(1,m) - 3.0*vn(1,m) + 6.0*ve(2,m) - vn(2,m) )/3.0
            !!dn(1,m) = (8.0*ve(2,m) - 15.0*vn(1,m) + 10.0*vn(2,m) - 3.0*vn(3,m) ) / 8.0
        end do
    end if

    if (nfe == nbc_inter_scheme) then
        ed = ni
    else if (nfe == nbc_intbc_scheme) then
    
        ed = ni - 3

        do m=nst,ned
            do i=ni-2,ni
                dn(i,m) = ( 2250.0*(ve(i+1,m) - ve(i  ,m)) - &
                             125.0*(ve(i+2,m) - ve(i-1,m)) + &
                               9.0*(ve(i+3,m) - ve(i-2,m)) ) / 1920.0
!!                dn(i,m) = ( 27.0*(ve(i+1,m) - ve(i  ,m)) - &
!!                             (ve(i+2,m) - ve(i-1,m)) ) / 24.0
            end do
        end do
    
    else
        ed = ni - 3

        do m=nst,ned
            dn(ni-2,m) = a2*(ve(ni-1,m) - ve(ni-2,m)) + &
                         b2*(vn(ni-1,m) - vn(ni-3,m)) + &
                         c2*(vn(ni  ,m) - vn(ni-4,m))
            dn(ni-1,m) = a1*(ve(ni  ,m) - ve(ni-1,m)) + &
                         b1*(vn(ni  ,m) - vn(ni-2,m))
            !!dn(ni  ,m) = -(-25.0*vn(ni,m) + 48.0*ve(ni,m) - 36.0*vn(ni-1,m) + 16.0*ve(ni-1,m) - 3.0*vn(ni-2,m) ) / 6.0
            dn(ni  ,m) = -( -2.0*ve(ni+1,m) - 3.0*vn(ni,m) + 6.0*ve(ni,m) - vn(ni-1,m) )/3.0
            !!dn(ni  ,m) = (8.0*ve(ni,m) + 9.0*vn(ni,m) - 22.0*vn(ni-1,m) + 5.0*vn(ni-2,m) ) / 8.0
        end do
    end if

    do m=nst,ned
        do i=st,ed
            dn(i,m) = ap*(ve(i+1,m) - ve(i  ,m)) + &
                      a3*(vn(i+1,m) - vn(i-1,m)) + &
                      b3*(vn(i+2,m) - vn(i-2,m)) + &
                      c3*(vn(i+3,m) - vn(i-3,m))
        end do
    end do
       
end subroutine dn_via_ehen6cc

subroutine dn_via_ehen6e(ni,nst,ned,ngn,vn,nge,ve,nfs,nfe,dn)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_inter_scheme
    use mod_variables, only : nghnode
    implicit none
    integer(kind_int),intent(in)    :: ni,nst,ned,ngn
    real(kind_real)  ,intent(in)    :: vn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)    :: nge
    real(kind_real)  ,intent(in)    :: ve(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)    :: nfs,nfe
    real(kind_real)  ,intent(inout) :: dn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int) :: i,m,st,ed
    real(kind_real)   :: a1,b1,a2,b2,c2,ap,a3,b3,c3

    if (ni < 4) then
        call dn_via_edge2(ni,nst,ned,ngn,vn,nge,ve,nfs,nfe,dn)
        return
    end if

    a1 = 4.D0/3.D0
    b1 = -1.D0/6.D0 

    a2 = 1.0D0 !!4.D0/3.D0
    b2 = (16.D0 - 15.D0*a2)/24.D0
    c2 = (3.D0*a2 - 4.D0)/48.D0

    ap = 1.D0 !!256.D0/175.D0
    a3 = (192.D0 - 175.D0*ap)/256.D0
    b3 = (35.D0*ap - 48.D0)/320.D0
    c3 = (64.D0 - 45.D0*ap)/3840.D0
    
    if (nfs == nbc_inter_scheme) then

        st = 1-ngn+nghnode
    
    else
        st = 4

        do m=nst,ned
            dn(3,m) = a2*(ve(3,m) - ve(2,m)) + &
                      b2*(vn(4,m) - vn(2,m)) + &
                      c2*(vn(5,m) - vn(1,m))
            dn(2,m) = a1*(ve(2,m) - ve(1,m)) + &
                      b1*(vn(3,m) - vn(1,m))
            !!dn(1,m) = (-25.0*vn(1,m) + 48.0*ve(1,m) - 36.0*vn(2,m) + 16.0*ve(2,m) - 3.0*vn(3,m) ) / 6.0
            dn(1,m) = (8.0*ve(1,m) - 15.0*vn(1,m) + 10.0*vn(2,m) - 3.0*vn(3,m) ) / 8.0
        end do
    end if

    if (nfe == nbc_inter_scheme) then

        ed = ni+ngn-nghnode
    
    else
        ed = ni - 3

        do m=nst,ned
            dn(ni-2,m) = a2*(ve(ni-2,m) - ve(ni-3,m)) + &
                         b2*(vn(ni-1,m) - vn(ni-3,m)) + &
                         c2*(vn(ni  ,m) - vn(ni-4,m))
            dn(ni-1,m) = a1*(ve(ni-1,m) - ve(ni-2,m)) + &
                         b1*(vn(ni  ,m) - vn(ni-2,m))
            !!dn(ni  ,m) = -(-25.0*vn(ni,m) + 48.0*ve(ni-1,m) - 36.0*vn(ni-1,m) + 16.0*ve(ni-2,m) - 3.0*vn(ni-2,m) ) / 6.0
            dn(ni  ,m) = (8.0*ve(ni-1,m) + 9.0*vn(ni,m) - 22.0*vn(ni-1,m) + 5.0*vn(ni-2,m) ) / 8.0
        end do
    end if

    do m=nst,ned
        do i=st,ed
            dn(i,m) = ap*(ve(i  ,m) - ve(i-1,m)) + &
                      a3*(vn(i+1,m) - vn(i-1,m)) + &
                      b3*(vn(i+2,m) - vn(i-2,m)) + &
                      c3*(vn(i+3,m) - vn(i-3,m))
        end do
    end do
       
end subroutine dn_via_ehen6e

subroutine dn_via_ehen6ecc(ni,nst,ned,ngn,vn,nge,ve,nfs,nfe,dn)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_inter_scheme
    use mod_variables, only : nghnode
    implicit none
    integer(kind_int),intent(in)    :: ni,nst,ned,ngn
    real(kind_real)  ,intent(in)    :: vn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)    :: nge
    real(kind_real)  ,intent(in)    :: ve(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)    :: nfs,nfe
    real(kind_real)  ,intent(inout) :: dn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int) :: i,m,st,ed
    real(kind_real)   :: a1,b1,a2,b2,c2,ap,a3,b3,c3

    if (ni < 4) then
        call dn_via_edge2cc(ni,nst,ned,ngn,vn,nge,ve,nfs,nfe,dn)
        return
    end if

    a1 = 4.D0/3.D0
    b1 = -1.D0/6.D0 

    a2 = 1.0D0 !!4.D0/3.D0
    b2 = (16.D0 - 15.D0*a2)/24.D0
    c2 = (3.D0*a2 - 4.D0)/48.D0

    ap = 1.D0 !!256.D0/175.D0
    a3 = (192.D0 - 175.D0*ap)/256.D0
    b3 = (35.D0*ap - 48.D0)/320.D0
    c3 = (64.D0 - 45.D0*ap)/3840.D0
    
    if (nfs == nbc_inter_scheme) then

        st = 1-ngn+nghnode
    
    else
        st = 4

        do m=nst,ned
            dn(3,m) = a2*(ve(4,m) - ve(3,m)) + &
                      b2*(vn(4,m) - vn(2,m)) + &
                      c2*(vn(5,m) - vn(1,m))
            dn(2,m) = a1*(ve(3,m) - ve(2,m)) + &
                      b1*(vn(3,m) - vn(1,m))
            !!dn(1,m) = (-25.0*vn(1,m) + 48.0*ve(2,m) - 36.0*vn(2,m) + 16.0*ve(3,m) - 3.0*vn(3,m) ) / 6.0
            dn(1,m) = ( -2.0*ve(1,m) - 3.0*vn(1,m) + 6.0*ve(2,m) - vn(2,m) )/3.0
            !!dn(1,m) = (8.0*ve(2,m) - 15.0*vn(1,m) + 10.0*vn(2,m) - 3.0*vn(3,m) ) / 8.0
        end do
    end if

    if (nfe == nbc_inter_scheme) then

        ed = ni+ngn-nghnode
    
    else
        ed = ni - 3

        do m=nst,ned
            dn(ni-2,m) = a2*(ve(ni-1,m) - ve(ni-2,m)) + &
                         b2*(vn(ni-1,m) - vn(ni-3,m)) + &
                         c2*(vn(ni  ,m) - vn(ni-4,m))
            dn(ni-1,m) = a1*(ve(ni  ,m) - ve(ni-1,m)) + &
                         b1*(vn(ni  ,m) - vn(ni-2,m))
            !!dn(ni  ,m) = -(-25.0*vn(ni,m) + 48.0*ve(ni,m) - 36.0*vn(ni-1,m) + 16.0*ve(ni-1,m) - 3.0*vn(ni-2,m) ) / 6.0
            dn(ni  ,m) = -( -2.0*ve(ni+1,m) - 3.0*vn(ni,m) + 6.0*ve(ni,m) - vn(ni-1,m) )/3.0
            !!dn(ni  ,m) = (8.0*ve(ni  ,m) + 9.0*vn(ni,m) - 22.0*vn(ni-1,m) + 5.0*vn(ni-2,m) ) / 8.0
        end do
    end if

    do m=nst,ned
        do i=st,ed
            dn(i,m) = ap*(ve(i+1,m) - ve(i  ,m)) + &
                      a3*(vn(i+1,m) - vn(i-1,m)) + &
                      b3*(vn(i+2,m) - vn(i-2,m)) + &
                      c3*(vn(i+3,m) - vn(i-3,m))
        end do
    end do
       
end subroutine dn_via_ehen6ecc

subroutine dn_via_ehen8(ni,nst,ned,ngn,vn,nge,ve,nfs,nfe,dn)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_inter_scheme,nbc_intbc_scheme
    implicit none
    integer(kind_int),intent(in)    :: ni,nst,ned,ngn
    real(kind_real)  ,intent(in)    :: vn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)    :: nge
    real(kind_real)  ,intent(in)    :: ve(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)    :: nfs,nfe
    real(kind_real)  ,intent(inout) :: dn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int) :: i,m,st,ed
    real(kind_real)   :: a1,b1,a2,b2,c2,ap,a3,b3,c3

    if (ni < 4) then
        call dn_via_edge2(ni,nst,ned,ngn,vn,nge,ve,nfs,nfe,dn)
        return
    end if

    a1 = 4.D0/3.D0
    b1 = -1.D0/6.D0 

    a2 = 64.d0/45.d0 !!1.0d0 !!
    b2 = (16.D0 - 15.D0*a2)/24.D0
    c2 = (3.D0*a2 - 4.D0)/48.D0

    ap = 256.d0/175.d0 !!1.d0 !!
    a3 = (192.D0 - 175.D0*ap)/256.D0
    b3 = (35.D0*ap - 48.D0)/320.D0
    c3 = (64.D0 - 45.D0*ap)/3840.D0
    
    if (nfs == nbc_inter_scheme) then
        st = 1
    else if (nfs == nbc_intbc_scheme) then
    
        st = 4

        do m=nst,ned
            do i=1,3
                dn(i,m) = ( 2250.0*(ve(i  ,m) - ve(i-1,m)) - &
                             125.0*(ve(i+1,m) - ve(i-2,m)) + &
                               9.0*(ve(i+2,m) - ve(i-3,m)) ) / 1920.0
!!                dn(i,m) = ( 27.0*(ve(i  ,m) - ve(i-1,m)) - &
!!                                 (ve(i+1,m) - ve(i-2,m)) ) / 24.0
            end do
        end do
    
    else
        st = 4

        do m=nst,ned
            dn(3,m) = a2*(ve(3,m) - ve(2,m)) + &
                      b2*(vn(4,m) - vn(2,m)) + &
                      c2*(vn(5,m) - vn(1,m))
            dn(2,m) = a1*(ve(2,m) - ve(1,m)) + &
                      b1*(vn(3,m) - vn(1,m))
            dn(1,m) = (-25.0*vn(1,m) + 48.0*ve(1,m) - 36.0*vn(2,m) + 16.0*ve(2,m) - 3.0*vn(3,m) ) / 6.0
            !!dn(1,m) = (8.0*ve(1,m) - 15.0*vn(1,m) + 10.0*vn(2,m) - 3.0*vn(3,m) ) / 8.0
        end do
    end if

    if (nfe == nbc_inter_scheme) then
        ed = ni
    else if (nfe == nbc_intbc_scheme) then
    
        ed = ni - 3

        do m=nst,ned
            do i=ni-2,ni
                dn(i,m) = ( 2250.0*(ve(i  ,m) - ve(i-1,m)) - &
                             125.0*(ve(i+1,m) - ve(i-2,m)) + &
                               9.0*(ve(i+2,m) - ve(i-3,m)) ) / 1920.0
!!                dn(i,m) = ( 27.0*(ve(i  ,m) - ve(i-1,m)) - &
!!                             (ve(i+1,m) - ve(i-2,m)) ) / 24.0
            end do
        end do
    
    else
        ed = ni - 3

        do m=nst,ned
            dn(ni-2,m) = a2*(ve(ni-2,m) - ve(ni-3,m)) + &
                         b2*(vn(ni-1,m) - vn(ni-3,m)) + &
                         c2*(vn(ni  ,m) - vn(ni-4,m))
            dn(ni-1,m) = a1*(ve(ni-1,m) - ve(ni-2,m)) + &
                         b1*(vn(ni  ,m) - vn(ni-2,m))
            dn(ni  ,m) = -(-25.0*vn(ni,m) + 48.0*ve(ni-1,m) - 36.0*vn(ni-1,m) + 16.0*ve(ni-2,m) - 3.0*vn(ni-2,m) ) / 6.0
            !!dn(ni  ,m) = (8.0*ve(ni-1,m) + 9.0*vn(ni,m) - 22.0*vn(ni-1,m) + 5.0*vn(ni-2,m) ) / 8.0
        end do
    end if

    do m=nst,ned
        do i=st,ed
            dn(i,m) = ap*(ve(i  ,m) - ve(i-1,m)) + &
                      a3*(vn(i+1,m) - vn(i-1,m)) + &
                      b3*(vn(i+2,m) - vn(i-2,m)) + &
                      c3*(vn(i+3,m) - vn(i-3,m))
        end do
    end do
       
end subroutine dn_via_ehen8

subroutine dn_via_ehen8cc(ni,nst,ned,ngn,vn,nge,ve,nfs,nfe,dn)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_inter_scheme,nbc_intbc_scheme
    implicit none
    integer(kind_int),intent(in)    :: ni,nst,ned,ngn
    real(kind_real)  ,intent(in)    :: vn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)    :: nge
    real(kind_real)  ,intent(in)    :: ve(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)    :: nfs,nfe
    real(kind_real)  ,intent(inout) :: dn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int) :: i,m,st,ed
    real(kind_real)   :: a1,b1,a2,b2,c2,ap,a3,b3,c3

    if (ni < 4) then
        call dn_via_edge2cc(ni,nst,ned,ngn,vn,nge,ve,nfs,nfe,dn)
        return
    end if

    a1 = 4.D0/3.D0
    b1 = -1.D0/6.D0 

    a2 = 64.d0/45.d0 !!1.0d0 !!
    b2 = (16.D0 - 15.D0*a2)/24.D0
    c2 = (3.D0*a2 - 4.D0)/48.D0

    ap = 256.d0/175.d0 !!1.d0 !!
    a3 = (192.D0 - 175.D0*ap)/256.D0
    b3 = (35.D0*ap - 48.D0)/320.D0
    c3 = (64.D0 - 45.D0*ap)/3840.D0
    
    if (nfs == nbc_inter_scheme) then
        st = 1
    else if (nfs == nbc_intbc_scheme) then
    
        st = 4

        do m=nst,ned
            do i=1,3
                dn(i,m) = ( 2250.0*(ve(i+1,m) - ve(i  ,m)) - &
                             125.0*(ve(i+2,m) - ve(i-1,m)) + &
                               9.0*(ve(i+3,m) - ve(i-2,m)) ) / 1920.0
!!                dn(i,m) = ( 27.0*(ve(i+1,m) - ve(i  ,m)) - &
!!                                 (ve(i+2,m) - ve(i-1,m)) ) / 24.0
            end do
        end do
    
    else
        st = 4

        do m=nst,ned
            dn(3,m) = a2*(ve(4,m) - ve(3,m)) + &
                      b2*(vn(4,m) - vn(2,m)) + &
                      c2*(vn(5,m) - vn(1,m))
            dn(2,m) = a1*(ve(3,m) - ve(2,m)) + &
                      b1*(vn(3,m) - vn(1,m))
            !!dn(1,m) = (-25.0*vn(1,m) + 48.0*ve(2,m) - 36.0*vn(2,m) + 16.0*ve(3,m) - 3.0*vn(3,m) ) / 6.0
            dn(1,m) = ( -2.0*ve(1,m) - 3.0*vn(1,m) + 6.0*ve(2,m) - vn(2,m) )/3.0
            !!dn(1,m) = (8.0*ve(2,m) - 15.0*vn(1,m) + 10.0*vn(2,m) - 3.0*vn(3,m) ) / 8.0
            !!dn(1,m) = ( -23.0*ve(1,m) +   21.0*ve(2,m) +    3.0*ve(3,m) -      ve(4,m) ) / 24.0              ! third-order
            !!dn(1,m) = ve(2,m) - ve(1,m)
        end do
    end if

    if (nfe == nbc_inter_scheme) then
        ed = ni
    else if (nfe == nbc_intbc_scheme) then
    
        ed = ni - 3

        do m=nst,ned
            do i=ni-2,ni
                dn(i,m) = ( 2250.0*(ve(i+1,m) - ve(i  ,m)) - &
                             125.0*(ve(i+2,m) - ve(i-1,m)) + &
                               9.0*(ve(i+3,m) - ve(i-2,m)) ) / 1920.0
!!                dn(i,m) = ( 27.0*(ve(i+1,m) - ve(i  ,m)) - &
!!                             (ve(i+2,m) - ve(i-1,m)) ) / 24.0
            end do
        end do
    
    else
        ed = ni - 3

        do m=nst,ned
            dn(ni-2,m) = a2*(ve(ni-1,m) - ve(ni-2,m)) + &
                         b2*(vn(ni-1,m) - vn(ni-3,m)) + &
                         c2*(vn(ni  ,m) - vn(ni-4,m))
            dn(ni-1,m) = a1*(ve(ni  ,m) - ve(ni-1,m)) + &
                         b1*(vn(ni  ,m) - vn(ni-2,m))
            !!dn(ni  ,m) = -(-25.0*vn(ni,m) + 48.0*ve(ni,m) - 36.0*vn(ni-1,m) + 16.0*ve(ni-1,m) - 3.0*vn(ni-2,m) ) / 6.0
            dn(ni,m) = -( -2.0*ve(ni+1,m) - 3.0*vn(ni,m) + 6.0*ve(ni,m) - vn(ni-1,m) )/3.0
            !!dn(ni  ,m) = (8.0*ve(ni,m) + 9.0*vn(ni,m) - 22.0*vn(ni-1,m) + 5.0*vn(ni-2,m) ) / 8.0
            !!dn(ni,m)   = -( -23.0*ve(ni+1,m) +   21.0*ve(ni,m) +    3.0*ve(ni-1,m) -      ve(ni-2,m) ) / 24.0      ! third-order
            !!dn(ni,m) = ve(ni+1,m) - ve(ni,m)
        end do
    end if

    do m=nst,ned
        do i=st,ed
            dn(i,m) = ap*(ve(i+1,m) - ve(i  ,m)) + &
                      a3*(vn(i+1,m) - vn(i-1,m)) + &
                      b3*(vn(i+2,m) - vn(i-2,m)) + &
                      c3*(vn(i+3,m) - vn(i-3,m))
        end do
    end do
       
end subroutine dn_via_ehen8cc

subroutine dn_via_ehen8e(ni,nst,ned,ngn,vn,nge,ve,nfs,nfe,dn)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_inter_scheme
    use mod_variables, only : nghnode
    implicit none
    integer(kind_int),intent(in)    :: ni,nst,ned,ngn
    real(kind_real)  ,intent(in)    :: vn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)    :: nge
    real(kind_real)  ,intent(in)    :: ve(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)    :: nfs,nfe
    real(kind_real)  ,intent(inout) :: dn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int) :: i,m,st,ed
    real(kind_real)   :: a1,b1,a2,b2,c2,ap,a3,b3,c3

    if (ni < 4) then
        call dn_via_edge2(ni,nst,ned,ngn,vn,nge,ve,nfs,nfe,dn)
        return
    end if

    a1 = 4.D0/3.D0
    b1 = -1.D0/6.D0 

    a2 = 64.d0/45.d0 !!1.0d0 !!
    b2 = (16.D0 - 15.D0*a2)/24.D0
    c2 = (3.D0*a2 - 4.D0)/48.D0

    ap = 256.d0/175.d0 !!1.d0 !!
    a3 = (192.D0 - 175.D0*ap)/256.D0
    b3 = (35.D0*ap - 48.D0)/320.D0
    c3 = (64.D0 - 45.D0*ap)/3840.D0
    
    if (nfs == nbc_inter_scheme) then

        st = 1-ngn+nghnode
    
    else
        st = 4

        do m=nst,ned
            dn(3,m) = a2*(ve(3,m) - ve(2,m)) + &
                      b2*(vn(4,m) - vn(2,m)) + &
                      c2*(vn(5,m) - vn(1,m))
            dn(2,m) = a1*(ve(2,m) - ve(1,m)) + &
                      b1*(vn(3,m) - vn(1,m))
            dn(1,m) = (-25.0*vn(1,m) + 48.0*ve(1,m) - 36.0*vn(2,m) + 16.0*ve(2,m) - 3.0*vn(3,m) ) / 6.0
            !!dn(1,m) = (8.0*ve(1,m) - 15.0*vn(1,m) + 10.0*vn(2,m) - 3.0*vn(3,m) ) / 8.0
        end do
    end if

    if (nfe == nbc_inter_scheme) then

        ed = ni+ngn-nghnode
    
    else
        ed = ni - 3

        do m=nst,ned
            dn(ni-2,m) = a2*(ve(ni-2,m) - ve(ni-3,m)) + &
                         b2*(vn(ni-1,m) - vn(ni-3,m)) + &
                         c2*(vn(ni  ,m) - vn(ni-4,m))
            dn(ni-1,m) = a1*(ve(ni-1,m) - ve(ni-2,m)) + &
                         b1*(vn(ni  ,m) - vn(ni-2,m))
            dn(ni  ,m) = -(-25.0*vn(ni,m) + 48.0*ve(ni-1,m) - 36.0*vn(ni-1,m) + 16.0*ve(ni-2,m) - 3.0*vn(ni-2,m) ) / 6.0
            !!dn(ni  ,m) = (8.0*ve(ni-1,m) + 9.0*vn(ni,m) - 22.0*vn(ni-1,m) + 5.0*vn(ni-2,m) ) / 8.0
        end do
    end if

    do m=nst,ned
        do i=st,ed
            dn(i,m) = ap*(ve(i  ,m) - ve(i-1,m)) + &
                      a3*(vn(i+1,m) - vn(i-1,m)) + &
                      b3*(vn(i+2,m) - vn(i-2,m)) + &
                      c3*(vn(i+3,m) - vn(i-3,m))
        end do
    end do
       
end subroutine dn_via_ehen8e

subroutine dn_via_ehen8ecc(ni,nst,ned,ngn,vn,nge,ve,nfs,nfe,dn)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_inter_scheme
    use mod_variables, only : nghnode
    implicit none
    integer(kind_int),intent(in)    :: ni,nst,ned,ngn
    real(kind_real)  ,intent(in)    :: vn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)    :: nge
    real(kind_real)  ,intent(in)    :: ve(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)    :: nfs,nfe
    real(kind_real)  ,intent(inout) :: dn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int) :: i,m,st,ed
    real(kind_real)   :: a1,b1,a2,b2,c2,ap,a3,b3,c3

    if (ni < 4) then
        call dn_via_edge2cc(ni,nst,ned,ngn,vn,nge,ve,nfs,nfe,dn)
        return
    end if

    a1 = 4.D0/3.D0
    b1 = -1.D0/6.D0 

    a2 = 64.d0/45.d0 !!1.0d0 !!
    b2 = (16.D0 - 15.D0*a2)/24.D0
    c2 = (3.D0*a2 - 4.D0)/48.D0

    ap = 256.d0/175.d0 !!1.d0 !!
    a3 = (192.D0 - 175.D0*ap)/256.D0
    b3 = (35.D0*ap - 48.D0)/320.D0
    c3 = (64.D0 - 45.D0*ap)/3840.D0
    
    if (nfs == nbc_inter_scheme) then

        st = 1-ngn+nghnode
    
    else
        st = 4

        do m=nst,ned
            dn(3,m) = a2*(ve(4,m) - ve(3,m)) + &
                      b2*(vn(4,m) - vn(2,m)) + &
                      c2*(vn(5,m) - vn(1,m))
            dn(2,m) = a1*(ve(3,m) - ve(2,m)) + &
                      b1*(vn(3,m) - vn(1,m))
            !!dn(1,m) = (-25.0*vn(1,m) + 48.0*ve(2,m) - 36.0*vn(2,m) + 16.0*ve(3,m) - 3.0*vn(3,m) ) / 6.0
            dn(1,m) = ( -2.0*ve(1,m) - 3.0*vn(1,m) + 6.0*ve(2,m) - vn(2,m) )/3.0
            !!dn(1,m) = (8.0*ve(2,m) - 15.0*vn(1,m) + 10.0*vn(2,m) - 3.0*vn(3,m) ) / 8.0
        end do
    end if

    if (nfe == nbc_inter_scheme) then

        ed = ni+ngn-nghnode
    
    else
        ed = ni - 3

        do m=nst,ned
            dn(ni-2,m) = a2*(ve(ni-1,m) - ve(ni-2,m)) + &
                         b2*(vn(ni-1,m) - vn(ni-3,m)) + &
                         c2*(vn(ni  ,m) - vn(ni-4,m))
            dn(ni-1,m) = a1*(ve(ni  ,m) - ve(ni-1,m)) + &
                         b1*(vn(ni  ,m) - vn(ni-2,m))
            !!dn(ni  ,m) = -(-25.0*vn(ni,m) + 48.0*ve(ni,m) - 36.0*vn(ni-1,m) + 16.0*ve(ni-1,m) - 3.0*vn(ni-2,m) ) / 6.0
            dn(ni  ,m) = -( -2.0*ve(ni+1,m) - 3.0*vn(ni,m) + 6.0*ve(ni,m) - vn(ni-1,m) )/3.0
            !!dn(ni  ,m) = (8.0*ve(ni,m) + 9.0*vn(ni,m) - 22.0*vn(ni-1,m) + 5.0*vn(ni-2,m) ) / 8.0
        end do
    end if

    do m=nst,ned
        do i=st,ed
            dn(i,m) = ap*(ve(i+1,m) - ve(i  ,m)) + &
                      a3*(vn(i+1,m) - vn(i-1,m)) + &
                      b3*(vn(i+2,m) - vn(i-2,m)) + &
                      c3*(vn(i+3,m) - vn(i-3,m))
        end do
    end do
       
end subroutine dn_via_ehen8ecc

subroutine dn_via_ehcs6(ni,nst,ned,ngn,vn,nge,ve,nfs,nfe,dn)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_inter_scheme,nbc_intbc_scheme
    implicit none
    integer(kind_int),intent(in)    :: ni,nst,ned,ngn
    real(kind_real)  ,intent(in)    :: vn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)    :: nge
    real(kind_real)  ,intent(in)    :: ve(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)    :: nfs,nfe
    real(kind_real)  ,intent(inout) :: dn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int) :: i,m,st,ed
    real(kind_real)   :: a1,b1,a2,b2,c2

    if (ni < 4) then
        call dn_via_edge2(ni,nst,ned,ngn,vn,nge,ve,nfs,nfe,dn)
        return
    end if

    a1 = 4.D0/3.D0
    b1 = -1.D0/6.D0 

    a2 = 64.d0/45.d0
    b2 = (16.D0 - 15.D0*a2)/24.D0
    c2 = (3.D0*a2 - 4.D0)/48.D0

    
    if (nfs == nbc_inter_scheme) then
        st = 1
    else if (nfs == nbc_intbc_scheme) then
    
        st = 3

        do m=nst,ned
            do i=1,2
                dn(i,m) = ( 27.0*(ve(i  ,m) - ve(i-1,m)) - &
                                 (ve(i+1,m) - ve(i-2,m)) ) / 24.0
            end do
        end do
    
    else
        st = 3

        do m=nst,ned
            dn(2,m) = a1*(ve(2,m) - ve(1,m)) + &
                      b1*(vn(3,m) - vn(1,m))
            !!dn(1,m) = (-25.0*vn(1,m) + 48.0*ve(1,m) - 36.0*vn(2,m) + 16.0*ve(2,m) - 3.0*vn(3,m) ) / 6.0
            dn(1,m) = (8.0*ve(1,m) - 15.0*vn(1,m) + 10.0*vn(2,m) - 3.0*vn(3,m) ) / 8.0
        end do
    end if

    if (nfe == nbc_inter_scheme) then
        ed = ni
    else if (nfe == nbc_intbc_scheme) then
    
        ed = ni - 2

        do m=nst,ned
            do i=ni-1,ni
                dn(i,m) = ( 27.0*(ve(i  ,m) - ve(i-1,m)) - &
                             (ve(i+1,m) - ve(i-2,m)) ) / 24.0
            end do
        end do
    
    else
        ed = ni - 2

        do m=nst,ned
            dn(ni-1,m) = a1*(ve(ni-1,m) - ve(ni-2,m)) + &
                         b1*(vn(ni  ,m) - vn(ni-2,m))
            !!dn(ni  ,m) = -(-25.0*vn(ni,m) + 48.0*ve(ni-1,m) - 36.0*vn(ni-1,m) + 16.0*ve(ni-2,m) - 3.0*vn(ni-2,m) ) / 6.0
            dn(ni  ,m) = (8.0*ve(ni-1,m) + 9.0*vn(ni,m) - 22.0*vn(ni-1,m) + 5.0*vn(ni-2,m) ) / 8.0
        end do
    end if

    do m=nst,ned
        do i=st,ed
            dn(i,m) = a2*(ve(i  ,m) - ve(i-1,m)) + &
                      b2*(vn(i+1,m) - vn(i-1,m)) + &
                      c2*(vn(i+2,m) - vn(i-2,m))
        end do
    end do
       
end subroutine dn_via_ehcs6

subroutine dn_via_ehcs6cc(ni,nst,ned,ngn,vn,nge,ve,nfs,nfe,dn)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_inter_scheme,nbc_intbc_scheme
    implicit none
    integer(kind_int),intent(in)    :: ni,nst,ned,ngn
    real(kind_real)  ,intent(in)    :: vn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)    :: nge
    real(kind_real)  ,intent(in)    :: ve(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)    :: nfs,nfe
    real(kind_real)  ,intent(inout) :: dn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int) :: i,m,st,ed
    real(kind_real)   :: a1,b1,a2,b2,c2

    if (ni < 4) then
        call dn_via_edge2cc(ni,nst,ned,ngn,vn,nge,ve,nfs,nfe,dn)
        return
    end if

    a1 = 4.D0/3.D0
    b1 = -1.D0/6.D0 

    a2 = 64.d0/45.d0
    b2 = (16.D0 - 15.D0*a2)/24.D0
    c2 = (3.D0*a2 - 4.D0)/48.D0

    
    if (nfs == nbc_inter_scheme) then
        st = 1
    else if (nfs == nbc_intbc_scheme) then
    
        st = 3

        do m=nst,ned
            do i=1,2
                dn(i,m) = ( 27.0*(ve(i+1,m) - ve(i  ,m)) - &
                                 (ve(i+2,m) - ve(i-1,m)) ) / 24.0
            end do
        end do
    
    else
        st = 3

        do m=nst,ned
            dn(2,m) = a1*(ve(3,m) - ve(2,m)) + &
                      b1*(vn(3,m) - vn(1,m))
            !!dn(1,m) = (-25.0*vn(1,m) + 48.0*ve(2,m) - 36.0*vn(2,m) + 16.0*ve(3,m) - 3.0*vn(3,m) ) / 6.0
            dn(1,m) = ( -2.0*ve(1,m) - 3.0*vn(1,m) + 6.0*ve(2,m) - vn(2,m) )/3.0
            !!dn(1,m) = (8.0*ve(2,m) - 15.0*vn(1,m) + 10.0*vn(2,m) - 3.0*vn(3,m) ) / 8.0
        end do
    end if

    if (nfe == nbc_inter_scheme) then
        ed = ni
    else if (nfe == nbc_intbc_scheme) then
    
        ed = ni - 2

        do m=nst,ned
            do i=ni-1,ni
                dn(i,m) = ( 27.0*(ve(i+1,m) - ve(i  ,m)) - &
                                 (ve(i+2,m) - ve(i-1,m)) ) / 24.0
            end do
        end do
    
    else
        ed = ni - 2

        do m=nst,ned
            dn(ni-1,m) = a1*(ve(ni  ,m) - ve(ni-1,m)) + &
                         b1*(vn(ni  ,m) - vn(ni-2,m))
            !!dn(ni  ,m) = -(-25.0*vn(ni,m) + 48.0*ve(ni,m) - 36.0*vn(ni-1,m) + 16.0*ve(ni-1,m) - 3.0*vn(ni-2,m) ) / 6.0
            dn(ni,m) = -( -2.0*ve(ni+1,m) - 3.0*vn(ni,m) + 6.0*ve(ni,m) - vn(ni-1,m) )/3.0
            !!dn(ni  ,m) = (8.0*ve(ni,m) + 9.0*vn(ni,m) - 22.0*vn(ni-1,m) + 5.0*vn(ni-2,m) ) / 8.0
        end do
    end if

    do m=nst,ned
        do i=st,ed
            dn(i,m) = a2*(ve(i+1,m) - ve(i  ,m)) + &
                      b2*(vn(i+1,m) - vn(i-1,m)) + &
                      c2*(vn(i+2,m) - vn(i-2,m))
        end do
    end do
       
end subroutine dn_via_ehcs6cc

subroutine dn_via_ehcs6e(ni,nst,ned,ngn,vn,nge,ve,nfs,nfe,dn)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_inter_scheme
    use mod_variables, only : nghnode
    implicit none
    integer(kind_int),intent(in)    :: ni,nst,ned,ngn
    real(kind_real)  ,intent(in)    :: vn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)    :: nge
    real(kind_real)  ,intent(in)    :: ve(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)    :: nfs,nfe
    real(kind_real)  ,intent(inout) :: dn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int) :: i,m,st,ed
    real(kind_real)   :: a1,b1,a2,b2,c2

    if (ni < 4) then
        call dn_via_edge2(ni,nst,ned,ngn,vn,nge,ve,nfs,nfe,dn)
        return
    end if

    a1 = 4.D0/3.D0
    b1 = -1.D0/6.D0 

    a2 = 64.d0/45.d0
    b2 = (16.D0 - 15.D0*a2)/24.D0
    c2 = (3.D0*a2 - 4.D0)/48.D0

    
    if (nfs == nbc_inter_scheme) then

        st = 1-ngn+nghnode
    
    else
        st = 3

        do m=nst,ned
            dn(2,m) = a1*(ve(2,m) - ve(1,m)) + &
                      b1*(vn(3,m) - vn(1,m))
            !!dn(1,m) = (-25.0*vn(1,m) + 48.0*ve(1,m) - 36.0*vn(2,m) + 16.0*ve(2,m) - 3.0*vn(3,m) ) / 6.0
            dn(1,m) = (8.0*ve(1,m) - 15.0*vn(1,m) + 10.0*vn(2,m) - 3.0*vn(3,m) ) / 8.0
        end do
    end if

    if (nfe == nbc_inter_scheme) then

        ed = ni+ngn-nghnode
    
    else
        ed = ni - 2

        do m=nst,ned
            dn(ni-1,m) = a1*(ve(ni-1,m) - ve(ni-2,m)) + &
                         b1*(vn(ni  ,m) - vn(ni-2,m))
            !!dn(ni  ,m) = -(-25.0*vn(ni,m) + 48.0*ve(ni-1,m) - 36.0*vn(ni-1,m) + 16.0*ve(ni-2,m) - 3.0*vn(ni-2,m) ) / 6.0
            dn(ni  ,m) = (8.0*ve(ni-1,m) + 9.0*vn(ni,m) - 22.0*vn(ni-1,m) + 5.0*vn(ni-2,m) ) / 8.0
        end do
    end if

    do m=nst,ned
        do i=st,ed
            dn(i,m) = a2*(ve(i  ,m) - ve(i-1,m)) + &
                      b2*(vn(i+1,m) - vn(i-1,m)) + &
                      c2*(vn(i+2,m) - vn(i-2,m))
        end do
    end do
       
end subroutine dn_via_ehcs6e

subroutine dn_via_ehcs6ecc(ni,nst,ned,ngn,vn,nge,ve,nfs,nfe,dn)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_inter_scheme
    use mod_variables, only : nghnode
    implicit none
    integer(kind_int),intent(in)    :: ni,nst,ned,ngn
    real(kind_real)  ,intent(in)    :: vn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)    :: nge
    real(kind_real)  ,intent(in)    :: ve(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)    :: nfs,nfe
    real(kind_real)  ,intent(inout) :: dn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int) :: i,m,st,ed
    real(kind_real)   :: a1,b1,a2,b2,c2

    if (ni < 4) then
        call dn_via_edge2cc(ni,nst,ned,ngn,vn,nge,ve,nfs,nfe,dn)
        return
    end if

    a1 = 4.D0/3.D0
    b1 = -1.D0/6.D0 

    a2 = 64.d0/45.d0
    b2 = (16.D0 - 15.D0*a2)/24.D0
    c2 = (3.D0*a2 - 4.D0)/48.D0

    
    if (nfs == nbc_inter_scheme) then

        st = 1-ngn+nghnode
    
    else
        st = 3

        do m=nst,ned
            dn(2,m) = a1*(ve(3,m) - ve(2,m)) + &
                      b1*(vn(3,m) - vn(1,m))
            !!dn(1,m) = (-25.0*vn(1,m) + 48.0*ve(2,m) - 36.0*vn(2,m) + 16.0*ve(3,m) - 3.0*vn(3,m) ) / 6.0
            dn(1,m) = ( -2.0*ve(1,m) - 3.0*vn(1,m) + 6.0*ve(2,m) - vn(2,m) )/3.0
            !!dn(1,m) = (8.0*ve(2,m) - 15.0*vn(1,m) + 10.0*vn(2,m) - 3.0*vn(3,m) ) / 8.0
        end do
    end if

    if (nfe == nbc_inter_scheme) then

        ed = ni+ngn-nghnode
    
    else
        ed = ni - 2

        do m=nst,ned
            dn(ni-1,m) = a1*(ve(ni  ,m) - ve(ni-1,m)) + &
                         b1*(vn(ni  ,m) - vn(ni-2,m))
            !!dn(ni  ,m) = -(-25.0*vn(ni,m) + 48.0*ve(ni,m) - 36.0*vn(ni-1,m) + 16.0*ve(ni-1,m) - 3.0*vn(ni-2,m) ) / 6.0
            dn(ni,m) = -( -2.0*ve(ni+1,m) - 3.0*vn(ni,m) + 6.0*ve(ni,m) - vn(ni-1,m) )/3.0
            !!dn(ni  ,m) = (8.0*ve(ni,m) + 9.0*vn(ni,m) - 22.0*vn(ni-1,m) + 5.0*vn(ni-2,m) ) / 8.0
        end do
    end if

    do m=nst,ned
        do i=st,ed
            dn(i,m) = a2*(ve(i+1,m) - ve(i  ,m)) + &
                      b2*(vn(i+1,m) - vn(i-1,m)) + &
                      c2*(vn(i+2,m) - vn(i-2,m))
        end do
    end do
       
end subroutine dn_via_ehcs6ecc

subroutine dn_via_scsl4(ni,nst,ned,ngn,vn,nge,ve,nfs,nfe,dn)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_inter_scheme,nbc_intbc_scheme
    implicit none
    integer(kind_int),intent(in)    :: ni,nst,ned,ngn
    real(kind_real)  ,intent(in)    :: vn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)    :: nge
    real(kind_real)  ,intent(in)    :: ve(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)    :: nfs,nfe
    real(kind_real)  ,intent(inout) :: dn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)                 :: dn3(1-ngn:ni+ngn,nst:ned)
    integer(kind_int) :: i,m,st,ed
    real(kind_real)   :: a1,b1
    
    a1 = 4.D0/3.D0
    b1 = -1.D0/6.D0  
    
    if (ni < 6) then
        call dn_via_edge2(ni,nst,ned,ngn,vn,nge,ve,nfs,nfe,dn)
        return
    end if   
    
    call dn3_via_E2(ni,nst,ned,ngn,vn,nge,ve,nfs,nfe,dn3)
    
    if (nfs == nbc_inter_scheme) then
        st = 1
        
    else if (nfs == nbc_intbc_scheme) then
    
        st = 3

        do m=nst,ned
            do i=1,2
                dn(i,m) = ( 27.0*(ve(i  ,m) - ve(i-1,m)) - &
                                 (ve(i+1,m) - ve(i-2,m)) ) / 24.0
            end do
        end do
    
    else
        st = 3

        do m=nst,ned
            dn(2,m) = a1*(ve(2,m) - ve(1,m)) + &
                      b1*(vn(3,m) - vn(1,m))
            !!dn(1,m) = (-25.0*vn(1,m) + 48.0*ve(1,m) - 36.0*vn(2,m) + 16.0*ve(2,m) - 3.0*vn(3,m) ) / 6.0
            dn(1,m) = (8.0*ve(1,m) - 15.0*vn(1,m) + 10.0*vn(2,m) - 3.0*vn(3,m) ) / 8.0
        end do
    end if

    if (nfe == nbc_inter_scheme) then
        ed = ni
    else if (nfe == nbc_intbc_scheme) then
    
        ed = ni - 2

        do m=nst,ned
            do i=ni-1,ni
                dn(i,m) = ( 27.0*(ve(i  ,m) - ve(i-1,m)) - &
                             (ve(i+1,m) - ve(i-2,m)) ) / 24.0
            end do
        end do
    
    else
        ed = ni - 2

        do m=nst,ned
            dn(ni-1,m) = a1*(ve(ni-1,m) - ve(ni-2,m)) + &
                         b1*(vn(ni  ,m) - vn(ni-2,m))
            !!dn(ni  ,m) = -(-25.0*vn(ni,m) + 48.0*ve(ni-1,m) - 36.0*vn(ni-1,m) + 16.0*ve(ni-2,m) - 3.0*vn(ni-2,m) ) / 6.0
            dn(ni  ,m) = (8.0*ve(ni-1,m) + 9.0*vn(ni,m) - 22.0*vn(ni-1,m) + 5.0*vn(ni-2,m) ) / 8.0
        end do
    end if

    do m=nst,ned
        do i=st,ed
            dn(i,m) = (ve(i  ,m) - ve(i-1,m)) - &
                      1.0/24.0*dn3(i,m)
        end do
    end do    

end subroutine dn_via_scsl4

subroutine dn_via_scsl4cc(ni,nst,ned,ngn,vn,nge,ve,nfs,nfe,dn)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_inter_scheme,nbc_intbc_scheme
    implicit none
    integer(kind_int),intent(in)    :: ni,nst,ned,ngn
    real(kind_real)  ,intent(in)    :: vn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)    :: nge
    real(kind_real)  ,intent(in)    :: ve(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)    :: nfs,nfe
    real(kind_real)  ,intent(inout) :: dn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)                 :: dn3(1-ngn:ni+ngn,nst:ned)
    integer(kind_int) :: i,m,st,ed
    real(kind_real)   :: a1,b1
    
    a1 = 4.D0/3.D0
    b1 = -1.D0/6.D0  
    
    if (ni < 6) then
        call dn_via_edge2cc(ni,nst,ned,ngn,vn,nge,ve,nfs,nfe,dn)
        return
    end if   
    
    call dn3_via_E2(ni,nst,ned,ngn,vn,nge,ve,nfs,nfe,dn3)
    
    if (nfs == nbc_inter_scheme) then
        st = 1
        
    else if (nfs == nbc_intbc_scheme) then
    
        st = 3

        do m=nst,ned
            do i=1,2
                dn(i,m) = ( 27.0*(ve(i+1,m) - ve(i  ,m)) - &
                                 (ve(i+2,m) - ve(i-1,m)) ) / 24.0
            end do
        end do
    
    else
        st = 3

        do m=nst,ned
            dn(2,m) = a1*(ve(3,m) - ve(2,m)) + &
                      b1*(vn(3,m) - vn(1,m))
            !!dn(1,m) = (-25.0*vn(1,m) + 48.0*ve(2,m) - 36.0*vn(2,m) + 16.0*ve(3,m) - 3.0*vn(3,m) ) / 6.0
            dn(1,m) = ( -2.0*ve(1,m) - 3.0*vn(1,m) + 6.0*ve(2,m) - vn(2,m) )/3.0
            !!dn(1,m) = (8.0*ve(2,m) - 15.0*vn(1,m) + 10.0*vn(2,m) - 3.0*vn(3,m) ) / 8.0
        end do
    end if

    if (nfe == nbc_inter_scheme) then
        ed = ni
    else if (nfe == nbc_intbc_scheme) then
    
        ed = ni - 2

        do m=nst,ned
            do i=ni-1,ni
                dn(i,m) = ( 27.0*(ve(i+1,m) - ve(i  ,m)) - &
                                 (ve(i+2,m) - ve(i-1,m)) ) / 24.0
            end do
        end do
    
    else
        ed = ni - 2

        do m=nst,ned
            dn(ni-1,m) = a1*(ve(ni  ,m) - ve(ni-1,m)) + &
                         b1*(vn(ni  ,m) - vn(ni-2,m))
            !!dn(ni  ,m) = -(-25.0*vn(ni,m) + 48.0*ve(ni,m) - 36.0*vn(ni-1,m) + 16.0*ve(ni-1,m) - 3.0*vn(ni-2,m) ) / 6.0
            dn(ni  ,m) = -( -2.0*ve(ni+1,m) - 3.0*vn(ni,m) + 6.0*ve(ni,m) - vn(ni-1,m) )/3.0
            !!dn(ni  ,m) = (8.0*ve(ni,m) + 9.0*vn(ni,m) - 22.0*vn(ni-1,m) + 5.0*vn(ni-2,m) ) / 8.0
        end do
    end if

    do m=nst,ned
        do i=st,ed
            dn(i,m) = (ve(i+1,m) - ve(i  ,m)) - &
                      1.0/24.0*dn3(i,m)
        end do
    end do    

end subroutine dn_via_scsl4cc

subroutine dn3_via_E2(ni,nst,ned,ngn,vn,nge,ve,nfs,nfe,dn3)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_inter_scheme,nbc_intbc_scheme
    implicit none
    integer(kind_int),intent(in)    :: ni,nst,ned,ngn
    real(kind_real)  ,intent(in)    :: vn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)    :: nge
    real(kind_real)  ,intent(in)    :: ve(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)    :: nfs,nfe
    real(kind_real)  ,intent(out)   :: dn3(1-ngn:ni+ngn,nst:ned)
    integer(kind_int) :: i,m,st,ed
    
    !>
    !!������׾��Ƚڵ����׿ռ䵼��
    if (nfs == nbc_inter_scheme) then
        
        st=1      
                                                     
    else
        
        st=3      

    end if
    !>

    !>
    if (nfe == nbc_inter_scheme) then

        ed=ni
        
    else

        ed=ni-2

    end if
    !>

    !>
    !!������׾��Ƚڵ����׿ռ䵼���ڵ�
    do m=nst,ned
        do i=st,ed
            dn3(i,m) =    -1.0*( vn(i+1,m) - vn(i-1,m) )          +  &
                           0.5*( vn(i+2,m) - vn(i-2,m) )     
        end do
    end do    

endsubroutine dn3_via_E2

subroutine dn_via_scsl4e(ni,nst,ned,ngn,vn,nge,ve,nfs,nfe,dn)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_inter_scheme
    use mod_variables, only : nghnode
    implicit none
    integer(kind_int),intent(in)    :: ni,nst,ned,ngn
    real(kind_real)  ,intent(in)    :: vn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)    :: nge
    real(kind_real)  ,intent(in)    :: ve(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)    :: nfs,nfe
    real(kind_real)  ,intent(inout) :: dn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)                 :: dn3(1-ngn:ni+ngn,nst:ned)
    integer(kind_int) :: i,m,st,ed
    real(kind_real)   :: a1,b1
    
    a1 = 4.D0/3.D0
    b1 = -1.D0/6.D0  
    
    if (ni < 6) then
        call dn_via_edge2(ni,nst,ned,ngn,vn,nge,ve,nfs,nfe,dn)
        return
    end if   
    
    call dn3_via_E2e(ni,nst,ned,ngn,vn,nge,ve,nfs,nfe,dn3)
    
    if (nfs == nbc_inter_scheme) then

        st = 1-ngn+nghnode
        
    else
        st = 3

        do m=nst,ned
            dn(2,m) = a1*(ve(2,m) - ve(1,m)) + &
                      b1*(vn(3,m) - vn(1,m))
            !!dn(1,m) = (-25.0*vn(1,m) + 48.0*ve(1,m) - 36.0*vn(2,m) + 16.0*ve(2,m) - 3.0*vn(3,m) ) / 6.0
            dn(1,m) = (8.0*ve(1,m) - 15.0*vn(1,m) + 10.0*vn(2,m) - 3.0*vn(3,m) ) / 8.0
        end do
    end if

    if (nfe == nbc_inter_scheme) then

        ed = ni+ngn-nghnode
    
    else
        ed = ni - 2

        do m=nst,ned
            dn(ni-1,m) = a1*(ve(ni-1,m) - ve(ni-2,m)) + &
                         b1*(vn(ni  ,m) - vn(ni-2,m))
            !!dn(ni  ,m) = -(-25.0*vn(ni,m) + 48.0*ve(ni-1,m) - 36.0*vn(ni-1,m) + 16.0*ve(ni-2,m) - 3.0*vn(ni-2,m) ) / 6.0
            dn(ni  ,m) = (8.0*ve(ni-1,m) + 9.0*vn(ni,m) - 22.0*vn(ni-1,m) + 5.0*vn(ni-2,m) ) / 8.0
        end do
    end if

    do m=nst,ned
        do i=st,ed
            dn(i,m) = (ve(i  ,m) - ve(i-1,m)) - &
                      1.0/24.0*dn3(i,m)
        end do
    end do    

end subroutine dn_via_scsl4e

subroutine dn_via_scsl4ecc(ni,nst,ned,ngn,vn,nge,ve,nfs,nfe,dn)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_inter_scheme
    use mod_variables, only : nghnode
    implicit none
    integer(kind_int),intent(in)    :: ni,nst,ned,ngn
    real(kind_real)  ,intent(in)    :: vn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)    :: nge
    real(kind_real)  ,intent(in)    :: ve(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)    :: nfs,nfe
    real(kind_real)  ,intent(inout) :: dn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)                 :: dn3(1-ngn:ni+ngn,nst:ned)
    integer(kind_int) :: i,m,st,ed
    real(kind_real)   :: a1,b1
    
    a1 = 4.D0/3.D0
    b1 = -1.D0/6.D0  
    
    if (ni < 6) then
        call dn_via_edge2cc(ni,nst,ned,ngn,vn,nge,ve,nfs,nfe,dn)
        return
    end if   
    
    call dn3_via_E2e(ni,nst,ned,ngn,vn,nge,ve,nfs,nfe,dn3)
    
    if (nfs == nbc_inter_scheme) then

        st = 1-ngn+nghnode
        
    else
        st = 3

        do m=nst,ned
            dn(2,m) = a1*(ve(3,m) - ve(2,m)) + &
                      b1*(vn(3,m) - vn(1,m))
            !!dn(1,m) = (-25.0*vn(1,m) + 48.0*ve(2,m) - 36.0*vn(2,m) + 16.0*ve(3,m) - 3.0*vn(3,m) ) / 6.0
            dn(1,m) = ( -2.0*ve(1,m) - 3.0*vn(1,m) + 6.0*ve(2,m) - vn(2,m) )/3.0
            !!dn(1,m) = (8.0*ve(2,m) - 15.0*vn(1,m) + 10.0*vn(2,m) - 3.0*vn(3,m) ) / 8.0
        end do
    end if

    if (nfe == nbc_inter_scheme) then

        ed = ni+ngn-nghnode
    
    else
        ed = ni - 2

        do m=nst,ned
            dn(ni-1,m) = a1*(ve(ni  ,m) - ve(ni-1,m)) + &
                         b1*(vn(ni  ,m) - vn(ni-2,m))
            !!dn(ni  ,m) = -(-25.0*vn(ni,m) + 48.0*ve(ni,m) - 36.0*vn(ni-1,m) + 16.0*ve(ni-1,m) - 3.0*vn(ni-2,m) ) / 6.0
            dn(ni,m) = -( -2.0*ve(ni+1,m) - 3.0*vn(ni,m) + 6.0*ve(ni,m) - vn(ni-1,m) )/3.0
            !!dn(ni  ,m) = (8.0*ve(ni,m) + 9.0*vn(ni,m) - 22.0*vn(ni-1,m) + 5.0*vn(ni-2,m) ) / 8.0
        end do
    end if

    do m=nst,ned
        do i=st,ed
            dn(i,m) = (ve(i+1,m) - ve(i  ,m)) - &
                      1.0/24.0*dn3(i,m)
        end do
    end do    

end subroutine dn_via_scsl4ecc

subroutine dn3_via_E2e(ni,nst,ned,ngn,vn,nge,ve,nfs,nfe,dn3)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_inter_scheme
    use mod_variables, only : nghnode
    implicit none
    integer(kind_int),intent(in)    :: ni,nst,ned,ngn
    real(kind_real)  ,intent(in)    :: vn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)    :: nge
    real(kind_real)  ,intent(in)    :: ve(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)    :: nfs,nfe
    real(kind_real)  ,intent(out)   :: dn3(1-ngn:ni+ngn,nst:ned)
    integer(kind_int) :: i,m,st,ed
    
    !>
    !!������׾��Ƚڵ����׿ռ䵼��
    if (nfs == nbc_inter_scheme) then
        
        st=1-ngn+nghnode
                                                     
    else
        
        st=3      

    end if
    !>

    !>
    if (nfe == nbc_inter_scheme) then

        ed=ni+ngn-nghnode
        
    else

        ed=ni-2

    end if
    !>

    !>
    !!������׾��Ƚڵ����׿ռ䵼���ڵ�
    do m=nst,ned
        do i=st,ed
            dn3(i,m) =    -1.0*( vn(i+1,m) - vn(i-1,m) )          +  &
                           0.5*( vn(i+2,m) - vn(i-2,m) )     
        end do
    end do    

endsubroutine dn3_via_E2e

subroutine dn_via_scsl6(ni,nst,ned,ngn,vn,nge,ve,nfs,nfe,dn)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_inter_scheme,nbc_intbc_scheme
    implicit none
    integer(kind_int),intent(in)    :: ni,nst,ned,ngn
    real(kind_real)  ,intent(in)    :: vn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)    :: nge
    real(kind_real)  ,intent(in)    :: ve(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)    :: nfs,nfe
    real(kind_real)  ,intent(inout) :: dn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)                 :: dn3(1-ngn:ni+ngn,nst:ned)
    real(kind_real)                 :: dn5(1-ngn:ni+ngn,nst:ned)
    integer(kind_int) :: i,m,st,ed
    real(kind_real)   :: a1,b1,a2,b2,c2
    
    a1 = 4.D0/3.D0
    b1 = -1.D0/6.D0
    
    a2 = 64.d0/45.d0 !!1.0d0 !!
    b2 = (16.D0 - 15.D0*a2)/24.D0
    c2 = (3.D0*a2 - 4.D0)/48.D0    
    
    if (ni < 8) then
        call dn_via_edge2(ni,nst,ned,ngn,vn,nge,ve,nfs,nfe,dn)
        return
    end if   
    
    call dn3_via_E4(ni,nst,ned,ngn,vn,nge,ve,nfs,nfe,dn3)
    
    call dn5_via_E2(ni,nst,ned,ngn,vn,nge,ve,nfs,nfe,dn5)
    
    if (nfs == nbc_inter_scheme) then
        
        st = 1
        
    else if (nfs == nbc_intbc_scheme) then
    
        st = 4

        do m=nst,ned
            do i=1,3
                
                dn(i,m) = ( 2250.0*(ve(i  ,m) - ve(i-1,m)) - &
                             125.0*(ve(i+1,m) - ve(i-2,m)) + &
                               9.0*(ve(i+2,m) - ve(i-3,m)) ) / 1920.0
                
            end do
        end do
    
    else
        st = 4

        do m=nst,ned
            
            dn(3,m) = a2*(ve(3,m) - ve(2,m)) + &
                      b2*(vn(4,m) - vn(2,m)) + &
                      c2*(vn(5,m) - vn(1,m))            
            dn(2,m) = a1*(ve(2,m) - ve(1,m)) + &
                      b1*(vn(3,m) - vn(1,m))
            !!dn(1,m) = (-25.0*vn(1,m) + 48.0*ve(1,m) - 36.0*vn(2,m) + 16.0*ve(2,m) - 3.0*vn(3,m) ) / 6.0
            dn(1,m) = (8.0*ve(1,m) - 15.0*vn(1,m) + 10.0*vn(2,m) - 3.0*vn(3,m) ) / 8.0
            
        end do
    end if

    if (nfe == nbc_inter_scheme) then
        
        ed = ni
        
    else if (nfe == nbc_intbc_scheme) then
    
        ed = ni - 3

        do m=nst,ned
            do i=ni-2,ni
                
                dn(i,m) = ( 2250.0*(ve(i  ,m) - ve(i-1,m)) - &
                             125.0*(ve(i+1,m) - ve(i-2,m)) + &
                               9.0*(ve(i+2,m) - ve(i-3,m)) ) / 1920.0
                
            end do
        end do
    
    else
        ed = ni - 3

        do m=nst,ned
            
            dn(ni-2,m) = a2*(ve(ni-2,m) - ve(ni-3,m)) + &
                         b2*(vn(ni-1,m) - vn(ni-3,m)) + &
                         c2*(vn(ni  ,m) - vn(ni-4,m))            
            dn(ni-1,m) = a1*(ve(ni-1,m) - ve(ni-2,m)) + &
                         b1*(vn(ni  ,m) - vn(ni-2,m))
            !!dn(ni  ,m) = -(-25.0*vn(ni,m) + 48.0*ve(ni-1,m) - 36.0*vn(ni-1,m) + 16.0*ve(ni-2,m) - 3.0*vn(ni-2,m) ) / 6.0
            dn(ni  ,m) = (8.0*ve(ni-1,m) + 9.0*vn(ni,m) - 22.0*vn(ni-1,m) + 5.0*vn(ni-2,m) ) / 8.0
            
        end do
    end if

    do m=nst,ned
        do i=st,ed
            dn(i,m) = (ve(i  ,m) - ve(i-1,m)) - &
                            1.0/24.0*dn3(i,m) - &
                          1.0/1920.0*dn5(i,m)
        end do
    end do    

end subroutine dn_via_scsl6

subroutine dn_via_scsl6cc(ni,nst,ned,ngn,vn,nge,ve,nfs,nfe,dn)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_inter_scheme,nbc_intbc_scheme
    implicit none
    integer(kind_int),intent(in)    :: ni,nst,ned,ngn
    real(kind_real)  ,intent(in)    :: vn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)    :: nge
    real(kind_real)  ,intent(in)    :: ve(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)    :: nfs,nfe
    real(kind_real)  ,intent(inout) :: dn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)                 :: dn3(1-ngn:ni+ngn,nst:ned)
    real(kind_real)                 :: dn5(1-ngn:ni+ngn,nst:ned)
    integer(kind_int) :: i,m,st,ed
    real(kind_real)   :: a1,b1,a2,b2,c2
    
    a1 = 4.D0/3.D0
    b1 = -1.D0/6.D0
    
    a2 = 64.d0/45.d0 !!1.0d0 !!
    b2 = (16.D0 - 15.D0*a2)/24.D0
    c2 = (3.D0*a2 - 4.D0)/48.D0    
    
    if (ni < 8) then
        call dn_via_edge2cc(ni,nst,ned,ngn,vn,nge,ve,nfs,nfe,dn)
        return
    end if   
    
    call dn3_via_E4(ni,nst,ned,ngn,vn,nge,ve,nfs,nfe,dn3)
    
    call dn5_via_E2(ni,nst,ned,ngn,vn,nge,ve,nfs,nfe,dn5)
    
    if (nfs == nbc_inter_scheme) then
        
        st = 1
        
    else if (nfs == nbc_intbc_scheme) then
    
        st = 4

        do m=nst,ned
            do i=1,3
                
                dn(i,m) = ( 2250.0*(ve(i+1,m) - ve(i  ,m)) - &
                             125.0*(ve(i+2,m) - ve(i-1,m)) + &
                               9.0*(ve(i+3,m) - ve(i-2,m)) ) / 1920.0
                
            end do
        end do
    
    else
        st = 4

        do m=nst,ned
            
            dn(3,m) = a2*(ve(4,m) - ve(3,m)) + &
                      b2*(vn(4,m) - vn(2,m)) + &
                      c2*(vn(5,m) - vn(1,m))            
            dn(2,m) = a1*(ve(3,m) - ve(2,m)) + &
                      b1*(vn(3,m) - vn(1,m))
            !!dn(1,m) = (-25.0*vn(1,m) + 48.0*ve(2,m) - 36.0*vn(2,m) + 16.0*ve(3,m) - 3.0*vn(3,m) ) / 6.0
            dn(1,m) = ( -2.0*ve(1,m) - 3.0*vn(1,m) + 6.0*ve(2,m) - vn(2,m) )/3.0
            !!dn(1,m) = (8.0*ve(2,m) - 15.0*vn(1,m) + 10.0*vn(2,m) - 3.0*vn(3,m) ) / 8.0
            
        end do
    end if

    if (nfe == nbc_inter_scheme) then
        
        ed = ni
        
    else if (nfe == nbc_intbc_scheme) then
    
        ed = ni - 3

        do m=nst,ned
            do i=ni-2,ni
                
                dn(i,m) = ( 2250.0*(ve(i+1,m) - ve(i  ,m)) - &
                             125.0*(ve(i+2,m) - ve(i-1,m)) + &
                               9.0*(ve(i+3,m) - ve(i-2,m)) ) / 1920.0
                
            end do
        end do
    
    else
        ed = ni - 3

        do m=nst,ned
            
            dn(ni-2,m) = a2*(ve(ni-1,m) - ve(ni-2,m)) + &
                         b2*(vn(ni-1,m) - vn(ni-3,m)) + &
                         c2*(vn(ni  ,m) - vn(ni-4,m))            
            dn(ni-1,m) = a1*(ve(ni  ,m) - ve(ni-1,m)) + &
                         b1*(vn(ni  ,m) - vn(ni-2,m))
            !!dn(ni  ,m) = -(-25.0*vn(ni,m) + 48.0*ve(ni,m) - 36.0*vn(ni-1,m) + 16.0*ve(ni-1,m) - 3.0*vn(ni-2,m) ) / 6.0
            dn(ni,m) = -( -2.0*ve(ni+1,m) - 3.0*vn(ni,m) + 6.0*ve(ni,m) - vn(ni-1,m) )/3.0
            !!dn(ni  ,m) = (8.0*ve(ni,m) + 9.0*vn(ni,m) - 22.0*vn(ni-1,m) + 5.0*vn(ni-2,m) ) / 8.0
            
        end do
    end if

    do m=nst,ned
        do i=st,ed
            dn(i,m) = (ve(i+1,m) - ve(i  ,m)) - &
                            1.0/24.0*dn3(i,m) - &
                          1.0/1920.0*dn5(i,m)
        end do
    end do    

end subroutine dn_via_scsl6cc

subroutine dn_via_scsl6e(ni,nst,ned,ngn,vn,nge,ve,nfs,nfe,dn)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_inter_scheme
    use mod_variables, only : nghnode
    implicit none
    integer(kind_int),intent(in)    :: ni,nst,ned,ngn
    real(kind_real)  ,intent(in)    :: vn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)    :: nge
    real(kind_real)  ,intent(in)    :: ve(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)    :: nfs,nfe
    real(kind_real)  ,intent(inout) :: dn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)                 :: dn3(1-ngn:ni+ngn,nst:ned)
    real(kind_real)                 :: dn5(1-ngn:ni+ngn,nst:ned)
    integer(kind_int) :: i,m,st,ed
    real(kind_real)   :: a1,b1,a2,b2,c2
    
    a1 = 4.D0/3.D0
    b1 = -1.D0/6.D0
    
    a2 = 64.d0/45.d0 !!1.0d0 !!
    b2 = (16.D0 - 15.D0*a2)/24.D0
    c2 = (3.D0*a2 - 4.D0)/48.D0    
    
    if (ni < 8) then
        call dn_via_edge2(ni,nst,ned,ngn,vn,nge,ve,nfs,nfe,dn)
        return
    end if   
    
    call dn3_via_E4e(ni,nst,ned,ngn,vn,nge,ve,nfs,nfe,dn3)
    
    call dn5_via_E2e(ni,nst,ned,ngn,vn,nge,ve,nfs,nfe,dn5)
    
    if (nfs == nbc_inter_scheme) then
        
        st = 1-ngn+nghnode
    
    else
        st = 4

        do m=nst,ned
            
            dn(3,m) = a2*(ve(3,m) - ve(2,m)) + &
                      b2*(vn(4,m) - vn(2,m)) + &
                      c2*(vn(5,m) - vn(1,m))            
            dn(2,m) = a1*(ve(2,m) - ve(1,m)) + &
                      b1*(vn(3,m) - vn(1,m))
            !!dn(1,m) = (-25.0*vn(1,m) + 48.0*ve(1,m) - 36.0*vn(2,m) + 16.0*ve(2,m) - 3.0*vn(3,m) ) / 6.0
            dn(1,m) = (8.0*ve(1,m) - 15.0*vn(1,m) + 10.0*vn(2,m) - 3.0*vn(3,m) ) / 8.0
            
        end do
    end if

    if (nfe == nbc_inter_scheme) then
        
        ed = ni+ngn-nghnode
    
    else
        ed = ni - 3

        do m=nst,ned
            
            dn(ni-2,m) = a2*(ve(ni-2,m) - ve(ni-3,m)) + &
                         b2*(vn(ni-1,m) - vn(ni-3,m)) + &
                         c2*(vn(ni  ,m) - vn(ni-4,m))            
            dn(ni-1,m) = a1*(ve(ni-1,m) - ve(ni-2,m)) + &
                         b1*(vn(ni  ,m) - vn(ni-2,m))
            !!dn(ni  ,m) = -(-25.0*vn(ni,m) + 48.0*ve(ni-1,m) - 36.0*vn(ni-1,m) + 16.0*ve(ni-2,m) - 3.0*vn(ni-2,m) ) / 6.0
            dn(ni  ,m) = (8.0*ve(ni-1,m) + 9.0*vn(ni,m) - 22.0*vn(ni-1,m) + 5.0*vn(ni-2,m) ) / 8.0
            
        end do
    end if

    do m=nst,ned
        do i=st,ed
            dn(i,m) = (ve(i  ,m) - ve(i-1,m)) - &
                            1.0/24.0*dn3(i,m) - &
                          1.0/1920.0*dn5(i,m)
        end do
    end do    

end subroutine dn_via_scsl6e

subroutine dn_via_scsl6ecc(ni,nst,ned,ngn,vn,nge,ve,nfs,nfe,dn)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_inter_scheme
    use mod_variables, only : nghnode
    implicit none
    integer(kind_int),intent(in)    :: ni,nst,ned,ngn
    real(kind_real)  ,intent(in)    :: vn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)    :: nge
    real(kind_real)  ,intent(in)    :: ve(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)    :: nfs,nfe
    real(kind_real)  ,intent(inout) :: dn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)                 :: dn3(1-ngn:ni+ngn,nst:ned)
    real(kind_real)                 :: dn5(1-ngn:ni+ngn,nst:ned)
    integer(kind_int) :: i,m,st,ed
    real(kind_real)   :: a1,b1,a2,b2,c2
    
    a1 = 4.D0/3.D0
    b1 = -1.D0/6.D0
    
    a2 = 64.d0/45.d0 !!1.0d0 !!
    b2 = (16.D0 - 15.D0*a2)/24.D0
    c2 = (3.D0*a2 - 4.D0)/48.D0    
    
    if (ni < 8) then
        call dn_via_edge2cc(ni,nst,ned,ngn,vn,nge,ve,nfs,nfe,dn)
        return
    end if   
    
    call dn3_via_E4e(ni,nst,ned,ngn,vn,nge,ve,nfs,nfe,dn3)
    
    call dn5_via_E2e(ni,nst,ned,ngn,vn,nge,ve,nfs,nfe,dn5)
    
    if (nfs == nbc_inter_scheme) then
        
        st = 1-ngn+nghnode
    
    else
        st = 4

        do m=nst,ned
            
            dn(3,m) = a2*(ve(4,m) - ve(3,m)) + &
                      b2*(vn(4,m) - vn(2,m)) + &
                      c2*(vn(5,m) - vn(1,m))            
            dn(2,m) = a1*(ve(3,m) - ve(2,m)) + &
                      b1*(vn(3,m) - vn(1,m))
            !!dn(1,m) = (-25.0*vn(1,m) + 48.0*ve(2,m) - 36.0*vn(2,m) + 16.0*ve(3,m) - 3.0*vn(3,m) ) / 6.0
            dn(1,m) = ( -2.0*ve(1,m) - 3.0*vn(1,m) + 6.0*ve(2,m) - vn(2,m) )/3.0
            !!dn(1,m) = (8.0*ve(2,m) - 15.0*vn(1,m) + 10.0*vn(2,m) - 3.0*vn(3,m) ) / 8.0
            
        end do
    end if

    if (nfe == nbc_inter_scheme) then
        
        ed = ni+ngn-nghnode
    
    else
        ed = ni - 3

        do m=nst,ned
            
            dn(ni-2,m) = a2*(ve(ni-1,m) - ve(ni-2,m)) + &
                         b2*(vn(ni-1,m) - vn(ni-3,m)) + &
                         c2*(vn(ni  ,m) - vn(ni-4,m))            
            dn(ni-1,m) = a1*(ve(ni  ,m) - ve(ni-1,m)) + &
                         b1*(vn(ni  ,m) - vn(ni-2,m))
            !!dn(ni  ,m) = -(-25.0*vn(ni,m) + 48.0*ve(ni,m) - 36.0*vn(ni-1,m) + 16.0*ve(ni-1,m) - 3.0*vn(ni-2,m) ) / 6.0
            dn(ni  ,m) = -( -2.0*ve(ni+1,m) - 3.0*vn(ni,m) + 6.0*ve(ni,m) - vn(ni-1,m) )/3.0
            !!dn(ni  ,m) = (8.0*ve(ni,m) + 9.0*vn(ni,m) - 22.0*vn(ni-1,m) + 5.0*vn(ni-2,m) ) / 8.0
            
        end do
    end if

    do m=nst,ned
        do i=st,ed
            dn(i,m) = (ve(i+1,m) - ve(i  ,m)) - &
                            1.0/24.0*dn3(i,m) - &
                          1.0/1920.0*dn5(i,m)
        end do
    end do    

end subroutine dn_via_scsl6ecc

subroutine dn3_via_E4(ni,nst,ned,ngn,vn,nge,ve,nfs,nfe,dn3)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_inter_scheme
    implicit none
    integer(kind_int),intent(in)    :: ni,nst,ned,ngn
    real(kind_real)  ,intent(in)    :: vn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)    :: nge
    real(kind_real)  ,intent(in)    :: ve(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)    :: nfs,nfe
    real(kind_real)  ,intent(out)   :: dn3(1-ngn:ni+ngn,nst:ned)
    integer(kind_int) :: i,m,st,ed
    
    !>
    !!������׾��Ƚڵ����׿ռ䵼��
    if (nfs == nbc_inter_scheme) then
        
        st=1      
                                                     
    else
        
        st=4      

    end if
    !>

    !>
    if (nfe == nbc_inter_scheme) then

        ed=ni
        
    else

        ed=ni-3

    end if
    !>

    !>
    !!������׾��Ƚڵ����׿ռ䵼���ڵ�
    do m=nst,ned
        do i=st,ed
            dn3(i,m) =    -13.0/8.0*( vn(i+1,m) - vn(i-1,m) )          +  &
                                    ( vn(i+2,m) - vn(i-2,m) )          -  &
                            1.0/8.0*( vn(i+3,m) - vn(i-3,m) )
        end do
    end do    

endsubroutine dn3_via_E4

subroutine dn3_via_E4e(ni,nst,ned,ngn,vn,nge,ve,nfs,nfe,dn3)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_inter_scheme
    use mod_variables, only : nghnode
    implicit none
    integer(kind_int),intent(in)    :: ni,nst,ned,ngn
    real(kind_real)  ,intent(in)    :: vn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)    :: nge
    real(kind_real)  ,intent(in)    :: ve(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)    :: nfs,nfe
    real(kind_real)  ,intent(out)   :: dn3(1-ngn:ni+ngn,nst:ned)
    integer(kind_int) :: i,m,st,ed
    
    !>
    !!������׾��Ƚڵ����׿ռ䵼��
    if (nfs == nbc_inter_scheme) then
        
        st=1-ngn+nghnode
                                                     
    else
        
        st=4      

    end if
    !>

    !>
    if (nfe == nbc_inter_scheme) then

        ed=ni+ngn-nghnode
        
    else

        ed=ni-3

    end if
    !>

    !>
    !!������׾��Ƚڵ����׿ռ䵼���ڵ�
    do m=nst,ned
        do i=st,ed
            dn3(i,m) =    -13.0/8.0*( vn(i+1,m) - vn(i-1,m) )          +  &
                                    ( vn(i+2,m) - vn(i-2,m) )          -  &
                            1.0/8.0*( vn(i+3,m) - vn(i-3,m) )
        end do
    end do    

endsubroutine dn3_via_E4e

subroutine dn5_via_E2(ni,nst,ned,ngn,vn,nge,ve,nfs,nfe,dn5)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_inter_scheme
    implicit none
    integer(kind_int),intent(in)    :: ni,nst,ned,ngn
    real(kind_real)  ,intent(in)    :: vn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)    :: nge
    real(kind_real)  ,intent(in)    :: ve(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)    :: nfs,nfe
    real(kind_real)  ,intent(out)   :: dn5(1-ngn:ni+ngn,nst:ned)
    integer(kind_int) :: i,m,st,ed
    
    !>
    !!������׾��Ƚڵ����׿ռ䵼��
    if (nfs == nbc_inter_scheme) then
        
        st=1      
                                                     
    else
        
        st=4      

    end if
    !>

    !>
    if (nfe == nbc_inter_scheme) then

        ed=ni
        
    else

        ed=ni-3

    end if
    !>

    !>
    !!������׾��Ƚڵ����׿ռ䵼���ڵ�
    do m=nst,ned
        do i=st,ed
            dn5(i,m) =     5.0/2.0*( vn(i+1,m) - vn(i-1,m) )          -  &
                               2.0*( vn(i+2,m) - vn(i-2,m) )          +  &
                           1.0/2.0*( vn(i+3,m) - vn(i-3,m) )
        end do
    end do    

endsubroutine dn5_via_E2

subroutine dn5_via_E2e(ni,nst,ned,ngn,vn,nge,ve,nfs,nfe,dn5)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_inter_scheme
    use mod_variables, only : nghnode
    implicit none
    integer(kind_int),intent(in)    :: ni,nst,ned,ngn
    real(kind_real)  ,intent(in)    :: vn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)    :: nge
    real(kind_real)  ,intent(in)    :: ve(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)    :: nfs,nfe
    real(kind_real)  ,intent(out)   :: dn5(1-ngn:ni+ngn,nst:ned)
    integer(kind_int) :: i,m,st,ed
    
    !>
    !!������׾��Ƚڵ����׿ռ䵼��
    if (nfs == nbc_inter_scheme) then
        
        st=1-ngn+nghnode
                                                     
    else
        
        st=4      

    end if
    !>

    !>
    if (nfe == nbc_inter_scheme) then

        ed=ni+ngn-nghnode
        
    else

        ed=ni-3

    end if
    !>

    !>
    !!������׾��Ƚڵ����׿ռ䵼���ڵ�
    do m=nst,ned
        do i=st,ed
            dn5(i,m) =     5.0/2.0*( vn(i+1,m) - vn(i-1,m) )          -  &
                               2.0*( vn(i+2,m) - vn(i-2,m) )          +  &
                           1.0/2.0*( vn(i+3,m) - vn(i-3,m) )
        end do
    end do    

endsubroutine dn5_via_E2e

subroutine muscl2pv(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,one,fourth,three
    use mod_constants, only : nbc_inter_scheme,nbc_intbc_scheme
    use mod_variables, only : rlimit,plimit
    use mod_variables, only : ckmuscl,cbmuscl
    implicit none
    integer(kind_int), intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  , intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(in)  :: sn(1-ngn:ni+ngn,1:3)
    integer(kind_int), intent(in)  :: nfs,nfe
    integer(kind_int), intent(in)  :: nge
    real(kind_real)  , intent(out) :: vl(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(out) :: vr(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , external    :: sub_limit
    integer(kind_int), intent(in)  :: chkid(nst:ned)
    real(kind_real),   intent(in)  :: chklim(nst:ned)
    integer(kind_int), intent(out) :: ndec
    real(kind_real)   :: ck1,ck2
    real(kind_real)   :: dv(1-ngn:ni+ngn)
    real(kind_real)   :: dvl(1-ngn:ni+ngn),dvr(1-ngn:ni+ngn)
    integer(kind_int) :: i,i1,m,st,ed,stc,edc

    ndec = 0

    ck1 = fourth*(one - ckmuscl)
    ck2 = fourth*(one + ckmuscl)

    if (nfs == nbc_inter_scheme) then
        st = -nge
        stc = st
    elseif (nfs == nbc_intbc_scheme) then
        st = -nge
        stc = st
    else
        st = 2
        stc = 0

        do m=nst,ned
            vl(0,m) = ( 3.0*vn(1,m) - vn(2,m) ) / 2.0    ! second-order
            !!vl(0,m) = vn(1,m)                            ! first-order

            vr(0,m) = vl(0,m)

            vl(1,m) = (     vn(1,m) + vn(2,m) ) / 2.0    ! second-order
            !!vl(1,m) = vn(1,m)                            ! first-order
        end do

        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vl(0,m) < chklim(m) ) then
                    vl(0,m) = vn(1,m)
                    ndec = ndec + 1
                end if
            end if
        end do
    end if

    if (nfe == nbc_inter_scheme) then
        ed = ni + nge
        edc = ed
    elseif (nfe == nbc_intbc_scheme) then
        ed = ni + nge
        edc = ed
    else
        ed = ni - 2
        edc = ni

        do m=nst,ned
            vr(ni,m) = ( 3.0*vn(ni,m) - vn(ni-1,m) ) / 2.0   ! second-order
            !!vr(ni,m) = vn(ni,m)                              ! first-order

            vl(ni,m) = vr(ni,m)

            vr(ni-1,m) = (   vn(ni,m) + vn(ni-1,m) ) / 2.0   ! second-order
            !!vr(ni-1,m) = vn(ni,m)                            ! first-order
        end do

        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vr(ni,m) < chklim(m) ) then
                    vr(ni,m) = vn(ni,m)
                    ndec = ndec + 1
                end if
            end if
        end do
    end if

    do m=nst,ned
       do i=st-1,ed+1
          dv(i) = vn(i+1,m) - vn(i,m)
       end do

       do i=st,ed+1
          i1 = i - 1
          dvl(i) = sub_limit(dv(i ), cbmuscl*dv(i1))
          dvr(i) = sub_limit(dv(i1), cbmuscl*dv(i ))
       end do

       do i=st,ed+1
          vl(i,m) = vn(i ,m) + ck1*dvr(i ) + ck2*dvl(i )
       end do

       do i=st-1,ed
          i1 = i + 1
          vr(i,m) = vn(i1,m) - ck1*dvl(i1) - ck2*dvr(i1)
       end do
    end do

    do i=stc,edc
        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vl(i,m) < chklim(m) ) then
                    vl(i,m) = vn(i,m)
                    ndec = ndec + 1
                end if

                if (vr(i,m) < chklim(m) ) then
                    vr(i,m) = vn(i+1,m)
                    ndec = ndec + 1
                end if
           end if
        end do
    end do

end subroutine muscl2pv

subroutine muscl2pvcc(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,one,fourth,three
    use mod_constants, only : nbc_inter_scheme,nbc_intbc_scheme
    use mod_variables, only : rlimit,plimit
    use mod_variables, only : ckmuscl,cbmuscl
    implicit none
    integer(kind_int), intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  , intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(in)  :: sn(1-ngn:ni+ngn,1:3)
    integer(kind_int), intent(in)  :: nfs,nfe
    integer(kind_int), intent(in)  :: nge
    real(kind_real)  , intent(out) :: vl(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(out) :: vr(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , external    :: sub_limit
    integer(kind_int), intent(in)  :: chkid(nst:ned)
    real(kind_real),   intent(in)  :: chklim(nst:ned)
    integer(kind_int), intent(out) :: ndec
    real(kind_real)   :: ck1,ck2
    real(kind_real)   :: dv(1-ngn:ni+ngn)
    real(kind_real)   :: dvl(1-ngn:ni+ngn),dvr(1-ngn:ni+ngn)
    integer(kind_int) :: i,i1,m,st,ed,stc,edc

    ndec = 0

    ck1 = fourth*(one - ckmuscl)
    ck2 = fourth*(one + ckmuscl)

    if (nfs == nbc_inter_scheme) then
        st = 1 - nge
        stc = st
    elseif (nfs == nbc_intbc_scheme) then
        st = 1 - nge
        stc = st
    else
        st = 3
        stc = 1

        do m=nst,ned
            vl(1,m) = ( 3.0*vn(1,m) - vn(2,m) ) / 2.0    ! second-order
            !!vl(1,m) = vn(1,m)                            ! first-order

            vr(1,m) = vl(1,m)

            vl(2,m) = (     vn(1,m) + vn(2,m) ) / 2.0    ! second-order
            !!vl(2,m) = vn(1,m)                            ! first-order
        end do

        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vl(1,m) < chklim(m) ) then
                    vl(1,m) = vn(1,m)
                    ndec = ndec + 1
                end if
            end if
        end do
    end if

    if (nfe == nbc_inter_scheme) then
        ed = ni + 1 + nge
        edc = ed
    elseif (nfe == nbc_intbc_scheme) then
        ed = ni + 1 + nge
        edc = ed
    else
        ed = ni - 1
        edc = ni + 1

        do m=nst,ned
            vr(ni+1,m) = ( 3.0*vn(ni,m) - vn(ni-1,m) ) / 2.0   ! second-order
            !!vr(ni+1,m) = vn(ni,m)                              ! first-order

            vl(ni+1,m) = vr(ni+1,m)

            vr(ni,m) = (   vn(ni,m) + vn(ni-1,m) ) / 2.0   ! second-order
            !!vr(ni,m) = vn(ni,m)                            ! first-order
        end do

        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vr(ni+1,m) < chklim(m) ) then
                    vr(ni+1,m) = vn(ni,m)
                    ndec = ndec + 1
                end if
            end if
        end do
    end if

    do m=nst,ned
       do i=st-2,ed
          dv(i+1) = vn(i+1,m) - vn(i,m)
       end do

       do i=st-1,ed
          i1 = i - 1
          dvl(i+1) = sub_limit(dv(i ), cbmuscl*dv(i1))
          dvr(i+1) = sub_limit(dv(i1), cbmuscl*dv(i ))
       end do

       do i=st-1,ed
          vl(i+1,m) = vn(i ,m) + ck1*dvr(i ) + ck2*dvl(i )
       end do

       do i=st-2,ed-1
          i1 = i + 1
          vr(i+1,m) = vn(i1,m) - ck1*dvl(i1) - ck2*dvr(i1)
       end do
    end do

    do i=stc,edc
        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vl(i,m) < chklim(m) ) then
                    vl(i,m) = vn(i,m)
                    ndec = ndec + 1
                end if

                if (vr(i,m) < chklim(m) ) then
                    vr(i,m) = vn(i+1,m)
                    ndec = ndec + 1
                end if
           end if
        end do
    end do

end subroutine muscl2pvcc

function nolimiter(x,y) result(lim)
    use mod_kndconsts, only : kind_real
    use mod_constants, only : zero
    implicit none
    real(kind_real) :: x,y
    real(kind_real) :: lim

    lim = zero !!x

end function nolimiter

function minmod(x,y) result(lim)
    use mod_kndconsts, only : kind_real
    use mod_constants, only : half,one
    implicit none
    real(kind_real) :: x,y
    real(kind_real) :: lim

    lim = half*(sign(one,x) + sign(one,y))* &
                min(abs(x),abs(y))

end function minmod

function vanleer(x,y) result(lim)
    use mod_kndconsts, only : kind_real
    use mod_constants, only : one
    implicit none
    real(kind_real) :: x,y
    real(kind_real) :: lim
    real(kind_real), parameter :: eps=1.0e-15

    lim = (sign(one,x) + sign(one,y))*x*y/ &
          (abs(x) + abs(y) + eps)

end function vanleer

function vanalbada(x,y) result(limit)
    use mod_kndconsts, only : kind_real
    use mod_constants, only : one
    implicit none
    real(kind_real) :: x,y
    real(kind_real) :: limit
    real(kind_real) :: x2,y2,oxy2
    real(kind_real), parameter :: eps=1.0e-6

    x2 = x*x + eps
    y2 = y*y + eps
    oxy2 = one/(x2 + y2)
    x2 = x2 * oxy2
    y2 = y2 * oxy2
    limit = x2 * y + y2 * x

end function vanalbada

subroutine wcns5pv(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,one
    use mod_constants, only : nbc_inter_scheme
    use mod_variables, only : rlimit,plimit
    implicit none
    integer(kind_int), intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  , intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(in)  :: sn(1-ngn:ni+ngn,1:3)
    integer(kind_int), intent(in)  :: nfs,nfe
    integer(kind_int), intent(in)  :: nge
    real(kind_real)  , intent(out) :: vl(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(out) :: vr(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , external    :: sub_limit
    integer(kind_int), intent(in)  :: chkid(nst:ned)
    real(kind_real),   intent(in)  :: chklim(nst:ned)
    integer(kind_int), intent(out) :: ndec
    real(kind_real), parameter :: cl(3) = (/1.D0/16.D0,10.D0/16.D0,5.D0/16.D0/)
    real(kind_real), parameter :: cr(3) = (/5.D0/16.D0,10.D0/16.D0,1.D0/16.D0/)
    real(kind_real), parameter :: eps   = 1.0e-6
    real(kind_real)            :: is,eis2,bs,bl(3),br(3)
    real(kind_real)            :: s(1:3,1-ngn:ni+ngn,nst:ned)
    real(kind_real)            :: g(1:3,1-ngn:ni+ngn,nst:ned)
    real(kind_real)            :: wl(1:3,1-ngn:ni+ngn,nst:ned)
    real(kind_real)            :: wr(1:3,1-ngn:ni+ngn,nst:ned)
    integer(kind_int)          :: i,i1,m,n,st,ed,stc,edc

    ndec = 0

    if (nfs == nbc_inter_scheme) then
        st = -nge
        stc = st
    else
        st = 3
        stc = 0

        do m=nst,ned
            vl(0,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                  ! third-order
            vl(1,m) = (  3.0*vn(1,m) +  6.0*vn(2,m) -      vn(3,m) ) / 8.0                  ! third-order

!!            vl(0,m) = (35.0*vn(1,m) - 35.0*vn(2,m) + 21.0*vn(3,m) - 5.0*vn(4,m))/16.0
!!            vl(1,m) = ( 5.0*vn(1,m) + 15.0*vn(2,m) -  5.0*vn(3,m) +     vn(4,m))/16.0
            vl(2,m) = (    -vn(1,m) +  9.0*vn(2,m) +  9.0*vn(3,m) -     vn(4,m))/16.0
            vr(0,m) = vl(0,m)
            vr(1,m) = vl(1,m)
        end do

        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vl(0,m) < chklim(m) ) then
                    vl(0,m) = vn(1,m)
                    ndec = ndec + 1
                end if
            end if
        end do
    end if

    if (nfe == nbc_inter_scheme) then
        ed = ni + nge
        edc = ed
    else
        ed = ni - 3
        edc = ni

        do m=nst,ned
            vr(ni  ,m) = ( 15.0*vn(ni,m) - 10.0*vn(ni-1,m) +  3.0*vn(ni-2,m) ) / 8.0                ! third-order
            vr(ni-1,m) = (  3.0*vn(ni,m) +  6.0*vn(ni-1,m) -      vn(ni-2,m) ) / 8.0                ! third-order

!!            vr(ni  ,m) = (35.0*vn(ni,m) - 35.0*vn(ni-1,m) + 21.0*vn(ni-2,m) - 5.0*vn(ni-3,m))/16.0
!!            vr(ni-1,m) = ( 5.0*vn(ni,m) + 15.0*vn(ni-1,m) -  5.0*vn(ni-2,m) +     vn(ni-3,m))/16.0
            vr(ni-2,m) = (    -vn(ni,m) +  9.0*vn(ni-1,m) +  9.0*vn(ni-2,m) -     vn(ni-3,m))/16.0
            vl(ni  ,m) = vr(ni  ,m)
            vl(ni-1,m) = vr(ni-1,m)
        end do

        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vr(ni,m) < chklim(m) ) then
                    vr(ni,m) = vn(ni,m)
                    ndec = ndec + 1
                end if
            end if
        end do
    end if

    do m=nst,ned
        i = st-1
        s(2,i,m) = vn(i-1,m) - 2.0*vn(i  ,m) + vn(i+1,m)
        s(3,i,m) = vn(i  ,m) - 2.0*vn(i+1,m) + vn(i+2,m)

        do i=st,ed+1
            g(1,i,m) = 0.5*(     vn(i-2,m) - 4.0*vn(i-1,m) + 3.0*vn(i,  m))
            g(2,i,m) = 0.5*(     vn(i+1,m) -     vn(i-1,m)               )
            g(3,i,m) = 0.5*(-3.0*vn(i,  m) + 4.0*vn(i+1,m) -     vn(i+2,m))

            s(1,i,m) = s(2,i-1,m) !vn(i-2,m) - 2.0*vn(i-1,m) + vn(i,  m)
            s(2,i,m) = s(3,i-1,m) !vn(i-1,m) - 2.0*vn(i  ,m) + vn(i+1,m) 
            s(3,i,m) = vn(i,m) - 2.0*vn(i+1,m) + vn(i+2,m)

            do n=1,3
                is = g(n,i,m)*g(n,i,m) + s(n,i,m)*s(n,i,m)
                eis2 = (eps+is)**2
                bl(n) = cl(n)/eis2
                br(n) = cr(n)/eis2
            end do

            bs = bl(1) + bl(2) + bl(3)
            do n=1,3
                wl(n,i,m) = bl(n)/bs
            end do

            bs = br(1) + br(2) + br(3)
            do n=1,3
                wr(n,i,m) = br(n)/bs
            end do
        end do
    end do

    do m=nst,ned
        do i=st,ed+1
            vl(i,m) = vn(i ,m) + 0.125*(wl(1,i ,m)*(s(1,i ,m)+4.0*g(1,i ,m)) + &
                                        wl(2,i ,m)*(s(2,i ,m)+4.0*g(2,i ,m)) + &
                                        wl(3,i ,m)*(s(3,i ,m)+4.0*g(3,i ,m)) )
        end do

        do i=st-1,ed
            i1 = i+1
            vr(i,m) = vn(i1,m) + 0.125*(wr(1,i1,m)*(s(1,i1,m)-4.0*g(1,i1,m)) + &
                                        wr(2,i1,m)*(s(2,i1,m)-4.0*g(2,i1,m)) + &
                                        wr(3,i1,m)*(s(3,i1,m)-4.0*g(3,i1,m)) )  
        end do
    end do

    do i=stc,edc
        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vl(i,m) < chklim(m) ) then
                    vl(i,m) = vn(i,m)
                    ndec = ndec + 1
                end if

                if (vr(i,m) < chklim(m) ) then
                    vr(i,m) = vn(i+1,m)
                    ndec = ndec + 1
                end if
           end if
        end do
    end do

end subroutine wcns5pv

subroutine wcns5pvcc(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,one
    use mod_constants, only : nbc_inter_scheme
    use mod_variables, only : rlimit,plimit
    implicit none
    integer(kind_int), intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  , intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(in)  :: sn(1-ngn:ni+ngn,1:3)
    integer(kind_int), intent(in)  :: nfs,nfe
    integer(kind_int), intent(in)  :: nge
    real(kind_real)  , intent(out) :: vl(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(out) :: vr(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , external    :: sub_limit
    integer(kind_int), intent(in)  :: chkid(nst:ned)
    real(kind_real),   intent(in)  :: chklim(nst:ned)
    integer(kind_int), intent(out) :: ndec
    real(kind_real), parameter :: cl(3) = (/1.D0/16.D0,10.D0/16.D0,5.D0/16.D0/)
    real(kind_real), parameter :: cr(3) = (/5.D0/16.D0,10.D0/16.D0,1.D0/16.D0/)
    real(kind_real), parameter :: eps   = 1.0e-6
    real(kind_real)            :: is,eis2,bs,bl(3),br(3)
    real(kind_real)            :: s(1:3,1-ngn:ni+ngn,nst:ned)
    real(kind_real)            :: g(1:3,1-ngn:ni+ngn,nst:ned)
    real(kind_real)            :: wl(1:3,1-ngn:ni+ngn,nst:ned)
    real(kind_real)            :: wr(1:3,1-ngn:ni+ngn,nst:ned)
    integer(kind_int)          :: i,i1,m,n,st,ed,stc,edc

    ndec = 0

    if (nfs == nbc_inter_scheme) then
        st = 1 - nge
        stc = st
    else
        st = 4
        stc = 1

        do m=nst,ned
            vl(1,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                  ! third-order
            vl(2,m) = (  3.0*vn(1,m) +  6.0*vn(2,m) -      vn(3,m) ) / 8.0                  ! third-order

!!            vl(1,m) = (35.0*vn(1,m) - 35.0*vn(2,m) + 21.0*vn(3,m) - 5.0*vn(4,m))/16.0
!!            vl(2,m) = ( 5.0*vn(1,m) + 15.0*vn(2,m) -  5.0*vn(3,m) +     vn(4,m))/16.0
            vl(3,m) = (    -vn(1,m) +  9.0*vn(2,m) +  9.0*vn(3,m) -     vn(4,m))/16.0
            vr(1,m) = vl(1,m)
            vr(2,m) = vl(2,m)
        end do

        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vl(1,m) < chklim(m) ) then
                    vl(1,m) = vn(1,m)
                    ndec = ndec + 1
                end if
            end if
        end do
    end if

    if (nfe == nbc_inter_scheme) then
        ed = ni + 1 + nge
        edc = ed
    else
        ed = ni - 2
        edc = ni + 1

        do m=nst,ned
            vr(ni+1,m) = ( 15.0*vn(ni,m) - 10.0*vn(ni-1,m) +  3.0*vn(ni-2,m) ) / 8.0                ! third-order
            vr(ni  ,m) = (  3.0*vn(ni,m) +  6.0*vn(ni-1,m) -      vn(ni-2,m) ) / 8.0                ! third-order

!!            vr(ni+1,m) = (35.0*vn(ni,m) - 35.0*vn(ni-1,m) + 21.0*vn(ni-2,m) - 5.0*vn(ni-3,m))/16.0
!!            vr(ni  ,m) = ( 5.0*vn(ni,m) + 15.0*vn(ni-1,m) -  5.0*vn(ni-2,m) +     vn(ni-3,m))/16.0
            vr(ni-1,m) = (    -vn(ni,m) +  9.0*vn(ni-1,m) +  9.0*vn(ni-2,m) -     vn(ni-3,m))/16.0
            vl(ni+1,m) = vr(ni+1,m)
            vl(ni  ,m) = vr(ni  ,m)
        end do

        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vr(ni+1,m) < chklim(m) ) then
                    vr(ni+1,m) = vn(ni,m)
                    ndec = ndec + 1
                end if
            end if
        end do
    end if

    do m=nst,ned
        i = st-2
        s(2,i,m) = vn(i-1,m) - 2.0*vn(i  ,m) + vn(i+1,m)
        s(3,i,m) = vn(i  ,m) - 2.0*vn(i+1,m) + vn(i+2,m)

        do i=st-1,ed
            g(1,i,m) = 0.5*(     vn(i-2,m) - 4.0*vn(i-1,m) + 3.0*vn(i,  m))
            g(2,i,m) = 0.5*(     vn(i+1,m) -     vn(i-1,m)               )
            g(3,i,m) = 0.5*(-3.0*vn(i,  m) + 4.0*vn(i+1,m) -     vn(i+2,m))

            s(1,i,m) = s(2,i-1,m) !vn(i-2,m) - 2.0*vn(i-1,m) + vn(i,  m)
            s(2,i,m) = s(3,i-1,m) !vn(i-1,m) - 2.0*vn(i  ,m) + vn(i+1,m) 
            s(3,i,m) = vn(i,m) - 2.0*vn(i+1,m) + vn(i+2,m)

            do n=1,3
                is = g(n,i,m)*g(n,i,m) + s(n,i,m)*s(n,i,m)
                eis2 = (eps+is)**2
                bl(n) = cl(n)/eis2
                br(n) = cr(n)/eis2
            end do

            bs = bl(1) + bl(2) + bl(3)
            do n=1,3
                wl(n,i,m) = bl(n)/bs
            end do

            bs = br(1) + br(2) + br(3)
            do n=1,3
                wr(n,i,m) = br(n)/bs
            end do
        end do
    end do

    do m=nst,ned
        do i=st-1,ed
            vl(i+1,m) = vn(i ,m) + 0.125*(wl(1,i ,m)*(s(1,i ,m)+4.0*g(1,i ,m)) + &
                                          wl(2,i ,m)*(s(2,i ,m)+4.0*g(2,i ,m)) + &
                                          wl(3,i ,m)*(s(3,i ,m)+4.0*g(3,i ,m)) )
        end do

        do i=st-2,ed-1
            i1 = i+1
            vr(i+1,m) = vn(i1,m) + 0.125*(wr(1,i1,m)*(s(1,i1,m)-4.0*g(1,i1,m)) + &
                                          wr(2,i1,m)*(s(2,i1,m)-4.0*g(2,i1,m)) + &
                                          wr(3,i1,m)*(s(3,i1,m)-4.0*g(3,i1,m)) )  
        end do
    end do

    do i=stc,edc
        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vl(i,m) < chklim(m) ) then
                    vl(i,m) = vn(i,m)
                    ndec = ndec + 1
                end if

                if (vr(i,m) < chklim(m) ) then
                    vr(i,m) = vn(i+1,m)
                    ndec = ndec + 1
                end if
           end if
        end do
    end do

end subroutine wcns5pvcc

subroutine wcns5cv(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,one
    use mod_constants, only : nbc_inter_scheme
    implicit none
    integer(kind_int), intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  , intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(in)  :: sn(1-ngn:ni+ngn,1:3)
    integer(kind_int), intent(in)  :: nfs,nfe
    integer(kind_int), intent(in)  :: nge
    real(kind_real)  , intent(out) :: vl(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(out) :: vr(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , external    :: sub_limit
    integer(kind_int), intent(in)  :: chkid(nst:ned)
    real(kind_real),   intent(in)  :: chklim(nst:ned)
    integer(kind_int), intent(out) :: ndec
    real(kind_real), parameter :: cl(3) = (/1.D0/16.D0,10.D0/16.D0,5.D0/16.D0/)
    real(kind_real), parameter :: cr(3) = (/5.D0/16.D0,10.D0/16.D0,1.D0/16.D0/)
    real(kind_real), parameter :: eps   = 1.0e-6
    real(kind_real)            :: is,eis2,bs,bl(3),br(3)
    real(kind_real)            :: nt,nx,ny,nz,fp(nst:ned)
    real(kind_real)            :: vp(nst:ned),vc(nst:ned)
    real(kind_real)            :: q(1-ngn:ni+ngn,nst:ned)
    real(kind_real)            :: s(1:3,1-ngn:ni+ngn,nst:ned)
    real(kind_real)            :: g(1:3,1-ngn:ni+ngn,nst:ned)
    real(kind_real)            :: wl(1:3,1-ngn:ni+ngn,nst:ned)
    real(kind_real)            :: wr(1:3,1-ngn:ni+ngn,nst:ned)
    integer(kind_int)          :: i,i1,m,n,st,ed,stc,edc

    ndec = 0

    if (nfs == nbc_inter_scheme) then
        st = -nge
        stc = st
    else
        st = 3
        stc = 0

        do m=nst,ned
            !!vl(0,m) = (35.0*vn(1,m) - 35.0*vn(2,m) + 21.0*vn(3,m) - 5.0*vn(4,m))/16.0  !cic
            !!vl(1,m) = ( 5.0*vn(1,m) + 15.0*vn(2,m) -  5.0*vn(3,m) +     vn(4,m))/16.0
            vl(2,m) = (    -vn(1,m) +  9.0*vn(2,m) +  9.0*vn(3,m) -     vn(4,m))/16.0
            vl(0,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                  ! third-order
            vl(1,m) = (  3.0*vn(1,m) +  6.0*vn(2,m) -      vn(3,m) ) / 8.0                  ! third-order
            !!vl(2,m) = ( 35.0*vn(1,m) - 42.0*vn(2,m) + 15.0*vn(3,m) ) / 8.0 

            vr(0,m) = vl(0,m)
            vr(1,m) = vl(1,m)
        end do

        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vl(0,m) < chklim(m) ) then
                    vl(0,m) = vn(1,m)
                    ndec = ndec + 1
                end if
            end if
        end do
    end if

    if (nfe == nbc_inter_scheme) then
        ed = ni + nge
        edc = ed
    else
        ed = ni - 3
        edc = ni

        do m=nst,ned
            !!vr(ni  ,m) = (35.0*vn(ni,m) - 35.0*vn(ni-1,m) + 21.0*vn(ni-2,m) - 5.0*vn(ni-3,m))/16.0 !cic
            !!vr(ni-1,m) = ( 5.0*vn(ni,m) + 15.0*vn(ni-1,m) -  5.0*vn(ni-2,m) +     vn(ni-3,m))/16.0
            vr(ni-2,m) = (    -vn(ni,m) +  9.0*vn(ni-1,m) +  9.0*vn(ni-2,m) -     vn(ni-3,m))/16.0
            vr(ni  ,m) = ( 15.0*vn(ni,m) - 10.0*vn(ni-1,m) +  3.0*vn(ni-2,m) ) / 8.0                ! third-order
            vr(ni-1,m) = (  3.0*vn(ni,m) +  6.0*vn(ni-1,m) -      vn(ni-2,m) ) / 8.0                ! third-order
            !!vr(ni-2,m) = ( 35.0*vn(ni,m) - 42.0*vn(ni-1,m) + 15.0*vn(ni-2,m) ) / 8.0 

            vl(ni  ,m) = vr(ni  ,m)
            vl(ni-1,m) = vr(ni-1,m)
        end do

        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vr(ni,m) < chklim(m) ) then
                    vr(ni,m) = vn(ni,m)
                    ndec = ndec + 1
                end if
            end if
        end do
    end if


    do m=nst,ned
        i = st-1
        s(2,i,m) = vn(i-1,m) - 2.0*vn(i  ,m) + vn(i+1,m)
        s(3,i,m) = vn(i  ,m) - 2.0*vn(i+1,m) + vn(i+2,m)

        do i=st,ed+1
            g(1,i,m) = 0.5*(     vn(i-2,m) - 4.0*vn(i-1,m) + 3.0*vn(i,  m))
            g(2,i,m) = 0.5*(     vn(i+1,m) -     vn(i-1,m)               )
            g(3,i,m) = 0.5*(-3.0*vn(i,  m) + 4.0*vn(i+1,m) -     vn(i+2,m))

            s(1,i,m) = s(2,i-1,m) !vn(i-2,m) - 2.0*vn(i-1,m) + vn(i,  m)
            s(2,i,m) = s(3,i-1,m) !vn(i-1,m) - 2.0*vn(i  ,m) + vn(i+1,m) 
            s(3,i,m) = vn(i,m) - 2.0*vn(i+1,m) + vn(i+2,m)
        end do
    end do

    do i=st,ed+1
        nx = sn(i,1)
        ny = sn(i,2)
        nz = sn(i,3)

        do m=nst,ned
            vp(m) = vn(i,m)
        end do
        call LxQ_p(nst,ned,vp,nt,nx,ny,nz,vp,vc) 
        do m=nst,ned
            q(i,m) = vc(m)
        end do
         
        do n=1,3
            do m=nst,ned
                fp(m) = g(n,i,m)
            end do
            call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
            do m=nst,ned
                g(n,i,m) = vc(m)
            end do

            do m=nst,ned
                fp(m) = s(n,i,m)
            end do
            call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
            do m=nst,ned
                s(n,i,m) = vc(m)
            end do
        end do
    end do

    do m=nst,ned
        do i=st,ed+1
            do n=1,3
                is = g(n,i,m)*g(n,i,m) + s(n,i,m)*s(n,i,m)
                eis2 = (eps+is)**2
                bl(n) = cl(n)/eis2
                br(n) = cr(n)/eis2
            end do


            bs = bl(1) + bl(2) + bl(3)
            do n=1,3
                wl(n,i,m) = bl(n)/bs
            end do

            bs = br(1) + br(2) + br(3)
            do n=1,3
                wr(n,i,m) = br(n)/bs
            end do
        end do
    end do

    do m=nst,ned
        do i=st,ed+1
            vl(i,m) = q(i ,m) + 0.125*(wl(1,i ,m)*(s(1,i ,m)+4.0*g(1,i ,m)) + &
                                       wl(2,i ,m)*(s(2,i ,m)+4.0*g(2,i ,m)) + &
                                       wl(3,i ,m)*(s(3,i ,m)+4.0*g(3,i ,m)) )
        end do

        do i=st-1,ed
            i1 = i+1
            vr(i,m) = q(i1,m) + 0.125*(wr(1,i1,m)*(s(1,i1,m)-4.0*g(1,i1,m)) + &
                                       wr(2,i1,m)*(s(2,i1,m)-4.0*g(2,i1,m)) + &
                                       wr(3,i1,m)*(s(3,i1,m)-4.0*g(3,i1,m)) )  
        end do
    end do

    do i=st,ed+1
        nx = sn(i,1)
        ny = sn(i,2)
        nz = sn(i,3)

        do m=nst,ned
            vp(m) = vn(i,m)
            vc(m) = vl(i,m)
        end do

        call RxQ_p(nst,ned,vp,nt,nx,ny,nz,vc,fp) 
        do m=nst,ned
            vl(i,m) = fp(m)
        end do
    end do

    do i=st-1,ed
        i1 = i+1
        nx = sn(i1,1)
        ny = sn(i1,2)
        nz = sn(i1,3)

        do m=nst,ned
            vp(m) = vn(i1,m)
            vc(m) = vr(i ,m)
        end do

        call RxQ_p(nst,ned,vp,nt,nx,ny,nz,vc,fp) 
        do m=nst,ned
            vr(i,m) = fp(m)
        end do
    end do

    do i=stc,edc
        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vl(i,m) < chklim(m) ) then
                    vl(i,m) = vn(i,m)
                    ndec = ndec + 1
                end if

                if (vr(i,m) < chklim(m) ) then
                    vr(i,m) = vn(i+1,m)
                    ndec = ndec + 1
                end if
           end if
        end do
    end do

end subroutine wcns5cv

subroutine wcns5cvcc(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,one
    use mod_constants, only : nbc_inter_scheme
    implicit none
    integer(kind_int), intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  , intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(in)  :: sn(1-ngn:ni+ngn,1:3)
    integer(kind_int), intent(in)  :: nfs,nfe
    integer(kind_int), intent(in)  :: nge
    real(kind_real)  , intent(out) :: vl(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(out) :: vr(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , external    :: sub_limit
    integer(kind_int), intent(in)  :: chkid(nst:ned)
    real(kind_real),   intent(in)  :: chklim(nst:ned)
    integer(kind_int), intent(out) :: ndec
    real(kind_real), parameter :: cl(3) = (/1.D0/16.D0,10.D0/16.D0,5.D0/16.D0/)
    real(kind_real), parameter :: cr(3) = (/5.D0/16.D0,10.D0/16.D0,1.D0/16.D0/)
    real(kind_real), parameter :: eps   = 1.0e-6
    real(kind_real)            :: is,eis2,bs,bl(3),br(3)
    real(kind_real)            :: nt,nx,ny,nz,fp(nst:ned)
    real(kind_real)            :: vp(nst:ned),vc(nst:ned)
    real(kind_real)            :: q(1-ngn:ni+ngn,nst:ned)
    real(kind_real)            :: s(1:3,1-ngn:ni+ngn,nst:ned)
    real(kind_real)            :: g(1:3,1-ngn:ni+ngn,nst:ned)
    real(kind_real)            :: wl(1:3,1-ngn:ni+ngn,nst:ned)
    real(kind_real)            :: wr(1:3,1-ngn:ni+ngn,nst:ned)
    integer(kind_int)          :: i,i1,m,n,st,ed,stc,edc

    ndec = 0

    if (nfs == nbc_inter_scheme) then
        st = 1 - nge
        stc = st
    else
        st = 4
        stc = 1

        do m=nst,ned
            !!vl(1,m) = (35.0*vn(1,m) - 35.0*vn(2,m) + 21.0*vn(3,m) - 5.0*vn(4,m))/16.0  !cic
            !!vl(2,m) = ( 5.0*vn(1,m) + 15.0*vn(2,m) -  5.0*vn(3,m) +     vn(4,m))/16.0
            vl(3,m) = (    -vn(1,m) +  9.0*vn(2,m) +  9.0*vn(3,m) -     vn(4,m))/16.0
            vl(1,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                  ! third-order
            vl(2,m) = (  3.0*vn(1,m) +  6.0*vn(2,m) -      vn(3,m) ) / 8.0                  ! third-order
            !!vl(3,m) = ( 35.0*vn(1,m) - 42.0*vn(2,m) + 15.0*vn(3,m) ) / 8.0 

            vr(1,m) = vl(1,m)
            vr(2,m) = vl(2,m)
        end do

        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vl(1,m) < chklim(m) ) then
                    vl(1,m) = vn(1,m)
                    ndec = ndec + 1
                end if
            end if
        end do
    end if

    if (nfe == nbc_inter_scheme) then
        ed = ni + 1 + nge
        edc = ed
    else
        ed = ni - 2
        edc = ni + 1

        do m=nst,ned
            !!vr(ni+1,m) = (35.0*vn(ni,m) - 35.0*vn(ni-1,m) + 21.0*vn(ni-2,m) - 5.0*vn(ni-3,m))/16.0 !cic
            !!vr(ni  ,m) = ( 5.0*vn(ni,m) + 15.0*vn(ni-1,m) -  5.0*vn(ni-2,m) +     vn(ni-3,m))/16.0
            vr(ni-1,m) = (    -vn(ni,m) +  9.0*vn(ni-1,m) +  9.0*vn(ni-2,m) -     vn(ni-3,m))/16.0
            vr(ni+1,m) = ( 15.0*vn(ni,m) - 10.0*vn(ni-1,m) +  3.0*vn(ni-2,m) ) / 8.0                ! third-order
            vr(ni  ,m) = (  3.0*vn(ni,m) +  6.0*vn(ni-1,m) -      vn(ni-2,m) ) / 8.0                ! third-order
            !!vr(ni-1,m) = ( 35.0*vn(ni,m) - 42.0*vn(ni-1,m) + 15.0*vn(ni-2,m) ) / 8.0 

            vl(ni+1,m) = vr(ni+1,m)
            vl(ni  ,m) = vr(ni  ,m)
        end do

        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vr(ni+1,m) < chklim(m) ) then
                    vr(ni+1,m) = vn(ni,m)
                    ndec = ndec + 1
                end if
            end if
        end do
    end if


    do m=nst,ned
        i = st-2
        s(2,i,m) = vn(i-1,m) - 2.0*vn(i  ,m) + vn(i+1,m)
        s(3,i,m) = vn(i  ,m) - 2.0*vn(i+1,m) + vn(i+2,m)

        do i=st-1,ed
            g(1,i,m) = 0.5*(     vn(i-2,m) - 4.0*vn(i-1,m) + 3.0*vn(i,  m))
            g(2,i,m) = 0.5*(     vn(i+1,m) -     vn(i-1,m)               )
            g(3,i,m) = 0.5*(-3.0*vn(i,  m) + 4.0*vn(i+1,m) -     vn(i+2,m))

            s(1,i,m) = s(2,i-1,m) !vn(i-2,m) - 2.0*vn(i-1,m) + vn(i,  m)
            s(2,i,m) = s(3,i-1,m) !vn(i-1,m) - 2.0*vn(i  ,m) + vn(i+1,m) 
            s(3,i,m) = vn(i,m) - 2.0*vn(i+1,m) + vn(i+2,m)
        end do
    end do

    do i=st-1,ed
        nx = sn(i,1)
        ny = sn(i,2)
        nz = sn(i,3)

        do m=nst,ned
            vp(m) = vn(i,m)
        end do
        call LxQ_p(nst,ned,vp,nt,nx,ny,nz,vp,vc) 
        do m=nst,ned
            q(i,m) = vc(m)
        end do
         
        do n=1,3
            do m=nst,ned
                fp(m) = g(n,i,m)
            end do
            call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
            do m=nst,ned
                g(n,i,m) = vc(m)
            end do

            do m=nst,ned
                fp(m) = s(n,i,m)
            end do
            call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
            do m=nst,ned
                s(n,i,m) = vc(m)
            end do
        end do
    end do

    do m=nst,ned
        do i=st-1,ed
            do n=1,3
                is = g(n,i,m)*g(n,i,m) + s(n,i,m)*s(n,i,m)
                eis2 = (eps+is)**2
                bl(n) = cl(n)/eis2
                br(n) = cr(n)/eis2
            end do


            bs = bl(1) + bl(2) + bl(3)
            do n=1,3
                wl(n,i,m) = bl(n)/bs
            end do

            bs = br(1) + br(2) + br(3)
            do n=1,3
                wr(n,i,m) = br(n)/bs
            end do
        end do
    end do

    do m=nst,ned
        do i=st-1,ed
            vl(i+1,m) = q(i ,m) + 0.125*(wl(1,i ,m)*(s(1,i ,m)+4.0*g(1,i ,m)) + &
                                         wl(2,i ,m)*(s(2,i ,m)+4.0*g(2,i ,m)) + &
                                         wl(3,i ,m)*(s(3,i ,m)+4.0*g(3,i ,m)) )
        end do

        do i=st-2,ed-1
            i1 = i+1
            vr(i+1,m) = q(i1,m) + 0.125*(wr(1,i1,m)*(s(1,i1,m)-4.0*g(1,i1,m)) + &
                                         wr(2,i1,m)*(s(2,i1,m)-4.0*g(2,i1,m)) + &
                                         wr(3,i1,m)*(s(3,i1,m)-4.0*g(3,i1,m)) )  
        end do
    end do

    do i=st,ed+1
        nx = sn(i,1)
        ny = sn(i,2)
        nz = sn(i,3)

        do m=nst,ned
            vp(m) = vn(i,m)
            vc(m) = vl(i,m)
        end do

        call RxQ_p(nst,ned,vp,nt,nx,ny,nz,vc,fp) 
        do m=nst,ned
            vl(i,m) = fp(m)
        end do
    end do

    do i=st-1,ed
        i1 = i+1
        nx = sn(i1,1)
        ny = sn(i1,2)
        nz = sn(i1,3)

        do m=nst,ned
            vp(m) = vn(i1,m)
            vc(m) = vr(i ,m)
        end do

        call RxQ_p(nst,ned,vp,nt,nx,ny,nz,vc,fp) 
        do m=nst,ned
            vr(i,m) = fp(m)
        end do
    end do

    do i=stc,edc
        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vl(i,m) < chklim(m) ) then
                    vl(i,m) = vn(i,m)
                    ndec = ndec + 1
                end if

                if (vr(i,m) < chklim(m) ) then
                    vr(i,m) = vn(i+1,m)
                    ndec = ndec + 1
                end if
           end if
        end do
    end do

end subroutine wcns5cvcc

subroutine wcns7pv(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,one
    use mod_constants, only : nbc_inter_scheme
    use mod_variables, only : rlimit,plimit
    implicit none
    integer(kind_int), intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  , intent(inout)  :: vn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(in)  :: sn(1-ngn:ni+ngn,1:3)
    integer(kind_int), intent(in)  :: nfs,nfe
    integer(kind_int), intent(in)  :: nge
    real(kind_real)  , intent(out) :: vl(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(out) :: vr(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , external    :: sub_limit
    integer(kind_int), intent(in)  :: chkid(nst:ned)
    real(kind_real),   intent(in)  :: chklim(nst:ned)
    integer(kind_int), intent(out) :: ndec
    real(kind_real), parameter :: cl(4) = (/1.D0/64.D0,21.D0/64.D0,35.D0/64.D0,7.D0/64.D0/)
    real(kind_real), parameter :: cr(4) = (/7.D0/64.D0,35.D0/64.D0,21.D0/64.D0,1.D0/64.D0/)
    real(kind_real), parameter :: eps   = 1.0e-6
    real(kind_real)            :: is,eis2,bs,bl(4),br(4)
    real(kind_real)            :: s(1:4,1-ngn:ni+ngn,nst:ned)
    real(kind_real)            :: g(1:4,1-ngn:ni+ngn,nst:ned)
    real(kind_real)            :: t(1:4,1-ngn:ni+ngn,nst:ned)
    real(kind_real)            :: wl(1:4,1-ngn:ni+ngn,nst:ned)
    real(kind_real)            :: wr(1:4,1-ngn:ni+ngn,nst:ned)
    integer(kind_int)          :: i,i1,m,n,st,ed,stc,edc

    ndec = 0

    if (nfs == nbc_inter_scheme) then
        st = -nge
        stc = st
    else
        st = 4
        stc = 0

        do m=nst,ned
            vl(0,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                  ! third-order
            vl(1,m) = (  3.0*vn(1,m) +  6.0*vn(2,m) -      vn(3,m) ) / 8.0                  ! third-order

!!            vl(0,m) = (35.0*vn(1,m) - 35.0*vn(2,m) + 21.0*vn(3,m) - 5.0*vn(4,m))/16.0
!!            vl(1,m) = ( 5.0*vn(1,m) + 15.0*vn(2,m) -  5.0*vn(3,m) +     vn(4,m))/16.0
            vl(2,m) = (    -vn(1,m) +  9.0*vn(2,m) +  9.0*vn(3,m) -      vn(4,m))/16.0
            vr(2,m) = (-5.0*vn(1,m) + 60.0*vn(2,m) + 90.0*vn(3,m) - 20.0*vn(4,m) + 3.0*vn(5,m) ) / 128.0

            vl(3,m) = ( 3.0*vn(1,m) - 20.0*vn(2,m) + 90.0*vn(3,m) + 60.0*vn(4,m) - 5.0*vn(5,m) ) / 128.0

            vr(0,m) = vl(0,m)
            vr(1,m) = vl(1,m)
        end do

        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vl(0,m) < chklim(m) ) then
                    vl(0,m) = vn(1,m)
                    ndec = ndec + 1
                end if
            end if
        end do
    end if

    if (nfe == nbc_inter_scheme) then
        ed = ni + nge
        edc = ed
    else
        ed = ni - 4
        edc = ni

        do m=nst,ned
            vr(ni  ,m) = ( 15.0*vn(ni,m) - 10.0*vn(ni-1,m) +  3.0*vn(ni-2,m) ) / 8.0                ! third-order
            vr(ni-1,m) = (  3.0*vn(ni,m) +  6.0*vn(ni-1,m) -      vn(ni-2,m) ) / 8.0                ! third-order

!!            vr(ni  ,m) = (35.0*vn(ni,m) - 35.0*vn(ni-1,m) + 21.0*vn(ni-2,m) - 5.0*vn(ni-3,m))/16.0
!!            vr(ni-1,m) = ( 5.0*vn(ni,m) + 15.0*vn(ni-1,m) -  5.0*vn(ni-2,m) +     vn(ni-3,m))/16.0
            vr(ni-2,m) = (    -vn(ni,m) +  9.0*vn(ni-1,m) +  9.0*vn(ni-2,m) -      vn(ni-3,m))/16.0
            vl(ni-2,m) = (-5.0*vn(ni,m) + 60.0*vn(ni-1,m) + 90.0*vn(ni-2,m) - 20.0*vn(ni-3,m) + 3.0*vn(ni-4,m) ) / 128.0   ! fifth-order

            vr(ni-3,m) = ( 3.0*vn(ni,m) - 20.0*vn(ni-1,m) + 90.0*vn(ni-2,m) + 60.0*vn(ni-3,m) - 5.0*vn(ni-4,m) ) / 128.0

            vl(ni  ,m) = vr(ni  ,m)
            vl(ni-1,m) = vr(ni-1,m)
        end do

        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vr(ni,m) < chklim(m) ) then
                    vr(ni,m) = vn(ni,m)
                    ndec = ndec + 1
                end if
            end if
        end do
    end if

    do m=nst,ned
        do i=st,ed+1
            g(1,i,m) = ( -2.0*vn(i-3,m) +  9.0*vn(i-2,m) - 18.0*vn(i-1,m) + 11.0*vn(i  ,m))/6.0
            g(2,i,m) = (      vn(i-2,m) -  6.0*vn(i-1,m) +  3.0*vn(i  ,m) +  2.0*vn(i+1,m))/6.0
            g(3,i,m) = ( -2.0*vn(i-1,m) -  3.0*vn(i  ,m) +  6.0*vn(i+1,m) -      vn(i+2,m))/6.0
            g(4,i,m) = (-11.0*vn(i  ,m) + 18.0*vn(i+1,m) -  9.0*vn(i+2,m) +  2.0*vn(i+3,m))/6.0

            s(1,i,m) =     -vn(i-3,m) + 4.0*vn(i-2,m) - 5.0*vn(i-1,m) + 2.0*vn(i  ,m)
            s(2,i,m) =                      vn(i-1,m) - 2.0*vn(i  ,m) +     vn(i+1,m) 
            s(3,i,m) =      vn(i-1,m) - 2.0*vn(i  ,m) +     vn(i+1,m)
            s(4,i,m) =  2.0*vn(i  ,m) - 5.0*vn(i+1,m) + 4.0*vn(i+2,m) -     vn(i+3,m)

            t(1,i,m) = -vn(i-3,m) + 3.0*vn(i-2,m) - 3.0*vn(i-1,m) + vn(i  ,m)
            t(2,i,m) = -vn(i-2,m) + 3.0*vn(i-1,m) - 3.0*vn(i  ,m) + vn(i+1,m)
            t(3,i,m) = -vn(i-1,m) + 3.0*vn(i  ,m) - 3.0*vn(i+1,m) + vn(i+2,m)
            t(4,i,m) = -vn(i  ,m) + 3.0*vn(i+1,m) - 3.0*vn(i+2,m) + vn(i+3,m)

            do n=1,4
                is = g(n,i,m)*g(n,i,m) + s(n,i,m)*s(n,i,m) + t(n,i,m)*t(n,i,m)
                eis2 = (eps+is)**2
                bl(n) = cl(n)/eis2
                br(n) = cr(n)/eis2
            end do

            bs = bl(1) + bl(2) + bl(3) + bl(4)
            do n=1,4
                wl(n,i,m) = bl(n)/bs
            end do

            bs = br(1) + br(2) + br(3) + br(4)
            do n=1,4
                wr(n,i,m) = br(n)/bs
            end do
        end do
    end do

    do m=nst,ned
        do i=st,ed+1
            vl(i,m) = vn(i ,m) + ( wl(1,i ,m)*(24.0*g(1,i ,m)+6.0*s(1,i ,m)+t(1,i ,m)) + &
                                   wl(2,i ,m)*(24.0*g(2,i ,m)+6.0*s(2,i ,m)+t(2,i ,m)) + &
                                   wl(3,i ,m)*(24.0*g(3,i ,m)+6.0*s(3,i ,m)+t(3,i ,m)) + &
                                   wl(4,i ,m)*(24.0*g(4,i ,m)+6.0*s(4,i ,m)+t(4,i ,m)) )/48.0
        end do

        do i=st-1,ed
            i1 = i+1
            vr(i,m) = vn(i1,m) + ( wr(1,i1,m)*(-24.0*g(1,i1,m)+6.0*s(1,i1,m)-t(1,i1,m)) + &
                                   wr(2,i1,m)*(-24.0*g(2,i1,m)+6.0*s(2,i1,m)-t(2,i1,m)) + &
                                   wr(3,i1,m)*(-24.0*g(3,i1,m)+6.0*s(3,i1,m)-t(3,i1,m)) + &
                                   wr(4,i1,m)*(-24.0*g(4,i1,m)+6.0*s(4,i1,m)-t(4,i1,m)) )/48.0
        end do
    end do

    do i=stc,edc
        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vl(i,m) < chklim(m) ) then
                    vl(i,m) = vn(i,m)
                    ndec = ndec + 1
                end if

                if (vr(i,m) < chklim(m) ) then
                    vr(i,m) = vn(i+1,m)
                    ndec = ndec + 1
                end if
           end if
        end do
    end do

end subroutine wcns7pv

subroutine wcns7pvcc(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,one
    use mod_constants, only : nbc_inter_scheme
    use mod_variables, only : rlimit,plimit
    implicit none
    integer(kind_int), intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  , intent(inout)  :: vn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(in)  :: sn(1-ngn:ni+ngn,1:3)
    integer(kind_int), intent(in)  :: nfs,nfe
    integer(kind_int), intent(in)  :: nge
    real(kind_real)  , intent(out) :: vl(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(out) :: vr(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , external    :: sub_limit
    integer(kind_int), intent(in)  :: chkid(nst:ned)
    real(kind_real),   intent(in)  :: chklim(nst:ned)
    integer(kind_int), intent(out) :: ndec
    real(kind_real), parameter :: cl(4) = (/1.D0/64.D0,21.D0/64.D0,35.D0/64.D0,7.D0/64.D0/)
    real(kind_real), parameter :: cr(4) = (/7.D0/64.D0,35.D0/64.D0,21.D0/64.D0,1.D0/64.D0/)
    real(kind_real), parameter :: eps   = 1.0e-6
    real(kind_real)            :: is,eis2,bs,bl(4),br(4)
    real(kind_real)            :: s(1:4,1-ngn:ni+ngn,nst:ned)
    real(kind_real)            :: g(1:4,1-ngn:ni+ngn,nst:ned)
    real(kind_real)            :: t(1:4,1-ngn:ni+ngn,nst:ned)
    real(kind_real)            :: wl(1:4,1-ngn:ni+ngn,nst:ned)
    real(kind_real)            :: wr(1:4,1-ngn:ni+ngn,nst:ned)
    integer(kind_int)          :: i,i1,m,n,st,ed,stc,edc

    ndec = 0

    if (nfs == nbc_inter_scheme) then
        st = 1 - nge
        stc = st
    else
        st = 5
        stc = 1

        do m=nst,ned
            vl(1,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                  ! third-order
            vl(2,m) = (  3.0*vn(1,m) +  6.0*vn(2,m) -      vn(3,m) ) / 8.0                  ! third-order

!!            vl(1,m) = (35.0*vn(1,m) - 35.0*vn(2,m) + 21.0*vn(3,m) - 5.0*vn(4,m))/16.0
!!            vl(2,m) = ( 5.0*vn(1,m) + 15.0*vn(2,m) -  5.0*vn(3,m) +     vn(4,m))/16.0
            vl(3,m) = (    -vn(1,m) +  9.0*vn(2,m) +  9.0*vn(3,m) -      vn(4,m))/16.0
            vr(3,m) = (-5.0*vn(1,m) + 60.0*vn(2,m) + 90.0*vn(3,m) - 20.0*vn(4,m) + 3.0*vn(5,m) ) / 128.0

            vl(4,m) = ( 3.0*vn(1,m) - 20.0*vn(2,m) + 90.0*vn(3,m) + 60.0*vn(4,m) - 5.0*vn(5,m) ) / 128.0

            vr(1,m) = vl(1,m)
            vr(2,m) = vl(2,m)
        end do

        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vl(1,m) < chklim(m) ) then
                    vl(1,m) = vn(1,m)
                    ndec = ndec + 1
                end if
            end if
        end do
    end if

    if (nfe == nbc_inter_scheme) then
        ed = ni + 1 + nge
        edc = ed
    else
        ed = ni - 3
        edc = ni + 1

        do m=nst,ned
            vr(ni+1,m) = ( 15.0*vn(ni,m) - 10.0*vn(ni-1,m) +  3.0*vn(ni-2,m) ) / 8.0                ! third-order
            vr(ni  ,m) = (  3.0*vn(ni,m) +  6.0*vn(ni-1,m) -      vn(ni-2,m) ) / 8.0                ! third-order

!!            vr(ni+1,m) = (35.0*vn(ni,m) - 35.0*vn(ni-1,m) + 21.0*vn(ni-2,m) - 5.0*vn(ni-3,m))/16.0
!!            vr(ni  ,m) = ( 5.0*vn(ni,m) + 15.0*vn(ni-1,m) -  5.0*vn(ni-2,m) +     vn(ni-3,m))/16.0
            vr(ni-1,m) = (    -vn(ni,m) +  9.0*vn(ni-1,m) +  9.0*vn(ni-2,m) -      vn(ni-3,m))/16.0
            vl(ni-1,m) = (-5.0*vn(ni,m) + 60.0*vn(ni-1,m) + 90.0*vn(ni-2,m) - 20.0*vn(ni-3,m) + 3.0*vn(ni-4,m) ) / 128.0   ! fifth-order

            vr(ni-2,m) = ( 3.0*vn(ni,m) - 20.0*vn(ni-1,m) + 90.0*vn(ni-2,m) + 60.0*vn(ni-3,m) - 5.0*vn(ni-4,m) ) / 128.0

            vl(ni+1,m) = vr(ni+1,m)
            vl(ni  ,m) = vr(ni  ,m)
        end do

        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vr(ni+1,m) < chklim(m) ) then
                    vr(ni+1,m) = vn(ni,m)
                    ndec = ndec + 1
                end if
            end if
        end do
    end if

    do m=nst,ned
        do i=st-1,ed
            g(1,i,m) = ( -2.0*vn(i-3,m) +  9.0*vn(i-2,m) - 18.0*vn(i-1,m) + 11.0*vn(i  ,m))/6.0
            g(2,i,m) = (      vn(i-2,m) -  6.0*vn(i-1,m) +  3.0*vn(i  ,m) +  2.0*vn(i+1,m))/6.0
            g(3,i,m) = ( -2.0*vn(i-1,m) -  3.0*vn(i  ,m) +  6.0*vn(i+1,m) -      vn(i+2,m))/6.0
            g(4,i,m) = (-11.0*vn(i  ,m) + 18.0*vn(i+1,m) -  9.0*vn(i+2,m) +  2.0*vn(i+3,m))/6.0

            s(1,i,m) =     -vn(i-3,m) + 4.0*vn(i-2,m) - 5.0*vn(i-1,m) + 2.0*vn(i  ,m)
            s(2,i,m) =                      vn(i-1,m) - 2.0*vn(i  ,m) +     vn(i+1,m) 
            s(3,i,m) =      vn(i-1,m) - 2.0*vn(i  ,m) +     vn(i+1,m)
            s(4,i,m) =  2.0*vn(i  ,m) - 5.0*vn(i+1,m) + 4.0*vn(i+2,m) -     vn(i+3,m)

            t(1,i,m) = -vn(i-3,m) + 3.0*vn(i-2,m) - 3.0*vn(i-1,m) + vn(i  ,m)
            t(2,i,m) = -vn(i-2,m) + 3.0*vn(i-1,m) - 3.0*vn(i  ,m) + vn(i+1,m)
            t(3,i,m) = -vn(i-1,m) + 3.0*vn(i  ,m) - 3.0*vn(i+1,m) + vn(i+2,m)
            t(4,i,m) = -vn(i  ,m) + 3.0*vn(i+1,m) - 3.0*vn(i+2,m) + vn(i+3,m)

            do n=1,4
                is = g(n,i,m)*g(n,i,m) + s(n,i,m)*s(n,i,m) + t(n,i,m)*t(n,i,m)
                eis2 = (eps+is)**2
                bl(n) = cl(n)/eis2
                br(n) = cr(n)/eis2
            end do

            bs = bl(1) + bl(2) + bl(3) + bl(4)
            do n=1,4
                wl(n,i,m) = bl(n)/bs
            end do

            bs = br(1) + br(2) + br(3) + br(4)
            do n=1,4
                wr(n,i,m) = br(n)/bs
            end do
        end do
    end do

    do m=nst,ned
        do i=st-1,ed
            vl(i+1,m) = vn(i ,m) + ( wl(1,i ,m)*(24.0*g(1,i ,m)+6.0*s(1,i ,m)+t(1,i ,m)) + &
                                     wl(2,i ,m)*(24.0*g(2,i ,m)+6.0*s(2,i ,m)+t(2,i ,m)) + &
                                     wl(3,i ,m)*(24.0*g(3,i ,m)+6.0*s(3,i ,m)+t(3,i ,m)) + &
                                     wl(4,i ,m)*(24.0*g(4,i ,m)+6.0*s(4,i ,m)+t(4,i ,m)) )/48.0
        end do

        do i=st-2,ed-1
            i1 = i+1
            vr(i+1,m) = vn(i1,m) + ( wr(1,i1,m)*(-24.0*g(1,i1,m)+6.0*s(1,i1,m)-t(1,i1,m)) + &
                                     wr(2,i1,m)*(-24.0*g(2,i1,m)+6.0*s(2,i1,m)-t(2,i1,m)) + &
                                     wr(3,i1,m)*(-24.0*g(3,i1,m)+6.0*s(3,i1,m)-t(3,i1,m)) + &
                                     wr(4,i1,m)*(-24.0*g(4,i1,m)+6.0*s(4,i1,m)-t(4,i1,m)) )/48.0
        end do
    end do

    do i=stc,edc
        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vl(i,m) < chklim(m) ) then
                    vl(i,m) = vn(i,m)
                    ndec = ndec + 1
                end if

                if (vr(i,m) < chklim(m) ) then
                    vr(i,m) = vn(i+1,m)
                    ndec = ndec + 1
                end if
           end if
        end do
    end do

end subroutine wcns7pvcc

subroutine wcns7cv(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,one
    use mod_constants, only : nbc_inter_scheme
    use mod_variables, only : rlimit,plimit
    implicit none
    integer(kind_int), intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  , intent(inout)  :: vn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(in)  :: sn(1-ngn:ni+ngn,1:3)
    integer(kind_int), intent(in)  :: nfs,nfe
    integer(kind_int), intent(in)  :: nge
    real(kind_real)  , intent(out) :: vl(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(out) :: vr(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , external    :: sub_limit
    integer(kind_int), intent(in)  :: chkid(nst:ned)
    real(kind_real),   intent(in)  :: chklim(nst:ned)
    integer(kind_int), intent(out) :: ndec
    real(kind_real), parameter :: cl(4) = (/1.D0/64.D0,21.D0/64.D0,35.D0/64.D0,7.D0/64.D0/)
    real(kind_real), parameter :: cr(4) = (/7.D0/64.D0,35.D0/64.D0,21.D0/64.D0,1.D0/64.D0/)
    real(kind_real), parameter :: eps   = 1.0e-6
    real(kind_real)            :: is,eis2,bs,bl(4),br(4)
    real(kind_real)            :: nt,nx,ny,nz,fp(nst:ned)
    real(kind_real)            :: vp(nst:ned),vc(nst:ned)
    real(kind_real)            :: q(1-ngn:ni+ngn,nst:ned)
    real(kind_real)            :: s(1:4,1-ngn:ni+ngn,nst:ned)
    real(kind_real)            :: g(1:4,1-ngn:ni+ngn,nst:ned)
    real(kind_real)            :: t(1:4,1-ngn:ni+ngn,nst:ned)
    real(kind_real)            :: wl(1:4,1-ngn:ni+ngn,nst:ned)
    real(kind_real)            :: wr(1:4,1-ngn:ni+ngn,nst:ned)
    integer(kind_int)          :: i,i1,m,n,st,ed,stc,edc

    ndec = 0

    if (nfs == nbc_inter_scheme) then
        st = -nge
        stc = st
    else
        st = 4
        stc = 0

        do m=nst,ned
            vl(0,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                  ! third-order
            vl(1,m) = (  3.0*vn(1,m) +  6.0*vn(2,m) -      vn(3,m) ) / 8.0                  ! third-order

!!            vl(0,m) = (35.0*vn(1,m) - 35.0*vn(2,m) + 21.0*vn(3,m) - 5.0*vn(4,m))/16.0
!!            vl(1,m) = ( 5.0*vn(1,m) + 15.0*vn(2,m) -  5.0*vn(3,m) +     vn(4,m))/16.0
            vl(2,m) = (    -vn(1,m) +  9.0*vn(2,m) +  9.0*vn(3,m) -      vn(4,m))/16.0
            vr(2,m) = (-5.0*vn(1,m) + 60.0*vn(2,m) + 90.0*vn(3,m) - 20.0*vn(4,m) + 3.0*vn(5,m) ) / 128.0

            vl(3,m) = ( 3.0*vn(1,m) - 20.0*vn(2,m) + 90.0*vn(3,m) + 60.0*vn(4,m) - 5.0*vn(5,m) ) / 128.0

            vr(0,m) = vl(0,m)
            vr(1,m) = vl(1,m)
        end do

        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vl(0,m) < chklim(m) ) then
                    vl(0,m) = vn(1,m)
                    ndec = ndec + 1
                end if
            end if
        end do
    end if

    if (nfe == nbc_inter_scheme) then
        ed = ni + nge
        edc = ed
    else
        ed = ni - 4
        edc = ni

        do m=nst,ned
            vr(ni  ,m) = ( 15.0*vn(ni,m) - 10.0*vn(ni-1,m) +  3.0*vn(ni-2,m) ) / 8.0                ! third-order
            vr(ni-1,m) = (  3.0*vn(ni,m) +  6.0*vn(ni-1,m) -      vn(ni-2,m) ) / 8.0                ! third-order

!!            vr(ni  ,m) = (35.0*vn(ni,m) - 35.0*vn(ni-1,m) + 21.0*vn(ni-2,m) - 5.0*vn(ni-3,m))/16.0
!!            vr(ni-1,m) = ( 5.0*vn(ni,m) + 15.0*vn(ni-1,m) -  5.0*vn(ni-2,m) +     vn(ni-3,m))/16.0
            vr(ni-2,m) = (    -vn(ni,m) +  9.0*vn(ni-1,m) +  9.0*vn(ni-2,m) -      vn(ni-3,m))/16.0
            vl(ni-2,m) = (-5.0*vn(ni,m) + 60.0*vn(ni-1,m) + 90.0*vn(ni-2,m) - 20.0*vn(ni-3,m) + 3.0*vn(ni-4,m) ) / 128.0   ! fifth-order

            vr(ni-3,m) = ( 3.0*vn(ni,m) - 20.0*vn(ni-1,m) + 90.0*vn(ni-2,m) + 60.0*vn(ni-3,m) - 5.0*vn(ni-4,m) ) / 128.0

            vl(ni  ,m) = vr(ni  ,m)
            vl(ni-1,m) = vr(ni-1,m)
        end do

        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vr(ni,m) < chklim(m) ) then
                    vr(ni,m) = vn(ni,m)
                    ndec = ndec + 1
                end if
            end if
        end do
    end if

    do m=nst,ned
        do i=st,ed+1
            g(1,i,m) = ( -2.0*vn(i-3,m) +  9.0*vn(i-2,m) - 18.0*vn(i-1,m) + 11.0*vn(i  ,m))/6.0
            g(2,i,m) = (      vn(i-2,m) -  6.0*vn(i-1,m) +  3.0*vn(i  ,m) +  2.0*vn(i+1,m))/6.0
            g(3,i,m) = ( -2.0*vn(i-1,m) -  3.0*vn(i  ,m) +  6.0*vn(i+1,m) -      vn(i+2,m))/6.0
            g(4,i,m) = (-11.0*vn(i  ,m) + 18.0*vn(i+1,m) -  9.0*vn(i+2,m) +  2.0*vn(i+3,m))/6.0

            s(1,i,m) =     -vn(i-3,m) + 4.0*vn(i-2,m) - 5.0*vn(i-1,m) + 2.0*vn(i  ,m)
            s(2,i,m) =                      vn(i-1,m) - 2.0*vn(i  ,m) +     vn(i+1,m) 
            s(3,i,m) =      vn(i-1,m) - 2.0*vn(i  ,m) +     vn(i+1,m)
            s(4,i,m) =  2.0*vn(i  ,m) - 5.0*vn(i+1,m) + 4.0*vn(i+2,m) -     vn(i+3,m)

            t(1,i,m) = -vn(i-3,m) + 3.0*vn(i-2,m) - 3.0*vn(i-1,m) + vn(i  ,m)
            t(2,i,m) = -vn(i-2,m) + 3.0*vn(i-1,m) - 3.0*vn(i  ,m) + vn(i+1,m)
            t(3,i,m) = -vn(i-1,m) + 3.0*vn(i  ,m) - 3.0*vn(i+1,m) + vn(i+2,m)
            t(4,i,m) = -vn(i  ,m) + 3.0*vn(i+1,m) - 3.0*vn(i+2,m) + vn(i+3,m)
        end do
    end do

    do i=st,ed+1
        nx = sn(i,1)
        ny = sn(i,2)
        nz = sn(i,3)

        do m=nst,ned
            vp(m) = vn(i,m)
        end do
        call LxQ_p(nst,ned,vp,nt,nx,ny,nz,vp,vc) 
        do m=nst,ned
            q(i,m) = vc(m)
        end do
         
        do n=1,4
            do m=nst,ned
                fp(m) = g(n,i,m)
            end do
            call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
            do m=nst,ned
                g(n,i,m) = vc(m)
            end do

            do m=nst,ned
                fp(m) = s(n,i,m)
            end do
            call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
            do m=nst,ned
                s(n,i,m) = vc(m)
            end do

            do m=nst,ned
                fp(m) = t(n,i,m)
            end do
            call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
            do m=nst,ned
                t(n,i,m) = vc(m)
            end do
        end do
    end do

    do m=nst,ned
        do i=st,ed+1
            do n=1,4
                is = g(n,i,m)*g(n,i,m) + s(n,i,m)*s(n,i,m) + t(n,i,m)*t(n,i,m)
                eis2 = (eps+is)**2
                bl(n) = cl(n)/eis2
                br(n) = cr(n)/eis2
            end do

            bs = bl(1) + bl(2) + bl(3) + bl(4)
            do n=1,4
                wl(n,i,m) = bl(n)/bs
            end do

            bs = br(1) + br(2) + br(3) + br(4)
            do n=1,4
                wr(n,i,m) = br(n)/bs
            end do
        end do
    end do

    do m=nst,ned
        do i=st,ed+1
            vl(i,m) = q(i ,m) + ( wl(1,i ,m)*(24.0*g(1,i ,m)+6.0*s(1,i ,m)+t(1,i ,m)) + &
                                  wl(2,i ,m)*(24.0*g(2,i ,m)+6.0*s(2,i ,m)+t(2,i ,m)) + &
                                  wl(3,i ,m)*(24.0*g(3,i ,m)+6.0*s(3,i ,m)+t(3,i ,m)) + &
                                  wl(4,i ,m)*(24.0*g(4,i ,m)+6.0*s(4,i ,m)+t(4,i ,m)) )/48.0
        end do

        do i=st-1,ed
            i1 = i+1
            vr(i,m) = q(i1,m) + ( wr(1,i1,m)*(-24.0*g(1,i1,m)+6.0*s(1,i1,m)-t(1,i1,m)) + &
                                  wr(2,i1,m)*(-24.0*g(2,i1,m)+6.0*s(2,i1,m)-t(2,i1,m)) + &
                                  wr(3,i1,m)*(-24.0*g(3,i1,m)+6.0*s(3,i1,m)-t(3,i1,m)) + &
                                  wr(4,i1,m)*(-24.0*g(4,i1,m)+6.0*s(4,i1,m)-t(4,i1,m)) )/48.0
        end do
    end do

    do i=st,ed+1
        nx = sn(i,1)
        ny = sn(i,2)
        nz = sn(i,3)

        do m=nst,ned
            vp(m) = vn(i,m)
            vc(m) = vl(i,m)
        end do

        call RxQ_p(nst,ned,vp,nt,nx,ny,nz,vc,fp) 
        do m=nst,ned
            vl(i,m) = fp(m)
        end do
    end do

    do i=st-1,ed
        i1 = i+1
        nx = sn(i1,1)
        ny = sn(i1,2)
        nz = sn(i1,3)

        do m=nst,ned
            vp(m) = vn(i1,m)
            vc(m) = vr(i ,m)
        end do

        call RxQ_p(nst,ned,vp,nt,nx,ny,nz,vc,fp) 
        do m=nst,ned
            vr(i,m) = fp(m)
        end do
    end do

    do i=stc,edc
        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vl(i,m) < chklim(m) ) then
                    vl(i,m) = vn(i,m)
                    ndec = ndec + 1
                end if

                if (vr(i,m) < chklim(m) ) then
                    vr(i,m) = vn(i+1,m)
                    ndec = ndec + 1
                end if
           end if
        end do
    end do

end subroutine wcns7cv

subroutine wcns7cvcc(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,one
    use mod_constants, only : nbc_inter_scheme
    use mod_variables, only : rlimit,plimit
    implicit none
    integer(kind_int), intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  , intent(inout)  :: vn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(in)  :: sn(1-ngn:ni+ngn,1:3)
    integer(kind_int), intent(in)  :: nfs,nfe
    integer(kind_int), intent(in)  :: nge
    real(kind_real)  , intent(out) :: vl(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(out) :: vr(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , external    :: sub_limit
    integer(kind_int), intent(in)  :: chkid(nst:ned)
    real(kind_real),   intent(in)  :: chklim(nst:ned)
    integer(kind_int), intent(out) :: ndec
    real(kind_real), parameter :: cl(4) = (/1.D0/64.D0,21.D0/64.D0,35.D0/64.D0,7.D0/64.D0/)
    real(kind_real), parameter :: cr(4) = (/7.D0/64.D0,35.D0/64.D0,21.D0/64.D0,1.D0/64.D0/)
    real(kind_real), parameter :: eps   = 1.0e-6
    real(kind_real)            :: is,eis2,bs,bl(4),br(4)
    real(kind_real)            :: nt,nx,ny,nz,fp(nst:ned)
    real(kind_real)            :: vp(nst:ned),vc(nst:ned)
    real(kind_real)            :: q(1-ngn:ni+ngn,nst:ned)
    real(kind_real)            :: s(1:4,1-ngn:ni+ngn,nst:ned)
    real(kind_real)            :: g(1:4,1-ngn:ni+ngn,nst:ned)
    real(kind_real)            :: t(1:4,1-ngn:ni+ngn,nst:ned)
    real(kind_real)            :: wl(1:4,1-ngn:ni+ngn,nst:ned)
    real(kind_real)            :: wr(1:4,1-ngn:ni+ngn,nst:ned)
    integer(kind_int)          :: i,i1,m,n,st,ed,stc,edc

    ndec = 0

    if (nfs == nbc_inter_scheme) then
        st = 1- nge
        stc = st
    else
        st = 5
        stc = 1

        do m=nst,ned
            vl(1,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                  ! third-order
            vl(2,m) = (  3.0*vn(1,m) +  6.0*vn(2,m) -      vn(3,m) ) / 8.0                  ! third-order

!!            vl(1,m) = (35.0*vn(1,m) - 35.0*vn(2,m) + 21.0*vn(3,m) - 5.0*vn(4,m))/16.0
!!            vl(2,m) = ( 5.0*vn(1,m) + 15.0*vn(2,m) -  5.0*vn(3,m) +     vn(4,m))/16.0
            vl(3,m) = (    -vn(1,m) +  9.0*vn(2,m) +  9.0*vn(3,m) -      vn(4,m))/16.0
            vr(3,m) = (-5.0*vn(1,m) + 60.0*vn(2,m) + 90.0*vn(3,m) - 20.0*vn(4,m) + 3.0*vn(5,m) ) / 128.0

            vl(4,m) = ( 3.0*vn(1,m) - 20.0*vn(2,m) + 90.0*vn(3,m) + 60.0*vn(4,m) - 5.0*vn(5,m) ) / 128.0

            vr(1,m) = vl(1,m)
            vr(2,m) = vl(2,m)
        end do

        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vl(1,m) < chklim(m) ) then
                    vl(1,m) = vn(1,m)
                    ndec = ndec + 1
                end if
            end if
        end do
    end if

    if (nfe == nbc_inter_scheme) then
        ed = ni + 1 + nge
        edc = ed
    else
        ed = ni - 3
        edc = ni + 1

        do m=nst,ned
            vr(ni+1,m) = ( 15.0*vn(ni,m) - 10.0*vn(ni-1,m) +  3.0*vn(ni-2,m) ) / 8.0                ! third-order
            vr(ni  ,m) = (  3.0*vn(ni,m) +  6.0*vn(ni-1,m) -      vn(ni-2,m) ) / 8.0                ! third-order

!!            vr(ni+1,m) = (35.0*vn(ni,m) - 35.0*vn(ni-1,m) + 21.0*vn(ni-2,m) - 5.0*vn(ni-3,m))/16.0
!!            vr(ni  ,m) = ( 5.0*vn(ni,m) + 15.0*vn(ni-1,m) -  5.0*vn(ni-2,m) +     vn(ni-3,m))/16.0
            vr(ni-1,m) = (    -vn(ni,m) +  9.0*vn(ni-1,m) +  9.0*vn(ni-2,m) -      vn(ni-3,m))/16.0
            vl(ni-1,m) = (-5.0*vn(ni,m) + 60.0*vn(ni-1,m) + 90.0*vn(ni-2,m) - 20.0*vn(ni-3,m) + 3.0*vn(ni-4,m) ) / 128.0   ! fifth-order

            vr(ni-2,m) = ( 3.0*vn(ni,m) - 20.0*vn(ni-1,m) + 90.0*vn(ni-2,m) + 60.0*vn(ni-3,m) - 5.0*vn(ni-4,m) ) / 128.0

            vl(ni+1,m) = vr(ni+1,m)
            vl(ni  ,m) = vr(ni  ,m)
        end do

        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vr(ni+1,m) < chklim(m) ) then
                    vr(ni+1,m) = vn(ni,m)
                    ndec = ndec + 1
                end if
            end if
        end do
    end if

    do m=nst,ned
        do i=st-1,ed
            g(1,i,m) = ( -2.0*vn(i-3,m) +  9.0*vn(i-2,m) - 18.0*vn(i-1,m) + 11.0*vn(i  ,m))/6.0
            g(2,i,m) = (      vn(i-2,m) -  6.0*vn(i-1,m) +  3.0*vn(i  ,m) +  2.0*vn(i+1,m))/6.0
            g(3,i,m) = ( -2.0*vn(i-1,m) -  3.0*vn(i  ,m) +  6.0*vn(i+1,m) -      vn(i+2,m))/6.0
            g(4,i,m) = (-11.0*vn(i  ,m) + 18.0*vn(i+1,m) -  9.0*vn(i+2,m) +  2.0*vn(i+3,m))/6.0

            s(1,i,m) =     -vn(i-3,m) + 4.0*vn(i-2,m) - 5.0*vn(i-1,m) + 2.0*vn(i  ,m)
            s(2,i,m) =                      vn(i-1,m) - 2.0*vn(i  ,m) +     vn(i+1,m) 
            s(3,i,m) =      vn(i-1,m) - 2.0*vn(i  ,m) +     vn(i+1,m)
            s(4,i,m) =  2.0*vn(i  ,m) - 5.0*vn(i+1,m) + 4.0*vn(i+2,m) -     vn(i+3,m)

            t(1,i,m) = -vn(i-3,m) + 3.0*vn(i-2,m) - 3.0*vn(i-1,m) + vn(i  ,m)
            t(2,i,m) = -vn(i-2,m) + 3.0*vn(i-1,m) - 3.0*vn(i  ,m) + vn(i+1,m)
            t(3,i,m) = -vn(i-1,m) + 3.0*vn(i  ,m) - 3.0*vn(i+1,m) + vn(i+2,m)
            t(4,i,m) = -vn(i  ,m) + 3.0*vn(i+1,m) - 3.0*vn(i+2,m) + vn(i+3,m)
        end do
    end do

    do i=st-1,ed
        nx = sn(i,1)
        ny = sn(i,2)
        nz = sn(i,3)

        do m=nst,ned
            vp(m) = vn(i,m)
        end do
        call LxQ_p(nst,ned,vp,nt,nx,ny,nz,vp,vc) 
        do m=nst,ned
            q(i,m) = vc(m)
        end do
         
        do n=1,4
            do m=nst,ned
                fp(m) = g(n,i,m)
            end do
            call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
            do m=nst,ned
                g(n,i,m) = vc(m)
            end do

            do m=nst,ned
                fp(m) = s(n,i,m)
            end do
            call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
            do m=nst,ned
                s(n,i,m) = vc(m)
            end do

            do m=nst,ned
                fp(m) = t(n,i,m)
            end do
            call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
            do m=nst,ned
                t(n,i,m) = vc(m)
            end do
        end do
    end do

    do m=nst,ned
        do i=st-1,ed
            do n=1,4
                is = g(n,i,m)*g(n,i,m) + s(n,i,m)*s(n,i,m) + t(n,i,m)*t(n,i,m)
                eis2 = (eps+is)**2
                bl(n) = cl(n)/eis2
                br(n) = cr(n)/eis2
            end do

            bs = bl(1) + bl(2) + bl(3) + bl(4)
            do n=1,4
                wl(n,i,m) = bl(n)/bs
            end do

            bs = br(1) + br(2) + br(3) + br(4)
            do n=1,4
                wr(n,i,m) = br(n)/bs
            end do
        end do
    end do

    do m=nst,ned
        do i=st-1,ed
            vl(i+1,m) = q(i ,m) + ( wl(1,i ,m)*(24.0*g(1,i ,m)+6.0*s(1,i ,m)+t(1,i ,m)) + &
                                    wl(2,i ,m)*(24.0*g(2,i ,m)+6.0*s(2,i ,m)+t(2,i ,m)) + &
                                    wl(3,i ,m)*(24.0*g(3,i ,m)+6.0*s(3,i ,m)+t(3,i ,m)) + &
                                    wl(4,i ,m)*(24.0*g(4,i ,m)+6.0*s(4,i ,m)+t(4,i ,m)) )/48.0
        end do

        do i=st-2,ed-1
            i1 = i+1
            vr(i+1,m) = q(i1,m) + ( wr(1,i1,m)*(-24.0*g(1,i1,m)+6.0*s(1,i1,m)-t(1,i1,m)) + &
                                    wr(2,i1,m)*(-24.0*g(2,i1,m)+6.0*s(2,i1,m)-t(2,i1,m)) + &
                                    wr(3,i1,m)*(-24.0*g(3,i1,m)+6.0*s(3,i1,m)-t(3,i1,m)) + &
                                    wr(4,i1,m)*(-24.0*g(4,i1,m)+6.0*s(4,i1,m)-t(4,i1,m)) )/48.0
        end do
    end do

    do i=st,ed+1
        nx = sn(i,1)
        ny = sn(i,2)
        nz = sn(i,3)

        do m=nst,ned
            vp(m) = vn(i,m)
            vc(m) = vl(i,m)
        end do

        call RxQ_p(nst,ned,vp,nt,nx,ny,nz,vc,fp) 
        do m=nst,ned
            vl(i,m) = fp(m)
        end do
    end do

    do i=st-1,ed
        i1 = i+1
        nx = sn(i1,1)
        ny = sn(i1,2)
        nz = sn(i1,3)

        do m=nst,ned
            vp(m) = vn(i1,m)
            vc(m) = vr(i ,m)
        end do

        call RxQ_p(nst,ned,vp,nt,nx,ny,nz,vc,fp) 
        do m=nst,ned
            vr(i,m) = fp(m)
        end do
    end do

    do i=stc,edc
        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vl(i,m) < chklim(m) ) then
                    vl(i,m) = vn(i,m)
                    ndec = ndec + 1
                end if

                if (vr(i,m) < chklim(m) ) then
                    vr(i,m) = vn(i+1,m)
                    ndec = ndec + 1
                end if
           end if
        end do
    end do

end subroutine wcns7cvcc

subroutine hdcs5ei(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_inter_scheme,nbc_intbc_scheme
    implicit none
    integer(kind_int), intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  , intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(in)  :: sn(1-ngn:ni+ngn,1:3)
    integer(kind_int), intent(in)  :: nfs,nfe
    integer(kind_int), intent(in)  :: nge
    real(kind_real)  , intent(out) :: vl(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(out) :: vr(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , external    :: sub_limit
    integer(kind_int), intent(in)  :: chkid(nst:ned)
    real(kind_real),   intent(in)  :: chklim(nst:ned)
    integer(kind_int), intent(out) :: ndec
    real(kind_real)   :: aa
    integer(kind_int) :: i,m,st,ed,stc,edc

    if (ni < 4) then
        call muscl2pv(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)
        return
    end if

    ndec = 0
    
    aa   = 0.01
   
    if (nfs == nbc_inter_scheme) then
        st  = -nge
        stc = st
    else if (nfs == nbc_intbc_scheme) then
        st  = -nge
        stc = st     
    else
        st  = 3
        stc = 0

        do m=nst,ned
            vl(0,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                  ! third-order
            vl(1,m) = (  3.0*vn(1,m) +  6.0*vn(2,m) -      vn(3,m) ) / 8.0   ! third-order 
            !!vl(1,m) = (  35.0*vn(1,m) + 140.0*vn(2,m) -  70.0*vn(3,m) +  28.0*vn(4,m) -  5.0*vn(5,m) ) / 128.0   ! fifth-order
            !!vl(2,m) = (  35.0*vn(2,m) + 140.0*vn(3,m) -  70.0*vn(4,m) +  28.0*vn(5,m) -  5.0*vn(6,m) ) / 128.0   ! fifth-order
            vl(2,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 + &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )
            vr(0,m) = vl(0,m)
            vr(1,m) = vl(1,m)
            !!vr(2,m) = vl(2,m)
            vr(2,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 - &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )            
        end do

    end if

    if (nfe == nbc_inter_scheme) then
        ed  = ni + nge
        edc = ed
    
    else if (nfe == nbc_intbc_scheme) then    
        ed  = ni + nge
        edc = ed    
    else
        ed  = ni - 3
        edc = ni

        do m=nst,ned
            vr(ni  ,m) = ( 15.0*vn(ni,m) - 10.0*vn(ni-1,m) +  3.0*vn(ni-2,m) ) / 8.0                ! third-order
            vl(ni-1,m) = (  3.0*vn(ni,m) +  6.0*vn(ni-1,m) -      vn(ni-2,m) ) / 8.0   ! third-order
            !!vl(ni-1,m) = (  35.0*vn(ni  ,  m) + 140.0*vn(ni-1,m) -  70.0*vn(ni-2,m) +  28.0*vn(ni-3,m) -  5.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            !!vl(ni-2,m) = (  35.0*vn(ni-1,  m) + 140.0*vn(ni-2,m) -  70.0*vn(ni-3,m) +  28.0*vn(ni-4,m) -  5.0*vn(ni-5,m) ) / 128.0   ! fifth-order
            vl(ni-2,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 + &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) )
            vl(ni  ,m) = vr(ni  ,m)
            vr(ni-1,m) = vl(ni-1,m)
            !!vr(ni-2,m) = vl(ni-2,m)
            vr(ni-2,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 - &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) )            
        end do

    end if



    do m=nst,ned
        do i=st,ed
            vl(i,m) =   ( 150.0*( vn(i+1,m) + vn(i,m  ) )          -  &
                           25.0*( vn(i+2,m) + vn(i-1,m) )          +  &
                            3.0*( vn(i+3,m) + vn(i-2,m) ) ) /256.0 -  &
                      aa*( 10.0*( vn(i+1,m) - vn(i,m  ) )          -  &
                            5.0*( vn(i+2,m) - vn(i-1,m) )          +  &
                                ( vn(i+3,m) - vn(i-2,m) ) ) 
            
            vr(i,m) =   ( 150.0*( vn(i+1,m) + vn(i,m  ) )          -  &
                           25.0*( vn(i+2,m) + vn(i-1,m) )          +  &
                            3.0*( vn(i+3,m) + vn(i-2,m) ) ) /256.0 +  &
                      aa*( 10.0*( vn(i+1,m) - vn(i,m  ) )          -  &
                            5.0*( vn(i+2,m) - vn(i-1,m) )          +  &
                                ( vn(i+3,m) - vn(i-2,m) ) ) 		                        
        end do
    end do

    do i=stc,edc
        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vl(i,m) < chklim(m) ) then
                    vl(i,m) = vn(i,m)
                    ndec = ndec + 1
                end if

                if (vr(i,m) < chklim(m) ) then
                    vr(i,m) = vn(i+1,m)
                    ndec = ndec + 1
                end if
           end if
        end do
    end do

end subroutine hdcs5ei

subroutine hdcs5eicc(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_inter_scheme,nbc_intbc_scheme
    implicit none
    integer(kind_int), intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  , intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(in)  :: sn(1-ngn:ni+ngn,1:3)
    integer(kind_int), intent(in)  :: nfs,nfe
    integer(kind_int), intent(in)  :: nge
    real(kind_real)  , intent(out) :: vl(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(out) :: vr(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , external    :: sub_limit
    integer(kind_int), intent(in)  :: chkid(nst:ned)
    real(kind_real),   intent(in)  :: chklim(nst:ned)
    integer(kind_int), intent(out) :: ndec
    real(kind_real)   :: aa
    integer(kind_int) :: i,m,st,ed,stc,edc

    if (ni < 4) then
        call muscl2pvcc(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)
        return
    end if

    ndec = 0
    
    aa   = 0.01
   
    if (nfs == nbc_inter_scheme) then
        st  = 1 - nge
        stc = st
    else if (nfs == nbc_intbc_scheme) then
        st  = 1 - nge
        stc = st     
    else
        st  = 4
        stc = 1

        do m=nst,ned
            vl(1,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                  ! third-order
            vl(2,m) = (  3.0*vn(1,m) +  6.0*vn(2,m) -      vn(3,m) ) / 8.0   ! third-order 
            !!vl(2,m) = (  35.0*vn(1,m) + 140.0*vn(2,m) -  70.0*vn(3,m) +  28.0*vn(4,m) -  5.0*vn(5,m) ) / 128.0   ! fifth-order
            !!vl(3,m) = (  35.0*vn(2,m) + 140.0*vn(3,m) -  70.0*vn(4,m) +  28.0*vn(5,m) -  5.0*vn(6,m) ) / 128.0   ! fifth-order
            vl(3,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 + &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )
            vr(1,m) = vl(1,m)
            vr(2,m) = vl(2,m)
            !!vr(3,m) = vl(3,m)
            vr(3,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 - &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )            
        end do

    end if

    if (nfe == nbc_inter_scheme) then
        ed  = ni + 1 + nge
        edc = ed
    
    else if (nfe == nbc_intbc_scheme) then    
        ed  = ni + 1 + nge
        edc = ed    
    else
        ed  = ni - 2
        edc = ni + 1

        do m=nst,ned
            vr(ni+1,m) = ( 15.0*vn(ni,m) - 10.0*vn(ni-1,m) +  3.0*vn(ni-2,m) ) / 8.0                ! third-order
            vl(ni  ,m) = (  3.0*vn(ni,m) +  6.0*vn(ni-1,m) -      vn(ni-2,m) ) / 8.0   ! third-order
            !!vl(ni  ,m) = (  35.0*vn(ni  ,  m) + 140.0*vn(ni-1,m) -  70.0*vn(ni-2,m) +  28.0*vn(ni-3,m) -  5.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            !!vl(ni-1,m) = (  35.0*vn(ni-1,  m) + 140.0*vn(ni-2,m) -  70.0*vn(ni-3,m) +  28.0*vn(ni-4,m) -  5.0*vn(ni-5,m) ) / 128.0   ! fifth-order
            vl(ni-1,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 + &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) )
            vl(ni+1,m) = vr(ni+1,m)
            vr(ni  ,m) = vl(ni  ,m)
            !!vr(ni-1,m) = vl(ni-1,m)
            vr(ni-1,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 - &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) )            
        end do

    end if



    do m=nst,ned
        do i=st-1,ed-1
            vl(i+1,m) =   ( 150.0*( vn(i+1,m) + vn(i,m  ) )          -  &
                             25.0*( vn(i+2,m) + vn(i-1,m) )          +  &
                              3.0*( vn(i+3,m) + vn(i-2,m) ) ) /256.0 -  &
                        aa*( 10.0*( vn(i+1,m) - vn(i,m  ) )          -  &
                              5.0*( vn(i+2,m) - vn(i-1,m) )          +  &
                                  ( vn(i+3,m) - vn(i-2,m) ) ) 
            
            vr(i+1,m) =   ( 150.0*( vn(i+1,m) + vn(i,m  ) )          -  &
                             25.0*( vn(i+2,m) + vn(i-1,m) )          +  &
                              3.0*( vn(i+3,m) + vn(i-2,m) ) ) /256.0 +  &
                        aa*( 10.0*( vn(i+1,m) - vn(i,m  ) )          -  &
                              5.0*( vn(i+2,m) - vn(i-1,m) )          +  &
                                  ( vn(i+3,m) - vn(i-2,m) ) ) 		                        
        end do
    end do

    do i=stc,edc
        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vl(i,m) < chklim(m) ) then
                    vl(i,m) = vn(i,m)
                    ndec = ndec + 1
                end if

                if (vr(i,m) < chklim(m) ) then
                    vr(i,m) = vn(i+1,m)
                    ndec = ndec + 1
                end if
           end if
        end do
    end do

end subroutine hdcs5eicc

subroutine hdcs7ci(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_bound_scheme
    implicit none
    integer(kind_int), intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  , intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(in)  :: sn(1-ngn:ni+ngn,1:3)
    integer(kind_int), intent(in)  :: nfs,nfe
    integer(kind_int), intent(in)  :: nge
    real(kind_real)  , intent(out) :: vl(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(out) :: vr(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , external    :: sub_limit
    integer(kind_int), intent(in)  :: chkid(nst:ned)
    real(kind_real),   intent(in)  :: chklim(nst:ned)
    integer(kind_int), intent(out) :: ndec
    real(kind_real)   :: ma(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵ����Aϵ��
    real(kind_real)   :: mb(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵ����Bϵ��
    real(kind_real)   :: mc(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵ����Cϵ��
    real(kind_real)   :: md(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵ��Ҳ�Dϵ��
    real(kind_real)   :: ms(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵļ�����     
    real(kind_real)   :: aa,bb,diss
    integer(kind_int) :: i,m,st,ed,stc,edc

    if (ni < 4) then
        call muscl2pv(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)
        return
    end if

    ndec = 0
    
    aa   = 0.005
       
    bb   = 0.3
    
    !>
    !!�����߽�HDCS��߽�������Ҳ��ֵ
    if (nfs == nbc_bound_scheme) then
        st = 3       !< ����nfs=1ʱ���ڵ��ֵ��ʽ�������ʼ����
        stc = 0      !< ����nfs=1ʱ����ֵ�������ʼ����
        ma(st-4)=0.0
        ma(st-3)=0.0
        ma(st-2)=0.0                           !< ���Խ�׷�ϵ����Aϵ������߽�����ֵ��
        ma(st-1)=0.0                           !< ���Խ�׷�ϵ����Aϵ�����ٽ���߽�����ֵ��
        mc(st-4)=0.0
        mc(st-3)=0.0
        mc(st-2)=0.0                           !< ���Խ�׷�ϵ����Cϵ������߽�����ֵ��
        mc(st-1)=0.0                           !< ���Խ�׷�ϵ����Cϵ�����ٽ���߽�����ֵ��

        do m=nst,ned
            vl(0,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                  ! third-order
            vl(1,m) = ( 3.0*vn(1,m) +   6.0*vn(2,m) -      vn(3,m) ) / 8.0                        !< ��������ֵ��߽�
                                                                                                  
            diss = aa*(    -vn(1,m) +   3.0*vn(2,m) -  3.0*vn(3,m) +      vn(4,m) )*10.0          !< ���׺�ɢ��
            vl(2,m) = (    -vn(1,m) +   9.0*vn(2,m) +  9.0*vn(3,m) -      vn(4,m) )/16.0 + diss   !< �Ľ�����ֵ+���׺�ɢ�����ٽ��߽磨i=2��
            vr(2,m) = (    -vn(1,m) +   9.0*vn(2,m) +  9.0*vn(3,m) -      vn(4,m) )/16.0 - diss   !< �Ľ��Ҳ��ֵ-���׺�ɢ�����ٽ��߽磨i=2��
                                                                                                  
            vr(0,m) = vl(0,m)
            vr(1,m) = vl(1,m)                                                                     !< �Ҳ��ֵ��߽�

        end do                                                                                   
    else
        st = -1    !< ����nfs=0������ȫ�ڵ�߽��ʽ����-1ʱ���ڵ��ֵ��ʽ�������ʼ����
        stc = -2   !< ����nfs=0������ȫ�ڵ�߽��ʽ����-1ʱ����ֵ�������ʼ����
        ma(st-1)=0.0
        mc(st-1)=0.0

        do m=nst,ned
            vl(stc,m) =   ( 1225.0*( vn(stc+1,m) + vn(stc,m  ) )          -  &
                             245.0*( vn(stc+2,m) + vn(stc-1,m) )          +  &
                              49.0*( vn(stc+3,m) + vn(stc-2,m) )          -  &
                               5.0*( vn(stc+4,m) + vn(stc-3,m) )) /2048.0 +  &
                        aa*( -35.0*( vn(stc+1,m) - vn(stc,m  ) )          +  &
                              21.0*( vn(stc+2,m) - vn(stc-1,m) )          -  &
                               7.0*( vn(stc+3,m) - vn(stc-2,m) )          +  &
                                   ( vn(stc+4,m) - vn(stc-3,m) )) 
            
            vr(stc,m) =   ( 1225.0*( vn(stc+1,m) + vn(stc,m  ) )          -  &
                             245.0*( vn(stc+2,m) + vn(stc-1,m) )          +  &
                              49.0*( vn(stc+3,m) + vn(stc-2,m) )          -  &
                               5.0*( vn(stc+4,m) + vn(stc-3,m) )) /2048.0 -  &
                        aa*( -35.0*( vn(stc+1,m) - vn(stc,m  ) )          +  &
                              21.0*( vn(stc+2,m) - vn(stc-1,m) )          -  &
                               7.0*( vn(stc+3,m) - vn(stc-2,m) )          +  &
                                   ( vn(stc+4,m) - vn(stc-3,m) ))    
        end do
    end if
    !>

    !>
    !!�����߽�HDCS�ұ߽�������Ҳ��ֵ
    if (nfe == nbc_bound_scheme) then
        ed = ni - 3      !< ����nfe=1ʱ���ڵ��ֵ��ʽ�������ֹ����
        edc = ni     !< ����nfe=1ʱ����ֵ�������ֹ����
        ma(ed+1)=0.0                           !< ���Խ�׷�ϵ����Aϵ�����ٽ��ұ߽�����ֵ��
        ma(ed+2)=0.0                           !< ���Խ�׷�ϵ����Aϵ�����ұ߽�����ֵ��
        ma(ed+3)=0.0
        ma(ed+4)=0.0
        mc(ed+1)=0.0                           !< ���Խ�׷�ϵ����Cϵ�����ٽ��ұ߽�����ֵ��
        mc(ed+2)=0.0                           !< ���Խ�׷�ϵ����Cϵ�����ұ߽�����ֵ��
        mc(ed+3)=0.0
        mc(ed+4)=0.0

        do m=nst,ned
            vr(ni  ,m) = ( 15.0*vn(ni,m) - 10.0*vn(ni-1,m) +  3.0*vn(ni-2,m) ) / 8.0                ! third-order
            vl(ni-1,m) = (  3.0*vn(ni,m) +   6.0*vn(ni-1,m) -      vn(ni-2,m) ) / 8.0   ! third-order         !< �����Ҳ��ֵ��߽�
                                                                                                              
            diss = aa*(         vn(ni,m) -   3.0*vn(ni-1,m) +  3.0*vn(ni-2,m) -      vn(ni-3,m) )*10.0        !< ���׺�ɢ��
            vl(ni-2,m) = (     -vn(ni,m) +   9.0*vn(ni-1,m) +  9.0*vn(ni-2,m) -      vn(ni-3,m) )/16.0 + diss !< �Ľ�����ֵ+���׺�ɢ�����ٽ��߽磨i=ni-2��
            vr(ni-2,m) = (     -vn(ni,m) +   9.0*vn(ni-1,m) +  9.0*vn(ni-2,m) -      vn(ni-3,m) )/16.0 - diss !< �Ľ��Ҳ��ֵ-���׺�ɢ�����ٽ��߽磨i=ni-2��
            
            vl(ni  ,m) = vr(ni  ,m)
            vr(ni-1,m) = vl(ni-1,m)                                                                           !< �Ҳ��ֵ�ұ߽�                               
        end do                                                                                               
    else
        ed = ni + 1    !< ����nfe=0������ȫ�ڵ�߽��ʽ����-1ʱ���ڵ��ֵ��ʽ�������ֹ����
        edc = ni + 2   !< ����nfe=0������ȫ�ڵ�߽��ʽ����-1ʱ����ֵ�������ֹ����
        ma(ed+1)=0.0
        mc(ed+1)=0.0

        do m=nst,ned
            vl(edc,m) =   ( 1225.0*( vn(edc+1,m) + vn(edc,m  ) )          -  &
                             245.0*( vn(edc+2,m) + vn(edc-1,m) )          +  &
                              49.0*( vn(edc+3,m) + vn(edc-2,m) )          -  &
                               5.0*( vn(edc+4,m) + vn(edc-3,m) )) /2048.0 +  &
                        aa*( -35.0*( vn(edc+1,m) - vn(edc,m  ) )          +  &
                              21.0*( vn(edc+2,m) - vn(edc-1,m) )          -  &
                               7.0*( vn(edc+3,m) - vn(edc-2,m) )          +  &
                                   ( vn(edc+4,m) - vn(edc-3,m) )) 
            
            vr(edc,m) =   ( 1225.0*( vn(edc+1,m) + vn(edc,m  ) )          -  &
                             245.0*( vn(edc+2,m) + vn(edc-1,m) )          +  &
                              49.0*( vn(edc+3,m) + vn(edc-2,m) )          -  &
                               5.0*( vn(edc+4,m) + vn(edc-3,m) )) /2048.0 -  &
                        aa*( -35.0*( vn(edc+1,m) - vn(edc,m  ) )          +  &
                              21.0*( vn(edc+2,m) - vn(edc-1,m) )          -  &
                               7.0*( vn(edc+3,m) - vn(edc-2,m) )          +  &
                                   ( vn(edc+4,m) - vn(edc-3,m) ))    
        end do
    end if
    !>

    !>
    !!�����߽�HDCS�����Ҳ��ڵ��ֵ
    do m=nst,ned
        do i=st,ed
            vl(i,m) =    ( 350.0 *( vn(i+1,m) + vn(i,m  ) )          +  &
                            35.0 *( vn(i+2,m) + vn(i-1,m) )          -  &
                                  ( vn(i+3,m) + vn(i-2,m) ) ) /448.0 -  &
                      bb*( 350.0 *( vn(i+1,m) - vn(i,m  ) )          +  &
                           105.0 *( vn(i+2,m) - vn(i-1,m) )          -  &
                             5.0 *( vn(i+3,m) - vn(i-2,m) ) ) /896.0            !< �߽�HDCS�ڵ��ֵ����ࣩ @attention vlֻ���߽�HDCS��ֵ��ʽ���Ҷ�  
            
            vr(i,m) =    ( 350.0 *( vn(i+1,m) + vn(i,m  ) )          +  &
                            35.0 *( vn(i+2,m) + vn(i-1,m) )          -  &
                                  ( vn(i+3,m) + vn(i-2,m) ) ) /448.0 +  &
                      bb*( 350.0 *( vn(i+1,m) - vn(i,m  ) )          +  &
                           105.0 *( vn(i+2,m) - vn(i-1,m) )          -  &
                             5.0 *( vn(i+3,m) - vn(i-2,m) ) ) /896.0            !< ���HDCS�ڵ��ֵ���Ҳࣩ @attention vrֻ���߽�HDCS��ֵ��ʽ���Ҷ�   	                        
        end do
    end do
    
    do i=stc,edc
        mb(i)=1.0                          !< ���Խ�׷�ϵ����Bϵ��
    enddo
        
  
    do i=st,ed
        ma(i)=5.0/14.0*(1.0+bb)            !< ���Խ�׷�ϵ����Aϵ�����ڵ�����ֵ��
    end do
            


    do i=st,ed                                                
        mc(i)=5.0/14.0*(1.0-bb)            !< ���Խ�׷�ϵ����Cϵ�����ڵ�����ֵ��
    end do                                               


    
    do m=nst,ned
        do i=stc,edc
            md(i) = vl(i,m)               !< ���Խ�׷�ϵ��Ҳ�Dϵ��������ֵ��
        end do
        call tri_pursue(ma,mb,mc,md,ms,stc,edc,1-ngn,ni+ngn) !< ���Խ�׷��
        do i=stc,edc
            vl(i,m) = ms(i)               !< �������Խ�׷�Ͻ⣨����ֵ��
        end do
    end do
    
       
    do i=st,ed
        ma(i)=5.0/14.0*(1.0-bb)           !< ���Խ�׷�ϵ����Aϵ�����ڵ��Ҳ��ֵ @attention �߽缰�ٽ��߽�ϵ�����䣩
    end do
     
    do i=st,ed
        mc(i)=5.0/14.0*(1.0+bb)           !< ���Խ�׷�ϵ����Cϵ�����ڵ��Ҳ��ֵ @attention �߽缰�ٽ��߽�ϵ�����䣩
    end do       

    do m=nst,ned        
        do i=stc,edc
            md(i) = vr(i,m)               !< ���Խ�׷�ϵ��Ҳ�Dϵ�����Ҳ��ֵ��
        end do
        call tri_pursue(ma,mb,mc,md,ms,stc,edc,1-ngn,ni+ngn) !< ���Խ�׷��
        do i=stc,edc
            vr(i,m) = ms(i)               !< �������Խ�׷�Ͻ⣨�Ҳ��ֵ��
        end do                      
    end do 

    !>
    !!��������Ҳ��ֵ�Ƿ����
    do i=stc,edc
        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vl(i,m) < chklim(m) ) then
                    vl(i,m) = vn(i,m)
                    ndec = ndec + 1
                end if

                if (vr(i,m) < chklim(m) ) then
                    vr(i,m) = vn(i+1,m)
                    ndec = ndec + 1
                end if
           end if
        end do
    end do
    !>
              
end subroutine hdcs7ci

subroutine hdcs7cicc(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_bound_scheme
    implicit none
    integer(kind_int), intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  , intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(in)  :: sn(1-ngn:ni+ngn,1:3)
    integer(kind_int), intent(in)  :: nfs,nfe
    integer(kind_int), intent(in)  :: nge
    real(kind_real)  , intent(out) :: vl(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(out) :: vr(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , external    :: sub_limit
    integer(kind_int), intent(in)  :: chkid(nst:ned)
    real(kind_real),   intent(in)  :: chklim(nst:ned)
    integer(kind_int), intent(out) :: ndec
    real(kind_real)   :: ma(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵ����Aϵ��
    real(kind_real)   :: mb(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵ����Bϵ��
    real(kind_real)   :: mc(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵ����Cϵ��
    real(kind_real)   :: md(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵ��Ҳ�Dϵ��
    real(kind_real)   :: ms(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵļ�����     
    real(kind_real)   :: aa,bb,diss
    integer(kind_int) :: i,m,st,ed,stc,edc

    if (ni < 4) then
        call muscl2pvcc(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)
        return
    end if

    ndec = 0
    
    aa   = 0.005
       
    bb   = 0.3
    
    !>
    !!�����߽�HDCS��߽�������Ҳ��ֵ
    if (nfs == nbc_bound_scheme) then
        st = 4       !< ����nfs=1ʱ���ڵ��ֵ��ʽ�������ʼ����
        stc = 1      !< ����nfs=1ʱ����ֵ�������ʼ����
        ma(st-4)=0.0
        ma(st-3)=0.0
        ma(st-2)=0.0                           !< ���Խ�׷�ϵ����Aϵ������߽�����ֵ��
        ma(st-1)=0.0                           !< ���Խ�׷�ϵ����Aϵ�����ٽ���߽�����ֵ��
        mc(st-4)=0.0
        mc(st-3)=0.0
        mc(st-2)=0.0                           !< ���Խ�׷�ϵ����Cϵ������߽�����ֵ��
        mc(st-1)=0.0                           !< ���Խ�׷�ϵ����Cϵ�����ٽ���߽�����ֵ��

        do m=nst,ned
            vl(1,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                  ! third-order
            vl(2,m) = ( 3.0*vn(1,m) +   6.0*vn(2,m) -      vn(3,m) ) / 8.0                        !< ��������ֵ��߽�
                                                                                                  
            diss = aa*(    -vn(1,m) +   3.0*vn(2,m) -  3.0*vn(3,m) +      vn(4,m) )*10.0          !< ���׺�ɢ��
            vl(3,m) = (    -vn(1,m) +   9.0*vn(2,m) +  9.0*vn(3,m) -      vn(4,m) )/16.0 + diss   !< �Ľ�����ֵ+���׺�ɢ�����ٽ��߽磨i=2��
            vr(3,m) = (    -vn(1,m) +   9.0*vn(2,m) +  9.0*vn(3,m) -      vn(4,m) )/16.0 - diss   !< �Ľ��Ҳ��ֵ-���׺�ɢ�����ٽ��߽磨i=2��
                                                                                                  
            vr(1,m) = vl(1,m)
            vr(2,m) = vl(2,m)                                                                     !< �Ҳ��ֵ��߽�

        end do                                                                                   
    else
        st = 2 - nge    !< ����nfs=0������ȫ�ڵ�߽��ʽ����-1ʱ���ڵ��ֵ��ʽ�������ʼ����
        stc = 1 - nge   !< ����nfs=0������ȫ�ڵ�߽��ʽ����-1ʱ����ֵ�������ʼ����
        ma(st-1)=0.0
        mc(st-1)=0.0

        stc = stc - 1
		do m=nst,ned
            vl(stc+1,m) =   ( 1225.0*( vn(stc+1,m) + vn(stc,m  ) )          -  &
                               245.0*( vn(stc+2,m) + vn(stc-1,m) )          +  &
                                49.0*( vn(stc+3,m) + vn(stc-2,m) )          -  &
                                 5.0*( vn(stc+4,m) + vn(stc-3,m) )) /2048.0 +  &
                          aa*( -35.0*( vn(stc+1,m) - vn(stc,m  ) )          +  &
                                21.0*( vn(stc+2,m) - vn(stc-1,m) )          -  &
                                 7.0*( vn(stc+3,m) - vn(stc-2,m) )          +  &
                                     ( vn(stc+4,m) - vn(stc-3,m) )) 
            
            vr(stc+1,m) =   ( 1225.0*( vn(stc+1,m) + vn(stc,m  ) )          -  &
                               245.0*( vn(stc+2,m) + vn(stc-1,m) )          +  &
                                49.0*( vn(stc+3,m) + vn(stc-2,m) )          -  &
                                 5.0*( vn(stc+4,m) + vn(stc-3,m) )) /2048.0 -  &
                          aa*( -35.0*( vn(stc+1,m) - vn(stc,m  ) )          +  &
                                21.0*( vn(stc+2,m) - vn(stc-1,m) )          -  &
                                 7.0*( vn(stc+3,m) - vn(stc-2,m) )          +  &
                                     ( vn(stc+4,m) - vn(stc-3,m) ))    
        end do
		stc = stc + 1
    end if
    !>

    !>
    !!�����߽�HDCS�ұ߽�������Ҳ��ֵ
    if (nfe == nbc_bound_scheme) then
        ed = ni - 2      !< ����nfe=1ʱ���ڵ��ֵ��ʽ�������ֹ����
        edc = ni + 1     !< ����nfe=1ʱ����ֵ�������ֹ����
        ma(ed+1)=0.0                           !< ���Խ�׷�ϵ����Aϵ�����ٽ��ұ߽�����ֵ��
        ma(ed+2)=0.0                           !< ���Խ�׷�ϵ����Aϵ�����ұ߽�����ֵ��
        ma(ed+3)=0.0
        ma(ed+4)=0.0
        mc(ed+1)=0.0                           !< ���Խ�׷�ϵ����Cϵ�����ٽ��ұ߽�����ֵ��
        mc(ed+2)=0.0                           !< ���Խ�׷�ϵ����Cϵ�����ұ߽�����ֵ��
        mc(ed+3)=0.0
        mc(ed+4)=0.0

        do m=nst,ned
            vr(ni+1,m) = ( 15.0*vn(ni,m) - 10.0*vn(ni-1,m) +  3.0*vn(ni-2,m) ) / 8.0                ! third-order
            vl(ni  ,m) = (  3.0*vn(ni,m) +   6.0*vn(ni-1,m) -      vn(ni-2,m) ) / 8.0   ! third-order         !< �����Ҳ��ֵ��߽�
                                                                                                              
            diss = aa*(         vn(ni,m) -   3.0*vn(ni-1,m) +  3.0*vn(ni-2,m) -      vn(ni-3,m) )*10.0        !< ���׺�ɢ��
            vl(ni-1,m) = (     -vn(ni,m) +   9.0*vn(ni-1,m) +  9.0*vn(ni-2,m) -      vn(ni-3,m) )/16.0 + diss !< �Ľ�����ֵ+���׺�ɢ�����ٽ��߽磨i=ni-2��
            vr(ni-1,m) = (     -vn(ni,m) +   9.0*vn(ni-1,m) +  9.0*vn(ni-2,m) -      vn(ni-3,m) )/16.0 - diss !< �Ľ��Ҳ��ֵ-���׺�ɢ�����ٽ��߽磨i=ni-2��
            
            vl(ni+1,m) = vr(ni+1,m)
            vr(ni  ,m) = vl(ni  ,m)                                                                           !< �Ҳ��ֵ�ұ߽�                               
        end do                                                                                               
    else
        ed = ni + nge    !< ����nfe=0������ȫ�ڵ�߽��ʽ����-1ʱ���ڵ��ֵ��ʽ�������ֹ����
        edc = ni + 1 + nge   !< ����nfe=0������ȫ�ڵ�߽��ʽ����-1ʱ����ֵ�������ֹ����
        ma(ed+1)=0.0
        mc(ed+1)=0.0

        edc = edc -1
		do m=nst,ned
            vl(edc+1,m) =   ( 1225.0*( vn(edc+1,m) + vn(edc,m  ) )          -  &
                               245.0*( vn(edc+2,m) + vn(edc-1,m) )          +  &
                                49.0*( vn(edc+3,m) + vn(edc-2,m) )          -  &
                                 5.0*( vn(edc+4,m) + vn(edc-3,m) )) /2048.0 +  &
                          aa*( -35.0*( vn(edc+1,m) - vn(edc,m  ) )          +  &
                                21.0*( vn(edc+2,m) - vn(edc-1,m) )          -  &
                                 7.0*( vn(edc+3,m) - vn(edc-2,m) )          +  &
                                     ( vn(edc+4,m) - vn(edc-3,m) )) 
            
            vr(edc+1,m) =   ( 1225.0*( vn(edc+1,m) + vn(edc,m  ) )          -  &
                               245.0*( vn(edc+2,m) + vn(edc-1,m) )          +  &
                                49.0*( vn(edc+3,m) + vn(edc-2,m) )          -  &
                                 5.0*( vn(edc+4,m) + vn(edc-3,m) )) /2048.0 -  &
                          aa*( -35.0*( vn(edc+1,m) - vn(edc,m  ) )          +  &
                                21.0*( vn(edc+2,m) - vn(edc-1,m) )          -  &
                                 7.0*( vn(edc+3,m) - vn(edc-2,m) )          +  &
                                     ( vn(edc+4,m) - vn(edc-3,m) ))    
        end do
		edc = edc + 1
    end if
    !>

    !>
    !!�����߽�HDCS�����Ҳ��ڵ��ֵ
    do m=nst,ned
        do i=st-1,ed-1
            vl(i+1,m) =    ( 350.0 *( vn(i+1,m) + vn(i,m  ) )          +  &
                              35.0 *( vn(i+2,m) + vn(i-1,m) )          -  &
                                    ( vn(i+3,m) + vn(i-2,m) ) ) /448.0 -  &
                        bb*( 350.0 *( vn(i+1,m) - vn(i,m  ) )          +  &
                             105.0 *( vn(i+2,m) - vn(i-1,m) )          -  &
                               5.0 *( vn(i+3,m) - vn(i-2,m) ) ) /896.0            !< �߽�HDCS�ڵ��ֵ����ࣩ @attention vlֻ���߽�HDCS��ֵ��ʽ���Ҷ�  
            
            vr(i+1,m) =    ( 350.0 *( vn(i+1,m) + vn(i,m  ) )          +  &
                              35.0 *( vn(i+2,m) + vn(i-1,m) )          -  &
                                    ( vn(i+3,m) + vn(i-2,m) ) ) /448.0 +  &
                        bb*( 350.0 *( vn(i+1,m) - vn(i,m  ) )          +  &
                             105.0 *( vn(i+2,m) - vn(i-1,m) )          -  &
                               5.0 *( vn(i+3,m) - vn(i-2,m) ) ) /896.0            !< ���HDCS�ڵ��ֵ���Ҳࣩ @attention vrֻ���߽�HDCS��ֵ��ʽ���Ҷ�   	                        
        end do
    end do
    
    do i=stc,edc
        mb(i)=1.0                          !< ���Խ�׷�ϵ����Bϵ��
    enddo
        
  
    do i=st,ed
        ma(i)=5.0/14.0*(1.0+bb)            !< ���Խ�׷�ϵ����Aϵ�����ڵ�����ֵ��
    end do
            


    do i=st,ed                                                
        mc(i)=5.0/14.0*(1.0-bb)            !< ���Խ�׷�ϵ����Cϵ�����ڵ�����ֵ��
    end do                                               


    
    do m=nst,ned
        do i=stc,edc
            md(i) = vl(i,m)               !< ���Խ�׷�ϵ��Ҳ�Dϵ��������ֵ��
        end do
        call tri_pursue(ma,mb,mc,md,ms,stc,edc,1-ngn,ni+ngn) !< ���Խ�׷��
        do i=stc,edc
            vl(i,m) = ms(i)               !< �������Խ�׷�Ͻ⣨����ֵ��
        end do
    end do
    
       
    do i=st,ed
        ma(i)=5.0/14.0*(1.0-bb)           !< ���Խ�׷�ϵ����Aϵ�����ڵ��Ҳ��ֵ @attention �߽缰�ٽ��߽�ϵ�����䣩
    end do
     
    do i=st,ed
        mc(i)=5.0/14.0*(1.0+bb)           !< ���Խ�׷�ϵ����Cϵ�����ڵ��Ҳ��ֵ @attention �߽缰�ٽ��߽�ϵ�����䣩
    end do       

    do m=nst,ned        
        do i=stc,edc
            md(i) = vr(i,m)               !< ���Խ�׷�ϵ��Ҳ�Dϵ�����Ҳ��ֵ��
        end do
        call tri_pursue(ma,mb,mc,md,ms,stc,edc,1-ngn,ni+ngn) !< ���Խ�׷��
        do i=stc,edc
            vr(i,m) = ms(i)               !< �������Խ�׷�Ͻ⣨�Ҳ��ֵ��
        end do                      
    end do 

    !>
    !!��������Ҳ��ֵ�Ƿ����
    do i=stc,edc
        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vl(i,m) < chklim(m) ) then
                    vl(i,m) = vn(i,m)
                    ndec = ndec + 1
                end if

                if (vr(i,m) < chklim(m) ) then
                    vr(i,m) = vn(i+1,m)
                    ndec = ndec + 1
                end if
           end if
        end do
    end do
    !>
              
end subroutine hdcs7cicc

subroutine hdcs7ciN(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec,st,ed)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_bound_scheme
    implicit none
    integer(kind_int), intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  , intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(in)  :: sn(1-ngn:ni+ngn,1:3)
    integer(kind_int), intent(in)  :: nfs,nfe
    integer(kind_int), intent(in)  :: nge,st,ed
    real(kind_real)  , intent(out) :: vl(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(out) :: vr(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , external    :: sub_limit
    integer(kind_int), intent(in)  :: chkid(nst:ned)
    real(kind_real),   intent(in)  :: chklim(nst:ned)
    integer(kind_int), intent(out) :: ndec
    real(kind_real)   :: ma(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵ����Aϵ��
    real(kind_real)   :: mb(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵ����Bϵ��
    real(kind_real)   :: mc(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵ����Cϵ��
    real(kind_real)   :: md(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵ��Ҳ�Dϵ��
    real(kind_real)   :: ms(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵļ�����     
    real(kind_real)   :: aa,bb
    integer(kind_int) :: i,m

    ndec = 0
    
    aa   = 0.005
       
    bb   = 0.3

    ma(st)=0.0
    mc(st)=0.0    
    do m=nst,ned
        vl(st,m) =   ( 1225.0*( vn(st+1,m) + vn(st,m  ) )          -  &
                        245.0*( vn(st+2,m) + vn(st-1,m) )          +  &
                         49.0*( vn(st+3,m) + vn(st-2,m) )          -  &
                          5.0*( vn(st+4,m) + vn(st-3,m) )) /2048.0 +  &
                   aa*( -35.0*( vn(st+1,m) - vn(st,m  ) )          +  &
                         21.0*( vn(st+2,m) - vn(st-1,m) )          -  &
                          7.0*( vn(st+3,m) - vn(st-2,m) )          +  &
                              ( vn(st+4,m) - vn(st-3,m) )) 
        
        vr(st,m) =   ( 1225.0*( vn(st+1,m) + vn(st,m  ) )          -  &
                        245.0*( vn(st+2,m) + vn(st-1,m) )          +  &
                         49.0*( vn(st+3,m) + vn(st-2,m) )          -  &
                          5.0*( vn(st+4,m) + vn(st-3,m) )) /2048.0 -  &
                   aa*( -35.0*( vn(st+1,m) - vn(st,m  ) )          +  &
                         21.0*( vn(st+2,m) - vn(st-1,m) )          -  &
                          7.0*( vn(st+3,m) - vn(st-2,m) )          +  &
                              ( vn(st+4,m) - vn(st-3,m) )) 
    end do

    ma(ed)=0.0
    mc(ed)=0.0    
    do m=nst,ned
        vl(ed,m) =   ( 1225.0*( vn(ed+1,m) + vn(ed,m  ) )          -  &
                        245.0*( vn(ed+2,m) + vn(ed-1,m) )          +  &
                         49.0*( vn(ed+3,m) + vn(ed-2,m) )          -  &
                          5.0*( vn(ed+4,m) + vn(ed-3,m) )) /2048.0 +  &
                   aa*( -35.0*( vn(ed+1,m) - vn(ed,m  ) )          +  &
                         21.0*( vn(ed+2,m) - vn(ed-1,m) )          -  &
                          7.0*( vn(ed+3,m) - vn(ed-2,m) )          +  &
                              ( vn(ed+4,m) - vn(ed-3,m) )) 
        
        vr(ed,m) =   ( 1225.0*( vn(ed+1,m) + vn(ed,m  ) )          -  &
                        245.0*( vn(ed+2,m) + vn(ed-1,m) )          +  &
                         49.0*( vn(ed+3,m) + vn(ed-2,m) )          -  &
                          5.0*( vn(ed+4,m) + vn(ed-3,m) )) /2048.0 -  &
                   aa*( -35.0*( vn(ed+1,m) - vn(ed,m  ) )          +  &
                         21.0*( vn(ed+2,m) - vn(ed-1,m) )          -  &
                          7.0*( vn(ed+3,m) - vn(ed-2,m) )          +  &
                              ( vn(ed+4,m) - vn(ed-3,m) )) 
    end do

    do m=nst,ned
        do i=st+1,ed-1
            vl(i,m) =    ( 350.0 *( vn(i+1,m) + vn(i,m  ) )          +  &
                            35.0 *( vn(i+2,m) + vn(i-1,m) )          -  &
                                  ( vn(i+3,m) + vn(i-2,m) ) ) /448.0 -  &
                      bb*( 350.0 *( vn(i+1,m) - vn(i,m  ) )          +  &
                           105.0 *( vn(i+2,m) - vn(i-1,m) )          -  &
                             5.0 *( vn(i+3,m) - vn(i-2,m) ) ) /896.0            !< �߽�HDCS�ڵ��ֵ����ࣩ @attention vlֻ���߽�HDCS��ֵ��ʽ���Ҷ�  
            
            vr(i,m) =    ( 350.0 *( vn(i+1,m) + vn(i,m  ) )          +  &
                            35.0 *( vn(i+2,m) + vn(i-1,m) )          -  &
                                  ( vn(i+3,m) + vn(i-2,m) ) ) /448.0 +  &
                      bb*( 350.0 *( vn(i+1,m) - vn(i,m  ) )          +  &
                           105.0 *( vn(i+2,m) - vn(i-1,m) )          -  &
                             5.0 *( vn(i+3,m) - vn(i-2,m) ) ) /896.0            !< ���HDCS�ڵ��ֵ���Ҳࣩ @attention vrֻ���߽�HDCS��ֵ��ʽ���Ҷ�   	                        
        end do
    end do
    
    do i=st,ed
        mb(i)=1.0                          !< ���Խ�׷�ϵ����Bϵ��
    enddo
        
  
    do i=st+1,ed-1
        ma(i)=5.0/14.0*(1.0+bb)            !< ���Խ�׷�ϵ����Aϵ�����ڵ�����ֵ��
    end do
            


    do i=st+1,ed-1
        mc(i)=5.0/14.0*(1.0-bb)            !< ���Խ�׷�ϵ����Cϵ�����ڵ�����ֵ��
    end do                                               


    
    do m=nst,ned
        do i=st,ed
            md(i) = vl(i,m)               !< ���Խ�׷�ϵ��Ҳ�Dϵ��������ֵ��
        end do
        call tri_pursue(ma,mb,mc,md,ms,st,ed,1-ngn,ni+ngn) !< ���Խ�׷��
        do i=st,ed
            vl(i,m) = ms(i)               !< �������Խ�׷�Ͻ⣨����ֵ��
        end do
    end do
    
       
    do i=st+1,ed-1
        ma(i)=5.0/14.0*(1.0-bb)           !< ���Խ�׷�ϵ����Aϵ�����ڵ��Ҳ��ֵ @attention �߽缰�ٽ��߽�ϵ�����䣩
    end do
     
    do i=st+1,ed-1
        mc(i)=5.0/14.0*(1.0+bb)           !< ���Խ�׷�ϵ����Cϵ�����ڵ��Ҳ��ֵ @attention �߽缰�ٽ��߽�ϵ�����䣩
    end do       

    do m=nst,ned        
        do i=st,ed
            md(i) = vr(i,m)               !< ���Խ�׷�ϵ��Ҳ�Dϵ�����Ҳ��ֵ��
        end do
        call tri_pursue(ma,mb,mc,md,ms,st,ed,1-ngn,ni+ngn) !< ���Խ�׷��
        do i=st,ed
            vr(i,m) = ms(i)               !< �������Խ�׷�Ͻ⣨�Ҳ��ֵ��
        end do                      
    end do 

    !>
    !!��������Ҳ��ֵ�Ƿ����
    do i=st,ed
        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vl(i,m) < chklim(m) ) then
                    vl(i,m) = vn(i,m)
                    ndec = ndec + 1
                end if

                if (vr(i,m) < chklim(m) ) then
                    vr(i,m) = vn(i+1,m)
                    ndec = ndec + 1
                end if
           end if
        end do
    end do
    !>
              
end subroutine hdcs7ciN

subroutine hdcs5ci(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_bound_scheme
    implicit none
    integer(kind_int), intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  , intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(in)  :: sn(1-ngn:ni+ngn,1:3)
    integer(kind_int), intent(in)  :: nfs,nfe
    integer(kind_int), intent(in)  :: nge
    real(kind_real)  , intent(out) :: vl(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(out) :: vr(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , external    :: sub_limit
    integer(kind_int), intent(in)  :: chkid(nst:ned)
    real(kind_real),   intent(in)  :: chklim(nst:ned)
    integer(kind_int), intent(out) :: ndec
    real(kind_real)   :: ma(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵ����Aϵ��
    real(kind_real)   :: mb(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵ����Bϵ��
    real(kind_real)   :: mc(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵ����Cϵ��
    real(kind_real)   :: md(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵ��Ҳ�Dϵ��
    real(kind_real)   :: ms(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵļ�����     
    real(kind_real)   :: aa,bb
    integer(kind_int) :: i,m,st,ed,stc,edc

    if (ni < 4) then
        call muscl2pv(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)
        return
    end if

    ndec = 0
    
    aa   = 0.01
       
    bb   = 0.3
    
    !>
    !!�����߽�HDCS��߽�������Ҳ��ֵ
    if (nfs == nbc_bound_scheme) then
        st = 2       !< ����nfs=1ʱ���ڵ��ֵ��ʽ�������ʼ����
        stc = 0      !< ����nfs=1ʱ����ֵ�������ʼ����
        ma(st-3)=0.0
        ma(st-2)=0.0                           !< ���Խ�׷�ϵ����Aϵ������߽�����ֵ��
        ma(st-1)=0.0                           !< ���Խ�׷�ϵ����Aϵ�����ٽ���߽�����ֵ��
        mc(st-3)=0.0
        mc(st-2)=0.0                           !< ���Խ�׷�ϵ����Cϵ������߽�����ֵ��
        mc(st-1)=0.0                           !< ���Խ�׷�ϵ����Cϵ�����ٽ���߽�����ֵ��

        do m=nst,ned
            vl(0,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                  ! third-order
            vl(1,m) = ( 3.0*vn(1,m) +   6.0*vn(2,m) -      vn(3,m) ) / 8.0                        !< ��������ֵ��߽�                                                                                                  
            vr(0,m) = vl(0,m)
            vr(1,m) = vl(1,m)                                                                     !< �Ҳ��ֵ��߽�
        end do                                                                                   
    else
        st = -2    !< ����nfs=0������ȫ�ڵ�߽��ʽ����-1ʱ���ڵ��ֵ��ʽ�������ʼ����
        stc = -3   !< ����nfs=0������ȫ�ڵ�߽��ʽ����-1ʱ����ֵ�������ʼ����
        ma(st-1)=0.0
        mc(st-1)=0.0

        do m=nst,ned
            vl(stc,m) =   ( 150.0*( vn(stc+1,m) + vn(stc,m  ) )          -  &
                             25.0*( vn(stc+2,m) + vn(stc-1,m) )          +  &
                              3.0*( vn(stc+3,m) + vn(stc-2,m) ) ) /256.0 -  &
                        aa*( 10.0*( vn(stc+1,m) - vn(stc,m  ) )          -  &
                              5.0*( vn(stc+2,m) - vn(stc-1,m) )          +  &
                                  ( vn(stc+3,m) - vn(stc-2,m) ) ) 
            
            vr(stc,m) =   ( 150.0*( vn(stc+1,m) + vn(stc,m  ) )          -  &
                             25.0*( vn(stc+2,m) + vn(stc-1,m) )          +  &
                              3.0*( vn(stc+3,m) + vn(stc-2,m) ) ) /256.0 +  &
                        aa*( 10.0*( vn(stc+1,m) - vn(stc,m  ) )          -  &
                              5.0*( vn(stc+2,m) - vn(stc-1,m) )          +  &
                                  ( vn(stc+3,m) - vn(stc-2,m) ) )      
        end do

    end if
    !>

    !>
    !!�����߽�HDCS�ұ߽�������Ҳ��ֵ
    if (nfe == nbc_bound_scheme) then
        ed = ni - 2      !< ����nfe=1ʱ���ڵ��ֵ��ʽ�������ֹ����
        edc = ni     !< ����nfe=1ʱ����ֵ�������ֹ����
        ma(ed+1)=0.0                           !< ���Խ�׷�ϵ����Aϵ�����ٽ��ұ߽�����ֵ��
        ma(ed+2)=0.0                           !< ���Խ�׷�ϵ����Aϵ�����ұ߽�����ֵ��
        ma(ed+3)=0.0
        mc(ed+1)=0.0                           !< ���Խ�׷�ϵ����Cϵ�����ٽ��ұ߽�����ֵ��
        mc(ed+2)=0.0                           !< ���Խ�׷�ϵ����Cϵ�����ұ߽�����ֵ��
        mc(ed+3)=0.0

        do m=nst,ned
            vr(ni  ,m) = ( 15.0*vn(ni,m) - 10.0*vn(ni-1,m) +  3.0*vn(ni-2,m) ) / 8.0                ! third-order
            vl(ni-1,m) = (  3.0*vn(ni,m) +   6.0*vn(ni-1,m) -      vn(ni-2,m) ) / 8.0   ! third-order         !< �����Ҳ��ֵ��߽�
            vl(ni  ,m) = vr(ni  ,m)
            vr(ni-1,m) = vl(ni-1,m)                                                                           !< �Ҳ��ֵ�ұ߽�                               
        end do                                                                                               
    else
        ed = ni + 2    !< ����nfe=0������ȫ�ڵ�߽��ʽ����-1ʱ���ڵ��ֵ��ʽ�������ֹ����
        edc = ni + 3   !< ����nfe=0������ȫ�ڵ�߽��ʽ����-1ʱ����ֵ�������ֹ����
        ma(ed+1)=0.0
        mc(ed+1)=0.0

        do m=nst,ned
            vl(edc,m) =   ( 150.0*( vn(edc+1,m) + vn(edc,m  ) )          -  &
                             25.0*( vn(edc+2,m) + vn(edc-1,m) )          +  &
                              3.0*( vn(edc+3,m) + vn(edc-2,m) ) ) /256.0 -  &
                        aa*( 10.0*( vn(edc+1,m) - vn(edc,m  ) )          -  &
                              5.0*( vn(edc+2,m) - vn(edc-1,m) )          +  &
                                  ( vn(edc+3,m) - vn(edc-2,m) ) ) 
            
            vr(edc,m) =   ( 150.0*( vn(edc+1,m) + vn(edc,m  ) )          -  &
                             25.0*( vn(edc+2,m) + vn(edc-1,m) )          +  &
                              3.0*( vn(edc+3,m) + vn(edc-2,m) ) ) /256.0 +  &
                        aa*( 10.0*( vn(edc+1,m) - vn(edc,m  ) )          -  &
                              5.0*( vn(edc+2,m) - vn(edc-1,m) )          +  &
                                  ( vn(edc+3,m) - vn(edc-2,m) ) )    
        end do
    end if
    !>

    !>
    !!�����߽�HDCS�����Ҳ��ڵ��ֵ
    do m=nst,ned
        do i=st,ed
            vl(i,m) =    (  15.0 *( vn(i+1,m) + vn(i,m  ) )          +  &
                                  ( vn(i+2,m) + vn(i-1,m) ) ) /20.0  -  &
                      bb*(  15.0 *( vn(i+1,m) - vn(i,m  ) )          +  &
                             3.0 *( vn(i+2,m) - vn(i-1,m) ) ) /40.0            !< �߽�HDCS�ڵ��ֵ����ࣩ @attention vlֻ���߽�HDCS��ֵ��ʽ���Ҷ�  
            
            vr(i,m) =    (  15.0 *( vn(i+1,m) + vn(i,m  ) )          +  &
                                  ( vn(i+2,m) + vn(i-1,m) ) ) /20.0  +  &
                      bb*(  15.0 *( vn(i+1,m) - vn(i,m  ) )          +  &
                             3.0 *( vn(i+2,m) - vn(i-1,m) ) ) /40.0            !< ���HDCS�ڵ��ֵ���Ҳࣩ @attention vrֻ���߽�HDCS��ֵ��ʽ���Ҷ�   	                        
        end do
    end do
    
    do i=stc,edc
        mb(i)=1.0                          !< ���Խ�׷�ϵ����Bϵ��
    enddo
        
  
    do i=st,ed
        ma(i)=3.0/10.0*(1.0+bb)            !< ���Խ�׷�ϵ����Aϵ�����ڵ�����ֵ��
    end do
            


    do i=st,ed                                                
        mc(i)=3.0/10.0*(1.0-bb)            !< ���Խ�׷�ϵ����Cϵ�����ڵ�����ֵ��
    end do                                               


    
    do m=nst,ned
        do i=stc,edc
            md(i) = vl(i,m)               !< ���Խ�׷�ϵ��Ҳ�Dϵ��������ֵ��
        end do
        call tri_pursue(ma,mb,mc,md,ms,stc,edc,1-ngn,ni+ngn) !< ���Խ�׷��
        do i=stc,edc
            vl(i,m) = ms(i)               !< �������Խ�׷�Ͻ⣨����ֵ��
        end do
    end do
    
       
    do i=st,ed
        ma(i)=3.0/10.0*(1.0-bb)           !< ���Խ�׷�ϵ����Aϵ�����ڵ��Ҳ��ֵ @attention �߽缰�ٽ��߽�ϵ�����䣩
    end do
     
    do i=st,ed
        mc(i)=3.0/10.0*(1.0+bb)           !< ���Խ�׷�ϵ����Cϵ�����ڵ��Ҳ��ֵ @attention �߽缰�ٽ��߽�ϵ�����䣩
    end do       

    do m=nst,ned        
        do i=stc,edc
            md(i) = vr(i,m)               !< ���Խ�׷�ϵ��Ҳ�Dϵ�����Ҳ��ֵ��
        end do
        call tri_pursue(ma,mb,mc,md,ms,stc,edc,1-ngn,ni+ngn) !< ���Խ�׷��
        do i=stc,edc
            vr(i,m) = ms(i)               !< �������Խ�׷�Ͻ⣨�Ҳ��ֵ��
        end do                      
    end do 

    !>
    !!��������Ҳ��ֵ�Ƿ����
    do i=stc,edc
        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vl(i,m) < chklim(m) ) then
                    vl(i,m) = vn(i,m)
                    ndec = ndec + 1
                end if

                if (vr(i,m) < chklim(m) ) then
                    vr(i,m) = vn(i+1,m)
                    ndec = ndec + 1
                end if
           end if
        end do
    end do
    !>
              
end subroutine hdcs5ci

subroutine hdcs5cicc(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_bound_scheme
    implicit none
    integer(kind_int), intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  , intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(in)  :: sn(1-ngn:ni+ngn,1:3)
    integer(kind_int), intent(in)  :: nfs,nfe
    integer(kind_int), intent(in)  :: nge
    real(kind_real)  , intent(out) :: vl(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(out) :: vr(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , external    :: sub_limit
    integer(kind_int), intent(in)  :: chkid(nst:ned)
    real(kind_real),   intent(in)  :: chklim(nst:ned)
    integer(kind_int), intent(out) :: ndec
    real(kind_real)   :: ma(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵ����Aϵ��
    real(kind_real)   :: mb(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵ����Bϵ��
    real(kind_real)   :: mc(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵ����Cϵ��
    real(kind_real)   :: md(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵ��Ҳ�Dϵ��
    real(kind_real)   :: ms(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵļ�����     
    real(kind_real)   :: aa,bb
    integer(kind_int) :: i,m,st,ed,stc,edc

    if (ni < 4) then
        call muscl2pvcc(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)
        return
    end if

    ndec = 0
    
    aa   = 0.01
       
    bb   = 0.3
    
    !>
    !!�����߽�HDCS��߽�������Ҳ��ֵ
    if (nfs == nbc_bound_scheme) then
        st = 3       !< ����nfs=1ʱ���ڵ��ֵ��ʽ�������ʼ����
        stc = 1      !< ����nfs=1ʱ����ֵ�������ʼ����
        ma(st-3)=0.0
        ma(st-2)=0.0                           !< ���Խ�׷�ϵ����Aϵ������߽�����ֵ��
        ma(st-1)=0.0                           !< ���Խ�׷�ϵ����Aϵ�����ٽ���߽�����ֵ��
        mc(st-3)=0.0
        mc(st-2)=0.0                           !< ���Խ�׷�ϵ����Cϵ������߽�����ֵ��
        mc(st-1)=0.0                           !< ���Խ�׷�ϵ����Cϵ�����ٽ���߽�����ֵ��

        do m=nst,ned
            vl(1,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                  ! third-order
            vl(2,m) = ( 3.0*vn(1,m) +   6.0*vn(2,m) -      vn(3,m) ) / 8.0                        !< ��������ֵ��߽�                                                                                                  
            vr(1,m) = vl(1,m)
            vr(2,m) = vl(2,m)                                                                     !< �Ҳ��ֵ��߽�
        end do                                                                                   
    else
        st = -1    !< ����nfs=0������ȫ�ڵ�߽��ʽ����-1ʱ���ڵ��ֵ��ʽ�������ʼ����
        stc = -2   !< ����nfs=0������ȫ�ڵ�߽��ʽ����-1ʱ����ֵ�������ʼ����
        ma(st-1)=0.0
        mc(st-1)=0.0

        stc=stc-1
		do m=nst,ned
            vl(stc+1,m) =   ( 150.0*( vn(stc+1,m) + vn(stc,m  ) )          -  &
                               25.0*( vn(stc+2,m) + vn(stc-1,m) )          +  &
                                3.0*( vn(stc+3,m) + vn(stc-2,m) ) ) /256.0 -  &
                          aa*( 10.0*( vn(stc+1,m) - vn(stc,m  ) )          -  &
                                5.0*( vn(stc+2,m) - vn(stc-1,m) )          +  &
                                    ( vn(stc+3,m) - vn(stc-2,m) ) ) 
            
            vr(stc+1,m) =   ( 150.0*( vn(stc+1,m) + vn(stc,m  ) )          -  &
                               25.0*( vn(stc+2,m) + vn(stc-1,m) )          +  &
                                3.0*( vn(stc+3,m) + vn(stc-2,m) ) ) /256.0 +  &
                          aa*( 10.0*( vn(stc+1,m) - vn(stc,m  ) )          -  &
                                5.0*( vn(stc+2,m) - vn(stc-1,m) )          +  &
                                    ( vn(stc+3,m) - vn(stc-2,m) ) )      
        end do
		stc=stc+1

    end if
    !>

    !>
    !!�����߽�HDCS�ұ߽�������Ҳ��ֵ
    if (nfe == nbc_bound_scheme) then
        ed = ni - 1      !< ����nfe=1ʱ���ڵ��ֵ��ʽ�������ֹ����
        edc = ni + 1     !< ����nfe=1ʱ����ֵ�������ֹ����
        ma(ed+1)=0.0                           !< ���Խ�׷�ϵ����Aϵ�����ٽ��ұ߽�����ֵ��
        ma(ed+2)=0.0                           !< ���Խ�׷�ϵ����Aϵ�����ұ߽�����ֵ��
        ma(ed+3)=0.0
        mc(ed+1)=0.0                           !< ���Խ�׷�ϵ����Cϵ�����ٽ��ұ߽�����ֵ��
        mc(ed+2)=0.0                           !< ���Խ�׷�ϵ����Cϵ�����ұ߽�����ֵ��
        mc(ed+3)=0.0

        do m=nst,ned
            vr(ni+1,m) = ( 15.0*vn(ni,m) - 10.0*vn(ni-1,m) +  3.0*vn(ni-2,m) ) / 8.0                ! third-order
            vl(ni  ,m) = (  3.0*vn(ni,m) +   6.0*vn(ni-1,m) -      vn(ni-2,m) ) / 8.0   ! third-order         !< �����Ҳ��ֵ��߽�
            vl(ni+1,m) = vr(ni+1,m)
            vr(ni  ,m) = vl(ni  ,m)                                                                           !< �Ҳ��ֵ�ұ߽�                               
        end do                                                                                               
    else
        ed = ni + 3    !< ����nfe=0������ȫ�ڵ�߽��ʽ����-1ʱ���ڵ��ֵ��ʽ�������ֹ����
        edc = ni + 4   !< ����nfe=0������ȫ�ڵ�߽��ʽ����-1ʱ����ֵ�������ֹ����
        ma(ed+1)=0.0
        mc(ed+1)=0.0

        edc=edc-1
		do m=nst,ned
            vl(edc+1,m) =   ( 150.0*( vn(edc+1,m) + vn(edc,m  ) )          -  &
                               25.0*( vn(edc+2,m) + vn(edc-1,m) )          +  &
                                3.0*( vn(edc+3,m) + vn(edc-2,m) ) ) /256.0 -  &
                          aa*( 10.0*( vn(edc+1,m) - vn(edc,m  ) )          -  &
                                5.0*( vn(edc+2,m) - vn(edc-1,m) )          +  &
                                    ( vn(edc+3,m) - vn(edc-2,m) ) ) 
            
            vr(edc+1,m) =   ( 150.0*( vn(edc+1,m) + vn(edc,m  ) )          -  &
                               25.0*( vn(edc+2,m) + vn(edc-1,m) )          +  &
                                3.0*( vn(edc+3,m) + vn(edc-2,m) ) ) /256.0 +  &
                          aa*( 10.0*( vn(edc+1,m) - vn(edc,m  ) )          -  &
                                5.0*( vn(edc+2,m) - vn(edc-1,m) )          +  &
                                    ( vn(edc+3,m) - vn(edc-2,m) ) )    
        end do
		edc=edc+1
    end if
    !>

    !>
    !!�����߽�HDCS�����Ҳ��ڵ��ֵ
    do m=nst,ned
        do i=st-1,ed-1
            vl(i+1,m) =    (  15.0 *( vn(i+1,m) + vn(i,m  ) )          +  &
                                    ( vn(i+2,m) + vn(i-1,m) ) ) /20.0  -  &
                        bb*(  15.0 *( vn(i+1,m) - vn(i,m  ) )          +  &
                               3.0 *( vn(i+2,m) - vn(i-1,m) ) ) /40.0            !< �߽�HDCS�ڵ��ֵ����ࣩ @attention vlֻ���߽�HDCS��ֵ��ʽ���Ҷ�  
            
            vr(i+1,m) =    (  15.0 *( vn(i+1,m) + vn(i,m  ) )          +  &
                                    ( vn(i+2,m) + vn(i-1,m) ) ) /20.0  +  &
                        bb*(  15.0 *( vn(i+1,m) - vn(i,m  ) )          +  &
                               3.0 *( vn(i+2,m) - vn(i-1,m) ) ) /40.0            !< ���HDCS�ڵ��ֵ���Ҳࣩ @attention vrֻ���߽�HDCS��ֵ��ʽ���Ҷ�   	                        
        end do
    end do
    
    do i=stc,edc
        mb(i)=1.0                          !< ���Խ�׷�ϵ����Bϵ��
    enddo
        
  
    do i=st,ed
        ma(i)=3.0/10.0*(1.0+bb)            !< ���Խ�׷�ϵ����Aϵ�����ڵ�����ֵ��
    end do
            


    do i=st,ed                                                
        mc(i)=3.0/10.0*(1.0-bb)            !< ���Խ�׷�ϵ����Cϵ�����ڵ�����ֵ��
    end do                                               


    
    do m=nst,ned
        do i=stc,edc
            md(i) = vl(i,m)               !< ���Խ�׷�ϵ��Ҳ�Dϵ��������ֵ��
        end do
        call tri_pursue(ma,mb,mc,md,ms,stc,edc,1-ngn,ni+ngn) !< ���Խ�׷��
        do i=stc,edc
            vl(i,m) = ms(i)               !< �������Խ�׷�Ͻ⣨����ֵ��
        end do
    end do
    
       
    do i=st,ed
        ma(i)=3.0/10.0*(1.0-bb)           !< ���Խ�׷�ϵ����Aϵ�����ڵ��Ҳ��ֵ @attention �߽缰�ٽ��߽�ϵ�����䣩
    end do
     
    do i=st,ed
        mc(i)=3.0/10.0*(1.0+bb)           !< ���Խ�׷�ϵ����Cϵ�����ڵ��Ҳ��ֵ @attention �߽缰�ٽ��߽�ϵ�����䣩
    end do       

    do m=nst,ned        
        do i=stc,edc
            md(i) = vr(i,m)               !< ���Խ�׷�ϵ��Ҳ�Dϵ�����Ҳ��ֵ��
        end do
        call tri_pursue(ma,mb,mc,md,ms,stc,edc,1-ngn,ni+ngn) !< ���Խ�׷��
        do i=stc,edc
            vr(i,m) = ms(i)               !< �������Խ�׷�Ͻ⣨�Ҳ��ֵ��
        end do                      
    end do 

    !>
    !!��������Ҳ��ֵ�Ƿ����
    do i=stc,edc
        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vl(i,m) < chklim(m) ) then
                    vl(i,m) = vn(i,m)
                    ndec = ndec + 1
                end if

                if (vr(i,m) < chklim(m) ) then
                    vr(i,m) = vn(i+1,m)
                    ndec = ndec + 1
                end if
           end if
        end do
    end do
    !>
              
end subroutine hdcs5cicc

subroutine hdcs5ciN(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec,st,ed)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_bound_scheme
    implicit none
    integer(kind_int), intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  , intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(in)  :: sn(1-ngn:ni+ngn,1:3)
    integer(kind_int), intent(in)  :: nfs,nfe
    integer(kind_int), intent(in)  :: nge,st,ed
    real(kind_real)  , intent(out) :: vl(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(out) :: vr(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , external    :: sub_limit
    integer(kind_int), intent(in)  :: chkid(nst:ned)
    real(kind_real),   intent(in)  :: chklim(nst:ned)
    integer(kind_int), intent(out) :: ndec
    real(kind_real)   :: ma(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵ����Aϵ��
    real(kind_real)   :: mb(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵ����Bϵ��
    real(kind_real)   :: mc(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵ����Cϵ��
    real(kind_real)   :: md(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵ��Ҳ�Dϵ��
    real(kind_real)   :: ms(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵļ�����     
    real(kind_real)   :: aa,bb
    integer(kind_int) :: i,m

    if (ni < 4) then
        call muscl2pv(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)
        return
    end if

    ndec = 0
    
    aa   = 0.01    
       
    bb   = 0.3
    
    ma(st)=0.0
    mc(st)=0.0

    do m=nst,ned
        vl(st,m) =   ( 150.0*( vn(st+1,m) + vn(st,m  ) )          -  &
                        25.0*( vn(st+2,m) + vn(st-1,m) )          +  &
                         3.0*( vn(st+3,m) + vn(st-2,m) ) ) /256.0 -  &
                   aa*( 10.0*( vn(st+1,m) - vn(st,m  ) )          -  &
                         5.0*( vn(st+2,m) - vn(st-1,m) )          +  &
                             ( vn(st+3,m) - vn(st-2,m) ) ) 
        
        vr(st,m) =   ( 150.0*( vn(st+1,m) + vn(st,m  ) )          -  &
                        25.0*( vn(st+2,m) + vn(st-1,m) )          +  &
                         3.0*( vn(st+3,m) + vn(st-2,m) ) ) /256.0 +  &
                   aa*( 10.0*( vn(st+1,m) - vn(st,m  ) )          -  &
                         5.0*( vn(st+2,m) - vn(st-1,m) )          +  &
                             ( vn(st+3,m) - vn(st-2,m) ) )      
    end do

    ma(ed)=0.0
    mc(ed)=0.0

    do m=nst,ned
        vl(ed,m) =   ( 150.0*( vn(ed+1,m) + vn(ed,m  ) )          -  &
                        25.0*( vn(ed+2,m) + vn(ed-1,m) )          +  &
                         3.0*( vn(ed+3,m) + vn(ed-2,m) ) ) /256.0 -  &
                   aa*( 10.0*( vn(ed+1,m) - vn(ed,m  ) )          -  &
                         5.0*( vn(ed+2,m) - vn(ed-1,m) )          +  &
                             ( vn(ed+3,m) - vn(ed-2,m) ) ) 
        
        vr(ed,m) =   ( 150.0*( vn(ed+1,m) + vn(ed,m  ) )          -  &
                        25.0*( vn(ed+2,m) + vn(ed-1,m) )          +  &
                         3.0*( vn(ed+3,m) + vn(ed-2,m) ) ) /256.0 +  &
                   aa*( 10.0*( vn(ed+1,m) - vn(ed,m  ) )          -  &
                         5.0*( vn(ed+2,m) - vn(ed-1,m) )          +  &
                             ( vn(ed+3,m) - vn(ed-2,m) ) )    
    end do

    !>
    !!�����߽�HDCS�����Ҳ��ڵ��ֵ
    do m=nst,ned
        do i=st+1,ed-1
            vl(i,m) =    (  15.0 *( vn(i+1,m) + vn(i,m  ) )          +  &
                                  ( vn(i+2,m) + vn(i-1,m) ) ) /20.0  -  &
                      bb*(  15.0 *( vn(i+1,m) - vn(i,m  ) )          +  &
                             3.0 *( vn(i+2,m) - vn(i-1,m) ) ) /40.0            !< �߽�HDCS�ڵ��ֵ����ࣩ @attention vlֻ���߽�HDCS��ֵ��ʽ���Ҷ�  
            
            vr(i,m) =    (  15.0 *( vn(i+1,m) + vn(i,m  ) )          +  &
                                  ( vn(i+2,m) + vn(i-1,m) ) ) /20.0  +  &
                      bb*(  15.0 *( vn(i+1,m) - vn(i,m  ) )          +  &
                             3.0 *( vn(i+2,m) - vn(i-1,m) ) ) /40.0            !< ���HDCS�ڵ��ֵ���Ҳࣩ @attention vrֻ���߽�HDCS��ֵ��ʽ���Ҷ�   	                        
        end do
    end do
    
    do i=st,ed
        mb(i)=1.0                          !< ���Խ�׷�ϵ����Bϵ��
    enddo
        
  
    do i=st+1,ed-1
        ma(i)=3.0/10.0*(1.0+bb)            !< ���Խ�׷�ϵ����Aϵ�����ڵ�����ֵ��
    end do
            


    do i=st+1,ed-1
        mc(i)=3.0/10.0*(1.0-bb)            !< ���Խ�׷�ϵ����Cϵ�����ڵ�����ֵ��
    end do                                               


    
    do m=nst,ned
        do i=st,ed
            md(i) = vl(i,m)               !< ���Խ�׷�ϵ��Ҳ�Dϵ��������ֵ��
        end do
        call tri_pursue(ma,mb,mc,md,ms,st,ed,1-ngn,ni+ngn) !< ���Խ�׷��
        do i=st,ed
            vl(i,m) = ms(i)               !< �������Խ�׷�Ͻ⣨����ֵ��
        end do
    end do
    
       
    do i=st+1,ed-1
        ma(i)=3.0/10.0*(1.0-bb)           !< ���Խ�׷�ϵ����Aϵ�����ڵ��Ҳ��ֵ @attention �߽缰�ٽ��߽�ϵ�����䣩
    end do
     
    do i=st+1,ed-1
        mc(i)=3.0/10.0*(1.0+bb)           !< ���Խ�׷�ϵ����Cϵ�����ڵ��Ҳ��ֵ @attention �߽缰�ٽ��߽�ϵ�����䣩
    end do       

    do m=nst,ned        
        do i=st,ed
            md(i) = vr(i,m)               !< ���Խ�׷�ϵ��Ҳ�Dϵ�����Ҳ��ֵ��
        end do
        call tri_pursue(ma,mb,mc,md,ms,st,ed,1-ngn,ni+ngn) !< ���Խ�׷��
        do i=st,ed
            vr(i,m) = ms(i)               !< �������Խ�׷�Ͻ⣨�Ҳ��ֵ��
        end do                      
    end do 

    !>
    !!��������Ҳ��ֵ�Ƿ����
    do i=st,ed
        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vl(i,m) < chklim(m) ) then
                    vl(i,m) = vn(i,m)
                    ndec = ndec + 1
                end if

                if (vr(i,m) < chklim(m) ) then
                    vr(i,m) = vn(i+1,m)
                    ndec = ndec + 1
                end if
           end if
        end do
    end do
    !>
              
end subroutine hdcs5ciN

subroutine scsl3ci(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_bound_scheme
    implicit none
    integer(kind_int), intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  , intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(in)  :: sn(1-ngn:ni+ngn,1:3)
    integer(kind_int), intent(in)  :: nfs,nfe
    integer(kind_int), intent(in)  :: nge
    real(kind_real)  , intent(out) :: vl(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(out) :: vr(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , external    :: sub_limit
    integer(kind_int), intent(in)  :: chkid(nst:ned)
    real(kind_real),   intent(in)  :: chklim(nst:ned)
    integer(kind_int), intent(out) :: ndec
    real(kind_real)   :: d3(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: ve(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: aa,aa_scs
    integer(kind_int) :: i,m,n,st,ed,stc,edc
    
    if (ni < 6) then
        call muscl2pv(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)
        return
    end if

    ndec = 0
    
    call d3_via_T2(ni,nst,ned,ngn,vn,nfs,nfe,nge,d3)
    
    call vec_via_scsl4(ni,nst,ned,ngn,vn,nfs,nfe,nge,ve)
    
    aa     = 0.01
    aa_scs = 0.005
           
    if (nfs == nbc_bound_scheme) then
        
        st  = 3
        stc = 0

        do m=nst,ned
            vl(0,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                  ! third-order
            vl(1,m) = (  3.0*vn(1,m) +  6.0*vn(2,m) -      vn(3,m) ) / 8.0   ! third-order 
            !!vl(1,m) = (  35.0*vn(1,m) + 140.0*vn(2,m) -  70.0*vn(3,m) +  28.0*vn(4,m) -  5.0*vn(5,m) ) / 128.0   ! fifth-order
            !!vl(2,m) = (  35.0*vn(2,m) + 140.0*vn(3,m) -  70.0*vn(4,m) +  28.0*vn(5,m) -  5.0*vn(6,m) ) / 128.0   ! fifth-order
            vl(2,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 + &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )
            vr(0,m) = vl(0,m)
            vr(1,m) = vl(1,m)
            !!vr(2,m) = vl(2,m)
            vr(2,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 - &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )            
        end do
        
    else
        
        st  = -1
        stc = st
     
    end if
     
    if (nfe == nbc_bound_scheme) then
        
        ed  = ni - 3
        edc = ni

        do m=nst,ned
            vr(ni  ,m) = ( 15.0*vn(ni,m) - 10.0*vn(ni-1,m) +  3.0*vn(ni-2,m) ) / 8.0                ! third-order
            vl(ni-1,m) = (  3.0*vn(ni,m) +  6.0*vn(ni-1,m) -      vn(ni-2,m) ) / 8.0   ! third-order
            !!vl(ni-1,m) = (  35.0*vn(ni  ,  m) + 140.0*vn(ni-1,m) -  70.0*vn(ni-2,m) +  28.0*vn(ni-3,m) -  5.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            !!vl(ni-2,m) = (  35.0*vn(ni-1,  m) + 140.0*vn(ni-2,m) -  70.0*vn(ni-3,m) +  28.0*vn(ni-4,m) -  5.0*vn(ni-5,m) ) / 128.0   ! fifth-order
            vl(ni-2,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 + &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) )
            vl(ni  ,m) = vr(ni  ,m)
            vr(ni-1,m) = vl(ni-1,m)
            !!vr(ni-2,m) = vl(ni-2,m)
            vr(ni-2,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 - &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) )            
        end do
        
    else
        
        ed  = ni + 1
        edc = ed
            
    end if



    do m=nst,ned
        do i=st,ed
            vl(i,m) =   ve(i,m) +  aa_scs*d3(i,m)
            
            vr(i,m) =   ve(i,m) -  aa_scs*d3(i,m)		                        
        end do
    end do

    do i=stc,edc
        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vl(i,m) < chklim(m) ) then
                    vl(i,m) = vn(i,m)
                    ndec = ndec + 1
                end if

                if (vr(i,m) < chklim(m) ) then
                    vr(i,m) = vn(i+1,m)
                    ndec = ndec + 1
                end if
           end if
        end do
    end do
    
end subroutine scsl3ci

subroutine scsl3cicc(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_bound_scheme
    implicit none
    integer(kind_int), intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  , intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(in)  :: sn(1-ngn:ni+ngn,1:3)
    integer(kind_int), intent(in)  :: nfs,nfe
    integer(kind_int), intent(in)  :: nge
    real(kind_real)  , intent(out) :: vl(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(out) :: vr(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , external    :: sub_limit
    integer(kind_int), intent(in)  :: chkid(nst:ned)
    real(kind_real),   intent(in)  :: chklim(nst:ned)
    integer(kind_int), intent(out) :: ndec
    real(kind_real)   :: d3(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: ve(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: aa,aa_scs
    integer(kind_int) :: i,m,n,st,ed,stc,edc
    
    if (ni < 6) then
        call muscl2pvcc(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)
        return
    end if

    ndec = 0
    
    call d3_via_T2(ni,nst,ned,ngn,vn,nfs,nfe,nge,d3)
    
    call vec_via_scsl4(ni,nst,ned,ngn,vn,nfs,nfe,nge,ve)
    
    aa     = 0.01
    aa_scs = 0.005
           
    if (nfs == nbc_bound_scheme) then
        
        st  = 4
        stc = 1

        do m=nst,ned
            vl(1,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                  ! third-order
            vl(2,m) = (  3.0*vn(1,m) +  6.0*vn(2,m) -      vn(3,m) ) / 8.0   ! third-order 
            !!vl(2,m) = (  35.0*vn(1,m) + 140.0*vn(2,m) -  70.0*vn(3,m) +  28.0*vn(4,m) -  5.0*vn(5,m) ) / 128.0   ! fifth-order
            !!vl(3,m) = (  35.0*vn(2,m) + 140.0*vn(3,m) -  70.0*vn(4,m) +  28.0*vn(5,m) -  5.0*vn(6,m) ) / 128.0   ! fifth-order
            vl(3,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 + &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )
            vr(1,m) = vl(1,m)
            vr(2,m) = vl(2,m)
            !!vr(2,m) = vl(2,m)
            vr(3,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 - &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )            
        end do
        
    else
        
        st  = 0
        stc = st
     
    end if
     
    if (nfe == nbc_bound_scheme) then
        
        ed  = ni - 2
        edc = ni + 1

        do m=nst,ned
            vr(ni+1,m) = ( 15.0*vn(ni,m) - 10.0*vn(ni-1,m) +  3.0*vn(ni-2,m) ) / 8.0                ! third-order
            vl(ni  ,m) = (  3.0*vn(ni,m) +  6.0*vn(ni-1,m) -      vn(ni-2,m) ) / 8.0   ! third-order
            !!vl(ni  ,m) = (  35.0*vn(ni  ,  m) + 140.0*vn(ni-1,m) -  70.0*vn(ni-2,m) +  28.0*vn(ni-3,m) -  5.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            !!vl(ni-1,m) = (  35.0*vn(ni-1,  m) + 140.0*vn(ni-2,m) -  70.0*vn(ni-3,m) +  28.0*vn(ni-4,m) -  5.0*vn(ni-5,m) ) / 128.0   ! fifth-order
            vl(ni-1,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 + &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) )
            vl(ni+1,m) = vr(ni+1,m)
            vr(ni  ,m) = vl(ni  ,m)
            !!vr(ni-2,m) = vl(ni-2,m)
            vr(ni-1,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 - &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) )            
        end do
        
    else
        
        ed  = ni + 2
        edc = ed
            
    end if



    do m=nst,ned
        do i=st-1,ed-1
            vl(i+1,m) =   ve(i,m) +  aa_scs*d3(i,m)
            
            vr(i+1,m) =   ve(i,m) -  aa_scs*d3(i,m)		                        
        end do
    end do

    do i=stc,edc
        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vl(i,m) < chklim(m) ) then
                    vl(i,m) = vn(i,m)
                    ndec = ndec + 1
                end if

                if (vr(i,m) < chklim(m) ) then
                    vr(i,m) = vn(i+1,m)
                    ndec = ndec + 1
                end if
           end if
        end do
    end do
    
end subroutine scsl3cicc

subroutine d3_via_T2(ni,nst,ned,ngn,vn,nfs,nfe,nge,d3)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_bound_scheme
    implicit none
    integer(kind_int),intent(in)  :: ni,nst,ned,ngn
    integer(kind_int),intent(in)  :: nfs,nfe
    integer(kind_int),intent(in)  :: nge    
    real(kind_real)  ,intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  ,intent(out) :: d3(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: ma(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵ����Aϵ��
    real(kind_real)   :: mb(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵ����Bϵ��
    real(kind_real)   :: mc(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵ����Cϵ��
    real(kind_real)   :: md(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵ��Ҳ�Dϵ��
    real(kind_real)   :: ms(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵļ�����
    real(kind_real)   :: x1                                        !< �ع���ʽϵ��
    integer(kind_int) :: i,m,st,ed,stc,edc
    
    x1=0.46
    
    !>
    !!������׾��Ȱ�����׿ռ䵼��
    if (nfs == nbc_bound_scheme) then

        st = 3       !< ����nfs=1ʱ���ڵ��ʽ�������ʼ����
        stc = 2      !< ����nfs=1ʱ���������ʼ����

        ma(st-1)=0.0
        mc(st-1)=0.0

        do m=nst,ned
            d3(2,m) = -3.0*(vn(3,m)-vn(2,m))+(vn(4,m)-vn(1,m))      !< �����ع���߽�
        end do 
        
    else
        
        st  = -1   !< ����nfs=0������ȫ�ڵ�߽��ʽ��ʱ���ڵ��ع���ʽ�������ʼ����
        stc = -2   !< ����nfs=0������ȫ�ڵ�߽��ʽ��ʱ���ع��������ʼ����
        
        ma(st-1)=0.0                           !< ���Խ�׷�ϵ����Aϵ������߽磩
        mc(st-1)=0.0                           !< ���Խ�׷�ϵ����Cϵ������߽磩
        
        do m=nst,ned
            d3(-2,m) = -3.0*(vn(-1,m)-vn(-2,m))+(vn(0,m)-vn(-3,m))   !< @attention �ع���δʹ�ã�ֻ��Ϊ���������Խ�׷��
        end do        
 
    end if
    !>

    if (nfe == nbc_bound_scheme) then
        
        ed = ni - 3      !< ����nfe=1ʱ���ڵ��ֵ��ʽ�������ֹ����
        edc = ni - 2     !< ����nfe=1ʱ����ֵ�������ֹ����

        ma(ed+1)=0.0                           !< ���Խ�׷�ϵ����Aϵ�����ұ߽磩
        mc(ed+1)=0.0                           !< ���Խ�׷�ϵ����Cϵ�����ұ߽磩

        do m=nst,ned
            d3(ni-2,m) = -3.0*(vn(ni-1,m)-vn(ni-2,m))+(vn(ni,m)-vn(ni-3,m))         !<  �����ع��ұ߽�                             
        end do
        
    else
        
        ed  = ni+1      !< ����nfe=0������ȫ�ڵ�߽��ʽ��ʱ���ڵ��ع���ʽ�������ֹ����
        edc = ni+2       !< ����nfe=0������ȫ�ڵ�߽��ʽ��ʱ���ع��������ֹ����
        
        ma(ed+1)=0.0
        mc(ed+1)=0.0

        do m=nst,ned
            d3(ni+2,m) = -3.0*(vn(ni+3,m)-vn(ni+2,m))+(vn(ni+4,m)-vn(ni+1,m))   !< @attention ��ֵ��δʹ�ã�ֻ��Ϊ���������Խ�׷��
        end do        
   
    end if
    !>

    !>
    !!������׾��Ȱ�����׿ռ䵼���ڵ�
    do m=nst,ned
        do i=st,ed
            d3(i,m) =    -3.0*(1.0+2.0*x1)*( vn(i+1,m) - vn(i,m  ) )          +  &
                              (1.0+2.0*x1)*( vn(i+2,m) - vn(i-1,m) )     
        end do
    end do
    
    do i=stc,edc
        mb(i)=1.0                          !< ���Խ�׷�ϵ�Bϵ��
    enddo
  
    do i=st,ed
        ma(i)=x1                           !< ���Խ�׷�ϵ�Aϵ�����ڵ㣩
    end do

    do i=st,ed                                                
        mc(i)=x1                           !< ���Խ�׷�ϵ�Cϵ�����ڵ㣩
    end do                                               


    
    do m=nst,ned
        do i=stc,edc
            md(i) = d3(i,m)               !< ���Խ�׷�ϵ��Ҳ�Dϵ��
        end do
        call tri_pursue(ma,mb,mc,md,ms,stc,edc,1-ngn,ni+ngn) !< ���Խ�׷��
        do i=stc,edc
            d3(i,m) = ms(i)               !< �������Խ�׷�Ͻ�
        end do
    end do    
    
end subroutine d3_via_T2

subroutine d3_via_T2N(ni,nst,ned,ngn,vn,nfs,nfe,nge,d3,st,ed)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_bound_scheme
    implicit none
    integer(kind_int),intent(in)  :: ni,nst,ned,ngn
    integer(kind_int),intent(in)  :: nfs,nfe
    integer(kind_int),intent(in)  :: nge,st,ed
    real(kind_real)  ,intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  ,intent(out) :: d3(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: ma(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵ����Aϵ��
    real(kind_real)   :: mb(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵ����Bϵ��
    real(kind_real)   :: mc(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵ����Cϵ��
    real(kind_real)   :: md(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵ��Ҳ�Dϵ��
    real(kind_real)   :: ms(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵļ�����
    real(kind_real)   :: x1                                        !< �ع���ʽϵ��
    integer(kind_int) :: i,m
    
    x1=0.46
    
    ma(st)=0.0
    mc(st)=0.0
        
    do m=nst,ned
        d3(st,m) =    -3.0*( vn(st+1,m) - vn(st,m  ) )          +  &
                           ( vn(st+2,m) - vn(st-1,m) ) 
    end do

    ma(ed)=0.0
    mc(ed)=0.0

    do m=nst,ned
        d3(ed,m) =    -3.0*( vn(ed+1,m) - vn(ed,m  ) )          +  &
                           ( vn(ed+2,m) - vn(ed-1,m) ) 
    end do

    !>
    !!������׾��Ȱ�����׿ռ䵼���ڵ�
    do m=nst,ned
        do i=st+1,ed-1
            d3(i,m) =    -3.0*(1.0+2.0*x1)*( vn(i+1,m) - vn(i,m  ) )          +  &
                              (1.0+2.0*x1)*( vn(i+2,m) - vn(i-1,m) )     
        end do
    end do
    
    do i=st,ed
        mb(i)=1.0                          !< ���Խ�׷�ϵ�Bϵ��
    enddo
  
    do i=st+1,ed-1
        ma(i)=x1                           !< ���Խ�׷�ϵ�Aϵ�����ڵ㣩
    end do

    do i=st+1,ed-1
        mc(i)=x1                           !< ���Խ�׷�ϵ�Cϵ�����ڵ㣩
    end do                                               


    
    do m=nst,ned
        do i=st,ed
            md(i) = d3(i,m)               !< ���Խ�׷�ϵ��Ҳ�Dϵ��
        end do
        call tri_pursue(ma,mb,mc,md,ms,st,ed,1-ngn,ni+ngn) !< ���Խ�׷��
        do i=st,ed
            d3(i,m) = ms(i)               !< �������Խ�׷�Ͻ�
        end do
    end do    
    
end subroutine d3_via_T2N

subroutine vec_via_scsl4(ni,nst,ned,ngn,vn,nfs,nfe,nge,ve)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_bound_scheme
    implicit none
    integer(kind_int),intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  ,intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)  :: nfs,nfe
    integer(kind_int),intent(in)  :: nge
    real(kind_real)  ,intent(out) :: ve(1-ngn:ni+ngn,nst:ned)
    real(kind_real)               :: ve2(1-ngn:ni+ngn,nst:ned)   
    integer(kind_int) :: i,m,st,ed,ierr
    
    if (ni < 6) then
        call ve_via_node2(ni,nst,ned,ngn,vn,nfs,nfe,nge,ve)
        return
    end if    
    
    call vec2_via_T2(ni,nst,ned,ngn,vn,nfs,nfe,nge,ve2)
    
    if (nfs == nbc_bound_scheme) then
        st = 3

        do m=nst,ned
            !!ve(2,m) = (      -vn(1,m) +   9.0*vn(2,m) +   9.0*vn(3,m) -       vn(4,m) ) / 16.0                   ! fourth-order
            !!ve(1,m) = (   5.0*vn(1,m) +  15.0*vn(2,m) -   5.0*vn(3,m) +       vn(4,m) ) / 16.0                   ! fourth-order
            !!ve(0,m) = (  35.0*vn(1,m) -  35.0*vn(2,m) +  21.0*vn(3,m) -   5.0*vn(4,m) ) / 16.0                   ! fourth-order
            ve(0,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                                       ! third-order
            ve(2,m) = (  -5.0*vn(1,m) +  60.0*vn(2,m) +  90.0*vn(3,m) -  20.0*vn(4,m) +  3.0*vn(5,m) ) / 128.0   ! fifth-order
            ve(1,m) = (  35.0*vn(1,m) + 140.0*vn(2,m) -  70.0*vn(3,m) +  28.0*vn(4,m) -  5.0*vn(5,m) ) / 128.0   ! fifth-order
            !!ve(0,m) = ( 315.0*vn(1,m) - 420.0*vn(2,m) + 378.0*vn(3,m) - 180.0*vn(4,m) + 35.0*vn(5,m) ) / 128.0   ! fifth-order
            !!ve(1,m) = ( 3.0*vn(1,m) + 6.0*vn(2,m) -     vn(3,m) ) / 8.0                                          ! third-order
            !!ve(2,m) = (    -vn(1,m) + 9.0*vn(2,m) + 9.0*vn(3,m) - vn(4,m) ) / 16.0                               ! fourth-order
        end do
        
    else
        
        st = -1
        
    end if
    
    if (nfe == nbc_bound_scheme) then
        
        ed = ni - 3

        do m=nst,ned
            !!ve(ni-2,m) = (      -vn(ni,m) +   9.0*vn(ni-1,m) +   9.0*vn(ni-2,m) -       vn(ni-3,m) ) / 16.0                      ! fourth-order
            !!ve(ni-1,m) = (   5.0*vn(ni,m) +  15.0*vn(ni-1,m) -   5.0*vn(ni-2,m) +       vn(ni-3,m) ) / 16.0                      ! fourth-order
            !!ve(ni  ,m) = (  35.0*vn(ni,m) -  35.0*vn(ni-1,m) +  21.0*vn(ni-2,m) -   5.0*vn(ni-3,m) ) / 16.0                      ! fourth-order
            ve(ni  ,m) = (  15.0*vn(ni,m) -  10.0*vn(ni-1,m) +   3.0*vn(ni-2,m) ) / 8.0                                          ! third-order
            ve(ni-2,m) = (  -5.0*vn(ni,m) +  60.0*vn(ni-1,m) +  90.0*vn(ni-2,m) -  20.0*vn(ni-3,m) +  3.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            ve(ni-1,m) = (  35.0*vn(ni,m) + 140.0*vn(ni-1,m) -  70.0*vn(ni-2,m) +  28.0*vn(ni-3,m) -  5.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            !!ve(ni  ,m) = ( 315.0*vn(ni,m) - 420.0*vn(ni-1,m) + 378.0*vn(ni-2,m) - 180.0*vn(ni-3,m) + 35.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            !!ve(ni-1,m) = (  3.0*vn(ni,m) +   6.0*vn(ni-1,m) -       vn(ni-2,m) ) / 8.0                                          ! third-order
            !!ve(ni-2,m) = (     -vn(ni,m) +   9.0*vn(ni-1,m) +   9.0*vn(ni-2,m) -       vn(ni-3,m) ) / 16.0                      ! fourth-order
        end do
        
    else
    
        ed = ni + 1 
      
    end if

    do m=nst,ned
        do i=st,ed
            ve(i,m) = 0.5*(vn(i+1,m) + vn(i  ,m)) - &
                      1.0/8.0*ve2(i,m)
        end do
    end do

end subroutine vec_via_scsl4

subroutine vec_via_scsl4cc(ni,nst,ned,ngn,vn,nfs,nfe,nge,ve)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_bound_scheme
    implicit none
    integer(kind_int),intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  ,intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)  :: nfs,nfe
    integer(kind_int),intent(in)  :: nge
    real(kind_real)  ,intent(out) :: ve(1-ngn:ni+ngn,nst:ned)
    real(kind_real)               :: ve2(1-ngn:ni+ngn,nst:ned)   
    integer(kind_int) :: i,m,st,ed,ierr
    
    if (ni < 6) then
        call ve_via_node2cc(ni,nst,ned,ngn,vn,nfs,nfe,nge,ve)
        return
    end if    
    
    call vec2_via_T2(ni,nst,ned,ngn,vn,nfs,nfe,nge,ve2)
    
    if (nfs == nbc_bound_scheme) then
        st = 4

        do m=nst,ned
            !!ve(3,m) = (      -vn(1,m) +   9.0*vn(2,m) +   9.0*vn(3,m) -       vn(4,m) ) / 16.0                   ! fourth-order
            !!ve(2,m) = (   5.0*vn(1,m) +  15.0*vn(2,m) -   5.0*vn(3,m) +       vn(4,m) ) / 16.0                   ! fourth-order
            !!ve(1,m) = (  35.0*vn(1,m) -  35.0*vn(2,m) +  21.0*vn(3,m) -   5.0*vn(4,m) ) / 16.0                   ! fourth-order
            ve(1,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                                       ! third-order
            ve(3,m) = (  -5.0*vn(1,m) +  60.0*vn(2,m) +  90.0*vn(3,m) -  20.0*vn(4,m) +  3.0*vn(5,m) ) / 128.0   ! fifth-order
            ve(2,m) = (  35.0*vn(1,m) + 140.0*vn(2,m) -  70.0*vn(3,m) +  28.0*vn(4,m) -  5.0*vn(5,m) ) / 128.0   ! fifth-order
            !!ve(1,m) = ( 315.0*vn(1,m) - 420.0*vn(2,m) + 378.0*vn(3,m) - 180.0*vn(4,m) + 35.0*vn(5,m) ) / 128.0   ! fifth-order
            !!ve(2,m) = ( 3.0*vn(1,m) + 6.0*vn(2,m) -     vn(3,m) ) / 8.0                                          ! third-order
            !!ve(3,m) = (    -vn(1,m) + 9.0*vn(2,m) + 9.0*vn(3,m) - vn(4,m) ) / 16.0                               ! fourth-order
        end do
        
    else
        
        st = 0
        
    end if
    
    if (nfe == nbc_bound_scheme) then
        
        ed = ni - 2

        do m=nst,ned
            !!ve(ni-1,m) = (      -vn(ni,m) +   9.0*vn(ni-1,m) +   9.0*vn(ni-2,m) -       vn(ni-3,m) ) / 16.0                      ! fourth-order
            !!ve(ni  ,m) = (   5.0*vn(ni,m) +  15.0*vn(ni-1,m) -   5.0*vn(ni-2,m) +       vn(ni-3,m) ) / 16.0                      ! fourth-order
            !!ve(ni+1,m) = (  35.0*vn(ni,m) -  35.0*vn(ni-1,m) +  21.0*vn(ni-2,m) -   5.0*vn(ni-3,m) ) / 16.0                      ! fourth-order
            ve(ni+1,m) = (  15.0*vn(ni,m) -  10.0*vn(ni-1,m) +   3.0*vn(ni-2,m) ) / 8.0                                          ! third-order
            ve(ni-1,m) = (  -5.0*vn(ni,m) +  60.0*vn(ni-1,m) +  90.0*vn(ni-2,m) -  20.0*vn(ni-3,m) +  3.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            ve(ni  ,m) = (  35.0*vn(ni,m) + 140.0*vn(ni-1,m) -  70.0*vn(ni-2,m) +  28.0*vn(ni-3,m) -  5.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            !!ve(ni+1,m) = ( 315.0*vn(ni,m) - 420.0*vn(ni-1,m) + 378.0*vn(ni-2,m) - 180.0*vn(ni-3,m) + 35.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            !!ve(ni  ,m) = (  3.0*vn(ni,m) +   6.0*vn(ni-1,m) -       vn(ni-2,m) ) / 8.0                                          ! third-order
            !!ve(ni-1,m) = (     -vn(ni,m) +   9.0*vn(ni-1,m) +   9.0*vn(ni-2,m) -       vn(ni-3,m) ) / 16.0                      ! fourth-order
        end do
        
    else
    
        ed = ni + 2 
      
    end if

    do m=nst,ned
        do i=st-1,ed-1
            ve(i+1,m) = 0.5*(vn(i+1,m) + vn(i  ,m)) - &
                        1.0/8.0*ve2(i,m)
        end do
    end do

end subroutine vec_via_scsl4cc

subroutine vec2_via_T2(ni,nst,ned,ngn,vn,nfs,nfe,nge,ve2)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_bound_scheme
    implicit none
    integer(kind_int),intent(in)  :: ni,nst,ned,ngn
    integer(kind_int),intent(in)  :: nfs,nfe
    integer(kind_int),intent(in)  :: nge    
    real(kind_real)  ,intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  ,intent(out) :: ve2(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: ma(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵ����Aϵ��
    real(kind_real)   :: mb(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵ����Bϵ��
    real(kind_real)   :: mc(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵ����Cϵ��
    real(kind_real)   :: md(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵ��Ҳ�Dϵ��
    real(kind_real)   :: ms(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵļ�����
    real(kind_real)   :: r1                                        !< �ع���ʽϵ��
    integer(kind_int) :: i,m,st,ed,stc,edc
    
    r1=0.44
    
    !>
    !!������׾��Ȱ����׿ռ䵼��                                                     
    if (nfs == nbc_bound_scheme) then

        st = 3       !< ����nfs=1ʱ���ڵ��ʽ�������ʼ����
        stc = 2      !< ����nfs=1ʱ���������ʼ����

        ma(st-1)=0.0
        mc(st-1)=0.0

        do m=nst,ned
            ve2(2,m) = -1.0/2.0*(vn(3,m)+vn(2,m))+1.0/2.0*(vn(4,m)+vn(1,m))      !< �����ع���߽�
        end do
        
    else
        
        st = -1     !< ����nfs=0������ȫ�ڵ�߽��ʽ��ʱ���ڵ��ع���ʽ�������ʼ����
        stc = -2    !< ����nfs=0������ȫ�ڵ�߽��ʽ��ʱ���ع��������ʼ����
        
        ma(st-1)=0.0                           !< ���Խ�׷�ϵ����Aϵ������߽磩
        mc(st-1)=0.0                           !< ���Խ�׷�ϵ����Cϵ������߽磩
        
        do m=nst,ned
            ve2(-2,m) = -1.0/2.0*(vn(-1,m)+vn(-2,m))+1.0/2.0*(vn(0,m)+vn(-3,m))   !< @attention �ع���δʹ�ã�ֻ��Ϊ���������Խ�׷��
        end do         
 
    end if
    !>

    !>        
    if (nfe == nbc_bound_scheme) then
        
        ed = ni - 3      !< ����nfe=1ʱ���ڵ��ֵ��ʽ�������ֹ����
        edc = ni - 2     !< ����nfe=1ʱ����ֵ�������ֹ����

        ma(ed+1)=0.0                           !< ���Խ�׷�ϵ����Aϵ�����ұ߽磩
        mc(ed+1)=0.0                           !< ���Խ�׷�ϵ����Cϵ�����ұ߽磩

        do m=nst,ned
            ve2(ni-2,m) = -1.0/2.0*(vn(ni-1,m)+vn(ni-2,m))+1.0/2.0*(vn(ni,m)+vn(ni-3,m))         !<  �����ع��ұ߽�                             
        end do
        
    else
        
        ed  = ni+1        !< ����nfe=0������ȫ�ڵ�߽��ʽ��ʱ���ڵ��ع���ʽ�������ֹ����
        edc = ni+2        !< ����nfe=0������ȫ�ڵ�߽��ʽ��ʱ���ع��������ֹ����
        
        ma(ed+1)=0.0
        mc(ed+1)=0.0

        do m=nst,ned
            ve2(ni+2,m) = -1.0/2.0*(vn(ni+3,m)+vn(ni+2,m))+1.0/2.0*(vn(ni+4,m)+vn(ni+1,m))   !< @attention ��ֵ��δʹ�ã�ֻ��Ϊ���������Խ�׷��
        end do         

    end if
    !>

    !>
    !!������׾��Ȱ����׿ռ䵼���ڵ�
    do m=nst,ned
        do i=st,ed
            ve2(i,m) =    (-0.5-r1)*( vn(i+1,m) + vn(i,m  ) )          +  &
                           (0.5+r1)*( vn(i+2,m) + vn(i-1,m) )     
        end do
    end do
    
    do i=stc,edc
        mb(i)=1.0                          !< ���Խ�׷�ϵ�Bϵ��
    enddo
  
    do i=st,ed
        ma(i)=r1                           !< ���Խ�׷�ϵ�Aϵ�����ڵ㣩
    end do

    do i=st,ed                                                
        mc(i)=r1                           !< ���Խ�׷�ϵ�Cϵ�����ڵ㣩
    end do                                               


    
    do m=nst,ned
        do i=stc,edc
            md(i) = ve2(i,m)               !< ���Խ�׷�ϵ��Ҳ�Dϵ��
        end do
        call tri_pursue(ma,mb,mc,md,ms,stc,edc,1-ngn,ni+ngn) !< ���Խ�׷��
        do i=stc,edc
            ve2(i,m) = ms(i)               !< �������Խ�׷�Ͻ�
        end do
    end do    
    
end subroutine vec2_via_T2

subroutine vec2_via_T2N(ni,nst,ned,ngn,vn,nfs,nfe,nge,ve2,st,ed)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_bound_scheme
    implicit none
    integer(kind_int),intent(in)  :: ni,nst,ned,ngn
    integer(kind_int),intent(in)  :: nfs,nfe
    integer(kind_int),intent(in)  :: nge,st,ed    
    real(kind_real)  ,intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  ,intent(out) :: ve2(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: ma(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵ����Aϵ��
    real(kind_real)   :: mb(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵ����Bϵ��
    real(kind_real)   :: mc(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵ����Cϵ��
    real(kind_real)   :: md(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵ��Ҳ�Dϵ��
    real(kind_real)   :: ms(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵļ�����
    real(kind_real)   :: r1                                        !< �ع���ʽϵ��
    integer(kind_int) :: i,m
    
    r1=0.44
    
    ma(st)=0.0
    mc(st)=0.0
        
    do m=nst,ned
        ve2(st,m) =    -0.5*( vn(st+1,m) + vn(st,m  ) )          +  &
                        0.5*( vn(st+2,m) + vn(st-1,m) )  
    end do

    ma(ed)=0.0
    mc(ed)=0.0

    do m=nst,ned
        ve2(ed,m) =    -0.5*( vn(ed+1,m) + vn(ed,m  ) )          +  &
                        0.5*( vn(ed+2,m) + vn(ed-1,m) )  
    end do

    do m=nst,ned
        do i=st+1,ed-1
            ve2(i,m) =    (-0.5-r1)*( vn(i+1,m) + vn(i,m  ) )          +  &
                           (0.5+r1)*( vn(i+2,m) + vn(i-1,m) )     
        end do
    end do
    
    do i=st,ed
        mb(i)=1.0                          !< ���Խ�׷�ϵ�Bϵ��
    enddo
  
    do i=st+1,ed-1
        ma(i)=r1                           !< ���Խ�׷�ϵ�Aϵ�����ڵ㣩
    end do

    do i=st+1,ed-1
        mc(i)=r1                           !< ���Խ�׷�ϵ�Cϵ�����ڵ㣩
    end do                                               


    
    do m=nst,ned
        do i=st,ed
            md(i) = ve2(i,m)               !< ���Խ�׷�ϵ��Ҳ�Dϵ��
        end do
        call tri_pursue(ma,mb,mc,md,ms,st,ed,1-ngn,ni+ngn) !< ���Խ�׷��
        do i=st,ed
            ve2(i,m) = ms(i)               !< �������Խ�׷�Ͻ�
        end do
    end do    
    
end subroutine vec2_via_T2N

subroutine scsh3pi2(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_bound_scheme
    implicit none
    integer(kind_int), intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  , intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(in)  :: sn(1-ngn:ni+ngn,1:3)
    integer(kind_int), intent(in)  :: nfs,nfe
    integer(kind_int), intent(in)  :: nge
    real(kind_real)  , intent(out) :: vl(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(out) :: vr(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , external    :: sub_limit
    integer(kind_int), intent(in)  :: chkid(nst:ned)
    real(kind_real),   intent(in)  :: chklim(nst:ned)
    integer(kind_int), intent(out) :: ndec  
    real(kind_real)   :: vln(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: vrn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: ve_lil(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: ve_lir(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: d2_li(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: d3_li(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: g(1:3,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: wl(1:3,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: wr(1:3,1-ngn:ni+ngn,nst:ned)    
    real(kind_real)   :: aa,aa_scs
    real(kind_real)   :: sigma(1-ngn:ni+ngn),imax,imin,indicator
    real(kind_real), parameter :: cl(3) = (/1.D0/4.D0,3.D0/4.D0,0.D0/4.D0/)
    real(kind_real), parameter :: cr(3) = (/0.D0/4.D0,3.D0/4.D0,1.D0/4.D0/)     
    real(kind_real), parameter :: eps        = 1.0e-30
    real(kind_real), parameter :: indicatorc = 0.1
    real(kind_real)            :: is,eis2(3),bs,tau3,bl(3),br(3)
    integer(kind_int)          :: i,m,n,j,st,ed,stc,edc
    integer(kind_int)          :: id,YN,ind(1:ni+2*ngn),ind_st,ind_ed,pre_ed
    integer(kind_int)          :: nsec,sec_st(1:ni+2*ngn),sec_ed(1:ni+2*ngn)
    integer(kind_int)          :: nsecn,sec_stn(1:ni+2*ngn),sec_edn(1:ni+2*ngn)
    
    if (ni < 8) then
        call muscl2pv(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)
        return
    end if
    
    ndec = 0
    
    aa     = 0.01
    aa_scs = 0.001
    
    wl(:,:,:) = 0.5
    wr(:,:,:) = 0.5
    vln(:,:)    = 0.0
    ve_lil(:,:) = 0.0
    vrn(:,:)    = 0.0
    ve_lir(:,:) = 0.0    
           
    if (nfs == nbc_bound_scheme) then
        
        st  = 4
        stc = 0

        do m=nst,ned
            vl(0,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                  ! third-order
            vl(1,m) = (  3.0*vn(1,m) +  6.0*vn(2,m) -      vn(3,m) ) / 8.0   ! third-order 
            vl(2,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 + &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )
            vr(0,m) = vl(0,m)
            vr(1,m) = vl(1,m)
            vr(2,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 - &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )
            vl(3,m) =   ( 150.0*( vn(4,m) + vn(3,m) )          -  &
                           25.0*( vn(5,m) + vn(2,m) )          +  &
                            3.0*( vn(6,m) + vn(1,m) ) ) /256.0 -  &
                      aa*( 10.0*( vn(4,m) - vn(3,m) )          -  &
                            5.0*( vn(5,m) - vn(2,m) )          +  &
                                ( vn(6,m) - vn(1,m) ) ) 
            
            vr(3,m) =   ( 150.0*( vn(4,m) + vn(3,m) )          -  &
                           25.0*( vn(5,m) + vn(2,m) )          +  &
                            3.0*( vn(6,m) + vn(1,m) ) ) /256.0 +  &
                      aa*( 10.0*( vn(4,m) - vn(3,m) )          -  &
                            5.0*( vn(5,m) - vn(2,m) )          +  &
                                ( vn(6,m) - vn(1,m) ) )          
        end do
        
    else
        
        st  = -2
        stc = st
     
    end if
     
    if (nfe == nbc_bound_scheme) then
        
        ed  = ni - 4
        edc = ni

        do m=nst,ned
            vr(ni  ,m) = ( 15.0*vn(ni,m) - 10.0*vn(ni-1,m) +  3.0*vn(ni-2,m) ) / 8.0                ! third-order
            vl(ni-1,m) = (  3.0*vn(ni,m) +  6.0*vn(ni-1,m) -      vn(ni-2,m) ) / 8.0   ! third-order
            vl(ni-2,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 + &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) )
            vl(ni  ,m) = vr(ni  ,m)
            vr(ni-1,m) = vl(ni-1,m)
            vr(ni-2,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 - &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) ) 
            vl(ni-3,m) =   ( 150.0*( vn(ni-2,m) + vn(ni-3,m) )          -  &
                              25.0*( vn(ni-1,m) + vn(ni-4,m) )          +  &
                               3.0*( vn(ni  ,m) + vn(ni-5,m) ) ) /256.0 -  &
                         aa*( 10.0*( vn(ni-2,m) - vn(ni-3,m) )          -  &
                               5.0*( vn(ni-1,m) - vn(ni-4,m) )          +  &
                                   ( vn(ni  ,m) - vn(ni-5,m) ) ) 
            
            vr(ni-3,m) =   ( 150.0*( vn(ni-2,m) + vn(ni-3,m) )          -  &
                              25.0*( vn(ni-1,m) + vn(ni-4,m) )          +  &
                               3.0*( vn(ni  ,m) + vn(ni-5,m) ) ) /256.0 +  &
                         aa*( 10.0*( vn(ni-2,m) - vn(ni-3,m) )          -  &
                               5.0*( vn(ni-1,m) - vn(ni-4,m) )          +  &
                                   ( vn(ni  ,m) - vn(ni-5,m) ) )              
        end do
        
    else
        
        ed  = ni + 2
        edc = ed
            
    end if
    
    do m=nst,ned
        do i=st-2,ed+2
            g(1,i,m) = vn(i  ,m) - vn(i-1,m)            
            g(2,i,m) = vn(i+1,m) - vn(i  ,m)
            g(3,i,m) = vn(i+2,m) - vn(i+1,m)            
        end do
    end do   
    
    do m=nst,ned
        do i=st-2,ed+2
            do n=1,3
                is = g(n,i,m)*g(n,i,m)
                eis2(n) = is**2
            end do

            tau3 = abs(eis2(2)-eis2(1))
            bl(1) = 1.0+tau3/(eis2(1)+eps)
            bl(2) = 1.0+tau3/(eis2(2)+eps)
            bs = bl(1) + bl(2)
            wl(2,i,m) = bl(2)/bs
            
            tau3 = abs(eis2(3)-eis2(2))
            br(2) = 1.0+tau3/(eis2(2)+eps)
            br(3) = 1.0+tau3/(eis2(3)+eps)
            bs = br(2) + br(3)
            wr(2,i,m) = br(2)/bs
        end do
    end do
    
    do i=st,ed
        imin=1.0
        imax=0.0
        do m=nst,ned
            imin         = min(wl(2,i-3,m),wl(2,i-2,m),wl(2,i-1,m),wl(2,i,m),wl(2,i+1,m),wl(2,i+2,m),wl(2,i+3,m), &
                               wr(2,i-3,m),wr(2,i-2,m),wr(2,i-1,m),wr(2,i,m),wr(2,i+1,m),wr(2,i+2,m),wr(2,i+3,m),imin)                                
            imax         = max(wl(2,i-3,m),wl(2,i-2,m),wl(2,i-1,m),wl(2,i,m),wl(2,i+1,m),wl(2,i+2,m),wl(2,i+3,m), &
                               wr(2,i-3,m),wr(2,i-2,m),wr(2,i-1,m),wr(2,i,m),wr(2,i+1,m),wr(2,i+2,m),wr(2,i+3,m),imax)
        end do
        indicator    = imin/imax
        if(indicator>indicatorc) then
            sigma(i)=1.0
        else
            sigma(i)=0.0
        endif        
    end do
    
    id=0
    do i=st,ed
        YN=0
        if(sigma(i)<0.5) then
            YN=YN+1
        end if
        if(YN>0) then
            id=id+1
            ind(id)=i           
        end if
    end do
       
    do i=1-ngn,ni+ngn
        sigma(i)=0.0
    end do
    if(ed-st<9) then
        ind_st=st
        ind_ed=ed
        nsec=0
        
        nsecn=1
        sec_stn(1)=st
        sec_edn(1)=ed
    else
        if(id==0) then
            nsec=1
            sec_st(1)=st
            sec_ed(1)=ed
            do n=sec_st(1),sec_ed(1)
                sigma(n)=1.0
            end do
            ind_st=st
            ind_ed=ed
            nsecn=0
        else
            ind_st=st
            ind_ed=ed
            nsec=0
            nsecn=0
            if(ind(1)-ind_st>9) then
                nsec=1
                sec_st(1)=ind_st
                sec_ed(1)=ind(1)-1
                pre_ed   =sec_ed(1)
                do n=sec_st(1),sec_ed(1)
                    sigma(n)=1.0
                end do
            else
                nsecn=1
                sec_stn(1)=ind_st
                sec_edn(1)=ind(1)
                pre_ed    =sec_edn(1)
            end if
            do i=2,id
                if(ind(i)-ind(i-1)>9) then
                    nsec=nsec+1
                    sec_st(nsec)=ind(i-1)+1
                    sec_ed(nsec)=ind(i)-1
                    do n=sec_st(nsec),sec_ed(nsec)
                        sigma(n)=1.0
                    end do
                    nsecn=nsecn+1
                    sec_stn(nsecn)=pre_ed+1
                    sec_edn(nsecn)=ind(i-1)
                    pre_ed        =sec_ed(nsec)
                else
                    if(ind(i)-ind(i-1)>1 .or. i==id) then
                        nsecn=nsecn+1
                        sec_stn(nsecn)=pre_ed+1
                        sec_edn(nsecn)=ind(i)
                        pre_ed        =sec_edn(nsecn)
                    end if
                end if
            end do
            if(ind_ed-ind(id)>9) then
                nsec=nsec+1
                sec_st(nsec)=pre_ed+1
                sec_ed(nsec)=ind_ed
                do n=sec_st(nsec),sec_ed(nsec)
                    sigma(n)=1.0
                end do
            else
                nsecn=nsecn+1
                sec_stn(nsecn)=pre_ed+1
                sec_edn(nsecn)=ind_ed
            end if
        end if
    end if
    
    do m=nst,ned
        do j=1,nsecn
            do i=sec_stn(j),sec_edn(j)
                g(1,i,m) = vn(i  ,m) - vn(i-1,m)
                g(2,i,m) = vn(i+1,m) - vn(i  ,m)
                g(3,i,m) = vn(i+2,m) - vn(i+1,m) 
            end do
        end do      
    end do            
            
    do m=nst,ned
        do j=1,nsecn
            do i=sec_stn(j),sec_edn(j)
                do n=1,3
                    is = g(n,i,m)*g(n,i,m)
                    eis2(n) = is**2
                end do
                
                tau3 = abs(eis2(2)-eis2(1))
                
                do n=1,3
                    bl(n) = cl(n)*(1.0+tau3/(eis2(n)+eps))
                end do
                
                bs = bl(1) + bl(2)
                do n=1,2
                    wl(n,i,m) = bl(n)/bs
                end do
          
                tau3 = abs(eis2(3)-eis2(2))           
                do n=2,3
                    br(n) = cr(n)*(1.0+tau3/(eis2(n)+eps))
                end do

                bs = br(2) + br(3)
                do n=2,3
                    wr(n,i,m) = br(n)/bs
                end do
            end do
        end do
    end do
    
!!    open(2,file='sigma.dat',form='formatted')
!!    write(2,*) 'variables="i" "sigmal" "sigmar"'
!!    do i=st,ed
!!
!!        write(2,200) i,sigmal(i,1),sigmar(i,1)
!!
!!    end do
!!
!!200 format(I4,1x,e12.5,1x,e12.5)
!!    close(2)
!!
!!    stop    
    
    do m=nst,ned
        do j=1,nsecn
            do i=sec_stn(j),sec_edn(j)
                vln(i,m) = vn(i  ,m) + ( wl(1,i,m)*( g(1,i,m)) + &
                                         wl(2,i,m)*( g(2,i,m)) )/2.0

                vrn(i,m) = vn(i+1,m) + ( wr(2,i,m)*(-g(2,i,m)) + &
                                         wr(3,i,m)*(-g(3,i,m)) )/2.0
            end do
        end do
    end do   
      
    do j=1,nsec
        call d3_via_T2N(ni,nst,ned,ngn,vn,nfs,nfe,nge,d3_li,sec_st(j),sec_ed(j))
        
        call vec2_via_T2N(ni,nst,ned,ngn,vn,nfs,nfe,nge,d2_li,sec_st(j),sec_ed(j))
    end do
    
    do m=nst,ned
        do j=1,nsec
            do i=sec_st(j),sec_ed(j)
                
                ve_lil(i,m) = 0.5*(vn(i+1,m) + vn(i  ,m)) - 1.0/8.0*d2_li(i,m) +  aa_scs*d3_li(i,m)
                
                ve_lir(i,m) = 0.5*(vn(i+1,m) + vn(i  ,m)) - 1.0/8.0*d2_li(i,m) -  aa_scs*d3_li(i,m)

            end do
        end do
    end do
        
    do m=nst,ned
        do i=st,ed
            
               vl(i,m) = sigma(i)*ve_lil(i,m) + vln(i,m)*(1.0-sigma(i))
                         
               vr(i,m) = sigma(i)*ve_lir(i,m) + vrn(i,m)*(1.0-sigma(i))            

        end do
    end do  

    do i=stc,edc
        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vl(i,m) < chklim(m) ) then
                    vl(i,m) = vn(i  ,m)
                    ndec = ndec + 1
                end if

                if (vr(i,m) < chklim(m) ) then
                    vr(i,m) = vn(i+1,m)
                    ndec = ndec + 1
                end if
           end if
        end do
    end do

end subroutine scsh3pi2

subroutine scsn2pi(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_bound_scheme
    implicit none
    integer(kind_int), intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  , intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(in)  :: sn(1-ngn:ni+ngn,1:3)
    integer(kind_int), intent(in)  :: nfs,nfe
    integer(kind_int), intent(in)  :: nge
    real(kind_real)  , intent(out) :: vl(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(out) :: vr(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , external    :: sub_limit
    integer(kind_int), intent(in)  :: chkid(nst:ned)
    real(kind_real),   intent(in)  :: chklim(nst:ned)
    integer(kind_int), intent(out) :: ndec
    real(kind_real)   :: g(1:3,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: wl(1:3,-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: wr(1:3,-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: aa
    real(kind_real), parameter :: cl(3) = (/1.D0/4.D0,3.D0/4.D0,0.D0/4.D0/)
    real(kind_real), parameter :: cr(3) = (/0.D0/4.D0,3.D0/4.D0,1.D0/4.D0/)    
    real(kind_real), parameter :: eps        = 1.0e-30
    real(kind_real)            :: is,eis2(3),bs,tau3,bl(3),br(3)    
    integer(kind_int)          :: i,m,n,st,ed,stc,edc
    
    if (ni < 6) then
        call muscl2pv(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)
        return
    end if

    ndec = 0 
    
    aa     = 0.01

    wl(:,:,:) = 0.5
    wr(:,:,:) = 0.5    
           
    if (nfs == nbc_bound_scheme) then
        
        st  = 3
        stc = 0

        do m=nst,ned
            vl(0,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                  ! third-order
            vl(1,m) = (  3.0*vn(1,m) +  6.0*vn(2,m) -      vn(3,m) ) / 8.0   ! third-order 

            vl(2,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 + &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )
            vr(0,m) = vl(0,m)
            vr(1,m) = vl(1,m)

            vr(2,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 - &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )         
        end do
        
    else
        
        st  = -1
        stc = st
     
    end if
     
    if (nfe == nbc_bound_scheme) then
        
        ed  = ni - 3
        edc = ni

        do m=nst,ned
            vr(ni  ,m) = ( 15.0*vn(ni,m) - 10.0*vn(ni-1,m) +  3.0*vn(ni-2,m) ) / 8.0                ! third-order
            vl(ni-1,m) = (  3.0*vn(ni,m) +  6.0*vn(ni-1,m) -      vn(ni-2,m) ) / 8.0   ! third-order
            vl(ni-2,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 + &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) )
            vl(ni  ,m) = vr(ni  ,m)
            vr(ni-1,m) = vl(ni-1,m)

            vr(ni-2,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 - &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) )           
        end do
        
    else
        
        ed  = ni + 1
        edc = ed
            
    end if

    do m=nst,ned
        do i=st-2,ed+2
            g(1,i,m) = vn(i  ,m) - vn(i-1,m) 
            g(2,i,m) = vn(i+1,m) - vn(i  ,m)
            g(3,i,m) = vn(i+2,m) - vn(i+1,m)
        end do         
    end do 
            
    do m=nst,ned
        do i=st-2,ed+2
            do n=1,3
                is = g(n,i,m)*g(n,i,m)
                eis2(n) = is**2
            end do
            
            tau3 = abs(eis2(2)-eis2(1))
            
            do n=1,3
                bl(n) = cl(n)*(1.0+tau3/(eis2(n)+eps))
            end do
            
            bs = bl(1) + bl(2)
            do n=1,2
                wl(n,i,m) = bl(n)/bs
            end do
          
            tau3 = abs(eis2(3)-eis2(2))           
            do n=2,3
                br(n) = cr(n)*(1.0+tau3/(eis2(n)+eps))
            end do

            bs = br(2) + br(3)
            do n=2,3
                wr(n,i,m) = br(n)/bs
            end do
        end do
    end do     
    
    do m=nst,ned
        do i=st,ed
            vl(i,m) = vn(i  ,m) + ( wl(1,i,m)*( g(1,i,m)) + &
                                    wl(2,i,m)*( g(2,i,m)) )/2.0
        end do

        do i=st,ed
            vr(i,m) = vn(i+1,m) + ( wr(2,i,m)*(-g(2,i,m)) + &
                                    wr(3,i,m)*(-g(3,i,m)) )/2.0
        end do
    end do

    do i=stc,edc
        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vl(i,m) < chklim(m) ) then
                    vl(i,m) = vn(i  ,m)
                    ndec = ndec + 1
                end if

                if (vr(i,m) < chklim(m) ) then
                    vr(i,m) = vn(i+1,m)
                    ndec = ndec + 1
                end if
           end if
        end do
    end do

end subroutine scsn2pi

subroutine scsn2picc(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_bound_scheme
    implicit none
    integer(kind_int), intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  , intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(in)  :: sn(1-ngn:ni+ngn,1:3)
    integer(kind_int), intent(in)  :: nfs,nfe
    integer(kind_int), intent(in)  :: nge
    real(kind_real)  , intent(out) :: vl(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(out) :: vr(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , external    :: sub_limit
    integer(kind_int), intent(in)  :: chkid(nst:ned)
    real(kind_real),   intent(in)  :: chklim(nst:ned)
    integer(kind_int), intent(out) :: ndec
    real(kind_real)   :: g(1:3,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: wl(1:3,-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: wr(1:3,-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: aa
    real(kind_real), parameter :: cl(3) = (/1.D0/4.D0,3.D0/4.D0,0.D0/4.D0/)
    real(kind_real), parameter :: cr(3) = (/0.D0/4.D0,3.D0/4.D0,1.D0/4.D0/)    
    real(kind_real), parameter :: eps        = 1.0e-30
    real(kind_real)            :: is,eis2(3),bs,tau3,bl(3),br(3)    
    integer(kind_int)          :: i,m,n,st,ed,stc,edc
    
    if (ni < 6) then
        call muscl2pvcc(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)
        return
    end if

    ndec = 0 
    
    aa     = 0.01

    wl(:,:,:) = 0.5
    wr(:,:,:) = 0.5    
           
    if (nfs == nbc_bound_scheme) then
        
        st  = 4
        stc = 1

        do m=nst,ned
            vl(1,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                  ! third-order
            vl(2,m) = (  3.0*vn(1,m) +  6.0*vn(2,m) -      vn(3,m) ) / 8.0   ! third-order 

            vl(3,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 + &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )
            vr(1,m) = vl(1,m)
            vr(2,m) = vl(2,m)

            vr(3,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 - &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )         
        end do
        
    else
        
        st  = 0
        stc = st
     
    end if
     
    if (nfe == nbc_bound_scheme) then
        
        ed  = ni - 2
        edc = ni + 1

        do m=nst,ned
            vr(ni+1,m) = ( 15.0*vn(ni,m) - 10.0*vn(ni-1,m) +  3.0*vn(ni-2,m) ) / 8.0                ! third-order
            vl(ni  ,m) = (  3.0*vn(ni,m) +  6.0*vn(ni-1,m) -      vn(ni-2,m) ) / 8.0   ! third-order
            vl(ni-1,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 + &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) )
            vl(ni+1,m) = vr(ni+1,m)
            vr(ni  ,m) = vl(ni  ,m)

            vr(ni-1,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 - &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) )           
        end do
        
    else
        
        ed  = ni + 2
        edc = ed
            
    end if

    do m=nst,ned
        do i=st-3,ed+1
            g(1,i,m) = vn(i  ,m) - vn(i-1,m) 
            g(2,i,m) = vn(i+1,m) - vn(i  ,m)
            g(3,i,m) = vn(i+2,m) - vn(i+1,m)
        end do         
    end do 
            
    do m=nst,ned
        do i=st-3,ed+1
            do n=1,3
                is = g(n,i,m)*g(n,i,m)
                eis2(n) = is**2
            end do
            
            tau3 = abs(eis2(2)-eis2(1))
            
            do n=1,3
                bl(n) = cl(n)*(1.0+tau3/(eis2(n)+eps))
            end do
            
            bs = bl(1) + bl(2)
            do n=1,2
                wl(n,i,m) = bl(n)/bs
            end do
          
            tau3 = abs(eis2(3)-eis2(2))           
            do n=2,3
                br(n) = cr(n)*(1.0+tau3/(eis2(n)+eps))
            end do

            bs = br(2) + br(3)
            do n=2,3
                wr(n,i,m) = br(n)/bs
            end do
        end do
    end do     
    
    do m=nst,ned
        do i=st-1,ed-1
            vl(i+1,m) = vn(i  ,m) + ( wl(1,i,m)*( g(1,i,m)) + &
                                      wl(2,i,m)*( g(2,i,m)) )/2.0
        end do

        do i=st-1,ed-1
            vr(i+1,m) = vn(i+1,m) + ( wr(2,i,m)*(-g(2,i,m)) + &
                                      wr(3,i,m)*(-g(3,i,m)) )/2.0
        end do
    end do

    do i=stc,edc
        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vl(i,m) < chklim(m) ) then
                    vl(i,m) = vn(i  ,m)
                    ndec = ndec + 1
                end if

                if (vr(i,m) < chklim(m) ) then
                    vr(i,m) = vn(i+1,m)
                    ndec = ndec + 1
                end if
           end if
        end do
    end do

end subroutine scsn2picc

subroutine scsh3ci2(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_bound_scheme
    implicit none
    integer(kind_int), intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  , intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(in)  :: sn(1-ngn:ni+ngn,1:3)
    integer(kind_int), intent(in)  :: nfs,nfe
    integer(kind_int), intent(in)  :: nge
    real(kind_real)  , intent(out) :: vl(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(out) :: vr(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , external    :: sub_limit
    integer(kind_int), intent(in)  :: chkid(nst:ned)
    real(kind_real),   intent(in)  :: chklim(nst:ned)
    integer(kind_int), intent(out) :: ndec  
    real(kind_real)   :: vln(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: vrn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: ve_lil(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: ve_lir(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: d2_li(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: d3_li(1-ngn:ni+ngn,nst:ned)  
    real(kind_real)   :: nt,nx,ny,nz,fp(nst:ned)
    real(kind_real)   :: vp(nst:ned),vc(nst:ned)
    real(kind_real)   :: q(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: q1(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: g(1:3,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: wl(1:3,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: wr(1:3,1-ngn:ni+ngn,nst:ned)    
    real(kind_real)   :: aa,aa_scs
    real(kind_real)   :: sigma(1-ngn:ni+ngn),imax,imin,indicator
    real(kind_real), parameter :: cl(3) = (/1.D0/4.D0,3.D0/4.D0,0.D0/4.D0/)
    real(kind_real), parameter :: cr(3) = (/0.D0/4.D0,3.D0/4.D0,1.D0/4.D0/)     
    real(kind_real), parameter :: eps        = 1.0e-30
    real(kind_real), parameter :: indicatorc = 0.1
    real(kind_real)            :: is,eis2(3),bs,tau3,bl(3),br(3)
    integer(kind_int)          :: i,m,n,j,st,ed,stc,edc
    integer(kind_int)          :: id,YN,ind(1:ni+2*ngn),ind_st,ind_ed,pre_ed
    integer(kind_int)          :: nsec,sec_st(1:ni+2*ngn),sec_ed(1:ni+2*ngn)
    integer(kind_int)          :: nsecn,sec_stn(1:ni+2*ngn),sec_edn(1:ni+2*ngn)
    
    if (ni < 8) then
        call muscl2pv(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)
        return
    end if
    
    ndec = 0
    
    aa     = 0.01
    aa_scs = 0.001
    
    wl(:,:,:) = 0.5
    wr(:,:,:) = 0.5
    vln(:,:)    = 0.0
    ve_lil(:,:) = 0.0
    vrn(:,:)    = 0.0
    ve_lir(:,:) = 0.0    
           
    if (nfs == nbc_bound_scheme) then
        
        st  = 4
        stc = 0

        do m=nst,ned
            vl(0,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                  ! third-order
            vl(1,m) = (  3.0*vn(1,m) +  6.0*vn(2,m) -      vn(3,m) ) / 8.0   ! third-order 
            vl(2,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 + &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )
            vr(0,m) = vl(0,m)
            vr(1,m) = vl(1,m)
            vr(2,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 - &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )
            vl(3,m) =   ( 150.0*( vn(4,m) + vn(3,m) )          -  &
                           25.0*( vn(5,m) + vn(2,m) )          +  &
                            3.0*( vn(6,m) + vn(1,m) ) ) /256.0 -  &
                      aa*( 10.0*( vn(4,m) - vn(3,m) )          -  &
                            5.0*( vn(5,m) - vn(2,m) )          +  &
                                ( vn(6,m) - vn(1,m) ) ) 
            
            vr(3,m) =   ( 150.0*( vn(4,m) + vn(3,m) )          -  &
                           25.0*( vn(5,m) + vn(2,m) )          +  &
                            3.0*( vn(6,m) + vn(1,m) ) ) /256.0 +  &
                      aa*( 10.0*( vn(4,m) - vn(3,m) )          -  &
                            5.0*( vn(5,m) - vn(2,m) )          +  &
                                ( vn(6,m) - vn(1,m) ) )          
        end do
        
    else
        
        st  = -2
        stc = st
     
    end if
     
    if (nfe == nbc_bound_scheme) then
        
        ed  = ni - 4
        edc = ni

        do m=nst,ned
            vr(ni  ,m) = ( 15.0*vn(ni,m) - 10.0*vn(ni-1,m) +  3.0*vn(ni-2,m) ) / 8.0                ! third-order
            vl(ni-1,m) = (  3.0*vn(ni,m) +  6.0*vn(ni-1,m) -      vn(ni-2,m) ) / 8.0   ! third-order
            vl(ni-2,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 + &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) )
            vl(ni  ,m) = vr(ni  ,m)
            vr(ni-1,m) = vl(ni-1,m)
            vr(ni-2,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 - &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) ) 
            vl(ni-3,m) =   ( 150.0*( vn(ni-2,m) + vn(ni-3,m) )          -  &
                              25.0*( vn(ni-1,m) + vn(ni-4,m) )          +  &
                               3.0*( vn(ni  ,m) + vn(ni-5,m) ) ) /256.0 -  &
                         aa*( 10.0*( vn(ni-2,m) - vn(ni-3,m) )          -  &
                               5.0*( vn(ni-1,m) - vn(ni-4,m) )          +  &
                                   ( vn(ni  ,m) - vn(ni-5,m) ) ) 
            
            vr(ni-3,m) =   ( 150.0*( vn(ni-2,m) + vn(ni-3,m) )          -  &
                              25.0*( vn(ni-1,m) + vn(ni-4,m) )          +  &
                               3.0*( vn(ni  ,m) + vn(ni-5,m) ) ) /256.0 +  &
                         aa*( 10.0*( vn(ni-2,m) - vn(ni-3,m) )          -  &
                               5.0*( vn(ni-1,m) - vn(ni-4,m) )          +  &
                                   ( vn(ni  ,m) - vn(ni-5,m) ) )              
        end do
        
    else
        
        ed  = ni + 2
        edc = ed
            
    end if
    
    do m=nst,ned
        do i=st-2,ed+2
            g(1,i,m) = vn(i  ,m) - vn(i-1,m)            
            g(2,i,m) = vn(i+1,m) - vn(i  ,m)
            g(3,i,m) = vn(i+2,m) - vn(i+1,m)            
        end do
    end do   
    
    do m=nst,ned
        do i=st-2,ed+2
            do n=1,3
                is = g(n,i,m)*g(n,i,m)
                eis2(n) = is**2
            end do

            tau3 = abs(eis2(2)-eis2(1))
            bl(1) = 1.0+tau3/(eis2(1)+eps)
            bl(2) = 1.0+tau3/(eis2(2)+eps)
            bs = bl(1) + bl(2)
            wl(2,i,m) = bl(2)/bs
            
            tau3 = abs(eis2(3)-eis2(2))
            br(2) = 1.0+tau3/(eis2(2)+eps)
            br(3) = 1.0+tau3/(eis2(3)+eps)
            bs = br(2) + br(3)
            wr(2,i,m) = br(2)/bs
        end do
    end do
    
    do i=st,ed
        imin=1.0
        imax=0.0
        do m=nst,ned
            imin         = min(wl(2,i-3,m),wl(2,i-2,m),wl(2,i-1,m),wl(2,i,m),wl(2,i+1,m),wl(2,i+2,m),wl(2,i+3,m), &
                               wr(2,i-3,m),wr(2,i-2,m),wr(2,i-1,m),wr(2,i,m),wr(2,i+1,m),wr(2,i+2,m),wr(2,i+3,m),imin)                                
            imax         = max(wl(2,i-3,m),wl(2,i-2,m),wl(2,i-1,m),wl(2,i,m),wl(2,i+1,m),wl(2,i+2,m),wl(2,i+3,m), &
                               wr(2,i-3,m),wr(2,i-2,m),wr(2,i-1,m),wr(2,i,m),wr(2,i+1,m),wr(2,i+2,m),wr(2,i+3,m),imax)
        end do
        indicator    = imin/imax
        if(indicator>indicatorc) then
            sigma(i)=1.0
        else
            sigma(i)=0.0
        endif        
    end do
    
    id=0
    do i=st,ed
        YN=0
        if(sigma(i)<0.5) then
            YN=YN+1
        end if
        if(YN>0) then
            id=id+1
            ind(id)=i           
        end if
    end do
       
    do i=1-ngn,ni+ngn
        sigma(i)=0.0
    end do
    if(ed-st<9) then
        ind_st=st
        ind_ed=ed
        nsec=0
        
        nsecn=1
        sec_stn(1)=st
        sec_edn(1)=ed
    else
        if(id==0) then
            nsec=1
            sec_st(1)=st
            sec_ed(1)=ed
            do n=sec_st(1),sec_ed(1)
                sigma(n)=1.0
            end do
            ind_st=st
            ind_ed=ed
            nsecn=0
        else
            ind_st=st
            ind_ed=ed
            nsec=0
            nsecn=0
            if(ind(1)-ind_st>9) then
                nsec=1
                sec_st(1)=ind_st
                sec_ed(1)=ind(1)-1
                pre_ed   =sec_ed(1)
                do n=sec_st(1),sec_ed(1)
                    sigma(n)=1.0
                end do
            else
                nsecn=1
                sec_stn(1)=ind_st
                sec_edn(1)=ind(1)
                pre_ed    =sec_edn(1)
            end if
            do i=2,id
                if(ind(i)-ind(i-1)>9) then
                    nsec=nsec+1
                    sec_st(nsec)=ind(i-1)+1
                    sec_ed(nsec)=ind(i)-1
                    do n=sec_st(nsec),sec_ed(nsec)
                        sigma(n)=1.0
                    end do
                    nsecn=nsecn+1
                    sec_stn(nsecn)=pre_ed+1
                    sec_edn(nsecn)=ind(i-1)
                    pre_ed        =sec_ed(nsec)
                else
                    if(ind(i)-ind(i-1)>1 .or. i==id) then
                        nsecn=nsecn+1
                        sec_stn(nsecn)=pre_ed+1
                        sec_edn(nsecn)=ind(i)
                        pre_ed        =sec_edn(nsecn)
                    end if
                end if
            end do
            if(ind_ed-ind(id)>9) then
                nsec=nsec+1
                sec_st(nsec)=pre_ed+1
                sec_ed(nsec)=ind_ed
                do n=sec_st(nsec),sec_ed(nsec)
                    sigma(n)=1.0
                end do
            else
                nsecn=nsecn+1
                sec_stn(nsecn)=pre_ed+1
                sec_edn(nsecn)=ind_ed
            end if
        end if
    end if
    
    do m=nst,ned
        do j=1,nsecn
            do i=sec_stn(j),sec_edn(j)
                g(1,i,m) = vn(i  ,m) - vn(i-1,m)
                g(2,i,m) = vn(i+1,m) - vn(i  ,m)
                g(3,i,m) = vn(i+2,m) - vn(i+1,m) 
            end do
        end do      
    end do
    
    do j=1,nsecn
        do i=sec_stn(j),sec_edn(j)
            nx = (sn(i,1) + sn(i+1,1))/2.0
            ny = (sn(i,2) + sn(i+1,2))/2.0
            nz = (sn(i,3) + sn(i+1,3))/2.0

            do m=nst,ned
                vp(m) = (vn(i,m) + vn(i+1,m))/2.0
                fp(m) = vn(i,m)
            end do
            call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
            do m=nst,ned
                q(i,m) = vc(m)
            end do
            
            do m=nst,ned
                fp(m) = vn(i+1,m)
            end do
            call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
            do m=nst,ned
                q1(i+1,m) = vc(m)
            end do        
             
            do n=1,3
                do m=nst,ned
                    fp(m) = g(n,i,m)
                end do
                call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
                do m=nst,ned
                    g(n,i,m) = vc(m)
                end do          
            end do
        end do
    end do   
            
            
    do m=nst,ned
        do j=1,nsecn
            do i=sec_stn(j),sec_edn(j)
                do n=1,3
                    is = g(n,i,m)*g(n,i,m)
                    eis2(n) = is**2
                end do
                
                tau3 = abs(eis2(2)-eis2(1))
                
                do n=1,3
                    bl(n) = cl(n)*(1.0+tau3/(eis2(n)+eps))
                end do
                
                bs = bl(1) + bl(2)
                do n=1,2
                    wl(n,i,m) = bl(n)/bs
                end do
          
                tau3 = abs(eis2(3)-eis2(2))           
                do n=2,3
                    br(n) = cr(n)*(1.0+tau3/(eis2(n)+eps))
                end do

                bs = br(2) + br(3)
                do n=2,3
                    wr(n,i,m) = br(n)/bs
                end do
            end do
        end do
    end do
    
!!    open(2,file='sigma.dat',form='formatted')
!!    write(2,*) 'variables="i" "sigmal" "sigmar"'
!!    do i=st,ed
!!
!!        write(2,200) i,sigmal(i,1),sigmar(i,1)
!!
!!    end do
!!
!!200 format(I4,1x,e12.5,1x,e12.5)
!!    close(2)
!!
!!    stop    
    
    do m=nst,ned
        do j=1,nsecn
            do i=sec_stn(j),sec_edn(j)
                vln(i,m) =  q(i  ,m) + ( wl(1,i,m)*( g(1,i,m)) + &
                                         wl(2,i,m)*( g(2,i,m)) )/2.0

                vrn(i,m) = q1(i+1,m) + ( wr(2,i,m)*(-g(2,i,m)) + &
                                         wr(3,i,m)*(-g(3,i,m)) )/2.0
            end do
        end do
    end do
    
    do j=1,nsecn
        do i=sec_stn(j),sec_edn(j)
            nx = (sn(i,1) + sn(i+1,1))/2.0
            ny = (sn(i,2) + sn(i+1,2))/2.0
            nz = (sn(i,3) + sn(i+1,3))/2.0

            do m=nst,ned
                vp(m) = (vn(i,m) + vn(i+1,m))/2.0
                vc(m) = vln(i,m)
            end do

            call RxQ_p(nst,ned,vp,nt,nx,ny,nz,vc,fp) 
            do m=nst,ned
                vln(i,m) = fp(m)
            end do
            
            do m=nst,ned
                vc(m) = vrn(i ,m)
            end do

            call RxQ_p(nst,ned,vp,nt,nx,ny,nz,vc,fp) 
            do m=nst,ned
                vrn(i,m) = fp(m)
            end do
        end do
    end do     
      
    do j=1,nsec
        call d3_via_T2N(ni,nst,ned,ngn,vn,nfs,nfe,nge,d3_li,sec_st(j),sec_ed(j))
        
        call vec2_via_T2N(ni,nst,ned,ngn,vn,nfs,nfe,nge,d2_li,sec_st(j),sec_ed(j))
    end do
    
    do m=nst,ned
        do j=1,nsec
            do i=sec_st(j),sec_ed(j)
                
                ve_lil(i,m) = 0.5*(vn(i+1,m) + vn(i  ,m)) - 1.0/8.0*d2_li(i,m) +  aa_scs*d3_li(i,m)
                
                ve_lir(i,m) = 0.5*(vn(i+1,m) + vn(i  ,m)) - 1.0/8.0*d2_li(i,m) -  aa_scs*d3_li(i,m)

            end do
        end do
    end do
        
    do m=nst,ned
        do i=st,ed
            
               vl(i,m) = sigma(i)*ve_lil(i,m) + vln(i,m)*(1.0-sigma(i))
                         
               vr(i,m) = sigma(i)*ve_lir(i,m) + vrn(i,m)*(1.0-sigma(i))            

        end do
    end do  

    do i=stc,edc
        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vl(i,m) < chklim(m) ) then
                    vl(i,m) = vn(i  ,m)
                    ndec = ndec + 1
                end if

                if (vr(i,m) < chklim(m) ) then
                    vr(i,m) = vn(i+1,m)
                    ndec = ndec + 1
                end if
           end if
        end do
    end do

end subroutine scsh3ci2

subroutine scsn2ci(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_bound_scheme
    implicit none
    integer(kind_int), intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  , intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(in)  :: sn(1-ngn:ni+ngn,1:3)
    integer(kind_int), intent(in)  :: nfs,nfe
    integer(kind_int), intent(in)  :: nge
    real(kind_real)  , intent(out) :: vl(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(out) :: vr(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , external    :: sub_limit
    integer(kind_int), intent(in)  :: chkid(nst:ned)
    real(kind_real),   intent(in)  :: chklim(nst:ned)
    integer(kind_int), intent(out) :: ndec
    real(kind_real)   :: nt,nx,ny,nz,fp(nst:ned)
    real(kind_real)   :: vp(nst:ned),vc(nst:ned)
    real(kind_real)   :: q(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: q1(1-ngn:ni+ngn,nst:ned)    
    real(kind_real)   :: g(1:3,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: wl(1:3,-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: wr(1:3,-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: aa
    real(kind_real), parameter :: cl(3) = (/1.D0/4.D0,3.D0/4.D0,0.D0/4.D0/)
    real(kind_real), parameter :: cr(3) = (/0.D0/4.D0,3.D0/4.D0,1.D0/4.D0/)    
    real(kind_real), parameter :: eps        = 1.0e-30
    real(kind_real)            :: is,eis2(3),bs,tau3,bl(3),br(3)    
    integer(kind_int)          :: i,m,n,st,ed,stc,edc
    
    if (ni < 6) then
        call muscl2pv(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)
        return
    end if

    ndec = 0  
    
    aa     = 0.01

    wl(:,:,:) = 0.5
    wr(:,:,:) = 0.5     
    
           
    if (nfs == nbc_bound_scheme) then
        
        st  = 3
        stc = 0

        do m=nst,ned
            vl(0,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                  ! third-order
            vl(1,m) = (  3.0*vn(1,m) +  6.0*vn(2,m) -      vn(3,m) ) / 8.0   ! third-order 

            vl(2,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 + &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )
            vr(0,m) = vl(0,m)
            vr(1,m) = vl(1,m)

            vr(2,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 - &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )         
        end do
        
    else
        
        st  = -1
        stc = st
     
    end if
     
    if (nfe == nbc_bound_scheme) then
        
        ed  = ni - 3
        edc = ni

        do m=nst,ned
            vr(ni  ,m) = ( 15.0*vn(ni,m) - 10.0*vn(ni-1,m) +  3.0*vn(ni-2,m) ) / 8.0                ! third-order
            vl(ni-1,m) = (  3.0*vn(ni,m) +  6.0*vn(ni-1,m) -      vn(ni-2,m) ) / 8.0   ! third-order
            vl(ni-2,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 + &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) )
            vl(ni  ,m) = vr(ni  ,m)
            vr(ni-1,m) = vl(ni-1,m)

            vr(ni-2,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 - &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) )           
        end do
        
    else
        
        ed  = ni + 1
        edc = ed
            
    end if

    do m=nst,ned
        do i=st-2,ed+2
            g(1,i,m) = vn(i  ,m) - vn(i-1,m)
            g(2,i,m) = vn(i+1,m) - vn(i  ,m)
            g(3,i,m) = vn(i+2,m) - vn(i+1,m)            
        end do         
    end do
    
    do i=st-2,ed+2
        nx = (sn(i,1) + sn(i+1,1))/2.0
        ny = (sn(i,2) + sn(i+1,2))/2.0
        nz = (sn(i,3) + sn(i+1,3))/2.0

        do m=nst,ned
            vp(m) = (vn(i,m) + vn(i+1,m))/2.0
            fp(m) = vn(i,m)
        end do
        call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc)  
        do m=nst,ned
            q(i,m) = vc(m)
        end do
        
        do m=nst,ned
            fp(m) = vn(i+1,m)
        end do
        call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
        do m=nst,ned
            q1(i+1,m) = vc(m)
        end do        
         
        do n=1,3
            do m=nst,ned
                fp(m) = g(n,i,m)
            end do
            call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
            do m=nst,ned
                g(n,i,m) = vc(m)
            end do
        end do
    end do 
            
    do m=nst,ned
        do i=st-2,ed+2
            do n=1,3
                is = g(n,i,m)*g(n,i,m)
                eis2(n) = is**2
            end do
            
            tau3 = abs(eis2(2)-eis2(1))
            
            do n=1,3
                bl(n) = cl(n)*(1.0+tau3/(eis2(n)+eps))
            end do
            
            bs = bl(1) + bl(2)
            do n=1,2
                wl(n,i,m) = bl(n)/bs
            end do
          
            tau3 = abs(eis2(3)-eis2(2))           
            do n=2,3
                br(n) = cr(n)*(1.0+tau3/(eis2(n)+eps))
            end do

            bs = br(2) + br(3)
            do n=2,3
                wr(n,i,m) = br(n)/bs
            end do
        end do
    end do         
    
    do m=nst,ned
        do i=st,ed
            vl(i,m) =  q(i  ,m) + ( wl(1,i,m)*( g(1,i,m)) + &
                                    wl(2,i,m)*( g(2,i,m)) )/2.0
        end do

        do i=st,ed
            vr(i,m) = q1(i+1,m) + ( wr(2,i,m)*(-g(2,i,m)) + &
                                    wr(3,i,m)*(-g(3,i,m)) )/2.0
        end do
    end do
    
    do i=st,ed
        nx = (sn(i,1) + sn(i+1,1))/2.0
        ny = (sn(i,2) + sn(i+1,2))/2.0
        nz = (sn(i,3) + sn(i+1,3))/2.0

        do m=nst,ned
            vp(m) = (vn(i,m) + vn(i+1,m))/2.0
            vc(m) = vl(i,m)
        end do

        call RxQ_p(nst,ned,vp,nt,nx,ny,nz,vc,fp) 
        do m=nst,ned
            vl(i,m) = fp(m)
        end do
        
        do m=nst,ned
            vc(m) = vr(i ,m)
        end do

        call RxQ_p(nst,ned,vp,nt,nx,ny,nz,vc,fp) 
        do m=nst,ned
            vr(i,m) = fp(m)
        end do        
    end do

    do i=stc,edc
        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vl(i,m) < chklim(m) ) then
                    vl(i,m) = vn(i  ,m)
                    ndec = ndec + 1
                end if

                if (vr(i,m) < chklim(m) ) then
                    vr(i,m) = vn(i+1,m)
                    ndec = ndec + 1
                end if
           end if
        end do
    end do
end subroutine scsn2ci

subroutine scsn2cicc(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_bound_scheme
    implicit none
    integer(kind_int), intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  , intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(in)  :: sn(1-ngn:ni+ngn,1:3)
    integer(kind_int), intent(in)  :: nfs,nfe
    integer(kind_int), intent(in)  :: nge
    real(kind_real)  , intent(out) :: vl(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(out) :: vr(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , external    :: sub_limit
    integer(kind_int), intent(in)  :: chkid(nst:ned)
    real(kind_real),   intent(in)  :: chklim(nst:ned)
    integer(kind_int), intent(out) :: ndec
    real(kind_real)   :: nt,nx,ny,nz,fp(nst:ned)
    real(kind_real)   :: vp(nst:ned),vc(nst:ned)
    real(kind_real)   :: q(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: q1(1-ngn:ni+ngn,nst:ned)    
    real(kind_real)   :: g(1:3,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: wl(1:3,-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: wr(1:3,-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: aa
    real(kind_real), parameter :: cl(3) = (/1.D0/4.D0,3.D0/4.D0,0.D0/4.D0/)
    real(kind_real), parameter :: cr(3) = (/0.D0/4.D0,3.D0/4.D0,1.D0/4.D0/)    
    real(kind_real), parameter :: eps        = 1.0e-30
    real(kind_real)            :: is,eis2(3),bs,tau3,bl(3),br(3)    
    integer(kind_int)          :: i,m,n,st,ed,stc,edc
    
    if (ni < 6) then
        call muscl2pvcc(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)
        return
    end if

    ndec = 0  
    
    aa     = 0.01

    wl(:,:,:) = 0.5
    wr(:,:,:) = 0.5     
    
           
    if (nfs == nbc_bound_scheme) then
        
        st  = 4
        stc = 1

        do m=nst,ned
            vl(1,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                  ! third-order
            vl(2,m) = (  3.0*vn(1,m) +  6.0*vn(2,m) -      vn(3,m) ) / 8.0   ! third-order 

            vl(3,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 + &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )
            vr(1,m) = vl(1,m)
            vr(2,m) = vl(2,m)

            vr(3,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 - &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )         
        end do
        
    else
        
        st  = 0
        stc = st
     
    end if
     
    if (nfe == nbc_bound_scheme) then
        
        ed  = ni - 2
        edc = ni + 1

        do m=nst,ned
            vr(ni+1,m) = ( 15.0*vn(ni,m) - 10.0*vn(ni-1,m) +  3.0*vn(ni-2,m) ) / 8.0                ! third-order
            vl(ni  ,m) = (  3.0*vn(ni,m) +  6.0*vn(ni-1,m) -      vn(ni-2,m) ) / 8.0   ! third-order
            vl(ni-1,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 + &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) )
            vl(ni+1,m) = vr(ni+1,m)
            vr(ni  ,m) = vl(ni  ,m)

            vr(ni-1,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 - &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) )           
        end do
        
    else
        
        ed  = ni + 2
        edc = ed
            
    end if

    do m=nst,ned
        do i=st-3,ed+1
            g(1,i,m) = vn(i  ,m) - vn(i-1,m)
            g(2,i,m) = vn(i+1,m) - vn(i  ,m)
            g(3,i,m) = vn(i+2,m) - vn(i+1,m)            
        end do         
    end do
    
    do i=st-3,ed+1
        nx = (sn(i,1) + sn(i+1,1))/2.0
        ny = (sn(i,2) + sn(i+1,2))/2.0
        nz = (sn(i,3) + sn(i+1,3))/2.0

        do m=nst,ned
            vp(m) = (vn(i,m) + vn(i+1,m))/2.0
            fp(m) = vn(i,m)
        end do
        call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc)  
        do m=nst,ned
            q(i,m) = vc(m)
        end do
        
        do m=nst,ned
            fp(m) = vn(i+1,m)
        end do
        call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
        do m=nst,ned
            q1(i+1,m) = vc(m)
        end do        
         
        do n=1,3
            do m=nst,ned
                fp(m) = g(n,i,m)
            end do
            call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
            do m=nst,ned
                g(n,i,m) = vc(m)
            end do
        end do
    end do 
            
    do m=nst,ned
        do i=st-3,ed+1
            do n=1,3
                is = g(n,i,m)*g(n,i,m)
                eis2(n) = is**2
            end do
            
            tau3 = abs(eis2(2)-eis2(1))
            
            do n=1,3
                bl(n) = cl(n)*(1.0+tau3/(eis2(n)+eps))
            end do
            
            bs = bl(1) + bl(2)
            do n=1,2
                wl(n,i,m) = bl(n)/bs
            end do
          
            tau3 = abs(eis2(3)-eis2(2))           
            do n=2,3
                br(n) = cr(n)*(1.0+tau3/(eis2(n)+eps))
            end do

            bs = br(2) + br(3)
            do n=2,3
                wr(n,i,m) = br(n)/bs
            end do
        end do
    end do         
    
    do m=nst,ned
        do i=st-1,ed-1
            vl(i+1,m) =  q(i  ,m) + ( wl(1,i,m)*( g(1,i,m)) + &
                                      wl(2,i,m)*( g(2,i,m)) )/2.0
        end do

        do i=st-1,ed-1
            vr(i+1,m) = q1(i+1,m) + ( wr(2,i,m)*(-g(2,i,m)) + &
                                      wr(3,i,m)*(-g(3,i,m)) )/2.0
        end do
    end do
    
    do i=st,ed
        nx = (sn(i,1) + sn(i+1,1))/2.0
        ny = (sn(i,2) + sn(i+1,2))/2.0
        nz = (sn(i,3) + sn(i+1,3))/2.0

        do m=nst,ned
            vp(m) = (vn(i,m) + vn(i+1,m))/2.0
            vc(m) = vl(i,m)
        end do

        call RxQ_p(nst,ned,vp,nt,nx,ny,nz,vc,fp) 
        do m=nst,ned
            vl(i,m) = fp(m)
        end do
        
        do m=nst,ned
            vc(m) = vr(i ,m)
        end do

        call RxQ_p(nst,ned,vp,nt,nx,ny,nz,vc,fp) 
        do m=nst,ned
            vr(i,m) = fp(m)
        end do        
    end do

    do i=stc,edc
        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vl(i,m) < chklim(m) ) then
                    vl(i,m) = vn(i  ,m)
                    ndec = ndec + 1
                end if

                if (vr(i,m) < chklim(m) ) then
                    vr(i,m) = vn(i+1,m)
                    ndec = ndec + 1
                end if
           end if
        end do
    end do
end subroutine scsn2cicc

subroutine scsh3pi(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_bound_scheme
    implicit none
    integer(kind_int), intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  , intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(in)  :: sn(1-ngn:ni+ngn,1:3)
    integer(kind_int), intent(in)  :: nfs,nfe
    integer(kind_int), intent(in)  :: nge
    real(kind_real)  , intent(out) :: vl(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(out) :: vr(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , external    :: sub_limit
    integer(kind_int), intent(in)  :: chkid(nst:ned)
    real(kind_real),   intent(in)  :: chklim(nst:ned)
    integer(kind_int), intent(out) :: ndec  
    real(kind_real)   :: vln(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: vrn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: ve_lil(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: ve_lir(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: d2_li(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: d3_li(1-ngn:ni+ngn,nst:ned)   
    real(kind_real)   :: s(1:4,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: g(1:4,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: wl(1:4,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: wr(1:4,1-ngn:ni+ngn,nst:ned)    
    real(kind_real)   :: aa,aa_scs
    real(kind_real)   :: sigma(1-ngn:ni+ngn),imax,imin,indicator
    real(kind_real), parameter :: cl(4) = (/1.D0/16.D0,10.D0/16.D0,5.D0/16.D0,0.D0/16.D0/)
    real(kind_real), parameter :: cr(4) = (/0.D0/16.D0,5.D0/16.D0,10.D0/16.D0,1.D0/16.D0/)      
    real(kind_real), parameter :: eps        = 1.0e-30
    real(kind_real), parameter :: indicatorc = 0.1
    real(kind_real)            :: is,eis2(4),bs,tau5,bl(4),br(4)
    integer(kind_int)          :: i,m,n,j,st,ed,stc,edc
    integer(kind_int)          :: id,YN,ind(1:ni+2*ngn),ind_st,ind_ed,pre_ed
    integer(kind_int)          :: nsec,sec_st(1:ni+2*ngn),sec_ed(1:ni+2*ngn)
    integer(kind_int)          :: nsecn,sec_stn(1:ni+2*ngn),sec_edn(1:ni+2*ngn)
    
    if (ni < 8) then
        call muscl2pv(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)
        return
    end if
    
    ndec = 0
    
    aa     = 0.01
    aa_scs = 0.001
    
    wl(:,:,:) = 0.3333
    wr(:,:,:) = 0.3333
    vln(:,:)    = 0.0
    ve_lil(:,:) = 0.0
    vrn(:,:)    = 0.0
    ve_lir(:,:) = 0.0    
           
    if (nfs == nbc_bound_scheme) then
        
        st  = 4
        stc = 0

        do m=nst,ned
            vl(0,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                  ! third-order
            vl(1,m) = (  3.0*vn(1,m) +  6.0*vn(2,m) -      vn(3,m) ) / 8.0   ! third-order 
            vl(2,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 + &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )
            vr(0,m) = vl(0,m)
            vr(1,m) = vl(1,m)
            vr(2,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 - &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )
            vl(3,m) =   ( 150.0*( vn(4,m) + vn(3,m) )          -  &
                           25.0*( vn(5,m) + vn(2,m) )          +  &
                            3.0*( vn(6,m) + vn(1,m) ) ) /256.0 -  &
                      aa*( 10.0*( vn(4,m) - vn(3,m) )          -  &
                            5.0*( vn(5,m) - vn(2,m) )          +  &
                                ( vn(6,m) - vn(1,m) ) ) 
            
            vr(3,m) =   ( 150.0*( vn(4,m) + vn(3,m) )          -  &
                           25.0*( vn(5,m) + vn(2,m) )          +  &
                            3.0*( vn(6,m) + vn(1,m) ) ) /256.0 +  &
                      aa*( 10.0*( vn(4,m) - vn(3,m) )          -  &
                            5.0*( vn(5,m) - vn(2,m) )          +  &
                                ( vn(6,m) - vn(1,m) ) )          
        end do
        
    else
        
        st  = -2
        stc = st
     
    end if
     
    if (nfe == nbc_bound_scheme) then
        
        ed  = ni - 4
        edc = ni

        do m=nst,ned
            vr(ni  ,m) = ( 15.0*vn(ni,m) - 10.0*vn(ni-1,m) +  3.0*vn(ni-2,m) ) / 8.0                ! third-order
            vl(ni-1,m) = (  3.0*vn(ni,m) +  6.0*vn(ni-1,m) -      vn(ni-2,m) ) / 8.0   ! third-order
            vl(ni-2,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 + &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) )
            vl(ni  ,m) = vr(ni  ,m)
            vr(ni-1,m) = vl(ni-1,m)
            vr(ni-2,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 - &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) ) 
            vl(ni-3,m) =   ( 150.0*( vn(ni-2,m) + vn(ni-3,m) )          -  &
                              25.0*( vn(ni-1,m) + vn(ni-4,m) )          +  &
                               3.0*( vn(ni  ,m) + vn(ni-5,m) ) ) /256.0 -  &
                         aa*( 10.0*( vn(ni-2,m) - vn(ni-3,m) )          -  &
                               5.0*( vn(ni-1,m) - vn(ni-4,m) )          +  &
                                   ( vn(ni  ,m) - vn(ni-5,m) ) ) 
            
            vr(ni-3,m) =   ( 150.0*( vn(ni-2,m) + vn(ni-3,m) )          -  &
                              25.0*( vn(ni-1,m) + vn(ni-4,m) )          +  &
                               3.0*( vn(ni  ,m) + vn(ni-5,m) ) ) /256.0 +  &
                         aa*( 10.0*( vn(ni-2,m) - vn(ni-3,m) )          -  &
                               5.0*( vn(ni-1,m) - vn(ni-4,m) )          +  &
                                   ( vn(ni  ,m) - vn(ni-5,m) ) )              
        end do
        
    else
        
        ed  = ni + 2
        edc = ed
            
    end if
    
    do m=nst,ned
        do i=st-2,ed+2
            g(1,i,m) = vn(i  ,m) - vn(i-1,m)            
            g(2,i,m) = vn(i+1,m) - vn(i  ,m)
            g(3,i,m) = vn(i+2,m) - vn(i+1,m)            
        end do
    end do   
    
    do m=nst,ned
        do i=st-2,ed+2
            do n=1,3
                is = g(n,i,m)*g(n,i,m)
                eis2(n) = is**2
            end do

            tau5 = abs(eis2(2)-eis2(1))
            bl(1) = 1.0+tau5/(eis2(1)+eps)
            bl(2) = 1.0+tau5/(eis2(2)+eps)
            bs = bl(1) + bl(2)
            wl(2,i,m) = bl(2)/bs
            
            tau5 = abs(eis2(3)-eis2(2))
            br(2) = 1.0+tau5/(eis2(2)+eps)
            br(3) = 1.0+tau5/(eis2(3)+eps)
            bs = br(2) + br(3)
            wr(2,i,m) = br(2)/bs
        end do
    end do
    
    do i=st,ed
        imin=1.0
        imax=0.0
        do m=nst,ned
            imin         = min(wl(2,i-3,m),wl(2,i-2,m),wl(2,i-1,m),wl(2,i,m),wl(2,i+1,m),wl(2,i+2,m),wl(2,i+3,m), &
                               wr(2,i-3,m),wr(2,i-2,m),wr(2,i-1,m),wr(2,i,m),wr(2,i+1,m),wr(2,i+2,m),wr(2,i+3,m),imin)                                
            imax         = max(wl(2,i-3,m),wl(2,i-2,m),wl(2,i-1,m),wl(2,i,m),wl(2,i+1,m),wl(2,i+2,m),wl(2,i+3,m), &
                               wr(2,i-3,m),wr(2,i-2,m),wr(2,i-1,m),wr(2,i,m),wr(2,i+1,m),wr(2,i+2,m),wr(2,i+3,m),imax)
        end do
        indicator    = imin/imax
        if(indicator>indicatorc) then
            sigma(i)=1.0
        else
            sigma(i)=0.0
        endif        
    end do
    
    id=0
    do i=st,ed
        YN=0
        if(sigma(i)<0.5) then
            YN=YN+1
        end if
        if(YN>0) then
            id=id+1
            ind(id)=i           
        end if
    end do
       
    do i=1-ngn,ni+ngn
        sigma(i)=0.0
    end do
    if(ed-st<9) then
        ind_st=st
        ind_ed=ed
        nsec=0
        
        nsecn=1
        sec_stn(1)=st
        sec_edn(1)=ed
    else
        if(id==0) then
            nsec=1
            sec_st(1)=st
            sec_ed(1)=ed
            do n=sec_st(1),sec_ed(1)
                sigma(n)=1.0
            end do
            ind_st=st
            ind_ed=ed
            nsecn=0
        else
            ind_st=st
            ind_ed=ed
            nsec=0
            nsecn=0
            if(ind(1)-ind_st>9) then
                nsec=1
                sec_st(1)=ind_st
                sec_ed(1)=ind(1)-1
                pre_ed   =sec_ed(1)
                do n=sec_st(1),sec_ed(1)
                    sigma(n)=1.0
                end do
            else
                nsecn=1
                sec_stn(1)=ind_st
                sec_edn(1)=ind(1)
                pre_ed    =sec_edn(1)
            end if
            do i=2,id
                if(ind(i)-ind(i-1)>9) then
                    nsec=nsec+1
                    sec_st(nsec)=ind(i-1)+1
                    sec_ed(nsec)=ind(i)-1
                    do n=sec_st(nsec),sec_ed(nsec)
                        sigma(n)=1.0
                    end do
                    nsecn=nsecn+1
                    sec_stn(nsecn)=pre_ed+1
                    sec_edn(nsecn)=ind(i-1)
                    pre_ed        =sec_ed(nsec)
                else
                    if(ind(i)-ind(i-1)>1 .or. i==id) then
                        nsecn=nsecn+1
                        sec_stn(nsecn)=pre_ed+1
                        sec_edn(nsecn)=ind(i)
                        pre_ed        =sec_edn(nsecn)
                    end if
                end if
            end do
            if(ind_ed-ind(id)>9) then
                nsec=nsec+1
                sec_st(nsec)=pre_ed+1
                sec_ed(nsec)=ind_ed
                do n=sec_st(nsec),sec_ed(nsec)
                    sigma(n)=1.0
                end do
            else
                nsecn=nsecn+1
                sec_stn(nsecn)=pre_ed+1
                sec_edn(nsecn)=ind_ed
            end if
        end if
    end if
    
    do m=nst,ned
        do j=1,nsecn
            do i=sec_stn(j),sec_edn(j)
                g(1,i,m) =      vn(i-2,m) - 3.0*vn(i-1,m) + 2.0*vn(i  ,m) 
                g(2,i,m) =     -vn(i  ,m) +     vn(i+1,m)
                g(3,i,m) =     -vn(i  ,m) +     vn(i+1,m)
                g(4,i,m) = -2.0*vn(i+1,m) + 3.0*vn(i+2,m) -     vn(i+3,m)

                s(1,i,m) = vn(i-2,m) - 2.0*vn(i-1,m) + vn(i,  m)
                s(2,i,m) = vn(i-1,m) - 2.0*vn(i  ,m) + vn(i+1,m) 
                s(3,i,m) = vn(i  ,m) - 2.0*vn(i+1,m) + vn(i+2,m)
                s(4,i,m) = vn(i+1,m) - 2.0*vn(i+2,m) + vn(i+3,m)
            end do
        end do      
    end do            
            
    do m=nst,ned
        do j=1,nsecn
            do i=sec_stn(j),sec_edn(j)
                do n=1,4
                    is = g(n,i,m)*g(n,i,m) + s(n,i,m)*s(n,i,m)
                    eis2(n) = is**2
                end do
                
                tau5 = abs(eis2(3)-eis2(1))
                
                do n=1,3
                    bl(n) = cl(n)*(1.0+tau5/(eis2(n)+eps))
                end do
                
                bs = bl(1) + bl(2) + bl(3)
                do n=1,3
                    wl(n,i,m) = bl(n)/bs
                end do
          
                tau5 = abs(eis2(4)-eis2(2))           
                do n=2,4
                    br(n) = cr(n)*(1.0+tau5/(eis2(n)+eps))
                end do

                bs = br(2) + br(3) + br(4)
                do n=2,4
                    wr(n,i,m) = br(n)/bs
                end do
            end do
        end do
    end do
    
!!    open(2,file='sigma.dat',form='formatted')
!!    write(2,*) 'variables="i" "sigmal" "sigmar"'
!!    do i=st,ed
!!
!!        write(2,200) i,sigmal(i,1),sigmar(i,1)
!!
!!    end do
!!
!!200 format(I4,1x,e12.5,1x,e12.5)
!!    close(2)
!!
!!    stop    
    
    do m=nst,ned
        do j=1,nsecn
            do i=sec_stn(j),sec_edn(j)
                vln(i,m) = vn(i  ,m) + ( wl(1,i,m)*( 4.0*g(1,i,m)-s(1,i,m)) + &
                                         wl(2,i,m)*( 4.0*g(2,i,m)-s(2,i,m)) + &
                                         wl(3,i,m)*( 4.0*g(3,i,m)-s(3,i,m)) )/8.0

                vrn(i,m) = vn(i+1,m) + ( wr(2,i,m)*(-4.0*g(2,i,m)-s(2,i,m)) + &
                                         wr(3,i,m)*(-4.0*g(3,i,m)-s(3,i,m)) + &
                                         wr(4,i,m)*(-4.0*g(4,i,m)-s(4,i,m)) )/8.0
            end do
        end do
    end do   
      
    do j=1,nsec
        call d3_via_T2N(ni,nst,ned,ngn,vn,nfs,nfe,nge,d3_li,sec_st(j),sec_ed(j))
        
        call vec2_via_T2N(ni,nst,ned,ngn,vn,nfs,nfe,nge,d2_li,sec_st(j),sec_ed(j))
    end do
    
    do m=nst,ned
        do j=1,nsec
            do i=sec_st(j),sec_ed(j)
                
                ve_lil(i,m) = 0.5*(vn(i+1,m) + vn(i  ,m)) - 1.0/8.0*d2_li(i,m) +  aa_scs*d3_li(i,m)
                
                ve_lir(i,m) = 0.5*(vn(i+1,m) + vn(i  ,m)) - 1.0/8.0*d2_li(i,m) -  aa_scs*d3_li(i,m)

            end do
        end do
    end do
        
    do m=nst,ned
        do i=st,ed
            
               vl(i,m) = sigma(i)*ve_lil(i,m) + vln(i,m)*(1.0-sigma(i))
                         
               vr(i,m) = sigma(i)*ve_lir(i,m) + vrn(i,m)*(1.0-sigma(i))            

        end do
    end do  

    do i=stc,edc
        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vl(i,m) < chklim(m) ) then
                    vl(i,m) = vn(i  ,m)
                    ndec = ndec + 1
                end if

                if (vr(i,m) < chklim(m) ) then
                    vr(i,m) = vn(i+1,m)
                    ndec = ndec + 1
                end if
           end if
        end do
    end do

end subroutine scsh3pi

subroutine scsh3picc(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_bound_scheme
    implicit none
    integer(kind_int), intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  , intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(in)  :: sn(1-ngn:ni+ngn,1:3)
    integer(kind_int), intent(in)  :: nfs,nfe
    integer(kind_int), intent(in)  :: nge
    real(kind_real)  , intent(out) :: vl(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(out) :: vr(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , external    :: sub_limit
    integer(kind_int), intent(in)  :: chkid(nst:ned)
    real(kind_real),   intent(in)  :: chklim(nst:ned)
    integer(kind_int), intent(out) :: ndec  
    real(kind_real)   :: vln(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: vrn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: ve_lil(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: ve_lir(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: d2_li(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: d3_li(1-ngn:ni+ngn,nst:ned)   
    real(kind_real)   :: s(1:4,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: g(1:4,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: wl(1:4,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: wr(1:4,1-ngn:ni+ngn,nst:ned)    
    real(kind_real)   :: aa,aa_scs
    real(kind_real)   :: sigma(1-ngn:ni+ngn),imax,imin,indicator
    real(kind_real), parameter :: cl(4) = (/1.D0/16.D0,10.D0/16.D0,5.D0/16.D0,0.D0/16.D0/)
    real(kind_real), parameter :: cr(4) = (/0.D0/16.D0,5.D0/16.D0,10.D0/16.D0,1.D0/16.D0/)      
    real(kind_real), parameter :: eps        = 1.0e-30
    real(kind_real), parameter :: indicatorc = 0.1
    real(kind_real)            :: is,eis2(4),bs,tau5,bl(4),br(4)
    integer(kind_int)          :: i,m,n,j,st,ed,stc,edc
    integer(kind_int)          :: id,YN,ind(1:ni+2*ngn),ind_st,ind_ed,pre_ed
    integer(kind_int)          :: nsec,sec_st(1:ni+2*ngn),sec_ed(1:ni+2*ngn)
    integer(kind_int)          :: nsecn,sec_stn(1:ni+2*ngn),sec_edn(1:ni+2*ngn)
    
    if (ni < 8) then
        call muscl2pvcc(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)
        return
    end if
    
    ndec = 0
    
    aa     = 0.01
    aa_scs = 0.001
    
    wl(:,:,:) = 0.3333
    wr(:,:,:) = 0.3333
    vln(:,:)    = 0.0
    ve_lil(:,:) = 0.0
    vrn(:,:)    = 0.0
    ve_lir(:,:) = 0.0    
           
    if (nfs == nbc_bound_scheme) then
        
        st  = 5
        stc = 1

        do m=nst,ned
            vl(1,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                  ! third-order
            vl(2,m) = (  3.0*vn(1,m) +  6.0*vn(2,m) -      vn(3,m) ) / 8.0   ! third-order 
            vl(3,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 + &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )
            vr(1,m) = vl(1,m)
            vr(2,m) = vl(2,m)
            vr(3,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 - &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )
            vl(4,m) =   ( 150.0*( vn(4,m) + vn(3,m) )          -  &
                           25.0*( vn(5,m) + vn(2,m) )          +  &
                            3.0*( vn(6,m) + vn(1,m) ) ) /256.0 -  &
                      aa*( 10.0*( vn(4,m) - vn(3,m) )          -  &
                            5.0*( vn(5,m) - vn(2,m) )          +  &
                                ( vn(6,m) - vn(1,m) ) ) 
            
            vr(4,m) =   ( 150.0*( vn(4,m) + vn(3,m) )          -  &
                           25.0*( vn(5,m) + vn(2,m) )          +  &
                            3.0*( vn(6,m) + vn(1,m) ) ) /256.0 +  &
                      aa*( 10.0*( vn(4,m) - vn(3,m) )          -  &
                            5.0*( vn(5,m) - vn(2,m) )          +  &
                                ( vn(6,m) - vn(1,m) ) )          
        end do
        
    else
        
        st  = -1
        stc = st
     
    end if
     
    if (nfe == nbc_bound_scheme) then
        
        ed  = ni - 3
        edc = ni + 1

        do m=nst,ned
            vr(ni+1,m) = ( 15.0*vn(ni,m) - 10.0*vn(ni-1,m) +  3.0*vn(ni-2,m) ) / 8.0                ! third-order
            vl(ni  ,m) = (  3.0*vn(ni,m) +  6.0*vn(ni-1,m) -      vn(ni-2,m) ) / 8.0   ! third-order
            vl(ni-1,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 + &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) )
            vl(ni+1,m) = vr(ni+1,m)
            vr(ni  ,m) = vl(ni  ,m)
            vr(ni-1,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 - &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) ) 
            vl(ni-2,m) =   ( 150.0*( vn(ni-2,m) + vn(ni-3,m) )          -  &
                              25.0*( vn(ni-1,m) + vn(ni-4,m) )          +  &
                               3.0*( vn(ni  ,m) + vn(ni-5,m) ) ) /256.0 -  &
                         aa*( 10.0*( vn(ni-2,m) - vn(ni-3,m) )          -  &
                               5.0*( vn(ni-1,m) - vn(ni-4,m) )          +  &
                                   ( vn(ni  ,m) - vn(ni-5,m) ) ) 
            
            vr(ni-2,m) =   ( 150.0*( vn(ni-2,m) + vn(ni-3,m) )          -  &
                              25.0*( vn(ni-1,m) + vn(ni-4,m) )          +  &
                               3.0*( vn(ni  ,m) + vn(ni-5,m) ) ) /256.0 +  &
                         aa*( 10.0*( vn(ni-2,m) - vn(ni-3,m) )          -  &
                               5.0*( vn(ni-1,m) - vn(ni-4,m) )          +  &
                                   ( vn(ni  ,m) - vn(ni-5,m) ) )              
        end do
        
    else
        
        ed  = ni + 3
        edc = ed
            
    end if
    
    do m=nst,ned
        do i=st-3,ed+1
            g(1,i,m) = vn(i  ,m) - vn(i-1,m)            
            g(2,i,m) = vn(i+1,m) - vn(i  ,m)
            g(3,i,m) = vn(i+2,m) - vn(i+1,m)            
        end do
    end do   
    
    do m=nst,ned
        do i=st-3,ed+1
            do n=1,3
                is = g(n,i,m)*g(n,i,m)
                eis2(n) = is**2
            end do

            tau5 = abs(eis2(2)-eis2(1))
            bl(1) = 1.0+tau5/(eis2(1)+eps)
            bl(2) = 1.0+tau5/(eis2(2)+eps)
            bs = bl(1) + bl(2)
            wl(2,i,m) = bl(2)/bs
            
            tau5 = abs(eis2(3)-eis2(2))
            br(2) = 1.0+tau5/(eis2(2)+eps)
            br(3) = 1.0+tau5/(eis2(3)+eps)
            bs = br(2) + br(3)
            wr(2,i,m) = br(2)/bs
        end do
    end do
    
    do i=st-1,ed-1
        imin=1.0
        imax=0.0
        do m=nst,ned
            imin         = min(wl(2,i-3,m),wl(2,i-2,m),wl(2,i-1,m),wl(2,i,m),wl(2,i+1,m),wl(2,i+2,m),wl(2,i+3,m), &
                               wr(2,i-3,m),wr(2,i-2,m),wr(2,i-1,m),wr(2,i,m),wr(2,i+1,m),wr(2,i+2,m),wr(2,i+3,m),imin)                                
            imax         = max(wl(2,i-3,m),wl(2,i-2,m),wl(2,i-1,m),wl(2,i,m),wl(2,i+1,m),wl(2,i+2,m),wl(2,i+3,m), &
                               wr(2,i-3,m),wr(2,i-2,m),wr(2,i-1,m),wr(2,i,m),wr(2,i+1,m),wr(2,i+2,m),wr(2,i+3,m),imax)
        end do
        indicator    = imin/imax
        if(indicator>indicatorc) then
            sigma(i)=1.0
        else
            sigma(i)=0.0
        endif        
    end do
    
    id=0
    do i=st-1,ed-1
        YN=0
        if(sigma(i)<0.5) then
            YN=YN+1
        end if
        if(YN>0) then
            id=id+1
            ind(id)=i           
        end if
    end do
       
    do i=1-ngn,ni+ngn
        sigma(i)=0.0
    end do
    if(ed-st<9) then
        ind_st=st-1
        ind_ed=ed-1
        nsec=0
        
        nsecn=1
        sec_stn(1)=st-1
        sec_edn(1)=ed-1
    else
        if(id==0) then
            nsec=1
            sec_st(1)=st-1
            sec_ed(1)=ed-1
            do n=sec_st(1),sec_ed(1)
                sigma(n)=1.0
            end do
            ind_st=st-1
            ind_ed=ed-1
            nsecn=0
        else
            ind_st=st-1
            ind_ed=ed-1
            nsec=0
            nsecn=0
            if(ind(1)-ind_st>9) then
                nsec=1
                sec_st(1)=ind_st
                sec_ed(1)=ind(1)-1
                pre_ed   =sec_ed(1)
                do n=sec_st(1),sec_ed(1)
                    sigma(n)=1.0
                end do
            else
                nsecn=1
                sec_stn(1)=ind_st
                sec_edn(1)=ind(1)
                pre_ed    =sec_edn(1)
            end if
            do i=2,id
                if(ind(i)-ind(i-1)>9) then
                    nsec=nsec+1
                    sec_st(nsec)=ind(i-1)+1
                    sec_ed(nsec)=ind(i)-1
                    do n=sec_st(nsec),sec_ed(nsec)
                        sigma(n)=1.0
                    end do
                    nsecn=nsecn+1
                    sec_stn(nsecn)=pre_ed+1
                    sec_edn(nsecn)=ind(i-1)
                    pre_ed        =sec_ed(nsec)
                else
                    if(ind(i)-ind(i-1)>1 .or. i==id) then
                        nsecn=nsecn+1
                        sec_stn(nsecn)=pre_ed+1
                        sec_edn(nsecn)=ind(i)
                        pre_ed        =sec_edn(nsecn)
                    end if
                end if
            end do
            if(ind_ed-ind(id)>9) then
                nsec=nsec+1
                sec_st(nsec)=pre_ed+1
                sec_ed(nsec)=ind_ed
                do n=sec_st(nsec),sec_ed(nsec)
                    sigma(n)=1.0
                end do
            else
                nsecn=nsecn+1
                sec_stn(nsecn)=pre_ed+1
                sec_edn(nsecn)=ind_ed
            end if
        end if
    end if
    
    do m=nst,ned
        do j=1,nsecn
            do i=sec_stn(j),sec_edn(j)
                g(1,i,m) =      vn(i-2,m) - 3.0*vn(i-1,m) + 2.0*vn(i  ,m) 
                g(2,i,m) =     -vn(i  ,m) +     vn(i+1,m)
                g(3,i,m) =     -vn(i  ,m) +     vn(i+1,m)
                g(4,i,m) = -2.0*vn(i+1,m) + 3.0*vn(i+2,m) -     vn(i+3,m)

                s(1,i,m) = vn(i-2,m) - 2.0*vn(i-1,m) + vn(i,  m)
                s(2,i,m) = vn(i-1,m) - 2.0*vn(i  ,m) + vn(i+1,m) 
                s(3,i,m) = vn(i  ,m) - 2.0*vn(i+1,m) + vn(i+2,m)
                s(4,i,m) = vn(i+1,m) - 2.0*vn(i+2,m) + vn(i+3,m)
            end do
        end do      
    end do            
            
    do m=nst,ned
        do j=1,nsecn
            do i=sec_stn(j),sec_edn(j)
                do n=1,4
                    is = g(n,i,m)*g(n,i,m) + s(n,i,m)*s(n,i,m)
                    eis2(n) = is**2
                end do
                
                tau5 = abs(eis2(3)-eis2(1))
                
                do n=1,3
                    bl(n) = cl(n)*(1.0+tau5/(eis2(n)+eps))
                end do
                
                bs = bl(1) + bl(2) + bl(3)
                do n=1,3
                    wl(n,i,m) = bl(n)/bs
                end do
          
                tau5 = abs(eis2(4)-eis2(2))           
                do n=2,4
                    br(n) = cr(n)*(1.0+tau5/(eis2(n)+eps))
                end do

                bs = br(2) + br(3) + br(4)
                do n=2,4
                    wr(n,i,m) = br(n)/bs
                end do
            end do
        end do
    end do
    
!!    open(2,file='sigma.dat',form='formatted')
!!    write(2,*) 'variables="i" "sigmal" "sigmar"'
!!    do i=st,ed
!!
!!        write(2,200) i,sigmal(i,1),sigmar(i,1)
!!
!!    end do
!!
!!200 format(I4,1x,e12.5,1x,e12.5)
!!    close(2)
!!
!!    stop    
    
    do m=nst,ned
        do j=1,nsecn
            do i=sec_stn(j),sec_edn(j)
                vln(i,m) = vn(i  ,m) + ( wl(1,i,m)*( 4.0*g(1,i,m)-s(1,i,m)) + &
                                         wl(2,i,m)*( 4.0*g(2,i,m)-s(2,i,m)) + &
                                         wl(3,i,m)*( 4.0*g(3,i,m)-s(3,i,m)) )/8.0

                vrn(i,m) = vn(i+1,m) + ( wr(2,i,m)*(-4.0*g(2,i,m)-s(2,i,m)) + &
                                         wr(3,i,m)*(-4.0*g(3,i,m)-s(3,i,m)) + &
                                         wr(4,i,m)*(-4.0*g(4,i,m)-s(4,i,m)) )/8.0
            end do
        end do
    end do   
      
    do j=1,nsec
        call d3_via_T2N(ni,nst,ned,ngn,vn,nfs,nfe,nge,d3_li,sec_st(j),sec_ed(j))
        
        call vec2_via_T2N(ni,nst,ned,ngn,vn,nfs,nfe,nge,d2_li,sec_st(j),sec_ed(j))
    end do
    
    do m=nst,ned
        do j=1,nsec
            do i=sec_st(j),sec_ed(j)
                
                ve_lil(i,m) = 0.5*(vn(i+1,m) + vn(i  ,m)) - 1.0/8.0*d2_li(i,m) +  aa_scs*d3_li(i,m)
                
                ve_lir(i,m) = 0.5*(vn(i+1,m) + vn(i  ,m)) - 1.0/8.0*d2_li(i,m) -  aa_scs*d3_li(i,m)

            end do
        end do
    end do
        
    do m=nst,ned
        do i=st-1,ed-1
            
               vl(i+1,m) = sigma(i)*ve_lil(i,m) + vln(i,m)*(1.0-sigma(i))
                         
               vr(i+1,m) = sigma(i)*ve_lir(i,m) + vrn(i,m)*(1.0-sigma(i))            

        end do
    end do  

    do i=stc,edc
        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vl(i,m) < chklim(m) ) then
                    vl(i,m) = vn(i  ,m)
                    ndec = ndec + 1
                end if

                if (vr(i,m) < chklim(m) ) then
                    vr(i,m) = vn(i+1,m)
                    ndec = ndec + 1
                end if
           end if
        end do
    end do

end subroutine scsh3picc

subroutine dcsh5pi(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_bound_scheme
    implicit none
    integer(kind_int), intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  , intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(in)  :: sn(1-ngn:ni+ngn,1:3)
    integer(kind_int), intent(in)  :: nfs,nfe
    integer(kind_int), intent(in)  :: nge
    real(kind_real)  , intent(out) :: vl(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(out) :: vr(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , external    :: sub_limit
    integer(kind_int), intent(in)  :: chkid(nst:ned)
    real(kind_real),   intent(in)  :: chklim(nst:ned)
    integer(kind_int), intent(out) :: ndec  
    real(kind_real)   :: vln(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: vrn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: ve_lil(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: ve_lir(1-ngn:ni+ngn,nst:ned)     
    real(kind_real)   :: s(1:4,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: g(1:4,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: wl(1:4,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: wr(1:4,1-ngn:ni+ngn,nst:ned)    
    real(kind_real)   :: aa,aa_scs
    real(kind_real)   :: sigma(1-ngn:ni+ngn),imax,imin,indicator
    real(kind_real), parameter :: cl(4) = (/1.D0/16.D0,10.D0/16.D0,5.D0/16.D0,0.D0/16.D0/)
    real(kind_real), parameter :: cr(4) = (/0.D0/16.D0,5.D0/16.D0,10.D0/16.D0,1.D0/16.D0/)  
    real(kind_real), parameter :: eps        = 1.0e-30
    real(kind_real), parameter :: indicatorc = 0.1
    real(kind_real)            :: is,eis2(4),bs,tau5,bl(4),br(4)
    integer(kind_int)          :: i,m,n,j,st,ed,stc,edc
    integer(kind_int)          :: id,YN,ind(1:ni+2*ngn),ind_st,ind_ed,pre_ed
    integer(kind_int)          :: nsec,sec_st(1:ni+2*ngn),sec_ed(1:ni+2*ngn)
    integer(kind_int)          :: nsecn,sec_stn(1:ni+2*ngn),sec_edn(1:ni+2*ngn)
    
    if (ni < 8) then
        call muscl2pv(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)
        return
    end if
    
    ndec = 0
    
    aa     = 0.01
    aa_scs = 0.001
    
    wl(:,:,:) = 0.3333
    wr(:,:,:) = 0.3333
    vln(:,:)    = 0.0
    ve_lil(:,:) = 0.0
    vrn(:,:)    = 0.0
    ve_lir(:,:) = 0.0    
           
    if (nfs == nbc_bound_scheme) then
        
        st  = 4
        stc = 0

        do m=nst,ned
            vl(0,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                  ! third-order
            vl(1,m) = (  3.0*vn(1,m) +  6.0*vn(2,m) -      vn(3,m) ) / 8.0   ! third-order 
            vl(2,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 + &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )
            vr(0,m) = vl(0,m)
            vr(1,m) = vl(1,m)
            vr(2,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 - &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )
            vl(3,m) =   ( 150.0*( vn(4,m) + vn(3,m) )          -  &
                           25.0*( vn(5,m) + vn(2,m) )          +  &
                            3.0*( vn(6,m) + vn(1,m) ) ) /256.0 -  &
                      aa*( 10.0*( vn(4,m) - vn(3,m) )          -  &
                            5.0*( vn(5,m) - vn(2,m) )          +  &
                                ( vn(6,m) - vn(1,m) ) ) 
            
            vr(3,m) =   ( 150.0*( vn(4,m) + vn(3,m) )          -  &
                           25.0*( vn(5,m) + vn(2,m) )          +  &
                            3.0*( vn(6,m) + vn(1,m) ) ) /256.0 +  &
                      aa*( 10.0*( vn(4,m) - vn(3,m) )          -  &
                            5.0*( vn(5,m) - vn(2,m) )          +  &
                                ( vn(6,m) - vn(1,m) ) )          
        end do
        
    else
        
        st  = -2
        stc = st
     
    end if
     
    if (nfe == nbc_bound_scheme) then
        
        ed  = ni - 4
        edc = ni

        do m=nst,ned
            vr(ni  ,m) = ( 15.0*vn(ni,m) - 10.0*vn(ni-1,m) +  3.0*vn(ni-2,m) ) / 8.0                ! third-order
            vl(ni-1,m) = (  3.0*vn(ni,m) +  6.0*vn(ni-1,m) -      vn(ni-2,m) ) / 8.0   ! third-order
            vl(ni-2,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 + &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) )
            vl(ni  ,m) = vr(ni  ,m)
            vr(ni-1,m) = vl(ni-1,m)
            vr(ni-2,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 - &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) ) 
            vl(ni-3,m) =   ( 150.0*( vn(ni-2,m) + vn(ni-3,m) )          -  &
                              25.0*( vn(ni-1,m) + vn(ni-4,m) )          +  &
                               3.0*( vn(ni  ,m) + vn(ni-5,m) ) ) /256.0 -  &
                         aa*( 10.0*( vn(ni-2,m) - vn(ni-3,m) )          -  &
                               5.0*( vn(ni-1,m) - vn(ni-4,m) )          +  &
                                   ( vn(ni  ,m) - vn(ni-5,m) ) ) 
            
            vr(ni-3,m) =   ( 150.0*( vn(ni-2,m) + vn(ni-3,m) )          -  &
                              25.0*( vn(ni-1,m) + vn(ni-4,m) )          +  &
                               3.0*( vn(ni  ,m) + vn(ni-5,m) ) ) /256.0 +  &
                         aa*( 10.0*( vn(ni-2,m) - vn(ni-3,m) )          -  &
                               5.0*( vn(ni-1,m) - vn(ni-4,m) )          +  &
                                   ( vn(ni  ,m) - vn(ni-5,m) ) )              
        end do
        
    else
        
        ed  = ni + 2
        edc = ed
            
    end if
    
    do m=nst,ned
        do i=st-2,ed+2
            g(1,i,m) = vn(i  ,m) - vn(i-1,m)            
            g(2,i,m) = vn(i+1,m) - vn(i  ,m)
            g(3,i,m) = vn(i+2,m) - vn(i+1,m)            
        end do
    end do   
    
    do m=nst,ned
        do i=st-2,ed+2
            do n=1,3
                is = g(n,i,m)*g(n,i,m)
                eis2(n) = is**2
            end do

            tau5 = abs(eis2(2)-eis2(1))
            bl(1) = 1.0+tau5/(eis2(1)+eps)
            bl(2) = 1.0+tau5/(eis2(2)+eps)
            bs = bl(1) + bl(2)
            wl(2,i,m) = bl(2)/bs
            
            tau5 = abs(eis2(3)-eis2(2))
            br(2) = 1.0+tau5/(eis2(2)+eps)
            br(3) = 1.0+tau5/(eis2(3)+eps)
            bs = br(2) + br(3)
            wr(2,i,m) = br(2)/bs
        end do
    end do
    
    do i=st,ed
        imin=1.0
        imax=0.0
        do m=nst,ned
            imin         = min(wl(2,i-3,m),wl(2,i-2,m),wl(2,i-1,m),wl(2,i,m),wl(2,i+1,m),wl(2,i+2,m),wl(2,i+3,m), &
                               wr(2,i-3,m),wr(2,i-2,m),wr(2,i-1,m),wr(2,i,m),wr(2,i+1,m),wr(2,i+2,m),wr(2,i+3,m),imin)                                
            imax         = max(wl(2,i-3,m),wl(2,i-2,m),wl(2,i-1,m),wl(2,i,m),wl(2,i+1,m),wl(2,i+2,m),wl(2,i+3,m), &
                               wr(2,i-3,m),wr(2,i-2,m),wr(2,i-1,m),wr(2,i,m),wr(2,i+1,m),wr(2,i+2,m),wr(2,i+3,m),imax)
        end do
        indicator    = imin/imax
        if(indicator>indicatorc) then
            sigma(i)=1.0
        else
            sigma(i)=0.0
        endif        
    end do
    
    id=0
    do i=st,ed
        YN=0
        if(sigma(i)<0.5) then
            YN=YN+1
        end if
        if(YN>0) then
            id=id+1
            ind(id)=i           
        end if
    end do
       
    do i=1-ngn,ni+ngn
        sigma(i)=0.0
    end do
    if(ed-st<9) then
        ind_st=st
        ind_ed=ed
        nsec=0
        
        nsecn=1
        sec_stn(1)=st
        sec_edn(1)=ed
    else
        if(id==0) then
            nsec=1
            sec_st(1)=st
            sec_ed(1)=ed
            do n=sec_st(1),sec_ed(1)
                sigma(n)=1.0
            end do
            ind_st=st
            ind_ed=ed
            nsecn=0
        else
            ind_st=st
            ind_ed=ed
            nsec=0
            nsecn=0
            if(ind(1)-ind_st>9) then
                nsec=1
                sec_st(1)=ind_st
                sec_ed(1)=ind(1)-1
                pre_ed   =sec_ed(1)
                do n=sec_st(1),sec_ed(1)
                    sigma(n)=1.0
                end do
            else
                nsecn=1
                sec_stn(1)=ind_st
                sec_edn(1)=ind(1)
                pre_ed    =sec_edn(1)
            end if
            do i=2,id
                if(ind(i)-ind(i-1)>9) then
                    nsec=nsec+1
                    sec_st(nsec)=ind(i-1)+1
                    sec_ed(nsec)=ind(i)-1
                    do n=sec_st(nsec),sec_ed(nsec)
                        sigma(n)=1.0
                    end do
                    nsecn=nsecn+1
                    sec_stn(nsecn)=pre_ed+1
                    sec_edn(nsecn)=ind(i-1)
                    pre_ed        =sec_ed(nsec)
                else
                    if(ind(i)-ind(i-1)>1 .or. i==id) then
                        nsecn=nsecn+1
                        sec_stn(nsecn)=pre_ed+1
                        sec_edn(nsecn)=ind(i)
                        pre_ed        =sec_edn(nsecn)
                    end if
                end if
            end do
            if(ind_ed-ind(id)>9) then
                nsec=nsec+1
                sec_st(nsec)=pre_ed+1
                sec_ed(nsec)=ind_ed
                do n=sec_st(nsec),sec_ed(nsec)
                    sigma(n)=1.0
                end do
            else
                nsecn=nsecn+1
                sec_stn(nsecn)=pre_ed+1
                sec_edn(nsecn)=ind_ed
            end if
        end if
    end if
    
    do m=nst,ned
        do j=1,nsecn
            do i=sec_stn(j),sec_edn(j)
                g(1,i,m) =      vn(i-2,m) - 3.0*vn(i-1,m) + 2.0*vn(i  ,m) 
                g(2,i,m) =     -vn(i  ,m) +     vn(i+1,m)
                g(3,i,m) =     -vn(i  ,m) +     vn(i+1,m)
                g(4,i,m) = -2.0*vn(i+1,m) + 3.0*vn(i+2,m) -     vn(i+3,m)

                s(1,i,m) = vn(i-2,m) - 2.0*vn(i-1,m) + vn(i,  m)
                s(2,i,m) = vn(i-1,m) - 2.0*vn(i  ,m) + vn(i+1,m) 
                s(3,i,m) = vn(i  ,m) - 2.0*vn(i+1,m) + vn(i+2,m)
                s(4,i,m) = vn(i+1,m) - 2.0*vn(i+2,m) + vn(i+3,m)
            end do
        end do      
    end do            
            
    do m=nst,ned
        do j=1,nsecn
            do i=sec_stn(j),sec_edn(j)
                do n=1,4
                    is = g(n,i,m)*g(n,i,m) + s(n,i,m)*s(n,i,m)
                    eis2(n) = is**2
                end do
                
                tau5 = abs(eis2(3)-eis2(1))
                
                do n=1,3
                    bl(n) = cl(n)*(1.0+tau5/(eis2(n)+eps))
                end do
                
                bs = bl(1) + bl(2) + bl(3)
                do n=1,3
                    wl(n,i,m) = bl(n)/bs
                end do
          
                tau5 = abs(eis2(4)-eis2(2))           
                do n=2,4
                    br(n) = cr(n)*(1.0+tau5/(eis2(n)+eps))
                end do

                bs = br(2) + br(3) + br(4)
                do n=2,4
                    wr(n,i,m) = br(n)/bs
                end do
            end do
        end do
    end do
    
!!    open(2,file='sigma.dat',form='formatted')
!!    write(2,*) 'variables="i" "sigmal" "sigmar"'
!!    do i=st,ed
!!
!!        write(2,200) i,sigmal(i,1),sigmar(i,1)
!!
!!    end do
!!
!!200 format(I4,1x,e12.5,1x,e12.5)
!!    close(2)
!!
!!    stop    
    
    do m=nst,ned
        do j=1,nsecn
            do i=sec_stn(j),sec_edn(j)
                vln(i,m) = vn(i  ,m) + ( wl(1,i,m)*( 4.0*g(1,i,m)-s(1,i,m)) + &
                                         wl(2,i,m)*( 4.0*g(2,i,m)-s(2,i,m)) + &
                                         wl(3,i,m)*( 4.0*g(3,i,m)-s(3,i,m)) )/8.0

                vrn(i,m) = vn(i+1,m) + ( wr(2,i,m)*(-4.0*g(2,i,m)-s(2,i,m)) + &
                                         wr(3,i,m)*(-4.0*g(3,i,m)-s(3,i,m)) + &
                                         wr(4,i,m)*(-4.0*g(4,i,m)-s(4,i,m)) )/8.0
            end do
        end do
    end do
    
    do j=1,nsec
        
        call hdcs5ciN(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,ve_lil,ve_lir,sub_limit,chkid,chklim,ndec,sec_st(j),sec_ed(j))        
        
    end do
        
    do m=nst,ned
        do i=st,ed
            
               vl(i,m) = sigma(i)*ve_lil(i,m) + vln(i,m)*(1.0-sigma(i))
                         
               vr(i,m) = sigma(i)*ve_lir(i,m) + vrn(i,m)*(1.0-sigma(i))            

        end do
    end do  

    do i=stc,edc
        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vl(i,m) < chklim(m) ) then
                    vl(i,m) = vn(i  ,m)
                    ndec = ndec + 1
                end if

                if (vr(i,m) < chklim(m) ) then
                    vr(i,m) = vn(i+1,m)
                    ndec = ndec + 1
                end if
           end if
        end do
    end do

end subroutine dcsh5pi

subroutine dcsh5picc(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_bound_scheme
    implicit none
    integer(kind_int), intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  , intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(in)  :: sn(1-ngn:ni+ngn,1:3)
    integer(kind_int), intent(in)  :: nfs,nfe
    integer(kind_int), intent(in)  :: nge
    real(kind_real)  , intent(out) :: vl(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(out) :: vr(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , external    :: sub_limit
    integer(kind_int), intent(in)  :: chkid(nst:ned)
    real(kind_real),   intent(in)  :: chklim(nst:ned)
    integer(kind_int), intent(out) :: ndec  
    real(kind_real)   :: vln(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: vrn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: ve_lil(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: ve_lir(1-ngn:ni+ngn,nst:ned)     
    real(kind_real)   :: s(1:4,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: g(1:4,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: wl(1:4,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: wr(1:4,1-ngn:ni+ngn,nst:ned)    
    real(kind_real)   :: aa,aa_scs
    real(kind_real)   :: sigma(1-ngn:ni+ngn),imax,imin,indicator
    real(kind_real), parameter :: cl(4) = (/1.D0/16.D0,10.D0/16.D0,5.D0/16.D0,0.D0/16.D0/)
    real(kind_real), parameter :: cr(4) = (/0.D0/16.D0,5.D0/16.D0,10.D0/16.D0,1.D0/16.D0/)  
    real(kind_real), parameter :: eps        = 1.0e-30
    real(kind_real), parameter :: indicatorc = 0.1
    real(kind_real)            :: is,eis2(4),bs,tau5,bl(4),br(4)
    integer(kind_int)          :: i,m,n,j,st,ed,stc,edc
    integer(kind_int)          :: id,YN,ind(1:ni+2*ngn),ind_st,ind_ed,pre_ed
    integer(kind_int)          :: nsec,sec_st(1:ni+2*ngn),sec_ed(1:ni+2*ngn)
    integer(kind_int)          :: nsecn,sec_stn(1:ni+2*ngn),sec_edn(1:ni+2*ngn)
    
    if (ni < 8) then
        call muscl2pvcc(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)
        return
    end if
    
    ndec = 0
    
    aa     = 0.01
    aa_scs = 0.001
    
    wl(:,:,:) = 0.3333
    wr(:,:,:) = 0.3333
    vln(:,:)    = 0.0
    ve_lil(:,:) = 0.0
    vrn(:,:)    = 0.0
    ve_lir(:,:) = 0.0    
           
    if (nfs == nbc_bound_scheme) then
        
        st  = 5
        stc = 1

        do m=nst,ned
            vl(1,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                  ! third-order
            vl(2,m) = (  3.0*vn(1,m) +  6.0*vn(2,m) -      vn(3,m) ) / 8.0   ! third-order 
            vl(3,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 + &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )
            vr(1,m) = vl(1,m)
            vr(2,m) = vl(2,m)
            vr(3,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 - &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )
            vl(4,m) =   ( 150.0*( vn(4,m) + vn(3,m) )          -  &
                           25.0*( vn(5,m) + vn(2,m) )          +  &
                            3.0*( vn(6,m) + vn(1,m) ) ) /256.0 -  &
                      aa*( 10.0*( vn(4,m) - vn(3,m) )          -  &
                            5.0*( vn(5,m) - vn(2,m) )          +  &
                                ( vn(6,m) - vn(1,m) ) ) 
            
            vr(4,m) =   ( 150.0*( vn(4,m) + vn(3,m) )          -  &
                           25.0*( vn(5,m) + vn(2,m) )          +  &
                            3.0*( vn(6,m) + vn(1,m) ) ) /256.0 +  &
                      aa*( 10.0*( vn(4,m) - vn(3,m) )          -  &
                            5.0*( vn(5,m) - vn(2,m) )          +  &
                                ( vn(6,m) - vn(1,m) ) )          
        end do
        
    else
        
        st  = -1
        stc = st
     
    end if
     
    if (nfe == nbc_bound_scheme) then
        
        ed  = ni - 3
        edc = ni + 1

        do m=nst,ned
            vr(ni+1,m) = ( 15.0*vn(ni,m) - 10.0*vn(ni-1,m) +  3.0*vn(ni-2,m) ) / 8.0                ! third-order
            vl(ni  ,m) = (  3.0*vn(ni,m) +  6.0*vn(ni-1,m) -      vn(ni-2,m) ) / 8.0   ! third-order
            vl(ni-1,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 + &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) )
            vl(ni+1,m) = vr(ni+1,m)
            vr(ni  ,m) = vl(ni  ,m)
            vr(ni-1,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 - &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) ) 
            vl(ni-2,m) =   ( 150.0*( vn(ni-2,m) + vn(ni-3,m) )          -  &
                              25.0*( vn(ni-1,m) + vn(ni-4,m) )          +  &
                               3.0*( vn(ni  ,m) + vn(ni-5,m) ) ) /256.0 -  &
                         aa*( 10.0*( vn(ni-2,m) - vn(ni-3,m) )          -  &
                               5.0*( vn(ni-1,m) - vn(ni-4,m) )          +  &
                                   ( vn(ni  ,m) - vn(ni-5,m) ) ) 
            
            vr(ni-2,m) =   ( 150.0*( vn(ni-2,m) + vn(ni-3,m) )          -  &
                              25.0*( vn(ni-1,m) + vn(ni-4,m) )          +  &
                               3.0*( vn(ni  ,m) + vn(ni-5,m) ) ) /256.0 +  &
                         aa*( 10.0*( vn(ni-2,m) - vn(ni-3,m) )          -  &
                               5.0*( vn(ni-1,m) - vn(ni-4,m) )          +  &
                                   ( vn(ni  ,m) - vn(ni-5,m) ) )              
        end do
        
    else
        
        ed  = ni + 3
        edc = ed
            
    end if
    
    do m=nst,ned
        do i=st-3,ed+1
            g(1,i,m) = vn(i  ,m) - vn(i-1,m)            
            g(2,i,m) = vn(i+1,m) - vn(i  ,m)
            g(3,i,m) = vn(i+2,m) - vn(i+1,m)            
        end do
    end do   
    
    do m=nst,ned
        do i=st-3,ed+1
            do n=1,3
                is = g(n,i,m)*g(n,i,m)
                eis2(n) = is**2
            end do

            tau5 = abs(eis2(2)-eis2(1))
            bl(1) = 1.0+tau5/(eis2(1)+eps)
            bl(2) = 1.0+tau5/(eis2(2)+eps)
            bs = bl(1) + bl(2)
            wl(2,i,m) = bl(2)/bs
            
            tau5 = abs(eis2(3)-eis2(2))
            br(2) = 1.0+tau5/(eis2(2)+eps)
            br(3) = 1.0+tau5/(eis2(3)+eps)
            bs = br(2) + br(3)
            wr(2,i,m) = br(2)/bs
        end do
    end do
    
    do i=st-1,ed-1
        imin=1.0
        imax=0.0
        do m=nst,ned
            imin         = min(wl(2,i-3,m),wl(2,i-2,m),wl(2,i-1,m),wl(2,i,m),wl(2,i+1,m),wl(2,i+2,m),wl(2,i+3,m), &
                               wr(2,i-3,m),wr(2,i-2,m),wr(2,i-1,m),wr(2,i,m),wr(2,i+1,m),wr(2,i+2,m),wr(2,i+3,m),imin)                                
            imax         = max(wl(2,i-3,m),wl(2,i-2,m),wl(2,i-1,m),wl(2,i,m),wl(2,i+1,m),wl(2,i+2,m),wl(2,i+3,m), &
                               wr(2,i-3,m),wr(2,i-2,m),wr(2,i-1,m),wr(2,i,m),wr(2,i+1,m),wr(2,i+2,m),wr(2,i+3,m),imax)
        end do
        indicator    = imin/imax
        if(indicator>indicatorc) then
            sigma(i)=1.0
        else
            sigma(i)=0.0
        endif        
    end do
    
    id=0
    do i=st-1,ed-1
        YN=0
        if(sigma(i)<0.5) then
            YN=YN+1
        end if
        if(YN>0) then
            id=id+1
            ind(id)=i           
        end if
    end do
       
    do i=1-ngn,ni+ngn
        sigma(i)=0.0
    end do
    if(ed-st<9) then
        ind_st=st
        ind_ed=ed
        nsec=0
        
        nsecn=1
        sec_stn(1)=st-1
        sec_edn(1)=ed-1
    else
        if(id==0) then
            nsec=1
            sec_st(1)=st-1
            sec_ed(1)=ed-1
            do n=sec_st(1),sec_ed(1)
                sigma(n)=1.0
            end do
            ind_st=st-1
            ind_ed=ed-1
            nsecn=0
        else
            ind_st=st-1
            ind_ed=ed-1
            nsec=0
            nsecn=0
            if(ind(1)-ind_st>9) then
                nsec=1
                sec_st(1)=ind_st-1
                sec_ed(1)=ind(1)-1
                pre_ed   =sec_ed(1)
                do n=sec_st(1),sec_ed(1)
                    sigma(n)=1.0
                end do
            else
                nsecn=1
                sec_stn(1)=ind_st
                sec_edn(1)=ind(1)
                pre_ed    =sec_edn(1)
            end if
            do i=2,id
                if(ind(i)-ind(i-1)>9) then
                    nsec=nsec+1
                    sec_st(nsec)=ind(i-1)+1
                    sec_ed(nsec)=ind(i)-1
                    do n=sec_st(nsec),sec_ed(nsec)
                        sigma(n)=1.0
                    end do
                    nsecn=nsecn+1
                    sec_stn(nsecn)=pre_ed+1
                    sec_edn(nsecn)=ind(i-1)
                    pre_ed        =sec_ed(nsec)
                else
                    if(ind(i)-ind(i-1)>1 .or. i==id) then
                        nsecn=nsecn+1
                        sec_stn(nsecn)=pre_ed+1
                        sec_edn(nsecn)=ind(i)
                        pre_ed        =sec_edn(nsecn)
                    end if
                end if
            end do
            if(ind_ed-ind(id)>9) then
                nsec=nsec+1
                sec_st(nsec)=pre_ed+1
                sec_ed(nsec)=ind_ed
                do n=sec_st(nsec),sec_ed(nsec)
                    sigma(n)=1.0
                end do
            else
                nsecn=nsecn+1
                sec_stn(nsecn)=pre_ed+1
                sec_edn(nsecn)=ind_ed
            end if
        end if
    end if
    
    do m=nst,ned
        do j=1,nsecn
            do i=sec_stn(j),sec_edn(j)
                g(1,i,m) =      vn(i-2,m) - 3.0*vn(i-1,m) + 2.0*vn(i  ,m) 
                g(2,i,m) =     -vn(i  ,m) +     vn(i+1,m)
                g(3,i,m) =     -vn(i  ,m) +     vn(i+1,m)
                g(4,i,m) = -2.0*vn(i+1,m) + 3.0*vn(i+2,m) -     vn(i+3,m)

                s(1,i,m) = vn(i-2,m) - 2.0*vn(i-1,m) + vn(i,  m)
                s(2,i,m) = vn(i-1,m) - 2.0*vn(i  ,m) + vn(i+1,m) 
                s(3,i,m) = vn(i  ,m) - 2.0*vn(i+1,m) + vn(i+2,m)
                s(4,i,m) = vn(i+1,m) - 2.0*vn(i+2,m) + vn(i+3,m)
            end do
        end do      
    end do            
            
    do m=nst,ned
        do j=1,nsecn
            do i=sec_stn(j),sec_edn(j)
                do n=1,4
                    is = g(n,i,m)*g(n,i,m) + s(n,i,m)*s(n,i,m)
                    eis2(n) = is**2
                end do
                
                tau5 = abs(eis2(3)-eis2(1))
                
                do n=1,3
                    bl(n) = cl(n)*(1.0+tau5/(eis2(n)+eps))
                end do
                
                bs = bl(1) + bl(2) + bl(3)
                do n=1,3
                    wl(n,i,m) = bl(n)/bs
                end do
          
                tau5 = abs(eis2(4)-eis2(2))           
                do n=2,4
                    br(n) = cr(n)*(1.0+tau5/(eis2(n)+eps))
                end do

                bs = br(2) + br(3) + br(4)
                do n=2,4
                    wr(n,i,m) = br(n)/bs
                end do
            end do
        end do
    end do
    
!!    open(2,file='sigma.dat',form='formatted')
!!    write(2,*) 'variables="i" "sigmal" "sigmar"'
!!    do i=st,ed
!!
!!        write(2,200) i,sigmal(i,1),sigmar(i,1)
!!
!!    end do
!!
!!200 format(I4,1x,e12.5,1x,e12.5)
!!    close(2)
!!
!!    stop    
    
    do m=nst,ned
        do j=1,nsecn
            do i=sec_stn(j),sec_edn(j)
                vln(i,m) = vn(i  ,m) + ( wl(1,i,m)*( 4.0*g(1,i,m)-s(1,i,m)) + &
                                         wl(2,i,m)*( 4.0*g(2,i,m)-s(2,i,m)) + &
                                         wl(3,i,m)*( 4.0*g(3,i,m)-s(3,i,m)) )/8.0

                vrn(i,m) = vn(i+1,m) + ( wr(2,i,m)*(-4.0*g(2,i,m)-s(2,i,m)) + &
                                         wr(3,i,m)*(-4.0*g(3,i,m)-s(3,i,m)) + &
                                         wr(4,i,m)*(-4.0*g(4,i,m)-s(4,i,m)) )/8.0
            end do
        end do
    end do
    
    do j=1,nsec
        
        call hdcs5ciN(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,ve_lil,ve_lir,sub_limit,chkid,chklim,ndec,sec_st(j),sec_ed(j))        
        
    end do
        
    do m=nst,ned
        do i=st-1,ed-1
            
               vl(i+1,m) = sigma(i)*ve_lil(i,m) + vln(i,m)*(1.0-sigma(i))
                         
               vr(i+1,m) = sigma(i)*ve_lir(i,m) + vrn(i,m)*(1.0-sigma(i))            

        end do
    end do  

    do i=stc,edc
        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vl(i,m) < chklim(m) ) then
                    vl(i,m) = vn(i  ,m)
                    ndec = ndec + 1
                end if

                if (vr(i,m) < chklim(m) ) then
                    vr(i,m) = vn(i+1,m)
                    ndec = ndec + 1
                end if
           end if
        end do
    end do

end subroutine dcsh5picc

subroutine scsn3pi(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_bound_scheme
    implicit none
    integer(kind_int), intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  , intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(in)  :: sn(1-ngn:ni+ngn,1:3)
    integer(kind_int), intent(in)  :: nfs,nfe
    integer(kind_int), intent(in)  :: nge
    real(kind_real)  , intent(out) :: vl(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(out) :: vr(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , external    :: sub_limit
    integer(kind_int), intent(in)  :: chkid(nst:ned)
    real(kind_real),   intent(in)  :: chklim(nst:ned)
    integer(kind_int), intent(out) :: ndec 
    real(kind_real)   :: s(1:4,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: g(1:4,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: wl(1:4,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: wr(1:4,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: dl3(-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: dr3(-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: aa
    real(kind_real), parameter :: cl(4) = (/1.D0/16.D0,10.D0/16.D0,5.D0/16.D0,0.D0/16.D0/)
    real(kind_real), parameter :: cr(4) = (/0.D0/16.D0,5.D0/16.D0,10.D0/16.D0,1.D0/16.D0/)
    real(kind_real), parameter :: eps        = 1.0e-30
    real(kind_real)            :: is,eis2(4),bs,tau5,bl(4),br(4)    
    integer(kind_int)          :: i,m,n,st,ed,stc,edc
    
    if (ni < 6) then
        call muscl2pv(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)
        return
    end if
    
    ndec = 0
    
    aa     = 0.01
    
    wl(:,:,:) = 0.33333
    wr(:,:,:) = 0.33333
           
    if (nfs == nbc_bound_scheme) then
        
        st  = 3
        stc = 0

        do m=nst,ned
            vl(0,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                  ! third-order
            vl(1,m) = (  3.0*vn(1,m) +  6.0*vn(2,m) -      vn(3,m) ) / 8.0   ! third-order 

            vl(2,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 + &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )
            vr(0,m) = vl(0,m)
            vr(1,m) = vl(1,m)

            vr(2,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 - &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )         
        end do
        
    else
        
        st  = -1
        stc = st
     
    end if
     
    if (nfe == nbc_bound_scheme) then
        
        ed  = ni - 3
        edc = ni

        do m=nst,ned
            vr(ni  ,m) = ( 15.0*vn(ni,m) - 10.0*vn(ni-1,m) +  3.0*vn(ni-2,m) ) / 8.0                ! third-order
            vl(ni-1,m) = (  3.0*vn(ni,m) +  6.0*vn(ni-1,m) -      vn(ni-2,m) ) / 8.0   ! third-order
            vl(ni-2,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 + &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) )
            vl(ni  ,m) = vr(ni  ,m)
            vr(ni-1,m) = vl(ni-1,m)

            vr(ni-2,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 - &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) )           
        end do
        
    else
        
        ed  = ni + 1
        edc = ed
            
    end if
    
    do m=nst,ned
        do i=st,ed
            g(1,i,m) =      vn(i-2,m) - 3.0*vn(i-1,m) + 2.0*vn(i  ,m) 
            g(2,i,m) =     -vn(i  ,m) +     vn(i+1,m)
            g(3,i,m) =     -vn(i  ,m) +     vn(i+1,m)
            g(4,i,m) = -2.0*vn(i+1,m) + 3.0*vn(i+2,m) -     vn(i+3,m)

            s(1,i,m) = vn(i-2,m) - 2.0*vn(i-1,m) + vn(i,  m)
            s(2,i,m) = vn(i-1,m) - 2.0*vn(i  ,m) + vn(i+1,m) 
            s(3,i,m) = vn(i  ,m) - 2.0*vn(i+1,m) + vn(i+2,m)
            s(4,i,m) = vn(i+1,m) - 2.0*vn(i+2,m) + vn(i+3,m)
        end do     
    end do            
            
    do m=nst,ned
        do i=st,ed
            do n=1,4
                is = g(n,i,m)*g(n,i,m) + s(n,i,m)*s(n,i,m)
                eis2(n) = is**2
            end do
            
            tau5 = abs(eis2(3)-eis2(1))
            
            do n=1,3
                bl(n) = cl(n)*(1.0+tau5/(eis2(n)+eps))
            end do
            
            bs = bl(1) + bl(2) + bl(3)
            do n=1,3
                wl(n,i,m) = bl(n)/bs
            end do
          
            tau5 = abs(eis2(4)-eis2(2))           
            do n=2,4
                br(n) = cr(n)*(1.0+tau5/(eis2(n)+eps))
            end do

            bs = br(2) + br(3) + br(4)
            do n=2,4
                wr(n,i,m) = br(n)/bs
            end do
        end do
    end do
      
    do m=nst,ned
        do i=st,ed
            vl(i,m) = vn(i  ,m) + ( wl(1,i,m)*( 4.0*g(1,i,m)-s(1,i,m)) + &
                                    wl(2,i,m)*( 4.0*g(2,i,m)-s(2,i,m)) + &
                                    wl(3,i,m)*( 4.0*g(3,i,m)-s(3,i,m)) )/8.0
        end do

        do i=st,ed
            vr(i,m) = vn(i+1,m) + ( wr(2,i,m)*(-4.0*g(2,i,m)-s(2,i,m)) + &
                                    wr(3,i,m)*(-4.0*g(3,i,m)-s(3,i,m)) + &
                                    wr(4,i,m)*(-4.0*g(4,i,m)-s(4,i,m)) )/8.0
        end do
    end do

    do i=stc,edc
        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vl(i,m) < chklim(m) ) then
                    vl(i,m) = vn(i  ,m)
                    ndec = ndec + 1
                end if

                if (vr(i,m) < chklim(m) ) then
                    vr(i,m) = vn(i+1,m)
                    ndec = ndec + 1
                end if
           end if
        end do
    end do

end subroutine scsn3pi

subroutine scsn3picc(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_bound_scheme
    implicit none
    integer(kind_int), intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  , intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(in)  :: sn(1-ngn:ni+ngn,1:3)
    integer(kind_int), intent(in)  :: nfs,nfe
    integer(kind_int), intent(in)  :: nge
    real(kind_real)  , intent(out) :: vl(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(out) :: vr(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , external    :: sub_limit
    integer(kind_int), intent(in)  :: chkid(nst:ned)
    real(kind_real),   intent(in)  :: chklim(nst:ned)
    integer(kind_int), intent(out) :: ndec 
    real(kind_real)   :: s(1:4,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: g(1:4,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: wl(1:4,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: wr(1:4,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: dl3(-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: dr3(-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: aa
    real(kind_real), parameter :: cl(4) = (/1.D0/16.D0,10.D0/16.D0,5.D0/16.D0,0.D0/16.D0/)
    real(kind_real), parameter :: cr(4) = (/0.D0/16.D0,5.D0/16.D0,10.D0/16.D0,1.D0/16.D0/)
    real(kind_real), parameter :: eps        = 1.0e-30
    real(kind_real)            :: is,eis2(4),bs,tau5,bl(4),br(4)    
    integer(kind_int)          :: i,m,n,st,ed,stc,edc
    
    if (ni < 6) then
        call muscl2pvcc(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)
        return
    end if
    
    ndec = 0
    
    aa     = 0.01
    
    wl(:,:,:) = 0.33333
    wr(:,:,:) = 0.33333
           
    if (nfs == nbc_bound_scheme) then
        
        st  = 4
        stc = 1

        do m=nst,ned
            vl(1,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                  ! third-order
            vl(2,m) = (  3.0*vn(1,m) +  6.0*vn(2,m) -      vn(3,m) ) / 8.0   ! third-order 

            vl(3,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 + &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )
            vr(1,m) = vl(1,m)
            vr(2,m) = vl(2,m)

            vr(3,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 - &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )         
        end do
        
    else
        
        st  = 0
        stc = st
     
    end if
     
    if (nfe == nbc_bound_scheme) then
        
        ed  = ni - 2
        edc = ni + 1

        do m=nst,ned
            vr(ni+1,m) = ( 15.0*vn(ni,m) - 10.0*vn(ni-1,m) +  3.0*vn(ni-2,m) ) / 8.0                ! third-order
            vl(ni  ,m) = (  3.0*vn(ni,m) +  6.0*vn(ni-1,m) -      vn(ni-2,m) ) / 8.0   ! third-order
            vl(ni-1,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 + &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) )
            vl(ni+1,m) = vr(ni+1,m)
            vr(ni  ,m) = vl(ni  ,m)

            vr(ni-1,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 - &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) )           
        end do
        
    else
        
        ed  = ni + 2
        edc = ed
            
    end if
    
    do m=nst,ned
        do i=st-1,ed-1
            g(1,i,m) =      vn(i-2,m) - 3.0*vn(i-1,m) + 2.0*vn(i  ,m) 
            g(2,i,m) =     -vn(i  ,m) +     vn(i+1,m)
            g(3,i,m) =     -vn(i  ,m) +     vn(i+1,m)
            g(4,i,m) = -2.0*vn(i+1,m) + 3.0*vn(i+2,m) -     vn(i+3,m)

            s(1,i,m) = vn(i-2,m) - 2.0*vn(i-1,m) + vn(i,  m)
            s(2,i,m) = vn(i-1,m) - 2.0*vn(i  ,m) + vn(i+1,m) 
            s(3,i,m) = vn(i  ,m) - 2.0*vn(i+1,m) + vn(i+2,m)
            s(4,i,m) = vn(i+1,m) - 2.0*vn(i+2,m) + vn(i+3,m)
        end do     
    end do            
            
    do m=nst,ned
        do i=st-1,ed-1
            do n=1,4
                is = g(n,i,m)*g(n,i,m) + s(n,i,m)*s(n,i,m)
                eis2(n) = is**2
            end do
            
            tau5 = abs(eis2(3)-eis2(1))
            
            do n=1,3
                bl(n) = cl(n)*(1.0+tau5/(eis2(n)+eps))
            end do
            
            bs = bl(1) + bl(2) + bl(3)
            do n=1,3
                wl(n,i,m) = bl(n)/bs
            end do
          
            tau5 = abs(eis2(4)-eis2(2))           
            do n=2,4
                br(n) = cr(n)*(1.0+tau5/(eis2(n)+eps))
            end do

            bs = br(2) + br(3) + br(4)
            do n=2,4
                wr(n,i,m) = br(n)/bs
            end do
        end do
    end do
      
    do m=nst,ned
        do i=st-1,ed-1
            vl(i+1,m) = vn(i  ,m) + ( wl(1,i,m)*( 4.0*g(1,i,m)-s(1,i,m)) + &
                                      wl(2,i,m)*( 4.0*g(2,i,m)-s(2,i,m)) + &
                                      wl(3,i,m)*( 4.0*g(3,i,m)-s(3,i,m)) )/8.0
        end do

        do i=st-1,ed-1
            vr(i+1,m) = vn(i+1,m) + ( wr(2,i,m)*(-4.0*g(2,i,m)-s(2,i,m)) + &
                                      wr(3,i,m)*(-4.0*g(3,i,m)-s(3,i,m)) + &
                                      wr(4,i,m)*(-4.0*g(4,i,m)-s(4,i,m)) )/8.0
        end do
    end do

    do i=stc,edc
        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vl(i,m) < chklim(m) ) then
                    vl(i,m) = vn(i  ,m)
                    ndec = ndec + 1
                end if

                if (vr(i,m) < chklim(m) ) then
                    vr(i,m) = vn(i+1,m)
                    ndec = ndec + 1
                end if
           end if
        end do
    end do

end subroutine scsn3picc

subroutine scsh3ci(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_bound_scheme
    implicit none
    integer(kind_int), intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  , intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(in)  :: sn(1-ngn:ni+ngn,1:3)
    integer(kind_int), intent(in)  :: nfs,nfe
    integer(kind_int), intent(in)  :: nge
    real(kind_real)  , intent(out) :: vl(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(out) :: vr(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , external    :: sub_limit
    integer(kind_int), intent(in)  :: chkid(nst:ned)
    real(kind_real),   intent(in)  :: chklim(nst:ned)
    integer(kind_int), intent(out) :: ndec  
    real(kind_real)   :: vln(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: vrn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: ve_lil(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: ve_lir(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: d2_li(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: d3_li(1-ngn:ni+ngn,nst:ned)  
    real(kind_real)   :: nt,nx,ny,nz,fp(nst:ned)
    real(kind_real)   :: vp(nst:ned),vc(nst:ned)
    real(kind_real)   :: q(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: q1(1-ngn:ni+ngn,nst:ned)    
    real(kind_real)   :: s(1:4,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: g(1:4,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: wl(1:4,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: wr(1:4,1-ngn:ni+ngn,nst:ned)    
    real(kind_real)   :: aa,aa_scs
    real(kind_real)   :: sigma(1-ngn:ni+ngn),imax,imin,indicator
    real(kind_real), parameter :: cl(4) = (/1.D0/16.D0,10.D0/16.D0,5.D0/16.D0,0.D0/16.D0/)
    real(kind_real), parameter :: cr(4) = (/0.D0/16.D0,5.D0/16.D0,10.D0/16.D0,1.D0/16.D0/)      
    real(kind_real), parameter :: eps        = 1.0e-30
    real(kind_real), parameter :: indicatorc = 0.1
    real(kind_real)            :: is,eis2(4),bs,tau5,bl(4),br(4)
    integer(kind_int)          :: i,m,n,j,st,ed,stc,edc
    integer(kind_int)          :: id,YN,ind(1:ni+2*ngn),ind_st,ind_ed,pre_ed
    integer(kind_int)          :: nsec,sec_st(1:ni+2*ngn),sec_ed(1:ni+2*ngn)
    integer(kind_int)          :: nsecn,sec_stn(1:ni+2*ngn),sec_edn(1:ni+2*ngn)
    
    if (ni < 8) then
        call muscl2pv(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)
        return
    end if
    
    ndec = 0
    
    aa     = 0.01
    aa_scs = 0.001
    
    wl(:,:,:) = 0.3333
    wr(:,:,:) = 0.3333
    vln(:,:)    = 0.0
    ve_lil(:,:) = 0.0
    vrn(:,:)    = 0.0
    ve_lir(:,:) = 0.0    
           
    if (nfs == nbc_bound_scheme) then
        
        st  = 4
        stc = 0

        do m=nst,ned
            vl(0,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                  ! third-order
            vl(1,m) = (  3.0*vn(1,m) +  6.0*vn(2,m) -      vn(3,m) ) / 8.0   ! third-order 
            vl(2,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 + &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )
            vr(0,m) = vl(0,m)
            vr(1,m) = vl(1,m)
            vr(2,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 - &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )
            vl(3,m) =   ( 150.0*( vn(4,m) + vn(3,m) )          -  &
                           25.0*( vn(5,m) + vn(2,m) )          +  &
                            3.0*( vn(6,m) + vn(1,m) ) ) /256.0 -  &
                      aa*( 10.0*( vn(4,m) - vn(3,m) )          -  &
                            5.0*( vn(5,m) - vn(2,m) )          +  &
                                ( vn(6,m) - vn(1,m) ) ) 
            
            vr(3,m) =   ( 150.0*( vn(4,m) + vn(3,m) )          -  &
                           25.0*( vn(5,m) + vn(2,m) )          +  &
                            3.0*( vn(6,m) + vn(1,m) ) ) /256.0 +  &
                      aa*( 10.0*( vn(4,m) - vn(3,m) )          -  &
                            5.0*( vn(5,m) - vn(2,m) )          +  &
                                ( vn(6,m) - vn(1,m) ) )          
        end do
        
    else
        
        st  = -2
        stc = st
     
    end if
     
    if (nfe == nbc_bound_scheme) then
        
        ed  = ni - 4
        edc = ni

        do m=nst,ned
            vr(ni  ,m) = ( 15.0*vn(ni,m) - 10.0*vn(ni-1,m) +  3.0*vn(ni-2,m) ) / 8.0                ! third-order
            vl(ni-1,m) = (  3.0*vn(ni,m) +  6.0*vn(ni-1,m) -      vn(ni-2,m) ) / 8.0   ! third-order
            vl(ni-2,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 + &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) )
            vl(ni  ,m) = vr(ni  ,m)
            vr(ni-1,m) = vl(ni-1,m)
            vr(ni-2,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 - &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) ) 
            vl(ni-3,m) =   ( 150.0*( vn(ni-2,m) + vn(ni-3,m) )          -  &
                              25.0*( vn(ni-1,m) + vn(ni-4,m) )          +  &
                               3.0*( vn(ni  ,m) + vn(ni-5,m) ) ) /256.0 -  &
                         aa*( 10.0*( vn(ni-2,m) - vn(ni-3,m) )          -  &
                               5.0*( vn(ni-1,m) - vn(ni-4,m) )          +  &
                                   ( vn(ni  ,m) - vn(ni-5,m) ) ) 
            
            vr(ni-3,m) =   ( 150.0*( vn(ni-2,m) + vn(ni-3,m) )          -  &
                              25.0*( vn(ni-1,m) + vn(ni-4,m) )          +  &
                               3.0*( vn(ni  ,m) + vn(ni-5,m) ) ) /256.0 +  &
                         aa*( 10.0*( vn(ni-2,m) - vn(ni-3,m) )          -  &
                               5.0*( vn(ni-1,m) - vn(ni-4,m) )          +  &
                                   ( vn(ni  ,m) - vn(ni-5,m) ) )              
        end do
        
    else
        
        ed  = ni + 2
        edc = ed
            
    end if
    
    do m=nst,ned
        do i=st-2,ed+2
            g(1,i,m) = vn(i  ,m) - vn(i-1,m)            
            g(2,i,m) = vn(i+1,m) - vn(i  ,m)
            g(3,i,m) = vn(i+2,m) - vn(i+1,m)            
        end do
    end do   
    
    do m=nst,ned
        do i=st-2,ed+2
            do n=1,3
                is = g(n,i,m)*g(n,i,m)
                eis2(n) = is**2
            end do

            tau5 = abs(eis2(2)-eis2(1))
            bl(1) = 1.0+tau5/(eis2(1)+eps)
            bl(2) = 1.0+tau5/(eis2(2)+eps)
            bs = bl(1) + bl(2)
            wl(2,i,m) = bl(2)/bs
            
            tau5 = abs(eis2(3)-eis2(2))
            br(2) = 1.0+tau5/(eis2(2)+eps)
            br(3) = 1.0+tau5/(eis2(3)+eps)
            bs = br(2) + br(3)
            wr(2,i,m) = br(2)/bs
        end do
    end do
    
    do i=st,ed
        imin=1.0
        imax=0.0
        do m=nst,ned
            imin         = min(wl(2,i-3,m),wl(2,i-2,m),wl(2,i-1,m),wl(2,i,m),wl(2,i+1,m),wl(2,i+2,m),wl(2,i+3,m), &
                               wr(2,i-3,m),wr(2,i-2,m),wr(2,i-1,m),wr(2,i,m),wr(2,i+1,m),wr(2,i+2,m),wr(2,i+3,m),imin)                                
            imax         = max(wl(2,i-3,m),wl(2,i-2,m),wl(2,i-1,m),wl(2,i,m),wl(2,i+1,m),wl(2,i+2,m),wl(2,i+3,m), &
                               wr(2,i-3,m),wr(2,i-2,m),wr(2,i-1,m),wr(2,i,m),wr(2,i+1,m),wr(2,i+2,m),wr(2,i+3,m),imax)
        end do
        indicator    = imin/imax
        if(indicator>indicatorc) then
            sigma(i)=1.0
        else
            sigma(i)=0.0
        endif        
    end do
    
    id=0
    do i=st,ed
        YN=0
        if(sigma(i)<0.5) then
            YN=YN+1
        end if
        if(YN>0) then
            id=id+1
            ind(id)=i           
        end if
    end do
       
    do i=1-ngn,ni+ngn
        sigma(i)=0.0
    end do
    if(ed-st<9) then
        ind_st=st
        ind_ed=ed
        nsec=0
        
        nsecn=1
        sec_stn(1)=st
        sec_edn(1)=ed
    else
        if(id==0) then
            nsec=1
            sec_st(1)=st
            sec_ed(1)=ed
            do n=sec_st(1),sec_ed(1)
                sigma(n)=1.0
            end do
            ind_st=st
            ind_ed=ed
            nsecn=0
        else
            ind_st=st
            ind_ed=ed
            nsec=0
            nsecn=0
            if(ind(1)-ind_st>9) then
                nsec=1
                sec_st(1)=ind_st
                sec_ed(1)=ind(1)-1
                pre_ed   =sec_ed(1)
                do n=sec_st(1),sec_ed(1)
                    sigma(n)=1.0
                end do
            else
                nsecn=1
                sec_stn(1)=ind_st
                sec_edn(1)=ind(1)
                pre_ed    =sec_edn(1)
            end if
            do i=2,id
                if(ind(i)-ind(i-1)>9) then
                    nsec=nsec+1
                    sec_st(nsec)=ind(i-1)+1
                    sec_ed(nsec)=ind(i)-1
                    do n=sec_st(nsec),sec_ed(nsec)
                        sigma(n)=1.0
                    end do
                    nsecn=nsecn+1
                    sec_stn(nsecn)=pre_ed+1
                    sec_edn(nsecn)=ind(i-1)
                    pre_ed        =sec_ed(nsec)
                else
                    if(ind(i)-ind(i-1)>1 .or. i==id) then
                        nsecn=nsecn+1
                        sec_stn(nsecn)=pre_ed+1
                        sec_edn(nsecn)=ind(i)
                        pre_ed        =sec_edn(nsecn)
                    end if
                end if
            end do
            if(ind_ed-ind(id)>9) then
                nsec=nsec+1
                sec_st(nsec)=pre_ed+1
                sec_ed(nsec)=ind_ed
                do n=sec_st(nsec),sec_ed(nsec)
                    sigma(n)=1.0
                end do
            else
                nsecn=nsecn+1
                sec_stn(nsecn)=pre_ed+1
                sec_edn(nsecn)=ind_ed
            end if
        end if
    end if
    
    do m=nst,ned
        do j=1,nsecn
            do i=sec_stn(j),sec_edn(j)
                g(1,i,m) =      vn(i-2,m) - 3.0*vn(i-1,m) + 2.0*vn(i  ,m) 
                g(2,i,m) =     -vn(i  ,m) +     vn(i+1,m)
                g(3,i,m) =     -vn(i  ,m) +     vn(i+1,m)
                g(4,i,m) = -2.0*vn(i+1,m) + 3.0*vn(i+2,m) -     vn(i+3,m)

                s(1,i,m) = vn(i-2,m) - 2.0*vn(i-1,m) + vn(i,  m)
                s(2,i,m) = vn(i-1,m) - 2.0*vn(i  ,m) + vn(i+1,m) 
                s(3,i,m) = vn(i  ,m) - 2.0*vn(i+1,m) + vn(i+2,m)
                s(4,i,m) = vn(i+1,m) - 2.0*vn(i+2,m) + vn(i+3,m)
            end do
        end do      
    end do
    
    do j=1,nsecn
        do i=sec_stn(j),sec_edn(j)
            nx = (sn(i,1) + sn(i+1,1))/2.0
            ny = (sn(i,2) + sn(i+1,2))/2.0
            nz = (sn(i,3) + sn(i+1,3))/2.0

            do m=nst,ned
                vp(m) = (vn(i,m) + vn(i+1,m))/2.0
                fp(m) = vn(i,m)
            end do
            call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
            do m=nst,ned
                q(i,m) = vc(m)
            end do
            
            do m=nst,ned
                fp(m) = vn(i+1,m)
            end do
            call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
            do m=nst,ned
                q1(i+1,m) = vc(m)
            end do        
             
            do n=1,4
                do m=nst,ned
                    fp(m) = g(n,i,m)
                end do
                call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
                do m=nst,ned
                    g(n,i,m) = vc(m)
                end do

                do m=nst,ned
                    fp(m) = s(n,i,m)
                end do
                call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
                do m=nst,ned
                    s(n,i,m) = vc(m)
                end do          
            end do
        end do
    end do   
            
            
    do m=nst,ned
        do j=1,nsecn
            do i=sec_stn(j),sec_edn(j)
                do n=1,4
                    is = g(n,i,m)*g(n,i,m) + s(n,i,m)*s(n,i,m)
                    eis2(n) = is**2
                end do
                
                tau5 = abs(eis2(3)-eis2(1))
                
                do n=1,3
                    bl(n) = cl(n)*(1.0+tau5/(eis2(n)+eps))
                end do
                
                bs = bl(1) + bl(2) + bl(3)
                do n=1,3
                    wl(n,i,m) = bl(n)/bs
                end do
          
                tau5 = abs(eis2(4)-eis2(2))           
                do n=2,4
                    br(n) = cr(n)*(1.0+tau5/(eis2(n)+eps))
                end do

                bs = br(2) + br(3) + br(4)
                do n=2,4
                    wr(n,i,m) = br(n)/bs
                end do
            end do
        end do
    end do
    
!!    open(2,file='sigma.dat',form='formatted')
!!    write(2,*) 'variables="i" "sigmal" "sigmar"'
!!    do i=st,ed
!!
!!        write(2,200) i,sigmal(i,1),sigmar(i,1)
!!
!!    end do
!!
!!200 format(I4,1x,e12.5,1x,e12.5)
!!    close(2)
!!
!!    stop    
    
    do m=nst,ned
        do j=1,nsecn
            do i=sec_stn(j),sec_edn(j)
                vln(i,m) =  q(i  ,m) + ( wl(1,i,m)*( 4.0*g(1,i,m)-s(1,i,m)) + &
                                         wl(2,i,m)*( 4.0*g(2,i,m)-s(2,i,m)) + &
                                         wl(3,i,m)*( 4.0*g(3,i,m)-s(3,i,m)) )/8.0

                vrn(i,m) = q1(i+1,m) + ( wr(2,i,m)*(-4.0*g(2,i,m)-s(2,i,m)) + &
                                         wr(3,i,m)*(-4.0*g(3,i,m)-s(3,i,m)) + &
                                         wr(4,i,m)*(-4.0*g(4,i,m)-s(4,i,m)) )/8.0
            end do
        end do
    end do
    
    do j=1,nsecn
        do i=sec_stn(j),sec_edn(j)
            nx = (sn(i,1) + sn(i+1,1))/2.0
            ny = (sn(i,2) + sn(i+1,2))/2.0
            nz = (sn(i,3) + sn(i+1,3))/2.0

            do m=nst,ned
                vp(m) = (vn(i,m) + vn(i+1,m))/2.0
                vc(m) = vln(i,m)
            end do

            call RxQ_p(nst,ned,vp,nt,nx,ny,nz,vc,fp) 
            do m=nst,ned
                vln(i,m) = fp(m)
            end do
            
            do m=nst,ned
                vc(m) = vrn(i ,m)
            end do

            call RxQ_p(nst,ned,vp,nt,nx,ny,nz,vc,fp) 
            do m=nst,ned
                vrn(i,m) = fp(m)
            end do
        end do
    end do     
      
    do j=1,nsec
        call d3_via_T2N(ni,nst,ned,ngn,vn,nfs,nfe,nge,d3_li,sec_st(j),sec_ed(j))
        
        call vec2_via_T2N(ni,nst,ned,ngn,vn,nfs,nfe,nge,d2_li,sec_st(j),sec_ed(j))
    end do
    
    do m=nst,ned
        do j=1,nsec
            do i=sec_st(j),sec_ed(j)
                
                ve_lil(i,m) = 0.5*(vn(i+1,m) + vn(i  ,m)) - 1.0/8.0*d2_li(i,m) +  aa_scs*d3_li(i,m)
                
                ve_lir(i,m) = 0.5*(vn(i+1,m) + vn(i  ,m)) - 1.0/8.0*d2_li(i,m) -  aa_scs*d3_li(i,m)

            end do
        end do
    end do
        
    do m=nst,ned
        do i=st,ed
            
               vl(i,m) = sigma(i)*ve_lil(i,m) + vln(i,m)*(1.0-sigma(i))
                         
               vr(i,m) = sigma(i)*ve_lir(i,m) + vrn(i,m)*(1.0-sigma(i))            

        end do
    end do  

    do i=stc,edc
        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vl(i,m) < chklim(m) ) then
                    vl(i,m) = vn(i  ,m)
                    ndec = ndec + 1
                end if

                if (vr(i,m) < chklim(m) ) then
                    vr(i,m) = vn(i+1,m)
                    ndec = ndec + 1
                end if
           end if
        end do
    end do

end subroutine scsh3ci

subroutine scsh3cicc(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_bound_scheme
    implicit none
    integer(kind_int), intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  , intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(in)  :: sn(1-ngn:ni+ngn,1:3)
    integer(kind_int), intent(in)  :: nfs,nfe
    integer(kind_int), intent(in)  :: nge
    real(kind_real)  , intent(out) :: vl(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(out) :: vr(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , external    :: sub_limit
    integer(kind_int), intent(in)  :: chkid(nst:ned)
    real(kind_real),   intent(in)  :: chklim(nst:ned)
    integer(kind_int), intent(out) :: ndec  
    real(kind_real)   :: vln(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: vrn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: ve_lil(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: ve_lir(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: d2_li(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: d3_li(1-ngn:ni+ngn,nst:ned)  
    real(kind_real)   :: nt,nx,ny,nz,fp(nst:ned)
    real(kind_real)   :: vp(nst:ned),vc(nst:ned)
    real(kind_real)   :: q(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: q1(1-ngn:ni+ngn,nst:ned)    
    real(kind_real)   :: s(1:4,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: g(1:4,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: wl(1:4,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: wr(1:4,1-ngn:ni+ngn,nst:ned)    
    real(kind_real)   :: aa,aa_scs
    real(kind_real)   :: sigma(1-ngn:ni+ngn),imax,imin,indicator
    real(kind_real), parameter :: cl(4) = (/1.D0/16.D0,10.D0/16.D0,5.D0/16.D0,0.D0/16.D0/)
    real(kind_real), parameter :: cr(4) = (/0.D0/16.D0,5.D0/16.D0,10.D0/16.D0,1.D0/16.D0/)      
    real(kind_real), parameter :: eps        = 1.0e-30
    real(kind_real), parameter :: indicatorc = 0.1
    real(kind_real)            :: is,eis2(4),bs,tau5,bl(4),br(4)
    integer(kind_int)          :: i,m,n,j,st,ed,stc,edc
    integer(kind_int)          :: id,YN,ind(1:ni+2*ngn),ind_st,ind_ed,pre_ed
    integer(kind_int)          :: nsec,sec_st(1:ni+2*ngn),sec_ed(1:ni+2*ngn)
    integer(kind_int)          :: nsecn,sec_stn(1:ni+2*ngn),sec_edn(1:ni+2*ngn)
    
    if (ni < 8) then
        call muscl2pvcc(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)
        return
    end if
    
    ndec = 0
    
    aa     = 0.01
    aa_scs = 0.001
    
    wl(:,:,:) = 0.3333
    wr(:,:,:) = 0.3333
    vln(:,:)    = 0.0
    ve_lil(:,:) = 0.0
    vrn(:,:)    = 0.0
    ve_lir(:,:) = 0.0    
           
    if (nfs == nbc_bound_scheme) then
        
        st  = 5
        stc = 1

        do m=nst,ned
            vl(1,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                  ! third-order
            vl(2,m) = (  3.0*vn(1,m) +  6.0*vn(2,m) -      vn(3,m) ) / 8.0   ! third-order 
            vl(3,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 + &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )
            vr(1,m) = vl(1,m)
            vr(2,m) = vl(2,m)
            vr(3,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 - &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )
            vl(4,m) =   ( 150.0*( vn(4,m) + vn(3,m) )          -  &
                           25.0*( vn(5,m) + vn(2,m) )          +  &
                            3.0*( vn(6,m) + vn(1,m) ) ) /256.0 -  &
                      aa*( 10.0*( vn(4,m) - vn(3,m) )          -  &
                            5.0*( vn(5,m) - vn(2,m) )          +  &
                                ( vn(6,m) - vn(1,m) ) ) 
            
            vr(4,m) =   ( 150.0*( vn(4,m) + vn(3,m) )          -  &
                           25.0*( vn(5,m) + vn(2,m) )          +  &
                            3.0*( vn(6,m) + vn(1,m) ) ) /256.0 +  &
                      aa*( 10.0*( vn(4,m) - vn(3,m) )          -  &
                            5.0*( vn(5,m) - vn(2,m) )          +  &
                                ( vn(6,m) - vn(1,m) ) )          
        end do
        
    else
        
        st  = -1
        stc = st
     
    end if
     
    if (nfe == nbc_bound_scheme) then
        
        ed  = ni - 3
        edc = ni + 1

        do m=nst,ned
            vr(ni+1,m) = ( 15.0*vn(ni,m) - 10.0*vn(ni-1,m) +  3.0*vn(ni-2,m) ) / 8.0                ! third-order
            vl(ni  ,m) = (  3.0*vn(ni,m) +  6.0*vn(ni-1,m) -      vn(ni-2,m) ) / 8.0   ! third-order
            vl(ni-1,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 + &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) )
            vl(ni+1,m) = vr(ni+1,m)
            vr(ni  ,m) = vl(ni  ,m)
            vr(ni-1,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 - &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) ) 
            vl(ni-2,m) =   ( 150.0*( vn(ni-2,m) + vn(ni-3,m) )          -  &
                              25.0*( vn(ni-1,m) + vn(ni-4,m) )          +  &
                               3.0*( vn(ni  ,m) + vn(ni-5,m) ) ) /256.0 -  &
                         aa*( 10.0*( vn(ni-2,m) - vn(ni-3,m) )          -  &
                               5.0*( vn(ni-1,m) - vn(ni-4,m) )          +  &
                                   ( vn(ni  ,m) - vn(ni-5,m) ) ) 
            
            vr(ni-2,m) =   ( 150.0*( vn(ni-2,m) + vn(ni-3,m) )          -  &
                              25.0*( vn(ni-1,m) + vn(ni-4,m) )          +  &
                               3.0*( vn(ni  ,m) + vn(ni-5,m) ) ) /256.0 +  &
                         aa*( 10.0*( vn(ni-2,m) - vn(ni-3,m) )          -  &
                               5.0*( vn(ni-1,m) - vn(ni-4,m) )          +  &
                                   ( vn(ni  ,m) - vn(ni-5,m) ) )              
        end do
        
    else
        
        ed  = ni + 3
        edc = ed
            
    end if
    
    do m=nst,ned
        do i=st-3,ed+1
            g(1,i,m) = vn(i  ,m) - vn(i-1,m)            
            g(2,i,m) = vn(i+1,m) - vn(i  ,m)
            g(3,i,m) = vn(i+2,m) - vn(i+1,m)            
        end do
    end do   
    
    do m=nst,ned
        do i=st-3,ed+1
            do n=1,3
                is = g(n,i,m)*g(n,i,m)
                eis2(n) = is**2
            end do

            tau5 = abs(eis2(2)-eis2(1))
            bl(1) = 1.0+tau5/(eis2(1)+eps)
            bl(2) = 1.0+tau5/(eis2(2)+eps)
            bs = bl(1) + bl(2)
            wl(2,i,m) = bl(2)/bs
            
            tau5 = abs(eis2(3)-eis2(2))
            br(2) = 1.0+tau5/(eis2(2)+eps)
            br(3) = 1.0+tau5/(eis2(3)+eps)
            bs = br(2) + br(3)
            wr(2,i,m) = br(2)/bs
        end do
    end do
    
    do i=st-1,ed-1
        imin=1.0
        imax=0.0
        do m=nst,ned
            imin         = min(wl(2,i-3,m),wl(2,i-2,m),wl(2,i-1,m),wl(2,i,m),wl(2,i+1,m),wl(2,i+2,m),wl(2,i+3,m), &
                               wr(2,i-3,m),wr(2,i-2,m),wr(2,i-1,m),wr(2,i,m),wr(2,i+1,m),wr(2,i+2,m),wr(2,i+3,m),imin)                                
            imax         = max(wl(2,i-3,m),wl(2,i-2,m),wl(2,i-1,m),wl(2,i,m),wl(2,i+1,m),wl(2,i+2,m),wl(2,i+3,m), &
                               wr(2,i-3,m),wr(2,i-2,m),wr(2,i-1,m),wr(2,i,m),wr(2,i+1,m),wr(2,i+2,m),wr(2,i+3,m),imax)
        end do
        indicator    = imin/imax
        if(indicator>indicatorc) then
            sigma(i)=1.0
        else
            sigma(i)=0.0
        endif        
    end do
    
    id=0
    do i=st-1,ed-1
        YN=0
        if(sigma(i)<0.5) then
            YN=YN+1
        end if
        if(YN>0) then
            id=id+1
            ind(id)=i           
        end if
    end do
       
    do i=1-ngn,ni+ngn
        sigma(i)=0.0
    end do
    if(ed-st<9) then
        ind_st=st-1
        ind_ed=ed-1
        nsec=0
        
        nsecn=1
        sec_stn(1)=st-1
        sec_edn(1)=ed-1
    else
        if(id==0) then
            nsec=1
            sec_st(1)=st-1
            sec_ed(1)=ed-1
            do n=sec_st(1),sec_ed(1)
                sigma(n)=1.0
            end do
            ind_st=st-1
            ind_ed=ed-1
            nsecn=0
        else
            ind_st=st-1
            ind_ed=ed-1
            nsec=0
            nsecn=0
            if(ind(1)-ind_st>9) then
                nsec=1
                sec_st(1)=ind_st
                sec_ed(1)=ind(1)-1
                pre_ed   =sec_ed(1)
                do n=sec_st(1),sec_ed(1)
                    sigma(n)=1.0
                end do
            else
                nsecn=1
                sec_stn(1)=ind_st
                sec_edn(1)=ind(1)
                pre_ed    =sec_edn(1)
            end if
            do i=2,id
                if(ind(i)-ind(i-1)>9) then
                    nsec=nsec+1
                    sec_st(nsec)=ind(i-1)+1
                    sec_ed(nsec)=ind(i)-1
                    do n=sec_st(nsec),sec_ed(nsec)
                        sigma(n)=1.0
                    end do
                    nsecn=nsecn+1
                    sec_stn(nsecn)=pre_ed+1
                    sec_edn(nsecn)=ind(i-1)
                    pre_ed        =sec_ed(nsec)
                else
                    if(ind(i)-ind(i-1)>1 .or. i==id) then
                        nsecn=nsecn+1
                        sec_stn(nsecn)=pre_ed+1
                        sec_edn(nsecn)=ind(i)
                        pre_ed        =sec_edn(nsecn)
                    end if
                end if
            end do
            if(ind_ed-ind(id)>9) then
                nsec=nsec+1
                sec_st(nsec)=pre_ed+1
                sec_ed(nsec)=ind_ed
                do n=sec_st(nsec),sec_ed(nsec)
                    sigma(n)=1.0
                end do
            else
                nsecn=nsecn+1
                sec_stn(nsecn)=pre_ed+1
                sec_edn(nsecn)=ind_ed
            end if
        end if
    end if
    
    do m=nst,ned
        do j=1,nsecn
            do i=sec_stn(j),sec_edn(j)
                g(1,i,m) =      vn(i-2,m) - 3.0*vn(i-1,m) + 2.0*vn(i  ,m) 
                g(2,i,m) =     -vn(i  ,m) +     vn(i+1,m)
                g(3,i,m) =     -vn(i  ,m) +     vn(i+1,m)
                g(4,i,m) = -2.0*vn(i+1,m) + 3.0*vn(i+2,m) -     vn(i+3,m)

                s(1,i,m) = vn(i-2,m) - 2.0*vn(i-1,m) + vn(i,  m)
                s(2,i,m) = vn(i-1,m) - 2.0*vn(i  ,m) + vn(i+1,m) 
                s(3,i,m) = vn(i  ,m) - 2.0*vn(i+1,m) + vn(i+2,m)
                s(4,i,m) = vn(i+1,m) - 2.0*vn(i+2,m) + vn(i+3,m)
            end do
        end do      
    end do
    
    do j=1,nsecn
        do i=sec_stn(j),sec_edn(j)
            nx = (sn(i,1) + sn(i+1,1))/2.0
            ny = (sn(i,2) + sn(i+1,2))/2.0
            nz = (sn(i,3) + sn(i+1,3))/2.0

            do m=nst,ned
                vp(m) = (vn(i,m) + vn(i+1,m))/2.0
                fp(m) = vn(i,m)
            end do
            call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
            do m=nst,ned
                q(i,m) = vc(m)
            end do
            
            do m=nst,ned
                fp(m) = vn(i+1,m)
            end do
            call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
            do m=nst,ned
                q1(i+1,m) = vc(m)
            end do        
             
            do n=1,4
                do m=nst,ned
                    fp(m) = g(n,i,m)
                end do
                call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
                do m=nst,ned
                    g(n,i,m) = vc(m)
                end do

                do m=nst,ned
                    fp(m) = s(n,i,m)
                end do
                call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
                do m=nst,ned
                    s(n,i,m) = vc(m)
                end do          
            end do
        end do
    end do   
            
            
    do m=nst,ned
        do j=1,nsecn
            do i=sec_stn(j),sec_edn(j)
                do n=1,4
                    is = g(n,i,m)*g(n,i,m) + s(n,i,m)*s(n,i,m)
                    eis2(n) = is**2
                end do
                
                tau5 = abs(eis2(3)-eis2(1))
                
                do n=1,3
                    bl(n) = cl(n)*(1.0+tau5/(eis2(n)+eps))
                end do
                
                bs = bl(1) + bl(2) + bl(3)
                do n=1,3
                    wl(n,i,m) = bl(n)/bs
                end do
          
                tau5 = abs(eis2(4)-eis2(2))           
                do n=2,4
                    br(n) = cr(n)*(1.0+tau5/(eis2(n)+eps))
                end do

                bs = br(2) + br(3) + br(4)
                do n=2,4
                    wr(n,i,m) = br(n)/bs
                end do
            end do
        end do
    end do
    
!!    open(2,file='sigma.dat',form='formatted')
!!    write(2,*) 'variables="i" "sigmal" "sigmar"'
!!    do i=st,ed
!!
!!        write(2,200) i,sigmal(i,1),sigmar(i,1)
!!
!!    end do
!!
!!200 format(I4,1x,e12.5,1x,e12.5)
!!    close(2)
!!
!!    stop    
    
    do m=nst,ned
        do j=1,nsecn
            do i=sec_stn(j),sec_edn(j)
                vln(i,m) =  q(i  ,m) + ( wl(1,i,m)*( 4.0*g(1,i,m)-s(1,i,m)) + &
                                         wl(2,i,m)*( 4.0*g(2,i,m)-s(2,i,m)) + &
                                         wl(3,i,m)*( 4.0*g(3,i,m)-s(3,i,m)) )/8.0

                vrn(i,m) = q1(i+1,m) + ( wr(2,i,m)*(-4.0*g(2,i,m)-s(2,i,m)) + &
                                         wr(3,i,m)*(-4.0*g(3,i,m)-s(3,i,m)) + &
                                         wr(4,i,m)*(-4.0*g(4,i,m)-s(4,i,m)) )/8.0
            end do
        end do
    end do
    
    do j=1,nsecn
        do i=sec_stn(j),sec_edn(j)
            nx = (sn(i,1) + sn(i+1,1))/2.0
            ny = (sn(i,2) + sn(i+1,2))/2.0
            nz = (sn(i,3) + sn(i+1,3))/2.0

            do m=nst,ned
                vp(m) = (vn(i,m) + vn(i+1,m))/2.0
                vc(m) = vln(i,m)
            end do

            call RxQ_p(nst,ned,vp,nt,nx,ny,nz,vc,fp) 
            do m=nst,ned
                vln(i,m) = fp(m)
            end do
            
            do m=nst,ned
                vc(m) = vrn(i ,m)
            end do

            call RxQ_p(nst,ned,vp,nt,nx,ny,nz,vc,fp) 
            do m=nst,ned
                vrn(i,m) = fp(m)
            end do
        end do
    end do     
      
    do j=1,nsec
        call d3_via_T2N(ni,nst,ned,ngn,vn,nfs,nfe,nge,d3_li,sec_st(j),sec_ed(j))
        
        call vec2_via_T2N(ni,nst,ned,ngn,vn,nfs,nfe,nge,d2_li,sec_st(j),sec_ed(j))
    end do
    
    do m=nst,ned
        do j=1,nsec
            do i=sec_st(j),sec_ed(j)
                
                ve_lil(i,m) = 0.5*(vn(i+1,m) + vn(i  ,m)) - 1.0/8.0*d2_li(i,m) +  aa_scs*d3_li(i,m)
                
                ve_lir(i,m) = 0.5*(vn(i+1,m) + vn(i  ,m)) - 1.0/8.0*d2_li(i,m) -  aa_scs*d3_li(i,m)

            end do
        end do
    end do
        
    do m=nst,ned
        do i=st-1,ed-1
            
               vl(i+1,m) = sigma(i)*ve_lil(i,m) + vln(i,m)*(1.0-sigma(i))
                         
               vr(i+1,m) = sigma(i)*ve_lir(i,m) + vrn(i,m)*(1.0-sigma(i))            

        end do
    end do  

    do i=stc,edc
        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vl(i,m) < chklim(m) ) then
                    vl(i,m) = vn(i  ,m)
                    ndec = ndec + 1
                end if

                if (vr(i,m) < chklim(m) ) then
                    vr(i,m) = vn(i+1,m)
                    ndec = ndec + 1
                end if
           end if
        end do
    end do

end subroutine scsh3cicc

subroutine dcsh5ci(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_bound_scheme
    implicit none
    integer(kind_int), intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  , intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(in)  :: sn(1-ngn:ni+ngn,1:3)
    integer(kind_int), intent(in)  :: nfs,nfe
    integer(kind_int), intent(in)  :: nge
    real(kind_real)  , intent(out) :: vl(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(out) :: vr(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , external    :: sub_limit
    integer(kind_int), intent(in)  :: chkid(nst:ned)
    real(kind_real),   intent(in)  :: chklim(nst:ned)
    integer(kind_int), intent(out) :: ndec  
    real(kind_real)   :: vln(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: vrn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: ve_lil(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: ve_lir(1-ngn:ni+ngn,nst:ned) 
    real(kind_real)   :: nt,nx,ny,nz,fp(nst:ned)
    real(kind_real)   :: vp(nst:ned),vc(nst:ned)
    real(kind_real)   :: q(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: q1(1-ngn:ni+ngn,nst:ned)    
    real(kind_real)   :: s(1:4,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: g(1:4,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: wl(1:4,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: wr(1:4,1-ngn:ni+ngn,nst:ned)    
    real(kind_real)   :: aa,aa_scs
    real(kind_real)   :: sigma(1-ngn:ni+ngn),imax,imin,indicator
    real(kind_real), parameter :: cl(4) = (/1.D0/16.D0,10.D0/16.D0,5.D0/16.D0,0.D0/16.D0/)
    real(kind_real), parameter :: cr(4) = (/0.D0/16.D0,5.D0/16.D0,10.D0/16.D0,1.D0/16.D0/)  
    real(kind_real), parameter :: eps        = 1.0e-30
    real(kind_real), parameter :: indicatorc = 0.1
    real(kind_real)            :: is,eis2(4),bs,tau5,bl(4),br(4)
    integer(kind_int)          :: i,m,n,j,st,ed,stc,edc
    integer(kind_int)          :: id,YN,ind(1:ni+2*ngn),ind_st,ind_ed,pre_ed
    integer(kind_int)          :: nsec,sec_st(1:ni+2*ngn),sec_ed(1:ni+2*ngn)
    integer(kind_int)          :: nsecn,sec_stn(1:ni+2*ngn),sec_edn(1:ni+2*ngn)
    
    if (ni < 8) then
        call muscl2pv(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)
        return
    end if
    
    ndec = 0
    
    aa     = 0.01
    aa_scs = 0.001
    
    wl(:,:,:) = 0.3333
    wr(:,:,:) = 0.3333
    vln(:,:)    = 0.0
    ve_lil(:,:) = 0.0
    vrn(:,:)    = 0.0
    ve_lir(:,:) = 0.0    
           
    if (nfs == nbc_bound_scheme) then
        
        st  = 4
        stc = 0

        do m=nst,ned
            vl(0,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                  ! third-order
            vl(1,m) = (  3.0*vn(1,m) +  6.0*vn(2,m) -      vn(3,m) ) / 8.0   ! third-order 
            vl(2,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 + &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )
            vr(0,m) = vl(0,m)
            vr(1,m) = vl(1,m)
            vr(2,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 - &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )
            vl(3,m) =   ( 150.0*( vn(4,m) + vn(3,m) )          -  &
                           25.0*( vn(5,m) + vn(2,m) )          +  &
                            3.0*( vn(6,m) + vn(1,m) ) ) /256.0 -  &
                      aa*( 10.0*( vn(4,m) - vn(3,m) )          -  &
                            5.0*( vn(5,m) - vn(2,m) )          +  &
                                ( vn(6,m) - vn(1,m) ) ) 
            
            vr(3,m) =   ( 150.0*( vn(4,m) + vn(3,m) )          -  &
                           25.0*( vn(5,m) + vn(2,m) )          +  &
                            3.0*( vn(6,m) + vn(1,m) ) ) /256.0 +  &
                      aa*( 10.0*( vn(4,m) - vn(3,m) )          -  &
                            5.0*( vn(5,m) - vn(2,m) )          +  &
                                ( vn(6,m) - vn(1,m) ) )          
        end do
        
    else
        
        st  = -2
        stc = st
     
    end if
     
    if (nfe == nbc_bound_scheme) then
        
        ed  = ni - 4
        edc = ni

        do m=nst,ned
            vr(ni  ,m) = ( 15.0*vn(ni,m) - 10.0*vn(ni-1,m) +  3.0*vn(ni-2,m) ) / 8.0                ! third-order
            vl(ni-1,m) = (  3.0*vn(ni,m) +  6.0*vn(ni-1,m) -      vn(ni-2,m) ) / 8.0   ! third-order
            vl(ni-2,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 + &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) )
            vl(ni  ,m) = vr(ni  ,m)
            vr(ni-1,m) = vl(ni-1,m)
            vr(ni-2,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 - &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) ) 
            vl(ni-3,m) =   ( 150.0*( vn(ni-2,m) + vn(ni-3,m) )          -  &
                              25.0*( vn(ni-1,m) + vn(ni-4,m) )          +  &
                               3.0*( vn(ni  ,m) + vn(ni-5,m) ) ) /256.0 -  &
                         aa*( 10.0*( vn(ni-2,m) - vn(ni-3,m) )          -  &
                               5.0*( vn(ni-1,m) - vn(ni-4,m) )          +  &
                                   ( vn(ni  ,m) - vn(ni-5,m) ) ) 
            
            vr(ni-3,m) =   ( 150.0*( vn(ni-2,m) + vn(ni-3,m) )          -  &
                              25.0*( vn(ni-1,m) + vn(ni-4,m) )          +  &
                               3.0*( vn(ni  ,m) + vn(ni-5,m) ) ) /256.0 +  &
                         aa*( 10.0*( vn(ni-2,m) - vn(ni-3,m) )          -  &
                               5.0*( vn(ni-1,m) - vn(ni-4,m) )          +  &
                                   ( vn(ni  ,m) - vn(ni-5,m) ) )              
        end do
        
    else
        
        ed  = ni + 2
        edc = ed
            
    end if
    
    do m=nst,ned
        do i=st-2,ed+2
            g(1,i,m) = vn(i  ,m) - vn(i-1,m)            
            g(2,i,m) = vn(i+1,m) - vn(i  ,m)
            g(3,i,m) = vn(i+2,m) - vn(i+1,m)            
        end do
    end do   
    
    do m=nst,ned
        do i=st-2,ed+2
            do n=1,3
                is = g(n,i,m)*g(n,i,m)
                eis2(n) = is**2
            end do

            tau5 = abs(eis2(2)-eis2(1))
            bl(1) = 1.0+tau5/(eis2(1)+eps)
            bl(2) = 1.0+tau5/(eis2(2)+eps)
            bs = bl(1) + bl(2)
            wl(2,i,m) = bl(2)/bs
            
            tau5 = abs(eis2(3)-eis2(2))
            br(2) = 1.0+tau5/(eis2(2)+eps)
            br(3) = 1.0+tau5/(eis2(3)+eps)
            bs = br(2) + br(3)
            wr(2,i,m) = br(2)/bs
        end do
    end do
    
    do i=st,ed
        imin=1.0
        imax=0.0
        do m=nst,ned
            imin         = min(wl(2,i-3,m),wl(2,i-2,m),wl(2,i-1,m),wl(2,i,m),wl(2,i+1,m),wl(2,i+2,m),wl(2,i+3,m), &
                               wr(2,i-3,m),wr(2,i-2,m),wr(2,i-1,m),wr(2,i,m),wr(2,i+1,m),wr(2,i+2,m),wr(2,i+3,m),imin)                                
            imax         = max(wl(2,i-3,m),wl(2,i-2,m),wl(2,i-1,m),wl(2,i,m),wl(2,i+1,m),wl(2,i+2,m),wl(2,i+3,m), &
                               wr(2,i-3,m),wr(2,i-2,m),wr(2,i-1,m),wr(2,i,m),wr(2,i+1,m),wr(2,i+2,m),wr(2,i+3,m),imax)
        end do
        indicator    = imin/imax
        if(indicator>indicatorc) then
            sigma(i)=1.0
        else
            sigma(i)=0.0
        endif        
    end do
    
    id=0
    do i=st,ed
        YN=0
        if(sigma(i)<0.5) then
            YN=YN+1
        end if
        if(YN>0) then
            id=id+1
            ind(id)=i           
        end if
    end do
       
    do i=1-ngn,ni+ngn
        sigma(i)=0.0
    end do
    if(ed-st<9) then
        ind_st=st
        ind_ed=ed
        nsec=0
        
        nsecn=1
        sec_stn(1)=st
        sec_edn(1)=ed
    else
        if(id==0) then
            nsec=1
            sec_st(1)=st
            sec_ed(1)=ed
            do n=sec_st(1),sec_ed(1)
                sigma(n)=1.0
            end do
            ind_st=st
            ind_ed=ed
            nsecn=0
        else
            ind_st=st
            ind_ed=ed
            nsec=0
            nsecn=0
            if(ind(1)-ind_st>9) then
                nsec=1
                sec_st(1)=ind_st
                sec_ed(1)=ind(1)-1
                pre_ed   =sec_ed(1)
                do n=sec_st(1),sec_ed(1)
                    sigma(n)=1.0
                end do
            else
                nsecn=1
                sec_stn(1)=ind_st
                sec_edn(1)=ind(1)
                pre_ed    =sec_edn(1)
            end if
            do i=2,id
                if(ind(i)-ind(i-1)>9) then
                    nsec=nsec+1
                    sec_st(nsec)=ind(i-1)+1
                    sec_ed(nsec)=ind(i)-1
                    do n=sec_st(nsec),sec_ed(nsec)
                        sigma(n)=1.0
                    end do
                    nsecn=nsecn+1
                    sec_stn(nsecn)=pre_ed+1
                    sec_edn(nsecn)=ind(i-1)
                    pre_ed        =sec_ed(nsec)
                else
                    if(ind(i)-ind(i-1)>1 .or. i==id) then
                        nsecn=nsecn+1
                        sec_stn(nsecn)=pre_ed+1
                        sec_edn(nsecn)=ind(i)
                        pre_ed        =sec_edn(nsecn)
                    end if
                end if
            end do
            if(ind_ed-ind(id)>9) then
                nsec=nsec+1
                sec_st(nsec)=pre_ed+1
                sec_ed(nsec)=ind_ed
                do n=sec_st(nsec),sec_ed(nsec)
                    sigma(n)=1.0
                end do
            else
                nsecn=nsecn+1
                sec_stn(nsecn)=pre_ed+1
                sec_edn(nsecn)=ind_ed
            end if
        end if
    end if
    
    do m=nst,ned
        do j=1,nsecn
            do i=sec_stn(j),sec_edn(j)
                g(1,i,m) =      vn(i-2,m) - 3.0*vn(i-1,m) + 2.0*vn(i  ,m) 
                g(2,i,m) =     -vn(i  ,m) +     vn(i+1,m)
                g(3,i,m) =     -vn(i  ,m) +     vn(i+1,m)
                g(4,i,m) = -2.0*vn(i+1,m) + 3.0*vn(i+2,m) -     vn(i+3,m)

                s(1,i,m) = vn(i-2,m) - 2.0*vn(i-1,m) + vn(i,  m)
                s(2,i,m) = vn(i-1,m) - 2.0*vn(i  ,m) + vn(i+1,m) 
                s(3,i,m) = vn(i  ,m) - 2.0*vn(i+1,m) + vn(i+2,m)
                s(4,i,m) = vn(i+1,m) - 2.0*vn(i+2,m) + vn(i+3,m)
            end do
        end do      
    end do
    
    do j=1,nsecn
        do i=sec_stn(j),sec_edn(j)
            nx = (sn(i,1) + sn(i+1,1))/2.0
            ny = (sn(i,2) + sn(i+1,2))/2.0
            nz = (sn(i,3) + sn(i+1,3))/2.0

            do m=nst,ned
                vp(m) = (vn(i,m) + vn(i+1,m))/2.0
                fp(m) = vn(i,m)
            end do
            call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
            do m=nst,ned
                q(i,m) = vc(m)
            end do
            
            do m=nst,ned
                fp(m) = vn(i+1,m)
            end do
            call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
            do m=nst,ned
                q1(i+1,m) = vc(m)
            end do        
             
            do n=1,4
                do m=nst,ned
                    fp(m) = g(n,i,m)
                end do
                call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
                do m=nst,ned
                    g(n,i,m) = vc(m)
                end do

                do m=nst,ned
                    fp(m) = s(n,i,m)
                end do
                call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
                do m=nst,ned
                    s(n,i,m) = vc(m)
                end do          
            end do
        end do
    end do   
            
            
    do m=nst,ned
        do j=1,nsecn
            do i=sec_stn(j),sec_edn(j)
                do n=1,4
                    is = g(n,i,m)*g(n,i,m) + s(n,i,m)*s(n,i,m)
                    eis2(n) = is**2
                end do
                
                tau5 = abs(eis2(3)-eis2(1))
                
                do n=1,3
                    bl(n) = cl(n)*(1.0+tau5/(eis2(n)+eps))
                end do
                
                bs = bl(1) + bl(2) + bl(3)
                do n=1,3
                    wl(n,i,m) = bl(n)/bs
                end do
          
                tau5 = abs(eis2(4)-eis2(2))           
                do n=2,4
                    br(n) = cr(n)*(1.0+tau5/(eis2(n)+eps))
                end do

                bs = br(2) + br(3) + br(4)
                do n=2,4
                    wr(n,i,m) = br(n)/bs
                end do
            end do
        end do
    end do
    
!!    open(2,file='sigma.dat',form='formatted')
!!    write(2,*) 'variables="i" "sigmal" "sigmar"'
!!    do i=st,ed
!!
!!        write(2,200) i,sigmal(i,1),sigmar(i,1)
!!
!!    end do
!!
!!200 format(I4,1x,e12.5,1x,e12.5)
!!    close(2)
!!
!!    stop    
    
    do m=nst,ned
        do j=1,nsecn
            do i=sec_stn(j),sec_edn(j)
                vln(i,m) =  q(i  ,m) + ( wl(1,i,m)*( 4.0*g(1,i,m)-s(1,i,m)) + &
                                         wl(2,i,m)*( 4.0*g(2,i,m)-s(2,i,m)) + &
                                         wl(3,i,m)*( 4.0*g(3,i,m)-s(3,i,m)) )/8.0

                vrn(i,m) = q1(i+1,m) + ( wr(2,i,m)*(-4.0*g(2,i,m)-s(2,i,m)) + &
                                         wr(3,i,m)*(-4.0*g(3,i,m)-s(3,i,m)) + &
                                         wr(4,i,m)*(-4.0*g(4,i,m)-s(4,i,m)) )/8.0
            end do
        end do
    end do
    
    do j=1,nsecn
        do i=sec_stn(j),sec_edn(j)
            nx = (sn(i,1) + sn(i+1,1))/2.0
            ny = (sn(i,2) + sn(i+1,2))/2.0
            nz = (sn(i,3) + sn(i+1,3))/2.0

            do m=nst,ned
                vp(m) = (vn(i,m) + vn(i+1,m))/2.0
                vc(m) = vln(i,m)
            end do

            call RxQ_p(nst,ned,vp,nt,nx,ny,nz,vc,fp) 
            do m=nst,ned
                vln(i,m) = fp(m)
            end do
            
            do m=nst,ned
                vc(m) = vrn(i ,m)
            end do

            call RxQ_p(nst,ned,vp,nt,nx,ny,nz,vc,fp) 
            do m=nst,ned
                vrn(i,m) = fp(m)
            end do
        end do
    end do     
      
    do j=1,nsec
        
        call hdcs5ciN(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,ve_lil,ve_lir,sub_limit,chkid,chklim,ndec,sec_st(j),sec_ed(j))        
        
    end do
        
    do m=nst,ned
        do i=st,ed
            
               vl(i,m) = sigma(i)*ve_lil(i,m) + vln(i,m)*(1.0-sigma(i))
                         
               vr(i,m) = sigma(i)*ve_lir(i,m) + vrn(i,m)*(1.0-sigma(i))            

        end do
    end do  

    do i=stc,edc
        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vl(i,m) < chklim(m) ) then
                    vl(i,m) = vn(i  ,m)
                    ndec = ndec + 1
                end if

                if (vr(i,m) < chklim(m) ) then
                    vr(i,m) = vn(i+1,m)
                    ndec = ndec + 1
                end if
           end if
        end do
    end do

end subroutine dcsh5ci

subroutine dcsh5cicc(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_bound_scheme
    implicit none
    integer(kind_int), intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  , intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(in)  :: sn(1-ngn:ni+ngn,1:3)
    integer(kind_int), intent(in)  :: nfs,nfe
    integer(kind_int), intent(in)  :: nge
    real(kind_real)  , intent(out) :: vl(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(out) :: vr(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , external    :: sub_limit
    integer(kind_int), intent(in)  :: chkid(nst:ned)
    real(kind_real),   intent(in)  :: chklim(nst:ned)
    integer(kind_int), intent(out) :: ndec  
    real(kind_real)   :: vln(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: vrn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: ve_lil(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: ve_lir(1-ngn:ni+ngn,nst:ned) 
    real(kind_real)   :: nt,nx,ny,nz,fp(nst:ned)
    real(kind_real)   :: vp(nst:ned),vc(nst:ned)
    real(kind_real)   :: q(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: q1(1-ngn:ni+ngn,nst:ned)    
    real(kind_real)   :: s(1:4,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: g(1:4,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: wl(1:4,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: wr(1:4,1-ngn:ni+ngn,nst:ned)    
    real(kind_real)   :: aa,aa_scs
    real(kind_real)   :: sigma(1-ngn:ni+ngn),imax,imin,indicator
    real(kind_real), parameter :: cl(4) = (/1.D0/16.D0,10.D0/16.D0,5.D0/16.D0,0.D0/16.D0/)
    real(kind_real), parameter :: cr(4) = (/0.D0/16.D0,5.D0/16.D0,10.D0/16.D0,1.D0/16.D0/)  
    real(kind_real), parameter :: eps        = 1.0e-30
    real(kind_real), parameter :: indicatorc = 0.1
    real(kind_real)            :: is,eis2(4),bs,tau5,bl(4),br(4)
    integer(kind_int)          :: i,m,n,j,st,ed,stc,edc
    integer(kind_int)          :: id,YN,ind(1:ni+2*ngn),ind_st,ind_ed,pre_ed
    integer(kind_int)          :: nsec,sec_st(1:ni+2*ngn),sec_ed(1:ni+2*ngn)
    integer(kind_int)          :: nsecn,sec_stn(1:ni+2*ngn),sec_edn(1:ni+2*ngn)
    
    if (ni < 8) then
        call muscl2pvcc(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)
        return
    end if
    
    ndec = 0
    
    aa     = 0.01
    aa_scs = 0.001
    
    wl(:,:,:) = 0.3333
    wr(:,:,:) = 0.3333
    vln(:,:)    = 0.0
    ve_lil(:,:) = 0.0
    vrn(:,:)    = 0.0
    ve_lir(:,:) = 0.0    
           
    if (nfs == nbc_bound_scheme) then
        
        st  = 5
        stc = 1

        do m=nst,ned
            vl(1,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                  ! third-order
            vl(2,m) = (  3.0*vn(1,m) +  6.0*vn(2,m) -      vn(3,m) ) / 8.0   ! third-order 
            vl(3,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 + &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )
            vr(1,m) = vl(1,m)
            vr(2,m) = vl(2,m)
            vr(3,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 - &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )
            vl(4,m) =   ( 150.0*( vn(4,m) + vn(3,m) )          -  &
                           25.0*( vn(5,m) + vn(2,m) )          +  &
                            3.0*( vn(6,m) + vn(1,m) ) ) /256.0 -  &
                      aa*( 10.0*( vn(4,m) - vn(3,m) )          -  &
                            5.0*( vn(5,m) - vn(2,m) )          +  &
                                ( vn(6,m) - vn(1,m) ) ) 
            
            vr(4,m) =   ( 150.0*( vn(4,m) + vn(3,m) )          -  &
                           25.0*( vn(5,m) + vn(2,m) )          +  &
                            3.0*( vn(6,m) + vn(1,m) ) ) /256.0 +  &
                      aa*( 10.0*( vn(4,m) - vn(3,m) )          -  &
                            5.0*( vn(5,m) - vn(2,m) )          +  &
                                ( vn(6,m) - vn(1,m) ) )          
        end do
        
    else
        
        st  = -1
        stc = st
     
    end if
     
    if (nfe == nbc_bound_scheme) then
        
        ed  = ni - 3
        edc = ni + 1

        do m=nst,ned
            vr(ni+1,m) = ( 15.0*vn(ni,m) - 10.0*vn(ni-1,m) +  3.0*vn(ni-2,m) ) / 8.0                ! third-order
            vl(ni  ,m) = (  3.0*vn(ni,m) +  6.0*vn(ni-1,m) -      vn(ni-2,m) ) / 8.0   ! third-order
            vl(ni-1,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 + &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) )
            vl(ni+1,m) = vr(ni+1,m)
            vr(ni  ,m) = vl(ni  ,m)
            vr(ni-1,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 - &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) ) 
            vl(ni-2,m) =   ( 150.0*( vn(ni-2,m) + vn(ni-3,m) )          -  &
                              25.0*( vn(ni-1,m) + vn(ni-4,m) )          +  &
                               3.0*( vn(ni  ,m) + vn(ni-5,m) ) ) /256.0 -  &
                         aa*( 10.0*( vn(ni-2,m) - vn(ni-3,m) )          -  &
                               5.0*( vn(ni-1,m) - vn(ni-4,m) )          +  &
                                   ( vn(ni  ,m) - vn(ni-5,m) ) ) 
            
            vr(ni-2,m) =   ( 150.0*( vn(ni-2,m) + vn(ni-3,m) )          -  &
                              25.0*( vn(ni-1,m) + vn(ni-4,m) )          +  &
                               3.0*( vn(ni  ,m) + vn(ni-5,m) ) ) /256.0 +  &
                         aa*( 10.0*( vn(ni-2,m) - vn(ni-3,m) )          -  &
                               5.0*( vn(ni-1,m) - vn(ni-4,m) )          +  &
                                   ( vn(ni  ,m) - vn(ni-5,m) ) )              
        end do
        
    else
        
        ed  = ni + 3
        edc = ed
            
    end if
    
    do m=nst,ned
        do i=st-3,ed+1
            g(1,i,m) = vn(i  ,m) - vn(i-1,m)            
            g(2,i,m) = vn(i+1,m) - vn(i  ,m)
            g(3,i,m) = vn(i+2,m) - vn(i+1,m)            
        end do
    end do   
    
    do m=nst,ned
        do i=st-3,ed+1
            do n=1,3
                is = g(n,i,m)*g(n,i,m)
                eis2(n) = is**2
            end do

            tau5 = abs(eis2(2)-eis2(1))
            bl(1) = 1.0+tau5/(eis2(1)+eps)
            bl(2) = 1.0+tau5/(eis2(2)+eps)
            bs = bl(1) + bl(2)
            wl(2,i,m) = bl(2)/bs
            
            tau5 = abs(eis2(3)-eis2(2))
            br(2) = 1.0+tau5/(eis2(2)+eps)
            br(3) = 1.0+tau5/(eis2(3)+eps)
            bs = br(2) + br(3)
            wr(2,i,m) = br(2)/bs
        end do
    end do
    
    do i=st-1,ed-1
        imin=1.0
        imax=0.0
        do m=nst,ned
            imin         = min(wl(2,i-3,m),wl(2,i-2,m),wl(2,i-1,m),wl(2,i,m),wl(2,i+1,m),wl(2,i+2,m),wl(2,i+3,m), &
                               wr(2,i-3,m),wr(2,i-2,m),wr(2,i-1,m),wr(2,i,m),wr(2,i+1,m),wr(2,i+2,m),wr(2,i+3,m),imin)                                
            imax         = max(wl(2,i-3,m),wl(2,i-2,m),wl(2,i-1,m),wl(2,i,m),wl(2,i+1,m),wl(2,i+2,m),wl(2,i+3,m), &
                               wr(2,i-3,m),wr(2,i-2,m),wr(2,i-1,m),wr(2,i,m),wr(2,i+1,m),wr(2,i+2,m),wr(2,i+3,m),imax)
        end do
        indicator    = imin/imax
        if(indicator>indicatorc) then
            sigma(i)=1.0
        else
            sigma(i)=0.0
        endif        
    end do
    
    id=0
    do i=st-1,ed-1
        YN=0
        if(sigma(i)<0.5) then
            YN=YN+1
        end if
        if(YN>0) then
            id=id+1
            ind(id)=i           
        end if
    end do
       
    do i=1-ngn,ni+ngn
        sigma(i)=0.0
    end do
    if(ed-st<9) then
        ind_st=st-1
        ind_ed=ed-1
        nsec=0
        
        nsecn=1
        sec_stn(1)=st-1
        sec_edn(1)=ed-1
    else
        if(id==0) then
            nsec=1
            sec_st(1)=st-1
            sec_ed(1)=ed-1
            do n=sec_st(1),sec_ed(1)
                sigma(n)=1.0
            end do
            ind_st=st-1
            ind_ed=ed-1
            nsecn=0
        else
            ind_st=st-1
            ind_ed=ed-1
            nsec=0
            nsecn=0
            if(ind(1)-ind_st>9) then
                nsec=1
                sec_st(1)=ind_st
                sec_ed(1)=ind(1)-1
                pre_ed   =sec_ed(1)
                do n=sec_st(1),sec_ed(1)
                    sigma(n)=1.0
                end do
            else
                nsecn=1
                sec_stn(1)=ind_st
                sec_edn(1)=ind(1)
                pre_ed    =sec_edn(1)
            end if
            do i=2,id
                if(ind(i)-ind(i-1)>9) then
                    nsec=nsec+1
                    sec_st(nsec)=ind(i-1)+1
                    sec_ed(nsec)=ind(i)-1
                    do n=sec_st(nsec),sec_ed(nsec)
                        sigma(n)=1.0
                    end do
                    nsecn=nsecn+1
                    sec_stn(nsecn)=pre_ed+1
                    sec_edn(nsecn)=ind(i-1)
                    pre_ed        =sec_ed(nsec)
                else
                    if(ind(i)-ind(i-1)>1 .or. i==id) then
                        nsecn=nsecn+1
                        sec_stn(nsecn)=pre_ed+1
                        sec_edn(nsecn)=ind(i)
                        pre_ed        =sec_edn(nsecn)
                    end if
                end if
            end do
            if(ind_ed-ind(id)>9) then
                nsec=nsec+1
                sec_st(nsec)=pre_ed+1
                sec_ed(nsec)=ind_ed
                do n=sec_st(nsec),sec_ed(nsec)
                    sigma(n)=1.0
                end do
            else
                nsecn=nsecn+1
                sec_stn(nsecn)=pre_ed+1
                sec_edn(nsecn)=ind_ed
            end if
        end if
    end if
    
    do m=nst,ned
        do j=1,nsecn
            do i=sec_stn(j),sec_edn(j)
                g(1,i,m) =      vn(i-2,m) - 3.0*vn(i-1,m) + 2.0*vn(i  ,m) 
                g(2,i,m) =     -vn(i  ,m) +     vn(i+1,m)
                g(3,i,m) =     -vn(i  ,m) +     vn(i+1,m)
                g(4,i,m) = -2.0*vn(i+1,m) + 3.0*vn(i+2,m) -     vn(i+3,m)

                s(1,i,m) = vn(i-2,m) - 2.0*vn(i-1,m) + vn(i,  m)
                s(2,i,m) = vn(i-1,m) - 2.0*vn(i  ,m) + vn(i+1,m) 
                s(3,i,m) = vn(i  ,m) - 2.0*vn(i+1,m) + vn(i+2,m)
                s(4,i,m) = vn(i+1,m) - 2.0*vn(i+2,m) + vn(i+3,m)
            end do
        end do      
    end do
    
    do j=1,nsecn
        do i=sec_stn(j),sec_edn(j)
            nx = (sn(i,1) + sn(i+1,1))/2.0
            ny = (sn(i,2) + sn(i+1,2))/2.0
            nz = (sn(i,3) + sn(i+1,3))/2.0

            do m=nst,ned
                vp(m) = (vn(i,m) + vn(i+1,m))/2.0
                fp(m) = vn(i,m)
            end do
            call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
            do m=nst,ned
                q(i,m) = vc(m)
            end do
            
            do m=nst,ned
                fp(m) = vn(i+1,m)
            end do
            call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
            do m=nst,ned
                q1(i+1,m) = vc(m)
            end do        
             
            do n=1,4
                do m=nst,ned
                    fp(m) = g(n,i,m)
                end do
                call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
                do m=nst,ned
                    g(n,i,m) = vc(m)
                end do

                do m=nst,ned
                    fp(m) = s(n,i,m)
                end do
                call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
                do m=nst,ned
                    s(n,i,m) = vc(m)
                end do          
            end do
        end do
    end do   
            
            
    do m=nst,ned
        do j=1,nsecn
            do i=sec_stn(j),sec_edn(j)
                do n=1,4
                    is = g(n,i,m)*g(n,i,m) + s(n,i,m)*s(n,i,m)
                    eis2(n) = is**2
                end do
                
                tau5 = abs(eis2(3)-eis2(1))
                
                do n=1,3
                    bl(n) = cl(n)*(1.0+tau5/(eis2(n)+eps))
                end do
                
                bs = bl(1) + bl(2) + bl(3)
                do n=1,3
                    wl(n,i,m) = bl(n)/bs
                end do
          
                tau5 = abs(eis2(4)-eis2(2))           
                do n=2,4
                    br(n) = cr(n)*(1.0+tau5/(eis2(n)+eps))
                end do

                bs = br(2) + br(3) + br(4)
                do n=2,4
                    wr(n,i,m) = br(n)/bs
                end do
            end do
        end do
    end do
    
!!    open(2,file='sigma.dat',form='formatted')
!!    write(2,*) 'variables="i" "sigmal" "sigmar"'
!!    do i=st,ed
!!
!!        write(2,200) i,sigmal(i,1),sigmar(i,1)
!!
!!    end do
!!
!!200 format(I4,1x,e12.5,1x,e12.5)
!!    close(2)
!!
!!    stop    
    
    do m=nst,ned
        do j=1,nsecn
            do i=sec_stn(j),sec_edn(j)
                vln(i,m) =  q(i  ,m) + ( wl(1,i,m)*( 4.0*g(1,i,m)-s(1,i,m)) + &
                                         wl(2,i,m)*( 4.0*g(2,i,m)-s(2,i,m)) + &
                                         wl(3,i,m)*( 4.0*g(3,i,m)-s(3,i,m)) )/8.0

                vrn(i,m) = q1(i+1,m) + ( wr(2,i,m)*(-4.0*g(2,i,m)-s(2,i,m)) + &
                                         wr(3,i,m)*(-4.0*g(3,i,m)-s(3,i,m)) + &
                                         wr(4,i,m)*(-4.0*g(4,i,m)-s(4,i,m)) )/8.0
            end do
        end do
    end do
    
    do j=1,nsecn
        do i=sec_stn(j),sec_edn(j)
            nx = (sn(i,1) + sn(i+1,1))/2.0
            ny = (sn(i,2) + sn(i+1,2))/2.0
            nz = (sn(i,3) + sn(i+1,3))/2.0

            do m=nst,ned
                vp(m) = (vn(i,m) + vn(i+1,m))/2.0
                vc(m) = vln(i,m)
            end do

            call RxQ_p(nst,ned,vp,nt,nx,ny,nz,vc,fp) 
            do m=nst,ned
                vln(i,m) = fp(m)
            end do
            
            do m=nst,ned
                vc(m) = vrn(i ,m)
            end do

            call RxQ_p(nst,ned,vp,nt,nx,ny,nz,vc,fp) 
            do m=nst,ned
                vrn(i,m) = fp(m)
            end do
        end do
    end do     
      
    do j=1,nsec
        
        call hdcs5ciN(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,ve_lil,ve_lir,sub_limit,chkid,chklim,ndec,sec_st(j),sec_ed(j))        
        
    end do
        
    do m=nst,ned
        do i=st-1,ed-1
            
               vl(i+1,m) = sigma(i)*ve_lil(i,m) + vln(i,m)*(1.0-sigma(i))
                         
               vr(i+1,m) = sigma(i)*ve_lir(i,m) + vrn(i,m)*(1.0-sigma(i))            

        end do
    end do  

    do i=stc,edc
        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vl(i,m) < chklim(m) ) then
                    vl(i,m) = vn(i  ,m)
                    ndec = ndec + 1
                end if

                if (vr(i,m) < chklim(m) ) then
                    vr(i,m) = vn(i+1,m)
                    ndec = ndec + 1
                end if
           end if
        end do
    end do

end subroutine dcsh5cicc

subroutine scsn3ci(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_bound_scheme
    implicit none
    integer(kind_int), intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  , intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(in)  :: sn(1-ngn:ni+ngn,1:3)
    integer(kind_int), intent(in)  :: nfs,nfe
    integer(kind_int), intent(in)  :: nge
    real(kind_real)  , intent(out) :: vl(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(out) :: vr(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , external    :: sub_limit
    integer(kind_int), intent(in)  :: chkid(nst:ned)
    real(kind_real),   intent(in)  :: chklim(nst:ned)
    integer(kind_int), intent(out) :: ndec   
    real(kind_real)   :: nt,nx,ny,nz,fp(nst:ned)
    real(kind_real)   :: vp(nst:ned),vc(nst:ned)
    real(kind_real)   :: q(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: q1(1-ngn:ni+ngn,nst:ned)    
    real(kind_real)   :: s(1:4,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: g(1:4,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: wl(1:4,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: wr(1:4,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: dl3(-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: dr3(-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: aa
    real(kind_real), parameter :: cl(4) = (/1.D0/16.D0,10.D0/16.D0,5.D0/16.D0,0.D0/16.D0/)
    real(kind_real), parameter :: cr(4) = (/0.D0/16.D0,5.D0/16.D0,10.D0/16.D0,1.D0/16.D0/)  
    real(kind_real), parameter :: eps        = 1.0e-30
    real(kind_real)            :: is,eis2(4),bs,tau5,bl(4),br(4)    
    integer(kind_int)          :: i,m,n,st,ed,stc,edc
    
    if (ni < 6) then
        call muscl2pv(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)
        return
    end if
    
    ndec = 0
    
    aa     = 0.01
    
    wl(:,:,:) = 0.33333
    wr(:,:,:) = 0.33333  
           
    if (nfs == nbc_bound_scheme) then
        
        st  = 3
        stc = 0

        do m=nst,ned
            vl(0,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                  ! third-order
            vl(1,m) = (  3.0*vn(1,m) +  6.0*vn(2,m) -      vn(3,m) ) / 8.0   ! third-order 

            vl(2,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 + &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )
            vr(0,m) = vl(0,m)
            vr(1,m) = vl(1,m)

            vr(2,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 - &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )         
        end do
        
    else
        
        st  = -1
        stc = st
     
    end if
     
    if (nfe == nbc_bound_scheme) then
        
        ed  = ni - 3
        edc = ni

        do m=nst,ned
            vr(ni  ,m) = ( 15.0*vn(ni,m) - 10.0*vn(ni-1,m) +  3.0*vn(ni-2,m) ) / 8.0                ! third-order
            vl(ni-1,m) = (  3.0*vn(ni,m) +  6.0*vn(ni-1,m) -      vn(ni-2,m) ) / 8.0   ! third-order
            vl(ni-2,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 + &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) )
            vl(ni  ,m) = vr(ni  ,m)
            vr(ni-1,m) = vl(ni-1,m)

            vr(ni-2,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 - &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) )           
        end do
        
    else
        
        ed  = ni + 1
        edc = ed
            
    end if
    
    do m=nst,ned
        do i=st,ed
            g(1,i,m) =      vn(i-2,m) - 3.0*vn(i-1,m) + 2.0*vn(i  ,m) 
            g(2,i,m) =     -vn(i  ,m) +     vn(i+1,m)
            g(3,i,m) =     -vn(i  ,m) +     vn(i+1,m)
            g(4,i,m) = -2.0*vn(i+1,m) + 3.0*vn(i+2,m) -     vn(i+3,m)

            s(1,i,m) = vn(i-2,m) - 2.0*vn(i-1,m) + vn(i,  m)
            s(2,i,m) = vn(i-1,m) - 2.0*vn(i  ,m) + vn(i+1,m) 
            s(3,i,m) = vn(i  ,m) - 2.0*vn(i+1,m) + vn(i+2,m)
            s(4,i,m) = vn(i+1,m) - 2.0*vn(i+2,m) + vn(i+3,m)
        end do       
    end do
    
    do i=st,ed
        nx = (sn(i,1) + sn(i+1,1))/2.0
        ny = (sn(i,2) + sn(i+1,2))/2.0
        nz = (sn(i,3) + sn(i+1,3))/2.0

        do m=nst,ned
            vp(m) = (vn(i,m) + vn(i+1,m))/2.0
            fp(m) = vn(i,m)
        end do
        call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc)  
        do m=nst,ned
            q(i,m) = vc(m)
        end do
        
        do m=nst,ned
            fp(m) = vn(i+1,m)
        end do
        call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
        do m=nst,ned
            q1(i+1,m) = vc(m)
        end do        
         
        do n=1,4
            do m=nst,ned
                fp(m) = g(n,i,m)
            end do
            call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
            do m=nst,ned
                g(n,i,m) = vc(m)
            end do

            do m=nst,ned
                fp(m) = s(n,i,m)
            end do
            call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
            do m=nst,ned
                s(n,i,m) = vc(m)
            end do
        end do
    end do
            
            
    do m=nst,ned
        do i=st,ed
            do n=1,4
                is = g(n,i,m)*g(n,i,m) + s(n,i,m)*s(n,i,m)
                eis2(n) = is**2
            end do
            
            tau5 = abs(eis2(3)-eis2(1))
            
            do n=1,3
                bl(n) = cl(n)*(1.0+tau5/(eis2(n)+eps))
            end do
            
            bs = bl(1) + bl(2) + bl(3)
            do n=1,3
                wl(n,i,m) = bl(n)/bs
            end do
          
            tau5 = abs(eis2(4)-eis2(2))           
            do n=2,4
                br(n) = cr(n)*(1.0+tau5/(eis2(n)+eps))
            end do

            bs = br(2) + br(3) + br(4)
            do n=2,4
                wr(n,i,m) = br(n)/bs
            end do
        end do
    end do
      
    do m=nst,ned
        do i=st,ed
            vl(i,m) =  q(i  ,m) + ( wl(1,i,m)*( 4.0*g(1,i,m)-s(1,i,m)) + &
                                    wl(2,i,m)*( 4.0*g(2,i,m)-s(2,i,m)) + &
                                    wl(3,i,m)*( 4.0*g(3,i,m)-s(3,i,m)) )/8.0
        end do

        do i=st,ed
            vr(i,m) = q1(i+1,m) + ( wr(2,i,m)*(-4.0*g(2,i,m)-s(2,i,m)) + &
                                    wr(3,i,m)*(-4.0*g(3,i,m)-s(3,i,m)) + &
                                    wr(4,i,m)*(-4.0*g(4,i,m)-s(4,i,m)) )/8.0
        end do
    end do
    
    do i=st,ed
        nx = (sn(i,1) + sn(i+1,1))/2.0
        ny = (sn(i,2) + sn(i+1,2))/2.0
        nz = (sn(i,3) + sn(i+1,3))/2.0

        do m=nst,ned
            vp(m) = (vn(i,m) + vn(i+1,m))/2.0
            vc(m) = vl(i,m)
        end do

        call RxQ_p(nst,ned,vp,nt,nx,ny,nz,vc,fp) 
        do m=nst,ned
            vl(i,m) = fp(m)
        end do
        
        do m=nst,ned
            vc(m) = vr(i ,m)
        end do

        call RxQ_p(nst,ned,vp,nt,nx,ny,nz,vc,fp) 
        do m=nst,ned
            vr(i,m) = fp(m)
        end do        
    end do   

    do i=stc,edc
        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vl(i,m) < chklim(m) ) then
                    vl(i,m) = vn(i  ,m)
                    ndec = ndec + 1
                end if

                if (vr(i,m) < chklim(m) ) then
                    vr(i,m) = vn(i+1,m)
                    ndec = ndec + 1
                end if
           end if
        end do
    end do

end subroutine scsn3ci

subroutine scsn3cicc(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_bound_scheme
    implicit none
    integer(kind_int), intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  , intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(in)  :: sn(1-ngn:ni+ngn,1:3)
    integer(kind_int), intent(in)  :: nfs,nfe
    integer(kind_int), intent(in)  :: nge
    real(kind_real)  , intent(out) :: vl(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(out) :: vr(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , external    :: sub_limit
    integer(kind_int), intent(in)  :: chkid(nst:ned)
    real(kind_real),   intent(in)  :: chklim(nst:ned)
    integer(kind_int), intent(out) :: ndec   
    real(kind_real)   :: nt,nx,ny,nz,fp(nst:ned)
    real(kind_real)   :: vp(nst:ned),vc(nst:ned)
    real(kind_real)   :: q(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: q1(1-ngn:ni+ngn,nst:ned)    
    real(kind_real)   :: s(1:4,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: g(1:4,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: wl(1:4,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: wr(1:4,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: dl3(-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: dr3(-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: aa
    real(kind_real), parameter :: cl(4) = (/1.D0/16.D0,10.D0/16.D0,5.D0/16.D0,0.D0/16.D0/)
    real(kind_real), parameter :: cr(4) = (/0.D0/16.D0,5.D0/16.D0,10.D0/16.D0,1.D0/16.D0/)  
    real(kind_real), parameter :: eps        = 1.0e-30
    real(kind_real)            :: is,eis2(4),bs,tau5,bl(4),br(4)    
    integer(kind_int)          :: i,m,n,st,ed,stc,edc
    
    if (ni < 6) then
        call muscl2pvcc(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)
        return
    end if
    
    ndec = 0
    
    aa     = 0.01
    
    wl(:,:,:) = 0.33333
    wr(:,:,:) = 0.33333  
           
    if (nfs == nbc_bound_scheme) then
        
        st  = 4
        stc = 1

        do m=nst,ned
            vl(1,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                  ! third-order
            vl(2,m) = (  3.0*vn(1,m) +  6.0*vn(2,m) -      vn(3,m) ) / 8.0   ! third-order 

            vl(3,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 + &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )
            vr(1,m) = vl(1,m)
            vr(2,m) = vl(2,m)

            vr(3,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 - &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )         
        end do
        
    else
        
        st  = 0
        stc = st
     
    end if
     
    if (nfe == nbc_bound_scheme) then
        
        ed  = ni - 2
        edc = ni + 1

        do m=nst,ned
            vr(ni+1,m) = ( 15.0*vn(ni,m) - 10.0*vn(ni-1,m) +  3.0*vn(ni-2,m) ) / 8.0                ! third-order
            vl(ni  ,m) = (  3.0*vn(ni,m) +  6.0*vn(ni-1,m) -      vn(ni-2,m) ) / 8.0   ! third-order
            vl(ni-1,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 + &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) )
            vl(ni+1,m) = vr(ni+1,m)
            vr(ni  ,m) = vl(ni  ,m)

            vr(ni-1,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 - &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) )           
        end do
        
    else
        
        ed  = ni + 2
        edc = ed
            
    end if
    
    do m=nst,ned
        do i=st-1,ed-1
            g(1,i,m) =      vn(i-2,m) - 3.0*vn(i-1,m) + 2.0*vn(i  ,m) 
            g(2,i,m) =     -vn(i  ,m) +     vn(i+1,m)
            g(3,i,m) =     -vn(i  ,m) +     vn(i+1,m)
            g(4,i,m) = -2.0*vn(i+1,m) + 3.0*vn(i+2,m) -     vn(i+3,m)

            s(1,i,m) = vn(i-2,m) - 2.0*vn(i-1,m) + vn(i,  m)
            s(2,i,m) = vn(i-1,m) - 2.0*vn(i  ,m) + vn(i+1,m) 
            s(3,i,m) = vn(i  ,m) - 2.0*vn(i+1,m) + vn(i+2,m)
            s(4,i,m) = vn(i+1,m) - 2.0*vn(i+2,m) + vn(i+3,m)
        end do       
    end do
    
    do i=st-1,ed-1
        nx = (sn(i,1) + sn(i+1,1))/2.0
        ny = (sn(i,2) + sn(i+1,2))/2.0
        nz = (sn(i,3) + sn(i+1,3))/2.0

        do m=nst,ned
            vp(m) = (vn(i,m) + vn(i+1,m))/2.0
            fp(m) = vn(i,m)
        end do
        call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc)  
        do m=nst,ned
            q(i,m) = vc(m)
        end do
        
        do m=nst,ned
            fp(m) = vn(i+1,m)
        end do
        call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
        do m=nst,ned
            q1(i+1,m) = vc(m)
        end do        
         
        do n=1,4
            do m=nst,ned
                fp(m) = g(n,i,m)
            end do
            call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
            do m=nst,ned
                g(n,i,m) = vc(m)
            end do

            do m=nst,ned
                fp(m) = s(n,i,m)
            end do
            call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
            do m=nst,ned
                s(n,i,m) = vc(m)
            end do
        end do
    end do
            
            
    do m=nst,ned
        do i=st-1,ed-1
            do n=1,4
                is = g(n,i,m)*g(n,i,m) + s(n,i,m)*s(n,i,m)
                eis2(n) = is**2
            end do
            
            tau5 = abs(eis2(3)-eis2(1))
            
            do n=1,3
                bl(n) = cl(n)*(1.0+tau5/(eis2(n)+eps))
            end do
            
            bs = bl(1) + bl(2) + bl(3)
            do n=1,3
                wl(n,i,m) = bl(n)/bs
            end do
          
            tau5 = abs(eis2(4)-eis2(2))           
            do n=2,4
                br(n) = cr(n)*(1.0+tau5/(eis2(n)+eps))
            end do

            bs = br(2) + br(3) + br(4)
            do n=2,4
                wr(n,i,m) = br(n)/bs
            end do
        end do
    end do
      
    do m=nst,ned
        do i=st-1,ed-1
            vl(i+1,m) =  q(i  ,m) + ( wl(1,i,m)*( 4.0*g(1,i,m)-s(1,i,m)) + &
                                      wl(2,i,m)*( 4.0*g(2,i,m)-s(2,i,m)) + &
                                      wl(3,i,m)*( 4.0*g(3,i,m)-s(3,i,m)) )/8.0
        end do

        do i=st-1,ed-1
            vr(i+1,m) = q1(i+1,m) + ( wr(2,i,m)*(-4.0*g(2,i,m)-s(2,i,m)) + &
                                      wr(3,i,m)*(-4.0*g(3,i,m)-s(3,i,m)) + &
                                      wr(4,i,m)*(-4.0*g(4,i,m)-s(4,i,m)) )/8.0
        end do
    end do
    
    do i=st-1,ed-1
        nx = (sn(i,1) + sn(i+1,1))/2.0
        ny = (sn(i,2) + sn(i+1,2))/2.0
        nz = (sn(i,3) + sn(i+1,3))/2.0

        do m=nst,ned
            vp(m) = (vn(i,m) + vn(i+1,m))/2.0
            vc(m) = vl(i,m)
        end do

        call RxQ_p(nst,ned,vp,nt,nx,ny,nz,vc,fp) 
        do m=nst,ned
            vl(i,m) = fp(m)
        end do
        
        do m=nst,ned
            vc(m) = vr(i ,m)
        end do

        call RxQ_p(nst,ned,vp,nt,nx,ny,nz,vc,fp) 
        do m=nst,ned
            vr(i,m) = fp(m)
        end do        
    end do   

    do i=stc,edc
        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vl(i,m) < chklim(m) ) then
                    vl(i,m) = vn(i  ,m)
                    ndec = ndec + 1
                end if

                if (vr(i,m) < chklim(m) ) then
                    vr(i,m) = vn(i+1,m)
                    ndec = ndec + 1
                end if
           end if
        end do
    end do

end subroutine scsn3cicc

subroutine scsl5ci(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_bound_scheme
    implicit none
    integer(kind_int), intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  , intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(in)  :: sn(1-ngn:ni+ngn,1:3)
    integer(kind_int), intent(in)  :: nfs,nfe
    integer(kind_int), intent(in)  :: nge
    real(kind_real)  , intent(out) :: vl(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(out) :: vr(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , external    :: sub_limit
    integer(kind_int), intent(in)  :: chkid(nst:ned)
    real(kind_real),   intent(in)  :: chklim(nst:ned)
    integer(kind_int), intent(out) :: ndec
    real(kind_real)   :: d5(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: ve(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: aa,aa_scs
    integer(kind_int) :: i,m,n,st,ed,stc,edc
    
    if (ni < 8) then
        call muscl2pv(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)
        return
    end if

    ndec = 0
    
    call d5_via_T2(ni,nst,ned,ngn,vn,nfs,nfe,nge,d5)
    
    call vec_via_scsl6(ni,nst,ned,ngn,vn,nfs,nfe,nge,ve)
    
    aa     = 0.01
    aa_scs = 0.001
           
    if (nfs == nbc_bound_scheme) then
        
        st  = 4
        stc = 0

        do m=nst,ned
            vl(0,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                  ! third-order
            vl(1,m) = (  3.0*vn(1,m) +  6.0*vn(2,m) -      vn(3,m) ) / 8.0   ! third-order 
            !!vl(1,m) = (  35.0*vn(1,m) + 140.0*vn(2,m) -  70.0*vn(3,m) +  28.0*vn(4,m) -  5.0*vn(5,m) ) / 128.0   ! fifth-order
            !!vl(2,m) = (  35.0*vn(2,m) + 140.0*vn(3,m) -  70.0*vn(4,m) +  28.0*vn(5,m) -  5.0*vn(6,m) ) / 128.0   ! fifth-order
            vl(2,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 + &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )
            vr(0,m) = vl(0,m)
            vr(1,m) = vl(1,m)
            !!vr(2,m) = vl(2,m)
            vr(2,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 - &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )
            vl(3,m) =   ( 150.0*( vn(4,m) + vn(3,m) )          -  &
                           25.0*( vn(5,m) + vn(2,m) )          +  &
                            3.0*( vn(6,m) + vn(1,m) ) ) /256.0 -  &
                      aa*( 10.0*( vn(4,m) - vn(3,m) )          -  &
                            5.0*( vn(5,m) - vn(2,m) )          +  &
                                ( vn(6,m) - vn(1,m) ) ) 
            
            vr(3,m) =   ( 150.0*( vn(4,m) + vn(3,m) )          -  &
                           25.0*( vn(5,m) + vn(2,m) )          +  &
                            3.0*( vn(6,m) + vn(1,m) ) ) /256.0 +  &
                      aa*( 10.0*( vn(4,m) - vn(3,m) )          -  &
                            5.0*( vn(5,m) - vn(2,m) )          +  &
                                ( vn(6,m) - vn(1,m) ) )             
        end do
        
    else
        
        st  = -2
        stc = st
     
    end if
     
    if (nfe == nbc_bound_scheme) then
        
        ed  = ni - 4
        edc = ni

        do m=nst,ned
            vr(ni  ,m) = ( 15.0*vn(ni,m) - 10.0*vn(ni-1,m) +  3.0*vn(ni-2,m) ) / 8.0                ! third-order
            vl(ni-1,m) = (  3.0*vn(ni,m) +  6.0*vn(ni-1,m) -      vn(ni-2,m) ) / 8.0   ! third-order
            !!vl(ni-1,m) = (  35.0*vn(ni  ,  m) + 140.0*vn(ni-1,m) -  70.0*vn(ni-2,m) +  28.0*vn(ni-3,m) -  5.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            !!vl(ni-2,m) = (  35.0*vn(ni-1,  m) + 140.0*vn(ni-2,m) -  70.0*vn(ni-3,m) +  28.0*vn(ni-4,m) -  5.0*vn(ni-5,m) ) / 128.0   ! fifth-order
            vl(ni-2,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 + &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) )
            vl(ni  ,m) = vr(ni  ,m)
            vr(ni-1,m) = vl(ni-1,m)
            !!vr(ni-2,m) = vl(ni-2,m)
            vr(ni-2,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 - &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) ) 
            vl(ni-3,m) =   ( 150.0*( vn(ni-2,m) + vn(ni-3,m) )          -  &
                              25.0*( vn(ni-1,m) + vn(ni-4,m) )          +  &
                               3.0*( vn(ni  ,m) + vn(ni-5,m) ) ) /256.0 -  &
                         aa*( 10.0*( vn(ni-2,m) - vn(ni-3,m) )          -  &
                               5.0*( vn(ni-1,m) - vn(ni-4,m) )          +  &
                                   ( vn(ni  ,m) - vn(ni-5,m) ) ) 
            
            vr(ni-3,m) =   ( 150.0*( vn(ni-2,m) + vn(ni-3,m) )          -  &
                              25.0*( vn(ni-1,m) + vn(ni-4,m) )          +  &
                               3.0*( vn(ni  ,m) + vn(ni-5,m) ) ) /256.0 +  &
                         aa*( 10.0*( vn(ni-2,m) - vn(ni-3,m) )          -  &
                               5.0*( vn(ni-1,m) - vn(ni-4,m) )          +  &
                                   ( vn(ni  ,m) - vn(ni-5,m) ) )             
        end do
        
    else
        
        ed  = ni + 2
        edc = ed
            
    end if

    do m=nst,ned
        do i=st,ed
            vl(i,m) =   ve(i,m) -  aa_scs*d5(i,m)
            
            vr(i,m) =   ve(i,m) +  aa_scs*d5(i,m)
        end do
    end do

    do i=stc,edc
        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vl(i,m) < chklim(m) ) then
                    vl(i,m) = vn(i,m)
                    ndec = ndec + 1
                end if

                if (vr(i,m) < chklim(m) ) then
                    vr(i,m) = vn(i+1,m)
                    ndec = ndec + 1
                end if
           end if
        end do
    end do
    
end subroutine scsl5ci

subroutine scsl5cicc(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_bound_scheme
    implicit none
    integer(kind_int), intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  , intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(in)  :: sn(1-ngn:ni+ngn,1:3)
    integer(kind_int), intent(in)  :: nfs,nfe
    integer(kind_int), intent(in)  :: nge
    real(kind_real)  , intent(out) :: vl(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(out) :: vr(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , external    :: sub_limit
    integer(kind_int), intent(in)  :: chkid(nst:ned)
    real(kind_real),   intent(in)  :: chklim(nst:ned)
    integer(kind_int), intent(out) :: ndec
    real(kind_real)   :: d5(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: ve(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: aa,aa_scs
    integer(kind_int) :: i,m,n,st,ed,stc,edc
    
    if (ni < 8) then
        call muscl2pvcc(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)
        return
    end if

    ndec = 0
    
    call d5_via_T2(ni,nst,ned,ngn,vn,nfs,nfe,nge,d5)
    
    call vec_via_scsl6(ni,nst,ned,ngn,vn,nfs,nfe,nge,ve)
    
    aa     = 0.01
    aa_scs = 0.001
           
    if (nfs == nbc_bound_scheme) then
        
        st  = 5
        stc = 1

        do m=nst,ned
            vl(1,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                  ! third-order
            vl(2,m) = (  3.0*vn(1,m) +  6.0*vn(2,m) -      vn(3,m) ) / 8.0   ! third-order 
            !!vl(2,m) = (  35.0*vn(1,m) + 140.0*vn(2,m) -  70.0*vn(3,m) +  28.0*vn(4,m) -  5.0*vn(5,m) ) / 128.0   ! fifth-order
            !!vl(3,m) = (  35.0*vn(2,m) + 140.0*vn(3,m) -  70.0*vn(4,m) +  28.0*vn(5,m) -  5.0*vn(6,m) ) / 128.0   ! fifth-order
            vl(3,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 + &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )
            vr(1,m) = vl(1,m)
            vr(2,m) = vl(2,m)
            !!vr(3,m) = vl(3,m)
            vr(3,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 - &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )
            vl(4,m) =   ( 150.0*( vn(4,m) + vn(3,m) )          -  &
                           25.0*( vn(5,m) + vn(2,m) )          +  &
                            3.0*( vn(6,m) + vn(1,m) ) ) /256.0 -  &
                      aa*( 10.0*( vn(4,m) - vn(3,m) )          -  &
                            5.0*( vn(5,m) - vn(2,m) )          +  &
                                ( vn(6,m) - vn(1,m) ) ) 
            
            vr(4,m) =   ( 150.0*( vn(4,m) + vn(3,m) )          -  &
                           25.0*( vn(5,m) + vn(2,m) )          +  &
                            3.0*( vn(6,m) + vn(1,m) ) ) /256.0 +  &
                      aa*( 10.0*( vn(4,m) - vn(3,m) )          -  &
                            5.0*( vn(5,m) - vn(2,m) )          +  &
                                ( vn(6,m) - vn(1,m) ) )             
        end do
        
    else
        
        st  = -1
        stc = st
     
    end if
     
    if (nfe == nbc_bound_scheme) then
        
        ed  = ni - 3
        edc = ni + 1

        do m=nst,ned
            vr(ni+1,m) = ( 15.0*vn(ni,m) - 10.0*vn(ni-1,m) +  3.0*vn(ni-2,m) ) / 8.0                ! third-order
            vl(ni  ,m) = (  3.0*vn(ni,m) +  6.0*vn(ni-1,m) -      vn(ni-2,m) ) / 8.0   ! third-order
            !!vl(ni  ,m) = (  35.0*vn(ni  ,  m) + 140.0*vn(ni-1,m) -  70.0*vn(ni-2,m) +  28.0*vn(ni-3,m) -  5.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            !!vl(ni-1,m) = (  35.0*vn(ni-1,  m) + 140.0*vn(ni-2,m) -  70.0*vn(ni-3,m) +  28.0*vn(ni-4,m) -  5.0*vn(ni-5,m) ) / 128.0   ! fifth-order
            vl(ni-1,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 + &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) )
            vl(ni+1,m) = vr(ni+1,m)
            vr(ni  ,m) = vl(ni  ,m)
            !!vr(ni-1,m) = vl(ni-1,m)
            vr(ni-1,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 - &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) ) 
            vl(ni-2,m) =   ( 150.0*( vn(ni-2,m) + vn(ni-3,m) )          -  &
                              25.0*( vn(ni-1,m) + vn(ni-4,m) )          +  &
                               3.0*( vn(ni  ,m) + vn(ni-5,m) ) ) /256.0 -  &
                         aa*( 10.0*( vn(ni-2,m) - vn(ni-3,m) )          -  &
                               5.0*( vn(ni-1,m) - vn(ni-4,m) )          +  &
                                   ( vn(ni  ,m) - vn(ni-5,m) ) ) 
            
            vr(ni-2,m) =   ( 150.0*( vn(ni-2,m) + vn(ni-3,m) )          -  &
                              25.0*( vn(ni-1,m) + vn(ni-4,m) )          +  &
                               3.0*( vn(ni  ,m) + vn(ni-5,m) ) ) /256.0 +  &
                         aa*( 10.0*( vn(ni-2,m) - vn(ni-3,m) )          -  &
                               5.0*( vn(ni-1,m) - vn(ni-4,m) )          +  &
                                   ( vn(ni  ,m) - vn(ni-5,m) ) )             
        end do
        
    else
        
        ed  = ni + 3
        edc = ed
            
    end if

    do m=nst,ned
        do i=st-1,ed-1
            vl(i,m) =   ve(i,m) -  aa_scs*d5(i,m)
            
            vr(i,m) =   ve(i,m) +  aa_scs*d5(i,m)
        end do
    end do

    do i=stc,edc
        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vl(i,m) < chklim(m) ) then
                    vl(i,m) = vn(i,m)
                    ndec = ndec + 1
                end if

                if (vr(i,m) < chklim(m) ) then
                    vr(i,m) = vn(i+1,m)
                    ndec = ndec + 1
                end if
           end if
        end do
    end do
    
end subroutine scsl5cicc

subroutine d5_via_T2(ni,nst,ned,ngn,vn,nfs,nfe,nge,d5)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_bound_scheme
    implicit none
    integer(kind_int),intent(in)  :: ni,nst,ned,ngn
    integer(kind_int),intent(in)  :: nfs,nfe
    integer(kind_int),intent(in)  :: nge    
    real(kind_real)  ,intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  ,intent(out) :: d5(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: ma(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵ����Aϵ��
    real(kind_real)   :: mb(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵ����Bϵ��
    real(kind_real)   :: mc(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵ����Cϵ��
    real(kind_real)   :: md(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵ��Ҳ�Dϵ��
    real(kind_real)   :: ms(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵļ�����
    real(kind_real)   :: x1                                        !< �ع���ʽϵ��
    integer(kind_int) :: i,m,n,st,ed,stc,edc
    
    x1=0.46
    
    !>
    !!������׾��Ȱ�����׿ռ䵼��
    if (nfs == nbc_bound_scheme) then

        st  = 4       !< ����nfs=1ʱ���ڵ��ʽ�������ʼ����
        stc = 3       !< ����nfs=1ʱ���������ʼ����

        ma(st-1)=0.0
        mc(st-1)=0.0

        do m=nst,ned
            d5(3,m) =   10.0*( vn(4,m) - vn(3,m) )          -  &
                         5.0*( vn(5,m) - vn(2,m) )          +  &
                             ( vn(6,m) - vn(1,m) )              !< �����ع���߽�
        end do 
        
    else
        
        st  = -2   !< ����nfs=0������ȫ�ڵ�߽��ʽ��ʱ���ڵ��ع���ʽ�������ʼ����
        stc = -3   !< ����nfs=0������ȫ�ڵ�߽��ʽ��ʱ���ع��������ʼ����
        
        ma(st-1)=0.0                           !< ���Խ�׷�ϵ����Aϵ������߽磩
        mc(st-1)=0.0                           !< ���Խ�׷�ϵ����Cϵ������߽磩
        
        do m=nst,ned
            d5(-3,m) =  10.0*( vn(-2,m) - vn(-3,m) )          -  &
                         5.0*( vn(-1,m) - vn(-4,m) )          +  &
                             ( vn( 0,m) - vn(-5,m) )              !< @attention �ع���δʹ�ã�ֻ��Ϊ���������Խ�׷��
        end do        
 
    end if
    !>

    if (nfe == nbc_bound_scheme) then
        
        ed  = ni - 4      !< ����nfe=1ʱ���ڵ��ֵ��ʽ�������ֹ����
        edc = ni - 3      !< ����nfe=1ʱ����ֵ�������ֹ����

        ma(ed+1)=0.0                           !< ���Խ�׷�ϵ����Aϵ�����ұ߽磩
        mc(ed+1)=0.0                           !< ���Խ�׷�ϵ����Cϵ�����ұ߽磩

        do m=nst,ned
            d5(ni-3,m) =  10.0*( vn(ni-2,m) - vn(ni-3,m) )          -  &
                           5.0*( vn(ni-1,m) - vn(ni-4,m) )          +  &
                               ( vn(ni  ,m) - vn(ni-5,m) )                      !<  �����ع��ұ߽�                             
        end do
        
    else
        
        ed  = ni+2      !< ����nfe=0������ȫ�ڵ�߽��ʽ��ʱ���ڵ��ع���ʽ�������ֹ����
        edc = ni+3       !< ����nfe=0������ȫ�ڵ�߽��ʽ��ʱ���ع��������ֹ����
        
        ma(ed+1)=0.0
        mc(ed+1)=0.0

        do m=nst,ned
            d5(ni+3,m) =  10.0*( vn(ni+4,m) - vn(ni+3,m) )          -  &
                           5.0*( vn(ni+5,m) - vn(ni+2,m) )          +  &
                               ( vn(ni+6,m) - vn(ni+1,m) )                         !< @attention ��ֵ��δʹ�ã�ֻ��Ϊ���������Խ�׷��
        end do        
   
    end if
    !>

    !>
    !!������׾��Ȱ�����׿ռ䵼���ڵ�
    do m=nst,ned
        do i=st,ed
            d5(i,m) =     10.0*(1.0+2.0*x1)*( vn(i+1,m) - vn(i,m  ) )          -  &
                           5.0*(1.0+2.0*x1)*( vn(i+2,m) - vn(i-1,m) )          +  &
                               (1.0+2.0*x1)*( vn(i+3,m) - vn(i-2,m) )
        end do
    end do
    
    do i=stc,edc
        mb(i)=1.0                          !< ���Խ�׷�ϵ�Bϵ��
    enddo
  
    do i=st,ed
        ma(i)=x1                           !< ���Խ�׷�ϵ�Aϵ�����ڵ㣩
    end do

    do i=st,ed                                                
        mc(i)=x1                           !< ���Խ�׷�ϵ�Cϵ�����ڵ㣩
    end do                                               


    
    do m=nst,ned
        do i=stc,edc
            md(i) = d5(i,m)               !< ���Խ�׷�ϵ��Ҳ�Dϵ��
        end do
        call tri_pursue(ma,mb,mc,md,ms,stc,edc,1-ngn,ni+ngn) !< ���Խ�׷��
        do i=stc,edc
            d5(i,m) = ms(i)               !< �������Խ�׷�Ͻ�
        end do
    end do    
    
end subroutine d5_via_T2

subroutine d5_via_T2N(ni,nst,ned,ngn,vn,nfs,nfe,nge,d5,st,ed)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_bound_scheme
    implicit none
    integer(kind_int),intent(in)  :: ni,nst,ned,ngn
    integer(kind_int),intent(in)  :: nfs,nfe
    integer(kind_int),intent(in)  :: nge,st,ed    
    real(kind_real)  ,intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  ,intent(out) :: d5(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: ma(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵ����Aϵ��
    real(kind_real)   :: mb(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵ����Bϵ��
    real(kind_real)   :: mc(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵ����Cϵ��
    real(kind_real)   :: md(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵ��Ҳ�Dϵ��
    real(kind_real)   :: ms(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵļ�����
    real(kind_real)   :: x1                                        !< �ع���ʽϵ��
    integer(kind_int) :: i,m
    
    x1=0.46

    ma(st)=0.0
    mc(st)=0.0

    do m=nst,ned
        d5(st,m) =   10.0*( vn(st+1,m) - vn(st  ,m) )          -  &
                      5.0*( vn(st+2,m) - vn(st-1,m) )          +  &
                          ( vn(st+3,m) - vn(st-2,m) )
    end do

    ma(ed)=0.0
    mc(ed)=0.0

    do m=nst,ned
        d5(ed,m) =  10.0*( vn(ed+1,m) - vn(ed  ,m) )          -  &
                     5.0*( vn(ed+2,m) - vn(ed-1,m) )          +  &
                         ( vn(ed+3,m) - vn(ed-2,m) )                     
    end do
        
    do m=nst,ned
        do i=st+1,ed-1
            d5(i,m) =     10.0*(1.0+2.0*x1)*( vn(i+1,m) - vn(i,m  ) )          -  &
                           5.0*(1.0+2.0*x1)*( vn(i+2,m) - vn(i-1,m) )          +  &
                               (1.0+2.0*x1)*( vn(i+3,m) - vn(i-2,m) )
        end do
    end do
    
    do i=st,ed
        mb(i)=1.0                          !< ���Խ�׷�ϵ�Bϵ��
    enddo
  
    do i=st+1,ed-1
        ma(i)=x1                           !< ���Խ�׷�ϵ�Aϵ�����ڵ㣩
    end do

    do i=st+1,ed-1
        mc(i)=x1                           !< ���Խ�׷�ϵ�Cϵ�����ڵ㣩
    end do                                               


    
    do m=nst,ned
        do i=st,ed
            md(i) = d5(i,m)               !< ���Խ�׷�ϵ��Ҳ�Dϵ��
        end do
        call tri_pursue(ma,mb,mc,md,ms,st,ed,1-ngn,ni+ngn) !< ���Խ�׷��
        do i=st,ed
            d5(i,m) = ms(i)               !< �������Խ�׷�Ͻ�
        end do
    end do    
    
end subroutine d5_via_T2N

subroutine vec_via_scsl6(ni,nst,ned,ngn,vn,nfs,nfe,nge,ve)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_bound_scheme
    implicit none
    integer(kind_int),intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  ,intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)  :: nfs,nfe
    integer(kind_int),intent(in)  :: nge
    real(kind_real)  ,intent(out) :: ve(1-ngn:ni+ngn,nst:ned)
    real(kind_real)               :: ve2(1-ngn:ni+ngn,nst:ned)
    real(kind_real)               :: ve4(1-ngn:ni+ngn,nst:ned)
    integer(kind_int) :: i,m,st,ed,ierr
    
    if (ni < 8) then
        call ve_via_node2(ni,nst,ned,ngn,vn,nfs,nfe,nge,ve)
        return
    end if    
    
    call vec2_via_T4(ni,nst,ned,ngn,vn,nfs,nfe,nge,ve2)
    
    call vec4_via_T2(ni,nst,ned,ngn,vn,nfs,nfe,nge,ve4)
    
    if (nfs == nbc_bound_scheme) then
        
        st = 4

        do m=nst,ned
            !!ve(2,m) = (      -vn(1,m) +   9.0*vn(2,m) +   9.0*vn(3,m) -       vn(4,m) ) / 16.0                   ! fourth-order
            !!ve(1,m) = (   5.0*vn(1,m) +  15.0*vn(2,m) -   5.0*vn(3,m) +       vn(4,m) ) / 16.0                   ! fourth-order
            !!ve(0,m) = (  35.0*vn(1,m) -  35.0*vn(2,m) +  21.0*vn(3,m) -   5.0*vn(4,m) ) / 16.0                   ! fourth-order
            ve(0,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                                       ! third-order
            ve(2,m) = (  -5.0*vn(1,m) +  60.0*vn(2,m) +  90.0*vn(3,m) -  20.0*vn(4,m) +  3.0*vn(5,m) ) / 128.0   ! fifth-order
            ve(1,m) = (  35.0*vn(1,m) + 140.0*vn(2,m) -  70.0*vn(3,m) +  28.0*vn(4,m) -  5.0*vn(5,m) ) / 128.0   ! fifth-order
            ve(3 ,m) = (150.0*(vn(4,m) + vn(3,m)) - 25.0*(vn(5,m) + vn(2,m)) + 3.0*(vn(6,m) + vn(1,m)) )/256.0   ! sixth-order            
            !!ve(0,m) = ( 315.0*vn(1,m) - 420.0*vn(2,m) + 378.0*vn(3,m) - 180.0*vn(4,m) + 35.0*vn(5,m) ) / 128.0   ! fifth-order
            !!ve(1,m) = ( 3.0*vn(1,m) + 6.0*vn(2,m) -     vn(3,m) ) / 8.0                                          ! third-order
            !!ve(2,m) = (    -vn(1,m) + 9.0*vn(2,m) + 9.0*vn(3,m) - vn(4,m) ) / 16.0                               ! fourth-order
        end do
        
    else
        
        st = -2
        
    end if
    
    if (nfe == nbc_bound_scheme) then
        
        ed = ni - 4

        do m=nst,ned
            !!ve(ni-2,m) = (      -vn(ni,m) +   9.0*vn(ni-1,m) +   9.0*vn(ni-2,m) -       vn(ni-3,m) ) / 16.0                      ! fourth-order
            !!ve(ni-1,m) = (   5.0*vn(ni,m) +  15.0*vn(ni-1,m) -   5.0*vn(ni-2,m) +       vn(ni-3,m) ) / 16.0                      ! fourth-order
            !!ve(ni  ,m) = (  35.0*vn(ni,m) -  35.0*vn(ni-1,m) +  21.0*vn(ni-2,m) -   5.0*vn(ni-3,m) ) / 16.0                      ! fourth-order
            ve(ni  ,m) = (  15.0*vn(ni,m) -  10.0*vn(ni-1,m) +   3.0*vn(ni-2,m) ) / 8.0                                            ! third-order
            ve(ni-2,m) = (  -5.0*vn(ni,m) +  60.0*vn(ni-1,m) +  90.0*vn(ni-2,m) -  20.0*vn(ni-3,m) +  3.0*vn(ni-4,m) ) / 128.0     ! fifth-order
            ve(ni-1,m) = (  35.0*vn(ni,m) + 140.0*vn(ni-1,m) -  70.0*vn(ni-2,m) +  28.0*vn(ni-3,m) -  5.0*vn(ni-4,m) ) / 128.0     ! fifth-order
            ve(ni-3,m) = ( 150.0*(vn(ni-2,m) + vn(ni-3,m)) - 25.0*(vn(ni-1,m) + vn(ni-4,m)) + 3.0*(vn(ni ,m) + vn(ni-5,m)) )/256.0 ! sixth-order            
            !!ve(ni  ,m) = ( 315.0*vn(ni,m) - 420.0*vn(ni-1,m) + 378.0*vn(ni-2,m) - 180.0*vn(ni-3,m) + 35.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            !!ve(ni-1,m) = (  3.0*vn(ni,m) +   6.0*vn(ni-1,m) -       vn(ni-2,m) ) / 8.0                                          ! third-order
            !!ve(ni-2,m) = (     -vn(ni,m) +   9.0*vn(ni-1,m) +   9.0*vn(ni-2,m) -       vn(ni-3,m) ) / 16.0                      ! fourth-order
        end do
        
    else
    
        ed = ni + 2 
      
    end if

    do m=nst,ned
        do i=st,ed
            ve(i,m) =    0.5*(vn(i+1,m) + vn(i,m)) - &
                                 1.0/8.0*ve2(i,m)  - &
                               1.0/384.0*ve4(i,m)
        end do
    end do

end subroutine vec_via_scsl6

subroutine vec_via_scsl6cc(ni,nst,ned,ngn,vn,nfs,nfe,nge,ve)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_bound_scheme
    implicit none
    integer(kind_int),intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  ,intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)  :: nfs,nfe
    integer(kind_int),intent(in)  :: nge
    real(kind_real)  ,intent(out) :: ve(1-ngn:ni+ngn,nst:ned)
    real(kind_real)               :: ve2(1-ngn:ni+ngn,nst:ned)
    real(kind_real)               :: ve4(1-ngn:ni+ngn,nst:ned)
    integer(kind_int) :: i,m,st,ed,ierr
    
    if (ni < 8) then
        call ve_via_node2cc(ni,nst,ned,ngn,vn,nfs,nfe,nge,ve)
        return
    end if    
    
    call vec2_via_T4(ni,nst,ned,ngn,vn,nfs,nfe,nge,ve2)
    
    call vec4_via_T2(ni,nst,ned,ngn,vn,nfs,nfe,nge,ve4)
    
    if (nfs == nbc_bound_scheme) then
        
        st = 5

        do m=nst,ned
            !!ve(3,m) = (      -vn(1,m) +   9.0*vn(2,m) +   9.0*vn(3,m) -       vn(4,m) ) / 16.0                   ! fourth-order
            !!ve(2,m) = (   5.0*vn(1,m) +  15.0*vn(2,m) -   5.0*vn(3,m) +       vn(4,m) ) / 16.0                   ! fourth-order
            !!ve(1,m) = (  35.0*vn(1,m) -  35.0*vn(2,m) +  21.0*vn(3,m) -   5.0*vn(4,m) ) / 16.0                   ! fourth-order
            ve(1,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                                       ! third-order
            ve(3,m) = (  -5.0*vn(1,m) +  60.0*vn(2,m) +  90.0*vn(3,m) -  20.0*vn(4,m) +  3.0*vn(5,m) ) / 128.0   ! fifth-order
            ve(2,m) = (  35.0*vn(1,m) + 140.0*vn(2,m) -  70.0*vn(3,m) +  28.0*vn(4,m) -  5.0*vn(5,m) ) / 128.0   ! fifth-order
            ve(4 ,m) = (150.0*(vn(4,m) + vn(3,m)) - 25.0*(vn(5,m) + vn(2,m)) + 3.0*(vn(6,m) + vn(1,m)) )/256.0   ! sixth-order            
            !!ve(1,m) = ( 315.0*vn(1,m) - 420.0*vn(2,m) + 378.0*vn(3,m) - 180.0*vn(4,m) + 35.0*vn(5,m) ) / 128.0   ! fifth-order
            !!ve(2,m) = ( 3.0*vn(1,m) + 6.0*vn(2,m) -     vn(3,m) ) / 8.0                                          ! third-order
            !!ve(3,m) = (    -vn(1,m) + 9.0*vn(2,m) + 9.0*vn(3,m) - vn(4,m) ) / 16.0                               ! fourth-order
        end do
        
    else
        
        st = -1
        
    end if
    
    if (nfe == nbc_bound_scheme) then
        
        ed = ni - 3

        do m=nst,ned
            !!ve(ni-1,m) = (      -vn(ni,m) +   9.0*vn(ni-1,m) +   9.0*vn(ni-2,m) -       vn(ni-3,m) ) / 16.0                      ! fourth-order
            !!ve(ni  ,m) = (   5.0*vn(ni,m) +  15.0*vn(ni-1,m) -   5.0*vn(ni-2,m) +       vn(ni-3,m) ) / 16.0                      ! fourth-order
            !!ve(ni+1,m) = (  35.0*vn(ni,m) -  35.0*vn(ni-1,m) +  21.0*vn(ni-2,m) -   5.0*vn(ni-3,m) ) / 16.0                      ! fourth-order
            ve(ni+1,m) = (  15.0*vn(ni,m) -  10.0*vn(ni-1,m) +   3.0*vn(ni-2,m) ) / 8.0                                            ! third-order
            ve(ni-1,m) = (  -5.0*vn(ni,m) +  60.0*vn(ni-1,m) +  90.0*vn(ni-2,m) -  20.0*vn(ni-3,m) +  3.0*vn(ni-4,m) ) / 128.0     ! fifth-order
            ve(ni  ,m) = (  35.0*vn(ni,m) + 140.0*vn(ni-1,m) -  70.0*vn(ni-2,m) +  28.0*vn(ni-3,m) -  5.0*vn(ni-4,m) ) / 128.0     ! fifth-order
            ve(ni-2,m) = ( 150.0*(vn(ni-2,m) + vn(ni-3,m)) - 25.0*(vn(ni-1,m) + vn(ni-4,m)) + 3.0*(vn(ni ,m) + vn(ni-5,m)) )/256.0 ! sixth-order            
            !!ve(ni+1,m) = ( 315.0*vn(ni,m) - 420.0*vn(ni-1,m) + 378.0*vn(ni-2,m) - 180.0*vn(ni-3,m) + 35.0*vn(ni-4,m) ) / 128.0   ! fifth-order
            !!ve(ni  ,m) = (  3.0*vn(ni,m) +   6.0*vn(ni-1,m) -       vn(ni-2,m) ) / 8.0                                          ! third-order
            !!ve(ni-1,m) = (     -vn(ni,m) +   9.0*vn(ni-1,m) +   9.0*vn(ni-2,m) -       vn(ni-3,m) ) / 16.0                      ! fourth-order
        end do
        
    else
    
        ed = ni + 3 
      
    end if

    do m=nst,ned
        do i=st-1,ed-1
            ve(i+1,m) =    0.5*(vn(i+1,m) + vn(i,m)) - &
                                   1.0/8.0*ve2(i,m)  - &
                                 1.0/384.0*ve4(i,m)
        end do
    end do

end subroutine vec_via_scsl6cc

subroutine vec2_via_T4(ni,nst,ned,ngn,vn,nfs,nfe,nge,ve2)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_bound_scheme
    implicit none
    integer(kind_int),intent(in)  :: ni,nst,ned,ngn
    integer(kind_int),intent(in)  :: nfs,nfe
    integer(kind_int),intent(in)  :: nge    
    real(kind_real)  ,intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  ,intent(out) :: ve2(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: ma(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵ����Aϵ��
    real(kind_real)   :: mb(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵ����Bϵ��
    real(kind_real)   :: mc(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵ����Cϵ��
    real(kind_real)   :: md(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵ��Ҳ�Dϵ��
    real(kind_real)   :: ms(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵļ�����
    real(kind_real)   :: r1                                        !< �ع���ʽϵ��
    integer(kind_int) :: i,m,n,st,ed,stc,edc
    
    r1=0.46
    
    !>
    !!������׾��Ȱ����׿ռ䵼��                                                     
    if (nfs == nbc_bound_scheme) then

        st  = 4       !< ����nfs=1ʱ���ڵ��ʽ�������ʼ����
        stc = 3       !< ����nfs=1ʱ���������ʼ����

        ma(st-1)=0.0
        mc(st-1)=0.0

        do m=nst,ned
            ve2(3,m) = -17.0/24.0*(vn(4,m)+vn(3,m)) + &
                        13.0/16.0*(vn(5,m)+vn(2,m)) - &
                         5.0/48.0*(vn(6,m)+vn(1,m))       !< �Ľ��ع���߽�
        end do
        
    else
        
        st  = -2     !< ����nfs=0������ȫ�ڵ�߽��ʽ��ʱ���ڵ��ع���ʽ�������ʼ����
        stc = -3     !< ����nfs=0������ȫ�ڵ�߽��ʽ��ʱ���ع��������ʼ����
        
        ma(st-1)=0.0                           !< ���Խ�׷�ϵ����Aϵ������߽磩
        mc(st-1)=0.0                           !< ���Խ�׷�ϵ����Cϵ������߽磩
        
        do m=nst,ned
            ve2(-3,m) = -17.0/24.0*(vn(-2,m)+vn(-3,m)) + &
                         13.0/16.0*(vn(-1,m)+vn(-4,m)) - &
                          5.0/48.0*(vn( 0,m)+vn(-5,m))          !< @attention �ع���δʹ�ã�ֻ��Ϊ���������Խ�׷��
        end do         
 
    end if
    !>

    !>        
    if (nfe == nbc_bound_scheme) then
        
        ed  = ni - 4      !< ����nfe=1ʱ���ڵ��ֵ��ʽ�������ֹ����
        edc = ni - 3      !< ����nfe=1ʱ����ֵ�������ֹ����

        ma(ed+1)=0.0                           !< ���Խ�׷�ϵ����Aϵ�����ұ߽磩
        mc(ed+1)=0.0                           !< ���Խ�׷�ϵ����Cϵ�����ұ߽磩

        do m=nst,ned
            ve2(ni-3,m) = -17.0/24.0*(vn(ni-2,m)+vn(ni-3,m)) + &
                           13.0/16.0*(vn(ni-1,m)+vn(ni-4,m)) - &
                            5.0/48.0*(vn(ni  ,m)+vn(ni-5,m))                   !<  �Ľ��ع��ұ߽�                             
        end do
        
    else
        
        ed  = ni+2        !< ����nfe=0������ȫ�ڵ�߽��ʽ��ʱ���ڵ��ع���ʽ�������ֹ����
        edc = ni+3        !< ����nfe=0������ȫ�ڵ�߽��ʽ��ʱ���ع��������ֹ����
        
        ma(ed+1)=0.0
        mc(ed+1)=0.0

        do m=nst,ned
            ve2(ni+3,m) = -17.0/24.0*(vn(ni+4,m)+vn(ni+3,m)) + &
                           13.0/16.0*(vn(ni+5,m)+vn(ni+2,m)) - &
                            5.0/48.0*(vn(ni+6,m)+vn(ni+1,m))    !< @attention ��ֵ��δʹ�ã�ֻ��Ϊ���������Խ�׷��
        end do         

    end if
    !>

    !>
    !!������׾��Ȱ����׿ռ䵼���ڵ�
    do m=nst,ned
        do i=st,ed
            ve2(i,m) =    (-17.0-10.0*r1)/24.0*( vn(i+1,m) + vn(i  ,m) )          +  &
                            (13.0+2.0*r1)/16.0*( vn(i+2,m) + vn(i-1,m) )          +  &
                           (-5.0+14.0*r1)/48.0*( vn(i+3,m) + vn(i-2,m) )
        end do
    end do
    
    do i=stc,edc
        mb(i)=1.0                          !< ���Խ�׷�ϵ�Bϵ��
    enddo
  
    do i=st,ed
        ma(i)=r1                           !< ���Խ�׷�ϵ�Aϵ�����ڵ㣩
    end do

    do i=st,ed                                                
        mc(i)=r1                           !< ���Խ�׷�ϵ�Cϵ�����ڵ㣩
    end do                                               


    
    do m=nst,ned
        do i=stc,edc
            md(i) = ve2(i,m)               !< ���Խ�׷�ϵ��Ҳ�Dϵ��
        end do
        call tri_pursue(ma,mb,mc,md,ms,stc,edc,1-ngn,ni+ngn) !< ���Խ�׷��
        do i=stc,edc
            ve2(i,m) = ms(i)               !< �������Խ�׷�Ͻ�
        end do
    end do    
    
end subroutine vec2_via_T4

subroutine vec2_via_T4N(ni,nst,ned,ngn,vn,nfs,nfe,nge,ve2,st,ed)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_bound_scheme
    implicit none
    integer(kind_int),intent(in)  :: ni,nst,ned,ngn
    integer(kind_int),intent(in)  :: nfs,nfe
    integer(kind_int),intent(in)  :: nge,st,ed    
    real(kind_real)  ,intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  ,intent(out) :: ve2(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: ma(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵ����Aϵ��
    real(kind_real)   :: mb(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵ����Bϵ��
    real(kind_real)   :: mc(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵ����Cϵ��
    real(kind_real)   :: md(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵ��Ҳ�Dϵ��
    real(kind_real)   :: ms(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵļ�����
    real(kind_real)   :: r1                                        !< �ع���ʽϵ��
    integer(kind_int) :: i,m
    
    r1=0.46

    ma(st)=0.0
    mc(st)=0.0
    
    do m=nst,ned
        ve2(st,m) = -17.0/24.0*(vn(st+1,m)+vn(st  ,m)) + &
                     13.0/16.0*(vn(st+2,m)+vn(st-1,m)) - &
                      5.0/48.0*(vn(st+3,m)+vn(st-2,m))       !< �Ľ��ع���߽�
    end do

    ma(ed)=0.0
    mc(ed)=0.0

    do m=nst,ned
        ve2(ed,m) = -17.0/24.0*(vn(ed+1,m)+vn(ed  ,m)) + &
                     13.0/16.0*(vn(ed+2,m)+vn(ed-1,m)) - &
                      5.0/48.0*(vn(ed+3,m)+vn(ed-2,m))
    end do
    
    do m=nst,ned
        do i=st+1,ed-1
            ve2(i,m) =    (-17.0-10.0*r1)/24.0*( vn(i+1,m) + vn(i  ,m) )          +  &
                            (13.0+2.0*r1)/16.0*( vn(i+2,m) + vn(i-1,m) )          +  &
                           (-5.0+14.0*r1)/48.0*( vn(i+3,m) + vn(i-2,m) )
        end do
    end do
    
    do i=st,ed
        mb(i)=1.0                          !< ���Խ�׷�ϵ�Bϵ��
    enddo
  
    do i=st+1,ed-1
        ma(i)=r1                           !< ���Խ�׷�ϵ�Aϵ�����ڵ㣩
    end do

    do i=st+1,ed-1
        mc(i)=r1                           !< ���Խ�׷�ϵ�Cϵ�����ڵ㣩
    end do                                               


    
    do m=nst,ned
        do i=st,ed
            md(i) = ve2(i,m)               !< ���Խ�׷�ϵ��Ҳ�Dϵ��
        end do
        call tri_pursue(ma,mb,mc,md,ms,st,ed,1-ngn,ni+ngn) !< ���Խ�׷��
        do i=st,ed
            ve2(i,m) = ms(i)               !< �������Խ�׷�Ͻ�
        end do
    end do    
    
end subroutine vec2_via_T4N

subroutine vec4_via_T2(ni,nst,ned,ngn,vn,nfs,nfe,nge,ve4)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_bound_scheme
    implicit none
    integer(kind_int),intent(in)  :: ni,nst,ned,ngn
    integer(kind_int),intent(in)  :: nfs,nfe
    integer(kind_int),intent(in)  :: nge    
    real(kind_real)  ,intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  ,intent(out) :: ve4(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: ma(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵ����Aϵ��
    real(kind_real)   :: mb(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵ����Bϵ��
    real(kind_real)   :: mc(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵ����Cϵ��
    real(kind_real)   :: md(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵ��Ҳ�Dϵ��
    real(kind_real)   :: ms(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵļ�����
    real(kind_real)   :: r1                                        !< �ع���ʽϵ��
    integer(kind_int) :: i,m,n,st,ed,stc,edc
    
    r1=0.24
    
    !>
    !!������׾��Ȱ����׿ռ䵼��                                                     
    if (nfs == nbc_bound_scheme) then

        st  = 4       !< ����nfs=1ʱ���ڵ��ʽ�������ʼ����
        stc = 3       !< ����nfs=1ʱ���������ʼ����

        ma(st-1)=0.0
        mc(st-1)=0.0

        do m=nst,ned
            ve4(3,m) = (vn(4,m)+vn(3,m))-3.0/2.0*(vn(5,m)+vn(2,m))+1.0/2.0*(vn(6,m)+vn(1,m))      !< �����ع���߽�
        end do
        
    else
        
        st  = -2    !< ����nfs=0������ȫ�ڵ�߽��ʽ��ʱ���ڵ��ع���ʽ�������ʼ����
        stc = -3    !< ����nfs=0������ȫ�ڵ�߽��ʽ��ʱ���ع��������ʼ����
        
        ma(st-1)=0.0                           !< ���Խ�׷�ϵ����Aϵ������߽磩
        mc(st-1)=0.0                           !< ���Խ�׷�ϵ����Cϵ������߽磩
        
        do m=nst,ned
            ve4(-3,m) = (vn(-2,m)+vn(-3,m))-3.0/2.0*(vn(-1,m)+vn(-4,m))+1.0/2.0*(vn(0,m)+vn(-5,m))   !< @attention �ع���δʹ�ã�ֻ��Ϊ���������Խ�׷��
        end do         
 
    end if
    !>

    !>        
    if (nfe == nbc_bound_scheme) then
        
        ed  = ni - 4     !< ����nfe=1ʱ���ڵ��ֵ��ʽ�������ֹ����
        edc = ni - 3     !< ����nfe=1ʱ����ֵ�������ֹ����

        ma(ed+1)=0.0                           !< ���Խ�׷�ϵ����Aϵ�����ұ߽磩
        mc(ed+1)=0.0                           !< ���Խ�׷�ϵ����Cϵ�����ұ߽磩

        do m=nst,ned
            ve4(ni-3,m) = (vn(ni-2,m)+vn(ni-3,m))-3.0/2.0*(vn(ni-1,m)+vn(ni-4,m))+1.0/2.0*(vn(ni,m)+vn(ni-5,m))         !<  �����ع��ұ߽�                             
        end do
        
    else
        
        ed  = ni+2        !< ����nfe=0������ȫ�ڵ�߽��ʽ��ʱ���ڵ��ع���ʽ�������ֹ����
        edc = ni+3        !< ����nfe=0������ȫ�ڵ�߽��ʽ��ʱ���ع��������ֹ����
        
        ma(ed+1)=0.0
        mc(ed+1)=0.0

        do m=nst,ned
            ve4(ni+3,m) = (vn(ni+4,m)+vn(ni+3,m))-3.0/2.0*(vn(ni+5,m)+vn(ni+2,m))+1.0/2.0*(vn(ni+6,m)+vn(ni+1,m))   !< @attention ��ֵ��δʹ�ã�ֻ��Ϊ���������Խ�׷��
        end do         

    end if
    !>

    !>
    !!������׾��Ȱ����׿ռ䵼���ڵ�
    do m=nst,ned
        do i=st,ed
            ve4(i,m) =          (1.0+2.0*r1)*( vn(i+1,m) + vn(i,m  ) )          -  &
                        3.0/2.0*(1.0+2.0*r1)*( vn(i+2,m) + vn(i-1,m) )          +  &
                        1.0/2.0*(1.0+2.0*r1)*( vn(i+3,m) + vn(i-2,m) )
        end do
    end do
    
    do i=stc,edc
        mb(i)=1.0                          !< ���Խ�׷�ϵ�Bϵ��
    enddo
  
    do i=st,ed
        ma(i)=r1                           !< ���Խ�׷�ϵ�Aϵ�����ڵ㣩
    end do

    do i=st,ed                                                
        mc(i)=r1                           !< ���Խ�׷�ϵ�Cϵ�����ڵ㣩
    end do                                               


    
    do m=nst,ned
        do i=stc,edc
            md(i) = ve4(i,m)               !< ���Խ�׷�ϵ��Ҳ�Dϵ��
        end do
        call tri_pursue(ma,mb,mc,md,ms,stc,edc,1-ngn,ni+ngn) !< ���Խ�׷��
        do i=stc,edc
            ve4(i,m) = ms(i)               !< �������Խ�׷�Ͻ�
        end do
    end do    
    
end subroutine vec4_via_T2

subroutine vec4_via_T2N(ni,nst,ned,ngn,vn,nfs,nfe,nge,ve4,st,ed)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_bound_scheme
    implicit none
    integer(kind_int),intent(in)  :: ni,nst,ned,ngn
    integer(kind_int),intent(in)  :: nfs,nfe
    integer(kind_int),intent(in)  :: nge,st,ed    
    real(kind_real)  ,intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  ,intent(out) :: ve4(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: ma(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵ����Aϵ��
    real(kind_real)   :: mb(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵ����Bϵ��
    real(kind_real)   :: mc(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵ����Cϵ��
    real(kind_real)   :: md(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵ��Ҳ�Dϵ��
    real(kind_real)   :: ms(1-ngn:ni+ngn)                          !< �������Խ�׷�ϵļ�����
    real(kind_real)   :: r1                                        !< �ع���ʽϵ��
    integer(kind_int) :: i,m
    
    r1=0.24

    ma(st)=0.0
    mc(st)=0.0

    do m=nst,ned
        ve4(st,m) = (vn(st+1,m)+vn(st,m))-3.0/2.0*(vn(st+2,m)+vn(st-1,m))+1.0/2.0*(vn(st+3,m)+vn(st-2,m))      !< �����ع���߽�
    end do

    ma(ed)=0.0                           !< ���Խ�׷�ϵ����Aϵ�����ұ߽磩
    mc(ed)=0.0                           !< ���Խ�׷�ϵ����Cϵ�����ұ߽磩

    do m=nst,ned
        ve4(ed,m) = (vn(ed+1,m)+vn(ed,m))-3.0/2.0*(vn(ed+2,m)+vn(ed-1,m))+1.0/2.0*(vn(ed+3,m)+vn(ed-2,m))         !<  �����ع��ұ߽�                             
    end do

    do m=nst,ned
        do i=st+1,ed-1
            ve4(i,m) =          (1.0+2.0*r1)*( vn(i+1,m) + vn(i,m  ) )          -  &
                        3.0/2.0*(1.0+2.0*r1)*( vn(i+2,m) + vn(i-1,m) )          +  &
                        1.0/2.0*(1.0+2.0*r1)*( vn(i+3,m) + vn(i-2,m) )
        end do
    end do
    
    do i=st,ed
        mb(i)=1.0                          !< ���Խ�׷�ϵ�Bϵ��
    enddo
  
    do i=st+1,ed-1
        ma(i)=r1                           !< ���Խ�׷�ϵ�Aϵ�����ڵ㣩
    end do

    do i=st+1,ed-1                                                
        mc(i)=r1                           !< ���Խ�׷�ϵ�Cϵ�����ڵ㣩
    end do                                               


    
    do m=nst,ned
        do i=st,ed
            md(i) = ve4(i,m)               !< ���Խ�׷�ϵ��Ҳ�Dϵ��
        end do
        call tri_pursue(ma,mb,mc,md,ms,st,ed,1-ngn,ni+ngn) !< ���Խ�׷��
        do i=st,ed
            ve4(i,m) = ms(i)               !< �������Խ�׷�Ͻ�
        end do
    end do    
    
end subroutine vec4_via_T2N

subroutine scsh5pi2(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_bound_scheme
    implicit none
    integer(kind_int), intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  , intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(in)  :: sn(1-ngn:ni+ngn,1:3)
    integer(kind_int), intent(in)  :: nfs,nfe
    integer(kind_int), intent(in)  :: nge
    real(kind_real)  , intent(out) :: vl(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(out) :: vr(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , external    :: sub_limit
    integer(kind_int), intent(in)  :: chkid(nst:ned)
    real(kind_real),   intent(in)  :: chklim(nst:ned)
    integer(kind_int), intent(out) :: ndec  
    real(kind_real)   :: vln(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: vrn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: ve_lil(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: ve_lir(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: d2_li(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: d4_li(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: d5_li(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: g(1:3,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: wl(1:3,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: wr(1:3,1-ngn:ni+ngn,nst:ned)    
    real(kind_real)   :: aa,aa_scs
    real(kind_real)   :: sigma(1-ngn:ni+ngn),imax,imin,indicator
    real(kind_real), parameter :: cl(3) = (/1.D0/4.D0,3.D0/4.D0,0.D0/4.D0/)
    real(kind_real), parameter :: cr(3) = (/0.D0/4.D0,3.D0/4.D0,1.D0/4.D0/)      
    real(kind_real), parameter :: eps        = 1.0e-30
    real(kind_real), parameter :: indicatorc = 0.1
    real(kind_real)            :: is,eis2(3),bs,tau3,bl(3),br(3)
    integer(kind_int)          :: i,m,n,j,st,ed,stc,edc
    integer(kind_int)          :: id,YN,ind(1:ni+2*ngn),ind_st,ind_ed,pre_ed
    integer(kind_int)          :: nsec,sec_st(1:ni+2*ngn),sec_ed(1:ni+2*ngn)
    integer(kind_int)          :: nsecn,sec_stn(1:ni+2*ngn),sec_edn(1:ni+2*ngn)
    
    if (ni < 8) then
        call muscl2pv(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)
        return
    end if
    
    ndec = 0
    
    aa     = 0.01
    aa_scs = 0.001
    
    wl(:,:,:) = 0.5
    wr(:,:,:) = 0.5
    vln(:,:)    = 0.0
    ve_lil(:,:) = 0.0
    vrn(:,:)    = 0.0
    ve_lir(:,:) = 0.0    
           
    if (nfs == nbc_bound_scheme) then
        
        st  = 4
        stc = 0

        do m=nst,ned
            vl(0,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                  ! third-order
            vl(1,m) = (  3.0*vn(1,m) +  6.0*vn(2,m) -      vn(3,m) ) / 8.0   ! third-order 
            vl(2,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 + &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )
            vr(0,m) = vl(0,m)
            vr(1,m) = vl(1,m)
            vr(2,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 - &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )
            vl(3,m) =   ( 150.0*( vn(4,m) + vn(3,m) )          -  &
                           25.0*( vn(5,m) + vn(2,m) )          +  &
                            3.0*( vn(6,m) + vn(1,m) ) ) /256.0 -  &
                      aa*( 10.0*( vn(4,m) - vn(3,m) )          -  &
                            5.0*( vn(5,m) - vn(2,m) )          +  &
                                ( vn(6,m) - vn(1,m) ) ) 
            
            vr(3,m) =   ( 150.0*( vn(4,m) + vn(3,m) )          -  &
                           25.0*( vn(5,m) + vn(2,m) )          +  &
                            3.0*( vn(6,m) + vn(1,m) ) ) /256.0 +  &
                      aa*( 10.0*( vn(4,m) - vn(3,m) )          -  &
                            5.0*( vn(5,m) - vn(2,m) )          +  &
                                ( vn(6,m) - vn(1,m) ) )          
        end do
        
    else
        
        st  = -2
        stc = st
     
    end if
     
    if (nfe == nbc_bound_scheme) then
        
        ed  = ni - 4
        edc = ni

        do m=nst,ned
            vr(ni  ,m) = ( 15.0*vn(ni,m) - 10.0*vn(ni-1,m) +  3.0*vn(ni-2,m) ) / 8.0                ! third-order
            vl(ni-1,m) = (  3.0*vn(ni,m) +  6.0*vn(ni-1,m) -      vn(ni-2,m) ) / 8.0   ! third-order
            vl(ni-2,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 + &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) )
            vl(ni  ,m) = vr(ni  ,m)
            vr(ni-1,m) = vl(ni-1,m)
            vr(ni-2,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 - &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) ) 
            vl(ni-3,m) =   ( 150.0*( vn(ni-2,m) + vn(ni-3,m) )          -  &
                              25.0*( vn(ni-1,m) + vn(ni-4,m) )          +  &
                               3.0*( vn(ni  ,m) + vn(ni-5,m) ) ) /256.0 -  &
                         aa*( 10.0*( vn(ni-2,m) - vn(ni-3,m) )          -  &
                               5.0*( vn(ni-1,m) - vn(ni-4,m) )          +  &
                                   ( vn(ni  ,m) - vn(ni-5,m) ) ) 
            
            vr(ni-3,m) =   ( 150.0*( vn(ni-2,m) + vn(ni-3,m) )          -  &
                              25.0*( vn(ni-1,m) + vn(ni-4,m) )          +  &
                               3.0*( vn(ni  ,m) + vn(ni-5,m) ) ) /256.0 +  &
                         aa*( 10.0*( vn(ni-2,m) - vn(ni-3,m) )          -  &
                               5.0*( vn(ni-1,m) - vn(ni-4,m) )          +  &
                                   ( vn(ni  ,m) - vn(ni-5,m) ) )              
        end do
        
    else
        
        ed  = ni + 2
        edc = ed
            
    end if
    
    do m=nst,ned
        do i=st-2,ed+2
            g(1,i,m) = vn(i  ,m) - vn(i-1,m)            
            g(2,i,m) = vn(i+1,m) - vn(i  ,m)
            g(3,i,m) = vn(i+2,m) - vn(i+1,m)            
        end do
    end do   
    
    do m=nst,ned
        do i=st-2,ed+2
            do n=1,3
                is = g(n,i,m)*g(n,i,m)
                eis2(n) = is**2
            end do

            tau3 = abs(eis2(2)-eis2(1))
            bl(1) = 1.0+tau3/(eis2(1)+eps)
            bl(2) = 1.0+tau3/(eis2(2)+eps)
            bs = bl(1) + bl(2)
            wl(2,i,m) = bl(2)/bs
            
            tau3 = abs(eis2(3)-eis2(2))
            br(2) = 1.0+tau3/(eis2(2)+eps)
            br(3) = 1.0+tau3/(eis2(3)+eps)
            bs = br(2) + br(3)
            wr(2,i,m) = br(2)/bs
        end do
    end do
    
    do i=st,ed
        imin=1.0
        imax=0.0
        do m=nst,ned
            imin         = min(wl(2,i-3,m),wl(2,i-2,m),wl(2,i-1,m),wl(2,i,m),wl(2,i+1,m),wl(2,i+2,m),wl(2,i+3,m), &
                               wr(2,i-3,m),wr(2,i-2,m),wr(2,i-1,m),wr(2,i,m),wr(2,i+1,m),wr(2,i+2,m),wr(2,i+3,m),imin)                                
            imax         = max(wl(2,i-3,m),wl(2,i-2,m),wl(2,i-1,m),wl(2,i,m),wl(2,i+1,m),wl(2,i+2,m),wl(2,i+3,m), &
                               wr(2,i-3,m),wr(2,i-2,m),wr(2,i-1,m),wr(2,i,m),wr(2,i+1,m),wr(2,i+2,m),wr(2,i+3,m),imax)
        end do
        indicator    = imin/imax
        if(indicator>indicatorc) then
            sigma(i)=1.0
        else
            sigma(i)=0.0
        endif        
    end do
    
    id=0
    do i=st,ed
        YN=0
        if(sigma(i)<0.5) then
            YN=YN+1
        end if
        if(YN>0) then
            id=id+1
            ind(id)=i           
        end if
    end do
       
    do i=1-ngn,ni+ngn
        sigma(i)=0.0
    end do
    if(ed-st<9) then
        ind_st=st
        ind_ed=ed
        nsec=0
        
        nsecn=1
        sec_stn(1)=st
        sec_edn(1)=ed
    else
        if(id==0) then
            nsec=1
            sec_st(1)=st
            sec_ed(1)=ed
            do n=sec_st(1),sec_ed(1)
                sigma(n)=1.0
            end do
            ind_st=st
            ind_ed=ed
            nsecn=0
        else
            ind_st=st
            ind_ed=ed
            nsec=0
            nsecn=0
            if(ind(1)-ind_st>9) then
                nsec=1
                sec_st(1)=ind_st
                sec_ed(1)=ind(1)-1
                pre_ed   =sec_ed(1)
                do n=sec_st(1),sec_ed(1)
                    sigma(n)=1.0
                end do
            else
                nsecn=1
                sec_stn(1)=ind_st
                sec_edn(1)=ind(1)
                pre_ed    =sec_edn(1)
            end if
            do i=2,id
                if(ind(i)-ind(i-1)>9) then
                    nsec=nsec+1
                    sec_st(nsec)=ind(i-1)+1
                    sec_ed(nsec)=ind(i)-1
                    do n=sec_st(nsec),sec_ed(nsec)
                        sigma(n)=1.0
                    end do
                    nsecn=nsecn+1
                    sec_stn(nsecn)=pre_ed+1
                    sec_edn(nsecn)=ind(i-1)
                    pre_ed        =sec_ed(nsec)
                else
                    if(ind(i)-ind(i-1)>1 .or. i==id) then
                        nsecn=nsecn+1
                        sec_stn(nsecn)=pre_ed+1
                        sec_edn(nsecn)=ind(i)
                        pre_ed        =sec_edn(nsecn)
                    end if
                end if
            end do
            if(ind_ed-ind(id)>9) then
                nsec=nsec+1
                sec_st(nsec)=pre_ed+1
                sec_ed(nsec)=ind_ed
                do n=sec_st(nsec),sec_ed(nsec)
                    sigma(n)=1.0
                end do
            else
                nsecn=nsecn+1
                sec_stn(nsecn)=pre_ed+1
                sec_edn(nsecn)=ind_ed
            end if
        end if
    end if
    
    do m=nst,ned
        do j=1,nsecn
            do i=sec_stn(j),sec_edn(j)
                g(1,i,m) = vn(i  ,m) - vn(i-1,m) 
                g(2,i,m) = vn(i+1,m) - vn(i  ,m)
                g(3,i,m) = vn(i+2,m) - vn(i+1,m)
            end do
        end do      
    end do            
            
    do m=nst,ned
        do j=1,nsecn
            do i=sec_stn(j),sec_edn(j)
                do n=1,3
                    is = g(n,i,m)*g(n,i,m)
                    eis2(n) = is**2
                end do
                
                tau3 = abs(eis2(2)-eis2(1))
                
                do n=1,3
                    bl(n) = cl(n)*(1.0+tau3/(eis2(n)+eps))
                end do
                
                bs = bl(1) + bl(2)
                do n=1,2
                    wl(n,i,m) = bl(n)/bs
                end do
          
                tau3 = abs(eis2(3)-eis2(2))           
                do n=2,3
                    br(n) = cr(n)*(1.0+tau3/(eis2(n)+eps))
                end do

                bs = br(2) + br(3)
                do n=2,3
                    wr(n,i,m) = br(n)/bs
                end do
            end do
        end do
    end do
    
!!    open(2,file='sigma.dat',form='formatted')
!!    write(2,*) 'variables="i" "sigmal" "sigmar"'
!!    do i=st,ed
!!
!!        write(2,200) i,sigmal(i,1),sigmar(i,1)
!!
!!    end do
!!
!!200 format(I4,1x,e12.5,1x,e12.5)
!!    close(2)
!!
!!    stop    
    
    do m=nst,ned
        do j=1,nsecn
            do i=sec_stn(j),sec_edn(j)
                vln(i,m) = vn(i  ,m) + ( wl(1,i,m)*( g(1,i,m)) + &
                                         wl(2,i,m)*( g(2,i,m)) )/2.0

                vrn(i,m) = vn(i+1,m) + ( wr(2,i,m)*(-g(2,i,m)) + &
                                         wr(3,i,m)*(-g(3,i,m)) )/2.0
            end do
        end do
    end do  
      
    do j=1,nsec
        call d5_via_T2N(ni,nst,ned,ngn,vn,nfs,nfe,nge,d5_li,sec_st(j),sec_ed(j))
        
        call vec2_via_T4N(ni,nst,ned,ngn,vn,nfs,nfe,nge,d2_li,sec_st(j),sec_ed(j))
        
        call vec4_via_T2N(ni,nst,ned,ngn,vn,nfs,nfe,nge,d4_li,sec_st(j),sec_ed(j))
    end do
    
    do m=nst,ned
        do j=1,nsec
            do i=sec_st(j),sec_ed(j)
                
                ve_lil(i,m) = 0.5*(vn(i+1,m) + vn(i  ,m)) - 1.0/8.0*d2_li(i,m) - 1.0/384.0*d4_li(i,m) -  aa_scs*d5_li(i,m)
                            
                ve_lir(i,m) = 0.5*(vn(i+1,m) + vn(i  ,m)) - 1.0/8.0*d2_li(i,m) - 1.0/384.0*d4_li(i,m) +  aa_scs*d5_li(i,m)

            end do
        end do
    end do
        
    do m=nst,ned
        do i=st,ed
            
               vl(i,m) = sigma(i)*ve_lil(i,m) + vln(i,m)*(1.0-sigma(i))
                         
               vr(i,m) = sigma(i)*ve_lir(i,m) + vrn(i,m)*(1.0-sigma(i))            

        end do
    end do  

    do i=stc,edc
        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vl(i,m) < chklim(m) ) then
                    vl(i,m) = vn(i  ,m)
                    ndec = ndec + 1
                end if

                if (vr(i,m) < chklim(m) ) then
                    vr(i,m) = vn(i+1,m)
                    ndec = ndec + 1
                end if
           end if
        end do
    end do

end subroutine scsh5pi2

subroutine scsh5ci2(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_bound_scheme
    implicit none
    integer(kind_int), intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  , intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(in)  :: sn(1-ngn:ni+ngn,1:3)
    integer(kind_int), intent(in)  :: nfs,nfe
    integer(kind_int), intent(in)  :: nge
    real(kind_real)  , intent(out) :: vl(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(out) :: vr(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , external    :: sub_limit
    integer(kind_int), intent(in)  :: chkid(nst:ned)
    real(kind_real),   intent(in)  :: chklim(nst:ned)
    integer(kind_int), intent(out) :: ndec  
    real(kind_real)   :: vln(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: vrn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: ve_lil(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: ve_lir(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: d2_li(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: d4_li(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: d5_li(1-ngn:ni+ngn,nst:ned)  
    real(kind_real)   :: nt,nx,ny,nz,fp(nst:ned)
    real(kind_real)   :: vp(nst:ned),vc(nst:ned)
    real(kind_real)   :: q(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: q1(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: g(1:3,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: wl(1:3,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: wr(1:3,1-ngn:ni+ngn,nst:ned)    
    real(kind_real)   :: aa,aa_scs
    real(kind_real)   :: sigma(1-ngn:ni+ngn),imax,imin,indicator
    real(kind_real), parameter :: cl(3) = (/1.D0/4.D0,3.D0/4.D0,0.D0/4.D0/)
    real(kind_real), parameter :: cr(3) = (/0.D0/4.D0,3.D0/4.D0,1.D0/4.D0/)      
    real(kind_real), parameter :: eps        = 1.0e-30
    real(kind_real), parameter :: indicatorc = 0.1
    real(kind_real)            :: is,eis2(3),bs,tau3,bl(3),br(3)
    integer(kind_int)          :: i,m,n,j,st,ed,stc,edc
    integer(kind_int)          :: id,YN,ind(1:ni+2*ngn),ind_st,ind_ed,pre_ed
    integer(kind_int)          :: nsec,sec_st(1:ni+2*ngn),sec_ed(1:ni+2*ngn)
    integer(kind_int)          :: nsecn,sec_stn(1:ni+2*ngn),sec_edn(1:ni+2*ngn)
    
    if (ni < 8) then
        call muscl2pv(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)
        return
    end if
    
    ndec = 0
    
    aa     = 0.01
    aa_scs = 0.001
    
    wl(:,:,:) = 0.5
    wr(:,:,:) = 0.5
    vln(:,:)    = 0.0
    ve_lil(:,:) = 0.0
    vrn(:,:)    = 0.0
    ve_lir(:,:) = 0.0    
           
    if (nfs == nbc_bound_scheme) then
        
        st  = 4
        stc = 0

        do m=nst,ned
            vl(0,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                  ! third-order
            vl(1,m) = (  3.0*vn(1,m) +  6.0*vn(2,m) -      vn(3,m) ) / 8.0   ! third-order 
            vl(2,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 + &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )
            vr(0,m) = vl(0,m)
            vr(1,m) = vl(1,m)
            vr(2,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 - &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )
            vl(3,m) =   ( 150.0*( vn(4,m) + vn(3,m) )          -  &
                           25.0*( vn(5,m) + vn(2,m) )          +  &
                            3.0*( vn(6,m) + vn(1,m) ) ) /256.0 -  &
                      aa*( 10.0*( vn(4,m) - vn(3,m) )          -  &
                            5.0*( vn(5,m) - vn(2,m) )          +  &
                                ( vn(6,m) - vn(1,m) ) ) 
            
            vr(3,m) =   ( 150.0*( vn(4,m) + vn(3,m) )          -  &
                           25.0*( vn(5,m) + vn(2,m) )          +  &
                            3.0*( vn(6,m) + vn(1,m) ) ) /256.0 +  &
                      aa*( 10.0*( vn(4,m) - vn(3,m) )          -  &
                            5.0*( vn(5,m) - vn(2,m) )          +  &
                                ( vn(6,m) - vn(1,m) ) )          
        end do
        
    else
        
        st  = -2
        stc = st
     
    end if
     
    if (nfe == nbc_bound_scheme) then
        
        ed  = ni - 4
        edc = ni

        do m=nst,ned
            vr(ni  ,m) = ( 15.0*vn(ni,m) - 10.0*vn(ni-1,m) +  3.0*vn(ni-2,m) ) / 8.0                ! third-order
            vl(ni-1,m) = (  3.0*vn(ni,m) +  6.0*vn(ni-1,m) -      vn(ni-2,m) ) / 8.0   ! third-order
            vl(ni-2,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 + &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) )
            vl(ni  ,m) = vr(ni  ,m)
            vr(ni-1,m) = vl(ni-1,m)
            vr(ni-2,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 - &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) ) 
            vl(ni-3,m) =   ( 150.0*( vn(ni-2,m) + vn(ni-3,m) )          -  &
                              25.0*( vn(ni-1,m) + vn(ni-4,m) )          +  &
                               3.0*( vn(ni  ,m) + vn(ni-5,m) ) ) /256.0 -  &
                         aa*( 10.0*( vn(ni-2,m) - vn(ni-3,m) )          -  &
                               5.0*( vn(ni-1,m) - vn(ni-4,m) )          +  &
                                   ( vn(ni  ,m) - vn(ni-5,m) ) ) 
            
            vr(ni-3,m) =   ( 150.0*( vn(ni-2,m) + vn(ni-3,m) )          -  &
                              25.0*( vn(ni-1,m) + vn(ni-4,m) )          +  &
                               3.0*( vn(ni  ,m) + vn(ni-5,m) ) ) /256.0 +  &
                         aa*( 10.0*( vn(ni-2,m) - vn(ni-3,m) )          -  &
                               5.0*( vn(ni-1,m) - vn(ni-4,m) )          +  &
                                   ( vn(ni  ,m) - vn(ni-5,m) ) )              
        end do
        
    else
        
        ed  = ni + 2
        edc = ed
            
    end if
    
    do m=nst,ned
        do i=st-2,ed+2
            g(1,i,m) = vn(i  ,m) - vn(i-1,m)            
            g(2,i,m) = vn(i+1,m) - vn(i  ,m)
            g(3,i,m) = vn(i+2,m) - vn(i+1,m)            
        end do
    end do   
    
    do m=nst,ned
        do i=st-2,ed+2
            do n=1,3
                is = g(n,i,m)*g(n,i,m)
                eis2(n) = is**2
            end do

            tau3 = abs(eis2(2)-eis2(1))
            bl(1) = 1.0+tau3/(eis2(1)+eps)
            bl(2) = 1.0+tau3/(eis2(2)+eps)
            bs = bl(1) + bl(2)
            wl(2,i,m) = bl(2)/bs
            
            tau3 = abs(eis2(3)-eis2(2))
            br(2) = 1.0+tau3/(eis2(2)+eps)
            br(3) = 1.0+tau3/(eis2(3)+eps)
            bs = br(2) + br(3)
            wr(2,i,m) = br(2)/bs
        end do
    end do
    
    do i=st,ed
        imin=1.0
        imax=0.0
        do m=nst,ned
            imin         = min(wl(2,i-3,m),wl(2,i-2,m),wl(2,i-1,m),wl(2,i,m),wl(2,i+1,m),wl(2,i+2,m),wl(2,i+3,m), &
                               wr(2,i-3,m),wr(2,i-2,m),wr(2,i-1,m),wr(2,i,m),wr(2,i+1,m),wr(2,i+2,m),wr(2,i+3,m),imin)                                
            imax         = max(wl(2,i-3,m),wl(2,i-2,m),wl(2,i-1,m),wl(2,i,m),wl(2,i+1,m),wl(2,i+2,m),wl(2,i+3,m), &
                               wr(2,i-3,m),wr(2,i-2,m),wr(2,i-1,m),wr(2,i,m),wr(2,i+1,m),wr(2,i+2,m),wr(2,i+3,m),imax)
        end do
        indicator    = imin/imax
        if(indicator>indicatorc) then
            sigma(i)=1.0
        else
            sigma(i)=0.0
        endif        
    end do
    
    id=0
    do i=st,ed
        YN=0
        if(sigma(i)<0.5) then
            YN=YN+1
        end if
        if(YN>0) then
            id=id+1
            ind(id)=i           
        end if
    end do
       
    do i=1-ngn,ni+ngn
        sigma(i)=0.0
    end do
    if(ed-st<9) then
        ind_st=st
        ind_ed=ed
        nsec=0
        
        nsecn=1
        sec_stn(1)=st
        sec_edn(1)=ed
    else
        if(id==0) then
            nsec=1
            sec_st(1)=st
            sec_ed(1)=ed
            do n=sec_st(1),sec_ed(1)
                sigma(n)=1.0
            end do
            ind_st=st
            ind_ed=ed
            nsecn=0
        else
            ind_st=st
            ind_ed=ed
            nsec=0
            nsecn=0
            if(ind(1)-ind_st>9) then
                nsec=1
                sec_st(1)=ind_st
                sec_ed(1)=ind(1)-1
                pre_ed   =sec_ed(1)
                do n=sec_st(1),sec_ed(1)
                    sigma(n)=1.0
                end do
            else
                nsecn=1
                sec_stn(1)=ind_st
                sec_edn(1)=ind(1)
                pre_ed    =sec_edn(1)
            end if
            do i=2,id
                if(ind(i)-ind(i-1)>9) then
                    nsec=nsec+1
                    sec_st(nsec)=ind(i-1)+1
                    sec_ed(nsec)=ind(i)-1
                    do n=sec_st(nsec),sec_ed(nsec)
                        sigma(n)=1.0
                    end do
                    nsecn=nsecn+1
                    sec_stn(nsecn)=pre_ed+1
                    sec_edn(nsecn)=ind(i-1)
                    pre_ed        =sec_ed(nsec)
                else
                    if(ind(i)-ind(i-1)>1 .or. i==id) then
                        nsecn=nsecn+1
                        sec_stn(nsecn)=pre_ed+1
                        sec_edn(nsecn)=ind(i)
                        pre_ed        =sec_edn(nsecn)
                    end if
                end if
            end do
            if(ind_ed-ind(id)>9) then
                nsec=nsec+1
                sec_st(nsec)=pre_ed+1
                sec_ed(nsec)=ind_ed
                do n=sec_st(nsec),sec_ed(nsec)
                    sigma(n)=1.0
                end do
            else
                nsecn=nsecn+1
                sec_stn(nsecn)=pre_ed+1
                sec_edn(nsecn)=ind_ed
            end if
        end if
    end if
    
    do m=nst,ned
        do j=1,nsecn
            do i=sec_stn(j),sec_edn(j)
                g(1,i,m) = vn(i  ,m) - vn(i-1,m) 
                g(2,i,m) = vn(i+1,m) - vn(i  ,m)
                g(3,i,m) = vn(i+2,m) - vn(i+1,m)
            end do
        end do      
    end do
    
    do j=1,nsecn
        do i=sec_stn(j),sec_edn(j)
            nx = (sn(i,1) + sn(i+1,1))/2.0
            ny = (sn(i,2) + sn(i+1,2))/2.0
            nz = (sn(i,3) + sn(i+1,3))/2.0

            do m=nst,ned
                vp(m) = (vn(i,m) + vn(i+1,m))/2.0
                fp(m) = vn(i,m)
            end do
            call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
            do m=nst,ned
                q(i,m) = vc(m)
            end do
            
            do m=nst,ned
                fp(m) = vn(i+1,m)
            end do
            call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
            do m=nst,ned
                q1(i+1,m) = vc(m)
            end do        
             
            do n=1,3
                do m=nst,ned
                    fp(m) = g(n,i,m)
                end do
                call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
                do m=nst,ned
                    g(n,i,m) = vc(m)
                end do        
            end do
        end do
    end do   
            
            
    do m=nst,ned
        do j=1,nsecn
            do i=sec_stn(j),sec_edn(j)
                do n=1,3
                    is = g(n,i,m)*g(n,i,m)
                    eis2(n) = is**2
                end do
                
                tau3 = abs(eis2(2)-eis2(1))
                
                do n=1,3
                    bl(n) = cl(n)*(1.0+tau3/(eis2(n)+eps))
                end do
                
                bs = bl(1) + bl(2)
                do n=1,2
                    wl(n,i,m) = bl(n)/bs
                end do
          
                tau3 = abs(eis2(3)-eis2(2))           
                do n=2,3
                    br(n) = cr(n)*(1.0+tau3/(eis2(n)+eps))
                end do

                bs = br(2) + br(3)
                do n=2,3
                    wr(n,i,m) = br(n)/bs
                end do
            end do
        end do
    end do
    
!!    open(2,file='sigma.dat',form='formatted')
!!    write(2,*) 'variables="i" "sigmal" "sigmar"'
!!    do i=st,ed
!!
!!        write(2,200) i,sigmal(i,1),sigmar(i,1)
!!
!!    end do
!!
!!200 format(I4,1x,e12.5,1x,e12.5)
!!    close(2)
!!
!!    stop    
    
    do m=nst,ned
        do j=1,nsecn
            do i=sec_stn(j),sec_edn(j)
                vln(i,m) =  q(i  ,m) + ( wl(1,i,m)*( g(1,i,m)) + &
                                         wl(2,i,m)*( g(2,i,m)) )/2.0

                vrn(i,m) = q1(i+1,m) + ( wr(2,i,m)*(-g(2,i,m)) + &
                                         wr(3,i,m)*(-g(3,i,m)) )/2.0
            end do
        end do
    end do
    
    do j=1,nsecn
        do i=sec_stn(j),sec_edn(j)
            nx = (sn(i,1) + sn(i+1,1))/2.0
            ny = (sn(i,2) + sn(i+1,2))/2.0
            nz = (sn(i,3) + sn(i+1,3))/2.0

            do m=nst,ned
                vp(m) = (vn(i,m) + vn(i+1,m))/2.0
                vc(m) = vln(i,m)
            end do

            call RxQ_p(nst,ned,vp,nt,nx,ny,nz,vc,fp) 
            do m=nst,ned
                vln(i,m) = fp(m)
            end do
            
            do m=nst,ned
                vc(m) = vrn(i ,m)
            end do

            call RxQ_p(nst,ned,vp,nt,nx,ny,nz,vc,fp) 
            do m=nst,ned
                vrn(i,m) = fp(m)
            end do
        end do
    end do     
      
    do j=1,nsec
        call d5_via_T2N(ni,nst,ned,ngn,vn,nfs,nfe,nge,d5_li,sec_st(j),sec_ed(j))
        
        call vec2_via_T4N(ni,nst,ned,ngn,vn,nfs,nfe,nge,d2_li,sec_st(j),sec_ed(j))
        
        call vec4_via_T2N(ni,nst,ned,ngn,vn,nfs,nfe,nge,d4_li,sec_st(j),sec_ed(j))
    end do
    
    do m=nst,ned
        do j=1,nsec
            do i=sec_st(j),sec_ed(j)
                
                ve_lil(i,m) = 0.5*(vn(i+1,m) + vn(i  ,m)) - 1.0/8.0*d2_li(i,m) - 1.0/384.0*d4_li(i,m) -  aa_scs*d5_li(i,m)
                            
                ve_lir(i,m) = 0.5*(vn(i+1,m) + vn(i  ,m)) - 1.0/8.0*d2_li(i,m) - 1.0/384.0*d4_li(i,m) +  aa_scs*d5_li(i,m)

            end do
        end do
    end do
        
    do m=nst,ned
        do i=st,ed
            
               vl(i,m) = sigma(i)*ve_lil(i,m) + vln(i,m)*(1.0-sigma(i))
                         
               vr(i,m) = sigma(i)*ve_lir(i,m) + vrn(i,m)*(1.0-sigma(i))            

        end do
    end do  

    do i=stc,edc
        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vl(i,m) < chklim(m) ) then
                    vl(i,m) = vn(i  ,m)
                    ndec = ndec + 1
                end if

                if (vr(i,m) < chklim(m) ) then
                    vr(i,m) = vn(i+1,m)
                    ndec = ndec + 1
                end if
           end if
        end do
    end do

end subroutine scsh5ci2

subroutine scsh5pi3(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_bound_scheme
    implicit none
    integer(kind_int), intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  , intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(in)  :: sn(1-ngn:ni+ngn,1:3)
    integer(kind_int), intent(in)  :: nfs,nfe
    integer(kind_int), intent(in)  :: nge
    real(kind_real)  , intent(out) :: vl(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(out) :: vr(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , external    :: sub_limit
    integer(kind_int), intent(in)  :: chkid(nst:ned)
    real(kind_real),   intent(in)  :: chklim(nst:ned)
    integer(kind_int), intent(out) :: ndec  
    real(kind_real)   :: vln(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: vrn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: ve_lil(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: ve_lir(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: d2_li(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: d4_li(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: d5_li(1-ngn:ni+ngn,nst:ned)   
    real(kind_real)   :: s(1:4,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: g(1:4,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: wl(1:4,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: wr(1:4,1-ngn:ni+ngn,nst:ned)    
    real(kind_real)   :: aa,aa_scs
    real(kind_real)   :: sigma(1-ngn:ni+ngn),imax,imin,indicator
    real(kind_real), parameter :: cl(4) = (/1.D0/16.D0,10.D0/16.D0,5.D0/16.D0,0.D0/16.D0/)
    real(kind_real), parameter :: cr(4) = (/0.D0/16.D0,5.D0/16.D0,10.D0/16.D0,1.D0/16.D0/)      
    real(kind_real), parameter :: eps        = 1.0e-30
    real(kind_real), parameter :: indicatorc = 0.1
    real(kind_real)            :: is,eis2(4),bs,tau5,bl(4),br(4)
    integer(kind_int)          :: i,m,n,j,st,ed,stc,edc
    integer(kind_int)          :: id,YN,ind(1:ni+2*ngn),ind_st,ind_ed,pre_ed
    integer(kind_int)          :: nsec,sec_st(1:ni+2*ngn),sec_ed(1:ni+2*ngn)
    integer(kind_int)          :: nsecn,sec_stn(1:ni+2*ngn),sec_edn(1:ni+2*ngn)
    
    if (ni < 8) then
        call muscl2pv(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)
        return
    end if
    
    ndec = 0
    
    aa     = 0.01
    aa_scs = 0.001
    
    wl(:,:,:) = 0.3333
    wr(:,:,:) = 0.3333
    vln(:,:)    = 0.0
    ve_lil(:,:) = 0.0
    vrn(:,:)    = 0.0
    ve_lir(:,:) = 0.0    
           
    if (nfs == nbc_bound_scheme) then
        
        st  = 4
        stc = 0

        do m=nst,ned
            vl(0,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                  ! third-order
            vl(1,m) = (  3.0*vn(1,m) +  6.0*vn(2,m) -      vn(3,m) ) / 8.0   ! third-order 
            vl(2,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 + &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )
            vr(0,m) = vl(0,m)
            vr(1,m) = vl(1,m)
            vr(2,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 - &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )
            vl(3,m) =   ( 150.0*( vn(4,m) + vn(3,m) )          -  &
                           25.0*( vn(5,m) + vn(2,m) )          +  &
                            3.0*( vn(6,m) + vn(1,m) ) ) /256.0 -  &
                      aa*( 10.0*( vn(4,m) - vn(3,m) )          -  &
                            5.0*( vn(5,m) - vn(2,m) )          +  &
                                ( vn(6,m) - vn(1,m) ) ) 
            
            vr(3,m) =   ( 150.0*( vn(4,m) + vn(3,m) )          -  &
                           25.0*( vn(5,m) + vn(2,m) )          +  &
                            3.0*( vn(6,m) + vn(1,m) ) ) /256.0 +  &
                      aa*( 10.0*( vn(4,m) - vn(3,m) )          -  &
                            5.0*( vn(5,m) - vn(2,m) )          +  &
                                ( vn(6,m) - vn(1,m) ) )          
        end do
        
    else
        
        st  = -2
        stc = st
     
    end if
     
    if (nfe == nbc_bound_scheme) then
        
        ed  = ni - 4
        edc = ni

        do m=nst,ned
            vr(ni  ,m) = ( 15.0*vn(ni,m) - 10.0*vn(ni-1,m) +  3.0*vn(ni-2,m) ) / 8.0                ! third-order
            vl(ni-1,m) = (  3.0*vn(ni,m) +  6.0*vn(ni-1,m) -      vn(ni-2,m) ) / 8.0   ! third-order
            vl(ni-2,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 + &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) )
            vl(ni  ,m) = vr(ni  ,m)
            vr(ni-1,m) = vl(ni-1,m)
            vr(ni-2,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 - &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) ) 
            vl(ni-3,m) =   ( 150.0*( vn(ni-2,m) + vn(ni-3,m) )          -  &
                              25.0*( vn(ni-1,m) + vn(ni-4,m) )          +  &
                               3.0*( vn(ni  ,m) + vn(ni-5,m) ) ) /256.0 -  &
                         aa*( 10.0*( vn(ni-2,m) - vn(ni-3,m) )          -  &
                               5.0*( vn(ni-1,m) - vn(ni-4,m) )          +  &
                                   ( vn(ni  ,m) - vn(ni-5,m) ) ) 
            
            vr(ni-3,m) =   ( 150.0*( vn(ni-2,m) + vn(ni-3,m) )          -  &
                              25.0*( vn(ni-1,m) + vn(ni-4,m) )          +  &
                               3.0*( vn(ni  ,m) + vn(ni-5,m) ) ) /256.0 +  &
                         aa*( 10.0*( vn(ni-2,m) - vn(ni-3,m) )          -  &
                               5.0*( vn(ni-1,m) - vn(ni-4,m) )          +  &
                                   ( vn(ni  ,m) - vn(ni-5,m) ) )              
        end do
        
    else
        
        ed  = ni + 2
        edc = ed
            
    end if
    
    do m=nst,ned
        do i=st-2,ed+2
            g(1,i,m) = vn(i  ,m) - vn(i-1,m)            
            g(2,i,m) = vn(i+1,m) - vn(i  ,m)
            g(3,i,m) = vn(i+2,m) - vn(i+1,m)            
        end do
    end do   
    
    do m=nst,ned
        do i=st-2,ed+2
            do n=1,3
                is = g(n,i,m)*g(n,i,m)
                eis2(n) = is**2
            end do

            tau5 = abs(eis2(2)-eis2(1))
            bl(1) = 1.0+tau5/(eis2(1)+eps)
            bl(2) = 1.0+tau5/(eis2(2)+eps)
            bs = bl(1) + bl(2)
            wl(2,i,m) = bl(2)/bs
            
            tau5 = abs(eis2(3)-eis2(2))
            br(2) = 1.0+tau5/(eis2(2)+eps)
            br(3) = 1.0+tau5/(eis2(3)+eps)
            bs = br(2) + br(3)
            wr(2,i,m) = br(2)/bs
        end do
    end do
    
    do i=st,ed
        imin=1.0
        imax=0.0
        do m=nst,ned
            imin         = min(wl(2,i-3,m),wl(2,i-2,m),wl(2,i-1,m),wl(2,i,m),wl(2,i+1,m),wl(2,i+2,m),wl(2,i+3,m), &
                               wr(2,i-3,m),wr(2,i-2,m),wr(2,i-1,m),wr(2,i,m),wr(2,i+1,m),wr(2,i+2,m),wr(2,i+3,m),imin)                                
            imax         = max(wl(2,i-3,m),wl(2,i-2,m),wl(2,i-1,m),wl(2,i,m),wl(2,i+1,m),wl(2,i+2,m),wl(2,i+3,m), &
                               wr(2,i-3,m),wr(2,i-2,m),wr(2,i-1,m),wr(2,i,m),wr(2,i+1,m),wr(2,i+2,m),wr(2,i+3,m),imax)
        end do
        indicator    = imin/imax
        if(indicator>indicatorc) then
            sigma(i)=1.0
        else
            sigma(i)=0.0
        endif        
    end do
    
    id=0
    do i=st,ed
        YN=0
        if(sigma(i)<0.5) then
            YN=YN+1
        end if
        if(YN>0) then
            id=id+1
            ind(id)=i           
        end if
    end do
       
    do i=1-ngn,ni+ngn
        sigma(i)=0.0
    end do
    if(ed-st<9) then
        ind_st=st
        ind_ed=ed
        nsec=0
        
        nsecn=1
        sec_stn(1)=st
        sec_edn(1)=ed
    else
        if(id==0) then
            nsec=1
            sec_st(1)=st
            sec_ed(1)=ed
            do n=sec_st(1),sec_ed(1)
                sigma(n)=1.0
            end do
            ind_st=st
            ind_ed=ed
            nsecn=0
        else
            ind_st=st
            ind_ed=ed
            nsec=0
            nsecn=0
            if(ind(1)-ind_st>9) then
                nsec=1
                sec_st(1)=ind_st
                sec_ed(1)=ind(1)-1
                pre_ed   =sec_ed(1)
                do n=sec_st(1),sec_ed(1)
                    sigma(n)=1.0
                end do
            else
                nsecn=1
                sec_stn(1)=ind_st
                sec_edn(1)=ind(1)
                pre_ed    =sec_edn(1)
            end if
            do i=2,id
                if(ind(i)-ind(i-1)>9) then
                    nsec=nsec+1
                    sec_st(nsec)=ind(i-1)+1
                    sec_ed(nsec)=ind(i)-1
                    do n=sec_st(nsec),sec_ed(nsec)
                        sigma(n)=1.0
                    end do
                    nsecn=nsecn+1
                    sec_stn(nsecn)=pre_ed+1
                    sec_edn(nsecn)=ind(i-1)
                    pre_ed        =sec_ed(nsec)
                else
                    if(ind(i)-ind(i-1)>1 .or. i==id) then
                        nsecn=nsecn+1
                        sec_stn(nsecn)=pre_ed+1
                        sec_edn(nsecn)=ind(i)
                        pre_ed        =sec_edn(nsecn)
                    end if
                end if
            end do
            if(ind_ed-ind(id)>9) then
                nsec=nsec+1
                sec_st(nsec)=pre_ed+1
                sec_ed(nsec)=ind_ed
                do n=sec_st(nsec),sec_ed(nsec)
                    sigma(n)=1.0
                end do
            else
                nsecn=nsecn+1
                sec_stn(nsecn)=pre_ed+1
                sec_edn(nsecn)=ind_ed
            end if
        end if
    end if
    
    do m=nst,ned
        do j=1,nsecn
            do i=sec_stn(j),sec_edn(j)
                g(1,i,m) =      vn(i-2,m) - 3.0*vn(i-1,m) + 2.0*vn(i  ,m) 
                g(2,i,m) =     -vn(i  ,m) +     vn(i+1,m)
                g(3,i,m) =     -vn(i  ,m) +     vn(i+1,m)
                g(4,i,m) = -2.0*vn(i+1,m) + 3.0*vn(i+2,m) -     vn(i+3,m)

                s(1,i,m) = vn(i-2,m) - 2.0*vn(i-1,m) + vn(i,  m)
                s(2,i,m) = vn(i-1,m) - 2.0*vn(i  ,m) + vn(i+1,m) 
                s(3,i,m) = vn(i  ,m) - 2.0*vn(i+1,m) + vn(i+2,m)
                s(4,i,m) = vn(i+1,m) - 2.0*vn(i+2,m) + vn(i+3,m)
            end do
        end do      
    end do            
            
    do m=nst,ned
        do j=1,nsecn
            do i=sec_stn(j),sec_edn(j)
                do n=1,4
                    is = g(n,i,m)*g(n,i,m) + s(n,i,m)*s(n,i,m)
                    eis2(n) = is**2
                end do
                
                tau5 = abs(eis2(3)-eis2(1))
                
                do n=1,3
                    bl(n) = cl(n)*(1.0+tau5/(eis2(n)+eps))
                end do
                
                bs = bl(1) + bl(2) + bl(3)
                do n=1,3
                    wl(n,i,m) = bl(n)/bs
                end do
          
                tau5 = abs(eis2(4)-eis2(2))           
                do n=2,4
                    br(n) = cr(n)*(1.0+tau5/(eis2(n)+eps))
                end do

                bs = br(2) + br(3) + br(4)
                do n=2,4
                    wr(n,i,m) = br(n)/bs
                end do
            end do
        end do
    end do
    
!!    open(2,file='sigma.dat',form='formatted')
!!    write(2,*) 'variables="i" "sigmal" "sigmar"'
!!    do i=st,ed
!!
!!        write(2,200) i,sigmal(i,1),sigmar(i,1)
!!
!!    end do
!!
!!200 format(I4,1x,e12.5,1x,e12.5)
!!    close(2)
!!
!!    stop    
    
    do m=nst,ned
        do j=1,nsecn
            do i=sec_stn(j),sec_edn(j)
                vln(i,m) = vn(i  ,m) + ( wl(1,i,m)*( 4.0*g(1,i,m)-s(1,i,m)) + &
                                         wl(2,i,m)*( 4.0*g(2,i,m)-s(2,i,m)) + &
                                         wl(3,i,m)*( 4.0*g(3,i,m)-s(3,i,m)) )/8.0

                vrn(i,m) = vn(i+1,m) + ( wr(2,i,m)*(-4.0*g(2,i,m)-s(2,i,m)) + &
                                         wr(3,i,m)*(-4.0*g(3,i,m)-s(3,i,m)) + &
                                         wr(4,i,m)*(-4.0*g(4,i,m)-s(4,i,m)) )/8.0
            end do
        end do
    end do
    
    do j=1,nsec
        call d5_via_T2N(ni,nst,ned,ngn,vn,nfs,nfe,nge,d5_li,sec_st(j),sec_ed(j))
        
        call vec2_via_T4N(ni,nst,ned,ngn,vn,nfs,nfe,nge,d2_li,sec_st(j),sec_ed(j))
        
        call vec4_via_T2N(ni,nst,ned,ngn,vn,nfs,nfe,nge,d4_li,sec_st(j),sec_ed(j))
    end do
    
    do m=nst,ned
        do j=1,nsec
            do i=sec_st(j),sec_ed(j)
                
                ve_lil(i,m) = 0.5*(vn(i+1,m) + vn(i  ,m)) - 1.0/8.0*d2_li(i,m) - 1.0/384.0*d4_li(i,m) -  aa_scs*d5_li(i,m)
                            
                ve_lir(i,m) = 0.5*(vn(i+1,m) + vn(i  ,m)) - 1.0/8.0*d2_li(i,m) - 1.0/384.0*d4_li(i,m) +  aa_scs*d5_li(i,m)

            end do
        end do
    end do
        
    do m=nst,ned
        do i=st,ed
            
               vl(i,m) = sigma(i)*ve_lil(i,m) + vln(i,m)*(1.0-sigma(i))
                         
               vr(i,m) = sigma(i)*ve_lir(i,m) + vrn(i,m)*(1.0-sigma(i))            

        end do
    end do  

    do i=stc,edc
        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vl(i,m) < chklim(m) ) then
                    vl(i,m) = vn(i  ,m)
                    ndec = ndec + 1
                end if

                if (vr(i,m) < chklim(m) ) then
                    vr(i,m) = vn(i+1,m)
                    ndec = ndec + 1
                end if
           end if
        end do
    end do

end subroutine scsh5pi3

subroutine scsh5ci3(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_bound_scheme
    implicit none
    integer(kind_int), intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  , intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(in)  :: sn(1-ngn:ni+ngn,1:3)
    integer(kind_int), intent(in)  :: nfs,nfe
    integer(kind_int), intent(in)  :: nge
    real(kind_real)  , intent(out) :: vl(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(out) :: vr(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , external    :: sub_limit
    integer(kind_int), intent(in)  :: chkid(nst:ned)
    real(kind_real),   intent(in)  :: chklim(nst:ned)
    integer(kind_int), intent(out) :: ndec  
    real(kind_real)   :: vln(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: vrn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: ve_lil(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: ve_lir(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: d2_li(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: d4_li(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: d5_li(1-ngn:ni+ngn,nst:ned)  
    real(kind_real)   :: nt,nx,ny,nz,fp(nst:ned)
    real(kind_real)   :: vp(nst:ned),vc(nst:ned)
    real(kind_real)   :: q(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: q1(1-ngn:ni+ngn,nst:ned)    
    real(kind_real)   :: s(1:4,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: g(1:4,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: wl(1:4,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: wr(1:4,1-ngn:ni+ngn,nst:ned)    
    real(kind_real)   :: aa,aa_scs
    real(kind_real)   :: sigma(1-ngn:ni+ngn),imax,imin,indicator
    real(kind_real), parameter :: cl(4) = (/1.D0/16.D0,10.D0/16.D0,5.D0/16.D0,0.D0/16.D0/)
    real(kind_real), parameter :: cr(4) = (/0.D0/16.D0,5.D0/16.D0,10.D0/16.D0,1.D0/16.D0/)      
    real(kind_real), parameter :: eps        = 1.0e-30
    real(kind_real), parameter :: indicatorc = 0.1
    real(kind_real)            :: is,eis2(4),bs,tau5,bl(4),br(4)
    integer(kind_int)          :: i,m,n,j,st,ed,stc,edc
    integer(kind_int)          :: id,YN,ind(1:ni+2*ngn),ind_st,ind_ed,pre_ed
    integer(kind_int)          :: nsec,sec_st(1:ni+2*ngn),sec_ed(1:ni+2*ngn)
    integer(kind_int)          :: nsecn,sec_stn(1:ni+2*ngn),sec_edn(1:ni+2*ngn)
    
    if (ni < 8) then
        call muscl2pv(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)
        return
    end if
    
    ndec = 0
    
    aa     = 0.01
    aa_scs = 0.001
    
    wl(:,:,:) = 0.3333
    wr(:,:,:) = 0.3333
    vln(:,:)    = 0.0
    ve_lil(:,:) = 0.0
    vrn(:,:)    = 0.0
    ve_lir(:,:) = 0.0    
           
    if (nfs == nbc_bound_scheme) then
        
        st  = 4
        stc = 0

        do m=nst,ned
            vl(0,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                  ! third-order
            vl(1,m) = (  3.0*vn(1,m) +  6.0*vn(2,m) -      vn(3,m) ) / 8.0   ! third-order 
            vl(2,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 + &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )
            vr(0,m) = vl(0,m)
            vr(1,m) = vl(1,m)
            vr(2,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 - &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )
            vl(3,m) =   ( 150.0*( vn(4,m) + vn(3,m) )          -  &
                           25.0*( vn(5,m) + vn(2,m) )          +  &
                            3.0*( vn(6,m) + vn(1,m) ) ) /256.0 -  &
                      aa*( 10.0*( vn(4,m) - vn(3,m) )          -  &
                            5.0*( vn(5,m) - vn(2,m) )          +  &
                                ( vn(6,m) - vn(1,m) ) ) 
            
            vr(3,m) =   ( 150.0*( vn(4,m) + vn(3,m) )          -  &
                           25.0*( vn(5,m) + vn(2,m) )          +  &
                            3.0*( vn(6,m) + vn(1,m) ) ) /256.0 +  &
                      aa*( 10.0*( vn(4,m) - vn(3,m) )          -  &
                            5.0*( vn(5,m) - vn(2,m) )          +  &
                                ( vn(6,m) - vn(1,m) ) )          
        end do
        
    else
        
        st  = -2
        stc = st
     
    end if
     
    if (nfe == nbc_bound_scheme) then
        
        ed  = ni - 4
        edc = ni

        do m=nst,ned
            vr(ni  ,m) = ( 15.0*vn(ni,m) - 10.0*vn(ni-1,m) +  3.0*vn(ni-2,m) ) / 8.0                ! third-order
            vl(ni-1,m) = (  3.0*vn(ni,m) +  6.0*vn(ni-1,m) -      vn(ni-2,m) ) / 8.0   ! third-order
            vl(ni-2,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 + &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) )
            vl(ni  ,m) = vr(ni  ,m)
            vr(ni-1,m) = vl(ni-1,m)
            vr(ni-2,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 - &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) ) 
            vl(ni-3,m) =   ( 150.0*( vn(ni-2,m) + vn(ni-3,m) )          -  &
                              25.0*( vn(ni-1,m) + vn(ni-4,m) )          +  &
                               3.0*( vn(ni  ,m) + vn(ni-5,m) ) ) /256.0 -  &
                         aa*( 10.0*( vn(ni-2,m) - vn(ni-3,m) )          -  &
                               5.0*( vn(ni-1,m) - vn(ni-4,m) )          +  &
                                   ( vn(ni  ,m) - vn(ni-5,m) ) ) 
            
            vr(ni-3,m) =   ( 150.0*( vn(ni-2,m) + vn(ni-3,m) )          -  &
                              25.0*( vn(ni-1,m) + vn(ni-4,m) )          +  &
                               3.0*( vn(ni  ,m) + vn(ni-5,m) ) ) /256.0 +  &
                         aa*( 10.0*( vn(ni-2,m) - vn(ni-3,m) )          -  &
                               5.0*( vn(ni-1,m) - vn(ni-4,m) )          +  &
                                   ( vn(ni  ,m) - vn(ni-5,m) ) )              
        end do
        
    else
        
        ed  = ni + 2
        edc = ed
            
    end if
    
    do m=nst,ned
        do i=st-2,ed+2
            g(1,i,m) = vn(i  ,m) - vn(i-1,m)            
            g(2,i,m) = vn(i+1,m) - vn(i  ,m)
            g(3,i,m) = vn(i+2,m) - vn(i+1,m)            
        end do
    end do   
    
    do m=nst,ned
        do i=st-2,ed+2
            do n=1,3
                is = g(n,i,m)*g(n,i,m)
                eis2(n) = is**2
            end do

            tau5 = abs(eis2(2)-eis2(1))
            bl(1) = 1.0+tau5/(eis2(1)+eps)
            bl(2) = 1.0+tau5/(eis2(2)+eps)
            bs = bl(1) + bl(2)
            wl(2,i,m) = bl(2)/bs
            
            tau5 = abs(eis2(3)-eis2(2))
            br(2) = 1.0+tau5/(eis2(2)+eps)
            br(3) = 1.0+tau5/(eis2(3)+eps)
            bs = br(2) + br(3)
            wr(2,i,m) = br(2)/bs
        end do
    end do
    
    do i=st,ed
        imin=1.0
        imax=0.0
        do m=nst,ned
            imin         = min(wl(2,i-3,m),wl(2,i-2,m),wl(2,i-1,m),wl(2,i,m),wl(2,i+1,m),wl(2,i+2,m),wl(2,i+3,m), &
                               wr(2,i-3,m),wr(2,i-2,m),wr(2,i-1,m),wr(2,i,m),wr(2,i+1,m),wr(2,i+2,m),wr(2,i+3,m),imin)                                
            imax         = max(wl(2,i-3,m),wl(2,i-2,m),wl(2,i-1,m),wl(2,i,m),wl(2,i+1,m),wl(2,i+2,m),wl(2,i+3,m), &
                               wr(2,i-3,m),wr(2,i-2,m),wr(2,i-1,m),wr(2,i,m),wr(2,i+1,m),wr(2,i+2,m),wr(2,i+3,m),imax)
        end do
        indicator    = imin/imax
        if(indicator>indicatorc) then
            sigma(i)=1.0
        else
            sigma(i)=0.0
        endif        
    end do
    
    id=0
    do i=st,ed
        YN=0
        if(sigma(i)<0.5) then
            YN=YN+1
        end if
        if(YN>0) then
            id=id+1
            ind(id)=i           
        end if
    end do
       
    do i=1-ngn,ni+ngn
        sigma(i)=0.0
    end do
    if(ed-st<9) then
        ind_st=st
        ind_ed=ed
        nsec=0
        
        nsecn=1
        sec_stn(1)=st
        sec_edn(1)=ed
    else
        if(id==0) then
            nsec=1
            sec_st(1)=st
            sec_ed(1)=ed
            do n=sec_st(1),sec_ed(1)
                sigma(n)=1.0
            end do
            ind_st=st
            ind_ed=ed
            nsecn=0
        else
            ind_st=st
            ind_ed=ed
            nsec=0
            nsecn=0
            if(ind(1)-ind_st>9) then
                nsec=1
                sec_st(1)=ind_st
                sec_ed(1)=ind(1)-1
                pre_ed   =sec_ed(1)
                do n=sec_st(1),sec_ed(1)
                    sigma(n)=1.0
                end do
            else
                nsecn=1
                sec_stn(1)=ind_st
                sec_edn(1)=ind(1)
                pre_ed    =sec_edn(1)
            end if
            do i=2,id
                if(ind(i)-ind(i-1)>9) then
                    nsec=nsec+1
                    sec_st(nsec)=ind(i-1)+1
                    sec_ed(nsec)=ind(i)-1
                    do n=sec_st(nsec),sec_ed(nsec)
                        sigma(n)=1.0
                    end do
                    nsecn=nsecn+1
                    sec_stn(nsecn)=pre_ed+1
                    sec_edn(nsecn)=ind(i-1)
                    pre_ed        =sec_ed(nsec)
                else
                    if(ind(i)-ind(i-1)>1 .or. i==id) then
                        nsecn=nsecn+1
                        sec_stn(nsecn)=pre_ed+1
                        sec_edn(nsecn)=ind(i)
                        pre_ed        =sec_edn(nsecn)
                    end if
                end if
            end do
            if(ind_ed-ind(id)>9) then
                nsec=nsec+1
                sec_st(nsec)=pre_ed+1
                sec_ed(nsec)=ind_ed
                do n=sec_st(nsec),sec_ed(nsec)
                    sigma(n)=1.0
                end do
            else
                nsecn=nsecn+1
                sec_stn(nsecn)=pre_ed+1
                sec_edn(nsecn)=ind_ed
            end if
        end if
    end if
    
    do m=nst,ned
        do j=1,nsecn
            do i=sec_stn(j),sec_edn(j)
                g(1,i,m) =      vn(i-2,m) - 3.0*vn(i-1,m) + 2.0*vn(i  ,m) 
                g(2,i,m) =     -vn(i  ,m) +     vn(i+1,m)
                g(3,i,m) =     -vn(i  ,m) +     vn(i+1,m)
                g(4,i,m) = -2.0*vn(i+1,m) + 3.0*vn(i+2,m) -     vn(i+3,m)

                s(1,i,m) = vn(i-2,m) - 2.0*vn(i-1,m) + vn(i,  m)
                s(2,i,m) = vn(i-1,m) - 2.0*vn(i  ,m) + vn(i+1,m) 
                s(3,i,m) = vn(i  ,m) - 2.0*vn(i+1,m) + vn(i+2,m)
                s(4,i,m) = vn(i+1,m) - 2.0*vn(i+2,m) + vn(i+3,m)
            end do
        end do      
    end do
    
    do j=1,nsecn
        do i=sec_stn(j),sec_edn(j)
            nx = (sn(i,1) + sn(i+1,1))/2.0
            ny = (sn(i,2) + sn(i+1,2))/2.0
            nz = (sn(i,3) + sn(i+1,3))/2.0

            do m=nst,ned
                vp(m) = (vn(i,m) + vn(i+1,m))/2.0
                fp(m) = vn(i,m)
            end do
            call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
            do m=nst,ned
                q(i,m) = vc(m)
            end do
            
            do m=nst,ned
                fp(m) = vn(i+1,m)
            end do
            call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
            do m=nst,ned
                q1(i+1,m) = vc(m)
            end do        
             
            do n=1,4
                do m=nst,ned
                    fp(m) = g(n,i,m)
                end do
                call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
                do m=nst,ned
                    g(n,i,m) = vc(m)
                end do

                do m=nst,ned
                    fp(m) = s(n,i,m)
                end do
                call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
                do m=nst,ned
                    s(n,i,m) = vc(m)
                end do          
            end do
        end do
    end do   
            
            
    do m=nst,ned
        do j=1,nsecn
            do i=sec_stn(j),sec_edn(j)
                do n=1,4
                    is = g(n,i,m)*g(n,i,m) + s(n,i,m)*s(n,i,m)
                    eis2(n) = is**2
                end do
                
                tau5 = abs(eis2(3)-eis2(1))
                
                do n=1,3
                    bl(n) = cl(n)*(1.0+tau5/(eis2(n)+eps))
                end do
                
                bs = bl(1) + bl(2) + bl(3)
                do n=1,3
                    wl(n,i,m) = bl(n)/bs
                end do
          
                tau5 = abs(eis2(4)-eis2(2))           
                do n=2,4
                    br(n) = cr(n)*(1.0+tau5/(eis2(n)+eps))
                end do

                bs = br(2) + br(3) + br(4)
                do n=2,4
                    wr(n,i,m) = br(n)/bs
                end do
            end do
        end do
    end do
    
!!    open(2,file='sigma.dat',form='formatted')
!!    write(2,*) 'variables="i" "sigmal" "sigmar"'
!!    do i=st,ed
!!
!!        write(2,200) i,sigmal(i,1),sigmar(i,1)
!!
!!    end do
!!
!!200 format(I4,1x,e12.5,1x,e12.5)
!!    close(2)
!!
!!    stop    
    
    do m=nst,ned
        do j=1,nsecn
            do i=sec_stn(j),sec_edn(j)
                vln(i,m) =  q(i  ,m) + ( wl(1,i,m)*( 4.0*g(1,i,m)-s(1,i,m)) + &
                                         wl(2,i,m)*( 4.0*g(2,i,m)-s(2,i,m)) + &
                                         wl(3,i,m)*( 4.0*g(3,i,m)-s(3,i,m)) )/8.0

                vrn(i,m) = q1(i+1,m) + ( wr(2,i,m)*(-4.0*g(2,i,m)-s(2,i,m)) + &
                                         wr(3,i,m)*(-4.0*g(3,i,m)-s(3,i,m)) + &
                                         wr(4,i,m)*(-4.0*g(4,i,m)-s(4,i,m)) )/8.0
            end do
        end do
    end do
    
    do j=1,nsecn
        do i=sec_stn(j),sec_edn(j)
            nx = (sn(i,1) + sn(i+1,1))/2.0
            ny = (sn(i,2) + sn(i+1,2))/2.0
            nz = (sn(i,3) + sn(i+1,3))/2.0

            do m=nst,ned
                vp(m) = (vn(i,m) + vn(i+1,m))/2.0
                vc(m) = vln(i,m)
            end do

            call RxQ_p(nst,ned,vp,nt,nx,ny,nz,vc,fp) 
            do m=nst,ned
                vln(i,m) = fp(m)
            end do
            
            do m=nst,ned
                vc(m) = vrn(i ,m)
            end do

            call RxQ_p(nst,ned,vp,nt,nx,ny,nz,vc,fp) 
            do m=nst,ned
                vrn(i,m) = fp(m)
            end do
        end do
    end do     
      
    do j=1,nsec
        call d5_via_T2N(ni,nst,ned,ngn,vn,nfs,nfe,nge,d5_li,sec_st(j),sec_ed(j))
        
        call vec2_via_T4N(ni,nst,ned,ngn,vn,nfs,nfe,nge,d2_li,sec_st(j),sec_ed(j))
        
        call vec4_via_T2N(ni,nst,ned,ngn,vn,nfs,nfe,nge,d4_li,sec_st(j),sec_ed(j))
    end do
    
    do m=nst,ned
        do j=1,nsec
            do i=sec_st(j),sec_ed(j)
                
                ve_lil(i,m) = 0.5*(vn(i+1,m) + vn(i  ,m)) - 1.0/8.0*d2_li(i,m) - 1.0/384.0*d4_li(i,m) -  aa_scs*d5_li(i,m)
                            
                ve_lir(i,m) = 0.5*(vn(i+1,m) + vn(i  ,m)) - 1.0/8.0*d2_li(i,m) - 1.0/384.0*d4_li(i,m) +  aa_scs*d5_li(i,m)

            end do
        end do
    end do
        
    do m=nst,ned
        do i=st,ed
            
               vl(i,m) = sigma(i)*ve_lil(i,m) + vln(i,m)*(1.0-sigma(i))
                         
               vr(i,m) = sigma(i)*ve_lir(i,m) + vrn(i,m)*(1.0-sigma(i))            

        end do
    end do  

    do i=stc,edc
        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vl(i,m) < chklim(m) ) then
                    vl(i,m) = vn(i  ,m)
                    ndec = ndec + 1
                end if

                if (vr(i,m) < chklim(m) ) then
                    vr(i,m) = vn(i+1,m)
                    ndec = ndec + 1
                end if
           end if
        end do
    end do

end subroutine scsh5ci3

subroutine scsh5pi(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_bound_scheme
    implicit none
    integer(kind_int), intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  , intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(in)  :: sn(1-ngn:ni+ngn,1:3)
    integer(kind_int), intent(in)  :: nfs,nfe
    integer(kind_int), intent(in)  :: nge
    real(kind_real)  , intent(out) :: vl(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(out) :: vr(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , external    :: sub_limit
    integer(kind_int), intent(in)  :: chkid(nst:ned)
    real(kind_real),   intent(in)  :: chklim(nst:ned)
    integer(kind_int), intent(out) :: ndec  
    real(kind_real)   :: vln(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: vrn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: ve_lil(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: ve_lir(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: d2_li(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: d4_li(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: d5_li(1-ngn:ni+ngn,nst:ned)   
    real(kind_real)   :: s(1:5,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: g(1:5,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: t(1:5,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: wl(1:5,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: wr(1:5,1-ngn:ni+ngn,nst:ned)    
    real(kind_real)   :: aa,aa_scs
    real(kind_real)   :: sigma(1-ngn:ni+ngn),imax,imin,indicator
    real(kind_real), parameter :: cl(5) = (/1.D0/64.D0,21.D0/64.D0,35.D0/64.D0,7.D0/64.D0,0.D0/64.D0/)
    real(kind_real), parameter :: cr(5) = (/0.D0/64.D0,7.D0/64.D0,35.D0/64.D0,21.D0/64.D0,1.D0/64.D0/)  
    real(kind_real), parameter :: eps        = 1.0e-30
    real(kind_real), parameter :: indicatorc = 0.1
    real(kind_real)            :: is,eis2(5),bs,tau7,bl(5),br(5)
    integer(kind_int)          :: i,m,n,j,st,ed,stc,edc
    integer(kind_int)          :: id,YN,ind(1:ni+2*ngn),ind_st,ind_ed,pre_ed
    integer(kind_int)          :: nsec,sec_st(1:ni+2*ngn),sec_ed(1:ni+2*ngn)
    integer(kind_int)          :: nsecn,sec_stn(1:ni+2*ngn),sec_edn(1:ni+2*ngn)
    
    if (ni < 8) then
        call muscl2pv(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)
        return
    end if
    
    ndec = 0
    
    aa     = 0.01
    aa_scs = 0.001
    
    wl(:,:,:) = 0.25
    wr(:,:,:) = 0.25
    vln(:,:)    = 0.0
    ve_lil(:,:) = 0.0
    vrn(:,:)    = 0.0
    ve_lir(:,:) = 0.0    
           
    if (nfs == nbc_bound_scheme) then
        
        st  = 4
        stc = 0

        do m=nst,ned
            vl(0,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                  ! third-order
            vl(1,m) = (  3.0*vn(1,m) +  6.0*vn(2,m) -      vn(3,m) ) / 8.0   ! third-order 
            vl(2,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 + &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )
            vr(0,m) = vl(0,m)
            vr(1,m) = vl(1,m)
            vr(2,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 - &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )
            vl(3,m) =   ( 150.0*( vn(4,m) + vn(3,m) )          -  &
                           25.0*( vn(5,m) + vn(2,m) )          +  &
                            3.0*( vn(6,m) + vn(1,m) ) ) /256.0 -  &
                      aa*( 10.0*( vn(4,m) - vn(3,m) )          -  &
                            5.0*( vn(5,m) - vn(2,m) )          +  &
                                ( vn(6,m) - vn(1,m) ) ) 
            
            vr(3,m) =   ( 150.0*( vn(4,m) + vn(3,m) )          -  &
                           25.0*( vn(5,m) + vn(2,m) )          +  &
                            3.0*( vn(6,m) + vn(1,m) ) ) /256.0 +  &
                      aa*( 10.0*( vn(4,m) - vn(3,m) )          -  &
                            5.0*( vn(5,m) - vn(2,m) )          +  &
                                ( vn(6,m) - vn(1,m) ) )          
        end do
        
    else
        
        st  = -2
        stc = st
     
    end if
     
    if (nfe == nbc_bound_scheme) then
        
        ed  = ni - 4
        edc = ni

        do m=nst,ned
            vr(ni  ,m) = ( 15.0*vn(ni,m) - 10.0*vn(ni-1,m) +  3.0*vn(ni-2,m) ) / 8.0                ! third-order
            vl(ni-1,m) = (  3.0*vn(ni,m) +  6.0*vn(ni-1,m) -      vn(ni-2,m) ) / 8.0   ! third-order
            vl(ni-2,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 + &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) )
            vl(ni  ,m) = vr(ni  ,m)
            vr(ni-1,m) = vl(ni-1,m)
            vr(ni-2,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 - &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) ) 
            vl(ni-3,m) =   ( 150.0*( vn(ni-2,m) + vn(ni-3,m) )          -  &
                              25.0*( vn(ni-1,m) + vn(ni-4,m) )          +  &
                               3.0*( vn(ni  ,m) + vn(ni-5,m) ) ) /256.0 -  &
                         aa*( 10.0*( vn(ni-2,m) - vn(ni-3,m) )          -  &
                               5.0*( vn(ni-1,m) - vn(ni-4,m) )          +  &
                                   ( vn(ni  ,m) - vn(ni-5,m) ) ) 
            
            vr(ni-3,m) =   ( 150.0*( vn(ni-2,m) + vn(ni-3,m) )          -  &
                              25.0*( vn(ni-1,m) + vn(ni-4,m) )          +  &
                               3.0*( vn(ni  ,m) + vn(ni-5,m) ) ) /256.0 +  &
                         aa*( 10.0*( vn(ni-2,m) - vn(ni-3,m) )          -  &
                               5.0*( vn(ni-1,m) - vn(ni-4,m) )          +  &
                                   ( vn(ni  ,m) - vn(ni-5,m) ) )              
        end do
        
    else
        
        ed  = ni + 2
        edc = ed
            
    end if
    
    do m=nst,ned
        do i=st-2,ed+2
            g(1,i,m) = vn(i  ,m) - vn(i-1,m)            
            g(2,i,m) = vn(i+1,m) - vn(i  ,m)
            g(3,i,m) = vn(i+2,m) - vn(i+1,m)            
        end do
    end do   
    
    do m=nst,ned
        do i=st-2,ed+2
            do n=1,3
                is = g(n,i,m)*g(n,i,m)
                eis2(n) = is**2
            end do

            tau7 = abs(eis2(2)-eis2(1))
            bl(1) = 1.0+tau7/(eis2(1)+eps)
            bl(2) = 1.0+tau7/(eis2(2)+eps)
            bs = bl(1) + bl(2)
            wl(2,i,m) = bl(2)/bs
            
            tau7 = abs(eis2(3)-eis2(2))
            br(2) = 1.0+tau7/(eis2(2)+eps)
            br(3) = 1.0+tau7/(eis2(3)+eps)
            bs = br(2) + br(3)
            wr(2,i,m) = br(2)/bs
        end do
    end do
    
    do i=st,ed
        imin=1.0
        imax=0.0
        do m=nst,ned
            imin         = min(wl(2,i-3,m),wl(2,i-2,m),wl(2,i-1,m),wl(2,i,m),wl(2,i+1,m),wl(2,i+2,m),wl(2,i+3,m), &
                               wr(2,i-3,m),wr(2,i-2,m),wr(2,i-1,m),wr(2,i,m),wr(2,i+1,m),wr(2,i+2,m),wr(2,i+3,m),imin)                                
            imax         = max(wl(2,i-3,m),wl(2,i-2,m),wl(2,i-1,m),wl(2,i,m),wl(2,i+1,m),wl(2,i+2,m),wl(2,i+3,m), &
                               wr(2,i-3,m),wr(2,i-2,m),wr(2,i-1,m),wr(2,i,m),wr(2,i+1,m),wr(2,i+2,m),wr(2,i+3,m),imax)
        end do
        indicator    = imin/imax
        if(indicator>indicatorc) then
            sigma(i)=1.0
        else
            sigma(i)=0.0
        endif        
    end do
    
    id=0
    do i=st,ed
        YN=0
        if(sigma(i)<0.5) then
            YN=YN+1
        end if
        if(YN>0) then
            id=id+1
            ind(id)=i           
        end if
    end do
       
    do i=1-ngn,ni+ngn
        sigma(i)=0.0
    end do
    if(ed-st<9) then
        ind_st=st
        ind_ed=ed
        nsec=0
        
        nsecn=1
        sec_stn(1)=st
        sec_edn(1)=ed
    else
        if(id==0) then
            nsec=1
            sec_st(1)=st
            sec_ed(1)=ed
            do n=sec_st(1),sec_ed(1)
                sigma(n)=1.0
            end do
            ind_st=st
            ind_ed=ed
            nsecn=0
        else
            ind_st=st
            ind_ed=ed
            nsec=0
            nsecn=0
            if(ind(1)-ind_st>9) then
                nsec=1
                sec_st(1)=ind_st
                sec_ed(1)=ind(1)-1
                pre_ed   =sec_ed(1)
                do n=sec_st(1),sec_ed(1)
                    sigma(n)=1.0
                end do
            else
                nsecn=1
                sec_stn(1)=ind_st
                sec_edn(1)=ind(1)
                pre_ed    =sec_edn(1)
            end if
            do i=2,id
                if(ind(i)-ind(i-1)>9) then
                    nsec=nsec+1
                    sec_st(nsec)=ind(i-1)+1
                    sec_ed(nsec)=ind(i)-1
                    do n=sec_st(nsec),sec_ed(nsec)
                        sigma(n)=1.0
                    end do
                    nsecn=nsecn+1
                    sec_stn(nsecn)=pre_ed+1
                    sec_edn(nsecn)=ind(i-1)
                    pre_ed        =sec_ed(nsec)
                else
                    if(ind(i)-ind(i-1)>1 .or. i==id) then
                        nsecn=nsecn+1
                        sec_stn(nsecn)=pre_ed+1
                        sec_edn(nsecn)=ind(i)
                        pre_ed        =sec_edn(nsecn)
                    end if
                end if
            end do
            if(ind_ed-ind(id)>9) then
                nsec=nsec+1
                sec_st(nsec)=pre_ed+1
                sec_ed(nsec)=ind_ed
                do n=sec_st(nsec),sec_ed(nsec)
                    sigma(n)=1.0
                end do
            else
                nsecn=nsecn+1
                sec_stn(nsecn)=pre_ed+1
                sec_edn(nsecn)=ind_ed
            end if
        end if
    end if
    
    do m=nst,ned
        do j=1,nsecn
            do i=sec_stn(j),sec_edn(j)
                g(1,i,m) = -23.0/24.0*vn(i-3,m) + 31.0/8.0*vn(i-2,m) - 47.0/8.0*vn(i-1,m) + 71.0/24.0*vn(i  ,m)
                g(2,i,m) =   1.0/24.0*vn(i-2,m) -  1.0/8.0*vn(i-1,m) -  7.0/8.0*vn(i  ,m) + 23.0/24.0*vn(i+1,m)
                g(3,i,m) =   1.0/24.0*vn(i-1,m) -  9.0/8.0*vn(i  ,m) +  9.0/8.0*vn(i+1,m) -  1.0/24.0*vn(i+2,m)
                g(4,i,m) = -23.0/24.0*vn(i  ,m) +  7.0/8.0*vn(i+1,m) +  1.0/8.0*vn(i+2,m) -  1.0/24.0*vn(i+3,m)
                g(5,i,m) = -71.0/24.0*vn(i+1,m) + 47.0/8.0*vn(i+2,m) - 31.0/8.0*vn(i+3,m) + 23.0/24.0*vn(i+4,m)
                
                s(1,i,m) =   -3.0/2.0*vn(i-3,m) + 11.0/2.0*vn(i-2,m) - 13.0/2.0*vn(i-1,m) +   5.0/2.0*vn(i  ,m)
                s(2,i,m) =   -1.0/2.0*vn(i-2,m) +  5.0/2.0*vn(i-1,m) -  7.0/2.0*vn(i  ,m) +   3.0/2.0*vn(i+1,m)
                s(3,i,m) =    1.0/2.0*vn(i-1,m) -  1.0/2.0*vn(i  ,m) -  1.0/2.0*vn(i+1,m) +   1.0/2.0*vn(i+2,m)
                s(4,i,m) =    3.0/2.0*vn(i  ,m) -  7.0/2.0*vn(i+1,m) +  5.0/2.0*vn(i+2,m) -   1.0/2.0*vn(i+3,m)
                s(5,i,m) =    5.0/2.0*vn(i+1,m) - 13.0/2.0*vn(i+2,m) + 11.0/2.0*vn(i+3,m) -   3.0/2.0*vn(i+4,m)
                
                t(1,i,m) =    -vn(i-3,m) + 3.0*vn(i-2,m) - 3.0*vn(i-1,m) + vn(i  ,m)
                t(2,i,m) =    -vn(i-2,m) + 3.0*vn(i-1,m) - 3.0*vn(i  ,m) + vn(i+1,m)
                t(3,i,m) =    -vn(i-1,m) + 3.0*vn(i  ,m) - 3.0*vn(i+1,m) + vn(i+2,m)
                t(4,i,m) =    -vn(i  ,m) + 3.0*vn(i+1,m) - 3.0*vn(i+2,m) + vn(i+3,m)
                t(5,i,m) =    -vn(i+1,m) + 3.0*vn(i+2,m) - 3.0*vn(i+3,m) + vn(i+4,m)
            end do
        end do      
    end do            
            
    do m=nst,ned
        do j=1,nsecn
            do i=sec_stn(j),sec_edn(j)
                do n=1,5
                    is = g(n,i,m)*g(n,i,m) + s(n,i,m)*s(n,i,m) + t(n,i,m)*t(n,i,m)              
                    eis2(n) = is**2
                end do
                
                tau7 = abs(eis2(4)-eis2(1))
                
                do n=1,4
                    bl(n) = cl(n)*(1.0+tau7/(eis2(n)+eps))
                end do
                
                bs = bl(1) + bl(2) + bl(3) + bl(4)
                do n=1,4
                    wl(n,i,m) = bl(n)/bs
                end do
          
                tau7 = abs(eis2(5)-eis2(2))           
                do n=2,5
                    br(n) = cr(n)*(1.0+tau7/(eis2(n)+eps))
                end do

                bs = br(2) + br(3) + br(4) + br(5)
                do n=2,5
                    wr(n,i,m) = br(n)/bs
                end do
            end do
        end do
    end do
    
!!    open(2,file='sigma.dat',form='formatted')
!!    write(2,*) 'variables="i" "sigmal" "sigmar"'
!!    do i=st,ed
!!
!!        write(2,200) i,sigmal(i,1),sigmar(i,1)
!!
!!    end do
!!
!!200 format(I4,1x,e12.5,1x,e12.5)
!!    close(2)
!!
!!    stop    
    
    do m=nst,ned
        do j=1,nsecn
            do i=sec_stn(j),sec_edn(j)
                vln(i,m) = vn(i  ,m) + ( wl(1,i,m)*( 24.0*g(1,i,m)-6.0*s(1,i,m)+t(1,i,m)) + &
                                         wl(2,i,m)*( 24.0*g(2,i,m)-6.0*s(2,i,m)+t(2,i,m)) + &
                                         wl(3,i,m)*( 24.0*g(3,i,m)-6.0*s(3,i,m)+t(3,i,m)) + &
                                         wl(4,i,m)*( 24.0*g(4,i,m)-6.0*s(4,i,m)+t(4,i,m)) )/48.0

                vrn(i,m) = vn(i+1,m) + ( wr(2,i,m)*(-24.0*g(2,i,m)-6.0*s(2,i,m)-t(2,i,m)) + &
                                         wr(3,i,m)*(-24.0*g(3,i,m)-6.0*s(3,i,m)-t(3,i,m)) + &
                                         wr(4,i,m)*(-24.0*g(4,i,m)-6.0*s(4,i,m)-t(4,i,m)) + &
                                         wr(5,i,m)*(-24.0*g(5,i,m)-6.0*s(5,i,m)-t(5,i,m)) )/48.0
            end do
        end do
    end do    
      
    do j=1,nsec
        call d5_via_T2N(ni,nst,ned,ngn,vn,nfs,nfe,nge,d5_li,sec_st(j),sec_ed(j))
        
        call vec2_via_T4N(ni,nst,ned,ngn,vn,nfs,nfe,nge,d2_li,sec_st(j),sec_ed(j))
        
        call vec4_via_T2N(ni,nst,ned,ngn,vn,nfs,nfe,nge,d4_li,sec_st(j),sec_ed(j))
    end do
    
    do m=nst,ned
        do j=1,nsec
            do i=sec_st(j),sec_ed(j)
                
                ve_lil(i,m) = 0.5*(vn(i+1,m) + vn(i  ,m)) - 1.0/8.0*d2_li(i,m) - 1.0/384.0*d4_li(i,m) -  aa_scs*d5_li(i,m)
                            
                ve_lir(i,m) = 0.5*(vn(i+1,m) + vn(i  ,m)) - 1.0/8.0*d2_li(i,m) - 1.0/384.0*d4_li(i,m) +  aa_scs*d5_li(i,m)

            end do
        end do
    end do
        
    do m=nst,ned
        do i=st,ed
            
               vl(i,m) = sigma(i)*ve_lil(i,m) + vln(i,m)*(1.0-sigma(i))
                         
               vr(i,m) = sigma(i)*ve_lir(i,m) + vrn(i,m)*(1.0-sigma(i))            

        end do
    end do  

    do i=stc,edc
        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vl(i,m) < chklim(m) ) then
                    vl(i,m) = vn(i  ,m)
                    ndec = ndec + 1
                end if

                if (vr(i,m) < chklim(m) ) then
                    vr(i,m) = vn(i+1,m)
                    ndec = ndec + 1
                end if
           end if
        end do
    end do

end subroutine scsh5pi

subroutine scsh5picc(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_bound_scheme
    implicit none
    integer(kind_int), intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  , intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(in)  :: sn(1-ngn:ni+ngn,1:3)
    integer(kind_int), intent(in)  :: nfs,nfe
    integer(kind_int), intent(in)  :: nge
    real(kind_real)  , intent(out) :: vl(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(out) :: vr(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , external    :: sub_limit
    integer(kind_int), intent(in)  :: chkid(nst:ned)
    real(kind_real),   intent(in)  :: chklim(nst:ned)
    integer(kind_int), intent(out) :: ndec  
    real(kind_real)   :: vln(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: vrn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: ve_lil(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: ve_lir(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: d2_li(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: d4_li(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: d5_li(1-ngn:ni+ngn,nst:ned)   
    real(kind_real)   :: s(1:5,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: g(1:5,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: t(1:5,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: wl(1:5,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: wr(1:5,1-ngn:ni+ngn,nst:ned)    
    real(kind_real)   :: aa,aa_scs
    real(kind_real)   :: sigma(1-ngn:ni+ngn),imax,imin,indicator
    real(kind_real), parameter :: cl(5) = (/1.D0/64.D0,21.D0/64.D0,35.D0/64.D0,7.D0/64.D0,0.D0/64.D0/)
    real(kind_real), parameter :: cr(5) = (/0.D0/64.D0,7.D0/64.D0,35.D0/64.D0,21.D0/64.D0,1.D0/64.D0/)  
    real(kind_real), parameter :: eps        = 1.0e-30
    real(kind_real), parameter :: indicatorc = 0.1
    real(kind_real)            :: is,eis2(5),bs,tau7,bl(5),br(5)
    integer(kind_int)          :: i,m,n,j,st,ed,stc,edc
    integer(kind_int)          :: id,YN,ind(1:ni+2*ngn),ind_st,ind_ed,pre_ed
    integer(kind_int)          :: nsec,sec_st(1:ni+2*ngn),sec_ed(1:ni+2*ngn)
    integer(kind_int)          :: nsecn,sec_stn(1:ni+2*ngn),sec_edn(1:ni+2*ngn)
    
    if (ni < 8) then
        call muscl2pvcc(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)
        return
    end if
    
    ndec = 0
    
    aa     = 0.01
    aa_scs = 0.001
    
    wl(:,:,:) = 0.25
    wr(:,:,:) = 0.25
    vln(:,:)    = 0.0
    ve_lil(:,:) = 0.0
    vrn(:,:)    = 0.0
    ve_lir(:,:) = 0.0    
           
    if (nfs == nbc_bound_scheme) then
        
        st  = 5
        stc = 1

        do m=nst,ned
            vl(1,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                  ! third-order
            vl(2,m) = (  3.0*vn(1,m) +  6.0*vn(2,m) -      vn(3,m) ) / 8.0   ! third-order 
            vl(3,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 + &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )
            vr(1,m) = vl(1,m)
            vr(2,m) = vl(2,m)
            vr(3,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 - &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )
            vl(4,m) =   ( 150.0*( vn(4,m) + vn(3,m) )          -  &
                           25.0*( vn(5,m) + vn(2,m) )          +  &
                            3.0*( vn(6,m) + vn(1,m) ) ) /256.0 -  &
                      aa*( 10.0*( vn(4,m) - vn(3,m) )          -  &
                            5.0*( vn(5,m) - vn(2,m) )          +  &
                                ( vn(6,m) - vn(1,m) ) ) 
            
            vr(4,m) =   ( 150.0*( vn(4,m) + vn(3,m) )          -  &
                           25.0*( vn(5,m) + vn(2,m) )          +  &
                            3.0*( vn(6,m) + vn(1,m) ) ) /256.0 +  &
                      aa*( 10.0*( vn(4,m) - vn(3,m) )          -  &
                            5.0*( vn(5,m) - vn(2,m) )          +  &
                                ( vn(6,m) - vn(1,m) ) )          
        end do
        
    else
        
        st  = -1
        stc = st
     
    end if
     
    if (nfe == nbc_bound_scheme) then
        
        ed  = ni - 3
        edc = ni + 1

        do m=nst,ned
            vr(ni+1,m) = ( 15.0*vn(ni,m) - 10.0*vn(ni-1,m) +  3.0*vn(ni-2,m) ) / 8.0                ! third-order
            vl(ni  ,m) = (  3.0*vn(ni,m) +  6.0*vn(ni-1,m) -      vn(ni-2,m) ) / 8.0   ! third-order
            vl(ni-1,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 + &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) )
            vl(ni+1,m) = vr(ni+1,m)
            vr(ni  ,m) = vl(ni  ,m)
            vr(ni-1,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 - &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) ) 
            vl(ni-2,m) =   ( 150.0*( vn(ni-2,m) + vn(ni-3,m) )          -  &
                              25.0*( vn(ni-1,m) + vn(ni-4,m) )          +  &
                               3.0*( vn(ni  ,m) + vn(ni-5,m) ) ) /256.0 -  &
                         aa*( 10.0*( vn(ni-2,m) - vn(ni-3,m) )          -  &
                               5.0*( vn(ni-1,m) - vn(ni-4,m) )          +  &
                                   ( vn(ni  ,m) - vn(ni-5,m) ) ) 
            
            vr(ni-2,m) =   ( 150.0*( vn(ni-2,m) + vn(ni-3,m) )          -  &
                              25.0*( vn(ni-1,m) + vn(ni-4,m) )          +  &
                               3.0*( vn(ni  ,m) + vn(ni-5,m) ) ) /256.0 +  &
                         aa*( 10.0*( vn(ni-2,m) - vn(ni-3,m) )          -  &
                               5.0*( vn(ni-1,m) - vn(ni-4,m) )          +  &
                                   ( vn(ni  ,m) - vn(ni-5,m) ) )              
        end do
        
    else
        
        ed  = ni + 3
        edc = ed
            
    end if
    
    do m=nst,ned
        do i=st-3,ed+1
            g(1,i,m) = vn(i  ,m) - vn(i-1,m)            
            g(2,i,m) = vn(i+1,m) - vn(i  ,m)
            g(3,i,m) = vn(i+2,m) - vn(i+1,m)            
        end do
    end do   
    
    do m=nst,ned
        do i=st-3,ed+1
            do n=1,3
                is = g(n,i,m)*g(n,i,m)
                eis2(n) = is**2
            end do

            tau7 = abs(eis2(2)-eis2(1))
            bl(1) = 1.0+tau7/(eis2(1)+eps)
            bl(2) = 1.0+tau7/(eis2(2)+eps)
            bs = bl(1) + bl(2)
            wl(2,i,m) = bl(2)/bs
            
            tau7 = abs(eis2(3)-eis2(2))
            br(2) = 1.0+tau7/(eis2(2)+eps)
            br(3) = 1.0+tau7/(eis2(3)+eps)
            bs = br(2) + br(3)
            wr(2,i,m) = br(2)/bs
        end do
    end do
    
    do i=st-1,ed-1
        imin=1.0
        imax=0.0
        do m=nst,ned
            imin         = min(wl(2,i-3,m),wl(2,i-2,m),wl(2,i-1,m),wl(2,i,m),wl(2,i+1,m),wl(2,i+2,m),wl(2,i+3,m), &
                               wr(2,i-3,m),wr(2,i-2,m),wr(2,i-1,m),wr(2,i,m),wr(2,i+1,m),wr(2,i+2,m),wr(2,i+3,m),imin)                                
            imax         = max(wl(2,i-3,m),wl(2,i-2,m),wl(2,i-1,m),wl(2,i,m),wl(2,i+1,m),wl(2,i+2,m),wl(2,i+3,m), &
                               wr(2,i-3,m),wr(2,i-2,m),wr(2,i-1,m),wr(2,i,m),wr(2,i+1,m),wr(2,i+2,m),wr(2,i+3,m),imax)
        end do
        indicator    = imin/imax
        if(indicator>indicatorc) then
            sigma(i)=1.0
        else
            sigma(i)=0.0
        endif        
    end do
    
    id=0
    do i=st-1,ed-1
        YN=0
        if(sigma(i)<0.5) then
            YN=YN+1
        end if
        if(YN>0) then
            id=id+1
            ind(id)=i           
        end if
    end do
       
    do i=1-ngn,ni+ngn
        sigma(i)=0.0
    end do
    if(ed-st<9) then
        ind_st=st-1
        ind_ed=ed-1
        nsec=0
        
        nsecn=1
        sec_stn(1)=st-1
        sec_edn(1)=ed-1
    else
        if(id==0) then
            nsec=1
            sec_st(1)=st-1
            sec_ed(1)=ed-1
            do n=sec_st(1),sec_ed(1)
                sigma(n)=1.0
            end do
            ind_st=st-1
            ind_ed=ed-1
            nsecn=0
        else
            ind_st=st-1
            ind_ed=ed-1
            nsec=0
            nsecn=0
            if(ind(1)-ind_st>9) then
                nsec=1
                sec_st(1)=ind_st
                sec_ed(1)=ind(1)-1
                pre_ed   =sec_ed(1)
                do n=sec_st(1),sec_ed(1)
                    sigma(n)=1.0
                end do
            else
                nsecn=1
                sec_stn(1)=ind_st
                sec_edn(1)=ind(1)
                pre_ed    =sec_edn(1)
            end if
            do i=2,id
                if(ind(i)-ind(i-1)>9) then
                    nsec=nsec+1
                    sec_st(nsec)=ind(i-1)+1
                    sec_ed(nsec)=ind(i)-1
                    do n=sec_st(nsec),sec_ed(nsec)
                        sigma(n)=1.0
                    end do
                    nsecn=nsecn+1
                    sec_stn(nsecn)=pre_ed+1
                    sec_edn(nsecn)=ind(i-1)
                    pre_ed        =sec_ed(nsec)
                else
                    if(ind(i)-ind(i-1)>1 .or. i==id) then
                        nsecn=nsecn+1
                        sec_stn(nsecn)=pre_ed+1
                        sec_edn(nsecn)=ind(i)
                        pre_ed        =sec_edn(nsecn)
                    end if
                end if
            end do
            if(ind_ed-ind(id)>9) then
                nsec=nsec+1
                sec_st(nsec)=pre_ed+1
                sec_ed(nsec)=ind_ed
                do n=sec_st(nsec),sec_ed(nsec)
                    sigma(n)=1.0
                end do
            else
                nsecn=nsecn+1
                sec_stn(nsecn)=pre_ed+1
                sec_edn(nsecn)=ind_ed
            end if
        end if
    end if
    
    do m=nst,ned
        do j=1,nsecn
            do i=sec_stn(j),sec_edn(j)
                g(1,i,m) = -23.0/24.0*vn(i-3,m) + 31.0/8.0*vn(i-2,m) - 47.0/8.0*vn(i-1,m) + 71.0/24.0*vn(i  ,m)
                g(2,i,m) =   1.0/24.0*vn(i-2,m) -  1.0/8.0*vn(i-1,m) -  7.0/8.0*vn(i  ,m) + 23.0/24.0*vn(i+1,m)
                g(3,i,m) =   1.0/24.0*vn(i-1,m) -  9.0/8.0*vn(i  ,m) +  9.0/8.0*vn(i+1,m) -  1.0/24.0*vn(i+2,m)
                g(4,i,m) = -23.0/24.0*vn(i  ,m) +  7.0/8.0*vn(i+1,m) +  1.0/8.0*vn(i+2,m) -  1.0/24.0*vn(i+3,m)
                g(5,i,m) = -71.0/24.0*vn(i+1,m) + 47.0/8.0*vn(i+2,m) - 31.0/8.0*vn(i+3,m) + 23.0/24.0*vn(i+4,m)
                
                s(1,i,m) =   -3.0/2.0*vn(i-3,m) + 11.0/2.0*vn(i-2,m) - 13.0/2.0*vn(i-1,m) +   5.0/2.0*vn(i  ,m)
                s(2,i,m) =   -1.0/2.0*vn(i-2,m) +  5.0/2.0*vn(i-1,m) -  7.0/2.0*vn(i  ,m) +   3.0/2.0*vn(i+1,m)
                s(3,i,m) =    1.0/2.0*vn(i-1,m) -  1.0/2.0*vn(i  ,m) -  1.0/2.0*vn(i+1,m) +   1.0/2.0*vn(i+2,m)
                s(4,i,m) =    3.0/2.0*vn(i  ,m) -  7.0/2.0*vn(i+1,m) +  5.0/2.0*vn(i+2,m) -   1.0/2.0*vn(i+3,m)
                s(5,i,m) =    5.0/2.0*vn(i+1,m) - 13.0/2.0*vn(i+2,m) + 11.0/2.0*vn(i+3,m) -   3.0/2.0*vn(i+4,m)
                
                t(1,i,m) =    -vn(i-3,m) + 3.0*vn(i-2,m) - 3.0*vn(i-1,m) + vn(i  ,m)
                t(2,i,m) =    -vn(i-2,m) + 3.0*vn(i-1,m) - 3.0*vn(i  ,m) + vn(i+1,m)
                t(3,i,m) =    -vn(i-1,m) + 3.0*vn(i  ,m) - 3.0*vn(i+1,m) + vn(i+2,m)
                t(4,i,m) =    -vn(i  ,m) + 3.0*vn(i+1,m) - 3.0*vn(i+2,m) + vn(i+3,m)
                t(5,i,m) =    -vn(i+1,m) + 3.0*vn(i+2,m) - 3.0*vn(i+3,m) + vn(i+4,m)
            end do
        end do      
    end do            
            
    do m=nst,ned
        do j=1,nsecn
            do i=sec_stn(j),sec_edn(j)
                do n=1,5
                    is = g(n,i,m)*g(n,i,m) + s(n,i,m)*s(n,i,m) + t(n,i,m)*t(n,i,m)              
                    eis2(n) = is**2
                end do
                
                tau7 = abs(eis2(4)-eis2(1))
                
                do n=1,4
                    bl(n) = cl(n)*(1.0+tau7/(eis2(n)+eps))
                end do
                
                bs = bl(1) + bl(2) + bl(3) + bl(4)
                do n=1,4
                    wl(n,i,m) = bl(n)/bs
                end do
          
                tau7 = abs(eis2(5)-eis2(2))           
                do n=2,5
                    br(n) = cr(n)*(1.0+tau7/(eis2(n)+eps))
                end do

                bs = br(2) + br(3) + br(4) + br(5)
                do n=2,5
                    wr(n,i,m) = br(n)/bs
                end do
            end do
        end do
    end do
    
!!    open(2,file='sigma.dat',form='formatted')
!!    write(2,*) 'variables="i" "sigmal" "sigmar"'
!!    do i=st,ed
!!
!!        write(2,200) i,sigmal(i,1),sigmar(i,1)
!!
!!    end do
!!
!!200 format(I4,1x,e12.5,1x,e12.5)
!!    close(2)
!!
!!    stop    
    
    do m=nst,ned
        do j=1,nsecn
            do i=sec_stn(j),sec_edn(j)
                vln(i,m) = vn(i  ,m) + ( wl(1,i,m)*( 24.0*g(1,i,m)-6.0*s(1,i,m)+t(1,i,m)) + &
                                         wl(2,i,m)*( 24.0*g(2,i,m)-6.0*s(2,i,m)+t(2,i,m)) + &
                                         wl(3,i,m)*( 24.0*g(3,i,m)-6.0*s(3,i,m)+t(3,i,m)) + &
                                         wl(4,i,m)*( 24.0*g(4,i,m)-6.0*s(4,i,m)+t(4,i,m)) )/48.0

                vrn(i,m) = vn(i+1,m) + ( wr(2,i,m)*(-24.0*g(2,i,m)-6.0*s(2,i,m)-t(2,i,m)) + &
                                         wr(3,i,m)*(-24.0*g(3,i,m)-6.0*s(3,i,m)-t(3,i,m)) + &
                                         wr(4,i,m)*(-24.0*g(4,i,m)-6.0*s(4,i,m)-t(4,i,m)) + &
                                         wr(5,i,m)*(-24.0*g(5,i,m)-6.0*s(5,i,m)-t(5,i,m)) )/48.0
            end do
        end do
    end do    
      
    do j=1,nsec
        call d5_via_T2N(ni,nst,ned,ngn,vn,nfs,nfe,nge,d5_li,sec_st(j),sec_ed(j))
        
        call vec2_via_T4N(ni,nst,ned,ngn,vn,nfs,nfe,nge,d2_li,sec_st(j),sec_ed(j))
        
        call vec4_via_T2N(ni,nst,ned,ngn,vn,nfs,nfe,nge,d4_li,sec_st(j),sec_ed(j))
    end do
    
    do m=nst,ned
        do j=1,nsec
            do i=sec_st(j),sec_ed(j)
                
                ve_lil(i,m) = 0.5*(vn(i+1,m) + vn(i  ,m)) - 1.0/8.0*d2_li(i,m) - 1.0/384.0*d4_li(i,m) -  aa_scs*d5_li(i,m)
                            
                ve_lir(i,m) = 0.5*(vn(i+1,m) + vn(i  ,m)) - 1.0/8.0*d2_li(i,m) - 1.0/384.0*d4_li(i,m) +  aa_scs*d5_li(i,m)

            end do
        end do
    end do
        
    do m=nst,ned
        do i=st-1,ed-1
            
               vl(i+1,m) = sigma(i)*ve_lil(i,m) + vln(i,m)*(1.0-sigma(i))
                         
               vr(i+1,m) = sigma(i)*ve_lir(i,m) + vrn(i,m)*(1.0-sigma(i))            

        end do
    end do  

    do i=stc,edc
        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vl(i,m) < chklim(m) ) then
                    vl(i,m) = vn(i  ,m)
                    ndec = ndec + 1
                end if

                if (vr(i,m) < chklim(m) ) then
                    vr(i,m) = vn(i+1,m)
                    ndec = ndec + 1
                end if
           end if
        end do
    end do

end subroutine scsh5picc

subroutine dcsh7pi(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_bound_scheme
    implicit none
    integer(kind_int), intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  , intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(in)  :: sn(1-ngn:ni+ngn,1:3)
    integer(kind_int), intent(in)  :: nfs,nfe
    integer(kind_int), intent(in)  :: nge
    real(kind_real)  , intent(out) :: vl(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(out) :: vr(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , external    :: sub_limit
    integer(kind_int), intent(in)  :: chkid(nst:ned)
    real(kind_real),   intent(in)  :: chklim(nst:ned)
    integer(kind_int), intent(out) :: ndec  
    real(kind_real)   :: vln(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: vrn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: ve_lil(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: ve_lir(1-ngn:ni+ngn,nst:ned)    
    real(kind_real)   :: s(1:5,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: g(1:5,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: t(1:5,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: wl(1:5,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: wr(1:5,1-ngn:ni+ngn,nst:ned)    
    real(kind_real)   :: aa,aa_scs
    real(kind_real)   :: sigma(1-ngn:ni+ngn),imax,imin,indicator
    real(kind_real), parameter :: cl(5) = (/1.D0/64.D0,21.D0/64.D0,35.D0/64.D0,7.D0/64.D0,0.D0/64.D0/)
    real(kind_real), parameter :: cr(5) = (/0.D0/64.D0,7.D0/64.D0,35.D0/64.D0,21.D0/64.D0,1.D0/64.D0/)  
    real(kind_real), parameter :: eps        = 1.0e-30
    real(kind_real), parameter :: indicatorc = 0.1
    real(kind_real)            :: is,eis2(5),bs,tau7,bl(5),br(5)
    integer(kind_int)          :: i,m,n,j,st,ed,stc,edc
    integer(kind_int)          :: id,YN,ind(1:ni+2*ngn),ind_st,ind_ed,pre_ed
    integer(kind_int)          :: nsec,sec_st(1:ni+2*ngn),sec_ed(1:ni+2*ngn)
    integer(kind_int)          :: nsecn,sec_stn(1:ni+2*ngn),sec_edn(1:ni+2*ngn)
    
    if (ni < 8) then
        call muscl2pv(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)
        return
    end if
    
    ndec = 0
    
    aa     = 0.01
    aa_scs = 0.001
    
    wl(:,:,:) = 0.25
    wr(:,:,:) = 0.25
    vln(:,:)    = 0.0
    ve_lil(:,:) = 0.0
    vrn(:,:)    = 0.0
    ve_lir(:,:) = 0.0    
           
    if (nfs == nbc_bound_scheme) then
        
        st  = 4
        stc = 0

        do m=nst,ned
            vl(0,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                  ! third-order
            vl(1,m) = (  3.0*vn(1,m) +  6.0*vn(2,m) -      vn(3,m) ) / 8.0   ! third-order 
            vl(2,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 + &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )
            vr(0,m) = vl(0,m)
            vr(1,m) = vl(1,m)
            vr(2,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 - &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )
            vl(3,m) =   ( 150.0*( vn(4,m) + vn(3,m) )          -  &
                           25.0*( vn(5,m) + vn(2,m) )          +  &
                            3.0*( vn(6,m) + vn(1,m) ) ) /256.0 -  &
                      aa*( 10.0*( vn(4,m) - vn(3,m) )          -  &
                            5.0*( vn(5,m) - vn(2,m) )          +  &
                                ( vn(6,m) - vn(1,m) ) ) 
            
            vr(3,m) =   ( 150.0*( vn(4,m) + vn(3,m) )          -  &
                           25.0*( vn(5,m) + vn(2,m) )          +  &
                            3.0*( vn(6,m) + vn(1,m) ) ) /256.0 +  &
                      aa*( 10.0*( vn(4,m) - vn(3,m) )          -  &
                            5.0*( vn(5,m) - vn(2,m) )          +  &
                                ( vn(6,m) - vn(1,m) ) )          
        end do
        
    else
        
        st  = -2
        stc = st
     
    end if
     
    if (nfe == nbc_bound_scheme) then
        
        ed  = ni - 4
        edc = ni

        do m=nst,ned
            vr(ni  ,m) = ( 15.0*vn(ni,m) - 10.0*vn(ni-1,m) +  3.0*vn(ni-2,m) ) / 8.0                ! third-order
            vl(ni-1,m) = (  3.0*vn(ni,m) +  6.0*vn(ni-1,m) -      vn(ni-2,m) ) / 8.0   ! third-order
            vl(ni-2,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 + &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) )
            vl(ni  ,m) = vr(ni  ,m)
            vr(ni-1,m) = vl(ni-1,m)
            vr(ni-2,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 - &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) ) 
            vl(ni-3,m) =   ( 150.0*( vn(ni-2,m) + vn(ni-3,m) )          -  &
                              25.0*( vn(ni-1,m) + vn(ni-4,m) )          +  &
                               3.0*( vn(ni  ,m) + vn(ni-5,m) ) ) /256.0 -  &
                         aa*( 10.0*( vn(ni-2,m) - vn(ni-3,m) )          -  &
                               5.0*( vn(ni-1,m) - vn(ni-4,m) )          +  &
                                   ( vn(ni  ,m) - vn(ni-5,m) ) ) 
            
            vr(ni-3,m) =   ( 150.0*( vn(ni-2,m) + vn(ni-3,m) )          -  &
                              25.0*( vn(ni-1,m) + vn(ni-4,m) )          +  &
                               3.0*( vn(ni  ,m) + vn(ni-5,m) ) ) /256.0 +  &
                         aa*( 10.0*( vn(ni-2,m) - vn(ni-3,m) )          -  &
                               5.0*( vn(ni-1,m) - vn(ni-4,m) )          +  &
                                   ( vn(ni  ,m) - vn(ni-5,m) ) )              
        end do
        
    else
        
        ed  = ni + 2
        edc = ed
            
    end if
    
    do m=nst,ned
        do i=st-2,ed+2
            g(1,i,m) = vn(i  ,m) - vn(i-1,m)            
            g(2,i,m) = vn(i+1,m) - vn(i  ,m)
            g(3,i,m) = vn(i+2,m) - vn(i+1,m)            
        end do
    end do   
    
    do m=nst,ned
        do i=st-2,ed+2
            do n=1,3
                is = g(n,i,m)*g(n,i,m)
                eis2(n) = is**2
            end do

            tau7 = abs(eis2(2)-eis2(1))
            bl(1) = 1.0+tau7/(eis2(1)+eps)
            bl(2) = 1.0+tau7/(eis2(2)+eps)
            bs = bl(1) + bl(2)
            wl(2,i,m) = bl(2)/bs
            
            tau7 = abs(eis2(3)-eis2(2))
            br(2) = 1.0+tau7/(eis2(2)+eps)
            br(3) = 1.0+tau7/(eis2(3)+eps)
            bs = br(2) + br(3)
            wr(2,i,m) = br(2)/bs
        end do
    end do
    
    do i=st,ed
        imin=1.0
        imax=0.0
        do m=nst,ned
            imin         = min(wl(2,i-3,m),wl(2,i-2,m),wl(2,i-1,m),wl(2,i,m),wl(2,i+1,m),wl(2,i+2,m),wl(2,i+3,m), &
                               wr(2,i-3,m),wr(2,i-2,m),wr(2,i-1,m),wr(2,i,m),wr(2,i+1,m),wr(2,i+2,m),wr(2,i+3,m),imin)                                
            imax         = max(wl(2,i-3,m),wl(2,i-2,m),wl(2,i-1,m),wl(2,i,m),wl(2,i+1,m),wl(2,i+2,m),wl(2,i+3,m), &
                               wr(2,i-3,m),wr(2,i-2,m),wr(2,i-1,m),wr(2,i,m),wr(2,i+1,m),wr(2,i+2,m),wr(2,i+3,m),imax)
        end do
        indicator    = imin/imax
        if(indicator>indicatorc) then
            sigma(i)=1.0
        else
            sigma(i)=0.0
        endif        
    end do
    
    id=0
    do i=st,ed
        YN=0
        if(sigma(i)<0.5) then
            YN=YN+1
        end if
        if(YN>0) then
            id=id+1
            ind(id)=i           
        end if
    end do
       
    do i=1-ngn,ni+ngn
        sigma(i)=0.0
    end do
    if(ed-st<9) then
        ind_st=st
        ind_ed=ed
        nsec=0
        
        nsecn=1
        sec_stn(1)=st
        sec_edn(1)=ed
    else
        if(id==0) then
            nsec=1
            sec_st(1)=st
            sec_ed(1)=ed
            do n=sec_st(1),sec_ed(1)
                sigma(n)=1.0
            end do
            ind_st=st
            ind_ed=ed
            nsecn=0
        else
            ind_st=st
            ind_ed=ed
            nsec=0
            nsecn=0
            if(ind(1)-ind_st>9) then
                nsec=1
                sec_st(1)=ind_st
                sec_ed(1)=ind(1)-1
                pre_ed   =sec_ed(1)
                do n=sec_st(1),sec_ed(1)
                    sigma(n)=1.0
                end do
            else
                nsecn=1
                sec_stn(1)=ind_st
                sec_edn(1)=ind(1)
                pre_ed    =sec_edn(1)
            end if
            do i=2,id
                if(ind(i)-ind(i-1)>9) then
                    nsec=nsec+1
                    sec_st(nsec)=ind(i-1)+1
                    sec_ed(nsec)=ind(i)-1
                    do n=sec_st(nsec),sec_ed(nsec)
                        sigma(n)=1.0
                    end do
                    nsecn=nsecn+1
                    sec_stn(nsecn)=pre_ed+1
                    sec_edn(nsecn)=ind(i-1)
                    pre_ed        =sec_ed(nsec)
                else
                    if(ind(i)-ind(i-1)>1 .or. i==id) then
                        nsecn=nsecn+1
                        sec_stn(nsecn)=pre_ed+1
                        sec_edn(nsecn)=ind(i)
                        pre_ed        =sec_edn(nsecn)
                    end if
                end if
            end do
            if(ind_ed-ind(id)>9) then
                nsec=nsec+1
                sec_st(nsec)=pre_ed+1
                sec_ed(nsec)=ind_ed
                do n=sec_st(nsec),sec_ed(nsec)
                    sigma(n)=1.0
                end do
            else
                nsecn=nsecn+1
                sec_stn(nsecn)=pre_ed+1
                sec_edn(nsecn)=ind_ed
            end if
        end if
    end if
    
    do m=nst,ned
        do j=1,nsecn
            do i=sec_stn(j),sec_edn(j)
                g(1,i,m) = -23.0/24.0*vn(i-3,m) + 31.0/8.0*vn(i-2,m) - 47.0/8.0*vn(i-1,m) + 71.0/24.0*vn(i  ,m)
                g(2,i,m) =   1.0/24.0*vn(i-2,m) -  1.0/8.0*vn(i-1,m) -  7.0/8.0*vn(i  ,m) + 23.0/24.0*vn(i+1,m)
                g(3,i,m) =   1.0/24.0*vn(i-1,m) -  9.0/8.0*vn(i  ,m) +  9.0/8.0*vn(i+1,m) -  1.0/24.0*vn(i+2,m)
                g(4,i,m) = -23.0/24.0*vn(i  ,m) +  7.0/8.0*vn(i+1,m) +  1.0/8.0*vn(i+2,m) -  1.0/24.0*vn(i+3,m)
                g(5,i,m) = -71.0/24.0*vn(i+1,m) + 47.0/8.0*vn(i+2,m) - 31.0/8.0*vn(i+3,m) + 23.0/24.0*vn(i+4,m)
                
                s(1,i,m) =   -3.0/2.0*vn(i-3,m) + 11.0/2.0*vn(i-2,m) - 13.0/2.0*vn(i-1,m) +   5.0/2.0*vn(i  ,m)
                s(2,i,m) =   -1.0/2.0*vn(i-2,m) +  5.0/2.0*vn(i-1,m) -  7.0/2.0*vn(i  ,m) +   3.0/2.0*vn(i+1,m)
                s(3,i,m) =    1.0/2.0*vn(i-1,m) -  1.0/2.0*vn(i  ,m) -  1.0/2.0*vn(i+1,m) +   1.0/2.0*vn(i+2,m)
                s(4,i,m) =    3.0/2.0*vn(i  ,m) -  7.0/2.0*vn(i+1,m) +  5.0/2.0*vn(i+2,m) -   1.0/2.0*vn(i+3,m)
                s(5,i,m) =    5.0/2.0*vn(i+1,m) - 13.0/2.0*vn(i+2,m) + 11.0/2.0*vn(i+3,m) -   3.0/2.0*vn(i+4,m)
                
                t(1,i,m) =    -vn(i-3,m) + 3.0*vn(i-2,m) - 3.0*vn(i-1,m) + vn(i  ,m)
                t(2,i,m) =    -vn(i-2,m) + 3.0*vn(i-1,m) - 3.0*vn(i  ,m) + vn(i+1,m)
                t(3,i,m) =    -vn(i-1,m) + 3.0*vn(i  ,m) - 3.0*vn(i+1,m) + vn(i+2,m)
                t(4,i,m) =    -vn(i  ,m) + 3.0*vn(i+1,m) - 3.0*vn(i+2,m) + vn(i+3,m)
                t(5,i,m) =    -vn(i+1,m) + 3.0*vn(i+2,m) - 3.0*vn(i+3,m) + vn(i+4,m)
            end do
        end do      
    end do            
            
    do m=nst,ned
        do j=1,nsecn
            do i=sec_stn(j),sec_edn(j)
                do n=1,5
                    is = g(n,i,m)*g(n,i,m) + s(n,i,m)*s(n,i,m) + t(n,i,m)*t(n,i,m)              
                    eis2(n) = is**2
                end do
                
                tau7 = abs(eis2(4)-eis2(1))
                
                do n=1,4
                    bl(n) = cl(n)*(1.0+tau7/(eis2(n)+eps))
                end do
                
                bs = bl(1) + bl(2) + bl(3) + bl(4)
                do n=1,4
                    wl(n,i,m) = bl(n)/bs
                end do
          
                tau7 = abs(eis2(5)-eis2(2))           
                do n=2,5
                    br(n) = cr(n)*(1.0+tau7/(eis2(n)+eps))
                end do

                bs = br(2) + br(3) + br(4) + br(5)
                do n=2,5
                    wr(n,i,m) = br(n)/bs
                end do
            end do
        end do
    end do
    
!!    open(2,file='sigma.dat',form='formatted')
!!    write(2,*) 'variables="i" "sigmal" "sigmar"'
!!    do i=st,ed
!!
!!        write(2,200) i,sigmal(i,1),sigmar(i,1)
!!
!!    end do
!!
!!200 format(I4,1x,e12.5,1x,e12.5)
!!    close(2)
!!
!!    stop    
    
    do m=nst,ned
        do j=1,nsecn
            do i=sec_stn(j),sec_edn(j)
                vln(i,m) = vn(i  ,m) + ( wl(1,i,m)*( 24.0*g(1,i,m)-6.0*s(1,i,m)+t(1,i,m)) + &
                                         wl(2,i,m)*( 24.0*g(2,i,m)-6.0*s(2,i,m)+t(2,i,m)) + &
                                         wl(3,i,m)*( 24.0*g(3,i,m)-6.0*s(3,i,m)+t(3,i,m)) + &
                                         wl(4,i,m)*( 24.0*g(4,i,m)-6.0*s(4,i,m)+t(4,i,m)) )/48.0

                vrn(i,m) = vn(i+1,m) + ( wr(2,i,m)*(-24.0*g(2,i,m)-6.0*s(2,i,m)-t(2,i,m)) + &
                                         wr(3,i,m)*(-24.0*g(3,i,m)-6.0*s(3,i,m)-t(3,i,m)) + &
                                         wr(4,i,m)*(-24.0*g(4,i,m)-6.0*s(4,i,m)-t(4,i,m)) + &
                                         wr(5,i,m)*(-24.0*g(5,i,m)-6.0*s(5,i,m)-t(5,i,m)) )/48.0
            end do
        end do
    end do
      
    do j=1,nsec
        
        call hdcs7ciN(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,ve_lil,ve_lir,sub_limit,chkid,chklim,ndec,sec_st(j),sec_ed(j))        
        
    end do
        
    do m=nst,ned
        do i=st,ed
            
               vl(i,m) = sigma(i)*ve_lil(i,m) + vln(i,m)*(1.0-sigma(i))
                         
               vr(i,m) = sigma(i)*ve_lir(i,m) + vrn(i,m)*(1.0-sigma(i))            

        end do
    end do  

    do i=stc,edc
        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vl(i,m) < chklim(m) ) then
                    vl(i,m) = vn(i  ,m)
                    ndec = ndec + 1
                end if

                if (vr(i,m) < chklim(m) ) then
                    vr(i,m) = vn(i+1,m)
                    ndec = ndec + 1
                end if
           end if
        end do
    end do

end subroutine dcsh7pi

subroutine dcsh7picc(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_bound_scheme
    implicit none
    integer(kind_int), intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  , intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(in)  :: sn(1-ngn:ni+ngn,1:3)
    integer(kind_int), intent(in)  :: nfs,nfe
    integer(kind_int), intent(in)  :: nge
    real(kind_real)  , intent(out) :: vl(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(out) :: vr(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , external    :: sub_limit
    integer(kind_int), intent(in)  :: chkid(nst:ned)
    real(kind_real),   intent(in)  :: chklim(nst:ned)
    integer(kind_int), intent(out) :: ndec  
    real(kind_real)   :: vln(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: vrn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: ve_lil(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: ve_lir(1-ngn:ni+ngn,nst:ned)    
    real(kind_real)   :: s(1:5,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: g(1:5,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: t(1:5,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: wl(1:5,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: wr(1:5,1-ngn:ni+ngn,nst:ned)    
    real(kind_real)   :: aa,aa_scs
    real(kind_real)   :: sigma(1-ngn:ni+ngn),imax,imin,indicator
    real(kind_real), parameter :: cl(5) = (/1.D0/64.D0,21.D0/64.D0,35.D0/64.D0,7.D0/64.D0,0.D0/64.D0/)
    real(kind_real), parameter :: cr(5) = (/0.D0/64.D0,7.D0/64.D0,35.D0/64.D0,21.D0/64.D0,1.D0/64.D0/)  
    real(kind_real), parameter :: eps        = 1.0e-30
    real(kind_real), parameter :: indicatorc = 0.1
    real(kind_real)            :: is,eis2(5),bs,tau7,bl(5),br(5)
    integer(kind_int)          :: i,m,n,j,st,ed,stc,edc
    integer(kind_int)          :: id,YN,ind(1:ni+2*ngn),ind_st,ind_ed,pre_ed
    integer(kind_int)          :: nsec,sec_st(1:ni+2*ngn),sec_ed(1:ni+2*ngn)
    integer(kind_int)          :: nsecn,sec_stn(1:ni+2*ngn),sec_edn(1:ni+2*ngn)
    
    if (ni < 8) then
        call muscl2pvcc(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)
        return
    end if
    
    ndec = 0
    
    aa     = 0.01
    aa_scs = 0.001
    
    wl(:,:,:) = 0.25
    wr(:,:,:) = 0.25
    vln(:,:)    = 0.0
    ve_lil(:,:) = 0.0
    vrn(:,:)    = 0.0
    ve_lir(:,:) = 0.0    
           
    if (nfs == nbc_bound_scheme) then
        
        st  = 5
        stc = 1

        do m=nst,ned
            vl(1,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                  ! third-order
            vl(2,m) = (  3.0*vn(1,m) +  6.0*vn(2,m) -      vn(3,m) ) / 8.0   ! third-order 
            vl(3,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 + &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )
            vr(1,m) = vl(1,m)
            vr(2,m) = vl(2,m)
            vr(3,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 - &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )
            vl(4,m) =   ( 150.0*( vn(4,m) + vn(3,m) )          -  &
                           25.0*( vn(5,m) + vn(2,m) )          +  &
                            3.0*( vn(6,m) + vn(1,m) ) ) /256.0 -  &
                      aa*( 10.0*( vn(4,m) - vn(3,m) )          -  &
                            5.0*( vn(5,m) - vn(2,m) )          +  &
                                ( vn(6,m) - vn(1,m) ) ) 
            
            vr(4,m) =   ( 150.0*( vn(4,m) + vn(3,m) )          -  &
                           25.0*( vn(5,m) + vn(2,m) )          +  &
                            3.0*( vn(6,m) + vn(1,m) ) ) /256.0 +  &
                      aa*( 10.0*( vn(4,m) - vn(3,m) )          -  &
                            5.0*( vn(5,m) - vn(2,m) )          +  &
                                ( vn(6,m) - vn(1,m) ) )          
        end do
        
    else
        
        st  = -1
        stc = st
     
    end if
     
    if (nfe == nbc_bound_scheme) then
        
        ed  = ni - 3
        edc = ni + 1

        do m=nst,ned
            vr(ni+1,m) = ( 15.0*vn(ni,m) - 10.0*vn(ni-1,m) +  3.0*vn(ni-2,m) ) / 8.0                ! third-order
            vl(ni  ,m) = (  3.0*vn(ni,m) +  6.0*vn(ni-1,m) -      vn(ni-2,m) ) / 8.0   ! third-order
            vl(ni-1,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 + &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) )
            vl(ni+1,m) = vr(ni+1,m)
            vr(ni  ,m) = vl(ni  ,m)
            vr(ni-1,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 - &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) ) 
            vl(ni-2,m) =   ( 150.0*( vn(ni-2,m) + vn(ni-3,m) )          -  &
                              25.0*( vn(ni-1,m) + vn(ni-4,m) )          +  &
                               3.0*( vn(ni  ,m) + vn(ni-5,m) ) ) /256.0 -  &
                         aa*( 10.0*( vn(ni-2,m) - vn(ni-3,m) )          -  &
                               5.0*( vn(ni-1,m) - vn(ni-4,m) )          +  &
                                   ( vn(ni  ,m) - vn(ni-5,m) ) ) 
            
            vr(ni-2,m) =   ( 150.0*( vn(ni-2,m) + vn(ni-3,m) )          -  &
                              25.0*( vn(ni-1,m) + vn(ni-4,m) )          +  &
                               3.0*( vn(ni  ,m) + vn(ni-5,m) ) ) /256.0 +  &
                         aa*( 10.0*( vn(ni-2,m) - vn(ni-3,m) )          -  &
                               5.0*( vn(ni-1,m) - vn(ni-4,m) )          +  &
                                   ( vn(ni  ,m) - vn(ni-5,m) ) )              
        end do
        
    else
        
        ed  = ni + 3
        edc = ed
            
    end if
    
    do m=nst,ned
        do i=st-3,ed+1
            g(1,i,m) = vn(i  ,m) - vn(i-1,m)            
            g(2,i,m) = vn(i+1,m) - vn(i  ,m)
            g(3,i,m) = vn(i+2,m) - vn(i+1,m)            
        end do
    end do   
    
    do m=nst,ned
        do i=st-3,ed+1
            do n=1,3
                is = g(n,i,m)*g(n,i,m)
                eis2(n) = is**2
            end do

            tau7 = abs(eis2(2)-eis2(1))
            bl(1) = 1.0+tau7/(eis2(1)+eps)
            bl(2) = 1.0+tau7/(eis2(2)+eps)
            bs = bl(1) + bl(2)
            wl(2,i,m) = bl(2)/bs
            
            tau7 = abs(eis2(3)-eis2(2))
            br(2) = 1.0+tau7/(eis2(2)+eps)
            br(3) = 1.0+tau7/(eis2(3)+eps)
            bs = br(2) + br(3)
            wr(2,i,m) = br(2)/bs
        end do
    end do
    
    do i=st-1,ed-1
        imin=1.0
        imax=0.0
        do m=nst,ned
            imin         = min(wl(2,i-3,m),wl(2,i-2,m),wl(2,i-1,m),wl(2,i,m),wl(2,i+1,m),wl(2,i+2,m),wl(2,i+3,m), &
                               wr(2,i-3,m),wr(2,i-2,m),wr(2,i-1,m),wr(2,i,m),wr(2,i+1,m),wr(2,i+2,m),wr(2,i+3,m),imin)                                
            imax         = max(wl(2,i-3,m),wl(2,i-2,m),wl(2,i-1,m),wl(2,i,m),wl(2,i+1,m),wl(2,i+2,m),wl(2,i+3,m), &
                               wr(2,i-3,m),wr(2,i-2,m),wr(2,i-1,m),wr(2,i,m),wr(2,i+1,m),wr(2,i+2,m),wr(2,i+3,m),imax)
        end do
        indicator    = imin/imax
        if(indicator>indicatorc) then
            sigma(i)=1.0
        else
            sigma(i)=0.0
        endif        
    end do
    
    id=0
    do i=st-1,ed-1
        YN=0
        if(sigma(i)<0.5) then
            YN=YN+1
        end if
        if(YN>0) then
            id=id+1
            ind(id)=i           
        end if
    end do
       
    do i=1-ngn,ni+ngn
        sigma(i)=0.0
    end do
    if(ed-st<9) then
        ind_st=st-1
        ind_ed=ed-1
        nsec=0
        
        nsecn=1
        sec_stn(1)=st-1
        sec_edn(1)=ed-1
    else
        if(id==0) then
            nsec=1
            sec_st(1)=st-1
            sec_ed(1)=ed-1
            do n=sec_st(1),sec_ed(1)
                sigma(n)=1.0
            end do
            ind_st=st-1
            ind_ed=ed-1
            nsecn=0
        else
            ind_st=st-1
            ind_ed=ed-1
            nsec=0
            nsecn=0
            if(ind(1)-ind_st>9) then
                nsec=1
                sec_st(1)=ind_st
                sec_ed(1)=ind(1)-1
                pre_ed   =sec_ed(1)
                do n=sec_st(1),sec_ed(1)
                    sigma(n)=1.0
                end do
            else
                nsecn=1
                sec_stn(1)=ind_st
                sec_edn(1)=ind(1)
                pre_ed    =sec_edn(1)
            end if
            do i=2,id
                if(ind(i)-ind(i-1)>9) then
                    nsec=nsec+1
                    sec_st(nsec)=ind(i-1)+1
                    sec_ed(nsec)=ind(i)-1
                    do n=sec_st(nsec),sec_ed(nsec)
                        sigma(n)=1.0
                    end do
                    nsecn=nsecn+1
                    sec_stn(nsecn)=pre_ed+1
                    sec_edn(nsecn)=ind(i-1)
                    pre_ed        =sec_ed(nsec)
                else
                    if(ind(i)-ind(i-1)>1 .or. i==id) then
                        nsecn=nsecn+1
                        sec_stn(nsecn)=pre_ed+1
                        sec_edn(nsecn)=ind(i)
                        pre_ed        =sec_edn(nsecn)
                    end if
                end if
            end do
            if(ind_ed-ind(id)>9) then
                nsec=nsec+1
                sec_st(nsec)=pre_ed+1
                sec_ed(nsec)=ind_ed
                do n=sec_st(nsec),sec_ed(nsec)
                    sigma(n)=1.0
                end do
            else
                nsecn=nsecn+1
                sec_stn(nsecn)=pre_ed+1
                sec_edn(nsecn)=ind_ed
            end if
        end if
    end if
    
    do m=nst,ned
        do j=1,nsecn
            do i=sec_stn(j),sec_edn(j)
                g(1,i,m) = -23.0/24.0*vn(i-3,m) + 31.0/8.0*vn(i-2,m) - 47.0/8.0*vn(i-1,m) + 71.0/24.0*vn(i  ,m)
                g(2,i,m) =   1.0/24.0*vn(i-2,m) -  1.0/8.0*vn(i-1,m) -  7.0/8.0*vn(i  ,m) + 23.0/24.0*vn(i+1,m)
                g(3,i,m) =   1.0/24.0*vn(i-1,m) -  9.0/8.0*vn(i  ,m) +  9.0/8.0*vn(i+1,m) -  1.0/24.0*vn(i+2,m)
                g(4,i,m) = -23.0/24.0*vn(i  ,m) +  7.0/8.0*vn(i+1,m) +  1.0/8.0*vn(i+2,m) -  1.0/24.0*vn(i+3,m)
                g(5,i,m) = -71.0/24.0*vn(i+1,m) + 47.0/8.0*vn(i+2,m) - 31.0/8.0*vn(i+3,m) + 23.0/24.0*vn(i+4,m)
                
                s(1,i,m) =   -3.0/2.0*vn(i-3,m) + 11.0/2.0*vn(i-2,m) - 13.0/2.0*vn(i-1,m) +   5.0/2.0*vn(i  ,m)
                s(2,i,m) =   -1.0/2.0*vn(i-2,m) +  5.0/2.0*vn(i-1,m) -  7.0/2.0*vn(i  ,m) +   3.0/2.0*vn(i+1,m)
                s(3,i,m) =    1.0/2.0*vn(i-1,m) -  1.0/2.0*vn(i  ,m) -  1.0/2.0*vn(i+1,m) +   1.0/2.0*vn(i+2,m)
                s(4,i,m) =    3.0/2.0*vn(i  ,m) -  7.0/2.0*vn(i+1,m) +  5.0/2.0*vn(i+2,m) -   1.0/2.0*vn(i+3,m)
                s(5,i,m) =    5.0/2.0*vn(i+1,m) - 13.0/2.0*vn(i+2,m) + 11.0/2.0*vn(i+3,m) -   3.0/2.0*vn(i+4,m)
                
                t(1,i,m) =    -vn(i-3,m) + 3.0*vn(i-2,m) - 3.0*vn(i-1,m) + vn(i  ,m)
                t(2,i,m) =    -vn(i-2,m) + 3.0*vn(i-1,m) - 3.0*vn(i  ,m) + vn(i+1,m)
                t(3,i,m) =    -vn(i-1,m) + 3.0*vn(i  ,m) - 3.0*vn(i+1,m) + vn(i+2,m)
                t(4,i,m) =    -vn(i  ,m) + 3.0*vn(i+1,m) - 3.0*vn(i+2,m) + vn(i+3,m)
                t(5,i,m) =    -vn(i+1,m) + 3.0*vn(i+2,m) - 3.0*vn(i+3,m) + vn(i+4,m)
            end do
        end do      
    end do            
            
    do m=nst,ned
        do j=1,nsecn
            do i=sec_stn(j),sec_edn(j)
                do n=1,5
                    is = g(n,i,m)*g(n,i,m) + s(n,i,m)*s(n,i,m) + t(n,i,m)*t(n,i,m)              
                    eis2(n) = is**2
                end do
                
                tau7 = abs(eis2(4)-eis2(1))
                
                do n=1,4
                    bl(n) = cl(n)*(1.0+tau7/(eis2(n)+eps))
                end do
                
                bs = bl(1) + bl(2) + bl(3) + bl(4)
                do n=1,4
                    wl(n,i,m) = bl(n)/bs
                end do
          
                tau7 = abs(eis2(5)-eis2(2))           
                do n=2,5
                    br(n) = cr(n)*(1.0+tau7/(eis2(n)+eps))
                end do

                bs = br(2) + br(3) + br(4) + br(5)
                do n=2,5
                    wr(n,i,m) = br(n)/bs
                end do
            end do
        end do
    end do
    
!!    open(2,file='sigma.dat',form='formatted')
!!    write(2,*) 'variables="i" "sigmal" "sigmar"'
!!    do i=st,ed
!!
!!        write(2,200) i,sigmal(i,1),sigmar(i,1)
!!
!!    end do
!!
!!200 format(I4,1x,e12.5,1x,e12.5)
!!    close(2)
!!
!!    stop    
    
    do m=nst,ned
        do j=1,nsecn
            do i=sec_stn(j),sec_edn(j)
                vln(i,m) = vn(i  ,m) + ( wl(1,i,m)*( 24.0*g(1,i,m)-6.0*s(1,i,m)+t(1,i,m)) + &
                                         wl(2,i,m)*( 24.0*g(2,i,m)-6.0*s(2,i,m)+t(2,i,m)) + &
                                         wl(3,i,m)*( 24.0*g(3,i,m)-6.0*s(3,i,m)+t(3,i,m)) + &
                                         wl(4,i,m)*( 24.0*g(4,i,m)-6.0*s(4,i,m)+t(4,i,m)) )/48.0

                vrn(i,m) = vn(i+1,m) + ( wr(2,i,m)*(-24.0*g(2,i,m)-6.0*s(2,i,m)-t(2,i,m)) + &
                                         wr(3,i,m)*(-24.0*g(3,i,m)-6.0*s(3,i,m)-t(3,i,m)) + &
                                         wr(4,i,m)*(-24.0*g(4,i,m)-6.0*s(4,i,m)-t(4,i,m)) + &
                                         wr(5,i,m)*(-24.0*g(5,i,m)-6.0*s(5,i,m)-t(5,i,m)) )/48.0
            end do
        end do
    end do
      
    do j=1,nsec
        
        call hdcs7ciN(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,ve_lil,ve_lir,sub_limit,chkid,chklim,ndec,sec_st(j),sec_ed(j))        
        
    end do
        
    do m=nst,ned
        do i=st-1,ed-1
            
               vl(i+1,m) = sigma(i)*ve_lil(i,m) + vln(i,m)*(1.0-sigma(i))
                         
               vr(i+1,m) = sigma(i)*ve_lir(i,m) + vrn(i,m)*(1.0-sigma(i))            

        end do
    end do  

    do i=stc,edc
        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vl(i,m) < chklim(m) ) then
                    vl(i,m) = vn(i  ,m)
                    ndec = ndec + 1
                end if

                if (vr(i,m) < chklim(m) ) then
                    vr(i,m) = vn(i+1,m)
                    ndec = ndec + 1
                end if
           end if
        end do
    end do

end subroutine dcsh7picc

subroutine scsn4pi(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_bound_scheme
    implicit none
    integer(kind_int), intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  , intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(in)  :: sn(1-ngn:ni+ngn,1:3)
    integer(kind_int), intent(in)  :: nfs,nfe
    integer(kind_int), intent(in)  :: nge
    real(kind_real)  , intent(out) :: vl(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(out) :: vr(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , external    :: sub_limit
    integer(kind_int), intent(in)  :: chkid(nst:ned)
    real(kind_real),   intent(in)  :: chklim(nst:ned)
    integer(kind_int), intent(out) :: ndec  
    real(kind_real)   :: s(1:5,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: g(1:5,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: t(1:5,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: wl(1:5,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: wr(1:5,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: dl3(-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: dr3(-ngn:ni+ngn,nst:ned)    
    real(kind_real)   :: aa
    real(kind_real), parameter :: cl(5) = (/1.D0/64.D0,21.D0/64.D0,35.D0/64.D0,7.D0/64.D0,0.D0/64.D0/)
    real(kind_real), parameter :: cr(5) = (/0.D0/64.D0,7.D0/64.D0,35.D0/64.D0,21.D0/64.D0,1.D0/64.D0/)
    real(kind_real), parameter :: eps        = 1.0e-30
    real(kind_real), parameter :: indicatorc = 0.1
    real(kind_real)            :: is,eis2(5),bs,tau7,bl(5),br(5)
    integer(kind_int)          :: i,m,n,st,ed,stc,edc
    
    if (ni < 8) then
        call muscl2pv(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)
        return
    end if
    
    ndec = 0
    
    aa     = 0.01
    
    wl(:,:,:) = 0.25
    wr(:,:,:) = 0.25  
           
    if (nfs == nbc_bound_scheme) then
        
        st  = 4
        stc = 0

        do m=nst,ned
            vl(0,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                  ! third-order
            vl(1,m) = (  3.0*vn(1,m) +  6.0*vn(2,m) -      vn(3,m) ) / 8.0   ! third-order 
            vl(2,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 + &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )
            vr(0,m) = vl(0,m)
            vr(1,m) = vl(1,m)
            vr(2,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 - &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )
            vl(3,m) =   ( 150.0*( vn(4,m) + vn(3,m) )          -  &
                           25.0*( vn(5,m) + vn(2,m) )          +  &
                            3.0*( vn(6,m) + vn(1,m) ) ) /256.0 -  &
                      aa*( 10.0*( vn(4,m) - vn(3,m) )          -  &
                            5.0*( vn(5,m) - vn(2,m) )          +  &
                                ( vn(6,m) - vn(1,m) ) ) 
            
            vr(3,m) =   ( 150.0*( vn(4,m) + vn(3,m) )          -  &
                           25.0*( vn(5,m) + vn(2,m) )          +  &
                            3.0*( vn(6,m) + vn(1,m) ) ) /256.0 +  &
                      aa*( 10.0*( vn(4,m) - vn(3,m) )          -  &
                            5.0*( vn(5,m) - vn(2,m) )          +  &
                                ( vn(6,m) - vn(1,m) ) )          
        end do
        
    else
        
        st  = -2
        stc = st
     
    end if
     
    if (nfe == nbc_bound_scheme) then
        
        ed  = ni - 4
        edc = ni

        do m=nst,ned
            vr(ni  ,m) = ( 15.0*vn(ni,m) - 10.0*vn(ni-1,m) +  3.0*vn(ni-2,m) ) / 8.0                ! third-order
            vl(ni-1,m) = (  3.0*vn(ni,m) +  6.0*vn(ni-1,m) -      vn(ni-2,m) ) / 8.0   ! third-order
            vl(ni-2,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 + &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) )
            vl(ni  ,m) = vr(ni  ,m)
            vr(ni-1,m) = vl(ni-1,m)
            vr(ni-2,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 - &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) ) 
            vl(ni-3,m) =   ( 150.0*( vn(ni-2,m) + vn(ni-3,m) )          -  &
                              25.0*( vn(ni-1,m) + vn(ni-4,m) )          +  &
                               3.0*( vn(ni  ,m) + vn(ni-5,m) ) ) /256.0 -  &
                         aa*( 10.0*( vn(ni-2,m) - vn(ni-3,m) )          -  &
                               5.0*( vn(ni-1,m) - vn(ni-4,m) )          +  &
                                   ( vn(ni  ,m) - vn(ni-5,m) ) ) 
            
            vr(ni-3,m) =   ( 150.0*( vn(ni-2,m) + vn(ni-3,m) )          -  &
                              25.0*( vn(ni-1,m) + vn(ni-4,m) )          +  &
                               3.0*( vn(ni  ,m) + vn(ni-5,m) ) ) /256.0 +  &
                         aa*( 10.0*( vn(ni-2,m) - vn(ni-3,m) )          -  &
                               5.0*( vn(ni-1,m) - vn(ni-4,m) )          +  &
                                   ( vn(ni  ,m) - vn(ni-5,m) ) )              
        end do
        
    else
        
        ed  = ni + 2
        edc = ed
            
    end if
    
    do m=nst,ned

        do i=st,ed
            g(1,i,m) = -23.0/24.0*vn(i-3,m) + 31.0/8.0*vn(i-2,m) - 47.0/8.0*vn(i-1,m) + 71.0/24.0*vn(i  ,m)
            g(2,i,m) =   1.0/24.0*vn(i-2,m) -  1.0/8.0*vn(i-1,m) -  7.0/8.0*vn(i  ,m) + 23.0/24.0*vn(i+1,m)
            g(3,i,m) =   1.0/24.0*vn(i-1,m) -  9.0/8.0*vn(i  ,m) +  9.0/8.0*vn(i+1,m) -  1.0/24.0*vn(i+2,m)
            g(4,i,m) = -23.0/24.0*vn(i  ,m) +  7.0/8.0*vn(i+1,m) +  1.0/8.0*vn(i+2,m) -  1.0/24.0*vn(i+3,m)
            g(5,i,m) = -71.0/24.0*vn(i+1,m) + 47.0/8.0*vn(i+2,m) - 31.0/8.0*vn(i+3,m) + 23.0/24.0*vn(i+4,m)
            
            s(1,i,m) =   -3.0/2.0*vn(i-3,m) + 11.0/2.0*vn(i-2,m) - 13.0/2.0*vn(i-1,m) +   5.0/2.0*vn(i  ,m)
            s(2,i,m) =   -1.0/2.0*vn(i-2,m) +  5.0/2.0*vn(i-1,m) -  7.0/2.0*vn(i  ,m) +   3.0/2.0*vn(i+1,m)
            s(3,i,m) =    1.0/2.0*vn(i-1,m) -  1.0/2.0*vn(i  ,m) -  1.0/2.0*vn(i+1,m) +   1.0/2.0*vn(i+2,m)
            s(4,i,m) =    3.0/2.0*vn(i  ,m) -  7.0/2.0*vn(i+1,m) +  5.0/2.0*vn(i+2,m) -   1.0/2.0*vn(i+3,m)
            s(5,i,m) =    5.0/2.0*vn(i+1,m) - 13.0/2.0*vn(i+2,m) + 11.0/2.0*vn(i+3,m) -   3.0/2.0*vn(i+4,m)
            
            t(1,i,m) =    -vn(i-3,m) + 3.0*vn(i-2,m) - 3.0*vn(i-1,m) + vn(i  ,m)
            t(2,i,m) =    -vn(i-2,m) + 3.0*vn(i-1,m) - 3.0*vn(i  ,m) + vn(i+1,m)
            t(3,i,m) =    -vn(i-1,m) + 3.0*vn(i  ,m) - 3.0*vn(i+1,m) + vn(i+2,m)
            t(4,i,m) =    -vn(i  ,m) + 3.0*vn(i+1,m) - 3.0*vn(i+2,m) + vn(i+3,m)
            t(5,i,m) =    -vn(i+1,m) + 3.0*vn(i+2,m) - 3.0*vn(i+3,m) + vn(i+4,m)
        end do       
    end do
            
    do m=nst,ned
        do i=st,ed
            do n=1,5
                is = g(n,i,m)*g(n,i,m) + s(n,i,m)*s(n,i,m) + t(n,i,m)*t(n,i,m)              
                eis2(n) = is**2
            end do
            
            tau7 = abs(eis2(4)-eis2(1))
            
            do n=1,4
                bl(n) = cl(n)*(1.0+tau7/(eis2(n)+eps))
            end do
            
            bs = bl(1) + bl(2) + bl(3) + bl(4)
            do n=1,4
                wl(n,i,m) = bl(n)/bs
            end do
          
            tau7 = abs(eis2(5)-eis2(2))           
            do n=2,5
                br(n) = cr(n)*(1.0+tau7/(eis2(n)+eps))
            end do

            bs = br(2) + br(3) + br(4) + br(5)
            do n=2,5
                wr(n,i,m) = br(n)/bs
            end do
        end do
    end do
    
    do m=nst,ned
        do i=st,ed
            vl(i,m) = vn(i  ,m) + ( wl(1,i,m)*( 24.0*g(1,i,m)-6.0*s(1,i,m)+t(1,i,m)) + &
                                    wl(2,i,m)*( 24.0*g(2,i,m)-6.0*s(2,i,m)+t(2,i,m)) + &
                                    wl(3,i,m)*( 24.0*g(3,i,m)-6.0*s(3,i,m)+t(3,i,m)) + &
                                    wl(4,i,m)*( 24.0*g(4,i,m)-6.0*s(4,i,m)+t(4,i,m)) )/48.0
        end do

        do i=st,ed
            vr(i,m) = vn(i+1,m) + ( wr(2,i,m)*(-24.0*g(2,i,m)-6.0*s(2,i,m)-t(2,i,m)) + &
                                    wr(3,i,m)*(-24.0*g(3,i,m)-6.0*s(3,i,m)-t(3,i,m)) + &
                                    wr(4,i,m)*(-24.0*g(4,i,m)-6.0*s(4,i,m)-t(4,i,m)) + &
                                    wr(5,i,m)*(-24.0*g(5,i,m)-6.0*s(5,i,m)-t(5,i,m)) )/48.0
        end do
    end do    

    do i=stc,edc
        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vl(i,m) < chklim(m) ) then
                    vl(i,m) = vn(i  ,m)
                    ndec = ndec + 1
                end if

                if (vr(i,m) < chklim(m) ) then
                    vr(i,m) = vn(i+1,m)
                    ndec = ndec + 1
                end if
           end if
        end do
    end do

end subroutine scsn4pi

subroutine scsn4picc(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_bound_scheme
    implicit none
    integer(kind_int), intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  , intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(in)  :: sn(1-ngn:ni+ngn,1:3)
    integer(kind_int), intent(in)  :: nfs,nfe
    integer(kind_int), intent(in)  :: nge
    real(kind_real)  , intent(out) :: vl(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(out) :: vr(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , external    :: sub_limit
    integer(kind_int), intent(in)  :: chkid(nst:ned)
    real(kind_real),   intent(in)  :: chklim(nst:ned)
    integer(kind_int), intent(out) :: ndec  
    real(kind_real)   :: s(1:5,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: g(1:5,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: t(1:5,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: wl(1:5,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: wr(1:5,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: dl3(-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: dr3(-ngn:ni+ngn,nst:ned)    
    real(kind_real)   :: aa
    real(kind_real), parameter :: cl(5) = (/1.D0/64.D0,21.D0/64.D0,35.D0/64.D0,7.D0/64.D0,0.D0/64.D0/)
    real(kind_real), parameter :: cr(5) = (/0.D0/64.D0,7.D0/64.D0,35.D0/64.D0,21.D0/64.D0,1.D0/64.D0/)
    real(kind_real), parameter :: eps        = 1.0e-30
    real(kind_real), parameter :: indicatorc = 0.1
    real(kind_real)            :: is,eis2(5),bs,tau7,bl(5),br(5)
    integer(kind_int)          :: i,m,n,st,ed,stc,edc
    
    if (ni < 8) then
        call muscl2pvcc(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)
        return
    end if
    
    ndec = 0
    
    aa     = 0.01
    
    wl(:,:,:) = 0.25
    wr(:,:,:) = 0.25  
           
    if (nfs == nbc_bound_scheme) then
        
        st  = 5
        stc = 1

        do m=nst,ned
            vl(1,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                  ! third-order
            vl(2,m) = (  3.0*vn(1,m) +  6.0*vn(2,m) -      vn(3,m) ) / 8.0   ! third-order 
            vl(3,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 + &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )
            vr(1,m) = vl(1,m)
            vr(2,m) = vl(2,m)
            vr(3,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 - &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )
            vl(4,m) =   ( 150.0*( vn(4,m) + vn(3,m) )          -  &
                           25.0*( vn(5,m) + vn(2,m) )          +  &
                            3.0*( vn(6,m) + vn(1,m) ) ) /256.0 -  &
                      aa*( 10.0*( vn(4,m) - vn(3,m) )          -  &
                            5.0*( vn(5,m) - vn(2,m) )          +  &
                                ( vn(6,m) - vn(1,m) ) ) 
            
            vr(4,m) =   ( 150.0*( vn(4,m) + vn(3,m) )          -  &
                           25.0*( vn(5,m) + vn(2,m) )          +  &
                            3.0*( vn(6,m) + vn(1,m) ) ) /256.0 +  &
                      aa*( 10.0*( vn(4,m) - vn(3,m) )          -  &
                            5.0*( vn(5,m) - vn(2,m) )          +  &
                                ( vn(6,m) - vn(1,m) ) )          
        end do
        
    else
        
        st  = -1
        stc = st
     
    end if
     
    if (nfe == nbc_bound_scheme) then
        
        ed  = ni - 3
        edc = ni + 1

        do m=nst,ned
            vr(ni+1,m) = ( 15.0*vn(ni,m) - 10.0*vn(ni-1,m) +  3.0*vn(ni-2,m) ) / 8.0                ! third-order
            vl(ni  ,m) = (  3.0*vn(ni,m) +  6.0*vn(ni-1,m) -      vn(ni-2,m) ) / 8.0   ! third-order
            vl(ni-1,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 + &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) )
            vl(ni+1,m) = vr(ni+1,m)
            vr(ni  ,m) = vl(ni  ,m)
            vr(ni-1,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 - &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) ) 
            vl(ni-2,m) =   ( 150.0*( vn(ni-2,m) + vn(ni-3,m) )          -  &
                              25.0*( vn(ni-1,m) + vn(ni-4,m) )          +  &
                               3.0*( vn(ni  ,m) + vn(ni-5,m) ) ) /256.0 -  &
                         aa*( 10.0*( vn(ni-2,m) - vn(ni-3,m) )          -  &
                               5.0*( vn(ni-1,m) - vn(ni-4,m) )          +  &
                                   ( vn(ni  ,m) - vn(ni-5,m) ) ) 
            
            vr(ni-2,m) =   ( 150.0*( vn(ni-2,m) + vn(ni-3,m) )          -  &
                              25.0*( vn(ni-1,m) + vn(ni-4,m) )          +  &
                               3.0*( vn(ni  ,m) + vn(ni-5,m) ) ) /256.0 +  &
                         aa*( 10.0*( vn(ni-2,m) - vn(ni-3,m) )          -  &
                               5.0*( vn(ni-1,m) - vn(ni-4,m) )          +  &
                                   ( vn(ni  ,m) - vn(ni-5,m) ) )              
        end do
        
    else
        
        ed  = ni + 3
        edc = ed
            
    end if
    
    do m=nst,ned

        do i=st-1,ed-1
            g(1,i,m) = -23.0/24.0*vn(i-3,m) + 31.0/8.0*vn(i-2,m) - 47.0/8.0*vn(i-1,m) + 71.0/24.0*vn(i  ,m)
            g(2,i,m) =   1.0/24.0*vn(i-2,m) -  1.0/8.0*vn(i-1,m) -  7.0/8.0*vn(i  ,m) + 23.0/24.0*vn(i+1,m)
            g(3,i,m) =   1.0/24.0*vn(i-1,m) -  9.0/8.0*vn(i  ,m) +  9.0/8.0*vn(i+1,m) -  1.0/24.0*vn(i+2,m)
            g(4,i,m) = -23.0/24.0*vn(i  ,m) +  7.0/8.0*vn(i+1,m) +  1.0/8.0*vn(i+2,m) -  1.0/24.0*vn(i+3,m)
            g(5,i,m) = -71.0/24.0*vn(i+1,m) + 47.0/8.0*vn(i+2,m) - 31.0/8.0*vn(i+3,m) + 23.0/24.0*vn(i+4,m)
            
            s(1,i,m) =   -3.0/2.0*vn(i-3,m) + 11.0/2.0*vn(i-2,m) - 13.0/2.0*vn(i-1,m) +   5.0/2.0*vn(i  ,m)
            s(2,i,m) =   -1.0/2.0*vn(i-2,m) +  5.0/2.0*vn(i-1,m) -  7.0/2.0*vn(i  ,m) +   3.0/2.0*vn(i+1,m)
            s(3,i,m) =    1.0/2.0*vn(i-1,m) -  1.0/2.0*vn(i  ,m) -  1.0/2.0*vn(i+1,m) +   1.0/2.0*vn(i+2,m)
            s(4,i,m) =    3.0/2.0*vn(i  ,m) -  7.0/2.0*vn(i+1,m) +  5.0/2.0*vn(i+2,m) -   1.0/2.0*vn(i+3,m)
            s(5,i,m) =    5.0/2.0*vn(i+1,m) - 13.0/2.0*vn(i+2,m) + 11.0/2.0*vn(i+3,m) -   3.0/2.0*vn(i+4,m)
            
            t(1,i,m) =    -vn(i-3,m) + 3.0*vn(i-2,m) - 3.0*vn(i-1,m) + vn(i  ,m)
            t(2,i,m) =    -vn(i-2,m) + 3.0*vn(i-1,m) - 3.0*vn(i  ,m) + vn(i+1,m)
            t(3,i,m) =    -vn(i-1,m) + 3.0*vn(i  ,m) - 3.0*vn(i+1,m) + vn(i+2,m)
            t(4,i,m) =    -vn(i  ,m) + 3.0*vn(i+1,m) - 3.0*vn(i+2,m) + vn(i+3,m)
            t(5,i,m) =    -vn(i+1,m) + 3.0*vn(i+2,m) - 3.0*vn(i+3,m) + vn(i+4,m)
        end do       
    end do
            
    do m=nst,ned
        do i=st-1,ed-1
            do n=1,5
                is = g(n,i,m)*g(n,i,m) + s(n,i,m)*s(n,i,m) + t(n,i,m)*t(n,i,m)              
                eis2(n) = is**2
            end do
            
            tau7 = abs(eis2(4)-eis2(1))
            
            do n=1,4
                bl(n) = cl(n)*(1.0+tau7/(eis2(n)+eps))
            end do
            
            bs = bl(1) + bl(2) + bl(3) + bl(4)
            do n=1,4
                wl(n,i,m) = bl(n)/bs
            end do
          
            tau7 = abs(eis2(5)-eis2(2))           
            do n=2,5
                br(n) = cr(n)*(1.0+tau7/(eis2(n)+eps))
            end do

            bs = br(2) + br(3) + br(4) + br(5)
            do n=2,5
                wr(n,i,m) = br(n)/bs
            end do
        end do
    end do
    
    do m=nst,ned
        do i=st-1,ed-1
            vl(i+1,m) = vn(i  ,m) + ( wl(1,i,m)*( 24.0*g(1,i,m)-6.0*s(1,i,m)+t(1,i,m)) + &
                                      wl(2,i,m)*( 24.0*g(2,i,m)-6.0*s(2,i,m)+t(2,i,m)) + &
                                      wl(3,i,m)*( 24.0*g(3,i,m)-6.0*s(3,i,m)+t(3,i,m)) + &
                                      wl(4,i,m)*( 24.0*g(4,i,m)-6.0*s(4,i,m)+t(4,i,m)) )/48.0
        end do

        do i=st,ed
            vr(i+1,m) = vn(i+1,m) + ( wr(2,i,m)*(-24.0*g(2,i,m)-6.0*s(2,i,m)-t(2,i,m)) + &
                                      wr(3,i,m)*(-24.0*g(3,i,m)-6.0*s(3,i,m)-t(3,i,m)) + &
                                      wr(4,i,m)*(-24.0*g(4,i,m)-6.0*s(4,i,m)-t(4,i,m)) + &
                                      wr(5,i,m)*(-24.0*g(5,i,m)-6.0*s(5,i,m)-t(5,i,m)) )/48.0
        end do
    end do    

    do i=stc,edc
        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vl(i,m) < chklim(m) ) then
                    vl(i,m) = vn(i  ,m)
                    ndec = ndec + 1
                end if

                if (vr(i,m) < chklim(m) ) then
                    vr(i,m) = vn(i+1,m)
                    ndec = ndec + 1
                end if
           end if
        end do
    end do

end subroutine scsn4picc

subroutine scsh5ci(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_bound_scheme
    implicit none
    integer(kind_int), intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  , intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(in)  :: sn(1-ngn:ni+ngn,1:3)
    integer(kind_int), intent(in)  :: nfs,nfe
    integer(kind_int), intent(in)  :: nge
    real(kind_real)  , intent(out) :: vl(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(out) :: vr(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , external    :: sub_limit
    integer(kind_int), intent(in)  :: chkid(nst:ned)
    real(kind_real),   intent(in)  :: chklim(nst:ned)
    integer(kind_int), intent(out) :: ndec  
    real(kind_real)   :: vln(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: vrn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: ve_lil(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: ve_lir(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: d2_li(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: d4_li(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: d5_li(1-ngn:ni+ngn,nst:ned)  
    real(kind_real)   :: nt,nx,ny,nz,fp(nst:ned)
    real(kind_real)   :: vp(nst:ned),vc(nst:ned)
    real(kind_real)   :: q(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: q1(1-ngn:ni+ngn,nst:ned)    
    real(kind_real)   :: s(1:5,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: g(1:5,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: t(1:5,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: wl(1:5,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: wr(1:5,1-ngn:ni+ngn,nst:ned)    
    real(kind_real)   :: aa,aa_scs
    real(kind_real)   :: sigma(1-ngn:ni+ngn),imax,imin,indicator
    real(kind_real), parameter :: cl(5) = (/1.D0/64.D0,21.D0/64.D0,35.D0/64.D0,7.D0/64.D0,0.D0/64.D0/)
    real(kind_real), parameter :: cr(5) = (/0.D0/64.D0,7.D0/64.D0,35.D0/64.D0,21.D0/64.D0,1.D0/64.D0/)  
    real(kind_real), parameter :: eps        = 1.0e-30
    real(kind_real), parameter :: indicatorc = 0.1
    real(kind_real)            :: is,eis2(5),bs,tau7,bl(5),br(5)
    integer(kind_int)          :: i,m,n,j,st,ed,stc,edc
    integer(kind_int)          :: id,YN,ind(1:ni+2*ngn),ind_st,ind_ed,pre_ed
    integer(kind_int)          :: nsec,sec_st(1:ni+2*ngn),sec_ed(1:ni+2*ngn)
    integer(kind_int)          :: nsecn,sec_stn(1:ni+2*ngn),sec_edn(1:ni+2*ngn)
    
    if (ni < 8) then
        call muscl2pv(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)
        return
    end if
    
    ndec = 0
    
    aa     = 0.01
    aa_scs = 0.001
    
    wl(:,:,:) = 0.25
    wr(:,:,:) = 0.25
    vln(:,:)    = 0.0
    ve_lil(:,:) = 0.0
    vrn(:,:)    = 0.0
    ve_lir(:,:) = 0.0    
           
    if (nfs == nbc_bound_scheme) then
        
        st  = 4
        stc = 0

        do m=nst,ned
            vl(0,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                  ! third-order
            vl(1,m) = (  3.0*vn(1,m) +  6.0*vn(2,m) -      vn(3,m) ) / 8.0   ! third-order 
            vl(2,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 + &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )
            vr(0,m) = vl(0,m)
            vr(1,m) = vl(1,m)
            vr(2,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 - &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )
            vl(3,m) =   ( 150.0*( vn(4,m) + vn(3,m) )          -  &
                           25.0*( vn(5,m) + vn(2,m) )          +  &
                            3.0*( vn(6,m) + vn(1,m) ) ) /256.0 -  &
                      aa*( 10.0*( vn(4,m) - vn(3,m) )          -  &
                            5.0*( vn(5,m) - vn(2,m) )          +  &
                                ( vn(6,m) - vn(1,m) ) ) 
            
            vr(3,m) =   ( 150.0*( vn(4,m) + vn(3,m) )          -  &
                           25.0*( vn(5,m) + vn(2,m) )          +  &
                            3.0*( vn(6,m) + vn(1,m) ) ) /256.0 +  &
                      aa*( 10.0*( vn(4,m) - vn(3,m) )          -  &
                            5.0*( vn(5,m) - vn(2,m) )          +  &
                                ( vn(6,m) - vn(1,m) ) )          
        end do
        
    else
        
        st  = -2
        stc = st
     
    end if
     
    if (nfe == nbc_bound_scheme) then
        
        ed  = ni - 4
        edc = ni

        do m=nst,ned
            vr(ni  ,m) = ( 15.0*vn(ni,m) - 10.0*vn(ni-1,m) +  3.0*vn(ni-2,m) ) / 8.0                ! third-order
            vl(ni-1,m) = (  3.0*vn(ni,m) +  6.0*vn(ni-1,m) -      vn(ni-2,m) ) / 8.0   ! third-order
            vl(ni-2,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 + &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) )
            vl(ni  ,m) = vr(ni  ,m)
            vr(ni-1,m) = vl(ni-1,m)
            vr(ni-2,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 - &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) ) 
            vl(ni-3,m) =   ( 150.0*( vn(ni-2,m) + vn(ni-3,m) )          -  &
                              25.0*( vn(ni-1,m) + vn(ni-4,m) )          +  &
                               3.0*( vn(ni  ,m) + vn(ni-5,m) ) ) /256.0 -  &
                         aa*( 10.0*( vn(ni-2,m) - vn(ni-3,m) )          -  &
                               5.0*( vn(ni-1,m) - vn(ni-4,m) )          +  &
                                   ( vn(ni  ,m) - vn(ni-5,m) ) ) 
            
            vr(ni-3,m) =   ( 150.0*( vn(ni-2,m) + vn(ni-3,m) )          -  &
                              25.0*( vn(ni-1,m) + vn(ni-4,m) )          +  &
                               3.0*( vn(ni  ,m) + vn(ni-5,m) ) ) /256.0 +  &
                         aa*( 10.0*( vn(ni-2,m) - vn(ni-3,m) )          -  &
                               5.0*( vn(ni-1,m) - vn(ni-4,m) )          +  &
                                   ( vn(ni  ,m) - vn(ni-5,m) ) )              
        end do
        
    else
        
        ed  = ni + 2
        edc = ed
            
    end if
    
    do m=nst,ned
        do i=st-2,ed+2
            g(1,i,m) = vn(i  ,m) - vn(i-1,m)            
            g(2,i,m) = vn(i+1,m) - vn(i  ,m)
            g(3,i,m) = vn(i+2,m) - vn(i+1,m)            
        end do
    end do   
    
    do m=nst,ned
        do i=st-2,ed+2
            do n=1,3
                is = g(n,i,m)*g(n,i,m)
                eis2(n) = is**2
            end do

            tau7 = abs(eis2(2)-eis2(1))
            bl(1) = 1.0+tau7/(eis2(1)+eps)
            bl(2) = 1.0+tau7/(eis2(2)+eps)
            bs = bl(1) + bl(2)
            wl(2,i,m) = bl(2)/bs
            
            tau7 = abs(eis2(3)-eis2(2))
            br(2) = 1.0+tau7/(eis2(2)+eps)
            br(3) = 1.0+tau7/(eis2(3)+eps)
            bs = br(2) + br(3)
            wr(2,i,m) = br(2)/bs
        end do
    end do
    
    do i=st,ed
        imin=1.0
        imax=0.0
        do m=nst,ned
            imin         = min(wl(2,i-3,m),wl(2,i-2,m),wl(2,i-1,m),wl(2,i,m),wl(2,i+1,m),wl(2,i+2,m),wl(2,i+3,m), &
                               wr(2,i-3,m),wr(2,i-2,m),wr(2,i-1,m),wr(2,i,m),wr(2,i+1,m),wr(2,i+2,m),wr(2,i+3,m),imin)                                
            imax         = max(wl(2,i-3,m),wl(2,i-2,m),wl(2,i-1,m),wl(2,i,m),wl(2,i+1,m),wl(2,i+2,m),wl(2,i+3,m), &
                               wr(2,i-3,m),wr(2,i-2,m),wr(2,i-1,m),wr(2,i,m),wr(2,i+1,m),wr(2,i+2,m),wr(2,i+3,m),imax)
        end do
        indicator    = imin/imax
        if(indicator>indicatorc) then
            sigma(i)=1.0
        else
            sigma(i)=0.0
        endif        
    end do
    
    id=0
    do i=st,ed
        YN=0
        if(sigma(i)<0.5) then
            YN=YN+1
        end if
        if(YN>0) then
            id=id+1
            ind(id)=i           
        end if
    end do
       
    do i=1-ngn,ni+ngn
        sigma(i)=0.0
    end do
    if(ed-st<9) then
        ind_st=st
        ind_ed=ed
        nsec=0
        
        nsecn=1
        sec_stn(1)=st
        sec_edn(1)=ed
    else
        if(id==0) then
            nsec=1
            sec_st(1)=st
            sec_ed(1)=ed
            do n=sec_st(1),sec_ed(1)
                sigma(n)=1.0
            end do
            ind_st=st
            ind_ed=ed
            nsecn=0
        else
            ind_st=st
            ind_ed=ed
            nsec=0
            nsecn=0
            if(ind(1)-ind_st>9) then
                nsec=1
                sec_st(1)=ind_st
                sec_ed(1)=ind(1)-1
                pre_ed   =sec_ed(1)
                do n=sec_st(1),sec_ed(1)
                    sigma(n)=1.0
                end do
            else
                nsecn=1
                sec_stn(1)=ind_st
                sec_edn(1)=ind(1)
                pre_ed    =sec_edn(1)
            end if
            do i=2,id
                if(ind(i)-ind(i-1)>9) then
                    nsec=nsec+1
                    sec_st(nsec)=ind(i-1)+1
                    sec_ed(nsec)=ind(i)-1
                    do n=sec_st(nsec),sec_ed(nsec)
                        sigma(n)=1.0
                    end do
                    nsecn=nsecn+1
                    sec_stn(nsecn)=pre_ed+1
                    sec_edn(nsecn)=ind(i-1)
                    pre_ed        =sec_ed(nsec)
                else
                    if(ind(i)-ind(i-1)>1 .or. i==id) then
                        nsecn=nsecn+1
                        sec_stn(nsecn)=pre_ed+1
                        sec_edn(nsecn)=ind(i)
                        pre_ed        =sec_edn(nsecn)
                    end if
                end if
            end do
            if(ind_ed-ind(id)>9) then
                nsec=nsec+1
                sec_st(nsec)=pre_ed+1
                sec_ed(nsec)=ind_ed
                do n=sec_st(nsec),sec_ed(nsec)
                    sigma(n)=1.0
                end do
            else
                nsecn=nsecn+1
                sec_stn(nsecn)=pre_ed+1
                sec_edn(nsecn)=ind_ed
            end if
        end if
    end if
    
    do m=nst,ned
        do j=1,nsecn
            do i=sec_stn(j),sec_edn(j)
                g(1,i,m) = -23.0/24.0*vn(i-3,m) + 31.0/8.0*vn(i-2,m) - 47.0/8.0*vn(i-1,m) + 71.0/24.0*vn(i  ,m)
                g(2,i,m) =   1.0/24.0*vn(i-2,m) -  1.0/8.0*vn(i-1,m) -  7.0/8.0*vn(i  ,m) + 23.0/24.0*vn(i+1,m)
                g(3,i,m) =   1.0/24.0*vn(i-1,m) -  9.0/8.0*vn(i  ,m) +  9.0/8.0*vn(i+1,m) -  1.0/24.0*vn(i+2,m)
                g(4,i,m) = -23.0/24.0*vn(i  ,m) +  7.0/8.0*vn(i+1,m) +  1.0/8.0*vn(i+2,m) -  1.0/24.0*vn(i+3,m)
                g(5,i,m) = -71.0/24.0*vn(i+1,m) + 47.0/8.0*vn(i+2,m) - 31.0/8.0*vn(i+3,m) + 23.0/24.0*vn(i+4,m)
                
                s(1,i,m) =   -3.0/2.0*vn(i-3,m) + 11.0/2.0*vn(i-2,m) - 13.0/2.0*vn(i-1,m) +   5.0/2.0*vn(i  ,m)
                s(2,i,m) =   -1.0/2.0*vn(i-2,m) +  5.0/2.0*vn(i-1,m) -  7.0/2.0*vn(i  ,m) +   3.0/2.0*vn(i+1,m)
                s(3,i,m) =    1.0/2.0*vn(i-1,m) -  1.0/2.0*vn(i  ,m) -  1.0/2.0*vn(i+1,m) +   1.0/2.0*vn(i+2,m)
                s(4,i,m) =    3.0/2.0*vn(i  ,m) -  7.0/2.0*vn(i+1,m) +  5.0/2.0*vn(i+2,m) -   1.0/2.0*vn(i+3,m)
                s(5,i,m) =    5.0/2.0*vn(i+1,m) - 13.0/2.0*vn(i+2,m) + 11.0/2.0*vn(i+3,m) -   3.0/2.0*vn(i+4,m)
                
                t(1,i,m) =    -vn(i-3,m) + 3.0*vn(i-2,m) - 3.0*vn(i-1,m) + vn(i  ,m)
                t(2,i,m) =    -vn(i-2,m) + 3.0*vn(i-1,m) - 3.0*vn(i  ,m) + vn(i+1,m)
                t(3,i,m) =    -vn(i-1,m) + 3.0*vn(i  ,m) - 3.0*vn(i+1,m) + vn(i+2,m)
                t(4,i,m) =    -vn(i  ,m) + 3.0*vn(i+1,m) - 3.0*vn(i+2,m) + vn(i+3,m)
                t(5,i,m) =    -vn(i+1,m) + 3.0*vn(i+2,m) - 3.0*vn(i+3,m) + vn(i+4,m)
            end do
        end do      
    end do
    
    do j=1,nsecn
        do i=sec_stn(j),sec_edn(j)
            nx = (sn(i,1) + sn(i+1,1))/2.0
            ny = (sn(i,2) + sn(i+1,2))/2.0
            nz = (sn(i,3) + sn(i+1,3))/2.0

            do m=nst,ned
                vp(m) = (vn(i,m) + vn(i+1,m))/2.0
                fp(m) = vn(i,m)
            end do
            call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
            do m=nst,ned
                q(i,m) = vc(m)
            end do
            
            do m=nst,ned
                fp(m) = vn(i+1,m)
            end do
            call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
            do m=nst,ned
                q1(i+1,m) = vc(m)
            end do        
             
            do n=1,5
                do m=nst,ned
                    fp(m) = g(n,i,m)
                end do
                call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
                do m=nst,ned
                    g(n,i,m) = vc(m)
                end do

                do m=nst,ned
                    fp(m) = s(n,i,m)
                end do
                call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
                do m=nst,ned
                    s(n,i,m) = vc(m)
                end do
                
                do m=nst,ned
                    fp(m) = t(n,i,m)
                end do
                call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
                do m=nst,ned
                    t(n,i,m) = vc(m)
                end do            
            end do
        end do
    end do   
            
            
    do m=nst,ned
        do j=1,nsecn
            do i=sec_stn(j),sec_edn(j)
                do n=1,5
                    is = g(n,i,m)*g(n,i,m) + s(n,i,m)*s(n,i,m) + t(n,i,m)*t(n,i,m)              
                    eis2(n) = is**2
                end do
                
                tau7 = abs(eis2(4)-eis2(1))
                
                do n=1,4
                    bl(n) = cl(n)*(1.0+tau7/(eis2(n)+eps))
                end do
                
                bs = bl(1) + bl(2) + bl(3) + bl(4)
                do n=1,4
                    wl(n,i,m) = bl(n)/bs
                end do
          
                tau7 = abs(eis2(5)-eis2(2))           
                do n=2,5
                    br(n) = cr(n)*(1.0+tau7/(eis2(n)+eps))
                end do

                bs = br(2) + br(3) + br(4) + br(5)
                do n=2,5
                    wr(n,i,m) = br(n)/bs
                end do
            end do
        end do
    end do
    
!!    open(2,file='sigma.dat',form='formatted')
!!    write(2,*) 'variables="i" "sigmal" "sigmar"'
!!    do i=st,ed
!!
!!        write(2,200) i,sigmal(i,1),sigmar(i,1)
!!
!!    end do
!!
!!200 format(I4,1x,e12.5,1x,e12.5)
!!    close(2)
!!
!!    stop    
    
    do m=nst,ned
        do j=1,nsecn
            do i=sec_stn(j),sec_edn(j)
                vln(i,m) =  q(i  ,m) + ( wl(1,i,m)*( 24.0*g(1,i,m)-6.0*s(1,i,m)+t(1,i,m)) + &
                                         wl(2,i,m)*( 24.0*g(2,i,m)-6.0*s(2,i,m)+t(2,i,m)) + &
                                         wl(3,i,m)*( 24.0*g(3,i,m)-6.0*s(3,i,m)+t(3,i,m)) + &
                                         wl(4,i,m)*( 24.0*g(4,i,m)-6.0*s(4,i,m)+t(4,i,m)) )/48.0

                vrn(i,m) = q1(i+1,m) + ( wr(2,i,m)*(-24.0*g(2,i,m)-6.0*s(2,i,m)-t(2,i,m)) + &
                                         wr(3,i,m)*(-24.0*g(3,i,m)-6.0*s(3,i,m)-t(3,i,m)) + &
                                         wr(4,i,m)*(-24.0*g(4,i,m)-6.0*s(4,i,m)-t(4,i,m)) + &
                                         wr(5,i,m)*(-24.0*g(5,i,m)-6.0*s(5,i,m)-t(5,i,m)) )/48.0
            end do
        end do
    end do
    
    do j=1,nsecn
        do i=sec_stn(j),sec_edn(j)
            nx = (sn(i,1) + sn(i+1,1))/2.0
            ny = (sn(i,2) + sn(i+1,2))/2.0
            nz = (sn(i,3) + sn(i+1,3))/2.0

            do m=nst,ned
                vp(m) = (vn(i,m) + vn(i+1,m))/2.0
                vc(m) = vln(i,m)
            end do

            call RxQ_p(nst,ned,vp,nt,nx,ny,nz,vc,fp) 
            do m=nst,ned
                vln(i,m) = fp(m)
            end do
            
            do m=nst,ned
                vc(m) = vrn(i ,m)
            end do

            call RxQ_p(nst,ned,vp,nt,nx,ny,nz,vc,fp) 
            do m=nst,ned
                vrn(i,m) = fp(m)
            end do
        end do
    end do     
      
    do j=1,nsec
        call d5_via_T2N(ni,nst,ned,ngn,vn,nfs,nfe,nge,d5_li,sec_st(j),sec_ed(j))
        
        call vec2_via_T4N(ni,nst,ned,ngn,vn,nfs,nfe,nge,d2_li,sec_st(j),sec_ed(j))
        
        call vec4_via_T2N(ni,nst,ned,ngn,vn,nfs,nfe,nge,d4_li,sec_st(j),sec_ed(j))
    end do
    
    do m=nst,ned
        do j=1,nsec
            do i=sec_st(j),sec_ed(j)
                
                ve_lil(i,m) = 0.5*(vn(i+1,m) + vn(i  ,m)) - 1.0/8.0*d2_li(i,m) - 1.0/384.0*d4_li(i,m) -  aa_scs*d5_li(i,m)
                            
                ve_lir(i,m) = 0.5*(vn(i+1,m) + vn(i  ,m)) - 1.0/8.0*d2_li(i,m) - 1.0/384.0*d4_li(i,m) +  aa_scs*d5_li(i,m)

            end do
        end do
    end do
        
    do m=nst,ned
        do i=st,ed
            
               vl(i,m) = sigma(i)*ve_lil(i,m) + vln(i,m)*(1.0-sigma(i))
                         
               vr(i,m) = sigma(i)*ve_lir(i,m) + vrn(i,m)*(1.0-sigma(i))            

        end do
    end do  

    do i=stc,edc
        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vl(i,m) < chklim(m) ) then
                    vl(i,m) = vn(i  ,m)
                    ndec = ndec + 1
                end if

                if (vr(i,m) < chklim(m) ) then
                    vr(i,m) = vn(i+1,m)
                    ndec = ndec + 1
                end if
           end if
        end do
    end do

end subroutine scsh5ci

subroutine scsh5cicc(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_bound_scheme
    implicit none
    integer(kind_int), intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  , intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(in)  :: sn(1-ngn:ni+ngn,1:3)
    integer(kind_int), intent(in)  :: nfs,nfe
    integer(kind_int), intent(in)  :: nge
    real(kind_real)  , intent(out) :: vl(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(out) :: vr(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , external    :: sub_limit
    integer(kind_int), intent(in)  :: chkid(nst:ned)
    real(kind_real),   intent(in)  :: chklim(nst:ned)
    integer(kind_int), intent(out) :: ndec  
    real(kind_real)   :: vln(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: vrn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: ve_lil(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: ve_lir(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: d2_li(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: d4_li(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: d5_li(1-ngn:ni+ngn,nst:ned)  
    real(kind_real)   :: nt,nx,ny,nz,fp(nst:ned)
    real(kind_real)   :: vp(nst:ned),vc(nst:ned)
    real(kind_real)   :: q(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: q1(1-ngn:ni+ngn,nst:ned)    
    real(kind_real)   :: s(1:5,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: g(1:5,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: t(1:5,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: wl(1:5,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: wr(1:5,1-ngn:ni+ngn,nst:ned)    
    real(kind_real)   :: aa,aa_scs
    real(kind_real)   :: sigma(1-ngn:ni+ngn),imax,imin,indicator
    real(kind_real), parameter :: cl(5) = (/1.D0/64.D0,21.D0/64.D0,35.D0/64.D0,7.D0/64.D0,0.D0/64.D0/)
    real(kind_real), parameter :: cr(5) = (/0.D0/64.D0,7.D0/64.D0,35.D0/64.D0,21.D0/64.D0,1.D0/64.D0/)  
    real(kind_real), parameter :: eps        = 1.0e-30
    real(kind_real), parameter :: indicatorc = 0.1
    real(kind_real)            :: is,eis2(5),bs,tau7,bl(5),br(5)
    integer(kind_int)          :: i,m,n,j,st,ed,stc,edc
    integer(kind_int)          :: id,YN,ind(1:ni+2*ngn),ind_st,ind_ed,pre_ed
    integer(kind_int)          :: nsec,sec_st(1:ni+2*ngn),sec_ed(1:ni+2*ngn)
    integer(kind_int)          :: nsecn,sec_stn(1:ni+2*ngn),sec_edn(1:ni+2*ngn)
    
    if (ni < 8) then
        call muscl2pvcc(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)
        return
    end if
    
    ndec = 0
    
    aa     = 0.01
    aa_scs = 0.001
    
    wl(:,:,:) = 0.25
    wr(:,:,:) = 0.25
    vln(:,:)    = 0.0
    ve_lil(:,:) = 0.0
    vrn(:,:)    = 0.0
    ve_lir(:,:) = 0.0    
           
    if (nfs == nbc_bound_scheme) then
        
        st  = 5
        stc = 1

        do m=nst,ned
            vl(1,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                  ! third-order
            vl(2,m) = (  3.0*vn(1,m) +  6.0*vn(2,m) -      vn(3,m) ) / 8.0   ! third-order 
            vl(3,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 + &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )
            vr(1,m) = vl(1,m)
            vr(2,m) = vl(2,m)
            vr(3,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 - &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )
            vl(4,m) =   ( 150.0*( vn(4,m) + vn(3,m) )          -  &
                           25.0*( vn(5,m) + vn(2,m) )          +  &
                            3.0*( vn(6,m) + vn(1,m) ) ) /256.0 -  &
                      aa*( 10.0*( vn(4,m) - vn(3,m) )          -  &
                            5.0*( vn(5,m) - vn(2,m) )          +  &
                                ( vn(6,m) - vn(1,m) ) ) 
            
            vr(4,m) =   ( 150.0*( vn(4,m) + vn(3,m) )          -  &
                           25.0*( vn(5,m) + vn(2,m) )          +  &
                            3.0*( vn(6,m) + vn(1,m) ) ) /256.0 +  &
                      aa*( 10.0*( vn(4,m) - vn(3,m) )          -  &
                            5.0*( vn(5,m) - vn(2,m) )          +  &
                                ( vn(6,m) - vn(1,m) ) )          
        end do
        
    else
        
        st  = -1
        stc = st
     
    end if
     
    if (nfe == nbc_bound_scheme) then
        
        ed  = ni - 3
        edc = ni + 1

        do m=nst,ned
            vr(ni+1,m) = ( 15.0*vn(ni,m) - 10.0*vn(ni-1,m) +  3.0*vn(ni-2,m) ) / 8.0                ! third-order
            vl(ni  ,m) = (  3.0*vn(ni,m) +  6.0*vn(ni-1,m) -      vn(ni-2,m) ) / 8.0   ! third-order
            vl(ni-1,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 + &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) )
            vl(ni+1,m) = vr(ni+1,m)
            vr(ni  ,m) = vl(ni  ,m)
            vr(ni-1,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 - &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) ) 
            vl(ni-2,m) =   ( 150.0*( vn(ni-2,m) + vn(ni-3,m) )          -  &
                              25.0*( vn(ni-1,m) + vn(ni-4,m) )          +  &
                               3.0*( vn(ni  ,m) + vn(ni-5,m) ) ) /256.0 -  &
                         aa*( 10.0*( vn(ni-2,m) - vn(ni-3,m) )          -  &
                               5.0*( vn(ni-1,m) - vn(ni-4,m) )          +  &
                                   ( vn(ni  ,m) - vn(ni-5,m) ) ) 
            
            vr(ni-2,m) =   ( 150.0*( vn(ni-2,m) + vn(ni-3,m) )          -  &
                              25.0*( vn(ni-1,m) + vn(ni-4,m) )          +  &
                               3.0*( vn(ni  ,m) + vn(ni-5,m) ) ) /256.0 +  &
                         aa*( 10.0*( vn(ni-2,m) - vn(ni-3,m) )          -  &
                               5.0*( vn(ni-1,m) - vn(ni-4,m) )          +  &
                                   ( vn(ni  ,m) - vn(ni-5,m) ) )              
        end do
        
    else
        
        ed  = ni + 3
        edc = ed
            
    end if
    
    do m=nst,ned
        do i=st-3,ed+1
            g(1,i,m) = vn(i  ,m) - vn(i-1,m)            
            g(2,i,m) = vn(i+1,m) - vn(i  ,m)
            g(3,i,m) = vn(i+2,m) - vn(i+1,m)            
        end do
    end do   
    
    do m=nst,ned
        do i=st-3,ed+1
            do n=1,3
                is = g(n,i,m)*g(n,i,m)
                eis2(n) = is**2
            end do

            tau7 = abs(eis2(2)-eis2(1))
            bl(1) = 1.0+tau7/(eis2(1)+eps)
            bl(2) = 1.0+tau7/(eis2(2)+eps)
            bs = bl(1) + bl(2)
            wl(2,i,m) = bl(2)/bs
            
            tau7 = abs(eis2(3)-eis2(2))
            br(2) = 1.0+tau7/(eis2(2)+eps)
            br(3) = 1.0+tau7/(eis2(3)+eps)
            bs = br(2) + br(3)
            wr(2,i,m) = br(2)/bs
        end do
    end do
    
    do i=st-1,ed-1
        imin=1.0
        imax=0.0
        do m=nst,ned
            imin         = min(wl(2,i-3,m),wl(2,i-2,m),wl(2,i-1,m),wl(2,i,m),wl(2,i+1,m),wl(2,i+2,m),wl(2,i+3,m), &
                               wr(2,i-3,m),wr(2,i-2,m),wr(2,i-1,m),wr(2,i,m),wr(2,i+1,m),wr(2,i+2,m),wr(2,i+3,m),imin)                                
            imax         = max(wl(2,i-3,m),wl(2,i-2,m),wl(2,i-1,m),wl(2,i,m),wl(2,i+1,m),wl(2,i+2,m),wl(2,i+3,m), &
                               wr(2,i-3,m),wr(2,i-2,m),wr(2,i-1,m),wr(2,i,m),wr(2,i+1,m),wr(2,i+2,m),wr(2,i+3,m),imax)
        end do
        indicator    = imin/imax
        if(indicator>indicatorc) then
            sigma(i)=1.0
        else
            sigma(i)=0.0
        endif        
    end do
    
    id=0
    do i=st-1,ed-1
        YN=0
        if(sigma(i)<0.5) then
            YN=YN+1
        end if
        if(YN>0) then
            id=id+1
            ind(id)=i           
        end if
    end do
       
    do i=1-ngn,ni+ngn
        sigma(i)=0.0
    end do
    if(ed-st<9) then
        ind_st=st-1
        ind_ed=ed-1
        nsec=0
        
        nsecn=1
        sec_stn(1)=st-1
        sec_edn(1)=ed-1
    else
        if(id==0) then
            nsec=1
            sec_st(1)=st-1
            sec_ed(1)=ed-1
            do n=sec_st(1),sec_ed(1)
                sigma(n)=1.0
            end do
            ind_st=st-1
            ind_ed=ed-1
            nsecn=0
        else
            ind_st=st-1
            ind_ed=ed-1
            nsec=0
            nsecn=0
            if(ind(1)-ind_st>9) then
                nsec=1
                sec_st(1)=ind_st
                sec_ed(1)=ind(1)-1
                pre_ed   =sec_ed(1)
                do n=sec_st(1),sec_ed(1)
                    sigma(n)=1.0
                end do
            else
                nsecn=1
                sec_stn(1)=ind_st
                sec_edn(1)=ind(1)
                pre_ed    =sec_edn(1)
            end if
            do i=2,id
                if(ind(i)-ind(i-1)>9) then
                    nsec=nsec+1
                    sec_st(nsec)=ind(i-1)+1
                    sec_ed(nsec)=ind(i)-1
                    do n=sec_st(nsec),sec_ed(nsec)
                        sigma(n)=1.0
                    end do
                    nsecn=nsecn+1
                    sec_stn(nsecn)=pre_ed+1
                    sec_edn(nsecn)=ind(i-1)
                    pre_ed        =sec_ed(nsec)
                else
                    if(ind(i)-ind(i-1)>1 .or. i==id) then
                        nsecn=nsecn+1
                        sec_stn(nsecn)=pre_ed+1
                        sec_edn(nsecn)=ind(i)
                        pre_ed        =sec_edn(nsecn)
                    end if
                end if
            end do
            if(ind_ed-ind(id)>9) then
                nsec=nsec+1
                sec_st(nsec)=pre_ed+1
                sec_ed(nsec)=ind_ed
                do n=sec_st(nsec),sec_ed(nsec)
                    sigma(n)=1.0
                end do
            else
                nsecn=nsecn+1
                sec_stn(nsecn)=pre_ed+1
                sec_edn(nsecn)=ind_ed
            end if
        end if
    end if
    
    do m=nst,ned
        do j=1,nsecn
            do i=sec_stn(j),sec_edn(j)
                g(1,i,m) = -23.0/24.0*vn(i-3,m) + 31.0/8.0*vn(i-2,m) - 47.0/8.0*vn(i-1,m) + 71.0/24.0*vn(i  ,m)
                g(2,i,m) =   1.0/24.0*vn(i-2,m) -  1.0/8.0*vn(i-1,m) -  7.0/8.0*vn(i  ,m) + 23.0/24.0*vn(i+1,m)
                g(3,i,m) =   1.0/24.0*vn(i-1,m) -  9.0/8.0*vn(i  ,m) +  9.0/8.0*vn(i+1,m) -  1.0/24.0*vn(i+2,m)
                g(4,i,m) = -23.0/24.0*vn(i  ,m) +  7.0/8.0*vn(i+1,m) +  1.0/8.0*vn(i+2,m) -  1.0/24.0*vn(i+3,m)
                g(5,i,m) = -71.0/24.0*vn(i+1,m) + 47.0/8.0*vn(i+2,m) - 31.0/8.0*vn(i+3,m) + 23.0/24.0*vn(i+4,m)
                
                s(1,i,m) =   -3.0/2.0*vn(i-3,m) + 11.0/2.0*vn(i-2,m) - 13.0/2.0*vn(i-1,m) +   5.0/2.0*vn(i  ,m)
                s(2,i,m) =   -1.0/2.0*vn(i-2,m) +  5.0/2.0*vn(i-1,m) -  7.0/2.0*vn(i  ,m) +   3.0/2.0*vn(i+1,m)
                s(3,i,m) =    1.0/2.0*vn(i-1,m) -  1.0/2.0*vn(i  ,m) -  1.0/2.0*vn(i+1,m) +   1.0/2.0*vn(i+2,m)
                s(4,i,m) =    3.0/2.0*vn(i  ,m) -  7.0/2.0*vn(i+1,m) +  5.0/2.0*vn(i+2,m) -   1.0/2.0*vn(i+3,m)
                s(5,i,m) =    5.0/2.0*vn(i+1,m) - 13.0/2.0*vn(i+2,m) + 11.0/2.0*vn(i+3,m) -   3.0/2.0*vn(i+4,m)
                
                t(1,i,m) =    -vn(i-3,m) + 3.0*vn(i-2,m) - 3.0*vn(i-1,m) + vn(i  ,m)
                t(2,i,m) =    -vn(i-2,m) + 3.0*vn(i-1,m) - 3.0*vn(i  ,m) + vn(i+1,m)
                t(3,i,m) =    -vn(i-1,m) + 3.0*vn(i  ,m) - 3.0*vn(i+1,m) + vn(i+2,m)
                t(4,i,m) =    -vn(i  ,m) + 3.0*vn(i+1,m) - 3.0*vn(i+2,m) + vn(i+3,m)
                t(5,i,m) =    -vn(i+1,m) + 3.0*vn(i+2,m) - 3.0*vn(i+3,m) + vn(i+4,m)
            end do
        end do      
    end do
    
    do j=1,nsecn
        do i=sec_stn(j),sec_edn(j)
            nx = (sn(i,1) + sn(i+1,1))/2.0
            ny = (sn(i,2) + sn(i+1,2))/2.0
            nz = (sn(i,3) + sn(i+1,3))/2.0

            do m=nst,ned
                vp(m) = (vn(i,m) + vn(i+1,m))/2.0
                fp(m) = vn(i,m)
            end do
            call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
            do m=nst,ned
                q(i,m) = vc(m)
            end do
            
            do m=nst,ned
                fp(m) = vn(i+1,m)
            end do
            call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
            do m=nst,ned
                q1(i+1,m) = vc(m)
            end do        
             
            do n=1,5
                do m=nst,ned
                    fp(m) = g(n,i,m)
                end do
                call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
                do m=nst,ned
                    g(n,i,m) = vc(m)
                end do

                do m=nst,ned
                    fp(m) = s(n,i,m)
                end do
                call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
                do m=nst,ned
                    s(n,i,m) = vc(m)
                end do
                
                do m=nst,ned
                    fp(m) = t(n,i,m)
                end do
                call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
                do m=nst,ned
                    t(n,i,m) = vc(m)
                end do            
            end do
        end do
    end do   
            
            
    do m=nst,ned
        do j=1,nsecn
            do i=sec_stn(j),sec_edn(j)
                do n=1,5
                    is = g(n,i,m)*g(n,i,m) + s(n,i,m)*s(n,i,m) + t(n,i,m)*t(n,i,m)              
                    eis2(n) = is**2
                end do
                
                tau7 = abs(eis2(4)-eis2(1))
                
                do n=1,4
                    bl(n) = cl(n)*(1.0+tau7/(eis2(n)+eps))
                end do
                
                bs = bl(1) + bl(2) + bl(3) + bl(4)
                do n=1,4
                    wl(n,i,m) = bl(n)/bs
                end do
          
                tau7 = abs(eis2(5)-eis2(2))           
                do n=2,5
                    br(n) = cr(n)*(1.0+tau7/(eis2(n)+eps))
                end do

                bs = br(2) + br(3) + br(4) + br(5)
                do n=2,5
                    wr(n,i,m) = br(n)/bs
                end do
            end do
        end do
    end do
    
!!    open(2,file='sigma.dat',form='formatted')
!!    write(2,*) 'variables="i" "sigmal" "sigmar"'
!!    do i=st,ed
!!
!!        write(2,200) i,sigmal(i,1),sigmar(i,1)
!!
!!    end do
!!
!!200 format(I4,1x,e12.5,1x,e12.5)
!!    close(2)
!!
!!    stop    
    
    do m=nst,ned
        do j=1,nsecn
            do i=sec_stn(j),sec_edn(j)
                vln(i,m) =  q(i  ,m) + ( wl(1,i,m)*( 24.0*g(1,i,m)-6.0*s(1,i,m)+t(1,i,m)) + &
                                         wl(2,i,m)*( 24.0*g(2,i,m)-6.0*s(2,i,m)+t(2,i,m)) + &
                                         wl(3,i,m)*( 24.0*g(3,i,m)-6.0*s(3,i,m)+t(3,i,m)) + &
                                         wl(4,i,m)*( 24.0*g(4,i,m)-6.0*s(4,i,m)+t(4,i,m)) )/48.0

                vrn(i,m) = q1(i+1,m) + ( wr(2,i,m)*(-24.0*g(2,i,m)-6.0*s(2,i,m)-t(2,i,m)) + &
                                         wr(3,i,m)*(-24.0*g(3,i,m)-6.0*s(3,i,m)-t(3,i,m)) + &
                                         wr(4,i,m)*(-24.0*g(4,i,m)-6.0*s(4,i,m)-t(4,i,m)) + &
                                         wr(5,i,m)*(-24.0*g(5,i,m)-6.0*s(5,i,m)-t(5,i,m)) )/48.0
            end do
        end do
    end do
    
    do j=1,nsecn
        do i=sec_stn(j),sec_edn(j)
            nx = (sn(i,1) + sn(i+1,1))/2.0
            ny = (sn(i,2) + sn(i+1,2))/2.0
            nz = (sn(i,3) + sn(i+1,3))/2.0

            do m=nst,ned
                vp(m) = (vn(i,m) + vn(i+1,m))/2.0
                vc(m) = vln(i,m)
            end do

            call RxQ_p(nst,ned,vp,nt,nx,ny,nz,vc,fp) 
            do m=nst,ned
                vln(i,m) = fp(m)
            end do
            
            do m=nst,ned
                vc(m) = vrn(i ,m)
            end do

            call RxQ_p(nst,ned,vp,nt,nx,ny,nz,vc,fp) 
            do m=nst,ned
                vrn(i,m) = fp(m)
            end do
        end do
    end do     
      
    do j=1,nsec
        call d5_via_T2N(ni,nst,ned,ngn,vn,nfs,nfe,nge,d5_li,sec_st(j),sec_ed(j))
        
        call vec2_via_T4N(ni,nst,ned,ngn,vn,nfs,nfe,nge,d2_li,sec_st(j),sec_ed(j))
        
        call vec4_via_T2N(ni,nst,ned,ngn,vn,nfs,nfe,nge,d4_li,sec_st(j),sec_ed(j))
    end do
    
    do m=nst,ned
        do j=1,nsec
            do i=sec_st(j),sec_ed(j)
                
                ve_lil(i,m) = 0.5*(vn(i+1,m) + vn(i  ,m)) - 1.0/8.0*d2_li(i,m) - 1.0/384.0*d4_li(i,m) -  aa_scs*d5_li(i,m)
                            
                ve_lir(i,m) = 0.5*(vn(i+1,m) + vn(i  ,m)) - 1.0/8.0*d2_li(i,m) - 1.0/384.0*d4_li(i,m) +  aa_scs*d5_li(i,m)

            end do
        end do
    end do
        
    do m=nst,ned
        do i=st-1,ed-1
            
               vl(i+1,m) = sigma(i)*ve_lil(i,m) + vln(i,m)*(1.0-sigma(i))
                         
               vr(i+1,m) = sigma(i)*ve_lir(i,m) + vrn(i,m)*(1.0-sigma(i))            

        end do
    end do  

    do i=stc,edc
        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vl(i,m) < chklim(m) ) then
                    vl(i,m) = vn(i  ,m)
                    ndec = ndec + 1
                end if

                if (vr(i,m) < chklim(m) ) then
                    vr(i,m) = vn(i+1,m)
                    ndec = ndec + 1
                end if
           end if
        end do
    end do

end subroutine scsh5cicc

subroutine dcsh7ci(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_bound_scheme
    implicit none
    integer(kind_int), intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  , intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(in)  :: sn(1-ngn:ni+ngn,1:3)
    integer(kind_int), intent(in)  :: nfs,nfe
    integer(kind_int), intent(in)  :: nge
    real(kind_real)  , intent(out) :: vl(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(out) :: vr(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , external    :: sub_limit
    integer(kind_int), intent(in)  :: chkid(nst:ned)
    real(kind_real),   intent(in)  :: chklim(nst:ned)
    integer(kind_int), intent(out) :: ndec  
    real(kind_real)   :: vln(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: vrn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: ve_lil(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: ve_lir(1-ngn:ni+ngn,nst:ned) 
    real(kind_real)   :: nt,nx,ny,nz,fp(nst:ned)
    real(kind_real)   :: vp(nst:ned),vc(nst:ned)
    real(kind_real)   :: q(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: q1(1-ngn:ni+ngn,nst:ned)    
    real(kind_real)   :: s(1:5,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: g(1:5,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: t(1:5,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: wl(1:5,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: wr(1:5,1-ngn:ni+ngn,nst:ned)    
    real(kind_real)   :: aa,aa_scs
    real(kind_real)   :: sigma(1-ngn:ni+ngn),imax,imin,indicator
    real(kind_real), parameter :: cl(5) = (/1.D0/64.D0,21.D0/64.D0,35.D0/64.D0,7.D0/64.D0,0.D0/64.D0/)
    real(kind_real), parameter :: cr(5) = (/0.D0/64.D0,7.D0/64.D0,35.D0/64.D0,21.D0/64.D0,1.D0/64.D0/)  
    real(kind_real), parameter :: eps        = 1.0e-30
    real(kind_real), parameter :: indicatorc = 0.1
    real(kind_real)            :: is,eis2(5),bs,tau7,bl(5),br(5)
    integer(kind_int)          :: i,m,n,j,st,ed,stc,edc
    integer(kind_int)          :: id,YN,ind(1:ni+2*ngn),ind_st,ind_ed,pre_ed
    integer(kind_int)          :: nsec,sec_st(1:ni+2*ngn),sec_ed(1:ni+2*ngn)
    integer(kind_int)          :: nsecn,sec_stn(1:ni+2*ngn),sec_edn(1:ni+2*ngn)
    
    if (ni < 8) then
        call muscl2pv(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)
        return
    end if
    
    ndec = 0
    
    aa     = 0.01
    aa_scs = 0.001
    
    wl(:,:,:) = 0.25
    wr(:,:,:) = 0.25
    vln(:,:)    = 0.0
    ve_lil(:,:) = 0.0
    vrn(:,:)    = 0.0
    ve_lir(:,:) = 0.0
           
    if (nfs == nbc_bound_scheme) then
        
        st  = 4
        stc = 0

        do m=nst,ned
            vl(0,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                  ! third-order
            vl(1,m) = (  3.0*vn(1,m) +  6.0*vn(2,m) -      vn(3,m) ) / 8.0   ! third-order 
            vl(2,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 + &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )
            vr(0,m) = vl(0,m)
            vr(1,m) = vl(1,m)
            vr(2,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 - &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )
            vl(3,m) =   ( 150.0*( vn(4,m) + vn(3,m) )          -  &
                           25.0*( vn(5,m) + vn(2,m) )          +  &
                            3.0*( vn(6,m) + vn(1,m) ) ) /256.0 -  &
                      aa*( 10.0*( vn(4,m) - vn(3,m) )          -  &
                            5.0*( vn(5,m) - vn(2,m) )          +  &
                                ( vn(6,m) - vn(1,m) ) ) 
            
            vr(3,m) =   ( 150.0*( vn(4,m) + vn(3,m) )          -  &
                           25.0*( vn(5,m) + vn(2,m) )          +  &
                            3.0*( vn(6,m) + vn(1,m) ) ) /256.0 +  &
                      aa*( 10.0*( vn(4,m) - vn(3,m) )          -  &
                            5.0*( vn(5,m) - vn(2,m) )          +  &
                                ( vn(6,m) - vn(1,m) ) )          
        end do
        
    else
        
        st  = -2
        stc = st
     
    end if
     
    if (nfe == nbc_bound_scheme) then
        
        ed  = ni - 4
        edc = ni

        do m=nst,ned
            vr(ni  ,m) = ( 15.0*vn(ni,m) - 10.0*vn(ni-1,m) +  3.0*vn(ni-2,m) ) / 8.0                ! third-order
            vl(ni-1,m) = (  3.0*vn(ni,m) +  6.0*vn(ni-1,m) -      vn(ni-2,m) ) / 8.0   ! third-order
            vl(ni-2,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 + &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) )
            vl(ni  ,m) = vr(ni  ,m)
            vr(ni-1,m) = vl(ni-1,m)
            vr(ni-2,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 - &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) ) 
            vl(ni-3,m) =   ( 150.0*( vn(ni-2,m) + vn(ni-3,m) )          -  &
                              25.0*( vn(ni-1,m) + vn(ni-4,m) )          +  &
                               3.0*( vn(ni  ,m) + vn(ni-5,m) ) ) /256.0 -  &
                         aa*( 10.0*( vn(ni-2,m) - vn(ni-3,m) )          -  &
                               5.0*( vn(ni-1,m) - vn(ni-4,m) )          +  &
                                   ( vn(ni  ,m) - vn(ni-5,m) ) ) 
            
            vr(ni-3,m) =   ( 150.0*( vn(ni-2,m) + vn(ni-3,m) )          -  &
                              25.0*( vn(ni-1,m) + vn(ni-4,m) )          +  &
                               3.0*( vn(ni  ,m) + vn(ni-5,m) ) ) /256.0 +  &
                         aa*( 10.0*( vn(ni-2,m) - vn(ni-3,m) )          -  &
                               5.0*( vn(ni-1,m) - vn(ni-4,m) )          +  &
                                   ( vn(ni  ,m) - vn(ni-5,m) ) )              
        end do
        
    else
        
        ed  = ni + 2
        edc = ed
            
    end if
    
    do m=nst,ned
        do i=st-2,ed+2
            g(1,i,m) = vn(i  ,m) - vn(i-1,m)            
            g(2,i,m) = vn(i+1,m) - vn(i  ,m)
            g(3,i,m) = vn(i+2,m) - vn(i+1,m)            
        end do
    end do   
    
    do m=nst,ned
        do i=st-2,ed+2
            do n=1,3
                is = g(n,i,m)*g(n,i,m)
                eis2(n) = is**2
            end do

            tau7 = abs(eis2(2)-eis2(1))
            bl(1) = 1.0+tau7/(eis2(1)+eps)
            bl(2) = 1.0+tau7/(eis2(2)+eps)
            bs = bl(1) + bl(2)
            wl(2,i,m) = bl(2)/bs
            
            tau7 = abs(eis2(3)-eis2(2))
            br(2) = 1.0+tau7/(eis2(2)+eps)
            br(3) = 1.0+tau7/(eis2(3)+eps)
            bs = br(2) + br(3)
            wr(2,i,m) = br(2)/bs
        end do
    end do
    
    do i=st,ed
        imin=1.0
        imax=0.0
        do m=nst,ned
            imin         = min(wl(2,i-3,m),wl(2,i-2,m),wl(2,i-1,m),wl(2,i,m),wl(2,i+1,m),wl(2,i+2,m),wl(2,i+3,m), &
                               wr(2,i-3,m),wr(2,i-2,m),wr(2,i-1,m),wr(2,i,m),wr(2,i+1,m),wr(2,i+2,m),wr(2,i+3,m),imin)                                
            imax         = max(wl(2,i-3,m),wl(2,i-2,m),wl(2,i-1,m),wl(2,i,m),wl(2,i+1,m),wl(2,i+2,m),wl(2,i+3,m), &
                               wr(2,i-3,m),wr(2,i-2,m),wr(2,i-1,m),wr(2,i,m),wr(2,i+1,m),wr(2,i+2,m),wr(2,i+3,m),imax)
        end do
        indicator    = imin/imax
        if(indicator>indicatorc) then
            sigma(i)=1.0
        else
            sigma(i)=0.0
        endif        
    end do
    
    id=0
    do i=st,ed
        YN=0
        if(sigma(i)<0.5) then
            YN=YN+1
        end if
        if(YN>0) then
            id=id+1
            ind(id)=i           
        end if
    end do
       
    do i=1-ngn,ni+ngn
        sigma(i)=0.0
    end do
    if(ed-st<9) then
        ind_st=st
        ind_ed=ed
        nsec=0
        
        nsecn=1
        sec_stn(1)=st
        sec_edn(1)=ed
    else
        if(id==0) then
            nsec=1
            sec_st(1)=st
            sec_ed(1)=ed
            do n=sec_st(1),sec_ed(1)
                sigma(n)=1.0
            end do
            ind_st=st
            ind_ed=ed
            nsecn=0
        else
            ind_st=st
            ind_ed=ed
            nsec=0
            nsecn=0
            if(ind(1)-ind_st>9) then
                nsec=1
                sec_st(1)=ind_st
                sec_ed(1)=ind(1)-1
                pre_ed   =sec_ed(1)
                do n=sec_st(1),sec_ed(1)
                    sigma(n)=1.0
                end do
            else
                nsecn=1
                sec_stn(1)=ind_st
                sec_edn(1)=ind(1)
                pre_ed    =sec_edn(1)
            end if
            do i=2,id
                if(ind(i)-ind(i-1)>9) then
                    nsec=nsec+1
                    sec_st(nsec)=ind(i-1)+1
                    sec_ed(nsec)=ind(i)-1
                    do n=sec_st(nsec),sec_ed(nsec)
                        sigma(n)=1.0
                    end do
                    nsecn=nsecn+1
                    sec_stn(nsecn)=pre_ed+1
                    sec_edn(nsecn)=ind(i-1)
                    pre_ed        =sec_ed(nsec)
                else
                    if(ind(i)-ind(i-1)>1 .or. i==id) then
                        nsecn=nsecn+1
                        sec_stn(nsecn)=pre_ed+1
                        sec_edn(nsecn)=ind(i)
                        pre_ed        =sec_edn(nsecn)
                    end if
                end if
            end do
            if(ind_ed-ind(id)>9) then
                nsec=nsec+1
                sec_st(nsec)=pre_ed+1
                sec_ed(nsec)=ind_ed
                do n=sec_st(nsec),sec_ed(nsec)
                    sigma(n)=1.0
                end do
            else
                nsecn=nsecn+1
                sec_stn(nsecn)=pre_ed+1
                sec_edn(nsecn)=ind_ed
            end if
        end if
    end if
    
    do m=nst,ned
        do j=1,nsecn
            do i=sec_stn(j),sec_edn(j)
                g(1,i,m) = -23.0/24.0*vn(i-3,m) + 31.0/8.0*vn(i-2,m) - 47.0/8.0*vn(i-1,m) + 71.0/24.0*vn(i  ,m)
                g(2,i,m) =   1.0/24.0*vn(i-2,m) -  1.0/8.0*vn(i-1,m) -  7.0/8.0*vn(i  ,m) + 23.0/24.0*vn(i+1,m)
                g(3,i,m) =   1.0/24.0*vn(i-1,m) -  9.0/8.0*vn(i  ,m) +  9.0/8.0*vn(i+1,m) -  1.0/24.0*vn(i+2,m)
                g(4,i,m) = -23.0/24.0*vn(i  ,m) +  7.0/8.0*vn(i+1,m) +  1.0/8.0*vn(i+2,m) -  1.0/24.0*vn(i+3,m)
                g(5,i,m) = -71.0/24.0*vn(i+1,m) + 47.0/8.0*vn(i+2,m) - 31.0/8.0*vn(i+3,m) + 23.0/24.0*vn(i+4,m)
                
                s(1,i,m) =   -3.0/2.0*vn(i-3,m) + 11.0/2.0*vn(i-2,m) - 13.0/2.0*vn(i-1,m) +   5.0/2.0*vn(i  ,m)
                s(2,i,m) =   -1.0/2.0*vn(i-2,m) +  5.0/2.0*vn(i-1,m) -  7.0/2.0*vn(i  ,m) +   3.0/2.0*vn(i+1,m)
                s(3,i,m) =    1.0/2.0*vn(i-1,m) -  1.0/2.0*vn(i  ,m) -  1.0/2.0*vn(i+1,m) +   1.0/2.0*vn(i+2,m)
                s(4,i,m) =    3.0/2.0*vn(i  ,m) -  7.0/2.0*vn(i+1,m) +  5.0/2.0*vn(i+2,m) -   1.0/2.0*vn(i+3,m)
                s(5,i,m) =    5.0/2.0*vn(i+1,m) - 13.0/2.0*vn(i+2,m) + 11.0/2.0*vn(i+3,m) -   3.0/2.0*vn(i+4,m)
                
                t(1,i,m) =    -vn(i-3,m) + 3.0*vn(i-2,m) - 3.0*vn(i-1,m) + vn(i  ,m)
                t(2,i,m) =    -vn(i-2,m) + 3.0*vn(i-1,m) - 3.0*vn(i  ,m) + vn(i+1,m)
                t(3,i,m) =    -vn(i-1,m) + 3.0*vn(i  ,m) - 3.0*vn(i+1,m) + vn(i+2,m)
                t(4,i,m) =    -vn(i  ,m) + 3.0*vn(i+1,m) - 3.0*vn(i+2,m) + vn(i+3,m)
                t(5,i,m) =    -vn(i+1,m) + 3.0*vn(i+2,m) - 3.0*vn(i+3,m) + vn(i+4,m)
            end do
        end do      
    end do
    
    do j=1,nsecn
        do i=sec_stn(j),sec_edn(j)
            nx = (sn(i,1) + sn(i+1,1))/2.0
            ny = (sn(i,2) + sn(i+1,2))/2.0
            nz = (sn(i,3) + sn(i+1,3))/2.0

            do m=nst,ned
                vp(m) = (vn(i,m) + vn(i+1,m))/2.0
                fp(m) = vn(i,m)
            end do
            call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
            do m=nst,ned
                q(i,m) = vc(m)
            end do
            
            do m=nst,ned
                fp(m) = vn(i+1,m)
            end do
            call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
            do m=nst,ned
                q1(i+1,m) = vc(m)
            end do        
             
            do n=1,5
                do m=nst,ned
                    fp(m) = g(n,i,m)
                end do
                call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
                do m=nst,ned
                    g(n,i,m) = vc(m)
                end do

                do m=nst,ned
                    fp(m) = s(n,i,m)
                end do
                call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
                do m=nst,ned
                    s(n,i,m) = vc(m)
                end do
                
                do m=nst,ned
                    fp(m) = t(n,i,m)
                end do
                call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
                do m=nst,ned
                    t(n,i,m) = vc(m)
                end do            
            end do
        end do
    end do   
            
            
    do m=nst,ned
        do j=1,nsecn
            do i=sec_stn(j),sec_edn(j)
                do n=1,5
                    is = g(n,i,m)*g(n,i,m) + s(n,i,m)*s(n,i,m) + t(n,i,m)*t(n,i,m)              
                    eis2(n) = is**2
                end do
                
                tau7 = abs(eis2(4)-eis2(1))
                
                do n=1,4
                    bl(n) = cl(n)*(1.0+tau7/(eis2(n)+eps))
                end do
                
                bs = bl(1) + bl(2) + bl(3) + bl(4)
                do n=1,4
                    wl(n,i,m) = bl(n)/bs
                end do
          
                tau7 = abs(eis2(5)-eis2(2))           
                do n=2,5
                    br(n) = cr(n)*(1.0+tau7/(eis2(n)+eps))
                end do

                bs = br(2) + br(3) + br(4) + br(5)
                do n=2,5
                    wr(n,i,m) = br(n)/bs
                end do
            end do
        end do
    end do
    
!!    open(2,file='sigma.dat',form='formatted')
!!    write(2,*) 'variables="i" "sigmal" "sigmar"'
!!    do i=st,ed
!!
!!        write(2,200) i,sigmal(i,1),sigmar(i,1)
!!
!!    end do
!!
!!200 format(I4,1x,e12.5,1x,e12.5)
!!    close(2)
!!
!!    stop    
    
    do m=nst,ned
        do j=1,nsecn
            do i=sec_stn(j),sec_edn(j)
                vln(i,m) =  q(i  ,m) + ( wl(1,i,m)*( 24.0*g(1,i,m)-6.0*s(1,i,m)+t(1,i,m)) + &
                                         wl(2,i,m)*( 24.0*g(2,i,m)-6.0*s(2,i,m)+t(2,i,m)) + &
                                         wl(3,i,m)*( 24.0*g(3,i,m)-6.0*s(3,i,m)+t(3,i,m)) + &
                                         wl(4,i,m)*( 24.0*g(4,i,m)-6.0*s(4,i,m)+t(4,i,m)) )/48.0

                vrn(i,m) = q1(i+1,m) + ( wr(2,i,m)*(-24.0*g(2,i,m)-6.0*s(2,i,m)-t(2,i,m)) + &
                                         wr(3,i,m)*(-24.0*g(3,i,m)-6.0*s(3,i,m)-t(3,i,m)) + &
                                         wr(4,i,m)*(-24.0*g(4,i,m)-6.0*s(4,i,m)-t(4,i,m)) + &
                                         wr(5,i,m)*(-24.0*g(5,i,m)-6.0*s(5,i,m)-t(5,i,m)) )/48.0
            end do
        end do
    end do
    
    do j=1,nsecn
        do i=sec_stn(j),sec_edn(j)
            nx = (sn(i,1) + sn(i+1,1))/2.0
            ny = (sn(i,2) + sn(i+1,2))/2.0
            nz = (sn(i,3) + sn(i+1,3))/2.0

            do m=nst,ned
                vp(m) = (vn(i,m) + vn(i+1,m))/2.0
                vc(m) = vln(i,m)
            end do

            call RxQ_p(nst,ned,vp,nt,nx,ny,nz,vc,fp) 
            do m=nst,ned
                vln(i,m) = fp(m)
            end do
            
            do m=nst,ned
                vc(m) = vrn(i ,m)
            end do

            call RxQ_p(nst,ned,vp,nt,nx,ny,nz,vc,fp) 
            do m=nst,ned
                vrn(i,m) = fp(m)
            end do
        end do
    end do     
      
    do j=1,nsec
        
        call hdcs7ciN(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,ve_lil,ve_lir,sub_limit,chkid,chklim,ndec,sec_st(j),sec_ed(j))        
        
    end do
        
    do m=nst,ned
        do i=st,ed
            
               vl(i,m) = sigma(i)*ve_lil(i,m) + vln(i,m)*(1.0-sigma(i))
                         
               vr(i,m) = sigma(i)*ve_lir(i,m) + vrn(i,m)*(1.0-sigma(i))            

        end do
    end do  

    do i=stc,edc
        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vl(i,m) < chklim(m) ) then
                    vl(i,m) = vn(i  ,m)
                    ndec = ndec + 1
                end if

                if (vr(i,m) < chklim(m) ) then
                    vr(i,m) = vn(i+1,m)
                    ndec = ndec + 1
                end if
           end if
        end do
    end do

end subroutine dcsh7ci

subroutine dcsh7cicc(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_bound_scheme
    implicit none
    integer(kind_int), intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  , intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(in)  :: sn(1-ngn:ni+ngn,1:3)
    integer(kind_int), intent(in)  :: nfs,nfe
    integer(kind_int), intent(in)  :: nge
    real(kind_real)  , intent(out) :: vl(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(out) :: vr(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , external    :: sub_limit
    integer(kind_int), intent(in)  :: chkid(nst:ned)
    real(kind_real),   intent(in)  :: chklim(nst:ned)
    integer(kind_int), intent(out) :: ndec  
    real(kind_real)   :: vln(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: vrn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: ve_lil(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: ve_lir(1-ngn:ni+ngn,nst:ned) 
    real(kind_real)   :: nt,nx,ny,nz,fp(nst:ned)
    real(kind_real)   :: vp(nst:ned),vc(nst:ned)
    real(kind_real)   :: q(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: q1(1-ngn:ni+ngn,nst:ned)    
    real(kind_real)   :: s(1:5,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: g(1:5,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: t(1:5,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: wl(1:5,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: wr(1:5,1-ngn:ni+ngn,nst:ned)    
    real(kind_real)   :: aa,aa_scs
    real(kind_real)   :: sigma(1-ngn:ni+ngn),imax,imin,indicator
    real(kind_real), parameter :: cl(5) = (/1.D0/64.D0,21.D0/64.D0,35.D0/64.D0,7.D0/64.D0,0.D0/64.D0/)
    real(kind_real), parameter :: cr(5) = (/0.D0/64.D0,7.D0/64.D0,35.D0/64.D0,21.D0/64.D0,1.D0/64.D0/)  
    real(kind_real), parameter :: eps        = 1.0e-30
    real(kind_real), parameter :: indicatorc = 0.1
    real(kind_real)            :: is,eis2(5),bs,tau7,bl(5),br(5)
    integer(kind_int)          :: i,m,n,j,st,ed,stc,edc
    integer(kind_int)          :: id,YN,ind(1:ni+2*ngn),ind_st,ind_ed,pre_ed
    integer(kind_int)          :: nsec,sec_st(1:ni+2*ngn),sec_ed(1:ni+2*ngn)
    integer(kind_int)          :: nsecn,sec_stn(1:ni+2*ngn),sec_edn(1:ni+2*ngn)
    
    if (ni < 8) then
        call muscl2pvcc(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)
        return
    end if
    
    ndec = 0
    
    aa     = 0.01
    aa_scs = 0.001
    
    wl(:,:,:) = 0.25
    wr(:,:,:) = 0.25
    vln(:,:)    = 0.0
    ve_lil(:,:) = 0.0
    vrn(:,:)    = 0.0
    ve_lir(:,:) = 0.0
           
    if (nfs == nbc_bound_scheme) then
        
        st  = 5
        stc = 1

        do m=nst,ned
            vl(1,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                  ! third-order
            vl(2,m) = (  3.0*vn(1,m) +  6.0*vn(2,m) -      vn(3,m) ) / 8.0   ! third-order 
            vl(3,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 + &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )
            vr(1,m) = vl(1,m)
            vr(2,m) = vl(2,m)
            vr(3,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 - &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )
            vl(4,m) =   ( 150.0*( vn(4,m) + vn(3,m) )          -  &
                           25.0*( vn(5,m) + vn(2,m) )          +  &
                            3.0*( vn(6,m) + vn(1,m) ) ) /256.0 -  &
                      aa*( 10.0*( vn(4,m) - vn(3,m) )          -  &
                            5.0*( vn(5,m) - vn(2,m) )          +  &
                                ( vn(6,m) - vn(1,m) ) ) 
            
            vr(4,m) =   ( 150.0*( vn(4,m) + vn(3,m) )          -  &
                           25.0*( vn(5,m) + vn(2,m) )          +  &
                            3.0*( vn(6,m) + vn(1,m) ) ) /256.0 +  &
                      aa*( 10.0*( vn(4,m) - vn(3,m) )          -  &
                            5.0*( vn(5,m) - vn(2,m) )          +  &
                                ( vn(6,m) - vn(1,m) ) )          
        end do
        
    else
        
        st  = -1
        stc = st
     
    end if
     
    if (nfe == nbc_bound_scheme) then
        
        ed  = ni - 3
        edc = ni + 1

        do m=nst,ned
            vr(ni+1,m) = ( 15.0*vn(ni,m) - 10.0*vn(ni-1,m) +  3.0*vn(ni-2,m) ) / 8.0                ! third-order
            vl(ni  ,m) = (  3.0*vn(ni,m) +  6.0*vn(ni-1,m) -      vn(ni-2,m) ) / 8.0   ! third-order
            vl(ni-1,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 + &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) )
            vl(ni+1,m) = vr(ni+1,m)
            vr(ni  ,m) = vl(ni  ,m)
            vr(ni-1,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 - &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) ) 
            vl(ni-2,m) =   ( 150.0*( vn(ni-2,m) + vn(ni-3,m) )          -  &
                              25.0*( vn(ni-1,m) + vn(ni-4,m) )          +  &
                               3.0*( vn(ni  ,m) + vn(ni-5,m) ) ) /256.0 -  &
                         aa*( 10.0*( vn(ni-2,m) - vn(ni-3,m) )          -  &
                               5.0*( vn(ni-1,m) - vn(ni-4,m) )          +  &
                                   ( vn(ni  ,m) - vn(ni-5,m) ) ) 
            
            vr(ni-2,m) =   ( 150.0*( vn(ni-2,m) + vn(ni-3,m) )          -  &
                              25.0*( vn(ni-1,m) + vn(ni-4,m) )          +  &
                               3.0*( vn(ni  ,m) + vn(ni-5,m) ) ) /256.0 +  &
                         aa*( 10.0*( vn(ni-2,m) - vn(ni-3,m) )          -  &
                               5.0*( vn(ni-1,m) - vn(ni-4,m) )          +  &
                                   ( vn(ni  ,m) - vn(ni-5,m) ) )              
        end do
        
    else
        
        ed  = ni + 3
        edc = ed
            
    end if
    
    do m=nst,ned
        do i=st-3,ed+1
            g(1,i,m) = vn(i  ,m) - vn(i-1,m)            
            g(2,i,m) = vn(i+1,m) - vn(i  ,m)
            g(3,i,m) = vn(i+2,m) - vn(i+1,m)            
        end do
    end do   
    
    do m=nst,ned
        do i=st-3,ed+1
            do n=1,3
                is = g(n,i,m)*g(n,i,m)
                eis2(n) = is**2
            end do

            tau7 = abs(eis2(2)-eis2(1))
            bl(1) = 1.0+tau7/(eis2(1)+eps)
            bl(2) = 1.0+tau7/(eis2(2)+eps)
            bs = bl(1) + bl(2)
            wl(2,i,m) = bl(2)/bs
            
            tau7 = abs(eis2(3)-eis2(2))
            br(2) = 1.0+tau7/(eis2(2)+eps)
            br(3) = 1.0+tau7/(eis2(3)+eps)
            bs = br(2) + br(3)
            wr(2,i,m) = br(2)/bs
        end do
    end do
    
    do i=st-1,ed-1
        imin=1.0
        imax=0.0
        do m=nst,ned
            imin         = min(wl(2,i-3,m),wl(2,i-2,m),wl(2,i-1,m),wl(2,i,m),wl(2,i+1,m),wl(2,i+2,m),wl(2,i+3,m), &
                               wr(2,i-3,m),wr(2,i-2,m),wr(2,i-1,m),wr(2,i,m),wr(2,i+1,m),wr(2,i+2,m),wr(2,i+3,m),imin)                                
            imax         = max(wl(2,i-3,m),wl(2,i-2,m),wl(2,i-1,m),wl(2,i,m),wl(2,i+1,m),wl(2,i+2,m),wl(2,i+3,m), &
                               wr(2,i-3,m),wr(2,i-2,m),wr(2,i-1,m),wr(2,i,m),wr(2,i+1,m),wr(2,i+2,m),wr(2,i+3,m),imax)
        end do
        indicator    = imin/imax
        if(indicator>indicatorc) then
            sigma(i)=1.0
        else
            sigma(i)=0.0
        endif        
    end do
    
    id=0
    do i=st-1,ed-1
        YN=0
        if(sigma(i)<0.5) then
            YN=YN+1
        end if
        if(YN>0) then
            id=id+1
            ind(id)=i           
        end if
    end do
       
    do i=1-ngn,ni+ngn
        sigma(i)=0.0
    end do
    if(ed-st<9) then
        ind_st=st-1
        ind_ed=ed-1
        nsec=0
        
        nsecn=1
        sec_stn(1)=st-1
        sec_edn(1)=ed-1
    else
        if(id==0) then
            nsec=1
            sec_st(1)=st-1
            sec_ed(1)=ed-1
            do n=sec_st(1),sec_ed(1)
                sigma(n)=1.0
            end do
            ind_st=st-1
            ind_ed=ed-1
            nsecn=0
        else
            ind_st=st-1
            ind_ed=ed-1
            nsec=0
            nsecn=0
            if(ind(1)-ind_st>9) then
                nsec=1
                sec_st(1)=ind_st
                sec_ed(1)=ind(1)-1
                pre_ed   =sec_ed(1)
                do n=sec_st(1),sec_ed(1)
                    sigma(n)=1.0
                end do
            else
                nsecn=1
                sec_stn(1)=ind_st
                sec_edn(1)=ind(1)
                pre_ed    =sec_edn(1)
            end if
            do i=2,id
                if(ind(i)-ind(i-1)>9) then
                    nsec=nsec+1
                    sec_st(nsec)=ind(i-1)+1
                    sec_ed(nsec)=ind(i)-1
                    do n=sec_st(nsec),sec_ed(nsec)
                        sigma(n)=1.0
                    end do
                    nsecn=nsecn+1
                    sec_stn(nsecn)=pre_ed+1
                    sec_edn(nsecn)=ind(i-1)
                    pre_ed        =sec_ed(nsec)
                else
                    if(ind(i)-ind(i-1)>1 .or. i==id) then
                        nsecn=nsecn+1
                        sec_stn(nsecn)=pre_ed+1
                        sec_edn(nsecn)=ind(i)
                        pre_ed        =sec_edn(nsecn)
                    end if
                end if
            end do
            if(ind_ed-ind(id)>9) then
                nsec=nsec+1
                sec_st(nsec)=pre_ed+1
                sec_ed(nsec)=ind_ed
                do n=sec_st(nsec),sec_ed(nsec)
                    sigma(n)=1.0
                end do
            else
                nsecn=nsecn+1
                sec_stn(nsecn)=pre_ed+1
                sec_edn(nsecn)=ind_ed
            end if
        end if
    end if
    
    do m=nst,ned
        do j=1,nsecn
            do i=sec_stn(j),sec_edn(j)
                g(1,i,m) = -23.0/24.0*vn(i-3,m) + 31.0/8.0*vn(i-2,m) - 47.0/8.0*vn(i-1,m) + 71.0/24.0*vn(i  ,m)
                g(2,i,m) =   1.0/24.0*vn(i-2,m) -  1.0/8.0*vn(i-1,m) -  7.0/8.0*vn(i  ,m) + 23.0/24.0*vn(i+1,m)
                g(3,i,m) =   1.0/24.0*vn(i-1,m) -  9.0/8.0*vn(i  ,m) +  9.0/8.0*vn(i+1,m) -  1.0/24.0*vn(i+2,m)
                g(4,i,m) = -23.0/24.0*vn(i  ,m) +  7.0/8.0*vn(i+1,m) +  1.0/8.0*vn(i+2,m) -  1.0/24.0*vn(i+3,m)
                g(5,i,m) = -71.0/24.0*vn(i+1,m) + 47.0/8.0*vn(i+2,m) - 31.0/8.0*vn(i+3,m) + 23.0/24.0*vn(i+4,m)
                
                s(1,i,m) =   -3.0/2.0*vn(i-3,m) + 11.0/2.0*vn(i-2,m) - 13.0/2.0*vn(i-1,m) +   5.0/2.0*vn(i  ,m)
                s(2,i,m) =   -1.0/2.0*vn(i-2,m) +  5.0/2.0*vn(i-1,m) -  7.0/2.0*vn(i  ,m) +   3.0/2.0*vn(i+1,m)
                s(3,i,m) =    1.0/2.0*vn(i-1,m) -  1.0/2.0*vn(i  ,m) -  1.0/2.0*vn(i+1,m) +   1.0/2.0*vn(i+2,m)
                s(4,i,m) =    3.0/2.0*vn(i  ,m) -  7.0/2.0*vn(i+1,m) +  5.0/2.0*vn(i+2,m) -   1.0/2.0*vn(i+3,m)
                s(5,i,m) =    5.0/2.0*vn(i+1,m) - 13.0/2.0*vn(i+2,m) + 11.0/2.0*vn(i+3,m) -   3.0/2.0*vn(i+4,m)
                
                t(1,i,m) =    -vn(i-3,m) + 3.0*vn(i-2,m) - 3.0*vn(i-1,m) + vn(i  ,m)
                t(2,i,m) =    -vn(i-2,m) + 3.0*vn(i-1,m) - 3.0*vn(i  ,m) + vn(i+1,m)
                t(3,i,m) =    -vn(i-1,m) + 3.0*vn(i  ,m) - 3.0*vn(i+1,m) + vn(i+2,m)
                t(4,i,m) =    -vn(i  ,m) + 3.0*vn(i+1,m) - 3.0*vn(i+2,m) + vn(i+3,m)
                t(5,i,m) =    -vn(i+1,m) + 3.0*vn(i+2,m) - 3.0*vn(i+3,m) + vn(i+4,m)
            end do
        end do      
    end do
    
    do j=1,nsecn
        do i=sec_stn(j),sec_edn(j)
            nx = (sn(i,1) + sn(i+1,1))/2.0
            ny = (sn(i,2) + sn(i+1,2))/2.0
            nz = (sn(i,3) + sn(i+1,3))/2.0

            do m=nst,ned
                vp(m) = (vn(i,m) + vn(i+1,m))/2.0
                fp(m) = vn(i,m)
            end do
            call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
            do m=nst,ned
                q(i,m) = vc(m)
            end do
            
            do m=nst,ned
                fp(m) = vn(i+1,m)
            end do
            call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
            do m=nst,ned
                q1(i+1,m) = vc(m)
            end do        
             
            do n=1,5
                do m=nst,ned
                    fp(m) = g(n,i,m)
                end do
                call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
                do m=nst,ned
                    g(n,i,m) = vc(m)
                end do

                do m=nst,ned
                    fp(m) = s(n,i,m)
                end do
                call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
                do m=nst,ned
                    s(n,i,m) = vc(m)
                end do
                
                do m=nst,ned
                    fp(m) = t(n,i,m)
                end do
                call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
                do m=nst,ned
                    t(n,i,m) = vc(m)
                end do            
            end do
        end do
    end do   
            
            
    do m=nst,ned
        do j=1,nsecn
            do i=sec_stn(j),sec_edn(j)
                do n=1,5
                    is = g(n,i,m)*g(n,i,m) + s(n,i,m)*s(n,i,m) + t(n,i,m)*t(n,i,m)              
                    eis2(n) = is**2
                end do
                
                tau7 = abs(eis2(4)-eis2(1))
                
                do n=1,4
                    bl(n) = cl(n)*(1.0+tau7/(eis2(n)+eps))
                end do
                
                bs = bl(1) + bl(2) + bl(3) + bl(4)
                do n=1,4
                    wl(n,i,m) = bl(n)/bs
                end do
          
                tau7 = abs(eis2(5)-eis2(2))           
                do n=2,5
                    br(n) = cr(n)*(1.0+tau7/(eis2(n)+eps))
                end do

                bs = br(2) + br(3) + br(4) + br(5)
                do n=2,5
                    wr(n,i,m) = br(n)/bs
                end do
            end do
        end do
    end do
    
!!    open(2,file='sigma.dat',form='formatted')
!!    write(2,*) 'variables="i" "sigmal" "sigmar"'
!!    do i=st,ed
!!
!!        write(2,200) i,sigmal(i,1),sigmar(i,1)
!!
!!    end do
!!
!!200 format(I4,1x,e12.5,1x,e12.5)
!!    close(2)
!!
!!    stop    
    
    do m=nst,ned
        do j=1,nsecn
            do i=sec_stn(j),sec_edn(j)
                vln(i,m) =  q(i  ,m) + ( wl(1,i,m)*( 24.0*g(1,i,m)-6.0*s(1,i,m)+t(1,i,m)) + &
                                         wl(2,i,m)*( 24.0*g(2,i,m)-6.0*s(2,i,m)+t(2,i,m)) + &
                                         wl(3,i,m)*( 24.0*g(3,i,m)-6.0*s(3,i,m)+t(3,i,m)) + &
                                         wl(4,i,m)*( 24.0*g(4,i,m)-6.0*s(4,i,m)+t(4,i,m)) )/48.0

                vrn(i,m) = q1(i+1,m) + ( wr(2,i,m)*(-24.0*g(2,i,m)-6.0*s(2,i,m)-t(2,i,m)) + &
                                         wr(3,i,m)*(-24.0*g(3,i,m)-6.0*s(3,i,m)-t(3,i,m)) + &
                                         wr(4,i,m)*(-24.0*g(4,i,m)-6.0*s(4,i,m)-t(4,i,m)) + &
                                         wr(5,i,m)*(-24.0*g(5,i,m)-6.0*s(5,i,m)-t(5,i,m)) )/48.0
            end do
        end do
    end do
    
    do j=1,nsecn
        do i=sec_stn(j),sec_edn(j)
            nx = (sn(i,1) + sn(i+1,1))/2.0
            ny = (sn(i,2) + sn(i+1,2))/2.0
            nz = (sn(i,3) + sn(i+1,3))/2.0

            do m=nst,ned
                vp(m) = (vn(i,m) + vn(i+1,m))/2.0
                vc(m) = vln(i,m)
            end do

            call RxQ_p(nst,ned,vp,nt,nx,ny,nz,vc,fp) 
            do m=nst,ned
                vln(i,m) = fp(m)
            end do
            
            do m=nst,ned
                vc(m) = vrn(i ,m)
            end do

            call RxQ_p(nst,ned,vp,nt,nx,ny,nz,vc,fp) 
            do m=nst,ned
                vrn(i,m) = fp(m)
            end do
        end do
    end do     
      
    do j=1,nsec
        
        call hdcs7ciN(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,ve_lil,ve_lir,sub_limit,chkid,chklim,ndec,sec_st(j),sec_ed(j))        
        
    end do
        
    do m=nst,ned
        do i=st-1,ed-1
            
               vl(i+1,m) = sigma(i)*ve_lil(i,m) + vln(i,m)*(1.0-sigma(i))
                         
               vr(i+1,m) = sigma(i)*ve_lir(i,m) + vrn(i,m)*(1.0-sigma(i))            

        end do
    end do  

    do i=stc,edc
        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vl(i,m) < chklim(m) ) then
                    vl(i,m) = vn(i  ,m)
                    ndec = ndec + 1
                end if

                if (vr(i,m) < chklim(m) ) then
                    vr(i,m) = vn(i+1,m)
                    ndec = ndec + 1
                end if
           end if
        end do
    end do

end subroutine dcsh7cicc

subroutine scsn4ci(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_bound_scheme
    implicit none
    integer(kind_int), intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  , intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(in)  :: sn(1-ngn:ni+ngn,1:3)
    integer(kind_int), intent(in)  :: nfs,nfe
    integer(kind_int), intent(in)  :: nge
    real(kind_real)  , intent(out) :: vl(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(out) :: vr(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , external    :: sub_limit
    integer(kind_int), intent(in)  :: chkid(nst:ned)
    real(kind_real),   intent(in)  :: chklim(nst:ned)
    integer(kind_int), intent(out) :: ndec   
    real(kind_real)   :: nt,nx,ny,nz,fp(nst:ned)
    real(kind_real)   :: vp(nst:ned),vc(nst:ned)
    real(kind_real)   :: q(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: q1(1-ngn:ni+ngn,nst:ned)    
    real(kind_real)   :: s(1:5,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: g(1:5,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: t(1:5,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: wl(1:5,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: wr(1:5,1-ngn:ni+ngn,nst:ned)   
    real(kind_real)   :: aa
    real(kind_real), parameter :: cl(5) = (/1.D0/64.D0,21.D0/64.D0,35.D0/64.D0,7.D0/64.D0,0.D0/64.D0/)
    real(kind_real), parameter :: cr(5) = (/0.D0/64.D0,7.D0/64.D0,35.D0/64.D0,21.D0/64.D0,1.D0/64.D0/)  
    real(kind_real), parameter :: eps        = 1.0e-30
    real(kind_real)            :: is,eis2(5),bs,tau7,bl(5),br(5)
    integer(kind_int)          :: i,m,n,st,ed,stc,edc
    
    if (ni < 8) then
        call muscl2pv(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)
        return
    end if
    
    ndec = 0
    
    aa     = 0.01

    wl(:,:,:) = 0.25
    wr(:,:,:) = 0.25

    if (nfs == nbc_bound_scheme) then
        
        st  = 4
        stc = 0

        do m=nst,ned
            vl(0,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                  ! third-order
            vl(1,m) = (  3.0*vn(1,m) +  6.0*vn(2,m) -      vn(3,m) ) / 8.0   ! third-order 
            vl(2,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 + &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )
            vr(0,m) = vl(0,m)
            vr(1,m) = vl(1,m)
            vr(2,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 - &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )
            vl(3,m) =   ( 150.0*( vn(4,m) + vn(3,m) )          -  &
                           25.0*( vn(5,m) + vn(2,m) )          +  &
                            3.0*( vn(6,m) + vn(1,m) ) ) /256.0 -  &
                      aa*( 10.0*( vn(4,m) - vn(3,m) )          -  &
                            5.0*( vn(5,m) - vn(2,m) )          +  &
                                ( vn(6,m) - vn(1,m) ) ) 
            
            vr(3,m) =   ( 150.0*( vn(4,m) + vn(3,m) )          -  &
                           25.0*( vn(5,m) + vn(2,m) )          +  &
                            3.0*( vn(6,m) + vn(1,m) ) ) /256.0 +  &
                      aa*( 10.0*( vn(4,m) - vn(3,m) )          -  &
                            5.0*( vn(5,m) - vn(2,m) )          +  &
                                ( vn(6,m) - vn(1,m) ) )          
        end do
        
    else
        
        st  = -2
        stc = st
     
    end if
     
    if (nfe == nbc_bound_scheme) then
        
        ed  = ni - 4
        edc = ni

        do m=nst,ned
            vr(ni  ,m) = ( 15.0*vn(ni,m) - 10.0*vn(ni-1,m) +  3.0*vn(ni-2,m) ) / 8.0                ! third-order
            vl(ni-1,m) = (  3.0*vn(ni,m) +  6.0*vn(ni-1,m) -      vn(ni-2,m) ) / 8.0   ! third-order
            vl(ni-2,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 + &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) )
            vl(ni  ,m) = vr(ni  ,m)
            vr(ni-1,m) = vl(ni-1,m)
            vr(ni-2,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 - &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) ) 
            vl(ni-3,m) =   ( 150.0*( vn(ni-2,m) + vn(ni-3,m) )          -  &
                              25.0*( vn(ni-1,m) + vn(ni-4,m) )          +  &
                               3.0*( vn(ni  ,m) + vn(ni-5,m) ) ) /256.0 -  &
                         aa*( 10.0*( vn(ni-2,m) - vn(ni-3,m) )          -  &
                               5.0*( vn(ni-1,m) - vn(ni-4,m) )          +  &
                                   ( vn(ni  ,m) - vn(ni-5,m) ) ) 
            
            vr(ni-3,m) =   ( 150.0*( vn(ni-2,m) + vn(ni-3,m) )          -  &
                              25.0*( vn(ni-1,m) + vn(ni-4,m) )          +  &
                               3.0*( vn(ni  ,m) + vn(ni-5,m) ) ) /256.0 +  &
                         aa*( 10.0*( vn(ni-2,m) - vn(ni-3,m) )          -  &
                               5.0*( vn(ni-1,m) - vn(ni-4,m) )          +  &
                                   ( vn(ni  ,m) - vn(ni-5,m) ) )              
        end do
        
    else
        
        ed  = ni + 2
        edc = ed
            
    end if
    
    do m=nst,ned

        do i=st,ed
            g(1,i,m) = -23.0/24.0*vn(i-3,m) + 31.0/8.0*vn(i-2,m) - 47.0/8.0*vn(i-1,m) + 71.0/24.0*vn(i  ,m)
            g(2,i,m) =   1.0/24.0*vn(i-2,m) -  1.0/8.0*vn(i-1,m) -  7.0/8.0*vn(i  ,m) + 23.0/24.0*vn(i+1,m)
            g(3,i,m) =   1.0/24.0*vn(i-1,m) -  9.0/8.0*vn(i  ,m) +  9.0/8.0*vn(i+1,m) -  1.0/24.0*vn(i+2,m)
            g(4,i,m) = -23.0/24.0*vn(i  ,m) +  7.0/8.0*vn(i+1,m) +  1.0/8.0*vn(i+2,m) -  1.0/24.0*vn(i+3,m)
            g(5,i,m) = -71.0/24.0*vn(i+1,m) + 47.0/8.0*vn(i+2,m) - 31.0/8.0*vn(i+3,m) + 23.0/24.0*vn(i+4,m)
            
            s(1,i,m) =   -3.0/2.0*vn(i-3,m) + 11.0/2.0*vn(i-2,m) - 13.0/2.0*vn(i-1,m) +   5.0/2.0*vn(i  ,m)
            s(2,i,m) =   -1.0/2.0*vn(i-2,m) +  5.0/2.0*vn(i-1,m) -  7.0/2.0*vn(i  ,m) +   3.0/2.0*vn(i+1,m)
            s(3,i,m) =    1.0/2.0*vn(i-1,m) -  1.0/2.0*vn(i  ,m) -  1.0/2.0*vn(i+1,m) +   1.0/2.0*vn(i+2,m)
            s(4,i,m) =    3.0/2.0*vn(i  ,m) -  7.0/2.0*vn(i+1,m) +  5.0/2.0*vn(i+2,m) -   1.0/2.0*vn(i+3,m)
            s(5,i,m) =    5.0/2.0*vn(i+1,m) - 13.0/2.0*vn(i+2,m) + 11.0/2.0*vn(i+3,m) -   3.0/2.0*vn(i+4,m)
            
            t(1,i,m) =    -vn(i-3,m) + 3.0*vn(i-2,m) - 3.0*vn(i-1,m) + vn(i  ,m)
            t(2,i,m) =    -vn(i-2,m) + 3.0*vn(i-1,m) - 3.0*vn(i  ,m) + vn(i+1,m)
            t(3,i,m) =    -vn(i-1,m) + 3.0*vn(i  ,m) - 3.0*vn(i+1,m) + vn(i+2,m)
            t(4,i,m) =    -vn(i  ,m) + 3.0*vn(i+1,m) - 3.0*vn(i+2,m) + vn(i+3,m)
            t(5,i,m) =    -vn(i+1,m) + 3.0*vn(i+2,m) - 3.0*vn(i+3,m) + vn(i+4,m)
        end do      
    end do
    
    do i=st,ed
        nx = (sn(i,1) + sn(i+1,1))/2.0
        ny = (sn(i,2) + sn(i+1,2))/2.0
        nz = (sn(i,3) + sn(i+1,3))/2.0

        do m=nst,ned
            vp(m) = (vn(i,m) + vn(i+1,m))/2.0
            fp(m) = vn(i,m)
        end do
        call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
        do m=nst,ned
            q(i,m) = vc(m)
        end do
        
        do m=nst,ned
            fp(m) = vn(i+1,m)
        end do
        call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
        do m=nst,ned
            q1(i+1,m) = vc(m)
        end do        
         
        do n=1,5
            do m=nst,ned
                fp(m) = g(n,i,m)
            end do
            call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
            do m=nst,ned
                g(n,i,m) = vc(m)
            end do

            do m=nst,ned
                fp(m) = s(n,i,m)
            end do
            call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
            do m=nst,ned
                s(n,i,m) = vc(m)
            end do
            
            do m=nst,ned
                fp(m) = t(n,i,m)
            end do
            call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
            do m=nst,ned
                t(n,i,m) = vc(m)
            end do            
        end do
    end do  
            
            
    do m=nst,ned
        do i=st,ed
            do n=1,5
                is = g(n,i,m)*g(n,i,m) + s(n,i,m)*s(n,i,m) + t(n,i,m)*t(n,i,m)              
                eis2(n) = is**2
            end do
            
            tau7 = abs(eis2(4)-eis2(1))
            
            do n=1,4
                bl(n) = cl(n)*(1.0+tau7/(eis2(n)+eps))
            end do
            
            bs = bl(1) + bl(2) + bl(3) + bl(4)
            do n=1,4
                wl(n,i,m) = bl(n)/bs
            end do
          
            tau7 = abs(eis2(5)-eis2(2))           
            do n=2,5
                br(n) = cr(n)*(1.0+tau7/(eis2(n)+eps))
            end do

            bs = br(2) + br(3) + br(4) + br(5)
            do n=2,5
                wr(n,i,m) = br(n)/bs
            end do
        end do
    end do
    
    do m=nst,ned
        do i=st,ed
            vl(i,m) =  q(i  ,m) + ( wl(1,i,m)*( 24.0*g(1,i,m)-6.0*s(1,i,m)+t(1,i,m)) + &
                                    wl(2,i,m)*( 24.0*g(2,i,m)-6.0*s(2,i,m)+t(2,i,m)) + &
                                    wl(3,i,m)*( 24.0*g(3,i,m)-6.0*s(3,i,m)+t(3,i,m)) + &
                                    wl(4,i,m)*( 24.0*g(4,i,m)-6.0*s(4,i,m)+t(4,i,m)) )/48.0
        end do

        do i=st,ed
            vr(i,m) = q1(i+1,m) + ( wr(2,i,m)*(-24.0*g(2,i,m)-6.0*s(2,i,m)-t(2,i,m)) + &
                                    wr(3,i,m)*(-24.0*g(3,i,m)-6.0*s(3,i,m)-t(3,i,m)) + &
                                    wr(4,i,m)*(-24.0*g(4,i,m)-6.0*s(4,i,m)-t(4,i,m)) + &
                                    wr(5,i,m)*(-24.0*g(5,i,m)-6.0*s(5,i,m)-t(5,i,m)) )/48.0
        end do
    end do   
    
    do i=st,ed
        nx = (sn(i,1) + sn(i+1,1))/2.0
        ny = (sn(i,2) + sn(i+1,2))/2.0
        nz = (sn(i,3) + sn(i+1,3))/2.0

        do m=nst,ned
            vp(m) = (vn(i,m) + vn(i+1,m))/2.0
            vc(m) = vl(i,m)
        end do

        call RxQ_p(nst,ned,vp,nt,nx,ny,nz,vc,fp) 
        do m=nst,ned
            vl(i,m) = fp(m)
        end do
        
        do m=nst,ned
            vc(m) = vr(i ,m)
        end do

        call RxQ_p(nst,ned,vp,nt,nx,ny,nz,vc,fp) 
        do m=nst,ned
            vr(i,m) = fp(m)
        end do        
    end do

    do i=stc,edc
        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vl(i,m) < chklim(m) ) then
                    vl(i,m) = vn(i  ,m)
                    ndec = ndec + 1
                end if

                if (vr(i,m) < chklim(m) ) then
                    vr(i,m) = vn(i+1,m)
                    ndec = ndec + 1
                end if
           end if
        end do
    end do

end subroutine scsn4ci

subroutine scsn4cicc(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_bound_scheme
    implicit none
    integer(kind_int), intent(in)  :: ni,nst,ned,ngn
    real(kind_real)  , intent(in)  :: vn(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(in)  :: sn(1-ngn:ni+ngn,1:3)
    integer(kind_int), intent(in)  :: nfs,nfe
    integer(kind_int), intent(in)  :: nge
    real(kind_real)  , intent(out) :: vl(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , intent(out) :: vr(1-ngn:ni+ngn,nst:ned)
    real(kind_real)  , external    :: sub_limit
    integer(kind_int), intent(in)  :: chkid(nst:ned)
    real(kind_real),   intent(in)  :: chklim(nst:ned)
    integer(kind_int), intent(out) :: ndec   
    real(kind_real)   :: nt,nx,ny,nz,fp(nst:ned)
    real(kind_real)   :: vp(nst:ned),vc(nst:ned)
    real(kind_real)   :: q(1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: q1(1-ngn:ni+ngn,nst:ned)    
    real(kind_real)   :: s(1:5,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: g(1:5,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: t(1:5,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: wl(1:5,1-ngn:ni+ngn,nst:ned)
    real(kind_real)   :: wr(1:5,1-ngn:ni+ngn,nst:ned)   
    real(kind_real)   :: aa
    real(kind_real), parameter :: cl(5) = (/1.D0/64.D0,21.D0/64.D0,35.D0/64.D0,7.D0/64.D0,0.D0/64.D0/)
    real(kind_real), parameter :: cr(5) = (/0.D0/64.D0,7.D0/64.D0,35.D0/64.D0,21.D0/64.D0,1.D0/64.D0/)  
    real(kind_real), parameter :: eps        = 1.0e-30
    real(kind_real)            :: is,eis2(5),bs,tau7,bl(5),br(5)
    integer(kind_int)          :: i,m,n,st,ed,stc,edc
    
    if (ni < 8) then
        call muscl2pvcc(ni,nst,ned,ngn,vn,sn,nfs,nfe,nge,vl,vr,sub_limit,chkid,chklim,ndec)
        return
    end if
    
    ndec = 0
    
    aa     = 0.01

    wl(:,:,:) = 0.25
    wr(:,:,:) = 0.25

    if (nfs == nbc_bound_scheme) then
        
        st  = 5
        stc = 1

        do m=nst,ned
            vl(1,m) = ( 15.0*vn(1,m) - 10.0*vn(2,m) +  3.0*vn(3,m) ) / 8.0                  ! third-order
            vl(2,m) = (  3.0*vn(1,m) +  6.0*vn(2,m) -      vn(3,m) ) / 8.0   ! third-order 
            vl(3,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 + &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )
            vr(1,m) = vl(1,m)
            vr(2,m) = vl(2,m)
            vr(3,m) = (    -vn(1,m) + 9.0*vn(2,m)   +   9.0*vn(3,m) -       vn(4,m) ) /16.0 - &
                      aa*10.0*(-3.0*( vn(3,m) - vn(2,m) )+ vn(4,m) - vn(1,m) )
            vl(4,m) =   ( 150.0*( vn(4,m) + vn(3,m) )          -  &
                           25.0*( vn(5,m) + vn(2,m) )          +  &
                            3.0*( vn(6,m) + vn(1,m) ) ) /256.0 -  &
                      aa*( 10.0*( vn(4,m) - vn(3,m) )          -  &
                            5.0*( vn(5,m) - vn(2,m) )          +  &
                                ( vn(6,m) - vn(1,m) ) ) 
            
            vr(4,m) =   ( 150.0*( vn(4,m) + vn(3,m) )          -  &
                           25.0*( vn(5,m) + vn(2,m) )          +  &
                            3.0*( vn(6,m) + vn(1,m) ) ) /256.0 +  &
                      aa*( 10.0*( vn(4,m) - vn(3,m) )          -  &
                            5.0*( vn(5,m) - vn(2,m) )          +  &
                                ( vn(6,m) - vn(1,m) ) )          
        end do
        
    else
        
        st  = -1
        stc = st
     
    end if
     
    if (nfe == nbc_bound_scheme) then
        
        ed  = ni - 3
        edc = ni + 1

        do m=nst,ned
            vr(ni+1,m) = ( 15.0*vn(ni,m) - 10.0*vn(ni-1,m) +  3.0*vn(ni-2,m) ) / 8.0                ! third-order
            vl(ni  ,m) = (  3.0*vn(ni,m) +  6.0*vn(ni-1,m) -      vn(ni-2,m) ) / 8.0   ! third-order
            vl(ni-1,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 + &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) )
            vl(ni+1,m) = vr(ni+1,m)
            vr(ni  ,m) = vl(ni  ,m)
            vr(ni-1,m) = (     -vn(ni-3,m) +   9.0*vn(ni-2,m) +   9.0*vn(ni-1,m)        -vn(ni,m) ) /16.0 - &
                         aa*10.0*(-3.0*( vn(ni-1,m) - vn(ni-2,m) )+ vn(ni,m) - vn(ni-3,m) ) 
            vl(ni-2,m) =   ( 150.0*( vn(ni-2,m) + vn(ni-3,m) )          -  &
                              25.0*( vn(ni-1,m) + vn(ni-4,m) )          +  &
                               3.0*( vn(ni  ,m) + vn(ni-5,m) ) ) /256.0 -  &
                         aa*( 10.0*( vn(ni-2,m) - vn(ni-3,m) )          -  &
                               5.0*( vn(ni-1,m) - vn(ni-4,m) )          +  &
                                   ( vn(ni  ,m) - vn(ni-5,m) ) ) 
            
            vr(ni-2,m) =   ( 150.0*( vn(ni-2,m) + vn(ni-3,m) )          -  &
                              25.0*( vn(ni-1,m) + vn(ni-4,m) )          +  &
                               3.0*( vn(ni  ,m) + vn(ni-5,m) ) ) /256.0 +  &
                         aa*( 10.0*( vn(ni-2,m) - vn(ni-3,m) )          -  &
                               5.0*( vn(ni-1,m) - vn(ni-4,m) )          +  &
                                   ( vn(ni  ,m) - vn(ni-5,m) ) )              
        end do
        
    else
        
        ed  = ni + 3
        edc = ed
            
    end if
    
    do m=nst,ned

        do i=st-1,ed-1
            g(1,i,m) = -23.0/24.0*vn(i-3,m) + 31.0/8.0*vn(i-2,m) - 47.0/8.0*vn(i-1,m) + 71.0/24.0*vn(i  ,m)
            g(2,i,m) =   1.0/24.0*vn(i-2,m) -  1.0/8.0*vn(i-1,m) -  7.0/8.0*vn(i  ,m) + 23.0/24.0*vn(i+1,m)
            g(3,i,m) =   1.0/24.0*vn(i-1,m) -  9.0/8.0*vn(i  ,m) +  9.0/8.0*vn(i+1,m) -  1.0/24.0*vn(i+2,m)
            g(4,i,m) = -23.0/24.0*vn(i  ,m) +  7.0/8.0*vn(i+1,m) +  1.0/8.0*vn(i+2,m) -  1.0/24.0*vn(i+3,m)
            g(5,i,m) = -71.0/24.0*vn(i+1,m) + 47.0/8.0*vn(i+2,m) - 31.0/8.0*vn(i+3,m) + 23.0/24.0*vn(i+4,m)
            
            s(1,i,m) =   -3.0/2.0*vn(i-3,m) + 11.0/2.0*vn(i-2,m) - 13.0/2.0*vn(i-1,m) +   5.0/2.0*vn(i  ,m)
            s(2,i,m) =   -1.0/2.0*vn(i-2,m) +  5.0/2.0*vn(i-1,m) -  7.0/2.0*vn(i  ,m) +   3.0/2.0*vn(i+1,m)
            s(3,i,m) =    1.0/2.0*vn(i-1,m) -  1.0/2.0*vn(i  ,m) -  1.0/2.0*vn(i+1,m) +   1.0/2.0*vn(i+2,m)
            s(4,i,m) =    3.0/2.0*vn(i  ,m) -  7.0/2.0*vn(i+1,m) +  5.0/2.0*vn(i+2,m) -   1.0/2.0*vn(i+3,m)
            s(5,i,m) =    5.0/2.0*vn(i+1,m) - 13.0/2.0*vn(i+2,m) + 11.0/2.0*vn(i+3,m) -   3.0/2.0*vn(i+4,m)
            
            t(1,i,m) =    -vn(i-3,m) + 3.0*vn(i-2,m) - 3.0*vn(i-1,m) + vn(i  ,m)
            t(2,i,m) =    -vn(i-2,m) + 3.0*vn(i-1,m) - 3.0*vn(i  ,m) + vn(i+1,m)
            t(3,i,m) =    -vn(i-1,m) + 3.0*vn(i  ,m) - 3.0*vn(i+1,m) + vn(i+2,m)
            t(4,i,m) =    -vn(i  ,m) + 3.0*vn(i+1,m) - 3.0*vn(i+2,m) + vn(i+3,m)
            t(5,i,m) =    -vn(i+1,m) + 3.0*vn(i+2,m) - 3.0*vn(i+3,m) + vn(i+4,m)
        end do      
    end do
    
    do i=st-1,ed-1
        nx = (sn(i,1) + sn(i+1,1))/2.0
        ny = (sn(i,2) + sn(i+1,2))/2.0
        nz = (sn(i,3) + sn(i+1,3))/2.0

        do m=nst,ned
            vp(m) = (vn(i,m) + vn(i+1,m))/2.0
            fp(m) = vn(i,m)
        end do
        call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
        do m=nst,ned
            q(i,m) = vc(m)
        end do
        
        do m=nst,ned
            fp(m) = vn(i+1,m)
        end do
        call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
        do m=nst,ned
            q1(i+1,m) = vc(m)
        end do        
         
        do n=1,5
            do m=nst,ned
                fp(m) = g(n,i,m)
            end do
            call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
            do m=nst,ned
                g(n,i,m) = vc(m)
            end do

            do m=nst,ned
                fp(m) = s(n,i,m)
            end do
            call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
            do m=nst,ned
                s(n,i,m) = vc(m)
            end do
            
            do m=nst,ned
                fp(m) = t(n,i,m)
            end do
            call LxQ_p(nst,ned,vp,nt,nx,ny,nz,fp,vc) 
            do m=nst,ned
                t(n,i,m) = vc(m)
            end do            
        end do
    end do  
            
            
    do m=nst,ned
        do i=st-1,ed-1
            do n=1,5
                is = g(n,i,m)*g(n,i,m) + s(n,i,m)*s(n,i,m) + t(n,i,m)*t(n,i,m)              
                eis2(n) = is**2
            end do
            
            tau7 = abs(eis2(4)-eis2(1))
            
            do n=1,4
                bl(n) = cl(n)*(1.0+tau7/(eis2(n)+eps))
            end do
            
            bs = bl(1) + bl(2) + bl(3) + bl(4)
            do n=1,4
                wl(n,i,m) = bl(n)/bs
            end do
          
            tau7 = abs(eis2(5)-eis2(2))           
            do n=2,5
                br(n) = cr(n)*(1.0+tau7/(eis2(n)+eps))
            end do

            bs = br(2) + br(3) + br(4) + br(5)
            do n=2,5
                wr(n,i,m) = br(n)/bs
            end do
        end do
    end do
    
    do m=nst,ned
        do i=st-1,ed-1
            vl(i+1,m) =  q(i  ,m) + ( wl(1,i,m)*( 24.0*g(1,i,m)-6.0*s(1,i,m)+t(1,i,m)) + &
                                      wl(2,i,m)*( 24.0*g(2,i,m)-6.0*s(2,i,m)+t(2,i,m)) + &
                                      wl(3,i,m)*( 24.0*g(3,i,m)-6.0*s(3,i,m)+t(3,i,m)) + &
                                      wl(4,i,m)*( 24.0*g(4,i,m)-6.0*s(4,i,m)+t(4,i,m)) )/48.0
        end do

        do i=st-1,ed-1
            vr(i+1,m) = q1(i+1,m) + ( wr(2,i,m)*(-24.0*g(2,i,m)-6.0*s(2,i,m)-t(2,i,m)) + &
                                      wr(3,i,m)*(-24.0*g(3,i,m)-6.0*s(3,i,m)-t(3,i,m)) + &
                                      wr(4,i,m)*(-24.0*g(4,i,m)-6.0*s(4,i,m)-t(4,i,m)) + &
                                      wr(5,i,m)*(-24.0*g(5,i,m)-6.0*s(5,i,m)-t(5,i,m)) )/48.0
        end do
    end do   
    
    do i=st,ed
        nx = (sn(i,1) + sn(i+1,1))/2.0
        ny = (sn(i,2) + sn(i+1,2))/2.0
        nz = (sn(i,3) + sn(i+1,3))/2.0

        do m=nst,ned
            vp(m) = (vn(i,m) + vn(i+1,m))/2.0
            vc(m) = vl(i,m)
        end do

        call RxQ_p(nst,ned,vp,nt,nx,ny,nz,vc,fp) 
        do m=nst,ned
            vl(i,m) = fp(m)
        end do
        
        do m=nst,ned
            vc(m) = vr(i ,m)
        end do

        call RxQ_p(nst,ned,vp,nt,nx,ny,nz,vc,fp) 
        do m=nst,ned
            vr(i,m) = fp(m)
        end do        
    end do

    do i=stc,edc
        do m=nst,ned
            if (chkid(m) /= 0) then
                if (vl(i,m) < chklim(m) ) then
                    vl(i,m) = vn(i  ,m)
                    ndec = ndec + 1
                end if

                if (vr(i,m) < chklim(m) ) then
                    vr(i,m) = vn(i+1,m)
                    ndec = ndec + 1
                end if
           end if
        end do
    end do

end subroutine scsn4cicc

subroutine dn_via_node6(ni,nst,ned,ngn,vn,nge,ve,nfs,nfe,dn)!jy
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_inter_scheme
    implicit none
    integer(kind_int),intent(in)    :: ni,nst,ned,ngn
    real(kind_real)  ,intent(in)    :: vn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)    :: nge
    real(kind_real)  ,intent(in)    :: ve(1-ngn:ni+ngn,nst:ned)
    integer(kind_int),intent(in)    :: nfs,nfe
    real(kind_real)  ,intent(inout) :: dn(1-ngn:ni+ngn,nst:ned)
    integer(kind_int) :: i,m,st,ed
    real(kind_real)   :: a,b,c

    a  = 192.0/256.0
    b  = -9.0/64.0
    c  = 64.0/3840.0
    
    if (nfs == nbc_inter_scheme) then
        st = 1
    else
        st = 4

        do m=nst,ned
            dn(3,m) = (  3.0*vn(1,m) - 30.0*vn(2,m) -  20.0*vn(3,m) + 60.0*vn(4,m) - 15.0*vn(5,m) + 2.0*vn(6,m)) / 60.0
            dn(2,m) = (-12.0*vn(1,m) - 65.0*vn(2,m) + 120.0*vn(3,m) - 60.0*vn(4,m) + 20.0*vn(5,m) - 3.0*vn(6,m)) / 60.0
            dn(1,m) = (-25.0*vn(1,m) + 48.0*vn(2,m) -  36.0*vn(3,m) + 16.0*vn(4,m) -  3.0*vn(5,m) ) / 12.0
        end do
    end if

    if (nfe == nbc_inter_scheme) then
        ed = ni
    else
        ed = ni - 3

        do m=nst,ned
            dn(ni-2,m) = &
            -(  3.0*vn(ni,m) - 30.0*vn(ni-1,m) -  20.0*vn(ni-2,m) + 60.0*vn(ni-3,m) - 15.0*vn(ni-4,m) + 2.0*vn(ni-5,m)) / 60.0
            dn(ni-1,m) = &
            -(-12.0*vn(ni,m) - 65.0*vn(ni-1,m) + 120.0*vn(ni-2,m) - 60.0*vn(ni-3,m) + 20.0*vn(ni-4,m) - 3.0*vn(ni-5,m)) / 60.0
            dn(ni  ,m) = &
            -(-25.0*vn(ni,m) + 48.0*vn(ni-1,m) -  36.0*vn(ni-2,m) + 16.0*vn(ni-3,m) -  3.0*vn(ni-4,m) ) / 12.0
        end do
    end if

    do m=nst,ned
        do i=st,ed
            dn(i,m) = a *(vn(i+1,m) - vn(i-1,m)) + &
                      b *(vn(i+2,m) - vn(i-2,m)) + &
                      c *(vn(i+3,m) - vn(i-3,m))
        end do
    end do
       
end subroutine dn_via_node6

subroutine Tri_pursue(A,B,C,D,S,us,ue,vst,ved)!jy
    use mod_kndconsts, only : kind_int,kind_real
    implicit none
    integer(kind_int), intent(in)  :: us,ue,vst,ved  
    real(kind_real),   intent(in)  :: A(vst:ved)
    real(kind_real),   intent(in)  :: B(vst:ved)
    real(kind_real),   intent(in)  :: C(vst:ved)
    real(kind_real),   intent(in)  :: D(vst:ved)
    real(kind_real),   intent(out) :: S(vst:ved)
    integer(kind_int) :: i   
    real(kind_real) :: M(vst:ved)
    real(kind_real) :: N(vst:ved)
    
    M(us) = -C(us)/B(us)
    N(us) = D(us)/B(us)
    do i=us+1,ue-1
        M(i) = -C(i)/(B(i)+A(i)*M(i-1))
        N(i) = (D(i)-A(i)*N(i-1))/(B(i)+A(i)*M(i-1))
    end do
    
    S(ue) = (D(ue)-A(ue)*N(ue-1))/(A(ue)*M(ue-1)+B(ue))
    do i=0,ue-us-1
        S(ue-1-i) = M(ue-1-i)*S(ue-i)+N(ue-1-i)
    end do

end subroutine Tri_pursue

