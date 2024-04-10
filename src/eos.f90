
function sutherland0(t,mu0,t0,ts) result(vis)
    use mod_kndconsts, only : kind_real
    implicit none
    real(kind_real) :: t,mu0,t0,ts
    real(kind_real) :: vis
    real(kind_real) :: rat

    rat = t/t0
    vis = mu0*rat*sqrt(rat)*(t0+ts)/(t+ts)

end function sutherland0

function sutherland(t,ts) result(vis)
    use mod_kndconsts, only : kind_real
    use mod_constants, only : one
    implicit none
    real(kind_real) :: t,ts
    real(kind_real) :: vis

    vis= t*sqrt(t)*(one+ts)/(t+ts)

end function sutherland

subroutine stdatm1976(alt, sigma, delta, theta)
    !   -------------------------------------------------------------------------
    ! purpose: compute the properties of the 1976 standard atmosphere to 86 km.
    ! author : Ralph Carmichael, Public Domain Aeronautical Software
    ! note   : if alt > 86, the values returned will not be correct, but they will
    !   not be too far removed from the correct values for density.
    !   the reference document does not use the terms pressure and temperature
    !   above 86 km.
    use mod_kndconsts, only : kind_int,kind_real=>kind_double
    implicit none
    real(kind_real), intent(in)  :: alt        ! geometric altitude, km.
    real(kind_real), intent(out) :: sigma      ! density/sea-level standard density
    real(kind_real), intent(out) :: delta      ! pressure/sea-level standard pressure
    real(kind_real), intent(out) :: theta      ! temperature/sea-level standard temperature
    ! local constants
    real(kind_real)  , parameter :: rearth = 6369.0     ! radius of the earth (km)
    real(kind_real)  , parameter :: gmr    = 34.163195  ! hydrostatic constant
    integer(kind_int), parameter :: ntab   = 8          ! number of entries in the defining tables
    ! local variables
    integer(kind_int) :: i,j,k             ! counters
    real(kind_real)   :: h                 ! geopotential altitude (km)
    real(kind_real)   :: tgrad, tbase      ! temperature gradient and base temp of this layer
    real(kind_real)   :: tlocal            ! local temperature
    real(kind_real)   :: deltah            ! height above base of this layer
    ! local arrays(1976 std. atmosphere)
    real(kind_real), parameter :: htab(ntab) = &
        (/0.0, 11.0, 20.0, 32.0, 47.0, 51.0, 71.0, 84.852/)
    real(kind_real), parameter :: ttab(ntab) = &
        (/288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65, 186.946/)
    real(kind_real), parameter :: ptab(ntab) = &
        (/1.0, 2.233611e-1, 5.403295e-2, 8.5666784e-3, 1.0945601e-3, &
                            6.6063531e-4, 3.9046834e-5, 3.68501e-6/)
    real(kind_real), parameter :: gtab(ntab) = &
        (/-6.5, 0.0, 1.0, 2.8, 0.0, -2.8, -2.0, 0.0/)

    h = alt*rearth/(alt+rearth)      ! convert geometric to geopotential altitude

    i = 1
    j = ntab                         ! setting up for binary search
    do
        k = (i+j)/2                  ! integer division
        if (h < htab(k)) then
            j=k
        else
            i=k
        end if
        if (j <= i+1) exit
    end do

    tgrad  = gtab(i)                 ! i will be in 1...ntab-1
    tbase  = ttab(i)
    deltah = h-htab(i)
    tlocal = tbase+tgrad*deltah
    theta  = tlocal/ttab(1)          ! temperature ratio

    if (tgrad == 0.0) then           ! pressure ratio
        delta = ptab(i)*exp(-gmr*deltah/tbase)
    else
        delta = ptab(i)*(tbase/tlocal)**(gmr/tgrad)
    end if

    sigma = delta/theta              ! density ratio

    sigma = sigma*1.225
    delta = delta*101325.0
    theta = theta*288.15

end subroutine stdatm1976