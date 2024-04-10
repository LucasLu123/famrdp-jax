
!-----------------------------------------------------------------------------
! convert primitive variable to conserved variable
!-----------------------------------------------------------------------------
subroutine prim2con(nsp,nep,prim,nsc,nec,con)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,half,one
    use mod_variables, only : gamma
    implicit none
    integer(kind_int), intent(in)  :: nsp,nep
    real(kind_real)  , intent(in)  :: prim(nsp:nep)
    integer(kind_int), intent(in)  :: nsc,nec
    real(kind_real)  , intent(out) :: con(nsc:nec)
    real(kind_real)   :: ro,v,v2,ae,e,gama
    integer(kind_int) :: i

    gama = gamma

    ro = prim(nsp)
    con(nsc) = ro

    v2 = zero
    do i=1,3
       v = prim(nsp+i)
       con(nsc+i) = ro*v
       v2 = v2 + v*v
    end do

    ae = gama - one
    e  = prim(nsp+4)/ae + half*ro*v2
    con(nsc+4) = e

end subroutine prim2con
!-----------------------------------------------------------------------------
! convert conserved variable to primitive variable
!-----------------------------------------------------------------------------
subroutine con2prim(nsc,nec,con,nsp,nep,prim)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,half,one
    use mod_variables, only : gamma
    implicit none
    integer(kind_int), intent(in)  :: nsc,nec
    real(kind_real)  , intent(in)  :: con(nsc:nec)
    integer(kind_int), intent(in)  :: nsp,nep
    real(kind_real)  , intent(out) :: prim(nsp:nep)
    real(kind_real)   :: ro,v,v2,ae,p,gama
    integer(kind_int) :: i

    gama = gamma

    ro = con(nsc)
    prim(nsp) = ro

    v2 = zero
    do i=1,3
        v = con(nsc+i)/ro
        prim(nsp+i) = v
        v2 = v2 + v*v
    end do

    ae = gama - one
    p  = ae*(con(nsc+4) - half*ro*v2)
    prim(nsp+4) = p

end subroutine con2prim


subroutine prim2tc(nsp,nep,prim,nst,net,tc)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_variables, only : gamma,refbeta
    implicit none
    integer(kind_int), intent(in)  :: nsp,nep
    real(kind_real)  , intent(in)  :: prim(nsp:nep)
    integer(kind_int), intent(in)  :: nst,net
    real(kind_real)  , intent(out) :: tc(nst:net)
    real(kind_real)   :: ro,p,gama,t,c

    gama = gamma

    ro = prim(nsp)
    p  = prim(nsp+1)
    t  = p/(refbeta*ro)
    c  = sqrt(gama*refbeta*t)
    tc(nst)   = t
    tc(nst+1) = c

end subroutine prim2tc


subroutine prim2tc_SCMP(nsp,nep,prim,nst,net,tc)    !> Done
    use mod_kndconsts, only : kind_int,kind_real
    use mod_variables, only : gamma,refbeta
    use mod_variables, only : moo,poo
    use mod_constants, only : SCMP_sigma    !> Done
    implicit none
    integer(kind_int), intent(in)  :: nsp,nep
    real(kind_real)  , intent(in)  :: prim(nsp:nep)
    integer(kind_int), intent(in)  :: nst,net
    real(kind_real)  , intent(out) :: tc(nst:net)
    real(kind_real)   :: ro,p,gama,t,c,pPrime
    
    gama = gamma

    p  = prim(nsp+1)
    pPrime = p - poo
    ro = (1.0 + gama*moo*moo*pPrime)**(1.0/SCMP_sigma)            !> Done
    t  = p/(refbeta*ro)
    c = sqrt( SCMP_sigma / (gama*moo*moo*ro**(1-SCMP_sigma)) )    !> Done
    tc(nst)   = t
    tc(nst+1) = c

end subroutine prim2tc_SCMP
    
subroutine consQ2consT(nsc,nec,consQ,nsp,nep,consT)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_variables, only : gamma,poo
    use mod_constants, only : zero,one,half
    implicit none
    integer(kind_int), intent(in)  :: nsc,nec
    real(kind_real)  , intent(in)  :: consQ(nsc:nec)
    integer(kind_int), intent(in)  :: nsp,nep
    real(kind_real)  , intent(out) :: consT(nsp:nep)
    real(kind_real)   :: p,pPrime
    real(kind_real)   :: ro,v,v2,ae,gama
    integer(kind_int) :: i

    gama = gamma

    ro = consQ(nsc) !> consQ = {��,��u,��v,��w,E}

    v2 = zero
    do i=1,3
        v = consQ(nsc+i)/ro
        v2 = v2 + v*v
    end do

    ae = gama - one
    p  = ae*(consQ(nsc+4) - half*ro*v2)
    pPrime = p - poo
    
    consT(nsp  ) = pPrime
    consT(nsp+4) = pPrime
    do i=1,3
        consT(nsp+i) = consQ(nsc+i) !> Done: consT = {p',��u,��v,��w,p'}
    end do

end subroutine consQ2consT
    
subroutine primQ2primT(nsc,nec,primQ,nsp,nep,primT)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_variables, only : poo
    implicit none
    integer(kind_int), intent(in)  :: nsc,nec
    real(kind_real)  , intent(in)  :: primQ(nsc:nec)
    integer(kind_int), intent(in)  :: nsp,nep
    real(kind_real)  , intent(out) :: primT(nsp:nep)
    real(kind_real)   :: p,pPrime
    integer(kind_int) :: i

    p  = primQ(nsc+4)   !> primQ = {��,u,v,w,p}
    pPrime = p - poo
    
    !> Done: primT = {��,u,v,w,p'}
    do i=0,3
        primT(nsp+i) = primQ(nsc+i)
    end do
    primT(nsp+4) = pPrime

end subroutine primQ2primT
    

subroutine calculatePseudoVelocityAndAcousticVelocity(prim,kx,ky,kz,viscosity,length,alpha,Ur,pseudo_c,pseudo_v)
	use mod_kndconsts, only : kind_int,kind_real
    use mod_variables, only : gamma,refbeta,nvis,reue,moo,uoo,voo,woo
    use mod_constants, only : nvis_euler
    implicit none
    real(kind_real)  , intent(in)  :: prim(1:5),kx,ky,kz,viscosity,length
    real(kind_real)  , intent(out)  :: pseudo_c,pseudo_v,alpha,Ur
    integer(kind_int) :: i,j
    real(kind_real)   :: ro,p,gama,T,c,vx,vy,vz,v,cp,H,uuinf
    real(kind_real)   :: drodp,drodT,theta,phi,beta
    real(kind_real)   :: contra_velocity, snk, ep
    real(kind_real), parameter ::epsilon = 1.0E-5

	gama = gamma

    ro = prim(1)
    p  = prim(5)
    vx = prim(2)
    vy = prim(3)
    vz = prim(4)
    v  = sqrt(vx*vx + vy*vy + vz*vz)
    
    T  = p/(refbeta*ro)
    cp = gama/(gama - 1)*refbeta
    H  = gama/(gama - 1)*p/ro + 0.5*(vx*vx + vy*vy + vz*vz)
    c  = sqrt(gama*refbeta*t)
    drodp = ro/p
    drodT = -ro/T
    
    
    !����Ur ws-95 ��ʽ(7) �����������
    if (v < epsilon*c) then
        Ur = epsilon*c
    else if (v > c) then
        Ur = c
    else
        Ur = v
    endif
    
    !ճ��������Ҫ��Urʩ������
    if (nvis > nvis_euler) then
        Ur = max(Ur, viscosity/(ro*length)/reue)
    end if
    
    !liaofei 2023 JCP
    !uuinf = sqrt(uoo*uoo + voo*voo + woo*woo)
    !Ur = min(max(v,2.*uuinf),c)
    !if (nvis > nvis_euler) then
    !    Ur = min(max(v,2.*uuinf,viscosity/(ro*length)/reue),c)
    !end if
    
    !����theta�������л�����drodp��ʾ ws-95 ��ʽ(6)
    theta = 1.0/(Ur*Ur) -drodT/(ro*cp)
    phi = drodT
   
    beta = drodp + drodT/(ro*cp)
    alpha = (1. - beta*Ur*Ur)*0.5
    snk = sqrt(kx**2 + ky**2 + kz**2)
    contra_velocity = (vx*kx + vy*ky + vz*kz)/snk
    pseudo_v = contra_velocity*(1.-alpha)
    pseudo_c = sqrt(alpha**2*contra_velocity**2+Ur**2)
    !pseudo_c = sqrt(alpha**2*contra_velocity**2+Ur**2)*snk

end subroutine calculatePseudoVelocityAndAcousticVelocity

    
    
subroutine calculate_precondition_eigenvalue_point(nsp,nep,sn,prim,viscosity,length,eigenvalue)
use mod_constants, only : nvis_euler,sml_ssf
    use mod_datatypes, only : fld_array_t
    use mod_kndconsts, only : kind_int,kind_real
    use mod_variables, only : gamma,refbeta,nvis,reue,moo,uoo,voo,woo
    implicit none
    integer(kind_int), intent(in)  :: nsp,nep
    real(kind_real)  , intent(in)  :: viscosity,length
    real(kind_real)  , intent(in)  :: prim(nsp:nep),sn(1:3,1:3)
    real(kind_real)  , intent(out) :: eigenvalue(1:3,nsp:nep)
    integer(kind_int) :: i,j
    real(kind_real)   :: ro,p,gama,T,c,vx,vy,vz,v,cp,H,uuinf
    real(kind_real)   :: drodp,drodT,theta,phi,beta,Ur,alpha
    real(kind_real)   :: contra_velocity, pseudo_c, snk, pseudo_v, ep
    real(kind_real)   :: machsqr,moosqr,vSqr
    real(kind_real)   :: roCpThetaPlusPhi,oroCpThetaPlusPhi,oro
    real(kind_real), parameter ::epsilon = 1.0E-5

    gama = gamma

    ro = prim(nsp)
    vx = prim(nsp+1)
    vy = prim(nsp+2)
    vz = prim(nsp+3)
    p  = prim(nsp+4)
    v  = sqrt(vx*vx + vy*vy + vz*vz)
    
    T  = p/(refbeta*ro)
    cp = gama/(gama - 1)*refbeta
    H  = gama/(gama - 1)*p/ro + 0.5*(vx*vx + vy*vy + vz*vz)
    c  = sqrt(gama*refbeta*t)
    drodp = ro/p
    drodT = -ro/T
    
    !����Ԥ��������
    !machsqr = (v/c)**2
    !moosqr = moo*moo
    !prec_coeff = min(max(machsqr,moosqr),1.0)
    
    !beta = prec_coeff
    
    !����Ur ws-95 ��ʽ(7) �����������
    if (v < epsilon*c) then
        Ur = epsilon*c
    else if (v > c) then
        Ur = c
    else
        Ur = v
    endif
    
    !ճ��������Ҫ��Urʩ������
    if (nvis > nvis_euler) then
        Ur = max(Ur, viscosity/(ro*length)/reue)
    end if
    
    !liaofei 2023 JCP
    !uuinf = sqrt(uoo*uoo + voo*voo + woo*woo)
    !Ur = min(max(v,2.*uuinf),c)
    !if (nvis > nvis_euler) then
    !    Ur = min(max(v,2.*uuinf,viscosity/(ro*length)/reue),c)
    !end if
    
    !����theta�������л�����drodp��ʾ ws-95 ��ʽ(6)
    theta = 1.0/(Ur*Ur) -drodT/(ro*cp)
    
    !Ф���� ��ʿ���� ��������
    !ep = beta/(1. + (gama - 1.)*beta)
    !theta = 1.0/(ep*c*c)
    phi = drodT
    
    beta = drodp + drodT/(ro*cp)
    alpha = (1. - beta*Ur*Ur)*0.5
    do i = 1,3
        contra_velocity = vx*sn(i,1) + vy*sn(i,2) + vz*sn(i,3)  
        snk = sqrt(sn(i,1)**2 + sn(i,2)**2 + sn(i,3)**2)
        eigenvalue(i,1:3) = contra_velocity
        contra_velocity = contra_velocity/max(snk,sml_ssf)
        pseudo_v = contra_velocity*(1.-alpha)
        pseudo_c = sqrt(alpha**2*contra_velocity**2+Ur**2)
        eigenValue(i,4) = (abs(pseudo_v) + pseudo_c)*snk
        eigenValue(i,5) = (abs(pseudo_v) - pseudo_c)*snk
    end do


end subroutine calculate_precondition_eigenvalue_point
    
    
!����Ԥ��������todo
!viscosity ճ��ϵ��
!length ��������߶�
!prim ԭʼ��������
!preMatrix Ԥ�������󣬴��غ������ԭʼ�����ĵ��������ݻ�����
!subroutine calculate_precondition_matrix(nsp,nep,kx,ky,kz,prec_coeff,viscosity,length,prim,eigenvalue,nonpreMatrix,preMatrix)
subroutine calculate_precondition_matrix_point(nsp,nep,prim,nonpreMatrix)
    use mod_constants, only : nvis_euler,sml_ssf
    use mod_datatypes, only : fld_array_t
    use mod_kndconsts, only : kind_int,kind_real
    use mod_variables, only : gamma,refbeta,nvis,reue,moo,uoo,voo,woo
    implicit none
    integer(kind_int), intent(in)  :: nsp,nep
    real(kind_real)  , intent(in)  :: prim(nsp:nep)
    real(kind_real)  , intent(out) :: nonpreMatrix(1:4)
    integer(kind_int) :: i,j
    real(kind_real)   :: ro,p,gama,T,c,vx,vy,vz,v,cp,H
    real(kind_real)   :: drodp,drodT,theta,phi,beta,Ur,alpha

    gama = gamma

    ro = prim(nsp)
    vx = prim(nsp+1)
    vy = prim(nsp+2)
    vz = prim(nsp+3)
    p  = prim(nsp+4)
    v  = sqrt(vx*vx + vy*vy + vz*vz)
    
    T  = p/(refbeta*ro)
    cp = gama/(gama - 1)*refbeta
    H  = gama/(gama - 1)*p/ro + 0.5*(vx*vx + vy*vy + vz*vz)
    c  = sqrt(gama*refbeta*t)
    drodp = ro/p
    drodT = -ro/T
    
    !ֻ������ά�������ά�����ѭ������Ҫ�ı�todo
    nonpreMatrix(:) = 0.
    
    !�����غ������ԭʼ�����ĵ�������
    !�غ���� rho, rhou, rhov, rhow, rhoE
    !ԭʼ���� p, u, v, w, T
    nonpreMatrix(1) = drodp
    nonpreMatrix(2) = drodT
    nonpreMatrix(3) = drodp*H - 1.
    nonpreMatrix(4) = ro*cp+H*drodT
    
end subroutine calculate_precondition_matrix_point
    
subroutine primitiveDeltaWToConservativeDeltaQ(point_nonprec_i,dqi,pvi,dUi)
    use mod_datatypes, only : kind_real,kind_int
    implicit none
    real(kind_real),intent(in)  :: point_nonprec_i(1:4),dqi(1:5),pvi(1:5)
    real(kind_real),intent(out) :: dUi(1:5)
    
    dUi(1) = point_nonprec_i(1)*dqi(1) + point_nonprec_i(2)*dqi(5)
    dUi(2) = pvi(2)*dUi(1) + pvi(1)*dqi(2)
    dUi(3) = pvi(3)*dUi(1) + pvi(1)*dqi(3)
    dUi(4) = pvi(4)*dUi(1) + pvi(1)*dqi(4)
    dUi(5) = point_nonprec_i(3)*dqi(1) + point_nonprec_i(4)*dqi(5)
    dUi(5) = dUi(5)+ pvi(1)*(pvi(2)*dqi(2) + pvi(3)*dqi(3) + pvi(4)*dqi(4))

end subroutine primitiveDeltaWToConservativeDeltaQ
    
subroutine calculateDiagMatrixInv(prim,dt,lambda,viscosity,length,diagInvMatrix)
    use mod_fieldvars, only : npvs
    use mod_datatypes, only : kind_real,kind_int
    use mod_variables, only : gamma,refbeta,nvis,reue,moo,vinf,uoo,voo,woo
    use mod_constants, only : nvis_euler
    implicit none
    real(kind_real),intent(in)    :: prim(1:npvs),dt,lambda,viscosity,length
    real(kind_real),intent(out)   :: diagInvMatrix(1:25)
    integer(kind_int) :: i,j
    real(kind_real)   :: ro,p,gama,T,c,vx,vy,vz,v,cp,H,uuinf
    real(kind_real)   :: drodp,drodT,theta,bbb,Ur
    !real(kind_real)   :: contra_velocity, pseudo_c, snk, pseudo_v, ep
    real(kind_real)   :: machsqr,moosqr,ep,prec_coeff
    real(kind_real)   :: alpha,beta,phi,psi,tau
    real(kind_real)   :: lineOneSame,lineTwoSame,lineThreeSame,lineFourSame,lineFiveSame,diagSame
    real(kind_real), parameter ::epsilon = 1.0E-5

    gama = gamma

    ro = prim(1)
    vx = prim(2)
    vy = prim(3)
    vz = prim(4)
    p  = prim(5)
    v  = sqrt(vx*vx + vy*vy + vz*vz)
    
    T  = p/(refbeta*ro)
    cp = gama/(gama - 1)*refbeta
    H  = gama/(gama - 1)*p/ro + 0.5*(vx*vx + vy*vy + vz*vz)
    c  = sqrt(gama*refbeta*t)
    drodp = ro/p
    drodT = -ro/T
    
    !����Ԥ��������
    machsqr = (v/c)**2
    moosqr = moo*moo
    prec_coeff = min(max(machsqr,moosqr),1.0)
    
    bbb = prec_coeff
    
    !����Ur ws-95 ��ʽ(7) �����������
    if (v < epsilon*c) then
        Ur = epsilon*c
    else if (v > c) then
        Ur = c
    else
        Ur = v
    endif
    
    !ճ��������Ҫ��Urʩ������
    if (nvis > nvis_euler) then
        Ur = max(Ur, viscosity/(ro*length)/reue)
    end if
    
    !liaofei 2023 JCP
    !uuinf = sqrt(uoo*uoo + voo*voo + woo*woo)
    !Ur = min(max(v,2.*uuinf),c)
    !if (nvis > nvis_euler) then
    !    Ur = min(max(v,2.*uuinf,viscosity/(ro*length)/reue),c)
    !end if
    
    !����theta�������л�����drodp��ʾ ws-95 ��ʽ(6)
    theta = 1.0/(Ur*Ur) -drodT/(ro*cp)
    
    !Ф���� ��ʿ���� ��������
    !ep = bbb/(1. + (gama - 1.)*bbb)
    !theta = 1.0/(ep*c*c)
    
    diagInvMatrix(:) = 0.

    tau = dt
    !�ؼ��͵ļ򻯹�ʽ
    alpha = lambda*tau+1
    beta = theta + lambda*tau*drodp
    phi = ro*cp*beta + alpha*drodT
    psi = vx*vx+vy*vy+vz*vz-H
    
    lineOneSame = tau*drodT/phi
    diagInvMatrix(1) = tau*(ro*cp-psi*drodT)/phi
    diagInvMatrix(2) = vx*lineOneSame
    diagInvMatrix(3) = vy*lineOneSame
    diagInvMatrix(4) = vz*lineOneSame
    diagInvMatrix(5) = -lineOneSame
    
    lineTwoSame = tau/ro/alpha
    diagInvMatrix(6) = -vx*lineTwoSame
    diagInvMatrix(7) = lineTwoSame
    
    lineThreeSame = lineTwoSame
    diagInvMatrix(11) = -vy*lineThreeSame
    diagInvMatrix(13) = lineThreeSame
    
    lineFourSame = lineTwoSame
    diagInvMatrix(16) = -vz*lineFourSame
    diagInvMatrix(19) = lineFourSame
    
    lineFiveSame = tau/alpha/phi
    diagInvMatrix(21) = (alpha+beta*psi)*lineFiveSame
    diagInvMatrix(22) = -vx*beta*lineFiveSame
    diagInvMatrix(23) = -vy*beta*lineFiveSame
    diagInvMatrix(24) = -vz*beta*lineFiveSame
    diagInvMatrix(25) = beta*lineFiveSame
    
end subroutine calculateDiagMatrixInv    
    
    
subroutine v1_eq_v2(ns1,ne1,v1,ns2,ne2,v2)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,half,one
    use mod_variables, only : gamma
    implicit none
    integer(kind_int), intent(in)  :: ns1,ne1
    real(kind_real)  , intent(in)  :: v1(ns1:ne1)
    integer(kind_int), intent(in)  :: ns2,ne2
    real(kind_real)  , intent(out) :: v2(ns2:ne2)
    integer(kind_int) :: i,imax

    do i=0,ne1-ns1
       v2(ns2+i) = v1(ns1+i)
    end do

end subroutine v1_eq_v2

!< todo ��ʱ�����prim����Ϊro,u,v,w,p'
!< todo SCM-P�Ķ���ͨ���������
subroutine flux_euler(nsp,nep,prim,nt,nx,ny,nz,nsf,nef,f)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : sml_ssf,one,two,half
    use mod_variables, only : gamma
    implicit none
    integer(kind_int), intent(in)  :: nsp,nep
    real(kind_real)  , intent(in)  :: prim(nsp:nep)
    real(kind_real)  , intent(in)  :: nt,nx,ny,nz
    integer(kind_int), intent(in)  :: nsf,nef
    real(kind_real)  , intent(out) :: f(nsf:nef)
    real(kind_real) :: ro,vx,vy,vz,ps
    real(kind_real) :: gama,ae,v2,h0,vn,rvn

    ro = prim(nsp)
    vx = prim(nsp+1)
    vy = prim(nsp+2)
    vz = prim(nsp+3)
    ps = prim(nsp+4)

    gama = gamma
    ae   = gama - one

    v2 = vx*vx + vy*vy + vz*vz
    h0 = gama*ps/(ro*ae) + half*v2

    vn = nx*vx + ny*vy + nz*vz

    rvn = ro*vn
    f(nsf  ) = rvn
    f(nsf+1) = rvn*vx + ps*nx
    f(nsf+2) = rvn*vy + ps*ny
    f(nsf+3) = rvn*vz + ps*nz
    f(nsf+4) = rvn*h0

end subroutine flux_euler

subroutine flux_steger(nsp,nep,priml,primr,nt,nx,ny,nz,nsf,nef,flr,efix)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : one
    implicit none
    integer(kind_int), intent(in)  :: nsp,nep
    real(kind_real)  , intent(in)  :: priml(nsp:nep),primr(nsp:nep)
    real(kind_real)  , intent(in)  :: nt,nx,ny,nz
    integer(kind_int), intent(in)  :: nsf,nef
    real(kind_real)  , intent(out) :: flr(nsf:nef)
    real(kind_real)  , intent(in)  :: efix
    integer(kind_int) :: m
    real(kind_real)   :: fl(nsf:nef),fr(nsf:nef)

    call flux_sw_pn(nsp,nep,priml,nt,nx,ny,nz,nsf,nef,fl,efix, one)
    call flux_sw_pn(nsp,nep,primr,nt,nx,ny,nz,nsf,nef,fr,efix,-one)

    do m=nsf,nef
        flr(m) = fl(m)+fr(m)
    end do

end subroutine flux_steger


subroutine flux_sw_pn(nsp,nep,prim,nt,nx,ny,nz,nsf,nef,f,efix,fsw)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : sml_ssf,one,two,half
    use mod_variables, only : gamma
    implicit none
    integer(kind_int), intent(in)  :: nsp,nep
    real(kind_real)  , intent(in)  :: prim(nsp:nep)
    real(kind_real)  , intent(in)  :: nt,nx,ny,nz
    integer(kind_int), intent(in)  :: nsf,nef
    real(kind_real)  , intent(out) :: f(nsf:nef)
    real(kind_real)  , intent(in)  :: efix,fsw
    real(kind_real)           :: ro,vx,vy,vz,ps,gama,ae
    real(kind_real)           :: c,c2,e,e0,h,h0,v2,vn,eps
    real(kind_real)           :: sn,osn,nxa,nya,nza,vna
    real(kind_real)           :: l1,l4,l5,x1,x2,u2,c2r
    real(kind_real)           :: al1,al4,al5,ma
    real(kind_real), external :: enfix_steger,enfix_harten
    integer(kind_int)         :: nfix=0

    ro = prim(nsp)
    vx = prim(nsp+1)
    vy = prim(nsp+2)
    vz = prim(nsp+3)
    ps = prim(nsp+4)

    gama = gamma
    ae   = gama - one

    v2 = vx*vx + vy*vy + vz*vz
    c2 = gama*ps/ro
    h  = c2/ae
    e  = h - ps/ro
    h0 = h + half*v2
    e0 = e + half*v2

    c  = sqrt(c2)
    vn = nx*vx + ny*vy + nz*vz

    sn = max(sqrt(nx*nx + ny*ny + nz*nz),sml_ssf)
    osn = one/sn
    nxa = nx * osn
    nya = ny * osn
    nza = nz * osn
    vna = vn * osn

    l1 = vn
    l4 = vn + sn*c
    l5 = vn - sn*c

    ma = sqrt(v2)/c
    eps = efix*sn*c*(one + ma)

    if (nfix == 0) then
        al1 = enfix_steger(l1,eps)
        al4 = enfix_steger(l4,eps)
        al5 = enfix_steger(l5,eps)
    else
        al1 = enfix_harten(l1,eps)
        al4 = enfix_harten(l4,eps)
        al5 = enfix_harten(l5,eps)
    end if

    l1 = half*(l1 + fsw*al1)
    l4 = half*(l4 + fsw*al4)
    l5 = half*(l5 + fsw*al5)

    c2r = c2 / gama
    x1 = c2r * ( two*l1 - l4 - l5 )/( two * c2 )
    x2 = c2r * ( l4 - l5 )/( two * c )

    f(nsf  ) = (l1    - x1             ) * ro
    f(nsf+1) = (l1*vx - x1*vx + nxa*x2 ) * ro
    f(nsf+2) = (l1*vy - x1*vy + nya*x2 ) * ro
    f(nsf+3) = (l1*vz - x1*vz + nza*x2 ) * ro
    f(nsf+4) = (l1*e0 - x1*h0 + vna*x2 ) * ro

end subroutine flux_sw_pn

function enfix_steger(l0,eps) result(l1)
    use mod_kndconsts, only : kind_real
    implicit none
    real(kind_real) :: l0,eps
    real(kind_real) :: l1

    l1 = sqrt(l0*l0 + eps*eps)

end function enfix_steger

function enfix_harten(l0,eps) result(l1)
    use mod_kndconsts, only : kind_real
    use mod_constants, only : half
    implicit none
    real(kind_real) :: l0,eps
    real(kind_real) :: l1
    real(kind_real) :: al0

    al0 = abs(l0)
    if (al0 < eps) then
        l1 = half*(l0*l0 + eps*eps)/eps
    else
        l1 = al0
    end if

end function enfix_harten

subroutine flux_sw_mod(nsp,nep,priml,primr,nt,nx,ny,nz,nsf,nef,flr,efix)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : half,one,two,id_ps
    implicit none
    integer(kind_int), intent(in)  :: nsp,nep
    real(kind_real)  , intent(in)  :: priml(nsp:nep),primr(nsp:nep)
    real(kind_real)  , intent(in)  :: nt,nx,ny,nz
    integer(kind_int), intent(in)  :: nsf,nef
    real(kind_real)  , intent(out) :: flr(nsf:nef)
    real(kind_real)  , intent(in)  :: efix
    integer(kind_int)          :: m
    real(kind_real)            :: pvl(nsp:nep),pvr(nsp:nep)
    real(kind_real)            :: ql(nsf:nef),qr(nsf:nef)
    real(kind_real)            :: fl(nsf:nef),fr(nsf:nef)
    real(kind_real)            :: pl,pr,gp,wp,rwp,we,rwe,efix0
    real(kind_real), parameter :: sig2=0.5

    call prim2con(nsp,nep,priml,nsf,nef,ql)
    call prim2con(nsp,nep,primr,nsf,nef,qr)

    pl = priml(id_ps)
    pr = primr(id_ps)

    gp = sig2*abs(pl - pr)/min(pl, pr)   !!!������,�д��Ľ�
    wp = half/(gp*gp + one)
    rwp = one - wp

    we = two*(half - wp)
    rwe = one - we
    !!efix0 = efix
    efix0 = we*efix + rwe*0.05

    pvl(:) = rwp*priml(:) + wp*primr(:)
    pvr(:) = rwp*primr(:) + wp*priml(:)

    call mxdq_sw(nsp,nep,pvl,nt,nx,ny,nz,nsf,nef,ql,fl,efix0, one)
    call mxdq_sw(nsp,nep,pvr,nt,nx,ny,nz,nsf,nef,qr,fr,efix0,-one)

    do m=nsf,nef
        flr(m) = fl(m)+fr(m)
    end do

end subroutine flux_sw_mod

subroutine flux_vanleer(nsp,nep,priml,primr,nt,nx,ny,nz,nsf,nef,flr,efix)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,one,two,half,fourth,sml_ssf
    use mod_variables, only : gamma
    implicit none
    integer(kind_int), intent(in)  :: nsp,nep
    real(kind_real)  , intent(in)  :: priml(nsp:nep),primr(nsp:nep)
    real(kind_real)  , intent(in)  :: nt,nx,ny,nz
    integer(kind_int), intent(in)  :: nsf,nef
    real(kind_real)  , intent(out) :: flr(nsf:nef)
    real(kind_real)  , intent(in)  :: efix
    real(kind_real) :: gama,ae,gmp1,sn,osn,plr
    real(kind_real) :: rl,ul,vl,wl,pl,v2l,c2l,cl,vnl
    real(kind_real) :: rr,ur,vr,wr,pr,v2r,c2r,cr,vnr
    real(kind_real) :: hl,ml,mla,mpl,ppl,fhl,mrcl
    real(kind_real) :: hr,mr,mra,mmr,pmr,fhr,mrcr

    gama = gamma
    ae   = gama - one
    gmp1 = gama + one

    rl = priml(nsp)
    ul = priml(nsp+1)
    vl = priml(nsp+2)
    wl = priml(nsp+3)
    pl = priml(nsp+4)

    rr = primr(nsp)
    ur = primr(nsp+1)
    vr = primr(nsp+2)
    wr = primr(nsp+3)
    pr = primr(nsp+4)

    v2l = ul*ul + vl*vl + wl*wl
    c2l = gama*pl/rl
    hl  = c2l/ae + half*v2l

    v2r = ur*ur + vr*vr + wr*wr
    c2r = gama*pr/rr
    hr  = c2r/ae + half*v2r

    sn = max(sqrt(nx*nx + ny*ny + nz*nz),sml_ssf)
    osn = one/sn

    cl  = sqrt(c2l)
    vnl = (nx*ul + ny*vl + nz*wl)*osn
    cr  = sqrt(c2r)
    vnr = (nx*ur + ny*vr + nz*wr)*osn

    ml  = vnl/cl
    mla = abs(ml)
    mr  = vnr/cr
    mra = abs(mr)

    if (mla >= one) then
        mpl = half*(ml + mla)
        ppl = half*(one + sign(one,ml))
        fhl = zero
    else
        mpl = fourth*(ml + one)**2
        ppl = fourth*(two - ml)*(ml + one)**2
        fhl = one
    end if

    if (mra >= one) then
        mmr = half*(mr - mra)
        pmr = half*(one - sign(one,mr))
        fhr = zero
    else
        mmr = -fourth*(mr - one)**2
        pmr =  fourth*(two + mr)*(mr - one)**2
        fhr = one
    end if

    ! original hr/hl,not ensure the exact constancy of
    ! the total enthaly for steady inviscid flows
    hl = hl - fhl*(vnl - cl)**2/gmp1
    hr = hr - fhr*(vnr + cr)**2/gmp1

    plr = ppl*pl + pmr*pr

    mrcl = mpl*rl*cl*sn
    mrcr = mmr*rr*cr*sn

    flr(nsf  ) = mrcl    + mrcr
    flr(nsf+1) = mrcl*ul + mrcr*ur + plr*nx
    flr(nsf+2) = mrcl*vl + mrcr*vr + plr*ny
    flr(nsf+3) = mrcl*wl + mrcr*wr + plr*nz
    flr(nsf+4) = mrcl*hl + mrcr*hr

end subroutine flux_vanleer

subroutine flux_roe(nsp,nep,priml,primr,nt,nx,ny,nz,nsf,nef,flr,efix)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,one,two,half,fourth,sml_ssf
    use mod_variables, only : gamma
    implicit none
    integer(kind_int), intent(in)  :: nsp,nep
    real(kind_real)  , intent(in)  :: priml(nsp:nep),primr(nsp:nep)
    real(kind_real)  , intent(in)  :: nt,nx,ny,nz
    integer(kind_int), intent(in)  :: nsf,nef
    real(kind_real)  , intent(out) :: flr(nsf:nef)
    real(kind_real)  , intent(in)  :: efix
    real(kind_real)           :: gama,ae,sn,osn,csn,plr,eps
    real(kind_real)           :: rl,ul,vl,wl,pl,v2l,c2l,cl,vnl
    real(kind_real)           :: rr,ur,vr,wr,pr,v2r,c2r,cr,vnr
    real(kind_real)           :: hl,rls,hr,rrs,oors,qi(5),adq(5)
    real(kind_real)           :: vni,vi2,vi2p2,ci2,ci,nxa,nya,nza
    real(kind_real)           :: dvn,dro,dvx,dvy,dvz,dps,dh,rvnl,rvnr
    real(kind_real)           :: l1,l4,l5,al1,al4,al5,rp1,rp4,rp5
    real(kind_real), external :: enfix_harten

    gama = gamma
    ae   = gama - one

    rl = priml(nsp)
    ul = priml(nsp+1)
    vl = priml(nsp+2)
    wl = priml(nsp+3)
    pl = priml(nsp+4)

    rr = primr(nsp)
    ur = primr(nsp+1)
    vr = primr(nsp+2)
    wr = primr(nsp+3)
    pr = primr(nsp+4)

    v2l = ul*ul + vl*vl + wl*wl
    c2l = gama*pl/rl
    hl  = c2l/ae + half*v2l

    v2r = ur*ur + vr*vr + wr*wr
    c2r = gama*pr/rr
    hr  = c2r/ae + half*v2r

    sn = max(sqrt(nx*nx + ny*ny + nz*nz),sml_ssf)
    osn = one/sn
    nxa = nx*osn
    nya = ny*osn
    nza = nz*osn

    cl  = sqrt(c2l)
    vnl = (nxa*ul + nya*vl + nza*wl)
    cr  = sqrt(c2r)
    vnr = (nxa*ur + nya*vr + nza*wr)

    rls = sqrt(rl)
    rrs = sqrt(rr)
    qi(1) = rls*rrs

    oors = one/(rls + rrs)
    rls = rls*oors
    rrs = rrs*oors
    qi(2) = rls*ul + rrs*ur
    qi(3) = rls*vl + rrs*vr
    qi(4) = rls*wl + rrs*wr
    qi(5) = rls*hl + rrs*hr

    vni = qi(2)*nxa + qi(3)*nya + qi(4)*nza
    vi2 = qi(2)**2  + qi(3)**2  + qi(4)**2
    vi2p2 = half*vi2
    ci2 = ae*(qi(5) - vi2p2)
    ci = sqrt(ci2)

    l1 = vni
    l4 = vni + ci
    l5 = vni - ci

    eps = efix*(abs(vni) + ci)

    al1 = enfix_harten(l1,eps)
    al4 = enfix_harten(l4,eps)
    al5 = enfix_harten(l5,eps)

    dvn = vnr - vnl
    dro = rr - rl
    dvx = ur - ul
    dvy = vr - vl
    dvz = wr - wl
    dps = pr - pl
    dh  = qi(2)*dvx + qi(3)*dvy + qi(4)*dvz

    rp1 = dro - dps/ci2
    adq(1) = al1*( rp1                              )
    adq(2) = al1*( rp1*qi(2) + qi(1)*(dvx - nxa*dvn) )
    adq(3) = al1*( rp1*qi(3) + qi(1)*(dvy - nya*dvn) )
    adq(4) = al1*( rp1*qi(4) + qi(1)*(dvz - nza*dvn) )
    adq(5) = al1*( rp1*vi2p2 + qi(1)*(dh  - vni*dvn) )

    rp4 = al4*(half*(dps + qi(1)*ci*dvn)/ci2)
    rp5 = al5*(half*(dps - qi(1)*ci*dvn)/ci2)
    adq(1) = adq(1) + rp4                    + rp5
    adq(2) = adq(2) + rp4*( qi(2) + nxa*ci ) + rp5*( qi(2) - nxa*ci )
    adq(3) = adq(3) + rp4*( qi(3) + nya*ci ) + rp5*( qi(3) - nya*ci )
    adq(4) = adq(4) + rp4*( qi(4) + nza*ci ) + rp5*( qi(4) - nza*ci )
    adq(5) = adq(5) + rp4*( qi(5) + vni*ci ) + rp5*( qi(5) - vni*ci )

    rvnl = rl*vnl
    rvnr = rr*vnr
    plr  = pl + pr

    csn = half*sn

    flr(nsf  ) = csn*( rvnl    + rvnr              - adq(1) )
    flr(nsf+1) = csn*( rvnl*ul + rvnr*ur + nxa*plr - adq(2) )
    flr(nsf+2) = csn*( rvnl*vl + rvnr*vr + nya*plr - adq(3) )
    flr(nsf+3) = csn*( rvnl*wl + rvnr*wr + nza*plr - adq(4) )
    flr(nsf+4) = csn*( rvnl*hl + rvnr*hr           - adq(5) )

end subroutine flux_roe
    
subroutine flux_roe_prec(nsp,nep,priml,primr,length,nt,nx,ny,nz,nsf,nef,flr,efix)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,one,two,half,fourth,sml_ssf
    use mod_variables, only : gamma,refbeta,tssnd
    implicit none
    integer(kind_int), intent(in)  :: nsp,nep
    real(kind_real)  , intent(in)  :: priml(nsp:nep),primr(nsp:nep)
    real(kind_real)  , intent(in)  :: nt,nx,ny,nz,length
    integer(kind_int), intent(in)  :: nsf,nef
    real(kind_real)  , intent(out) :: flr(nsf:nef)
    real(kind_real)  , intent(in)  :: efix
    real(kind_real)           :: gama,ae,csn,plr,eps
    real(kind_real)           :: rl,ul,vl,wl,pl,v2l,c2l,cl,vnl,El
    real(kind_real)           :: rr,ur,vr,wr,pr,v2r,c2r,cr,vnr,Er
    real(kind_real)           :: sn,osn,nxa,nya,nza
    real(kind_real)           :: Hl,rls,Hr,rrs,oors,qi(5),adq(5)
    real(kind_real)           :: prim(1:5),viSqr
    real(kind_real)           :: l1,l4,l5,al1,al4,al5,rp1,rp4,rp5
    real(kind_real)           :: pseudo_c,pseudo_v,alpha,Urr
    real(kind_real)           :: uu,du,deltaP,deltaU,temp1,temp2,cStar,MStar
    real(kind_real)           :: dissipation(nsf:nef),vis,temperature
    real(kind_real), external :: enfix_harten,sutherland

    gama = gamma
    ae   = gama - one

    rl = priml(nsp)
    ul = priml(nsp+1)
    vl = priml(nsp+2)
    wl = priml(nsp+3)
    pl = priml(nsp+4)

    rr = primr(nsp)
    ur = primr(nsp+1)
    vr = primr(nsp+2)
    wr = primr(nsp+3)
    pr = primr(nsp+4)

    v2l = ul*ul + vl*vl + wl*wl
    c2l = gama*pl/rl
    Hl  = c2l/ae + half*v2l
    El  = pl/rl/ae + half*v2l

    v2r = ur*ur + vr*vr + wr*wr
    c2r = gama*pr/rr
    Hr  = c2r/ae + half*v2r
    Er  = pr/rr/ae + half*v2r

    sn = max(sqrt(nx*nx + ny*ny + nz*nz),sml_ssf)
    osn = one/sn
    nxa = nx*osn
    nya = ny*osn
    nza = nz*osn

    cl  = sqrt(c2l)
    vnl = (nxa*ul + nya*vl + nza*wl)
    cr  = sqrt(c2r)
    vnr = (nxa*ur + nya*vr + nza*wr)
    
    
    flr(nsf)   = rl*vnl
    flr(nsf+1) = rl*ul*vnl + pl*nxa
    flr(nsf+2) = rl*vl*vnl + pl*nya
    flr(nsf+3) = rl*wl*vnl + pl*nza
    flr(nsf+4) = rl*Hl*vnl
    
    flr(nsf)   = flr(nsf) + rr*vnr
    flr(nsf+1) = flr(nsf+1) + rr*ur*vnr + pr*nxa
    flr(nsf+2) = flr(nsf+2) + rr*vr*vnr + pr*nya
    flr(nsf+3) = flr(nsf+3) + rr*wr*vnr + pr*nza
    flr(nsf+4) = flr(nsf+4) + rr*Hr*vnr

    rls = sqrt(rl)
    rrs = sqrt(rr)
    qi(1) = rls*rrs!Roeƽ���ܶ�

    oors = one/(rls + rrs)
    rls = rls*oors
    rrs = rrs*oors
    qi(2) = rls*ul + rrs*ur!<Roeƽ���ٶ�u
    qi(3) = rls*vl + rrs*vr!<Roeƽ���ٶ�v
    qi(4) = rls*wl + rrs*wr!<Roeƽ���ٶ�w
    qi(5) = rls*Hl + rrs*Hr!<Roeƽ������H
    
    prim(:) = qi(:)
    viSqr = qi(2)*qi(2)+qi(3)*qi(3)+qi(4)*qi(4)
    prim(5) = (qi(5)-0.5*viSqr)*ae/gama*qi(1)
    temperature = prim(5)/(qi(1)*refbeta)
    vis = sutherland(temperature,tssnd)

    !vni = qi(2)*nxa + qi(3)*nya + qi(4)*nza!<Roeƽ������ٶ�
    !vi2 = qi(2)**2  + qi(3)**2  + qi(4)**2 !<Roeƽ���ٶ�ƽ��
    !vi2p2 = half*vi2
    !ci2 = ae*(qi(5) - vi2p2)
    !ci = sqrt(ci2)                         !<Roeƽ������
    
    call calculatePseudoVelocityAndAcousticVelocity(prim,nx,ny,nz,vis,length,alpha,Urr,pseudo_c,pseudo_v)
    
    temp1 = abs(pseudo_v + pseudo_c)
	temp2 = abs(pseudo_v - pseudo_c)
	cStar = 0.5*(temp1 + temp2)
	MStar = 0.5*(temp1 - temp2)/pseudo_c

	du = (ur-ul)*nxa + (vr-vl)*nya + (wr-wl)*nza
	uu = qi(2)*nxa + qi(3)*nya + qi(4)*nza
	deltaP = MStar*(pr - pl) + (cStar - abs(uu) + alpha*uu*MStar)*qi(1)*du
	deltaU = MStar*du + (cStar - (1-2.*alpha)*abs(uu) - alpha*uu*MStar)*(pr - pl)/(qi(1)*Urr*Urr)
    
    dissipation(nsf)   = abs(uu)*(rr - rl) + deltaU*qi(1)
    dissipation(nsf+1) = abs(uu)*(rr*ur-rl*ul) + deltaU*qi(1)*qi(2) + deltaP*nxa
    dissipation(nsf+2) = abs(uu)*(rr*vr-rl*vl) + deltaU*qi(1)*qi(3) + deltaP*nya
    dissipation(nsf+3) = abs(uu)*(rr*wr-rl*wl) + deltaU*qi(1)*qi(4) + deltaP*nza
    dissipation(nsf+4) = abs(uu)*(rr*Er-rl*El) + deltaU*qi(1)*qi(5) + deltaP*uu
    
    flr(:) = 0.5*(flr(:) - dissipation(:))*sn
    !flr(:) = 0.5*flr(:)*sn

end subroutine flux_roe_prec
    
subroutine flux_roe_scmp(nsp,nep,priml,primr,nt,nx,ny,nz,nsf,nef,flr,efix)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,one,two,half,fourth,sml_ssf,scmp_sigma
    use mod_variables, only : gamma,refbeta,tssnd,moo
    implicit none
    integer(kind_int), intent(in)  :: nsp,nep
    real(kind_real)  , intent(in)  :: priml(nsp:nep),primr(nsp:nep)
    real(kind_real)  , intent(in)  :: nt,nx,ny,nz
    integer(kind_int), intent(in)  :: nsf,nef
    real(kind_real)  , intent(out) :: flr(nsf:nef)
    real(kind_real)  , intent(in)  :: efix
    integer(kind_int)              :: m
    real(kind_real)           :: gama,csn,plr,eps
    real(kind_real)           :: rl,ul,vl,wl,pl,cl,vnl
    real(kind_real)           :: rr,ur,vr,wr,pr,cr,vnr
    real(kind_real)           :: sn,osn,nxa,nya,nza
    real(kind_real)           :: rls,rrs,oors,qi(5)
    real(kind_real)           :: prim(1:5),viSqr
    real(kind_real)           :: l1,l4,l5,al1,al4,al5,rp1,rp4,rp5
    real(kind_real)           :: omega,delta,oasqr,alpha
    real(kind_real)           :: dtilde1,dtilde2,dtilde3,dtilde4,dtilde5,dtilde6
    real(kind_real)           :: lambda1,lambda2,lambda3,lambda4
    real(kind_real)           :: delta_rhou,delta_rhov,delta_rhow,delta_p
    real(kind_real)           :: dissipation(nsf:nef)
    real(kind_real), external :: enfix_harten,sutherland

    gama = gamma
    
    !< Ĭ��nt=0
    !< �ܶ�ͨ��״̬���̼���õ�
    !< �������ro,u,v,w,p' 
    ul = priml(2)
    vl = priml(3)
    wl = priml(4)
    pl = priml(5) 
    rl = (1.0 + gama*moo*moo*pl)**(1.0/scmp_sigma)

    ur = primr(2)
    vr = primr(3)
    wr = primr(4)
    pr = primr(5)
    rr = (1.0 + gama*moo*moo*pr)**(1.0/scmp_sigma)

    vnl = (nx*ul + ny*vl + nz*wl)
    vnr = (nx*ur + ny*vr + nz*wr)   
    
    !< �������ͨ�����Ĳ���
    flr(1) = rl*vnl
    flr(2) = rl*ul*vnl + pl*nx
    flr(3) = rl*vl*vnl + pl*ny
    flr(4) = rl*wl*vnl + pl*nz
    
    flr(1) = flr(1) + rr*vnr
    flr(2) = flr(2) + rr*ur*vnr + pr*nx
    flr(3) = flr(3) + rr*vr*vnr + pr*ny
    flr(4) = flr(4) + rr*wr*vnr + pr*nz

    !< ����Roeƽ������
    rls = sqrt(rl)
    rrs = sqrt(rr)
    qi(1) = rls*rrs         !<Roeƽ���ܶȦ�
    
    oors = one/(rls + rrs)
    rls = rls*oors
    rrs = rrs*oors
    qi(2) = rls*ul + rrs*ur !<Roeƽ���ٶ�u
    qi(3) = rls*vl + rrs*vr !<Roeƽ���ٶ�v
    qi(4) = rls*wl + rrs*wr !<Roeƽ���ٶ�w
    qi(5) = rls*pl + rrs*pr !<Roeƽ��ѹ��p'

    !< ��������ֵ theta = 0
    omega = qi(2)*nx + qi(3)*ny + qi(4)*nz
    sn = nx*nx + ny*ny + nz*nz
    oasqr = gama*moo*moo/scmp_sigma*qi(1)**(1.0 - scmp_sigma)
    delta = sqrt(omega*omega - oasqr*omega*omega + sn)
    lambda1 = abs(omega)
    lambda2 = abs(omega)
    lambda3 = abs(omega + delta)
    lambda4 = abs(omega - delta)
    
    !< ��������ֵ֮�ע���ʱ�ı���˳��Ϊp',u,v,w
    delta_rhou = rr*ur - rl*ul
    delta_rhov = rr*vr - rl*vl
    delta_rhow = rr*wr - rl*wl
    delta_p    = pr - pl
    alpha = delta_rhou*nx + delta_rhov*ny + delta_rhow*nz
    
    !< �������ͨ����ɢ���� Vk = omega
    dtilde1 = 0.5*(lambda3 + lambda4)/delta
    dtilde2 = 0.5*(lambda3 - lambda4)/delta
    dtilde3 = (lambda1 - delta*dtilde1)/(delta*delta)
    dtilde4 = dtilde2 + dtilde3*(oasqr - 1.0)*omega
    dtilde5 = dtilde2 + dtilde3*omega
    dtilde6 = oasqr*dtilde2*omega + dtilde3*sn
    
    dissipation(1) = (dtilde1*delta + dtilde2*(0. - omega))*delta_p + dtilde2*alpha
    dissipation(2) = lambda1*delta_rhou + (dtilde5*nx - qi(2)*dtilde6)*delta_p + (dtilde4*qi(2) - dtilde3*nx)*alpha
    dissipation(3) = lambda1*delta_rhov + (dtilde5*ny - qi(3)*dtilde6)*delta_p + (dtilde4*qi(3) - dtilde3*ny)*alpha
    dissipation(4) = lambda1*delta_rhow + (dtilde5*nz - qi(4)*dtilde6)*delta_p + (dtilde4*qi(4) - dtilde3*nz)*alpha

    do m = 1,4
        flr(m) = 0.5*(flr(m) - dissipation(m))
    end do

end subroutine flux_roe_scmp
    
subroutine flux_slau(nsp,nep,priml,primr,nt,nx,ny,nz,nsf,nef,flr,efix)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,one,two,half,fourth,sml_ssf
    use mod_variables, only : gamma
    implicit none
    integer(kind_int), intent(in)  :: nsp,nep
    real(kind_real)  , intent(in)  :: priml(nsp:nep),primr(nsp:nep)
    real(kind_real)  , intent(in)  :: nt,nx,ny,nz
    integer(kind_int), intent(in)  :: nsf,nef
    real(kind_real)  , intent(out) :: flr(nsf:nef)
    real(kind_real)  , intent(in)  :: efix
    real(kind_real), external :: pressure_function
    
    integer(kind_int) :: flag
    real(kind_real) :: gama,ae
    real(kind_real) :: rl,ul,vl,wl,pl,rr,ur,vr,wr,pr
    real(kind_real) :: vnl,cl2,cl,Ml,vnr,cr2,cr,Mr
    real(kind_real) :: sn,osn,nxa,nya,nza
    real(kind_real) :: cBar,vnBar,g,vnBarPlus,vnBarMinus
    real(kind_real) :: mach,qta,mdot,pTilde
    real(kind_real) :: alpha,fpl,fpr
    real(kind_real) :: psiL(1:5),psiR(1:5),normal(1:5)
    
    gama = gamma
    ae = gama - 1.
    
    rl = priml(1)
    ul = priml(2)
    vl = priml(3)
    wl = priml(4)
    pl = priml(5)
    
    rr = primr(1)
    ur = primr(2)
    vr = primr(3)
    wr = primr(4)
    pr = primr(5)
    
    sn = sqrt(nx*nx + ny*ny + nz*nz)
    osn = 1.0/max(sn,sml_ssf)
    nxa = nx*osn
    nya = ny*osn
    nza = nz*osn
    
    vnl = ul*nxa + vl*nya + wl*nza
    cl2 = gama*pl/rl
    cl  = sqrt(cl2)
    !Ml  = vnl/cl
    
    vnr = ur*nxa + vr*nya + wr*nza
    cr2 = gama*pr/rr
    cr  = sqrt(cr2)
    !Mr  = vnr/cr
    
    cBar = 0.5*(cl + cr)
    Ml = vnl/cBar
    Mr = vnr/cBar
    vnBar = (rl*abs(vnl) + rr*abs(vnr))/(rl + rr)
    
    g = -max(min(Ml, 0.), -1.)*min(max(Mr, 0.), 1.)
    
    vnBarPlus  = (1. - g)*vnBar + g*abs(vnl)
    vnBarMinus = (1. - g)*vnBar + g*abs(vnr)
    
    mach = min(1.0, sqrt(0.5*(ul*ul + vl*vl + wl*wl + ur*ur + vr*vr + wr*wr))/cBar)
    qta = (1-mach)**2
    mdot = 0.5*(rl*(vnl + vnBarPlus) + rr*(vnr - vnBarMinus) - qta/cBar*(pr - pl))
    !mdot = 0.5*(rl*vnl + rr*vnr - vnBar*(rr - rl))*(1 - g) - 0.5*qta/cBar*(pr - pl)
    
    alpha = 0.
    flag = 1
    fpl = pressure_function(Ml, alpha, flag)
    flag = 0
    fpr = pressure_function(Mr, alpha, flag)

    pTilde = 0.5*(pl + pr) + 0.5*(fpl - fpr)*(pl - pr) + (1 - qta)*(fpl + fpr - 1.)*(pl + pr)*0.5
    
    psiL(1) = 1
    psiL(2) = ul
    psiL(3) = vl
    psiL(4) = wl
    psiL(5) = cl2/ae + 0.5*(ul*ul + vl*vl + wl*wl)
    
    psiR(1) = 1
    psiR(2) = ur
    psiR(3) = vr
    psiR(4) = wr
    psiR(5) = cr2/ae + 0.5*(ur*ur + vr*vr + wr*wr)
    
    normal(1) = 0
    normal(2) = nxa
    normal(3) = nya
    normal(4) = nza
    normal(5) = 0
    
    flr(:) = 0.5*(mdot + abs(mdot))*psiL(:) + 0.5*(mdot - abs(mdot))*psiR(:) + pTilde*normal(:)
    
    flr(:) = flr(:)*sn

end subroutine flux_slau
    
function pressure_function(mach,alpha,flag) result(fp)
    use mod_kndconsts, only : kind_int,kind_real
    implicit none
    integer(kind_int) :: flag
    real(kind_real)   :: mach,alpha
    real(kind_real)   :: fp_plus,fp_minus,fp
    
    if (abs(mach) > 1.) then
        !fp_plus = 0.5*(1+sign(1., mach))
        !fp_minus = 0.5*(1+sign(1., -mach))
    else
        fp_plus  = 0.25*(mach+1.)**2*(2-mach)+alpha*mach*(mach**2-1)**2  
        fp_minus = 0.25*(mach-1.)**2*(2+mach)-alpha*mach*(mach**2-1)**2  
    end if
    
    if (flag == 1) then
        fp = fp_plus
    else
        fp = fp_minus
    end if
    
end function

subroutine mxdq_sw(nsp,nep,prim,nt,nx,ny,nz,nsf,nef,dq,df,efix,fsw)
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
    real(kind_real),   intent(in)  :: efix,fsw
    real(kind_real)           :: ro,vx,vy,vz,ps,gama,ae
    real(kind_real)           :: c,c2,h,h0,v2,vn,af,eps,eps2
    real(kind_real)           :: sn,osn,nxa,nya,nza,vna
    real(kind_real)           :: l1,l4,l5,x1,x2,dc,dh,c2dc
    real(kind_real)           :: al1,al4,al5,ma
    real(kind_real), external :: enfix_steger,enfix_harten
    integer(kind_int)         :: nfix=0

    ro = prim(nsp)
    vx = prim(nsp+1)
    vy = prim(nsp+2)
    vz = prim(nsp+3)
    ps = prim(nsp+4)

    gama = gamma
    ae   = gama - one

    v2 = vx*vx + vy*vy + vz*vz
    c2 = gama*ps/ro
    h  = c2/ae
    h0 = h + half*v2

    c  = sqrt(c2)
    vn = nx*vx + ny*vy + nz*vz

    sn = max(sqrt(nx*nx + ny*ny + nz*nz),sml_ssf)
    osn = one/sn
    nxa = nx * osn
    nya = ny * osn
    nza = nz * osn
    vna = vn * osn

    l1 = vn
    l4 = vn + sn*c
    l5 = vn - sn*c

    ma = sqrt(v2)/c
    eps = efix*sn*c*(one + ma)

    if (nfix == 0) then
        al1 = enfix_steger(l1,eps)
        al4 = enfix_steger(l4,eps)
        al5 = enfix_steger(l5,eps)
    else
        al1 = enfix_harten(l1,eps)
        al4 = enfix_harten(l4,eps)
        al5 = enfix_harten(l5,eps)
    end if

    l1 = half*(l1 + fsw*al1)
    l4 = half*(l4 + fsw*al4)
    l5 = half*(l5 + fsw*al5)

    x1 = ( two*l1 - l4 - l5 )/( two * c2 )
    x2 = ( l4 - l5 )/( two * c )


    af = half*ae*v2

    dc = vna * dq(nsf) -        nxa * dq(nsf+1) - nya * dq(nsf+2) - nza * dq(nsf+3)
    dh = af  * dq(nsf) - ae * (  vx * dq(nsf+1) +  vy * dq(nsf+2) +  vz * dq(nsf+3) - dq(nsf+4) )

    c2dc = c2 * dc

    df(nsf  ) = l1 * dq(nsf  )              -    dh   * x1            -    dc   * x2
    df(nsf+1) = l1 * dq(nsf+1) + ( nxa*c2dc - vx*dh ) * x1 + ( nxa*dh - vx*dc ) * x2
    df(nsf+2) = l1 * dq(nsf+2) + ( nya*c2dc - vy*dh ) * x1 + ( nya*dh - vy*dc ) * x2
    df(nsf+3) = l1 * dq(nsf+3) + ( nza*c2dc - vz*dh ) * x1 + ( nza*dh - vz*dc ) * x2
    df(nsf+4) = l1 * dq(nsf+4) + ( vna*c2dc - h0*dh ) * x1 + ( vna*dh - h0*dc ) * x2

end subroutine mxdq_sw

subroutine mxdq_std(nsp,nep,prim,nt,nx,ny,nz,nsf,nef,dq,df,srad,fsw)
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
    real(kind_real),   intent(in)  :: srad
    real(kind_real),   intent(in)  :: fsw
    real(kind_real) :: ro,vx,vy,vz,ps,gama,ae
    real(kind_real) :: c,c2,h,h0,v2,vn,af
    real(kind_real) :: sn,osn,nxa,nya,nza,vna
    real(kind_real) :: l1,l4,l5,x1,x2,y1,y2,rcl,sl

    ro = prim(nsp)
    vx = prim(nsp+1)
    vy = prim(nsp+2)
    vz = prim(nsp+3)
    ps = prim(nsp+4)

    gama = gamma
    ae   = gama - one

    v2 = vx*vx + vy*vy + vz*vz
    c2 = gama*ps/ro
    h  = c2/ae
    h0 = h + half*v2

    c  = sqrt(c2)
    vn = nx*vx + ny*vy + nz*vz

    sn  = max(sqrt(nx*nx + ny*ny + nz*nz),sml_ssf)
    osn = one/sn
    nxa = nx * osn
    nya = ny * osn
    nza = nz * osn
    vna = vn * osn

    l1 = vn
    l4 = vn + sn*c
    l5 = vn - sn*c

    l1 = half*(l1 + fsw*srad)!�����װ뾶����
    l4 = half*(l4 + fsw*srad)
    l5 = half*(l5 + fsw*srad)

    x1 = ( two*l1 - l4 - l5 )/( two * c2 )
    x2 = ( l4 - l5 )/( two * c )

    af = half*ae*v2

    y1 = vna * dq(nsf) -      ( nxa * dq(nsf+1) + nya * dq(nsf+2) + nza * dq(nsf+3) )
    y2 = af  * dq(nsf) - ae * ( vx  * dq(nsf+1) +  vy * dq(nsf+2) +  vz * dq(nsf+3) - dq(nsf+4) )

    rcl = c2*x1*y1 + x2*y2
    sl  = x1*y2 + x2*y1

    df(nsf  ) = l1 * dq(nsf  ) -    sl
    df(nsf+1) = l1 * dq(nsf+1) - vx*sl + nxa*rcl
    df(nsf+2) = l1 * dq(nsf+2) - vy*sl + nya*rcl
    df(nsf+3) = l1 * dq(nsf+3) - vz*sl + nza*rcl
    df(nsf+4) = l1 * dq(nsf+4) - h0*sl + vna*rcl

end subroutine mxdq_std
    
subroutine mxdq_scmp(nsp,nep,prim,nt,nx,ny,nz,nsf,nef,dq,df,srad,fsw)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : one,half,two,sml_ssf,SCMP_sigma
    use mod_variables, only : gamma,moo
    implicit none
    integer(kind_int), intent(in)  :: nsp,nep
    real(kind_real),   intent(in)  :: prim(nsp:nep)
    real(kind_real),   intent(in)  :: nt,nx,ny,nz
    integer(kind_int), intent(in)  :: nsf,nef
    real(kind_real),   intent(in)  :: dq(nsf:nef)
    real(kind_real),   intent(out) :: df(nsf:nef)
    real(kind_real),   intent(in)  :: srad
    real(kind_real),   intent(in)  :: fsw
    real(kind_real) :: ro,vx,vy,vz,ps,gama,ae
    real(kind_real) :: c,c2,oc2,h,h0,v2,vn,af
    real(kind_real) :: sn,osn,sn2,nxa,nya,nza,vna
    real(kind_real) :: l1,l3,l4,x1,x2,y1,y2,rcl,sl
    real(kind_real) :: Delta
    real(kind_real) :: D1,D2,D3,D4,D5,D6

    ro = prim(nsp)
    vx = prim(nsp+1)
    vy = prim(nsp+2)
    vz = prim(nsp+3)
    ps = prim(nsp+4)

    gama = gamma

    c2 = SCMP_sigma / (gama*moo*moo*ro**(1-SCMP_sigma))    !> Done
    oc2 = 1.0/c2

    sn  = max(sqrt(nx*nx + ny*ny + nz*nz),sml_ssf)
    sn2 = sn*sn
    osn = one/sn

    vn = nx*vx + ny*vy + nz*vz
    Delta = sqrt( vn*vn*(1 - oc2) + sn2 )     !> Done

    l1 = vn
    l3 = vn + Delta
    l4 = vn - Delta

    l1 = half*(l1 + fsw*srad)!�����װ뾶����
    l3 = half*(l3 + fsw*srad)
    l4 = half*(l4 + fsw*srad)
    
    D1 = (l3+l4)/(2*Delta)
    D2 = (l3-l4)/(2*Delta)
    D3 = (l1-Delta*D1)/(Delta*Delta)
    D4 = D2 + D3*(oc2-1)*vn
    D5 = D2 + D3*vn
    D6 = D2*vn*oc2 + D3*sn2
    
    df(nsf  ) = dq(1)*(D1*Delta - D2*vn) + dq(2)*(D2*nx)                   + &
                dq(3)*(D2*ny)                   + dq(4)*(D2*nz)
    df(nsf+1) = dq(1)*(D5*nx    - vx*D6) + dq(2)*(l1 + (D4*vx - D3*nx)*nx) + &
                dq(3)*(     (D4*vx - D3*nx)*ny) + dq(4)*(     (D4*vx - D3*nx)*nz)
    df(nsf+2) = dq(1)*(D5*ny    - vy*D6) + dq(2)*(     (D4*vy - D3*ny)*nx) + &
                dq(3)*(l1 + (D4*vy - D3*ny)*ny) + dq(4)*(     (D4*vy - D3*ny)*nz)
    df(nsf+3) = dq(1)*(D5*nz    - vz*D6) + dq(2)*(     (D4*vz - D3*nz)*nx) + &
                dq(3)*(     (D4*vz - D3*nz)*ny) + dq(4)*(l1 + (D4*vz - D3*nz)*nz)

end subroutine mxdq_scmp

subroutine aflux_inv(nsp,nep,prim,nt,nx,ny,nz,nsf,nef,af,efix)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : one,half,two,sml_ssf
    use mod_variables, only : gamma
    implicit none
    integer(kind_int), intent(in)  :: nsp,nep
    real(kind_real),   intent(in)  :: prim(nsp:nep)
    real(kind_real),   intent(in)  :: nt,nx,ny,nz
    integer(kind_int), intent(in)  :: nsf,nef
    real(kind_real),   intent(out) :: af(nsf:nef,nsf:nef)
    real(kind_real),   intent(in)  :: efix
    real(kind_real)           :: ro,vx,vy,vz,ps,gama,ae
    real(kind_real)           :: a,a2,h,h0,v2,vn,eps
    real(kind_real)           :: sn,osn,nxa,nya,nza,vna,almax
    real(kind_real)           :: l1,l4,l5,al1,al4,al5,x1,x2,x3,y1
    real(kind_real)           :: c1,d1,c2,d2,c3,d3,c4,d4,c5,d5
    real(kind_real), external :: enfix_steger,enfix_harten
    integer(kind_int)         :: nfix=0

    ro = prim(nsp)
    vx = prim(nsp+1)
    vy = prim(nsp+2)
    vz = prim(nsp+3)
    ps = prim(nsp+4)

    gama = gamma
    ae   = gama - one

    v2 = vx*vx + vy*vy + vz*vz
    a2 = gama*ps/ro
    h  = a2/ae
    h0 = h + half*v2

    a  = sqrt(a2)
    vn = nx*vx + ny*vy + nz*vz

    sn  = max(sqrt(nx*nx + ny*ny + nz*nz),sml_ssf)
    osn = one/sn
    nxa = nx * osn
    nya = ny * osn
    nza = nz * osn
    vna = vn * osn

    l1 = vn
    l4 = vn + sn*a
    l5 = vn - sn*a

    almax = abs(vn) + sn*a

    eps = efix*sn
    !!eps = efix*almax
    if (nfix == 0) then
        al1 = enfix_steger(l1,eps)
        al4 = enfix_steger(l4,eps)
        al5 = enfix_steger(l5,eps)
    else
        al1 = enfix_harten(l1,eps)
        al4 = enfix_harten(l4,eps)
        al5 = enfix_harten(l5,eps)
    end if

    l1 = al1
    l4 = al4
    l5 = al5

    x1 = half*(l4 + l5)
    x2 = half*(l4 - l5)
    x3 = x1 - l1

    y1 = half*v2

    c1 = ae*x3/a2
    d1 = x2/a
    af(nsf  ,nsf  ) =  c1*y1 - d1*vna + l1
    af(nsf  ,nsf+1) = -c1*vx + d1*nxa
    af(nsf  ,nsf+2) = -c1*vy + d1*nya
    af(nsf  ,nsf+3) = -c1*vz + d1*nza
    af(nsf  ,nsf+4) =  c1

    c2 = c1*vx  + d1*nxa*ae
    d2 = x3*nxa + d1*vx
    af(nsf+1,nsf  ) =  c2*y1 - d2*vna
    af(nsf+1,nsf+1) = -c2*vx + d2*nxa + l1
    af(nsf+1,nsf+2) = -c2*vy + d2*nya
    af(nsf+1,nsf+3) = -c2*vz + d2*nza
    af(nsf+1,nsf+4) =  c2

    c3 = c1*vy  + d1*nya*ae
    d3 = x3*nya + d1*vy
    af(nsf+2,nsf  ) =  c3*y1 - d3*vna
    af(nsf+2,nsf+1) = -c3*vx + d3*nxa
    af(nsf+2,nsf+2) = -c3*vy + d3*nya + l1
    af(nsf+2,nsf+3) = -c3*vz + d3*nza
    af(nsf+2,nsf+4) =  c3

    c4 = c1*vz  + d1*nza*ae
    d4 = x3*nza + d1*vz
    af(nsf+3,nsf  ) =  c4*y1 - d4*vna
    af(nsf+3,nsf+1) = -c4*vx + d4*nxa
    af(nsf+3,nsf+2) = -c4*vy + d4*nya
    af(nsf+3,nsf+3) = -c4*vz + d4*nza + l1
    af(nsf+3,nsf+4) =  c4

    c5 = c1*h0  + d1*vna*ae
    d5 = x3*vna + d1*h0
    af(nsf+4,nsf  ) =  c5*y1 - d5*vna
    af(nsf+4,nsf+1) = -c5*vx + d5*nxa
    af(nsf+4,nsf+2) = -c5*vy + d5*nya
    af(nsf+4,nsf+3) = -c5*vz + d5*nza
    af(nsf+4,nsf+4) =  c5             + l1

end subroutine aflux_inv

subroutine aflux_vis(nsp,nep,prim,nt,nx,ny,nz,vol,vsl,vst,nsf,nef,af)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,third,half,one,two,sml_ssf
    use mod_variables, only : refbeta,gamma,prlam,prtur
    use mod_variables, only : reue,csrvis
    implicit none
    integer(kind_int), intent(in)  :: nsp,nep
    real(kind_real),   intent(in)  :: prim(nsp:nep)
    real(kind_real),   intent(in)  :: nt,nx,ny,nz,vol
    real(kind_real),   intent(in)  :: vsl,vst
    integer(kind_int), intent(in)  :: nsf,nef
    real(kind_real),   intent(out) :: af(nsf:nef,nsf:nef)
    integer(kind_int) :: m,n
    real(kind_real)   :: ro,vx,vy,vz,ps,vn,v2,e0
    real(kind_real)   :: sn2,sn,osn,nxa,nya,nza,vna
    real(kind_real)   :: gama,ae,cv,cp,cp_prl,cp_prt
    real(kind_real)   :: nx3,ny3,nz3,kcp,vis,c5,wc5,cav

    gama = gamma
    ae = gama - one

    cv = refbeta/ae
    cp = gama*cv
    cp_prl = cp/prlam
    cp_prt = cp/prtur

    ro = prim(nsp)
    vx = prim(nsp+1)
    vy = prim(nsp+2)
    vz = prim(nsp+3)
    ps = prim(nsp+4)

    vn = nx*vx + ny*vy + nz*vz
    v2 = vx*vx + vy*vy + vz*vz
    e0 = ps/(ae*ro) + half*v2

    sn2 = nx*nx + ny*ny + nz*nz
    sn  = max(sqrt(sn2),sml_ssf)
    osn = one/sn
    nxa = nx * osn
    nya = ny * osn
    nza = nz * osn
    vna = vn * osn

    nx3 = third*nxa
    ny3 = third*nya
    nz3 = third*nza

! CHEyg   ������μ�������Ԫ�ط���˳����Ż�һЩ�� ���Ż�����openmpӰ��
!#define AFLUX_VIS_OPT  ��makefile�ж���

#ifndef AFLUX_VIS_OPT
! ԭ���Ĵ���
    af(nsf  ,nsf  ) = zero
    af(nsf  ,nsf+1) = zero
    af(nsf  ,nsf+2) = zero
    af(nsf  ,nsf+3) = zero
    af(nsf  ,nsf+4) = zero

    af(nsf+1,nsf  ) = -(nx3*vna + vx)
    af(nsf+1,nsf+1) = nx3*nxa + one
    af(nsf+1,nsf+2) = nx3*nya
    af(nsf+1,nsf+3) = nx3*nza
    af(nsf+1,nsf+4) = zero

    af(nsf+2,nsf  ) = -(ny3*vna + vy)
    af(nsf+2,nsf+1) = ny3*nxa
    af(nsf+2,nsf+2) = ny3*nya + one
    af(nsf+2,nsf+3) = ny3*nza
    af(nsf+2,nsf+4) = zero

    af(nsf+3,nsf  ) = -(nz3*vna + vz)
    af(nsf+3,nsf+1) = nz3*nxa
    af(nsf+3,nsf+2) = nz3*nya
    af(nsf+3,nsf+3) = nz3*nza + one
    af(nsf+3,nsf+4) = zero

#else
! �Ż��Ĵ���
    af(nsf  ,nsf  ) = zero
    af(nsf+1,nsf  ) = -(nx3*vna + vx)
    af(nsf+2,nsf  ) = -(ny3*vna + vy)
    af(nsf+3,nsf  ) = -(nz3*vna + vz)
                
    af(nsf  ,nsf+1) = zero
    af(nsf+1,nsf+1) = nx3*nxa + one    
    af(nsf+2,nsf+1) = ny3*nxa
    af(nsf+3,nsf+1) = nz3*nxa
    
            
    af(nsf  ,nsf+2) = zero
    af(nsf+1,nsf+2) = nx3*nya
    af(nsf+2,nsf+2) = ny3*nya + one
    af(nsf+3,nsf+2) = nz3*nya    
            
    af(nsf  ,nsf+3) = zero
    af(nsf+1,nsf+3) = nx3*nza
    af(nsf+2,nsf+3) = ny3*nza
    af(nsf+3,nsf+3) = nz3*nza + one    
            
    af(nsf  ,nsf+4) = zero
    af(nsf+1,nsf+4) = zero
    af(nsf+2,nsf+4) = zero
    af(nsf+3,nsf+4) = zero

#endif    

    vis = vsl + vst
    kcp = vsl*cp_prl + vst*cp_prt

    c5 = kcp/(cv*vis)
    wc5 = one - c5
    af(nsf+4,nsf  ) = -(third*vna*vna + wc5*v2 + c5*e0)
    af(nsf+4,nsf+1) = nx3*vna + vx*wc5
    af(nsf+4,nsf+2) = ny3*vna + vy*wc5
    af(nsf+4,nsf+3) = nz3*vna + vz*wc5
    af(nsf+4,nsf+4) = c5

    cav = csrvis*two*vis*sn2/(reue*ro*vol)
    do m=nsf,nef
    do n=nsf,nef
        af(m,n) = af(m,n)*cav
    end do
    end do

end subroutine aflux_vis

subroutine LxQ_p(nsp,nep,prim,nt,nx,ny,nz,q,f)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : one,sml_ssf
    use mod_variables, only : gamma
    implicit none
    integer(kind_int), intent(in)  :: nsp,nep
    real(kind_real),   intent(in)  :: prim(nsp:nep)
    real(kind_real),   intent(in)  :: nt,nx,ny,nz
    real(kind_real),   intent(in)  :: q(nsp:nep)
    real(kind_real),   intent(out) :: f(nsp:nep)
    real(kind_real) :: gama,qn,wp1,wp2
    real(kind_real) :: nta,nxa,nya,nza,osn
    real(kind_real) :: ro,vx,vy,vz,ps,c2,c

    gama = gamma

    osn = one/max(sqrt(nx*nx + ny*ny + nz*nz),sml_ssf)
    nta = nt*osn
    nxa = nx*osn
    nya = ny*osn
    nza = nz*osn

    ro = prim(nsp  )
    vx = prim(nsp+1)
    vy = prim(nsp+2)
    vz = prim(nsp+3)
    ps = prim(nsp+4)

    c2 = gama*ps/ro
    c  = sqrt(c2)

    qn  = nxa*q(2) + nya*q(3) + nza*q(4)
    wp1 = q(1) - q(5)/c2
    wp2 = q(5)/(ro*c)

    f(nsp  ) = nza*q(3) - nya*q(4) + nxa*wp1
    f(nsp+1) = nxa*q(4) - nza*q(2) + nya*wp1
    f(nsp+2) = nya*q(2) - nxa*q(3) + nza*wp1
    f(nsp+3) = wp2+qn
    f(nsp+4) = wp2-qn

end subroutine LxQ_p

subroutine RxQ_p(nsp,nep,prim,nt,nx,ny,nz,q,f)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : one,half,sml_ssf
    use mod_variables, only : gamma
    implicit none
    integer(kind_int), intent(in)  :: nsp,nep
    real(kind_real),   intent(in)  :: prim(nsp:nep)
    real(kind_real),   intent(in)  :: nt,nx,ny,nz
    real(kind_real),   intent(in)  :: q(nsp:nep)
    real(kind_real),   intent(out) :: f(nsp:nep)
    real(kind_real) :: gama,qn,wp1,wp2
    real(kind_real) :: nta,nxa,nya,nza,osn
    real(kind_real) :: ro,vx,vy,vz,ps,c2,c

    gama = gamma

    osn = one/max(sqrt(nx*nx + ny*ny + nz*nz),sml_ssf)
    nta = nt*osn
    nxa = nx*osn
    nya = ny*osn
    nza = nz*osn

    ro = prim(nsp  )
    vx = prim(nsp+1)
    vy = prim(nsp+2)
    vz = prim(nsp+3)
    ps = prim(nsp+4)

    c2 = gama*ps/ro
    c  = sqrt(c2)

    qn  = nxa*q(1) + nya*q(2) + nza*q(3)
    wp1 = half*ro*(q(4) + q(5))
    wp2 = half*   (q(4) - q(5))

    f(nsp  ) =     qn              + wp1/c
    f(nsp+1) = nya*q(3) - nza*q(2) + nxa*wp2
    f(nsp+2) = nza*q(1) - nxa*q(3) + nya*wp2
    f(nsp+3) = nxa*q(2) - nya*q(1) + nza*wp2
    f(nsp+4) =                       wp1*c

end subroutine RxQ_p

subroutine t2visl(nst,net,tem,nsv,nev,visl)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_variables, only : tssnd
    implicit none
    integer(kind_int), intent(in)  :: nst,net
    real(kind_real)  , intent(in)  :: tem(nst:net)
    integer(kind_int), intent(in)  :: nsv,nev
    real(kind_real)  , intent(out) :: visl(nsv:nev)
    real(kind_real), external :: sutherland

    visl(nsv) = sutherland(tem(nst),tssnd)

end subroutine t2visl

subroutine t2visl_SCMP(nst,net,tem,nsv,nev,visl)        !> Done
    use mod_kndconsts, only : kind_int,kind_real
    use mod_variables, only : tssnd
    implicit none
    integer(kind_int), intent(in)  :: nst,net
    real(kind_real)  , intent(in)  :: tem(nst:net)
    integer(kind_int), intent(in)  :: nsv,nev
    real(kind_real)  , intent(out) :: visl(nsv:nev)

    visl(nsv) = 1.0

end subroutine t2visl_SCMP
    

subroutine pv2cv_SCMP ( npvs, nx,ny,nz, primRef, prim, chrt )
    use mod_kndconsts, only : kind_int,kind_real
    use mod_variables, only : gamma,poo,moo
    use mod_constants, only : scmp_sigma,sml_ssf
    implicit none
    integer(kind_int), intent(in)  :: npvs
    real(kind_real)  , intent(in)  :: nx,ny,nz
    real(kind_real)  , intent(in)  :: primRef(npvs)
    real(kind_real)  , intent(in)  :: prim(npvs)
    real(kind_real)  , intent(out) :: chrt(4)
    real(kind_real)   :: ro,vx,vy,vz,p,pPrime,gama
    real(kind_real)   :: vecV(3), vecN(3), vecL(3), vecM(3)
    real(kind_real)   :: Vn, Vl,lx,ly,lz, Vm,mx,my,mz
    real(kind_real)   :: d1,d2,d3,d4, Delta, oa2, sn2

    gama = gamma

    vx = primRef(2)
    vy = primRef(3)
    vz = primRef(4)
    p  = primRef(5)
    pPrime = p - poo
    ro  = (1.0 + gama*moo*moo*pPrime)**(1.0/SCMP_sigma)
    oa2 = (gama*moo*moo*ro**(1-SCMP_sigma)) / SCMP_sigma
    
    vecV(1) = vx
    vecV(2) = vy
    vecV(3) = vz
    
    vecN(1) = nx
    vecN(2) = ny
    vecN(3) = nz
    !call normalizeVector_dim3 ( vecN )
    !nx = vecN(1)
    !ny = vecN(2)
    !nz = vecN(3)
    
    call cross_Product_dim3 ( vecN, vecV, vecL )    !> l = k x v
    call cross_Product_dim3 ( vecN, vecL, vecM )    !> m = k x l
    call normalizeVector_dim3 ( vecL )
    call normalizeVector_dim3 ( vecM )
    lx = vecL(1)
    ly = vecL(2)
    lz = vecL(3)
    mx = vecM(1)
    my = vecM(2)
    mz = vecM(3)
    
    Vn = vx*nx + vy*ny + vz*nz
    Vl = vx*lx + vy*ly + vz*lz
    Vm = vx*mx + vy*my + vz*mz
    
    sn2 = nx*nx + ny*ny + nz*nz
    sn2 = max(sn2, sml_ssf)
    Delta = sqrt( Vn*Vn*(1 - oa2) + sn2 )     !> Done
    
    d1 = (oa2 - 1.)*Vn*Vl/(Delta**2)
    d2 = (oa2 - 1.)*Vn*Vm/(Delta**2)
    d3 = 0.5/(Delta*Delta)
    d4 = d3
    
    vx = prim(2)
    vy = prim(3)
    vz = prim(4)
    p  = prim(5)
    pPrime = p - poo
    ro  = (1.0 + gama*moo*moo*pPrime)**(1.0/SCMP_sigma)
    
    chrt(1) = -(Vl + d1*Vn)    *pPrime + (lx + d1*nx)*ro*vx + (ly + d1*ny)*ro*vy + (lz + d1*nz)*ro*vz
    chrt(2) = -(Vm + d2*Vn)    *pPrime + (mx + d2*nx)*ro*vx + (my + d2*ny)*ro*vy + (mz + d2*nz)*ro*vz
    chrt(3) = ( Delta - Vn)*d3 *pPrime +       d3*nx *ro*vx +       d3*ny *ro*vy +       d3*nz *ro*vz
    chrt(4) = (-Delta - Vn)*d4 *pPrime +       d4*nx *ro*vx +       d4*ny *ro*vy +       d4*nz *ro*vz

end subroutine pv2cv_SCMP


subroutine cv2pv_SCMP ( npvs, nx,ny,nz, primRef, chrt, primOut )
    use mod_kndconsts, only : kind_int,kind_real
    use mod_variables, only : gamma,poo,moo
    use mod_constants, only : scmp_sigma,sml_ssf
    implicit none
    integer(kind_int), intent(in)  :: npvs
    real(kind_real)  , intent(in)  :: nx,ny,nz
    real(kind_real)  , intent(in)  :: primRef(npvs)
    real(kind_real)  , intent(in)  :: chrt(4)
    real(kind_real)  , intent(out) :: primOut(npvs)
    real(kind_real)   :: ro,vx,vy,vz,p,pPrime,gama
    real(kind_real)   :: vecV(3), vecN(3), vecL(3), vecM(3)
    real(kind_real)   :: Vn, Vl,lx,ly,lz, Vm,mx,my,mz
    real(kind_real)   :: t3,t4, Delta, oa2, sn2
    real(kind_real)   :: consT(4)

    gama = gamma

    vx = primRef(2)
    vy = primRef(3)
    vz = primRef(4)
    p  = primRef(5)
    pPrime = p - poo
    ro  = (1.0 + gama*moo*moo*pPrime)**(1.0/SCMP_sigma)
    oa2 = (gama*moo*moo*ro**(1-SCMP_sigma)) / SCMP_sigma
    
    vecV(1) = vx
    vecV(2) = vy
    vecV(3) = vz
    
    vecN(1) = nx
    vecN(2) = ny
    vecN(3) = nz
    !call normalizeVector_dim3 ( vecN )
    !nx = vecN(1)
    !ny = vecN(2)
    !nz = vecN(3)
    
    call cross_Product_dim3 ( vecN, vecV, vecL )    !> l = k x v
    call cross_Product_dim3 ( vecN, vecL, vecM )    !> m = k x l
    call normalizeVector_dim3 ( vecL )
    call normalizeVector_dim3 ( vecM )
    lx = vecL(1)
    ly = vecL(2)
    lz = vecL(3)
    mx = vecM(1)
    my = vecM(2)
    mz = vecM(3)
    
    Vn = vx*nx + vy*ny + vz*nz
    Vl = vx*lx + vy*ly + vz*lz
    Vm = vx*mx + vy*my + vz*mz
    
    sn2 = nx*nx + ny*ny + nz*nz
    sn2 = max(sn2, sml_ssf)
    Delta = sqrt( Vn*Vn*(1 - oa2) + sn2 )     !> Done
    
    t3 = (1 - oa2)*Vn + Delta
    t4 = (1 - oa2)*Vn - Delta
    
    consT(1) =                                 Delta *chrt(3) -       Delta *chrt(4)
    consT(2) = lx*chrt(1) + mx*chrt(2) + (nx + t3*vx)*chrt(3) + (nx + t4*vx)*chrt(4)
    consT(3) = ly*chrt(1) + my*chrt(2) + (ny + t3*vy)*chrt(3) + (ny + t4*vy)*chrt(4)
    consT(4) = lz*chrt(1) + mz*chrt(2) + (nz + t3*vz)*chrt(3) + (nz + t4*vz)*chrt(4)
    
    pPrime = consT(1)
    p  = pPrime + poo
    ro = (1.0 + gama*moo*moo*pPrime)**(1.0/SCMP_sigma)
    
    primOut(1) = ro
    primOut(2) = consT(2)/ro
    primOut(3) = consT(3)/ro
    primOut(4) = consT(4)/ro
    primOut(5) = p

end subroutine cv2pv_SCMP


!> ���
subroutine cross_Product_dim3 ( vecA, vecB, vecC )
    use mod_kndconsts, only : kind_real
    implicit none
    real(kind_real)  , intent(in)  :: vecA(3)
    real(kind_real)  , intent(in)  :: vecB(3)
    real(kind_real)  , intent(out) :: vecC(3)
    real(kind_real)   :: a1,a2,a3, b1,b2,b3

    a1 = vecA(1)
    a2 = vecA(2)
    a3 = vecA(3)
    
    b1 = vecB(1)
    b2 = vecB(2)
    b3 = vecB(3)
    
    vecC(1) = a2*b3 - a3*b2
    vecC(2) = a3*b1 - a1*b3
    vecC(3) = a1*b2 - a2*b1

end subroutine cross_Product_dim3

!> ������λ��
subroutine normalizeVector_dim3 ( vec )
    use mod_kndconsts, only : kind_real
    implicit none
    real(kind_real)  , intent(inout) :: vec(3)
    real(kind_real)   :: mod

    mod = sqrt( vec(1)**2 + vec(2)**2 + vec(3)**2 )
    
    vec(1) = vec(1)/mod
    vec(2) = vec(2)/mod
    vec(3) = vec(3)/mod

end subroutine normalizeVector_dim3

subroutine flux_vis(vx,vy,vz,nt,nx,ny,nz,kcp,vis,nsd,ned,der,nsf,nef,f)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,one,two,two3rd
    use mod_variables, only : gamma
    implicit none
    real(kind_real)  , intent(in)  :: vx,vy,vz
    real(kind_real)  , intent(in)  :: nt,nx,ny,nz
    real(kind_real)  , intent(in)  :: kcp,vis
    integer(kind_int), intent(in)  :: nsd,ned
    real(kind_real)  , intent(in)  :: der(nsd:ned)
    integer(kind_int), intent(in)  :: nsf,nef
    real(kind_real)  , intent(out) :: f(nsf:nef)
    real(kind_real) :: dux,duy,duz,dvx,dvy,dvz
    real(kind_real) :: dwx,dwy,dwz,dtx,dty,dtz
    real(kind_real) :: vis2p3,tauxx,tauyy,tauzz
    real(kind_real) :: tauxy,tauxz,tauyz

    dux = der(nsd  )
    duy = der(nsd+1)
    duz = der(nsd+2)

    dvx = der(nsd+3)
    dvy = der(nsd+4)
    dvz = der(nsd+5)

    dwx = der(nsd+6)
    dwy = der(nsd+7)
    dwz = der(nsd+8)

    dtx = der(nsd+9)
    dty = der(nsd+10)
    dtz = der(nsd+11)

    vis2p3 = two3rd*vis
    tauxx = vis2p3 * ( two*dux - dvy - dwz )
    tauyy = vis2p3 * ( two*dvy - dwz - dux )
    tauzz = vis2p3 * ( two*dwz - dux - dvy )
    tauxy = vis * ( duy + dvx )
    tauxz = vis * ( duz + dwx )
    tauyz = vis * ( dvz + dwy )

    f(nsf  ) = zero
    f(nsf+1) = nx*tauxx + ny*tauxy + nz*tauxz
    f(nsf+2) = nx*tauxy + ny*tauyy + nz*tauyz
    f(nsf+3) = nx*tauxz + ny*tauyz + nz*tauzz
    f(nsf+4) = vx*f(nsf+1) + vy*f(nsf+2) + vz*f(nsf+3) + &
               kcp*( nx*dtx + ny*dty + nz*dtz )

end subroutine flux_vis

subroutine set_patched_offsets
    use mod_kndconsts, only : kind_int
    use mod_constants, only : nbc_inter_buf_dyn,nsgl_aver_art
    use mod_variables, only : nghnode
    use mod_fieldvars, only : mb_xyz
    use mod_interface, only : pre_exchange_bc_var,post_exchange_bc_var
    implicit none
        
    call create_buff_offsets        
    
    call calc_patched_indexes
    
    call backup_patched_coordinates
    
    call pre_exchange_bc_var(mb_xyz,1,3,nghnode,nbc_inter_buf_dyn,nsgl_aver_art)
    
    call post_exchange_bc_var(mb_xyz,1,3,nghnode,nbc_inter_buf_dyn,nsgl_aver_art)    
    
    call recover_patched_coordinates

end subroutine set_patched_offsets

subroutine backup_patched_coordinates
    use mod_kndconsts, only : kind_int,kind_real
    use mod_datatypes, only : var_block_t
    use mod_variables, only : nghnode
    use mod_fieldvars, only : mb_top
    use mod_fieldvars, only : ninters,ninterf
    use mod_fieldvars, only : bc_seq
    use mod_fieldvars, only : inters,patched
    use mod_fieldvars, only : mb_xyz
    use mod_parallels
    implicit none
    integer(kind_int)             :: bc_id,nb
    integer(kind_int)             :: s_st(3),s_ed(3)
    integer(kind_int)             :: s_nd,s_lr
    integer(kind_int)             :: st(3),ed(3)
    integer(kind_int)             :: id_src
    integer(kind_int)             :: i,j,k,nint,m

    do nint=ninterf+1,ninters
        s_st(:)  = inters(nint)%bc%s_st(:)
        s_ed(:)  = inters(nint)%bc%s_ed(:)
        s_nd     = inters(nint)%bc%s_nd
        s_lr     = inters(nint)%bc%s_lr

        bc_id   = nint - ninterf
        nb      = bc_seq(bc_id,1)

        call bc_extend_outward(s_st,s_ed,s_nd,s_lr,nghnode,st,ed)

#ifdef PARALLEL
        id_src = mb_top(nb)%pid
        if (myid == id_src) then
#endif

            do m=1,3
            do k=st(3),ed(3)
            do j=st(2),ed(2)
            do i=st(1),ed(1)

                    patched(bc_id)%dats(i,j,k,m) = mb_xyz(nb)%fld(m)%r3d(i,j,k)
        
            end do
            end do
            end do
            end do           

#ifdef PARALLEL
        end if
#endif
    end do

end subroutine backup_patched_coordinates

subroutine recover_patched_coordinates
    use mod_kndconsts, only : kind_int,kind_real
    use mod_datatypes, only : var_block_t
    use mod_variables, only : nghnode
    use mod_fieldvars, only : mb_top
    use mod_fieldvars, only : ninters,ninterf
    use mod_fieldvars, only : bc_seq
    use mod_fieldvars, only : inters,patched
    use mod_fieldvars, only : mb_xyz
    use mod_parallels
    implicit none
    integer(kind_int)             :: bc_id,nb
    integer(kind_int)             :: s_st(3),s_ed(3)
    integer(kind_int)             :: s_nd,s_lr
    integer(kind_int)             :: st(3),ed(3)
    integer(kind_int)             :: id_src
    integer(kind_int)             :: i,j,k,nint,m

    do nint=ninterf+1,ninters
        s_st(:)  = inters(nint)%bc%s_st(:)
        s_ed(:)  = inters(nint)%bc%s_ed(:)
        s_nd     = inters(nint)%bc%s_nd
        s_lr     = inters(nint)%bc%s_lr

        bc_id   = nint - ninterf
        nb      = bc_seq(bc_id,1)

        call bc_extend_outward(s_st,s_ed,s_nd,s_lr,nghnode,st,ed)

#ifdef PARALLEL
        id_src = mb_top(nb)%pid
        if (myid == id_src) then
#endif

            do m=1,3
            do k=st(3),ed(3)
            do j=st(2),ed(2)
            do i=st(1),ed(1)

                    patched(bc_id)%datt(i,j,k,m) = mb_xyz(nb)%fld(m)%r3d(i,j,k)
                    mb_xyz(nb)%fld(m)%r3d(i,j,k) = patched(bc_id)%dats(i,j,k,m)
        
            end do
            end do
            end do
            end do           

#ifdef PARALLEL
        end if
#endif
    end do

end subroutine recover_patched_coordinates

subroutine calc_patched_indexes
    use mod_kndconsts, only : kind_int,kind_real
    use mod_fieldvars, only : mb_top
    use mod_fieldvars, only : ninters,ninterf
    use mod_fieldvars, only : bc_seq,scal
    use mod_fieldvars, only : inters,patched
    use mod_parallels
    implicit none
    integer(kind_int), parameter :: nghp = 5
    integer(kind_int)             :: bc_id,nb,nr
    integer(kind_int)             :: nbs,nbt,nrt,nrs
    integer(kind_int)             :: s_st(3),s_ed(3),t_st(3),t_ed(3)
    integer(kind_int)             :: s_lr,t_lr,s_nd,t_nd
    integer(kind_int)             :: t_dir(3)
    integer(kind_int)             :: id_src
    integer(kind_int)             :: i,j,k,nint,m,n
    integer(kind_int)             :: P1,P2,P1t,P2t
    integer(kind_int)             :: il,jl,kl
    integer(kind_int)             :: ilt,jlt,klt
    integer(kind_int)             :: ijks(3),ijkt(3)
    integer(kind_int)             :: is,js,ks,it,jt,kt

    do nint=ninterf+1,ninters
        nbs      = inters(nint)%bc%nbs
        nrs      = inters(nint)%bc%nrs
        s_st(:)  = inters(nint)%bc%s_st(:)
        s_ed(:)  = inters(nint)%bc%s_ed(:)
        s_nd     = inters(nint)%bc%s_nd
        s_lr     = inters(nint)%bc%s_lr
        
        nbt      = inters(nint)%bc%nbt
        nrt      = inters(nint)%bc%nrt
        t_st(:)  = inters(nint)%bc%t_st(:)
        t_ed(:)  = inters(nint)%bc%t_ed(:)
        t_nd     = inters(nint)%bc%t_nd
        t_lr     = inters(nint)%bc%t_lr
        t_dir(:) = inters(nint)%bc%t_dir(:)

        bc_id   = nint - ninterf
        nb      = bc_seq(bc_id,1)
        nr      = bc_seq(bc_id,2)

#ifdef PARALLEL
        id_src = mb_top(nbs)%pid
        if (myid == id_src) then
#endif
        if (s_nd == 1) then

            do i =s_st(1),s_ed(1)
                do k=s_st(3),s_ed(3)
                do j=s_st(2),s_ed(2)
                    if(scal(bc_id,2) >= 1.0) then
                        if(j<s_st(2)+nghp/2+1) then
                            jl      = -1 + s_st(2) - j
                            P1      = -jl
                        else if(j>s_ed(2)-nghp/2-1) then
                            jl      = -nghp + s_ed(2) - j
                            P1      = -jl
                        else
                            jl      = -nghp/2-1
                            P1      = -jl
                        end if

                        do m=1,nghp
                            patched(bc_id)%id1(i,j,k,m) = j+jl+m
                        end do
                        patched(bc_id)%id1(i,j,k,nghp+1) = j
                    else
                        ijks(:) = (/i,j,k/)
                        ijkt(:) = mb_top(nb)%bcs(nr)%mapijk(i,j,k,:)
                        
                        jt = ijkt(t_dir(2))
                        if(jt<t_st(t_dir(2))+nghp/2+1) then
                            jlt      = -1 + t_st(t_dir(2)) - jt
                            P1t      = -jlt
                        else if(jt>t_ed(t_dir(2))-nghp/2-1) then
                            jlt      = -nghp + t_ed(t_dir(2)) - jt
                            P1t      = -jlt
                        else
                            jlt      = -nghp/2-1
                            P1t      = -jlt
                        end if

                        do m=1,nghp
                            ijkt(t_dir(2)) = jt+jlt+m
                            ijks(:)        = mb_top(nbt)%bcs(nrt)%mapijk(ijkt(1),ijkt(2),ijkt(3),:)

                            patched(bc_id)%id1(i,j,k,m) = ijks(2)
                        end do
                        patched(bc_id)%id1(i,j,k,nghp+1) = j
                    end if
        
                    if(scal(bc_id,3) >= 1.0) then
                        if(k<s_st(3)+nghp/2+1) then
                            kl      = -1 + s_st(3) - k
                            P2      = -kl
                        else if(k>s_ed(3)-nghp/2-1) then
                            kl      = -nghp + s_ed(3) - k
                            P2      = -kl
                        else
                            kl      = -nghp/2-1
                            P2      = -kl
                        end if

                        do m=1,nghp
                            patched(bc_id)%id2(i,j,k,m) = k+kl+m
                        end do
                        patched(bc_id)%id2(i,j,k,nghp+1) = k
                    else
                        ijks(:) = (/i,j,k/)
                        ijkt(:) = mb_top(nb)%bcs(nr)%mapijk(i,j,k,:)

                        kt = ijkt(t_dir(3))
                        if(kt<t_st(t_dir(3))+nghp/2+1) then
                            klt      = -1 + t_st(t_dir(3)) - kt
                            P2t      = -klt
                        else if(kt>t_ed(t_dir(3))-nghp/2-1) then
                            klt      = -nghp + t_ed(t_dir(3)) - kt
                            P2t      = -klt
                        else
                            klt      = -nghp/2-1
                            P2t      = -klt
                        end if

                        do m=1,nghp
                            ijkt(t_dir(3)) = kt+klt+m
                            ijks(:)        = mb_top(nbt)%bcs(nrt)%mapijk(ijkt(1),ijkt(2),ijkt(3),:)

                            patched(bc_id)%id2(i,j,k,m) = ijks(3)
                        end do
                        patched(bc_id)%id2(i,j,k,nghp+1) = k
                    end if
        
                end do
                end do

            end do
        else if (s_nd == 2) then
        
            do j =s_st(2),s_ed(2)
                do k=s_st(3),s_ed(3)
                do i=s_st(1),s_ed(1)
                    if(scal(bc_id,1) >= 1.0) then
                        if(i<s_st(1)+nghp/2+1) then
                            il      = -1 + s_st(1) - i
                            P1      = -il
                        else if(i>s_ed(1)-nghp/2-1) then
                            il      = -nghp + s_ed(1) - i
                            P1      = -il
                        else
                            il      = -nghp/2-1
                            P1      = -il
                        end if

                        do m=1,nghp
                            patched(bc_id)%id1(i,j,k,m) = i+il+m
                        end do
                        patched(bc_id)%id1(i,j,k,nghp+1) = i
                    else
                        ijks(:) = (/i,j,k/)
                        ijkt(:) = mb_top(nb)%bcs(nr)%mapijk(i,j,k,:)
                        
                        it = ijkt(t_dir(2))
                        if(it<t_st(t_dir(2))+nghp/2+1) then
                            ilt      = -1 + t_st(t_dir(2)) - it
                            P1t      = -ilt
                        else if(it>t_ed(t_dir(2))-nghp/2-1) then
                            ilt      = -nghp + t_ed(t_dir(2)) - it
                            P1t      = -ilt
                        else
                            ilt      = -nghp/2-1
                            P1t      = -ilt
                        end if

                        do m=1,nghp
                            ijkt(t_dir(2)) = it+ilt+m
                            ijks(:)        = mb_top(nbt)%bcs(nrt)%mapijk(ijkt(1),ijkt(2),ijkt(3),:)

                            patched(bc_id)%id1(i,j,k,m) = ijks(1)
                        end do
                        patched(bc_id)%id1(i,j,k,nghp+1) = i
                    end if
        
                    if(scal(bc_id,3) >= 1.0) then
                        if(k<s_st(3)+nghp/2+1) then
                            kl      = -1 + s_st(3) - k
                            P2      = -kl
                        else if(k>s_ed(3)-nghp/2-1) then
                            kl      = -nghp + s_ed(3) - k
                            P2      = -kl
                        else
                            kl      = -nghp/2-1
                            P2      = -kl
                        end if

                        do m=1,nghp
                            patched(bc_id)%id2(i,j,k,m) = k+kl+m
                        end do
                        patched(bc_id)%id2(i,j,k,nghp+1) = k
                    else
                        ijks(:) = (/i,j,k/)
                        ijkt(:) = mb_top(nb)%bcs(nr)%mapijk(i,j,k,:)

                        kt = ijkt(t_dir(3))
                        if(kt<t_st(t_dir(3))+nghp/2+1) then
                            klt      = -1 + t_st(t_dir(3)) - kt
                            P2t      = -klt
                        else if(kt>t_ed(t_dir(3))-nghp/2-1) then
                            klt      = -nghp + t_ed(t_dir(3)) - kt
                            P2t      = -klt
                        else
                            klt      = -nghp/2-1
                            P2t      = -klt
                        end if

                        do m=1,nghp
                            ijkt(t_dir(3)) = kt+klt+m
                            ijks(:)        = mb_top(nbt)%bcs(nrt)%mapijk(ijkt(1),ijkt(2),ijkt(3),:)

                            patched(bc_id)%id2(i,j,k,m) = ijks(3)
                        end do
                        patched(bc_id)%id2(i,j,k,nghp+1) = k
                    end if
        
                end do
                end do

            end do
        
        else if (s_nd == 3) then

            do k =s_st(3),s_ed(3)
                do j=s_st(2),s_ed(2)
                do i=s_st(1),s_ed(1)
                    if(scal(bc_id,1) >= 1.0) then
                        if(i<s_st(1)+nghp/2+1) then
                            il      = -1 + s_st(1) - i
                            P1      = -il
                        else if(i>s_ed(1)-nghp/2-1) then
                            il      = -nghp + s_ed(1) - i
                            P1      = -il
                        else
                            il      = -nghp/2-1
                            P1      = -il
                        end if

                        do m=1,nghp
                            patched(bc_id)%id1(i,j,k,m) = i+il+m
                        end do
                        patched(bc_id)%id1(i,j,k,nghp+1) = i
                    else
                        ijks(:) = (/i,j,k/)
                        ijkt(:) = mb_top(nb)%bcs(nr)%mapijk(i,j,k,:)
                        
                        it = ijkt(t_dir(2))
                        if(it<t_st(t_dir(2))+nghp/2+1) then
                            ilt      = -1 + t_st(t_dir(2)) - it
                            P1t      = -ilt
                        else if(it>t_ed(t_dir(2))-nghp/2-1) then
                            ilt      = -nghp + t_ed(t_dir(2)) - it
                            P1t      = -ilt
                        else
                            ilt      = -nghp/2-1
                            P1t      = -ilt
                        end if

                        do m=1,nghp
                            ijkt(t_dir(2)) = it+ilt+m
                            ijks(:)        = mb_top(nbt)%bcs(nrt)%mapijk(ijkt(1),ijkt(2),ijkt(3),:)

                            patched(bc_id)%id1(i,j,k,m) = ijks(1)
                        end do
                        patched(bc_id)%id1(i,j,k,nghp+1) = i
                    end if
        
                    if(scal(bc_id,2) >= 1.0) then
                        if(j<s_st(2)+nghp/2+1) then
                            jl      = -1 + s_st(2) - j
                            P2      = -jl
                        else if(j>s_ed(2)-nghp/2-1) then
                            jl      = -nghp + s_ed(2) - j
                            P2      = -jl
                        else
                            jl      = -nghp/2-1
                            P2      = -jl
                        end if

                        do m=1,nghp
                            patched(bc_id)%id2(i,j,k,m) = j+jl+m
                        end do
                        patched(bc_id)%id2(i,j,k,nghp+1) = j
                    else
                        ijks(:) = (/i,j,k/)
                        ijkt(:) = mb_top(nb)%bcs(nr)%mapijk(i,j,k,:)

                        jt = ijkt(t_dir(3))
                        if(jt<t_st(t_dir(3))+nghp/2+1) then
                            jlt      = -1 + t_st(t_dir(3)) - jt
                            P2t      = -jlt
                        else if(jt>t_ed(t_dir(3))-nghp/2-1) then
                            jlt      = -nghp + t_ed(t_dir(3)) - jt
                            P2t      = -jlt
                        else
                            jlt      = -nghp/2-1
                            P2t      = -jlt
                        end if

                        do m=1,nghp
                            ijkt(t_dir(3)) = jt+jlt+m
                            ijks(:)        = mb_top(nbt)%bcs(nrt)%mapijk(ijkt(1),ijkt(2),ijkt(3),:)

                            patched(bc_id)%id2(i,j,k,m) = ijks(2)
                        end do
                        patched(bc_id)%id2(i,j,k,nghp+1) = j
                    end if
        
                end do
                end do

            end do
        
        end if
#ifdef PARALLEL
        end if
#endif
    end do

end subroutine calc_patched_indexes

subroutine extending_ghost_points(mb_var,nst,ned,ngh)
    use mod_kndconsts, only : kind_int
    use mod_datatypes, only : var_block_t
    use mod_fieldvars, only : nblkcoms,blkcoms,mb_top
    implicit none
    type(var_block_t), pointer, intent(in) :: mb_var(:)
    integer(kind_int),          intent(in) :: nst,ned,ngh
    integer(kind_int)          :: nc,nb,st(1:2,1:3),ierr
    integer(kind_int)          :: ni,nj,nk
    integer(kind_int)          :: i,j,k,m,n,n1,n2,lr(2),offsets,offsett,offsetc
    integer(kind_int)          :: offsetit,offsetic,offsetjt,offsetjc,offsetkt,offsetkc


    lr(1) = -1
    lr(2) = 1
    do nc=1,nblkcoms
        nb  =  blkcoms(nc)%nb

        ni = mb_top(nb)%nijk(1)
        nj = mb_top(nb)%nijk(2)
        nk = mb_top(nb)%nijk(3)
        
        st(1,1) = 1
        st(2,1) = ni
        st(1,2) = 1
        st(2,2) = nj
        st(1,3) = 1
        st(2,3) = nk
                  
        do m=nst,ned
            do n=1,2
            do i=1,ngh
            offsetc = st(n,1)+(i-1)*lr(n)
            offsets = offsetc-lr(n)
            offsett = offsetc+lr(n)
            do k=1,nk
            do j=1,nj
        
                mb_var(nb)%fld(m)%r3d(offsett,j,k) = mb_var(nb)%fld(m)%r3d(offsetc,j,k)*2.0 -  &
                                                     mb_var(nb)%fld(m)%r3d(offsets,j,k)
        
            end do
            end do
            end do
            end do
        end do
        
        do m=nst,ned
            do n=1,2
            do j=1,ngh
            offsetc = st(n,2)+(j-1)*lr(n)
            offsets = offsetc-lr(n)
            offsett = offsetc+lr(n)
            do k=1,nk
            do i=1,ni
        
                mb_var(nb)%fld(m)%r3d(i,offsett,k) = mb_var(nb)%fld(m)%r3d(i,offsetc,k)*2.0 -  &
                                                     mb_var(nb)%fld(m)%r3d(i,offsets,k)
        
            end do
            end do
            end do
            end do
        end do
        
        
        do m=nst,ned
            do n=1,2
            do k=1,ngh
            offsetc = st(n,3)+(k-1)*lr(n)
            offsets = offsetc-lr(n)
            offsett = offsetc+lr(n)
            do i=1,ni
            do j=1,nj
        
                mb_var(nb)%fld(m)%r3d(i,j,offsett) = mb_var(nb)%fld(m)%r3d(i,j,offsetc)*2.0 -  &
                                                     mb_var(nb)%fld(m)%r3d(i,j,offsets)
        
            end do
            end do
            end do
            end do
        end do
!---------------------------------------------------------------

        do m=nst,ned
            do n=1,2
            do i=1,ngh
            offsetic = st(n,1)+(i-1)*lr(n)
            offsetit = offsetic+lr(n)
            do n1=1,2
            do k=1,ngh
            offsetc = st(n1,3) + (k-1)*lr(n1)
            offsets = offsetc-lr(n1)
            offsett = offsetc+lr(n1)
            do j=1,nj
        
                mb_var(nb)%fld(m)%r3d(offsetit,j,offsett) = mb_var(nb)%fld(m)%r3d(offsetit,j,offsetc)*2.0 -  &
                                                            mb_var(nb)%fld(m)%r3d(offsetit,j,offsets)
        
            end do
            end do
            end do
            end do
            end do
        end do

        do m=nst,ned
            do n=1,2
            do i=1,ngh
            offsetic = st(n,1)+(i-1)*lr(n)
            offsetit = offsetic+lr(n)
            do n1=1,2
            do j=1,ngh
            offsetc = st(n1,2) + (j-1)*lr(n1)
            offsets = offsetc-lr(n1)
            offsett = offsetc+lr(n1)
            do k=1,nk
        
                mb_var(nb)%fld(m)%r3d(offsetit,offsett,k) = mb_var(nb)%fld(m)%r3d(offsetit,offsetc,k)*2.0 -  &
                                                            mb_var(nb)%fld(m)%r3d(offsetit,offsets,k)
        
            end do
            end do
            end do
            end do
            end do
        end do
        
        do m=nst,ned
            do n=1,2
            do j=1,ngh
            offsetjc = st(n,2)+(j-1)*lr(n)
            offsetjt = offsetjc+lr(n)
            do n1=1,2
            do k=1,ngh
            offsetc = st(n1,3)+(k-1)*lr(n1)
            offsets = offsetc-lr(n1)
            offsett = offsetc+lr(n1)
            do i=1,ni
        
                mb_var(nb)%fld(m)%r3d(i,offsetjt,offsett) = mb_var(nb)%fld(m)%r3d(i,offsetjt,offsetc)*2.0 -  &
                                                            mb_var(nb)%fld(m)%r3d(i,offsetjt,offsets)
        
            end do
            end do
            end do
            end do
            end do
        end do

        do m=nst,ned
            do n=1,2
            do j=1,ngh
            offsetjc = st(n,2)+(j-1)*lr(n)
            offsetjt = offsetjc+lr(n)
            do n1=1,2
            do i=1,ngh
            offsetc = st(n1,1)+(i-1)*lr(n1)
            offsets = offsetc-lr(n1)
            offsett = offsetc+lr(n1)
            do k=1,nk
        
                mb_var(nb)%fld(m)%r3d(offsett,offsetjt,k) = mb_var(nb)%fld(m)%r3d(offsetc,offsetjt,k)*2.0 -  &
                                                            mb_var(nb)%fld(m)%r3d(offsets,offsetjt,k)
        
            end do
            end do
            end do
            end do
            end do
        end do        
        
        do m=nst,ned
            do n=1,2
            do k=1,ngh
            offsetkc = st(n,3) + (k-1)*lr(n)
            offsetkt = offsetkc+lr(n)
            do n1=1,2
            do i=1,ngh
            offsetc = st(n1,1)+(i-1)*lr(n1)
            offsets = offsetc-lr(n1)
            offsett = offsetc+lr(n1)
            do j=1,nj
        
                mb_var(nb)%fld(m)%r3d(offsett,j,offsetkt) = mb_var(nb)%fld(m)%r3d(offsetc,j,offsetkt)*2.0 -  &
                                                            mb_var(nb)%fld(m)%r3d(offsets,j,offsetkt)
        
            end do
            end do
            end do
            end do
            end do
        end do

        do m=nst,ned
            do n=1,2
            do k=1,ngh
            offsetkc = st(n,3) + (k-1)*lr(n)
            offsetkt = offsetkc+lr(n)
            do n1=1,2
            do j=1,ngh
            offsetc = st(n1,2)+(j-1)*lr(n1)
            offsets = offsetc-lr(n1)
            offsett = offsetc+lr(n1)
            do i=1,ni
        
                mb_var(nb)%fld(m)%r3d(i,offsett,offsetkt) = mb_var(nb)%fld(m)%r3d(i,offsetc,offsetkt)*2.0 -  &
                                                            mb_var(nb)%fld(m)%r3d(i,offsets,offsetkt)
        
            end do
            end do
            end do
            end do
            end do
        end do
!---------------------------------------------------------------

        do m=nst,ned
            do n=1,2
            do i=1,ngh
            offsetic = st(n,1)+(i-1)*lr(n)
            offsetit = offsetic+lr(n)
            do n1=1,2
            do k=1,ngh
            offsetkc = st(n1,3) + (k-1)*lr(n1)
            offsetkt = offsetkc+lr(n1)
            do n2=1,2
            do j=1,ngh
            offsetc = st(n2,2)+(j-1)*lr(n2)
            offsets = offsetc-lr(n2)
            offsett = offsetc+lr(n2)
        
                mb_var(nb)%fld(m)%r3d(offsetit,offsett,offsetkt) = mb_var(nb)%fld(m)%r3d(offsetit,offsetc,offsetkt)*2.0 -  &
                                                                   mb_var(nb)%fld(m)%r3d(offsetit,offsets,offsetkt)
        
            end do
            end do
            end do
            end do
            end do
            end do
        end do

        do m=nst,ned
            do n=1,2
            do i=1,ngh
            offsetic = st(n,1)+(i-1)*lr(n)
            offsetit = offsetic+lr(n)
            do n1=1,2
            do j=1,ngh
            offsetjc = st(n1,2) + (j-1)*lr(n1)
            offsetjt = offsetjc+lr(n1)
            do n2=1,2
            do k=1,ngh
            offsetc = st(n2,3)+(k-1)*lr(n2)
            offsets = offsetc-lr(n2)
            offsett = offsetc+lr(n2)
        
                mb_var(nb)%fld(m)%r3d(offsetit,offsetjt,offsett) = mb_var(nb)%fld(m)%r3d(offsetit,offsetjt,offsetc)*2.0 -  &
                                                                   mb_var(nb)%fld(m)%r3d(offsetit,offsetjt,offsets)
        
            end do
            end do
            end do
            end do
            end do
            end do
        end do
        
        do m=nst,ned
            do n=1,2
            do j=1,ngh
            offsetjc = st(n,2)+(j-1)*lr(n)
            offsetjt = offsetjc+lr(n)
            do n1=1,2
            do k=1,ngh
            offsetkc = st(n1,3)+(k-1)*lr(n1)
            offsetkt = offsetkc+lr(n1)
            do n2=1,2
            do i=1,ngh
            offsetc = st(n2,1)+(i-1)*lr(n2)
            offsets = offsetc-lr(n2)
            offsett = offsetc+lr(n2)
        
                mb_var(nb)%fld(m)%r3d(offsett,offsetjt,offsetkt) = mb_var(nb)%fld(m)%r3d(offsetc,offsetjt,offsetkt)*2.0 -  &
                                                                   mb_var(nb)%fld(m)%r3d(offsets,offsetjt,offsetkt)
        
            end do
            end do
            end do
            end do
            end do
            end do
        end do

        do m=nst,ned
            do n=1,2
            do j=1,ngh
            offsetjc = st(n,2)+(j-1)*lr(n)
            offsetjt = offsetjc+lr(n)
            do n1=1,2
            do i=1,ngh
            offsetic = st(n1,1)+(i-1)*lr(n1)
            offsetit = offsetic+lr(n1)
            do n2=1,2
            do k=1,ngh
            offsetc = st(n2,3)+(k-1)*lr(n2)
            offsets = offsetc-lr(n2)
            offsett = offsetc+lr(n2)
        
                mb_var(nb)%fld(m)%r3d(offsetit,offsetjt,offsett) = mb_var(nb)%fld(m)%r3d(offsetit,offsetjt,offsetc)*2.0 -  &
                                                                   mb_var(nb)%fld(m)%r3d(offsetit,offsetjt,offsets)
        
            end do
            end do
            end do
            end do
            end do
            end do
        end do        
        
        do m=nst,ned
            do n=1,2
            do k=1,ngh
            offsetkc = st(n,3) + (k-1)*lr(n)
            offsetkt = offsetkc+lr(n)
            do n1=1,2
            do i=1,ngh
            offsetic = st(n1,1)+(i-1)*lr(n1)
            offsetit = offsetic+lr(n1)
            do n2=1,2
            do j=1,ngh
            offsetc = st(n2,2)+(j-1)*lr(n2)
            offsets = offsetc-lr(n2)
            offsett = offsetc+lr(n2)
        
                mb_var(nb)%fld(m)%r3d(offsetit,offsett,offsetkt) = mb_var(nb)%fld(m)%r3d(offsetit,offsetc,offsetkt)*2.0 -  &
                                                                   mb_var(nb)%fld(m)%r3d(offsetit,offsets,offsetkt)
        
            end do
            end do
            end do
            end do
            end do
            end do
        end do

        do m=nst,ned
            do n=1,2
            do k=1,ngh
            offsetkc = st(n,3) + (k-1)*lr(n)
            offsetkt = offsetkc+lr(n)
            do n1=1,2
            do j=1,ngh
            offsetjc = st(n1,2)+(j-1)*lr(n1)
            offsetjt = offsetjc+lr(n1)
            do n2=1,2
            do i=1,ngh
            offsetc = st(n2,1)+(i-1)*lr(n2)
            offsets = offsetc-lr(n2)
            offsett = offsetc+lr(n2)
        
                mb_var(nb)%fld(m)%r3d(offsett,offsetjt,offsetkt) = mb_var(nb)%fld(m)%r3d(offsetc,offsetjt,offsetkt)*2.0 -  &
                                                                   mb_var(nb)%fld(m)%r3d(offsets,offsetjt,offsetkt)
        
            end do
            end do
            end do
            end do
            end do
            end do
        end do
!---------------------------------------------------------------

    end do

end subroutine extending_ghost_points

!!subroutine patched_ghost_points(mb_var,nst,ned,ngh)
!!    use mod_kndconsts, only : kind_int,kind_real
!!    use mod_datatypes, only : var_block_t
!!    use mod_fieldvars, only : mb_top
!!    use mod_fieldvars, only : ninters,ninterf
!!    use mod_fieldvars, only : bc_seq
!!    use mod_fieldvars, only : inters,patched
!!    use mod_parallels
!!    implicit none
!!    type(var_block_t), pointer, intent(in) :: mb_var(:)
!!    integer(kind_int),          intent(in) :: nst,ned,ngh
!!    integer(kind_int), parameter :: nghp = 5
!!    integer(kind_int)             :: bc_id,nb,nr
!!    integer(kind_int)             :: nbs,nrs
!!    integer(kind_int)             :: s_st(3),s_ed(3)
!!    integer(kind_int)             :: s_lr,s_nd
!!    integer(kind_int)             :: id_src
!!    integer(kind_int)             :: i,j,k,nint,m,n,r,s
!!    integer(kind_int)             :: offsets,P1,P2,indx,indy,indz
!!    real(kind_real)                :: var_L,Ix(1:nghp),Iy(1:nghp),Iz(1:nghp),Px,Py,Pz,F(1:nghp,1:nghp)
!!
!!    do nint=ninterf+1,ninters
!!        nbs      = inters(nint)%bc%nbs
!!        nrs      = inters(nint)%bc%nrs
!!        s_st(:)  = inters(nint)%bc%s_st(:)
!!        s_ed(:)  = inters(nint)%bc%s_ed(:)
!!        s_nd     = inters(nint)%bc%s_nd
!!        s_lr     = inters(nint)%bc%s_lr
!!
!!        bc_id   = nint - ninterf
!!        nb      = bc_seq(bc_id,1)
!!        nr      = bc_seq(bc_id,2)
!!
!!#ifdef PARALLEL
!!        id_src = mb_top(nbs)%pid
!!        if (myid == id_src) then
!!#endif
!!        if (s_nd == 1) then
!!
!!            do m=nst,ned
!!            do i =1,ngh
!!            
!!                offsets = s_st(1)+s_lr*i
!!
!!                do k=s_st(3),s_ed(3)
!!                do j=s_st(2),s_ed(2)
!!
!!                    P1 = patched(bc_id)%id1(s_st(1),j,k,nghp+1)
!!                    P2 = patched(bc_id)%id2(s_st(1),j,k,nghp+1)
!!
!!                    Py = P1*1.0
!!                    Pz = P2*1.0
!!                    do n =1,nghp
!!                        Iy(n) = patched(bc_id)%id1(s_st(1),j,k,n)*1.0
!!                        Iz(n) = patched(bc_id)%id2(s_st(1),j,k,n)*1.0
!!                    end do
!!                    
!!                    do r=1,nghp
!!                    do s=1,nghp
!!                        indy   = patched(bc_id)%id1(s_st(1),j,k,r)
!!                        indz   = patched(bc_id)%id2(s_st(1),j,k,s)
!!                        F(r,s) = mb_var(nb)%fld(m)%r3d(offsets,indy,indz)
!!                    end do
!!                    end do
!!                    
!!                    call lagrange2(Iy,Iz,F,nghp,nghp,Py,Pz,var_L)
!!                    mb_var(nb)%fld(m)%r3d(offsets,j,k) = var_L
!!        
!!                end do
!!                end do
!!
!!            end do
!!            end do
!!        else if (s_nd == 2) then
!!        
!!            do m=nst,ned
!!            do j = 1,ngh
!!            
!!                offsets = s_st(2)+s_lr*j
!!
!!                do k=s_st(3),s_ed(3)
!!                do i=s_st(1),s_ed(1)
!!
!!                    P1 = patched(bc_id)%id1(i,s_st(2),k,nghp+1)
!!                    P2 = patched(bc_id)%id2(i,s_st(2),k,nghp+1)
!!
!!                    Px = P1*1.0
!!                    Pz = P2*1.0
!!                    do n =1,nghp
!!                        Ix(n) = patched(bc_id)%id1(i,s_st(2),k,n)*1.0
!!                        Iz(n) = patched(bc_id)%id2(i,s_st(2),k,n)*1.0
!!                    end do
!!                    
!!                    do r=1,nghp
!!                    do s=1,nghp
!!                        indx   = patched(bc_id)%id1(i,s_st(2),k,r)
!!                        indz   = patched(bc_id)%id2(i,s_st(2),k,s)
!!                        F(r,s) = mb_var(nb)%fld(m)%r3d(indx,offsets,indz)
!!                    end do
!!                    end do
!!                    
!!                    call lagrange2(Ix,Iz,F,nghp,nghp,Px,Pz,var_L)
!!                    mb_var(nb)%fld(m)%r3d(i,offsets,k) = var_L
!!        
!!                end do
!!                end do
!!
!!            end do
!!            end do
!!        
!!        else if (s_nd == 3) then
!!
!!            do m=nst,ned
!!            do k =1,ngh
!!            
!!                offsets = s_st(3)+s_lr*k
!!
!!                do j=s_st(2),s_ed(2)
!!                do i=s_st(1),s_ed(1)
!!
!!                    P1 = patched(bc_id)%id1(i,j,s_st(3),nghp+1)
!!                    P2 = patched(bc_id)%id2(i,j,s_st(3),nghp+1)
!!
!!                    Px = P1*1.0
!!                    Py = P2*1.0
!!                    do n =1,nghp
!!                        Ix(n) = patched(bc_id)%id1(i,j,s_st(3),n)*1.0
!!                        Iy(n) = patched(bc_id)%id2(i,j,s_st(3),n)*1.0
!!                    end do
!!                    
!!                    do r=1,nghp
!!                    do s=1,nghp
!!                        indx   = patched(bc_id)%id1(i,j,s_st(3),r)
!!                        indy   = patched(bc_id)%id2(i,j,s_st(3),s)
!!                        F(r,s) = mb_var(nb)%fld(m)%r3d(indx,indy,offsets)
!!                    end do
!!                    end do
!!                    
!!                    call lagrange2(Ix,Iy,F,nghp,nghp,Px,Py,var_L)
!!                    mb_var(nb)%fld(m)%r3d(i,j,offsets) = var_L
!!        
!!                end do
!!                end do
!!
!!            end do
!!            end do
!!        
!!        end if
!!#ifdef PARALLEL
!!        end if
!!#endif
!!    end do
!!
!!end subroutine patched_ghost_points

subroutine patched_ghost_points(mb_var,nst,ned,ngh,nswave)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : m3x3,one,nsgl_aver_art    
    use mod_datatypes, only : var_block_t,fld_array_t
    use mod_fieldvars, only : mb_top
    use mod_fieldvars, only : ninters,ninterf
    use mod_fieldvars, only : bc_seq,mb_sxyz
    use mod_fieldvars, only : inters,patched,mb_xyz
    use mod_parallels
    implicit none
    type(var_block_t), pointer, intent(inout) :: mb_var(:)
    integer(kind_int),          intent(in)    :: nst,ned,ngh,nswave
    integer(kind_int), parameter :: nghp = 5
    integer(kind_int)             :: bc_id,nb,nr
    integer(kind_int)             :: nbs,nrs
    integer(kind_int)             :: s_st(3),s_ed(3)
    integer(kind_int)             :: s_lr,s_nd
    integer(kind_int)             :: id_src
    integer(kind_int)             :: i,j,k,nint,m,n,r,s
    integer(kind_int)             :: m11,m12,m13,m21,m22,m23    
    integer(kind_int)             :: offsets,P1,P2,indx,indy,indz
    real(kind_real)                :: var_L,Ix(1:nghp),Iy(1:nghp),Iz(1:nghp),Px,Py,Pz,F(1:nghp,1:nghp)
    type(fld_array_t), pointer    :: sxyz(:)    

    do nint=ninterf+1,ninters
        nbs      = inters(nint)%bc%nbs
        nrs      = inters(nint)%bc%nrs
        s_st(:)  = inters(nint)%bc%s_st(:)
        s_ed(:)  = inters(nint)%bc%s_ed(:)
        s_nd     = inters(nint)%bc%s_nd
        s_lr     = inters(nint)%bc%s_lr

        bc_id   = nint - ninterf
        nb      = bc_seq(bc_id,1)
        nr      = bc_seq(bc_id,2)

#ifdef PARALLEL
        id_src = mb_top(nbs)%pid
        if (myid == id_src) then
#endif

        sxyz => mb_sxyz(nbs)%fld        
        
        if (s_nd == 1) then
            
            m11 = m3x3(1,2)
            m12 = m3x3(2,2)
            m13 = m3x3(3,2)
            m21 = m3x3(1,3)
            m22 = m3x3(2,3)
            m23 = m3x3(3,3)           

            do m=nst,ned
            do i =1,ngh
            
                offsets = s_st(1)+s_lr*i

                do k=s_st(3),s_ed(3)
                do j=s_st(2),s_ed(2)

                    P1 = patched(bc_id)%id1(s_st(1),j,k,nghp+1)
                    P2 = patched(bc_id)%id2(s_st(1),j,k,nghp+1)

                    Py = P1*1.0 !mb_xyz(nb)%fld(1)%r3d(offsets,P1,k)*sxyz(m11)%r3d(offsets,P1,k) + &
                                !mb_xyz(nb)%fld(2)%r3d(offsets,P1,k)*sxyz(m12)%r3d(offsets,P1,k) + &
                                !mb_xyz(nb)%fld(3)%r3d(offsets,P1,k)*sxyz(m13)%r3d(offsets,P1,k)
                    Pz = P2*1.0 !mb_xyz(nb)%fld(1)%r3d(offsets,j,P2)*sxyz(m21)%r3d(offsets,j,P2) + &
                                !mb_xyz(nb)%fld(2)%r3d(offsets,j,P2)*sxyz(m22)%r3d(offsets,j,P2) + &
                                !mb_xyz(nb)%fld(3)%r3d(offsets,j,P2)*sxyz(m23)%r3d(offsets,j,P2)
                    do n =1,nghp
                        indy  = patched(bc_id)%id1(s_st(1),j,k,n)
                        indz  = patched(bc_id)%id2(s_st(1),j,k,n)                       
                        Iy(n) = indy*1.0 !patched(bc_id)%datt(offsets,indy,k,1)*sxyz(m11)%r3d(offsets,P1,k) + &
                                         !patched(bc_id)%datt(offsets,indy,k,2)*sxyz(m12)%r3d(offsets,P1,k) + &
                                         !patched(bc_id)%datt(offsets,indy,k,3)*sxyz(m13)%r3d(offsets,P1,k)
                        Iz(n) = indz*1.0 !patched(bc_id)%datt(offsets,j,indz,1)*sxyz(m21)%r3d(offsets,j,P2) + &
                                         !patched(bc_id)%datt(offsets,j,indz,2)*sxyz(m22)%r3d(offsets,j,P2) + &
                                         !patched(bc_id)%datt(offsets,j,indz,3)*sxyz(m23)%r3d(offsets,j,P2)
                    end do
                    
                    do r=1,nghp
                    do s=1,nghp
                        indy   = patched(bc_id)%id1(s_st(1),j,k,r)
                        indz   = patched(bc_id)%id2(s_st(1),j,k,s)                     
                        F(r,s) = mb_var(nb)%fld(m)%r3d(offsets,indy,indz)
                    end do
                    end do
                    
                    call lagrange2(Iy,Iz,F,nghp,nghp,Py,Pz,var_L)
                    mb_var(nb)%fld(m)%r3d(offsets,j,k) = var_L
        
                end do
                end do

            end do
            end do
        else if (s_nd == 2) then
            
            m11 = m3x3(1,1)
            m12 = m3x3(2,1)
            m13 = m3x3(3,1)
            m21 = m3x3(1,3)
            m22 = m3x3(2,3)
            m23 = m3x3(3,3)            
        
            do m=nst,ned
            do j = 1,ngh
            
                offsets = s_st(2)+s_lr*j

                do k=s_st(3),s_ed(3)
                do i=s_st(1),s_ed(1)

                    P1 = patched(bc_id)%id1(i,s_st(2),k,nghp+1)
                    P2 = patched(bc_id)%id2(i,s_st(2),k,nghp+1)

                    Px = P1*1.0 !mb_xyz(nb)%fld(1)%r3d(P1,offsets,k)*sxyz(m11)%r3d(P1,offsets,k) + &
                                !mb_xyz(nb)%fld(2)%r3d(P1,offsets,k)*sxyz(m12)%r3d(P1,offsets,k) + &
                                !mb_xyz(nb)%fld(3)%r3d(P1,offsets,k)*sxyz(m13)%r3d(P1,offsets,k)
                    Pz = P2*1.0 !mb_xyz(nb)%fld(1)%r3d(i,offsets,P2)*sxyz(m21)%r3d(i,offsets,P2) + &
                                !mb_xyz(nb)%fld(2)%r3d(i,offsets,P2)*sxyz(m22)%r3d(i,offsets,P2) + &
                                !mb_xyz(nb)%fld(3)%r3d(i,offsets,P2)*sxyz(m23)%r3d(i,offsets,P2)
                    do n =1,nghp
                        indx  = patched(bc_id)%id1(i,s_st(2),k,n)
                        indz  = patched(bc_id)%id2(i,s_st(2),k,n)                        
                        Ix(n) = indx*1.0 !patched(bc_id)%datt(indx,offsets,k,1)*sxyz(m11)%r3d(P1,offsets,k) + &
                                         !patched(bc_id)%datt(indx,offsets,k,2)*sxyz(m12)%r3d(P1,offsets,k) + &
                                         !patched(bc_id)%datt(indx,offsets,k,3)*sxyz(m13)%r3d(P1,offsets,k)
                        Iz(n) = indz*1.0 !patched(bc_id)%datt(i,offsets,indz,1)*sxyz(m21)%r3d(i,offsets,P2) + &
                                         !patched(bc_id)%datt(i,offsets,indz,2)*sxyz(m22)%r3d(i,offsets,P2) + &
                                         !patched(bc_id)%datt(i,offsets,indz,3)*sxyz(m23)%r3d(i,offsets,P2)
                    end do
                    
                    do r=1,nghp
                    do s=1,nghp
                        indx   = patched(bc_id)%id1(i,s_st(2),k,r)
                        indz   = patched(bc_id)%id2(i,s_st(2),k,s)                        
                        F(r,s) = mb_var(nb)%fld(m)%r3d(indx,offsets,indz)
                    end do
                    end do
                    
                    call lagrange2(Ix,Iz,F,nghp,nghp,Px,Pz,var_L)                   
                    mb_var(nb)%fld(m)%r3d(i,offsets,k) = var_L
        
                end do
                end do

            end do
            end do
        
        else if (s_nd == 3) then
            
            m11 = m3x3(1,1)
            m12 = m3x3(2,1)
            m13 = m3x3(3,1)
            m21 = m3x3(1,2)
            m22 = m3x3(2,2)
            m23 = m3x3(3,2)            

            do m=nst,ned
            do k =1,ngh
            
                offsets = s_st(3)+s_lr*k

                do j=s_st(2),s_ed(2)
                do i=s_st(1),s_ed(1)

                    P1 = patched(bc_id)%id1(i,j,s_st(3),nghp+1)
                    P2 = patched(bc_id)%id2(i,j,s_st(3),nghp+1)

                    Px = P1*1.0 !mb_xyz(nb)%fld(1)%r3d(P1,j,offsets)*sxyz(m11)%r3d(P1,j,offsets) + &
                                !mb_xyz(nb)%fld(2)%r3d(P1,j,offsets)*sxyz(m12)%r3d(P1,j,offsets) + &
                                !mb_xyz(nb)%fld(3)%r3d(P1,j,offsets)*sxyz(m13)%r3d(P1,j,offsets)
                    Py = P2*1.0 !mb_xyz(nb)%fld(1)%r3d(i,P2,offsets)*sxyz(m21)%r3d(i,P2,offsets) + &
                                !mb_xyz(nb)%fld(2)%r3d(i,P2,offsets)*sxyz(m22)%r3d(i,P2,offsets) + &
                                !mb_xyz(nb)%fld(3)%r3d(i,P2,offsets)*sxyz(m23)%r3d(i,P2,offsets)
                    do n =1,nghp
                        indx  = patched(bc_id)%id1(i,j,s_st(3),n)
                        indy  = patched(bc_id)%id2(i,j,s_st(3),n)                        
                        Ix(n) = indx*1.0 !patched(bc_id)%datt(indx,j,offsets,1)*sxyz(m11)%r3d(P1,j,offsets) + &
                                         !patched(bc_id)%datt(indx,j,offsets,2)*sxyz(m12)%r3d(P1,j,offsets) + &
                                         !patched(bc_id)%datt(indx,j,offsets,3)*sxyz(m13)%r3d(P1,j,offsets)
                        Iy(n) = indy*1.0 !patched(bc_id)%datt(i,indy,offsets,1)*sxyz(m21)%r3d(i,P2,offsets) + &
                                         !patched(bc_id)%datt(i,indy,offsets,2)*sxyz(m22)%r3d(i,P2,offsets) + &
                                         !patched(bc_id)%datt(i,indy,offsets,3)*sxyz(m23)%r3d(i,P2,offsets)
                    end do
                    
                    do r=1,nghp
                    do s=1,nghp
                        indx   = patched(bc_id)%id1(i,j,s_st(3),r)
                        indy   = patched(bc_id)%id2(i,j,s_st(3),s)                       
                        F(r,s) = mb_var(nb)%fld(m)%r3d(indx,indy,offsets)
                    end do
                    end do
                    
                    call lagrange2(Ix,Iy,F,nghp,nghp,Px,Py,var_L)                   
                    mb_var(nb)%fld(m)%r3d(i,j,offsets) = var_L
        
                end do
                end do

            end do
            end do
        
        end if
#ifdef PARALLEL
        end if
#endif
    end do

end subroutine patched_ghost_points

subroutine lagrange(xa,ya,n,x,y)
!n      ���ͱ���������������ڵ����
!xa     n��Ԫ�ص�һάʵ�������飬�������������Ա�����ֵ�ڵ�xi��i=1,��,n��
!ya     n��Ԫ�ص�һάʵ�������飬�������,��ź���ֵ��y1,��,yn��T
!x      ʵ�ͱ����������������ֵ�Ա���
!y      ʵ�ͱ������������������ֵ
!------------------------------------------------------------------------------
    use mod_kndconsts, only : kind_int,kind_real
    implicit none
    integer(kind_int), intent(in)  :: n
    integer(kind_int)               :: i,j
!    integer(kind_int), parameter  :: nmax=10
    real(kind_real)  , intent(in)  :: x,xa(n),ya(n)
    real(kind_real)  , intent(out) :: y
    real(kind_real)                 :: l(n),dy

    l(1)=1

    do j=2,n
        l(1)=l(1)*(x-xa(j))/(xa(1)-xa(j))            !����l1(x)
    end do

    do i=2,n-1
        l(i)=1
        do j=1,i-1
            l(i)=l(i)*(x-xa(j))/(xa(i)-xa(j)) 
        end do

        do j=i+1,n
            l(i)=l(i)*(x-xa(j))/(xa(i)-xa(j))        !����li(x),1<i<n
        end do
    end do

    l(n)=1
    do j=1,n-1
        l(n)=l(n)*(x-xa(j))/(xa(n)-xa(j))        !����ln(x)
    end do

    y=0
    do i=1,n
        y=y+l(i)*ya(i)                 !����y=���li(x)*ya(i)
    end do
end subroutine lagrange

subroutine lagrange2(x1a,x2a,ya,m,n,x1,x2,y)
!m      ���ͱ��������������x�Ա����ڵ����
!n      ���ͱ��������������y�Ա����ڵ����
!x1a     m��Ԫ�ص�һάʵ�������飬������������x�Ա�����ֵ�ڵ�xi��i=1,��,m��
!x2a     n��Ԫ�ص�һάʵ�������飬������������y�Ա�����ֵ�ڵ�yj��j=1,��,n��
!x1     ʵ�ͱ����������������ֵx�Ա���
!x2     ʵ�ͱ����������������ֵy�Ա���
!ya     m��n��Ԫ�صĶ�άʵ�������飬�������,��ţ�xi,yj����i=1,��,m��j=1,��,n������ֵ��y1,��,yn��T
!y      ʵ�ͱ�������������������ֵ���
!------------------------------------------------------------------------------------------------------
    use mod_kndconsts, only : kind_int,kind_real
    implicit none
    integer(kind_int), intent(in)  :: m,n
    integer(kind_int)               :: i,j
    integer(kind_int), parameter  :: nmax=100,mmax=100
    real(kind_real)  , intent(in)  :: x1,x2,x1a(m),x2a(n),ya(m,n)
    real(kind_real)  , intent(out) :: y
    real(kind_real)                 :: ym(mmax),yn(nmax)

    do i=1,m

        do j=1,n
            yn(j)=ya(i,j)     !��ÿһ��xi,�ԣ�yj,zij����Ϊ��ֵ�ڵ�
        end do

        call lagrange(x2a,yn,n,x2,ym(i))   !��ã�xi,y���ĺ���ֵui

    end do

    call lagrange(x1a,ym,m,x1,y)   !�ԣ�xi,ui����ֵ����x,y������ֵ
end subroutine lagrange2