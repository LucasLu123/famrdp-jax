
subroutine init_inflow
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,one,deg2rad,runiv
    use mod_variables, only : moo,reno,nscmp,nprec,nflux,nlhs
    use mod_variables, only : alt,attack,ayaw
    use mod_variables, only : rinf,vinf,pinf,tinf
    use mod_variables, only : reflen,reflgrd,twall
    use mod_variables, only : wgas,gamma,mu0sth,t0sth,tssth
    use mod_variables, only : rgas,reftime,reue,tssnd,twlnd
    use mod_variables, only : refden,refvel,refpres,reftem
    use mod_variables, only : refengy,refcp,refbeta,refvis
    use mod_variables, only : refqw,roo,uoo,voo,woo,poo,visoo
    implicit none
    integer(kind_int)         :: n(3),m
    real(kind_real), external :: sutherland0
    real(kind_real)           :: ainf,muinf
    external                  :: inflow_check

    wgas = 1.0e-3*wgas
    rgas = runiv/wgas

    if(nscmp>0) then
        if(nprec>0 .or. nflux/=7 .or. (nlhs/=1 .and. nlhs/=12 .and. nlhs/=13)) then
            call error_check(1, "Parameters are invalid for SCM model")
        end if
    end if
    
    if(nprec>0) then
        if(nscmp>0 .or. (nflux/=5 .and. nflux/=6) .or. (nlhs/=1 .and. nlhs/=11)) then
            call error_check(1, "Parameters are invalid for Precondition model")
        end if
    end if
    
    if(nscmp<=0 .and. nprec<=0) then
        if(nflux==5 .or. nflux==6 .or. nflux==7 .or. nlhs==11 .or. nlhs==12 .or. nlhs==13) then
            call error_check(1, "Parameters are invalid for simulation model")
        end if
    end if
    
    if (alt < zero) then
        n(:) = 0
        if (pinf > zero) n(1) = 1
        if (rinf > zero) n(2) = 1
        if (tinf > zero) n(3) = 1
        m = sum(n(:))
        if (m == 2) then
            if (n(1) == 0) then
                pinf = rinf*rgas*tinf
            else if (n(2) == 0) then
                rinf = pinf/(rgas*tinf)
            else !if (n(3) == 0) then
                tinf = pinf/(rgas*rinf)
            end if
        else
            call error_check(1, "(pinf,rinf,tinf) is invalid")
        end if
    else
        call stdatm1976(alt,rinf,pinf,tinf)
    end if


    ainf = sqrt(gamma*rgas*tinf)
    if (vinf > zero) then
        moo = vinf/ainf
    else
        vinf = moo*ainf
    end if

    refden  = rinf
    refvel  = vinf
    refpres = rinf*vinf*vinf
    reftem  = tinf
    refengy = refpres/refden

    refcp   = refengy/reftem               ! Cp and R have the same dimension
    refbeta = rgas/refcp                   ! rgas/refcp=Roo(nondimensional), p=beta*ro*t

    refqw   = rinf*vinf**3

    tssnd = tssth/reftem

    muinf = sutherland0(tinf,mu0sth,t0sth,tssth)

    if (alt >= zero) then
        reno = rinf*vinf*reflen/muinf
    end if

    refvis = muinf

    reue = reno*reflgrd/reflen!粘性项前面的雷诺数

    reftime = reflgrd/refvel

    attack = attack*deg2rad
    ayaw   = ayaw  *deg2rad

    roo = one
    uoo = cos(attack)*cos(ayaw)
    voo = sin(attack)*cos(ayaw)
    woo = sin(ayaw)
    poo = pinf/refpres
    visoo = one

    twlnd = twall/reftem

    call run_seq_and_master(inflow_check)

end subroutine init_inflow

subroutine inflow_check
    use mod_kndconsts, only : kind_real
    use mod_constants, only : epsil
    use mod_variables, only : moo,rgas,refvis,refvel
    use mod_variables, only : reftem,refpres,refden
    use mod_variables, only : reno,reflen,reue,reflgrd
    implicit none
    real(kind_real) :: roc,psc,rec,lcc

    lcc = reno*refvis/(refden*refvel)
    roc = reno*refvis/(reflen*refvel)
    psc = roc*reftem*rgas
    rec = refden*refvel*reflen/refvis
    if (abs(lcc-reflen)/reflen > epsil) then
        call msg_prompt("Optional freestream and reference parameters")
        write(*,*) "  moo    =",moo
        write(*,*) "  reflen =",lcc
        write(*,*) "  rinf   =",roc
        write(*,*) "  pinf   =",psc
        write(*,*) "  reno   =",rec

        call error_check(1, "inflow is invalid")
    end if

    write(*,*) "========================================================="
    write(*,*) "         The freestream and reference parameters"
    write(*,*) "---------------------------------------------------------"
    write(*,'(3x,3a12  )') "Mach","Viscosity","Velocity"
    write(*,'(3x,3e12.5)') moo,refvis,refvel
    write(*,'(3x,4a12  )') "Temperature","Pressure","Density","Sound"
    write(*,'(3x,4e12.5)') reftem,refpres,refden,refvel/moo
    write(*,'(3x,4a12  )') "Reynolds","RefLength","RenoGrid","RefLGrid"
    write(*,'(3x,4e12.5)') reno,reflen,reue,reflgrd
    write(*,*) "========================================================="

end subroutine inflow_check


subroutine init_fld_variables
    use mod_constants, only : nvis_ns_lam
    use mod_variables, only : nvis
    implicit none

    if (nvis > nvis_ns_lam) then
        call init_tur_field
    end if
    call init_flow_field    

end subroutine init_fld_variables


subroutine init_flow_field
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,deg2rad,nvis_euler
    use mod_constants, only : nincst_close,nincst_inter,nincst_intbc
    use mod_constants, only : nlhs_prsgs_ust_mat,nlhs_prsgs_ust_sca,nlhs_prsgs_ust_sca_scmp
    use mod_variables, only : solfile,nrestrt,nincst,nghnode,nprms,nacous
    use mod_variables, only : mincst,alpcst,betcst,pincst,tincst
    use mod_variables, only : roo,uoo,voo,woo,poo,nvis,nlhs,nstepsav
    use mod_variables, only : refpres,reftem,refbeta,gamma,nsub,nsubp,nsubp1,nsubmax
    use mod_fieldvars, only : npvs,mb_pv,neqn,mb_qc
    use mod_fieldvars, only : mb_vst,mb_q0,ntime,npoints
    use mod_fieldvars, only : mb_frms,mb_fmean,nsteprms,nstepmean
    use mod_interface, only : assign_mb_var_uniform
    use mod_interface, only : calc_mb_var_via_sub
    implicit none
    real(kind_real) :: pvoo(1:npvs),gama,fsoo(1:npvs)
    real(kind_real) :: a,va,tm,ro,vx,vy,vz,ps
    external        :: prim2con,v1_eq_v2

    gama = gamma

    select case(nincst)
    case(:nincst_close)
        pvoo(:) = (/roo,uoo,voo,woo,poo/)
    case(nincst_inter,nincst_intbc)
        alpcst = alpcst*deg2rad
        betcst = betcst*deg2rad
        ps = pincst/refpres
        tm = tincst/reftem
        ro = ps/(refbeta*tm)

        a = sqrt(gama*ps/ro)
        va = mincst*a
        vx = va*cos(alpcst)*cos(betcst)
        vy = va*sin(alpcst)*cos(betcst)
        vz = va*sin(betcst)

        pvoo(:) = (/ro,vx,vy,vz,ps/)

        if (nincst == nincst_intbc) then
            roo = ro
            uoo = vx
            voo = vy
            woo = vz
            poo = ps
        end if
    end select

    if (nacous > 0) then
        ntime     = 0
        npoints   = 0
        call calculate_points_acoustic
    end if

    fsoo(:) = (/0.0,0.0,0.0,0.0,0.0/)
    if (nprms > 0) then
        nstepmean = 0
        nsteprms  = 0
        call assign_mb_var_uniform(mb_fmean,1,npvs,nghnode,fsoo)        
        call assign_mb_var_uniform(mb_frms ,1,npvs,nghnode,fsoo)
    end if

    call input_sol

    call calc_mb_var_via_sub(mb_pv,1,npvs,prim2con,mb_qc,1,neqn,nghnode)

    select case(nlhs)
    case(nlhs_prsgs_ust_mat,nlhs_prsgs_ust_sca,nlhs_prsgs_ust_sca_scmp)
        call calc_mb_var_via_sub(mb_qc,1,neqn,v1_eq_v2,mb_q0,1,neqn,nghnode)
        if (nrestrt == 0 .or. nstepsav == 0) then
            nsub   = nsubmax
            nsubp  = nsubmax
            nsubp1 = nsubmax
        end if
    case default
        nsub   = nsubmax
        nsubp  = nsubmax
        nsubp1 = nsubmax
    end select

    if (nvis > nvis_euler) then
        call assign_mb_var_uniform(mb_vst,1,1,nghnode,(/zero/))
    end if


end subroutine init_flow_field

subroutine init_fld_variables_sp
    use mod_constants, only : nvis_ns_lam
    use mod_variables, only : nvis
    implicit none

    !todo if (nvis > nvis_ns_lam) then
    !todo     call init_tur_field_sp
    !todo end if
    call init_flow_field_sp    

end subroutine init_fld_variables_sp


subroutine init_flow_field_sp
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,deg2rad,nvis_euler
    use mod_constants, only : nincst_close,nincst_inter,nincst_intbc
    use mod_constants, only : nlhs_prsgs_ust_mat,nlhs_prsgs_ust_sca,nlhs_prsgs_ust_sca_scmp
    use mod_variables, only : solfile,nrestrt,nincst,nghnode,nprms,nacous
    use mod_variables, only : mincst,alpcst,betcst,pincst,tincst
    use mod_variables, only : roo,uoo,voo,woo,poo,nvis,nlhs,nstepsav
    use mod_variables, only : refpres,reftem,refbeta,gamma,nsub,nsubp,nsubp1,nsubmax
    use mod_fieldvars, only : npvs,mb_pv,neqn,mb_qc
    use mod_fieldvars, only : mb_vst,mb_q0,ntime,npoints
    use mod_fieldvars, only : mb_frms,mb_fmean,nsteprms,nstepmean
    use mod_interface, only : assign_mb_var_uniform_sp
    use mod_interface, only : calc_mb_var_via_sub_sp
    implicit none
    real(kind_real) :: pvoo(1:npvs),gama,fsoo(1:npvs)
    real(kind_real) :: a,va,tm,ro,vx,vy,vz,ps
    external        :: prim2con,v1_eq_v2

    gama = gamma

    select case(nincst)
    case(:nincst_close)
        pvoo(:) = (/roo,uoo,voo,woo,poo/)
    case(nincst_inter,nincst_intbc)
        alpcst = alpcst*deg2rad
        betcst = betcst*deg2rad
        ps = pincst/refpres
        tm = tincst/reftem
        ro = ps/(refbeta*tm)

        a = sqrt(gama*ps/ro)
        va = mincst*a
        vx = va*cos(alpcst)*cos(betcst)
        vy = va*sin(alpcst)*cos(betcst)
        vz = va*sin(betcst)

        pvoo(:) = (/ro,vx,vy,vz,ps/)

        if (nincst == nincst_intbc) then
            roo = ro
            uoo = vx
            voo = vy
            woo = vz
            poo = ps
        end if
    end select

    if (nacous > 0) then
        ntime     = 0
        npoints   = 0
        call calculate_points_acoustic_sp
    end if

    fsoo(:) = (/0.0,0.0,0.0,0.0,0.0/)
    if (nprms > 0) then
        nstepmean = 0
        nsteprms  = 0
        call assign_mb_var_uniform_sp(mb_fmean,1,npvs,nghnode,fsoo)        
        call assign_mb_var_uniform_sp(mb_frms ,1,npvs,nghnode,fsoo)
    end if

    call input_sol_sp

    call calc_mb_var_via_sub_sp(mb_pv,1,npvs,prim2con,mb_qc,1,neqn,nghnode)

    select case(nlhs)
    case(nlhs_prsgs_ust_mat,nlhs_prsgs_ust_sca,nlhs_prsgs_ust_sca_scmp)
        call calc_mb_var_via_sub_sp(mb_qc,1,neqn,v1_eq_v2,mb_q0,1,neqn,nghnode)
        if (nrestrt == 0 .or. nstepsav == 0) then
            nsub   = nsubmax
            nsubp  = nsubmax
            nsubp1 = nsubmax
        end if
    case default
        nsub   = nsubmax
        nsubp  = nsubmax
        nsubp1 = nsubmax
    end select

    if (nvis > nvis_euler) then
        call assign_mb_var_uniform_sp(mb_vst,1,1,nghnode,(/zero/))
    end if


end subroutine init_flow_field_sp

subroutine Rest_flow_field
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,deg2rad,nvis_euler
    use mod_constants, only : nincst_close,nincst_inter,nincst_intbc
    use mod_constants, only : nlhs_prsgs_ust_mat,nlhs_prsgs_ust_sca
    use mod_variables, only : solfile,nincst,nrestrt,nghnode,nprms,nacous
    use mod_variables, only : mincst,alpcst,betcst,pincst,tincst
    use mod_variables, only : roo,uoo,voo,woo,poo,nvis,nlhs,nstepsav
    use mod_variables, only : refpres,reftem,refbeta,gamma,nsub,nsubp,nsubp1,nsubmax
    use mod_fieldvars, only : npvs,mb_pv,neqn,mb_qc
    use mod_fieldvars, only : mb_vst,mb_q0,ntime,npoints
    use mod_fieldvars, only : mb_frms,mb_fmean,nsteprms,nstepmean
    use mod_interface, only : assign_mb_var_uniform
    implicit none
    real(kind_real) :: pvoo(1:npvs),gama,fsoo(1:npvs)
    real(kind_real) :: a,va,tm,ro,vx,vy,vz,ps
    external        :: prim2con,v1_eq_v2

    gama    = gamma
    nrestrt = 0    

    select case(nincst)
    case(:nincst_close)
        pvoo(:) = (/roo,uoo,voo,woo,poo/)
    case(nincst_inter,nincst_intbc)
        alpcst = alpcst*deg2rad
        betcst = betcst*deg2rad
        ps = pincst/refpres
        tm = tincst/reftem
        ro = ps/(refbeta*tm)

        a = sqrt(gama*ps/ro)
        va = mincst*a
        vx = va*cos(alpcst)*cos(betcst)
        vy = va*sin(alpcst)*cos(betcst)
        vz = va*sin(betcst)

        pvoo(:) = (/ro,vx,vy,vz,ps/)

        if (nincst == nincst_intbc) then
            roo = ro
            uoo = vx
            voo = vy
            woo = vz
            poo = ps
        end if
    end select

    call assign_mb_var_uniform(mb_pv,1,npvs,nghnode,pvoo)
    !!call assign_mb_var_uniform1(mb_pv,1,npvs,nghnode,pvoo)

    if (nacous > 0) then
        ntime     = 0
        npoints   = 0
        call calculate_points_acoustic
    end if

    fsoo(:) = (/0.0,0.0,0.0,0.0,0.0/)
    if (nprms > 0) then
        nstepmean = 0
        nsteprms  = 0
        call assign_mb_var_uniform(mb_fmean,1,npvs,nghnode,fsoo)        
        call assign_mb_var_uniform(mb_frms ,1,npvs,nghnode,fsoo)
    end if

end subroutine Rest_flow_field

subroutine Rest_flow_field_sp
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,deg2rad,nvis_euler
    use mod_constants, only : nincst_close,nincst_inter,nincst_intbc
    use mod_constants, only : nlhs_prsgs_ust_mat,nlhs_prsgs_ust_sca
    use mod_variables, only : solfile,nincst,nrestrt,nghnode,nprms,nacous
    use mod_variables, only : mincst,alpcst,betcst,pincst,tincst
    use mod_variables, only : roo,uoo,voo,woo,poo,nvis,nlhs,nstepsav
    use mod_variables, only : refpres,reftem,refbeta,gamma,nsub,nsubp,nsubp1,nsubmax
    use mod_fieldvars, only : npvs,mb_pv,mb_pvfp,neqn,mb_qc
    use mod_fieldvars, only : mb_vst,mb_q0,ntime,npoints
    use mod_fieldvars, only : mb_frms,mb_fmean,nsteprms,nstepmean
    use mod_interface, only : assign_mb_var_uniform_sp,assign_mb_var_uniform
    implicit none
    real(kind_real) :: pvoo(1:npvs),gama,fsoo(1:npvs)
    real(kind_real) :: a,va,tm,ro,vx,vy,vz,ps
    external        :: prim2con,v1_eq_v2

    gama    = gamma
    nrestrt = 0    

    select case(nincst)
    case(:nincst_close)
        pvoo(:) = (/roo,uoo,voo,woo,poo/)
    case(nincst_inter,nincst_intbc)
        alpcst = alpcst*deg2rad
        betcst = betcst*deg2rad
        ps = pincst/refpres
        tm = tincst/reftem
        ro = ps/(refbeta*tm)

        a = sqrt(gama*ps/ro)
        va = mincst*a
        vx = va*cos(alpcst)*cos(betcst)
        vy = va*sin(alpcst)*cos(betcst)
        vz = va*sin(betcst)

        pvoo(:) = (/ro,vx,vy,vz,ps/)

        if (nincst == nincst_intbc) then
            roo = ro
            uoo = vx
            voo = vy
            woo = vz
            poo = ps
        end if
    end select

    call assign_mb_var_uniform_sp(mb_pv,1,npvs,nghnode,pvoo)
    call assign_mb_var_uniform(mb_pvfp,1,npvs,nghnode,pvoo)
    !!call assign_mb_var_uniform1(mb_pv,1,npvs,nghnode,pvoo)

    if (nacous > 0) then
        ntime     = 0
        npoints   = 0
        call calculate_points_acoustic_sp
    end if

    fsoo(:) = (/0.0,0.0,0.0,0.0,0.0/)
    if (nprms > 0) then
        nstepmean = 0
        nsteprms  = 0
        call assign_mb_var_uniform_sp(mb_fmean,1,npvs,nghnode,fsoo)        
        call assign_mb_var_uniform_sp(mb_frms ,1,npvs,nghnode,fsoo)
    end if

end subroutine Rest_flow_field_sp

subroutine init_tur_field
    use mod_kndconsts, only : kind_real
    use mod_constants, only : thr2nd,ten,nrestrt_cont_tur
    use mod_constants, only : nvis_tur_sa,nvis_tur_sst,nvis_tur_hst
    use mod_variables, only : nvis,vfstur,kfstur
    use mod_variables, only : roo,uoo,voo,woo,visoo
    use mod_variables, only : reflen,reflgrd,reue
    use mod_variables, only : nrestrt,nghnode
    use mod_variables, only : edvisoo,nuoo,tkeoo,omeoo
    use mod_fieldvars, only : mb_qt,neqt,mb_vst
    use mod_interface, only : assign_mb_var_uniform
    implicit none
    real(kind_real) :: v2oo
    real(kind_real) :: qtoo(1:neqt)

    edvisoo = vfstur*visoo
    select case(nvis)
    case(nvis_tur_sa)
       nuoo = kfstur*visoo/roo
       qtoo(:) = (/nuoo/)
    case(nvis_tur_sst,nvis_tur_hst)
       ! low(intensity=1%,viscosity ratio=1)
       ! medium(intensity=5%,viscosity ratio=10)
       ! high(intensity=10%,viscosity ratio=100)
       v2oo  = uoo*uoo+voo*voo+woo*woo
       tkeoo = thr2nd*kfstur*kfstur*v2oo
       omeoo = roo*tkeoo*reue/edvisoo

       !omeoo = ten/(reflen/reflgrd)
       !tkeoo = edvisoo*omeoo/(roo*reue)

       qtoo(:) = (/tkeoo,omeoo/)
    end select

    call assign_mb_var_uniform(mb_qt ,1,neqt,nghnode,qtoo)
    call assign_mb_var_uniform(mb_vst,1,   1,nghnode,(/edvisoo/))

    if (nrestrt == nrestrt_cont_tur) then
        call input_tur_sol
    end if

end subroutine init_tur_field

subroutine init_tur_field_sp
    use mod_kndconsts, only : kind_real
    use mod_constants, only : thr2nd,ten,nrestrt_cont_tur
    use mod_constants, only : nvis_tur_sa,nvis_tur_sst,nvis_tur_hst
    use mod_variables, only : nvis,vfstur,kfstur
    use mod_variables, only : roo,uoo,voo,woo,visoo
    use mod_variables, only : reflen,reflgrd,reue
    use mod_variables, only : nrestrt,nghnode
    use mod_variables, only : edvisoo,nuoo,tkeoo,omeoo
    use mod_fieldvars, only : mb_qt,neqt,mb_vst
    use mod_interface, only : assign_mb_var_uniform_sp
    implicit none
    real(kind_real) :: v2oo
    real(kind_real) :: qtoo(1:neqt)

    edvisoo = vfstur*visoo
    select case(nvis)
    case(nvis_tur_sa)
       nuoo = kfstur*visoo/roo
       qtoo(:) = (/nuoo/)
    case(nvis_tur_sst,nvis_tur_hst)
       ! low(intensity=1%,viscosity ratio=1)
       ! medium(intensity=5%,viscosity ratio=10)
       ! high(intensity=10%,viscosity ratio=100)
       v2oo  = uoo*uoo+voo*voo+woo*woo
       tkeoo = thr2nd*kfstur*kfstur*v2oo
       omeoo = roo*tkeoo*reue/edvisoo

       !omeoo = ten/(reflen/reflgrd)
       !tkeoo = edvisoo*omeoo/(roo*reue)

       qtoo(:) = (/tkeoo,omeoo/)
    end select

    call assign_mb_var_uniform_sp(mb_qt ,1,neqt,nghnode,qtoo)
    call assign_mb_var_uniform_sp(mb_vst,1,   1,nghnode,(/edvisoo/))

    if (nrestrt == nrestrt_cont_tur) then
        call input_tur_sol_sp
    end if

end subroutine init_tur_field_sp





