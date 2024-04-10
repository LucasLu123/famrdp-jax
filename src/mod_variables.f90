
module mod_variables
    use mod_kndconsts, only : kind_int,kind_real,len_char_name,len_char_file
    implicit none

    character(len_char_file) :: parfile

    character(len_char_name) :: casname
    character(len_char_file) :: indir,outdir
    namelist /general/ casname,indir,outdir

    real(kind_real) :: moo,reno
    real(kind_real) :: alt,attack,ayaw
    real(kind_real) :: rinf,vinf,pinf,tinf
    namelist /inflow/ moo,reno,alt,attack,ayaw, &
                      rinf,vinf,pinf,tinf

    real(kind_real) :: reflen,reflgrd
    real(kind_real) :: freflen,frefsc,frefxc,frefyc,frefzc
    namelist /reference/ reflen,reflgrd, &
                         freflen,frefsc,frefxc,frefyc,frefzc

    integer(kind_int)        :: topsty,grdsty,pltsty
    character(len_char_file) :: topfile,grdfile
    character(len_char_file) :: resfile,fcefile
    character(len_char_file) :: solfile,pltfile
    namelist /filename/ topfile,topsty,grdfile,grdsty, &
                        resfile,fcefile,solfile, &
                        pltfile,pltsty

    integer(kind_int) :: nincst
    real(kind_real)   :: mincst
    real(kind_real)   :: alpcst,betcst
    real(kind_real)   :: pincst,tincst
    namelist /initconst/ nincst,mincst,alpcst,betcst, &
                         pincst,tincst

    integer(kind_int) :: nrestrt,nstepmx
    integer(kind_int) :: nressav,nsolsav
    integer(kind_int) :: nfcesav,npltsav
    namelist /control/ nrestrt,nstepmx,nressav, &
                       nfcesav,nsolsav,npltsav

    integer(kind_int) :: ndim,nvis
    real(kind_real)   :: twall
    namelist /physical/ ndim,nvis,twall

    integer(kind_int) :: nturlhs,ntursub
    integer(kind_int) :: ntursch,nturint
    integer(kind_int) :: nrddst,nrddsp
    integer(kind_int) :: ntke2s,nearsm
    real(kind_real)   :: kfstur,kmaxtur
    real(kind_real)   :: vfstur,vmaxtur
    namelist /turbulent/ nturlhs,ntursub,ntursch,nturint,nrddst, &
                         ntke2s,nearsm,kfstur,kmaxtur,vfstur,vmaxtur

    !todo 预处理相关
    integer(kind_int) :: nprec!预处理标识
    integer(kind_int) :: nscmp!SCM-P标识
    integer(kind_int) :: nlhs,nsubmax
    real(kind_real)   :: tolsub
    integer(kind_int) :: nunst
    real(kind_real)   :: cflst,cfled,cflfac
    real(kind_real)   :: cdtmax,relaxs,relaxp,dtau,dtime
    namelist /timestep/ nprec,nscmp,nlhs,nsubmax,tolsub, &
                        nunst,cflst,cfled,cflfac, &
                        cdtmax,relaxs,relaxp,dtau,dtime

    integer(kind_int) :: nscheme,nflux,nintnon,nlimit
    real(kind_real)   :: enfix,ckmuscl,cbmuscl,csrvis,csafeup
    real(kind_real)   :: rlimit(2),plimit(2),tlimit(2)
    namelist /method/ nscheme,nflux,nintnon,nlimit, &
                      enfix,ckmuscl,cbmuscl,csrvis,csafeup, &
                      rlimit,plimit,tlimit

    integer(kind_int) :: ncutpol,ncelcet
    integer(kind_int) :: nghnode,nghedge
    integer(kind_int) :: nsponge
    real(kind_real)   :: sigspg
    integer(kind_int) :: nacous,nprms,ns_mean,ns_prms
    namelist /technic/ ncutpol,ncelcet,nghnode,nghedge,nsponge,sigspg,nrddsp,nacous,nprms,ns_mean,ns_prms

    real(kind_real) :: wgas,gamma
    real(kind_real) :: mu0sth,t0sth,tssth
    real(kind_real) :: prlam,prtur,sclam,sctur
    namelist /gasmodel/ wgas,gamma,mu0sth,t0sth,tssth, &
                        prlam,prtur,sclam,sctur
    
    integer(kind_int) :: nmoni,nbmoni(30),nimoni(30),njmoni(30),nkmoni(30)
    namelist /monitor/ nmoni,nbmoni,nimoni,njmoni,nkmoni

    character(len_char_file) :: splfile

    real(kind_real) :: rgas,mgas

    real(kind_real) :: reue
    real(kind_real) :: refden,refvel,refpres,reftem
    real(kind_real) :: refengy,refcp,refbeta,refvis
    real(kind_real) :: reftime,refqw
    real(kind_real) :: tssnd
    real(kind_real) :: roo,uoo,voo,woo,poo
    real(kind_real) :: visoo

    real(kind_real) :: twlnd

    integer(kind_int) :: nsw_kdir
    real(kind_real)   :: fsw_kdir

    integer(kind_int) :: nd1der_con,nd1int_con
    integer(kind_int) :: nd2der_vis,nd2int_vis,nd3der_vis,nd3int_vis
    integer(kind_int) :: nd2der_grd,nd2int_grd,nd3der_grd,nd3int_grd

    real(kind_real) :: cfl,cflmin,cflmax
    real(kind_real) :: dtaumin,dtaumax
    real(kind_real) :: rcflmin,rcflmax

    integer(kind_int) :: nstep,nstepsav
    integer(kind_int) :: nstepst,nsteped

    integer(kind_int) :: respos(5)
    real(kind_real)   :: resmax,restot
    real(kind_real)   :: restp0,restpn,restpm
    real(kind_real)   :: resdts

    real(kind_real)   :: turconsts(32)
    integer(kind_int) :: nsw_dst,nsw_bld

    integer(kind_int) :: nd1der_tur_con,nd1int_tur_con
    integer(kind_int) :: nd2der_tur_vis,nd2int_tur_vis
    integer(kind_int) :: nd3der_tur_vis,nd3int_tur_vis

    real(kind_real) :: cfx,cfy,cfz,cmx,cmy,cmz

    real(kind_real) :: edvisoo
    real(kind_real) :: nuoo
    real(kind_real) :: tkeoo,omeoo

    real(kind_real) :: dspmin,dspmax
    
    integer(kind_int) :: niter,nsub,nsubp,nsubp1

end module mod_variables
