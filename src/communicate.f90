
subroutine broadcast_par
#ifdef PARALLEL
    use mod_kndconsts, only : kind_long,kind_int
    use mod_variables
    use mod_parallels
    implicit none
    integer(kind_int),parameter :: packsize=10000
    integer(kind_long)          :: packbuf(packsize)
    integer(kind_int)           :: position,ierr

    if (myid == master) then
        position = 0

        call MPI_PACK(casname,len_char_name,kind_char_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(indir  ,len_char_file,kind_char_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(outdir ,len_char_file,kind_char_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)

        call MPI_PACK(moo   ,1,kind_real_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(reno  ,1,kind_real_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(alt   ,1,kind_real_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(attack,1,kind_real_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(ayaw  ,1,kind_real_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(rinf  ,1,kind_real_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(vinf  ,1,kind_real_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(pinf  ,1,kind_real_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(tinf  ,1,kind_real_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)

        call MPI_PACK(reflen ,1,kind_real_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(reflgrd,1,kind_real_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(freflen,1,kind_real_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(frefsc ,1,kind_real_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(frefxc ,1,kind_real_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(frefyc ,1,kind_real_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(frefzc ,1,kind_real_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)

        call MPI_PACK(topsty ,1,kind_int_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(grdsty ,1,kind_int_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(pltsty ,1,kind_int_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(topfile,len_char_file,kind_char_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(grdfile,len_char_file,kind_char_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(resfile,len_char_file,kind_char_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(fcefile,len_char_file,kind_char_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(solfile,len_char_file,kind_char_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(pltfile,len_char_file,kind_char_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)

        call MPI_PACK(nincst,1,kind_int_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(mincst,1,kind_real_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(alpcst,1,kind_real_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(betcst,1,kind_real_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(pincst,1,kind_real_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(tincst,1,kind_real_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)

        call MPI_PACK(nrestrt,1,kind_int_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(nstepmx,1,kind_int_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(nressav,1,kind_int_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(nsolsav,1,kind_int_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(nfcesav,1,kind_int_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(npltsav,1,kind_int_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)

        call MPI_PACK(ndim ,1,kind_int_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(nvis ,1,kind_int_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(twall,1,kind_real_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)

        call MPI_PACK(nturlhs,1,kind_int_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(ntursub,1,kind_int_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(ntursch,1,kind_int_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(nturint,1,kind_int_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(nrddst ,1,kind_int_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(ntke2s ,1,kind_int_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(nearsm ,1,kind_int_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(kfstur ,1,kind_real_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(kmaxtur,1,kind_real_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(vfstur ,1,kind_real_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(vmaxtur,1,kind_real_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)

        call MPI_PACK(nprec  ,1,kind_int_mpi ,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(nscmp  ,1,kind_int_mpi ,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(nlhs   ,1,kind_int_mpi ,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(nsubmax,1,kind_int_mpi ,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(tolsub ,1,kind_real_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(nunst  ,1,kind_int_mpi ,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(cflst  ,1,kind_real_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(cfled  ,1,kind_real_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(cflfac ,1,kind_real_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(cdtmax ,1,kind_real_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(relaxs ,1,kind_real_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(relaxp ,1,kind_real_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(dtau   ,1,kind_real_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(dtime  ,1,kind_real_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)

        call MPI_PACK(nscheme,1,kind_int_mpi ,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(nflux  ,1,kind_int_mpi ,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(nintnon,1,kind_int_mpi ,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(nlimit ,1,kind_int_mpi ,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(enfix  ,1,kind_real_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(ckmuscl,1,kind_real_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(cbmuscl,1,kind_real_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(csrvis ,1,kind_real_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(csafeup,1,kind_real_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(rlimit ,2,kind_real_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(plimit ,2,kind_real_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(tlimit ,2,kind_real_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)

        call MPI_PACK(ncutpol,1,kind_int_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(ncelcet,1,kind_int_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(nghnode,1,kind_int_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(nghedge,1,kind_int_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(nsponge,1,kind_int_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(sigspg ,1,kind_real_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(nrddsp ,1,kind_int_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(nacous ,1,kind_int_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(nprms  ,1,kind_int_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(ns_mean,1,kind_int_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(ns_prms,1,kind_int_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)

        call MPI_PACK(wgas  ,1,kind_real_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(gamma ,1,kind_real_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(mu0sth,1,kind_real_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(t0sth ,1,kind_real_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(tssth ,1,kind_real_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(prlam ,1,kind_real_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(prtur ,1,kind_real_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(sclam ,1,kind_real_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(sctur ,1,kind_real_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        
        call MPI_PACK(nmoni ,1,kind_int_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(nbmoni,30,kind_int_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(nimoni,30,kind_int_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(njmoni,30,kind_int_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        call MPI_PACK(nkmoni,30,kind_int_mpi,packbuf,packsize,position,MPI_COMM_WORLD,ierr)
    end if

    call MPI_BCAST(position,1,kind_int_mpi,master,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(packbuf,position,MPI_PACKED,master,MPI_COMM_WORLD,ierr)

    if (myid /= master) then
        position = 0

        call MPI_UNPACK(packbuf,packsize,position,casname,len_char_name,kind_char_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,indir  ,len_char_file,kind_char_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,outdir ,len_char_file,kind_char_mpi,MPI_COMM_WORLD,ierr)

        call MPI_UNPACK(packbuf,packsize,position,moo   ,1,kind_real_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,reno  ,1,kind_real_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,alt   ,1,kind_real_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,attack,1,kind_real_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,ayaw  ,1,kind_real_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,rinf  ,1,kind_real_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,vinf  ,1,kind_real_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,pinf  ,1,kind_real_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,tinf  ,1,kind_real_mpi,MPI_COMM_WORLD,ierr)

        call MPI_UNPACK(packbuf,packsize,position,reflen ,1,kind_real_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,reflgrd,1,kind_real_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,freflen,1,kind_real_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,frefsc ,1,kind_real_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,frefxc ,1,kind_real_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,frefyc ,1,kind_real_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,frefzc ,1,kind_real_mpi,MPI_COMM_WORLD,ierr)

        call MPI_UNPACK(packbuf,packsize,position,topsty ,1,kind_int_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,grdsty ,1,kind_int_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,pltsty ,1,kind_int_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,topfile,len_char_file,kind_char_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,grdfile,len_char_file,kind_char_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,resfile,len_char_file,kind_char_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,fcefile,len_char_file,kind_char_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,solfile,len_char_file,kind_char_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,pltfile,len_char_file,kind_char_mpi,MPI_COMM_WORLD,ierr)

        call MPI_UNPACK(packbuf,packsize,position,nincst,1,kind_int_mpi ,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,mincst,1,kind_real_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,alpcst,1,kind_real_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,betcst,1,kind_real_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,pincst,1,kind_real_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,tincst,1,kind_real_mpi,MPI_COMM_WORLD,ierr)

        call MPI_UNPACK(packbuf,packsize,position,nrestrt,1,kind_int_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,nstepmx,1,kind_int_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,nressav,1,kind_int_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,nsolsav,1,kind_int_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,nfcesav,1,kind_int_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,npltsav,1,kind_int_mpi,MPI_COMM_WORLD,ierr)

        call MPI_UNPACK(packbuf,packsize,position,ndim ,1,kind_int_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,nvis ,1,kind_int_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,twall,1,kind_real_mpi,MPI_COMM_WORLD,ierr)

        call MPI_UNPACK(packbuf,packsize,position,nturlhs,1,kind_int_mpi ,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,ntursub,1,kind_int_mpi ,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,ntursch,1,kind_int_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,nturint,1,kind_int_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,nrddst ,1,kind_int_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,ntke2s ,1,kind_int_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,nearsm ,1,kind_int_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,kfstur ,1,kind_real_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,kmaxtur,1,kind_real_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,vfstur ,1,kind_real_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,vmaxtur,1,kind_real_mpi,MPI_COMM_WORLD,ierr)

        call MPI_UNPACK(packbuf,packsize,position,nprec  ,1,kind_int_mpi ,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,nscmp  ,1,kind_int_mpi ,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,nlhs   ,1,kind_int_mpi ,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,nsubmax,1,kind_int_mpi ,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,tolsub ,1,kind_real_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,nunst  ,1,kind_int_mpi ,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,cflst  ,1,kind_real_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,cfled  ,1,kind_real_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,cflfac ,1,kind_real_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,cdtmax ,1,kind_real_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,relaxs ,1,kind_real_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,relaxp ,1,kind_real_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,dtau   ,1,kind_real_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,dtime  ,1,kind_real_mpi,MPI_COMM_WORLD,ierr)

        call MPI_UNPACK(packbuf,packsize,position,nscheme,1,kind_int_mpi ,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,nflux  ,1,kind_int_mpi ,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,nintnon,1,kind_int_mpi ,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,nlimit ,1,kind_int_mpi ,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,enfix  ,1,kind_real_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,ckmuscl,1,kind_real_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,cbmuscl,1,kind_real_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,csrvis ,1,kind_real_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,csafeup,1,kind_real_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,rlimit ,2,kind_real_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,plimit ,2,kind_real_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,tlimit ,2,kind_real_mpi,MPI_COMM_WORLD,ierr)

        call MPI_UNPACK(packbuf,packsize,position,ncutpol,1,kind_int_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,ncelcet,1,kind_int_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,nghnode,1,kind_int_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,nghedge,1,kind_int_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,nsponge,1,kind_int_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,sigspg ,1,kind_real_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,nrddsp ,1,kind_int_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,nacous  ,1,kind_int_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,nprms  ,1,kind_int_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,ns_mean,1,kind_int_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,ns_prms,1,kind_int_mpi,MPI_COMM_WORLD,ierr)

        call MPI_UNPACK(packbuf,packsize,position,wgas  ,1,kind_real_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,gamma ,1,kind_real_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,mu0sth,1,kind_real_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,t0sth ,1,kind_real_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,tssth ,1,kind_real_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,prlam ,1,kind_real_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,prtur ,1,kind_real_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,sclam ,1,kind_real_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,sctur ,1,kind_real_mpi,MPI_COMM_WORLD,ierr)
        
        call MPI_UNPACK(packbuf,packsize,position,nmoni ,1,kind_int_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,nbmoni,30,kind_int_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,nimoni,30,kind_int_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,njmoni,30,kind_int_mpi,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position,nkmoni,30,kind_int_mpi,MPI_COMM_WORLD,ierr)
    end if

#endif
end subroutine broadcast_par


subroutine broadcast_top
#ifdef PARALLEL
    use mod_kndconsts, only : kind_long,kind_int,len_char_name
    use mod_datatypes, only : bc_region_t,top_block_t
    use mod_fieldvars, only : nblocks,nprocs,mb_top
    use mod_parallels
    implicit none
    integer(kind_int),parameter :: packsize=10000
    integer(kind_long)          :: packbuf(packsize)
    integer(kind_int)           :: position
    integer(kind_int)           :: nb,nr,ntms,pid,ierr
    type(top_block_t), pointer  :: top
    type(bc_region_t), pointer  :: reg

    call MPI_BCAST(nblocks,1,kind_int_mpi,master,MPI_COMM_WORLD,ierr)

    if (myid == master) then
        if (nprocs > numprocs) then
            ntms = nprocs/numprocs
            if (mod(nprocs, numprocs) > 0) then
                ntms = ntms + 1
            end if
        else if (nprocs == numprocs) then
            ntms = -1
        else
            ntms = -1
            ierr = 1
            call error_check(ierr, &
                     "nprocs < numprocs in subroutine broadcast_top")
        end if

        nprocs = numprocs

        do nb=1,nblocks
            top => mb_top(nb)

            pid = top%pid
            if (ntms > 0) then
                pid = pid / ntms
            end if

            top%pid = pid
        end do
    else
        allocate(mb_top(nblocks),stat=ierr)
    end if

    call MPI_BCAST(nprocs ,1,kind_int_mpi,master,MPI_COMM_WORLD,ierr)

    ! distribute the info. of blocks<*>
    do nb=1,nblocks
        top => mb_top(nb)
        if (myid == master) then
            position = 0

            call MPI_PACK(top%pid,1,kind_int_mpi, &
                          packbuf,packsize,position,MPI_COMM_WORLD,ierr)
            call MPI_PACK(top%pnb,1,kind_int_mpi, &
                          packbuf,packsize,position,MPI_COMM_WORLD,ierr)
            call MPI_PACK(top%pst,3,kind_int_mpi, &
                          packbuf,packsize,position,MPI_COMM_WORLD,ierr)
            call MPI_PACK(top%ped,3,kind_int_mpi, &
                          packbuf,packsize,position,MPI_COMM_WORLD,ierr)
            call MPI_PACK(top%nupd,1,kind_int_mpi, &
                          packbuf,packsize,position,MPI_COMM_WORLD,ierr)
            call MPI_PACK(top%nijk,3,kind_int_mpi, &
                          packbuf,packsize,position,MPI_COMM_WORLD,ierr)
            call MPI_PACK(top%name,len_char_name,kind_char_mpi, &
                          packbuf,packsize,position,MPI_COMM_WORLD,ierr)
            call MPI_PACK(top%nregions,1,kind_int_mpi, &
                          packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        end if

        ! broadcast the info. of blocks to each processor
        ! from master processor(myid=0)
        call MPI_BCAST(position,1,kind_int_mpi,master,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(packbuf,position,MPI_PACKED,master,MPI_COMM_WORLD,ierr)

        if (myid /= master) then
            position = 0

            call MPI_UNPACK(packbuf,packsize,position, &
                            top%pid,1,kind_int_mpi,MPI_COMM_WORLD,ierr)
            call MPI_UNPACK(packbuf,packsize,position, &
                            top%pnb,1,kind_int_mpi,MPI_COMM_WORLD,ierr)
            call MPI_UNPACK(packbuf,packsize,position, &
                            top%pst,3,kind_int_mpi,MPI_COMM_WORLD,ierr)
            call MPI_UNPACK(packbuf,packsize,position, &
                            top%ped,3,kind_int_mpi,MPI_COMM_WORLD,ierr)
            call MPI_UNPACK(packbuf,packsize,position, &
                            top%nupd,1,kind_int_mpi,MPI_COMM_WORLD,ierr)
            call MPI_UNPACK(packbuf,packsize,position, &
                            top%nijk,3,kind_int_mpi,MPI_COMM_WORLD,ierr)
            call MPI_UNPACK(packbuf,packsize,position, &
                            top%name,len_char_name,kind_char_mpi,MPI_COMM_WORLD,ierr)
            call MPI_UNPACK(packbuf,packsize,position, &
                            top%nregions,1,kind_int_mpi,MPI_COMM_WORLD,ierr)
            allocate(top%bcs(top%nregions),stat=ierr)
        end if

    end do

    do nb=1,nblocks
        top => mb_top(nb)
        do nr=1,top%nregions
            reg => top%bcs(nr)
            if (myid == master) then
                position = 0

                call MPI_PACK(reg%bctype,1,kind_int_mpi, &
                              packbuf,packsize,position,MPI_COMM_WORLD,ierr)
                call MPI_PACK(reg%s_st  ,3,kind_int_mpi, &
                              packbuf,packsize,position,MPI_COMM_WORLD,ierr)
                call MPI_PACK(reg%s_ed  ,3,kind_int_mpi, &
                              packbuf,packsize,position,MPI_COMM_WORLD,ierr)

                call MPI_PACK(reg%nbt   ,1,kind_int_mpi, &
                              packbuf,packsize,position,MPI_COMM_WORLD,ierr)
                call MPI_PACK(reg%t_st  ,3,kind_int_mpi, &
                              packbuf,packsize,position,MPI_COMM_WORLD,ierr)
                call MPI_PACK(reg%t_ed  ,3,kind_int_mpi, &
                              packbuf,packsize,position,MPI_COMM_WORLD,ierr)

                call MPI_PACK(reg%subtype  ,1,kind_int_mpi, &
                              packbuf,packsize,position,MPI_COMM_WORLD,ierr)
                call MPI_PACK(reg%subtype_A,1,kind_int_mpi, &
                              packbuf,packsize,position,MPI_COMM_WORLD,ierr)
            end if

            ! broadcast the info. of BCs to each processor
            ! from master processor(myid=0)
            call MPI_BCAST(position,1,kind_int_mpi,master,MPI_COMM_WORLD,ierr)
            call MPI_BCAST(packbuf,position,MPI_PACKED,master,MPI_COMM_WORLD,ierr)

            if (myid /= master) then
                position = 0

                call MPI_UNPACK(packbuf,packsize,position, &
                                reg%bctype,1,kind_int_mpi,MPI_COMM_WORLD,ierr)
                call MPI_UNPACK(packbuf,packsize,position, &
                                reg%s_st  ,3,kind_int_mpi,MPI_COMM_WORLD,ierr)
                call MPI_UNPACK(packbuf,packsize,position, &
                                reg%s_ed  ,3,kind_int_mpi,MPI_COMM_WORLD,ierr)

                call MPI_UNPACK(packbuf,packsize,position, &
                                reg%nbt   ,1,kind_int_mpi,MPI_COMM_WORLD,ierr)
                call MPI_UNPACK(packbuf,packsize,position, &
                                reg%t_st  ,3,kind_int_mpi,MPI_COMM_WORLD,ierr)
                call MPI_UNPACK(packbuf,packsize,position, &
                                reg%t_ed  ,3,kind_int_mpi,MPI_COMM_WORLD,ierr)

                call MPI_UNPACK(packbuf,packsize,position, &
                                reg%subtype  ,1,kind_int_mpi,MPI_COMM_WORLD,ierr)
                call MPI_UNPACK(packbuf,packsize,position, &
                                reg%subtype_A,1,kind_int_mpi,MPI_COMM_WORLD,ierr)
            end if

        end do
    end do
#endif
end subroutine broadcast_top

subroutine broadcast_top_cc
#ifdef PARALLEL
    use mod_kndconsts, only : kind_long,kind_int,len_char_name
    use mod_datatypes, only : bc_region_t,top_block_t
    use mod_fieldvars, only : nblocks,nprocs,mb_top,mb_topc,mb_topsp
    use mod_parallels
    implicit none
    integer(kind_int),parameter :: packsize=10000
    integer(kind_long)          :: packbuf(packsize)
    integer(kind_int)           :: position
    integer(kind_int)           :: nb,nr,ntms,pid,ierr
    type(top_block_t), pointer  :: top,topc,topsp
    type(bc_region_t), pointer  :: reg,regc,regsp

    call MPI_BCAST(nblocks,1,kind_int_mpi,master,MPI_COMM_WORLD,ierr)

    if (myid == master) then
        if (nprocs > numprocs) then
            ntms = nprocs/numprocs
            if (mod(nprocs, numprocs) > 0) then
                ntms = ntms + 1
            end if
        else if (nprocs == numprocs) then
            ntms = -1
        else
            ntms = -1
            ierr = 1
            call error_check(ierr, &
                     "nprocs < numprocs in subroutine broadcast_top_cc")
        end if

        nprocs = numprocs

        do nb=1,nblocks
            top => mb_top(nb)

            pid = top%pid
            if (ntms > 0) then
                pid = pid / ntms
            end if

            top%pid = pid
        end do
        
        do nb=1,nblocks
            topc => mb_topc(nb)

            pid = topc%pid
            if (ntms > 0) then
                pid = pid / ntms
            end if

            topc%pid = pid
        end do 
        
        do nb=1,nblocks
            topsp => mb_topsp(nb)

            pid = topsp%pid
            if (ntms > 0) then
                pid = pid / ntms
            end if

            topsp%pid = pid
        end do
    else
        allocate(mb_top(nblocks),stat=ierr)
        allocate(mb_topc(nblocks),stat=ierr)
        allocate(mb_topsp(nblocks),stat=ierr)
    end if  

    call MPI_BCAST(nprocs ,1,kind_int_mpi,master,MPI_COMM_WORLD,ierr)

    ! distribute the info. of blocks<*>
    do nb=1,nblocks
        top => mb_top(nb)
        if (myid == master) then
            position = 0

            call MPI_PACK(top%pid,1,kind_int_mpi, &
                          packbuf,packsize,position,MPI_COMM_WORLD,ierr)
            call MPI_PACK(top%pnb,1,kind_int_mpi, &
                          packbuf,packsize,position,MPI_COMM_WORLD,ierr)
            call MPI_PACK(top%pst,3,kind_int_mpi, &
                          packbuf,packsize,position,MPI_COMM_WORLD,ierr)
            call MPI_PACK(top%ped,3,kind_int_mpi, &
                          packbuf,packsize,position,MPI_COMM_WORLD,ierr)
            call MPI_PACK(top%nupd,1,kind_int_mpi, &
                          packbuf,packsize,position,MPI_COMM_WORLD,ierr)
            call MPI_PACK(top%nijk,3,kind_int_mpi, &
                          packbuf,packsize,position,MPI_COMM_WORLD,ierr)
            call MPI_PACK(top%name,len_char_name,kind_char_mpi, &
                          packbuf,packsize,position,MPI_COMM_WORLD,ierr)
            call MPI_PACK(top%nregions,1,kind_int_mpi, &
                          packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        end if

        ! broadcast the info. of blocks to each processor
        ! from master processor(myid=0)
        call MPI_BCAST(position,1,kind_int_mpi,master,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(packbuf,position,MPI_PACKED,master,MPI_COMM_WORLD,ierr)

        if (myid /= master) then
            position = 0

            call MPI_UNPACK(packbuf,packsize,position, &
                            top%pid,1,kind_int_mpi,MPI_COMM_WORLD,ierr)
            call MPI_UNPACK(packbuf,packsize,position, &
                            top%pnb,1,kind_int_mpi,MPI_COMM_WORLD,ierr)
            call MPI_UNPACK(packbuf,packsize,position, &
                            top%pst,3,kind_int_mpi,MPI_COMM_WORLD,ierr)
            call MPI_UNPACK(packbuf,packsize,position, &
                            top%ped,3,kind_int_mpi,MPI_COMM_WORLD,ierr)
            call MPI_UNPACK(packbuf,packsize,position, &
                            top%nupd,1,kind_int_mpi,MPI_COMM_WORLD,ierr)
            call MPI_UNPACK(packbuf,packsize,position, &
                            top%nijk,3,kind_int_mpi,MPI_COMM_WORLD,ierr)
            call MPI_UNPACK(packbuf,packsize,position, &
                            top%name,len_char_name,kind_char_mpi,MPI_COMM_WORLD,ierr)
            call MPI_UNPACK(packbuf,packsize,position, &
                            top%nregions,1,kind_int_mpi,MPI_COMM_WORLD,ierr)
            allocate(top%bcs(top%nregions),stat=ierr)
        end if

    end do

    do nb=1,nblocks
        top => mb_top(nb)
        do nr=1,top%nregions
            reg => top%bcs(nr)
            if (myid == master) then
                position = 0

                call MPI_PACK(reg%bctype,1,kind_int_mpi, &
                              packbuf,packsize,position,MPI_COMM_WORLD,ierr)
                call MPI_PACK(reg%s_st  ,3,kind_int_mpi, &
                              packbuf,packsize,position,MPI_COMM_WORLD,ierr)
                call MPI_PACK(reg%s_ed  ,3,kind_int_mpi, &
                              packbuf,packsize,position,MPI_COMM_WORLD,ierr)

                call MPI_PACK(reg%nbt   ,1,kind_int_mpi, &
                              packbuf,packsize,position,MPI_COMM_WORLD,ierr)
                call MPI_PACK(reg%t_st  ,3,kind_int_mpi, &
                              packbuf,packsize,position,MPI_COMM_WORLD,ierr)
                call MPI_PACK(reg%t_ed  ,3,kind_int_mpi, &
                              packbuf,packsize,position,MPI_COMM_WORLD,ierr)

                call MPI_PACK(reg%subtype  ,1,kind_int_mpi, &
                              packbuf,packsize,position,MPI_COMM_WORLD,ierr)
                call MPI_PACK(reg%subtype_A,1,kind_int_mpi, &
                              packbuf,packsize,position,MPI_COMM_WORLD,ierr)
            end if

            ! broadcast the info. of BCs to each processor
            ! from master processor(myid=0)
            call MPI_BCAST(position,1,kind_int_mpi,master,MPI_COMM_WORLD,ierr)
            call MPI_BCAST(packbuf,position,MPI_PACKED,master,MPI_COMM_WORLD,ierr)

            if (myid /= master) then
                position = 0

                call MPI_UNPACK(packbuf,packsize,position, &
                                reg%bctype,1,kind_int_mpi,MPI_COMM_WORLD,ierr)
                call MPI_UNPACK(packbuf,packsize,position, &
                                reg%s_st  ,3,kind_int_mpi,MPI_COMM_WORLD,ierr)
                call MPI_UNPACK(packbuf,packsize,position, &
                                reg%s_ed  ,3,kind_int_mpi,MPI_COMM_WORLD,ierr)

                call MPI_UNPACK(packbuf,packsize,position, &
                                reg%nbt   ,1,kind_int_mpi,MPI_COMM_WORLD,ierr)
                call MPI_UNPACK(packbuf,packsize,position, &
                                reg%t_st  ,3,kind_int_mpi,MPI_COMM_WORLD,ierr)
                call MPI_UNPACK(packbuf,packsize,position, &
                                reg%t_ed  ,3,kind_int_mpi,MPI_COMM_WORLD,ierr)

                call MPI_UNPACK(packbuf,packsize,position, &
                                reg%subtype  ,1,kind_int_mpi,MPI_COMM_WORLD,ierr)
                call MPI_UNPACK(packbuf,packsize,position, &
                                reg%subtype_A,1,kind_int_mpi,MPI_COMM_WORLD,ierr)
            end if

        end do
    end do
    
    do nb=1,nblocks
        topc => mb_topc(nb)
        if (myid == master) then
            position = 0

            call MPI_PACK(topc%pid,1,kind_int_mpi, &
                          packbuf,packsize,position,MPI_COMM_WORLD,ierr)
            call MPI_PACK(topc%pnb,1,kind_int_mpi, &
                          packbuf,packsize,position,MPI_COMM_WORLD,ierr)
            call MPI_PACK(topc%pst,3,kind_int_mpi, &
                          packbuf,packsize,position,MPI_COMM_WORLD,ierr)
            call MPI_PACK(topc%ped,3,kind_int_mpi, &
                          packbuf,packsize,position,MPI_COMM_WORLD,ierr)
            call MPI_PACK(topc%nupd,1,kind_int_mpi, &
                          packbuf,packsize,position,MPI_COMM_WORLD,ierr)
            call MPI_PACK(topc%nijk,3,kind_int_mpi, &
                          packbuf,packsize,position,MPI_COMM_WORLD,ierr)
            call MPI_PACK(topc%name,len_char_name,kind_char_mpi, &
                          packbuf,packsize,position,MPI_COMM_WORLD,ierr)
            call MPI_PACK(topc%nregions,1,kind_int_mpi, &
                          packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        end if

        ! broadcast the info. of blocks to each processor
        ! from master processor(myid=0)
        call MPI_BCAST(position,1,kind_int_mpi,master,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(packbuf,position,MPI_PACKED,master,MPI_COMM_WORLD,ierr)

        if (myid /= master) then
            position = 0

            call MPI_UNPACK(packbuf,packsize,position, &
                            topc%pid,1,kind_int_mpi,MPI_COMM_WORLD,ierr)
            call MPI_UNPACK(packbuf,packsize,position, &
                            topc%pnb,1,kind_int_mpi,MPI_COMM_WORLD,ierr)
            call MPI_UNPACK(packbuf,packsize,position, &
                            topc%pst,3,kind_int_mpi,MPI_COMM_WORLD,ierr)
            call MPI_UNPACK(packbuf,packsize,position, &
                            topc%ped,3,kind_int_mpi,MPI_COMM_WORLD,ierr)
            call MPI_UNPACK(packbuf,packsize,position, &
                            topc%nupd,1,kind_int_mpi,MPI_COMM_WORLD,ierr)
            call MPI_UNPACK(packbuf,packsize,position, &
                            topc%nijk,3,kind_int_mpi,MPI_COMM_WORLD,ierr)
            call MPI_UNPACK(packbuf,packsize,position, &
                            topc%name,len_char_name,kind_char_mpi,MPI_COMM_WORLD,ierr)
            call MPI_UNPACK(packbuf,packsize,position, &
                            topc%nregions,1,kind_int_mpi,MPI_COMM_WORLD,ierr)
            allocate(topc%bcs(topc%nregions),stat=ierr)
        end if

    end do

    do nb=1,nblocks
        topc => mb_topc(nb)
        do nr=1,topc%nregions
            regc => topc%bcs(nr)
            if (myid == master) then
                position = 0

                call MPI_PACK(regc%bctype,1,kind_int_mpi, &
                              packbuf,packsize,position,MPI_COMM_WORLD,ierr)
                call MPI_PACK(regc%s_st  ,3,kind_int_mpi, &
                              packbuf,packsize,position,MPI_COMM_WORLD,ierr)
                call MPI_PACK(regc%s_ed  ,3,kind_int_mpi, &
                              packbuf,packsize,position,MPI_COMM_WORLD,ierr)

                call MPI_PACK(regc%nbt   ,1,kind_int_mpi, &
                              packbuf,packsize,position,MPI_COMM_WORLD,ierr)
                call MPI_PACK(regc%t_st  ,3,kind_int_mpi, &
                              packbuf,packsize,position,MPI_COMM_WORLD,ierr)
                call MPI_PACK(regc%t_ed  ,3,kind_int_mpi, &
                              packbuf,packsize,position,MPI_COMM_WORLD,ierr)

                call MPI_PACK(regc%subtype  ,1,kind_int_mpi, &
                              packbuf,packsize,position,MPI_COMM_WORLD,ierr)
                call MPI_PACK(regc%subtype_A,1,kind_int_mpi, &
                              packbuf,packsize,position,MPI_COMM_WORLD,ierr)
            end if

            ! broadcast the info. of BCs to each processor
            ! from master processor(myid=0)
            call MPI_BCAST(position,1,kind_int_mpi,master,MPI_COMM_WORLD,ierr)
            call MPI_BCAST(packbuf,position,MPI_PACKED,master,MPI_COMM_WORLD,ierr)

            if (myid /= master) then
                position = 0

                call MPI_UNPACK(packbuf,packsize,position, &
                                regc%bctype,1,kind_int_mpi,MPI_COMM_WORLD,ierr)
                call MPI_UNPACK(packbuf,packsize,position, &
                                regc%s_st  ,3,kind_int_mpi,MPI_COMM_WORLD,ierr)
                call MPI_UNPACK(packbuf,packsize,position, &
                                regc%s_ed  ,3,kind_int_mpi,MPI_COMM_WORLD,ierr)

                call MPI_UNPACK(packbuf,packsize,position, &
                                regc%nbt   ,1,kind_int_mpi,MPI_COMM_WORLD,ierr)
                call MPI_UNPACK(packbuf,packsize,position, &
                                regc%t_st  ,3,kind_int_mpi,MPI_COMM_WORLD,ierr)
                call MPI_UNPACK(packbuf,packsize,position, &
                                regc%t_ed  ,3,kind_int_mpi,MPI_COMM_WORLD,ierr)

                call MPI_UNPACK(packbuf,packsize,position, &
                                regc%subtype  ,1,kind_int_mpi,MPI_COMM_WORLD,ierr)
                call MPI_UNPACK(packbuf,packsize,position, &
                                regc%subtype_A,1,kind_int_mpi,MPI_COMM_WORLD,ierr)
            end if

        end do
    end do
    
    do nb=1,nblocks
        topsp => mb_topsp(nb)
        if (myid == master) then
            position = 0

            call MPI_PACK(topsp%pid,1,kind_int_mpi, &
                          packbuf,packsize,position,MPI_COMM_WORLD,ierr)
            call MPI_PACK(topsp%pnb,1,kind_int_mpi, &
                          packbuf,packsize,position,MPI_COMM_WORLD,ierr)
            call MPI_PACK(topsp%pst,3,kind_int_mpi, &
                          packbuf,packsize,position,MPI_COMM_WORLD,ierr)
            call MPI_PACK(topsp%ped,3,kind_int_mpi, &
                          packbuf,packsize,position,MPI_COMM_WORLD,ierr)
            call MPI_PACK(topsp%nupd,1,kind_int_mpi, &
                          packbuf,packsize,position,MPI_COMM_WORLD,ierr)
            call MPI_PACK(topsp%nijk,3,kind_int_mpi, &
                          packbuf,packsize,position,MPI_COMM_WORLD,ierr)
            call MPI_PACK(topsp%name,len_char_name,kind_char_mpi, &
                          packbuf,packsize,position,MPI_COMM_WORLD,ierr)
            call MPI_PACK(topsp%nregions,1,kind_int_mpi, &
                          packbuf,packsize,position,MPI_COMM_WORLD,ierr)
        end if

        ! broadcast the info. of blocks to each processor
        ! from master processor(myid=0)
        call MPI_BCAST(position,1,kind_int_mpi,master,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(packbuf,position,MPI_PACKED,master,MPI_COMM_WORLD,ierr)

        if (myid /= master) then
            position = 0

            call MPI_UNPACK(packbuf,packsize,position, &
                            topsp%pid,1,kind_int_mpi,MPI_COMM_WORLD,ierr)
            call MPI_UNPACK(packbuf,packsize,position, &
                            topsp%pnb,1,kind_int_mpi,MPI_COMM_WORLD,ierr)
            call MPI_UNPACK(packbuf,packsize,position, &
                            topsp%pst,3,kind_int_mpi,MPI_COMM_WORLD,ierr)
            call MPI_UNPACK(packbuf,packsize,position, &
                            topsp%ped,3,kind_int_mpi,MPI_COMM_WORLD,ierr)
            call MPI_UNPACK(packbuf,packsize,position, &
                            topsp%nupd,1,kind_int_mpi,MPI_COMM_WORLD,ierr)
            call MPI_UNPACK(packbuf,packsize,position, &
                            topsp%nijk,3,kind_int_mpi,MPI_COMM_WORLD,ierr)
            call MPI_UNPACK(packbuf,packsize,position, &
                            topsp%name,len_char_name,kind_char_mpi,MPI_COMM_WORLD,ierr)
            call MPI_UNPACK(packbuf,packsize,position, &
                            topsp%nregions,1,kind_int_mpi,MPI_COMM_WORLD,ierr)
            allocate(topsp%bcs(topsp%nregions),stat=ierr)
        end if

    end do

    do nb=1,nblocks
        topsp => mb_topsp(nb)
        do nr=1,topsp%nregions
            regsp => topsp%bcs(nr)
            if (myid == master) then
                position = 0

                call MPI_PACK(regsp%bctype,1,kind_int_mpi, &
                              packbuf,packsize,position,MPI_COMM_WORLD,ierr)
                call MPI_PACK(regsp%s_st  ,3,kind_int_mpi, &
                              packbuf,packsize,position,MPI_COMM_WORLD,ierr)
                call MPI_PACK(regsp%s_ed  ,3,kind_int_mpi, &
                              packbuf,packsize,position,MPI_COMM_WORLD,ierr)

                call MPI_PACK(regsp%nbt   ,1,kind_int_mpi, &
                              packbuf,packsize,position,MPI_COMM_WORLD,ierr)
                call MPI_PACK(regsp%t_st  ,3,kind_int_mpi, &
                              packbuf,packsize,position,MPI_COMM_WORLD,ierr)
                call MPI_PACK(regsp%t_ed  ,3,kind_int_mpi, &
                              packbuf,packsize,position,MPI_COMM_WORLD,ierr)

                call MPI_PACK(regsp%subtype  ,1,kind_int_mpi, &
                              packbuf,packsize,position,MPI_COMM_WORLD,ierr)
                call MPI_PACK(regsp%subtype_A,1,kind_int_mpi, &
                              packbuf,packsize,position,MPI_COMM_WORLD,ierr)
            end if

            ! broadcast the info. of BCs to each processor
            ! from master processor(myid=0)
            call MPI_BCAST(position,1,kind_int_mpi,master,MPI_COMM_WORLD,ierr)
            call MPI_BCAST(packbuf,position,MPI_PACKED,master,MPI_COMM_WORLD,ierr)

            if (myid /= master) then
                position = 0

                call MPI_UNPACK(packbuf,packsize,position, &
                                regsp%bctype,1,kind_int_mpi,MPI_COMM_WORLD,ierr)
                call MPI_UNPACK(packbuf,packsize,position, &
                                regsp%s_st  ,3,kind_int_mpi,MPI_COMM_WORLD,ierr)
                call MPI_UNPACK(packbuf,packsize,position, &
                                regsp%s_ed  ,3,kind_int_mpi,MPI_COMM_WORLD,ierr)

                call MPI_UNPACK(packbuf,packsize,position, &
                                regsp%nbt   ,1,kind_int_mpi,MPI_COMM_WORLD,ierr)
                call MPI_UNPACK(packbuf,packsize,position, &
                                regsp%t_st  ,3,kind_int_mpi,MPI_COMM_WORLD,ierr)
                call MPI_UNPACK(packbuf,packsize,position, &
                                regsp%t_ed  ,3,kind_int_mpi,MPI_COMM_WORLD,ierr)

                call MPI_UNPACK(packbuf,packsize,position, &
                                regsp%subtype  ,1,kind_int_mpi,MPI_COMM_WORLD,ierr)
                call MPI_UNPACK(packbuf,packsize,position, &
                                regsp%subtype_A,1,kind_int_mpi,MPI_COMM_WORLD,ierr)
            end if

        end do
    end do    
#endif
end subroutine broadcast_top_cc

subroutine scatter_input_mb_var(io_unit,mb_var,nst,ned,nfun)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : fldtype_r3d,nplot3d_std
    use mod_datatypes, only : var_block_t,fld_array_t,top_block_t
    use mod_fieldvars, only : nblocks,mb_top
    use mod_parallels
    implicit none
    integer(kind_int),          intent(in) :: io_unit
    type(var_block_t), pointer, intent(in) :: mb_var(:)
    integer(kind_int),          intent(in) :: nst,ned,nfun

    ! local storage
    integer(kind_int)          :: nb,m,n,nvar,nv,ierr
    integer(kind_int)          :: nblkvars,ni,nj,nk
    integer(kind_int)          :: i,j,k,st(3),ed(3)
    type(top_block_t), pointer :: top
    type(fld_array_t), pointer :: var(:)
    real(kind_real),   pointer :: packbuf(:,:,:,:)

#ifdef PARALLEL
    integer(kind_int) :: packsize,pid
    integer(kind_int) :: status(MPI_STATUS_SIZE)
#endif

    nvar = ned - nst + 1

#ifdef PARALLEL
    if (myid == master) then
#endif
        read(io_unit, iostat=ierr) nblkvars
        call error_check(nblkvars-nblocks, &
                         "The number of blocks is incorrect in subroutine scatter_input_mb_var")

        do nb=1,nblocks
            top => mb_top(nb)

            if (nfun == nplot3d_std) then
                read(io_unit, iostat=ierr) ni, nj, nk
                n = 0
            else
                read(io_unit, iostat=ierr) ni, nj, nk, nv
                n = nvar - nv
            end if

            m = abs(top%nijk(1) - ni) + &
                abs(top%nijk(2) - nj) + &
                abs(top%nijk(3) - nk) + abs(n)
            call error_check(m,"The dimension of blocks is incorrect in subroutine scatter_input_mb_var")
        end do
#ifdef PARALLEL
    end if
#endif

    do nb=1,nblocks
        top => mb_top(nb)
        var => mb_var(nb)%fld

        ni = top%nijk(1)
        nj = top%nijk(2)
        nk = top%nijk(3)

#ifdef PARALLEL
        packsize = product(top%nijk(:))*nvar

        pid = top%pid
        if (myid == master) then
#endif
            allocate(packbuf(1:ni, &
                             1:nj, &
                             1:nk, &
                             nst:ned), stat=ierr)
            read(io_unit, iostat=ierr) ((((packbuf(i,j,k,m), i=1,ni), &
                                                             j=1,nj), &
                                                             k=1,nk), &
                                                             m=nst,ned)
#ifdef PARALLEL
        end if

        if (pid /= master) then
            if (myid == master) then
                call MPI_SEND(packbuf,packsize,kind_real_mpi, &
                              pid,nb,MPI_COMM_WORLD,ierr)
            end if
            if (myid == pid) then
                allocate(packbuf(1:ni, &
                                 1:nj, &
                                 1:nk, &
                                 nst:ned), stat=ierr)
                call MPI_RECV(packbuf,packsize,kind_real_mpi, &
                              master,nb,MPI_COMM_WORLD,status,ierr)
            end if

        end if

        if (myid == pid) then
#endif
            do m=nst,ned
                do k=1,nk
                do j=1,nj
                do i=1,ni
                    var(m)%r3d(i,j,k) = packbuf(i,j,k,m)
                end do
                end do
                end do
            end do

            deallocate(packbuf, stat=ierr)
#ifdef PARALLEL
        end if

        if (pid /= master .and. myid == master) then
            deallocate(packbuf, stat=ierr)
        end if
#endif
    end do

end subroutine scatter_input_mb_var

subroutine scatter_input_mb_var_cc(io_unit,mb_var,nst,ned,nfun)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : fldtype_r3d,nplot3d_std
    use mod_datatypes, only : var_block_t,fld_array_t,top_block_t
    use mod_fieldvars, only : nblocks,mb_topc
    use mod_parallels
    implicit none
    integer(kind_int),          intent(in) :: io_unit
    type(var_block_t), pointer, intent(in) :: mb_var(:)
    integer(kind_int),          intent(in) :: nst,ned,nfun

    ! local storage
    integer(kind_int)          :: nb,m,n,nvar,nv,ierr
    integer(kind_int)          :: nblkvars,ni,nj,nk
    integer(kind_int)          :: i,j,k,st(3),ed(3)
    type(top_block_t), pointer :: top
    type(fld_array_t), pointer :: var(:)
    real(kind_real),   pointer :: packbuf(:,:,:,:)

#ifdef PARALLEL
    integer(kind_int) :: packsize,pid
    integer(kind_int) :: status(MPI_STATUS_SIZE)
#endif

    nvar = ned - nst + 1

#ifdef PARALLEL
    if (myid == master) then
#endif
        read(io_unit, iostat=ierr) nblkvars
        call error_check(nblkvars-nblocks, &
                         "The number of blocks is incorrect in subroutine scatter_input_mb_var")

        do nb=1,nblocks
            top => mb_topc(nb)

            if (nfun == nplot3d_std) then
                read(io_unit, iostat=ierr) ni, nj, nk
                n = 0
            else
                read(io_unit, iostat=ierr) ni, nj, nk, nv
                n = nvar - nv
            end if

            m = abs(top%nijk(1) - ni) + &
                abs(top%nijk(2) - nj) + &
                abs(top%nijk(3) - nk) + abs(n)
            call error_check(m,"The dimension of blocks is incorrect in subroutine scatter_input_mb_var")
        end do
#ifdef PARALLEL
    end if
#endif

    do nb=1,nblocks
        top => mb_topc(nb)
        var => mb_var(nb)%fld

        ni = top%nijk(1)
        nj = top%nijk(2)
        nk = top%nijk(3)

#ifdef PARALLEL
        packsize = product(top%nijk(:))*nvar

        pid = top%pid
        if (myid == master) then
#endif
            allocate(packbuf(1:ni, &
                             1:nj, &
                             1:nk, &
                             nst:ned), stat=ierr)
            read(io_unit, iostat=ierr) ((((packbuf(i,j,k,m), i=1,ni), &
                                                             j=1,nj), &
                                                             k=1,nk), &
                                                             m=nst,ned)
#ifdef PARALLEL
        end if

        if (pid /= master) then
            if (myid == master) then
                call MPI_SEND(packbuf,packsize,kind_real_mpi, &
                              pid,nb,MPI_COMM_WORLD,ierr)
            end if
            if (myid == pid) then
                allocate(packbuf(1:ni, &
                                 1:nj, &
                                 1:nk, &
                                 nst:ned), stat=ierr)
                call MPI_RECV(packbuf,packsize,kind_real_mpi, &
                              master,nb,MPI_COMM_WORLD,status,ierr)
            end if

        end if

        if (myid == pid) then
#endif
            do m=nst,ned
                do k=1,nk
                do j=1,nj
                do i=1,ni
                    var(m)%r3d(i,j,k) = packbuf(i,j,k,m)
                end do
                end do
                end do
            end do

            deallocate(packbuf, stat=ierr)
#ifdef PARALLEL
        end if

        if (pid /= master .and. myid == master) then
            deallocate(packbuf, stat=ierr)
        end if
#endif
    end do

end subroutine scatter_input_mb_var_cc

subroutine scatter_input_mb_var_sp(io_unit,mb_var,nst,ned,nfun)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : fldtype_r3d,nplot3d_std
    use mod_datatypes, only : var_block_t,fld_array_t,top_block_t
    use mod_fieldvars, only : nblocks,mb_topsp
    use mod_parallels
    implicit none
    integer(kind_int),          intent(in) :: io_unit
    type(var_block_t), pointer, intent(in) :: mb_var(:)
    integer(kind_int),          intent(in) :: nst,ned,nfun

    ! local storage
    integer(kind_int)          :: nb,m,n,nvar,nv,ierr
    integer(kind_int)          :: nblkvars,ni,nj,nk
    integer(kind_int)          :: i,j,k,st(3),ed(3)
    type(top_block_t), pointer :: top
    type(fld_array_t), pointer :: var(:)
    real(kind_real),   pointer :: packbuf(:,:,:,:)

#ifdef PARALLEL
    integer(kind_int) :: packsize,pid
    integer(kind_int) :: status(MPI_STATUS_SIZE)
#endif

    nvar = ned - nst + 1

#ifdef PARALLEL
    if (myid == master) then
#endif
        read(io_unit, iostat=ierr) nblkvars
        call error_check(nblkvars-nblocks, &
                         "The number of blocks is incorrect in subroutine scatter_input_mb_var")

        do nb=1,nblocks
            top => mb_topsp(nb)

            if (nfun == nplot3d_std) then
                read(io_unit, iostat=ierr) ni, nj, nk
                n = 0
            else
                read(io_unit, iostat=ierr) ni, nj, nk, nv
                n = nvar - nv
            end if

            m = abs(top%nijk(1) - ni) + &
                abs(top%nijk(2) - nj) + &
                abs(top%nijk(3) - nk) + abs(n)
            call error_check(m,"The dimension of blocks is incorrect in subroutine scatter_input_mb_var")
        end do
#ifdef PARALLEL
    end if
#endif

    do nb=1,nblocks
        top => mb_topsp(nb)
        var => mb_var(nb)%fld

        ni = top%nijk(1)
        nj = top%nijk(2)
        nk = top%nijk(3)

#ifdef PARALLEL
        packsize = product(top%nijk(:))*nvar

        pid = top%pid
        if (myid == master) then
#endif
            allocate(packbuf(1:ni, &
                             1:nj, &
                             1:nk, &
                             nst:ned), stat=ierr)
            read(io_unit, iostat=ierr) ((((packbuf(i,j,k,m), i=1,ni), &
                                                             j=1,nj), &
                                                             k=1,nk), &
                                                             m=nst,ned)
#ifdef PARALLEL
        end if

        if (pid /= master) then
            if (myid == master) then
                call MPI_SEND(packbuf,packsize,kind_real_mpi, &
                              pid,nb,MPI_COMM_WORLD,ierr)
            end if
            if (myid == pid) then
                allocate(packbuf(1:ni, &
                                 1:nj, &
                                 1:nk, &
                                 nst:ned), stat=ierr)
                call MPI_RECV(packbuf,packsize,kind_real_mpi, &
                              master,nb,MPI_COMM_WORLD,status,ierr)
            end if

        end if

        if (myid == pid) then
#endif
            do m=nst,ned
                do k=1,nk
                do j=1,nj
                do i=1,ni
                    var(m)%r3d(i,j,k) = packbuf(i,j,k,m)
                end do
                end do
                end do
            end do

            deallocate(packbuf, stat=ierr)
#ifdef PARALLEL
        end if

        if (pid /= master .and. myid == master) then
            deallocate(packbuf, stat=ierr)
        end if
#endif
    end do

end subroutine scatter_input_mb_var_sp

subroutine broadcast_input_aux(io_unit)
    use mod_kndconsts, only : kind_int
    use mod_variables, only : nstepsav,cfl
    use mod_variables, only : nsub,nsubp,nsubp1
    use mod_fieldvars, only : nstepmean,nsteprms,ntime
    use mod_runtimers, only : time_cpu_sav,time_wall_sav
    use mod_parallels
    implicit none
    integer(kind_int), intent(in) :: io_unit
    integer(kind_int) :: ierr

#ifdef PARALLEL
    if (myid == master) then
#endif
        read(io_unit,*) nstepsav,nstepmean,nsteprms,ntime,cfl
        read(io_unit,*) nsub,nsubp,nsubp1,time_cpu_sav,time_wall_sav
#ifdef PARALLEL
    end if

    call MPI_BCAST(nstepsav     ,1,kind_int_mpi,master,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(nstepmean    ,1,kind_int_mpi,master,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(nsteprms     ,1,kind_int_mpi,master,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(ntime        ,1,kind_int_mpi,master,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(cfl          ,1,kind_real_mpi,master,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(time_cpu_sav ,1,kind_real_mpi,master,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(time_wall_sav,1,kind_real_mpi,master,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(nsub         ,1,kind_int_mpi,master,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(nsubp        ,1,kind_int_mpi,master,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(nsubp1       ,1,kind_int_mpi,master,MPI_COMM_WORLD,ierr)    
#endif

end subroutine broadcast_input_aux

subroutine broadcast_timestep
    use mod_kndconsts, only : kind_int,kind_real
    use mod_variables, only : dtaumin,dtaumax
    use mod_parallels
    implicit none
    integer(kind_int) :: ierr
    real(kind_real)   :: dtaumingb,dtaumaxgb

#ifdef PARALLEL
    call MPI_REDUCE(dtaumin,dtaumingb,1,kind_real_mpi,MPI_MIN,master,MPI_COMM_WORLD,ierr)
    call MPI_REDUCE(dtaumax,dtaumaxgb,1,kind_real_mpi,MPI_MAX,master,MPI_COMM_WORLD,ierr)

    call MPI_BCAST(dtaumingb,1,kind_real_mpi,master,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(dtaumaxgb,1,kind_real_mpi,master,MPI_COMM_WORLD,ierr)

    dtaumin = dtaumingb
    dtaumax = dtaumaxgb
#endif

end subroutine broadcast_timestep

subroutine broadcast_sponge
    use mod_kndconsts, only : kind_int,kind_real
    use mod_variables, only : dspmin,dspmax
    use mod_parallels
    implicit none
    integer(kind_int) :: ierr
    real(kind_real)   :: dspmingb,dspmaxgb

#ifdef PARALLEL
    call MPI_REDUCE(dspmin,dspmingb,1,kind_real_mpi,MPI_MIN,master,MPI_COMM_WORLD,ierr)
    call MPI_REDUCE(dspmax,dspmaxgb,1,kind_real_mpi,MPI_MAX,master,MPI_COMM_WORLD,ierr)

    call MPI_BCAST(dspmingb,1,kind_real_mpi,master,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(dspmaxgb,1,kind_real_mpi,master,MPI_COMM_WORLD,ierr)

    dspmin = dspmingb
    dspmax = dspmaxgb
#endif

end subroutine broadcast_sponge

subroutine gather_output_mb_var(io_unit,mb_var,nst,ned,nfun)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : fldtype_r3d,nplot3d_std
    use mod_datatypes, only : var_block_t,fld_array_t,top_block_t
    use mod_fieldvars, only : nblocks,mb_top
    use mod_parallels
    implicit none
    integer(kind_int),          intent(in) :: io_unit
    type(var_block_t), pointer, intent(in) :: mb_var(:)
    integer(kind_int),          intent(in) :: nst,ned,nfun

    ! local storage
    integer(kind_int)          :: nb,m,nvar,ierr
    integer(kind_int)          :: nblkvars,ni,nj,nk
    integer(kind_int)          :: i,j,k,st(3),ed(3)
    type(top_block_t), pointer :: top
    type(fld_array_t), pointer :: var(:)
    real(kind_real),   pointer :: packbuf(:,:,:,:)

#ifdef PARALLEL
    integer(kind_int) :: packsize,pid
    integer(kind_int) :: status(MPI_STATUS_SIZE)
#endif

    nvar = ned - nst + 1

#ifdef PARALLEL
    if (myid == master) then
#endif
        write(io_unit, iostat=ierr) nblocks

        do nb=1,nblocks
            top => mb_top(nb)

            if (nfun == nplot3d_std) then
                write(io_unit, iostat=ierr) top%nijk(:)
            else
                write(io_unit, iostat=ierr) top%nijk(:),nvar
            end if
        end do
#ifdef PARALLEL
    end if
#endif

    do nb=1,nblocks
        top => mb_top(nb)
        var => mb_var(nb)%fld

        ni = top%nijk(1)
        nj = top%nijk(2)
        nk = top%nijk(3)

#ifdef PARALLEL
        packsize = product(top%nijk(:))*nvar

        pid = top%pid
        if (myid == pid) then
#endif
            allocate(packbuf(1:ni, &
                             1:nj, &
                             1:nk, &
                             nst:ned), stat=ierr)
            do m=nst,ned
                do k=1,nk
                do j=1,nj
                do i=1,ni
                    packbuf(i,j,k,m) = var(m)%r3d(i,j,k)
                end do
                end do
                end do
            end do
#ifdef PARALLEL
        end if

        if (pid /= master) then
            if (myid == master) then
                allocate(packbuf(1:ni, &
                                 1:nj, &
                                 1:nk, &
                                 nst:ned), stat=ierr)
                call MPI_RECV(packbuf,packsize,kind_real_mpi, &
                              pid,nb,MPI_COMM_WORLD,status,ierr)
            end if
            if (myid == pid) then
                call MPI_SEND(packbuf,packsize,kind_real_mpi, &
                              master,nb,MPI_COMM_WORLD,ierr)
            end if

        end if


        if (myid == master) then
#endif
            write(io_unit, iostat=ierr) ((((packbuf(i,j,k,m), i=1,ni), &
                                                              j=1,nj), &
                                                              k=1,nk), &
                                                              m=nst,ned)
            deallocate(packbuf, stat=ierr)
#ifdef PARALLEL
        end if

        if (pid /= master .and. myid == pid) then
            deallocate(packbuf, stat=ierr)
        end if
#endif
    end do

end subroutine gather_output_mb_var

subroutine gather_output_mb_var_sp(io_unit,mb_var,nst,ned,nfun)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : fldtype_r3d,nplot3d_std
    use mod_datatypes, only : var_block_t,fld_array_t,top_block_t
    use mod_fieldvars, only : nblocks,mb_topsp
    use mod_parallels
    implicit none
    integer(kind_int),          intent(in) :: io_unit
    type(var_block_t), pointer, intent(in) :: mb_var(:)
    integer(kind_int),          intent(in) :: nst,ned,nfun

    ! local storage
    integer(kind_int)          :: nb,m,nvar,ierr
    integer(kind_int)          :: nblkvars,ni,nj,nk
    integer(kind_int)          :: i,j,k,st(3),ed(3)
    type(top_block_t), pointer :: top
    type(fld_array_t), pointer :: var(:)
    real(kind_real),   pointer :: packbuf(:,:,:,:)

#ifdef PARALLEL
    integer(kind_int) :: packsize,pid
    integer(kind_int) :: status(MPI_STATUS_SIZE)
#endif

    nvar = ned - nst + 1

#ifdef PARALLEL
    if (myid == master) then
#endif
        write(io_unit, iostat=ierr) nblocks

        do nb=1,nblocks
            top => mb_topsp(nb)

            if (nfun == nplot3d_std) then
                write(io_unit, iostat=ierr) top%nijk(:)
            else
                write(io_unit, iostat=ierr) top%nijk(:),nvar
            end if
        end do
#ifdef PARALLEL
    end if
#endif

    do nb=1,nblocks
        top => mb_topsp(nb)
        var => mb_var(nb)%fld

        ni = top%nijk(1)
        nj = top%nijk(2)
        nk = top%nijk(3)

#ifdef PARALLEL
        packsize = product(top%nijk(:))*nvar

        pid = top%pid
        if (myid == pid) then
#endif
            allocate(packbuf(1:ni, &
                             1:nj, &
                             1:nk, &
                             nst:ned), stat=ierr)
            do m=nst,ned
                do k=1,nk
                do j=1,nj
                do i=1,ni
                    packbuf(i,j,k,m) = var(m)%r3d(i,j,k)
                end do
                end do
                end do
            end do
#ifdef PARALLEL
        end if

        if (pid /= master) then
            if (myid == master) then
                allocate(packbuf(1:ni, &
                                 1:nj, &
                                 1:nk, &
                                 nst:ned), stat=ierr)
                call MPI_RECV(packbuf,packsize,kind_real_mpi, &
                              pid,nb,MPI_COMM_WORLD,status,ierr)
            end if
            if (myid == pid) then
                call MPI_SEND(packbuf,packsize,kind_real_mpi, &
                              master,nb,MPI_COMM_WORLD,ierr)
            end if

        end if


        if (myid == master) then
#endif
            write(io_unit, iostat=ierr) ((((packbuf(i,j,k,m), i=1,ni), &
                                                              j=1,nj), &
                                                              k=1,nk), &
                                                              m=nst,ned)
            deallocate(packbuf, stat=ierr)
#ifdef PARALLEL
        end if

        if (pid /= master .and. myid == pid) then
            deallocate(packbuf, stat=ierr)
        end if
#endif
    end do

end subroutine gather_output_mb_var_sp

!>
!! @brief 交换边界流场信息函数
!! @details 以对接边界为基本单位，逐一交换流场信息
!! @param[in] mb_var  场变量指针
!! @param[in] nst     场变量起始索引
!! @param[in] ned     场变量终止索引
!! @param[in] ngh     边界交换信息层数
!! @param[in] nswmem  通讯内存分配方式
!! @param[in] nswave  边界信息平均方式
!! @author 刘化勇
!! @date 2011年05月31日
!! @author 王勇献
!! @remark 2012年02月14日 添加非阻塞式通信
!! @remark 2012年04月05日 TODO: 待删除
!>
subroutine exchange_bc_var(mb_var,nst,ned,ngh,nswmem,nswave)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_inter_buf_dyn
    use mod_constants, only : one,nsgl_aver_art
    use mod_datatypes, only : var_block_t
    use mod_fieldvars, only : mb_top,mb_vol
    use mod_fieldvars, only : ninters,inters
    use mod_parallels
    implicit none
    type(var_block_t), pointer, intent(in) :: mb_var(:)
    integer(kind_int),          intent(in) :: nst,ned,ngh
    integer(kind_int),          intent(in) :: nswmem
    integer(kind_int),          intent(in) :: nswave
    integer(kind_int) :: nint,ierr
    integer(kind_int) :: nbs,nrs,s_st(3),s_ed(3),s_nd,s_lr
    integer(kind_int) :: nbt,nrt,t_st(3),t_ed(3),t_nd,t_lr
    integer(kind_int) :: st(3),ed(3)
    integer(kind_int) :: ijks(3),ijkt(3),ijkgs(3),ijkgt(3)
    integer(kind_int) :: i,j,k,m,n,is,js,ks,it,jt,kt
    real(kind_real)   :: vol
#ifdef PARALLEL
    integer(kind_int) :: id_src,id_des,tag,packsize
    integer(kind_int) :: index_msg
#endif

#ifdef PARALLEL
    index_msg = 0
#endif

    do nint=1,ninters
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

        call bc_extend_inward(s_st,s_ed,s_nd,s_lr,ngh,st,ed)

#ifdef PARALLEL
        packsize = product(ed(:)-st(:)+1)*(ned-nst+1)

        id_src = mb_top(nbs)%pid
        id_des = mb_top(nbt)%pid
        tag    = nint + ninters*3


        if (id_src /= id_des) then
            if (myid == id_src) then
                if (nswmem == nbc_inter_buf_dyn) then
                    allocate(inters(nint)%dat(st(1):ed(1), &
                                              st(2):ed(2), &
                                              st(3):ed(3), &
                                              nst:ned),stat=ierr)
                else
                    inters(nint)%dat => inters(nint)%buf(nswmem)%dat
                end if

                do k=st(3),ed(3)
                do j=st(2),ed(2)
                do i=st(1),ed(1)
                    if (nswave == nsgl_aver_art) then
                        vol = one
                    else
                        vol = mb_vol(nbs)%fld(1)%r3d(i,j,k)
                    end if

                    do m=nst,ned
                        inters(nint)%dat(i,j,k,m) = mb_var(nbs)%fld(m)%r3d(i,j,k)/vol
                    end do
                end do
                end do
                end do

                index_msg = index_msg + 1
                call MPI_ISEND(inters(nint)%dat,packsize,kind_real_mpi, &
                              id_des,tag,MPI_COMM_WORLD,request_int(index_msg),ierr)

            end if

            if (myid == id_des) then
                if (nswmem == nbc_inter_buf_dyn) then
                    allocate(inters(nint)%dat(st(1):ed(1), &
                                              st(2):ed(2), &
                                              st(3):ed(3), &
                                              nst:ned),stat=ierr)
                else
                    inters(nint)%dat => inters(nint)%buf(nswmem)%dat
                end if

                index_msg = index_msg + 1
                call MPI_IRECV(inters(nint)%dat,packsize,kind_real_mpi, &
                              id_src,tag,MPI_COMM_WORLD,request_int(index_msg),ierr)

            end if

        else    ! i.e. (id_src == id_des)
            if (myid == id_src) then
#endif
                if (nswmem == nbc_inter_buf_dyn) then
                    allocate(inters(nint)%dat(st(1):ed(1), &
                                              st(2):ed(2), &
                                              st(3):ed(3), &
                                              nst:ned),stat=ierr)
                else
                    inters(nint)%dat => inters(nint)%buf(nswmem)%dat
                end if

                do k=st(3),ed(3)
                do j=st(2),ed(2)
                do i=st(1),ed(1)
                    if (nswave == nsgl_aver_art) then
                        vol = one
                    else
                        vol = mb_vol(nbs)%fld(1)%r3d(i,j,k)
                    end if

                    do m=nst,ned
                        inters(nint)%dat(i,j,k,m) = mb_var(nbs)%fld(m)%r3d(i,j,k)/vol
                    end do
                end do
                end do
                end do

                do k=t_st(3),t_ed(3)
                do j=t_st(2),t_ed(2)
                do i=t_st(1),t_ed(1)
                    ijkt(:) = (/i,j,k/)
                    ijks(:) = mb_top(nbt)%bcs(nrt)%mapijk(i,j,k,:)
                    ijkgt(:) = ijkt(:)
                    ijkgs(:) = ijks(:)

                    do n=1,ngh
                        ijkgs(s_nd) = ijks(s_nd) - n*s_lr
                        is = ijkgs(1)
                        js = ijkgs(2)
                        ks = ijkgs(3)

                        ijkgt(t_nd) = ijkt(t_nd) + n*t_lr
                        it = ijkgt(1)
                        jt = ijkgt(2)
                        kt = ijkgt(3)

                        do m=nst,ned
                            mb_var(nbt)%fld(m)%r3d(it,jt,kt) = inters(nint)%dat(is,js,ks,m)
                        end do
                    end do
                end do
                end do
                end do

                if (nswmem == nbc_inter_buf_dyn) then
                    deallocate(inters(nint)%dat,stat=ierr)
                else
                    nullify(inters(nint)%dat)
                end if
#ifdef PARALLEL
            end if
        end if
#endif
    end do


#ifdef PARALLEL
    if (nmsg_int > 0) then
        call MPI_WAITALL(nmsg_int, request_int(1:nmsg_int), status_int, ierr)

        do nint=1,ninters
            nbs     = inters(nint)%bc%nbs
            nbt     = inters(nint)%bc%nbt
            id_src = mb_top(nbs)%pid
            id_des = mb_top(nbt)%pid

            s_nd    = inters(nint)%bc%s_nd
            s_lr    = inters(nint)%bc%s_lr

            nrt     = inters(nint)%bc%nrt
            t_st(:) = inters(nint)%bc%t_st(:)
            t_ed(:) = inters(nint)%bc%t_ed(:)
            t_nd    = inters(nint)%bc%t_nd
            t_lr    = inters(nint)%bc%t_lr

            if (id_src /= id_des) then
                if (myid == id_des) then
                    do k=t_st(3),t_ed(3)
                    do j=t_st(2),t_ed(2)
                    do i=t_st(1),t_ed(1)
                        ijkt(:) = (/i,j,k/)
                        ijks(:) = mb_top(nbt)%bcs(nrt)%mapijk(i,j,k,:)
                        ijkgt(:) = ijkt(:)
                        ijkgs(:) = ijks(:)

                        do n=1,ngh
                            ijkgs(s_nd) = ijks(s_nd) - n*s_lr
                            is = ijkgs(1)
                            js = ijkgs(2)
                            ks = ijkgs(3)

                            ijkgt(t_nd) = ijkt(t_nd) + n*t_lr
                            it = ijkgt(1)
                            jt = ijkgt(2)
                            kt = ijkgt(3)

                            do m=nst,ned
                                mb_var(nbt)%fld(m)%r3d(it,jt,kt) = inters(nint)%dat(is,js,ks,m)
                            end do
                        end do
                    end do
                    end do
                    end do
                end if

                if (myid == id_src .or. myid == id_des) then
                    if (nswmem == nbc_inter_buf_dyn) then
                        deallocate(inters(nint)%dat,stat=ierr)
                    else
                        nullify(inters(nint)%dat)
                    end if
                end if

            end if
        end do

    end if

#endif

end subroutine exchange_bc_var

!>
!! @brief 交换边界流场信息函数
!! @details 以对接边界为基本单位，逐一交换流场信息
!! @param[in] mb_var  场变量指针
!! @param[in] nst     场变量起始索引
!! @param[in] ned     场变量终止索引
!! @param[in] ngh     边界交换信息层数
!! @param[in] nswmem  通讯内存分配方式
!! @param[in] nswave  边界信息平均方式
!! @author 王勇献
!! @remark 2012年04月05日 修改为非阻塞式通信（投递消息部分）
!>
subroutine pre_exchange_bc_var(mb_var,nst,ned,ngh,nswmem,nswave)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_inter_buf_dyn
    use mod_constants, only : one,nsgl_aver_art
    use mod_datatypes, only : var_block_t
    use mod_fieldvars, only : mb_top,mb_vol
    use mod_fieldvars, only : ninters,inters
    use mod_parallels
    implicit none
    type(var_block_t), pointer             :: mb_var(:)
    integer(kind_int),          intent(in) :: nst,ned,ngh
    integer(kind_int),          intent(in) :: nswmem
    integer(kind_int),          intent(in) :: nswave
    integer(kind_int) :: nint,ierr
    integer(kind_int) :: nbs,s_st(3),s_ed(3),s_nd,s_lr
    integer(kind_int) :: nbt,nrt,t_st(3),t_ed(3),t_nd,t_lr
    integer(kind_int) :: st(3),ed(3)
    integer(kind_int) :: ijks(3),ijkt(3),ijkgs(3),ijkgt(3)
    integer(kind_int) :: i,j,k,m,n,is,js,ks,it,jt,kt
    real(kind_real)   :: vol
#ifdef PARALLEL
    integer(kind_int) :: id_src,id_des,tag,packsize
    integer(kind_int) :: index_msg
#endif

#ifdef PARALLEL
    index_msg = 0
#endif

    do nint=1,ninters
        nbs     = inters(nint)%bc%nbs
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

        call bc_extend_inward(s_st,s_ed,s_nd,s_lr,ngh,st,ed)

#ifdef PARALLEL
        packsize = product(ed(:)-st(:)+1)*(ned-nst+1)

        id_src = mb_top(nbs)%pid
        id_des = mb_top(nbt)%pid
        tag    = nint + ninters*3


        if (id_src /= id_des) then
            if (myid == id_src) then
                if (nswmem == nbc_inter_buf_dyn) then
                    allocate(inters(nint)%dat(st(1):ed(1), &
                                              st(2):ed(2), &
                                              st(3):ed(3), &
                                              nst:ned),stat=ierr)
                else
                    inters(nint)%dat => inters(nint)%buf(nswmem)%dat
                end if

                do k=st(3),ed(3)
                do j=st(2),ed(2)
                do i=st(1),ed(1)
                    if (nswave == nsgl_aver_art) then
                        vol = one
                    else
                        vol = mb_vol(nbs)%fld(1)%r3d(i,j,k)
                    end if

                    do m=nst,ned
                        inters(nint)%dat(i,j,k,m) = mb_var(nbs)%fld(m)%r3d(i,j,k)/vol
                    end do
                end do
                end do
                end do

                index_msg = index_msg + 1
                call MPI_ISEND(inters(nint)%dat,packsize,kind_real_mpi, &
                              id_des,tag,MPI_COMM_WORLD,request_int(index_msg),ierr)

            end if

            if (myid == id_des) then
                if (nswmem == nbc_inter_buf_dyn) then
                    allocate(inters(nint)%dat(st(1):ed(1), &
                                              st(2):ed(2), &
                                              st(3):ed(3), &
                                              nst:ned),stat=ierr)
                else
                    inters(nint)%dat => inters(nint)%buf(nswmem)%dat
                end if

                index_msg = index_msg + 1
                call MPI_IRECV(inters(nint)%dat,packsize,kind_real_mpi, &
                              id_src,tag,MPI_COMM_WORLD,request_int(index_msg),ierr)

            end if

        else    ! i.e. (id_src == id_des)
            if (myid == id_src) then
#endif
                if (nswmem == nbc_inter_buf_dyn) then
                    allocate(inters(nint)%dat(st(1):ed(1), &
                                              st(2):ed(2), &
                                              st(3):ed(3), &
                                              nst:ned),stat=ierr)
                else
                    inters(nint)%dat => inters(nint)%buf(nswmem)%dat
                end if

                do k=st(3),ed(3)
                do j=st(2),ed(2)
                do i=st(1),ed(1)
                    if (nswave == nsgl_aver_art) then
                        vol = one
                    else
                        vol = mb_vol(nbs)%fld(1)%r3d(i,j,k)
                    end if

                    do m=nst,ned
                        inters(nint)%dat(i,j,k,m) = mb_var(nbs)%fld(m)%r3d(i,j,k)/vol
                    end do
                end do
                end do
                end do

                do k=t_st(3),t_ed(3)
                do j=t_st(2),t_ed(2)
                do i=t_st(1),t_ed(1)
                    ijkt(:) = (/i,j,k/)
                    ijks(:) = mb_top(nbt)%bcs(nrt)%mapijk(i,j,k,:)
                    ijkgt(:) = ijkt(:)
                    ijkgs(:) = ijks(:)

                    do n=1,ngh
                        ijkgs(s_nd) = ijks(s_nd) - n*s_lr
                        is = ijkgs(1)
                        js = ijkgs(2)
                        ks = ijkgs(3)

                        ijkgt(t_nd) = ijkt(t_nd) + n*t_lr
                        it = ijkgt(1)
                        jt = ijkgt(2)
                        kt = ijkgt(3)

                        do m=nst,ned
                            mb_var(nbt)%fld(m)%r3d(it,jt,kt) = inters(nint)%dat(is,js,ks,m)
                        end do
                    end do
                end do
                end do
                end do

                if (nswmem == nbc_inter_buf_dyn) then
                    deallocate(inters(nint)%dat,stat=ierr)
                else
                    nullify(inters(nint)%dat)
                end if
#ifdef PARALLEL
            end if
        end if
#endif
    end do

end subroutine pre_exchange_bc_var

subroutine pre_exchange_cc_bc_var(mb_var,nst,ned,ngh,nswmem,nswave)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_inter_buf_dyn
    use mod_constants, only : one,nsgl_aver_art
    use mod_datatypes, only : var_block_t
    use mod_fieldvars, only : mb_topc,mb_volcc
    use mod_fieldvars, only : ninters,interscc
    use mod_parallels
    implicit none
    type(var_block_t), pointer             :: mb_var(:)
    integer(kind_int),          intent(in) :: nst,ned,ngh
    integer(kind_int),          intent(in) :: nswmem
    integer(kind_int),          intent(in) :: nswave
    integer(kind_int) :: nint,ierr
    integer(kind_int) :: nbs,s_st(3),s_ed(3),s_nd,s_lr
    integer(kind_int) :: nbt,nrt,t_st(3),t_ed(3),t_nd,t_lr
    integer(kind_int) :: st(3),ed(3)
    integer(kind_int) :: ijks(3),ijkt(3),ijkgs(3),ijkgt(3)
    integer(kind_int) :: i,j,k,m,n,is,js,ks,it,jt,kt
    real(kind_real)   :: vol
#ifdef PARALLEL
    integer(kind_int) :: id_src,id_des,tag,packsize
    integer(kind_int) :: index_msg
#endif

#ifdef PARALLEL
    index_msg = 0
#endif

    do nint=1,ninters
        nbs     = interscc(nint)%bc%nbs
        s_st(:) = interscc(nint)%bc%s_st(:)
        s_ed(:) = interscc(nint)%bc%s_ed(:)
        s_nd    = interscc(nint)%bc%s_nd
        s_lr    = interscc(nint)%bc%s_lr

        nbt     = interscc(nint)%bc%nbt
        nrt     = interscc(nint)%bc%nrt
        t_st(:) = interscc(nint)%bc%t_st(:)
        t_ed(:) = interscc(nint)%bc%t_ed(:)
        t_nd    = interscc(nint)%bc%t_nd
        t_lr    = interscc(nint)%bc%t_lr

        call bc_extend_inward(s_st,s_ed,s_nd,s_lr,ngh,st,ed)

#ifdef PARALLEL
        packsize = product(ed(:)-st(:)+1)*(ned-nst+1)

        id_src = mb_topc(nbs)%pid
        id_des = mb_topc(nbt)%pid
        tag    = nint + ninters*3


        if (id_src /= id_des) then
            if (myid == id_src) then
                if (nswmem == nbc_inter_buf_dyn) then
                    allocate(interscc(nint)%dat(st(1):ed(1), &
                                                st(2):ed(2), &
                                                st(3):ed(3), &
                                                nst:ned),stat=ierr)
                else
                    interscc(nint)%dat => interscc(nint)%buf(nswmem)%dat
                end if

                do k=st(3),ed(3)
                do j=st(2),ed(2)
                do i=st(1),ed(1)
                    if (nswave == nsgl_aver_art) then
                        vol = one
                    else
                        vol = mb_volcc(nbs)%fld(1)%r3d(i,j,k)
                    end if

                    do m=nst,ned
                        interscc(nint)%dat(i,j,k,m) = mb_var(nbs)%fld(m)%r3d(i,j,k)/vol
                    end do
                end do
                end do
                end do

                index_msg = index_msg + 1
                call MPI_ISEND(interscc(nint)%dat,packsize,kind_real_mpi, &
                              id_des,tag,MPI_COMM_WORLD,request_int(index_msg),ierr)

            end if

            if (myid == id_des) then
                if (nswmem == nbc_inter_buf_dyn) then
                    allocate(interscc(nint)%dat(st(1):ed(1), &
                                              st(2):ed(2), &
                                              st(3):ed(3), &
                                              nst:ned),stat=ierr)
                else
                    interscc(nint)%dat => interscc(nint)%buf(nswmem)%dat
                end if

                index_msg = index_msg + 1
                call MPI_IRECV(interscc(nint)%dat,packsize,kind_real_mpi, &
                              id_src,tag,MPI_COMM_WORLD,request_int(index_msg),ierr)

            end if

        else    ! i.e. (id_src == id_des)
            if (myid == id_src) then
#endif
                if (nswmem == nbc_inter_buf_dyn) then
                    allocate(interscc(nint)%dat(st(1):ed(1), &
                                                st(2):ed(2), &
                                                st(3):ed(3), &
                                                nst:ned),stat=ierr)
                else
                    interscc(nint)%dat => interscc(nint)%buf(nswmem)%dat
                end if

                do k=st(3),ed(3)
                do j=st(2),ed(2)
                do i=st(1),ed(1)
                    if (nswave == nsgl_aver_art) then
                        vol = one
                    else
                        vol = mb_volcc(nbs)%fld(1)%r3d(i,j,k)
                    end if

                    do m=nst,ned
                        interscc(nint)%dat(i,j,k,m) = mb_var(nbs)%fld(m)%r3d(i,j,k)/vol
                    end do
                end do
                end do
                end do

                do k=t_st(3),t_ed(3)
                do j=t_st(2),t_ed(2)
                do i=t_st(1),t_ed(1)
                    ijkt(:) = (/i,j,k/)
                    ijks(:) = mb_topc(nbt)%bcs(nrt)%mapijk(i,j,k,:)
                    ijkgt(:) = ijkt(:)
                    ijkgs(:) = ijks(:)

                    do n=0,ngh
                        ijkgs(s_nd) = ijks(s_nd) - n*s_lr
                        is = ijkgs(1)
                        js = ijkgs(2)
                        ks = ijkgs(3)

                        ijkgt(t_nd) = ijkt(t_nd) + n*t_lr
                        it = ijkgt(1)
                        jt = ijkgt(2)
                        kt = ijkgt(3)

                        do m=nst,ned
                            mb_var(nbt)%fld(m)%r3d(it,jt,kt) = interscc(nint)%dat(is,js,ks,m)
                        end do
                    end do
                end do
                end do
                end do

                if (nswmem == nbc_inter_buf_dyn) then
                    deallocate(interscc(nint)%dat,stat=ierr)
                else
                    nullify(interscc(nint)%dat)
                end if
#ifdef PARALLEL
            end if
        end if
#endif
    end do

end subroutine pre_exchange_cc_bc_var

!>
!! @brief 交换边界流场信息函数
!! @details 以对接边界为基本单位，逐一交换流场信息
!! @param[in] mb_var  场变量指针
!! @param[in] nst     场变量起始索引
!! @param[in] ned     场变量终止索引
!! @param[in] ngh     边界交换信息层数
!! @param[in] nswmem  通讯内存分配方式
!! @param[in] nswave  边界信息平均方式
!! @author 王勇献
!! @remark 2012年04月05日 修改为非阻塞式通信（投递消息部分）
!>
subroutine pre_exchange_bc_var_sp(mb_var,nst,ned,ngh,nswmem,nswave)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_inter_buf_dyn
    use mod_constants, only : one,nsgl_aver_art
    use mod_datatypes, only : var_block_t
    use mod_fieldvars, only : mb_topsp,mb_volsp
    use mod_fieldvars, only : ninters,interssp
    use mod_parallels
    implicit none
    type(var_block_t), pointer             :: mb_var(:)
    integer(kind_int),          intent(in) :: nst,ned,ngh
    integer(kind_int),          intent(in) :: nswmem
    integer(kind_int),          intent(in) :: nswave
    integer(kind_int) :: nint,ierr
    integer(kind_int) :: nbs,s_st(3),s_ed(3),s_nd,s_lr
    integer(kind_int) :: nbt,nrt,t_st(3),t_ed(3),t_nd,t_lr
    integer(kind_int) :: st(3),ed(3),s_nd3d(3),t_nd3d(3)
    integer(kind_int) :: ijks(3),ijkt(3),ijkgs(3),ijkgt(3)
    integer(kind_int) :: i,j,k,m,n,is,js,ks,it,jt,kt
    real(kind_real)   :: vol
#ifdef PARALLEL
    integer(kind_int) :: id_src,id_des,tag,packsize
    integer(kind_int) :: index_msg
#endif

#ifdef PARALLEL
    index_msg = 0
#endif

    do nint=1,ninters
        nbs     = interssp(nint)%bc%nbs
        s_st(:) = interssp(nint)%bc%s_st(:)
        s_ed(:) = interssp(nint)%bc%s_ed(:)
        s_nd    = interssp(nint)%bc%s_nd
        s_lr    = interssp(nint)%bc%s_lr

        nbt     = interssp(nint)%bc%nbt
        nrt     = interssp(nint)%bc%nrt
        t_st(:) = interssp(nint)%bc%t_st(:)
        t_ed(:) = interssp(nint)%bc%t_ed(:)
        t_nd    = interssp(nint)%bc%t_nd
        t_lr    = interssp(nint)%bc%t_lr

        call bc_extend_inward(s_st,s_ed,s_nd,s_lr,ngh,st,ed)

#ifdef PARALLEL
        packsize = product(ed(:)-st(:)+1)*(ned-nst+1)

        id_src = mb_topsp(nbs)%pid
        id_des = mb_topsp(nbt)%pid
        tag    = nint + ninters*3


        if (id_src /= id_des) then
            if (myid == id_src) then
                if (nswmem == nbc_inter_buf_dyn) then
                    allocate(interssp(nint)%dat(st(1):ed(1), &
                                                st(2):ed(2), &
                                                st(3):ed(3), &
                                                nst:ned),stat=ierr)
                else
                    interssp(nint)%dat => interssp(nint)%buf(nswmem)%dat
                end if

                do k=st(3),ed(3)
                do j=st(2),ed(2)
                do i=st(1),ed(1)
                    if (nswave == nsgl_aver_art) then
                        vol = one
                    else
                        vol = mb_volsp(nbs)%fld(1)%r3d(i,j,k)
                    end if

                    do m=nst,ned
                        interssp(nint)%dat(i,j,k,m) = mb_var(nbs)%fld(m)%r3d(i,j,k)/vol
                    end do
                end do
                end do
                end do

                index_msg = index_msg + 1
                call MPI_ISEND(interssp(nint)%dat,packsize,kind_real_mpi, &
                              id_des,tag,MPI_COMM_WORLD,request_int(index_msg),ierr)

            end if

            if (myid == id_des) then
                if (nswmem == nbc_inter_buf_dyn) then
                    allocate(interssp(nint)%dat(st(1):ed(1), &
                                                st(2):ed(2), &
                                                st(3):ed(3), &
                                                nst:ned),stat=ierr)
                else
                    interssp(nint)%dat => interssp(nint)%buf(nswmem)%dat
                end if

                index_msg = index_msg + 1
                call MPI_IRECV(interssp(nint)%dat,packsize,kind_real_mpi, &
                              id_src,tag,MPI_COMM_WORLD,request_int(index_msg),ierr)

            end if

        else    ! i.e. (id_src == id_des)
            if (myid == id_src) then
#endif
                if (nswmem == nbc_inter_buf_dyn) then
                    allocate(interssp(nint)%dat(st(1):ed(1), &
                                              st(2):ed(2), &
                                              st(3):ed(3), &
                                              nst:ned),stat=ierr)
                else
                    interssp(nint)%dat => interssp(nint)%buf(nswmem)%dat
                end if

                do k=st(3),ed(3)
                do j=st(2),ed(2)
                do i=st(1),ed(1)
                    if (nswave == nsgl_aver_art) then
                        vol = one
                    else
                        vol = mb_volsp(nbs)%fld(1)%r3d(i,j,k)
                    end if

                    do m=nst,ned
                        interssp(nint)%dat(i,j,k,m) = mb_var(nbs)%fld(m)%r3d(i,j,k)/vol
                    end do
                end do
                end do
                end do

                do k=t_st(3),t_ed(3)
                do j=t_st(2),t_ed(2)
                do i=t_st(1),t_ed(1)
                    ijkt(:) = (/i,j,k/)
                    ijks(:) = mb_topsp(nbt)%bcs(nrt)%mapijk(i,j,k,:)
                    ijkgt(:) = ijkt(:)
                    ijkgs(:) = ijks(:)

                    do n=1,ngh
                        ijkgs(s_nd) = ijks(s_nd) - (n-1)*s_lr
                        is = ijkgs(1)
                        js = ijkgs(2)
                        ks = ijkgs(3)

                        ijkgt(t_nd) = ijkt(t_nd) + n*t_lr
                        it = ijkgt(1)
                        jt = ijkgt(2)
                        kt = ijkgt(3)

                        do m=nst,ned
                            mb_var(nbt)%fld(m)%r3d(it,jt,kt) = interssp(nint)%dat(is,js,ks,m)
                        end do
                    end do
                end do
                end do
                end do

                if (nswmem == nbc_inter_buf_dyn) then
                    deallocate(interssp(nint)%dat,stat=ierr)
                else
                    nullify(interssp(nint)%dat)
                end if
#ifdef PARALLEL
            end if
        end if
#endif
    end do

end subroutine pre_exchange_bc_var_sp

!>
!! @brief 交换边界流场信息函数
!! @details 以对接边界为基本单位，逐一交换流场信息
!! @param[in] mb_var  场变量指针
!! @param[in] nst     场变量起始索引
!! @param[in] ned     场变量终止索引
!! @param[in] ngh     边界交换信息层数
!! @param[in] nswmem  通讯内存分配方式
!! @param[in] nswave  边界信息平均方式
!! @author 王勇献
!! @remark 2012年04月05日 修改为非阻塞式通信（接受消息后的处理）
!>
subroutine post_exchange_bc_var(mb_var,nst,ned,ngh,nswmem,nswave)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_inter_buf_dyn
    use mod_constants, only : one,nsgl_aver_art
    use mod_datatypes, only : var_block_t
    use mod_fieldvars, only : mb_top,mb_vol
    use mod_fieldvars, only : ninters,inters
    use mod_parallels
    implicit none
    type(var_block_t), pointer, intent(in) :: mb_var(:)
    integer(kind_int),          intent(in) :: nst,ned,ngh
    integer(kind_int),          intent(in) :: nswmem
    integer(kind_int),          intent(in) :: nswave
    integer(kind_int) :: nint,ierr
    integer(kind_int) :: nbs,s_nd,s_lr
    integer(kind_int) :: nbt,nrt,t_st(3),t_ed(3),t_nd,t_lr
    integer(kind_int) :: ijks(3),ijkt(3),ijkgs(3),ijkgt(3)
    integer(kind_int) :: i,j,k,m,n,is,js,ks,it,jt,kt
#ifdef PARALLEL
    integer(kind_int) :: id_src,id_des
#endif

#ifdef PARALLEL
    if (nmsg_int > 0) then
        call MPI_WAITALL(nmsg_int, request_int(1:nmsg_int), status_int, ierr)

        do nint=1,ninters
            nbs     = inters(nint)%bc%nbs
            nbt     = inters(nint)%bc%nbt
            id_src = mb_top(nbs)%pid
            id_des = mb_top(nbt)%pid

            s_nd    = inters(nint)%bc%s_nd
            s_lr    = inters(nint)%bc%s_lr

            nrt     = inters(nint)%bc%nrt
            t_st(:) = inters(nint)%bc%t_st(:)
            t_ed(:) = inters(nint)%bc%t_ed(:)
            t_nd    = inters(nint)%bc%t_nd
            t_lr    = inters(nint)%bc%t_lr

            if (id_src /= id_des) then
                if (myid == id_des) then
                    do k=t_st(3),t_ed(3)
                    do j=t_st(2),t_ed(2)
                    do i=t_st(1),t_ed(1)
                        ijkt(:) = (/i,j,k/)                        
                        ijks(:) = mb_top(nbt)%bcs(nrt)%mapijk(i,j,k,:)
                        ijkgt(:) = ijkt(:)
                        ijkgs(:) = ijks(:)

                        do n=1,ngh
                            ijkgs(s_nd) = ijks(s_nd) - n*s_lr
                            is = ijkgs(1)
                            js = ijkgs(2)
                            ks = ijkgs(3)

                            ijkgt(t_nd) = ijkt(t_nd) + n*t_lr
                            it = ijkgt(1)
                            jt = ijkgt(2)
                            kt = ijkgt(3)

                            do m=nst,ned
                                mb_var(nbt)%fld(m)%r3d(it,jt,kt) = inters(nint)%dat(is,js,ks,m)
                            end do
                        end do
                    end do
                    end do
                    end do
                end if

                if (myid == id_src .or. myid == id_des) then
                    if (nswmem == nbc_inter_buf_dyn) then
                        deallocate(inters(nint)%dat,stat=ierr)
                    else
                        nullify(inters(nint)%dat)
                    end if
                end if

            end if
        end do

    end if

#endif

end subroutine post_exchange_bc_var

subroutine post_exchange_cc_bc_var(mb_var,nst,ned,ngh,nswmem,nswave)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_inter_buf_dyn
    use mod_constants, only : one,nsgl_aver_art
    use mod_datatypes, only : var_block_t
    use mod_fieldvars, only : mb_topc,mb_volcc
    use mod_fieldvars, only : ninters,interscc
    use mod_parallels
    implicit none
    type(var_block_t), pointer, intent(in) :: mb_var(:)
    integer(kind_int),          intent(in) :: nst,ned,ngh
    integer(kind_int),          intent(in) :: nswmem
    integer(kind_int),          intent(in) :: nswave
    integer(kind_int) :: nint,ierr
    integer(kind_int) :: nbs,s_nd,s_lr
    integer(kind_int) :: nbt,nrt,t_st(3),t_ed(3),t_nd,t_lr
    integer(kind_int) :: ijks(3),ijkt(3),ijkgs(3),ijkgt(3)
    integer(kind_int) :: i,j,k,m,n,is,js,ks,it,jt,kt
#ifdef PARALLEL
    integer(kind_int) :: id_src,id_des
#endif

#ifdef PARALLEL
    if (nmsg_int > 0) then
        call MPI_WAITALL(nmsg_int, request_int(1:nmsg_int), status_int, ierr)

        do nint=1,ninters
            nbs     = interscc(nint)%bc%nbs
            nbt     = interscc(nint)%bc%nbt
            id_src = mb_topc(nbs)%pid
            id_des = mb_topc(nbt)%pid

            s_nd    = interscc(nint)%bc%s_nd
            s_lr    = interscc(nint)%bc%s_lr

            nrt     = interscc(nint)%bc%nrt
            t_st(:) = interscc(nint)%bc%t_st(:)
            t_ed(:) = interscc(nint)%bc%t_ed(:)
            t_nd    = interscc(nint)%bc%t_nd
            t_lr    = interscc(nint)%bc%t_lr

            if (id_src /= id_des) then
                if (myid == id_des) then
                    do k=t_st(3),t_ed(3)
                    do j=t_st(2),t_ed(2)
                    do i=t_st(1),t_ed(1)
                        ijkt(:) = (/i,j,k/)                        
                        ijks(:) = mb_topc(nbt)%bcs(nrt)%mapijk(i,j,k,:)
                        ijkgt(:) = ijkt(:)
                        ijkgs(:) = ijks(:)

                        do n=0,ngh
                            ijkgs(s_nd) = ijks(s_nd) - n*s_lr
                            is = ijkgs(1)
                            js = ijkgs(2)
                            ks = ijkgs(3)

                            ijkgt(t_nd) = ijkt(t_nd) + n*t_lr
                            it = ijkgt(1)
                            jt = ijkgt(2)
                            kt = ijkgt(3)

                            do m=nst,ned
                                mb_var(nbt)%fld(m)%r3d(it,jt,kt) = interscc(nint)%dat(is,js,ks,m)
                            end do
                        end do
                    end do
                    end do
                    end do
                end if

                if (myid == id_src .or. myid == id_des) then
                    if (nswmem == nbc_inter_buf_dyn) then
                        deallocate(interscc(nint)%dat,stat=ierr)
                    else
                        nullify(interscc(nint)%dat)
                    end if
                end if

            end if
        end do

    end if

#endif

end subroutine post_exchange_cc_bc_var

subroutine post_exchange_bc_var_sp(mb_var,nst,ned,ngh,nswmem,nswave)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_inter_buf_dyn
    use mod_constants, only : one,nsgl_aver_art
    use mod_datatypes, only : var_block_t
    use mod_fieldvars, only : mb_topsp,mb_volsp
    use mod_fieldvars, only : ninters,interssp
    use mod_parallels
    implicit none
    type(var_block_t), pointer, intent(in) :: mb_var(:)
    integer(kind_int),          intent(in) :: nst,ned,ngh
    integer(kind_int),          intent(in) :: nswmem
    integer(kind_int),          intent(in) :: nswave
    integer(kind_int) :: nint,ierr
    integer(kind_int) :: nbs,s_nd,s_lr
    integer(kind_int) :: nbt,nrt,t_st(3),t_ed(3),t_nd,t_lr
    integer(kind_int) :: ijks(3),ijkt(3),ijkgs(3),ijkgt(3)
    integer(kind_int) :: i,j,k,m,n,is,js,ks,it,jt,kt
#ifdef PARALLEL
    integer(kind_int) :: id_src,id_des
#endif

#ifdef PARALLEL
    if (nmsg_int > 0) then
        call MPI_WAITALL(nmsg_int, request_int(1:nmsg_int), status_int, ierr)

        do nint=1,ninters
            nbs     = interssp(nint)%bc%nbs
            nbt     = interssp(nint)%bc%nbt
            id_src = mb_topsp(nbs)%pid
            id_des = mb_topsp(nbt)%pid

            s_nd    = interssp(nint)%bc%s_nd
            s_lr    = interssp(nint)%bc%s_lr

            nrt     = interssp(nint)%bc%nrt
            t_st(:) = interssp(nint)%bc%t_st(:)
            t_ed(:) = interssp(nint)%bc%t_ed(:)
            t_nd    = interssp(nint)%bc%t_nd
            t_lr    = interssp(nint)%bc%t_lr

            if (id_src /= id_des) then
                if (myid == id_des) then
                    do k=t_st(3),t_ed(3)
                    do j=t_st(2),t_ed(2)
                    do i=t_st(1),t_ed(1)
                        ijkt(:) = (/i,j,k/)                        
                        ijks(:) = mb_topsp(nbt)%bcs(nrt)%mapijk(i,j,k,:)
                        ijkgt(:) = ijkt(:)
                        ijkgs(:) = ijks(:)

                        do n=1,ngh
                            ijkgs(s_nd) = ijks(s_nd) - (n-1)*s_lr
                            is = ijkgs(1)
                            js = ijkgs(2)
                            ks = ijkgs(3)

                            ijkgt(t_nd) = ijkt(t_nd) + n*t_lr
                            it = ijkgt(1)
                            jt = ijkgt(2)
                            kt = ijkgt(3)

                            do m=nst,ned
                                mb_var(nbt)%fld(m)%r3d(it,jt,kt) = interssp(nint)%dat(is,js,ks,m)
                            end do
                        end do
                    end do
                    end do
                    end do
                end if

                if (myid == id_src .or. myid == id_des) then
                    if (nswmem == nbc_inter_buf_dyn) then
                        deallocate(interssp(nint)%dat,stat=ierr)
                    else
                        nullify(interssp(nint)%dat)
                    end if
                end if

            end if
        end do

    end if

#endif

end subroutine post_exchange_bc_var_sp

subroutine average_bc_var(mb_var,nst,ned,nswmem,nswave)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_inter_buf_dyn,nbc_inter_buf_max
    use mod_constants, only : half,one,nsgl_aver_art
    use mod_datatypes, only : var_block_t
    use mod_fieldvars, only : mb_top,mb_vol
    use mod_fieldvars, only : ninterf,inters
    use mod_parallels
    implicit none
    type(var_block_t), pointer, intent(in) :: mb_var(:)
    integer(kind_int),          intent(in) :: nst,ned
    integer(kind_int),          intent(in) :: nswmem
    integer(kind_int),          intent(in) :: nswave
    integer(kind_int) :: nint,ierr
    integer(kind_int) :: nbs,nrs,s_st(3),s_ed(3),s_nd,s_lr
    integer(kind_int) :: nbt,nrt,t_st(3),t_ed(3),t_nd,t_lr
    integer(kind_int) :: ijks(3),ijkt(3),id_src,id_des
    integer(kind_int) :: i,j,k,m,is,js,ks,it,jt,kt
    real(kind_real)   :: vol

    if (nswmem >  nbc_inter_buf_dyn .and. &
        nswmem <= nbc_inter_buf_max) then
        ierr = 0
    else
        ierr = 1
    end if
    call error_check(ierr, &
                     "nsw isn't in the range in subroutine average_bc_var")

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
            inters(nint)%dat => inters(nint)%buf(nswmem)%dat

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

                if (nswave == nsgl_aver_art) then
                    vol = one
                else
                    vol = mb_vol(nbt)%fld(1)%r3d(it,jt,kt)
                end if

                do m=nst,ned
                    mb_var(nbt)%fld(m)%r3d(it,jt,kt) = &
                        half*(mb_var(nbt)%fld(m)%r3d(it,jt,kt) + &
                              vol*inters(nint)%dat(is,js,ks,m))
                end do
            end do
            end do
            end do

            nullify(inters(nint)%dat)
#ifdef PARALLEL
        end if
#endif
    end do

end subroutine average_bc_var

subroutine average_bc_cc_var(mb_var,nst,ned,nswmem,nswave)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nbc_inter_buf_dyn,nbc_inter_buf_max
    use mod_constants, only : half,one,nsgl_aver_art
    use mod_datatypes, only : var_block_t
    use mod_fieldvars, only : mb_topc,mb_volcc
    use mod_fieldvars, only : ninterf,interscc
    use mod_parallels
    implicit none
    type(var_block_t), pointer, intent(in) :: mb_var(:)
    integer(kind_int),          intent(in) :: nst,ned
    integer(kind_int),          intent(in) :: nswmem
    integer(kind_int),          intent(in) :: nswave
    integer(kind_int) :: nint,ierr
    integer(kind_int) :: nbs,nrs,s_st(3),s_ed(3),s_nd,s_lr
    integer(kind_int) :: nbt,nrt,t_st(3),t_ed(3),t_nd,t_lr
    integer(kind_int) :: ijks(3),ijkt(3),id_src,id_des
    integer(kind_int) :: i,j,k,m,is,js,ks,it,jt,kt
    real(kind_real)   :: vol

    if (nswmem >  nbc_inter_buf_dyn .and. &
        nswmem <= nbc_inter_buf_max) then
        ierr = 0
    else
        ierr = 1
    end if
    call error_check(ierr, &
                     "nsw isn't in the range in subroutine average_bc_var")

    do nint=1,ninterf
        nbs     = interscc(nint)%bc%nbs
        nrs     = interscc(nint)%bc%nrs
        s_st(:) = interscc(nint)%bc%s_st(:)
        s_ed(:) = interscc(nint)%bc%s_ed(:)
        s_nd    = interscc(nint)%bc%s_nd
        s_lr    = interscc(nint)%bc%s_lr

        nbt     = interscc(nint)%bc%nbt
        nrt     = interscc(nint)%bc%nrt
        t_st(:) = interscc(nint)%bc%t_st(:)
        t_ed(:) = interscc(nint)%bc%t_ed(:)
        t_nd    = interscc(nint)%bc%t_nd
        t_lr    = interscc(nint)%bc%t_lr

#ifdef PARALLEL
        id_src = mb_topc(nbs)%pid
        id_des = mb_topc(nbt)%pid

        if (myid == id_des) then
#endif
            interscc(nint)%dat => interscc(nint)%buf(nswmem)%dat

            do k=t_st(3),t_ed(3)
            do j=t_st(2),t_ed(2)
            do i=t_st(1),t_ed(1)
                ijkt(:) = (/i,j,k/)
                ijks(:) = mb_topc(nbt)%bcs(nrt)%mapijk(i,j,k,:)

                is = ijks(1)
                js = ijks(2)
                ks = ijks(3)
                it = ijkt(1)
                jt = ijkt(2)
                kt = ijkt(3)

                if (nswave == nsgl_aver_art) then
                    vol = one
                else
                    vol = mb_volcc(nbt)%fld(1)%r3d(it,jt,kt)
                end if

                do m=nst,ned
                    mb_var(nbt)%fld(m)%r3d(it,jt,kt) = &
                        half*(mb_var(nbt)%fld(m)%r3d(it,jt,kt) + &
                              vol*interscc(nint)%dat(is,js,ks,m))
                end do
            end do
            end do
            end do

            nullify(interscc(nint)%dat)
#ifdef PARALLEL
        end if
#endif
    end do

end subroutine average_bc_cc_var

!>
!! @brief 交换边界流场梯度函数
!! @details 以对接边界为基本单位，逐一交换流场梯度
!! @param[in] mb_var  场变量指针
!! @param[in] nst     场变量起始索引
!! @param[in] ned     场变量终止索引
!! @param[in] ngh     边界交换信息层数
!! @param[in] nswmem  通讯内存分配方式
!! @author 刘化勇
!! @date 2011年05月31日
!! @author 王勇献
!! @remark 2012年02月14日 修改为非阻塞式通信
!>
subroutine exchange_bc_der(mb_der,nst,ned,ngh,nswmem)
    use mod_kndconsts, only : kind_int
    use mod_constants, only : nbc_inter_buf_dyn
    use mod_datatypes, only : var_block_t
    use mod_fieldvars, only : nblocks,mb_top,ninters,inters
    use mod_parallels
    implicit none
    type(var_block_t), pointer, intent(in) :: mb_der(:)
    integer(kind_int),          intent(in) :: nst,ned
    integer(kind_int),          intent(in) :: ngh
    integer(kind_int),          intent(in) :: nswmem
    integer(kind_int) :: nint,ierr
    integer(kind_int) :: nbs,nrs,s_st(3),s_ed(3),s_nd,s_lr
    integer(kind_int) :: nbt,nrt,t_st(3),t_ed(3),t_nd,t_lr
    integer(kind_int) :: s_ord(3),s_sgn(3),t_ord(3),t_sgn(3)
    integer(kind_int) :: st(3),ed(3),sgn(nst:ned),ord(nst:ned)
    integer(kind_int) :: ijks(3),ijkt(3),ijkgs(3),ijkgt(3)
    integer(kind_int) :: i,j,k,m,n,is,js,ks,it,jt,kt,sd,td
#ifdef PARALLEL
    integer(kind_int) :: id_src,id_des,tag,packsize
    integer(kind_int) :: index_msg
#endif

    ierr = mod(ned-nst+1,3)
    call error_check(ierr, &
                     "The size of array isn't times of 3 in subroutine exchange_bc_der")

#ifdef PARALLEL
    index_msg = 0
#endif

    do nint=1,ninters
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

        s_ord(:) = mb_top(nbt)%bcs(nrt)%s_ord(:)
        s_sgn(:) = mb_top(nbt)%bcs(nrt)%s_sgn(:)
        t_ord(:) = mb_top(nbt)%bcs(nrt)%t_ord(:)
        t_sgn(:) = mb_top(nbt)%bcs(nrt)%t_sgn(:)

        do m=nst,ned
           sd = m-nst+1 - (m-nst)/3*3
           ord(m) = t_ord(sd) + (m-nst)/3*3
           sgn(m) = t_sgn(sd)
        end do

        call bc_extend_inward(s_st,s_ed,s_nd,s_lr,ngh,st,ed)

#ifdef PARALLEL
        packsize = product(ed(:)-st(:)+1)*(ned-nst+1)

        id_src = mb_top(nbs)%pid
        id_des = mb_top(nbt)%pid
        tag    = nint + ninters*3


        if (id_src /= id_des) then
            if (myid == id_src) then
                if (nswmem == nbc_inter_buf_dyn) then
                    allocate(inters(nint)%dat(st(1):ed(1), &
                                              st(2):ed(2), &
                                              st(3):ed(3), &
                                              nst:ned),stat=ierr)
                else
                    inters(nint)%dat => inters(nint)%buf(nswmem)%dat
                end if

                do m=nst,ned
                    do k=st(3),ed(3)
                    do j=st(2),ed(2)
                    do i=st(1),ed(1)
                        inters(nint)%dat(i,j,k,m) = mb_der(nbs)%fld(m)%r3d(i,j,k)
                    end do
                    end do
                    end do
                end do

                index_msg = index_msg + 1
                call MPI_ISEND(inters(nint)%dat,packsize,kind_real_mpi, &
                              id_des,tag,MPI_COMM_WORLD,request_int(index_msg),ierr)

            end if

            if (myid == id_des) then
                if (nswmem == nbc_inter_buf_dyn) then
                    allocate(inters(nint)%dat(st(1):ed(1), &
                                              st(2):ed(2), &
                                              st(3):ed(3), &
                                              nst:ned),stat=ierr)
                else
                    inters(nint)%dat => inters(nint)%buf(nswmem)%dat
                end if

                index_msg = index_msg + 1
                call MPI_IRECV(inters(nint)%dat,packsize,kind_real_mpi, &
                              id_src,tag,MPI_COMM_WORLD,request_int(index_msg),ierr)

            end if

        else
            if (myid == id_src) then
#endif
                if (nswmem == nbc_inter_buf_dyn) then
                    allocate(inters(nint)%dat(st(1):ed(1), &
                                              st(2):ed(2), &
                                              st(3):ed(3), &
                                              nst:ned),stat=ierr)
                else
                    inters(nint)%dat => inters(nint)%buf(nswmem)%dat
                end if

                do m=nst,ned
                    do k=st(3),ed(3)
                    do j=st(2),ed(2)
                    do i=st(1),ed(1)
                        inters(nint)%dat(i,j,k,m) = mb_der(nbs)%fld(m)%r3d(i,j,k)
                    end do
                    end do
                    end do
                end do

                do k=t_st(3),t_ed(3)
                do j=t_st(2),t_ed(2)
                do i=t_st(1),t_ed(1)
                    ijkt(:) = (/i,j,k/)
                    ijks(:) = mb_top(nbt)%bcs(nrt)%mapijk(i,j,k,:)
                    ijkgt(:) = ijkt(:)
                    ijkgs(:) = ijks(:)
                    do n=1,ngh
                        ijkgs(s_nd) = ijks(s_nd) - n*s_lr
                        is = ijkgs(1)
                        js = ijkgs(2)
                        ks = ijkgs(3)

                        ijkgt(t_nd) = ijkt(t_nd) + n*t_lr
                        it = ijkgt(1)
                        jt = ijkgt(2)
                        kt = ijkgt(3)

                        do m=nst,ned
                            td = ord(m)
                            mb_der(nbt)%fld(m)%r3d(it,jt,kt) = sgn(m)*inters(nint)%dat(is,js,ks,td)
                        end do
                    end do
                end do
                end do
                end do

                if (nswmem == nbc_inter_buf_dyn) then
                    deallocate(inters(nint)%dat,stat=ierr)
                else
                    nullify(inters(nint)%dat)
                end if
#ifdef PARALLEL
            end if
        end if
#endif
    end do

#ifdef PARALLEL
    if (nmsg_int > 0) then
        call MPI_WAITALL(nmsg_int, request_int(1:nmsg_int), status_int, ierr)

        do nint=1,ninters
            nbs     = inters(nint)%bc%nbs
            nbt     = inters(nint)%bc%nbt
            id_src = mb_top(nbs)%pid
            id_des = mb_top(nbt)%pid

            s_nd    = inters(nint)%bc%s_nd
            s_lr    = inters(nint)%bc%s_lr

            nrt     = inters(nint)%bc%nrt
            t_st(:) = inters(nint)%bc%t_st(:)
            t_ed(:) = inters(nint)%bc%t_ed(:)
            t_nd    = inters(nint)%bc%t_nd
            t_lr    = inters(nint)%bc%t_lr

            if (id_src /= id_des) then

                if (myid == id_des) then
                    do k=t_st(3),t_ed(3)
                    do j=t_st(2),t_ed(2)
                    do i=t_st(1),t_ed(1)
                        ijkt(:) = (/i,j,k/)
                        ijks(:) = mb_top(nbt)%bcs(nrt)%mapijk(i,j,k,:)
                        ijkgt(:) = ijkt(:)
                        ijkgs(:) = ijks(:)
                        do n=1,ngh
                            ijkgs(s_nd) = ijks(s_nd) - n*s_lr
                            is = ijkgs(1)
                            js = ijkgs(2)
                            ks = ijkgs(3)

                            ijkgt(t_nd) = ijkt(t_nd) + n*t_lr
                            it = ijkgt(1)
                            jt = ijkgt(2)
                            kt = ijkgt(3)

                            do m=nst,ned
                                td = ord(m)
                                mb_der(nbt)%fld(m)%r3d(it,jt,kt) = sgn(m)*inters(nint)%dat(is,js,ks,td)
                            end do
                        end do
                    end do
                    end do
                    end do

                end if

                if (myid == id_src .or. myid == id_des) then
                    if (nswmem == nbc_inter_buf_dyn) then
                        deallocate(inters(nint)%dat,stat=ierr)
                    else
                        nullify(inters(nint)%dat)
                    end if
                end if

            end if
        end do

    end if
#endif
end subroutine exchange_bc_der

subroutine exchange_cc_bc_der(mb_der,nst,ned,ngh,nswmem)
    use mod_kndconsts, only : kind_int
    use mod_constants, only : nbc_inter_buf_dyn
    use mod_datatypes, only : var_block_t
    use mod_fieldvars, only : nblocks,mb_topc,ninters,interscc
    use mod_parallels
    implicit none
    type(var_block_t), pointer, intent(in) :: mb_der(:)
    integer(kind_int),          intent(in) :: nst,ned
    integer(kind_int),          intent(in) :: ngh
    integer(kind_int),          intent(in) :: nswmem
    integer(kind_int) :: nint,ierr
    integer(kind_int) :: nbs,nrs,s_st(3),s_ed(3),s_nd,s_lr
    integer(kind_int) :: nbt,nrt,t_st(3),t_ed(3),t_nd,t_lr
    integer(kind_int) :: s_ord(3),s_sgn(3),t_ord(3),t_sgn(3)
    integer(kind_int) :: st(3),ed(3),sgn(nst:ned),ord(nst:ned)
    integer(kind_int) :: ijks(3),ijkt(3),ijkgs(3),ijkgt(3)
    integer(kind_int) :: i,j,k,m,n,is,js,ks,it,jt,kt,sd,td
#ifdef PARALLEL
    integer(kind_int) :: id_src,id_des,tag,packsize
    integer(kind_int) :: index_msg
#endif

    ierr = mod(ned-nst+1,3)
    call error_check(ierr, &
                     "The size of array isn't times of 3 in subroutine exchange_bc_der")

#ifdef PARALLEL
    index_msg = 0
#endif

    do nint=1,ninters
        nbs     = interscc(nint)%bc%nbs
        nrs     = interscc(nint)%bc%nrs
        s_st(:) = interscc(nint)%bc%s_st(:)
        s_ed(:) = interscc(nint)%bc%s_ed(:)
        s_nd    = interscc(nint)%bc%s_nd
        s_lr    = interscc(nint)%bc%s_lr

        nbt     = interscc(nint)%bc%nbt
        nrt     = interscc(nint)%bc%nrt
        t_st(:) = interscc(nint)%bc%t_st(:)
        t_ed(:) = interscc(nint)%bc%t_ed(:)
        t_nd    = interscc(nint)%bc%t_nd
        t_lr    = interscc(nint)%bc%t_lr

        s_ord(:) = mb_topc(nbt)%bcs(nrt)%s_ord(:)
        s_sgn(:) = mb_topc(nbt)%bcs(nrt)%s_sgn(:)
        t_ord(:) = mb_topc(nbt)%bcs(nrt)%t_ord(:)
        t_sgn(:) = mb_topc(nbt)%bcs(nrt)%t_sgn(:)

        do m=nst,ned
           sd = m-nst+1 - (m-nst)/3*3
           ord(m) = t_ord(sd) + (m-nst)/3*3
           sgn(m) = t_sgn(sd)
        end do

        call bc_extend_inward(s_st,s_ed,s_nd,s_lr,ngh,st,ed)

#ifdef PARALLEL
        packsize = product(ed(:)-st(:)+1)*(ned-nst+1)

        id_src = mb_topc(nbs)%pid
        id_des = mb_topc(nbt)%pid
        tag    = nint + ninters*3


        if (id_src /= id_des) then
            if (myid == id_src) then
                if (nswmem == nbc_inter_buf_dyn) then
                    allocate(interscc(nint)%dat(st(1):ed(1), &
                                                st(2):ed(2), &
                                                st(3):ed(3), &
                                                nst:ned),stat=ierr)
                else
                    interscc(nint)%dat => interscc(nint)%buf(nswmem)%dat
                end if

                do m=nst,ned
                    do k=st(3),ed(3)
                    do j=st(2),ed(2)
                    do i=st(1),ed(1)
                        interscc(nint)%dat(i,j,k,m) = mb_der(nbs)%fld(m)%r3d(i,j,k)
                    end do
                    end do
                    end do
                end do

                index_msg = index_msg + 1
                call MPI_ISEND(interscc(nint)%dat,packsize,kind_real_mpi, &
                              id_des,tag,MPI_COMM_WORLD,request_int(index_msg),ierr)

            end if

            if (myid == id_des) then
                if (nswmem == nbc_inter_buf_dyn) then
                    allocate(interscc(nint)%dat(st(1):ed(1), &
                                                st(2):ed(2), &
                                                st(3):ed(3), &
                                                nst:ned),stat=ierr)
                else
                    interscc(nint)%dat => interscc(nint)%buf(nswmem)%dat
                end if

                index_msg = index_msg + 1
                call MPI_IRECV(interscc(nint)%dat,packsize,kind_real_mpi, &
                              id_src,tag,MPI_COMM_WORLD,request_int(index_msg),ierr)

            end if

        else
            if (myid == id_src) then
#endif
                if (nswmem == nbc_inter_buf_dyn) then
                    allocate(interscc(nint)%dat(st(1):ed(1), &
                                                st(2):ed(2), &
                                                st(3):ed(3), &
                                                nst:ned),stat=ierr)
                else
                    interscc(nint)%dat => interscc(nint)%buf(nswmem)%dat
                end if

                do m=nst,ned
                    do k=st(3),ed(3)
                    do j=st(2),ed(2)
                    do i=st(1),ed(1)
                        interscc(nint)%dat(i,j,k,m) = mb_der(nbs)%fld(m)%r3d(i,j,k)
                    end do
                    end do
                    end do
                end do

                do k=t_st(3),t_ed(3)
                do j=t_st(2),t_ed(2)
                do i=t_st(1),t_ed(1)
                    ijkt(:) = (/i,j,k/)
                    ijks(:) = mb_topc(nbt)%bcs(nrt)%mapijk(i,j,k,:)
                    ijkgt(:) = ijkt(:)
                    ijkgs(:) = ijks(:)
                    do n=1,ngh
                        ijkgs(s_nd) = ijks(s_nd) - n*s_lr
                        is = ijkgs(1)
                        js = ijkgs(2)
                        ks = ijkgs(3)

                        ijkgt(t_nd) = ijkt(t_nd) + n*t_lr
                        it = ijkgt(1)
                        jt = ijkgt(2)
                        kt = ijkgt(3)

                        do m=nst,ned
                            td = ord(m)
                            mb_der(nbt)%fld(m)%r3d(it,jt,kt) = sgn(m)*interscc(nint)%dat(is,js,ks,td)
                        end do
                    end do
                end do
                end do
                end do

                if (nswmem == nbc_inter_buf_dyn) then
                    deallocate(interscc(nint)%dat,stat=ierr)
                else
                    nullify(interscc(nint)%dat)
                end if
#ifdef PARALLEL
            end if
        end if
#endif
    end do

#ifdef PARALLEL
    if (nmsg_int > 0) then
        call MPI_WAITALL(nmsg_int, request_int(1:nmsg_int), status_int, ierr)

        do nint=1,ninters
            nbs     = interscc(nint)%bc%nbs
            nbt     = interscc(nint)%bc%nbt
            id_src = mb_topc(nbs)%pid
            id_des = mb_topc(nbt)%pid

            s_nd    = interscc(nint)%bc%s_nd
            s_lr    = interscc(nint)%bc%s_lr

            nrt     = interscc(nint)%bc%nrt
            t_st(:) = interscc(nint)%bc%t_st(:)
            t_ed(:) = interscc(nint)%bc%t_ed(:)
            t_nd    = interscc(nint)%bc%t_nd
            t_lr    = interscc(nint)%bc%t_lr

            if (id_src /= id_des) then

                if (myid == id_des) then
                    do k=t_st(3),t_ed(3)
                    do j=t_st(2),t_ed(2)
                    do i=t_st(1),t_ed(1)
                        ijkt(:) = (/i,j,k/)
                        ijks(:) = mb_topc(nbt)%bcs(nrt)%mapijk(i,j,k,:)
                        ijkgt(:) = ijkt(:)
                        ijkgs(:) = ijks(:)
                        do n=1,ngh
                            ijkgs(s_nd) = ijks(s_nd) - n*s_lr
                            is = ijkgs(1)
                            js = ijkgs(2)
                            ks = ijkgs(3)

                            ijkgt(t_nd) = ijkt(t_nd) + n*t_lr
                            it = ijkgt(1)
                            jt = ijkgt(2)
                            kt = ijkgt(3)

                            do m=nst,ned
                                td = ord(m)
                                mb_der(nbt)%fld(m)%r3d(it,jt,kt) = sgn(m)*interscc(nint)%dat(is,js,ks,td)
                            end do
                        end do
                    end do
                    end do
                    end do

                end if

                if (myid == id_src .or. myid == id_des) then
                    if (nswmem == nbc_inter_buf_dyn) then
                        deallocate(interscc(nint)%dat,stat=ierr)
                    else
                        nullify(interscc(nint)%dat)
                    end if
                end if

            end if
        end do

    end if
#endif
end subroutine exchange_cc_bc_der

subroutine exchange_sp_bc_der(mb_der,nst,ned,ngh,nswmem)
    use mod_kndconsts, only : kind_int
    use mod_constants, only : nbc_inter_buf_dyn
    use mod_datatypes, only : var_block_t
    use mod_fieldvars, only : nblocks,mb_topsp,ninters,interssp
    use mod_parallels
    implicit none
    type(var_block_t), pointer, intent(in) :: mb_der(:)
    integer(kind_int),          intent(in) :: nst,ned
    integer(kind_int),          intent(in) :: ngh
    integer(kind_int),          intent(in) :: nswmem
    integer(kind_int) :: nint,ierr
    integer(kind_int) :: nbs,nrs,s_st(3),s_ed(3),s_nd,s_lr
    integer(kind_int) :: nbt,nrt,t_st(3),t_ed(3),t_nd,t_lr
    integer(kind_int) :: s_ord(3),s_sgn(3),t_ord(3),t_sgn(3)
    integer(kind_int) :: st(3),ed(3),sgn(nst:ned),ord(nst:ned)
    integer(kind_int) :: ijks(3),ijkt(3),ijkgs(3),ijkgt(3)
    integer(kind_int) :: i,j,k,m,n,is,js,ks,it,jt,kt,sd,td
#ifdef PARALLEL
    integer(kind_int) :: id_src,id_des,tag,packsize
    integer(kind_int) :: index_msg
#endif

    ierr = mod(ned-nst+1,3)
    call error_check(ierr, &
                     "The size of array isn't times of 3 in subroutine exchange_bc_der")

#ifdef PARALLEL
    index_msg = 0
#endif

    do nint=1,ninters
        nbs     = interssp(nint)%bc%nbs
        nrs     = interssp(nint)%bc%nrs
        s_st(:) = interssp(nint)%bc%s_st(:)
        s_ed(:) = interssp(nint)%bc%s_ed(:)
        s_nd    = interssp(nint)%bc%s_nd
        s_lr    = interssp(nint)%bc%s_lr

        nbt     = interssp(nint)%bc%nbt
        nrt     = interssp(nint)%bc%nrt
        t_st(:) = interssp(nint)%bc%t_st(:)
        t_ed(:) = interssp(nint)%bc%t_ed(:)
        t_nd    = interssp(nint)%bc%t_nd
        t_lr    = interssp(nint)%bc%t_lr

        s_ord(:) = mb_topsp(nbt)%bcs(nrt)%s_ord(:)
        s_sgn(:) = mb_topsp(nbt)%bcs(nrt)%s_sgn(:)
        t_ord(:) = mb_topsp(nbt)%bcs(nrt)%t_ord(:)
        t_sgn(:) = mb_topsp(nbt)%bcs(nrt)%t_sgn(:)

        do m=nst,ned
           sd = m-nst+1 - (m-nst)/3*3
           ord(m) = t_ord(sd) + (m-nst)/3*3
           sgn(m) = t_sgn(sd)
        end do

        call bc_extend_inward(s_st,s_ed,s_nd,s_lr,ngh,st,ed)

#ifdef PARALLEL
        packsize = product(ed(:)-st(:)+1)*(ned-nst+1)

        id_src = mb_topsp(nbs)%pid
        id_des = mb_topsp(nbt)%pid
        tag    = nint + ninters*3


        if (id_src /= id_des) then
            if (myid == id_src) then
                if (nswmem == nbc_inter_buf_dyn) then
                    allocate(interssp(nint)%dat(st(1):ed(1), &
                                                st(2):ed(2), &
                                                st(3):ed(3), &
                                                nst:ned),stat=ierr)
                else
                    interssp(nint)%dat => interssp(nint)%buf(nswmem)%dat
                end if

                do m=nst,ned
                    do k=st(3),ed(3)
                    do j=st(2),ed(2)
                    do i=st(1),ed(1)
                        interssp(nint)%dat(i,j,k,m) = mb_der(nbs)%fld(m)%r3d(i,j,k)
                    end do
                    end do
                    end do
                end do

                index_msg = index_msg + 1
                call MPI_ISEND(interssp(nint)%dat,packsize,kind_real_mpi, &
                              id_des,tag,MPI_COMM_WORLD,request_int(index_msg),ierr)

            end if

            if (myid == id_des) then
                if (nswmem == nbc_inter_buf_dyn) then
                    allocate(interssp(nint)%dat(st(1):ed(1), &
                                              st(2):ed(2), &
                                              st(3):ed(3), &
                                              nst:ned),stat=ierr)
                else
                    interssp(nint)%dat => interssp(nint)%buf(nswmem)%dat
                end if

                index_msg = index_msg + 1
                call MPI_IRECV(interssp(nint)%dat,packsize,kind_real_mpi, &
                              id_src,tag,MPI_COMM_WORLD,request_int(index_msg),ierr)

            end if

        else
            if (myid == id_src) then
#endif
                if (nswmem == nbc_inter_buf_dyn) then
                    allocate(interssp(nint)%dat(st(1):ed(1), &
                                                st(2):ed(2), &
                                                st(3):ed(3), &
                                                nst:ned),stat=ierr)
                else
                    interssp(nint)%dat => interssp(nint)%buf(nswmem)%dat
                end if

                do m=nst,ned
                    do k=st(3),ed(3)
                    do j=st(2),ed(2)
                    do i=st(1),ed(1)
                        interssp(nint)%dat(i,j,k,m) = mb_der(nbs)%fld(m)%r3d(i,j,k)
                    end do
                    end do
                    end do
                end do

                do k=t_st(3),t_ed(3)
                do j=t_st(2),t_ed(2)
                do i=t_st(1),t_ed(1)
                    ijkt(:) = (/i,j,k/)
                    ijks(:) = mb_topsp(nbt)%bcs(nrt)%mapijk(i,j,k,:)
                    ijkgt(:) = ijkt(:)
                    ijkgs(:) = ijks(:)
                    do n=1,ngh
                        ijkgs(s_nd) = ijks(s_nd) - (n-1)*s_lr
                        is = ijkgs(1)
                        js = ijkgs(2)
                        ks = ijkgs(3)

                        ijkgt(t_nd) = ijkt(t_nd) + n*t_lr
                        it = ijkgt(1)
                        jt = ijkgt(2)
                        kt = ijkgt(3)

                        do m=nst,ned
                            td = ord(m)
                            mb_der(nbt)%fld(m)%r3d(it,jt,kt) = sgn(m)*interssp(nint)%dat(is,js,ks,td)
                        end do
                    end do
                end do
                end do
                end do

                if (nswmem == nbc_inter_buf_dyn) then
                    deallocate(interssp(nint)%dat,stat=ierr)
                else
                    nullify(interssp(nint)%dat)
                end if
#ifdef PARALLEL
            end if
        end if
#endif
    end do

#ifdef PARALLEL
    if (nmsg_int > 0) then
        call MPI_WAITALL(nmsg_int, request_int(1:nmsg_int), status_int, ierr)

        do nint=1,ninters
            nbs     = interssp(nint)%bc%nbs
            nbt     = interssp(nint)%bc%nbt
            id_src = mb_topsp(nbs)%pid
            id_des = mb_topsp(nbt)%pid

            s_nd    = interssp(nint)%bc%s_nd
            s_lr    = interssp(nint)%bc%s_lr

            nrt     = interssp(nint)%bc%nrt
            t_st(:) = interssp(nint)%bc%t_st(:)
            t_ed(:) = interssp(nint)%bc%t_ed(:)
            t_nd    = interssp(nint)%bc%t_nd
            t_lr    = interssp(nint)%bc%t_lr

            if (id_src /= id_des) then

                if (myid == id_des) then
                    do k=t_st(3),t_ed(3)
                    do j=t_st(2),t_ed(2)
                    do i=t_st(1),t_ed(1)
                        ijkt(:) = (/i,j,k/)
                        ijks(:) = mb_topsp(nbt)%bcs(nrt)%mapijk(i,j,k,:)
                        ijkgt(:) = ijkt(:)
                        ijkgs(:) = ijks(:)
                        do n=1,ngh
                            ijkgs(s_nd) = ijks(s_nd) - (n-1)*s_lr
                            is = ijkgs(1)
                            js = ijkgs(2)
                            ks = ijkgs(3)

                            ijkgt(t_nd) = ijkt(t_nd) + n*t_lr
                            it = ijkgt(1)
                            jt = ijkgt(2)
                            kt = ijkgt(3)

                            do m=nst,ned
                                td = ord(m)
                                mb_der(nbt)%fld(m)%r3d(it,jt,kt) = sgn(m)*interssp(nint)%dat(is,js,ks,td)
                            end do
                        end do
                    end do
                    end do
                    end do

                end if

                if (myid == id_src .or. myid == id_des) then
                    if (nswmem == nbc_inter_buf_dyn) then
                        deallocate(interssp(nint)%dat,stat=ierr)
                    else
                        nullify(interssp(nint)%dat)
                    end if
                end if

            end if
        end do

    end if
#endif
end subroutine exchange_sp_bc_der

!>
!! @brief 交换边界网格导数函数
!! @details 以对接边界为基本单位，逐一交换网格导数
!! @param[in] mb_sxyz 网格导数变量指针
!! @param[in] nst     网格导数变量起始索引
!! @param[in] ned     网格导数变量终止索引
!! @param[in] ngh     边界网格导数层数
!! @author 刘化勇
!! @date 2011年05月31日
!! @author 王勇献
!! @remark 2012年02月14日 修改为非阻塞式通信
!>
subroutine exchange_bc_sxyz(mb_sxyz,nst,ned,ngh)
    use mod_kndconsts, only : kind_int
    use mod_datatypes, only : var_block_t
    use mod_fieldvars, only : nblocks,mb_top,ninters,inters
    use mod_parallels
    implicit none
    type(var_block_t), pointer, intent(in) :: mb_sxyz(:)
    integer(kind_int),          intent(in) :: nst,ned,ngh
    integer(kind_int) :: nint,ierr
    integer(kind_int) :: nbs,nrs,s_st(3),s_ed(3),s_nd,s_lr
    integer(kind_int) :: nbt,nrt,t_st(3),t_ed(3),t_nd,t_lr
    integer(kind_int) :: s_ord(3),s_sgn(3),t_ord(3),t_sgn(3)
    integer(kind_int) :: st(3),ed(3),sgn(nst:ned),ord(nst:ned)
    integer(kind_int) :: ijks(3),ijkt(3),ijkgs(3),ijkgt(3)
    integer(kind_int) :: i,j,k,m,n,is,js,ks,it,jt,kt,sd,td
#ifdef PARALLEL
    integer(kind_int) :: id_src,id_des,tag,packsize
    integer(kind_int) :: index_msg
#endif

    ierr = ned - nst - 8
    call error_check(ierr, &
                     "The size of array isn't 9 in subroutine exchange_bc_sxyz")

#ifdef PARALLEL
    index_msg = 0
#endif

    do nint=1,ninters
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

        s_ord(:) = mb_top(nbt)%bcs(nrt)%s_ord(:)
        s_sgn(:) = mb_top(nbt)%bcs(nrt)%s_sgn(:)
        t_ord(:) = mb_top(nbt)%bcs(nrt)%t_ord(:)
        t_sgn(:) = mb_top(nbt)%bcs(nrt)%t_sgn(:)

        do m=nst,ned
           sd = (m-nst)/3 + 1
           ord(m) = (t_ord(sd) - 1)*3 + (m-nst+1 - (m-nst)/3*3)
           sgn(m) = t_sgn(sd)
        end do

        call bc_extend_inward(s_st,s_ed,s_nd,s_lr,ngh,st,ed)


#ifdef PARALLEL
        packsize = product(ed(:)-st(:)+1)*(ned-nst+1)

        id_src = mb_top(nbs)%pid
        id_des = mb_top(nbt)%pid
        tag    = nint + ninters*3


        if (id_src /= id_des) then
            if (myid == id_src) then
                allocate(inters(nint)%dat(st(1):ed(1), &
                                          st(2):ed(2), &
                                          st(3):ed(3), &
                                          nst:ned),stat=ierr)

                do m=nst,ned
                    do k=st(3),ed(3)
                    do j=st(2),ed(2)
                    do i=st(1),ed(1)
                        inters(nint)%dat(i,j,k,m) = mb_sxyz(nbs)%fld(m)%r3d(i,j,k)
                    end do
                    end do
                    end do
                end do

                index_msg = index_msg + 1
                call MPI_ISEND(inters(nint)%dat,packsize,kind_real_mpi, &
                              id_des,tag,MPI_COMM_WORLD,request_int(index_msg),ierr)
            end if

            if (myid == id_des) then
                allocate(inters(nint)%dat(st(1):ed(1), &
                                          st(2):ed(2), &
                                          st(3):ed(3), &
                                          nst:ned),stat=ierr)

                index_msg = index_msg + 1
                call MPI_IRECV(inters(nint)%dat,packsize,kind_real_mpi, &
                              id_src,tag,MPI_COMM_WORLD,request_int(index_msg),ierr)
            end if

        else
            if (myid == id_src) then
#endif
                allocate(inters(nint)%dat(st(1):ed(1), &
                                          st(2):ed(2), &
                                          st(3):ed(3), &
                                          nst:ned),stat=ierr)

                do m=nst,ned
                    do k=st(3),ed(3)
                    do j=st(2),ed(2)
                    do i=st(1),ed(1)
                        inters(nint)%dat(i,j,k,m) = mb_sxyz(nbs)%fld(m)%r3d(i,j,k)
                    end do
                    end do
                    end do
                end do


                do k=t_st(3),t_ed(3)
                do j=t_st(2),t_ed(2)
                do i=t_st(1),t_ed(1)
                    ijkt(:) = (/i,j,k/)
                    ijks(:) = mb_top(nbt)%bcs(nrt)%mapijk(i,j,k,:)
                    ijkgt(:) = ijkt(:)
                    ijkgs(:) = ijks(:)
                    do n=1,ngh
                        ijkgs(s_nd) = ijks(s_nd) - n*s_lr
                        is = ijkgs(1)
                        js = ijkgs(2)
                        ks = ijkgs(3)

                        ijkgt(t_nd) = ijkt(t_nd) + n*t_lr
                        it = ijkgt(1)
                        jt = ijkgt(2)
                        kt = ijkgt(3)

                        do m=nst,ned
                            td = ord(m)
                            mb_sxyz(nbt)%fld(m)%r3d(it,jt,kt) = sgn(m)*inters(nint)%dat(is,js,ks,td)
                        end do
                    end do
                end do
                end do
                end do

                deallocate(inters(nint)%dat,stat=ierr)
#ifdef PARALLEL
            end if
        end if
#endif
    end do

#ifdef PARALLEL
    if (nmsg_int > 0) then
        call MPI_WAITALL(nmsg_int, request_int(1:nmsg_int), status_int, ierr)

        do nint=1,ninters
            nbs     = inters(nint)%bc%nbs
            nbt     = inters(nint)%bc%nbt
            id_src = mb_top(nbs)%pid
            id_des = mb_top(nbt)%pid

            if (id_src /= id_des) then
                if (myid == id_src) then
                    deallocate(inters(nint)%dat,stat=ierr)
                end if

                if (myid == id_des) then
                    s_nd    = inters(nint)%bc%s_nd
                    s_lr    = inters(nint)%bc%s_lr

                    nrt     = inters(nint)%bc%nrt
                    t_st(:) = inters(nint)%bc%t_st(:)
                    t_ed(:) = inters(nint)%bc%t_ed(:)
                    t_nd    = inters(nint)%bc%t_nd
                    t_lr    = inters(nint)%bc%t_lr

                    t_ord(:) = mb_top(nbt)%bcs(nrt)%t_ord(:)
                    t_sgn(:) = mb_top(nbt)%bcs(nrt)%t_sgn(:)

                    do m=nst,ned
                       sd = (m-nst)/3 + 1
                       ord(m) = (t_ord(sd) - 1)*3 + (m-nst+1 - (m-nst)/3*3)
                       sgn(m) = t_sgn(sd)
                    end do

                    do k=t_st(3),t_ed(3)
                    do j=t_st(2),t_ed(2)
                    do i=t_st(1),t_ed(1)
                        ijkt(:) = (/i,j,k/)
                        ijks(:) = mb_top(nbt)%bcs(nrt)%mapijk(i,j,k,:)
                        ijkgt(:) = ijkt(:)
                        ijkgs(:) = ijks(:)
                        do n=1,ngh
                            ijkgs(s_nd) = ijks(s_nd) - n*s_lr
                            is = ijkgs(1)
                            js = ijkgs(2)
                            ks = ijkgs(3)

                            ijkgt(t_nd) = ijkt(t_nd) + n*t_lr
                            it = ijkgt(1)
                            jt = ijkgt(2)
                            kt = ijkgt(3)

                            do m=nst,ned
                                td = ord(m)
                                mb_sxyz(nbt)%fld(m)%r3d(it,jt,kt) = sgn(m)*inters(nint)%dat(is,js,ks,td)
                            end do
                        end do
                    end do
                    end do
                    end do

                    deallocate(inters(nint)%dat,stat=ierr)
                end if
            end if

        end do

    end if
#endif
end subroutine exchange_bc_sxyz

subroutine exchange_cc_bc_sxyz(mb_sxyzcc,nst,ned,ngh)
    use mod_kndconsts, only : kind_int
    use mod_datatypes, only : var_block_t
    use mod_fieldvars, only : nblocks,mb_topc,ninters,interscc
    use mod_parallels
    implicit none
    type(var_block_t), pointer, intent(in) :: mb_sxyzcc(:)
    integer(kind_int),          intent(in) :: nst,ned,ngh
    integer(kind_int) :: nint,ierr
    integer(kind_int) :: nbs,nrs,s_st(3),s_ed(3),s_nd,s_lr
    integer(kind_int) :: nbt,nrt,t_st(3),t_ed(3),t_nd,t_lr
    integer(kind_int) :: s_ord(3),s_sgn(3),t_ord(3),t_sgn(3)
    integer(kind_int) :: st(3),ed(3),sgn(nst:ned),ord(nst:ned)
    integer(kind_int) :: ijks(3),ijkt(3),ijkgs(3),ijkgt(3)
    integer(kind_int) :: i,j,k,m,n,is,js,ks,it,jt,kt,sd,td
#ifdef PARALLEL
    integer(kind_int) :: id_src,id_des,tag,packsize
    integer(kind_int) :: index_msg
#endif

    ierr = ned - nst - 8
    call error_check(ierr, &
                     "The size of array isn't 9 in subroutine exchange_bc_sxyz")

#ifdef PARALLEL
    index_msg = 0
#endif

    do nint=1,ninters
        nbs     = interscc(nint)%bc%nbs
        nrs     = interscc(nint)%bc%nrs
        s_st(:) = interscc(nint)%bc%s_st(:)
        s_ed(:) = interscc(nint)%bc%s_ed(:)
        s_nd    = interscc(nint)%bc%s_nd
        s_lr    = interscc(nint)%bc%s_lr

        nbt     = interscc(nint)%bc%nbt
        nrt     = interscc(nint)%bc%nrt
        t_st(:) = interscc(nint)%bc%t_st(:)
        t_ed(:) = interscc(nint)%bc%t_ed(:)
        t_nd    = interscc(nint)%bc%t_nd
        t_lr    = interscc(nint)%bc%t_lr

        s_ord(:) = mb_topc(nbt)%bcs(nrt)%s_ord(:)
        s_sgn(:) = mb_topc(nbt)%bcs(nrt)%s_sgn(:)
        t_ord(:) = mb_topc(nbt)%bcs(nrt)%t_ord(:)
        t_sgn(:) = mb_topc(nbt)%bcs(nrt)%t_sgn(:)

        do m=nst,ned
           sd = (m-nst)/3 + 1
           ord(m) = (t_ord(sd) - 1)*3 + (m-nst+1 - (m-nst)/3*3)
           sgn(m) = t_sgn(sd)
        end do

        call bc_extend_inward(s_st,s_ed,s_nd,s_lr,ngh,st,ed)


#ifdef PARALLEL
        packsize = product(ed(:)-st(:)+1)*(ned-nst+1)

        id_src = mb_topc(nbs)%pid
        id_des = mb_topc(nbt)%pid
        tag    = nint + ninters*3


        if (id_src /= id_des) then
            if (myid == id_src) then
                allocate(interscc(nint)%dat(st(1):ed(1), &
                                            st(2):ed(2), &
                                            st(3):ed(3), &
                                            nst:ned),stat=ierr)

                do m=nst,ned
                    do k=st(3),ed(3)
                    do j=st(2),ed(2)
                    do i=st(1),ed(1)
                        interscc(nint)%dat(i,j,k,m) = mb_sxyzcc(nbs)%fld(m)%r3d(i,j,k)
                    end do
                    end do
                    end do
                end do

                index_msg = index_msg + 1
                call MPI_ISEND(interscc(nint)%dat,packsize,kind_real_mpi, &
                              id_des,tag,MPI_COMM_WORLD,request_int(index_msg),ierr)
            end if

            if (myid == id_des) then
                allocate(interscc(nint)%dat(st(1):ed(1), &
                                            st(2):ed(2), &
                                            st(3):ed(3), &
                                            nst:ned),stat=ierr)

                index_msg = index_msg + 1
                call MPI_IRECV(interscc(nint)%dat,packsize,kind_real_mpi, &
                              id_src,tag,MPI_COMM_WORLD,request_int(index_msg),ierr)
            end if

        else
            if (myid == id_src) then
#endif
                allocate(interscc(nint)%dat(st(1):ed(1), &
                                            st(2):ed(2), &
                                            st(3):ed(3), &
                                            nst:ned),stat=ierr)

                do m=nst,ned
                    do k=st(3),ed(3)
                    do j=st(2),ed(2)
                    do i=st(1),ed(1)
                        interscc(nint)%dat(i,j,k,m) = mb_sxyzcc(nbs)%fld(m)%r3d(i,j,k)
                    end do
                    end do
                    end do
                end do


                do k=t_st(3),t_ed(3)
                do j=t_st(2),t_ed(2)
                do i=t_st(1),t_ed(1)
                    ijkt(:) = (/i,j,k/)
                    ijks(:) = mb_topc(nbt)%bcs(nrt)%mapijk(i,j,k,:)
                    ijkgt(:) = ijkt(:)
                    ijkgs(:) = ijks(:)
                    do n=1,ngh
                        ijkgs(s_nd) = ijks(s_nd) - n*s_lr
                        is = ijkgs(1)
                        js = ijkgs(2)
                        ks = ijkgs(3)

                        ijkgt(t_nd) = ijkt(t_nd) + n*t_lr
                        it = ijkgt(1)
                        jt = ijkgt(2)
                        kt = ijkgt(3)

                        do m=nst,ned
                            td = ord(m)
                            mb_sxyzcc(nbt)%fld(m)%r3d(it,jt,kt) = sgn(m)*interscc(nint)%dat(is,js,ks,td)
                        end do
                    end do
                end do
                end do
                end do

                deallocate(interscc(nint)%dat,stat=ierr)
#ifdef PARALLEL
            end if
        end if
#endif
    end do

#ifdef PARALLEL
    if (nmsg_int > 0) then
        call MPI_WAITALL(nmsg_int, request_int(1:nmsg_int), status_int, ierr)

        do nint=1,ninters
            nbs     = interscc(nint)%bc%nbs
            nbt     = interscc(nint)%bc%nbt
            id_src = mb_topc(nbs)%pid
            id_des = mb_topc(nbt)%pid

            if (id_src /= id_des) then
                if (myid == id_src) then
                    deallocate(interscc(nint)%dat,stat=ierr)
                end if

                if (myid == id_des) then
                    s_nd    = interscc(nint)%bc%s_nd
                    s_lr    = interscc(nint)%bc%s_lr

                    nrt     = interscc(nint)%bc%nrt
                    t_st(:) = interscc(nint)%bc%t_st(:)
                    t_ed(:) = interscc(nint)%bc%t_ed(:)
                    t_nd    = interscc(nint)%bc%t_nd
                    t_lr    = interscc(nint)%bc%t_lr

                    t_ord(:) = mb_topc(nbt)%bcs(nrt)%t_ord(:)
                    t_sgn(:) = mb_topc(nbt)%bcs(nrt)%t_sgn(:)

                    do m=nst,ned
                       sd = (m-nst)/3 + 1
                       ord(m) = (t_ord(sd) - 1)*3 + (m-nst+1 - (m-nst)/3*3)
                       sgn(m) = t_sgn(sd)
                    end do

                    do k=t_st(3),t_ed(3)
                    do j=t_st(2),t_ed(2)
                    do i=t_st(1),t_ed(1)
                        ijkt(:) = (/i,j,k/)
                        ijks(:) = mb_topc(nbt)%bcs(nrt)%mapijk(i,j,k,:)
                        ijkgt(:) = ijkt(:)
                        ijkgs(:) = ijks(:)
                        do n=1,ngh
                            ijkgs(s_nd) = ijks(s_nd) - n*s_lr
                            is = ijkgs(1)
                            js = ijkgs(2)
                            ks = ijkgs(3)

                            ijkgt(t_nd) = ijkt(t_nd) + n*t_lr
                            it = ijkgt(1)
                            jt = ijkgt(2)
                            kt = ijkgt(3)

                            do m=nst,ned
                                td = ord(m)
                                mb_sxyzcc(nbt)%fld(m)%r3d(it,jt,kt) = sgn(m)*interscc(nint)%dat(is,js,ks,td)
                            end do
                        end do
                    end do
                    end do
                    end do

                    deallocate(interscc(nint)%dat,stat=ierr)
                end if
            end if

        end do

    end if
#endif
end subroutine exchange_cc_bc_sxyz

subroutine exchange_sp_bc_sxyz(mb_sxyzsp,nst,ned,ngh)
    use mod_kndconsts, only : kind_int
    use mod_datatypes, only : var_block_t
    use mod_fieldvars, only : nblocks,mb_topsp,ninters,interssp
    use mod_parallels
    implicit none
    type(var_block_t), pointer, intent(in) :: mb_sxyzsp(:)
    integer(kind_int),          intent(in) :: nst,ned,ngh
    integer(kind_int) :: nint,ierr
    integer(kind_int) :: nbs,nrs,s_st(3),s_ed(3),s_nd,s_lr
    integer(kind_int) :: nbt,nrt,t_st(3),t_ed(3),t_nd,t_lr
    integer(kind_int) :: s_ord(3),s_sgn(3),t_ord(3),t_sgn(3)
    integer(kind_int) :: st(3),ed(3),sgn(nst:ned),ord(nst:ned)
    integer(kind_int) :: ijks(3),ijkt(3),ijkgs(3),ijkgt(3)
    integer(kind_int) :: i,j,k,m,n,is,js,ks,it,jt,kt,sd,td
#ifdef PARALLEL
    integer(kind_int) :: id_src,id_des,tag,packsize
    integer(kind_int) :: index_msg
#endif

    ierr = ned - nst - 8
    call error_check(ierr, &
                     "The size of array isn't 9 in subroutine exchange_bc_sxyz")

#ifdef PARALLEL
    index_msg = 0
#endif

    do nint=1,ninters
        nbs     = interssp(nint)%bc%nbs
        nrs     = interssp(nint)%bc%nrs
        s_st(:) = interssp(nint)%bc%s_st(:)
        s_ed(:) = interssp(nint)%bc%s_ed(:)
        s_nd    = interssp(nint)%bc%s_nd
        s_lr    = interssp(nint)%bc%s_lr

        nbt     = interssp(nint)%bc%nbt
        nrt     = interssp(nint)%bc%nrt
        t_st(:) = interssp(nint)%bc%t_st(:)
        t_ed(:) = interssp(nint)%bc%t_ed(:)
        t_nd    = interssp(nint)%bc%t_nd
        t_lr    = interssp(nint)%bc%t_lr

        s_ord(:) = mb_topsp(nbt)%bcs(nrt)%s_ord(:)
        s_sgn(:) = mb_topsp(nbt)%bcs(nrt)%s_sgn(:)
        t_ord(:) = mb_topsp(nbt)%bcs(nrt)%t_ord(:)
        t_sgn(:) = mb_topsp(nbt)%bcs(nrt)%t_sgn(:)

        do m=nst,ned
           sd = (m-nst)/3 + 1
           ord(m) = (t_ord(sd) - 1)*3 + (m-nst+1 - (m-nst)/3*3)
           sgn(m) = t_sgn(sd)
        end do

        call bc_extend_inward(s_st,s_ed,s_nd,s_lr,ngh,st,ed)


#ifdef PARALLEL
        packsize = product(ed(:)-st(:)+1)*(ned-nst+1)

        id_src = mb_topsp(nbs)%pid
        id_des = mb_topsp(nbt)%pid
        tag    = nint + ninters*3


        if (id_src /= id_des) then
            if (myid == id_src) then
                allocate(interssp(nint)%dat(st(1):ed(1), &
                                            st(2):ed(2), &
                                            st(3):ed(3), &
                                            nst:ned),stat=ierr)

                do m=nst,ned
                    do k=st(3),ed(3)
                    do j=st(2),ed(2)
                    do i=st(1),ed(1)
                        interssp(nint)%dat(i,j,k,m) = mb_sxyzsp(nbs)%fld(m)%r3d(i,j,k)
                    end do
                    end do
                    end do
                end do

                index_msg = index_msg + 1
                call MPI_ISEND(interssp(nint)%dat,packsize,kind_real_mpi, &
                              id_des,tag,MPI_COMM_WORLD,request_int(index_msg),ierr)
            end if

            if (myid == id_des) then
                allocate(interssp(nint)%dat(st(1):ed(1), &
                                            st(2):ed(2), &
                                            st(3):ed(3), &
                                            nst:ned),stat=ierr)

                index_msg = index_msg + 1
                call MPI_IRECV(interssp(nint)%dat,packsize,kind_real_mpi, &
                              id_src,tag,MPI_COMM_WORLD,request_int(index_msg),ierr)
            end if

        else
            if (myid == id_src) then
#endif
                allocate(interssp(nint)%dat(st(1):ed(1), &
                                            st(2):ed(2), &
                                            st(3):ed(3), &
                                            nst:ned),stat=ierr)

                do m=nst,ned
                    do k=st(3),ed(3)
                    do j=st(2),ed(2)
                    do i=st(1),ed(1)
                        interssp(nint)%dat(i,j,k,m) = mb_sxyzsp(nbs)%fld(m)%r3d(i,j,k)
                    end do
                    end do
                    end do
                end do


                do k=t_st(3),t_ed(3)
                do j=t_st(2),t_ed(2)
                do i=t_st(1),t_ed(1)
                    ijkt(:) = (/i,j,k/)
                    ijks(:) = mb_topsp(nbt)%bcs(nrt)%mapijk(i,j,k,:)
                    ijkgt(:) = ijkt(:)
                    ijkgs(:) = ijks(:)
                    do n=1,ngh
                        ijkgs(s_nd) = ijks(s_nd) - (n-1)*s_lr
                        is = ijkgs(1)
                        js = ijkgs(2)
                        ks = ijkgs(3)

                        ijkgt(t_nd) = ijkt(t_nd) + n*t_lr
                        it = ijkgt(1)
                        jt = ijkgt(2)
                        kt = ijkgt(3)

                        do m=nst,ned
                            td = ord(m)
                            mb_sxyzsp(nbt)%fld(m)%r3d(it,jt,kt) = sgn(m)*interssp(nint)%dat(is,js,ks,td)
                        end do
                    end do
                end do
                end do
                end do

                deallocate(interssp(nint)%dat,stat=ierr)
#ifdef PARALLEL
            end if
        end if
#endif
    end do

#ifdef PARALLEL
    if (nmsg_int > 0) then
        call MPI_WAITALL(nmsg_int, request_int(1:nmsg_int), status_int, ierr)

        do nint=1,ninters
            nbs     = interssp(nint)%bc%nbs
            nbt     = interssp(nint)%bc%nbt
            id_src = mb_topsp(nbs)%pid
            id_des = mb_topsp(nbt)%pid

            if (id_src /= id_des) then
                if (myid == id_src) then
                    deallocate(interssp(nint)%dat,stat=ierr)
                end if

                if (myid == id_des) then
                    s_nd    = interssp(nint)%bc%s_nd
                    s_lr    = interssp(nint)%bc%s_lr

                    nrt     = interssp(nint)%bc%nrt
                    t_st(:) = interssp(nint)%bc%t_st(:)
                    t_ed(:) = interssp(nint)%bc%t_ed(:)
                    t_nd    = interssp(nint)%bc%t_nd
                    t_lr    = interssp(nint)%bc%t_lr

                    t_ord(:) = mb_topsp(nbt)%bcs(nrt)%t_ord(:)
                    t_sgn(:) = mb_topsp(nbt)%bcs(nrt)%t_sgn(:)

                    do m=nst,ned
                       sd = (m-nst)/3 + 1
                       ord(m) = (t_ord(sd) - 1)*3 + (m-nst+1 - (m-nst)/3*3)
                       sgn(m) = t_sgn(sd)
                    end do

                    do k=t_st(3),t_ed(3)
                    do j=t_st(2),t_ed(2)
                    do i=t_st(1),t_ed(1)
                        ijkt(:) = (/i,j,k/)
                        ijks(:) = mb_topsp(nbt)%bcs(nrt)%mapijk(i,j,k,:)
                        ijkgt(:) = ijkt(:)
                        ijkgs(:) = ijks(:)
                        do n=1,ngh
                            ijkgs(s_nd) = ijks(s_nd) - (n-1)*s_lr
                            is = ijkgs(1)
                            js = ijkgs(2)
                            ks = ijkgs(3)

                            ijkgt(t_nd) = ijkt(t_nd) + n*t_lr
                            it = ijkgt(1)
                            jt = ijkgt(2)
                            kt = ijkgt(3)

                            do m=nst,ned
                                td = ord(m)
                                mb_sxyzsp(nbt)%fld(m)%r3d(it,jt,kt) = sgn(m)*interssp(nint)%dat(is,js,ks,td)
                            end do
                        end do
                    end do
                    end do
                    end do

                    deallocate(interssp(nint)%dat,stat=ierr)
                end if
            end if

        end do

    end if
#endif
end subroutine exchange_sp_bc_sxyz

subroutine exchange_singulars(mb_var,nst,ned,nswmem,nswave)
    use mod_kndconsts, only : kind_int
    use mod_constants, only : one,nsgl_aver_art
    use mod_datatypes, only : var_block_t
    use mod_fieldvars, only : mb_top,mb_vol
    use mod_fieldvars, only : nprocs
    use mod_singulars
    use mod_parallels
    implicit none
    type(var_block_t), pointer, intent(in) :: mb_var(:)
    integer(kind_int),          intent(in) :: nst,ned
    integer(kind_int),          intent(in) :: nswmem
    integer(kind_int),          intent(in) :: nswave
    integer(kind_int) :: i,j,pid,m,nvar,ierr
    integer(kind_int) :: npnts,nb,is,js,ks,ip,nmsg
    real(kind_real)   :: vol
#ifdef PARALLEL
    integer(kind_int) :: msg_count, msg_size, msg_tag
#endif

    nmsg = sum(nmsg_proc(:, 2))           ! send_recv 消息总数
    if (nmsg == 0) return ! 没有奇点的进程，直接返回
    nvar = simp_buf(nswmem, 0)%nvar


    ierr = ned - nst + 1 - nvar
    call error_check(ierr, &
                     "The size of array isn't nvar in subroutine exchange_singulars")


    ! step 1: 收集所有要发送的数据: simp_arr --> send_buf
    do i = 1, nsimp_arr
        npnts = simp_arr(i)%npnts
        do j=1,npnts
            nb  = simp_arr(i)%pnts(1,j)
            pid = mb_top(nb)%pid
#ifdef PARALLEL
            if (pid == myid) then
#endif
                is = simp_arr(i)%pnts(2,j)
                js = simp_arr(i)%pnts(3,j)
                ks = simp_arr(i)%pnts(4,j)

                ip = simp_arr(i)%pnts(5,j)
                if (nswave == nsgl_aver_art) then
                    vol = one
                else
                    vol = mb_vol(nb)%fld(1)%r3d(is,js,ks)
                end if
                do m=1,nvar
                    simp_buf(nswmem, pid)%buf(m,ip) = mb_var(nb)%fld(nst+m-1)%r3d(is,js,ks)/vol
                end do
#ifdef PARALLEL
            end if
#endif
        end do
    end do

    ! step 2: 发送消息: send_buf --> recv_buf


#ifdef PARALLEL
    msg_count = 0

    do i = 0, nprocs-1
        if (nmsg_proc(i, 2) > 0) then
            if (myid /= i) then
                msg_count = msg_count + 1
                msg_size  = nmsg_proc(i, 1) * nvar
                msg_tag   = myid
                call mpi_irecv( simp_buf(nswmem, i)%buf, msg_size, kind_real_mpi, i, msg_tag, &
                                MPI_COMM_WORLD, request_sin(msg_count), ierr)
            end if
        end if
    end do

    npnts = nmsg_proc(myid, 1)
    do i = 0, nprocs-1
        if (nmsg_proc(i, 2) > 0) then
            if (myid /= i) then
                msg_count = msg_count + 1
                msg_size  = npnts * nvar
                msg_tag   = i   ! 取接收者PID作为tag
                call mpi_isend( simp_buf(nswmem, myid)%buf, msg_size, kind_real_mpi, i, msg_tag, &
                                MPI_COMM_WORLD, request_sin(msg_count), ierr)
            end if
        end if
    end do

    if (msg_count /= nmsg_sin) then
        write(*, *) "!! <ERROR> !! msg_count /= nmsg_sin : ", msg_count, nmsg_sin
    end if

#endif

end subroutine exchange_singulars

subroutine exchange_singulars_cc(mb_var,nst,ned,nswmem,nswave)
    use mod_kndconsts, only : kind_int
    use mod_constants, only : one,nsgl_aver_art
    use mod_datatypes, only : var_block_t
    use mod_fieldvars, only : mb_topc,mb_volcc
    use mod_fieldvars, only : nprocs
    use mod_singulars
    use mod_parallels
    implicit none
    type(var_block_t), pointer, intent(in) :: mb_var(:)
    integer(kind_int),          intent(in) :: nst,ned
    integer(kind_int),          intent(in) :: nswmem
    integer(kind_int),          intent(in) :: nswave
    integer(kind_int) :: i,j,pid,m,nvar,ierr
    integer(kind_int) :: npnts,nb,is,js,ks,ip,nmsg
    real(kind_real)   :: vol
#ifdef PARALLEL
    integer(kind_int) :: msg_count, msg_size, msg_tag
#endif

    nmsg = sum(nmsg_proc(:, 2))           ! send_recv 消息总数
    if (nmsg == 0) return ! 没有奇点的进程，直接返回
    nvar = simp_buf(nswmem, 0)%nvar


    ierr = ned - nst + 1 - nvar
    call error_check(ierr, &
                     "The size of array isn't nvar in subroutine exchange_singulars")


    ! step 1: 收集所有要发送的数据: simp_arr --> send_buf
    do i = 1, nsimp_arr
        npnts = simp_arr(i)%npnts
        do j=1,npnts
            nb  = simp_arr(i)%pnts(1,j)
            pid = mb_topc(nb)%pid
#ifdef PARALLEL
            if (pid == myid) then
#endif
                is = simp_arr(i)%pnts(2,j)
                js = simp_arr(i)%pnts(3,j)
                ks = simp_arr(i)%pnts(4,j)

                ip = simp_arr(i)%pnts(5,j)
                if (nswave == nsgl_aver_art) then
                    vol = one
                else
                    vol = mb_volcc(nb)%fld(1)%r3d(is,js,ks)
                end if
                do m=1,nvar
                    simp_buf(nswmem, pid)%buf(m,ip) = mb_var(nb)%fld(nst+m-1)%r3d(is,js,ks)/vol
                end do
#ifdef PARALLEL
            end if
#endif
        end do
    end do

    ! step 2: 发送消息: send_buf --> recv_buf


#ifdef PARALLEL
    msg_count = 0

    do i = 0, nprocs-1
        if (nmsg_proc(i, 2) > 0) then
            if (myid /= i) then
                msg_count = msg_count + 1
                msg_size  = nmsg_proc(i, 1) * nvar
                msg_tag   = myid
                call mpi_irecv( simp_buf(nswmem, i)%buf, msg_size, kind_real_mpi, i, msg_tag, &
                                MPI_COMM_WORLD, request_sin(msg_count), ierr)
            end if
        end if
    end do

    npnts = nmsg_proc(myid, 1)
    do i = 0, nprocs-1
        if (nmsg_proc(i, 2) > 0) then
            if (myid /= i) then
                msg_count = msg_count + 1
                msg_size  = npnts * nvar
                msg_tag   = i   ! 取接收者PID作为tag
                call mpi_isend( simp_buf(nswmem, myid)%buf, msg_size, kind_real_mpi, i, msg_tag, &
                                MPI_COMM_WORLD, request_sin(msg_count), ierr)
            end if
        end if
    end do

    if (msg_count /= nmsg_sin) then
        write(*, *) "!! <ERROR> !! msg_count /= nmsg_sin : ", msg_count, nmsg_sin
    end if

#endif

end subroutine exchange_singulars_cc

subroutine average_singulars(mb_var,nst,ned,nswmem,nswave)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,one,nsgl_aver_art
    use mod_datatypes, only : var_block_t
    use mod_fieldvars, only : mb_top,mb_vol
    use mod_singulars
    use mod_parallels
    implicit none
    type(var_block_t), pointer, intent(in) :: mb_var(:)
    integer(kind_int),          intent(in) :: nst,ned
    integer(kind_int),          intent(in) :: nswmem
    integer(kind_int),          intent(in) :: nswave
    integer(kind_int)        :: i,j,k,pid,m,nvar,ierr
    integer(kind_int)        :: npnts,nb,nbt,is,js,ks,pdt,ipl
    real(kind_real)          :: vol,fsw,sumfsw
    real(kind_real), pointer :: var(:)

    nvar = simp_buf(nswmem, 0)%nvar

    ierr = ned - nst + 1 - nvar
    call error_check(ierr, &
                     "The size of array isn't nvar in subroutine exchange_singulars")

#ifdef PARALLEL
    ! step 0: 完成通信
    if (nmsg_sin > 0) then
        call mpi_waitall(nmsg_sin, request_sin(1:nmsg_sin), status_sin, ierr)
    end if
#endif

    allocate(var(nvar), stat=ierr)

    do i=1,nsimp_arr
        npnts = simp_arr(i)%npnts
        do j=1,npnts
            nb = simp_arr(i)%pnts(1,j)
            pid = mb_top(nb)%pid
#ifdef PARALLEL
            if (pid == myid) then
#endif
                is = simp_arr(i)%pnts(2,j)
                js = simp_arr(i)%pnts(3,j)
                ks = simp_arr(i)%pnts(4,j)

                if (nswave == nsgl_aver_art) then
                    vol = one
                else
                    vol = mb_vol(nb)%fld(1)%r3d(is,js,ks)
                end if

                sumfsw = zero
                var(:) = zero
                do k=1,npnts
                    nbt = simp_arr(i)%pnts(1,k)
                    pdt = mb_top(nbt)%pid
                    ipl = simp_arr(i)%pnts(5,k)
                    fsw = real( simp_arr(i)%pnts(6,k) )

                    sumfsw = sumfsw + fsw
                    var(:) = var(:) + fsw * simp_buf(nswmem, pdt)%buf(:, ipl)

                end do

                do m=1,nvar
                    var(m) = (var(m) / sumfsw) * vol
                    mb_var(nb)%fld(nst+m-1)%r3d(is,js,ks) = var(m)
                end do
#ifdef PARALLEL
            end if
#endif
        end do
    end do

    deallocate(var, stat=ierr)

end subroutine average_singulars

subroutine average_singulars_cc(mb_var,nst,ned,nswmem,nswave)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,one,nsgl_aver_art
    use mod_datatypes, only : var_block_t
    use mod_fieldvars, only : mb_topc,mb_volcc
    use mod_singulars
    use mod_parallels
    implicit none
    type(var_block_t), pointer, intent(in) :: mb_var(:)
    integer(kind_int),          intent(in) :: nst,ned
    integer(kind_int),          intent(in) :: nswmem
    integer(kind_int),          intent(in) :: nswave
    integer(kind_int)        :: i,j,k,pid,m,nvar,ierr
    integer(kind_int)        :: npnts,nb,nbt,is,js,ks,pdt,ipl
    real(kind_real)          :: vol,fsw,sumfsw
    real(kind_real), pointer :: var(:)

    nvar = simp_buf(nswmem, 0)%nvar

    ierr = ned - nst + 1 - nvar
    call error_check(ierr, &
                     "The size of array isn't nvar in subroutine exchange_singulars")

#ifdef PARALLEL
    ! step 0: 完成通信
    if (nmsg_sin > 0) then
        call mpi_waitall(nmsg_sin, request_sin(1:nmsg_sin), status_sin, ierr)
    end if
#endif

    allocate(var(nvar), stat=ierr)

    do i=1,nsimp_arr
        npnts = simp_arr(i)%npnts
        do j=1,npnts
            nb = simp_arr(i)%pnts(1,j)
            pid = mb_topc(nb)%pid
#ifdef PARALLEL
            if (pid == myid) then
#endif
                is = simp_arr(i)%pnts(2,j)
                js = simp_arr(i)%pnts(3,j)
                ks = simp_arr(i)%pnts(4,j)

                if (nswave == nsgl_aver_art) then
                    vol = one
                else
                    vol = mb_volcc(nb)%fld(1)%r3d(is,js,ks)
                end if

                sumfsw = zero
                var(:) = zero
                do k=1,npnts
                    nbt = simp_arr(i)%pnts(1,k)
                    pdt = mb_topc(nbt)%pid
                    ipl = simp_arr(i)%pnts(5,k)
                    fsw = real( simp_arr(i)%pnts(6,k) )

                    sumfsw = sumfsw + fsw
                    var(:) = var(:) + fsw * simp_buf(nswmem, pdt)%buf(:, ipl)

                end do

                do m=1,nvar
                    var(m) = (var(m) / sumfsw) * vol
                    mb_var(nb)%fld(nst+m-1)%r3d(is,js,ks) = var(m)
                end do
#ifdef PARALLEL
            end if
#endif
        end do
    end do

    deallocate(var, stat=ierr)

end subroutine average_singulars_cc

subroutine exchange_wall_pnts(bc_target,mark)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : tenth
    use mod_datatypes, only : fld_array_t,bc_region_t
    use mod_fieldvars, only : nblocks,nprocs,mb_top,mb_xyz
    use mod_fieldvars, only : nblkcoms,blkcoms
    use mod_fieldvars, only : nwallpnts,wallpnts
    use mod_parallels
    implicit none
    integer(kind_int), intent(in) :: bc_target,mark
    integer(kind_int)          :: nc,nb,pid,nr,i,j,k,m,n,ip,ierr
    integer(kind_int)          :: nwps,nwpnts,s_st(3),s_ed(3)
    integer(kind_int)          :: nregs,s_nd,s_lr,bctype,bcmark
    integer(kind_int), pointer :: nwpnts_local(:),nwists_local(:)
    integer(kind_int),parameter:: nwdat_max=7
    real(kind_real)  , pointer :: wpnts_local(:,:)
    real(kind_real)  , pointer :: wpnts_total(:,:)
    type(bc_region_t), pointer :: reg
    type(fld_array_t), pointer :: xyz(:)
    logical                    :: istar
#ifdef PARALLEL
    integer(kind_int) :: sendcount,recvcounts(numprocs)
    integer(kind_int) :: displs(numprocs)
#endif

    allocate(nwpnts_local(0:nprocs-1),stat=ierr)

    nwpnts_local(:) = 0

    nwallpnts = 0
    do nb=1,nblocks
        pid = mb_top(nb)%pid

        nregs = mb_top(nb)%nregions
        do nr=1,nregs
            reg => mb_top(nb)%bcs(nr)

            bctype = reg%bctype
            istar = (bctype == bc_target)
            if (mark > 0) then
                bcmark = reg%mark
                istar = istar .and. (bcmark == mark)
            end if

            if (istar) then
                s_st(:) = reg%s_st(:)
                s_ed(:) = reg%s_ed(:)
                s_nd    = reg%s_nd
                s_lr    = reg%s_lr
                nwps = product(s_ed(:)-s_st(:)+1)

                nwpnts_local(pid) = nwpnts_local(pid) + nwps

                nwallpnts = nwallpnts + nwps
            end if
        end do
    end do

    allocate(nwists_local(0:nprocs-1), stat=ierr)
    nwists_local(0) = 0
    do n=1,nprocs-1
       nwists_local(n) = nwists_local(n-1) + nwpnts_local(n-1)
    end do

    nwpnts = 0
    do nc=1,nblkcoms
        nb = blkcoms(nc)%nb

        nregs = mb_top(nb)%nregions
        do nr=1,nregs
            reg => mb_top(nb)%bcs(nr)

            bctype = reg%bctype
            istar = (bctype == bc_target)
            if (mark > 0) then
                bcmark = reg%mark
                istar = istar .and. (bcmark == mark)
            end if

            if (istar) then
                s_st(:) = reg%s_st(:)
                s_ed(:) = reg%s_ed(:)
                s_nd    = reg%s_nd
                s_lr    = reg%s_lr
                nwps = product(s_ed(:)-s_st(:)+1)

                nwpnts = nwpnts + nwps
            end if
        end do
    end do

    allocate(wpnts_local(nwdat_max,nwpnts),stat=ierr)

    ip = 0
    do nc=1,nblkcoms
        nb = blkcoms(nc)%nb

        xyz => mb_xyz(nb)%fld

        nregs = mb_top(nb)%nregions
        do nr=1,nregs
            reg => mb_top(nb)%bcs(nr)

            bctype = reg%bctype
            istar = (bctype == bc_target)
            if (mark > 0) then
                bcmark = reg%mark
                istar = istar .and. (bcmark == mark)
            end if

            if (istar) then
                s_st(:) = reg%s_st(:)
                s_ed(:) = reg%s_ed(:)
                s_nd    = reg%s_nd
                s_lr    = reg%s_lr

                do k=s_st(3),s_ed(3)
                do j=s_st(2),s_ed(2)
                do i=s_st(1),s_ed(1)
                    ip = ip + 1
                    wpnts_local(1,ip) = nb + tenth
                    wpnts_local(2,ip) = i  + tenth
                    wpnts_local(3,ip) = j  + tenth
                    wpnts_local(4,ip) = k  + tenth
                    wpnts_local(5,ip) = xyz(1)%r3d(i,j,k)
                    wpnts_local(6,ip) = xyz(2)%r3d(i,j,k)
                    wpnts_local(7,ip) = xyz(3)%r3d(i,j,k)
                end do
                end do
                end do
            end if
        end do
    end do

    allocate(wpnts_total(nwdat_max,nwallpnts),stat=ierr)

#ifdef PARALLEL
    if (nwallpnts > 0) then
        sendcount  = nwpnts_local(myid)*nwdat_max
        recvcounts = nwpnts_local(:)*nwdat_max
        displs(:)  = nwists_local(:)*nwdat_max

        call MPI_ALLGATHERV(wpnts_local,sendcount,kind_real_mpi, &
                            wpnts_total,recvcounts,displs,kind_real_mpi,MPI_COMM_WORLD,ierr)
    end if
#else
    do ip=1,nwallpnts
        do m=1,nwdat_max
            wpnts_total(m,ip) = wpnts_local(m,ip)
        end do
    end do
#endif
    deallocate(wpnts_local,stat=ierr)

    allocate(wallpnts(nwallpnts),stat=ierr)

    do ip=1,nwallpnts
        wallpnts(ip)%nb = floor( wpnts_total(1,ip) )
        wallpnts(ip)%i  = floor( wpnts_total(2,ip) )
        wallpnts(ip)%j  = floor( wpnts_total(3,ip) )
        wallpnts(ip)%k  = floor( wpnts_total(4,ip) )
        wallpnts(ip)%xw = wpnts_total(5,ip)
        wallpnts(ip)%yw = wpnts_total(6,ip)
        wallpnts(ip)%zw = wpnts_total(7,ip)
    end do

    deallocate(wpnts_total,stat=ierr)

end subroutine exchange_wall_pnts
