
subroutine reset_file_status
    use mod_constants, only : nrestrt_restart
    use mod_constants, only : io_unit_res,io_unit_fce
    use mod_variables, only : nrestrt,resfile,fcefile
    implicit none
    character(len=128) :: header

    if (nrestrt == nrestrt_restart) then
        header = 'variables="NSTEP","CFL","DT","RESAVE","RESMAX","NB","I","J","K","NV","CPU Time","Wall Time"'
        call openfile(io_unit_res,resfile,"unknown","formatted","sequential")
        write(io_unit_res,'(a128)') adjustl(header)
        call closefile(io_unit_res,resfile)

        header = 'variables="NSTEP","Cfx","Cfy","Cfz","Cmx","Cmy","Cmz","Cd","Cl","Xcp","CPU Time","Wall Time"'
        call openfile(io_unit_fce,fcefile,"unknown","formatted","sequential")
        write(io_unit_fce,'(a128)') adjustl(header)
        call closefile(io_unit_fce,fcefile)
    end if

end subroutine reset_file_status


subroutine output_res
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,io_unit_res
    use mod_variables, only : cflmax,dtaumax
    use mod_variables, only : nstep,nressav,resfile
    use mod_variables, only : respos,resmax,restot
    use mod_fieldvars, only : ntotpts
    use mod_runtimers, only : time_cpu,time_wall
    use mod_parallels
    implicit none
    integer(kind_int) :: m,ierr
    real(kind_real)   :: resave,cflmax0,dtaumax0

#ifdef PARALLEL
    call MPI_REDUCE(cflmax,cflmax0,1,kind_real_mpi,MPI_MAX,master,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(cflmax0,1,kind_real_mpi,master,MPI_COMM_WORLD,ierr)
    call MPI_REDUCE(dtaumax,dtaumax0,1,kind_real_mpi,MPI_MAX,master,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(dtaumax0,1,kind_real_mpi,master,MPI_COMM_WORLD,ierr)
    cflmax  = cflmax0
    dtaumax = dtaumax0

    if (myid == master) then
#endif
        resave = sqrt(restot/ntotpts)

        if (mod(nstep-nressav,10*nressav) == 0) then
            write(*,*)
            write(*,10)"NSTEP","CFL","DT","RESAVE","RESMAX","NB","I","J","K","NV"
        end if

        call openfile(io_unit_res,resfile,"unknown","formatted","append")
        write(*          ,20)nstep,cflmax,dtaumax,resave,resmax,(respos(m),m=1,5)
        write(io_unit_res,30)nstep,cflmax,dtaumax,resave,resmax,(respos(m),m=1,5),time_cpu,time_wall
        call closefile(io_unit_res,resfile)
#ifdef PARALLEL
    end if
#endif

10  format(a7,2(1x,a9  ),2(1x,a12  ),5(1x,a4))
20  format(i7,2(1x,e9.2),2(1x,e12.5),5(1x,i4))
30  format(i7,2(1x,e9.2),2(1x,e12.5),5(1x,i4),2(1x,e12.5))

end subroutine output_res

subroutine output_fce
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : one,small,io_unit_fce
    use mod_variables, only : nstep,nfcesav,fcefile
    use mod_variables, only : reflgrd,frefsc,freflen
    use mod_variables, only : attack,ayaw
    use mod_variables, only : cfx,cfy,cfz,cmx,cmy,cmz
    use mod_runtimers, only : time_cpu,time_wall
    use mod_parallels
    implicit none
    integer(kind_int) :: ierr
    real(kind_real)   :: forces(6),pforces(6),cl,cd,xcp
    real(kind_real)   :: scf,vcm,sina,cosa,sinb,cosb

#ifdef PARALLEL
    pforces(1) = cfx
    pforces(2) = cfy
    pforces(3) = cfz
    pforces(4) = cmx
    pforces(5) = cmy
    pforces(6) = cmz
    call MPI_REDUCE(pforces,forces,6,kind_real_mpi,MPI_SUM,master,MPI_COMM_WORLD,ierr)
    cfx = forces(1)
    cfy = forces(2)
    cfz = forces(3)
    cmx = forces(4)
    cmy = forces(5)
    cmz = forces(6)
#endif

#ifdef PARALLEL
    if (myid == master) then
#endif
        scf = reflgrd*reflgrd/frefsc
        vcm = reflgrd**3/(frefsc*freflen)

        cfx = cfx * scf
        cfy = cfy * scf
        cfz = cfz * scf
        cmx = cmx * vcm
        cmy = cmy * vcm
        cmz = cmz * vcm

        sina = sin(attack)
        cosa = cos(attack)
        sinb = sin(ayaw)
        cosb = cos(ayaw)

        cl = cfy * cosa - cfx * sina
        cd = cfy * sina * cosb + cfx * cosa * cosb + cfz * sinb
        xcp = sign(one,cfy)*cmz/(abs(cfy) + small)

        if (mod(nstep,100*nfcesav) == 0) then
            write(*,*)
            write(*,'(1x,3(1x,a12  ))') "cfx","cfy","cfz"
            write(*,'(1x,3(1x,e12.5))') cfx,cfy,cfz
            write(*,'(1x,3(1x,a12  ))') "cmx","cmy","cmz"
            write(*,'(1x,3(1x,e12.5))') cmx,cmy,cmz
            write(*,'(1x,3(1x,a12  ))') "cd","cl","xcp"
            write(*,'(1x,3(1x,e12.5))') cd,cl,xcp
        end if

        call openfile(io_unit_fce,fcefile,"unknown","formatted","append")
        write(io_unit_fce,10) nstep,cfx,cfy,cfz,cmx,cmy,cmz,cd,cl,xcp,time_cpu,time_wall
        call closefile(io_unit_fce,fcefile)
#ifdef PARALLEL
    end if
#endif

10  format(1x,i6,11(1x,e12.5))

end subroutine output_fce

subroutine output_sol
    use mod_kndconsts, only : kind_int
    use mod_constants, only : io_unit_fld,io_unit_tur
    use mod_constants, only : io_unit_aux,nvis_ns_lam 
    use mod_variables, only : solfile,cfl,nvis,nstep,nrestrt,nsolsav
    use mod_variables, only : nsub,nsubp,nsubp1
    use mod_fieldvars, only : nstepmean,nsteprms,ntime
    use mod_runtimers, only : time_cpu,time_wall
    use mod_parallels
    implicit none
    integer(kind_int) :: error,ierr,check

    call openfile_check(io_unit_fld,trim(solfile),"unknown","unformatted","stream",ierr)
    
#ifdef PARALLEL
    call MPI_BCAST(ierr,1,kind_int_mpi,master,MPI_COMM_WORLD,error)
#endif

    if(ierr == 0) then
        check = check
    else
        check = check + 1
    end if
    
    call openfile_check(io_unit_aux,trim(solfile)//".aux","unknown","unformatted","stream",ierr)
    
#ifdef PARALLEL
    call MPI_BCAST(ierr,1,kind_int_mpi,master,MPI_COMM_WORLD,error)
#endif

    if(ierr == 0) then
        check = check
    else
        check = check + 1
    end if    
    
    if (check == 0) then
        call openfile(io_unit_fld,trim(solfile)//".bak","unknown","unformatted","stream")
        call deletefile(io_unit_fld,trim(solfile)//".bak")
        call openfile(io_unit_aux,trim(solfile)//".aux"//".bak","unknown","formatted","sequential")
        call deletefile(io_unit_aux,trim(solfile)//".aux"//".bak")
        
        call renamefile(io_unit_fld,solfile)
        call renamefile(io_unit_aux,trim(solfile)//".aux")
    end if
    
    call openfile(io_unit_fld,solfile,"unknown","unformatted","stream")
    call deletefile(io_unit_fld,solfile)
    
    call openfile(io_unit_fld,solfile,"unknown","unformatted","stream")
    call output_sol_plot3d(io_unit_fld)
    call closefile(io_unit_fld,solfile) 
    
    call openfile(io_unit_aux,trim(solfile)//".aux","unknown","formatted","sequential")
    call run_seq_and_master(output_sol_aux)
    call closefile(io_unit_aux,trim(solfile)//".aux")

    if (nvis > nvis_ns_lam) then
        call openfile_check(io_unit_tur,trim(solfile)//".tur","unknown","unformatted","stream",ierr)
    
#ifdef PARALLEL
        call MPI_BCAST(ierr,1,kind_int_mpi,master,MPI_COMM_WORLD,error)
#endif        
        if (ierr == 0) then
            call openfile(io_unit_tur,trim(solfile)//".tur"//".bak","unknown","unformatted","stream")
            call deletefile(io_unit_tur,trim(solfile)//".tur"//".bak")
            
            call renamefile(io_unit_tur,trim(solfile)//".tur")  
        end if
       
        call openfile(io_unit_tur,trim(solfile)//".tur","unknown","unformatted","stream")
        call output_tur_sol_plot3d(io_unit_tur)
        call closefile(io_unit_tur,trim(solfile)//".tur")
    end if

    call msg_seq_and_master("output the flow field successfully")

    contains

    subroutine output_sol_aux
        write(io_unit_aux,*) nstep,nstepmean,nsteprms,ntime,cfl
        write(io_unit_aux,*) nsub,nsubp,nsubp1,time_cpu,time_wall
    end subroutine output_sol_aux

end subroutine output_sol

subroutine output_sol_sp
    use mod_kndconsts, only : kind_int
    use mod_constants, only : io_unit_fld,io_unit_tur
    use mod_constants, only : io_unit_aux,nvis_ns_lam 
    use mod_variables, only : solfile,cfl,nvis,nstep,nrestrt,nsolsav
    use mod_variables, only : nsub,nsubp,nsubp1
    use mod_fieldvars, only : nstepmean,nsteprms,ntime
    use mod_runtimers, only : time_cpu,time_wall
    use mod_parallels
    implicit none
    integer(kind_int) :: error,ierr,check

    call openfile_check(io_unit_fld,trim(solfile),"unknown","unformatted","stream",ierr)
    
#ifdef PARALLEL
    call MPI_BCAST(ierr,1,kind_int_mpi,master,MPI_COMM_WORLD,error)
#endif

    if(ierr == 0) then
        check = check
    else
        check = check + 1
    end if
    
    call openfile_check(io_unit_aux,trim(solfile)//".aux","unknown","unformatted","stream",ierr)
    
#ifdef PARALLEL
    call MPI_BCAST(ierr,1,kind_int_mpi,master,MPI_COMM_WORLD,error)
#endif

    if(ierr == 0) then
        check = check
    else
        check = check + 1
    end if    
    
    if (check == 0) then
        call openfile(io_unit_fld,trim(solfile)//".bak","unknown","unformatted","stream")
        call deletefile(io_unit_fld,trim(solfile)//".bak")
        call openfile(io_unit_aux,trim(solfile)//".aux"//".bak","unknown","formatted","sequential")
        call deletefile(io_unit_aux,trim(solfile)//".aux"//".bak")
        
        call renamefile(io_unit_fld,solfile)
        call renamefile(io_unit_aux,trim(solfile)//".aux")
    end if
    
    call openfile(io_unit_fld,solfile,"unknown","unformatted","stream")
    call deletefile(io_unit_fld,solfile)
    
    call openfile(io_unit_fld,solfile,"unknown","unformatted","stream")
    call output_sol_plot3d_sp(io_unit_fld)
    call closefile(io_unit_fld,solfile) 
    
    call openfile(io_unit_aux,trim(solfile)//".aux","unknown","formatted","sequential")
    call run_seq_and_master(output_sol_aux)
    call closefile(io_unit_aux,trim(solfile)//".aux")

    if (nvis > nvis_ns_lam) then
        call openfile_check(io_unit_tur,trim(solfile)//".tur","unknown","unformatted","stream",ierr)
    
#ifdef PARALLEL
        call MPI_BCAST(ierr,1,kind_int_mpi,master,MPI_COMM_WORLD,error)
#endif        
        if (ierr == 0) then
            call openfile(io_unit_tur,trim(solfile)//".tur"//".bak","unknown","unformatted","stream")
            call deletefile(io_unit_tur,trim(solfile)//".tur"//".bak")
            
            call renamefile(io_unit_tur,trim(solfile)//".tur")  
        end if
       
        call openfile(io_unit_tur,trim(solfile)//".tur","unknown","unformatted","stream")
        call output_tur_sol_plot3d_sp(io_unit_tur)
        call closefile(io_unit_tur,trim(solfile)//".tur")
    end if

    call msg_seq_and_master("output the flow field successfully")

    contains

    subroutine output_sol_aux
        write(io_unit_aux,*) nstep,nstepmean,nsteprms,ntime,cfl
        write(io_unit_aux,*) nsub,nsubp,nsubp1,time_cpu,time_wall
    end subroutine output_sol_aux

end subroutine output_sol_sp

subroutine output_plt
    use mod_variables, only : pltsty
    implicit none

    select case(pltsty)
    case(0)
       call output_surface
    case(1)
       call output_bcwall
    case(2)
       call output_surface
       call output_bcwall
    case(3)
       call output_surface
       call output_bcwall
       call output_flowfield
    case(4)
       call output_surface
       call output_bcwall
       call output_flowfield_unsteady
    case(5)
       call output_surface
       call output_bcwall
       call output_vortfield_unsteady
    end select

    call msg_seq_and_master("output the tecout flow successfully")

end subroutine output_plt

subroutine output_plt_sp
    use mod_variables, only : pltsty
    implicit none

    select case(pltsty)
    case(0)
       call output_surface_sp
    case(1)
       call output_bcwall
    case(2)
       call output_surface
       call output_bcwall
    case(3)
       call output_surface
       call output_bcwall
       call output_flowfield
    case(4)
       call output_surface
       call output_bcwall
       call output_flowfield_unsteady
    case(5)
       call output_surface
       call output_bcwall
       call output_vortfield_unsteady
    end select

    call msg_seq_and_master("output the tecout flow successfully")

end subroutine output_plt_sp

subroutine output_dst
    use mod_constants, only : io_unit_dst
    use mod_variables, only : grdfile
    implicit none

    call openfile(io_unit_dst,trim(grdfile)//".dst","unknown","unformatted","stream")
    call output_dst_plot3d(io_unit_dst)
    call closefile(io_unit_dst,trim(grdfile)//".dst")

end subroutine output_dst

subroutine output_dsp
    use mod_constants, only : io_unit_dsp
    use mod_variables, only : grdfile
    implicit none

    call openfile(io_unit_dsp,trim(grdfile)//".dsp","unknown","unformatted","stream")
    call output_dsp_plot3d(io_unit_dsp)
    call closefile(io_unit_dsp,trim(grdfile)//".dsp")

end subroutine output_dsp

subroutine output_sol_plot3d(io_unit)
    use mod_kndconsts, only : kind_int
    use mod_constants, only : nplot3d_fun
    use mod_fieldvars, only : mb_pv,npvs
    use mod_interface, only : gather_output_mb_var
    implicit none
    integer(kind_int), intent(in) :: io_unit

    call gather_output_mb_var(io_unit,mb_pv,1,npvs,nplot3d_fun)

end subroutine output_sol_plot3d

subroutine output_sol_plot3d_sp(io_unit)
    use mod_kndconsts, only : kind_int
    use mod_constants, only : nplot3d_fun
    use mod_fieldvars, only : mb_pv,npvs
    use mod_interface, only : gather_output_mb_var_sp
    implicit none
    integer(kind_int), intent(in) :: io_unit

    call gather_output_mb_var_sp(io_unit,mb_pv,1,npvs,nplot3d_fun)

end subroutine output_sol_plot3d_sp

subroutine output_tur_sol_plot3d(io_unit)
    use mod_kndconsts, only : kind_int
    use mod_constants, only : nplot3d_fun
    use mod_fieldvars, only : mb_qvst,nqvst
    use mod_interface, only : gather_output_mb_var
    implicit none
    integer(kind_int), intent(in) :: io_unit

    call gather_output_mb_var(io_unit,mb_qvst,1,nqvst,nplot3d_fun)

end subroutine output_tur_sol_plot3d

subroutine output_tur_sol_plot3d_sp(io_unit)
    use mod_kndconsts, only : kind_int
    use mod_constants, only : nplot3d_fun
    use mod_fieldvars, only : mb_qvst,nqvst
    use mod_interface, only : gather_output_mb_var_sp
    implicit none
    integer(kind_int), intent(in) :: io_unit

    call gather_output_mb_var_sp(io_unit,mb_qvst,1,nqvst,nplot3d_fun)

end subroutine output_tur_sol_plot3d_sp

subroutine output_dst_plot3d(io_unit)
    use mod_kndconsts, only : kind_int
    use mod_constants, only : nplot3d_fun
    use mod_fieldvars, only : mb_dst
    use mod_interface, only : gather_output_mb_var
    implicit none
    integer(kind_int), intent(in) :: io_unit

    call gather_output_mb_var(io_unit,mb_dst,1,2,nplot3d_fun)

end subroutine output_dst_plot3d

subroutine output_dsp_plot3d(io_unit)
    use mod_kndconsts, only : kind_int
    use mod_constants, only : nplot3d_fun
    use mod_fieldvars, only : mb_dsp
    use mod_interface, only : gather_output_mb_var
    implicit none
    integer(kind_int), intent(in) :: io_unit

    call gather_output_mb_var(io_unit,mb_dsp,1,2,nplot3d_fun)

end subroutine output_dsp_plot3d

subroutine output_surface
    use mod_constants, only : io_unit_tec
    use mod_variables, only : pltfile
    implicit none

    call openfile(io_unit_tec,pltfile,"unknown","unformatted","stream")
    call deletefile(io_unit_tec,pltfile)
    
    call openfile(io_unit_tec,pltfile,"unknown","unformatted","stream")
    call tecout_bin_surface(io_unit_tec)
    call closefile(io_unit_tec,pltfile)

end subroutine output_surface

subroutine output_surface_sp
    use mod_constants, only : io_unit_tec
    use mod_variables, only : pltfile
    implicit none

    call openfile(io_unit_tec,pltfile,"unknown","unformatted","stream")
    call deletefile(io_unit_tec,pltfile)
    
    call openfile(io_unit_tec,pltfile,"unknown","unformatted","stream")
    call tecout_bin_surface_sp(io_unit_tec)
    call closefile(io_unit_tec,pltfile)

end subroutine output_surface_sp

subroutine output_surface_mean
    use mod_kndconsts, only : kind_int,len_char_file
    use mod_constants, only : io_unit_tec
    use mod_variables, only : pltfile
    implicit none
    integer(kind_int)        :: m
    character(len_char_file) :: tecfile

    m = index(pltfile,'.',back=.TRUE.)
    if (m == 0) then
       tecfile = trim(adjustl(pltfile)) // "_mean.plt"
    else
       tecfile = trim(adjustl(pltfile(1:m-1))) // "_mean.plt"
    end if
    
    call openfile(io_unit_tec,tecfile,"unknown","unformatted","stream")
    call deletefile(io_unit_tec,tecfile)
    
    call openfile(io_unit_tec,tecfile,"unknown","unformatted","stream")
    call tecout_bin_surface(io_unit_tec)
    call closefile(io_unit_tec,tecfile)

end subroutine output_surface_mean

subroutine output_surface_rms
    use mod_kndconsts, only : kind_int,len_char_file
    use mod_constants, only : io_unit_tec
    use mod_variables, only : pltfile
    implicit none
    integer(kind_int)        :: m
    character(len_char_file) :: tecfile

    m = index(pltfile,'.',back=.TRUE.)
    if (m == 0) then
       tecfile = trim(adjustl(pltfile)) // "_rms.plt"
    else
       tecfile = trim(adjustl(pltfile(1:m-1))) // "_rms.plt"
    end if
    
    call openfile(io_unit_tec,tecfile,"unknown","unformatted","stream")
    call deletefile(io_unit_tec,tecfile)
    
    call openfile(io_unit_tec,tecfile,"unknown","unformatted","stream")
    call tecout_bin_surface_rms(io_unit_tec)
    call closefile(io_unit_tec,tecfile)

end subroutine output_surface_rms

subroutine tecout_bin_surface(io_unit)
    use mod_kndconsts, only : kind_int,kind_real,kind_single
    use mod_constants, only : sml_ssf,one,m3x3,bc_wall
    use mod_constants, only : nvis_ns_lam,nvis_tur_sa
    use mod_constants, only : nvis_tur_sst,nvis_tur_hst
    use mod_datatypes, only : var_block_t,fld_array_t
    use mod_datatypes, only : top_block_t,bc_region_t
    use mod_variables, only : gamma,refbeta,pltfile,nvis
    use mod_fieldvars, only : nblocks,mb_top,mb_xyz,mb_sxyz
    use mod_fieldvars, only : npvs,mb_pv,nqvst,mb_qvst
    use mod_parallels
    use mod_tecplotio
    implicit none
    integer(kind_int), intent(in) :: io_unit
    integer(kind_int)          :: nb,nr,i,j,k,m,m1,m2,m3,ierr
    integer(kind_int)          :: s_st(3),s_ed(3),st(3),ed(3)
    integer(kind_int)          :: nregs,s_nd,s_lr,bctype,subtype
    integer(kind_int)          :: ni,nj,nk,ijk2(3),is,js,ks
    real(kind_real)            :: gama,ro,u,v,w,p,c2,t,ma
    real(kind_real)            :: nx,ny,nz,on,vn
    type(top_block_t), pointer :: top
    type(bc_region_t), pointer :: reg
    type(fld_array_t), pointer :: xyz(:),sxyz(:),pv(:),qvst(:)
    real(kind_real),   pointer :: vars(:,:,:,:)
    integer(kind_int)          :: ntecvars,ndata
    character(len=256)         :: tecvarnames
    character(len=128)         :: zonename
    character(len=32 )         :: str1,str2,str3
#ifdef PARALLEL
    integer(kind_int) :: pid,packsize,tag
    integer(kind_int) :: status(MPI_STATUS_SIZE)
#endif

    gama = gamma

    ntecvars = 10
    tecvarnames = "X,Y,Z,Ro,U,V,W,P,T,Ma"

    if (nvis > nvis_ns_lam) then
        ntecvars = ntecvars + nqvst
        select case(nvis)
        case(nvis_tur_sa)
            tecvarnames = "X,Y,Z,Ro,U,V,W,P,T,Ma,Nu,Vist"
        case(nvis_tur_sst,nvis_tur_hst)
            tecvarnames = "X,Y,Z,Ro,U,V,W,P,T,Ma,Tke,Ome,Vist"
        end select
    end if

    if (kind_real == kind_single) then
        ndata = 1
    else
        ndata = 2
    end if

#ifdef PARALLEL
    if (myid == master) then
#endif
        call tecio_ini(io_unit,"HOSTA",ntecvars,tecvarnames)

        do nb=1,nblocks
            top => mb_top(nb)

            nregs = top%nregions
            do nr=1,nregs
                reg => top%bcs(nr)

                bctype  = reg%bctype
                subtype = reg%subtype
                if ( isoutbc(bctype,subtype) ) then
                    s_st(:) = reg%s_st(:)
                    s_ed(:) = reg%s_ed(:)
                    s_nd    = reg%s_nd
                    s_lr    = reg%s_lr

                    st(:) = s_st(:)
                    ed(:) = s_ed(:)

                    ni = ed(1) - st(1) + 1
                    nj = ed(2) - st(2) + 1
                    nk = ed(3) - st(3) + 1

                    write(str1,'(i6)') nb
                    write(str2,'(i6)') nr
                    write(str3,'(i6)') bctype

                    zonename = 'BLK'//trim(adjustl(str1))// &
                               'BC'//trim(adjustl(str2))// &
                               'T'//trim(adjustl(str3))

                    call tecio_zone(io_unit,trim(zonename),0,ni,nj,nk)
                end if
            end do
        end do

        call tecio_eohmark(io_unit)
#ifdef PARALLEL
    end if
#endif

    do nb=1,nblocks
        top => mb_top(nb)

        xyz  => mb_xyz(nb)%fld
        sxyz => mb_sxyz(nb)%fld
        pv   => mb_pv(nb)%fld

        if (nvis > nvis_ns_lam) then
            qvst => mb_qvst(nb)%fld
        end if

        nregs = top%nregions
        do nr=1,nregs
            reg => top%bcs(nr)

            bctype  = reg%bctype
            subtype = reg%subtype
            if ( isoutbc(bctype,subtype) ) then
                s_st(:) = reg%s_st(:)
                s_ed(:) = reg%s_ed(:)
                s_nd    = reg%s_nd
                s_lr    = reg%s_lr

                st(:) = s_st(:)
                ed(:) = s_ed(:)

                m1 = m3x3(1,s_nd)
                m2 = m3x3(2,s_nd)
                m3 = m3x3(3,s_nd)

#ifdef PARALLEL
                packsize = product(ed(:)-st(:)+1)*ntecvars

                tag = nb*1000+nr

                pid = mb_top(nb)%pid
                if ((pid /= master .and. myid == pid) .or. myid == master) then
#endif
                    allocate(vars(st(1):ed(1),st(2):ed(2),st(3):ed(3),ntecvars), stat=ierr)
#ifdef PARALLEL
                end if

                if (myid == pid) then
#endif
                    do k=st(3),ed(3)
                    do j=st(2),ed(2)
                    do i=st(1),ed(1)
                        ro = pv(1)%r3d(i,j,k)
                        u  = pv(2)%r3d(i,j,k)
                        v  = pv(3)%r3d(i,j,k)
                        w  = pv(4)%r3d(i,j,k)
                        p  = pv(5)%r3d(i,j,k)
                        t  = p/(refbeta*ro)
                        c2 = gama*p/ro
                        ma = sqrt((u*u+v*v+w*w)/c2)

                        if (bctype == bc_wall) then
                            nx = sxyz(m1)%r3d(i,j,k)
                            ny = sxyz(m2)%r3d(i,j,k)
                            nz = sxyz(m3)%r3d(i,j,k)
                            on = one/max(sqrt(nx*nx + ny*ny + nz*nz),sml_ssf)
                            nx = nx*on
                            ny = ny*on
                            nz = nz*on
                            ijk2(:) = (/i,j,k/)
                            ijk2(s_nd) = ijk2(s_nd) - s_lr
                            is = ijk2(1)
                            js = ijk2(2)
                            ks = ijk2(3)
                            u  = pv(2)%r3d(is,js,ks)
                            v  = pv(3)%r3d(is,js,ks)
                            w  = pv(4)%r3d(is,js,ks)
                            vn = nx*u + ny*v + nz*w
                            u  = u - nx*vn
                            v  = v - ny*vn
                            w  = w - nz*vn
                        end if

                        vars(i,j,k, 1) = xyz(1)%r3d(i,j,k)
                        vars(i,j,k, 2) = xyz(2)%r3d(i,j,k)
                        vars(i,j,k, 3) = xyz(3)%r3d(i,j,k)
                        vars(i,j,k, 4) = ro
                        vars(i,j,k, 5) = u
                        vars(i,j,k, 6) = v
                        vars(i,j,k, 7) = w
                        vars(i,j,k, 8) = p
                        vars(i,j,k, 9) = t
                        vars(i,j,k,10) = ma

                        if (nvis > nvis_ns_lam) then
                            do m=1,nqvst
                                vars(i,j,k,m+10) = qvst(m)%r3d(i,j,k)
                            end do
                        end if
                    end do
                    end do
                    end do
#ifdef PARALLEL
                end if

                if (pid /= master) then
                    if (myid == pid) then
                        call MPI_SEND(vars,packsize,kind_real_mpi, &
                                      master,tag,MPI_COMM_WORLD,ierr)
                    end if

                    if (myid == master) then
                        call MPI_RECV(vars,packsize,kind_real_mpi, &
                                      pid,tag,MPI_COMM_WORLD,status,ierr)
                    end if
                end if

                if (myid == master) then
#endif
                   call tecio_data(io_unit,ndata)

                   write(io_unit)((((vars(i,j,k,m),i=st(1),ed(1)), &
                                                   j=st(2),ed(2)), &
                                                   k=st(3),ed(3)), &
                                                   m=1,ntecvars)
#ifdef PARALLEL
                end if

                if ((pid /= master .and. myid == pid) .or. myid == master) then
#endif
                    deallocate(vars, stat=ierr)
#ifdef PARALLEL
                end if
#endif
            end if
        end do
    end do

    contains

    function isoutbc(bct,subct) result(isout)
        use mod_constants, only : bc_inflow,bc_farfield,bc_cut1to1
        implicit none
        integer(kind_int) :: bct,subct
        logical           :: isout

        if ( bct /= bc_inflow   .and. &
             bct /= bc_farfield .and. &
             bct /= bc_cut1to1 ) then
            isout = .true.
        else
            isout = .false.
        end if

    end function isoutbc

end subroutine tecout_bin_surface

subroutine tecout_bin_surface_sp(io_unit)
    use mod_kndconsts, only : kind_int,kind_real,kind_single
    use mod_constants, only : sml_ssf,one,m3x3,bc_wall
    use mod_constants, only : nvis_ns_lam,nvis_tur_sa
    use mod_constants, only : nvis_tur_sst,nvis_tur_hst
    use mod_datatypes, only : var_block_t,fld_array_t
    use mod_datatypes, only : top_block_t,bc_region_t
    use mod_variables, only : gamma,refbeta,pltfile,nvis
    use mod_fieldvars, only : nblocks,mb_topsp,mb_xyzsp,mb_sxyzsp !,mb_topc,mb_xyzcc,mb_sxyzcc!
    use mod_fieldvars, only : npvs,mb_pv,nqvst,mb_qvst
    use mod_parallels
    use mod_tecplotio
    implicit none
    integer(kind_int), intent(in) :: io_unit
    integer(kind_int)          :: nb,nr,i,j,k,m,m1,m2,m3,ierr
    integer(kind_int)          :: s_st(3),s_ed(3),st(3),ed(3)
    integer(kind_int)          :: nregs,s_nd,s_lr,bctype,subtype
    integer(kind_int)          :: ni,nj,nk,ijk2(3),is,js,ks
    real(kind_real)            :: gama,ro,u,v,w,p,c2,t,ma
    real(kind_real)            :: nx,ny,nz,on,vn
    type(top_block_t), pointer :: top
    type(bc_region_t), pointer :: reg
    type(fld_array_t), pointer :: xyz(:),sxyz(:),pv(:),qvst(:)
    real(kind_real),   pointer :: vars(:,:,:,:)
    integer(kind_int)          :: ntecvars,ndata
    character(len=256)         :: tecvarnames
    character(len=128)         :: zonename
    character(len=32 )         :: str1,str2,str3
#ifdef PARALLEL
    integer(kind_int) :: pid,packsize,tag
    integer(kind_int) :: status(MPI_STATUS_SIZE)
#endif

    gama = gamma

    ntecvars = 10
    tecvarnames = "X,Y,Z,Ro,U,V,W,P,T,Ma"

    if (nvis > nvis_ns_lam) then
        ntecvars = ntecvars + nqvst
        select case(nvis)
        case(nvis_tur_sa)
            tecvarnames = "X,Y,Z,Ro,U,V,W,P,T,Ma,Nu,Vist"
        case(nvis_tur_sst,nvis_tur_hst)
            tecvarnames = "X,Y,Z,Ro,U,V,W,P,T,Ma,Tke,Ome,Vist"
        end select
    end if

    if (kind_real == kind_single) then
        ndata = 1
    else
        ndata = 2
    end if

#ifdef PARALLEL
    if (myid == master) then
#endif
        call tecio_ini(io_unit,"HOSTA",ntecvars,tecvarnames)

        do nb=1,nblocks
            top => mb_topsp(nb)

            nregs = top%nregions
            do nr=1,nregs
                reg => top%bcs(nr)

                bctype  = reg%bctype
                subtype = reg%subtype
                if ( isoutbcsp(bctype,subtype) ) then
                    s_st(:) = reg%s_st(:)
                    s_ed(:) = reg%s_ed(:)
                    s_nd    = reg%s_nd
                    s_lr    = reg%s_lr

                    st(:) = s_st(:)
                    ed(:) = s_ed(:)

                    ni = ed(1) - st(1) + 1
                    nj = ed(2) - st(2) + 1
                    nk = ed(3) - st(3) + 1

                    write(str1,'(i6)') nb
                    write(str2,'(i6)') nr
                    write(str3,'(i6)') bctype

                    zonename = 'BLK'//trim(adjustl(str1))// &
                               'BC'//trim(adjustl(str2))// &
                               'T'//trim(adjustl(str3))

                    call tecio_zone(io_unit,trim(zonename),0,ni,nj,nk)
                end if
            end do
        end do

        call tecio_eohmark(io_unit)
#ifdef PARALLEL
    end if
#endif

    do nb=1,nblocks
        top => mb_topsp(nb)

        xyz  => mb_xyzsp(nb)%fld
        sxyz => mb_sxyzsp(nb)%fld
        pv   => mb_pv(nb)%fld

        if (nvis > nvis_ns_lam) then
            qvst => mb_qvst(nb)%fld
        end if

        nregs = top%nregions
        do nr=1,nregs
            reg => top%bcs(nr)

            bctype  = reg%bctype
            subtype = reg%subtype
            if ( isoutbcsp(bctype,subtype) ) then
                s_st(:) = reg%s_st(:)
                s_ed(:) = reg%s_ed(:)
                s_nd    = reg%s_nd
                s_lr    = reg%s_lr

                st(:) = s_st(:)
                ed(:) = s_ed(:)

                m1 = m3x3(1,s_nd)
                m2 = m3x3(2,s_nd)
                m3 = m3x3(3,s_nd)

#ifdef PARALLEL
                packsize = product(ed(:)-st(:)+1)*ntecvars

                tag = nb*1000+nr

                pid = mb_topsp(nb)%pid
                if ((pid /= master .and. myid == pid) .or. myid == master) then
#endif
                    allocate(vars(st(1):ed(1),st(2):ed(2),st(3):ed(3),ntecvars), stat=ierr)
#ifdef PARALLEL
                end if

                if (myid == pid) then
#endif
                    do k=st(3),ed(3)
                    do j=st(2),ed(2)
                    do i=st(1),ed(1)
                        ro = pv(1)%r3d(i,j,k)
                        u  = pv(2)%r3d(i,j,k)
                        v  = pv(3)%r3d(i,j,k)
                        w  = pv(4)%r3d(i,j,k)
                        p  = pv(5)%r3d(i,j,k)
                        t  = p/(refbeta*ro)
                        c2 = gama*p/ro
                        ma = sqrt((u*u+v*v+w*w)/c2)

                        !if (bctype == bc_wall) then
                        !    nx = sxyz(m1)%r3d(i,j,k)
                        !    ny = sxyz(m2)%r3d(i,j,k)
                        !    nz = sxyz(m3)%r3d(i,j,k)
                        !    on = one/max(sqrt(nx*nx + ny*ny + nz*nz),sml_ssf)
                        !    nx = nx*on
                        !    ny = ny*on
                        !    nz = nz*on
                        !    ijk2(:) = (/i,j,k/)
                        !    ijk2(s_nd) = ijk2(s_nd) - s_lr
                        !    is = ijk2(1)
                        !    js = ijk2(2)
                        !    ks = ijk2(3)
                        !    u  = pv(2)%r3d(is,js,ks)
                        !    v  = pv(3)%r3d(is,js,ks)
                        !    w  = pv(4)%r3d(is,js,ks)
                        !    vn = nx*u + ny*v + nz*w
                        !    u  = u - nx*vn
                        !    v  = v - ny*vn
                        !    w  = w - nz*vn
                        !end if

                        vars(i,j,k, 1) = xyz(1)%r3d(i,j,k)
                        vars(i,j,k, 2) = xyz(2)%r3d(i,j,k)
                        vars(i,j,k, 3) = xyz(3)%r3d(i,j,k)
                        vars(i,j,k, 4) = ro
                        vars(i,j,k, 5) = u
                        vars(i,j,k, 6) = v
                        vars(i,j,k, 7) = w
                        vars(i,j,k, 8) = p
                        vars(i,j,k, 9) = t
                        vars(i,j,k,10) = ma

                        if (nvis > nvis_ns_lam) then
                            do m=1,nqvst
                                vars(i,j,k,m+10) = qvst(m)%r3d(i,j,k)
                            end do
                        end if
                    end do
                    end do
                    end do
#ifdef PARALLEL
                end if

                if (pid /= master) then
                    if (myid == pid) then
                        call MPI_SEND(vars,packsize,kind_real_mpi, &
                                      master,tag,MPI_COMM_WORLD,ierr)
                    end if

                    if (myid == master) then
                        call MPI_RECV(vars,packsize,kind_real_mpi, &
                                      pid,tag,MPI_COMM_WORLD,status,ierr)
                    end if
                end if

                if (myid == master) then
#endif
                   call tecio_data(io_unit,ndata)

                   write(io_unit)((((vars(i,j,k,m),i=st(1),ed(1)), &
                                                   j=st(2),ed(2)), &
                                                   k=st(3),ed(3)), &
                                                   m=1,ntecvars)
#ifdef PARALLEL
                end if

                if ((pid /= master .and. myid == pid) .or. myid == master) then
#endif
                    deallocate(vars, stat=ierr)
#ifdef PARALLEL
                end if
#endif
            end if
        end do
    end do

    contains

    function isoutbcsp(bct,subct) result(isout)
        use mod_constants, only : bc_inflow,bc_farfield,bc_cut1to1
        implicit none
        integer(kind_int) :: bct,subct
        logical           :: isout

        if ( bct /= bc_inflow   .and. &
             bct /= bc_farfield .and. &
             bct /= bc_cut1to1 ) then
            isout = .true.
        else
            isout = .false.
        end if

    end function isoutbcsp

end subroutine tecout_bin_surface_sp

subroutine tecout_bin_surface_sp1(io_unit)
    use mod_kndconsts, only : kind_int,kind_real,kind_single
    use mod_constants, only : sml_ssf,one,m3x3,bc_wall
    use mod_constants, only : nvis_ns_lam,nvis_tur_sa
    use mod_constants, only : nvis_tur_sst,nvis_tur_hst
    use mod_datatypes, only : var_block_t,fld_array_t
    use mod_datatypes, only : top_block_t,bc_region_t
    use mod_variables, only : gamma,refbeta,pltfile,nvis
    use mod_fieldvars, only : nblocks,mb_topc,mb_xyzcc,mb_sxyzcc!,mb_topsp,mb_xyzsp,mb_sxyzsp !
    use mod_fieldvars, only : npvs,mb_pvfp,nqvst,mb_qvst
    use mod_parallels
    use mod_tecplotio
    implicit none
    integer(kind_int), intent(in) :: io_unit
    integer(kind_int)          :: nb,nr,i,j,k,m,m1,m2,m3,ierr
    integer(kind_int)          :: s_st(3),s_ed(3),st(3),ed(3)
    integer(kind_int)          :: nregs,s_nd,s_lr,bctype,subtype
    integer(kind_int)          :: ni,nj,nk,ijk2(3),is,js,ks
    real(kind_real)            :: gama,ro,u,v,w,p,c2,t,ma
    real(kind_real)            :: nx,ny,nz,on,vn
    type(top_block_t), pointer :: top
    type(bc_region_t), pointer :: reg
    type(fld_array_t), pointer :: xyz(:),sxyz(:),pv(:),qvst(:)
    real(kind_real),   pointer :: vars(:,:,:,:)
    integer(kind_int)          :: ntecvars,ndata
    character(len=256)         :: tecvarnames
    character(len=128)         :: zonename
    character(len=32 )         :: str1,str2,str3
#ifdef PARALLEL
    integer(kind_int) :: pid,packsize,tag
    integer(kind_int) :: status(MPI_STATUS_SIZE)
#endif

    gama = gamma

    ntecvars = 10
    tecvarnames = "X,Y,Z,Ro,U,V,W,P,T,Ma"

    if (nvis > nvis_ns_lam) then
        ntecvars = ntecvars + nqvst
        select case(nvis)
        case(nvis_tur_sa)
            tecvarnames = "X,Y,Z,Ro,U,V,W,P,T,Ma,Nu,Vist"
        case(nvis_tur_sst,nvis_tur_hst)
            tecvarnames = "X,Y,Z,Ro,U,V,W,P,T,Ma,Tke,Ome,Vist"
        end select
    end if

    if (kind_real == kind_single) then
        ndata = 1
    else
        ndata = 2
    end if
    
    call map_sptocp

#ifdef PARALLEL
    if (myid == master) then
#endif
        call tecio_ini(io_unit,"HOSTA",ntecvars,tecvarnames)

        do nb=1,nblocks
            top => mb_topc(nb)

            nregs = top%nregions
            do nr=1,nregs
                reg => top%bcs(nr)

                bctype  = reg%bctype
                subtype = reg%subtype
                if ( isoutbcsp(bctype,subtype) ) then
                    s_st(:) = reg%s_st(:)
                    s_ed(:) = reg%s_ed(:)
                    s_nd    = reg%s_nd
                    s_lr    = reg%s_lr

                    st(:) = s_st(:)
                    ed(:) = s_ed(:)

                    ni = ed(1) - st(1) + 1
                    nj = ed(2) - st(2) + 1
                    nk = ed(3) - st(3) + 1

                    write(str1,'(i6)') nb
                    write(str2,'(i6)') nr
                    write(str3,'(i6)') bctype

                    zonename = 'BLK'//trim(adjustl(str1))// &
                               'BC'//trim(adjustl(str2))// &
                               'T'//trim(adjustl(str3))

                    call tecio_zone(io_unit,trim(zonename),0,ni,nj,nk)
                end if
            end do
        end do

        call tecio_eohmark(io_unit)
#ifdef PARALLEL
    end if
#endif

    do nb=1,nblocks
        top => mb_topc(nb)

        xyz  => mb_xyzcc(nb)%fld
        sxyz => mb_sxyzcc(nb)%fld
        pv   => mb_pvfp(nb)%fld

        if (nvis > nvis_ns_lam) then
            qvst => mb_qvst(nb)%fld
        end if

        nregs = top%nregions
        do nr=1,nregs
            reg => top%bcs(nr)

            bctype  = reg%bctype
            subtype = reg%subtype
            if ( isoutbcsp(bctype,subtype) ) then
                s_st(:) = reg%s_st(:)
                s_ed(:) = reg%s_ed(:)
                s_nd    = reg%s_nd
                s_lr    = reg%s_lr

                st(:) = s_st(:)
                ed(:) = s_ed(:)

                m1 = m3x3(1,s_nd)
                m2 = m3x3(2,s_nd)
                m3 = m3x3(3,s_nd)

#ifdef PARALLEL
                packsize = product(ed(:)-st(:)+1)*ntecvars

                tag = nb*1000+nr

                pid = mb_topc(nb)%pid
                if ((pid /= master .and. myid == pid) .or. myid == master) then
#endif
                    allocate(vars(st(1):ed(1),st(2):ed(2),st(3):ed(3),ntecvars), stat=ierr)
#ifdef PARALLEL
                end if

                if (myid == pid) then
#endif
                    do k=st(3),ed(3)
                    do j=st(2),ed(2)
                    do i=st(1),ed(1)
                        ro = pv(1)%r3d(i,j,k)
                        u  = pv(2)%r3d(i,j,k)
                        v  = pv(3)%r3d(i,j,k)
                        w  = pv(4)%r3d(i,j,k)
                        p  = pv(5)%r3d(i,j,k)
                        t  = p/(refbeta*ro)
                        c2 = gama*p/ro
                        ma = sqrt((u*u+v*v+w*w)/c2)

                        !if (bctype == bc_wall) then
                        !    nx = sxyz(m1)%r3d(i,j,k)
                        !    ny = sxyz(m2)%r3d(i,j,k)
                        !    nz = sxyz(m3)%r3d(i,j,k)
                        !    on = one/max(sqrt(nx*nx + ny*ny + nz*nz),sml_ssf)
                        !    nx = nx*on
                        !    ny = ny*on
                        !    nz = nz*on
                        !    ijk2(:) = (/i,j,k/)
                        !    ijk2(s_nd) = ijk2(s_nd) - s_lr
                        !    is = ijk2(1)
                        !    js = ijk2(2)
                        !    ks = ijk2(3)
                        !    u  = pv(2)%r3d(is,js,ks)
                        !    v  = pv(3)%r3d(is,js,ks)
                        !    w  = pv(4)%r3d(is,js,ks)
                        !    vn = nx*u + ny*v + nz*w
                        !    u  = u - nx*vn
                        !    v  = v - ny*vn
                        !    w  = w - nz*vn
                        !end if

                        vars(i,j,k, 1) = xyz(1)%r3d(i,j,k)
                        vars(i,j,k, 2) = xyz(2)%r3d(i,j,k)
                        vars(i,j,k, 3) = xyz(3)%r3d(i,j,k)
                        vars(i,j,k, 4) = ro
                        vars(i,j,k, 5) = u
                        vars(i,j,k, 6) = v
                        vars(i,j,k, 7) = w
                        vars(i,j,k, 8) = p
                        vars(i,j,k, 9) = t
                        vars(i,j,k,10) = ma

                        if (nvis > nvis_ns_lam) then
                            do m=1,nqvst
                                vars(i,j,k,m+10) = qvst(m)%r3d(i,j,k)
                            end do
                        end if
                    end do
                    end do
                    end do
#ifdef PARALLEL
                end if

                if (pid /= master) then
                    if (myid == pid) then
                        call MPI_SEND(vars,packsize,kind_real_mpi, &
                                      master,tag,MPI_COMM_WORLD,ierr)
                    end if

                    if (myid == master) then
                        call MPI_RECV(vars,packsize,kind_real_mpi, &
                                      pid,tag,MPI_COMM_WORLD,status,ierr)
                    end if
                end if

                if (myid == master) then
#endif
                   call tecio_data(io_unit,ndata)

                   write(io_unit)((((vars(i,j,k,m),i=st(1),ed(1)), &
                                                   j=st(2),ed(2)), &
                                                   k=st(3),ed(3)), &
                                                   m=1,ntecvars)
#ifdef PARALLEL
                end if

                if ((pid /= master .and. myid == pid) .or. myid == master) then
#endif
                    deallocate(vars, stat=ierr)
#ifdef PARALLEL
                end if
#endif
            end if
        end do
    end do

    contains

    function isoutbcsp(bct,subct) result(isout)
        use mod_constants, only : bc_inflow,bc_farfield,bc_cut1to1
        implicit none
        integer(kind_int) :: bct,subct
        logical           :: isout

        if ( bct /= bc_inflow   .and. &
             bct /= bc_farfield .and. &
             bct /= bc_cut1to1 ) then
            isout = .true.
        else
            isout = .false.
        end if

    end function isoutbcsp

end subroutine tecout_bin_surface_sp1
    
subroutine tecout_bin_surface_rms(io_unit)
    use mod_kndconsts, only : kind_int,kind_real,kind_single
    use mod_constants, only : sml_ssf,one,m3x3,bc_wall
    use mod_constants, only : nvis_ns_lam,nvis_tur_sa
    use mod_constants, only : nvis_tur_sst,nvis_tur_hst
    use mod_datatypes, only : var_block_t,fld_array_t
    use mod_datatypes, only : top_block_t,bc_region_t
    use mod_variables, only : gamma,refbeta,pltfile,nvis
    use mod_fieldvars, only : nblocks,mb_top,mb_xyz,mb_sxyz
    use mod_fieldvars, only : npvs,mb_pv,nqvst,mb_qvst
    use mod_parallels
    use mod_tecplotio
    implicit none
    integer(kind_int), intent(in) :: io_unit
    integer(kind_int)          :: nb,nr,i,j,k,m,m1,m2,m3,ierr
    integer(kind_int)          :: s_st(3),s_ed(3),st(3),ed(3)
    integer(kind_int)          :: nregs,s_nd,s_lr,bctype,subtype
    integer(kind_int)          :: ni,nj,nk,ijk2(3),is,js,ks
    real(kind_real)            :: ro,u,v,w,p
    real(kind_real)            :: nx,ny,nz,on,vn
    type(top_block_t), pointer :: top
    type(bc_region_t), pointer :: reg
    type(fld_array_t), pointer :: xyz(:),pv(:)
    real(kind_real),   pointer :: vars(:,:,:,:)
    integer(kind_int)          :: ntecvars,ndata
    character(len=256)         :: tecvarnames
    character(len=128)         :: zonename
    character(len=32 )         :: str1,str2,str3
#ifdef PARALLEL
    integer(kind_int) :: pid,packsize,tag
    integer(kind_int) :: status(MPI_STATUS_SIZE)
#endif

    ntecvars = 8
    tecvarnames = "X,Y,Z,Ro,U,V,W,P"

    if (kind_real == kind_single) then
        ndata = 1
    else
        ndata = 2
    end if

#ifdef PARALLEL
    if (myid == master) then
#endif
        call tecio_ini(io_unit,"HOSTA",ntecvars,tecvarnames)

        do nb=1,nblocks
            top => mb_top(nb)

            nregs = top%nregions
            do nr=1,nregs
                reg => top%bcs(nr)

                bctype  = reg%bctype
                subtype = reg%subtype
                if ( isoutbc(bctype,subtype) ) then
                    s_st(:) = reg%s_st(:)
                    s_ed(:) = reg%s_ed(:)
                    s_nd    = reg%s_nd
                    s_lr    = reg%s_lr

                    st(:) = s_st(:)
                    ed(:) = s_ed(:)

                    ni = ed(1) - st(1) + 1
                    nj = ed(2) - st(2) + 1
                    nk = ed(3) - st(3) + 1

                    write(str1,'(i6)') nb
                    write(str2,'(i6)') nr
                    write(str3,'(i6)') bctype

                    zonename = 'BLK'//trim(adjustl(str1))// &
                               'BC'//trim(adjustl(str2))// &
                               'T'//trim(adjustl(str3))

                    call tecio_zone(io_unit,trim(zonename),0,ni,nj,nk)
                end if
            end do
        end do

        call tecio_eohmark(io_unit)
#ifdef PARALLEL
    end if
#endif

    do nb=1,nblocks
        top => mb_top(nb)

        xyz  => mb_xyz(nb)%fld
        pv   => mb_pv(nb)%fld

        nregs = top%nregions
        do nr=1,nregs
            reg => top%bcs(nr)

            bctype  = reg%bctype
            subtype = reg%subtype
            if ( isoutbc(bctype,subtype) ) then
                s_st(:) = reg%s_st(:)
                s_ed(:) = reg%s_ed(:)
                s_nd    = reg%s_nd
                s_lr    = reg%s_lr

                st(:) = s_st(:)
                ed(:) = s_ed(:)

                m1 = m3x3(1,s_nd)
                m2 = m3x3(2,s_nd)
                m3 = m3x3(3,s_nd)

#ifdef PARALLEL
                packsize = product(ed(:)-st(:)+1)*ntecvars

                tag = nb*1000+nr

                pid = mb_top(nb)%pid
                if ((pid /= master .and. myid == pid) .or. myid == master) then
#endif
                    allocate(vars(st(1):ed(1),st(2):ed(2),st(3):ed(3),ntecvars), stat=ierr)
#ifdef PARALLEL
                end if

                if (myid == pid) then
#endif
                    do k=st(3),ed(3)
                    do j=st(2),ed(2)
                    do i=st(1),ed(1)
                        ro = pv(1)%r3d(i,j,k)
                        u  = pv(2)%r3d(i,j,k)
                        v  = pv(3)%r3d(i,j,k)
                        w  = pv(4)%r3d(i,j,k)
                        p  = pv(5)%r3d(i,j,k)

                        vars(i,j,k, 1) = xyz(1)%r3d(i,j,k)
                        vars(i,j,k, 2) = xyz(2)%r3d(i,j,k)
                        vars(i,j,k, 3) = xyz(3)%r3d(i,j,k)
                        vars(i,j,k, 4) = ro
                        vars(i,j,k, 5) = u
                        vars(i,j,k, 6) = v
                        vars(i,j,k, 7) = w
                        vars(i,j,k, 8) = p
                    end do
                    end do
                    end do
#ifdef PARALLEL
                end if

                if (pid /= master) then
                    if (myid == pid) then
                        call MPI_SEND(vars,packsize,kind_real_mpi, &
                                      master,tag,MPI_COMM_WORLD,ierr)
                    end if

                    if (myid == master) then
                        call MPI_RECV(vars,packsize,kind_real_mpi, &
                                      pid,tag,MPI_COMM_WORLD,status,ierr)
                    end if
                end if

                if (myid == master) then
#endif
                   call tecio_data(io_unit,ndata)

                   write(io_unit)((((vars(i,j,k,m),i=st(1),ed(1)), &
                                                   j=st(2),ed(2)), &
                                                   k=st(3),ed(3)), &
                                                   m=1,ntecvars)
#ifdef PARALLEL
                end if

                if ((pid /= master .and. myid == pid) .or. myid == master) then
#endif
                    deallocate(vars, stat=ierr)
#ifdef PARALLEL
                end if
#endif
            end if
        end do
    end do

    contains

    function isoutbc(bct,subct) result(isout)
        use mod_constants, only : bc_inflow,bc_farfield,bc_cut1to1
        implicit none
        integer(kind_int) :: bct,subct
        logical           :: isout

        if ( bct /= bc_inflow   .and. &
             bct /= bc_farfield .and. &
             bct /= bc_cut1to1 ) then
            isout = .true.
        else
            isout = .false.
        end if

    end function isoutbc

end subroutine tecout_bin_surface_rms 

subroutine output_bcwall
    use mod_kndconsts, only : kind_int,len_char_file
    use mod_constants, only : io_unit_wfl,nvis_euler
    use mod_variables, only : nvis,pltfile
    implicit none
    integer(kind_int)        :: m
    character(len_char_file) :: tecfile

    if (nvis > nvis_euler) then
        m = index(pltfile,'.',back=.TRUE.)
        if (m == 0) then
           tecfile = trim(adjustl(pltfile)) // "_wall.plt"
        else
           tecfile = trim(adjustl(pltfile(1:m-1))) // "_wall.plt"
        end if

        call openfile(io_unit_wfl,tecfile,"unknown","unformatted","stream")
        call deletefile(io_unit_wfl,tecfile)
        
        call openfile(io_unit_wfl,tecfile,"unknown","unformatted","stream")
        call tecout_bin_bcwall(io_unit_wfl)
        call closefile(io_unit_wfl,tecfile)
    end if

end subroutine output_bcwall

subroutine output_bcwall_mean
    use mod_kndconsts, only : kind_int,len_char_file
    use mod_constants, only : io_unit_wfl,nvis_euler
    use mod_variables, only : nvis,pltfile
    implicit none
    integer(kind_int)        :: m
    character(len_char_file) :: tecfile
    
    if (nvis > nvis_euler) then
        m = index(pltfile,'.',back=.TRUE.)
        if (m == 0) then
           tecfile = trim(adjustl(pltfile)) // "_wall_mean.plt"
        else
           tecfile = trim(adjustl(pltfile(1:m-1))) // "_wall_mean.plt"
        end if

        call openfile(io_unit_wfl,tecfile,"unknown","unformatted","stream")
        call deletefile(io_unit_wfl,tecfile)
        
        call openfile(io_unit_wfl,tecfile,"unknown","unformatted","stream")
        call tecout_bin_bcwall(io_unit_wfl)
        call closefile(io_unit_wfl,tecfile)
    end if    

end subroutine output_bcwall_mean

subroutine tecout_bin_bcwall(io_unit)
    use mod_kndconsts, only : kind_int,kind_real,kind_single
    use mod_constants, only : half,one,two,two3rd,m3x3
    use mod_constants, only : sml_ssf,bc_wall
    use mod_datatypes, only : var_block_t,fld_array_t
    use mod_datatypes, only : top_block_t,bc_region_t
    use mod_variables, only : gamma,uoo,voo,woo,poo,reue
    use mod_variables, only : refbeta,refqw,prlam,prtur,pltfile
    use mod_fieldvars, only : nblocks,mb_top,mb_xyz,mb_sxyz
    use mod_fieldvars, only : mb_pv,mb_dpv,mb_vsl,mb_vst
    use mod_parallels
    use mod_tecplotio
    implicit none
    integer(kind_int), intent(in) :: io_unit
    integer(kind_int)          :: nb,nr,i,j,k,m,m1,m2,m3,ierr
    integer(kind_int)          :: s_st(3),s_ed(3),st(3),ed(3)
    integer(kind_int)          :: nregs,s_nd,s_lr,bctype,ni,nj,nk
    real(kind_real)            :: gama,ae,cp,cp_prl,cp_prt,re
    real(kind_real)            :: fp,vis,kcp,clr,clrre,nx,ny,nz,on
    real(kind_real)            :: dux,duy,duz,dvx,dvy,dvz
    real(kind_real)            :: dwx,dwy,dwz,dtx,dty,dtz
    real(kind_real)            :: vis2p3,tauxx,tauyy,tauzz,cf,sf
    real(kind_real)            :: tauxy,tauxz,tauyz,fvx,fvy,fvz,qw
    type(top_block_t), pointer :: top
    type(bc_region_t), pointer :: reg
    type(fld_array_t), pointer :: xyz(:),sxyz(:),dpv(:)
    type(fld_array_t), pointer :: vsl(:),vst(:),pv(:)
    real(kind_real),   pointer :: vars(:,:,:,:)
    integer(kind_int)          :: ntecvars,ndata
    character(len=256)         :: tecvarnames
    character(len=128)         :: zonename
    character(len=32 )         :: str1,str2,str3
#ifdef PARALLEL
    integer(kind_int) :: pid,packsize,tag
    integer(kind_int) :: status(MPI_STATUS_SIZE)
#endif

    gama = gamma
    ae = gama - one

    cp = gama*refbeta/ae
    cp_prl = cp/prlam
    cp_prt = cp/prtur

    re = one/reue

    ntecvars = 9
    tecvarnames = "X,Y,Z,cfx,cfy,cfz,cf,qw,cp"

    if (kind_real == kind_single) then
        ndata = 1
    else
        ndata = 2
    end if


#ifdef PARALLEL
    if (myid == master) then
#endif
        call tecio_ini(io_unit,"HOSTA",ntecvars,tecvarnames)

        do nb=1,nblocks
            top => mb_top(nb)

            nregs = top%nregions
            do nr=1,nregs
                reg => top%bcs(nr)

                bctype  = reg%bctype
                if (bctype == bc_wall) then
                    s_st(:) = reg%s_st(:)
                    s_ed(:) = reg%s_ed(:)
                    s_nd    = reg%s_nd
                    s_lr    = reg%s_lr

                    st(:) = s_st(:)
                    ed(:) = s_ed(:)

                    ni = ed(1) - st(1) + 1
                    nj = ed(2) - st(2) + 1
                    nk = ed(3) - st(3) + 1

                    write(str1,'(i6)') nb
                    write(str2,'(i6)') nr
                    write(str3,'(i6)') bctype

                    zonename = 'BLK'//trim(adjustl(str1))// &
                               'BC'//trim(adjustl(str2))// &
                               'T'//trim(adjustl(str3))

                    call tecio_zone(io_unit,trim(zonename),0,ni,nj,nk)
                end if
            end do
        end do

        call tecio_eohmark(io_unit)
#ifdef PARALLEL
    end if
#endif

    do nb=1,nblocks
        top => mb_top(nb)

        xyz  => mb_xyz(nb)%fld
        sxyz => mb_sxyz(nb)%fld
        vsl  => mb_vsl(nb)%fld
        vst  => mb_vst(nb)%fld
        pv   => mb_pv(nb)%fld
        dpv  => mb_dpv(nb)%fld

        nregs = top%nregions
        do nr=1,nregs
            reg => top%bcs(nr)

            bctype  = reg%bctype
            if (bctype == bc_wall) then
                s_st(:) = reg%s_st(:)
                s_ed(:) = reg%s_ed(:)
                s_nd    = reg%s_nd
                s_lr    = reg%s_lr

                st(:) = s_st(:)
                ed(:) = s_ed(:)

                m1 = m3x3(1,s_nd)
                m2 = m3x3(2,s_nd)
                m3 = m3x3(3,s_nd)

                clr   = two*s_lr
                clrre = -clr*re
#ifdef PARALLEL
                packsize = product(ed(:)-st(:)+1)*ntecvars

                tag = nb*1000+nr

                pid = mb_top(nb)%pid
                if ((pid /= master .and. myid == pid) .or. myid == master) then
#endif
                    allocate(vars(st(1):ed(1),st(2):ed(2),st(3):ed(3),ntecvars), stat=ierr)
#ifdef PARALLEL
                end if

                if (myid == pid) then
#endif
                    do k=st(3),ed(3)
                    do j=st(2),ed(2)
                    do i=st(1),ed(1)
                        dux = dpv(1)%r3d(i,j,k)
                        duy = dpv(2)%r3d(i,j,k)
                        duz = dpv(3)%r3d(i,j,k)

                        dvx = dpv(4)%r3d(i,j,k)
                        dvy = dpv(5)%r3d(i,j,k)
                        dvz = dpv(6)%r3d(i,j,k)

                        dwx = dpv(7)%r3d(i,j,k)
                        dwy = dpv(8)%r3d(i,j,k)
                        dwz = dpv(9)%r3d(i,j,k)

                        dtx = dpv(10)%r3d(i,j,k)
                        dty = dpv(11)%r3d(i,j,k)
                        dtz = dpv(12)%r3d(i,j,k)

                        vis = vsl(1)%r3d(i,j,k)
                        kcp = vis*cp_prl

                        vis = vis + vst(1)%r3d(i,j,k)
                        kcp = kcp + vst(1)%r3d(i,j,k)*cp_prt

                        nx = sxyz(m1)%r3d(i,j,k)
                        ny = sxyz(m2)%r3d(i,j,k)
                        nz = sxyz(m3)%r3d(i,j,k)
                        on = one/max(sqrt(nx*nx + ny*ny + nz*nz),sml_ssf)
                        nx = nx*on
                        ny = ny*on
                        nz = nz*on

!!                        vis2p3 = two3rd*vis
!!                        tauxx = vis2p3 * ( two*dux - dvy - dwz )
!!                        tauyy = vis2p3 * ( two*dvy - dwz - dux )
!!                        tauzz = vis2p3 * ( two*dwz - dux - dvy )
!!                        tauxy = vis * ( duy + dvx )
!!                        tauxz = vis * ( duz + dwx )
!!                        tauyz = vis * ( dvz + dwy )
!!                        fvx = clrre*(nx*tauxx + ny*tauxy + nz*tauxz)
!!                        fvy = clrre*(nx*tauxy + ny*tauyy + nz*tauyz)
!!                        fvz = clrre*(nx*tauxz + ny*tauyz + nz*tauzz)

                        tauxx = vis*(nx*dux + ny*duy + nz*duz)
                        tauyy = vis*(nx*dvx + ny*dvy + nz*dvz)
                        tauzz = vis*(nx*dwx + ny*dwy + nz*dwz)

                        tauxy = nx*tauxx + ny*tauyy + nz*tauzz

                        fvx = clrre*(tauxx - nx*tauxy)
                        fvy = clrre*(tauyy - ny*tauxy)
                        fvz = clrre*(tauzz - nz*tauxy)

                        !!sf = sign(one,fvx)
                        sf = sign(one,uoo*fvx + voo*fvy + woo*fvz)
                        cf = sf*sqrt(fvx*fvx + fvy*fvy + fvz*fvz)

                        qw = half*clrre*kcp*(nx*dtx + ny*dty + nz*dtz)

                        fp = two*(pv(5)%r3d(i,j,k) - poo)

                        vars(i,j,k,1) = xyz(1)%r3d(i,j,k)
                        vars(i,j,k,2) = xyz(2)%r3d(i,j,k)
                        vars(i,j,k,3) = xyz(3)%r3d(i,j,k)
                        vars(i,j,k,4) = fvx
                        vars(i,j,k,5) = fvy
                        vars(i,j,k,6) = fvz
                        vars(i,j,k,7) = cf
                        vars(i,j,k,8) = qw    !!*refqw
                        vars(i,j,k,9) = fp
                    end do
                    end do
                    end do
#ifdef PARALLEL
                end if

                if (pid /= master) then
                    if (myid == pid) then
                        call MPI_SEND(vars,packsize,kind_real_mpi, &
                                      master,tag,MPI_COMM_WORLD,ierr)
                    end if

                    if (myid == master) then
                        call MPI_RECV(vars,packsize,kind_real_mpi, &
                                      pid,tag,MPI_COMM_WORLD,status,ierr)
                    end if
                end if

                if (myid == master) then
#endif
                   call tecio_data(io_unit,ndata)

                   write(io_unit)((((vars(i,j,k,m),i=st(1),ed(1)), &
                                                   j=st(2),ed(2)), &
                                                   k=st(3),ed(3)), &
                                                   m=1,ntecvars)
#ifdef PARALLEL
                end if

                if ((pid /= master .and. myid == pid) .or. myid == master) then
#endif
                    deallocate(vars, stat=ierr)
#ifdef PARALLEL
                end if
#endif
            end if
        end do
    end do

end subroutine tecout_bin_bcwall

subroutine output_flowfield
    use mod_kndconsts, only : kind_int,len_char_file
    use mod_constants, only : io_unit_tfld
    use mod_variables, only : pltfile
    implicit none
    integer(kind_int)        :: m    
    character(len_char_file) :: tecfile
    
    m = index(pltfile,'.',back=.TRUE.)
    if (m == 0) then
       tecfile = trim(adjustl(pltfile)) // "_flow.plt"
    else
       tecfile = trim(adjustl(pltfile(1:m-1))) // "_flow.plt"
    end if
    
    call openfile(io_unit_tfld,tecfile,"unknown","unformatted","stream")
    call deletefile(io_unit_tfld,tecfile)
    
    call openfile(io_unit_tfld,tecfile,"unknown","unformatted","stream")
    call tecout_bin_flowfield(io_unit_tfld)
    call closefile(io_unit_tfld,tecfile)

end subroutine output_flowfield

subroutine output_flowfield_unsteady
    use mod_kndconsts, only : kind_int,len_char_file
    use mod_constants, only : io_unit_tfld
    use mod_variables, only : pltfile,nstep,ns_prms,nprms
    implicit none
    integer(kind_int)        :: m    
    character(len_char_file) :: tecfile
    character(len=64 )       :: str1
    
    write(str1,'(i6)') nstep    
    m = index(pltfile,'.',back=.TRUE.)    
    if (m == 0) then
        if(nprms>0 .and. nstep>ns_prms) then
            tecfile = trim(adjustl(pltfile)) //trim(adjustl(str1)) // "_flow.plt"
        else
            tecfile = trim(adjustl(pltfile)) // "_flow.plt"
        end if
    else
        if(nprms>0 .and. nstep>ns_prms) then
            tecfile = trim(adjustl(pltfile(1:m-1))) //trim(adjustl(str1)) // "_flow.plt"
        else
            tecfile = trim(adjustl(pltfile(1:m-1))) // "_flow.plt"
        endif
    end if
    
    call openfile(io_unit_tfld,tecfile,"unknown","unformatted","stream")
    call deletefile(io_unit_tfld,tecfile)
    
    call openfile(io_unit_tfld,tecfile,"unknown","unformatted","stream")
    call tecout_bin_flowfield(io_unit_tfld)
    call closefile(io_unit_tfld,tecfile)

end subroutine output_flowfield_unsteady

subroutine output_monitor_unsteady(n)
    use mod_kndconsts, only : kind_int,len_char_file
    use mod_constants, only : io_unit_tfld
    use mod_variables, only : nstep,ns_mean,nprms
    implicit none
    integer(kind_int), intent(in)  :: n
    integer(kind_int)        :: m
    character(len_char_file) :: tecfile
    character(len=64 )       :: str1,str2,monfile
    
    monfile="monitor.plt"
    write(str1,'(i6)') nstep
    write(str2,'(i6)') n
    m = index(monfile,'.',back=.TRUE.)    
    
    tecfile = 'monitor/'//trim(adjustl(monfile(1:m-1)))//trim(adjustl(str2)) //"_" //trim(adjustl(str1)) // ".plt"
    
    call openfile(io_unit_tfld,tecfile,"unknown","unformatted","stream")
    call deletefile(io_unit_tfld,tecfile)
    
    call openfile(io_unit_tfld,tecfile,"unknown","unformatted","stream")
    call tecout_bin_slice(io_unit_tfld,n)
    call closefile(io_unit_tfld,tecfile)

end subroutine output_monitor_unsteady

subroutine output_vortfield_unsteady
    use mod_kndconsts, only : kind_int,len_char_file
    use mod_constants, only : io_unit_tfld
    use mod_variables, only : pltfile,nstep,ns_prms,nprms
    implicit none
    integer(kind_int)        :: m    
    character(len_char_file) :: tecfile
    character(len=64 )       :: str1
    
    write(str1,'(i6)') nstep    
    m = index(pltfile,'.',back=.TRUE.)    
    if (m == 0) then
        if(nprms>0 .and. nstep>ns_prms) then
            tecfile = trim(adjustl(pltfile)) //trim(adjustl(str1)) // "_vort.plt"
        else
            tecfile = trim(adjustl(pltfile)) // "_vort.plt"
        end if
    else
        if(nprms>0 .and. nstep>ns_prms) then
            tecfile = trim(adjustl(pltfile(1:m-1))) //trim(adjustl(str1)) // "_vort.plt"
        else
            tecfile = trim(adjustl(pltfile(1:m-1))) // "_vort.plt"
        endif
    end if
    
    call openfile(io_unit_tfld,tecfile,"unknown","unformatted","stream")
    call deletefile(io_unit_tfld,tecfile)
    
    call openfile(io_unit_tfld,tecfile,"unknown","unformatted","stream")
    call tecout_bin_vortfield(io_unit_tfld)
    call closefile(io_unit_tfld,tecfile)

end subroutine output_vortfield_unsteady

subroutine output_flowfield_mean
    use mod_kndconsts, only : kind_int,len_char_file
    use mod_constants, only : io_unit_tfld
    use mod_variables, only : pltfile
    implicit none
    integer(kind_int)        :: m    
    character(len_char_file) :: tecfile
    
    m = index(pltfile,'.',back=.TRUE.)
    if (m == 0) then
       tecfile = trim(adjustl(pltfile)) // "_flow_mean.plt"
    else
       tecfile = trim(adjustl(pltfile(1:m-1))) // "_flow_mean.plt"
    end if
    
    call openfile(io_unit_tfld,tecfile,"unknown","unformatted","stream")
    call deletefile(io_unit_tfld,tecfile)
    
    call openfile(io_unit_tfld,tecfile,"unknown","unformatted","stream")
    call tecout_bin_flowfield(io_unit_tfld)
    call closefile(io_unit_tfld,tecfile)

end subroutine output_flowfield_mean

subroutine output_flowfield_frms
    use mod_kndconsts, only : kind_int,len_char_file
    use mod_constants, only : io_unit_tfld
    use mod_variables, only : pltfile
    implicit none
    integer(kind_int)        :: m    
    character(len_char_file) :: tecfile
    
    m = index(pltfile,'.',back=.TRUE.)
    if (m == 0) then
       tecfile = trim(adjustl(pltfile)) // "_flow_frms.plt"
    else
       tecfile = trim(adjustl(pltfile(1:m-1))) // "_flow_frms.plt"
    end if
    
    call openfile(io_unit_tfld,tecfile,"unknown","unformatted","stream")
    call deletefile(io_unit_tfld,tecfile)
    
    call openfile(io_unit_tfld,tecfile,"unknown","unformatted","stream")
    call tecout_bin_flowfield_frms(io_unit_tfld)
    call closefile(io_unit_tfld,tecfile)

end subroutine output_flowfield_frms

subroutine tecout_bin_flowfield(io_unit)
    use mod_kndconsts, only : kind_int,kind_real,kind_single
    use mod_datatypes, only : fld_array_t,top_block_t
    use mod_fieldvars, only : nblocks,mb_top,mb_xyz,mb_pv,mb_dpv
    use mod_constants, only : nvis_euler
    use mod_variables, only : nvis
    use mod_parallels
    use mod_tecplotio
    implicit none
    integer(kind_int), intent(in) :: io_unit
    integer(kind_int)          :: nb,i,j,k,m,ierr
    integer(kind_int)          :: ni,nj,nk
    type(top_block_t), pointer :: top
    type(fld_array_t), pointer :: xyz(:),pv(:),dpv(:)
    real(kind_real),   pointer :: vars(:,:,:,:)
    integer(kind_int)          :: ntecvars,ndata
    character(len=256)         :: tecvarnames
    character(len=128)         :: zonename
    character(len=32 )         :: str1
#ifdef PARALLEL
    integer(kind_int) :: pid,packsize,tag
    integer(kind_int) :: status(MPI_STATUS_SIZE)
#endif

    if (nvis > nvis_euler) then
        ntecvars = 10
        tecvarnames = "X,Y,Z,Ro,U,V,W,P,Q,DV"
    else
        ntecvars = 8
        tecvarnames = "X,Y,Z,Ro,U,V,W,P"        
    end if

    if (kind_real == kind_single) then
        ndata = 1
    else
        ndata = 2
    end if

#ifdef PARALLEL
    if (myid == master) then
#endif
        call tecio_ini(io_unit,"HOSTA",ntecvars,tecvarnames)

        do nb=1,nblocks
            top => mb_top(nb)

            ni = top%nijk(1)
            nj = top%nijk(2)
            nk = top%nijk(3)

            write(str1,'(i6)') nb

            zonename = 'BLK'//trim(adjustl(str1))

            call tecio_zone(io_unit,trim(zonename),0,ni,nj,nk)
        end do

        call tecio_eohmark(io_unit)
#ifdef PARALLEL
    end if
#endif

    do nb=1,nblocks
        top => mb_top(nb)
        
        ni = top%nijk(1)
        nj = top%nijk(2)
        nk = top%nijk(3)        

        xyz  => mb_xyz(nb)%fld
        pv   => mb_pv(nb)%fld
        if (nvis > nvis_euler) then
            dpv  => mb_dpv(nb)%fld
        end if

#ifdef PARALLEL
        packsize = ni*nj*nk*ntecvars

        tag = nb

        pid = mb_top(nb)%pid
        if ((pid /= master .and. myid == pid) .or. myid == master) then
#endif
            allocate(vars(1:ni,1:nj,1:nk,ntecvars), stat=ierr)
#ifdef PARALLEL
        end if

        if (myid == pid) then
#endif
            do k=1,nk
            do j=1,nj
            do i=1,ni
                vars(i,j,k, 1) = xyz(1)%r3d(i,j,k)
                vars(i,j,k, 2) = xyz(2)%r3d(i,j,k)
                vars(i,j,k, 3) = xyz(3)%r3d(i,j,k)
                vars(i,j,k, 4) = pv(1)%r3d(i,j,k)
                vars(i,j,k, 5) = pv(2)%r3d(i,j,k)
                vars(i,j,k, 6) = pv(3)%r3d(i,j,k)
                vars(i,j,k, 7) = pv(4)%r3d(i,j,k)
                vars(i,j,k, 8) = pv(5)%r3d(i,j,k)
                if (nvis > nvis_euler) then
                    vars(i,j,k, 9) = dpv(5)%r3d(i,j,k)*dpv(9)%r3d(i,j,k) - dpv(6)%r3d(i,j,k)*dpv(8)%r3d(i,j,k) + &
                                     dpv(1)%r3d(i,j,k)*dpv(9)%r3d(i,j,k) - dpv(3)%r3d(i,j,k)*dpv(7)%r3d(i,j,k) + &
                                     dpv(1)%r3d(i,j,k)*dpv(5)%r3d(i,j,k) - dpv(2)%r3d(i,j,k)*dpv(4)%r3d(i,j,k)
                    vars(i,j,k,10) = dpv(1)%r3d(i,j,k) + dpv(5)%r3d(i,j,k) + dpv(9)%r3d(i,j,k)
                end if
            end do
            end do
            end do
#ifdef PARALLEL
        end if

        if (pid /= master) then
            if (myid == pid) then
                call MPI_SEND(vars,packsize,kind_real_mpi, &
                              master,tag,MPI_COMM_WORLD,ierr)
            end if

            if (myid == master) then
                call MPI_RECV(vars,packsize,kind_real_mpi, &
                              pid,tag,MPI_COMM_WORLD,status,ierr)
            end if
        end if

        if (myid == master) then
#endif
            call tecio_data(io_unit,ndata)

            write(io_unit)((((vars(i,j,k,m),i=1,ni), &
                                            j=1,nj), &
                                            k=1,nk), &
                                            m=1,ntecvars)
#ifdef PARALLEL
        end if

        if ((pid /= master .and. myid == pid) .or. myid == master) then
#endif
            deallocate(vars, stat=ierr)
#ifdef PARALLEL
        end if
#endif
    end do

end subroutine tecout_bin_flowfield

subroutine tecout_bin_slice(io_unit,n)
    use mod_kndconsts, only : kind_int,kind_real,kind_single
    use mod_datatypes, only : fld_array_t,top_block_t
    use mod_fieldvars, only : mb_top,mb_xyz,mb_pv
    use mod_constants, only : nvis_euler
    use mod_variables, only : nbmoni,nimoni,njmoni,nkmoni
    use mod_variables, only : nvis
    use mod_parallels
    use mod_tecplotio
    implicit none
    integer(kind_int), intent(in) :: io_unit,n
    integer(kind_int)          :: nb,i,j,k,m,ierr
    integer(kind_int)          :: ni,nj,nk,st(3),ed(3)
    type(top_block_t), pointer :: top
    type(fld_array_t), pointer :: xyz(:),pv(:)
    real(kind_real),   pointer :: vars(:,:,:,:)
    integer(kind_int)          :: ntecvars,ndata
    character(len=256)         :: tecvarnames
    character(len=128)         :: zonename
    character(len=32 )         :: str1
#ifdef PARALLEL
    integer(kind_int) :: pid,packsize,tag
    integer(kind_int) :: status(MPI_STATUS_SIZE)
#endif

    ntecvars = 8
    tecvarnames = "X,Y,Z,Ro,U,V,W,P"

    if (kind_real == kind_single) then
        ndata = 1
    else
        ndata = 2
    end if

#ifdef PARALLEL
    if (myid == master) then
#endif
        call tecio_ini(io_unit,"HOSTA",ntecvars,tecvarnames)

        nb=nbmoni(n)
        
        top => mb_top(nb)

        if(nimoni(n)>0) then
            ni    = 1
            st(1) = nimoni(n)
            ed(1) = nimoni(n)
        else
            ni = top%nijk(1)
            st(1) = 1
            ed(1) = top%nijk(1)
        end if
        
        if(njmoni(n)>0) then
            nj    = 1
            st(2) = njmoni(n)
            ed(2) = njmoni(n)
        else
            nj = top%nijk(2)
            st(2) = 1
            ed(2) = top%nijk(2)
        end if
        
        if(nkmoni(n)>0) then
            nk    = 1
            st(3) = nkmoni(n)
            ed(3) = nkmoni(n)
        else
            nk = top%nijk(3)
            st(3) = 1
            ed(3) = top%nijk(3)
        end if        

        write(str1,'(i6)') nb

        zonename = 'BLK'//trim(adjustl(str1))

        call tecio_zone(io_unit,trim(zonename),0,ni,nj,nk)

        call tecio_eohmark(io_unit)
#ifdef PARALLEL
    end if
#endif

    nb=nbmoni(n)
    
    top => mb_top(nb)

    if(nimoni(n)>0) then
        ni    = 1
        st(1) = nimoni(n)
        ed(1) = nimoni(n)
    else
        ni = top%nijk(1)
        st(1) = 1
        ed(1) = top%nijk(1)
    end if
    
    if(njmoni(n)>0) then
        nj    = 1
        st(2) = njmoni(n)
        ed(2) = njmoni(n)
    else
        nj = top%nijk(2)
        st(2) = 1
        ed(2) = top%nijk(2)
    end if
    
    if(nkmoni(n)>0) then
        nk    = 1
        st(3) = nkmoni(n)
        ed(3) = nkmoni(n)
    else
        nk = top%nijk(3)
        st(3) = 1
        ed(3) = top%nijk(3)
    end if        

    pv   => mb_pv(nb)%fld
    xyz  => mb_xyz(nb)%fld

#ifdef PARALLEL
    packsize = ni*nj*nk*ntecvars

    tag = nb

    pid = mb_top(nb)%pid
    if ((pid /= master .and. myid == pid) .or. myid == master) then
#endif
        allocate(vars(st(1):ed(1),st(2):ed(2),st(3):ed(3),ntecvars), stat=ierr)
#ifdef PARALLEL
    end if

    if (myid == pid) then
#endif
        do k=st(3),ed(3)
        do j=st(2),ed(2)
        do i=st(1),ed(1)
            vars(i,j,k, 1) = xyz(1)%r3d(i,j,k)
            vars(i,j,k, 2) = xyz(2)%r3d(i,j,k)
            vars(i,j,k, 3) = xyz(3)%r3d(i,j,k)
            vars(i,j,k, 4) = pv(1)%r3d(i,j,k)
            vars(i,j,k, 5) = pv(2)%r3d(i,j,k)
            vars(i,j,k, 6) = pv(3)%r3d(i,j,k)
            vars(i,j,k, 7) = pv(4)%r3d(i,j,k)
            vars(i,j,k, 8) = pv(5)%r3d(i,j,k)
        end do
        end do
        end do
#ifdef PARALLEL
    end if

    if (pid /= master) then
        if (myid == pid) then
            call MPI_SEND(vars,packsize,kind_real_mpi, &
                          master,tag,MPI_COMM_WORLD,ierr)
        end if

        if (myid == master) then
            call MPI_RECV(vars,packsize,kind_real_mpi, &
                          pid,tag,MPI_COMM_WORLD,status,ierr)
        end if
    end if

    if (myid == master) then
#endif
        call tecio_data(io_unit,ndata)

        write(io_unit)((((vars(i,j,k,m),i=st(1),ed(1)), &
                                        j=st(2),ed(2)), &
                                        k=st(3),ed(3)), &
                                        m=1,ntecvars)
#ifdef PARALLEL
    end if

    if ((pid /= master .and. myid == pid) .or. myid == master) then
#endif
        deallocate(vars, stat=ierr)
#ifdef PARALLEL
    end if
#endif

end subroutine tecout_bin_slice

subroutine tecout_bin_vortfield(io_unit)
    use mod_kndconsts, only : kind_int,kind_real,kind_single
    use mod_datatypes, only : fld_array_t,top_block_t
    use mod_fieldvars, only : nblocks,mb_top,mb_xyz,mb_pv,mb_dpv
    use mod_fieldvars, only : mb_ws,mb_dws,mb_t
    use mod_constants, only : nvis_euler
    use mod_variables, only : nvis
    use mod_parallels
    use mod_tecplotio
    implicit none
    integer(kind_int), intent(in) :: io_unit
    integer(kind_int)          :: nb,i,j,k,m,ierr
    integer(kind_int)          :: ni,nj,nk
    type(top_block_t), pointer :: top
    type(fld_array_t), pointer :: xyz(:),pv(:),dpv(:)
    type(fld_array_t), pointer :: ws(:),dws(:),t(:)
    real(kind_real),   pointer :: vars(:,:,:,:)
    integer(kind_int)          :: ntecvars,ndata
    character(len=256)         :: tecvarnames
    character(len=128)         :: zonename
    character(len=32 )         :: str1
#ifdef PARALLEL
    integer(kind_int) :: pid,packsize,tag
    integer(kind_int) :: status(MPI_STATUS_SIZE)
#endif

    if (nvis > nvis_euler) then
        ntecvars = 16
        tecvarnames = "X,Y,Z,Ro,U,V,W,P,Q,DV,LCx,LCy,LCz,LD1,LD2,LD3"
    else
        ntecvars = 8
        tecvarnames = "X,Y,Z,Ro,U,V,W,P"        
    end if
    
    call solve_Lamb

    if (kind_real == kind_single) then
        ndata = 1
    else
        ndata = 2
    end if

#ifdef PARALLEL
    if (myid == master) then
#endif
        call tecio_ini(io_unit,"HOSTA",ntecvars,tecvarnames)

        do nb=1,nblocks
            top => mb_top(nb)

            ni = top%nijk(1)
            nj = top%nijk(2)
            nk = top%nijk(3)

            write(str1,'(i6)') nb

            zonename = 'BLK'//trim(adjustl(str1))

            call tecio_zone(io_unit,trim(zonename),0,ni,nj,nk)
        end do

        call tecio_eohmark(io_unit)
#ifdef PARALLEL
    end if
#endif

    do nb=1,nblocks
        top => mb_top(nb)
        
        ni = top%nijk(1)
        nj = top%nijk(2)
        nk = top%nijk(3)        

        xyz  => mb_xyz(nb)%fld
        pv   => mb_pv(nb)%fld
        if (nvis > nvis_euler) then
            dpv  => mb_dpv(nb)%fld
            ws   => mb_ws(nb)%fld
            dws  => mb_dws(nb)%fld
            t    => mb_t(nb)%fld
        end if

#ifdef PARALLEL
        packsize = ni*nj*nk*ntecvars

        tag = nb

        pid = mb_top(nb)%pid
        if ((pid /= master .and. myid == pid) .or. myid == master) then
#endif
            allocate(vars(1:ni,1:nj,1:nk,ntecvars), stat=ierr)
#ifdef PARALLEL
        end if

        if (myid == pid) then
#endif
            do k=1,nk
            do j=1,nj
            do i=1,ni
                vars(i,j,k, 1) = xyz(1)%r3d(i,j,k)
                vars(i,j,k, 2) = xyz(2)%r3d(i,j,k)
                vars(i,j,k, 3) = xyz(3)%r3d(i,j,k)
                vars(i,j,k, 4) = pv(1)%r3d(i,j,k)
                vars(i,j,k, 5) = pv(2)%r3d(i,j,k)
                vars(i,j,k, 6) = pv(3)%r3d(i,j,k)
                vars(i,j,k, 7) = pv(4)%r3d(i,j,k)
                vars(i,j,k, 8) = pv(5)%r3d(i,j,k)
                if (nvis > nvis_euler) then
                    vars(i,j,k, 9) = dpv(5)%r3d(i,j,k)*dpv(9)%r3d(i,j,k) - dpv(6)%r3d(i,j,k)*dpv(8)%r3d(i,j,k) + &
                                     dpv(1)%r3d(i,j,k)*dpv(9)%r3d(i,j,k) - dpv(3)%r3d(i,j,k)*dpv(7)%r3d(i,j,k) + &
                                     dpv(1)%r3d(i,j,k)*dpv(5)%r3d(i,j,k) - dpv(2)%r3d(i,j,k)*dpv(4)%r3d(i,j,k)
                    vars(i,j,k,10) = dpv(1)%r3d(i,j,k) + dpv(5)%r3d(i,j,k) + dpv(9)%r3d(i,j,k)
                    vars(i,j,k,11) = dws(11)%r3d(i,j,k) - dws(9)%r3d(i,j,k) ! ws(5)%r3d(i,j,k) !
                    vars(i,j,k,12) = dws(6)%r3d(i,j,k) - dws(10)%r3d(i,j,k) ! ws(6)%r3d(i,j,k) !
                    vars(i,j,k,13) = dws(7)%r3d(i,j,k) - dws(5)%r3d(i,j,k)  ! ws(7)%r3d(i,j,k) !
                    vars(i,j,k,14) = dws(13)%r3d(i,j,k) + dws(17)%r3d(i,j,k) + dws(21)%r3d(i,j,k)
                    vars(i,j,k,15) = dpv(10)%r3d(i,j,k)*dws(1)%r3d(i,j,k) + &
                                     dpv(11)%r3d(i,j,k)*dws(2)%r3d(i,j,k) + &
                                     dpv(12)%r3d(i,j,k)*dws(3)%r3d(i,j,k)
                    vars(i,j,k,16) = t(1)%r3d(i,j,k)*dws(22)%r3d(i,j,k) + &
                                     t(1)%r3d(i,j,k)*dws(26)%r3d(i,j,k) + &
                                     t(1)%r3d(i,j,k)*dws(30)%r3d(i,j,k)
                end if
            end do
            end do
            end do
#ifdef PARALLEL
        end if

        if (pid /= master) then
            if (myid == pid) then
                call MPI_SEND(vars,packsize,kind_real_mpi, &
                              master,tag,MPI_COMM_WORLD,ierr)
            end if

            if (myid == master) then
                call MPI_RECV(vars,packsize,kind_real_mpi, &
                              pid,tag,MPI_COMM_WORLD,status,ierr)
            end if
        end if

        if (myid == master) then
#endif
            call tecio_data(io_unit,ndata)

            write(io_unit)((((vars(i,j,k,m),i=1,ni), &
                                            j=1,nj), &
                                            k=1,nk), &
                                            m=1,ntecvars)
#ifdef PARALLEL
        end if

        if ((pid /= master .and. myid == pid) .or. myid == master) then
#endif
            deallocate(vars, stat=ierr)
#ifdef PARALLEL
        end if
#endif
    end do

end subroutine tecout_bin_vortfield

subroutine tecout_bin_flowfield_frms(io_unit)
    use mod_kndconsts, only : kind_int,kind_real,kind_single
    use mod_datatypes, only : fld_array_t,top_block_t
    use mod_fieldvars, only : nblocks,mb_top,mb_xyz,mb_pv,mb_dpv
    use mod_constants, only : nvis_euler
    use mod_variables, only : nvis
    use mod_parallels
    use mod_tecplotio
    implicit none
    integer(kind_int), intent(in) :: io_unit
    integer(kind_int)          :: nb,i,j,k,m,ierr
    integer(kind_int)          :: ni,nj,nk
    type(top_block_t), pointer :: top
    type(fld_array_t), pointer :: xyz(:),pv(:),dpv(:)
    real(kind_real),   pointer :: vars(:,:,:,:)
    integer(kind_int)          :: ntecvars,ndata
    character(len=256)         :: tecvarnames
    character(len=128)         :: zonename
    character(len=32 )         :: str1
#ifdef PARALLEL
    integer(kind_int) :: pid,packsize,tag
    integer(kind_int) :: status(MPI_STATUS_SIZE)
#endif

    ntecvars = 8
    tecvarnames = "X,Y,Z,Ro,U,V,W,P"

    if (kind_real == kind_single) then
        ndata = 1
    else
        ndata = 2
    end if

#ifdef PARALLEL
    if (myid == master) then
#endif
        call tecio_ini(io_unit,"HOSTA",ntecvars,tecvarnames)

        do nb=1,nblocks
            top => mb_top(nb)

            ni = top%nijk(1)
            nj = top%nijk(2)
            nk = top%nijk(3)

            write(str1,'(i6)') nb

            zonename = 'BLK'//trim(adjustl(str1))

            call tecio_zone(io_unit,trim(zonename),0,ni,nj,nk)
        end do

        call tecio_eohmark(io_unit)
#ifdef PARALLEL
    end if
#endif

    do nb=1,nblocks
        top => mb_top(nb)
        
        ni = top%nijk(1)
        nj = top%nijk(2)
        nk = top%nijk(3)        

        xyz  => mb_xyz(nb)%fld
        pv   => mb_pv(nb)%fld

#ifdef PARALLEL
        packsize = ni*nj*nk*ntecvars

        tag = nb

        pid = mb_top(nb)%pid
        if ((pid /= master .and. myid == pid) .or. myid == master) then
#endif
            allocate(vars(1:ni,1:nj,1:nk,ntecvars), stat=ierr)
#ifdef PARALLEL
        end if

        if (myid == pid) then
#endif
            do k=1,nk
            do j=1,nj
            do i=1,ni
                vars(i,j,k, 1) = xyz(1)%r3d(i,j,k)
                vars(i,j,k, 2) = xyz(2)%r3d(i,j,k)
                vars(i,j,k, 3) = xyz(3)%r3d(i,j,k)
                vars(i,j,k, 4) = pv(1)%r3d(i,j,k)
                vars(i,j,k, 5) = pv(2)%r3d(i,j,k)
                vars(i,j,k, 6) = pv(3)%r3d(i,j,k)
                vars(i,j,k, 7) = pv(4)%r3d(i,j,k)
                vars(i,j,k, 8) = pv(5)%r3d(i,j,k)
            end do
            end do
            end do
#ifdef PARALLEL
        end if

        if (pid /= master) then
            if (myid == pid) then
                call MPI_SEND(vars,packsize,kind_real_mpi, &
                              master,tag,MPI_COMM_WORLD,ierr)
            end if

            if (myid == master) then
                call MPI_RECV(vars,packsize,kind_real_mpi, &
                              pid,tag,MPI_COMM_WORLD,status,ierr)
            end if
        end if

        if (myid == master) then
#endif
            call tecio_data(io_unit,ndata)

            write(io_unit)((((vars(i,j,k,m),i=1,ni), &
                                            j=1,nj), &
                                            k=1,nk), &
                                            m=1,ntecvars)
#ifdef PARALLEL
        end if

        if ((pid /= master .and. myid == pid) .or. myid == master) then
#endif
            deallocate(vars, stat=ierr)
#ifdef PARALLEL
        end if
#endif
    end do

end subroutine tecout_bin_flowfield_frms

subroutine output_mean
    use mod_kndconsts, only : kind_int
    use mod_constants, only : io_unit_men,nplot3d_fun
    use mod_variables, only : solfile
    use mod_variables, only : nstep,ns_mean,nsolsav    
    use mod_fieldvars, only : mb_fmean,npvs
    use mod_interface, only : gather_output_mb_var
    use mod_parallels
    implicit none
    integer(kind_int) :: error,ierr
    
    call openfile_check(io_unit_men,trim(solfile)//".mean","unknown","unformatted","stream",ierr)
    
#ifdef PARALLEL
    call MPI_BCAST(ierr,1,kind_int_mpi,master,MPI_COMM_WORLD,error)
#endif    
    
    if (ierr == 0) then
        call openfile(io_unit_men,trim(solfile)//".mean"//".bak","unknown","unformatted","stream")
        call deletefile(io_unit_men,trim(solfile)//".mean"//".bak")
        
        call renamefile(io_unit_men,trim(solfile)//".mean")
    end if
    
    call openfile(io_unit_men,trim(solfile)//".mean","unknown","unformatted","stream")
    call deletefile(io_unit_men,trim(solfile)//".mean")
    
    call openfile(io_unit_men,trim(solfile)//".mean","unknown","unformatted","stream")
    call gather_output_mb_var(io_unit_men,mb_fmean,1,npvs,nplot3d_fun)
    call closefile(io_unit_men,trim(solfile)//".mean") 

end subroutine output_mean

subroutine output_rms
    use mod_kndconsts, only : kind_int,len_char_file
    use mod_constants, only : io_unit_rms,nplot3d_fun
    use mod_variables, only : solfile
    use mod_variables, only : nstep,nsolsav,ns_prms     
    use mod_fieldvars, only : mb_frms,npvs
    use mod_interface, only : gather_output_mb_var
    use mod_parallels
    implicit none
    integer(kind_int) :: error,ierr
    
    call openfile_check(io_unit_rms,trim(solfile)//".rms","unknown","unformatted","stream",ierr)
    
#ifdef PARALLEL
    call MPI_BCAST(ierr,1,kind_int_mpi,master,MPI_COMM_WORLD,error)
#endif

    if (ierr == 0) then
        call openfile(io_unit_rms,trim(solfile)//".rms"//".bak","unknown","unformatted","stream")
        call deletefile(io_unit_rms,trim(solfile)//".rms"//".bak")
        
        call renamefile(io_unit_rms,trim(solfile)//".rms")
    end if
    
    call openfile(io_unit_rms,trim(solfile)//".rms","unknown","unformatted","stream")
    call deletefile(io_unit_rms,trim(solfile)//".rms")
    
    call openfile(io_unit_rms,trim(solfile)//".rms","unknown","unformatted","stream")
    call gather_output_mb_var(io_unit_rms,mb_frms,1,npvs,nplot3d_fun)
    call closefile(io_unit_rms,trim(solfile)//".rms")

end subroutine output_rms

subroutine output_fmean
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : zero,nvis_euler
    use mod_constants, only : nlhs_prsgs_ust_mat,nlhs_prsgs_ust_sca,nlhs_prsgs_ust_sca_scmp
    use mod_datatypes, only : var_block_t    
    use mod_variables, only : nghnode,nvis,nlhs
    use mod_fieldvars, only : nblkcoms,blkcoms,npvs
    use mod_fieldvars, only : mb_fmean,mb_frms,mb_pv
    use mod_fieldvars, only : neqn,mb_qc
    use mod_fieldvars, only : mb_vst,mb_q0
    use mod_interface, only : assign_mb_var_uniform
    use mod_interface, only : calc_mb_var_via_sub    
    use mod_openmp
    implicit none
    integer(kind_int) :: nc,nb,i,j,k,m
    integer(kind_int) :: st(3),ed(3)
    external        :: prim2con,v1_eq_v2    

    do nc=1,nblkcoms
        nb = blkcoms(nc)%nb

        st(:) = 1 - nghnode
        ed(:) = blkcoms(nc)%top%nijk(:) + nghnode
        do m=1,npvs
!$OMP parallel do
            do k=st(3),ed(3)
            do j=st(2),ed(2)
            do i=st(1),ed(1)
                mb_frms(nb)%fld(m)%r3d(i,j,k) = mb_pv(nb)%fld(m)%r3d(i,j,k)
            end do
            end do
            end do
!$OMP end parallel do
        end do
    end do
    
    do nc=1,nblkcoms
        nb = blkcoms(nc)%nb

        st(:) = 1
        ed(:) = blkcoms(nc)%top%nijk(:)
        do m=1,npvs
!$OMP parallel do
            do k=st(3),ed(3)
            do j=st(2),ed(2)
            do i=st(1),ed(1)
                mb_pv(nb)%fld(m)%r3d(i,j,k) = mb_fmean(nb)%fld(m)%r3d(i,j,k)
            end do
            end do
            end do
!$OMP end parallel do
        end do
    end do    
    
    call prepare_linear_system
    call output_surface_mean
    call output_bcwall_mean
    call output_flowfield_mean
    
    do nc=1,nblkcoms
        nb = blkcoms(nc)%nb

        st(:) = 1 - nghnode
        ed(:) = blkcoms(nc)%top%nijk(:) + nghnode
        do m=1,npvs
!$OMP parallel do
            do k=st(3),ed(3)
            do j=st(2),ed(2)
            do i=st(1),ed(1)
                mb_pv(nb)%fld(m)%r3d(i,j,k) = mb_frms(nb)%fld(m)%r3d(i,j,k)
            end do
            end do
            end do
!$OMP end parallel do
        end do
    end do
    
    call calc_mb_var_via_sub(mb_pv,1,npvs,prim2con,mb_qc,1,neqn,nghnode)

    select case(nlhs)
    case(nlhs_prsgs_ust_mat,nlhs_prsgs_ust_sca,nlhs_prsgs_ust_sca_scmp)
        call calc_mb_var_via_sub(mb_qc,1,neqn,v1_eq_v2,mb_q0,1,neqn,nghnode)
    end select

    if (nvis > nvis_euler) then
        call assign_mb_var_uniform(mb_vst,1,1,nghnode,(/zero/))
    end if    

end subroutine output_fmean

subroutine output_frms
    use mod_kndconsts, only : kind_int,kind_real
    use mod_datatypes, only : var_block_t    
    use mod_variables, only : nghnode
    use mod_fieldvars, only : nblkcoms,blkcoms,npvs
    use mod_fieldvars, only : mb_fmean,mb_frms,mb_pv,npvs    
    use mod_openmp
    implicit none
    integer(kind_int) :: nc,nb,i,j,k,m
    integer(kind_int) :: st(3),ed(3)
    
    do nc=1,nblkcoms
        nb = blkcoms(nc)%nb

        st(:) = 1 - nghnode
        ed(:) = blkcoms(nc)%top%nijk(:) + nghnode
        do m=1,npvs
!$OMP parallel do
            do k=st(3),ed(3)
            do j=st(2),ed(2)
            do i=st(1),ed(1)
                mb_fmean(nb)%fld(m)%r3d(i,j,k) = mb_pv(nb)%fld(m)%r3d(i,j,k)
            end do
            end do
            end do
!$OMP end parallel do
        end do
    end do    

    do nc=1,nblkcoms
        nb = blkcoms(nc)%nb

        st(:) = 1 - nghnode
        ed(:) = blkcoms(nc)%top%nijk(:) + nghnode
        do m=1,npvs
!$OMP parallel do
            do k=st(3),ed(3)
            do j=st(2),ed(2)
            do i=st(1),ed(1)
                mb_pv(nb)%fld(m)%r3d(i,j,k) = mb_frms(nb)%fld(m)%r3d(i,j,k)
            end do
            end do
            end do
!$OMP end parallel do
        end do
    end do   
    
    call output_flowfield_frms
    call output_surface_rms
    
    do nc=1,nblkcoms
        nb = blkcoms(nc)%nb

        st(:) = 1 - nghnode
        ed(:) = blkcoms(nc)%top%nijk(:) + nghnode
        do m=1,npvs
!$OMP parallel do
            do k=st(3),ed(3)
            do j=st(2),ed(2)
            do i=st(1),ed(1)
                mb_pv(nb)%fld(m)%r3d(i,j,k) = mb_fmean(nb)%fld(m)%r3d(i,j,k)
            end do
            end do
            end do
!$OMP end parallel do
        end do
    end do    

end subroutine output_frms

subroutine solve_Lamb
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : one,two
    use mod_constants, only : nsgl_buffer_pvs,nsgl_buffer_dpv
    use mod_constants, only : nsgl_aver_art,nbc_inter_buf_pvs,nbc_inter_buf_dpv
    use mod_datatypes, only : fld_array_t,top_block_t
    use mod_variables, only : gamma,poo,nghnode
    use mod_fieldvars, only : npvs
    use mod_fieldvars, only : nblkcoms,blkcoms
    use mod_fieldvars, only : mb_pv,mb_dpv,mb_t
    use mod_fieldvars, only : mb_ws,mb_dws
    use mod_interface, only : pre_exchange_bc_var,post_exchange_bc_var,average_bc_var
    use mod_interface, only : patched_ghost_points,exchange_singulars,average_singulars
    implicit none
    integer(kind_int)          :: nc,nb,i,j,k,ierr
    integer(kind_int)          :: m,ni,nj,nk
    real(kind_real)            :: gama,ae,fp
    type(top_block_t), pointer :: top
    type(fld_array_t), pointer :: dpv(:),pv(:),t(:)
    type(fld_array_t), pointer :: ws(:),dws(:)

    gama = gamma
    ae = gama - one
    
    do nc=1,nblkcoms
        nb  =  blkcoms(nc)%nb
        top => blkcoms(nc)%top

        ni = top%nijk(1)
        nj = top%nijk(2)
        nk = top%nijk(3)

        t     => mb_t(nb)%fld
        pv    => mb_pv(nb)%fld
        dpv   => mb_dpv(nb)%fld
        ws    => mb_ws(nb)%fld

        do k=1-nghnode,nk+nghnode
        do j=1-nghnode,nj+nghnode
        do i=1-nghnode,ni+nghnode
            fp = two*(pv(5)%r3d(i,j,k) - poo)
            ws(1)%r3d(i,j,k) = dpv(8)%r3d(i,j,k) - dpv(6)%r3d(i,j,k)
            ws(2)%r3d(i,j,k) = dpv(3)%r3d(i,j,k) - dpv(7)%r3d(i,j,k)
            ws(3)%r3d(i,j,k) = dpv(4)%r3d(i,j,k) - dpv(2)%r3d(i,j,k)
            ws(4)%r3d(i,j,k) = fp*log(t(1)%r3d(i,j,k)/pv(5)%r3d(i,j,k)**(ae/gama))
        end do
        end do
        end do
    end do
    
    call pre_exchange_bc_var(mb_ws,1,4,nghnode,nbc_inter_buf_dpv,nsgl_aver_art)
    call exchange_singulars(mb_ws,1,5,nsgl_buffer_pvs,nsgl_aver_art)
    call post_exchange_bc_var(mb_ws,1,4,nghnode,nbc_inter_buf_dpv,nsgl_aver_art)
    call patched_ghost_points(mb_ws,1,4,nghnode)
    call average_bc_var(mb_ws,1,4,nbc_inter_buf_dpv,nsgl_aver_art)
    call average_singulars(mb_ws,1,5,nsgl_buffer_pvs,nsgl_aver_art)        
    
    call s_gradient
    
    call pre_exchange_bc_var(mb_dws,1,3,nghnode,nbc_inter_buf_dpv,nsgl_aver_art)
    call exchange_singulars(mb_dws,1,5,nsgl_buffer_pvs,nsgl_aver_art)
    call post_exchange_bc_var(mb_dws,1,3,nghnode,nbc_inter_buf_dpv,nsgl_aver_art)
    call patched_ghost_points(mb_dws,1,3,nghnode)
    call average_bc_var(mb_dws,1,3,nbc_inter_buf_dpv,nsgl_aver_art)
    call average_singulars(mb_dws,1,5,nsgl_buffer_pvs,nsgl_aver_art)
    
    do nc=1,nblkcoms
        nb  =  blkcoms(nc)%nb
        top => blkcoms(nc)%top

        ni = top%nijk(1)
        nj = top%nijk(2)
        nk = top%nijk(3)

        t     => mb_t(nb)%fld
        pv    => mb_pv(nb)%fld
        ws    => mb_ws(nb)%fld
        dws   => mb_dws(nb)%fld        

        do k=1-nghnode,nk+nghnode
        do j=1-nghnode,nj+nghnode
        do i=1-nghnode,ni+nghnode
            ws(5)%r3d(i,j,k) = ws(2)%r3d(i,j,k)*pv(4)%r3d(i,j,k) - ws(3)%r3d(i,j,k)*pv(3)%r3d(i,j,k) - &
                               t(1)%r3d(i,j,k)*dws(1)%r3d(i,j,k)
            ws(6)%r3d(i,j,k) = ws(3)%r3d(i,j,k)*pv(2)%r3d(i,j,k) - ws(1)%r3d(i,j,k)*pv(4)%r3d(i,j,k) - &
                               t(1)%r3d(i,j,k)*dws(2)%r3d(i,j,k)
            ws(7)%r3d(i,j,k) = ws(1)%r3d(i,j,k)*pv(3)%r3d(i,j,k) - ws(2)%r3d(i,j,k)*pv(2)%r3d(i,j,k) - &
                               t(1)%r3d(i,j,k)*dws(3)%r3d(i,j,k)
        end do
        end do
        end do
    end do
    
    call pre_exchange_bc_var(mb_ws,1,7,nghnode,nbc_inter_buf_dpv,nsgl_aver_art)
    call exchange_singulars(mb_ws,3,7,nsgl_buffer_pvs,nsgl_aver_art)
    call post_exchange_bc_var(mb_ws,1,7,nghnode,nbc_inter_buf_dpv,nsgl_aver_art)
    call patched_ghost_points(mb_ws,1,7,nghnode)
    call average_bc_var(mb_ws,1,7,nbc_inter_buf_dpv,nsgl_aver_art)
    call average_singulars(mb_ws,3,7,nsgl_buffer_pvs,nsgl_aver_art)
    
    do nc=1,nblkcoms
        nb  =  blkcoms(nc)%nb
        top => blkcoms(nc)%top

        ni = top%nijk(1)
        nj = top%nijk(2)
        nk = top%nijk(3)

        pv    => mb_pv(nb)%fld
        ws    => mb_ws(nb)%fld      

        do k=1-nghnode,nk+nghnode
        do j=1-nghnode,nj+nghnode
        do i=1-nghnode,ni+nghnode
            ws(8)%r3d(i,j,k)  = ws(3)%r3d(i,j,k)*pv(3)%r3d(i,j,k) - ws(2)%r3d(i,j,k)*pv(4)%r3d(i,j,k)
            ws(9)%r3d(i,j,k)  = ws(1)%r3d(i,j,k)*pv(4)%r3d(i,j,k) - ws(3)%r3d(i,j,k)*pv(2)%r3d(i,j,k)
            ws(10)%r3d(i,j,k) = ws(2)%r3d(i,j,k)*pv(2)%r3d(i,j,k) - ws(1)%r3d(i,j,k)*pv(3)%r3d(i,j,k)
        end do
        end do
        end do
    end do
    
    call pre_exchange_bc_var(mb_ws,1,10,nghnode,nbc_inter_buf_dpv,nsgl_aver_art)
    call exchange_singulars(mb_ws,6,10,nsgl_buffer_pvs,nsgl_aver_art)
    call post_exchange_bc_var(mb_ws,1,10,nghnode,nbc_inter_buf_dpv,nsgl_aver_art)
    call patched_ghost_points(mb_ws,1,10,nghnode)
    call average_bc_var(mb_ws,1,10,nbc_inter_buf_dpv,nsgl_aver_art)
    call average_singulars(mb_ws,6,10,nsgl_buffer_pvs,nsgl_aver_art)
    
    call l_gradient

end subroutine solve_Lamb

subroutine s_gradient
    use mod_kndconsts, only : kind_int
    use mod_constants, only : nintplt_node2,nintplt_node4
    use mod_constants, only : nintplt_node6,nintplt_node6w,nintplt_node8
    use mod_constants, only : nintplt_scsl4,nintplt_scsl6
    use mod_constants, only : nintplt_node6e,nintplt_node8e
    use mod_constants, only : nintplt_scsl4e,nintplt_scsl6e
    use mod_variables, only : nd3int_vis
    implicit none
    external :: ve_via_node2,ve_via_node4
    external :: ve_via_node6,ve_via_node6w,ve_via_node8
    external :: ve_via_scsl4,ve_via_scsl6
    external :: ve_via_node6e,ve_via_node8e
    external :: ve_via_scsl4e,ve_via_scsl6e
    
    select case(nd3int_vis)
    case(nintplt_node2)
        call grad_s_dnvis(ve_via_node2)
    case(nintplt_node4)
        call grad_s_dnvis(ve_via_node4)
    case(nintplt_node6)
        call grad_s_dnvis(ve_via_node6)
    case(nintplt_node6w)
        call grad_s_dnvis(ve_via_node6w)
    case(nintplt_node6e)
        call grad_s_dnvis(ve_via_node6e)
    case(nintplt_node8)
        call grad_s_dnvis(ve_via_node8)
    case(nintplt_node8e)
        call grad_s_dnvis(ve_via_node8e)
    case(nintplt_scsl4)
        call grad_s_dnvis(ve_via_scsl4)
    case(nintplt_scsl4e)
        call grad_s_dnvis(ve_via_scsl4e)
    case(nintplt_scsl6)
        call grad_s_dnvis(ve_via_scsl6)
    case(nintplt_scsl6e)
        call grad_s_dnvis(ve_via_scsl6e)                
    case default

    end select

end subroutine s_gradient


subroutine grad_s_dnvis(sub_intvis)
    use mod_constants, only : nderive_edge6,nderive_edge4,nderive_edge2
    use mod_constants, only : nderive_ehen6,nderive_ehen4,nderive_ehen8,nderive_ehcs6
    use mod_constants, only : nderive_scsl4,nderive_scsl6
    use mod_constants, only : nderive_ehen6e,nderive_ehen8e,nderive_ehcs6e
    use mod_constants, only : nderive_scsl4e,nderive_scsl6e
    use mod_variables, only : nd3der_vis,nghnode,nghedge
    implicit none
    external :: sub_intvis
    external :: dn_via_edge6,dn_via_edge4,dn_via_edge2
    external :: dn_via_ehen6,dn_via_ehen4,dn_via_ehen8,dn_via_ehcs6
    external :: dn_via_scsl4,dn_via_scsl6
    external :: dn_via_ehen6e,dn_via_ehen8e,dn_via_ehcs6e
    external :: dn_via_scsl4e,dn_via_scsl6e

    select case(nd3der_vis)
    case(nderive_edge2)
        call calc_grad_s(nghnode,nghedge,sub_intvis,dn_via_edge2)
    case(nderive_edge4)
        call calc_grad_s(nghnode,nghedge,sub_intvis,dn_via_edge4)
    case(nderive_edge6)
        call calc_grad_s(nghnode,nghedge,sub_intvis,dn_via_edge6)
    case(nderive_ehen4)
        call calc_grad_s(nghnode,nghedge,sub_intvis,dn_via_ehen4)
    case(nderive_ehen6)
        call calc_grad_s(nghnode,nghedge,sub_intvis,dn_via_ehen6)
    case(nderive_ehen6e)
        call calc_grad_s(nghnode,nghedge,sub_intvis,dn_via_ehen6e)
    case(nderive_ehcs6)
        call calc_grad_s(nghnode,nghedge,sub_intvis,dn_via_ehcs6)
    case(nderive_ehcs6e)
        call calc_grad_s(nghnode,nghedge,sub_intvis,dn_via_ehcs6e)
    case(nderive_ehen8)
        call calc_grad_s(nghnode,nghedge,sub_intvis,dn_via_ehen8)
    case(nderive_ehen8e)
        call calc_grad_s(nghnode,nghedge,sub_intvis,dn_via_ehen8e)
    case(nderive_scsl4)
        call calc_grad_s(nghnode,nghedge,sub_intvis,dn_via_scsl4)
    case(nderive_scsl4e)
        call calc_grad_s(nghnode,nghedge,sub_intvis,dn_via_scsl4e)
    case(nderive_scsl6)
        call calc_grad_s(nghnode,nghedge,sub_intvis,dn_via_scsl6)
    case(nderive_scsl6e)
        call calc_grad_s(nghnode,nghedge,sub_intvis,dn_via_scsl6e)                
    case default

    end select

end subroutine grad_s_dnvis


subroutine calc_grad_s(ngn,nge,sub_intvis,sub_dnvis)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : one
    use mod_constants, only : nfsf_vis_d3int,nfsf_vis_d3der
    use mod_datatypes, only : fld_array_t,var_block_t,top_block_t
    use mod_fieldvars, only : nblkcoms,blkcoms
    use mod_fieldvars, only : mb_ws,mb_sxyz,mb_vol,mb_dws
    use mod_interface, only : mb_var_pointer_create,mb_var_pointer_assign
    use mod_interface, only : mb_var_pointer_delete,calc_mb_dn_via_node3
    use mod_openmp
    implicit none
    integer(kind_int), intent(in) :: ngn,nge
    external                      :: sub_intvis,sub_dnvis
    integer(kind_int)          :: nc,nb,i,j,k,m,m1,m2,m3
    integer(kind_int)          :: st(3),ed(3)
    real(kind_real)            :: kx,ky,kz,ex,ey,ez,cx,cy,cz
    real(kind_real)            :: ovol,der(12)
    type(top_block_t), pointer :: top
    type(fld_array_t), pointer :: sxyz(:),vol(:),dws(:)
    type(var_block_t), pointer :: mb_uvwt(:)

    call mb_var_pointer_create(mb_uvwt,1,1)
    call mb_var_pointer_assign(mb_uvwt,1,1,mb_ws,4,4)
    do m=1,1
        m1 = 3*m - 2
        m2 = m1 + 1
        m3 = m2 + 1

        call calc_mb_dn_via_node3(mb_uvwt,m,ngn,sub_intvis,nge,sub_dnvis,nfsf_vis_d3int,nfsf_vis_d3der,mb_dws,m1,1)
        call calc_mb_dn_via_node3(mb_uvwt,m,ngn,sub_intvis,nge,sub_dnvis,nfsf_vis_d3int,nfsf_vis_d3der,mb_dws,m2,2)
        call calc_mb_dn_via_node3(mb_uvwt,m,ngn,sub_intvis,nge,sub_dnvis,nfsf_vis_d3int,nfsf_vis_d3der,mb_dws,m3,3)
    end do
    call mb_var_pointer_delete(mb_uvwt)

    do nc=1,nblkcoms
        nb  =  blkcoms(nc)%nb
        top => blkcoms(nc)%top

        sxyz => mb_sxyz(nb)%fld
        vol  => mb_vol(nb)%fld
        dws  => mb_dws(nb)%fld

        st(:) = 1
        ed(:) = top%nijk(:)
!$OMP parallel do private(k,j,i,kx,ky,kz,ex,ey,ez,cx,cy,cz,ovol,m1,m2,m3,der)
        do k=st(3),ed(3)
        do j=st(2),ed(2)
        do i=st(1),ed(1)
            kx = sxyz(1)%r3d(i,j,k)
            ky = sxyz(2)%r3d(i,j,k)
            kz = sxyz(3)%r3d(i,j,k)

            ex = sxyz(4)%r3d(i,j,k)
            ey = sxyz(5)%r3d(i,j,k)
            ez = sxyz(6)%r3d(i,j,k)

            cx = sxyz(7)%r3d(i,j,k)
            cy = sxyz(8)%r3d(i,j,k)
            cz = sxyz(9)%r3d(i,j,k)

            ovol = one/vol(1)%r3d(i,j,k)

            do m=1,3
               der(m) = dws(m)%r3d(i,j,k)
            end do

            do m=1,1
                m1 = 3*m - 2
                m2 = m1 + 1
                m3 = m2 + 1

                dws(m1)%r3d(i,j,k) = ( kx*der(m1) + ex*der(m2) + cx*der(m3) ) * ovol
                dws(m2)%r3d(i,j,k) = ( ky*der(m1) + ey*der(m2) + cy*der(m3) ) * ovol
                dws(m3)%r3d(i,j,k) = ( kz*der(m1) + ez*der(m2) + cz*der(m3) ) * ovol
            end do
        end do
        end do
        end do
!$OMP end parallel do
    end do


end subroutine calc_grad_s

subroutine l_gradient
    use mod_kndconsts, only : kind_int
    use mod_constants, only : nintplt_node2,nintplt_node4
    use mod_constants, only : nintplt_node6,nintplt_node6w,nintplt_node8
    use mod_constants, only : nintplt_scsl4,nintplt_scsl6
    use mod_constants, only : nintplt_node6e,nintplt_node8e
    use mod_constants, only : nintplt_scsl4e,nintplt_scsl6e
    use mod_variables, only : nd3int_vis
    implicit none
    external :: ve_via_node2,ve_via_node4
    external :: ve_via_node6,ve_via_node6w,ve_via_node8
    external :: ve_via_scsl4,ve_via_scsl6
    external :: ve_via_node6e,ve_via_node8e
    external :: ve_via_scsl4e,ve_via_scsl6e
    
    select case(nd3int_vis)
    case(nintplt_node2)
        call grad_l_dnvis(ve_via_node2)
    case(nintplt_node4)
        call grad_l_dnvis(ve_via_node4)
    case(nintplt_node6)
        call grad_l_dnvis(ve_via_node6)
    case(nintplt_node6w)
        call grad_l_dnvis(ve_via_node6w)
    case(nintplt_node6e)
        call grad_l_dnvis(ve_via_node6e)
    case(nintplt_node8)
        call grad_l_dnvis(ve_via_node8)
    case(nintplt_node8e)
        call grad_l_dnvis(ve_via_node8e)
    case(nintplt_scsl4)
        call grad_l_dnvis(ve_via_scsl4)
    case(nintplt_scsl4e)
        call grad_l_dnvis(ve_via_scsl4e)
    case(nintplt_scsl6)
        call grad_l_dnvis(ve_via_scsl6)
    case(nintplt_scsl6e)
        call grad_l_dnvis(ve_via_scsl6e)                
    case default

    end select

end subroutine l_gradient


subroutine grad_l_dnvis(sub_intvis)
    use mod_constants, only : nderive_edge6,nderive_edge4,nderive_edge2
    use mod_constants, only : nderive_ehen6,nderive_ehen4,nderive_ehen8,nderive_ehcs6
    use mod_constants, only : nderive_scsl4,nderive_scsl6
    use mod_constants, only : nderive_ehen6e,nderive_ehen8e,nderive_ehcs6e
    use mod_constants, only : nderive_scsl4e,nderive_scsl6e
    use mod_variables, only : nd3der_vis,nghnode,nghedge
    implicit none
    external :: sub_intvis
    external :: dn_via_edge6,dn_via_edge4,dn_via_edge2
    external :: dn_via_ehen6,dn_via_ehen4,dn_via_ehen8,dn_via_ehcs6
    external :: dn_via_scsl4,dn_via_scsl6
    external :: dn_via_ehen6e,dn_via_ehen8e,dn_via_ehcs6e
    external :: dn_via_scsl4e,dn_via_scsl6e

    select case(nd3der_vis)
    case(nderive_edge2)
        call calc_grad_l(nghnode,nghedge,sub_intvis,dn_via_edge2)
    case(nderive_edge4)
        call calc_grad_l(nghnode,nghedge,sub_intvis,dn_via_edge4)
    case(nderive_edge6)
        call calc_grad_l(nghnode,nghedge,sub_intvis,dn_via_edge6)
    case(nderive_ehen4)
        call calc_grad_l(nghnode,nghedge,sub_intvis,dn_via_ehen4)
    case(nderive_ehen6)
        call calc_grad_l(nghnode,nghedge,sub_intvis,dn_via_ehen6)
    case(nderive_ehen6e)
        call calc_grad_l(nghnode,nghedge,sub_intvis,dn_via_ehen6e)
    case(nderive_ehcs6)
        call calc_grad_l(nghnode,nghedge,sub_intvis,dn_via_ehcs6)
    case(nderive_ehcs6e)
        call calc_grad_l(nghnode,nghedge,sub_intvis,dn_via_ehcs6e)
    case(nderive_ehen8)
        call calc_grad_l(nghnode,nghedge,sub_intvis,dn_via_ehen8)
    case(nderive_ehen8e)
        call calc_grad_l(nghnode,nghedge,sub_intvis,dn_via_ehen8e)
    case(nderive_scsl4)
        call calc_grad_l(nghnode,nghedge,sub_intvis,dn_via_scsl4)
    case(nderive_scsl4e)
        call calc_grad_l(nghnode,nghedge,sub_intvis,dn_via_scsl4e)
    case(nderive_scsl6)
        call calc_grad_l(nghnode,nghedge,sub_intvis,dn_via_scsl6)
    case(nderive_scsl6e)
        call calc_grad_l(nghnode,nghedge,sub_intvis,dn_via_scsl6e)                
    case default

    end select

end subroutine grad_l_dnvis


subroutine calc_grad_l(ngn,nge,sub_intvis,sub_dnvis)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : one
    use mod_constants, only : nfsf_vis_d3int,nfsf_vis_d3der
    use mod_datatypes, only : fld_array_t,var_block_t,top_block_t
    use mod_fieldvars, only : nblkcoms,blkcoms
    use mod_fieldvars, only : mb_ws,mb_sxyz,mb_vol,mb_dws
    use mod_interface, only : mb_var_pointer_create,mb_var_pointer_assign
    use mod_interface, only : mb_var_pointer_delete,calc_mb_dn_via_node3
    use mod_openmp
    implicit none
    integer(kind_int), intent(in) :: ngn,nge
    external                      :: sub_intvis,sub_dnvis
    integer(kind_int)          :: nc,nb,i,j,k,m,m1,m2,m3
    integer(kind_int)          :: st(3),ed(3)
    real(kind_real)            :: kx,ky,kz,ex,ey,ez,cx,cy,cz
    real(kind_real)            :: ovol,der(30)
    type(top_block_t), pointer :: top
    type(fld_array_t), pointer :: sxyz(:),vol(:),dws(:)
    type(var_block_t), pointer :: mb_uvwt(:)

    call mb_var_pointer_create(mb_uvwt,1,10)
    call mb_var_pointer_assign(mb_uvwt,2,7,mb_ws,5,10)
    call mb_var_pointer_assign(mb_uvwt,8,10,mb_dws,1,3)
    do m=2,10
        m1 = 3*m - 2
        m2 = m1 + 1
        m3 = m2 + 1

        call calc_mb_dn_via_node3(mb_uvwt,m,ngn,sub_intvis,nge,sub_dnvis,nfsf_vis_d3int,nfsf_vis_d3der,mb_dws,m1,1)
        call calc_mb_dn_via_node3(mb_uvwt,m,ngn,sub_intvis,nge,sub_dnvis,nfsf_vis_d3int,nfsf_vis_d3der,mb_dws,m2,2)
        call calc_mb_dn_via_node3(mb_uvwt,m,ngn,sub_intvis,nge,sub_dnvis,nfsf_vis_d3int,nfsf_vis_d3der,mb_dws,m3,3)
    end do
    call mb_var_pointer_delete(mb_uvwt)

    do nc=1,nblkcoms
        nb  =  blkcoms(nc)%nb
        top => blkcoms(nc)%top

        sxyz => mb_sxyz(nb)%fld
        vol  => mb_vol(nb)%fld
        dws  => mb_dws(nb)%fld

        st(:) = 1
        ed(:) = top%nijk(:)
!$OMP parallel do private(k,j,i,kx,ky,kz,ex,ey,ez,cx,cy,cz,ovol,m1,m2,m3,der)
        do k=st(3),ed(3)
        do j=st(2),ed(2)
        do i=st(1),ed(1)
            kx = sxyz(1)%r3d(i,j,k)
            ky = sxyz(2)%r3d(i,j,k)
            kz = sxyz(3)%r3d(i,j,k)

            ex = sxyz(4)%r3d(i,j,k)
            ey = sxyz(5)%r3d(i,j,k)
            ez = sxyz(6)%r3d(i,j,k)

            cx = sxyz(7)%r3d(i,j,k)
            cy = sxyz(8)%r3d(i,j,k)
            cz = sxyz(9)%r3d(i,j,k)

            ovol = one/vol(1)%r3d(i,j,k)

            do m=4,30
               der(m) = dws(m)%r3d(i,j,k)
            end do

            do m=2,10
                m1 = 3*m - 2
                m2 = m1 + 1
                m3 = m2 + 1

                dws(m1)%r3d(i,j,k) = ( kx*der(m1) + ex*der(m2) + cx*der(m3) ) * ovol
                dws(m2)%r3d(i,j,k) = ( ky*der(m1) + ey*der(m2) + cy*der(m3) ) * ovol
                dws(m3)%r3d(i,j,k) = ( kz*der(m1) + ez*der(m2) + cz*der(m3) ) * ovol
            end do
        end do
        end do
        end do
!$OMP end parallel do
    end do


end subroutine calc_grad_l

subroutine solve_Q
    use mod_kndconsts, only : kind_int,kind_real
    use mod_datatypes, only : fld_array_t,top_block_t
    use mod_fieldvars, only : npvs
    use mod_fieldvars, only : nblkcoms,blkcoms
    use mod_fieldvars, only : mb_dpv,mb_fmean
    implicit none
    integer(kind_int)          :: nc,nb,i,j,k,ierr
    integer(kind_int)          :: m,ni,nj,nk
    type(top_block_t), pointer :: top
    real(kind_real)            :: temp(5)
    type(fld_array_t), pointer :: dpv(:),fmean(:)

    do nc=1,nblkcoms
        nb  =  blkcoms(nc)%nb
        top => blkcoms(nc)%top

        ni = top%nijk(1)
        nj = top%nijk(2)
        nk = top%nijk(3)

        dpv   => mb_dpv(nb)%fld
        fmean => mb_fmean(nb)%fld

        do k=1,nk
        do j=1,nj
        do i=1,ni
            temp(1) = dpv(8)%r3d(i,j,k) - dpv(6)%r3d(i,j,k)
            temp(2) = dpv(3)%r3d(i,j,k) - dpv(7)%r3d(i,j,k)
            temp(3) = dpv(4)%r3d(i,j,k) - dpv(2)%r3d(i,j,k)
            temp(4) = sqrt( temp(1)**2.0 + temp(2)**2.0 + temp(3)**2.0 )
            temp(5) = dpv(5)%r3d(i,j,k)*dpv(9)%r3d(i,j,k) - dpv(6)%r3d(i,j,k)*dpv(8)%r3d(i,j,k) + &
                      dpv(1)%r3d(i,j,k)*dpv(9)%r3d(i,j,k) - dpv(3)%r3d(i,j,k)*dpv(7)%r3d(i,j,k) + &
                      dpv(1)%r3d(i,j,k)*dpv(5)%r3d(i,j,k) - dpv(2)%r3d(i,j,k)*dpv(4)%r3d(i,j,k)
        do m=1,npvs
            fmean(m)%r3d(i,j,k) = temp(m)
        end do
        end do
        end do
        end do

    end do

end subroutine solve_Q

subroutine map_sptocp
    use mod_kndconsts, only : kind_int
    use mod_constants, only : nintnon_muscl2pv
    use mod_constants, only : nintnon_wcns5pv,nintnon_wcns5cv
    use mod_constants, only : nintnon_wcns7pv,nintnon_wcns7cv
    use mod_constants, only : nintnon_hdcs5ei,nintnon_hdcs7ci
    use mod_constants, only : nintnon_hdcs5ci,nintnon_scsl3ci
    use mod_constants, only : nintnon_scsl5ci
    use mod_constants, only : nintnon_scsh3ci,nintnon_scsh5ci
    use mod_constants, only : nintnon_scsn2ci,nintnon_scsn3ci
    use mod_constants, only : nintnon_scsn4ci,nintnon_scsh3pi
    use mod_constants, only : nintnon_scsh5pi,nintnon_scsn2pi
    use mod_constants, only : nintnon_scsn3pi,nintnon_scsn4pi
    use mod_constants, only : nintnon_dcsh5ci,nintnon_dcsh7ci
    use mod_constants, only : nintnon_dcsh5pi,nintnon_dcsh7pi
    use mod_variables, only : nintnon,nghnode,nghedge
    use mod_fieldvars, only : nblkcoms,blkcomssp,npvs,mb_pv,neqn
    use mod_interface, only : cal_sptocp    
    implicit none
    integer(kind_int)         :: nc,nb
    external :: muscl2pvcc,wcns5pvcc,wcns5cvcc
    external :: wcns7pvcc,wcns7cvcc
    external :: hdcs5eicc,hdcs7cicc
    external :: hdcs5cicc,scsl3cicc
    external :: scsl5cicc
    external :: scsh3cicc,scsh5cicc
    external :: scsn2cicc,scsn3cicc
    external :: scsn4cicc,scsh3picc
    external :: scsh5picc,scsn2picc
    external :: scsn3picc,scsn4picc
    external :: dcsh5cicc,dcsh7cicc
    external :: dcsh5picc,dcsh7picc

    do nc=1,nblkcoms
        nb = blkcomssp(nc)%nb
        select case(nintnon)
        case(nintnon_muscl2pv)
           call cal_sptocp(nb,nghnode,nghedge,npvs,mb_pv,neqn,muscl2pvcc)
        case(nintnon_wcns5pv)
           call cal_sptocp(nb,nghnode,nghedge,npvs,mb_pv,neqn,wcns5pvcc)
        case(nintnon_wcns5cv)
           call cal_sptocp(nb,nghnode,nghedge,npvs,mb_pv,neqn,wcns5cvcc)
        case(nintnon_wcns7pv)
           call cal_sptocp(nb,nghnode,nghedge,npvs,mb_pv,neqn,wcns7pvcc)
        case(nintnon_wcns7cv)
           call cal_sptocp(nb,nghnode,nghedge,npvs,mb_pv,neqn,wcns7cvcc)
        case(nintnon_hdcs5ei)
           call cal_sptocp(nb,nghnode,nghedge,npvs,mb_pv,neqn,hdcs5eicc)
        case(nintnon_hdcs7ci)
           call cal_sptocp(nb,nghnode,nghedge,npvs,mb_pv,neqn,hdcs7cicc)
        case(nintnon_hdcs5ci)
           call cal_sptocp(nb,nghnode,nghedge,npvs,mb_pv,neqn,hdcs5cicc)
        case(nintnon_scsl3ci)
           call cal_sptocp(nb,nghnode,nghedge,npvs,mb_pv,neqn,scsl3cicc)
        case(nintnon_scsl5ci)
           call cal_sptocp(nb,nghnode,nghedge,npvs,mb_pv,neqn,scsl5cicc)
        case(nintnon_scsh3ci)
           call cal_sptocp(nb,nghnode,nghedge,npvs,mb_pv,neqn,scsh3cicc)
        case(nintnon_scsh5ci)
           call cal_sptocp(nb,nghnode,nghedge,npvs,mb_pv,neqn,scsh5cicc)
        case(nintnon_scsn2ci)
           call cal_sptocp(nb,nghnode,nghedge,npvs,mb_pv,neqn,scsn2cicc)       
        case(nintnon_scsn3ci)
           call cal_sptocp(nb,nghnode,nghedge,npvs,mb_pv,neqn,scsn3cicc)
        case(nintnon_scsn4ci)
           call cal_sptocp(nb,nghnode,nghedge,npvs,mb_pv,neqn,scsn4cicc) 
        case(nintnon_scsh3pi)
           call cal_sptocp(nb,nghnode,nghedge,npvs,mb_pv,neqn,scsh3picc)
        case(nintnon_scsh5pi)
           call cal_sptocp(nb,nghnode,nghedge,npvs,mb_pv,neqn,scsh5picc)
        case(nintnon_scsn2pi)
           call cal_sptocp(nb,nghnode,nghedge,npvs,mb_pv,neqn,scsn2picc)       
        case(nintnon_scsn3pi)
           call cal_sptocp(nb,nghnode,nghedge,npvs,mb_pv,neqn,scsn3picc)
        case(nintnon_scsn4pi)
           call cal_sptocp(nb,nghnode,nghedge,npvs,mb_pv,neqn,scsn4picc)
        case(nintnon_dcsh5ci)
           call cal_sptocp(nb,nghnode,nghedge,npvs,mb_pv,neqn,dcsh5cicc)
        case(nintnon_dcsh7ci)
           call cal_sptocp(nb,nghnode,nghedge,npvs,mb_pv,neqn,dcsh7cicc)
        case(nintnon_dcsh5pi)
           call cal_sptocp(nb,nghnode,nghedge,npvs,mb_pv,neqn,dcsh5picc)
        case(nintnon_dcsh7pi)
           call cal_sptocp(nb,nghnode,nghedge,npvs,mb_pv,neqn,dcsh7picc)
        case default        

        end select
    end do

end subroutine map_sptocp

subroutine cal_sptocp(nb,ngn,nge,npvs,mb_pv,neqn,sub_intnon)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : nsgl_buffer_pvs
    use mod_constants, only : nsgl_aver_art,nbc_inter_buf_pvs    
    use mod_constants, only : nfsf_con_d1nint,nsw_dir_close,nscmp_non
    use mod_constants, only : nbc_inter_scheme,nbc_intbc_scheme
    use mod_datatypes, only : var_block_t,fld_array_t
    use mod_variables, only : nsw_kdir,nscmp,poo,nghnode
    use mod_fieldvars, only : mb_topsp,mb_fsfsp,mb_pvfp,mb_sxyzsp,mb_vsl,mb_vst
    use mod_interface, only : pre_exchange_bc_var_sp,post_exchange_bc_var_sp
    use mod_interface, only : pre_exchange_cc_bc_var,post_exchange_cc_bc_var
    use mod_interface, only : average_bc_cc_var,exchange_singulars_cc,average_singulars_cc
    implicit none
    integer(kind_int), intent(in) :: nb
    integer(kind_int), intent(in) :: ngn,nge
    integer(kind_int), intent(in) :: npvs
    type(var_block_t),    pointer :: mb_pv(:)
    integer(kind_int), intent(in) :: neqn
    external                      :: sub_intnon
    real(kind_real), external  :: nolimiter
    integer(kind_int)          :: i,j,k,m,nfs,nfe,ierr
    integer(kind_int)          :: ni,nj,nk,stn,edn,ste,ede
    integer(kind_int)          :: stn0,edn0,ste0,ede0
    integer(kind_int)          :: nkst,nked,njed,ndec
    integer(kind_int)          :: chkid(npvs)
    real(kind_real)            :: chklim(npvs)
    real(kind_real)            :: vl(1:npvs),vr(1:npvs)
    type(fld_array_t), pointer :: sxyzsp(:),pv(:),pvfp(:),vsl(:),vst(:)
	type(fld_array_t), pointer :: fsfs(:),fsfe(:)
    real(kind_real),   pointer :: pvn(:,:),pvl(:,:),pvr(:,:),sn(:,:)
    
    !call pre_exchange_bc_var_sp(mb_pv,1,npvs,nghnode,nbc_inter_buf_pvs,nsgl_aver_art)
    !call post_exchange_bc_var_sp(mb_pv,1,npvs,nghnode,nbc_inter_buf_pvs,nsgl_aver_art)  
    
    chkid(:) = 0
    chklim(:) = 0.0
    ndec = 0  

    ni = mb_topsp(nb)%nijk(1)
    nj = mb_topsp(nb)%nijk(2)
    nk = mb_topsp(nb)%nijk(3)

    nkst = mb_topsp(nb)%ndst(3)
    nked = mb_topsp(nb)%nded(3)

    sxyzsp => mb_sxyzsp(nb)%fld    
    pv   => mb_pv(nb)%fld
    pvfp => mb_pvfp(nb)%fld

    ! I-direction

    fsfs => mb_fsfsp(nb,1)%fld
    fsfe => mb_fsfsp(nb,2)%fld    
	stn = 1 - ngn
    edn = ni + ngn
    ste = 1 - nge
    ede = ni + 1 + nge

    allocate(pvn(stn:edn,1:npvs), stat=ierr)
    allocate(pvl(stn:edn,1:npvs), stat=ierr)
    allocate(pvr(stn:edn,1:npvs), stat=ierr)
    allocate(sn(stn:edn,1:3), stat=ierr)
    
    if (nscmp > nscmp_non) then
        do k=1-ngn,nk+ngn
        do j=1-ngn,nj+ngn
        do i=1-ngn,ni+ngn
            pv(5)%r3d(i,j,k) = pv(5)%r3d(i,j,k) - poo
        end do
        end do
        end do 
    end if    
    
    do k=1,nk
    do j=1,nj
        do m=1,npvs
            do i=stn,edn
                pvn(i,m) = pv(m)%r3d(i,j,k)
            end do
        end do       

        do m=1,3
            do i=stn,edn
               sn(i,m) = sxyzsp(m)%r3d(i,j,k)
            end do
        end do        
        nfs = 1!fsfs(nfsf_con_d1nint)%i3d(1 ,j,k)
        nfe = 1!fsfe(nfsf_con_d1nint)%i3d(ni,j,k)
        call sub_intnon(ni,1,npvs,ngn,pvn,sn, &
                        nfs,nfe, &
                        nge,pvl,pvr,nolimiter, &
                        chkid,chklim,ndec)
       do m=1,npvs
            do i=1,ni+1
                pvfp(m)%r3d(i,j,k) = 0.5*(pvl(i,m) + pvr(i,m))
            end do
       end do
    end do
    end do
    
    !call pre_exchange_cc_bc_var(mb_pvfp,1,npvs,nghnode,nbc_inter_buf_pvs,nsgl_aver_art)
    !call exchange_singulars_cc(mb_pvfp,1,npvs,nsgl_buffer_pvs,nsgl_aver_art)
    !call post_exchange_cc_bc_var(mb_pvfp,1,npvs,nghnode,nbc_inter_buf_pvs,nsgl_aver_art)
    !call average_bc_cc_var(mb_pvfp,1,npvs,nbc_inter_buf_pvs,nsgl_aver_art)
	!call average_singulars_cc(mb_pvfp,1,npvs,nsgl_buffer_pvs,nsgl_aver_art)

    deallocate(sn, stat=ierr)
    deallocate(pvr, stat=ierr)
    deallocate(pvl, stat=ierr)
    deallocate(pvn, stat=ierr)

    ! J-direction
 
    fsfs => mb_fsfsp(nb,3)%fld
    fsfe => mb_fsfsp(nb,4)%fld    
	stn = 1 - ngn
    edn = nj + ngn
    ste = -nge
    ede = nj + 1 + nge

    allocate(pvn(stn:edn,1:npvs), stat=ierr)
    allocate(pvl(stn:edn,1:npvs), stat=ierr)
    allocate(pvr(stn:edn,1:npvs), stat=ierr)
    allocate(sn(stn:edn,1:3), stat=ierr)

    do k=1,nk
    do i=1,ni+1
        do m=1,npvs
            do j=stn,edn
                pvn(j,m) = pvfp(m)%r3d(i,j,k)
            end do
        end do     

        do m=1,3
            do j=stn,edn
               sn(j,m) = sxyzsp(m+3)%r3d(i,j,k)
            end do
        end do      

        nfs = 1!fsfs(nfsf_con_d1nint)%i3d(i, 1,k)
        nfe = 1!fsfe(nfsf_con_d1nint)%i3d(i,nj,k)
        call sub_intnon(nj,1,npvs,ngn,pvn,sn, &
                        nfs,nfe, &
                        nge,pvl,pvr,nolimiter, &
                        chkid,chklim,ndec)
       do m=1,npvs
            do j=1,nj+1
                pvfp(m)%r3d(i,j,k) = 0.5*(pvl(j,m) + pvr(j,m))
            end do
       end do
    end do
    end do
    
    !call pre_exchange_cc_bc_var(mb_pvfp,1,npvs,nghnode,nbc_inter_buf_pvs,nsgl_aver_art)
    !call exchange_singulars_cc(mb_pvfp,1,npvs,nsgl_buffer_pvs,nsgl_aver_art)
    !call post_exchange_cc_bc_var(mb_pvfp,1,npvs,nghnode,nbc_inter_buf_pvs,nsgl_aver_art)
    !call average_bc_cc_var(mb_pvfp,1,npvs,nbc_inter_buf_pvs,nsgl_aver_art)
	!call average_singulars_cc(mb_pvfp,1,npvs,nsgl_buffer_pvs,nsgl_aver_art)

    deallocate(sn, stat=ierr)
    deallocate(pvr, stat=ierr)
    deallocate(pvl, stat=ierr)
    deallocate(pvn, stat=ierr)

    ! K-direction

    fsfs => mb_fsfsp(nb,5)%fld
    fsfe => mb_fsfsp(nb,6)%fld    
	stn = 1 - ngn
    edn = nk + ngn
    ste = -nge
    ede = nk + 1 + nge

    allocate(pvn(stn:edn,1:npvs), stat=ierr)
    allocate(pvl(stn:edn,1:npvs), stat=ierr)
    allocate(pvr(stn:edn,1:npvs), stat=ierr)
    allocate(sn(stn:edn,1:3), stat=ierr)

    do j=1,nj+1
    do i=1,ni+1
        do m=1,npvs
            do k=stn,edn
                pvn(k,m) = pvfp(m)%r3d(i,j,k)
            end do
        end do      

        do m=1,3
            do k=stn,edn
               sn(k,m) = sxyzsp(m+6)%r3d(i,j,k)
            end do
        end do       

        nfs = 1!fsfs(nfsf_con_d1nint)%i3d(i,j, 1)
        nfe = 1!fsfe(nfsf_con_d1nint)%i3d(i,j,nk)
        call sub_intnon(nk,1,npvs,ngn,pvn,sn, &
                        nfs,nfe, &
                        nge,pvl,pvr,nolimiter, &
                        chkid,chklim,ndec)
       do m=1,npvs
            do k=1,nk+1
                pvfp(m)%r3d(i,j,k) = 0.5*(pvl(k,m) + pvr(k,m))
            end do
       end do
    end do
    end do

    deallocate(sn, stat=ierr)
    deallocate(pvr, stat=ierr)
    deallocate(pvl, stat=ierr)
    deallocate(pvn, stat=ierr)
    
    call pre_exchange_cc_bc_var(mb_pvfp,1,npvs,nghnode,nbc_inter_buf_pvs,nsgl_aver_art)
    call exchange_singulars_cc(mb_pvfp,1,npvs,nsgl_buffer_pvs,nsgl_aver_art)
    call post_exchange_cc_bc_var(mb_pvfp,1,npvs,nghnode,nbc_inter_buf_pvs,nsgl_aver_art)
    call average_bc_cc_var(mb_pvfp,1,npvs,nbc_inter_buf_pvs,nsgl_aver_art)
	call average_singulars_cc(mb_pvfp,1,npvs,nsgl_buffer_pvs,nsgl_aver_art)	    
    
    if (nscmp > nscmp_non) then
        do k=1,nk+1
        do j=1,nj+1
        do i=1,ni+1
            pvfp(5)%r3d(i,j,k) = pvfp(5)%r3d(i,j,k) + poo
        end do
        end do
        end do 
    end if
    
    if (nscmp > nscmp_non) then
        do k=1-ngn,nk+ngn
        do j=1-ngn,nj+ngn
        do i=1-ngn,ni+ngn
            pv(5)%r3d(i,j,k) = pv(5)%r3d(i,j,k) + poo
        end do
        end do
        end do 
    end if    

    if (nsw_kdir == nsw_dir_close) then
        do m=1,neqn
            do k=1,nk+1
            do j=1,nj+1
            do i=1,ni+1
                pvfp(m)%r3d(i,j,k) = pvfp(m)%r3d(i,j,nkst)
            end do
            end do
            end do
        end do
    end if
		
end subroutine cal_sptocp


