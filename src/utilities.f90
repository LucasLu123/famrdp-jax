

subroutine env_initialize
    use mod_kndconsts
    use mod_variables, only : parfile
    use mod_parallels
    implicit none
    integer(kind_int) :: ierr

#ifdef PARALLEL
    if (kind_int == kind_long) then
        kind_int_mpi = MPI_INTEGER
    else if (kind_int == kind_short) then
        kind_int_mpi = MPI_INTEGER2
    end if
    
    if (kind_real == kind_single) then
        kind_real_mpi = MPI_REAL
    else if (kind_real == kind_double) then
        kind_real_mpi = MPI_DOUBLE_PRECISION
    end if
    
    kind_char_mpi = MPI_CHARACTER

    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)
    call MPI_GET_PROCESSOR_NAME(procname,len_proc_name,ierr)
#endif

#ifdef PARALLEL
    if (myid == master) then
#endif
        write(*,*) "*********************************************************"
        write(*,*) "*    <HOSTA> - High-Order SimulaTor for Aerodynamics    *"
        write(*,*) "*     Copyright (c) 2005-2011 Fox Liu, Song Li.         *"
        write(*,*) "*     All Rights Reserved.                              *"
        write(*,*) "*********************************************************"
#ifdef PARALLEL
    end if
#endif

    call get_processors_name

    parfile = 'param.inp'
    
end subroutine env_initialize

subroutine get_processors_name
#ifdef PARALLEL
    use mod_kndconsts, only : kind_int
    use mod_parallels
    implicit none
    character(len=MPI_MAX_PROCESSOR_NAME) :: name
    integer(kind_int) :: np,ierr
    character(len=5)  :: strint
    
    allocate(procnames(0:numprocs-1),stat=ierr)
    do np=0,numprocs-1
        if (myid == np) then
           name = procname
        end if
        call MPI_BCAST(name,MPI_MAX_PROCESSOR_NAME,kind_char_mpi, &
                       np,MPI_COMM_WORLD,ierr)
        procnames(np) = name
    end do

    if (myid == master) then
       write(strint,'(i5)') numprocs
       write(*,'(3a)') "   The communicator starts ",trim(adjustl(strint))," processes."
       write(*,*) "========================================================="
       write(*,*) "               The processors name and pid"
       write(*,*) "---------------------------------------------------------"
       do np=0,numprocs-1
          if (mod(np,3) == 0) write(*,'(a)') "   "
          write(strint,'(i5)') np
          write(*, '(a18,1x)') trim(procnames(np))//'('//trim(adjustl(strint))//')'
          if (mod(np+1,3) == 0 .or. np==numprocs-1) write(*,*) ''
       end do
       write(*,*) "---------------------------------------------------------"
       write(*,*) "  The master pid:",master
       write(*,*) "========================================================="
    end if

    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif

end subroutine get_processors_name

subroutine env_finalize
    use mod_kndconsts, only : kind_int
    use mod_parallels
    implicit none
    integer(kind_int) :: ierr
 
#ifdef PARALLEL
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)

    if (myid == master) then
#endif
        write(*,*) "*********************************************************"
        write(*,*) "*    <HOSTA> Exit successfully.                         *"
        write(*,*) "*********************************************************"
#ifdef PARALLEL
    end if
#endif

#ifdef PARALLEL
    call MPI_FINALIZE(ierr)
#endif

end subroutine env_finalize

subroutine msg_seq_and_master(msg)
    use mod_parallels
    implicit none
    character(len=*) :: msg
    
#ifdef PARALLEL
    if (myid == master) then
#endif
        call msg_one_line(msg)
#ifdef PARALLEL
    end if
#endif

end subroutine msg_seq_and_master

subroutine msg_on_master(msg)
    use mod_parallels
    implicit none
    character(len=*) :: msg
    
#ifdef PARALLEL
    if (myid == master) then
        call msg_one_line(msg)
    end if
#endif

end subroutine msg_on_master

subroutine msg_one_line(msg)
    implicit none
    character(len=*) :: msg
    
    write(*,*) "HOSTA> "//trim(msg)//"."
        
end subroutine msg_one_line    

subroutine msg_prompt(msg)
    implicit none
    character(len=*) :: msg
    
    write(*,*) "HOSTA> "//trim(msg)//":"
        
end subroutine msg_prompt    

subroutine run_select(seqprocess,parprocess)
    use mod_parallels
    implicit none
    external :: seqprocess,parprocess
    
#ifdef PARALLEL
    call parprocess
#else
    call seqprocess
#endif

end subroutine run_select

subroutine run_parallel(subprocess)
    use mod_parallels
    implicit none
    external :: subprocess
    
#ifdef PARALLEL
    call subprocess
#endif

end subroutine run_parallel

subroutine run_seq_and_master(subprocess)
    use mod_parallels
    implicit none
    external :: subprocess
    
#ifdef PARALLEL
    if (myid == master) then
#endif
        call subprocess
#ifdef PARALLEL
    end if
#endif

end subroutine run_seq_and_master

subroutine run_on_master(subprocess)
    use mod_parallels
    implicit none
    external :: subprocess
    
#ifdef PARALLEL
    if (myid == master) then
        call subprocess
    end if
#endif

end subroutine run_on_master


subroutine run_on_others(subprocess)
    use mod_parallels
    implicit none
    external :: subprocess
    
#ifdef PARALLEL
    if (myid /= master) then
        call subprocess
    end if
#endif

end subroutine run_on_others

subroutine error_check(code,msg)
    use mod_kndconsts, only : kind_int
    use mod_parallels
    implicit none
    integer(kind_int) :: code
    character(len=*)  :: msg
    integer(kind_int) :: ierr
    
    if (code /= 0) then
        write(*,'(1x,a,i4)') "ERROR: "//trim(adjustl(msg))//",",code
#ifdef PARALLEL
        call MPI_ABORT(MPI_COMM_WORLD, code, ierr)
#else
        stop
#endif
    end if

end subroutine error_check


subroutine openfile(io_unit,name,sta,fm,acc)
    use mod_kndconsts, only : kind_int
    use mod_parallels
    implicit none
    integer(kind_int), intent(in) :: io_unit
    character(len=*) , intent(in) :: name,sta,fm,acc
    integer(kind_int) :: ierr
    
    ! sta: 'old','new','unknown'
    ! fm : 'formatted','binary'
    ! acc: 'append','sequential' 
#ifdef PARALLEL
    if (myid == master) then
#endif
        open(io_unit,file=name,status=sta, &
             access=acc,form=fm,iostat=ierr)
         
        call error_check(ierr,trim(adjustl(name))//" file cannot be opened")
#ifdef PARALLEL
    end if
#endif

end subroutine openfile

subroutine openfile_check(io_unit,name,sta,fm,acc,ierr)
    use mod_kndconsts, only : kind_int
    use mod_parallels
    implicit none
    integer(kind_int), intent(in) :: io_unit
    character(len=*) , intent(in) :: name,sta,fm,acc
    integer(kind_int), intent(out) :: ierr
    
    ! sta: 'old','new','unknown'
    ! fm : 'formatted','binary'
    ! acc: 'append','sequential' 
#ifdef PARALLEL
    if (myid == master) then
#endif
        open(io_unit,file=name,status=sta, &
             access=acc,form=fm,iostat=ierr)
#ifdef PARALLEL
    end if
#endif

end subroutine openfile_check

subroutine closefile(io_unit,name)
    use mod_kndconsts, only : kind_int
    use mod_parallels
    implicit none
    integer(kind_int), intent(in) :: io_unit
    character(len=*) , intent(in) :: name
    integer(kind_int) :: ierr

#ifdef PARALLEL
    if (myid == master) then
#endif
        close(io_unit,iostat=ierr)
    
        call error_check(ierr,trim(adjustl(name))//" file cannot be closed")
#ifdef PARALLEL
    end if
#endif

end subroutine closefile

subroutine deletefile(io_unit,name)
    use mod_kndconsts, only : kind_int
    use mod_parallels
    implicit none
    integer(kind_int), intent(in) :: io_unit
    character(len=*) , intent(in) :: name
    integer(kind_int) :: ierr

#ifdef PARALLEL
    if (myid == master) then
#endif
        close(io_unit,status="delete",iostat=ierr)
        
        call error_check(ierr,trim(adjustl(name))//" file cannot be deleted")        
#ifdef PARALLEL
    end if
#endif

end subroutine deletefile

subroutine renamefile(io_unit,name)
    use mod_kndconsts, only : kind_int
    use mod_parallels
    implicit none
    integer(kind_int), intent(in) :: io_unit
    character(len=*) , intent(in) :: name
    integer(kind_int) :: ierr,rename

#ifdef PARALLEL
    if (myid == master) then
#endif
        ierr=rename(name,trim(name)//".bak")
        
        call error_check(ierr,trim(adjustl(name))//" file cannot be renamed")        
#ifdef PARALLEL
    end if
#endif

end subroutine renamefile

subroutine fld_array_create(fldtype,array,st1,ed1,st2,ed2,st3,ed3,ierr)
    use mod_kndconsts, only : kind_int
    use mod_constants, only : fldtype_none,fldtype_i3d,fldtype_r3d
    use mod_datatypes, only : fld_array_t
    implicit none
    integer(kind_int), intent(in)  :: fldtype
    type(fld_array_t)              :: array
    integer(kind_int), intent(in)  :: st1,ed1
    integer(kind_int), intent(in)  :: st2,ed2
    integer(kind_int), intent(in)  :: st3,ed3
    integer(kind_int), intent(out) :: ierr

    select case(fldtype)
    case(fldtype_i3d)
        allocate(array%i3d(st1:ed1,st2:ed2,st3:ed3), stat=ierr)
        array%fldtype = fldtype_i3d
        nullify(array%r3d)
    case(fldtype_r3d)
        allocate(array%r3d(st1:ed1,st2:ed2,st3:ed3), stat=ierr)
        array%fldtype = fldtype_r3d
        nullify(array%i3d)
    case default
        ierr = 1
        array%fldtype = fldtype_none
    end select
    
end subroutine fld_array_create

subroutine fld_array_delete(array,ierr)
    use mod_kndconsts, only : kind_int
    use mod_constants, only : fldtype_none,fldtype_i3d,fldtype_r3d
    use mod_datatypes, only : fld_array_t
    implicit none
    type(fld_array_t)              :: array
    integer(kind_int), intent(out) :: ierr
    integer(kind_int) :: fldtype
    
    fldtype = array%fldtype
    
    select case(fldtype)
    case(fldtype_i3d)
        deallocate(array%i3d, stat=ierr)
        array%fldtype = fldtype_none
    case(fldtype_r3d)
        deallocate(array%r3d, stat=ierr)
        array%fldtype = fldtype_none
    case default
        ierr = 1
    end select
    
end subroutine fld_array_delete


subroutine var_array_create(array,nst,ned,fldtype,st1,ed1,st2,ed2,st3,ed3)
    use mod_kndconsts, only : kind_int
    use mod_datatypes, only : var_block_t
    implicit none
    type(var_block_t)             :: array
    integer(kind_int), intent(in) :: nst,ned
    integer(kind_int), intent(in) :: fldtype
    integer(kind_int), intent(in) :: st1,ed1
    integer(kind_int), intent(in) :: st2,ed2
    integer(kind_int), intent(in) :: st3,ed3
    integer(kind_int) :: m,ierr

    allocate(array%fld(nst:ned), stat=ierr)

    do m=nst,ned
       call fld_array_create(fldtype,array%fld(m),st1,ed1,st2,ed2,st3,ed3,ierr)
    end do
    
    call error_check(ierr,"var_array can not be created!")
    
end subroutine var_array_create

subroutine var_array_delete(array)
    use mod_kndconsts, only : kind_int
    use mod_datatypes, only : var_block_t
    implicit none
    type(var_block_t)             :: array
    integer(kind_int) :: nst,ned,m,ierr

    nst = lbound(array%fld,1)
    ned = ubound(array%fld,1)
    do m=nst,ned
       call fld_array_delete(array%fld(m),ierr)
    end do

    deallocate(array%fld, stat=ierr)

    call error_check(ierr,"var_array can not be deleted!")
    
end subroutine var_array_delete

subroutine bc_extend_inward(s_st,s_ed,s_nd,s_lr,nghost,st,ed)
    use mod_kndconsts,only : kind_int
    implicit none
    integer(kind_int), intent(in)  :: s_st(3),s_ed(3)
    integer(kind_int), intent(in)  :: s_nd,s_lr,nghost
    integer(kind_int), intent(out) :: st(3),ed(3)
    
    st(:) = s_st(:)
    ed(:) = s_ed(:)
    if (s_lr > 0) then
        st(s_nd) = st(s_nd) - nghost 
    else
        ed(s_nd) = ed(s_nd) + nghost 
    end if

end subroutine bc_extend_inward

subroutine bc_extend_outward(s_st,s_ed,s_nd,s_lr,nghost,st,ed)
    use mod_kndconsts,only : kind_int
    implicit none
    integer(kind_int), intent(in)  :: s_st(3),s_ed(3)
    integer(kind_int), intent(in)  :: s_nd,s_lr,nghost
    integer(kind_int), intent(out) :: st(3),ed(3)
    
    st(:) = s_st(:)
    ed(:) = s_ed(:)
    if (s_lr > 0) then
        ed(s_nd) = ed(s_nd) + nghost 
    else
        st(s_nd) = st(s_nd) - nghost 
    end if

end subroutine bc_extend_outward

subroutine bc_extend_inoutward(s_st,s_ed,s_nd,s_lr,nghost,st,ed)
    use mod_kndconsts,only : kind_int
    implicit none
    integer(kind_int), intent(in)  :: s_st(3),s_ed(3)
    integer(kind_int), intent(in)  :: s_nd,s_lr,nghost
    integer(kind_int), intent(out) :: st(3),ed(3)
    
    st(:) = s_st(:)
    ed(:) = s_ed(:)
    st(s_nd) = st(s_nd) - nghost 
    ed(s_nd) = ed(s_nd) + nghost 

end subroutine bc_extend_inoutward

subroutine bc_extend_ghost(s_st,s_ed,s_nd,s_lr,ngst,nged,st,ed)
    use mod_kndconsts,only : kind_int
    implicit none
    integer(kind_int), intent(in)  :: s_st(3),s_ed(3)
    integer(kind_int), intent(in)  :: s_nd,s_lr,ngst,nged
    integer(kind_int), intent(out) :: st(3),ed(3)
    
    st(:) = s_st(:)
    ed(:) = s_ed(:)
    if (s_lr > 0) then
        st(s_nd) = st(s_nd) + ngst
        ed(s_nd) = ed(s_nd) + nged 
    else
        st(s_nd) = st(s_nd) - nged 
        ed(s_nd) = ed(s_nd) - ngst 
    end if

end subroutine bc_extend_ghost


subroutine set_bc_inters
    use mod_kndconsts, only : kind_int
    use mod_constants, only : bc_cut1to1
    use mod_constants, only : nbc_inter_buf_max,nbc_inter_buf_pvs
    use mod_constants, only : nbc_inter_buf_dqc,nbc_inter_buf_dpv
    use mod_constants, only : nbc_inter_buf_qts,nbc_inter_buf_dqt
    use mod_constants, only : nbc_inter_buf_dtur,nvis_ns_lam
    use mod_constants, only : subc_cut_patched
    use mod_datatypes, only : bc_region_t,top_block_t
    use mod_variables, only : nghnode,nvis
    use mod_fieldvars, only : nblocks,mb_top,npvs,neqt,nqvst
    use mod_fieldvars, only : ninters,ninterf,inters,nmsg_intf
    use mod_parallels
    implicit none
    integer(kind_int) :: nint,nst,ned,bctype,subtype,ierr
    integer(kind_int) :: nbs,nrs,s_st(3),s_ed(3),s_nd,s_lr
    integer(kind_int) :: nbt,nrt,t_st(3),t_ed(3),t_nd,t_lr
    integer(kind_int) :: nregs,st(3),ed(3),id_src,id_des
    integer(kind_int) :: nbc_inter_bufs
    integer(kind_int) :: m,nbufpars(3,nbc_inter_buf_max)
    type(top_block_t), pointer :: top
    type(bc_region_t), pointer :: reg

    nint = 0
    do nbs=1,nblocks
        top => mb_top(nbs)

        nregs = top%nregions
        do nrs=1,nregs
            reg => top%bcs(nrs)
            
            bctype  = reg%bctype
            subtype = reg%subtype
            if (bctype == bc_cut1to1 .and. subtype /= subc_cut_patched) then
               nint = nint + 1
               reg%nint = nint
            end if      
        end do
    end do
    ninterf = nint

    do nbs=1,nblocks
        top => mb_top(nbs)

        nregs = top%nregions
        do nrs=1,nregs
            reg => top%bcs(nrs)
            
            bctype  = reg%bctype
            subtype = reg%subtype
            if (bctype == bc_cut1to1 .and. subtype == subc_cut_patched) then
               nint = nint + 1
               reg%nint = nint
            end if      
        end do
    end do
    
    ninters = nint
    allocate(inters(ninters),stat=ierr)
    
    nint = 0
    do nbs=1,nblocks
        top => mb_top(nbs)

        nregs = top%nregions
        do nrs=1,nregs
            reg => top%bcs(nrs)
            
            nbt = reg%nbt
            nrt = reg%nrt
            
            bctype  = reg%bctype
            subtype = reg%subtype
            if (bctype == bc_cut1to1  .and. subtype /= subc_cut_patched) then
               nint = nint + 1
               inters(nint)%bc => reg
               inters(nint)%nint = mb_top(nbt)%bcs(nrt)%nint
               nullify(inters(nint)%dat)
            end if      
        end do
    end do

    ! 为对接信息通讯分配内存
    nbufpars(:,nbc_inter_buf_pvs ) = (/nghnode,1,npvs/)
    nbufpars(:,nbc_inter_buf_dqc ) = (/1      ,1,npvs/)
    nbufpars(:,nbc_inter_buf_dpv ) = (/nghnode,1,12  /)
    nbufpars(:,nbc_inter_buf_qts ) = (/nghnode,1,nqvst/)
    nbufpars(:,nbc_inter_buf_dqt ) = (/1      ,1,neqt  /)
    nbufpars(:,nbc_inter_buf_dtur) = (/nghnode,1,neqt*3/)

    if (nvis > nvis_ns_lam) then
        nbc_inter_bufs = nbc_inter_buf_max
    else
        nbc_inter_bufs = nbc_inter_buf_dpv
    end if

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
        
        allocate(inters(nint)%buf(nbc_inter_bufs),stat=ierr)

        do m=1,nbc_inter_bufs
#ifdef PARALLEL
            id_src = mb_top(nbs)%pid
            id_des = mb_top(nbt)%pid
            if (id_src /= id_des) then                         
                if (myid == id_src .or. myid == id_des) then
                    call bc_extend_inward(s_st,s_ed,s_nd,s_lr,nbufpars(1,m),st,ed)
                    nst = nbufpars(2,m)
                    ned = nbufpars(3,m)
                    allocate(inters(nint)%buf(m)%dat(st(1):ed(1), &
                                                     st(2):ed(2), &
                                                     st(3):ed(3), &
                                                     nst:ned),stat=ierr)
                end if
            else
                if (myid == id_src) then
#endif
                    call bc_extend_inward(s_st,s_ed,s_nd,s_lr,nbufpars(1,m),st,ed)
                    nst = nbufpars(2,m)
                    ned = nbufpars(3,m)
                    allocate(inters(nint)%buf(m)%dat(st(1):ed(1), &
                                                     st(2):ed(2), &
                                                     st(3):ed(3), &
                                                     nst:ned),stat=ierr)
#ifdef PARALLEL
                end if
            end if
#endif
        end do
    end do

#ifdef PARALLEL
    nmsg_int = 0
    do nint=1,ninterf
        nbs    = inters(nint)%bc%nbs
        nbt    = inters(nint)%bc%nbt
        id_src = mb_top(nbs)%pid
        id_des = mb_top(nbt)%pid

        if (id_src /= id_des) then  ! 仅当源进程与目标进程不同时，才需要MPI通信
            if (myid == id_src .or. myid == id_des) then
                nmsg_int = nmsg_int + 1
            end if
        end if
    end do

    nmsg_intf = nmsg_int
    if (nmsg_int > 0) then
        allocate(request_int(nmsg_int),stat=ierr)
        allocate(status_int(MPI_STATUS_SIZE, nmsg_int),stat=ierr)
    end if
#endif
end subroutine set_bc_inters

subroutine set_bc_inters_cc
    use mod_kndconsts, only : kind_int
    use mod_constants, only : bc_cut1to1
    use mod_constants, only : nbc_inter_buf_max,nbc_inter_buf_pvs
    use mod_constants, only : nbc_inter_buf_dqc,nbc_inter_buf_dpv
    use mod_constants, only : nbc_inter_buf_qts,nbc_inter_buf_dqt
    use mod_constants, only : nbc_inter_buf_dtur,nvis_ns_lam
    use mod_datatypes, only : bc_region_t,top_block_t
    use mod_variables, only : nghnode,nvis
    use mod_fieldvars, only : nblocks,mb_top,mb_topc,mb_topsp,npvs,neqt,nqvst
    use mod_fieldvars, only : ninters,ninterf,inters,interscc,interssp,nmsg_intf
    use mod_parallels
    implicit none
    integer(kind_int) :: nint,nst,ned,bctype,subtype,ierr
    integer(kind_int) :: nbs,nrs,s_st(3),s_ed(3),s_nd,s_lr
    integer(kind_int) :: nbt,nrt,t_st(3),t_ed(3),t_nd,t_lr
    integer(kind_int) :: nregs,st(3),ed(3),id_src,id_des
    integer(kind_int) :: nbc_inter_bufs
    integer(kind_int) :: m,nbufpars(3,nbc_inter_buf_max)
    type(top_block_t), pointer :: top
    type(bc_region_t), pointer :: reg

    nint = 0
    do nbs=1,nblocks
        top => mb_top(nbs)

        nregs = top%nregions
        do nrs=1,nregs
            reg => top%bcs(nrs)
            
            bctype  = reg%bctype
            subtype = reg%subtype
            if (bctype == bc_cut1to1) then
               nint = nint + 1
               reg%nint = nint
            end if      
        end do
    end do
    ninterf = nint
    ninters = nint
    
    allocate(inters(ninterf),stat=ierr)
    allocate(interscc(ninterf),stat=ierr)
    allocate(interssp(ninterf),stat=ierr)
    
    nint = 0
    do nbs=1,nblocks
        top => mb_top(nbs)

        nregs = top%nregions
        do nrs=1,nregs
            reg => top%bcs(nrs)
            
            nbt = reg%nbt
            nrt = reg%nrt
            
            bctype  = reg%bctype
            if (bctype == bc_cut1to1) then
               nint = nint + 1
               inters(nint)%bc => reg
               inters(nint)%nint = mb_top(nbt)%bcs(nrt)%nint
               nullify(inters(nint)%dat)
            end if      
        end do
    end do    

    ! 为对接信息通讯分配内存
    nbufpars(:,nbc_inter_buf_pvs ) = (/nghnode,1,npvs/)
    nbufpars(:,nbc_inter_buf_dqc ) = (/1      ,1,npvs/)
    nbufpars(:,nbc_inter_buf_dpv ) = (/nghnode,1,12  /)
    nbufpars(:,nbc_inter_buf_qts ) = (/nghnode,1,nqvst/)
    nbufpars(:,nbc_inter_buf_dqt ) = (/1      ,1,neqt  /)
    nbufpars(:,nbc_inter_buf_dtur) = (/nghnode,1,neqt*3/)

    if (nvis > nvis_ns_lam) then
        nbc_inter_bufs = nbc_inter_buf_max
    else
        nbc_inter_bufs = nbc_inter_buf_dpv
    end if

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
        
        allocate(inters(nint)%buf(nbc_inter_bufs),stat=ierr)

        do m=1,nbc_inter_bufs
#ifdef PARALLEL
            id_src = mb_top(nbs)%pid
            id_des = mb_top(nbt)%pid
            if (id_src /= id_des) then                         
                if (myid == id_src .or. myid == id_des) then
                    call bc_extend_inward(s_st,s_ed,s_nd,s_lr,nbufpars(1,m),st,ed)
                    nst = nbufpars(2,m)
                    ned = nbufpars(3,m)
                    allocate(inters(nint)%buf(m)%dat(st(1):ed(1), &
                                                     st(2):ed(2), &
                                                     st(3):ed(3), &
                                                     nst:ned),stat=ierr)
                end if
            else
                if (myid == id_src) then
#endif
                    call bc_extend_inward(s_st,s_ed,s_nd,s_lr,nbufpars(1,m),st,ed)
                    nst = nbufpars(2,m)
                    ned = nbufpars(3,m)
                    allocate(inters(nint)%buf(m)%dat(st(1):ed(1), &
                                                     st(2):ed(2), &
                                                     st(3):ed(3), &
                                                     nst:ned),stat=ierr)
#ifdef PARALLEL
                end if
            end if
#endif
        end do
    end do    
    
    nint = 0
    do nbs=1,nblocks
        top => mb_topc(nbs)

        nregs = top%nregions
        do nrs=1,nregs
            reg => top%bcs(nrs)
            
            nbt = reg%nbt
            nrt = reg%nrt
            
            bctype  = reg%bctype
            if (bctype == bc_cut1to1) then
               nint = nint + 1
               interscc(nint)%bc => reg
               interscc(nint)%nint = mb_topc(nbt)%bcs(nrt)%nint
               nullify(interscc(nint)%dat)
            end if      
        end do
    end do    

    if (nvis > nvis_ns_lam) then
        nbc_inter_bufs = nbc_inter_buf_max
    else
        nbc_inter_bufs = nbc_inter_buf_dpv
    end if

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
        
        allocate(interscc(nint)%buf(nbc_inter_bufs),stat=ierr)

        do m=1,nbc_inter_bufs
#ifdef PARALLEL
            id_src = mb_topc(nbs)%pid
            id_des = mb_topc(nbt)%pid
            if (id_src /= id_des) then                         
                if (myid == id_src .or. myid == id_des) then
                    call bc_extend_inward(s_st,s_ed,s_nd,s_lr,nbufpars(1,m),st,ed)
                    nst = nbufpars(2,m)
                    ned = nbufpars(3,m)
                    allocate(interscc(nint)%buf(m)%dat(st(1):ed(1), &
                                                       st(2):ed(2), &
                                                       st(3):ed(3), &
                                                       nst:ned),stat=ierr)
                end if
            else
                if (myid == id_src) then
#endif
                    call bc_extend_inward(s_st,s_ed,s_nd,s_lr,nbufpars(1,m),st,ed)
                    nst = nbufpars(2,m)
                    ned = nbufpars(3,m)
                    allocate(interscc(nint)%buf(m)%dat(st(1):ed(1), &
                                                       st(2):ed(2), &
                                                       st(3):ed(3), &
                                                       nst:ned),stat=ierr)
#ifdef PARALLEL
                end if
            end if
#endif
        end do
    end do
    
    nint = 0
    do nbs=1,nblocks
        top => mb_topsp(nbs)

        nregs = top%nregions
        do nrs=1,nregs
            reg => top%bcs(nrs)
            
            nbt = reg%nbt
            nrt = reg%nrt
            
            bctype  = reg%bctype
            if (bctype == bc_cut1to1) then
               nint = nint + 1
               interssp(nint)%bc => reg
               interssp(nint)%nint = mb_topsp(nbt)%bcs(nrt)%nint
               nullify(interssp(nint)%dat)
            end if      
        end do
    end do
    
    do nint=1,ninterf
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
        
        allocate(interssp(nint)%buf(nbc_inter_bufs),stat=ierr)

        do m=1,nbc_inter_bufs
#ifdef PARALLEL
            id_src = mb_topsp(nbs)%pid
            id_des = mb_topsp(nbt)%pid
            if (id_src /= id_des) then                         
                if (myid == id_src .or. myid == id_des) then
                    call bc_extend_inward(s_st,s_ed,s_nd,s_lr,nbufpars(1,m),st,ed)
                    nst = nbufpars(2,m)
                    ned = nbufpars(3,m)
                    allocate(interssp(nint)%buf(m)%dat(st(1):ed(1), &
                                                       st(2):ed(2), &
                                                       st(3):ed(3), &
                                                       nst:ned),stat=ierr)
                end if
            else
                if (myid == id_src) then
#endif
                    call bc_extend_inward(s_st,s_ed,s_nd,s_lr,nbufpars(1,m),st,ed)
                    nst = nbufpars(2,m)
                    ned = nbufpars(3,m)
                    allocate(interssp(nint)%buf(m)%dat(st(1):ed(1), &
                                                       st(2):ed(2), &
                                                       st(3):ed(3), &
                                                       nst:ned),stat=ierr)
#ifdef PARALLEL
                end if
            end if
#endif
        end do
    end do
    
#ifdef PARALLEL
    nmsg_int = 0
    do nint=1,ninterf
        nbs    = interscc(nint)%bc%nbs
        nbt    = interscc(nint)%bc%nbt
        id_src = mb_topc(nbs)%pid
        id_des = mb_topc(nbt)%pid

        if (id_src /= id_des) then  ! 仅当源进程与目标进程不同时，才需要MPI通信
            if (myid == id_src .or. myid == id_des) then
                nmsg_int = nmsg_int + 1
            end if
        end if
    end do

    nmsg_intf = nmsg_int
    if (nmsg_int > 0) then
        allocate(request_int(nmsg_int),stat=ierr)
        allocate(status_int(MPI_STATUS_SIZE, nmsg_int),stat=ierr)
    end if
#endif    
end subroutine set_bc_inters_cc

subroutine set_bc_inters_patched
    use mod_kndconsts, only : kind_int
    use mod_constants, only : bc_cut1to1
    use mod_constants, only : nbc_inter_buf_max,nbc_inter_buf_pvs
    use mod_constants, only : nbc_inter_buf_dqc,nbc_inter_buf_dpv
    use mod_constants, only : nbc_inter_buf_qts,nbc_inter_buf_dqt
    use mod_constants, only : nbc_inter_buf_dtur,nvis_ns_lam
    use mod_constants, only : subc_cut_patched
    use mod_datatypes, only : bc_region_t,top_block_t
    use mod_variables, only : nghnode,nvis
    use mod_fieldvars, only : nblocks,mb_top,npvs,neqt,nqvst
    use mod_fieldvars, only : ninters,ninterf,inters,nmsg_intf
    use mod_parallels
    implicit none
    integer(kind_int) :: nint,nst,ned,bctype,subtype,ierr
    integer(kind_int) :: nbs,nrs,s_st(3),s_ed(3),s_nd,s_lr
    integer(kind_int) :: nbt,nrt,t_st(3),t_ed(3),t_nd,t_lr
    integer(kind_int) :: nregs,st(3),ed(3),id_src,id_des
    integer(kind_int) :: nbc_inter_bufs
    integer(kind_int) :: m,nbufpars(3,nbc_inter_buf_max)
    type(top_block_t), pointer :: top
    type(bc_region_t), pointer :: reg

    nint = ninterf
    do nbs=1,nblocks
        top => mb_top(nbs)

        nregs = top%nregions
        do nrs=1,nregs
            reg => top%bcs(nrs)
            
            nbt = reg%nbt
            nrt = reg%nrt
            
            bctype  = reg%bctype
            subtype = reg%subtype
            if (bctype == bc_cut1to1  .and. subtype == subc_cut_patched) then
               nint = nint + 1
               inters(nint)%bc => reg
               inters(nint)%nint = mb_top(nbt)%bcs(nrt)%nint
               nullify(inters(nint)%dat)
            end if      
        end do
    end do

    ! 为对接信息通讯分配内存
    nbufpars(:,nbc_inter_buf_pvs ) = (/nghnode,1,npvs/)
    nbufpars(:,nbc_inter_buf_dqc ) = (/1      ,1,npvs/)
    nbufpars(:,nbc_inter_buf_dpv ) = (/nghnode,1,12  /)
    nbufpars(:,nbc_inter_buf_qts ) = (/nghnode,1,nqvst/)
    nbufpars(:,nbc_inter_buf_dqt ) = (/1      ,1,neqt  /)
    nbufpars(:,nbc_inter_buf_dtur) = (/nghnode,1,neqt*3/)

    if (nvis > nvis_ns_lam) then
        nbc_inter_bufs = nbc_inter_buf_max
    else
        nbc_inter_bufs = nbc_inter_buf_dpv
    end if

    do nint=ninterf+1,ninters
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
        
        allocate(inters(nint)%buf(nbc_inter_bufs),stat=ierr)

        do m=1,nbc_inter_bufs
#ifdef PARALLEL
            id_src = mb_top(nbs)%pid
            id_des = mb_top(nbt)%pid
            if (id_src /= id_des) then                         
                if (myid == id_src .or. myid == id_des) then
                    call bc_extend_inward(s_st,s_ed,s_nd,s_lr,nbufpars(1,m),st,ed)
                    nst = nbufpars(2,m)
                    ned = nbufpars(3,m)
                    allocate(inters(nint)%buf(m)%dat(st(1):ed(1), &
                                                     st(2):ed(2), &
                                                     st(3):ed(3), &
                                                     nst:ned),stat=ierr)
                end if
            else
                if (myid == id_src) then
#endif
                    call bc_extend_inward(s_st,s_ed,s_nd,s_lr,nbufpars(1,m),st,ed)
                    nst = nbufpars(2,m)
                    ned = nbufpars(3,m)
                    allocate(inters(nint)%buf(m)%dat(st(1):ed(1), &
                                                     st(2):ed(2), &
                                                     st(3):ed(3), &
                                                     nst:ned),stat=ierr)
#ifdef PARALLEL
                end if
            end if
#endif
        end do
    end do

#ifdef PARALLEL
    nmsg_int = nmsg_intf
    do nint=ninterf+1,ninters
        nbs    = inters(nint)%bc%nbs
        nbt    = inters(nint)%bc%nbt
        id_src = mb_top(nbs)%pid
        id_des = mb_top(nbt)%pid

        if (id_src /= id_des) then  ! 仅当源进程与目标进程不同时，才需要MPI通信
            if (myid == id_src .or. myid == id_des) then
                nmsg_int = nmsg_int + 1
            end if
        end if        
    end do

    if (nmsg_int > 0) then
        allocate(request_int(nmsg_int),stat=ierr)
        allocate(status_int(MPI_STATUS_SIZE, nmsg_int),stat=ierr)
    end if
#endif
end subroutine set_bc_inters_patched

subroutine create_buff_offsets
    use mod_kndconsts, only : kind_int
    use mod_variables, only : nghnode
    use mod_fieldvars, only : mb_top
    use mod_fieldvars, only : ninters,ninterf
    use mod_fieldvars, only : inters,patched
    use mod_parallels
    implicit none
    integer(kind_int) :: nint,bc_id,nst,ned,ierr
    integer(kind_int) :: nbs,nrs,s_st(3),s_ed(3),s_nd,s_lr
    integer(kind_int) :: st(3),ed(3),id_src
    integer(kind_int), parameter :: nghp = 5


    allocate(patched(1:ninters-ninterf),stat=ierr)

    do nint=ninterf+1,ninters
        nbs     = inters(nint)%bc%nbs
        nrs     = inters(nint)%bc%nrs
        s_st(:) = inters(nint)%bc%s_st(:)
        s_ed(:) = inters(nint)%bc%s_ed(:)
        s_nd    = inters(nint)%bc%s_nd
        s_lr    = inters(nint)%bc%s_lr

        bc_id   = nint - ninterf

#ifdef PARALLEL
        id_src = mb_top(nbs)%pid
        if (myid == id_src) then
#endif
            nst = 1
            ned = nghp+1
            allocate(patched(bc_id)%id1(s_st(1):s_ed(1), &
                                         s_st(2):s_ed(2), &
                                         s_st(3):s_ed(3), &
                                         nst:ned),stat=ierr)
            allocate(patched(bc_id)%id2(s_st(1):s_ed(1), &
                                         s_st(2):s_ed(2), &
                                         s_st(3):s_ed(3), &
                                         nst:ned),stat=ierr)
            call bc_extend_outward(s_st,s_ed,s_nd,s_lr,nghnode,st,ed)
            allocate(patched(bc_id)%dats( st(1):ed(1), &
                                          st(2):ed(2), &
                                          st(3):ed(3), &
                                          1:3),stat=ierr)
            allocate(patched(bc_id)%datt( st(1):ed(1), &
                                          st(2):ed(2), &
                                          st(3):ed(3), &
                                          1:3),stat=ierr)             
            
#ifdef PARALLEL
        end if
#endif
    end do
end subroutine create_buff_offsets

subroutine set_block_coms
    use mod_kndconsts, only : kind_int
    use mod_datatypes, only : top_block_t
    use mod_fieldvars, only : nblocks,mb_top
    use mod_fieldvars, only : nblkcoms,blkcoms,ntotpts
    use mod_parallels
    implicit none
    integer(kind_int) :: nb,pid,nc,ierr
    type(top_block_t), pointer :: top

    nc = 0
    do nb=1,nblocks
        top => mb_top(nb)

        pid = top%pid
        
#ifdef PARALLEL
        if (myid == pid) then
#endif       
            nc = nc + 1
#ifdef PARALLEL
        end if
#endif       
    end do
    
    nblkcoms = nc
    allocate(blkcoms(nblkcoms),stat=ierr)
    
    nc = 0
    do nb=1,nblocks
        top => mb_top(nb)

        pid = top%pid
        
#ifdef PARALLEL
        if (myid == pid) then
#endif       
            nc = nc + 1
            blkcoms(nc)%nb  = nb
            blkcoms(nc)%top => top
#ifdef PARALLEL
        end if
#endif       
    end do

    ntotpts = 0
    do nb=1,nblocks
        top => mb_top(nb)
        ntotpts = ntotpts + product(top%nijk(:)) 
    end do

end subroutine set_block_coms

subroutine set_block_coms_cc
    use mod_kndconsts, only : kind_int
    use mod_datatypes, only : top_block_t
    use mod_fieldvars, only : nblocks,mb_top,mb_topc,mb_topsp
    use mod_fieldvars, only : nblkcoms,blkcoms,ntotpts,blkcomscc,ntotptscc,blkcomssp,ntotptssp
    use mod_parallels
    implicit none
    integer(kind_int) :: nb,pid,nc,ierr
    type(top_block_t), pointer :: top

    nc = 0
    do nb=1,nblocks
        top => mb_top(nb)

        pid = top%pid
        
#ifdef PARALLEL
        if (myid == pid) then
#endif       
            nc = nc + 1
#ifdef PARALLEL
        end if
#endif       
    end do
    
    nblkcoms = nc
    allocate(blkcoms(nblkcoms),stat=ierr)
    
    nc = 0
    do nb=1,nblocks
        top => mb_top(nb)

        pid = top%pid
        
#ifdef PARALLEL
        if (myid == pid) then
#endif       
            nc = nc + 1
            blkcoms(nc)%nb  = nb
            blkcoms(nc)%top => top
#ifdef PARALLEL
        end if
#endif       
    end do

    ntotpts = 0
    do nb=1,nblocks
        top => mb_top(nb)
        ntotpts = ntotpts + product(top%nijk(:)) 
    end do
!----------------------    

    nc = 0
    do nb=1,nblocks
        top => mb_topc(nb)

        pid = top%pid
        
#ifdef PARALLEL
        if (myid == pid) then
#endif       
            nc = nc + 1
#ifdef PARALLEL
        end if
#endif       
    end do
    
    nblkcoms = nc
    allocate(blkcomscc(nblkcoms),stat=ierr)
    
    nc = 0
    do nb=1,nblocks
        top => mb_topc(nb)

        pid = top%pid
        
#ifdef PARALLEL
        if (myid == pid) then
#endif       
            nc = nc + 1
            blkcomscc(nc)%nb  = nb
            blkcomscc(nc)%top => top
#ifdef PARALLEL
        end if
#endif       
    end do

    ntotpts = 0
    do nb=1,nblocks
        top => mb_topc(nb)
        ntotpts = ntotpts + product(top%nijk(:)) 
    end do  
    
!----------------------    

    nc = 0
    do nb=1,nblocks
        top => mb_topsp(nb)

        pid = top%pid
        
#ifdef PARALLEL
        if (myid == pid) then
#endif       
            nc = nc + 1
#ifdef PARALLEL
        end if
#endif       
    end do
    
    nblkcoms = nc
    allocate(blkcomssp(nblkcoms),stat=ierr)
    
    nc = 0
    do nb=1,nblocks
        top => mb_topsp(nb)

        pid = top%pid
        
#ifdef PARALLEL
        if (myid == pid) then
#endif       
            nc = nc + 1
            blkcomssp(nc)%nb  = nb
            blkcomssp(nc)%top => top
#ifdef PARALLEL
        end if
#endif       
    end do

    ntotpts = 0
    do nb=1,nblocks
        top => mb_topsp(nb)
        ntotpts = ntotpts + product(top%nijk(:)) 
    end do    

end subroutine set_block_coms_cc

!subroutine set_block_coms_cc
!    use mod_kndconsts, only : kind_int
!    use mod_datatypes, only : top_block_t
!    use mod_fieldvars, only : nblocks,mb_top,mb_topc,mb_topsp
!    use mod_fieldvars, only : nblkcoms,blkcoms,ntotpts,blkcomscc,ntotptscc,blkcomssp,ntotptssp
!    use mod_parallels
!    implicit none
!    integer(kind_int) :: nb,pid,nc,ierr
!    type(top_block_t), pointer :: top
!
!    nc = 0
!    do nb=1,nblocks
!        top => mb_topc(nb)
!
!        pid = top%pid
!        
!#ifdef PARALLEL
!        if (myid == pid) then
!#endif       
!            nc = nc + 1
!#ifdef PARALLEL
!        end if
!#endif       
!    end do
!    
!    nblkcoms = nc
!    allocate(blkcoms(nblkcoms),stat=ierr)
!    allocate(blkcomscc(nblkcoms),stat=ierr)
!    allocate(blkcomssp(nblkcoms),stat=ierr)
!    
!    nc = 0
!    do nb=1,nblocks
!        top => mb_top(nb)
!
!        pid = top%pid
!        
!#ifdef PARALLEL
!        if (myid == pid) then
!#endif       
!            nc = nc + 1
!            blkcoms(nc)%nb  = nb
!            blkcoms(nc)%top => top
!#ifdef PARALLEL
!        end if
!#endif       
!    end do
!
!    ntotpts = 0
!    do nb=1,nblocks
!        top => mb_top(nb)
!        ntotpts = ntotpts + product(top%nijk(:)) 
!    end do    
!    
!    nc = 0
!    do nb=1,nblocks
!        top => mb_topc(nb)
!
!        pid = top%pid
!        
!#ifdef PARALLEL
!        if (myid == pid) then
!#endif       
!            nc = nc + 1
!            blkcomscc(nc)%nb  = nb
!            blkcomscc(nc)%top => top
!#ifdef PARALLEL
!        end if
!#endif       
!    end do
!
!    ntotptscc = 0
!    do nb=1,nblocks
!        top => mb_topc(nb)
!        ntotptscc = ntotptscc + product(top%nijk(:)) 
!    end do
!    
!    nc = 0
!    do nb=1,nblocks
!        top => mb_topsp(nb)
!
!        pid = top%pid
!        
!#ifdef PARALLEL
!        if (myid == pid) then
!#endif       
!            nc = nc + 1
!            blkcomssp(nc)%nb  = nb
!            blkcomssp(nc)%top => top
!#ifdef PARALLEL
!        end if
!#endif       
!    end do
!
!    ntotptssp = 0
!    do nb=1,nblocks
!        top => mb_topsp(nb)
!        ntotptssp = ntotptssp + product(top%nijk(:)) 
!    end do    
!
!end subroutine set_block_coms_cc

subroutine run_on_blkcoms(subprocess)
    use mod_kndconsts, only : kind_int
    use mod_fieldvars, only : nblkcoms,blkcoms
    implicit none
    external :: subprocess
    integer(kind_int) :: nc,nb

    do nc=1,nblkcoms
        nb = blkcoms(nc)%nb

        call subprocess(nb)
    end do

end subroutine run_on_blkcoms

subroutine mb_var_create(mb_var,nst,ned,ngh)
    use mod_kndconsts, only : kind_int
    use mod_constants, only : fldtype_r3d
    use mod_datatypes, only : var_block_t
    use mod_fieldvars, only : nblocks,nblkcoms,blkcoms
    implicit none
    type(var_block_t), pointer, intent(inout) :: mb_var(:)
    integer(kind_int)         , intent(in)    :: nst,ned,ngh
    integer(kind_int) :: nc,nb,ierr
    integer(kind_int) :: st(3),ed(3)
    
    allocate(mb_var(nblocks), stat=ierr)
    
    do nc=1,nblkcoms
        nb = blkcoms(nc)%nb

        st(:) = 1 - ngh
        ed(:) = blkcoms(nc)%top%nijk(:) + ngh
        
        call var_array_create(mb_var(nb),nst,ned, &
                              fldtype_r3d, &
                              st(1),ed(1), &
                              st(2),ed(2), &
                              st(3),ed(3))
    end do

end subroutine mb_var_create

subroutine mb_var_dg_create(mb_var,nst,ned,ngh)
    use mod_kndconsts, only : kind_int
    use mod_constants, only : fldtype_r3d
    use mod_datatypes, only : var_block_t
    use mod_fieldvars, only : nblocks,nblkcoms,blkcoms
    implicit none
    type(var_block_t), pointer, intent(inout) :: mb_var(:)
    integer(kind_int)         , intent(in)    :: nst,ned,ngh
    integer(kind_int) :: nc,nb,ierr
    integer(kind_int) :: st(3),ed(3)
    
    allocate(mb_var(nblocks), stat=ierr)
    
    do nc=1,nblkcoms
        nb = blkcoms(nc)%nb

        st(:) = 1 - ngh*2
        ed(:) = blkcoms(nc)%top%nijk(:)*2 -1 + ngh*2
        
        call var_array_create(mb_var(nb),nst,ned, &
                              fldtype_r3d, &
                              st(1),ed(1), &
                              st(2),ed(2), &
                              st(3),ed(3))
    end do

end subroutine mb_var_dg_create

subroutine mb_var_cc_create(mb_var,nst,ned,ngh)
    use mod_kndconsts, only : kind_int
    use mod_constants, only : fldtype_r3d
    use mod_datatypes, only : var_block_t
    use mod_fieldvars, only : nblocks,nblkcoms,blkcomscc
    implicit none
    type(var_block_t), pointer, intent(inout) :: mb_var(:)
    integer(kind_int)         , intent(in)    :: nst,ned,ngh
    integer(kind_int) :: nc,nb,ierr
    integer(kind_int) :: st(3),ed(3)
    
    allocate(mb_var(nblocks), stat=ierr)
    
    do nc=1,nblkcoms
        nb = blkcomscc(nc)%nb

        st(:) = 1 - ngh
        ed(:) = blkcomscc(nc)%top%nijk(:) + ngh
        
        call var_array_create(mb_var(nb),nst,ned, &
                              fldtype_r3d, &
                              st(1),ed(1), &
                              st(2),ed(2), &
                              st(3),ed(3))
    end do

end subroutine mb_var_cc_create    

subroutine mb_var_create_sp(mb_var,nst,ned,ngh)
    use mod_kndconsts, only : kind_int
    use mod_constants, only : fldtype_r3d
    use mod_datatypes, only : var_block_t
    use mod_fieldvars, only : nblocks,nblkcoms,blkcomssp
    implicit none
    type(var_block_t), pointer, intent(inout) :: mb_var(:)
    integer(kind_int)         , intent(in)    :: nst,ned,ngh
    integer(kind_int) :: nc,nb,ierr
    integer(kind_int) :: st(3),ed(3)
    
    allocate(mb_var(nblocks), stat=ierr)
    
    do nc=1,nblkcoms
        nb = blkcomssp(nc)%nb

        st(:) = 1 - ngh
        ed(:) = blkcomssp(nc)%top%nijk(:) + ngh
        
        call var_array_create(mb_var(nb),nst,ned, &
                              fldtype_r3d, &
                              st(1),ed(1), &
                              st(2),ed(2), &
                              st(3),ed(3))
    end do

end subroutine mb_var_create_sp

subroutine mb_var_delete(mb_var)
    use mod_kndconsts, only : kind_int
    use mod_datatypes, only : var_block_t
    use mod_fieldvars, only : nblocks,nblkcoms,blkcoms
    implicit none
    type(var_block_t), pointer, intent(inout) :: mb_var(:)
    integer(kind_int) :: nc,nb,ierr
    
    do nc=1,nblkcoms
        nb = blkcoms(nc)%nb

        call var_array_delete(mb_var(nb))
    end do
    
    deallocate(mb_var, stat=ierr)

end subroutine mb_var_delete

subroutine mb_var_pointer_create(mb_var,nst,ned)
    use mod_kndconsts, only : kind_int
    use mod_datatypes, only : var_block_t
    use mod_fieldvars, only : nblocks,nblkcoms,blkcoms
    implicit none
    type(var_block_t), pointer, intent(inout) :: mb_var(:)
    integer(kind_int)         , intent(in)    :: nst,ned
    integer(kind_int) :: nc,nb,ierr
    
    allocate(mb_var(nblocks), stat=ierr)
    
    do nc=1,nblkcoms
        nb = blkcoms(nc)%nb

        allocate(mb_var(nb)%fld(nst:ned), stat=ierr)
    end do

end subroutine mb_var_pointer_create

subroutine mb_var_pointer_delete(mb_var)
    use mod_kndconsts, only : kind_int
    use mod_datatypes, only : var_block_t
    use mod_fieldvars, only : nblocks,nblkcoms,blkcoms
    implicit none
    type(var_block_t), pointer, intent(inout) :: mb_var(:)
    integer(kind_int) :: nc,nb,ierr
    
    do nc=1,nblkcoms
        nb = blkcoms(nc)%nb

        deallocate(mb_var(nb)%fld, stat=ierr)
    end do

    deallocate(mb_var, stat=ierr)

end subroutine mb_var_pointer_delete

subroutine mb_var_pointer_assign(mb_vars,nvs_st,nvs_ed, &
                                 mb_vard,nvd_st,nvd_ed)
    use mod_kndconsts, only : kind_int
    use mod_datatypes, only : var_block_t
    use mod_fieldvars, only : nblocks,nblkcoms,blkcoms
    implicit none
    type(var_block_t), pointer, intent(inout) :: mb_vars(:)
    integer(kind_int)         , intent(in)    :: nvs_st,nvs_ed
    type(var_block_t), pointer, intent(inout) :: mb_vard(:)
    integer(kind_int)         , intent(in)    :: nvd_st,nvd_ed
    integer(kind_int) :: nc,nb,m,ierr

    ierr = (nvd_ed-nvd_st) - (nvs_ed-nvs_st)
    call error_check(ierr, &
                     "The size of array isn't equal in subroutine mb_var_pointer_assign")
    
    do nc=1,nblkcoms
        nb = blkcoms(nc)%nb

        do m=0,nvs_ed-nvs_st
            mb_vars(nb)%fld(nvs_st+m)%r3d => mb_vard(nb)%fld(nvd_st+m)%r3d
        end do
    end do

end subroutine mb_var_pointer_assign

