
subroutine input_par
    use mod_constants, only : io_unit_par
    use mod_variables, only : parfile
    implicit none

    call openfile(io_unit_par,parfile,"old","formatted","sequential")
    call input_par_primary(io_unit_par)
    call closefile(io_unit_par,parfile)

end subroutine input_par

subroutine input_top
    use mod_constants, only : io_unit_top
    use mod_variables, only : topfile,splfile
    implicit none

    call openfile(io_unit_top,topfile,"old","formatted","sequential")
    call input_top_gridgen(io_unit_top)
    call closefile(io_unit_top,topfile)

    call init_top_split
    call run_on_master(input_split)

    contains

    subroutine input_split
        implicit none
        logical :: exists
        splfile = trim(topfile)//".splt"
        inquire(file=splfile,exist=exists)
        if (exists) then
            call openfile(io_unit_top,splfile,"old","formatted","sequential")
            call input_top_split(io_unit_top)
            call closefile(io_unit_top,splfile)
        else
            call set_top_split
        end if
    end subroutine input_split
end subroutine input_top

subroutine input_top_cc
    use mod_constants, only : io_unit_top
    use mod_variables, only : topfile,splfile
    implicit none

    call openfile(io_unit_top,topfile,"old","formatted","sequential")
    call input_top_gridgen_cc(io_unit_top)
    call closefile(io_unit_top,topfile)

    call init_top_split_cc
    call run_on_master(input_split_cc)

    contains

    subroutine input_split_cc
        implicit none
        logical :: exists
        splfile = trim(topfile)//".splt"
        inquire(file=splfile,exist=exists)
        if (exists) then
            call openfile(io_unit_top,splfile,"old","formatted","sequential")
            call input_top_split_cc(io_unit_top)
            call closefile(io_unit_top,splfile)
        else
            call set_top_split_cc
        end if
    end subroutine input_split_cc
end subroutine input_top_cc

subroutine input_grd
    use mod_kndconsts, only : kind_int
    use mod_constants, only : io_unit_grd
    use mod_variables, only : ncutpol,ncelcet
    use mod_variables, only : grdfile,nghnode
    use mod_fieldvars, only : nblocks,mb_xyz
    use mod_interface, only : mb_var_create
    implicit none

    if(ncutpol == 5) then
        call mb_var_create(mb_xyz,1,3,nghnode*3)
    else
        call mb_var_create(mb_xyz,1,3,nghnode)
    end if

    call openfile(io_unit_grd,grdfile,"old","unformatted","stream")
    call input_grd_plot3d(io_unit_grd)
    call closefile(io_unit_grd,grdfile)

end subroutine input_grd

subroutine input_grd_cc
    use mod_kndconsts, only : kind_int
    use mod_constants, only : io_unit_grd
    use mod_variables, only : ncutpol,ncelcet
    use mod_variables, only : grdfile,nghnode
    use mod_fieldvars, only : nblocks,mb_xyzcc
    use mod_interface, only : mb_var_create
    implicit none

    if(ncutpol == 5) then
        call msg_seq_and_master("Patched grid is not supported in dual-mesh")
    else
        call mb_var_create(mb_xyzcc,1,3,nghnode)
    end if

    call openfile(io_unit_grd,grdfile,"old","unformatted","stream")
    call input_grd_plot3d_cc(io_unit_grd)
    call closefile(io_unit_grd,grdfile)
    
    call reset_coordinate 

end subroutine input_grd_cc    
    
subroutine reset_coordinate
    use mod_kndconsts, only : kind_int
    use mod_datatypes, only : var_block_t,fld_array_t
    use mod_variables, only : nghnode
    use mod_fieldvars, only : nblocks,nblkcoms,blkcoms,blkcomscc,blkcomssp
    use mod_fieldvars, only : mb_xyzcc,mb_xyzsp
    use mod_interface, only : mb_var_create,mb_var_cc_create,mb_var_delete
    use mod_interface, only : mb_var_create_sp,mb_var_dg_create
    implicit none
    integer(kind_int) :: nc,nb,i,j,k,icc,jcc,kcc
    integer(kind_int) :: st(3),ed(3)    
    type(var_block_t), pointer :: mb_xyz_bak(:),mb_xyzdg(:)
    type(fld_array_t), pointer :: xyz(:),xyzcc(:),xyzsp(:)
    type(fld_array_t), pointer :: xyz_bak(:)
    
    call mb_var_create(mb_xyz_bak,1,3,nghnode)
    do nc=1,nblkcoms
        nb = blkcoms(nc)%nb
        xyz => mb_xyzcc(nb)%fld
        xyz_bak => mb_xyz_bak(nb)%fld
        
        st(:) = 1
        ed(:) = blkcoms(nc)%top%nijk(:)
        
        do i = st(1),ed(1)
            do j = st(2),ed(2)
                do k = st(3),ed(3)
                    xyz_bak(1)%r3d(i,j,k) = xyz(1)%r3d(i,j,k)
                    xyz_bak(2)%r3d(i,j,k) = xyz(2)%r3d(i,j,k)
                    xyz_bak(3)%r3d(i,j,k) = xyz(3)%r3d(i,j,k)                    
                end do
            end do
        end do
    end do

    call mb_var_dg_create(mb_xyzdg,1,3,nghnode)
    call mb_var_create_sp(mb_xyzsp,1,3,nghnode)

    do nc=1,nblkcoms
        nb = blkcoms(nc)%nb
        
        xyz => mb_xyzdg(nb)%fld
        xyz_bak => mb_xyz_bak(nb)%fld

        st(:) = 1
        ed(:) = blkcoms(nc)%top%nijk(:)*2
        
        do i = st(1),ed(1)-1,2
            do j = st(2),ed(2)-1,2
                do k = st(3),ed(3)-1,2
                    icc = (i+1)/2
                    jcc = (j+1)/2
                    kcc = (k+1)/2
                    xyz(1)%r3d(i,j,k) = xyz_bak(1)%r3d(icc,jcc,kcc)
                    xyz(2)%r3d(i,j,k) = xyz_bak(2)%r3d(icc,jcc,kcc)
                    xyz(3)%r3d(i,j,k) = xyz_bak(3)%r3d(icc,jcc,kcc)
                end do
            enddo
        end do
        
        do i = 2,ed(1)-2,2
            do j = st(2),ed(2)-1,2
                do k = st(3),ed(3)-1,2
                    xyz(1)%r3d(i,j,k) = 0.5*(xyz(1)%r3d(i-1,j,k) + xyz(1)%r3d(i+1,j,k))
                    xyz(2)%r3d(i,j,k) = 0.5*(xyz(2)%r3d(i-1,j,k) + xyz(2)%r3d(i+1,j,k))
                    xyz(3)%r3d(i,j,k) = 0.5*(xyz(3)%r3d(i-1,j,k) + xyz(3)%r3d(i+1,j,k))
                end do
            enddo
        end do 
        
        do i = st(1),ed(1)-1,2
            do j = 2,ed(2)-2,2
                do k = st(3),ed(3)-1,2
                    xyz(1)%r3d(i,j,k) = 0.5*(xyz(1)%r3d(i,j-1,k) + xyz(1)%r3d(i,j+1,k))
                    xyz(2)%r3d(i,j,k) = 0.5*(xyz(2)%r3d(i,j-1,k) + xyz(2)%r3d(i,j+1,k))
                    xyz(3)%r3d(i,j,k) = 0.5*(xyz(3)%r3d(i,j-1,k) + xyz(3)%r3d(i,j+1,k))
                end do
            enddo
        end do
        
        do i = st(1),ed(1)-1,2
            do j = st(2),ed(2)-1,2
                do k = 2,ed(3)-2,2
                    xyz(1)%r3d(i,j,k) = 0.5*(xyz(1)%r3d(i,j,k-1) + xyz(1)%r3d(i,j,k+1))
                    xyz(2)%r3d(i,j,k) = 0.5*(xyz(2)%r3d(i,j,k-1) + xyz(2)%r3d(i,j,k+1))
                    xyz(3)%r3d(i,j,k) = 0.5*(xyz(3)%r3d(i,j,k-1) + xyz(3)%r3d(i,j,k+1))
                end do
            enddo
        end do
        
        do i = st(1),ed(1)-1,2
            do j = 2,ed(2)-2,2
                do k = 2,ed(3)-2,2
                    xyz(1)%r3d(i,j,k) = 0.5*(xyz(1)%r3d(i,j,k-1) + xyz(1)%r3d(i,j,k+1))
                    xyz(2)%r3d(i,j,k) = 0.5*(xyz(2)%r3d(i,j,k-1) + xyz(2)%r3d(i,j,k+1))
                    xyz(3)%r3d(i,j,k) = 0.5*(xyz(3)%r3d(i,j,k-1) + xyz(3)%r3d(i,j,k+1))
                end do
            enddo
        end do   
        
        do i = 2,ed(1)-2,2
            do j = st(2),ed(2)-1,2
                do k = 2,ed(3)-2,2
                    xyz(1)%r3d(i,j,k) = 0.5*(xyz(1)%r3d(i,j,k-1) + xyz(1)%r3d(i,j,k+1))
                    xyz(2)%r3d(i,j,k) = 0.5*(xyz(2)%r3d(i,j,k-1) + xyz(2)%r3d(i,j,k+1))
                    xyz(3)%r3d(i,j,k) = 0.5*(xyz(3)%r3d(i,j,k-1) + xyz(3)%r3d(i,j,k+1))
                end do
            enddo
        end do  
        
        do i = 2,ed(1)-2,2
            do j = 2,ed(2)-2,2
                do k = st(3),ed(3)-1,2
                    xyz(1)%r3d(i,j,k) = 0.5*(xyz(1)%r3d(i,j-1,k) + xyz(1)%r3d(i,j+1,k))
                    xyz(2)%r3d(i,j,k) = 0.5*(xyz(2)%r3d(i,j-1,k) + xyz(2)%r3d(i,j+1,k))
                    xyz(3)%r3d(i,j,k) = 0.5*(xyz(3)%r3d(i,j-1,k) + xyz(3)%r3d(i,j+1,k))
                end do
            enddo
        end do 
        
        do i = 2,ed(1)-2,2
            do j = 2,ed(2)-2,2
                do k = 2,ed(3)-2,2
                    xyz(1)%r3d(i,j,k) = 0.5*(xyz(1)%r3d(i,j,k-1) + xyz(1)%r3d(i,j,k+1))
                    xyz(2)%r3d(i,j,k) = 0.5*(xyz(2)%r3d(i,j,k-1) + xyz(2)%r3d(i,j,k+1))
                    xyz(3)%r3d(i,j,k) = 0.5*(xyz(3)%r3d(i,j,k-1) + xyz(3)%r3d(i,j,k+1))
                end do
            enddo
        end do        
            
    end do
    
    do nc=1,nblkcoms
        nb = blkcomscc(nc)%nb
        
        xyz => mb_xyzdg(nb)%fld
        xyzcc => mb_xyzcc(nb)%fld
        xyzsp => mb_xyzsp(nb)%fld

        st(:) = 1
        ed(:) = blkcomscc(nc)%top%nijk(:)
        
        do i = st(1),ed(1)
            do j = st(2),ed(2)
                do k = st(3),ed(3)
                    icc = 2*i-1
                    jcc = 2*j-1
                    kcc = 2*k-1
                    xyzcc(1)%r3d(i,j,k) = xyz(1)%r3d(icc,jcc,kcc)
                    xyzcc(2)%r3d(i,j,k) = xyz(2)%r3d(icc,jcc,kcc)
                    xyzcc(3)%r3d(i,j,k) = xyz(3)%r3d(icc,jcc,kcc)
                end do
            enddo
        end do
        
        st(:) = 1
        ed(:) = blkcomssp(nc)%top%nijk(:)
        
        do i = st(1),ed(1)
            do j = st(2),ed(2)
                do k = st(3),ed(3)
                    icc = 2*i
                    jcc = 2*j
                    kcc = 2*k
                    xyzsp(1)%r3d(i,j,k) = xyz(1)%r3d(icc,jcc,kcc)
                    xyzsp(2)%r3d(i,j,k) = xyz(2)%r3d(icc,jcc,kcc)
                    xyzsp(3)%r3d(i,j,k) = xyz(3)%r3d(icc,jcc,kcc)
                end do
            enddo
        end do
            
    end do    

    call mb_var_delete(mb_xyzdg)
    call mb_var_delete(mb_xyz_bak)
    
end subroutine reset_coordinate      


subroutine input_par_primary(io_unit)
    use mod_kndconsts, only : kind_int
    use mod_variables, only : general,inflow,reference
    use mod_variables, only : filename,initconst,control
    use mod_variables, only : physical,turbulent,timestep
    use mod_variables, only : method,technic,gasmodel,monitor
    implicit none
    integer(kind_int), intent(in) :: io_unit
    integer(kind_int) :: ierr

    read(io_unit,nml=general  , iostat=ierr)
    read(io_unit,nml=inflow   , iostat=ierr)
    read(io_unit,nml=reference, iostat=ierr)
    read(io_unit,nml=filename , iostat=ierr)
    read(io_unit,nml=initconst, iostat=ierr)
    read(io_unit,nml=control  , iostat=ierr)

    read(io_unit,nml=physical , iostat=ierr)
    read(io_unit,nml=turbulent, iostat=ierr)
    read(io_unit,nml=timestep , iostat=ierr)
    read(io_unit,nml=method   , iostat=ierr)
    read(io_unit,nml=technic  , iostat=ierr)

    read(io_unit,nml=gasmodel , iostat=ierr)
    read(io_unit,nml=monitor  , iostat=ierr)

end subroutine input_par_primary

subroutine input_top_gridgen(io_unit)
    use mod_kndconsts, only : kind_int
    use mod_constants, only : bc_cut1to1,bc_cut1to1A,bc_cut1to1_split,bc_cut1to1_splitA
    use mod_constants, only : subc_none,subc_cut_original,bc_wallA,bc_wall,bc_patched,bc_patchedA
    use mod_constants, only : subc_cut_spliting,subc_cut_patched
    use mod_constants, only : bc_pole,bc_pole_i,bc_pole_j,bc_pole_k
    use mod_constants, only : subc_pole_i,subc_pole_j,subc_pole_k
    use mod_datatypes, only : bc_region_t,top_block_t
    use mod_fieldvars, only : nblocks,nprocs,mb_top
    implicit none
    integer(kind_int), intent(in) :: io_unit

    ! local storage
    integer(kind_int)          :: solverid
    type(top_block_t), pointer :: top
    type(bc_region_t), pointer :: reg
    integer(kind_int)          :: nb,nr,m,bctype,subtype,ierr

    read(io_unit, *, iostat=ierr) solverid

    read(io_unit, *, iostat=ierr) nblocks
    allocate(mb_top(nblocks), stat=ierr)

    do nb=1,nblocks
        top => mb_top(nb)

        read(io_unit, *, iostat=ierr) (top%nijk(m),m=1,3)
        read(io_unit, *, iostat=ierr)  top%name

        read(io_unit, *, iostat=ierr)  top%nregions
        allocate(top%bcs(top%nregions), stat=ierr)
        do nr=1,top%nregions
            reg => top%bcs(nr)
            read(io_unit, *, iostat=ierr) (reg%s_st(m), &
                                           reg%s_ed(m),m=1,3), &
                                           reg%bctype

            bctype = reg%bctype
            select case(bctype)
            case(bc_wallA)
                reg%bctype    = bc_wall
                reg%subtype_A = 1            
            case(bc_cut1to1)
                reg%subtype   = subc_cut_original
                reg%subtype_A = 0
            case(bc_cut1to1A)
                reg%bctype    = bc_cut1to1
                reg%subtype   = subc_cut_original
                reg%subtype_A = 1
            case(bc_cut1to1_split)   !! 剖分出来的对接面
                reg%bctype    = bc_cut1to1
                reg%subtype   = subc_cut_spliting
                reg%subtype_A = 0
            case(bc_cut1to1_splitA)  !! 剖分出来的对接面
                reg%bctype    = bc_cut1to1
                reg%subtype   = subc_cut_patched
                reg%subtype_A = 1
            case(bc_patched)
                reg%bctype    = bc_cut1to1
                reg%subtype   = subc_cut_patched
                reg%subtype_A = 0
            case(bc_patchedA)
                reg%bctype    = bc_cut1to1
                reg%subtype   = subc_cut_patched
                reg%subtype_A = 1
            case(bc_pole_i)
                reg%bctype    = bc_pole
                reg%subtype   = subc_pole_i
                reg%subtype_A = 0
            case(bc_pole_j)
                reg%bctype    = bc_pole
                reg%subtype   = subc_pole_j
                reg%subtype_A = 0
            case(bc_pole_k)
                reg%bctype    = bc_pole
                reg%subtype   = subc_pole_k
                reg%subtype_A = 0
            case default
                reg%subtype   = subc_none
                reg%subtype_A = 0
            end select

            bctype  = reg%bctype
            subtype = reg%subtype
            if (bctype == bc_cut1to1 .and. subtype /= subc_cut_patched) then
                read(io_unit, *, iostat=ierr) (reg%t_st(m), &
                                               reg%t_ed(m),m=1,3), &
                                               reg%nbt
            end if
       end do
    end do

    call error_check(ierr,"BC file is incorrect!")

end subroutine input_top_gridgen

subroutine input_top_gridgen_cc(io_unit)
    use mod_kndconsts, only : kind_int
    use mod_constants, only : bc_cut1to1,bc_cut1to1A,bc_cut1to1_split,bc_cut1to1_splitA
    use mod_constants, only : subc_none,subc_cut_original,bc_wallA,bc_wall,bc_patched,bc_patchedA
    use mod_constants, only : subc_cut_spliting,subc_cut_patched
    use mod_constants, only : bc_pole,bc_pole_i,bc_pole_j,bc_pole_k
    use mod_constants, only : subc_pole_i,subc_pole_j,subc_pole_k
    use mod_datatypes, only : bc_region_t,top_block_t
    use mod_fieldvars, only : nblocks,nprocs,mb_top,mb_topc,mb_topsp
    implicit none
    integer(kind_int), intent(in) :: io_unit

    ! local storage
    integer(kind_int)          :: solverid
    type(top_block_t), pointer :: top
    type(bc_region_t), pointer :: reg
    type(top_block_t), pointer :: topc
    type(bc_region_t), pointer :: regc
    type(top_block_t), pointer :: topsp
    type(bc_region_t), pointer :: regsp    
    integer(kind_int)          :: nb,nr,m,bctype,subtype,ierr

    read(io_unit, *, iostat=ierr) solverid

    read(io_unit, *, iostat=ierr) nblocks
    allocate(mb_top(nblocks), stat=ierr)

    do nb=1,nblocks
        top => mb_top(nb)

        read(io_unit, *, iostat=ierr) (top%nijk(m),m=1,3)
        read(io_unit, *, iostat=ierr)  top%name

        read(io_unit, *, iostat=ierr)  top%nregions
        allocate(top%bcs(top%nregions), stat=ierr)
        do nr=1,top%nregions
            reg => top%bcs(nr)
            read(io_unit, *, iostat=ierr) (reg%s_st(m), &
                                           reg%s_ed(m),m=1,3), &
                                           reg%bctype

            bctype = reg%bctype
            select case(bctype)
            case(bc_wallA)
                reg%bctype    = bc_wall
                reg%subtype_A = 1            
            case(bc_cut1to1)
                reg%subtype   = subc_cut_original
                reg%subtype_A = 0
            case(bc_cut1to1A)
                reg%bctype    = bc_cut1to1
                reg%subtype   = subc_cut_original
                reg%subtype_A = 1
            case(bc_cut1to1_split)   !! 剖分出来的对接面
                reg%bctype    = bc_cut1to1
                reg%subtype   = subc_cut_spliting
                reg%subtype_A = 0
            case(bc_cut1to1_splitA)  !! 剖分出来的对接面
                reg%bctype    = bc_cut1to1
                reg%subtype   = subc_cut_patched
                reg%subtype_A = 1
            case(bc_patched)
                reg%bctype    = bc_cut1to1
                reg%subtype   = subc_cut_patched
                reg%subtype_A = 0
            case(bc_patchedA)
                reg%bctype    = bc_cut1to1
                reg%subtype   = subc_cut_patched
                reg%subtype_A = 1
            case(bc_pole_i)
                reg%bctype    = bc_pole
                reg%subtype   = subc_pole_i
                reg%subtype_A = 0
            case(bc_pole_j)
                reg%bctype    = bc_pole
                reg%subtype   = subc_pole_j
                reg%subtype_A = 0
            case(bc_pole_k)
                reg%bctype    = bc_pole
                reg%subtype   = subc_pole_k
                reg%subtype_A = 0
            case default
                reg%subtype   = subc_none
                reg%subtype_A = 0
            end select

            bctype  = reg%bctype
            subtype = reg%subtype
            if (bctype == bc_cut1to1 .and. subtype /= subc_cut_patched) then
                read(io_unit, *, iostat=ierr) (reg%t_st(m), &
                                               reg%t_ed(m),m=1,3), &
                                               reg%nbt
            end if
       end do
    end do

    allocate(mb_topc(nblocks), stat=ierr)
    allocate(mb_topsp(nblocks), stat=ierr)
    
    do nb=1,nblocks
        top => mb_top(nb)
        topc => mb_topc(nb)
        topsp => mb_topsp(nb)

        topc%nijk(:) = top%nijk(:)
        topc%name = top%name
        topsp%nijk(:) = top%nijk(:) - 1
        topsp%name = top%name

        topc%nregions = top%nregions
        topsp%nregions = top%nregions
        allocate(topc%bcs(topc%nregions), stat=ierr)
        allocate(topsp%bcs(topsp%nregions), stat=ierr)
        do nr=1,topc%nregions
            reg => top%bcs(nr)
            regc => topc%bcs(nr)

            do m=1,3
                regc%s_st(m) = reg%s_st(m)
                regc%s_ed(m) = reg%s_ed(m)     
            end do
            regc%bctype = reg%bctype
            bctype = regc%bctype
            select case(bctype)
            case(bc_wallA)
                regc%bctype    = bc_wall
                regc%subtype_A = 1            
            case(bc_cut1to1)
                regc%subtype   = subc_cut_original
                regc%subtype_A = 0
            case(bc_cut1to1A)
                regc%bctype    = bc_cut1to1
                regc%subtype   = subc_cut_original
                regc%subtype_A = 1
            case(bc_cut1to1_split)   !! 剖分出来的对接面
                regc%bctype    = bc_cut1to1
                regc%subtype   = subc_cut_spliting
                regc%subtype_A = 0
            case(bc_cut1to1_splitA)  !! 剖分出来的对接面
                regc%bctype    = bc_cut1to1
                regc%subtype   = subc_cut_patched
                regc%subtype_A = 1
            case(bc_pole_i)
                regc%bctype    = bc_pole
                regc%subtype   = subc_pole_i
                regc%subtype_A = 0
            case(bc_pole_j)
                regc%bctype    = bc_pole
                regc%subtype   = subc_pole_j
                regc%subtype_A = 0
            case(bc_pole_k)
                regc%bctype    = bc_pole
                regc%subtype   = subc_pole_k
                regc%subtype_A = 0
            case default
                regc%subtype   = subc_none
                regc%subtype_A = 0
            end select

            bctype  = regc%bctype
            if (bctype == bc_cut1to1) then
                do m=1,3
                    regc%t_st(m) = reg%t_st(m)
                    regc%t_ed(m) = reg%t_ed(m)
                end do
                regc%nbt = reg%nbt
            end if
        end do
        
        do nr=1,topsp%nregions
            reg => top%bcs(nr)
            regsp => topsp%bcs(nr)

            do m=1,3
                if(abs(reg%s_st(m)) > abs(reg%s_ed(m))) then
                    regsp%s_st(m) = sign(abs(reg%s_st(m)) - 1,reg%s_st(m))
                    regsp%s_ed(m) = reg%s_ed(m)
                elseif(abs(reg%s_st(m)) == abs(reg%s_ed(m)) .and. abs(reg%s_st(m)) == 1) then
                    regsp%s_st(m) = reg%s_st(m)
                    regsp%s_ed(m) = reg%s_ed(m)
                elseif(abs(reg%s_st(m)) == abs(reg%s_ed(m)) .and. abs(reg%s_st(m)) /= 1) then
                    regsp%s_st(m) = sign(abs(reg%s_st(m)) - 1,reg%s_st(m))
                    regsp%s_ed(m) = sign(abs(reg%s_ed(m)) - 1,reg%s_ed(m))
                else
                    regsp%s_st(m) = reg%s_st(m)
                    regsp%s_ed(m) = sign(abs(reg%s_ed(m)) - 1,reg%s_ed(m))                        
                end if    
            end do
            regsp%bctype = reg%bctype
            bctype = regsp%bctype
            select case(bctype)
            case(bc_wallA)
                regsp%bctype    = bc_wall
                regsp%subtype_A = 1            
            case(bc_cut1to1)
                regsp%subtype   = subc_cut_original
                regsp%subtype_A = 0
            case(bc_cut1to1A)
                regsp%bctype    = bc_cut1to1
                regsp%subtype   = subc_cut_original
                regsp%subtype_A = 1
            case(bc_cut1to1_split)   !! 剖分出来的对接面
                regsp%bctype    = bc_cut1to1
                regsp%subtype   = subc_cut_spliting
                regsp%subtype_A = 0
            case(bc_cut1to1_splitA)  !! 剖分出来的对接面
                regsp%bctype    = bc_cut1to1
                regsp%subtype   = subc_cut_patched
                regsp%subtype_A = 1
            case(bc_pole_i)
                regsp%bctype    = bc_pole
                regsp%subtype   = subc_pole_i
                regsp%subtype_A = 0
            case(bc_pole_j)
                regsp%bctype    = bc_pole
                regsp%subtype   = subc_pole_j
                regsp%subtype_A = 0
            case(bc_pole_k)
                regsp%bctype    = bc_pole
                regsp%subtype   = subc_pole_k
                regsp%subtype_A = 0
            case default
                regsp%subtype   = subc_none
                regsp%subtype_A = 0
            end select

            bctype  = regsp%bctype
            if (bctype == bc_cut1to1) then
                do m=1,3
                    if(abs(reg%t_st(m)) > abs(reg%t_ed(m))) then
                        regsp%t_st(m) = sign(abs(reg%t_st(m)) - 1,reg%t_st(m))
                        regsp%t_ed(m) = reg%t_ed(m)
                    elseif(abs(reg%t_st(m)) == abs(reg%t_ed(m)) .and. abs(reg%t_st(m)) == 1) then
                        regsp%t_st(m) = reg%t_st(m)
                        regsp%t_ed(m) = reg%t_ed(m)
                    elseif(abs(reg%t_st(m)) == abs(reg%t_ed(m)) .and. abs(reg%t_st(m)) /= 1) then
                        regsp%t_st(m) = sign(abs(reg%t_st(m)) - 1,reg%t_st(m))
                        regsp%t_ed(m) = sign(abs(reg%t_ed(m)) - 1,reg%t_ed(m))    
                    else    
                        regsp%t_st(m) = reg%t_st(m)
                        regsp%t_ed(m) = sign(abs(reg%t_ed(m)) - 1,reg%t_ed(m))                        
                    end if    
                end do
                regsp%nbt = reg%nbt
            end if
       end do        
    end do   

    call error_check(ierr,"BC file is incorrect!")

end subroutine input_top_gridgen_cc

subroutine init_top_split
    use mod_kndconsts, only : kind_int
    use mod_datatypes, only : top_block_t
    use mod_fieldvars, only : nblocks,nprocs,mb_top
    implicit none

    ! local storage
    type(top_block_t), pointer :: top
    integer(kind_int)          :: nb

    nprocs = 1
    do nb=1,nblocks
        top => mb_top(nb)

        top%pid = 0
        top%pnb = nb
        top%pst(:) = 1
        top%ped(:) = top%nijk(:)
        top%nupd = 1
    end do

end subroutine init_top_split

subroutine init_top_split_cc
    use mod_kndconsts, only : kind_int
    use mod_datatypes, only : top_block_t
    use mod_fieldvars, only : nblocks,nprocs,mb_top,mb_topc,mb_topsp  
    implicit none

    ! local storage
    type(top_block_t), pointer :: top
    integer(kind_int)          :: nb

    nprocs = 1
    do nb=1,nblocks
        top => mb_top(nb)

        top%pid = 0
        top%pnb = nb
        top%pst(:) = 1
        top%ped(:) = top%nijk(:)
        top%nupd = 1
    end do
   
    do nb=1,nblocks
        top => mb_topc(nb)

        top%pid = 0
        top%pnb = nb
        top%pst(:) = 1
        top%ped(:) = top%nijk(:)
        top%nupd = 1
    end do  
    
    do nb=1,nblocks
        top => mb_topsp(nb)

        top%pid = 0
        top%pnb = nb
        top%pst(:) = 1
        top%ped(:) = top%nijk(:)
        top%nupd = 1
    end do    

end subroutine init_top_split_cc

subroutine set_top_split
    use mod_kndconsts, only : kind_int
    use mod_datatypes, only : top_block_t
    use mod_fieldvars, only : nblocks,nprocs,mb_top
    implicit none

    ! local storage
    type(top_block_t), pointer :: top
    integer(kind_int)          :: nb

    nprocs = nblocks
    do nb=1,nblocks
        top => mb_top(nb)

        top%pid = nb - 1
        top%pnb = nb
        top%pst(:) = 1
        top%ped(:) = top%nijk(:)
        top%nupd = 1
    end do

end subroutine set_top_split

subroutine set_top_split_cc
    use mod_kndconsts, only : kind_int
    use mod_datatypes, only : top_block_t
    use mod_fieldvars, only : nblocks,nprocs,mb_top,mb_topc,mb_topsp  
    implicit none

    ! local storage
    type(top_block_t), pointer :: top
    integer(kind_int)          :: nb

    nprocs = nblocks
    do nb=1,nblocks
        top => mb_top(nb)

        top%pid = nb - 1
        top%pnb = nb
        top%pst(:) = 1
        top%ped(:) = top%nijk(:)
        top%nupd = 1
    end do

    do nb=1,nblocks
        top => mb_topc(nb)

        top%pid = nb - 1
        top%pnb = nb
        top%pst(:) = 1
        top%ped(:) = top%nijk(:)
        top%nupd = 1
    end do  
    
    do nb=1,nblocks
        top => mb_topsp(nb)

        top%pid = nb - 1
        top%pnb = nb
        top%pst(:) = 1
        top%ped(:) = top%nijk(:)
        top%nupd = 1
    end do    


end subroutine set_top_split_cc

subroutine input_top_split(io_unit)
    use mod_kndconsts, only : kind_int
    use mod_datatypes, only : top_block_t
    use mod_fieldvars, only : nblocks,nprocs,mb_top
    implicit none
    integer(kind_int), intent(in) :: io_unit

    ! local storage
    type(top_block_t), pointer :: top
    integer(kind_int)          :: nblkspart,nbpart
    integer(kind_int)          :: nb,m,ierr

    read(io_unit, *, iostat=ierr) nblkspart,nprocs
    call error_check(abs(ierr)+abs(nblocks-nblkspart),"Split file is incorrect!")
    do nb=1,nblocks
        top => mb_top(nb)
        read(io_unit, *, iostat=ierr) nbpart,top%pid,top%pnb,(top%pst(m),top%ped(m),m=1,3),top%nupd
        call error_check(abs(ierr)+abs(nb-nbpart),"Split file is incorrect!")
    end do

end subroutine input_top_split

subroutine input_top_split_cc(io_unit)
    use mod_kndconsts, only : kind_int
    use mod_datatypes, only : top_block_t
    use mod_fieldvars, only : nblocks,nprocs,mb_top,mb_topc,mb_topsp
    implicit none
    integer(kind_int), intent(in) :: io_unit

    ! local storage
    type(top_block_t), pointer :: top
    type(top_block_t), pointer :: topc
    type(top_block_t), pointer :: topsp
    integer(kind_int)          :: nblkspart,nbpart
    integer(kind_int)          :: nb,m,ierr

    read(io_unit, *, iostat=ierr) nblkspart,nprocs
    call error_check(abs(ierr)+abs(nblocks-nblkspart),"Split file is incorrect!")
    do nb=1,nblocks
        top => mb_top(nb)
        read(io_unit, *, iostat=ierr) nbpart,top%pid,top%pnb,(top%pst(m),top%ped(m),m=1,3),top%nupd
        call error_check(abs(ierr)+abs(nb-nbpart),"Split file is incorrect!")
    end do
   
    do nb=1,nblocks
        topc => mb_topc(nb)
        topsp => mb_topsp(nb)
        top => mb_top(nb)
    
        topc%pid = top%pid
        topc%pnb = top%pnb
        topsp%pid = top%pid
        topsp%pnb = top%pnb        
        do m =1,3
            topc%pst(m) = top%pst(m)
            topc%ped(m) = top%ped(m)
            topsp%pst(m) = top%pst(m)
            topsp%ped(m) = top%ped(m)            
        end do
        topc%nupd = top%nupd
        topsp%nupd = top%nupd
    end do                
end subroutine input_top_split_cc

subroutine input_grd_plot3d(io_unit)
    use mod_kndconsts, only : kind_int
    use mod_constants, only : nplot3d_std
    use mod_fieldvars, only : mb_xyz
    use mod_interface, only : scatter_input_mb_var
    implicit none
    integer(kind_int), intent(in) :: io_unit

    call scatter_input_mb_var(io_unit,mb_xyz,1,3,nplot3d_std)

end subroutine input_grd_plot3d

subroutine input_grd_plot3d_cc(io_unit)
    use mod_kndconsts, only : kind_int
    use mod_constants, only : nplot3d_std
    use mod_fieldvars, only : mb_xyzcc
    use mod_interface, only : scatter_input_mb_var_cc
    implicit none
    integer(kind_int), intent(in) :: io_unit

    call scatter_input_mb_var_cc(io_unit,mb_xyzcc,1,3,nplot3d_std)

end subroutine input_grd_plot3d_cc

subroutine input_sol
    use mod_kndconsts, only : kind_int
    use mod_constants, only : io_unit_fld,io_unit_aux
    use mod_constants, only : io_unit_men,io_unit_rms,io_unit_fmen,nvis_ns_lam
    use mod_variables, only : solfile,nvis,nrestrt,nstepsav,ns_mean,ns_prms,nprms
    use mod_parallels
    implicit none
    integer(kind_int) :: error,ierr,check1,check2

    check1 = 0
    check2 = 0
    ierr   = 0

    call openfile_check(io_unit_fld,solfile,"old","unformatted","stream",ierr)    
    
#ifdef PARALLEL
    call MPI_BCAST(ierr,1,kind_int_mpi,master,MPI_COMM_WORLD,error)
#endif
    
    if(ierr == 0) then
        call input_sol_plot3d(io_unit_fld)
        call closefile(io_unit_fld,solfile)        
    else        
        check1 = check1 + 1       
    end if
    
    call openfile_check(io_unit_aux,trim(solfile)//".aux","old","formatted","sequential",ierr)
    
#ifdef PARALLEL
    call MPI_BCAST(ierr,1,kind_int_mpi,master,MPI_COMM_WORLD,error)
#endif
    
    if(ierr == 0 .and. check1 == 0) then
        call broadcast_input_aux(io_unit_aux)
        call closefile(io_unit_aux,trim(solfile)//".aux")
        
        if (nvis > nvis_ns_lam) then
            nrestrt = 2
        else
            nrestrt = 1
        endif
        call msg_seq_and_master("input the flow field successfully")
    else
        check1 = check1 + 1
    end if
    
    if (nprms > 0) then
    if (nstepsav > ns_mean .and. nstepsav < ns_prms) then
        call openfile_check(io_unit_men,trim(solfile)//".mean","old","unformatted","stream",ierr)
        
#ifdef PARALLEL
        call MPI_BCAST(ierr,1,kind_int_mpi,master,MPI_COMM_WORLD,error)
#endif

        if(ierr == 0 .and. check1 == 0) then
            call input_mean_plot3d(io_unit_men)
            call closefile(io_unit_men,trim(solfile)//".mean")
            call msg_seq_and_master("input the statistical flow field successfully")
        else
            check1 = check1 + 1
        end if
    end if
    
    if (nstepsav == ns_prms) then
        call openfile_check(io_unit_fmen,trim(solfile)//".mean","old","unformatted","stream",ierr)
        
#ifdef PARALLEL
        call MPI_BCAST(ierr,1,kind_int_mpi,master,MPI_COMM_WORLD,error)
#endif

        if(ierr == 0 .and. check1 == 0) then
            call input_mean_plot3d(io_unit_fmen)
            call closefile(io_unit_fmen,trim(solfile)//".mean")
            call msg_seq_and_master("input the statistical flow field successfully")
        else
            check1 = check1 + 1
        end if
    end if
    
    if (nstepsav > ns_prms) then
        call openfile_check(io_unit_fmen,trim(solfile)//".mean","old","unformatted","stream",ierr)
        
#ifdef PARALLEL
        call MPI_BCAST(ierr,1,kind_int_mpi,master,MPI_COMM_WORLD,error)
#endif            
        
        if(ierr == 0 .and. check1 == 0) then
            call input_mean_plot3d(io_unit_fmen)
            call closefile(io_unit_fmen,trim(solfile)//".mean")
        else
            check1 = check1 + 1
        end if
        
        call openfile_check(io_unit_rms,trim(solfile)//".rms","old","unformatted","stream",ierr)
        
#ifdef PARALLEL
        call MPI_BCAST(ierr,1,kind_int_mpi,master,MPI_COMM_WORLD,error)
#endif             
        
        if(ierr == 0 .and. check1 == 0) then
            call input_rms_plot3d(io_unit_rms)
            call closefile(io_unit_rms,trim(solfile)//".rms")
            call msg_seq_and_master("input the statistical flow field successfully")
        else
            check1 = check1 + 1
        end if
    end if
    end if

    if(check1 /= 0) then    
        call input_sol_backup(check2)
    end if   
    
    if(check2 /= 0) then
        call Rest_flow_field
    end if

end subroutine input_sol

subroutine input_sol_backup(check)
    use mod_kndconsts, only : kind_int
    use mod_constants, only : io_unit_fld,io_unit_aux
    use mod_constants, only : io_unit_men,io_unit_rms,io_unit_fmen,nvis_ns_lam
    use mod_variables, only : solfile,nrestrt,nvis,nstepsav,ns_mean,ns_prms,nprms
    use mod_parallels    
    implicit none
    integer(kind_int), intent(out) :: check    
    integer(kind_int) :: ierr,error
    
    check = 0
    ierr  = 0

    call openfile_check(io_unit_fld,trim(solfile)//".bak","old","unformatted","stream",ierr)
    
#ifdef PARALLEL
    call MPI_BCAST(ierr,1,kind_int_mpi,master,MPI_COMM_WORLD,error)
#endif       

    if(ierr == 0) then
        call input_sol_plot3d(io_unit_fld)
        call closefile(io_unit_fld,trim(solfile)//".bak")        
    else       
        check = check + 1
    end if    
    
    call openfile_check(io_unit_aux,trim(solfile)//".aux"//".bak","old","formatted","sequential",ierr)
    
#ifdef PARALLEL
    call MPI_BCAST(ierr,1,kind_int_mpi,master,MPI_COMM_WORLD,error)
#endif     
    
    if(ierr == 0 .and. check == 0) then
        call broadcast_input_aux(io_unit_aux)
        call closefile(io_unit_aux,trim(solfile)//".aux"//".bak")
        
        if (nvis > nvis_ns_lam) then
            nrestrt = 2
        else
            nrestrt = 1
        endif
        call msg_seq_and_master("input the flow field successfully")   
    else
        check = check + 1
    end if
    
    if (nprms > 0) then
    if (nstepsav > ns_mean .and. nstepsav < ns_prms) then
        call openfile_check(io_unit_men,trim(solfile)//".mean"//".bak","old","unformatted","stream",ierr)
        
#ifdef PARALLEL
        call MPI_BCAST(ierr,1,kind_int_mpi,master,MPI_COMM_WORLD,error)
#endif        
        
        if(ierr == 0 .and. check == 0) then
            call input_mean_plot3d(io_unit_men)
            call closefile(io_unit_men,trim(solfile)//".mean"//".bak")
            call msg_seq_and_master("input the statistical flow field successfully")
        else
            check = check + 1
        end if
    end if
    
    if (nstepsav == ns_prms) then
        call openfile_check(io_unit_fmen,trim(solfile)//".mean"//".bak","old","unformatted","stream",ierr)
        
#ifdef PARALLEL
        call MPI_BCAST(ierr,1,kind_int_mpi,master,MPI_COMM_WORLD,error)
#endif        
        
        if(ierr == 0 .and. check == 0) then
            call input_mean_plot3d(io_unit_fmen)
            call closefile(io_unit_fmen,trim(solfile)//".mean"//".bak")
            call msg_seq_and_master("input the statistical flow field successfully")
        else
            check = check + 1
        end if
    end if
    
    if (nstepsav > ns_prms) then
        call openfile_check(io_unit_fmen,trim(solfile)//".mean"//".bak","old","unformatted","stream",ierr)
        
#ifdef PARALLEL
        call MPI_BCAST(ierr,1,kind_int_mpi,master,MPI_COMM_WORLD,error)
#endif        
        
        if(ierr == 0 .and. check == 0) then
            call input_mean_plot3d(io_unit_fmen)
            call closefile(io_unit_fmen,trim(solfile)//".mean"//".bak")
        else
            check = check + 1
        end if
    
        call openfile_check(io_unit_rms,trim(solfile)//".rms"//".bak","old","unformatted","stream",ierr)
        
#ifdef PARALLEL
        call MPI_BCAST(ierr,1,kind_int_mpi,master,MPI_COMM_WORLD,error)
#endif        
        
        if(ierr == 0 .and. check == 0) then
            call input_rms_plot3d(io_unit_rms)
            call closefile(io_unit_rms,trim(solfile)//".rms"//".bak")
            call msg_seq_and_master("input the statistical flow field successfully")
        else
            check = check + 1
        end if
    end if
    end if    
    
end subroutine input_sol_backup

subroutine input_sol_sp
    use mod_kndconsts, only : kind_int
    use mod_constants, only : io_unit_fld,io_unit_aux
    use mod_constants, only : io_unit_men,io_unit_rms,io_unit_fmen,nvis_ns_lam
    use mod_variables, only : solfile,nvis,nrestrt,nstepsav,ns_mean,ns_prms,nprms
    use mod_parallels
    implicit none
    integer(kind_int) :: error,ierr,check1,check2

    check1 = 0
    check2 = 0
    ierr   = 0

    call openfile_check(io_unit_fld,solfile,"old","unformatted","stream",ierr)    
    
#ifdef PARALLEL
    call MPI_BCAST(ierr,1,kind_int_mpi,master,MPI_COMM_WORLD,error)
#endif
    
    if(ierr == 0) then
        call input_sol_plot3d_sp(io_unit_fld)
        call closefile(io_unit_fld,solfile)        
    else        
        check1 = check1 + 1       
    end if
    
    call openfile_check(io_unit_aux,trim(solfile)//".aux","old","formatted","sequential",ierr)
    
#ifdef PARALLEL
    call MPI_BCAST(ierr,1,kind_int_mpi,master,MPI_COMM_WORLD,error)
#endif
    
    if(ierr == 0 .and. check1 == 0) then
        call broadcast_input_aux(io_unit_aux)
        call closefile(io_unit_aux,trim(solfile)//".aux")
        
        if (nvis > nvis_ns_lam) then
            nrestrt = 2
        else
            nrestrt = 1
        endif
        call msg_seq_and_master("input the flow field successfully")
    else
        check1 = check1 + 1
    end if
    
    if (nprms > 0) then
    if (nstepsav > ns_mean .and. nstepsav < ns_prms) then
        call openfile_check(io_unit_men,trim(solfile)//".mean","old","unformatted","stream",ierr)
        
#ifdef PARALLEL
        call MPI_BCAST(ierr,1,kind_int_mpi,master,MPI_COMM_WORLD,error)
#endif

        if(ierr == 0 .and. check1 == 0) then
            call input_mean_plot3d_sp(io_unit_men)
            call closefile(io_unit_men,trim(solfile)//".mean")
            call msg_seq_and_master("input the statistical flow field successfully")
        else
            check1 = check1 + 1
        end if
    end if
    
    if (nstepsav == ns_prms) then
        call openfile_check(io_unit_fmen,trim(solfile)//".mean","old","unformatted","stream",ierr)
        
#ifdef PARALLEL
        call MPI_BCAST(ierr,1,kind_int_mpi,master,MPI_COMM_WORLD,error)
#endif

        if(ierr == 0 .and. check1 == 0) then
            call input_mean_plot3d_sp(io_unit_fmen)
            call closefile(io_unit_fmen,trim(solfile)//".mean")
            call msg_seq_and_master("input the statistical flow field successfully")
        else
            check1 = check1 + 1
        end if
    end if
    
    if (nstepsav > ns_prms) then
        call openfile_check(io_unit_fmen,trim(solfile)//".mean","old","unformatted","stream",ierr)
        
#ifdef PARALLEL
        call MPI_BCAST(ierr,1,kind_int_mpi,master,MPI_COMM_WORLD,error)
#endif            
        
        if(ierr == 0 .and. check1 == 0) then
            call input_mean_plot3d_sp(io_unit_fmen)
            call closefile(io_unit_fmen,trim(solfile)//".mean")
        else
            check1 = check1 + 1
        end if
        
        call openfile_check(io_unit_rms,trim(solfile)//".rms","old","unformatted","stream",ierr)
        
#ifdef PARALLEL
        call MPI_BCAST(ierr,1,kind_int_mpi,master,MPI_COMM_WORLD,error)
#endif             
        
        if(ierr == 0 .and. check1 == 0) then
            call input_rms_plot3d_sp(io_unit_rms)
            call closefile(io_unit_rms,trim(solfile)//".rms")
            call msg_seq_and_master("input the statistical flow field successfully")
        else
            check1 = check1 + 1
        end if
    end if
    end if

    if(check1 /= 0) then    
        call input_sol_backup_sp(check2)
    end if   
    
    if(check2 /= 0) then
        call Rest_flow_field_sp
    end if

end subroutine input_sol_sp

subroutine input_sol_backup_sp(check)
    use mod_kndconsts, only : kind_int
    use mod_constants, only : io_unit_fld,io_unit_aux
    use mod_constants, only : io_unit_men,io_unit_rms,io_unit_fmen,nvis_ns_lam
    use mod_variables, only : solfile,nrestrt,nvis,nstepsav,ns_mean,ns_prms,nprms
    use mod_parallels    
    implicit none
    integer(kind_int), intent(out) :: check    
    integer(kind_int) :: ierr,error
    
    check = 0
    ierr  = 0

    call openfile_check(io_unit_fld,trim(solfile)//".bak","old","unformatted","stream",ierr)
    
#ifdef PARALLEL
    call MPI_BCAST(ierr,1,kind_int_mpi,master,MPI_COMM_WORLD,error)
#endif       

    if(ierr == 0) then
        call input_sol_plot3d_sp(io_unit_fld)
        call closefile(io_unit_fld,trim(solfile)//".bak")        
    else       
        check = check + 1
    end if    
    
    call openfile_check(io_unit_aux,trim(solfile)//".aux"//".bak","old","formatted","sequential",ierr)
    
#ifdef PARALLEL
    call MPI_BCAST(ierr,1,kind_int_mpi,master,MPI_COMM_WORLD,error)
#endif     
    
    if(ierr == 0 .and. check == 0) then
        call broadcast_input_aux(io_unit_aux)
        call closefile(io_unit_aux,trim(solfile)//".aux"//".bak")
        
        if (nvis > nvis_ns_lam) then
            nrestrt = 2
        else
            nrestrt = 1
        endif
        call msg_seq_and_master("input the flow field successfully")   
    else
        check = check + 1
    end if
    
    if (nprms > 0) then
    if (nstepsav > ns_mean .and. nstepsav < ns_prms) then
        call openfile_check(io_unit_men,trim(solfile)//".mean"//".bak","old","unformatted","stream",ierr)
        
#ifdef PARALLEL
        call MPI_BCAST(ierr,1,kind_int_mpi,master,MPI_COMM_WORLD,error)
#endif        
        
        if(ierr == 0 .and. check == 0) then
            call input_mean_plot3d_sp(io_unit_men)
            call closefile(io_unit_men,trim(solfile)//".mean"//".bak")
            call msg_seq_and_master("input the statistical flow field successfully")
        else
            check = check + 1
        end if
    end if
    
    if (nstepsav == ns_prms) then
        call openfile_check(io_unit_fmen,trim(solfile)//".mean"//".bak","old","unformatted","stream",ierr)
        
#ifdef PARALLEL
        call MPI_BCAST(ierr,1,kind_int_mpi,master,MPI_COMM_WORLD,error)
#endif        
        
        if(ierr == 0 .and. check == 0) then
            call input_mean_plot3d_sp(io_unit_fmen)
            call closefile(io_unit_fmen,trim(solfile)//".mean"//".bak")
            call msg_seq_and_master("input the statistical flow field successfully")
        else
            check = check + 1
        end if
    end if
    
    if (nstepsav > ns_prms) then
        call openfile_check(io_unit_fmen,trim(solfile)//".mean"//".bak","old","unformatted","stream",ierr)
        
#ifdef PARALLEL
        call MPI_BCAST(ierr,1,kind_int_mpi,master,MPI_COMM_WORLD,error)
#endif        
        
        if(ierr == 0 .and. check == 0) then
            call input_mean_plot3d_sp(io_unit_fmen)
            call closefile(io_unit_fmen,trim(solfile)//".mean"//".bak")
        else
            check = check + 1
        end if
    
        call openfile_check(io_unit_rms,trim(solfile)//".rms"//".bak","old","unformatted","stream",ierr)
        
#ifdef PARALLEL
        call MPI_BCAST(ierr,1,kind_int_mpi,master,MPI_COMM_WORLD,error)
#endif        
        
        if(ierr == 0 .and. check == 0) then
            call input_rms_plot3d_sp(io_unit_rms)
            call closefile(io_unit_rms,trim(solfile)//".rms"//".bak")
            call msg_seq_and_master("input the statistical flow field successfully")
        else
            check = check + 1
        end if
    end if
    end if    
    
end subroutine input_sol_backup_sp

subroutine input_tur_sol
    use mod_kndconsts, only : kind_int
    use mod_constants, only : io_unit_tur
    use mod_variables, only : solfile
    implicit none

    call openfile(io_unit_tur,trim(solfile)//".tur","old","unformatted","stream")
    call input_tur_sol_plot3d(io_unit_tur)
    call closefile(io_unit_tur,trim(solfile)//".tur")

    call msg_seq_and_master("input the turbulence flow field successfully")

end subroutine input_tur_sol

subroutine input_tur_sol_sp
    use mod_kndconsts, only : kind_int
    use mod_constants, only : io_unit_tur
    use mod_variables, only : solfile
    implicit none

    call openfile(io_unit_tur,trim(solfile)//".tur","old","unformatted","stream")
    call input_tur_sol_plot3d_sp(io_unit_tur)
    call closefile(io_unit_tur,trim(solfile)//".tur")

    call msg_seq_and_master("input the turbulence flow field successfully")

end subroutine input_tur_sol_sp

subroutine input_dst
    use mod_kndconsts, only : kind_int
    use mod_constants, only : io_unit_dst
    use mod_variables, only : grdfile
    implicit none

    call openfile(io_unit_dst,trim(grdfile)//".dst","old","unformatted","stream")
    call input_dst_plot3d(io_unit_dst)
    call closefile(io_unit_dst,trim(grdfile)//".dst")

end subroutine input_dst

subroutine input_dsp
    use mod_kndconsts, only : kind_int
    use mod_constants, only : io_unit_dsp
    use mod_variables, only : grdfile
    implicit none

    call openfile(io_unit_dsp,trim(grdfile)//".dsp","old","unformatted","stream")
    call input_dsp_plot3d(io_unit_dsp)
    call closefile(io_unit_dsp,trim(grdfile)//".dsp")

end subroutine input_dsp

subroutine input_sol_plot3d(io_unit)
    use mod_kndconsts, only : kind_int
    use mod_constants, only : nplot3d_fun
    use mod_fieldvars, only : mb_pv,npvs
    use mod_interface, only : scatter_input_mb_var
    implicit none
    integer(kind_int), intent(in) :: io_unit

    call scatter_input_mb_var(io_unit,mb_pv,1,npvs,nplot3d_fun)

end subroutine input_sol_plot3d

subroutine input_sol_plot3d_sp(io_unit)
    use mod_kndconsts, only : kind_int
    use mod_constants, only : nplot3d_fun
    use mod_fieldvars, only : mb_pv,npvs
    use mod_interface, only : scatter_input_mb_var_sp
    implicit none
    integer(kind_int), intent(in) :: io_unit

    call scatter_input_mb_var_sp(io_unit,mb_pv,1,npvs,nplot3d_fun)

end subroutine input_sol_plot3d_sp

subroutine input_mean_plot3d(io_unit)
    use mod_kndconsts, only : kind_int
    use mod_constants, only : nplot3d_fun
    use mod_fieldvars, only : mb_fmean,npvs
    use mod_interface, only : scatter_input_mb_var
    implicit none
    integer(kind_int), intent(in) :: io_unit

    call scatter_input_mb_var(io_unit,mb_fmean,1,npvs,nplot3d_fun)

end subroutine input_mean_plot3d

subroutine input_mean_plot3d_sp(io_unit)
    use mod_kndconsts, only : kind_int
    use mod_constants, only : nplot3d_fun
    use mod_fieldvars, only : mb_fmean,npvs
    use mod_interface, only : scatter_input_mb_var_sp
    implicit none
    integer(kind_int), intent(in) :: io_unit

    call scatter_input_mb_var_sp(io_unit,mb_fmean,1,npvs,nplot3d_fun)

end subroutine input_mean_plot3d_sp

subroutine input_rms_plot3d(io_unit)
    use mod_kndconsts, only : kind_int
    use mod_constants, only : nplot3d_fun
    use mod_fieldvars, only : mb_frms,npvs
    use mod_interface, only : scatter_input_mb_var
    implicit none
    integer(kind_int), intent(in) :: io_unit

    call scatter_input_mb_var(io_unit,mb_frms,1,npvs,nplot3d_fun)

end subroutine input_rms_plot3d

subroutine input_rms_plot3d_sp(io_unit)
    use mod_kndconsts, only : kind_int
    use mod_constants, only : nplot3d_fun
    use mod_fieldvars, only : mb_frms,npvs
    use mod_interface, only : scatter_input_mb_var_sp
    implicit none
    integer(kind_int), intent(in) :: io_unit

    call scatter_input_mb_var_sp(io_unit,mb_frms,1,npvs,nplot3d_fun)

end subroutine input_rms_plot3d_sp

subroutine input_tur_sol_plot3d(io_unit)
    use mod_kndconsts, only : kind_int
    use mod_constants, only : nplot3d_fun
    use mod_fieldvars, only : mb_qvst,nqvst
    use mod_interface, only : scatter_input_mb_var
    implicit none
    integer(kind_int), intent(in) :: io_unit

    call scatter_input_mb_var(io_unit,mb_qvst,1,nqvst,nplot3d_fun)

end subroutine input_tur_sol_plot3d

subroutine input_tur_sol_plot3d_sp(io_unit)
    use mod_kndconsts, only : kind_int
    use mod_constants, only : nplot3d_fun
    use mod_fieldvars, only : mb_qvst,nqvst
    use mod_interface, only : scatter_input_mb_var_sp
    implicit none
    integer(kind_int), intent(in) :: io_unit

    call scatter_input_mb_var_sp(io_unit,mb_qvst,1,nqvst,nplot3d_fun)

end subroutine input_tur_sol_plot3d_sp

subroutine input_dst_plot3d(io_unit)
    use mod_kndconsts, only : kind_int
    use mod_constants, only : nplot3d_fun
    use mod_fieldvars, only : mb_dst
    use mod_interface, only : scatter_input_mb_var
    implicit none
    integer(kind_int), intent(in) :: io_unit

    call scatter_input_mb_var(io_unit,mb_dst,1,2,nplot3d_fun)

end subroutine input_dst_plot3d

subroutine input_dsp_plot3d(io_unit)
    use mod_kndconsts, only : kind_int
    use mod_constants, only : nplot3d_fun
    use mod_fieldvars, only : mb_dsp
    use mod_interface, only : scatter_input_mb_var
    implicit none
    integer(kind_int), intent(in) :: io_unit

    call scatter_input_mb_var(io_unit,mb_dsp,1,2,nplot3d_fun)

end subroutine input_dsp_plot3d
