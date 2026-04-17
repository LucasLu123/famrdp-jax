
subroutine analyze_top
    use mod_kndconsts, only : kind_int
    use mod_constants, only : bc_cut1to1,nsw_dir_close
    use mod_constants, only : subc_cut_patched
    use mod_datatypes, only : bc_region_t,top_block_t
    use mod_variables, only : nsw_kdir
    use mod_fieldvars, only : nblocks,mb_top
    use mod_parallels
    implicit none
    integer(kind_int) :: nb,nr,m,n,i,j,k,ierr
    integer(kind_int) :: nregs,bctype,subtype
    integer(kind_int) :: st(3),ed(3),nk,nkst
    integer(kind_int) :: ijk(3),mapmat(3,4)
    type(top_block_t), pointer :: top
    type(bc_region_t), pointer :: reg

    call interface_patched
        
    do nb=1,nblocks
        top => mb_top(nb)

        top%ndst(:) = 1
        top%nded(:) = top%nijk(:)
        if (nsw_kdir == nsw_dir_close) then
            nk = top%nijk(3)
            nkst = (nk+1)/2
            top%ndst(3) = nkst
            top%nded(3) = nkst
        end if

        nregs = top%nregions
        do nr=1,nregs
            reg => top%bcs(nr)

            bctype  = reg%bctype
            subtype = reg%subtype
            if (bctype == bc_cut1to1 .and.  subtype /= subc_cut_patched) then
                call get_bc_map(nb,nr, &
                                reg%s_st ,reg%s_ed , &
                                reg%t_st ,reg%t_ed , &
                                reg%s_nd ,reg%s_lr , &
                                reg%s_ord,reg%s_sgn, &
                                reg%t_nd ,reg%t_lr , &
                                reg%t_ord,reg%t_sgn, &
                                reg%s_dir,reg%t_dir, &
                                reg%mapmat)

                st(:) = reg%s_st
                ed(:) = reg%s_ed
                allocate(reg%mapijk(st(1):ed(1), &
                                    st(2):ed(2), &
                                    st(3):ed(3), &
                                    1:3),stat=ierr)

                mapmat(:,:) = reg%mapmat(:,:)
                do k=st(3),ed(3)
                do j=st(2),ed(2)
                do i=st(1),ed(1)
                    ijk(1) = i
                    ijk(2) = j
                    ijk(3) = k
                    reg%mapijk(i,j,k,1) = sum(ijk(1:3)*mapmat(1,1:3)) + mapmat(1,4)
                    reg%mapijk(i,j,k,2) = sum(ijk(1:3)*mapmat(2,1:3)) + mapmat(2,4)
                    reg%mapijk(i,j,k,3) = sum(ijk(1:3)*mapmat(3,1:3)) + mapmat(3,4)
                end do
                end do
                end do
            else
                reg%s_st(:) = abs(reg%s_st(:))
                reg%s_ed(:) = abs(reg%s_ed(:))
                do m=1,3
                    if (reg%s_st(m) > reg%s_ed(m)) then
                        n = reg%s_ed(m)
                        reg%s_ed(m) = reg%s_st(m)
                        reg%s_st(m) = n
                    else if (reg%s_st(m) == reg%s_ed(m)) then
                        reg%s_nd = m
                        if (reg%s_st(m) == 1) then
                            reg%s_lr = -1
                        else
                            reg%s_lr = 1
                        end if
                    end if
                end do

                reg%s_ord(:) = (/1,2,3/)
                reg%s_sgn(:) = (/1,1,1/)
            end if
        end do
    end do

    do nb=1,nblocks
        top => mb_top(nb)

        nregs = top%nregions
        do nr=1,nregs
            reg => top%bcs(nr)

            reg%nbs = nb
            reg%nrs = nr

            bctype  = reg%bctype
            subtype = reg%subtype
            if (bctype == bc_cut1to1 .and. subtype /= subc_cut_patched) then
                call get_inter_nr(reg%nbs,reg%nrs, &
                                  reg%nbt,reg%nrt)
            end if
        end do
    end do

end subroutine analyze_top

subroutine analyze_top_cc
    use mod_kndconsts, only : kind_int
    use mod_constants, only : bc_cut1to1,nsw_dir_close
    use mod_datatypes, only : bc_region_t,top_block_t
    use mod_variables, only : nsw_kdir
    use mod_fieldvars, only : nblocks,mb_top,mb_topc,mb_topsp
    use mod_parallels
    implicit none
    integer(kind_int) :: nb,nr,m,n,i,j,k,ierr
    integer(kind_int) :: nregs,bctype
    integer(kind_int) :: st(3),ed(3),nk,nkst
    integer(kind_int) :: ijk(3),mapmat(3,4)
    type(top_block_t), pointer :: top
    type(bc_region_t), pointer :: reg
    
    do nb=1,nblocks
        top => mb_top(nb)

        top%ndst(:) = 1
        top%nded(:) = top%nijk(:)
        if (nsw_kdir == nsw_dir_close) then
            nk = top%nijk(3)
            nkst = (nk+1)/2
            top%ndst(3) = nkst
            top%nded(3) = nkst
        end if

        nregs = top%nregions
        do nr=1,nregs
            reg => top%bcs(nr)

            bctype  = reg%bctype
            if (bctype == bc_cut1to1) then
                call get_bc_map(nb,nr, &
                                reg%s_st ,reg%s_ed , &
                                reg%t_st ,reg%t_ed , &
                                reg%s_nd ,reg%s_lr , &
                                reg%s_ord,reg%s_sgn, &
                                reg%t_nd ,reg%t_lr , &
                                reg%t_ord,reg%t_sgn, &
                                reg%s_dir,reg%t_dir, &
                                reg%mapmat)

                st(:) = reg%s_st
                ed(:) = reg%s_ed
                allocate(reg%mapijk(st(1):ed(1), &
                                    st(2):ed(2), &
                                    st(3):ed(3), &
                                    1:3),stat=ierr)

                mapmat(:,:) = reg%mapmat(:,:)
                do k=st(3),ed(3)
                do j=st(2),ed(2)
                do i=st(1),ed(1)
                    ijk(1) = i
                    ijk(2) = j
                    ijk(3) = k
                    reg%mapijk(i,j,k,1) = sum(ijk(1:3)*mapmat(1,1:3)) + mapmat(1,4)
                    reg%mapijk(i,j,k,2) = sum(ijk(1:3)*mapmat(2,1:3)) + mapmat(2,4)
                    reg%mapijk(i,j,k,3) = sum(ijk(1:3)*mapmat(3,1:3)) + mapmat(3,4)
                end do
                end do
                end do
            else
                reg%s_st(:) = abs(reg%s_st(:))
                reg%s_ed(:) = abs(reg%s_ed(:))
                do m=1,3
                    if (reg%s_st(m) > reg%s_ed(m)) then
                        n = reg%s_ed(m)
                        reg%s_ed(m) = reg%s_st(m)
                        reg%s_st(m) = n
                    else if (reg%s_st(m) == reg%s_ed(m)) then
                        reg%s_nd = m
                        if (reg%s_st(m) == 1) then
                            reg%s_lr = -1
                        else
                            reg%s_lr = 1
                        end if
                    end if
                end do

                reg%s_ord(:) = (/1,2,3/)
                reg%s_sgn(:) = (/1,1,1/)
            end if
        end do
    end do

    do nb=1,nblocks
        top => mb_top(nb)

        nregs = top%nregions
        do nr=1,nregs
            reg => top%bcs(nr)

            reg%nbs = nb
            reg%nrs = nr

            bctype  = reg%bctype
            if (bctype == bc_cut1to1) then
                call get_inter_nr(reg%nbs,reg%nrs, &
                                  reg%nbt,reg%nrt)
            end if
        end do
    end do    
        
    do nb=1,nblocks
        top => mb_topc(nb)

        top%ndst(:) = 1
        top%nded(:) = top%nijk(:)
        if (nsw_kdir == nsw_dir_close) then
            nk = top%nijk(3)
            nkst = (nk+1)/2
            top%ndst(3) = nkst
            top%nded(3) = nkst
        end if

        nregs = top%nregions
        do nr=1,nregs
            reg => top%bcs(nr)

            bctype  = reg%bctype
            if (bctype == bc_cut1to1) then
                call get_bc_map(nb,nr, &
                                reg%s_st ,reg%s_ed , &
                                reg%t_st ,reg%t_ed , &
                                reg%s_nd ,reg%s_lr , &
                                reg%s_ord,reg%s_sgn, &
                                reg%t_nd ,reg%t_lr , &
                                reg%t_ord,reg%t_sgn, &
                                reg%s_dir,reg%t_dir, &
                                reg%mapmat)

                st(:) = reg%s_st
                ed(:) = reg%s_ed
                allocate(reg%mapijk(st(1):ed(1), &
                                    st(2):ed(2), &
                                    st(3):ed(3), &
                                    1:3),stat=ierr)

                mapmat(:,:) = reg%mapmat(:,:)
                do k=st(3),ed(3)
                do j=st(2),ed(2)
                do i=st(1),ed(1)
                    ijk(1) = i
                    ijk(2) = j
                    ijk(3) = k
                    reg%mapijk(i,j,k,1) = sum(ijk(1:3)*mapmat(1,1:3)) + mapmat(1,4)
                    reg%mapijk(i,j,k,2) = sum(ijk(1:3)*mapmat(2,1:3)) + mapmat(2,4)
                    reg%mapijk(i,j,k,3) = sum(ijk(1:3)*mapmat(3,1:3)) + mapmat(3,4)
                end do
                end do
                end do
            else
                reg%s_st(:) = abs(reg%s_st(:))
                reg%s_ed(:) = abs(reg%s_ed(:))
                do m=1,3
                    if (reg%s_st(m) > reg%s_ed(m)) then
                        n = reg%s_ed(m)
                        reg%s_ed(m) = reg%s_st(m)
                        reg%s_st(m) = n
                    else if (reg%s_st(m) == reg%s_ed(m)) then
                        reg%s_nd = m
                        if (reg%s_st(m) == 1) then
                            reg%s_lr = -1
                        else
                            reg%s_lr = 1
                        end if
                    end if
                end do

                reg%s_ord(:) = (/1,2,3/)
                reg%s_sgn(:) = (/1,1,1/)
            end if
        end do
    end do

    do nb=1,nblocks
        top => mb_topc(nb)

        nregs = top%nregions
        do nr=1,nregs
            reg => top%bcs(nr)

            reg%nbs = nb
            reg%nrs = nr

            bctype  = reg%bctype
            if (bctype == bc_cut1to1) then
                call get_inter_nr_cc(reg%nbs,reg%nrs, &
                                  reg%nbt,reg%nrt)
            end if
        end do
    end do
    
    do nb=1,nblocks
        top => mb_topsp(nb)

        top%ndst(:) = 1
        top%nded(:) = top%nijk(:)
        if (nsw_kdir == nsw_dir_close) then
            nk = top%nijk(3)
            nkst = (nk+1)/2
            top%ndst(3) = nkst
            top%nded(3) = nkst
        end if

        nregs = top%nregions
        do nr=1,nregs
            reg => top%bcs(nr)

            bctype  = reg%bctype
            if (bctype == bc_cut1to1) then
                call get_bc_map(nb,nr, &
                                reg%s_st ,reg%s_ed , &
                                reg%t_st ,reg%t_ed , &
                                reg%s_nd ,reg%s_lr , &
                                reg%s_ord,reg%s_sgn, &
                                reg%t_nd ,reg%t_lr , &
                                reg%t_ord,reg%t_sgn, &
                                reg%s_dir,reg%t_dir, &
                                reg%mapmat)

                st(:) = reg%s_st
                ed(:) = reg%s_ed
                allocate(reg%mapijk(st(1):ed(1), &
                                    st(2):ed(2), &
                                    st(3):ed(3), &
                                    1:3),stat=ierr)

                mapmat(:,:) = reg%mapmat(:,:)
                do k=st(3),ed(3)
                do j=st(2),ed(2)
                do i=st(1),ed(1)
                    ijk(1) = i
                    ijk(2) = j
                    ijk(3) = k
                    reg%mapijk(i,j,k,1) = sum(ijk(1:3)*mapmat(1,1:3)) + mapmat(1,4)
                    reg%mapijk(i,j,k,2) = sum(ijk(1:3)*mapmat(2,1:3)) + mapmat(2,4)
                    reg%mapijk(i,j,k,3) = sum(ijk(1:3)*mapmat(3,1:3)) + mapmat(3,4)
                end do
                end do
                end do
            else
                reg%s_st(:) = abs(reg%s_st(:))
                reg%s_ed(:) = abs(reg%s_ed(:))
                do m=1,3
                    if (reg%s_st(m) > reg%s_ed(m)) then
                        n = reg%s_ed(m)
                        reg%s_ed(m) = reg%s_st(m)
                        reg%s_st(m) = n
                    else if (reg%s_st(m) == reg%s_ed(m)) then
                        reg%s_nd = m
                        if (reg%s_st(m) == 1) then
                            reg%s_lr = -1
                        else
                            reg%s_lr = 1
                        end if
                    end if
                end do

                reg%s_ord(:) = (/1,2,3/)
                reg%s_sgn(:) = (/1,1,1/)
            end if
        end do
    end do

    do nb=1,nblocks
        top => mb_topsp(nb)

        nregs = top%nregions
        do nr=1,nregs
            reg => top%bcs(nr)

            reg%nbs = nb
            reg%nrs = nr

            bctype  = reg%bctype
            if (bctype == bc_cut1to1) then
                call get_inter_nr_sp(reg%nbs,reg%nrs, &
                                  reg%nbt,reg%nrt)
            end if
        end do
    end do    

end subroutine analyze_top_cc

subroutine analyze_top_patched
    use mod_kndconsts, only : kind_int,kind_real,kind_long
    use mod_constants, only : bc_cut1to1,nsw_dir_close
    use mod_constants, only : subc_cut_patched
    use mod_datatypes, only : bc_region_t,top_block_t
    use mod_variables, only : nsw_kdir
    use mod_fieldvars, only : nblocks,mb_top,scal
    use mod_fieldvars, only : nblkcoms,blkcoms
    use mod_fieldvars, only : ninters,ninterf,bc_seq,surf_pxyz,surf_check
    use mod_fieldvars, only : s_dir_all,L_check_all,P_check_all,scal,scalt
    use mod_parallels
    implicit none
    integer(kind_int) :: nb,nc,nr,m,n,i,j,k,ierr
    integer(kind_int) :: nregs,bctype,subtype
    integer(kind_int) :: st(3),ed(3),nk,nkst
    integer(kind_int) :: t_st(3),t_ed(3)
    integer(kind_int) :: ijk(3),mapmat(3,4),tempijk(3)
    type(top_block_t), pointer :: top
    type(bc_region_t), pointer :: reg
    integer(kind_int) :: t_lr,td
    integer(kind_int) :: s_dir(3),t_dir(3)
    integer(kind_int) :: t_sd,t_td,bc_id,s_nd
    real(kind_real)    :: L_check(2),P_check(4)
#ifdef PARALLEL
    integer(kind_int) :: pid
    integer(kind_int),parameter :: packsize=10000
    integer(kind_long)          :: packbuf(packsize)
    integer(kind_int)           :: position
#endif

    allocate(s_dir_all(1:ninters-ninterf,1:3))
    allocate(L_check_all(1:ninters-ninterf,1:2),P_check_all(1:ninters-ninterf,1:4))
    allocate(scal(1:ninters-ninterf,1:3),scalt(1:ninters-ninterf,1:3))
        
    do nb=1,nblocks
        top => mb_top(nb)
        nregs = top%nregions
        do nr=1,nregs
            reg => top%bcs(nr)
            bctype  = reg%bctype
            subtype = reg%subtype
#ifdef PARALLEL
            pid  = mb_top(nb)%pid
#endif
            if (bctype == bc_cut1to1 .and.  subtype == subc_cut_patched) then
                do m=1,ninters-ninterf                    
                    if(nb == bc_seq(m,1) .and. nr == bc_seq(m,2)) then
                        bc_id = m                        
                    end if
                end do
#ifdef PARALLEL
                if (myid == pid) then
#endif
                    call get_bc_dir_patched_s(nb,nr,bc_id,s_dir,L_check,P_check)

#ifdef PARALLEL
                    position = 0
                    call MPI_PACK(s_dir,3,kind_int_mpi, &
                          packbuf,packsize,position,MPI_COMM_WORLD,ierr)
                    call MPI_PACK(L_check,2,kind_real_mpi, &
                          packbuf,packsize,position,MPI_COMM_WORLD,ierr)
                    call MPI_PACK(P_check,4,kind_real_mpi, &
                          packbuf,packsize,position,MPI_COMM_WORLD,ierr)
                end if

                call MPI_BCAST(position,1,kind_int_mpi,pid,MPI_COMM_WORLD,ierr)
                call MPI_BCAST(packbuf,position,MPI_PACKED,pid,MPI_COMM_WORLD,ierr)

                if (myid /= pid) then
                    position = 0

                    call MPI_UNPACK(packbuf,packsize,position, &
                                    s_dir,3,kind_int_mpi,MPI_COMM_WORLD,ierr)
                    call MPI_UNPACK(packbuf,packsize,position, &
                                    L_check,2,kind_real_mpi,MPI_COMM_WORLD,ierr)
                    call MPI_UNPACK(packbuf,packsize,position, &
                                    P_check,4,kind_real_mpi,MPI_COMM_WORLD,ierr)
                end if
#endif
                s_dir_all(bc_id,:)   = s_dir(:)
                L_check_all(bc_id,:) = L_check(:)
                P_check_all(bc_id,:) = P_check(:)
            end if
        end do
    end do

    do nb=1,nblocks
        top => mb_top(nb)

        top%ndst(:) = 1
        top%nded(:) = top%nijk(:)
        if (nsw_kdir == nsw_dir_close) then
            nk = top%nijk(3)
            nkst = (nk+1)/2
            top%ndst(3) = nkst
            top%nded(3) = nkst
        end if

        nregs = top%nregions
        do nr=1,nregs
            reg => top%bcs(nr)

            st(:) = reg%s_st
            ed(:) = reg%s_ed
            s_nd  = reg%s_nd

            bctype  = reg%bctype
            subtype = reg%subtype
            if (bctype == bc_cut1to1 .and.  subtype == subc_cut_patched) then
                do m=1,ninters-ninterf                    
                    if(nb == bc_seq(m,1) .and. nr == bc_seq(m,2)) then
                        bc_id = m                        
                    end if
                end do

                s_dir(:) = s_dir_all(bc_id,:)
                L_check(:) = L_check_all(bc_id,:)
                P_check(:) = P_check_all(bc_id,:)
                call get_bc_dir_patched_t(bc_id,reg%t_st,reg%t_ed,reg%t_lr,t_dir,L_check,P_check)
                call get_bc_map_patched(bc_id, &
                                        s_dir,t_dir, &
                                        reg%s_ord,reg%s_sgn, &
                                        reg%t_ord,reg%t_sgn, &
                                        reg%mapmat)

                allocate(reg%mapijk(st(1):ed(1), &
                                    st(2):ed(2), &
                                    st(3):ed(3), &
                                    1:3),stat=ierr)

                t_st(:) = reg%t_st
                t_ed(:) = reg%t_ed                
                mapmat(:,:) = reg%mapmat(:,:)
                do k=st(3),ed(3)
                do j=st(2),ed(2)
                do i=st(1),ed(1)
                    ijk(1) = i
                    ijk(2) = j
                    ijk(3) = k
                    tempijk(1) = int(scalt(bc_id,1)*(sum(ijk(1:3)*mapmat(1,1:3)) + mapmat(1,4) - t_st(1)) +  t_st(1))
                    tempijk(2) = int(scalt(bc_id,2)*(sum(ijk(1:3)*mapmat(2,1:3)) + mapmat(2,4) - t_st(2)) +  t_st(2))
                    tempijk(3) = int(scalt(bc_id,3)*(sum(ijk(1:3)*mapmat(3,1:3)) + mapmat(3,4) - t_st(3)) +  t_st(3))

                    do m=1,3
                        if(tempijk(m)<min(t_st(m),t_ed(m))) then
                            tempijk(m) = min(t_st(m),t_ed(m))
                        end if
                    end do

                    do m=1,3
                        if(tempijk(m)>max(t_st(m),t_ed(m))) then
                            tempijk(m) = max(t_st(m),t_ed(m))
                        end if
                    end do

                    reg%mapijk(i,j,k,1) = tempijk(1)
                    reg%mapijk(i,j,k,2) = tempijk(2)
                    reg%mapijk(i,j,k,3) = tempijk(3)
                end do
                end do
                end do

                do m=2,3
                    td = t_dir(m)
                    if (t_st(td) > t_ed(td)) then
                        n = t_ed(td)
                        t_ed(td) = t_st(td)
                        t_st(td) = n
                    end if
                end do

                reg%t_st = t_st(:)
                reg%t_ed = t_ed(:)

            end if
        end do
    end do

    deallocate(surf_pxyz,surf_check,s_dir_all,L_check_all,P_check_all)

    do nb=1,nblocks
        top => mb_top(nb)

        nregs = top%nregions
        do nr=1,nregs
            reg => top%bcs(nr)

            reg%nbs = nb
            reg%nrs = nr

            bctype  = reg%bctype
            subtype = reg%subtype
            if (bctype == bc_cut1to1 .and. subtype == subc_cut_patched) then
                call get_inter_nr(reg%nbs,reg%nrs, &
                                  reg%nbt,reg%nrt)
            end if
        end do
    end do

end subroutine analyze_top_patched

subroutine interface_patched
    use mod_kndconsts, only : kind_int
    use mod_constants, only : bc_cut1to1
    use mod_constants, only : subc_cut_patched
    use mod_datatypes, only : bc_region_t,top_block_t
    use mod_fieldvars, only : nblocks,mb_top
    implicit none
    integer(kind_int) :: nb,nr,m,n,i,j,k,ierr
    integer(kind_int) :: nregs,bctype,subtype
    integer(kind_int) :: nbt
    type(top_block_t), pointer :: top
    type(bc_region_t), pointer :: reg

    do nb=1,nblocks
        top => mb_top(nb)
        nregs = top%nregions
        do nr=1,nregs
            reg => top%bcs(nr)
            subtype = reg%subtype
            if (subtype == subc_cut_patched) then

                do m=1,3
                    reg%t_st(m) = -10000000
                    reg%t_ed(m) = -10000000
                end do
                reg%nbt = -10000000
                reg%nrs = nr

                reg%s_st(:) = abs(reg%s_st(:))
                reg%s_ed(:) = abs(reg%s_ed(:))
                do m=1,3
                    if (reg%s_st(m) > reg%s_ed(m)) then
                        n = reg%s_ed(m)
                        reg%s_ed(m) = reg%s_st(m)
                        reg%s_st(m) = n
                    else if (reg%s_st(m) == reg%s_ed(m)) then
                        reg%s_nd = m
                        if (reg%s_st(m) == 1) then
                            reg%s_lr = -1
                        else
                            reg%s_lr = 1
                        end if
                    end if
                end do
            end if
        end do
    end do

end subroutine interface_patched

subroutine set_bc_patched
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : bc_cut1to1
    use mod_constants, only : subc_cut_patched
    use mod_datatypes, only : bc_region_t,top_block_t
    use mod_fieldvars, only : nblocks,mb_top
    use mod_fieldvars, only : ninters,ninterf
    use mod_fieldvars, only : bc_seq,surf_pxyz,surf_check
    use mod_parallels
    implicit none
    integer(kind_int) :: nb,nr,id,ierr
    integer(kind_int) :: s_nd,st(3),ed(3)
    integer(kind_int) :: nregs,bctype,subtype
    integer(kind_int) :: t_nb,t_nr,t_s_nd,t_nregs,t_subtype
    integer(kind_int) :: t_st(3),t_ed(3)
    integer(kind_int) :: nbt,m,n,l,bc_id
    real(kind_real)    :: pxyz(1:3,1:4),check,t_check
    real(kind_real)    :: psurf_check(1:ninters-ninterf)
    real(kind_real)    :: P_pxyz(1:(ninters-ninterf)*12),S_pxyz(1:(ninters-ninterf)*12)
    type(top_block_t), pointer :: top
    type(bc_region_t), pointer :: reg
    type(top_block_t), pointer :: t_top
    type(bc_region_t), pointer :: t_reg

    allocate(bc_seq(1:ninters-ninterf,1:3),surf_pxyz(1:ninters-ninterf,1:3,1:4),surf_check(1:ninters-ninterf))

    bc_id = 0
    do nb=1,nblocks
        top => mb_top(nb)
        nregs = top%nregions
        do nr=1,nregs
            reg => top%bcs(nr)
            bctype  = reg%bctype
            subtype = reg%subtype
            if (bctype == bc_cut1to1 .and.  subtype == subc_cut_patched) then
                bc_id = bc_id + 1
                bc_seq(bc_id,1) = nb
                bc_seq(bc_id,2) = nr
                bc_seq(bc_id,3) = -10000
            end if
        end do
    end do

    surf_check(:) = 0.0
    surf_pxyz(:,:,:) = 0.0
    do nb=1,nblocks
        top => mb_top(nb)
        nregs = top%nregions
        do nr=1,nregs
            reg => top%bcs(nr)
            bctype  = reg%bctype
            subtype = reg%subtype
            if (bctype == bc_cut1to1) then
                nbt     = reg%nbt
                if (subtype == subc_cut_patched .and. nbt < 0) then
                    do m=1,ninters-ninterf
                        if(nb == bc_seq(m,1) .and. nr == bc_seq(m,2)) then
                            bc_id = m                        
                        end if
                    end do
                    call find_patched_check(nb,nr,check,pxyz)
                    psurf_check(bc_id)   = check
                    do m=1,3
                        do n=1,4
                            surf_pxyz(bc_id,m,n) = pxyz(m,n)
                        end do
                    end do                    
                end if
            end if
        end do
    end do

    surf_check(:)    = psurf_check(:)

    do m=1,ninters-ninterf
        do n=1,3
            do l=1,4
                id = l + 4*(n-1) + 12*(m-1)
                P_pxyz(id) = surf_pxyz(m,n,l)
                S_pxyz(id) = P_pxyz(id)
            end do
        end do
    end do

#ifdef PARALLEL

    call MPI_REDUCE(psurf_check,surf_check,ninters-ninterf,kind_real_mpi,MPI_SUM,master,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(surf_check,ninters-ninterf,kind_real_mpi,master,MPI_COMM_WORLD,ierr)
    call MPI_REDUCE(P_pxyz,S_pxyz,(ninters-ninterf)*12,kind_real_mpi,MPI_SUM,master,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(S_pxyz,(ninters-ninterf)*12,kind_real_mpi,master,MPI_COMM_WORLD,ierr)

#endif

    do m=1,ninters-ninterf
        do n=1,3
            do l=1,4
                id = l + 4*(n-1) + 12*(m-1)
                surf_pxyz(m,n,l) = S_pxyz(id)
            end do
        end do
    end do

    do m=1,ninters-ninterf
        nb = bc_seq(m,1)
        nr = bc_seq(m,2)
        top => mb_top(nb)
        reg => top%bcs(nr)
        st(:) = reg%s_st(:)
        ed(:) = reg%s_ed(:)
        s_nd  = reg%s_nd
        check = surf_check(m)
        do n = 1,ninters-ninterf
            if(m /=n) then
                t_nb = bc_seq(n,1)
                t_nr = bc_seq(n,2)
                t_top => mb_top(t_nb)
                t_reg => t_top%bcs(t_nr)
                t_st(:)   = t_reg%s_st(:)
                t_ed(:)   = t_reg%s_ed(:)
                t_s_nd    = t_reg%s_nd
                t_check   = surf_check(n)
                if(abs(1.0-t_check**2.0/(check**2.0+1.0e-30))<1.0e-3) then
                    bc_seq(m,3) = n

                    mb_top(nb)%bcs(nr)%nbt  = t_nb
                    mb_top(nb)%bcs(nr)%nrt  = t_nr
                    mb_top(nb)%bcs(nr)%t_nd = t_s_nd
                    mb_top(nb)%bcs(nr)%t_st = t_st(:)
                    mb_top(nb)%bcs(nr)%t_ed = t_ed(:)

                    mb_top(t_nb)%bcs(t_nr)%nbt   = nb
                    mb_top(t_nb)%bcs(t_nr)%nrt   = nr
                    mb_top(t_nb)%bcs(t_nr)%t_nd  = s_nd
                    mb_top(t_nb)%bcs(t_nr)%t_st  = st(:)
                    mb_top(t_nb)%bcs(t_nr)%t_ed  = ed(:)
                end if
            end if
        end do
    end do

    ierr = 0
    do nb=1,nblocks
        top => mb_top(nb)
        nregs = top%nregions
        do nr=1,nregs
            reg => top%bcs(nr)
            bctype  = reg%bctype
            subtype = reg%subtype
            if (bctype == bc_cut1to1) then
                nbt     = reg%nbt
                if (subtype == subc_cut_patched .and. nbt < 0) then
                    ierr = nb
                    goto 10
                end if
            end if
        end do
    end do

10  call error_check(ierr,"Some boundaries are not patched")

end subroutine set_bc_patched

subroutine find_patched_check(nb,nr,check,pxyz)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : subc_cut_patched
    use mod_datatypes, only : bc_region_t,top_block_t
    use mod_fieldvars, only : nblocks,mb_top,mb_xyz
    use mod_parallels
    implicit none
    integer(kind_int), intent(in) :: nb,nr
    real(kind_real),   intent(out) :: pxyz(1:3,1:4),check
    integer(kind_int)              :: s_nd,st(3),ed(3)
    integer(kind_int)              :: m,n,ierr
    type(top_block_t), pointer    :: top
    type(bc_region_t), pointer    :: reg
#ifdef PARALLEL
    integer(kind_int) :: pid
#endif

    top => mb_top(nb)
    reg => top%bcs(nr)
    st(:) = reg%s_st(:)
    ed(:) = reg%s_ed(:)
    s_nd  = reg%s_nd

    check = 0.0
    pxyz(:,:) = 0.0

#ifdef PARALLEL
    pid  = mb_top(nb)%pid
    if (myid == pid) then
#endif

    if (s_nd == 1) then
        do m=1,3
            pxyz(m,1) = mb_xyz(nb)%fld(m)%r3d(st(1),st(2),st(3))
            pxyz(m,2) = mb_xyz(nb)%fld(m)%r3d(st(1),ed(2),st(3))
            pxyz(m,3) = mb_xyz(nb)%fld(m)%r3d(st(1),st(2),ed(3))
            pxyz(m,4) = mb_xyz(nb)%fld(m)%r3d(st(1),ed(2),ed(3))
        end do
    else if (s_nd == 2) then
        do m=1,3
            pxyz(m,1) = mb_xyz(nb)%fld(m)%r3d(st(1),st(2),st(3))
            pxyz(m,2) = mb_xyz(nb)%fld(m)%r3d(ed(1),st(2),st(3))
            pxyz(m,3) = mb_xyz(nb)%fld(m)%r3d(st(1),st(2),ed(3))
            pxyz(m,4) = mb_xyz(nb)%fld(m)%r3d(ed(1),st(2),ed(3))
        end do
    else
        do m=1,3
            pxyz(m,1) = mb_xyz(nb)%fld(m)%r3d(st(1),st(2),st(3))
            pxyz(m,2) = mb_xyz(nb)%fld(m)%r3d(ed(1),st(2),st(3))
            pxyz(m,3) = mb_xyz(nb)%fld(m)%r3d(st(1),ed(2),st(3))
            pxyz(m,4) = mb_xyz(nb)%fld(m)%r3d(ed(1),ed(2),st(3))
        end do
    end if

    do m=1,3
        do n=1,4
            check = check +pxyz(m,n)
        end do
    end do

#ifdef PARALLEL
    end if
#endif

end subroutine find_patched_check

subroutine get_bc_map(nb,nr,s_st,s_ed,t_st,t_ed, &
                      s_nd,s_lr,s_ord,s_sgn, &
                      t_nd,t_lr,t_ord,t_sgn, &
                      s_dir,t_dir, &
                      mapmat)
    use mod_kndconsts, only : kind_int
    use mod_constants, only : mide
    implicit none
    integer(kind_int),intent(in) :: nb,nr
    integer(kind_int),intent(inout) :: s_st(3),s_ed(3)
    integer(kind_int),intent(inout) :: t_st(3),t_ed(3)
    integer(kind_int),intent(inout) :: s_dir(3),t_dir(3)
    integer(kind_int),intent(out)   :: s_nd,s_lr
    integer(kind_int),intent(out)   :: t_nd,t_lr
    integer(kind_int),intent(out)   :: s_ord(3),s_sgn(3)
    integer(kind_int),intent(out)   :: t_ord(3),t_sgn(3)
    integer(kind_int),intent(out)   :: mapmat(3,4)
    integer(kind_int) :: m,n,sd,td

    call get_bc_dir(nb,nr,s_st,s_ed,s_lr,s_dir)
    call get_bc_dir(nb,nr,t_st,t_ed,t_lr,t_dir)

    s_nd = s_dir(1)
    t_nd = t_dir(1)

    do m=1,3
        s_ord(t_dir(m)) = s_dir(m)
        t_ord(s_dir(m)) = t_dir(m)
    end do

    do m=2,3
        sd = s_dir(m)
        td = t_dir(m)
        if (s_st(sd) > s_ed(sd)) then
            n = s_ed(sd)
            s_ed(sd) = s_st(sd)
            s_st(sd) = n

            n = t_ed(td)
            t_ed(td) = t_st(td)
            t_st(td) = n
        end if
    end do

    s_sgn(:) = (/1,1,1/)
    do m=1,3
        td = t_ord(m)
        t_sgn(m) = sign(1,t_ed(td)-t_st(td))
    end do

    do m=1,3
       mapmat(m,4) = t_st(m)
       sd = s_ord(m)
       do n=1,3
          mapmat(m,n) = t_sgn(n)*mide(n,sd)
          mapmat(m,4) = mapmat(m,4) - mapmat(m,n)*s_st(n)
       end do
    end do

    t_sgn(s_nd) = -s_lr*t_lr  !����������

    do m=2,3
        td = t_dir(m)
        if (t_st(td) > t_ed(td)) then
            n = t_ed(td)
            t_ed(td) = t_st(td)
            t_st(td) = n
        end if
    end do

end subroutine get_bc_map

subroutine get_bc_map_patched(bc_id,      &
                               s_dir,t_dir, &
                               s_ord,s_sgn, &
                               t_ord,t_sgn, &
                               mapmat)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : mide
    use mod_fieldvars, only : ninters,ninterf
    use mod_fieldvars, only : nblocks,mb_top,bc_seq,scal,scalt
    implicit none
    integer(kind_int),intent(in) :: bc_id
    integer(kind_int),intent(in) :: s_dir(3),t_dir(3)
    integer(kind_int),intent(out)   :: s_ord(3),s_sgn(3)
    integer(kind_int),intent(out)   :: t_ord(3),t_sgn(3)
    integer(kind_int),intent(out)   :: mapmat(3,4)
    integer(kind_int) :: nb,nr,s_st(3),s_ed(3),t_st(3),t_ed(3),s_nd,s_lr,t_lr
    integer(kind_int) :: m,n,sd,td,lengths,lengtht
    real(kind_real)    :: L_check(2),P_check(1:4)
    
    nb = bc_seq(bc_id,1)
    nr = bc_seq(bc_id,2)

    mb_top(nb)%bcs(nr)%s_dir(:) = s_dir(:)
    mb_top(nb)%bcs(nr)%t_dir(:) = t_dir(:)

    s_st(:) = mb_top(nb)%bcs(nr)%s_st(:)
    s_ed(:) = mb_top(nb)%bcs(nr)%s_ed(:)
    t_st(:) = mb_top(nb)%bcs(nr)%t_st(:)
    t_ed(:) = mb_top(nb)%bcs(nr)%t_ed(:)
    s_nd    = mb_top(nb)%bcs(nr)%s_nd
    s_lr    = mb_top(nb)%bcs(nr)%s_lr
    t_lr    = mb_top(nb)%bcs(nr)%t_lr
    
    do m=1,3
        s_ord(t_dir(m)) = s_dir(m)
        t_ord(s_dir(m)) = t_dir(m)
    end do

    s_sgn(:) = (/1,1,1/)
    do m=1,3
        td = t_ord(m)
        t_sgn(m) = sign(1,t_ed(td)-t_st(td))
    end do

    do m=1,3
       mapmat(m,4) = t_st(m)
       sd = s_ord(m)
       do n=1,3
          mapmat(m,n) = t_sgn(n)*mide(n,sd)
          mapmat(m,4) = mapmat(m,4) - mapmat(m,n)*s_st(n)
       end do        
    end do

    do m=1,3
        lengths = abs(s_ed(s_dir(m))-s_st(s_dir(m)))
        lengtht = abs(t_ed(t_dir(m))-t_st(t_dir(m)))
        if(lengths == lengtht) then
            scal(bc_id,s_dir(m)) = 1.0            
        else
            scal(bc_id,s_dir(m)) = real((lengtht+1.0e-30)/(lengths+1.0e-10))
        end if
        scalt(bc_id,t_dir(m)) = scal(bc_id,s_dir(m))
    end do

    t_sgn(s_nd) = -s_lr*t_lr  !����������

end subroutine get_bc_map_patched

subroutine get_bc_dir(nb,nr,st,ed,lr,dir)
    use mod_kndconsts, only : kind_int
    implicit none
    integer(kind_int), intent(in) :: nb,nr
    integer(kind_int), intent(inout) :: st(3),ed(3)
    integer(kind_int), intent(out) :: lr,dir(3)
    integer(kind_int) :: m,ierr

    ierr = 0

    dir(:) = 0
    do m=1,3
        if (st(m) == ed(m)) then
            dir(1) = m
            if (st(m) == 1) then
                lr = -1
            else
                lr = 1
            end if
        else if (st(m) < 0) then
            dir(2) = m
        end if
    end do
    dir(3) = 6 - dir(1) - dir(2)

    do m=1,3
        if (dir(m) == 0) then
            ierr = 1
            write(*,'(2(a,i5))') " # BLOCK:",nb,", BC:",nr
            exit
        end if
    end do

    call error_check(ierr,"BC file is incorrect in subroutine get_bc_dir!")

    st(:) = abs(st(:))
    ed(:) = abs(ed(:))

end subroutine get_bc_dir

subroutine get_bc_dir_patched_s(nb,nr,bc_id,dir,dir_check,p_check)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_fieldvars, only : nblocks,mb_top,mb_xyz
    use mod_fieldvars, only : s_dir_all,L_check_all,P_check_all
    use mod_parallels
    implicit none
    integer(kind_int), intent(in)  :: nb,nr,bc_id
    integer(kind_int), intent(out) :: dir(1:3)
    integer(kind_int)               :: m,s_nd,st(3),ed(3)
    real(kind_real), intent(out)    :: dir_check(1:2),p_check(1:4)
    real(kind_real)               :: pxyz(1:3,1:4)

    dir(:)       = 0
    dir_check(:) = 0.0

    s_nd   = mb_top(nb)%bcs(nr)%s_nd
    st(:)  = mb_top(nb)%bcs(nr)%s_st(:)
    ed(:)  = mb_top(nb)%bcs(nr)%s_ed(:)
    dir(1) = s_nd

    if (s_nd == 1) then
        do m=1,3
            pxyz(m,1) = mb_xyz(nb)%fld(m)%r3d(st(1),st(2),st(3))
            pxyz(m,2) = mb_xyz(nb)%fld(m)%r3d(st(1),ed(2),st(3))
            pxyz(m,3) = mb_xyz(nb)%fld(m)%r3d(st(1),st(2),ed(3))
            pxyz(m,4) = mb_xyz(nb)%fld(m)%r3d(st(1),ed(2),ed(3))
        end do
        dir(2) = 2
        dir(3) = 3
    else if (s_nd == 2) then
        do m=1,3
            pxyz(m,1) = mb_xyz(nb)%fld(m)%r3d(st(1),st(2),st(3))
            pxyz(m,2) = mb_xyz(nb)%fld(m)%r3d(ed(1),st(2),st(3))
            pxyz(m,3) = mb_xyz(nb)%fld(m)%r3d(st(1),st(2),ed(3))
            pxyz(m,4) = mb_xyz(nb)%fld(m)%r3d(ed(1),st(2),ed(3))
        end do
        dir(2) = 1
        dir(3) = 3
    else
        do m=1,3
            pxyz(m,1) = mb_xyz(nb)%fld(m)%r3d(st(1),st(2),st(3))
            pxyz(m,2) = mb_xyz(nb)%fld(m)%r3d(ed(1),st(2),st(3))
            pxyz(m,3) = mb_xyz(nb)%fld(m)%r3d(st(1),ed(2),st(3))
            pxyz(m,4) = mb_xyz(nb)%fld(m)%r3d(ed(1),ed(2),st(3))
        end do
        dir(2) = 1
        dir(3) = 2
    end if

    do m=1,3
        dir_check(1) = dir_check(1) + (pxyz(m,1)- pxyz(m,2))**2.0 + (pxyz(m,3)- pxyz(m,4))**2.0
        dir_check(2) = dir_check(2) + (pxyz(m,1)- pxyz(m,3))**2.0 + (pxyz(m,2)- pxyz(m,4))**2.0
    end do

    do m=1,4
        p_check(1) = pxyz(1,1) + pxyz(2,1) + pxyz(3,1)
        p_check(2) = pxyz(1,2) + pxyz(2,2) + pxyz(3,2)
        p_check(3) = pxyz(1,3) + pxyz(2,3) + pxyz(3,3)
        p_check(4) = pxyz(1,4) + pxyz(2,4) + pxyz(3,4)
    end do

end subroutine get_bc_dir_patched_s

subroutine get_bc_dir_patched_t(bc_id,st,ed,lr,dir,dir_check,p_check)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_fieldvars, only : nblocks,mb_top
    use mod_fieldvars, only : bc_seq,surf_pxyz
    use mod_parallels
    implicit none
    integer(kind_int), intent(in) :: bc_id
    integer(kind_int), intent(inout) :: st(3),ed(3)
    integer(kind_int), intent(out)   :: lr,dir(3)
    integer(kind_int)                 :: nb,nr,m,n,t_bc_id,td,nbt,nrt,t_nd
    integer(kind_int)                 :: direction(2)
    real(kind_real),    intent(in)    :: dir_check(2),p_check(1:4)
    real(kind_real)                    :: pxyz(1:3,1:4),check(1:2),tp_check(1:4)

    direction(:) = 0

    dir(:)      = 0
    check(:)    = 0.0
    tp_check(:) = 0.0

    nb = bc_seq(bc_id,1)
    nr = bc_seq(bc_id,2)
    t_bc_id = bc_seq(bc_id,3)
        
    nbt    = mb_top(nb)%bcs(nr)%nbt
    nrt    = mb_top(nb)%bcs(nr)%nrt
    t_nd   = mb_top(nb)%bcs(nr)%t_nd
    dir(1) = t_nd
    st(:)  = mb_top(nbt)%bcs(nrt)%s_st(:)
    ed(:)  = mb_top(nbt)%bcs(nrt)%s_ed(:)

    if (t_nd == 1) then
        if (st(1) == 1) then
            lr = -1
        else
            lr = 1
        end if
               
        do m=1,3
            pxyz(m,1) = surf_pxyz(t_bc_id,m,1)
            pxyz(m,2) = surf_pxyz(t_bc_id,m,2)
            pxyz(m,3) = surf_pxyz(t_bc_id,m,3)
            pxyz(m,4) = surf_pxyz(t_bc_id,m,4)
        end do

        do m=1,3
            check(1) = check(1) + (pxyz(m,1)- pxyz(m,2))**2.0 + (pxyz(m,3)- pxyz(m,4))**2.0
            check(2) = check(2) + (pxyz(m,1)- pxyz(m,3))**2.0 + (pxyz(m,2)- pxyz(m,4))**2.0
        end do

        if(abs(1.0-check(1)**2.0/(dir_check(1)**2.0+1.0e-30))<1.0e-3 .and. &
           abs(1.0-check(2)**2.0/(dir_check(2)**2.0+1.0e-30))<1.0e-3) then
            dir(2)    = 2
            dir(3)    = 3
            do m=1,4
                tp_check(1) = pxyz(1,1) + pxyz(2,1) + pxyz(3,1)
                tp_check(2) = pxyz(1,2) + pxyz(2,2) + pxyz(3,2)
                tp_check(3) = pxyz(1,3) + pxyz(2,3) + pxyz(3,3)
                tp_check(4) = pxyz(1,4) + pxyz(2,4) + pxyz(3,4)
            end do
        else
            dir(2)    = 3
            dir(3)    = 2
            do m=1,4
                tp_check(1) = pxyz(1,1) + pxyz(2,1) + pxyz(3,1)
                tp_check(2) = pxyz(1,3) + pxyz(2,3) + pxyz(3,3)
                tp_check(3) = pxyz(1,2) + pxyz(2,2) + pxyz(3,2)
                tp_check(4) = pxyz(1,4) + pxyz(2,4) + pxyz(3,4)
            end do
        end if

        if(abs(1.0-tp_check(1)**2.0/(p_check(1)**2.0+1.0e-30))<1.0e-3 .and. &
           abs(1.0-tp_check(2)**2.0/(p_check(2)**2.0+1.0e-30))<1.0e-3) then
            direction(1) = 1
        else
            direction(1) = 2
        end if
        if(abs(1.0-tp_check(1)**2.0/(p_check(1)**2.0+1.0e-30))<1.0e-3 .and. &
           abs(1.0-tp_check(3)**2.0/(p_check(3)**2.0+1.0e-30))<1.0e-3) then
            direction(2) = 1
        else
            direction(2) = 2
        end if

        do m=1,2
            if(direction(m) == 2) then
                td = dir(m+1)
                n = ed(td)
                ed(td) = st(td)
                st(td) = n
            end if
        end do

    else if (t_nd == 2) then
        if (st(2) == 1) then
            lr = -1
        else
            lr = 1
        end if

        do m=1,3
            pxyz(m,1) = surf_pxyz(t_bc_id,m,1)
            pxyz(m,2) = surf_pxyz(t_bc_id,m,2)
            pxyz(m,3) = surf_pxyz(t_bc_id,m,3)
            pxyz(m,4) = surf_pxyz(t_bc_id,m,4)
        end do

        do m=1,3
            check(1) = check(1) + (pxyz(m,1)- pxyz(m,2))**2.0 + (pxyz(m,3)- pxyz(m,4))**2.0
            check(2) = check(2) + (pxyz(m,1)- pxyz(m,3))**2.0 + (pxyz(m,2)- pxyz(m,4))**2.0
        end do

        if(abs(1.0-check(1)**2.0/(dir_check(1)**2.0+1.0e-30))<1.0e-3 .and. &
           abs(1.0-check(2)**2.0/(dir_check(2)**2.0+1.0e-30))<1.0e-3) then
            dir(2)    = 1
            dir(3)    = 3
            do m=1,4
                tp_check(1) = pxyz(1,1) + pxyz(2,1) + pxyz(3,1)
                tp_check(2) = pxyz(1,2) + pxyz(2,2) + pxyz(3,2)
                tp_check(3) = pxyz(1,3) + pxyz(2,3) + pxyz(3,3)
                tp_check(4) = pxyz(1,4) + pxyz(2,4) + pxyz(3,4)
            end do
        else
            dir(2)    = 3
            dir(3)    = 1
            do m=1,4
                tp_check(1) = pxyz(1,1) + pxyz(2,1) + pxyz(3,1)
                tp_check(2) = pxyz(1,3) + pxyz(2,3) + pxyz(3,3)
                tp_check(3) = pxyz(1,2) + pxyz(2,2) + pxyz(3,2)
                tp_check(4) = pxyz(1,4) + pxyz(2,4) + pxyz(3,4)
            end do
        end if

        if(abs(1.0-tp_check(1)**2.0/(p_check(1)**2.0+1.0e-30))<1.0e-3 .and. &
           abs(1.0-tp_check(2)**2.0/(p_check(2)**2.0+1.0e-30))<1.0e-3) then
            direction(1) = 1
        else
            direction(1) = 2
        end if
        if(abs(1.0-tp_check(1)**2.0/(p_check(1)**2.0+1.0e-30))<1.0e-3 .and. &
           abs(1.0-tp_check(3)**2.0/(p_check(3)**2.0+1.0e-30))<1.0e-3) then
            direction(2) = 1
        else
            direction(2) = 2
        end if

        do m=1,2
            if(direction(m) == 2) then
                td = dir(m+1)
                n = ed(td)
                ed(td) = st(td)
                st(td) = n
            end if
        end do

    else
        if (st(3) == 1) then
            lr = -1
        else
            lr = 1
        end if

        do m=1,3
            pxyz(m,1) = surf_pxyz(t_bc_id,m,1)
            pxyz(m,2) = surf_pxyz(t_bc_id,m,2)
            pxyz(m,3) = surf_pxyz(t_bc_id,m,3)
            pxyz(m,4) = surf_pxyz(t_bc_id,m,4)
        end do

        do m=1,3
            check(1) = check(1) + (pxyz(m,1)- pxyz(m,2))**2.0 + (pxyz(m,3)- pxyz(m,4))**2.0
            check(2) = check(2) + (pxyz(m,1)- pxyz(m,3))**2.0 + (pxyz(m,2)- pxyz(m,4))**2.0
        end do

        if(abs(1.0-check(1)**2.0/(dir_check(1)**2.0+1.0e-30))<1.0e-3 .and. &
           abs(1.0-check(2)**2.0/(dir_check(2)**2.0+1.0e-30))<1.0e-3) then
            dir(2)    = 1
            dir(3)    = 2
            do m=1,4
                tp_check(1) = pxyz(1,1) + pxyz(2,1) + pxyz(3,1)
                tp_check(2) = pxyz(1,2) + pxyz(2,2) + pxyz(3,2)
                tp_check(3) = pxyz(1,3) + pxyz(2,3) + pxyz(3,3)
                tp_check(4) = pxyz(1,4) + pxyz(2,4) + pxyz(3,4)
            end do
        else
            dir(2)    = 2
            dir(3)    = 1
            do m=1,4
                tp_check(1) = pxyz(1,1) + pxyz(2,1) + pxyz(3,1)
                tp_check(2) = pxyz(1,3) + pxyz(2,3) + pxyz(3,3)
                tp_check(3) = pxyz(1,2) + pxyz(2,2) + pxyz(3,2)
                tp_check(4) = pxyz(1,4) + pxyz(2,4) + pxyz(3,4)
            end do
        end if

        if(abs(1.0-tp_check(1)**2.0/(p_check(1)**2.0+1.0e-30))<1.0e-3 .and. &
           abs(1.0-tp_check(2)**2.0/(p_check(2)**2.0+1.0e-30))<1.0e-3) then
            direction(1) = 1
        else
            direction(1) = 2
        end if
        if(abs(1.0-tp_check(1)**2.0/(p_check(1)**2.0+1.0e-30))<1.0e-3 .and. &
           abs(1.0-tp_check(3)**2.0/(p_check(3)**2.0+1.0e-30))<1.0e-3) then
            direction(2) = 1
        else
            direction(2) = 2
        end if

        do m=1,2
            if(direction(m) == 2) then
                td = dir(m+1)
                n = ed(td)
                ed(td) = st(td)
                st(td) = n
            end if
        end do
    end if

    do m=1,3
        if (dir(m) == 0) then  !.or. ierr>0
            write(*,'(2(a,i5))') " # BLOCK:",nb,", BC:",nr
            call error_check(m,"Patched boundary is incorrect in subroutine get_bc_dir!")
            exit
        end if
    end do
end subroutine get_bc_dir_patched_t

subroutine get_inter_nr(nbs,nrs,nbt,nrt)
    use mod_kndconsts, only : kind_int
    use mod_constants, only : bc_cut1to1
    use mod_fieldvars, only : mb_top
    implicit none
    integer(kind_int), intent(in)  :: nbs,nrs,nbt
    integer(kind_int), intent(out) :: nrt
    integer(kind_int) :: s_st(3),s_ed(3)
    integer(kind_int) :: t_st(3),t_ed(3)
    integer(kind_int) :: bctype,d(3)
    integer(kind_int) :: nr,m,n,ierr

    s_st(:) = mb_top(nbs)%bcs(nrs)%s_st(:)
    s_ed(:) = mb_top(nbs)%bcs(nrs)%s_ed(:)
    t_st(:) = mb_top(nbs)%bcs(nrs)%t_st(:)
    t_ed(:) = mb_top(nbs)%bcs(nrs)%t_ed(:)

    n = 0
    do nr=1,mb_top(nbt)%nregions
        bctype = mb_top(nbt)%bcs(nr)%bctype
        if (bctype == bc_cut1to1) then
            d(:) =        abs(s_st(:) - mb_top(nbt)%bcs(nr)%t_st(:))
            d(:) = d(:) + abs(s_ed(:) - mb_top(nbt)%bcs(nr)%t_ed(:))
            d(:) = d(:) + abs(t_st(:) - mb_top(nbt)%bcs(nr)%s_st(:))
            d(:) = d(:) + abs(t_ed(:) - mb_top(nbt)%bcs(nr)%s_ed(:))
            m = d(1) + d(2) + d(3)
            if (m == 0) then
               nrt = nr
               n = n + 1
            end if
        end if
    end do

    ierr = n - 1

    call error_check(ierr,"BC file is incorrect in subroutine get_inter_nr!")

end subroutine get_inter_nr

subroutine get_inter_nr_cc(nbs,nrs,nbt,nrt)
    use mod_kndconsts, only : kind_int
    use mod_constants, only : bc_cut1to1
    use mod_fieldvars, only : mb_topc
    implicit none
    integer(kind_int), intent(in)  :: nbs,nrs,nbt
    integer(kind_int), intent(out) :: nrt
    integer(kind_int) :: s_st(3),s_ed(3)
    integer(kind_int) :: t_st(3),t_ed(3)
    integer(kind_int) :: bctype,d(3)
    integer(kind_int) :: nr,m,n,ierr

    s_st(:) = mb_topc(nbs)%bcs(nrs)%s_st(:)
    s_ed(:) = mb_topc(nbs)%bcs(nrs)%s_ed(:)
    t_st(:) = mb_topc(nbs)%bcs(nrs)%t_st(:)
    t_ed(:) = mb_topc(nbs)%bcs(nrs)%t_ed(:)

    n = 0
    do nr=1,mb_topc(nbt)%nregions
        bctype = mb_topc(nbt)%bcs(nr)%bctype
        if (bctype == bc_cut1to1) then
            d(:) =        abs(s_st(:) - mb_topc(nbt)%bcs(nr)%t_st(:))
            d(:) = d(:) + abs(s_ed(:) - mb_topc(nbt)%bcs(nr)%t_ed(:))
            d(:) = d(:) + abs(t_st(:) - mb_topc(nbt)%bcs(nr)%s_st(:))
            d(:) = d(:) + abs(t_ed(:) - mb_topc(nbt)%bcs(nr)%s_ed(:))
            m = d(1) + d(2) + d(3)
            if (m == 0) then
               nrt = nr
               n = n + 1
            end if
        end if
    end do

    ierr = n - 1

    call error_check(ierr,"BC file is incorrect in subroutine get_inter_nr_cc!")

end subroutine get_inter_nr_cc

subroutine get_inter_nr_sp(nbs,nrs,nbt,nrt)
    use mod_kndconsts, only : kind_int
    use mod_constants, only : bc_cut1to1
    use mod_fieldvars, only : mb_topsp
    implicit none
    integer(kind_int), intent(in)  :: nbs,nrs,nbt
    integer(kind_int), intent(out) :: nrt
    integer(kind_int) :: s_st(3),s_ed(3)
    integer(kind_int) :: t_st(3),t_ed(3)
    integer(kind_int) :: bctype,d(3)
    integer(kind_int) :: nr,m,n,ierr

    s_st(:) = mb_topsp(nbs)%bcs(nrs)%s_st(:)
    s_ed(:) = mb_topsp(nbs)%bcs(nrs)%s_ed(:)
    t_st(:) = mb_topsp(nbs)%bcs(nrs)%t_st(:)
    t_ed(:) = mb_topsp(nbs)%bcs(nrs)%t_ed(:)

    n = 0
    do nr=1,mb_topsp(nbt)%nregions
        bctype = mb_topsp(nbt)%bcs(nr)%bctype
        if (bctype == bc_cut1to1) then
            d(:) =        abs(s_st(:) - mb_topsp(nbt)%bcs(nr)%t_st(:))
            d(:) = d(:) + abs(s_ed(:) - mb_topsp(nbt)%bcs(nr)%t_ed(:))
            d(:) = d(:) + abs(t_st(:) - mb_topsp(nbt)%bcs(nr)%s_st(:))
            d(:) = d(:) + abs(t_ed(:) - mb_topsp(nbt)%bcs(nr)%s_ed(:))
            m = d(1) + d(2) + d(3)
            if (m == 0) then
               nrt = nr
               n = n + 1
            end if
        end if
    end do

    ierr = n - 1

    call error_check(ierr,"BC file is incorrect in subroutine get_inter_nr_sp!")

end subroutine get_inter_nr_sp

subroutine set_bc_sequence
    use mod_kndconsts, only : kind_int
    use mod_constants, only : nbclistmax,bclists
    use mod_datatypes, only : bc_region_t,top_block_t
    use mod_fieldvars, only : nblocks,mb_top
    implicit none
    integer(kind_int) :: nb,nl,nr,m,ierr
    integer(kind_int) :: nregs,bctype,bcid
    type(top_block_t), pointer :: top
    type(bc_region_t), pointer :: reg

    do nb=1,nblocks
        top => mb_top(nb)

        nregs = top%nregions

        allocate(top%bcnrs(nregs), stat=ierr)

        m = 0
        do nl=1,nbclistmax
            bcid = bclists(nl)

            do nr=1,nregs
                reg => top%bcs(nr)

                bctype = reg%bctype
                if (bctype == bcid) then
                    m = m + 1
                    top%bcnrs(m) =  nr
                end if
            end do
        end do

        call error_check(m-nregs,"Some BC types are not in our lists")
    end do

end subroutine set_bc_sequence

subroutine set_bc_sequence_cc
    use mod_kndconsts, only : kind_int
    use mod_constants, only : nbclistmax,bclists
    use mod_datatypes, only : bc_region_t,top_block_t
    use mod_fieldvars, only : nblocks,mb_top,mb_topc,mb_topsp
    implicit none
    integer(kind_int) :: nb,nl,nr,m,ierr
    integer(kind_int) :: nregs,bctype,bcid
    type(top_block_t), pointer :: top
    type(bc_region_t), pointer :: reg
    
    do nb=1,nblocks
        top => mb_top(nb)

        nregs = top%nregions

        allocate(top%bcnrs(nregs), stat=ierr)

        m = 0
        do nl=1,nbclistmax
            bcid = bclists(nl)

            do nr=1,nregs
                reg => top%bcs(nr)

                bctype = reg%bctype
                if (bctype == bcid) then
                    m = m + 1
                    top%bcnrs(m) =  nr
                end if
            end do
        end do

        call error_check(m-nregs,"Some BC types are not in our lists")
    end do    

    do nb=1,nblocks
        top => mb_topc(nb)

        nregs = top%nregions

        allocate(top%bcnrs(nregs), stat=ierr)

        m = 0
        do nl=1,nbclistmax
            bcid = bclists(nl)

            do nr=1,nregs
                reg => top%bcs(nr)

                bctype = reg%bctype
                if (bctype == bcid) then
                    m = m + 1
                    top%bcnrs(m) =  nr
                end if
            end do
        end do

        call error_check(m-nregs,"Some BC types are not in our lists")
    end do
    
    do nb=1,nblocks
        top => mb_topsp(nb)

        nregs = top%nregions

        allocate(top%bcnrs(nregs), stat=ierr)

        m = 0
        do nl=1,nbclistmax
            bcid = bclists(nl)

            do nr=1,nregs
                reg => top%bcs(nr)

                bctype = reg%bctype
                if (bctype == bcid) then
                    m = m + 1
                    top%bcnrs(m) =  nr
                end if
            end do
        end do

        call error_check(m-nregs,"Some BC types are not in our lists")
    end do    

end subroutine set_bc_sequence_cc

subroutine set_flag_surface
    use mod_kndconsts, only : kind_int
    use mod_constants, only : nfsf_id_max,nfsf_con_d1der
    use mod_constants, only : nfsf_con_d1int,nfsf_con_d1nint
    use mod_constants, only : nfsf_vis_d2der,nfsf_vis_d2int
    use mod_constants, only : nfsf_vis_d3der,nfsf_vis_d3int
    use mod_constants, only : nfsf_grd_d2der,nfsf_grd_d2int
    use mod_constants, only : nfsf_grd_d3der,nfsf_grd_d3int
    use mod_constants, only : nfsf_bc_policys
    use mod_constants, only : nbc_policy_fullext,nbc_policy_partext
    use mod_constants, only : nbc_policy_partint,nbc_policy_fullint
    use mod_constants, only : fldtype_i3d,bc_wall,bc_symmetry
    use mod_constants, only : bc_cut1to1,subc_cut_spliting,bc_pole
    use mod_constants, only : zero,nscheme_policy_edge2
    use mod_datatypes, only : top_block_t,bc_region_t
    use mod_variables, only : ncutpol,nscheme,twlnd
    use mod_fieldvars, only : nblocks,nblkcoms,blkcoms
    use mod_fieldvars, only : mb_top,mb_fsf
    implicit none
    integer(kind_int)          :: nc,nb,ns,nr,i,j,k,m,n,ierr
    integer(kind_int)          :: st(3),ed(3),nst,ned
    integer(kind_int)          :: nregs,bctype,nd,lr,subtype
    integer(kind_int)          :: nbcpol_wall,nbcpol_symm
    type(top_block_t), pointer :: top
    type(bc_region_t), pointer :: reg

    select case(nscheme)
    case(nscheme_policy_edge2)
        if (twlnd > zero) then
            nbcpol_wall = nbc_policy_fullext
        else
            nbcpol_wall = nbc_policy_partext
        end if
    case default
        nbcpol_wall = nbc_policy_fullext
    end select

    !!nbcpol_wall = nbc_policy_partext
    !!nbcpol_symm = nbc_policy_partext
    nbcpol_symm = nbc_policy_fullext


    allocate(mb_fsf(nblocks,6), stat=ierr)

    nst = 1
    ned = nfsf_id_max

    do nc=1,nblkcoms
        nb  =  blkcoms(nc)%nb
        top => blkcoms(nc)%top

        st(:) = 1
        ed(:) = top%nijk(:)

        call var_array_create(mb_fsf(nb,1),nst,ned, &
                              fldtype_i3d, &
                              st(1),st(1), &
                              st(2),ed(2), &
                              st(3),ed(3))
        call var_array_create(mb_fsf(nb,2),nst,ned, &
                              fldtype_i3d, &
                              ed(1),ed(1), &
                              st(2),ed(2), &
                              st(3),ed(3))

        call var_array_create(mb_fsf(nb,3),nst,ned, &
                              fldtype_i3d, &
                              st(1),ed(1), &
                              st(2),st(2), &
                              st(3),ed(3))
        call var_array_create(mb_fsf(nb,4),nst,ned, &
                              fldtype_i3d, &
                              st(1),ed(1), &
                              ed(2),ed(2), &
                              st(3),ed(3))

        call var_array_create(mb_fsf(nb,5),nst,ned, &
                              fldtype_i3d, &
                              st(1),ed(1), &
                              st(2),ed(2), &
                              st(3),st(3))
        call var_array_create(mb_fsf(nb,6),nst,ned, &
                              fldtype_i3d, &
                              st(1),ed(1), &
                              st(2),ed(2), &
                              ed(3),ed(3))

        nregs = top%nregions
        do ns=1,nregs
            nr  =  top%bcnrs(ns)
            reg => top%bcs(nr)

            bctype  = reg%bctype
            st(:)   = reg%s_st(:)
            ed(:)   = reg%s_ed(:)
            nd      = reg%s_nd
            lr      = reg%s_lr
            subtype = reg%subtype

            m = 2*nd + (lr-1)/2

            select case(bctype)
            case(bc_cut1to1)
                if (subtype == subc_cut_spliting) then
                    do k=st(3),ed(3)
                    do j=st(2),ed(2)
                    do i=st(1),ed(1)
                        do n=nst,ned
                            mb_fsf(nb,m)%fld(n)%i3d(i,j,k) = &
                                nfsf_bc_policys(n,nbc_policy_fullint)
                        end do
                    end do
                    end do
                    end do
                else
                    do k=st(3),ed(3)
                    do j=st(2),ed(2)
                    do i=st(1),ed(1)
                        do n=nst,ned
                            mb_fsf(nb,m)%fld(n)%i3d(i,j,k) = &
                                nfsf_bc_policys(n,ncutpol)
                        end do
                    end do
                    end do
                    end do
                end if
            case(bc_symmetry)
                do k=st(3),ed(3)
                do j=st(2),ed(2)
                do i=st(1),ed(1)
                    do n=nst,ned
                        mb_fsf(nb,m)%fld(n)%i3d(i,j,k) = &
                            nfsf_bc_policys(n,nbcpol_symm)
                    end do
                end do
                end do
                end do
            case(bc_wall)
                do k=st(3),ed(3)
                do j=st(2),ed(2)
                do i=st(1),ed(1)
                    do n=nst,ned
                        mb_fsf(nb,m)%fld(n)%i3d(i,j,k) = &
                            nfsf_bc_policys(n,nbcpol_wall)
                    end do
                end do
                end do
                end do
            case(bc_pole)
                do k=st(3),ed(3)
                do j=st(2),ed(2)
                do i=st(1),ed(1)
                    do n=nst,ned
                        mb_fsf(nb,m)%fld(n)%i3d(i,j,k) = &
                            nfsf_bc_policys(n,nbc_policy_partext)
                    end do
                end do
                end do
                end do
            case default
                do k=st(3),ed(3)
                do j=st(2),ed(2)
                do i=st(1),ed(1)
                    do n=nst,ned
                        mb_fsf(nb,m)%fld(n)%i3d(i,j,k) = &
                            nfsf_bc_policys(n,nbc_policy_fullext)
                    end do
                end do
                end do
                end do
            end select
        end do
    end do

end subroutine set_flag_surface

subroutine set_flag_surface_cc
    use mod_kndconsts, only : kind_int
    use mod_constants, only : nfsf_id_max,nfsf_con_d1der
    use mod_constants, only : nfsf_con_d1int,nfsf_con_d1nint
    use mod_constants, only : nfsf_vis_d2der,nfsf_vis_d2int
    use mod_constants, only : nfsf_vis_d3der,nfsf_vis_d3int
    use mod_constants, only : nfsf_grd_d2der,nfsf_grd_d2int
    use mod_constants, only : nfsf_grd_d3der,nfsf_grd_d3int
    use mod_constants, only : nfsf_bc_policys
    use mod_constants, only : nbc_policy_fullext,nbc_policy_partext
    use mod_constants, only : nbc_policy_partint,nbc_policy_fullint
    use mod_constants, only : fldtype_i3d,fldtype_r3d,bc_wall,bc_symmetry
    use mod_constants, only : bc_cut1to1,bc_farfield,subc_cut_spliting,bc_pole
    use mod_constants, only : zero,nscheme_policy_edge2,bc_inflow,bc_outflow
    use mod_datatypes, only : top_block_t,bc_region_t
    use mod_variables, only : ncutpol,nscheme,twlnd
    use mod_fieldvars, only : nblocks,nblkcoms,blkcoms,blkcomscc,blkcomssp
    use mod_fieldvars, only : mb_topc,mb_fsfcc,mb_topsp,mb_fsfsp,mb_top,mb_fsffp,mb_dpvfp
    implicit none
    integer(kind_int)          :: nc,nb,ns,nr,i,j,k,m,n,ierr
    integer(kind_int)          :: st(3),ed(3),nst,ned
    integer(kind_int)          :: nregs,bctype,nd,lr,subtype
    integer(kind_int)          :: nbcpol_wall,nbcpol_symm
    type(top_block_t), pointer :: top
    type(bc_region_t), pointer :: reg

    select case(nscheme)
    case(nscheme_policy_edge2)
        if (twlnd > zero) then
            nbcpol_wall = nbc_policy_fullext
        else
            nbcpol_wall = nbc_policy_partext
        end if
    case default
        nbcpol_wall = nbc_policy_fullext
    end select

    !!nbcpol_wall = nbc_policy_partext
    !!nbcpol_symm = nbc_policy_partext
    nbcpol_symm = nbc_policy_fullext


    allocate(mb_fsfcc(nblocks,6), stat=ierr)
    allocate(mb_fsfsp(nblocks,6), stat=ierr)
    allocate(mb_fsffp(nblocks,6), stat=ierr)
    allocate(mb_dpvfp(nblocks,6), stat=ierr)

    nst = 1
    ned = nfsf_id_max

    do nc=1,nblkcoms
        nb  =  blkcomscc(nc)%nb
        top => blkcomscc(nc)%top

        st(:) = 1
        ed(:) = top%nijk(:)

        call var_array_create(mb_fsfcc(nb,1),nst,ned, &
                              fldtype_i3d, &
                              st(1),st(1), &
                              st(2),ed(2), &
                              st(3),ed(3))
        call var_array_create(mb_fsfcc(nb,2),nst,ned, &
                              fldtype_i3d, &
                              ed(1),ed(1), &
                              st(2),ed(2), &
                              st(3),ed(3))

        call var_array_create(mb_fsfcc(nb,3),nst,ned, &
                              fldtype_i3d, &
                              st(1),ed(1), &
                              st(2),st(2), &
                              st(3),ed(3))
        call var_array_create(mb_fsfcc(nb,4),nst,ned, &
                              fldtype_i3d, &
                              st(1),ed(1), &
                              ed(2),ed(2), &
                              st(3),ed(3))

        call var_array_create(mb_fsfcc(nb,5),nst,ned, &
                              fldtype_i3d, &
                              st(1),ed(1), &
                              st(2),ed(2), &
                              st(3),st(3))
        call var_array_create(mb_fsfcc(nb,6),nst,ned, &
                              fldtype_i3d, &
                              st(1),ed(1), &
                              st(2),ed(2), &
                              ed(3),ed(3))

        nregs = top%nregions
        do ns=1,nregs
            nr  =  top%bcnrs(ns)
            reg => top%bcs(nr)

            bctype  = reg%bctype
            st(:)   = reg%s_st(:)
            ed(:)   = reg%s_ed(:)
            nd      = reg%s_nd
            lr      = reg%s_lr
            subtype = reg%subtype

            m = 2*nd + (lr-1)/2

            select case(bctype)
            case(bc_cut1to1)
                if (subtype == subc_cut_spliting) then
                    do k=st(3),ed(3)
                    do j=st(2),ed(2)
                    do i=st(1),ed(1)
                        do n=nst,ned
                            mb_fsfcc(nb,m)%fld(n)%i3d(i,j,k) = &
                                nfsf_bc_policys(n,nbc_policy_fullint)
                        end do
                    end do
                    end do
                    end do
                else
                    do k=st(3),ed(3)
                    do j=st(2),ed(2)
                    do i=st(1),ed(1)
                        do n=nst,ned
                            mb_fsfcc(nb,m)%fld(n)%i3d(i,j,k) = &
                                nfsf_bc_policys(n,ncutpol)
                        end do
                    end do
                    end do
                    end do
                end if
            case(bc_symmetry)
                do k=st(3),ed(3)
                do j=st(2),ed(2)
                do i=st(1),ed(1)
                    do n=nst,ned
                        mb_fsfcc(nb,m)%fld(n)%i3d(i,j,k) = &
                            nfsf_bc_policys(n,nbcpol_symm)
                    end do
                end do
                end do
                end do
            case(bc_wall)
                do k=st(3),ed(3)
                do j=st(2),ed(2)
                do i=st(1),ed(1)
                    do n=nst,ned
                        mb_fsfcc(nb,m)%fld(n)%i3d(i,j,k) = &
                            nfsf_bc_policys(n,nbcpol_wall)
                    end do
                end do
                end do
                end do
            case(bc_pole)
                do k=st(3),ed(3)
                do j=st(2),ed(2)
                do i=st(1),ed(1)
                    do n=nst,ned
                        mb_fsfcc(nb,m)%fld(n)%i3d(i,j,k) = &
                            nfsf_bc_policys(n,nbc_policy_partext)
                    end do
                end do
                end do
                end do
            case default
                do k=st(3),ed(3)
                do j=st(2),ed(2)
                do i=st(1),ed(1)
                    do n=nst,ned
                        mb_fsfcc(nb,m)%fld(n)%i3d(i,j,k) = &
                            nfsf_bc_policys(n,nbc_policy_fullext)
                    end do
                end do
                end do
                end do
            end select
        end do
    end do
    
    do nc=1,nblkcoms
        nb  =  blkcomssp(nc)%nb
        top => blkcomssp(nc)%top

        st(:) = 1
        ed(:) = top%nijk(:)

        call var_array_create(mb_fsfsp(nb,1),nst,ned, &
                              fldtype_i3d, &
                              st(1),st(1), &
                              st(2),ed(2), &
                              st(3),ed(3))
        call var_array_create(mb_fsfsp(nb,2),nst,ned, &
                              fldtype_i3d, &
                              ed(1),ed(1), &
                              st(2),ed(2), &
                              st(3),ed(3))

        call var_array_create(mb_fsfsp(nb,3),nst,ned, &
                              fldtype_i3d, &
                              st(1),ed(1), &
                              st(2),st(2), &
                              st(3),ed(3))
        call var_array_create(mb_fsfsp(nb,4),nst,ned, &
                              fldtype_i3d, &
                              st(1),ed(1), &
                              ed(2),ed(2), &
                              st(3),ed(3))

        call var_array_create(mb_fsfsp(nb,5),nst,ned, &
                              fldtype_i3d, &
                              st(1),ed(1), &
                              st(2),ed(2), &
                              st(3),st(3))
        call var_array_create(mb_fsfsp(nb,6),nst,ned, &
                              fldtype_i3d, &
                              st(1),ed(1), &
                              st(2),ed(2), &
                              ed(3),ed(3))

        nregs = top%nregions
        do ns=1,nregs
            nr  =  top%bcnrs(ns)
            reg => top%bcs(nr)

            bctype  = reg%bctype
            st(:)   = reg%s_st(:)
            ed(:)   = reg%s_ed(:)
            nd      = reg%s_nd
            lr      = reg%s_lr
            subtype = reg%subtype

            m = 2*nd + (lr-1)/2

            select case(bctype)
            case(bc_cut1to1)
                if (subtype == subc_cut_spliting) then
                    do k=st(3),ed(3)
                    do j=st(2),ed(2)
                    do i=st(1),ed(1)
                        do n=nst,ned
                            mb_fsfsp(nb,m)%fld(n)%i3d(i,j,k) = &
                                nfsf_bc_policys(n,nbc_policy_fullint)
                        end do
                    end do
                    end do
                    end do
                else
                    do k=st(3),ed(3)
                    do j=st(2),ed(2)
                    do i=st(1),ed(1)
                        do n=nst,ned
                            mb_fsfsp(nb,m)%fld(n)%i3d(i,j,k) = &
                                nfsf_bc_policys(n,ncutpol)
                        end do
                    end do
                    end do
                    end do
                end if
            case(bc_symmetry)
                do k=st(3),ed(3)
                do j=st(2),ed(2)
                do i=st(1),ed(1)
                    do n=nst,ned
                        mb_fsfsp(nb,m)%fld(n)%i3d(i,j,k) = &
                            nfsf_bc_policys(n,nbcpol_symm)
                    end do
                end do
                end do
                end do
            case(bc_wall)
                do k=st(3),ed(3)
                do j=st(2),ed(2)
                do i=st(1),ed(1)
                    do n=nst,ned
                        mb_fsfsp(nb,m)%fld(n)%i3d(i,j,k) = &
                            nfsf_bc_policys(n,nbcpol_wall)
                    end do
                end do
                end do
                end do
            case(bc_pole)
                do k=st(3),ed(3)
                do j=st(2),ed(2)
                do i=st(1),ed(1)
                    do n=nst,ned
                        mb_fsfsp(nb,m)%fld(n)%i3d(i,j,k) = &
                            nfsf_bc_policys(n,nbc_policy_partext)
                    end do
                end do
                end do
                end do
            case default
                do k=st(3),ed(3)
                do j=st(2),ed(2)
                do i=st(1),ed(1)
                    do n=nst,ned
                        mb_fsfsp(nb,m)%fld(n)%i3d(i,j,k) = &
                            nfsf_bc_policys(n,nbc_policy_fullext)
                    end do
                end do
                end do
                end do
            end select
        end do
    end do
    
    do nc=1,nblkcoms
        nb  =  blkcoms(nc)%nb
        top => blkcoms(nc)%top

        st(:) = 1
        ed(:) = top%nijk(:)

        call var_array_create(mb_fsffp(nb,1),1,2, &
                              fldtype_i3d, &
                              st(1),st(1), &
                              st(2),ed(2), &
                              st(3),ed(3))
        call var_array_create(mb_fsffp(nb,2),1,2, &
                              fldtype_i3d, &
                              ed(1),ed(1), &
                              st(2),ed(2), &
                              st(3),ed(3))

        call var_array_create(mb_fsffp(nb,3),1,2, &
                              fldtype_i3d, &
                              st(1),ed(1), &
                              st(2),st(2), &
                              st(3),ed(3))
        call var_array_create(mb_fsffp(nb,4),1,2, &
                              fldtype_i3d, &
                              st(1),ed(1), &
                              ed(2),ed(2), &
                              st(3),ed(3))

        call var_array_create(mb_fsffp(nb,5),1,2, &
                              fldtype_i3d, &
                              st(1),ed(1), &
                              st(2),ed(2), &
                              st(3),st(3))
        call var_array_create(mb_fsffp(nb,6),1,2, &
                              fldtype_i3d, &
                              st(1),ed(1), &
                              st(2),ed(2), &
                              ed(3),ed(3))
        
        call var_array_create(mb_dpvfp(nb,1),1,12, &
                              fldtype_r3d, &
                              st(1),st(1), &
                              st(2),ed(2), &
                              st(3),ed(3))
        call var_array_create(mb_dpvfp(nb,2),1,12, &
                              fldtype_r3d, &
                              ed(1),ed(1), &
                              st(2),ed(2), &
                              st(3),ed(3))

        call var_array_create(mb_dpvfp(nb,3),1,12, &
                              fldtype_r3d, &
                              st(1),ed(1), &
                              st(2),st(2), &
                              st(3),ed(3))
        call var_array_create(mb_dpvfp(nb,4),1,12, &
                              fldtype_r3d, &
                              st(1),ed(1), &
                              ed(2),ed(2), &
                              st(3),ed(3))

        call var_array_create(mb_dpvfp(nb,5),1,12, &
                              fldtype_r3d, &
                              st(1),ed(1), &
                              st(2),ed(2), &
                              st(3),st(3))
        call var_array_create(mb_dpvfp(nb,6),1,12, &
                              fldtype_r3d, &
                              st(1),ed(1), &
                              st(2),ed(2), &
                              ed(3),ed(3))        

        nregs = top%nregions
        do ns=1,nregs
            nr  =  top%bcnrs(ns)
            reg => top%bcs(nr)

            bctype  = reg%bctype
            st(:)   = reg%s_st(:)
            ed(:)   = reg%s_ed(:)
            nd      = reg%s_nd
            lr      = reg%s_lr

            m = 2*nd + (lr-1)/2
            
            do k=st(3),ed(3)
            do j=st(2),ed(2)
            do i=st(1),ed(1)
                do n=1,12
                    mb_dpvfp(nb,m)%fld(n)%r3d(i,j,k) = 0.0
                end do
            end do
            end do
            end do

            select case(bctype)
            case(bc_cut1to1)
                do k=st(3),ed(3)
                do j=st(2),ed(2)
                do i=st(1),ed(1)
                    mb_fsffp(nb,m)%fld(1)%i3d(i,j,k) = nr
                    mb_fsffp(nb,m)%fld(2)%i3d(i,j,k) = bc_cut1to1
                end do
                end do
                end do
            case(bc_wall)
                do k=st(3),ed(3)
                do j=st(2),ed(2)
                do i=st(1),ed(1)
                    mb_fsffp(nb,m)%fld(1)%i3d(i,j,k) = nr
                    mb_fsffp(nb,m)%fld(2)%i3d(i,j,k) = bc_wall
                end do
                end do
                end do
            case(bc_symmetry)
                do k=st(3),ed(3)
                do j=st(2),ed(2)
                do i=st(1),ed(1)
                    mb_fsffp(nb,m)%fld(1)%i3d(i,j,k) = nr
                    mb_fsffp(nb,m)%fld(2)%i3d(i,j,k) = bc_symmetry
                end do
                end do
                end do                
            case(bc_farfield)
                do k=st(3),ed(3)
                do j=st(2),ed(2)
                do i=st(1),ed(1)
                    mb_fsffp(nb,m)%fld(1)%i3d(i,j,k) = nr
                    mb_fsffp(nb,m)%fld(2)%i3d(i,j,k) = bc_farfield
                end do
                end do
                end do
            case(bc_inflow)
                do k=st(3),ed(3)
                do j=st(2),ed(2)
                do i=st(1),ed(1)
                    mb_fsffp(nb,m)%fld(1)%i3d(i,j,k) = nr
                    mb_fsffp(nb,m)%fld(2)%i3d(i,j,k) = bc_inflow
                end do
                end do
                end do
            case(bc_outflow)
                do k=st(3),ed(3)
                do j=st(2),ed(2)
                do i=st(1),ed(1)
                    mb_fsffp(nb,m)%fld(1)%i3d(i,j,k) = nr
                    mb_fsffp(nb,m)%fld(2)%i3d(i,j,k) = bc_outflow
                end do
                end do
                end do                
            case(bc_pole)
                do k=st(3),ed(3)
                do j=st(2),ed(2)
                do i=st(1),ed(1)
                    mb_fsffp(nb,m)%fld(1)%i3d(i,j,k) = nr
                    mb_fsffp(nb,m)%fld(2)%i3d(i,j,k) = bc_pole
                end do
                end do
                end do
            case default
                do k=st(3),ed(3)
                do j=st(2),ed(2)
                do i=st(1),ed(1)
                    mb_fsffp(nb,m)%fld(1)%i3d(i,j,k) = 0
                    mb_fsffp(nb,m)%fld(2)%i3d(i,j,k) = 0
                end do
                end do
                end do
            end select
        end do
    end do    

end subroutine set_flag_surface_cc

!>
!! @brief ����ȫ�������������ӵ㣨�����㣩�ϵ�������Ϣ
!! @author ������
!! @remark 2012��04��01�� �������ṹ������㣬��ʱ�俪�������ɺ��Բ���
!>
subroutine search_singulars()
    use mod_kndconsts, only : kind_int
    use mod_fieldvars, only : nblocks,nprocs
    use mod_singulars
    use mod_parallels
    implicit none

    integer(kind_int)  :: nvars,nmsg
    integer(kind_int)  :: i,j,m,nbs,ierr
    integer(kind_int)  :: nsgl_buffers
    logical            :: changed
    type(box_node_ptr_t),target  :: box_group(nblocks)
    type(box_node_ptr_t),pointer :: bgrp_root(:)

    bgrp_root => box_group

    ! step 1: �����ʼ�߼�
    call build_initial_segment_set(bgrp_root)

    ! step 2: ͨ�����������������ռ����й����ı�
    changed = .true.
    do while (changed)
        call divide_into_disjoint_sets(bgrp_root, changed)
        call collect_all_singular_segment(bgrp_root, changed)
    end do

    ! step 3: �����ȼ۹�ϵ
    call build_equiv_relation(bgrp_root)

    ! step 4: ɸѡ���˵�����Ҫ�ĵȼۼ������˻���������2�ĵȼۼ���
    do nbs = 1, nblocks
        call tree_node_filter(box_group(nbs)%ptr)
    end do

    ! step 5: �γɵ���ȼۼ���
    call build_segment_set(bgrp_root)

    do nbs=1,nblocks
        call tree_destroy(box_group(nbs)%ptr)
    end do

    ! step 6: ��� nsimp_pnts_local, �Լ� simp_arr �� gids ��

    allocate(nmsg_proc(0:nprocs-1, 3), stat=ierr)
    nmsg_proc(:, :) = 0

    call build_nsimp_pnts_local()

    ! step 7: ���simp_arr �� fsw ��
    call fill_simp_fsw_field()

    ! step 8: ׼����Ϣӳ���, ��Ϣ������
    call build_msg_table()

contains

    subroutine build_initial_segment_set(box_group)
        use mod_kndconsts, only : kind_int
        use mod_fieldvars, only : nblocks,nprocs,mb_top
        use mod_datatypes, only : top_block_t,bc_region_t
        use mod_constants, only : bc_cut1to1
        use mod_singulars
        use mod_parallels
        implicit none
        type(box_node_ptr_t), pointer, intent(inout) :: box_group(:)

        integer(kind_int)           :: nbs, nrs, n, ierr
        logical                     :: succ
        type(top_block_t), pointer  :: top
        type(bc_region_t), pointer  :: reg
        type(box_node_t), pointer   :: bn_tree
        type(box_t)                 :: init_edge_list(12)
        type(box_t)                 :: box_a, box
        type(box_node_t), pointer   :: box_node

        do nbs=1,nblocks
            top => mb_top(nbs)

            ! ������������12����
            call edge_box_create(init_edge_list, top%nijk)
            nullify(bn_tree)
            do n = 1, 12
            do nrs=1,top%nregions
                reg => top%bcs(nrs)

                if ( reg%bctype == bc_cut1to1 ) then
                    call box_create( box_a, reg%s_st, reg%s_ed )
                    call box_intersection( init_edge_list(n), box_a, box )
                    if ( .not. box_is_empty(box) ) then
                        if (box_dim(box) > 1) then
                            write(*, *) "!!! warning: box_dim > 1:"
                        end if
                        call box_node_create(box_node, box, nbs)
                        call tree_insert(bn_tree, box_node, stat=succ)
                        call box_node_destroy(box_node, .not. succ)
                    end if
                end if
            end do ! nrs=1,nregs
            end do ! edge = 1, 12

            box_group(nbs)%ptr => bn_tree
            box_group(nbs)%flag = .true. ! flag = changed
        end do !  nbs=1,nblocks

    end subroutine build_initial_segment_set

    subroutine build_msg_table()
        use mod_kndconsts, only : kind_int
        use mod_constants, only : zero,one,mcyc
        use mod_variables, only : nvis
        use mod_constants, only : nvis_ns_lam
        use mod_constants, only : nsgl_buffer_max,nsgl_buffer_pvs
        use mod_constants, only : nsgl_buffer_qts,nsgl_buffer_dqt
        use mod_constants, only : nsgl_buffer_dpv,nsgl_buffer_dtur
        use mod_fieldvars, only : neqn,neqt,nqvst
        use mod_fieldvars, only : nblocks,nprocs,mb_top
        use mod_singulars
        use mod_parallels
#ifdef PARALLEL
        use mpi
#endif
        implicit none

        type(item_list_t), pointer :: simp_this
        integer(kind_int)  :: npnts, nbs, nbt, pds, pdt, pid
        integer(kind_int)  :: n
        integer(kind_int)  :: nbufpars(1,nsgl_buffer_max)

#ifdef PARALLEL
    pid = myid
#else
    pid = 0
#endif

        loop_group: do n = 1, nsimp_arr
            simp_this => simp_arr(n)
            npnts     = simp_this%npnts

            loop_this: do i = 1, npnts
                nbs = simp_this%pnts(1, i)
                pds = mb_top(nbs)%pid
                if (pds == pid) then
                    do j = 1, npnts
                        if (i == j) cycle
                        nbt = simp_this%pnts(1, j)
                        pdt = mb_top(nbt)%pid

                        nmsg_proc(pdt, 2) = 1
                    end do
                    cycle loop_group
                end if
            end do loop_this
        end do loop_group


#ifdef PARALLEL
        nmsg = sum(nmsg_proc(:, 2))
        nmsg_sin = ( nmsg - nmsg_proc(myid, 2) ) * 2
        if (nmsg_sin > 0) then
            allocate(request_sin(nmsg_sin), stat=ierr)
            allocate(status_sin(MPI_STATUS_SIZE, nmsg_sin), stat=ierr)
        end if
#endif

        allocate(simp_buf(1:nsgl_buffer_max, 0:nprocs-1), stat=ierr)
        nbufpars(:,nsgl_buffer_pvs ) = (/neqn/)
        nbufpars(:,nsgl_buffer_dpv ) = (/4*3/)
        nbufpars(:,nsgl_buffer_qts ) = (/nqvst/)
        nbufpars(:,nsgl_buffer_dtur) = (/neqt*3/)
        nbufpars(:,nsgl_buffer_dqt ) = (/neqt/)

        if (nvis > nvis_ns_lam) then
            nsgl_buffers = nsgl_buffer_max
        else
            nsgl_buffers = nsgl_buffer_dpv
        end if

        do m = 1, nsgl_buffers
            nvars = nbufpars(1,m)
            do i = 0, nprocs-1
                npnts = nmsg_proc(i, 1)

                simp_buf(m, i)%nvar = nvars
                simp_buf(m, i)%npnt = npnts
                nullify(simp_buf(m, i)%buf)

!                if (npnts > 0) then
                    if ( (nmsg_proc(i, 2) > 0) .or. (pid == i)) then
                        allocate(simp_buf(m, i)%buf(nvars, npnts), stat=ierr)
                    end if
!                end if
            end do
        end do

    end subroutine build_msg_table

    subroutine divide_into_disjoint_sets(box_group, changed)
        use mod_kndconsts, only : kind_int
        use mod_fieldvars, only : nblocks
        use mod_singulars
        implicit none
        type(box_node_ptr_t), pointer, intent(inout) :: box_group(:)
        logical, intent(out) :: changed

        type(box_node_t), pointer :: tree_root

        type(box_node_ptr_t), target  :: bnp_head
        type(box_node_ptr_t), pointer :: bnp_tail
        type(box_node_ptr_t), pointer :: bnp_this, bnp_that, bnp_temp
        type(box_node_t),     pointer :: bn_this, bn_that

        type(box_node_t), pointer   :: new_tree

        type(box_t) :: box_this, box_that
        type(box_t) :: outbox(5)    ! �����߶����ֽ��5�������Ӽ�
        integer(kind_int) :: nbox
        integer(kind_int) :: nbs, n, ierr

        logical :: block_changed
        logical :: succ

        changed = .false.
        loop_block : do nbs=1,nblocks
            tree_root => box_group(nbs)%ptr

            if (.not. box_group(nbs)%flag)   cycle loop_block

            if (.not. associated(tree_root)) then
                write(*, *) "!!! error: box node tree is NULL in block ", nbs
                cycle loop_block
            end if
            if (nbs /= tree_root%nblk) then
                write(*, *) "   !!! error: nbs /= tree_root%nblk [nbs, nblk]: ", nbs, tree_root%nblk
            end if
            block_changed = .false.

            ! ��������ѭ������, bnp_head ��ָ�����õ�ͷ���, bnp_tail ��ָ��������β���.
            ! ���б�Ϊ��ʱ, head = tail = head%next = tail%next
            bnp_head%next => bnp_head
            bnp_tail => bnp_head
            call list_append(bnp_tail, tree_root)

            ! ������������������������������Ƿ��ཻ
            bnp_this => bnp_head
            loop_this: do while (.not. associated(bnp_this%next, bnp_head))
                if ( .not. bnp_this%next%flag ) then
                    call box_copy(box_this, bnp_this%next%ptr%box)

                    bnp_that => bnp_this%next
                    loop_that: do while (.not. associated(bnp_that%next, bnp_head))
                        succ = .not. bnp_that%next%flag

                        if (succ) then
                            call box_copy(box_that, bnp_that%next%ptr%box)
                            call box_divide(box_this, box_that, outbox, nbox)
                            if (nbox > 0) then

                                block_changed = .true.
                                changed = .true.

                                ! ɾ�� this, that, �������� box ����β
                                bnp_this%next%flag = .true.
                                bnp_that%next%flag = .true.
                                do n = 1, nbox
                                    nullify(bn_this)
                                    call box_node_create(bn_this, outbox(n), nbs)
                                    call tree_insert(tree_root, bn_this, stat=succ, pos=bn_that)
                                    call box_node_destroy(bn_this, .not. succ)
                                    if (succ .or. &
                                        associated(bn_that, bnp_this%next%ptr) .or. &
                                        associated(bn_that, bnp_that%next%ptr)) then

                                        nullify(bnp_temp)
                                        allocate(bnp_temp, stat=ierr)
                                        bnp_temp%flag = .false.
                                        bnp_temp%ptr => bn_that
                                        bnp_temp%next => bnp_tail%next
                                        bnp_tail%next => bnp_temp
                                        bnp_tail => bnp_tail%next

                                    end if
                                end do

                                exit loop_that
                            end if
                        end if

                        bnp_that => bnp_that%next
                    end do loop_that
                end if

                bnp_this => bnp_this%next
            end do loop_this

            if (block_changed) then
                ! �½� box_node ��������������н��以���ཻ
                nullify(new_tree)
                bnp_this => bnp_head
                do while (.not. associated(bnp_this%next, bnp_head))
                    if ( .not. bnp_this%next%flag ) then
                        nullify(bn_this)
                        call box_node_create(bn_this, bnp_this%next%ptr%box, nbs)
                        call tree_insert(new_tree, bn_this, stat=succ)
                        call box_node_destroy(bn_this, .not. succ)
                    end if
                    bnp_this => bnp_this%next
                end do

                call tree_destroy(tree_root)

                box_group(nbs)%ptr => new_tree
            else
            end if

            box_group(nbs)%flag = block_changed

            ! �ͷŴ洢
            call list_destroy(bnp_head)

        end do loop_block
    end subroutine divide_into_disjoint_sets

    subroutine collect_all_singular_segment(box_group, changed)
        use mod_kndconsts, only : kind_int
        use mod_constants, only : bc_cut1to1,bc_wall
        use mod_fieldvars, only : nblocks,mb_top
        use mod_datatypes, only : top_block_t,bc_region_t
        use mod_singulars
        implicit none
        type(box_node_ptr_t), pointer, intent(inout) :: box_group(:)
        logical, intent(out) :: changed

        type(top_block_t), pointer  :: top
        type(bc_region_t), pointer  :: reg
        integer(kind_int)           :: bctype, nbs, nbt, nrs, ierr
        integer(kind_int)           :: s_st(3),s_ed(3)
        integer(kind_int)           :: d_st(3),d_ed(3)
        logical :: succ

        type(box_t)                 :: box_this, box_that
        type(box_node_t), pointer   :: bn_this
        type(box_node_ptr_t),target :: bnp_head
        type(box_node_ptr_t),pointer:: bnp_tail
        type(box_node_ptr_t),pointer:: bnp_temp

        changed = .false.
        bnp_head%next => bnp_head
        bnp_tail => bnp_head
        do nbs=1,nblocks
            if (box_group(nbs)%flag) then
                call list_append(bnp_tail, box_group(nbs)%ptr)
                box_group(nbs)%flag = .false.
            end if
        end do

        do while (.not. associated(bnp_head%next, bnp_head))
            ! ȡ����ͷ����box(�������ͷſռ�)
            bnp_temp      => bnp_head%next
            bnp_head%next => bnp_temp%next

            bn_this => bnp_temp%ptr
            nbs = bn_this%nblk
            call box_copy(box_this, bn_this%box)
            s_st = box_this%st
            s_ed = box_this%ed
            nullify(bn_this)
            deallocate(bnp_temp, stat=ierr)

            ! ���� box �������ھӣ��ǽ�����һһ����ϲ�������

            top => mb_top(nbs)
            do nrs = 1, top%nregions
                reg => top%bcs(nrs)
                nbt = reg%nbt
                bctype = reg%bctype

                if (bctype == bc_cut1to1) then
                    call box_create(box_that, reg%s_st, reg%s_ed)
                    if (box_contains_box(box_that, box_this)) then
                        d_st = reg%mapijk(s_st(1), s_st(2), s_st(3), :)
                        d_ed = reg%mapijk(s_ed(1), s_ed(2), s_ed(3), :)
                        call box_make(box_that, d_st, d_ed)

                        call box_node_create(bn_this, box_that, nbt)
                        call tree_insert(box_group(nbt)%ptr, bn_this, stat=succ)
                        call box_node_destroy(bn_this, .not. succ)
                        if (succ) then
                            box_group(nbt)%flag = .true.
                            changed = .true.
                        end if
                    end if
                end if
            end do ! nrs = 1, nregs

        end do
        ! ���ˣ��Ѿ����������е����������Ϣ�������˼ӷ������Ժ�ֻ��Ҫ������

    end subroutine collect_all_singular_segment

    subroutine build_equiv_relation(box_group)
        use mod_kndconsts, only : kind_int
        use mod_constants, only : bc_cut1to1,bc_wall
        use mod_fieldvars, only : nblocks,mb_top
        use mod_datatypes, only : top_block_t,bc_region_t
        use mod_singulars
        implicit none
        type(box_node_ptr_t), pointer, intent(inout) :: box_group(:)

        type(top_block_t), pointer  :: top
        type(bc_region_t), pointer  :: reg
        integer(kind_int)           :: bctype, nbs, nbt, nrs
        integer(kind_int)           :: dir, dir_desired
        integer(kind_int)           :: s_st(3),s_ed(3)
        integer(kind_int)           :: d_st(3),d_ed(3)
        logical :: succ, found

        type(box_t)                 :: box_this, box_that
        type(box_node_t), pointer   :: bn_this, bn_that, bn_pos
        type(box_node_ptr_t),target :: bnp_head
        type(box_node_ptr_t),pointer:: bnp_tail
        type(box_node_ptr_t),pointer:: bnp_this

        bnp_head%next => bnp_head
        bnp_tail => bnp_head
        do nbs=1,nblocks
            call list_append(bnp_tail, box_group(nbs)%ptr)
        end do

        bnp_this => bnp_head
        do while (.not. associated(bnp_this%next, bnp_head))
            bn_this => bnp_this%next%ptr
            nbs = bn_this%nblk
            dir = bn_this%dir
            call box_copy(box_this, bn_this%box)
            s_st = box_this%st
            s_ed = box_this%ed

            ! ���� box �������ھӣ��ǽ�����һһ����ϲ�������

            top => mb_top(nbs)
            loop_that: do nrs = 1, top%nregions
                reg => top%bcs(nrs)
                bctype = reg%bctype

                if (bctype == bc_cut1to1) then
                    nbt = reg%nbt

                    call box_create(box_that, reg%s_st, reg%s_ed)
                    if (box_contains_box(box_that, box_this)) then

                        d_st = reg%mapijk(s_st(1), s_st(2), s_st(3), :)
                        d_ed = reg%mapijk(s_ed(1), s_ed(2), s_ed(3), :)
                        call box_make(box_that, d_st, d_ed)
                        if (point_equal(box_that%st, d_st)) then
                            dir_desired = dir
                        else if (point_equal(box_that%st, d_ed)) then
                            dir_desired = - dir
                        else
                            write(*,*) "!! error: impossible situation!"
                        end if

                        call box_node_create(bn_that, box_that, nbt)
                        call tree_insert(box_group(nbt)%ptr, bn_that, stat=succ, pos=bn_pos)
                        if (.not. succ) then ! ���Ѵ���������, ���һ���ж����� bn_this �Ƿ�ȼ�
                            call box_node_destroy(bn_that)
                            bn_that => bn_pos
                            call box_node_chain_search(bn_that, bn_this, found)

                            if (.not. found) then
                                ! ���Ѵ��ڵķ����Ƿ��롰�ڴ��ķ�����ͬ? ��ͬ��ת
                                if (bn_that%dir /= dir_desired) then
                                    call box_node_chain_flip_dir(bn_that)
                                end if

                                ! �ϲ� bn_this �� bn_that �����ȼۼ�
                                bn_pos            => bn_this%next
                                bn_this%next      => bn_that
                                bn_pos%prev       => bn_that%prev
                                bn_that%prev%next => bn_pos
                                bn_that%prev      => bn_this
                            else
                                if (bn_that%dir /= dir_desired) then
                                    write(*,*) "!! warn: direction conflicked in a same chain!"
                                end if
                            end if
                        else
                            write(*,*) "!! error: impossible situation! [find a new neighboring seg]"
                        end if
                    end if
                end if
            end do loop_that

            bnp_this => bnp_this%next
        end do

        call list_destroy(bnp_head)
    end subroutine build_equiv_relation

    subroutine build_segment_set(box_group)
        use mod_kndconsts, only : kind_int
        use mod_fieldvars, only : nblocks,mb_top
        use mod_singulars
        implicit none
        type(box_node_ptr_t), pointer, intent(inout) :: box_group(:)

        type(box_node_ptr_t),target :: bnp_head
        type(box_node_ptr_t),pointer:: bnp_tail
        type(box_node_ptr_t),pointer:: bnp_this

        integer(kind_int) :: nsegs_total, npnts_total

        ! ���켯���б�
        bnp_head%next   => bnp_head
        bnp_tail        => bnp_head
        do nbs=1, nblocks
            call list_append_by_set_member(bnp_tail, box_group(nbs)%ptr)
        end do

        ! ͳ�Ƶ㼯�������������
        nsegs_total = 0
        npnts_total = 0
        bnp_this => bnp_head%next
        do while (.not. associated(bnp_this, bnp_head))
            nsegs_total = nsegs_total + 1
            npnts_total = npnts_total + bnp_this%nlen_segmt

            bnp_this => bnp_this%next
        end do

        ! ����ȫ������ sim_top_tmp (������ÿһ��Ϊһ���㼯��)
        call build_simp_array(bnp_head, npnts_total)
    end subroutine build_segment_set

    subroutine build_simp_array(bnp_head, npnts_total)
        use mod_kndconsts, only : kind_int
        use mod_singulars
        implicit none
        type(box_node_ptr_t), intent(in), target :: bnp_head
        integer(kind_int),    intent(in) :: npnts_total

        type(box_node_ptr_t), pointer :: bnp_this
        type(box_node_t),     pointer :: bn_this
        integer(kind_int) :: pnt(3)
        integer(kind_int) :: nlen_chain, nlen_segmt, igroup
        integer(kind_int) :: i, j, ierr

        type(item_list_t), pointer :: item_this

        allocate(simp_arr(npnts_total), stat=ierr)
        nsimp_arr = npnts_total

        igroup = 0
        bnp_this => bnp_head%next
        do while (.not. associated(bnp_this, bnp_head))
            nlen_segmt   = bnp_this%nlen_segmt
            nlen_chain   = bnp_this%nlen_chain

            do j = 1, nlen_segmt
                igroup = igroup + 1
                item_this => simp_arr(igroup)

                item_this%npnts = nlen_chain
                allocate(item_this%pnts(6, nlen_chain), stat=ierr)

                bn_this => bnp_this%ptr
                do i = 1, nlen_chain
                    call box_node_get_point(bn_this, j, pnt)
                    item_this%pnts(1,   i) = bn_this%nblk
                    item_this%pnts(2:4, i) = pnt(1:3)

                    bn_this => bn_this%next
                end do

            end do

            bnp_this => bnp_this%next
        end do

    end subroutine build_simp_array

    subroutine build_nsimp_pnts_local()
        use mod_kndconsts, only : kind_int
        use mod_fieldvars, only : nblocks,nprocs,mb_top
        use mod_singulars
        implicit none
        integer(kind_int) :: i, j, npnts, nbs, pid


        do i = 1, nsimp_arr
            npnts = simp_arr(i)%npnts
            do j = 1, npnts
                nbs = simp_arr(i)%pnts(1, j)
                pid = mb_top(nbs)%pid
                nmsg_proc(pid, 1) = nmsg_proc(pid, 1) + 1
                simp_arr(i)%pnts(5, j) = nmsg_proc(pid, 1)
            end do
        end do
    end subroutine build_nsimp_pnts_local

    subroutine fill_simp_fsw_field()
        use mod_kndconsts, only : kind_int
        use mod_fieldvars, only : nblocks,nprocs,mb_top
        use mod_datatypes, only : top_block_t,bc_region_t
        use mod_constants, only : bc_wall
        use mod_singulars
        implicit none

        type(top_block_t), pointer :: top
        type(bc_region_t), pointer :: reg
        type(box_t) :: box
        integer(kind_int) :: i, j, npnts, nbs, nrs, bctype
        logical :: succ
        integer(kind_int) :: pnt(3)

        do i = 1, nsimp_arr    !! ����������ԭ��(rhsҲ��Ҫ�𣿣���)
            npnts = simp_arr(i)%npnts

            succ = .false.  ! ��ʾû���ĸ���λ��������
            simp_arr(i)%pnts(6, :) = 0
            do j = 1, npnts
                nbs = simp_arr(i)%pnts(1, j)
                pnt = simp_arr(i)%pnts(2:4, j)
                top => mb_top(nbs)
                loop_all_bc: do nrs=1,top%nregions
                    reg => top%bcs(nrs)
                    bctype = reg%bctype
                    if ( bctype == bc_wall ) then
                        call box_create(box, reg%s_st, reg%s_ed)
                        if ( box_contains_pnt(box, pnt) ) then
                            succ = .true.
                            simp_arr(i)%pnts(6, j) = 1
                            exit loop_all_bc
                        end if
                    end if
                end do loop_all_bc
            end do

            if (succ) then
                ! do nothing
            else
                simp_arr(i)%pnts(6, :) = 1
            end if
        end do
    end subroutine fill_simp_fsw_field

    recursive   &
    subroutine  list_append_by_set_member(tail, subtree)
        use mod_kndconsts, only : kind_int
        use mod_singulars
        implicit none
        type(box_node_ptr_t), pointer :: tail
        type(box_node_t),     pointer, intent(in)    :: subtree

        type(box_node_ptr_t), pointer :: node
        integer :: ierr

        if ( associated(subtree) ) then
            call list_append_by_set_member(tail, subtree%lnode)

            if ( .not. (subtree%visited .or. subtree%del) ) then
                allocate(node, stat=ierr)
                node%ptr  => subtree
                node%next => tail%next
                node%nlen_segmt = box_volume(subtree%box)
                node%nlen_chain = box_node_chain_count(subtree)
                node%flag = .false.
                tail%next => node

                tail => tail%next
                call box_node_chain_set_visited(subtree)
            end if

            call list_append_by_set_member(tail, subtree%rnode)
        end if
    end subroutine list_append_by_set_member

end subroutine search_singulars

subroutine grid_derivative
    use mod_constants, only : nintplt_node2,nintplt_node4
    use mod_constants, only : nintplt_node6,nintplt_node6w
    use mod_constants, only : nintplt_node8
    use mod_constants, only : nintplt_scsl4,nintplt_scsl6
    use mod_variables, only : nd3int_grd
    implicit none
    external :: ve_via_node2,ve_via_node4
    external :: ve_via_node6,ve_via_node6w,ve_via_node8
    external :: ve_via_scsl4,ve_via_scsl6

    select case(nd3int_grd)
    case(nintplt_node2)
        call grid_derivative_dergrd(ve_via_node2)
    case(nintplt_node4)
        call grid_derivative_dergrd(ve_via_node4)
    case(nintplt_node6)
        call grid_derivative_dergrd(ve_via_node6)
    case(nintplt_node6w)
        call grid_derivative_dergrd(ve_via_node6w)
    case(nintplt_node8)
        call grid_derivative_dergrd(ve_via_node8)
    case(nintplt_scsl4)
        call grid_derivative_dergrd(ve_via_scsl4)
    case(nintplt_scsl6)
        call grid_derivative_dergrd(ve_via_scsl6)        
    case default

    end select

end subroutine grid_derivative

subroutine grid_derivative_exp
    use mod_constants, only : nintplt_node6e
    use mod_constants, only : nintplt_node8e
    use mod_constants, only : nintplt_scsl4e,nintplt_scsl6e
    use mod_variables, only : nd3int_grd
    implicit none
    external :: ve_via_node6e,ve_via_node8e
    external :: ve_via_scsl4e,ve_via_scsl6e

    select case(nd3int_grd)
    case(nintplt_node6e)
        call grid_derivative_dergrd_exp(ve_via_node6e)
    case(nintplt_node8e)
        call grid_derivative_dergrd_exp(ve_via_node8e)
    case(nintplt_scsl4e)
        call grid_derivative_dergrd_exp(ve_via_scsl4e)
    case(nintplt_scsl6e)
        call grid_derivative_dergrd_exp(ve_via_scsl6e)        
    case default

    end select

end subroutine grid_derivative_exp

subroutine grid_derivative_cc
    use mod_constants, only : nintplt_node2,nintplt_node4
    use mod_constants, only : nintplt_node6,nintplt_node6w
    use mod_constants, only : nintplt_node8
    use mod_constants, only : nintplt_scsl4,nintplt_scsl6
    use mod_variables, only : nd3int_grd
    implicit none
    external :: ve_via_node2cc,ve_via_node4cc
    external :: ve_via_node6cc,ve_via_node6wcc,ve_via_node8cc
    external :: ve_via_scsl4cc,ve_via_scsl6cc

    select case(nd3int_grd)
    case(nintplt_node2)
        call grid_derivative_dergrd_cc(ve_via_node2cc)
    case(nintplt_node4)
        call grid_derivative_dergrd_cc(ve_via_node4cc)
    case(nintplt_node6)
        call grid_derivative_dergrd_cc(ve_via_node6cc)
    case(nintplt_node6w)
        call grid_derivative_dergrd_cc(ve_via_node6wcc)
    case(nintplt_node8)
        call grid_derivative_dergrd_cc(ve_via_node8cc)
    case(nintplt_scsl4)
        call grid_derivative_dergrd_cc(ve_via_scsl4cc)
    case(nintplt_scsl6)
        call grid_derivative_dergrd_cc(ve_via_scsl6cc)        
    case default

    end select

end subroutine grid_derivative_cc

subroutine grid_derivative_dergrd(sub_ve3)
    use mod_constants, only : nderive_edge6,nderive_edge4
    use mod_constants, only : nderive_edge2,nderive_ehen4
    use mod_constants, only : nderive_ehen6,nderive_ehen8
    use mod_constants, only : nderive_ehcs6,nderive_scsl4
    use mod_constants, only : nderive_scsl6
    use mod_variables, only : nd3der_grd,nghnode,nghedge
    implicit none
    external :: sub_ve3
    external :: dn_via_edge6,dn_via_edge4,dn_via_edge2
    external :: dn_via_ehen6,dn_via_ehen4,dn_via_ehen8,dn_via_ehcs6
    external :: dn_via_scsl4,dn_via_scsl6
    
    select case(nd3der_grd)
    case(nderive_edge2)
       call grid_derivative_intplt(sub_ve3,dn_via_edge2)
    case(nderive_edge4)
       call grid_derivative_intplt(sub_ve3,dn_via_edge4)
    case(nderive_edge6)
       call grid_derivative_intplt(sub_ve3,dn_via_edge6)
    case(nderive_ehen4)
       call grid_derivative_intplt(sub_ve3,dn_via_ehen4)
    case(nderive_ehen6)
       call grid_derivative_intplt(sub_ve3,dn_via_ehen6)
    case(nderive_ehen8)
       call grid_derivative_intplt(sub_ve3,dn_via_ehen8)
    case(nderive_ehcs6)
       call grid_derivative_intplt(sub_ve3,dn_via_ehcs6)
    case(nderive_scsl4)
       call grid_derivative_intplt(sub_ve3,dn_via_scsl4)
    case(nderive_scsl6)
       call grid_derivative_intplt(sub_ve3,dn_via_scsl6)       
    case default

    end select

end subroutine grid_derivative_dergrd

subroutine grid_derivative_dergrd_exp(sub_ve3)
    use mod_constants, only : nderive_ehen6e,nderive_ehen8e
    use mod_constants, only : nderive_ehcs6e,nderive_scsl4e
    use mod_constants, only : nderive_scsl6e
    use mod_variables, only : nd3der_grd,nghnode,nghedge
    implicit none
    external :: sub_ve3
    external :: dn_via_ehen6e,dn_via_ehen8e,dn_via_ehcs6e
    external :: dn_via_scsl4e,dn_via_scsl6e
    
    select case(nd3der_grd)
    case(nderive_ehen6e)
       call grid_derivative_intplt_exp(sub_ve3,dn_via_ehen6e)
    case(nderive_ehen8e)
       call grid_derivative_intplt_exp(sub_ve3,dn_via_ehen8e)
    case(nderive_ehcs6e)
       call grid_derivative_intplt_exp(sub_ve3,dn_via_ehcs6e)
    case(nderive_scsl4e)
       call grid_derivative_intplt_exp(sub_ve3,dn_via_scsl4e)
    case(nderive_scsl6e)
       call grid_derivative_intplt_exp(sub_ve3,dn_via_scsl6e)       
    case default

    end select

end subroutine grid_derivative_dergrd_exp

subroutine grid_derivative_dergrd_cc(sub_ve3)
    use mod_constants, only : nderive_edge6,nderive_edge4
    use mod_constants, only : nderive_edge2,nderive_ehen4
    use mod_constants, only : nderive_ehen6,nderive_ehen8
    use mod_constants, only : nderive_ehcs6,nderive_scsl4
    use mod_constants, only : nderive_scsl6
    use mod_variables, only : nd3der_grd,nghnode,nghedge
    implicit none
    external :: sub_ve3
    external :: dn_via_edge6cc,dn_via_edge4cc,dn_via_edge2cc
    external :: dn_via_ehen6cc,dn_via_ehen4cc,dn_via_ehen8cc,dn_via_ehcs6cc
    external :: dn_via_scsl4cc,dn_via_scsl6cc
    
    select case(nd3der_grd)
    case(nderive_edge2)
       call grid_derivative_intplt_cc(sub_ve3,dn_via_edge2cc)
    case(nderive_edge4)
       call grid_derivative_intplt_cc(sub_ve3,dn_via_edge4cc)
    case(nderive_edge6)
       call grid_derivative_intplt_cc(sub_ve3,dn_via_edge6cc)
    case(nderive_ehen4)
       call grid_derivative_intplt_cc(sub_ve3,dn_via_ehen4cc)
    case(nderive_ehen6)
       call grid_derivative_intplt_cc(sub_ve3,dn_via_ehen6cc)
    case(nderive_ehen8)
       call grid_derivative_intplt_cc(sub_ve3,dn_via_ehen8cc)
    case(nderive_ehcs6)
       call grid_derivative_intplt_cc(sub_ve3,dn_via_ehcs6cc)
    case(nderive_scsl4)
       call grid_derivative_intplt_cc(sub_ve3,dn_via_scsl4cc)
    case(nderive_scsl6)
       call grid_derivative_intplt_cc(sub_ve3,dn_via_scsl6cc)       
    case default

    end select

end subroutine grid_derivative_dergrd_cc

subroutine grid_derivative_intplt(sub_ve3,sub_dn3)
    use mod_constants, only : nintplt_node2,nintplt_node4
    use mod_constants, only : nintplt_node6,nintplt_node6w
    use mod_constants, only : nintplt_node8
    use mod_constants, only : nintplt_scsl4,nintplt_scsl6
    use mod_variables, only : nd2int_grd
    implicit none
    external :: sub_ve3,sub_dn3
    external :: ve_via_node2,ve_via_node4
    external :: ve_via_node6,ve_via_node6w,ve_via_node8
    external :: ve_via_scsl4,ve_via_scsl6


    select case(nd2int_grd)
    case(nintplt_node2)
        call grid_derivative_scheme(sub_ve3,sub_dn3,ve_via_node2)
    case(nintplt_node4)
        call grid_derivative_scheme(sub_ve3,sub_dn3,ve_via_node4)
    case(nintplt_node6)
        call grid_derivative_scheme(sub_ve3,sub_dn3,ve_via_node6)
    case(nintplt_node6w)
        call grid_derivative_scheme(sub_ve3,sub_dn3,ve_via_node6w)
    case(nintplt_node8)
        call grid_derivative_scheme(sub_ve3,sub_dn3,ve_via_node8)
    case(nintplt_scsl4)
        call grid_derivative_scheme(sub_ve3,sub_dn3,ve_via_scsl4)
    case(nintplt_scsl6)
        call grid_derivative_scheme(sub_ve3,sub_dn3,ve_via_scsl6)         
    case default

    end select

end subroutine grid_derivative_intplt

subroutine grid_derivative_intplt_exp(sub_ve3,sub_dn3)
    use mod_constants, only : nintplt_node6e
    use mod_constants, only : nintplt_node8e
    use mod_constants, only : nintplt_scsl4e,nintplt_scsl6e
    use mod_variables, only : nd2int_grd
    implicit none
    external :: sub_ve3,sub_dn3
    external :: ve_via_node6e,ve_via_node8e
    external :: ve_via_scsl4e,ve_via_scsl6e


    select case(nd2int_grd)
    case(nintplt_node6e)
        call grid_derivative_scheme_exp(sub_ve3,sub_dn3,ve_via_node6e)
    case(nintplt_node8e)
        call grid_derivative_scheme_exp(sub_ve3,sub_dn3,ve_via_node8e)
    case(nintplt_scsl4e)
        call grid_derivative_scheme_exp(sub_ve3,sub_dn3,ve_via_scsl4e)
    case(nintplt_scsl6e)
        call grid_derivative_scheme_exp(sub_ve3,sub_dn3,ve_via_scsl6e)         
    case default

    end select

end subroutine grid_derivative_intplt_exp

subroutine grid_derivative_intplt_cc(sub_ve3,sub_dn3)
    use mod_constants, only : nintplt_node2,nintplt_node4
    use mod_constants, only : nintplt_node6,nintplt_node6w
    use mod_constants, only : nintplt_node8
    use mod_constants, only : nintplt_scsl4,nintplt_scsl6
    use mod_variables, only : nd2int_grd
    implicit none
    external :: sub_ve3,sub_dn3
    external :: ve_via_node2cc,ve_via_node4cc
    external :: ve_via_node6cc,ve_via_node6wcc,ve_via_node8cc
    external :: ve_via_scsl4cc,ve_via_scsl6cc


    select case(nd2int_grd)
    case(nintplt_node2)
        call grid_derivative_scheme_cc(sub_ve3,sub_dn3,ve_via_node2cc)
    case(nintplt_node4)
        call grid_derivative_scheme_cc(sub_ve3,sub_dn3,ve_via_node4cc)
    case(nintplt_node6)
        call grid_derivative_scheme_cc(sub_ve3,sub_dn3,ve_via_node6cc)
    case(nintplt_node6w)
        call grid_derivative_scheme_cc(sub_ve3,sub_dn3,ve_via_node6wcc)
    case(nintplt_node8)
        call grid_derivative_scheme_cc(sub_ve3,sub_dn3,ve_via_node8cc)
    case(nintplt_scsl4)
        call grid_derivative_scheme_cc(sub_ve3,sub_dn3,ve_via_scsl4cc)
    case(nintplt_scsl6)
        call grid_derivative_scheme_cc(sub_ve3,sub_dn3,ve_via_scsl6cc)         
    case default

    end select

end subroutine grid_derivative_intplt_cc

subroutine grid_derivative_scheme(sub_ve3,sub_dn3,sub_ve2)
    use mod_constants, only : nderive_edge6,nderive_edge4
    use mod_constants, only : nderive_edge2,nderive_ehen4
    use mod_constants, only : nderive_ehen6,nderive_ehen8
    use mod_constants, only : nderive_ehcs6,nderive_scsl4
    use mod_constants, only : nderive_scsl6    
    use mod_variables, only : nd2der_grd,nghnode,nghedge
    implicit none
    external :: sub_ve3,sub_dn3,sub_ve2
    external :: dn_via_edge6,dn_via_edge4,dn_via_edge2
    external :: dn_via_ehen6,dn_via_ehen4,dn_via_ehen8,dn_via_ehcs6
    external :: dn_via_scsl4,dn_via_scsl6

    select case(nd2der_grd)
    case(nderive_edge2)
       call calc_grid_derivative(nghnode,nghedge,sub_ve3,sub_dn3,sub_ve2,dn_via_edge2)
    case(nderive_edge4)
       call calc_grid_derivative(nghnode,nghedge,sub_ve3,sub_dn3,sub_ve2,dn_via_edge4)
    case(nderive_edge6)
       call calc_grid_derivative(nghnode,nghedge,sub_ve3,sub_dn3,sub_ve2,dn_via_edge6)
    case(nderive_ehen4)
       call calc_grid_derivative(nghnode,nghedge,sub_ve3,sub_dn3,sub_ve2,dn_via_ehen4)
    case(nderive_ehen6)
       call calc_grid_derivative(nghnode,nghedge,sub_ve3,sub_dn3,sub_ve2,dn_via_ehen6)
    case(nderive_ehen8)
       call calc_grid_derivative(nghnode,nghedge,sub_ve3,sub_dn3,sub_ve2,dn_via_ehen8)
    case(nderive_ehcs6)
       call calc_grid_derivative(nghnode,nghedge,sub_ve3,sub_dn3,sub_ve2,dn_via_ehcs6)
    case(nderive_scsl4)
       call calc_grid_derivative(nghnode,nghedge,sub_ve3,sub_dn3,sub_ve2,dn_via_scsl4)
    case(nderive_scsl6)
       call calc_grid_derivative(nghnode,nghedge,sub_ve3,sub_dn3,sub_ve2,dn_via_scsl6)       
    case default

    end select

end subroutine grid_derivative_scheme

subroutine grid_derivative_scheme_exp(sub_ve3,sub_dn3,sub_ve2)
    use mod_constants, only : nderive_ehen6e,nderive_ehen8e
    use mod_constants, only : nderive_ehcs6e,nderive_scsl4e
    use mod_constants, only : nderive_scsl6e
    use mod_variables, only : nd2der_grd,nghnode,nghedge
    implicit none
    external :: sub_ve3,sub_dn3,sub_ve2
    external :: dn_via_ehen6e,dn_via_ehen8e,dn_via_ehcs6e
    external :: dn_via_scsl4e,dn_via_scsl6e

    select case(nd2der_grd)
    case(nderive_ehen6e)
       call calc_grid_derivative_exp(nghnode,nghedge,sub_ve3,sub_dn3,sub_ve2,dn_via_ehen6e)    
    case(nderive_ehcs6e)
       call calc_grid_derivative_exp(nghnode,nghedge,sub_ve3,sub_dn3,sub_ve2,dn_via_ehcs6e)
    case(nderive_ehen8e)
       call calc_grid_derivative_exp(nghnode,nghedge,sub_ve3,sub_dn3,sub_ve2,dn_via_ehen8e)
    case(nderive_scsl4e)
       call calc_grid_derivative_exp(nghnode,nghedge,sub_ve3,sub_dn3,sub_ve2,dn_via_scsl4e)
    case(nderive_scsl6e)
       call calc_grid_derivative_exp(nghnode,nghedge,sub_ve3,sub_dn3,sub_ve2,dn_via_scsl6e)       
    case default

    end select

end subroutine grid_derivative_scheme_exp

subroutine grid_derivative_scheme_cc(sub_ve3,sub_dn3,sub_ve2)
    use mod_constants, only : nderive_edge6,nderive_edge4
    use mod_constants, only : nderive_edge2,nderive_ehen4
    use mod_constants, only : nderive_ehen6,nderive_ehen8
    use mod_constants, only : nderive_ehcs6,nderive_scsl4
    use mod_constants, only : nderive_scsl6    
    use mod_variables, only : nd2der_grd,nghnode,nghedge
    implicit none
    external :: sub_ve3,sub_dn3,sub_ve2
    external :: dn_via_edge6cc,dn_via_edge4cc,dn_via_edge2cc
    external :: dn_via_ehen6cc,dn_via_ehen4cc,dn_via_ehen8cc,dn_via_ehcs6cc
    external :: dn_via_scsl4cc,dn_via_scsl6cc

    select case(nd2der_grd)
    case(nderive_edge2)
       call calc_grid_derivative_cc(nghnode,nghedge,sub_ve3,sub_dn3,sub_ve2,dn_via_edge2cc)
    case(nderive_edge4)
       call calc_grid_derivative_cc(nghnode,nghedge,sub_ve3,sub_dn3,sub_ve2,dn_via_edge4cc)
    case(nderive_edge6)
       call calc_grid_derivative_cc(nghnode,nghedge,sub_ve3,sub_dn3,sub_ve2,dn_via_edge6cc)
    case(nderive_ehen4)
       call calc_grid_derivative_cc(nghnode,nghedge,sub_ve3,sub_dn3,sub_ve2,dn_via_ehen4cc)
    case(nderive_ehen6)
       call calc_grid_derivative_cc(nghnode,nghedge,sub_ve3,sub_dn3,sub_ve2,dn_via_ehen6cc)
    case(nderive_ehen8)
       call calc_grid_derivative_cc(nghnode,nghedge,sub_ve3,sub_dn3,sub_ve2,dn_via_ehen8cc)
    case(nderive_ehcs6)
       call calc_grid_derivative_cc(nghnode,nghedge,sub_ve3,sub_dn3,sub_ve2,dn_via_ehcs6cc)
    case(nderive_scsl4)
       call calc_grid_derivative_cc(nghnode,nghedge,sub_ve3,sub_dn3,sub_ve2,dn_via_scsl4cc)
    case(nderive_scsl6)
       call calc_grid_derivative_cc(nghnode,nghedge,sub_ve3,sub_dn3,sub_ve2,dn_via_scsl6cc)       
    case default

    end select

end subroutine grid_derivative_scheme_cc

subroutine calc_grid_derivative(nghn,nghe,sub_ve3,sub_dn3,sub_ve2,sub_dn2)
#include "validation/fortran_probes/probes.h"
    use mod_kndconsts, only : kind_int
    use mod_constants, only : nfsf_grd_d3int,nfsf_grd_d3der
    use mod_constants, only : nfsf_grd_d2int,nfsf_grd_d2der
    use mod_constants, only : nbc_inter_buf_dyn,nsgl_aver_art
    use mod_datatypes, only : var_block_t
    use mod_fieldvars, only : mb_xyz,mb_sxyz,mb_vol,nblocks
    use mod_interface, only : pre_exchange_bc_var,post_exchange_bc_var,exchange_bc_der
    use mod_interface, only : mb_var_create,mb_var_delete
    use mod_interface, only : calc_mb_dn_via_node3,calc_mb_dn_via_node2
    use mod_interface, only : calc_mb_sxyz_vol,exchange_bc_sxyz
    use mod_interface, only : calc_mb_dn_via_vec_node2,calc_mb_vol
    use mod_interface, only : ghost_bc_var_all,gather_output_mb_var
    use probe_utils, only : probe_write_jac
    implicit none
    integer(kind_int), intent(in) :: nghn,nghe
    external                      :: sub_ve3,sub_dn3
    external                      :: sub_ve2,sub_dn2
    type(var_block_t), pointer    :: mb_der3(:),mb_der2(:),mb_der1(:)
    integer(kind_int)             :: nb

    !
    ! +----*----+----*----+----*- ... ... -*----+----*----+----*----+
    ! 0    1    1    2    2       ... ...          ni-1  ni-1  ni   ni
    ! "+" stand for cell-edge, "*" stand for node(cell-center)
    ! nghn: number of ghost nodes based on (1,ni)
    ! nghc: number of ghost edges based on (0,ni)
    ! nfsf : index of bc flag

    call mb_var_create(mb_vol ,1,1,nghn)
    call mb_var_create(mb_sxyz,1,9,nghn)

    call pre_exchange_bc_var(mb_xyz,1,3,nghn,nbc_inter_buf_dyn,nsgl_aver_art)
    call post_exchange_bc_var(mb_xyz,1,3,nghn,nbc_inter_buf_dyn,nsgl_aver_art)

    call mb_var_create(mb_der3,1, 9,nghn)
    call mb_var_create(mb_der2,1,36,nghn)

    call calc_mb_dn_via_node3(mb_xyz,1,nghn,sub_ve3,nghe,sub_dn3,nfsf_grd_d3int,nfsf_grd_d3der,mb_der3,1,1)   ! x_kc
    call calc_mb_dn_via_node3(mb_xyz,1,nghn,sub_ve3,nghe,sub_dn3,nfsf_grd_d3int,nfsf_grd_d3der,mb_der3,2,2)   ! x_et
    call calc_mb_dn_via_node3(mb_xyz,1,nghn,sub_ve3,nghe,sub_dn3,nfsf_grd_d3int,nfsf_grd_d3der,mb_der3,3,3)   ! x_ct
    call calc_mb_dn_via_node3(mb_xyz,2,nghn,sub_ve3,nghe,sub_dn3,nfsf_grd_d3int,nfsf_grd_d3der,mb_der3,4,1)   ! y_kc
    call calc_mb_dn_via_node3(mb_xyz,2,nghn,sub_ve3,nghe,sub_dn3,nfsf_grd_d3int,nfsf_grd_d3der,mb_der3,5,2)   ! y_et
    call calc_mb_dn_via_node3(mb_xyz,2,nghn,sub_ve3,nghe,sub_dn3,nfsf_grd_d3int,nfsf_grd_d3der,mb_der3,6,3)   ! y_ct
    call calc_mb_dn_via_node3(mb_xyz,3,nghn,sub_ve3,nghe,sub_dn3,nfsf_grd_d3int,nfsf_grd_d3der,mb_der3,7,1)   ! z_kc
    call calc_mb_dn_via_node3(mb_xyz,3,nghn,sub_ve3,nghe,sub_dn3,nfsf_grd_d3int,nfsf_grd_d3der,mb_der3,8,2)   ! z_et
    call calc_mb_dn_via_node3(mb_xyz,3,nghn,sub_ve3,nghe,sub_dn3,nfsf_grd_d3int,nfsf_grd_d3der,mb_der3,9,3)   ! z_ct

    call exchange_bc_der(mb_der3,1,9,nghn,nbc_inter_buf_dyn)

    call calc_mb_dn_via_node2(mb_der3,1,mb_xyz,2,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,1,2)  ! (x_kc*y)_et
    call calc_mb_dn_via_node2(mb_der3,2,mb_xyz,2,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,2,3)  ! (x_et*y)_ct
    call calc_mb_dn_via_node2(mb_der3,3,mb_xyz,2,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,3,1)  ! (x_ct*y)_kc
    call calc_mb_dn_via_node2(mb_der3,4,mb_xyz,3,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,4,2)  ! (y_kc*z)_et
    call calc_mb_dn_via_node2(mb_der3,5,mb_xyz,3,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,5,3)  ! (y_et*z)_ct
    call calc_mb_dn_via_node2(mb_der3,6,mb_xyz,3,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,6,1)  ! (y_ct*z)_kc
    call calc_mb_dn_via_node2(mb_der3,7,mb_xyz,1,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,7,2)  ! (z_kc*x)_et
    call calc_mb_dn_via_node2(mb_der3,8,mb_xyz,1,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,8,3)  ! (z_et*x)_ct
    call calc_mb_dn_via_node2(mb_der3,9,mb_xyz,1,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,9,1)  ! (z_ct*x)_kc

    call calc_mb_dn_via_node2(mb_der3,1,mb_xyz,2,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,10,3) ! (x_kc*y)_ct
    call calc_mb_dn_via_node2(mb_der3,2,mb_xyz,2,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,11,1) ! (x_et*y)_kc
    call calc_mb_dn_via_node2(mb_der3,3,mb_xyz,2,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,12,2) ! (x_ct*y)_et
    call calc_mb_dn_via_node2(mb_der3,4,mb_xyz,3,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,13,3) ! (y_kc*z)_ct
    call calc_mb_dn_via_node2(mb_der3,5,mb_xyz,3,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,14,1) ! (y_et*z)_kc
    call calc_mb_dn_via_node2(mb_der3,6,mb_xyz,3,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,15,2) ! (y_ct*z)_et
    call calc_mb_dn_via_node2(mb_der3,7,mb_xyz,1,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,16,3) ! (z_kc*x)_ct
    call calc_mb_dn_via_node2(mb_der3,8,mb_xyz,1,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,17,1) ! (z_et*x)_kc
    call calc_mb_dn_via_node2(mb_der3,9,mb_xyz,1,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,18,2) ! (z_ct*x)_et


    call calc_mb_dn_via_node2(mb_der3,5,mb_xyz,1,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,19,1) ! (x*y_et)_kc
    call calc_mb_dn_via_node2(mb_der3,6,mb_xyz,1,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,20,2) ! (x*y_ct)_et
    call calc_mb_dn_via_node2(mb_der3,4,mb_xyz,1,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,21,3) ! (x*y_kc)_ct
    call calc_mb_dn_via_node2(mb_der3,8,mb_xyz,2,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,22,1) ! (y*z_et)_kc
    call calc_mb_dn_via_node2(mb_der3,9,mb_xyz,2,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,23,2) ! (y*z_ct)_et
    call calc_mb_dn_via_node2(mb_der3,7,mb_xyz,2,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,24,3) ! (y*z_kc)_ct
    call calc_mb_dn_via_node2(mb_der3,2,mb_xyz,3,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,25,1) ! (z*x_et)_kc
    call calc_mb_dn_via_node2(mb_der3,3,mb_xyz,3,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,26,2) ! (z*x_ct)_et
    call calc_mb_dn_via_node2(mb_der3,1,mb_xyz,3,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,27,3) ! (z*x_kc)_ct

    call calc_mb_dn_via_node2(mb_der3,6,mb_xyz,1,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,28,1) ! (x*y_ct)_kc
    call calc_mb_dn_via_node2(mb_der3,4,mb_xyz,1,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,29,2) ! (x*y_kc)_et
    call calc_mb_dn_via_node2(mb_der3,5,mb_xyz,1,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,30,3) ! (x*y_et)_ct
    call calc_mb_dn_via_node2(mb_der3,9,mb_xyz,2,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,31,1) ! (y*z_ct)_kc
    call calc_mb_dn_via_node2(mb_der3,7,mb_xyz,2,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,32,2) ! (y*z_kc)_et
    call calc_mb_dn_via_node2(mb_der3,8,mb_xyz,2,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,33,3) ! (y*z_et)_ct
    call calc_mb_dn_via_node2(mb_der3,3,mb_xyz,3,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,34,1) ! (z*x_ct)_kc
    call calc_mb_dn_via_node2(mb_der3,1,mb_xyz,3,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,35,2) ! (z*x_kc)_et
    call calc_mb_dn_via_node2(mb_der3,2,mb_xyz,3,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,36,3) ! (z*x_et)_ct

    call calc_mb_sxyz_vol(mb_der3,1,9,mb_der2,1,36,mb_sxyz,1,9,mb_vol,1)

!!    call openfile(101,"fort.101","unknown","unformatted","stream")
!!    call gather_output_mb_var(101,mb_sxyz,1,9,1)
!!    call closefile(101,"fort.101")
!!    stop

    call mb_var_delete(mb_der2)
    call mb_var_delete(mb_der3)

    call exchange_bc_sxyz(mb_sxyz,1,9,nghn)

!!------------------------------------------------------------------------------------------------------------------------------------
!    call mb_var_create(mb_der3,1, 9,nghn)
!
!    call calc_mb_dn_via_node3(mb_sxyz,1,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d3int,nfsf_grd_d3der,mb_der3,1,1)   ! kcx_kc
!    call calc_mb_dn_via_node3(mb_sxyz,4,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d3int,nfsf_grd_d3der,mb_der3,2,2)   ! etx_et
!    call calc_mb_dn_via_node3(mb_sxyz,7,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d3int,nfsf_grd_d3der,mb_der3,3,3)   ! ctx_ct
!    call calc_mb_dn_via_node3(mb_sxyz,2,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d3int,nfsf_grd_d3der,mb_der3,4,1)   ! kcy_kc
!    call calc_mb_dn_via_node3(mb_sxyz,5,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d3int,nfsf_grd_d3der,mb_der3,5,2)   ! ety_et
!    call calc_mb_dn_via_node3(mb_sxyz,8,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d3int,nfsf_grd_d3der,mb_der3,6,3)   ! cty_ct
!    call calc_mb_dn_via_node3(mb_sxyz,3,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d3int,nfsf_grd_d3der,mb_der3,7,1)   ! kcz_kc
!    call calc_mb_dn_via_node3(mb_sxyz,6,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d3int,nfsf_grd_d3der,mb_der3,8,2)   ! etz_et
!    call calc_mb_dn_via_node3(mb_sxyz,9,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d3int,nfsf_grd_d3der,mb_der3,9,3)   ! ctz_ct
!
!    call openfile(101,"fort.101","unknown","unformatted","stream")
!    call gather_output_mb_var(101,mb_der3,1,9,1)
!    call closefile(101,"fort.101")
!
!    call mb_var_delete(mb_der3)
!!------------------------------------------------------------------------------------------------------------------------------------

    call mb_var_create(mb_der1,1,3,nghn)
    call calc_mb_dn_via_vec_node2(mb_sxyz,1,3,mb_xyz,1,3,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der1,1,1)
    call calc_mb_dn_via_vec_node2(mb_sxyz,4,6,mb_xyz,1,3,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der1,2,2)
    call calc_mb_dn_via_vec_node2(mb_sxyz,7,9,mb_xyz,1,3,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der1,3,3)
    call calc_mb_vol(mb_der1,1,3,mb_vol,1)
    call mb_var_delete(mb_der1)

    call pre_exchange_bc_var(mb_vol,1,1,nghn,nbc_inter_buf_dyn,nsgl_aver_art)
    call post_exchange_bc_var(mb_vol,1,1,nghn,nbc_inter_buf_dyn,nsgl_aver_art)


    call ghost_bc_var_all(mb_sxyz,1,9)
    call ghost_bc_var_all(mb_vol,1,1)

    ! Probe: dump jac for JAX validation (enabled by -DPROBE_METRIC)
    do nb = 1, nblocks
        PROBE_DUMP_JAC(nb, mb_vol(nb)%fld(1)%r3d)
    end do

end subroutine calc_grid_derivative

subroutine calc_grid_derivative_exp(nghn,nghe,sub_ve3,sub_dn3,sub_ve2,sub_dn2)
    use mod_kndconsts, only : kind_int
    use mod_constants, only : nfsf_grd_d3int,nfsf_grd_d3der
    use mod_constants, only : nfsf_grd_d2int,nfsf_grd_d2der
    use mod_constants, only : nbc_inter_buf_dyn,nsgl_aver_art
    use mod_datatypes, only : var_block_t
    use mod_fieldvars, only : mb_xyz,mb_sxyz,mb_vol
    use mod_interface, only : pre_exchange_bc_var,post_exchange_bc_var
    use mod_interface, only : mb_var_create,mb_var_delete
    use mod_interface, only : extending_ghost_points
    use mod_interface, only : calc_mb_dn_via_node3_exp,calc_mb_dn_via_node2_exp !,calc_mb_dn_via_node3
    use mod_interface, only : calc_mb_sxyz_vol_exp
    use mod_interface, only : calc_mb_dn_via_vec_node2,calc_mb_vol
    use mod_interface, only : ghost_bc_var_all !,gather_output_mb_var
    implicit none
    integer(kind_int), intent(in) :: nghn,nghe
    external                      :: sub_ve3,sub_dn3
    external                      :: sub_ve2,sub_dn2
    type(var_block_t), pointer    :: mb_der3(:),mb_der2(:),mb_der1(:)

    !
    ! +----*----+----*----+----*- ... ... -*----+----*----+----*----+
    ! 0    1    1    2    2       ... ...          ni-1  ni-1  ni   ni
    ! "+" stands for cell-edge, "*" stands for node(cell-center)
    ! nghn: number of ghost nodes based on (1,ni)
    ! nghc: number of ghost edges based on (0,ni)
    ! nfsf : index of bc flag

    call mb_var_create(mb_vol ,1,1,nghn)
    call mb_var_create(mb_sxyz,1,9,nghn)

    call extending_ghost_points(mb_xyz,1,3,3*nghn)

    call mb_var_create(mb_der3,1, 9,nghn*2)
    call mb_var_create(mb_der2,1,36,nghn)

    call calc_mb_dn_via_node3_exp(mb_xyz,1,nghn*3,sub_ve3,nghn*2+nghe,sub_dn3,nfsf_grd_d3int,nfsf_grd_d3der,mb_der3,1,1)   ! x_kc
    call calc_mb_dn_via_node3_exp(mb_xyz,1,nghn*3,sub_ve3,nghn*2+nghe,sub_dn3,nfsf_grd_d3int,nfsf_grd_d3der,mb_der3,2,2)   ! x_et
    call calc_mb_dn_via_node3_exp(mb_xyz,1,nghn*3,sub_ve3,nghn*2+nghe,sub_dn3,nfsf_grd_d3int,nfsf_grd_d3der,mb_der3,3,3)   ! x_ct
    call calc_mb_dn_via_node3_exp(mb_xyz,2,nghn*3,sub_ve3,nghn*2+nghe,sub_dn3,nfsf_grd_d3int,nfsf_grd_d3der,mb_der3,4,1)   ! y_kc
    call calc_mb_dn_via_node3_exp(mb_xyz,2,nghn*3,sub_ve3,nghn*2+nghe,sub_dn3,nfsf_grd_d3int,nfsf_grd_d3der,mb_der3,5,2)   ! y_et
    call calc_mb_dn_via_node3_exp(mb_xyz,2,nghn*3,sub_ve3,nghn*2+nghe,sub_dn3,nfsf_grd_d3int,nfsf_grd_d3der,mb_der3,6,3)   ! y_ct
    call calc_mb_dn_via_node3_exp(mb_xyz,3,nghn*3,sub_ve3,nghn*2+nghe,sub_dn3,nfsf_grd_d3int,nfsf_grd_d3der,mb_der3,7,1)   ! z_kc
    call calc_mb_dn_via_node3_exp(mb_xyz,3,nghn*3,sub_ve3,nghn*2+nghe,sub_dn3,nfsf_grd_d3int,nfsf_grd_d3der,mb_der3,8,2)   ! z_et
    call calc_mb_dn_via_node3_exp(mb_xyz,3,nghn*3,sub_ve3,nghn*2+nghe,sub_dn3,nfsf_grd_d3int,nfsf_grd_d3der,mb_der3,9,3)   ! z_ct

    call calc_mb_dn_via_node2_exp(mb_der3,1,mb_xyz,2,nghn*2,sub_ve2,nghn+nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,1,2)  ! (x_kc*y)_et
    call calc_mb_dn_via_node2_exp(mb_der3,2,mb_xyz,2,nghn*2,sub_ve2,nghn+nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,2,3)  ! (x_et*y)_ct
    call calc_mb_dn_via_node2_exp(mb_der3,3,mb_xyz,2,nghn*2,sub_ve2,nghn+nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,3,1)  ! (x_ct*y)_kc
    call calc_mb_dn_via_node2_exp(mb_der3,4,mb_xyz,3,nghn*2,sub_ve2,nghn+nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,4,2)  ! (y_kc*z)_et
    call calc_mb_dn_via_node2_exp(mb_der3,5,mb_xyz,3,nghn*2,sub_ve2,nghn+nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,5,3)  ! (y_et*z)_ct
    call calc_mb_dn_via_node2_exp(mb_der3,6,mb_xyz,3,nghn*2,sub_ve2,nghn+nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,6,1)  ! (y_ct*z)_kc
    call calc_mb_dn_via_node2_exp(mb_der3,7,mb_xyz,1,nghn*2,sub_ve2,nghn+nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,7,2)  ! (z_kc*x)_et
    call calc_mb_dn_via_node2_exp(mb_der3,8,mb_xyz,1,nghn*2,sub_ve2,nghn+nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,8,3)  ! (z_et*x)_ct
    call calc_mb_dn_via_node2_exp(mb_der3,9,mb_xyz,1,nghn*2,sub_ve2,nghn+nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,9,1)  ! (z_ct*x)_kc

    call calc_mb_dn_via_node2_exp(mb_der3,1,mb_xyz,2,nghn*2,sub_ve2,nghn+nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,10,3) ! (x_kc*y)_ct
    call calc_mb_dn_via_node2_exp(mb_der3,2,mb_xyz,2,nghn*2,sub_ve2,nghn+nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,11,1) ! (x_et*y)_kc
    call calc_mb_dn_via_node2_exp(mb_der3,3,mb_xyz,2,nghn*2,sub_ve2,nghn+nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,12,2) ! (x_ct*y)_et
    call calc_mb_dn_via_node2_exp(mb_der3,4,mb_xyz,3,nghn*2,sub_ve2,nghn+nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,13,3) ! (y_kc*z)_ct
    call calc_mb_dn_via_node2_exp(mb_der3,5,mb_xyz,3,nghn*2,sub_ve2,nghn+nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,14,1) ! (y_et*z)_kc
    call calc_mb_dn_via_node2_exp(mb_der3,6,mb_xyz,3,nghn*2,sub_ve2,nghn+nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,15,2) ! (y_ct*z)_et
    call calc_mb_dn_via_node2_exp(mb_der3,7,mb_xyz,1,nghn*2,sub_ve2,nghn+nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,16,3) ! (z_kc*x)_ct
    call calc_mb_dn_via_node2_exp(mb_der3,8,mb_xyz,1,nghn*2,sub_ve2,nghn+nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,17,1) ! (z_et*x)_kc
    call calc_mb_dn_via_node2_exp(mb_der3,9,mb_xyz,1,nghn*2,sub_ve2,nghn+nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,18,2) ! (z_ct*x)_et


    call calc_mb_dn_via_node2_exp(mb_der3,5,mb_xyz,1,nghn*2,sub_ve2,nghn+nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,19,1) ! (x*y_et)_kc
    call calc_mb_dn_via_node2_exp(mb_der3,6,mb_xyz,1,nghn*2,sub_ve2,nghn+nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,20,2) ! (x*y_ct)_et
    call calc_mb_dn_via_node2_exp(mb_der3,4,mb_xyz,1,nghn*2,sub_ve2,nghn+nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,21,3) ! (x*y_kc)_ct
    call calc_mb_dn_via_node2_exp(mb_der3,8,mb_xyz,2,nghn*2,sub_ve2,nghn+nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,22,1) ! (y*z_et)_kc
    call calc_mb_dn_via_node2_exp(mb_der3,9,mb_xyz,2,nghn*2,sub_ve2,nghn+nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,23,2) ! (y*z_ct)_et
    call calc_mb_dn_via_node2_exp(mb_der3,7,mb_xyz,2,nghn*2,sub_ve2,nghn+nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,24,3) ! (y*z_kc)_ct
    call calc_mb_dn_via_node2_exp(mb_der3,2,mb_xyz,3,nghn*2,sub_ve2,nghn+nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,25,1) ! (z*x_et)_kc
    call calc_mb_dn_via_node2_exp(mb_der3,3,mb_xyz,3,nghn*2,sub_ve2,nghn+nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,26,2) ! (z*x_ct)_et
    call calc_mb_dn_via_node2_exp(mb_der3,1,mb_xyz,3,nghn*2,sub_ve2,nghn+nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,27,3) ! (z*x_kc)_ct

    call calc_mb_dn_via_node2_exp(mb_der3,6,mb_xyz,1,nghn*2,sub_ve2,nghn+nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,28,1) ! (x*y_ct)_kc
    call calc_mb_dn_via_node2_exp(mb_der3,4,mb_xyz,1,nghn*2,sub_ve2,nghn+nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,29,2) ! (x*y_kc)_et
    call calc_mb_dn_via_node2_exp(mb_der3,5,mb_xyz,1,nghn*2,sub_ve2,nghn+nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,30,3) ! (x*y_et)_ct
    call calc_mb_dn_via_node2_exp(mb_der3,9,mb_xyz,2,nghn*2,sub_ve2,nghn+nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,31,1) ! (y*z_ct)_kc
    call calc_mb_dn_via_node2_exp(mb_der3,7,mb_xyz,2,nghn*2,sub_ve2,nghn+nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,32,2) ! (y*z_kc)_et
    call calc_mb_dn_via_node2_exp(mb_der3,8,mb_xyz,2,nghn*2,sub_ve2,nghn+nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,33,3) ! (y*z_et)_ct
    call calc_mb_dn_via_node2_exp(mb_der3,3,mb_xyz,3,nghn*2,sub_ve2,nghn+nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,34,1) ! (z*x_ct)_kc
    call calc_mb_dn_via_node2_exp(mb_der3,1,mb_xyz,3,nghn*2,sub_ve2,nghn+nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,35,2) ! (z*x_kc)_et
    call calc_mb_dn_via_node2_exp(mb_der3,2,mb_xyz,3,nghn*2,sub_ve2,nghn+nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,36,3) ! (z*x_et)_ct

    call calc_mb_sxyz_vol_exp(mb_der3,1,9,mb_der2,1,36,mb_sxyz,1,9,mb_vol,1)

    call mb_var_delete(mb_der2)
    call mb_var_delete(mb_der3)

!!------------------------------------------------------------------------------------------------------------------------------------
!    call mb_var_create(mb_der3,1, 9,nghn)
!
!    call calc_mb_dn_via_node3(mb_sxyz,1,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d3int,nfsf_grd_d3der,mb_der3,1,1)   ! kcx_kc
!    call calc_mb_dn_via_node3(mb_sxyz,4,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d3int,nfsf_grd_d3der,mb_der3,2,2)   ! etx_et
!    call calc_mb_dn_via_node3(mb_sxyz,7,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d3int,nfsf_grd_d3der,mb_der3,3,3)   ! ctx_ct
!    call calc_mb_dn_via_node3(mb_sxyz,2,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d3int,nfsf_grd_d3der,mb_der3,4,1)   ! kcy_kc
!    call calc_mb_dn_via_node3(mb_sxyz,5,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d3int,nfsf_grd_d3der,mb_der3,5,2)   ! ety_et
!    call calc_mb_dn_via_node3(mb_sxyz,8,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d3int,nfsf_grd_d3der,mb_der3,6,3)   ! cty_ct
!    call calc_mb_dn_via_node3(mb_sxyz,3,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d3int,nfsf_grd_d3der,mb_der3,7,1)   ! kcz_kc
!    call calc_mb_dn_via_node3(mb_sxyz,6,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d3int,nfsf_grd_d3der,mb_der3,8,2)   ! etz_et
!    call calc_mb_dn_via_node3(mb_sxyz,9,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d3int,nfsf_grd_d3der,mb_der3,9,3)   ! ctz_ct
!
!    call openfile(101,"fort.101","unknown","unformatted","stream")
!    call gather_output_mb_var(101,mb_der3,1,9,1)
!    call closefile(101,"fort.101")
!
!    call mb_var_delete(mb_der3)
!!------------------------------------------------------------------------------------------------------------------------------------

    call mb_var_create(mb_der1,1,3,nghn)
    call calc_mb_dn_via_vec_node2(mb_sxyz,1,3,mb_xyz,1,3,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der1,1,1)
    call calc_mb_dn_via_vec_node2(mb_sxyz,4,6,mb_xyz,1,3,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der1,2,2)
    call calc_mb_dn_via_vec_node2(mb_sxyz,7,9,mb_xyz,1,3,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der1,3,3)
    call calc_mb_vol(mb_der1,1,3,mb_vol,1)
    call mb_var_delete(mb_der1)

    call pre_exchange_bc_var(mb_vol,1,1,nghn,nbc_inter_buf_dyn,nsgl_aver_art)
    call post_exchange_bc_var(mb_vol,1,1,nghn,nbc_inter_buf_dyn,nsgl_aver_art)


    call ghost_bc_var_all(mb_sxyz,1,9)
    call ghost_bc_var_all(mb_vol,1,1)

end subroutine calc_grid_derivative_exp

subroutine calc_grid_derivative_cc(nghn,nghe,sub_ve3,sub_dn3,sub_ve2,sub_dn2)
    use mod_kndconsts, only : kind_int
    use mod_constants, only : nfsf_grd_d3int,nfsf_grd_d3der
    use mod_constants, only : nfsf_grd_d2int,nfsf_grd_d2der
    use mod_constants, only : nbc_inter_buf_dyn,nsgl_aver_art
    use mod_datatypes, only : var_block_t
    use mod_fieldvars, only : mb_xyzcc,mb_sxyzcc,mb_volcc
    use mod_fieldvars, only : mb_xyzsp,mb_sxyzsp,mb_volsp
    use mod_interface, only : pre_exchange_cc_bc_var,post_exchange_cc_bc_var,exchange_cc_bc_der
    use mod_interface, only : pre_exchange_bc_var_sp,post_exchange_bc_var_sp,exchange_sp_bc_der
    use mod_interface, only : mb_var_cc_create,mb_var_create_sp,mb_var_delete
    use mod_interface, only : calc_cc_mb_dn_via_node3,calc_cc_mb_dn_via_node2
    use mod_interface, only : calc_sp_mb_dn_via_node3,calc_sp_mb_dn_via_node2
    use mod_interface, only : calc_cc_mb_sxyz_vol,exchange_cc_bc_sxyz
    use mod_interface, only : calc_sp_mb_sxyz_vol,exchange_sp_bc_sxyz
    use mod_interface, only : calc_cc_mb_dn_via_vec_node2,calc_cc_mb_vol
    use mod_interface, only : calc_sp_mb_dn_via_vec_node2,calc_sp_mb_vol
    use mod_interface, only : ghost_bc_var_all_cc,ghost_bc_var_all_sp !,gather_output_mb_var
    implicit none
    integer(kind_int), intent(in) :: nghn,nghe
    external                      :: sub_ve3,sub_dn3
    external                      :: sub_ve2,sub_dn2
    type(var_block_t), pointer    :: mb_der3(:),mb_der2(:),mb_der1(:)

    !
    ! +----*----+----*----+----*- ... ... -*----+----*------+-----*-----+
    ! 0    1    1    2    2       ... ...         2*ni-1 2*ni-1 2*ni  2*ni
    ! "+" stand for cell-edge, "*" stand for node(cell-center)
    ! nghn: number of ghost nodes based on (1,2*ni)
    ! nghc: number of ghost edges based on (0,2*ni)
    ! nfsf : index of bc flag

    call mb_var_cc_create(mb_volcc ,1,1,nghn)
    call mb_var_cc_create(mb_sxyzcc,1,9,nghn)

    call pre_exchange_cc_bc_var(mb_xyzcc,1,3,nghn,nbc_inter_buf_dyn,nsgl_aver_art)
    call post_exchange_cc_bc_var(mb_xyzcc,1,3,nghn,nbc_inter_buf_dyn,nsgl_aver_art)

    call mb_var_cc_create(mb_der3,1, 9,nghn)
    call mb_var_cc_create(mb_der2,1,36,nghn)

    call calc_cc_mb_dn_via_node3(mb_xyzcc,1,nghn,sub_ve3,nghe,sub_dn3,nfsf_grd_d3int,nfsf_grd_d3der,mb_der3,1,1)   ! x_kc
    call calc_cc_mb_dn_via_node3(mb_xyzcc,1,nghn,sub_ve3,nghe,sub_dn3,nfsf_grd_d3int,nfsf_grd_d3der,mb_der3,2,2)   ! x_et
    call calc_cc_mb_dn_via_node3(mb_xyzcc,1,nghn,sub_ve3,nghe,sub_dn3,nfsf_grd_d3int,nfsf_grd_d3der,mb_der3,3,3)   ! x_ct
    call calc_cc_mb_dn_via_node3(mb_xyzcc,2,nghn,sub_ve3,nghe,sub_dn3,nfsf_grd_d3int,nfsf_grd_d3der,mb_der3,4,1)   ! y_kc
    call calc_cc_mb_dn_via_node3(mb_xyzcc,2,nghn,sub_ve3,nghe,sub_dn3,nfsf_grd_d3int,nfsf_grd_d3der,mb_der3,5,2)   ! y_et
    call calc_cc_mb_dn_via_node3(mb_xyzcc,2,nghn,sub_ve3,nghe,sub_dn3,nfsf_grd_d3int,nfsf_grd_d3der,mb_der3,6,3)   ! y_ct
    call calc_cc_mb_dn_via_node3(mb_xyzcc,3,nghn,sub_ve3,nghe,sub_dn3,nfsf_grd_d3int,nfsf_grd_d3der,mb_der3,7,1)   ! z_kc
    call calc_cc_mb_dn_via_node3(mb_xyzcc,3,nghn,sub_ve3,nghe,sub_dn3,nfsf_grd_d3int,nfsf_grd_d3der,mb_der3,8,2)   ! z_et
    call calc_cc_mb_dn_via_node3(mb_xyzcc,3,nghn,sub_ve3,nghe,sub_dn3,nfsf_grd_d3int,nfsf_grd_d3der,mb_der3,9,3)   ! z_ct

    call exchange_cc_bc_der(mb_der3,1,9,nghn,nbc_inter_buf_dyn)

    call calc_cc_mb_dn_via_node2(mb_der3,1,mb_xyzcc,2,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,1,2)  ! (x_kc*y)_et
    call calc_cc_mb_dn_via_node2(mb_der3,2,mb_xyzcc,2,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,2,3)  ! (x_et*y)_ct
    call calc_cc_mb_dn_via_node2(mb_der3,3,mb_xyzcc,2,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,3,1)  ! (x_ct*y)_kc
    call calc_cc_mb_dn_via_node2(mb_der3,4,mb_xyzcc,3,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,4,2)  ! (y_kc*z)_et
    call calc_cc_mb_dn_via_node2(mb_der3,5,mb_xyzcc,3,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,5,3)  ! (y_et*z)_ct
    call calc_cc_mb_dn_via_node2(mb_der3,6,mb_xyzcc,3,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,6,1)  ! (y_ct*z)_kc
    call calc_cc_mb_dn_via_node2(mb_der3,7,mb_xyzcc,1,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,7,2)  ! (z_kc*x)_et
    call calc_cc_mb_dn_via_node2(mb_der3,8,mb_xyzcc,1,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,8,3)  ! (z_et*x)_ct
    call calc_cc_mb_dn_via_node2(mb_der3,9,mb_xyzcc,1,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,9,1)  ! (z_ct*x)_kc

    call calc_cc_mb_dn_via_node2(mb_der3,1,mb_xyzcc,2,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,10,3) ! (x_kc*y)_ct
    call calc_cc_mb_dn_via_node2(mb_der3,2,mb_xyzcc,2,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,11,1) ! (x_et*y)_kc
    call calc_cc_mb_dn_via_node2(mb_der3,3,mb_xyzcc,2,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,12,2) ! (x_ct*y)_et
    call calc_cc_mb_dn_via_node2(mb_der3,4,mb_xyzcc,3,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,13,3) ! (y_kc*z)_ct
    call calc_cc_mb_dn_via_node2(mb_der3,5,mb_xyzcc,3,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,14,1) ! (y_et*z)_kc
    call calc_cc_mb_dn_via_node2(mb_der3,6,mb_xyzcc,3,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,15,2) ! (y_ct*z)_et
    call calc_cc_mb_dn_via_node2(mb_der3,7,mb_xyzcc,1,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,16,3) ! (z_kc*x)_ct
    call calc_cc_mb_dn_via_node2(mb_der3,8,mb_xyzcc,1,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,17,1) ! (z_et*x)_kc
    call calc_cc_mb_dn_via_node2(mb_der3,9,mb_xyzcc,1,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,18,2) ! (z_ct*x)_et


    call calc_cc_mb_dn_via_node2(mb_der3,5,mb_xyzcc,1,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,19,1) ! (x*y_et)_kc
    call calc_cc_mb_dn_via_node2(mb_der3,6,mb_xyzcc,1,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,20,2) ! (x*y_ct)_et
    call calc_cc_mb_dn_via_node2(mb_der3,4,mb_xyzcc,1,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,21,3) ! (x*y_kc)_ct
    call calc_cc_mb_dn_via_node2(mb_der3,8,mb_xyzcc,2,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,22,1) ! (y*z_et)_kc
    call calc_cc_mb_dn_via_node2(mb_der3,9,mb_xyzcc,2,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,23,2) ! (y*z_ct)_et
    call calc_cc_mb_dn_via_node2(mb_der3,7,mb_xyzcc,2,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,24,3) ! (y*z_kc)_ct
    call calc_cc_mb_dn_via_node2(mb_der3,2,mb_xyzcc,3,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,25,1) ! (z*x_et)_kc
    call calc_cc_mb_dn_via_node2(mb_der3,3,mb_xyzcc,3,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,26,2) ! (z*x_ct)_et
    call calc_cc_mb_dn_via_node2(mb_der3,1,mb_xyzcc,3,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,27,3) ! (z*x_kc)_ct

    call calc_cc_mb_dn_via_node2(mb_der3,6,mb_xyzcc,1,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,28,1) ! (x*y_ct)_kc
    call calc_cc_mb_dn_via_node2(mb_der3,4,mb_xyzcc,1,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,29,2) ! (x*y_kc)_et
    call calc_cc_mb_dn_via_node2(mb_der3,5,mb_xyzcc,1,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,30,3) ! (x*y_et)_ct
    call calc_cc_mb_dn_via_node2(mb_der3,9,mb_xyzcc,2,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,31,1) ! (y*z_ct)_kc
    call calc_cc_mb_dn_via_node2(mb_der3,7,mb_xyzcc,2,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,32,2) ! (y*z_kc)_et
    call calc_cc_mb_dn_via_node2(mb_der3,8,mb_xyzcc,2,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,33,3) ! (y*z_et)_ct
    call calc_cc_mb_dn_via_node2(mb_der3,3,mb_xyzcc,3,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,34,1) ! (z*x_ct)_kc
    call calc_cc_mb_dn_via_node2(mb_der3,1,mb_xyzcc,3,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,35,2) ! (z*x_kc)_et
    call calc_cc_mb_dn_via_node2(mb_der3,2,mb_xyzcc,3,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,36,3) ! (z*x_et)_ct

    call calc_cc_mb_sxyz_vol(mb_der3,1,9,mb_der2,1,36,mb_sxyzcc,1,9,mb_volcc,1)

!!    call openfile(101,"fort.101","unknown","unformatted","stream")
!!    call gather_output_mb_var(101,mb_sxyzcc,1,9,1)
!!    call closefile(101,"fort.101")
!!    stop

    call mb_var_delete(mb_der2)
    call mb_var_delete(mb_der3)

    call exchange_cc_bc_sxyz(mb_sxyzcc,1,9,nghn)

!!------------------------------------------------------------------------------------------------------------------------------------
!    call mb_var_create(mb_der3,1, 9,nghn)
!
!    call calc_mb_dn_via_node3(mb_sxyzcc,1,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d3int,nfsf_grd_d3der,mb_der3,1,1)   ! kcx_kc
!    call calc_mb_dn_via_node3(mb_sxyzcc,4,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d3int,nfsf_grd_d3der,mb_der3,2,2)   ! etx_et
!    call calc_mb_dn_via_node3(mb_sxyzcc,7,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d3int,nfsf_grd_d3der,mb_der3,3,3)   ! ctx_ct
!    call calc_mb_dn_via_node3(mb_sxyzcc,2,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d3int,nfsf_grd_d3der,mb_der3,4,1)   ! kcy_kc
!    call calc_mb_dn_via_node3(mb_sxyzcc,5,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d3int,nfsf_grd_d3der,mb_der3,5,2)   ! ety_et
!    call calc_mb_dn_via_node3(mb_sxyzcc,8,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d3int,nfsf_grd_d3der,mb_der3,6,3)   ! cty_ct
!    call calc_mb_dn_via_node3(mb_sxyzcc,3,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d3int,nfsf_grd_d3der,mb_der3,7,1)   ! kcz_kc
!    call calc_mb_dn_via_node3(mb_sxyzcc,6,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d3int,nfsf_grd_d3der,mb_der3,8,2)   ! etz_et
!    call calc_mb_dn_via_node3(mb_sxyzcc,9,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d3int,nfsf_grd_d3der,mb_der3,9,3)   ! ctz_ct
!
!    call openfile(101,"fort.101","unknown","unformatted","stream")
!    call gather_output_mb_var(101,mb_der3,1,9,1)
!    call closefile(101,"fort.101")
!
!    call mb_var_delete(mb_der3)
!!------------------------------------------------------------------------------------------------------------------------------------

    call mb_var_cc_create(mb_der1,1,3,nghn)
    call calc_cc_mb_dn_via_vec_node2(mb_sxyzcc,1,3,mb_xyzcc,1,3,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der1,1,1)
    call calc_cc_mb_dn_via_vec_node2(mb_sxyzcc,4,6,mb_xyzcc,1,3,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der1,2,2)
    call calc_cc_mb_dn_via_vec_node2(mb_sxyzcc,7,9,mb_xyzcc,1,3,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der1,3,3)
    call calc_cc_mb_vol(mb_der1,1,3,mb_volcc,1)
    call mb_var_delete(mb_der1)

    call pre_exchange_cc_bc_var(mb_volcc,1,1,nghn,nbc_inter_buf_dyn,nsgl_aver_art)
    call post_exchange_cc_bc_var(mb_volcc,1,1,nghn,nbc_inter_buf_dyn,nsgl_aver_art)


    call ghost_bc_var_all_cc(mb_sxyzcc,1,9)
    call ghost_bc_var_all_cc(mb_volcc,1,1)

!-----------------------------------------------------------------SP    
    call mb_var_create_sp(mb_volsp ,1,1,nghn)
    call mb_var_create_sp(mb_sxyzsp,1,9,nghn)

    call pre_exchange_bc_var_sp(mb_xyzsp,1,3,nghn,nbc_inter_buf_dyn,nsgl_aver_art)
    call post_exchange_bc_var_sp(mb_xyzsp,1,3,nghn,nbc_inter_buf_dyn,nsgl_aver_art)

    call mb_var_create_sp(mb_der3,1, 9,nghn)
    call mb_var_create_sp(mb_der2,1,36,nghn)

    call calc_sp_mb_dn_via_node3(mb_xyzsp,1,nghn,sub_ve3,nghe,sub_dn3,nfsf_grd_d3int,nfsf_grd_d3der,mb_der3,1,1)   ! x_kc
    call calc_sp_mb_dn_via_node3(mb_xyzsp,1,nghn,sub_ve3,nghe,sub_dn3,nfsf_grd_d3int,nfsf_grd_d3der,mb_der3,2,2)   ! x_et
    call calc_sp_mb_dn_via_node3(mb_xyzsp,1,nghn,sub_ve3,nghe,sub_dn3,nfsf_grd_d3int,nfsf_grd_d3der,mb_der3,3,3)   ! x_ct
    call calc_sp_mb_dn_via_node3(mb_xyzsp,2,nghn,sub_ve3,nghe,sub_dn3,nfsf_grd_d3int,nfsf_grd_d3der,mb_der3,4,1)   ! y_kc
    call calc_sp_mb_dn_via_node3(mb_xyzsp,2,nghn,sub_ve3,nghe,sub_dn3,nfsf_grd_d3int,nfsf_grd_d3der,mb_der3,5,2)   ! y_et
    call calc_sp_mb_dn_via_node3(mb_xyzsp,2,nghn,sub_ve3,nghe,sub_dn3,nfsf_grd_d3int,nfsf_grd_d3der,mb_der3,6,3)   ! y_ct
    call calc_sp_mb_dn_via_node3(mb_xyzsp,3,nghn,sub_ve3,nghe,sub_dn3,nfsf_grd_d3int,nfsf_grd_d3der,mb_der3,7,1)   ! z_kc
    call calc_sp_mb_dn_via_node3(mb_xyzsp,3,nghn,sub_ve3,nghe,sub_dn3,nfsf_grd_d3int,nfsf_grd_d3der,mb_der3,8,2)   ! z_et
    call calc_sp_mb_dn_via_node3(mb_xyzsp,3,nghn,sub_ve3,nghe,sub_dn3,nfsf_grd_d3int,nfsf_grd_d3der,mb_der3,9,3)   ! z_ct

    call exchange_sp_bc_der(mb_der3,1,9,nghn,nbc_inter_buf_dyn)

    call calc_sp_mb_dn_via_node2(mb_der3,1,mb_xyzsp,2,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,1,2)  ! (x_kc*y)_et
    call calc_sp_mb_dn_via_node2(mb_der3,2,mb_xyzsp,2,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,2,3)  ! (x_et*y)_ct
    call calc_sp_mb_dn_via_node2(mb_der3,3,mb_xyzsp,2,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,3,1)  ! (x_ct*y)_kc
    call calc_sp_mb_dn_via_node2(mb_der3,4,mb_xyzsp,3,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,4,2)  ! (y_kc*z)_et
    call calc_sp_mb_dn_via_node2(mb_der3,5,mb_xyzsp,3,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,5,3)  ! (y_et*z)_ct
    call calc_sp_mb_dn_via_node2(mb_der3,6,mb_xyzsp,3,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,6,1)  ! (y_ct*z)_kc
    call calc_sp_mb_dn_via_node2(mb_der3,7,mb_xyzsp,1,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,7,2)  ! (z_kc*x)_et
    call calc_sp_mb_dn_via_node2(mb_der3,8,mb_xyzsp,1,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,8,3)  ! (z_et*x)_ct
    call calc_sp_mb_dn_via_node2(mb_der3,9,mb_xyzsp,1,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,9,1)  ! (z_ct*x)_kc

    call calc_sp_mb_dn_via_node2(mb_der3,1,mb_xyzsp,2,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,10,3) ! (x_kc*y)_ct
    call calc_sp_mb_dn_via_node2(mb_der3,2,mb_xyzsp,2,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,11,1) ! (x_et*y)_kc
    call calc_sp_mb_dn_via_node2(mb_der3,3,mb_xyzsp,2,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,12,2) ! (x_ct*y)_et
    call calc_sp_mb_dn_via_node2(mb_der3,4,mb_xyzsp,3,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,13,3) ! (y_kc*z)_ct
    call calc_sp_mb_dn_via_node2(mb_der3,5,mb_xyzsp,3,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,14,1) ! (y_et*z)_kc
    call calc_sp_mb_dn_via_node2(mb_der3,6,mb_xyzsp,3,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,15,2) ! (y_ct*z)_et
    call calc_sp_mb_dn_via_node2(mb_der3,7,mb_xyzsp,1,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,16,3) ! (z_kc*x)_ct
    call calc_sp_mb_dn_via_node2(mb_der3,8,mb_xyzsp,1,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,17,1) ! (z_et*x)_kc
    call calc_sp_mb_dn_via_node2(mb_der3,9,mb_xyzsp,1,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,18,2) ! (z_ct*x)_et


    call calc_sp_mb_dn_via_node2(mb_der3,5,mb_xyzsp,1,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,19,1) ! (x*y_et)_kc
    call calc_sp_mb_dn_via_node2(mb_der3,6,mb_xyzsp,1,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,20,2) ! (x*y_ct)_et
    call calc_sp_mb_dn_via_node2(mb_der3,4,mb_xyzsp,1,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,21,3) ! (x*y_kc)_ct
    call calc_sp_mb_dn_via_node2(mb_der3,8,mb_xyzsp,2,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,22,1) ! (y*z_et)_kc
    call calc_sp_mb_dn_via_node2(mb_der3,9,mb_xyzsp,2,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,23,2) ! (y*z_ct)_et
    call calc_sp_mb_dn_via_node2(mb_der3,7,mb_xyzsp,2,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,24,3) ! (y*z_kc)_ct
    call calc_sp_mb_dn_via_node2(mb_der3,2,mb_xyzsp,3,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,25,1) ! (z*x_et)_kc
    call calc_sp_mb_dn_via_node2(mb_der3,3,mb_xyzsp,3,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,26,2) ! (z*x_ct)_et
    call calc_sp_mb_dn_via_node2(mb_der3,1,mb_xyzsp,3,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,27,3) ! (z*x_kc)_ct

    call calc_sp_mb_dn_via_node2(mb_der3,6,mb_xyzsp,1,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,28,1) ! (x*y_ct)_kc
    call calc_sp_mb_dn_via_node2(mb_der3,4,mb_xyzsp,1,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,29,2) ! (x*y_kc)_et
    call calc_sp_mb_dn_via_node2(mb_der3,5,mb_xyzsp,1,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,30,3) ! (x*y_et)_ct
    call calc_sp_mb_dn_via_node2(mb_der3,9,mb_xyzsp,2,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,31,1) ! (y*z_ct)_kc
    call calc_sp_mb_dn_via_node2(mb_der3,7,mb_xyzsp,2,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,32,2) ! (y*z_kc)_et
    call calc_sp_mb_dn_via_node2(mb_der3,8,mb_xyzsp,2,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,33,3) ! (y*z_et)_ct
    call calc_sp_mb_dn_via_node2(mb_der3,3,mb_xyzsp,3,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,34,1) ! (z*x_ct)_kc
    call calc_sp_mb_dn_via_node2(mb_der3,1,mb_xyzsp,3,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,35,2) ! (z*x_kc)_et
    call calc_sp_mb_dn_via_node2(mb_der3,2,mb_xyzsp,3,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der2,36,3) ! (z*x_et)_ct

    call calc_sp_mb_sxyz_vol(mb_der3,1,9,mb_der2,1,36,mb_sxyzsp,1,9,mb_volsp,1)

!!    call openfile(101,"fort.101","unknown","unformatted","stream")
!!    call gather_output_mb_var(101,mb_sxyzsp,1,9,1)
!!    call closefile(101,"fort.101")
!!    stop

    call mb_var_delete(mb_der2)
    call mb_var_delete(mb_der3)

    call exchange_sp_bc_sxyz(mb_sxyzsp,1,9,nghn)

!!------------------------------------------------------------------------------------------------------------------------------------
!    call mb_var_create(mb_der3,1, 9,nghn)
!
!    call calc_mb_dn_via_node3(mb_sxyzsp,1,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d3int,nfsf_grd_d3der,mb_der3,1,1)   ! kcx_kc
!    call calc_mb_dn_via_node3(mb_sxyzsp,4,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d3int,nfsf_grd_d3der,mb_der3,2,2)   ! etx_et
!    call calc_mb_dn_via_node3(mb_sxyzsp,7,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d3int,nfsf_grd_d3der,mb_der3,3,3)   ! ctx_ct
!    call calc_mb_dn_via_node3(mb_sxyzsp,2,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d3int,nfsf_grd_d3der,mb_der3,4,1)   ! kcy_kc
!    call calc_mb_dn_via_node3(mb_sxyzsp,5,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d3int,nfsf_grd_d3der,mb_der3,5,2)   ! ety_et
!    call calc_mb_dn_via_node3(mb_sxyzsp,8,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d3int,nfsf_grd_d3der,mb_der3,6,3)   ! cty_ct
!    call calc_mb_dn_via_node3(mb_sxyzsp,3,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d3int,nfsf_grd_d3der,mb_der3,7,1)   ! kcz_kc
!    call calc_mb_dn_via_node3(mb_sxyzsp,6,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d3int,nfsf_grd_d3der,mb_der3,8,2)   ! etz_et
!    call calc_mb_dn_via_node3(mb_sxyzsp,9,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d3int,nfsf_grd_d3der,mb_der3,9,3)   ! ctz_ct
!
!    call openfile(101,"fort.101","unknown","unformatted","stream")
!    call gather_output_mb_var(101,mb_der3,1,9,1)
!    call closefile(101,"fort.101")
!
!    call mb_var_delete(mb_der3)
!!------------------------------------------------------------------------------------------------------------------------------------

    call mb_var_create_sp(mb_der1,1,3,nghn)
    call calc_sp_mb_dn_via_vec_node2(mb_sxyzsp,1,3,mb_xyzsp,1,3,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der1,1,1)
    call calc_sp_mb_dn_via_vec_node2(mb_sxyzsp,4,6,mb_xyzsp,1,3,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der1,2,2)
    call calc_sp_mb_dn_via_vec_node2(mb_sxyzsp,7,9,mb_xyzsp,1,3,nghn,sub_ve2,nghe,sub_dn2,nfsf_grd_d2int,nfsf_grd_d2der,mb_der1,3,3)
    call calc_sp_mb_vol(mb_der1,1,3,mb_volsp,1)
    call mb_var_delete(mb_der1)

    call pre_exchange_bc_var_sp(mb_volsp,1,1,nghn,nbc_inter_buf_dyn,nsgl_aver_art)
    call post_exchange_bc_var_sp(mb_volsp,1,1,nghn,nbc_inter_buf_dyn,nsgl_aver_art)


    call ghost_bc_var_all_sp(mb_sxyzsp,1,9)
    call ghost_bc_var_all_sp(mb_volsp,1,1)    

end subroutine calc_grid_derivative_cc

subroutine calc_mb_sxyz_vol(mb_der3,nst3,ned3, &
                            mb_der2,nst2,ned2, &
                            mb_sxyz,nst,ned, &
                            mb_vol,nvol)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : sml_vol,half,third
    use mod_datatypes, only : var_block_t,fld_array_t
    use mod_fieldvars, only : nblkcoms,blkcoms
    implicit none
    type(var_block_t), pointer, intent(in) :: mb_der3(:)
    integer(kind_int),          intent(in) :: nst3,ned3
    type(var_block_t), pointer, intent(in) :: mb_der2(:)
    integer(kind_int),          intent(in) :: nst2,ned2
    type(var_block_t), pointer, intent(in) :: mb_sxyz(:)
    integer(kind_int),          intent(in) :: nst,ned
    type(var_block_t), pointer, intent(in) :: mb_vol(:)
    integer(kind_int),          intent(in) :: nvol
    integer(kind_int)           :: nc,nb,ni,nj,nk,i,j,k,m,ierr
    real(kind_real)             :: d3(9),d2(36),sf(9)
    real(kind_real)             :: vkc,vet,vct,vol0
    real(kind_real), external   :: det_mat3x3
    type(fld_array_t), pointer  :: der3(:),der2(:),sxyz(:),vol(:)

    ierr = ned3 - nst3 - 8
    call error_check(ierr, &
                     "The size of array isn't equal to 9 in subroutine calc_mb_sxyz_vol")

    ierr = ned2 - nst2 - 35
    call error_check(ierr, &
                     "The size of array isn't equal to 9 in subroutine calc_mb_sxyz_vol")

    ierr = ned  - nst  - 8
    call error_check(ierr, &
                     "The size of array isn't equal to 9 in subroutine calc_mb_sxyz_vol")


    do nc=1,nblkcoms
        nb = blkcoms(nc)%nb

        ni = blkcoms(nc)%top%nijk(1)
        nj = blkcoms(nc)%top%nijk(2)
        nk = blkcoms(nc)%top%nijk(3)

        der3 => mb_der3(nb)%fld
        der2 => mb_der2(nb)%fld
        sxyz => mb_sxyz(nb)%fld
        vol  => mb_vol(nb)%fld

        do k=1,nk
        do j=1,nj
        do i=1,ni
            do m=nst3,ned3
                d3(m-nst3+1) = der3(m)%r3d(i,j,k)
            end do
            do m=nst2,ned2
                d2(m-nst2+1) = der2(m)%r3d(i,j,k)
            end do

            sf(1) = d2(5) - d2(15) + d2(23) - d2(33)
            sf(4) = d2(6) - d2(13) + d2(24) - d2(31)
            sf(7) = d2(4) - d2(14) + d2(22) - d2(32)

            sf(2) = d2(8) - d2(18) + d2(26) - d2(36)
            sf(5) = d2(9) - d2(16) + d2(27) - d2(34)
            sf(8) = d2(7) - d2(17) + d2(25) - d2(35)

            sf(3) = d2(2) - d2(12) + d2(20) - d2(30)
            sf(6) = d2(3) - d2(10) + d2(21) - d2(28)
            sf(9) = d2(1) - d2(11) + d2(19) - d2(29)

            sf(:) = half*sf(:)

            do m=nst,ned
                sxyz(m)%r3d(i,j,k) = sf(m-nst+1)
            end do

            vkc = sf(1)*d3(1) + sf(4)*d3(2) + sf(7)*d3(3)
            vet = sf(2)*d3(4) + sf(5)*d3(5) + sf(8)*d3(6)
            vct = sf(3)*d3(7) + sf(6)*d3(8) + sf(9)*d3(9)
        
            !!vol0 = det_mat3x3(der3)
            vol0 = third*(vkc + vet + vct)
        
            vol(nvol)%r3d(i,j,k) = vol0
        end do
        end do
        end do
    end do

end subroutine calc_mb_sxyz_vol

subroutine calc_cc_mb_sxyz_vol(mb_der3,nst3,ned3, &
                               mb_der2,nst2,ned2, &
                               mb_sxyzcc,nst,ned, &
                               mb_volcc,nvol)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : sml_vol,half,third
    use mod_datatypes, only : var_block_t,fld_array_t
    use mod_fieldvars, only : nblkcoms,blkcomscc
    implicit none
    type(var_block_t), pointer, intent(in) :: mb_der3(:)
    integer(kind_int),          intent(in) :: nst3,ned3
    type(var_block_t), pointer, intent(in) :: mb_der2(:)
    integer(kind_int),          intent(in) :: nst2,ned2
    type(var_block_t), pointer, intent(in) :: mb_sxyzcc(:)
    integer(kind_int),          intent(in) :: nst,ned
    type(var_block_t), pointer, intent(in) :: mb_volcc(:)
    integer(kind_int),          intent(in) :: nvol
    integer(kind_int)           :: nc,nb,ni,nj,nk,i,j,k,m,ierr
    real(kind_real)             :: d3(9),d2(36),sf(9)
    real(kind_real)             :: vkc,vet,vct,vol0
    real(kind_real), external   :: det_mat3x3
    type(fld_array_t), pointer  :: der3(:),der2(:),sxyz(:),vol(:)

    ierr = ned3 - nst3 - 8
    call error_check(ierr, &
                     "The size of array isn't equal to 9 in subroutine calc_cc_mb_sxyz_vol")

    ierr = ned2 - nst2 - 35
    call error_check(ierr, &
                     "The size of array isn't equal to 9 in subroutine calc_cc_mb_sxyz_vol")

    ierr = ned  - nst  - 8
    call error_check(ierr, &
                     "The size of array isn't equal to 9 in subroutine calc_cc_mb_sxyz_vol")


    do nc=1,nblkcoms
        nb = blkcomscc(nc)%nb

        ni = blkcomscc(nc)%top%nijk(1)
        nj = blkcomscc(nc)%top%nijk(2)
        nk = blkcomscc(nc)%top%nijk(3)

        der3 => mb_der3(nb)%fld
        der2 => mb_der2(nb)%fld
        sxyz => mb_sxyzcc(nb)%fld
        vol  => mb_volcc(nb)%fld

        do k=1,nk
        do j=1,nj
        do i=1,ni
            do m=nst3,ned3
                d3(m-nst3+1) = der3(m)%r3d(i,j,k)
            end do
            do m=nst2,ned2
                d2(m-nst2+1) = der2(m)%r3d(i,j,k)
            end do

            sf(1) = d2(5) - d2(15) + d2(23) - d2(33)
            sf(4) = d2(6) - d2(13) + d2(24) - d2(31)
            sf(7) = d2(4) - d2(14) + d2(22) - d2(32)

            sf(2) = d2(8) - d2(18) + d2(26) - d2(36)
            sf(5) = d2(9) - d2(16) + d2(27) - d2(34)
            sf(8) = d2(7) - d2(17) + d2(25) - d2(35)

            sf(3) = d2(2) - d2(12) + d2(20) - d2(30)
            sf(6) = d2(3) - d2(10) + d2(21) - d2(28)
            sf(9) = d2(1) - d2(11) + d2(19) - d2(29)

            sf(:) = half*sf(:)

            do m=nst,ned
                sxyz(m)%r3d(i,j,k) = sf(m-nst+1)
            end do

            vkc = sf(1)*d3(1) + sf(4)*d3(2) + sf(7)*d3(3)
            vet = sf(2)*d3(4) + sf(5)*d3(5) + sf(8)*d3(6)
            vct = sf(3)*d3(7) + sf(6)*d3(8) + sf(9)*d3(9)
        
            !!vol0 = det_mat3x3(der3)
            vol0 = third*(vkc + vet + vct)
        
            vol(nvol)%r3d(i,j,k) = vol0
        end do
        end do
        end do
    end do

end subroutine calc_cc_mb_sxyz_vol

subroutine calc_sp_mb_sxyz_vol(mb_der3,nst3,ned3, &
                               mb_der2,nst2,ned2, &
                               mb_sxyzsp,nst,ned, &
                               mb_volsp,nvol)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : sml_vol,half,third
    use mod_datatypes, only : var_block_t,fld_array_t
    use mod_fieldvars, only : nblkcoms,blkcomscc
    implicit none
    type(var_block_t), pointer, intent(in) :: mb_der3(:)
    integer(kind_int),          intent(in) :: nst3,ned3
    type(var_block_t), pointer, intent(in) :: mb_der2(:)
    integer(kind_int),          intent(in) :: nst2,ned2
    type(var_block_t), pointer, intent(in) :: mb_sxyzsp(:)
    integer(kind_int),          intent(in) :: nst,ned
    type(var_block_t), pointer, intent(in) :: mb_volsp(:)
    integer(kind_int),          intent(in) :: nvol
    integer(kind_int)           :: nc,nb,ni,nj,nk,i,j,k,m,ierr
    real(kind_real)             :: d3(9),d2(36),sf(9)
    real(kind_real)             :: vkc,vet,vct,vol0
    real(kind_real), external   :: det_mat3x3
    type(fld_array_t), pointer  :: der3(:),der2(:),sxyz(:),vol(:)

    ierr = ned3 - nst3 - 8
    call error_check(ierr, &
                     "The size of array isn't equal to 9 in subroutine calc_sp_mb_sxyz_vol")

    ierr = ned2 - nst2 - 35
    call error_check(ierr, &
                     "The size of array isn't equal to 9 in subroutine calc_sp_mb_sxyz_vol")

    ierr = ned  - nst  - 8
    call error_check(ierr, &
                     "The size of array isn't equal to 9 in subroutine calc_sp_mb_sxyz_vol")


    do nc=1,nblkcoms
        nb = blkcomscc(nc)%nb

        ni = blkcomscc(nc)%top%nijk(1)
        nj = blkcomscc(nc)%top%nijk(2)
        nk = blkcomscc(nc)%top%nijk(3)

        der3 => mb_der3(nb)%fld
        der2 => mb_der2(nb)%fld
        sxyz => mb_sxyzsp(nb)%fld
        vol  => mb_volsp(nb)%fld

        do k=1,nk
        do j=1,nj
        do i=1,ni
            do m=nst3,ned3
                d3(m-nst3+1) = der3(m)%r3d(i,j,k)
            end do
            do m=nst2,ned2
                d2(m-nst2+1) = der2(m)%r3d(i,j,k)
            end do

            sf(1) = d2(5) - d2(15) + d2(23) - d2(33)
            sf(4) = d2(6) - d2(13) + d2(24) - d2(31)
            sf(7) = d2(4) - d2(14) + d2(22) - d2(32)

            sf(2) = d2(8) - d2(18) + d2(26) - d2(36)
            sf(5) = d2(9) - d2(16) + d2(27) - d2(34)
            sf(8) = d2(7) - d2(17) + d2(25) - d2(35)

            sf(3) = d2(2) - d2(12) + d2(20) - d2(30)
            sf(6) = d2(3) - d2(10) + d2(21) - d2(28)
            sf(9) = d2(1) - d2(11) + d2(19) - d2(29)

            sf(:) = half*sf(:)

            do m=nst,ned
                sxyz(m)%r3d(i,j,k) = sf(m-nst+1)
            end do

            vkc = sf(1)*d3(1) + sf(4)*d3(2) + sf(7)*d3(3)
            vet = sf(2)*d3(4) + sf(5)*d3(5) + sf(8)*d3(6)
            vct = sf(3)*d3(7) + sf(6)*d3(8) + sf(9)*d3(9)
        
            !!vol0 = det_mat3x3(der3)
            vol0 = third*(vkc + vet + vct)
        
            vol(nvol)%r3d(i,j,k) = vol0
        end do
        end do
        end do
    end do

end subroutine calc_sp_mb_sxyz_vol

subroutine calc_mb_sxyz_vol_exp(mb_der3,nst3,ned3, &
                                 mb_der2,nst2,ned2, &
                                 mb_sxyz,nst,ned, &
                                 mb_vol,nvol)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : sml_vol,half,third
    use mod_datatypes, only : var_block_t,fld_array_t
    use mod_variables, only : nghnode
    use mod_fieldvars, only : nblkcoms,blkcoms
    implicit none
    type(var_block_t), pointer, intent(in) :: mb_der3(:)
    integer(kind_int),          intent(in) :: nst3,ned3
    type(var_block_t), pointer, intent(in) :: mb_der2(:)
    integer(kind_int),          intent(in) :: nst2,ned2
    type(var_block_t), pointer, intent(in) :: mb_sxyz(:)
    integer(kind_int),          intent(in) :: nst,ned
    type(var_block_t), pointer, intent(in) :: mb_vol(:)
    integer(kind_int),          intent(in) :: nvol
    integer(kind_int)           :: nc,nb,ni,nj,nk,i,j,k,m,ierr
    real(kind_real)             :: d3(9),d2(36),sf(9)
    real(kind_real)             :: vkc,vet,vct,vol0
    real(kind_real), external   :: det_mat3x3
    type(fld_array_t), pointer  :: der3(:),der2(:),sxyz(:),vol(:)

    ierr = ned3 - nst3 - 8
    call error_check(ierr, &
                     "The size of array isn't equal to 9 in subroutine calc_mb_sxyz_vol")

    ierr = ned2 - nst2 - 35
    call error_check(ierr, &
                     "The size of array isn't equal to 9 in subroutine calc_mb_sxyz_vol")

    ierr = ned  - nst  - 8
    call error_check(ierr, &
                     "The size of array isn't equal to 9 in subroutine calc_mb_sxyz_vol")


    do nc=1,nblkcoms
        nb = blkcoms(nc)%nb

        ni = blkcoms(nc)%top%nijk(1)
        nj = blkcoms(nc)%top%nijk(2)
        nk = blkcoms(nc)%top%nijk(3)

        der3 => mb_der3(nb)%fld
        der2 => mb_der2(nb)%fld
        sxyz => mb_sxyz(nb)%fld
        vol  => mb_vol(nb)%fld

        do k=1-nghnode,nk+nghnode
        do j=1-nghnode,nj+nghnode
        do i=1-nghnode,ni+nghnode
            do m=nst3,ned3
                d3(m-nst3+1) = der3(m)%r3d(i,j,k)
            end do
            do m=nst2,ned2
                d2(m-nst2+1) = der2(m)%r3d(i,j,k)
            end do

            sf(1) = d2(5) - d2(15) + d2(23) - d2(33)
            sf(4) = d2(6) - d2(13) + d2(24) - d2(31)
            sf(7) = d2(4) - d2(14) + d2(22) - d2(32)

            sf(2) = d2(8) - d2(18) + d2(26) - d2(36)
            sf(5) = d2(9) - d2(16) + d2(27) - d2(34)
            sf(8) = d2(7) - d2(17) + d2(25) - d2(35)

            sf(3) = d2(2) - d2(12) + d2(20) - d2(30)
            sf(6) = d2(3) - d2(10) + d2(21) - d2(28)
            sf(9) = d2(1) - d2(11) + d2(19) - d2(29)

            sf(:) = half*sf(:)

            do m=nst,ned
                sxyz(m)%r3d(i,j,k) = sf(m-nst+1)
            end do

            vkc = sf(1)*d3(1) + sf(4)*d3(2) + sf(7)*d3(3)
            vet = sf(2)*d3(4) + sf(5)*d3(5) + sf(8)*d3(6)
            vct = sf(3)*d3(7) + sf(6)*d3(8) + sf(9)*d3(9)
        
            !!vol0 = det_mat3x3(der3)
            vol0 = third*(vkc + vet + vct)
        
            vol(nvol)%r3d(i,j,k) = vol0
        end do
        end do
        end do
    end do

end subroutine calc_mb_sxyz_vol_exp

function det_mat3x3(mat) result(det)
    use mod_kndconsts, only : kind_real
    implicit none
    real(kind_real), intent(in)  :: mat(9)
    real(kind_real)              :: det

    ! | mat(1) mat(2) mat(3) |
    ! | mat(4) mat(5) mat(6) |
    ! | mat(7) mat(8) mat(9) |
    det = mat(1)*mat(5)*mat(9) + mat(2)*mat(6)*mat(7) + &
          mat(3)*mat(8)*mat(4) - mat(3)*mat(5)*mat(7) - &
          mat(2)*mat(4)*mat(9) - mat(1)*mat(8)*mat(6)

end function det_mat3x3

subroutine calc_mb_vol(mb_der1,nst1,ned1,mb_vol,nvol)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : sml_vol,third
    use mod_datatypes, only : var_block_t,fld_array_t
    use mod_fieldvars, only : nblkcoms,blkcoms
    implicit none
    type(var_block_t), pointer, intent(in) :: mb_der1(:)
    integer(kind_int),          intent(in) :: nst1,ned1
    type(var_block_t), pointer, intent(in) :: mb_vol(:)
    integer(kind_int),          intent(in) :: nvol
    integer(kind_int)           :: nc,nb,ni,nj,nk,i,j,k,m,ierr
    real(kind_real)             :: d1(3),vol0
    type(fld_array_t), pointer  :: der1(:),vol(:)

    ierr = ned1 - nst1 - 2
    call error_check(ierr, &
                     "The size of array isn't equal to 3 in subroutine calc_mb_vol")

    do nc=1,nblkcoms
        nb = blkcoms(nc)%nb

        ni = blkcoms(nc)%top%nijk(1)
        nj = blkcoms(nc)%top%nijk(2)
        nk = blkcoms(nc)%top%nijk(3)

        der1 => mb_der1(nb)%fld
        vol  => mb_vol(nb)%fld

        do k=1,nk
        do j=1,nj
        do i=1,ni
            do m=nst1,ned1
                d1(m-nst1+1) = der1(m)%r3d(i,j,k)
            end do

            vol0 = third*(d1(1) + d1(2) + d1(3))

            if ( vol0 < sml_vol ) then
                print '(a5,4(1x,i5),1x,e12.5)',"V3<0:",nb,i,j,k,vol0
            end if

            vol(nvol)%r3d(i,j,k) = vol0
        end do
        end do
        end do
    end do

end subroutine calc_mb_vol

subroutine calc_cc_mb_vol(mb_der1,nst1,ned1,mb_volcc,nvol)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : sml_vol,third
    use mod_datatypes, only : var_block_t,fld_array_t
    use mod_fieldvars, only : nblkcoms,blkcomscc
    implicit none
    type(var_block_t), pointer, intent(in) :: mb_der1(:)
    integer(kind_int),          intent(in) :: nst1,ned1
    type(var_block_t), pointer, intent(in) :: mb_volcc(:)
    integer(kind_int),          intent(in) :: nvol
    integer(kind_int)           :: nc,nb,ni,nj,nk,i,j,k,m,ierr
    real(kind_real)             :: d1(3),vol0
    type(fld_array_t), pointer  :: der1(:),vol(:)

    ierr = ned1 - nst1 - 2
    call error_check(ierr, &
                     "The size of array isn't equal to 3 in subroutine calc_cc_mb_vol")

    do nc=1,nblkcoms
        nb = blkcomscc(nc)%nb

        ni = blkcomscc(nc)%top%nijk(1)
        nj = blkcomscc(nc)%top%nijk(2)
        nk = blkcomscc(nc)%top%nijk(3)

        der1 => mb_der1(nb)%fld
        vol  => mb_volcc(nb)%fld

        do k=1,nk
        do j=1,nj
        do i=1,ni
            do m=nst1,ned1
                d1(m-nst1+1) = der1(m)%r3d(i,j,k)
            end do

            vol0 = third*(d1(1) + d1(2) + d1(3))

            if ( vol0 < sml_vol ) then
                print '(a5,4(1x,i5),1x,e12.5)',"V3<0:",nb,i,j,k,vol0
            end if

            vol(nvol)%r3d(i,j,k) = vol0
        end do
        end do
        end do
    end do

end subroutine calc_cc_mb_vol

subroutine calc_sp_mb_vol(mb_der1,nst1,ned1,mb_volsp,nvol)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : sml_vol,third
    use mod_datatypes, only : var_block_t,fld_array_t
    use mod_fieldvars, only : nblkcoms,blkcomssp
    implicit none
    type(var_block_t), pointer, intent(in) :: mb_der1(:)
    integer(kind_int),          intent(in) :: nst1,ned1
    type(var_block_t), pointer, intent(in) :: mb_volsp(:)
    integer(kind_int),          intent(in) :: nvol
    integer(kind_int)           :: nc,nb,ni,nj,nk,i,j,k,m,ierr
    real(kind_real)             :: d1(3),vol0
    type(fld_array_t), pointer  :: der1(:),vol(:)

    ierr = ned1 - nst1 - 2
    call error_check(ierr, &
                     "The size of array isn't equal to 3 in subroutine calc_sp_mb_vol")

    do nc=1,nblkcoms
        nb = blkcomssp(nc)%nb

        ni = blkcomssp(nc)%top%nijk(1)
        nj = blkcomssp(nc)%top%nijk(2)
        nk = blkcomssp(nc)%top%nijk(3)

        der1 => mb_der1(nb)%fld
        vol  => mb_volsp(nb)%fld

        do k=1,nk
        do j=1,nj
        do i=1,ni
            do m=nst1,ned1
                d1(m-nst1+1) = der1(m)%r3d(i,j,k)
            end do

            vol0 = third*(d1(1) + d1(2) + d1(3))

            if ( vol0 < sml_vol ) then
                print '(a5,4(1x,i5),1x,e12.5)',"V3<0:",nb,i,j,k,vol0
            end if

            vol(nvol)%r3d(i,j,k) = vol0
        end do
        end do
        end do
    end do

end subroutine calc_sp_mb_vol

subroutine reset_grid_derivative
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : small,large    
    use mod_constants, only : mide,sml_vol,sml_ssf
    use mod_constants, only : pole_ssf,pole_vol,bc_pole
    use mod_datatypes, only : top_block_t,bc_region_t,fld_array_t
    use mod_fieldvars, only : nblkcoms,blkcoms
    use mod_fieldvars, only : mb_sxyz,mb_vol
    implicit none
    integer(kind_int)          :: nc,nb,ni,nj,nk,i,j,k,m,ierr
    integer(kind_int)          :: nr,nregs,bctype,m1,m2,m3
    integer(kind_int)          :: s_st(3),s_ed(3),s_nd,s_lr
    integer(kind_int)          :: s_lr3d(3),ijkb(3),ijk2(3)
    real(kind_real)            :: sf(9),vol0,nx,ny,nz,on
    type(top_block_t), pointer :: top
    type(bc_region_t), pointer :: reg
    type(fld_array_t), pointer :: sxyz(:),vol(:) 

    do nc=1,nblkcoms
        nb  =  blkcoms(nc)%nb
        top => blkcoms(nc)%top

        ni = top%nijk(1)
        nj = top%nijk(2)
        nk = top%nijk(3)

        sxyz => mb_sxyz(nb)%fld
        vol  => mb_vol(nb)%fld

        nregs = top%nregions
        do nr=1,nregs
            reg => top%bcs(nr)
            s_st(:) = reg%s_st(:)
            s_ed(:) = reg%s_ed(:)
            s_nd    = reg%s_nd
            s_lr    = reg%s_lr

            s_lr3d(:) = s_lr*mide(:,s_nd)

            bctype = reg%bctype
            if (bctype == bc_pole) then
                do k=s_st(3),s_ed(3)
                do j=s_st(2),s_ed(2)
                do i=s_st(1),s_ed(1)
                    ijkb(:) = (/i,j,k/)
                    ijk2(:) = ijkb(:) - s_lr3d(:)

                    do m=1,3
                        m1 = 3*m - 2
                        m2 = m1 + 1
                        m3 = m2 + 1

                        nx = sxyz(m1)%r3d(ijk2(1),ijk2(2),ijk2(3))
                        ny = sxyz(m2)%r3d(ijk2(1),ijk2(2),ijk2(3))
                        nz = sxyz(m3)%r3d(ijk2(1),ijk2(2),ijk2(3))
                        !!on = pole_ssf/max(sqrt(nx*nx + ny*ny + nz*nz),sml_ssf)
                        on = pole_ssf
                        nx = nx*on
                        ny = ny*on
                        nz = nz*on
                        sxyz(m1)%r3d(i,j,k) = nx
                        sxyz(m2)%r3d(i,j,k) = ny
                        sxyz(m3)%r3d(i,j,k) = nz
                    end do

                    vol0 = vol(1)%r3d(ijk2(1),ijk2(2),ijk2(3))

                    !!vol(1)%r3d(i,j,k) = pole_vol
                    vol(1)%r3d(i,j,k) = pole_vol*vol0
                end do
                end do
                end do
            end if
        end do

        do k=1,nk
        do j=1,nj
        do i=1,ni
            vol0 = vol(1)%r3d(i,j,k)

            if ( vol0 < sml_vol ) then
                print '(a5,4(1x,i5),1x,e12.5)',"V2<0:",nb,i,j,k,vol0
            end if
        end do
        end do
        end do       
    end do

end subroutine reset_grid_derivative

subroutine reset_grid_derivative_cc
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : small,large    
    use mod_constants, only : mide,sml_vol,sml_ssf
    use mod_constants, only : pole_ssf,pole_vol,bc_pole
    use mod_datatypes, only : top_block_t,bc_region_t,fld_array_t
    use mod_fieldvars, only : nblkcoms,blkcomscc,blkcomssp
    use mod_fieldvars, only : mb_sxyzcc,mb_volcc,mb_sxyzsp,mb_volsp
    implicit none
    integer(kind_int)          :: nc,nb,ni,nj,nk,i,j,k,m,ierr
    integer(kind_int)          :: nr,nregs,bctype,m1,m2,m3
    integer(kind_int)          :: s_st(3),s_ed(3),s_nd,s_lr
    integer(kind_int)          :: s_lr3d(3),ijkb(3),ijk2(3)
    real(kind_real)            :: sf(9),vol0,nx,ny,nz,on
    type(top_block_t), pointer :: top
    type(bc_region_t), pointer :: reg
    type(fld_array_t), pointer :: sxyz(:),vol(:) 

    do nc=1,nblkcoms
        nb  =  blkcomscc(nc)%nb
        top => blkcomscc(nc)%top

        ni = top%nijk(1)
        nj = top%nijk(2)
        nk = top%nijk(3)

        sxyz => mb_sxyzcc(nb)%fld
        vol  => mb_volcc(nb)%fld

        nregs = top%nregions
        do nr=1,nregs
            reg => top%bcs(nr)
            s_st(:) = reg%s_st(:)
            s_ed(:) = reg%s_ed(:)
            s_nd    = reg%s_nd
            s_lr    = reg%s_lr

            s_lr3d(:) = s_lr*mide(:,s_nd)

            bctype = reg%bctype
            if (bctype == bc_pole) then
                do k=s_st(3),s_ed(3)
                do j=s_st(2),s_ed(2)
                do i=s_st(1),s_ed(1)
                    ijkb(:) = (/i,j,k/)
                    ijk2(:) = ijkb(:) - s_lr3d(:)

                    do m=1,3
                        m1 = 3*m - 2
                        m2 = m1 + 1
                        m3 = m2 + 1

                        nx = sxyz(m1)%r3d(ijk2(1),ijk2(2),ijk2(3))
                        ny = sxyz(m2)%r3d(ijk2(1),ijk2(2),ijk2(3))
                        nz = sxyz(m3)%r3d(ijk2(1),ijk2(2),ijk2(3))
                        !!on = pole_ssf/max(sqrt(nx*nx + ny*ny + nz*nz),sml_ssf)
                        on = pole_ssf
                        nx = nx*on
                        ny = ny*on
                        nz = nz*on
                        sxyz(m1)%r3d(i,j,k) = nx
                        sxyz(m2)%r3d(i,j,k) = ny
                        sxyz(m3)%r3d(i,j,k) = nz
                    end do

                    vol0 = vol(1)%r3d(ijk2(1),ijk2(2),ijk2(3))

                    !!vol(1)%r3d(i,j,k) = pole_vol
                    vol(1)%r3d(i,j,k) = pole_vol*vol0
                end do
                end do
                end do
            end if
        end do

        do k=1,nk
        do j=1,nj
        do i=1,ni
            vol0 = vol(1)%r3d(i,j,k)

            if ( vol0 < sml_vol ) then
                print '(a5,4(1x,i5),1x,e12.5)',"Cell Edge Based V2<0:",nb,i,j,k,vol0
            end if
        end do
        end do
        end do       
    end do
    
    do nc=1,nblkcoms
        nb  =  blkcomssp(nc)%nb
        top => blkcomssp(nc)%top

        ni = top%nijk(1)
        nj = top%nijk(2)
        nk = top%nijk(3)

        sxyz => mb_sxyzsp(nb)%fld
        vol  => mb_volsp(nb)%fld

        nregs = top%nregions
        do nr=1,nregs
            reg => top%bcs(nr)
            s_st(:) = reg%s_st(:)
            s_ed(:) = reg%s_ed(:)
            s_nd    = reg%s_nd
            s_lr    = reg%s_lr

            s_lr3d(:) = s_lr*mide(:,s_nd)

            bctype = reg%bctype
            if (bctype == bc_pole) then
                do k=s_st(3),s_ed(3)
                do j=s_st(2),s_ed(2)
                do i=s_st(1),s_ed(1)
                    ijkb(:) = (/i,j,k/)
                    ijk2(:) = ijkb(:) - s_lr3d(:)

                    do m=1,3
                        m1 = 3*m - 2
                        m2 = m1 + 1
                        m3 = m2 + 1

                        nx = sxyz(m1)%r3d(ijk2(1),ijk2(2),ijk2(3))
                        ny = sxyz(m2)%r3d(ijk2(1),ijk2(2),ijk2(3))
                        nz = sxyz(m3)%r3d(ijk2(1),ijk2(2),ijk2(3))
                        !!on = pole_ssf/max(sqrt(nx*nx + ny*ny + nz*nz),sml_ssf)
                        on = pole_ssf
                        nx = nx*on
                        ny = ny*on
                        nz = nz*on
                        sxyz(m1)%r3d(i,j,k) = nx
                        sxyz(m2)%r3d(i,j,k) = ny
                        sxyz(m3)%r3d(i,j,k) = nz
                    end do

                    vol0 = vol(1)%r3d(ijk2(1),ijk2(2),ijk2(3))

                    !!vol(1)%r3d(i,j,k) = pole_vol
                    vol(1)%r3d(i,j,k) = pole_vol*vol0
                end do
                end do
                end do
            end if
        end do

        do k=1,nk
        do j=1,nj
        do i=1,ni
            vol0 = vol(1)%r3d(i,j,k)

            if ( vol0 < sml_vol ) then
                print '(a5,4(1x,i5),1x,e12.5)',"Cell Center Based V2<0:",nb,i,j,k,vol0
            end if
        end do
        end do
        end do       
    end do    

end subroutine reset_grid_derivative_cc

subroutine set_wall_dist
    use mod_kndconsts, only : kind_int
    use mod_constants, only : nbc_inter_buf_dyn
    use mod_constants, only : nsgl_aver_art,bc_wall
    use mod_constants, only : blk_mark_dstmin
    use mod_variables, only : nrddst,nghnode,nsw_dst
    use mod_fieldvars, only : nblocks,mb_top
    use mod_fieldvars, only : mb_dst,wallpnts
    use mod_interface, only : calc_wall_dist
    use mod_interface, only : pre_exchange_bc_var,post_exchange_bc_var
    implicit none
    integer(kind_int) :: nb,ierr

    if (nsw_dst > 0) then
        do nb=1,nblocks
            mb_top(nb)%mark = blk_mark_dstmin
        end do

        call exchange_wall_pnts(bc_wall,-1)

        if (nrddst == 0) then
            call calc_wall_dist(mb_dst,bc_wall,-1)
            call output_dst
            call msg_seq_and_master("Compute the wall distances successfully")
        else
            call input_dst
            call msg_seq_and_master("Input the wall distances successfully")
        end if

        call pre_exchange_bc_var(mb_dst,1,2,nghnode,nbc_inter_buf_dyn,nsgl_aver_art)
        call post_exchange_bc_var(mb_dst,1,2,nghnode,nbc_inter_buf_dyn,nsgl_aver_art)

        deallocate(wallpnts, stat=ierr)

    end if

end subroutine set_wall_dist

subroutine calc_wall_dist(mb_dst,bc_target,mark)
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : small,large,tenth
    use mod_datatypes, only : fld_array_t,var_block_t
    use mod_datatypes, only : top_block_t,bc_region_t
    use mod_fieldvars, only : mb_top,mb_xyz
    use mod_fieldvars, only : nblkcoms,blkcoms
    use mod_fieldvars, only : nwallpnts,wallpnts
    implicit none
    type(var_block_t), pointer, intent(in) :: mb_dst(:)
    integer(kind_int),          intent(in) :: bc_target,mark
    integer(kind_int)          :: nc,nb,ni,nj,nk,i,j,k,ip,ierr
    integer(kind_int)          :: nr,nregs,s_st(3),s_ed(3)
    integer(kind_int)          :: s_nd,s_lr,bctype,bcmark
    real(kind_real)            :: xw,yw,zw,dx,dy,dz,dist
    logical                    :: istar
    type(top_block_t), pointer :: top
    type(bc_region_t), pointer :: reg
    type(fld_array_t), pointer :: xyz(:),dst(:)

    do nc=1,nblkcoms
        nb  =  blkcoms(nc)%nb
        top => blkcoms(nc)%top

!!        if (top%mark > 0) then
            ni = top%nijk(1)
            nj = top%nijk(2)
            nk = top%nijk(3)

            dst => mb_dst(nb)%fld

            do k=1,nk
            do j=1,nj
            do i=1,ni
               dst(1)%r3d(i,j,k) = large
               dst(2)%r3d(i,j,k) = tenth
            end do
            end do
            end do
!!        end if
    end do

    do ip=1,nwallpnts
        xw = wallpnts(ip)%xw
        yw = wallpnts(ip)%yw
        zw = wallpnts(ip)%zw

        do nc=1,nblkcoms
            nb  =  blkcoms(nc)%nb
            top => blkcoms(nc)%top

            if (top%mark > 0) then
                ni = top%nijk(1)
                nj = top%nijk(2)
                nk = top%nijk(3)

                xyz => mb_xyz(nb)%fld
                dst => mb_dst(nb)%fld

                do k=1,nk
                do j=1,nj
                do i=1,ni
                    dx = xyz(1)%r3d(i,j,k) - xw
                    dy = xyz(2)%r3d(i,j,k) - yw
                    dz = xyz(3)%r3d(i,j,k) - zw
                    dist = dx*dx + dy*dy + dz*dz

                    if (dst(1)%r3d(i,j,k) >= dist ) then
                       dst(1)%r3d(i,j,k) = dist
                       dst(2)%r3d(i,j,k) = ip + tenth
                    end if
                end do
                end do
                end do
            end if
        end do
    end do

    do nc=1,nblkcoms
        nb  =  blkcoms(nc)%nb
        top => blkcoms(nc)%top

        if (top%mark > 0) then
            dst => mb_dst(nb)%fld

            nregs = top%nregions
            do nr=1,nregs
                reg => top%bcs(nr)

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
                        dst(1)%r3d(i,j,k) = small
                    end do
                    end do
                    end do
                end if
            end do
        end if
    end do

    do nc=1,nblkcoms
        nb  =  blkcoms(nc)%nb
        top => blkcoms(nc)%top

        if (top%mark > 0) then
            ni = top%nijk(1)
            nj = top%nijk(2)
            nk = top%nijk(3)

            dst => mb_dst(nb)%fld

            do k=1,nk
            do j=1,nj
            do i=1,ni
                dist = dst(1)%r3d(i,j,k)
                dst(1)%r3d(i,j,k) = sqrt(dist)
            end do
            end do
            end do
        end if
    end do

end subroutine calc_wall_dist

subroutine set_sponge_dist
    use mod_kndconsts, only : kind_int
    use mod_constants, only : bc_farfield,bc_cut1to1
    use mod_constants, only : blk_mark_sponge,bc_mark_sponge
    use mod_datatypes, only : top_block_t,bc_region_t
    use mod_variables, only : nrddsp,nghnode,nsponge
    use mod_fieldvars, only : nblocks,mb_top
    use mod_fieldvars, only : mb_dsp,wallpnts
    use mod_interface, only : calc_wall_dist
    implicit none
    integer(kind_int)          :: nb,nr,nbt
    integer(kind_int)          :: nregs,bctype,ierr
    type(top_block_t), pointer :: top
    type(bc_region_t), pointer :: reg

    if (nsponge > 0) then
        do nb=1,nblocks
            mb_top(nb)%mark = 0

            nregs = mb_top(nb)%nregions
            do nr=1,nregs
                reg => mb_top(nb)%bcs(nr)

                bctype = reg%bctype
                if (bctype == bc_farfield) then
                    mb_top(nb)%mark = blk_mark_sponge
                    exit
                end if
            end do
        end do

        do nb=1,nblocks
            nregs = mb_top(nb)%nregions
            do nr=1,nregs
                reg => mb_top(nb)%bcs(nr)

                nbt = reg%nbt

                bctype = reg%bctype
                if (bctype == bc_cut1to1) then
                    if (mb_top(nb)%mark == blk_mark_sponge .and. &
                        mb_top(nbt)%mark == 0) then
                        reg%mark = bc_mark_sponge
                    else
                        reg%mark = 0
                    end if
                else
                    reg%mark = 0
                end if
            end do
        end do

        call exchange_wall_pnts(bc_cut1to1,bc_mark_sponge)

        if (nrddsp == 0) then
            call calc_wall_dist(mb_dsp,bc_cut1to1,bc_mark_sponge)
            call scale_sponge_dist
            call output_dsp
            call msg_seq_and_master("Compute the sponge distances successfully")
        else
            call input_dsp
            call msg_seq_and_master("Input the sponge distances successfully")
        endif

        deallocate(wallpnts, stat=ierr)

    end if

end subroutine set_sponge_dist

subroutine scale_sponge_dist
    use mod_kndconsts, only : kind_int,kind_real
    use mod_constants, only : small,large,bc_farfield
    use mod_constants, only : blk_mark_sponge,io_unit_dsp
    use mod_datatypes, only : fld_array_t,var_block_t
    use mod_datatypes, only : top_block_t,bc_region_t
    use mod_variables, only : dspmin,dspmax,grdfile
    use mod_fieldvars, only : mb_top,mb_xyz,mb_dsp
    use mod_fieldvars, only : nblkcoms,blkcoms
    use mod_fieldvars, only : nwallpnts,wallpnts
    implicit none
    integer(kind_int)          :: nc,nb,ni,nj,nk,i,j,k,ierr
    integer(kind_int)          :: nr,nregs,s_st(3),s_ed(3)
    integer(kind_int)          :: s_nd,s_lr,bctype
    real(kind_real)            :: disp
    type(top_block_t), pointer :: top
    type(bc_region_t), pointer :: reg
    type(fld_array_t), pointer :: dsp(:)

    dspmin = large
    dspmax = small
    do nc=1,nblkcoms
        nb  =  blkcoms(nc)%nb
        top => blkcoms(nc)%top

        if (top%mark == blk_mark_sponge) then
            dsp => mb_dsp(nb)%fld

            nregs = top%nregions
            do nr=1,nregs
                reg => top%bcs(nr)

                bctype = reg%bctype
                if (bctype == bc_farfield) then
                    s_st(:) = reg%s_st(:)
                    s_ed(:) = reg%s_ed(:)
                    s_nd    = reg%s_nd
                    s_lr    = reg%s_lr

                    do k=s_st(3),s_ed(3)
                    do j=s_st(2),s_ed(2)
                    do i=s_st(1),s_ed(1)
                        disp = dsp(1)%r3d(i,j,k)
                        dspmin = min(dspmin,disp)
                        dspmax = max(dspmax,disp)
                    end do
                    end do
                    end do
                end if
            end do
        end if
    end do

    call broadcast_sponge

    do nc=1,nblkcoms
        nb  =  blkcoms(nc)%nb
        top => blkcoms(nc)%top

        if (top%mark == blk_mark_sponge) then
            ni = top%nijk(1)
            nj = top%nijk(2)
            nk = top%nijk(3)

            dsp => mb_dsp(nb)%fld

            do k=1,nk
            do j=1,nj
            do i=1,ni
                disp = dsp(1)%r3d(i,j,k)
                dsp(1)%r3d(i,j,k) = min(disp/dspmin,3.0)
            end do
            end do
            end do

        end if
    end do

end subroutine scale_sponge_dist


