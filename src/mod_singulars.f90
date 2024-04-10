
module mod_singulars
    use mod_kndconsts, only : kind_int,kind_real

    type box_t
        integer(kind_int) :: st(3)
        integer(kind_int) :: ed(3)
    end type box_t

    type box_node_t
        type(box_t)                 :: box
        integer(kind_int)           :: nblk
        logical                     :: del
        logical                     :: visited
        integer(kind_int)           :: dir
        type(box_node_t), pointer   :: lnode    ! 二叉树的左右子结点
        type(box_node_t), pointer   :: rnode
        type(box_node_t), pointer   :: prev     ! 双向列表的前后继结点
        type(box_node_t), pointer   :: next
    end type box_node_t

    type box_node_ptr_t
        type(box_node_t), pointer     :: ptr
        logical                       :: flag
        integer(kind_int)             :: nlen_segmt
        integer(kind_int)             :: nlen_chain
        type(box_node_ptr_t), pointer :: next
    end type box_node_ptr_t

    ! 物理空间上重合的一组点的计算坐标信息
    type item_list_t
        integer(kind_int)           :: npnts
        ! pnts(1:npnts, 6), 其中 pnts(i, :) = (/nb, si, sj, sk, seqno, valid_flag/)
        integer(kind_int), pointer  :: pnts(:,:)
    end type item_list_t

    ! 奇点数据交换用的缓冲区
    type simp_buf_t
        integer(kind_int)        :: nvar
        integer(kind_int)        :: npnt
        real(kind_real), pointer :: buf(:,:)
    end type simp_buf_t

    integer(kind_int)          :: nsimp_arr
    type(item_list_t), pointer :: simp_arr(:)
    type(simp_buf_t),  pointer :: simp_buf(:, :)
    integer(kind_int), pointer :: nmsg_proc(:, :)  ! nproc * 3, 1=进程内奇点数, 2=当前进程是否需要与各进程交换数据

contains

    ! point class

    function    point_equal(pnt1, pnt2) result(ret)
        implicit none
        integer(kind_int), intent(in) :: pnt1(:), pnt2(:)
        logical :: ret

        ret = all(pnt1(1:3) == pnt2(1:3))
    end function point_equal

    ! box class

    function    box_compare(a, b) result(ret)
        implicit none
        type(box_t) :: a, b
        integer(kind_int) :: ret, i

        do i = 1, 3
            if ( a%st(i) < b%st(i) ) then
                ret = -1
                return
            else if ( a%st(i) > b%st(i) ) then
                ret = 1
                return
            else
                ! do nothing
            end if
        end do
        do i = 1, 3
            if ( a%ed(i) < b%ed(i) ) then
                ret = -1
                return
            else if ( a%ed(i) > b%ed(i) ) then
                ret = 1
                return
            else
                ! do nothing
            end if
        end do
        ret = 0
    end function box_compare

    function    box_volume(box) result(vol)
        implicit none
        type(box_t) :: box
        integer(kind_int) :: vol

        vol = product(box%ed - box%st + 1)
    end function box_volume

    function    box_dim(box) result(ret)
        implicit none
        type(box_t) :: box
        integer(kind_int) :: ret
        integer(kind_int) :: i

        ret = 0
        do i = 1, 3
            if (box%st(i) /= box%ed(i)) then
                ret = ret + 1
            end if
        end do
    end function box_dim

    subroutine  box_make(box, first, second)
        implicit none
        type(box_t) :: box
        integer(kind_int) :: first(3)
        integer(kind_int) :: second(3)
        integer           :: i

        do i = 1, 3
            box%st(i) = min(first(i), second(i))
            box%ed(i) = max(first(i), second(i))
        end do
    end subroutine box_make

    subroutine  edge_box_create(edge_boxes, maxijk)
        implicit none
        type(box_t) :: edge_boxes(12)
        integer(kind_int) :: maxijk(3)

        integer(kind_int) :: imin, jmin, kmin, imax, jmax, kmax

        imin = 1;     imax = maxijk(1)
        jmin = 1;     jmax = maxijk(2)
        kmin = 1;     kmax = maxijk(3)

        call box_create( edge_boxes(11), (/imin, jmin, kmin/), (/imin, jmin, kmax/) )
        call box_create( edge_boxes( 9), (/imin, jmin, kmin/), (/imin, jmax, kmin/) )
        call box_create( edge_boxes( 7), (/imin, jmin, kmax/), (/imin, jmax, kmax/) )
        call box_create( edge_boxes( 5), (/imin, jmax, kmin/), (/imin, jmax, kmax/) )

        call box_create( edge_boxes( 3), (/imin, jmin, kmin/), (/imax, jmin, kmin/) )
        call box_create( edge_boxes( 1), (/imin, jmin, kmax/), (/imax, jmin, kmax/) )
        call box_create( edge_boxes( 2), (/imin, jmax, kmin/), (/imax, jmax, kmin/) )
        call box_create( edge_boxes( 4), (/imin, jmax, kmax/), (/imax, jmax, kmax/) )

        call box_create( edge_boxes( 6), (/imax, jmin, kmin/), (/imax, jmin, kmax/) )
        call box_create( edge_boxes( 8), (/imax, jmin, kmin/), (/imax, jmax, kmin/) )
        call box_create( edge_boxes(10), (/imax, jmin, kmax/), (/imax, jmax, kmax/) )
        call box_create( edge_boxes(12), (/imax, jmax, kmin/), (/imax, jmax, kmax/) )

    end subroutine edge_box_create

    subroutine  box_create(box, st, ed)
        implicit none
        integer(kind_int), intent(in):: st(:), ed(:)
        type(box_t), intent(out):: box

        box%st = st(1:3)
        box%ed = ed(1:3)
    end subroutine box_create

    subroutine  box_copy(box, box_another)
        implicit none
        type(box_t), intent(out):: box
        type(box_t), intent(in) :: box_another

        call box_create(box, box_another%st, box_another%ed)
    end subroutine box_copy

    function    box_dir(box) result(ret)
        implicit none
        type(box_t) :: box
        integer(kind_int) :: ret
        integer(kind_int) :: i

        ret = 0
        do i = 1, 3
            if (box%st(i) /= box%ed(i)) then
                ret = i
                exit
            end if
        end do
    end function box_dir

    ! 两个线段盒子拆分成互不相交的子集之和
    subroutine  box_divide(box1, box2, outbox, nout)
        implicit none
        type(box_t), intent(in)             :: box1, box2
        type(box_t), intent(out)            :: outbox(:)
        integer(kind_int), intent(out)      :: nout

        type(box_t) :: box_int
        logical     :: lft, rht
        integer(kind_int) :: nret

        call box_intersection(box1, box2, box_int)

        if (box_is_empty(box_int)) then
            ! case 1: A, B 无交集
            nout = -2
            call box_copy(outbox(1), box1)
            call box_copy(outbox(2), box2)
            return
        end if

        call box_copy(outbox(1), box_int)
        nout = 1

        lft = box_contains_box(box1, box2)
        rht = box_contains_box(box2, box1)
        if (lft .and. rht) then
            ! case 2: A = B
            return
        end if

        if (lft) then
            ! add A\B to the output set
            call box_diff_subbox(box1, box_int, outbox(nout+1:nout+2), nret)
            nout = nout + nret

        else if (rht) then
            ! add B\A to the output set
            call box_diff_subbox(box2, box_int, outbox(nout+1:nout+2), nret)
            nout = nout + nret

        else
            call box_diff_subbox(box1, box_int, outbox(nout+1:nout+2), nret)
            nout = nout + nret

            call box_diff_subbox(box2, box_int, outbox(nout+1:nout+2), nret)
            nout = nout + nret
        end if
    end subroutine box_divide

    subroutine  box_diff_subbox(box, subbox, outbox, nout)
        implicit none
        type(box_t), intent(in)         :: box, subbox
        type(box_t), intent(out)        :: outbox(:)
        integer(kind_int), intent(out)  :: nout

        logical           :: left, right
        integer(kind_int) :: dir

        dir  = box_dir(box)
        left = point_equal(box%st, subbox%st)
        right= point_equal(box%ed, subbox%ed)

        if (dir == 0) then
            write(*, *) " !! Error: dir = 0 in <! box_diff_subbox !>"
        end if

        if (left .and. right) then
            nout = 0
        else if(left) then
            nout = 1
            call box_create(outbox(1), subbox%ed, box%ed)
            outbox(1)%st(dir) = outbox(1)%st(dir) + 1
        else if(right) then
            nout = 1
            call box_create(outbox(1), box%st, subbox%st)
            outbox(1)%ed(dir) = outbox(1)%ed(dir) - 1
        else
            nout = 2
            call box_create(outbox(1), box%st, subbox%st)
            outbox(1)%ed(dir) = outbox(1)%ed(dir) - 1

            call box_create(outbox(2), subbox%ed, box%ed)
            outbox(2)%st(dir) = outbox(2)%st(dir) + 1
        end if
    end subroutine box_diff_subbox

    subroutine  box_intersection(box1, box2, box)
        implicit none
        type(box_t), intent(in) :: box1, box2
        type(box_t), intent(out):: box
        integer(kind_int) :: i

        do i = 1, 3
            box%st(i) = max(box1%st(i), box2%st(i))
            box%ed(i) = min(box1%ed(i), box2%ed(i))
        end do
    end subroutine box_intersection

    function    box_is_empty(box) result(is_empty)
        implicit none
        type(box_t), intent(in):: box
        logical :: is_empty
        integer(kind_int) :: i

        is_empty = .false.
        do i = 1, 3
            if (box%st(i) > box%ed(i)) then
                is_empty = .true.
                exit
            end if
        end do
    end function box_is_empty

    function    box_contains_pnt(box, pnt) result(ret)
        implicit none
        type(box_t), intent(in):: box
        integer(kind_int), intent(in):: pnt(3)
        logical :: ret
        integer(kind_int) :: i

        ret = .true.
        do i = 1, 3
            if (pnt(i) < box%st(i) .or. pnt(i) > box%ed(i)) then
                ret = .false.
                exit
            end if
        end do
    end function box_contains_pnt

    function    box_contains_box(box, box2) result(ret)
        implicit none
        type(box_t), intent(in):: box, box2
        logical :: ret

        ret = (box_contains_pnt(box, box2%st) .and. box_contains_pnt(box, box2%ed))
    end function box_contains_box

    subroutine  box_get_point(box, idx_pnt, pnt)
        implicit none
        type(box_t) :: box
        integer(kind_int) :: idx_pnt
        integer(kind_int), intent(out):: pnt(:)

        integer(kind_int) :: n(3), stride(3), idx

        n = box%ed - box%st + 1
        stride(1) = 1
        stride(2) = n(1) * stride(1)
        stride(3) = n(2) * stride(2)

        idx = idx_pnt - 1
        pnt(3)  = int(idx / stride(3))
        idx     = mod(idx, stride(3))

        pnt(2)  = int(idx / stride(2))
        idx     = mod(idx, stride(2))

        pnt(1)  = idx

        pnt = pnt + box%st
    end subroutine  box_get_point

    ! box node class

    subroutine  box_node_create(node, box, nblk)
        implicit none
        type(box_node_t), pointer :: node
        type(box_t)               :: box
        integer(kind_int)         :: nblk
        integer(kind_int)         :: ierr

        if ( associated(node) ) then
            nullify(node)
        end if
        allocate(node, stat = ierr)

        call box_create(node%box, box%st, box%ed)
        node%nblk = nblk
        node%del = .false.
        node%visited = .false.
        node%dir = 1
        nullify(node%lnode)
        nullify(node%rnode)
        node%prev => node
        node%next => node
    end subroutine box_node_create

    subroutine  box_node_destroy(node, istrue)
        implicit none
        type(box_node_t), pointer :: node
        logical, optional :: istrue

        logical :: ok
        integer(kind_int) :: ierr

        if (present(istrue)) then
            ok = istrue
        else
            ok = .true.
        end if

        if (ok .and. associated(node)) then
            deallocate(node, stat=ierr)
            nullify(node)
        end if
    end subroutine box_node_destroy

    subroutine  box_node_chain_search(head, node, found)
        implicit none
        type(box_node_t), pointer :: head
        type(box_node_t), pointer :: node
        logical, intent(out) :: found

        type(box_node_t), pointer :: curr

        found = associated(head, node)
        if (found) return

        curr => head
        do while (.not. associated(curr%next, head))
            if (associated(curr%next, node)) then
                found = .true.
                return
            end if
            curr => curr%next
        end do
        found = .false.
    end subroutine box_node_chain_search

    subroutine  box_node_chain_flip_dir(head)
        implicit none
        type(box_node_t), pointer :: head
        type(box_node_t), pointer :: curr

        head%dir = - head%dir

        curr => head
        do while (.not. associated(curr%next, head))
            curr%next%dir = - curr%next%dir
            curr => curr%next
        end do
    end subroutine box_node_chain_flip_dir

    subroutine  box_node_chain_set_del(head, stat)
        implicit none
        type(box_node_t), pointer :: head
        logical, optional :: stat

        type(box_node_t), pointer :: curr
        logical :: del_stat

        if (present(stat)) then
            del_stat = stat
        else
            del_stat = .true.
        end if

        if (head%del .eqv. del_stat) return

        curr => head
        curr%del = del_stat
        do while (.not. associated(curr%next, head))
            curr%next%del = del_stat
            curr => curr%next
        end do
    end subroutine box_node_chain_set_del

    subroutine  box_node_chain_set_visited(head, stat)
        implicit none
        type(box_node_t), pointer :: head
        logical, optional :: stat

        type(box_node_t), pointer :: curr
        logical :: visited_stat

        if (present(stat)) then
            visited_stat = stat
        else
            visited_stat = .true.
        end if

        if (head%visited .eqv. visited_stat) return

        curr => head
        curr%visited = visited_stat
        do while (.not. associated(curr%next, head))
            curr%next%visited = visited_stat
            curr => curr%next
        end do
    end subroutine box_node_chain_set_visited

    function    box_node_chain_count(head) result(n)
        implicit none
        type(box_node_t), pointer :: head
        integer(kind_int) :: n

        type(box_node_t), pointer :: curr

        n = 1
        curr => head
        do while (.not. associated(curr%next, head))
            n = n + 1
            curr => curr%next
        end do
    end function box_node_chain_count

    function    box_node_compare(a, b) result(ret)
        implicit none
        type(box_node_t) :: a, b
        integer(kind_int) :: ret

        if (a%nblk < b%nblk) then
            ret = -1
        else if(a%nblk > b%nblk) then
            ret = 1
        else
            ret = box_compare(a%box, b%box)
        end if
    end function box_node_compare

    subroutine  box_node_get_point(box_node, idx_pnt, pnt)
        implicit none
        type(box_node_t), pointer :: box_node
        integer(kind_int) :: idx_pnt
        integer(kind_int), intent(out):: pnt(:)

        integer(kind_int) :: dir, vol, idx

        dir = box_node%dir
        vol = box_volume( box_node%box )
        idx = idx_pnt

        if (dir < 0) idx = vol + 1 - idx

        call box_get_point(box_node%box, idx, pnt)
    end subroutine  box_node_get_point

    ! tree of box node

    recursive   &
    subroutine  tree_insert(root, node, stat, pos)
        implicit none
        type(box_node_t), pointer   :: root
        type(box_node_t), pointer,  intent(in)      :: node
        logical,                    intent(out), optional :: stat
        type(box_node_t), pointer,  intent(out), optional :: pos

        logical :: success
        type(box_node_t), pointer :: ptr_found

        nullify(ptr_found)
        success = .false.
        if ( associated(node) ) then
            if ( associated(root) ) then
                if ( box_node_compare(node, root) < 0 ) then
                    call tree_insert(root%lnode, node, success, ptr_found)
                else if ( box_node_compare(node, root) > 0 ) then
                    call tree_insert(root%rnode, node, success, ptr_found)
                else     ! node 已经存在于树中
                    ptr_found => root
                    success = .false.
                end if
            else
                root        => node
                ptr_found   => root
                success     = .true.
            end if
        end if

        if (present(stat)) then
            stat = success
        end if
        if (present(pos)) then
            pos => ptr_found
        end if
    end subroutine tree_insert

    recursive   &
    subroutine  tree_node_filter(root)
        implicit none
        type(box_node_t), pointer :: root

        integer(kind_int) :: n

        if ( associated(root) ) then
            call tree_node_filter(root%lnode)

            n = box_node_chain_count(root)
            if (n <= 2) then
                call box_node_chain_set_del(root)
            end if

            call tree_node_filter(root%rnode)
        end if
    end subroutine tree_node_filter

    recursive   &
    subroutine  tree_destroy(root)
        implicit none
        type(box_node_t), pointer :: root
        integer :: ierr

        if ( associated(root) ) then
            call tree_destroy(root%lnode)
            call tree_destroy(root%rnode)

            deallocate(root, stat=ierr)
        end if
    end subroutine tree_destroy

    subroutine  list_destroy(bnp_head)
        implicit none
        type(box_node_ptr_t), target :: bnp_head ! 专用头结点，自身不能释放

        type(box_node_ptr_t), pointer :: bnp_this
        integer(kind_int) :: ierr

        do while (.not. associated(bnp_head%next, bnp_head))
            bnp_this => bnp_head%next
            bnp_head%next => bnp_this%next
            deallocate(bnp_this, stat=ierr)
        end do
    end subroutine  list_destroy

    ! 遍历 subtree 中所有结点，创建相应 box_node_ptr 包装结点并追加在 tail 之后
    recursive   &
    subroutine  list_append(tail, subtree)
        implicit none
        type(box_node_ptr_t), pointer :: tail
        type(box_node_t), pointer, intent(in) :: subtree

        type(box_node_ptr_t), pointer :: node
        integer :: ierr

        if ( associated(subtree) ) then
            call list_append(tail, subtree%lnode)

            allocate(node, stat=ierr)
            node%ptr => subtree
            node%next => tail%next
            node%flag = .false.
            tail%next => node

            tail => tail%next

            call list_append(tail, subtree%rnode)
        end if

    end subroutine list_append

end module mod_singulars

