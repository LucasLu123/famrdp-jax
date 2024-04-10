module mod_opt_vars
!   CHEyg    只为优化mb_sxyz (在子程序中被sxyz指向)，它包含9个数组，
!   连续访问9个数组相同下标元素，但地址完全不连续，如下：
!     kx = sxyz(1)%r3d(i,j,k)  ...
!     cz = sxyz(9)%r3d(i,j,k)
!   因为 mb_sxyz 在不同次迭代中不变，故只要在网格导数计算之后将其存储到
!   另一个数组 r4d (9,I,J,K) 中，在子程序中使用 r4d 即可。另外因为每个
!   进程算多块，所以要存多块的 r4d， 多块指针为 sxyz1。

    use mod_kndconsts, only : kind_int, kind_real
    implicit none
    
    !---------------------------------------------------------------
    !CHEyg added
    type dim4_fld_array_t
        ! integer(kind_int) :: fldtype
        ! integer(kind_int), pointer :: i4d(:,:,:)
        
        integer(kind_int)        :: nb      ! 全局的块号
        real(kind_real), pointer :: r4d(:,:,:,:)
    end type dim4_fld_array_t

    type(dim4_fld_array_t), pointer :: sxyz1(:)

end module mod_opt_vars

subroutine save_mod_opt_vars()
!   将本进程各块的 mb_sxyz 缓存在 sxyz1 中
 
 
#ifdef PARALLEL
    use mod_variables, only : nghnode    !虚网格层数
	use mod_parallels, only : myid
#endif	
	
	
    use mod_kndconsts, only : kind_int, kind_real
    use mod_fieldvars, only : nblkcoms, blkcoms   ! CHEyg 本地块数和块指针
    use mod_fieldvars, only : mb_top,mb_sxyz    
    use mod_opt_vars
    implicit none
    
    integer(kind_int)  :: nerror, i,j,k,m,nb,i1,i2,j1,j2,k1,k2
    integer(kind_int)  :: nb1, nc, vgrid_ser 
    integer(kind_int)  :: st(3),ed(3)  

    INTEGER vec_shape(3)
    
    allocate(sxyz1( nblkcoms ), stat=nerror)

!    write(*,'(a6,I4,a12,I4)') "myid=",myid,"nblkcoms=",nblkcoms  
	   
	vgrid_ser = 1
	
    do nc=1, nblkcoms
    
       nb1   = blkcoms(nc)%nb        !全局块编号
       st(:) = mb_top(nb1)%ndst(:)
       ed(:) = mb_top(nb1)%nded(:)


       i1 = st(1)
       i2 = ed(1)
       
       j1 = st(2)
       j2 = ed(2)  
               
       k1 = st(3) 
       k2 = ed(3)  


#ifdef PARALLEL
!CHEyg  虚网格最多6层
      i1 = i1 - nghnode
      i2 = i2 + nghnode
	  
      j1 = j1 - nghnode
      j2 = j2 + nghnode
	  
      k1 = k1 - nghnode
      k2 = k2 + nghnode	  
#else
      i1 = i1 - vgrid_ser
      i2 = i2 + vgrid_ser      

      j1 = j1 - vgrid_ser
      j2 = j2 + vgrid_ser      

      k1 = k1 - vgrid_ser
      k2 = k2 + vgrid_ser                 
      
#endif

       allocate(sxyz1(nc)%r4d(9, i1:i2, j1:j2, k1:k2), stat=nerror)
       
       if(nerror .NE. 0) then
          write(*,*) "Allocating memory for sxyz1 wrong. Stop. !"
          stop         
       endif
	   
!       write(*,'(a6,I4,a10,I4)') "myid=",myid,"nghnode=",nghnode  

!      将 mb_sxyz(nb1)%fld的值暂存在 sxyz1(nc)%r4d中

!       vec_shape = SHAPE( mb_sxyz(nb1)%fld(m)%r3d )
!       write(*,*) "myid=",myid,"nghnode=",nghnode,"(i2-i1+1),(j2-j1+1),(k2-k1+1)=", (i2-i1+1),(j2-j1+1),(k2-k1+1),"shape=",vec_shape
	   
	   
       do k = k1,k2
	   ! write(*,*)"myid=",myid,"k=",k
       do j = j1,j2
       do i = i1,i2 
          do m = 1,9
             sxyz1(nc)%r4d(m,i,j,k) = mb_sxyz(nb1)%fld(m)%r3d(i,j,k)
          enddo
       
       enddo
       enddo
       enddo

       sxyz1(nc)%nb = nb1  ! 记录全局的块号，暂无用
       
    enddo
    
    ! write(*,'(a48,6I4)') "save_mod_opt_vars finished. i1,i2,j1,j2,k1,k2=",i1,i2,j1,j2,k1,k2


end subroutine save_mod_opt_vars  
