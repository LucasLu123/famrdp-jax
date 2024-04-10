
module mod_datatypes
    use mod_kndconsts, only : kind_int,kind_real,len_char_name
    implicit none
    
    type fld_array_t
        integer(kind_int) :: fldtype
        
        integer(kind_int), pointer :: i3d(:,:,:)
        real(kind_real)  , pointer :: r3d(:,:,:)
    end type fld_array_t

    type var_block_t
        character(len=len_char_name) :: varname
        
        type(fld_array_t), pointer :: fld(:)
    end type var_block_t
    
    type bc_region_t
        integer(kind_int) :: bctype
        integer(kind_int) :: s_st(3),s_ed(3)
        integer(kind_int) :: nbt
        integer(kind_int) :: t_st(3),t_ed(3)
        
        integer(kind_int) :: nbs
        integer(kind_int) :: nrs
        integer(kind_int) :: s_nd,s_lr
        integer(kind_int) :: nrt
        integer(kind_int) :: t_nd,t_lr
        
        integer(kind_int) :: s_dir(3),s_ord(3),s_sgn(3)
        integer(kind_int) :: t_dir(3),t_ord(3),t_sgn(3)
        
        integer(kind_int) :: mapmat(3,4)
        integer(kind_int), pointer :: mapijk(:,:,:,:)

        integer(kind_int), pointer :: bcts(:,:,:)  !!每个点上的bctype
        
        integer(kind_int) :: subtype

        integer(kind_int) :: subtype_A
        
        integer(kind_int) :: nint

        integer(kind_int) :: mark
    end type bc_region_t
    
    type top_block_t
        integer(kind_int) :: pid

        integer(kind_int) :: pnb
        integer(kind_int) :: pst(3),ped(3)

        integer(kind_int) :: nupd

        integer(kind_int) :: nijk(3)

        integer(kind_int) :: ndst(3)
        integer(kind_int) :: nded(3)

        character(len_char_name) :: name
        
        integer(kind_int) :: nregions
        type(bc_region_t), pointer :: bcs(:)
        
        integer(kind_int), pointer :: bcnrs(:)

        integer(kind_int) :: mark
    end type top_block_t

    type bc_buffer_t
        real(kind_real), pointer :: dat(:,:,:,:)
    end type bc_buffer_t

    type bc_inter_t
        integer(kind_int)          :: nint
        type(bc_region_t), pointer :: bc
        real(kind_real)  , pointer :: dat(:,:,:,:)

        integer(kind_int)          :: nbuf
        type(bc_buffer_t), pointer :: buf(:)
    end type bc_inter_t

    type block_compute_t
        integer(kind_int) :: nb
        
        type(top_block_t), pointer :: top
    end type block_compute_t

    type wll_point_t
       real(kind_real)   :: xw,yw,zw
       integer(kind_int) :: nb
       integer(kind_int) :: i,j,k
    end type wll_point_t

    type patched_t
        real(kind_real),   pointer :: dats(:,:,:,:)
        real(kind_real),   pointer :: datt(:,:,:,:)        
        integer(kind_int), pointer :: id1(:,:,:,:)
        integer(kind_int), pointer :: id2(:,:,:,:)
    end type patched_t

end module mod_datatypes
