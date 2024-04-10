
module mod_parallels
#ifdef PARALLEL
    use mod_kndconsts, only : kind_int,kind_real
    use mod_datatypes, only : bc_region_t
    use mpi
    implicit none

    integer(kind_int), parameter :: master = 0
    integer(kind_int) :: kind_int_mpi,kind_real_mpi,kind_char_mpi

    integer(kind_int) :: myid,numprocs
    character(len=MPI_MAX_PROCESSOR_NAME) :: procname
    integer(kind_int) :: len_proc_name

    character(len=MPI_MAX_PROCESSOR_NAME), pointer :: procnames(:)

    integer(kind_int)           :: nmsg_sin
    integer(kind_int), pointer  :: request_sin(:), status_sin(:, :)

    integer(kind_int)           :: nmsg_int
    integer(kind_int), pointer  :: request_int(:), status_int(:, :)

#endif

end module mod_parallels

module mod_openmp
#ifdef OMP_IMP
    use mod_kndconsts, only : kind_int,kind_real
    use omp_lib
    implicit none

    integer(kind_int), parameter :: omp_max_num_threads = 128
#endif
end module mod_openmp
