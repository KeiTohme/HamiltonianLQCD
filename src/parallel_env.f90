module parallel_env
  use kinds
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private
  integer :: world_rank = 0
  integer :: world_size = 1
  logical :: mpi_enabled = .false.
  logical :: openmp_enabled = .false.

  public :: init_parallel, finalize_parallel, get_rank, get_size, barrier
  public :: is_mpi_enabled, is_openmp_enabled

contains

  subroutine init_parallel()
    implicit none
#ifdef USE_MPI
    integer :: ierr, provided
    call MPI_Initialized(mpi_enabled, ierr)
    if (.not. mpi_enabled) then
      call MPI_Init_thread(MPI_THREAD_FUNNELED, provided, ierr)
      mpi_enabled = .true.
    end if
    call MPI_Comm_rank(MPI_COMM_WORLD, world_rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, world_size, ierr)
#else
    mpi_enabled = .false.
    world_rank = 0
    world_size = 1
#endif
#ifdef USE_OPENMP
    openmp_enabled = .true.
#else
    openmp_enabled = .false.
#endif
  end subroutine init_parallel

  subroutine finalize_parallel()
    implicit none
#ifdef USE_MPI
    integer :: ierr, finalized
    call MPI_Finalized(finalized, ierr)
    if (.not. finalized) then
      call MPI_Finalize(ierr)
    end if
#endif
  end subroutine finalize_parallel

  integer function get_rank()
    implicit none
    get_rank = world_rank
  end function get_rank

  integer function get_size()
    implicit none
    get_size = world_size
  end function get_size

  subroutine barrier()
    implicit none
#ifdef USE_MPI
    integer :: ierr
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
  end subroutine barrier

  logical function is_mpi_enabled()
    implicit none
    is_mpi_enabled = mpi_enabled
  end function is_mpi_enabled

  logical function is_openmp_enabled()
    implicit none
    is_openmp_enabled = openmp_enabled
  end function is_openmp_enabled

end module parallel_env
