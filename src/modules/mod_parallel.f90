!***********************************************************************
! Module: mod_parallel
! Purpose: Implement OpenMP and MPI parallelization
!***********************************************************************
module mod_parallel
  use mod_parameters
  implicit none
  
contains

  !=====================================================================
  ! Subroutine: initialize_parallel
  ! Purpose: Initialize parallel environment (OpenMP/MPI)
  !=====================================================================
  subroutine initialize_parallel()
    implicit none
    integer :: ierr
    
    ! Initialize MPI if requested
    if (use_mpi) then
#ifdef USE_MPI
      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, mpi_rank, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_size, ierr)
      
      if (mpi_rank == 0) then
        write(*,'(A,I4,A)') ' MPI initialized with ', mpi_size, ' processes'
      end if
#else
      if (mpi_rank == 0) then
        write(*,*) 'Warning: MPI requested but not compiled with MPI support'
        use_mpi = .false.
      end if
#endif
    end if
    
    ! Initialize OpenMP if requested
    if (use_openmp) then
#ifdef _OPENMP
      call omp_set_num_threads(n_threads)
      write(*,'(A,I4,A)') ' OpenMP initialized with ', n_threads, ' threads'
#else
      write(*,*) 'Warning: OpenMP requested but not compiled with OpenMP support'
      use_openmp = .false.
#endif
    end if
    
    if (.not. use_mpi .and. .not. use_openmp) then
      write(*,'(A)') ' Running in serial mode'
    end if
    
  end subroutine initialize_parallel
  
  !=====================================================================
  ! Subroutine: finalize_parallel
  ! Purpose: Finalize parallel environment
  !=====================================================================
  subroutine finalize_parallel()
    implicit none
    integer :: ierr
    
    if (use_mpi) then
#ifdef USE_MPI
      call MPI_FINALIZE(ierr)
#endif
    end if
    
  end subroutine finalize_parallel
  
  !=====================================================================
  ! Subroutine: parallel_sum_real
  ! Purpose: Perform parallel sum reduction for real array
  !=====================================================================
  subroutine parallel_sum_real(local_val, global_val)
    implicit none
    real(dp), intent(in) :: local_val
    real(dp), intent(out) :: global_val
    integer :: ierr
    
    if (use_mpi) then
#ifdef USE_MPI
      call MPI_ALLREDUCE(local_val, global_val, 1, MPI_DOUBLE_PRECISION, &
                         MPI_SUM, MPI_COMM_WORLD, ierr)
#else
      global_val = local_val
#endif
    else
      global_val = local_val
    end if
    
  end subroutine parallel_sum_real
  
  !=====================================================================
  ! Subroutine: parallel_barrier
  ! Purpose: Synchronization barrier
  !=====================================================================
  subroutine parallel_barrier()
    implicit none
    integer :: ierr
    
    if (use_mpi) then
#ifdef USE_MPI
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif
    end if
    
  end subroutine parallel_barrier

end module mod_parallel
