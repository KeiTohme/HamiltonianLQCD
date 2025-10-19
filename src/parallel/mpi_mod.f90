! ==============================================================================
! Module: mpi_mod
! Description: MPI parallelization support
! ==============================================================================
module mpi_mod
    use parameters_mod
    implicit none
    
    ! MPI parameters
    integer :: mpi_rank, mpi_size
    integer :: mpi_ierr
    logical :: mpi_initialized = .false.
    
contains

    ! ==========================================================================
    ! Subroutine: initialize_mpi
    ! Description: Initialize MPI environment
    ! ==========================================================================
    subroutine initialize_mpi()
        implicit none
        
        ! MPI support not enabled in this build
        mpi_rank = 0
        mpi_size = 1
        
        if (use_mpi) then
            print *, 'Warning: MPI requested but not available in this build'
            print *, 'Compile with MPI support (make mpi) to enable MPI'
            use_mpi = .false.
        end if
        
    end subroutine initialize_mpi
    
    ! ==========================================================================
    ! Subroutine: finalize_mpi
    ! Description: Finalize MPI environment
    ! ==========================================================================
    subroutine finalize_mpi()
        implicit none
        
        ! Nothing to do when MPI is not enabled
        
    end subroutine finalize_mpi
    
    ! ==========================================================================
    ! Subroutine: mpi_barrier_all
    ! Description: Synchronize all MPI processes
    ! ==========================================================================
    subroutine mpi_barrier_all()
        implicit none
        
        ! Nothing to do when MPI is not enabled
        
    end subroutine mpi_barrier_all
    
    ! ==========================================================================
    ! Subroutine: mpi_broadcast_parameters
    ! Description: Broadcast parameters from rank 0 to all processes
    ! ==========================================================================
    subroutine mpi_broadcast_parameters()
        implicit none
        
        ! Nothing to do when MPI is not enabled
        
    end subroutine mpi_broadcast_parameters
    
    ! ==========================================================================
    ! Function: is_master
    ! Description: Check if current process is master (rank 0)
    ! ==========================================================================
    function is_master() result(master)
        implicit none
        logical :: master
        
        master = (mpi_rank == 0)
        
    end function is_master

end module mpi_mod
