! ==============================================================================
! Module: openmp_mod
! Description: OpenMP parallelization support
! ==============================================================================
module openmp_mod
    use parameters_mod
    implicit none
    
contains

    ! ==========================================================================
    ! Subroutine: initialize_openmp
    ! Description: Initialize OpenMP environment
    ! ==========================================================================
    subroutine initialize_openmp()
        implicit none
        
        if (use_openmp) then
            print *, 'Warning: OpenMP requested but not available in this build'
            print *, 'Compile with OpenMP flag (make openmp) to enable OpenMP'
            use_openmp = .false.
        end if
        
    end subroutine initialize_openmp
    
    ! ==========================================================================
    ! Subroutine: parallel_loop_example
    ! Description: Example of parallelized loop over lattice sites
    ! ==========================================================================
    subroutine parallel_loop_example()
        implicit none
        integer :: site
        real(dp) :: local_sum
        
        local_sum = 0.0_dp
        do site = 1, total_sites
            ! Example computation on each site
            local_sum = local_sum + 1.0_dp
        end do
        
    end subroutine parallel_loop_example

end module openmp_mod
