module parallel_module
    implicit none
    
    ! MPI related variables
    integer :: mpi_rank, mpi_size
    integer :: mpi_comm_world
    logical :: mpi_initialized = .false.
    
    ! OpenMP related variables
    integer :: omp_num_threads
    logical :: omp_initialized = .false.
    
    ! Domain decomposition
    integer, allocatable :: local_sites(:)
    integer :: n_local_sites
    integer :: site_offset
    
contains
    
    subroutine initialize_parallel(use_mpi, use_openmp, ierr)
        logical, intent(in) :: use_mpi, use_openmp
        integer, intent(out) :: ierr
        
        ierr = 0
        
        if (use_mpi) then
            call initialize_mpi(ierr)
            if (ierr /= 0) return
        else
            mpi_rank = 0
            mpi_size = 1
        endif
        
        if (use_openmp) then
            call initialize_openmp()
        endif
        
    end subroutine initialize_parallel
    
    subroutine initialize_mpi(ierr)
        integer, intent(out) :: ierr
        
#ifdef USE_MPI
        include 'mpif.h'
        
        call MPI_Init(ierr)
        if (ierr /= 0) then
            print *, "Error initializing MPI"
            return
        endif
        
        mpi_comm_world = MPI_COMM_WORLD
        call MPI_Comm_rank(mpi_comm_world, mpi_rank, ierr)
        call MPI_Comm_size(mpi_comm_world, mpi_size, ierr)
        
        mpi_initialized = .true.
        
        if (mpi_rank == 0) then
            print *, "MPI initialized with", mpi_size, "processes"
        endif
#else
        ierr = 0
        mpi_rank = 0
        mpi_size = 1
        if (mpi_rank == 0) then
            print *, "MPI support not compiled in"
        endif
#endif
        
    end subroutine initialize_mpi
    
    subroutine initialize_openmp()
        
#ifdef _OPENMP
        use omp_lib
        
        omp_num_threads = omp_get_max_threads()
        omp_initialized = .true.
        
        if (mpi_rank == 0) then
            print *, "OpenMP initialized with", omp_num_threads, "threads"
        endif
#else
        omp_num_threads = 1
        if (mpi_rank == 0) then
            print *, "OpenMP support not compiled in"
        endif
#endif
        
    end subroutine initialize_openmp
    
    subroutine setup_domain_decomposition(n_sites, n_local, local_list, offset)
        integer, intent(in) :: n_sites
        integer, intent(out) :: n_local
        integer, allocatable, intent(out) :: local_list(:)
        integer, intent(out) :: offset
        
        integer :: sites_per_proc, remainder
        integer :: i, start_site, end_site
        
        ! Simple 1D decomposition
        sites_per_proc = n_sites / mpi_size
        remainder = mod(n_sites, mpi_size)
        
        ! Calculate local number of sites
        if (mpi_rank < remainder) then
            n_local = sites_per_proc + 1
            offset = mpi_rank * (sites_per_proc + 1)
        else
            n_local = sites_per_proc
            offset = remainder * (sites_per_proc + 1) + &
                    (mpi_rank - remainder) * sites_per_proc
        endif
        
        ! Create list of local sites
        allocate(local_list(n_local))
        start_site = offset + 1
        end_site = offset + n_local
        
        do i = 1, n_local
            local_list(i) = start_site + i - 1
        end do
        
    end subroutine setup_domain_decomposition
    
    subroutine parallel_sum_real(local_val, global_val)
        real(8), intent(in) :: local_val
        real(8), intent(out) :: global_val
        
#ifdef USE_MPI
        include 'mpif.h'
        integer :: ierr
        
        if (mpi_initialized) then
            call MPI_Allreduce(local_val, global_val, 1, MPI_DOUBLE_PRECISION, &
                              MPI_SUM, mpi_comm_world, ierr)
        else
            global_val = local_val
        endif
#else
        global_val = local_val
#endif
        
    end subroutine parallel_sum_real
    
    subroutine parallel_sum_complex(local_val, global_val)
        complex(8), intent(in) :: local_val
        complex(8), intent(out) :: global_val
        
#ifdef USE_MPI
        include 'mpif.h'
        integer :: ierr
        
        if (mpi_initialized) then
            call MPI_Allreduce(local_val, global_val, 1, MPI_DOUBLE_COMPLEX, &
                              MPI_SUM, mpi_comm_world, ierr)
        else
            global_val = local_val
        endif
#else
        global_val = local_val
#endif
        
    end subroutine parallel_sum_complex
    
    subroutine parallel_max_real(local_val, global_val)
        real(8), intent(in) :: local_val
        real(8), intent(out) :: global_val
        
#ifdef USE_MPI
        include 'mpif.h'
        integer :: ierr
        
        if (mpi_initialized) then
            call MPI_Allreduce(local_val, global_val, 1, MPI_DOUBLE_PRECISION, &
                              MPI_MAX, mpi_comm_world, ierr)
        else
            global_val = local_val
        endif
#else
        global_val = local_val
#endif
        
    end subroutine parallel_max_real
    
    subroutine sync_gauge_field(gauge, lat, params)
        use gauge_field_module
        use lattice_module
        use parameter_module
        type(gauge_field), intent(inout) :: gauge
        type(lattice), intent(in) :: lat
        type(parameters), intent(in) :: params
        
#ifdef USE_MPI
        include 'mpif.h'
        integer :: ierr
        integer :: buffer_size
        
        if (mpi_initialized .and. mpi_size > 1) then
            ! Synchronize gauge field across processes
            buffer_size = params%Nc * params%Nc * lat%n_sites * gauge%n_links
            
            ! Broadcast gauge links
            call MPI_Bcast(gauge%U, buffer_size, MPI_DOUBLE_COMPLEX, 0, &
                          mpi_comm_world, ierr)
            
            ! Broadcast electric field
            call MPI_Bcast(gauge%E, buffer_size, MPI_DOUBLE_COMPLEX, 0, &
                          mpi_comm_world, ierr)
        endif
#endif
        
    end subroutine sync_gauge_field
    
    subroutine sync_fermion_field(fermions, lat, params)
        use fermion_module
        use lattice_module
        use parameter_module
        type(fermion_field), intent(inout) :: fermions
        type(lattice), intent(in) :: lat
        type(parameters), intent(in) :: params
        
#ifdef USE_MPI
        include 'mpif.h'
        integer :: ierr
        integer :: buffer_size
        
        if (mpi_initialized .and. mpi_size > 1) then
            ! Synchronize fermion fields across processes
            buffer_size = params%Nc * fermions%n_spin * lat%n_sites * fermions%n_flavors
            
            ! Broadcast fermion fields
            call MPI_Bcast(fermions%psi, buffer_size, MPI_DOUBLE_COMPLEX, 0, &
                          mpi_comm_world, ierr)
            
            call MPI_Bcast(fermions%pi_psi, buffer_size, MPI_DOUBLE_COMPLEX, 0, &
                          mpi_comm_world, ierr)
        endif
#endif
        
    end subroutine sync_fermion_field
    
    subroutine parallel_wilson_dirac_operator(psi_out, psi_in, gauge, fermions, lat, params, dag)
        use fermion_module
        use gauge_field_module
        use lattice_module
        use parameter_module
        
        complex(8), intent(out) :: psi_out(:,:,:,:)
        complex(8), intent(in) :: psi_in(:,:,:,:)
        type(gauge_field), intent(in) :: gauge
        type(fermion_field), intent(in) :: fermions
        type(lattice), intent(in) :: lat
        type(parameters), intent(in) :: params
        logical, intent(in) :: dag
        
        integer :: site, local_site, color, spin, flavor
        
        ! Initialize output
        psi_out = (0.0d0, 0.0d0)
        
#ifdef _OPENMP
        !$omp parallel do private(local_site, site, color, spin, flavor) &
        !$omp& schedule(static)
#endif
        do local_site = 1, n_local_sites
            site = local_sites(local_site)
            
            ! Apply Wilson operator to this site
            ! (Implementation would go here - simplified for brevity)
            
        end do
#ifdef _OPENMP
        !$omp end parallel do
#endif
        
        ! Communication for boundary sites if using MPI
        if (mpi_initialized .and. mpi_size > 1) then
            call communicate_fermion_boundaries(psi_out, lat, params)
        endif
        
    end subroutine parallel_wilson_dirac_operator
    
    subroutine communicate_fermion_boundaries(psi, lat, params)
        use lattice_module
        use parameter_module
        
        complex(8), intent(inout) :: psi(:,:,:,:)
        type(lattice), intent(in) :: lat
        type(parameters), intent(in) :: params
        
#ifdef USE_MPI
        include 'mpif.h'
        integer :: ierr
        
        ! Exchange boundary data between neighboring processes
        ! (Implementation would depend on domain decomposition strategy)
#endif
        
    end subroutine communicate_fermion_boundaries
    
    subroutine parallel_plaquette_sum(local_plaq, global_plaq, n_plaq_local, n_plaq_global)
        real(8), intent(in) :: local_plaq
        real(8), intent(out) :: global_plaq
        integer, intent(in) :: n_plaq_local
        integer, intent(out) :: n_plaq_global
        
#ifdef USE_MPI
        include 'mpif.h'
        integer :: ierr
        
        if (mpi_initialized) then
            call MPI_Allreduce(local_plaq, global_plaq, 1, MPI_DOUBLE_PRECISION, &
                              MPI_SUM, mpi_comm_world, ierr)
            call MPI_Allreduce(n_plaq_local, n_plaq_global, 1, MPI_INTEGER, &
                              MPI_SUM, mpi_comm_world, ierr)
        else
            global_plaq = local_plaq
            n_plaq_global = n_plaq_local
        endif
#else
        global_plaq = local_plaq
        n_plaq_global = n_plaq_local
#endif
        
    end subroutine parallel_plaquette_sum
    
    subroutine finalize_parallel(use_mpi)
        logical, intent(in) :: use_mpi
        
#ifdef USE_MPI
        include 'mpif.h'
        integer :: ierr
        
        if (use_mpi .and. mpi_initialized) then
            call MPI_Finalize(ierr)
            mpi_initialized = .false.
        endif
#endif
        
        if (allocated(local_sites)) deallocate(local_sites)
        
    end subroutine finalize_parallel
    
    logical function is_master_process()
        is_master_process = (mpi_rank == 0)
    end function is_master_process
    
    subroutine parallel_barrier()
#ifdef USE_MPI
        include 'mpif.h'
        integer :: ierr
        
        if (mpi_initialized) then
            call MPI_Barrier(mpi_comm_world, ierr)
        endif
#endif
    end subroutine parallel_barrier
    
end module parallel_module