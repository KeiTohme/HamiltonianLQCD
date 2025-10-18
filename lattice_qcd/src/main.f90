program lattice_qcd_hamiltonian
    use parameter_module
    use lattice_module
    use gauge_field_module
    use fermion_module
    use hamiltonian_module
    use time_evolution_module
    use parallel_module
    implicit none
    
    ! Main program variables
    type(parameters) :: params
    type(lattice) :: lat
    type(gauge_field) :: gauge
    type(fermion_field) :: fermions
    type(hamiltonian) :: H
    character(len=256) :: input_file
    logical :: use_parallel
    integer :: ierr
    
    ! Initialize MPI if requested
    call get_command_argument(1, input_file)
    if (len_trim(input_file) == 0) then
        print *, "Usage: lattice_qcd_hamiltonian <input_parameter_file>"
        stop
    endif
    
    ! Read parameters
    call read_parameters(input_file, params)
    
    ! Initialize parallel environment if requested
    use_parallel = params%use_mpi .or. params%use_openmp
    if (use_parallel) then
        call initialize_parallel(params%use_mpi, params%use_openmp, ierr)
        if (ierr /= 0) then
            print *, "Error initializing parallel environment"
            stop
        endif
    endif
    
    ! Initialize lattice
    call initialize_lattice(lat, params)
    
    ! Initialize gauge fields
    call initialize_gauge_field(gauge, lat, params)
    
    ! Initialize fermion fields
    call initialize_fermion_field(fermions, lat, params)
    
    ! Construct Hamiltonian
    call construct_hamiltonian(H, gauge, fermions, lat, params)
    
    ! Perform time evolution
    call time_evolution_main(H, gauge, fermions, lat, params)
    
    ! Cleanup
    call cleanup_hamiltonian(H)
    call cleanup_fermion_field(fermions)
    call cleanup_gauge_field(gauge)
    call cleanup_lattice(lat)
    
    if (use_parallel) then
        call finalize_parallel(params%use_mpi)
    endif
    
end program lattice_qcd_hamiltonian