! ==============================================================================
! Program: Hamiltonian LQCD
! Description: Hamiltonian formulation of Lattice QCD with real-time dynamics
! ==============================================================================
program hamiltonian_lqcd
    use parameters_mod
    use lattice_mod
    use su_n_mod
    use wilson_fermion_mod
    use staggered_fermion_mod
    use hamiltonian_mod
    use time_evolution_mod
    use openmp_mod
    use mpi_mod
    implicit none
    
    character(len=256) :: input_file, output_file
    integer :: nargs
    
    ! Print banner
    call print_banner()
    
    ! Initialize MPI (if enabled)
    call initialize_mpi()
    
    ! Get command line arguments
    nargs = command_argument_count()
    if (nargs >= 1) then
        call get_command_argument(1, input_file)
    else
        input_file = 'input.dat'
    end if
    
    ! Read parameters from input file
    if (is_master()) then
        print *, 'Reading parameters from: ', trim(input_file)
    end if
    call read_parameters(input_file)
    
    ! Broadcast parameters to all MPI processes
    call mpi_broadcast_parameters()
    
    ! Initialize OpenMP (if enabled)
    call initialize_openmp()
    
    ! Initialize lattice structure
    if (is_master()) then
        print *, ''
        print *, 'Initializing lattice structure...'
    end if
    call initialize_lattice()
    
    ! Initialize gauge fields
    if (is_master()) then
        print *, ''
        print *, 'Initializing gauge fields...'
    end if
    call initialize_gauge_fields()
    
    ! Initialize fermion fields based on type
    if (is_master()) then
        print *, ''
        print *, 'Initializing fermion fields...'
    end if
    select case (trim(fermion_type))
        case ('wilson')
            call initialize_wilson_fermions()
        case ('staggered')
            call initialize_staggered_fermions()
        case default
            if (is_master()) then
                print *, 'Warning: Unknown fermion type. Using Wilson fermions.'
            end if
            call initialize_wilson_fermions()
    end select
    
    ! Initialize time evolution
    if (is_master()) then
        print *, ''
        print *, 'Initializing time evolution...'
    end if
    call initialize_time_evolution()
    
    ! Compute initial Hamiltonian
    if (is_master()) then
        print *, ''
        print *, 'Computing initial Hamiltonian...'
    end if
    call compute_hamiltonian()
    if (is_master()) then
        call print_hamiltonian_info()
    end if
    
    ! Perform time evolution (gluon real-time dynamics)
    if (is_master()) then
        print *, ''
        print *, '================================================'
        print *, 'Starting Gluon Real-Time Dynamics Simulation'
        print *, '================================================'
    end if
    call mpi_barrier_all()
    call evolve_gauge_fields()
    
    ! Save results
    if (is_master()) then
        if (nargs >= 2) then
            call get_command_argument(2, output_file)
        else
            output_file = 'time_evolution.dat'
        end if
        
        print *, ''
        print *, 'Saving results...'
        call save_results(output_file)
    end if
    
    ! Print computational method info
    if (is_master()) then
        print *, ''
        print *, '============================================='
        print *, 'Computational method: ', trim(comp_method)
        select case (trim(comp_method))
            case ('tensor_network')
                print *, 'Note: Tensor network methods can be used for'
                print *, '      efficient representation of quantum states'
            case ('quantum')
                print *, 'Note: Quantum computer implementation would use'
                print *, '      quantum gates for time evolution'
            case ('standard')
                print *, 'Using standard classical simulation'
        end select
        print *, '============================================='
    end if
    
    ! Cleanup
    if (is_master()) then
        print *, ''
        print *, 'Cleaning up...'
    end if
    call cleanup_time_evolution()
    
    select case (trim(fermion_type))
        case ('wilson')
            call cleanup_wilson_fermions()
        case ('staggered')
            call cleanup_staggered_fermions()
    end select
    
    call cleanup_gauge_fields()
    call cleanup_lattice()
    call cleanup_parameters()
    
    ! Finalize MPI
    call finalize_mpi()
    
    if (is_master()) then
        print *, ''
        print *, '================================================'
        print *, 'Hamiltonian LQCD simulation completed!'
        print *, '================================================'
        print *, ''
    end if

contains

    ! ==========================================================================
    ! Subroutine: print_banner
    ! Description: Print program banner
    ! ==========================================================================
    subroutine print_banner()
        implicit none
        
        print *, ''
        print *, '================================================'
        print *, '    Hamiltonian Lattice QCD Simulation'
        print *, '    Real-Time Gluon Dynamics'
        print *, '================================================'
        print *, ''
        
    end subroutine print_banner

end program hamiltonian_lqcd
