module parameter_module
    implicit none
    
    type parameters
        ! Physical parameters
        real(8) :: g                    ! Gauge coupling constant
        real(8) :: m_fermion           ! Fermion mass
        integer :: Nc                   ! Number of colors
        integer :: D_euclid            ! Euclidean spacetime dimension
        
        ! Lattice parameters
        integer, allocatable :: Nsize(:)  ! Lattice size in each dimension
        character(len=20) :: lattice_type ! "square" or "hexagonal"
        
        ! Fermion parameters
        character(len=20) :: fermion_type ! "wilson" or "staggered"
        
        ! Simulation parameters
        real(8) :: dt                   ! Time step
        integer :: n_steps              ! Number of time steps
        integer :: measure_interval     ! Measurement interval
        
        ! Parallelization options
        logical :: use_mpi
        logical :: use_openmp
        integer :: n_threads           ! For OpenMP
        
        ! Advanced options
        logical :: use_tensor_network
        logical :: use_quantum_computer
        character(len=256) :: tn_backend    ! Tensor network backend
        character(len=256) :: qc_backend    ! Quantum computer backend
        
        ! Output options
        character(len=256) :: output_dir
        logical :: save_configurations
        integer :: save_interval
    end type parameters
    
contains
    
    subroutine read_parameters(filename, params)
        character(len=*), intent(in) :: filename
        type(parameters), intent(out) :: params
        integer :: unit_num, iostat, i
        integer :: D_euclid_temp
        character(len=256) :: line, keyword, value
        
        unit_num = 10
        open(unit=unit_num, file=filename, status='old', action='read', iostat=iostat)
        if (iostat /= 0) then
            print *, "Error: Cannot open parameter file: ", trim(filename)
            stop
        endif
        
        ! Set default values
        params%g = 1.0d0
        params%m_fermion = 0.1d0
        params%Nc = 3
        params%D_euclid = 4
        params%lattice_type = "square"
        params%fermion_type = "wilson"
        params%dt = 0.01d0
        params%n_steps = 1000
        params%measure_interval = 10
        params%use_mpi = .false.
        params%use_openmp = .false.
        params%n_threads = 1
        params%use_tensor_network = .false.
        params%use_quantum_computer = .false.
        params%tn_backend = "none"
        params%qc_backend = "none"
        params%output_dir = "./output"
        params%save_configurations = .true.
        params%save_interval = 100
        
        ! Read parameters line by line
        do
            read(unit_num, '(A)', iostat=iostat) line
            if (iostat /= 0) exit
            
            ! Skip comments and empty lines
            line = adjustl(line)
            if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
            
            ! Parse keyword = value
            i = index(line, '=')
            if (i == 0) cycle
            
            keyword = adjustl(line(1:i-1))
            value = adjustl(line(i+1:))
            
            select case(trim(keyword))
                case('g', 'gauge_coupling')
                    read(value, *) params%g
                case('m', 'm_fermion', 'fermion_mass')
                    read(value, *) params%m_fermion
                case('Nc', 'n_colors')
                    read(value, *) params%Nc
                case('D', 'D_euclid', 'dimension')
                    read(value, *) D_euclid_temp
                    params%D_euclid = D_euclid_temp
                    allocate(params%Nsize(D_euclid_temp))
                case('lattice_type')
                    params%lattice_type = trim(value)
                case('fermion_type')
                    params%fermion_type = trim(value)
                case('dt', 'time_step')
                    read(value, *) params%dt
                case('n_steps')
                    read(value, *) params%n_steps
                case('measure_interval')
                    read(value, *) params%measure_interval
                case('use_mpi')
                    read(value, *) params%use_mpi
                case('use_openmp')
                    read(value, *) params%use_openmp
                case('n_threads')
                    read(value, *) params%n_threads
                case('use_tensor_network')
                    read(value, *) params%use_tensor_network
                case('use_quantum_computer')
                    read(value, *) params%use_quantum_computer
                case('output_dir')
                    params%output_dir = trim(value)
            end select
            
            ! Special handling for Nsize array
            if (index(keyword, 'Nsize') > 0) then
                if (allocated(params%Nsize)) then
                    read(value, *) params%Nsize
                endif
            endif
        end do
        
        close(unit_num)
        
        ! Validate parameters
        call validate_parameters(params)
        
    end subroutine read_parameters
    
    subroutine validate_parameters(params)
        type(parameters), intent(in) :: params
        
        if (params%g <= 0.0d0) then
            print *, "Error: Gauge coupling must be positive"
            stop
        endif
        
        if (params%Nc < 2) then
            print *, "Error: Number of colors must be at least 2"
            stop
        endif
        
        if (params%D_euclid < 2 .or. params%D_euclid > 4) then
            print *, "Error: Dimension must be between 2 and 4"
            stop
        endif
        
        if (trim(params%lattice_type) /= "square" .and. &
            trim(params%lattice_type) /= "hexagonal") then
            print *, "Error: Lattice type must be 'square' or 'hexagonal'"
            stop
        endif
        
        if (trim(params%fermion_type) /= "wilson" .and. &
            trim(params%fermion_type) /= "staggered") then
            print *, "Error: Fermion type must be 'wilson' or 'staggered'"
            stop
        endif
        
    end subroutine validate_parameters
    
    subroutine print_parameters(params)
        type(parameters), intent(in) :: params
        integer :: i
        
        print *, "===== Lattice QCD Parameters ====="
        print *, "Physical parameters:"
        print *, "  Gauge coupling g =", params%g
        print *, "  Fermion mass m =", params%m_fermion
        print *, "  Number of colors Nc =", params%Nc
        print *, "  Dimension D =", params%D_euclid
        print *, "Lattice parameters:"
        print *, "  Lattice type:", trim(params%lattice_type)
        print *, "  Lattice size:", (params%Nsize(i), i=1,params%D_euclid)
        print *, "  Fermion type:", trim(params%fermion_type)
        print *, "Simulation parameters:"
        print *, "  Time step dt =", params%dt
        print *, "  Number of steps =", params%n_steps
        print *, "Parallelization:"
        print *, "  Use MPI:", params%use_mpi
        print *, "  Use OpenMP:", params%use_openmp
        if (params%use_openmp) print *, "  Number of threads:", params%n_threads
        print *, "Advanced options:"
        print *, "  Use tensor network:", params%use_tensor_network
        print *, "  Use quantum computer:", params%use_quantum_computer
        print *, "=================================="
        
    end subroutine print_parameters
    
end module parameter_module