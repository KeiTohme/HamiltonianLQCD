! ==============================================================================
! Module: parameters_mod
! Description: Parameter management for Hamiltonian LQCD
! ==============================================================================
module parameters_mod
    implicit none
    
    ! Precision parameters
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer, parameter :: sp = selected_real_kind(6, 37)
    
    ! Mathematical constants
    real(dp), parameter :: pi = 3.141592653589793238462643383279502884197_dp
    
    ! Physical parameters
    real(dp) :: gauge_coupling       ! Gauge coupling constant g
    real(dp) :: fermion_mass         ! Fermion mass m
    integer :: num_colors            ! Number of colors Nc (typically 3 for QCD)
    
    ! Lattice parameters
    integer :: spacetime_dim         ! Space-time dimension D_euclid
    integer, allocatable :: lattice_size(:)  ! Lattice size Nsize[D_euclid]
    integer :: total_sites           ! Total number of lattice sites
    character(len=20) :: lattice_type ! 'square' or 'hexagonal'
    
    ! Fermion parameters
    character(len=20) :: fermion_type ! 'wilson' or 'staggered'
    
    ! Parallel parameters
    logical :: use_openmp            ! Use OpenMP parallelization
    logical :: use_mpi               ! Use MPI parallelization
    integer :: num_threads           ! Number of OpenMP threads
    
    ! Time evolution parameters
    real(dp) :: time_step            ! Time step dt
    integer :: num_time_steps        ! Number of time steps
    real(dp) :: total_time           ! Total evolution time
    
    ! Computational method
    character(len=20) :: comp_method ! 'standard', 'tensor_network', 'quantum'
    
contains

    ! ==========================================================================
    ! Subroutine: read_parameters
    ! Description: Read parameters from input file
    ! ==========================================================================
    subroutine read_parameters(filename)
        implicit none
        character(len=*), intent(in) :: filename
        integer :: ios, unit_num
        character(len=100) :: line, key, value
        integer :: eq_pos, i
        
        unit_num = 10
        
        ! Default values
        gauge_coupling = 1.0_dp
        fermion_mass = 0.1_dp
        num_colors = 3
        spacetime_dim = 4
        lattice_type = 'square'
        fermion_type = 'wilson'
        use_openmp = .false.
        use_mpi = .false.
        num_threads = 1
        time_step = 0.01_dp
        num_time_steps = 100
        comp_method = 'standard'
        
        ! Open input file
        open(unit=unit_num, file=filename, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            print *, 'Warning: Could not open parameter file: ', trim(filename)
            print *, 'Using default parameters'
            call set_default_lattice_size()
            return
        end if
        
        ! Read parameters line by line
        do
            read(unit_num, '(A)', iostat=ios) line
            if (ios /= 0) exit
            
            ! Skip empty lines and comments
            line = adjustl(line)
            if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
            
            ! Find '=' position
            eq_pos = index(line, '=')
            if (eq_pos == 0) cycle
            
            key = adjustl(line(1:eq_pos-1))
            value = adjustl(line(eq_pos+1:))
            
            ! Remove inline comments from value
            eq_pos = index(value, '#')
            if (eq_pos > 0) then
                value = value(1:eq_pos-1)
            end if
            value = adjustl(value)
            
            ! Parse parameters
            select case (trim(key))
                case ('gauge_coupling', 'g')
                    read(value, *) gauge_coupling
                case ('fermion_mass', 'm')
                    read(value, *) fermion_mass
                case ('num_colors', 'Nc')
                    read(value, *) num_colors
                case ('spacetime_dim', 'D_euclid')
                    read(value, *) spacetime_dim
                    if (allocated(lattice_size)) deallocate(lattice_size)
                    allocate(lattice_size(spacetime_dim))
                case ('lattice_size', 'Nsize')
                    call parse_lattice_size(value)
                case ('lattice_type')
                    read(value, '(A)') lattice_type
                case ('fermion_type')
                    read(value, '(A)') fermion_type
                case ('use_openmp')
                    read(value, *) use_openmp
                case ('use_mpi')
                    read(value, *) use_mpi
                case ('num_threads')
                    read(value, *) num_threads
                case ('time_step', 'dt')
                    read(value, *) time_step
                case ('num_time_steps')
                    read(value, *) num_time_steps
                case ('comp_method')
                    read(value, '(A)') comp_method
            end select
        end do
        
        close(unit_num)
        
        ! Calculate derived parameters
        total_time = time_step * num_time_steps
        call calculate_total_sites()
        
        ! Print parameters
        call print_parameters()
        
    end subroutine read_parameters
    
    ! ==========================================================================
    ! Subroutine: parse_lattice_size
    ! Description: Parse lattice size array from string
    ! ==========================================================================
    subroutine parse_lattice_size(str)
        implicit none
        character(len=*), intent(in) :: str
        character(len=100) :: temp_str
        integer :: i, start, finish, count
        
        if (.not. allocated(lattice_size)) then
            allocate(lattice_size(spacetime_dim))
        end if
        
        temp_str = adjustl(str)
        count = 1
        start = 1
        
        do i = 1, len_trim(temp_str)
            if (temp_str(i:i) == ',' .or. i == len_trim(temp_str)) then
                if (i == len_trim(temp_str) .and. temp_str(i:i) /= ',') then
                    finish = i
                else
                    finish = i - 1
                end if
                read(temp_str(start:finish), *) lattice_size(count)
                count = count + 1
                start = i + 1
                if (count > spacetime_dim) exit
            end if
        end do
        
    end subroutine parse_lattice_size
    
    ! ==========================================================================
    ! Subroutine: set_default_lattice_size
    ! Description: Set default lattice size
    ! ==========================================================================
    subroutine set_default_lattice_size()
        implicit none
        integer :: i
        
        if (.not. allocated(lattice_size)) then
            allocate(lattice_size(spacetime_dim))
        end if
        
        do i = 1, spacetime_dim
            lattice_size(i) = 8
        end do
        
        call calculate_total_sites()
        
    end subroutine set_default_lattice_size
    
    ! ==========================================================================
    ! Subroutine: calculate_total_sites
    ! Description: Calculate total number of lattice sites
    ! ==========================================================================
    subroutine calculate_total_sites()
        implicit none
        integer :: i
        
        if (.not. allocated(lattice_size)) then
            call set_default_lattice_size()
            return
        end if
        
        total_sites = 1
        do i = 1, spacetime_dim
            total_sites = total_sites * lattice_size(i)
        end do
        
    end subroutine calculate_total_sites
    
    ! ==========================================================================
    ! Subroutine: print_parameters
    ! Description: Print all parameters
    ! ==========================================================================
    subroutine print_parameters()
        implicit none
        integer :: i
        
        print *, '============================================='
        print *, 'Hamiltonian LQCD Parameters'
        print *, '============================================='
        print *, 'Physical parameters:'
        print *, '  Gauge coupling g     = ', gauge_coupling
        print *, '  Fermion mass m       = ', fermion_mass
        print *, '  Number of colors Nc  = ', num_colors
        print *, ''
        print *, 'Lattice parameters:'
        print *, '  Lattice type         = ', trim(lattice_type)
        print *, '  Space-time dimension = ', spacetime_dim
        if (allocated(lattice_size)) then
            print *, '  Lattice size         = ', (lattice_size(i), i=1,spacetime_dim)
        end if
        print *, '  Total sites          = ', total_sites
        print *, ''
        print *, 'Fermion parameters:'
        print *, '  Fermion type         = ', trim(fermion_type)
        print *, ''
        print *, 'Time evolution:'
        print *, '  Time step dt         = ', time_step
        print *, '  Number of steps      = ', num_time_steps
        print *, '  Total time           = ', total_time
        print *, ''
        print *, 'Parallel computing:'
        print *, '  Use OpenMP           = ', use_openmp
        print *, '  Use MPI              = ', use_mpi
        print *, '  Number of threads    = ', num_threads
        print *, ''
        print *, 'Computational method   = ', trim(comp_method)
        print *, '============================================='
        print *, ''
        
    end subroutine print_parameters
    
    ! ==========================================================================
    ! Subroutine: cleanup_parameters
    ! Description: Cleanup allocated memory
    ! ==========================================================================
    subroutine cleanup_parameters()
        implicit none
        
        if (allocated(lattice_size)) deallocate(lattice_size)
        
    end subroutine cleanup_parameters

end module parameters_mod
