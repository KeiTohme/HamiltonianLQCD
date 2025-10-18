!***********************************************************************
! Module: mod_parameters
! Purpose: Manage all input parameters for Lattice QCD simulation
!***********************************************************************
module mod_parameters
  implicit none
  
  ! Precision
  integer, parameter :: dp = selected_real_kind(15, 307)
  integer, parameter :: sp = selected_real_kind(6, 37)
  
  ! Lattice configuration
  character(len=20) :: lattice_type
  integer :: d_euclid
  integer, allocatable :: nsize(:)
  integer :: total_sites
  
  ! Physical parameters
  real(dp) :: gauge_coupling
  real(dp) :: fermion_mass
  integer :: n_colors
  real(dp) :: beta
  
  ! Fermion configuration
  character(len=20) :: fermion_type
  integer :: n_flavors
  
  ! Time evolution
  real(dp) :: time_step
  integer :: n_steps
  integer :: output_interval
  
  ! Parallelization
  logical :: use_openmp
  logical :: use_mpi
  integer :: n_threads
  
  ! Computational method
  character(len=20) :: computation_method
  integer :: tensor_chi
  character(len=20) :: quantum_backend
  
  ! MPI parameters (if used)
  integer :: mpi_rank = 0
  integer :: mpi_size = 1
  
contains

  !=====================================================================
  ! Subroutine: read_parameters
  ! Purpose: Read parameters from input file
  !=====================================================================
  subroutine read_parameters(filename)
    implicit none
    character(len=*), intent(in) :: filename
    integer :: unit_num, ios, i
    
    namelist /lattice_config/ lattice_type, d_euclid, nsize
    namelist /physical_params/ gauge_coupling, fermion_mass, n_colors, beta
    namelist /fermion_config/ fermion_type, n_flavors
    namelist /time_evolution/ time_step, n_steps, output_interval
    namelist /parallel_config/ use_openmp, use_mpi, n_threads
    namelist /computation_method/ computation_method, tensor_chi, quantum_backend
    
    ! Default values
    lattice_type = 'square'
    d_euclid = 4
    gauge_coupling = 1.0_dp
    fermion_mass = 0.1_dp
    n_colors = 3
    beta = 6.0_dp
    fermion_type = 'wilson'
    n_flavors = 2
    time_step = 0.01_dp
    n_steps = 1000
    output_interval = 10
    use_openmp = .false.
    use_mpi = .false.
    n_threads = 1
    computation_method = 'classical'
    tensor_chi = 50
    quantum_backend = 'qiskit'
    
    ! Allocate nsize with default dimension
    if (.not. allocated(nsize)) allocate(nsize(4))
    nsize = 8
    
    ! Open and read file
    open(newunit=unit_num, file=filename, status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'Warning: Could not open parameter file: ', trim(filename)
      write(*,*) 'Using default parameters'
      return
    end if
    
    read(unit_num, nml=lattice_config, iostat=ios)
    if (ios /= 0) write(*,*) 'Warning: Error reading lattice_config'
    
    ! Reallocate nsize if dimension changed
    if (size(nsize) /= d_euclid) then
      deallocate(nsize)
      allocate(nsize(d_euclid))
      nsize = 8
      rewind(unit_num)
      read(unit_num, nml=lattice_config, iostat=ios)
    end if
    
    rewind(unit_num)
    read(unit_num, nml=physical_params, iostat=ios)
    if (ios /= 0) write(*,*) 'Warning: Error reading physical_params'
    
    rewind(unit_num)
    read(unit_num, nml=fermion_config, iostat=ios)
    if (ios /= 0) write(*,*) 'Warning: Error reading fermion_config'
    
    rewind(unit_num)
    read(unit_num, nml=time_evolution, iostat=ios)
    if (ios /= 0) write(*,*) 'Warning: Error reading time_evolution'
    
    rewind(unit_num)
    read(unit_num, nml=parallel_config, iostat=ios)
    if (ios /= 0) write(*,*) 'Warning: Error reading parallel_config'
    
    rewind(unit_num)
    read(unit_num, nml=computation_method, iostat=ios)
    if (ios /= 0) write(*,*) 'Warning: Error reading computation_method'
    
    close(unit_num)
    
    ! Calculate total sites
    total_sites = 1
    do i = 1, d_euclid
      total_sites = total_sites * nsize(i)
    end do
    
    call print_parameters()
    
  end subroutine read_parameters
  
  !=====================================================================
  ! Subroutine: print_parameters
  ! Purpose: Display current parameters
  !=====================================================================
  subroutine print_parameters()
    implicit none
    integer :: i
    
    write(*,'(A)') '====================================================================='
    write(*,'(A)') '           Lattice QCD Hamiltonian Formalism'
    write(*,'(A)') '====================================================================='
    write(*,'(A,A)') ' Lattice type:        ', trim(lattice_type)
    write(*,'(A,I3)') ' Spacetime dimension: ', d_euclid
    write(*,'(A)', advance='no') ' Lattice size:        '
    do i = 1, d_euclid
      write(*,'(I4)', advance='no') nsize(i)
      if (i < d_euclid) write(*,'(A)', advance='no') ' x'
    end do
    write(*,*)
    write(*,'(A,I10)') ' Total sites:         ', total_sites
    write(*,'(A)') '---------------------------------------------------------------------'
    write(*,'(A,F10.4)') ' Gauge coupling g:    ', gauge_coupling
    write(*,'(A,F10.4)') ' Fermion mass m:      ', fermion_mass
    write(*,'(A,I3)') ' Number of colors:    ', n_colors
    write(*,'(A,F10.4)') ' Beta (2Nc/g²):       ', beta
    write(*,'(A)') '---------------------------------------------------------------------'
    write(*,'(A,A)') ' Fermion type:        ', trim(fermion_type)
    write(*,'(A,I3)') ' Number of flavors:   ', n_flavors
    write(*,'(A)') '---------------------------------------------------------------------'
    write(*,'(A,ES12.4)') ' Time step dt:        ', time_step
    write(*,'(A,I10)') ' Number of steps:     ', n_steps
    write(*,'(A,I10)') ' Output interval:     ', output_interval
    write(*,'(A)') '---------------------------------------------------------------------'
    write(*,'(A,L3)') ' Use OpenMP:          ', use_openmp
    write(*,'(A,L3)') ' Use MPI:             ', use_mpi
    if (use_openmp) write(*,'(A,I3)') ' Number of threads:   ', n_threads
    write(*,'(A)') '---------------------------------------------------------------------'
    write(*,'(A,A)') ' Computation method:  ', trim(computation_method)
    if (trim(computation_method) == 'tensor_network') then
      write(*,'(A,I5)') ' Bond dimension χ:    ', tensor_chi
    else if (trim(computation_method) == 'quantum') then
      write(*,'(A,A)') ' Quantum backend:     ', trim(quantum_backend)
    end if
    write(*,'(A)') '====================================================================='
    write(*,*)
    
  end subroutine print_parameters
  
  !=====================================================================
  ! Subroutine: cleanup_parameters
  ! Purpose: Deallocate arrays
  !=====================================================================
  subroutine cleanup_parameters()
    implicit none
    if (allocated(nsize)) deallocate(nsize)
  end subroutine cleanup_parameters

end module mod_parameters
