!***********************************************************************
! Program: Lattice QCD Hamiltonian Formalism
! Purpose: Main program for real-time gluon dynamics simulation
!          using Hamiltonian formulation of Lattice QCD
!
! Features:
!   - Square and hexagonal lattice structures
!   - SU(Nc) gauge theory with arbitrary Nc
!   - Wilson and Staggered fermion formulations
!   - Real-time Hamiltonian evolution
!   - OpenMP/MPI parallelization
!   - Tensor Network and Quantum Computing interfaces
!***********************************************************************
program lattice_qcd_hamiltonian
  use mod_parameters
  use mod_su_algebra
  use mod_lattice
  use mod_gauge_field
  use mod_wilson_fermion
  use mod_staggered_fermion
  use mod_hamiltonian
  use mod_time_evolution
  use mod_parallel
  use mod_tensor_network
  use mod_quantum_interface
  implicit none
  
  character(len=256) :: input_file
  character(len=256) :: output_file
  real(dp) :: start_time, end_time
  real(dp) :: initial_energy, final_energy
  
  ! Print header
  call print_header()
  
  ! Get input file from command line or use default
  if (command_argument_count() >= 1) then
    call get_command_argument(1, input_file)
  else
    input_file = 'input/parameters.inp'
  end if
  
  ! ====================================================================
  ! Initialization Phase
  ! ====================================================================
  
  call cpu_time(start_time)
  
  ! Read parameters
  call read_parameters(input_file)
  
  ! Initialize parallel environment
  call initialize_parallel()
  
  ! Initialize SU(Nc) algebra
  call initialize_su_algebra()
  
  ! Initialize lattice structure
  call initialize_lattice()
  
  ! Initialize gauge field
  call initialize_gauge_field('cold')  ! or 'hot' for random start
  
  ! Initialize fermion fields
  if (trim(fermion_type) == 'wilson') then
    call initialize_wilson_fermion()
  else if (trim(fermion_type) == 'staggered') then
    call initialize_staggered_fermion()
  end if
  
  ! Calculate initial observables
  initial_energy = calculate_hamiltonian()
  write(*,'(A,ES16.8)') ' Initial Hamiltonian: ', initial_energy
  write(*,'(A,ES16.8)') '   Electric energy:   ', H_electric
  write(*,'(A,ES16.8)') '   Magnetic energy:   ', H_magnetic
  write(*,'(A,ES16.8)') '   Fermion energy:    ', H_fermion
  write(*,'(A,ES16.8)') ' Initial plaquette:   ', plaquette_average()
  write(*,*)
  
  ! ====================================================================
  ! Select Computation Method
  ! ====================================================================
  
  select case (trim(computation_method))
    
  case ('classical')
    ! Classical Hamiltonian evolution
    write(*,'(A)') ' Using classical Hamiltonian evolution'
    write(*,*)
    
    call initialize_time_evolution()
    
    ! Choose integrator
    call evolve_leapfrog()  ! Symplectic integrator
    ! call evolve_rk4()     ! Alternative: 4th order Runge-Kutta
    
    ! Write results
    call write_observables('output/observables.dat')
    call write_configuration('output/final_config.dat')
    
  case ('tensor_network')
    ! Tensor Network method
    write(*,'(A)') ' Using Tensor Network method'
    write(*,*)
    
    call initialize_tensor_network()
    call gauge_to_tensor()
    
    ! Time evolution using TEBD or similar
    write(*,'(A)') ' Note: Full tensor network evolution requires external library'
    write(*,'(A)') ' (e.g., ITensor, TeNPy, etc.)'
    write(*,*)
    
  case ('quantum')
    ! Quantum computing method
    write(*,'(A)') ' Using Quantum Computing method'
    write(*,*)
    
    call initialize_quantum_interface()
    call gauge_to_quantum()
    
    write(*,'(A)') ' Note: Quantum execution requires external quantum backend'
    write(*,'(A)') ' (e.g., Qiskit, Cirq, Q#, etc.)'
    write(*,*)
    
  case default
    write(*,*) 'Error: Unknown computation method: ', trim(computation_method)
    stop
  end select
  
  ! ====================================================================
  ! Analysis and Output
  ! ====================================================================
  
  ! Calculate final observables
  final_energy = calculate_hamiltonian()
  write(*,*)
  write(*,'(A)') ' ============================================================='
  write(*,'(A)') '                     Final Results'
  write(*,'(A)') ' ============================================================='
  write(*,'(A,ES16.8)') ' Final Hamiltonian:   ', final_energy
  write(*,'(A,ES16.8)') '   Electric energy:   ', H_electric
  write(*,'(A,ES16.8)') '   Magnetic energy:   ', H_magnetic
  write(*,'(A,ES16.8)') '   Fermion energy:    ', H_fermion
  write(*,'(A,ES16.8)') ' Final plaquette:     ', plaquette_average()
  write(*,'(A,ES16.8)') ' Energy change:       ', final_energy - initial_energy
  write(*,'(A,ES16.8)') ' Relative change:     ', &
    (final_energy - initial_energy) / abs(initial_energy)
  
  ! Fermion observables
  if (trim(fermion_type) == 'wilson') then
    write(*,'(A,2ES16.8)') ' Fermion bilinear ψ̄ψ:', fermion_bilinear(5)
  else if (trim(fermion_type) == 'staggered') then
    write(*,'(A,2ES16.8)') ' Chiral condensate:   ', staggered_condensate()
  end if
  
  write(*,'(A)') ' ============================================================='
  
  ! ====================================================================
  ! Cleanup
  ! ====================================================================
  
  call cleanup_time_evolution()
  call cleanup_gauge_field()
  if (trim(fermion_type) == 'wilson') then
    call cleanup_wilson_fermion()
  else if (trim(fermion_type) == 'staggered') then
    call cleanup_staggered_fermion()
  end if
  call cleanup_lattice()
  call cleanup_su_algebra()
  call cleanup_parameters()
  
  if (trim(computation_method) == 'tensor_network') then
    call cleanup_tensor_network()
  else if (trim(computation_method) == 'quantum') then
    call cleanup_quantum_interface()
  end if
  
  call finalize_parallel()
  
  ! Print timing
  call cpu_time(end_time)
  write(*,*)
  write(*,'(A,F12.3,A)') ' Total CPU time: ', end_time - start_time, ' seconds'
  write(*,*)
  
  write(*,'(A)') ' Program completed successfully'
  
contains

  !=====================================================================
  ! Subroutine: print_header
  ! Purpose: Print program header
  !=====================================================================
  subroutine print_header()
    implicit none
    
    write(*,*)
    write(*,'(A)') ' ============================================================='
    write(*,'(A)') '        Lattice QCD - Hamiltonian Formalism                  '
    write(*,'(A)') '        Real-Time Gluon Dynamics Simulation                   '
    write(*,'(A)') ' ============================================================='
    write(*,'(A)') '  Features:                                                   '
    write(*,'(A)') '    - Square and Hexagonal Lattices                           '
    write(*,'(A)') '    - SU(Nc) Gauge Theory                                     '
    write(*,'(A)') '    - Wilson and Staggered Fermions                           '
    write(*,'(A)') '    - Real-Time Hamiltonian Evolution                         '
    write(*,'(A)') '    - OpenMP/MPI Parallelization                              '
    write(*,'(A)') '    - Tensor Network Interface                                '
    write(*,'(A)') '    - Quantum Computing Interface                             '
    write(*,'(A)') ' ============================================================='
    write(*,*)
    
  end subroutine print_header

end program lattice_qcd_hamiltonian
