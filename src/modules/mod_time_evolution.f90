!***********************************************************************
! Module: mod_time_evolution
! Purpose: Implement real-time evolution for gluon dynamics
!          Using various integration schemes (leapfrog, RK4, etc.)
!***********************************************************************
module mod_time_evolution
  use mod_parameters
  use mod_lattice
  use mod_gauge_field
  use mod_hamiltonian
  use mod_su_algebra
  implicit none
  
  ! Time evolution variables
  real(dp) :: current_time = 0.0_dp
  integer :: current_step = 0
  
  ! Observables history
  real(dp), allocatable :: energy_history(:)
  real(dp), allocatable :: plaquette_history(:)
  real(dp), allocatable :: time_array(:)
  
contains

  !=====================================================================
  ! Subroutine: initialize_time_evolution
  ! Purpose: Initialize time evolution arrays and observables
  !=====================================================================
  subroutine initialize_time_evolution()
    implicit none
    
    write(*,'(A)') ' Initializing time evolution...'
    
    ! Allocate history arrays
    if (allocated(energy_history)) deallocate(energy_history)
    if (allocated(plaquette_history)) deallocate(plaquette_history)
    if (allocated(time_array)) deallocate(time_array)
    
    allocate(energy_history(n_steps+1))
    allocate(plaquette_history(n_steps+1))
    allocate(time_array(n_steps+1))
    
    energy_history = 0.0_dp
    plaquette_history = 0.0_dp
    time_array = 0.0_dp
    
    current_time = 0.0_dp
    current_step = 0
    
    ! Record initial values
    call record_observables()
    
    write(*,'(A,ES12.4)') '   Time step dt = ', time_step
    write(*,'(A,I8)') '   Number of steps = ', n_steps
    write(*,'(A,ES12.4)') '   Total time = ', time_step * n_steps
    
  end subroutine initialize_time_evolution
  
  !=====================================================================
  ! Subroutine: evolve_leapfrog
  ! Purpose: Evolve system using leapfrog (symplectic) integrator
  !          Preserves energy and phase space volume
  !=====================================================================
  subroutine evolve_leapfrog()
    implicit none
    integer :: istep
    real(dp) :: dt_half
    complex(dp) :: dU_dt(n_colors, n_colors, n_links)
    complex(dp) :: dE_dt(n_colors**2-1, n_links)
    
    dt_half = 0.5_dp * time_step
    
    write(*,'(A)') ' Starting leapfrog time evolution...'
    write(*,'(A)') ' Step      Time         Energy          Plaquette'
    write(*,'(A)') ' --------------------------------------------------------'
    
    do istep = 1, n_steps
      ! Half step for E: E(t+dt/2) = E(t) + (dt/2) dE/dt
      call calculate_equations_of_motion(dU_dt, dE_dt)
      electric_field = electric_field + dt_half * dE_dt
      
      ! Full step for U: U(t+dt) = U(t) + dt dU/dt
      call calculate_equations_of_motion(dU_dt, dE_dt)
      link_matrices = link_matrices + time_step * dU_dt
      
      ! Project back to SU(Nc)
      call project_links_to_sun()
      
      ! Half step for E: E(t+dt) = E(t+dt/2) + (dt/2) dE/dt
      call calculate_equations_of_motion(dU_dt, dE_dt)
      electric_field = electric_field + dt_half * dE_dt
      
      ! Update time
      current_step = istep
      current_time = current_time + time_step
      
      ! Record observables
      call record_observables()
      
      ! Output progress
      if (mod(istep, output_interval) == 0) then
        write(*,'(I6,3ES16.6)') istep, current_time, &
          energy_history(istep+1), plaquette_history(istep+1)
      end if
    end do
    
    write(*,'(A)') ' --------------------------------------------------------'
    write(*,'(A)') ' Time evolution completed'
    
  end subroutine evolve_leapfrog
  
  !=====================================================================
  ! Subroutine: evolve_rk4
  ! Purpose: Evolve system using 4th-order Runge-Kutta
  !=====================================================================
  subroutine evolve_rk4()
    implicit none
    integer :: istep
    complex(dp) :: U0(n_colors, n_colors, n_links)
    complex(dp) :: E0(n_colors**2-1, n_links)
    complex(dp) :: k1_U(n_colors, n_colors, n_links)
    complex(dp) :: k2_U(n_colors, n_colors, n_links)
    complex(dp) :: k3_U(n_colors, n_colors, n_links)
    complex(dp) :: k4_U(n_colors, n_colors, n_links)
    complex(dp) :: k1_E(n_colors**2-1, n_links)
    complex(dp) :: k2_E(n_colors**2-1, n_links)
    complex(dp) :: k3_E(n_colors**2-1, n_links)
    complex(dp) :: k4_E(n_colors**2-1, n_links)
    
    write(*,'(A)') ' Starting RK4 time evolution...'
    write(*,'(A)') ' Step      Time         Energy          Plaquette'
    write(*,'(A)') ' --------------------------------------------------------'
    
    do istep = 1, n_steps
      ! Save initial values
      U0 = link_matrices
      E0 = electric_field
      
      ! k1
      call calculate_equations_of_motion(k1_U, k1_E)
      
      ! k2
      link_matrices = U0 + 0.5_dp * time_step * k1_U
      electric_field = E0 + 0.5_dp * time_step * k1_E
      call calculate_equations_of_motion(k2_U, k2_E)
      
      ! k3
      link_matrices = U0 + 0.5_dp * time_step * k2_U
      electric_field = E0 + 0.5_dp * time_step * k2_E
      call calculate_equations_of_motion(k3_U, k3_E)
      
      ! k4
      link_matrices = U0 + time_step * k3_U
      electric_field = E0 + time_step * k3_E
      call calculate_equations_of_motion(k4_U, k4_E)
      
      ! Final update
      link_matrices = U0 + time_step * (k1_U + 2.0_dp*k2_U + 2.0_dp*k3_U + k4_U) / 6.0_dp
      electric_field = E0 + time_step * (k1_E + 2.0_dp*k2_E + 2.0_dp*k3_E + k4_E) / 6.0_dp
      
      ! Project back to SU(Nc)
      call project_links_to_sun()
      
      ! Update time
      current_step = istep
      current_time = current_time + time_step
      
      ! Record observables
      call record_observables()
      
      ! Output progress
      if (mod(istep, output_interval) == 0) then
        write(*,'(I6,3ES16.6)') istep, current_time, &
          energy_history(istep+1), plaquette_history(istep+1)
      end if
    end do
    
    write(*,'(A)') ' --------------------------------------------------------'
    write(*,'(A)') ' Time evolution completed'
    
  end subroutine evolve_rk4
  
  !=====================================================================
  ! Subroutine: project_links_to_sun
  ! Purpose: Project all link matrices back to SU(Nc)
  !=====================================================================
  subroutine project_links_to_sun()
    implicit none
    integer :: ilink
    
    do ilink = 1, n_links
      call project_to_sun(link_matrices(:,:,ilink))
    end do
    
  end subroutine project_links_to_sun
  
  !=====================================================================
  ! Subroutine: record_observables
  ! Purpose: Record current observables for analysis
  !=====================================================================
  subroutine record_observables()
    implicit none
    real(dp) :: H
    
    time_array(current_step+1) = current_time
    
    ! Calculate Hamiltonian
    H = calculate_hamiltonian()
    energy_history(current_step+1) = H
    
    ! Calculate average plaquette
    plaquette_history(current_step+1) = plaquette_average()
    
  end subroutine record_observables
  
  !=====================================================================
  ! Subroutine: write_observables
  ! Purpose: Write observables to file
  !=====================================================================
  subroutine write_observables(filename)
    implicit none
    character(len=*), intent(in) :: filename
    integer :: unit_num, i
    
    open(newunit=unit_num, file=filename, status='replace', action='write')
    
    write(unit_num, '(A)') '# Time evolution data'
    write(unit_num, '(A)') '# Columns: Step, Time, Energy, Plaquette'
    
    do i = 1, current_step + 1
      write(unit_num, '(I8,3ES20.10)') i-1, time_array(i), &
        energy_history(i), plaquette_history(i)
    end do
    
    close(unit_num)
    
    write(*,'(A,A)') ' Observables written to: ', trim(filename)
    
  end subroutine write_observables
  
  !=====================================================================
  ! Subroutine: write_configuration
  ! Purpose: Write current gauge configuration to file
  !=====================================================================
  subroutine write_configuration(filename)
    implicit none
    character(len=*), intent(in) :: filename
    integer :: unit_num, ilink, i, j
    
    open(newunit=unit_num, file=filename, status='replace', &
         action='write', form='unformatted')
    
    ! Write header
    write(unit_num) n_colors, n_links, current_step, current_time
    
    ! Write link matrices
    do ilink = 1, n_links
      do i = 1, n_colors
        do j = 1, n_colors
          write(unit_num) link_matrices(i, j, ilink)
        end do
      end do
    end do
    
    ! Write electric field
    write(unit_num) electric_field
    
    close(unit_num)
    
    write(*,'(A,A)') ' Configuration written to: ', trim(filename)
    
  end subroutine write_configuration
  
  !=====================================================================
  ! Subroutine: read_configuration
  ! Purpose: Read gauge configuration from file
  !=====================================================================
  subroutine read_configuration(filename)
    implicit none
    character(len=*), intent(in) :: filename
    integer :: unit_num, ilink, i, j
    integer :: nc_file, nl_file, step_file
    real(dp) :: time_file
    
    open(newunit=unit_num, file=filename, status='old', &
         action='read', form='unformatted')
    
    ! Read header
    read(unit_num) nc_file, nl_file, step_file, time_file
    
    if (nc_file /= n_colors .or. nl_file /= n_links) then
      write(*,*) 'Error: Configuration file dimensions do not match'
      close(unit_num)
      return
    end if
    
    current_step = step_file
    current_time = time_file
    
    ! Read link matrices
    do ilink = 1, n_links
      do i = 1, n_colors
        do j = 1, n_colors
          read(unit_num) link_matrices(i, j, ilink)
        end do
      end do
    end do
    
    ! Read electric field
    read(unit_num) electric_field
    
    close(unit_num)
    
    write(*,'(A,A)') ' Configuration read from: ', trim(filename)
    write(*,'(A,I8,A,F12.4)') ' Restored step ', current_step, ' at time ', current_time
    
  end subroutine read_configuration
  
  !=====================================================================
  ! Subroutine: cleanup_time_evolution
  ! Purpose: Deallocate time evolution arrays
  !=====================================================================
  subroutine cleanup_time_evolution()
    implicit none
    if (allocated(energy_history)) deallocate(energy_history)
    if (allocated(plaquette_history)) deallocate(plaquette_history)
    if (allocated(time_array)) deallocate(time_array)
  end subroutine cleanup_time_evolution

end module mod_time_evolution
