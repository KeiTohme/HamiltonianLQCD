module hamiltonian
  use kinds
  use params
  use lattice_geometry
  use gauge_field
  use parallel_env
  implicit none
  private
  public :: run_simulation
contains

  subroutine run_simulation(p)
    type(simulation_params), intent(in) :: p
    type(lattice_t) :: lat
    type(gauge_t) :: G
    integer :: step
    real(dp) :: Ee, Eb

    call build_lattice(lat, trim(p%lattice_type), p%D_euclid, p%Nsize)
    call init_gauge(G, lat, p%Nc)

    if (get_rank() == 0) then
      write(*,'(A, I0, A, I0, A)') 'Initialized lattice: volume = ', lat%volume, ', dirs = ', lat%num_dirs, ' '
    end if

    do step = 1, p%steps
      call step_leapfrog(G, p%g, p%dt)
      if (mod(step, max(1,p%save_interval)) == 0 .or. step == 1 .or. step == p%steps) then
        call measure_energies(G, p%g, Ee, Eb)
        if (get_rank() == 0) then
          call append_energies(step, p%dt, Ee, Eb)
          write(*,'(A, I6, A, 2ES14.5)') 'Step ', step, ' energies Ee,Eb = ', Ee, Eb
        end if
      end if
    end do
  end subroutine run_simulation

  subroutine append_energies(step, dt, Ee, Eb)
    integer, intent(in) :: step
    real(dp), intent(in) :: dt, Ee, Eb
    integer :: unit, ios
    character(len=*), parameter :: fname = 'output/energies.dat'
    inquire(file=fname, exist=ios)
    if (step == 1) then
      open(newunit=unit, file=fname, status='replace', action='write', iostat=ios)
      write(unit,'(A)') '# step  t  Ee  Eb'
    else
      open(newunit=unit, file=fname, status='old', position='append', action='write', iostat=ios)
    end if
    write(unit,'(I8, 3X, F12.6, 3X, ES14.5, 3X, ES14.5)') step, dt*real(step-1,dp), Ee, Eb
    close(unit)
  end subroutine append_energies

end module hamiltonian
