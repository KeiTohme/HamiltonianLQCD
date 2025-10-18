program lqcd_main
  use kinds
  use params
  use parallel_env
  use hamiltonian
  implicit none

  type(simulation_params) :: p

  call init_parallel()

  if (get_rank() == 0) then
    call read_params('params.in', p)
  end if
  call broadcast_params(p)

  if (get_rank() == 0) then
    write(*,'(A, F8.4, A, F8.4, A, I0, A, I0)') 'g=', p%g, ' m=', p%m, ' D=', p%D_euclid, ' Nc=', p%Nc
    write(*,'(A, A)') 'lattice=', trim(p%lattice_type)
  end if

  call run_simulation(p)

  call finalize_parallel()
end program lqcd_main
