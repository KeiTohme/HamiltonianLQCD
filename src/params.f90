module params
  use kinds
  use parallel_env
  implicit none
  private

  type :: simulation_params
     real(dp) :: g = 1.0_dp
     real(dp) :: m = 0.0_dp
     integer  :: D_euclid = 2
     integer  :: Nc = 2
     integer  :: steps = 10
     real(dp) :: dt = 0.01_dp
     character(len=16) :: lattice_type = 'square'  ! 'square' or 'hex'
     integer, allocatable :: Nsize(:)              ! length D_euclid
     logical :: enable_mpi = .false.
     logical :: enable_openmp = .false.
     integer :: gauge_seed = 12345
     character(len=16) :: fermion_type = 'Wilson'  ! 'Wilson' or 'Staggered'
     integer :: save_interval = 10
  end type simulation_params

  public :: simulation_params, read_params, broadcast_params

contains

  subroutine read_params(filename, p)
    character(len=*), intent(in) :: filename
    type(simulation_params), intent(inout) :: p

    integer :: ios
    integer :: D_euclid_in
    integer :: Nc_in
    integer :: steps_in
    integer :: Nsize_in(6)
    real(dp) :: g_in, m_in, dt_in
    character(len=16) :: lattice_type_in, fermion_type_in
    logical :: enable_mpi_in, enable_openmp_in
    integer :: gauge_seed_in, save_interval_in

    ! Defaults
    D_euclid_in = p%D_euclid
    Nc_in = p%Nc
    steps_in = p%steps
    Nsize_in = 1
    Nsize_in(1) = 8
    Nsize_in(2) = 8
    g_in = p%g
    m_in = p%m
    dt_in = p%dt
    lattice_type_in = p%lattice_type
    fermion_type_in = p%fermion_type
    enable_mpi_in = p%enable_mpi
    enable_openmp_in = p%enable_openmp
    gauge_seed_in = p%gauge_seed
    save_interval_in = p%save_interval

    namelist /simulation/ g_in, m_in, D_euclid_in, Nc_in, steps_in, dt_in, &
         lattice_type_in, Nsize_in, enable_mpi_in, enable_openmp_in, &
         gauge_seed_in, fermion_type_in, save_interval_in

    open(unit=10, file=filename, status='old', action='read', iostat=ios)
    if (ios /= 0) then
      if (get_rank() == 0) then
        write(*,*) 'Could not open params file: ', trim(filename), ' using defaults.'
      end if
    else
      read(10, nml=simulation, iostat=ios)
      if (ios /= 0) then
        if (get_rank() == 0) then
          write(*,*) 'Error reading params file; using defaults.'
        end if
      end if
      close(10)
    end if

    p%D_euclid = D_euclid_in
    p%Nc = Nc_in
    p%steps = steps_in
    p%g = g_in
    p%m = m_in
    p%dt = dt_in
    p%lattice_type = lattice_type_in
    p%fermion_type = fermion_type_in
    p%enable_mpi = enable_mpi_in
    p%enable_openmp = enable_openmp_in
    p%gauge_seed = gauge_seed_in
    p%save_interval = save_interval_in

    if (allocated(p%Nsize)) deallocate(p%Nsize)
    allocate(p%Nsize(p%D_euclid))
    p%Nsize = Nsize_in(1:p%D_euclid)

  end subroutine read_params

  subroutine broadcast_params(p)
    type(simulation_params), intent(inout) :: p
#ifdef USE_MPI
    integer :: ierr
    integer :: root
    root = 0
    call MPI_Bcast(p%g, 1, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(p%m, 1, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(p%dt, 1, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(p%D_euclid, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(p%Nc, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(p%steps, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(p%gauge_seed, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(p%save_interval, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
    ! Broadcast char arrays as bytes
    call MPI_Bcast(p%lattice_type, len(p%lattice_type), MPI_CHARACTER, root, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(p%fermion_type, len(p%fermion_type), MPI_CHARACTER, root, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(p%enable_mpi, 1, MPI_LOGICAL, root, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(p%enable_openmp, 1, MPI_LOGICAL, root, MPI_COMM_WORLD, ierr)
    if (allocated(p%Nsize)) then
      call MPI_Bcast(p%Nsize, size(p%Nsize), MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
    end if
#endif
  end subroutine broadcast_params

end module params
