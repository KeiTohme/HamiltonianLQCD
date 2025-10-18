module lattice_geometry
  use kinds
  use parallel_env
  implicit none
  private

  type :: lattice_t
     integer :: ndim
     character(len=16) :: geometry ! 'square' or 'hex'
     integer, allocatable :: size(:)
     integer :: volume
     integer :: num_dirs  ! number of positive directions per site
     logical :: is_honeycomb
  end type lattice_t

  public :: lattice_t, build_lattice, site_index, index_to_coord, neighbor_site

contains

  subroutine build_lattice(lat, geometry, ndim, size)
    type(lattice_t), intent(out) :: lat
    character(len=*), intent(in) :: geometry
    integer, intent(in) :: ndim
    integer, intent(in) :: size(:)

    integer :: d

    lat%geometry = geometry
    lat%ndim = ndim
    if (allocated(lat%size)) deallocate(lat%size)
    allocate(lat%size(ndim))
    lat%size = size

    lat%is_honeycomb = (trim(geometry) == 'hex')

    if (lat%is_honeycomb) then
      if (ndim /= 2) then
        if (get_rank() == 0) write(*,*) 'Honeycomb lattice only supports D=2. Aborting.'
        stop 1
      end if
      ! Honeycomb: two sublattices per (x,y)
      lat%num_dirs = 3
      lat%volume = 2 * product(lat%size)
    else
      lat%num_dirs = ndim
      lat%volume = product(lat%size)
    end if
  end subroutine build_lattice

  integer function site_index(lat, coord, sublattice)
    type(lattice_t), intent(in) :: lat
    integer, intent(in) :: coord(:)
    integer, intent(in), optional :: sublattice ! only for honeycomb: 0 or 1
    integer :: i, idx, stride

    if (lat%is_honeycomb) then
      if (.not. present(sublattice)) stop 'sublattice index required for honeycomb.'
      idx = coord(1)
      stride = lat%size(1)
      idx = idx + coord(2) * stride
      site_index = 2 * (mod(idx, lat%size(1)*lat%size(2))) + sublattice
    else
      idx = 0
      stride = 1
      do i = 1, lat%ndim
        idx = idx + mod(coord(i), lat%size(i)) * stride
        stride = stride * lat%size(i)
      end do
      site_index = idx + 1
    end if
  end function site_index

  subroutine index_to_coord(lat, idx, coord, sublattice)
    type(lattice_t), intent(in) :: lat
    integer, intent(in) :: idx
    integer, intent(out) :: coord(:)
    integer, intent(out), optional :: sublattice
    integer :: i, t, tmp

    if (lat%is_honeycomb) then
      if (.not. present(sublattice)) stop 'sublattice intent(out) required for honeycomb.'
      tmp = idx - 1
      sublattice = mod(tmp, 2)
      t = tmp / 2
      coord(1) = mod(t, lat%size(1))
      coord(2) = mod(t / lat%size(1), lat%size(2))
    else
      t = idx - 1
      do i = 1, lat%ndim
        coord(i) = mod(t, lat%size(i))
        t = t / lat%size(i)
      end do
    end if
  end subroutine index_to_coord

  integer function neighbor_site(lat, idx, mu, dir)
    type(lattice_t), intent(in) :: lat
    integer, intent(in) :: idx   ! 1..volume (or 1..2*Nx*Ny for honeycomb)
    integer, intent(in) :: mu    ! direction index: 1..ndim (square) or 1..3 (hex)
    integer, intent(in) :: dir   ! +1 forward, -1 backward

    integer :: coord(4), subl, nx, ny, x, y
    integer :: nxp, nyp, subl2

    if (.not. lat%is_honeycomb) then
      integer :: t, i
      t = idx - 1
      do i = 1, lat%ndim
        coord(i) = mod(t, lat%size(i))
        t = t / lat%size(i)
      end do
      coord(mu) = mod(coord(mu) + (dir>0) - (dir<0) + lat%size(mu), lat%size(mu))
      neighbor_site = 1
      do i = lat%ndim, 1, -1
        neighbor_site = (neighbor_site - 1) * lat%size(i) + coord(i) + 1
      end do
    else
      ! Honeycomb 2D: each site has 3 neighbors across sublattices
      call index_to_coord(lat, idx, coord(1:2), subl)
      nx = lat%size(1); ny = lat%size(2)
      x = coord(1); y = coord(2)
      if (subl == 0) then
        ! neighbors (to subl=1): (x,y), (x-1,y), (x,y-1)
        select case (mu)
        case (1)
          nxp = x;   nyp = y
        case (2)
          nxp = mod(x-1+nx, nx); nyp = y
        case (3)
          nxp = x;   nyp = mod(y-1+ny, ny)
        end select
        subl2 = 1
      else
        ! sublattice 1 neighbors to subl=0: (x,y), (x+1,y), (x,y+1)
        select case (mu)
        case (1)
          nxp = x;   nyp = y
        case (2)
          nxp = mod(x+1, nx); nyp = y
        case (3)
          nxp = x;   nyp = mod(y+1, ny)
        end select
        subl2 = 0
      end if
      if (dir < 0) then
        ! reverse direction means swap roles
        subl2 = 1 - subl
        ! reverse mapping: choose the unique mu' such that neighbor_site forward equals idx
        ! For simplicity, approximate by moving the opposite way using complementary rule
        if (subl == 0) then
          select case (mu)
          case (1)
            nxp = x;   nyp = y
          case (2)
            nxp = mod(x+1, nx); nyp = y
          case (3)
            nxp = x;   nyp = mod(y+1, ny)
          end select
        else
          select case (mu)
          case (1)
            nxp = x;   nyp = y
          case (2)
            nxp = mod(x-1+nx, nx); nyp = y
          case (3)
            nxp = x;   nyp = mod(y-1+ny, ny)
          end select
        end if
      end if
      neighbor_site = 2 * (nxp + nyp * nx) + subl2 + 1
    end if
  end function neighbor_site

end module lattice_geometry
