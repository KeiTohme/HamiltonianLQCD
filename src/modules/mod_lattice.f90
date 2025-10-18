!***********************************************************************
! Module: mod_lattice
! Purpose: Implement lattice structure (square and hexagonal)
!          Including site indexing, neighbor finding, and link structure
!***********************************************************************
module mod_lattice
  use mod_parameters
  implicit none
  
  ! Site and link data structures
  type :: site_type
    integer, allocatable :: coords(:)     ! Coordinates in each dimension
    integer, allocatable :: neighbors(:,:) ! Neighbor indices (dir, ±)
    integer :: index                      ! Linear index
  end type site_type
  
  type :: link_type
    integer :: site1, site2              ! Connected sites
    integer :: direction                 ! Direction index
    logical :: is_boundary = .false.     ! Boundary condition flag
  end type link_type
  
  type(site_type), allocatable :: lattice_sites(:)
  type(link_type), allocatable :: lattice_links(:)
  integer :: n_links
  integer :: n_directions
  
contains

  !=====================================================================
  ! Subroutine: initialize_lattice
  ! Purpose: Initialize lattice structure based on type
  !=====================================================================
  subroutine initialize_lattice()
    implicit none
    
    write(*,'(A)') ' Initializing lattice structure...'
    
    select case (trim(lattice_type))
    case ('square')
      call initialize_square_lattice()
    case ('hexagonal')
      call initialize_hexagonal_lattice()
    case default
      write(*,*) 'Error: Unknown lattice type: ', trim(lattice_type)
      stop
    end select
    
    write(*,'(A,I8,A)') ' Created lattice with ', total_sites, ' sites'
    write(*,'(A,I8,A)') ' and ', n_links, ' links'
    
  end subroutine initialize_lattice
  
  !=====================================================================
  ! Subroutine: initialize_square_lattice
  ! Purpose: Initialize hypercubic lattice structure
  !=====================================================================
  subroutine initialize_square_lattice()
    implicit none
    integer :: isite, idir, i, j, k, l
    integer :: coords(d_euclid)
    integer :: neighbor_coords(d_euclid)
    integer :: link_count
    
    n_directions = d_euclid
    
    ! Allocate sites
    allocate(lattice_sites(total_sites))
    
    ! Initialize sites
    isite = 0
    select case (d_euclid)
    case (2)
      do i = 0, nsize(1)-1
        do j = 0, nsize(2)-1
          isite = isite + 1
          allocate(lattice_sites(isite)%coords(d_euclid))
          allocate(lattice_sites(isite)%neighbors(d_euclid, 2))
          lattice_sites(isite)%coords = [i, j]
          lattice_sites(isite)%index = isite
        end do
      end do
      
    case (3)
      do i = 0, nsize(1)-1
        do j = 0, nsize(2)-1
          do k = 0, nsize(3)-1
            isite = isite + 1
            allocate(lattice_sites(isite)%coords(d_euclid))
            allocate(lattice_sites(isite)%neighbors(d_euclid, 2))
            lattice_sites(isite)%coords = [i, j, k]
            lattice_sites(isite)%index = isite
          end do
        end do
      end do
      
    case (4)
      do i = 0, nsize(1)-1
        do j = 0, nsize(2)-1
          do k = 0, nsize(3)-1
            do l = 0, nsize(4)-1
              isite = isite + 1
              allocate(lattice_sites(isite)%coords(d_euclid))
              allocate(lattice_sites(isite)%neighbors(d_euclid, 2))
              lattice_sites(isite)%coords = [i, j, k, l]
              lattice_sites(isite)%index = isite
            end do
          end do
        end do
      end do
      
    case default
      ! General case for arbitrary dimension
      call generate_sites_recursive(1, coords, isite)
    end select
    
    ! Set up neighbors with periodic boundary conditions
    do isite = 1, total_sites
      do idir = 1, d_euclid
        ! Forward neighbor
        neighbor_coords = lattice_sites(isite)%coords
        neighbor_coords(idir) = mod(neighbor_coords(idir) + 1, nsize(idir))
        lattice_sites(isite)%neighbors(idir, 1) = coords_to_index(neighbor_coords)
        
        ! Backward neighbor
        neighbor_coords = lattice_sites(isite)%coords
        neighbor_coords(idir) = mod(neighbor_coords(idir) - 1 + nsize(idir), nsize(idir))
        lattice_sites(isite)%neighbors(idir, 2) = coords_to_index(neighbor_coords)
      end do
    end do
    
    ! Create links
    n_links = total_sites * d_euclid
    allocate(lattice_links(n_links))
    
    link_count = 0
    do isite = 1, total_sites
      do idir = 1, d_euclid
        link_count = link_count + 1
        lattice_links(link_count)%site1 = isite
        lattice_links(link_count)%site2 = lattice_sites(isite)%neighbors(idir, 1)
        lattice_links(link_count)%direction = idir
      end do
    end do
    
  end subroutine initialize_square_lattice
  
  !=====================================================================
  ! Subroutine: initialize_hexagonal_lattice
  ! Purpose: Initialize hexagonal lattice structure (2D/3D)
  !=====================================================================
  subroutine initialize_hexagonal_lattice()
    implicit none
    integer :: isite, i, j, k
    integer :: neighbor_i, neighbor_j
    integer :: link_count
    
    if (d_euclid < 2 .or. d_euclid > 3) then
      write(*,*) 'Error: Hexagonal lattice only supported for 2D and 3D'
      stop
    end if
    
    if (d_euclid == 2) then
      n_directions = 3  ! Three primary directions in hexagonal lattice
    else
      n_directions = 4  ! For 3D hexagonal (layered)
    end if
    
    ! Allocate sites
    allocate(lattice_sites(total_sites))
    
    ! Initialize 2D hexagonal lattice
    if (d_euclid == 2) then
      isite = 0
      do i = 0, nsize(1)-1
        do j = 0, nsize(2)-1
          isite = isite + 1
          allocate(lattice_sites(isite)%coords(2))
          allocate(lattice_sites(isite)%neighbors(3, 2))
          lattice_sites(isite)%coords = [i, j]
          lattice_sites(isite)%index = isite
          
          ! Hexagonal neighbors
          ! Direction 1: right
          neighbor_i = mod(i + 1, nsize(1))
          neighbor_j = j
          lattice_sites(isite)%neighbors(1, 1) = neighbor_i * nsize(2) + neighbor_j + 1
          
          ! Direction 2: upper-right
          neighbor_i = i
          neighbor_j = mod(j + 1, nsize(2))
          lattice_sites(isite)%neighbors(2, 1) = neighbor_i * nsize(2) + neighbor_j + 1
          
          ! Direction 3: upper-left
          if (mod(j, 2) == 0) then
            neighbor_i = mod(i - 1 + nsize(1), nsize(1))
          else
            neighbor_i = i
          end if
          neighbor_j = mod(j + 1, nsize(2))
          lattice_sites(isite)%neighbors(3, 1) = neighbor_i * nsize(2) + neighbor_j + 1
          
          ! Backward neighbors (opposite directions)
          neighbor_i = mod(i - 1 + nsize(1), nsize(1))
          neighbor_j = j
          lattice_sites(isite)%neighbors(1, 2) = neighbor_i * nsize(2) + neighbor_j + 1
          
          neighbor_i = i
          neighbor_j = mod(j - 1 + nsize(2), nsize(2))
          lattice_sites(isite)%neighbors(2, 2) = neighbor_i * nsize(2) + neighbor_j + 1
          
          if (mod(j, 2) == 1) then
            neighbor_i = mod(i + 1, nsize(1))
          else
            neighbor_i = i
          end if
          neighbor_j = mod(j - 1 + nsize(2), nsize(2))
          lattice_sites(isite)%neighbors(3, 2) = neighbor_i * nsize(2) + neighbor_j + 1
        end do
      end do
    end if
    
    ! Create links
    n_links = total_sites * n_directions
    allocate(lattice_links(n_links))
    
    link_count = 0
    do isite = 1, total_sites
      do i = 1, n_directions
        link_count = link_count + 1
        lattice_links(link_count)%site1 = isite
        lattice_links(link_count)%site2 = lattice_sites(isite)%neighbors(i, 1)
        lattice_links(link_count)%direction = i
      end do
    end do
    
  end subroutine initialize_hexagonal_lattice
  
  !=====================================================================
  ! Recursive subroutine: generate_sites_recursive
  ! Purpose: Generate sites for arbitrary dimensions
  !=====================================================================
  recursive subroutine generate_sites_recursive(dim, coords, isite)
    implicit none
    integer, intent(in) :: dim
    integer, intent(inout) :: coords(:)
    integer, intent(inout) :: isite
    integer :: i
    
    if (dim > d_euclid) then
      isite = isite + 1
      allocate(lattice_sites(isite)%coords(d_euclid))
      allocate(lattice_sites(isite)%neighbors(d_euclid, 2))
      lattice_sites(isite)%coords = coords
      lattice_sites(isite)%index = isite
      return
    end if
    
    do i = 0, nsize(dim)-1
      coords(dim) = i
      call generate_sites_recursive(dim + 1, coords, isite)
    end do
    
  end subroutine generate_sites_recursive
  
  !=====================================================================
  ! Function: coords_to_index
  ! Purpose: Convert lattice coordinates to linear index
  !=====================================================================
  function coords_to_index(coords) result(idx)
    implicit none
    integer, intent(in) :: coords(:)
    integer :: idx
    integer :: i, stride
    
    idx = coords(1) + 1
    stride = nsize(1)
    
    do i = 2, d_euclid
      idx = idx + coords(i) * stride
      stride = stride * nsize(i)
    end do
    
  end function coords_to_index
  
  !=====================================================================
  ! Function: index_to_coords
  ! Purpose: Convert linear index to lattice coordinates
  !=====================================================================
  function index_to_coords(idx) result(coords)
    implicit none
    integer, intent(in) :: idx
    integer :: coords(d_euclid)
    integer :: i, temp, stride
    
    temp = idx - 1
    stride = 1
    
    do i = 1, d_euclid
      coords(i) = mod(temp / stride, nsize(i))
      stride = stride * nsize(i)
    end do
    
  end function index_to_coords
  
  !=====================================================================
  ! Function: get_plaquette_sites
  ! Purpose: Get the four sites of a plaquette in directions mu, nu
  !=====================================================================
  subroutine get_plaquette_sites(isite, mu, nu, plaq_sites)
    implicit none
    integer, intent(in) :: isite, mu, nu
    integer, intent(out) :: plaq_sites(4)
    
    plaq_sites(1) = isite
    plaq_sites(2) = lattice_sites(isite)%neighbors(mu, 1)
    plaq_sites(3) = lattice_sites(plaq_sites(2))%neighbors(nu, 1)
    plaq_sites(4) = lattice_sites(isite)%neighbors(nu, 1)
    
  end subroutine get_plaquette_sites
  
  !=====================================================================
  ! Subroutine: cleanup_lattice
  ! Purpose: Deallocate lattice arrays
  !=====================================================================
  subroutine cleanup_lattice()
    implicit none
    integer :: i
    
    if (allocated(lattice_sites)) then
      do i = 1, total_sites
        if (allocated(lattice_sites(i)%coords)) deallocate(lattice_sites(i)%coords)
        if (allocated(lattice_sites(i)%neighbors)) deallocate(lattice_sites(i)%neighbors)
      end do
      deallocate(lattice_sites)
    end if
    
    if (allocated(lattice_links)) deallocate(lattice_links)
    
  end subroutine cleanup_lattice

end module mod_lattice
