module lattice_module
    use parameter_module
    implicit none
    
    type lattice
        integer :: n_sites              ! Total number of lattice sites
        integer :: n_links              ! Total number of links
        integer :: n_dim                ! Number of dimensions
        integer :: coordination_number  ! Number of neighbors per site
        character(len=20) :: lattice_type
        
        ! Site information
        integer, allocatable :: site_coords(:,:)  ! Coordinates of each site
        integer, allocatable :: site_index(:,:,:,:) ! Reverse lookup: coords -> site index
        
        ! Neighbor information
        integer, allocatable :: neighbors(:,:)    ! neighbors(site, direction)
        real(8), allocatable :: link_vectors(:,:) ! Unit vectors for each link direction
        
        ! For hexagonal lattice
        integer :: n_sublattices        ! Number of sublattices (2 for hexagonal)
        integer, allocatable :: sublattice(:)     ! Sublattice index for each site
        
        ! Lattice spacing and volume
        real(8) :: a                    ! Lattice spacing
        real(8) :: volume              ! Physical volume
        
        ! Boundary conditions
        logical :: periodic_bc          ! Use periodic boundary conditions
    end type lattice
    
contains
    
    subroutine initialize_lattice(lat, params)
        type(lattice), intent(out) :: lat
        type(parameters), intent(in) :: params
        
        lat%n_dim = params%D_euclid
        lat%lattice_type = params%lattice_type
        lat%a = 1.0d0  ! Set lattice spacing to 1
        lat%periodic_bc = .true.  ! Default to periodic BC
        
        select case(trim(params%lattice_type))
            case("square")
                call initialize_square_lattice(lat, params)
            case("hexagonal")
                if (params%D_euclid /= 2) then
                    print *, "Error: Hexagonal lattice only supported in 2D"
                    stop
                endif
                call initialize_hexagonal_lattice(lat, params)
            case default
                print *, "Error: Unknown lattice type"
                stop
        end select
        
    end subroutine initialize_lattice
    
    subroutine initialize_square_lattice(lat, params)
        type(lattice), intent(inout) :: lat
        type(parameters), intent(in) :: params
        integer :: i, j, k, l, site_idx, d
        integer :: coords(params%D_euclid)
        integer :: neighbor_coords(params%D_euclid)
        
        ! Calculate total number of sites
        lat%n_sites = 1
        do i = 1, params%D_euclid
            lat%n_sites = lat%n_sites * params%Nsize(i)
        end do
        
        lat%coordination_number = 2 * params%D_euclid
        lat%n_links = lat%n_sites * params%D_euclid
        lat%n_sublattices = 1
        
        ! Allocate arrays
        allocate(lat%site_coords(lat%n_sites, params%D_euclid))
        allocate(lat%neighbors(lat%n_sites, lat%coordination_number))
        allocate(lat%link_vectors(params%D_euclid, lat%coordination_number))
        allocate(lat%sublattice(lat%n_sites))
        
        ! Allocate site index array based on dimensions
        select case(params%D_euclid)
            case(2)
                allocate(lat%site_index(params%Nsize(1), params%Nsize(2), 1, 1))
            case(3)
                allocate(lat%site_index(params%Nsize(1), params%Nsize(2), params%Nsize(3), 1))
            case(4)
                allocate(lat%site_index(params%Nsize(1), params%Nsize(2), &
                                      params%Nsize(3), params%Nsize(4)))
        end select
        
        ! Initialize site coordinates and index
        site_idx = 0
        do i = 1, params%Nsize(1)
            coords(1) = i - 1
            do j = 1, params%Nsize(2)
                coords(2) = j - 1
                if (params%D_euclid >= 3) then
                    do k = 1, params%Nsize(3)
                        coords(3) = k - 1
                        if (params%D_euclid == 4) then
                            do l = 1, params%Nsize(4)
                                coords(4) = l - 1
                                site_idx = site_idx + 1
                                lat%site_coords(site_idx, :) = coords
                                lat%site_index(i,j,k,l) = site_idx
                                lat%sublattice(site_idx) = 1
                            end do
                        else
                            site_idx = site_idx + 1
                            lat%site_coords(site_idx, :) = coords
                            lat%site_index(i,j,k,1) = site_idx
                            lat%sublattice(site_idx) = 1
                        endif
                    end do
                else
                    site_idx = site_idx + 1
                    lat%site_coords(site_idx, :) = coords
                    lat%site_index(i,j,1,1) = site_idx
                    lat%sublattice(site_idx) = 1
                endif
            end do
        end do
        
        ! Set up neighbors with periodic boundary conditions
        do site_idx = 1, lat%n_sites
            coords = lat%site_coords(site_idx, :)
            
            ! Forward and backward in each dimension
            do d = 1, params%D_euclid
                ! Forward neighbor
                neighbor_coords = coords
                neighbor_coords(d) = modulo(coords(d) + 1, params%Nsize(d))
                lat%neighbors(site_idx, 2*d-1) = coords_to_site(neighbor_coords, params%Nsize, params%D_euclid)
                
                ! Backward neighbor
                neighbor_coords = coords
                neighbor_coords(d) = modulo(coords(d) - 1 + params%Nsize(d), params%Nsize(d))
                lat%neighbors(site_idx, 2*d) = coords_to_site(neighbor_coords, params%Nsize, params%D_euclid)
            end do
        end do
        
        ! Set up link vectors
        lat%link_vectors = 0.0d0
        do d = 1, params%D_euclid
            lat%link_vectors(d, 2*d-1) = 1.0d0   ! Forward
            lat%link_vectors(d, 2*d) = -1.0d0    ! Backward
        end do
        
        ! Calculate volume
        lat%volume = lat%a**params%D_euclid
        do i = 1, params%D_euclid
            lat%volume = lat%volume * params%Nsize(i)
        end do
        
    end subroutine initialize_square_lattice
    
    subroutine initialize_hexagonal_lattice(lat, params)
        type(lattice), intent(inout) :: lat
        type(parameters), intent(in) :: params
        integer :: i, j, site_idx
        real(8) :: a1(2), a2(2), r(2)
        integer :: n1, n2
        
        ! Hexagonal lattice is 2D only
        if (params%D_euclid /= 2) then
            print *, "Error: Hexagonal lattice only in 2D"
            stop
        endif
        
        n1 = params%Nsize(1)
        n2 = params%Nsize(2)
        
        ! Two sublattices for hexagonal
        lat%n_sublattices = 2
        lat%n_sites = 2 * n1 * n2
        lat%coordination_number = 3  ! Each site has 3 neighbors
        lat%n_links = lat%n_sites * 3 / 2  ! Each link shared by 2 sites
        
        ! Allocate arrays
        allocate(lat%site_coords(lat%n_sites, 2))
        allocate(lat%neighbors(lat%n_sites, 3))
        allocate(lat%link_vectors(2, 3))
        allocate(lat%sublattice(lat%n_sites))
        allocate(lat%site_index(n1, n2, 2, 1))  ! Extra dimension for sublattice
        
        ! Hexagonal lattice vectors
        a1 = [1.0d0, 0.0d0]
        a2 = [0.5d0, sqrt(3.0d0)/2.0d0]
        
        ! Link vectors for hexagonal lattice (120 degrees apart)
        lat%link_vectors(:, 1) = [1.0d0, 0.0d0]
        lat%link_vectors(:, 2) = [-0.5d0, sqrt(3.0d0)/2.0d0]
        lat%link_vectors(:, 3) = [-0.5d0, -sqrt(3.0d0)/2.0d0]
        
        ! Initialize sites
        site_idx = 0
        do i = 1, n1
            do j = 1, n2
                ! Sublattice A
                site_idx = site_idx + 1
                r = (i-1)*a1 + (j-1)*a2
                lat%site_coords(site_idx, :) = [i-1, j-1]
                lat%sublattice(site_idx) = 1
                lat%site_index(i, j, 1, 1) = site_idx
                
                ! Sublattice B (shifted)
                site_idx = site_idx + 1
                lat%site_coords(site_idx, :) = [i-1, j-1]
                lat%sublattice(site_idx) = 2
                lat%site_index(i, j, 2, 1) = site_idx
            end do
        end do
        
        ! Set up neighbors for hexagonal lattice
        call setup_hexagonal_neighbors(lat, n1, n2)
        
        ! Calculate volume (area for 2D)
        lat%volume = sqrt(3.0d0)/2.0d0 * lat%a**2 * n1 * n2
        
    end subroutine initialize_hexagonal_lattice
    
    subroutine setup_hexagonal_neighbors(lat, n1, n2)
        type(lattice), intent(inout) :: lat
        integer, intent(in) :: n1, n2
        integer :: i, j, site_A, site_B
        integer :: i_next, j_next, i_prev, j_prev
        
        do i = 1, n1
            do j = 1, n2
                site_A = lat%site_index(i, j, 1, 1)
                site_B = lat%site_index(i, j, 2, 1)
                
                ! Periodic boundary conditions
                i_next = modulo(i, n1) + 1
                i_prev = modulo(i-2+n1, n1) + 1
                j_next = modulo(j, n2) + 1
                j_prev = modulo(j-2+n2, n2) + 1
                
                ! Neighbors for sublattice A
                lat%neighbors(site_A, 1) = lat%site_index(i, j, 2, 1)         ! Same unit cell B
                lat%neighbors(site_A, 2) = lat%site_index(i_prev, j, 2, 1)    ! Left B
                lat%neighbors(site_A, 3) = lat%site_index(i_prev, j_prev, 2, 1) ! Lower-left B
                
                ! Neighbors for sublattice B
                lat%neighbors(site_B, 1) = lat%site_index(i, j, 1, 1)         ! Same unit cell A
                lat%neighbors(site_B, 2) = lat%site_index(i_next, j, 1, 1)    ! Right A
                lat%neighbors(site_B, 3) = lat%site_index(i_next, j_next, 1, 1) ! Upper-right A
            end do
        end do
        
    end subroutine setup_hexagonal_neighbors
    
    integer function coords_to_site(coords, Nsize, n_dim)
        integer, intent(in) :: coords(:), Nsize(:), n_dim
        integer :: i, stride
        
        coords_to_site = 1
        stride = 1
        do i = 1, n_dim
            coords_to_site = coords_to_site + coords(i) * stride
            stride = stride * Nsize(i)
        end do
        
    end function coords_to_site
    
    subroutine cleanup_lattice(lat)
        type(lattice), intent(inout) :: lat
        
        if (allocated(lat%site_coords)) deallocate(lat%site_coords)
        if (allocated(lat%site_index)) deallocate(lat%site_index)
        if (allocated(lat%neighbors)) deallocate(lat%neighbors)
        if (allocated(lat%link_vectors)) deallocate(lat%link_vectors)
        if (allocated(lat%sublattice)) deallocate(lat%sublattice)
        
    end subroutine cleanup_lattice
    
end module lattice_module