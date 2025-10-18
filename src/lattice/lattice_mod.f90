! ==============================================================================
! Module: lattice_mod
! Description: Lattice structure and neighbor relations
! ==============================================================================
module lattice_mod
    use parameters_mod
    implicit none
    
    ! Lattice site structure
    type :: lattice_site
        integer :: site_index
        integer, allocatable :: coords(:)
        integer, allocatable :: neighbors(:,:)  ! (direction, forward/backward)
    end type lattice_site
    
    ! Lattice structure
    type(lattice_site), allocatable :: lattice(:)
    
    ! Lattice type flag for efficient checking
    logical :: is_square_lattice
    
contains

    ! ==========================================================================
    ! Subroutine: initialize_lattice
    ! Description: Initialize lattice structure based on type
    ! ==========================================================================
    subroutine initialize_lattice()
        implicit none
        
        select case (trim(lattice_type))
            case ('square')
                call initialize_square_lattice()
                is_square_lattice = .true.
            case ('hexagonal')
                call initialize_hexagonal_lattice()
                is_square_lattice = .false.
            case default
                print *, 'Error: Unknown lattice type: ', trim(lattice_type)
                stop
        end select
        
    end subroutine initialize_lattice
    
    ! ==========================================================================
    ! Subroutine: initialize_square_lattice
    ! Description: Initialize square/hypercubic lattice
    ! ==========================================================================
    subroutine initialize_square_lattice()
        implicit none
        integer :: i, j, mu, site_idx
        integer :: coords(spacetime_dim)
        integer :: neighbor_coords(spacetime_dim)
        
        ! Allocate lattice sites
        allocate(lattice(total_sites))
        
        ! Initialize each site
        do site_idx = 1, total_sites
            lattice(site_idx)%site_index = site_idx
            allocate(lattice(site_idx)%coords(spacetime_dim))
            allocate(lattice(site_idx)%neighbors(spacetime_dim, 2))  ! 2 for forward/backward
            
            ! Convert linear index to coordinates
            call index_to_coords(site_idx, lattice(site_idx)%coords)
            
            ! Find neighbors
            do mu = 1, spacetime_dim
                ! Forward neighbor
                neighbor_coords = lattice(site_idx)%coords
                neighbor_coords(mu) = mod(neighbor_coords(mu), lattice_size(mu)) + 1
                lattice(site_idx)%neighbors(mu, 1) = coords_to_index(neighbor_coords)
                
                ! Backward neighbor
                neighbor_coords = lattice(site_idx)%coords
                neighbor_coords(mu) = neighbor_coords(mu) - 1
                if (neighbor_coords(mu) < 1) neighbor_coords(mu) = lattice_size(mu)
                lattice(site_idx)%neighbors(mu, 2) = coords_to_index(neighbor_coords)
            end do
        end do
        
        print *, 'Square lattice initialized with ', total_sites, ' sites'
        
    end subroutine initialize_square_lattice
    
    ! ==========================================================================
    ! Subroutine: initialize_hexagonal_lattice
    ! Description: Initialize hexagonal lattice (2D case)
    ! ==========================================================================
    subroutine initialize_hexagonal_lattice()
        implicit none
        integer :: i, j, site_idx
        
        ! For simplicity, we implement a 2D hexagonal lattice
        ! In hexagonal lattice, each site has 6 nearest neighbors
        
        if (spacetime_dim /= 2) then
            print *, 'Warning: Hexagonal lattice currently only supports 2D'
            print *, 'Falling back to square lattice'
            call initialize_square_lattice()
            is_square_lattice = .true.  ! Update flag since we're using square lattice
            return
        end if
        
        ! Allocate lattice sites
        allocate(lattice(total_sites))
        
        ! Initialize each site with 6 neighbors for hexagonal structure
        do site_idx = 1, total_sites
            lattice(site_idx)%site_index = site_idx
            allocate(lattice(site_idx)%coords(spacetime_dim))
            allocate(lattice(site_idx)%neighbors(6, 1))  ! 6 neighbors
            
            ! Convert linear index to coordinates
            call index_to_coords(site_idx, lattice(site_idx)%coords)
            
            ! Set up hexagonal neighbors (simplified version)
            call setup_hexagonal_neighbors(site_idx)
        end do
        
        print *, 'Hexagonal lattice initialized with ', total_sites, ' sites'
        
    end subroutine initialize_hexagonal_lattice
    
    ! ==========================================================================
    ! Subroutine: setup_hexagonal_neighbors
    ! Description: Setup neighbors for hexagonal lattice
    ! ==========================================================================
    subroutine setup_hexagonal_neighbors(site_idx)
        implicit none
        integer, intent(in) :: site_idx
        integer :: i, j
        integer :: neighbor_coords(2)
        integer :: parity
        
        i = lattice(site_idx)%coords(1)
        j = lattice(site_idx)%coords(2)
        parity = mod(i, 2)
        
        ! Six neighbors in hexagonal lattice
        ! Right
        neighbor_coords = [mod(i, lattice_size(1)) + 1, j]
        lattice(site_idx)%neighbors(1, 1) = coords_to_index(neighbor_coords)
        
        ! Left
        neighbor_coords(1) = i - 1
        if (neighbor_coords(1) < 1) neighbor_coords(1) = lattice_size(1)
        neighbor_coords(2) = j
        lattice(site_idx)%neighbors(2, 1) = coords_to_index(neighbor_coords)
        
        ! Upper right
        neighbor_coords(1) = i
        if (parity == 0) then
            neighbor_coords(1) = mod(i, lattice_size(1)) + 1
        end if
        neighbor_coords(2) = mod(j, lattice_size(2)) + 1
        lattice(site_idx)%neighbors(3, 1) = coords_to_index(neighbor_coords)
        
        ! Upper left
        neighbor_coords(1) = i
        if (parity == 1) then
            neighbor_coords(1) = i - 1
            if (neighbor_coords(1) < 1) neighbor_coords(1) = lattice_size(1)
        end if
        neighbor_coords(2) = mod(j, lattice_size(2)) + 1
        lattice(site_idx)%neighbors(4, 1) = coords_to_index(neighbor_coords)
        
        ! Lower right
        neighbor_coords(1) = i
        if (parity == 0) then
            neighbor_coords(1) = mod(i, lattice_size(1)) + 1
        end if
        neighbor_coords(2) = j - 1
        if (neighbor_coords(2) < 1) neighbor_coords(2) = lattice_size(2)
        lattice(site_idx)%neighbors(5, 1) = coords_to_index(neighbor_coords)
        
        ! Lower left
        neighbor_coords(1) = i
        if (parity == 1) then
            neighbor_coords(1) = i - 1
            if (neighbor_coords(1) < 1) neighbor_coords(1) = lattice_size(1)
        end if
        neighbor_coords(2) = j - 1
        if (neighbor_coords(2) < 1) neighbor_coords(2) = lattice_size(2)
        lattice(site_idx)%neighbors(6, 1) = coords_to_index(neighbor_coords)
        
    end subroutine setup_hexagonal_neighbors
    
    ! ==========================================================================
    ! Function: index_to_coords
    ! Description: Convert linear index to lattice coordinates
    ! ==========================================================================
    subroutine index_to_coords(idx, coords)
        implicit none
        integer, intent(in) :: idx
        integer, intent(out) :: coords(:)
        integer :: i, temp
        
        temp = idx - 1
        do i = 1, spacetime_dim
            coords(i) = mod(temp, lattice_size(i)) + 1
            temp = temp / lattice_size(i)
        end do
        
    end subroutine index_to_coords
    
    ! ==========================================================================
    ! Function: coords_to_index
    ! Description: Convert lattice coordinates to linear index
    ! ==========================================================================
    function coords_to_index(coords) result(idx)
        implicit none
        integer, intent(in) :: coords(:)
        integer :: idx
        integer :: i, multiplier
        
        idx = coords(1)
        multiplier = lattice_size(1)
        
        do i = 2, spacetime_dim
            idx = idx + (coords(i) - 1) * multiplier
            multiplier = multiplier * lattice_size(i)
        end do
        
    end function coords_to_index
    
    ! ==========================================================================
    ! Subroutine: cleanup_lattice
    ! Description: Cleanup lattice memory
    ! ==========================================================================
    subroutine cleanup_lattice()
        implicit none
        integer :: i
        
        if (allocated(lattice)) then
            do i = 1, total_sites
                if (allocated(lattice(i)%coords)) deallocate(lattice(i)%coords)
                if (allocated(lattice(i)%neighbors)) deallocate(lattice(i)%neighbors)
            end do
            deallocate(lattice)
        end if
        
    end subroutine cleanup_lattice

end module lattice_mod
