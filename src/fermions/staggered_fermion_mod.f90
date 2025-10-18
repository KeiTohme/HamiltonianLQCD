! ==============================================================================
! Module: staggered_fermion_mod
! Description: Staggered (Kogut-Susskind) fermion formulation on the lattice
! ==============================================================================
module staggered_fermion_mod
    use parameters_mod
    use lattice_mod
    use su_n_mod
    implicit none
    
    ! Staggered fermion field (one component per site and color)
    ! chi(site, color)
    complex(dp), allocatable :: chi(:,:)
    complex(dp), allocatable :: chi_bar(:,:)
    
contains

    ! ==========================================================================
    ! Subroutine: initialize_staggered_fermions
    ! Description: Initialize staggered fermion fields
    ! ==========================================================================
    subroutine initialize_staggered_fermions()
        implicit none
        
        ! Allocate fermion fields (one component per site, color index)
        allocate(chi(total_sites, num_colors))
        allocate(chi_bar(total_sites, num_colors))
        
        ! Initialize to zero
        chi = cmplx(0.0_dp, 0.0_dp, dp)
        chi_bar = cmplx(0.0_dp, 0.0_dp, dp)
        
        print *, 'Staggered fermion fields initialized'
        
    end subroutine initialize_staggered_fermions
    
    ! ==========================================================================
    ! Subroutine: apply_staggered_dirac_operator
    ! Description: Apply staggered Dirac operator to fermion field
    ! D_staggered = m * chi(x) + 
    !               sum_mu eta_mu(x) [U_mu(x) chi(x+mu) - U_mu^dagger(x-mu) chi(x-mu)] / 2
    ! ==========================================================================
    subroutine apply_staggered_dirac_operator(chi_in, chi_out)
        implicit none
        complex(dp), intent(in) :: chi_in(:,:)
        complex(dp), intent(out) :: chi_out(:,:)
        integer :: site, mu, neighbor_fwd, neighbor_bwd
        real(dp) :: eta_mu
        complex(dp) :: temp_fwd(num_colors), temp_bwd(num_colors)
        
        ! Initialize output
        chi_out = cmplx(0.0_dp, 0.0_dp, dp)
        
        ! Mass term
        do site = 1, total_sites
            chi_out(site, :) = fermion_mass * chi_in(site, :)
        end do
        
        ! Hopping term with staggered phases
        do site = 1, total_sites
            do mu = 1, spacetime_dim
                ! Calculate staggered phase eta_mu(x)
                eta_mu = staggered_phase(site, mu)
                
                ! Forward hopping: eta_mu(x) * U_mu(x) * chi(x+mu)
                neighbor_fwd = lattice(site)%neighbors(mu, 1)
                temp_fwd = matmul(gauge_links(site, mu)%elements, chi_in(neighbor_fwd, :))
                chi_out(site, :) = chi_out(site, :) + eta_mu * temp_fwd / 2.0_dp
                
                ! Backward hopping: -eta_mu(x) * U_mu^dagger(x-mu) * chi(x-mu)
                neighbor_bwd = lattice(site)%neighbors(mu, 2)
                temp_bwd = matmul(conjg(transpose(gauge_links(neighbor_bwd, mu)%elements)), &
                                 chi_in(neighbor_bwd, :))
                chi_out(site, :) = chi_out(site, :) - eta_mu * temp_bwd / 2.0_dp
            end do
        end do
        
    end subroutine apply_staggered_dirac_operator
    
    ! ==========================================================================
    ! Function: staggered_phase
    ! Description: Calculate staggered phase factor eta_mu(x)
    ! eta_mu(x) = (-1)^(x_1 + x_2 + ... + x_(mu-1))
    ! ==========================================================================
    function staggered_phase(site, mu) result(eta)
        implicit none
        integer, intent(in) :: site, mu
        real(dp) :: eta
        integer :: nu, sum_coords
        
        sum_coords = 0
        do nu = 1, mu - 1
            sum_coords = sum_coords + lattice(site)%coords(nu) - 1  ! coords start from 1
        end do
        
        ! eta_mu = (-1)^sum_coords
        if (mod(sum_coords, 2) == 0) then
            eta = 1.0_dp
        else
            eta = -1.0_dp
        end if
        
    end function staggered_phase
    
    ! ==========================================================================
    ! Subroutine: apply_staggered_improved_operator
    ! Description: Apply improved staggered Dirac operator (Fat link + Naik)
    ! This is a placeholder for more advanced staggered actions
    ! ==========================================================================
    subroutine apply_staggered_improved_operator(chi_in, chi_out)
        implicit none
        complex(dp), intent(in) :: chi_in(:,:)
        complex(dp), intent(out) :: chi_out(:,:)
        
        ! For now, use standard staggered operator
        ! Future improvements: fat links, Naik term, HISQ action, etc.
        call apply_staggered_dirac_operator(chi_in, chi_out)
        
    end subroutine apply_staggered_improved_operator
    
    ! ==========================================================================
    ! Function: compute_staggered_fermion_matrix
    ! Description: Compute fermion matrix element for specific sites
    ! ==========================================================================
    function compute_staggered_fermion_matrix(site1, site2, color1, color2) result(matrix_elem)
        implicit none
        integer, intent(in) :: site1, site2, color1, color2
        complex(dp) :: matrix_elem
        integer :: mu, neighbor_fwd, neighbor_bwd
        real(dp) :: eta_mu
        
        matrix_elem = cmplx(0.0_dp, 0.0_dp, dp)
        
        ! Mass term (diagonal)
        if (site1 == site2 .and. color1 == color2) then
            matrix_elem = matrix_elem + cmplx(fermion_mass, 0.0_dp, dp)
        end if
        
        ! Hopping terms
        do mu = 1, spacetime_dim
            eta_mu = staggered_phase(site1, mu)
            neighbor_fwd = lattice(site1)%neighbors(mu, 1)
            neighbor_bwd = lattice(site1)%neighbors(mu, 2)
            
            ! Forward hopping
            if (site2 == neighbor_fwd) then
                matrix_elem = matrix_elem + &
                    eta_mu * gauge_links(site1, mu)%elements(color1, color2) / 2.0_dp
            end if
            
            ! Backward hopping
            if (site2 == neighbor_bwd) then
                matrix_elem = matrix_elem - &
                    eta_mu * conjg(gauge_links(site2, mu)%elements(color2, color1)) / 2.0_dp
            end if
        end do
        
    end function compute_staggered_fermion_matrix
    
    ! ==========================================================================
    ! Subroutine: compute_taste_matrix
    ! Description: Compute taste matrix for staggered fermions
    ! Staggered fermions have residual 4-fold taste degeneracy
    ! ==========================================================================
    subroutine compute_taste_matrix(taste_index, xi_matrix)
        implicit none
        integer, intent(in) :: taste_index
        complex(dp), intent(out) :: xi_matrix(total_sites, total_sites)
        integer :: site1, site2, mu, nu
        real(dp) :: phase_product
        
        ! Taste matrices xi_mu are constructed from products of eta phases
        xi_matrix = cmplx(0.0_dp, 0.0_dp, dp)
        
        do site1 = 1, total_sites
            do site2 = 1, total_sites
                if (site1 == site2) then
                    ! Diagonal elements
                    phase_product = 1.0_dp
                    do mu = 1, min(taste_index, spacetime_dim)
                        phase_product = phase_product * staggered_phase(site1, mu)
                    end do
                    xi_matrix(site1, site2) = cmplx(phase_product, 0.0_dp, dp)
                end if
            end do
        end do
        
    end subroutine compute_taste_matrix
    
    ! ==========================================================================
    ! Subroutine: taste_splitting_analysis
    ! Description: Analyze taste symmetry breaking (important for staggered)
    ! ==========================================================================
    subroutine taste_splitting_analysis()
        implicit none
        
        print *, 'Taste splitting analysis for staggered fermions'
        print *, 'Note: Staggered fermions have 4-fold taste degeneracy'
        print *, 'Taste splitting can be reduced with improved actions (e.g., ASQTAD, HISQ)'
        
    end subroutine taste_splitting_analysis
    
    ! ==========================================================================
    ! Subroutine: cleanup_staggered_fermions
    ! Description: Cleanup staggered fermion memory
    ! ==========================================================================
    subroutine cleanup_staggered_fermions()
        implicit none
        
        if (allocated(chi)) deallocate(chi)
        if (allocated(chi_bar)) deallocate(chi_bar)
        
    end subroutine cleanup_staggered_fermions

end module staggered_fermion_mod
