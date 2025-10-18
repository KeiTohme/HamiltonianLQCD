! ==============================================================================
! Module: wilson_fermion_mod
! Description: Wilson fermion formulation on the lattice
! ==============================================================================
module wilson_fermion_mod
    use parameters_mod
    use lattice_mod
    use su_n_mod
    implicit none
    
    ! Wilson fermion field (complex spinor with color index)
    ! psi(site, spin, color)
    complex(dp), allocatable :: psi(:,:,:)
    complex(dp), allocatable :: psi_bar(:,:,:)
    
    ! Wilson parameter (usually r=1)
    real(dp), parameter :: wilson_r = 1.0_dp
    
contains

    ! ==========================================================================
    ! Subroutine: initialize_wilson_fermions
    ! Description: Initialize Wilson fermion fields
    ! ==========================================================================
    subroutine initialize_wilson_fermions()
        implicit none
        integer :: num_spin
        
        ! Number of spin components (4 for Dirac spinor in 4D)
        num_spin = 4
        
        ! Allocate fermion fields
        allocate(psi(total_sites, num_spin, num_colors))
        allocate(psi_bar(total_sites, num_spin, num_colors))
        
        ! Initialize to zero
        psi = cmplx(0.0_dp, 0.0_dp, dp)
        psi_bar = cmplx(0.0_dp, 0.0_dp, dp)
        
        print *, 'Wilson fermion fields initialized'
        
    end subroutine initialize_wilson_fermions
    
    ! ==========================================================================
    ! Subroutine: apply_wilson_dirac_operator
    ! Description: Apply Wilson-Dirac operator to fermion field
    ! D_wilson = m + sum_mu [(1-gamma_mu) U_mu(x) psi(x+mu) + 
    !                         (1+gamma_mu) U_mu^dagger(x-mu) psi(x-mu)] / 2
    ! ==========================================================================
    subroutine apply_wilson_dirac_operator(psi_in, psi_out)
        implicit none
        complex(dp), intent(in) :: psi_in(:,:,:)
        complex(dp), intent(out) :: psi_out(:,:,:)
        integer :: site, mu, neighbor_fwd, neighbor_bwd
        complex(dp) :: temp_spinor(4, num_colors)
        complex(dp) :: temp_fwd(4, num_colors), temp_bwd(4, num_colors)
        
        ! Initialize output
        psi_out = cmplx(0.0_dp, 0.0_dp, dp)
        
        ! Mass term
        do site = 1, total_sites
            psi_out(site, :, :) = fermion_mass * psi_in(site, :, :)
        end do
        
        ! Hopping term
        do site = 1, total_sites
            do mu = 1, min(spacetime_dim, 4)  ! Only first 4 dimensions for Dirac
                ! Forward hopping
                neighbor_fwd = lattice(site)%neighbors(mu, 1)
                temp_fwd = psi_in(neighbor_fwd, :, :)
                
                ! Apply (1 - gamma_mu) projector
                call apply_projector_minus(mu, temp_fwd)
                
                ! Apply gauge link U_mu(x)
                call apply_gauge_link(site, mu, temp_fwd)
                
                psi_out(site, :, :) = psi_out(site, :, :) + temp_fwd / 2.0_dp
                
                ! Backward hopping
                neighbor_bwd = lattice(site)%neighbors(mu, 2)
                temp_bwd = psi_in(neighbor_bwd, :, :)
                
                ! Apply (1 + gamma_mu) projector
                call apply_projector_plus(mu, temp_bwd)
                
                ! Apply gauge link U_mu^dagger(x-mu)
                call apply_gauge_link_dagger(neighbor_bwd, mu, temp_bwd)
                
                psi_out(site, :, :) = psi_out(site, :, :) + temp_bwd / 2.0_dp
            end do
        end do
        
        ! Wilson term for improved locality
        call apply_wilson_term(psi_in, psi_out)
        
    end subroutine apply_wilson_dirac_operator
    
    ! ==========================================================================
    ! Subroutine: apply_wilson_term
    ! Description: Apply Wilson term -r/2 * sum_mu [nabla_mu^2 psi]
    ! ==========================================================================
    subroutine apply_wilson_term(psi_in, psi_out)
        implicit none
        complex(dp), intent(in) :: psi_in(:,:,:)
        complex(dp), intent(inout) :: psi_out(:,:,:)
        integer :: site, mu, neighbor_fwd, neighbor_bwd
        complex(dp) :: laplacian(4, num_colors)
        
        do site = 1, total_sites
            laplacian = cmplx(0.0_dp, 0.0_dp, dp)
            
            do mu = 1, spacetime_dim
                neighbor_fwd = lattice(site)%neighbors(mu, 1)
                neighbor_bwd = lattice(site)%neighbors(mu, 2)
                
                ! Discrete laplacian: psi(x+mu) + psi(x-mu) - 2*psi(x)
                laplacian = laplacian + psi_in(neighbor_fwd, :, :) + &
                           psi_in(neighbor_bwd, :, :) - 2.0_dp * psi_in(site, :, :)
            end do
            
            psi_out(site, :, :) = psi_out(site, :, :) - wilson_r * laplacian / 2.0_dp
        end do
        
    end subroutine apply_wilson_term
    
    ! ==========================================================================
    ! Subroutine: apply_projector_minus
    ! Description: Apply (1 - gamma_mu) projector
    ! ==========================================================================
    subroutine apply_projector_minus(mu, spinor)
        implicit none
        integer, intent(in) :: mu
        complex(dp), intent(inout) :: spinor(4, num_colors)
        complex(dp) :: temp(4, num_colors)
        
        temp = spinor
        
        ! Simplified gamma matrix action (Dirac basis)
        select case(mu)
            case(1)  ! gamma_1
                spinor(1,:) = temp(1,:) - cmplx(0.0_dp, 1.0_dp, dp) * temp(4,:)
                spinor(2,:) = temp(2,:) - cmplx(0.0_dp, 1.0_dp, dp) * temp(3,:)
                spinor(3,:) = temp(3,:) + cmplx(0.0_dp, 1.0_dp, dp) * temp(2,:)
                spinor(4,:) = temp(4,:) + cmplx(0.0_dp, 1.0_dp, dp) * temp(1,:)
            case(2)  ! gamma_2
                spinor(1,:) = temp(1,:) - temp(4,:)
                spinor(2,:) = temp(2,:) + temp(3,:)
                spinor(3,:) = temp(3,:) + temp(2,:)
                spinor(4,:) = temp(4,:) - temp(1,:)
            case(3)  ! gamma_3
                spinor(1,:) = temp(1,:) - cmplx(0.0_dp, 1.0_dp, dp) * temp(3,:)
                spinor(2,:) = temp(2,:) + cmplx(0.0_dp, 1.0_dp, dp) * temp(4,:)
                spinor(3,:) = temp(3,:) + cmplx(0.0_dp, 1.0_dp, dp) * temp(1,:)
                spinor(4,:) = temp(4,:) - cmplx(0.0_dp, 1.0_dp, dp) * temp(2,:)
            case(4)  ! gamma_0 (or gamma_4)
                spinor(1,:) = temp(1,:) - temp(3,:)
                spinor(2,:) = temp(2,:) - temp(4,:)
                spinor(3,:) = temp(3,:) - temp(1,:)
                spinor(4,:) = temp(4,:) - temp(2,:)
        end select
        
    end subroutine apply_projector_minus
    
    ! ==========================================================================
    ! Subroutine: apply_projector_plus
    ! Description: Apply (1 + gamma_mu) projector
    ! ==========================================================================
    subroutine apply_projector_plus(mu, spinor)
        implicit none
        integer, intent(in) :: mu
        complex(dp), intent(inout) :: spinor(4, num_colors)
        complex(dp) :: temp(4, num_colors)
        
        temp = spinor
        
        ! Simplified gamma matrix action (Dirac basis)
        select case(mu)
            case(1)  ! gamma_1
                spinor(1,:) = temp(1,:) + cmplx(0.0_dp, 1.0_dp, dp) * temp(4,:)
                spinor(2,:) = temp(2,:) + cmplx(0.0_dp, 1.0_dp, dp) * temp(3,:)
                spinor(3,:) = temp(3,:) - cmplx(0.0_dp, 1.0_dp, dp) * temp(2,:)
                spinor(4,:) = temp(4,:) - cmplx(0.0_dp, 1.0_dp, dp) * temp(1,:)
            case(2)  ! gamma_2
                spinor(1,:) = temp(1,:) + temp(4,:)
                spinor(2,:) = temp(2,:) - temp(3,:)
                spinor(3,:) = temp(3,:) - temp(2,:)
                spinor(4,:) = temp(4,:) + temp(1,:)
            case(3)  ! gamma_3
                spinor(1,:) = temp(1,:) + cmplx(0.0_dp, 1.0_dp, dp) * temp(3,:)
                spinor(2,:) = temp(2,:) - cmplx(0.0_dp, 1.0_dp, dp) * temp(4,:)
                spinor(3,:) = temp(3,:) - cmplx(0.0_dp, 1.0_dp, dp) * temp(1,:)
                spinor(4,:) = temp(4,:) + cmplx(0.0_dp, 1.0_dp, dp) * temp(2,:)
            case(4)  ! gamma_0 (or gamma_4)
                spinor(1,:) = temp(1,:) + temp(3,:)
                spinor(2,:) = temp(2,:) + temp(4,:)
                spinor(3,:) = temp(3,:) + temp(1,:)
                spinor(4,:) = temp(4,:) + temp(2,:)
        end select
        
    end subroutine apply_projector_plus
    
    ! ==========================================================================
    ! Subroutine: apply_gauge_link
    ! Description: Apply gauge link matrix to spinor
    ! ==========================================================================
    subroutine apply_gauge_link(site, mu, spinor)
        implicit none
        integer, intent(in) :: site, mu
        complex(dp), intent(inout) :: spinor(4, num_colors)
        integer :: s
        complex(dp) :: temp(num_colors)
        
        do s = 1, 4
            temp = matmul(gauge_links(site, mu)%elements, spinor(s, :))
            spinor(s, :) = temp
        end do
        
    end subroutine apply_gauge_link
    
    ! ==========================================================================
    ! Subroutine: apply_gauge_link_dagger
    ! Description: Apply conjugate transpose of gauge link to spinor
    ! ==========================================================================
    subroutine apply_gauge_link_dagger(site, mu, spinor)
        implicit none
        integer, intent(in) :: site, mu
        complex(dp), intent(inout) :: spinor(4, num_colors)
        integer :: s
        complex(dp) :: temp(num_colors)
        complex(dp) :: U_dag(num_colors, num_colors)
        
        U_dag = conjg(transpose(gauge_links(site, mu)%elements))
        
        do s = 1, 4
            temp = matmul(U_dag, spinor(s, :))
            spinor(s, :) = temp
        end do
        
    end subroutine apply_gauge_link_dagger
    
    ! ==========================================================================
    ! Subroutine: cleanup_wilson_fermions
    ! Description: Cleanup Wilson fermion memory
    ! ==========================================================================
    subroutine cleanup_wilson_fermions()
        implicit none
        
        if (allocated(psi)) deallocate(psi)
        if (allocated(psi_bar)) deallocate(psi_bar)
        
    end subroutine cleanup_wilson_fermions

end module wilson_fermion_mod
