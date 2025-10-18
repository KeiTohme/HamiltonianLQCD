! ==============================================================================
! Module: hamiltonian_mod
! Description: Hamiltonian formulation for lattice gauge theory
! ==============================================================================
module hamiltonian_mod
    use parameters_mod
    use lattice_mod
    use su_n_mod
    implicit none
    
    ! Electric and magnetic field energies
    real(dp) :: electric_energy
    real(dp) :: magnetic_energy
    real(dp) :: total_energy
    
contains

    ! ==========================================================================
    ! Subroutine: compute_hamiltonian
    ! Description: Compute the gauge field Hamiltonian
    ! H = (g^2/2) * sum_x sum_i E_i^a(x)^2 - (1/g^2) * sum_x sum_plaq Re[Tr(U_plaq)]
    ! ==========================================================================
    subroutine compute_hamiltonian()
        implicit none
        
        ! Compute electric field energy (kinetic term)
        call compute_electric_energy()
        
        ! Compute magnetic field energy (plaquette term)
        call compute_magnetic_energy()
        
        ! Total Hamiltonian
        total_energy = electric_energy + magnetic_energy
        
    end subroutine compute_hamiltonian
    
    ! ==========================================================================
    ! Subroutine: compute_electric_energy
    ! Description: Compute electric field energy from gauge momenta
    ! E_i^a = P_i^a (conjugate momentum to gauge links)
    ! ==========================================================================
    subroutine compute_electric_energy()
        implicit none
        integer :: site, mu, a, b
        real(dp) :: E_squared
        complex(dp) :: E_trace
        
        electric_energy = 0.0_dp
        
        ! Sum over all sites and directions
        do site = 1, total_sites
            do mu = 1, spacetime_dim
                ! Compute Tr(E^2) where E is represented by gauge_momenta
                E_trace = cmplx(0.0_dp, 0.0_dp, dp)
                
                do a = 1, num_colors
                    do b = 1, num_colors
                        E_trace = E_trace + gauge_momenta(site, mu)%elements(a,b) * &
                                           conjg(gauge_momenta(site, mu)%elements(a,b))
                    end do
                end do
                
                electric_energy = electric_energy + real(E_trace, dp)
            end do
        end do
        
        ! Multiply by gauge coupling factor
        electric_energy = electric_energy * gauge_coupling**2 / 2.0_dp
        
    end subroutine compute_electric_energy
    
    ! ==========================================================================
    ! Subroutine: compute_magnetic_energy
    ! Description: Compute magnetic field energy from plaquettes
    ! B = (1/ig) [U_mu(x) U_nu(x+mu) U_mu^dag(x+nu) U_nu^dag(x) - h.c.]
    ! ==========================================================================
    subroutine compute_magnetic_energy()
        implicit none
        integer :: site, mu, nu
        real(dp) :: plaquette_sum
        
        magnetic_energy = 0.0_dp
        plaquette_sum = 0.0_dp
        
        ! Sum over all plaquettes
        do site = 1, total_sites
            do mu = 1, spacetime_dim - 1
                do nu = mu + 1, spacetime_dim
                    plaquette_sum = plaquette_sum + compute_plaquette(site, mu, nu)
                end do
            end do
        end do
        
        ! Magnetic energy (Wilson action)
        magnetic_energy = -plaquette_sum / gauge_coupling**2
        
    end subroutine compute_magnetic_energy
    
    ! ==========================================================================
    ! Function: compute_plaquette
    ! Description: Compute plaquette at given site and directions
    ! U_plaq = U_mu(x) U_nu(x+mu) U_mu^dag(x+nu) U_nu^dag(x)
    ! ==========================================================================
    function compute_plaquette(site, mu, nu) result(plaq_value)
        implicit none
        integer, intent(in) :: site, mu, nu
        real(dp) :: plaq_value
        integer :: site_plus_mu, site_plus_nu
        type(su_n_matrix) :: U_plaq, temp
        complex(dp) :: trace_plaq
        
        ! Get neighbor sites (check if lattice has forward/backward structure)
        if (size(lattice(site)%neighbors, 2) < 2) then
            ! Hexagonal or other non-standard lattice - use simplified plaquette
            plaq_value = 1.0_dp
            return
        end if
        
        site_plus_mu = lattice(site)%neighbors(mu, 1)
        site_plus_nu = lattice(site)%neighbors(nu, 1)
        
        ! Construct plaquette: U_mu(x) * U_nu(x+mu)
        U_plaq = matrix_multiply(gauge_links(site, mu), gauge_links(site_plus_mu, nu))
        
        ! Multiply by U_mu^dag(x+nu)
        temp = matrix_dagger(gauge_links(site_plus_nu, mu))
        U_plaq = matrix_multiply(U_plaq, temp)
        if (allocated(temp%elements)) deallocate(temp%elements)
        
        ! Multiply by U_nu^dag(x)
        temp = matrix_dagger(gauge_links(site, nu))
        U_plaq = matrix_multiply(U_plaq, temp)
        if (allocated(temp%elements)) deallocate(temp%elements)
        
        ! Take real part of trace
        trace_plaq = trace(U_plaq)
        plaq_value = real(trace_plaq, dp)
        
        ! Cleanup
        if (allocated(U_plaq%elements)) deallocate(U_plaq%elements)
        
    end function compute_plaquette
    
    ! ==========================================================================
    ! Subroutine: compute_field_strength
    ! Description: Compute gauge field strength tensor F_mu_nu
    ! ==========================================================================
    subroutine compute_field_strength(site, mu, nu, F_mu_nu)
        implicit none
        integer, intent(in) :: site, mu, nu
        type(su_n_matrix), intent(out) :: F_mu_nu
        integer :: site_plus_mu, site_plus_nu
        type(su_n_matrix) :: U_plaq, temp
        complex(dp) :: i_unit
        
        i_unit = cmplx(0.0_dp, 1.0_dp, dp)
        
        ! Get neighbor sites
        site_plus_mu = lattice(site)%neighbors(mu, 1)
        site_plus_nu = lattice(site)%neighbors(nu, 1)
        
        ! Construct plaquette
        U_plaq = matrix_multiply(gauge_links(site, mu), gauge_links(site_plus_mu, nu))
        temp = matrix_dagger(gauge_links(site_plus_nu, mu))
        U_plaq = matrix_multiply(U_plaq, temp)
        if (allocated(temp%elements)) deallocate(temp%elements)
        temp = matrix_dagger(gauge_links(site, nu))
        U_plaq = matrix_multiply(U_plaq, temp)
        if (allocated(temp%elements)) deallocate(temp%elements)
        
        ! Field strength: F = (i/g) (U_plaq - U_plaq^dag)
        allocate(F_mu_nu%elements(num_colors, num_colors))
        temp = matrix_dagger(U_plaq)
        F_mu_nu%elements = (U_plaq%elements - temp%elements) * i_unit / gauge_coupling
        
        ! Cleanup
        if (allocated(U_plaq%elements)) deallocate(U_plaq%elements)
        if (allocated(temp%elements)) deallocate(temp%elements)
        
    end subroutine compute_field_strength
    
    ! ==========================================================================
    ! Subroutine: compute_wilson_loops
    ! Description: Compute Wilson loops of various sizes for confinement study
    ! ==========================================================================
    subroutine compute_wilson_loops(R, T, wilson_loop_value)
        implicit none
        integer, intent(in) :: R, T  ! Spatial and temporal extent
        real(dp), intent(out) :: wilson_loop_value
        integer :: site, i, j
        integer :: current_site, next_site
        type(su_n_matrix) :: W_loop, temp
        complex(dp) :: trace_W
        integer :: num_loops
        
        wilson_loop_value = 0.0_dp
        num_loops = 0
        
        ! Loop over all starting positions
        do site = 1, total_sites
            ! Construct R x T Wilson loop starting from this site
            current_site = site
            
            ! Initialize to identity
            allocate(W_loop%elements(num_colors, num_colors))
            call set_to_identity(W_loop)
            
            ! Move in spatial direction R steps
            do i = 1, R
                temp = matrix_multiply(W_loop, gauge_links(current_site, 1))
                if (allocated(W_loop%elements)) deallocate(W_loop%elements)
                W_loop = temp
                current_site = lattice(current_site)%neighbors(1, 1)
            end do
            
            ! Move in temporal direction T steps
            do j = 1, T
                temp = matrix_multiply(W_loop, gauge_links(current_site, spacetime_dim))
                if (allocated(W_loop%elements)) deallocate(W_loop%elements)
                W_loop = temp
                current_site = lattice(current_site)%neighbors(spacetime_dim, 1)
            end do
            
            ! Move back in spatial direction R steps
            do i = 1, R
                current_site = lattice(current_site)%neighbors(1, 2)
                temp = matrix_dagger(gauge_links(current_site, 1))
                W_loop = matrix_multiply(W_loop, temp)
                if (allocated(temp%elements)) deallocate(temp%elements)
            end do
            
            ! Move back in temporal direction T steps
            do j = 1, T
                current_site = lattice(current_site)%neighbors(spacetime_dim, 2)
                temp = matrix_dagger(gauge_links(current_site, spacetime_dim))
                W_loop = matrix_multiply(W_loop, temp)
                if (allocated(temp%elements)) deallocate(temp%elements)
            end do
            
            ! Take trace
            trace_W = trace(W_loop)
            wilson_loop_value = wilson_loop_value + real(trace_W, dp) / num_colors
            num_loops = num_loops + 1
            
            if (allocated(W_loop%elements)) deallocate(W_loop%elements)
            
            ! For efficiency, only compute a subset of Wilson loops
            if (num_loops > 100) exit
        end do
        
        ! Average over all loops
        if (num_loops > 0) then
            wilson_loop_value = wilson_loop_value / num_loops
        end if
        
    end subroutine compute_wilson_loops
    
    ! ==========================================================================
    ! Subroutine: print_hamiltonian_info
    ! Description: Print Hamiltonian energy components
    ! ==========================================================================
    subroutine print_hamiltonian_info()
        implicit none
        
        print *, '----------------------------------------'
        print *, 'Hamiltonian Energy Components:'
        print *, '  Electric energy  = ', electric_energy
        print *, '  Magnetic energy  = ', magnetic_energy
        print *, '  Total energy     = ', total_energy
        print *, '----------------------------------------'
        
    end subroutine print_hamiltonian_info

end module hamiltonian_mod
