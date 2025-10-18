! ==============================================================================
! Module: time_evolution_mod
! Description: Real-time evolution of gauge fields (Hamiltonian formalism)
! ==============================================================================
module time_evolution_mod
    use parameters_mod
    use lattice_mod
    use su_n_mod
    use hamiltonian_mod
    implicit none
    
    ! Time evolution parameters
    integer :: current_time_step
    real(dp) :: current_time
    
    ! Observable history
    real(dp), allocatable :: energy_history(:)
    real(dp), allocatable :: plaquette_history(:)
    
contains

    ! ==========================================================================
    ! Subroutine: initialize_time_evolution
    ! Description: Initialize time evolution arrays
    ! ==========================================================================
    subroutine initialize_time_evolution()
        implicit none
        
        current_time_step = 0
        current_time = 0.0_dp
        
        ! Allocate history arrays
        allocate(energy_history(num_time_steps + 1))
        allocate(plaquette_history(num_time_steps + 1))
        
        energy_history = 0.0_dp
        plaquette_history = 0.0_dp
        
        print *, 'Time evolution initialized'
        print *, '  Time step dt    = ', time_step
        print *, '  Number of steps = ', num_time_steps
        print *, '  Total time      = ', total_time
        
    end subroutine initialize_time_evolution
    
    ! ==========================================================================
    ! Subroutine: evolve_gauge_fields
    ! Description: Evolve gauge fields in real time
    ! Uses symplectic integrator (leapfrog/Verlet method)
    ! ==========================================================================
    subroutine evolve_gauge_fields()
        implicit none
        integer :: step
        real(dp) :: avg_plaquette
        
        print *, ''
        print *, 'Starting real-time evolution of gauge fields...'
        print *, ''
        
        ! Initial measurement
        call compute_hamiltonian()
        call measure_observables(0, energy_history(1), plaquette_history(1))
        
        ! Time evolution loop
        do step = 1, num_time_steps
            current_time_step = step
            current_time = step * time_step
            
            ! Perform one time step using symplectic integrator
            call symplectic_step()
            
            ! Measure observables
            call compute_hamiltonian()
            call measure_observables(step, energy_history(step+1), plaquette_history(step+1))
            
            ! Print progress
            if (mod(step, max(1, num_time_steps/10)) == 0) then
                print '(A,I6,A,I6,A,F10.4)', '  Step ', step, ' / ', num_time_steps, &
                                              '  Time = ', current_time
                print '(A,E14.6)', '    Energy = ', total_energy
                print '(A,F10.6)', '    Avg plaquette = ', plaquette_history(step+1)
            end if
        end do
        
        print *, ''
        print *, 'Time evolution completed!'
        print *, ''
        
    end subroutine evolve_gauge_fields
    
    ! ==========================================================================
    ! Subroutine: symplectic_step
    ! Description: Perform one symplectic time step (leapfrog method)
    ! Hamilton's equations: dU/dt = g^2 * E, dE/dt = -dS/dU
    ! ==========================================================================
    subroutine symplectic_step()
        implicit none
        real(dp) :: half_dt
        
        half_dt = time_step / 2.0_dp
        
        ! Leapfrog integration:
        ! 1. Half step for momenta (E field)
        call update_momenta(half_dt)
        
        ! 2. Full step for positions (U field)
        call update_links(time_step)
        
        ! 3. Half step for momenta (E field)
        call update_momenta(half_dt)
        
    end subroutine symplectic_step
    
    ! ==========================================================================
    ! Subroutine: update_momenta
    ! Description: Update gauge momenta (electric fields)
    ! dE/dt = -dS/dU = force from plaquette action
    ! ==========================================================================
    subroutine update_momenta(dt)
        implicit none
        real(dp), intent(in) :: dt
        integer :: site, mu
        type(su_n_matrix) :: force
        
        ! Update momenta based on force from gauge action
        do site = 1, total_sites
            do mu = 1, spacetime_dim
                ! Compute force from plaquette action
                call compute_gauge_force(site, mu, force)
                
                ! Update momentum: E(t+dt) = E(t) - dt * dS/dU
                gauge_momenta(site, mu)%elements = gauge_momenta(site, mu)%elements - &
                                                   dt * force%elements
                
                if (allocated(force%elements)) deallocate(force%elements)
            end do
        end do
        
    end subroutine update_momenta
    
    ! ==========================================================================
    ! Subroutine: update_links
    ! Description: Update gauge links
    ! dU/dt = i * g^2 * E * U (matrix exponential)
    ! ==========================================================================
    subroutine update_links(dt)
        implicit none
        real(dp), intent(in) :: dt
        integer :: site, mu
        type(su_n_matrix) :: delta_U, exp_matrix
        complex(dp) :: i_unit
        
        i_unit = cmplx(0.0_dp, 1.0_dp, dp)
        
        ! Update links based on momenta (electric fields)
        do site = 1, total_sites
            do mu = 1, spacetime_dim
                ! Compute delta_U = exp(i * g^2 * dt * E) * U
                ! Simplified: use first-order approximation
                allocate(delta_U%elements(num_colors, num_colors))
                
                ! delta_U = (1 + i * g^2 * dt * E) * U
                delta_U%elements = gauge_momenta(site, mu)%elements * &
                                  i_unit * gauge_coupling**2 * dt
                
                ! Add identity
                call add_identity(delta_U)
                
                ! Multiply by current U
                delta_U = matrix_multiply(delta_U, gauge_links(site, mu))
                
                ! Update link
                gauge_links(site, mu)%elements = delta_U%elements
                
                ! Reunitarize to maintain SU(N) structure
                call reunitarize(gauge_links(site, mu))
                
                if (allocated(delta_U%elements)) deallocate(delta_U%elements)
            end do
        end do
        
    end subroutine update_links
    
    ! ==========================================================================
    ! Subroutine: compute_gauge_force
    ! Description: Compute force on gauge link from plaquette action
    ! ==========================================================================
    subroutine compute_gauge_force(site, mu, force)
        implicit none
        integer, intent(in) :: site, mu
        type(su_n_matrix), intent(out) :: force
        integer :: nu, site_minus_nu, site_plus_mu
        type(su_n_matrix) :: staple, temp1, temp2
        complex(dp) :: i_unit
        
        i_unit = cmplx(0.0_dp, 1.0_dp, dp)
        
        allocate(force%elements(num_colors, num_colors))
        force%elements = cmplx(0.0_dp, 0.0_dp, dp)
        
        allocate(staple%elements(num_colors, num_colors))
        
        ! Sum over staples in all orthogonal directions
        do nu = 1, spacetime_dim
            if (nu == mu) cycle
            
            ! Forward staple
            site_plus_mu = lattice(site)%neighbors(mu, 1)
            
            temp1 = matrix_multiply(gauge_links(site_plus_mu, nu), &
                                   matrix_dagger(gauge_links(lattice(site)%neighbors(nu, 1), mu)))
            temp2 = matrix_multiply(temp1, matrix_dagger(gauge_links(site, nu)))
            staple = temp2
            
            if (allocated(temp1%elements)) deallocate(temp1%elements)
            if (allocated(temp2%elements)) deallocate(temp2%elements)
            
            force%elements = force%elements + staple%elements
            
            ! Backward staple
            site_minus_nu = lattice(site)%neighbors(nu, 2)
            
            temp1 = matrix_dagger(gauge_links(lattice(site_minus_nu)%neighbors(mu, 1), nu))
            temp2 = matrix_multiply(temp1, matrix_dagger(gauge_links(site_minus_nu, mu)))
            temp1 = matrix_multiply(temp2, gauge_links(site_minus_nu, nu))
            
            if (allocated(temp2%elements)) deallocate(temp2%elements)
            
            force%elements = force%elements + temp1%elements
            
            if (allocated(temp1%elements)) deallocate(temp1%elements)
        end do
        
        ! Multiply by U(x,mu) and project to Lie algebra
        temp1 = matrix_multiply(gauge_links(site, mu), staple)
        call project_to_algebra(temp1, force)
        
        ! Multiply by gauge coupling factor
        force%elements = force%elements / gauge_coupling**2
        
        if (allocated(staple%elements)) deallocate(staple%elements)
        if (allocated(temp1%elements)) deallocate(temp1%elements)
        
    end subroutine compute_gauge_force
    
    ! ==========================================================================
    ! Subroutine: add_identity
    ! Description: Add identity matrix to SU(N) matrix
    ! ==========================================================================
    subroutine add_identity(matrix)
        implicit none
        type(su_n_matrix), intent(inout) :: matrix
        integer :: i
        
        do i = 1, num_colors
            matrix%elements(i,i) = matrix%elements(i,i) + cmplx(1.0_dp, 0.0_dp, dp)
        end do
        
    end subroutine add_identity
    
    ! ==========================================================================
    ! Subroutine: project_to_algebra
    ! Description: Project matrix to Lie algebra (traceless anti-Hermitian)
    ! ==========================================================================
    subroutine project_to_algebra(matrix_in, matrix_out)
        implicit none
        type(su_n_matrix), intent(in) :: matrix_in
        type(su_n_matrix), intent(inout) :: matrix_out
        complex(dp) :: tr
        integer :: i
        
        if (.not. allocated(matrix_out%elements)) then
            allocate(matrix_out%elements(num_colors, num_colors))
        end if
        
        ! Anti-Hermitian part: (M - M^dagger) / 2
        matrix_out%elements = (matrix_in%elements - &
                              conjg(transpose(matrix_in%elements))) / 2.0_dp
        
        ! Subtract trace to make traceless
        tr = cmplx(0.0_dp, 0.0_dp, dp)
        do i = 1, num_colors
            tr = tr + matrix_out%elements(i,i)
        end do
        tr = tr / num_colors
        
        do i = 1, num_colors
            matrix_out%elements(i,i) = matrix_out%elements(i,i) - tr
        end do
        
    end subroutine project_to_algebra
    
    ! ==========================================================================
    ! Subroutine: reunitarize
    ! Description: Restore unitarity of gauge link (Gram-Schmidt)
    ! ==========================================================================
    subroutine reunitarize(U)
        implicit none
        type(su_n_matrix), intent(inout) :: U
        complex(dp) :: v(num_colors), norm
        integer :: i, j, k
        complex(dp) :: U_new(num_colors, num_colors)
        
        ! Gram-Schmidt orthogonalization
        do i = 1, num_colors
            v = U%elements(:, i)
            
            ! Subtract projections onto previous columns
            do j = 1, i-1
                v = v - dot_product(U_new(:,j), U%elements(:,i)) * U_new(:,j)
            end do
            
            ! Normalize
            norm = sqrt(dot_product(v, v))
            if (abs(norm) > 1.0e-10_dp) then
                U_new(:, i) = v / norm
            else
                U_new(:, i) = v
            end if
        end do
        
        U%elements = U_new
        
    end subroutine reunitarize
    
    ! ==========================================================================
    ! Subroutine: measure_observables
    ! Description: Measure and store observables
    ! ==========================================================================
    subroutine measure_observables(step, energy, avg_plaq)
        implicit none
        integer, intent(in) :: step
        real(dp), intent(out) :: energy, avg_plaq
        integer :: site, mu, nu
        real(dp) :: plaq_sum
        integer :: num_plaquettes
        
        ! Energy
        energy = total_energy
        
        ! Average plaquette
        plaq_sum = 0.0_dp
        num_plaquettes = 0
        
        do site = 1, total_sites
            do mu = 1, spacetime_dim - 1
                do nu = mu + 1, spacetime_dim
                    plaq_sum = plaq_sum + compute_plaquette(site, mu, nu)
                    num_plaquettes = num_plaquettes + 1
                end do
            end do
        end do
        
        avg_plaq = plaq_sum / (num_plaquettes * num_colors)
        
    end subroutine measure_observables
    
    ! ==========================================================================
    ! Subroutine: save_results
    ! Description: Save time evolution results to file
    ! ==========================================================================
    subroutine save_results(filename)
        implicit none
        character(len=*), intent(in) :: filename
        integer :: unit_num, i
        
        unit_num = 20
        open(unit=unit_num, file=filename, status='replace', action='write')
        
        write(unit_num, '(A)') '# Time evolution results'
        write(unit_num, '(A)') '# Step, Time, Energy, Avg_Plaquette'
        
        do i = 1, num_time_steps + 1
            write(unit_num, '(I6,3E16.8)') i-1, (i-1)*time_step, &
                energy_history(i), plaquette_history(i)
        end do
        
        close(unit_num)
        
        print *, 'Results saved to: ', trim(filename)
        
    end subroutine save_results
    
    ! ==========================================================================
    ! Subroutine: cleanup_time_evolution
    ! Description: Cleanup time evolution memory
    ! ==========================================================================
    subroutine cleanup_time_evolution()
        implicit none
        
        if (allocated(energy_history)) deallocate(energy_history)
        if (allocated(plaquette_history)) deallocate(plaquette_history)
        
    end subroutine cleanup_time_evolution

end module time_evolution_mod
