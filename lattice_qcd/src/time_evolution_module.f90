module time_evolution_module
    use parameter_module
    use lattice_module
    use gauge_field_module
    use fermion_module
    use hamiltonian_module
    implicit none
    
    ! Time evolution methods
    integer, parameter :: LEAPFROG = 1
    integer, parameter :: RUNGE_KUTTA_4 = 2
    integer, parameter :: SYMPLECTIC_4 = 3
    
contains
    
    subroutine time_evolution_main(H, gauge, fermions, lat, params)
        type(hamiltonian), intent(inout) :: H
        type(gauge_field), intent(inout) :: gauge
        type(fermion_field), intent(inout) :: fermions
        type(lattice), intent(in) :: lat
        type(parameters), intent(in) :: params
        
        integer :: step, measure_step
        real(8) :: time
        integer :: evolution_method
        character(len=256) :: filename
        integer :: unit_num
        
        ! Select evolution method
        evolution_method = LEAPFROG  ! Default
        
        ! Open output files
        unit_num = 20
        write(filename, '(A,A)') trim(params%output_dir), '/energy.dat'
        open(unit=unit_num, file=filename, status='replace')
        write(unit_num, '(A)') '# time, E_total, E_gauge, E_fermion'
        
        print *, "Starting time evolution..."
        print *, "Method: Leapfrog"
        print *, "Time step dt =", params%dt
        print *, "Number of steps =", params%n_steps
        
        ! Initial measurements
        time = 0.0d0
        call measure_observables(gauge, fermions, H, lat, params, time, unit_num)
        
        ! Time evolution loop
        do step = 1, params%n_steps
            time = step * params%dt
            
            ! Evolve one time step
            select case(evolution_method)
                case(LEAPFROG)
                    call leapfrog_step(gauge, fermions, lat, params)
                case(RUNGE_KUTTA_4)
                    call runge_kutta_4_step(gauge, fermions, lat, params)
                case(SYMPLECTIC_4)
                    call symplectic_4_step(gauge, fermions, lat, params)
            end select
            
            ! Measurements
            if (mod(step, params%measure_interval) == 0) then
                ! Recompute Hamiltonian
                call construct_hamiltonian(H, gauge, fermions, lat, params)
                call measure_observables(gauge, fermions, H, lat, params, time, unit_num)
                
                print *, "Step", step, "of", params%n_steps, &
                        " t =", time, " E =", H%H_total
            endif
            
            ! Save configuration
            if (params%save_configurations .and. mod(step, params%save_interval) == 0) then
                call save_configuration(gauge, fermions, lat, params, step)
            endif
        end do
        
        close(unit_num)
        print *, "Time evolution completed"
        
    end subroutine time_evolution_main
    
    subroutine leapfrog_step(gauge, fermions, lat, params)
        type(gauge_field), intent(inout) :: gauge
        type(fermion_field), intent(inout) :: fermions
        type(lattice), intent(in) :: lat
        type(parameters), intent(in) :: params
        
        ! Leapfrog integration: 
        ! 1. Update momenta by dt/2
        ! 2. Update positions by dt
        ! 3. Update momenta by dt/2
        
        ! Update gauge field momenta (electric field)
        call update_gauge_momentum(gauge, fermions, lat, params, params%dt/2.0d0)
        
        ! Update fermion momenta
        call update_fermion_momentum(fermions, gauge, lat, params, params%dt/2.0d0)
        
        ! Update gauge links
        call update_gauge_links(gauge, lat, params, params%dt)
        
        ! Update fermion fields
        call update_fermion_fields(fermions, gauge, lat, params, params%dt)
        
        ! Update momenta again
        call update_gauge_momentum(gauge, fermions, lat, params, params%dt/2.0d0)
        call update_fermion_momentum(fermions, gauge, lat, params, params%dt/2.0d0)
        
    end subroutine leapfrog_step
    
    subroutine update_gauge_momentum(gauge, fermions, lat, params, dt)
        type(gauge_field), intent(inout) :: gauge
        type(fermion_field), intent(in) :: fermions
        type(lattice), intent(in) :: lat
        type(parameters), intent(in) :: params
        real(8), intent(in) :: dt
        
        integer :: site, mu, nu, i, j
        complex(8) :: F_U(params%Nc, params%Nc)  ! Force term
        complex(8) :: staple(params%Nc, params%Nc)
        
        ! Update E_mu(x) using Hamilton's equation: dE/dt = -dH/dA
        do mu = 1, gauge%n_links
            do site = 1, lat%n_sites
                ! Compute gauge force (staple term)
                call compute_gauge_force(F_U, gauge, lat, site, mu, params)
                
                ! Add fermion contribution to force
                call add_fermion_force(F_U, fermions, gauge, lat, site, mu, params)
                
                ! Update electric field: E_mu(x) -> E_mu(x) - dt * F_U
                gauge%E(:,:,site,mu) = gauge%E(:,:,site,mu) - dt * F_U
                
                ! Project to Lie algebra (traceless anti-Hermitian)
                call project_to_algebra(gauge%E(:,:,site,mu), params%Nc)
            end do
        end do
        
    end subroutine update_gauge_momentum
    
    subroutine update_gauge_links(gauge, lat, params, dt)
        type(gauge_field), intent(inout) :: gauge
        type(lattice), intent(in) :: lat
        type(parameters), intent(in) :: params
        real(8), intent(in) :: dt
        
        integer :: site, mu
        complex(8) :: exp_E(params%Nc, params%Nc)
        
        ! Update U_mu(x) using U_mu(x) -> exp(igt*E_mu(x)) * U_mu(x)
        do mu = 1, gauge%n_links
            do site = 1, lat%n_sites
                ! Compute exp(igt*E)
                call matrix_exponential(ci * params%g * dt * gauge%E(:,:,site,mu), &
                                      exp_E, params%Nc)
                
                ! Update link
                gauge%U(:,:,site,mu) = matmul(exp_E, gauge%U(:,:,site,mu))
                
                ! Reunitarize to maintain SU(N)
                call reunitarize(gauge%U(:,:,site,mu), params%Nc)
            end do
        end do
        
    end subroutine update_gauge_links
    
    subroutine update_fermion_momentum(fermions, gauge, lat, params, dt)
        type(fermion_field), intent(inout) :: fermions
        type(gauge_field), intent(in) :: gauge
        type(lattice), intent(in) :: lat
        type(parameters), intent(in) :: params
        real(8), intent(in) :: dt
        
        complex(8), allocatable :: D_psi(:,:,:,:)
        
        allocate(D_psi(params%Nc, fermions%n_spin, lat%n_sites, fermions%n_flavors))
        
        ! Update pi_psi using dpi/dt = -dH/dpsi = -D^\dagger psi
        select case(trim(fermions%fermion_type))
            case("wilson")
                call wilson_dirac_operator(D_psi, fermions%psi, gauge, fermions, lat, params, .true.)
            case("staggered")
                call staggered_dirac_operator(D_psi, fermions%psi, gauge, fermions, lat, params, .true.)
        end select
        
        ! Update momentum
        fermions%pi_psi = fermions%pi_psi - dt * D_psi
        
        ! Similarly for psi_bar
        select case(trim(fermions%fermion_type))
            case("wilson")
                call wilson_dirac_operator(D_psi, fermions%psi_bar, gauge, fermions, lat, params, .false.)
            case("staggered")
                call staggered_dirac_operator(D_psi, fermions%psi_bar, gauge, fermions, lat, params, .false.)
        end select
        
        fermions%pi_psi_bar = fermions%pi_psi_bar - dt * D_psi
        
        deallocate(D_psi)
        
    end subroutine update_fermion_momentum
    
    subroutine update_fermion_fields(fermions, gauge, lat, params, dt)
        type(fermion_field), intent(inout) :: fermions
        type(gauge_field), intent(in) :: gauge
        type(lattice), intent(in) :: lat
        type(parameters), intent(in) :: params
        real(8), intent(in) :: dt
        
        ! Update psi using dpsi/dt = dH/dpi = pi
        fermions%psi = fermions%psi + dt * fermions%pi_psi
        fermions%psi_bar = fermions%psi_bar + dt * fermions%pi_psi_bar
        
    end subroutine update_fermion_fields
    
    subroutine compute_gauge_force(F_U, gauge, lat, site, mu, params)
        complex(8), intent(out) :: F_U(:,:)
        type(gauge_field), intent(in) :: gauge
        type(lattice), intent(in) :: lat
        integer, intent(in) :: site, mu
        type(parameters), intent(in) :: params
        
        integer :: nu
        complex(8) :: staple(params%Nc, params%Nc), total_staple(params%Nc, params%Nc)
        complex(8) :: temp(params%Nc, params%Nc)
        
        ! Initialize
        total_staple = (0.0d0, 0.0d0)
        
        ! Sum over staples
        do nu = 1, lat%n_dim
            if (nu /= mu) then
                call compute_staple(staple, gauge, lat, site, mu, nu, params%Nc)
                total_staple = total_staple + staple
            endif
        end do
        
        ! Force is proportional to U * staple^\dagger
        temp = matmul(gauge%U(:,:,site,mu), conjg(transpose(total_staple)))
        
        ! Project to algebra: F = (temp - temp^\dagger) / 2i
        F_U = (temp - conjg(transpose(temp))) / (2.0d0 * ci)
        
        ! Make traceless
        call make_traceless(F_U, params%Nc)
        
        ! Scale by coupling
        F_U = F_U * 2.0d0 / params%g**2
        
    end subroutine compute_gauge_force
    
    subroutine compute_staple(staple, gauge, lat, site, mu, nu, Nc)
        complex(8), intent(out) :: staple(:,:)
        type(gauge_field), intent(in) :: gauge
        type(lattice), intent(in) :: lat
        integer, intent(in) :: site, mu, nu, Nc
        
        integer :: site_mu, site_nu, site_mu_nu
        integer :: site_neg_nu, site_mu_neg_nu
        complex(8) :: temp1(Nc,Nc), temp2(Nc,Nc)
        
        ! Forward staple
        site_mu = lat%neighbors(site, 2*mu-1)
        site_nu = lat%neighbors(site, 2*nu-1)
        site_mu_nu = lat%neighbors(site_mu, 2*nu-1)
        
        temp1 = matmul(gauge%U(:,:,site,nu), gauge%U(:,:,site_nu,mu))
        staple = matmul(temp1, conjg(transpose(gauge%U(:,:,site_mu,nu))))
        
        ! Backward staple
        site_neg_nu = lat%neighbors(site, 2*nu)
        site_mu_neg_nu = lat%neighbors(site_mu, 2*nu)
        
        temp1 = matmul(conjg(transpose(gauge%U(:,:,site_neg_nu,nu))), &
                      gauge%U(:,:,site_neg_nu,mu))
        temp2 = matmul(temp1, gauge%U(:,:,site_mu_neg_nu,nu))
        
        staple = staple + temp2
        
    end subroutine compute_staple
    
    subroutine add_fermion_force(F_U, fermions, gauge, lat, site, mu, params)
        complex(8), intent(inout) :: F_U(:,:)
        type(fermion_field), intent(in) :: fermions
        type(gauge_field), intent(in) :: gauge
        type(lattice), intent(in) :: lat
        integer, intent(in) :: site, mu
        type(parameters), intent(in) :: params
        
        ! Add fermion contribution to gauge force
        ! This involves computing derivatives of the fermion action with respect to U
        
        ! Placeholder - full implementation would compute
        ! F_U += -i * psi^\dagger * T^a * (fermion hopping terms)
        
    end subroutine add_fermion_force
    
    subroutine project_to_algebra(M, Nc)
        complex(8), intent(inout) :: M(:,:)
        integer, intent(in) :: Nc
        
        ! Project to su(N) Lie algebra (traceless anti-Hermitian)
        call make_anti_hermitian(M, Nc)
        call make_traceless(M, Nc)
        
    end subroutine project_to_algebra
    
    subroutine make_anti_hermitian(M, Nc)
        complex(8), intent(inout) :: M(:,:)
        integer, intent(in) :: Nc
        
        M = 0.5d0 * (M - conjg(transpose(M)))
        
    end subroutine make_anti_hermitian
    
    subroutine make_traceless(M, Nc)
        complex(8), intent(inout) :: M(:,:)
        integer, intent(in) :: Nc
        complex(8) :: tr
        integer :: i
        
        tr = (0.0d0, 0.0d0)
        do i = 1, Nc
            tr = tr + M(i,i)
        end do
        
        do i = 1, Nc
            M(i,i) = M(i,i) - tr / Nc
        end do
        
    end subroutine make_traceless
    
    subroutine reunitarize(U, Nc)
        complex(8), intent(inout) :: U(:,:)
        integer, intent(in) :: Nc
        
        ! Gram-Schmidt orthogonalization to project back to SU(N)
        integer :: i, j, k
        complex(8) :: norm, dot_product
        complex(8) :: V(Nc,Nc)
        
        V = U
        
        ! Gram-Schmidt process
        do i = 1, Nc
            ! Normalize column i
            norm = (0.0d0, 0.0d0)
            do j = 1, Nc
                norm = norm + conjg(V(j,i)) * V(j,i)
            end do
            norm = sqrt(norm)
            V(:,i) = V(:,i) / norm
            
            ! Orthogonalize remaining columns
            do j = i+1, Nc
                dot_product = (0.0d0, 0.0d0)
                do k = 1, Nc
                    dot_product = dot_product + conjg(V(k,i)) * V(k,j)
                end do
                V(:,j) = V(:,j) - dot_product * V(:,i)
            end do
        end do
        
        ! Ensure det(U) = 1
        call ensure_unit_determinant(V, Nc)
        
        U = V
        
    end subroutine reunitarize
    
    subroutine ensure_unit_determinant(U, Nc)
        complex(8), intent(inout) :: U(:,:)
        integer, intent(in) :: Nc
        complex(8) :: det
        real(8) :: phase
        
        ! Calculate determinant
        call calculate_determinant(U, Nc, det)
        
        ! Extract phase
        phase = atan2(aimag(det), real(det, 8))
        
        ! Divide last column by exp(i*phase/Nc)
        U(:,Nc) = U(:,Nc) * exp(-ci * phase / Nc)
        
    end subroutine ensure_unit_determinant
    
    subroutine calculate_determinant(A, n, det)
        complex(8), intent(in) :: A(:,:)
        integer, intent(in) :: n
        complex(8), intent(out) :: det
        
        ! Simple determinant calculation for small matrices
        select case(n)
            case(1)
                det = A(1,1)
            case(2)
                det = A(1,1)*A(2,2) - A(1,2)*A(2,1)
            case(3)
                det = A(1,1)*(A(2,2)*A(3,3) - A(2,3)*A(3,2)) - &
                      A(1,2)*(A(2,1)*A(3,3) - A(2,3)*A(3,1)) + &
                      A(1,3)*(A(2,1)*A(3,2) - A(2,2)*A(3,1))
            case default
                ! For larger matrices, use LU decomposition
                det = (1.0d0, 0.0d0)  ! Placeholder
        end select
        
    end subroutine calculate_determinant
    
    subroutine runge_kutta_4_step(gauge, fermions, lat, params)
        type(gauge_field), intent(inout) :: gauge
        type(fermion_field), intent(inout) :: fermions
        type(lattice), intent(in) :: lat
        type(parameters), intent(in) :: params
        
        ! 4th order Runge-Kutta - placeholder
        print *, "RK4 method not yet implemented, using leapfrog instead"
        call leapfrog_step(gauge, fermions, lat, params)
        
    end subroutine runge_kutta_4_step
    
    subroutine symplectic_4_step(gauge, fermions, lat, params)
        type(gauge_field), intent(inout) :: gauge
        type(fermion_field), intent(inout) :: fermions
        type(lattice), intent(in) :: lat
        type(parameters), intent(in) :: params
        
        ! 4th order symplectic integrator - placeholder
        print *, "Symplectic-4 method not yet implemented, using leapfrog instead"
        call leapfrog_step(gauge, fermions, lat, params)
        
    end subroutine symplectic_4_step
    
    subroutine measure_observables(gauge, fermions, H, lat, params, time, unit_num)
        type(gauge_field), intent(in) :: gauge
        type(fermion_field), intent(in) :: fermions
        type(hamiltonian), intent(in) :: H
        type(lattice), intent(in) :: lat
        type(parameters), intent(in) :: params
        real(8), intent(in) :: time
        integer, intent(in) :: unit_num
        
        real(8) :: plaq, polyakov, chiral_condensate
        real(8) :: gauss_violation
        
        ! Measure plaquette
        call measure_plaquette(plaq, gauge, lat, params)
        
        ! Measure Polyakov loop
        call measure_polyakov_loop(polyakov, gauge, lat, params)
        
        ! Measure chiral condensate (for fermions)
        call measure_chiral_condensate(chiral_condensate, fermions, lat, params)
        
        ! Check Gauss law constraint
        call check_gauss_law(gauss_violation, gauge, fermions, lat, params)
        
        ! Write to file
        write(unit_num, '(5E16.8)') time, H%H_total, H%H_gauge, H%H_fermion, gauss_violation
        
        ! Print to screen
        if (gauss_violation > 1.0d-6) then
            print *, "Warning: Gauss law violation =", gauss_violation
        endif
        
    end subroutine measure_observables
    
    subroutine measure_plaquette(plaq, gauge, lat, params)
        real(8), intent(out) :: plaq
        type(gauge_field), intent(in) :: gauge
        type(lattice), intent(in) :: lat
        type(parameters), intent(in) :: params
        
        integer :: site, mu, nu, i
        complex(8) :: P(params%Nc, params%Nc), tr_P
        real(8) :: sum_plaq
        integer :: n_plaq
        
        sum_plaq = 0.0d0
        n_plaq = 0
        
        do site = 1, lat%n_sites
            do mu = 1, lat%n_dim - 1
                do nu = mu + 1, lat%n_dim
                    call plaquette(gauge, lat, site, mu, nu, params%Nc, P)
                    
                    tr_P = (0.0d0, 0.0d0)
                    do i = 1, params%Nc
                        tr_P = tr_P + P(i,i)
                    end do
                    
                    sum_plaq = sum_plaq + real(tr_P, 8)
                    n_plaq = n_plaq + 1
                end do
            end do
        end do
        
        plaq = sum_plaq / (n_plaq * params%Nc)
        
    end subroutine measure_plaquette
    
    subroutine measure_polyakov_loop(polyakov, gauge, lat, params)
        real(8), intent(out) :: polyakov
        type(gauge_field), intent(in) :: gauge
        type(lattice), intent(in) :: lat
        type(parameters), intent(in) :: params
        
        ! Placeholder - would compute Polyakov loop
        polyakov = 0.0d0
        
    end subroutine measure_polyakov_loop
    
    subroutine measure_chiral_condensate(condensate, fermions, lat, params)
        real(8), intent(out) :: condensate
        type(fermion_field), intent(in) :: fermions
        type(lattice), intent(in) :: lat
        type(parameters), intent(in) :: params
        
        ! Placeholder - would compute <psi_bar psi>
        condensate = 0.0d0
        
    end subroutine measure_chiral_condensate
    
    subroutine check_gauss_law(violation, gauge, fermions, lat, params)
        real(8), intent(out) :: violation
        type(gauge_field), intent(in) :: gauge
        type(fermion_field), intent(in) :: fermions
        type(lattice), intent(in) :: lat
        type(parameters), intent(in) :: params
        
        real(8), allocatable :: G_a(:)
        integer :: site, a
        real(8) :: sum_violation
        
        allocate(G_a(gauge%n_generators))
        
        sum_violation = 0.0d0
        
        do site = 1, lat%n_sites
            call compute_gauss_law_constraint(G_a, gauge, fermions, lat, params, site)
            
            do a = 1, gauge%n_generators
                sum_violation = sum_violation + G_a(a)**2
            end do
        end do
        
        violation = sqrt(sum_violation / (lat%n_sites * gauge%n_generators))
        
        deallocate(G_a)
        
    end subroutine check_gauss_law
    
    subroutine save_configuration(gauge, fermions, lat, params, step)
        type(gauge_field), intent(in) :: gauge
        type(fermion_field), intent(in) :: fermions
        type(lattice), intent(in) :: lat
        type(parameters), intent(in) :: params
        integer, intent(in) :: step
        
        character(len=256) :: filename
        integer :: unit_num
        
        unit_num = 30
        
        ! Save gauge configuration
        write(filename, '(A,A,I6.6,A)') trim(params%output_dir), '/gauge_config_', step, '.dat'
        open(unit=unit_num, file=filename, form='unformatted', status='replace')
        write(unit_num) gauge%U
        write(unit_num) gauge%E
        close(unit_num)
        
        ! Save fermion configuration
        write(filename, '(A,A,I6.6,A)') trim(params%output_dir), '/fermion_config_', step, '.dat'
        open(unit=unit_num, file=filename, form='unformatted', status='replace')
        write(unit_num) fermions%psi
        write(unit_num) fermions%pi_psi
        close(unit_num)
        
    end subroutine save_configuration
    
end module time_evolution_module