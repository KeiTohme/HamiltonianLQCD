module hamiltonian_module
    use parameter_module
    use lattice_module
    use gauge_field_module
    use fermion_module
    implicit none
    
    type hamiltonian
        ! Hamiltonian components
        real(8) :: H_total              ! Total energy
        real(8) :: H_gauge              ! Pure gauge energy
        real(8) :: H_fermion            ! Fermion kinetic energy
        real(8) :: H_interaction        ! Gauge-fermion interaction
        
        ! For efficient computation
        logical :: is_constructed
        
        ! Tensor network representation (if used)
        logical :: use_tn
        integer :: tn_bond_dim          ! Bond dimension for tensor network
        
        ! Quantum computer representation (if used)
        logical :: use_qc
        integer :: n_qubits             ! Number of qubits needed
    end type hamiltonian
    
contains
    
    subroutine construct_hamiltonian(H, gauge, fermions, lat, params)
        type(hamiltonian), intent(out) :: H
        type(gauge_field), intent(in) :: gauge
        type(fermion_field), intent(in) :: fermions
        type(lattice), intent(in) :: lat
        type(parameters), intent(in) :: params
        
        H%is_constructed = .false.
        H%use_tn = params%use_tensor_network
        H%use_qc = params%use_quantum_computer
        
        ! Initialize energies
        H%H_total = 0.0d0
        H%H_gauge = 0.0d0
        H%H_fermion = 0.0d0
        H%H_interaction = 0.0d0
        
        ! Compute gauge field Hamiltonian (electric + magnetic)
        call compute_gauge_hamiltonian(H%H_gauge, gauge, lat, params)
        
        ! Compute fermion Hamiltonian
        call compute_fermion_hamiltonian(H%H_fermion, fermions, gauge, lat, params)
        
        ! Compute interaction terms (included in fermion part for QCD)
        H%H_interaction = 0.0d0  ! Already included in fermion Hamiltonian
        
        ! Total Hamiltonian
        H%H_total = H%H_gauge + H%H_fermion + H%H_interaction
        
        H%is_constructed = .true.
        
        ! Setup tensor network representation if requested
        if (H%use_tn) then
            call setup_tensor_network_hamiltonian(H, gauge, fermions, lat, params)
        endif
        
        ! Setup quantum computer representation if requested
        if (H%use_qc) then
            call setup_quantum_computer_hamiltonian(H, gauge, fermions, lat, params)
        endif
        
    end subroutine construct_hamiltonian
    
    subroutine compute_gauge_hamiltonian(H_gauge, gauge, lat, params)
        real(8), intent(out) :: H_gauge
        type(gauge_field), intent(in) :: gauge
        type(lattice), intent(in) :: lat
        type(parameters), intent(in) :: params
        
        real(8) :: H_electric, H_magnetic
        integer :: site, link, mu, nu
        integer :: i, j
        complex(8) :: E2, tr_E2
        complex(8) :: P(params%Nc, params%Nc), tr_P
        real(8) :: g2
        
        g2 = params%g**2
        H_electric = 0.0d0
        H_magnetic = 0.0d0
        
        ! Electric field energy: H_E = (g^2/2) * sum_{x,mu} Tr[E^a_mu(x) E^a_mu(x)]
        do link = 1, gauge%n_links
            do site = 1, lat%n_sites
                ! Compute Tr[E^2]
                E2 = matmul(gauge%E(:,:,site,link), gauge%E(:,:,site,link))
                tr_E2 = (0.0d0, 0.0d0)
                do i = 1, params%Nc
                    tr_E2 = tr_E2 + E2(i,i)
                end do
                H_electric = H_electric + real(tr_E2, 8)
            end do
        end do
        H_electric = H_electric * g2 / 2.0d0
        
        ! Magnetic field energy: H_B = (2/g^2) * sum_{x,mu<nu} Re[Tr(1 - P_{mu,nu}(x))]
        ! where P_{mu,nu} is the plaquette
        do site = 1, lat%n_sites
            do mu = 1, lat%n_dim - 1
                do nu = mu + 1, lat%n_dim
                    ! Calculate plaquette
                    call plaquette(gauge, lat, site, mu, nu, params%Nc, P)
                    
                    ! Compute Re[Tr(1 - P)]
                    tr_P = (0.0d0, 0.0d0)
                    do i = 1, params%Nc
                        tr_P = tr_P + P(i,i)
                    end do
                    
                    H_magnetic = H_magnetic + params%Nc - real(tr_P, 8)
                end do
            end do
        end do
        H_magnetic = H_magnetic * 2.0d0 / g2
        
        H_gauge = H_electric + H_magnetic
        
    end subroutine compute_gauge_hamiltonian
    
    subroutine compute_fermion_hamiltonian(H_fermion, fermions, gauge, lat, params)
        real(8), intent(out) :: H_fermion
        type(fermion_field), intent(in) :: fermions
        type(gauge_field), intent(in) :: gauge
        type(lattice), intent(in) :: lat
        type(parameters), intent(in) :: params
        
        complex(8), allocatable :: D_psi(:,:,:,:)
        complex(8) :: inner_product
        integer :: site, color, spin, flavor
        
        H_fermion = 0.0d0
        
        ! Allocate temporary array for D*psi
        allocate(D_psi(params%Nc, fermions%n_spin, lat%n_sites, fermions%n_flavors))
        
        select case(trim(fermions%fermion_type))
            case("wilson")
                ! Apply Wilson Dirac operator
                call wilson_dirac_operator(D_psi, fermions%psi, gauge, fermions, lat, params, .false.)
                
                ! Compute psi^\dagger D psi
                do flavor = 1, fermions%n_flavors
                    do site = 1, lat%n_sites
                        do spin = 1, 4
                            do color = 1, params%Nc
                                inner_product = conjg(fermions%psi(color,spin,site,flavor)) * &
                                               D_psi(color,spin,site,flavor)
                                H_fermion = H_fermion + real(inner_product, 8)
                            end do
                        end do
                    end do
                end do
                
            case("staggered")
                ! Apply staggered Dirac operator
                call staggered_dirac_operator(D_psi, fermions%psi, gauge, fermions, lat, params, .false.)
                
                ! Compute chi^\dagger D chi
                do flavor = 1, fermions%n_flavors
                    do site = 1, lat%n_sites
                        do color = 1, params%Nc
                            inner_product = conjg(fermions%psi(color,1,site,flavor)) * &
                                           D_psi(color,1,site,flavor)
                            H_fermion = H_fermion + real(inner_product, 8)
                        end do
                    end do
                end do
        end select
        
        ! Add contribution from canonical momenta if needed
        ! In the Hamiltonian formalism: H_f = pi^\dagger pi + psi^\dagger D psi
        do flavor = 1, fermions%n_flavors
            do site = 1, lat%n_sites
                do spin = 1, fermions%n_spin
                    do color = 1, params%Nc
                        inner_product = conjg(fermions%pi_psi(color,spin,site,flavor)) * &
                                       fermions%pi_psi(color,spin,site,flavor)
                        H_fermion = H_fermion + real(inner_product, 8)
                    end do
                end do
            end do
        end do
        
        deallocate(D_psi)
        
    end subroutine compute_fermion_hamiltonian
    
    subroutine setup_tensor_network_hamiltonian(H, gauge, fermions, lat, params)
        type(hamiltonian), intent(inout) :: H
        type(gauge_field), intent(in) :: gauge
        type(fermion_field), intent(in) :: fermions
        type(lattice), intent(in) :: lat
        type(parameters), intent(in) :: params
        
        ! Set default bond dimension
        H%tn_bond_dim = 100  ! Can be adjusted based on accuracy requirements
        
        ! Tensor network setup would go here
        ! This is a placeholder for the actual implementation
        print *, "Setting up tensor network representation..."
        print *, "Bond dimension:", H%tn_bond_dim
        print *, "Backend:", trim(params%tn_backend)
        
    end subroutine setup_tensor_network_hamiltonian
    
    subroutine setup_quantum_computer_hamiltonian(H, gauge, fermions, lat, params)
        type(hamiltonian), intent(inout) :: H
        type(gauge_field), intent(in) :: gauge
        type(fermion_field), intent(in) :: fermions
        type(lattice), intent(in) :: lat
        type(parameters), intent(in) :: params
        
        integer :: n_gauge_qubits, n_fermion_qubits
        
        ! Estimate number of qubits needed
        ! For gauge fields: ~log2(Nc) * n_links * truncation_level
        ! For fermions: n_sites * n_colors * n_spin (for Wilson)
        
        n_gauge_qubits = ceiling(log(real(params%Nc, 8))/log(2.0d0)) * gauge%n_links * 4
        
        select case(trim(fermions%fermion_type))
            case("wilson")
                n_fermion_qubits = lat%n_sites * params%Nc * 4 * 2  ! Factor of 2 for particle/hole
            case("staggered")
                n_fermion_qubits = lat%n_sites * params%Nc * 2
        end select
        
        H%n_qubits = n_gauge_qubits + n_fermion_qubits
        
        print *, "Setting up quantum computer representation..."
        print *, "Total qubits needed:", H%n_qubits
        print *, "  Gauge qubits:", n_gauge_qubits
        print *, "  Fermion qubits:", n_fermion_qubits
        print *, "Backend:", trim(params%qc_backend)
        
    end subroutine setup_quantum_computer_hamiltonian
    
    subroutine hamiltonian_matrix_vector_product(H_psi, psi, H, gauge, fermions, lat, params)
        ! Apply Hamiltonian to a state vector
        complex(8), intent(out) :: H_psi(:)
        complex(8), intent(in) :: psi(:)
        type(hamiltonian), intent(in) :: H
        type(gauge_field), intent(inout) :: gauge
        type(fermion_field), intent(inout) :: fermions
        type(lattice), intent(in) :: lat
        type(parameters), intent(in) :: params
        
        ! This would implement the action of H on a state vector
        ! For now, this is a placeholder
        H_psi = (0.0d0, 0.0d0)
        
        print *, "Hamiltonian matrix-vector product not yet fully implemented"
        
    end subroutine hamiltonian_matrix_vector_product
    
    subroutine compute_gauss_law_constraint(G_a, gauge, fermions, lat, params, site)
        ! Compute Gauss law constraint G^a(x) = 0
        real(8), intent(out) :: G_a(:)  ! G_a(generator_index)
        type(gauge_field), intent(in) :: gauge
        type(fermion_field), intent(in) :: fermions
        type(lattice), intent(in) :: lat
        type(parameters), intent(in) :: params
        integer, intent(in) :: site
        
        integer :: mu, a, neighbor_bwd
        complex(8) :: div_E
        
        G_a = 0.0d0
        
        ! Divergence of electric field: sum_mu [E^a_mu(x) - E^a_mu(x-mu)]
        do a = 1, gauge%n_generators
            div_E = (0.0d0, 0.0d0)
            
            do mu = 1, lat%n_dim
                ! E^a_mu(x)
                div_E = div_E + trace_with_generator(gauge%E(:,:,site,mu), &
                                                    gauge%generators(:,:,a), params%Nc)
                
                ! -E^a_mu(x-mu)
                neighbor_bwd = lat%neighbors(site, 2*mu)
                div_E = div_E - trace_with_generator(gauge%E(:,:,neighbor_bwd,mu), &
                                                    gauge%generators(:,:,a), params%Nc)
            end do
            
            G_a(a) = real(div_E, 8)
        end do
        
        ! Add fermion color charge density
        ! This would require computing psi^\dagger T^a psi at site
        
    end subroutine compute_gauss_law_constraint
    
    complex(8) function trace_with_generator(M, T_a, n)
        complex(8), intent(in) :: M(:,:), T_a(:,:)
        integer, intent(in) :: n
        integer :: i
        
        trace_with_generator = (0.0d0, 0.0d0)
        do i = 1, n
            trace_with_generator = trace_with_generator + sum(M(i,:) * T_a(:,i))
        end do
        
    end function trace_with_generator
    
    subroutine cleanup_hamiltonian(H)
        type(hamiltonian), intent(inout) :: H
        
        H%is_constructed = .false.
        
    end subroutine cleanup_hamiltonian
    
end module hamiltonian_module