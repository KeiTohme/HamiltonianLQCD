module fermion_module
    use parameter_module
    use lattice_module
    use gauge_field_module
    implicit none
    
    type fermion_field
        character(len=20) :: fermion_type  ! "wilson" or "staggered"
        integer :: n_sites
        integer :: n_colors
        integer :: n_flavors
        integer :: n_spin                  ! 4 for Wilson, 1 for staggered
        
        ! Fermion fields: psi(color, spin, site, flavor)
        complex(8), allocatable :: psi(:,:,:,:)
        complex(8), allocatable :: psi_bar(:,:,:,:)
        
        ! For Hamiltonian formalism - canonical momentum
        complex(8), allocatable :: pi_psi(:,:,:,:)
        complex(8), allocatable :: pi_psi_bar(:,:,:,:)
        
        ! Mass parameters
        real(8), allocatable :: mass(:)    ! Mass for each flavor
        
        ! Wilson fermion specific
        real(8) :: r_wilson               ! Wilson parameter (usually 1)
        
        ! Staggered fermion specific
        integer, allocatable :: staggered_phase(:,:)  ! eta(x,mu)
    end type fermion_field
    
    ! Gamma matrices for Dirac algebra
    complex(8), dimension(4,4,5) :: gamma_matrices  ! gamma_0 to gamma_3, plus gamma_5
    
contains
    
    subroutine initialize_fermion_field(fermions, lat, params)
        type(fermion_field), intent(out) :: fermions
        type(lattice), intent(in) :: lat
        type(parameters), intent(in) :: params
        
        fermions%fermion_type = params%fermion_type
        fermions%n_sites = lat%n_sites
        fermions%n_colors = params%Nc
        fermions%n_flavors = 1  ! Can be extended
        
        ! Initialize gamma matrices
        call initialize_gamma_matrices()
        
        select case(trim(params%fermion_type))
            case("wilson")
                call initialize_wilson_fermion(fermions, lat, params)
            case("staggered")
                call initialize_staggered_fermion(fermions, lat, params)
            case default
                print *, "Error: Unknown fermion type"
                stop
        end select
        
    end subroutine initialize_fermion_field
    
    subroutine initialize_gamma_matrices()
        ! Initialize gamma matrices in the Dirac representation
        integer :: i, j
        
        ! Initialize all to zero
        gamma_matrices = (0.0d0, 0.0d0)
        
        ! gamma^0 (timelike)
        do i = 1, 2
            gamma_matrices(i,i,1) = (1.0d0, 0.0d0)
            gamma_matrices(i+2,i+2,1) = (-1.0d0, 0.0d0)
        end do
        
        ! gamma^1 (Pauli sigma_1 structure)
        gamma_matrices(1,4,2) = (1.0d0, 0.0d0)
        gamma_matrices(2,3,2) = (1.0d0, 0.0d0)
        gamma_matrices(3,2,2) = (-1.0d0, 0.0d0)
        gamma_matrices(4,1,2) = (-1.0d0, 0.0d0)
        
        ! gamma^2 (Pauli sigma_2 structure)
        gamma_matrices(1,4,3) = (0.0d0, -1.0d0)
        gamma_matrices(2,3,3) = (0.0d0, 1.0d0)
        gamma_matrices(3,2,3) = (0.0d0, 1.0d0)
        gamma_matrices(4,1,3) = (0.0d0, -1.0d0)
        
        ! gamma^3 (Pauli sigma_3 structure)
        gamma_matrices(1,3,4) = (1.0d0, 0.0d0)
        gamma_matrices(2,4,4) = (-1.0d0, 0.0d0)
        gamma_matrices(3,1,4) = (-1.0d0, 0.0d0)
        gamma_matrices(4,2,4) = (1.0d0, 0.0d0)
        
        ! gamma^5 = i * gamma^0 * gamma^1 * gamma^2 * gamma^3
        do i = 1, 2
            gamma_matrices(i,i+2,5) = (1.0d0, 0.0d0)
            gamma_matrices(i+2,i,5) = (1.0d0, 0.0d0)
        end do
        
    end subroutine initialize_gamma_matrices
    
    subroutine initialize_wilson_fermion(fermions, lat, params)
        type(fermion_field), intent(inout) :: fermions
        type(lattice), intent(in) :: lat
        type(parameters), intent(in) :: params
        integer :: site, color, spin, flavor
        
        fermions%n_spin = 4  ! Dirac spinor
        fermions%r_wilson = 1.0d0  ! Standard Wilson parameter
        
        ! Allocate fermion fields
        allocate(fermions%psi(params%Nc, 4, lat%n_sites, fermions%n_flavors))
        allocate(fermions%psi_bar(params%Nc, 4, lat%n_sites, fermions%n_flavors))
        allocate(fermions%pi_psi(params%Nc, 4, lat%n_sites, fermions%n_flavors))
        allocate(fermions%pi_psi_bar(params%Nc, 4, lat%n_sites, fermions%n_flavors))
        allocate(fermions%mass(fermions%n_flavors))
        
        ! Set mass
        fermions%mass(1) = params%m_fermion
        
        ! Initialize fermion fields to zero (vacuum)
        fermions%psi = (0.0d0, 0.0d0)
        fermions%psi_bar = (0.0d0, 0.0d0)
        fermions%pi_psi = (0.0d0, 0.0d0)
        fermions%pi_psi_bar = (0.0d0, 0.0d0)
        
        ! Option: Initialize with random small values
        if (.false.) then
            call random_fermion_field(fermions)
        endif
        
    end subroutine initialize_wilson_fermion
    
    subroutine initialize_staggered_fermion(fermions, lat, params)
        type(fermion_field), intent(inout) :: fermions
        type(lattice), intent(in) :: lat
        type(parameters), intent(in) :: params
        integer :: site, mu, i
        integer :: coords(params%D_euclid)
        
        fermions%n_spin = 1  ! Staggered fermions have no explicit spin
        
        ! Allocate fermion fields
        allocate(fermions%psi(params%Nc, 1, lat%n_sites, fermions%n_flavors))
        allocate(fermions%psi_bar(params%Nc, 1, lat%n_sites, fermions%n_flavors))
        allocate(fermions%pi_psi(params%Nc, 1, lat%n_sites, fermions%n_flavors))
        allocate(fermions%pi_psi_bar(params%Nc, 1, lat%n_sites, fermions%n_flavors))
        allocate(fermions%mass(fermions%n_flavors))
        allocate(fermions%staggered_phase(lat%n_sites, params%D_euclid))
        
        ! Set mass
        fermions%mass(1) = params%m_fermion
        
        ! Initialize fermion fields to zero
        fermions%psi = (0.0d0, 0.0d0)
        fermions%psi_bar = (0.0d0, 0.0d0)
        fermions%pi_psi = (0.0d0, 0.0d0)
        fermions%pi_psi_bar = (0.0d0, 0.0d0)
        
        ! Calculate staggered phases eta(x,mu) = (-1)^(x_1+...+x_{mu-1})
        do site = 1, lat%n_sites
            coords = lat%site_coords(site, :)
            do mu = 1, params%D_euclid
                fermions%staggered_phase(site, mu) = 1
                do i = 1, mu-1
                    if (mod(coords(i), 2) == 1) then
                        fermions%staggered_phase(site, mu) = -fermions%staggered_phase(site, mu)
                    endif
                end do
            end do
        end do
        
    end subroutine initialize_staggered_fermion
    
    subroutine wilson_dirac_operator(psi_out, psi_in, gauge, fermions, lat, params, dag)
        complex(8), intent(out) :: psi_out(:,:,:,:)
        complex(8), intent(in) :: psi_in(:,:,:,:)
        type(gauge_field), intent(in) :: gauge
        type(fermion_field), intent(in) :: fermions
        type(lattice), intent(in) :: lat
        type(parameters), intent(in) :: params
        logical, intent(in) :: dag  ! If true, apply D^\dagger
        
        integer :: site, mu, neighbor_fwd, neighbor_bwd
        integer :: color, spin, flavor
        complex(8) :: U_mu(params%Nc, params%Nc)
        complex(8) :: term(params%Nc, 4)
        real(8) :: kappa, sign
        
        ! Wilson fermion parameters
        kappa = 1.0d0 / (2.0d0 * (fermions%mass(1) + fermions%r_wilson * params%D_euclid))
        
        ! Initialize output to zero
        psi_out = (0.0d0, 0.0d0)
        
        do flavor = 1, fermions%n_flavors
            do site = 1, lat%n_sites
                ! Diagonal term (mass + r-term)
                psi_out(:,:,site,flavor) = psi_in(:,:,site,flavor)
                
                ! Hopping terms
                do mu = 1, params%D_euclid
                    neighbor_fwd = lat%neighbors(site, 2*mu-1)
                    neighbor_bwd = lat%neighbors(site, 2*mu)
                    
                    ! Forward hopping: -kappa * sum_mu (1-gamma_mu) U_mu(x) psi(x+mu)
                    U_mu = gauge%U(:,:,site,mu)
                    term = matmul(U_mu, psi_in(:,:,neighbor_fwd,flavor))
                    
                    if (dag) then
                        ! For D^\dagger, use (1+gamma_mu)
                        call apply_one_plus_gamma(term, mu, params%Nc)
                    else
                        ! For D, use (1-gamma_mu)
                        call apply_one_minus_gamma(term, mu, params%Nc)
                    endif
                    
                    psi_out(:,:,site,flavor) = psi_out(:,:,site,flavor) - kappa * term
                    
                    ! Backward hopping: -kappa * sum_mu (1+gamma_mu) U_mu^\dagger(x-mu) psi(x-mu)
                    U_mu = conjg(transpose(gauge%U(:,:,neighbor_bwd,mu)))
                    term = matmul(U_mu, psi_in(:,:,neighbor_bwd,flavor))
                    
                    if (dag) then
                        ! For D^\dagger, use (1-gamma_mu)
                        call apply_one_minus_gamma(term, mu, params%Nc)
                    else
                        ! For D, use (1+gamma_mu)
                        call apply_one_plus_gamma(term, mu, params%Nc)
                    endif
                    
                    psi_out(:,:,site,flavor) = psi_out(:,:,site,flavor) - kappa * term
                end do
            end do
        end do
        
    end subroutine wilson_dirac_operator
    
    subroutine staggered_dirac_operator(chi_out, chi_in, gauge, fermions, lat, params, dag)
        complex(8), intent(out) :: chi_out(:,:,:,:)
        complex(8), intent(in) :: chi_in(:,:,:,:)
        type(gauge_field), intent(in) :: gauge
        type(fermion_field), intent(in) :: fermions
        type(lattice), intent(in) :: lat
        type(parameters), intent(in) :: params
        logical, intent(in) :: dag
        
        integer :: site, mu, neighbor_fwd, neighbor_bwd, flavor
        complex(8) :: U_mu(params%Nc, params%Nc)
        real(8) :: eta_mu
        
        ! Initialize output
        chi_out = (0.0d0, 0.0d0)
        
        do flavor = 1, fermions%n_flavors
            do site = 1, lat%n_sites
                ! Mass term
                chi_out(:,1,site,flavor) = 2.0d0 * fermions%mass(flavor) * chi_in(:,1,site,flavor)
                
                ! Hopping terms with staggered phases
                do mu = 1, params%D_euclid
                    neighbor_fwd = lat%neighbors(site, 2*mu-1)
                    neighbor_bwd = lat%neighbors(site, 2*mu)
                    eta_mu = real(fermions%staggered_phase(site, mu), 8)
                    
                    ! Forward hopping: eta_mu(x) U_mu(x) chi(x+mu)
                    U_mu = gauge%U(:,:,site,mu)
                    chi_out(:,1,site,flavor) = chi_out(:,1,site,flavor) + &
                        eta_mu * matmul(U_mu, chi_in(:,1,neighbor_fwd,flavor))
                    
                    ! Backward hopping: -eta_mu(x-mu) U_mu^\dagger(x-mu) chi(x-mu)
                    U_mu = conjg(transpose(gauge%U(:,:,neighbor_bwd,mu)))
                    eta_mu = real(fermions%staggered_phase(neighbor_bwd, mu), 8)
                    chi_out(:,1,site,flavor) = chi_out(:,1,site,flavor) - &
                        eta_mu * matmul(U_mu, chi_in(:,1,neighbor_bwd,flavor))
                end do
            end do
        end do
        
        ! For staggered, D^\dagger = D* (just complex conjugate)
        if (dag) then
            chi_out = conjg(chi_out)
        endif
        
    end subroutine staggered_dirac_operator
    
    subroutine apply_one_minus_gamma(psi, mu, Nc)
        complex(8), intent(inout) :: psi(:,:)  ! psi(color, spin)
        integer, intent(in) :: mu, Nc
        complex(8) :: temp(Nc, 4)
        integer :: spin
        
        temp = psi
        
        ! Apply (1 - gamma_mu)
        do spin = 1, 4
            psi(:,spin) = temp(:,spin) - matmul(temp(:,:), gamma_matrices(:,spin,mu+1))
        end do
        
    end subroutine apply_one_minus_gamma
    
    subroutine apply_one_plus_gamma(psi, mu, Nc)
        complex(8), intent(inout) :: psi(:,:)  ! psi(color, spin)
        integer, intent(in) :: mu, Nc
        complex(8) :: temp(Nc, 4)
        integer :: spin
        
        temp = psi
        
        ! Apply (1 + gamma_mu)
        do spin = 1, 4
            psi(:,spin) = temp(:,spin) + matmul(temp(:,:), gamma_matrices(:,spin,mu+1))
        end do
        
    end subroutine apply_one_plus_gamma
    
    subroutine random_fermion_field(fermions)
        type(fermion_field), intent(inout) :: fermions
        real(8) :: rand_real, rand_imag
        integer :: i, j, k, l
        
        do l = 1, fermions%n_flavors
            do k = 1, fermions%n_sites
                do j = 1, fermions%n_spin
                    do i = 1, fermions%n_colors
                        call random_number(rand_real)
                        call random_number(rand_imag)
                        fermions%psi(i,j,k,l) = cmplx(rand_real - 0.5d0, rand_imag - 0.5d0, kind=8) * 0.1d0
                        
                        call random_number(rand_real)
                        call random_number(rand_imag)
                        fermions%psi_bar(i,j,k,l) = cmplx(rand_real - 0.5d0, rand_imag - 0.5d0, kind=8) * 0.1d0
                    end do
                end do
            end do
        end do
        
    end subroutine random_fermion_field
    
    subroutine cleanup_fermion_field(fermions)
        type(fermion_field), intent(inout) :: fermions
        
        if (allocated(fermions%psi)) deallocate(fermions%psi)
        if (allocated(fermions%psi_bar)) deallocate(fermions%psi_bar)
        if (allocated(fermions%pi_psi)) deallocate(fermions%pi_psi)
        if (allocated(fermions%pi_psi_bar)) deallocate(fermions%pi_psi_bar)
        if (allocated(fermions%mass)) deallocate(fermions%mass)
        if (allocated(fermions%staggered_phase)) deallocate(fermions%staggered_phase)
        
    end subroutine cleanup_fermion_field
    
end module fermion_module