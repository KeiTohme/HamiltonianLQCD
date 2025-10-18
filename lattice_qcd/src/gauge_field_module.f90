module gauge_field_module
    use parameter_module
    use lattice_module
    implicit none
    
    ! Constants
    real(8), parameter :: PI = 3.141592653589793d0
    complex(8), parameter :: ci = (0.0d0, 1.0d0)
    
    type gauge_field
        integer :: Nc                   ! Number of colors
        integer :: n_sites             ! Number of lattice sites
        integer :: n_links             ! Number of links per site
        complex(8), allocatable :: U(:,:,:,:)  ! U(i,j,site,link) - link matrices
        complex(8), allocatable :: E(:,:,:,:)  ! E(i,j,site,link) - electric field
        
        ! For Hamiltonian formalism
        complex(8), allocatable :: A(:,:,:,:)  ! A(i,j,site,link) - gauge potential
        complex(8), allocatable :: Pi(:,:,:,:) ! Pi(i,j,site,link) - conjugate momentum
        
        ! SU(N) generators
        complex(8), allocatable :: generators(:,:,:)  ! T^a matrices
        integer :: n_generators         ! Number of generators (Nc^2 - 1)
    end type gauge_field
    
contains
    
    subroutine initialize_gauge_field(gauge, lat, params)
        type(gauge_field), intent(out) :: gauge
        type(lattice), intent(in) :: lat
        type(parameters), intent(in) :: params
        integer :: i, j, site, link
        
        gauge%Nc = params%Nc
        gauge%n_sites = lat%n_sites
        gauge%n_links = lat%coordination_number / 2  ! Only positive directions
        gauge%n_generators = params%Nc**2 - 1
        
        ! Allocate gauge field arrays
        allocate(gauge%U(params%Nc, params%Nc, lat%n_sites, gauge%n_links))
        allocate(gauge%E(params%Nc, params%Nc, lat%n_sites, gauge%n_links))
        allocate(gauge%A(params%Nc, params%Nc, lat%n_sites, gauge%n_links))
        allocate(gauge%Pi(params%Nc, params%Nc, lat%n_sites, gauge%n_links))
        
        ! Allocate and initialize SU(N) generators
        allocate(gauge%generators(params%Nc, params%Nc, gauge%n_generators))
        call initialize_su_n_generators(gauge%generators, params%Nc)
        
        ! Initialize gauge fields to identity (cold start)
        do link = 1, gauge%n_links
            do site = 1, lat%n_sites
                do i = 1, params%Nc
                    do j = 1, params%Nc
                        if (i == j) then
                            gauge%U(i,j,site,link) = (1.0d0, 0.0d0)
                            gauge%A(i,j,site,link) = (0.0d0, 0.0d0)
                        else
                            gauge%U(i,j,site,link) = (0.0d0, 0.0d0)
                            gauge%A(i,j,site,link) = (0.0d0, 0.0d0)
                        endif
                    end do
                end do
                ! Initialize electric field and momentum to zero
                gauge%E(:,:,site,link) = (0.0d0, 0.0d0)
                gauge%Pi(:,:,site,link) = (0.0d0, 0.0d0)
            end do
        end do
        
        ! Option for hot start or specific initialization
        if (.false.) then  ! Change to .true. for hot start
            call hot_start_gauge_field(gauge, lat, params)
        endif
        
    end subroutine initialize_gauge_field
    
    subroutine initialize_su_n_generators(generators, Nc)
        complex(8), intent(out) :: generators(:,:,:)
        integer, intent(in) :: Nc
        integer :: a, i, j, k, l, gen_idx
        real(8) :: norm
        
        ! Initialize generators based on Nc
        select case(Nc)
            case(2)
                ! Pauli matrices for SU(2)
                ! sigma_1
                generators(1,1,1) = (0.0d0, 0.0d0)
                generators(1,2,1) = (1.0d0, 0.0d0)
                generators(2,1,1) = (1.0d0, 0.0d0)
                generators(2,2,1) = (0.0d0, 0.0d0)
                
                ! sigma_2
                generators(1,1,2) = (0.0d0, 0.0d0)
                generators(1,2,2) = (0.0d0, -1.0d0)
                generators(2,1,2) = (0.0d0, 1.0d0)
                generators(2,2,2) = (0.0d0, 0.0d0)
                
                ! sigma_3
                generators(1,1,3) = (1.0d0, 0.0d0)
                generators(1,2,3) = (0.0d0, 0.0d0)
                generators(2,1,3) = (0.0d0, 0.0d0)
                generators(2,2,3) = (-1.0d0, 0.0d0)
                
                ! Normalize to T^a = sigma^a/2
                generators = generators * 0.5d0
                
            case(3)
                ! Gell-Mann matrices for SU(3)
                call initialize_gell_mann_matrices(generators)
                
            case default
                ! General SU(N) construction
                call initialize_general_su_n_generators(generators, Nc)
        end select
        
    end subroutine initialize_su_n_generators
    
    subroutine initialize_gell_mann_matrices(generators)
        complex(8), intent(out) :: generators(:,:,:)
        
        ! Lambda_1
        generators(:,:,1) = 0.0d0
        generators(1,2,1) = (1.0d0, 0.0d0)
        generators(2,1,1) = (1.0d0, 0.0d0)
        
        ! Lambda_2
        generators(:,:,2) = 0.0d0
        generators(1,2,2) = (0.0d0, -1.0d0)
        generators(2,1,2) = (0.0d0, 1.0d0)
        
        ! Lambda_3
        generators(:,:,3) = 0.0d0
        generators(1,1,3) = (1.0d0, 0.0d0)
        generators(2,2,3) = (-1.0d0, 0.0d0)
        
        ! Lambda_4
        generators(:,:,4) = 0.0d0
        generators(1,3,4) = (1.0d0, 0.0d0)
        generators(3,1,4) = (1.0d0, 0.0d0)
        
        ! Lambda_5
        generators(:,:,5) = 0.0d0
        generators(1,3,5) = (0.0d0, -1.0d0)
        generators(3,1,5) = (0.0d0, 1.0d0)
        
        ! Lambda_6
        generators(:,:,6) = 0.0d0
        generators(2,3,6) = (1.0d0, 0.0d0)
        generators(3,2,6) = (1.0d0, 0.0d0)
        
        ! Lambda_7
        generators(:,:,7) = 0.0d0
        generators(2,3,7) = (0.0d0, -1.0d0)
        generators(3,2,7) = (0.0d0, 1.0d0)
        
        ! Lambda_8
        generators(:,:,8) = 0.0d0
        generators(1,1,8) = (1.0d0, 0.0d0)
        generators(2,2,8) = (1.0d0, 0.0d0)
        generators(3,3,8) = (-2.0d0, 0.0d0)
        generators(:,:,8) = generators(:,:,8) / sqrt(3.0d0)
        
        ! Normalize
        generators = generators * 0.5d0
        
    end subroutine initialize_gell_mann_matrices
    
    subroutine initialize_general_su_n_generators(generators, Nc)
        complex(8), intent(out) :: generators(:,:,:)
        integer, intent(in) :: Nc
        integer :: a, i, j, k, l, gen_idx
        real(8) :: norm
        
        gen_idx = 0
        
        ! Off-diagonal generators (symmetric and antisymmetric)
        do i = 1, Nc
            do j = i+1, Nc
                ! Symmetric generator
                gen_idx = gen_idx + 1
                generators(:,:,gen_idx) = 0.0d0
                generators(i,j,gen_idx) = (1.0d0, 0.0d0)
                generators(j,i,gen_idx) = (1.0d0, 0.0d0)
                
                ! Antisymmetric generator
                gen_idx = gen_idx + 1
                generators(:,:,gen_idx) = 0.0d0
                generators(i,j,gen_idx) = (0.0d0, -1.0d0)
                generators(j,i,gen_idx) = (0.0d0, 1.0d0)
            end do
        end do
        
        ! Diagonal generators (Cartan subalgebra)
        do l = 1, Nc-1
            gen_idx = gen_idx + 1
            generators(:,:,gen_idx) = 0.0d0
            
            do k = 1, l
                generators(k,k,gen_idx) = (1.0d0, 0.0d0)
            end do
            generators(l+1,l+1,gen_idx) = cmplx(-l, 0.0d0, kind=8)
            
            ! Normalize
            norm = sqrt(2.0d0 / (l * (l + 1.0d0)))
            generators(:,:,gen_idx) = generators(:,:,gen_idx) * norm
        end do
        
        ! Final normalization for proper Tr(T^a T^b) = delta^{ab}/2
        generators = generators * 0.5d0
        
    end subroutine initialize_general_su_n_generators
    
    subroutine hot_start_gauge_field(gauge, lat, params)
        type(gauge_field), intent(inout) :: gauge
        type(lattice), intent(in) :: lat
        type(parameters), intent(in) :: params
        integer :: site, link, a
        real(8) :: theta(gauge%n_generators)
        real(8) :: epsilon
        complex(8) :: temp_matrix(params%Nc, params%Nc)
        
        epsilon = 0.1d0  ! Small parameter for near-identity initialization
        
        do link = 1, gauge%n_links
            do site = 1, lat%n_sites
                ! Generate random angles
                call random_number(theta)
                theta = (theta - 0.5d0) * 2.0d0 * epsilon
                
                ! Construct U = exp(i * sum_a theta^a T^a)
                temp_matrix = 0.0d0
                do a = 1, gauge%n_generators
                    temp_matrix = temp_matrix + theta(a) * gauge%generators(:,:,a)
                end do
                
                ! Exponentiate
                call matrix_exponential(ci * temp_matrix, gauge%U(:,:,site,link), params%Nc)
                
                ! Initialize A from U = exp(igA)
                gauge%A(:,:,site,link) = temp_matrix / params%g
            end do
        end do
        
    end subroutine hot_start_gauge_field
    
    subroutine matrix_exponential(A, expA, n)
        complex(8), intent(in) :: A(:,:)
        complex(8), intent(out) :: expA(:,:)
        integer, intent(in) :: n
        complex(8) :: temp(n,n), An(n,n)
        integer :: i, j, k
        real(8) :: factorial
        integer :: max_terms
        real(8) :: tolerance
        
        max_terms = 20
        tolerance = 1.0d-12
        
        ! Initialize: expA = I
        expA = 0.0d0
        do i = 1, n
            expA(i,i) = (1.0d0, 0.0d0)
        end do
        
        ! Initialize: An = I, factorial = 1
        An = 0.0d0
        do i = 1, n
            An(i,i) = (1.0d0, 0.0d0)
        end do
        factorial = 1.0d0
        
        ! Taylor series: exp(A) = I + A + A^2/2! + A^3/3! + ...
        do k = 1, max_terms
            ! An = An * A
            temp = matmul(An, A)
            An = temp
            
            factorial = factorial * k
            
            ! Add term
            expA = expA + An / factorial
            
            ! Check convergence
            if (maxval(abs(An)) / factorial < tolerance) exit
        end do
        
    end subroutine matrix_exponential
    
    subroutine plaquette(gauge, lat, site, mu, nu, Nc, P)
        type(gauge_field), intent(in) :: gauge
        type(lattice), intent(in) :: lat
        integer, intent(in) :: site, mu, nu, Nc
        complex(8), intent(out) :: P(Nc, Nc)
        integer :: site_mu, site_nu, site_mu_nu
        complex(8) :: U_mu(Nc,Nc), U_nu(Nc,Nc), U_mu_dag(Nc,Nc), U_nu_dag(Nc,Nc)
        complex(8) :: temp1(Nc,Nc), temp2(Nc,Nc)
        
        ! Get neighboring sites
        site_mu = lat%neighbors(site, 2*mu-1)      ! site + mu
        site_nu = lat%neighbors(site, 2*nu-1)      ! site + nu
        site_mu_nu = lat%neighbors(site_mu, 2*nu-1) ! site + mu + nu
        
        ! Get link matrices
        U_mu = gauge%U(:,:,site,mu)
        U_nu = gauge%U(:,:,site,nu)
        
        ! Dagger of U matrices at shifted sites
        U_mu_dag = conjg(transpose(gauge%U(:,:,site_nu,mu)))
        U_nu_dag = conjg(transpose(gauge%U(:,:,site_mu,nu)))
        
        ! Calculate plaquette: P = U_mu(x) U_nu(x+mu) U_mu^dag(x+nu) U_nu^dag(x)
        temp1 = matmul(U_mu, gauge%U(:,:,site_mu,nu))
        temp2 = matmul(temp1, U_mu_dag)
        P = matmul(temp2, U_nu_dag)
        
    end subroutine plaquette
    
    subroutine cleanup_gauge_field(gauge)
        type(gauge_field), intent(inout) :: gauge
        
        if (allocated(gauge%U)) deallocate(gauge%U)
        if (allocated(gauge%E)) deallocate(gauge%E)
        if (allocated(gauge%A)) deallocate(gauge%A)
        if (allocated(gauge%Pi)) deallocate(gauge%Pi)
        if (allocated(gauge%generators)) deallocate(gauge%generators)
        
    end subroutine cleanup_gauge_field
    
end module gauge_field_module