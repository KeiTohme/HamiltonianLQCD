! ==============================================================================
! Module: su_n_mod
! Description: SU(N) gauge group operations and representations
! ==============================================================================
module su_n_mod
    use parameters_mod
    implicit none
    
    ! SU(N) matrix type (complex)
    type :: su_n_matrix
        complex(dp), allocatable :: elements(:,:)
    end type su_n_matrix
    
    ! Gauge link variables (one for each direction at each site)
    type(su_n_matrix), allocatable :: gauge_links(:,:)  ! (site, direction)
    
    ! Gauge momenta (conjugate to gauge links) for Hamiltonian formalism
    type(su_n_matrix), allocatable :: gauge_momenta(:,:)
    
    ! SU(N) generators (Lie algebra)
    type(su_n_matrix), allocatable :: su_n_generators(:)
    
contains

    ! ==========================================================================
    ! Subroutine: initialize_gauge_fields
    ! Description: Initialize gauge link variables and momenta
    ! ==========================================================================
    subroutine initialize_gauge_fields()
        implicit none
        integer :: i, mu, a
        
        ! Allocate gauge links and momenta
        allocate(gauge_links(total_sites, spacetime_dim))
        allocate(gauge_momenta(total_sites, spacetime_dim))
        
        ! Initialize each link and momentum
        do i = 1, total_sites
            do mu = 1, spacetime_dim
                allocate(gauge_links(i,mu)%elements(num_colors, num_colors))
                allocate(gauge_momenta(i,mu)%elements(num_colors, num_colors))
                
                ! Initialize links to identity (cold start)
                call set_to_identity(gauge_links(i,mu))
                
                ! Initialize momenta to zero
                gauge_momenta(i,mu)%elements = cmplx(0.0_dp, 0.0_dp, dp)
            end do
        end do
        
        ! Initialize SU(N) generators
        call initialize_su_n_generators()
        
        print *, 'Gauge fields initialized for SU(', num_colors, ')'
        
    end subroutine initialize_gauge_fields
    
    ! ==========================================================================
    ! Subroutine: initialize_su_n_generators
    ! Description: Initialize SU(N) Lie algebra generators
    ! ==========================================================================
    subroutine initialize_su_n_generators()
        implicit none
        integer :: num_generators, a
        
        ! Number of generators for SU(N) is N^2 - 1
        num_generators = num_colors * num_colors - 1
        allocate(su_n_generators(num_generators))
        
        do a = 1, num_generators
            allocate(su_n_generators(a)%elements(num_colors, num_colors))
            su_n_generators(a)%elements = cmplx(0.0_dp, 0.0_dp, dp)
        end do
        
        ! For SU(2), use Pauli matrices
        if (num_colors == 2) then
            call initialize_su2_generators()
        ! For SU(3), use Gell-Mann matrices
        else if (num_colors == 3) then
            call initialize_su3_generators()
        else
            ! Generic SU(N) generators
            call initialize_generic_su_n_generators()
        end if
        
    end subroutine initialize_su_n_generators
    
    ! ==========================================================================
    ! Subroutine: initialize_su2_generators
    ! Description: Initialize SU(2) generators (Pauli matrices / 2)
    ! ==========================================================================
    subroutine initialize_su2_generators()
        implicit none
        complex(dp) :: i_unit
        
        i_unit = cmplx(0.0_dp, 1.0_dp, dp)
        
        ! sigma_1 / 2
        su_n_generators(1)%elements(1,1) = cmplx(0.0_dp, 0.0_dp, dp)
        su_n_generators(1)%elements(1,2) = cmplx(0.5_dp, 0.0_dp, dp)
        su_n_generators(1)%elements(2,1) = cmplx(0.5_dp, 0.0_dp, dp)
        su_n_generators(1)%elements(2,2) = cmplx(0.0_dp, 0.0_dp, dp)
        
        ! sigma_2 / 2
        su_n_generators(2)%elements(1,1) = cmplx(0.0_dp, 0.0_dp, dp)
        su_n_generators(2)%elements(1,2) = cmplx(0.0_dp, -0.5_dp, dp)
        su_n_generators(2)%elements(2,1) = cmplx(0.0_dp, 0.5_dp, dp)
        su_n_generators(2)%elements(2,2) = cmplx(0.0_dp, 0.0_dp, dp)
        
        ! sigma_3 / 2
        su_n_generators(3)%elements(1,1) = cmplx(0.5_dp, 0.0_dp, dp)
        su_n_generators(3)%elements(1,2) = cmplx(0.0_dp, 0.0_dp, dp)
        su_n_generators(3)%elements(2,1) = cmplx(0.0_dp, 0.0_dp, dp)
        su_n_generators(3)%elements(2,2) = cmplx(-0.5_dp, 0.0_dp, dp)
        
    end subroutine initialize_su2_generators
    
    ! ==========================================================================
    ! Subroutine: initialize_su3_generators
    ! Description: Initialize SU(3) generators (Gell-Mann matrices / 2)
    ! ==========================================================================
    subroutine initialize_su3_generators()
        implicit none
        complex(dp) :: i_unit
        real(dp) :: sqrt3_inv
        
        i_unit = cmplx(0.0_dp, 1.0_dp, dp)
        sqrt3_inv = 1.0_dp / sqrt(3.0_dp)
        
        ! Lambda_1
        su_n_generators(1)%elements = cmplx(0.0_dp, 0.0_dp, dp)
        su_n_generators(1)%elements(1,2) = cmplx(0.5_dp, 0.0_dp, dp)
        su_n_generators(1)%elements(2,1) = cmplx(0.5_dp, 0.0_dp, dp)
        
        ! Lambda_2
        su_n_generators(2)%elements = cmplx(0.0_dp, 0.0_dp, dp)
        su_n_generators(2)%elements(1,2) = cmplx(0.0_dp, -0.5_dp, dp)
        su_n_generators(2)%elements(2,1) = cmplx(0.0_dp, 0.5_dp, dp)
        
        ! Lambda_3
        su_n_generators(3)%elements = cmplx(0.0_dp, 0.0_dp, dp)
        su_n_generators(3)%elements(1,1) = cmplx(0.5_dp, 0.0_dp, dp)
        su_n_generators(3)%elements(2,2) = cmplx(-0.5_dp, 0.0_dp, dp)
        
        ! Lambda_4
        su_n_generators(4)%elements = cmplx(0.0_dp, 0.0_dp, dp)
        su_n_generators(4)%elements(1,3) = cmplx(0.5_dp, 0.0_dp, dp)
        su_n_generators(4)%elements(3,1) = cmplx(0.5_dp, 0.0_dp, dp)
        
        ! Lambda_5
        su_n_generators(5)%elements = cmplx(0.0_dp, 0.0_dp, dp)
        su_n_generators(5)%elements(1,3) = cmplx(0.0_dp, -0.5_dp, dp)
        su_n_generators(5)%elements(3,1) = cmplx(0.0_dp, 0.5_dp, dp)
        
        ! Lambda_6
        su_n_generators(6)%elements = cmplx(0.0_dp, 0.0_dp, dp)
        su_n_generators(6)%elements(2,3) = cmplx(0.5_dp, 0.0_dp, dp)
        su_n_generators(6)%elements(3,2) = cmplx(0.5_dp, 0.0_dp, dp)
        
        ! Lambda_7
        su_n_generators(7)%elements = cmplx(0.0_dp, 0.0_dp, dp)
        su_n_generators(7)%elements(2,3) = cmplx(0.0_dp, -0.5_dp, dp)
        su_n_generators(7)%elements(3,2) = cmplx(0.0_dp, 0.5_dp, dp)
        
        ! Lambda_8
        su_n_generators(8)%elements = cmplx(0.0_dp, 0.0_dp, dp)
        su_n_generators(8)%elements(1,1) = cmplx(0.5_dp * sqrt3_inv, 0.0_dp, dp)
        su_n_generators(8)%elements(2,2) = cmplx(0.5_dp * sqrt3_inv, 0.0_dp, dp)
        su_n_generators(8)%elements(3,3) = cmplx(-sqrt3_inv, 0.0_dp, dp)
        
    end subroutine initialize_su3_generators
    
    ! ==========================================================================
    ! Subroutine: initialize_generic_su_n_generators
    ! Description: Initialize generic SU(N) generators
    ! ==========================================================================
    subroutine initialize_generic_su_n_generators()
        implicit none
        integer :: a, i, j, k
        real(dp) :: norm
        
        ! Simplified version: create orthogonal traceless Hermitian matrices
        a = 1
        
        ! Symmetric generators
        do i = 1, num_colors
            do j = i+1, num_colors
                su_n_generators(a)%elements = cmplx(0.0_dp, 0.0_dp, dp)
                su_n_generators(a)%elements(i,j) = cmplx(0.5_dp, 0.0_dp, dp)
                su_n_generators(a)%elements(j,i) = cmplx(0.5_dp, 0.0_dp, dp)
                a = a + 1
            end do
        end do
        
        ! Antisymmetric generators
        do i = 1, num_colors
            do j = i+1, num_colors
                su_n_generators(a)%elements = cmplx(0.0_dp, 0.0_dp, dp)
                su_n_generators(a)%elements(i,j) = cmplx(0.0_dp, -0.5_dp, dp)
                su_n_generators(a)%elements(j,i) = cmplx(0.0_dp, 0.5_dp, dp)
                a = a + 1
            end do
        end do
        
        ! Diagonal generators (traceless)
        do k = 1, num_colors - 1
            norm = sqrt(2.0_dp / real(k*(k+1), dp))
            su_n_generators(a)%elements = cmplx(0.0_dp, 0.0_dp, dp)
            do i = 1, k
                su_n_generators(a)%elements(i,i) = cmplx(0.5_dp * norm, 0.0_dp, dp)
            end do
            su_n_generators(a)%elements(k+1,k+1) = cmplx(-0.5_dp * k * norm, 0.0_dp, dp)
            a = a + 1
        end do
        
    end subroutine initialize_generic_su_n_generators
    
    ! ==========================================================================
    ! Subroutine: set_to_identity
    ! Description: Set SU(N) matrix to identity
    ! ==========================================================================
    subroutine set_to_identity(matrix)
        implicit none
        type(su_n_matrix), intent(inout) :: matrix
        integer :: i
        
        matrix%elements = cmplx(0.0_dp, 0.0_dp, dp)
        do i = 1, num_colors
            matrix%elements(i,i) = cmplx(1.0_dp, 0.0_dp, dp)
        end do
        
    end subroutine set_to_identity
    
    ! ==========================================================================
    ! Function: matrix_multiply
    ! Description: Multiply two SU(N) matrices
    ! ==========================================================================
    function matrix_multiply(A, B) result(C)
        implicit none
        type(su_n_matrix), intent(in) :: A, B
        type(su_n_matrix) :: C
        
        allocate(C%elements(num_colors, num_colors))
        C%elements = matmul(A%elements, B%elements)
        
    end function matrix_multiply
    
    ! ==========================================================================
    ! Function: matrix_dagger
    ! Description: Hermitian conjugate of SU(N) matrix
    ! ==========================================================================
    function matrix_dagger(A) result(A_dag)
        implicit none
        type(su_n_matrix), intent(in) :: A
        type(su_n_matrix) :: A_dag
        
        allocate(A_dag%elements(num_colors, num_colors))
        A_dag%elements = conjg(transpose(A%elements))
        
    end function matrix_dagger
    
    ! ==========================================================================
    ! Function: trace
    ! Description: Trace of SU(N) matrix
    ! ==========================================================================
    function trace(A) result(tr)
        implicit none
        type(su_n_matrix), intent(in) :: A
        complex(dp) :: tr
        integer :: i
        
        tr = cmplx(0.0_dp, 0.0_dp, dp)
        do i = 1, num_colors
            tr = tr + A%elements(i,i)
        end do
        
    end function trace
    
    ! ==========================================================================
    ! Subroutine: cleanup_gauge_fields
    ! Description: Cleanup gauge field memory
    ! ==========================================================================
    subroutine cleanup_gauge_fields()
        implicit none
        integer :: i, mu, a
        
        if (allocated(gauge_links)) then
            do i = 1, total_sites
                do mu = 1, spacetime_dim
                    if (allocated(gauge_links(i,mu)%elements)) &
                        deallocate(gauge_links(i,mu)%elements)
                end do
            end do
            deallocate(gauge_links)
        end if
        
        if (allocated(gauge_momenta)) then
            do i = 1, total_sites
                do mu = 1, spacetime_dim
                    if (allocated(gauge_momenta(i,mu)%elements)) &
                        deallocate(gauge_momenta(i,mu)%elements)
                end do
            end do
            deallocate(gauge_momenta)
        end if
        
        if (allocated(su_n_generators)) then
            do a = 1, size(su_n_generators)
                if (allocated(su_n_generators(a)%elements)) &
                    deallocate(su_n_generators(a)%elements)
            end do
            deallocate(su_n_generators)
        end if
        
    end subroutine cleanup_gauge_fields

end module su_n_mod
