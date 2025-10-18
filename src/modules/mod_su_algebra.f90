!***********************************************************************
! Module: mod_su_algebra
! Purpose: Implement SU(Nc) group algebra and operations
!          Including generators, structure constants, group operations
!***********************************************************************
module mod_su_algebra
  use mod_parameters, only: dp, n_colors
  implicit none
  
  ! Constants
  complex(dp), parameter :: zi = (0.0_dp, 1.0_dp)  ! imaginary unit
  real(dp), parameter :: pi = 3.141592653589793_dp
  
  ! Generators (Gell-Mann matrices for SU(Nc))
  complex(dp), allocatable :: generators(:,:,:)  ! (Nc, Nc, Nc²-1)
  real(dp), allocatable :: structure_const(:,:,:) ! f^abc structure constants
  
contains

  !=====================================================================
  ! Subroutine: initialize_su_algebra
  ! Purpose: Initialize SU(Nc) algebra generators and structure constants
  !=====================================================================
  subroutine initialize_su_algebra()
    implicit none
    integer :: n_gen
    
    n_gen = n_colors**2 - 1
    
    if (allocated(generators)) deallocate(generators)
    if (allocated(structure_const)) deallocate(structure_const)
    
    allocate(generators(n_colors, n_colors, n_gen))
    allocate(structure_const(n_gen, n_gen, n_gen))
    
    generators = (0.0_dp, 0.0_dp)
    structure_const = 0.0_dp
    
    ! Generate SU(Nc) generators
    if (n_colors == 2) then
      call generate_su2_matrices()
    else if (n_colors == 3) then
      call generate_su3_matrices()
    else
      call generate_sun_matrices()
    end if
    
    ! Calculate structure constants
    call calculate_structure_constants()
    
    write(*,'(A,I2,A)') ' Initialized SU(', n_colors, ') algebra'
    
  end subroutine initialize_su_algebra
  
  !=====================================================================
  ! Subroutine: generate_su2_matrices
  ! Purpose: Generate Pauli matrices (generators of SU(2))
  !=====================================================================
  subroutine generate_su2_matrices()
    implicit none
    
    ! σ₁ / 2
    generators(1,2,1) = 0.5_dp
    generators(2,1,1) = 0.5_dp
    
    ! σ₂ / 2
    generators(1,2,2) = -0.5_dp * zi
    generators(2,1,2) = 0.5_dp * zi
    
    ! σ₃ / 2
    generators(1,1,3) = 0.5_dp
    generators(2,2,3) = -0.5_dp
    
  end subroutine generate_su2_matrices
  
  !=====================================================================
  ! Subroutine: generate_su3_matrices
  ! Purpose: Generate Gell-Mann matrices (generators of SU(3))
  !=====================================================================
  subroutine generate_su3_matrices()
    implicit none
    real(dp) :: norm
    
    ! λ₁
    generators(1,2,1) = 0.5_dp
    generators(2,1,1) = 0.5_dp
    
    ! λ₂
    generators(1,2,2) = -0.5_dp * zi
    generators(2,1,2) = 0.5_dp * zi
    
    ! λ₃
    generators(1,1,3) = 0.5_dp
    generators(2,2,3) = -0.5_dp
    
    ! λ₄
    generators(1,3,4) = 0.5_dp
    generators(3,1,4) = 0.5_dp
    
    ! λ₅
    generators(1,3,5) = -0.5_dp * zi
    generators(3,1,5) = 0.5_dp * zi
    
    ! λ₆
    generators(2,3,6) = 0.5_dp
    generators(3,2,6) = 0.5_dp
    
    ! λ₇
    generators(2,3,7) = -0.5_dp * zi
    generators(3,2,7) = 0.5_dp * zi
    
    ! λ₈
    norm = 1.0_dp / (2.0_dp * sqrt(3.0_dp))
    generators(1,1,8) = norm
    generators(2,2,8) = norm
    generators(3,3,8) = -2.0_dp * norm
    
  end subroutine generate_su3_matrices
  
  !=====================================================================
  ! Subroutine: generate_sun_matrices
  ! Purpose: Generate generators for general SU(N)
  !=====================================================================
  subroutine generate_sun_matrices()
    implicit none
    integer :: i, j, k, idx
    real(dp) :: norm
    
    idx = 0
    
    ! Symmetric generators (i < j)
    do i = 1, n_colors
      do j = i+1, n_colors
        idx = idx + 1
        generators(i,j,idx) = 0.5_dp
        generators(j,i,idx) = 0.5_dp
      end do
    end do
    
    ! Antisymmetric generators (i < j)
    do i = 1, n_colors
      do j = i+1, n_colors
        idx = idx + 1
        generators(i,j,idx) = -0.5_dp * zi
        generators(j,i,idx) = 0.5_dp * zi
      end do
    end do
    
    ! Diagonal generators
    do k = 1, n_colors-1
      idx = idx + 1
      norm = sqrt(2.0_dp / (k * (k + 1.0_dp)))
      do i = 1, k
        generators(i,i,idx) = 0.5_dp * norm
      end do
      generators(k+1,k+1,idx) = -0.5_dp * k * norm
    end do
    
  end subroutine generate_sun_matrices
  
  !=====================================================================
  ! Subroutine: calculate_structure_constants
  ! Purpose: Calculate structure constants f^abc from [T^a, T^b] = i f^abc T^c
  !=====================================================================
  subroutine calculate_structure_constants()
    implicit none
    integer :: a, b, c, i, j, k
    integer :: n_gen
    complex(dp) :: commutator(n_colors, n_colors)
    complex(dp) :: trace_val
    
    n_gen = n_colors**2 - 1
    structure_const = 0.0_dp
    
    do a = 1, n_gen
      do b = 1, n_gen
        ! Compute [T^a, T^b]
        commutator = matmul(generators(:,:,a), generators(:,:,b)) - &
                     matmul(generators(:,:,b), generators(:,:,a))
        
        ! Project onto generators: Tr([T^a,T^b] T^c) = i/2 f^abc
        do c = 1, n_gen
          trace_val = (0.0_dp, 0.0_dp)
          do i = 1, n_colors
            do j = 1, n_colors
              trace_val = trace_val + commutator(i,j) * conjg(generators(j,i,c))
            end do
          end do
          ! f^abc = -2i Tr([T^a,T^b] T^c)
          structure_const(a,b,c) = real(-2.0_dp * zi * trace_val, dp)
        end do
      end do
    end do
    
  end subroutine calculate_structure_constants
  
  !=====================================================================
  ! Function: su_matrix_exp
  ! Purpose: Compute matrix exponential exp(i α T^a) for SU(Nc) element
  !=====================================================================
  function su_matrix_exp(alpha, generator_idx) result(U)
    implicit none
    real(dp), intent(in) :: alpha
    integer, intent(in) :: generator_idx
    complex(dp) :: U(n_colors, n_colors)
    complex(dp) :: M(n_colors, n_colors)
    integer :: n, i
    real(dp) :: fact
    
    ! M = i α T^a
    M = zi * alpha * generators(:,:,generator_idx)
    
    ! Taylor expansion: exp(M) = Σ M^n / n!
    U = 0.0_dp
    do i = 1, n_colors
      U(i,i) = 1.0_dp
    end do
    
    fact = 1.0_dp
    do n = 1, 20  ! Sufficient for most cases
      fact = fact / real(n, dp)
      U = U + (M**n) * fact
      if (maxval(abs(M**n)) * fact < 1.0e-12_dp) exit
    end do
    
  end function su_matrix_exp
  
  !=====================================================================
  ! Function: su_trace
  ! Purpose: Compute trace of a matrix
  !=====================================================================
  function su_trace(M) result(tr)
    implicit none
    complex(dp), intent(in) :: M(:,:)
    complex(dp) :: tr
    integer :: i
    
    tr = (0.0_dp, 0.0_dp)
    do i = 1, size(M,1)
      tr = tr + M(i,i)
    end do
    
  end function su_trace
  
  !=====================================================================
  ! Function: su_adjoint
  ! Purpose: Compute Hermitian conjugate (adjoint) of a matrix
  !=====================================================================
  function su_adjoint(M) result(Madj)
    implicit none
    complex(dp), intent(in) :: M(:,:)
    complex(dp) :: Madj(size(M,1), size(M,2))
    
    Madj = transpose(conjg(M))
    
  end function su_adjoint
  
  !=====================================================================
  ! Subroutine: cleanup_su_algebra
  ! Purpose: Deallocate arrays
  !=====================================================================
  subroutine cleanup_su_algebra()
    implicit none
    if (allocated(generators)) deallocate(generators)
    if (allocated(structure_const)) deallocate(structure_const)
  end subroutine cleanup_su_algebra

end module mod_su_algebra
