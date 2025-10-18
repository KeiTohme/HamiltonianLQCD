!***********************************************************************
! Module: mod_gauge_field
! Purpose: Implement gauge field and link variables for Lattice QCD
!          Including electric and magnetic field operators
!***********************************************************************
module mod_gauge_field
  use mod_parameters
  use mod_lattice
  use mod_su_algebra
  implicit none
  
  ! Gauge field configuration
  complex(dp), allocatable :: link_matrices(:,:,:)    ! U_μ(x): (Nc, Nc, n_links)
  complex(dp), allocatable :: electric_field(:,:,:)   ! E^a_μ(x): (Nc²-1, n_links)
  real(dp), allocatable :: plaquette_values(:)        ! For measuring field strength
  
  ! Wilson loops for observables
  complex(dp), allocatable :: wilson_loops(:)
  
contains

  !=====================================================================
  ! Subroutine: initialize_gauge_field
  ! Purpose: Initialize gauge field configuration
  !=====================================================================
  subroutine initialize_gauge_field(init_type)
    implicit none
    character(len=*), intent(in), optional :: init_type
    character(len=20) :: init_method
    integer :: ilink, a
    integer :: n_gen
    real(dp) :: random_val
    
    n_gen = n_colors**2 - 1
    
    if (present(init_type)) then
      init_method = init_type
    else
      init_method = 'cold'
    end if
    
    write(*,'(A)') ' Initializing gauge field...'
    
    ! Allocate arrays
    if (allocated(link_matrices)) deallocate(link_matrices)
    if (allocated(electric_field)) deallocate(electric_field)
    
    allocate(link_matrices(n_colors, n_colors, n_links))
    allocate(electric_field(n_gen, n_links))
    
    ! Initialize link matrices
    select case (trim(init_method))
    case ('cold')
      ! Cold start: all links to identity
      link_matrices = (0.0_dp, 0.0_dp)
      do ilink = 1, n_links
        do a = 1, n_colors
          link_matrices(a, a, ilink) = (1.0_dp, 0.0_dp)
        end do
      end do
      write(*,'(A)') '   Cold start: links initialized to identity'
      
    case ('hot')
      ! Hot start: random SU(Nc) matrices
      do ilink = 1, n_links
        call random_su_matrix(link_matrices(:,:,ilink))
      end do
      write(*,'(A)') '   Hot start: links initialized randomly'
      
    case default
      write(*,*) 'Warning: Unknown initialization type. Using cold start.'
      link_matrices = (0.0_dp, 0.0_dp)
      do ilink = 1, n_links
        do a = 1, n_colors
          link_matrices(a, a, ilink) = (1.0_dp, 0.0_dp)
        end do
      end do
    end select
    
    ! Initialize electric field (canonical momentum)
    electric_field = 0.0_dp
    
    write(*,'(A,I8,A)') '   Initialized ', n_links, ' link variables'
    
  end subroutine initialize_gauge_field
  
  !=====================================================================
  ! Subroutine: random_su_matrix
  ! Purpose: Generate random SU(Nc) matrix
  !=====================================================================
  subroutine random_su_matrix(U)
    implicit none
    complex(dp), intent(out) :: U(n_colors, n_colors)
    complex(dp) :: H(n_colors, n_colors)
    complex(dp) :: V(n_colors, n_colors)
    real(dp) :: random_vals(2)
    integer :: i, j, n_gen, a
    real(dp) :: alpha
    
    n_gen = n_colors**2 - 1
    
    ! Generate random Hermitian matrix from generators
    H = (0.0_dp, 0.0_dp)
    do a = 1, n_gen
      call random_number(alpha)
      alpha = (alpha - 0.5_dp) * 2.0_dp * pi
      H = H + alpha * generators(:,:,a)
    end do
    
    ! Exponentiate: U = exp(i H)
    U = matrix_exp(zi * H)
    
    ! Project to SU(Nc) by normalizing determinant
    call project_to_sun(U)
    
  end subroutine random_su_matrix
  
  !=====================================================================
  ! Function: matrix_exp
  ! Purpose: Matrix exponential via Taylor series
  !=====================================================================
  function matrix_exp(M) result(expM)
    implicit none
    complex(dp), intent(in) :: M(:,:)
    complex(dp) :: expM(size(M,1), size(M,2))
    complex(dp) :: term(size(M,1), size(M,2))
    complex(dp) :: M_power(size(M,1), size(M,2))
    integer :: n, i
    real(dp) :: fact
    
    ! Initialize identity
    expM = (0.0_dp, 0.0_dp)
    do i = 1, size(M,1)
      expM(i,i) = (1.0_dp, 0.0_dp)
    end do
    
    M_power = expM
    fact = 1.0_dp
    
    do n = 1, 20
      M_power = matmul(M_power, M)
      fact = fact * real(n, dp)
      term = M_power / fact
      expM = expM + term
      
      if (maxval(abs(term)) < 1.0e-12_dp) exit
    end do
    
  end function matrix_exp
  
  !=====================================================================
  ! Subroutine: project_to_sun
  ! Purpose: Project matrix to SU(N) by Gram-Schmidt
  !=====================================================================
  subroutine project_to_sun(U)
    implicit none
    complex(dp), intent(inout) :: U(:,:)
    complex(dp) :: Q(size(U,1), size(U,2))
    complex(dp) :: det
    integer :: i, j, k
    real(dp) :: norm
    
    ! Gram-Schmidt orthogonalization
    Q = U
    do j = 1, n_colors
      ! Subtract projections on previous columns
      do k = 1, j-1
        Q(:,j) = Q(:,j) - dot_product(Q(:,k), U(:,j)) * Q(:,k)
      end do
      ! Normalize
      norm = sqrt(real(dot_product(Q(:,j), Q(:,j)), dp))
      if (norm > 1.0e-14_dp) then
        Q(:,j) = Q(:,j) / norm
      end if
    end do
    
    ! Ensure det = 1
    det = matrix_det(Q)
    Q(:,n_colors) = Q(:,n_colors) / det**(1.0_dp/n_colors)
    
    U = Q
    
  end subroutine project_to_sun
  
  !=====================================================================
  ! Function: matrix_det
  ! Purpose: Calculate determinant (simple implementation)
  !=====================================================================
  function matrix_det(M) result(det)
    implicit none
    complex(dp), intent(in) :: M(:,:)
    complex(dp) :: det
    complex(dp) :: L(size(M,1), size(M,2)), U_mat(size(M,1), size(M,2))
    integer :: i, j, k, n
    
    n = size(M, 1)
    
    if (n == 1) then
      det = M(1,1)
      return
    else if (n == 2) then
      det = M(1,1)*M(2,2) - M(1,2)*M(2,1)
      return
    else if (n == 3) then
      det = M(1,1)*(M(2,2)*M(3,3) - M(2,3)*M(3,2)) - &
            M(1,2)*(M(2,1)*M(3,3) - M(2,3)*M(3,1)) + &
            M(1,3)*(M(2,1)*M(3,2) - M(2,2)*M(3,1))
      return
    end if
    
    ! For larger matrices, use product of diagonal (approximate)
    det = (1.0_dp, 0.0_dp)
    do i = 1, n
      det = det * M(i,i)
    end do
    
  end function matrix_det
  
  !=====================================================================
  ! Function: calculate_plaquette
  ! Purpose: Calculate plaquette (elementary magnetic field)
  !=====================================================================
  function calculate_plaquette(isite, mu, nu) result(plaq)
    implicit none
    integer, intent(in) :: isite, mu, nu
    complex(dp) :: plaq(n_colors, n_colors)
    integer :: plaq_sites(4)
    integer :: link1, link2, link3, link4
    
    ! Get the four corners of the plaquette
    call get_plaquette_sites(isite, mu, nu, plaq_sites)
    
    ! Find corresponding links
    link1 = (plaq_sites(1) - 1) * n_directions + mu
    link2 = (plaq_sites(2) - 1) * n_directions + nu
    link3 = (plaq_sites(4) - 1) * n_directions + mu
    link4 = (plaq_sites(1) - 1) * n_directions + nu
    
    ! Plaquette = U_μ(x) U_ν(x+μ) U†_μ(x+ν) U†_ν(x)
    plaq = matmul(link_matrices(:,:,link1), link_matrices(:,:,link2))
    plaq = matmul(plaq, su_adjoint(link_matrices(:,:,link3)))
    plaq = matmul(plaq, su_adjoint(link_matrices(:,:,link4)))
    
  end function calculate_plaquette
  
  !=====================================================================
  ! Function: plaquette_average
  ! Purpose: Calculate average plaquette value
  !=====================================================================
  function plaquette_average() result(avg_plaq)
    implicit none
    real(dp) :: avg_plaq
    complex(dp) :: plaq(n_colors, n_colors)
    integer :: isite, mu, nu
    real(dp) :: sum_plaq
    integer :: n_plaq
    
    sum_plaq = 0.0_dp
    n_plaq = 0
    
    do isite = 1, total_sites
      do mu = 1, n_directions
        do nu = mu+1, n_directions
          plaq = calculate_plaquette(isite, mu, nu)
          sum_plaq = sum_plaq + real(su_trace(plaq), dp) / n_colors
          n_plaq = n_plaq + 1
        end do
      end do
    end do
    
    avg_plaq = sum_plaq / n_plaq
    
  end function plaquette_average
  
  !=====================================================================
  ! Function: field_strength_tensor
  ! Purpose: Calculate field strength F_μν at a site
  !=====================================================================
  function field_strength_tensor(isite, mu, nu) result(F)
    implicit none
    integer, intent(in) :: isite, mu, nu
    complex(dp) :: F(n_colors, n_colors)
    complex(dp) :: plaq(n_colors, n_colors)
    integer :: a
    
    ! F_μν = (1/ig) [U_P - U_P†] where U_P is the plaquette
    plaq = calculate_plaquette(isite, mu, nu)
    
    ! Extract anti-Hermitian traceless part
    F = (plaq - su_adjoint(plaq)) / (2.0_dp * zi)
    
    ! Remove trace
    F = F - su_trace(F) / n_colors
    
  end function field_strength_tensor
  
  !=====================================================================
  ! Subroutine: cleanup_gauge_field
  ! Purpose: Deallocate gauge field arrays
  !=====================================================================
  subroutine cleanup_gauge_field()
    implicit none
    if (allocated(link_matrices)) deallocate(link_matrices)
    if (allocated(electric_field)) deallocate(electric_field)
    if (allocated(plaquette_values)) deallocate(plaquette_values)
    if (allocated(wilson_loops)) deallocate(wilson_loops)
  end subroutine cleanup_gauge_field

end module mod_gauge_field
