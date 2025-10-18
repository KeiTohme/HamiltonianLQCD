!***********************************************************************
! Module: mod_wilson_fermion
! Purpose: Implement Wilson fermion formulation on the lattice
!          Including Dirac operator and fermion action
!***********************************************************************
module mod_wilson_fermion
  use mod_parameters
  use mod_lattice
  use mod_gauge_field
  use mod_su_algebra
  implicit none
  
  ! Wilson fermion field: ψ(x, spin, color)
  complex(dp), allocatable :: fermion_field(:,:,:)    ! (total_sites, 4, Nc)
  complex(dp), allocatable :: fermion_mom(:,:,:)      ! Conjugate momentum
  
  ! Gamma matrices (Dirac matrices)
  complex(dp), allocatable :: gamma(:,:,:)            ! (4, 4, 0:5)
  ! γ⁰, γ¹, γ², γ³, γ⁵
  
  ! Wilson parameter
  real(dp) :: wilson_r = 1.0_dp
  
contains

  !=====================================================================
  ! Subroutine: initialize_wilson_fermion
  ! Purpose: Initialize Wilson fermion fields and gamma matrices
  !=====================================================================
  subroutine initialize_wilson_fermion()
    implicit none
    integer :: isite, ispin, icolor
    real(dp) :: rand_re, rand_im
    
    write(*,'(A)') ' Initializing Wilson fermion...'
    
    ! Allocate fermion fields
    if (allocated(fermion_field)) deallocate(fermion_field)
    if (allocated(fermion_mom)) deallocate(fermion_mom)
    
    allocate(fermion_field(total_sites, 4, n_colors))
    allocate(fermion_mom(total_sites, 4, n_colors))
    
    ! Initialize with random values (or zero)
    do isite = 1, total_sites
      do ispin = 1, 4
        do icolor = 1, n_colors
          call random_number(rand_re)
          call random_number(rand_im)
          fermion_field(isite, ispin, icolor) = &
            cmplx(rand_re - 0.5_dp, rand_im - 0.5_dp, dp)
        end do
      end do
    end do
    
    fermion_mom = (0.0_dp, 0.0_dp)
    
    ! Initialize gamma matrices
    call setup_gamma_matrices()
    
    write(*,'(A,I8,A)') '   Initialized Wilson fermion field with ', &
                        total_sites * 4 * n_colors, ' components'
    
  end subroutine initialize_wilson_fermion
  
  !=====================================================================
  ! Subroutine: setup_gamma_matrices
  ! Purpose: Define Dirac gamma matrices (Dirac representation)
  !=====================================================================
  subroutine setup_gamma_matrices()
    implicit none
    
    if (allocated(gamma)) deallocate(gamma)
    allocate(gamma(4, 4, 0:5))
    
    gamma = (0.0_dp, 0.0_dp)
    
    ! γ⁰ (Dirac representation)
    gamma(1,1,0) = (1.0_dp, 0.0_dp)
    gamma(2,2,0) = (1.0_dp, 0.0_dp)
    gamma(3,3,0) = (-1.0_dp, 0.0_dp)
    gamma(4,4,0) = (-1.0_dp, 0.0_dp)
    
    ! γ¹
    gamma(1,4,1) = zi
    gamma(2,3,1) = zi
    gamma(3,2,1) = -zi
    gamma(4,1,1) = -zi
    
    ! γ²
    gamma(1,4,2) = (1.0_dp, 0.0_dp)
    gamma(2,3,2) = (-1.0_dp, 0.0_dp)
    gamma(3,2,2) = (-1.0_dp, 0.0_dp)
    gamma(4,1,2) = (1.0_dp, 0.0_dp)
    
    ! γ³
    gamma(1,3,3) = zi
    gamma(2,4,3) = -zi
    gamma(3,1,3) = -zi
    gamma(4,2,3) = zi
    
    ! γ⁵ = γ⁰γ¹γ²γ³
    gamma(1,1,5) = (1.0_dp, 0.0_dp)
    gamma(2,2,5) = (1.0_dp, 0.0_dp)
    gamma(3,3,5) = (-1.0_dp, 0.0_dp)
    gamma(4,4,5) = (-1.0_dp, 0.0_dp)
    
  end subroutine setup_gamma_matrices
  
  !=====================================================================
  ! Function: wilson_dirac_operator
  ! Purpose: Apply Wilson-Dirac operator to fermion field
  !          D_W = m + (4 + m) - 1/2 Σ_μ [(1-γ_μ)U_μ(x)δ_{x+μ,y} + 
  !                                         (1+γ_μ)U†_μ(x-μ)δ_{x-μ,y}]
  !=====================================================================
  function wilson_dirac_operator(psi) result(D_psi)
    implicit none
    complex(dp), intent(in) :: psi(:,:,:)  ! (sites, spin, color)
    complex(dp) :: D_psi(size(psi,1), size(psi,2), size(psi,3))
    integer :: isite, mu, neighbor_fwd, neighbor_bwd
    integer :: link_fwd, link_bwd
    integer :: ispin, jspin, icolor, jcolor
    complex(dp) :: U_fwd(n_colors, n_colors)
    complex(dp) :: U_bwd(n_colors, n_colors)
    complex(dp) :: psi_fwd(4, n_colors), psi_bwd(4, n_colors)
    complex(dp) :: temp_spin(4), temp_color(n_colors)
    real(dp) :: kappa
    
    ! Hopping parameter κ = 1/(2m + 8)
    kappa = 1.0_dp / (2.0_dp * fermion_mass + 8.0_dp)
    
    ! Mass term and diagonal
    D_psi = fermion_mass * psi
    
    ! Hopping terms
    do isite = 1, total_sites
      do mu = 1, min(n_directions, 4)  ! Only spacetime directions
        ! Forward hopping
        neighbor_fwd = lattice_sites(isite)%neighbors(mu, 1)
        link_fwd = (isite - 1) * n_directions + mu
        U_fwd = link_matrices(:,:,link_fwd)
        psi_fwd = psi(neighbor_fwd, :, :)
        
        ! Backward hopping
        neighbor_bwd = lattice_sites(isite)%neighbors(mu, 2)
        link_bwd = (neighbor_bwd - 1) * n_directions + mu
        U_bwd = su_adjoint(link_matrices(:,:,link_bwd))
        psi_bwd = psi(neighbor_bwd, :, :)
        
        ! Apply (1 - γ_μ) and (1 + γ_μ) projectors
        do icolor = 1, n_colors
          ! Forward: (1 - γ_μ) U_μ ψ(x+μ)
          temp_spin = (0.0_dp, 0.0_dp)
          do jspin = 1, 4
            do jcolor = 1, n_colors
              temp_spin(jspin) = temp_spin(jspin) + &
                U_fwd(icolor, jcolor) * psi_fwd(jspin, jcolor)
            end do
          end do
          
          do ispin = 1, 4
            D_psi(isite, ispin, icolor) = D_psi(isite, ispin, icolor) + &
              kappa * (temp_spin(ispin) - &
                      sum(gamma(ispin, :, mu) * temp_spin(:)))
          end do
          
          ! Backward: (1 + γ_μ) U†_μ ψ(x-μ)
          temp_spin = (0.0_dp, 0.0_dp)
          do jspin = 1, 4
            do jcolor = 1, n_colors
              temp_spin(jspin) = temp_spin(jspin) + &
                U_bwd(icolor, jcolor) * psi_bwd(jspin, jcolor)
            end do
          end do
          
          do ispin = 1, 4
            D_psi(isite, ispin, icolor) = D_psi(isite, ispin, icolor) + &
              kappa * (temp_spin(ispin) + &
                      sum(gamma(ispin, :, mu) * temp_spin(:)))
          end do
        end do
      end do
    end do
    
  end function wilson_dirac_operator
  
  !=====================================================================
  ! Function: wilson_fermion_action
  ! Purpose: Calculate Wilson fermion action S_F = ψ̄ D_W ψ
  !=====================================================================
  function wilson_fermion_action() result(action)
    implicit none
    real(dp) :: action
    complex(dp) :: D_psi(total_sites, 4, n_colors)
    complex(dp) :: sum_val
    integer :: isite, ispin, icolor
    
    ! Apply Dirac operator
    D_psi = wilson_dirac_operator(fermion_field)
    
    ! Calculate ψ̄ D_W ψ
    sum_val = (0.0_dp, 0.0_dp)
    do isite = 1, total_sites
      do ispin = 1, 4
        do icolor = 1, n_colors
          sum_val = sum_val + &
            conjg(fermion_field(isite, ispin, icolor)) * &
            D_psi(isite, ispin, icolor)
        end do
      end do
    end do
    
    action = real(sum_val, dp)
    
  end function wilson_fermion_action
  
  !=====================================================================
  ! Function: fermion_bilinear
  ! Purpose: Calculate fermion bilinear ψ̄ Γ ψ for operator Γ
  !=====================================================================
  function fermion_bilinear(gamma_idx) result(bilinear)
    implicit none
    integer, intent(in) :: gamma_idx  ! 0-5 for γ⁰-γ³, γ⁵
    complex(dp) :: bilinear
    integer :: isite, ispin, jspin, icolor
    complex(dp) :: temp
    
    bilinear = (0.0_dp, 0.0_dp)
    
    do isite = 1, total_sites
      do icolor = 1, n_colors
        temp = (0.0_dp, 0.0_dp)
        do ispin = 1, 4
          do jspin = 1, 4
            temp = temp + &
              conjg(fermion_field(isite, ispin, icolor)) * &
              gamma(ispin, jspin, gamma_idx) * &
              fermion_field(isite, jspin, icolor)
          end do
        end do
        bilinear = bilinear + temp
      end do
    end do
    
    bilinear = bilinear / (total_sites * n_colors)
    
  end function fermion_bilinear
  
  !=====================================================================
  ! Subroutine: cleanup_wilson_fermion
  ! Purpose: Deallocate Wilson fermion arrays
  !=====================================================================
  subroutine cleanup_wilson_fermion()
    implicit none
    if (allocated(fermion_field)) deallocate(fermion_field)
    if (allocated(fermion_mom)) deallocate(fermion_mom)
    if (allocated(gamma)) deallocate(gamma)
  end subroutine cleanup_wilson_fermion

end module mod_wilson_fermion
