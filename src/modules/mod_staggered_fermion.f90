!***********************************************************************
! Module: mod_staggered_fermion
! Purpose: Implement Staggered (Kogut-Susskind) fermion formulation
!          Reduced degrees of freedom: one component per site
!***********************************************************************
module mod_staggered_fermion
  use mod_parameters
  use mod_lattice
  use mod_gauge_field
  use mod_su_algebra
  implicit none
  
  ! Staggered fermion field: χ(x, color) - one spin component per site
  complex(dp), allocatable :: staggered_field(:,:)      ! (total_sites, Nc)
  complex(dp), allocatable :: staggered_mom(:,:)        ! Conjugate momentum
  
  ! Staggered phase factors η_μ(x)
  real(dp), allocatable :: eta_phase(:,:)               ! (total_sites, directions)
  
contains

  !=====================================================================
  ! Subroutine: initialize_staggered_fermion
  ! Purpose: Initialize staggered fermion fields and phase factors
  !=====================================================================
  subroutine initialize_staggered_fermion()
    implicit none
    integer :: isite, icolor, mu
    real(dp) :: rand_re, rand_im
    
    write(*,'(A)') ' Initializing Staggered fermion...'
    
    ! Allocate fermion fields
    if (allocated(staggered_field)) deallocate(staggered_field)
    if (allocated(staggered_mom)) deallocate(staggered_mom)
    if (allocated(eta_phase)) deallocate(eta_phase)
    
    allocate(staggered_field(total_sites, n_colors))
    allocate(staggered_mom(total_sites, n_colors))
    allocate(eta_phase(total_sites, n_directions))
    
    ! Initialize with random values
    do isite = 1, total_sites
      do icolor = 1, n_colors
        call random_number(rand_re)
        call random_number(rand_im)
        staggered_field(isite, icolor) = &
          cmplx(rand_re - 0.5_dp, rand_im - 0.5_dp, dp)
      end do
    end do
    
    staggered_mom = (0.0_dp, 0.0_dp)
    
    ! Calculate staggered phase factors
    call setup_staggered_phases()
    
    write(*,'(A,I8,A)') '   Initialized Staggered fermion field with ', &
                        total_sites * n_colors, ' components'
    
  end subroutine initialize_staggered_fermion
  
  !=====================================================================
  ! Subroutine: setup_staggered_phases
  ! Purpose: Calculate staggered phase factors η_μ(x) = (-1)^(x₁+...+x_{μ-1})
  !=====================================================================
  subroutine setup_staggered_phases()
    implicit none
    integer :: isite, mu, nu
    integer :: phase_sum
    integer, allocatable :: coords(:)
    
    allocate(coords(d_euclid))
    
    do isite = 1, total_sites
      coords = lattice_sites(isite)%coords
      
      do mu = 1, n_directions
        phase_sum = 0
        do nu = 1, mu - 1
          if (nu <= d_euclid) then
            phase_sum = phase_sum + coords(nu)
          end if
        end do
        
        ! η_μ(x) = (-1)^(sum)
        if (mod(phase_sum, 2) == 0) then
          eta_phase(isite, mu) = 1.0_dp
        else
          eta_phase(isite, mu) = -1.0_dp
        end if
      end do
    end do
    
    deallocate(coords)
    
  end subroutine setup_staggered_phases
  
  !=====================================================================
  ! Function: staggered_dirac_operator
  ! Purpose: Apply staggered Dirac operator
  !          D_stag = m χ(x) + 1/2 Σ_μ η_μ(x) [U_μ(x)χ(x+μ) - U†_μ(x-μ)χ(x-μ)]
  !=====================================================================
  function staggered_dirac_operator(chi) result(D_chi)
    implicit none
    complex(dp), intent(in) :: chi(:,:)  ! (sites, color)
    complex(dp) :: D_chi(size(chi,1), size(chi,2))
    integer :: isite, mu, neighbor_fwd, neighbor_bwd
    integer :: link_fwd, link_bwd
    integer :: icolor, jcolor
    complex(dp) :: U_fwd(n_colors, n_colors)
    complex(dp) :: U_bwd(n_colors, n_colors)
    complex(dp) :: chi_fwd(n_colors), chi_bwd(n_colors)
    real(dp) :: eta
    
    ! Mass term
    D_chi = fermion_mass * chi
    
    ! Hopping terms with staggered phases
    do isite = 1, total_sites
      do mu = 1, min(n_directions, 4)  ! Spacetime directions only
        eta = eta_phase(isite, mu)
        
        ! Forward neighbor
        neighbor_fwd = lattice_sites(isite)%neighbors(mu, 1)
        link_fwd = (isite - 1) * n_directions + mu
        U_fwd = link_matrices(:,:,link_fwd)
        chi_fwd = chi(neighbor_fwd, :)
        
        ! Backward neighbor
        neighbor_bwd = lattice_sites(isite)%neighbors(mu, 2)
        link_bwd = (neighbor_bwd - 1) * n_directions + mu
        U_bwd = su_adjoint(link_matrices(:,:,link_bwd))
        chi_bwd = chi(neighbor_bwd, :)
        
        ! Add hopping contribution
        do icolor = 1, n_colors
          ! Forward hopping: η_μ(x) U_μ(x) χ(x+μ)
          do jcolor = 1, n_colors
            D_chi(isite, icolor) = D_chi(isite, icolor) + &
              0.5_dp * eta * U_fwd(icolor, jcolor) * chi_fwd(jcolor)
          end do
          
          ! Backward hopping: -η_μ(x) U†_μ(x-μ) χ(x-μ)
          do jcolor = 1, n_colors
            D_chi(isite, icolor) = D_chi(isite, icolor) - &
              0.5_dp * eta * U_bwd(icolor, jcolor) * chi_bwd(jcolor)
          end do
        end do
      end do
    end do
    
  end function staggered_dirac_operator
  
  !=====================================================================
  ! Function: staggered_fermion_action
  ! Purpose: Calculate staggered fermion action S_F = χ̄ D_stag χ
  !=====================================================================
  function staggered_fermion_action() result(action)
    implicit none
    real(dp) :: action
    complex(dp) :: D_chi(total_sites, n_colors)
    complex(dp) :: sum_val
    integer :: isite, icolor
    
    ! Apply staggered Dirac operator
    D_chi = staggered_dirac_operator(staggered_field)
    
    ! Calculate χ̄ D_stag χ
    sum_val = (0.0_dp, 0.0_dp)
    do isite = 1, total_sites
      do icolor = 1, n_colors
        sum_val = sum_val + &
          conjg(staggered_field(isite, icolor)) * D_chi(isite, icolor)
      end do
    end do
    
    action = real(sum_val, dp)
    
  end function staggered_fermion_action
  
  !=====================================================================
  ! Function: staggered_condensate
  ! Purpose: Calculate chiral condensate ⟨χ̄χ⟩
  !=====================================================================
  function staggered_condensate() result(condensate)
    implicit none
    complex(dp) :: condensate
    integer :: isite, icolor
    
    condensate = (0.0_dp, 0.0_dp)
    
    do isite = 1, total_sites
      do icolor = 1, n_colors
        condensate = condensate + &
          conjg(staggered_field(isite, icolor)) * &
          staggered_field(isite, icolor)
      end do
    end do
    
    condensate = condensate / (total_sites * n_colors)
    
  end function staggered_condensate
  
  !=====================================================================
  ! Function: staggered_pion_correlator
  ! Purpose: Calculate pion correlator for staggered fermions
  !=====================================================================
  function staggered_pion_correlator(separation) result(correlator)
    implicit none
    integer, intent(in) :: separation
    complex(dp) :: correlator
    integer :: isite, jsite, icolor
    real(dp) :: distance
    integer :: count
    
    correlator = (0.0_dp, 0.0_dp)
    count = 0
    
    ! Sum over all site pairs with given separation
    do isite = 1, total_sites
      do jsite = 1, total_sites
        ! Check if separation matches (simplified - just use one direction)
        if (abs(lattice_sites(isite)%coords(1) - &
                lattice_sites(jsite)%coords(1)) == separation) then
          do icolor = 1, n_colors
            correlator = correlator + &
              conjg(staggered_field(isite, icolor)) * &
              staggered_field(jsite, icolor)
          end do
          count = count + 1
        end if
      end do
    end do
    
    if (count > 0) then
      correlator = correlator / count
    end if
    
  end function staggered_pion_correlator
  
  !=====================================================================
  ! Subroutine: cleanup_staggered_fermion
  ! Purpose: Deallocate staggered fermion arrays
  !=====================================================================
  subroutine cleanup_staggered_fermion()
    implicit none
    if (allocated(staggered_field)) deallocate(staggered_field)
    if (allocated(staggered_mom)) deallocate(staggered_mom)
    if (allocated(eta_phase)) deallocate(eta_phase)
  end subroutine cleanup_staggered_fermion

end module mod_staggered_fermion
