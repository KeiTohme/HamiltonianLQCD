!***********************************************************************
! Module: mod_hamiltonian
! Purpose: Implement Hamiltonian for Lattice QCD
!          H = H_gauge + H_fermion
!          Gauge part: Electric field (E²) and Magnetic field (B²) terms
!***********************************************************************
module mod_hamiltonian
  use mod_parameters
  use mod_lattice
  use mod_gauge_field
  use mod_su_algebra
  use mod_wilson_fermion
  use mod_staggered_fermion
  implicit none
  
  ! Hamiltonian components
  real(dp) :: H_electric, H_magnetic, H_fermion, H_total
  
contains

  !=====================================================================
  ! Function: calculate_hamiltonian
  ! Purpose: Calculate total Hamiltonian H = H_E + H_B + H_F
  !=====================================================================
  function calculate_hamiltonian() result(H)
    implicit none
    real(dp) :: H
    
    ! Electric field energy (canonical momentum)
    H_electric = electric_field_energy()
    
    ! Magnetic field energy (from plaquettes)
    H_magnetic = magnetic_field_energy()
    
    ! Fermion contribution
    if (trim(fermion_type) == 'wilson') then
      H_fermion = wilson_fermion_action()
    else if (trim(fermion_type) == 'staggered') then
      H_fermion = staggered_fermion_action()
    else
      H_fermion = 0.0_dp
    end if
    
    H_total = H_electric + H_magnetic + H_fermion
    H = H_total
    
  end function calculate_hamiltonian
  
  !=====================================================================
  ! Function: electric_field_energy
  ! Purpose: Calculate electric field energy H_E = (g²/2) Σ_x,μ,a [E^a_μ(x)]²
  !=====================================================================
  function electric_field_energy() result(H_E)
    implicit none
    real(dp) :: H_E
    integer :: ilink, a
    integer :: n_gen
    
    n_gen = n_colors**2 - 1
    H_E = 0.0_dp
    
    ! Sum over all links and color components
    do ilink = 1, n_links
      do a = 1, n_gen
        H_E = H_E + abs(electric_field(a, ilink))**2
      end do
    end do
    
    ! Multiply by coupling factor
    H_E = 0.5_dp * gauge_coupling**2 * H_E
    
  end function electric_field_energy
  
  !=====================================================================
  ! Function: magnetic_field_energy
  ! Purpose: Calculate magnetic field energy H_B = (1/g²) Σ_x,μ<ν Tr[1-Re(U_P)]
  !          where U_P is the plaquette
  !=====================================================================
  function magnetic_field_energy() result(H_B)
    implicit none
    real(dp) :: H_B
    complex(dp) :: plaq(n_colors, n_colors)
    integer :: isite, mu, nu
    real(dp) :: plaq_sum
    
    H_B = 0.0_dp
    
    ! Sum over all plaquettes
    do isite = 1, total_sites
      do mu = 1, min(n_directions, 4)
        do nu = mu+1, min(n_directions, 4)
          plaq = calculate_plaquette(isite, mu, nu)
          
          ! Contribution: 1 - Re(Tr[U_P])/Nc
          plaq_sum = 1.0_dp - real(su_trace(plaq), dp) / n_colors
          H_B = H_B + plaq_sum
        end do
      end do
    end do
    
    ! Multiply by coupling factor
    H_B = (2.0_dp * n_colors / gauge_coupling**2) * H_B
    
  end function magnetic_field_energy
  
  !=====================================================================
  ! Subroutine: calculate_equations_of_motion
  ! Purpose: Calculate dU/dt and dE/dt from Hamilton's equations
  !          dU/dt = δH/δE,  dE/dt = -δH/δU
  !=====================================================================
  subroutine calculate_equations_of_motion(dU_dt, dE_dt)
    implicit none
    complex(dp), intent(out) :: dU_dt(n_colors, n_colors, n_links)
    complex(dp), intent(out) :: dE_dt(n_colors**2-1, n_links)
    integer :: ilink, a
    integer :: n_gen
    
    n_gen = n_colors**2 - 1
    
    ! dU/dt = i g² [E^a T^a, U]
    do ilink = 1, n_links
      dU_dt(:,:,ilink) = (0.0_dp, 0.0_dp)
      do a = 1, n_gen
        dU_dt(:,:,ilink) = dU_dt(:,:,ilink) + &
          zi * gauge_coupling**2 * electric_field(a, ilink) * &
          (matmul(generators(:,:,a), link_matrices(:,:,ilink)) - &
           matmul(link_matrices(:,:,ilink), generators(:,:,a)))
      end do
    end do
    
    ! dE/dt = -δH_magnetic/δU (force from magnetic field)
    call calculate_gauge_force(dE_dt)
    
  end subroutine calculate_equations_of_motion
  
  !=====================================================================
  ! Subroutine: calculate_gauge_force
  ! Purpose: Calculate force on electric field from magnetic energy
  !          F^a_μ(x) = -δH_B/δU^a_μ(x)
  !=====================================================================
  subroutine calculate_gauge_force(force)
    implicit none
    complex(dp), intent(out) :: force(n_colors**2-1, n_links)
    integer :: isite, mu, nu, ilink, a
    integer :: neighbor
    integer :: n_gen
    complex(dp) :: staple(n_colors, n_colors)
    complex(dp) :: U(n_colors, n_colors)
    complex(dp) :: force_matrix(n_colors, n_colors)
    
    n_gen = n_colors**2 - 1
    force = (0.0_dp, 0.0_dp)
    
    do isite = 1, total_sites
      do mu = 1, min(n_directions, 4)
        ilink = (isite - 1) * n_directions + mu
        U = link_matrices(:,:,ilink)
        
        ! Calculate staples in all perpendicular directions
        staple = (0.0_dp, 0.0_dp)
        do nu = 1, min(n_directions, 4)
          if (nu /= mu) then
            staple = staple + calculate_staple(isite, mu, nu)
          end if
        end do
        
        ! Force matrix: dS/dU = -β/Nc [staple U†]
        force_matrix = -beta / n_colors * &
                      matmul(staple, su_adjoint(U))
        
        ! Project onto generators
        do a = 1, n_gen
          force(a, ilink) = 2.0_dp * real( &
            su_trace(matmul(generators(:,:,a), force_matrix)), dp)
        end do
      end do
    end do
    
  end subroutine calculate_gauge_force
  
  !=====================================================================
  ! Function: calculate_staple
  ! Purpose: Calculate staple contributions for force calculation
  !          Staple = U_ν(x+μ)U†_μ(x+ν)U†_ν(x) + U†_ν(x+μ-ν)U†_μ(x-ν)U_ν(x-ν)
  !=====================================================================
  function calculate_staple(isite, mu, nu) result(staple)
    implicit none
    integer, intent(in) :: isite, mu, nu
    complex(dp) :: staple(n_colors, n_colors)
    integer :: site_pmu, site_pnu, site_mnu
    integer :: link1, link2, link3, link4, link5, link6
    complex(dp) :: term1(n_colors, n_colors)
    complex(dp) :: term2(n_colors, n_colors)
    
    ! Forward staple: U_ν(x+μ) U†_μ(x+ν) U†_ν(x)
    site_pmu = lattice_sites(isite)%neighbors(mu, 1)
    site_pnu = lattice_sites(isite)%neighbors(nu, 1)
    
    link1 = (site_pmu - 1) * n_directions + nu
    link2 = (site_pnu - 1) * n_directions + mu
    link3 = (isite - 1) * n_directions + nu
    
    term1 = matmul(link_matrices(:,:,link1), su_adjoint(link_matrices(:,:,link2)))
    term1 = matmul(term1, su_adjoint(link_matrices(:,:,link3)))
    
    ! Backward staple: U†_ν(x+μ-ν) U†_μ(x-ν) U_ν(x-ν)
    site_mnu = lattice_sites(isite)%neighbors(nu, 2)
    site_pmu = lattice_sites(site_mnu)%neighbors(mu, 1)
    
    link4 = (site_pmu - 1) * n_directions + nu
    link5 = (site_mnu - 1) * n_directions + mu
    link6 = (site_mnu - 1) * n_directions + nu
    
    term2 = matmul(su_adjoint(link_matrices(:,:,link4)), &
                   su_adjoint(link_matrices(:,:,link5)))
    term2 = matmul(term2, link_matrices(:,:,link6))
    
    staple = term1 + term2
    
  end function calculate_staple
  
  !=====================================================================
  ! Function: calculate_energy_density
  ! Purpose: Calculate energy density at each site
  !=====================================================================
  function calculate_energy_density() result(energy_density)
    implicit none
    real(dp) :: energy_density(total_sites)
    integer :: isite, mu, nu, ilink, a
    integer :: n_gen
    real(dp) :: E_local, B_local
    complex(dp) :: plaq(n_colors, n_colors)
    
    n_gen = n_colors**2 - 1
    energy_density = 0.0_dp
    
    do isite = 1, total_sites
      E_local = 0.0_dp
      B_local = 0.0_dp
      
      ! Electric field contribution
      do mu = 1, min(n_directions, 4)
        ilink = (isite - 1) * n_directions + mu
        do a = 1, n_gen
          E_local = E_local + abs(electric_field(a, ilink))**2
        end do
      end do
      E_local = 0.5_dp * gauge_coupling**2 * E_local
      
      ! Magnetic field contribution (plaquettes starting from this site)
      do mu = 1, min(n_directions, 4)
        do nu = mu+1, min(n_directions, 4)
          plaq = calculate_plaquette(isite, mu, nu)
          B_local = B_local + (1.0_dp - real(su_trace(plaq), dp) / n_colors)
        end do
      end do
      B_local = (2.0_dp * n_colors / gauge_coupling**2) * B_local
      
      energy_density(isite) = E_local + B_local
    end do
    
  end function calculate_energy_density
  
  !=====================================================================
  ! Function: calculate_action
  ! Purpose: Calculate total action S = S_gauge + S_fermion
  !=====================================================================
  function calculate_action() result(action)
    implicit none
    real(dp) :: action
    real(dp) :: S_gauge, S_fermion
    
    ! Wilson gauge action
    S_gauge = -beta * plaquette_average() * total_sites * &
              n_directions * (n_directions - 1) / 2
    
    ! Fermion action
    if (trim(fermion_type) == 'wilson') then
      S_fermion = wilson_fermion_action()
    else if (trim(fermion_type) == 'staggered') then
      S_fermion = staggered_fermion_action()
    else
      S_fermion = 0.0_dp
    end if
    
    action = S_gauge + S_fermion
    
  end function calculate_action

end module mod_hamiltonian
