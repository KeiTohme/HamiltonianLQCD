module gauge_field
  use kinds
  use lattice_geometry
  use sun_group
  implicit none
  private

  type :: gauge_t
     integer :: Nc
     type(lattice_t) :: lat
     complex(dp), allocatable :: U(:,:,:,:) ! (Nc,Nc,ndir,volume)
     complex(dp), allocatable :: E(:,:,:,:) ! Hermitian traceless algebra element (Nc,Nc,ndir,volume)
  end type gauge_t

  public :: gauge_t, init_gauge, cold_start, hot_start, step_leapfrog, measure_energies

contains

  subroutine init_gauge(G, lat, Nc)
    type(gauge_t), intent(out) :: G
    type(lattice_t), intent(in) :: lat
    integer, intent(in) :: Nc
    integer :: v

    G%Nc = Nc
    G%lat = lat
    allocate(G%U(Nc,Nc,lat%num_dirs,lat%volume))
    allocate(G%E(Nc,Nc,lat%num_dirs,lat%volume))

    call cold_start(G)
  end subroutine init_gauge

  subroutine cold_start(G)
    type(gauge_t), intent(inout) :: G
    integer :: v, mu
    do v = 1, G%lat%volume
      do mu = 1, G%lat%num_dirs
        call set_identity(G%U(:,:,mu,v))
        G%E(:,:,mu,v) = (0.0_dp, 0.0_dp)
      end do
    end do
  end subroutine cold_start

  subroutine hot_start(G, seed)
    type(gauge_t), intent(inout) :: G
    integer, intent(in) :: seed
    integer :: v, mu, i, j
    real(dp) :: r
    call random_seed()
    do v = 1, G%lat%volume
      do mu = 1, G%lat%num_dirs
        do i = 1, G%Nc
          do j = 1, G%Nc
            call random_number(r)
            G%U(i,j,mu,v) = cmplx(2.0_dp*r-1.0_dp, 0.0_dp, kind=dp)
          end do
        end do
        call project_to_su(G%U(:,:,mu,v))
        G%E(:,:,mu,v) = (0.0_dp, 0.0_dp)
      end do
    end do
  end subroutine hot_start

  subroutine update_links(G, dt)
    type(gauge_t), intent(inout) :: G
    real(dp), intent(in) :: dt
    integer :: v, mu
    complex(dp), allocatable :: Ustep(:,:)
    allocate(Ustep(G%Nc, G%Nc))

    do v = 1, G%lat%volume
      do mu = 1, G%lat%num_dirs
        call exp_iH(real(G%E(:,:,mu,v)), dt, Ustep)  ! H should be Hermitian; real() applied elementwise to ensure hermitian part
        ! update: U <- exp(i*dt*E) U
        G%U(:,:,mu,v) = matmul(Ustep, G%U(:,:,mu,v))
        call project_to_su(G%U(:,:,mu,v))
      end do
    end do
    deallocate(Ustep)
  end subroutine update_links

  subroutine update_electric_hypercubic(G, g)
    type(gauge_t), intent(inout) :: G
    real(dp), intent(in) :: g
    integer :: v, mu, nu, vp, vpp
    complex(dp), allocatable :: staple(:,:), tmp(:,:)
    integer :: Nc

    if (G%lat%is_honeycomb) then
      stop 'Magnetic force update for honeycomb not implemented.'
    end if

    Nc = G%Nc
    allocate(staple(Nc,Nc), tmp(Nc,Nc))

    do v = 1, G%lat%volume
      do mu = 1, G%lat%num_dirs
        staple = (0.0_dp, 0.0_dp)
        do nu = 1, G%lat%num_dirs
          if (nu == mu) cycle
          ! forward staple: U_mu(x) U_nu(x+mu) U_mu^\dagger(x+nu) U_nu^\dagger(x)
          vp  = neighbor_site(G%lat, v, mu, +1)
          vpp = neighbor_site(G%lat, v, nu, +1)
          tmp = matmul(G%U(:,:,nu,vp), dagger(G%U(:,:,mu,vpp)))
          tmp = matmul(tmp, dagger(G%U(:,:,nu,v)))
          staple = staple + tmp
          ! backward staple: U_mu(x) U_nu^\dagger(x+mu-nu) U_mu^\dagger(x-nu) U_nu(x-nu)
          vp  = neighbor_site(G%lat, v, nu, -1)
          vpp = neighbor_site(G%lat, vp, mu, +1)
          tmp = matmul(dagger(G%U(:,:,nu,vpp)), dagger(G%U(:,:,mu,vp)))
          tmp = matmul(tmp, G%U(:,:,nu,vp))
          staple = staple + tmp
        end do
        ! Force ~ antihermitian traceless( staple * U_mu^\dagger(x) ) / (g*g)
        call antihermitian_traceless(matmul(staple, dagger(G%U(:,:,mu,v))), tmp)
        G%E(:,:,mu,v) = G%E(:,:,mu,v) - (1.0_dp/(g*g)) * tmp
      end do
    end do

    deallocate(staple, tmp)
  end subroutine update_electric_hypercubic

  subroutine measure_energies(G, g, E_electric, E_magnetic)
    type(gauge_t), intent(in) :: G
    real(dp), intent(in) :: g
    real(dp), intent(out) :: E_electric, E_magnetic
    integer :: v, mu, nu
    real(dp) :: tr_real
    complex(dp) :: trc
    integer :: Nc

    Nc = G%Nc
    E_electric = 0.0_dp
    do v = 1, G%lat%volume
      do mu = 1, G%lat%num_dirs
        trc = (0.0_dp, 0.0_dp)
        trc = 0.0_dp
        trc = trace_matrix(matmul(dagger(G%E(:,:,mu,v)), G%E(:,:,mu,v)))
        E_electric = E_electric + real(trc, dp)
      end do
    end do
    E_electric = 0.5_dp * E_electric

    E_magnetic = 0.0_dp
    if (.not. G%lat%is_honeycomb) then
      do v = 1, G%lat%volume
        do mu = 1, G%lat%num_dirs
          do nu = mu+1, G%lat%num_dirs
            ! plaquette U_mu(x) U_nu(x+mu) U_mu^\dagger(x+nu) U_nu^\dagger(x)
            integer :: vp, vpp
            complex(dp) :: Uplaq(G%Nc,G%Nc)
            vp = neighbor_site(G%lat, v, mu, +1)
            vpp = neighbor_site(G%lat, v, nu, +1)
            Uplaq = matmul(G%U(:,:,mu,v), G%U(:,:,nu,vp))
            Uplaq = matmul(Uplaq, dagger(G%U(:,:,mu,vpp)))
            Uplaq = matmul(Uplaq, dagger(G%U(:,:,nu,v)))
            trc = trace_matrix(Uplaq)
            E_magnetic = E_magnetic + (1.0_dp - real(trc,dp)/real(Nc,dp))
          end do
        end do
      end do
      E_magnetic = (1.0_dp/(g*g)) * E_magnetic
    else
      ! Placeholder: honeycomb magnetic energy not implemented
      E_magnetic = -1.0_dp
    end if
  end subroutine measure_energies

  function trace_matrix(A) result(tr)
    complex(dp), intent(in) :: A(:,:)
    complex(dp) :: tr
    integer :: i, n
    n = size(A,1)
    tr = (0.0_dp, 0.0_dp)
    do i = 1, n
      tr = tr + A(i,i)
    end do
  end function trace_matrix

  subroutine step_leapfrog(G, g, dt)
    type(gauge_t), intent(inout) :: G
    real(dp), intent(in) :: g, dt
    ! Half-step E
    call update_electric_hypercubic(G, g)
    G%E = G%E * (dt * 0.5_dp)
    ! Full-step U
    call update_links(G, dt)
    ! Half-step E
    call update_electric_hypercubic(G, g)
    G%E = G%E * (dt * 0.5_dp)
  end subroutine step_leapfrog

end module gauge_field
