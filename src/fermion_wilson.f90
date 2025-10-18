module fermion_wilson
  use kinds
  use lattice_geometry
  implicit none
  private
  public :: apply_wilson_dirac
contains

  subroutine apply_wilson_dirac(psi_in, psi_out, U, lat, m)
    complex(dp), intent(in) :: psi_in(:,:)     ! (Nc*Ns, volume)
    complex(dp), intent(out) :: psi_out(:,:)   ! (Nc*Ns, volume)
    complex(dp), intent(in) :: U(:,:,:,:)      ! (Nc,Nc,ndir,volume)
    type(lattice_t), intent(in) :: lat
    real(dp), intent(in) :: m
    integer :: v, mu, Nc, Ns, vol

    Nc = size(U,1)
    vol = size(U,4)
    Ns = size(psi_in,1) / Nc

    psi_out = (m + 4.0_dp) * psi_in  ! naive mass + Wilson r=1 diagonal

    if (lat%is_honeycomb) then
      stop 'Wilson fermion operator not implemented for honeycomb.'
    end if

    do v = 1, vol
      do mu = 1, lat%num_dirs
        integer :: vp
        complex(dp) :: psi_v(:), psi_p(:)
        complex(dp), allocatable :: tmp(:)
        vp = neighbor_site(lat, v, mu, +1)
        psi_v => psi_in(:, v)
        psi_p => psi_in(:, vp)
        allocate(tmp(Nc*Ns))
        tmp = psi_p
        ! gauge transport on color indices only; spin left untouched (Ns=4)
        call apply_link_on_color(tmp, U(:,:,mu,v), Nc, Ns)
        psi_out(:,v) = psi_out(:,v) - 0.5_dp * ( (1.0_dp,0.0_dp) * tmp )
        deallocate(tmp)

        ! backward hop
        vp = neighbor_site(lat, v, mu, -1)
        psi_p => psi_in(:, vp)
        allocate(tmp(Nc*Ns))
        tmp = psi_p
        call apply_link_dagger_on_color(tmp, U(:,:,mu,vp), Nc, Ns)
        psi_out(:,v) = psi_out(:,v) - 0.5_dp * ( (1.0_dp,0.0_dp) * tmp )
        deallocate(tmp)
      end do
    end do
  end subroutine apply_wilson_dirac

  subroutine apply_link_on_color(psi, Ulink, Nc, Ns)
    complex(dp), intent(inout) :: psi(:)     ! Nc*Ns vector
    complex(dp), intent(in) :: Ulink(:,:)
    integer, intent(in) :: Nc, Ns
    complex(dp), allocatable :: tmp(:)
    integer :: s, i

    allocate(tmp(Nc))
    do s = 1, Ns
      tmp = matmul(Ulink, psi((s-1)*Nc+1:s*Nc))
      psi((s-1)*Nc+1:s*Nc) = tmp
    end do
    deallocate(tmp)
  end subroutine apply_link_on_color

  subroutine apply_link_dagger_on_color(psi, Ulink, Nc, Ns)
    complex(dp), intent(inout) :: psi(:)
    complex(dp), intent(in) :: Ulink(:,:)
    integer, intent(in) :: Nc, Ns
    complex(dp), allocatable :: tmp(:)
    integer :: s

    allocate(tmp(Nc))
    do s = 1, Ns
      tmp = matmul(transpose(conjg(Ulink)), psi((s-1)*Nc+1:s*Nc))
      psi((s-1)*Nc+1:s*Nc) = tmp
    end do
    deallocate(tmp)
  end subroutine apply_link_dagger_on_color

end module fermion_wilson
