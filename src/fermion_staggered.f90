module fermion_staggered
  use kinds
  use lattice_geometry
  implicit none
  private
  public :: apply_staggered
contains

  subroutine apply_staggered(chi_in, chi_out, U, lat, m)
    complex(dp), intent(in) :: chi_in(:,:)     ! (Nc, volume)
    complex(dp), intent(out) :: chi_out(:,:)   ! (Nc, volume)
    complex(dp), intent(in) :: U(:,:,:,:)      ! (Nc,Nc,ndir,volume)
    type(lattice_t), intent(in) :: lat
    real(dp), intent(in) :: m
    integer :: v, mu, Nc, vol

    Nc = size(U,1)
    vol = size(U,4)

    chi_out = m * chi_in

    do v = 1, vol
      do mu = 1, lat%num_dirs
        integer :: vp
        complex(dp) :: phase
        complex(dp), allocatable :: tmp(:)
        vp = neighbor_site(lat, v, mu, +1)
        phase = eta_phase(lat, v, mu)
        allocate(tmp(Nc))
        tmp = matmul(U(:,:,mu,v), chi_in(:,vp))
        chi_out(:,v) = chi_out(:,v) + 0.5_dp * phase * tmp
        deallocate(tmp)

        vp = neighbor_site(lat, v, mu, -1)
        phase = conjg(eta_phase(lat, vp, mu))
        allocate(tmp(Nc))
        tmp = matmul(transpose(conjg(U(:,:,mu,vp))), chi_in(:,vp))
        chi_out(:,v) = chi_out(:,v) - 0.5_dp * phase * tmp
        deallocate(tmp)
      end do
    end do
  end subroutine apply_staggered

  complex(dp) function eta_phase(lat, v, mu) result(eta)
    type(lattice_t), intent(in) :: lat
    integer, intent(in) :: v, mu
    integer :: coord(4), s, i
    integer :: subl

    if (lat%is_honeycomb) then
      ! Simple choice: phases all 1 on honeycomb as placeholder
      eta = (1.0_dp, 0.0_dp)
      return
    end if

    ! Hypercubic: eta_mu(x) = (-1)^{sum_{nu<mu} x_nu}
    call index_to_coord(lat, v, coord(1:lat%ndim))
    s = 0
    do i = 1, mu-1
      s = s + coord(i)
    end do
    if (mod(s,2) == 0) then
      eta = (1.0_dp, 0.0_dp)
    else
      eta = (-1.0_dp, 0.0_dp)
    end if
  end function eta_phase

end module fermion_staggered
