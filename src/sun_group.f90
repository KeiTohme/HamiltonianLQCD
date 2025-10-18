module sun_group
  use kinds
  implicit none
  private
  public :: set_identity, dagger, antihermitian_traceless, project_to_su, exp_iH

contains

  subroutine set_identity(U)
    complex(dp), intent(inout) :: U(:,:)
    integer :: i, n
    n = size(U,1)
    U = (0.0_dp, 0.0_dp)
    do i = 1, n
      U(i,i) = (1.0_dp, 0.0_dp)
    end do
  end subroutine set_identity

  pure function dagger(A) result(B)
    complex(dp), intent(in) :: A(:,:)
    complex(dp) :: B(size(A,2), size(A,1))
    B = transpose(conjg(A))
  end function dagger

  subroutine antihermitian_traceless(X, Y)
    complex(dp), intent(in) :: X(:,:)
    complex(dp), intent(out) :: Y(size(X,1), size(X,2))
    complex(dp) :: A(size(X,1), size(X,2))
    complex(dp) :: tr
    integer :: n, i
    n = size(X,1)
    A = 0.5_dp * (X - dagger(X))
    tr = (0.0_dp, 0.0_dp)
    do i = 1, n
      tr = tr + A(i,i)
    end do
    A = A - (tr / real(n,dp)) * identity_matrix(n)
    Y = A
  end subroutine antihermitian_traceless

  function identity_matrix(n) result(I)
    integer, intent(in) :: n
    complex(dp) :: I(n,n)
    integer :: i
    I = (0.0_dp, 0.0_dp)
    do i = 1, n
      I(i,i) = (1.0_dp, 0.0_dp)
    end do
  end function identity_matrix

  subroutine project_to_su(U)
    complex(dp), intent(inout) :: U(:,:)
    integer :: n, i, j, k
    complex(dp) :: Q(size(U,1), size(U,2))
    complex(dp) :: v(size(U,1)), uvec(size(U,1))
    complex(dp) :: detU, phase
    real(dp) :: norm

    n = size(U,1)
    Q = (0.0_dp, 0.0_dp)

    do j = 1, n
      v = U(:,j)
      do k = 1, j-1
        v = v - dot_product(conjg(Q(:,k)), U(:,j)) * Q(:,k)
      end do
      norm = sqrt(real(dot_product(conjg(v), v), dp))
      if (norm > 0.0_dp) then
        Q(:,j) = v / norm
      else
        Q(:,j) = (0.0_dp, 0.0_dp)
        Q(j,j) = (1.0_dp, 0.0_dp)
      end if
    end do

    U = Q

    detU = determinant(U)
    if (abs(detU) > 0.0_dp) then
      phase = detU / abs(detU)
    else
      phase = (1.0_dp, 0.0_dp)
    end if
    phase = phase ** (1.0_dp/real(n,dp))
    do i = 1, n
      U(i,i) = U(i,i) / phase
    end do
  end subroutine project_to_su

  function determinant(A) result(det)
    complex(dp), intent(in) :: A(:,:)
    complex(dp) :: det
    complex(dp), allocatable :: M(:,:)
    integer :: n, i, j, k, piv
    complex(dp) :: factor
    real(dp) :: maxabs
    integer :: maxi

    n = size(A,1)
    allocate(M(n,n))
    M = A
    det = (1.0_dp, 0.0_dp)

    do k = 1, n
      ! partial pivoting
      maxabs = 0.0_dp; maxi = k
      do i = k, n
        if (abs(M(i,k)) > maxabs) then
          maxabs = abs(M(i,k)); maxi = i
        end if
      end do
      if (maxabs == 0.0_dp) then
        det = (0.0_dp, 0.0_dp)
        return
      end if
      if (maxi /= k) then
        M([k,maxi],:) = M([maxi,k],:)
        det = -det
      end if
      det = det * M(k,k)
      do i = k+1, n
        factor = M(i,k) / M(k,k)
        M(i,k:n) = M(i,k:n) - factor * M(k,k:n)
      end do
    end do

    deallocate(M)
  end function determinant

  subroutine matmul_nn(A, B, C)
    complex(dp), intent(in) :: A(:,:), B(:,:)
    complex(dp), intent(out) :: C(size(A,1), size(B,2))
    integer :: i,j,k, n, m, p
    n = size(A,1); m = size(A,2); p = size(B,2)
    C = (0.0_dp, 0.0_dp)
    do i = 1, n
      do k = 1, m
        do j = 1, p
          C(i,j) = C(i,j) + A(i,k) * B(k,j)
        end do
      end do
    end do
  end subroutine matmul_nn

  subroutine exp_iH(H, a, Uexp)
    complex(dp), intent(in) :: H(:,:)
    real(dp), intent(in) :: a
    complex(dp), intent(out) :: Uexp(size(H,1), size(H,2))
    integer :: n, k
    complex(dp), allocatable :: X(:,:), term(:,:)
    real(dp) :: tol

    n = size(H,1)
    allocate(X(n,n), term(n,n))
    ! Compute X = i * a * H
    X = (0.0_dp, 0.0_dp)
    X = (0.0_dp, 1.0_dp) * a * H
    Uexp = identity_matrix(n)
    term = identity_matrix(n)

    tol = 1.0e-12_dp
    do k = 1, 24
      call matmul_nn(term, X, term)
      term = term / real(k, dp)
      Uexp = Uexp + term
      if (maxval(abs(term)) < tol) exit
    end do
    call project_to_su(Uexp)
    deallocate(X, term)
  end subroutine exp_iH

end module sun_group
