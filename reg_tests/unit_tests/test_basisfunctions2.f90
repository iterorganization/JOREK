!> Module testing accuracy of calculating basis functions with 2 different routines
module test_basisfunctions2
use fruit
implicit none

real*8, parameter :: tol = 1d-14
contains

!> Check that the first-derivative and second-derivative versions of the basisfunctions match at the gaussian points
subroutine test_basisfunctions_2D_1_basisfunctions_2D_2
  use gauss
  use mod_basisfunctions
  use mod_parameters, only: n_vertex_max, n_order

  integer :: k, l, m, n
  real*8  :: s, t
  real*8, dimension(n_vertex_max,n_order+1) :: H, H_s, H_t, H_st, H_ss, H_tt, G, G_s, G_t, G_st, F, F_s, F_t

  ! Verify against other expression
  do k=1,n_gauss
    s = xgauss(k)
    do l=1,n_gauss
      t = xgauss(l)
      call basisfunctions(s,t,H, H_s, H_t, H_st, H_ss, H_tt)
      call basisfunctions(s,t,G, G_s, G_t, G_st)
      call basisfunctions(s,t,F, F_s, F_t)

      do m=1,n_vertex_max
        do n=1,n_order+1
          call assert_equals(H(m,n), G(m,n), tol, 'value match')
          call assert_equals(H(m,n), F(m,n), tol, 'value match')
          call assert_equals(H_s(m,n), G_s(m,n), tol, 's-derivative match')
          call assert_equals(H_s(m,n), F_s(m,n), tol, 's-derivative match')
          call assert_equals(H_t(m,n), G_t(m,n), tol, 't-derivative match')
          call assert_equals(H_t(m,n), F_t(m,n), tol, 't-derivative match')
          call assert_equals(H_st(m,n), G_st(m,n), tol, 'st-derivative match')
        enddo
      enddo
    enddo
  enddo
end subroutine test_basisfunctions_2D_1_basisfunctions_2D_2

!> Check that the transposed version of the basisfunctions is equal to the transpose of the normal versino
subroutine test_basisfunctions_2D_vs_transpose
  use gauss
  use mod_basisfunctions
  use mod_parameters, only: n_vertex_max, n_order

  integer :: k, l, m, n
  real*8  :: s, t
  real*8, dimension(n_vertex_max,n_order+1) :: H, H_s, H_t, G, G_s, G_t

  do k=1,n_gauss
    s = xgauss(k)
    do l=1,n_gauss
      t = xgauss(l)
      call basisfunctions(s,t,H, H_s, H_t)
      call basisfunctions_T(s,t,G, G_s, G_t)

      do m=1,n_vertex_max
        do n=1,n_order+1
          call assert_equals(H(m,n), G(n,m), tol, 'value match')
          call assert_equals(H_s(m,n), G_s(n,m), tol, 's-derivative match')
          call assert_equals(H_t(m,n), G_t(n,m), tol, 't-derivative match')
        enddo
      enddo
    enddo
  enddo
end subroutine test_basisfunctions_2D_vs_transpose


!> Test if the basisfunctions are 1 on their node
subroutine test_basisfunction_properties
  use mod_basisfunctions
  use mod_parameters, only: n_vertex_max, n_order

  integer :: j, k
  real*8  :: s, t
  real*8, dimension(n_vertex_max,n_order+1) :: H

  do k=1,n_vertex_max
  
    select case (k)
    case (1)
      s = 0
      t = 0
    case (2)
      s = 1
      t = 0
    case (3)
      s = 1
      t = 1
    case (4)
      s = 0
      t = 1
    end select
    call basisfunctions(s,t,H)

    call assert_equals(1.d0,H(k,1), tol, 'value 1 at right node')
    do j=0,2
      call assert_equals(0.d0,H(mod(k+j,4)+1,1), tol, 'value 0 at other node')
    enddo
  enddo
end subroutine test_basisfunction_properties
end module test_basisfunctions2
