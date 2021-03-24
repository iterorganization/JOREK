!> Module testing hermite_birkhoff interpolation
module test_hermite_birkhoff
use fruit
use mod_hermite_birkhoff
contains

subroutine test_interp_hermite_birkhoff_cos
  integer, parameter   :: N_tests = 10
  real*8, parameter :: t0 = 0.d0, t1 = 1.d0
  real*8, parameter :: ftol = 3d-3, dftol = 1d-2
  real*8, dimension(1) :: y0, y1, y, dy0, dy1
  integer :: i
  real*8 :: t
  y0 = cos(t0)
  y1 = cos(t1)
  dy0 = -sin(t0)
  dy1 = -sin(t1)

  do i=1,N_tests
    t = real(i-1)/real(N_tests-1) ! include endpoints
    call HB_interp(t0, t1, 1, y0, y1, dy0, dy1, t, y)
    call assert_equals(cos(t), y(1), ftol, "Interp function test")
    call HB_interp_dt(t0, t1, 1, y0, y1, dy0, dy1, t, y)
    call assert_equals(-sin(t), y(1), dftol, "Interp derivative test")
  end do
end subroutine test_interp_hermite_birkhoff_cos
end module test_hermite_birkhoff
