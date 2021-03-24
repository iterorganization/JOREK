!> Test root-finding algorithms
module rootfinding_spec
use mod_rootfinding
use fruit
implicit none

real*8, parameter :: tolerance = 1d-8

contains

!> Test whether the minimum and maximum values are between 0 and 1
subroutine test_newtons_method_tanh
  integer :: i, ierr
  real*8  :: y, x, ref
  do i=2,9
    y   = real((i-1)/10.d0,8) ! solve tanh(x) = y
    ref = atanh(y) ! reference solution
    ! guess at ref-0.08, ref+0.08
    call newtons_method(f, df, y, ref-0.08, x, ierr)
    call assert_equals(ierr, 0, "must find a solution")
    call assert_equals(ref, x, tolerance, "must find the right solution (1)")
    call newtons_method(f, df, y, ref+0.08, x, ierr)
    call assert_equals(ierr, 0, "must find a solution")
    call assert_equals(ref, x, tolerance, "must find the right solution (2)")
  end do
end subroutine test_newtons_method_tanh

!> Test whether the minimum and maximum values are between 0 and 1
subroutine test_halleys_method_sqrt
  integer :: i, ierr
  real*8  :: y, x, ref
  do i=2,9
    y   = real((i-1)/10.d0,8) ! solve sqrt(x) = y
    ref = atanh(y) ! reference solution
    ! guess at ref-0.08, ref+0.08
    call halleys_method(f, df, ddf, y, ref-0.08, x, ierr)
    call assert_equals(ierr, 0, "must find a solution")
    call assert_equals(ref, x, tolerance, "must find the right solution (1)")
    call halleys_method(f, df, ddf, y, ref+0.08, x, ierr)
    call assert_equals(ierr, 0, "must find a solution")
    call assert_equals(ref, x, tolerance, "must find the right solution (2)")
  end do
end subroutine test_halleys_method_sqrt

pure function f(x)
  real*8, intent(in) :: x
  real*8 :: f
  f = tanh(x)
end function f

pure function df(x)
  real*8, intent(in) :: x
  real*8 :: df
  df = 1.d0 - tanh(x)*tanh(x)
end function df

pure function ddf(x)
  real*8, intent(in) :: x
  real*8 :: ddf
  ddf = -2.d0*tanh(x)/cosh(x)/cosh(x)
end function ddf
end module rootfinding_spec
