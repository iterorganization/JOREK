!> Module testing our quicksort implementation
module test_quicksort
use fruit
use mod_quicksort
implicit none
contains
subroutine test_quicksort_1
  real*8 :: a(1)
  a(1) = 1.d0
  call quicksort(a)
  call assert_equals(1.d0, a(1), 'does nothing')
end subroutine test_quicksort_1

subroutine test_quicksort_100
  integer, parameter :: n = 100
  real*8 :: a(n)
  integer :: i
  do i=1,n
    a(i) = real((n-i)**2)
  end do
  call quicksort(a)
  do i=1,n
    call assert_equals(real((i-1)**2), a(i), 'element i right')
  end do
end subroutine test_quicksort_100
end module test_quicksort
