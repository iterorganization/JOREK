!> This module contains some testcases for the sobol series implementation
module sample_gaussian_spec
use mod_pcg32_rng
use mod_sampling
use fruit
implicit none

contains

!> Test whether the minimum and maximum values are between 0 and 1
subroutine test_gaussian_mean
  type(pcg32_rng) :: rng
  integer :: ifail, i
  integer, parameter :: n_tries = 1000
  real*8 :: x(n_tries)
  call rng%initialize(n_dims=1, seed=1231789264, n_streams=1, i_stream=1, ierr=ifail)
  call assert_equals(0, ifail)
  do i=1,n_tries
    call rng%next(x(i:i))
    x(i) = sample_gaussian(x(i))
  end do

  call assert_equals(sum(x)/n_tries, 0.d0, 7d-2, "Mean should be 0")
end subroutine test_gaussian_mean


!> Test whether the minimum and maximum values are between 0 and 1
subroutine test_gaussian_variance
  type(pcg32_rng) :: rng
  integer :: ifail, i
  integer, parameter :: n_tries = 1000
  real*8 :: x(n_tries)
  call rng%initialize(n_dims=1, seed=1231789264, n_streams=1, i_stream=1, ierr=ifail)
  call assert_equals(0, ifail)
  do i=1,n_tries
    call rng%next(x(i:i))
    x(i) = sample_gaussian(x(i))
  end do

  call assert_equals(dot_product(x,x)/n_tries, 1.d0, 7d-2, "Variance should be 0")
end subroutine test_gaussian_variance
end module sample_gaussian_spec
