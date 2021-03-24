!> This module contains some testcases for the sobol series implementation
module sobseq_spec
use mod_sobseq_rng
use fruit
implicit none

contains

!> Test whether the minimum and maximum values are between 0 and 1
subroutine test_sobseq_minmax
  type(sobseq_rng) :: rng
  integer :: ifail, i
  integer, parameter :: n_tries = 1000
  real*8 :: x(n_tries)
  call rng%initialize(n_dims=1, seed=0, n_streams=1, i_stream=1, ierr=ifail)
  call assert_equals(0, ifail)
  do i=1,n_tries
    call rng%next(x(i:i))
  end do
  call assert_true(1.d0 > maxval(x), "no values can be above 1")
  call assert_true(0.d0 < minval(x), "no values can be below 0")
end subroutine test_sobseq_minmax

!> Test the mean of the series in sequential mode
subroutine test_sobseq_mean
  call assert_equals(0d0, sobseq_mean(1, 100), 1d-2, "1 stream 100 points")
  call assert_equals(0d0, sobseq_mean(1, 1000), 1d-3, "1 stream 1000 points")
  call assert_equals(0d0, sobseq_mean(1, 10000), 1d-4, "1 stream 10000 points")
  call assert_equals(0d0, sobseq_mean(2, 100), 1d-2, "2 streams 100 points")
  call assert_equals(0d0, sobseq_mean(2, 1000), 1d-3, "2 streams 1000 points")
  call assert_equals(0d0, sobseq_mean(2, 10000), 1d-4, "2 streams 10000 points")
  call assert_equals(0d0, sobseq_mean(8, 10000), 1.03d-4, "8 streams 10000 points")
end subroutine test_sobseq_mean

!> Calculate error in mean of a uniform sobol sequence in 1d with n_streams and n_points
function sobseq_mean(n_streams, n_points) result(mean)
  type(sobseq_rng) :: rng(n_streams)
  integer :: ifail, i, n_streams, n_points
  real*8 :: x(n_points), mean
  do i=1,n_streams
    call rng(i)%initialize(n_dims=1, seed=0, n_streams=n_streams, i_stream=i, ierr=ifail)
  end do
  do i=1,n_points
    call rng(mod(i,n_streams)+1)%next(x(i:i))
  end do
  mean = (sum(x)/n_points) - 0.5d0
end function sobseq_mean

subroutine test_sobseq_errors
  type(sobseq_rng) :: rng
  integer :: ifail
  call rng%initialize(n_dims=0, seed=0, n_streams=1, i_stream=1, ierr=ifail)
  call assert_equals(1, ifail, "no dims produces error 1")
  call rng%initialize(n_dims=1000, seed=0, n_streams=1, i_stream=1, ierr=ifail)
  call assert_equals(6, ifail, "too many dims produces error 6")
  call rng%initialize(n_dims=1, seed=0, n_streams=0, i_stream=0, ierr=ifail)
  call assert_equals(3, ifail, "too few streams")
  call rng%initialize(n_dims=1, seed=0, n_streams=1, i_stream=0, ierr=ifail)
  call assert_equals(4, ifail, "i_stream < 1")
  call rng%initialize(n_dims=1, seed=0, n_streams=1, i_stream=2, ierr=ifail)
  call assert_equals(5, ifail, "i_stream > n_streams")
end subroutine test_sobseq_errors

subroutine test_sobseq_known_values
  type(sobseq_rng) :: rng, rng_strided(2)
  integer :: ifail, i, j
  real*8 :: x(8)
  real*8, dimension(10,12) :: known_values !< taken from http://michaelcarteronline.com/MCM/LDSequences/SobolExample.pdf
  character(len=5) :: s
  known_values(1,:) = 1.d0/2.d0
  known_values(2,:) = real([1,1,1,3,3,1,3,3,3,3,3,1],8)/4.d0
  known_values(3,:) = real([3,3,3,1,1,3,1,1,1,1,1,3],8)/4.d0
  known_values(4,:) = real([3,5,7,3,1,3,7,7,5,7,3,3],8)/8.d0
  known_values(5,:) = real([7,1,3,7,5,7,3,3,1,3,7,7],8)/8.d0
  known_values(6,:) = real([1,7,5,5,7,1,1,1,3,1,5,1],8)/8.d0
  known_values(7,:) = real([5,3,1,1,3,5,5,5,7,5,1,5],8)/8.d0
  known_values(8,:) = real([5,15,7,9,5,7,15,15,5,11,1,15],8)/16.d0
  known_values(9,:) = real([13,7,15,1,13,15,7,7,13,3,9,7],8)/16.d0
  call rng%initialize(n_dims=8, seed=0, n_streams=1, i_stream=1, ierr=ifail)
  call assert_equals(0, ifail, "must seed normal without error")
  call rng_strided(1)%initialize(n_dims=8, seed=0, n_streams=2, i_stream=1, ierr=ifail)
  call assert_equals(0, ifail, "must seed strided 1 without error")
  call rng_strided(2)%initialize(n_dims=8, seed=0, n_streams=2, i_stream=2, ierr=ifail)
  call assert_equals(0, ifail, "must seed strided 2 without error")
  ! rng actually starts at position 1
  do i=2,9
    call rng%next(x)
    do j=1,8
      write(s,"(A1,I1,A1,I1,A1)") "(", i, ",", j, ")"
      call assert_equals(known_values(i,j), x(j), s)
    end do
  end do
  do i=3,8,2
    call rng_strided(1)%next(x)
    do j=1,8
      write(s,"(A1,I1,A1,I1,A1)") "(", i, ",", j, ")"
      call assert_equals(known_values(i,j), x(j), s)
    end do
    call rng_strided(2)%next(x)
    do j=1,8
      write(s,"(A1,I1,A1,I1,A1)") "(", i+1, ",", j, ")"
      call assert_equals(known_values(i+1,j), x(j), s)
    end do
  end do
end subroutine test_sobseq_known_values

subroutine test_sobseq_1D
  call sobseq_1D_integral(10)
  call sobseq_1D_integral(100)
  call sobseq_1D_integral(1000)
end subroutine test_sobseq_1D
subroutine test_sobseq_1D_stride_4
  call sobseq_1D_integral_strided(10, 4)
  call sobseq_1D_integral_strided(100, 4)
  call sobseq_1D_integral_strided(1000, 4)
end subroutine test_sobseq_1D_stride_4
subroutine test_sobseq_1D_stride_8
  call sobseq_1D_integral_strided(10, 8)
  call sobseq_1D_integral_strided(100, 8)
  call sobseq_1D_integral_strided(1000, 8)
end subroutine test_sobseq_1D_stride_8
subroutine test_sobseq_1D_stride_16
  call sobseq_1D_integral_strided(100, 16)
  call sobseq_1D_integral_strided(1000, 16)
  call sobseq_1D_integral_strided(10000, 16)
end subroutine test_sobseq_1D_stride_16
subroutine test_sobseq_2D
  call sobseq_2D_integral(10)
  call sobseq_2D_integral(100)
  call sobseq_2D_integral(1000)
end subroutine test_sobseq_2D
subroutine test_sobseq_1D_stride_3_error
  type(sobseq_rng) :: rng
  integer :: ifail
  call rng%initialize(n_dims=1, seed=0, n_streams=3, i_stream=1, ierr=ifail)
  call assert_true(ifail .gt. 0, "using n_streams .ne. 2^n is an error")
end subroutine test_sobseq_1D_stride_3_error


subroutine test_sobseq_1D_strided_equivalence
  type(sobseq_rng) :: rng, rng_strided(2)
  integer :: ifail
  real*8, dimension(1) :: x, x_strided

  call            rng%initialize(n_dims=1, seed=0, n_streams=1, i_stream=1, ierr=ifail)
  call assert_equals(0, ifail)
  call rng_strided(1)%initialize(n_dims=1, seed=0, n_streams=2, i_stream=1, ierr=ifail)
  call assert_equals(0, ifail)
  call rng_strided(2)%initialize(n_dims=1, seed=0, n_streams=2, i_stream=2, ierr=ifail)
  call assert_equals(0, ifail)

  ! one skip to get synced
  call rng%next(x)

  call rng%next(x)
  call rng_strided(1)%next(x_strided)
  call assert_equals(x(1), x_strided(1), "element 1 should be same")
  call rng%next(x)
  call rng_strided(2)%next(x_strided)
  call assert_equals(x(1), x_strided(1), "element 2 should be same")
  call rng%next(x)
  call rng_strided(1)%next(x_strided)
  call assert_equals(x(1), x_strided(1), "element 3 should be same")
  call rng%next(x)
  call rng_strided(2)%next(x_strided)
  call assert_equals(x(1), x_strided(1), "element 4 should be same")
  call rng%next(x)
  call rng_strided(1)%next(x_strided)
  call assert_equals(x(1), x_strided(1), "element 5 should be same")
  call rng%next(x)
  call rng_strided(2)%next(x_strided)
  call assert_equals(x(1), x_strided(1), "element 6 should be same")
end subroutine test_sobseq_1D_strided_equivalence

subroutine test_ilog2_b_ceil
  call assert_equals(-1, ilog2_b_ceil(-1), "negative values not allowed, return -1")
  call assert_equals(-1, ilog2_b_ceil(0), "negative values not allowed, return -1")
  call assert_equals(0, ilog2_b_ceil(1), "2^0 >= 1")
  call assert_equals(1, ilog2_b_ceil(2), "2^1 >= 2")
  call assert_equals(2, ilog2_b_ceil(3), "2^2 >= 3")
  call assert_equals(2, ilog2_b_ceil(4), "2^2 >= 4")
  call assert_equals(3, ilog2_b_ceil(5), "2^3 >= 5")
  call assert_equals(3, ilog2_b_ceil(6), "2^3 >= 6")
  call assert_equals(3, ilog2_b_ceil(7), "2^3 >= 7")
  call assert_equals(3, ilog2_b_ceil(8), "2^3 >= 8")
end subroutine test_ilog2_b_ceil


!> Test the sobol series by integrating \(8r^3-6r+2\) over 0 to 1 (=1)
subroutine sobseq_1D_integral(n)
  type(sobseq_rng) :: rng
  integer, intent(in) :: n
  integer :: ifail, i
  real*8 :: x(1), integral(1), tolerance
  character(len=6) :: n_s

  call rng%initialize(n_dims=1, seed=0, n_streams=1, i_stream=1, ierr=ifail)
  call assert_equals(0, ifail)
  if (ifail .gt. 0) return

  x(1) = 0.d0
  integral = 0.d0
  do i=1,n
    call rng%next(x)
    integral = integral + 8.d0*(x**3) - 6*x + 2
  end do
  integral = integral / real(n,8)

  tolerance = 2.1d0/real(n,8)
  write(n_s,"(i6)") n
  call assert_equals(1.d0, integral(1), tolerance, "integral 8r^3-6r+2 = 1" // n_s)
end subroutine sobseq_1D_integral

!> Test the sobol series by integrating \(8r^3-6r+2\) over 0 to 1 (=1)
subroutine sobseq_1D_integral_strided(n, num_streams)
  integer, intent(in) :: n, num_streams
  type(sobseq_rng) :: rng(num_streams)
  integer :: ifail, i, j
  real*8 :: x(1), integral(1), tolerance
  character(len=6) :: n_s

  do i=1,num_streams
    call rng(i)%initialize(n_dims=1, seed=0, n_streams=num_streams, i_stream=i, ierr=ifail)
    call assert_equals(0, ifail)
    if (ifail .gt. 0) return
  end do

  x(1) = 0.d0
  integral = 0.d0
  do i=1,n/num_streams
    do j=1,num_streams
      call rng(j)%next(x)
      integral = integral + 8.d0*(x**3) - 6*x + 2
    end do
  end do
  integral = integral / real(n,8)

  tolerance = 7d0/real(n,8) ! could be 2.1 for strides 1-4
  write(n_s,"(i6)") n
  call assert_equals(1.d0, integral(1), tolerance, "integral 8r^3-6r+2 = 1" // n_s)
end subroutine sobseq_1D_integral_strided

!> Test the sobol series by integrating \(1-|x|-|y|\) from 0 to 1 (result = 0)
subroutine sobseq_2D_integral(n)
  type(sobseq_rng) :: rng
  integer, intent(in) :: n
  integer :: ifail, i
  real*8 :: x(2), integral(2), tolerance
  character(len=6) :: n_s

  call rng%initialize(n_dims=2, seed=0, n_streams=1, i_stream=1, ierr=ifail)
  call assert_equals(0, ifail)

  integral = 0.d0
  do i=1,n
    call rng%next(x)
    integral = integral + (1.d0 - abs(x(1)) - abs(x(2)))
  end do
  integral = integral / real(n,8)

  write(n_s,"(i6)") n
  tolerance = 1.d0/real(n,8)
  call assert_equals(0.d0, integral(1), tolerance, "integral 1-|x|-|y| = 0")
end subroutine sobseq_2D_integral
end module sobseq_spec
