!> Module testing our gaussian integration routine
module test_gauss
use fruit
use gauss
implicit none
contains

!> Test $\int_0^1 x^n dx$, integrated from 0 to 1
subroutine test_gauss_p
  integer :: i, j
  real*8 :: a
  character(len=2) :: n
  character(len=20) :: s
  real*8 :: tol
  
  do j=1,14
    a = 0.d0
    do i=1,n_gauss
      a = a + Xgauss(i)**j*Wgauss(i)
    end do

    if (j .lt. 2*n_gauss) then
      tol = n_gauss/4d-17 ! integration is exact up to 2*n_gauss. allow for some floating-point errors
    else
      if (n_gauss .eq. 4) then
        tol = 2d-3 ! below that it is not so accurate (this value is for n_gauss=4)
      else
        tol = 0.d0 ! for testing other schemes
      end if
    end if

    write(n,'(i2)') j
    write(s,'(g20.12)') 1.d0/real(j+1)
    call assert_equals(1.d0/real(j+1), a, tol, 'integral x^'//trim(adjustl(n))//' must be '//s)
  end do
end subroutine test_gauss_p
end module test_gauss
