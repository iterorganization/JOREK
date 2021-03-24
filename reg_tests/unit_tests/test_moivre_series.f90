!> Module testing accuracy of calculating sines by applying De Moivre's formula
module test_moivre_series
use fruit
use mod_interp
use mod_parameters, only: n_tor, n_period
implicit none
contains

subroutine test_moivre_sincos
  integer, parameter :: N_max = 40, N_phi = 4
  real*8, parameter :: tol = 1d-14
  real*8, parameter :: phi0(N_phi) = [0.1d0, 0.5d0, 1.d0, 12.d0]
  character(len=3) :: i_s
  integer :: i, j
  real*8 :: e1r, e1i, enr, eni, enrtmp

  do i=1,N_phi
    e1r = cos(phi0(i))
    e1i = sin(phi0(i))
    enr = 1.d0
    eni = 0.d0
    do j=1,N_max
      enrtmp = enr*e1r - eni*e1i
      eni = eni*e1r + enr*e1i
      enr = enrtmp
      write(i_s,'(i3)') j
      call assert_equals(cos(real(j)*phi0(i)), enr, tol, 'Cosine match '//i_s)
      call assert_equals(sin(real(j)*phi0(i)), eni, tol, 'Sine match '//i_s)
    end do
  end do
end subroutine test_moivre_sincos

! Note that these tests should really be run with n_tor > 1 (and maybe n_period > 1)
subroutine test_mode_moivre
  integer, parameter :: N_max = 40, N_phi = 4
  real*8, parameter :: tol = 1d-14
  real*8, parameter :: phi0(N_phi) = [0.1d0, 0.5d0, 1.d0, 12.d0]
  character(len=3) :: i_s
  integer :: i, j
  real*8 :: HZ(n_tor)

  do j=1,N_phi
    call mode_moivre(phi0(j), HZ)
    call assert_equals(1.d0, HZ(1), tol, 'static match')
    do i=1,(n_tor-1)/2
      write(i_s,'(i3)') i
      call assert_equals(cos(n_period*i*phi0(j)), HZ(2*i), tol, 'Cosine match '//i_s)
      call assert_equals(sin(n_period*i*phi0(j)), HZ(2*i+1), tol, 'Sine match '//i_s)
    end do
  end do
end subroutine test_mode_moivre

subroutine test_sincosperiod_moivre
  integer, parameter :: N_max = 40, N_phi = 4
  real*8, parameter :: tol = 1d-14
  real*8, parameter :: phi0(N_phi) = [0.1d0, 0.5d0, 1.d0, 12.d0]
  character(len=3) :: i_s
  integer :: i, j
  real*8 :: HZ(n_tor), dHZ(n_tor)

  do j=1,N_phi
    call sincosperiod_moivre(phi0(j), HZ, dHZ)
    call assert_equals(1.d0, HZ(1), tol, 'static match')
    call assert_equals(0.d0, dHZ(1), tol, 'static match')
    do i=1,(n_tor-1)/2
      write(i_s,'(i3)') i
      call assert_equals(cos(real(i)*phi0(j)), HZ(2*i), tol, 'Cosine match '//i_s)
      call assert_equals(sin(real(i)*phi0(j)), HZ(2*i+1), tol, 'Sine match '//i_s)

      call assert_equals(-n_period*i*sin(n_period*i*phi0(j)), dHZ(2*i), tol, 'd Cosine match '//i_s)
      call assert_equals(n_period*i*cos(n_period*i*phi0(j)), dHZ(2*i+1), tol, 'Sine match '//i_s)
    end do
  end do
end subroutine test_sincosperiod_moivre
end module test_moivre_series
