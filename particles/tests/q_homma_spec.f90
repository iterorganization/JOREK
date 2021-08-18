!> This module contains some testcases for calculating the heat flux with classical
!> transport coefficients
module q_homma_spec
use mod_collisions
use constants
use fruit
implicit none

!> Global parameters
real*8, parameter  :: n_b = 1d20 !< [m^-3]
real*8, parameter  :: m_b = 1 !< [u]
integer*1, parameter :: q_b = 1 !< [e]
real*8, parameter :: B(3) = [0.d0, 0.d0, 0.1d0] !< [T]
real*8, parameter :: kT = 50.d0*EL_CHG !< [J]
logical :: fixed_coulomb_log = .true.
contains

!> Calculate analytically the heat flux vector and compare our implementation
!> Test for a dia heatflux, i.e. we'll have diamagnetic drift
subroutine test_q_homma2013_perp
  real*8, parameter :: grad_kT(3) = [100.d0, 0.d0, 0.d0]*EL_CHG ! [J/m]
  real*8, dimension(3) :: q
  real*8, parameter :: tol = 1.d0

  q = q_homma2013(kT, grad_kT, B, n_b, m_b, q_b)
  ! our gradient is in x direction
  if (fixed_coulomb_log) then
    call assert_equals(-33910.4d0, q(1), tol, 'perp component 33 kW/m^2') ! -33910.4d0 calculated with ln_lambda = 15
  else
    call assert_equals(-25820d0, q(1), tol, 'perp component 25 kW/m^2')
    ! -25820d0 calculated with temperature-dependent coulomb log
  end if
  ! Since our gradient is in x direction the dia direction is y
  call assert_equals(2.002720d6, q(2), tol, 'dia component 2 MW/m^2')
  call assert_equals(0.d0, q(3), tol, 'parallel component 0')
end subroutine test_q_homma2013_perp

!> Test for a parallel heatflux
subroutine test_q_homma2013_par
  real*8, parameter :: grad_kT(3) = [0.d0, 0.d0, 100.d0]*EL_CHG ! [J/m]
  real*8, dimension(3) :: q
  real*8, parameter :: tol = 1.d0

  q = q_homma2013(kT, grad_kT, B, n_b, m_b, q_b)
  call assert_equals(0.d0, q(1), tol, 'x component 0')
  call assert_equals(0.d0, q(2), tol, 'y component 0')
  if (fixed_coulomb_log) then
    call assert_equals(-1.47612167d8, q(3), tol, 'parallel component 147 MW/m^2')
  else
    call assert_equals(-1.93861d8, q(3), tol, 'parallel component 147 MW/m^2')
  end if
end subroutine test_q_homma2013_par

end module q_homma_spec
