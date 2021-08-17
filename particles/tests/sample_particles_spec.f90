!> This module contains some testcases for the binary collision model implementation
!> (sampling part only)
module sample_particles_spec
use mod_pcg32_rng
use mod_random_seed
use mod_collisions
use constants
use fruit
implicit none

!> Global parameters
real*8, parameter  :: n_b = 1d20 !< [m^-3]
real*8, parameter  :: m_b = 1 !< [u]
integer*1, parameter :: q_b = 1 !< [e]
real*8, parameter  :: v_b(3) = 0
real*8, parameter  :: coulomb_log = 15

contains

subroutine test_averages_100
  call averages_for_n(100)
end subroutine test_averages_100
subroutine test_averages_1000
  call averages_for_n(1000)
end subroutine test_averages_1000
subroutine test_averages_10000
  call averages_for_n(10000)
end subroutine test_averages_10000
subroutine test_averages_100000
  call averages_for_n(100000)
end subroutine test_averages_100000
subroutine test_averages_1000000
  call averages_for_n(1000000)
end subroutine test_averages_1000000

!> Verify that the moments of the distribution function are correct (no grad T)
subroutine averages_for_n(n)
type(pcg32_rng) :: rng
integer, intent(in) :: n
real*8 :: w_1(3,n), u(6,n)
real*8, parameter :: ktb = 50.d0, grad_ktb(3) = 0.d0
integer :: i

call rng%initialize(6, 1234897543, 1, 1)
do i=1,n
  call rng%next(u(:,i))
end do
call sample_velocity_dist_unmagnetized(n, u, coulomb_log, ktb*EL_CHG, grad_ktb, n_b, m_b, q_b, v_b, w_1)

! Verify moments of distribution function
! Momentum (scaling is correct but absolute values are large)
call assert_equals(0.d0, sum(m_b*ATOMIC_MASS_UNIT*w_1(1,:))/real(n), 18d4/sqrt(real(n)), "x momentum zero")
call assert_equals(0.d0, sum(m_b*ATOMIC_MASS_UNIT*w_1(2,:))/real(n), 18d4/sqrt(real(n)), "y momentum zero")
call assert_equals(0.d0, sum(m_b*ATOMIC_MASS_UNIT*w_1(3,:))/real(n), 18d4/sqrt(real(n)), "z momentum zero")

! Mean speed
call assert_equals(sqrt(8.d0*ktb*el_chg/(PI*m_b*ATOMIC_MASS_UNIT)), sum(norm2(w_1,dim=1))/real(n), &
    8d4/sqrt(real(n)), 'mean speed correct')
! RMS speed
call assert_equals(sqrt(3.d0*ktb*el_chg/(m_b*ATOMIC_MASS_UNIT)), sqrt(sum(norm2(w_1,dim=1)**2)/real(n)), &
    8d4/sqrt(real(n)), 'RMS speed correct')

end subroutine averages_for_n
end module sample_particles_spec
