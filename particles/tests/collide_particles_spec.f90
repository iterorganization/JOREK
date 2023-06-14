!> This module contains some testcases for the binary collision model implementation
module collide_particles_spec
use mod_pcg32_rng
use mod_random_seed
use mod_collisions
use fruit
implicit none

!> Global parameters
real*8, parameter  :: n_b = 1d20 !< [m^-3]
real*8, parameter  :: coulomb_log = 15
real*8, parameter  :: dt = 7.03d-10 !< [s]

contains

!> Test a like-particle collision against conservation of energy and momentum
subroutine test_conservation_laws
type(pcg32_rng) :: rng
integer, parameter :: n_tries = 100
real*8  :: m_a, m_b, r(13), v_a(3), v_b(3)
real*8  :: v_a0(3), v_b0(3)
integer*1 :: q_a, q_b, i

call rng%initialize(13, random_seed(), 1, 1)

do i=1,n_tries
  call rng%next(r)
  ! If these are the same the test works!
  r(7) = 0.5d0
  r(9) = 1.d0
  ! setup velocities
  v_a = r(1:3)*1d3
  v_b = r(4:6)*1d3
  v_a0 = v_a
  v_b0 = v_b
  m_a = r(7)*100.d0
  q_a = nint(r(8)*10.d0,1)
  m_b = r(9)*100.d0
  q_b = nint(r(10)*10.d0,1)

  call collide_particles(r(11:13), q_a, m_a, v_a, q_b, m_b, v_b, n_b, coulomb_log, dt)

  ! Test energy conservation
  call assert_equals(m_a*sq(v_a0)+m_b*sq(v_b0), m_a*sq(v_a)+m_b*sq(v_b), 1d-6, 'energy must remain same') ! absolute tolerance, so ok
  call assert_equals(0.d0, norm2(m_a*v_a+m_b*v_b - (m_a*v_a0+m_b*v_b0)), 1d-10, 'momentum must remain same')

end do
end subroutine test_conservation_laws

!> Test that a collision with heavy and light particles produces the right result
!> Calculate a single case manually and then perform it with `collide_particles`
!> Use as test a W3+ particle in a H+ background, with opposite velocities
!> Preselect some random numbers to calculate the impact parameters
subroutine test_heavy_light_head_on
use mod_sampling
use constants
real*8, parameter :: m_a = 183.84d0, m_b=1.d0
real*8, parameter :: v_a(3) = [0.d0,0.d0,1.d4], v_b(3) = [0.d0,0.d0,-1.d5]
integer*1, parameter :: q_a = 3, q_b = 1
real*8 :: r(3) = [0.1d0, 0.4d0, 0.2d0] ! 'random' numbers for this test
real*8, parameter :: tol = 2d-11
real*8 :: v_a_test(3), v_b_test(3), delta, theta, tmp(2), mab, v_rel(3), delta_v(3), u, phi
v_a_test = v_a
v_b_test = v_b

tmp = boxmueller_transform(r(1:2))
mab = m_a*m_b/(m_a+m_b)
v_rel = v_a-v_b
u     = norm2(v_rel)
delta = sqrt(q_a**2*q_b**2*EL_CHG**4*n_b*coulomb_log*dt/(8.d0*PI*eps_zero**2*mab**2*ATOMIC_MASS_UNIT**2*u**3))
theta = 2.d0*atan(delta*tmp(1))
phi   = r(3)*TWOPI
delta_v = [u*sin(theta)*cos(phi), u*sin(theta)*sin(phi), -u*(1-cos(theta))]
write(*,*) "headon", norm2(delta_v), delta

call collide_particles(r, q_a, m_a, v_a_test, q_b, m_b, v_b_test, n_b, coulomb_log, dt)
! Test all components
call assert_equals(v_a(1)+mab/m_a*delta_v(1), v_a_test(1), tol, 'x-component of v_a')
call assert_equals(v_a(2)+mab/m_a*delta_v(2), v_a_test(2), tol, 'y-component of v_a')
call assert_equals(v_a(3)+mab/m_a*delta_v(3), v_a_test(3), tol, 'z-component of v_a')
call assert_equals(v_b(1)-mab/m_b*delta_v(1), v_b_test(1), tol, 'x-component of v_b')
call assert_equals(v_b(2)-mab/m_b*delta_v(2), v_b_test(2), tol, 'y-component of v_b')
call assert_equals(v_b(3)-mab/m_b*delta_v(3), v_b_test(3), tol, 'z-component of v_b')
write(*,*) v_b(3)-mab/m_b*delta_v(3)- v_b_test(3)
end subroutine test_heavy_light_head_on

!> Test that a collision with heavy and light particles produces the right result
!> Calculate a single case manually and then perform it with `collide_particles`
!> Use as test a W3+ particle in a H+ background, with velocities under an angle
!> Preselect some random numbers to calculate the impact parameters
subroutine test_heavy_light_angle
use mod_sampling
use constants
real*8, parameter :: m_a = 183.84d0, m_b=1.d0
real*8, parameter :: v_a(3) = [0.d0,0.d0,1.d4], v_b(3) = [0.d0,1.d5,0.d0]
integer*1, parameter :: q_a = 3, q_b = 1
real*8 :: r(3) = [0.7d0, 0.2d0, 0.5d0] ! 'random' numbers for this test
real*8, parameter :: tol = 1d-12 !< WARNING: very high tolerance!
real*8 :: v_a_test(3), v_b_test(3), delta, theta, tmp(2), mab, v_rel(3), delta_v(3), u, phi
v_a_test = v_a
v_b_test = v_b

tmp = boxmueller_transform(r(1:2))
mab = m_a*m_b/(m_a+m_b)
v_rel = v_a-v_b
u     = norm2(v_rel)
delta = sqrt(q_a**2*q_b**2*EL_CHG**4*n_b*coulomb_log*dt/(8.d0*PI*eps_zero**2*mab**2*ATOMIC_MASS_UNIT**2*u**3))
theta = 2.d0*atan(delta*tmp(1))
phi   = r(3)*TWOPI
delta_v = [v_rel(1)/norm2(v_rel(1:2)) * v_rel(3)*sin(theta)*cos(phi) &
             - v_rel(2)/norm2(v_rel(1:2))*u*sin(theta)*sin(phi) - v_rel(1)*(1.d0-cos(theta)), &
           v_rel(2)/norm2(v_rel(1:2)) * v_rel(3)*sin(theta)*cos(phi) &
             + v_rel(1)/norm2(v_rel(1:2))*u*sin(theta)*sin(phi) - v_rel(2)*(1.d0-cos(theta)), &
           -norm2(v_rel(1:2))*sin(theta)*cos(phi) - v_rel(3)*(1.d0-cos(theta))]
write(*,*) "angle", norm2(delta_v), delta

call collide_particles(r, q_a, m_a, v_a_test, q_b, m_b, v_b_test, n_b, coulomb_log, dt)
! Test all components
call assert_equals(v_a(1)+mab/m_a*delta_v(1), v_a_test(1), tol, 'x-component of v_a')
call assert_equals(v_a(2)+mab/m_a*delta_v(2), v_a_test(2), tol, 'y-component of v_a')
call assert_equals(v_a(3)+mab/m_a*delta_v(3), v_a_test(3), tol, 'z-component of v_a')
call assert_equals(v_b(1)-mab/m_b*delta_v(1), v_b_test(1), tol, 'x-component of v_b')
call assert_equals(v_b(2)-mab/m_b*delta_v(2), v_b_test(2), tol, 'y-component of v_b')
call assert_equals(v_b(3)-mab/m_b*delta_v(3), v_b_test(3), tol, 'z-component of v_b')


end subroutine test_heavy_light_angle

pure function sq(v) result(n)
  real*8, intent(in), dimension(:) :: v
  real*8 :: n
  n = dot_product(v,v)
end function sq
end module collide_particles_spec
