!> Module containing test routines for collisions with the BCM (Homma, Hatayama JCP 2012,2013)
module mod_collision_cases_homma2012
  use mod_collision_case
  use constants
  implicit none
  private
  public :: homma_2012_collision_case, collision_case
  public :: n_cases, case_numbers, f_theo

  ! Cases from Homma JCP 2012
  ! All particles launched in the z-direction, from the origin
  ! No gradients in y-direction
  integer, parameter :: n_cases = 9 
  real*8, parameter  :: kT(n_cases) = 50.d0*EL_CHG !< [J]
  real*8, parameter  :: grad_KT_x(n_cases) = &!< [J]
  ! 0-0, 1-1, 1-2, 2-1, 2-2, 3-1, 3-2, 4-1, 4-2
  [ 0d0, 0d0, 0d0, 3d0, 5d0, 0d0, 0d0, 0d0, 0d0]*EL_CHG
  real*8, parameter  :: grad_KT_z(n_cases) = &!< [J]
  [ 0d0, 3d0, 5d0, 0d0, 0d0, 0d0, 1d1, 5d0, 5d0]*EL_CHG
  real*8, parameter  :: v_normalized(n_cases) = &!< [dimensionless], Homma 2012 eq 7
  [0.1278,0.1278,0.1278,0.1278,0.1278,1.80,1.80,0.1278,0.1278] ! v_z only
  real*8, parameter  :: n_b(n_cases) = & !< [m^-3]
  [1d20,1d20,1d20,1d20,1d20,1d20,1d20,0.85d20,2d20]
  real*8, parameter  :: E(3) = 0.d0, B(3) = 0.d0 !< No electric or magnetic fields
  real*8, parameter  :: coulomb_log = 15d0
  integer*1, parameter :: q_a = 3 !< [e]
  integer*1, parameter :: q_b = 1 !< [e]
  real*8, parameter  :: m_a = 183.84 !< [u] (tungsten)
  character(len=3), parameter :: case_numbers(n_cases) = ["0-0", "1-1", "1-2", "2-1", "2-2", &
                                                          "3-1", "3-2", "4-1", "4-2"]
  real*8, parameter :: f_theo(n_cases) = [-3.83d0, -2.7d0, -1.95d0, 1.14d0, 1.91d0, -22.5d0, -24.2d0, -1.38d0, -5.78d0] ! Table 3
  ! from 2012 paper. Units of 10^-17 N

  real*8, parameter  :: m_b = 1.0 ! < [u] (hydrogen)
contains

!> Select a single case from the lists above
function homma_2012_collision_case(n) result(c)
  integer :: n
  type(collision_case) :: c
  if (n .gt. n_cases .or. n .lt. 1) then
    write(*,*) "Invalid N given!"
    call exit(1)
  end if
  c%kT = kT(n)
  c%grad_kT = [grad_kT_x(n), 0.d0, grad_kT_z(n)]
  c%v_0 = [0.d0, 0.d0, v_normalized(n)*sqrt(c%kT/(m_b*ATOMIC_MASS_UNIT))] ! note that the reference velocity is the background thermal velocity
  c%n_b = n_b(n)
  c%E = E
  c%B = B
  c%q_a = q_a
  c%q_b = q_b
  c%m_a = m_a
  c%m_b = m_b
  c%x_0 = 0.d0
  c%coulomb_log = coulomb_log

  c%slowdown_time = (1.d0/( &
                     (1.d0+m_a/m_b)* &
                     mu((m_b*ATOMIC_MASS_UNIT)*c%v_0(3)**2/(2.d0*c%kT)))) * &
      (4.d0*PI*(EPS_ZERO**2) * ((m_a*ATOMIC_MASS_UNIT)**2) * (c%v_0(3)**3)) / &
      ((real(q_a)**2) * (real(q_b)**2) * EL_CHG**4 * c%n_b * c%coulomb_log)
end function homma_2012_collision_case

!> Integral of 2/sqrt(pi) int_0^x exp(-t)sqrt(t)dt
pure function mu(x)
  real*8, intent(in) :: x
  real*8 :: mu
  mu = erf(sqrt(x)) - sqrt(2.d0/PI)*exp(-x)*sqrt(2.d0*x)
end function mu
end module mod_collision_cases_homma2012
