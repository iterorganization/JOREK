!> Module containing test routines for collisions with the BCM (Homma, Hatayama JCP 2012,2013)
module mod_collision_case
  use constants, only: EL_CHG
  implicit none
  type :: collision_case
    real*8 :: E(3) = [0.d0,0.d0,0.d0], B(3) = [0.d0,0.d0,0.d0] ! [V/m], [T]
    real*8 :: kT = 50.d0*EL_CHG
    real*8 :: grad_kT(3) = 0.d0
    real*8 :: v_0(3), x_0(3) = [0.d0,0.d0,0.d0] ! [m/s], [m]
    real*8 :: n_b = 1d20 ! [m^-3]
    real*8 :: coulomb_log = 15.d0
    integer*1 :: q_a = 3 ! [e] (tungsten)
    integer*1 :: q_b = 1 ! [e] (deuterium)
    real*8 :: m_a = 183.84 ! [u] (tungsten)
    real*8 :: m_b = 2.0 ! [u] (deuterium)

    real*8 :: slowdown_time ! Reference timescale
  end type collision_case
end module mod_collision_case
