!> Particle pusher module with simple forward euler stepping
!> This module contains routines for pushing particles in the RZPhi (cylindrical)
!> or Cartesian XYZ coordinate systems. The Cartesian version is included mostly
!> for performance comparisons and testing and should not be used in production
!> (as all diagnostics assume RZPhi coordinates)
module mod_fieldline_euler
  use mod_particle_types
  use constants, only: EL_CHG, ATOMIC_MASS_UNIT
  implicit none
  private

  public fieldline_euler_push_cylindrical, fieldline_euler_push_cartesian
  public fieldline_adams_bashforth_push_cylindrical, fieldline_adams_bashforth_push_cartesian
  public gc_to_fieldline, v_fieldline_adams_bashforth_push_cylindrical
contains

!> Follow a fieldline for a single timestep with forward euler
!> This routine works in RZPhi coordinates
pure subroutine fieldline_euler_push_cylindrical(particle, B, dt)
  type(particle_fieldline), intent(inout) :: particle
  real*8, dimension(3), intent(in) :: B
  real*8, intent(in) :: dt
  real*8 :: R, Rphi
  real*8 :: B_hat(3)
  B_hat = B / norm2(B)
  R    = particle%x(1) + B_hat(1)*particle%v * dt
  RPhi = particle%v*B_hat(3) * dt

  ! Calculate the new R, Phi, Z
  particle%x(1) = sqrt(R**2 + RPhi**2)
  particle%x(2) = particle%x(2) + dt * particle%v*B_hat(2)
  particle%x(3) = particle%x(3) + asin(RPhi / particle%x(1))
end subroutine fieldline_euler_push_cylindrical

!> Follow a fieldline for a single timestep with a Two-step Adams-Bashfort method
!> This routine works in RZPhi coordinates.
!> B_hat_prev must be set in the particle or the first step will be inaccurate
pure subroutine fieldline_adams_bashforth_push_cylindrical(particle, B, dt)
  type(particle_fieldline), intent(inout) :: particle
  real*8, dimension(3), intent(in) :: B
  real*8, intent(in) :: dt
  real*8 :: R, Rphi
  real*8 :: B_hat(3)
  B_hat = B / norm2(B)

  ! No cylindrical correction! works better because adams-bashforth needs linear steps
  particle%x(3) = (particle%x(3)*particle%x(1)  + (B_hat(3)*1.5d0 - particle%B_hat_prev(3)*0.5d0) * particle%v * dt)/particle%x(1)
  particle%x(1) = particle%x(1)                 + (B_hat(1)*1.5d0 - particle%B_hat_prev(1)*0.5d0) * particle%v * dt
  particle%x(2) = particle%x(2) + dt * particle%v*(B_hat(2)*1.5d0 - particle%B_hat_prev(2)*0.5d0)
  particle%B_hat_prev = B_hat
end subroutine fieldline_adams_bashforth_push_cylindrical

!> Follow a VELOCITY fieldline for a single timestep with a Two-step Adams-Bashfort method
!> This routine works in RZPhi coordinates.
!> v_hat_prev must be set in the particle or the first step will be inaccurate
pure subroutine v_fieldline_adams_bashforth_push_cylindrical(particle, V, dt, V_add)
  type(particle_marker), intent(inout) :: particle
  real*8, dimension(3), intent(in) :: V
  real*8, intent(in) :: dt
  real*8, dimension(3), intent(in), optional :: V_add ! Additional velocity given externally apart from the fluid velocity
  real*8 :: R, Rphi, V_value, V_value_prev
  real*8 :: V_hat(3), V_hat_prev(3), V_ext(3)

  V_ext        = 0.0
  if (present(V_add)) V_ext = V_add

  V_value      = norm2(V+V_ext)
  V_value_prev = norm2(particle%V_prev)
  V_hat        = (V+V_ext) / V_value
  V_hat_prev   = particle%V_prev / V_value_prev

  ! No cylindrical correction! works better because adams-bashforth needs linear steps
  particle%x(3) = (particle%x(3)*particle%x(1) + (V_hat(3)*1.5d0 - V_hat_prev(3)*0.5d0) * V_value * dt)/particle%x(1)
  particle%x(1) = particle%x(1)                + (V_hat(1)*1.5d0 - V_hat_prev(1)*0.5d0) * V_value * dt
  particle%x(2) = particle%x(2)                + (V_hat(2)*1.5d0 - V_hat_prev(2)*0.5d0) * V_value * dt

  particle%V_prev = V + V_ext
end subroutine v_fieldline_adams_bashforth_push_cylindrical

!> Follow a fieldline for a single timestep with forward euler
!> This routine works in RZPhi coordinates
pure subroutine fieldline_euler_push_cartesian(particle, B, dt)
  type(particle_fieldline), intent(inout) :: particle
  real*8, dimension(3), intent(in) :: B
  real*8, intent(in) :: dt
  real*8 :: B_hat(3)
  B_hat = B / norm2(B)
  particle%x = particle%x + particle%v*B_hat*dt
end subroutine fieldline_euler_push_cartesian

!> Follow a fieldline for a single timestep with a Two-step Adams-Bashfort method
!> This routine works in RZPhi coordinates.
!> B_hat_prev must be set in the particle or the first step will be inaccurate
pure subroutine fieldline_adams_bashforth_push_cartesian(particle, B, dt)
  type(particle_fieldline), intent(inout) :: particle
  real*8, dimension(3), intent(in) :: B
  real*8, intent(in) :: dt
  real*8 :: B_hat(3)
  B_hat = B / norm2(B)
  particle%x = particle%x + particle%v*dt*(B_hat*1.5d0-particle%B_hat_prev*0.5d0)
  particle%B_hat_prev = B_hat
end subroutine fieldline_adams_bashforth_push_cartesian

!> Take a particle_gc and get the fieldline particle.
function gc_to_fieldline(in) result(out)
  use constants
  use data_structure
  type(particle_gc), intent(in)       :: in
  type(particle_fieldline)     :: out !< Particle of which to update x and v
  out = in
end function gc_to_fieldline
end module mod_fieldline_euler
