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

  public fieldline_euler_push_cylindrical, fieldline_euler_push_cartesian, field_line_runge_kutta_fixed_dt_push_jorek
  public fieldline_adams_bashforth_push_cylindrical, fieldline_adams_bashforth_push_cartesian
  public gc_to_fieldline
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

!> Calculate RHS of the field line tracer
!> This just evaluates the magnetic field direction in cylindrical coordinates along which
!> the marker is pushed. Intented to be used with the RK4 solver.
subroutine compute_field_line_rhs(fields,n_variables, &
     n_int_parameters,n_real_parameters,t,solution_old,solution,            &
     int_parameters,real_parameters,derivatives,i_elm_new)
  
  !> load modules
  use mod_fields, only: fields_base
  use mod_find_rz_nearby
  implicit none
  
  !> declare input variables
  class(fields_base), intent(in)                         :: fields
  integer, intent(in)                                    :: n_variables, n_int_parameters, n_real_parameters
  integer, dimension(n_variables), intent(in)            :: int_parameters
  real(kind=8), intent(in)                               :: t
  real(kind=8), dimension(n_variables), intent(in)       :: solution, solution_old
  real(kind=8), dimension(n_real_parameters), intent(in) :: real_parameters
  
  !> declare output variables
  integer, intent(out)                              :: i_elm_new
  real(kind=8), dimension(n_variables), intent(out) :: derivatives

  real(kind=8) :: E(3), B(3), psi, U, st_new(2)
  integer :: ifail

  call find_RZ_nearby(fields%node_list,fields%element_list,solution_old(1), &
       solution_old(2),real_parameters(1),real_parameters(2),                  &
       int_parameters(1),solution(1),solution(2),st_new(1),                    &
       st_new(2),i_elm_new,ifail)

  !> compute required fields
  if(i_elm_new .ne. 0) then
     call fields%calc_EBpsiU(t, i_elm_new, st_new, solution(3), E, B, psi, U)
  end if

  B = B / norm2(B)
  derivatives(1) = B(1)
  derivatives(2) = B(2)
  derivatives(3) = B(3) / solution(1)

end subroutine compute_field_line_rhs

!> Push field line tracer with RK4
subroutine field_line_runge_kutta_fixed_dt_push_jorek(fields, particle, t, dt)
  !> modules
  use mod_fields, only: fields_base
  use mod_find_rz_nearby
  use mod_runge_kutta, only: runge_kutta_fixed_dt
  implicit none

  !> input/output variables
  class(fields_base),       intent(in)    :: fields
  type(particle_fieldline), intent(inout) :: particle
  real(kind=8),             intent(in)    :: t
  real(kind=8),             intent(in)    :: dt

  !> internal variables
  integer                    :: ifail, i_elm_new 
  real(kind=8), dimension(2) :: st_new
  real(kind=8), dimension(3) :: solution_new !< Coordinates 1:R, 2:Z, 3:phi

  !> compute Runge-Kutta differentials
  call runge_kutta_fixed_dt(compute_field_line_rhs, &
       fields,3,1,2,t,dt,[particle%x(1),particle%x(2),particle%x(3)],  &
       [particle%i_elm],[particle%st(1),particle%st(2)], solution_new, i_elm_new)

  !> compute the new local coordinates
  if(i_elm_new.ne.0) call find_rz_nearby(fields%node_list, &
       fields%element_list,particle%x(1),particle%x(2),       &
       particle%st(1),particle%st(2),particle%i_elm,          &
       solution_new(1),solution_new(2),st_new(1),st_new(2),   &
       i_elm_new,ifail)

  !> overwrite GC fields
  particle%x     = solution_new
  particle%st    = st_new
  particle%i_elm = i_elm_new
  
end subroutine field_line_runge_kutta_fixed_dt_push_jorek


end module mod_fieldline_euler
