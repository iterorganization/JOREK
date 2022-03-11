!> Particle pusher module with the Boris scheme.
!> This module contains routines for pushing particles in the RZPhi (cylindrical)
!> or Cartesian XYZ coordinate systems. The Cartesian version is included mostly
!> for performance comparisons and testing and should not be used in production
!> (as all diagnostics assume RZPhi coordinates)
module mod_boris
  use mod_particle_types
  use constants, only: EL_CHG, ATOMIC_MASS_UNIT
  implicit none
  private

  public boris_push_cylindrical, boris_push_cartesian, boris_initial_half_step_backwards_XYZ
  public boris_initial_half_step_backwards_RZPhi
  public boris_all_initial_half_step_backwards_RZPhi
  public gc_to_kinetic, kinetic_to_gc
  public kinetic_to_kinetic_leapfrog, kinetic_leapfrog_to_kinetic
  public kinetic_leapfrog_to_gc, gc_to_kinetic_leapfrog
contains

!> Push a single particle for some timesteps with the boris method
!> See G.L. Delzanno, E. Camporeale / JCP 253 (2013) 259-277 for details.
!> This routine works in RZPhi coordinates
pure subroutine boris_push_cylindrical(particle, m, E, B, dt)
  use mod_math_operators, only: cross_product
  type(particle_kinetic_leapfrog), intent(inout)  :: particle
  real*8, intent(in) :: m
  real*8, dimension(3), intent(in) :: E, B
  real*8, intent(in) :: dt
  real*8 :: R, Rphi
  real*8 :: fE, fB, eom
  real*8 :: B2, Bnorm
  eom = EL_CHG / (m * ATOMIC_MASS_UNIT)

  B2    = dot_product(B,B)
  Bnorm = sqrt(B2)

  ! update the velocity from v^(n-1/2) to v^(n+1/2)
  ! Calculate the geometric factor f = tan(q/m delta_t/2 |B|)/|B|
  fE =     particle%q*eom * dt * 0.5d0
  fB = tan(particle%q*eom * dt * 0.5d0 * Bnorm) / Bnorm

  ! Calculate the electric field update (v^n-1/2 -> v-) with the Boris method
  particle%v = particle%v + fE * E
  ! Calculate the rotation
  particle%v = (particle%v + 2.d0*fB/(1.d0+fB*fB*B2)*( &
    cross_product(particle%v,B) &
    - fB * particle%v * B2 &
    + fB * B * dot_product(particle%v,B)))
  ! Calculate the next electric field update (v+ -> v^n+1/2)
  particle%v = particle%v + fE * E
  ! update the position from v^n to v^(n+1)
  ! Calculate the new R and RPhi
  R    = particle%x(1) + particle%v(1) * dt
  RPhi = particle%v(3) * dt

  ! Calculate the new R, Phi, Z
  particle%x(1) = sqrt(R**2 + RPhi**2)
  particle%x(2) = particle%x(2) + dt * particle%v(2)
  particle%x(3) = particle%x(3) + asin(RPhi / particle%x(1))

  ! Adjust R and Phi velocities (component 1 and 3) to the new reference frame
  particle%v(1:3:2) = [R     * particle%v(1) + RPhi * particle%v(3), &
                       -RPhi * particle%v(1) + R    * particle%v(3)] / particle%x(1)
end subroutine boris_push_cylindrical

!> Push a single particle for some timesteps with the boris method
!> See G.L. Delzanno, E. Camporeale / JCP 253 (2013) 259-277 for details
!> This routine works in a left-handed cartesian coordinate system (XYZ)
pure subroutine boris_push_cartesian(particle, m, E, B, dt)
  use mod_math_operators, only: cross_product
  class(particle_kinetic_leapfrog), intent(inout)  :: particle
  real*8, intent(in) :: m
  real*8, dimension(3), intent(in) :: E, B
  real*8, intent(in) :: dt
  real*8  :: fE, fB, eom
  real*8  :: B2, Bnorm
  eom = EL_CHG / (m * ATOMIC_MASS_UNIT)

  B2    = dot_product(B,B)
  Bnorm = sqrt(B2)

  ! update the velocity from v^(n-1/2) to v^(n+1/2)
  ! Calculate the geometric factor f = tan(q/m delta_t/2 |B|)/|B|
  fE =     particle%q*eom * dt * 0.5d0
  fB = tan(particle%q*eom * dt * 0.5d0 * Bnorm) / Bnorm

  ! Calculate the electric field update (v^n-1/2 -> v-) with the Boris method
  particle%v = particle%v + fE * E
  ! Calculate the rotation
  particle%v = (particle%v + 2.d0*fB/(1.d0+fB*fB*B2)*( &
    cross_product(particle%v,B) &
    - fB * particle%v * B2 &
    + fB * B * dot_product(particle%v,B)))
  ! Calculate the next electric field update (v+ -> v^n+1/2)
  particle%v = particle%v + fE * E
  particle%x = particle%x + particle%v * dt
end subroutine boris_push_cartesian



!> Given a particle with position x and velocity v at time t=0 (t^0), calculate v^-1/2
!> (see G.L. Delzanno, E. Camporeale / JCP 253 (2013) 259-277 for details)
pure subroutine boris_initial_half_step_backwards_XYZ(particle, m, E, B, dt)
  use constants, only: EL_CHG, ATOMIC_MASS_UNIT
  use mod_math_operators, only: cross_product
  class(particle_kinetic_leapfrog), intent(inout) :: particle
  real*8, intent(in) :: m
  real*8, dimension(3), intent(in) :: E, B
  real*8, intent(in) :: dt
  real*8, dimension(3) :: v !< for calculating the initial half-step
  real*8 :: f, B2
  f = - (EL_CHG * real(particle%q)) / (ATOMIC_MASS_UNIT * m) * dt * 0.25d0
  B2 = dot_product(B, B)
  v = particle%v + f*E
  v = (v + 2.d0*f/(1.d0+f**2*B2) &
      * (cross_product(v, B) - f*v*B2 + f*B*dot_product(v,B)))
  v = v + f*E
  particle%v = v
end subroutine boris_initial_half_step_backwards_XYZ

!> Given a particle with position x and velocity v at time t=0 (t^0), calculate v^-1/2
!> in cylindrical coordinates this is the same as in cartesian coordinates
!> (see G.L. Delzanno, E. Camporeale / JCP 253 (2013) 259-277 for details)
pure subroutine boris_initial_half_step_backwards_RZPhi(particle, m, E, B, dt)
  use constants, only: EL_CHG, ATOMIC_MASS_UNIT
  use mod_math_operators, only: cross_product
  class(particle_kinetic_leapfrog), intent(inout) :: particle
  real*8, intent(in) :: m
  real*8, dimension(3), intent(in) :: E, B
  real*8, intent(in) :: dt
  real*8, dimension(3) :: v !< for calculating the initial half-step
  real*8 :: f, B2
  f = - (EL_CHG * real(particle%q)) / (ATOMIC_MASS_UNIT * m) * dt * 0.25d0
  B2 = dot_product(B, B)
  v = particle%v + f*E
  v = (v + 2.d0*f/(1.d0+f**2*B2) &
      * (cross_product(v, B) - f*v*B2 + f*B*dot_product(v,B)))
  v = v + f*E
  particle%v = v
end subroutine boris_initial_half_step_backwards_RZPhi

!> Given a list of particles and the fields, perfom half-steps backwards for all
subroutine boris_all_initial_half_step_backwards_RZPhi(particles, m, fields, t, dt)
  use mod_fields
  class(particle_kinetic_leapfrog), intent(inout), dimension(:) :: particles
  real*8, intent(in) :: m
  class(fields_base), intent(in) :: fields
  real*8, intent(in) :: t !< Time of the current simulation
  real*8, intent(in) :: dt !< Timestep

  integer :: i
  real*8  :: psi, U, E(3), B(3)

#ifndef __NVCOMPILER
  !$omp parallel do default(none)   &
  !$omp private(E, B, psi, U, i)    &
  !$omp shared(particles, fields, dt, m, t)
#endif
  do i=1,size(particles,1)
    if (particles(i)%i_elm .le. 0) cycle
    call fields%calc_EBpsiU(t, particles(i)%i_elm, particles(i)%st, particles(i)%x(3), E, B, psi, U)
    call boris_initial_half_step_backwards_RZPhi(particles(i), m, E, B, dt)
  end do
#ifndef __NVCOMPILER
  !$omp end parallel do
#endif
end subroutine boris_all_initial_half_step_backwards_RZPhi

!> Shortcut functions for converting between kinetic leapfrog and GC,
!> by converting first to kinetic as an intermediate step.
function kinetic_leapfrog_to_gc(node_list, element_list, in, E, B, mass, dt) result(out)
  use data_structure
  type(type_node_list), intent(in)            :: node_list
  type(type_element_list), intent(in)         :: element_list
  type(particle_kinetic_leapfrog), intent(in) :: in
  real*8, dimension(3), intent(in)            :: E !< Electric field at kinetic position [V/m]
  real*8, dimension(3), intent(in)            :: B !< Magnetic field at kinetic position [T]
  real*8, intent(in)                          :: mass !< Mass of the particle [amu]
  real*8, intent(in)                          :: dt !< Timestep [s]
  type(particle_gc)             :: out
  out = kinetic_to_gc(node_list, element_list, kinetic_leapfrog_to_kinetic(in, E, B, mass, dt),B, mass)
end function kinetic_leapfrog_to_gc

!> Shortcut functions for converting between GC and kinetic leapfrog
!> by converting first to kinetic as an intermediate step.
function gc_to_kinetic_leapfrog(in, node_list, element_list, chi, E, B, mass, dt) result(out)
  use data_structure
  type(particle_gc), intent(in)               :: in
  real*8, dimension(3), intent(in)            :: E !< Electric field at kinetic position [V/m]
  real*8, dimension(3), intent(in)            :: B !< Magnetic field at kinetic position [T]
  real*8, intent(in)                          :: chi !< Gyrophase to select from the guiding center ring
  type(type_node_list), intent(in)            :: node_list
  type(type_element_list), intent(in)         :: element_list
  real*8, intent(in)                          :: mass !< Mass of the particle [amu]
  real*8, intent(in)                          :: dt !< Timestep [s]
  type(particle_kinetic_leapfrog)             :: out
  out = kinetic_to_kinetic_leapfrog(gc_to_kinetic(node_list, element_list, in, chi, B, mass), E, B, mass, dt)
end function gc_to_kinetic_leapfrog


!> Convert a particle_kinetic_leapfrog to particle_kinetic by performing a step forwards for the velocity only and averaging
function kinetic_to_kinetic_leapfrog(in, E, B, mass, dt) result(out)
  use data_structure
  type(particle_kinetic), intent(in)          :: in
  real*8, dimension(3), intent(in)            :: B !< Magnetic field at kinetic position [T]
  real*8, dimension(3), intent(in)            :: E !< Electric field at kinetic position [V/m]
  real*8, intent(in)                          :: mass !< Mass of the particle [amu]
  real*8, intent(in)                          :: dt !< Timestep,[s]
  type(particle_kinetic_leapfrog)             :: out
  out = in
  out%q = in%q
  out%v = in%v
  call boris_initial_half_step_backwards_RZPhi(out, mass, E, B, dt)
end function kinetic_to_kinetic_leapfrog

!> Convert a particle_kinetic to particle_kinetic_leapfrog by performing a half step backwards
function kinetic_leapfrog_to_kinetic(in, E, B, mass, dt) result(out)
  use data_structure
  type(particle_kinetic_leapfrog), intent(in) :: in
  real*8, dimension(3), intent(in)            :: B !< Magnetic field at kinetic position [T]
  real*8, dimension(3), intent(in)            :: E !< Electric field at kinetic position [V/m]
  real*8, intent(in)                          :: mass !< Mass of the particle [amu]
  real*8, intent(in)                          :: dt !< Timestep [s]
  type(particle_kinetic)                      :: out
  type(particle_kinetic_leapfrog) :: tmp ! needed because we cannot alter `in`
  tmp = in
  call boris_initial_half_step_backwards_RZPhi(tmp, mass, E, B, -dt)
  out = in
  out%q = in%q
  out%v = tmp%v
end function kinetic_leapfrog_to_kinetic


!> Take a particle_kinetic and B and return a particle_guiding_centre.
!> This is an approximation based on the magnetic field at the kinetic position.
!> gc_to_kinetic is not a true inverse of this function! Use with care. (only works if B uniform)
!> This should actually use the current timestep size and fields to be a bit more accurate.
!>
!> See for instance http://iter.rma.ac.be/Stufftodownload/Texts/GuidingCenterMotion.pdf
!> 
!> The kinetic position is given by
!> \[ x = x_{gc} - \frac{m}{q B^2} v \times B \]
!> Which can easily be adjusted to obtain the gc position from a kinetic position and velocity.
function kinetic_to_gc(node_list, element_list, in, B, mass) result(out)
  use data_structure
  use mod_math_operators, only: cross_product
  use mod_find_rz_nearby
  type(type_node_list), intent(in)            :: node_list
  type(type_element_list), intent(in)         :: element_list
  type(particle_kinetic), intent(in)          :: in
  real*8, dimension(3), intent(in)            :: B !< Magnetic field at kinetic position [T]
  real*8, intent(in)                          :: mass !< Mass of the particle [amu]
  type(particle_gc)                           :: out
  real*8  :: B_hat(3), B_norm, v_par, v2
  integer :: ifail

!  call copy_particle_base(in, out)
  out = in
  out%q      = in%q

  B_norm = norm2(B)
  B_hat = B/B_norm
  v_par = dot_product(in%v,B_hat)
  v2    = dot_product(in%v,in%v)
  ! Calculate GC position
  if (out%q .ne. 0) then
    out%x = in%x + (mass*ATOMIC_MASS_UNIT*cross_product(in%v,B_hat))/(in%q*EL_CHG*B_norm)
  else
    out%x = in%x
  end if

  ! Calculate velocity-related variables
  out%E  = mass * ATOMIC_MASS_UNIT * 0.5d0 * v2 / EL_CHG ! [eV]
  ! Perhaps we could also use the magnetic field in the guiding center to calculate
  ! v_perp, but don't do it for now.
  out%mu = sign(mass * ATOMIC_MASS_UNIT * 0.5d0 * (v2 - v_par**2)&
      /B_norm/EL_CHG, v_par) ! [eV/T]

  ! Calculate new st and i_elm
  call find_RZ_nearby(node_list, element_list, in%x(1), in%x(2), in%st(1), in%st(2), in%i_elm, &
      out%x(1), out%x(2), out%st(1), out%st(2), out%i_elm, ifail)
end function kinetic_to_gc

!> Take a particle_gc and get the kinetic particle.
!> Only updates performed are a transform of E and mu (and chi) to v, and of x_gc to x
!> Calls find_RZ_nearby to find st
!> See for reference http://iter.rma.ac.be/Stufftodownload/Texts/GuidingCenterMotion.pdf
function gc_to_kinetic(node_list, element_list, in, chi, B, mass) result(out)
  use constants
  use data_structure
  use mod_pusher_tools, only: get_orthonormals
  use mod_math_operators, only: cross_product
  use mod_find_rz_nearby
  type(type_node_list), intent(in)    :: node_list
  type(type_element_list), intent(in) :: element_list
  type(particle_gc), intent(in)       :: in
  type(particle_kinetic)     :: out !< Particle of which to update x and v
  real*8, intent(in) :: chi !< Gyrophase [0,2pi]
  real*8, intent(in) :: B(3) !< Magnetic field at GC position [T]
  real*8, intent(in) :: mass !< Mass of the particle [amu]
  real*8 :: B_norm, v_perp, v_par, B_hat(3), e1(3), e2(3)
  integer :: ifail

  out = in
  out%q      = in%q

  B_norm = norm2(B)
  B_hat  = B/B_norm
  ! mu [eV/T] * B_norm [T] * EL_CHG [C] = E [J]
  v_perp = sqrt(2.d0*abs(in%mu*B_norm*EL_CHG)/(mass*ATOMIC_MASS_UNIT)) ! [m/s]
  v_par  = sign(sqrt(2.d0*(in%E-abs(in%mu)*B_norm)*EL_CHG/(mass*ATOMIC_MASS_UNIT)),in%mu)

  ! Define chi as the angle of the velocity vector with b x r
  call get_orthonormals(B_hat, e1, e2)
  out%v  = v_par * B_hat + v_perp * (cos(chi) * e1 + sin(chi) * e2)

  if (out%q .ne. 0) then
    out%x = in%x - (mass*ATOMIC_MASS_UNIT*cross_product(out%v,B_hat))/(real(out%q,8)*EL_CHG*B_norm)
  else
    out%x = in%x
  end if

  ! Calculate new st and i_elm
  call find_RZ_nearby(node_list, element_list, in%x(1), in%x(2), in%st(1), in%st(2), in%i_elm, &
      out%x(1), out%x(2), out%st(1), out%st(2), out%i_elm, ifail)
end function gc_to_kinetic
end module mod_boris
