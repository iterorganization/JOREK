!> Particle pusher module for integrating full orbits of relativistic
!> particles. Available integrators:
!> -5-steps Volume Preserving Algorithm (VPA): R. Zhang et al., Phys. Plasmas 22 (2015) 044501 
!>                                       (see also C. Sommariva et al., Nucl. Fusion 58 (2018) 016043)
module mod_kinetic_relativistic
use mod_particle_types
use mod_radreactforce
use constants, only: EL_CHG,ATOMIC_MASS_UNIT,SPEED_OF_LIGHT

implicit none

private
public volume_preserving_push_cartesian,volume_preserving_push_jorek
public relativistic_kinetic_to_particle
public gc_to_relativistic_kinetic
public relativistic_kinetic_to_relativistic_gc
public runge_kutta_fixed_dt_relativistic_particle_push
public runge_kutta_fixed_dt_relativistic_particle_push_jorek
public volume_preserving_push_analytical

contains

!---------------------------------------------------------------------------
!> This subroutine computes the first half-step of the JOREK-specific VPA
!> inputs:
!>   half_position: real(8)(3) initial particle position in cartesian coordinates  
!>   particle: (particle_kinetic_relativistic) relativistic particle
!>   mass:     (real8) particle mass in [AMU]
!>   dt:       (real8) time step in [s]
!> outputs:
!>   scaling_factor: (real8)(3) scaling factor to be used in subsequent steps
!>   half_position:  (real8)(2) particle position after half-step 
!>			        in cartesian coordinates
pure subroutine volume_preserving_first_half_step_jorek(particle,half_position,&
       mass,dt,scaling_factor)
  ! input variables
  real(kind=8), intent(in) :: mass, dt !< mass and time step
  ! input/output variables
  real(kind=8),dimension(3),intent(inout) :: half_position
  class(particle_kinetic_relativistic), intent(inout) :: particle !< relativistic particle
  ! output variables
  real(kind=8), intent(out) :: scaling_factor

  scaling_factor = 5.d-1*dt*particle%q*EL_CHG/(ATOMIC_MASS_UNIT*mass*SPEED_OF_LIGHT)
  ! compute dimensionless momentum
  particle%p = particle%p/(mass*SPEED_OF_LIGHT)
  ! compute coordinates at half-step
  half_position = half_position + (5.d-1*dt*SPEED_OF_LIGHT*particle%p) &
                  /(sqrt(1.d0+dot_product(particle%p,particle%p)))
end subroutine volume_preserving_first_half_step_jorek

!---------------------------------------------------------------------------
!> This subroutine computes the second half-step of the JOREK-specific VPA
!> inputs:
!>   particle:       (particle_kinetic_relativistic) relativistic particle
!>   half_position:  (real8)(2) half-position in cartesian coordinates
!>   scaling_factor: (real8) scaling factor for computing momentum
!>   B:		     (real8)(3) magnetic field in cartesian coordinates
!>   E:		     (real8)(3) electric field in cartesian coordinates
!>   mass:	     (real8) particle mass in [AMU]
!>   dt:	     (real8) time step in [s]
!> outputs:
!>   particle: (particle_kinetic_relativistic) relativistic particle type
pure subroutine volume_preserving_second_half_step_jorek(particle,&
                half_position,scaling_factor,E,B,mass,dt)
  ! load methods
  use mod_pusher_tools, only: cayley_transform !< use full Cayley transform
  ! define input/output variables
  class(particle_kinetic_relativistic), intent(inout) :: particle !< relativistic particle
  real(kind=8),dimension(3),intent(inout) :: half_position
  ! define input variables
  real(kind=8),dimension(3),intent(in) :: B, E
  real(kind=8),intent(in) :: scaling_factor, mass, dt

  ! compute momentum at t_(i+1/2)
  particle%p = particle%p + scaling_factor*E
  ! rotate momentum with respect to the magnetic field
  particle%p = matmul(cayley_transform(SPEED_OF_LIGHT*scaling_factor/&
    (sqrt(1.d0+dot_product(particle%p,particle%p))),B),particle%p)
  ! compute momentum at t_(i+1)
  particle%p = particle%p + scaling_factor*E
  ! update position at t_(i+1)
  half_position = half_position + (5.d-1*dt*SPEED_OF_LIGHT*particle%p)/&
    (sqrt(1.d0+dot_product(particle%p,particle%p)))
  ! compute dimensional momentum
  particle%p = particle%p*mass*SPEED_OF_LIGHT
end subroutine volume_preserving_second_half_step_jorek

!---------------------------------------------------------------------------
!> This subroutine integrates a relativistic particle trajectory in JOREK
!> fields using the Volume Preserving Algorithm (VPA)
subroutine volume_preserving_push_jorek(particle,fields,mass,time,timestep,ifail)
  ! load functions
  use mod_coordinate_transforms, only: cartesian_to_cylindrical
  use mod_coordinate_transforms, only: cylindrical_to_cartesian
  use mod_coordinate_transforms, only: vector_cylindrical_to_cartesian
  use mod_fields
  use mod_find_rz_nearby
  ! declare input/output variables
  integer(kind=4),intent(inout) :: ifail
  class(particle_kinetic_relativistic), intent(inout) :: particle !< relativistic particle
  ! declare input variables
  real(kind=8),intent(in) :: mass, time, timestep
  class(fields_base), intent(in) :: fields
  ! declare internal variables
  real(kind=8) :: psi, U
  real(kind=8),dimension(3) :: B, E
  ! half_position coordinates: 1:x, 2:y, 3:z, 4:R, 5:Z, 6:phi 
  real(kind=8),dimension(6) :: half_position
  real(kind=8) :: scaling_factor !< in [s^2*C/(kg*m)]

  ! check if the particle is valid
  if(particle%i_elm.eq.0) return
  ! transform the particle position from cylindrical to cartesian coordinates
  half_position(1:3) = cylindrical_to_cartesian(particle%x)
  ! compute first half-step
  call volume_preserving_first_half_step_jorek(particle,half_position(1:3),&
       mass,timestep,scaling_factor)
  ! calculate cylindrical coordinates from cartesian ones
  half_position(4:6) = cartesian_to_cylindrical(half_position(1:3))
  ! find the (i_elm,s,t) coordinates
  call find_RZ_nearby(fields%node_list,fields%element_list,particle%x(1),&
       particle%x(2),particle%st(1),particle%st(2),particle%i_elm,&
       half_position(4),half_position(5),particle%st(1),particle%st(2),&
       particle%i_elm,ifail)
  ! check if the particle is lost, exit if it is the case
  if(particle%i_elm.eq.0) return
  ! copy RZPHI coordinates in particles
  particle%x = half_position(4:6)
  ! compute magnetic and electric fields
  call fields%calc_EBpsiU(time+5.d-1*timestep,particle%i_elm,&
       particle%st,particle%x(3),E,B,psi,U)
  ! compute the second half-step  
  call volume_preserving_second_half_step_jorek(particle, &
    half_position(1:3),scaling_factor,                    &
    vector_cylindrical_to_cartesian(particle%x(3),E),     &
    vector_cylindrical_to_cartesian(particle%x(3),B),     &
    mass,timestep)
  ! transform back from cartesian to cylindrical coordinates
  half_position(4:6) = cartesian_to_cylindrical(half_position(1:3))
  ! find the (i_elm,s,t) coordinates
  call find_RZ_nearby(fields%node_list,fields%element_list,particle%x(1), &
    particle%x(2),particle%st(1),particle%st(2),particle%i_elm,           &
    half_position(4),half_position(5),particle%st(1),particle%st(2),      &
    particle%i_elm,ifail)
  ! copy new RZPHI position into particle
  particle%x = half_position(4:6)
end subroutine volume_preserving_push_jorek


!---------------------------------------------------------------------------
!> This subroutine integrates a relativistic particle trajectory in JOREK
!> fields using the Volume Preserving Algorithm (VPA) while also including
!> the effect from the radiation reaction force
subroutine volume_preserving_radiation_push_jorek(particle,fields,mass,time,timestep,ifail)
  ! load functions
  use mod_coordinate_transforms, only: cartesian_to_cylindrical
  use mod_coordinate_transforms, only: cylindrical_to_cartesian
  use mod_coordinate_transforms, only: vector_cylindrical_to_cartesian
  use mod_fields
  use mod_find_rz_nearby
  ! declare input/output variables
  integer(kind=4),intent(inout) :: ifail
  class(particle_kinetic_relativistic), intent(inout) :: particle !< relativistic particle
  ! declare input variables
  real(kind=8),intent(in) :: mass, time, timestep
  class(fields_base), intent(in) :: fields
  ! declare internal variables
  real(kind=8) :: psi, U
  real(kind=8),dimension(3) :: B, E
  ! half_position coordinates: 1:x, 2:y, 3:z, 4:R, 5:Z, 6:phi 
  real(kind=8),dimension(6) :: half_position
  real(kind=8) :: scaling_factor !< in [s^2*C/(kg*m)]

  ! check if the particle is valid
  if(particle%i_elm.eq.0) return
  ! transform the particle position from cylindrical to cartesian coordinates
  half_position(1:3) = cylindrical_to_cartesian(particle%x)
  ! compute first half-step
  call volume_preserving_first_half_step_jorek(particle,half_position(1:3),&
       mass,timestep,scaling_factor)
  ! calculate cylindrical coordinates from cartesian ones
  half_position(4:6) = cartesian_to_cylindrical(half_position(1:3))
  ! find the (i_elm,s,t) coordinates
  call find_RZ_nearby(fields%node_list,fields%element_list,particle%x(1),&
       particle%x(2),particle%st(1),particle%st(2),particle%i_elm,&
       half_position(4),half_position(5),particle%st(1),particle%st(2),&
       particle%i_elm,ifail)
  ! check if the particle is lost, exit if it is the case
  if(particle%i_elm.eq.0) return
  ! copy RZPHI coordinates in particles
  particle%x = half_position(4:6)
  ! compute magnetic and electric fields
  call fields%calc_EBpsiU(time+5.d-1*timestep,particle%i_elm,&
       particle%st,particle%x(3),E,B,psi,U)
  ! compute the second half-step  
  call volume_preserving_second_half_step_jorek(particle, &
    half_position(1:3),scaling_factor,                    &
    vector_cylindrical_to_cartesian(particle%x(3),E),     &
    vector_cylindrical_to_cartesian(particle%x(3),B),     &
    mass,timestep)
  call radreactforce_kinetic(B, timestep, mass, particle)

  ! transform back from cartesian to cylindrical coordinates
  half_position(4:6) = cartesian_to_cylindrical(half_position(1:3))
  ! find the (i_elm,s,t) coordinates
  call find_RZ_nearby(fields%node_list,fields%element_list,particle%x(1), &
    particle%x(2),particle%st(1),particle%st(2),particle%i_elm,           &
    half_position(4),half_position(5),particle%st(1),particle%st(2),      &
    particle%i_elm,ifail)
  ! copy new RZPHI position into particle
  particle%x = half_position(4:6)
end subroutine volume_preserving_radiation_push_jorek

!--------------------------------------------------------------------------

!> This procedure implements a test version of the VPA for analytical
!> magnetic and electric fields. Not for production.
!> inputs:
!>   particle: (particle_kinetic_relativistic) particle to integrate
!>   fields:   (fields_base) jorek fields
!>   mass:     (real8) particle mass
!>   time:     (real8) time integration 
!>   timestep: (real8) time step
!> outputs:
!>   particle: (particle_kinetic_relativistic) integrated particle
pure subroutine volume_preserving_push_analytical(particle,fields,&
     mass,time,timestep)
  use mod_fields
  use mod_coordinate_transforms, only: cartesian_to_cylindrical
  use mod_coordinate_transforms, only: cylindrical_to_cartesian
  use mod_coordinate_transforms, only: vector_cylindrical_to_cartesian
  implicit none
  !> declare input/output varibales
  type(particle_kinetic_relativistic),intent(inout) :: particle
  !> declare inputs
  class(fields_base),intent(in) :: fields
  real(kind=8),intent(in) :: mass,time,timestep
  !> declare variables
  real(kind=8) :: scaling_factor,psi,U
  real(kind=8),dimension(3) :: half_position,E,B

  !> transform the particle position from cylindrical to cartesian
  half_position = cylindrical_to_cartesian(particle%x)
  !> apply the first vpa time step
  call volume_preserving_first_half_step_jorek(particle,half_position,&
       mass,timestep,scaling_factor)
  !> convert back the position to cylindrical coordinates
  particle%x = cartesian_to_cylindrical(half_position)
  !> compute the analytical electromagnetic fields
  call fields%calc_analytical_EBpsiU(particle%x(1:2),E,B,psi,U)
  !> apply the secon VPA step
  call volume_preserving_second_half_step_jorek(particle,&
       half_position,scaling_factor,&
       vector_cylindrical_to_cartesian(particle%x(3),E),&
       vector_cylindrical_to_cartesian(particle%x(3),B),&
       mass,timestep)
  !> convert back the position to cylindrical coordinates
  particle%x = cartesian_to_cylindrical(half_position)
  
end subroutine volume_preserving_push_analytical

!---------------------------------------------------------------------------

!> This subroutine implements a test version of the VPA assuming uniform B and E.
!> This subroutine is to be used for tests and comparisons, not for production.
!> inputs:
!>   particle: (particle_kinetic_relativistic) relativistic particle
!>   mass:     (real8) particle mass in [AMU]
!>   E:        (real8)(3) electric field in [V/m]
!>   B:        (real8)(3) magnetic field in [T]
!>   dt:       (real8) time step in [s]
!> outputs:
!>   particle: (particle_kinetic_relativistic) relativistic particle type
pure subroutine volume_preserving_push_cartesian(particle,mass,E,B,dt)
  use mod_pusher_tools, only: cayley_transform !< use full Cayley transform
  ! define input output variables
  class(particle_kinetic_relativistic), intent(inout) :: particle !< relativistic particle
  ! define input variables
  real(kind=8), intent(in) :: mass, dt !< mass and time step
  real(kind=8), dimension(3), intent(in) :: E, B !< electric and magnetic fields
  ! internal variable
  real(kind=8) :: scaling_factor !< in [s^2*C/(kg*m)]

  scaling_factor = 5.d-1*dt*particle%q*EL_CHG/(ATOMIC_MASS_UNIT*mass*SPEED_OF_LIGHT)

  ! compute dimensionless momentum
  particle%p = particle%p/(mass*SPEED_OF_LIGHT)

  ! compute position at t_(i+1/2)
  particle%x = particle%x + (5.d-1*dt*SPEED_OF_LIGHT*particle%p)/&
    (sqrt(1.d0+dot_product(particle%p,particle%p)))
  
  ! compute momentum at t_(i+1/2)
  particle%p = particle%p + scaling_factor*E
  
  ! rotate momentum with respect to the magnetic field
  particle%p = matmul(cayley_transform(SPEED_OF_LIGHT*scaling_factor/&
    (sqrt(1.d0+dot_product(particle%p,particle%p))),B),particle%p)

  ! compute momentum at t_(i+1)
  particle%p = particle%p + scaling_factor*E

  ! compute position at t_(i+1)
  particle%x = particle%x + (5.d-1*dt*SPEED_OF_LIGHT*particle%p)/&
    (sqrt(1.d0+dot_product(particle%p,particle%p)))

  ! compute dimensional momentum
  particle%p = particle%p*mass*SPEED_OF_LIGHT  
end subroutine volume_preserving_push_cartesian

!----------------------------------------------------------------------

!> This subroutine integrates a relativistic particle trajectory
!> in JOREK fields using the Runge-Kutta integrator
!> inputs
!>  fields:    (fields_base) JOREK fields
!>  t:         (real8) current time
!>  dt:        (real8) time step
!>  mass:      (real8) particle mass
!>  particle:  (particle_kinetic_relativistic) particle to be pushed
!> outputs:
!>   particle: (particle_kinetic_relativistic) pushed particle
subroutine runge_kutta_fixed_dt_relativistic_particle_push_jorek( &
  fields,t,dt,mass,particle)
  !> load modules
  use mod_find_rz_nearby
  use mod_fields, only: fields_base
  use mod_runge_kutta, only: runge_kutta_fixed_dt
  implicit none
  !> declare input/output variables
  type(particle_kinetic_relativistic), intent(inout) :: particle
  !> declare input variables
  class(fields_base), intent(in) :: fields
  real(kind=8), intent(in)       :: t, dt, mass
  !> declare internal variables
  integer                    :: i_elm_new, ifail
  real(kind=8), dimension(2) :: st_new
  real(kind=8), dimension(6) :: solution_new

  !> apply the RK scheme
  call runge_kutta_fixed_dt(                                    &
    compute_relativistic_particle_derivatives_jorek,fields,     &
    6,2,3,t,dt,[particle%x(1),particle%x(2),particle%x(3),      &
    particle%p(1),particle%p(2),particle%p(3)],[particle%i_elm, &
    int(particle%p)],[particle%st(1),particle%st(2),mass],      &
    solution_new,i_elm_new)
  
  !> find the particle in the JOREK mesh in case of success
  if(i_elm_new.ne.0) call find_rz_nearby(fields%node_list,&
       fields%element_list,particle%x(1),particle%x(2),&
       particle%st(1),particle%st(2),particle%i_elm,solution_new(1),&
       solution_new(2),st_new(1),st_new(2),i_elm_new,ifail)

  !> overwrite the particle state
  particle%x     = solution_new(1:3)
  particle%p     = solution_new(4:6)
  particle%st    = st_new
  particle%i_elm = i_elm_new
  
end subroutine runge_kutta_fixed_dt_relativistic_particle_push_jorek

!-----------------------------------------------------------------------

!> This subroutine integrates a relativistic particle trajectory
!> in analytical fields using the Runge-Kutta integrator
!> inputs:
!>   fields:   (fields_base) analytical fields
!>   t:        (real8) current time
!>   dt:       (real8) time step
!>   mass:     (real8) particle mass
!>   particle: (particle_kinetic_relativistic) particle to be pushed
!> outputs:
!>   t:        (real8) new time
!>   particle: (particle_kinetic_relativistic) pushed particle
subroutine runge_kutta_fixed_dt_relativistic_particle_push(fields,t,dt, &
  mass,particle)
  !> load modules
  use mod_fields, only: fields_base
  use mod_runge_kutta, only: runge_kutta_fixed_dt
  implicit none
  !> declare input/output variables
  type(particle_kinetic_relativistic),intent(inout) :: particle
  !> delcare input variables
  class(fields_base), intent(in) :: fields
  real(kind=8), intent(in)       :: t, dt, mass
  !> declare internal variables
  integer                    :: ifail
  real(kind=8), dimension(6) :: solution_new

  call runge_kutta_fixed_dt(compute_relativistic_particle_derivatives, &
    fields,6,1,1,t,dt,[particle%x(1),particle%x(2),particle%x(3),      &
    particle%p(1),particle%p(2),particle%p(3)],[int(particle%q)],      &
    [mass],solution_new,ifail)

  !> overwrite the particle state
  particle%x = solution_new(1:3)
  particle%p = solution_new(4:6)
  
end subroutine runge_kutta_fixed_dt_relativistic_particle_push

!-----------------------------------------------------------------------

!> This procedure computes derivatives for the full orbit particles
!> to be used by the Runge-Kutta integrator in JOREK fields
!> inputs:
!>   fields:            (fields_base) JOREK fields
!>   n_variables:       (integer) number of variables describing the particle = 6
!>   n_int_parameters:  (integer) number of integer parameters = 2
!>   n_real_parameters: (integer) number of real parameters = 3
!>   t:                 (real8) time of current RK step
!>   solution_old:      (real8)(n_variables) particle description at previous time step
!>   solution:          (real8)(n_variables) particle description at current RK step
!>   int_parameters:    (integer)(n_int_parameters) parameters:
!>                      1:i_elm old, 2:charge
!>   real_parameters:   (real8)(n_real_parameters) parameter
!>                      1:s old, 2:t old, 3:mass
!> outputs:
!>   derivatives: (real8)(6) derivatives for RK integrator
!>   ifail:       (integer) if 0 integration failed
subroutine compute_relativistic_particle_derivatives_jorek(fields, &
  n_variables,n_int_parameters,n_real_parameters,t,solution_old,   &
  solution,int_parameters,real_parameters,derivatives,ifail)
  use mod_coordinate_transforms, only: vector_cylindrical_to_cartesian
  use mod_find_rz_nearby
  use mod_fields, only: fields_base
  implicit none
  !> input variables
  class(fields_base), intent(in)                         :: fields
  integer, intent(in)                                    :: n_variables, n_int_parameters, n_real_parameters
  real(kind=8), intent(in)                               :: t
  real(kind=8), dimension(n_variables), intent(in)       :: solution_old, solution
  integer, dimension(n_int_parameters), intent(in)       :: int_parameters
  real(kind=8), dimension(n_real_parameters), intent(in) :: real_parameters
  !> output variables
  real(kind=8), dimension(n_variables), intent(out) :: derivatives
  integer, intent(out)                              :: ifail
  !> internal variables
  integer                    :: ierr
  real(kind=8)               :: psi, U
  real(kind=8), dimension(2) :: st_new !< new local coordinates
  real(kind=8), dimension(3) :: B, E

  !> find the particle local position at current RK step
  call find_rz_nearby(fields%node_list,fields%element_list,solution_old(1),&
       solution_old(2),real_parameters(1),real_parameters(2),&
       int_parameters(1),solution(1),solution(2),st_new(1),&
       st_new(2),ifail,ierr)
       
  !> compute the electromagnetic fields
  if(ifail.ne.0) call fields%calc_EBpsiU(t,ifail,st_new,solution(3),&
       E,B,psi,U)
       
  !> compute RHS terms of relativistic equations of motion needed for RK integration
  derivatives = compute_relativistic_particle_rhs(int_parameters(2), &
    real_parameters(3),solution,                                     &
    vector_cylindrical_to_cartesian(solution(3),E),                  &
    vector_cylindrical_to_cartesian(solution(3),B))
  
end subroutine compute_relativistic_particle_derivatives_jorek

!------------------------------------------------------------------------

!> This procedure computes derivatives for the Runge-Kutta full orbit integrator
!> of relativistic particles in analytical fields.
!> inputs:
!>   fields:           (fields_base) analytical fields
!>   n_variables:      (integer) number of variables describing the particle = 6
!>   n_int_parameters: (integer) number of integer parameters = 1
!>   n_real_paramters: (integer) number of real parameters = 1
!>   t:                (real8) time of current RK step
!>   solution_old:     (real8)(n_variables) particle description at previous time step:
!>                     1:x, 2:y, 3:z, 4:px, 5:py, 6:pz
!>   solution:         (real8)(n_variables) particle description at current RK step
!>   int_parameters:   (integer)(n_int_parameters) integer parameters
!>                     1:charge
!>   real_parameters:  (real8)(n_real_parameters) real parameters
!>                     1:mass
!> outputs:
!>   derivatives:      (real8)(n_variables) solution derivatives
!>   ifail:            (integer) if 0 integration failed
subroutine  compute_relativistic_particle_derivatives(fields,n_variables, &
  n_int_parameters,n_real_parameters,t,solution_old,solution,             &
  int_parameters,real_parameters,derivatives,ifail)
  !> load modules
  use mod_coordinate_transforms, only: vector_cylindrical_to_cartesian
  use mod_fields, only: fields_base
  implicit none
  !> declare input variables
  class(fields_base), intent(in)                         :: fields
  integer, intent(in)                                    :: n_variables, n_int_parameters, n_real_parameters
  real(kind=8), intent(in)                               :: t
  real(kind=8), dimension(n_variables), intent(in)       :: solution, solution_old
  integer, dimension(n_int_parameters), intent(in)       :: int_parameters
  real(kind=8), dimension(n_real_parameters), intent(in) :: real_parameters
  !> declare output variables
  real(kind=8), dimension(n_variables), intent(out) :: derivatives
  integer, intent(out)                              :: ifail
  !> internal variables
  real(kind=8)               :: psi, U
  real(kind=8), dimension(3) :: E, B

  !> compute electromagnetic fields at current RK step
  call fields%calc_analytical_EBpsiU(solution(1:2),E,B,psi,U)

  !> compute RHS terms of relativistic equations of motion needed for RK integration
  derivatives = compute_relativistic_particle_rhs(int_parameters(1), &
    real_parameters(1),solution,                                     &
    vector_cylindrical_to_cartesian(solution(3),E),                  &
    vector_cylindrical_to_cartesian(solution(3),B))
    
  ifail = 1 !< set ifail to true
  
end subroutine compute_relativistic_particle_derivatives

!-------------------------------------------------------------------------

!> This procedure computes the RHS of the full orbit equations of motion
!> inputs:
!>   charge:    (integer) particle charge
!>   mass:      (real8) particle mass
!>   solutions: (real8)(6) particle solution
!>   E:         (real8)(3) electric field [V/m]
!>   B:         (real8)(3) magnetic field [T]
!> outputs:
!>   rhs: (real8)(6) RHS of full orbit equations of motion
pure function compute_relativistic_particle_rhs(charge,mass,solution,E,B) &
  result(rhs)
  !> load modules
  use constants, only: EL_CHG, SPEED_OF_LIGHT, ATOMIC_MASS_UNIT
  use mod_math_operators, only: cross_product
  use mod_coordinate_transforms, only: vector_cartesian_to_cylindrical
  implicit none
  !> declare intput variables
  integer, intent(in)                    :: charge
  real(kind=8), intent(in)               :: mass
  real(kind=8), dimension(6), intent(in) :: solution
  real(kind=8), dimension(3), intent(in) :: E, B
  !> declare output variables
  real(kind=8), dimension(6) :: rhs
  !> internal variables
  real(kind=8) :: gamma !< relativistic factor

  !> compute the relativistic factor
  gamma = sqrt(mass*mass + ((solution(4)*solution(4) + &
    solution(5)*solution(5)+solution(6)*solution(6))/  &
    (SPEED_OF_LIGHT*SPEED_OF_LIGHT)))
    
  !> compute the spatial derivatives
  rhs(1:3) = vector_cartesian_to_cylindrical(solution(3),solution(4:6))/gamma
  rhs(3) = rhs(3)/solution(1)

  !> compute the momenta derivatives
  rhs(4:6) = EL_CHG*real(charge,kind=8)*(      &
    E + cross_product(solution(4:6),B)/gamma)/ &
    ATOMIC_MASS_UNIT
  
end function compute_relativistic_particle_rhs

!--------------------------------------------------------------------------

!> This procedure transforms a particle_kinetic_relativistic into a different
!> particle type
!> inputs:
!>   node_list:    (type_node_list) JOREK node list
!>   element_list: (type_element_list) JOREK element list
!>   particle_in:  (particle_kinetic_relativistic) a relativistic particle
!>   mass:
!>   B:            (real8)(3) magnetic field
!> outputs:
!>   particle_out: (particle_base) output particle
subroutine relativistic_kinetic_to_particle(node_list,element_list,particle_in, &
  particle_out,mass,B)
  !> load modules
  use data_structure
  implicit none
  !> declare input variables
  type(type_node_list), intent(in)                :: node_list
  type(type_element_list), intent(in)             :: element_list
  type(particle_kinetic_relativistic), intent(in) :: particle_in
  real(kind=8), intent(in)                        :: mass
  real(kind=8), dimension(3), intent(in)          :: B
  !> declare output variables
  class(particle_base), intent(out)               :: particle_out

  select type (particle_out)
  type is (particle_gc)
     particle_out = relativistic_kinetic_to_gc(node_list,element_list, &
       particle_in,mass,B)
  type is (particle_gc_relativistic)
     particle_out = relativistic_kinetic_to_relativistic_gc(node_list, &
       element_list,particle_in,mass,B)
  end select
  
end subroutine relativistic_kinetic_to_particle

!--------------------------------------------------------------------------

!> This procedure transforms a particle_kinetic_relativistic 
!> into a particle_gc_relativistic
!> inputs:
!>   node_list:    (type_node_list) JOREK node list
!>   element_list: (type_element_list) JOREK element list
!>   in:           (particle_kinetic_relativistic) relativistic particle
!>   mass:         (real8) particle mass
!>   B:            (real8)(3) magnetic field
!> outputs:
!>   out: (particle_gc_relativistic) a relativistic gc 
function relativistic_kinetic_to_relativistic_gc(node_list,element_list, &
  in,mass,B) result(out)
  !> load modules
  use data_structure
  use mod_coordinate_transforms, only: vector_cartesian_to_cylindrical
  use mod_pusher_tools, only: particle_position_to_gc
  implicit none
  !> declare input variables
  type(type_node_list), intent(in)                :: node_list
  type(type_element_list), intent(in)             :: element_list
  type(particle_kinetic_relativistic), intent(in) :: in
  real(kind=8), intent(in)                        :: mass
  real(kind=8), dimension(3), intent(in)          :: B
  !> declare output variables
  type(particle_gc_relativistic)                  :: out
  !> delcare internal variables
  real(kind=8)                                    :: norm_B
  real(kind=8), dimension(3)                      :: B_hat, p_perp

  !> compute magnetic field direction and intensity
  norm_B = sqrt(B(1)*B(1)+B(2)*B(2)+B(3)*B(3))
  B_hat = B/norm_B

  !> copy base particle fields
  out = in

  !> copy charge
  out%q = in%q

  !> extract momenta in cylindrical coordinates
  p_perp = vector_cartesian_to_cylindrical(in%x(3),in%p)

  !> compute parallel momentum
  out%p(1) = p_perp(1)*B_hat(1)+p_perp(2)*B_hat(2)+p_perp(3)*B_hat(3)

  !> compute perpendicular momenta
  p_perp = p_perp - out%p(1)*B_hat

  !> compute magnetic moment
  out%p(2) = (p_perp(1)*p_perp(1)+p_perp(2)*p_perp(2)+&
       p_perp(3)*p_perp(3))/(2.d0*norm_B*mass)

  !> compute GC position
  if(out%q.ne.0) then
     call particle_position_to_gc(node_list,element_list,in%x,&
          in%st,in%i_elm,vector_cartesian_to_cylindrical(in%x(3),in%p),&
          in%q,B_hat,norm_B,out%x,out%st,out%i_elm)
  endif
  
end function relativistic_kinetic_to_relativistic_gc

!---------------------------------------------------------------------------

!> This function transforms a particle_kinetic_relativistic into a
!> particle_gc. Phase space particle coordinates are transformed
!> into GC position, energy in [eV] and magnetic moment in [eV/T]
!> inputs:
!>   node_list:    (type_node_list) JOREK node list
!>   element_list: (type_element_list) JOREK element list
!>   in:           (particle_relativistic_kinetic) a relativistic particle
!>   mass:         (real8) particle mass in AMU
!>   B:            (real8)(3)(optional) magnetic field [T]
!> outputs:
!>   out: (particle_gc) a guiding center particle
function relativistic_kinetic_to_gc(node_list,element_list,in,mass,B) result(out)
  ! load modules
  use data_structure
  use mod_coordinate_transforms, only: vector_cartesian_to_cylindrical
  use mod_pusher_tools, only: particle_position_to_gc
  ! declare input variables
  type(type_node_list), intent(in)                :: node_list
  type(type_element_list), intent(in)             :: element_list
  type(particle_kinetic_relativistic), intent(in) :: in
  real(kind=8), intent(in) :: mass
  real(kind=8), dimension(3), intent(in) :: B
  ! declare output variables
  type(particle_gc) :: out
  ! declare internal variables
  real(kind=8)               :: B_norm, p_par !< magnetic intensity, parallel momentum
  real(kind=8), dimension(3) :: p_perp, B_hat !< perp. momentum, magnetic field direction

  ! initialise default variables for particle_gc
  out = in
  
  ! copy base particle fields
  out%q = in%q

  ! compute the guiding center total (i.e. rest+kinetic) energy in [eV]
  out%E = ATOMIC_MASS_UNIT*SPEED_OF_LIGHT*                     &
    sqrt((mass*SPEED_OF_LIGHT)*(mass*SPEED_OF_LIGHT)+          &
    (in%p(1)*in%p(1)+in%p(2)*in%p(2)+in%p(3)*in%p(3)))/EL_CHG  

  ! compute magnetic field intensity and direction
  B_norm = sqrt(B(1)*B(1)+B(2)*B(2)+B(3)*B(3)) !< intensity
  B_hat = B/B_norm  !< direction
  
  ! compute the parallel and perpendicular momenta
  p_perp = vector_cartesian_to_cylindrical(in%x(3),in%p)
  p_par  = p_perp(1)*B_hat(1) + p_perp(2)*B_hat(2) + p_perp(3)*B_hat(3)
  p_perp = p_perp - p_par*B_hat

  ! compute the magnetic moment p_perp^2/(2*B) in [eV/T]
  ! the sign is given by the particle parallel momentum 
  out%mu = sign((ATOMIC_MASS_UNIT*(p_perp(1)*p_perp(1)+ &
    p_perp(2)*p_perp(2)+p_perp(3)*p_perp(3))/           &
    (2.d0*B_norm*mass*EL_CHG)),p_par)

  ! compute the GC position
  if(out%q.ne.0) then 
    call particle_position_to_gc(node_list,element_list,  &
      in%x,in%st,in%i_elm,                                &
      vector_cartesian_to_cylindrical(in%x(3),in%p),in%q, &
      B_hat,B_norm,out%x,out%st,out%i_elm)
  endif  
end function relativistic_kinetic_to_gc

!---------------------------------------------------------------------------

!> This function transforms a particle_gc into a particle_kinetic_relativistic.
!> The GC energy should be in [eV] and its magnetic moment in [eV/T].
!> inputs:
!>   node_list:    (type_node_list) JOREK node list
!>   element_list: (type_element_list) JORKE element list
!>   in:           (particle_gc) a guiding center particle
!>   chi:          (real8) gyro-angle in [rad]
!>   time:         (real8) current time
!>   mass:         (real8) particle mass
!>   B:	           (real8)(3)(optional) magnetic field [T]
!> outputs:
!>   out: (particle_kinetic_relativistic) a kinetic relativistic particle
function gc_to_relativistic_kinetic(node_list,element_list,in,time,mass,chi,B) result(out)
  use data_structure
  use mod_coordinate_transforms, only: vector_cylindrical_to_cartesian
  use mod_pusher_tools, only: get_orthonormals
  use mod_pusher_tools, only: gc_position_to_particle
  ! declare input variables
  type(type_node_list), intent(in)       :: node_list
  type(type_element_list), intent(in)    :: element_list
  type(particle_gc), intent(in)          :: in
  real(kind=8), intent(in)               :: time, mass, chi !< time, particle mass [AMU], gyroangle [rad]
  real(kind=8), dimension(3), intent(in) :: B !< magnetic field in [T]
  ! declare output variables
  type(particle_kinetic_relativistic) :: out
  ! declare internal variables
  real(kind=8)               :: B_norm, p_perp, p_par   !< magnetic intensity, perpendicular/parallel momenta
  real(kind=8), dimension(3) :: B_hat, e1, e2 !< B-field-aligned cartesian vector basis

  ! copy base particle fields
  out = in
  
  ! copy particle charge
  out%q = in%q
  
  ! compute the magnetic field intensity and direction
  B_norm = norm2(B)
  B_hat = B/B_norm
  ! compute a B-field-aligned orthonormal vector basis
  call get_orthonormals(B_hat,e1,e2)

  ! compute the perpendicular momentum squared in (AMU*m/s)^2
  ! (put it in p_perp temporarily to save one variable)
  p_perp = (EL_CHG*2.d0*mass*B_norm*abs(in%mu))/ATOMIC_MASS_UNIT

  ! compute the parallel momentum in (AMU*m/s)
  p_par = sign(sqrt((((in%E*EL_CHG)*(in%E*EL_CHG))/                         &
    ((ATOMIC_MASS_UNIT*SPEED_OF_LIGHT)*(ATOMIC_MASS_UNIT*SPEED_OF_LIGHT)))- &
    (mass*SPEED_OF_LIGHT)*(mass*SPEED_OF_LIGHT)-p_perp),in%mu)

  ! compute the perpendicular momentum in (AMU*m/s)
  p_perp = sqrt(p_perp)

  ! computing the particle momentum
  out%p = p_par*B_hat + p_perp*(e1*cos(chi)+e2*sin(chi)) 

  ! compute the particle position in (R,Z,phi) coordinates
  if(out%q.ne.0) then
    call gc_position_to_particle(node_list,element_list,in%x,in%st, &
      in%i_elm,out%p,in%q,B_hat,B_norm,out%x,out%st,out%i_elm)
  endif
 
  !> transform the momentum in cartesian coodinates
  out%p = vector_cylindrical_to_cartesian(out%x(3),out%p)
end function gc_to_relativistic_kinetic


end module mod_kinetic_relativistic
