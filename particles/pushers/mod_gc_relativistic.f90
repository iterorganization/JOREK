!> this module contains procedures for pushing and transforming
!> relativistic guiding centers
module mod_gc_relativistic
  !> load modules
  use mod_particle_types
  use mod_radreactforce
  
  implicit none
  
  !> declare default private
  private
  !> declare public procedures and variables
  public relativistic_gc_to_particle
  public gc_to_relativistic_gc
  public compute_relativistic_factor
  public relativistic_gc_to_relativistic_kinetic
  public relativistic_gc_momenta_from_E_cospitch
  public runge_kutta_fixed_dt_gc_push_jorek
  public runge_kutta_fixed_dt_gc_push_jorek_radreact
  public runge_kutta_fixed_dt_gc_push
  public runge_kutta_adapt_dt_gc_push_jorek
  public runge_kutta_error_control_dt_gc_push_jorek

contains

  !> This subroutine pushes a relativistic guiding center in JOREK fields
  !> using a Runge-Kutta integrator in which the time step is adjusted to control 
  !> the integration error. At the moment, the error metric is hard-coded, 
  !> but a user-defined metric can be easily added as last argument to the 
  !> procedure runge_kutta_order_error_control_dt.
  !> inputs:
  !>   fields:     (fields_base) JOREK fields
  !>   tolerances: (real8)(4) tolerances on the integration error for
  !>                          1:R, 2:Z, 3:phi, 4:p_parallel
  !>   t:          (real8) current time
  !>   dt:         (real8) first time step that will be tried
  !>   t_stop:     (real8) time of next event
  !>   mass:       (real8) particle mass in AMU
  !>   particle:   (particle_gc_relativistic) GC to be pushed
  !> outputs:
  !>   dt:       (real8) time step used for push
  !>   dt_new:   (real8) time step suggested for the next push
  !>   particle: (particle_gc_relativistic) pushed GC
  !>   error:    (real8) final Runge-Kutta error (if defined TEST only)
#ifdef TEST
  subroutine runge_kutta_error_control_dt_gc_push_jorek(fields,tolerances, &
    t,dt,t_stop,mass,dt_new,particle,error)
#else
  subroutine runge_kutta_error_control_dt_gc_push_jorek(fields,tolerances, &
    t,dt,t_stop,mass,dt_new,particle)
#endif
    !> modules
    use mod_fields, only: fields_base
    use mod_find_rz_nearby
    use mod_runge_kutta, only: runge_kutta_order_error_control_dt
    implicit none
    !> input/output variables
    type(particle_gc_relativistic), intent(inout) :: particle
    real(kind=8), intent(inout)                   :: dt
    !> input variables
    class(fields_base), intent(in)         :: fields
    real(kind=8), intent(in)               :: t, t_stop, mass
    real(kind=8), dimension(4), intent(in) :: tolerances
    !> output variables
    real(kind=8), intent(out) :: dt_new
#ifdef TEST
    real(kind=8), intent(out) :: error
#endif
    !> internal variables
    integer :: ifail, i_elm_new
    real(kind=8), dimension(2) :: st_new
    real(kind=8), dimension(4) :: solution_new

    !> compute Runge-Kutta solution and new time step
#ifdef TEST
    call runge_kutta_order_error_control_dt(compute_relativistic_gc_derivatives_jorek, &
      fields,4,2,4,t,t_stop,dt,[particle%x(1),particle%x(2),particle%x(3),             &
      particle%p(1)],[particle%i_elm,int(particle%q)],[particle%st(1),                 &
      particle%st(2),mass,particle%p(2)],                                              &
      tolerances,dt_new,solution_new,i_elm_new,error)
#else
    call runge_kutta_order_error_control_dt(compute_relativistic_gc_derivatives_jorek, &
      fields,4,2,4,t,t_stop,dt,[particle%x(1),particle%x(2),particle%x(3),             &
      particle%p(1)],[particle%i_elm,int(particle%q)],[particle%st(1),                 &
      particle%st(2),mass,particle%p(2)],                                              &
      tolerances,dt_new,solution_new,i_elm_new)
#endif

    !> compute the new local coordinates
    particle%i_elm = i_elm_new
    if(i_elm_new.gt.0) then 
      call find_rz_nearby(fields%node_list,                  &
        fields%element_list,particle%x(1),particle%x(2),     &
        particle%st(1),particle%st(2),particle%i_elm,        &
        solution_new(1),solution_new(2),st_new(1),st_new(2), &
        i_elm_new,ifail)
      particle%st = st_new       
    endif

    !> overwrite particle fields
    particle%x    = solution_new(1:3)
    particle%p(1) = solution_new(4)
    
  end subroutine runge_kutta_error_control_dt_gc_push_jorek

  !> This subroutine pushes a relativistic guiding center in JOREK fields
  !> using a Runge-Kutta integrator with an adaptation of the time step
  !> before the push based on local gradients.
  !> inputs:
  !>   fields:   (fields_base) JOREK fields
  !>   t:        (real8) current time
  !>   dt:       (real8) "suggested" time step
  !>   t_stop:   (real8) time of next event
  !>   mass:     (real8) particle mass in AMU
  !>   particle: (particle_gc_relativistic) GC to be pushed
  !> outputs:
  !>   t:        (real8) new time
  !>   dt:       (real8) time step used for push
  !>   particle: (particle_gc_relativistic) pushed GC
  subroutine runge_kutta_adapt_dt_gc_push_jorek(fields,t,dt, &
    t_stop,mass,particle)
    !> modules
    use mod_fields, only: fields_base
    use mod_find_rz_nearby
    use mod_runge_kutta, only: runge_kutta_adaptative_dt
    implicit none
    !> input/output variables
    type(particle_gc_relativistic), intent(inout) :: particle
    real(kind=8), intent(inout) :: t, dt
    !> input variables
    class(fields_base), intent(in) :: fields
    real(kind=8), intent(in)       :: t_stop,mass
    !> internal variables
    integer :: ifail, i_elm_new
    real(kind=8), dimension(2) :: st_new
    real(kind=8), dimension(4) :: solution_new

    !> compute Runge-Kutta solution and new time step
    call runge_kutta_adaptative_dt(compute_relativistic_gc_derivatives_jorek,&
         adapt_time_step_gradB_curlb_dbdt,fields,4,2,4,&
         t,t_stop,dt,[particle%x(1),particle%x(2),particle%x(3),&
         particle%p(1)],[particle%i_elm,int(particle%q)],[particle%st(1),&
         particle%st(2),mass,particle%p(2)],solution_new,i_elm_new)
    
    !> compute new local coordinates
    if(i_elm_new.ne.0) call find_rz_nearby(fields%node_list,&
         fields%element_list,particle%x(1),particle%x(2),&
         particle%st(1),particle%st(2),particle%i_elm,&
         solution_new(1),solution_new(2),st_new(1),st_new(2),&
         i_elm_new,ifail)
    
    !> overwrite particle fields
    particle%x     = solution_new(1:3)
    particle%p(1)  = solution_new(4)
    particle%st    = st_new
    particle%i_elm = i_elm_new
    
  end subroutine runge_kutta_adapt_dt_gc_push_jorek


   !> This function adapts the time step as a function as a function of the
   !> gradient and curvature length
   !> inputs:
   !>   fields:            (field_base)(optional) fields for particle pushing
   !>   n_variables:       (integer) number of variables describing the particle = 4
   !>   n_int_parameters:  (integer) number of integer parameters = 2
   !>   n_real_parameters: (integer) number of real parameters = 4
   !>   t:                 (real8) time
   !>   dt:                (real8) "suggested" time step
   !>   solution:          (real8)(n_variables) solution 1:R, 2:Z, 3:phi, 4:p_parallel
   !>   int_parameters:    (integer)(2) integer parameters 1:i_elm, 2:q
   !>   real_parameters:   (real8)(4) real parameters 1:s, 2:t, 3:mass, 4:magnetic moment
   !> outputs:
   !>   dt_new: (real8) adapted time step
  pure function adapt_time_step_gradB_curlb_dbdt(fields,n_variables,  &
    n_int_parameters, n_real_parameters,t,dt,solution,int_parameters, &
    real_parameters) result(dt_new)
     !> modules
     use constants, only: EL_CHG, ATOMIC_MASS_UNIT, TWOPI, SPEED_OF_LIGHT
     use mod_fields, only: fields_base
     implicit none
     !> parameters
     real(kind=8), parameter :: inv_number_steps=0.1 !< 10 steps per characteristic time
     real(kind=8), parameter :: minimum_n_gyroperiod=1.d0
     !> input variables
     class(fields_base), intent(in)                         :: fields
     integer, intent(in)                                    :: n_variables, n_int_parameters, n_real_parameters
     real(kind=8), intent(in)                               :: t, dt
     real(kind=8), dimension(n_variables), intent(in)       :: solution
     integer, dimension(n_int_parameters), intent(in)       :: int_parameters
     real(kind=8), dimension(n_real_parameters), intent(in) :: real_parameters
     !> output variables
     real(kind=8) :: dt_new
     !> internal variables
     real(kind=8)               :: velocity, relativistic_factor, normB
     real(kind=8), dimension(3) :: E, b, gradB, curlb, dbdt 

     !> compute fields
     call fields%calc_EBNormBGradBCurlbDbdt(t,int_parameters(1),    &
       real_parameters(1:2),solution(3),E,b,normB,gradB,curlb,dbdt)
     
     !> computing the particle velocity
     velocity = (solution(4)*solution(4) +                &
       2.d0*real_parameters(3)*real_parameters(4)*normB)/ &
       (real_parameters(3)*real_parameters(3))
     relativistic_factor = sqrt(1.d0 + velocity/(SPEED_OF_LIGHT*SPEED_OF_LIGHT))
     velocity = sqrt(velocity)/relativistic_factor

     !> compute the characteristic times
     gradB = normB/(velocity*gradB)
     curlb = 1.d0/(velocity*curlb)
     dbdt = 1.d0/dbdt

     !> compute the adapted time step
     dt_new = inv_number_steps*minval(abs([gradB(1),gradB(2),gradB(3), &
       curlb(1),curlb(2),curlb(3),dbdt(1),dbdt(2),dbdt(3)]))

     !> limit the time step to a fixed number of gyroperiods
     dt_new = max(dt_new,                                  &
       minimum_n_gyroperiod*TWOPI*abs((real_parameters(3)* &
       ATOMIC_MASS_UNIT*relativistic_factor)/              &
       (real(int_parameters(2),kind=8)*normB)))

   end function adapt_time_step_gradB_curlb_dbdt

  !> This subroutine pushes a relativistic guiding center in JOREK fields
  !> using a standard Runge-Kutta integrator without time step control.
  !> inputs:
  !>   fields:   (fields_base) JOREK fields
  !>   t:        (real8) current time
  !>   dt:       (real8) time step
  !>   mass:     (real8) GC mass in AMU
  !>   particle: (particle_gc_relativistic) GC to be pushed
  !> outputs:
  !>   t:        (real8) new time 
  !>   particle: (particle_gc_relativistic) pushed GC
  subroutine runge_kutta_fixed_dt_gc_push_jorek(fields,t,dt,mass,particle)
    !> modules
    use mod_fields, only: fields_base
    use mod_find_rz_nearby
    use mod_runge_kutta, only: runge_kutta_fixed_dt
    implicit none
    !> input/output variables
    type(particle_gc_relativistic), intent(inout) :: particle
    real(kind=8), intent(inout)                   :: t
    !> input variables
    class(fields_base), intent(in) :: fields
    real(kind=8), intent(in)       :: dt, mass
    !> internal variables
    integer                    :: ifail, i_elm_new 
    real(kind=8), dimension(2) :: st_new
    !> global coordinates used during RK integration: 1:R, 2:Z, 3:phi, 4:p_parallel
    real(kind=8), dimension(4) :: solution_new
    !> global coordinates when RR is included includes also 5:mu
    real(kind=8), dimension(5) :: solution_new_rr

    call runge_kutta_fixed_dt(compute_relativistic_gc_derivatives_jorek, &
         fields,4,2,4,t,dt,[particle%x(1),particle%x(2),particle%x(3),      &
         particle%p(1)],[particle%i_elm,int(particle%q)],[particle%st(1),   &
         particle%st(2),mass,particle%p(2)],solution_new,i_elm_new)

    !> compute the new local coordinates
    if(i_elm_new.ne.0) call find_rz_nearby(fields%node_list, &
         fields%element_list,particle%x(1),particle%x(2),       &
         particle%st(1),particle%st(2),particle%i_elm,          &
         solution_new(1),solution_new(2),st_new(1),st_new(2),   &
         i_elm_new,ifail)
    
    !> overwrite GC fields
    particle%x     = solution_new(1:3)
    particle%p(1)  = solution_new(4)
    particle%st    = st_new
    particle%i_elm = i_elm_new
    
  end subroutine runge_kutta_fixed_dt_gc_push_jorek


  !> This subroutine pushes a relativistic guiding center in JOREK fields
  !> using a standard Runge-Kutta integrator without time step control. This pusher
  !> takes into account the radiation reaction force.
  !> inputs:
  !>   fields:   (fields_base) JOREK fields
  !>   t:        (real8) current time
  !>   dt:       (real8) time step
  !>   mass:     (real8) GC mass in AMU
  !>   particle: (particle_gc_relativistic) GC to be pushed
  !> outputs:
  !>   t:        (real8) new time 
  !>   particle: (particle_gc_relativistic) pushed GC
  subroutine runge_kutta_fixed_dt_gc_push_jorek_radreact(fields,t,dt,mass,particle)
    !> modules
    use mod_fields, only: fields_base
    use mod_find_rz_nearby
    use mod_runge_kutta, only: runge_kutta_fixed_dt
    implicit none
    !> input/output variables
    type(particle_gc_relativistic), intent(inout) :: particle
    real(kind=8), intent(inout)                   :: t
    !> input variables
    class(fields_base), intent(in) :: fields
    real(kind=8), intent(in)       :: dt, mass
    !> internal variables
    integer                    :: ifail, i_elm_new 
    real(kind=8), dimension(2) :: st_new
    !> global coordinates used during RK integration: 1:R, 2:Z, 3:phi, 4:p_parallel
    real(kind=8), dimension(4) :: solution_new
    !> global coordinates when RR is included includes also 5:mu
    real(kind=8), dimension(5) :: solution_new_rr

    !> compute Runge-Kutta differentials
    call runge_kutta_fixed_dt(compute_relativistic_gc_derivatives_jorek_radreactionforce, &
         fields,5,2,3,t,dt,[particle%x(1),particle%x(2),particle%x(3),      &
         particle%p(1),particle%p(2)],[particle%i_elm,int(particle%q)],[particle%st(1),   &
         particle%st(2),mass],solution_new_rr,i_elm_new)
       
    !> compute the new local coordinates
    if(i_elm_new.ne.0) call find_rz_nearby(fields%node_list, &
         fields%element_list,particle%x(1),particle%x(2),       &
         particle%st(1),particle%st(2),particle%i_elm,          &
         solution_new_rr(1),solution_new_rr(2),st_new(1),st_new(2),   &
         i_elm_new,ifail)
    
    !> overwrite GC fields
    particle%x     = solution_new_rr(1:3)
    particle%p(1)  = solution_new_rr(4)
    particle%p(2)  = solution_new_rr(5)
    particle%st    = st_new
    particle%i_elm = i_elm_new
    
  end subroutine runge_kutta_fixed_dt_gc_push_jorek_radreact

  !> This procedure pushes a relativistic guiding center in analytical fields
  !> using a standard Runge-Kutta integrator without time step control
  !> inputs:
  !>   fields:   (fields_base) analytical fields
  !>   t:        (real8) current time
  !>   dt:       (real8) time step
  !>   mass:     (real8) particle mass in AMU
  !>   particle: (particle_gc_relativistic) GC to be pushed
  !> outputs:
  !>   t:        (real8) updated time
  !>   particle: (particle_gc_relativistic) pushed particle
  subroutine runge_kutta_fixed_dt_gc_push(fields,t,dt,mass,particle)
    !> load modules
    use mod_fields, only: fields_base
    use mod_runge_kutta, only: runge_kutta_fixed_dt
    implicit none
    !> declare input/output varibales
    real(kind=8),intent(inout) :: t
    type(particle_gc_relativistic), intent(inout) :: particle
    !> declare input variables
    class(fields_base), intent(in) :: fields
    real(kind=8), intent(in)       :: dt, mass
    !> declare internal variables
    integer                    :: ifail
    real(kind=8), dimension(4) :: solution_new

    !> push GC
    call runge_kutta_fixed_dt(compute_relativistic_gc_derivatives,&
         fields,4,1,2,t,dt,[particle%x(1),particle%x(2),particle%x(3),&
         particle%p(1)],[int(particle%q)],[mass,particle%p(2)],&
         solution_new,ifail)
	 
    !> overwrite GC fields
    particle%x    = solution_new(1:3)
    particle%p(1) = solution_new(4)
    
  end subroutine runge_kutta_fixed_dt_gc_push
  
  !> This procedure computes the derivatives at a given Runge-Kutta step required  
  !> for the Runge-Kutta integration of a relativistic guiding center in JOREK fields.
  !> inputs:
  !>   fields:            (fields_base) JOREK fields
  !>   n_variables:       (integer) number of variables describing the GC = 4
  !>   n_int_parameters:  (integer) number of integer parameters = 2
  !>   n_real_parameters: (integer) number of real parameters = 4
  !>   t:                 (real8) time at the current RK step
  !>   solution_old:      (real8)(n_variables) GC at previous time step:
  !>                      1:R, 2:Z, 3:phi, 4:parallel momentum
  !>   solution:          (real8)(n_variables) GC at current Runge-Kutta step
  !>   int_parameters:    (integer)(n_int_parameters) integer parameters
  !>                      1:old i_elm, 2:charge
  !>   real_parameters:   (real8)(n_real_parameters) real parameters:
  !>                      1:s_old, 2:t_old, 3:mass, 4:magnetic moment
  !> outputs:
  !>   derivatives:       (real8)(n_variables) Runge-Kutta derivatives
  !>   ifail:             (integer) if 0 calculation failed
  subroutine compute_relativistic_gc_derivatives_jorek(fields,n_variables, &
    n_int_parameters,n_real_parameters,t,solution_old,solution,            &
    int_parameters,real_parameters,derivatives,ifail)
    !> load modules
    use constants, only: SPEED_OF_LIGHT, EL_CHG, ATOMIC_MASS_UNIT
    use mod_fields, only: fields_base
    use mod_math_operators, only: cross_product
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
    integer, intent(out)                              :: ifail
    real(kind=8), dimension(n_variables), intent(out) :: derivatives
    !> internal variables
    integer                    :: ierr !< error for find_rz nearby
    real(kind=8)               :: normB, gamma !< magnetic field intensity and relativistic factor
    real(kind=8), dimension(2) :: st_new !< local GC postion at current RK step
    !> fields required to push the relativistic GC
    real(kind=8), dimension(3) :: E, b, gradB, curlb, dbdt, B_star

    !> find GC at current RK step
    call find_RZ_nearby(fields%node_list,fields%element_list,solution_old(1), &
      solution_old(2),real_parameters(1),real_parameters(2),                  &
      int_parameters(1),solution(1),solution(2),st_new(1),                    &
      st_new(2),ifail,ierr)

    !> compute required fields at GC position
    if(ifail.ne.0) call fields%calc_EBNormBGradBCurlbDbdt(t,ifail,st_new, &
      solution(3),E,b,normB,gradB,curlb,dbdt)

    !> compute RHS of GC evolution (Cary-Brizard) equations
    derivatives = compute_relativistic_gc_rhs(int_parameters(2), &
      real_parameters(3),real_parameters(4),solution(1),         &
      solution(4),normB,E,b,gradB,curlb,dbdt)

  end subroutine compute_relativistic_gc_derivatives_jorek

  !> This procedure computes the derivatives at a given Runge-Kutta step required  
  !> for the Runge-Kutta integration of a relativistic guiding center in JOREK fields
  !> when the radiation reaction force is included. Note that here magnetic moment
  !> is no longer a parameter but variable.
  !>
  !> inputs:
  !>   fields:            (fields_base) JOREK fields
  !>   n_variables:       (integer) number of variables describing the GC = 5
  !>   n_int_parameters:  (integer) number of integer parameters = 2
  !>   n_real_parameters: (integer) number of real parameters = 3
  !>   t:                 (real8) time at the current RK step
  !>   solution_old:      (real8)(n_variables) GC at previous time step:
  !>                      1:R, 2:Z, 3:phi, 4:parallel momentum, 5:magnetic moment
  !>   solution:          (real8)(n_variables) GC at current Runge-Kutta step
  !>   int_parameters:    (integer)(n_int_parameters) integer parameters
  !>                      1:old i_elm, 2:charge
  !>   real_parameters:   (real8)(n_real_parameters) real parameters:
  !>                      1:s_old, 2:t_old, 3:mass
  !> outputs:
  !>   derivatives:       (real8)(n_variables) Runge-Kutta derivatives
  !>   ifail:             (integer) if 0 calculation failed
  subroutine compute_relativistic_gc_derivatives_jorek_radreactionforce(fields,n_variables, &
    n_int_parameters,n_real_parameters,t,solution_old,solution,            &
    int_parameters,real_parameters,derivatives,ifail)
    !> load modules
    use constants, only: SPEED_OF_LIGHT, EL_CHG, ATOMIC_MASS_UNIT
    use mod_fields, only: fields_base
    use mod_math_operators, only: cross_product
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
    integer, intent(out)                              :: ifail
    real(kind=8), dimension(n_variables), intent(out) :: derivatives
    !> internal variables
    integer                    :: ierr !< error for find_rz nearby
    real(kind=8)               :: normB, gamma !< magnetic field intensity and relativistic factor
    real(kind=8), dimension(2) :: st_new !< local GC postion at current RK step
    real(kind=8), dimension(2) :: derivatives_rr !< derivatives from radiation reaction force (ppardot, mudot)
    real(kind=8), dimension(4) :: derivatives_h  !< derivatives from Hamiltonian motion (Rdot,zdor,phidot,ppardot)
    !> fields required to push the relativistic GC
    real(kind=8), dimension(3) :: E, b, gradB, curlb, dbdt, B_star

    !> find GC at current RK step
    call find_RZ_nearby(fields%node_list,fields%element_list,solution_old(1), &
      solution_old(2),real_parameters(1),real_parameters(2),                  &
      int_parameters(1),solution(1),solution(2),st_new(1),                    &
      st_new(2),ifail,ierr)

    !> compute required fields at GC position
    if(ifail.ne.0) call fields%calc_EBNormBGradBCurlbDbdt(t,ifail,st_new, &
      solution(3),E,b,normB,gradB,curlb,dbdt)

    !> compute RHS of GC evolution (Cary-Brizard) equations
    derivatives_h = compute_relativistic_gc_rhs(int_parameters(2), &
      real_parameters(3),solution(5),solution(1),         &
      solution(4),normB,E,b,gradB,curlb,dbdt)

    !> compute RHS of radiation reaction force and add it to the Hamiltonian  motion
    call radreactforce_gc_rhs(B*normB, real_parameters(3), int_parameters(2), (/solution(4), solution(5)/), derivatives_rr)

    derivatives(1:3) = derivatives_h(1:3)
    derivatives(4)   = derivatives_h(4) + derivatives_rr(1)
    derivatives(5)   = derivatives_rr(2)

  end subroutine compute_relativistic_gc_derivatives_jorek_radreactionforce

  !> This procedure computes the derivatives at a given Runge-Kutta step required  
  !> for the Runge-Kutta integration of a relativistic guiding center in analytical fields.
  !> This is mainly used for testing models.
  !> inputs:
  !>   fields:            (fields_base) analytical fields
  !>   n_variables:       (n_variables) number of variables describing the GC
  !>   n_int_parameters:  (integer) number of integer parameters = 1
  !>   n_real_parameters: (integer) number of real parameters = 2
  !>   t:                 (real8) time at the current RK step
  !>   solution_old:      (real8)(n_variables) GC at previous time step
  !>   solution:          (real8)(n_variables) GC at current Runge-Kutta step
  !>   int_parameters:    (integer)(n_int_parameters) 1: charge
  !>   real_parameters:   (real8)(n_real_parameters) 1:mass, 2:magnetic moment
  !> outputs:
  !>   ifail:       (integer) if 0 the integration failed
  !>   derivatives: (real8)(n_variables) Runge-Kutta derivatives
  pure subroutine compute_relativistic_gc_derivatives(fields,n_variables, &
    n_int_parameters,n_real_parameters,t,solution_old,solution,           &
    int_parameters,real_parameters,derivatives,ifail)
    !> load modules
    use mod_fields, only: fields_base
    implicit none
    !> declare input variables
    class(fields_base), intent(in)                         :: fields
    integer, intent(in)                                    :: n_variables, n_int_parameters, n_real_parameters
    real(kind=8), intent(in)                               :: t
    real(kind=8), dimension(n_variables), intent(in)       :: solution_old, solution
    integer, dimension(n_int_parameters), intent(in)       :: int_parameters
    real(kind=8), dimension(n_real_parameters), intent(in) :: real_parameters
    !> declare output variables
    integer, intent(out)                              :: ifail
    real(kind=8), dimension(n_variables), intent(out) :: derivatives
    !> internal variables
    real(kind=8)               :: normB
    real(kind=8), dimension(3) :: E, b, gradB, curlb, dbdt

    !> compute required fields at GC position
    call fields%calc_analytical_EBNormBGradBCurlbDbdt(solution(1:2),E,b,normB, &
         gradB,curlb,dbdt)
	 
    !> compute RHS of GC evolution (Cary-Brizard) equations
    derivatives = compute_relativistic_gc_rhs(int_parameters(1),real_parameters(1), &
         real_parameters(2),solution(1),solution(4),normB,E,b,gradB,curlb,dbdt)

    !> set ifail to true
    ifail = 1
    
  end subroutine compute_relativistic_gc_derivatives
  
  !> This procedure computes the RHS of the relativistic GC evolution equations 
  !> as given in J.R Cary, A.J. Brizard, Rev. Mod. Phys, vol.81, p.693, 2009
  !> inputs:
  !> outputs:
  !>   derivatives: (real8)(4) gc right hand side: 1:R_dot, 2:Z_dot,
  !>   3:phi_dot, 4:p_parallel_dot
  pure function compute_relativistic_gc_rhs(charge,mass,magnetic_moment, &
    R,p_parallel,normB,E,b,gradB,curlb,dbdt) result(derivatives)
    use constants, only: EL_CHG, SPEED_OF_LIGHT, ATOMIC_MASS_UNIT
    use mod_math_operators, only: cross_product
    implicit none
    !> declare input variables
    integer,intent(in)                     :: charge
    real(kind=8),intent(in)                :: mass, R, p_parallel, magnetic_moment, normB
    real(kind=8), dimension(3), intent(in) :: b, E, gradB, curlb, dbdt
    !> declare output variable
    real(kind=8), dimension(4) :: derivatives
    !> declare input variables
    real(kind=8)               :: gamma !< relativistic factor
    real(kind=8), dimension(3) :: B_star, E_star

    gamma = sqrt(mass*mass*SPEED_OF_LIGHT*SPEED_OF_LIGHT +    &
      p_parallel*p_parallel+2.d0*mass*normB*magnetic_moment)/ &
      (mass*SPEED_OF_LIGHT)

    B_star = p_parallel*curlb +                               &
      ((EL_CHG*real(charge,kind=8)*normB)/ATOMIC_MASS_UNIT)*b

    E_star = (EL_CHG*real(charge,kind=8)*E/ATOMIC_MASS_UNIT) - &
      p_parallel*dbdt - ((magnetic_moment*gradB)/gamma)
    
    derivatives(1:3) = cross_product(E_star,b) + ((p_parallel*B_star)/(mass*gamma))
    derivatives(3) = derivatives(3)/R
    derivatives(4) = B_star(1)*E_star(1) + B_star(2)*E_star(2) + B_star(3)*E_star(3)
    derivatives = derivatives/(B_star(1)*b(1)+B_star(2)*b(2)+B_star(3)*b(3))
    
  end function compute_relativistic_gc_rhs
  
  !> This procedure transforms a relativistic GC into a different type
  !> inputs:
  !>   node_list:       (type_node_list) JOREK node list
  !>   element_list:    (type_element_list) JOREK element list  
  !>   relativistic_gc: (particle_gc_relativistic) a relativistic GC
  !>   time:            (real8) current time
  !>   mass:            (real8) particle mass
  !>   B:               (real8)(3) magnetic field
  !>   gyro_angle:      (real8)(optional) the gyro-angle, defaut=0
  !> outputs:
  !>   particle_out: (particle_base) the output particle
  !>   gyro_angle:   (real8)(optional) the gyro-angle, default=0
  subroutine relativistic_gc_to_particle(node_list,element_list, &
    relativistic_gc,particle_out,mass,B,gyro_angle)
    !> load modules
    use data_structure
    implicit none
    !> declare input variables
    type(type_node_list), intent(in)           :: node_list
    type(type_element_list), intent(in)        :: element_list
    type(particle_gc_relativistic), intent(in) :: relativistic_gc
    real(kind=8), intent(in)                   :: mass
    real(kind=8), dimension(3), intent(in)     :: B
    real(kind=8), intent(in), optional         :: gyro_angle
    !> declare output variables
    class(particle_base), intent(out)          :: particle_out
    !> declare internal variables
    real(kind=8)                               :: gyro_angle_local

    gyro_angle_local = 0.d0
    if(present(gyro_angle)) gyro_angle_local = gyro_angle
        
    !> select particle type
    select type (particle_out)
    type is (particle_gc)
       particle_out = relativistic_gc_to_gc(relativistic_gc,mass,B)
    type is (particle_kinetic_relativistic)
       particle_out = relativistic_gc_to_relativistic_kinetic(              &
          node_list,element_list,relativistic_gc,mass,B,gyro_angle_local)
    end select
    
  end subroutine relativistic_gc_to_particle

  !> This procedure transforms a particle_gc_relativistic 
  !> into a particle_kinetic_relativistic
  !> inputs:
  !>   node_list:       (type_node_list) JOREK node list
  !>   element_list:    (type_element_list) JOREK element list
  !>   relativistic_gc: (particle_gc_relativistic) a relativistic GC
  !>   mass:            (real8) particle mass 
  !>   B:               (real8)(3) magnetic field
  !>   gyro_angle:      (real8) gyro-angle for initialising the particle
  !> outputs:
  !>   relativistic_particle: (particle_kinetic_relativistic) a relativistic kinetic particle
  function relativistic_gc_to_relativistic_kinetic(node_list,element_list, &
    in,mass,B,gyro_angle) result(out)
    !> load modules
    use data_structure
    use mod_coordinate_transforms, only: vector_cylindrical_to_cartesian
    use mod_pusher_tools, only: gc_position_to_particle,get_orthonormals
    implicit none
    !> declare input varibales
    type(type_node_list), intent(in)           :: node_list
    type(type_element_list), intent(in)        :: element_list
    type(particle_gc_relativistic), intent(in) :: in
    real(kind=8), intent(in)                   :: mass, gyro_angle
    real(kind=8), dimension(3), intent(in)     :: B
    !> declare output variables
    type(particle_kinetic_relativistic)        :: out
    !> declare internal variables
    real(kind=8)                               :: B_norm, p_perp
    real(kind=8), dimension(3)                 :: B_hat, e1, e2
    
    out = in !< copy fields common to both particle types
    out%q = in%q !< copy the charge
    
    !> compute the magnetic field intensity and direction
    B_norm = norm2(B)
    B_hat = B/B_norm
    !> construct orthogonal basis
    call get_orthonormals(B_hat,e1,e2)
    !> compute the perpendicular momentum
    p_perp = sqrt(2.d0*mass*B_norm*in%p(2))
    !> compute particle momenta
    out%p = in%p(1)*B_hat + p_perp*(e1*cos(gyro_angle)+e2*sin(gyro_angle))
    if(out%q.ne.0) then
       call gc_position_to_particle(node_list,element_list,in%x, &
         in%st,in%i_elm,out%p,out%q,B_hat,B_norm,out%x,          &
         out%st,out%i_elm)
    endif
    !> transform the momenta to cartesian coordinates
    out%p = vector_cylindrical_to_cartesian(out%x(3),out%p)
  end function relativistic_gc_to_relativistic_kinetic 

  !> This procedure transforms a particle_gc_relativistic into a particle_gc
  !> inputs:
  !>   relativistic_gc: (particle_gc_relativistic) relativistic GC
  !>   mass:            (real8) particle mass
  !>   B:               (real8)(3) magnetic field
  !> outputs:
  !>   particle_out:    (particle_gc) GC
  pure function relativistic_gc_to_gc(relativistic_gc,mass,B) &
    result(particle_out)
    !> load modules
    use constants, only: EL_CHG, ATOMIC_MASS_UNIT, SPEED_OF_LIGHT
    implicit none
    !> declare inputs:
    type(particle_gc_relativistic), intent(in) :: relativistic_gc
    real(kind=8), intent(in)                   :: mass
    real(kind=8), dimension(3), intent(in)     :: B
    !> declare output variable
    type(particle_gc)                          :: particle_out

    particle_out%x = relativistic_gc%x
    particle_out%q = relativistic_gc%q
    particle_out%i_elm = relativistic_gc%i_elm
    particle_out%st = relativistic_gc%st
    
    !> copy the magnetic moment in eV/T with p_parallel sign
    particle_out%mu = sign(ATOMIC_MASS_UNIT*relativistic_gc%p(2)/EL_CHG, &
         relativistic_gc%p(1))
	 
    !> compute the guiding center energy in eV
    particle_out%E = ATOMIC_MASS_UNIT*SPEED_OF_LIGHT          &
         * sqrt(mass*mass*SPEED_OF_LIGHT*SPEED_OF_LIGHT +     &
                relativistic_gc%p(1)*relativistic_gc%p(1) +   &
                2.d0*mass*norm2(B)*relativistic_gc%p(2)) /EL_CHG
  end function relativistic_gc_to_gc
  
  !> This procedure transforms a particle_gc into particle_gc_relativistic
  !> inputs:
  !>   gc_in: (particle_gc) GC
  !>   mass:  (real8) particle mass
  !>   B:     (real8)(3) magnetic field
  !> outputs:
  !>   relativistic_gc: (particle_gc_relativistic) relativistic gc
  pure function gc_to_relativistic_gc(gc_in,mass,B) result(relativistic_gc)
    !> load modules
    use constants, only: EL_CHG, ATOMIC_MASS_UNIT, SPEED_OF_LIGHT
    implicit none
    !> declare input variables
    type(particle_gc), intent(in)          :: gc_in
    real(kind=8), intent(in)               :: mass
    real(kind=8), dimension(3), intent(in) :: B
    !> declare output variables
    type(particle_gc_relativistic)         :: relativistic_gc
    
    relativistic_gc%x = gc_in%x
    relativistic_gc%q = gc_in%q
    relativistic_gc%i_elm = gc_in%i_elm
    relativistic_gc%st = gc_in%st
    relativistic_gc%p(2) = abs(EL_CHG*gc_in%mu/ATOMIC_MASS_UNIT)
    relativistic_gc%p(1) = sign(sqrt(((gc_in%E*gc_in%E*EL_CHG*EL_CHG)/         &
         (ATOMIC_MASS_UNIT*ATOMIC_MASS_UNIT*SPEED_OF_LIGHT*SPEED_OF_LIGHT)) -  &
         mass*mass*SPEED_OF_LIGHT*SPEED_OF_LIGHT -                             &
         2.d0*mass*norm2(B)*relativistic_gc%p(2)),gc_in%mu)
  end function gc_to_relativistic_gc

  !> This function fills in the p(1) and p(2) fields of a particle_gc_relativistic
  !> from its energy and (cosine of) pitch-angle
  function relativistic_gc_momenta_from_E_cospitch(rel_gc_in,energy,ksi,mass,fields,time) result(rel_gc_out)
    use constants, only: ATOMIC_MASS_UNIT, EL_CHG, SPEED_OF_LIGHT
    use mod_fields

    implicit none

    !> input variables
    type(particle_gc_relativistic), intent(in) :: rel_gc_in
    real*8, intent(in)                         :: energy !< Particle energy in eV (including rest energy)
    real*8, intent(in)                         :: ksi    !< Cosine of particle pitch-angle 
    real*8, intent(in)                         :: mass   !< Particle mass in AMU
    class(fields_base), intent(in)             :: fields
    real*8, intent(in)                         :: time

    !> output variables
    type(particle_gc_relativistic) :: rel_gc_out

    !> internal variables
    real*8               :: p_norm_sq, energy_at_rest
    real*8               :: psi, U, B_norm
    real*8, dimension(3) :: E, B

    !> Copy existing fields from input to output
    rel_gc_out%q     = rel_gc_in%q
    rel_gc_out%x     = rel_gc_in%x
    rel_gc_out%i_elm = rel_gc_in%i_elm
    rel_gc_out%st    = rel_gc_in%st
    
    if (abs(ksi)>1) then
      write(*,*) 'Error in relativistic_gc_momenta_from_E_cospitch: &
        abs(cos(pitch-angle)) should be <1, exiting'
      stop
    end if
    
    energy_at_rest = ATOMIC_MASS_UNIT*mass*SPEED_OF_LIGHT**2
    
    p_norm_sq = ((energy*EL_CHG)**2-energy_at_rest**2)/SPEED_OF_LIGHT**2

    !> mu [AMU (m/s)**2 /T]
    call fields%calc_EBpsiU(time, rel_gc_out%i_elm, rel_gc_out%st, rel_gc_out%x(3), E, B, psi, U)
    B_norm = norm2(B)
    rel_gc_out%p(2) = p_norm_sq*(1.d0-ksi**2)/(2.d0*B_norm*mass*ATOMIC_MASS_UNIT**2)

    !> Parallel momentum [AMU m/s]
    rel_gc_out%p(1) = sqrt(p_norm_sq)*ksi/ATOMIC_MASS_UNIT		       
  end function relativistic_gc_momenta_from_E_cospitch

  !> Compute relativistic factor from guiding center momentum
  function compute_relativistic_factor(gc, mass, Bnorm) result(gamma)
    type(particle_gc_relativistic), intent(in) :: gc
    real*8, intent(in) :: mass, Bnorm !< mass [amu], B-field magnitude [T]
    real*8 :: gamma

    gamma = sqrt( 1.d0 + ( gc%p(1) / ( mass * SPEED_OF_LIGHT ) )**2 &
                       + 2.d0 * Bnorm * gc%p(2) / ( mass * SPEED_OF_LIGHT**2 ) )
  end function compute_relativistic_factor


end module mod_gc_relativistic

