!> Mod_particle_type contains the default particle type, type_particle
!> If you wish to add a new particle type, create one here, inheriting from
!> [[particle_base]], and update [[mod_particle_io]] (search for `particle_kinetic`
!> and add your particle at each spot)
module mod_particle_types
  implicit none
  private
  public particle_base, particle_kinetic, particle_kinetic_leapfrog, particle_gc, particle_fieldline
  public particle_kinetic_relativistic, particle_gc_relativistic
  public particle_get_q
  public copy_particle
  public copy_particle_base
  public copy_particle_kinetic_leapfrog


  !> The base type for all other particles. Includes only the position and weight elements
  !> Integration in a 2D finite element method is included in the form of 2 coordinates
  !> and an element index.
  type, abstract :: particle_base
    real*8    :: x(3)             !< particle position in real space
    real*8    :: st(2)            !< particle position in the element
    real*4    :: weight = 1.0     !< weight (i.e. number of particles)
    integer*4 :: i_elm = 0        !< index in element_list. Negative indices indicate lost particles on the edge of - that element.
    integer*4 :: i_life = 0       !< particle lifetime index (i.e. is this still the same particle?)
    real*4    :: t_birth = 0.0    !< birth time of this particle
    !< zero means lost without location specification.
  contains
    procedure :: copy => copy_particle
    generic :: assignment(=) => copy
  end type particle_base

  !> A simple type just for fieldline tracing in two-step methods (Adams Bashforth) or for forward euler
  type, extends(particle_base) :: particle_fieldline
    real*8    :: B_hat_prev(3) = 0.d0 !< Field direction at previous timestep
    real*8    :: v = 0.d0 !< Parallel velocity along the fieldline
  end type particle_fieldline

  !> A simple guiding-center particle type.
  type, extends(particle_base) :: particle_gc
    real*8    :: E = 0.d0 !< The particle energy [eV]
    real*8    :: mu = 0.d0 !< The magnetic moment [eV/T]. Sign determines sign of v_par
    integer*1 :: q = 0_1 !< Charge [e]
  end type particle_gc

  !> For most kinetic methods the velocity is required at time \(t\)
  type, extends(particle_base) :: particle_kinetic
    real*8, dimension(3) :: v = 0.d0 !< Velocity [m/s]
    integer*1            :: q = 0_1 !< charge [e]
  end type particle_kinetic

  !> Leapfrog methods define the particle velocity at time \(t^{n-1/2}\)
  !> and are therefore incompatible with normal kinetic methods (but a conversion
  !> function should not be too difficult)
  type, extends(particle_base) :: particle_kinetic_leapfrog
    real*8, dimension(3) :: v = 0.d0 !< Velocity [m/s] at t=t^(n-1/2) (where the position is known at t^n)
    integer*1            :: q = 0_1 !< charge [e]
  end type particle_kinetic_leapfrog

  !> This particle type is used for computing the full orbit trajectory
  !> of a relativistic particle. Particle position and momentum are given at time \(t\)
 type, extends(particle_base) :: particle_kinetic_relativistic
    real(kind=8),dimension(3) :: p !< Momentum in Cartesian coordinates (p_x,p_y,p_z) in [AMU*m/s]
    integer(kind=1)           :: q !< charge [e]
 end type particle_kinetic_relativistic

 !> This particle type is used for computing the guiding center trajectory of a
 !> relativistic particle. GC position and momentum are given at the time (\t\)
 type, extends(particle_base) :: particle_gc_relativistic
    real(kind=8), dimension(2) :: p  !< 1: parallel momentum [AMU m/s], 2: magnetic moment [(AMU*m**2)/(T*s**2)]
    integer(kind=1) :: q !< charge [e]
 end type particle_gc_relativistic
contains
  !> Convenience function to obtain q if it exists, or 0 otherwise
  !> Here also because of https://gcc.gnu.org/bugzilla/show_bug.cgi?id=82064
  !> which means that we cannot use the same derived type in too many modules which will be
  !> imported in the main program (roughly)
  pure function particle_get_q(in) result(q)
    class(particle_base), intent(in) :: in
    integer*1 :: q
    select type (p => in)
    type is (particle_kinetic)
      q = p%q
    type is (particle_kinetic_leapfrog)
      q = p%q
    type is (particle_gc)
      q = p%q
    type is (particle_kinetic_relativistic)
      q = p%q
    type is (particle_gc_relativistic)
      q=p%q
    class default
      q = 0
    end select
  end function particle_get_q
  !> Copy the base variables from one particle of class(particle_base) to another
  
  pure subroutine copy_particle_base(in, out)
    class(particle_base), intent(in)    :: in
    class(particle_base), intent(inout) :: out
    out%x      = in%x
    out%st     = in%st
    out%weight = in%weight
    out%i_elm  = in%i_elm
    out%i_life = in%i_life
    out%t_birth= in%t_birth
  end subroutine copy_particle_base

  !> Copy one particle of a type kinetic_leapfrog to another
  pure subroutine copy_particle_kinetic_leapfrog(in, out)
    type(particle_kinetic_leapfrog), intent(in)    :: in
    type(particle_kinetic_leapfrog), intent(inout) :: out
    out%x       = in%x
    out%st      = in%st
    out%weight  = in%weight
    out%i_elm   = in%i_elm
    out%i_life  = in%i_life
    out%t_birth = in%t_birth
    out%v       = in%v
    out%q       = in%q
  end subroutine copy_particle_kinetic_leapfrog

  !> Copy a descendant of particle_base into another descendant of particle_base
  !> as requested by the types of the input and output parameters.
  !> We do not transform any of the quantities between eachother, just copy them
  !> if the particles are of the same type.
  !> To do a proper transform we typically need additional parameters, like the
  !> mass or magnetic field. See [[mod_boris.f90]] for some examples.
  !> See https://stackoverflow.com/a/19082934
  subroutine copy_particle(particle_out, particle_in)
    class(particle_base), intent(out) :: particle_out !< Particle to copy attributes into
    class(particle_base), intent(in)  :: particle_in  !< Particle to copy attributes from

    particle_out%x      = particle_in%x
    particle_out%st     = particle_in%st
    particle_out%weight = particle_in%weight
    particle_out%i_elm  = particle_in%i_elm

    select type (p_out => particle_out)
    type is (particle_fieldline)
      select type (p_in => particle_in)
      type is (particle_fieldline) ! this is a straight copy, simple
        p_out%B_hat_prev = p_in%B_hat_prev
        p_out%v = p_in%v
      class default
        ! Maybe we should warn here instead of just putting zeros,
        ! or put nothing so the compiler can catch the uninitialized value
        p_out%B_hat_prev = [0.d0, 0.d0, 0.d0]
        p_out%v = 0.d0
      end select
    type is (particle_gc)
      select type (p_in => particle_in)
      type is (particle_gc)
        p_out%E  = p_in%E
        p_out%mu = p_in%mu
        p_out%q  = p_in%q
      class default
        p_out%E  = 0.d0
        p_out%mu = 0.d0
        p_out%q  = 0
      end select
    type is (particle_kinetic)
      select type (p_in => particle_in)
      type is (particle_kinetic)
        p_out%v  = p_in%v
        p_out%q  = p_in%q
      class default
        p_out%v  = [0.d0, 0.d0, 0.d0]
        p_out%q  = 0
      end select
    type is (particle_kinetic_leapfrog)
      select type (p_in => particle_in)
      type is (particle_kinetic_leapfrog)
        p_out%v  = p_in%v
        p_out%q  = p_in%q
      class default
        ! the transformation from kinetic to kinetic_leapfrog could be done with a small error here
        p_out%v  = [0.d0, 0.d0, 0.d0]
        p_out%q  = 0
      end select
    type is (particle_kinetic_relativistic)
      select type (p_in => particle_in)
      type is (particle_kinetic_relativistic)
        p_out%p  = p_in%p
        p_out%q  = p_in%q
      class default
        p_out%p = [0.d0,0.d0,0.d0]
        p_out%q = 0
      end select     
     type is (particle_gc_relativistic)
       select type (p_in => particle_in)
       type is (particle_gc_relativistic)
         p_out%p = p_in%p
         p_out%q = p_in%q
     class default
        p_out%p  = [0.d0, 0.d0]
        p_out%q  = 0
      end select
    end select
  end subroutine copy_particle

end module mod_particle_types
