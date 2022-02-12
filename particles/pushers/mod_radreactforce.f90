!< Radiation reaction force a.k.a. synchrotron loss operator
!!
!! This module contains operators for relativistic particle and guiding center
!! to take into account the radiation reaction force. This is only relevant
!! for relativistic electrons and, hence, this module assumes the test particle is
!! an electron.
!!
!! The implemented operators are taken from here:
!! https://arxiv.org/pdf/1412.1966.pdf
module mod_radreactforce
  use constants, only: EL_CHG, SPEED_OF_LIGHT, ATOMIC_MASS_UNIT, PI, EPS_ZERO
  use mod_particle_types

  implicit none
  
contains

  !< Calculate the characteristic time for the radiation reaction force
  function radreactforce_chartime(B, gamma, mass, charge)
    implicit none
    real*8 :: radreactforce_chartime
    real*8, intent(in) :: B     !< Magnetic field strength [T]
    real*8, intent(in) :: gamma !< Particle Lorentz factor
    real*8, intent(in) :: mass  !< Particle mass [kg]
    integer, intent(in) :: charge !< Particle charge [e]

    radreactforce_chartime = 6 * PI * EPS_ZERO * gamma * ( mass * SPEED_OF_LIGHT )**3 &
                           / ( (charge * EL_CHG)**4 * B**2 )

  end function radreactforce_chartime

  !< Apply radiation reaction force to a particle
  !!
  !! Radiation reaction force is applied with forward Euler method. The term
  !! containing explicit magnetic field time-dependence is omitted for now
  subroutine radreactforce_kinetic(B, dt, mass, particle)

    use mod_coordinate_transforms, only: vector_cylindrical_to_cartesian

    implicit none
    
    real(kind=8), intent(in)       :: B(3)   !< Magnetic field vector
    real(kind=8), intent(in)       :: dt     !< Time step
    real(kind=8), intent(in)       :: mass   !< Particle mass [AMU]
    type(particle_kinetic_relativistic), intent(inout) :: particle !< Particle to be operated on

    real*8 :: pperp(3), p(3) !< Perpendicular and total momentum in SI units and in cartesian coordinates
    real*8 :: Bxyz(3), Bnorm !< Magnetic field vector in cartesian coordinates and magnetic field norm
    real*8 :: tau, gamma, m  !< Characteristic time, Lorentz factor, mass in SI units

    ! Convert to SI units
    m = mass * ATOMIC_MASS_UNIT
    p = particle%p * ATOMIC_MASS_UNIT

    Bxyz  = vector_cylindrical_to_cartesian(particle%x(3), B)
    Bnorm = norm2(B)

    ! Calculate the characteristic time and pperp
    gamma = sqrt(1.0 + ( norm2(p) / ( m * SPEED_OF_LIGHT ) )**2 )
    tau   = radreactforce_chartime(Bnorm, gamma, m, int(particle%q))
    pperp = p - Bxyz * dot_product(Bxyz,p) / Bnorm**2

    ! Apply RR-force and convert back to JOREK units
    p = p - dt * ( pperp + p * norm2(pperp)**2 / ( m * SPEED_OF_LIGHT )**2 ) / tau
    particle%p = p / ATOMIC_MASS_UNIT

  end subroutine radreactforce_kinetic

  !< Apply radiation reaction force to a guiding center
  !!
  !! Radiation reaction force is applied with forward Euler method. Only the zeroth
  !! order terms are included for now (magnetic field non-uniformity is neglegted)
  subroutine radreactforce_gc(B, dt, mass, particle)
    implicit none

    real(kind=8), intent(in)       :: B(3)   !< Magnetic field vector
    real(kind=8), intent(in)       :: dt     !< Time step
    real(kind=8), intent(in)       :: mass   !< Particle mass [AMU]
    type(particle_gc_relativistic), intent(inout) :: particle !< Particle to be operated

    real*8 :: tau, Bnorm, ppar, mu, gamma, m

    ! Convert to SI units
    m    = mass * ATOMIC_MASS_UNIT
    ppar = particle%p(1) * ATOMIC_MASS_UNIT
    mu   = particle%p(2) * ATOMIC_MASS_UNIT

    ! Calculate the characteristic time
    Bnorm = norm2(B)
    gamma = sqrt( 1.0 + ( ppar / ( m * SPEED_OF_LIGHT ) )**2 + 2 * Bnorm * m * mu / ( m * SPEED_OF_LIGHT )**2 )
    tau   = radreactforce_chartime(Bnorm, gamma, m, int(particle%q))

    ! Apply RR-force
    particle%p(1) = ppar - dt * 2 * ppar * mu * Bnorm / ( m * SPEED_OF_LIGHT**2 * tau )
    particle%p(2) = mu   - dt * 2 * mu * ( 1.0 + 2 * mu * Bnorm / ( m * SPEED_OF_LIGHT**2 ) ) / tau

    ! Convert back to JOREK units
    particle%p(1) = particle%p(1) / ATOMIC_MASS_UNIT
    particle%p(2) = particle%p(2) / ATOMIC_MASS_UNIT

  end subroutine radreactforce_gc

  !< RHS of the guiding center radiation reaction force to be used in RK scheme
  subroutine radreactforce_gc_rhs(B, mass, charge, p, yout)
    implicit none

    real(kind=8), intent(in)       :: B(3)    !< Magnetic field vector
    real(kind=8), intent(in)       :: mass    !< Particle mass [AMU]
    integer,      intent(in)       :: charge  !< Particle charge [e]
    real(kind=8), intent(in)       :: p(2)    !< Guiding center ppar and mu
    real(kind=8), intent(inout)    :: yout(2) !< [dppar / dt, dmu / dt]

    real*8 :: tau, Bnorm, ppar, mu, gamma, m

    ! Convert to SI units
    m    = mass * ATOMIC_MASS_UNIT
    ppar = p(1) * ATOMIC_MASS_UNIT
    mu   = p(2) * ATOMIC_MASS_UNIT

    ! Calculate the characteristic time
    Bnorm = norm2(B)
    gamma = sqrt( 1.0 + ( ppar / ( m * SPEED_OF_LIGHT ) )**2 + 2 * Bnorm * m * mu / ( m * SPEED_OF_LIGHT )**2 )
    tau   = radreactforce_chartime(Bnorm, gamma, m, charge)

    ! Apply RR-force
    yout(1) = - 2 * ppar * mu * Bnorm / ( m * SPEED_OF_LIGHT**2 * tau )
    yout(2) = - 2 * mu * ( 1.0 + 2 * mu * Bnorm / ( m * SPEED_OF_LIGHT**2 ) ) / tau

    ! Convert back to JOREK units
    yout = yout / ATOMIC_MASS_UNIT

  end subroutine radreactforce_gc_rhs

end module mod_radreactforce
