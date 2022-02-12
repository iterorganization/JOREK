!> Unit test for verifying radiation reaction force operator has not changed
module radreactforce
  use constants, only: EL_CHG, SPEED_OF_LIGHT, ATOMIC_MASS_UNIT, MASS_ELECTRON
  use mod_particle_types
  use mod_radreactforce
  use fruit

  implicit none

contains

  subroutine test_radreactforce
    implicit none

    type(particle_kinetic_relativistic) :: prt
    type(particle_gc_relativistic) :: gc

    real(kind=8) :: mass = MASS_ELECTRON / ATOMIC_MASS_UNIT
    real(kind=8) :: dt   = 1.0e-3, energy, pitch, predefval
    real(kind=8) :: B(3), p0_prt(3), ptemp_prt(3), p0_gc, p_gc, ptemp_gc, yout(2)

    prt%q = -1.0
    gc%q  = -1.0

    prt%x = (/1.d0, 0.d0, 0.d0/)
    gc%x  = (/1.d0, 0.d0, 0.d0/)

    ! Radiation reaction force is zero when p_par == p_tot and p_perp == 0. Check!
    B = (/0.d0, 0.d0, 1.d0/)
    energy = 1e6
    pitch  = 1.d0
    call momentum_prt_and_gc(energy, pitch, mass, norm2(B), prt, gc)
    p0_prt = prt%p
    p0_gc  = gc2momentum(gc, mass, norm2(B))

    call radreactforce_kinetic(B, dt, mass, prt)
    call radreactforce_gc(B, dt, mass, gc)
    p_gc = gc2momentum(gc, mass, norm2(B))

    call assert_equals(0.d0, norm2(prt%p - p0_prt) / norm2(p0_prt), 1.d-6, "p - p0 = 0 (prt)")
    call assert_equals(0.d0, abs(p_gc - p0_gc) / p0_gc, 1d-6, "p - p0 = 0 (gc)")

    ! Check B^-2 scaling
    B = (/0.d0, 0.d0, 1.d0/)
    energy = 1e6
    pitch  = 0.5d0
    call momentum_prt_and_gc(energy, pitch, mass, norm2(B), prt, gc)
    p0_prt = prt%p
    p0_gc  = gc%p(1)

    call radreactforce_kinetic(B, dt, mass, prt)
    call radreactforce_gc(B, dt, mass, gc)
    ptemp_prt = prt%p
    ptemp_gc  = gc%p(1)

    B = (/0.d0, 0.d0, 1.d0/)*4 ! Scale by a factor of 4
    call momentum_prt_and_gc(energy, pitch, mass, norm2(B), prt, gc)
    call radreactforce_kinetic(B, dt, mass, prt)
    call radreactforce_gc(B, dt, mass, gc)
    p_gc  = gc%p(1)

    call assert_equals( 4.d0**2, norm2(prt%p - p0_prt) / norm2(ptemp_prt - p0_prt), 1.d-6, "Delta p propto B^2")
    call assert_equals( 4.d0**2, (p_gc - p0_gc) / (ptemp_gc - p0_gc), 1.d-6, "Delta p_par propto B^2")

    ! Running a proper physics test (bump-on-tail) is too expensive so we only
    ! check that the operator has not changed since the physics tests were done.
    ! We do this by comparing the result to a precomputed value (evaluated
    ! when the operator was confirmed to pass the physics test). Also
    ! check that prt == gc (when the field is uniform).

    B = (/0.d0, 0.d0, 1.d0/)
    energy = 1e6
    pitch  = 0.5d0
    call momentum_prt_and_gc(energy, pitch, mass, norm2(B), prt, gc)
    call radreactforce_kinetic(B, dt, mass, prt)
    call radreactforce_gc(B, dt, mass, gc)
    p0_prt = prt%p
    p0_gc  = gc2momentum(gc, mass, norm2(B))
    
    ! Check RHS routine as well
    call momentum_prt_and_gc(energy, pitch, mass, norm2(B), prt, gc)
    call radreactforce_gc_rhs(B, mass, int(gc%q), gc%p, yout)
    gc%p = gc%p + dt * yout
    p_gc = gc2momentum(gc, mass, norm2(B))

    predefval = 457450.d0
    call assert_equals(0.d0, (p0_gc - norm2(p0_prt)) / p0_gc, 1.d-6, "GC == PRT")
    call assert_equals(0.d0, (p0_gc - p_gc) / p0_gc, 1.d-6, "GC == GC_RHS")
    call assert_equals(0.d0, (p0_gc - predefval) / p0_gc, 1.d-5, "Regression test")
    

  end subroutine test_radreactforce

  !> Simple function to initialize particle and gc momentum coordiantes from energy and pitch
  subroutine momentum_prt_and_gc(energy, pitch, mass, Bnorm, prt, gc)
    implicit none

    real*8, intent(in) :: energy, pitch, mass, Bnorm !< Ekin [eV], ppar/ptot, mass [amu], B-field magnitude [T]
    type(particle_kinetic_relativistic), intent(inout) :: prt
    type(particle_gc_relativistic), intent(inout) :: gc

    real*8 :: p

    p = sqrt( (energy * EL_CHG / ( mass * SPEED_OF_LIGHT**2 * ATOMIC_MASS_UNIT ) + 1.d0)**2 - 1.d0  ) * mass * SPEED_OF_LIGHT

    prt%p = (/sqrt( 1.d0 - pitch**2 ) * p, pitch * p, 0.d0/)

    gc%p =(/ pitch * p, ( 1.d0 - pitch**2 ) * p**2 / ( 2.d0 * mass * Bnorm )/)

  end subroutine momentum_prt_and_gc

  !> Evaluate momentum (norm) from gc momentum coordinates
  function gc2momentum(gc, mass, Bnorm) result(p)
    type(particle_gc_relativistic), intent(in) :: gc
    real*8, intent(in) :: mass, Bnorm !< mass [amu], B-field magnitude [T]
    real*8 :: p

    p = sqrt( gc%p(1)**2 + 2.d0 * mass * Bnorm * gc%p(2) )
  end function gc2momentum

end module radreactforce
