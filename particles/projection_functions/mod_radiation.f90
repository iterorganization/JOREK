!> module takes the OPEN-ADAS data to calculate the radiated power of a single atom
module mod_radiation
use mod_parameters
use mod_openadas
use data_structure
use mod_particle_sim
use mod_particle_types, only: particle_get_q, particle_base, particle_marker
use mod_ionisation_recombination, only: fields_interp_ne_Te
implicit none
private
public proj_Lz, proj_Lz_equil, get_Lz

contains

!> Project the particle radiated power density L_z (W/atom)
!> This function is not pure due to the SPLINE call in %interp in adas
!> which is a bit risky
function proj_Lz(sim, group, particle)
  type(particle_sim), intent(in) :: sim
  integer, intent(in) :: group
  class(particle_base), intent(in) :: particle
  real*8 :: proj_Lz
  real*8 :: n_e, T_e, log_T_e, log_n_e
  real*8 :: prb, plt, prc
  integer :: q
#ifdef WITH_Neutrals
  real*8 :: n_n
  real*8, dimension(1) :: P, P_s, P_t, P_phi, P_time
  real*8 :: R, R_s, R_t, Z, Z_s, Z_t
#endif

  ! Calculate local temperature, density
  call fields_interp_ne_Te(sim%fields, sim%time, particle%st(1), particle%st(2), &
      particle%x(3), particle%i_elm, n_e, T_e)
  log_T_e = log10(T_e)
  log_n_e = log10(n_e)

#ifdef WITH_Neutrals
  ! Calculate neutral_density if model5XX (model501 has n_imp in 8)
  call sim%fields%interp_PRZ(sim%time,particle%i_elm,[var_rhon],1,particle%st(1), &
      particle%st(2),particle%x(3),P,P_s,P_t,P_phi,P_time,R,R_s,R_t,Z,Z_s,Z_t)
  n_n = P(1)
#endif

  q = particle_get_q(particle)
  ! From here on out we have a q
  call sim%groups(group)%ad%PRB%interp_linear(q, log_n_e, log_T_e, prb)
  call sim%groups(group)%ad%PLT%interp_linear(q, log_n_e, log_T_e, plt)
  proj_Lz      = (prb + plt) * n_e
#ifdef WITH_Neutrals
  call sim%groups(group)%ad%PRC%interp_linear(q, log_n_e, log_T_e, prc)
  proj_Lz      = proj_Lz + prc * n_n
#endif
end function proj_Lz
!function proj_Lz(sim, group, particle)
!  type(particle_sim), intent(in) :: sim
!  integer, intent(in) :: group
!  class(particle_base), intent(in) :: particle
!  real*8 :: proj_Lz
!  real*8 :: n_e, T_e, log_T_e, log_n_e
!  real*8 :: prb, plt, prc
!  real*8, dimension(0:sim%groups(group)%ad%n_Z) :: rad, rad_RC
!  integer :: q, iZ
!#if (JOREK_MODEL == 500) || (JOREK_MODEL == 555)
!  real*8 :: n_n
!  real*8, dimension(1) :: P, P_s, P_t, P_phi, P_time
!  real*8 :: R, R_s, R_t, Z, Z_s, Z_t
!#endif
!
!  ! Calculate local temperature, density
!  call fields_interp_ne_Te(sim%fields, sim%time, particle%st(1), particle%st(2), &
!      particle%x(3), particle%i_elm, n_e, T_e)
!  log_T_e = log(T_e)
!  log_n_e = log(n_e)
!
!#if (JOREK_MODEL == 500) || (JOREK_MODEL == 555)
!  ! Calculate neutral_density if model5XX (model501 has n_imp in 8)
!  call sim%fields%interp_PRZ(sim%time,particle%i_elm,[var_rhon],1,particle%st(1), &
!      particle%st(2),particle%x(3),P,P_s,P_t,P_phi,P_time,R,R_s,R_t,Z,Z_s,Z_t)
!  n_n = P(1)
!#endif
!
!  select type(p=>particle)
!  type is (particle_marker)
!    do iZ=0,sim%groups(group)%ad%n_Z
!      call sim%groups(group)%ad%PRB%interp_linear(iZ, log_n_e, log_T_e, prb)
!      call sim%groups(group)%ad%PLT%interp_linear(iZ, log_n_e, log_T_e, plt)
!      call sim%groups(group)%ad%PRC%interp_linear(iZ, log_n_e, log_T_e, prc)
!      rad(iZ) = prb + plt
!    enddo
!    proj_Lz      = dot_product(p%P_imp(0:sim%groups(group)%ad%n_Z), rad*n_e)
!#if (JOREK_MODEL == 500) || (JOREK_MODEL == 555)
!    proj_Lz      = proj_Lz + dot_product(p%P_imp(0:sim%groups(group)%ad%n_Z), prc*n_n)
!#endif
!  class default
!    q = particle_get_q(particle)
!    ! From here on out we have a q
!    call sim%groups(group)%ad%PRB%interp_linear(q, log_n_e, log_T_e, prb)
!    call sim%groups(group)%ad%PLT%interp_linear(q, log_n_e, log_T_e, plt)
!    proj_Lz      = (prb + plt) * n_e
!#if (JOREK_MODEL == 500) || (JOREK_MODEL == 555)
!    call sim%groups(group)%ad%PRC%interp_linear(q, log_n_e, log_T_e, prc)
!    proj_Lz      = proj_Lz + prc * n_n
!#endif
!  end select
!end function proj_Lz

function get_Lz(sim, group, iz, n_e, T_e, n_n)
  type(particle_sim), intent(in) :: sim
  integer, intent(in) :: group
  integer, intent(in) :: iz
  real*8, intent(in)  :: n_e, T_e
  real*8, intent(in), optional :: n_n
  real*8 :: get_Lz
  real*8 :: log_T_e, log_n_e
  real*8 :: prb, plt, prc
  integer :: q

  log_T_e = log10(T_e)
  log_n_e = log10(n_e)

  ! From here on out we have a q
  call sim%groups(group)%ad%PRB%interp_linear(iz, log_n_e, log_T_e, prb)
  call sim%groups(group)%ad%PLT%interp_linear(iz, log_n_e, log_T_e, plt)
  get_Lz      = (prb + plt) * n_e
#ifdef WITH_Neutrals
  call sim%groups(group)%ad%PRC%interp_linear(iz, log_n_e, log_T_e, prc)
  get_Lz      = get_Lz + prc * n_n
#endif
end function get_Lz

function get_PRB(sim, group, iz, n_e, T_e, n_n)
  type(particle_sim), intent(in) :: sim
  integer, intent(in) :: group
  integer, intent(in) :: iz
  real*8, intent(in)  :: n_e, T_e
  real*8, intent(in), optional :: n_n
  real*8 :: get_PRB
  real*8 :: log_T_e, log_n_e
  real*8 :: prb
  integer :: q

  log_T_e = log10(T_e)
  log_n_e = log10(n_e)

  ! From here on out we have a q
  call sim%groups(group)%ad%PRB%interp_linear(iz, log_n_e, log_T_e, prb)
  get_PRB     = prb * n_e
end function get_PRB

function get_PLT(sim, group, iz, n_e, T_e, n_n)
  type(particle_sim), intent(in) :: sim
  integer, intent(in) :: group
  integer, intent(in) :: iz
  real*8, intent(in)  :: n_e, T_e
  real*8, intent(in), optional :: n_n
  real*8 :: get_PRB
  real*8 :: log_T_e, log_n_e
  real*8 :: plt
  integer :: q

  log_T_e = log10(T_e)
  log_n_e = log10(n_e)

  ! From here on out we have a q
  call sim%groups(group)%ad%PLT%interp_linear(iz, log_n_e, log_T_e, plt)
  get_PLT     = plt * n_e
end function get_PLT

!> Project the particle radiated power density L_z (W/atom) in the equilibrium
!> calculation without neutrals
function proj_Lz_equil(sim, group, particle) result(P_rad)
  use mod_coronal
  type(particle_sim), intent(in) :: sim
  integer, intent(in) :: group
  class(particle_base), intent(in) :: particle
  real*8 :: P_rad
  real*8 :: n_e, T_e, log_T_e, log_n_e
  real*8, allocatable :: fractions(:)

  ! Calculate local temperature, density
  call sim%fields%calc_NeTe(sim%time, particle%i_elm, particle%st, particle%x(3), n_e, T_e)
  log_T_e = log10(T_e)
  log_n_e = log10(n_e)

  !call sim%groups(group)%cor%interp_linear(log_n_e, log_T_e, rad=P_rad)

  ! This is more expensive but maybe correct
  allocate(fractions(0:sim%groups(group)%ad%n_Z))
  fractions = specific_coronal_equilibrium(sim%groups(group)%ad, log_n_e, log_T_e, .true.)

  P_rad = coronal_Prad(sim%groups(group)%ad, log_n_e, log_T_e, fractions)
end function proj_Lz_equil

!subroutine Lz_and_ionization(sim, group, particle, Lz, ion, rcb)
!  type(particle_sim), intent(in) :: sim
!  integer, intent(in) :: group
!  real*8, intent(out) :: Lz
!  real*8, dimension(0:ad%n_Z), intent(out), optional :: ion, rcb
!  class(particle_base), intent(in) :: particle
!  real*8 :: n_e, T_e, log_T_e, log_n_e
!  real*8 :: prb, plt, prc
!  real*8, dimension(0:ad%n_Z) :: rad, rad_RC
!  integer :: q, iZ
!#if (JOREK_MODEL == 500) || (JOREK_MODEL == 555)
!  real*8 :: n_n
!  real*8, dimension(1) :: P, P_s, P_t, P_phi, P_time
!  real*8 :: R, R_s, R_t, Z, Z_s, Z_t
!#endif
!
!  ! Calculate local temperature, density
!  call fields_interp_ne_Te(sim%fields, sim%time, particle%st(1), particle%st(2), &
!      particle%x(3), particle%i_elm, n_e, T_e)
!  log_T_e = log(T_e)
!  log_n_e = log(n_e)
!
!#if (JOREK_MODEL == 500) || (JOREK_MODEL == 555)
!  ! Calculate neutral_density if model5XX (model501 has n_imp in 8)
!  call sim%fields%interp_PRZ(sim%time,particle%i_elm,[var_rhon],1,particle%st(1), &
!      particle%st(2),particle%x(3),P,P_s,P_t,P_phi,P_time,R,R_s,R_t,Z,Z_s,Z_t)
!  n_n = P(1)
!#endif
!
!  select type(p=>particle)
!  type is (particle_marker)
!    do iZ=0,sim%groups(group)%ad%n_Z
!      call sim%groups(group)%ad%PRB%interp_linear(iZ, log_n_e, log_T_e, prb)
!      call sim%groups(group)%ad%PLT%interp_linear(iZ, log_n_e, log_T_e, plt)
!      call sim%groups(group)%ad%PRC%interp_linear(iZ, log_n_e, log_T_e, prc)
!      if (present(ion)) call sim%groups(group)%ad%SCD%interp_linear(iZ, log_n_e, log_T_e, ion(iZ))
!      if (present(rcb)) call sim%groups(group)%ad%ACD%interp_linear(iZ, log_n_e, log_T_e, rcb(iZ))
!      
!      rad(iZ) = prb + plt
!    enddo
!    Lz      = dot_product(p%P_imp(0:sim%groups(group)%ad%n_Z), rad*n_e)
!#if (JOREK_MODEL == 500) || (JOREK_MODEL == 555)
!    Lz      = Lz + dot_product(p%P_imp(0:sim%groups(group)%ad%n_Z), prc*n_n)
!#endif
!  class default
!    q = particle_get_q(particle)
!    ! From here on out we have a q
!    call sim%groups(group)%ad%PRB%interp_linear(q, log_n_e, log_T_e, prb)
!    call sim%groups(group)%ad%PLT%interp_linear(q, log_n_e, log_T_e, plt)
!    Lz      = (prb + plt) * n_e
!#if (JOREK_MODEL == 500) || (JOREK_MODEL == 555)
!    call sim%groups(group)%ad%PRC%interp_linear(q, log_n_e, log_T_e, prc)
!    Lz      = Lz + prc * n_n
!#endif
!  end select
!end subroutine Lz_and_ionization
end module mod_radiation
