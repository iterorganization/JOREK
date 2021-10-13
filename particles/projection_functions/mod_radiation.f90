!> module takes the OPEN-ADAS data to calculate the radiated power of a single atom
module mod_radiation
use mod_openadas
use data_structure
use mod_particle_sim
use mod_particle_types, only: particle_get_q, particle_base
implicit none
private
public proj_Lz, proj_Lz_equil

contains

!> Project the particle radiated power density L_z (W/atom)
!> This function is not pure due to the SPLINE call in %interp in adas
!> which is a bit risky
function proj_Lz(sim, group, particle)
  type(particle_sim), intent(in) :: sim
  integer, intent(in) :: group
  class(particle_base), intent(in) :: particle
  real*8 :: proj_Lz, PRB, PLT, PRC
  real*8 :: n_e, T_e, log_T_e, log_n_e
  integer :: q
#if (JOREK_MODEL == 500) || (JOREK_MODEL == 555)
  real*8 :: n_n
  real*8, dimension(1) :: P, P_s, P_t, P_phi, P_time
  real*8 :: R, R_s, R_t, Z, Z_s, Z_t
#endif

  ! Calculate local temperature, density
  call sim%fields%calc_NeTe(sim%time, particle%i_elm, particle%st, particle%x(3), n_e, T_e)
  log_T_e = log10(T_e)
  log_n_e = log10(n_e)

#if (JOREK_MODEL == 500) || (JOREK_MODEL == 555)
  ! Calculate neutral_density if model5XX (model501 has n_imp in 8)
  call sim%fields%interp_PRZ(sim%time,particle%i_elm,[8],1,particle%st(1), &
      particle%st(2),particle%x(3),P,P_s,P_t,P_phi,P_time,R,R_s,R_t,Z,Z_s,Z_t)
  n_n = P(1)
#endif

  q = particle_get_q(particle)
  ! From here on out we have a q
  PRB = sim%groups(group)%ad%PRB%interp_linear(q, log_n_e, log_T_e)
  PLT = sim%groups(group)%ad%PLT%interp_linear(q, log_n_e, log_T_e)
  proj_Lz      = (PRB + PLT) * n_e
#if (JOREK_MODEL == 500) || (JOREK_MODEL == 555)
  PRC = sim%groups(group)%ad%PRC%interp_linear(q, log_n_e, log_T_e)
  proj_Lz      = proj_Lz + PRC * n_n
#endif
end function proj_Lz



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
  fractions = specific_coronal_equilibrium(sim%groups(group)%ad, log_n_e, log_T_e)

  P_rad = coronal_Prad(sim%groups(group)%ad, log_n_e, log_T_e, fractions)
end function proj_Lz_equil
end module mod_radiation
