!> module takes the OPEN-ADAS data to calculate the radiated power of a single atom
module mod_radiation
use mod_parameters
use mod_openadas
use data_structure
use mod_particle_sim
use mod_particle_types, only: particle_get_q, particle_base
use mod_ionisation_recombination, only: fields_interp_ne_Te
implicit none
private
public proj_Lz

contains

!> Project the particle radiated power density L_z (W/atom)
function proj_Lz(sim, group, particle)
  type(particle_sim), intent(in) :: sim
  integer, intent(in) :: group
  class(particle_base), intent(in) :: particle
  real*8 :: proj_Lz
  real*8 :: n_e, T_e, log_T_e, log_n_e
  real*8 :: prb, plt, prc
  real*8, dimension(0:ad%n_Z) :: rad, rad_RC
  integer :: q, iZ
#if (JOREK_MODEL == 500) || (JOREK_MODEL == 555)
  real*8 :: n_n
  real*8, dimension(1) :: P, P_s, P_t, P_phi, P_time
  real*8 :: R, R_s, R_t, Z, Z_s, Z_t
#endif

  ! Calculate local temperature, density
  call fields_interp_ne_Te(sim%fields, sim%time, particle%st(1), particle%st(2), &
      particle%x(3), particle%i_elm, n_e, T_e)
  log_T_e = log(T_e)
  log_n_e = log(n_e)

#if (JOREK_MODEL == 500) || (JOREK_MODEL == 555)
  ! Calculate neutral_density if model5XX (model501 has n_imp in 8)
  call sim%fields%interp_PRZ(sim%time,particle%i_elm,[var_rhon],1,particle%st(1), &
      particle%st(2),particle%x(3),P,P_s,P_t,P_phi,P_time,R,R_s,R_t,Z,Z_s,Z_t)
  n_n = P(1)
#endif

  select type(p=>particle)
  type is (particle_marker)
    do iZ=0,sim%groups(group)%ad%n_Z
      call sim%groups(group)%ad%PRB%interp_linear(iZ, log_n_e, log_T_e, prb)
      call sim%groups(group)%ad%PLT%interp_linear(iZ, log_n_e, log_T_e, plt)
      call sim%groups(group)%ad%PRC%interp_linear(iZ, log_n_e, log_T_e, prc)
      rad(iZ) = prb + plt
    enddo
    proj_Lz      = dot_product(p%P_imp(0:sim%groups(group)%ad%n_Z), rad*n_e)
#if (JOREK_MODEL == 500) || (JOREK_MODEL == 555)
    proj_Lz      = proj_Lz + dot_product(p%P_imp(0:sim%groups(group)%ad%n_Z), prc*n_n)
#endif
  class default
    q = particle_get_q(particle)
    ! From here on out we have a q
    call sim%groups(group)%ad%PRB%interp_linear(q, log_n_e, log_T_e, prb)
    call sim%groups(group)%ad%PLT%interp_linear(q, log_n_e, log_T_e, plt)
    proj_Lz      = (prb + plt) * n_e
#if (JOREK_MODEL == 500) || (JOREK_MODEL == 555)
    call sim%groups(group)%ad%PRC%interp_linear(q, log_n_e, log_T_e, prc)
    proj_Lz      = proj_Lz + prc * n_n
#endif
  end select
end function proj_Lz
end module mod_radiation
