!> module contains functions to project certain quantities of a particle to a jorek field.
module mod_projection_functions
use mod_openadas
use data_structure
use mod_particle_sim
use mod_particle_types
use phys_module, only: CENTRAL_MASS, CENTRAL_DENSITY
use constants,   only: MU_ZERO, MASS_PROTON, ATOMIC_MASS_UNIT, K_BOLTZ

implicit none

private
public proj_f_ion_source,          proj_f_ionization_energy,    proj_f_ionization_par_momentum,  &
       proj_f_CX_source,           proj_f_CX_energy,            proj_f_CX_par_momentum,          &
       proj_f_combined_density_SI, proj_f_combined_energy_SI,   proj_f_combined_par_momentum_SI, &
       proj_f_combined_density,    proj_f_combined_energy,      proj_f_combined_par_momentum

contains

function proj_f_combined_density_SI(sim, i_group,  particle) result(combined_density_SI)
  type(particle_sim), intent(in)    :: sim
  class(particle_base), intent(in)  :: particle
  integer, intent(in)               :: i_group
  real*8                            :: combined_density_SI

  select type (pa => particle)
  type is (particle_kinetic_leapfrog)
    combined_density_SI = proj_f_ion_source(sim, i_group, particle)
  class default
    combined_density_SI = 0.d0
  end select
end function proj_f_combined_density_SI

function proj_f_combined_density(sim, i_group,  particle) result(combined_density)
  type(particle_sim), intent(in)    :: sim
  class(particle_base), intent(in)  :: particle
  integer, intent(in)               :: i_group
  real*8                            :: combined_density

  select type (pa => particle)
  type is (particle_kinetic_leapfrog)
    combined_density =  proj_f_combined_density_SI(sim, i_group, particle) & 
                      * 1.d0  *sqrt((MU_ZERO * CENTRAL_MASS * MASS_PROTON) / (CENTRAL_DENSITY * 1.d20))
  class default
    combined_density = 0.d0
  end select
end function proj_f_combined_density

function proj_f_combined_energy_SI(sim, i_group,  particle) result(combined_energy_SI)
  type(particle_sim), intent(in)    :: sim
  class(particle_base), intent(in)  :: particle
  integer, intent(in)               :: i_group
  real*8                            :: combined_energy_SI

  select type (pa => particle)
  type is (particle_kinetic_leapfrog)
    combined_energy_SI = proj_f_ionization_energy(sim, i_group, particle) &
                    + proj_f_CX_energy(sim, i_group, particle)
  class default
    combined_energy_SI = 0.d0
  end select
end function proj_f_combined_energy_SI
  
function proj_f_combined_energy(sim, i_group,  particle) result(combined_energy)
  type(particle_sim), intent(in)    :: sim
  class(particle_base), intent(in)  :: particle
  integer, intent(in)               :: i_group
  real*8                            :: combined_energy
  
  real*8    :: rho(1), R, Z
  real*8    :: rho_s(1), rho_t(1), rho_phi(1), rho_time(1), R_s, R_t, Z_s, Z_t

  select type (pa => particle)
  type is (particle_kinetic_leapfrog)
    
    call sim%fields%interp_PRZ(sim%time, particle%i_elm, [5], 1, particle%st(1), particle%st(2), particle%x(3), &
        rho, rho_s, rho_t, rho_phi, rho_time, R, R_s, R_t, Z, Z_s, Z_t)
    combined_energy = rho(1) * MU_ZERO**(3.d0/2.d0) * sqrt(CENTRAL_MASS * CENTRAL_DENSITY * 1.d20 * MASS_PROTON) &
                  * proj_f_combined_energy_SI(sim, i_group, particle) 
  class default
    combined_energy = 0.d0
  end select
end function proj_f_combined_energy

function proj_f_combined_par_momentum_SI(sim, i_group,  particle) result(combined_par_momentum_SI)
  type(particle_sim), intent(in)    :: sim
  class(particle_base), intent(in)  :: particle
  integer, intent(in)               :: i_group
  real*8                            :: combined_par_momentum_SI

  select type (pa => particle)
  type is (particle_kinetic_leapfrog)
    combined_par_momentum_SI = proj_f_ionization_par_momentum(sim, i_group, particle) &
                            +  proj_f_CX_par_momentum(sim, i_group, particle)
  class default
    combined_par_momentum_SI = 0.d0
  end select
end function proj_f_combined_par_momentum_SI

function proj_f_combined_par_momentum(sim, i_group,  particle) result(combined_par_momentum)
  type(particle_sim), intent(in)    :: sim
  class(particle_base), intent(in)  :: particle
  integer, intent(in)               :: i_group
  real*8                            :: combined_par_momentum

  real*8    :: rho(1), R, Z
  real*8    :: rho_s(1), rho_t(1), rho_phi(1), rho_time(1), R_s, R_t, Z_s, Z_t
  real*8    :: E(3), B(3), psi, U, B_norm
 
  select type (pa => particle)
  type is (particle_kinetic_leapfrog)
    call sim%fields%interp_PRZ(sim%time, particle%i_elm, [5], 1, particle%st(1), particle%st(2), particle%x(3), &
        rho, rho_s, rho_t, rho_phi, rho_time, R, R_s, R_t, Z, Z_s, Z_t)
    call sim%fields%calc_EBpsiU(sim%time, particle%i_elm, &
            particle%st, particle%x(3), E, B, psi, U)
    B_norm = norm2(B)
    combined_par_momentum = MU_ZERO / (B_norm * CENTRAL_MASS * MASS_PROTON * rho(1)) &
                         * proj_f_combined_par_momentum_SI(sim, i_group, particle) 
  class default
    combined_par_momentum = 0.d0
  end select
end function proj_f_combined_par_momentum

function proj_f_ion_source(sim, i_group,  particle) result(ion_source)
  type(particle_sim), intent(in)    :: sim
  class(particle_base), intent(in)  :: particle
  integer, intent(in)               :: i_group
  real*8                            :: ion_source

  real*8                            :: n_e, T_e, ion_rate

  select type (pa => particle)
  type is (particle_kinetic_leapfrog)
    call sim%fields%calc_NeTe(sim%time, pa%i_elm, pa%st, pa%x(3), n_e, T_e)
    call sim%groups(i_group)%ad%SCD%interp(pa%q+1, log10(n_e), log10(T_e), ion_rate) 
    ion_source = ion_rate * n_e
  class default
    ion_source = 0d0
  end select
end function proj_f_ion_source

function proj_f_ionization_energy(sim, i_group, particle) result(ion_energy)
  type(particle_sim), intent(in)   :: sim
  class(particle_base), intent(in) :: particle
  integer, intent(in)              :: i_group
  real*8                           :: ion_energy

  real*8, parameter :: binding_energy = 2.18d-18 ! ionization energy of a hydrogen atom [J] (= 13.6 eV)
  real*8            :: ion_source, kinetic_energy
  
  select type (pa => particle)
  type is (particle_kinetic_leapfrog)   
    ion_source     = proj_f_ion_source(sim, i_group, pa)
    kinetic_energy = sim%groups(i_group)%mass * ATOMIC_MASS_UNIT * dot_product(pa%v, pa%v)/2.0d0
    ion_energy     = (kinetic_energy - binding_energy) * ion_source
  class default
    ion_energy = 0.d0
  end select
end function proj_f_ionization_energy

function proj_f_ionization_par_momentum(sim, i_group, particle) result(ion_par_moment)
  type(particle_sim), intent(in)   :: sim
  class(particle_base), intent(in) :: particle
  integer, intent(in)              :: i_group
  real*8                           :: ion_par_moment

  real*8               :: psi, U, ion_source
  real*8, dimension(3) :: E, B, b_hat

  select type (pa => particle)
  type is (particle_kinetic_leapfrog)
    call sim%fields%calc_EBpsiU(sim%time, pa%i_elm, pa%st, pa%x(3), E, B, psi, U)
    b_hat          = B/norm2(B)
    ion_source     = proj_f_ion_source(sim, i_group, pa) 
    ion_par_moment = sim%groups(i_group)%mass * ATOMIC_MASS_UNIT * &
        dot_product(pa%v, b_hat) * ion_source
  class default
    ion_par_moment = 0.d0
  end select
end function proj_f_ionization_par_momentum
  
function proj_f_CX_source(sim, i_group,  particle) result(CX_source)
  type(particle_sim), intent(in)   :: sim
  class(particle_base), intent(in) :: particle
  integer, intent(in)              :: i_group
  real*8                           :: CX_source

  real*8 :: n_e, T_e, CCD

  select type (pa => particle)
  type is (particle_kinetic_leapfrog)
    call sim%fields%calc_NeTe(sim%time, pa%i_elm, pa%st, pa%x(3), n_e, T_e)
    call sim%groups(i_group)%ad%CCD%interp(pa%q+1, log10(n_e), log10(T_e), CCD) ! [m^3/s]
    CX_source = CCD * CENTRAL_DENSITY
  class default
    CX_source = 0d0
  end select
end function proj_f_CX_source

function proj_f_CX_energy(sim, i_group, particle) result(CX_energy)
  type(particle_sim), intent(in)    :: sim
  class(particle_base), intent(in)  :: particle
  integer, intent(in)               :: i_group
  real*8                            :: CX_energy

  real*8               :: n_e, T_e, v_local_squared, delta_v_squared
  
  select type (pa => particle)
  type is (particle_kinetic_leapfrog)
    call sim%fields%calc_NeTe(sim%time, pa%i_elm, pa%st, pa%x(3), n_e, T_e)
    v_local_squared = 2.d0*K_BOLTZ*T_e/(sim%groups(i_group)%mass * ATOMIC_MASS_UNIT)
    delta_v_squared = v_local_squared - dot_product(pa%v, pa%v)
    CX_energy       = 0.5d0 * sim%groups(i_group)%mass * ATOMIC_MASS_UNIT * &
                      delta_v_squared * proj_f_CX_source(sim, i_group, particle)
  class default
    CX_energy = 0.d0
  end select
end function proj_f_CX_energy

function proj_f_CX_par_momentum(sim, i_group, particle) result(CX_par_moment)
  type(particle_sim), intent(in)    :: sim
  class(particle_base), intent(in)  :: particle
  integer, intent(in)               :: i_group
  real*8                            :: CX_par_moment

  real*8               :: psi, U, n_e, T_e
  real*8, dimension(3) :: E, B, b_hat

  ! <STIJN> Note that v-paralel can be gotten from JOREK
  ! <STIJN> Note that maybe different units should be used (v = v*B)
  ! <STIJN> v_perpendicular = 0 (or cross_product(v, B))
  select type (pa => particle)
  type is (particle_kinetic_leapfrog)
    call sim%fields%calc_EBpsiU(sim%time, pa%i_elm, pa%st, pa%x(3), E, B, psi, U)
    b_hat = B/norm2(B)
    call sim%fields%calc_NeTe(sim%time, pa%i_elm, pa%st, pa%x(3), n_e, T_e)
    ! <STIJN> It is assumed that on average the velocity of the ions is 0, in that way delta_v = v_particle
    !v_local       = sqrt(2*K_BOLTZ*T_e/(sim%groups(i_group)%mass * ATOMIC_MASS_UNIT))/sqrt(3.d0)
    CX_par_moment = sim%groups(i_group)%mass * ATOMIC_MASS_UNIT * &
        dot_product(-pa%v, b_hat) * proj_f_CX_source(sim, i_group, particle)
  class default
    CX_par_moment = 0.d0
  end select
end function proj_f_CX_par_momentum

end module mod_projection_functions
