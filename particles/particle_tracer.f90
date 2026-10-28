module particle_tracer
! Base section
use constants
use mod_particle_sim
use mod_particle_types
use mod_event
use mod_initialise_particles
use initialisers_base

! IO
use mod_io_actions

! Pushers
use mod_boris
use mod_fieldline_euler
use mod_kinetic_relativistic

! Fields
use mod_fields
use mod_fields_linear
use mod_fields_hermite_birkhoff

! ADAS
use mod_openadas
use mod_coronal

! RNGs
use mod_pcg32_rng
use mod_sobseq_rng

! Diagnostics
use mod_diag_print_kinetic_energy
use mod_project_particles
use mod_rhs_projections

! JOREK
use data_structure
use mod_find_rz_nearby

! Default variables
implicit none
type(particle_sim) :: sim
type(event), dimension(:), allocatable, target :: events

public
end module particle_tracer
