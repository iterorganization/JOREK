!> Master module which handles the creation of the initial distribution of
!> superparticles for specific particle groups the for kinetic simulations.
!> This module mainly handles the interface between the initialisers and the simulation
!> The function of the specific initialiser subroutines are contained in the 
!> "initialisers_*.f90" files
!>
!> init_function  -> 1st dispatch, choses global initialiser. 
!>                   Could be only phase-space (e.g. maxwell) or both real and phase-space (e.g. experimental)
!> init_pdf       -> 2nd dispatch. Ignored if init_function includes real-space init (e.g. experimental)
!>                   Optional real-space distribution function. Implemented via rejection sampling (see mod_rej_f).
!> particle_type  -> Determines valid init_functions
!>
!> Valid init_functions:
!>   'maxwell'      -> Maxwellian Energy distribution, uniform pitch. Spatial dist set by init_pdf (or uniform if no init_pdf specified)
!>                     Temperature specified by T_maxwell in input file. Only valid for particle_kinetic_leapfrog
!>   're_gaussian'  -> Monoenergetic or Gaussian energy distribution, delta function in pitch.
!>                     Spatial dist set by init_pdf (or uniform if no init_pdf specified)
!>                     re_energy, re_std_energy, re_pitch specified in input file. Only valid for particle_kinetic_relativistic
!>   'experimental' -> Samples from external F(R, Z, energy, pitch) distribution, expects 'experimental_dist.h5' hdf5 file with the correct format
!>                     real-space dist is implicit in file, init_pdf is ignored. Only valid for particle_kinetic_leapfrog
!>
!>  Any other init_function produces an explicit error at initialisation time.
module mod_initialise_particles
  use mod_rej_f
  use initialisers_RE
  use initialisers_base
  use mod_import_experimental_dist
  use phys_module, only: part_group_configs, type_part_group_config, n_part_groups
  use mod_particle_group_id, only: matching_part_config_indices
  use mpi, only: MPI_ABORT, MPI_COMM_WORLD


  implicit none
  private

  public :: initialise_particles_for_sim

  contains

  !> Top-level initialisation loop over all particle groups
  !> atm skips ics/ncs groups (assume particles born via other methods)
  !> delegates everything else to initialise_group
  subroutine initialise_particles_for_sim(sim)
    implicit none
    class(particle_sim), intent(inout)       :: sim
    integer                                  :: i, j

    do i=1, n_part_groups ! loop over part_groups_in_use
      j = matching_part_config_indices(i) ! get the matching part_group_config index

      select case(trim(part_group_configs(j)%coupling_scheme))
      case('ics', 'ncs', 'non')
        if (sim%my_id == 0) then
          write(*,*) "Particle initialisation skipped for ics/ncs, particles are born elsewhere"
        endif
      case default
        if (sim%my_id == 0) then
          write(*,*) "----- Initialising group '", part_group_configs(j)%id, "' -----"
          write(*,*) "  init_function : '", trim(part_group_configs(j)%init_function), "'"
          write(*,*) "  init_pdf      : '", trim(part_group_configs(j)%init_pdf), "'"
        endif
        call initialise_group(sim, i)
      endselect
    enddo
  end subroutine initialise_particles_for_sim

  !> Initialise a single particle group using namelist defined init_function and init_pdf
  subroutine initialise_group(sim, group_num)
    use mod_pcg32_rng

    implicit none

    !> i/o vars
    class(particle_sim), intent(inout) :: sim
    integer,             intent(in)    :: group_num

    !> internal vars
    type(type_part_group_config) :: config
    type(spatial_pdf)            :: space_pdf
    logical                      :: exists
    integer                      :: ierr

    config = part_group_configs(matching_part_config_indices(group_num))

    !> Compatability check : init_function x particle_type
    select type(particles => sim%groups(group_num)%particles)
    type is (particle_kinetic_relativistic)
      if (trim(config%init_function) /= 're_gaussian') then
        write(*,*) "ERROR : particle_kinetic_relativistic requred init_function='re_gaussian'"
        write(*,*) "        got init_function='",trim(config%init_function),"' for group '", trim(config%id), "'"
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
      endif
    type is (particle_kinetic_leapfrog)
      if (trim(config%init_function) == 're_gaussian') then
        write(*,*) "ERROR (initialise_group): particle_kinetic_leapfrog incompatible with init_function='re_gaussian'"
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
      endif
    end select

    !> Set real-space pdf
    space_pdf = spatial_pdf_from_name(trim(config%init_pdf), sim%fields%node_list, sim%fields%element_list)

    !> Select phase-space initialiser and sample in both real- and phase-space
    select case(trim(config%init_function))

    !> EPs : Maxwellian distribution
    case ('maxwell')
      if (sim%my_id == 0) then
        write(*,*) "  Sampler: initialise_particles_H_mu_psi_phiplanes"
        write(*, '(A,ES12.3,A)') "  T_maxwell    : ", config%T_maxwell, " [eV]"
        write(*, '(A,I0)')       "  n_phi_planes : ", config%n_phi_planes
      endif

      call initialise_particles_H_mu_psi_phiplanes(               &
        sim%groups(group_num)%particles, sim%fields, pcg32_rng(), &
        sim%groups(group_num)%mass, uniform_space=.true.,         & 
        space_pdf=space_pdf, charge=1,                            &
        T_maxwell=config%T_maxwell, n_phi_planes_in=config%n_phi_planes)

    !> EPs : experimental distribution - f(R,Z,energy,pitch), expects a .h5 file - see mod_import_experimental_dist
    case ('experimental')
      if (sim%my_id == 0) write(*,*) "  Sampler: import_particles - experimental distribution from 'experimental_dist.h5' file"

      !> check for existence of experimental_dist.h5 file
      inquire(file="experimental_dist.h5", exist=exists)
      if (.not. exists) then
        if (sim%my_id == 0) write(*,*) "ERROR (initialise_group): init_function='experimental' requires 'experimental_dist.h5' file in working directory"
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
      endif

      call import_particles(sim%groups(group_num)%particles, sim%fields, &
        "experimental_dist.h5", pcg32_rng(), sim%groups(group_num)%mass, &
        n_phi_planes_in=config%n_phi_planes, fraction_phi_planes=1.d0)

    !> REs : Monoenergetic or Gaussian energy distribution
    case('re_gaussian')
      if (sim%my_id == 0) then
        write(*,*) "  Sampler: re_gaussian_initialisation"
        write(*,'(A,ES12.3,A)') "  Energy : ", config%re_energy, " [eV]"
        write(*,'(A,ES12.3,A)') "  std    : ", config%re_std_energy, " [eV]"
        write(*,'(A,ES12.3)')   "  Pitch  : ", config%re_pitch
      endif

      call initialise_re_gaussian(sim, group_num, pcg32_rng(), &
        space_pdf, config%re_energy, config%re_pitch, config%re_std_energy)

    !> Unknown init functions
    case default
      if (sim%my_id == 0) then
        write(*,*) "ERROR (initialise_group): unknown init_function '", &
          trim(config%init_function), "' for group '", trim(config%id), "'"
        write(*,*) "  Valid init_functions: 'maxwell', 're_gaussian', 'experimental'"
      endif
      call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
    
    end select

    !> Finalize weights and charge
    select type(particles => sim%groups(group_num)%particles)
    type is (particle_kinetic_relativistic)
      particles(:)%q      = -1
   end select
  
    !> Finalize weights
    sim%groups(group_num)%particles(:)%weight = 1.d0
    call adjust_particle_weights(sim%groups(group_num)%particles, config%n_particles_total)

    !> Summary print
    if (sim%my_id == 0) then
      write(*,'(A,ES12.3)') "  Number of super particles  : ", config%n_particles
      write(*,'(A,ES12.3)') "  Number of actual particles : ", config%n_particles_total
      write(*,'(A,ES12.3)') "  Particle weight            : ", sim%groups(group_num)%particles(1)%weight
      write(*,*) "----- Finished initialisation for group '", trim(config%id), "'-----"
      write(*,*) ""
    endif

  end subroutine initialise_group

end module mod_initialise_particles