!> Initialise a kinetic relativistic particle group
!> The particle position is initialised using the 
!> normalised plasma density using an acceptance 
!> rejection method. Requires a jorek restart 
!> file named: test_jorek_3d_fields_restart.h5
!> and its input file (to be piped at execution)
program seed_relativistic_particles
  use constants,                             only: PI,TWOPI
  use mod_pcg32_rng,                         only: pcg32_rng
  use particle_tracer
  use mod_particle_io
  use mod_accept_reject_funct,               only: accept_lower_values_rand
  use mod_initialise_relativistic_particles, only: init_p_gc_relativistic_from_fluid_energypitchgyro
  implicit none
  !> variables:
  type(event) :: read_field
  character(len=20),parameter :: fields_name="test_jorek_3d_fields"
  character(len=33),parameter :: particle_name="test_part_relativistic_init_3d.h5"
  integer*1,dimension(2)           :: q_box
  integer                          :: ii,n_fields,n_groups,n_particles_per_task
  integer                          :: n_trials,ifail
  integer,dimension(:),allocatable :: field_ids
  real*8                           :: mass_e
  real*8,dimension(2)              :: phi_box,energy_kin_box,pitch_box,gyro_box

  !> initialise the simulation ---------------------------------------------------
  n_fields=1; allocate(field_ids(n_fields)); field_ids=(/5/) !< density
  phi_box=(/0.d0,TWOPI/); pitch_box=(/0.d0,PI/); gyro_box=(/0.d0,TWOPI/);
  energy_kin_box=(/1.d5,5.d7/); n_groups=1; n_particles_per_task=1000000;
  mass_e = 5.48579909065d-4; q_box=(/-1,-1/); n_trials=1000000;
  !> set initialisation variables ------------------------------------------------
  call sim%initialize(n_groups,.true.)

  !> basic group setup -----------------------------------------------------------
  do ii=1,n_groups
    sim%groups(ii)%mass = mass_e
    allocate(particle_kinetic_relativistic::sim%groups(ii)%particles(n_particles_per_task))
  enddo

  !> logs ------------------------------------------------------------------------
  write(*,*) "--------------------------------------------------------------------"
  write(*,*) " "
  write(*,*) "Example: initialise relativistic particles from fluid density profile"
  write(*,*) "Number of groups: ",size(sim%groups)
  write(*,*) "Number of particles per task: ",(/(size(sim%groups(ii)%particles),ii=1,size(sim%groups))/)  
  write(*,*) " " 
  write(*,*) "--------------------------------------------------------------------"

  !> load the jorek fields -------------------------------------------------------
  read_field = event(read_jorek_fields_interp_linear(basename=trim(fields_name),i=-1))
  call with(sim,read_field)

  !> initialise the particles ----------------------------------------------------
  call init_p_gc_relativistic_from_fluid_energypitchgyro(&
  sim%groups,sim%fields,sim%time,pcg32_rng(),n_fields,field_ids,&
  accept_lower_values_rand,phi_box,energy_kin_box,pitch_box,&
  gyro_box,q_box,sim%my_id,sim%n_cpu,ifail,n_trials)

  !> write particle in restart file ----------------------------------------------
  call write_simulation_hdf5(sim,trim(particle_name))

  !> cleanup ----------------------------------------------------------------------
  deallocate(field_ids)
end program seed_relativistic_particles

