!> Initialise a kinetic relativistic particle group
!> The particle position is initialised using the 
!> normalised plasma density using an acceptance 
!> rejection method. Requires a jorek restart 
!> file named: test_jorek_3d_fields_restart.h5
!> and its input file (to be piped at execution)
program seed_relativistic_particles
  use mpi
  use constants,                             only: PI,TWOPI
  use mod_settings
  use phys_module,                           only: filter_perp, filter_hyper, filter_par
  use phys_module,                           only: filter_perp_n0,filter_hyper_n0,filter_par_n0
  use mod_rng,                               only: type_rng
  use mod_pcg32_rng,                         only: pcg32_rng
  use mod_event
  use mod_particle_types,                    only: particle_kinetic_relativistic
  use mod_particle_sim,                      only: particle_sim
  use mod_fields_linear,                     only: read_jorek_fields_interp_linear
  use mod_project_particles
  use mod_particle_io
  use mod_accept_reject_funct,               only: accept_larger_values_rand
  use mod_initialise_relativistic_particles, only: init_p_gc_relativistic_from_fluid_energypitchgyro
  use mod_fields_minmax,                     only: field_minmax_monte_carlo
  implicit none
  !> variables:
  class(type_rng),dimension(:),allocatable :: rngs
  type(event)                 :: read_field,project_density_event
  type(projection),target     :: project_density
  type(particle_sim)          :: sim
  character(len=20),parameter :: fields_name="test_jorek_3d_fields"
  character(len=33),parameter :: particle_name="test_part_relativistic_init_3d.h5"
  integer*1,dimension(2)              :: q_box
  integer                             :: ii,n_fields,n_groups,n_particles_per_task
  integer                             :: n_trials,ifail
  integer,dimension(:),allocatable    :: fluid_field_ids,particle_field_ids
  real*8                              :: mass_e,max_error
  real*8,dimension(2)                 :: phi_box,energy_kin_box,pitch_box,gyro_box
  real*8,dimension(:,:),allocatable   :: minmax_global_particle,minmax_global_fluid
  real*8,dimension(:,:,:),allocatable :: minmax_list

  !> initialise the simulation ---------------------------------------------------
  n_fields=1; allocate(fluid_field_ids(n_fields)); fluid_field_ids=(/5/); !< fluid density
  allocate(particle_field_ids(n_fields)); particle_field_ids=(/1/) !< particledensity
  phi_box=(/0.d0,TWOPI/); pitch_box=(/0.d0,PI/); gyro_box=(/0.d0,TWOPI/);
  energy_kin_box=(/1.d5,5.d7/); n_groups=1; n_particles_per_task=10000000;
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
  call with(sim,read_field);

  !> initialise the particles ----------------------------------------------------
  call init_p_gc_relativistic_from_fluid_energypitchgyro(&
  sim%groups,sim%fields,sim%time,pcg32_rng(),n_fields,fluid_field_ids,&
  accept_larger_values_rand,phi_box,energy_kin_box,pitch_box,&
  gyro_box,q_box,sim%my_id,sim%n_cpu,ifail,n_trials)

  !> write particle in restart file ----------------------------------------------
  call write_simulation_hdf5(sim,trim(particle_name))

  !> project particle to vtk readable file ----------------------------------------
  project_density = new_projection(sim%fields%node_list,sim%fields%element_list, &
  filter=filter_perp,filter_hyper=filter_hyper,filter_parallel=filter_par, &
  filter_n0=filter_perp,filter_hyper_n0=filter_hyper,filter_parallel_n0=filter_par_n0, &
  f=[proj_f(proj_one,group=1)],fractional_digits = 16,calc_integrals=.true.,to_vtk=.true.,&
  to_h5=.false.,basename='initial_relativistic_particle_density',nsub=5)
  call with(sim,project_density)

  !> verify congruency between normalised particle and fluid profiles -------------
  allocate(minmax_list(n_fields,2,project_density%element_list%n_elements))
  allocate(minmax_global_particle(n_fields,2)); allocate(minmax_global_fluid(n_fields,2));
  call field_minmax_monte_carlo(sim%fields%node_list,sim%fields%element_list,n_fields,&
  fluid_field_ids,n_trials,phi_box,pcg32_rng(),rngs,minmax_list,minmax_global_fluid,&
  sim%my_id,sim%n_cpu,ifail)
  call field_minmax_monte_carlo(project_density%node_list,project_density%element_list,&
  n_fields,particle_field_ids,n_trials,phi_box,pcg32_rng(),rngs,minmax_list,&
  minmax_global_particle,sim%my_id,sim%n_cpu,ifail)
  call error_fluid_particle_profiles(sim%fields%node_list,project_density%node_list,&
  sim%fields%element_list,project_density%element_list,pcg32_rng(),rngs,&
  minmax_global_fluid,minmax_global_particle,n_fields,fluid_field_ids,&
  particle_field_ids,phi_box,max_error,n_trials)
  write(*,*) "Maximum error between fluid and particle profiles: ",max_error

  !> cleanup ----------------------------------------------------------------------
  deallocate(particle_field_ids); deallocate(fluid_field_ids); deallocate(minmax_list); 
  deallocate(minmax_global_particle); deallocate(minmax_global_fluid); call sim%finalize

  !> Program procedure ------------------------------------------------------------
  contains

  !> find maximum error between fluid and particle profiles
  subroutine error_fluid_particle_profiles(node_list_fluid,node_list_part,&
  element_list_fluid,element_list_part,rng_type,rngs,minmax_global_fluid,&
  minmax_global_part,n_fields,fluid_field_ids,part_field_ids,phi_box,&
  max_error,n_trials_in)
    use mod_rng,        only: type_rng,setup_shared_rngs
    use data_structure, only: type_node_list,type_element_list
    use mod_interp,     only: interp_prz
    !$ use omp_lib
    implicit none
    !> inputs-outputs:
    class(type_rng),dimension(:),allocatable,intent(inout) :: rngs
    !> inputs:
    type(type_node_list),intent(in)    :: node_list_fluid,node_list_part
    type(type_element_list),intent(in) :: element_list_fluid,element_list_part
    class(type_rng),intent(in)         :: rng_type
    integer,intent(in)                 :: n_fields
    integer,dimension(n_fields),intent(in)  :: fluid_field_ids,part_field_ids
    real*8,dimension(2),intent(in)          :: phi_box
    real*8,dimension(n_fields,2),intent(in) :: minmax_global_fluid
    real*8,dimension(n_fields,2),intent(in) :: minmax_global_part
    integer,intent(in),optional             :: n_trials_in
    !> outputs:
    real*8,intent(out) :: max_error
    !> variables:
    integer :: ii,n_trials,i_elm,thread_id
    real*8  :: phi,R,Z
    real*8,dimension(4) :: rands
    real*8,dimension(n_fields) :: fluid_fields,part_fields
    !> initialisation
    thread_id=1; max_error=-1.d21; n_trials=1000000;
    if(present(n_trials_in)) n_trials = n_trials_in
    call setup_shared_rngs(size(rands),rng_type,rngs)
    !> loop for checking
    !$omp parallel default(private) firstprivate(n_trials,thread_id,phi_box,&
    !$omp n_fields) shared(node_list_part,node_list_fluid,element_list_part,&
    !$omp element_list_fluid,rngs,minmax_global_part,minmax_global_fluid,&
    !$omp fluid_field_ids,part_field_ids) reduction(max:max_error)
    !$ thread_id = omp_get_thread_num()+1
    !$omp do
    do ii=1,n_trials
      call rngs(thread_id)%next(rands)
      i_elm = 1+floor(real(element_list_fluid%n_elements,kind=8)*rands(1))
      phi = phi_box(1) + (phi_box(2)-phi_box(1))*rands(2)
      call interp_prz(node_list_fluid,element_list_fluid,i_elm,&
      fluid_field_ids,n_fields,rands(3),rands(4),phi,fluid_fields,R,Z)
      call interp_prz(node_list_part,element_list_part,i_elm,&
      part_field_ids,n_fields,rands(3),rands(4),phi,part_fields,R,Z)
      fluid_fields = (fluid_fields-minmax_global_fluid(:,1))/&
      (minmax_global_fluid(:,2) - minmax_global_fluid(:,1))
      part_fields = (part_fields-minmax_global_part(:,1))/&
      (minmax_global_part(:,2) - minmax_global_part(:,1))
      max_error = max(max_error,maxval(abs(fluid_fields-part_fields))) 
    enddo
    !$omp end do
    !$omp end parallel
    !> cleanup
    if(allocated(rngs)) deallocate(rngs)
  end subroutine error_fluid_particle_profiles

end program seed_relativistic_particles

