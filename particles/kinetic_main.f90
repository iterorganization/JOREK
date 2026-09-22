!> The equivalent of jorek2_main.f90 for kinetic simulations
!>
!> OpenADAS is used for atomic physics
!> The wall interaction rates are based on external files with Eckstein coefficients such as 
!> y_DD.dat (for yield of a D -> D reaction, which would be wall recombination) and 
!> ye_DD.dat (for the resulting energy of the same reaction).
!>
!> To use a particle restart file: use restart_particles=.t. in the input file.
!>
!> To change to 3D recombination:
!> Change subroutine do1particlerecombination: use mod_integrate_recomb3D, only : integrate_recombination
!> Particle puffing is axisymmetric by default.

program kinetic_main

use particle_tracer
use mod_particle_diagnostics
use mpi
use mod_interp
use mod_atomic_elements
use mod_particle_evolution
use mod_particle_recomb
use mod_particle_conservation
use mod_particle_io
use mod_particle_allocation
use mod_event
use mod_project_particles
use mod_jorek_timestepping
use mod_random_seed
use mod_basisfunctions
use nodes_elements
use constants,   only: MU_ZERO, ATOMIC_MASS_UNIT, K_BOLTZ, EL_CHG
use mod_particle_wall_interaction
use mod_neutral_collision, only: neutral_collisions_from_config, type_neutral_collision, gcd_neutral_collisions
use mod_projection_functions, only: proj_f_combined_density, proj_f_combined_energy, proj_f_combined_par_momentum
use mod_particle_puffing
use mod_edge_domain
use mod_edge_elements, only: edge_elements
use mod_atomic_coeff_deuterium, only: ad_deuterium 
use data_structure, only: type_bnd_element_list, type_bnd_node_list 
use mod_boundary,   only: boundary_from_grid
use mod_coupling_settings, only: use_kin_recomb_global
use mod_initialise_particles
use equil_info
use mod_output_file_routines, only: write_to_outputfile
use mod_neutral_density, only: get_neutral_density
use mod_impurity,        only: init_imp_adas

use phys_module, only: index_now
use phys_module, only: tstep,tstep_n,restart_particles, restart, t_start, nout
use phys_module, only: CENTRAL_MASS, CENTRAL_DENSITY, xcase, xpoint
use phys_module, only: n_part_groups, n_aux_var, n_valves_max
use phys_module, only: nstep_particles, nsubstep_particles, tstep_particles, nout_particles
use phys_module, only: deuterium_adas,sqrt_mu0_over_rho0
use phys_module, only: filter_perp, filter_hyper, filter_par, filter_perp_n0, filter_hyper_n0, filter_par_n0
use phys_module, only: apply_dirichlet_proj, part_group_configs, init_particles_only
use phys_module, only: use_manual_random_seed, manual_seed
use phys_module, only: use_imp_adas, nimp_bg 

use mod_particle_group_id, only: matching_part_config_indices

use mod_pcg32_rng, only: pcg32_rng
use mod_rng, only: type_rng, setup_shared_rngs
!$ use omp_lib

implicit none

type(event)                                       :: fieldreader, partreader
type(event)                                       :: gas_puff_event, gas_puff2_event
type(event), target                               :: project_jorek_feedback, jorek_stepper_event, project_neutral_density
type(pcg32_rng), dimension(:), allocatable        :: rng
type(count_action)                                :: counter
type(projection), target                          :: jorek_feedback
type(projection), target                          :: neutral_density_proj
type(jorek_timestep_action), target               :: jorek_stepper
type(type_edge_domain), allocatable, dimension(:) :: edge_domains
type(edge_elements)                               :: edge_elm_template
type(particle_puffing)                            :: gas_puff, gas_puff2
character(len=50)                                 :: rst_part_file

real*8    :: rho_norm, t_norm, n_norm

integer   :: n_reflect
integer   :: i, j, istep_inner_loop, group_num, config_num, valve_num, n_lcm_blocks, inner_stepsize
integer   :: seed, i_rng, n_stream
logical   :: last_step !< whether it is the last step of the inner particle loop
logical   :: puff_from_namelist, puff_from_file

!> For keeping track of groups requiring specific physics (e.g. wall actions, recombination, puffing...)
integer   :: n_wall_act_groups = 0
integer   :: recomb_counter  = 0
integer   :: puff_counter    = 0

integer,                      dimension(:), allocatable :: recomb_groups
type(particle_puffing),       dimension(:), allocatable :: puff_actions   
type(wall_act_group),         dimension(:), allocatable :: wall_act_groups
type(type_neutral_collision), dimension(:), allocatable :: neutral_collisions

!tmp
class(type_rng), dimension(:), allocatable :: wall_rng

integer :: n_particles_local

character(len=100) :: header_line

!***********************************************************************
!*                            initialisation                            *
!***********************************************************************

! Start up MPI, jorek
call sim%initialize()

! Loading the jorek fields
if (restart) then
  fieldreader = event(read_jorek_fields_interp_linear(basename='jorek', i=-1))
  call with(sim, fieldreader)
else
  if (sim%my_id == 0) write(*,*) 'ERROR: using this program without restarting from a jorek field is not possible. Please set restart=.t. in the namelist and provide a jorek_restart.h5 file'
  stop
end if

!> update and broadcast the equilibrium state
if (sim%my_id .eq. 0) call boundary_from_grid(sim%fields%node_list, sim%fields%element_list, bnd_node_list, bnd_elm_list, .false.)
call broadcast_boundary(sim%my_id, bnd_elm_list, bnd_node_list)
call update_equil_state(sim%my_id, sim%fields%node_list, sim%fields%element_list, bnd_elm_list, xpoint, xcase )
call broadcast_equil_state(sim%my_id)

! setting up the particles
if (restart_particles) then
  ! reading the particles from a file
  if (sim%my_id == 0) write(*,*) 'INFO: READING PARTICLES RESTART FILE =========='
  partreader = event(read_action(filename='part_restart.h5'))
  call with(sim, partreader) !<defines sim%groups and the corresponding particles

  !TODO? Sven: We should make an option to use partreader but increase n_particles; may be similar to phi_zero_whrite to a sim_in and sim_out but with different allocation size.
else
  if (sim%my_id == 0) write(*,*) 'INFO: INITIALIZING PARTICLES', sim%n_mpi, " mpi's "

  call allocate_particles_for_sim(sim) ! populate the particle arrays in the particle groups

  call initialise_particles_for_sim(sim) !< initialise the particle positions for groups requiring it (e.g. REs)

  call write_simulation_hdf5(sim, 'part_restart_init.h5')

  if (init_particles_only) then
    call sim%finalize
    stop
  endif

  
endif ! (restart_particles)


! Read Open ADAS data for plasma fluid
if (deuterium_adas .and. use_kin_recomb_global) ad_deuterium =  read_adf11(sim%my_id,'96_h') !< move to core (jorek2_main for particles)

#if (defined WITH_Impurities)
  ! --- Read ADAS data and generate coronal equilibrium if needed
  call init_imp_adas(sim%my_id)
#else
  if (use_imp_adas .and. (nimp_bg(1) > 0.d0)) then
    call init_imp_adas(sim%my_id)
  endif
#endif


! --- Setting up random numbers for ionisation probability
seed = random_seed()
n_stream = 1
!$ n_stream = omp_get_max_threads()
write(*,*) "id, n_mpi, n_stream",sim%my_id, sim%n_mpi, n_stream
allocate(rng(n_stream))
do i=1,n_stream
  call rng(i)%initialize(1, seed, n_stream, i)
end do

! --- Check if tstep_particles is positive
if (tstep_particles <= 0.d0) then
  if (sim%my_id == 0) write(*,*) "ERROR: tstep_particles <= 0 which is not allowed. Stopping now."
  stop
end if

!***********************************************************************
!*                         setting up physics                          *
!***********************************************************************

! --- Calculating normalisation constants
n_norm    = CENTRAL_DENSITY * 1.d20                              ! (number) density normalisation
rho_norm  = CENTRAL_MASS * ATOMIC_MASS_UNIT * n_norm             ! rho_SI = rho_norm * rho
t_norm    = sqrt((MU_ZERO * rho_norm))                           ! t_SI   = t_norm * t_jorek 

! --- Setting up wall actions
! edge domains for wall actions
call find_edge_domains(node_list,element_list, edge_domains)!, discont_corner=.true.)
if (sim%my_id .eq. 0) write(*,*) "n_domains = ", size(edge_domains,1)
call edge_elm_template%prepare(node_list, element_list, edge_domains, nsub=6, nsub_toroidal=1)!,wall_albedo=wall_albedo)

! getting the wall actions from the input
wall_act_groups =  wall_actions_from_config(sim, edge_elm_template)
n_wall_act_groups = size(wall_act_groups,1)

! --- Setting up particle events
allocate(recomb_groups(n_part_groups), puff_actions(n_part_groups*n_valves_max)) 

do group_num=1, n_part_groups
  config_num = matching_part_config_indices(group_num)
  
  ! puffing
  if (sim%groups(group_num)%use_kin_puffing) then
    do valve_num=1, n_valves_max
      !> check for the puff_ctrl objects that have been set
      
      puff_from_namelist = count(part_group_configs(group_num)%puff_ctrl(valve_num)%rates /= -1) + count(part_group_configs(group_num)%puff_ctrl(valve_num)%times /= -1) > 0
      puff_from_file = trim(part_group_configs(group_num)%puff_ctrl(valve_num)%from_file) /= "none"
      if (puff_from_namelist .or. puff_from_file) then
        puff_counter = puff_counter + 1   ! increase the number of puffing events required per fluid time step
        puff_actions(puff_counter) = particle_puffing(sim, group_num, valve_num, seed=seed)  
      endif
    enddo ! n_valves_max
  endif

  ! recombination
  if (sim%groups(group_num)%use_kin_recombination) then
    recomb_counter = recomb_counter + 1   ! increase the number of groups that requires recombination
    recomb_groups(recomb_counter) = group_num
  endif

enddo ! n_part_groups

! neutral collisions
neutral_collisions = neutral_collisions_from_config(sim)

! --- Set up feedback to the plasma (does not currently include recombination)
jorek_feedback = new_projection(sim%fields%node_list, sim%fields%element_list, &
                                filter_n0 = filter_perp_n0, filter_hyper_n0 = filter_hyper_n0, filter_parallel_n0=filter_par_n0,            &
                                filter = filter_perp, filter_hyper = filter_hyper, filter_parallel=filter_par, fractional_digits = 9,       &
                                do_zonal = .false., calc_integrals=.false., to_vtk=.false., to_h5 = .false., basename='projections', nsub=2, &
                                do_dirichlet=apply_dirichlet_proj)
aux_node_list => jorek_feedback%node_list

!> define feedback size dependent on the number of variables required for coupling
allocate(jorek_feedback%rhs(n_order+1, n_vertex_max, sim%fields%element_list%n_elements, n_tor, n_aux_var))

jorek_feedback%rhs = 0.d0

! --- Setting up jorek timestepper
! For proper timestepping, the projections need to be defined before the jorek timestepper
jorek_stepper = new_jorek_timestep_action(jorek_feedback%node_list)

project_jorek_feedback = new_event_ptr(jorek_feedback,   start = sim%time)
jorek_stepper_event    = new_event_ptr(jorek_stepper,    start = sim%time)

! setting up projection for neutral density
if(size(neutral_collisions) > 0) then
  neutral_density_proj = new_projection(sim%fields%node_list, sim%fields%element_list, &
                                  filter_n0 = filter_perp_n0, filter_hyper_n0 = filter_hyper_n0, filter_parallel_n0=filter_par_n0,            &
                                  filter = filter_perp, filter_hyper = filter_hyper, filter_parallel=filter_par, fractional_digits = 9,       &
                                  do_zonal = .false., calc_integrals=.false., to_vtk=.false., to_h5 = .false., basename='projections', nsub=2, &
                                  do_dirichlet=apply_dirichlet_proj)

  ! define feedback size dependent on the number of variables required for coupling
  allocate(neutral_density_proj%rhs(n_order+1, n_vertex_max, sim%fields%element_list%n_elements, n_tor, 1))
endif

!if no %each_nstep_part was set, the least common multiple of all %each_nstep_part is 1 (meaning all particle-particle actions will happen only once each fluid timestep)
if(sim%lcm_inner_loop == -9999991) sim%lcm_inner_loop = 1

!***********************************************************************
!*                           main loop                                 *
!***********************************************************************

sim%istep_fluid = 0
call write_to_outputfile(sim,"Starting main loop",next_block_write_conserv=.false.,next_block_write_timing=.false.) ! next_block_write_...=.false. because this is only a header, not the announcement of some action, so we don't want to time or write particle conservation for the "content" of this block as there is no content
do while (.not. sim%stop_now)
  sim%istep_fluid = sim%istep_fluid + 1
  write(header_line,'(A,I6)') "Starting main loop iteration sim%istep_fluid = ",sim%istep_fluid
  call write_to_outputfile(sim,header_line,next_block_write_conserv=.false.,next_block_write_timing=.false.,call_main_loop_timer=.true.) ! call_main_loop_timer=.true. because we want to write the total time taken by the previous main loop iteration before starting the next main loop iteration

  ! --- Determining the time stepping for this fluid step
  tstep = get_tstep_n(sim%istep_fluid)     ! tstep is also set in stepper, but tstep is already used in the calls before the stepper
  sim%tstep_fluid_si = tstep*t_norm
  sim%time = sim%time + sim%tstep_fluid_si ! carries the time at the end of the current main loop

  if (sim%my_id .eq. 0) then
    write(*,*) "sim%time       : ",sim%time
    write(*,*) "tstep_fluid_si : ",sim%tstep_fluid_si
  endif

  ! --- Interactions that happen on the fluid timestep (creating kinetic particles)

  if (n_wall_act_groups > 0) then
    ! here only do the wall actions that create particles
    call write_to_outputfile(sim, "Particle creating wall_actions")
    do i=1, n_wall_act_groups
      call wall_act_groups(i)%do(sim,.false.)
    enddo
  endif

  if (recomb_counter > 0) then
    call write_to_outputfile(sim, "Volume recombination") !< as opposed to wall recombination which is part of wall_actions
    do i=1, recomb_counter
      call do_1particle_recombination(element_list,node_list, recomb_groups(i), jorek_stepper,rng, sim%tstep_fluid_si) 
    enddo
  endif
    
  if (puff_counter > 0) then 
    call write_to_outputfile(sim, "Puffing")
    if (sim%my_id == 0) write(*,"(A,G13.6,A)") "====== Puffing details for time t=", sim%time, " s ======"
    do i=1, puff_counter
      call puff_actions(i)%do(sim)
    enddo
  endif

  ! --- Inner particle loop, with interactions that happen on the particle timesteps
  
  call write_to_outputfile(sim,"Starting inner particle loop",next_block_write_conserv=.false.,next_block_write_timing=.false.)
  
  ! We divide the fluid timestep into an integer number of particle timesteps such that all actions with their own %each_nstep_part fit an integer amount 
  ! of times in the fluid timestep. We ensure this by fitting the least common multiple (lcm) of all %each_nstep_part into the fluid timestep an integer 
  ! amount of times. 
  n_lcm_blocks             = ceiling(sim%tstep_fluid_si / (sim%lcm_inner_loop * tstep_particles)) ! integer number of lcm blocks within one fluid step
  sim%nstep_inner_loop     = n_lcm_blocks*sim%lcm_inner_loop                                      ! number of inner loop steps
  sim%tstep_part_adj       = sim%tstep_fluid_si / sim%nstep_inner_loop                            ! slightly smaller tstep_particles to fit an exact integer amount in one fluid timestep
  
  ! if gcd > 1, we don't have to run each particle loop separately, but we can bunch them up into steps of size gcd
  ! if no %each_nstep_part is set, gcd is the maximum possible value of sim%nstep_inner_loop
  inner_stepsize = sim%gcd_inner_loop
  if(sim%gcd_inner_loop == -9999991) inner_stepsize = sim%nstep_inner_loop

  if (sim%my_id .eq. 0) then
    if(sim%tstep_fluid_si < tstep_particles) then
      write(*,"(A)") "WARNING: tstep < tstep_particles which is currently not supported. Effectively tstep_part_adj = tstep will be used,"
      write(*,"(A)") "or, if the least common multiple (lcm) > 1, tstep_part_adj < tstep will be used."
    else if(sim%lcm_inner_loop * tstep_particles > sim%tstep_fluid_si) then
      write(*,"(A)") "WARNING: the least common multiple (lcm) of all specified %each_nstep_part makes the particle_timestep smaller than it"
      write(*,"(A)") "needs to be this fluid tstep. Consider changing your %each_nstep_part to be more compatible with eachother (e.g. avoid   "
      write(*,"(A)") "co-primes), and smaller in general to keep the lcm small."
    end if
    write(*,*) "sim%time                   : ",sim%time
    write(*,*) "tstep_fluid_si             : ",sim%tstep_fluid_si
    write(*,*) "aim tstep_particles        : ",tstep_particles
    write(*,*) "used tstep_part_adj        : ",sim%tstep_part_adj
    write(*,*) "nstep_inner_loop           : ",sim%nstep_inner_loop
    write(*,*) "lcm, n_lcm_blocks          : ",sim%lcm_inner_loop, n_lcm_blocks
    write(*,*) "gcd, inner_stepsize        : ",sim%gcd_inner_loop, inner_stepsize
    write(*,*) "n_part*dt_part - dt_fluid  : ",sim%nstep_inner_loop*sim%tstep_part_adj - sim%tstep_fluid_si
  endif
  
  !> As we call multiple kinetic loops multiple times and only want to use 1 %rhs,
  !> we should set it to zero here before the inner loop starts
  jorek_feedback%rhs = 0.d0

  do istep_inner_loop=inner_stepsize,sim%nstep_inner_loop,inner_stepsize
    write(header_line,'(A,I6,A,I6)') "Starting inner particle loop iteration getting us to istep_inner_loop=",istep_inner_loop," out of ",sim%nstep_inner_loop
    call write_to_outputfile(sim,header_line,next_block_write_conserv=.false.,next_block_write_timing=.false.)

    !updating inner loop steps
    sim%istep_inner_loop = istep_inner_loop
    if(istep_inner_loop == sim%nstep_inner_loop) then
      last_step = .true.
    else
      last_step = .false.
    endif

    call write_to_outputfile(sim, "Particle evolution loop")
    
    !> Evolution loop: calculating interaction between particles in a group and its environment (i.e. CX, ionisation) 
    !> + pushing the particles + calculating the feedback
    !> Currently supports: (Work in progress)
    !>  - kinetic neutrals
    !>  - kinetic impurities

    ! evolution loop is called every inner particle step, for inner_stepsize number of steps at once
    do group_num=1, n_part_groups
      call evolve_particle_group(sim, group_num, jorek_feedback, rng, sim%tstep_part_adj, inner_stepsize)
    enddo  

    ! --- Handling the particles that left the domain

    ! At the moment the inner particle loop does not include any particle generation routines (such as puffing). If such routines would be moved 
    ! to the inner particle loop, it would reintroduce the problem of overwriting %ielm<0 particles (particles that hit the wall but have not been 
    ! reintroduced into the domain yet during the below particle-particle wall_actions) because the particle-particle wall_actions don't necessarily 
    ! happen directly after the evolve_particle_group due to the flexible timestepping.
    if (n_wall_act_groups > 0 .and. (mod(istep_inner_loop,gcd_wall_acts)==0 .or. last_step)) then
      call write_to_outputfile(sim, "Particle particle wall_actions")
      do i=1, n_wall_act_groups
        call wall_act_groups(i)%do(sim,.true.)
      enddo
    endif
  
    !neutral self collisions, which need the projected neutral density
    if (size(neutral_collisions) > 0 .and. (mod(istep_inner_loop,gcd_neutral_collisions)==0 .or. last_step)) then
      call write_to_outputfile(sim, "Neutral self collisions")
  
      ! update the neutral density which are necessary for the neutral collisions
      call get_neutral_density(sim,neutral_density_proj)

      ! do the neutral collisions
      do i=1, size(neutral_collisions)
        call neutral_collisions(i)%do(sim,neutral_density_proj%node_list,neutral_density_proj%element_list)
      enddo
    endif
    
  enddo ! inner particle loop
  sim%istep_inner_loop=-1 ! set to -1 to notify that sim is outside of inner particle loop

  ! --- Update the fluid
  
  !> Project the collected feedback from the particles onto the finite element grid so that the MHD solver can use it
  !> Also writes the projection.vtk file which contains the interaction terms (particle, energy and momentum exchange to the fluid) and neutral density
  call write_to_outputfile(sim, "Projecting feedback from particles to fluid FE grid",next_block_write_conserv=.false.) ! next_block_write_conserv=.false. since the particle array is not changed during this action, so we don't want to log the global #particles, momentum and energy here
  call with(sim, project_jorek_feedback)
  
  !> Calls the MHD solver which timesteps the MHD fluid based on the fluid itself using the projected
  !> feedback of the particles as sources and sinks in the MHD equations 
  !> Also writes .h5 file, updates tstep and sets sim%stop_now = .true. if all fluid steps are done
  call write_to_outputfile(sim, "Fluid stepper",next_block_write_conserv=.false.)
  call with(sim, jorek_stepper_event) 
  

  ! -- Finalising the fluid timestep
  
  !Writing interim particle restart files every nout fluid steps done. Overwrites previous restart file to save space
  if ( nout_particles .eq. 9999999 ) then
    if ( mod(index_now,nout) .eq. 0 ) then
      call write_to_outputfile(sim, "Writing interim_part_restart.h5",next_block_write_conserv=.false.)
    endif
  else
    !Writing regular particle restart files every nout_particles.
    if ( mod(index_now,nout_particles) .eq. 0 ) then
      write(rst_part_file, '(a,i6.6,a)') 'part_restart_', index_now, '.h5'
      write(header_line, "(2A)") "Writing ",trim(rst_part_file)
      call write_to_outputfile(sim, header_line,next_block_write_conserv=.false.)
      call write_simulation_hdf5(sim, trim(rst_part_file))
    endif
  endif

  ! Writing some conservation checks to the ouput file
  call write_to_outputfile(sim, "Conservation checks")
  call conservation_checks(sim)

end do ! while

!***********************************************************************
!*                          end of simulation                          *
!***********************************************************************
call write_to_outputfile(sim, "End of simulation",next_block_write_conserv=.false.,next_block_write_timing=.false.,call_main_loop_timer=.true.)

call write_to_outputfile(sim, "Writing final part_restart.h5",next_block_write_conserv=.false.)
call write_simulation_hdf5(sim, 'part_restart.h5')

call write_to_outputfile(sim, "Finalizing",next_block_write_conserv=.false., next_block_write_timing=.false.)
deallocate(recomb_groups)
call sim%finalize()

!***********************************************************************
!*                          end of main program                        *
!***********************************************************************

contains

pure function f_adapted(n, P, grad_P) result(f)
  integer, intent(in) :: n
  real*8, intent(in) :: P(n), grad_P(3,n)
  real*8 ::s, coeff(0:4)
  real*4 :: f

  coeff(0)=0.53
  coeff(1)=0.3
  coeff(2)=0.2
  coeff(3)=0.52
  coeff(4)=0.26

  s = max((P(1) - ES%Psi_axis) / ( ES%Psi_bnd - ES%Psi_axis),0.d0)

  f = coeff(3)*exp(-coeff(2)/coeff(1)*(tanh((sqrt(s)-coeff(0))/coeff(2))))

  f = (f - coeff(4)) / (1.d0 - coeff(4))

end function f_adapted

pure function f_original(n, P, grad_P) result(f)
  integer, intent(in) :: n
  real*8, intent(in) :: P(n), grad_P(3,n)
  real*8 ::s, coeff(0:3)
  real*4 :: f

  ! central densiy should be 1.44131x10^17

  coeff(0)=0.49123
  coeff(1)=0.298228
  coeff(2)=0.198739
  coeff(3)=0.521298

  s = max((P(1) - ES%Psi_axis) / ( ES%Psi_bnd - ES%Psi_axis),0.d0)

  f = coeff(3)*exp(-coeff(2)/coeff(1)*(tanh((sqrt(s)-coeff(0))/coeff(2))))

end function f_original


end program kinetic_main
