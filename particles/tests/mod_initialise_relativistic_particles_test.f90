!> test the initialisation routines for relativistic particles
module mod_initialise_relativistic_particles_test
use fruit
use constants,             only: PI,TWOPI
use mod_particle_sim,      only: particle_group
use mod_particle_types,    only: particle_kinetic_relativistic_id
use mod_particle_types,    only: particle_gc_relativistic_id
use mod_particle_types,    only: particle_kinetic_id
use mod_particle_types,    only: particle_kinetic_leapfrog_id
use mod_particle_types,    only: particle_gc_vpar_id
use mod_particle_types,    only: particle_kinetic_relativistic
use mod_particle_types,    only: particle_gc_relativistic
use mod_particle_types,    only: particle_kinetic
use mod_particle_types,    only: particle_kinetic_leapfrog
use mod_particle_types,    only: particle_gc_vpar
use mod_fields_linear,     only: jorek_fields_interp_linear
use mod_fields_analytical, only: fields_analytical
implicit none

private
public :: run_fruit_initialise_relativistic_particles

!> Variables and datatypes ------------------------------------
character(len=20),parameter :: fields_eq_name_sol="test_jorek_3d_fields"
integer,parameter  :: n_v=1
integer,parameter  :: n_groups=2
integer,parameter  :: n_groups_2=5
integer,parameter  :: n_samples=343
integer,parameter  :: n_fields_sol=1
integer,parameter  :: n_active_groups_2_sol=2
integer,parameter  :: n_trials_sol=1000
integer,parameter,dimension(n_v)           :: i_psi=(/1/)
integer,dimension(n_groups),parameter      :: n_particles=(/123,234/)
integer,dimension(n_groups),parameter      :: p_types_sol=(/&
                      particle_kinetic_relativistic_id,&
                      particle_gc_relativistic_id/)
integer,dimension(n_groups_2),parameter    :: n_particles_2=(/1,1,1,1,1/)
integer,dimension(n_groups_2),parameter    :: p_types_2_sol=(/&
                      particle_kinetic_id,&
                      particle_kinetic_relativistic_id,&
                      particle_kinetic_leapfrog_id,&
                      particle_gc_relativistic_id,&
                      particle_gc_vpar_id/)
integer,dimension(n_fields_sol),parameter  :: field_ids_sol=(/5/) !< index of the density
integer,dimension(n_groups_2),parameter    :: active_group_ids_2_sol=(/2,4,0,0,0/)
real*8,parameter                           :: tol_real8=5.d-13
real*8,parameter                           :: tol_interp_real8=7.5d-10
real*8,parameter                           :: time_sol=0.d0
real*8,parameter                           :: p_neg=-2.d0
real*8,dimension(2),parameter              :: psi_big_sol=(/-2.d1,3.d2/)
real*8,dimension(2),parameter              :: psi_small_sol=(/3.d-1,6.d-1/)
real*8,dimension(2),parameter              :: psi_axisbnd=(/-1.d0,2.5d0/)
real*8,dimension(2),parameter              :: R_big_sol=(/1.d-1,3.d2/)
real*8,dimension(2),parameter              :: R_small_sol=(/1.d0,2.5d0/)
real*8,dimension(2),parameter              :: Z_big_sol=(/-3.d1,3.d1/)
real*8,dimension(2),parameter              :: Z_small_sol=(/-1.d0,2.5d0/)
real*8,dimension(2),parameter              :: theta_big_sol=(/-3.d1,3.d1/)
real*8,dimension(2),parameter              :: theta_small_sol=(/PI/6.d0,TWOPI/3.d0/)
real*8,dimension(2),parameter              :: phi_big_sol=(/-3.d1,3.d1/)
real*8,dimension(2),parameter              :: phi_small_sol=(/PI/6.d0,TWOPI/3.d0/)
real*8,dimension(2),parameter              :: pitch_small_sol=(/PI/6.d0,2.d0*PI/3.d0/)
real*8,dimension(2),parameter              :: gyro_small_sol=(/PI/6.d0,3.d0*PI/2.d0/)
real*8,dimension(2),parameter              :: pitch_big_sol=(/-PI/2.d0,TWOPI/)
real*8,dimension(2),parameter              :: gyro_big_sol=(/-PI/2.d0,3.d0*PI/)
real*8,dimension(2),parameter              :: psi_minmax_sol=(/-5.d0,1.d1/)
real*8,dimension(2),parameter              :: R_minmax=(/5.d-1,1.d1/)
real*8,dimension(2),parameter              :: Z_minmax=(/-1.d1,1.d1/)
real*8,dimension(2),parameter              :: theta_minmax=(/0.d0,TWOPI/)
real*8,dimension(2),parameter              :: phi_minmax=(/0.d0,TWOPI/)
real*8,dimension(2),parameter              :: R_box_fraction_sol=(/0.22,0.34/)
real*8,dimension(2),parameter              :: Z_box_fraction_sol=(/0.41,0.73/)
class(particle_group),dimension(:),allocatable :: groups_sol
type(particle_group),dimension(n_groups_2)     :: groups_2_sol
type(fields_analytical)                        :: fields_sol
type(jorek_fields_interp_linear)               :: fields_linear_sol
real*8,dimension(2)                            :: R_box_jorek_sol
real*8,dimension(2)                            :: Z_box_jorek_sol
real*8,dimension(2)                            :: R_minmax_jorek_sol
real*8,dimension(2)                            :: Z_minmax_jorek_sol
real*8,dimension(2)                            :: psi_minmax_global_jorek_sol
real*8,dimension(2)                            :: field_minmax_global_sol
real*8,dimension(:,:),allocatable              :: psi_minmax_list_jorek_sol
real*8,dimension(:,:,:),allocatable            :: field_minmax_list_sol
!> Interfaces--------------------------------------------------
contains
!> Fruit basket -----------------------------------------------
!> fruit test basket: set-up, run and tear-down tests
subroutine run_fruit_initialise_relativistic_particles()
  implicit none
  write(*,'(/A)') "  ... setting-up: initialise relativistic particles tests"
  call setup
  call setup_jorek_simulation
  call setup_init_intervals
  call setup_field_minmax
  write(*,'(/A)') "  ... running: initialise relativistic particles tests"
  call test_dummy
  call test_find_relativistic_kinetic_gc_groups
  call test_find_RZPhi_minmax_global
  call test_find_psithetaphi_minmax_global_list
  call test_sample_position_uniformly_psi_theta_phi
  call test_sample_position_uniformly_cylinder
  call test_sample_position_acceptreject_from_fluid_profiles
  call test_init_p_gc_relativistic_psithetaphi_energypitchgyro
  call test_init_p_gc_relativistic_RZPhi_energypitchgyro
  call test_init_p_gc_relativistic_from_fluid_energypitchgyro
  call test_sampling_cartesian_p_kinetic_relativistic
  call test_sampling_cartesian_gc_kinetic_relativistic
  call test_sampling_uniform_ppitchgyro_kinetic_relativistic
  call test_sampling_uniform_ppitchgyro_gc_relativistic
  call test_particle_base_init_to_zero
  call test_sampling_uniform_charge
  call test_check_thetapsiphi_interval
  call test_check_RZPhi_interval
  call test_check_energykinpitchgyro_interval
  write(*,'(/A)') "  ... tearing-down: initialise relativistic particles tests"
  call teardown
  call teardown_field_minmax
end subroutine run_fruit_initialise_relativistic_particles

!> Set-up and tear-down ---------------------------------------
!> set-up relativistic particle initialisation test features
subroutine setup()
  use mod_gnu_rng,                    only: gnu_rng_interval
  use mod_particle_common_test_tools, only: fill_mass_RE
  use mod_particle_common_test_tools, only: rng_seed_interval
  use mod_particle_common_test_tools, only: allocate_one_particle_list_type
  use mod_particle_common_test_tools, only: RZ0_lowbnd,RZ0_uppbnd
  use mod_particle_common_test_tools, only: BE0_lowbnd,BE0_uppbnd
  implicit none
  !> variables:
  integer :: ifail
  integer,dimension(0) :: int_param
  real*8,dimension(4)  :: real_param 
  !> allocate groups
  allocate(groups_sol(n_groups))
  !> allocate particle lists
  call allocate_one_particle_list_type(n_groups,n_particles,&
  p_types_sol,groups_sol,ifail)
  call allocate_one_particle_list_type(n_groups_2,n_particles_2,&
  p_types_2_sol,groups_2_sol,ifail)
  !> initialise the group masses as runaway
  call fill_mass_RE(n_groups,groups_sol)
  !> initialise particle fields
  call gnu_rng_interval(2,RZ0_lowbnd,RZ0_uppbnd,real_param(1:2))
  call gnu_rng_interval(2,BE0_lowbnd,BE0_uppbnd,real_param(3:4))
  call fields_sol%init_fields(0,4,int_param,real_param)
end subroutine setup

!> set-up a jorek simulation reading directly from a input
!> and restart file. No mpi is expected for this tests hence
!> the mpi rank is set to 0
subroutine setup_jorek_simulation()
  use mod_fields_linear, only: read_jorek_fields_interp_linear
  use mod_particle_sim,  only: particle_sim
  use mod_event,         only: event,with
  use phys_module,       only: xpoint,xcase
  use data_structure,    only: type_bnd_node_list
  use data_structure,    only: type_bnd_element_list
  use mod_boundary,      only: boundary_from_grid
  use equil_info,        only: update_equil_state
  implicit none
  !> variables
  type(type_bnd_node_list)        :: bnd_node_list
  type(type_bnd_element_list)     :: bnd_element_list
  type(particle_sim)              :: sim_particles
  type(event),dimension(1),target :: events
  !> initialise particle simulation and all jorek fields
  call sim_particles%initialize(0,.true.); sim_particles%time=time_sol;
  !> set-up the read event
  events = [event(read_jorek_fields_interp_linear(basename=fields_eq_name_sol,i=-1))]
  !> read the fields from jorek restart file
  call with(sim_particles,events,at=0.d0)
  !> copy particle field to external structure
  select type(fields=>sim_particles%fields)
  type is (jorek_fields_interp_linear)
    fields_linear_sol = fields
  end select
  !> find the boundary nodes and elements
  call boundary_from_grid(fields_linear_sol%node_list,fields_linear_sol%element_list,&
  bnd_node_list,bnd_element_list,.false.)
  !> compute the equilibrium
  call update_equil_state(0,fields_linear_sol%node_list,&
  fields_linear_sol%element_list,bnd_element_list,xpoint,xcase)
end subroutine setup_jorek_simulation

!> setup initialisation intervals
subroutine setup_init_intervals()
  implicit none
  !> variables:
  integer :: ii
  real*8 :: R_min,R_max,Z_min,Z_max
  !> initialisations
  allocate(psi_minmax_list_jorek_sol(fields_linear_sol%element_list%n_elements,2))
  !> compute minimum and maximum RZ jorek mesh
  R_minmax_jorek_sol = (/1.d21,-1.d21/); 
  Z_minmax_jorek_sol = (/1.d21,-1.d21/);
  psi_minmax_global_jorek_sol = (/1.d21,-1.d21/);
  psi_minmax_list_jorek_sol(:,1) = 1.d21
  psi_minmax_list_jorek_sol(:,2) = -1.d21
  do ii=1,fields_linear_sol%element_list%n_elements
    call RZ_minmax(fields_linear_sol%node_list,&
    fields_linear_sol%element_list,ii,R_min,R_max,Z_min,Z_max)
    R_minmax_jorek_sol(1) = min(R_min,R_minmax_jorek_sol(1))
    R_minmax_jorek_sol(2) = max(R_max,R_minmax_jorek_sol(2))
    Z_minmax_jorek_sol(1) = min(Z_min,Z_minmax_jorek_sol(1))
    Z_minmax_jorek_sol(2) = max(Z_max,Z_minmax_jorek_sol(2))
    call psi_minmax(fields_linear_sol%node_list,fields_linear_sol%element_list,&
    ii,psi_minmax_list_jorek_sol(ii,1),psi_minmax_list_jorek_sol(ii,2))
    psi_minmax_global_jorek_sol(1) = min(psi_minmax_global_jorek_sol(1),&
    psi_minmax_list_jorek_sol(ii,1))
    psi_minmax_global_jorek_sol(2) = max(psi_minmax_global_jorek_sol(2),&
    psi_minmax_list_jorek_sol(ii,2))
  enddo
  R_box_jorek_sol = R_minmax_jorek_sol*R_box_fraction_sol 
  R_box_jorek_sol = 5.d-1*sum(R_minmax_jorek_sol) + &
  5.d-1*(R_box_jorek_sol(2)-R_box_jorek_sol(1))*(/-1.d0,1.d0/)
  Z_box_jorek_sol = Z_minmax_jorek_sol*Z_box_fraction_sol
  Z_box_jorek_sol = 5.d-1*sum(Z_minmax_jorek_sol) + &
  5.d-1*(Z_box_jorek_sol(2)-Z_box_jorek_sol(1))*(/-1.d0,1.d0/)
end subroutine setup_init_intervals

!> find minimum and maximum of the requested fields
subroutine setup_field_minmax()
  use mod_rng,           only: type_rng
  use mod_pcg32_rng,     only: pcg32_rng
  use mod_fields_minmax, only: field_minmax_monte_carlo
  implicit none
  !> variables
  class(type_rng),dimension(:),allocatable :: rngs
  integer :: rank,n_tasks,ifail
  !> initialisation
  rank = 0; n_tasks=1;
  !> allocate and compute field minmax
  allocate(field_minmax_list_sol(n_fields_sol,2,fields_linear_sol%element_list%n_elements))
  call field_minmax_monte_carlo(fields_linear_sol%node_list,fields_linear_sol%element_list,&
  n_fields_sol,field_ids_sol,n_trials_sol,phi_small_sol,pcg32_rng(),rngs,&
  field_minmax_list_sol,field_minmax_global_sol,rank,n_tasks,ifail)
end subroutine setup_field_minmax

!> tear-down test feature
subroutine teardown()
  implicit none
  !> variables
  integer :: ii
  !> clean particle group
  if(allocated(groups_sol)) deallocate(groups_sol)
  !> cleanup fields
  call fields_sol%deallocate_fields()
  if(allocated(psi_minmax_list_jorek_sol)) deallocate(psi_minmax_list_jorek_sol)
end subroutine teardown

!> cleanup field minamx
subroutine teardown_field_minmax()
  implicit none
  if(allocated(field_minmax_list_sol)) deallocate(field_minmax_list_sol)
end subroutine teardown_field_minmax

!> Tests ------------------------------------------------------
!> test routine used for finding relativistic particles
subroutine test_find_relativistic_kinetic_gc_groups()
  use mod_initialise_relativistic_particles, only: find_relativistic_kinetic_gc_groups 
  implicit none
  !> variables
  integer :: n_active_groups_2_test
  integer,dimension(n_groups_2) :: active_group_ids_2_test
  !> extract relativitic groups
  call find_relativistic_kinetic_gc_groups(n_groups_2,groups_2_sol,&
  n_active_groups_2_test,active_group_ids_2_test)
  !> checks
  call assert_equals(n_active_groups_2_test,n_active_groups_2_sol,&
  "Error find relativistic particle groups: N# active groups mismatch!")
  call assert_equals(active_group_ids_2_test,active_group_ids_2_sol,n_groups_2,&
  "Error find relativistic particle groups: active group ids mismatch!")
end subroutine test_find_relativistic_kinetic_gc_groups

!> test find R,Z min max
subroutine test_find_RZPhi_minmax_global()
  use mod_initialise_relativistic_particles, only: find_RZPhi_minmax_global
  implicit none
  real*8,dimension(2) :: R_minmax_test,Z_minmax_test,phi_minmax_test
  !> test
  call find_RZPhi_minmax_global(R_minmax_test,Z_minmax_test,&
  phi_minmax_test,fields_linear_sol)
  !> checks
  call assert_equals(R_minmax_test,R_minmax_jorek_sol,2,&
  "Error find RZPhi minmax global: R minmax mismatch!")
  call assert_equals(Z_minmax_test,Z_minmax_jorek_sol,2,&
  "Error find RZPhi minmax global: Z minmax mismatch!")
  call assert_equals(phi_minmax_test,phi_minmax,2,&
  "Error find RZPhi minmax global: phi minmax mismatch!")
end subroutine test_find_RZPhi_minmax_global

!> test find poloidal flux, angle and toroidal angle minmax
subroutine test_find_psithetaphi_minmax_global_list()
  use mod_initialise_relativistic_particles, only: find_psithetaphi_minmax_global_list
  implicit none
  real*8,dimension(2) :: psi_minmax_global_test,theta_minmax_test,phi_minmax_test
  real*8,dimension(fields_linear_sol%element_list%n_elements,2) :: psi_minmax_list_test
  !> compute maximum and minimum
  call find_psithetaphi_minmax_global_list(psi_minmax_global_test,&
  psi_minmax_list_test,theta_minmax_test,phi_minmax_test,fields_linear_sol)
  !> checks
  call assert_equals(psi_minmax_global_test,psi_minmax_global_jorek_sol,2,&
  "Error find psi,theta,phi minmax global: psi global minmax mismatch!")
  call assert_equals(psi_minmax_list_test,psi_minmax_list_jorek_sol,&
  fields_linear_sol%element_list%n_elements,2,&
  "Error find psi,theta,phi minmax global: psi list minmax mismatch!")
  call assert_equals(theta_minmax_test,theta_minmax,2,&
  "Error find psi,theta,phi minmax global: theta minmax mismatch!")
  call assert_equals(phi_minmax_test,phi_minmax,2,&
  "Error find psi,theta,phi minmax global: phi minmax mismatch!")
end subroutine test_find_psithetaphi_minmax_global_list

!> test initialise particles in R,Z,phi coordinates
subroutine test_sample_position_uniformly_cylinder()
  use mod_interp, only: interp_RZ
  use mod_initialise_relativistic_particles, only: sample_position_uniformly_cylinder
  implicit none
  !> variables
  integer :: ii,jj,it,maxit
  real*8 :: R_test,Z_test
  real*8,dimension(2) :: R2_box
  real*8,dimension(3) :: rands
  type(particle_kinetic_relativistic) :: p_test
  logical,dimension(:),allocatable    :: success
  real*8,dimension(:),allocatable     :: errors,zeros
  !> initialisation
  maxit=100; R2_box = R_box_jorek_sol*R_box_jorek_sol
  !> loop for initialising particles
  do ii=1,n_groups
    allocate(success(n_particles(ii))); allocate(errors(n_particles(ii)));
    allocate(zeros(n_particles(ii))); errors=1.d16; zeros=0.d0;
    do jj=1,n_particles(ii)
      p_test%i_elm=0; it=0
      do while((p_test%i_elm.le.0).and.(it.le.maxit))
        call random_number(rands)
        call sample_position_uniformly_cylinder(p_test,&
        fields_linear_sol%node_list,fields_linear_sol%element_list,&
        rands,R2_box,Z_box_jorek_sol,phi_small_sol)
        it = it+1
      enddo
      call test_position_in_RZPhi_box(p_test,fields_linear_sol,&
      R_box_jorek_sol,Z_box_jorek_sol,phi_small_sol,success(jj),errors(jj))
    enddo
    call assert_true(all(success),&
    "Error sample position uniformly cylinder: R,Z,Phi not in bound!")
    call assert_equals(errors,zeros,n_particles(ii),tol_interp_real8,&
    "Error sample position uniformly cylinder: R,Z mismatch!")
    deallocate(success); deallocate(errors); deallocate(zeros);
  enddo
end subroutine test_sample_position_uniformly_cylinder

!> test initialise particles within psi,theta,phi bounds
subroutine test_sample_position_uniformly_psi_theta_phi()
  use mod_interp, only: interp_PRZ
  use equil_info, only: ES
  use mod_initialise_relativistic_particles, only: sample_position_uniformly_psi_theta_phi
  implicit none
  !> variables
  logical                             :: ifail
  integer                             :: ii,jj,it,maxit
  real*8                              :: R_test,Z_test
  real*8,dimension(n_v)               :: psi_test
  real*8,dimension(3)                 :: rands
  real*8,dimension(3,2)               :: thetapsiphi_box
  type(particle_kinetic_relativistic) :: p_test
  logical,dimension(:),allocatable    :: success
  real*8,dimension(:),allocatable     :: errors,zeros
  !> initialisation
  maxit = 100000
  !> denormalise the psi_box and store everything in the thetapsiphi variable
  thetapsiphi_box(1,:) = theta_small_sol
  thetapsiphi_box(2,:) = ES%psi_axis + (ES%psi_bnd-ES%psi_axis)*psi_small_sol
  thetapsiphi_box(3,:) = phi_small_sol
  !> loop for computing the particle position and check the bound
  do ii=1,n_groups
    allocate(success(n_particles(ii))); allocate(errors(n_particles(ii)));
    allocate(zeros(n_particles(ii)));   zeros = 0.d0;
    do jj=1,n_particles(ii)
      p_test%i_elm = 0; it = 0; ifail = .true.;
      do while(ifail.and.(p_test%i_elm.le.0).and.(it.le.maxit))
        call random_number(rands)
        call sample_position_uniformly_psi_theta_phi(p_test,&
        fields_linear_sol%node_list,fields_linear_sol%element_list,rands,&
        thetapsiphi_box,psi_minmax_list_jorek_sol,(/ES%R_axis,ES%Z_axis/),ifail)
        it = it+1
      enddo
      !> check end errors
      call test_position_in_psi_phi_box(p_test,n_v,i_psi,fields_linear_sol,&
      ES,psi_small_sol,phi_small_sol,success(jj),errors(jj))
    enddo
    call assert_true(all(success),&
    "Error sample position uniformly psi-theta-phi: psi,phi not in bound!")
    call assert_equals(errors,zeros,n_particles(ii),tol_interp_real8,&
    "Error sample position uniformly psi-theta-phi: R,Z mismatch!")
    deallocate(success); deallocate(errors); deallocate(zeros);
  enddo
end subroutine test_sample_position_uniformly_psi_theta_phi

!> test the sampling w.r.t. fluid field profile via accept-reject method
subroutine test_sample_position_acceptreject_from_fluid_profiles()
  use mod_rng,                               only: type_rng
  use mod_rng,                               only: setup_shared_rngs
  use mod_pcg32_rng,                         only: pcg32_rng
  use mod_accept_reject_funct,               only: accept_lower_values_rand 
  use mod_particle_types,                    only: particle_kinetic_relativistic
  use mod_initialise_relativistic_particles, only: sample_position_acceptreject_from_fluid_profiles
  implicit none
  !> variables:
  class(type_rng),dimension(:),allocatable :: rngs
  type(particle_kinetic_relativistic)      :: p_test
  integer                                  :: ii,jj,n_rngs,rng_id
  logical,dimension(:),allocatable         :: success
  !> initialisation
  call setup_shared_rngs(4+n_fields_sol,pcg32_rng(),rngs); 
  n_rngs=size(rngs); rng_id=1;
  !> sampled the particle positions
  do ii=1,n_groups
     allocate(success(n_particles(ii)))
     do jj=1,n_particles(ii)
       !> generate a particle
       call sample_position_acceptreject_from_fluid_profiles(p_test,&
       fields_linear_sol,n_rngs,n_fields_sol,rng_id,field_ids_sol,&
       field_minmax_global_sol,phi_small_sol,rngs,accept_lower_values_rand)
       !> check if the particle is in the st element box
       call test_position_in_stelement_box(p_test,&
       fields_linear_sol%element_list%n_elements,phi_small_sol,success(jj))
     enddo
     call assert_true(all(success),&
     "Error sample position accept-reject from fluid: s,t,i_elm not in bound!")
     deallocate(success)
  enddo
  !> cleanup
  if(allocated(rngs)) deallocate(rngs)
end subroutine test_sample_position_acceptreject_from_fluid_profiles

!> test particle kinetic and gc relativistic psithetaphi and energypitchgyro intervals
subroutine test_init_p_gc_relativistic_psithetaphi_energypitchgyro()
  use constants, only: EL_CHG,ATOMIC_MASS_UNIT,SPEED_OF_LIGHT
  use mod_coordinate_transforms,             only: vector_cylindrical_to_cartesian
  use mod_pcg32_rng,                         only: pcg32_rng
  use mod_interp,                            only: interp_PRZ
  use equil_info,                            only: ES
  use mod_particle_types,                    only: particle_kinetic_relativistic
  use mod_particle_types,                    only: particle_gc_relativistic
  use mod_pusher_tools,                      only: get_orthonormals
  use mod_particle_common_test_tools,        only: EThetaChi_RE_lowbnd,EThetaChi_RE_uppbnd
  use mod_particle_common_test_tools,        only: q1_posneg_interval
  use mod_initialise_relativistic_particles, only: init_particle_kinetic_relativistic_to_zero
  use mod_initialise_relativistic_particles, only: init_particle_gc_relativistic_to_zero
  use mod_initialise_relativistic_particles, only: init_p_gc_relativistic_psithetaphi_energypitchgyro
  implicit none
  !> variables:
  integer :: ii,jj
  real*8 :: psi_test,B_norm,psi_2,U
  real*8,dimension(2)   :: energy_box,pitch_box,gyro_box,p_box
  real*8,dimension(3)   :: B,e2,e3,E
  real*8,dimension(3,2) :: thetapsiphi_norm_box
  real*8,dimension(:),allocatable  :: errors,zeros
  logical,dimension(:),allocatable :: success_pos,success_vel,success_q
  !> initialisation:
  thetapsiphi_norm_box(1,:) = theta_small_sol
  thetapsiphi_norm_box(2,:) = psi_small_sol
  thetapsiphi_norm_box(3,:) = phi_small_sol
  energy_box = (/EThetaChi_RE_lowbnd(1),EThetaChi_RE_uppbnd(1)/)
  pitch_box  = (/EThetaChi_RE_lowbnd(2),EThetaChi_RE_uppbnd(2)/)
  gyro_box   = (/EThetaChi_RE_lowbnd(3),EThetaChi_RE_uppbnd(3)/)
  !> initialise particle groups
  call init_p_gc_relativistic_psithetaphi_energypitchgyro(groups_sol,&
  fields_linear_sol,ES,time_sol,pcg32_rng(),thetapsiphi_norm_box,&
  energy_box,pitch_box,gyro_box,q1_posneg_interval)
  !> perform checks
  do ii=1,n_groups
    select type (p_list=>groups_sol(ii)%particles)
    type is (particle_kinetic_relativistic)
      p_box = (EL_CHG*(/EThetaChi_RE_lowbnd(1),EThetaChi_RE_uppbnd(1)/))/&
      (ATOMIC_MASS_UNIT*groups_sol(ii)%mass*SPEED_OF_LIGHT*SPEED_OF_LIGHT)
      p_box = (p_box+1.d0)*(p_box+1.0); p_box = sqrt(p_box-1.d0);
      p_box = SPEED_OF_LIGHT*groups_sol(ii)%mass*p_box;
    type is (particle_gc_relativistic)
      p_box = (EL_CHG*(/EThetaChi_RE_lowbnd(1),EThetaChi_RE_uppbnd(1)/))/&
      (ATOMIC_MASS_UNIT*groups_sol(ii)%mass*SPEED_OF_LIGHT*SPEED_OF_LIGHT)
      p_box = (p_box+1.d0)*(p_box+1.0); p_box = sqrt(p_box-1.d0);
      p_box = SPEED_OF_LIGHT*groups_sol(ii)%mass*p_box;
    end select
    allocate(success_pos(n_particles(ii))); success_pos = .false.;
    allocate(success_vel(n_particles(ii))); success_vel = .false.;
    allocate(success_q(n_particles(ii))); success_q = .false.;
    allocate(errors(n_particles(ii))); errors = 1.d21;
    allocate(zeros(n_particles(ii)));  zeros = 0.d0;
    do jj=1,n_particles(ii)
      if(groups_sol(ii)%particles(jj)%i_elm.gt.0) then
        !> check the position initialisation
        call test_position_in_psi_phi_box(groups_sol(ii)%particles(jj),&
        n_v,i_psi,fields_linear_sol,ES,psi_small_sol,phi_small_sol,&
        success_pos(jj),errors(jj))
        !> compute magnetic coordinate system
        call fields_linear_sol%calc_EBpsiU(time_sol,&
        groups_sol(ii)%particles(jj)%i_elm,groups_sol(ii)%particles(jj)%st,&
        groups_sol(ii)%particles(jj)%x(3),E,B,psi_2,U)
        B_norm = norm2(B); B = B/B_norm;
        B = vector_cylindrical_to_cartesian(groups_sol(ii)%particles(jj)%x(3),B)
        call get_orthonormals(B,e2,e3)
        !> check momentum and charge solutions
        select type(p_test=>groups_sol(ii)%particles(jj))
        type is (particle_kinetic_relativistic)
          call test_momentum_in_ppitchgyro_box_p(p_test,B,e2,e3,&
          p_box,pitch_box,gyro_box,success_vel(jj))
          success_q(jj) = (p_test%q.ge.q1_posneg_interval(1)).and.&
          (p_test%q.le.q1_posneg_interval(2))
          call init_particle_kinetic_relativistic_to_zero(p_test) !< cleanup
        type is (particle_gc_relativistic)
          call test_momentum_in_ppitch_box_gc(p_test,groups_sol(ii)%mass,&
          B_norm,p_box,pitch_box,success_vel(jj))
          success_q(jj) = (p_test%q.ge.q1_posneg_interval(1)).and.&
          (p_test%q.le.q1_posneg_interval(2))
          call init_particle_gc_relativistic_to_zero(p_test) !< cleanup
        end select
      endif
    enddo
    call assert_true(all(success_pos),&
    "Error initialise relativistic kinetic p gc psi/theta/phi - E/pitch/gyro: psi,phi not in bound!")
    call assert_true(all(success_vel),&
    "Error initialise relativistic kinetic p gc psi/theta/phi - E/pitch/gyro: E,pitch,gyro not in bound!")
    call assert_true(all(success_q),&
    "Error initialise relativistic kinetic p gc psi/theta/phi - E/pitch/gyro: charge not in bound!")
    call assert_equals(errors,zeros,n_particles(ii),tol_interp_real8,&
    "Error initialise relativistic kinetic p gc psi/theta/phi - E/pitch/gyro: R,Z mismatch!")
    deallocate(errors); deallocate(zeros); 
    deallocate(success_pos); deallocate(success_vel); deallocate(success_q);
  enddo
end subroutine test_init_p_gc_relativistic_psithetaphi_energypitchgyro

!> test the initialisation of particles in RZPhi energy/pitch/gyro boxes
subroutine test_init_p_gc_relativistic_RZPhi_energypitchgyro()
  use constants, only: EL_CHG,ATOMIC_MASS_UNIT,SPEED_OF_LIGHT
  use mod_coordinate_transforms,             only: vector_cylindrical_to_cartesian
  use mod_pcg32_rng,                         only: pcg32_rng
  use mod_interp,                            only: interp_PRZ
  use equil_info,                            only: ES
  use mod_particle_types,                    only: particle_kinetic_relativistic
  use mod_particle_types,                    only: particle_gc_relativistic
  use mod_pusher_tools,                      only: get_orthonormals
  use mod_particle_common_test_tools,        only: EThetaChi_RE_lowbnd,EThetaChi_RE_uppbnd
  use mod_particle_common_test_tools,        only: q1_posneg_interval
  use mod_initialise_relativistic_particles, only: init_particle_kinetic_relativistic_to_zero
  use mod_initialise_relativistic_particles, only: init_particle_gc_relativistic_to_zero
  use mod_initialise_relativistic_particles, only: init_p_gc_relativistic_RZPhi_energypitchgyro
  implicit none
  !> variables:
  integer :: ii,jj
  real*8 :: psi_test,B_norm,psi_2,U
  real*8,dimension(2)   :: R_box_test,Z_box_test,phi_box_test
  real*8,dimension(2)   :: energy_box,pitch_box,gyro_box,p_box
  real*8,dimension(3)   :: B,e2,e3,E
  real*8,dimension(:),allocatable  :: errors,zeros
  logical,dimension(:),allocatable :: success_pos,success_vel,success_q
  !> initialisation:
  R_box_test=R_box_jorek_sol; Z_box_test=Z_box_jorek_sol; phi_box_test=phi_small_sol;
  energy_box = (/EThetaChi_RE_lowbnd(1),EThetaChi_RE_uppbnd(1)/)
  pitch_box  = (/EThetaChi_RE_lowbnd(2),EThetaChi_RE_uppbnd(2)/)
  gyro_box   = (/EThetaChi_RE_lowbnd(3),EThetaChi_RE_uppbnd(3)/)
  !> initialise particle groups
  call init_p_gc_relativistic_RZPhi_energypitchgyro(groups_sol,&
  fields_linear_sol,time_sol,pcg32_rng(),R_box_test,Z_box_test,&
  phi_box_test,energy_box,pitch_box,gyro_box,q1_posneg_interval) 
  !> perform checks
  do ii=1,n_groups
    select type (p_list=>groups_sol(ii)%particles)
    type is (particle_kinetic_relativistic)
      p_box = (EL_CHG*(/EThetaChi_RE_lowbnd(1),EThetaChi_RE_uppbnd(1)/))/&
      (ATOMIC_MASS_UNIT*groups_sol(ii)%mass*SPEED_OF_LIGHT*SPEED_OF_LIGHT)
      p_box = (p_box+1.d0)*(p_box+1.0); p_box = sqrt(p_box-1.d0);
      p_box = SPEED_OF_LIGHT*groups_sol(ii)%mass*p_box;
    type is (particle_gc_relativistic)
      p_box = (EL_CHG*(/EThetaChi_RE_lowbnd(1),EThetaChi_RE_uppbnd(1)/))/&
      (ATOMIC_MASS_UNIT*groups_sol(ii)%mass*SPEED_OF_LIGHT*SPEED_OF_LIGHT)
      p_box = (p_box+1.d0)*(p_box+1.0); p_box = sqrt(p_box-1.d0);
      p_box = SPEED_OF_LIGHT*groups_sol(ii)%mass*p_box;
    end select
    allocate(success_pos(n_particles(ii))); success_pos = .false.;
    allocate(success_vel(n_particles(ii))); success_vel = .false.;
    allocate(success_q(n_particles(ii))); success_q = .false.;
    allocate(errors(n_particles(ii))); errors = 1.d21;
    allocate(zeros(n_particles(ii)));  zeros = 0.d0;
    do jj=1,n_particles(ii)
      if(groups_sol(ii)%particles(jj)%i_elm.gt.0) then
        !> test position
        call test_position_in_RZPhi_box(groups_sol(ii)%particles(jj),&
        fields_linear_sol,R_box_jorek_sol,Z_box_jorek_sol,phi_small_sol,&
        success_pos(jj),errors(jj))
        !> compute magnetic coordinate system
        call fields_linear_sol%calc_EBpsiU(time_sol,&
        groups_sol(ii)%particles(jj)%i_elm,groups_sol(ii)%particles(jj)%st,&
        groups_sol(ii)%particles(jj)%x(3),E,B,psi_2,U)
        B_norm = norm2(B); B = B/B_norm;
        B = vector_cylindrical_to_cartesian(groups_sol(ii)%particles(jj)%x(3),B)
        call get_orthonormals(B,e2,e3)
        !> check momentum and charge solutions
        select type(p_test=>groups_sol(ii)%particles(jj))
        type is (particle_kinetic_relativistic)
          call test_momentum_in_ppitchgyro_box_p(p_test,B,e2,e3,&
          p_box,pitch_box,gyro_box,success_vel(jj))
          success_q(jj) = (p_test%q.ge.q1_posneg_interval(1)).and.&
          (p_test%q.le.q1_posneg_interval(2))
          call init_particle_kinetic_relativistic_to_zero(p_test) !< cleanup
        type is (particle_gc_relativistic)
          call test_momentum_in_ppitch_box_gc(p_test,groups_sol(ii)%mass,&
          B_norm,p_box,pitch_box,success_vel(jj))
          success_q(jj) = (p_test%q.ge.q1_posneg_interval(1)).and.&
          (p_test%q.le.q1_posneg_interval(2))
          call init_particle_gc_relativistic_to_zero(p_test) !< cleanup
        end select
      endif
    enddo
    call assert_true(all(success_pos),&
    "Error initialise relativistic kinetic p gc RZPhi - E/pitch/gyro: RZPhi not in bound!")
    call assert_true(all(success_vel),&
    "Error initialise relativistic kinetic p gc RZPhi - E/pitch/gyro: E,pitch,gyro not in bound!")
    call assert_true(all(success_q),&
    "Error initialise relativistic kinetic p gc RZPhi - E/pitch/gyro: charge not in bound!")
    call assert_equals(errors,zeros,n_particles(ii),tol_interp_real8,&
    "Error initialise relativistic kinetic p gc RZPhi - E/pitch/gyro: R,Z mismatch!")
    deallocate(errors); deallocate(zeros); 
    deallocate(success_pos); deallocate(success_vel); deallocate(success_q);
  enddo
end subroutine test_init_p_gc_relativistic_RZPhi_energypitchgyro

!> test the initialisaiton via accpet-reject sampling using
!> normalised fluid profiles as distribution
subroutine test_init_p_gc_relativistic_from_fluid_energypitchgyro()
  use constants, only: EL_CHG,ATOMIC_MASS_UNIT,SPEED_OF_LIGHT
  use mod_coordinate_transforms,             only: vector_cylindrical_to_cartesian
  use mod_pcg32_rng,                         only: pcg32_rng
  use mod_interp,                            only: interp_PRZ
  use mod_particle_types,                    only: particle_kinetic_relativistic
  use mod_particle_types,                    only: particle_gc_relativistic
  use mod_pusher_tools,                      only: get_orthonormals
  use mod_particle_common_test_tools,        only: EThetaChi_RE_lowbnd,EThetaChi_RE_uppbnd
  use mod_particle_common_test_tools,        only: q1_posneg_interval
  use mod_accept_reject_funct,               only: accept_lower_values_rand
  use mod_initialise_relativistic_particles, only: init_particle_kinetic_relativistic_to_zero
  use mod_initialise_relativistic_particles, only: init_particle_gc_relativistic_to_zero
  use mod_initialise_relativistic_particles, only: init_p_gc_relativistic_from_fluid_energypitchgyro
  implicit none
  !> variables:
  integer                          :: ii,jj,rank,n_tasks,ifail
  real*8                           :: psi_test,B_norm,psi_2,U
  real*8,dimension(2)              :: phi_box_test,p_box,energy_box
  real*8,dimension(2)              :: pitch_box,gyro_box
  real*8,dimension(3)              :: B,e2,e3,E
  logical,dimension(:),allocatable :: success_pos,success_vel,success_q
  !> initialisations
  rank = 0; n_tasks = 1; phi_box_test = phi_small_sol;
  energy_box = (/EThetaChi_RE_lowbnd(1),EThetaChi_RE_uppbnd(1)/)
  pitch_box  = (/EThetaChi_RE_lowbnd(2),EThetaChi_RE_uppbnd(2)/)
  gyro_box   = (/EThetaChi_RE_lowbnd(3),EThetaChi_RE_uppbnd(3)/)
  !> initialise particle groups
  call init_p_gc_relativistic_from_fluid_energypitchgyro(groups_sol,&
  fields_linear_sol,time_sol,pcg32_rng(),n_fields_sol,field_ids_sol,&
  accept_lower_values_rand,phi_box_test,energy_box,pitch_box,gyro_box,&
  q1_posneg_interval,rank,n_tasks,ifail,n_trials_sol)
  !> perform checks
  do ii=1,n_groups
    select type (p_list=>groups_sol(ii)%particles)
    type is (particle_kinetic_relativistic)
      p_box = (EL_CHG*(/EThetaChi_RE_lowbnd(1),EThetaChi_RE_uppbnd(1)/))/&
      (ATOMIC_MASS_UNIT*groups_sol(ii)%mass*SPEED_OF_LIGHT*SPEED_OF_LIGHT)
      p_box = (p_box+1.d0)*(p_box+1.0); p_box = sqrt(p_box-1.d0);
      p_box = SPEED_OF_LIGHT*groups_sol(ii)%mass*p_box;
    type is (particle_gc_relativistic)
      p_box = (EL_CHG*(/EThetaChi_RE_lowbnd(1),EThetaChi_RE_uppbnd(1)/))/&
      (ATOMIC_MASS_UNIT*groups_sol(ii)%mass*SPEED_OF_LIGHT*SPEED_OF_LIGHT)
      p_box = (p_box+1.d0)*(p_box+1.0); p_box = sqrt(p_box-1.d0);
      p_box = SPEED_OF_LIGHT*groups_sol(ii)%mass*p_box;
    end select
    allocate(success_pos(n_particles(ii))); success_pos = .false.;
    allocate(success_vel(n_particles(ii))); success_vel = .false.;
    allocate(success_q(n_particles(ii))); success_q = .false.;
    do jj=1,n_particles(ii)
      if(groups_sol(ii)%particles(jj)%i_elm.gt.0) then
        !> TODO: TEST SPATIAL DISTRIBUTION
        !> check if the particle is in the st element box
        call test_position_in_stelement_box(groups_sol(ii)%particles(jj),&
        fields_linear_sol%element_list%n_elements,phi_small_sol,success_pos(jj))
        !> compute magnetic coordinate system
        call fields_linear_sol%calc_EBpsiU(time_sol,&
        groups_sol(ii)%particles(jj)%i_elm,groups_sol(ii)%particles(jj)%st,&
        groups_sol(ii)%particles(jj)%x(3),E,B,psi_2,U)
        B_norm = norm2(B); B = B/B_norm;
        B = vector_cylindrical_to_cartesian(groups_sol(ii)%particles(jj)%x(3),B)
        call get_orthonormals(B,e2,e3)
        !> check momentum and charge solutions
        select type(p_test=>groups_sol(ii)%particles(jj))
        type is (particle_kinetic_relativistic)
          call test_momentum_in_ppitchgyro_box_p(p_test,B,e2,e3,&
          p_box,pitch_box,gyro_box,success_vel(jj))
          success_q(jj) = (p_test%q.ge.q1_posneg_interval(1)).and.&
          (p_test%q.le.q1_posneg_interval(2))
          call init_particle_kinetic_relativistic_to_zero(p_test) !< cleanup
        type is (particle_gc_relativistic)
          call test_momentum_in_ppitch_box_gc(p_test,groups_sol(ii)%mass,&
          B_norm,p_box,pitch_box,success_vel(jj))
          success_q(jj) = (p_test%q.ge.q1_posneg_interval(1)).and.&
          (p_test%q.le.q1_posneg_interval(2))
          call init_particle_gc_relativistic_to_zero(p_test) !< cleanup
        end select
      endif
    enddo
    call assert_true(all(success_pos),&
    "Error initialise relativistic kinetic p gc fluid prof. - E/pitch/gyro: s,t,i_elm not in bound!")
    call assert_true(all(success_vel),&
    "Error initialise relativistic kinetic p gc fluid prof. - E/pitch/gyro: E,pitch,gyro not in bound!")
    call assert_true(all(success_q),&
    "Error initialise relativistic kinetic p gc fluid prof. - E/pitch/gyro: charge not in bound!")
    deallocate(success_pos); deallocate(success_vel); deallocate(success_q);
  enddo
end subroutine test_init_p_gc_relativistic_from_fluid_energypitchgyro

!> test initialisation particle momentum between limits
subroutine test_sampling_cartesian_p_kinetic_relativistic()
  use mod_gnu_rng,                           only: gnu_rng_interval
  use mod_particle_common_test_tools,        only: vp3d_lowbnd,vp3d_uppbnd
  use mod_initialise_relativistic_particles, only: sampling_cartesian_p_kinetic_relativistic
  implicit none
  !> variables:
  integer :: ii,jj
  real*8,dimension(3)  :: rand
  real*8,dimension(3,2) :: pxpypz_int
  logical,dimension(:),allocatable :: success
  !> initialisation
  pxpypz_int(:,1) = vp3d_lowbnd; pxpypz_int(:,2) = vp3d_uppbnd;
  !> the the random kinetic particle initialisation 
  do ii=1,n_groups
    select type (p_list=>groups_sol(ii)%particles)
    type is (particle_kinetic_relativistic)
      allocate(success(n_particles(ii))); success = .false.;
      do jj=1,n_particles(ii)
        call random_number(rand)
        call sampling_cartesian_p_kinetic_relativistic(p_list(jj),rand,pxpypz_int)
        success(jj) = &
        ((p_list(jj)%p(1).ge.pxpypz_int(1,1)).and.(p_list(jj)%p(1).le.pxpypz_int(1,2))).and.&
        ((p_list(jj)%p(2).ge.pxpypz_int(2,1)).and.(p_list(jj)%p(2).le.pxpypz_int(2,2))).and.&
        ((p_list(jj)%p(3).ge.pxpypz_int(3,1)).and.(p_list(jj)%p(3).le.pxpypz_int(3,2)))
        p_list(jj)%p =0.d0
      enddo
      call assert_true(all(success),"Error sampling cart. p kinetic relat.: momenta not in bound!")
      deallocate(success)
    end select
  enddo
end subroutine test_sampling_cartesian_p_kinetic_relativistic

!> test initialisation gc momentum between limits
subroutine test_sampling_cartesian_gc_kinetic_relativistic()
  use mod_gnu_rng,                           only: gnu_rng_interval
  use mod_coordinate_transforms,             only: vector_cylindrical_to_cartesian
  use mod_particle_common_test_tools,        only: RZPhi_lowbnd,RZPhi_uppbnd
  use mod_particle_common_test_tools,        only: vp3d_lowbnd,vp3d_uppbnd
  use mod_initialise_relativistic_particles, only: sampling_cartesian_p_gc_relativistic
  implicit none
  !> variables:
  integer :: ii,jj
  real*8 :: psi,U,B_norm
  real*8,dimension(2)   :: p_gc_sol
  real*8,dimension(3)   :: rand,B,E,RZPhi,p_kin_sol
  real*8,dimension(3,2) :: pxpypz_int
  logical,dimension(:),allocatable :: success_ppar,success_mu
  !> initialisation
  pxpypz_int(:,1) = vp3d_lowbnd; pxpypz_int(:,2) = vp3d_uppbnd;
  call gnu_rng_interval(3,RZPhi_lowbnd,RZPhi_uppbnd,RZPhi)
  call fields_sol%calc_EBPsiU(time_sol,0,RZPhi(1:2),RZPhi(3),&
  E,B,psi,U); B_norm = norm2(B); B = B/B_norm; 
  B = vector_cylindrical_to_cartesian(RZPhi(3),B)
  !> test the sampling of gc from cartesian momentum
  do ii=1,n_groups
    select type (p_list=>groups_sol(ii)%particles)
    type is (particle_gc_relativistic)
      allocate(success_ppar(n_particles(ii))); allocate(success_mu(n_particles(ii)));
      !> loops on the particles
      do jj=1,n_particles(ii)
        !> initialising particles
        p_list(jj)%x=RZPhi; p_list(jj)%st=RZPhi(1:2); p_list(jj)%i_elm=0;
        !> compute momentum
        call random_number(rand)
        call sampling_cartesian_p_gc_relativistic(p_list(jj),fields_sol,&
        rand,time_sol,groups_sol(ii)%mass,pxpypz_int)
        !> compute new momentum
        p_kin_sol = pxpypz_int(:,1)+(pxpypz_int(:,2)-pxpypz_int(:,1))*rand
        p_gc_sol(1) = dot_product(p_kin_sol,B);
        p_kin_sol = p_kin_sol - p_gc_sol(1)*B
        p_gc_sol(2) = dot_product(p_kin_sol,p_kin_sol)/(2.d0*groups_sol(ii)%mass*B_norm)
        !> checks
        success_ppar(jj) = .true.; success_mu(jj) = .true.;
        if(p_gc_sol(1).ne.0.d0) then
          success_ppar(jj) = (abs((p_gc_sol(1)-p_list(jj)%p(1))/p_gc_sol(1))).le.tol_real8
        else
          if(p_list(jj)%p(1).ne.0.d0) success_ppar(1) = .false.
        endif
        if(p_gc_sol(2).ne.0.d0) then 
          success_mu(jj) = (abs((p_gc_sol(2)-p_list(jj)%p(2))/p_gc_sol(2))).le.tol_real8
        else
          if(p_list(jj)%p(2).ne.0.d0) success_mu(jj) = .false.
        endif
        !> cleaning particles
        p_list(jj)%x=0.d0; p_list(jj)%st=0.d0; p_list(jj)%p=0.d0;
      enddo
      call assert_true(all(success_ppar),&
      "Error sampling cart. gc kinetic relat.: parallel momentum not in bound!")
      call assert_true(all(success_mu),&
      "Error sampling cart. gc kinetic relat.: magnetic moment not in bound!")
      deallocate(success_ppar); deallocate(success_mu);
    end select
  enddo
end subroutine test_sampling_cartesian_gc_kinetic_relativistic

!> test initialisation relativistic kinetic particle momentum from
!> momentum intensity, pitch and gyro angles
subroutine test_sampling_uniform_ppitchgyro_kinetic_relativistic()
  use constants,                             only: PI,SPEED_OF_LIGHT,ATOMIC_MASS_UNIT,EL_CHG
  use mod_coordinate_transforms,             only: vector_cylindrical_to_cartesian
  use mod_gnu_rng,                           only: gnu_rng_interval
  use mod_pusher_tools,                      only: get_orthonormals
  use mod_particle_common_test_tools,        only: RZPhi_lowbnd,RZPhi_uppbnd
  use mod_particle_common_test_tools,        only: EThetaChi_RE_lowbnd,EThetaChi_RE_uppbnd
  use mod_initialise_relativistic_particles, only: sampling_uniform_ppitchgyro_kinetic_relativistic
  implicit none
  !> variabes:
  integer :: ii,jj
  real*8 :: psi,U,B_norm,p_norm_test,theta_test,gyro_test
  real*8,dimension(2) :: p_int,costheta_int,gyro_int
  real*8,dimension(3) :: rand,RZPhi,B,E,e2,e3
  real*8,dimension(3,2) :: EThetaChi
  logical,dimension(:),allocatable :: success
  !> initialisations
  costheta_int= (/cos(EThetaChi_RE_lowbnd(2)),cos(EThetaChi_RE_uppbnd(2))/); 
  gyro_int = (/EThetaChi_RE_lowbnd(3),EThetaChi_RE_uppbnd(3)/);
  call gnu_rng_interval(2,RZPhi_lowbnd,RZPhi_uppbnd,RZPhi)
  call fields_sol%calc_EBPsiU(time_sol,0,RZPhi(1:2),RZPhi(3),&
  E,B,psi,U); B_norm = norm2(B); B = B/B_norm; 
  B = vector_cylindrical_to_cartesian(RZPhi(3),B);
  !> compute orthogonal coordinate system
  call get_orthonormals(B,e2,e3)
  !> loop on the particles
  do ii=1,n_groups
    select type (p_list=>groups_sol(ii)%particles)
    type is (particle_kinetic_relativistic)
      p_int = (EL_CHG*(/EThetaChi_RE_lowbnd(1),EThetaChi_RE_uppbnd(1)/))/&
      (ATOMIC_MASS_UNIT*groups_sol(ii)%mass*SPEED_OF_LIGHT*SPEED_OF_LIGHT)
      p_int = (p_int+1.d0)*(p_int+1.0); p_int = sqrt(p_int-1.d0);
      p_int = SPEED_OF_LIGHT*groups_sol(ii)%mass*p_int;
      allocate(success(n_particles(ii)))
      !> loop on the particles
      do jj=1,n_particles(ii)
        !> initialising particles
        call random_number(rand)
        p_list(jj)%x=RZPhi; p_list(jj)%st=RZPhi(1:2); p_list(jj)%i_elm=0;
        call sampling_uniform_ppitchgyro_kinetic_relativistic(p_list(jj),&
        fields_sol,rand,time_sol,p_int**3.d0,costheta_int,gyro_int)
        !> check solotion
        call test_momentum_in_ppitchgyro_box_p(p_list(jj),B,e2,e3,&
        p_int,(/EThetaChi_RE_lowbnd(2),EThetaChi_RE_uppbnd(2)/),&
        (/EThetaChi_RE_lowbnd(3),EThetaChi_RE_uppbnd(3)/),success(jj))
        !> cleaning particles
        p_list(jj)%x=0.d0; p_list(jj)%st=0.d0; p_list(jj)%p=0.d0;
      enddo
      call assert_true(all(success),&
      "Error sampling particle uniform p, pitch, gyro kinetic relat.: momenta not in bound!")
      deallocate(success)
    end select
  enddo
end subroutine test_sampling_uniform_ppitchgyro_kinetic_relativistic

!> test the initialisation of relativistic kinetic gcs from relativistic kinetic
!> particles using spherical coordinates
subroutine test_sampling_uniform_ppitchgyro_gc_relativistic()
  use constants,                             only: PI,SPEED_OF_LIGHT,ATOMIC_MASS_UNIT,EL_CHG
  use mod_coordinate_transforms,             only: vector_cylindrical_to_cartesian
  use mod_gnu_rng,                           only: gnu_rng_interval
  use mod_particle_common_test_tools,        only: RZPhi_lowbnd,RZPhi_uppbnd
  use mod_particle_common_test_tools,        only: EThetaChi_RE_lowbnd,EThetaChi_RE_uppbnd
  use mod_initialise_relativistic_particles, only: sampling_uniform_ppitchgyro_gc_relativistic
  implicit none
  !> variabes:
  integer :: ii,jj
  real*8 :: psi,U,B_norm,p_norm_test,theta_test,theta_test_2
  real*8,dimension(2) :: p_int,costheta_int,gyro_int
  real*8,dimension(3) :: rand,RZPhi,B,E,e2,e3
  real*8,dimension(3,2) :: EThetaChi
  logical,dimension(:),allocatable :: success
  !> initialisation
  costheta_int= (/cos(EThetaChi_RE_lowbnd(2)),cos(EThetaChi_RE_uppbnd(2))/); 
  gyro_int = (/EThetaChi_RE_lowbnd(3),EThetaChi_RE_uppbnd(3)/);
  call gnu_rng_interval(2,RZPhi_lowbnd,RZPhi_uppbnd,RZPhi)
  call fields_sol%calc_EBPsiU(time_sol,0,RZPhi(1:2),RZPhi(3),&
  E,B,psi,U); B_norm = norm2(B); B = B/B_norm; 
  B = vector_cylindrical_to_cartesian(RZPhi(3),B);
  do ii=1,n_groups
    select type (p_list=>groups_sol(ii)%particles)
    type is (particle_gc_relativistic)
      p_int = (EL_CHG*(/EThetaChi_RE_lowbnd(1),EThetaChi_RE_uppbnd(1)/))/&
      (ATOMIC_MASS_UNIT*groups_sol(ii)%mass*SPEED_OF_LIGHT*SPEED_OF_LIGHT)
      p_int = (p_int+1.d0)*(p_int+1.0); p_int = sqrt(p_int-1.d0);
      p_int = SPEED_OF_LIGHT*groups_sol(ii)%mass*p_int;
      allocate(success(n_particles(ii)))
      do jj=1,n_particles(ii)
        !> initialising particles
        call random_number(rand)
        p_list(jj)%x=RZPhi; p_list(jj)%st=RZPhi(1:2); p_list(jj)%i_elm=0;
        call sampling_uniform_ppitchgyro_gc_relativistic(p_list(jj),fields_sol,&
        rand,groups_sol(ii)%mass,time_sol,p_int**3.d0,costheta_int,gyro_int)
        !> checks
        call test_momentum_in_ppitch_box_gc(p_list(jj),groups_sol(ii)%mass,&
        B_norm,p_int,(/EThetaChi_RE_lowbnd(2),EThetaChi_RE_uppbnd(2)/),success(jj))
      enddo
      call assert_true(all(success),&
      "Error sampling gc uniform p, pitch, gyro kinetic relat.: momenta not in bound!")
      deallocate(success)
    end select
  enddo
end subroutine test_sampling_uniform_ppitchgyro_gc_relativistic

!> test initialisation of particles to zero
subroutine test_particle_base_init_to_zero()
  use mod_initialise_relativistic_particles, only: init_particle_base_to_zero
  implicit none
  !> variables:
  integer :: ii,jj
  integer,dimension(:),allocatable  :: i_elm_sol,i_life_sol,i_elm_test,i_life_test
  real*4,dimension(:),allocatable   :: t_birth_sol,t_birth_test
  real*8,dimension(:),allocatable   :: weight_sol,weight_test
  real*8,dimension(:,:),allocatable :: x_sol,st_sol,x_test,st_test
  !> test initialisation to zero
  do ii=1,n_groups
    allocate(i_elm_sol(n_particles(ii)));   allocate(i_elm_test(n_particles(ii)));
    allocate(i_life_sol(n_particles(ii)));  allocate(i_life_test(n_particles(ii)));
    allocate(t_birth_sol(n_particles(ii))); allocate(t_birth_test(n_particles(ii)));
    allocate(weight_sol(n_particles(ii)));  allocate(weight_test(n_particles(ii)));
    allocate(x_sol(3,n_particles(ii)));     allocate(x_test(3,n_particles(ii)));
    allocate(st_sol(2,n_particles(ii)));      allocate(st_test(2,n_particles(ii)));
    i_elm_sol=0; i_life_sol=0; t_birth_sol=0.d0; weight_sol=1.d0; x_sol=0.d0; st_sol=0.d0;
    !> initialise particle to zero
    do jj=1,n_particles(ii)
      call init_particle_base_to_zero(groups_sol(ii)%particles(jj))
      i_elm_test(jj)   = groups_sol(ii)%particles(jj)%i_elm
      i_life_test(jj)  = groups_sol(ii)%particles(jj)%i_life
      t_birth_test(jj) = groups_sol(ii)%particles(jj)%t_birth
      weight_test(jj)  = groups_sol(ii)%particles(jj)%weight
      x_test(:,jj)     = groups_sol(ii)%particles(jj)%x
      st_test(:,jj)    = groups_sol(ii)%particles(jj)%st
    enddo
    !> checks
    call assert_equals(i_elm_test,i_elm_sol,n_particles(ii),&
    "initialise particle base to zero: mismatch!")
    call assert_equals(i_life_test,i_life_sol,n_particles(ii),&
    "initialise particle base to zero: mismatch!")
    call assert_equals(t_birth_test,t_birth_sol,n_particles(ii),&
    "initialise particle base to zero: mismatch!")
    call assert_equals(weight_test,weight_sol,n_particles(ii),&
    "initialise particle base to zero: mismatch!")
    call assert_equals(x_test,x_sol,3,n_particles(ii),&
    "initialise particle base to zero: x mismatch!")
    call assert_equals(st_test,st_sol,2,n_particles(ii),&
    "initialise particle base to zero: st mismatch!")
    !> cleanups
    deallocate(i_elm_sol);   deallocate(i_elm_test);
    deallocate(i_life_sol);  deallocate(i_life_test);
    deallocate(t_birth_sol); deallocate(t_birth_test);
    deallocate(weight_sol);  deallocate(weight_test);
    deallocate(x_sol);       deallocate(x_test);
    deallocate(st_sol);      deallocate(st_test);
  enddo
end subroutine test_particle_base_init_to_zero

!> test uniform random sampling
subroutine test_sampling_uniform_charge()
  use mod_particle_common_test_tools,        only: q1_pos_interval
  use mod_particle_common_test_tools,        only: q1_neg_interval
  use mod_particle_common_test_tools,        only: q1_posneg_interval
  use mod_initialise_relativistic_particles, only: sampling_uniform_charge
  implicit none
  !> variables
  integer                      :: ii,jj
  integer*1                    :: q1_test
  integer*1,dimension(2,5)     :: q1_intervals
  logical,dimension(n_samples) :: success
  real*8                       :: rand
  !> initialisation
  q1_intervals(:,1)=(/1,1/); q1_intervals(:,2)=(/-1,-1/);
  q1_intervals(:,3)=q1_pos_interval; q1_intervals(:,4)=q1_neg_interval;
  q1_intervals(:,5) = q1_posneg_interval;
  !> test charge sampling
  do ii=1,size(q1_intervals,2)
    do jj=1,n_samples
      call random_number(rand)
      q1_test = sampling_uniform_charge(rand,q1_intervals(:,ii))
      success(jj) = (q1_test.ge.q1_intervals(1,ii)).and.(q1_test.le.q1_intervals(2,ii))
    enddo
    !> checks
    call assert_true(all(success),"Error sampling uniform charge: charges not in interval!")
  enddo
end subroutine test_sampling_uniform_charge

!> test min max and renormalization of psi bound
subroutine test_check_thetapsiphi_interval()
  use mod_initialise_relativistic_particles, only: check_thetapsiphi_interval
  implicit none
  !> variables:
  real*8,dimension(2) :: psi_local_sol,psi_test,theta_test,phi_test
  !> check complete inflow
  psi_local_sol = psi_axisbnd(1)+(psi_axisbnd(2)-psi_axisbnd(1))*psi_small_sol
  psi_test=psi_small_sol; theta_test=theta_small_sol; phi_test=phi_small_sol;
  call check_thetapsiphi_interval(theta_test,psi_test,phi_test,psi_axisbnd(1),&
  psi_axisbnd(2),psi_minmax_sol,theta_minmax,phi_minmax)
  call assert_equals(psi_test,psi_local_sol,2,tol_real8,&
  "Error check psi-theta-phi interval: psi inflow test mismatch!")
  call assert_equals(theta_test,theta_small_sol,2,tol_real8,&
  "Error check psi-theta-phi interval: theta inflow test mismatch!")
  call assert_equals(phi_test,phi_small_sol,2,tol_real8,&
  "Error check psi-theta-phi interval: phi inflow test mismatch!")
  !> check complete overflow
  psi_test=psi_big_sol; theta_test=theta_big_sol; phi_test=phi_big_sol;
  call check_thetapsiphi_interval(theta_test,psi_test,phi_test,psi_axisbnd(1),&
  psi_axisbnd(2),psi_minmax_sol,theta_minmax,phi_minmax)
  call assert_equals(psi_test,psi_minmax_sol,2,tol_real8,&
  "Error check psi-theta-phi interval: psi overflow test mismatch!")
  call assert_equals(theta_test,theta_minmax,2,tol_real8,&
  "Error check psi-theta-phi interval: theta overflow test mismatch!")
  call assert_equals(phi_test,phi_minmax,2,tol_real8,&
  "Error check psi-theta-phi interval: phi overflow test mismatch!")
  !> checke underflow min and inflow max
  psi_test=(/psi_big_sol(1),psi_small_sol(2)/)
  theta_test=(/theta_big_sol(1),theta_small_sol(2)/)
  phi_test=(/phi_big_sol(1),phi_small_sol(2)/)
  call check_thetapsiphi_interval(theta_test,psi_test,phi_test,&
  psi_axisbnd(1),psi_axisbnd(2),psi_minmax_sol,theta_minmax,phi_minmax)
  call assert_equals(psi_test,(/psi_minmax_sol(1),psi_local_sol(2)/),2,&
  tol_real8,"Error check psi-theta-phi interval: psi min underflow test mismatch!")
  call assert_equals(theta_test,(/theta_minmax(1),theta_small_sol(2)/),2,&
  tol_real8,"Error check psi-theta-phi interval: theta min underflow test mismatch!")
  call assert_equals(phi_test,(/phi_minmax(1),phi_small_sol(2)/),2,&
  tol_real8,"Error check psi-theta-phi interval: phi min underflow test mismatch!")
  !> check inflow min and overflow max
  psi_test=(/psi_small_sol(1),psi_big_sol(2)/)
  theta_test=(/theta_small_sol(1),theta_big_sol(2)/)
  phi_test=(/phi_small_sol(1),phi_big_sol(2)/)
  call check_thetapsiphi_interval(theta_test,psi_test,phi_test,&
  psi_axisbnd(1),psi_axisbnd(2),psi_minmax_sol,theta_minmax,phi_minmax)
  call assert_equals(psi_test,(/psi_local_sol(1),psi_minmax_sol(2)/),2,&
  tol_real8,"Error check psi-theta-phi interval: psi max overflow test mismatch!")
  call assert_equals(theta_test,(/theta_small_sol(1),theta_minmax(2)/),2,&
  tol_real8,"Error check psi-theta-phi interval: theta max overflow test mismatch!")
  call assert_equals(phi_test,(/phi_small_sol(1),phi_minmax(2)/),2,&
  tol_real8,"Error check psi-theta-phi interval: phi max overflow test mismatch!")
end subroutine test_check_thetapsiphi_interval

!> test RZPhi min max bounding
subroutine test_check_RZPhi_interval()
  use mod_initialise_relativistic_particles, only: check_RZPhi_interval
  implicit none
  !> variables
  real*8,dimension(2) :: R_test,Z_test,phi_test
  !> check complete inflow
  R_test=R_small_sol;Z_test=Z_small_sol;phi_test=phi_small_sol;
  call check_RZPhi_interval(R_test,Z_test,phi_test,&
  R_minmax,Z_minmax,phi_minmax)
  call assert_equals(R_test,R_small_sol,2,&
  "Error check RZPhi interval: R inflow test mismatch!")
  call assert_equals(Z_test,Z_small_sol,2,&
  "Error check RZPhi interval: Z inflow test mismatch!")
  call assert_equals(phi_test,phi_small_sol,2,&
  "Error check RZPhi interval: phi inflow test mismatch!")
  !> check complete overflow
  R_test=R_big_sol;Z_test=Z_big_sol;phi_test=phi_big_sol;
  call check_RZPhi_interval(R_test,Z_test,phi_test,&
  R_minmax,Z_minmax,phi_minmax)
  call assert_equals(R_test,R_minmax,2,&
  "Error check RZPhi interval: R overflow test mismatch!")
  call assert_equals(Z_test,Z_minmax,2,&
  "Error check RZPhi interval: Z overflow test mismatch!")
  call assert_equals(phi_test,phi_minmax,2,&
  "Error check RZPhi interval: phi overflow test mismatch!")
  !> check min underflow and max inflow
  R_test=(/R_big_sol(1),R_small_sol(2)/)
  Z_test=(/Z_big_sol(1),Z_small_sol(2)/)
  phi_test=(/phi_big_sol(1),phi_small_sol(2)/)
  call check_RZPhi_interval(R_test,Z_test,phi_test,&
  R_minmax,Z_minmax,phi_minmax)
  call assert_equals(R_test,(/R_minmax(1),R_small_sol(2)/),2,&
  "Error check RZPhi interval: R underflow test mismatch!")
  call assert_equals(Z_test,(/Z_minmax(1),Z_small_sol(2)/),2,&
  "Error check RZPhi interval: Z underflow test mismatch!")
  call assert_equals(phi_test,(/phi_minmax(1),phi_small_sol(2)/),2,&
  "Error check RZPhi interval: phi underflow test mismatch!")
  !> check min inflow and max overflow
  R_test=(/R_small_sol(1),R_big_sol(2)/)
  Z_test=(/Z_small_sol(1),Z_big_sol(2)/)
  phi_test=(/phi_small_sol(1),phi_big_sol(2)/)
  call check_RZPhi_interval(R_test,Z_test,phi_test,&
  R_minmax,Z_minmax,phi_minmax)
  call assert_equals(R_test,(/R_small_sol(1),R_minmax(2)/),2,&
  "Error check RZPhi interval: R overflow test mismatch!")
  call assert_equals(Z_test,(/Z_small_sol(1),Z_minmax(2)/),2,&
  "Error check RZPhi interval: Z overflow test mismatch!")
  call assert_equals(phi_test,(/phi_small_sol(1),phi_minmax(2)/),2,&
  "Error check RZPhi interval: phi overflow test mismatch!")
end subroutine test_check_RZPhi_interval

!> test the momentum, pitch and gyro angles bounding box tests
subroutine test_check_energykinpitchgyro_interval()
  use constants, only: ATOMIC_MASS_UNIT,SPEED_OF_LIGHT,EL_CHG
  use mod_particle_common_test_tools, only: EThetaChi_RE_lowbnd,EThetaChi_RE_uppbnd
  use mod_initialise_relativistic_particles, only: check_energykinpitchgyro_interval
  implicit none
  !> variables:
  real*8 :: E0
  real*8,dimension(2) :: p_small_sol,e_test,p_test,pitch_test,gyro_test
  !> computes momentum for solution
  E0 = (groups_sol(1)%mass*ATOMIC_MASS_UNIT*&
  SPEED_OF_LIGHT*SPEED_OF_LIGHT)/EL_CHG
  p_small_sol = groups_sol(1)%mass*SPEED_OF_LIGHT*sqrt(&
  (/((EThetaChi_RE_lowbnd(1)/E0 + 1.d0)**2.d0)-1.d0,&
  ((EThetaChi_RE_uppbnd(1)/E0 + 1.d0)**2.d0)-1.d0/))
  !> check complete inflow
  e_test=(/EThetaChi_RE_lowbnd(1),EThetaChi_RE_uppbnd(1)/);
  pitch_test=pitch_small_sol;gyro_test=gyro_small_sol;
  call check_energykinpitchgyro_interval(e_test,p_test,&
  pitch_test,gyro_test,groups_sol(1)%mass)
  call assert_equals(abs((p_test-p_small_sol)/p_small_sol),(/0.d0,0.d0/),&
  2,tol_real8,"Error check p-pitch-gyro interval: p inflow test mismatch!")
  call assert_equals(pitch_test,theta_small_sol,2,&
  "Error check p-pitch-gyro interval: theta inflow test mismatch!")
  call assert_equals(gyro_test,gyro_small_sol,2,&
  "Error check p-pitch-gyro interval: gyro inflow test mismatch!")
  !> check complete overflow
  e_test=(/p_neg,p_neg/);pitch_test=pitch_big_sol;gyro_test=gyro_big_sol;
  call check_energykinpitchgyro_interval(e_test,p_test,&
  pitch_test,gyro_test,groups_sol(1)%mass)
  call assert_true(all(p_test.gt.0.d0),&
  "Error check p-pitch-gyro interval: overflow negative momentum found!")
  call assert_equals(pitch_test,(/0.d0,PI/),2,&
  "Error check p-pitch-gyro interval: pitch overflow test mismatch!")
  call assert_equals(gyro_test,(/0.d0,TWOPI/),2,&
  "Error check p-pitch-gyro interval: gyro overflow test mismatch!")
  !> check min underflow and max inflow
  e_test=(/p_neg,EThetaChi_RE_uppbnd(1)/)
  pitch_test=(/pitch_big_sol(1),pitch_small_sol(2)/)
  gyro_test=(/gyro_big_sol(1),gyro_small_sol(2)/)
  call check_energykinpitchgyro_interval(e_test,p_test,&
  pitch_test,gyro_test,groups_sol(1)%mass)
  call assert_true(all(p_test.gt.0.d0),&
  "Error check p-pitch-gyro interval: underflow negative momentum found!")
  call assert_equals(p_test(2),p_small_sol(2),&
  "Error check p-pitch-gyro interval: underflow momenum mismatch!")
  call assert_equals(pitch_test,(/0.d0,pitch_small_sol(2)/),2,&
  "Error check p-pitch-gyro interval: pitch underflow test mismatch!")
  call assert_equals(gyro_test,(/0.d0,gyro_small_sol(2)/),2,&
  "Error check p-pitch-gyro interval: gyro underflow test mismatch!")
  !> check min inflow and max overflow
  e_test=(/EThetaChi_RE_lowbnd(1),p_neg/)
  pitch_test=(/pitch_small_sol(1),pitch_big_sol(2)/)
  gyro_test=(/gyro_small_sol(1),gyro_big_sol(2)/)
  call check_energykinpitchgyro_interval(e_test,p_test,&
  pitch_test,gyro_test,groups_sol(1)%mass)
  call assert_true(all(p_test.gt.0.d0),&
  "Error check p-pitch-gyro interval: half overflow negative momentum found!")
  call assert_equals(abs((p_test(1)-p_small_sol(1))/p_small_sol(1)),0.d0,tol_real8,&
  "Error check p-pitch-gyro interval: half overflow momenum mismatch!")
  call assert_equals(pitch_test,(/pitch_small_sol(1),PI/),2,&
  "Error check p-pitch-gyro interval: half pitch overflow test mismatch!")
  call assert_equals(gyro_test,(/gyro_small_sol(1),phi_minmax(2)/),2,&
  "Error check p-pitch-gyro interval: half gyro overflow test mismatch!")
end subroutine test_check_energykinpitchgyro_interval

!> dummy test procedure
subroutine test_dummy()
  use mod_initialise_relativistic_particles
  use mod_fields_analytical, only: fields_analytical
  implicit none
  logical :: success
  success=.true.
  call assert_true(success,"Dummy test")
end subroutine test_dummy

!> Tools ------------------------------------------------------
!> test if a particle position is in the R,Z,Phi box
!> inputs:
!>   p_test:  (particle_base) test particle
!>   fields:  (fields_base) jorek fields
!>   R_box:   (real8)(2) major radius box
!>   Z_box:   (real8)(2) vertical position box
!>   phi_box: (real8)(2) toroidal angle box
!> outputs:
!>   success: (logical) true if particle in box
!>   error:   (real8) error in the R,Z interpolation
subroutine test_position_in_RZPhi_box(p_test,fields,R_box,&
Z_box,phi_box,success,error)
  use mod_interp,         only: interp_RZ
  use mod_fields,         only: fields_base
  use mod_particle_types, only: particle_base
  implicit none
  !> inputs:
  class(particle_base),intent(in) :: p_test
  class(fields_base),intent(in)   :: fields
  real*8,dimension(2),intent(in)  :: R_box,Z_box,phi_box
  !> outputs:
  logical,intent(out) :: success
  real*8,intent(out)  :: error
  !> variables:
  real*8 :: R_test,Z_test
  !> test correctness
  success = (p_test%x(1).ge.R_box(1)).and.(p_test%x(1).le.R_box(2)).and.&
  (p_test%x(2).ge.Z_box(1)).and.(p_test%x(2).le.Z_box(2)).and.&
  (p_test%x(3).ge.phi_box(1)).and.(p_test%x(3).le.phi_box(2))
  if(p_test%i_elm.gt.0) &
  call interp_RZ(fields%node_list,fields%element_list,p_test%i_elm,&
  p_test%st(1),p_test%st(2),R_test,Z_test)
  error = max(abs(p_test%x(1)-R_test),abs(p_test%x(2)-Z_test))
end subroutine test_position_in_RZPhi_box

!> test if a particle position is in the poloidal flux, 
!> toroidal angle box
!> inputs:
!>   p_test:  (particle_base) particle to test
!>   n_v_loc: (integer) number of varibales to interpolate
!>   i_v_loc: (integer)(n_v_loc) index of the variables to interpolate
!>   fields:  (fields_base) jorek fields
!>   eq_info: (t_equil_info) equilibrium object
!>   psi_box: (real8)(2) poloidal flux box
!>   phi_box: (real8)(2) toroidal angle box
!> outputs:
!>   success: (logical) true if in box
!>   error:   (real8) maximum error of the R,Z calculation
subroutine test_position_in_psi_phi_box(p_test,n_v_loc,i_v_loc,&
fields,eq_info,psi_box,phi_box,success,error)
  use mod_interp,         only: interp_PRZ
  use mod_particle_types, only: particle_base
  use mod_fields,         only: fields_base
  use equil_info,         only: t_equil_state
  implicit none
  !> inputs:
  class(particle_base),intent(in)       :: p_test
  class(fields_base),intent(in)         :: fields
  type(t_equil_state),intent(in)        :: eq_info
  integer,intent(in)                    :: n_v_loc
  integer,dimension(n_v_loc),intent(in) :: i_v_loc
  real*8,dimension(2),intent(in)        :: psi_box,phi_box
  !> outputs:
  logical,intent(out) :: success
  real*8,intent(out)  :: error
  !> variables
  real*8 :: R_test,Z_test
  real*8,dimension(n_v_loc) :: psi_test
  !> initialisation
  success = .false.; error = 1.d21;
  !> check end errors
  if(p_test%i_elm.gt.0) then
    call interp_PRZ(fields%node_list,fields%element_list,p_test%i_elm,&
    i_v_loc,n_v_loc,p_test%st(1),p_test%st(2),p_test%x(3),psi_test,R_test,Z_test)
    psi_test = (psi_test-eq_info%psi_axis)/(eq_info%psi_bnd-eq_info%psi_axis)
    success = ((psi_test(1).ge.psi_box(1)).and.(psi_test(1).le.psi_box(2))).and.&
    ((p_test%x(3).ge.phi_box(1)).and.(p_test%x(3).le.phi_box(2)))
    error = max(abs(p_test%x(1)-R_test),abs(p_test%x(2)-Z_test))
  endif
end subroutine test_position_in_psi_phi_box

!> test if a particle local coordinate system is well defined
!> inputs:
!>   p_test:     (particle_base) particle to test
!>   n_elements: (integer) maximum number of elements
!>   phi_box:    (real8)(2) toroidal angle interval
!> outputs:
!>   success: (logical) true if values in bound
subroutine test_position_in_stelement_box(p_test,n_elements,phi_box,success)
  use mod_particle_types, only: particle_base
  implicit none
  !> inputs:
  class(particle_base),intent(in) :: p_test
  integer,intent(in)              :: n_elements
  real*8,dimension(2),intent(in)  :: phi_box
  !> outputs:
  logical,intent(out) :: success
  !> check if particle values are in bound
  success = (all(p_test%st.ge.0.d0).and.all(p_test%st.le.1.d0)).and.&
            ((p_test%i_elm.ge.1).and.(p_test%i_elm.le.n_elements)).and.&
            ((p_test%x(3).ge.phi_box(1)).and.(p_test%x(3).le.phi_box(2)))
end subroutine test_position_in_stelement_box

!> test if the kinetic particle momentum is within the momentum,
!> pitch and gryo angle boxes
!> inputs:
!>   p_test: (particle_kinetic_relativistic) particle to test
!>   b:         (real8)(3) magnetic field direction
!>   e2:        (real8)(3) second direction
!>   e3:        (real8)(3) third direction
!>   p_box:     (real8)(2) normalised momentum box
!>   pitch_box: (real8)(2) pitch angle box
!>   gyro_box:  (real8)(2) gyro angle box
!> outputs:
!>   success: (logical) true if test passes
subroutine test_momentum_in_ppitchgyro_box_p(p_test,b,e2,e3,&
p_box,pitch_box,gyro_box,success)
  use constants,          only: PI
  use mod_particle_types, only: particle_kinetic_relativistic
  implicit none
  !> inputs:
  type(particle_kinetic_relativistic),intent(in) :: p_test
  real*8,dimension(3),intent(in)                 :: b,e2,e3
  real*8,dimension(2),intent(in)                 :: p_box,pitch_box,gyro_box
  !> outputs:
  logical,intent(out) :: success
  !> variables:
  real*8 :: p_norm_test,pitch_test,gyro_test
  !> check if momentum in box
  p_norm_test = norm2(p_test%p);
  pitch_test = acos(dot_product(p_test%p,b)/p_norm_test)
  gyro_test  = PI+atan2(dot_product(p_test%p,e3),dot_product(p_test%p,e2))
  success =  ((p_norm_test.ge.p_box(1)).and.(p_norm_test.le.p_box(2))).and.&
             ((pitch_test.ge.(pitch_box(1))).and.(pitch_test.le.pitch_box(2))).and.&
             ((gyro_test.ge.gyro_box(1)).and.(gyro_test.le.gyro_box(2)))
end subroutine test_momentum_in_ppitchgyro_box_p

!> test if the kinetic gc momentum is within the momentum,
!> and pitch boxes
!> inputs:
!>   p_test: (particle_kinetic_relativistic) particle to test
!>   mass:      (real8) gc mass
!>   B_norm:    (real8) magnetic field intensity
!>   p_box:     (real8)(2) normalised momentum box
!>   pitch_box: (real8)(2) pitch angle box
!> outputs:
!>   success: (logical) true if test passes
subroutine test_momentum_in_ppitch_box_gc(p_test,mass,&
B_norm,p_box,pitch_box,success)
  implicit none
  !> inputs:
  type(particle_gc_relativistic),intent(in) :: p_test
  real*8,intent(in)                         :: mass,B_norm
  real*8,dimension(2),intent(in)            :: p_box,pitch_box
  !> outputs:
  logical,intent(out) :: success
  !> variables:
  real*8 :: p_norm_test,theta_test,theta_test_2
  !> checks
  p_norm_test = sqrt(p_test%p(1)*p_test%p(1) + p_test%p(2)*2.d0*B_norm*mass) 
  theta_test = acos(p_test%p(1)/p_norm_test)
  theta_test_2 = asin(sqrt(p_test%p(2)*2.d0*B_norm*mass)/p_norm_test)
  success = ((p_norm_test.ge.p_box(1)).and.(p_norm_test.le.p_box(2))).and.&
            ((theta_test.ge.(pitch_box(1))).and.(theta_test.le.pitch_box(2))).and.&
            ((theta_test_2.ge.pitch_box(1)).and.(theta_test_2.le.pitch_box(2)))
end subroutine test_momentum_in_ppitch_box_gc
!>-------------------------------------------------------------
end module mod_initialise_relativistic_particles_test
