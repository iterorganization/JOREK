!> Module for the initialization of multiple relativistic
!> particles in the phase space (3D spatial 2D/3D 
!> velocity coordinates. Acceptance-rejection method
!> is used for the Monte-Carlo initialisation
module mod_initialise_relativistic_particles
implicit none
private
public :: accept_function
public :: init_p_gc_relativistic_psithetaphi_energypitchgyro
public :: init_p_gc_relativistic_RZPhi_energypitchgyro
public :: init_p_gc_relativistic_from_fluid_energypitchgyro
#ifdef UNIT_TESTS
public :: find_relativistic_kinetic_gc_groups
public :: sample_position_uniformly_cylinder
public :: sample_position_uniformly_psi_theta_phi
public :: sample_position_acceptreject_from_fluid_profiles
public :: sampling_cartesian_p_kinetic_relativistic
public :: sampling_cartesian_p_gc_relativistic
public :: sampling_uniform_ppitchgyro_kinetic_relativistic
public :: sampling_uniform_ppitchgyro_gc_relativistic
public :: sampling_uniform_charge
public :: init_particle_base_to_zero
public :: find_RZPhi_minmax_global
public :: check_RZPhi_interval
public :: find_psithetaphi_minmax_global_list
public :: check_thetapsiphi_interval
public :: check_energykinpitchgyro_interval
public :: init_particle_kinetic_relativistic_to_zero
public :: init_particle_gc_relativistic_to_zero
#endif

!> Variables and datatypes --------------------------------
!> Interfaces ---------------------------------------------
abstract interface
  !> function for acceptance rejection method it accepts 
  !> a list of critieria between [intervals(1),intervals(2)]
  !>and a random number between [0,1] and return true is 
  !> the value is accpeted
  !> inputs:
  !>   n_values:  (integer) number of values to test
  !>   values:    (real8)(n_values) values to compare
  !>   intervals: (real8)(n_values,2) minimum and maximum 
  !>              values for all intervals
  !>   rands:     (real8)(n_values) random numbers in [0,1]
  !> outputs:
  !>   success:    (logical) true if accepted
  function accept_function(n_values,values,intervals,rands) result(success)
    implicit none
    integer,intent(in)                      :: n_values
    real*8,dimension(n_values),intent(in)   :: rands,values
    real*8,dimension(n_values,2),intent(in) :: intervals
    logical                                 :: success
  end function accept_function
end interface

interface
  !> interface for field_minmax_monte_carlo
  subroutine field_minmax_monte_carlo(node_list,element_list,&
  n_fields,field_ids,n_trials,phi_int,rng_type,rngs,minmax_list,&
  minmax_global,my_id,n_tasks,ifail)
    use data_structure, only: type_node_list,type_element_list
    use mod_interp, only: interp_prz
    use mod_rng,    only: type_rng
    implicit none
    !> inputs-outputs:
    integer,intent(inout) :: ifail
    !> inputs:
    type(type_node_list),intent(in)                        :: node_list
    type(type_element_list),intent(in)                     :: element_list
    class(type_rng),dimension(:),allocatable,intent(inout) :: rngs
    class(type_rng),intent(in) :: rng_type
    integer,intent(in)                                     :: n_fields,n_trials
    integer,intent(in)                                     :: my_id,n_tasks
    integer,dimension(n_fields),intent(in)                 :: field_ids
    real*8,dimension(2),intent(in)                         :: phi_int
    !> outputs:
    real*8,dimension(n_fields,2),intent(out)               :: minmax_global
    real*8,dimension(n_fields,2,element_list%n_elements),intent(out) :: minmax_list
  end subroutine field_minmax_monte_carlo

  !> interface of the find RZ procedure
  subroutine find_RZ(node_list,element_list,R_find,Z_find,&
  R_out,Z_out,ielm_out,s_out,t_out,ifail)
    use data_structure
    use mod_element_rtree, only: elements_containing_point
    implicit none
    type(type_node_list),intent(in)    :: node_list
    type(type_element_list),intent(in) :: element_list
    real*8,intent(in)                  :: R_find,Z_find
    real*8,intent(out)                 :: R_out,Z_out,s_out,t_out
    integer,intent(inout)              :: ielm_out
    integer,intent(out)                :: ifail
  end subroutine find_RZ

  !> interface of the find theta psi procedure
  subroutine find_theta_psi(node_list,element_list,psi_minmax,&
  theta,psi,phi,R_axis,Z_axis,i_elm,s,t,R,Z,verbose_in)
    use data_structure
    implicit none
    !> inputs:
    type(type_node_list),intent(in)    :: node_list
    type(type_element_list),intent(in) :: element_list
    real*8,dimension(element_list%n_elements,2),intent(in) :: psi_minmax
    real*8,intent(in)                  :: theta,psi,phi
    real*8,intent(in)                  :: R_axis,Z_axis
    integer,intent(in)                 :: i_elm
    logical,intent(in),optional        :: verbose_in
    !> outputs:
    real*8,intent(out)                 :: s,t,R,Z 
  end subroutine find_theta_psi

  !> find minim and maximum value of the poloidal flux (psi)
  !> for each mesh element
  subroutine psi_minmax(node_list,element_list,i_elm,psimin,psimax)
    use data_structure
    implicit none
    type(type_node_list)    :: node_list
    type(type_element_list) :: element_list
    real*8                  :: psimin,psimax
    integer                 :: i_elm
  end subroutine psi_minmax

  !> find maximum and minimum values of the major radius and 
  !> vertical position
  subroutine RZ_minmax(node_list,element_list,i_elm,Rmin,Rmax,Zmin,Zmax)
    use data_structure
    implicit none
    type(type_node_list),intent(in)    :: node_list
    type(type_element_list),intent(in) :: element_list
    integer,intent(in)                 :: i_elm
    real*8,intent(out)                 :: Rmin,Rmax,Zmin,Zmax
  end subroutine RZ_minmax
end interface
contains
!> Procedures ---------------------------------------------
!> initialise kinetic particle relativistic in poloidal
!> flux, poloidal and toroidal angle while in the
!> velocity space the energy, pitch and gyro angles
!> coordinates are use
!> inputs:
!>   groups:          (particle_group)(:)(allocatable) group to initialise
!>   fields:          (fields_base) jorek fields base type
!>   eq_info:         (t_equil_state) information on the equilibrium
!>   rng_type:        (type_rng) type of random number generator
!>   eq_info:         (t_equil_info) information on the equilibrium
!>   time:            (real8) time of the interpolation
!>   thetapsiphi_box: (real8)(3,2) normalised (1,:) poloidal angle box,
!>                    (2,:) poloidal flux and (3,:) toroidal angle boxes
!>   energy_kin_box:  (real8)(2) kinetic energy box for sampling in eV
!>   pitch_box:       (real8)(2) pitch angle box for sampling 
!>   gyro_box:        (real8)(2) gyro angle box for sampling
!>   q_box:           (integer1)(2) minimum and maximum particle chatge
!> outputs:
!>   groups:          (particle_group)(:)(allocatable) initialised group
!>   fields:          (fields_base) jorek fields base type
!>   thetapsiphi_box: (real8)(3,2) normalised (1,:) poloidal angle box,
!>                    (2,:) poloidal flux and (3,:) toroidal angle boxes
!>   energy_box:      (real8)(2) energy box for sampling in eV
!>   pitch_box:       (real8)(2) pitch angle box for sampling 
!>   gyro_box:        (real8)(2) gyro angle box for sampling
subroutine init_p_gc_relativistic_psithetaphi_energypitchgyro(&
groups,fields,eq_info,time,rng_type,thetapsiphi_box,&
energy_kin_box,pitch_box,gyro_box,q_box)
  use mod_rng,            only: type_rng
  use mod_rng,            only: setup_shared_rngs
  use equil_info,         only: t_equil_state
  use mod_fields,         only: fields_base
  use mod_particle_sim,   only: particle_group
  use mod_particle_types, only: particle_kinetic_relativistic
  use mod_particle_types, only: particle_gc_relativistic
  !$ use omp_lib
  implicit none
  !> inputs-outputs:
  type(particle_group),dimension(:),allocatable,intent(inout) :: groups
  class(fields_base),intent(inout)                             :: fields
  real*8,dimension(3,2),intent(inout) :: thetapsiphi_box
  real*8,dimension(2),intent(inout)   :: energy_kin_box,pitch_box,gyro_box
  !> inputs:
  type(t_equil_state),intent(in)    :: eq_info
  class(type_rng),intent(in)        :: rng_type
  integer*1,dimension(2),intent(in) :: q_box
  real*8,intent(in)                 :: time
  !> variables:
  class(type_rng),dimension(:),allocatable :: rngs
  logical :: ifail
  integer :: ii,jj,kk,n_groups,thread_id,maxit,n_active_groups
  integer,dimension(:),allocatable :: active_group_ids
  real*8,dimension(2)   :: cospitch_box,psi_minmax_global
  real*8,dimension(2)   :: theta_minmax,phi_minmax
  real*8,dimension(7)   :: rands
  real*8,dimension(:,:),allocatable :: momentum_box
  real*8,dimension(fields%element_list%n_elements,2) :: psi_minmax_list

  !> initialisation
  maxit=1000000; n_groups=size(groups); n_active_groups=0; 
  allocate(active_group_ids(n_groups)); allocate(momentum_box(2,n_groups));
  call find_relativistic_kinetic_gc_groups(n_groups,groups,&
  n_active_groups,active_group_ids)
  !> extract bounding boxes
  call find_psithetaphi_minmax_global_list(psi_minmax_global,&
  psi_minmax_list,theta_minmax,phi_minmax,fields)
  !> check bounding boxes in physical and velocity spaces
  call check_thetapsiphi_interval(thetapsiphi_box(1,:),&
  thetapsiphi_box(2,:),thetapsiphi_box(3,:),eq_info%Psi_axis,&
  eq_info%Psi_bnd,psi_minmax_global,theta_minmax,phi_minmax)
  do ii=1,n_active_groups
    call check_energykinpitchgyro_interval(energy_kin_box,momentum_box(:,ii),&
    pitch_box,gyro_box,groups(active_group_ids(ii))%mass)
  enddo
  !> compute the cube of the momentum and the cosinus of the pitch angle
  momentum_box = momentum_box**3; cospitch_box = cos(pitch_box);
  !> initialise random number generator
  call setup_shared_rngs(7,rng_type,rngs) 
  
  !$omp parallel default(private) firstprivate(n_active_groups,&
  !$omp active_group_ids,maxit) shared(groups,rngs,fields,&
  !$omp thetapsiphi_box,psi_minmax_list,eq_info,time,momentum_box,&
  !$omp cospitch_box,gyro_box,q_box)
  !> initialise the the relativistic and particle lists
  thread_id = 1;
  !$ thread_id = omp_get_thread_num()+1
  do ii=1,n_active_groups
    !$omp do
    do jj=1,size(groups(active_group_ids(ii))%particles)
      !> initialise particle to 0
      call init_particle_base_to_zero(groups(active_group_ids(ii))%particles(jj))
      !> while loop until a valid element is not found
      kk = 0; ifail = .true.;
      do while(ifail.and.(kk.le.maxit))
        call rngs(thread_id)%next(rands)
        call sample_position_uniformly_psi_theta_phi(&
             groups(active_group_ids(ii))%particles(jj),fields%node_list,&
             fields%element_list,rands(1:3),thetapsiphi_box,&
             psi_minmax_list,(/eq_info%R_axis,eq_info%Z_axis/),ifail)
          kk = kk + 1
        enddo
        !> if valid, initialise both velocity space and charge
        select type (particle=>groups(active_group_ids(ii))%particles(jj))
        type is (particle_kinetic_relativistic)
          call sampling_uniform_ppitchgyro_kinetic_relativistic(particle,&
          fields,rands(4:6),time,momentum_box,cospitch_box,gyro_box)
          particle%q = sampling_uniform_charge(rands(7),q_box)
        type is (particle_gc_relativistic)
          call sampling_uniform_ppitchgyro_gc_relativistic(particle,fields,&
          rands(4:6),groups(active_group_ids(ii))%mass,time,&
          momentum_box,cospitch_box,gyro_box)
          particle%q = sampling_uniform_charge(rands(7),q_box)
        end select 
      enddo
      !$omp end do
  enddo
  !$omp end parallel

  !> cleanup
  deallocate(active_group_ids); deallocate(momentum_box);
end subroutine init_p_gc_relativistic_psithetaphi_energypitchgyro

!> initialise kinetic particle relativistic in major radius
!> verticla coordinate and toroidal angle while in the
!> velocity space the energy, pitch and gyro angles
!> coordinates are use
!> inputs:
!>   groups:         (particle_group)(:)(allocatable) group to initialise
!>   fields:         (fields_base) jorek fields base type
!>   eq_info:        (t_equil_state) information on the equilibrium
!>   rng_type:       (type_rng) type of random number generator
!>   time:           (real8) time of the interpolation
!>   R_box:          (real8)(2) major radius box
!>   Z_box:          (real8)(2) vertical coordinate box
!>   phi_box         (real8)(2) toroidal angles box
!>   energy_kin_box: (real8)(2) kinetic energy box for sampling in eV
!>   pitch_box:      (real8)(2) pitch angle box for sampling 
!>   gyro_box:       (real8)(2) gyro angle box for sampling
!>   q_box:          (integer1)(2) minimum and maximum particle chatge
!> outputs:
!>   groups:         (particle_group)(:)(allocatable) initialised group
!>   fields:         (fields_base) jorek fields base type
!>   R_box:          (real8)(2) major radius box
!>   Z_box:          (real8)(2) vertical coordinate box
!>   phi_box         (real8)(2) toroidal angles box
!>   energy_box:     (real8)(2) energy box for sampling in eV
!>   pitch_box:      (real8)(2) pitch angle box for sampling 
!>   gyro_box:       (real8)(2) gyro angle box for sampling
subroutine init_p_gc_relativistic_RZPhi_energypitchgyro(&
groups,fields,time,rng_type,R_box,Z_box,phi_box,&
energy_kin_box,pitch_box,gyro_box,q_box)
  use mod_rng,            only: type_rng
  use mod_rng,            only: setup_shared_rngs
  use mod_fields,         only: fields_base
  use mod_particle_sim,   only: particle_group
  use mod_particle_types, only: particle_kinetic_relativistic
  use mod_particle_types, only: particle_gc_relativistic
  !$ use omp_lib
  implicit none
  !> inputs-outputs:
  type(particle_group),dimension(:),allocatable,intent(inout) :: groups
  class(fields_base),intent(inout)                             :: fields
  real*8,dimension(2),intent(inout) :: R_box,Z_box,phi_box
  real*8,dimension(2),intent(inout) :: energy_kin_box,pitch_box,gyro_box
  !> inputs:
  class(type_rng),intent(in)        :: rng_type
  integer*1,dimension(2),intent(in) :: q_box
  real*8,intent(in)                 :: time
  !> variables:
  class(type_rng),dimension(:),allocatable :: rngs
  integer :: ii,jj,kk,n_groups,thread_id,maxit,n_active_groups
  integer,dimension(:),allocatable :: active_group_ids 
  real*8,dimension(2)   :: cospitch_box,R_minmax,R2_box
  real*8,dimension(2)   :: Z_minmax,phi_minmax
  real*8,dimension(7)   :: rands
  real*8,dimension(:,:),allocatable :: momentum_box

  !> initialisation
  maxit=1000000; n_groups=size(groups); n_active_groups=0; 
  allocate(active_group_ids(n_groups)); allocate(momentum_box(2,n_groups));
  call find_relativistic_kinetic_gc_groups(n_groups,groups,&
  n_active_groups,active_group_ids)
  !> extract bounding maximum and minimum boxes
  call find_RZPhi_minmax_global(R_minmax,Z_minmax,phi_minmax,fields)
  !> check bounding boxes
  call check_RZPhi_interval(R_box,Z_box,phi_box,R_minmax,Z_minmax,phi_minmax)
  do ii=1,n_active_groups
    call check_energykinpitchgyro_interval(energy_kin_box,momentum_box(:,ii),&
    pitch_box,gyro_box,groups(active_group_ids(ii))%mass)
  enddo
  !> compute the cube of the momentum and the cosinus of the pitch angle
  R2_box=R_box**2; momentum_box=momentum_box**3; cospitch_box=cos(pitch_box);
  !> initialise random number generator
  call setup_shared_rngs(7,rng_type,rngs)

  !> fill particle list
  !$omp parallel default(private) firstprivate(n_active_groups,&
  !$omp maxit,active_group_ids) shared(groups,fields,rngs,R2_box,&
  !$omp Z_box,phi_box,time,momentum_box,cospitch_box,gyro_box,q_box)
  thread_id = 1;
  !$ thread_id = omp_get_thread_num()+1
  do ii=1,n_active_groups
    !$omp do
    do jj=1,size(groups(active_group_ids(ii))%particles)
      !> initialise particle to 0
      call init_particle_base_to_zero(groups(active_group_ids(ii))%particles(jj))
      !> while loop until a valid element is not found
      kk = 0;
      do while((groups(active_group_ids(ii))%particles(jj)%i_elm.le.0).and.(kk.le.maxit))
        call rngs(thread_id)%next(rands)
        call sample_position_uniformly_cylinder(&
        groups(active_group_ids(ii))%particles(jj),&
        fields%node_list,fields%element_list,rands(1:3),&
        R2_box,Z_box,phi_box)
        kk = kk + 1
      enddo
      !> if valid, initialise both velocity space and charge
      if(groups(active_group_ids(ii))%particles(jj)%i_elm.gt.0) then
        select type (particle=>groups(active_group_ids(ii))%particles(jj))
        type is (particle_kinetic_relativistic)
          call sampling_uniform_ppitchgyro_kinetic_relativistic(particle,&
          fields,rands(4:6),time,momentum_box,cospitch_box,gyro_box)
          particle%q = sampling_uniform_charge(rands(7),q_box)
        type is (particle_gc_relativistic)
          call sampling_uniform_ppitchgyro_gc_relativistic(particle,fields,&
          rands(4:6),groups(active_group_ids(ii))%mass,time,&
          momentum_box,cospitch_box,gyro_box)
          particle%q = sampling_uniform_charge(rands(7),q_box)
        end select 
      endif
    enddo
    !$omp end do
  enddo
  !$omp end parallel

  !> cleanup
  deallocate(active_group_ids); deallocate(momentum_box);
end subroutine init_p_gc_relativistic_RZPhi_energypitchgyro

!> method for spacially initialising particles as a function
!> of the fluid fields using the accept-rejection method.
!> initialisation in velocity space is uniform on the sphere
!> defined by the particle energy, pitch and gyro angles
!> TODO: further optimise the initialisation so that 
!>       particles are initialisated sequentially in
!>       w.r.t. the mesh element and not randomly
!>       for increasing data coherency
!> inputs:
!>   groups:              (particle_group) particle groups
!>   fields:              (fields_base) jorek fields
!>   time:                (real8) simulation time
!>   rng_type:            (type_rng) type of the RNG to use
!>   n_profiles:          (integer) number of fluid fields to use
!>   prof_ids:            (integer)(n_fields) node indices of the
!>                        fluid fields to use
!>   accept:              (accept_function) function for accepting
!>                        a Monte-Carlo solution
!>   phi_box:             (real8)(2) toroidal angle sampling box
!>   energy_kin_box:      (real8)(2) kinetic energy (eV) sampling box
!>   pitch_box:           (real8)(2) pitch angle sampling box
!>   gyro_box:            (real8)(2) gyro angle sampling box
!>   q_box:               (integer1)(2) charge sampling box
!>   my_id:               (integer) id of the MPI task
!>   n_tasks:             (integer) total number of MPI tasks
!>   ifail:               (integer) MPI error
!>   n_trials_in:         (integer)(optional) N# trials for finding
!>                        the minimum and maximum values of fluid fields
!>   n_min_prof_to_set:   (integer) number of the field minima to override
!>   n_max_prof_to_set:   (integer) number of the field maxima to override
!>   min_prof_to_set_ids: (integer)(n_min_prof_to_set) indices of
!>                        the field minima to override
!>   max_prof_to_set_ids: (integer)(n_min_prof_to_set) indices of
!>                        the field maxima to override
!>   min_prof_to_set_val: (real8)(n_min_prof_to_set) values of the field
!>                        minima to override
!>   min_prof_to_set_val: (real8)(n_min_prof_to_set) values of the field
!>                        maxima to override
!> outputs:
!>   groups:         (particle_group) particle groups
!>   fields:         (fields_base) jorek fields
!>   accept:         (accept_function) function for accepting
!>                   a Monte-Carlo solution
!>   phi_box:        (real8)(2) toroidal angle sampling box
!>   energy_kin_box: (real8)(2) kinetic energy (eV) sampling box
!>   pitch_box:      (real8)(2) pitch angle sampling box
!>   gyro_box:       (real8)(2) gyro angle sampling box
!>   q_box:          (integer1)(2) charge sampling box
!>   ifail:          (integer) MPI error
subroutine init_p_gc_relativistic_from_fluid_energypitchgyro(&
groups,fields,time,rng_type,n_profiles,prof_ids,accept,&
R_box,Z_box,phi_box,energy_kin_box,pitch_box,gyro_box,&
q_box,my_id,n_tasks,ifail,n_trials_in,n_min_prof_to_set,&
n_max_prof_to_set,min_prof_to_set_ids,max_prof_to_set_ids,&
min_prof_to_set_val,max_prof_to_set_val)
  use constants,          only: TWOPI
  use mod_array_handlers, only: set_values_in_array
  use mod_rng,            only: type_rng
  use mod_rng,            only: setup_shared_rngs
  use mod_fields,         only: fields_base
  use mod_particle_sim,   only: particle_group
  use mod_particle_types, only: particle_kinetic_relativistic
  use mod_particle_types, only: particle_gc_relativistic
  use mod_fields_minmax,  only: field_minmax_monte_carlo
  !$ use omp_lib
  implicit none
  !> inputs-outputs:
  type(particle_group),dimension(:),allocatable,intent(inout) :: groups
  class(fields_base),intent(inout)                            :: fields
  procedure(accept_function)                                  :: accept
  integer,intent(inout)             :: ifail
  real*8,dimension(2),intent(inout) :: R_box,Z_box,phi_box
  real*8,dimension(2),intent(inout) :: energy_kin_box,pitch_box,gyro_box
  !> inputs:
  class(type_rng),intent(in)                :: rng_type
  integer*1,dimension(2),intent(in)         :: q_box
  integer,intent(in)                        :: n_profiles,my_id,n_tasks
  integer,intent(in),optional               :: n_trials_in
  integer,dimension(n_profiles),intent(in)  :: prof_ids 
  integer,intent(in),optional               :: n_min_prof_to_set,n_max_prof_to_set
  integer,dimension(:),intent(in),optional :: min_prof_to_set_ids
  integer,dimension(:),intent(in),optional :: max_prof_to_set_ids
  real*8,intent(in)                         :: time
  real*8,dimension(:),intent(in),optional :: min_prof_to_set_val
  real*8,dimension(:),intent(in),optional :: max_prof_to_set_val
  !> variables:
  class(type_rng),dimension(:),allocatable :: rngs
  integer :: ii,jj,kk,n_rngs,n_groups,thread_id,maxit,n_active_groups
  integer :: n_trials
  integer,dimension(:),allocatable  :: active_group_ids
  real*8,dimension(2)               :: R2_box,cospitch_box
  real*8,dimension(2)               :: R_minmax,Z_minmax,phi_minmax
  real*8,dimension(3+n_profiles)    :: rands
  real*8,dimension(n_profiles,2)    :: prof_norms
  real*8,dimension(:,:),allocatable :: momentum_box
  real*8,dimension(n_profiles,2,fields%element_list%n_elements) :: prof_norms_list
  !> initialisations
  n_trials=1000000; if(present(n_trials_in)) n_trials=n_trials_in;
  n_groups=size(groups); n_active_groups=0; 
  allocate(active_group_ids(n_groups)); allocate(momentum_box(2,n_groups));
  call find_relativistic_kinetic_gc_groups(n_groups,groups,&
  n_active_groups,active_group_ids)
  !> extract bounding maximum and minimum boxes
  call find_RZPhi_minmax_global(R_minmax,Z_minmax,phi_minmax,fields)
  !> check bounding boxes
  call check_RZPhi_interval(R_box,Z_box,phi_box,R_minmax,Z_minmax,phi_minmax)
  !> check bounding boxes
  call check_phi_interval(phi_box,(/0.d0,TWOPI/))
  do ii=1,n_active_groups
    call check_energykinpitchgyro_interval(energy_kin_box,momentum_box(:,ii),&
    pitch_box,gyro_box,groups(active_group_ids(ii))%mass)
  enddo
  !> compute the cube of the momentum and the cosinus of the pitch angle
  R2_box = R_box**2; momentum_box=momentum_box**3; cospitch_box=cos(pitch_box);
  !> estimate minimum and maximum values of the fluid fields
  call field_minmax_monte_carlo(fields%node_list,fields%element_list,&
  n_profiles,prof_ids,n_trials,phi_box,rng_type,rngs,prof_norms_list,&
  prof_norms,my_id,n_tasks,ifail)
  !> set minimum and maximum profiles if present
  if(present(n_min_prof_to_set).and.present(min_prof_to_set_ids).and.present(min_prof_to_set_val)) &
  call set_values_in_array(n_min_prof_to_set,n_profiles,min_prof_to_set_ids,min_prof_to_set_val,prof_norms(:,1))
  if(present(n_max_prof_to_set).and.present(max_prof_to_set_ids).and.present(max_prof_to_set_val)) &
  call set_values_in_array(n_max_prof_to_set,n_profiles,max_prof_to_set_ids,max_prof_to_set_val,prof_norms(:,2))
  !> initialise random number generator
  call setup_shared_rngs(4+n_profiles,rng_type,rngs)
  n_rngs = size(rngs)

  !> initialise particle groups
  !$omp parallel default(private) firstprivate(n_rngs,n_active_groups,&
  !$omp maxit,active_group_ids,n_profiles) shared(groups,fields,rngs,&
  !$omp prof_ids,prof_norms,R2_box,Z_box,phi_box,time,momentum_box,&
  !$omp cospitch_box,gyro_box,q_box)
  thread_id = 1
  !$ thread_id = omp_get_thread_num()+1
  do ii=1,n_active_groups
    !$omp do
    do jj=1,size(groups(active_group_ids(ii))%particles)
      !> initialise particle to 0
      call init_particle_base_to_zero(groups(active_group_ids(ii))%particles(jj))
      !> initialise particle in space from plasma profiles
      call sample_position_acceptreject_from_fluid_profiles(&
      groups(active_group_ids(ii))%particles(jj),fields,n_rngs,n_profiles,&
      thread_id,prof_ids,prof_norms,R2_box,Z_box,phi_box,rngs,accept)
      !> if valid, initialise both velocity space and charge
      if(groups(active_group_ids(ii))%particles(jj)%i_elm.gt.0) then
        call rngs(thread_id)%next(rands)
        select type (particle=>groups(active_group_ids(ii))%particles(jj))
        type is (particle_kinetic_relativistic)
          call sampling_uniform_ppitchgyro_kinetic_relativistic(particle,&
          fields,rands(1:3),time,momentum_box,cospitch_box,gyro_box)
          particle%q = sampling_uniform_charge(rands(4),q_box)
        type is (particle_gc_relativistic)
          call sampling_uniform_ppitchgyro_gc_relativistic(particle,fields,&
          rands(1:3),groups(active_group_ids(ii))%mass,time,&
          momentum_box,cospitch_box,gyro_box)
          particle%q = sampling_uniform_charge(rands(4),q_box)
        end select 
      endif
    enddo
    !$omp end do
  enddo
  !$omp end parallel
end subroutine init_p_gc_relativistic_from_fluid_energypitchgyro

!> find number and id of groups being reltivistic
!> kinetic or relativistic gc
!> inputs:
!>   n_groups: (integer) number of groups
!>   groups:   (particle_groups) particle groups
!> outputs:
!>   n_active_groups:  (integer) number of relativistic groups
!>   active_group_ids: (integer)(n_groups) id of relativitic group
subroutine find_relativistic_kinetic_gc_groups(n_groups,&
groups,n_active_groups,active_group_ids)
  use mod_particle_sim,   only: particle_group
  use mod_particle_types, only: particle_kinetic_relativistic
  use mod_particle_types, only: particle_gc_relativistic
  implicit none
  !> inputs:
  class(particle_group),dimension(n_groups),intent(in) :: groups
  integer,intent(in) :: n_groups
  !> outputs:
  integer,intent(out)                     :: n_active_groups
  integer,dimension(n_groups),intent(out) :: active_group_ids
  !> variables
  integer :: ii
  !> initialisation
  n_active_groups = 0; active_group_ids = 0;
  !> find relativistic particles
  do ii=1,n_groups
    select type (p_list=>groups(ii)%particles)
    type is (particle_kinetic_relativistic)
      n_active_groups = n_active_groups + 1
      active_group_ids(n_active_groups) = ii
    type is (particle_gc_relativistic)
      n_active_groups = n_active_groups + 1
      active_group_ids(n_active_groups) = ii
    end select
  enddo
end subroutine find_relativistic_kinetic_gc_groups

!> sample particles in the physical space. Particles
!> are distributed uniformely in the physical space
!> using the cylindrical coordinates.
!> inputs:
!>   particle:     (particle_base) particle to be sampled
!>   node_list:    (type_node_list) jorek element nodes
!>   element_list: (type_element_list) jorek element list
!>   rand:         (real8)(3) random numbers
!>   R2_bounds:    (real8)(2) squared of the major radius 
!>                            bounding box
!>   Z_bound:      (real8)(2) vertical position bounding box
!>   phi_bound:    (real8)(3) toroidal angle bounding box
!> outputs:
!>   particle: (particle_base) sampled particle
subroutine sample_position_uniformly_cylinder(particle,&
node_list,element_list,rand,R2_bound,Z_bound,phi_bound)
  use data_structure,     only: type_node_list
  use data_structure,     only: type_element_list
  use mod_particle_types, only: particle_base
  use mod_sampling,       only: transform_uniform_cylindrical
  implicit none
  !> inputs-outputs:
  class(particle_base),intent(inout) :: particle
  !> inputs:
  type(type_node_list),intent(in)    :: node_list
  type(type_element_list),intent(in) :: element_list
  real*8,dimension(2),intent(in)     :: R2_bound,Z_bound,phi_bound
  real*8,dimension(3),intent(in)     :: rand
  !> variables:
  integer :: ifail
  !> compute the particle position in global and local coords
  call transform_uniform_cylindrical(rand,R2_bound,Z_bound,&
  phi_bound,particle%x)
  call find_RZ(node_list,element_list,particle%x(1),particle%x(2),&
  particle%x(1),particle%x(2),particle%i_elm,particle%st(1),&
  particle%st(2),ifail)
end subroutine sample_position_uniformly_cylinder

!> sample particles in physical space. Particles are 
!> distributed uniformly in the poloidal flux and 
!> poloidal angle coordinates which implies that 
!> the distribution is not uniformly sampled in in 
!> physical space
!> inputs:
!>   particle:           (particle_base) particle to be sampled
!>   node_list:          (type_node_list) jorek element nodes
!>   element_list:       (type_element_list) jorek element list
!>   rand:               (real8)(3) random numbers
!>   thetapsiphi_bound:  (real8)(3,2) sampling box boundaries
!>                       :,1) minimum theta,psi,phi boundary
!>                       :,2) maximum theta,psi,phi boundary
!>   psi_element_minmax: (real8)(n_elements,2) minimum and 
!>                       maximum value of psi per each element
!>   RZ_axis:            (real8)(2) 1)-R 2)-Z magnetic axis coords.
!> outputs: 
!>   particle:           (particle_base) sampled particle
!>   ifail:              (logical) if 1 psi not in bound
subroutine sample_position_uniformly_psi_theta_phi(&
particle,node_list,element_list,rand,thetapsiphi_bound,&
psi_element_minmax,RZ_axis,ifail)
  use data_structure,     only: type_node_list
  use data_structure,     only: type_element_list
  use mod_interp,         only: interp_PRZ
  use mod_particle_types, only: particle_base
  implicit none
  !> inputs-outputs:
  class(particle_base),intent(inout) :: particle
  !> inputs:
  type(type_node_list),intent(in)    :: node_list
  type(type_element_list),intent(in) :: element_list
  real*8,dimension(2),intent(in)     :: RZ_axis
  real*8,dimension(3),intent(in)     :: rand
  real*8,dimension(3,2),intent(in)   :: thetapsiphi_bound
  real*8,dimension(element_list%n_elements,2),intent(in) :: psi_element_minmax
  !> outputs:
  logical,intent(out) :: ifail
  !> varibales
  real*8              :: R_test,Z_test
  real*8,dimension(1) :: psi_test
  !> initialisation
  ifail = .false.;
  !> compute particle position in global coordinates
  particle%x = thetapsiphi_bound(:,1)+(thetapsiphi_bound(:,2)-&
  thetapsiphi_bound(:,1))*rand
  call find_theta_psi(node_list,element_list,psi_element_minmax,&
  particle%x(1),particle%x(2),particle%x(3),RZ_axis(1),RZ_axis(2),&
  particle%i_elm,particle%st(1),particle%st(2),particle%x(1),&
  particle%x(2),.false.)
  !> check if psi is in bound
  call interp_PRZ(node_list,element_list,particle%i_elm,(/1/),1,&
  particle%st(1),particle%st(2),particle%x(3),psi_test,R_test,Z_test) 
  !> it does not work for psi box values of different sign!
  psi_test = (psi_test-thetapsiphi_bound(2,1))/(thetapsiphi_bound(2,2)-thetapsiphi_bound(2,1))
  if((psi_test(1).lt.0.d0).or.(psi_test(1).gt.1.d0)) ifail=.true.;
end subroutine sample_position_uniformly_psi_theta_phi

!> Sample position from fluid profiles using the accept-reject method.
!> Proper profile normalization must be provided so that the 
!> normalised profile is within [0,1]
!> TODO: further optimise the initialisation so that 
!>       particles are initialisated sequentially in
!>       w.r.t. the mesh element and not randomly
!>       for increasing data coherency
!> inputs:
!>   particle:   (particle_base) particle to be initialised
!>   fields:     (fields_base) jorek fields object
!>   n_rngs:     (integer) number of random number generator
!>   n_profiles: (integer) number of profiles to use
!>   rng_id:     (integer) id of the random number ot be used
!>   rngs:       (n_rngs) random number generators
!>   prof_ids:   (integer)(n_profiles) index of the profile in the node
!>   prof_norms: (real8)(n_profiles,2) profile normalization
!>               1: minimum value, 2: extension: maximum-minimum
!>   R2_box:     (real8)(2) box of the squared major radius
!>   Z_box:      (real8)(2) vertical coord. box
!>   phi_box:    (real8)(2) toroidal angle box
!>   accept:     (accept_function) accept function: takes as arguments
!>               a random number and a set of values within [0,1] and
!>               returns true if the variable is accepted
!> outputs:
!>   particle: (particle_base) initialised particle
!>   rngs:     (n_rngs) random number generators
!>   accept:   (accept_function) accept function: takes as arguments
!>             a random number and a set of values within [0,1] and
!>             returns true if the variable is accepted
subroutine sample_position_acceptreject_from_fluid_profiles(particle,&
fields,n_rngs,n_profiles,rng_id,prof_ids,prof_norms,R2_box,Z_box,&
phi_box,rngs,accept)
  use mod_interp,         only: interp_PRZ
  use mod_fields,         only: fields_base
  use mod_rng,            only: type_rng
  use mod_particle_types, only: particle_base
  implicit none
  !> inputs:
  class(fields_base),intent(in)                   :: fields
  integer,intent(in)                              :: n_rngs,n_profiles,rng_id
  integer,dimension(n_profiles),intent(in)        :: prof_ids
  real*8,dimension(2),intent(in)                  :: R2_box,Z_box,phi_box
  real*8,dimension(n_profiles,2),intent(in)       :: prof_norms
  !> inputs-outputs:
  class(particle_base),intent(inout)              :: particle
  class(type_rng),dimension(n_rngs),intent(inout) :: rngs
  procedure(accept_function)                      :: accept
  !> varibales
  logical                        :: fail
  integer                        :: maxit,it
  real*8,dimension(3+n_profiles) :: rands
  real*8,dimension(n_profiles)   :: profiles
  !> initialisations
  maxit = 100000; fail=.true.; it = -1;

  !> while loop until a particle is accepted
  do while(fail.and.(it.le.maxit))
    it = it+1 !< update the counter
    call rngs(rng_id)%next(rands) !< generate random numbers
    !> sampling the particle position uniformly in the cylinder
    call sample_position_uniformly_cylinder(particle,fields%node_list,&
    fields%element_list,rands(1:3),R2_box,Z_box,phi_box)
    if(particle%i_elm.le.0) cycle
    !> interpolate the profiles
    call interp_PRZ(fields%node_list,fields%element_list,particle%i_elm,&
    prof_ids,n_profiles,particle%st(1),particle%st(2),particle%x(3),&
    profiles,particle%x(1),particle%x(2))
    !> check for failures
    fail = .not.accept(n_profiles,profiles,prof_norms,rands(4:n_profiles+3))
  enddo
  !> if the method failed just set the particle element equal to 0
  if(fail) particle%i_elm = 0
end subroutine sample_position_acceptreject_from_fluid_profiles

!> uniform sampling of the cartesian momentum
!> for kinetic relativistic particles
!> inputs:
!>   particle:        (particle_kinetic_relativistic) particle to be sampled
!>   rand:            (real8)(2) random numbers
!>   pxpypz_interval: (real8)(3,2) momenta interval (in units of AMU)
!>                    1- px,py,pz momentum lower intervals
!>                    2- px,py,pz momentum upper intervals
!> outputs:
!>   particle: (particle_base) particle with sample momenta
subroutine sampling_cartesian_p_kinetic_relativistic(particle,&
rand,pxpypz_interval)
  use mod_particle_types, only: particle_kinetic_relativistic
  implicit none
  !> inputs-outputs:
  type(particle_kinetic_relativistic),intent(inout) :: particle
  !> inputs:
  real*8,dimension(3),intent(in)   :: rand
  real*8,dimension(3,2),intent(in) :: pxpypz_interval
  !> sample the cartesian momentum
  !> WARNING: the momentum must be in unit of AMU
  particle%p = pxpypz_interval(:,1)+(pxpypz_interval(:,2)-pxpypz_interval(:,1))*rand
end subroutine sampling_cartesian_p_kinetic_relativistic

!> uniform sampling of the cartesiam momentum 
!> for gc relativistic particles
!> inputs:
!>   gc:              (particle_gc_relativistic) relativistic gc to be sampled
!>   fields:          (fields_base) MHD fields type object
!>   rand:            (real8)(3) random numbers
!>   time:            (real8) time at which perform the interpolation
!>   mass:            (real8) particle mass in AMU
!>   pxpypz_interval: (real8)(3,2) momenta interval (in units of AMU)
!>                    1- px,py,pz momentum lower intervals
!>                    2- px,py,pz momentum upper intervals
!> outputs:
!>   gc_inout: (particle_gc_relativistic) relativistic gc with sampled momentum
subroutine sampling_cartesian_p_gc_relativistic(gc,fields,&
rand,time,mass,pxpypz_interval)
  use mod_particle_types,       only: particle_base
  use mod_particle_types,       only: particle_gc_relativistic
  use mod_fields,               only: fields_base
  use mod_kinetic_relativistic, only: momentum_relativistic_kinetic_to_relativistic_gc
  implicit none
  !> inputs-outputs:
  type(particle_gc_relativistic),intent(inout) :: gc
  class(fields_base),intent(inout) :: fields
  !> inputs:
  real*8,intent(in)                :: time,mass
  real*8,dimension(3),intent(in)   :: rand
  real*8,dimension(3,2),intent(in) :: pxpypz_interval
  !> variables:
  real*8 :: psi,U,B_norm 
  real*8,dimension(3) :: momentum,B,E
  !> sample the relativistic gc
  !> WARNING: the momentum must be in units of AMU
  momentum = pxpypz_interval(:,1)+(pxpypz_interval(:,2)-pxpypz_interval(:,1))*rand
  call fields%calc_EBpsiU(time,gc%i_elm,gc%st,gc%x(3),&
  E,B,psi,U)
  B_norm = sqrt(B(1)**2+B(2)**2+B(3)**2); B = B/B_norm;
  gc%p = momentum_relativistic_kinetic_to_relativistic_gc(mass,&
  momentum,gc%x(3),B_norm,B)
end subroutine sampling_cartesian_p_gc_relativistic

!> uniform sampling of the energy, pitch angle
!> and gyro angle for kinetic relativistic particles
!> inputs:
!>   particle:       (particle_kinetic_relativistic) particle to be sampled
!>   fields:         (fields_base) MHD fields type object
!>   rand:           (real8)(3) random numbers
!>   time:           (real8) time of the initialisation
!>   p_cube_int:     (real8)(2) cube of the particle momentum sampling interval
!>   cos_theta_int:  (real8)(2) cosinus of the pitch angle sampling interval
!>   gyro_int:       (real8)(2) gyro angle sampling interval
!> outputs:
!>   particle: (particle_kinetic_relativistic) particle with sampled momentum
subroutine sampling_uniform_ppitchgyro_kinetic_relativistic(particle,&
fields,rand,time,p_cube_int,cos_theta_int,gyro_int)
  use mod_sampling,             only: sample_uniform_sphere_corona_rcosphi
  use mod_particle_types,       only: particle_kinetic_relativistic
  use mod_fields,               only: fields_base
  use mod_kinetic_relativistic, only: kinetic_relativistic_momentum_spherical_to_cart
  implicit none
  !> inputs-outputs:
  type(particle_kinetic_relativistic),intent(inout) :: particle
  class(fields_base),intent(inout) :: fields
  !> inputs:
  real*8,intent(in)               :: time
  real*8,dimension(2),intent(in)  :: p_cube_int,cos_theta_int,gyro_int
  real*8,dimension(3),intent(in)  :: rand
  !> variables:
  real*8              :: psi,U
  real*8,dimension(3) :: B,E
  !> compute the magnetic field direction
  call fields%calc_EBpsiU(time,particle%i_elm,particle%st,&
  particle%x(3),E,B,psi,U)
  B = B/sqrt(B(1)**2+B(2)**2+B(3)**2)
  !> sample the momentum in spherical coordinates
  particle%p = sample_uniform_sphere_corona_rcosphi(p_cube_int,&
  cos_theta_int,gyro_int,rand)
  !> compute the particle kinetic momentum
  particle%p = kinetic_relativistic_momentum_spherical_to_cart(&
  particle%x(3),particle%p,B)
end subroutine sampling_uniform_ppitchgyro_kinetic_relativistic

!> uniform sampling of the energy, pitch angle
!> and gyro angle for gc relativistic particles
!> inputs:
!>   particle:       (particle_gc_relativistic) particle to be sampled
!>   fields:         (fields_base) MHD fields type object
!>   rand:           (real8)(3) random numbers
!>   mass:           (real8) particle mass in AMU
!>   time:           (real8) time of the initialisation
!>   p_cube_int:     (real8)(2) cube of the gc momentum sampling interval
!>   cos_theta_int:  (real8)(2) cosinus of the pitch angle sampling interval
!>   gyro_int:       (real8)(2) gyro angle sampling interval
!> outputs:
!>   particle: (particle_gc_relativistic) particle with sampled momentum
subroutine sampling_uniform_ppitchgyro_gc_relativistic(gc,fields,&
rand,mass,time,p_cube_int,cos_theta_int,gyro_int)
  use mod_sampling,             only: sample_uniform_sphere_corona_rcosphi
  use mod_particle_types,       only: particle_gc_relativistic
  use mod_fields,               only: fields_base
  use mod_kinetic_relativistic, only: kinetic_relativistic_momentum_spherical_to_cart
  use mod_kinetic_relativistic, only: momentum_relativistic_kinetic_to_relativistic_gc
  implicit none
  !> inputs-outputs:
  type(particle_gc_relativistic),intent(inout) :: gc
  class(fields_base),intent(inout) :: fields
  !> inputs:
  real*8,intent(in)               :: time,mass
  real*8,dimension(2),intent(in)  :: p_cube_int,cos_theta_int,gyro_int
  real*8,dimension(3),intent(in)  :: rand
  !> variables
  real*8              :: B_norm,psi,U
  real*8,dimension(3) :: p_kin,B,E
  !> compute magnetic field direction
  call fields%calc_EBpsiU(time,gc%i_elm,gc%st,gc%x(3),&
  E,B,psi,U)
  B_norm = sqrt(B(1)**2+B(2)**2+B(3)**2); B = B/B_norm;
  !> sample the momentum in spherical coordinates
  p_kin = sample_uniform_sphere_corona_rcosphi(p_cube_int,&
  cos_theta_int,gyro_int,rand)
  !> compute the particle kinetic momentum
  p_kin = kinetic_relativistic_momentum_spherical_to_cart(gc%x(3),p_kin,B)
  !> transform in guiding center coords
  gc%p = momentum_relativistic_kinetic_to_relativistic_gc(mass,&
  p_kin,gc%x(3),B_norm,B)
end subroutine sampling_uniform_ppitchgyro_gc_relativistic

!> sampling particle charge
!> inputs:
!>   rand:  (real8) random number
!>   q_int: (integer1)(2) charge interval
!> outputs:
!>   q:     (integer1) selected charge
function sampling_uniform_charge(rand,q_int) result(q)
  implicit none
  !> inputs:
  integer*1,dimension(2),intent(in) :: q_int
  real*8,intent(in)                 :: rand
  !> outputs:
  integer*1 :: q
  q = floor(real(q_int(1),kind=8)+(real(q_int(2),kind=8)-&
  real(q_int(1),kind=8)+1.d0)*rand,kind=1)
end function sampling_uniform_charge

!> find the bounding box in RZ
!> inputs:
!>   fields: (fields_base) MHD fields object
!> outputs:
!>   R_minmax:   (real8)(2) minimum and maximum major radius (global)
!>   Z_minmax:   (real8)(2) minimum and maximum vertical coord. (global)
!>   phi_minmax: (real8)(2) minimum and maximum toroidal angle (global)
subroutine find_RZPhi_minmax_global(R_minmax,Z_minmax,phi_minmax,fields)
  use constants,  only: TWOPI
  use mod_fields, only: fields_base
  implicit none
  !> inputs:
  class(fields_base),intent(in)   :: fields
  !> outputs:
  real*8,dimension(2),intent(out) :: R_minmax,Z_minmax,phi_minmax
  !> variables:
  integer :: ii
  real*8 :: R_min_elem,R_max_elem,Z_min_elem,Z_max_elem
  !> find bounding box
  R_minmax=(/1.d16,-1.d16/); Z_minmax=R_minmax; phi_minmax=(/0.d0,TWOPI/);
  do ii=1,fields%element_list%n_elements
    call RZ_minmax(fields%node_list,fields%element_list,ii,&
    R_min_elem,R_max_elem,Z_min_elem,Z_max_elem)
    R_minmax(1) = min(R_min_elem,R_minmax(1))
    R_minmax(2) = max(R_max_elem,R_minmax(2))
    Z_minmax(1) = min(Z_min_elem,Z_minmax(1))
    Z_minmax(2) = max(Z_max_elem,Z_minmax(2))
  enddo
end subroutine find_RZPhi_minmax_global

!> check RZPhi interval in bounding box
!> inputs:
!>  R_box:      (real8)(2) major radius box for uniform sampling
!>  Z_box:      (real8)(2) vertical coord. box for uniform sampling
!>  phi_box:    (real8)(2) toroidal angle box for uniform sampling
!>  R_minmax:   (real8)(2) minimum and maximum major radius (global)
!>  Z_minmax:   (real8)(2) minimum and maximum vertical coord. (global)
!>  phi_minmax: (real8)(2) minimum and maximum toroidal angle (global)
!> outputs:
!>  R_box:      (real8)(2) major radius box for uniform sampling
!>  Z_box:      (real8)(2) vertical coord. box for uniform sampling
!>  phi_box:    (real8)(2) toroidal angle box for uniform sampling
subroutine check_RZPhi_interval(R_box,Z_box,Phi_box,&
R_minmax,Z_minmax,Phi_minmax)
  implicit none
  !> input-outputs:
  real*8,dimension(2),intent(inout) :: R_box,Z_box,phi_box
  !> inputs:
  real*8,dimension(2),intent(in)    :: R_minmax,Z_minmax,phi_minmax
  !> check boundin boxes limits
  R_box(1)=max(R_box(1),R_minmax(1)); R_box(2)=min(R_box(2),R_minmax(2));
  Z_box(1)=max(Z_box(1),Z_minmax(1)); Z_box(2)=min(Z_box(2),Z_minmax(2));
  call check_phi_interval(phi_box,phi_minmax)
end subroutine check_RZPhi_interval

!> check the toroidal angle interval
!> inputs:
!>  phi_box:    (real8)(2) toroidal angle box for uniform sampling
!>  phi_minmax: (real8)(2) minimum and maximum toroidal angle (global)
!> outputs:
!>  phi_box:    (real8)(2) toroidal angle box for uniform sampling
subroutine check_phi_interval(phi_box,phi_minmax)
  implicit none
  !> inputs-outputs:
  real*8,dimension(2),intent(inout) :: phi_box
  !> inputs:
  real*8,dimension(2),intent(in)    :: phi_minmax
  !> check in bound
  phi_box(1) = max(phi_box(1),phi_minmax(1))
  phi_box(2) = min(phi_box(2),phi_minmax(2))
end subroutine check_phi_interval

!> find maximum and minimum poloidal flux: global and per element
!> inputs:
!>   fields: (fields_base) MHD fields object
!> outputs:
!>   psi_minmax_global: (real8)(2) maximum and minimum poloidal fluxes
!>   psi_minmax_list:   (real8)(n_elements,2) maxmum and minimum
!>                      poloidal fluxes for each mesh element
!>   theta_minmax:      (real8)(2) maximum and minimum poloidal angle
!>   phi_minmax:        (real8)(2) maximum and minimum toroidal angle
subroutine find_psithetaphi_minmax_global_list(psi_minmax_global,&
psi_minmax_list,theta_minmax,phi_minmax,fields)
  use constants,  only: TWOPI
  use mod_fields, only: fields_base
  implicit none
  !> inputs:
  class(fields_base),intent(in) :: fields
  !> outputs:
  real*8,dimension(2),intent(out) :: psi_minmax_global,theta_minmax
  real*8,dimension(2),intent(out) :: phi_minmax
  real*8,dimension(fields%element_list%n_elements,2),intent(out) :: psi_minmax_list
  !> variables
  integer :: ii
  !> compute maximum and minimum of the poloidal flux
  theta_minmax = (/0.d0,TWOPI/); phi_minmax = (/0.d0,TWOPI/);
  psi_minmax_global=(/1.d16,-1.d16/)
  do ii=1,fields%element_list%n_elements
    call psi_minmax(fields%node_list,fields%element_list,ii,&
    psi_minmax_list(ii,1),psi_minmax_list(ii,2))
    psi_minmax_global(1) = min(psi_minmax_global(1),psi_minmax_list(ii,1))
    psi_minmax_global(2) = max(psi_minmax_global(2),psi_minmax_list(ii,2))
  enddo
end subroutine find_psithetaphi_minmax_global_list

!> denormalise and check the poloidal flux interval
!> inputs:
!>   theta_box:    (real8)(2) poloidal angle bounding box
!>   psi_box:      (real8)(2) poloidal flux bounding box
!>   phi_box:      (real8)(2) toroidal angle bounding box
!>   psi_axis:     (real8) poloidal flux at the magnetic axis
!>   psi_boundary: (real8) poloidal flux at the plasma boundary
!>   psi_minmax:   (real8)(2) minimum-maximum poloidal flux
!>   theta_minmax: (real8)(2) minimum-maxumum poloidal angle
!>   phi_minmax:   (real8)(2) minimum-maximum toroidal angle
!> outputs:
!>   theta_box:    (real8)(2) poloidal angle bounding box
!>   psi_box:      (real8)(2) poloidal flux bounding box
!>   phi_box:      (real8)(2) toroidal angle bounding box
subroutine check_thetapsiphi_interval(theta_box,psi_box,phi_box,&
psi_axis,psi_boundary,psi_minmax,theta_minmax,phi_minmax)
  implicit none
  !> inputs-outputs:
  real*8,dimension(2),intent(inout) :: theta_box,psi_box,phi_box
  !> inputs:
  real*8,intent(in)              :: psi_axis,psi_boundary
  real*8,dimension(2),intent(in) :: psi_minmax,theta_minmax,phi_minmax
  !> denormalize and check boundaries
  psi_box      = psi_box*(psi_boundary-psi_axis)+psi_axis
  psi_box(1)   = max(psi_box(1),psi_minmax(1))
  psi_box(2)   = min(psi_box(2),psi_minmax(2))
  theta_box(1) = max(theta_box(1),theta_minmax(1))
  theta_box(2) = min(theta_box(2),theta_minmax(2))
  call check_phi_interval(phi_box,phi_minmax)
end subroutine check_thetapsiphi_interval

!> compute the momentum from particle energy and check bounds
!> inputs:
!>   e_interval:     (real8)(2) kinetic energy interval in eV
!>   pitch_interval: (real8)(2) pitch angle interval
!>   gyro_interval:  (real8)(2) gyro angle interval
!>   mass:           (real8) particle mass in AMU
!> outputs:
!>   p_interval:     (real8)(2) momentum interval in AMU*m/s
!>   pitch_interval: (real8)(2) pitch angle interval
!>   gyro_interval:  (real8)(2) gyro angle interval
subroutine check_energykinpitchgyro_interval(e_interval,p_interval,&
pitch_interval,gyro_interval,mass)
  use constants, only: PI,TWOPI
  use constants, only: EL_CHG,ATOMIC_MASS_UNIT,SPEED_OF_LIGHT
  implicit none
  !> inputs-outputs:
  real*8,dimension(2),intent(inout) :: pitch_interval
  real*8,dimension(2),intent(inout) :: gyro_interval
  !> inputs:
  real*8,intent(in)   :: mass
  real*8,dimension(2) :: e_interval
  !> outpus:
  real*8,dimension(2) :: p_interval
  !> variables:
  real*8 :: E0
  !> check if the energy is not smaller thant the rest energy in eV
  E0 = (ATOMIC_MASS_UNIT*mass*SPEED_OF_LIGHT*SPEED_OF_LIGHT)/EL_CHG
  p_interval = (e_interval/E0) + 1.d0 ;
  if(p_interval(1).le.1.d0) p_interval(1) = 1.d0+1.d-13
  if(p_interval(2).le.1.d0) p_interval(2) = 1.d0+1.d-13
  !> transform the energy in momentum intensity
  p_interval = mass*SPEED_OF_LIGHT*sqrt(p_interval*p_interval-1.d0)
  !> check the pitch angle
  pitch_interval(1) = max(pitch_interval(1),0.d0)
  pitch_interval(2) = min(pitch_interval(2),PI)
  !> check the pitch angle
  gyro_interval(1) = max(gyro_interval(1),0.d0)
  gyro_interval(2) = min(gyro_interval(2),TWOPI)
end subroutine check_energykinpitchgyro_interval

!> initialise particle base fields with unit particle weight
!> inputs:
!>   p_inout: (particle_base) particle to be initialised to zero
!> outputs:
!>   p_inout: (particlle_base) initialised particle to zero
subroutine init_particle_base_to_zero(p_inout)
  use mod_particle_types, only: particle_base
  implicit none
  class(particle_base),intent(inout) :: p_inout
  !> reset the particle base fields
  p_inout%x=0.d0; p_inout%st=0.d0; p_inout%weight=1.d0;
  p_inout%i_elm=0; p_inout%i_life=0; p_inout%t_birth=0.d0;
end subroutine init_particle_base_to_zero

!> initialise particle kinetic relativistic fields with unit particle weight
!> inputs:
!>   p_inout: (particle_kinetic_relativistic) particle to be initialised to zero
!> outputs:
!>   p_inout: (particlle_kinetic_relativistic) initialised particle to zero
subroutine init_particle_kinetic_relativistic_to_zero(p_inout)
  use mod_particle_types, only: particle_kinetic_relativistic
  implicit none
  type(particle_kinetic_relativistic),intent(inout) :: p_inout
  !> reset particle kinetic relativistic
  call init_particle_base_to_zero(p_inout)
  p_inout%p = 0.d0; p_inout%q = 0;
end subroutine init_particle_kinetic_relativistic_to_zero

!> initialise particle gc relativistic fields with unit particle weight
!> inputs:
!>   p_inout: (particle_gc_relativistic) particle to be initialised to zero
!> outputs:
!>   p_inout: (particlle_gc_relativistic) initialised particle to zero
subroutine init_particle_gc_relativistic_to_zero(p_inout)
  use mod_particle_types, only: particle_gc_relativistic
  implicit none
  type(particle_gc_relativistic),intent(inout) :: p_inout
  !> reset particle gc relativistic
  call init_particle_base_to_zero(p_inout)
  p_inout%p = 0.d0; p_inout%q = 0;
end subroutine init_particle_gc_relativistic_to_zero

!>---------------------------------------------------------
end module mod_initialise_relativistic_particles
