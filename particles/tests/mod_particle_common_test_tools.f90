!> module mod_particle_type_test_tools contains variables and
!> procedure used for initialising and finalising particle data
!> for unit testing
module mod_particle_common_test_tools
use constants, only: PI
implicit none

private
public :: n_particle_types,rng_seed_interval
public :: q_interval,i_elm_interval,i_life_interval,sim_time_interval
public :: t_birth_interval,st_interval,mass_interval,v_interval
public :: Ekin_interval,mu_interval,Bnorm_interval,weight_interval
public :: x_lowbnd,x_uppbnd,vp3d_lowbnd,vp3d_uppbnd,ABE_lowbnd,ABE_uppbnd
public :: RZPhi_lowbnd,RZPhi_uppbnd,BE0_lowbnd,BE0_uppbnd
public :: RZ0_lowbnd,RZ0_uppbnd,EThetaChi_RE_lowbnd,EThetaChi_RE_uppbnd
public :: q1_pos_interval,q1_neg_interval,q1_posneg_interval
public :: fill_particles,invalidate_particles,obtain_active_particle_ids
public :: fill_groups,fill_particle_base,fill_particle_fieldline
public :: fill_particle_gc,fill_particle_gc_vpar,fill_particle_gc_Qin
public :: fill_particle_kinetic,fill_particle_kinetic_leapfrog
public :: fill_particle_kinetic_relativistic,fill_particle_gc_relativistic
public :: obtain_particle_charges,allocate_one_particle_list_type
public :: copy_group_fieldline_B_hat_prev
public :: compute_test_E_B_fields
public :: fill_particle_kinetic_relativistic_RE,fill_particle_gc_relativistic_RE
public :: fill_particles_tokamak,fill_mass_RE
public :: fill_particle_simulations_no_init

!> Variables --------------------------------------------------
integer,parameter :: n_particle_types=8
integer*1,dimension(2),parameter :: q1_pos_interval=(/1,14/)
integer*1,dimension(2),parameter :: q1_neg_interval=(/-12,-1/)
integer*1,dimension(2),parameter :: q1_posneg_interval=(/-23,43/)
integer,dimension(2),parameter   :: rng_seed_interval=(/-1234,9876/)
integer,dimension(2),parameter   :: q_interval=(/1,100/)
integer,dimension(2),parameter   :: i_elm_interval=(/1,1000000/)
integer,dimension(2),parameter   :: i_life_interval=(/1,1000000/)
real*8,dimension(2),parameter    :: sim_time_interval=(/0.d0,1.d3/)
real*8,dimension(2),parameter    :: t_birth_interval=(/0.d0,3.45d4/)
real*8,dimension(2),parameter    :: st_interval=(/0.d0,1.d0/)
real*8,dimension(2),parameter    :: mass_interval=(/5.485d-4,124.d0/)
real*8,dimension(2),parameter    :: v_interval=(/-6.75d3,8.45d3/)
real*8,dimension(2),parameter    :: Ekin_interval=(/0.d0,1.d7/)
real*8,dimension(2),parameter    :: mu_interval=(/0.d0,1.d-5/)
real*8,dimension(2),parameter    :: Bnorm_interval=(/0.d0,1.4d1/)
real*8,dimension(2),parameter    :: weight_interval=(/0.d0,1.d3/)
real*8,dimension(3),parameter    :: x_lowbnd=(/-5.d2,-1.d2,1.d0/)
real*8,dimension(3),parameter    :: x_uppbnd=(/7.d2,2.d2,4.d2/)
real*8,dimension(3),parameter    :: vp3d_lowbnd=(/-1.25d3,-7.5d2,-8.d1/)
real*8,dimension(3),parameter    :: vp3d_uppbnd=(/7.5d1,2.35d2,4.85d3/)
real*8,dimension(3),parameter    :: ABE_lowbnd=(/-2.67d0,-9.85d0,0.35d0/)
real*8,dimension(3),parameter    :: ABE_uppbnd=(/0.78d0,2.35d0,5.67d0/)
real*8,dimension(3),parameter    :: RZPhi_lowbnd=(/0.d0,-1.5d0,0.d0/)
real*8,dimension(3),parameter    :: RZPhi_uppbnd=(/1.d0,1.5d0,2.d0*PI/)
real*8,dimension(2),parameter    :: Rminmax=(/2.5d0,3.8d0/)
!> parameter for tets electric and magnetic fields
real*8,parameter :: B0=3.5d0 !< toroidal magnetic field on axis
real*8,parameter :: R0=3.d0  !< axis major radius
real*8,parameter :: Z0=1.d-1 !< axis vertical position
real*8,parameter :: E0=5.3d0 !< toroidal electric field on axis
real*8,dimension(2),parameter :: RZ0_lowbnd=(/5.d-1,-1.5d0/)
real*8,dimension(2),parameter :: RZ0_uppbnd=(/6.d0,1.5d0/)
real*8,dimension(2),parameter :: BE0_lowbnd=(/1.d-1,-5.2d0/)
real*8,dimension(2),parameter :: BE0_uppbnd=(/6.d0,5.2d0/)
real*8,parameter,dimension(3) :: EThetaChi_RE_lowbnd=(/1.d5,PI/3.d0,0.d0/) !< for RE
real*8,parameter,dimension(3) :: EThetaChi_RE_uppbnd=(/5.d7,2.d0*PI/3.d0,2.d0*PI/) !< for RE
!> Interfaces -------------------------------------------------
interface allocate_one_particle_list_type
  module procedure allocate_one_particle_list_all_types
  module procedure allocate_one_particle_list_from_types
end interface allocate_one_particle_list_type

interface fill_groups
  module procedure fill_groups_seq
  module procedure fill_groups_mpi
end interface fill_groups

interface compute_test_E_B_fields
  module procedure compute_test_E_B_fields_parabolic
end interface compute_test_E_B_fields

interface compute_test_gradpsi
  module procedure compute_test_gradpsi_parabolic
end interface compute_test_gradpsi

contains
!> Procedures -------------------------------------------------
!> compute test electric and magnetic fields using
!> a parabolic poloidal flux:
!> psi = B0*((R-R0)**2+(Z-Z0)**2)/2
!> and a radially dependent toroidal magnetic/electric fields
!> inputs:
!>   x: (real8)(3) position in cylindrical coord. (R,Z,phi)
!> outpus:
!>   E_field: (real8)(3) electric field cylindrical coord.
!>   B_field: (real8)(3) magnetic field cylindrical coord.
subroutine compute_test_E_B_fields_parabolic(x,E_field,B_field)
  implicit none
  !> inputs:
  real*8,dimension(3),intent(in) :: x
  !> outputs:
  real*8,dimension(3),intent(out) :: E_field,B_field
  !> compute magnetic field
  !> compute the electric field
  B_field = B0*(/x(2)-Z0,R0-x(1),R0/)/x(1)
  E_field = (/0.d0,0.d0,-E0*R0/x(1)/)
end subroutine compute_test_E_B_fields_parabolic

!> compute test gradient of the poloidal flux
!> a parabolic poloidal flux:
!> psi = B0*((R-R0)**2+(Z-Z0)**2)/2
!> inputs:
!>   x: (real8)(3) position in cylindrical coord. (R,Z,phi)
!> outpus:
!>   gradpsi: (real8)(3) grad poloidal flux
subroutine compute_test_gradpsi_parabolic(x,gradpsi)
  implicit none
  !> inputs:
  real*8,dimension(3),intent(in) :: x
  !> outputs:
  real*8,dimension(3),intent(out) :: gradpsi
  !> compute the grad psi
  gradpsi = B0*(/x(1)-R0,x(2)-Z0,0.d0/)
end subroutine compute_test_gradpsi_parabolic

!> fill particle simulations without direct initialisation
!> if fill_particle_in = 1 -> cylindrical initialisation is used
!>                       cartesian otherwise
subroutine fill_particle_simulations_no_init(sims_particles,n_times,n_groups,&
n_particles,p_types,ifail,rank_in,n_tasks_in,fill_particle_in)
  use mpi
  use mod_gnu_rng,      only: gnu_rng_interval
  use mod_particle_sim, only: particle_sim
  implicit none
  !> inputs-outputs:
  type(particle_sim),dimension(n_times),intent(inout) :: sims_particles
  integer,intent(inout) :: ifail
  !> inputs:
  integer,intent(in) :: n_times,n_groups
  integer,dimension(n_groups),intent(in) :: n_particles,p_types
  integer,intent(in),optional :: rank_in,n_tasks_in,fill_particle_in
  !> variables:
  integer :: ii,jj,rank,n_tasks,fill_particle
  real*8,dimension(n_times*(n_groups+1)) :: sims_real8
  !> initialisation
  rank = 0; if(present(rank_in)) rank = rank_in;
  n_tasks = 1; if(present(n_tasks_in)) n_tasks = n_tasks_in;
  fill_particle = 1; if(present(fill_particle_in)) fill_particle = fill_particle_in;
  !> fill the simulations without fields
  if(rank.eq.0) then 
    call gnu_rng_interval(n_times,sim_time_interval,sims_real8(1:n_times))
    call gnu_rng_interval(n_groups*n_times,mass_interval,sims_real8(n_times+1:n_times*(n_groups+1)))
  endif
  call MPI_Bcast(sims_real8,n_times*(n_groups+1),MPI_DOUBLE,0,MPI_COMM_WORLD,ifail)
  do ii=1,n_times
    sims_particles(ii)%my_id = rank; sims_particles(ii)%n_cpu = n_tasks;
    sims_particles(ii)%wtime_start = MPI_Wtime()
    sims_particles(ii)%time = sims_real8(ii)
    call sims_particles(ii)%set_t_norm
    !> allocate and initialise groups
    call sims_particles(ii)%allocate_groups(n_groups)
    do jj=1,n_groups
      sims_particles(ii)%groups(jj)%mass = sims_real8((ii-1)*n_groups+jj+n_times)
    enddo
    call fill_mass_RE(n_groups,sims_particles(ii)%groups)
    !> allocate and initialise particles
    call allocate_one_particle_list_type(n_groups,n_particles,p_types,&
    sims_particles(ii)%groups,ifail)
    call fill_particles_tokamak(n_groups,sims_particles(ii)%groups,fill_particle,rank)
  enddo
end subroutine fill_particle_simulations_no_init

!> allocate particle list as a function of the particle type
subroutine allocate_one_particle_list_from_types(n_groups,n_particles,p_types,groups,ifail)
  use mod_particle_sim,   only: particle_group
  use mod_particle_types, only: particle_kinetic,particle_kinetic_leapfrog
  use mod_particle_types, only: particle_gc,particle_fieldline
  use mod_particle_types, only: particle_kinetic_relativistic
  use mod_particle_types, only: particle_gc_relativistic
  use mod_particle_types, only: particle_gc_vpar,particle_gc_Qin 
  use mod_particle_types, only: particle_kinetic_id,particle_kinetic_leapfrog_id
  use mod_particle_types, only: particle_gc_id,particle_fieldline_id
  use mod_particle_types, only: particle_kinetic_relativistic_id
  use mod_particle_types, only: particle_gc_relativistic_id
  use mod_particle_types, only: particle_gc_vpar_id,particle_gc_Qin_id
  implicit none
  !> inputs
  integer,intent(in) :: n_groups
  integer,dimension(n_groups),intent(in) :: p_types,n_particles
  !> inputs-outputs
  integer,intent(inout) :: ifail
  type(particle_group),dimension(n_groups),intent(inout) :: groups
  !> variables
  integer :: ii
  
  !> initialize the particle as a function of the type
  do ii=1,n_groups
    if(p_types(ii).eq.particle_fieldline_id) then
      allocate(particle_fieldline::groups(ii)%particles(n_particles(ii)))
    elseif(p_types(ii).eq.particle_gc_id) then
      allocate(particle_gc::groups(ii)%particles(n_particles(ii)))
    elseif(p_types(ii).eq.particle_gc_vpar_id) then
      allocate(particle_gc_vpar::groups(ii)%particles(n_particles(ii)))
    elseif(p_types(ii).eq.particle_kinetic_id) then
      allocate(particle_kinetic::groups(ii)%particles(n_particles(ii)))
    elseif(p_types(ii).eq.particle_kinetic_leapfrog_id) then
      allocate(particle_kinetic_leapfrog::groups(ii)%particles(n_particles(ii)))
    elseif(p_types(ii).eq.particle_kinetic_relativistic_id) then
      allocate(particle_kinetic_relativistic::groups(ii)%particles(n_particles(ii)))
    elseif(p_types(ii).eq.particle_gc_relativistic_id) then 
      allocate(particle_gc_relativistic::groups(ii)%particles(n_particles(ii)))
    elseif(p_types(ii).eq.particle_gc_Qin_id) then
      allocate(particle_gc_Qin::groups(ii)%particles(n_particles(ii)))
    else
      write(*,'(/A)') "allocate one particle list from type failed!"
      ifail = 1
    endif 
  enddo  
end subroutine allocate_one_particle_list_from_types

!> allocate one particle list per type
subroutine allocate_one_particle_list_all_types(n_groups,n_particles,groups,ifail)
  use mod_particle_sim,   only: particle_group
  use mod_particle_types, only: particle_kinetic,particle_kinetic_leapfrog
  use mod_particle_types, only: particle_gc,particle_fieldline
  use mod_particle_types, only: particle_kinetic_relativistic
  use mod_particle_types, only: particle_gc_relativistic
  use mod_particle_types, only: particle_gc_vpar,particle_gc_Qin 
  implicit none
  !> inputs
  integer,intent(in) :: n_groups,n_particles
  !> inputs-outputs
  integer,intent(inout) :: ifail
  type(particle_group),dimension(n_groups),intent(inout) :: groups
  if(n_groups.lt.n_particle_types) then
    write(*,'(/A)') "allocate one particle list for all types failed!"
    ifail = 1
    return
  endif
  allocate(particle_fieldline::groups(1)%particles(n_particles))
  allocate(particle_gc::groups(2)%particles(n_particles))
  allocate(particle_gc_vpar::groups(3)%particles(n_particles))
  allocate(particle_kinetic::groups(4)%particles(n_particles))
  allocate(particle_kinetic_leapfrog::groups(5)%particles(n_particles))
  allocate(particle_kinetic_relativistic::groups(6)%particles(n_particles))
  allocate(particle_gc_relativistic::groups(7)%particles(n_particles))
  allocate(particle_gc_Qin::groups(8)%particles(n_particles))
end subroutine allocate_one_particle_list_all_types

!> obtain charges from all particles in a simulation 
subroutine obtain_particle_charges(n_groups,n_particles_max,charge_list,groups)
  use mod_particle_sim,   only: particle_group
  use mod_particle_types, only: particle_kinetic,particle_kinetic_leapfrog
  use mod_particle_types, only: particle_gc
  use mod_particle_types, only: particle_kinetic_relativistic
  use mod_particle_types, only: particle_gc_relativistic
  use mod_particle_types, only: particle_gc_vpar,particle_gc_Qin 
  implicit none
  !> inputs
  integer,intent(in) :: n_groups,n_particles_max
  !> inputs-outputs
  type(particle_group),dimension(n_groups),intent(inout) :: groups
  !> outputs
  integer*1,dimension(n_particles_max,n_groups),intent(out) :: charge_list
  !> variables
  integer :: ii,jj

  !> initialisation
  charge_list = 0
  !> extract charges
  !$omp parallel default(shared) firstprivate(n_groups) &
  !$omp private(ii,jj)
  do jj=1,n_groups
    !$omp do
    do ii=1,size(groups(jj)%particles)
      select type (p=>groups(jj)%particles(ii))
        type is (particle_kinetic)
        charge_list(ii,jj) = p%q
        type is (particle_kinetic_leapfrog)
        charge_list(ii,jj) = p%q
        type is (particle_gc)
        charge_list(ii,jj) = p%q
        type is (particle_kinetic_relativistic)
        charge_list(ii,jj) = p%q
        type is (particle_gc_relativistic)
        charge_list(ii,jj) = p%q
        type is (particle_gc_vpar)
        charge_list(ii,jj) = p%q
        type is (particle_gc_Qin)
        charge_list(ii,jj) = p%q
      end select
    enddo
    !$omp end do
  enddo
  !$omp end parallel
end subroutine obtain_particle_charges

!> extract the index of valid particles in a simulatiion
subroutine obtain_active_particle_ids(n_groups,n_particles_max,&
active_particle_ids,groups)
  use mod_particle_sim, only: particle_group
  implicit none
  !> inputs
  integer,intent(in) :: n_groups,n_particles_max
  !> inputs-outputs
  type(particle_group),dimension(n_groups),intent(inout) :: groups
  !> outputs:
  integer,dimension(n_particles_max,n_groups),intent(out) :: active_particle_ids
  !> variables
  integer :: ii,jj,counter
  do jj=1,n_groups
    counter = 0
    do ii=1,size(groups(jj)%particles)
      if(groups(jj)%particles(ii)%i_elm.gt.0) then
        counter = counter + 1
        active_particle_ids(counter,jj) = ii
      endif
    enddo
  enddo
end subroutine obtain_active_particle_ids

!> invalidate some of the particles in the particle groups of a
!> simulation as a function of a survival probability
subroutine invalidate_particles(n_groups,n_particles_max,survival_threshold_in,&
n_active_particles,groups,rank_in)
  !$ use omp_lib
  use mod_particle_sim, only: particle_group
  use mod_gnu_rng, only: gnu_rng_interval
  use mod_gnu_rng, only: set_seed_sys_time 
  implicit none
  !> inputs
  integer,intent(in) :: n_groups,n_particles_max
  real*8,intent(in)  :: survival_threshold_in
  integer,intent(in),optional :: rank_in
  !> inputs-outputs
  type(particle_group),dimension(n_groups),intent(inout) :: groups
  !> outputs
  integer,dimension(n_groups),intent(out) :: n_active_particles
  !> variables
  integer :: ii,jj,rank,thread_id
  real*8  :: survival_threshold
  real*8,dimension(n_particles_max,n_groups) :: survival_array
  integer,dimension(n_groups) :: n_particles
  !> initialisation
  rank = 1; n_particles=0; if(present(rank_in)) rank=rank_in;
  survival_threshold = abs(survival_threshold_in)
  if(abs(survival_threshold).gt.1.d0) survival_threshold = survival_threshold - floor(survival_threshold)
  n_active_particles = 0
  !$omp parallel default(shared) private(thread_id)
  !$ thread_id = omp_get_thread_num()
  call set_seed_sys_time(rng_seed_interval,rank,thread_id)
  !$omp end parallel
  
  !> generate survival probability
  call gnu_rng_interval(n_particles_max,n_groups,(/0.d0,1.d0/),survival_array)
  do ii=1,n_groups
    n_particles(ii) = size(groups(ii)%particles)
    n_active_particles(ii) = count(survival_array(1:n_particles(ii),ii).ge.survival_threshold)
  enddo
  !> invalidate particles
  !$omp parallel default(shared) firstprivate(n_groups,n_particles,survival_threshold) &
  !$omp private(ii,jj)
  do jj=1,n_groups
    !$omp do
    do ii=1,n_particles(jj)
      if(survival_array(ii,jj).lt.survival_threshold) groups(jj)%particles(ii)%i_elm = 0
    enddo
    !$omp end do
  enddo 
  !$omp end parallel
end subroutine invalidate_particles

!> copy fieldlines B_hat between two simulations 
!> used for IO because it is not stored in hdf5
subroutine copy_group_fieldline_B_hat_prev(n_groups,groups_in,groups_out)
  use mod_particle_types, only: particle_fieldline
  use mod_particle_sim,   only: particle_group
  implicit none
  integer,intent(in) :: n_groups
  type(particle_group),dimension(n_groups),intent(in)    :: groups_in
  type(particle_group),dimension(n_groups),intent(inout) :: groups_out
  integer :: ii,jj
  !$omp parallel default(private) firstprivate(n_groups) &
  !$omp shared(groups_in,groups_out)
  do ii=1,n_groups
    !$omp do
    do jj=1,size(groups_out(ii)%particles)
      select type (p_out=>groups_out(ii)%particles(jj))
      type is (particle_fieldline)
        select type (p_in=>groups_in(ii)%particles(jj))
          type is (particle_fieldline)
          p_out%B_hat_prev = p_in%B_hat_prev
        end select
      end select
    enddo
    !$omp end do
  enddo
  !$omp end parallel
end subroutine copy_group_fieldline_B_hat_prev

!> generate random values for filling the groups type
!> Sequential version
subroutine fill_groups_seq(n_groups,groups)
  use mod_particle_sim, only: particle_group
  use mod_gnu_rng, only: gnu_rng_interval
  use mod_gnu_rng, only: set_seed_sys_time
  implicit none
  integer,intent(in) :: n_groups
  type(particle_group),dimension(n_groups),intent(inout) :: groups
  !> variables
  integer :: ii
  real*8 :: rn_real
  call set_seed_sys_time(rng_seed_interval,1)
  do ii=1,n_groups
    call gnu_rng_interval(mass_interval,rn_real)
    groups(ii)%mass = rn_real
  enddo
end subroutine fill_groups_seq

!> generate random values for filling the groups type
!> we do not create random mass for all groups because
!> it seems that only rank 0 is saved in hdf5
subroutine fill_groups_mpi(n_groups,groups,rank,ifail)
  use mpi
  use mod_particle_sim, only: particle_group
  use mod_gnu_rng, only: gnu_rng_interval
  use mod_gnu_rng, only: set_seed_sys_time
  implicit none
  !> inputs
  integer,intent(in) :: n_groups,rank
  type(particle_group),dimension(n_groups),intent(inout) :: groups
  !> input-outputs
  integer,intent(inout) :: ifail
  !> variables
  integer :: ii
  real*8 :: rn_real
  real*8,dimension(n_groups) :: val_to_bcst
  if(rank.eq.0) then
    call set_seed_sys_time(rng_seed_interval,rank)
    do ii=1,n_groups
      call gnu_rng_interval(mass_interval,rn_real)
      groups(ii)%mass = rn_real
      val_to_bcst(ii) = groups(ii)%mass
    enddo
  endif
  call MPI_Bcast(val_to_bcst,n_groups,MPI_REAL8,0,MPI_COMM_WORLD,ifail)
  if(rank.ne.0) then
    do ii=1,n_groups
      groups(ii)%mass = val_to_bcst(ii)
    enddo
  endif
end subroutine fill_groups_mpi

!> fill the group mass with the electron mass in AMU is particle
!> is of a relativistic type
subroutine fill_mass_RE(n_groups,groups)
  use mod_particle_sim,   only: particle_group
  use mod_particle_types, only: particle_kinetic_relativistic
  use mod_particle_types, only: particle_gc_relativistic
  implicit none
  !> parameters
  real*8,parameter   :: mass_RE=5.48579909065d-4
  !> inputs
  integer,intent(in) :: n_groups
  !> inputs-outputs
  type(particle_group),dimension(n_groups),intent(inout) :: groups
  !> variables
  integer :: ii
  !> fill particle mass
  do ii=1,n_groups
    select type(p_list=>groups(ii)%particles)
      type is(particle_kinetic_relativistic)
      groups(ii)%mass = mass_re
      type is(particle_gc_relativistic)
      groups(ii)%mass = mass_re
    end select
  enddo
end subroutine fill_mass_RE

!> fill particle list with runaways electron as relativistic partices
!> or using the standard particle fill for all other types unless
!> a specific type implementation is provided
subroutine fill_particles_tokamak(n_groups,groups,fill_particle_in,rank_in)
  use mod_particle_sim,   only: particle_group
  use mod_particle_types, only: particle_kinetic,particle_kinetic_leapfrog
  use mod_particle_types, only: particle_gc,particle_fieldline
  use mod_particle_types, only: particle_kinetic_relativistic
  use mod_particle_types, only: particle_gc_relativistic
  use mod_particle_types, only: particle_gc_vpar,particle_gc_Qin
  implicit none
  integer,intent(in) :: n_groups
  type(particle_group),dimension(n_groups),intent(inout) :: groups
  integer,intent(in),optional :: fill_particle_in,rank_in
  integer :: fill_particle,rank,ii,n_particles
  fill_particle = 1; rank = 1;
  if(present(fill_particle_in)) fill_particle=fill_particle_in
  if(present(rank_in)) rank = rank_in
  !> fill particle type base
  call fill_particle_base(n_groups,fill_particle,groups,rank)
  !> fiil particle types 
  do ii=1,n_groups
    n_particles = size(groups(ii)%particles)
    select type (p_list=>groups(ii)%particles)
    type is(particle_fieldline)
    call fill_particle_fieldline(n_particles,p_list,rank)
    type is(particle_gc)
    call fill_particle_gc(n_particles,p_list,rank)
    type is(particle_gc_vpar)
    call fill_particle_gc_vpar(n_particles,p_list,rank)
    type is(particle_gc_Qin)
    call fill_particle_gc_Qin(n_particles,p_list,rank)
    type is(particle_kinetic)
    call fill_particle_kinetic(n_particles,p_list,rank)
    type is(particle_kinetic_leapfrog)
    call fill_particle_kinetic_leapfrog(n_particles,p_list,rank)
    type is(particle_kinetic_relativistic)
    call fill_particle_kinetic_relativistic_RE(n_particles,p_list,rank)
    type is(particle_gc_relativistic)
    call fill_particle_gc_relativistic_RE(n_particles,p_list,rank)
    end select
  enddo
end subroutine fill_particles_tokamak

!> fills particle list of each groups with random data
subroutine fill_particles(n_groups,groups,fill_particle_in,rank_in)
  use mod_particle_sim,   only: particle_group
  use mod_particle_types, only: particle_kinetic,particle_kinetic_leapfrog
  use mod_particle_types, only: particle_gc,particle_fieldline
  use mod_particle_types, only: particle_kinetic_relativistic
  use mod_particle_types, only: particle_gc_relativistic
  use mod_particle_types, only: particle_gc_vpar,particle_gc_Qin
  implicit none
  integer,intent(in) :: n_groups
  type(particle_group),dimension(n_groups),intent(inout) :: groups
  integer,intent(in),optional :: fill_particle_in,rank_in
  integer :: fill_particle,rank,jj,n_particles
  fill_particle = -1; rank = 1;
  if(present(rank_in)) rank = rank_in
  if(present(fill_particle_in)) fill_particle = fill_particle_in
  !> fill particle basic type
  call fill_particle_base(n_groups,fill_particle,groups,rank)
  !> fill particle specific types
  do jj=1,n_groups
    n_particles = size(groups(jj)%particles)
    select type (p_list=>groups(jj)%particles)
    type is(particle_fieldline)
    call fill_particle_fieldline(n_particles,p_list,rank)
    type is(particle_gc)
    call fill_particle_gc(n_particles,p_list,rank)
    type is(particle_gc_vpar)
    call fill_particle_gc_vpar(n_particles,p_list,rank)
    type is(particle_gc_Qin)
    call fill_particle_gc_Qin(n_particles,p_list,rank)
    type is(particle_kinetic)
    call fill_particle_kinetic(n_particles,p_list,rank)
    type is(particle_kinetic_leapfrog)
    call fill_particle_kinetic_leapfrog(n_particles,p_list,rank)
    type is(particle_kinetic_relativistic)
    call fill_particle_kinetic_relativistic(n_particles,p_list,rank)
    type is(particle_gc_relativistic)
    call fill_particle_gc_relativistic(n_particles,p_list,rank)
    end select
  enddo
end subroutine fill_particles

!> generate random values for filling the particle base type
!> parameter fill_type: default: fill particle in cartesian coordinates
!>                      1) fill particle base using cylindrical coordinates
subroutine fill_particle_base(n_groups,fill_type,groups,rank_in)
  use mod_particle_sim, only: particle_group
  use mod_gnu_rng, only: gnu_rng_interval
  use mod_gnu_rng, only: set_seed_sys_time
  !$ use omp_lib
  implicit none
  !> inputs
  integer,intent(in) :: n_groups,fill_type
  integer,intent(in),optional :: rank_in
  !> inputs-outputs:
  type(particle_group),dimension(n_groups),intent(inout) :: groups
  !> variables
  integer :: ii,jj,rank
  !$ integer :: thread_id
  rank = 1
  if(present(rank_in)) rank=rank_in
  thread_id = 0
  !> fill-up the particle_base variables for all particles and all groups
  !$omp parallel default(shared) firstprivate(n_groups) &
  !$omp private(ii,jj,rank,thread_id)
  !$ thread_id = omp_get_thread_num()
  call set_seed_sys_time(rng_seed_interval,rank,thread_id)
  if(fill_type.eq.1) then
    do jj=1,n_groups
      !$omp do
      do ii=1,size(groups(jj)%particles)
        call fill_particle_base_cyl(groups(jj)%particles(ii))
      enddo
      !$omp end do
    enddo
  else
    do jj=1,n_groups
      !$omp do
      do ii=1,size(groups(jj)%particles)
        call fill_particle_base_cart(groups(jj)%particles(ii))
      enddo
      !$omp end do
    enddo
  endif
  !$omp end parallel
end subroutine fill_particle_base

!> generate particle position in cylindrical coordinates
subroutine fill_particle_base_cyl(particle)
  use mod_particle_types, only: particle_base
  use mod_gnu_rng,        only: gnu_rng_interval
  implicit none
  class(particle_base),intent(inout) :: particle
  !> variables
  real*8,dimension(3) :: rn_real_size3
  !> generate random numbers
  call gnu_rng_interval(3,RZPhi_lowbnd,RZPhi_uppbnd,rn_real_size3)
  !> compute number major radius
  rn_real_size3(1) = sqrt(Rminmax(1)**2 + (Rminmax(2)**2 - &
  Rminmax(1)**2)*rn_real_size3(1))
  particle%x = rn_real_size3
  call fill_particle_base_nopos(particle)
end subroutine fill_particle_base_cyl

!> generate particle base in cartesian coordinates
subroutine fill_particle_base_cart(particle)
  use mod_particle_types, only: particle_base
  use mod_gnu_rng,        only: gnu_rng_interval
  implicit none
  !> inputs-outpus
  class(particle_base),intent(inout) :: particle
  !> variables
  real*8,dimension(3) :: rn_real_size3
  call gnu_rng_interval(3,x_lowbnd,x_uppbnd,rn_real_size3)
  particle%x = rn_real_size3
  call fill_particle_base_nopos(particle)
end subroutine fill_particle_base_cart

!> generate random values for particles without position
subroutine fill_particle_base_nopos(particle)
  use mod_particle_types, only: particle_base
  use mod_gnu_rng,        only: gnu_rng_interval
  implicit none
  !> inputs-outpus
  class(particle_base),intent(inout) :: particle
  !> variables
  integer             :: rn_integer
  real*8              :: rn_real
  real*8,dimension(2) :: rn_real_size2
  !> fill particle fields
  call gnu_rng_interval(weight_interval,rn_real)
  call gnu_rng_interval(2,st_interval,rn_real_size2)
  particle%st      = rn_real_size2
  particle%weight  = rn_real
  call gnu_rng_interval(t_birth_interval,rn_real)
  particle%t_birth = real(rn_real,kind=4)
  call gnu_rng_interval(i_elm_interval,rn_integer)
  particle%i_elm   = rn_integer
  call gnu_rng_interval(i_life_interval,rn_integer)
  particle%i_life  = rn_integer  
end subroutine fill_particle_base_nopos

!> fill up particle_fieldline with random numbers
subroutine fill_particle_fieldline(n_particles,particles,rank_in)
  use mod_particle_types, only: particle_fieldline
  use mod_gnu_rng,        only: gnu_rng_interval
  use mod_gnu_rng,        only: set_seed_sys_time
  !$ use omp_lib
  implicit none
  !> inputs
  integer,intent(in) :: n_particles
  integer,intent(in),optional :: rank_in
  type(particle_fieldline),dimension(n_particles),intent(inout) :: particles
  !> variables
  integer :: ii,rank
  !$ integer :: thread_id
  real*8              :: rn_real
  real*8,dimension(3) :: rn_real_size3
  rank = 1
  if(present(rank_in)) rank=rank_in
  thread_id = 0
  !$omp parallel default(private) firstprivate(n_particles,particles)
  !$ thread_id = omp_get_thread_num()
  call set_seed_sys_time(rng_seed_interval,rank,thread_id)
  !$omp do
  do ii=1,n_particles
    call gnu_rng_interval(3,ABE_lowbnd,ABE_uppbnd,rn_real_size3)
    call gnu_rng_interval(v_interval,rn_real)
    particles(ii)%B_hat_prev = rn_real_size3
    particles(ii)%v = rn_real
  enddo
  !$omp end do
  !$omp end parallel
end subroutine fill_particle_fieldline

!> fill up particle_gc with random numbers
subroutine fill_particle_gc(n_particles,particles,rank_in)
  use mod_particle_types, only: particle_gc
  use mod_gnu_rng,        only: gnu_rng_interval
  use mod_gnu_rng,        only: set_seed_sys_time
  !$ use omp_lib
  implicit none
  !> inputs
  integer,intent(in) :: n_particles
  integer,intent(in),optional :: rank_in
  type(particle_gc),dimension(n_particles),intent(inout) :: particles
  !> variables
  integer :: ii,rn_integer,rank
  !$ integer :: thread_id
  real*8 :: rn_real
  rank = 1
  if(present(rank_in)) rank=rank_in
  thread_id = 0
  !$omp parallel default(private) firstprivate(n_particles) shared(particles)
  !$ thread_id = omp_get_thread_num()
  call set_seed_sys_time(rng_seed_interval,rank,thread_id)
  !$omp do
  do ii=1,n_particles
    call gnu_rng_interval(Ekin_interval,rn_real)   
    particles(ii)%E = rn_real
    call gnu_rng_interval(mu_interval,rn_real)
    particles(ii)%mu = rn_real
    call gnu_rng_interval(q_interval,rn_integer)
    particles(ii)%q = int(rn_integer,kind=1)
  enddo
  !$omp end do
  !$omp end parallel
end subroutine fill_particle_gc

!> fill up particle_gc_vpar with random numbers
subroutine fill_particle_gc_vpar(n_particles,particles,rank_in)
  use mod_particle_types, only: particle_gc_vpar
  use mod_gnu_rng,        only: gnu_rng_interval
  use mod_gnu_rng,        only: set_seed_sys_time
  !$ use omp_lib
  implicit none
  !> inputs
  integer,intent(in) :: n_particles
  integer,intent(in),optional :: rank_in
  type(particle_gc_vpar),dimension(n_particles),intent(inout) :: particles
  !> variables
  integer :: ii,rn_integer,rank
  !$ integer :: thread_id
  real*8 :: rn_real
  rank = 1
  if(present(rank_in)) rank=rank_in
  !$omp parallel default(private) firstprivate(n_particles) shared(particles)
  thread_id = 0
  !$ thread_id = omp_get_thread_num()
  call set_seed_sys_time(rng_seed_interval,rank,thread_id)
  !$omp do
  do ii=1,n_particles
    call gnu_rng_interval(v_interval,rn_real)   
    particles(ii)%vpar = rn_real
    call gnu_rng_interval(mu_interval,rn_real)
    particles(ii)%mu = rn_real
    call gnu_rng_interval(Bnorm_interval,rn_real)
    particles(ii)%B_norm = rn_real
    call gnu_rng_interval(q_interval,rn_integer)
    particles(ii)%q = int(rn_integer,kind=1)
  enddo
  !$omp end do
  !$omp end parallel
end subroutine fill_particle_gc_vpar

!> fill up particle_gc_Qin with random numbers
subroutine fill_particle_gc_Qin(n_particles,particles,rank_in)
  use mod_particle_types, only: particle_gc_Qin
  use mod_gnu_rng,        only: gnu_rng_interval
  use mod_gnu_rng,        only: set_seed_sys_time
  !$ use omp_lib
  implicit none
  !> inputs
  integer,intent(in) :: n_particles
  integer,intent(in),optional :: rank_in
  type(particle_gc_Qin),dimension(n_particles),intent(inout) :: particles
  !> variables
  integer :: ii,rn_integer,rank
  !$ integer :: thread_id
  real*8              :: rn_real
  real*8,dimension(3) :: rn_real_size3
  real*8,dimension(3,3) :: rn_real_size33
  rank = 1
  if(present(rank_in)) rank=rank_in
  !$omp parallel default(private) firstprivate(n_particles) shared(particles)
  thread_id = 0
  !$ thread_id = omp_get_thread_num()
  call set_seed_sys_time(rng_seed_interval,rank,thread_id)
  !$omp do
  do ii=1,n_particles
    call gnu_rng_interval(3,x_lowbnd,x_uppbnd,rn_real_size3)
    particles(ii)%x_m = rn_real_size3       !< x_m
    call gnu_rng_interval(3,ABE_lowbnd,ABE_uppbnd,rn_real_size3)
    particles(ii)%Astar_m = rn_real_size3   !< Astar_m
    call gnu_rng_interval(3,ABE_lowbnd,ABE_uppbnd,rn_real_size3)
    particles(ii)%Astar_k = rn_real_size3   !< Astar_k
    call gnu_rng_interval(3,ABE_lowbnd,ABE_uppbnd,rn_real_size33(:,1))
    call gnu_rng_interval(3,ABE_lowbnd,ABE_uppbnd,rn_real_size33(:,2))
    call gnu_rng_interval(3,ABE_lowbnd,ABE_uppbnd,rn_real_size33(:,3))
    particles(ii)%dAstar_k = rn_real_size33 !< dAstar_k  
    call gnu_rng_interval(3,ABE_lowbnd,ABE_uppbnd,rn_real_size3)
    particles(ii)%dBn_k = rn_real_size3     !< dBn_k
    call gnu_rng_interval(3,ABE_lowbnd,ABE_uppbnd,rn_real_size3)
    particles(ii)%Bnorm_k = rn_real_size3   !< Bnorm_l
    call gnu_rng_interval(3,ABE_lowbnd,ABE_uppbnd,rn_real_size3)
    particles(ii)%E_k = rn_real_size3       !< E_k
    call gnu_rng_interval(v_interval,rn_real)
    particles(ii)%vpar_m = rn_real          !< vpar_m
    call gnu_rng_interval(Bnorm_interval,rn_real)
    particles(ii)%Bn_k = rn_real            !< Bn_k
    call gnu_rng_interval(v_interval,rn_real)
    particles(ii)%vpar = rn_real            !< vpar
    call gnu_rng_interval(v_interval,rn_real)
    particles(ii)%mu = rn_real              !< mu
    call gnu_rng_interval(q_interval,rn_integer)
    particles(ii)%q = int(rn_integer,kind=1)
  enddo
  !$omp end do
  !$omp end parallel
end subroutine fill_particle_gc_Qin

!> fill up particle_kinetic with random numbers
subroutine fill_particle_kinetic(n_particles,particles,rank_in)
  use mod_particle_types, only: particle_kinetic
  use mod_gnu_rng,        only: gnu_rng_interval
  use mod_gnu_rng,        only: set_seed_sys_time
  !$ use omp_lib
  implicit none
  !> inputs
  integer,intent(in) :: n_particles
  integer,intent(in),optional :: rank_in
  type(particle_kinetic),dimension(n_particles),intent(inout) :: particles
  !> variables
  integer :: ii,rn_integer,rank
  !$ integer :: thread_id
  real*8,dimension(3) :: rn_real_size3
  rank = 1
  if(present(rank_in)) rank=rank_in
  !$omp parallel default(private) firstprivate(n_particles) shared(particles)
  thread_id = 0
  !$ thread_id = omp_get_thread_num()
  call set_seed_sys_time(rng_seed_interval,rank,thread_id)
  !$omp do
  do ii=1,n_particles
    call gnu_rng_interval(3,vp3d_lowbnd,vp3d_uppbnd,rn_real_size3)
    particles(ii)%v = rn_real_size3
    call gnu_rng_interval(q_interval,rn_integer)
    particles(ii)%q = int(rn_integer,kind=1)
  enddo
  !$omp end do
  !$omp end parallel
end subroutine fill_particle_kinetic

!> fill up particle_kinetic_leapfrog with random numbers
subroutine fill_particle_kinetic_leapfrog(n_particles,particles,rank_in)
  use mod_particle_types, only: particle_kinetic_leapfrog
  use mod_gnu_rng,        only: gnu_rng_interval
  use mod_gnu_rng,        only: set_seed_sys_time
  !$ use omp_lib
  implicit none
  !> inputs
  integer,intent(in) :: n_particles
  integer,intent(in),optional :: rank_in
  type(particle_kinetic_leapfrog),dimension(n_particles),intent(inout) :: particles
  !> variables
  integer :: ii,rn_integer,rank
  !$ integer :: thread_id
  real*8,dimension(3) :: rn_real_size3
  rank = 1
  if(present(rank_in)) rank=rank_in
  !$omp parallel default(private) firstprivate(n_particles) shared(particles)
  thread_id = 0
  !$ thread_id = omp_get_thread_num()
  call set_seed_sys_time(rng_seed_interval,rank,thread_id)
  !$omp do
  do ii=1,n_particles
    call gnu_rng_interval(3,vp3d_lowbnd,vp3d_uppbnd,rn_real_size3)
    particles(ii)%v = rn_real_size3
    call gnu_rng_interval(q_interval,rn_integer)
    particles(ii)%q = int(rn_integer,kind=1)
  enddo
  !$omp end do
  !$omp end parallel
end subroutine fill_particle_kinetic_leapfrog

!> fill up particle_kinetic_relativistic with random numbers
subroutine fill_particle_kinetic_relativistic(n_particles,particles,rank_in)
  use mod_particle_types, only: particle_kinetic_relativistic
  use mod_gnu_rng,        only: gnu_rng_interval
  use mod_gnu_rng,        only: set_seed_sys_time
  !$ use omp_lib
  implicit none
  !> inputs
  integer,intent(in) :: n_particles
  integer,intent(in),optional :: rank_in
  type(particle_kinetic_relativistic),dimension(n_particles),intent(inout) :: particles
  !> variables
  integer :: ii,rn_integer,rank
  !$ integer :: thread_id
  real*8,dimension(3) :: rn_real_size3
  rank = 1
  if(present(rank_in)) rank=rank_in
  !$omp parallel default(private) firstprivate(n_particles) shared(particles)
  thread_id = 0
  !$ thread_id = omp_get_thread_num()
  call set_seed_sys_time(rng_seed_interval,rank,thread_id)
  !$omp do
  do ii=1,n_particles
    call gnu_rng_interval(3,vp3d_lowbnd,vp3d_uppbnd,rn_real_size3)
    particles(ii)%p = rn_real_size3
    call gnu_rng_interval(q_interval,rn_integer)
    particles(ii)%q = int(rn_integer,kind=1)
  enddo
  !$omp end do
  !$omp end parallel
end subroutine fill_particle_kinetic_relativistic

!> fill up particle_gc_relativistic with random numbers
subroutine fill_particle_gc_relativistic(n_particles,particles,rank_in)
  use mod_particle_types, only: particle_gc_relativistic
  use mod_gnu_rng,        only: gnu_rng_interval
  use mod_gnu_rng,        only: set_seed_sys_time
  !$ use omp_lib
  implicit none
  !> inputs
  integer,intent(in) :: n_particles
  integer,intent(in),optional :: rank_in
  type(particle_gc_relativistic),dimension(n_particles),intent(inout) :: particles
  !> variables
  integer :: ii,rn_integer,rank
  !$ integer :: thread_id
  real*8,dimension(2) :: rn_real_size2
  rank = 1
  if(present(rank_in)) rank=rank_in
  !$omp parallel default(private) firstprivate(n_particles) shared(particles)
  thread_id = 0
  !$ thread_id = omp_get_thread_num()
  call set_seed_sys_time(rng_seed_interval,rank,thread_id)
  !$omp do
  do ii=1,n_particles
    call gnu_rng_interval(2,vp3d_lowbnd(1:2),vp3d_uppbnd(1:2),rn_real_size2)
    particles(ii)%p = rn_real_size2
    call gnu_rng_interval(q_interval,rn_integer)
    particles(ii)%q = int(rn_integer,kind=1)
  enddo
  !$omp end do
  !$omp end parallel
end subroutine fill_particle_gc_relativistic

!> fill up particle_kinetic_relativistic as runaway electrons
!> the particle position must be provided in cylindrical coord.
subroutine fill_particle_kinetic_relativistic_RE(n_particles,particles,rank_in)
  use mod_particle_types,        only: particle_kinetic_relativistic
  use mod_gnu_rng,               only: gnu_rng_interval
  use mod_gnu_rng,               only: set_seed_sys_time
  use mod_coordinate_transforms, only: vector_cylindrical_to_cartesian
  use mod_math_operators,        only: cross_product
  use constants,                 only: SPEED_OF_LIGHT
  !$ use omp_lib
  implicit none
  !> parameters
  real*8,parameter :: E0=5.1099895d5 !< electron rest energy in MV
  real*8,parameter :: mass=5.48579909065d-4 !< electron mass in AMU
  !> inputs
  integer,intent(in) :: n_particles
  integer,intent(in),optional :: rank_in
  type(particle_kinetic_relativistic),dimension(n_particles),intent(inout) :: particles
  !> variables
  integer :: ii,rank
  !$ integer :: thread_id
  real*8 :: p_tot
  real*8,dimension(3) :: rn_real_size3,B_field,E_field,gradpsi,binorm
  rank = 1; if(present(rank_in)) rank=rank_in;
  !$omp parallel default(private) firstprivate(n_particles) shared(particles)
  thread_id = 0
  !$ thread_id = omp_get_thread_num()
  call set_seed_sys_time(rng_seed_interval,rank,thread_id) 
  !$omp do
  do ii=1,n_particles
    !> create orthonormal basis
    call compute_test_E_B_fields(particles(ii)%x,B_field,E_field); B_field = B_field/norm2(B_field) 
    call compute_test_gradpsi(particles(ii)%x,gradpsi); gradpsi = gradpsi/norm2(gradpsi)
    binorm = cross_product(B_field,gradpsi); binorm = binorm/norm2(binorm)
    !> compute particle momentum cartesian coord
    call gnu_rng_interval(3,EThetaChi_RE_lowbnd,EThetaChi_RE_uppbnd,rn_real_size3)
    p_tot = mass*SPEED_OF_LIGHT*sqrt(((rn_real_size3(1)/E0)+1.d0)**2.d0 - 1.d0)
    particles(ii)%p = p_tot*(B_field*cos(rn_real_size3(2))+sin(rn_real_size3(2))*(&
    cos(rn_real_size3(3))*gradpsi+sin(rn_real_size3(3))*binorm))
    particles(ii)%p = vector_cylindrical_to_cartesian(particles(ii)%x(3),particles(ii)%p)
    !> store electric charge
    particles(ii)%q = int(-1,kind=1)
  enddo 
  !$omp end do
  !$omp end parallel
end subroutine fill_particle_kinetic_relativistic_RE
 
!> fill up particle_kinetic_relativistic as runaway electrons
subroutine fill_particle_gc_relativistic_RE(n_particles,particles,rank_in)
  use mod_particle_types,        only: particle_gc_relativistic
  use mod_gnu_rng,               only: gnu_rng_interval
  use mod_gnu_rng,               only: set_seed_sys_time
  use constants,                 only: SPEED_OF_LIGHT
  !$ use omp_lib
  implicit none
  !> parameters
  real*8,parameter :: E0=5.1099895d5 !< electron rest energy in MV
  real*8,parameter :: mass=5.48579909065d-4 !< electron mass in AMU
  !> inputs
  integer,intent(in) :: n_particles
  integer,intent(in),optional :: rank_in
  type(particle_gc_relativistic),dimension(n_particles),intent(inout) :: particles
  !> variables
  integer :: ii,rn_integer,rank
  !$ integer :: thread_id
  real*8              :: p_tot,B_intensity
  real*8,dimension(2) :: rn_real_size2
  real*8,dimension(3) :: B_field,E_field
  rank = 1
  if(present(rank_in)) rank=rank_in
  !$omp parallel default(private) firstprivate(n_particles) shared(particles)
  thread_id = 0
  !$ thread_id = omp_get_thread_num()
  call set_seed_sys_time(rng_seed_interval,rank,thread_id)
  !$omp do
  do ii=1,n_particles
    !> computing the magnetic field
    call compute_test_E_B_fields(particles(ii)%x,B_field,E_field); 
    B_intensity = norm2(B_field); B_field = B_field/B_intensity;
    !> computing and storing the parallel momentum and magnetic moment
    call gnu_rng_interval(2,EThetaChi_RE_lowbnd(1:2),EThetaChi_RE_uppbnd(1:2),rn_real_size2)
    p_tot = mass*SPEED_OF_LIGHT*sqrt(((rn_real_size2(1)/E0)+1.d0)**2.d0 - 1.d0)
    particles(ii)%p(1) = p_tot*cos(rn_real_size2(2))
    particles(ii)%p(2) = (p_tot*p_tot-particles(ii)%p(1)*particles(ii)%p(1))/(2.d0*mass*B_intensity)
    !> store the electron charge    
    particles(ii)%q = int(-1,kind=1)
  enddo
  !$omp end do
  !$omp end parallel
end subroutine fill_particle_gc_relativistic_RE 

!>-------------------------------------------------------------
end module mod_particle_common_test_tools
