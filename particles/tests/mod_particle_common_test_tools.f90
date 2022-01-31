!> module mod_particle_type_test_tools contains variables and
!> procedure used for initialising and finalising particle data
!> for unit testing
module mod_particle_common_test_tools
implicit none

private
public :: n_particle_types
public :: q_interval,i_elm_interval,i_life_interval,sim_time_interval
public :: t_birth_interval,st_interval,mass_interval,v_interval
public :: Ekin_interval,mu_interval,Bnorm_interval,weight_interval
public :: x_lowbnd,x_uppbnd,vp3d_lowbnd,vp3d_uppbnd,ABE_lowbnd,ABE_uppbnd
public :: fill_particles,invalidate_particles,obtain_active_particle_ids
public :: fill_groups,fill_particle_base,fill_particle_fieldline
public :: fill_particle_gc,fill_particle_gc_vpar,fill_particle_gc_Qin
public :: fill_particle_kinetic,fill_particle_kinetic_leapfrog
public :: fill_particle_kinetic_relativistic,fill_particle_gc_relativistic
public :: obtain_particle_charges,allocate_one_particle_list_type
public :: copy_group_fieldline_B_hat_prev

!> Variables --------------------------------------------------
integer,parameter :: n_particle_types=8
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

!> Interfaces -------------------------------------------------
interface fill_groups
  module procedure fill_groups_seq
  module procedure fill_groups_mpi
end interface fill_groups
contains
!> Procedures -------------------------------------------------
!> allocate one particle list per type
subroutine allocate_one_particle_list_type(n_groups,n_particles,groups,ifail)
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
    write(*,'(/A)') "allocate one particle list per type failed!"
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
end subroutine allocate_one_particle_list_type

!> obtain charges from all particles in a simulation 
subroutine obtain_particle_charges(n_groups,n_particles,charge_list,groups)
  use mod_particle_sim,   only: particle_group
  use mod_particle_types, only: particle_kinetic,particle_kinetic_leapfrog
  use mod_particle_types, only: particle_gc
  use mod_particle_types, only: particle_kinetic_relativistic
  use mod_particle_types, only: particle_gc_relativistic
  use mod_particle_types, only: particle_gc_vpar,particle_gc_Qin 
  implicit none
  !> inputs
  integer,intent(in) :: n_groups,n_particles
  !> inputs-outputs
  type(particle_group),dimension(n_groups),intent(inout) :: groups
  !> outputs
  integer*1,dimension(n_particles,n_groups),intent(out) :: charge_list
  !> variables
  integer :: ii,jj

  !> initialisation
  charge_list = 0
  !> extract charges
  !$omp parallel do default(shared) firstprivate(n_groups,n_particles) &
  !$omp private(ii,jj) collapse(2)
  do jj=1,n_groups
    do ii=1,n_particles
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
  enddo
  !$omp end parallel do
end subroutine obtain_particle_charges

!> extract the index of valid particles in a simulatiion
subroutine obtain_active_particle_ids(n_groups,n_particles,&
active_particle_ids,groups)
  use mod_particle_sim, only: particle_group
  implicit none
  !> inputs
  integer,intent(in) :: n_groups,n_particles
  !> inputs-outputs
  type(particle_group),dimension(n_groups),intent(inout) :: groups
  !> outputs:
  integer,dimension(n_particles,n_groups),intent(out) :: active_particle_ids
  !> variables
  integer :: ii,jj,counter
  do jj=1,n_groups
    counter = 0
    do ii=1,n_particles
      if(groups(jj)%particles(ii)%i_elm.gt.0) then
        counter = counter + 1
        active_particle_ids(counter,jj) = ii
      endif
    enddo
  enddo
end subroutine obtain_active_particle_ids

!> invalidate some of the particles in the particle groups of a
!> simulation as a function of a survival probability
subroutine invalidate_particles(n_groups,n_particles,&
  survival_prob_in,n_active_particles,groups,rank_in)
  !$ use omp_lib
  use mod_particle_sim, only: particle_group
  use mod_gnu_rng, only: gnu_rng_interval
  use mod_gnu_rng, only: set_seed_sys_time 
  implicit none
  !> inputs
  integer,intent(in) :: n_groups,n_particles
  real*8,intent(in)  :: survival_prob_in
  integer,intent(in),optional :: rank_in
  !> inputs-outputs
  type(particle_group),dimension(n_groups),intent(inout) :: groups
  !> outputs
  integer,dimension(n_groups),intent(out) :: n_active_particles
  !> variables
  integer :: ii,jj,rank,thread_id
  real*8  :: survival_prob
  real*8,dimension(n_particles,n_groups) :: survival_array
  !> initialisation
  rank = 1; if(present(rank_in)) rank=rank_in;
  survival_prob = abs(survival_prob_in)
  if(abs(survival_prob).gt.1.d0) survival_prob = survival_prob - floor(survival_prob)
  n_active_particles = 0
  !$omp parallel default(shared) private(thread_id)
  !$ thread_id = omp_get_thread_num()
  call set_seed_sys_time(rng_seed_interval,rank,thread_id)
  !$omp end parallel
  
  !> generate survival probability
  call gnu_rng_interval(n_particles,n_groups,(/0.d0,1.d0/),survival_array)
  do ii=1,n_groups
    n_active_particles(ii) = count(survival_array(:,ii).ge.survival_prob)
  enddo
  !> invalidate particles
  !$omp parallel do default(shared) firstprivate(n_particles,n_groups,survival_prob) &
  !$omp private(ii,jj) collapse(2)
  do jj=1,n_groups
    do ii=1,n_particles
      if(survival_array(ii,jj).lt.survival_prob) groups(jj)%particles(ii)%i_elm = 0
    enddo
  enddo 
  !$omp end parallel do
end subroutine invalidate_particles

!> copy fieldlines B_hat between two simulations 
!> used for IO because it is not stored in hdf5
subroutine copy_group_fieldline_B_hat_prev(n_groups,n_particles,groups_in,groups_out)
  use mod_particle_types, only: particle_fieldline
  use mod_particle_sim,   only: particle_group
  implicit none
  integer,intent(in) :: n_groups,n_particles
  type(particle_group),dimension(n_groups),intent(in)    :: groups_in
  type(particle_group),dimension(n_groups),intent(inout) :: groups_out
  integer :: ii,jj
  !$omp parallel do default(private) firstprivate(n_groups,n_particles) &
  !$omp shared(groups_in,groups_out) collapse(2)
  do ii=1,n_groups
    do jj=1,n_particles
      select type (p_out=>groups_out(ii)%particles(jj))
      type is (particle_fieldline)
        select type (p_in=>groups_in(ii)%particles(jj))
          type is (particle_fieldline)
          p_out%B_hat_prev = p_in%B_hat_prev
        end select
      end select
    enddo
  enddo
  !$omp end parallel do
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

!> fills particle list of each groups with random data
subroutine fill_particles(n_groups,n_particles,groups,rank_in)
  use mod_particle_sim,   only: particle_group
  use mod_particle_types, only: particle_kinetic,particle_kinetic_leapfrog
  use mod_particle_types, only: particle_gc,particle_fieldline
  use mod_particle_types, only: particle_kinetic_relativistic
  use mod_particle_types, only: particle_gc_relativistic
  use mod_particle_types, only: particle_gc_vpar,particle_gc_Qin
  implicit none
  integer,intent(in) :: n_groups,n_particles
  type(particle_group),dimension(n_groups),intent(inout) :: groups
  integer,intent(in),optional :: rank_in
  integer :: rank,jj
  rank = 1
  if(present(rank_in)) rank = rank_in
  !> fill particle basic type
  call fill_particle_base(n_groups,n_particles,groups,rank)
  !> fill particle specific types
  do jj=1,n_groups
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
subroutine fill_particle_base(n_groups,n_particles,groups,rank_in)
  use mod_particle_sim, only: particle_group
  use mod_gnu_rng, only: gnu_rng_interval
  use mod_gnu_rng, only: set_seed_sys_time
  !$ use omp_lib
  implicit none
  !> inputs
  integer,intent(in) :: n_groups,n_particles
  integer,intent(in),optional :: rank_in
  !> inputs-outputs:
  type(particle_group),dimension(n_groups),intent(inout) :: groups
  !> variables
  integer :: ii,jj,rn_integer,rank
  !$ integer :: thread_id
  real*8              :: rn_real
  real*8,dimension(2) :: rn_real_size2
  real*8,dimension(3) :: rn_real_size3
  rank = 1
  if(present(rank_in)) rank=rank_in
  thread_id = 0
  !> fill-up the particle_base variables for all particles and all groups
  !$omp parallel default(shared) firstprivate(n_groups,n_particles) &
  !$omp private(ii,jj,rank,rn_integer,rn_real,rn_real_size2,rn_real_size3,thread_id)
  !$ thread_id = omp_get_thread_num()
  call set_seed_sys_time(rng_seed_interval,rank,thread_id)
  !$omp do collapse(2)
  do jj=1,n_groups
    do ii=1,n_particles
      call gnu_rng_interval(weight_interval,rn_real)
      call gnu_rng_interval(2,st_interval,rn_real_size2)
      call gnu_rng_interval(3,x_lowbnd,x_uppbnd,rn_real_size3)
      groups(jj)%particles(ii)%x       = rn_real_size3
      groups(jj)%particles(ii)%st      = rn_real_size2
      groups(jj)%particles(ii)%weight  = rn_real
      call gnu_rng_interval(t_birth_interval,rn_real)
      groups(jj)%particles(ii)%t_birth = real(rn_real,kind=4)
      call gnu_rng_interval(i_elm_interval,rn_integer)
      groups(jj)%particles(ii)%i_elm   = rn_integer
      call gnu_rng_interval(i_life_interval,rn_integer)
      groups(jj)%particles(ii)%i_life  = rn_integer 
    enddo
  enddo
  !$omp end do
  !$omp end parallel
end subroutine fill_particle_base

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

!>-------------------------------------------------------------
end module mod_particle_common_test_tools
