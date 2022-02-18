!> test the initialisation routines for relativistic particles
module mod_initialise_relativistic_particles_test
use fruit
use mod_particle_sim,      only: particle_group
use mod_particle_types,    only: particle_kinetic_relativistic_id
use mod_particle_types,    only: particle_gc_relativistic_id
use mod_particle_types,    only: particle_kinetic_relativistic
use mod_particle_types,    only: particle_gc_relativistic
use mod_fields_analytical, only: fields_analytical
implicit none

private
public :: run_fruit_initialise_relativistic_particles

!> Variables and datatypes ------------------------------------
integer,parameter  :: n_groups=2
integer,dimension(n_groups),parameter  :: n_particles=(/123,234/)
integer,dimension(n_groups),parameter  :: p_types_sol=(/&
                      particle_kinetic_relativistic_id,&
                      particle_gc_relativistic_id/)
type(particle_group),dimension(n_groups) :: groups_sol
type(fields_analytical)                  :: fields_sol
!> Interfaces--------------------------------------------------
contains
!> Fruit basket -----------------------------------------------
!> fruit test basket: set-up, run and tear-down tests
subroutine run_fruit_initialise_relativistic_particles()
  implicit none
  write(*,'(/A)') "  ... setting-up: initialise relativistic particles tests"
  call setup
  write(*,'(/A)') "  ... running: initialise relativistic particles tests"
  call test_sampling_cartesian_p_kinetic_relativistic
  write(*,'(/A)') "  ... tearing-down: initialise relativistic particles tests"
  call teardown
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
  !> initialise the group masses as runaway
  call fill_mass_RE(n_groups,groups_sol)
  !> allocate particle lists
  call allocate_one_particle_list_type(n_groups,n_particles,&
  p_types_sol,groups_sol,ifail)
  !> initialise particle fields
  call gnu_rng_interval(2,RZ0_lowbnd,RZ0_uppbnd,real_param(1:2))
  call gnu_rng_interval(2,BE0_lowbnd,BE0_uppbnd,real_param(3:4))
  call fields_sol%init_fields(0,4,int_param,real_param)
end subroutine setup

!> tear-down test feature
subroutine teardown()
  implicit none
  !> variables
  integer :: ii
  !> clean particle group
  do ii=1,n_groups; groups_sol(ii)%mass=0.d0; enddo
  !> cleanup fields
  call fields_sol%deallocate_fields()
end subroutine teardown

!> Tests ------------------------------------------------------
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
      call assert_true(all(success),"Errror sampling cart. p kinetic relat.: momenta not in bound!")
      deallocate(success)
    end select
  enddo
end subroutine test_sampling_cartesian_p_kinetic_relativistic

!> dummy test procedure
subroutine test_dummy()
  use mod_initialise_relativistic_particles
  use mod_fields_analytical, only: fields_analytical
  implicit none
  type(fields_analytical) :: fields
  real*8 :: psi,U
  real*8,dimension(3) :: B,E
  write(*,'(/A)') "initialise relativistic particle dummy test"
  call fields%calc_EBPsiU(0.d0,0,(/0.d0,0.d0/),0.d0,B,E,psi,U)
  write(*,*) "fields analytical B: ",B
  write(*,*) "fields analytical E: ",E
  write(*,*) "fields analytical psi: ",psi
  write(*,*) "fields analytical U: ",U
end subroutine test_dummy

!> Tools ------------------------------------------------------
!>-------------------------------------------------------------
end module mod_initialise_relativistic_particles_test
