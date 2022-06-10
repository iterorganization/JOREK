!> mod_particle_types_test contains variables and 
!> procedure for testing the methods coded in
!> mod_particle_types
module mod_particle_types_test
use fruit
use mod_particle_types
use mod_particle_sim, only: particle_group
use mod_particle_common_test_tools, only: n_particle_types
implicit none

private
public :: run_fruit_particle_types

!> Variables ------------------------------------
integer,parameter :: n_particles=1000
real*8,parameter  :: survival_prob=6.25d-1
integer*1,dimension(n_particles,n_particle_types) :: particle_charge_list_sol
type(particle_group),dimension(n_particle_types)  :: groups_sol
!> Interfaces -----------------------------------

contains

!> Fruit basket ---------------------------------
!> fruit basket containing all set-up, test and
!> tear-down procedures
subroutine run_fruit_particle_types()
  implicit none
  write(*,'(/A)') "  ... setting-up: particle types tests"
  call setup
  write(*,'(/A)') "  ... running: particle types tests"
  call test_particle_copy
  call test_individual_particle_copy
  call test_particle_get_q
  write(*,'(/A)') "  ... tearing-up: particle types tests"
end subroutine run_fruit_particle_types

!> Set-up and tear-down -------------------------
!> set-up particle types test features
subroutine setup()
  use mod_particle_common_test_tools, only: fill_particles
  use mod_particle_common_test_tools, only: fill_groups
  use mod_particle_common_test_tools, only: obtain_particle_charges
  use mod_particle_common_test_tools, only: allocate_one_particle_list_type
  implicit none
  !> variables
  integer :: ifail
  
  !> allocate the particle lists
  ifail = 0
  call allocate_one_particle_list_type(n_particle_types,n_particles,groups_sol,ifail)
  call assert_true(ifail.eq.0,"Error particle_types test setup: particle list not allocated!")

  !> fill-up the group and particle base variables
  call fill_groups(n_particle_types,groups_sol)
  call fill_particles(n_particle_types,n_particles,groups_sol)
  
  !> get particle charges
  call obtain_particle_charges(n_particle_types,n_particles,particle_charge_list_sol,groups_sol)
end subroutine setup

!> Tests ----------------------------------------
!> test particle copy function
subroutine test_particle_copy()
  use mod_particle_sim, only: particle_group
  use mod_particle_assert_equal, only: assert_equal_particle
  use mod_particle_common_test_tools, only: allocate_one_particle_list_type
  implicit none
  !> variables
  type(particle_group),dimension(n_particle_types) :: group_particles
  integer :: ii,jj,ifail
  !> allocate particle lists
  ifail = 0
  call allocate_one_particle_list_type(n_particle_types,n_particles,group_particles,ifail)
  call assert_true(ifail.eq.0,"Error particle_types test copy: particle list not allocated!")
  !> copy particles
  !$omp parallel do default(shared) private(ii,jj) collapse(2)
  do jj=1,n_particle_types
    do ii=1,n_particles
      group_particles(jj)%particles(ii) = groups_sol(jj)%particles(ii)
    enddo
  enddo
  !$omp end parallel do
  !> compare particle list
  do ii=1,n_particle_types
    call assert_equal_particle(n_particles,group_particles(ii)%particles,groups_sol(ii)%particles)
  enddo
end subroutine test_particle_copy

!> test indiviual particle copy functions
!> it assumes that the standard particle
!> copy function works
subroutine test_individual_particle_copy()
  use mod_particle_assert_equal, only: assert_equal_particle
  implicit none
  type(particle_kinetic_leapfrog) :: p_leapfrog_in,p_leapfrog_out
  integer :: ii

  !> initialise particle in
  !> copy particle types assuming that the standard particle copy function works
  do ii=1,n_particle_types
    select type (p_in=>groups_sol(ii)%particles(1))
    type is (particle_kinetic_leapfrog)
      p_leapfrog_in = p_in
    end select
  enddo
  
  !> apply individual copy function
  call copy_particle_kinetic_leapfrog(p_leapfrog_in,p_leapfrog_out)

  !> test equality
  call assert_equal_particle(p_leapfrog_in,p_leapfrog_out)
end subroutine test_individual_particle_copy

!> test return charge function
subroutine test_particle_get_q()
  implicit none
  !> variables
  integer :: ii,jj
  integer*1,dimension(n_particles,n_particle_types) :: q_array
  !> extract particle charge and test
  !$omp parallel do default(shared) private(ii,jj) collapse(2)
  do jj=1,n_particle_types
    do ii=1,n_particles
      q_array(ii,jj) = particle_get_q(groups_sol(jj)%particles(ii))
    enddo
  enddo
  !$omp end parallel do
  call assert_equals(int(q_array),int(particle_charge_list_sol),n_particles,&
  n_particle_types,"Errpr particle_type test get q: particle charge mismatch!")
end subroutine test_particle_get_q

!>-----------------------------------------------

end module mod_particle_types_test
