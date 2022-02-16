!> test the initialisation routines for relativistic particles
module mod_initialise_relativistic_particles_test
use fruit
implicit none
private
public :: run_fruit_initialise_relativistic_particles

!> Variables and datatypes ------------------------------------
!> Interfaces--------------------------------------------------
contains
!> Fruit basket -----------------------------------------------
!> fruit test basket: set-up, run and tear-down tests
subroutine run_fruit_initialise_relativistic_particles()
  implicit none
  write(*,'(/A)') "  ... setting-up: initialise relativistic particles tests"
  write(*,'(/A)') "  ... running: initialise relativistic particles tests"
  call test_dummy
  write(*,'(/A)') "  ... tearing-down: initialise relativistic particles tests"
end subroutine run_fruit_initialise_relativistic_particles

!> Set-up and tear-down ---------------------------------------
!> Tests ------------------------------------------------------
subroutine test_dummy()
  use mod_initialise_relativistic_particles
  implicit none
  write(*,'(/A)') "initialise relativistic particle dummy test"
end subroutine test_dummy
!> Tools ------------------------------------------------------
!>-------------------------------------------------------------
end module mod_initialise_relativistic_particles_test
