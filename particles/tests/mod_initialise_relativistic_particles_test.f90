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
