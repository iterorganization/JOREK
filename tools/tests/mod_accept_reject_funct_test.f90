!> mod_accept_reject_funct_test contains variables
!> and procedures used for testing the accept-rejection
!> functions implemented in mod_accept_reject_funct
module mod_accept_reject_funct_test
use fruit
implicit none
private 
public :: run_fruit_accept_reject_funct

!> Variables and datatypes ------------------------------
integer,parameter            :: n_trials=4
real*8,dimension(n_trials,2) :: positive_interval=reshape(&
  (/-5.d1,-3.d0,5.d2,1.d-1,-1.d0,6.d1,1.d3,5.d-1/),(/n_trials,2/))
real*8,dimension(n_trials,2) :: negative_interval=reshape(&
  (/-5.d1,6.d0,5.d-1,4.d1,-1.d2,-2.d0,7.d-2,1.1d0/),(/n_trials,2/))
real*8,dimension(n_trials)   :: rands 
!> Interfaces -------------------------------------------
contains
!> Fruit basket -----------------------------------------
!> initialise, run and tear down all test features
subroutine run_fruit_accept_reject_funct()
  implicit none
  write(*,'(/A)') "  ... setting-up: accept reject funct tests"
  call setup()
  write(*,'(/A)') "  ... running: accept reject funct tests"
  call test_accept_lower_value_rand
  write(*,'(/A)') "  ... tearing-down: accept reject funct tests"
  call teardown()
end subroutine run_fruit_accept_reject_funct

!> Set-up and tear-down ---------------------------------
!> set-up test features
subroutine setup()
  implicit none
  call random_number(rands)
end subroutine setup

!> tear-down test features
subroutine teardown()
  implicit none
  rands = 0.d0;
end subroutine teardown

!> Tests ------------------------------------------------
!> test acceptance if value smaller than rand
subroutine test_accept_lower_value_rand()
  use mod_accept_reject_funct, only: accept_lower_values_rand
  implicit none
  !> variables
  real*8 :: rand_success,rand_fail
  real*8,dimension(n_trials) :: values
  logical :: test
  !> initialisation
  rand_success=maxval(rands); rand_success=max(1.d0,1.01*rand_success);
  rand_fail=minval(rands); rand_fail=max(0.d0,0.99*rand_fail);
  !> test positive interval
  values = positive_interval(:,1)+(positive_interval(:,2)-&
  positive_interval(:,1))*rands
  test = accept_lower_values_rand(n_trials,values,positive_interval,rand_success)
  call assert_true(test,"Error accept-reject accept lower value: positive interval not success!")
  test = .not.accept_lower_values_rand(n_trials,values,positive_interval,rand_fail)
  call assert_true(test,"Error accept-reject accept lower value: positive interval not fail!")
  !> test negative interval
  values = negative_interval(:,1)+(negative_interval(:,2)-&
  negative_interval(:,1))*rands
  test = accept_lower_values_rand(n_trials,values,negative_interval,rand_success)
  call assert_true(test,"Error accept-reject accept lower value: negative interval not success!")
  test = .not.accept_lower_values_rand(n_trials,values,negative_interval,rand_fail)
  call assert_true(test,"Error accept-reject accept lower value: negative interval not fail!")
end subroutine test_accept_lower_value_rand

!> Tools ------------------------------------------------
!>-------------------------------------------------------
end module mod_accept_reject_funct_test
