!> method used for testing the array handlers implemented 
!> in mod_array_handlers
module mod_array_handlers_test
use fruit
implicit none
private
public :: run_fruit_array_handlers

!> Variables and datatypes ------------------------------
integer,parameter             :: n_values=100
integer,parameter             :: array_size=1000
real*8,dimension(2),parameter :: value_int_1d_r8=(/-2.d2,1.d5/)
real*8,dimension(2),parameter :: array_int_1d_r8=(/-4.d2,3.d2/)
integer,dimension(n_values)   :: value_ids_1d
real*8,dimension(n_values)    :: values_1d_r8
real*8,dimension(array_size)  :: init_array_1d_r8,array_sol_1d_r8
!> Interfaces -------------------------------------------
contains
!> Fruit basket -----------------------------------------
!> running all set-up, test and tear-down procedures
subroutine run_fruit_array_handlers()
  implicit none
  write(*,'(/A)') "  ... setting-up: array handlers tests"
  call setup
  write(*,'(/A)') "  ... running: array handlers tests"
  call test_set_values_in_array_1d_r8
  write(*,'(/A)') "  ... tearing-down: array handlers tests"
end subroutine run_fruit_array_handlers

!> Set-up and tear-down ---------------------------------
!> set-up procedure
subroutine setup()
  use mod_gnu_rng,only: gnu_rng_interval
  implicit none
  !> generate a random set of ids
  call gnu_rng_interval(n_values,(/1,n_values/),value_ids_1d)
  call gnu_rng_interval(n_values,value_int_1d_r8,values_1d_r8)
  call gnu_rng_interval(array_size,array_int_1d_r8,init_array_1d_r8)
  array_sol_1d_r8 = init_array_1d_r8
end subroutine setup

!> Tests ------------------------------------------------
!> test the set_values_in_array for 1d r8 array
subroutine test_set_values_in_array_1d_r8()
  use mod_array_handlers, only: set_values_in_array
  implicit none
  !> variables
  integer :: ii
  real*8,dimension(array_size) :: test_array
  !> prepare test and solution arrays
  do ii=1,n_values
    array_sol_1d_r8(value_ids_1d(ii))=values_1d_r8(ii)
  enddo
  test_array = init_array_1d_r8
  !> test
  call set_values_in_array(n_values,array_size,value_ids_1d,&
  values_1d_r8,test_array)
  call assert_equals(test_array,array_sol_1d_r8,array_size,&
  "Error test set values in array 1d real8: array mismatch!")
  !> cleanup solution 
  array_sol_1d_r8 = init_array_1d_r8
end subroutine test_set_values_in_array_1d_r8

!> Tools ------------------------------------------------
!>-------------------------------------------------------
end module mod_array_handlers_test
