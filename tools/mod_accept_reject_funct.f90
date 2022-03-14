!> module containing functions used for
!> acceptance-rejection procedure
module mod_accept_reject_funct
implicit none
private
public :: accept_lower_values_rand
!> Variables and datatypes ----------------------------------------------
!> Interfaces -----------------------------------------------------------
contains
!> Procedures -----------------------------------------------------------
!> accept a test if the computed value (normalised to be withing [0,1])
!> is smaller or equal to a random number in [0,1]
function accept_lower_values_rand(n_values,values,intervals,rand) &
result(success)
  implicit none
  !> inputs:
  integer,intent(in) :: n_values
  real*8,intent(in)  :: rand
  real*8,dimension(n_values),intent(in)   :: values
  real*8,dimension(n_values,2),intent(in) :: intervals
  !> outputs
  logical :: success
  !> check if value accepted
  success = all(((intervals(:,2)-intervals(:,1))*(values-intervals(:,1))).le.&
  (((intervals(:,2)-intervals(:,1))**2)*rand))
end function accept_lower_values_rand
!>-----------------------------------------------------------------------
end module mod_accept_reject_funct
