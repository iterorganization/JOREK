!> mod_array_handlers is a collection of procedures used 
!> for manipulating arrays
module mod_array_handlers
implicit none
private
public :: set_values_in_array

!> Variables and datatypes -------------------------------
!> Interfaces --------------------------------------------
interface set_values_in_array
  module procedure set_values_in_array_1d_r8
end interface set_values_in_array
contains
!> Procedures --------------------------------------------
!> set values in a 1d real8 array to specific indices.
!> WARNING: no consistency check is performed on the values
!> and array sizes
!> inputs:
!>   n_values:         (integer) number of values to set
!>   array_size:       (integer) size of the array
!>   value_to_set_ids: (integer)(n_values) indices of the values
!>                     to set in array
!>   values_to_set:    (real8)(n_values) values to set in array
!>   array:            (real8)(array_size) array with values to
!>                     be replaced
!> outputs:
!>   array:            (real8)(array_size) array with replaced values
subroutine set_values_in_array_1d_r8(n_values_to_set,array_size,&
value_to_set_ids,values_to_set,array)
  implicit none
  !> inputs:
  integer,intent(in) :: n_values_to_set,array_size
  integer,dimension(n_values_to_set),intent(in) :: value_to_set_ids
  real*8,dimension(n_values_to_set),intent(in)  :: values_to_set
  !> inputs-outputs:
  real*8,dimension(array_size),intent(inout) :: array
  !> variables:
  integer :: ii
  !> copy values
  do ii=1,n_values_to_set
    array(value_to_set_ids(ii)) = values_to_set(ii)
  enddo
end subroutine set_values_in_array_1d_r8
!> -------------------------------------------------------
end module mod_array_handlers

