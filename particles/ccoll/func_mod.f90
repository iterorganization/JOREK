!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
!
!> @brief
!> A module for defining function pointer types and interfaces
!<
module func_mod
  implicit none
  
  !> Interface for passing double precision 1D functions as arguments.
  !> To use these, declare
  !>
  !> procedure(func_real8_1d), pointer :: my_function_pointer => null()
  !> 
  !> and then set the my_function_pointer to point to the function
  !> by declaring
  !>
  !> my_function_pointer => my_function
  !>
  !<
  abstract interface
     function func_real8_1D (x)
       import
       real*8 :: func_real8_1d
       real*8, intent (in) :: x
     end function func_real8_1D
  end interface

  public :: func_real8_1d

contains 

end module func_mod
