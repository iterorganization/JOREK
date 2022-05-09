module mod_newton
!    inexact newton method for adaptive time stepping
!    Details may be found on p. 24 in: 
!    Franck et al., Energy conservation and numerical stability for the 
!    reduced MHD models of the non-linear JOREK code, 2014, arXiv.1408:2099v3


  use construct_matrix_mod, only: construct_matrix
  use mod_clock
  use global_distributed_matrix, only: A_glob, RHS_glob, ndof_glob, nz_glob, deltas
#ifdef USE_BICGSTAB
  use mod_bicgstab, only: bicgstab_driver, bicgstab_finalize
#else
  use mod_gmres, only: gmres_driver
  use data_structure
#endif
  implicit none
!> solve Jk delta_k_1 = -J_n delta_k_n+rhs_glob using inexact newton
!! x
  private
  public :: inexact_newton
  contains  
  subroutine inexact_newton(n_step_newton, element_list, node_list)
    implicit none
    !integer(kind=C_INT), pointer, intent(in) :: nz_glob(:) 
    real*8, allocatable        :: A_tmp(:)    ! to copy A_glob at step n
    real*8, allocatable        :: rhs_tmp(:)  ! to copy rhs_glob at step n
    type(clcktype)             :: t0,t1
    real*8                     :: tsecond
    integer(kind=4)                    :: n_step_newton
    type(type_element_list),intent(in)    :: element_list
    type(type_node_list),intent(in)       :: node_list
    !n_step_newton = 5

    call clck_time(t0)
    allocate(A_tmp(1:nz_glob))
    allocate(rhs_tmp(1:ndof_glob))
    A_tmp(1:nz_glob) = A_glob(1:nz_glob)
    rhs_tmp(1:ndof_glob) = rhs_glob(1:ndof_glob)

    call update_values(n_step_newton,element_list,node_list,deltas)         ! add solution to node values
    call update_deltas(n_step_newton,node_list)


    write(*,*)
    write(*,*) '>>>> SAVE COPY nz_glob, ndof_glob, my_id = ',nz_glob, ndof_glob, n_step_newton
    write(*,*)
    !open(1, file='A_tmp.dat', status='new') ! can't override file, need to check
    !write(1,*) A_tmp
    !close(1)
    call clck_time_barrier(t1)
    call clck_ldiff(t0,t1,tsecond)
    write(*,*) 'Elapsed time save copy: ', tsecond
    deallocate(A_tmp)
    deallocate(rhs_tmp)
  end subroutine inexact_newton

end module mod_newton
