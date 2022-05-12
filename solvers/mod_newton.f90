module mod_newton
!    Inexact newton method, details may be found on p. 24 in: 
!    Franck et al., Energy conservation and numerical stability for the 
!    reduced MHD models of the non-linear JOREK code, 2014, arXiv:1408.2099v3


  use construct_matrix_mod, only: construct_matrix
  use mod_clock
!  use global_distributed_matrix, only: A_glob, RHS_glob, ndof_glob, nz_glob, deltas
#ifdef USE_BICGSTAB
  use mod_bicgstab, only: bicgstab_driver, bicgstab_finalize
#else
  use mod_gmres, only: gmres_driver
#endif
  use data_structure
  implicit none
!> solve Jk delta_k_1 = -J_n delta_k_n+rhs_glob using inexact newton
!! x
  private
  public :: inexact_newton
  contains  
  subroutine inexact_newton(val, x, b, max_it, tol, comm_glob, comm_n, comm_master,                             &
    ! end of additional arguments of bicgstab_driver 
                            iter_gmres,                                                                         &
    ! end of additional arguments of gmres_driver         
                            element_list, node_list,                                                            &
    ! end of additional arguments
                            my_id, MPI_COMM_N, my_id_n, MPI_COMM_MASTER, my_id_master, local_elms, n_local_elms,& 
                            index_min, index_max, xpoint2, xcase2, R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint, &
                            Z_xpoint, psi_xpoint,i_tor_min, i_tor_max, n, nz, ndof, n_matrix_block_size, A_mat, &
                            rhs, irn, jcn, ijA_index, ijA_size, irn_jcn, harmonic_matrix)
    ! end of arguments of construct_matrix


    ! inexact_newton(arguments of bicgstab_driver,  &
    !                arguments of gmres_driver,     &
    !                additional arguments,          & 
    !                arguments of construct_matrix)                   



    implicit none
    !--------------------- input variables -----------------------------------------
    !--- definitions from bicgstab_driver, gmres_driver
    integer,                      intent(inout)       :: my_id, my_id_n, my_id_master
    !integer(kind=C_INT), pointer, intent(in)       :: irn(:), jcn(:)
    integer,                      intent(inout)       :: MPI_COMM_N, MPI_COMM_MASTER
    real(kind=C_DOUBLE)                            :: tol   
    type(type_element_list)                        :: element_list
    type(type_node_list)                           :: node_list
    real(kind=C_DOUBLE), pointer, intent(in)       :: val(:)
    real(kind=C_DOUBLE), allocatable               :: b(:)
    real(kind=C_DOUBLE), allocatable               :: x(:)
    integer, intent(in)                            :: comm_glob, comm_n, comm_master
    integer, intent(inout)                         :: max_it, iter_gmres
    !--- definitions from module global_distributed_matrix
    !real*8,                allocatable, target     :: A_glob(:)    
    !real*8,                allocatable, target     :: rhs_glob(:)  
    !integer(kind=int_all), allocatable, target     :: irn_glob(:)  
    !integer(kind=int_all), allocatable, target     :: jcn_glob(:)  
    !real*8,                allocatable             :: deltas(:)    
    !integer(kind=int_all)                          :: ndof_glob, n_glob, nz_glob
    !--- definitions subroutine construct_matrix
    integer,               intent(inout)              :: local_elms(*)
    integer,               intent(inout)              :: n_local_elms
    integer,               intent(inout)              :: index_min
    integer,               intent(inout)              :: index_max
    integer,               intent(inout)              :: xcase2
    real*8,                intent(inout)              :: R_axis
    real*8,                intent(inout)              :: Z_axis
    real*8,                intent(inout)              :: psi_axis
    real*8,                intent(inout)              :: psi_bnd
    real*8,                intent(inout)              :: R_xpoint(2)
    real*8,                intent(inout)              :: Z_xpoint(2)
    real*8,                intent(inout)              :: psi_xpoint(2)
    logical,               intent(inout)              :: xpoint2
    integer,               intent(in)              :: i_tor_min
    integer,               intent(in)              :: i_tor_max
    integer(kind=int_all), intent(inout)              :: n, nz, ndof
    integer(kind=int_all), intent(inout)              :: n_matrix_block_size
    logical,               intent(in)              :: harmonic_matrix
    real*8,                intent(inout), allocatable :: A_mat(:)
    real*8,                intent(inout), allocatable :: rhs(:)
    integer(kind=int_all), intent(inout), allocatable :: irn(:)
    integer(kind=int_all), intent(inout), allocatable :: jcn(:)   
    integer(kind=int_all), intent(inout),    allocatable :: ijA_index(:,:), ijA_size(:), irn_jcn(:,:)
    !--------------------- end of input variables ----------------------------------


    !--------------------- routine variables ---------------------------------------
    real*8, allocatable                            :: A_n(:)          ! to copy A_glob at step n
    real*8, allocatable                            :: rhs_n(:)        ! to copy rhs_glob at step n
    real*8, allocatable                            :: delta_k_n(:)    ! U(k)-U(n)
    real*8, allocatable                            :: delta_k(:)      ! initial guess and sol
    real*8, allocatable                            :: rhs_k(:)
    type(clcktype)                                 :: t0,t1
    real*8                                         :: tsecond
    type(type_element_list)                        :: element_list_temp
    type(type_node_list)                           :: node_list_temp
    integer                                        :: newton_iter
    !n_step_newton = 5
    !--------------------- end of routine variables --------------------------------

    call clck_time(t0)
    allocate(A_n(1:nz),rhs_n(1:ndof),rhs_k(1:ndof),delta_k_n(1:ndof),delta_k(1:ndof))
  
    delta_k_n            = 0.0d0
    delta_k(1:ndof)      = x(1:ndof)
    A_n(1:nz)            = val(1:nz)
    rhs_n(1:ndof)        = b(1:ndof)
    element_list_temp    = element_list
    node_list_temp       = node_list 
    call clck_time_barrier(t1)
    call clck_ldiff(t0,t1,tsecond)
    !write(*,*) 'Elapsed time copy element_list, node_list, A_glob, rhs_glob', tsecond
    write(*,*)   '>>>>>>>>>>>>>>>>>>>>>> START NEWTON LOOP <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
    ! start loop
    newton_loop: do newton_iter = 1, 3
      write(*,'(a,i2.2,a)') '###################### NEWTON ITERATION STEP ', newton_iter,' ##########################'
      write(*,*)
      !write(*,*) 'Newton loop, newton_iter=',newton_iter

      delta_k_n(1:ndof) = delta_k_n(1:ndof) + delta_k(1:ndof)
      write(*,*) 'Set delta_k_n=delta_k_n+delta_k'

      call update_values(my_id,element_list,node_list,delta_k_n)
      call update_deltas(my_id,node_list)
      write(*,*) 'update_values, update_deltas'

      !call construct_matrix(my_id, MPI_COMM_N, my_id_n, MPI_COMM_MASTER, my_id_master, local_elms, n_local_elms, index_min, index_max,& 
      !                      xpoint2, xcase2, R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint, Z_xpoint, psi_xpoint,i_tor_min, i_tor_max,  &
      !                      n, nz, ndof, n_matrix_block_size, A_mat, rhs, irn, jcn, ijA_index, ijA_size, irn_jcn, harmonic_matrix)
      write(*,*) 'construct_matrix'

      rhs_k = rhs_n ! -A_n.dU, maybe it has to '+'
      !write(*,*) 'newton_iter = ', newton_iter
#ifdef USE_BICGSTAB
      call bicgstab_driver(irn, jcn, val, delta_k, rhs_k, max_it, tol, comm_glob, comm_n, comm_master)
#else
      call gmres_driver(my_id,my_id_n,MPI_COMM_N,MPI_COMM_MASTER,iter_gmres)
#endif
    enddo newton_loop
    write(*,*) '>>>>>>>>>>>>>>>>>>>>>> END NEWTON LOOP <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
    element_list    = element_list_temp
    node_list       = node_list_temp 

    !call update_values(my_id,element_list,node_list,x)  ! add solution to node values
    !call update_deltas(my_id,node_list)
    !write(*,*)
    !write(*,*) '>>>> SAVE COPY nz_glob, ndof_glob, my_id = ',nz, ndof, my_id
    !write(*,*)
    !call clck_time(t0)
    !open(1, file='A_n.dat', status='replace') 
    !write(1,*) A_n
    !close(1)
    !call clck_time_barrier(t1)
    !call clck_ldiff(t0,t1,tsecond)
    !write(*,*) 'Elapsed time write A_n to file: ', tsecond
    deallocate(A_n,rhs_n,rhs_k,delta_k_n,delta_k)
  end subroutine inexact_newton

end module mod_newton
