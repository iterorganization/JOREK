module mod_newton
!##################################################################################    
!#    Inexact newton method, details may be found on p. 24 in:                    #
!#    Franck et al., Energy conservation and numerical stability for the          # 
!#    reduced MHD models of the non-linear JOREK code, 2014, arXiv:1408.2099v3    #
!##################################################################################

  !----------------------- LOADING MODULES, SUBROUTINES, TYPES --------------------
  use construct_matrix_mod, only: construct_matrix
  use mod_clock
#ifdef USE_BICGSTAB
  use mod_bicgstab, only: bicgstab_driver, bicgstab_finalize
#else
  use mod_gmres, only: gmres_driver
#endif
  use mod_gmres, only: gmres_matrix_vector
  use data_structure, only: type_element_list, type_node_list, thread_struct, new_thread_buffers, del_thread_buffers
  use mod_integer_types
  !----------------------- END OF LOADING -----------------------------------------
 
  implicit none
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
    !--------------------- INPUT VARIABLES -----------------------------------------
    !--- definitions for subroutine bicgstab_driver, gmres_driver
    integer,               intent(inout)              :: my_id, my_id_n, my_id_master
    integer,               intent(in)              :: MPI_COMM_N, MPI_COMM_MASTER
    real(kind=C_DOUBLE)                               :: tol   
    type(type_element_list)                           :: element_list
    type(type_node_list)                              :: node_list
    real(kind=C_DOUBLE), pointer, intent(in)          :: val(:)
    real(kind=C_DOUBLE), allocatable                  :: b(:)
    real(kind=C_DOUBLE), allocatable                  :: x(:)
    integer, intent(in)                               :: comm_glob, comm_n, comm_master
    integer, intent(inout)                            :: max_it, iter_gmres
    !--- definitions for subroutine construct_matrix
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
    integer,               intent(in)                 :: i_tor_min
    integer,               intent(in)                 :: i_tor_max
    integer(kind=int_all), intent(inout)              :: n, nz, ndof
    integer(kind=int_all), intent(inout)              :: n_matrix_block_size
    logical,               intent(in)                 :: harmonic_matrix
    real*8,                intent(inout), allocatable :: A_mat(:)
    real*8,                intent(inout), allocatable :: rhs(:)
#ifdef USE_BICGSTAB 
    integer(kind=C_INT),   target, intent(inout), allocatable        :: irn(:), jcn(:)
#else 
    integer(kind=int_all), intent(inout), allocatable :: irn(:), jcn(:)   
#endif
    integer(kind=int_all), intent(inout), allocatable :: ijA_index(:,:), ijA_size(:), irn_jcn(:,:)
    !--------------------- END OF INPUT VARIABLES ----------------------------------


    !--------------------- ROUTINE VARIABLES ---------------------------------------
    real*8, allocatable                               :: A_n(:)          ! to copy A_glob at step n
    real*8, allocatable                               :: rhs_n(:)        ! to copy rhs_glob at step n
    real*8, allocatable                               :: delta_k_n(:)    ! U(k)-U(n)
    real*8, allocatable                               :: delta_k(:)      ! initial guess and sol
    real*8, allocatable                               :: rhs_k(:)
    real*8, allocatable                               :: deltas_temp(:)
    type(clcktype)                                    :: t0,t1
    real*8                                            :: tsecond
    type(type_element_list)                           :: element_list_temp
    type(type_node_list)                              :: node_list_temp
    integer                                           :: newton_iter
    real*8, allocatable                               :: mat_vec_prod(:)
    real*8, allocatable                               :: res_vec(:)
    !--------------------- END OF ROUTINE VARIABLES --------------------------------


    !--------------------- ALLOCATE, ASSIGN VALUES --------------------------------- 
    allocate(A_n(1:nz),rhs_n(1:ndof),rhs_k(1:ndof),delta_k_n(1:ndof),&
             delta_k(1:ndof),deltas_temp(1:ndof),mat_vec_prod(1:ndof),res_vec(1:ndof))
    delta_k_n            = 0.0d0
    delta_k(1:ndof)      = x(1:ndof) * 1.0E-1
    A_n(1:nz)            = val(1:nz)
    rhs_n(1:ndof)        = b(1:ndof)
    element_list_temp    = element_list
    node_list_temp       = node_list
    deltas_temp(1:ndof)  = x(1:ndof)
    !--------------------- END OF ALLOCATION ---------------------------------------


    !--------------------- START OF NEWTON ITERATION -------------------------------
    write(*,*) '>>>>>>>>>>>>>>>>>>>>>> START NEWTON LOOP <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
    call clck_time(t0)    
    newton_loop: do newton_iter = 1, 10
      write(*,'(a,i2.2,a)') '###################### NEWTON ITERATION STEP ', newton_iter,' ##########################'
      write(*,*)
      !--- iterate commulative delta, i.e. delta_k_n = U(n+1)-U(n) if converged
      delta_k_n(1:ndof) = delta_k_n(1:ndof) + delta_k(1:ndof)
      write(*,*) 'ITERATING DELTA_K_N'

      !--- element_list, node_list now local iterates, global values in element_list_temp, node_list_temp
      call update_values(my_id,element_list,node_list,delta_k)
      call update_deltas(my_id,node_list)
      write(*,*) 'UPDATING GLOBAL ELEMENT_LIST, NODE_LIST'

      !--- call construct_matrix, need 'call new_thread_buffers()' for global 'thread_struct', otherwise:
      ! forrtl: severe (408): fort: (7): Attempt to use pointer THREAD_STRUCT when it is not associated with a target
      call new_thread_buffers() 
      call construct_matrix(my_id, MPI_COMM_N, my_id_n, MPI_COMM_MASTER, my_id_master, local_elms, n_local_elms, index_min, index_max,& 
                            xpoint2, xcase2, R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint, Z_xpoint, psi_xpoint,i_tor_min, i_tor_max,  &
                            n, nz, ndof, n_matrix_block_size, A_mat, rhs, irn, jcn, ijA_index, ijA_size, irn_jcn, harmonic_matrix)
      call del_thread_buffers()
      write(*,*) 'CONSTRUCTED A_GLOB AT U_K'

      !call matv(delta_k_n,mat_vec_prod)
      call gmres_matrix_vector(ndof,delta_k_n,ndof,mat_vec_prod,my_id)
      rhs_k(1:ndof) = rhs_n(1:ndof) - mat_vec_prod(1:ndof)
      write(*,*) 'COMPUTED MATRIX-VECTOR PRODUCT'
 
      x(1:ndof) = delta_k(1:ndof)
      b(1:ndof) = rhs_k(1:ndof) 
      !tol = tol*1000
      !--- call solvers
#ifdef USE_BICGSTAB
      call bicgstab_driver(irn, jcn, val, x, b, max_it, tol, comm_glob, comm_n, comm_master)
#else
      ! gmres implicitly assumes A_glob * deltas = RHS_glob
      call gmres_driver(my_id,my_id_n,MPI_COMM_N,MPI_COMM_MASTER,iter_gmres)
#endif

      delta_k(1:ndof) = x(1:ndof)

      !--- check convergence
      !--- |Jk.delta_k-rhs_k|/|rhs_k|
      res_vec = 0.0d0
      call gmres_matrix_vector(ndof,delta_k,ndof,res_vec,my_id)
      res_vec(1:ndof) = res_vec(1:ndof) - rhs_k(1:ndof)
      write(*,*) 'NEWTON RESIDUAL: ', DSQRT(DOT_PRODUCT(res_vec,res_vec)), DSQRT(DOT_PRODUCT(rhs_k,rhs_k))

    enddo newton_loop
 
    call clck_time_barrier(t1)
    call clck_ldiff(t0,t1,tsecond)
    write(*,*) 'ELAPSED TIME IN NEWTON LOOP: ', tsecond
    write(*,*) '>>>>>>>>>>>>>>>>>>>>>> END NEWTON LOOP <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
    !--------------------- END OF NEWTON ITERATION ---------------------------------


    element_list    = element_list_temp
    node_list       = node_list_temp
    x               = delta_k_n+delta_k!deltas_temp
    !x               = delta_k_n+x  

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
    deallocate(A_n,rhs_n,rhs_k,delta_k_n,delta_k,deltas_temp,res_vec)
  end subroutine inexact_newton

end module mod_newton
