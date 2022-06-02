module mod_newton
!##################################################################################    
!#    Inexact newton method, details may be found on p. 24 in:                    #
!# [1]Franck et al., Energy conservation and numerical stability for the          # 
!#    reduced MHD models of the non-linear JOREK code, 2014, arXiv:1408.2099v3    #
!##################################################################################

  !----------------------- LOADING MODULES, SUBROUTINES, TYPES --------------------
  use construct_matrix_mod, only: construct_matrix
  use mod_clock
  use mpi
#ifdef USE_BICGSTAB
  use mod_bicgstab, only: bicgstab_driver, bicgstab_finalize
#else
  use mod_gmres, only: gmres_driver
#endif
  use mod_gmres, only: gmres_matrix_vector
  use data_structure, only: type_element_list, type_node_list, thread_struct, new_thread_buffers, del_thread_buffers
  use mod_integer_types
  use mod_parameters, only : n_tor, n_var
  use global_distributed_matrix, only: local_index_start, local_index_end
  !----------------------- END OF LOADING -----------------------------------------
 
  implicit none
  integer :: MPI_GLOB, n_cpu
  private
  public :: inexact_newton
  contains  

  subroutine inexact_newton(val, x, b, max_it, tol, comm_glob, comm_n, comm_master,                             &
    ! end of additional arguments of bicgstab_driver 
                            iter_gmres,                                                                         &
    ! end of additional arguments of gmres_driver         
                            element_list, node_list, index_now,                                                 &
    ! end of additional arguments
                            my_id, MPI_COMM_N, my_id_n, MPI_COMM_MASTER, my_id_master, local_elms, n_local_elms,& 
                            index_min, index_max, xpoint2, xcase2, R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint, &
                            Z_xpoint, psi_xpoint,i_tor_min, i_tor_max, n, nz, ndof, n_matrix_block_size, A_mat, &
                            rhs, irn, jcn, ijA_index, ijA_size, irn_jcn, harmonic_matrix                        )
    ! end of arguments of construct_matrix


    ! inexact_newton(arguments of bicgstab_driver,  &
    !                arguments of gmres_driver,     &
    !                additional arguments,          & 
    !                arguments of construct_matrix)                   

    use phys_module, only:  newton_start, newton_gamma,    newton_alpha, newton_eps_a, &
                            newton_eps_r, newton_max_iter, newton_eps_gmres     
    implicit none
    !--------------------- INPUT VARIABLES -----------------------------------------
    !--- definitions for subroutine bicgstab_driver, gmres_driver
    integer,               intent(inout)              :: my_id, my_id_n, my_id_master
    integer,               intent(in)                 :: MPI_COMM_N, MPI_COMM_MASTER
    real(kind=C_DOUBLE)                               :: tol   
    type(type_element_list)                           :: element_list
    type(type_node_list)                              :: node_list
    integer                                           :: index_now, ierr
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
    real(kind=C_DOUBLE), allocatable                  :: rhs_n(:)        ! to copy rhs_glob at step n
    real(kind=C_DOUBLE), allocatable                  :: delta_k_n(:)    ! U(k)-U(n)
    real(kind=C_DOUBLE), allocatable                  :: delta_k(:)      ! initial guess and sol
    real(kind=C_DOUBLE), allocatable                  :: rhs_k(:), rhs_prev(:)
    real(kind=C_DOUBLE), allocatable                  :: deltas_temp(:)
    type(clcktype)                                    :: t0,t1
    real*8                                            :: tsecond
    type(type_element_list)                           :: element_list_temp
    type(type_node_list)                              :: node_list_temp
    integer                                           :: newton_i
    real(kind=C_DOUBLE), allocatable                  :: mat_vec_prod(:)
    real(kind=C_DOUBLE), allocatable                  :: resi(:)
    real(kind=C_DOUBLE), allocatable                  :: resi_n(:)
    real*8                                            :: sqrt_resi_n_2, sqrt_resi_2
    real*8                                            :: sqrt_rhs_n_2, sqrt_rhs_k_2, sqrt_rhs_prev_2
    !--------------------- END OF ROUTINE VARIABLES --------------------------------


    !--------------------- MPI STUFF -----------------------------------------------
    MPI_GLOB = comm_glob
    !call MPI_COMM_RANK(MPI_GLOB,   my_id, ierr)
    !call MPI_COMM_SIZE(MPI_GLOB,   n_cpu, ierr)
    !call MPI_COMM_RANK(MPI_COMM_N, my_id, ierr)
    !--------------------- END OF MPI STUFF ----------------------------------------

    !--------------------- ALLOCATE, ASSIGN VALUES --------------------------------- 
    !call new_thread_buffers() 
    allocate(A_n(1:nz),rhs_n(1:ndof),rhs_k(1:ndof),delta_k_n(1:ndof),  &
             delta_k(1:ndof),deltas_temp(1:ndof),mat_vec_prod(1:ndof), &
             resi(1:ndof), resi_n(1:ndof),rhs_prev(1:ndof))
      
    delta_k_n            = 0.0d0
    delta_k(1:ndof)      = x(1:ndof) * newton_start  ! 1.0E-10
    A_n(1:nz)            = val(1:nz)
    rhs_n(1:ndof)        = b(1:ndof)
    element_list_temp    = element_list
    node_list_temp       = node_list
    deltas_temp(1:ndof)  = x(1:ndof)
    rhs_prev(1:ndof)     = rhs_n(1:ndof)
    !--------------------- END OF ALLOCATION ---------------------------------------

     
    !--------------------- START OF NEWTON ITERATION -------------------------------
    if (my_id.eq.0) write(*,*) '>>>>>>>>>>>>>>>>>>>>>> START NEWTON LOOP <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
    call clck_time(t0)    
    call new_thread_buffers()

    newton_loop: do newton_i = 1, newton_max_iter+1
      if (newton_i.eq.newton_max_iter) then
        if (my_id.eq.0) write(*,*) 'MAXIMUM NUMBER OF NEWTON ITERATIONS. ABORT.'
        stop
      endif 
      if (my_id.eq.0) write(*,'(a,i2.2,a)') '###################### NEWTON ITERATION STEP ', newton_i,' ##########################'
      if (my_id.eq.0) write(*,*)
      !--- iterate commulative delta, i.e. delta_k_n = U(n+1)-U(n) if converged
      delta_k_n(1:ndof) = delta_k_n(1:ndof) + delta_k(1:ndof)
      if (my_id.eq.0) write(*,*) 'ITERATING DELTA_K_N'
      !--- element_list, node_list now local iterates, global values in element_list_temp, node_list_temp
      call update_values(my_id,element_list,node_list,delta_k)
      call update_deltas(my_id,node_list)
      if (my_id.eq.0)  write(*,*) 'UPDATING GLOBAL ELEMENT_LIST, NODE_LIST'
      !--- call construct_matrix, need 'call new_thread_buffers()' for global 'thread_struct', otherwise:
      ! forrtl: severe (408): fort: (7): Attempt to use pointer THREAD_STRUCT when it is not associated with a target
      !call new_thread_buffers() 
      call construct_matrix(my_id, MPI_COMM_N, my_id_n, MPI_COMM_MASTER, my_id_master, local_elms, n_local_elms, index_min, index_max,& 
                            xpoint2, xcase2, R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint, Z_xpoint, psi_xpoint,i_tor_min, i_tor_max,  &
                            n, nz, ndof, n_matrix_block_size, A_mat, rhs, irn, jcn, ijA_index, ijA_size, irn_jcn, harmonic_matrix)
      !call del_thread_buffers()
      if (my_id.eq.0) write(*,*) 'CONSTRUCTED A_GLOB AT U_K'
      !call MPI_BARRIER(MPI_GLOB, ierr)
      call gmres_matrix_vector(ndof,delta_k_n,ndof,mat_vec_prod,my_id)
      call MPI_Bcast(mat_vec_prod,ndof,MPI_DOUBLE_PRECISION,0,MPI_GLOB,ierr)
      rhs_k(1:ndof) = rhs_n(1:ndof) - mat_vec_prod(1:ndof)

      sqrt_rhs_k_2    = DSQRT(DOT_PRODUCT(rhs_k,rhs_k))
      sqrt_rhs_n_2    = DSQRT(DOT_PRODUCT(rhs_n,rhs_n))
      sqrt_rhs_prev_2 = DSQRT(DOT_PRODUCT(rhs_prev,rhs_prev))

      !--- check if initial guess already satisfies convergence criterion
      !if (sqrt_rhs_k_2<1.d-7+1.d-7*sqrt_rhs_n_2) then 
      !  if (my_id.eq.0) write(*,*) 'EXITING NEWTON LOOP AT ', newton_i
      !  if (my_id.eq.0) write(*,'(A8,i3,A13,i2,A18,1E14.4,A13,1E14.4)')           &
      !                          ' t_step=',index_now,' newton_i=', newton_i,      &
      !                          ' R(U_k,U^n) ', sqrt_rhs_k_2, &
      !                          ' R(U^n) '    , sqrt_rhs_n_2
      !  exit newton_loop
      !endif
 

      if (my_id.eq.0) write(*,*) 'COMPUTED MATRIX-VECTOR PRODUCT'
      x(1:ndof) = delta_k(1:ndof)
      b(1:ndof) = rhs_k(1:ndof) 
 
      if(my_id.eq.0) write(*,'(A19,E9.2)') '|rhs_k|/|rhs_prev|', sqrt_rhs_k_2/sqrt_rhs_prev_2     
      if(my_id.eq.0) write(*,'(A8,E9.2)') '|rhs_k|', sqrt_rhs_k_2 
     

      !-- tol is eps_k from [1]
      if(my_id.eq.0) write(*,'(A24,E9.2,E9.2)') 'eps_k, newton_eps_gmres',newton_gamma * (sqrt_rhs_k_2/sqrt_rhs_prev_2)**newton_alpha, newton_eps_gmres
      tol = MIN( newton_gamma * (sqrt_rhs_k_2/sqrt_rhs_prev_2)**newton_alpha, newton_eps_gmres)
      !--- call solvers
#ifdef USE_BICGSTAB
      call bicgstab_driver(irn, jcn, val, x, b, max_it, tol, comm_glob, comm_n, comm_master)
#else
      ! gmres implicitly assumes A_glob * deltas = RHS_glob
      call gmres_driver(my_id,my_id_n,MPI_COMM_N,MPI_COMM_MASTER,iter_gmres)
#endif
 
      delta_k(1:ndof)  = x(1:ndof)


      !--- check convergence
      !--- |Jk.delta_k-rhs_k|/|rhs_k|
      resi   = 0.0d0
      resi_n = 0.0d0
      call gmres_matrix_vector(ndof,delta_k,ndof,resi,my_id)
      call MPI_Bcast(resi,ndof,MPI_DOUBLE_PRECISION,0,MPI_GLOB,ierr)
      resi(1:ndof)   = resi(1:ndof) - rhs_k(1:ndof) 
      !if (my_id.eq.0) write(*,'(A32,i2,2E14.3)') 'newton_i,resi(1),resi(ndof)',newton_i,resi(1),resi(ndof)
      
      A_mat(1:nz) = A_n(1:nz)
      call gmres_matrix_vector(ndof,delta_k+delta_k_n,ndof,resi_n,my_id)
      call MPI_Bcast(resi_n,ndof,MPI_DOUBLE_PRECISION,0,MPI_GLOB,ierr)
      resi_n(1:ndof) = resi_n(1:ndof) - rhs_n(1:ndof)
      !if (my_id.eq.0) write(*,'(A32,i2,2E14.3)') 'newton_i,resi_n(1),resi_n(ndof)',newton_i,resi_n(1),resi_n(ndof)
      
      sqrt_resi_2   = DSQRT(DOT_PRODUCT(resi,resi))
      sqrt_resi_n_2 = DSQRT(DOT_PRODUCT(resi_n,resi_n))
         
      if (my_id.eq.0) write(*,'(A8,i3,A13,i2,A18,2E14.4,A13,2E14.4)')           &
                              ' t_step=',index_now,' newton_i=', newton_i,      &
                              ' r_k, R(U_k,U^n) ', sqrt_resi_2  , sqrt_rhs_k_2, &
                              ' r^n, R(U^n) '    , sqrt_resi_n_2, sqrt_rhs_n_2
      rhs_prev(1:ndof) = rhs_k(1:ndof)
      if (sqrt_rhs_k_2<newton_eps_a+newton_eps_r*sqrt_rhs_n_2) then 
        if (my_id.eq.0) write(*,*) 'EXITING NEWTON LOOP AT ', newton_i
        exit newton_loop
      endif
    end do newton_loop

    call del_thread_buffers() 
    call clck_time_barrier(t1)
    call clck_ldiff(t0,t1,tsecond)
    if (my_id.eq.0) write(*,*) 'ELAPSED TIME IN NEWTON LOOP: ', tsecond
    if (my_id.eq.0) write(*,*) '>>>>>>>>>>>>>>>>>>>>>> END NEWTON LOOP <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
    !--------------------- END OF NEWTON ITERATION ---------------------------------

   
    element_list    = element_list_temp
    node_list       = node_list_temp
    x               = delta_k_n+delta_k!deltas_temp
    if (my_id.eq.0) write(*,*) 'RESET ELEMENT_LIST and NODE_LIST, SAVED DELTA'
    deallocate(A_n,rhs_n,rhs_k,delta_k_n,delta_k,deltas_temp,mat_vec_prod,resi,resi_n,rhs_prev)
  end subroutine inexact_newton

end module mod_newton
