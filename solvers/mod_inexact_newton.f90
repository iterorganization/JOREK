module mod_inexact_newton
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
  use mod_gmres, only: gmres_driver, rinfo
#endif
  use mod_gmres, only: gmres_matrix_vector
  use data_structure, only: type_element_list, type_node_list, thread_struct, new_thread_buffers, del_thread_buffers
  use mod_integer_types
  use mod_parameters, only : n_tor, n_var
  use global_distributed_matrix, only: local_index_start, local_index_end
  !----------------------- END OF LOADING -----------------------------------------
 
  implicit none
  integer :: newton_adapt_time_code = 0 ! report exit status of newton loop for tstep 
  !logical :: newton_adapt_time = .true. ! use adaptive time stepping or not
  private
  public :: inexact_newton, adaptive_tstep, newton_adapt_time_code
  contains  

  subroutine inexact_newton(val, x, b, max_it, tol, MPI_GLOB, comm_n, comm_master,                             &
    ! end of additional arguments of bicgstab_driver 
                            iter_gmres,                                                                        &
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
                            newton_eps_r, newton_max_iter, newton_eps_gmres,           &
                            newton_eps_0, gmres_max_iter, gmres, tstep, newton_adapt_time     
    implicit none
    !--------------------- INPUT VARIABLES -----------------------------------------
    !--- definitions for subroutine bicgstab_driver, gmres_driver
    integer,               intent(inout)              :: my_id, my_id_n, my_id_master
    integer,               intent(in)                 :: MPI_GLOB, MPI_COMM_N, MPI_COMM_MASTER
    real(kind=C_DOUBLE)                               :: tol   
    type(type_element_list)                           :: element_list
    type(type_node_list)                              :: node_list
    integer                                           :: index_now, ierr
    real(kind=C_DOUBLE), pointer, intent(in)          :: val(:)
    real(kind=C_DOUBLE), allocatable                  :: b(:)
    real(kind=C_DOUBLE), allocatable                  :: x(:)
    integer, intent(in)                               :: comm_n, comm_master
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
    real(kind=C_DOUBLE),     allocatable              :: rhs_n(:)        ! to copy rhs_glob at step n
    real(kind=C_DOUBLE),     allocatable              :: delta_k_n(:)    ! U(k)-U(n)
    real(kind=C_DOUBLE),     allocatable              :: delta_k(:)      ! initial guess and sol
    real(kind=C_DOUBLE),     allocatable              :: rhs_k(:)
    real(kind=C_DOUBLE),     allocatable              :: deltas_temp(:)
    type(clcktype)                                    :: t0,t1
    real*8                                            :: tsecond
    type(type_element_list), allocatable              :: element_list_temp
    type(type_node_list),    allocatable              :: node_list_temp
    integer                                           :: newton_i, iter_prev
    real(kind=C_DOUBLE),     allocatable              :: mat_vec_prod(:)
    real(kind=C_DOUBLE),     allocatable              :: resi(:)
    real*8                                            :: normRESk, eps_k
    real*8                                            :: normRHSn, normRHSk, normRHSprev, rhs_ratio
    real*8                                            :: FT, FT_0  ! forcing term [1]
    real(kind=C_DOUBLE),     allocatable              :: tol_array(:), normRHSk_array(:), rhs_ratio_array(:), normRESk_array(:)
    real(kind=C_DOUBLE),     allocatable              :: rinfo1(:), rinfo2(:)
    integer                                           :: i
    integer,                 allocatable              :: iter_array(:)
    integer                                           :: exit_status  ! 0: convergence achieved, 1: max newton iter, 2: max gmres iter
    integer                                           :: n_RHSk_increase = 0
    !--------------------- END OF ROUTINE VARIABLES --------------------------------


    !--------------------- ALLOCATE, ASSIGN VALUES --------------------------------- 
    !call new_thread_buffers() 
    allocate(rhs_n(1:ndof),rhs_k(1:ndof),delta_k_n(1:ndof),delta_k(1:ndof),deltas_temp(1:ndof),mat_vec_prod(1:ndof),resi(1:ndof))
    allocate(element_list_temp, source=element_list)
    allocate(node_list_temp   , source=node_list   )
    allocate(tol_array(1:newton_max_iter),       normRHSk_array(1:newton_max_iter), &
             rhs_ratio_array(1:newton_max_iter), normRESk_array(1:newton_max_iter), iter_array(1:newton_max_iter))
    allocate(rinfo1(1:newton_max_iter), rinfo2(1:newton_max_iter))
    delta_k_n            = 0.0d0
    delta_k(1:ndof)      = x(1:ndof) * newton_start
    rhs_n(1:ndof)        = b(1:ndof)
    element_list_temp    = element_list
    node_list_temp       = node_list
    deltas_temp(1:ndof)  = x(1:ndof)
    normRHSn             = DSQRT(DOT_PRODUCT(rhs_n,rhs_n)) 
    normRHSprev          = normRHSn
    iter_array           = 0.0d0
    !--------------------- END OF ALLOCATION ---------------------------------------

     
    !--------------------- START OF NEWTON ITERATION -------------------------------
    if (my_id.eq.0) write(*,'(A27,A22,A27)') REPEAT('-',27), ' START OF NEWTON LOOP ', REPEAT('-',27)
    call clck_time(t0)    
    call new_thread_buffers()

    newton_loop: do newton_i = 1, newton_max_iter
 
      !--- abort if maximum number of newton iterations is reached
      if (newton_i.eq.newton_max_iter) then
        exit_status = 1
        exit newton_loop
      endif 
      
      !--- reset iter_gmres
      if (gmres) then
        iter_prev = iter_gmres
        iter_gmres = gmres_max_iter
      endif
      
      if (my_id.eq.0) write(*,'(A24,A23,i3,A2,A24)') REPEAT('-',24),' NEWTON ITERATION STEP ', newton_i, '  ',REPEAT('-',24)
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
      call construct_matrix(my_id, MPI_COMM_N, my_id_n, MPI_COMM_MASTER, my_id_master, local_elms, n_local_elms, index_min, index_max,& 
                            xpoint2, xcase2, R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint, Z_xpoint, psi_xpoint,i_tor_min, i_tor_max,  &
                            n, nz, ndof, n_matrix_block_size, A_mat, rhs, irn, jcn, ijA_index, ijA_size, irn_jcn, harmonic_matrix)
      if (my_id.eq.0) write(*,*) 'CONSTRUCTED A_GLOB AT U_K'

      !--- compute J_k.delta_k_n to compute rhs_k
      call gmres_matrix_vector(ndof,delta_k_n,ndof,mat_vec_prod,my_id)
      call MPI_Bcast(mat_vec_prod,ndof,MPI_DOUBLE_PRECISION,0,MPI_GLOB,ierr)
      rhs_k(1:ndof) = rhs_n(1:ndof) - mat_vec_prod(1:ndof)
      if (my_id.eq.0) write(*,*) 'COMPUTED MATRIX-VECTOR PRODUCT'

      !--- compute norms and eps_k
      normRHSk        = DSQRT(DOT_PRODUCT(rhs_k,rhs_k))
      rhs_ratio       = normRHSk/normRHSprev
      normRHSprev     = normRHSk

      !--- abort if rhs_ratio >= 1 for more than 2 iterations
      if ((rhs_ratio.ge.1).or.(isnan(rhs_ratio))) then
        n_RHSk_increase = n_RHSk_increase + 1
        if ((n_RHSk_increase.gt.2).or.(isnan(rhs_ratio))) then
          iter_gmres  = iter_prev
          exit_status = 3
          exit newton_loop
        endif
      else
        n_RHSk_increase = 0
      endif 

      !--- compute eps analogously to [1]
      FT_0    = newton_eps_0
      if (newton_i.eq.1) then
        FT    = FT_0
        eps_k = FT
      else
        FT    = eps_k
        if (newton_gamma*FT**newton_alpha.gt.1.d-1) then
          FT  = MIN(MAX(newton_gamma*rhs_ratio**newton_alpha,newton_gamma*FT**newton_alpha),FT_0)
        else
          FT  = MIN(newton_gamma*rhs_ratio**newton_alpha,FT_0)
        endif
        eps_k = MAX(MIN(0.5*eps_k,FT), newton_eps_gmres)
      endif
 
      !--- set variables for the solvers, tol is eps_k from [1]
      x(1:ndof)   = delta_k(1:ndof)
      b(1:ndof)   = rhs_k(1:ndof) 
      tol         = MAX(eps_k, newton_eps_gmres)  ! set tol>=newton_eps_gmres, just as in standard gmres

      !--- call solvers
#ifdef USE_BICGSTAB
      call bicgstab_driver(irn, jcn, val, x, b, max_it, tol, comm_glob, comm_n, comm_master)
#else
      ! gmres implicitly assumes A_glob * deltas = RHS_glob
      call gmres_driver(my_id,my_id_n,MPI_COMM_N,MPI_COMM_MASTER,iter_gmres)
#endif

      !--- solution of solvers is stored in x
      delta_k(1:ndof)  = x(1:ndof)

      !--- check convergence
      !--- |Jk.delta_k-rhs_k|/|rhs_k|
      resi   = 0.0d0
      call gmres_matrix_vector(ndof,delta_k,ndof,resi,my_id)
      call MPI_Bcast(resi,ndof,MPI_DOUBLE_PRECISION,0,MPI_GLOB,ierr)
      resi(1:ndof) = resi(1:ndof) - rhs_k(1:ndof)       
      normRESk = DSQRT(DOT_PRODUCT(resi,resi))  ! B.E. on unpreconditioned system, i.e. rinfo(2)
         
      !--- save to array for printing
      tol_array(newton_i)       = tol 
      normRHSk_array(newton_i)  = normRHSk  
      rhs_ratio_array(newton_i) = rhs_ratio  
      normRESk_array(newton_i)  = normRESk
      iter_array(newton_i)      = iter_gmres
      rinfo1(newton_i)          = rinfo(1)
      rinfo2(newton_i)          = rinfo(2)
  
      !--- check number of GMRES iterations
      if (iter_gmres.eq.gmres_max_iter) then
        exit_status = 2
        iter_gmres = iter_gmres - 1  ! prevent stop in jorek2main.f90
        exit newton_loop
      endif

      !--- check convergence, residual <= newton_eps_gmres
      if (rinfo(1).le.newton_eps_gmres) then  ! (normRHSk<newton_eps_a+newton_eps_r*normRHSn) then 
        exit_status = 0
        exit newton_loop
      endif
    end do newton_loop

    call del_thread_buffers() 
    call clck_time_barrier(t1)
    call clck_ldiff(t0,t1,tsecond)
    !--------------------- END OF NEWTON ITERATION ---------------------------------

    !--- print information about newton loop
    if(my_id.eq.0) then
      write(*,'(A76)') REPEAT('=',76)
      write(*,'(A3,A23,A50)') REPEAT('-',3), ' Inexact Newton Method ', REPEAT('-',50)
      if (exit_status.eq.0) write(*,'(A4,A30)') '','Convergence has been achieved.'
      if (exit_status.eq.2) write(*,'(A4,A27,i4,A12)') '','No GMRES convergence after ', iter_gmres, ' iterations.'
      if (exit_status.eq.1) write(*,'(A4,A44)') '','Maximum number of Newton iterations reached.'
      if (exit_status.eq.3) write(*,'(A4,A48)') '','Bad behaviour: |R_k|>|R_k-1|, need to recompute.'
      write(*,'(A40,i4,A2,i4,A1)')'Number of Newton (GMRES) iterations:', newton_i,' (',SUM(iter_array),')'
      write(*,'(A40,1f10.2,A1)')  'Elapsed time in inexact Newton loop:', tsecond, 's' 
      write(*,*)
      write(*,'(A3,A18,A55)') REPEAT('-',3), ' Input Parameters ', REPEAT('-',55)
      write(*,'(A11,i9,A17,i4,A7,1f5.2,A7,1E9.2)')'n_step:', index_now, 'newton_iter_max:',newton_max_iter, 'alpha:', newton_alpha,'tol_0:', newton_eps_0
      write(*,'(A11,1f9.3,A17,i4,A7,1f5.2,A7,1E9.2)')'t_step:', tstep, 'gmres_iter_max:', gmres_max_iter, 'gamma:', newton_gamma, 'tol_f:', newton_eps_gmres
      write(*,*)
      write(*,'(A3,A19,A54)') REPEAT('-',3), ' Iteration history ', REPEAT('-',54)
      write(*,'(A5,A5,6A11)') 'i_n', 'i_g', 'gmres_tol', '|R_k|', 'R_k/R_k-1', '|Ax-b|', 'rinfo(1)', 'rinfo(2)'
      if ( ( (exit_status.eq.1).and.(newton_i>1) ) .or. ( (exit_status.eq.3).and.(newton_i>1) )  ) newton_i=newton_i-1
      do i=1,newton_i
        write(*,'(i5,i5,6E11.2)') i,iter_array(i),tol_array(i),normRHSk_array(i),rhs_ratio_array(i),normRESk_array(i),rinfo1(i),rinfo2(i)
      end do
      write(*,'(A76)') REPEAT('=',76)  
    end if
    
    !--- if newton_adapt_time,  define newton_adapt_time_code
    if (newton_adapt_time) then
      if ( exit_status.eq.0                                             ) newton_adapt_time_code = +1  ! fast convergence, increase next tstep
      if ((exit_status.eq.0).and.( SUM(iter_array).ge.gmres_max_iter/2) ) newton_adapt_time_code = -1  ! slow convergence, reduce next tstep
      if ((exit_status.eq.3).or.(exit_status.eq.1).or.(exit_status.eq.2)) newton_adapt_time_code = -2  ! no gmres/newton convergence or bad residue behaviour
    else
    !--- if NOT newton_adapt_time, no convergence for exit_status=1,2,3
      if (exit_status /= 0) then
        if (my_id.eq.0) write(*,*) 'No convergence in Newton loop. Aborting.'
        call MPI_Finalize(ierr)
        stop
      endif
    endif

    !--- reset element_list and node_list, result is stored in deltas, i.e. x and propagated in jorek2_main.f90
    element_list    = element_list_temp
    node_list       = node_list_temp
    !--- reset deltas in case of non-convergence
    if (exit_status.eq.0) then
      x(1:ndof)     = delta_k_n(1:ndof)+delta_k(1:ndof)
    else
      x(1:ndof)     = deltas_temp(1:ndof)
    endif
  
    if (my_id.eq.0) write(*,*) 'RESET ELEMENT_LIST and NODE_LIST, SAVED DELTA'
    if (my_id.eq.0) write(*,'(A28,A20,A28)') REPEAT('-',28), ' END OF NEWTON LOOP ', REPEAT('-',28)
    deallocate(rhs_n,rhs_k,delta_k_n,delta_k,deltas_temp,mat_vec_prod,resi)
    deallocate(element_list_temp, node_list_temp, rinfo1, rinfo2)
    deallocate(tol_array, normRHSk_array, rhs_ratio_array, normRESk_array, iter_array)
  end subroutine inexact_newton


  subroutine adaptive_tstep(newton_adapt_time_code, tstep)
    use phys_module, only: newton_alpha_dt, newton_beta_dt, newton_gamma_dt 
    implicit none
    real*8 ::  t_min, tstep
    integer :: newton_adapt_time_code
    !newton_alpha_dt= 1.05
    !newton_beta_dt = 0.95
    !newton_gamma_dt= 0.75
    t_min          = 0.001   
    if (newton_adapt_time_code.eq.+1) tstep = tstep*newton_alpha_dt              ! fast convergence, increase next tstep
    if (newton_adapt_time_code.eq.-1) tstep = MAX(tstep*newton_beta_dt , t_min)  ! slow convergence, reduce next tstep
    if (newton_adapt_time_code.eq.-2) tstep = MAX(tstep*newton_gamma_dt, t_min)  ! no gmres/newton convergence or bad residue behaviour
  end subroutine adaptive_tstep


end module mod_inexact_newton
