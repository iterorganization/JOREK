!> Solve step of the local matrices for each toroidal harmonic (preconditioner for gmres)
subroutine gmres_precondition(x,y,my_id,my_id_n,MPI_COMM_MASTER,MPI_COMM_N)

#ifdef USE_COMPLEX_PRECOND
  use real2complex_mod
#endif
  use tr_module
  use mod_parameters
  use mumps_module
  use pastix_module
  use wsmp_module
  use global_distributed_matrix
  use mpi_mod
  use mod_clock
  use phys_module, only: use_pastix, use_mumps, use_strumpack
  use preconditioner_module, only: my_row_index, my_row_factor
  use mod_integer_types

#ifdef USE_PASTIX6
! -- For PaStiX solver version 6.x
  use pastixf
  use pastix_enums
  use spmf
#elif USE_PASTIX
#include "pastix_fortran.h"
#endif

#ifdef USE_STRUMPACK
  use strumpack_module
#endif

  implicit none

  integer               :: my_id, my_id_n, MPI_COMM_MASTER, MPI_COMM_N
  real, intent (in)     :: x(*)
  real, intent (inout)  :: y(*)

  real, allocatable     :: y_tmp(:)
  integer               :: i, ierr
  real, dimension(:), allocatable :: y_dum
  integer(kind=int_all), parameter   :: Int1=1

  real*8                :: DUMMY_REAL(1:1)
  integer(kind=int_all) :: DUMMY_INT (1:1)

#ifdef USE_PASTIX6
! -- For PaStiX solver version 6.x
  integer(c_int)     :: pastix_info
  type(c_ptr)        :: pastix_rhs_ptr
#endif

  if (associated(mumps_par%rhs)) call tr_deallocatep(mumps_par%rhs,"mumps_par%rhs",CAT_DMATRIX)
  call tr_allocatep(mumps_par%rhs,Int1,mumps_par%n,"mumps_par%rhs",CAT_DMATRIX)

  if (my_id_n.eq.0) then
    call MPI_BCAST(x,ndof_glob,MPI_DOUBLE_PRECISION,0,MPI_COMM_MASTER,ierr)
! construct local RHS from global vector x
    do i = 1, mumps_par%n
      mumps_par%rhs(i) = x(my_row_index(i))
    enddo
  endif

#ifdef USE_MUMPS
  if (use_mumps) then
    mumps_par%JOB = 3 ! Solve
    call DMUMPS(mumps_par)
  endif
#endif

#if defined(USE_PASTIX)||defined(USE_PASTIX6)
  if (use_pastix) then
    if ( (.not. pastix_smp_only) .or. (pastix_smp_only .and. (my_id_n .eq.0)) ) then

      if (.not. associated(mumps_par%rhs)) then
        call tr_allocatep(mumps_par%rhs,Int1,mumps_par%n,"mumps_par%rhs",CAT_DMATRIX)
      endif

      if (.not. pastix_smp_only) call MPI_BCAST(mumps_par%rhs,mumps_par%n,MPI_DOUBLE_PRECISION,0,MPI_COMM_N,ierr)

#ifdef USE_COMPLEX_PRECOND
        !-- converting RHS from real to complex
        call real2complex_rhs(my_id, my_id_n, rhs_cmplx_sol)
        if(my_id_n .gt. 0) then
          if (allocated(rhs_cmplx_sol))  deallocate(rhs_cmplx_sol)
          allocate(rhs_cmplx_sol(1:n_cmplx))
        endif
        call MPI_BCAST(rhs_cmplx_sol,n_cmplx,MPI_DOUBLE_COMPLEX,0,MPI_COMM_N,ierr)
#endif
        ! pastix input parameters working in Pastix5 and Pastix6
        pastix_iparm(IPARM_ITERMAX)               = pastix_iter                ! refinement : max number of iterations
        pastix_iparm(IPARM_FACTORIZATION)         = pastix_facto
        pastix_iparm(IPARM_THREAD_NBR)            = pastix_nthrd               ! number of threads
        pastix_iparm(IPARM_INCOMPLETE)            = pastix_ricar
        pastix_iparm(IPARM_LEVEL_OF_FILL)         = pastix_iluk
        pastix_dparm(DPARM_EPSILON_REFINEMENT)    = pastix_epsilon             ! error level refinement
        pastix_dparm(DPARM_EPSILON_MAGN_CTRL)     = pastix_pivot               ! pivot threshold
#ifdef USE_BLOCK
        pastix_iparm(IPARM_DOF_NBR)               = block_size                 ! block size
#else
        pastix_iparm(IPARM_DOF_NBR)               = 1
#endif


#ifndef USE_PASTIX6
        ! -- For PaStiX solver before version 6.x
        pastix_iparm(IPARM_START_TASK)            = API_TASK_SOLVE
        pastix_iparm(IPARM_END_TASK)              = pastix_endsolve
        pastix_iparm(IPARM_RHS_MAKING)            = pastix_rhs                 ! right hand side (0 : use RHS)
        pastix_iparm(IPARM_SYM)                   = pastix_sym
        pastix_iparm(IPARM_AMALGAMATION_LEVEL)    = pastix_amalg

#else
        ! -- For PaStiX solver version 6.x
        pastix_iparm(IPARM_MTX_TYPE)              = pastix_sym
        pastix_iparm(IPARM_AMALGAMATION_LVLCBLK)  = pastix_amalg
#endif

#ifndef USE_PASTIX6
        ! -- For PaStiX solver before version 6.x
#ifdef USE_BLOCK

#ifndef USE_COMPLEX_PRECOND
        call pastix_fortran(pastix_data,MPI_COMM_N, n_block,                        &
             !mumps_par%jcn,mumps_par%irn,mumps_par%A, &
                   DUMMY_INT, DUMMY_INT, DUMMY_REAL, &
                      pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,Int1,pastix_iparm,pastix_dparm)
#else
        call pastix_fortran(pastix_data,MPI_COMM_N, n_block,                        &
             !mumps_par%jcn,mumps_par%irn,mumps_par%A, &
                   DUMMY_INT, DUMMY_INT, DUMMY_REAL, &
                      pastix_perm_vars,pastix_iperm_vars,rhs_cmplx_sol,Int1,pastix_iparm,pastix_dparm)
#endif

#else

#ifndef USE_COMPLEX_PRECOND
        call pastix_fortran(pastix_data,MPI_COMM_N,mumps_par%n, DUMMY_INT, DUMMY_INT, DUMMY_REAL, &
             pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,Int1,pastix_iparm,pastix_dparm)
#else
        call pastix_fortran(pastix_data,MPI_COMM_N,n_cmplx, DUMMY_INT, DUMMY_INT, DUMMY_REAL, &
             pastix_perm_vars,pastix_iperm_vars,rhs_cmplx_sol,Int1,pastix_iparm,pastix_dparm)
#endif

#endif

#else
         ! -- For PaStiX solver version 6.x
         pastix_rhs_ptr = c_loc(mumps_par%rhs)
#ifdef USE_BLOCK
         call pastix_task_solve(pastix_data,Int1,pastix_rhs_ptr,n_block,pastix_info)
#else
         call pastix_task_solve(pastix_data,Int1,pastix_rhs_ptr,mumps_par%n,pastix_info)
#endif

#endif
    endif
  endif ! use_pastix
#endif /* defined(USE_PASTIX)||defined(USE_PASTIX6) */

#ifdef USE_WSMP
  if (use_wsmp) then
    call MPI_BCAST(mumps_par%rhs,mumps_par%n,MPI_DOUBLE_PRECISION,0,MPI_COMM_N,ierr)
    call PWGSMP__back_substitution(mumps_par%rhs, my_id_n)
  endif
#endif

#ifdef USE_STRUMPACK
  if (use_strumpack) then
    call MPI_BCAST(mumps_par%rhs,mumps_par%n,MPI_DOUBLE_PRECISION,0,MPI_COMM_N,ierr)
    call strumpack_solve(mumps_par%n,mumps_par%rhs,MPI_COMM_N)
  endif
#endif

  if (my_id_n .eq. 0) then

    if (.not.use_strumpack) then
!------------------------------------------ undo column scaling
#ifdef USE_COMPLEX_PRECOND
!-- converting RHS from complex to real
      do i=1,n_cmplx
        if(my_id .eq. 0) then
          mumps_par%rhs(i) = REAL(rhs_cmplx_sol(i))
        else
          mumps_par%rhs(2*i-1) = REAL(rhs_cmplx_sol(i))
          mumps_par%rhs(2*i) = AIMAG(rhs_cmplx_sol(i))
        endif
      enddo
#endif
      do i=1,mumps_par%n
        mumps_par%rhs(i) =  mumps_par%rhs(i) / column_scaling(i)
      enddo
    endif

    allocate(y_dum(ndof_glob))
    y_dum = 0.d0

    ! put local solution into global vector y
    do i = 1, mumps_par%n
      y_dum(my_row_index(i)) = mumps_par%rhs(i)*my_row_factor
    enddo

    call MPI_AllReduce(MPI_IN_PLACE,y_dum,ndof_glob,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_MASTER,ierr)
    if (my_id.eq.0) y(1:ndof_glob) = y_dum(1:ndof_glob)

    deallocate(y_dum)

  endif

  return

end subroutine gmres_precondition
