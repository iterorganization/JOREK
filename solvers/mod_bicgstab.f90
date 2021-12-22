module mod_bicgstab
!     Details of this algorithm are described in "Templates for the
!     Solution of Linear Systems: Building Blocks for Iterative
!     Methods", Barrett, Berry, Chan, Demmel, Donato, Dongarra,
!     Eijkhout, Pozo, Romine, and van der Vorst, SIAM Publications,
!     1993. (ftp netlib2.cs.utk.edu; cd linalg; get templates.ps).
!     http://www.netlib.org/templates/matlab/bicgstab.m
#ifdef USE_BICGSTAB
  use iso_c_binding
  use mpi
  use mumps_module, only: mumps_par
  use phys_module, only: use_pastix, use_mumps, use_strumpack

  implicit none

  type SPARSE_MATRIX_T
    integer(kind=C_INT), pointer :: irn(:), jcn(:)
    real(kind=C_DOUBLE), pointer :: val(:)
    integer                      :: indexing
    integer(kind=C_INT)          :: n
    integer(kind=C_INT)          :: nnz
  end type SPARSE_MATRIX_T

  type(SPARSE_MATRIX_T)          :: cooA

  ! MPI related
  integer                        :: my_id, my_id_n, n_cpu
  integer                        :: MPI_GLOB, MPI_COMM_N, MPI_COMM_MASTER

  ! Global-matrix related
  integer                           :: blocksize, blocksize2
  integer                           :: n_blocks, n_glob, nnz, index_offset, n_local
  real(kind=C_DOUBLE), allocatable  :: b_tmp(:)
  integer, allocatable              :: rcv_c(:), rcv_d(:)

  logical                        :: bicgstab_initialized = .false.


  private
  public :: bicgstab_driver, bicgstab_finalize

  contains

!> solve Ax = b using iterative preconditioned BiCGStab method
!! x=sol contains the initial guess
!! max_it - maximum number of iterations
!! tol - iteration tolerance
  subroutine bicgstab_driver(irn, jcn, val, x, b, max_it, tol, comm_glob, comm_n, comm_master)
    implicit none

    integer(kind=C_INT), pointer, intent(in)    :: irn(:), jcn(:)
    real(kind=C_DOUBLE), pointer, intent(in)    :: val(:)
    real(kind=C_DOUBLE), allocatable            :: b(:)
    real(kind=C_DOUBLE), allocatable            :: x(:)
    integer, intent(in)                         :: comm_glob, comm_n, comm_master
    integer, intent(inout)                      :: max_it

    real(kind=C_DOUBLE), allocatable, target :: r(:), r_tld(:), s(:), s_hat(:), tmp(:), p(:), p_hat(:), v(:), t(:)
    integer                          :: iter, flag
    integer                          :: ierr
    real(kind=C_DOUBLE)              :: tol, error, alpha, beta, omega, bnrm2, rho, rho_1, resid, snrm2

    real(kind=C_DOUBLE), external :: dnrm2, ddot ! 2-norm and dot product functions from BLAS
    real :: t0, t1, t2

    MPI_GLOB = comm_glob
    MPI_COMM_N = comm_n
    MPI_COMM_MASTER = comm_master
    call MPI_COMM_RANK(MPI_GLOB, my_id, ierr)
    call MPI_COMM_SIZE(MPI_GLOB, n_cpu, ierr)
    call MPI_COMM_RANK(MPI_COMM_N, my_id_n, ierr)

    if (.not.bicgstab_initialized) then
      call bicgstab_init()
      bicgstab_initialized = .true.
    endif

    cooA%irn => irn
    cooA%jcn => jcn
    cooA%val => val
    cooA%indexing = 1
    cooA%n = n_glob
    cooA%nnz = nnz

    call MPI_Bcast(b,n_glob,MPI_DOUBLE_PRECISION,0,MPI_GLOB,ierr)
    call MPI_Bcast(x,n_glob,MPI_DOUBLE_PRECISION,0,MPI_GLOB,ierr)

    iter = 0
    flag = 0

    allocate(r(n_glob), r_tld(n_glob), s(n_glob), s_hat(n_glob), &
             tmp(n_glob), p(n_glob), p_hat(n_glob), v(n_glob), t(n_glob))

    bnrm2 = dnrm2(n_glob, b, 1)

    if (bnrm2 == 0.0) bnrm2 = 1.0

    call matv(x,tmp)
    r = b - tmp

    error = dnrm2(n_glob, r, 1)/bnrm2;
    if (error <= tol) then
      !if (my_id.eq.0) write(*,*) "bicgstab exiting, initial relative error:", error
      max_it = 0
    endif

    omega  = 1.0
    r_tld = r

    t1 = 0; t2 = 0;

    if (my_id.eq.0) write(*,*) "bicgstab initial relative error:", error

    do iter = 1, max_it
      rho = ddot(n_glob,r_tld,1,r,1) ! direction vector

      if (rho == 0.0) exit

      if (iter > 1) then
        beta  = (rho/rho_1)*(alpha/omega)
        p = r + beta*(p - omega*v)
      else
        p = r
      endif

      t0 = get_time()
      call prec(p,p_hat)
      t1 = t1 + get_time() - t0
      t0 = get_time()
      call matv(p_hat,v)
      t2 = t2 + get_time() - t0

      alpha = rho/ddot(n_glob,r_tld,1,v,1)
      s = r - alpha*v
      !snrm2 = dnrm2(n_glob, s, 1)

      !if (snrm2 < tol) then
      !  x = x + alpha*p_hat
      !  resid = snrm2/bnrm2
      !  exit
      !endif

      t0 = get_time()
      call prec(s,s_hat) ! stabilizer
      t1 = t1 + get_time() - t0
      t0 = get_time()
      call matv(s_hat,t)
      t2 = t2 + get_time() - t0

      omega = ddot(n_glob,t,1,s,1)/ddot(n_glob,t,1,t,1)
      x = x + alpha*p_hat + omega*s_hat ! update approximation
      r = s - omega*t
      error = dnrm2(n_glob, r, 1)/bnrm2

      if (error <= tol) exit
      if (omega == 0.0) exit

      rho_1 = rho
    enddo

    max_it = iter! - 1 ! actual number of iterations

    if (error <= tol) then ! converged
     flag = 0
     if (my_id.eq.0) write(*,*) "bicgstab completed successfully with n_iter: ", max_it
    elseif (omega == 0.0) then ! breakdown
     flag = -2
     if (my_id.eq.0) write(*,*) "bicgstab fails, flag: ", flag
    elseif (rho == 0.0) then
     flag = -1
     write(*,*) "bicgstab fails, flag: ", flag
    else ! no convergence
     flag = 1
     if (my_id.eq.0) write(*,*) "bicgstab failed to converge, relative error: ", error
    endif

    if (my_id.eq.0) write(*,*) "bicgstab times (prec, matv):", t1, t2
    if (my_id.eq.0) write(*,*) "bicgstab final relative error:", error


    deallocate(r, r_tld, s, s_hat, tmp, p, p_hat, v, t)

  end subroutine bicgstab_driver

!> get matrix-vector product b=Ax
  subroutine matv(x,b)
    implicit none

    real(kind=C_DOUBLE), allocatable  :: x(:), b(:)
    integer                           :: i, j, ir, jc
    integer                           :: ierr
    integer                           :: iA_start, ix_start, iy_start
    real(kind=C_DOUBLE)               :: b_tmp_block(blocksize)

    !b = 0.d0

    !do i=1,cooA%nnz
    !  ir = cooA%irn(i)
    !  jc = cooA%jcn(i)
    !  b(ir) = b(ir) + cooA%val(i) * x(jc)
    !enddo
    !call MPI_AllReduce(MPI_IN_PLACE,b,n_glob,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_GLOB,ierr)

    b_tmp = 0.d0

!$omp parallel default(none)                                      &
!$omp shared(cooA,x,n_blocks,blocksize,blocksize2,index_offset)   &
!$omp private(i,iA_start,ix_start,iy_start,b_tmp_block)           &
!$omp reduction(+:b_tmp)
!$omp do schedule(guided)
    do i = 1, n_blocks

      iA_start = (i - 1)*blocksize2
      ix_start = cooA%jcn(iA_start + 1)
      iy_start = cooA%irn(iA_start + 1) - index_offset

      call dgemv('T',blocksize,blocksize,1.d0,cooA%val(iA_start + 1),blocksize,x(ix_start),1,0.d0,b_tmp_block,1)

      b_tmp(iy_start:iy_start + blocksize - 1) = b_tmp(iy_start:iy_start + blocksize - 1) + b_tmp_block(1:blocksize)

    enddo
!$omp end do
!$omp end parallel

    call MPI_Allgatherv(b_tmp,n_local,MPI_DOUBLE_PRECISION,b,rcv_c,rcv_d,MPI_DOUBLE_PRECISION,MPI_GLOB,ierr)

  end subroutine matv

!> apply preconditioner b = M\x
  subroutine prec(x,b)
    use preconditioner_module, only: my_row_index, my_row_factor
#ifdef USE_STRUMPACK
    use strumpack_module, only: strumpack_solve
#endif
    implicit none

    real(kind=C_DOUBLE), allocatable, target :: x(:), b(:)
    integer :: i
    integer :: ierr
    real :: t0, t1, t2

    !t0 = get_time()

    if (.not.associated(mumps_par%rhs)) allocate(mumps_par%rhs(mumps_par%n))

    if (my_id_n.eq.0) then
      do i = 1, mumps_par%n
        mumps_par%rhs(i) = x(my_row_index(i))
      enddo
    endif
    call MPI_Bcast(mumps_par%rhs,mumps_par%n,MPI_DOUBLE_PRECISION,0,MPI_COMM_N,ierr)

    !t1 = get_time()
    if (use_strumpack) then
#ifdef USE_STRUMPACK
      call strumpack_solve(mumps_par%n,mumps_par%rhs,MPI_COMM_N)
#endif
    elseif (use_pastix) then
#ifdef USE_PASTIX
      call pastix_solve(mumps_par%n,mumps_par%rhs,MPI_COMM_N)
#endif
    endif

    !if (my_id_n.eq.0) write(*,*) my_id, "bicgstab pc solve time", get_time() - t1

    b = 0.d0
    if (my_id_n.eq.0) then
      do i = 1, mumps_par%n
        b(my_row_index(i)) = mumps_par%rhs(i)*my_row_factor
      enddo
    endif
    call MPI_BARRIER(MPI_GLOB,ierr)
    call MPI_AllReduce(MPI_IN_PLACE,b,n_glob,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_GLOB,ierr)
    ! now all ranks have the global solution vector

    !if (my_id_n.eq.0) write(*,*) my_id, "bicgstab pc time", get_time() - t0

  end subroutine prec

!> initialize local module variables
  subroutine bicgstab_init()
    use mod_parameters, only : n_tor, n_var
    use global_distributed_matrix, only: ndof_glob, nz_glob, local_index_start, local_index_end
    implicit none

    integer :: i

    ! set module values
    n_glob = ndof_glob ! rank of global sparse matrix
    nnz = nz_glob ! number of nonzero entries in the local piece of global sparse matrix

    blocksize = n_tor*n_var

    blocksize2 = blocksize*blocksize
    n_blocks = nz_glob/blocksize2

    index_offset = (local_index_start(my_id + 1) - 1)*blocksize
    n_local   = (local_index_end(my_id + 1) - local_index_start(my_id + 1) + 1)*blocksize

    allocate(b_tmp(n_local))
    allocate(rcv_c(n_cpu),rcv_d(n_cpu))

    do i = 1, n_cpu
      rcv_c(i) = (local_index_end(i) - local_index_start(i) + 1)*blocksize
    enddo

    rcv_d(1) = 0
    do i = 2, n_cpu
      rcv_d(i) = rcv_d(i-1) + rcv_c(i-1)
    enddo

  end subroutine bicgstab_init

!> clean-up local module variables
  subroutine bicgstab_finalize()
    implicit none

    deallocate(rcv_c,rcv_d)
    deallocate(b_tmp)
  end subroutine bicgstab_finalize

  real function get_time()
    implicit none
    integer :: cc, cr

    call system_clock(count=cc, count_rate=cr)
    get_time =  real(cc)/cr
  end function get_time

#ifdef USE_PASTIX
!> call PaStiX solver
  subroutine pastix_solve(n, rhs, comm)
#include "pastix_fortran.h"
    use pastix_module
    use global_distributed_matrix, only: column_scaling
    implicit none

    real(kind=C_DOUBLE), dimension(:), pointer :: rhs
    integer(kind=C_INT_ALL), intent(in) :: n
    integer, intent(in) :: comm

    integer :: ierr
    integer :: i, n_blocks_loc, blocksize_loc
    real*8  :: DUMMY_REAL(1:1)
    integer :: DUMMY_INT (1:1)

    pastix_iparm(IPARM_START_TASK) = API_TASK_SOLVE
    pastix_iparm(IPARM_END_TASK)   = pastix_endsolve
    pastix_iparm(IPARM_RHS_MAKING) = pastix_rhs

    blocksize_loc = pastix_iparm(IPARM_DOF_NBR)
    n_blocks_loc = n/blocksize_loc

    call pastix_fortran(pastix_data,comm, n_blocks_loc, DUMMY_INT, DUMMY_INT, DUMMY_REAL, &
                        pastix_perm_vars, pastix_iperm_vars, rhs, 1, pastix_iparm, pastix_dparm)

    if (my_id_n.eq.0) then
      do i = 1, n
        rhs(i) = rhs(i)/column_scaling(i)
      enddo
    endif

    call MPI_BARRIER(comm,ierr)

  end subroutine pastix_solve
#endif

#endif
end module mod_bicgstab
