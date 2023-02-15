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
  use mod_sparse_data, only: pastix, mumps, strumpack

  use mod_integer_types

  implicit none

  ! MPI related
  integer                        :: my_id, my_id_n, n_cpu
  integer                        :: MPI_GLOB, MPI_COMM_N, MPI_COMM_MASTER

  ! Global-matrix related
  integer                           :: blocksize, blocksize2
  integer                           :: n_blocks, n_glob, nnz, index_offset, n_local
  real(kind=C_DOUBLE), allocatable  :: b_tmp(:)
  integer, allocatable              :: rcv_c(:), rcv_d(:)

  logical                           :: bicgstab_initialized = .false.


  private
  public :: bicgstab_driver, bicgstab_finalize

  contains

!> solve Ax = b using iterative preconditioned BiCGStab method
  subroutine bicgstab_driver(a_mat, rhs_vec, sol_vec, solver)

    use data_structure,  only: type_SP_MATRIX, type_RHS
    use mod_sparse_data, only: type_SP_SOLVER
    implicit none

    type(type_SP_MATRIX)                    :: a_mat
    type(type_SP_SOLVER)                    :: solver
    type(type_RHS)                          :: sol_vec, rhs_vec

    integer                                 :: comm_glob, comm_n, comm_master, ierr

    real(kind=C_DOUBLE), pointer :: r(:), r_tld(:), s(:), s_hat(:), tmp(:), p(:), p_hat(:), v(:), t(:)
    integer                          :: iter, flag
    real(kind=C_DOUBLE)              :: error, alpha, beta, omega, bnrm2, rho, rho_1, resid, snrm2

    real(kind=C_DOUBLE), external :: dnrm2, ddot ! 2-norm and dot product functions from BLAS
    real :: t0, t1, t2

    MPI_GLOB   = solver%pc%comm
    MPI_COMM_N = solver%pc%MPI_COMM_N
    MPI_COMM_MASTER = solver%pc%MPI_COMM_MASTER

    call MPI_COMM_RANK(MPI_GLOB, my_id, ierr)
    call MPI_COMM_SIZE(MPI_GLOB, n_cpu, ierr)
    call MPI_COMM_RANK(MPI_COMM_N, my_id_n, ierr)

    if (.not.bicgstab_initialized) then
      call bicgstab_init(a_mat)
      bicgstab_initialized = .true.
    endif

    iter = 0
    flag = 0
    t1 = 0; t2 = 0;

    allocate(r(n_glob), r_tld(n_glob), s(n_glob), s_hat(n_glob), &
             tmp(n_glob), p(n_glob), p_hat(n_glob), v(n_glob), t(n_glob))

    bnrm2 = dnrm2(n_glob, rhs_vec%val, 1)

    if (bnrm2 == 0.0) bnrm2 = 1.0

    call matv(a_mat, sol_vec%val, tmp)
    r = rhs_vec%val - tmp

    error = dnrm2(n_glob, r, 1)/bnrm2;

    if (my_id.eq.0) write(*,*) "bicgstab initial relative error:", error

    if (error <= solver%iter_tol) then

      solver%iter_gmres = 1

    else ! go into iterations

      omega  = 1.0
      r_tld = r

      do iter = 1, solver%iter_max
        rho = ddot(n_glob,r_tld,1,r,1) ! direction vector

        if (rho == 0.0) exit

        if (iter > 1) then
          beta  = (rho/rho_1)*(alpha/omega)
          p = r + beta*(p - omega*v)
        else
          p = r
        endif

        t0 = get_time()
        call prec(solver, p, p_hat)
        t1 = t1 + get_time() - t0
        t0 = get_time()
        call matv(a_mat, p_hat, v)
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
        call prec(solver, s, s_hat) ! stabilizer
        t1 = t1 + get_time() - t0
        t0 = get_time()
        call matv(a_mat, s_hat, t)
        t2 = t2 + get_time() - t0

        omega = ddot(n_glob,t,1,s,1)/ddot(n_glob,t,1,t,1)
        sol_vec%val = sol_vec%val + alpha*p_hat + omega*s_hat ! update approximation
        r = s - omega*t
        error = dnrm2(n_glob, r, 1)/bnrm2

        if (my_id.eq.0) write(*,'(A12,1X,I3,5X,A6,1X,E10.3)') "bicgstab it:", iter, "error:", error

        if (error <= solver%iter_tol) exit
        if (omega == 0.0) exit

        rho_1 = rho
      enddo

      solver%iter_gmres = iter

    endif

    if (error <= solver%iter_tol) then ! converged
     flag = 0
     if (my_id.eq.0) write(*,*) "bicgstab completed successfully with n_iter: ", solver%iter_gmres
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
  subroutine matv(a_mat, x, b)
    use data_structure,  only: type_SP_MATRIX

    implicit none

    type(type_SP_MATRIX)                    :: a_mat
    real(kind=8), pointer  :: x(:), b(:)
    integer                           :: i, j, ir, jc
    integer                           :: ierr
    integer                           :: iA_start, ix_start, iy_start
    real(kind=8), allocatable         :: b_tmp_block(:)

    allocate(b_tmp_block(a_mat%block_size))

    !b = 0.d0

    !do i=1,a_mat%nnz
    !  ir = a_mat%irn(i)
    !  jc = a_mat%jcn(i)
    !  b(ir) = b(ir) + a_mat%val(i) * x(jc)
    !enddo
    !call MPI_AllReduce(MPI_IN_PLACE,b,n_glob,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_GLOB,ierr)

    b_tmp = 0.d0

!$omp parallel                                    &
!$omp private(i,iA_start,ix_start,iy_start,b_tmp_block)           &
!$omp reduction(+:b_tmp)
!$omp do schedule(guided)
    do i = 1, n_blocks

      iA_start = (i - 1)*blocksize2
      ix_start = a_mat%jcn(iA_start + 1)
      iy_start = a_mat%irn(iA_start + 1) - index_offset

      call dgemv('T',blocksize,blocksize,1.d0,a_mat%val(iA_start + 1),blocksize,x(ix_start),1,0.d0,b_tmp_block,1)

      b_tmp(iy_start:iy_start + blocksize - 1) = b_tmp(iy_start:iy_start + blocksize - 1) + b_tmp_block(1:blocksize)

    enddo
!$omp end do
!$omp end parallel

    call MPI_Allgatherv(b_tmp,n_local,MPI_DOUBLE_PRECISION,b,rcv_c,rcv_d,MPI_DOUBLE_PRECISION,MPI_GLOB,ierr)
    deallocate(b_tmp_block)

  end subroutine matv

!> apply preconditioner b = M\x
  subroutine prec(solver,x,b)
    use mod_sparse_data, only: type_SP_SOLVER
#ifdef USE_STRUMPACK
    use mod_strumpack, only: strumpack_solve
#endif
#ifdef USE_PASTIX
    use mod_pastix, only: pastix_solve
#endif
#ifdef USE_MUMPS
    use mod_mumps, only: mumps_solve
#endif
    implicit none

    type(type_SP_SOLVER)         :: solver

    real(kind=8), pointer :: x(:), b(:)
    integer :: i
    integer :: ierr
    real :: t0, t1, t2

    !t0 = get_time()
    do i = 1, solver%pc%rhs%n
      solver%pc%rhs%val(i) = x(solver%pc%row_index(i))
    enddo

    !t1 = get_time()
    if (solver%library.eq.strumpack) then
#ifdef USE_STRUMPACK
      call strumpack_solve(solver%spss, solver%pc%rhs)
#endif
    elseif (solver%library.eq.pastix) then
#ifdef USE_PASTIX
      call pastix_solve(solver%ptss, solver%pc%rhs)
#endif
    elseif (solver%library.eq.mumps) then
#ifdef USE_MUMPS
      call mumps_solve(solver%mmss, solver%pc%rhs)
#endif
    endif

    !if (my_id_n.eq.0) write(*,*) my_id, "bicgstab pc solve time", get_time() - t1

    b = 0.d0
    if (solver%pc%my_id_n.eq.0) then
      do i = 1, solver%pc%rhs%n
        b(solver%pc%row_index(i)) = solver%pc%rhs%val(i)*solver%pc%row_factor
      enddo
    endif
    call MPI_BARRIER(MPI_GLOB,ierr)
    call MPI_AllReduce(MPI_IN_PLACE,b,n_glob,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_GLOB,ierr)
    ! now all ranks have the global solution vector

    !if (my_id_n.eq.0) write(*,*) my_id, "bicgstab pc time", get_time() - t0

  end subroutine prec

!> initialize local module variables
  subroutine bicgstab_init(a_mat)
    use data_structure,  only: type_SP_MATRIX

    implicit none

    type(type_SP_MATRIX)  :: a_mat

    integer :: i
    integer(kind=int_all) ::nnz

    ! set module values
    n_glob = a_mat%ng ! rank of global sparse matrix
    nnz = a_mat%nnz ! number of nonzero entries in the local piece of global sparse matrix

    blocksize = a_mat%block_size

    blocksize2 = blocksize*blocksize
    n_blocks = nnz/blocksize2

    index_offset = (a_mat%index_min(my_id + 1) - 1)*blocksize
    n_local   = (a_mat%index_max(my_id + 1) - a_mat%index_min(my_id + 1) + 1)*blocksize

    allocate(b_tmp(n_local))
    allocate(rcv_c(n_cpu),rcv_d(n_cpu))

    do i = 1, n_cpu
      rcv_c(i) = (a_mat%index_max(i) - a_mat%index_min(i) + 1)*blocksize
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

#endif
end module mod_bicgstab
