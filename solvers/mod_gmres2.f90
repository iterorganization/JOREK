module mod_gmres2
!#ifdef USE_GMRES
  use iso_c_binding
  use mpi_mod
  use mod_sparse_data, only: pastix, mumps, strumpack
  use mod_integer_types
  use mod_matv

  private
  public :: gmres2_driver
contains

!> solve a_mat x=b using iterative GMRES method with right (or left) preconditioning
subroutine gmres2_driver(a_mat,b,x,n,solver,converged)
  use mod_sparse_data, only: type_SP_SOLVER
  use data_structure,  only: type_SP_MATRIX

  implicit none
  type(type_SP_MATRIX)                :: a_mat
  real(kind=8), dimension(:), pointer :: x, b
  integer :: n
  type(type_SP_SOLVER)  :: solver
  logical, intent(out)  :: converged

  real(kind=8) :: atol, rtol, gamma, delta, rho, rho0, scale_norm, res_target
  real(kind=8) :: norm_p, norm_gs, norm_gs_new
  integer :: totit, maxit, restart, nrit, it, ldh, k, j
  integer :: n_ortho
  logical :: no_conv, breakdown, hard_breakdown
  logical :: GSC=.false., GSM=.false., GSCI=.true., GSMI=.false.
  !> .true.  : right preconditioning, A M^-1 y = b -> minimises the TRUE residual ||b-Ax||
  !> .false. : left  preconditioning, M^-1 A x = M^-1 b -> minimises ||M^-1 (b-Ax)||
  logical :: right_prec
  !> Relative threshold below which h_{j+1,j} is treated as a (happy) breakdown
  real(kind=8), parameter :: bd_tol = 1.d-14
  !> Floor on the relative residual ||r||/||b||. The convergence target is never set below this.
  real(kind=8), parameter :: res_floor = 1.d-14
  !> V : Arnoldi basis (restart+1 columns), Z : M^-1 V (restart columns, right prec. only),
  !> w_ : scratch vector of length n
  real(kind=8), dimension(:), allocatable, target :: givens_c, givens_s, hess, V, Z, w_, b_, s_

  integer :: my_id, my_id_n, n_cpu, ierr
  integer :: MPI_GLOB, MPI_COMM_N

  external :: dcopy, daxpby, daxpy, dgemv, dscal, dtrsv
  real(kind=8), external :: dnrm2, ddot

  MPI_GLOB   = a_mat%comm
  MPI_COMM_N = solver%pc%MPI_COMM_N

  call MPI_COMM_RANK(MPI_GLOB, my_id, ierr)
  call MPI_COMM_SIZE(MPI_GLOB, n_cpu, ierr)
  call MPI_COMM_RANK(MPI_COMM_N, my_id_n, ierr)

  rtol = solver%iter_tol
  atol = 1.d-36
  maxit = solver%iter_max
  restart = solver%gmres_m
  right_prec = solver%right_prec
  n_ortho = 3 ! Number of iterations for the orthogonalization methos (in case of iterative methods)
  if (restart > maxit) restart = maxit

  allocate(givens_c(restart),givens_s(restart),b_(restart+1),hess((restart+1)*restart),V(n*(restart+1)),w_(n))
  if (right_prec) allocate(Z(n*restart))
  if (GSCI .or. GSMI) allocate(s_(restart))

  givens_c(1:restart) = 0.d0
  givens_s(1:restart) = 0.d0

  ldh = restart+1

  ! --- Reference norm for the relative stopping criterion ---
  if (right_prec) then
    scale_norm = dnrm2(n, b(1:n), 1)
  else
    call prec(solver, b(1:n), w_(1:n), n, MPI_GLOB, MPI_COMM_N)
    scale_norm = dnrm2(n, w_(1:n), 1)
  endif
  if (scale_norm .le. 0.d0) scale_norm = 1.d0

  no_conv        = .true.
  breakdown      = .false.
  hard_breakdown = .false.
  converged      = .false.
  totit = 0

  rho0       = 0.d0
  res_target = 0.d0

  do while (no_conv)
    ! --- v_1 = A * x ---
    call bcsr_matv(a_mat, x, V(1:n))

    ! --- v_1 = b - A * x (True residual) ---
    call daxpby(n, 1.d0, b(1:n), 1, -1.d0, V(1:n), 1)

    ! --- Left preconditioning: v_1 = M^-1 (b - A * x) ---
    if (.not. right_prec) then
      call dcopy(n, V(1), 1, w_(1), 1)
      call prec(solver, w_(1:n), V(1:n), n, MPI_GLOB, MPI_COMM_N)
    endif

    ! --- rho = ||v_1||_2 ---
    rho = dnrm2(n, V(1:n), 1)

    ! --- Convergence target ---
    if (totit .eq. 0) then
      rho0 = rho
      res_target = max(rtol*rho0, res_floor*scale_norm)
    endif

    if ((rho .le. res_target) .or. (rho .lt. atol)) then
      no_conv   = .false.
      converged = .true.
      exit
    endif
    !--- v_1 = v_1 / rho
    call dscal(n, 1.d0/rho, V(1:n), 1)
    b_(1) = rho
    b_(2:restart+1) = 0.d0
    nrit = restart-1
    if (my_id.eq.0) then
      write(*, "(A, X, I0, X, A, X, ES14.6, X, A, X, ES14.6, X, A, X, ES14.6)") &
        "[GMRES] iteration", totit, "res =", rho, "res/res0 =", rho/rho0, "res/rhs =", rho/scale_norm
      write(*, "(A)") "[GMRES] --- Restart ---"
    endif

    do it = 1, restart
      totit = totit +1
      if (right_prec) then
        ! --- z_j = M^-1 * v_j (kept, so that the solution update needs no extra solve) ---
        call prec(solver, V((it-1)*n+1:it*n), Z((it-1)*n+1:it*n), n, MPI_GLOB, MPI_COMM_N)

        ! --- v_j+1 = A * z_j = A * M^-1 * v_j ---
        call bcsr_matv(a_mat, Z((it-1)*n+1:it*n), V(it*n+1:(it+1)*n))
      else
        ! --- w_ = A * v_j ---
        call bcsr_matv(a_mat, V((it-1)*n+1:it*n), w_(1:n))

        ! --- v_j+1 = M^-1 * w_ = M^-1 * A * v_j ---
        call prec(solver, w_(1:n), V(it*n+1:(it+1)*n), n, MPI_GLOB, MPI_COMM_N)
      endif

      ! --- Norm of the un-orthogonalized vector (reference for the breakdown test). It must
      !     NOT be overwritten by the iterative Gram-Schmidt loops below, which track their
      !     own norm: after heavy cancellation their norm collapses with the vector and the
      !     breakdown test would then never trigger. ---
      norm_p = dnrm2(n, V(it*n+1), 1)

      ! --- Orthogonalization ---
      if (GSC) then ! Gram-Schmidt Classical
        call dgemv('C', n, it, 1.d0, V(1), n, V(it*n+1), 1, 0.d0, hess((it-1)*ldh+1), 1)
        call dgemv('N', n, it, -1.d0, V(1), n, hess((it-1)*ldh+1), 1, 1.d0, V(it*n+1), 1)
      elseif (GSM) then ! Gram-Schmidt Modified
        do k=1,it
          hess(k+(it-1)*ldh) = ddot(n, V((k-1)*n+1), 1, V(it*n+1), 1)
          call daxpy(n, -hess(k+(it-1)*ldh), V((k-1)*n+1), 1, V(it*n+1), 1)
        enddo
      elseif (GSCI) then ! Gram-Schmidt Classical Iterative
        hess((it-1)*ldh+1:(it-1)*ldh+1+it) = 0.d0
        norm_gs = norm_p
        do j=1,n_ortho
          call dgemv('C', n, it, 1.d0, V(1), n, V(it*n+1), 1, 0.d0, s_(1), 1)
          call dgemv('N', n, it, -1.d0, V(1), n, s_(1), 1, 1.d0, V(it*n+1), 1)
          call daxpy(it, 1.d0, s_(1), 1, hess((it-1)*ldh+1), 1)
          norm_gs_new = dnrm2(n, V(it*n+1), 1)
          if (2.d0 * norm_gs_new .gt. norm_gs) exit ! Stopping criterion for iterative GS methods
          norm_gs = norm_gs_new
        enddo
      elseif (GSMI) then ! Gram-Schmidt Modified Iterative
        hess((it-1)*ldh+1:(it-1)*ldh+1+it) = 0.d0
        norm_gs = norm_p
        do j=1,n_ortho
          do k=1,it
            s_(k) = ddot(n, V((k-1)*n+1), 1, V(it*n+1), 1)
            call daxpy(n, -s_(k), V((k-1)*n+1), 1, V(it*n+1), 1)
          enddo
          call daxpy(it, 1.d0, s_(1), 1, hess((it-1)*ldh+1), 1)
          norm_gs_new = dnrm2(n, V(it*n+1), 1)
          if (2.d0 * norm_gs_new .gt. norm_gs) exit ! Stopping criterion for iterative GS methods
          norm_gs = norm_gs_new
        enddo
      endif
      ! --- h_j+1,j = ||v_j+1||_2 ---
      hess(it+(it-1)*ldh+1) = dnrm2(n, V(it*n+1), 1)
      ! --- Breakdown test ---
      if (hess(it+(it-1)*ldh+1) .le. bd_tol*max(norm_p, atol)) then
        hess(it+(it-1)*ldh+1) = 0.d0
        V(it*n+1:(it+1)*n) = 0.d0
        breakdown = .true.
      else
        ! --- v_j+1 = v_j+1 / h_j+1,j ---
        call dscal(n, 1.d0/hess(it+(it-1)*ldh+1), V(it*n+1), 1)
      endif
      ! --- Givens Rotation ---
      do k = 1, it-1
        gamma = givens_c(k)*hess(k+(it-1)*ldh) + givens_s(k)*hess(k+(it-1)*ldh+1)
        hess(k+(it-1)*ldh+1) = -givens_s(k)*hess(k+(it-1)*ldh) + givens_c(k)*hess(k+(it-1)*ldh+1)
        hess(k+(it-1)*ldh) = gamma
      enddo
      delta = sqrt(abs(hess(it+(it-1)*ldh))**2 + hess(it+(it-1)*ldh+1)**2);
      if (delta .le. 0.d0) then
        ! Hard breakdown: the whole column vanished, no further progress is possible.
        givens_c(it) = 1.d0
        givens_s(it) = 0.d0
        breakdown      = .true.
        hard_breakdown = .true.
      else
        givens_c(it) = hess(it+(it-1)*ldh) / delta
        givens_s(it) = hess(it+(it-1)*ldh+1) / delta
      endif
      hess(it+(it-1)*ldh) = givens_c(it)*hess(it+(it-1)*ldh) + givens_s(it)*hess(it+(it-1)*ldh+1)
      b_(it+1) = -givens_s(it)*b_(it)
      b_(it) = givens_c(it)*b_(it)
      rho = abs(b_(it+1))
      if (my_id.eq.0) write(*, "(A, X, I0, X, A, X, ES14.6, X, A, X, ES14.6, X, A, X, ES14.6)") &
        "[GMRES] iteration", totit, "res =", rho, "res/res0 =", rho/rho0, "res/rhs =", rho/scale_norm
      if ((rho .le. res_target).or.(rho .lt. atol).or.(totit .ge. maxit).or.breakdown) then
        no_conv = .false.
        ! A happy breakdown leaves rho = 0 -> convergence.
        ! A hard breakdown zeroes rho as well -> no convergence.
        converged = ((rho .le. res_target) .or. (rho .lt. atol)) .and. (.not. hard_breakdown)
        nrit = it-1
        ! The last column of the triangular factor is zero after a hard breakdown: drop it
        ! and keep the iterate of the previous step instead of dividing by zero.
        if (hard_breakdown) nrit = it-2
        exit
      endif

    enddo
    ! --- Solve upper triangular system b_ = H \ b_ ---
    call dtrsv('U', 'N', 'N', nrit+1, hess, ldh, b_, 1)
    if (right_prec) then
      ! --- Update solution x = x + Z * b_ = x + M^-1 * V * b_ ---
      call dgemv('N', n, nrit+1, 1.d0, Z(1), n, b_(1), 1, 1.d0, x(1:n), 1)
    else
      ! --- Update solution x = x + V * b_ ---
      call dgemv('N', n, nrit+1, 1.d0, V(1), n, b_(1), 1, 1.d0, x(1:n), 1)
    endif

  enddo

  solver%iter_gmres = totit

  deallocate(givens_c,givens_s,b_,hess,V,w_)
  if (right_prec) deallocate(Z)
  if (GSCI .or. GSMI) deallocate(s_)

end subroutine gmres2_driver


!> apply preconditioner b = M^-1 x
  subroutine prec(solver, x, b, n_glob, MPI_GLOB, MPI_COMM_N)
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

    !real(kind=8), pointer :: x(:), b(:)
    real(kind=8), dimension(:), intent(inout) :: x, b
    integer :: n_glob
    integer :: i
    integer :: ierr, MPI_GLOB, MPI_COMM_N, my_id_n, my_id
    real :: t0, t1, t2
    real(kind=8), external :: dnrm2

    call MPI_COMM_RANK(MPI_COMM_N, my_id_n, ierr)
    call MPI_COMM_RANK(MPI_GLOB, my_id, ierr)

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

    !if (my_id_n.eq.0) write(*,*) my_id, "gmres pc solve time", get_time() - t1

    b = 0.d0
    if (my_id_n.eq.0) then
      do i = 1, solver%pc%rhs%n
        b(solver%pc%row_index(i)) = solver%pc%rhs%val(i)*solver%pc%row_factor
      enddo
    endif
    !call MPI_BARRIER(MPI_GLOB,ierr)
    call MPI_AllReduce(MPI_IN_PLACE,b,n_glob,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_GLOB,ierr)
    ! now all ranks have the global solution vector

    !if (my_id_n.eq.0) write(*,*) my_id, "gmres pc time", get_time() - t0

  end subroutine prec

end module mod_gmres2
