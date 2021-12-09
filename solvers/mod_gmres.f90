module mod_gmres

  use tr_module
  use mod_integer_types
  use mod_clock
  use mumps_module, only: mumps_par
  use mod_parameters, only : n_tor, n_var
  use global_distributed_matrix, only: A_glob, IRN_glob, JCN_glob, RHS_glob, deltas, ndof_glob, nz_glob, &
                                       local_index_start, local_index_end, column_scaling

  private
  public :: gmres_driver


  contains

!> Driver for the reverse communication GMRES routine from dPackgmres (CERFACS)
  subroutine gmres_driver(my_id,my_id_n,MPI_COMM_N,MPI_COMM_MASTER,iter_gmres)

    use phys_module,  only : gmres_4, gmres_tol, gmres_m, index_now
    use mod_dpackgmres, only: init_dgmres, drive_dgmres

    implicit none
#include "r3_info.h"

    integer               :: j, my_id, my_id_n, my_id_master, MPI_COMM_N, MPI_COMM_MASTER
    integer               :: revcom, iter_gmres
    integer               :: icntl(8), info(3)
    integer(kind=int_all) :: i, m, colx, coly, colz, nbscal, lwork, n_dof, irc(5)

    integer(kind=int_all), parameter   :: Int1=1

    integer :: matvec, precondLeft, precondRight, dotProd, ierr
    real*8  :: cntl(5), rinfo(2), sum, err, Bnorm, Xnorm
    real*8, allocatable :: work(:), work_ndof(:), work_ndof2(:)
    type(clcktype)           :: t0, t1
    real*8                   :: tsecond
    REAL*8, ALLOCATABLE      :: rhs_tmp(:)
    real*8 ::ZERO, ONE
    CHARACTER(LEN=128) :: fname

    parameter (ZERO = 0.0d0, ONE = 1.0d0)
    parameter (matvec=1, precondLeft=2, precondRight=3, dotProd=4)

    !write(*,*) ' GMRES DRIVER : ',my_id,my_id_n
    call r3_info_begin (r3_info_index_0, 'gmres_driver')  ! timing
    call clck_time(t0)
    call init_dgmres(icntl,cntl)

    icntl(3) = 6            ! output unit
    icntl(7) = iter_gmres   ! Maximum number of iterations
    icntl(4) = 1            ! preconditioner (1) = left preconditioner
    icntl(5) = 3            ! orthogonalization scheme
    icntl(6) = 1            ! initial guess  (1) = user supplied guess
    icntl(8) = 1            ! residual calculation strategy at restart

    cntl(1) = gmres_tol         ! stopping tolerance
    cntl(2) = 1.d0
    cntl(3) = 1.d0
    cntl(4) = gmres_4          ! 1.d0
    cntl(5) = 1.d0

    m = gmres_m

    n_dof = ndof_glob

    lwork = m*m + m*(n_dof+5) + 6*n_dof + m + 1

    call tr_allocate(work,Int1,lwork,"work",CAT_GMRES)
    call tr_allocate(work_ndof,Int1,n_dof,"work_ndof",CAT_GMRES)
    call tr_allocate(work_ndof2,Int1,n_dof,"work_ndof2",CAT_GMRES)

    work(1:n_dof)         = deltas(1:n_dof)                     ! the initial guess
    work(n_dof+1:2*n_dof) = RHS_glob(1:n_dof)                   ! the right hand side

    work_ndof(Int1:n_dof)  = work(Int1:n_dof)
    work_ndof2(Int1:n_dof) = work(2*n_dof+Int1:3*n_dof)

    call gmres_matrix_vector(n_dof,work_ndof,n_dof,work_ndof2,my_id)

    work(Int1:n_dof)           = work_ndof(Int1:n_dof)
    work(2*n_dof+Int1:3*n_dof) = work_ndof2(Int1:n_dof)

    if (my_id .eq. 0) then
      sum = 0.d0
      err = -1.d20
      Bnorm = 0.d0
      Xnorm = 0.d0
      do i=1,n_dof
        sum = sum + (work(2*n_dof+i)-work(n_dof+i))**2
        err = max(err,abs(work(2*n_dof+i)-work(n_dof+i)))
        Bnorm = Bnorm + RHS_glob(i)**2
        Xnorm = Xnorm + deltas(i)**2
      enddo

      write(*,'(A,4e16.8)') ' residu test before : ',sqrt(sum),err,sqrt(Bnorm),sqrt(Xnorm)
    endif
    if (my_id .eq. 0) then
       write(fname,'(A,I6.6,A1,I6.6)')"product_before",index_now
       call tr_vdump(fname, work(2*n_dof+1:3*n_dof), n_dof)
       write(fname,'(A,I6.6,A1,I6.6)')"X_before",index_now
       call tr_vdump(fname, RHS_glob, n_dof)
    end if

    !*****************************************
    !** Reverse communication implementation
    !*****************************************

    10     call MPI_barrier(MPI_COMM_WORLD,ierr)

           if (my_id .eq. 0) then
             call drive_dgmres(n_dof,n_dof,m,lwork,work,irc,icntl,cntl,info,rinfo)
           endif

           call MPI_BCAST(irc,5,MPI_INTEGER_ALL,0,MPI_COMM_WORLD,ierr)
           revcom = irc(1)
           colx   = irc(2)
           coly   = irc(3)
           colz   = irc(4)
           nbscal = irc(5)

           if (revcom.eq.matvec) then                  ! perform the matrix vector product
                                                       ! work(colz) <-- A * work(colx)

             work_ndof(Int1:n_dof) = work(colx:colx+n_dof-Int1)
             work_ndof2(Int1:n_dof) = work(colz:colz+n_dof-Int1)
             call gmres_matrix_vector(n_dof,work_ndof,n_dof,work_ndof2,my_id)
             work(colx:colx+n_dof-Int1) = work_ndof(Int1:n_dof)
             work(colz:colz+n_dof-Int1) = work_ndof2(Int1:n_dof)

             goto 10

           else if (revcom.eq.precondLeft) then        ! perform the left preconditioning
                                                       ! work(colz) <-- M^{-1} * work(colx)
             call gmres_precondition(work(colx),work(colz),my_id,my_id_n,MPI_COMM_MASTER,MPI_COMM_N)
             goto 10

           else if (revcom.eq.precondRight) then       ! perform the right preconditioning

             if (my_id .eq. 0) call dcopy(n_dof,work(colx),Int1,work(colz),Int1)

             goto 10

           else if (revcom.eq.dotProd) then            ! perform the scalar product
                                                       ! work(colz) <-- work(colx) work(coly)

             if (my_id .eq. 0) call dgemv('C',n_dof,nbscal,ONE, work(colx),n_dof,work(coly),Int1,ZERO,work(colz),Int1)

             goto 10

           endif

    !******************************** end of GMRES reverse communication

    if (my_id .eq. 0) write(*,*) 'check work : ',n_dof,maxval(abs(work(1:n_dof)))

    if (my_id .eq. 0) deltas(1:n_dof) = work(1:n_dof)

    work_ndof(Int1:n_dof) = deltas(Int1:n_dof)
    work_ndof2(Int1:n_dof) = work(n_dof+Int1:2*n_dof)
    call gmres_matrix_vector(n_dof,work_ndof,n_dof,work_ndof2,my_id)
    deltas(Int1:n_dof) = work_ndof(Int1:n_dof)
    work(n_dof+Int1:2*n_dof) = work_ndof2(Int1:n_dof)

    if (my_id .eq. 0) then
      sum = 0.d0
      err = -1.d20
      Bnorm = 0.d0
      Xnorm = 0.d0
      do i=1,n_dof
        sum   = sum + (work(n_dof+i)-RHS_glob(i))**2
        err   = max(err,abs(work(n_dof+i)-RHS_glob(i)))
        Bnorm = Bnorm + RHS_glob(i)**2
        Xnorm = Xnorm + deltas(i)**2
      enddo
      write(*,'(A,4e16.8)') ' residu test after : ',sqrt(sum),err,sqrt(Bnorm),sqrt(Xnorm)
      call tr_locvnorms("gmres_rhs",rhs_glob,n_dof)
    endif

    iter_gmres = info(2) ! Actual number of iterations
    call tr_debug_writei("nbiter_gmres", iter_gmres)
    call MPI_BCAST(iter_gmres,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

    call tr_deallocate(work,"work",CAT_GMRES)
    call tr_deallocate(work_ndof,"work_ndof",CAT_GMRES)
    call tr_deallocate(work_ndof2,"work_ndof2",CAT_GMRES)

    call clck_time(t1)
    call clck_ldiff(t0,t1,tsecond)
    !if (my_id .eq. 0)  write(*,FMT_TIMING) my_id, '## Elapsed time gmres :', tsecond

    call r3_info_end (r3_info_index_0)  ! timing
    return

  end subroutine gmres_driver


!> Solve step of the local matrices for each toroidal harmonic (preconditioner for gmres)
  subroutine gmres_precondition(x,y,my_id,my_id_n,MPI_COMM_MASTER,MPI_COMM_N)

#ifdef USE_COMPLEX_PRECOND
    use real2complex_mod
    use global_distributed_matrix, only: rhs_cmplx_sol, n_cmplx
#endif
    use pastix_module
    use wsmp_module
    use phys_module, only: use_pastix, use_mumps, use_strumpack
    use preconditioner_module, only: my_row_index, my_row_factor

#if USE_PASTIX
#include "pastix_fortran.h"
#endif

#ifdef USE_STRUMPACK
    use strumpack_module, only: strumpack_solve
#endif
#ifdef USE_PASTIX6
    use mod_pastix, only: pastix_solve
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

#if defined(USE_PASTIX)
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

        pastix_iparm(IPARM_START_TASK)            = API_TASK_SOLVE
        pastix_iparm(IPARM_END_TASK)              = pastix_endsolve
        pastix_iparm(IPARM_RHS_MAKING)            = pastix_rhs                 ! right hand side (0 : use RHS)
        pastix_iparm(IPARM_SYM)                   = pastix_sym
        pastix_iparm(IPARM_AMALGAMATION_LEVEL)    = pastix_amalg

#ifdef USE_BLOCK
#if !defined(USE_COMPLEX_PRECOND)
        call pastix_fortran(pastix_data,MPI_COMM_N, n_block,                        &
                     DUMMY_INT, DUMMY_INT, DUMMY_REAL, &
                     pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,Int1,pastix_iparm,pastix_dparm)
#else
        call pastix_fortran(pastix_data,MPI_COMM_N, n_block,                        &
                     DUMMY_INT, DUMMY_INT, DUMMY_REAL, &
                     pastix_perm_vars,pastix_iperm_vars,rhs_cmplx_sol,Int1,pastix_iparm,pastix_dparm)
#endif /* !defined(USE_COMPLEX_PRECOND)  */

#else /* USE_BLOCK */
#if !defined(USE_COMPLEX_PRECOND)
        call pastix_fortran(pastix_data,MPI_COMM_N,mumps_par%n, DUMMY_INT, DUMMY_INT, DUMMY_REAL, &
               pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,Int1,pastix_iparm,pastix_dparm)
#else
        call pastix_fortran(pastix_data,MPI_COMM_N,n_cmplx, DUMMY_INT, DUMMY_INT, DUMMY_REAL, &
               pastix_perm_vars,pastix_iperm_vars,rhs_cmplx_sol,Int1,pastix_iparm,pastix_dparm)
#endif /* !defined(USE_COMPLEX_PRECOND)  */
#endif /* USE_BLOCK */

      endif
    endif ! use_pastix
#endif /* defined(USE_PASTIX) */

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

#ifdef USE_PASTIX6
    if (use_pastix) then
      call MPI_BCAST(mumps_par%rhs,mumps_par%n,MPI_DOUBLE_PRECISION,0,MPI_COMM_N,ierr)
      call pastix_solve(mumps_par%n,mumps_par%rhs)
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
#if !defined(USE_PASTIX6)
        do i=1,mumps_par%n
          mumps_par%rhs(i) =  mumps_par%rhs(i) / column_scaling(i)
        enddo
#endif
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

!> Sparse matrix vector product using coordinate scheme to be called only on all of MPI_COMM_WORLD
!!
!! result (y) is only known on id=0
  subroutine gmres_matrix_vector(size_x, x, size_y, y, my_id)

    implicit none

    logical, parameter :: PRINT_TIMING_INFO = .false.

    ! --- Routine parameters
    integer(kind=int_all), intent(in)  :: size_x, size_y
    real*8,                intent(in)  :: x(size_x)
    real*8,                intent(out) :: y(size_y)
    integer,               intent(in)  :: my_id

    ! --- Local variables
    real*8                :: t1, t2, t3, t4, t5
    real*8, allocatable   :: y_tmp(:)
    real*8                :: y_tmp_block(n_tor*n_var)
    integer,allocatable   :: recv_counts(:), recv_disp(:)
    integer               :: n, i, ir, jc, ierr, n_cpu
    integer(kind=int_all) :: n_blocksize, n_blocks, iA_start, ix_start, iy_start, ndof_local, index_offset, Int_tmp
    integer               :: ndof_glob_short
    integer index_ytmp_min,index_ytmp_max
    integer(kind=int_all), parameter :: Int1 = 1

    if ( PRINT_TIMING_INFO ) then
      call cpu_time(t1)
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
      call cpu_time(t2)
    end if

    call MPI_COMM_SIZE(MPI_COMM_WORLD, n_cpu, ierr) ! number of MPI ranks

    ndof_glob_short = ndof_glob
    call MPI_BCAST(x,ndof_glob_short,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

    if ( PRINT_TIMING_INFO ) call cpu_time(t3)

    y(1:size_y)  = 0.d0
    n_blocksize  = n_tor * n_var
    n_blocks     = nz_glob/n_blocksize**2
    ndof_local   = (local_index_end(my_id+1) - local_index_start(my_id+1) + 1) * n_blocksize
    index_offset = (local_index_start(my_id+1)-1) * n_blocksize
    call tr_allocate(y_tmp,Int1,ndof_local,"y_tmp",CAT_GMRES)
    y_tmp        = 0.d0

! --- The actual matrix vector multiplication uses dense matrix-vector products for the small
!     dense blocks within our sparse matrix. The size of these blocks depends on n_tor. Depending on
!     this block size (so depending on n_tor), two slightly different kernels are implemented.
!     Warning: this simply seg-faults with long-integers.
#ifndef INTSIZE64
    if ( n_tor <= 7 ) then

!$omp parallel default(none)                                                                     &
!$omp   shared(A_glob, jcn_glob, irn_glob, x, n_blocks, n_blocksize, index_offset)               &
!$omp   private(i,iA_start,ix_start, iy_start, ir, jc, y_tmp_block )                             &
!$omp   reduction(+:y_tmp)

      !$omp do schedule(guided)
      do i = 1, n_blocks

        iA_start = (i-1) * n_blocksize**2
        ix_start = jcn_glob(iA_start+1)
        iy_start = irn_glob(iA_start+1) - index_offset

        call dgemv('T',n_blocksize,n_blocksize,1.d0,A_glob(iA_start+1),n_blocksize,x(ix_start),Int1,0.d0,y_tmp_block,Int1)

        y_tmp(iy_start:iy_start+n_blocksize-1) = y_tmp(iy_start:iy_start+n_blocksize-1) + y_tmp_block(1:n_blocksize)

      end do
!$omp end do
!$omp end parallel

    else ! ... so in case n_tor is larger than 7
#endif

!$omp parallel default(none)                                                                     &
!$omp   shared(y_tmp, A_glob, jcn_glob, irn_glob, x, n_blocks, n_blocksize, index_offset)        &
!$omp   private(i,iA_start,ix_start, iy_start, ir, jc, y_tmp_block)

!$omp do schedule(guided)
      do i = 1, n_blocks

        iA_start = (i-1) * n_blocksize**2
        ix_start = jcn_glob(iA_start+1)
        iy_start = irn_glob(iA_start+1) - index_offset

        call dgemv('T',n_blocksize,n_blocksize,1.d0,A_glob(iA_start+1),n_blocksize,x(ix_start),Int1,0.d0,y_tmp_block,Int1)

!$omp critical
        y_tmp(iy_start:iy_start+n_blocksize-1) = y_tmp(iy_start:iy_start+n_blocksize-1) + y_tmp_block(1:n_blocksize)
        !$omp end critical

      end do
!$omp end do
!$omp end parallel

#ifndef INTSIZE64
    end if
! --- End: Two different kernels for matrix-vector multiplication depending on n_tor
#endif

! --- The unparallelized and unoptimized alternative for reference
    !do i=1,nz_glob
    !  ir = irn_glob(i)
    !  jc = jcn_glob(i)
    !  y_tmp(ir) = y_tmp(ir) + A_glob(i) * x(jc)
    !enddo

    if ( PRINT_TIMING_INFO ) call cpu_time(t4)

    y(1:size_y) = 0.d0

    call tr_allocate(recv_counts,1,n_cpu,"recv_counts",CAT_GMRES)
    call tr_allocate(recv_disp,1,n_cpu,"recv_disp",CAT_GMRES)

    do i = 1, n_cpu
       Int_tmp = (local_index_end(i) - local_index_start(i) + Int1) * n_blocksize
       recv_counts(i) = Int_tmp
    enddo

    recv_disp(1) = 0
    do i = 2, n_cpu
       recv_disp(i) = recv_disp(i-1) + recv_counts(i-1)
    enddo

    call mpi_gatherv(y_tmp,ndof_local,MPI_DOUBLE_PRECISION,y,recv_counts,recv_disp,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call tr_deallocate(y_tmp,"y_tmp",CAT_GMRES)
    call tr_deallocate(recv_counts,"recv_counts",CAT_GMRES)
    call tr_deallocate(recv_disp,"recv_disp",CAT_GMRES)

    if ( PRINT_TIMING_INFO ) then
      call cpu_time(t5)
      write(*,'(A,i3,3f14.6)') ' M-V timing  barrier: ',my_id,t2-t1
      write(*,'(A,i3,3f14.6)') ' M-V timing  bcast  : ',my_id,t3-t2
      write(*,'(A,i3,3f14.6)') ' M-V timing  dgemv  : ',my_id,t4-t3
      write(*,'(A,i3,3f14.6)') ' M-V timing  reduce : ',my_id,t5-t4
      write(*,'(A,i3,3f14.6)') ' M-V timing  TOTAL  : ',my_id,t5-t1
    end if

  end subroutine gmres_matrix_vector



end module mod_gmres
