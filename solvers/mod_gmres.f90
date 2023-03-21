module mod_gmres
#ifdef USE_GMRES

  use mod_integer_types
  use mod_clock

  use data_structure, only: type_SP_MATRIX, type_RHS
  use phys_module,    only: gmres_4, gmres_m, index_now
  

  private
  public :: gmres_driver


  contains

!> Driver for the reverse communication GMRES routine from dPackgmres (CERFACS)
  subroutine gmres_driver(a_mat, rhs_vec, sol_vec, solver)
    
    use mod_dpackgmres, only: init_dgmres, drive_dgmres
    use tr_module,      only: tr_locvnorms, tr_vdump
    use mod_sparse_data,     only: type_SP_SOLVER

    implicit none
#include "r3_info.h"
    type(type_SP_MATRIX)  :: a_mat
    type(type_RHS)        :: rhs_vec, sol_vec
    type(type_SP_SOLVER)  :: solver

    integer               :: j
    integer               :: my_id, comm, ierr, counts
    integer               :: revcom, iter_gmres
    integer               :: icntl(8), info(3)
    integer(kind=int_all) :: i, m, colx, coly, colz, nbscal, lwork, irc(5), n_dof

    integer :: matvec, precondLeft, precondRight, dotProd
    real*8  :: cntl(5), rinfo(2), sum, err, Bnorm, Xnorm
    real*8, allocatable          :: work(:), work_ndof(:), work_ndof2(:)
    type(clcktype)               :: t0, t1
    real*8                       :: tsecond
    real*8, allocatable          :: rhs_tmp(:)
    real*8 ::ZERO, ONE
    CHARACTER(LEN=128) :: fname
    
    parameter (ZERO = 0.0d0, ONE = 1.0d0)
    parameter (matvec=1, precondLeft=2, precondRight=3, dotProd=4)
    
    my_id      = solver%pc%my_id
    comm       = solver%pc%comm
    iter_gmres = solver%iter_gmres

    !write(*,*) ' GMRES DRIVER : ',my_id,my_id_n
    call r3_info_begin (r3_info_index_0, 'gmres_driver')  ! timing
    call clck_time(t0)
    call init_dgmres(icntl,cntl)
    

    icntl(3) = 0            ! output unit
    if (my_id.eq.0) icntl(3) = 6
    icntl(7) = iter_gmres   ! Maximum number of iterations
    icntl(4) = 1            ! preconditioner (1) = left preconditioner
    icntl(5) = 3            ! orthogonalization scheme
    icntl(6) = 1            ! initial guess  (1) = user supplied guess
    icntl(8) = 1            ! residual calculation strategy at restart

    cntl(1) = solver%iter_tol         ! stopping tolerance
    cntl(2) = 1.d0
    cntl(3) = 1.d0
    cntl(4) = gmres_4          ! 1.d0
    cntl(5) = 1.d0

    m = gmres_m

    n_dof = rhs_vec%n

    lwork = m*m + m*(n_dof+5) + 6*n_dof + m + 1

    allocate(work(lwork))
    allocate(work_ndof(n_dof))
    allocate(work_ndof2(n_dof))

    work(1:n_dof)         = sol_vec%val(1:n_dof)                     ! the initial guess
    work(n_dof+1:2*n_dof) = rhs_vec%val(1:n_dof)                   ! the right hand side

    work_ndof(Int1:n_dof)  = work(Int1:n_dof)
    work_ndof2(Int1:n_dof) = work(2*n_dof+Int1:3*n_dof)

    call gmres_matrix_vector(n_dof,work_ndof,n_dof,work_ndof2,a_mat)

    work(Int1:n_dof)           = work_ndof(Int1:n_dof)
    work(2*n_dof+Int1:3*n_dof) = work_ndof2(Int1:n_dof)


    sum = 0.d0
    err = -1.d20
    Bnorm = 0.d0
    Xnorm = 0.d0
    do i=1,n_dof
      sum = sum + (work(2*n_dof+i)-work(n_dof+i))**2
      err = max(err,abs(work(2*n_dof+i)-work(n_dof+i)))
      Bnorm = Bnorm + rhs_vec%val(i)**2
      Xnorm = Xnorm + sol_vec%val(i)**2
    enddo

    if (my_id.eq.0) write(*,'(A,4e16.8)') ' residu test before : ',sqrt(sum),err,sqrt(Bnorm),sqrt(Xnorm)

    if (my_id .eq. 0) then
       write(fname,'(A,I6.6,A1,I6.6)')"product_before",index_now
       call tr_vdump(fname, work(2*n_dof+1:3*n_dof), n_dof)
       write(fname,'(A,I6.6,A1,I6.6)')"X_before",index_now
       call tr_vdump(fname, rhs_vec%val, n_dof)
    end if

    !*****************************************
    !** Reverse communication implementation
    !*****************************************

    10     call MPI_barrier(MPI_COMM_WORLD,ierr)

           call drive_dgmres(n_dof,n_dof,m,lwork,work,irc,icntl,cntl,info,rinfo)

           call MPI_BCAST(irc,5,MPI_INTEGER_ALL,0,a_mat%comm,ierr)
           revcom = irc(1)
           colx   = irc(2)
           coly   = irc(3)
           colz   = irc(4)
           nbscal = irc(5)

           if (revcom.eq.matvec) then                  ! perform the matrix vector product
                                                       ! work(colz) <-- A * work(colx)

             work_ndof(Int1:n_dof) = work(colx:colx+n_dof-Int1)
             work_ndof2(Int1:n_dof) = work(colz:colz+n_dof-Int1)
             call gmres_matrix_vector(n_dof,work_ndof,n_dof,work_ndof2,a_mat)
             work(colx:colx+n_dof-Int1) = work_ndof(Int1:n_dof)
             work(colz:colz+n_dof-Int1) = work_ndof2(Int1:n_dof)

             goto 10

           else if (revcom.eq.precondLeft) then        ! perform the left preconditioning
                                                       ! work(colz) <-- M^{-1} * work(colx)
             call gmres_precondition(work(colx), work(colz), n_dof, solver)
             goto 10

           else if (revcom.eq.precondRight) then       ! perform the right preconditioning

             call dcopy(n_dof,work(colx),Int1,work(colz),Int1)

             goto 10

           else if (revcom.eq.dotProd) then            ! perform the scalar product
                                                       ! work(colz) <-- work(colx) work(coly)

             call dgemv('C',n_dof,nbscal,ONE, work(colx),n_dof,work(coly),Int1,ZERO,work(colz),Int1)

             goto 10

           endif

    !******************************** end of GMRES reverse communication

    if (my_id .eq. 0) write(*,*) 'check work : ',n_dof,maxval(abs(work(1:n_dof)))

    sol_vec%val(1:n_dof) = work(1:n_dof)

    work_ndof(Int1:n_dof) = sol_vec%val(Int1:n_dof)
    work_ndof2(Int1:n_dof) = work(n_dof+Int1:2*n_dof)
    call gmres_matrix_vector(n_dof,work_ndof,n_dof,work_ndof2, a_mat)
    sol_vec%val(Int1:n_dof) = work_ndof(Int1:n_dof)
    work(n_dof+Int1:2*n_dof) = work_ndof2(Int1:n_dof)

    sum = 0.d0
    err = -1.d20
    Bnorm = 0.d0
    Xnorm = 0.d0
    do i=1,n_dof
      sum   = sum + (work(n_dof+i)-rhs_vec%val(i))**2
      err   = max(err,abs(work(n_dof+i)-rhs_vec%val(i)))
      Bnorm = Bnorm + rhs_vec%val(i)**2
      Xnorm = Xnorm + sol_vec%val(i)**2
    enddo
    if (my_id.eq.0) write(*,'(A,4e16.8)') ' residu test after : ',sqrt(sum),err,sqrt(Bnorm),sqrt(Xnorm)
    call tr_locvnorms("gmres_rhs",rhs_vec%val,n_dof)

    iter_gmres = info(2) ! Actual number of iterations
    call MPI_BCAST(iter_gmres,1,MPI_INTEGER,0,a_mat%comm,ierr)
    
    solver%iter_gmres = info(2)

    deallocate(work)
    deallocate(work_ndof)
    deallocate(work_ndof2)

    call clck_time(t1)
    call clck_ldiff(t0,t1,tsecond)
    !if (my_id .eq. 0)  write(*,FMT_TIMING) my_id, '## Elapsed time gmres :', tsecond

    call r3_info_end (r3_info_index_0)  ! timing
    return

  end subroutine gmres_driver

!> Sparse matrix vector product using coordinate scheme to be called only on all of MPI_COMM_WORLD
!! result (y) is only known on id=0
  subroutine gmres_matrix_vector(size_x, x, size_y, y, a_mat)
    use mod_parameters, only: n_tor

    implicit none
    
    type(type_SP_MATRIX)  :: a_mat

    logical, parameter :: PRINT_TIMING_INFO = .false.

    ! --- Routine parameters
    integer(kind=int_all), intent(in)  :: size_x, size_y
    real*8,                intent(in)  :: x(size_x)
    real*8,                intent(out) :: y(size_y)

    ! --- Local variables
    real*8                :: t1, t2, t3, t4, t5
    real*8, allocatable   :: y_tmp(:)
    real*8                :: y_tmp_block(a_mat%block_size)
    integer,allocatable   :: recv_counts(:), recv_disp(:)
    integer               :: n, i, ir, jc
    integer               :: my_id, n_cpu, ierr, counts
    integer(kind=int_all) :: n_blocksize, n_blocks, iA_start, ix_start, iy_start, index_offset, Int_tmp
    integer               :: ndof_local
    
    if ( PRINT_TIMING_INFO ) then
      call cpu_time(t1)
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
      call cpu_time(t2)
    end if

    call MPI_COMM_SIZE(a_mat%comm, n_cpu, ierr)
    call MPI_COMM_RANK(a_mat%comm, my_id, ierr) 

    counts = size_x
    call MPI_BCAST(x,counts,MPI_DOUBLE_PRECISION,0,a_mat%comm,ierr)

    if ( PRINT_TIMING_INFO ) call cpu_time(t3)

    y(1:size_y)  = 0.d0
    n_blocksize  = a_mat%block_size
    n_blocks     = a_mat%nnz/n_blocksize**2
    ndof_local   = (a_mat%index_max(my_id + 1) - a_mat%index_min(my_id + 1) + 1)*n_blocksize
    
    index_offset = (a_mat%index_min(my_id + 1) - 1)*n_blocksize
    
    allocate(y_tmp(ndof_local))
    y_tmp        = 0.d0
    
! --- The actual matrix vector multiplication uses dense matrix-vector products for the small
!     dense blocks within our sparse matrix. The size of these blocks depends on n_tor. Depending on
!     this block size (so depending on n_tor), two slightly different kernels are implemented.
!     Warning: this simply seg-faults with long-integers.
#ifndef INTSIZE64
    if ( n_tor <= 7 ) then

!$omp parallel default(none) &
!$omp   shared(a_mat, x, n_blocks, n_blocksize, index_offset) &
!$omp   private(i, iA_start, ix_start, iy_start, ir, jc, y_tmp_block ) &
!$omp   reduction(+:y_tmp)
!$omp do schedule(guided)
      do i = 1, n_blocks

        iA_start = (i-1) * n_blocksize**2
        ix_start = a_mat%jcn(iA_start+1)
        iy_start = a_mat%irn(iA_start+1) - index_offset

        call dgemv('T', n_blocksize, n_blocksize, 1.d0, a_mat%val(iA_start+1), n_blocksize, x(ix_start), Int1, 0.d0, y_tmp_block, Int1)

        y_tmp(iy_start:iy_start+n_blocksize-1) = y_tmp(iy_start:iy_start+n_blocksize-1) + y_tmp_block(1:n_blocksize)

      end do
!$omp end do
!$omp end parallel
    else ! ... so in case n_tor is larger than 7
#endif

!$omp parallel default(none) &
!$omp   shared(y_tmp, a_mat, x, n_blocks, n_blocksize, index_offset)  &
!$omp   private(i,iA_start,ix_start, iy_start, ir, jc, y_tmp_block)
!$omp do schedule(guided)
      do i = 1, n_blocks

        iA_start = (i-1) * n_blocksize**2
        ix_start = a_mat%jcn(iA_start+1)
        iy_start = a_mat%irn(iA_start+1) - index_offset

        call dgemv('T',n_blocksize,n_blocksize,1.d0,a_mat%val(iA_start+1),n_blocksize,x(ix_start),Int1,0.d0,y_tmp_block,Int1)

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

    allocate(recv_counts(n_cpu))
    allocate(recv_disp(n_cpu))

    do i = 1, n_cpu
       int_tmp = (a_mat%index_max(i) - a_mat%index_min(i) + 1)*n_blocksize       
       recv_counts(i) = int_tmp
    enddo

    recv_disp(1) = 0
    do i = 2, n_cpu
       recv_disp(i) = recv_disp(i-1) + recv_counts(i-1)
    enddo

    call mpi_allgatherv(y_tmp, ndof_local, MPI_DOUBLE_PRECISION, y, recv_counts, recv_disp, MPI_DOUBLE_PRECISION, a_mat%comm,ierr)
    
    deallocate(y_tmp)
    deallocate(recv_counts)
    deallocate(recv_disp)

    if ( PRINT_TIMING_INFO ) then
      call cpu_time(t5)
      write(*,'(A,i3,3f14.6)') ' M-V timing  barrier: ',my_id,t2-t1
      write(*,'(A,i3,3f14.6)') ' M-V timing  bcast  : ',my_id,t3-t2
      write(*,'(A,i3,3f14.6)') ' M-V timing  dgemv  : ',my_id,t4-t3
      write(*,'(A,i3,3f14.6)') ' M-V timing  reduce : ',my_id,t5-t4
      write(*,'(A,i3,3f14.6)') ' M-V timing  TOTAL  : ',my_id,t5-t1
    end if
    
    return

  end subroutine gmres_matrix_vector
  
!> Solve step of the local matrices for each toroidal harmonic (preconditioner for gmres)
  subroutine gmres_precondition(x, y, n_dof, solver)

    use phys_module, only: use_pastix, use_mumps, use_strumpack
    use mod_sparse_data,     only: type_SP_SOLVER

#ifdef USE_STRUMPACK
    use mod_strumpack, only: strumpack_solve
#endif
#if defined(USE_PASTIX)||defined(USE_PASTIX6)
    use mod_pastix, only: pastix_solve
#endif
#ifdef USE_MUMPS
    use mod_mumps, only: mumps_solve
#endif

    implicit none
    
    type(type_SP_SOLVER)  :: solver
    
    real, intent (in)     :: x(*)
    real, intent (inout)  :: y(*)

    real, allocatable     :: y_tmp(:)
    integer               :: i, ierr, counts
    real, dimension(:), allocatable :: y_dum

    real*8                :: DUMMY_REAL(1:1)
    integer(kind=int_all) :: DUMMY_INT (1:1)
    integer(kind=int_all) :: n_dof
    
    counts = n_dof
    call MPI_BCAST(x,counts,MPI_DOUBLE_PRECISION,0,solver%pc%MPI_COMM_N,ierr)
   
    do i = 1, solver%pc%rhs%n
      solver%pc%rhs%val(i) = x(solver%pc%row_index(i))
    enddo

#ifdef USE_MUMPS
    if (use_mumps) then
    endif
#endif

#if defined(USE_PASTIX)||defined(USE_PASTIX6)
    if (use_pastix) then
      call pastix_solve(solver%ptss, solver%pc%rhs)
    endif ! use_pastix
#endif /* defined(USE_PASTIX) */

#ifdef USE_STRUMPACK
    if (use_strumpack) then
      call strumpack_solve(solver%spss, solver%pc%rhs)
    endif
#endif

#ifdef USE_MUMPS
    if (use_mumps) then
      call mumps_solve(solver%mmss, solver%pc%rhs)
    endif
#endif

    allocate(y_dum(n_dof))
    y_dum = 0.d0

    if (solver%pc%my_id_n .eq. 0) then      
    ! put local solution into global vector y
      do i = 1, solver%pc%rhs%n
        y_dum(solver%pc%row_index(i)) = solver%pc%rhs%val(i)*solver%pc%row_factor
      enddo
    endif

    counts = n_dof
    call MPI_AllReduce(MPI_IN_PLACE,y_dum,counts,MPI_DOUBLE_PRECISION,MPI_SUM,solver%pc%comm,ierr)
    y(1:n_dof) = y_dum(1:n_dof)

    deallocate(y_dum)


    return
  end subroutine gmres_precondition


#endif
end module mod_gmres
