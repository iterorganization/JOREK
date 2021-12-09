!> Module for the harmonic matrix systems
!!
!! * This library provides routines for:
!!   - Initialize the solvers, analyzes the matrix, factorizes the matrix
!!   - Solve the matrix system to obtain initial conditions for the time step
!!   - The LU factorized matrix stored internally by the solvers is re-used in gmres_precondition
!! * Supports the following solver libraries:
!!   - PaSTiX 5.x (real or complex)
!!   - MUMPS
!!   - WSMP (not tested for a while)
!!   - STRUMPACK
!!   - Typical choices for production are PaStiX 5.x or STRUMPACK
module solve_mat_n
  use phys_module, only: use_mumps, use_pastix, use_strumpack, use_wsmp
  use preconditioner_module, only : my_row_index, my_row_factor, my_mode_set_n
  use matio_module, only: timestamp
  use mod_integer_types

  implicit none        



contains

  !> Routine for the binding of threads to cores in PaStiX
  subroutine pastix_bind_threads(my_id)

    use pastix_module
   
    implicit none

#ifdef USE_PASTIX
#include "pastix_fortran.h"
#else
#include "no_pastix_fortran.h"
#endif


    integer, intent(in) :: my_id

    integer*4, dimension(1:pastix_nthrd) :: thread_map

    integer*4 k, packsize, procpernode

#if defined(WORLDWAR2) && defined(CORES_PER_NODE)
    procpernode = CORES_PER_NODE/pastix_nthrd
    packsize = CORES_PER_NODE/procpernode 
!    if (my_id .eq. 0) print *, "packsize", packsize, "procpernode", procpernode
    do k = 1, pastix_nthrd
      thread_map(k) = mod(my_id * packsize,CORES_PER_NODE) + k-1
    end do

    call pastix_fortran_bindthreads(pastix_data, pastix_nthrd, thread_map(1:))

#endif
  end subroutine pastix_bind_threads
  
  
#if defined(USE_PASTIX) || defined(USE_MUMPS)    
  !> Solves the system of equation for each harmonic using mumps, pastix, or wsmp
  subroutine solve_matrix_n(my_id,MPI_COMM_N,MPI_COMM_MASTER,solve_only)

#ifdef USE_COMPLEX_PRECOND
    use real2complex_mod
#endif
    use tr_module
    use mod_parameters
    use mumps_module
    use wsmp_module
    use pastix_module
    use global_distributed_matrix
    use mpi_mod 
    use mod_clock
    use phys_module, only : index_now, use_BLR_compression, epsilon_BLR, just_in_time_BLR, pastix_blr_abs_tol
    use mod_coicsr
    use mod_integer_types
   
    implicit none

#ifdef USE_PASTIX
#include "pastix_fortran.h"
#endif

#include "r3_info.h"
    
    ! --- Routine parameters
    integer, intent(in) :: my_id
    integer, intent(in) :: MPI_COMM_N, MPI_COMM_MASTER
    logical, intent(in) :: solve_only
   
    ! --- Local variables
    integer               :: my_id_n, n_cpu_n, ierr
    integer(kind=int_all) :: i, j, k
    type(clcktype)        :: t_itstart, t0, t1, t2, t3
    real*8                :: tsecond
    real*8, allocatable   :: RHS_tmp(:)

    !Split broadcast
    character*8           :: type
    real*8                :: DUMMY_REAL(1:1)
    integer(kind=int_all) :: DUMMY_INT (1:1)
    integer(kind=int_all), parameter   :: Int1=1
    CHARACTER(LEN=128) :: fname

    call r3_info_begin (r3_info_index_0, 'solve_matrix_n')                  ! timing
    call tr_print_memsize("BeforeSolveN")
    call tr_debug_write("smn_A_mumps_par%n",mumps_par%n)

    if (my_id .eq. 0) then
      write(*,*) my_id,'*********************************'
      write(*,*) my_id,'*      solve local matrix  (n)  *'
      write(*,*) my_id,'*********************************'

      if (use_mumps)  write(*,*) my_id,'*       using solver MUMPS      *'
#ifndef USE_COMPLEX_PRECOND
      if (use_pastix) write(*,*) my_id,'*       using solver PastiX     *'
#else 
      if (use_pastix) write(*,*) my_id,'*  using complex solver PastiX  *'
#endif
      if (use_wsmp)   write(*,*) my_id,'*       using solver WSMP       *'

      write(*,*) my_id,'*********************************'
    endif

    call MPI_COMM_RANK(MPI_COMM_N, my_id_n, ierr)     ! the id of each cpu
    call MPI_COMM_SIZE(MPI_COMM_N, n_cpu_n, ierr)     ! the number of cpus

    NOTSOLVEONLY: if (.not. solve_only) then
      ! This part of the code is needed when the preconditioning matrix is updated. Otherwise, the
      ! the LU decomposed matrix of a previous time step is re-used.

      ! --- Column scaling -------------------------------------------------------------------------
      if (my_id_n .eq. 0) then
        if (allocated(column_scaling))  call tr_deallocate(column_scaling,"column_scaling",CAT_DMATRIX)
        call tr_allocate(column_scaling,Int1,mumps_par%N,"column_scaling",CAT_DMATRIX)

        column_scaling = 1.d-20
        do k=1,mumps_par%nz
          j = mumps_par%jcn(k)
          column_scaling(j) = min(max(column_scaling(j),abs(mumps_par%A(k))),1d20)
        enddo
        if (my_id .eq. 0) then
           write(fname,'(A,I6.6)')"column_scaling",index_now
           call tr_vdump(fname,column_scaling,mumps_par%N)
        end if
        !CALL MPI_Abort(MPI_COMM_WORLD, 1, ierr)
        write(*,'(2i4,A,2e12.4)') my_id,my_id_n,' COLUMN SCALING : ',minval(column_scaling),maxval(column_scaling)
        do k=1,mumps_par%nz
          j = mumps_par%jcn(k)
          mumps_par%A(k) = mumps_par%A(k) / column_scaling(j)
        enddo
        
        ! elapsed time analysis start
        call MPI_Barrier(MPI_COMM_MASTER,ierr)
        call clck_time(t0)
      endif
      ! --- End column scaling ---------------------------------------------------------------------


      if (use_mumps) then
#ifdef USE_MUMPS

        mumps_par%JOB = 1                                  ! Analysis, only needed when grid has changed

        mumps_par%icntl(7)  = mumps_ordering               ! ordering option (7:automatic, 3:Scotch, 4:PORD, 5:METIS), default: 7
        mumps_par%icntl(8)  = 7                            ! row and column scaling
        mumps_par%icntl(14) = 30                           ! MAXS
        mumps_par%icntl(18) = 0
        if (use_BLR_compression) then
          mumps_par%icntl(35) = 1                          ! Block-low-rank (BLR) compression. 0: off (default), 1: automatic, 2: factorisation and solution, 3: only factorisation
          mumps_par%cntl(7)   = epsilon_BLR                ! Accuracy of BLR approximation
        endif

        call DMUMPS(mumps_par)

#endif /* ifdef USE_MUMPS */
      else ! .not. use_mumps --> use_pastix or use_wsmp

        if (my_id_n .eq. 0) then           

          ! Timing
          call MPI_Barrier(MPI_COMM_MASTER,ierr)
          call clck_time(t2)

#ifdef USE_BLOCK
          ! --- Preparation for USE_BLOCK ----------------------------------------------------------
          ! - reduce IRN,JCN to make use of blocksize ntor*nvar
          !   (temporary solution before using blocks everywhere)
          ! - convert matrix and rhs to complex if necessary
          ! - convert row/column sparse matrix to CSR format
#ifndef USE_COMPLEX_PRECOND
          block_size = n_var*my_mode_set_n
          block_size2 = block_size**2
          n_block   = mumps_par%n  / block_size
          nnz_block = mumps_par%nz / block_size2

          do i=1,nnz_block  
            mumps_par%irn(i) = (mumps_par%irn((i-1)*block_size2+1) - 1) / block_size + 1 
            mumps_par%jcn(i) = (mumps_par%jcn((i-1)*block_size2+1) - 1) / block_size + 1 
          enddo

          if (allocated(sparskit_work)) deallocate(sparskit_work)
          allocate(sparskit_work(n_block+1))
          call coicsr2(n_block,nnz_block,mumps_par%A,mumps_par%IRN(1:nnz_block),mumps_par%JCN(1:nnz_block),block_size,sparskit_work)
#else /* ifndef USE_COMPLEX_PRECOND */
          block_size  = n_var
          block_size2 = block_size**2

          !-- converting real harmonic blocks into the complex ones
          call real2complex_a(my_id, my_id_n) 
          !-- converting RHS into the complex form
          call real2complex_rhs(my_id, my_id_n, rhs_cmplx) 
  
          n_block   = n_cmplx  / block_size
          nnz_block = nz_cmplx / block_size2
           
          do i=1,nnz_block  
            irn_cmplx(i) = (irn_cmplx((i-1)*block_size2+1) - 1) / block_size + 1 
            jcn_cmplx(i) = (jcn_cmplx((i-1)*block_size2+1) - 1) / block_size + 1 
          enddo

          if (allocated(sparskit_work)) deallocate(sparskit_work)
          allocate(sparskit_work(n_block+1))
          call coicsr2_cmplx(n_block,nnz_block,A_cmplx,irn_cmplx(1:nnz_block),jcn_cmplx(1:nnz_block),block_size,sparskit_work)
#endif /* ifndef USE_COMPLEX_PRECOND */

          ! WARNING:  USE_BLOCK does not (yet) work with WSMP!!!
          if (use_wsmp) then
#ifdef USE_WSMP
            call PWGSMP__allocate(n_block, nnz_block, my_id_n)
            call PWGSMP__initialize_matrix(n_block, nnz_block,                                     &
              mumps_par%a, mumps_par%jcn, mumps_par%irn, my_id_n )
#endif
          endif
          ! --- End Preparation for USE_BLOCK ------------------------------------------------------

#else /* ifdef USE_BLOCK */

          ! --- Preparation without USE_BLOCK ------------------------------------------------------
          ! - convert matrix and rhs to complex if necessary
          ! - convert row/column sparse matrix to CSR format
#ifndef USE_COMPLEX_PRECOND
          if (allocated(sparskit_work)) deallocate(sparskit_work)
          allocate(sparskit_work(mumps_par%N + 1))

          call coicsr(mumps_par%N,mumps_par%NZ,1,mumps_par%A,mumps_par%IRN,mumps_par%JCN,sparskit_work)
#else
          !-- converting real harmonic blocks into the complex ones
          call real2complex_a(my_id, my_id_n) 
          !-- converting RHS into the complex form
          call real2complex_rhs(my_id, my_id_n, rhs_cmplx) 

          if (allocated(sparskit_work)) deallocate(sparskit_work)
          allocate(sparskit_work(n_cmplx + 1))

          call coicsr_cmplx(n_cmplx,nz_cmplx,1,A_cmplx,irn_cmplx,jcn_cmplx,sparskit_work)
#endif /* ifndef USE_COMPLEX_PRECOND */
          if (use_wsmp) then
#ifdef USE_WSMP
            call PWGSMP__allocate(mumps_par%N, mumps_par%NZ, my_id_n)
            call PWGSMP__initialize_matrix(mumps_par%N, mumps_par%NZ,                              &
              mumps_par%a, mumps_par%jcn, mumps_par%irn, my_id_n )
#endif
          endif
          ! --- End Preparation without USE_BLOCK --------------------------------------------------
#endif /* ifdef USE_BLOCK */

          if (allocated(sparskit_work)) deallocate(sparskit_work)

          ! Timing
          call MPI_Barrier(MPI_COMM_MASTER,ierr) 
          call clck_time(t3)
          call clck_ldiff(t2,t3,tsecond)
          write(*,FMT_TIMING) my_id, '### Elapsed time coicsr :', tsecond

        else  ! (my_id_n > 0) below
#ifdef USE_WSMP
          if (use_wsmp) call PWGSMP__allocate(0, 0, my_id_n)
#endif
        endif ! end (my_id_n .eq. 0)

        ! --- Distribute data to the MPI "slave" tasks (>0) ----------------------------------------
        !     (When using WSMP, this is *not necessary* in 0-master mode!)
        if ((.not. use_wsmp).and.(.not. pastix_smp_only)) then
!          call pastix_init_num_threads(my_id)

#ifndef USE_COMPLEX_PRECOND
          call MPI_BCAST(mumps_par%n,1,MPI_INTEGER_ALL,0,MPI_COMM_N,ierr)
          call MPI_BCAST(mumps_par%nz,1,MPI_INTEGER_ALL,0,MPI_COMM_N,ierr)
#else
          call MPI_BCAST(n_cmplx,1,MPI_INTEGER_ALL,0,MPI_COMM_N,ierr)
          call MPI_BCAST(nz_cmplx,1,MPI_INTEGER_ALL,0,MPI_COMM_N,ierr)
#endif 

#ifdef USE_BLOCK
          call MPI_BCAST(block_size,1,MPI_INTEGER_ALL,0,MPI_COMM_N,ierr)
          call MPI_BCAST(n_block,1,MPI_INTEGER_ALL,0,MPI_COMM_N,ierr)
          call MPI_BCAST(nnz_block,1,MPI_INTEGER_ALL,0,MPI_COMM_N,ierr)
#endif
          if (my_id_n .gt. 0) then
#ifndef USE_COMPLEX_PRECOND
            if (associated(mumps_par%irn)) call tr_deallocatep(mumps_par%irn,"mumps_par%irn",CAT_DMATRIX)
            if (associated(mumps_par%jcn)) call tr_deallocatep(mumps_par%jcn,"mumps_par%jcn",CAT_DMATRIX)
            if (associated(mumps_par%A))   call tr_deallocatep(mumps_par%A,"mumps_par%A",CAT_DMATRIX)
            if (associated(mumps_par%rhs)) call tr_deallocatep(mumps_par%rhs,"mumps_par%rhs",CAT_DMATRIX)
            call tr_allocatep(mumps_par%irn,Int1,mumps_par%nz,"mumps_par%irn",CAT_DMATRIX)
            call tr_allocatep(mumps_par%jcn,Int1,mumps_par%nz,"mumps_par%jcn",CAT_DMATRIX)
            call tr_allocatep(mumps_par%a,Int1,mumps_par%nz,"mumps_par%a",CAT_DMATRIX)
            call tr_allocatep(mumps_par%rhs,Int1,mumps_par%n,"mumps_par%rhs",CAT_DMATRIX)
#else
            if (allocated(A_cmplx))  deallocate(A_cmplx)
            if (allocated(rhs_cmplx))  deallocate(rhs_cmplx)
            if (allocated(irn_cmplx))deallocate(irn_cmplx)
            if (allocated(jcn_cmplx))deallocate(jcn_cmplx) 
            allocate(rhs_cmplx(1:n_cmplx))
            allocate(A_cmplx(1:nz_cmplx))
            allocate(irn_cmplx(1:nz_cmplx))
            allocate(jcn_cmplx(1:nz_cmplx))
#endif /* ifndef USE_COMPLEX_PRECOND */
          endif

#ifdef USE_COMPLEX_PRECOND
          call MPI_BCAST(irn_cmplx,nz_cmplx,MPI_INTEGER_ALL,0,MPI_COMM_N,ierr)
          call MPI_BCAST(jcn_cmplx,nz_cmplx,MPI_INTEGER_ALL,0,MPI_COMM_N,ierr)
          call MPI_BCAST(A_cmplx,nz_cmplx,MPI_DOUBLE_COMPLEX,0,MPI_COMM_N,ierr)
          call MPI_BCAST(rhs_cmplx,n_cmplx,MPI_DOUBLE_COMPLEX,0,MPI_COMM_N,ierr)
#else
          ! Split MPI_BCAST if MPI buffer beyond 2Go
          type='intIRN'
          call split_broadcast(type,MPI_COMM_N)
          type='intJCN'
          call split_broadcast(type,MPI_COMM_N)
          type='double'
          call split_broadcast(type,MPI_COMM_N)
#endif /* ifdef USE_COMPLEX_PRECOND */

        endif
        ! --- End distribute data to the MPI "slave" tasks (>0) ------------------------------------
        
        ! --- Initialize the PaStiX solver ---------------------------------------------------------
        ! - set parameters
        ! - call initialization routine
        if  (.not. pastix_initialised)  then

          if ( (.not. pastix_smp_only) .or. (pastix_smp_only .and. (my_id_n .eq.0)) ) then

            if (use_pastix) then

              call pastix_init_num_threads(my_id)

              pastix_iparm(IPARM_MODIFY_PARAMETER) = API_NO         ! insert default values
              pastix_iparm(IPARM_START_TASK)       = API_TASK_INIT  ! initializse
              pastix_iparm(IPARM_END_TASK)         = API_TASK_INIT

#ifndef USE_COMPLEX_PRECOND 
              if (.not. pastix_smp_only) call MPI_BCAST(mumps_par%n,1,MPI_INTEGER_ALL,0,MPI_COMM_N,ierr)
#else
              if (.not. pastix_smp_only) call MPI_BCAST(n_cmplx,1,MPI_INTEGER_ALL,0,MPI_COMM_N,ierr)
#endif

#ifdef USE_BLOCK
              call tr_allocate(pastix_perm_vars,Int1,n_block,"pastix_perm_vars",CAT_UNKNOWN)
              call tr_allocate(pastix_iperm_vars,Int1,n_block,"pastix_iperm_vars",CAT_UNKNOWN)

#ifndef USE_COMPLEX_PRECOND
              call pastix_fortran(pastix_data,MPI_COMM_N,n_block,mumps_par%jcn,mumps_par%irn,mumps_par%A, &
                pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,Int1,pastix_iparm,pastix_dparm)
#else /* ifndef USE_COMPLEX_PRECOND */
              call pastix_fortran(pastix_data,MPI_COMM_N,n_block,jcn_cmplx,irn_cmplx,A_cmplx, &
                pastix_perm_vars,pastix_iperm_vars,rhs_cmplx,Int1,pastix_iparm,pastix_dparm)
#endif /* ifndef USE_COMPLEX_PRECOND */
#else /* ifdef USE_BLOCK */
#ifndef USE_COMPLEX_PRECOND
              call tr_allocate(pastix_perm_vars,Int1 ,mumps_par%n,"pastix_perm_vars",CAT_UNKNOWN)
              call tr_allocate(pastix_iperm_vars,Int1,mumps_par%n,"pastix_iperm_vars",CAT_UNKNOWN)

              call pastix_fortran(pastix_data,MPI_COMM_N,mumps_par%n,mumps_par%jcn,mumps_par%irn,mumps_par%A, &
                pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,1,pastix_iparm,pastix_dparm)
#else /* ifndef USE_COMPLEX_PRECOND */
              call tr_allocate(pastix_perm_vars,Int1 ,n_cmplx,"pastix_perm_vars",CAT_UNKNOWN)
              call tr_allocate(pastix_iperm_vars,Int1,n_cmplx,"pastix_iperm_vars",CAT_UNKNOWN)

              call pastix_fortran(pastix_data,MPI_COMM_N,n_cmplx,jcn_cmplx,irn_cmplx,A_cmplx, &
                pastix_perm_vars,pastix_iperm_vars,rhs_cmplx,Int1,pastix_iparm,pastix_dparm)
#endif /* ifndef USE_COMPLEX_PRECOND */

#endif /* ifdef USE_BLOCK */


              ! pastix input parameters working in Pastix5 and Pastix6
              pastix_iparm(IPARM_VERBOSE)               = pastix_verb              
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

              pastix_iparm(IPARM_RHS_MAKING)            = pastix_rhs                 ! right hand side (0 : use RHS)
              pastix_iparm(IPARM_SYM)                   = pastix_sym
              pastix_iparm(IPARM_AMALGAMATION_LEVEL)    = pastix_amalg
#ifdef WORLDWAR2
#ifdef FUNNELED
              pastix_iparm(IPARM_THREAD_COMM_MODE)      = API_THREAD_FUNNELED
#else
              pastix_iparm(IPARM_THREAD_COMM_MODE)      = API_THREAD_MULTIPLE
#endif
#endif /* ifdef WORLDWAR2 */


             else if (use_wsmp) then
#ifdef USE_WSMP
              call PWGSMP__initialize_solver(my_id_n, MPI_COMM_N)
#endif
            end if

           pastix_initialised = .true.

          endif

        endif !.not. pastix_initialised
        ! --- End Initialize the PaStiX solver -----------------------------------------------------

        ! --- Analyze the matrix -------------------------------------------------------------------
        if (.not. pastix_analysed) then
          if ( (.not. pastix_smp_only) .or. (pastix_smp_only .and. (my_id_n .eq.0)) ) then

            if (use_pastix) then

              pastix_iparm(IPARM_THREAD_NBR) = pastix_nthrd
              pastix_iparm(IPARM_START_TASK) = API_TASK_ORDERING
              pastix_iparm(IPARM_END_TASK)   = API_TASK_ANALYSE
!              pastix_iparm(IPARM_BINDTHRD)   = API_NO
              !if (my_id_n.eq.0) call timestamp("Reorder",my_id)
#ifdef USE_BLOCK
#ifndef USE_COMPLEX_PRECOND
              call pastix_fortran(pastix_data,MPI_COMM_N, n_block, &
                mumps_par%jcn(1:n_block+1), mumps_par%irn(1:nnz_block), mumps_par%A, &
                pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,Int1,pastix_iparm,pastix_dparm)
#else /* ifndef USE_COMPLEX_PRECOND */
              call pastix_fortran(pastix_data,MPI_COMM_N, n_block, &
                jcn_cmplx(1:n_block+1), irn_cmplx(1:nnz_block), A_cmplx, &
                pastix_perm_vars,pastix_iperm_vars,rhs_cmplx,Int1,pastix_iparm,pastix_dparm)
#endif /* ifndef USE_COMPLEX_PRECOND */
#else /* ifdef USE_BLOCK */
#ifndef USE_COMPLEX_PRECOND
              call pastix_fortran(pastix_data,MPI_COMM_N,mumps_par%n,mumps_par%jcn,mumps_par%irn,mumps_par%A, &
                pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,Int1,pastix_iparm,pastix_dparm)
#else /* ifndef USE_COMPLEX_PRECOND */
              call pastix_fortran(pastix_data,MPI_COMM_N,n_cmplx,jcn_cmplx,irn_cmplx,A_cmplx, &
                pastix_perm_vars,pastix_iperm_vars,rhs_cmplx,Int1,pastix_iparm,pastix_dparm)
#endif /* ifndef USE_COMPLEX_PRECOND */
#endif /* ifdef USE_BLOCK */

            else if (use_wsmp) then
              ! do nothing
            endif

            pastix_analysed = .true.

          endif
        endif ! .not. pastix_analysed
        ! --- End analyze the matrix -------------------------------------------------------------------
        
      endif   ! (else, use_mumps)

      if (my_id_n .eq.0) then                            ! elapsed time analysis end
         call MPI_Barrier(MPI_COMM_MASTER,ierr)
         call clck_time(t1)
         call clck_ldiff(t0,t1,tsecond)
         write(*, FMT_TIMING) my_id,' ## Elapsed time, analysis :',tsecond
         call clck_time(t0)                              ! elapsed time facto start 
      endif

      ! --- Factorize the matrix -------------------------------------------------------------------
      if (use_mumps) then
#ifdef USE_MUMPS

        mumps_par%JOB = 2                                   ! factorisation

        call DMUMPS(mumps_par)

        if (my_id_n .eq.0)   write(*,'(i3,A,i8)')    my_id,' MUMPS, mem       : ',mumps_par%info(16)
#endif
      elseif ( (.not. pastix_smp_only) .or. (pastix_smp_only .and. (my_id_n .eq.0)) ) then

        if (use_pastix) then

          pastix_iparm(IPARM_THREAD_NBR) = pastix_nthrd
          pastix_iparm(IPARM_START_TASK) = API_TASK_NUMFACT
          pastix_iparm(IPARM_END_TASK)   = API_TASK_NUMFACT
          !if (my_id_n.eq.0) call timestamp("Factorize",my_id)
#if defined(WORLDWAR2) && defined(CORES_PER_NODE)
          pastix_iparm(IPARM_BINDTHRD)   = API_BIND_TAB
#endif
          call pastix_bind_threads(my_id)

#ifdef USE_BLOCK

#ifndef USE_COMPLEX_PRECOND
          call pastix_fortran(pastix_data,MPI_COMM_N, n_block, &
            mumps_par%jcn, mumps_par%irn, mumps_par%A, &
            pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,Int1,pastix_iparm,pastix_dparm)
#else /* ifndef USE_COMPLEX_PRECOND */
          call pastix_fortran(pastix_data,MPI_COMM_N, n_block, &
            jcn_cmplx, irn_cmplx, A_cmplx, &
            pastix_perm_vars,pastix_iperm_vars,rhs_cmplx,Int1,pastix_iparm,pastix_dparm)
#endif /* ifndef USE_COMPLEX_PRECOND */

#else /* ifdef USE_BLOCK */

#ifndef USE_COMPLEX_PRECOND
          call pastix_fortran(pastix_data,MPI_COMM_N,mumps_par%n,mumps_par%jcn,mumps_par%irn,mumps_par%A, &
            pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,Int1,pastix_iparm,pastix_dparm)
#else /* ifndef USE_COMPLEX_PRECOND */
          call pastix_fortran(pastix_data,MPI_COMM_N,n_cmplx,jcn_cmplx,irn_cmplx,A_cmplx, &
            pastix_perm_vars,pastix_iperm_vars,rhs_cmplx,Int1,pastix_iparm,pastix_dparm)
#endif /* ifndef USE_COMPLEX_PRECOND */

#endif /* ifdef USE_BLOCK */

        else if (use_wsmp) then
#ifdef USE_WSMP
          call PWGSMP__LU_factorization(my_id_n)
#endif
        end if

      endif
     
      if (my_id_n .eq.0) then                            ! elapsed time facto end
         call MPI_Barrier(MPI_COMM_MASTER,ierr)
         call clck_time(t1)
         call clck_ldiff(t0,t1,tsecond)
         write(*, FMT_TIMING) my_id,' ## Elapsed time, facto :',tsecond
         call clck_time(t0)
      end if
      ! --- End Factorize the matrix ---------------------------------------------------------------
      
   endif NOTSOLVEONLY
   call tr_debug_write("smn_B_mumps_par%n",mumps_par%n)

#ifdef USE_COMPLEX_PRECOND
   if (allocated(rhs_cmplx_guess))  deallocate(rhs_cmplx_guess)
   allocate(rhs_cmplx_guess(1:n_cmplx))
   rhs_cmplx_guess(:) = rhs_cmplx(:) 
#endif

   if (my_id_n .eq. 0) then                          ! elapsed time solve start
      call MPI_Barrier(MPI_COMM_MASTER,ierr)
      call clck_time(t0)
   endif
   
   ! --- Solve the matrix system -------------------------------------------------------------------
   if (use_mumps) then
#ifdef USE_MUMPS
      mumps_par%JOB = 3                                   ! Solve

      call DMUMPS(mumps_par)

#endif
    elseif ( (.not. pastix_smp_only) .or. (pastix_smp_only .and. (my_id_n .eq.0)) ) then

      if (use_pastix) then
        if (.not. pastix_smp_only) then
           call tr_debug_write("smn_C_mumps_par%n",mumps_par%n)
#ifndef USE_COMPLEX_PRECOND
           call MPI_BCAST(mumps_par%rhs,mumps_par%n,MPI_DOUBLE_PRECISION,0,MPI_COMM_N,ierr)
#else
           call MPI_BCAST(rhs_cmplx_guess,n_cmplx,MPI_DOUBLE_COMPLEX,0,MPI_COMM_N,ierr)
#endif
        end if

#ifndef USE_COMPLEX_PRECOND
        if (associated(mumps_par%irn)) call tr_deallocatep(mumps_par%irn,"mumps_par%irn",CAT_DMATRIX)
        if (associated(mumps_par%jcn)) call tr_deallocatep(mumps_par%jcn,"mumps_par%jcn",CAT_DMATRIX)
        if (associated(mumps_par%A))   call tr_deallocatep(mumps_par%A,"mumps_par%A",CAT_DMATRIX)
#else
        if (allocated(A_cmplx))   deallocate(A_cmplx)
        if (allocated(irn_cmplx)) deallocate(irn_cmplx)
        if (allocated(jcn_cmplx)) deallocate(jcn_cmplx)
#endif

#ifndef USE_COMPLEX_PRECOND
        call tr_locvnorms("smn_rhs",mumps_par%rhs,mumps_par%n)
#else
        call tr_locvnorms_cmplx("smn_rhs",rhs_cmplx_guess,n_cmplx)
#endif

        pastix_iparm(IPARM_THREAD_NBR) = pastix_nthrd
        pastix_iparm(IPARM_START_TASK) = API_TASK_SOLVE
        pastix_iparm(IPARM_END_TASK)   = pastix_endsolve
!        pastix_iparm(IPARM_BINDTHRD)   = API_NO
       !if (my_id_n.eq.0) call timestamp("Solve",my_id)
#ifdef USE_BLOCK
#ifndef USE_COMPLEX_PRECOND
        call pastix_fortran(pastix_data,MPI_COMM_N, n_block,                &
!             mumps_par%jcn,mumps_par%irn,mumps_par%A, &
             DUMMY_INT,DUMMY_INT,DUMMY_REAL, &
             pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,Int1,pastix_iparm,pastix_dparm)
#else /* ifndef USE_COMPLEX_PRECOND */
        call pastix_fortran(pastix_data,MPI_COMM_N, n_block,                &
             DUMMY_INT,DUMMY_INT,DUMMY_REAL, &
             pastix_perm_vars,pastix_iperm_vars,rhs_cmplx_guess,Int1,pastix_iparm,pastix_dparm)
#endif /* ifndef USE_COMPLEX_PRECOND */

#else /* ifdef USE_BLOCK */
#ifndef USE_COMPLEX_PRECOND
        call pastix_fortran(pastix_data,MPI_COMM_N,mumps_par%n,&
!             mumps_par%jcn,mumps_par%irn,mumps_par%A, &
          DUMMY_INT,DUMMY_INT,DUMMY_REAL, &
          pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,Int1,pastix_iparm,pastix_dparm)
#else /* ifndef USE_COMPLEX_PRECOND */
        call pastix_fortran(pastix_data,MPI_COMM_N,n_cmplx,&
          DUMMY_INT,DUMMY_INT,DUMMY_REAL, &
          pastix_perm_vars,pastix_iperm_vars,rhs_cmplx_guess,Int1,pastix_iparm,pastix_dparm)
#endif /* ifndef USE_COMPLEX_PRECOND */
#endif /* ifdef USE_BLOCK */

      else if (use_wsmp) then
#ifdef USE_WSMP
        call PWGSMP__back_substitution(mumps_par%rhs, my_id_n)
#endif
     end if
   endif

   if (my_id_n .eq.0) then                            ! elapsed time solve end
     call MPI_Barrier(MPI_COMM_MASTER,ierr)
     call clck_time(t1)
     call clck_ldiff(t0,t1,tsecond)
     write(*, FMT_TIMING) my_id,' ## Elapsed time, solve :',tsecond
     call clck_time(t0)
   end if

#ifdef USE_COMPLEX_PRECOND
   do i = 1, n_cmplx
     if(my_id .eq. 0) then
       mumps_par%rhs(i) = REAL(rhs_cmplx_guess(i)) 
     else
       mumps_par%rhs(2*i-1) = real(rhs_cmplx_guess(i))
       mumps_par%rhs(2*i) = aimag(rhs_cmplx_guess(i))
     endif
   enddo 
#endif /* ifdef USE_COMPLEX_PRECOND */
   ! --- End solve the matrix system ---------------------------------------------------------------

    if (my_id_n .eq. 0) then

      ! --- Undo the column scaling ----------------------------------------------------------------
      do k=1,mumps_par%n
        mumps_par%rhs(k) =  mumps_par%rhs(k) / column_scaling(k)
      enddo

      if (allocated(deltas)) call tr_deallocate(deltas,"deltas",CAT_PRECOND)
      call tr_allocate(deltas,Int1,ndof_glob,"deltas",CAT_PRECOND)
      deltas = 0.d0

      call tr_allocate(rhs_tmp,Int1,ndof_glob,"rhs_tmp",CAT_PRECOND)

      rhs_tmp = 0.d0

      do i = 1, mumps_par%n
        rhs_tmp(my_row_index(i)) = mumps_par%rhs(i)*my_row_factor
      enddo

      ! --- End undo the column scaling ------------------------------------------------------------
      
      ! --- Collect the RHSs from all harmonic matrices --------------------------------------------
      call MPI_AllReduce(RHS_tmp,deltas,ndof_glob,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_MASTER,ierr)
      call tr_deallocate(rhs_tmp,"rhs_tmp",CAT_PRECOND)
      ! --- End collect the RHSs from all harmonic matrices ----------------------------------------

      call tr_locvnorms("smn_res",mumps_par%rhs,mumps_par%n)
      call tr_locvnorms("smn_delta",deltas,ndof_glob)
    endif

    ! -- For PaStiX solver before version 6.x
    call tr_set_precondmem(pastix_dparm(2)) ! DPARM_MEM_MAX DEPRECATED IN PASTIX6: how to change this?
    ! ############### This should be looked at at some point

    call tr_print_memsize("AfterSolveN")
    call r3_info_end (r3_info_index_0)         ! timing
    return
 
  end subroutine solve_matrix_n
#endif /* defined(USE_PASTIX) */

! > Solve the harmonic matrix system using STRUMPACK
#ifdef USE_STRUMPACK  
subroutine solve_matrix_n_spk(my_id,MPI_COMM_N,MPI_COMM_MASTER,solve_only)
    use tr_module
    use iso_c_binding
    use mod_parameters
    use mumps_module
    use global_distributed_matrix
    use mpi_mod 
    use mod_clock
    use phys_module, only : index_now, centralize_harm_mat

    use strumpack_module
    use matio_module, only :  save_mat_h5
    use mod_integer_types
  
    implicit none

#include "r3_info.h"

    integer, intent(in) :: my_id
    integer, intent(in) :: MPI_COMM_N, MPI_COMM_MASTER
    logical, intent(in) :: solve_only

    integer               :: my_id_n, n_cpu_n, ierr, block_size
    integer(kind=int_all) :: i, j, k
    type(clcktype)        :: t_itstart, t0, t1, t2, t3
    real*8                :: tsecond
    real*8, allocatable   :: RHS_tmp(:)

    !Split broadcast
    character*8 :: type

    integer(kind=int_all), parameter   :: Int1=1
    
    integer(kind=C_INT_ALL) :: n, nnz

    call r3_info_begin (r3_info_index_0, 'solve_matrix_n')                  ! timing
    call tr_print_memsize("BeforeSolveN")
    call tr_debug_write("smn_A_mumps_par%n",mumps_par%n)

    if (my_id .eq. 0) then
      write(*,*) my_id,'*********************************'
      write(*,*) my_id,'*      solve local matrix  (n)  *'
      write(*,*) my_id,'*********************************'
      write(*,*) my_id,'*     using solver STRUMPACK    *'
      write(*,*) my_id,'*********************************'
    endif

    call MPI_COMM_RANK(MPI_COMM_N, my_id_n, ierr)     ! the id of each cpu
    call MPI_COMM_SIZE(MPI_COMM_N, n_cpu_n, ierr)     ! the number of cpus

    if (centralize_harm_mat) then 
      call MPI_BCAST(mumps_par%n,1,MPI_INTEGER_ALL,0,MPI_COMM_N,ierr)
      call MPI_BCAST(mumps_par%nz,1,MPI_INTEGER_ALL,0,MPI_COMM_N,ierr)
    endif
    
    n = mumps_par%n
    nnz = mumps_par%nz
    block_size = n_var*my_mode_set_n
    
    if (.not. solve_only) then
      
      if (.not. spss_initialized) then
        call strumpack_init(MPI_COMM_N)
        spss_initialized = .true.
      endif     

      if (centralize_harm_mat) then
        ! broadcast centralized matrix
        if (my_id_n.gt.0) then
          if (associated(mumps_par%irn)) call tr_deallocatep(mumps_par%irn,"mumps_par%irn",CAT_DMATRIX)
          if (associated(mumps_par%jcn)) call tr_deallocatep(mumps_par%jcn,"mumps_par%jcn",CAT_DMATRIX)
          if (associated(mumps_par%A))   call tr_deallocatep(mumps_par%a,"mumps_par%A",CAT_DMATRIX)
          call tr_allocatep(mumps_par%irn,Int1,nnz,"mumps_par%irn",CAT_DMATRIX)
          call tr_allocatep(mumps_par%jcn,Int1,nnz,"mumps_par%jcn",CAT_DMATRIX)
          call tr_allocatep(mumps_par%a,Int1,nnz,"mumps_par%a",CAT_DMATRIX)
        endif  
  
        ! Split MPI_BCAST if MPI buffer beyond 2Go
        type='intIRN'
        call split_broadcast(type,MPI_COMM_N)
        type='intJCN'
        call split_broadcast(type,MPI_COMM_N)
        type='double'
        call split_broadcast(type,MPI_COMM_N)

      endif
      
      !if (my_id_n.eq.0) call save_mat_h5(my_id,n,nnz,mumps_par%irn,mumps_par%jcn,mumps_par%a,mumps_par%rhs)      
      !if (my_id_n.eq.0) call timestamp("Set mat",my_id)
      call strumpack_set_mat(mumps_par%n,mumps_par%nz,mumps_par%irn,mumps_par%jcn,mumps_par%a,block_size,&
                MPI_COMM_N,UPDATE=spss_analyzed,DISTRIBUTED=.not.centralize_harm_mat,EQUILIBRIUM=.false.)      

      if (associated(mumps_par%irn)) call tr_deallocatep(mumps_par%irn,"mumps_par%irn",CAT_DMATRIX)
      if (associated(mumps_par%jcn)) call tr_deallocatep(mumps_par%jcn,"mumps_par%jcn",CAT_DMATRIX)
      if (associated(mumps_par%A))   call tr_deallocatep(mumps_par%A,"mumps_par%A",CAT_DMATRIX)

      if (.not. spss_analyzed) then
        if (my_id_n.eq. 0) then                  ! elapsed time reorder start
          call MPI_Barrier(MPI_COMM_MASTER,ierr)
          call clck_time(t0)
        endif
              
        !if (my_id_n.eq.0) call timestamp("Reorder",my_id)
        call strumpack_analyze(MPI_COMM_N)
        spss_analyzed = .true.
        if (my_id_n .eq.0) then                  ! elapsed time reorder end
          call MPI_Barrier(MPI_COMM_MASTER,ierr)
          call clck_time(t1)
          call clck_ldiff(t0,t1,tsecond)
          write(*, FMT_TIMING) my_id,' ## Elapsed time, analysis :',tsecond
        endif        
      endif
      
      if (my_id_n.eq. 0) then                   ! elapsed time factorization start
      call MPI_Barrier(MPI_COMM_MASTER,ierr)
          call clck_time(t0)
        endif      
      !if (my_id_n.eq.0) call timestamp("Factorize",my_id)
      call strumpack_factorize(MPI_COMM_N)
      
      if (my_id_n.eq.0) then                   ! elapsed time facto end
        call MPI_Barrier(MPI_COMM_MASTER,ierr)
        call clck_time(t1)
        call clck_ldiff(t0,t1,tsecond)
        write(*, FMT_TIMING) my_id,' ## Elapsed time, factorization :',tsecond
      endif       

    endif ! .not. solve_only
    
    if (my_id_n.gt.0) then
      if (associated(mumps_par%rhs)) call tr_deallocatep(mumps_par%rhs,"mumps_par%rhs",CAT_DMATRIX)
      call tr_allocatep(mumps_par%rhs,Int1,n,"mumps_par%rhs",CAT_DMATRIX)
    endif
    
    call MPI_BCAST(mumps_par%rhs,n,MPI_DOUBLE_PRECISION,0,MPI_COMM_N,ierr)
    
    if (my_id_n .eq. 0) then                          ! elapsed time solve start
      call MPI_Barrier(MPI_COMM_MASTER,ierr)
      call clck_time(t0)
    endif    
    
    call MPI_Barrier(MPI_COMM_N,ierr)
    !if (my_id_n.eq.0) call timestamp("Solve",my_id)
    call strumpack_solve(n,mumps_par%rhs,MPI_COMM_N)
    
    if (my_id_n .eq.0) then                            ! elapsed time solve end
       call MPI_Barrier(MPI_COMM_MASTER,ierr)
       call clck_time(t1)
       call clck_ldiff(t0,t1,tsecond)
       write(*, FMT_TIMING) my_id,' ## Elapsed time, solve :',tsecond
       call clck_time(t0)
    endif    

    if (my_id_n .eq. 0) then

      if (allocated(deltas)) call tr_deallocate(deltas,"deltas",CAT_PRECOND)
      call tr_allocate(deltas,Int1,ndof_glob,"deltas",CAT_PRECOND)
      deltas = 0.d0

      call tr_allocate(rhs_tmp,Int1,ndof_glob,"rhs_tmp",CAT_PRECOND)

      rhs_tmp = 0.d0
      do i = 1, mumps_par%n
        rhs_tmp(my_row_index(i)) = mumps_par%rhs(i)*my_row_factor
      enddo


      call MPI_AllReduce(RHS_tmp,deltas,ndof_glob,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_MASTER,ierr)
      call tr_deallocate(rhs_tmp,"rhs_tmp",CAT_PRECOND)

      call tr_locvnorms("smn_res",mumps_par%rhs,mumps_par%n)
      call tr_locvnorms("smn_delta",deltas,ndof_glob)
    endif  
    
    call tr_print_memsize("AfterSolveN")
    call r3_info_end (r3_info_index_0)         ! timing

    return

  end subroutine solve_matrix_n_spk
#endif

#if defined(USE_PASTIX6)
! > Solve the preconditioner system using PaStiX 6.2
subroutine solve_matrix_n_ptx(my_id,MPI_COMM_N,MPI_COMM_MASTER,solve_only)
    use tr_module
    use iso_c_binding
    use mod_parameters
    use mumps_module
    use global_distributed_matrix, only: ndof_glob, column_scaling, deltas
    use mpi_mod 
    use mod_clock
    use phys_module, only : index_now, centralize_harm_mat

    use mod_pastix, only: spm_initialized, spm_analyzed, pastix_init, pastix_set_mat, pastix_analyze, &
            pastix_factorize, pastix_solve, pastix_finalize
    use matio_module, only :  save_mat_h5, timestamp
    use mod_integer_types
  
    implicit none

#include "r3_info.h"

    integer, intent(in) :: my_id
    integer, intent(in) :: MPI_COMM_N, MPI_COMM_MASTER
    logical, intent(in) :: solve_only

    integer               :: my_id_n, n_cpu_n, ierr, block_size, indexing=1
    integer(kind=int_all) :: i, j, k
    type(clcktype)        :: t_itstart, t0, t1, t2, t3
    real*8                :: tsecond
    real*8, allocatable   :: RHS_tmp(:)

    !Split broadcast
    character*8 :: type

    integer(kind=int_all), parameter   :: Int1=1
    
    integer(kind=C_INT_ALL) :: n, nnz

    call r3_info_begin (r3_info_index_0, 'solve_matrix_n')                  ! timing
    call tr_print_memsize("BeforeSolveN")
    call tr_debug_write("smn_A_mumps_par%n",mumps_par%n)

    if (my_id .eq. 0) then
      write(*,*) my_id,'*********************************'
      write(*,*) my_id,'*      solve local matrix  (n)  *'
      write(*,*) my_id,'*********************************'
      write(*,*) my_id,'*     using solver PaStiX v6.2  *'
      write(*,*) my_id,'*********************************'
    endif

    call MPI_COMM_RANK(MPI_COMM_N, my_id_n, ierr)     ! the id of each cpu
    call MPI_COMM_SIZE(MPI_COMM_N, n_cpu_n, ierr)     ! the number of cpus

    if (centralize_harm_mat) then 
      call MPI_BCAST(mumps_par%n,1,MPI_INTEGER_ALL,0,MPI_COMM_N,ierr)
      call MPI_BCAST(mumps_par%nz,1,MPI_INTEGER_ALL,0,MPI_COMM_N,ierr)
    endif
    
    n = mumps_par%n
    nnz = mumps_par%nz
    block_size = n_var*my_mode_set_n
    
    if (.not. solve_only) then
      
      if (.not. spm_initialized) then
        call pastix_init(MPI_COMM_N)
        spm_initialized = .true.
      endif

      if (centralize_harm_mat) then
        ! broadcast centralized matrix
        if (my_id_n.gt.0) then
          if (associated(mumps_par%irn)) call tr_deallocatep(mumps_par%irn,"mumps_par%irn",CAT_DMATRIX)
          if (associated(mumps_par%jcn)) call tr_deallocatep(mumps_par%jcn,"mumps_par%jcn",CAT_DMATRIX)
          if (associated(mumps_par%A))   call tr_deallocatep(mumps_par%a,"mumps_par%A",CAT_DMATRIX)
          call tr_allocatep(mumps_par%irn,Int1,nnz,"mumps_par%irn",CAT_DMATRIX)
          call tr_allocatep(mumps_par%jcn,Int1,nnz,"mumps_par%jcn",CAT_DMATRIX)
          call tr_allocatep(mumps_par%a,Int1,nnz,"mumps_par%a",CAT_DMATRIX)
        endif  
  
        ! Split MPI_BCAST if MPI buffer beyond 2Go
        type='intIRN'
        call split_broadcast(type,MPI_COMM_N)
        type='intJCN'
        call split_broadcast(type,MPI_COMM_N)
        type='double'
        call split_broadcast(type,MPI_COMM_N)
      endif
      
      !if (my_id_n.eq.0) call save_mat_h5(my_id,n,nnz,mumps_par%irn,mumps_par%jcn,mumps_par%a,mumps_par%rhs)
      !if (my_id_n.eq.0) call timestamp("Set mat",my_id)
      call pastix_set_mat(n, nnz, mumps_par%irn, mumps_par%jcn, mumps_par%a, block_size, MPI_COMM_N,&
                UPDATE=spm_analyzed, DISTRIBUTED=.not.centralize_harm_mat, EQUILIBRIUM=.false.)
                
      if (.not. spm_analyzed) then
        if (my_id_n.eq. 0) then
          call MPI_Barrier(MPI_COMM_MASTER,ierr)
          call clck_time(t0)
        endif
              
        !if (my_id_n.eq.0) call timestamp("Reorder",my_id)
        call pastix_analyze()
        spm_analyzed = .true.
        if (my_id_n .eq.0) then                  ! elapsed time reorder end
          call MPI_Barrier(MPI_COMM_MASTER,ierr)
          call clck_time(t1)
          call clck_ldiff(t0,t1,tsecond)
          write(*, FMT_TIMING) my_id,' ## Elapsed time, analysis :',tsecond
        endif        
      endif
      
      if (my_id_n.eq. 0) then                   ! elapsed time factorization start
      call MPI_Barrier(MPI_COMM_MASTER,ierr)
          call clck_time(t0)
        endif      
      !if (my_id_n.eq.0) call timestamp("Factorize",my_id)
      call pastix_factorize()
      
      if (my_id_n.eq.0) then                   ! elapsed time facto end
        call MPI_Barrier(MPI_COMM_MASTER,ierr)
        call clck_time(t1)
        call clck_ldiff(t0,t1,tsecond)
        write(*, FMT_TIMING) my_id,' ## Elapsed time, factorization :',tsecond
      endif       

    endif ! .not. solve_only
    
    if (my_id_n.gt.0) then
      if (associated(mumps_par%rhs)) call tr_deallocatep(mumps_par%rhs,"mumps_par%rhs",CAT_DMATRIX)
      call tr_allocatep(mumps_par%rhs,Int1,n,"mumps_par%rhs",CAT_DMATRIX)
    endif
    
    call MPI_BCAST(mumps_par%rhs,n,MPI_DOUBLE_PRECISION,0,MPI_COMM_N,ierr)
    
    if (my_id_n .eq. 0) then                          ! elapsed time solve start
      call MPI_Barrier(MPI_COMM_MASTER,ierr)
      call clck_time(t0)
    endif    
    
    call MPI_Barrier(MPI_COMM_N,ierr)
    !if (my_id_n.eq.0) call timestamp("Solve",my_id)
    call pastix_solve(mumps_par%n, mumps_par%rhs)
    
    if (my_id_n .eq.0) then                            ! elapsed time solve end
       call MPI_Barrier(MPI_COMM_MASTER,ierr)
       call clck_time(t1)
       call clck_ldiff(t0,t1,tsecond)
       write(*, FMT_TIMING) my_id,' ## Elapsed time, solve :',tsecond
       call clck_time(t0)
    endif    

    if (my_id_n .eq. 0) then

      if (allocated(deltas)) call tr_deallocate(deltas,"deltas",CAT_PRECOND)
      call tr_allocate(deltas,Int1,ndof_glob,"deltas",CAT_PRECOND)
      deltas = 0.d0

      call tr_allocate(rhs_tmp,Int1,ndof_glob,"rhs_tmp",CAT_PRECOND)
     
      rhs_tmp = 0.d0
      do i = 1, mumps_par%n
        rhs_tmp(my_row_index(i)) = mumps_par%rhs(i)*my_row_factor
      enddo
      write(*,*) my_id, "solution", mumps_par%rhs(1), mumps_par%rhs(mumps_par%n)

      call MPI_AllReduce(rhs_tmp,deltas,ndof_glob,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_MASTER,ierr)
      call tr_deallocate(rhs_tmp,"rhs_tmp",CAT_PRECOND)

      call tr_locvnorms("smn_res",mumps_par%rhs,mumps_par%n)
      call tr_locvnorms("smn_delta",deltas,ndof_glob)
    endif  
    
    call tr_print_memsize("AfterSolveN")
    call r3_info_end (r3_info_index_0)         ! timing

    return

  end subroutine solve_matrix_n_ptx

#endif

end module solve_mat_n
