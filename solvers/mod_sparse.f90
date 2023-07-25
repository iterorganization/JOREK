module mod_sparse
  use mod_sparse_data

  private
  public :: solve_sparse_system

  contains

!> solve Ax = rhs
!! sol_vec contains the initial guess
!! sol_vec, rhs_vec are broadcasted
!! solve_type - type of system, e.g. GS equilibrium, MHD system with preconditioner, etc.
  subroutine solve_sparse_system(a_mat, rhs_vec, sol_vec, solver, mhd_sim)
    use mod_integer_types
    use mod_clock
    use data_structure, only: type_SP_MATRIX, type_PRECOND, type_RHS
    use mod_simulation_data, only: type_MHD_SIM
    use mod_sparse_data, only: type_SP_SOLVER, mumps, pastix, strumpack
    use mod_preconditioner, only: initialize_preconditioner, reset_preconditioner, update_pc_rhs, gather_solution
#ifdef DIRECT_CONSTRUCTION
    use mod_direct_construction, only: update_pc_mat
#else
    use mod_distribute_preconditioner, only: update_pc_mat
#endif
#ifdef USE_STRUMPACK
    use mod_strumpack, only: spk_delete_factors
#endif
#ifdef USE_BICGSTAB
    use mod_bicgstab, only: bicgstab_driver
#else
    use mod_gmres, only: gmres_driver
#endif
    use matio_module, only: save_mat_h5

    implicit none

    type(type_SP_SOLVER)          :: solver
    type(type_SP_MATRIX)          :: a_mat
    type(type_RHS)                :: rhs_vec, sol_vec
    type(type_MHD_SIM), optional  :: mhd_sim

    integer                  :: my_id, n_cpu, ierr
    type(clcktype)           :: t_itstart, t0, t1, t2, t3
    real*8                   :: tsecond
    integer(kind=int_all)    :: i
    logical                  :: verbose = .false.
    integer                  :: tag = -1   !< tag for log file output
    character(len=10)        :: fname


    external :: solve_mumps_all, solve_pastix_all, solve_strumpack_all

    call MPI_COMM_SIZE(a_mat%comm, n_cpu, ierr)
    call MPI_COMM_RANK(a_mat%comm, my_id, ierr)

    sol_vec%n = rhs_vec%n

    verbose = solver%verbose.and.(my_id.eq.0)
#ifdef SAVEMATRIX
    write(fname,'(A5,I2.2,A3)') "matA_",my_id,".h5"
    call save_mat_h5_ext(fname, a_mat%ng, a_mat%ng, a_mat%nnz, &
                           a_mat%irn, a_mat%jcn, a_mat%val, rhs=rhs_vec%val, &
                           ind_min=a_mat%index_min(my_id+1),ind_max=a_mat%index_max(my_id+1), &
                           block_size=a_mat%block_size)
#endif

    if (.not.solver%iterative) then

      if (verbose) tag = 0

      if (solver%equilibrium) then
        if (verbose) write(*,*) "Solving MHD equilibrium system"
      else
        if (verbose) write(*,*) "Solving MHD system using direct solver"
      endif

      if (solver%library.eq.mumps) then
#ifdef USE_MUMPS
        if (verbose) write(*,*) "Using MUMPS solver"
        solver%mmss%equilibrium = solver%equilibrium
        call solve_mumps_all(solver%mmss, a_mat, rhs_vec, solver%solve_only, tag)
#endif
      elseif (solver%library.eq.strumpack) then
#ifdef USE_STRUMPACK
        if (verbose) write(*,*) "Using STRUMPACK solver"
        solver%spss%equilibrium = solver%equilibrium
        call solve_strumpack_all(solver%spss, a_mat, rhs_vec, solver%solve_only, tag)
#endif
      elseif (solver%library.eq.pastix) then
#if (defined USE_PASTIX) || (defined USE_PASTIX6)
        if (verbose) write(*,*) "Using PaStiX solver"
        solver%ptss%equilibrium = solver%equilibrium
        solver%ptss%refine = .true.
        call solve_pastix_all(solver%ptss, a_mat, rhs_vec, solver%solve_only, tag)
#endif
      endif

      do i=1,rhs_vec%n
        sol_vec%val(i) =  rhs_vec%val(i)
      enddo

      solver%step_success = .true.

    elseif (solver%iterative) then

      if (verbose) write(*,*) "Solving MHD system using iterative solver"

      if (solver%verbose) tag = my_id

      ! condition for no PC update
      solver%solve_only = (solver%istep.gt.1) .and. ((solver%iter_gmres + solver%iter_prev <= 2*solver%iter_precon) &
                                              .and.  (solver%n_since_update < solver%max_steps_noUpdate))
      if (solver%newton%it.lt.2) then ! no counter within Newton loop
        if (solver%solve_only) then
          solver%n_since_update = solver%n_since_update + 1
        else
          solver%n_since_update = 0
        endif
      endif
      solver%solve_only = (solver%solve_only).or.(solver%newton%it.gt.1) ! no PC update within Newton loop

      if (.not.solver%pc%initialized) then
        call initialize_preconditioner(solver%pc,a_mat%comm)
        ! set whether to distribute pc matrix when constructing by communication
        if (solver%library.eq.strumpack) solver%pc%mat%row_distributed = .true.
#if (defined USE_PASTIX6)
        if (solver%library.eq.pastix) solver%pc%mat%col_distributed = .true.
#endif
      endif

! Finding PC solution
      if (.not.solver%solve_only) then
#ifdef USE_STRUMPACK
        if ((solver%library.eq.strumpack).and.(solver%spss%analyzed)) call spk_delete_factors(solver%spss%sscp)
#endif
        call update_pc_mat(solver%pc,a_mat,mhd_sim)
      endif

      call update_pc_rhs(solver%pc,rhs_vec)
#ifdef SAVEMATRIX
      if (.not.solver%solve_only) then
        write(fname,'(A3,I2.2,A3)') "pc_",my_id,".h5"
        call save_mat_h5_ext(fname, a_mat%ng, solver%pc%mat%ng,solver%pc%mat%nnz, &
                                       solver%pc%mat%irn, solver%pc%mat%jcn,solver%pc%mat%val, &
                                       l2g=solver%pc%row_index,rhs=solver%pc%rhs%val, block_size=solver%pc%mat%block_size)
      endif
#endif
      if (solver%library.eq.mumps) then
#ifdef USE_MUMPS
        call solve_mumps_all(solver%mmss, solver%pc%mat, solver%pc%rhs, solver%solve_only, tag)
#endif
      elseif (solver%library.eq.strumpack) then
#ifdef USE_STRUMPACK
        call solve_strumpack_all(solver%spss, solver%pc%mat, solver%pc%rhs, solver%solve_only, tag)
#endif
      elseif (solver%library.eq.pastix) then
#if (defined USE_PASTIX) || (defined USE_PASTIX6)
        call solve_pastix_all(solver%ptss, solver%pc%mat, solver%pc%rhs, solver%solve_only, tag)
#endif
      endif

      call MPI_Barrier(a_mat%comm, ierr)

      call gather_solution(solver%pc,sol_vec)

! iterative part
      solver%iter_prev  = solver%iter_gmres
      solver%iter_gmres = solver%iter_max

#ifdef USE_BICGSTAB
      call bicgstab_driver(a_mat, rhs_vec, sol_vec, solver)
#else
      call gmres_driver(a_mat, rhs_vec, sol_vec, solver)
#endif

      if (verbose) write(*,'(A32,I5)') 'Number of iterations: ', solver%iter_gmres

      solver%step_success = (solver%iter_gmres .lt. solver%iter_max)

    endif

  end subroutine solve_sparse_system


end module mod_sparse
