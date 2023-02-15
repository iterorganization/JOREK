#if defined(USE_PASTIX)||defined(USE_PASTIX6)
!> subroutine solves the complete system of equation using pastix with
!  distributed matrix ad_mat on the main group mpi_comm_world.
!  For pastix5 solver matrix is centralized into ac_mat
subroutine solve_pastix_all(ptss, ad_mat, rhs_vec, solve_only, tag)
  use tr_module
  use mod_parameters, only: n_tor, n_var
  use mpi_mod
  use mod_clock

  use mod_integer_types
  use data_structure, only: type_SP_MATRIX, type_RHS
  use mod_pastix, only:     type_PASTIX_SOLVER, pastix_set_mat, pastix_solve, pastix_factorize, pastix_analyze, pastix_initialize

  !use matio_module, only: save_mat_h5

  implicit none

! --- Input variables
  type(type_SP_MATRIX)              :: ad_mat
  type(type_RHS)                    :: rhs_vec
  type(type_PASTIX_SOLVER)          :: ptss
  logical                           :: solve_only
  integer                           :: tag

! --- Local variables
  type(clcktype)                    :: t_itstart, t0, t1, t2, t3
  real*8                            :: tsecond
  integer                           :: n_cpu, my_id, ierr, comm
  type(type_SP_MATRIX)              :: ac_mat
  logical                           :: verbose = .false.

  comm = ad_mat%comm

  call MPI_COMM_RANK(comm, my_id, ierr)
  call MPI_COMM_SIZE(comm, n_cpu, ierr)
  
  if ((tag.ge.0).and.(my_id.eq.0)) verbose = .true.

  if (.not.solve_only) then

    ptss%rhs_val => rhs_vec%val

#ifdef USE_PASTIX
    call pastix_set_mat(ptss, ad_mat, ac_mat, tag)

    if (.not. ptss%initialized) then
      ptss%comm = ad_mat%comm
      call pastix_initialize(ptss)
    endif
#elif USE_PASTIX6
    if (.not. ptss%initialized) then
      ptss%comm = ad_mat%comm
      call pastix_initialize(ptss)
    endif

    call pastix_set_mat(ptss, ad_mat, ac_mat, tag)
#endif

    !call save_mat_h5(tag,ac_mat%ng,ac_mat%nnz,ac_mat%irn,ac_mat%jcn,ac_mat%val,rhs_vec%val)

    if (.not. ptss%analyzed) then
      if (verbose) write(*,*) "PaStiX: analyzing matrix"
      call clck_time(t0)

      call pastix_analyze(ptss)

      call clck_time(t1); call clck_ldiff(t0,t1,tsecond)
      if (verbose)  write(*,FMT_TIMING) tag, '## Elapsed time analysis:', tsecond
    endif

    if (verbose) write(*,*) "PaStiX: factorizing matrix"
    call clck_time(t0)

    call pastix_factorize(ptss)

    call clck_time(t1); call clck_ldiff(t0,t1,tsecond)
    if (verbose)  write(*,FMT_TIMING) tag, '## Elapsed time factorization:', tsecond

! matrix isn't needed after factorization
    if (associated(ac_mat%irn)) call tr_deallocatep(ac_mat%irn,"irn",CAT_DMATRIX)
    if (associated(ac_mat%jcn)) call tr_deallocatep(ac_mat%jcn,"jcn",CAT_DMATRIX)
    if (associated(ac_mat%val)) call tr_deallocatep(ac_mat%val,"val",CAT_DMATRIX)
    if (associated(ad_mat%irn)) call tr_deallocatep(ad_mat%irn,"irn",CAT_DMATRIX)
    if (associated(ad_mat%jcn)) call tr_deallocatep(ad_mat%jcn,"jcn",CAT_DMATRIX)
    if (associated(ad_mat%val)) call tr_deallocatep(ad_mat%val,"val",CAT_DMATRIX)

  endif

  call clck_time(t0)

  call pastix_solve(ptss,rhs_vec)

  call clck_time(t1); call clck_ldiff(t0,t1,tsecond)
  if (verbose)  write(*,FMT_TIMING) tag, '## Elapsed time solve:', tsecond

  return
end
#endif
