#ifdef USE_STRUMPACK
!> subroutine solves the complete system of equation using STRUMPACK
! takes distributed matrix ad_mat, centralize it and solve, placing the solution into the rhs_vec
subroutine solve_strumpack_all(spss, ad_mat, rhs_vec, solve_only, tag)
  use mod_strumpack

  use tr_module
  use mpi_mod
  use mod_clock
  use mod_integer_types
  use data_structure, only: type_SP_MATRIX, type_RHS

  implicit none

! --- Input variables
  type(type_SP_MATRIX)        :: ad_mat
  type(type_RHS)              :: rhs_vec
  type(type_STRUMPACK_SOLVER) :: spss
  logical                     :: solve_only
  integer                     :: tag

! --- Local variables
  type(type_SP_MATRIX)        :: ac_mat
  type(clcktype)              :: t_itstart, t0, t1, t2, t3
  real*8                      :: tsecond
  integer                     :: my_id, n_cpu, comm, ierr
  logical                     :: centralize = .true.
  logical                     :: verbose = .false.

  external matrix_split_reduce

  comm = ad_mat%comm

  call MPI_COMM_RANK(comm, my_id, ierr)
  call MPI_COMM_SIZE(comm, n_cpu, ierr)

  if ((tag.ge.0).and.(my_id.eq.0)) verbose = .true.

  if (.not.solve_only) then

    centralize = (n_cpu.gt.1).and.(.not.ad_mat%row_distributed)

    if (centralize) then

      call clck_time(t0)

      call matrix_split_reduce(ad_mat,ac_mat)

      call clck_time(t1); call clck_ldiff(t0,t1,tsecond)
      if (verbose)  write(*,FMT_TIMING) tag, '## Elapsed time mpi_gather :', tsecond

    else
      call ad_mat%move_to(ac_mat, with_data=.true.)
    endif

    if (.not. spss%initialized) then
      call strumpack_init(spss, comm)
      spss%initialized = .true.
    endif

    call strumpack_set_mat(spss, ac_mat)

    ! strumpack copies matrix values into permuted sparse matrix, thus the original matrix can be deallocated
    if (associated(ac_mat%irn)) call tr_deallocatep(ac_mat%irn,"irn",CAT_DMATRIX)
    if (associated(ac_mat%jcn)) call tr_deallocatep(ac_mat%jcn,"jcn",CAT_DMATRIX)
    if (associated(ac_mat%val)) call tr_deallocatep(ac_mat%val,"val",CAT_DMATRIX)
    if (associated(ad_mat%irn)) call tr_deallocatep(ad_mat%irn,"irn",CAT_DMATRIX)
    if (associated(ad_mat%jcn)) call tr_deallocatep(ad_mat%jcn,"jcn",CAT_DMATRIX)
    if (associated(ad_mat%val)) call tr_deallocatep(ad_mat%val,"val",CAT_DMATRIX)

    if (.not. spss%analyzed) then
      call clck_time(t0)

      call strumpack_analyze(spss)

      call clck_time(t1); call clck_ldiff(t0,t1,tsecond)
      if (verbose)  write(*,FMT_TIMING) tag, '## Elapsed time analysis:', tsecond
    endif

    call clck_time(t0)

    call strumpack_factorize(spss)

    call clck_time(t1); call clck_ldiff(t0,t1,tsecond)
    if (verbose)  write(*,FMT_TIMING) tag, '## Elapsed time factorize:', tsecond

  endif ! .not.solve_only

  call clck_time(t0)

  call strumpack_solve(spss, rhs_vec)

  call clck_time(t1); call clck_ldiff(t0,t1,tsecond)
  if (verbose)  write(*,FMT_TIMING) tag, '## Elapsed time solve:', tsecond



  return
end subroutine solve_strumpack_all
#endif
