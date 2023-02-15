#ifdef USE_MUMPS
subroutine solve_mumps_all(mmss, a_mat, rhs_vec, solve_only, tag)
  use tr_module
  use mpi_mod
  use mod_clock
  use data_structure, only: type_SP_MATRIX, type_RHS
  use mod_mumps,      only: type_MUMPS_SOLVER, mumps_initialize, mumps_analyze, mumps_factorize, mumps_solve

  implicit none

  type(type_SP_MATRIX)     :: a_mat
  type(type_RHS)           :: rhs_vec
  type(type_MUMPS_SOLVER)  :: mmss
  logical                  :: solve_only
  integer                  :: tag

  real*8,allocatable       :: column_local(:)
  real*8                   :: tsecond, t_analysis_0, t_analysis_1, t_fact_0, t_fact_1
  type(clcktype)           :: t0, t1
  integer                  :: k, j
  integer                  :: my_id, n_cpu, comm, ierr
  integer(kind=int_all)    :: i

  logical                  :: verbose = .false.

  comm = a_mat%comm

  call MPI_COMM_RANK(comm, my_id, ierr)
  call MPI_COMM_SIZE(comm, n_cpu, ierr)

  if ((tag.ge.0).and.(my_id.eq.0)) verbose = .true.

  if (.not.solve_only) then

    call scale_by_cols(a_mat)

    if (associated(mmss%solution_scaling)) then
      deallocate(mmss%solution_scaling); mmss%solution_scaling => Null()
    endif
    allocate(mmss%solution_scaling(a_mat%ng))
    do i = 1, a_mat%ng
      mmss%solution_scaling(i) = a_mat%column_scaling(i)
    enddo
    mmss%scaled = .true.

    if (.not. mmss%initialized) then
      call mumps_initialize(mmss,comm)
    endif

    if (.not. mmss%analyzed) then
      call clck_time(t0)

      call mumps_analyze(mmss,a_mat)

      call clck_time(t1); call clck_ldiff(t0,t1,tsecond)
      if (verbose) write(*,FMT_TIMING) tag,  '## Elapsed time analysis:', tsecond

    endif

    call clck_time(t0)

    call mumps_factorize(mmss,a_mat)

    call clck_time(t1); call clck_ldiff(t0,t1,tsecond)
    if (verbose) write(*,FMT_TIMING) tag,  '## Elapsed time factorize:', tsecond

  endif

  call clck_time(t0)

  call mumps_solve(mmss,rhs_vec)

  call clck_time(t1); call clck_ldiff(t0,t1,tsecond)
  if (verbose) write(*,FMT_TIMING) tag,  '## Elapsed time solve:', tsecond

  return

end
#endif