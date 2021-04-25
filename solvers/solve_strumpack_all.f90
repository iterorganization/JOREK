#ifdef USE_STRUMPACK      
!> subroutine solves the complete system of equation using STRUMPACK
subroutine solve_strumpack_all(n_cpu,my_id,index_min,index_max)
  use strumpack_module

  use tr_module 
  use mod_parameters
  use mumps_module
  use global_distributed_matrix
  use mpi_mod
  use mod_clock
  use mod_integer_types

!$ use omp_lib

  implicit none

! --- Routine parameters
  integer,               intent(in) :: n_cpu, my_id
  integer,               intent(in) :: index_min, index_max

! --- Local variables
  type(clcktype)                    :: t_itstart, t0, t1, t2, t3
  real*8                            :: tsecond
  integer                           :: i, k, j, ierr
  integer(kind=int_all)             :: m_loc
  integer(kind=int_all),allocatable :: counts(:), displacements(:)
  integer(kind=int_all), parameter  :: Int1=1
  
  integer(kind=C_INT_ALL) :: n, nnz

!write(*,*) my_id,'*********************************'
!write(*,*) my_id,'*  solve global matrix using STRUMPACK *'
!write(*,*) my_id,'*********************************'

  m_loc = (index_max - index_min + 1) * n_tor * n_var
  mumps_par%nz_loc = nz_glob

  call MPI_Allreduce(m_loc,mumps_par%n,1,MPI_INTEGER_ALL,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_Allreduce(mumps_par%nz_loc,mumps_par%nz,1,MPI_INTEGER_ALL,MPI_SUM,MPI_COMM_WORLD,ierr)
  
  n = mumps_par%n
  nnz = mumps_par%nz

  
  call clck_time(t0)

!------------------------------------------------------ collect the distributed matrix onto all procs
  if (allocated(counts))        call tr_deallocate(counts,"counts",CAT_DMATRIX)
  if (allocated(displacements)) call tr_deallocate(displacements,"displacements",CAT_DMATRIX)

  call tr_allocate(counts,1,n_cpu,"counts",CAT_DMATRIX)
  call tr_allocate(displacements,1,n_cpu,"displacements",CAT_DMATRIX)

  call MPI_Allgather(mumps_par%nz_loc,1,MPI_INTEGER_ALL,counts,1,MPI_INTEGER_ALL,MPI_COMM_WORLD,ierr)

  displacements(1) = 0
  do i=2,n_cpu
    displacements(i) = displacements(i-1) + counts(i-1)
  enddo

  if (associated(mumps_par%irn)) call tr_deallocatep(mumps_par%irn,"mumps_par%IRN",CAT_DMATRIX)
  if (associated(mumps_par%jcn)) call tr_deallocatep(mumps_par%jcn,"mumps_par%JCN",CAT_DMATRIX)
  if (associated(mumps_par%a) )  call tr_deallocatep(mumps_par%a,"mumps_par%A",CAT_DMATRIX)
  if (associated(mumps_par%rhs)) call tr_deallocatep(mumps_par%rhs,"mumps_par%rhs",CAT_DMATRIX)

  call tr_allocatep(mumps_par%irn,Int1,nnz,"mumps_par%IRN",CAT_DMATRIX)
  call tr_allocatep(mumps_par%jcn,Int1,nnz,"mumps_par%JCN",CAT_DMATRIX)
  call tr_allocatep(mumps_par%a,Int1,nnz,"mumps_par%A",CAT_DMATRIX)
  call tr_allocatep(mumps_par%rhs,Int1,n,"mumps_par%rhs",CAT_DMATRIX)

  call split_allgathersolve(n_cpu,my_id,counts,displacements)

  call MPI_AllReduce(rhs_glob,mumps_par%rhs,mumps_par%n,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  
  call clck_time(t1)
  call clck_ldiff(t0,t1,tsecond)
  if (my_id .eq. 0)  write(*,FMT_TIMING) my_id, '## Elapsed time mpi_gather :', tsecond
  
  if (.not. spss_initialized) then
    call strumpack_init(MPI_COMM_WORLD)
    spss_initialized = .true.
  endif

  call strumpack_set_mat(n,nnz,mumps_par%irn,mumps_par%jcn,mumps_par%a,1,MPI_COMM_WORLD,&
                         UPDATE=spss_analyzed,DISTRIBUTED=.false.,EQUILIBRIUM=.false.)
  call tr_deallocatep(mumps_par%irn,"mumps_par%IRN",CAT_DMATRIX)
  call tr_deallocatep(mumps_par%jcn,"mumps_par%JCN",CAT_DMATRIX)
  call tr_deallocatep(mumps_par%a,"mumps_par%A",CAT_DMATRIX)
  
  if (.not. spss_analyzed) then
    call clck_time(t0)
    call strumpack_analyze(MPI_COMM_WORLD)    
    spss_analyzed = .true.
    call clck_time(t1)
    call clck_ldiff(t0,t1,tsecond)
    if (my_id .eq. 0)  write(*,FMT_TIMING) my_id, '## Elapsed time analysis :', tsecond    
  endif
  
  call clck_time(t0)

  call strumpack_factorize(MPI_COMM_WORLD)   
  call strumpack_solve(n,mumps_par%rhs,MPI_COMM_WORLD)
 
  call clck_time(t1)
  call clck_ldiff(t0,t1,tsecond)
  if (my_id .eq. 0)  write(*,FMT_TIMING) my_id, '## Elapsed facto/solve :', tsecond  
 
  do k=1,n
    deltas(k) =  mumps_par%rhs(k)
  enddo

  if (allocated(counts))        call tr_deallocate(counts,"counts",CAT_DMATRIX)
  if (allocated(displacements)) call tr_deallocate(displacements,"displacements",CAT_DMATRIX)

  return
end subroutine solve_strumpack_all
#endif  
