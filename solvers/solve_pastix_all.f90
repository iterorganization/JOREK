#if defined(USE_PASTIX)
subroutine solve_pastix_all(n_cpu,my_id,index_min,index_max)
!---------------------------------------------------------------------
! subroutine solves the complete system of equation using pastix with
! distributed matrix on the main group mpi_comm_world
!---------------------------------------------------------------------
use tr_module 
use mod_parameters
use mumps_module
use pastix_module
use global_distributed_matrix
use mpi_mod
use mod_clock
use mod_coicsr
use phys_module, only: use_BLR_compression, epsilon_BLR, just_in_time_BLR, pastix_blr_abs_tol
use mod_integer_types
 
implicit none

#ifdef USE_PASTIX
#include "pastix_fortran.h"
#else
#include "no_pastix_fortran.h"
#endif

! --- Routine variables
integer                           :: n_cpu, my_id, index_min, index_max       ! global index_min, index_max for this cpu
! --- Local variables
real*8,               allocatable :: column_local(:)
type(clcktype)                    :: t_itstart, t0, t1, t2, t3
real*8                            :: tsecond
integer                           :: i, k, j, ierr
integer(kind=int_all)             :: m_loc
integer(kind=int_all),allocatable :: counts(:), displacements(:)
integer(kind=int_all), parameter  :: Int1=1

!write(*,*) my_id,'*********************************'
!write(*,*) my_id,'*  solve global matrix (PastiX) *'
!write(*,*) my_id,'*********************************'

m_loc = (index_max - index_min + 1) * n_tor * n_var
mumps_par%nz_loc = nz_glob

call MPI_Allreduce(m_loc,mumps_par%N,1,MPI_INTEGER_ALL,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_Allreduce(mumps_par%NZ_loc,mumps_par%nz,1,MPI_INTEGER_ALL,MPI_SUM,MPI_COMM_WORLD,ierr)

!------------------------------------------------------- colunm scaling of global distributed matrix
call clck_time(t0)

if (allocated(column_scaling))  call tr_deallocate(column_scaling,"column_scaling",CAT_DMATRIX)
if (allocated(column_local))    call tr_deallocate(column_local,"column_local",CAT_DMATRIX)
call tr_allocate(column_scaling,Int1,mumps_par%N,"column_scaling",CAT_DMATRIX)
call tr_allocate(column_local,Int1,mumps_par%N,"column_local",CAT_DMATRIX)

column_local = 1.d-20;   column_scaling = 1.d-20
do k=1,nz_glob
  j = jcn_glob(k)
  column_local(j) = max(column_local(j),abs(A_glob(k)))
enddo

call MPI_AllReduce(column_local,column_scaling,mumps_par%N,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)

do k = 1, nz_glob
  j = jcn_glob(k)
  A_glob(k) = A_glob(k) / column_scaling(j)
enddo

call clck_time(t1)
call clck_ldiff(t0,t1,tsecond)
if (my_id .eq. 0)  write(*,FMT_TIMING) my_id, '## Elapsed time scale :', tsecond

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

if (associated(mumps_par%IRN)) call tr_deallocatep(mumps_par%IRN,"mumps_par%IRN",CAT_DMATRIX)
if (associated(mumps_par%JCN)) call tr_deallocatep(mumps_par%JCN,"mumps_par%JCN",CAT_DMATRIX)
if (associated(mumps_par%A) )  call tr_deallocatep(mumps_par%A,"mumps_par%A",CAT_DMATRIX)
if (associated(mumps_par%rhs)) call tr_deallocatep(mumps_par%rhs,"mumps_par%rhs",CAT_DMATRIX)

call tr_allocatep(mumps_par%IRN,Int1,mumps_par%nz,"mumps_par%IRN",CAT_DMATRIX)
call tr_allocatep(mumps_par%JCN,Int1,mumps_par%nz,"mumps_par%JCN",CAT_DMATRIX)
call tr_allocatep(mumps_par%A,Int1,mumps_par%nz,"mumps_par%A",CAT_DMATRIX)
call tr_allocatep(mumps_par%rhs,Int1,mumps_par%n,"mumps_par%rhs",CAT_DMATRIX)

call split_allgathersolve(n_cpu,my_id,counts,displacements)

call MPI_AllReduce(RHS_glob,mumps_par%RHS,mumps_par%N,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

call clck_time(t1)
call clck_ldiff(t0,t1,tsecond)
if (my_id .eq. 0)  write(*,FMT_TIMING) my_id, '## Elapsed time mpi_gather :', tsecond
call clck_time(t0)

#ifdef USE_BLOCK
!---------------------------- reduce IRN,JCN to make use of blocksize ntor*nvar
!                             temporary solution before using blocks everywhere

block_size  = n_tor * n_var
block_size2 = block_size**2
!---------------------------- reduce IRN,JCN to make use of blocksize ntor*nvar
n_block   = mumps_par%n  / block_size
nnz_block = mumps_par%nz / block_size2

do i=1,nnz_block  
  mumps_par%irn(i) = (mumps_par%irn((i-1)*block_size2+1) - 1) / block_size + 1 
  mumps_par%jcn(i) = (mumps_par%jcn((i-1)*block_size2+1) - 1) / block_size + 1 
enddo

if (allocated(sparskit_work)) deallocate(sparskit_work)
allocate(sparskit_work(n_block+1))

call coicsr2(n_block,nnz_block,mumps_par%A,mumps_par%IRN(1:nnz_block),mumps_par%JCN(1:nnz_block),block_size,sparskit_work)

! -- For PaStiX solver before version 6.x
if (.not. allocated(pastix_perm_vars))  call tr_allocate(pastix_perm_vars,Int1,n_block,"pastix_perm_vars",CAT_UNKNOWN)
if (.not. allocated(pastix_iperm_vars)) call tr_allocate(pastix_iperm_vars,Int1,n_block,"pastix_iperm_vars",CAT_UNKNOWN)

#else /* USE_BLOCK */

if (allocated(sparskit_work)) deallocate(sparskit_work)
allocate(sparskit_work(mumps_par%N + 1))

call coicsr(mumps_par%N,mumps_par%NZ,1,mumps_par%A,mumps_par%IRN,mumps_par%JCN,sparskit_work)
#endif /* USE_BLOCK */

if (allocated(sparskit_work)) deallocate(sparskit_work)
call clck_time(t1)
call clck_ldiff(t0,t1,tsecond)
if (my_id .eq. 0)  write(*,FMT_TIMING) my_id, '## Elapsed time coicsr :', tsecond

! -- For PaStiX solver before version 6.x
if (.not. allocated(pastix_perm_vars))  call tr_allocate(pastix_perm_vars,Int1,mumps_par%n,"pastix_perm_vars",CAT_UNKNOWN)
if (.not. allocated(pastix_iperm_vars)) call tr_allocate(pastix_iperm_vars,Int1,mumps_par%n,"pastix_iperm_vars",CAT_UNKNOWN)

call pastix_init_num_threads(my_id)

if (.not. pastix_initialised) then


  ! -- For PaStiX solver before version 6.x
  pastix_iparm(IPARM_MODIFY_PARAMETER)  = API_NO          ! insert default values
  pastix_iparm(IPARM_START_TASK)        = API_TASK_INIT   ! initializse
  pastix_iparm(IPARM_END_TASK)          = API_TASK_INIT

  if (my_id .eq. 0) then
    write(*,*) '***********************************'
    write(*,*) '* initialise PastiX               *'
    write(*,*) '***********************************'
  endif

  ! -- For PaStiX solver before version 6.x
#ifdef USE_BLOCK
  call pastix_fortran(pastix_data,MPI_COMM_WORLD,n_block,mumps_par%jcn,mumps_par%irn,mumps_par%A, &
                        pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,Int1,pastix_iparm,pastix_dparm)
#else
  call pastix_fortran(pastix_data,MPI_COMM_WORLD,mumps_par%n,mumps_par%jcn,mumps_par%irn,mumps_par%A, &
                        pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,Int1,pastix_iparm,pastix_dparm)
#endif

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

  ! -- For PaStiX solver before version 6.x
  pastix_iparm(IPARM_RHS_MAKING)            = pastix_rhs                 ! right hand side (0 : use RHS)
  pastix_iparm(IPARM_SYM)                   = pastix_sym
  pastix_iparm(IPARM_AMALGAMATION_LEVEL)    = pastix_amalg
#ifdef WORLDWAR2
#ifdef FUNNELED
  pastix_iparm(IPARM_THREAD_COMM_MODE)      = API_THREAD_FUNNELED
#else
  pastix_iparm(IPARM_THREAD_COMM_MODE)      = API_THREAD_MULTIPLE
#endif
#endif

  pastix_initialised = .true.
endif

if (.not. pastix_analysed) then

  ! -- For PaStiX solver before version 6.x
  pastix_iparm(IPARM_THREAD_NBR) = pastix_nthrd
  pastix_iparm(IPARM_START_TASK) = API_TASK_ORDERING
  pastix_iparm(IPARM_END_TASK)   = API_TASK_ANALYSE

  call clck_time(t0)
  
  if (my_id .eq. 0) then
    write(*,*) '***********************************'
    write(*,*) '* analyse PastiX                  *'
    write(*,*) '***********************************'
  endif

  ! -- For PaStiX solver before version 6.x
#ifdef USE_BLOCK
  
  call pastix_fortran(pastix_data,MPI_COMM_WORLD, n_block, &
                      mumps_par%jcn(1:n_block+1), mumps_par%irn(1:nnz_block), mumps_par%A(1:mumps_par%nz), &
                      pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,Int1,pastix_iparm,pastix_dparm)

#else

  call pastix_fortran(pastix_data,MPI_COMM_WORLD, mumps_par%n, &
                      mumps_par%jcn(1:mumps_par%n+1), mumps_par%irn(1:mumps_par%nz), mumps_par%A(1:mumps_par%nz), &
                      pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,Int1,pastix_iparm,pastix_dparm)

#endif
 
  call clck_time(t1)
  call clck_ldiff(t0,t1,tsecond)
  if (my_id .eq. 0)  write(*,FMT_TIMING) my_id, '## Elapsed time analysis :', tsecond

  pastix_analysed = .true.

endif ! .not. pastix_analysed

call clck_time(t0)


! -- For PaStiX solver before version 6.x
pastix_iparm(IPARM_THREAD_NBR) = pastix_nthrd
pastix_iparm(IPARM_START_TASK) = API_TASK_NUMFACT
pastix_iparm(IPARM_END_TASK)   = pastix_endsolve
#ifdef WORLDWAR2
#ifdef FUNNELED
  pastix_iparm(IPARM_THREAD_COMM_MODE)  = API_THREAD_FUNNELED
#else
  pastix_iparm(IPARM_THREAD_COMM_MODE)  = API_THREAD_MULTIPLE
#endif
#endif

if (my_id .eq. 0) then
  write(*,*) '***********************************'
  write(*,*) '* call PastiX                     *'
  write(*,*) '***********************************'
endif


! -- For PaStiX solver before version 6.x
#ifdef USE_BLOCK
pastix_iparm(IPARM_DOF_NBR)            = block_size

call pastix_fortran(pastix_data,MPI_COMM_WORLD, n_block,                                                 &
                    mumps_par%jcn(1:n_block+1), mumps_par%irn(1:nnz_block), mumps_par%A(1:mumps_par%nz), &
                    pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,Int1,pastix_iparm,pastix_dparm)
#else

call pastix_fortran(pastix_data,MPI_COMM_WORLD, mumps_par%n, mumps_par%jcn, mumps_par%irn, mumps_par%A, &
                    pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,Int1,pastix_iparm,pastix_dparm)
#endif

call clck_time(t1)
call clck_ldiff(t0,t1,tsecond)
if (my_id .eq. 0)  write(*,FMT_TIMING) my_id, '## Elapsed facto/solve :', tsecond

do k=1,mumps_par%n
  deltas(k) =  mumps_par%rhs(k)  / column_scaling(k)
enddo

if (allocated(column_local))  call tr_deallocate(column_local,"column_local",CAT_DMATRIX)
if (allocated(counts))        call tr_deallocate(counts,"counts",CAT_DMATRIX)
if (allocated(displacements)) call tr_deallocate(displacements,"displacements",CAT_DMATRIX)

return
end
#endif


#ifdef USE_PASTIX6
!> subroutine solves the complete system of equation using PaStiX 6.2
subroutine solve_pastix_all(n_cpu,my_id,index_min,index_max)
  use mod_pastix

  use tr_module 
  use mod_parameters
  use mumps_module
  use global_distributed_matrix, only: column_scaling, deltas, nz_glob, &
                                       irn_glob, jcn_glob, a_glob, rhs_glob  
  use mpi_mod
  use mod_clock
  use mod_integer_types

!$ use omp_lib

  implicit none

! --- Routine parameters
  integer,               intent(in) :: n_cpu, my_id
  integer,               intent(in) :: index_min, index_max

! --- Local variables
  real*8,               allocatable :: column_local(:)
  type(clcktype)                    :: t_itstart, t0, t1, t2, t3
  real*8                            :: tsecond
  integer                           :: i, k, j, ierr
  integer(kind=int_all)             :: m_loc
  integer(kind=int_all),allocatable :: counts(:), displacements(:)
  integer(kind=int_all), parameter  :: Int1=1
  
!write(*,*) my_id,'*****************************************'
!write(*,*) my_id,'*  solve global matrix using PaStiX 6.2 *'
!write(*,*) my_id,'*****************************************'

  m_loc = (index_max - index_min + 1) * n_tor * n_var
  mumps_par%nz_loc = nz_glob

  call MPI_Allreduce(m_loc,mumps_par%n,1,MPI_INTEGER_ALL,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_Allreduce(mumps_par%nz_loc,mumps_par%nz,1,MPI_INTEGER_ALL,MPI_SUM,MPI_COMM_WORLD,ierr)
  
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

  call tr_allocatep(mumps_par%irn,Int1,mumps_par%nz,"mumps_par%IRN",CAT_DMATRIX)
  call tr_allocatep(mumps_par%jcn,Int1,mumps_par%nz,"mumps_par%JCN",CAT_DMATRIX)
  call tr_allocatep(mumps_par%a,Int1,mumps_par%nz,"mumps_par%A",CAT_DMATRIX)
  call tr_allocatep(mumps_par%rhs,Int1,mumps_par%n,"mumps_par%rhs",CAT_DMATRIX)

  call split_allgathersolve(n_cpu,my_id,counts,displacements)

  call MPI_AllReduce(rhs_glob,mumps_par%rhs,mumps_par%n,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  
  call clck_time(t1)
  call clck_ldiff(t0,t1,tsecond)
  if (my_id .eq. 0)  write(*,FMT_TIMING) my_id, '## Elapsed time mpi_gather :', tsecond
  
  call clck_time(t0)
  
  if (.not. spm_initialized) then
    call pastix_init(MPI_COMM_WORLD)
    spm_initialized = .true.
  endif

  call pastix_set_mat(mumps_par%n,mumps_par%nz,mumps_par%irn,mumps_par%jcn,mumps_par%a,1,MPI_COMM_WORLD,&
                         UPDATE=spm_analyzed,DISTRIBUTED=.false.,EQUILIBRIUM=.false.)
  
  if (.not. spm_analyzed) then
    call clck_time(t0)
    call pastix_analyze()    
    spm_analyzed = .true.
    call clck_time(t1)
    call clck_ldiff(t0,t1,tsecond)
    if (my_id .eq. 0)  write(*,FMT_TIMING) my_id, '## Elapsed time analysis :', tsecond    
  endif
  
  call clck_time(t0)

  call pastix_factorize()   
  call pastix_solve(mumps_par%n,mumps_par%rhs,REFINE=.true.)
 
  call clck_time(t1)
  call clck_ldiff(t0,t1,tsecond)
  if (my_id .eq. 0)  write(*,FMT_TIMING) my_id, '## Elapsed facto/solve :', tsecond

  deltas(1:mumps_par%n) =  mumps_par%rhs(1:mumps_par%n)

  if (allocated(counts))        call tr_deallocate(counts,"counts",CAT_DMATRIX)
  if (allocated(displacements)) call tr_deallocate(displacements,"displacements",CAT_DMATRIX)  
  
  return
end subroutine solve_pastix_all
#endif
