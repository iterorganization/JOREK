!> Sparse matrix vector product using coordinate scheme to be called only on all of MPI_COMM_WORLD
!!
!! result (y) is only known on id=0
subroutine gmres_matrix_vector(size_x, x, size_y, y, my_id)

use tr_module
use mod_parameters
use global_distributed_matrix
use mpi_mod
use mod_integer_types

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

end if
! --- End: Two different kernels for matrix-vector multiplication depending on n_tor

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
