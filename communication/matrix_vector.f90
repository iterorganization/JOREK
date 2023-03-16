!> Sparse matrix vector product
!! y = Ax
subroutine matrix_vector(x_vec, a_mat, y_vec)
  use mpi_mod
  use mod_integer_types
  use data_structure, only: type_SP_MATRIX, type_RHS

  implicit none
  
  type(type_SP_MATRIX)  :: a_mat
  type(type_RHS)        :: x_vec
  type(type_RHS)        :: y_vec

  ! --- Local variables
  real*8                :: t1, t2, t3, t4, t5
  real*8, allocatable   :: y_tmp(:), y_tmp_block(:)
  integer,allocatable   :: recv_counts(:), recv_disp(:)
  integer               :: n, i, ir, jc
  integer               :: my_id, n_cpu, ierr, counts
  integer(kind=int_all) :: n_blocksize, n_blocks, iA_start, ix_start, iy_start, index_offset, Int_tmp
  integer               :: ndof_local
  external              :: dgemv
  
  call MPI_COMM_SIZE(a_mat%comm, n_cpu, ierr)
  call MPI_COMM_RANK(a_mat%comm, my_id, ierr) 

  counts = x_vec%n
  if (.not.associated(y_vec%val)) allocate(y_vec%val(y_vec%n))
  y_vec%val(1:y_vec%n)  = 0.d0
  n_blocksize  = a_mat%block_size
  n_blocks     = a_mat%nnz/n_blocksize**2
  ndof_local   = (a_mat%index_max(my_id + 1) - a_mat%index_min(my_id + 1) + 1)*n_blocksize
  
  index_offset = (a_mat%index_min(my_id + 1) - 1)*n_blocksize
  
  allocate(y_tmp(ndof_local))
  y_tmp        = 0.d0
  allocate(y_tmp_block(a_mat%block_size))
  
! --- The actual matrix vector multiplication uses dense matrix-vector products for the small
!     dense blocks within our sparse matrix. The size of these blocks depends on n_tor. Depending on
!     this block size (so depending on n_tor), two slightly different kernels are implemented.
!     Warning: this simply seg-faults with long-integers.

!$omp parallel default(none) &
!$omp   shared(y_tmp, a_mat, x_vec, n_blocks, n_blocksize, index_offset)  &
!$omp   private(i,iA_start,ix_start, iy_start, ir, jc, y_tmp_block)
!$omp do schedule(guided)
  do i = 1, n_blocks

    iA_start = (i-1) * n_blocksize**2
    ix_start = a_mat%jcn(iA_start+1)
    iy_start = a_mat%irn(iA_start+1) - index_offset

    call dgemv('T',n_blocksize,n_blocksize,1.d0,a_mat%val(iA_start+1),n_blocksize,x_vec%val(ix_start),Int1,0.d0,y_tmp_block,Int1)

!$omp critical
    y_tmp(iy_start:iy_start+n_blocksize-1) = y_tmp(iy_start:iy_start+n_blocksize-1) + y_tmp_block(1:n_blocksize)
!$omp end critical

  enddo
!$omp end do
!$omp end parallel

  y_vec%val(1:y_vec%n) = 0.d0

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

  call mpi_allgatherv(y_tmp, ndof_local, MPI_DOUBLE_PRECISION, y_vec%val, recv_counts, recv_disp, MPI_DOUBLE_PRECISION, a_mat%comm,ierr)
  
  deallocate(y_tmp,y_tmp_block)
  deallocate(recv_counts)
  deallocate(recv_disp)
 
  return

end subroutine matrix_vector