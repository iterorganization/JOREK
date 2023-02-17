subroutine matrix_split_bcast(ad_mat, ac_mat)

  use mpi_mod
  use mod_integer_types
  use data_structure, only: type_SP_MATRIX
  use tr_module
  implicit none

  type(type_SP_MATRIX), intent(inout)  :: ad_mat
  type(type_SP_MATRIX), intent(inout)  :: ac_mat

  integer               :: my_id, n_cpu, comm, ierr
  integer(kind=int_all) :: is, ie, nnz
  integer               :: buff_max, nz_split_end, n_split, i

  comm = ad_mat%comm
  call MPI_COMM_RANK(comm, my_id, ierr)
  call MPI_COMM_SIZE(comm, n_cpu, ierr)

  nnz = ad_mat%nnz
  call MPI_Bcast(nnz,1,MPI_INTEGER_ALL,0,comm,ierr)

  if (my_id.eq.0) then
    call ad_mat%move_to(ac_mat, with_data=.true.) ! move ad_mat to ac_mat for my_id.eq.0
  else
    call ad_mat%move_to(ac_mat, with_data=.false.) ! copy matrix parameters except arrays
    ac_mat%nnz = nnz
    if (associated(ac_mat%irn)) call tr_deallocatep(ac_mat%irn,"irn",CAT_DMATRIX)
    if (associated(ac_mat%jcn)) call tr_deallocatep(ac_mat%jcn,"jcn",CAT_DMATRIX)
    if (associated(ac_mat%val)) call tr_deallocatep(ac_mat%val,"val",CAT_DMATRIX)
    call tr_allocatep(ac_mat%irn,Int1,ac_mat%nnz,"irn",CAT_DMATRIX)
    call tr_allocatep(ac_mat%jcn,Int1,ac_mat%nnz,"jcn",CAT_DMATRIX)
    call tr_allocatep(ac_mat%val,Int1,ac_mat%nnz,"val",CAT_DMATRIX)
  endif

  buff_max = INT_MAX

  if (nnz > buff_max) then
    n_split = nnz/buff_max
    do i=1,n_split
      is = (i-1)*buff_max+1
      ie = i*buff_max
        call MPI_BCAST(ac_mat%irn(is:ie),buff_max,MPI_INTEGER_ALL,0,comm,ierr)
        call MPI_BCAST(ac_mat%jcn(is:ie),buff_max,MPI_INTEGER_ALL,0,comm,ierr)
        call MPI_BCAST(ac_mat%val(is:ie),buff_max,MPI_DOUBLE_PRECISION,0,comm,ierr)
    enddo

    nz_split_end = mod(nnz,buff_max)
    is = n_split*buff_max + 1
    ie = nnz
    call MPI_BCAST(ac_mat%irn(is:ie),nz_split_end,MPI_INTEGER_ALL,0,comm,ierr)
    call MPI_BCAST(ac_mat%jcn(is:ie),nz_split_end,MPI_INTEGER_ALL,0,comm,ierr)
    call MPI_BCAST(ac_mat%val(is:ie),nz_split_end,MPI_DOUBLE_PRECISION,0,comm,ierr)

  else
    call MPI_BCAST(ac_mat%irn(1:nnz),nnz,MPI_INTEGER_ALL,0,comm,ierr)
    call MPI_BCAST(ac_mat%jcn(1:nnz),nnz,MPI_INTEGER_ALL,0,comm,ierr)
    call MPI_BCAST(ac_mat%val(1:nnz),nnz,MPI_DOUBLE_PRECISION,0,comm,ierr)
  endif

end subroutine matrix_split_bcast
