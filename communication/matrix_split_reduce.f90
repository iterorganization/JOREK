!> Reduces distributed matrix ad_mat into centralized matrix ac_mat
! using split MPI_Allgatherv if MPI counts beyond INT_MAX
subroutine matrix_split_reduce(ad_mat, ac_mat)

  use mpi_mod
  use mod_integer_types
  use data_structure, only: type_SP_MATRIX
  use tr_module
  implicit none

  type(type_SP_MATRIX), intent(inout)  :: ad_mat
  type(type_SP_MATRIX), intent(inout)  :: ac_mat


  integer               :: my_id, n_cpu, comm
  integer               :: i, i_cpu, ierr
  integer(kind=int_all) :: is, ie
  integer, allocatable  :: counts_short(:,:), displs_short(:,:)
  integer, allocatable  :: index_buffer(:,:), index_target(:,:)
  integer(kind=int_all),allocatable :: counts_long(:), displs_long(:)

  integer(kind=int_all) :: i_long

  logical                :: need_to_split
  integer                :: n_split, i_split
  integer(kind=int_all)  :: count_split


  comm = ad_mat%comm
  call MPI_COMM_RANK(comm, my_id, ierr)
  call MPI_COMM_SIZE(comm, n_cpu, ierr)

  allocate(counts_long(n_cpu))
  allocate(displs_long(n_cpu))

  call MPI_Allgather(ad_mat%nnz,1,MPI_INTEGER_ALL,counts_long,1,MPI_INTEGER_ALL,comm,ierr)

  displs_long(1) = 0
  do i=2,n_cpu
    displs_long(i) = displs_long(i-1) + counts_long(i-1)
  enddo

  call ad_mat%move_to(ac_mat, with_data=.false.) ! copy matrix parameters except arrays
  ac_mat%reduced = .true.

  ! Allocate centralized matrix
  if (associated(ac_mat%irn)) call tr_deallocatep(ac_mat%irn,"irn",CAT_DMATRIX)
  if (associated(ac_mat%jcn)) call tr_deallocatep(ac_mat%jcn,"jcn",CAT_DMATRIX)
  if (associated(ac_mat%val)) call tr_deallocatep(ac_mat%val,"val",CAT_DMATRIX)

  ac_mat%nnz = sum(counts_long(1:n_cpu))

  ! --- Check if we need to split
  if (ac_mat%nnz .gt. INT_MAX) then
    need_to_split = .true.
    n_split = ac_mat%nnz / INT_MAX + 1
  else
    need_to_split = .false.
    n_split = 1
  endif

  ! --- Allocate short-integer counts and displacements for MPI calls
  ! --- Counts still need to be copied because MPI count types are always short ints
  allocate(counts_short(n_split,n_cpu))
  allocate(displs_short(n_split,n_cpu))


  ! --- Split MPI calls
  if (need_to_split) then
    ! --- Split respective to the max send/recv

    if (my_id .eq. 0) write(*,*) 'Warning: splitting matrix MPI centralisation', n_split

    allocate(index_buffer(n_split,n_cpu))
    allocate(index_target(n_split,n_cpu))

    do i_split=1,n_split

      ! --- Split counts for each MPI chunk
      do i_cpu=1,n_cpu
        count_split = counts_long(i_cpu) / n_split
        if (i_split .gt. 1) then
          index_buffer(i_split,i_cpu) = (i_split-1)*count_split
        else
          index_buffer(i_split,i_cpu) = 0
        endif
        counts_short(i_split,i_cpu) = count_split
        if (i_split .eq. n_split) then
          counts_short(i_split,i_cpu) = counts_long(i_cpu) - (n_split-1)*count_split
        endif
      enddo

      ! --- Split displacements for each MPI chunk
      displs_short(i_split,1:n_cpu) = 0
      index_target(i_split,1:n_cpu) = 0
      count_split = counts_long(1) / n_split
      index_target(i_split,1) = displs_long(1) + (i_split-1)*count_split
      do i_cpu=2,n_cpu
        displs_short(i_split,i_cpu) = displs_short(i_split,i_cpu-1) + counts_short(i_split,i_cpu-1)
        count_split = counts_long(i_cpu) / n_split
        index_target(i_split,i_cpu) = displs_long(i_cpu) + (i_split-1)*count_split
      enddo

    enddo

! do irn, jcn, val one by one to save memory
    i_cpu = my_id+1

    call tr_allocatep(ac_mat%irn,Int1,ac_mat%nnz,"irn",CAT_DMATRIX)
    do i_split=1,n_split
      is = index_buffer(i_split,i_cpu) + 1
      ie = index_buffer(i_split,i_cpu) + counts_short(i_split,i_cpu)
      call MPI_AllgatherV(ad_mat%irn(is:ie),counts_short(i_split,i_cpu),MPI_INTEGER_ALL,ac_mat%irn, &
                          counts_short(i_split,1:n_cpu),index_target(i_split,1:n_cpu),MPI_INTEGER_ALL,comm,ierr)
    enddo
    call tr_deallocatep(ad_mat%irn,"irn",CAT_DMATRIX)

    call tr_allocatep(ac_mat%jcn,Int1,ac_mat%nnz,"jcn",CAT_DMATRIX)
    do i_split=1,n_split
      is = index_buffer(i_split,i_cpu) + 1
      ie = index_buffer(i_split,i_cpu) + counts_short(i_split,i_cpu)
      call MPI_AllgatherV(ad_mat%jcn(is:ie),counts_short(i_split,i_cpu),MPI_INTEGER_ALL,ac_mat%jcn, &
                          counts_short(i_split,1:n_cpu),index_target(i_split,1:n_cpu),MPI_INTEGER_ALL,comm,ierr)
    enddo
    call tr_deallocatep(ad_mat%jcn,"jcn",CAT_DMATRIX)

    call tr_allocatep(ac_mat%val,Int1,ac_mat%nnz,"val",CAT_DMATRIX)
    do i_split=1,n_split
      is = index_buffer(i_split,i_cpu) + 1
      ie = index_buffer(i_split,i_cpu) + counts_short(i_split,i_cpu)
      call MPI_AllgatherV(ad_mat%val(is:ie),counts_short(i_split,i_cpu),MPI_DOUBLE_PRECISION,ac_mat%val, &
                          counts_short(i_split,1:n_cpu),index_target(i_split,1:n_cpu),MPI_DOUBLE_PRECISION,comm,ierr)
    enddo
    call tr_deallocatep(ad_mat%val,"val",CAT_DMATRIX)


    deallocate(index_buffer)
    deallocate(index_target)

  ! --- Don't split MPI calls
  else

    ! --- Counts still need to be copied because MPI count types are always short ints
    counts_short(1,1:n_cpu) = counts_long(1:n_cpu)
    displs_short(1,1:n_cpu) = displs_long(1:n_cpu)

    call tr_allocatep(ac_mat%irn,Int1,ac_mat%nnz,"irn",CAT_DMATRIX)
    call MPI_AllgatherV(ad_mat%irn,counts_short(1,my_id+1),MPI_INTEGER_ALL,ac_mat%irn, &
                        counts_short(1,1:n_cpu),displs_short(1,1:n_cpu),MPI_INTEGER_ALL,comm,ierr)
    call tr_deallocatep(ad_mat%irn,"irn",CAT_DMATRIX)

    call tr_allocatep(ac_mat%jcn,Int1,ac_mat%nnz,"jcn",CAT_DMATRIX)
    call MPI_AllgatherV(ad_mat%jcn,counts_short(1,my_id+1),MPI_INTEGER_ALL,ac_mat%jcn, &
                        counts_short(1,1:n_cpu),displs_short(1,1:n_cpu),MPI_INTEGER_ALL,comm,ierr)
    call tr_deallocatep(ad_mat%jcn,"jcn",CAT_DMATRIX)

    call tr_allocatep(ac_mat%val,Int1,ac_mat%nnz,"val",CAT_DMATRIX)
    call MPI_AllgatherV(ad_mat%val,counts_short(1,my_id+1),MPI_DOUBLE_PRECISION,ac_mat%val, &
                        counts_short(1,1:n_cpu),displs_short(1,1:n_cpu),MPI_DOUBLE_PRECISION,comm,ierr)
    call tr_deallocatep(ad_mat%val,"val",CAT_DMATRIX)

  endif

  deallocate(counts_short)
  deallocate(displs_short)
  deallocate(counts_long)
  deallocate(displs_long)

  return
end subroutine matrix_split_reduce