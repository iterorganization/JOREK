module setup_mpi
contains
subroutine setup
  use mpi
  integer :: ierr, provided
  call MPI_init_thread(MPI_THREAD_SINGLE, provided, ierr)
  write(*,*) "MPI setup: ", provided, ierr
end subroutine setup

subroutine teardown
  use mpi
  integer :: ierr
  call MPI_finalize(ierr)
end subroutine teardown
end module setup_mpi
