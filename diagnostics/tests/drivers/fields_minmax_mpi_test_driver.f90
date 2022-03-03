! This program is the driver of the fields minmax tests
program fields_minmax_mpi_test_driver
use fruit
use fruit_mpi
use mod_mpi_tools, only: init_mpi_threads,finalize_mpi_threads
use mod_fields_minmax_mpi_test, only: run_fruit_fields_minmax_mpi
  implicit none
  integer :: rank,n_tasks,ifail

  !> initialise MPI
  call init_mpi_threads(rank,n_tasks,ifail)

  ! init fruit suite
  call fruit_init_mpi_xml(rank)

  ! run the fields minmax test basket
  call run_fruit_fields_minmax_mpi(rank,n_tasks,ifail)

  ! write test summary and finilize test suit
  call fruit_summary_mpi(n_tasks,rank)
  call fruit_summary_mpi_xml(n_tasks,rank)
  call fruit_finalize_mpi(n_tasks,rank)

  !> finalize mpi
  call finalize_mpi_threads(ifail)

end program fields_minmax_mpi_test_driver
