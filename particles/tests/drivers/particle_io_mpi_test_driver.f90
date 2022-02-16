! This program is the driver of the particle_io_mpi tests
program particle_io_mpi_test_driver
use fruit_mpi
use mod_mpi_tools,            only: init_mpi_threads
use mod_mpi_tools,            only: finalize_mpi_threads
use mod_particle_io_mpi_test, only: run_fruit_particle_io_mpi
  implicit none
  integer :: rank,n_tasks,ifail

  ! initialize the mpi comm_world
  call init_mpi_threads(rank,n_tasks,ifail)

  ! init fruit suite
  call fruit_init_mpi_xml(rank)

  ! run the particle type openmp test basket
  call run_fruit_particle_io_mpi(rank,n_tasks,ifail)

  ! write test summary and finilize test suit
  call fruit_summary_mpi(n_tasks,rank)
  call fruit_summary_mpi_xml(n_tasks,rank)
  call fruit_finalize_mpi(n_tasks,rank)

  ! close the mpi communicator
  call finalize_mpi_threads(ifail)
end program particle_io_mpi_test_driver
