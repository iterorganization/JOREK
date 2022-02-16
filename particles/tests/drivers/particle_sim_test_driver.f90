! This program is the driver of the particle_sim tests (NO MPI)
program particle_sim_test_driver
use fruit
use mod_particle_sim_test, only: run_fruit_particle_sim
  implicit none

  ! init fruit suite
  call init_fruit
  call init_fruit_xml

  ! run the particle sim test basket
  call run_fruit_particle_sim

  ! write test summary and finilize test suit
  call fruit_summary
  call fruit_summary_xml
  call fruit_finalize

end program particle_sim_test_driver
