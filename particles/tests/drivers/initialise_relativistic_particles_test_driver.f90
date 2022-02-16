! This program is the driver of the initialise relativistic particles tests
program initialise_relativistic_particles_test_driver
use fruit
use mod_initialise_relativistic_particles_test, only: run_fruit_initialise_relativistic_particles
  implicit none

  ! init fruit suite
  call init_fruit
  call init_fruit_xml

  ! run the initialise relativistic particles test basket
  call run_fruit_initialise_relativistic_particles

  ! write test summary and finilize test suit
  call fruit_summary
  call fruit_summary_xml
  call fruit_finalize

end program initialise_relativistic_particles_test_driver

