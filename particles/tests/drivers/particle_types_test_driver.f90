! This program is the driver of the particle_types tests
program particle_types_test_driver
use fruit
use mod_particle_types_test, only: run_fruit_particle_types
  implicit none

  ! init fruit suite
  call init_fruit
  call init_fruit_xml

  ! run the particle types test basket
  call run_fruit_particle_types

  ! write test summary and finilize test suit
  call fruit_summary
  call fruit_summary_xml
  call fruit_finalize

end program particle_types_test_driver
