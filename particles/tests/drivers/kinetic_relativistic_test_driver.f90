! This program is the driver of the kinetic_relativistic tests
program particle_kinetic_relativistic_test_driver
use fruit
use mod_kinetic_relativistic_test, only: run_fruit_kinetic_relativistic
  implicit none

  ! init fruit suite
  call init_fruit
  call init_fruit_xml

  ! run the kinetic relativistic test basket
  call run_fruit_kinetic_relativistic

  ! write test summary and finilize test suit
  call fruit_summary
  call fruit_summary_xml
  call fruit_finalize

end program particle_kinetic_relativistic_test_driver
