! This program is the driver of the root finding tests
program rootfinding_test_driver
use fruit
use mod_rootfinding_test, only: run_fruit_rootfinding
  implicit none

  ! init fruit suite
  call init_fruit
  call init_fruit_xml

  ! run the root finding test basket
  call run_fruit_rootfinding

  ! write test summary and finilize test suit
  call fruit_summary
  call fruit_summary_xml
  call fruit_finalize

end program rootfinding_test_driver
