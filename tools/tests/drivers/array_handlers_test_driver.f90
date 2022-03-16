!> This program is the driver of the array handlers tests
program array_handlers_test_driver
use fruit
use mod_array_handlers_test, only: run_fruit_array_handlers
  implicit none

  ! init fruit suite
  call init_fruit
  call init_fruit_xml

  ! run the array handlers test basket
  call run_fruit_array_handlers

  ! write test summary and finilize test suit
  call fruit_summary
  call fruit_summary_xml
  call fruit_finalize

end program array_handlers_test_driver
