! This program is the driver of the accept reject function tests
program accept_reject_funct_test_driver
use fruit
use mod_accept_reject_funct_test, only: run_fruit_accept_reject_funct
  implicit none

  ! init fruit suite
  call init_fruit
  call init_fruit_xml

  ! run the accept reject funct test basket
  call run_fruit_accept_reject_funct

  ! write test summary and finilize test suit
  call fruit_summary
  call fruit_summary_xml
  call fruit_finalize

end program accept_reject_funct_test_driver
