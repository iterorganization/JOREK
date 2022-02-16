! This program is the driver of the particle_type_openmp tests
program particle_type_openmp_test_driver
use fruit
use mod_particle_type_openmp, only: run_fruit_particle_type_openmp
  implicit none

  ! init fruit suite
  call init_fruit
  call init_fruit_xml

  ! run the particle type openmp test basket
  call run_fruit_particle_type_openmp

  ! write test summary and finilize test suit
  call fruit_summary
  call fruit_summary_xml
  call fruit_finalize

end program particle_type_openmp_test_driver
