#ifdef CUDA_KERNELS
program particle_kernels_test_driver
use fruit
use mod_particle_kernels_test, only: run_fruit_particle_kernels
  implicit none

  ! init fruit suite
  call init_fruit
!  call init_fruit_xml

  ! run the particle sim test basket
  call run_fruit_particle_kernels

  ! write test summary and finalize test suit
  call fruit_summary
!  call fruit_summary_xml
  call fruit_finalize

end program particle_kernels_test_driver
#endif
