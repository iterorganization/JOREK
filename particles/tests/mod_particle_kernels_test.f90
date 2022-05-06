#ifdef CUDA_KERNELS
module mod_particle_kernels_test
  use fruit
  use cudafor
  use mpi
  use mod_particle_kernels
  implicit none
  
  private
  public :: run_fruit_particle_kernels


contains

  subroutine run_fruit_particle_kernels()
      write(*,'(/A)') "  ... setting-up: particle kernels tests"
      call setup
      write(*,'(/A)') "  ... running: particle kernels tests"
      call test_particle_kinetic_leapfrog_loop
      write(*,'(/A)') "  ... tearing-down: particle kernels tests"
      call teardown
  end subroutine run_fruit_particle_kernels

  !> Set-up and tear-down -------------------------
  !> set-up particle types test features
  subroutine setup()
  end subroutine setup

  subroutine teardown()
  end subroutine teardown

  subroutine test_particle_kinetic_leapfrog_loop
  end subroutine test_particle_kinetic_leapfrog_loop    
  
end module mod_particle_kernels_test
#endif
