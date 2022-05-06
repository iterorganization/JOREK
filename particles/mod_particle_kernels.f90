#ifdef CUDA_KERNELS
!> Contains routines for offloading particle computation to GPU devices using CUDA kernels
!> At present only the particle_kinetic_leapfrog type is supported

module mod_particle_kernels

  use mod_particle_sim
  use mod_particle_types, only: particle_kinetic_leapfrog
  use cudafor
  use mpi
  
  implicit none

  private
  public particle_kinetic_leapfrog_loop
  
contains

  subroutine particle_kinetic_leapfrog_loop

  end subroutine particle_kinetic_leapfrog_loop
  
end module mod_particle_kernels
#endif
