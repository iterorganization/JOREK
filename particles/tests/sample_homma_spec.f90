!> This module contains some testcases for the binary collision model implementation
module sample_homma_spec
use mod_pcg32_rng
use mod_random_seed
use mod_collisions
use constants
use fruit
implicit none

!> Global parameters
real*8, parameter :: n_b = 1d20 !< [m^-3]
real*8, parameter :: kTb = 50.d0*EL_CHG !< [J] (50 eV)
real*8, parameter :: zero3(3) = [0d0, 0d0, 0d0]
real*8, parameter :: grad_kTb = 30.d0*EL_CHG !< [J/m] (30 eV/m)
real*8, parameter :: m_b = 2d0 !< deuterium
integer*1, parameter :: q_b = 1
real*8, parameter :: coulomb_log = 15

contains

!> Verify that the sampling without any temperature gradient produces a symmetrical distribution
subroutine test_sampling_symmetry
type(pcg32_rng) :: rng
integer, parameter :: n_samples = 1000000
real*8  :: v_b(3,n_samples)
real*8  :: u(6,n_samples)
integer :: i

call rng%initialize(6, random_seed(), 1, 1)
do i=1,n_samples
  call rng%next(u(:,i))
end do

call sample_velocity_dist_unmagnetized(n_samples, u, coulomb_log, kTb, zero3, n_b, m_b, q_b, zero3, v_b)

! Check the average of v
call assert_equals(0.d0, sum(v_b(1,:))/n_samples, 200d0, 'x-dir zero')
call assert_equals(0.d0, sum(v_b(2,:))/n_samples, 200d0, 'y-dir zero')
call assert_equals(0.d0, sum(v_b(3,:))/n_samples, 200d0, 'z-dir zero')
end subroutine test_sampling_symmetry


!> Verify that sampling with a temperature gradient produces the right flux of particles in that direction
subroutine test_sampling_flux
use constants, only: PI, TWOPI
use mod_sampling, only: normal_vectors
use mod_sobseq_rng
type(sobseq_rng) :: rng
integer, parameter :: n_samples = 100000
integer, parameter :: n_phi = 10, n_theta=10 ! Try different polar and azimuthal angles
real*8  :: v_b(3,n_samples)
real*8  :: u(6,n_samples)
real*8  :: grad_kT(3)
integer :: i, i_phi, i_theta
real*8 :: phi, theta

real*8 :: v_perp1, v_perp2, v_par
real*8, dimension(3) :: perp1, perp2, par1


call rng%initialize(6, random_seed(), 1, 1)
do i=1,n_samples
  call rng%next(u(:,i))
end do

do i_phi=1,n_phi
  phi = (i_phi-1)*TWOPI/(n_phi-1)
  do i_theta=1,n_theta
    theta = (i_theta-1)*PI/(n_theta-1)
    grad_kT = grad_kTb * [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)]

    call sample_velocity_dist_unmagnetized(n_samples, u, coulomb_log, kTb, grad_kT, n_b, m_b, q_b, zero3, v_b)

    ! Calculate parallel and perpendicular fluxes
    v_perp1 = 0.d0
    v_perp2 = 0.d0
    v_par = 0.d0
    par1 = grad_kT/norm2(grad_kT)
    call normal_vectors(par1, perp1, perp2)
    do i=1,n_samples
      v_perp1 = v_perp1 + dot_product(v_b(:,i), perp1)
      v_perp2 = v_perp2 + dot_product(v_b(:,i), perp2)
      v_par   = v_par   + dot_product(v_b(:,i), par1)
    end do
      
    ! Check there is no mean flux in the directions perpendicular to grad_kT
    ! 2 is hard-coded for 100k sobol' samples
    ! reasonable compared to the values themselves, which are order few tens km/s
    call assert_equals(0.d0, v_perp1/n_samples, 4.d0, 'perp1 zero')
    call assert_equals(0.d0, v_perp2/n_samples, 4.d0, 'perp2 zero')

    ! TODO Verify the mean flux parallel to grad_kT
  end do
end do
end subroutine test_sampling_flux
end module sample_homma_spec
