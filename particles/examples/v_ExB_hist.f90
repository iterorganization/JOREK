!> Calculate the radial ExB drift velocity at many locations uniformly throughout the plasma
!> Find the closest file to the time given on the commandline and output a histogram ordered by normalized psi.
!>
!> The histogram is calculated for radial bins of size 0.01.
!> Velocity information is histogrammed from +- 6000 psi_N/s in 2000 bins
!>
!> We need to read an input file to obtain F0 for the magnetic field.
!>
!> Not really suitable for all types of cases, optimized for single xpoint case
program v_ExB_hist
use mod_particle_sim
use mod_fields_linear, only: last_file_before_time
use mod_fields_hermite_birkhoff
use phys_module, only: xpoint, xcase
use mod_sobseq_rng, only: sobseq_rng
use mod_random_seed, only: random_seed
use constants, only: TWOPI
use mod_initialise_particles, only: domain_bounding_box
use mod_math_operators, only: cross_product
use domains
!$ use omp_lib
implicit none

integer, parameter :: n_points = 100000000, nx = 440, nv = 1000
real*8, parameter :: x_min = 0.d0, x_max = 1.1d0, dx = (x_max-x_min)/real(nx,8)
real*8, parameter :: v_min = -6d3, v_max = 6d3, dv = (v_max-v_min)/real(nv,8)

integer :: counts(nv,nx)

type(particle_sim) :: sim
type(sobseq_rng) :: rng
type(read_jorek_fields_interp_hermite_birkhoff) :: fieldreader
integer :: i_elm, ifail, i, seed, n_threads, i_thread, unit
integer :: i_v, i_x
real*8 :: s, t, x(3), psi_N, R, Z, grad_Psi(2)
real*8 :: Rmin, Rmax, Zmin, Zmax, Z_xpoint
real*8 :: psi_axis, psi_xpoint(2)
real*8 :: E(3), B(3), psi, U, v(3)
character(len=20) :: time_s

! Start up MPI, jorek
call sim%initialize(num_groups=0)
call get_command_argument(1, time_s)
read(time_s,*) sim%time

fieldreader = read_jorek_fields_interp_hermite_birkhoff(i=last_file_before_time(sim%time))
call fieldreader%do(sim)

! Calculate normalization coefficients psi_axis and psi_xpoint(1)
call find_axis(0, sim%fields%node_list, sim%fields%element_list, psi_axis, x(1), x(2), i_elm, s, t, ifail)
call find_xpoint(0, sim%fields%node_list, sim%fields%element_list, psi_xpoint, x(1), Z_xpoint, i_elm, s, t, xcase, ifail)
write(*,*) "Psi_xpoint: ", psi_xpoint(1), " Psi_axis: ", psi_axis
! Get the size of the domain
call domain_bounding_box(sim%fields%node_list, sim%fields%element_list, Rmin, Rmax, Zmin, Zmax)

seed = random_seed()
n_threads = 1
i_thread = 1
counts = 0
!$omp parallel default(private) shared(sim, seed, Rmin, Rmax, Zmin, Zmax, psi_xpoint, psi_axis, Z_xpoint) reduction(+:counts)
!$ n_threads = omp_get_num_threads()
!$ i_thread = omp_get_thread_num()+1
call rng%initialize(3, seed, n_threads, i_thread)
!$omp do
do i=1,n_points
  ! Generate a random position
  call rng%next(x)
  ! Find corresponding s, t and i_elm
  call find_RZ(sim%fields%node_list, sim%fields%element_list, &
      Rmin + x(1)*(Rmax - Rmin), Zmin + x(2)*(Zmax - Zmin), &
      R, Z, i_elm, s, t, ifail)
  if (ifail .ne. 0 .or. i_elm .le. 0) cycle ! Particle out of domain, try a new position in the next loop

  ! Calculate the region of the plasma the particle is in
  if (Z .lt. Z_xpoint) cycle ! skip positions below the xpoint

  ! Calculate E and B fields
  call sim%fields%calc_EBpsiU(sim%time, i_elm, [s, t], TWOPI*x(3), E, B, psi, U)

  v = cross_product(E, B)/dot_product(B,B) ! in RZPhi coordinates, m/s

  ! Calculate the normal vector to the flux surface, Grad Psi / |Grad Psi|
  ! Assume no variation of Psi in the toroidal direction
  ! 
  grad_psi = -([-B(2),B(1)])*R/(psi_xpoint(1) - psi_axis) ! in poloidal plane, normalize to grad_psi_N
  ! Convert v_r from m/s to psi/s by multiplying by |grad_psi|
  i_v = ceiling((dot_product(v(1:2), grad_psi) - v_min)/dv)
  psi_N = (psi - psi_axis)/(psi_xpoint(1) - psi_axis)
  i_x = ceiling((psi_N - x_min)/dx)
  if (i_v .ge. lbound(counts,1) .and. i_x .ge. lbound(counts,2) .and. &
      i_v .le. ubound(counts,1) .and. i_x .le. ubound(counts,2)) then
    counts(i_v, i_x) = counts(i_v, i_x) + 1
  end if
end do
!$omp end do
!$omp end parallel

! Normalize counts per layer

! Write out the counts
open(newunit=unit, file=trim(time_s)//'_v_ExB.dat', status='replace', access='stream', action='write')
write(unit) counts
close(unit)
write(*,*) 'output written to ', trim(time_s), '_v_ExB.dat'

call sim%finalize
end program v_ExB_hist
