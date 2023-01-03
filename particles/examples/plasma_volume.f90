!> Calculate the plasma volume by sampling many particles uniformly
!> Find the closest file to the time given on the commandline and output psi values of particles
!>
!> We need to read an input file to obtain F0 for the magnetic field.
!>
!> Not really suitable for all types of cases, optimized for single xpoint case
program plasma_volume
use mod_particle_sim
use mod_fields_linear, only: last_file_before_time
use mod_fields_hermite_birkhoff
use phys_module, only: xpoint, xcase
use mod_random_seed, only: random_seed
use mod_sobseq_rng, only: sobseq_rng
use constants, only: TWOPI, PI
use mod_initialise_particles, only: domain_bounding_box
use domains
use hdf5_io_module
use mod_sampling, only: transform_uniform_cylindrical
!$ use omp_lib
implicit none

integer, parameter :: n_points = 10000000
logical, parameter :: above_xpoint_only = .true.

type(particle_sim) :: sim
type(sobseq_rng) :: rng
type(read_jorek_fields_interp_hermite_birkhoff) :: fieldreader
integer :: i_elm, ifail, i, seed, n_threads, i_thread, unit
integer :: i_v, i_x
real*8 :: s, t, x(3), psi_N, R, Z, Phi, grad_Psi(2)
real*8 :: Rmin, Rmax, Zmin, Zmax, Z_xpoint
real*8 :: psi_axis, psi_xpoint(2)
real*8 :: E(3), B(3), psi, U, v(3)
character(len=20) :: time_s

integer(HID_T) :: file_id

real*4, allocatable, dimension(:) :: psi_norm

! Start up MPI, jorek
call sim%initialize(num_groups=0)
call get_command_argument(1, time_s)
read(time_s,*) sim%time

fieldreader = read_jorek_fields_interp_hermite_birkhoff(i=last_file_before_time(sim%time))
call fieldreader%do(sim)

! Calculate normalization coefficients psi_axis and psi_xpoint(1)
call find_axis(0, sim%fields%node_list, sim%fields%element_list, psi_axis, x(1), x(2), i_elm, s, t, ifail)
call find_xpoint(0, sim%fields%node_list, sim%fields%element_list, psi_xpoint, x(1), Z_xpoint, i_elm, s, t, xcase, ifail)
! Get the size of the domain
call domain_bounding_box(sim%fields%node_list, sim%fields%element_list, Rmin, Rmax, Zmin, Zmax)
write(*,*) "Box size: R=[", Rmin, ",", Rmax, "], Z=[", Zmin, ",", Zmax, "]"
! Calculate the density from the surface of an annulus (0 - 2pi) times the height
write(*,*) "Point density: ", real(n_points)/(PI*(Rmax**2 - Rmin**2)*(Zmax - Zmin)), " [m^-3]"

allocate(psi_norm(n_points))

seed = random_seed()
n_threads = 1
i_thread = 1
!$omp parallel default(private) shared(sim, seed, Rmin, Rmax, Zmin, Zmax, psi_xpoint, psi_axis, Z_xpoint, psi_norm)
!$ n_threads = omp_get_num_threads()
!$ i_thread = omp_get_thread_num()+1
call rng%initialize(3, seed, n_threads, i_thread)
!$omp do
do i=1,n_points
  ! Generate a random position
  call rng%next(x)
  call transform_uniform_cylindrical(x, [Rmin, Rmax], [Zmin, Zmax], [0.d0,TWOPI], R, Z, Phi)
  ! Find corresponding s, t and i_elm
  call find_RZ(sim%fields%node_list, sim%fields%element_list, &
      R, Z, & ! inputs
      R, Z, i_elm, s, t, ifail) ! outputs
  if (ifail .ne. 0 .or. i_elm .le. 0) cycle ! Particle out of domain, try a new position in the next loop

  ! Calculate the region of the plasma the particle is in
  if (Z .lt. Z_xpoint .and. above_xpoint_only) cycle ! skip positions below the xpoint

  ! Calculate E and B fields
  call sim%fields%calc_EBpsiU(sim%time, i_elm, [s, t], phi, E, B, psi, U)

  psi_norm(i) = real((psi - psi_axis)/(psi_xpoint(1) - psi_axis), 4)
end do
!$omp end do
!$omp end parallel

call HDF5_create('psi.h5', file_id)
call HDF5_array1D_saving_r4(file_id, psi_norm, n_points, '/psi_norm')
call HDF5_close(file_id)

call sim%finalize
end program plasma_volume
