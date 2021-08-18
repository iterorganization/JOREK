!> particles/benchmarks/collisions/diamagnetic_thermal.f90
!> Simulate the temperature screening force on test particles on timescales
!> longer than the gyromotion.
!> This corresponds to case II from Homma JCP 2013.
!>
!> Use a magnetic field of 1 T and 10^4 test particles.
!>
!> We do not use the full simulation capabilities here, we just fake the temperature gradient,
!> plasma density and most other variables.
!>
!> Parameters for this case
!> Initial position: 0
!> Initial velocity: sampled from 50 eV maxwell distribution
!> Background ion species: H+
!> Background ion density: 10^20 m^-3
!> Background flow velocity: 0
!> Temperature at origin: 50 eV
!> Coulomb logarithm: 15
!> Particle: W^{3+}
!>
!> We cheat a bit here and store the xyz position in particle%x instead of the rzphi
!> and use the cartesian boris method, without using JOREK elements.
program temperature_screening
use particle_tracer
use constants
use mod_pcg32_rng
use mod_random_seed
use mod_collisions
use mod_sampling
!$ use omp_lib
implicit none

real*8, parameter :: timestep = 2.03d-10 ! Timestep for BCM, use same step for boris method
real*8, parameter :: omega_A = 3d0*EL_CHG * 1d0 / (183.84 * ATOMIC_MASS_UNIT)
real*8, parameter :: t_end = 650d0/omega_A
real*8, parameter :: dt_out = t_end/40.d0
real*8, parameter :: E(3) = [0.d0, 0.d0, 0.d0]
real*8, parameter :: B(3) = [0.d0, 0.d0, 1d0]

real*8, parameter :: kTb = 50.d0*EL_CHG! [J]
!real*8, parameter :: grad_kTb(3) = [100.d0, 0.d0, 0.d0]*EL_CHG !< [J/m] (II-1)
real*8, parameter :: grad_kTb(3) = [300.d0, 0.d0, 0.d0]*EL_CHG !< [J/m] (II-2)

integer, parameter :: n_particles = 10000

real*8, parameter :: n_b = 1d20 !< [m^-3]
real*8, parameter :: m_b = 1.d0 !< [u]
integer*1, parameter :: q_b = 1 !< [e]
real*8, parameter :: v_bg(3) = [0.d0, 0.d0, 0.d0] !< [m^-3]

integer*1, parameter :: q_a = 3 !< [e]
real*8, parameter :: m_a = 183.84 !< atomic mass units
real*8, parameter :: V_thermal = sqrt(50.d0 * EL_CHG / (m_a * ATOMIC_MASS_UNIT))


type(particle_kinetic_leapfrog), dimension(:), allocatable :: particles

real*8 :: u(6), coulomb_log, x_av(3), t, v_tmp(4)

integer :: i, j, k, seed, n_threads, i_thread, unit
type(pcg32_rng) :: rng
real*8 :: v_b(3) !< Sampled background particle velocity
real*8 :: q(3) !< Heat flux vector [W m^-2]

n_threads = 1
i_thread = 1
seed = random_seed()
! RNG for sampling initial velocities
call rng%initialize(6, seed, n_threads, i_thread)

! Set up W particles
do i=1,1
  allocate(particles(n_particles))
  particles(:)%q = q_a
  do j=1,size(particles,1)
    particles(j)%x = [0.d0, 0.d0, 0.d0]
    call rng%next(u)
    v_tmp = boxmueller_transform(u(1:4))*V_thermal ! 2 gaussian distributed random numbers
    particles(j)%v = v_tmp(1:3)
  end do
end do

coulomb_log = coulomb_logarithm(kTb, n_b, 3_1, q_b, m_a, m_b)
write(*,*) 'Coulomb log background ion-W: ', coulomb_log
write(*,*) 'Coulomb log background ion-ion: ', coulomb_logarithm(kTb, n_b, q_b, q_b, m_b, m_b)
write(*,*) 'Forcing coulomb log to 15 everywhere'
coulomb_log = 15.d0

seed = random_seed()
!seed = -1904078196 ! crashes after ~10 minutes
seed = -367353379 ! crashes quickly with R(5) = 0
write(*,*) seed ! store for reproducibility


call rm('II.txt')
open(newunit=unit, file='II.txt', status='new')

! We need private here for gfortran which does not handle derived types correctly
!$omp parallel default(private) shared(seed, unit, particles, coulomb_log) &
!$omp private(n_threads, i_thread, rng, i, j, k, q, u, x_av, v_b, t)
n_threads = 1
i_thread = 1
!$ n_threads = omp_get_num_threads()
!$ i_thread = omp_get_thread_num()
call rng%initialize(6, seed, n_threads, i_thread)

q = q_homma2013(kTb, grad_kTb, B, n_b, m_b, q_b)

do i=1,nint(t_end/dt_out) ! number of output steps (roughly, not exactly at dt_out)
  !$omp do
  do j=1,size(particles,1)
    do k=1,nint(dt_out/timestep) ! number of timesteps * timestep size describe real time
      ! Calculate collision
      call rng%next(u)
      call sample_velocity_dist_magnetized(1, u, kTb, q, n_b, m_b, q_b, v_bg, v_b)

      call rng%next(u)
      call collide_particles(u(1:3), particles(j)%q, m_a, particles(j)%v, &
          q_b, m_b, v_b, n_b, coulomb_log, timestep)
      call boris_push_cartesian(particles(j), m_a, E, B, timestep)
    end do ! steps
  end do ! particles
  !$omp end do

  !$omp barrier
  ! Calculate average particle velocity in y direction
  !$omp single
  write(6,'(A)',advance='no') '.'
  call flush(6)
  ! calculate average velocity
  do j=1,3
    x_av(j) = sum(particles(:)%x(j))/real(size(particles,1),8)
  end do
  t = timestep * real(nint(dt_out/timestep)*i, 8)
  write(unit,'(4g16.7)') omega_a*t, x_av
  !$omp end single
  !!$omp barrier ! is implied by omp single
end do ! output times
!$omp end parallel

close(unit)
contains
subroutine rm(file)
  character(len=*), intent(in) :: file
  integer :: u, stat
  open(newunit=u, iostat=stat, file=file, status='old')
  if (stat .eq. 0) close(u, status='delete')
end subroutine rm
end program temperature_screening
