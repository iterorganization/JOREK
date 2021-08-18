!> particles/benchmarks/collisions/diamagnetic_thermal.f90
!> Simulate the diamagnetic thermal force on test particles on very short timescale.
!> This corresponds to case I from Homma JCP 2013.
!>
!> Use a weak magnetic field of 0.1 T and 10^6 test particles for a very short time.
!> Even on a scale of less than a gyroperiod we can discern the thermal force
!> given enough particles to get a statistically accurate result.
!>
!> We do not use the full simulation capabilities here, we just fake the temperature gradient,
!> plasma density and most other variables.
!>
!> Parameters for this case
!> Initial position: 0
!> Initial velocity: 8.84e3 m/s e_x (50 eV)
!> Normalized velocity: 0.09 e_x
!> Background ion species: H+
!> Background ion density: 10^20 m^-3
!> Background flow velocity: 0
!> Temperature at origin: 50 eV
!> Coulomb logarithm: 15 (we compare our place-dependent value against this)
!> Particle: W^{3+}
!>
!> We cheat a bit here and store the xyz position in particle%x instead of the rzphi
!> and use the cartesian boris method, without using JOREK elements.
program diamagnetic_thermal
use particle_tracer
use constants
use mod_pcg32_rng
use mod_random_seed
use mod_collisions
!$ use omp_lib
implicit none

real*8, parameter :: timestep = 2.03d-10 ! Timestep for BCM, use same step for boris method
real*8, parameter :: omega_A = 3d0*EL_CHG * 0.1d0 / (183.84 * ATOMIC_MASS_UNIT)
real*8, parameter :: t_end = 0.01d0/omega_A
real*8, parameter :: dt_out = t_end/40.d0
real*8, parameter :: E(3) = [0.d0, 0.d0, 0.d0]
real*8, parameter :: B(3) = [0.d0, 0.d0, 0.1d0]

real*8, parameter :: kTb = 50.d0*EL_CHG ! [J]
real*8, parameter :: grad_kTb(3) = [100.d0, 0.d0, 0.d0]*EL_CHG !< [J/m] (I-1)
!real*8, parameter :: grad_kTb(3) = [300.d0, 0.d0, 0.d0]*EL_CHG !< [J/m] (I-2)
!real*8, parameter :: grad_kTb(3) = [300.d0, 0.d0, 5.d0]*EL_CHG !< [J/m] (I-3)

integer, parameter :: n_particles = 1000000

real*8, parameter :: n_b = 1d20 !< [m^-3]
real*8, parameter :: m_b = 1.d0 !< [u]
integer*1, parameter :: q_b = 1 !< [e]
real*8, parameter :: v_bg(3) = [0.d0, 0.d0, 0.d0] !< [m^-3]

integer*1, parameter :: q_a = 3 !< [e]
real*8, parameter :: v_a = 8.849e3 !< [m/s]
real*8, parameter :: m_a = 183.84 !< atomic mass units
type(particle_kinetic_leapfrog), dimension(:), allocatable :: particles

real*8 :: u(6), coulomb_log, v_av(3), t

integer :: i, j, k, seed, n_threads, i_thread, unit
type(pcg32_rng) :: rng
real*8 :: v_b(3) !< Sampled background particle velocity
real*8 :: q(3) !< Heat flux vector [W m^-2]

! Set up W particles
do i=1,1
  allocate(particles(n_particles))
  particles(:)%q = q_a
  do j=1,size(particles,1)
    particles(j)%x = [0.d0, 0.d0, 0.d0]
    particles(j)%v = [v_a, 0.d0, 0.d0]
  end do
end do

coulomb_log = coulomb_logarithm(kTb, n_b, 3_1, q_b, m_a, m_b)
write(*,*) 'Coulomb log background ion-W: ', coulomb_log
write(*,*) 'Coulomb log background ion-ion: ', coulomb_logarithm(kTb, n_b, q_b, q_b, m_b, m_b)
write(*,*) 'Forcing coulomb log to 15 everywhere'
coulomb_log = 15.d0

seed = random_seed()
!seed = -1289433091 ! this one failed twice, fixed now
write(*,*) seed ! store for reproducibility


call rm('I.txt')
open(newunit=unit, file='I.txt', status='new')

! We need private here for gfortran which does not handle derived types correctly
!$omp parallel default(private) shared(seed, unit, particles, coulomb_log) &
!$omp private(n_threads, i_thread, rng, i, j, k, q, u, v_av, v_b, t)
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
    v_av(j) = sum(particles(:)%v(j))/real(size(particles,1),8)
  end do
  t = timestep * real(nint(dt_out/timestep)*i, 8)
  write(unit,'(4g16.7)') omega_a*t, v_av - v_a * [cos(omega_a * t), -sin(omega_a * t), 0.d0]
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
end program diamagnetic_thermal
