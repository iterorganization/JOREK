!> Module containing test routines for collisions with the BCM (Homma, Hatayama JCP 2012,2013)
module mod_test_collisions
  implicit none
contains
subroutine run_collision_test(testcase, n_particles, dt, n_step, n_out, v_out, rng)
  use mod_collision_case
  use mod_particle_types
  use mod_boris
  use mod_rng
  use mod_pcg32_rng
  use mod_random_seed
  use mod_collisions
  !$ use omp_lib
  type(collision_case), intent(in)  :: testcase
  integer, intent(in)               :: n_particles
  real*8, intent(in)                :: dt
  integer, intent(in)               :: n_step ! number of total steps
  integer, intent(in)               :: n_out ! number of outputs
  real*8, dimension(3,n_out), intent(out) :: v_out
  class(type_rng), intent(in), optional :: rng

  class(type_rng), allocatable :: my_rng
  integer :: ifail
  class(particle_kinetic_leapfrog), allocatable, dimension(:) :: particles
  integer :: i, j, k, m, i_thread, n_threads, seed
  real*8 :: u(6), v_p(3)

  if (present(rng)) then
    allocate(my_rng, source=rng)
  else
    allocate(my_rng, source=pcg32_rng())
  end if

  write(*,*) "Running collision test with n=", n_particles, "particles"
  write(*,*) "dt =", dt, "n_step=",n_step
  write(*,*) "v_0=", testcase%v_0
  write(*,*) "q_a=", testcase%q_a, "m_a=", testcase%m_a
  write(*,*) "q_b=", testcase%q_b, "m_b=", testcase%m_b
  write(*,*) "v_b=", [0.d0,0.d0,0.d0]
  write(*,*) "kT =", testcase%kT/EL_CHG
  write(*,*) "Grad kT=", testcase%grad_kT/EL_CHG
  write(*,*) "n_b=", testcase%n_b
  write(*,*) "E=", testcase%E
  write(*,*) "B=", testcase%B
  write(*,*) "tau_s=", testcase%slowdown_time
  write(*,*) "coulomb logarithm=", testcase%coulomb_log

  ! Preparation
  allocate(particles(n_particles))
  do m=1,n_particles
    particles(m)%v = testcase%v_0
    particles(m)%x = testcase%x_0
    particles(m)%q = testcase%q_a
  end do
  seed=random_seed()
  !$omp parallel default(none) shared(particles, testcase, dt, seed, n_step, n_out, n_particles, v_out) &
  !$omp private(j,k,m,i_thread,n_threads,ifail,u,v_p) firstprivate(my_rng)
  i_thread=1
  n_threads=1
  !$ i_thread = omp_get_thread_num()
  !$ n_threads = omp_get_num_threads()
  call my_rng%initialize(6, seed, n_threads, i_thread, ifail)
  do i=1,n_out
    !$omp do
    do m=1,n_particles
      do j=1,n_step/n_out
        call my_rng%next(u)
        call sample_velocity_dist_unmagnetized(1, u, testcase%coulomb_log, &
            testcase%kT, testcase%grad_kT, testcase%n_b, testcase%m_b, testcase%q_b, v_b=[0.d0,0.d0,0.d0], w_1=v_p)
        call my_rng%next(u)
        call collide_particles(u(1:3), testcase%q_a, testcase%m_a, particles(m)%v, &
            testcase%q_b, testcase%m_b, v_p, testcase%n_b, testcase%coulomb_log, dt)

        call boris_push_cartesian(particles(m), testcase%m_a, testcase%E, testcase%B, dt)
      end do
    end do
    !$omp end do
    !$omp barrier
    !$omp master
    v_out(:,i) = [(sum(particles(:)%v(k)),k=1,3)]
    write(*,'(A)',advance='no') '.'
    !$omp end master
    !$omp barrier
  end do
  !$omp end parallel
  v_out = v_out/(norm2(testcase%v_0)*real(n_particles))
  write(*,*) "Done"
end subroutine run_collision_test

end module mod_test_collisions
