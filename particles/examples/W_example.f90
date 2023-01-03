!> particles/examples/W_example.f90
!> Requires: acd50_w.dat, scd50_w.dat and the JOREK input file
!> Compile with `j2p W_example` and run with `W_example < in1`
!> Read fields from jorek00721.h5 and onwards and push W particles
program W_example
use particle_tracer
use mod_sobseq_rng
use mod_particle_diagnostics
implicit none

type(adf11_all) :: adas
type(coronal)   :: cor
type(write_particle_diagnostics) :: diag

real*8 :: timesteps(1) = [2d-9]
integer :: i, j, k, n_steps, n_lost, i_elm_old, ifail
real*8 :: target_time, t
real*8 :: E(3), B(3), rz_old(2), st_old(2), psi, U
real*8, parameter :: alpha = 0.d0 !< if 0, sample from maxwellian
integer :: ierr, my_id
type(event) :: fieldreader

! Start up MPI, jorek
call sim%initialize(num_groups=1)

! Set up the field reader
fieldreader = event(read_jorek_fields_interp_hermite_birkhoff(&
    basename='jorek', i=721, stop_at_end=.false.))
call with(sim, fieldreader)

! Prepare the coronal equilibrium
call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)
adas = read_adf11(my_id, '50_w')
cor  = coronal(adas)

! Set up particles
sim%groups(:)%Z    = 74
sim%groups(:)%mass = 183.84 !< atomic mass units
do i=1,1
  allocate(particle_kinetic_leapfrog::sim%groups(i)%particles(10000))
  ! For every particle accept or reject it with probability f_psi_inside(psi)
  call initialise_particles_H_mu_psi(sim%groups(i)%particles, &
      sim%fields, sobseq_rng(), &
      sim%groups(i)%mass, &
      uniform_space=.true., uniform_space_rej_f=f_psi_inside, &
      uniform_space_rej_vars=[1], cor=cor)

  select type (p => sim%groups(i)%particles)
  type is (particle_kinetic_leapfrog)
    p(:)%weight = 1.d0
    call boris_all_initial_half_step_backwards_RZPhi(p, sim%groups(i)%mass, &
        sim%fields, sim%time, timesteps(i))
  end select
end do

! Set up the diagnostics output
diag = write_particle_diagnostics(filename='diag.h5')

! Perform a simple equilibrium simulation
events = [fieldreader, &
          event(write_action(),   step=5d-6), &
          event(projection(sim%fields%node_list, sim%fields%element_list, &
            filter=1d-4, to_h5=.true.), step=5d-6), &
          event(diag, step=5d-6), &
          event(stop_action(), start=5d-3)]
call check_and_fix_timesteps(timesteps, events)
call with(sim, events, at=0.d0)

do while (.not. sim%stop_now)
  target_time = next_event_at(sim, events)
  do i=1,size(sim%groups,1)
    n_steps = nint((target_time - sim%time)/timesteps(i))
    select type (particles => sim%groups(i)%particles)
    type is (particle_kinetic_leapfrog)
      !$omp parallel do default(private) &
      !$omp shared(sim, n_steps, timesteps, i)
      do j=1,size(particles,1)
        do k=1,n_steps
          if (particles(j)%i_elm .le. 0) exit
          t = sim%time + k*timesteps(i)
          call sim%fields%calc_EBpsiU(t, particles(j)%i_elm, &
              particles(j)%st, particles(j)%x(3), E, B, psi, U)
          rz_old    = particles(j)%x(1:2)
          st_old    = particles(j)%st
          i_elm_old = particles(j)%i_elm

          call boris_push_cylindrical(particles(j), sim%groups(i)%mass, E, B, timesteps(i))
          call find_RZ_nearby(sim%fields%node_list, sim%fields%element_list, rz_old(1), rz_old(2), st_old(1), st_old(2), i_elm_old, &
              particles(j)%x(1), particles(j)%x(2), particles(j)%st(1), particles(j)%st(2), particles(j)%i_elm, ifail)
        end do ! steps
      end do ! particles
      !$omp end parallel do
    end select
  end do ! groups
  sim%time = target_time
  call with(sim, events, at=sim%time)
end do
call sim%finalize
contains
pure function f_psi_inside(n, P, grad_P) result(f)
  integer, intent(in) :: n
  real*8, intent(in) :: P(n), grad_P(3,n)
  real*4 :: f
  f = 5e-2
  if (P(1) .lt. -0.26 .and. P(1) .gt. -0.35) f = 1e0
end function f_psi_inside
end program W_example
