!> Calculate the particle trajectories with the Boris method in
!> static JOREK fields, file jorek_restart.h5 or jorek_restart.rst
program ex1_jorek
use particle_tracer
implicit none

real*8 :: timesteps(2) = [1d-7,1d-7] !< seconds

type(adf11_all) :: adas
type(coronal)   :: cor
integer :: i, j, k, n_steps, i_elm_old, ifail, n_lost
real*8 :: target_time, t
real*8 :: E(3), B(3), rz_old(2), st_old(2), psi, U
type(read_jorek_fields_interp_linear) :: fieldreader

! Start up MPI, jorek
call sim%initialize(num_groups=2)

! Set up the field reader
fieldreader = read_jorek_fields_interp_linear(basename='jorek_restart', i=-1)
call with(sim, fieldreader)

! Set up particles
allocate(particle_kinetic_leapfrog::sim%groups(1)%particles(100000))
allocate(particle_kinetic_leapfrog::sim%groups(2)%particles(100000))
sim%groups(:)%Z    = 74
sim%groups(:)%mass = 183.84 !< atomic mass units

! Prepare the coronal equilibrium
adas = read_adf11(0,'50_w')
cor  = coronal(adas)

! Distribute particles uniformly throughout the domain
call initialise_particles(sim%groups(1)%particles, &
    sim%fields%node_list, sim%fields%element_list, sobseq_rng())
call set_velocity_from_T(sim%groups(1)%particles, sim%groups(1)%mass, &
    sim%fields%node_list, sim%fields%element_list, cor, v_par=.false.)
call initialise_particles(sim%groups(2)%particles, &
    sim%fields%node_list, sim%fields%element_list, pcg32_rng())
call set_velocity_from_T(sim%groups(2)%particles, sim%groups(2)%mass, &
    sim%fields%node_list, sim%fields%element_list, cor, v_par=.false.)
call adjust_particle_weights(sim%groups(1)%particles, num_atoms_total=1d23)
call adjust_particle_weights(sim%groups(2)%particles, num_atoms_total=1d23)

events = [event(write_action(basename='test'),   step=1d-4), &
          !event(diag_print_kinetic_energy(),     step=1d-6), &
          event(projection(sim%fields%node_list, sim%fields%element_list, f=[proj_f(proj_one, 1)], filter=1d-3, basename='proj', to_vtk=.true.), step=1d-5), &
          event(stop_action(), start=1d-3)]
call check_and_fix_timesteps(timesteps, events)
call with(sim, events, at=0.d0)

do while (.not. sim%stop_now)
  target_time = next_event_at(sim, events)

  do i=1,size(sim%groups,1)
    n_steps = nint((target_time - sim%time)/timesteps(i))
    n_lost = 0

    select type (particles => sim%groups(i)%particles)
    type is (particle_kinetic_leapfrog)
      !$omp parallel do default(private) &
      !$omp shared(sim, n_steps, timesteps, i) &
      !$omp reduction(+:n_lost)
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
          if (particles(j)%i_elm .le. 0) n_lost = n_lost + 1
        end do ! steps
      end do ! particles
      !$omp end parallel do
    end select
    write(*,*) "number of lost particles: ", n_lost
  end do ! groups
  sim%time = target_time
  call with(sim, events, at=sim%time)
end do
call sim%finalize
end program ex1_jorek
