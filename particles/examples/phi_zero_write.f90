!> Read a set of particles and static fields and calculate the fieldline trace
!> Command-line arguments: part_restart.h5
program calc_phi_zero
use constants
use mod_particle_sim
use mod_particle_types
use mod_io_actions
use mod_fields
use mod_fields_linear
use mpi
use mod_event
use mod_boris
use mod_fieldline_euler
use mod_find_rz_nearby
!$ use omp_lib
implicit none

type(particle_sim) :: sim_in, sim_out
real*8, parameter :: dt = 1d-8
real*8, parameter :: v = 1d6
type(event) :: fieldreader, partreader
type(event), dimension(:), allocatable :: events
integer :: i, j, k, n_steps, n_lost, n_lost_all, i_elm_old, ifail, ierr
real*8 :: target_time
!$ real*8 :: w0, w1
real*8 :: E(3), B(3), rz_old(2), st_old(2), psi, U, mmm(3)
character(len=500) :: part_file
part_file = ''

! Start up MPI, jorek
call sim_out%initialize(num_groups=0)
call get_command_argument(1, part_file)

partreader = event(read_action(filename=trim(part_file)))
call with(sim_in, partreader) 
fieldreader = event(read_jorek_fields_interp_linear(i=-1))
call with(sim_in, fieldreader)
! For some reason the code will segfault if you use the same event for field and particle reading
! It is also important to read the fields second. Something tricky is going on with the particle reader
! which ruins the fields that have been read because on its own this works fine.
! This is a very bad sign but I have no time to fix it now. Leave the ordering above as-is.

! Create new groups for particles in sim_out
! messy due to all the select-types, refactor?
deallocate(sim_out%groups)
allocate(sim_out%groups(size(sim_in%groups)))
do i=1,size(sim_in%groups)
  allocate(particle_fieldline::sim_out%groups(i)%particles(size(sim_in%groups(i)%particles)))
  select type (particles => sim_in%groups(i)%particles)
  type is (particle_kinetic_leapfrog)
    do j=1,size(particles)
      if (particles(j)%i_elm .le. 0) then
        sim_out%groups(i)%particles(j) = particles(j)
      else
        ! Calculate magnetic field to get GC coordinate
        call sim_in%fields%calc_EBpsiU(0.d0, particles(j)%i_elm, &
            particles(j)%st, particles(j)%x(3), E, B, psi, U)
        sim_out%groups(i)%particles(j) = kinetic_leapfrog_to_gc(sim_in%fields%node_list, sim_in%fields%element_list, particles(j), E, B, sim_in%groups(i)%mass, dt=0.d0)
        ! dt above is not per-se the right dt for this particle (since we read it from a file). Use 0 instead
      end if
    end do
  type is (particle_gc)
    do j=1,size(particles)
      sim_out%groups(i)%particles(j) = particles(j)
    end do
  end select ! Add types here as needed
  select type (p => sim_out%groups(i)%particles)
  type is (particle_fieldline)
    do j=1,size(p)
      p(j)%v = v
    end do
  end select
end do


! Perform a simple equilibrium simulation
events = [event(write_action(basename='phi0_'), step=1d-7)]
!call check_and_fix_timesteps(timesteps, events)
call with(sim_out, events, at=0.d0)
! Trace each fieldline particle until we get to the outboard midplane (theta = 0)
do while (.not. sim_out%stop_now)
  target_time = next_event_at(sim_out, events)

  !$ w0 = omp_get_wtime()
  do i=1,size(sim_out%groups,1)
    n_steps = nint((target_time - sim_out%time)/dt)
    n_lost = 0

    select type (particles => sim_out%groups(i)%particles)
    type is (particle_fieldline)
      !$omp parallel default(private) &
      !$omp shared(sim_in, sim_out, n_steps, i) &
      !$omp reduction(+:n_lost)

      !$omp do
      do j=1,size(particles,1)
        do k=1,n_steps
          if (particles(j)%i_elm .le. 0) exit
          call sim_in%fields%calc_EBpsiU(0.d0, particles(j)%i_elm, &
              particles(j)%st, particles(j)%x(3), E, B, psi, U)
          rz_old    = particles(j)%x(1:2)
          st_old    = particles(j)%st
          i_elm_old = particles(j)%i_elm

          call fieldline_euler_push_cylindrical(particles(j), B, dt)
          call find_RZ_nearby(sim_in%fields%node_list, sim_in%fields%element_list, rz_old(1), rz_old(2), st_old(1), st_old(2), i_elm_old, &
              particles(j)%x(1), particles(j)%x(2), particles(j)%st(1), particles(j)%st(2), particles(j)%i_elm, ifail)
          if (particles(j)%i_elm .le. 0) n_lost = n_lost + 1
        end do ! steps
      end do ! particles
      !$omp end do
      !$omp end parallel
    end select
    call MPI_REDUCE(n_lost, n_lost_all, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    if (sim_out%my_id .eq. 0) write(*,*) i, "number of lost particles: ", n_lost_all
  end do ! groups
  !$ w1 = omp_get_wtime()
  !$ mmm = mpi_minmeanmax(w1-w0)
  !$ if (sim_out%my_id .eq. 0) write(*,"(A,3f9.5,A)") "Particle stepping complete in ", mmm, "s"
  sim_out%time = target_time
  if (sim_out%my_id .eq. 0) write(*,*) "sim%time=", sim_out%time
  call with(sim_out, events, at=sim_out%time)
end do
call sim_out%finalize
end program calc_phi_zero
