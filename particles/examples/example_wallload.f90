!> Assess wall loads using relativistic guiding centers.
!>
!> This script requires that markers are given as an input in a part_restart.h5 file. You will
!> also need the wall input. See the wiki page "How to assess wall loads with the particle tracker?"
!> for details.
!>
Program  example_wallload

use particle_tracer
use mod_particle_io
use mod_particle_diagnostics
use mod_fields_linear   
use mod_fields_hermite_birkhoff
use mod_gc_relativistic
use mod_wall_collision
                  
implicit none

! Set up the simulation variables
real(kind=8)                      :: timesteps(1) = [1.d-10] 
real(kind=8)                      :: target_time, t
integer(kind=4)                   :: n_part, i, j, k, l, n_steps, ifail, max_depth, wall_id
type(write_particle_diagnostics)  :: diag
real(kind=8),dimension(3)         :: pos_prev, wall_pos

type(octree_node) :: wall

max_depth = 6
call mod_wall_collision_init('wall.h5',max_depth,wall)

call sim%initialize(num_groups=1)
call read_simulation_hdf5(sim, 'part_restart.h5')
n_part = size(sim%groups(1)%particles)

! Set up the diagnostics output
!diag = write_particle_diagnostics(filename='diag.h5',only=[1,2,6,12,13,14,15]) ! store total and kinetic energies, p_phi, ielm, phi, R, Z

! Set events to write output data and stop the simulation.
! One can use read_jorek_fields_interp_linear or read_jorek_fields_interp_hermite_birkhoff,
! and i=-1 (to read jorek_restart.h5 and keep this field at all time) or i=last_file_before_time(sim%time)
! (to read a sequel of jorekXXXXX.h5 files and use time-evolving fields)
events = [event(read_jorek_fields_interp_linear(i=last_file_before_time(sim%time))), & 
     !event(diag,start=sim%time,step=1d-8),         &
     event(stop_action(),start=sim%time+7.d-4)]

! Run first event to read the JOREK fields
call with(sim, events, at=0.d0)

do i=1,size(sim%groups(1)%particles)
   select type (p=>sim%groups(1)%particles(i))
   type is (particle_gc_relativistic)
      call find_RZ(sim%fields%node_list, sim%fields%element_list, &
           p%x(1), p%x(2), &
           p%x(1), p%x(2), p%i_elm, p%st(1), p%st(2), ifail)
   end select
end do

call check_and_fix_timesteps(timesteps, events)

do while (.not. sim%stop_now)
  target_time = next_event_at(sim, events)

  do i=1,1
    n_steps = nint((target_time - sim%time)/timesteps(i))

    select type (particles => sim%groups(i)%particles)
    type is (particle_gc_relativistic)	
      !$omp parallel do default(private) &
      !$omp shared (i, n_steps, timesteps, sim, wall)
       do j=1,size(particles,1)
          do k=1,n_steps
             if(particles(j)%i_elm .le. 0) exit

             pos_prev = particles(j)%x

             call runge_kutta_fixed_dt_gc_push_jorek(sim%fields,sim%time,timesteps(i), &
                  sim%groups(i)%mass,particles(j))
             if(particles(j)%i_elm .le. 0) then
                particles(j)%i_elm = 0 ! This marker was lost outside the grid without hitting the 3D triangular mesh wall
                exit
             end if

             call mod_wall_collision_check(pos_prev, particles(j)%x, wall, wall_id, wall_pos)
             if(wall_id .gt. 0) then
                particles(j)%x      = wall_pos
                particles(j)%i_elm  = -wall_id
             end if
          end do
       end do
       !$omp end parallel do
    end select

  enddo

  sim%time = target_time
  call with(sim, events, at=sim%time)
enddo

call mod_wall_collision_free(wall)

call write_simulation_hdf5(sim, 'part_out.h5')

call mod_wall_collision_export(sim, 'wallload.h5')

! Finalize the simulation
call sim%finalize

end program example_wallload
