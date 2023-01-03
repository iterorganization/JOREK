!> Push relativistic guiding centre(s) with error control in JOREK fields
!>
!> Compile with `make ex8_jorek`
!> Run with `./ex8_jorek < JOREK_namelist`

program ex8_jorek

use particle_tracer
use mod_particle_io
use mod_particle_diagnostics
use mod_fields_linear   
use mod_fields_hermite_birkhoff
use mod_gc_relativistic
                  
implicit none

! Set up the simulation variables
logical,parameter                   :: write_timestep=.true.
!> error control runge kutta error tolerances: 1:R, 2:Z, 3:phi, 4:p_parallel
real(kind=8),dimension(4),parameter :: tolerances = [1.d-2,1.d-2,1.d-2,1.d3]
real(kind=8)                        :: timesteps(1) = [3.5723d-3]
real(kind=8)                        :: target_time, t
real(kind=8)                        :: energy !< Initial energy in eV (including rest energy)
real(kind=8)                        :: ksi    !< Initial cosine of pitch-angle
real(kind=8)                        :: time_local,dt_local !< local time and time step
integer(kind=4)                     :: n_part, i, j, k, n_steps, ifail, n_lost
logical                             :: restart
type(diag_print_kinetic_energy)     :: print_kinetic_energy
type(write_particle_diagnostics)    :: diag
type(particle_gc_relativistic)      :: particle_in, particle_out
type(particle_gc)                   :: particle_out_gc
real(kind=8)                        :: psi, U, gyro_angle,error,dt_try
real(kind=8),dimension(3)           :: E, B

call sim%initialize(num_groups=1)

! Set up the diagnostics output
diag = write_particle_diagnostics(filename='diag.h5',only=[1,2,6,12,13,14,15]) ! store total and kinetic energies, p_phi, ielm, phi, R, Z

restart = .false. !

if (.not. restart) then
  sim%time = 2.5d-3 ! 1.d-7 !   ! start time 

  ! Allocate a group and particle(s) of type particle_gc_relativistic
  n_part = 1
  allocate(particle_gc_relativistic::sim%groups(1)%particles(n_part))
  sim%groups(1)%mass = 5.4857990907016d-4 !< particle mass in AMU

  ! Set events to write output data and stop the simulation.
  ! One can use read_jorek_fields_interp_linear or read_jorek_fields_interp_hermite_birkhoff,
  ! and i=-1 (to read jorek_restart.h5 and keep this field at all time) or i=last_file_before_time(sim%time)
  ! (to read a sequel of jorekXXXXX.h5 files and use time-evolving fields)
  !events = [event(read_jorek_fields_interp_linear(i=last_file_before_time(sim%time))), & 
  !          event(diag,start=sim%time,step=1d-7),         &
  !          event(stop_action(),start=sim%time+3.d-6)]

  events = [event(read_jorek_fields_interp_linear(i=-1)),&
       event(diag,start=sim%time,step=1d-7),&
       event(stop_action(),start=sim%time+5d-6)]

  ! Run first event to read the JOREK fields
  call with(sim, events, at=0.d0)

  select type (p=>sim%groups(1)%particles(1))
  type is (particle_gc_relativistic)
    p%q = -1
    p%x = [3.1d0,0.d0,0.d0]
    call find_RZ(sim%fields%node_list, sim%fields%element_list, &
                 p%x(1), p%x(2), & ! inputs
                 p%x(1), p%x(2), p%i_elm, p%st(1), p%st(2), ifail) ! outputs
    !p%p = [1.d7,0.]
    energy = 5.12d5 ! 1.d7 ! Particle energy, including rest energy
    ksi    = cos(0.785398) ! Cosine of pitch-angle
    particle_in = p
    particle_out = relativistic_gc_momenta_from_E_cospitch(particle_in,energy,ksi,sim%groups(1)%mass,sim%fields,sim%time)
    p%p = particle_out%p
    write(*,*) 'Initial p(1), p(2)', p%p(1), p%p(2)
  end select

else

  write(*,*) '*******************************'
  write(*,*) 'Reading particle restart file'
  call read_simulation_hdf5(sim, 'part_restart.h5')
  write(*,*) 'Time = ', sim%time
  write(*,*) '*******************************' 
  ! Set events to write output data and stop the simulation.
  ! One can use read_jorek_fields_interp_linear or read_jorek_fields_interp_hermite_birkhoff,
  ! and i=-1 (to read jorek_restart.h5 and keep this field at all time) or i=last_file_before_time(sim%time)
  ! (to read a sequel of jorekXXXXX.h5 files and use time-evolving fields)
  events = [event(read_jorek_fields_interp_linear(i=-1)), & 
            event(diag,start=sim%time,step=1d-8),         &
	    event(stop_action(),start=sim%time+1.d-8)]

  ! Run first event to read the JOREK fields
  call with(sim, events, at=0.d0)   

endif

! Check all events conform to the requested timestep
call check_and_fix_timesteps(timesteps, events)

! Set dpsi/dt=0 (useful e.g. to check the conservation of the total particle energy) 
call sim%fields%set_flag_dpsidt(.true.)

! Open a file where to write some fields at a given position to test time interpolation routines
open(22,file='timestep_vs_time.dat')

! Loop until the simulation is stopped
do while (.not. sim%stop_now)
  ! Extract the next event time
  target_time = next_event_at(sim, events)
  ! Loop over all particle groups
  do i=1,1
    n_lost = 0

    select type (particles => sim%groups(i)%particles)
    type is (particle_gc_relativistic)	
!      !$omp parallel do default(private) &
!      !$omp shared (i, n_steps, timesteps, sim) &
!      !$omp reduction(+:n_lost)	
       do j=1,size(particles,1)
          time_local = sim%time !< copy simulation time in local variable
          dt_local = timesteps(i)!< time step in local variable
          do while((time_local.lt.target_time) .and. (particles(j)%i_elm.gt.0)) !< continue until we reach the target time
#ifdef TEST
          call runge_kutta_error_control_dt_gc_push_jorek(sim%fields,tolerances,&
               time_local,dt_local,target_time,sim%groups(i)%mass,dt_try,particles(j),error) !< push in jorek 
#else   
          call runge_kutta_error_control_dt_gc_push_jorek(sim%fields,tolerances,&
               time_local,dt_local,target_time,sim%groups(i)%mass,dt_try,particles(j)) !< push in jorek fields
#endif
          time_local = time_local + dt_local !< time update
          !> write time step profile if enables
#ifdef TEST
          if(write_timestep) write(22,'(i6,4e26.16)') j,time_local,dt_local,dt_try,error
#else
          if(write_timestep) write(22,'(i6,3e26.16)') j,time_local,dt_local,dt_try
#endif
          dt_local = dt_try !< update the time step

          if (particles(j)%i_elm .le. 0) n_lost = n_lost + 1		
       end do !< time steps
      end do !< particles
!      !$omp end parallel do
    end select
    write(*,*) "number of lost particles: ", n_lost	  
  enddo !< groups

  ! Update current time and run events
  sim%time = target_time
  call with(sim, events, at=sim%time)
enddo !< event

! Print particle information
write(*,*) 'Final R, Z, phi: ', sim%groups(1)%particles(1)%x(1), sim%groups(1)%particles(1)%x(2), sim%groups(1)%particles(1)%x(3)

! Convert particle to particle_gc to get energy and mu
gyro_angle = 0.
call sim%fields%calc_EBpsiU(sim%time, sim%groups(1)%particles(1)%i_elm, sim%groups(1)%particles(1)%st, &              
                            sim%groups(1)%particles(1)%x(3), E, B, psi, U)
particle_in = sim%groups(1)%particles(1)
call relativistic_gc_to_particle(sim%fields%node_list, sim%fields%element_list, &
                                 particle_in, particle_out_gc, sim%groups(1)%mass, &
				 B, gyro_angle)
write(*,*) 'Final E (eV), mu (eV/T): ', particle_out_gc%E, particle_out_gc%mu

call write_simulation_hdf5(sim, 'part_restart.h5')

! Finalize the simulation
call sim%finalize

end program ex8_jorek

