!> Push relativistic particles in JOREK fields
!>
!> Compile with `make ex6_jorek`
!> Run with `./ex6_jorek < JOREK_namelist`

program ex6_jorek

use particle_tracer
use mod_particle_io
use mod_particle_diagnostics
use mod_fields_linear   
use mod_fields_hermite_birkhoff 
use mod_kinetic_relativistic
use mod_gc_relativistic
                 
implicit none

! Set up the simulation variables
real*8                              :: timesteps(1) = [3.5723d-13]
real*8                              :: target_time, t
real*8                              :: energy     !< Initial energy in eV (including rest energy)
real*8                              :: ksi        !< Initial cosine of pitch-angle
real*8                              :: gyro_angle !< Initial gyro-angle
real*8                              :: E(3), B(3), rz_old(2), st_old(2), psi, U
integer*4                           :: n_part, i, j, k, n_steps, i_elm_old, ifail, n_lost
logical                             :: restart
type(diag_print_kinetic_energy)     :: print_kinetic_energy
type(write_particle_diagnostics)    :: diag
type(particle_kinetic_relativistic) :: particle_kin_rel
type(particle_gc_relativistic)      :: particle_gc_rel, particle_gc_rel2
type(particle_gc)                   :: particle_gc_out
! For interp_PRZ
real*8 :: P(1), P_s(1), P_t(1), P_phi(1), P_time(1)
real*8 :: R, R_s, R_t, Z, Z_s, Z_t
! For CPU time
real*8 :: t0, t1

call cpu_time(t0)

call sim%initialize(num_groups=1)
  
! Set up the diagnostics output
diag = write_particle_diagnostics(filename='diag.h5',only=[1,2,6,12,13,14,15]) ! store total and kinetic energies, p_phi, phi, R, Z

restart = .false. !

if (.not. restart) then

  sim%time = 2.5d-3 ! 1.d-7 !  ! start time 

  n_part = 1
  allocate(particle_kinetic_relativistic::sim%groups(1)%particles(n_part))
  sim%groups(1)%mass = 5.4857990907016d-4 ! particle mass in AMU
  
  ! Set events to write output data and stop the simulation.
  ! One can use read_jorek_fields_interp_linear or read_jorek_fields_interp_hermite_birkhoff,
  ! and i=-1 (to read jorek_restart.h5 and keep this field at all time) or i=last_file_before_time(sim%time)
  ! (to read a sequel of jorekXXXXX.h5 files and use time-evolving fields)
  events = [event(read_jorek_fields_interp_linear(i=last_file_before_time(sim%time))), & 
            event(diag,start=sim%time,step=1.d-7),         &
	    event(stop_action(),start=sim%time+1.d-5)]

  ! Run first event to read the JOREK fields
  call with(sim, events, at=0.d0)

  do i=1,n_part

    select type (p=>sim%groups(1)%particles(i))

    type is (particle_kinetic_relativistic)

      p%q = -1

      p%x = [3.68d0,0.d0,0.d0] ! in (R,Z,phi) coordinates

      call find_RZ(sim%fields%node_list, sim%fields%element_list,    &
                   p%x(1), p%x(2),                                   & ! inputs
                   p%x(1), p%x(2), p%i_elm, p%st(1), p%st(2), ifail)   ! outputs

      !p%p = [0.d0,1.d+7,0.d0] ! in (X,Y,Z) coordinates

      ! Initializing momentum from energy, cosine of pitch-angle and gyro-angle
      ! To do this, we will convert the particle from particle_kinetic_relativistic
      ! into particle_gc_relativistic, then set the p_parallel and mu fields for the latter,
      ! and then convert it back into particle_kinetic_relativistic

      energy     = 1.d7 ! 5.12d5 ! Particle energy, including rest energy
      ksi        = 0. !1.d0      ! Cosine of pitch-angle
      gyro_angle = 0.

      call sim%fields%calc_EBpsiU(sim%time, p%i_elm, p%st, p%x(3), E, B, psi, U)
      
      ! Convert to particle_gc_relativistic
      particle_gc_rel = relativistic_kinetic_to_relativistic_gc(sim%fields%node_list, &
        sim%fields%element_list, p, sim%groups(1)%mass, B)

      ! Set p_parallel and mu
      particle_gc_rel2 = relativistic_gc_momenta_from_E_cospitch(particle_gc_rel, energy, &
        ksi, sim%groups(1)%mass, sim%fields, sim%time)

      ! Convert back to particle_kinetic_relativistic
      particle_kin_rel = relativistic_gc_to_relativistic_kinetic(sim%fields%node_list, &
        sim%fields%element_list, particle_gc_rel2, sim%groups(1)%mass, B, gyro_angle)

      p%x = particle_kin_rel%x
      p%p = particle_kin_rel%p

      ! Checks
      ! 1) Convert particle to particle_gc to get energy and mu
      call relativistic_kinetic_to_particle(sim%fields%node_list, sim%fields%element_list, &
                                      p, particle_gc_out, sim%groups(1)%mass, B)
      write(*,*) 'Initial E (eV), mu (eV/T): ', particle_gc_out%E, particle_gc_out%mu

      ! 2) Convert to relativistic GC as another way to get mu
      particle_gc_rel = relativistic_kinetic_to_relativistic_gc(sim%fields%node_list, &
        sim%fields%element_list, p, sim%groups(1)%mass, B)
      write(*,*) 'Initial mu from relativistic GC: ', particle_gc_rel%p(2)

    end select
  end do

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
	    event(stop_action(),start=sim%time+5.d-8)]

  ! Run first event to read the JOREK fields
  call with(sim, events, at=0.d0)   

endif

call check_and_fix_timesteps(timesteps, events)

! Set dpsi/dt=0 (useful e.g. to check the conservation of the total particle energy) 
!call sim%fields%set_flag_dpsidt(.true.)

! Open a file where to write some fields at a given position to test time interpolation routines
!open(22,file='field_vs_t.dat')

! Loop until the simulation is stopped
do while (.not. sim%stop_now)
  ! Extract the next event time
  target_time = next_event_at(sim, events)
  ! Loop over all particle groups
  do i=1,1
    ! Compute the number of time steps
    n_steps = nint((target_time - sim%time)/timesteps(i))
    write(*,*) 'n_steps', n_steps
    n_lost = 0
    select type (particles => sim%groups(i)%particles)
    type is (particle_kinetic_relativistic)	
      ! Loop on particles whithin the i-th group
!      !$omp parallel do default(private) &
!      !$omp shared (i, n_steps, timesteps, sim) &
!      !$omp reduction(+:n_lost)	
      do j=1,size(particles,1)
        do k=1,n_steps
          if (particles(j)%i_elm .le. 0) exit
	  sim%time = sim%time + timesteps(i)

          !call volume_preserving_push_analytical(particles(j),sim%fields,sim%groups(i)%mass,sim%time,timesteps(i))
          call volume_preserving_push_jorek(particles(j),sim%fields,sim%groups(i)%mass,sim%time,timesteps(i),ifail)

!	  if (modulo(k-1,10000)==0) then	                
!	    call sim%fields%interp_PRZ(sim%time, 1000, [1], 1, 0.5, 0.5, 0.5, P, P_s, P_t, P_phi, P_time, R, R_s, R_t, Z, Z_s, Z_t)
!	    write(22,'(7e26.16)') sim%time, P, P_time, R, Z
!	  end if

	  if (particles(j)%i_elm .le. 0) n_lost = n_lost + 1	
        end do !< time steps
      end do !< particles
!      !$omp end parallel do
    end select
    write(*,*) "number of lost particles: ", n_lost	  
  end do !< groups

  ! Update current time and run events
  sim%time = target_time
  call with(sim, events, at=sim%time)
enddo

! Print final particle information
!write(*,*) 'Final x,y,z: ', sim%groups(1)%particles(1)%x(1), sim%groups(1)%particles(1)%x(2), sim%groups(1)%particles(1)%x(3)
!call print_kinetic_energy%do(sim,events(1))  !< print particle kinetic energy

! Convert particle to particle_gc to get energy and mu
particle_kin_rel = sim%groups(1)%particles(1)
call sim%fields%calc_EBpsiU(sim%time, sim%groups(1)%particles(1)%i_elm, sim%groups(1)%particles(1)%st, &              
                            sim%groups(1)%particles(1)%x(3), E, B, psi, U)
call relativistic_kinetic_to_particle(sim%fields%node_list, sim%fields%element_list, &
                                      particle_kin_rel, particle_gc_out, sim%groups(1)%mass, B)
write(*,*) 'Final E (eV), mu (eV/T): ', particle_gc_out%E, particle_gc_out%mu

! Convert to relativistic GC as another way to get mu
particle_gc_rel = relativistic_kinetic_to_relativistic_gc(sim%fields%node_list, &
        sim%fields%element_list, particle_kin_rel, sim%groups(1)%mass, B)
write(*,*) 'Final mu from relativistic GC: ', particle_gc_rel%p(2)

call cpu_time(t1)
write(*,*) 'CPU time: ', t1-t0
 
call write_simulation_hdf5(sim, 'part_restart.h5')

! Finalize the simulation
call sim%finalize

end program ex6_jorek

