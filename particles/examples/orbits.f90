!> Comparing full-orbit and guiding-centre particles in JOREK
program orbits

use particle_tracer
use mpi
use mod_atomic_elements
use mod_particle_io
use mod_event
use mod_particle_loop
use nodes_elements
use mod_random_seed
use mod_math_operators, only: cross_product
use mod_pusher_tools, only: particle_position_to_gc
use mod_gc_variational
use mod_interp, only: mode_moivre, interp_RZ, interp_0
use phys_module, only: F0, tstep, nstep, nout, restart
use phys_module, only: CENTRAL_MASS, CENTRAL_DENSITY, xcase, xpoint
use phys_module, only: n_particles, nstep_particles, nsubstep_particles, tstep_particles
use phys_module, only: xtime, mode
use constants,   only: MU_ZERO, MASS_PROTON, ATOMIC_MASS_UNIT, K_BOLTZ, EL_CHG
use mod_export_restart

use mod_edge_domain
use mod_edge_elements
use data_structure, only: type_bnd_element_list, type_bnd_node_list 
use equil_info
use mod_boundary, only: boundary_from_grid

!$ use omp_lib

implicit none

type(event)                                       :: fieldreader, partreader, partwriter
type(count_action)                                :: counter
type(type_edge_domain), allocatable, dimension(:) :: edge_domains
type(edge_elements)                               :: D_edge
!type(write_particle_diagnostics)                  :: diag
!type(type_bnd_element_list) :: bnd_elm_list !< List of boundary elements
!type(type_bnd_node_list)    :: bnd_node_list !< List of boundary nodes.

real*8    :: target_time, tstep_keep
real*8    :: physical_particles, weight
real*8    :: oldtime, step_rest_time, particle_step_time, particle_start_time, timestep_gc
real*8    :: rho_norm, t_norm, v_norm, E_norm, M_norm, N_norm, tstep_si, timesteps
real*8    :: v_kin_temp, E(3), B(3), psi, U, B_norm, psi_prev
real*8    :: ran(6), T_scale_factor
real*8    :: p_in(3), B_hat(3), v_par, v_par2, qom, larmor_radius, v_perp
real*8    :: A(3), dA(3,3), dB(3,3), bn, dbn(3), Bnorm(3), dBnorm(3,3)
real*8    :: r_out, z_out, s_out, t_out
real*8    :: p_phi_gc, energy_gc, p_phi_gc_start, energy_gc_start 
real*8    :: p_phi_lf, energy_lf, p_phi_lf_start, energy_lf_start
real*8    :: p_phi_qin, energy_qin, p_phi_qin_start, energy_qin_start
real*8,allocatable :: error_W_rk4(:), error_W_qin(:), error_W_lf(:)
real*8,allocatable :: error_P_rk4(:), error_P_qin(:), error_P_lf(:)
integer   :: i_elm_out
!$ real*8 :: w0, w1, mmm(3)

integer   :: ifail, n_part_phi, n_particles_local
integer   :: i, j, n_steps, n_phases
integer   :: seed, i_rng, n_stream, ierr, n_particle_out
character*14 :: fileout, filepart
type(particle_gc_vpar) :: p_gc, p_check_gc
type(particle_gc)      :: p_Emu
type(particle_kinetic_leapfrog), allocatable :: p_orbit(:) 

call sim%initialize(num_groups=3)

tstep_keep  = tstep       ! fluid time step
n_particles_local = 1  
timesteps   = tstep_particles
nstep       = nstep_particles
n_steps     = nsubstep_particles

allocate(error_W_rk4(nstep), error_W_qin(nstep), error_W_lf(nstep))
allocate(error_P_rk4(nstep), error_P_qin(nstep), error_P_lf(nstep))

n_phases = 16  ! number of phase angles for orbit reconstruction from gc

open(111,file='orbits.txt')
open(112,file='orbits_full.txt')

fieldreader = event(read_jorek_fields_interp_linear(basename='jorek', i=-1))
call with(sim, fieldreader)

tstep     = tstep_keep

call det_modes()

if (.not. restart) then
  do j=1, sim%fields%node_list%n_nodes
    sim%fields%node_list%node(j)%values(:,:,2) = 0.d0
  enddo
endif

n_norm   = CENTRAL_DENSITY * 1.d20                              ! (number) density normalisation
rho_norm = CENTRAL_MASS * MASS_PROTON * n_norm                  ! rho_SI = rho_norm * rho
t_norm   = sqrt((MU_ZERO * rho_norm))                           ! t_SI   = t_norm * t_jorek

if (sim%my_id .eq. 0) call boundary_from_grid(sim%fields%node_list, sim%fields%element_list, bnd_node_list, bnd_elm_list, .false.)
call broadcast_boundary(sim%my_id, bnd_elm_list, bnd_node_list)
call update_equil_state(sim%my_id, sim%fields%node_list, sim%fields%element_list, bnd_elm_list, xpoint, xcase)

  
sim%groups(1)%Z    = 1
sim%groups(1)%mass = atomic_weights(1) !< atomic mass units, -2 for deuterium
sim%groups(2)%Z    = 1
sim%groups(2)%mass = atomic_weights(1) !< atomic mass units, -2 for deuterium
sim%groups(3)%Z    = 1
sim%groups(3)%mass = atomic_weights(1) !< atomic mass units, -2 for deuterium
  
allocate(particle_kinetic_leapfrog::sim%groups(1)%particles(n_particles_local))
allocate(particle_gc_vpar::sim%groups(2)%particles(n_particles_local))
allocate(particle_gc_Qin::sim%groups(3)%particles(n_particles_local))

allocate(p_orbit(n_phases))

select type (p_gc => sim%groups(2)%particles)
type is (particle_gc_vpar)
  select type (p_lf => sim%groups(1)%particles)
  type is (particle_kinetic_leapfrog)
  
      do i=1,n_particles_local

        p_lf(1)%q = 1
        p_gc(1)%q = 1

        qom =  p_lf(i)%q * EL_CHG / (sim%groups(1)%mass * ATOMIC_MASS_UNIT)

        p_gc(1)%vpar = 0.4122032306d+05          ! from Qin paper, PoP 16, 042510 (2009) (P_phi (normalised) = -1.077d-3)
        p_gc(1)%mu   = 0.2061748464d+11          ! mu (normalised) = 2.25d-6
        p_gc(1)%x    = (/ 1.05d0, 0.d0, 0.d0 /)

        r_out = p_gc(1)%x(1)
        z_out = p_gc(1)%x(2)
        call find_RZ_nearby(node_list, element_list, p_gc(1)%x(1), p_gc(1)%x(2), p_gc(1)%st(1), p_gc(1)%st(2), -1, &
                            r_out, z_out, s_out, t_out, i_elm_out, ifail)
        write(*,'(A,i3,8e16.8)') ' find_RZ_nearby : ',ifail,p_gc(1)%x(1:2),r_out,z_out
        p_gc(1)%i_elm = i_elm_out
        p_gc(1)%st    = (/ s_out, t_out /) 

        call sim%fields%calc_EBpsiU(sim%time, p_gc(i)%i_elm, p_gc(i)%st, p_gc(i)%x(3), E, B, psi, U)
        call sim%fields%calc_RK4(sim%time, p_gc(1)%i_elm, p_gc(1)%st, p_gc(1)%x(3), A, dA, B, dB, Bnorm, dBnorm, bn, dbn, E)
      !  call sim%fields%calc_RK4_analytic(p_gc(1)%x(1), p_gc(1)%x(2), p_gc(1)%x(3), A, dA, B, dB, Bnorm, dBnorm, bn, dbn, E)

        p_gc(1)%B_norm = bn

        v_perp = sqrt(2.d0*norm2(B)*p_gc(1)%mu)
        larmor_radius = v_perp / norm2(B) /qom
        write(*,'(A,3e18.10)') ' perpendicular velocity at (R,Z)=(1.05,0) : ',v_perp
        write(*,'(A,3e18.10)') ' Larmor radius             (R,Z)=(1.05,0) : ',larmor_radius

        p_lf(1)%x = p_gc(1)%x + (/ 0.d0, -larmor_radius, 0.d0 /)

        r_out = p_lf(1)%x(1)
        z_out = p_lf(1)%x(2)
        call find_RZ_nearby(node_list, element_list, p_lf(1)%x(1), p_lf(1)%x(2), p_lf(1)%st(1), p_lf(1)%st(2), -1, &
                            r_out, z_out, s_out, t_out, i_elm_out, ifail)
        p_lf(1)%i_elm = i_elm_out
        p_lf(1)%st    = (/ s_out, t_out /) 

        call sim%fields%calc_EBpsiU(sim%time, p_lf(i)%i_elm, p_lf(i)%st, p_lf(i)%x(3), E, B, psi, U)

        p_lf(1)%v = p_gc(1)%vpar * B / norm2(B) + v_perp * cross(B,(/0.d0,0.d0,1.d0/)) / norm2(cross(B,(/0.d0,0.d0,1.d0/)))

        write(*,*) 'leapfrog from gc'
        write(*,'(A,3e18.10)') 'Position     [m/s] : ',p_lf(1)%x
        write(*,'(A,3e18.10)') 'Velocity     [m/s] : ',p_lf(1)%v
        write(*,'(A,e18.10)')  'CHECK v_par  [m/s] : ',dot_product(p_lf(1)%v, B) / norm2(B)
        write(*,'(A,3e18.10)') 'CHECK v_perp [m/s] : ',norm2(cross(cross(p_lf(1)%v, B), B)) / norm2(B)**2

        call convert_gc_vpar_to_kinetic(sim%fields%node_list, sim%fields%element_list, p_gc(1), B, sim%groups(2)%mass, 1, p_orbit, ifail)

        write(*,*) ' CHECK convert gc to kinetic'
        write(*,'(A,3e18.10)') 'Position     [m/s] : ',p_orbit(1)%x
        write(*,'(A,3e18.10)') 'Velocity     [m/s] : ',p_orbit(1)%v
        write(*,'(A,e18.10)')  'CHECK v_par  [m/s] : ',dot_product(p_orbit(1)%v, B) / norm2(B)
        write(*,'(A,3e18.10)') 'CHECK v_perp [m/s] : ',norm2(cross(cross(p_orbit(1)%v, B), B)) / norm2(B)**2

        call convert_leapfrog_to_gc_vpar(sim%fields%node_list, sim%fields%element_list, p_lf(i), B, sim%groups(1)%mass, p_check_gc)
        write(*,*) ' CHECK convert leapfrog to gc'
        write(*,'(A,8e18.10)') ' check : x         : ',p_check_gc%x
        write(*,'(A,8e18.10)') ' check : v_par, mu : ',p_check_gc%vpar, p_check_gc%mu

        p_Emu = kinetic_leapfrog_to_gc(sim%fields%node_list, sim%fields%element_list, p_lf(i), E, B, sim%groups(1)%mass, 0.d0)
        write(*,'(A,8e18.10)') 'EMu x  : ',p_Emu%x
        write(*,'(A,8e18.10)') 'EMu E  : ',p_Emu%E
        write(*,'(A,8e18.10)') 'EMu mu / (q/m) : ',p_Emu%mu

!        B_norm = norm2(B)
!        B_hat  = B / B_norm
!        p_in   = sim%groups(1)%mass * p_lf(i)%v
!        p_gc(i)%q = p_lf(i)%q
!        p_gc(i)%x = p_gc(i)%x + ATOMIC_MASS_UNIT*cross_product(p_in,B_hat)/(EL_CHG*real(p_gc(i)%q,8)*B_norm)

!        call particle_position_to_gc(sim%fields%node_list,sim%fields%element_list,            &
!                                     p_lf(i)%x, p_lf(i)%st, p_lf(i)%i_elm, p_in, p_lf(i)%q, B_hat, B_norm, &
!                                     p_gc(i)%x, p_gc(i)%st, p_gc(i)%i_elm)     
                                     
        call sim%fields%calc_RK4(sim%time, p_gc(1)%i_elm, p_gc(1)%st, p_gc(1)%x(3), A, dA, B, dB, Bnorm, dBnorm, bn, dbn, E)
!        call sim%fields%calc_RK4_analytic(p_gc(1)%x(1), p_gc(1)%x(2), p_gc(1)%x(3), A, dA, B, dB, Bnorm, dBnorm, bn, dbn, E)

        p_phi_gc_start  = p_gc(1)%x(1) * ( p_gc(1)%vpar * Bnorm(3) + qom * A(3))
!        energy_gc_start = 0.5d0 * p_gc(1)%vpar**2 + p_gc(1)%mu * bn
        energy_gc_start = 0.5d0 * p_gc(1)%vpar**2 + p_gc(1)%mu * p_gc(1)%B_norm

        write(*,'(A,12e18.10)') ' RK4 start : ',p_gc(1)%x,p_gc(1)%vpar, p_phi_gc_start, energy_gc_start
       
      enddo
    end select
  end select

select type (p_gc => sim%groups(2)%particles)
type is (particle_gc_vpar)
  select type (p_qin => sim%groups(3)%particles)
  type is (particle_gc_qin)
    do i=1, n_particles_local
      
      call copy_particle_gc_Vpar_to_Qin(p_gc(i),p_qin(i))
      
      call initialise_gc_Qin(sim%fields, p_qin(i), sim%groups(3)%mass, timesteps)
      
      call sim%fields%calc_Qin(sim%time, p_qin(1)%i_elm, p_qin(1)%st, p_qin(1)%x(3), A, dA, B, dB, Bnorm, dBnorm, bn, dbn, E)
      !call sim%fields%calc_Qin_analytic(p_qin(1)%x(1), p_qin(1)%x(2), p_qin(1)%x(3), A, dA, B, dB, Bnorm, dBnorm, bn, dbn, E)
      
      p_phi_qin_start  = p_Qin(1)%vpar * Bnorm(3) + qom * A(3)
      energy_qin_start = 0.5d0 * p_Qin(1)%vpar**2 + p_Qin(1)%mu * bn

      write(*,'(A,12e18.10)') ' Qin start : ',p_Qin(1)%x,p_Qin(1)%vpar, p_phi_Qin_start, energy_qin_start

    enddo
  end select
end select

select type (p_lf => sim%groups(1)%particles)
type is (particle_kinetic_leapfrog)  
  call boris_all_initial_half_step_backwards_RZPhi(p_lf, sim%groups(1)%mass, sim%fields, sim%time, timesteps)
  call sim%fields%calc_EBpsiU(sim%time, p_lf(1)%i_elm, p_lf(1)%st, p_lf(1)%x(3), E, B, psi_prev, U)
  call loop_particle_kinetic_leapfrog(sim, timesteps, 1, particle_start_time)
  call sim%fields%calc_EBpsiU(sim%time, p_lf(1)%i_elm, p_lf(1)%st, p_lf(1)%x(3), E, B, psi, U)
  p_phi_lf_start  = p_lf(1)%x(1) * p_lf(1)%v(3) + 0.5d0 * qom * (psi + psi_prev)
  energy_lf_start = 0.5d0 * dot_product(p_lf(1)%v,p_lf(1)%v) 
  write(*,'(A,12e16.8)') ' LF start : ',p_lf(1)%x, p_phi_lf_start, energy_lf_start
end select

timestep_gc = n_steps * timesteps

write(*,*) 'total time : ',timesteps * nstep_particles, timesteps * nstep_particles/ 0.1044656224E-07

do i=1, nstep_particles

  particle_start_time = sim%time
  
  call loop_particle_kinetic_leapfrog(sim, timesteps, n_steps, particle_start_time)

  select type (p_gc => sim%groups(2)%particles)
  type is (particle_gc_vpar)	

    do j=1, n_particles 
      call push_gc_rk4(sim%fields, p_gc(j), sim%groups(2)%mass, timesteps, n_steps, 0)
    enddo

    call sim%fields%calc_RK4(sim%time, p_gc(1)%i_elm, p_gc(1)%st, p_gc(1)%x(3), A, dA, B, dB, Bnorm, dBnorm, bn, dbn, E)
    !call sim%fields%calc_RK4_analytic(p_gc(1)%x(1), p_gc(1)%x(2), p_gc(1)%x(3), A, dA, B, dB, Bnorm, dBnorm, bn, dbn, E)
    
 !   call convert_gc_vpar_to_kinetic(sim%fields%node_list, sim%fields%element_list, p_gc(1), B, sim%groups(2)%mass, 0, p_orbit, ifail)
 !   write(112,'(3e18.10)') (p_orbit(j)%x,j=1,n_phases)

    p_phi_gc  = p_gc(1)%x(1) * ( p_gc(1)%vpar * Bnorm(3) + qom * A(3))
    energy_gc = 0.5d0 * p_gc(1)%vpar**2 + p_gc(1)%mu * p_gc(1)%B_norm

    error_W_rk4(i) = abs((energy_gc - energy_gc_start)/energy_gc_start)
    error_P_rk4(i) = abs((P_phi_gc  - P_phi_gc_start) /P_phi_gc_start)

    write(114,'(A,i6,12e18.10)') 'RK4: ',p_gc(1)%i_elm,p_gc(1)%mu,p_gc(1)%x,p_gc(1)%vpar,(p_phi_gc-p_phi_gc_start)/p_phi_gc_start,(energy_gc-energy_gc_start)/energy_gc_start

  end select

  select type (p_Qin =>sim%groups(3)%particles)
  type is (particle_gc_Qin)	

    do j=1, n_particles 
      call push_gc_Qin(sim%fields, p_Qin(j), sim%groups(3)%mass, timesteps, n_steps)
    enddo

    call sim%fields%calc_Qin(sim%time, p_qin(1)%i_elm, p_qin(1)%st, p_qin(1)%x(3), A, dA, B, dB, Bnorm, dBnorm, bn, dbn, E)
    !call sim%fields%calc_Qin_analytic(p_qin(1)%x(1), p_qin(1)%x(2), p_qin(1)%x(3), A, dA, B, dB, Bnorm, dBnorm, bn, dbn, E)
    
    p_phi_qin  = p_Qin(1)%vpar * Bnorm(3) + qom * A(3)
    energy_qin = 0.5d0 * p_Qin(1)%vpar**2 + p_Qin(1)%mu * bn

    error_W_qin(i) = abs((energy_qin - energy_qin_start)/energy_qin_start)
    error_P_qin(i) = abs((P_phi_qin - P_phi_qin_start)/P_phi_qin_start)

    write(113,'(A,i6,12e18.10)') 'QIN: ',p_Qin(1)%i_elm,p_Qin(1)%mu,p_Qin(1)%x,p_Qin(1)%vpar,(p_phi_Qin-p_phi_qin_start)/p_phi_Qin_start,(energy_Qin-energy_Qin_start)/energy_qin_start

  end select
  
  select type (p_lf => sim%groups(1)%particles)
    type is (particle_kinetic_leapfrog)
  
    psi_prev = psi
    call sim%fields%calc_EBpsiU(sim%time, p_lf(1)%i_elm, p_lf(1)%st, p_lf(1)%x(3), E, B, psi, U)
  
    p_phi_lf  = p_lf(1)%x(1) * p_lf(1)%v(3) + 0.5d0 * qom * (psi + psi_prev)
    energy_lf = 0.5d0 * dot_product(p_lf(1)%v,p_lf(1)%v) 

    error_W_lf(i) = abs((energy_lf - energy_lf_start)/energy_lf_start)
    error_P_lf(i) = abs((P_phi_lf- P_phi_lf_start)/P_phi_lf_start)

  end select
  
  !write(111,'(12e16.8)') sim%groups(1)%particles(1)%x, sim%groups(2)%particles(1)%x, sim%groups(3)%particles(1)%x

enddo

write(*,'(A,12e18.10)') 'LF  : P_phi, energy : ',sim%time, P_phi_lf_start, energy_lf_start, &
                      P_phi_lf, energy_lf, (P_phi_lf-P_phi_lf_start)/P_phi_lf_start,(energy_lf-energy_lf_start)/energy_lf_start
write(*,'(A,12e18.10)') 'RK4 : P_phi, energy : ',sim%time,  P_phi_gc_start, energy_gc_start,&
                      P_phi_gc, energy_gc, (P_phi_gc-P_phi_gc_start)/P_phi_gc_start,(energy_gc-energy_gc_start)/energy_gc_start
write(*,'(A,12e18.10)') 'Qin : P_phi, energy : ',sim%time,  P_phi_qin_start, energy_qin_start, &
                      P_phi_qin, energy_qin, (P_phi_qin-P_phi_qin_start)/P_phi_qin_start,(energy_qin-energy_qin_start)/energy_qin_start
write(*,'(A,8e14.6)') 'LF  max(error) : ',maxval(error_P_lf),maxval(error_W_lf),sum(error_P_lf)/real(nstep,8),sum(error_W_lf)/real(nstep,8)
write(*,'(A,8e14.6)') 'RK4 max(error) : ',maxval(error_P_rk4),maxval(error_W_rk4),sum(error_P_rk4)/real(nstep,8),sum(error_W_rk4)/real(nstep,8)
write(*,'(A,8e14.6)') 'Qin max(error) : ',maxval(error_P_qin),maxval(error_W_qin),sum(error_P_qin)/real(nstep,8),sum(error_W_qin)/real(nstep,8)

call sim%finalize

close(111)
close(112)

contains

    
subroutine loop_particle_kinetic_leapfrog(sim, timesteps, n_steps, particle_start_time)
use mod_random_seed
use mod_interp, only: sincosperiod_moivre
use mod_particle_types, only: copy_particle_kinetic_leapfrog

implicit none

class(particle_sim), target, intent(inout)                :: sim
type(count_action)                                        :: counter
type(particle_kinetic_leapfrog)                           :: particle_tmp

real*8, intent(in)     :: timesteps, particle_start_time 
integer, intent(in)    :: n_steps

real*8    :: n_norm, rho_norm, t_norm, v_norm, E_norm, M_norm
real*8    :: t, E(3), B(3), psi, U, rz_old(2), st_old(2)
real*8    :: v_temp(3)
!$ real*8 :: w0, w1, mmm(3)

integer   :: i, j, k, l, m, i_elm_old, i_elm 
integer   :: seed, i_rng, n_stream, ierr, nthreads
integer   :: i_tor, index_lm, i_elm_temp
integer   :: ifail

!$ w0 = omp_get_wtime()

n_norm   = CENTRAL_DENSITY * 1.d20                              ! (number) density normalisation
rho_norm = CENTRAL_MASS * MASS_PROTON * n_norm                  ! rho_SI = rho_norm * rho
t_norm   = sqrt((MU_ZERO * rho_norm))                           ! t_SI   = t_norm * t_jorek
v_norm   = 1.d0 / t_norm                                        ! V_SI   = v_norm * v_jorek
E_norm   = 1.5d0 / MU_ZERO                                      ! E_SI   = E_norm * E_jorek
M_norm   = rho_norm * v_norm                                    ! momentum normalisation


select type (particles => sim%groups(1)%particles)
type is (particle_kinetic_leapfrog)
#ifdef __GFORTRAN__
   !$omp parallel do default(shared) & 
#else
   !$omp parallel do default(none) &
   !$omp shared(sim, particles, n_steps, timesteps, particle_start_time,          &
   !$omp rho_norm, t_norm, v_norm, E_norm, M_norm, N_norm,                        &
   !$omp central_density, central_mass)                                           &
#endif
   !$omp private(particle_tmp, i_rng, j, k, t, E, B, psi, U, rz_old, st_old,      &
   !$omp i_elm_old, i_elm, ifail)                                                 &
   !$omp schedule(dynamic,10)                                                     
   do j=1,size(particles,1)

      call copy_particle_kinetic_leapfrog(particles(j),particle_tmp)
!      i_rng = 1
    !$ i_rng = omp_get_thread_num()+1

      do k=1,n_steps

        if (particle_tmp%i_elm .le. 0) exit

        t = particle_start_time + (k-1)*timesteps

        call sim%fields%calc_EBpsiU(t, particle_tmp%i_elm, particle_tmp%st, particle_tmp%x(3), E, B, psi, U)

        rz_old    = particle_tmp%x(1:2)
        st_old    = particle_tmp%st
        i_elm_old = particle_tmp%i_elm
                
        call boris_push_cylindrical(particle_tmp, sim%groups(1)%mass, E, B, timesteps)

        call find_RZ_nearby(sim%fields%node_list, sim%fields%element_list, rz_old(1), rz_old(2), st_old(1), st_old(2), i_elm_old, &
                            particle_tmp%x(1), particle_tmp%x(2), particle_tmp%st(1), particle_tmp%st(2), particle_tmp%i_elm, ifail)
        
      end do ! steps

      call copy_particle_kinetic_leapfrog(particle_tmp, particles(j))
    
    end do   ! particles
    !$omp end parallel do
      
  end select

end subroutine


end program orbits
