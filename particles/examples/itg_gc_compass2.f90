!> Testing the coupling of the projections of particles to JOREK
module domain
  real*8 :: psi_start, psi_end, psi_n_start, psi_n_end
end module

program itg_loop

use particle_tracer
use mod_pusher_tools, only : particle_position_to_gc
use mod_particle_diagnostics
use mpi
use mod_atomic_elements
use mod_particle_io
use mod_event
use mod_project_particles
use mod_particle_loop
use mod_gc_variational
use nodes_elements
use mod_jorek_timestepping
use mod_random_seed
use mod_interp, only: mode_moivre, interp_RZ, interp_0
use mod_basisfunctions
use basis_at_gaussian
use phys_module, only: F0, tstep, nstep, nout, restart, rho_0, rho_1, rho_coef
use phys_module, only: CENTRAL_MASS, CENTRAL_DENSITY, xcase, xpoint, index_now, index_start
use phys_module, only: n_particles, nstep_particles, nsubstep_particles, tstep_particles
use phys_module, only: filter_perp, filter_hyper, filter_par, filter_perp_n0, filter_hyper_n0, filter_par_n0
use phys_module, only: xtime, energies, mode, restart_particles
use constants,   only: MU_ZERO, MASS_PROTON, ATOMIC_MASS_UNIT, K_BOLTZ, EL_CHG
use mod_export_restart
use live_data

use mod_edge_domain
use mod_edge_elements
use data_structure, only: type_bnd_element_list, type_bnd_node_list 
use equil_info
use mod_boundary, only: boundary_from_grid
use domain

!$ use omp_lib

implicit none

type(event)                                       :: fieldreader, partreader, partwriter
type(count_action)                                :: counter
type(projection), target                          :: jorek_feedback, project_profiles
type(jorek_timestep_action), target               :: jorek_stepper
type(type_edge_domain), allocatable, dimension(:) :: edge_domains
type(edge_elements)                               :: D_edge
type(write_particle_diagnostics)                  :: diag
!type(type_bnd_element_list) :: bnd_elm_list !< List of boundary elements
!type(type_bnd_node_list)    :: bnd_node_list !< List of boundary nodes.

real*8    :: target_time, tstep_keep
real*8    :: physical_particles, weight
real*8    :: oldtime, step_rest_time, particle_step_time, particle_start_time, diag_time
real*8    :: rho_norm, t_norm, v_norm, E_norm, M_norm, zn_norm, tstep_si, timesteps
real*8    :: v_kin_temp, E(3), B(3), psi, U, B_norm
real*8    :: rescale_coef, T_axis(1), E_axis, E_hot, rho_part, v2, delta_phi
real*8    :: W_mag(n_tor), W_kin(n_tor), ran(6), T_scale_factor, Te0_eV, alfa, growth_kin
real*8    :: psi_n, dpsi_n, ss, psi_axis, psi_bnd, average_potential, previous_potential
real*8    :: w_alive, w_alive_total, w_alive_previous, W_thermal_local, W_thermal_total
real*8    :: sum_rhs(4), sum_rhs_local(4), sum_weights(4)
real*8    :: fraction_particles, total_particles, total_volume, p_in(3), B_hat(3), v_par, v_par2
real*8    :: zne0, zne0_s, zne0_st, zne0_t, Te0, Tev_norm, Te0_s, Te0_st, Te0_t

!$ real*8 :: w0, w1, mmm(3)

integer   :: ifail, n_part_phi, n_particles_local, node_start, node_end
integer   :: i, j, k, l, m, in, jn, inode, i_elm, n_steps, i_elm_old, index_rhs, i_diagno(3)
integer   :: seed, i_rng, n_stream, ierr, n_particle_out, i_tor
logical   :: compensate_n0_field
character*14 :: fileout, filepart

real*8, allocatable :: rhs_nodes(:,:), rhs_nodes_local(:,:)


psi_n_start = 0.d0          ! limit domain to psi_start:psi_end (in normalised psi)
psi_n_end   = 1.d0 

T_scale_factor = 1.d0 !1.d4                ! for GENE benchmark

n_particle_out = -1

compensate_n0_field = .false.

call sim%initialize(num_groups=1)

if (sim%my_id .eq. 0) write(*,*) 'RESTART = ',restart

tstep_keep  = tstep       ! fluid time step
n_particles_local = int(n_particles/sim%n_cpu) 
timesteps   = tstep_particles
nstep       = nstep_particles
n_steps     = nsubstep_particles

write(*,*) sim%my_id,' number of particles : ',n_particles_local,int(n_particles)
open(111,file='diagno.txt')

if (sim%my_id .eq. 0) call init_live_data()

if (restart_particles) then
  deallocate(sim%groups)
  allocate(sim%groups(0))

  partreader = event(read_action(filename='part_restart.h5'))
  call with(sim, partreader)
endif

fieldreader = event(read_jorek_fields_interp_linear(basename='jorek', i=-1))
call with(sim, fieldreader)

tstep     = tstep_keep
index_now = index_start

i_diagno(1) =   sim%fields%node_list%n_nodes / 4 + 16 
i_diagno(2) = 2*sim%fields%node_list%n_nodes / 4 + 16 
i_diagno(3) = 3*sim%fields%node_list%n_nodes / 4 + 16 
if (sim%my_id .eq. 0) write(*,'(A,6f8.4)') ' probe at : ',sim%fields%node_list%node(i_diagno(1))%x(1,1,1:2), &
                                                          sim%fields%node_list%node(i_diagno(2))%x(1,1,1:2), &
                                                          sim%fields%node_list%node(i_diagno(3))%x(1,1,1:2)

call det_modes()

if (.not. restart) then
  do j=1, sim%fields%node_list%n_nodes
    sim%fields%node_list%node(j)%values(:,:,2) = 0.d0
    sim%fields%node_list%node(j)%values(:,:,6) = T_scale_factor * sim%fields%node_list%node(j)%values(:,:,6)
  enddo
endif

if (sim%my_id .eq. 0) then
  do i=1, index_start
    call write_live_data(i)
  enddo 
endif

zn_norm   = CENTRAL_DENSITY * 1.d20                              ! (number) density normalisation
rho_norm  = CENTRAL_MASS * MASS_PROTON * zn_norm                 ! rho_SI = rho_norm * rho
t_norm    = sqrt((MU_ZERO * rho_norm))                           ! t_SI   = t_norm * t_jorek
Tev_norm  = 1.d0 / (2.d0 * EL_CHG * MU_ZERO * zn_norm)           ! T_ev   = Tev_norm * T (factor 2 for electron temperature)

if (sim%my_id .eq. 0) write(*,*) ' n_norm : ',zn_norm

if (sim%my_id .eq. 0) call boundary_from_grid(sim%fields%node_list, sim%fields%element_list, bnd_node_list, bnd_elm_list, .false.)

call broadcast_boundary(sim%my_id, bnd_elm_list, bnd_node_list)
call update_equil_state(sim%my_id, sim%fields%node_list, sim%fields%element_list, bnd_elm_list, xpoint, xcase)

psi_axis  = ES%psi_axis
psi_bnd   = ES%psi_bnd
psi_start = psi_axis + (psi_bnd - psi_axis) * psi_n_start
psi_end   = psi_axis + (psi_bnd - psi_axis) * psi_n_end

if (sim%my_id .eq. 0) write(*,'(A,6e14.6)') 'PSI_AXIS, PSI_BND : ',psi_axis, psi_bnd, psi_start, psi_end, psi_n_start, psi_n_end


if (.not. restart_particles) then
  
  sim%groups(1)%Z    = 1
  sim%groups(1)%mass = atomic_weights(-2) !< atomic mass units, -2 for deuterium
  
  if (sim%my_id .eq. 0) write(*,'(A,e12.4)') ' ion mass : ',sim%groups(1)%mass

  allocate(particle_gc_vpar::sim%groups(1)%particles(n_particles_local))

  call initialise_particles_H_mu_psi(sim%groups(1)%particles, sim%fields, sobseq_rng(),sim%groups(1)%mass, &
                                     uniform_space=.true., uniform_space_rej_f=f_density, &
                                     uniform_space_rej_vars=[1], charge = 1)

  call density_integral(sim%fields%node_list,sim%fields%element_list,f_density,min(psi_axis,psi_bnd),max(psi_axis,psi_bnd),total_particles,total_volume) 
  call density_integral(sim%fields%node_list,sim%fields%element_list,f_density,min(psi_axis,psi_bnd),max(psi_axis,psi_bnd),fraction_particles,total_volume)

  total_particles    = total_particles    * zn_norm
  fraction_particles = fraction_particles * zn_norm
  rho_part           = fraction_particles
  if (sim%my_id .eq. 0) write(*,'(A,8e14.6)') ' total/fraction particles, volume : ',total_particles, fraction_particles, total_volume

  call with(sim, counter)

endif ! not restart

! SHOULD MAKE PROFILES USING THE MULTIPLE KINETIC POINTS PER GYROCENTRE

!project_profiles = new_projection(sim%fields%node_list, sim%fields%element_list, &
!                      filter    = filter_perp,    filter_hyper    = filter_hyper,    filter_parallel    = filter_par, &
!                      filter_n0 = filter_perp_n0, filter_hyper_n0 = filter_hyper_n0, filter_parallel_n0 = filter_par_n0, &
!                      f=[proj_f(proj_one, group = 1),proj_f(proj_Pressure, group = 1)], fractional_digits = 9,  &
!                      do_zonal = .false., calc_integrals=.true., to_vtk=.true., to_h5=.false., basename='profiles', nsub=5)
!call with(sim, project_profiles)

call MPI_BARRIER(MPI_COMM_WORLD, ierr)

if (nstep .gt. 0) then
  jorek_feedback = new_projection(sim%fields%node_list, sim%fields%element_list, &
                     filter    = filter_perp,    filter_hyper    = filter_hyper,    filter_parallel    = filter_par,    &
                     filter_n0 = filter_perp_n0, filter_hyper_n0 = filter_hyper_n0, filter_parallel_n0 = filter_par_n0, &
                     fractional_digits = 9, &
                     do_zonal = .true., calc_integrals=.false., to_vtk=.false., to_h5 = .false., basename='projections')

  allocate(jorek_feedback%rhs(n_order+1, n_vertex_max, sim%fields%element_list%n_elements, n_tor, 1))
  jorek_feedback%rhs = 0.d0
  
  aux_node_list => jorek_feedback%node_list

  events = [new_event_ptr(jorek_feedback,  start = sim%time), event(stop_action(), start=1d12)  ]

! Call events at sim%time once to help event scheduler, before entering particle loop
  step_rest_time = 0.d0
  call with(sim, events, at=sim%time)

  jorek_feedback%scaling_integral_weights = 1.00d0 ! subtract 1.d0 from the rhs of n=0

endif

allocate(rhs_nodes(4,sim%fields%node_list%n_nodes),rhs_nodes_local(4,sim%fields%node_list%n_nodes))

node_start = 1
node_end   = sim%fields%node_list%n_nodes !- 3*128

if (sim%my_id .eq. 0) then
  
  do inode = 1, sim%fields%node_list%n_nodes
    if ((inode .lt. node_start) .or. (inode .gt. node_end)) then
      do jn=1,n_order+1
        index_rhs = 2*(sim%fields%node_list%node(inode)%index(jn)-1) + 1
        jorek_feedback%integral_weights(index_rhs) = 0.d0
      enddo
    endif
  enddo

  sum_weights  = 0.d0
  do inode = 1, sim%fields%node_list%n_nodes
    do jn=1, 4
      index_rhs      = 2*(sim%fields%node_list%node(inode)%index(jn)-1) + 1
      sum_weights(jn) = sum_weights(jn) + jorek_feedback%integral_weights(index_rhs)
    enddo
  enddo
  
endif

!===============================================================================================

do i=1, nstep_particles

  particle_start_time = sim%time

  index_now = index_now + 1

  call loop_particle_gc_local(sim, jorek_feedback, timesteps, n_steps, particle_start_time, .true.)
 
  sim%time = particle_start_time + timesteps * n_steps

  if (sim%my_id == 0) xtime(index_now) = sim%time

  do i_elm=1,sim%fields%element_list%n_elements    
    do in=1,n_vertex_max      
      inode = sim%fields%element_list%element(i_elm)%vertex(in)
      if ( (inode .lt. node_start) .or. (inode .gt. node_end) ) then
        jorek_feedback%rhs(:, in, i_elm, 1, 1)    = 0.d0
      endif
    enddo
  enddo

 
  call with(sim, jorek_feedback)

  do j=1, jorek_feedback%node_list%n_nodes

    if (sim%fields%node_list%node(j)%boundary .eq. 0) then

      do i_tor = 2, n_tor

        zne0    = max(0.01, sim%fields%node_list%node(j)%values(1,1,5)) * zn_norm  ! in m^3
        zne0_s  = sim%fields%node_list%node(j)%values(1,2,5) * zn_norm  
        zne0_t  = sim%fields%node_list%node(j)%values(1,3,5) * zn_norm
        zne0_st = sim%fields%node_list%node(j)%values(1,4,5) * zn_norm
 
        Te0_eV  = sim%fields%node_list%node(j)%values(1,1,6) * Tev_norm   ! in eV
        Te0_s   = sim%fields%node_list%node(j)%values(1,2,6) * Tev_norm
        Te0_t   = sim%fields%node_list%node(j)%values(1,3,6) * Tev_norm
        Te0_st  = sim%fields%node_list%node(j)%values(1,4,6) * Tev_norm

        sim%fields%node_list%node(j)%values(i_tor,1,2) = jorek_feedback%node_list%node(j)%values(i_tor,1,1) * Te0_eV

        sim%fields%node_list%node(j)%values(i_tor,2,2) = jorek_feedback%node_list%node(j)%values(i_tor,2,1) * Te0_eV &
                                                       + jorek_feedback%node_list%node(j)%values(i_tor,1,1) * Te0_s  &
                                                       - jorek_feedback%node_list%node(j)%values(i_tor,1,1) * Te0_eV * zne0_s / zne0 

        sim%fields%node_list%node(j)%values(i_tor,3,2) = jorek_feedback%node_list%node(j)%values(i_tor,3,1) * Te0_eV  &
                                                       + jorek_feedback%node_list%node(j)%values(i_tor,1,1) * Te0_t   &
                                                       - jorek_feedback%node_list%node(j)%values(i_tor,1,1) * Te0_eV * zne0_t / zne0     

        sim%fields%node_list%node(j)%values(i_tor,4,2) = jorek_feedback%node_list%node(j)%values(i_tor,4,1) * Te0_eV &
                                                       + jorek_feedback%node_list%node(j)%values(i_tor,1,1) * Te0_st &
                                                       + jorek_feedback%node_list%node(j)%values(i_tor,2,1) * Te0_t  &
                                                       + jorek_feedback%node_list%node(j)%values(i_tor,3,1) * Te0_s  &

                                                       - jorek_feedback%node_list%node(j)%values(i_tor,1,1) * Te0_t  * zne0_s / zne0 &
                                                       - jorek_feedback%node_list%node(j)%values(i_tor,1,1) * Te0_s  * zne0_t / zne0 &
                                                       - jorek_feedback%node_list%node(j)%values(i_tor,3,1) * Te0_eV * zne0_s / zne0 &
                                                       - jorek_feedback%node_list%node(j)%values(i_tor,2,1) * Te0_eV * zne0_t / zne0 &

                                                       + jorek_feedback%node_list%node(j)%values(i_tor,1,1) / (zne0**2) &
                                                           * (2.d0 * zne0_s  * zne0_t - zne0 * zne0_st) * Te0_eV  

        sim%fields%node_list%node(j)%values(i_tor,:,2) = sim%fields%node_list%node(j)%values(i_tor,:,2) / (F0 * zne0) * t_norm
      
      enddo

    else

      sim%fields%node_list%node(j)%values(:,:,2) = 0.d0
    
    endif

  enddo

  if (sim%my_id == 0) write(111,'(12e20.12)') sim%time, sim%fields%node_list%node(i_diagno(1))%values(1:n_tor,1,2), &
                                                        sim%fields%node_list%node(i_diagno(2))%values(1:n_tor,1,2), &
                                                        sim%fields%node_list%node(i_diagno(3))%values(1:n_tor,1,2)  


  if (sim%my_id .eq. 0) then
    call itg_energy(node_list,element_list,min(psi_start,psi_end),max(psi_start,psi_end),W_kin)
    energies(:,1,index_now) = W_kin(:)
    call write_live_data(index_now)
    if (index_now > index_start+1) then
      growth_kin = 0.d0
      if (energies(n_tor,1,index_now-1) .gt. 0.d0) then
        growth_kin = 0.5d0*log(abs(energies(n_tor,1,index_now)/energies(n_tor,1,index_now-1)))/ timesteps
      endif
    endif
    write(*,'(A,32e14.6)') 'energies   : ',sim%time, W_kin(2:n_tor)
    write(*,'(A,32e14.6)') 'growth rate: ',sim%time, growth_kin
  endif

  if ( (sim%my_id .eq. 0) .and. (mod(index_now,nout).eq. 0) ) then

    write(fileout,'(A5,i5.5)') 'jorek',index_now
  
    call export_restart(sim%fields%node_list, sim%fields%element_list, fileout)

!    call with(sim, project_profiles)

  endif

  if (n_particle_out .gt. 0) then
    if (mod(index_now,n_particle_out) == 0) then
      write(filepart,'(A4,i5.5,A3)') 'part',index_now,'.h5'
      partwriter = event(write_action(filename=filepart))
      call with(sim, partwriter)
    endif
  endif

end do

if (nstep_particles .gt. 0) then
  
  close(111)

  if (sim%my_id .eq. 0) call export_restart(sim%fields%node_list, sim%fields%element_list, 'restart_jorek')

  partwriter = event(write_action(filename='part_restart.h5'))
  call with(sim, partwriter)

  if ( sim%my_id .eq. 0 ) call finalize_live_data()
  
  call sim%finalize

endif

contains

function f_density(n, P, grad_P) result(f)
  use domain
  integer, intent(in) :: n
  real*8,  intent(in) :: P(n), grad_P(3,n)
  real*8 :: psi, psi_axis, psi_bnd, Z, Z_xpoint, psi_n
  real*8 :: zn,dn_dpsi,dn_dz,dn_dpsi2,dn_dz2,dn_dpsi_dz,dn_dpsi3,dn_dpsi_dz2, dn_dpsi2_dz
  real*4 :: f
  logical :: xpoint
  integer :: xcase

  
  psi      = P(1)
  psi_axis = ES%psi_axis
  psi_bnd  = ES%psi_bnd
  Z        = 0.d0
  Z_xpoint = -99.d0
  xpoint   = .false.
  xcase    = 1

  call density(xpoint, xcase, Z, Z_xpoint, psi, psi_axis, psi_bnd,&
               zn,dn_dpsi,dn_dz,dn_dpsi2,dn_dz2,dn_dpsi_dz,dn_dpsi3,dn_dpsi_dz2, dn_dpsi2_dz)

  psi_n = (psi - psi_axis)/(psi_bnd - psi_axis)

  if ((psi_n .ge. psi_n_start) .and. (psi_n .le. psi_n_end)) then
    f = zn
  else
    f = 0.d0
  endif

end function f_density


pure function f_itg(n, P, grad_P) result(f)
  integer, intent(in) :: n
  real*8, intent(in) :: P(n), grad_P(3,n)
  real*8 :: psi_n, prof
  real*4 :: f
  
  psi_n = (P(1) - ES%Psi_axis) / ( ES%Psi_bnd - ES%Psi_axis)

  prof = (rho_0-rho_1)*(1.d0 + rho_coef(1) * psi_n + rho_coef(2) * psi_n**2 + rho_coef(3) * psi_n**3)

  prof = prof * (0.5d0 - 0.5d0*tanh((psi_n - rho_coef(5))/rho_coef(4)))

  f = prof + rho_1

end function f_itg


subroutine loop_particle_gc_local(sim, jorek_feedback, timesteps, n_steps, particle_start_time, update)
use mod_project_particles
use mod_random_seed
use mod_interp, only: sincosperiod_moivre
use mod_basisfunctions
use mod_particle_types, only: copy_particle_kinetic_leapfrog

implicit none

class(particle_sim), target, intent(inout)    :: sim
type(projection), target, intent(inout)       :: jorek_feedback
type(count_action)                            :: counter
type(particle_gc_vpar)                        :: particle_tmp
type(particle_kinetic_leapfrog), allocatable  :: p_orbit(:) 

real*8, intent(in)     :: timesteps, particle_start_time 
logical, intent(in)    :: update

real*8    :: n_norm, rho_norm, t_norm, v_norm, E_norm, M_norm
real*8    :: t, E(3), B(3), psi, U, zne0, Te0, rz_old(2), st_old(2)
real*8    :: v_temp(3), T_eV, K_eV
real*8    :: v, v_s, v_t, v_R, v_Z
real*8    :: R_g, Z_g, R_s, R_t, Z_s, Z_t, xjac
real*8    :: HHZ(n_tor), HHZ_p(n_tor), HH(4,4), HH_s(4,4), HH_t(4,4)
!$ real*8 :: w0, w1, mmm(3)

integer, intent(in)   :: n_steps
integer   :: i, j, k, l, m, i_elm_old, i_elm 
integer   :: seed, i_rng, n_stream, ierr, nthreads
integer   :: i_tor, index_lm, i_elm_temp, n_phases
integer   :: ifail

!$ w0 = omp_get_wtime()

n_norm   = CENTRAL_DENSITY * 1.d20                              ! (number) density normalisation
rho_norm = CENTRAL_MASS * MASS_PROTON * n_norm                  ! rho_SI = rho_norm * rho
t_norm   = sqrt((MU_ZERO * rho_norm))                           ! t_SI   = t_norm * t_jorek
v_norm   = 1.d0 / t_norm                                        ! V_SI   = v_norm * v_jorek
E_norm   = 1.5d0 / MU_ZERO                                      ! E_SI   = E_norm * E_jorek
M_norm   = rho_norm * v_norm                                    ! momentum normalisation

jorek_feedback%rhs_gather_time = jorek_feedback%rhs_gather_time + n_steps * timesteps

jorek_feedback%rhs = 0.d0

n_phases = 8
allocate(p_orbit(n_phases))

!call with(sim, counter)

select type (particles => sim%groups(1)%particles)
type is (particle_gc_vpar)

if (sim%my_id .eq. 0) write(*,*) 'starting loop gc : ',size(particles,1),n_phases

#ifdef __GFORTRAN__
   !$omp parallel do default(shared) & ! workaround for Error: �__vtab_mod_pcg32_rng_Pcg32_rng� not specified in enclosing �parallel�
#else
   !$omp parallel do default(none) &
#endif
   !$omp schedule(dynamic,10)                                                     &
   !$omp shared(sim, particles, n_steps, timesteps, particle_start_time, update,  &
   !$omp rho_norm, t_norm, v_norm, E_norm, M_norm, N_norm, n_phases,              &
   !$omp jorek_feedback, central_density, central_mass)                           &
   !$omp private(particle_tmp, i_rng, i,j,k,l,m, t, E, B, psi, U, rz_old, st_old, &
   !$omp i_elm_old, i_elm, zne0, Te0, p_orbit,                                    & 
   !$omp R_g, R_s, R_t, Z_g, Z_s, Z_t, xjac, HH, HH_s, HH_t, index_lm,            &
   !$omp ifail, v, v_s, v_t, v_R, v_Z, HHZ, HHZ_p)
   do j=1,size(particles,1)

    call copy_particle_gc_vpar(particles(j),particle_tmp)
    !  i_rng = 1
    !$ i_rng = omp_get_thread_num()+1

      call push_gc_rk4(sim%fields, particle_tmp, sim%groups(1)%mass, timesteps, n_steps, n_phases) ! add B to output of push_gc_rk4

      call copy_particle_gc_vpar(particle_tmp, particles(j))

      if (particle_tmp%i_elm .gt. 0) then

        call sim%fields%calc_EBpsiU(sim%time, particle_tmp%i_elm, particle_tmp%st, particle_tmp%x(3), E, B, psi, U)

        call convert_gc_vpar_to_kinetic(sim%fields%node_list, sim%fields%element_list, particle_tmp, B, sim%groups(1)%mass, n_phases, p_orbit, ifail)
        
        if (ifail .ne. 0) cycle

        do i=1, n_phases

          if (p_orbit(i)%i_elm .gt. 0) then

            i_elm = p_orbit(i)%i_elm

            call basisfunctions(p_orbit(i)%st(1), p_orbit(i)%st(2), HH, HH_s, HH_t)
  
            call mode_moivre(p_orbit(i)%x(3), HHZ)

!           zne0 = f_density(1, [psi], [0.d0,0.d0,0.d0]) * central_density * 1d20
!           zne0 = central_density * 1d20
!           call sim%fields%calc_NeTe(t, p_orbit(i)%i_elm, p_orbit%st, p_orbit%x(3), zne0, Te0)
!           zne0   = 4.66e19  ! for GENE benchmark
  
            do l=1,n_vertex_max
              do m=1,n_order+1
  
                index_lm = (l-1)*(n_order+1) + m
  
                v   = HH(l,m) * sim%fields%element_list%element(i_elm)%size(l,m) * particle_tmp%weight

                do i_tor=1,n_tor
                  jorek_feedback%rhs(m,l,i_elm,i_tor,1) = jorek_feedback%rhs(m,l,i_elm,i_tor,1) + HHZ(i_tor) * v  
                enddo
    
              enddo   !< order
            enddo     !< vertex
          endif
        enddo       !< phases of gc orbit

      endif
  
    end do   ! particles
    !$omp end parallel do
    
  end select

  t = particle_start_time + n_steps*timesteps

  jorek_feedback%rhs = jorek_feedback%rhs /real(n_phases,8)
   
  !  write(*,*) 'CAREFUL: averaging over n_steps : ',n_steps
  !  jorek_feedback%rhs = jorek_feedback%rhs / real(n_steps,8)
 
  !$ w1 = omp_get_wtime()
  !$ mmm = mpi_minmeanmax(w1-w0)
  !$ if (sim%my_id .eq. 0) write(*,"(A,3f9.4,A)") " Particle stepping complete !in ", mmm, "s"

  if (sim%my_id .eq. 0) write(*,*) 'done loop_particle_gc_local'
  
end subroutine


subroutine density_integral(node_list,element_list,f_density,psi_start,psi_end,total_particles,total_volume)
  !---------------------------------------------------------------
  !
  !---------------------------------------------------------------
  use data_structure
  use gauss
  use basis_at_gaussian
  use phys_module
  
  implicit none
  
  type (type_node_list)    :: node_list
  type (type_element_list) :: element_list
  type (type_element)      :: element
  type (type_node)         :: nodes(n_vertex_max)
  real*8, intent(in)       :: psi_start, psi_end
  real*8, intent(out)      :: total_particles, total_volume
  real*4, external         :: f_density
  
  real*8     :: x_g(n_gauss,n_gauss),  x_s(n_gauss,n_gauss),  x_t(n_gauss,n_gauss)
  real*8     :: y_s(n_gauss,n_gauss),  y_t(n_gauss,n_gauss)
  real*8     :: eq_g(n_gauss,n_gauss), psi_eq(n_gauss,n_gauss)
  integer    :: i, j, ms, mt, iv, inode, ife, n_elements
  real*8     :: xjac, BigR, wst
  real*8     :: sum_particles, sum_volume
  integer, parameter :: ivar_psi=1, ivar_rho=5
  
 sum_particles = 0.d0
 sum_volume    = 0.d0
  
!$omp parallel do default(none) &
!$omp schedule(dynamic,10)                                                         &
!$omp shared(element_list, node_list, psi_start, psi_end, H, H_s, H_t)             &
!$omp private(ife, element, nodes, iv, inode, x_g, x_s, x_t, y_s, y_t, eq_g,       &
!$omp         i, j, ms, mt, wst, bigR, xjac, psi_eq )              &
!$omp reduction(+:sum_particles, sum_volume)
  do ife =1,  element_list%n_elements
  
    element = element_list%element(ife)
  
    do iv = 1, n_vertex_max
      inode     = element%vertex(iv)
      nodes(iv) = node_list%node(inode)
    enddo
  
    x_g(:,:)  = 0.d0;   x_s(:,:)    = 0.d0;   x_t(:,:) = 0.d0;
    y_s(:,:)  = 0.d0;   y_t(:,:)    = 0.d0;
    eq_g(:,:) = 0.d0;   psi_eq(:,:) = 0.d0
  
    do i=1,n_vertex_max
      do j=1,n_order+1
        do ms=1, n_gauss
          do mt=1, n_gauss
  
            x_g(ms,mt) = x_g(ms,mt) + nodes(i)%x(1,j,1) * element%size(i,j) * H(i,j,ms,mt)
            x_s(ms,mt) = x_s(ms,mt) + nodes(i)%x(1,j,1) * element%size(i,j) * H_s(i,j,ms,mt)
            x_t(ms,mt) = x_t(ms,mt) + nodes(i)%x(1,j,1) * element%size(i,j) * H_t(i,j,ms,mt)
            y_s(ms,mt) = y_s(ms,mt) + nodes(i)%x(1,j,2) * element%size(i,j) * H_s(i,j,ms,mt)
            y_t(ms,mt) = y_t(ms,mt) + nodes(i)%x(1,j,2) * element%size(i,j) * H_t(i,j,ms,mt)
  
            psi_eq(ms,mt) = psi_eq(ms,mt) + nodes(i)%values(1,j,ivar_psi) * element%size(i,j) * H(i,j,ms,mt)
            eq_g(ms,mt)   = eq_g(ms,mt)   + nodes(i)%values(1,j,ivar_rho) * element%size(i,j) * H(i,j,ms,mt)
            
          enddo
        enddo
      enddo
    enddo

    do ms=1, n_gauss
      do mt=1, n_gauss
            
!          if (.not. ((psi_eq(ms,mt) .gt. psi_start) .and. (psi_eq(ms,mt) .lt. psi_end))) cycle
  
        wst = wgauss(ms)*wgauss(mt)
  
        xjac = x_s(ms,mt)*y_t(ms,mt) - x_t(ms,mt)*y_s(ms,mt)
        BigR = x_g(ms,mt)

        sum_particles = sum_particles + BigR * xjac * wst * f_density(1, psi_eq(ms,mt), (/0.d0, 0.d0,0.d0/))    
        sum_volume    = sum_volume    + BigR * xjac * wst   
  
      enddo
    enddo

  enddo    ! elements loop
  !$omp end parallel do
  
  total_volume    = sum_volume    * TWOPI
  total_particles = sum_particles * TWOPI

  !write(*,'(A,4e14.6)') 'DENSITY INTEGRAL : ',total_particles, total_volume
  return
  end
  

subroutine itg_energy(node_list,element_list,psi_start,psi_end,W_kin)
!---------------------------------------------------------------
!
!---------------------------------------------------------------
use data_structure
use gauss
use basis_at_gaussian
use phys_module

implicit none

type (type_node_list)    :: node_list
type (type_element_list) :: element_list
type (type_element)      :: element
type (type_node)         :: nodes(n_vertex_max)

real*8     :: x_g(n_gauss,n_gauss),  x_s(n_gauss,n_gauss),  x_t(n_gauss,n_gauss)
real*8     :: y_s(n_gauss,n_gauss),  y_t(n_gauss,n_gauss)
real*8     :: eq_s(n_gauss,n_gauss), eq_t(n_gauss,n_gauss), psi_eq(n_gauss,n_gauss)
integer    :: i, j, k, in, ms, mt, iv, inode, ife, n_elements, i_tor
real*8     :: W_kin(n_tor), xjac, BigR, wst
real*8     :: u0_x, u0_y, psi_start, psi_end, psi_norm, psi_axis, psi_bnd
integer, parameter :: ivar_psi=1, ivar_u=2, ivar_rho=5

W_kin = 0.d0

#ifdef __GFORTRAN__
   !$omp parallel do default(shared) & ! workaround for Error: �__vtab_mod_pcg32_rng_Pcg32_rng� not specified in enclosing �parallel�
#else
   !$omp parallel do default(none) &
#endif
   !$omp schedule(dynamic,10)                                                         &
   !$omp shared(element_list, node_list, psi_start, psi_end, H, H_s, H_t)             &
   !$omp private(ife, element, nodes, iv, inode, x_g, x_s, x_t, y_s, y_t, eq_s, eq_t, &
   !$omp         i, j, ms, mt, in, wst, bigR, xjac, u0_x, u0_y, psi_eq )              &
   !$omp reduction(+:W_kin)
do ife =1,  element_list%n_elements

  element = element_list%element(ife)

  do iv = 1, n_vertex_max
    inode     = element%vertex(iv)
    nodes(iv) = node_list%node(inode)
  enddo

  x_g(:,:)  = 0.d0; x_s(:,:) = 0.d0; x_t(:,:) = 0.d0;
  y_s(:,:)  = 0.d0; y_t(:,:) = 0.d0;
  eq_s(:,:) = 0.d0; eq_t(:,:) = 0.d0; psi_eq(:,:) = 0.d0

  do i=1,n_vertex_max
    do j=1,n_order+1
      do ms=1, n_gauss
        do mt=1, n_gauss

          x_g(ms,mt) = x_g(ms,mt) + nodes(i)%x(1,j,1) * element%size(i,j) * H(i,j,ms,mt)
          x_s(ms,mt) = x_s(ms,mt) + nodes(i)%x(1,j,1) * element%size(i,j) * H_s(i,j,ms,mt)
          x_t(ms,mt) = x_t(ms,mt) + nodes(i)%x(1,j,1) * element%size(i,j) * H_t(i,j,ms,mt)
          y_s(ms,mt) = y_s(ms,mt) + nodes(i)%x(1,j,2) * element%size(i,j) * H_s(i,j,ms,mt)
          y_t(ms,mt) = y_t(ms,mt) + nodes(i)%x(1,j,2) * element%size(i,j) * H_t(i,j,ms,mt)

          psi_eq(ms,mt) = psi_eq(ms,mt) + nodes(i)%values(1,j,ivar_psi) * element%size(i,j) * H(i,j,ms,mt)
          
        enddo
      enddo
    enddo
  enddo

  do in=2,n_tor

    eq_s(:,:) = 0.d0; eq_t(:,:) = 0.d0

    do ms=1, n_gauss
      do mt=1, n_gauss

        do i=1,n_vertex_max
          do j=1,n_order+1

            eq_s(ms,mt)  = eq_s(ms,mt)  + nodes(i)%values(in,j,ivar_u) * element%size(i,j) * H_s(i,j,ms,mt)
            eq_t(ms,mt)  = eq_t(ms,mt)  + nodes(i)%values(in,j,ivar_u) * element%size(i,j) * H_t(i,j,ms,mt)
        
          enddo
        enddo

      enddo
    enddo

    do ms=1, n_gauss
      do mt=1, n_gauss
          
        if (.not. ((psi_eq(ms,mt) .gt. psi_start) .and. (psi_eq(ms,mt) .lt. psi_end))) cycle

        wst = wgauss(ms)*wgauss(mt)

        xjac = x_s(ms,mt)*y_t(ms,mt) - x_t(ms,mt)*y_s(ms,mt)
        BigR = x_g(ms,mt)

        u0_x  = (   y_t(ms,mt) * eq_s(ms,mt) - y_s(ms,mt) * eq_t(ms,mt) ) / xjac
        u0_y  = ( - x_t(ms,mt) * eq_s(ms,mt) + x_s(ms,mt) * eq_t(ms,mt) ) / xjac

        W_kin(in) = W_kin(in) + (u0_x*u0_x   + u0_y*u0_y) * BigR**3 * xjac * wst        

      enddo
    enddo

  enddo  ! n_tor loop

enddo    ! elements loop
!$omp end parallel do

do in=1,n_tor
  if (mode(in) .ne. 0) then
    W_kin(in) = 0.5d0 * W_kin(in)
  endif
enddo

return
end


end program itg_loop
