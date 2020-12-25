!> Testing the coupling of the projections of particles to JOREK

program marker_502
use particle_tracer
use mod_particle_diagnostics
use mpi
use mod_atomic_elements
use mod_particle_io
use mod_event
use mod_project_particles
use mod_jorek_timestepping
use mod_random_seed
use mod_interp, only: mode_moivre, interp_RZ
use mod_basisfunctions
use nodes_elements
use phys_module, only: n_particles, nstep_particles, nsubstep_particles, tstep_particles, use_ncs, use_pcs, use_ccs
use phys_module, only: filter_perp, filter_hyper, filter_par, filter_perp_n0, filter_hyper_n0, filter_par_n0
use phys_module, only: tstep, gas_type, imp_adas, imp_cor, use_marker
use phys_module, only: CENTRAL_MASS, CENTRAL_DENSITY
use constants,   only: MU_ZERO, MASS_PROTON, ATOMIC_MASS_UNIT, K_BOLTZ, EL_CHG

use mod_particle_sputtering, only: particle_sputter, sample_fluid_particle_energy
use mod_projection_functions, only: proj_f_combined_density, &
                                    proj_f_combined_energy, proj_f_combined_par_momentum
use mod_edge_domain
use mod_edge_elements
use mod_impurity, only: init_imp_adas
!$ use omp_lib

implicit none

type(event)                                       :: fieldreader
type(pcg32_rng), dimension(:), allocatable        :: rng
type(count_action)                                :: counter
type(projection), target                          :: jorek_feedback, project_density
type(jorek_timestep_action), target               :: jorek_stepper
type(particle_sputter)                            :: D_sputter_source
type(type_edge_domain), allocatable, dimension(:) :: edge_domains
type(edge_elements)                               :: D_edge


real*8    :: timesteps, tstep_si, t_norm, rho_norm, n_norm
real*8    :: target_time, t, E(3), B(3), psi, U, n_rho, T_e, rz_old(2), st_old(2)
real*8    :: diag_time 
real*8    :: temp(3), T_eV, K_eV, v_kin_temp, B_norm(3)
real*8    :: physical_particles, weight
integer   :: n_particles_local, n_steps, ifail
integer   :: n_reflect
integer   :: j, seed, i_rng, n_stream

! For live updating the rhs of the projection
real*8  :: R_g, Z_g, R_s, R_t, Z_s, Z_t, xjac, HZ(n_tor), HH(4,4), HH_s(4,4), HH_t(4,4)
integer :: i_tor, index_lm, i_elm_temp
logical :: use_cx, use_sputtering

! Start up MPI, jorek
call sim%initialize(num_groups=1)

n_particles_local = int(n_particles/sim%n_cpu) 
timesteps         = tstep_particles

use_cx         = .false.
use_sputtering = .false.

! --- Read ADAS data and generate coronal equilibrium is needed
call init_imp_adas(sim%my_id)

! Set up the field reader
fieldreader = event(read_jorek_fields_interp_linear(basename='jorek', i=-1))
call with(sim, fieldreader)

if (use_sputtering) then  
  n_reflect = int(n_particles * 2.d-3)
  D_sputter_source = initialise_sputtering(sim%fields%node_list, sim%fields%element_list, n_reflect)
endif

n_norm   = CENTRAL_DENSITY * 1.d20                              ! (number) density normalisation
rho_norm = CENTRAL_MASS * MASS_PROTON * n_norm                  ! rho_SI = rho_norm * rho
t_norm   = sqrt((MU_ZERO * rho_norm))                           ! t_SI   = t_norm * t_jorek

tstep_si  = tstep * t_norm
n_steps   = floor(tstep_si / timesteps)
timesteps = tstep_si / n_steps  ! Timestep in SI unit!
n_steps   = tstep_si / timesteps

if (sim%my_id .eq.0) then
  write(*,*) ' adapt time step to be multiple of jorek time step'
  write(*,*) ' tstep = ', tstep_si, n_steps, timesteps
  write(*,*) ' check : ', n_steps, tstep_si - n_steps*timesteps
endif

! Set up particles

select case ( trim(gas_type) )
  case('D2')
    sim%groups(1)%Z    = -2
    sim%groups(1)%mass = atomic_weights(-2) !< atomic mass units
    sim%groups(1)%ad   = imp_adas(1)
    sim%groups(1)%cor  = imp_cor(1)
  case('Ar')
    sim%groups(1)%Z    = 18
    sim%groups(1)%mass = atomic_weights(18) !< atomic mass units
    sim%groups(1)%ad   = imp_adas(1)
    sim%groups(1)%cor  = imp_cor(1)
  case('Ne')
    sim%groups(1)%Z    = 10
    sim%groups(1)%mass = atomic_weights(10) !< atomic mass units
    sim%groups(1)%ad   = imp_adas(1)
    sim%groups(1)%cor  = imp_cor(1)
  case default
    write(*,*) '!! Gas type "', trim(gas_type), '" unknown (in marker_502) !!'
    write(*,*) 'Exiting NOW!!!'
    stop
end select

allocate(particle_kinetic_leapfrog::sim%groups(1)%particles(n_particles_local))

call initialise_particles_H_mu_psi(sim%groups(1)%particles, sim%fields, pcg32_rng(), sim%groups(1)%mass, &
           uniform_space=.true., uniform_space_rej_f=f_psi_inside, uniform_space_rej_vars=[1], charge = 0)

physical_particles = 1.d21
weight = physical_particles/n_particles

v_kin_temp = sqrt( (2.d0 * 1d5) / (sim%groups(1)%mass* ATOMIC_MASS_UNIT) / physical_particles)

select type (p => sim%groups(1)%particles)
type is (particle_kinetic_leapfrog)
 
  p(:)%q      = 0
  p(:)%weight = weight

  do j=1,size(p,1)
    call sim%fields%calc_EBpsiU(sim%time , p(j)%i_elm, p(j)%st, p(j)%x(3), E, B, psi, U)
    B_norm = B/norm2(B)
    p(j)%v(1)  = v_kin_temp * B_norm(1)
    p(j)%v(2)  = v_kin_temp * B_norm(2)
    p(j)%v(3)  = v_kin_temp * B_norm(3)
  end do
  call boris_all_initial_half_step_backwards_RZPhi(p, sim%groups(1)%mass, sim%fields, sim%time, timesteps)
end select

! Set up feedback
jorek_feedback = new_projection(sim%fields%node_list, sim%fields%element_list, &
                     filter    = filter_perp,    filter_hyper    = filter_hyper,    filter_parallel    = filter_par, &
                     filter_n0 = filter_perp_n0, filter_hyper_n0 = filter_hyper_n0, filter_parallel_n0 = filter_par_n0, &
                     fractional_digits = 9,  to_vtk=.FALSE., to_h5 = .FALSE., basename='projections')

aux_node_list => jorek_feedback%node_list

if (use_ncs) then
  allocate(jorek_feedback%rhs(n_order+1, n_vertex_max, sim%fields%element_list%n_elements, n_tor, 3))
elseif (use_pcs) then  ! not implemented yet!
  allocate(jorek_feedback%rhs(n_order+1, n_vertex_max, sim%fields%element_list%n_elements, n_tor, 1))
elseif (use_ccs) then  ! not implemented yet!
  allocate(jorek_feedback%rhs(n_order+1, n_vertex_max, sim%fields%element_list%n_elements, n_tor, 4))
elseif (use_marker) then
  allocate(jorek_feedback%rhs(n_order+1, n_vertex_max, sim%fields%element_list%n_elements, n_tor, 5))
else
  stop 'define use_ncs, use_pcs, use_marker or use_ccs'
endif

jorek_feedback%rhs = 0.d0

!project_density = new_projection(sim%fields%node_list, sim%fields%element_list, &
!                      filter    = filter_perp,    filter_hyper    = filter_hyper,    filter_parallel    = filter_par, &
!                      filter_n0 = filter_perp_n0, filter_hyper_n0 = filter_hyper_n0, filter_parallel_n0 = filter_par_n0, &
!                      f=[proj_f(proj_one, group = 1)], &
!                      fractional_digits = 9,  to_vtk=.TRUE., to_h5=.FALSE., basename='density', nsub=5)
!
!call with(sim, project_density)

! For proper timestepping, the projections need to be defined before the jorek timestepper
jorek_stepper = new_jorek_timestep_action(jorek_feedback%node_list)

diag_time = 1.0d12
events = [ new_event_ptr(jorek_feedback,   start = sim%time),            &
           new_event_ptr(jorek_stepper,    start = sim%time),            &
!          new_event_ptr(D_sputter_source, start = sim%time, step=5d-6), &
!          event(count_action(),           start = sim%time, step=1d-5), &
!          event(write_particle_diagnostics(filename='diag.h5'), step=diag_time), &
!          event(write_action(), step=diag_time),                        &
!           new_event_ptr(project_density, step=1.d-5),                  &
           event(stop_action(), start=1d12)                              &
        ]

jorek_stepper%extra_event => events(1)

call main_particle_loop(jorek_stepper, jorek_feedback, project_density, timesteps, &
                        use_cx, use_sputtering)

call sim%finalize

contains

!================================================================================================
subroutine main_particle_loop(jorek_stepper, jorek_feedback, project_density, timesteps, &
                              use_cx, use_sputtering)
!================================================================================================
use particle_tracer
use mod_particle_diagnostics
use mpi
use mod_atomic_elements
use mod_particle_io
use mod_event
use mod_project_particles
use mod_random_seed
use mod_interp, only: mode_moivre, interp_RZ
use mod_jorek_timestepping
use mod_basisfunctions
use phys_module, only: tstep, use_ncs, use_pcs, use_ccs, use_marker
use phys_module, only: CENTRAL_MASS, CENTRAL_DENSITY, GAMMA
use constants,   only: MU_ZERO, MASS_PROTON, ATOMIC_MASS_UNIT, K_BOLTZ, EL_CHG
use mod_integrals3D, only: int3d_new
use mod_radiation, only: proj_Lz

implicit none
real*8, parameter  :: binding_energy = 2.18d-18 ! ionization energy of a hydrogen atom [J] (= 13.6 eV)

real*8, intent(in) :: timesteps
logical            :: use_cx, use_sputtering

type(projection), target                          :: jorek_feedback, project_density
type(jorek_timestep_action), target               :: jorek_stepper

real*8,allocatable :: feedback_rhs(:,:,:,:,:)
real*8    :: oldtime, step_rest_time, particle_step_time, particle_start_time, diag_time
real*8    :: rho_norm, t_norm, v_norm, E_norm, M_norm, N_norm, tstep_si
real*8    :: kinetic_energy, ion_energy
real*8    :: E_lost_ion, E_lost_ion_all, E_lost_rad, E_lost_rad_all
!$ real*8 :: w0, w1, mmm(3)
integer   :: i, j, k, l, m, n_steps, i_elm_old
integer   :: seed, i_rng, n_stream, ierr, nthreads, myid
real*8    :: ion_rate, ion_source, ion_prob, ion_rec_ran(2), cx_ran(7), cx_source, cx_energy
real*8    :: rec_rate, dEion_dT, Z_imp, Z_eff, N_imp, Lrad, rad_sink
real*8    :: cx_prob, CX_rate
real*8    :: particle_source, velocity_par_source, energy_source
real*8    :: v_temp(3), T_eV, K_eV, v_kin_temp, B_norm(3), v_1, v_2, v_3, v_4, v_5
real*8    :: density_tot, density_in, density_out,  pressure, pressure_in, pressure_out
real*8    :: mom_par_tot, mom_par_in, mom_par_out, kin_par_tot, kin_par_out, kin_par_in
real*8    :: particles_remaining, momentum_remaining, energy_remaining, all_particles, all_momentum, all_energy

n_norm   = CENTRAL_DENSITY * 1.d20                              ! (number) density normalisation
rho_norm = CENTRAL_MASS * MASS_PROTON * n_norm                  ! rho_SI = rho_norm * rho
t_norm   = sqrt((MU_ZERO * rho_norm))                           ! t_SI   = t_norm * t_jorek
v_norm   = 1.d0 / t_norm                                        ! V_SI   = v_norm * v_jorek
E_norm   = 1.d0 / (MU_ZERO * (GAMMA-1.))                       ! E_SI   = E_norm * E_jorek
M_norm   = rho_norm * v_norm                                    ! momentum normalisation

particle_source = 0.0
velocity_par_source = 0.0
energy_source = 0.0
dEion_dT = 0.0
Z_imp    = 0.0; Z_eff = 0.0; N_imp = 0.0

v_1 = 0.0; v_2 = 0.0; v_3 = 0.0; v_4 = 0.0; v_5 = 0.0

if (sim%my_id .eq. 0) then
  if (use_cx) then
    write(*,*) ' including charge exchange'
  else
    write(*,*) ' NOT including charge exchange'
  endif
  write(*,*)
  write(*,'(A,e14.6)') ' N_norm   : ',N_norm
  write(*,'(A,e14.6)') ' rho_norm : ',rho_norm
  write(*,'(A,e14.6)') ' t_norm   : ',t_norm
  write(*,'(A,e14.6)') ' V_norm   : ',v_norm
  write(*,'(A,e14.6)') ' M_norm   : ',M_norm
  write(*,'(A,e14.6)') ' E_norm   : ',E_norm
  write(*,*)
  write(*,'(A,e14.6)') ' tstep    : ', tstep
endif

! Set up random numbers for ionisation probability
seed = random_seed()
n_stream = 1
!$ n_stream = omp_get_max_threads()
allocate(rng(n_stream))
do i=1,n_stream
  call rng(i)%initialize(1, seed, n_stream, i)
end do

! Call events at sim%time once to help event scheduler, before entering particle loop
step_rest_time = 0.d0

call with(sim, events, at=sim%time)

do while (.not. sim%stop_now)

  target_time = next_event_at(sim, events) 

  if (sim%my_id .eq. 0) write(*,'(A,3e16.8)') "PARTICLE : target_time : ",target_time

!$ w0 = omp_get_wtime()

  ! step_rest_time is the amount of time that needs to be covered in addition to target_time - sim%time,
  ! i.e. (sim%time - step_rest_time(i)) is the 'real' particle time before the next loop
  particle_start_time = (sim%time - step_rest_time)
  particle_step_time  = target_time - particle_start_time
  n_steps             = particle_step_time/timesteps
  step_rest_time      = particle_step_time - real(n_steps,8) * timesteps

  if (sim%my_id .eq. 0) then
    if (n_steps < 10) write(*,*) 'low n_steps,', n_steps
    write(*,*) 'Time difference between particles and jorek: ', step_rest_time

    write(*,*) "PARTICLE : sim%time            : ",sim%time
    write(*,*) "PARTICLE : particle_start_time : ",particle_start_time
    write(*,*) "PARTICLE : particle_step_time  : ",particle_step_time
    write(*,*) "PARTICLE : n_steps             : ",n_steps
    write(*,*) "PARTICLE : step_rest_time      : ",step_rest_time
  endif

  E_lost_ion = 0.d0
  E_lost_ion_all = 0.d0
  E_lost_rad = 0.d0
  E_lost_rad_all = 0.d0

!  jorek_feedback%rhs_gather_time = jorek_feedback%rhs_gather_time + n_steps * timesteps
  jorek_feedback%rhs_gather_time = n_steps * timesteps

  allocate(feedback_rhs,source=jorek_feedback%rhs)

  jorek_feedback%rhs = 0.d0
  feedback_rhs       = 0.d0
  
  call with(sim, counter)
 
  select type (particles => sim%groups(1)%particles)
  type is (particle_kinetic_leapfrog)

#ifdef __GFORTRAN__
    !$omp parallel do default(shared) & ! workaround for Error: '__vtab_mod_openadas_Adf11' not specified in enclosing ‘parallel’
#else
    !$omp parallel do default(none) &
#endif
    !$omp schedule(dynamic,10)      &
#ifdef __GFORTRAN__
    !$omp shared(sim, n_particles, n_steps, timesteps, rng, particle_start_time, & ! This is to work around the GNU compiler error: ASSOCIATE name '__tmp_type_particle_kinetic_leapfrog' in SHARED clause
#else
    !$omp shared(sim, particles, n_particles, n_steps, timesteps, rng, particle_start_time, &
#endif
    !$omp rho_norm, t_norm, v_norm, E_norm, M_norm, N_norm, &
    !$omp use_cx, use_sputtering, use_ncs, use_marker,                                      &
    !$omp CENTRAL_DENSITY, CENTRAL_MASS)                    &
    !$omp private(i_rng, i,j,k,l,m, t, E, B, psi, U, rz_old, st_old,                        &
    !$omp i_elm_old, n_rho, T_e, ion_rate, ion_prob, ion_source, ion_energy, kinetic_energy,& 
    !$omp rec_rate, dEion_dt, ion_rec_ran, Z_imp, Z_eff, N_imp, Lrad, rad_sink, & 
    !$omp R_g, R_s, R_t, Z_g, Z_s, Z_t, xjac, HH, HH_s, HH_t, HZ, index_lm,                 &
    !$omp ifail, CX_rate, CX_prob, CX_source, CX_energy, v_1, v_2, v_3, v_4, v_5,           &
    !$omp particle_source, velocity_par_source, energy_source, v_temp, K_eV, T_eV, cx_ran)  &
    !$omp reduction(+:feedback_rhs, E_lost_ion, E_lost_rad)
    do j=1,size(particles,1)

!      i_rng = 1
    !$ i_rng = omp_get_thread_num()+1

      do k=1,n_steps

        if (particles(j)%i_elm .le. 0) exit

        t = particle_start_time + (k-1)*timesteps

        call sim%fields%calc_EBpsiU(t, particles(j)%i_elm, particles(j)%st, particles(j)%x(3), E, B, psi, U)
        rz_old    = particles(j)%x(1:2)
        st_old    = particles(j)%st
        i_elm_old = particles(j)%i_elm

        call sim%fields%calc_NeTe(t, particles(j)%i_elm, particles(j)%st, particles(j)%x(3), n_rho, T_e)

        ion_source = 0.d0
        ion_energy = 0.d0

        ! Do the time average of the the impurity number density as well as average charge
        Z_imp = real(particles(j)%weight,8) * real(particles(j)%q,8) * timesteps
        Z_eff = real(particles(j)%weight,8) * real(particles(j)%q,8)**2 * timesteps
        N_imp = real(particles(j)%weight,8) * timesteps 

        ! Do the particle radiation
        ! We take away the n_rho here since it is not really the electron density
        Lrad = proj_Lz(sim ,1, particles(j)) / n_rho
        rad_sink = real(particles(j)%weight,8) * Lrad * timesteps 
        E_lost_rad = E_lost_rad + real(particles(j)%weight,8) * Lrad * timesteps * n_rho

        ! Do the ionization and recombination
        if (particles(j)%q .eq. 0) then
          call sim%groups(1)%ad%SCD%interp_linear(int(particles(j)%q), log10(n_rho), log10(T_e), ion_rate) ! [m^3/s]
          rec_rate = 0.
        elseif(particles(j)%q .eq. sim%groups(1)%ad%n_Z) then
          call sim%groups(1)%ad%ACD%interp_linear(int(particles(j)%q), log10(n_rho), log10(T_e), rec_rate) ! [m^3/s]
          ion_rate = 0.
        else             
          call sim%groups(1)%ad%SCD%interp_linear(int(particles(j)%q), log10(n_rho), log10(T_e), ion_rate) ! [m^3/s]
          call sim%groups(1)%ad%ACD%interp_linear(int(particles(j)%q), log10(n_rho), log10(T_e), rec_rate) ! [m^3/s]
        endif

        ion_prob = 1.d0 - exp(-(ion_rate+rec_rate) * n_rho * timesteps) ! [0] poisson point process, exponential 

        ! If the weight is small ionize the particle with the probability, else stop the code

        if (particles(j)%weight .le. 1.0d7) then

          call rng(i_rng)%next(ion_rec_ran)

          ! Take ionization energy from the electron when ionizing, give the ionization energy back to be
          ! cancelled with the recombination radiation power density to avoid double counting.
          if (ion_rec_ran(1) .le. ion_prob) then
            if ((ion_rec_ran(2) .le. (ion_rate/(ion_rate+rec_rate))) .and. (ion_rec_ran(2) .ne. 0.d0)) then
              particles(j)%q  = particles(j)%q + 1
              ion_source = real(particles(j)%weight,8)
              dEion_dt = ion_source * sim%groups(1)%ad%ionisation_energy(int(particles(j)%q))
            else
              particles(j)%q  = particles(j)%q - 1
              ion_source = -real(particles(j)%weight,8)
              dEion_dt = ion_source * sim%groups(1)%ad%ionisation_energy(int(particles(j)%q+1))
            endif
          endif
          
          dEion_dt = dEion_dt * EL_CHG ! Convert from eV to Joule
          E_lost_ion = E_lost_ion + dEion_dt ! Accumulated ionization energy

        else 
          write(*,*) "The particle weight is too large! Not supported for now, EXITING!!!"
          stop
        endif 

        kinetic_energy = dot_product(particles(j)%v,particles(j)%v) *sim%groups(1)%mass * ATOMIC_MASS_UNIT /2.d0

        ion_energy     = kinetic_energy !- binding_energy


        ! Charge Exchange
        ! It is assumed that we will have a exchange between hydrogen isotopes

        v_temp    = particles(j)%v
        cx_source = 0.d0
        cx_energy = 0.d0

        if (use_cx) then
  
          call sim%groups(1)%ad%CCD%interp(int(particles(j)%q+1), log10(n_rho), log10(T_e), CX_rate) ! [m^3/s]

          CX_prob = 1.d0 - exp(-CX_rate * n_rho * timesteps)

          call rng(i_rng)%next(cx_ran)
 
          if (cx_ran(1) .le. CX_prob) then

            ! sample boltzman, randomize velocity
            T_eV = T_e * K_BOLTZ / EL_CHG

            ! sample from main plasma (should this not be a shifted Maxwellian?)

            call sample_fluid_particle_energy(T_eV, cx_ran(2:4), 1, K_eV) ! K_eV in eV. 

!THIS IS WRONG: use sample distorted maxwellian (or box-Mueller transform)
            v_temp    = sqrt(2.d0* K_eV *EL_CHG/(sim%groups(1)%mass * ATOMIC_MASS_UNIT)) * cx_ran(5:7)

            CX_source = particles(j)%weight

            CX_energy   = 0.5d0 * sim%groups(1)%mass * ATOMIC_MASS_UNIT *  (dot_product(particles(j)%v,particles(j)%v) - dot_product(v_temp,v_temp))

          endif

        endif

        if (use_ncs) then
          energy_source       = ion_source * ion_energy + cx_source * cx_energy

          particle_source     = ion_source * sim%groups(1)%mass * ATOMIC_MASS_UNIT
       
          velocity_par_source = ion_source * dot_product(B, particles(j)%v)          * sim%groups(1)%mass * ATOMIC_MASS_UNIT &
                            
                              + CX_source  * dot_product(B, particles(j)%v - v_temp) * sim%groups(1)%mass * ATOMIC_MASS_UNIT 
        endif        
      
        particles(j)%v = v_temp 

        ! Calculate the projection of the ion source in real-time

        call basisfunctions(particles(j)%st(1), particles(j)%st(2), HH, HH_s, HH_t)
        call mode_moivre(particles(j)%x(3), HZ)
              
        do l=1,n_vertex_max
          do m=1,n_order+1

            index_lm = (l-1)*(n_order+1) + m

            if (use_ncs) then
              v_1 = HH(l,m) * sim%fields%element_list%element(i_elm_old)%size(l,m) * particle_source     * t_norm / rho_norm
              v_2 = HH(l,m) * sim%fields%element_list%element(i_elm_old)%size(l,m) * energy_source       * t_norm / E_norm
              v_3 = HH(l,m) * sim%fields%element_list%element(i_elm_old)%size(l,m) * velocity_par_source * t_norm / m_norm
            else if (use_marker) then
              v_1 = HH(l,m) * sim%fields%element_list%element(i_elm_old)%size(l,m) * dEion_dt * t_norm / E_norm
              v_2 = HH(l,m) * sim%fields%element_list%element(i_elm_old)%size(l,m) * rad_sink * t_norm * n_norm / E_norm
              v_3 = HH(l,m) * sim%fields%element_list%element(i_elm_old)%size(l,m) * Z_eff / n_norm
              v_4 = HH(l,m) * sim%fields%element_list%element(i_elm_old)%size(l,m) * Z_imp / n_norm
              v_5 = HH(l,m) * sim%fields%element_list%element(i_elm_old)%size(l,m) * N_imp / n_norm
            endif

            do i_tor=1,n_tor
              feedback_rhs(m,l,i_elm_old,i_tor,1) = feedback_rhs(m,l,i_elm_old,i_tor,1) + HZ(i_tor) * v_1
              feedback_rhs(m,l,i_elm_old,i_tor,2) = feedback_rhs(m,l,i_elm_old,i_tor,2) + HZ(i_tor) * v_2
              feedback_rhs(m,l,i_elm_old,i_tor,3) = feedback_rhs(m,l,i_elm_old,i_tor,3) + HZ(i_tor) * v_3
              feedback_rhs(m,l,i_elm_old,i_tor,4) = feedback_rhs(m,l,i_elm_old,i_tor,4) + HZ(i_tor) * v_4
              feedback_rhs(m,l,i_elm_old,i_tor,5) = feedback_rhs(m,l,i_elm_old,i_tor,5) + HZ(i_tor) * v_5
            enddo

          enddo
        enddo

        if (particles(j)%i_elm .gt. 0) then
          ! Push the particle and determine it's new location.
          call boris_push_cylindrical(particles(j), sim%groups(1)%mass, E, B, timesteps)

          call find_RZ_nearby(sim%fields%node_list, sim%fields%element_list, rz_old(1), rz_old(2), st_old(1), st_old(2), i_elm_old, &
                              particles(j)%x(1), particles(j)%x(2), particles(j)%st(1), particles(j)%st(2), particles(j)%i_elm, ifail)
        end if

      end do ! steps
    end do   ! particles
    !$omp end parallel do

  end select

  if (use_ncs) then
    write(*,*) 'GATHER TIME : ',jorek_feedback%rhs_gather_time
    jorek_feedback%rhs = feedback_rhs / jorek_feedback%rhs_gather_time !* TWOPI
    jorek_feedback%rhs_gather_time = 0.d0
  else if (use_marker) then
    write(*,*) 'GATHER TIME : ',jorek_feedback%rhs_gather_time
    jorek_feedback%rhs(:,:,:,:,1) = feedback_rhs(:,:,:,:,1) / jorek_feedback%rhs_gather_time
    jorek_feedback%rhs(:,:,:,:,2) = feedback_rhs(:,:,:,:,2) / jorek_feedback%rhs_gather_time
    jorek_feedback%rhs(:,:,:,:,3) = feedback_rhs(:,:,:,:,3) / (feedback_rhs(:,:,:,:,5) + jorek_feedback%rhs_gather_time) !  In case of feedback_rhs(:,:,:,:,4) == 0
    jorek_feedback%rhs(:,:,:,:,4) = feedback_rhs(:,:,:,:,4) / (feedback_rhs(:,:,:,:,5) + jorek_feedback%rhs_gather_time) !  In case of feedback_rhs(:,:,:,:,4) == 0
    jorek_feedback%rhs(:,:,:,:,5) = feedback_rhs(:,:,:,:,5) / jorek_feedback%rhs_gather_time
    jorek_feedback%rhs_gather_time = 0.d0
  else
    jorek_feedback%rhs = feedback_rhs 
  endif

  deallocate(feedback_rhs)

!  if (minval(jorek_feedback%rhs(:,:,:,:,2)) < 0.0) then
!    write(*,*) "SOMETHING WRONG in rad feedbacks,", minval(jorek_feedback%rhs(:,:,:,:,2))
!    stop
!  endif
!
!  if ((maxval(jorek_feedback%rhs(:,:,:,:,4)) > real(sim%groups(1)%ad%n_Z,8)) .or. &
!        (minval(jorek_feedback%rhs(:,:,:,:,4)) < 0.0)) then
!    write(*,*) "SOMETHING WRONG in mean charge feedbacks,", maxval(jorek_feedback%rhs(:,:,:,:,3)), minval(jorek_feedback%rhs(:,:,:,:,3))
!    stop
!  endif
!
!  if (minval(jorek_feedback%rhs(:,:,:,:,3)) < 0.0) then
!    write(*,*) "SOMETHING WRONG in effective charge feedbacks,", minval(jorek_feedback%rhs(:,:,:,:,4))
!    stop
!  endif
!
!  if (minval(jorek_feedback%rhs(:,:,:,:,5)) < 0.0) then
!    write(*,*) "SOMETHING WRONG in density feedbacks,", minval(jorek_feedback%rhs(:,:,:,:,4))
!    stop
!  endif

  call MPI_AllReduce(E_lost_ion,E_lost_ion_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_AllReduce(E_lost_rad,E_lost_rad_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  
  if (sim%my_id .eq. 0) write(*,*) " Lost energy at t due to ionisation: ", sim%time, E_lost_ion_all
  if (sim%my_id .eq. 0) write(*,*) " Lost energy at t due to radiation: ",  sim%time, E_lost_rad_all
  !$ w1 = omp_get_wtime()
  !$ mmm = mpi_minmeanmax(w1-w0)
  !$ if (sim%my_id .eq. 0) write(*,"(f10.7,A,3f9.4,A)") sim%time, " Particle stepping complete in ", mmm, "s"

!===================================================  
  sim%time = target_time 

  call with(sim, events, at=sim%time)
!===================================================

!  call Integrals_3D(sim%my_id, sim%fields%node_list, sim%fields%element_list, density_tot, density_in, density_out, &
!                    pressure, pressure_in, pressure_out, kin_par_tot, kin_par_in, kin_par_out, mom_par_tot, mom_par_in, mom_par_out)

  particles_remaining = 0.d0
  momentum_remaining  = 0.d0
  energy_remaining    = 0.d0

  select type (particles => sim%groups(1)%particles)
  type is (particle_kinetic_leapfrog)

#ifdef __GFORTRAN__
    !$omp parallel do default(shared) & ! workaround for Error: '__vtab_mod_openadas_Adf11' not specified in enclosing 'parallel'
#else
    !$omp parallel do default(none) & 
#endif
    !$omp reduction(+:particles_remaining, momentum_remaining, energy_remaining) &
    !$omp shared(sim) &
    !$omp private(j, E, B, psi, U, B_norm)
    do j=1,size(particles,1)

      if (particles(j)%i_elm .le. 0) cycle

      call sim%fields%calc_EBpsiU(sim%time , particles(j)%i_elm, particles(j)%st, particles(j)%x(3), E, B, psi, U)
      B_norm = B/norm2(B)

      particles_remaining = particles_remaining + particles(j)%weight
      momentum_remaining  = momentum_remaining  + particles(j)%weight * dot_product(B_norm,particles(j)%v) *sim%groups(1)%mass * ATOMIC_MASS_UNIT
      energy_remaining    = energy_remaining    + particles(j)%weight * dot_product(particles(j)%v,particles(j)%v) *sim%groups(1)%mass * ATOMIC_MASS_UNIT /2.d0
!      energy_remaining    = energy_remaining    + particles(j)%weight * 2.18d-15

    enddo
    !omp end parallel do
  end select

  call MPI_REDUCE(particles_remaining, all_particles, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_REDUCE(momentum_remaining,  all_momentum,  1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_REDUCE(energy_remaining,    all_energy,    1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  
  if (sim%my_id .eq. 0) then
    write(*,'(A,3e16.8)') 'REMAINING (START) : ',all_particles, all_momentum, all_energy

    write(*,'(A,126e16.8)') ' TOTAL : ',sim%time,density_tot+all_particles/1.d20, density_tot, all_particles/1.d20, &
                                        mom_par_tot+all_momentum, mom_par_tot, all_momentum, &
                                        pressure+kin_par_tot+all_energy, pressure, all_energy, kin_par_tot
  endif

end do


end subroutine

function initialise_sputtering(node_list, element_list, n_reflect) result(D_sputter_source)

  use mod_edge_domain
  use mod_edge_elements

  type(type_node_list), intent(in)    :: node_list
  type(type_element_list)             :: element_list
  type(particle_sputter)              :: D_sputter_source
  integer                             :: n_reflect
  type(type_edge_domain), allocatable, dimension(:) :: edge_domains

  ! number of particles to sputter per species (should be renormalized to yield)

  call find_edge_domains(node_list,element_list, edge_domains)

  call D_edge%prepare(node_list, element_list, edge_domains, nsub=6, nsub_toroidal=1)

  ! target group, number of particles per mpi task, densities, Zs, basename
  D_sputter_source = particle_sputter(D_edge, 1, n_reflect, basename='D_reflect')
  D_sputter_source%use_Yn_func = .false.

end function


pure function f_psi_inside(n, P, grad_P) result(f)
  integer, intent(in) :: n
  real*8, intent(in) :: P(n), grad_P(3,n)
  real*4 :: f, psi_norm

  psi_norm = (p(1) + 0.4463)/0.4463
 
  f = max((1. - psi_norm/0.2), 0.e0)**2

end function f_psi_inside

end program marker_502
