!> ITPA TAE case using perpendicular and parallel pressures, CGL-type.
!> Very similar to tae_loop example.
program tae_loop_full_isotrop

  use particle_tracer
  use mod_particle_diagnostics
  use mpi
  use mod_atomic_elements
  use mod_particle_io
  use mod_event
  use mod_project_particles
  use mod_particle_loop
  use mod_jorek_timestepping
  use mod_random_seed
  use mod_interp, only: mode_moivre, interp_RZ, interp_0
  use mod_basisfunctions
  use nodes_elements
  use phys_module, only: tstep, restart, t_start, restart_particles
  use phys_module, only: CENTRAL_MASS, CENTRAL_DENSITY, xcase, xpoint
  use phys_module, only: n_particles, nstep_particles, nsubstep_particles, tstep_particles
  use phys_module, only: filter_perp, filter_hyper, filter_par, filter_perp_n0, filter_hyper_n0, filter_par_n0
  use phys_module, only: n_mode_families

  use constants,   only: MU_ZERO, MASS_PROTON, ATOMIC_MASS_UNIT, K_BOLTZ, EL_CHG

  use mod_particle_sputtering, only: particle_sputter, sample_fluid_particle_energy
  use mod_projection_functions, only: proj_f_combined_density, &
                                    proj_f_combined_energy, proj_f_combined_par_momentum
  use mod_edge_domain
  use mod_edge_elements
  use data_structure, only: type_bnd_element_list, type_bnd_node_list
  use equil_info
  use mod_boundary, only: boundary_from_grid


!$ use omp_lib

  implicit none

  type(event)                                       :: fieldreader, partreader, partwriter
!type(adf11_all)                                   :: adas
  type(pcg32_rng), dimension(:), allocatable        :: rng
  type(count_action)                                :: counter
  type(projection), target                          :: jorek_feedback, project_density, project_current
  type(jorek_timestep_action), target               :: jorek_stepper
  type(type_edge_domain), allocatable, dimension(:) :: edge_domains
  type(edge_elements)                               :: D_edge
  type(write_particle_diagnostics)                  :: diag

  real*8, parameter  :: binding_energy = 2.18d-18 ! ionization energy of a hydrogen atom [J] (= 13.6 eV)
  real*8    :: target_time
  real*8    :: physical_particles, weight
  real*8    :: oldtime, step_rest_time, particle_step_time, particle_start_time, diag_time
  real*8    :: rho_norm, t_norm, v_norm, E_norm, M_norm, N_norm, tstep_si, timesteps
  real*8    :: v_kin_temp, E(3), B(3), psi, U, B_norm(3)
  real*8    :: rescale_coef, T_axis(1), E_axis, E_hot, rho_part, v2, tstart_jorek
!$ real*8 :: w0, w1, mmm(3)

  real*8 :: test_x(3), test_E(3), test_B(3), test_psi, test_U


  integer   :: n_particles_local,n_reflect,ifail
  integer   :: i, j, k, l, m, n_steps, i_elm_old, i_diagno
  integer   :: seed, i_rng, n_stream,i_tor

! Start up MPI, jorek
  call sim%initialize(num_groups=1)

  rho_part    = 1.195d19 !(corrected value to obtain density=1.441e17 (as in benchmark, for original profile with toroidal flux)

  n_particles_local = int(n_particles/sim%n_cpu)
  timesteps         = tstep_particles

! Set up the field reader for importing restart file. The imported mode is not scaled.
  fieldreader = event(read_jorek_fields_interp_linear(basename='jorek', i=-1,mode_divisor=1))
  call with(sim, fieldreader)

  write(*,*) 'main : t_start = ',t_start

  if (sim%my_id .eq. 0) call boundary_from_grid(sim%fields%node_list, sim%fields%element_list, bnd_node_list, bnd_elm_list, .false.)

  call broadcast_boundary(sim%my_id, bnd_elm_list, bnd_node_list)

  call update_equil_state(sim%my_id,sim%fields%node_list, sim%fields%element_list, bnd_elm_list, xpoint, xcase)

  n_norm   = CENTRAL_DENSITY * 1.d20                              ! (number) density normalisation
  rho_norm = CENTRAL_MASS * MASS_PROTON * n_norm                  ! rho_SI = rho_norm * rho
  t_norm   = sqrt((MU_ZERO * rho_norm))                           ! t_SI   = t_norm * t_jorek

  tstep_si  = tstep * t_norm
  n_steps   = floor(tstep_si / timesteps)
  timesteps = tstep_si / n_steps
  n_steps   = tstep_si / timesteps

  if (sim%my_id .eq.0) then
    write(*,*) ' adapt time step to be multiple of jorek time step'
    write(*,*) "tstep = ", tstep_si, n_steps, timesteps
    write(*,*) "check :", n_steps, tstep_si - n_steps*timesteps

  endif

  if (.not. restart_particles) then
    ! Set up particles
    sim%groups(1)%Z    = 1
    sim%groups(1)%mass = atomic_weights(-2) !< atomic mass units


    allocate(particle_kinetic_leapfrog::sim%groups(1)%particles(n_particles_local))



    !< Not initialised on phi planes in order to obtain a tae mode without having to step 2000 times. This does introduce
    !< a lot of noise and produces less accurate growth rates and more noisy eigenfunctions. However, it is needed
    !< to obtain a starting point for the plane-initialized simulations.
    call initialise_particles_H_mu_psi_phiplanes(sim%groups(1)%particles, sim%fields, pcg32_rng(),sim%groups(1)%mass, &
                                     uniform_space=.true., uniform_space_rej_f=f_toroidal_flux, &
                                     uniform_space_rej_vars=[1], charge = 1, T_maxwell = 4d5)

    call adjust_particle_weights(sim%groups(1)%particles, rho_part)
    if (sim%my_id .eq. 0) write(*,*) "Particle density was adjusted to:", rho_part, sim%groups(1)%particles(1:10)%weight

    select type (p => sim%groups(1)%particles)
      type is (particle_kinetic_leapfrog)

        call boris_all_initial_half_step_backwards_RZPhi(p, sim%groups(1)%mass, sim%fields, sim%time, timesteps)

    end select

  else  ! restarting particles

    if (sim%my_id .eq. 0) write(*,*) 'restarting particles: reading part_restart.h5'

    deallocate(sim%groups)
    allocate(sim%groups(0))

    partreader = event(read_action(filename='part_restart.h5'))
    call with(sim, partreader)

  endif


  jorek_feedback = new_projection(sim%fields%node_list, sim%fields%element_list, &
                                filter    = filter_perp, filter_hyper    = filter_hyper, filter_parallel    = filter_par, &
                                filter_n0 = filter_perp, filter_hyper_n0 = filter_hyper, filter_parallel_n0 = filter_par_n0, &
                                calc_integrals=.false., to_vtk=.false., to_h5 = .false., basename='projections')


!Full tensor has 6 elements due to symmetry

  allocate(jorek_feedback%rhs(n_degrees, n_vertex_max, sim%fields%element_list%n_elements, n_tor, 6))

  jorek_feedback%rhs = 0.d0

  aux_node_list => jorek_feedback%node_list



! For proper timestepping, the projections need to be defined before the jorek timestepper
  jorek_stepper = new_jorek_timestep_action(jorek_feedback%node_list)

  diag = write_particle_diagnostics(filename='diag.h5', append=.true.)

  if (restart) then
    tstart_jorek = sim%time + tstep_si
  else
    tstart_jorek = sim%time
  endif

  if (sim%my_id .eq. 0) write(*,*) 'tstart_jorek : ',tstart_jorek

  diag_time = timesteps
  events = [new_event_ptr(jorek_feedback,  start = tstart_jorek),  &
          new_event_ptr(jorek_stepper,   start = tstart_jorek), &
          event(stop_action(), start=1d12)  ]

  jorek_stepper%extra_event => events(1)


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
  if (.not. restart_particles) call with(sim, events, at=sim%time)



  partwriter = event(write_action(filename='part_restart_init.h5'))
  call with(sim, partwriter)

  do while (.not. sim%stop_now)

    target_time = next_event_at(sim, events)
    particle_start_time = (sim%time - step_rest_time)
    particle_step_time  = target_time - particle_start_time
    n_steps             = particle_step_time/timesteps
    step_rest_time      = particle_step_time - real(n_steps,8) * timesteps

    if (sim%my_id .eq. 0) then
      if (n_steps < 10) write(*,*) 'low n_steps,', n_steps
      write(*,*) 'Time difference between particles and jorek: ', step_rest_time
      write(*,*) "PARTICLE : target time         : ",target_time
      write(*,*) "PARTICLE : timesteps           : ",timesteps
      write(*,*) "PARTICLE : sim%time            : ",sim%time
      write(*,*) "PARTICLE : particle_start_time : ",particle_start_time
      write(*,*) "PARTICLE : particle_step_time  : ",particle_step_time
      write(*,*) "PARTICLE : n_steps             : ",n_steps
      write(*,*) "PARTICLE : step_rest_time      : ",step_rest_time
    endif

    call loop_particle_kinetic_local(sim, jorek_feedback, rng, timesteps, n_steps, particle_start_time)

    sim%time = target_time

    call with(sim, events, at=sim%time)

  end do



  partwriter = event(write_action(filename='part_restart.h5'))
  call with(sim, partwriter)

  call sim%finalize

contains


subroutine loop_particle_kinetic_local(sim, jorek_feedback, rng, timesteps, n_steps, particle_start_time)
  use mod_project_particles
  use mod_random_seed
  use mod_interp, only: mode_moivre
  use mod_basisfunctions
  use mod_particle_types, only: copy_particle_kinetic_leapfrog

  implicit none

  class(particle_sim), target, intent(inout)                :: sim
  type(projection), target, intent(inout)                   :: jorek_feedback
  type(count_action)                                        :: counter
  type(pcg32_rng), dimension(:), allocatable, intent(inout) :: rng
  type(particle_kinetic_leapfrog)                           :: particle_tmp

  real*8, intent(in)     :: timesteps, particle_start_time
  real*8    :: n_norm, rho_norm, t_norm, v_norm, E_norm, M_norm
  real*8    :: t, E(3), B(3), psi, U, n_e, T_e, rz_old(2), st_old(2)
  real*8    :: v_temp(3), T_eV, K_eV, v_kin_temp, B_norm(3), v
  real*8    :: R_g, Z_g, R_s, R_t, Z_s, Z_t, xjac, HZ(n_tor), HH(4,4), HH_s(4,4), HH_t(4,4)
  real*8    :: p_isotrop
!$ real*8 :: w0, w1, mmm(3)

  integer, intent(in)   :: n_steps
  integer   :: i, j, k, l, m, i_elm_old, i_elm
  integer   :: seed, i_rng, n_stream, ierr, nthreads
  integer   :: i_tor, index_lm, i_elm_temp
  integer   :: n_particles, ifail
  real*8,allocatable :: feedback_rhs(:,:,:,:,:)

!$ w0 = omp_get_wtime()

  n_norm   = CENTRAL_DENSITY * 1.d20                              ! (number) density normalisation
  rho_norm = CENTRAL_MASS * MASS_PROTON * n_norm                  ! rho_SI = rho_norm * rho
  t_norm   = sqrt((MU_ZERO * rho_norm))                           ! t_SI   = t_norm * t_jorek
  v_norm   = 1.d0 / t_norm                                        ! V_SI   = v_norm * v_jorek
  E_norm   = 1.5d0 / MU_ZERO                                      ! E_SI   = E_norm * E_jorek
  M_norm   = rho_norm * v_norm                                    ! momentum normalisation

  jorek_feedback%rhs_gather_time = jorek_feedback%rhs_gather_time + n_steps * timesteps
  allocate(feedback_rhs,source=jorek_feedback%rhs)

  jorek_feedback%rhs = 0.d0
  feedback_rhs       = 0.d0

  call with(sim, counter)

  select type (particles => sim%groups(1)%particles)
    type is (particle_kinetic_leapfrog)
#ifdef __GFORTRAN__
      !$omp parallel do default(shared) &
#else
      !$omp parallel do default(none) &
      !$omp shared(sim, particles, n_steps, timesteps, rng, particle_start_time,        &
      !$omp rho_norm, t_norm, v_norm, E_norm, M_norm, N_norm,                           &
      !$omp jorek_feedback, CENTRAL_DENSITY, CENTRAL_MASS)                              &
#endif
      !$omp private(particle_tmp, i_rng, i,j,k,l,m, t, E, B, psi, U, rz_old, st_old,    &
      !$omp i_elm_old, i_elm, n_e, T_e, p_isotrop,                                      &
      !$omp R_g, R_s, R_t, Z_g, Z_s, Z_t, xjac, HH, HH_s, HH_t, HZ, index_lm, ifail, v) &
      !$omp schedule(dynamic,10) &
      !$omp reduction(+:feedback_rhs)
      do j=1,size(particles,1)

        call copy_particle_kinetic_leapfrog(particles(j),particle_tmp)

!      i_rng = 1
!$      i_rng = omp_get_thread_num()+1

        do k=1,n_steps

          if (particle_tmp%i_elm .le. 0) exit

          t = particle_start_time + (k-1)*timesteps

          call sim%fields%calc_EBpsiU(t, particle_tmp%i_elm, particle_tmp%st, particle_tmp%x(3), E, B, psi, U)

          rz_old    = particle_tmp%x(1:2)
          st_old    = particle_tmp%st
          i_elm_old = particle_tmp%i_elm

          if (particle_tmp%i_elm .gt. 0) then
            ! Push the particle and determine its new location.
            call boris_push_cylindrical(particle_tmp, sim%groups(1)%mass, E, B, timesteps)

            call find_RZ_nearby(sim%fields%node_list, sim%fields%element_list, rz_old(1), rz_old(2), st_old(1), st_old(2), i_elm_old, &
             particle_tmp%x(1), particle_tmp%x(2), particle_tmp%st(1), particle_tmp%st(2), particle_tmp%i_elm, ifail)
            if(particle_tmp%i_elm .gt.0) then
              call basisfunctions(particle_tmp%st(1), particle_tmp%st(2), HH, HH_s, HH_t)
              call mode_moivre(particle_tmp%x(3), HZ)
              i_elm=particle_tmp%i_elm
              p_isotrop=1.d0/3.d0*(particle_tmp%v(1)**2+particle_tmp%v(2)**2+particle_tmp%v(3)**2)

              do l=1,n_vertex_max
                do m=1,n_degrees

                  index_lm = (l-1)*n_degrees + m

                  v = HH(l,m) * sim%fields%element_list%element(i_elm)%size(l,m)

                  do i_tor=1,n_tor

                    feedback_rhs(m,l,i_elm,i_tor,1) = feedback_rhs(m,l,i_elm,i_tor,1) &

                                                      + HZ(i_tor) * v * particle_tmp%weight * sim%groups(1)%mass * mass_proton &

                                                       * ( p_isotrop ) * mu_zero !Pi_RR
                    feedback_rhs(m,l,i_elm,i_tor,2) = feedback_rhs(m,l,i_elm,i_tor,2) &

                                                       + HZ(i_tor) * v * particle_tmp%weight * sim%groups(1)%mass * mass_proton &

                                                       * ( p_isotrop ) * mu_zero                       !Pi_ZZ
                    feedback_rhs(m,l,i_elm,i_tor,3) = feedback_rhs(m,l,i_elm,i_tor,3) &

                                                     + HZ(i_tor) * v * particle_tmp%weight * sim%groups(1)%mass * mass_proton &

                                                       * ( p_isotrop ) * mu_zero !Pi_phiphi
                    feedback_rhs(m,l,i_elm,i_tor,4) = feedback_rhs(m,l,i_elm,i_tor,4) &

                                                      + HZ(i_tor) * v * particle_tmp%weight * sim%groups(1)%mass * mass_proton &

                                                      * ( 0.d0 ) * mu_zero                            !Pi_RZ
                    feedback_rhs(m,l,i_elm,i_tor,5) = feedback_rhs(m,l,i_elm,i_tor,5) &

                                                      + HZ(i_tor) * v * particle_tmp%weight * sim%groups(1)%mass * mass_proton &

                                                       * ( 0.d0 ) * mu_zero !Pi_Rphi
                    feedback_rhs(m,l,i_elm,i_tor,6) = feedback_rhs(m,l,i_elm,i_tor,6) &

                                                       + HZ(i_tor) * v * particle_tmp%weight * sim%groups(1)%mass * mass_proton &

                                                       *( 0.d0 ) * mu_zero   !Pi_Zphi

                  enddo! <tor harmonic

                enddo   !< order
              enddo     !< vertex

            end if !<particle_temp%i_elm gt 0 after pushing
          end if !<particle_temp%i_elm gt 0

        end do ! steps

        call copy_particle_kinetic_leapfrog(particle_tmp, particles(j))

      end do   ! particles
      !$omp end parallel do

      if (sim%my_id .eq. 0) write(*,*) "End of the particle loop"

  end select

  jorek_feedback%rhs = feedback_rhs/n_steps

  deallocate(feedback_rhs)

end subroutine


pure function f_adapted(n, P, grad_P) result(f)
  integer, intent(in) :: n
  real*8, intent(in) :: P(n), grad_P(3,n)
  real*8 ::s, coeff(0:4)
  real*4 :: f

  coeff(0)=0.53
  coeff(1)=0.3
  coeff(2)=0.2
  coeff(3)=0.52
  coeff(4)=0.26

  s = max((P(1) - ES%Psi_axis) / ( ES%Psi_bnd - ES%Psi_axis),0.d0)

  f = coeff(3)*exp(-coeff(2)/coeff(1)*(tanh((sqrt(s)-coeff(0))/coeff(2))))

  f = (f - coeff(4)) / (1.d0 - coeff(4))

end function f_adapted

pure function f_original(n, P, grad_P) result(f)
  integer, intent(in) :: n
  real*8, intent(in) :: P(n), grad_P(3,n)
  real*8 ::s, coeff(0:3)
  real*4 :: f

  ! central densiy should be 1.44131x10^17

  coeff(0)=0.49123
  coeff(1)=0.298228
  coeff(2)=0.198739
  coeff(3)=0.521298

  s = max((P(1) - ES%Psi_axis) / ( ES%Psi_bnd - ES%Psi_axis),0.d0)

  f = coeff(3)*exp(-coeff(2)/coeff(1)*(tanh((sqrt(s)-coeff(0))/coeff(2))))

end function f_original

pure function f_toroidal_flux(n, P, grad_P) result(f)
  integer, intent(in) :: n
  real*8, intent(in) :: P(n), grad_P(3,n)
  real*8 :: s, psi_norm, coeff(0:3)
  real*4 :: f

  ! central densiy should be 1.44131x10^17

  coeff(0)=0.49123
  coeff(1)=0.298228
  coeff(2)=0.198739
  coeff(3)=0.521298

  psi_norm = max((P(1) - ES%Psi_axis) / ( ES%Psi_bnd - ES%Psi_axis),0.d0)

  s = 0.957 * psi_norm + 0.043 * psi_norm**2

  f = coeff(3)*exp(-coeff(2)/coeff(1)*(tanh((sqrt(s)-coeff(0))/coeff(2))))

end function f_toroidal_flux

end program tae_loop_full_isotrop
