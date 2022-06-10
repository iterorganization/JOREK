! Example for using phase space projections
program tae_phase_space_project

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
  use phys_module, only: CENTRAL_MASS, CENTRAL_DENSITY, xcase, xpoint,F0
  use phys_module, only: n_particles, nstep_particles, nsubstep_particles, tstep_particles
  use phys_module, only: filter_perp, filter_hyper, filter_par, filter_perp_n0, filter_hyper_n0, filter_par_n0
  use phys_module, only: n_mode_families, nout

  use constants,   only: MU_ZERO, MASS_PROTON, ATOMIC_MASS_UNIT, K_BOLTZ, EL_CHG
  use mod_math_operators, only: cross_product
  use mod_particle_sputtering, only: particle_sputter, sample_fluid_particle_energy
  use mod_projection_functions, only: proj_f_combined_density, &
      proj_f_combined_energy, proj_f_combined_par_momentum
  use mod_edge_domain
  use mod_edge_elements
  use data_structure, only: type_bnd_element_list, type_bnd_node_list
  use equil_info
  use mod_boundary, only: boundary_from_grid
  use mod_phase_space_project

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
  type(phase_space_projection)                      :: fourD_dist
  type(phase_space_projection)                      :: power_exchange_vpar_mu
  real*8, parameter  :: binding_energy = 2.18d-18 ! ionization energy of a hydrogen atom [J] (= 13.6 eV)
  real*8    :: target_time
  real*8    :: physical_particles, weight
  real*8    :: oldtime, step_rest_time, particle_step_time, particle_start_time, diag_time
  real*8    :: rho_norm, t_norm, v_norm, E_norm, M_norm, N_norm, tstep_si, timesteps
  real*8    :: v_kin_temp, E(3), B(3), psi, U, B_norm(3)
  real*8    :: rescale_coef, T_axis(1), E_axis, E_hot, rho_part, v2, tstart_jorek
!$ real*8 :: w0, w1, mmm(3)

  real*8 :: test_x(3), test_E(3), test_B(3), test_psi, test_U
  !real*8 :: b_norm_r, b_norm_z, b_norm_phi

  integer   :: n_particles_local,n_reflect,ifail, ino
  integer   :: i, j, k, l, m, n_steps, i_elm_old, i_diagno
  integer   :: seed, i_rng, n_stream,i_tor,ierr

  ! Start up MPI, jorek
  call sim%initialize(num_groups=1)

  rho_part    = 1.195d19 !(corrected value to obtain density=1.441e17 (as in benchmark, for original profile with toroidal flux)
  n_particles_local = int(n_particles/sim%n_cpu)
  timesteps         = tstep_particles

  ! Set up the field reader
  fieldreader = event(read_jorek_fields_interp_linear(basename='jorek', i=-1,mode_divisor=100))
  call with(sim, fieldreader)

  if(sim%my_id .eq. 0 ) write(*,*) 'main : t_start = ',t_start


  if (sim%my_id .eq. 0) call boundary_from_grid(sim%fields%node_list, sim%fields%element_list, bnd_node_list, bnd_elm_list, .false.)

  call broadcast_boundary(sim%my_id, bnd_elm_list, bnd_node_list)

  call update_equil_state(sim%my_id,sim%fields%node_list, sim%fields%element_list, bnd_elm_list, xpoint, xcase)

  n_norm   = CENTRAL_DENSITY * 1.d20                              ! (number) density normalisation
  rho_norm = CENTRAL_MASS * MASS_PROTON * n_norm                  ! rho_SI = rho_norm * rho
  t_norm   = sqrt((MU_ZERO * rho_norm))                           ! t_SI   = t_norm * t_jorek
  write(*,*) 'Alfven =', F0/10.d0/t_norm
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
    !sim%groups(1)%ad   = adas

    allocate(particle_kinetic_leapfrog::sim%groups(1)%particles(n_particles_local))

    !< If projecting phase space, it is vital to use a by construction phi independent initial distribution, i.e. phi planes.
    !< This lowers the sampling in the other coordinates naturally, so for the example it's not done 
    !< Example of usage is here in any case
    !call initialise_particles_H_mu_psi_phiplanes(sim%groups(1)%particles, sim%fields, pcg32_rng(),sim%groups(1)%mass, &
    !     uniform_space=.true., uniform_space_rej_f=f_toroidal_flux, &
    !     uniform_space_rej_vars=[1], charge = 1, T_maxwell = 4d5)   
    call initialise_particles_H_mu_psi_phiplanes(sim%groups(1)%particles, sim%fields, pcg32_rng(),sim%groups(1)%mass, &
         uniform_space=.true., uniform_space_rej_f=f_toroidal_flux, &
         uniform_space_rej_vars=[1], charge = 1, T_maxwell = 4d5,n_phi_planes_in=12)   

    call adjust_particle_weights(sim%groups(1)%particles, rho_part)
    if (sim%my_id .eq. 0) write(*,*) "Particle density was adjusted to:", rho_part, sim%groups(1)%particles(1)%weight
    diag = write_particle_diagnostics(filename='diag.h5', append=.false.)
    call with(sim,diag)
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

  ! Be aware of OMP_STACKSIZE limits here as well. Each OMP process will create a private copy of the array on the stack (for reduction later), so for large arrays
  ! such as this 4D array, it will be expensive. 90*95*50*75*8 ~ 256.5e6 = 256.5 MegaByte. If you run 
  !    export OMP_STACKSIZE=256M 
  ! the OMP stack size of each worker will be 256 MebiByte = 268 MegaByte, which barely fits this array. (also, be aware that intel OMP can overwrite this with KMP_STACKSIZE)
  ! The default on MPCDF machines is 256M, so it will work there without further action.
  ! Comment out these three lines if it is not needed.

  fourD_dist = new_phase_space_projection(ndim=4,res=[90,95,50,75],start=[9.d0,-1.d0,-1.1d0, -20.d0],end=[11.0,1.d0,1.1d0,1700.d0], f_proj = proj_f(proj_one,1), f_grids=proj_ndim_f(f=proj_RZPE, group=1),basename='fourD_dist')
  call with(sim,fourD_dist)
  call output_phase_project(fourD_dist,0,output_grids_in=.true.)

  ! The output is only done on the root process. To prevent a too big imbalance, barrier here.
  call MPI_BARRIER(MPI_COMM_WORLD,ifail)


  ! Example for power exchange with custom bandwidths in particle loop. Will be more useful with n_phi_planes initialisation and a restarted mode structure.
  ! However, then during the first timesteps the energy difference is so small floating point limitations limit the accuracy of the diagnostic. It does immediately
  ! work properly if you don't use n_phi_planes which you can use to test the diagnostic script.
  power_exchange_vpar_mu = new_phase_space_projection(ndim=2,res=[300,300],start=[-2.d7,-0.1d6],end=[2.d7,1.5d6],basename="power_exchange",bandwidths=[0.15d7,0.08d6])
  
  ! Full tensor + density (for density flattening was the idea)
  allocate(jorek_feedback%rhs(n_order+1, n_vertex_max, sim%fields%element_list%n_elements, n_tor, 7))

  jorek_feedback%rhs = 0.d0

  aux_node_list => jorek_feedback%node_list

  ! For proper timestepping, the projections need to be defined before the jorek timestepper
  jorek_stepper = new_jorek_timestep_action(jorek_feedback%node_list)

  if (restart) then
    tstart_jorek = sim%time + tstep_si
  else
    tstart_jorek = sim%time
  endif

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
  ino = 0
  if (.not. restart_particles) call with(sim, events, at=sim%time)

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
    power_exchange_vpar_mu%values = 0.d0
    call loop_particle_kinetic_local(sim, jorek_feedback, rng, timesteps, n_steps , particle_start_time,power_exchange_vpar_mu)
    ! Will output every time step. A comulative sum of these is needed, but for now this allows you to do anything with the resulting files.
    ! You can also move the power_exchange_vpar_mu%values=0.d0 line outside the loop, at which point it will do the cumulative sum itself.
    if(ino .eq.0) then
       call output_phase_project(power_exchange_vpar_mu,ino,output_grids_in=.true.)
    else
       call output_phase_project(power_exchange_vpar_mu,ino,output_grids_in=.false.)
    endif
    ! Output a 4D distribution function each nout so we can analyuze it later
    if(mod(ino,nout) .eq. 0) then 
       call with(sim,fourD_dist)
       call output_phase_project(fourD_dist,ino+1,output_grids_in=.false.)

       ! The output is only done on the root process. To prevent a too big imbalance, barrier here.
       call MPI_BARRIER(MPI_COMM_WORLD,ifail)
    endif
    ! But it can be useful to have only the power exchange between steps n and n+1.
    ino = ino +1
    sim%time = target_time
    call with(sim,events, at=sim%time)
  end do



  partwriter = event(write_action(filename='part_restart.h5'))
  call with(sim, partwriter)

  call sim%finalize



contains


subroutine loop_particle_kinetic_local(sim, jorek_feedback, rng, timesteps, n_steps, particle_start_time,test_phase)
  use mod_project_particles
  use mod_random_seed
  use mod_interp, only: mode_moivre
  use mod_basisfunctions
  use mod_particle_types, only: copy_particle_kinetic_leapfrog

  implicit none

  class(particle_sim), target, intent(inout)                :: sim
  type(projection), target, intent(inout)                   :: jorek_feedback
  type(phase_space_projection), target, intent(inout)       :: test_phase
  type(count_action)                                        :: counter
  type(pcg32_rng), dimension(:), allocatable, intent(inout) :: rng
  type(particle_kinetic_leapfrog)                           :: particle_tmp

  real*8, intent(in)     :: timesteps, particle_start_time
  real*8    :: n_norm, rho_norm, t_norm, v_norm, E_norm, M_norm
  real*8    :: t, E(3), B(3), psi, U, n_e, T_e, rz_old(2), st_old(2),rzp_old(3),vcart_old(3),vcart_new(3)
  real*8    :: v_temp(3), T_eV, K_eV, v_kin_temp, B_norm(3), v,E_diff
  real*8    :: R_g, Z_g, R_s, R_t, Z_s, Z_t, xjac, HZ(n_tor), HH(4,4), HH_s(4,4), HH_t(4,4), E_tot,E_tot_red, E_after, E_after_red
  real*8    :: b_norm_r, b_norm_z, b_norm_phi, vr_tilde,vz_tilde,v_par, p_par, p_perp, p_atrop, iterations
  real*8    :: mucontainer(n_steps), vparcontainer(n_steps), bcontainer(n_steps),xcontainer(n_steps), fitsinevpar(4), fitsinemu(4), omega,mumid,vparmid
!$ real*8 :: w0, w1, mmm(3)


  integer, intent(in)   :: n_steps
  integer   :: i, j, k, l, m, i_elm_old, i_elm,i_phase
  integer   :: seed, i_rng, n_stream, ierr, nthreads
  integer   :: i_tor, index_lm, i_elm_temp
  integer   :: n_particles, ifail,proj_factor
  real*8,allocatable :: feedback_rhs(:,:,:,:,:)
  real*8, allocatable :: phase_proj(:)

!$ w0 = omp_get_wtime()

  n_norm   = CENTRAL_DENSITY * 1.d20                              ! (number) density normalisation
  rho_norm = CENTRAL_MASS * MASS_PROTON * n_norm                  ! rho_SI = rho_norm * rho
  t_norm   = sqrt((MU_ZERO * rho_norm))                           ! t_SI   = t_norm * t_jorek
  v_norm   = 1.d0 / t_norm                                        ! V_SI   = v_norm * v_jorek
  E_norm   = 1.5d0 / MU_ZERO                                      ! E_SI   = E_norm * E_jorek
  M_norm   = rho_norm * v_norm                                    ! momentum normalisation

  jorek_feedback%rhs_gather_time = jorek_feedback%rhs_gather_time + n_steps * timesteps
  allocate(feedback_rhs,source=jorek_feedback%rhs)
  allocate(phase_proj,source=test_phase%values)

  jorek_feedback%rhs = 0.d0
  feedback_rhs       = 0.d0
  proj_factor = 10
  call with(sim, counter)
  E_tot = 0.d0
  E_tot_red=0.d0
  E_after =0.d0
  E_after_red = 0.d0
  ! For energy conservation checks we sum op all the particle energies.
  select type (particles => sim%groups(1)%particles)
    type is (particle_kinetic_leapfrog)
      do j=1,size(particles,1)
        E_tot = E_tot + 0.5d0*particles(j)%weight*sim%groups(1)%mass*mass_proton*dot_product(particles(j)%v, particles(j)%v)
      enddo
  end select
  call MPI_REDUCE(E_tot,E_tot_red,1,MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  
  if(sim%my_id .eq. 0 ) write(*,"(A,E16.8)") " PARTICLE LOOP: Total energy before:", E_tot_red
  
 ! For least squares fitting
  mucontainer=0.d0
  vparcontainer = 0.d0
  bcontainer = 0.d0
  xcontainer = (/(I,I=1,n_steps)/)
  omega = 0.d0
  
  
  select type (particles => sim%groups(1)%particles)
    type is (particle_kinetic_leapfrog)
#ifdef __GFORTRAN__
      !$omp parallel do default(shared) &
#else
      !$omp parallel do default(none) &
      !$omp shared(sim, particles, n_steps, timesteps, rng, particle_start_time,        &
      !$omp rho_norm, t_norm, v_norm, E_norm, M_norm, N_norm, proj_factor,                           &
      !$omp jorek_feedback, CENTRAL_DENSITY, CENTRAL_MASS,test_phase, xcontainer)                              &
#endif
      !$omp private(particle_tmp, i_rng, i,j,k,l,m, t, E, B, psi, U, rz_old, st_old, fitsinevpar,fitsinemu,  &
      !$omp i_elm_old, i_elm, n_e, T_e, b_norm_r, b_norm_z,b_norm_phi, vr_tilde, vz_tilde,v_par, p_par, p_perp, p_atrop, omega,vparmid,mumid,&
      !$omp R_g, R_s, R_t, Z_g, Z_s, Z_t, xjac, HH, HH_s, HH_t, HZ, index_lm, ifail, v,rzp_old,vcart_old,vcart_new,E_diff, mucontainer,vparcontainer,bcontainer) &
      !$omp schedule(dynamic,10) &
      !$omp reduction(+:feedback_rhs)&
      !$omp reduction(+:phase_proj)
      do j=1,size(particles,1)
        
        call copy_particle_kinetic_leapfrog(particles(j),particle_tmp)

        !      i_rng = 1
!$      i_rng = omp_get_thread_num()+1
        E_diff = particle_tmp%weight*0.5d0*sim%groups(1)%mass*MASS_PROTON*dot_product(particle_tmp%v,particle_tmp%v)
        do k=1,n_steps

          if (particle_tmp%i_elm .le. 0) exit

          t = particle_start_time + (k-1)*timesteps

          call sim%fields%calc_EBpsiU(t, particle_tmp%i_elm, particle_tmp%st, particle_tmp%x(3), E, B, psi, U)

          rz_old    = particle_tmp%x(1:2)
          rzp_old   = particle_tmp%x
          st_old    = particle_tmp%st
          i_elm_old = particle_tmp%i_elm

          ! For least squares fitting, arrays with all the mu & vpar
          mucontainer(k)   = 0.5d0*sim%groups(1)%mass*ATOMIC_MASS_UNIT*norm2(cross_product(particle_tmp%v,B/norm2(B)))**2/norm2(B)/el_chg ! [eV/T]
          vparcontainer(k) = dot_product(particle_tmp%v,B)/norm2(B)
          bcontainer(k)    = norm2(B)
          if (particle_tmp%i_elm .gt. 0) then
            ! Do phase space projection before pushing
            ! For completeness, this is the nearest neighbour implementation.
            ! call calc_index_val_phaseproj(test_phase,particle_tmp,index_phase_tmp2,sim)
            ! if(index_phase_tmp2 > 0) then
            !  phase_proj(index_phase_tmp2)=phase_proj(index_phase_tmp2)+1.d0/n_steps*particle_tmp%weight*dot_product(particle_tmp%v,E)!
            ! endif

            ! Push the particle and determine its new location.
            
            call boris_push_cylindrical(particle_tmp, sim%groups(1)%mass, E, B, timesteps)

            call find_RZ_nearby(sim%fields%node_list, sim%fields%element_list, rz_old(1), rz_old(2), st_old(1), st_old(2), i_elm_old, &
                     particle_tmp%x(1), particle_tmp%x(2), particle_tmp%st(1), particle_tmp%st(2), particle_tmp%i_elm, ifail)
            
            ! Particle projection of energy difference each timestep. 4D
          endif 
          if(particle_tmp%i_elm .gt.0 .and. mod(k,proj_factor).eq. 0) then

            !Particle perpendicular & parallel pressure averaging
            !Normalized b vector
            b_norm_r= B(1)/sqrt(B(1)**2+B(2)**2+B(3)**2)
            b_norm_z= B(2)/sqrt(B(1)**2+B(2)**2+B(3)**2)
            b_norm_phi= B(3)/sqrt(B(1)**2+B(2)**2+B(3)**2)
            !Orthonormal (compared to magnetic field and eachother)
            vr_tilde=(-(b_norm_r)*particle_tmp%v(3)+b_norm_phi*particle_tmp%v(1))/sqrt(b_norm_phi**2+b_norm_r**2)
            vz_tilde = (particle_tmp%v(2)-b_norm_z*(b_norm_phi*particle_tmp%v(3)+b_norm_r*particle_tmp%v(1)+b_norm_z*particle_tmp%v(2)))/sqrt(b_norm_phi**2+b_norm_r**2)
            v_par= b_norm_phi*particle_tmp%v(3)+b_norm_r*particle_tmp%v(1)+b_norm_z*particle_tmp%v(2)
  
            !Parallel & perpendicular pressures
            p_perp = 1.d0/2.d0*(vr_tilde**2+vz_tilde**2)
            p_par = v_par**2
            p_atrop = p_par-p_perp
  
            call basisfunctions(particle_tmp%st(1), particle_tmp%st(2), HH, HH_s, HH_t)
            call mode_moivre(particle_tmp%x(3), HZ)
  
            i_elm=particle_tmp%i_elm
  
            do l=1,n_vertex_max
              do m=1,n_order+1
  
                index_lm = (l-1)*(n_order+1) + m
  
                v = HH(l,m) * sim%fields%element_list%element(i_elm)%size(l,m)
  
                do i_tor=1,n_tor
  
  
                  feedback_rhs(m,l,i_elm,i_tor,1) = feedback_rhs(m,l,i_elm,i_tor,1) &
  
                                                             + HZ(i_tor) * v * particle_tmp%weight * sim%groups(1)%mass * mass_proton &
  
                                                              * (p_perp+b_norm_r**2*p_atrop) * mu_zero !PI_RR
                  feedback_rhs(m,l,i_elm,i_tor,2) = feedback_rhs(m,l,i_elm,i_tor,2) &
  
                                                              + HZ(i_tor) * v * particle_tmp%weight * sim%groups(1)%mass * mass_proton &
  
                                                              * ( p_perp+b_norm_z**2*p_atrop ) * mu_zero                       !PI_ZZ
                  feedback_rhs(m,l,i_elm,i_tor,3) = feedback_rhs(m,l,i_elm,i_tor,3) &
  
                                                              + HZ(i_tor) * v * particle_tmp%weight * sim%groups(1)%mass * mass_proton &
  
                                                              * (p_perp+b_norm_phi**2*p_atrop) * mu_zero !PI_PHIPHI
                  feedback_rhs(m,l,i_elm,i_tor,4) = feedback_rhs(m,l,i_elm,i_tor,4) &
  
                                                              + HZ(i_tor) * v * particle_tmp%weight * sim%groups(1)%mass * mass_proton &
  
                                                              * ( b_norm_r*b_norm_z*p_atrop ) * mu_zero                            !PI_RZ
                  feedback_rhs(m,l,i_elm,i_tor,5) = feedback_rhs(m,l,i_elm,i_tor,5) &
  
                                                             + HZ(i_tor) * v * particle_tmp%weight * sim%groups(1)%mass * mass_proton &
  
                                                              * (b_norm_r*b_norm_phi*p_atrop ) * mu_zero !PI_RPHI
                  feedback_rhs(m,l,i_elm,i_tor,6) = feedback_rhs(m,l,i_elm,i_tor,6) &
  
                                                             + HZ(i_tor) * v * particle_tmp%weight * sim%groups(1)%mass * mass_proton &
  
                                                              *(b_norm_z*b_norm_phi*p_atrop ) * mu_zero !PI_ZPHI
                  feedback_rhs(m,l,i_elm,i_tor,7) = feedback_rhs(m,l,i_elm,i_tor,7) &
  
                                                             + HZ(i_tor) * v * particle_tmp%weight  !Density
  
  
                enddo! <tor harmonic
  
              enddo   !< order
            enddo     !< vertex
  
          end if !<particle_temp%i_elm gt 0 after pushing and projection needed 
        end do ! steps
        
        E_diff =-E_diff+ particle_tmp%weight*0.5d0*sim%groups(1)%mass*mass_proton*dot_product(particle_tmp%v,particle_tmp%v)

        ! For the fitting, it is important to take a gyrofrequency that makes sense.
        ! Although the average B is not expected to give the exact gyrofrequency, it is often good enough.

        ! If you want to check that this is working, simply add if(sim%my_id .eq. 0 .and. j.eq.1 )statements with 
        ! outputs into files of the vparcontainer & mucontainer & fits.
        omega = el_chg*sum(bcontainer)/size(bcontainer)/sim%groups(1)%mass/atomic_mass_unit 
        fitsinevpar = l2_sinfit(n_steps,xcontainer*timesteps,vparcontainer,omega)
        fitsinemu   = l2_sinfit(n_steps,xcontainer*timesteps,mucontainer,  omega)

        ! Projecting at the linearized part at the middle. As long as these quantities
        ! do not vary too much over one particle loop, this is good enough.
        ! You could improve it by also using an Econtainer and depositing at every (or every proj_factor) timestep
        ! of the linear part of the sine fits. It is not recommended to fit this Econtainer, as this might break conservation
        ! of energy. 
        vparmid =  fitsinevpar(1)+fitsinevpar(2)*real(n_steps/2,8)*timesteps
        mumid   = fitsinemu(1)+fitsinemu(2)*real(n_steps/2,8)*timesteps
        call project_single_particle_x(test_phase,[vparmid,mumid],phase_proj,E_diff)

        call copy_particle_kinetic_leapfrog(particle_tmp, particles(j))
      end do   ! particles
      !$omp end parallel do

  end select
  select type(particles => sim%groups(1)%particles)
    type is (particle_kinetic_leapfrog)
      do j=1,size(particles,1)
        E_after = E_after + 0.5d0*particles(j)%weight*sim%groups(1)%mass*mass_proton*dot_product(particles(j)%v, particles(j)%v)
      enddo
  end select
  call MPI_REDUCE(E_after,E_after_red,1,MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  
  if(sim%my_id .eq. 0 ) then 
    write(*,"(A,E16.8)") " PARTICLE LOOP: Total energy after:", E_after_red
    write(*,"(A,E16.8)") " PARTICLE LOOP: Energy difference: ", E_after_red - E_tot_red
  endif
  iterations = real(n_steps/proj_factor,8)

  if(sim%my_id .eq.0 )write(*,"(A,I4,A,F4.0,A)") " PARTICLE LOOP: Projecting every",proj_factor," timesteps, for a total of ", iterations, " times."
  if (sim%my_id .eq. 0) write(*,*) "PARTICLE LOOP: End of the particle loop"
  jorek_feedback%rhs = feedback_rhs/iterations
  test_phase%values=test_phase%values+phase_proj
  deallocate(feedback_rhs)
  deallocate(phase_proj)

end subroutine

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

pure function proj_RZPE(ndim, sim,group, particle)
  use mod_import_experimental_dist, only: calculate_B
  type(particle_sim),           intent(in) :: sim
  integer,                      intent(in) :: ndim,group
  class(particle_base),         intent(in) :: particle
  real*8                                   :: B_tmp(3)
  real*8                                   :: proj_RZPE(ndim)
  B_tmp= calculate_B(sim%fields, particle%i_elm,particle%st(1),particle%st(2),particle%x(3))
  select type( p => particle)
    type is (particle_kinetic_leapfrog)
      proj_RZPE = [particle%x(1),particle%x(2),dot_product(p%v,B_tmp)/norm2(B_tmp)/norm2(p%v),0.5d0*sim%groups(1)%mass*mass_proton*dot_product(p%v, p%v)/1000.d0/el_CHG]
  end select  
end function proj_RZPE




function l2_sinfit(n,x,y,gyro_frequency) result(xfit)
  integer, intent(in) :: n
  real*8, intent(in)  :: x(n),y(n),gyro_frequency !gyro_frequency in rad/timestep
  real*8              :: xfit(4)
  real*8              :: a11,a12,a13,a14,a22,a23,a24,a33,a34,a44, r3(n),r4(n), determinant
  real*8              :: cofactor(4,4),XTY(4)
  integer             :: it
  r3 = cos(gyro_frequency*x)
  r4 = sin(gyro_frequency*x)
  a11=n
  a12=sum(x)
  a13=sum(r3)
  a14=sum(r4)
  a22=dot_product(x,x)
  a23=dot_product(x,r3)
  a24=dot_product(x,r4)
  a33=dot_product(r3,r3)
  a34=dot_product(r3,r4)
  a44=dot_product(r4,r4)
  determinant = a14**2*a23**2-2*a13*a14*a23*a24+a13**2*a24**2-a14**2*a22*a33+2*a12*a14*a24*a33-a11*a24**2*a33+2*a13*a14*a22*a34-2*a12*a14*a23*a34-2*a12*a13*a24*a34+2*a11*a23*a24*a34+a12**2*a34**2-a11*a22*a34**2-a13**2*a22*a44+2*a12*a13*a23*a44-a11*a23**2*a44-a12**2*a33*a44+a11*a22*a33*a44
  XTY(1) = sum(y)
  XTY(2) = dot_product(x,y)
  XTY(3) = dot_product(r3,y)
  XTY(4) = dot_product(r4,y)
  cofactor(1,1) = -a24**2*a33 + 2*a23*a24*a34 - a22*a34**2 - a23**2*a44 + a22*a33*a44
  cofactor(1,2) = a14*a24*a33 - a14*a23*a34 - a13*a24*a34 + a12*a34**2 + a13*a23*a44 - a12*a33*a44
  cofactor(1,3) = -a14*a23*a24 + a13*a24**2 + a14*a22*a34 - a12*a24*a34 - a13*a22*a44 + a12*a23*a44
  cofactor(1,4) = a14*a23**2 - a13*a23*a24 - a14*a22*a33 + a12*a24*a33 + a13*a22*a34 - a12*a23*a34
  cofactor(2,1) = cofactor(1,2)
  cofactor(2,2) = -a14**2*a33 + 2*a13*a14*a34 - a11*a34**2 - a13**2*a44 + a11*a33*a44
  cofactor(2,3) = a14**2*a23 - a13*a14*a24 - a12*a14*a34 + a11*a24*a34 + a12*a13*a44 - a11*a23*a44
  cofactor(2,4) = -a13*a14*a23 + a13**2*a24 + a12*a14*a33 - a11*a24*a33 - a12*a13*a34 + a11*a23*a34
  cofactor(3,1) = cofactor(1,3)
  cofactor(3,2) = cofactor(2,3)
  cofactor(3,3) = -a14**2*a22 + 2*a12*a14*a24 - a11*a24**2 - a12**2*a44 + a11*a22*a44
  cofactor(3,4) = a13*a14*a22 - a12*a14*a23 - a12*a13*a24 + a11*a23*a24 + a12**2*a34 - a11*a22*a34
  cofactor(4,1) = cofactor(1,4)
  cofactor(4,2) = cofactor(2,4)
  cofactor(4,3) = cofactor(3,4)
  cofactor(4,4) = -a13**2*a22 + 2*a12*a13*a23 - a11*a23**2 - a12**2*a33 + a11*a22*a33
  
  cofactor = cofactor/determinant
  !do it=1,4
  !  write(*,"(4E16.8)") cofactor(it,:)
  !enddo 
  xfit = matmul(cofactor, XTY)
  !write(*,*) "OMEGA ", gyro_frequency
  
  !write(*,"(A,4E16.8)") "XFIT:", xfit
  !write(*,"(A,4E16.8)") "XTY:", XTY
  
end function l2_sinfit
end program tae_phase_space_project

