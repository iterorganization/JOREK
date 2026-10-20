!> this module contains the functionality related to performing physics on particles that needs to happen each particle step
!> the main functionality is all contained in evolve_particle_group, so below for more explanation
module mod_particle_evolution
    use mod_model_settings
    use particle_tracer
    use phys_module, only: CENTRAL_MASS, CENTRAL_DENSITY
    use phys_module, only: use_manual_random_seed, n_aux_var, part_kill_ratio, proj_collection_period
    use mod_coupling_settings
    use coupling_variables
    use mod_project_particles
    use mod_random_seed
    use mod_interp, only: mode_moivre, interp_0
    use mod_basisfunctions
    use mod_particle_types, only: copy_particle_kinetic_leapfrog
    use mod_sampling, only: boxmueller_transform,sample_chi_squared_3
    use mod_coordinate_transforms, only: vector_cartesian_to_cylindrical
    !$ use omp_lib

    implicit none
    private
    public :: evolve_particle_group, evolve_REs
contains

  !> For each particle group, this function does the following:
  !> - performs coupling scheme specific physics (e.g. ionisation, radiation... etc)
  !> - creates the feedback rhs required for projections of kinetic variables
  !> - evolves the location and velocity of the particles based on the background plasma (pushing)
  subroutine evolve_particle_group(sim, group_num, jorek_feedback, rng, tstep_part_adj, nstep_part_adj)
    use mod_project_particles
    use mod_random_seed
    use mod_interp, only: mode_moivre
    use mod_basisfunctions
    use mod_particle_types, only: copy_particle_kinetic_leapfrog
    use mod_sampling, only: boxmueller_transform,sample_chi_squared_3

    implicit none
    class(particle_sim), target, intent(inout)                :: sim
    integer, intent(in)                                       :: group_num
    type(projection), target, intent(inout)                   :: jorek_feedback
    type(count_action)                                        :: counter
    type(pcg32_rng), dimension(:), allocatable, intent(inout) :: rng
    real*8,  intent(in)                                       :: tstep_part_adj
    integer, intent(in)                                       :: nstep_part_adj
    
    real*8,allocatable :: feedback_rhs(:,:,:,:,:)
    type (type_node_list),         pointer :: feedback_nodelist
    type (type_element_list),      pointer :: feedback_element_list
    type (particle_group),         pointer :: part_group
    character(len=3) :: cs

    !> Coupling scheme specific
    integer :: imp_q_idx

    !> ================================ INITIALISATION =======================================
    part_group => sim%groups(group_num)
    if (sim%my_id .eq. 0) write(*,*) '---------- Evolving particle group: ', part_group%id, " ----------"

    !> if ics, determine index for impurity charge projection
    if (part_group%coupling_scheme == 'ics') then
      imp_q_idx = ics_indices_kin(part_group%ics_group_idx)
    endif

    !> Set up storage of feedback
    allocate(feedback_rhs,source=jorek_feedback%rhs)
    feedback_nodelist => jorek_feedback%node_list
    feedback_element_list => jorek_feedback%element_list
    feedback_rhs       = 0.d0
    
    !> count number of particles in system, and update sim%groups(...)%average_weight
    call with(sim, counter)
    
    !> ================================ COUPLING SPECIFIC LOOPS =======================================
    !> gathers feedback rhs per particle per tstep_part_adj and pushes particle
    !> this is where coupling specific physics such as ionisation, charge exchange... etc happens

    select case (part_group%coupling_scheme)
      case ('ncs')
        call evolve_ncs_ics(sim, group_num, feedback_rhs, feedback_nodelist, feedback_element_list, rng, tstep_part_adj, nstep_part_adj)
      case ('ics')
        call evolve_ncs_ics(sim, group_num, feedback_rhs, feedback_nodelist, feedback_element_list, rng, tstep_part_adj, nstep_part_adj, imp_q_idx)
      case ('rep')
        call evolve_REs(sim, group_num, feedback_rhs, rng, tstep_part_adj, nstep_part_adj)
      case ('epf')
        call evolve_epf(sim, group_num, feedback_rhs, rng, tstep_part_adj, nstep_part_adj)
      case default
        write(*,*) "ERROR: Unknown coupling scheme: '", part_group%coupling_scheme, "' found for group '", part_group%id, "'"
        stop 1
    end select
    
    ! ================================= CONSTRUCT PROJECTION RHS =======================================
    !> enter gathered rhs into jorek_feedback
    if (part_group%coupling_scheme == 'ncs' .or. part_group%coupling_scheme == 'ics') then
      ! To get rates in the feedback, we need to divide the change by the time. Since we keep adding changes each evolve_particle_group call until 
      ! the rhs is reset to 0 when the jorek_feedback is projected (each fluid tstep), we should divide by sim%tstep_fluid_si
      jorek_feedback%rhs(:,:,:,:,mom_par_idx_kin) = jorek_feedback%rhs(:,:,:,:,mom_par_idx_kin) + feedback_rhs(:,:,:,:,mom_par_idx_kin) / sim%tstep_fluid_si
#ifdef WITH_TiTe
      jorek_feedback%rhs(:,:,:,:,E_Te_idx_kin)    = jorek_feedback%rhs(:,:,:,:,E_Te_idx_kin)    + feedback_rhs(:,:,:,:,E_Te_idx_kin)    / sim%tstep_fluid_si
      jorek_feedback%rhs(:,:,:,:,E_Ti_idx_kin)    = jorek_feedback%rhs(:,:,:,:,E_Ti_idx_kin)    + feedback_rhs(:,:,:,:,E_Ti_idx_kin)    / sim%tstep_fluid_si
#else
      jorek_feedback%rhs(:,:,:,:,E_idx_kin)       = jorek_feedback%rhs(:,:,:,:,E_idx_kin)       + feedback_rhs(:,:,:,:,E_idx_kin)       / sim%tstep_fluid_si
#endif
      !> ncs specific projections
      if (part_group%coupling_scheme == 'ncs') then
        jorek_feedback%rhs(:,:,:,:,rho_idx_kin)   = jorek_feedback%rhs(:,:,:,:,rho_idx_kin)     + feedback_rhs(:,:,:,:,rho_idx_kin)     / sim%tstep_fluid_si
        ! for the density, we should divide by the amount of times we will double count the same particle (=the number of particle steps in a fluid step)
        jorek_feedback%rhs(:,:,:,:,6)             = jorek_feedback%rhs(:,:,:,:,6)               + feedback_rhs(:,:,:,:,6)               / (sim%tstep_fluid_si/tstep_part_adj) !< extra diagnostic projection (density) 
      endif

      !> ics specific projections
      if (part_group%coupling_scheme == 'ics') then
        jorek_feedback%rhs(:,:,:,:,imp_q_idx)     = jorek_feedback%rhs(:,:,:,:,imp_q_idx)       + feedback_rhs(:,:,:,:,imp_q_idx)       / (sim%tstep_fluid_si/tstep_part_adj)
        jorek_feedback%rhs(:,:,:,:,7)             = jorek_feedback%rhs(:,:,:,:,7)               + feedback_rhs(:,:,:,:,7)               / sim%tstep_fluid_si                  !< extra projection (impurity radiated power)
        jorek_feedback%rhs(:,:,:,:,8)             = jorek_feedback%rhs(:,:,:,:,8)               + feedback_rhs(:,:,:,:,8)               / (sim%tstep_fluid_si/tstep_part_adj) !< extra projection (impurity density)
      endif
    endif

    !> rep specific projections
    if (part_group%coupling_scheme == 'rep') then
      feedback_rhs = feedback_rhs / real(nstep_part_adj,8) 
      jorek_feedback%rhs(:,:,:,:,P_par_idx_kin)   = jorek_feedback%rhs(:,:,:,:,P_par_idx_kin)   + feedback_rhs(:,:,:,:,P_par_idx_kin)
      jorek_feedback%rhs(:,:,:,:,P_perp_idx_kin)  = jorek_feedback%rhs(:,:,:,:,P_perp_idx_kin)  + feedback_rhs(:,:,:,:,P_perp_idx_kin)
      jorek_feedback%rhs(:,:,:,:,j_Phi_idx_kin)   = jorek_feedback%rhs(:,:,:,:,j_Phi_idx_kin)   + feedback_rhs(:,:,:,:,j_Phi_idx_kin)
    endif

    !> epf specific projection
    if (part_group%coupling_scheme == 'epf') then
      jorek_feedback%rhs(:,:,:,:,PI_RR_idx_kin)     = jorek_feedback%rhs(:,:,:,:,PI_RR_idx_kin)     + feedback_rhs(:,:,:,:,PI_RR_idx_kin)
      jorek_feedback%rhs(:,:,:,:,PI_ZZ_idx_kin)     = jorek_feedback%rhs(:,:,:,:,PI_ZZ_idx_kin)     + feedback_rhs(:,:,:,:,PI_ZZ_idx_kin)
      jorek_feedback%rhs(:,:,:,:,PI_PHIPHI_idx_kin) = jorek_feedback%rhs(:,:,:,:,PI_PHIPHI_idx_kin) + feedback_rhs(:,:,:,:,PI_PHIPHI_idx_kin)
      jorek_feedback%rhs(:,:,:,:,PI_RZ_idx_kin)     = jorek_feedback%rhs(:,:,:,:,PI_RZ_idx_kin)     + feedback_rhs(:,:,:,:,PI_RZ_idx_kin)
      jorek_feedback%rhs(:,:,:,:,PI_RPHI_idx_kin)   = jorek_feedback%rhs(:,:,:,:,PI_RPHI_idx_kin)   + feedback_rhs(:,:,:,:,PI_RPHI_idx_kin)
      jorek_feedback%rhs(:,:,:,:,PI_ZPHI_idx_kin)   = jorek_feedback%rhs(:,:,:,:,PI_ZPHI_idx_kin)   + feedback_rhs(:,:,:,:,PI_ZPHI_idx_kin)
      jorek_feedback%rhs(:,:,:,:,rho_ep_idx_kin)    = jorek_feedback%rhs(:,:,:,:,rho_ep_idx_kin)    + feedback_rhs(:,:,:,:,rho_ep_idx_kin)
    endif

    jorek_feedback%rhs_gather_time = 0.d0
    deallocate(feedback_rhs)
    
    if (sim%my_id .eq. 0) write(*,*) '---------- Finished evolving group: ', part_group%id, " ----------"
    
  end subroutine evolve_particle_group

  subroutine evolve_REs(sim, group_num, feedback_rhs, rng, tstep_part_adj, nstep_part_adj)
    use mod_project_particles
    use mod_random_seed
    use mod_interp, only: mode_moivre
    use mod_basisfunctions
    use mod_particle_types, only: copy_particle_kinetic_leapfrog
    use mod_sampling, only: boxmueller_transform,sample_chi_squared_3
    
    implicit none
    class(particle_sim), target, intent(inout)                :: sim
    integer, intent(in)                                       :: group_num
    real*8,allocatable,          intent(inout)                :: feedback_rhs(:,:,:,:,:)
    type(count_action)                                        :: counter
    type(pcg32_rng), dimension(:), allocatable, intent(inout) :: rng
    real*8,  intent(in)                                       :: tstep_part_adj
    integer, intent(in)                                       :: nstep_part_adj

    character(len=3) :: cs

    !> RE specific variables
    real*8    :: n_norm, rho_norm
    real*8    :: HZ(n_tor), HH(4,4), HH_s(4,4), HH_t(4,4)
    real*8    :: v_jR, v_jZ, v_jPhi, v_Ppar, v_Pperp

    real*8    :: E(3), B(3), B_norm2(3), psi, U
    real*8    :: cylindrical_velocity(3), cylindrical_momentum(3)
    real*8    :: v_par, v_perp, gamma_m, proj_factor
    integer   :: j, k, m, n, ifail, i_tor, n_lost

    n_norm   = CENTRAL_DENSITY * 1.d20                              ! (number) density normalisation
    rho_norm = CENTRAL_MASS * ATOMIC_MASS_UNIT * n_norm                  ! rho_SI = rho_norm * rho

    ! Loop over all particle groups
    n_lost = 0
    select type (particles => sim%groups(group_num)%particles)
    type is (particle_kinetic_relativistic)
      if(use_manual_random_seed) then
        !$ call omp_set_schedule(omp_sched_static,10)
      else
        !$ call omp_set_schedule(omp_sched_dynamic,10)
      end if  
      !$omp parallel do default(none) &
      !$omp schedule(runtime)         &
      !$omp private(j, k, m, n, HZ, HH, HH_s, HH_t, E, B, psi, U,  &
      !$omp B_norm2, proj_factor, v_Ppar, v_Pperp, v_jPhi, i_tor, ifail, &
      !$omp cylindrical_velocity, cylindrical_momentum, v_par, v_perp, gamma_m ) &
      !$omp shared (nstep_part_adj, tstep_part_adj, sim, group_num, rho_norm, &
      !$omp P_par_idx_kin, P_perp_idx_kin, j_phi_idx_kin) &
      !$omp reduction(+:feedback_rhs)
  
      do j=1,size(particles,1)
        do k=1,nstep_part_adj
          if (particles(j)%i_elm .le. 0) exit
  
          call basisfunctions(particles(j)%st(1), particles(j)%st(2), HH, HH_s, HH_t)
          call mode_moivre(particles(j)%x(3), HZ)
  
          ! Determines velocity in cylindrical coordinates
          cylindrical_momentum = vector_cartesian_to_cylindrical(particles(j)%x(3), particles(j)%p)
          cylindrical_velocity = cylindrical_momentum / &
                                  sqrt(dot_product(cylindrical_momentum,cylindrical_momentum)/SPEED_OF_LIGHT**2 + sim%groups(group_num)%mass**2)
  
          ! Uncomment for PCS
          call sim%fields%calc_EBpsiU(sim%time, particles(j)%i_elm, particles(j)%st, particles(j)%x(3), E, B, psi, U)
          B_norm2 = B/norm2(B)
  
          v_par   = dot_product(cylindrical_velocity, B_norm2)
          v_perp  = norm2(cylindrical_velocity - v_par * B_norm2)
          gamma_m = sqrt(MASS_ELECTRON**2 + dot_product(cylindrical_momentum,cylindrical_momentum)*ATOMIC_MASS_UNIT**2/SPEED_OF_LIGHT**2)
  
          do n=1,n_degrees
            do m=1,n_vertex_max
  
              proj_factor = HH(m,n) * sim%fields%element_list%element(particles(j)%i_elm)%size(m,n) * particles(j)%weight
  

              ! PCS for REs
              v_Ppar  = proj_factor * gamma_m * v_par**2 * MU_ZERO
              v_Pperp = proj_factor * gamma_m * v_perp**2 / 2.d0 * MU_ZERO
  
              v_jPhi  = - proj_factor * real(particles(j)%q, 8) * EL_CHG * cylindrical_velocity(3) * particles(j)%x(1) * MU_ZERO
              

              do i_tor = 1,n_tor
                feedback_rhs(n,m,particles(j)%i_elm,i_tor,P_par_idx_kin)  = feedback_rhs(n,m,particles(j)%i_elm,i_tor,P_par_idx_kin)  + HZ(i_tor)*v_Ppar
                feedback_rhs(n,m,particles(j)%i_elm,i_tor,P_perp_idx_kin) = feedback_rhs(n,m,particles(j)%i_elm,i_tor,P_perp_idx_kin) + HZ(i_tor)*v_Pperp
                feedback_rhs(n,m,particles(j)%i_elm,i_tor,j_Phi_idx_kin)  = feedback_rhs(n,m,particles(j)%i_elm,i_tor,j_Phi_idx_kin)  + HZ(i_tor)*v_jPhi

                !> inverse implementation from Hannes (TODO: reverse data structure of the projections generally)
                ! feedback_rhs_inv(P_par_idx_kin,i_tor,particles(j)%i_elm,m,n) = feedback_rhs_inv(P_par_idx_kin,i_tor,particles(j)%i_elm,m,n) + HZ(i_tor)*v_Ppar
                ! feedback_rhs_inv(P_perp_idx_kin,i_tor,particles(j)%i_elm,m,n) = feedback_rhs_inv(P_perp_idx_kin,i_tor,particles(j)%i_elm,m,n) + HZ(i_tor)*v_Pperp
                ! feedback_rhs_inv(j_Phi_idx_kin,i_tor,particles(j)%i_elm,m,n) = feedback_rhs_inv(j_Phi_idx_kin,i_tor,particles(j)%i_elm,m,n) + HZ(i_tor)*v_jPhi

              enddo
            enddo
          enddo
  
          call volume_preserving_push_jorek(particles(j),sim%fields,sim%groups(group_num)%mass,sim%time,tstep_part_adj,ifail)
  
  
        end do !< steps
      end do !< particles
      !$omp end parallel do 
  
    end select
    
  end subroutine evolve_REs



  !> Internal function for gathering the feedback rhs values when using the ncs or ics coupling scheme
  !> The two coupling schemes are handled by the same function due to large degree of overlap in the physics 
  !> experienced by neutrals and impurities. 
  !> The pushing of the particle is also done here
  subroutine evolve_ncs_ics(sim, group_num, feedback_rhs, feedback_nodelist, feedback_element_list, rng, tstep_part_adj, nstep_part_adj, imp_q_idx)
    use mod_collisions
    use mod_ionisation_recombination

    implicit none
    class(particle_sim),       target, intent(inout)          :: sim
    integer,                           intent(in)             :: group_num
    real*8,allocatable,                intent(inout)          :: feedback_rhs(:,:,:,:,:)
    type (type_node_list),    pointer, intent(in)             :: feedback_nodelist
    type (type_element_list), pointer, intent(in)             :: feedback_element_list
    type(pcg32_rng), dimension(:), allocatable, intent(inout) :: rng
    real*8,                            intent(in)             :: tstep_part_adj
    integer,                           intent(in)             :: nstep_part_adj
    integer, optional,                 intent(in)             :: imp_q_idx


    real*8, parameter  :: H_binding_energy = 2.18d-18

    !> Diagnostics ----------------------------------- 
    real*8    :: n_lost_ion, n_lost_ion_all, p_lost_ion, p_lost_ion_all 
    real*8    :: p_lost_cx, p_lost_cx_all, p_lost_plt, p_lost_plt_all 
    integer   :: n_super_ionized, n_super_ionized_all

    !> Coupling --------------------------------------
    real*8    :: density_source, mom_par_source, energy_source
    real*8    :: energy_source_Te, energy_source_Ti
    real*8    :: density_fb, mom_par_fb, E_fb, imp_q_fb, imp_density_fb, imp_P_rad_fb, extra_proj, imp_P_line_rad_fb
    real*8    :: E_fb_Te, E_fb_Ti
    real*8    :: v_old(3), v_new(3), T_eV, B_norm(3)
    real*8    :: vvector(3), ran_norm(4)

    logical   :: limits, limits_coll
    real*8    :: T_e_raw, T_i_raw,  n_e_raw
    real*8    :: ionize_rate, ionize_energy, ionize_source, ionize_prob
    real*8    :: cx_rate, cx_energy,  cx_source, cx_prob
    real*8    :: PLT, PRB, Srec, line_rad_energy
    real*8    :: ionize_ran(1), ionize_ran_imp(2), cx_ran(8), st_ran(2) 
    real*8    :: kinetic_energy, radiation_energy, binding_energy
    real*8    :: imp_charge_density ! impurity charge density in units of [e]

    !> impurity collision
    integer(kind=1) :: q_b
    integer, parameter :: n_coll=20
    real*8    :: ran(6), ran2(6,n_coll), q(3), m_b
    real*8    :: coulomb_log, kTb, n_b, v_b(3,n_coll) 
    real*8, dimension(1) :: P, P_s, P_t, P_phi, P_time
    real*8    :: delta_E_kin

    !> System variables ------------------------------
    type(particle_kinetic_leapfrog) :: particle_tmp

    real*8    :: n_norm, rho_norm, t_norm, v_norm, E_norm, M_norm
    real*8    :: t, E(3), B(3), psi, U, n_i, n_e, T_e, T_i, grad_T_i(3), rz_old(2), st_old(2)
    real*8    :: R_g, Z_g, R_s, R_t, Z_s, Z_t, xjac, R, Z
    real*8    :: HZ(n_tor), HH(4,4), HH_s(4,4), HH_t(4,4)

    integer   :: i, j, k, l, m, n, i_elm_old, i_elm, q_old 
    integer   :: seed, i_rng, n_stream, ierr, nthreads
    integer   :: i_tor, i_elm_temp
    integer   :: n_particles, ifail 
    integer   :: imp_q_idx_temp

    !> Initialise diagnostics to 0
    n_lost_ion     = 0.d0; n_lost_ion_all = 0.d0; p_lost_ion     = 0.d0; p_lost_ion_all  = 0.d0; p_lost_plt          = 0.d0; 
    p_lost_plt_all = 0.d0; p_lost_cx      = 0.d0; p_lost_cx_all  = 0.d0; n_super_ionized = 0;    n_super_ionized_all = 0;

    !> Normalise variables 
    n_norm   = CENTRAL_DENSITY * 1.d20                              ! (number) density normalisation
    rho_norm = CENTRAL_MASS * ATOMIC_MASS_UNIT * n_norm             ! rho_SI = rho_norm * rho
    t_norm   = sqrt((MU_ZERO * rho_norm))                           ! t_SI   = t_norm * t_jorek
    v_norm   = 1.d0 / t_norm                                        ! V_SI   = v_norm * v_jorek
    E_norm   = 1.5d0 / MU_ZERO                                      ! E_SI   = E_norm * E_jorek
    M_norm   = rho_norm * v_norm                                    ! momentum normalisation

    if(use_manual_random_seed) then
      !$ call omp_set_schedule(omp_sched_static,10)
    else
      !$ call omp_set_schedule(omp_sched_dynamic,10)
    end if  
#ifdef __GFORTRAN__
    !$omp parallel do default(shared) & ! workaround for Error: «__vtab_mod_pcg32_rng_Pcg32_rng» not specified in enclosing «parallel»
#else
    !$omp parallel do default(none)                                                                       &
#endif
    !$omp schedule(runtime)                                                                               &
    !$omp shared(sim, group_num, nstep_part_adj, tstep_part_adj, rng,                                    &
    !$omp rho_norm, t_norm, v_norm, E_norm, M_norm, N_norm, part_kill_ratio,                              &    
    !$omp rho_idx_kin, mom_par_idx_kin,                                                                   &
#ifdef WITH_TiTe
    !$omp E_Te_idx_kin, E_Ti_idx_kin,                                                                     &
#else
    !$omp E_idx_kin,                                                                                      &
#endif
    !$omp imp_q_idx, ics_indices_kin,                                                                     &
    !$omp CENTRAL_DENSITY, CENTRAL_MASS, feedback_nodelist, feedback_element_list)                        &
    !$omp private(particle_tmp, i_rng, i, j, k, l, m, t, E, B, psi, U, rz_old, st_old,                    &
    !$omp i_elm_old, i_elm, n_i, n_e, T_e, T_i, imp_charge_density, PLT, PRB, Srec, q_old,                &
    !$omp ionize_rate, ionize_prob, ionize_ran, ionize_ran_imp, ionize_source, ionize_energy,             &
    !$omp cx_rate, cx_prob, cx_source, cx_energy, cx_ran, grad_T_i,                                       &
    !$omp kinetic_energy, line_rad_energy, radiation_energy, binding_energy,                              &  
    !$omp R_s, R_t, Z_g, Z_s, Z_t, R, Z, xjac, HH, HH_s, HH_t, HZ, ifail,                                 &
    !$omp density_fb, E_fb, mom_par_fb,extra_proj, imp_q_fb, imp_density_fb, imp_P_rad_fb,                &
    !$omp density_source, mom_par_source, energy_source, v_old, v_new, T_eV, imp_P_line_rad_fb,           &
    !$omp m_b, kTb,coulomb_log ,n_b,v_b, ran, ran2, q_b, q, E_fb_Te, E_fb_Ti,                             &
    !$omp P, P_s, P_t, P_phi, P_time, limits, limits_coll, energy_source_Te, energy_source_Ti,            &
    !$omp vvector, ran_norm, imp_q_idx_temp, T_e_raw, T_i_raw, n_e_raw, delta_E_kin, i_tor, n)            &
    !$omp reduction(+:feedback_rhs,n_lost_ion,p_lost_plt,p_lost_cx,p_lost_ion,n_super_ionized)
    do j=1,size(sim%groups(group_num)%particles,1)
      particle_tmp = sim%groups(group_num)%particles(j)
      !$ i_rng = omp_get_thread_num()+1
      do k=1,nstep_part_adj

        !> exit evolution loop if particle is outside domain
        if (particle_tmp%i_elm .le. 0) exit
        
        t = sim%time + (k-1)*tstep_part_adj

        ionize_source = 0.d0
        ionize_energy = 0.d0
        line_rad_energy = 0.d0
        radiation_energy = 0.d0
        cx_source = 0.d0
        cx_energy = 0.d0
        delta_E_kin = 0.d0

        !> calculate local fields
        call sim%fields%calc_EBpsiU(t, particle_tmp%i_elm, particle_tmp%st, particle_tmp%x(3), E, B, psi, U)
        rz_old    = particle_tmp%x(1:2)
        st_old    = particle_tmp%st
        i_elm_old = particle_tmp%i_elm
        q_old     = particle_tmp%q 
        v_old     = particle_tmp%v

        v_new = v_old
        
        !> calculate ion density and electron temperature (jorek model assumption: n_e = n_i)     
#ifdef WITH_TiTe
          call sim%fields%calc_NeTeTi(t, particle_tmp%i_elm, particle_tmp%st, particle_tmp%x(3),n_e=n_i, T_e=T_e, T_i=T_i, n_e_raw=n_e_raw, &
                            T_e_raw=T_e_raw, T_i_raw=T_i_raw, grad_T_i=grad_T_i)
          limits = (n_e_raw .le. 1e14) .or. (T_e_raw * K_BOLTZ / EL_CHG .le. 1.d0) .or. (T_i_raw * K_BOLTZ / EL_CHG .le. 1.d0) !ADAS limits
          limits_coll = T_i_raw * K_BOLTZ / EL_CHG < 1.d0 !< limits for collisions
#else
          call sim%fields%calc_NeTeTi(t, particle_tmp%i_elm, particle_tmp%st, particle_tmp%x(3), n_e=n_i, T_e=T_e, n_e_raw=n_e_raw, T_e_raw=T_e_raw, grad_T_e=grad_T_i)
          limits = (n_e_raw .le. 1e14) .or. (T_e_raw * K_BOLTZ / EL_CHG .le. 1.d0)
          limits_coll = T_e_raw * K_BOLTZ / EL_CHG < 1.d0 !< limits for collisions
#endif

        !> loop over impurities groups and calculate their contribution to electron density
        imp_charge_density = 0.d0
        do n=1, size(sim%groups)
          if (trim(sim%groups(n)%coupling_scheme) == 'ics') then
            imp_q_idx_temp = ics_indices_kin(sim%groups(n)%ics_group_idx)
            call interp_0(feedback_nodelist, feedback_element_list, particle_tmp%i_elm, [imp_q_idx_temp], 1 , particle_tmp%st(1), particle_tmp%st(2), particle_tmp%x(3), P)
            imp_charge_density = imp_charge_density + P(1) ! charge density of impurities in units of [e]
          endif
        enddo
        
        !> adjust n_e based on impurity charge
        n_e = n_i + max(0.d0, imp_charge_density)
        
        !> check that particle weight is non negative
        if (particle_tmp%weight .lt. 0.0d0) write(*,*) "Negative particle weight p(j)%w=", particle_tmp%weight
        
        ! ============================================ NCS SPECIFIC PHYSICS ===========================================
        if (sim%groups(group_num)%coupling_scheme == 'ncs') then
          !> calculate fluid flow velocity [v_R, v_Z, v_phi] m/s
          call sim%fields%calc_vvector(t, particle_tmp%i_elm, particle_tmp%st, particle_tmp%x(3), vvector)

          !> RADIATION (only line radiation* for ncs)
          !>   (*Line radiation due to collisional excitation of the neutral's bound electrons with the 
          !>    electrons in the background plasma. The spectrum of these radiations are discrete.)
          if (sim%groups(group_num)%use_kin_radiation .and. .not. limits) then
            call sim%groups(group_num)%ad%PLT%interp(int(particle_tmp%q), log10(n_e), log10(T_e), PLT) ! [J m^3/s]
            line_rad_energy = n_e * particle_tmp%weight * PLT * tstep_part_adj
          endif ! RADIATION
          
          !> IONISATION (Neutrals)
          if (sim%groups(group_num)%use_kin_ionisation .and. .not. limits) then
            call sim%groups(group_num)%ad%SCD%interp(int(particle_tmp%q), log10(n_e), log10(T_e), ionize_rate) ! [m^3/s]
            ionize_prob = 1.d0 - exp(-ionize_rate * n_e * tstep_part_adj) ! [0] poisson point process, exponential 

            ! If the weight is to small throw away the particle with the probability, else reduce weight with ionising probability
            ionize_source = 0.d0
    
            if (particle_tmp%weight .le. sim%groups(group_num)%average_weight * part_kill_ratio) then
              call rng(i_rng)%next(ionize_ran)
              if (ionize_ran(1) .le. ionize_prob) then
                particle_tmp%i_elm  = 0
                ionize_source = particle_tmp%weight
                particle_tmp%weight = 0.d0
                !superparticles ionized
                n_super_ionized = n_super_ionized +1
              else
                ionize_source = 0.d0
              endif
            else 
              ionize_source = particle_tmp%weight * ionize_prob
              particle_tmp%weight = particle_tmp%weight * (1.d0 - ionize_prob)
            endif 
    
            kinetic_energy = dot_product(v_old,v_old) * sim%groups(group_num)%mass * ATOMIC_MASS_UNIT / 2.d0
            ionize_energy     = kinetic_energy - H_binding_energy
            !<including binding energy will make ionize_energy negative, so it becomes a sink for the plasma
          endif ! IONISATION
          
          !> CHARGE EXCHANGE
          ! It is assumed that we will have a exchange between hydrogen isotopes
          if (sim%groups(group_num)%use_kin_cx  .and. .not. limits .and. particle_tmp%weight .gt. 0.d0) then !< CX uses adas as well. Te limit could be lower.
            call sim%groups(group_num)%ad%CCD%interp(int(particle_tmp%q+1), log10(n_e), log10(T_e), cx_rate) ! [m^3/s]
            CX_prob = 1.d0 - exp(-cx_rate * n_e * tstep_part_adj)
    
            call rng(i_rng)%next(cx_ran)
            if (cx_ran(1) .le. CX_prob) then
              ! sample boltzman, randomize velocity
  
              !> ----- NEW CX PARTICLE ---------
              !Box-Mueller sample velocities with st.dev=1
              ran_norm = boxmueller_transform(cx_ran(2:5))
              !> v_new = sqrt(kT/m) * ran_norm
#ifdef WITH_TiTe
              v_new = sqrt(T_i * K_BOLTZ/(sim%groups(group_num)%mass * ATOMIC_MASS_UNIT))*ran_norm(2:4)
#else
              v_new = sqrt(T_e * K_BOLTZ/(sim%groups(group_num)%mass * ATOMIC_MASS_UNIT))*ran_norm(2:4)
#endif
              !>add bulk fluid flow
              v_new = v_new + vvector 

              CX_source = particle_tmp%weight
              !> Compute kinetic energy change to plasma due to velocity change in this substep
              delta_E_kin = 0.5d0 * sim%groups(group_num)%mass * ATOMIC_MASS_UNIT * particle_tmp%weight &
                        * (dot_product(v_old,v_old) - dot_product(v_new,v_new))
              cx_energy   = delta_E_kin  !< for diagnostics only
            endif ! cx_ran
          endif ! CHARGE EXCHANGE

          
          !> check that the energy feedback is valid
          if (isnan(ionize_source * ionize_energy + delta_E_kin - line_rad_energy)) then
            write(*,*) "ionize_energy", ionize_energy
            write(*,*) "delta_E_kin", delta_E_kin
            write(*,*) "line_rad_energy", line_rad_energy
            particle_tmp%i_elm  = 0
            CYCLE !< don't feed this particle into the feedback
          endif

          !> ----- CONSTRUCT FEEDBACK -----
          !> the feedback per particle per time step is accumulated which is then divided by gather time later
          density_source = ionize_source * sim%groups(group_num)%mass * ATOMIC_MASS_UNIT !< mass source in SI
          mom_par_source = ionize_source * dot_product(B, v_old) * sim%groups(group_num)%mass * ATOMIC_MASS_UNIT & 
                + CX_source  * dot_product(B, v_old - v_new) * sim%groups(group_num)%mass * ATOMIC_MASS_UNIT 
#ifdef WITH_TiTe
          energy_source_Te = -ionize_source * H_binding_energy - line_rad_energy
          energy_source_Ti = ionize_source * kinetic_energy + delta_E_kin
#else
          energy_source  = ionize_source * ionize_energy + delta_E_kin - line_rad_energy
#endif
          n_lost_ion = n_lost_ion + ionize_source	!< local sum #particles lost due to ionisation
          p_lost_ion = p_lost_ion + ionize_source * ionize_energy
          p_lost_plt = p_lost_plt + line_rad_energy
          p_lost_cx  = p_lost_cx + cx_source * cx_energy


          !> Calculate the projection of the ion source in real-time
          call basisfunctions(particle_tmp%st(1), particle_tmp%st(2), HH, HH_s, HH_t)
          call mode_moivre(particle_tmp%x(3), HZ)
                
          do l=1,n_vertex_max
            do m=1,n_order+1
  
              density_fb = HH(l,m) * sim%fields%element_list%element(i_elm_old)%size(l,m) * density_source   * t_norm / rho_norm
              mom_par_fb = HH(l,m) * sim%fields%element_list%element(i_elm_old)%size(l,m) * mom_par_source   * t_norm / m_norm
#ifdef WITH_TiTe
              E_fb_Te    = HH(l,m) * sim%fields%element_list%element(i_elm_old)%size(l,m) * energy_source_Te * t_norm / E_norm
              E_fb_Ti    = HH(l,m) * sim%fields%element_list%element(i_elm_old)%size(l,m) * energy_source_Ti * t_norm / E_norm
#else
              E_fb       = HH(l,m) * sim%fields%element_list%element(i_elm_old)%size(l,m) * energy_source    * t_norm / E_norm
#endif
              extra_proj = HH(l,m) * sim%fields%element_list%element(i_elm_old)%size(l,m) * particle_tmp%weight
  
              do i_tor=1,n_tor
                feedback_rhs(m,l,i_elm_old,i_tor,rho_idx_kin)     = feedback_rhs(m,l,i_elm_old,i_tor,rho_idx_kin)  + HZ(i_tor) * density_fb
#ifdef WITH_TiTe
                feedback_rhs(m,l,i_elm_old,i_tor,E_Te_idx_kin)    = feedback_rhs(m,l,i_elm_old,i_tor,E_Te_idx_kin) + HZ(i_tor) * E_fb_Te
                feedback_rhs(m,l,i_elm_old,i_tor,E_Ti_idx_kin)    = feedback_rhs(m,l,i_elm_old,i_tor,E_Ti_idx_kin) + HZ(i_tor) * E_fb_Ti
#else
                feedback_rhs(m,l,i_elm_old,i_tor,E_idx_kin)       = feedback_rhs(m,l,i_elm_old,i_tor,E_idx_kin)    + HZ(i_tor) * E_fb
#endif
                feedback_rhs(m,l,i_elm_old,i_tor,mom_par_idx_kin) = feedback_rhs(m,l,i_elm_old,i_tor,mom_par_idx_kin) + HZ(i_tor) * mom_par_fb
                feedback_rhs(m,l,i_elm_old,i_tor,6)               = feedback_rhs(m,l,i_elm_old,i_tor,6) + HZ(i_tor) * extra_proj
              enddo
            enddo
          enddo

        endif ! END OF NCS SPECIFIC PHYSICS

        ! ============================================ ICS SPECIFIC PHYSICS ===========================================

        if (sim%groups(group_num)%coupling_scheme == 'ics') then
          line_rad_energy = 0.d0
          delta_E_kin = 0.d0
          !> IONISATION & RECOMBINATION (Impurities)
          if (sim%groups(group_num)%use_kin_ionisation .and. .not. limits) then
            call rng(i_rng)%next(ionize_ran_imp)

            !> determines the new charge of the particle using ionisation/recombination rate coefficients
            particle_tmp%q = int(new_charge(int(q_old,4), sim%groups(group_num)%ad, log10(n_e), log10(T_e), tstep_part_adj, ionize_ran_imp(1:2)),1)
            
            if (particle_tmp%q .gt. q_old) then
              binding_energy = sim%groups(group_num)%ad%ionisation_energy(q_old+1) * EL_CHG
              ionize_energy     =  - binding_energy * particle_tmp%weight
              !< including binding energy will make ionize_energy negative, so it becomes a sink for the plasma
            endif
          endif ! IONISATION

          !> RADIATION (Line radiation + Bremsstrahlung + Recombination radiation)
          if (sim%groups(group_num)%use_kin_radiation .and. .not. limits) then !< before or after Ionisation and CX ??
            call sim%groups(group_num)%ad%PLT%interp(int(particle_tmp%q), log10(n_e), log10(T_e), PLT)  ! [J m^3/s] Line radiation
            call sim%groups(group_num)%ad%PRB%interp(int(particle_tmp%q), log10(n_e), log10(T_e), PRB)  ! [J m^3/s] Bremsstrahlung
            call sim%groups(group_num)%ad%ACD%interp(int(particle_tmp%q), log10(n_e), log10(T_e), Srec) ! [J m^3/s] Recomb radiation 
            binding_energy = sim%groups(group_num)%ad%ionisation_energy(particle_tmp%q) * EL_CHG
            radiation_energy = - n_e * particle_tmp%weight * (PLT +PRB-Srec*binding_energy)* tstep_part_adj
          endif ! RADIATION
          
          !> COLLISIONS WITH THE BACKGROUND PLASMA (Neoclassical collisions)
          if (sim%groups(group_num)%use_kin_bg_collisions .and. .not. limits_coll) then
            if (particle_tmp%q .gt. 0) then
              ! Calculate collisions
#ifdef WITH_TiTe
              kTb = T_i*K_BOLTZ
#else
              kTb = T_e*K_BOLTZ ! assume T_e == T_i
#endif
              n_b = n_i
              ! Assumes deuterium background
              q_b = 1
              m_b = central_mass

              if (sim%groups(group_num)%kin_bg_coll_type .eq. 'Homma2013') then
                !> Homma use temperature in [J] (kb [j/K]* T_e [K] or e [J/eV] * Te_eV [eV])
                q = q_homma2013(kTb, grad_T_i*K_BOLTZ, B, n_b, m_b, q_b)
              elseif (sim%groups(group_num)%kin_bg_coll_type .eq. 'Homma2020') then
                q = q_homma2020(kTb, grad_T_i*K_BOLTZ, B, n_b, m_b, q_b, sim%groups(group_num)%homma2020_alpha)
              endif

              !> Calculate coulomb logarithm and limit it to reasonable values
              coulomb_log = coulomb_logarithm(kTb, n_b, particle_tmp%q, q_b, sim%groups(group_num)%mass, m_b)
              coulomb_log = max(10.d0, coulomb_log)
              coulomb_log = min(20.d0, coulomb_log)

              !> Get parallel flow velocity
              call sim%fields%interp_PRZ(t, particle_tmp%i_elm, [var_Vpar], 1, particle_tmp%st(1), particle_tmp%st(2), &
                  particle_tmp%x(3), P, P_s, P_t, P_phi, P_time, R, R_s, R_t, Z, Z_s, Z_t)
              
              do l=1,n_coll
                call rng(i_rng)%next(ran2(:,l))
              end do

              call sample_velocity_dist_magnetized(n_coll, ran2(1:6,:), kTb, q, n_b, m_b, q_b, P(1)*B/sim%t_norm, v_b)
  
              do l=1,n_coll
                call rng(i_rng)%next(ran)
                call collide_particles(ran(1:3), particle_tmp%q, sim%groups(group_num)%mass, v_new, &
                    q_b, m_b, v_b(:,l), n_b, coulomb_log, tstep_part_adj/real(n_coll,8))
              end do
            end if
          endif ! COLLISIONS
          
          !> Kinetic energy change transferred to plasma fluid this substep
          delta_E_kin = 0.5d0 * particle_tmp%weight * sim%groups(group_num)%mass * ATOMIC_MASS_UNIT &
                        * (dot_product(v_old, v_old) - dot_product(v_new, v_new))

          !> check that the particle energy sources are valid
          if (isnan(imp_charge_density + ionize_energy + radiation_energy + delta_E_kin)) then
            write(*,*) "imp_charge_density", imp_charge_density
            write(*,*) "ionize_energy", ionize_energy
            write(*,*) "rad_energy", radiation_energy
            write(*,*) "delta_E_kin", delta_E_kin
            particle_tmp%i_elm  = 0
            CYCLE !< don't feed this particle into the feedback
          endif
      
          !> ----- CONSTRUCT FEEDBACK -----
          !> the feedback per particle per time step is accumulated which is then divided by gather time later
#ifdef WITH_TiTe
          energy_source_Te = ionize_energy + radiation_energy
          energy_source_Ti = delta_E_kin
#else
          energy_source    = ionize_energy + radiation_energy + delta_E_kin
#endif
          mom_par_source   = particle_tmp%weight * dot_product(B, (v_old - v_new)) * sim%groups(group_num)%mass * ATOMIC_MASS_UNIT ! parallel momentum given to plasma

          n_lost_ion = n_lost_ion
          p_lost_ion = p_lost_ion + ionize_energy
          p_lost_plt = p_lost_plt + radiation_energy
          p_lost_cx  = p_lost_cx + cx_source * cx_energy

          !> Calculate the projection of the ion source in real-time
          call basisfunctions(particle_tmp%st(1), particle_tmp%st(2), HH, HH_s, HH_t)
          call mode_moivre(particle_tmp%x(3), HZ)
          
          do l=1,n_vertex_max
            do m=1,n_order+1
  
              mom_par_fb     = HH(l,m) * sim%fields%element_list%element(i_elm_old)%size(l,m) * mom_par_source    * t_norm / m_norm
#ifdef WITH_TiTe
              E_fb_Te        = HH(l,m) * sim%fields%element_list%element(i_elm_old)%size(l,m) * energy_source_Te  * t_norm / E_norm
              E_fb_Ti        = HH(l,m) * sim%fields%element_list%element(i_elm_old)%size(l,m) * energy_source_Ti  * t_norm / E_norm
#else
              E_fb           = HH(l,m) * sim%fields%element_list%element(i_elm_old)%size(l,m) * energy_source     * t_norm / E_norm
#endif
              imp_q_fb       = HH(l,m) * sim%fields%element_list%element(i_elm_old)%size(l,m) * particle_tmp%weight * particle_tmp%q
              imp_density_fb = HH(l,m) * sim%fields%element_list%element(i_elm_old)%size(l,m) * particle_tmp%weight
              imp_P_rad_fb   = HH(l,m) * sim%fields%element_list%element(i_elm_old)%size(l,m) * radiation_energy / tstep_part_adj
              do i_tor=1,n_tor
#ifdef WITH_TiTe
                feedback_rhs(m,l,i_elm_old,i_tor,E_Te_idx_kin)    = feedback_rhs(m,l,i_elm_old,i_tor,E_Te_idx_kin)    + HZ(i_tor) * E_fb_Te
                feedback_rhs(m,l,i_elm_old,i_tor,E_Ti_idx_kin)    = feedback_rhs(m,l,i_elm_old,i_tor,E_Ti_idx_kin)    + HZ(i_tor) * E_fb_Ti
#else
                feedback_rhs(m,l,i_elm_old,i_tor,E_idx_kin)       = feedback_rhs(m,l,i_elm_old,i_tor,E_idx_kin)       + HZ(i_tor) * E_fb
#endif
                feedback_rhs(m,l,i_elm_old,i_tor,mom_par_idx_kin) = feedback_rhs(m,l,i_elm_old,i_tor,mom_par_idx_kin) + HZ(i_tor) * mom_par_fb
                feedback_rhs(m,l,i_elm_old,i_tor,imp_q_idx)       = feedback_rhs(m,l,i_elm_old,i_tor,imp_q_idx)       + HZ(i_tor) * imp_q_fb       ! impurity charge density
                feedback_rhs(m,l,i_elm_old,i_tor,7)               = feedback_rhs(m,l,i_elm_old,i_tor,7)               + HZ(i_tor) * imp_P_rad_fb   ! impurity radiated power [to be moved to diag feedback]
                feedback_rhs(m,l,i_elm_old,i_tor,8)               = feedback_rhs(m,l,i_elm_old,i_tor,8)               + HZ(i_tor) * imp_density_fb ! impurity density [to be moved to diag feedback]
              enddo
            enddo
          enddo

        endif ! END OF ICS SPECIFIC PHYSICS

        !> explicitly store updated velocity
        particle_tmp%v = v_new

        !> =============================== PUSH PARTICLE ====================================
        if (particle_tmp%i_elm .gt. 0) then
          call boris_push_cylindrical(particle_tmp, sim%groups(group_num)%mass, E, B, tstep_part_adj)
          call find_RZ_nearby(sim%fields%node_list, sim%fields%element_list, rz_old(1), rz_old(2), st_old(1), st_old(2), i_elm_old, &
                              particle_tmp%x(1), particle_tmp%x(2), particle_tmp%st(1), particle_tmp%st(2), particle_tmp%i_elm, ifail)
        end if
    
      end do ! steps 

      sim%groups(group_num)%particles(j) = particle_tmp
    end do   ! particles
    !$omp end parallel do
      
    !> ===================================== DIAGNOSTICS ===================================
    call MPI_REDUCE(n_lost_ion, n_lost_ion_all, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call MPI_REDUCE(p_lost_ion, p_lost_ion_all, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call MPI_REDUCE(p_lost_plt, p_lost_plt_all, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call MPI_REDUCE(p_lost_cx, p_lost_cx_all, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call MPI_REDUCE(n_super_ionized, n_super_ionized_all, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    
    if (sim%my_id .eq. 0) write(*,'(A46,E14.6,I6)') "Lost superparticles at t due to ionisation: ", sim%time, n_super_ionized_all
    if (sim%my_id .eq. 0) write(*,'(A46,2E14.6)') "Lost particles at t due to ionisation: ", sim%time, n_lost_ion_all
    p_lost_ion_all = p_lost_ion_all / (tstep_part_adj * nstep_part_adj)
    p_lost_plt_all = p_lost_plt_all / (tstep_part_adj * nstep_part_adj)
    p_lost_cx_all  = p_lost_cx_all / (tstep_part_adj * nstep_part_adj)
    if (sim%my_id .eq. 0) write(*,'(A46,2E14.6)') "Lost energy [W] at t due to ionisation: ", sim%time, p_lost_ion_all
    if (sim%my_id .eq. 0) write(*,'(A46,2E14.6)') "Lost energy [W] at t due to radiation: ", sim%time, p_lost_plt_all
    if (sim%my_id .eq. 0) write(*,'(A46,2E14.6)') "Lost energy [W] at t due to CX radiation: ", sim%time, p_lost_cx_all
    if (sim%my_id .eq. 0) write(*,'(A46,2E14.6)') "Total energy exchange to plasma [W]: ", sim%time, p_lost_ion_all -p_lost_plt_all+ p_lost_cx_all
  end subroutine evolve_ncs_ics


  subroutine evolve_epf(sim, group_num, feedback_rhs, rng, tstep_part_adj, nstep_part_adj)
    use mod_project_particles
    use mod_random_seed
    use mod_interp, only: mode_moivre
    use mod_basisfunctions

    implicit none
    class(particle_sim), target, intent(inout)                :: sim
    integer, intent(in)                                       :: group_num
    real*8, allocatable, intent(inout)                        :: feedback_rhs(:,:,:,:,:)
    type(count_action)                                        :: counter
    type(pcg32_rng), dimension(:), allocatable, intent(inout) :: rng
    real*8, intent(in)                                        :: tstep_part_adj
    integer, intent(in)                                       :: nstep_part_adj

    !> local variables
    real*8   :: t, E(3), B(3), B_norm(3), psi, U
    real*8   :: rzp_old(3), st_old(2)
    real*8   :: v_tilde_r, v_tilde_z, v_par, p_par, v, p_perp, p_atrop, base
    real*8   :: HZ(n_tor), HH(4,4), HH_s(4,4), HH_t(4,4) !> Bezier basis functions
    integer  :: iterations                               !> number of tsteps for which projection quantities are collected

    integer  :: i, j, k, l, m, ifail, i_tor
    integer  :: i_elm, i_elm_old

    if (nstep_part_adj < proj_collection_period) then
      proj_collection_period = nstep_part_adj
    endif

    !> loop over all particles in group(group_num)
    select type (particles => sim%groups(group_num)%particles)
    type is (particle_kinetic_leapfrog)
      if (use_manual_random_seed) then
        !$ call omp_set_schedule(omp_sched_static,10)
      else
        !$ call omp_set_schedule(omp_sched_dynamic,10)
      end if
      !$omp parallel do default(none) &
      !$omp schedule(runtime) &
      !$omp private(j, k, l, m, HZ, HH, HH_s, HH_t, E, B, psi, U, rzp_old, st_old, i_elm_old, &
      !$omp i_elm, B_norm, v_tilde_r, v_tilde_z, v_par, p_perp, p_par, p_atrop, base, v, i_tor, ifail) &
      !$omp shared(nstep_part_adj, tstep_part_adj, sim, group_num, proj_collection_period, &
      !$omp PI_RR_idx_kin, PI_ZZ_idx_kin, PI_PHIPHI_idx_kin, PI_RZ_idx_kin, PI_RPHI_idx_kin, PI_ZPHI_idx_kin, rho_ep_idx_kin) &
      !$omp reduction(+:feedback_rhs)

      do j=1,size(particles,1)
        do k=1,nstep_part_adj

          !> if particle is lost, skip
          if (particles(j)%i_elm .le. 0) exit

          !> Determine E, B at particles location
          call sim%fields%calc_EBpsiU(sim%time, particles(j)%i_elm, particles(j)%st, particles(j)%x(3), E, B, psi, U)

          rzp_old   = particles(j)%x
          st_old    = particles(j)%st
          i_elm_old = particles(j)%i_elm

          !> push particles and find (R,Z)
          if (particles(j)%i_elm .gt. 0) then
            call boris_push_cylindrical(particles(j), sim%groups(group_num)%mass, E, B, tstep_part_adj)
            call find_RZ_nearby(sim%fields%node_list, sim%fields%element_list, rzp_old(1), rzp_old(2), st_old(1), st_old(2), i_elm_old, particles(j)%x(1), particles(j)%x(2), particles(j)%st(1), particles(j)%st(2), particles(j)%i_elm, ifail)
          endif

          !> may have pushed particle out of domain, check again if it is lost
          if (particles(j)%i_elm .le. 0) exit

          !> only collect projections every proj_collection_period number of timesteps
          if (mod(k,proj_collection_period) .ne. 0) cycle

          !> calc normalised B and orthonormal v cmpts
          B_norm    = B / norm2(B)
          v_par     = dot_product(B_norm, particles(j)%v)
          v_tilde_r = (-(B_norm(1))*particles(j)%v(3)+B_norm(3)*particles(j)%v(1))/sqrt(B_norm(1)**2+B_norm(3)**2)
          v_tilde_z = (particles(j)%v(2)-B_norm(2)*v_par)/sqrt(B_norm(1)**2+B_norm(3)**2)

          !> calc parallel and perpendicular pressures
          p_perp  = 1.d0/2.d0 * (v_tilde_r**2 + v_tilde_z**2)
          p_par   = v_par**2
          p_atrop = p_par - p_perp

          !> calc FEM basis functions
          call basisfunctions(particles(j)%st(1), particles(j)%st(2), HH, HH_s, HH_t)
          call mode_moivre(particles(j)%x(3), HZ)

          !> Gather particle pressures for jorek feedback
          i_elm = particles(j)%i_elm
          do l = 1,n_vertex_max
            do m = 1,n_order+1

              v = HH(l,m) * sim%fields%element_list%element(i_elm)%size(l,m)

              do i_tor = 1,n_tor
                !> expression used in all of the following
                base = HZ(i_tor)*v*particles(j)%weight*sim%groups(group_num)%mass*ATOMIC_MASS_UNIT*MU_ZERO

                !> PI_RR
                feedback_rhs(m,l,i_elm,i_tor,PI_RR_idx_kin) = feedback_rhs(m,l,i_elm,i_tor,PI_RR_idx_kin) + base*(p_perp+B_norm(1)**2*p_atrop)

                !> PI_ZZ
                feedback_rhs(m,l,i_elm,i_tor,PI_ZZ_idx_kin) = feedback_rhs(m,l,i_elm,i_tor,PI_ZZ_idx_kin) + base*(p_perp+B_norm(2)**2*p_atrop)

                !> PI_PHIPHI
                feedback_rhs(m,l,i_elm,i_tor,PI_PHIPHI_idx_kin) = feedback_rhs(m,l,i_elm,i_tor,PI_PHIPHI_idx_kin) + base*(p_perp+B_norm(3)**2*p_atrop)

                !> PI_RZ
                feedback_rhs(m,l,i_elm,i_tor,PI_RZ_idx_kin) = feedback_rhs(m,l,i_elm,i_tor,PI_RZ_idx_kin) + base*(B_norm(1)*B_norm(2)*p_atrop)

                !> PI_RPHI
                feedback_rhs(m,l,i_elm,i_tor,PI_RPHI_idx_kin) = feedback_rhs(m,l,i_elm,i_tor,PI_RPHI_idx_kin) + base*(B_norm(1)*B_norm(3)*p_atrop)

                !> PI_ZPHI
                feedback_rhs(m,l,i_elm,i_tor,PI_ZPHI_idx_kin) = feedback_rhs(m,l,i_elm,i_tor,PI_ZPHI_idx_kin) + base*(B_norm(2)*B_norm(3)*p_atrop)

                !> density
                feedback_rhs(m,l,i_elm,i_tor,rho_ep_idx_kin) = feedback_rhs(m,l,i_elm,i_tor,rho_ep_idx_kin) + HZ(i_tor)*v*particles(j)%weight

              end do !> toroidal harmonics
            end do   !> FEM basis order
          end do     !> verticies of element
        end do       !> particle timesteps
      end do         !> particles


      !$omp end parallel do
    end select

    !> Renormalise by number of timesteps the projection quantities were collected for. Need integer division
    !> Ie if nstep = 201, proj_col_period = 10, then iterations = 20 (we only collected projections 20 times)
    iterations   = nstep_part_adj / proj_collection_period
    feedback_rhs = feedback_rhs/iterations

  end subroutine evolve_epf

end module mod_particle_evolution
