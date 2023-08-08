!> Initialise input parameters and read the input namelist
subroutine initialise_parameters(my_id, filename)

use tr_module
use phys_module
use vacuum
use live_data
use pellet_module

implicit none

! --- Routine parameters
integer,                      intent(in) :: my_id
character(len=*),             intent(in) :: filename

! --- Local variables
integer :: ierr, err, i

! --- Namelist with input parameters.
namelist /in1/  tstep, nstep, tstep_n, nstep_n,                     &
                rst_hdf5, rst_hdf5_version, keep_current_prof,      &
                restart, regrid, write_ps, time_evol_theta,         &
                time_evol_zeta, force_horizontal_Xline,             &
                Mach1_openBC, thermalization, Mach1_fix_B,          &
                eta_ARAZ_const, eta_ARAZ_on, eta_ARAZ_simple,       & 
                tauIC_ARAZ_on,                                      &
                n_tor_fft_thresh, fix_axis_nodes,                   &
                n_R, n_Z, n_radial, n_pol, n_tht, n_flux,           &
                n_open, n_private, n_leg,                           &
                n_outer, n_inner, n_up_priv, n_up_leg,              &
                psi_axis_init, XR_r, SIG_r, XR_tht, SIG_tht,        &
                SIG_closed, SIG_open, SIG_private, SIG_theta,       &
                SIG_leg_0, SIG_leg_1, dPSI_open, dPSI_private,      &
                SIG_up_leg_0, SIG_up_leg_1, SIG_up_priv,            &
                SIG_outer, SIG_inner, SIG_theta_up,                 &
                dPSI_outer, dPSI_inner, dPSI_up_priv,               &
                nout, nout_projection, xr1, sig1, xr2, sig2,        &
                R_begin, R_end, Z_begin, Z_end,                     &
                R_geo, Z_geo, amin, mf, fbnd, fpsi, mode,           &
                R_Z_psi_bnd_file,                                   &
                R_boundary, Z_boundary, psi_boundary, n_boundary,   &
                extend_existing_grid, no_mach1_bc,                  &
                grid_to_wall, RZ_grid_inside_wall, eqdsk_psi_fact,  &
                RZ_grid_jump_thres, n_tht_equidistant,              &
                n_ext_equidistant,                                  &
                n_wall_blocks, n_ext_block, corner_block,           &
                n_block_points_left,  n_block_points_right,         &
                R_block_points_left,  R_block_points_right,         &
                Z_block_points_left,  Z_block_points_right,         &
                use_simple_bnd_types,                               &
                tokamak_device, manipulate_psi_map,                 &
                F0, gamma, gamma_stangeby,                          &
                gamma_sheath_i, gamma_sheath_e,                     &
                zjz_0, zjz_1, zj_coef,                              &
                rho_0, rho_1, rho_coef, rho_min,                    &
                T_0,   T_1,   T_coef, T_min,                        &
				T_min_neg,rho_min_neg,                              &
				corr_neg_temp_coef, corr_neg_dens_coef,             &
                Ti_0,  Ti_1,  Ti_coef,                              &
                Te_0,  Te_1,  Te_coef,                              &
                FF_0,  FF_1,  FF_coef,                              &
                V_0, V_1, V_coef,                                   &
                ZK_par, ZK_i_par, ZK_e_par, ZK_par_max,             &
                ZK_perp, ZK_i_perp, ZK_e_perp, D_par, D_perp,       &
                heatsource_e, heatsource_i,                         &
                eta, visco, visco_par, ZK_i_perp_num, ZK_e_perp_num,&
                eta_num, visco_num, visco_par_num, D_perp_num,      &
                Dn_perp_num,                                        &
                heatsource_psin, heatsource_sig,                    &
                particlesource_psin, particlesource_sig,            &
                edgeparticlesource, edgeparticlesource_psin,        &
                edgeparticlesource_sig,                             &
                particlesource_gauss, heatsource_gauss,             &
                heatsource_gauss_i, heatsource_gauss_e,             &
                heatsource_gauss_psin, heatsource_gauss_sig,        &
                particlesource_gauss_psin, particlesource_gauss_sig,&
                particlesource, heatsource, tauIC,                  &
                pellet_amplitude, pellet_R, pellet_Z, pellet_phi,   &
                pellet_radius, pellet_sig, pellet_length,           &
                pellet_psi, pellet_delta_psi,                       &
                central_density, central_mass,                      &
                pellet_particles, use_pellet,                       &
                ellip,tria_u,tria_l,quad_u,quad_l,                  &
                xampl,xwidth,xsig,xtheta,xshift,xleft, xpoint,      &
                xcase, time_evol_scheme,                            &
                freeboundary_equil,                                 &
                rho_file, T_file, ffprime_file, Fprofile_file,      &
                Ti_file, Te_file,                                   &
                bc_natural_open, bc_natural_flux, gamma_sheath,     &
                freeboundary, resistive_wall, freeb_change_indices, &
                use_mumps_eq, use_pastix_eq, use_strumpack_eq,      &
                use_mumps, mumps_ordering, use_strumpack,           &
                use_BLR_compression, epsilon_BLR, just_in_time_BLR, &
                use_pastix, use_murge, use_murge_element,           &
                refinement, grid_to_wall,                           &
                fix_axis_nodes,                                     &
                adaptive_time, equil, bench_without_plot,           &
                eta_T_dependent, visco_T_dependent,ZKpar_T_dependent,&
                heatsource_psin, heatsource_sig,                    &
                particlesource_psin, particlesource_sig,            &
                edgeparticlesource, edgeparticlesource_psin,        &
                edgeparticlesource_sig,                             &
                neutral_line_source,                                &
                neutral_line_R_start, neutral_line_Z_start,         &
                neutral_line_R_end,   neutral_line_Z_end,           &
                produce_live_data, gmres, gmres_max_iter,           &
                iter_precon, gmres_4, gmres_m, gmres_tol,           &
                max_steps_noUpdate,                                 &
                keep_n0_const, linear_run, export_for_nemec,        &
                output_bnd_elements,                                &
                wall_file, Z_xpoint_limit,                          &
                first_target_point, last_target_point,              &
                n_limiter, R_limiter, Z_limiter,                    &
                spi_tor_rot, tor_frequency, spi_num_vol,            &
                NEO, neo_file, aki_neo_const, amu_neo_const,        &
                D_neutral_x, D_neutral_y, D_neutral_p,              &
                neutral_reflection,                                 &
                ns_deltaphi, ns_delta_minor_rad, ksi_ion, spi_rnd_seed, &
                ns_amplitude, ns_R, ns_Z, ns_phi, ns_radius,        &
                spi_Vel_Rref,spi_Vel_Zref, using_spi, n_spi, n_inj, &
                spi_Vel_RxZref, spi_quantity, spi_abl_model,        &
                ns_radius_ratio, ns_radius_min, spi_angle,          &
                spi_L_inj, spi_L_inj_diff,                          &
                drift_distance, energy_teleported,                  &
                K_Dmv, A_Dmv, L_tube, V_Dmv, P_Dmv,                 &
                spi_Vel_diff, t_ns, JET_MGI, ASDEX_MGI,             &
                delta_n_convection, nimp_bg,                        &
                amix, amix_freeb, equil_accuracy,                   &
                equil_accuracy_freeb, current_ref, FB_Ip_position,  &
                FB_Ip_integral, Z_axis_ref, FB_Zaxis_position,      &
                FB_Zaxis_derivative,FB_Zaxis_integral, start_VFB,   &
                n_feedback_current, n_feedback_vertical,            &
                n_iter_freeb, n_pf_coils, pf_coils, R_axis_ref,     &
                axis_srch_radius, rad_FB_amp,                       &
                starwall_equil_coils, freeb_equil_iterate_area,     &
                psi_offset_freeb, diag_coils, rmp_coils,            &
                voltage_coils, vert_FB_amp, find_pf_coil_currents,  &
                delta_psi_GS, newton_GS_fixbnd, newton_GS_freebnd,  &
                pastix_maxthrd, centralize_harm_mat, treat_axis,    &
                n_particles, tstep_particles, nstep_particles,      & !Particles extension
                nsubstep_particles, restart_particles,              &
                filter_perp,    filter_hyper,    filter_par,        &
                filter_perp_n0, filter_hyper_n0, filter_par_n0,     &
                use_cx, use_sputtering, use_ionisation,             &
                use_ncs, use_pcs, use_ccs, cte_current_FB_fact,     &
                eta_ohmic, force_central_node, D_par_imp,           &
                D_perp_imp, imp_type, adas_dir,                     &
                spi_shard_file, spi_tor_rot, pellet_density,        &
                pellet_density_bg, cte_current_FB_fact, Te_eV_min,  &
                spi_plume_file, n_adas, index_main_imp,             &
                output_prad_phi, spi_quantity_bg, CARIDDI_mode,     &
                autodistribute_modes, modes_per_family,             &
                mode_families_modes, n_mode_families,               &
                weights_per_family, autodistribute_ranks,           &
                ranks_per_family,                                   &                
                use_sc, add_sources_in_sc, visco_sc_num,            &
                D_perp_sc_num, D_par_sc_num, ZK_perp_sc_num,        &
                ZK_par_sc_num, ZK_i_perp_sc_num, ZK_i_par_sc_num,   &
                ZK_e_perp_sc_num, ZK_e_par_sc_num, visco_par_sc_num,&
                Dn_pol_sc_num, Dn_p_sc_num, D_perp_imp_sc_num,      &
                D_par_imp_sc_num, use_vms,                          &
                vms_coeff_AR, vms_coeff_AZ, vms_coeff_A3,           &
                vms_coeff_UR, vms_coeff_UZ, vms_coeff_Up,           &
                vms_coeff_rho, vms_coeff_Ti, vms_coeff_Te,          &
                vms_coeff_T, vms_coeff_rhon, vms_coeff_rhoimp,      &
                vacuum_min, strumpack_matching

if (my_id .eq. 0) then

  ! --- Preset input parameters to reasonable default values.
  call preset_parameters()
  call vacuum_preset(my_id, freeboundary_equil, freeboundary, resistive_wall)

  ! --- Model-specific presets
  ! -none-
  
  ! --- Read input parameters from namelist.
  if (trim(filename) .ne. "__NO_FILENAME__" ) then
    open(42, file=filename, status='old', action='read', iostat=ierr)
    if ( ierr /= 0 ) then
      write(*,*) 'ERROR: COULD NOT OPEN NAMELIST FILE "', trim(filename), '".'
      stop
    end if
    read(42,in1)
    close(42)
  else
    read(5,in1)
  end if

  !==============================R_Z_psi_bnd==========================
  if ( (n_boundary.ne.0) .and. (R_Z_psi_bnd_file /= 'none') ) then
    ! --- Open the file.
    OPEN(UNIT=243, FILE=R_Z_psi_bnd_file, FORM='FORMATTED', STATUS='OLD', ACTION='READ', IOSTAT=err)
    if ( err /= 0 ) then
      write(*,*) 'ERROR in initialise_parameters: Cannot open file '//TRIM(R_Z_psi_bnd_file)//'.'
      write(*,*) 'Assuming data is in main input file '//TRIM(filename)//'.'
    else
      write(*,'(A)') ' boundary info from R_Z_psi_bnd_file: R_boundary, Z_boundary, psi_boundary ' 
      do i=1,n_boundary
        read(243,*) R_boundary(i),Z_boundary(i),psi_boundary(i)
        write(*,*) R_boundary(i),Z_boundary(i),psi_boundary(i)  
      enddo
    endif    
    CLOSE(243)
  endif
  !=========================================
  
  !==============================Limiter==========================
  if (n_limiter.ne.0) then
    ! --- Open the file.
    OPEN(UNIT=244, FILE=wall_file, FORM='FORMATTED', STATUS='OLD', ACTION='READ', IOSTAT=err)
    if ( err /= 0 ) then
      write(*,*) 'ERROR in initialise_parameters: Cannot open file '//TRIM(wall_file)//'.'
      write(*,*) 'Assuming data is in main input file '//TRIM(filename)//'.'
    else
      write(*,'(A)') ' wall info from wall_file: R_wall, Z_wall ' 
      do i=1,n_limiter
        read(244,*) R_limiter(i),Z_limiter(i)
        write(*,*)  R_limiter(i),Z_limiter(i)
      enddo
    endif    
    CLOSE(244)
  endif
  
  ! --- Calculate normalisation factor for MGI source (related to its toroidal shape)
  ns_tor_norm = ns_deltaphi * PI**0.5 * ERF(PI/ns_deltaphi)

  ! --- Calculate JOREK gamma_sheath from gamma_stangeby if provided (otherwise the other way around)
  if (gamma_e_stangeby > -1.d89) then
    gamma_sheath_e = (gamma-1.d0) * (gamma_e_stangeby - 1.d0)
  else
    gamma_e_stangeby = gamma_sheath_e / (gamma-1.d0) + 1.d0
  end if
  if (gamma_i_stangeby > -1.d89) then
    gamma_sheath_i = (gamma-1.d0) * (gamma_i_stangeby - 1.d0 - gamma)
  else
    gamma_i_stangeby = gamma_sheath_i / (gamma-1.d0) + 1.d0 + gamma
  end if

  if (sum(nstep_n) .gt. 0) then
    nstep = sum(nstep_n)
  else
    tstep_n    = 0.d0
    tstep_n(1) = tstep
    nstep_n    = 0
    nstep_n(1) = nstep
  endif

  ! --- Checking consistency of eta_ARAZ parameters
  if (eta_ARAZ_on .eqv. .true.) then
     if (eta_ARAZ_const .ne. 0) then
        write(*,*) 'One should not use both eta_ARAZ_on and eta_ARAZ_const simultaneously, to avoid double-counting. Please use eta_ARAZ_on = .t. with eta_ARAZ_const = 0.d0, or eta_ARAZ_on = .f. with eta_ARAZ_const .ne. 0'
        stop
     endif
  endif

  ! --- Fill the same ablation model to others if not specified to keep the old behavior
  do i = 2,n_inj
    if (spi_abl_model(i) < 0) then
      spi_abl_model(i) = spi_abl_model(1)
    end if
  end do
  
  call allocate_live_data()

endif

keep_n0_const  = ( keep_n0_const .or. linear_run )
! --- Read numerical profiles for rho, T, and ff'.
call read_num_profiles(my_id)

! --- Determine the derivatives of the numerical input profiles.
call derive_num_profiles(my_id)
 
! --- Initialize the shattered pellet position

#if (defined WITH_Neutrals) && (!defined WITH_Impurities)
spi_quantity_bg   = spi_quantity
pellet_density_bg = pellet_density
#endif

if ( my_id == 0 ) then
  if (2*PI/(n_tor*n_period) >= ns_deltaphi) then
    write(*,*) "WARNING! ns_deltaphi too small for the n_tor, BEWARE!"
    if (t_now > minval(t_ns)) then
      write(*,*) "EXITING NOW!!!"
      stop
    end if
  end if


  if (n_inj > n_inj_max .or. n_inj < 1) then
    write(*,*) "ERROR! Do not support n_inj larger than n_inj_max or smaller than 1, EXITING!"
    stop
  end if

    do i = 1, n_inj_max
    if (n_spi(i)/=0 .and. i > n_inj) then
      write(*,*) "ERROR! Something wrong with n_inj, double check, EXITING!", n_spi, n_inj
      stop
    end if
  end do

  if (n_adas > n_imp_max) then
    write(*,*) "ERROR: n_adas should be no larger than n_imp_max, EXITING!"
    stop
  end if

  if (n_adas > 1 .and. (.not. use_imp_adas)) then
    write(*,*) "ERROR: Only support ADAS data for more than one impurities, through setting use_imp_adas to true, EXITING!"
    stop
  end if

  if (index_main_imp < 0 .or. index_main_imp > n_adas) then
    write(*,*) "ERROR: Illegal value of index_main_imp, EXITING!"
    write(*,*) "ERROR: index_main_imp:", index_main_imp
    stop
  end if

  do i = 1,n_inj
    if (drift_distance(i) < 0.d0 .or. energy_teleported(i) < 0.d0) then
      write(*,*) "ERROR: drift_distance and energy_teleported should be 0 or positive as signs already handled in codes, EXITING!"
      stop
    end if
  end do

  if (using_spi) call init_spi_all()
  ! ---- This is an ad hoc way to stop duplication of neutral source and ionized main species source, should find a permanenty solution later

  if (with_impurities .and. with_neutrals .and. using_spi .and. (maxval(spi_quantity_bg)>0.0)) then
    write(*,*) "WARNNING! Currently do not support both with_neutrals and with impurities enabled at the same time"
    write(*,*)  "while injecting background species! Should be fixed soon. EXITING!"
    stop
  endif

endif

return
end subroutine initialise_parameters
