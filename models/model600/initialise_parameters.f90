!> Initialise input parameters and read the input namelist
subroutine initialise_parameters(my_id, filename)

use tr_module
use phys_module
use mumps_module,  only: no_zeros_mumps, mumps_ordering
use pastix_module, only: no_zeros_pastix, pastix_smp_only, pastix_pivot, &
    pastix_maxthrd
use vacuum
use pellet_module
use live_data

implicit none

! --- Routine parameters
integer,                      intent(in) :: my_id
character(len=*),             intent(in) :: filename
real*8 :: vacuum_fraction, b_over_a, a_over_b

! --- Local variables
integer :: ierr,err,i

! --- Namelist with input parameters.
namelist /in1/  tstep, nstep, tstep_n, nstep_n,                     &
                rst_hdf5, rst_hdf5_version, keep_current_prof,      &
                eta, visco, visco_par,                              &
                restart, rst_format, regrid, bootstrap, write_ps,   &
                n_R, n_Z, n_radial, n_pol, n_tht, n_flux,           &
                n_open, n_private, n_leg, n_leg_out, n_ext,         &
                n_outer, n_inner, n_up_priv, n_up_leg, n_up_leg_out,&
                n_tht_equidistant,                                  &
                psi_axis_init, XR_r, SIG_r, XR_tht, SIG_tht,        &
                XR_z, SIG_z, bgf_r, bgf_z,                          &
                SIG_closed, SIG_open, SIG_private, SIG_theta,       &
                SIG_leg_0, SIG_leg_1, dPSI_open, dPSI_private,      &
                SIG_up_leg_0, SIG_up_leg_1, SIG_up_priv,            &
                SIG_outer, SIG_inner,                               &
                dPSI_outer, dPSI_inner, dPSI_up_priv,               &
                nout, xr1, sig1, xr2, sig2,                         &
                R_begin, R_end, Z_begin, Z_end,                     &
                rect_grid_vac_psi,                                  &
                R_geo, Z_geo, amin, mf, fbnd, fpsi, mode,           &
                R_boundary, Z_boundary, psi_boundary, n_boundary,   &
                R_Z_psi_bnd_file,                                   &
                force_horizontal_Xline,                             &
                n_pfc, manipulate_psi_map,                          &
                Rmin_pfc, Rmax_pfc, Zmin_pfc, Zmax_pfc, current_pfc,&
                n_jropes,                                           &
                R_jropes, Z_jropes, w_jropes, current_jropes,       &
                extend_existing_grid, no_mach1_bc,                  &
                grid_to_wall, RZ_grid_inside_wall, eqdsk_psi_fact,  &
                RZ_grid_jump_thres,                                 &
                n_wall_blocks, n_ext_block, corner_block,           &
                n_ext_equidistant,                                  &
                n_block_points_left,  n_block_points_right,         &
                R_block_points_left,  R_block_points_right,         &
                Z_block_points_left,  Z_block_points_right,         &
                use_simple_bnd_types,                               &
                tokamak_device, thermalization,                     &
                F0,                                                 &
                gamma_stangeby,gamma_i_stangeby,gamma_e_stangeby,   &
                gamma_sheath, gamma_sheath_i, gamma_sheath_e,       &
                deuterium_adas, deuterium_adas_1e20,                &
                old_deuterium_atomic,                               &
                density_reflection,                                 &
                mach_one_bnd_integral, Vpar_smoothing,              &
                Vpar_smoothing_coef,                                &
                zjz_0, zjz_1, zj_coef,                              &
                rho_0, rho_1, rho_coef,                             &
                rhon_0, rhon_1, rhon_coef,                          &
                T_0,   T_1,   T_coef,                               &
                Ti_0,  Ti_1,  Ti_coef,                              &
                Te_0,  Te_1,  Te_coef,                              &
                FF_0,  FF_1,  FF_coef,                              &
                ZK_par, ZK_i_par, ZK_e_par, ZK_par_max,             &
                ZK_perp, ZK_i_perp, ZK_e_perp, D_par, D_perp,       &
                heatsource_e, heatsource_i, heatsource,             &
                particlesource, tauIC, Wdia,                        &
                eta_num, visco_num, visco_par_num, D_perp_num,      &
                eta_num_T_dependent, visco_num_T_dependent,         &
                ZK_perp_num, Dn_perp_num, time_evol_scheme,         &
                pellet_amplitude, pellet_R, pellet_Z, pellet_phi,   &
                pellet_radius, pellet_sig, pellet_length,           &
                pellet_psi, pellet_delta_psi, pellet_density,       &
                pellet_velocity_R, pellet_velocity_Z, pellet_theta, &
                pellet_ellipse,                                     &
                central_density, central_mass,                      &
                pellet_particles, use_pellet,                       &
                ellip,tria_u,tria_l,quad_u,quad_l,                  &
                xampl,xwidth,xsig,xtheta,xshift,xleft, xpoint,      &
                xcase, SDN_threshold, D_perp_file, ZK_perp_file,    &
                ZK_e_perp_file, ZK_i_perp_file,                     &
                rho_file, T_file, Ti_file, Te_file, ffprime_file,   &
                rot_file, normalized_velocity_profile,              &
                freeboundary_equil, freeboundary,  freeb_change_indices, &
                resistive_wall,                                     &
                wall_resistivity, wall_resistivity_fact,            &
                bc_natural_open,                                    &
                use_mumps_eq, use_pastix_eq, use_strumpack_eq,      &
                use_mumps, mumps_ordering,                          &
                use_BLR_compression, epsilon_BLR, just_in_time_BLR, &
                use_pastix, use_murge, use_murge_element, use_wsmp, &
                n_tor_fft_thresh, use_strumpack,                    &
                pastix_smp_only, refinement, force_central_node,    &
                fix_axis_nodes,                                     &
                adaptive_time, equil, bench_without_plot,           &
                no_zeros_pastix, no_zeros_mumps,                    &
                eta_T_dependent, visco_T_dependent,                 &
                zkpar_T_dependent, T_max_eta, T_max_eta_ohm,        & 
                heatsource_psin, heatsource_sig,                    &
                heatsource_e_psin, heatsource_e_sig,                &
                heatsource_i_psin, heatsource_i_sig,                &
                particlesource_psin, particlesource_sig,            &
                edgeparticlesource, edgeparticlesource_psin,        &
                edgeparticlesource_sig,                             &
                particlesource_gauss,    heatsource_gauss,          &
                heatsource_gauss_i,      heatsource_gauss_e,        &
                heatsource_gauss_psin,   heatsource_gauss_sig,      &
                heatsource_gauss_i_psin, heatsource_gauss_i_sig,    &
                heatsource_gauss_e_psin, heatsource_gauss_e_sig,    &
                particlesource_gauss_psin, particlesource_gauss_sig,&
                neutral_line_source,                                &
                neutral_line_R_start, neutral_line_Z_start,         &
                neutral_line_R_end,   neutral_line_Z_end,           &
                produce_live_data, gmres, gmres_max_iter,           &
                gmres_m, gmres_4, gmres_tol, iter_precon,           &
                pastix_pivot, max_steps_noUpdate,                   &
                keep_n0_const, linear_run, export_for_nemec,        &
                RMP_on, RMP_har_cos,RMP_har_sin,                    &
                RMP_growth_rate, RMP_ramp_up_time,                  &
                RMP_psi_cos_file, RMP_psi_sin_file,                 &
                V_0,V_1,V_coef, output_bnd_elements,                &
                n_limiter, R_limiter, Z_limiter,                    &
                first_target_point, last_target_point,		    &
                R_Z_psi_bnd_file, wall_file,time_evol_scheme,       &
                spi_tor_rot, tor_frequency, spi_num_vol,            &
                NEO, neo_file, aki_neo_const, amu_neo_const,        &
                D_prof_neg_thresh, ZK_prof_neg_thresh, T_min,       &
                ne_SI_min, Te_eV_min, rn0_min,                      &
                D_neutral_x, D_neutral_y, D_neutral_p,              &
                neutral_reflection, rho_min,                        &
                corr_neg_temp_coef,                                 &
                corr_neg_dens_coef, D_prof_neg, ZK_prof_neg,        &  
                ZK_par_neg,                                         & 
                ns_deltaphi, ns_deltaminrad, ksi_ion, spi_rnd_seed, &
                ns_amplitude, ns_R, ns_Z, ns_phi, ns_radius,        &
                spi_Vel_Rref,spi_Vel_Zref, using_spi, n_spi, n_inj, &
                spi_Vel_RxZref, spi_quantity, spi_abl_model,        &
                ng_radius_ratio, ng_radius_min, spi_angle,          &
                spi_L_inj, spi_L_inj_diff,                          &
                K_Dmv, A_Dmv, L_tube, V_Dmv, P_Dmv,                 &
                spi_Vel_diff, t_ns, JET_MGI, ASDEX_MGI,             &
                delta_n_convection, nimp_bg, output_prad_phi,       &
                RMP_on, RMP_har_cos,RMP_har_sin, spi_shard_file,    &
                spi_plume_file, spi_plume_hdf5,                     &
                RMP_growth_rate, RMP_ramp_up_time,                  &
                RMP_psi_cos_file, RMP_psi_sin_file,                 &
                Number_RMP_harmonics,RMP_har_cos_spectrum,          &
                RMP_har_sin_spectrum, imp_type, adas_dir, n_adas,   &
                index_main_imp,                                     &
                amix, amix_freeb, equil_accuracy, use_imp_adas,     &
                equil_accuracy_freeb, current_ref, FB_Ip_position,  &
                FB_Ip_integral, Z_axis_ref, FB_Zaxis_position,      &
                FB_Zaxis_derivative,FB_Zaxis_integral, start_VFB,   &
                n_feedback_current, n_feedback_vertical,            &
                n_iter_freeb, n_pf_coils, pf_coils, R_axis_ref,     &
                axis_srch_radius,                                   &
                starwall_equil_coils, freeb_equil_iterate_area,     &
                psi_offset_freeb, diag_coils, rmp_coils,            &
                voltage_coils, vert_FB_amp, find_pf_coil_currents,  &
                delta_psi_GS, newton_GS_fixbnd, newton_GS_freebnd,  &
                pastix_maxthrd, eta_ohmic, centralize_harm_mat,     &
                vert_FB_amp_ts, vert_FB_gain, vert_pos_file,        & 
                vert_FB_tact, start_VFB_ts, I_coils_max, rad_FB_amp,&
                autodistribute_modes, modes_per_family,             &
                mode_families_modes, n_mode_families,               &
                weights_per_family, autodistribute_ranks,           &
                ranks_per_family,                                   &
                tgnum_psi, tgnum_u, tgnum_zj, tgnum_w, tgnum_rho,   &
                tgnum_T, tgnum_Ti, tgnum_Te, tgnum_vpar, tgnum_rhon,&
                tgnum_nre, tgnum_AR, tgnum_AZ, tgnum_A3,            &
                n_particles, tstep_particles, nstep_particles,      &
                nsubstep_particles,                                 &
                filter_perp,    filter_hyper,    filter_par,        &
                filter_perp_n0, filter_hyper_n0, filter_par_n0,     &
                use_cx, use_sputtering, use_ionisation,             &
                use_ncs, use_pcs, use_ccs,                          &
                min_sheath_angle, bcs,                              &
                use_sc, add_sources_in_sc, visco_sc_num,            &
                D_perp_sc_num, D_par_sc_num, ZK_perp_sc_num,        &
                ZK_par_sc_num, ZK_i_perp_sc_num, ZK_i_par_sc_num,   &
                ZK_e_perp_sc_num, ZK_e_par_sc_num, visco_par_sc_num,&
                Dn_pol_sc_num, Dn_p_sc_num, cte_current_FB_fact
      
if (my_id .eq. 0) then

  ! --- Preset input parameters to reasonable default values.
  call preset_parameters()

  call vacuum_preset(my_id, freeboundary_equil, freeboundary, resistive_wall)
  
  ! --- Model-specific presets
  particlesource_psin = 100.d0
  
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
  endif

  ! --- Calculate normalisation factor for MGI source (related to its toroidal shape)
  ns_tor_norm = ns_deltaphi * PI**0.5 * ERF(PI/ns_deltaphi)

  if (trim(R_Z_psi_bnd_file) .ne. 'none') then
    ! --- Open the file.
    OPEN(UNIT=243, FILE=R_Z_psi_bnd_file, FORM='FORMATTED', STATUS='OLD', ACTION='READ', IOSTAT=err)
    if ( err /= 0 ) then
      write(*,*) 'ERROR in initialise_parameters: Cannot open file '//TRIM(R_Z_psi_bnd_file)//'.'
      stop
    endif
    write(*,'(A)') ' boundary info from R_Z_psi_bnd_file: R_boundary, Z_boundary, psi_boundary '
    do i=1,n_boundary
      read(243,*) R_boundary(i),Z_boundary(i),psi_boundary(i)
      write(*,*) R_boundary(i),Z_boundary(i),psi_boundary(i)
    enddo
  endif

  ! --- Calculate JOREK gamma_sheath from gamma_stangeby if provided (otherwise the other way around)
  if ( with_TiTe ) then
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
  else
    if (gamma_stangeby > -1.d89) then
      gamma_sheath = (gamma-1.d0) * (0.5d0*gamma_stangeby - 1.d0 - 0.5d0*gamma)
    else
      gamma_stangeby = 2.d0 * ( gamma_sheath / (gamma-1.d0) + 1.d0 + 0.5d0 * gamma )
    end if
  end if

  if (sum(nstep_n) .gt. 0) then
    nstep = sum(nstep_n)
  else
    tstep_n    = 0.d0
    tstep_n(1) = tstep
    nstep_n    = 0
    nstep_n(1) = nstep
  endif

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

  if (using_spi) call init_spi_all()

end if

return
end subroutine initialise_parameters
