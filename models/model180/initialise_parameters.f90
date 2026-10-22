!> Initialise input parameters and read the input namelist
subroutine initialise_parameters(my_id, filename)

use tr_module
use phys_module
use pellet_module
use vacuum
use live_data

implicit none

! --- Routine parameters
integer,                      intent(in) :: my_id
character(len=*),             intent(in) :: filename

! --- Local variables
integer :: ierr,err,i

! --- Namelist with input parameters.
namelist /in1/  tstep, nstep, tstep_n, nstep_n,                     &
                rst_hdf5, rst_hdf5_version, keep_current_prof,      &
                init_current_prof, eta, visco, visco_par,           &
                restart, rst_format, regrid, bootstrap, write_ps,   &
                bootstrap_psin_cutoff,                              &
                regrid_from_rz, force_horizontal_Xline,             &
                n_R, n_Z, n_radial, n_pol, n_tht, n_flux,           &
                n_open, n_private, n_leg, n_ext, i_plane_rtree,     &
                n_outer, n_inner, n_up_priv, n_up_leg, m_pol_bc,    &
                SDN_threshold,                                      &
                psi_axis_init, XR_r, SIG_r, XR_tht, SIG_tht,        &
                xr_closed,                                          &
                SIG_closed, SIG_open, SIG_private, SIG_theta,       &
                SIG_leg_0, SIG_leg_1, dPSI_open, dPSI_private,      &
                SIG_up_leg_0, SIG_up_leg_1, SIG_up_priv,            &
                SIG_outer, SIG_inner,                               &
                dPSI_outer, dPSI_inner, dPSI_up_priv,               &
                nout, xr1, sig1, xr2, sig2,                         &
                R_begin, R_end, Z_begin, Z_end,                     &
                R_geo, Z_geo, amin, mf, fbnd, fpsi, mode,           &
                R_Z_psi_bnd_file,                                   &
                R_boundary, Z_boundary, psi_boundary, n_boundary,   &
                n_pfc, n_tor_fft_thresh, manipulate_psi_map,        &
                Rmin_pfc, Rmax_pfc, Zmin_pfc, Zmax_pfc, current_pfc,&
                tokamak_device, gvec_grid_import, extended_boundary,&
                j_cutoff_rcoord, j_cutoff_sig, bloating_factor,     &
                F0, gamma_sheath, density_reflection,               &
                zjz_0, zjz_1, zj_coef,                              &
                rho_0, rho_1, rho_coef,                             &
                T_0,   T_1,   T_coef,                               &
                Ti_0, Ti_1, Ti_coef,                                &
                Te_0, Te_1, Te_coef, TiTe_ratio,                    &
                FF_0,  FF_1,  FF_coef,                              &
                Phi_0,  Phi_1,  Phi_coef,                           &
                ZK_par, ZK_par_max, ZK_perp, D_par, D_perp,         &
                particlesource, tauIC,                              &
                heatsource, heatsource_i, heatsource_e,             &
                eta_num, visco_num, visco_par_num, D_perp_num,      &
                ZK_perp_num,                                        &
                pellet_amplitude, pellet_R, pellet_Z, pellet_phi,   &
                pellet_radius, pellet_sig, pellet_length,           &
                pellet_psi, pellet_delta_psi, pellet_density,       &
                pellet_velocity_R, pellet_velocity_Z, pellet_theta, &
                pellet_ellipse,                                     &
                central_density, central_mass,                      &
                pellet_particles, use_pellet,                       &
                ellip,tria_u,tria_l,quad_u,quad_l,                  &
                xampl,xwidth,xsig,xtheta,xshift,xleft, xpoint,      &
                xcase, D_perp_file, ZK_perp_file,                   &
                rho_file, T_file, ffprime_file,                     &
                rot_file, domm_file, phi_file,                      &
                normalized_velocity_profile,                        &
                freeboundary_equil, freeboundary,  freeb_change_indices, &
                resistive_wall,                                     &
                wall_resistivity, wall_resistivity_fact,            &
                bc_natural_open,                                    &
                use_mumps_eq, use_pastix_eq, use_strumpack_eq,      &
                use_mumps_prj, use_pastix_prj, use_strumpack_prj,   &
                use_mumps, mumps_ordering,                          &
                use_BLR_compression, epsilon_BLR, just_in_time_BLR, &
                use_pastix, use_murge, use_murge_element, use_wsmp, &
                refinement, force_central_node,    &
                fix_axis_nodes, treat_axis,                         &
                grid_to_wall, use_strumpack,                        &
                adaptive_time, equil, bench_without_plot,           &
                eta_T_dependent, visco_T_dependent,                 &
                zkpar_T_dependent,                                  & 
                heatsource_psin, heatsource_i_psin,                 &
                heatsource_e_psin, heatsource_sig,                  &
                heatsource_i_sig, heatsource_e_sig,                 &
                particlesource_psin, particlesource_sig,            &
                edgeparticlesource, edgeparticlesource_psin,        &
                edgeparticlesource_sig,                             &
                particlesource_gauss, heatsource_gauss,             &
                heatsource_gauss_i, heatsource_gauss_e,             &
                heatsource_gauss_psin, heatsource_gauss_i_psin,     &
                heatsource_gauss_e_psin, heatsource_gauss_sig,      &
                heatsource_gauss_i_sig, heatsource_gauss_e_sig,     &
                particlesource_gauss_psin, particlesource_gauss_sig,&
                produce_live_data, gmres, gmres_max_iter,           &
                gmres_m, gmres_4, gmres_tol, iter_precon,           &
                tgnum,  pastix_pivot, max_steps_noUpdate,           &
                export_for_nemec,                                   &
                RMP_on, RMP_har_cos,RMP_har_sin,                    &
                RMP_growth_rate, RMP_ramp_up_time,                  &
                RMP_psi_cos_file, RMP_psi_sin_file,                 &
                V_0,V_1,V_coef, output_bnd_elements,                &
                wall_file,                                          &
                n_limiter, R_limiter, Z_limiter,                    &
                first_target_point, last_target_point,              &
                NEO, neo_file, aki_neo_const, amu_neo_const,        &
                time_evol_scheme, corr_neg_temp_coef,               &
                corr_neg_dens_coef, D_prof_neg, ZK_prof_neg,        &
                D_prof_neg_thresh, ZK_prof_neg_thresh, T_min,       &
                amix, amix_freeb, equil_accuracy,                   &
                equil_accuracy_freeb, current_ref, FB_Ip_position,  &
                FB_Ip_integral, Z_axis_ref, FB_Zaxis_position,      &
                FB_Zaxis_derivative,FB_Zaxis_integral, start_VFB,   &
                n_feedback_current, n_feedback_vertical,            &
                n_iter_freeb, n_pf_coils, pf_coils,                 &
                axis_srch_radius,                                   &
                starwall_equil_coils, freeb_equil_iterate_area,     &
                psi_offset_freeb, diag_coils, rmp_coils,            &
                voltage_coils, vert_FB_amp, find_pf_coil_currents,  &
                pastix_maxthrd, eta_ohmic, centralize_harm_mat,     &
                vert_FB_amp_ts, vert_FB_gain, axis_pos_file,        & 
                rad_FB_amp_ts, rad_FB_gain,                         &
                pos_FB_tact,  start_PFB_ts, I_coils_max,            &
                autodistribute_modes, modes_per_family,             &
                mode_families_modes, n_mode_families,               &
                weights_per_family, autodistribute_ranks,           &
                ranks_per_family,                                   &
                use_manual_random_seed, manual_seed,                &
                use_fixed_rng_value, fixed_rng_value                                

                
namelist /dommcoef/  R_domm, dcoef

if (my_id .eq. 0) then

  ! --- Preset input parameters to reasonable default values.
  call preset_parameters()
  call vacuum_preset(my_id, freeboundary_equil, freeboundary, resistive_wall)
  
  ! --- Model-specific presets
  particlesource_psin = 100.d0
  
  ! Set false because temperature profiles are not valid in stellarator models
  write_ps = .false.

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


  if ( ( n_tor .eq. 1 ) .and. freeboundary .and. (.not. freeboundary_equil) ) then
    write(*,*) 'WARNING: The parameter freeboundary is automatically changed to .false. since n_tor==1 and freeboundary_equil is .false.'
    freeboundary= .false.
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
 !=========================================
  
  if (sum(nstep_n) .gt. 0) then
    nstep = sum(nstep_n)
    tstep = tstep_n(1)
  else
    tstep_n    = 0.d0
    tstep_n(1) = tstep
    nstep_n    = 0
    nstep_n(1) = nstep
  endif

  if (i_plane_rtree .gt. n_plane .or. i_plane_rtree .lt. 1) then
    write(*,*) 'ERROR: The variable i_plane_rtree must be between 1 and the total number of poloidal planes'
    write(*,'(A,I4,A,I4)') 'i_plane_rtree = ', i_plane_rtree, '; n_plane = ', n_plane
    stop
  end if

  call allocate_live_data()

endif

! --- Read numerical profiles for rho, T, ff', toroidal rotation and neoclassical coefficients.
call read_num_profiles(my_id)

! --- Determine the derivatives of the numerical input profiles.
call derive_num_profiles(my_id)

domm = ( domm_file /= 'none' )
if (domm .and. my_id .eq. 0 ) then
  open(43, file=domm_file, status='old', action='read', iostat=ierr)
  if (ierr /= 0) then
    write(*,*) 'ERROR: COULD NOT OPEN FILE "', trim(domm_file), '".'
    stop
  end if
  read(43,dommcoef)
  close(43)
  if (R_domm .le. 0.d0) then
    write(*,*) 'ERROR: THE VARIABLE R_DOMM IS NOT SPECIFIED IN THE FILE "', trim(domm_file), '".'
    write(*,*) 'A POSITIVE VALUE FOR R_DOMM MUST BE SPECIFIED'
    stop
  end if
end if

#ifdef USE_DOMM
! Runtime check: If compiled with USE_DOMM, vacuum field representation uses ONLY Dommaschk potentials (no FE correction)
! This requires domm_file to provide dcoef array. Without it, dcoef=0 causes NaN in field calculations.
if (.not. domm .and. my_id .eq. 0) then
  write(*,*) '**************************************************************************'
  write(*,*) 'WARNING: Compiled with USE_DOMM=1 (Dommaschk-only vacuum field)'
  write(*,*) '         but domm_file="', trim(domm_file), '"'
  write(*,*) '         Vacuum field will be ZERO (dcoef array uninitialized)!'
  write(*,*) '         This will cause NaN in field line tracing and Poincare plots.'
  write(*,*) ''
  write(*,*) 'SOLUTION: Either:'
  write(*,*) '  1. Provide domm_file with Dommaschk coefficients in namelist'
  write(*,*) '  2. Recompile with USE_DOMM=0 to use GVEC import + FE correction'
  write(*,*) '**************************************************************************'
  stop
end if
#endif
  
! --- TiTe_ratio has to be between 1 and 0
if ((with_TiTe) .and. (my_id .eq. 0)) then
  if ((TiTe_ratio .ge. 1.0) .or. (TiTe_ratio<0.d0)) then
    write(*,*) 'ERROR: The temperature ratio coefficient should be 0<=TiTe_ratio<1.0 but is ', TiTe_ratio
    stop
  end if
end if

return
end subroutine initialise_parameters
