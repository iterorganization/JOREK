!> Initialise input parameters and read the input namelist
subroutine initialise_parameters(my_id, filename)

use tr_module
use phys_module
use mumps_module,  only: no_zeros_mumps, mumps_ordering
use pastix_module, only: no_zeros_pastix, pastix_smp_only, &
    pastix_maxthrd
use vacuum
use live_data

implicit none

! --- Routine parameters
integer,                      intent(in) :: my_id
character(len=*),             intent(in) :: filename

! --- Local variables
integer :: ierr, err, i

! --- Namelist with input parameters.
namelist /in1/  tstep, nstep, tstep_n, nstep_n,                     &
                rst_hdf5, rst_hdf5_version, keep_current_prof,      &
                eta, visco, restart, regrid, write_ps,              &
                force_horizontal_Xline,                             &
                n_R, n_Z, n_radial, n_pol, n_tht, n_flux,           &
                n_open, n_private, n_leg, n_ext,                    &
                n_outer, n_inner, n_up_priv, n_up_leg,              &
                psi_axis_init, XR_r, SIG_r, XR_tht, SIG_tht,        &
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
                tokamak_device, manipulate_psi_map,                 &
                F0, time_evol_scheme,                               &
				T_min_neg,rho_min_neg,                              &
                zjz_0, zjz_1, zj_coef,                              &
                rho_0, rho_1, rho_coef,                             &
                T_0,   T_1,   T_coef,                               &
                FF_0,  FF_1,  FF_coef,                              &
                ZK_par, ZK_perp, D_par, D_perp,                     &
                particlesource, heatsource,                         &
                eta_num, visco_num, tgnum,                          &
                ellip, tria_u, tria_l, quad_u, quad_l,              &
                xampl, xwidth, xsig, xtheta, xshift, xleft, xpoint, &
                xcase, D_perp_file, ZK_perp_file,                   &
                rho_file, T_file, ffprime_file,                     &
                freeboundary, resistive_wall, freeb_change_indices, &
                wall_resistivity, wall_resistivity_fact,            &
                use_mumps_eq, use_pastix_eq, use_strumpack_eq,      &
                use_mumps, mumps_ordering,                          &
                use_BLR_compression, epsilon_BLR, just_in_time_BLR, &
                use_pastix, use_murge, use_murge_element, use_wsmp, &
                n_tor_fft_thresh,                                   &
                pastix_smp_only, refinement, force_central_node,    &
                fix_axis_nodes,                                     &
                grid_to_wall, use_strumpack,                        &
                adaptive_time, equil, bench_without_plot,           &
                no_zeros_pastix, no_zeros_mumps,                    &
                eta_T_dependent, visco_T_dependent,                 &
                heatsource_psin, heatsource_sig,                    &
                particlesource_psin, particlesource_sig,            &
                edgeparticlesource, edgeparticlesource_psin,        &
                edgeparticlesource_sig,                             &
                produce_live_data, gmres, gmres_max_iter,           &
                wall_file,                                          &
                n_limiter, R_limiter, Z_limiter,                    &
                first_target_point, last_target_point,		    &
                keep_n0_const, linear_run, export_for_nemec,        &
                amix, amix_freeb, equil_accuracy,                   &
                equil_accuracy_freeb, current_ref, FB_Ip_position,  &
                FB_Ip_integral, Z_axis_ref, FB_Zaxis_position,      &
                FB_Zaxis_derivative,FB_Zaxis_integral, start_VFB,   &
                n_feedback_current, n_feedback_vertical,            &
                n_iter_freeb, n_pf_coils, pf_coils, R_axis_ref,     &
                axis_srch_radius,                                   &
                starwall_equil_coils, freeb_equil_iterate_area,     &
                psi_offset_freeb, diag_coils, rmp_coils,            &
                voltage_coils, vert_FB_amp, find_pf_coil_currents,  &
                pastix_maxthrd, centralize_harm_mat,                & 
                vert_FB_amp_ts, vert_FB_gain, vert_pos_file,        & 
                delta_psi_GS, newton_GS_fixbnd, newton_GS_freebnd,  &
                vert_FB_tact, start_VFB_ts, I_coils_max, rad_FB_amp,&
                autodistribute_modes, modes_per_family,             &
                mode_families_modes, n_mode_families,               &
                weights_per_family, autodistribute_ranks,           &
                ranks_per_family, cte_current_FB_fact


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
  !=========================================

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

return
end subroutine initialise_parameters
