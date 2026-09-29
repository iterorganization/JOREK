!> Presets input parameters to reasonable default values.
!!
!! The model-specific routines initialise_parameters may overwrite
!! these defaults according to the requirements of the respective
!! model.
subroutine preset_parameters
  
  use phys_module
  implicit none

  integer :: i, j ! for iterations
  character(len=20) :: tmp !< for temporary name tag

  time_evol_scheme = 'Crank-Nicholson'
  
  n_tor_fft_thresh = 2
  if(jorek_model == 305 .or. jorek_model == 306) n_tor_fft_thresh = 99
  
  ! --- DoubleNull flag
  xcase = LOWER_XPOINT
  
  tstep    = 1.d0
  tstep_n  = 1.d0
  nstep    = 0
  nstep_n  = 0
  
  eta_T_dependent   = .true.
  eta_coul_log_dep  = .true.
  visco_T_dependent = .true.
  ZKpar_T_dependent = .true.

  eta_num_T_dependent   = .false.
  visco_num_T_dependent = .false.
  add_sources_in_sc     = .false.

  eta_num_psin_dependent= .false.
  eta_num_prof = 0.d0
  eta_num_prof(1) = 0.8d0
  eta_num_prof(2) = 0.03d0

  eta           = 1.d-5
  T_max_eta     = 1.d99
  eta_ohmic     = 0.d0
  T_max_eta_ohm = 1.d99
  
  TiTe_ratio    = 0.5d0

  visco = 1.d-5
  T_max_visco   = 1.d99
  visco_par = 1.d-5
  visco_par_par = 0.d0  
  visco_heating     = 0.d0
  visco_par_heating = 0.d0
  visco_old_setup   = .false.
  
  central_density = 1.d0            ! the central density in units 10^20 m^-3
  central_mass    = 2.01410177811d0 ! the central average mass (atomic mass of deuterium, including electron)

  n_tor_restart= 0
  restart      = .false.
  import_equil = .false.
  regrid       = .false.
  regrid_from_rz = .false.
  rst_format   = 0             ! use 'old' format for restart import
  write_ps     = .true.           ! write postscript file at the end of the run 
  gvec_grid_import = .false.
  extended_boundary = .false.
  j_cutoff_rcoord = 99.0
  j_cutoff_sig = 0.025
  bloating_factor = 1.0d0

  freeboundary_equil = .false. ! use free or fixed boundary equilibrium
  freeboundary       = .false. ! use free or fixed boundary?
  resistive_wall     = .false. ! use a resistive or ideal wall?    (freeboundary only)
  freeb_equil_iterate_area = .false.
  freeb_change_indices = .true. ! exchange grid node indices to parallelize boundary integral

  bc_natural_flux    = .false.! boundary conditions for flux surface boundaries (2 and 3)
  bc_natural_open    = .false. ! use sheath (Bohm) boundary conditions

  gamma_sheath       = 4.5d0  ! sheath transmission factor (single fluid) in the JOREK definition
  gamma_stangeby     = -1.d99 ! sheath transmission factor (single fluid) given by Stangeby
  gamma_sheath_e     = 3.00d0 ! sheath transmission factor (electron fluid) in the JOREK definition
  gamma_e_stangeby   = -1.d99 ! sheath transmission factor (electron fluid) given by Stangeby
  gamma_sheath_i     = -1.11d-1! sheath transmission factor (ion fluid) in the JOREK definition
  gamma_i_stangeby   = -1.d99 ! sheath transmission factor (ion fluid) given by Stangeby
  density_reflection = 0.d0   ! reflection coefficient for outgoing density
  neutral_reflection = 0.d0   ! reflection coefficient for (fluid) neutrals
  imp_reflection     = 0.d0   ! reflection coefficient for (fluid) impurities
  
  deuterium_adas        = .false. 
  deuterium_adas_1e20   = .false. 
  old_deuterium_atomic  = .false. 
  mach_one_bnd_integral = .false. ! implement Mach one condition as boundary integral
  Vpar_smoothing        = .false. ! smooth the transitions of Vpar positive/negavtive at B.n
  Vpar_smoothing_coef   = (/0.01d0, 0.d0, 0.d0 /) !(/ 0.01d0, 0.016d0, 0.00575446347d0/)
  min_sheath_angle      = 1.d0   ! 1 degree (not in radians)

  amix                 = 0.d0
  amix_freeb           = 0.85d0
  equil_accuracy       = 1.d-6
  equil_accuracy_freeb = 1.d-6
  axis_srch_radius     = 99.d0
  delta_psi_GS         = 10000.d0
  newton_GS_fixbnd     = .false.
  newton_GS_freebnd    = .true.
  
  n_R          = 0
  n_Z          = 0

  n_radial     = 11
  n_pol        = 16

  n_flux       = 11
  n_tht        = 16
  n_tht_equidistant = .false.
  m_pol_bc     = 1
  i_plane_rtree = 1
  
  n_open       = 5
  n_outer      = 0
  n_inner      = 0
  n_leg        = 5
  n_leg_out    = 0
  n_private    = 5
  n_up_leg     = 0
  n_up_leg_out = 0
  n_up_priv    = 0
  
  n_ext        = 0

  export_polar_boundary = .false.

  psi_axis_init = -0.1d0
  XR_r(:)       = 999.d0
  SIG_r(:)      = 999.d0
  XR_tht(:)     = 999.d0
  SIG_tht(:)    = 999.d0
  XR_z(:)       = 999.d0
  SIG_z(:)      = 999.d0
  bgf_r         = 0.7
  bgf_z         = 0.7
  bgf_rpolar    = 0.6
  bgf_tht       = 0.6

  xr_closed   = (/   1.0d0, 9999.d0, 9999.d0 /)
  SIG_closed  = (/   0.1d0, 9999.d0, 0.1d0   /)
  SIG_open    = 0.1d0
  SIG_outer   = 0.1d0
  SIG_inner   = 0.1d0
  SIG_private = 0.1d0
  SIG_up_priv = 0.1d0
  SIG_theta   = 0.03d0
  SIG_theta_up= 999.d0
  SIG_leg_0   = 0.05d0
  SIG_leg_1   = 0.2d0
  SIG_up_leg_0= 0.05d0
  SIG_up_leg_1= 0.2d0
  
  dPSI_open    = 0.11
  dPSI_outer   = 0.11
  dPSI_inner   = 0.11
  dPSI_private = 0.03
  dPSI_up_priv = 0.03

  forceSDN = .false.
  SDN_threshold = 1.d-4
  
  R_geo     = 10.d0
  Z_geo     = 0.d0
  amin      = 1.d0

  R_domm        = -10.d0

  F0        = 10.d0
  GAMMA     = 5.d0 / 3.d0

  mf        = 2
  fbnd      = 0.d0;   fbnd(1)  = 2.d0

  R_boundary   = 0.d0
  Z_boundary   = 0.d0
  psi_boundary = 0.d0
  n_boundary   = 0

  n_pfc       = 0
  Rmin_pfc    = 0.d0
  Rmax_pfc    = 0.d0
  Zmin_pfc    = 0.d0
  Zmax_pfc    = 0.d0
  current_pfc = 0.d0

  n_jropes       = 0
  R_jropes       = 0.d0
  Z_jropes       = 0.d0
  w_jropes       = 0.d0
  current_jropes = 0.d0
  rho_jropes     = 0.d0
  T_jropes       = 0.d0

  bootstrap = .false.
  bootstrap_psin_cutoff = 0.9995

  ellip  = 1.d0
  tria_u = 0.d0
  tria_l = 0.d0
  quad_l = 0.d0
  quad_u = 0.d0

  xampl  = 0.d0
  xwidth = 0.d0
  xsig   = 1.d0
  xtheta = 0.d0
  xshift = 0.d0
  xleft  = 0.d0
  xpoint = .false.
  force_horizontal_Xline = .false.
  Z_xpoint_limit(1) = -0.4d0
  Z_xpoint_limit(2) =  0.4d0
  xpoint_search_tries = 500

  xr1  = 9999.d0
  sig1 = 9999.d0
  xr2  = 99999.d0
  sig2 = 99999.d0

  R_begin = -0.1d0
  R_end   =  0.1d0
  Z_begin = -0.1d0
  Z_end   = 0.1d0

  rect_grid_vac_psi = 0.d0
  
  maintain_profiles = .false.
  ZK_perp(1:5)   = (/ 1.d-5, 0.d0, 0.d0, 99.d0, 99.d0 /)
  ZK_i_perp(1:5) = (/ 1.d-5, 0.d0, 0.d0, 99.d0, 99.d0 /)
  ZK_e_perp(1:5) = (/ 1.d-5, 0.d0, 0.d0, 99.d0, 99.d0 /)
  ZK_par       = 1.d0
  ZK_i_par     = 1.d0
  ZK_e_par     = 1.d0
  ZK_par_max   = 1.d20
  D_perp(1:5)  = (/ 1.d-5, 0.d0, 0.d0, 99.d0, 99.d0 /)
  D_par        = 0.d0
  V_pinch_gauss = 0.d0
  V_pinch_psin  = 0.d0
  V_pinch_sig   = 1.d0
  D_perp_imp(1:5)  = (/ 1.d-5, 0.d0, 0.d0, 99.d0, 99.d0 /)
  D_par_imp        = 0.d0

  D_prof_neg         = 1.d-5
  D_prof_neg_thresh  = 0.d0 ! default is zero for keeping the old behavior
  D_prof_imp_neg_thresh  = -1.d3 ! disabled by default to avoid convergence issues
  D_prof_tot_neg_thresh  = 0.d0 ! default is zero for keeping the old behavior

  D_imp_extra_R = 0.d0
  D_imp_extra_Z = 0.d0
  D_imp_extra_p = 0.d0
  D_imp_extra_neg      = 1.d-6
  D_imp_extra_neg_thresh  = -1.d3 ! to disable by default enhanced impurities diffusion in model501
  ZK_prof_neg        = 1.d-5
  ZK_par_neg         = 1.d-3
  ZK_prof_neg_thresh = 0.d0 ! default is zero for keeping the old behavior
  ZK_par_neg_thresh  = 0.d0
  ZK_e_prof_neg        = 1.d-5
  ZK_e_par_neg         = 1.d-3
  ZK_e_prof_neg_thresh = 0.d0 ! default is zero for keeping the old behavior
  ZK_e_par_neg_thresh  = 0.d0
  ZK_i_prof_neg        = 1.d-5
  ZK_i_par_neg         = 1.d-3
  ZK_i_prof_neg_thresh = 0.d0 ! default is zero for keeping the old behavior
  ZK_i_par_neg_thresh  = 0.d0

  HW_coef    = 0.d0
  HW_coef(1) = 1.d-6
  HW_coef(2) = 1.d0

  ne_SI_min          = 1.d18
  Te_eV_min          = 5.
  rn0_min            = 1.d-8
  T_min              = 1.0d-20  !-1.0d20
  rho_min            = 1.0d-20  !-1.0d20
  T_min_neg          = -1.d12 !< only used if T_min_neg>0 , 2.01d-5*central_density*Tmin_ev (cd = 1, 20 eV)
  T_min_ZKpar        = -1.d12 
  Ti_min_ZKpar       = -1.d12 
  Te_min_ZKpar       = -1.d12 
  rho_min_neg        = -1.d12
  
  implicit_heat_source = 0.d0
  
  corr_neg_temp_coef(:) = (/ 0.5, 0.5 /)
  corr_neg_dens_coef(:) = (/ 0.5, 0.5 /)

  eta_num            = 0.d0
  visco_num          = 0.d0
  visco_par_num      = 0.d0
  D_perp_num         = 0.d0
  D_perp_num_tanh    = 0.d0; D_perp_num_tanh_psin    = 3.d-1; D_perp_num_tanh_sig    = 1.d-1
  ZK_perp_num        = 0.d0
  ZK_perp_num_tanh   = 0.d0; ZK_perp_num_tanh_psin   = 3.d-1; ZK_perp_num_tanh_sig   = 1.d-1
  ZK_i_perp_num      = 0.d0
  ZK_i_perp_num_tanh = 0.d0; ZK_i_perp_num_tanh_psin = 3.d-1; ZK_i_perp_num_tanh_sig = 1.d-1
  ZK_e_perp_num      = 0.d0
  ZK_e_perp_num_tanh = 0.d0; ZK_e_perp_num_tanh_psin = 3.d-1; ZK_e_perp_num_tanh_sig = 1.d-1
  Dn_perp_num        = 0.d0

  use_sc = .false.
  visco_sc_num     = 0.d0
  D_perp_sc_num    = 0.d0
  D_par_sc_num     = 0.d0
  ZK_perp_sc_num   = 0.d0
  ZK_par_sc_num    = 0.d0
  ZK_i_perp_sc_num = 0.d0
  ZK_i_par_sc_num  = 0.d0
  ZK_e_perp_sc_num = 0.d0
  ZK_e_par_sc_num  = 0.d0
  visco_par_sc_num = 0.d0
  Dn_pol_sc_num    = 0.d0
  Dn_p_sc_num      = 0.d0
  D_perp_imp_sc_num= 0.d0
  D_par_imp_sc_num = 0.d0
 
  use_vms = .false.
  vms_coeff_AR     = 0.d0
  vms_coeff_AZ     = 0.d0
  vms_coeff_A3     = 0.d0
  vms_coeff_UR     = 0.d0
  vms_coeff_UZ     = 0.d0
  vms_coeff_Up     = 0.d0
  vms_coeff_rho    = 0.d0 
  vms_coeff_T      = 0.d0
  vms_coeff_Te     = 0.d0
  vms_coeff_Ti     = 0.d0  
  vms_coeff_rhon   = 0.d0 
  vms_coeff_rhoimp = 0.d0

  heatsource          = 1.e-7
  heatsource_e        = 0.5e-7
  heatsource_i        = 0.5e-7
  heatsource_psin     = 1.0d0
  heatsource_e_psin   = 1.0d0
  heatsource_i_psin   = 1.0d0
  heatsource_sig      = 0.1d0
  heatsource_e_sig    = 0.1d0
  heatsource_i_sig    = 0.1d0
  particlesource      = 1.e-5
  particlesource_psin = 1.0d0
  particlesource_sig  = 0.1d0
  edgeparticlesource      = 0.d0
  edgeparticlesource_psin = 0.98
  edgeparticlesource_sig  = 0.01
  heatsource_gauss          = 0.d0
  heatsource_gauss_e        = 0.d0
  heatsource_gauss_i        = 0.d0
  heatsource_gauss_psin     = 0.9d0
  heatsource_gauss_e_psin   = 0.9d0
  heatsource_gauss_i_psin   = 0.9d0
  heatsource_gauss_sig      = 0.1d0
  heatsource_gauss_e_sig    = 0.1d0
  heatsource_gauss_i_sig    = 0.1d0
  particlesource_gauss      = 0.d0
  particlesource_gauss_psin = 0.9d0
  particlesource_gauss_sig  = 0.1d0
  neutral_line_source       = 0.d0
  neutral_line_R_start      = 1.d20
  neutral_line_Z_start      = 1.d20
  neutral_line_R_end        = 2.d20
  neutral_line_Z_end        = 2.d20

  ! ------------------------------------------
  ! --- Default boundary conditions ----------
  ! ------------------------------------------
  loop_voltage = 0.d0
  ! --- Dirichlet
  bcs(:)%dirichlet%psi     = .true.
  bcs(:)%dirichlet%u       = .true.
  bcs(:)%dirichlet%zj      = .true.
  bcs(:)%dirichlet%w       = .true.
  bcs(:)%dirichlet%rho     = .true.
  bcs(:)%dirichlet%T       = .true.
  bcs(:)%dirichlet%Ti      = .true.
  bcs(:)%dirichlet%Te      = .true.
  bcs(:)%dirichlet%Vpar    = .true.
  bcs(:)%dirichlet%rhon    = .true.
  bcs(:)%dirichlet%rho_imp = .true.
  bcs(:)%dirichlet%nre     = .true.
  bcs(:)%dirichlet%AR      = .true.
  bcs(:)%dirichlet%AZ      = .true.
  bcs(:)%dirichlet%A3      = .true.

  bcs(  1)%dirichlet%rho   = .false.
  bcs(4:5)%dirichlet%rho   = .false.
  bcs(  9)%dirichlet%rho   = .false.
  bcs( 11)%dirichlet%rho   = .false.
  bcs( 15)%dirichlet%rho   = .false.
  bcs( 19)%dirichlet%rho   = .false.

  bcs(  1)%dirichlet%T     = .false.
  bcs(4:5)%dirichlet%T     = .false.
  bcs(  9)%dirichlet%T     = .false.
  bcs( 11)%dirichlet%T     = .false.
  bcs( 15)%dirichlet%T     = .false.
  bcs( 19)%dirichlet%T     = .false.

  bcs(  1)%dirichlet%Te    = .false.
  bcs(4:5)%dirichlet%Te    = .false.
  bcs(  9)%dirichlet%Te    = .false.
  bcs( 11)%dirichlet%Te    = .false.
  bcs( 15)%dirichlet%Te    = .false.
  bcs( 19)%dirichlet%Te    = .false.

  bcs(  1)%dirichlet%Ti    = .false.
  bcs(4:5)%dirichlet%Ti    = .false.
  bcs(  9)%dirichlet%Ti    = .false.
  bcs( 11)%dirichlet%Ti    = .false.
  bcs( 15)%dirichlet%Ti    = .false.
  bcs( 19)%dirichlet%Ti    = .false.

  bcs(  1)%dirichlet%vpar  = .false.
  bcs(4:5)%dirichlet%vpar  = .false.
  bcs(  9)%dirichlet%vpar  = .false.
  bcs( 11)%dirichlet%vpar  = .false.
  bcs( 15)%dirichlet%vpar  = .false.
  bcs( 19)%dirichlet%vpar  = .false.

  bcs(  1)%dirichlet%rhon  = .false.
  bcs(4:5)%dirichlet%rhon  = .false.
  bcs(  9)%dirichlet%rhon  = .false.
  bcs( 11)%dirichlet%rhon  = .false.
  bcs( 15)%dirichlet%rhon  = .false.
  bcs( 19)%dirichlet%rhon  = .false.

  bcs(  1)%dirichlet%rho_imp  = .false.
  bcs(4:5)%dirichlet%rho_imp  = .false.
  bcs(  9)%dirichlet%rho_imp  = .false.
  bcs( 11)%dirichlet%rho_imp  = .false.
  bcs( 15)%dirichlet%rho_imp  = .false.
  bcs( 19)%dirichlet%rho_imp  = .false.

  ! --- Mach 1
  bcs(:)%mach1   = .false.

  bcs(  1)%mach1 = .true.
  bcs(3:5)%mach1 = .true.
  bcs(  9)%mach1 = .true.
  bcs( 11)%mach1 = .true.
  bcs( 15)%mach1 = .true.
  bcs( 19)%mach1 = .true.

  ! --- Natural BCs
  bcs(:)%natural%rho     = .false.
  bcs(:)%natural%T       = .false.
  bcs(:)%natural%Ti      = .false.
  bcs(:)%natural%Te      = .false.
  bcs(:)%natural%Vpar    = .false.
  bcs(:)%natural%rhon    = .false.

  bcs(  1)%natural%rho   = .true.
  bcs(4:5)%natural%rho   = .true.
  bcs(  9)%natural%rho   = .true.
  bcs( 11)%natural%rho   = .true.
  bcs( 15)%natural%rho   = .true.
  bcs( 19)%natural%rho   = .true.

  bcs(  1)%natural%T     = .true.
  bcs(4:5)%natural%T     = .true.
  bcs(  9)%natural%T     = .true.
  bcs( 11)%natural%T     = .true.
  bcs( 15)%natural%T     = .true.
  bcs( 19)%natural%T     = .true.

  bcs(  1)%natural%Te    = .true.
  bcs(4:5)%natural%Te    = .true.
  bcs(  9)%natural%Te    = .true.
  bcs( 11)%natural%Te    = .true.
  bcs( 15)%natural%Te    = .true.
  bcs( 19)%natural%Te    = .true.

  bcs(  1)%natural%Ti    = .true.
  bcs(4:5)%natural%Ti    = .true.
  bcs(  9)%natural%Ti    = .true.
  bcs( 11)%natural%Ti    = .true.
  bcs( 15)%natural%Ti    = .true.
  bcs( 19)%natural%Ti    = .true.

  bcs(  1)%natural%vpar  = .true.
  bcs(4:5)%natural%vpar  = .true.
  bcs(  9)%natural%vpar  = .true.
  bcs( 11)%natural%vpar  = .true.
  bcs( 15)%natural%vpar  = .true.
  bcs( 19)%natural%vpar  = .true.

  bcs(  1)%natural%rhon  = .true.
  bcs(4:5)%natural%rhon  = .true.
  bcs(  9)%natural%rhon  = .true.
  bcs( 11)%natural%rhon  = .true.
  bcs( 15)%natural%rhon  = .true.
  bcs( 19)%natural%rhon  = .true.

  bcs(  1)%natural%rho_imp  = .true.
  bcs(4:5)%natural%rho_imp  = .true.
  bcs(  9)%natural%rho_imp  = .true.
  bcs( 11)%natural%rho_imp  = .true.
  bcs( 15)%natural%rho_imp  = .true.
  bcs( 19)%natural%rho_imp  = .true.
  ! -------------------------------------------

  
  U_sheath = .false.
  renormalise = .false.
  tauIC = 0.d0
  Wdia  = .false.

  zjz_0 =  0.1173d0   
  zjz_1 =  0.0d0   

  T_0   =  1.d-6  
  Ti_0  =  5.d-7
  Te_0  =  5.d-7

  T_1   =  1.d-8  
  Te_1  =  5.d-9
  Ti_1  =  5.d-9

  rho_0 =  1.d0   
  rho_1 =  1.d0   
  FF_0  =  1.d0
  FF_1  =  0.d0
  phi_0 =  0.d0
  phi_1 =  0.d0

  zj_coef     = 0.d0;  zj_coef(1)  = -1.d0
  T_coef      = 0.d0;  T_coef(1)   = -1.d0
  Te_coef     = 0.d0;  Te_coef(1)  = -1.d0
  Ti_coef     = 0.d0;  Ti_coef(1)  = -1.d0
  rho_coef    = 0.d0;  rho_coef(1) =  0.d0
  FF_coef     = 0.d0;  FF_coef(1)  = -1.d0
  dcoef       = 0.d0

  phi_coef    = 0.d0;  phi_coef(1) =  0.d0; phi_coef(4) = 1.d0
  nu_phi_source = 0.d0

  rhon_0 =  0.d0
  rhon_1 =  0.d0
  rhon_coef    = 0.d0
  rhon_coef(4) = 0.01
  rhon_coef(8) = 0.01

  pellet_amplitude  = 0.d0
  pellet_R          = 3.8d0
  pellet_Z          = 0.0d0
  pellet_phi        = 1.57d0
  pellet_theta      = 0.d0
  pellet_ellipse    = 5.d0
  pellet_radius     = 0.08d0
  pellet_sig        = 0.02
  pellet_length     = 0.785
  pellet_psi        = 1.0d0
  pellet_delta_psi  = 999.d0
  pellet_velocity_R = 0.d0
  pellet_velocity_Z = 0.d0
  pellet_particles  = 0.d0
  pellet_density    = 5.985d8       ! pellet density (in units 10^20 m^-3)
  pellet_density_bg = 5.958d8
  use_pellet        = .false.
  
  tstep_rst   = 0.d0
  t_now       = 0.d0
  t_start     = 0.d0
  index_start = 0

  nout = 9999999
  nout_projection = -1
  nout_particles  = 9999999

  rst_hdf5 = 1   ! =0,restart with binary files; =1, with HDF5 files

  !> Write out newest HDF5 restart file version this code supports, writing
  !! out an older version is possible by changing rst_hdf5_verison via the
  !! namelist input file
  rst_hdf5_version   = rst_hdf5_version_supported

  tokamak_device     = 'none'
  rho_file           = 'none'
  rhon_file          = 'none'
  T_file             = 'none'
  Te_file            = 'none'
  Ti_file            = 'none'
  phi_file           = 'none'
  Fprofile_file      = 'none'
  ffprime_file       = 'none'
  d_perp_file        = 'none'
  d_perp_imp_file    = 'none'
  zk_perp_file       = 'none'
  zk_e_perp_file     = 'none'
  zk_i_perp_file     = 'none'
  v_pinch_file       = 'none'
  R_Z_psi_bnd_file   = 'none'
  wall_file          = 'wall.txt'
  rot_file           = 'none'
  domm_file          = 'none'
  normalized_velocity_profile = .true.

  n_Fprofile_internal = 300 ! model710 only: size of internal numerical F-profile
  Fprofile_psi_max    = 1.5 ! model710 only: max-psi_norm of internal numerical F-profile
  Fprofile_tolerance  = 1.0 ! model710 only: tolerance of average different between final FFprime and requested FFprime

  produce_live_data  = .true.
  
  keep_n0_const      = .false.
  linear_run         = .false.

  use_zkperp_times_density = .false.
  zkperp_density_floor = 1.d-2
  
  export_for_nemec   = .false.
  export_aux_node_list = .true.
  
  ! Use iterative solver by default if n_tor>1.
  if ( n_tor == 1) then
    gmres            = .false.
  else
    gmres            = .true.
  end if
  
  gmres_max_iter     = 200                  ! Max number of GMRES iterations
  gmres_tol          = 1.d-8                ! converge tolerance GMRES
  gmres_4            = 1.d3                 ! error estimate GMRES (ratio preconditioned versus non-preconditioned error
  gmres_m            = 20                   ! gmres restart parameter
  iter_precon        = 10                   ! redo preconditioner when gmres iterations > iter_precon
  max_steps_noUpdate = 10000000             ! redo preconditioner when steps without preconditioning matrix update > max_steps_noUpdate
  centralize_harm_mat= .true.              ! centralize harmonic matrices on toroidal master rank 
  
  ! --- deprecated, code will stop if these parameters are set to .true. ---
  use_murge          = .false.
  use_murge_element  = .false.
  ! ------------------------------------------------------------------------

  tgnum              = 0.d0                 ! Taylor-Galerkin Stabilisation coefficients (0.d0 == TG not used)
  tgnum_psi          = 0.d0          
  tgnum_u            = 0.d0
  tgnum_zj           = 0.d0
  tgnum_w            = 0.d0
  tgnum_rho          = 0.d0
  tgnum_T            = 0.d0
  tgnum_Ti           = 0.d0
  tgnum_Te           = 0.d0
  tgnum_vpar         = 0.d0
  tgnum_rhon         = 0.d0
  tgnum_rhoimp       = 0.d0
  tgnum_nre          = 0.d0
  tgnum_AR           = 0.d0
  tgnum_AZ           = 0.d0
  tgnum_A3           = 0.d0

  keep_current_prof  = .true.               ! Keep the current_source term
  init_current_prof  = .false.
  current_prof_initialized = .false.
  
  use_mumps          = .false.              ! Use MUMPS solver
  use_pastix         = .false.              ! Use PASTIX solver
  use_strumpack      = .true.               ! Use STRUMPACK solver  
  use_wsmp           = .false.              ! Use WSMP solver (use with care, still in development!)
  
  use_mumps_eq       = .false.              ! Use MUMPS equilibrium solver
  use_pastix_eq      = .false.              ! Use PASTIX equilibrium solver
  use_strumpack_eq   = .false.              ! Use STRUMPACK equilibrium olver  
  
  use_mumps_prj      = .true.               ! Use MUMPS projection solver
  use_pastix_prj     = .false.              ! Use PASTIX projection solver
  use_strumpack_prj  = .false.              ! Use STRUMPACK projection olver  

  refinement         = .false.              ! enable mesh refinement
  force_central_node = .true.               ! force all nodes in the grid center to have the same values in flux surface aligned grids
  fix_axis_nodes     = .false.              ! Fix t-derivative and cross st-derivative on axis to avoid noise
  
  grid_to_wall       = .false.              ! extend the grid to a physical wall
  RZ_grid_inside_wall= .false.              ! build the rectangular grid inside first wall
  RZ_grid_jump_thres = 0.85                 ! threshold for jump of R-resolution as RZ-grid gets squeezed by limiter contour
  
  ! --- Option to manipulate psi_boundary, switched off by default
  manipulate_psi_map(:,1) = 0.
  manipulate_psi_map(:,2) = 99.
  manipulate_psi_map(:,3) = 99.
  manipulate_psi_map(:,4) = 0.1
  manipulate_psi_map(:,5) = 0.1
  
  adaptive_time      = .false.              ! requires no_mpi for Pastix library
  
  equil              = .true.               ! compute equilibrium
  
  no_mach1_bc        = .false.              ! Never apply Mach-1 BCs

  Mach1_openBC       = .true.               ! Full-MHD: Apply Mach-1 BCs inside mod_boundary_matrix_open.f90 (or mod_boundary_conditions.f90)
  Mach1_fix_B        = .true.               !< Full-MHD: Use the initial magnetic field for Mach1 BCs on targets, ie. without AR and AZ variations

  eta_ARAZ_const     = 0.d0                 !< Use uniform resistivity for AR and AZ equations, used only if eta_ARAZ_on=.false.
  eta_ARAZ_on        = .true.               !< Full-MHD: to switch on/off resistive   terms for AR and AZ equations
  eta_ARAZ_simple    = .false.              !< Full-MHD: remove the Fprof dependence of Bphi in the resistive terms for AR and AZ (which should be compensated by current source anyway)
  tauIC_ARAZ_on      = .true.               !< Full-MHD: to switch on/off diamagnetic terms for AR and AZ equations

  bench_without_plot = .false.              ! .true. for benchmark (mesuring elapsed time without plot phases) 

  mumps_ordering     = 7                    ! MUMPS ordering option (7:automatic, 3:Scotch, 4:PORD, 5:METIS)
  use_BLR_compression = .false.             ! Use MUMPS / PaStiX 6 solver with Block-low-rank (BLR) compression
  pastix_blr_abs_tol = .true.               ! Use absolute tolerance
  epsilon_BLR        = 0.                   ! Accuracy of BLR compression (0. = lossless)
  just_in_time_BLR   = .true.               ! Use Just-in-time strategy for BLR compression (.false. = memory-optimal)
  pastix_maxthrd     = 1024
  use_newton         = .false.
  maxNewton          = 20
  gamma_Newton       = 0.5
  alpha_Newton       = 2.d0
  strumpack_matching = .false.

  
!==== RMP parameters =====
  RMP_on             = .false.              ! .true. to activate RMPs (changes boundary conditions)
  RMP_psi_cos_file   = 'none'
  RMP_psi_sin_file   = 'none'
  RMP_growth_rate    = 0.011 ! RMP_growth_rate * RMP_ramp_up_time must be ~cst
  RMP_ramp_up_time   = 1000  ! in JOREK times
  output_bnd_elements = .false.  ! writes bnd nodes and elements in output files (boundary_nodes.dat and boundary_elements.dat)
  RMP_har_cos = 2
  RMP_har_sin = 3
  Number_RMP_harmonics = 1
  RMP_har_cos_spectrum(1) = RMP_har_cos ! 2 if only one harmonic (ntor=3) and this harmonic is RMP 
  RMP_har_sin_spectrum(1) = RMP_har_sin ! 3 if only one harmonic (ntor=3) and this harmonic is RMP 

! ===== Neoclassical parameters ======
  NEO = .false.
  neo_file ='none'
  amu_neo_const = 0.
  aki_neo_const = 0.
  

  n_limiter = 0
  R_limiter = 0.d0
  Z_limiter = 0.d0
  
  surface_cross_tol = 1.005
  eqdsk_psi_fact = 1.d0
  extend_existing_grid = .false.
  n_wall_blocks        = 0
  corner_block         = 0
  n_ext_block          = 0
  n_ext_equidistant    = .false.
  n_block_points_left  = 0
  R_block_points_left  = 0.d0
  Z_block_points_left  = 0.d0
  n_block_points_right = 0
  R_block_points_right = 0.d0
  Z_block_points_right = 0.d0
  use_simple_bnd_types = .false.
 
 !======================MB rotation profile
  V_0 = 0.d0
  V_1 = 0.d0
  V_coef(1:5) = (/ 0.d0, 0.d0, 0.d0, 0.1d0, 1.0d0 /)
!======================MB

  JET_MGI = .false.
  ASDEX_MGI = .false.
  t_ns  = 2.d3
  ns_amplitude = 0.d0
  ns_R      = 3.2d0
  ns_Z      =  1.5d0
  ns_phi    = 1.57d0
  ns_radius =   0.08d0
  ns_deltaphi =  0.5
  ns_delta_minor_rad = 0.d0
  ns_tor_norm = 1.
  ksi_ion = 1.84d-24
  D_neutral_x = 1.d-5
  D_neutral_y = 1.d-5
  D_neutral_p = 1.d-5
  delta_n_convection = 0
  nimp_bg = 0.

  n_adas = 1
  adas_dir = ' '
  imp_type = ' '
  index_main_imp = 0
  if (with_impurities) index_main_imp = 1
  use_imp_adas = .true. ! Directly use adas for impurity radiation; hard-coded one only implemented for argon
  drift_distance = 0.d0 ! No artificial plasmoid drift by default
  energy_teleported = 0.d0 
  constant_imp_source = 0.d0

  L_tube = 0. ! Needed to ensure injection starts at t_ns when JET_MGI=ASDEX_MGI=.false.

  !====== JET DMV-2 parameters
  !L_tube = 2.4d0
  K_Dmv = 4.d-2
  A_Dmv = 1.77d-2
  V_Dmv = 9.75d-4

  !======= Additional parameters for SPI =======
  spi_Vel_Rref    = 0.0d0
  spi_Vel_Zref    = 0.0d0
  spi_Vel_RxZref  = 0.0d0
  spi_quantity    = 0.0
  spi_quantity_bg = 0.0
  ns_radius_ratio = 1.4d0
  ns_radius_min   = 8.d-2
  spi_Vel_diff    = 0.0
  spi_angle       = 0.0
  spi_L_inj       = 0.25
  spi_L_inj_diff  = 0.0
  ns_phi_rotate   = 0.0
  tor_frequency   = 0.0
  n_spi           = 0
  n_spi(1)        = 1
  n_inj           = 1
  spi_rnd_seed    = 0
  spi_abl_model   = -1
  spi_shard_file(:) = 'none'
  spi_plume_file(:) = 'none'
  spi_plume_hdf5  = .false.
  spi_abl_mag_reduction  = .false.
  spi_tor_rot     = .false.
  spi_num_vol     = .true.
  using_spi       = .false.
  spi_abl_history_old = .false.

  output_prad_phi = .false.

!======================JP ECCD injection parameters
  nu_jec_fast=1.d1
  nu_jec1_fast=1.d1
  nu_jec2_fast=1.d1
  JJ_par=0.d1
  jecamp=1.d1
  jec_pos1=0.6d0
  jec_pos2=0.6d0
  jec_pos3=0.6d0
  jec_pos4=0.6d0
  jec_width=0.5d0
  jec_width2=0.5d0
  jw1=5.d-1 ! inner cut-off
  jw2=1.d0  ! outer cut-off
  jw3=1.d0  ! outer cut-off
  
  !> @name Mode families preconditioner parameters
  n_mode_families      = (n_tor + 1)/2
  autodistribute_modes = .true.
  mode_families_modes  = 0
  modes_per_family     = 0
  weights_per_family   = 1.0
  autodistribute_ranks = .true.
  ranks_per_family     = 0

!===================== Thermalization flag========

  thermalization = .true.

!===================== Polar axis treatment flag========

  treat_axis = .false.
  
!===================== not used?
  Q_bar = 0.d0
  Sigma = 0.d0

!===================== particle input values
use_particles      = .false.
nstep_particles    = 0
nsubstep_particles = 1
tstep_particles    = 1d-9
filter_perp        = 0.d0
filter_hyper       = 1.d-10
filter_par         = 0.d0
filter_perp_n0     = 0.d0
filter_hyper_n0    = 1.d-10
filter_par_n0      = 0.d0
restart_particles  = .false.
use_marker         = .false.
apply_dirichlet_proj = .false.
init_particles_only = .false.
find_RZ_nearby_iter = 16
find_RZ_nearby_tol  = 1.d-22

!--------------- valves -------------------------
valves(:)%type = 'none'
valves(:)%r_valve = -1.d0
valves(:)%R_valve_loc = -1.d0
valves(:)%Z_valve_loc = -1.d0
valves(:)%phi = -1.d0
do i=1, n_valves_max
  valves(i)%poly_R = 0.d0
  valves(i)%poly_Z = 0.d0
enddo

! -------------- particle groups ---------------
n_part_groups = 0
part_groups_in_use(:) = 'non'
proj_collection_period = 1

part_group_configs(:)%Z                 = 1
part_group_configs(:)%mass              = 0.d0
part_group_configs(:)%coupling_scheme   = 'non'
part_group_configs(:)%n_particles       = 0.d0
part_group_configs(:)%type              = 'none'
part_group_configs(:)%id                = 'non'
part_group_configs(:)%init_function     = 'none'
part_group_configs(:)%init_pdf          = 'none'
part_group_configs(:)%do_conservation_checks = .false.

!----- specific to ics and ncs 
part_group_configs(:)%atom_data_suffix      = ''
part_group_configs(:)%use_kin_puffing        = .false.
part_group_configs(:)%use_kin_radiation      = .false.
part_group_configs(:)%use_kin_ionisation     = .false.
! --- ncs only
part_group_configs(:)%use_kin_recombination  = .false.
part_group_configs(:)%use_kin_cx             = .false.
part_group_configs(:)%use_kin_neutral_coll   = .false.
do i=1,n_part_groups_max
  part_group_configs(i)%neutral_coll_dTw(:)  = -1.d99
end do
part_group_configs(:)%ncoll_each_nstep_part  = -9999991
! --- ics only
part_group_configs(:)%use_kin_bg_collisions  = .false.
part_group_configs(:)%kin_bg_coll_type       = 'Homma2020'
part_group_configs(:)%homma2020_alpha        = 1.5d0
part_group_configs(:)%ics_group_idx          = -1

!----- specific to rep 
part_group_configs(:)%num_re                 = 0.d0
part_group_configs(:)%re_energy              = 0.d0
part_group_configs(:)%re_std_energy          = 0.d0
part_group_configs(:)%re_pitch               = 0.d0

!----- specific to epf
part_group_configs(:)%T_maxwell              = 0.d0
part_group_configs(:)%n_phi_planes           = 0
part_group_configs(:)%n_particles_total      = 0.d0



do i=1, n_part_groups_max
  do j=1, n_valves_max
    part_group_configs(i)%puff_ctrl(j)%supers_num_puff    = -1
    part_group_configs(i)%puff_ctrl(j)%supers_weight_puff = -1.d0
    part_group_configs(i)%puff_ctrl(j)%supers_ratio_puff  = -1.d0   
    !< if none of these three above options are set, the supers_ratio_puff method
    !< will be used, with its default value being set by supers_ratio_puff_default in mod_particle_puffing.f90
    !< which overrides the default value of supers_ratio_puff set here
    part_group_configs(i)%puff_ctrl(j)%times = -1.d0
    part_group_configs(i)%puff_ctrl(j)%rates = -1.d0
    part_group_configs(i)%puff_ctrl(j)%from_file = "none"
  enddo
enddo

do i=1, n_part_groups_max
  do j=1, n_part_groups_max
    part_group_configs(i)%wall_act_configs(j)%type            = "none"
    part_group_configs(i)%wall_act_configs(j)%target_group_id = "non"
    part_group_configs(i)%wall_act_configs(j)%weight_factor   = 1.d0
    part_group_configs(i)%wall_act_configs(j)%only_in_polygon = .false.
    part_group_configs(i)%wall_act_configs(j)%poly_R          = -1.d99
    part_group_configs(i)%wall_act_configs(j)%poly_Z          = -1.d99
    write(tmp,"(I5)") j
    write(part_group_configs(i)%wall_act_configs(j)%nametag,"(A)") adjustl(trim(tmp)) ! sets the default nametag to the wall_act_configs(j) number j

    part_group_configs(i)%wall_act_configs(j)%supers_num_wall    = -1
    part_group_configs(i)%wall_act_configs(j)%supers_weight_wall = -1.d0
    part_group_configs(i)%wall_act_configs(j)%supers_ratio_wall  = -1.d0   
    !< if none of these three above options are set, the supers_ratio_wall method
    !< will be used, with its default value being set by supers_ratio_wall_default in mod_particle_wall_interaction.f90
  enddo
  part_group_configs(i)%wall_act_each_nstep_part = -9999991
enddo

part_kill_ratio = 1.d-3

! --- fluid groups
n_fluid_groups = 0

fluid_configs(:)%Z = -999
fluid_configs(:)%density_fraction = -1.d99
fluid_configs(1)%density_fraction = 1.d0 !< the first one should have default 1, if more then one fluid is used, the user should specify the distribution

do i=1, n_fluid_groups_max
  do j=1, n_part_groups_max
    fluid_configs(i)%wall_act_configs(j)%type            = "none"
    fluid_configs(i)%wall_act_configs(j)%target_group_id = "non"
    fluid_configs(i)%wall_act_configs(j)%weight_factor   = 1.d0
    fluid_configs(i)%wall_act_configs(j)%only_in_polygon = .false.
    fluid_configs(i)%wall_act_configs(j)%poly_R          = -1.d99
    fluid_configs(i)%wall_act_configs(j)%poly_Z          = -1.d99
    write(tmp,"(I5)") j
    write(fluid_configs(i)%wall_act_configs(j)%nametag,"(A)") adjustl(trim(tmp))

    fluid_configs(i)%wall_act_configs(j)%supers_num_wall    = -1
    fluid_configs(i)%wall_act_configs(j)%supers_weight_wall = -1.d0
    fluid_configs(i)%wall_act_configs(j)%supers_ratio_wall  = -1.d0
  enddo
enddo
!-----------------------------------------------

use_manual_random_seed = .false.
manual_seed = 498932990          !< chosen arbitarily
use_fixed_rng_value = .false.
fixed_rng_value = 0.5

end subroutine preset_parameters
