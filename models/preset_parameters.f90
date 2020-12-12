!> Presets input parameters to reasonable default values.
!!
!! The model-specific routines initialise_parameters may overwrite
!! these defaults according to the requirements of the respective
!! model.
subroutine preset_parameters
  
  use phys_module
  use mumps_module,  only: no_zeros_mumps, mumps_ordering
  use pastix_module, only: no_zeros_pastix, pastix_smp_only
  
  implicit none
  
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
  visco_T_dependent = .true.
  ZKpar_T_dependent = .true.

  eta_num_T_dependent   = .false.
  visco_num_T_dependent = .false.

  eta           = 1.d-5
  T_max_eta     = 1.d3
  eta_ohmic     = 0.d0
  T_max_eta_ohm = 1.d3

  visco = 1.d-5
  visco_par = 1.d-5
  
  central_density = 1.d0        ! the central density in units 10^20 m^-3
  central_mass    = 2.d0        ! the central average ion mass (D)

  restart      = .false.
  import_equil = .false.
  regrid       = .false.
  rst_format   = 0             ! use 'old' format for restart import
  write_ps     = .true.           ! write postscript file at the end of the run 
  
  freeboundary_equil = .false. ! use free or fixed boundary equilibrium
  freeboundary       = .false. ! use free or fixed boundary?
  resistive_wall     = .false. ! use a resistive or ideal wall?    (freeboundary only)
  freeb_equil_iterate_area = .false.
  freeb_change_indices = .true. ! exchange grid node indices to parallelize boundary integral

  bc_natural_flux    = .false.! boundary conditions for flux surface boundaries (2 and 3)
  bc_natural_open    = .false. ! use sheath (Bohm) boundary conditions
  gamma_sheath       = 4.5d0  ! sheath transmission factor (single fluid) in the JOREK definition
  gamma_stangeby     = -1.d99 ! sheath transmission factor (single fluid) given by Stangeby
  density_reflection = 0.d0   ! reflection coefficient for outgoing density
  neutral_reflection = 0.d0   ! reflection coefficient for (fluid) neutrals
  imp_reflection     = 0.d0   ! reflection coefficient for (fluid) impurities
  
  deuterium_adas        = .false. 
  old_deuterium_atomic  = .false. 
  mach_one_bnd_integral = .false. ! implement Mach one condition as boundary integral
  Vpar_smoothing        = .false. ! smooth the transitions of Vpar positive/negavtive at B.n
  Vpar_smoothing_coef   = (/0.01d0, 0.d0, 0.d0 /) !(/ 0.01d0, 0.016d0, 0.00575446347d0/)

  amix                 = 0.d0
  amix_freeb           = 0.85d0
  equil_accuracy       = 1.d-6
  equil_accuracy_freeb = 1.d-6
  axis_srch_radius     = 99.d0
  
  n_R          = 0
  n_Z          = 0

  n_radial     = 11
  n_pol        = 16

  n_flux       = 11
  n_tht        = 16
  n_tht_equidistant = .false.
  
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

  psi_axis_init = -0.1d0
  XR_r(:)       = 999.d0
  SIG_r(:)      = 999.d0
  XR_tht(:)     = 999.d0
  SIG_tht(:)    = 999.d0

  SIG_closed  = 0.1d0
  SIG_open    = 0.1d0
  SIG_outer   = 0.1d0
  SIG_inner   = 0.1d0
  SIG_private = 0.1d0
  SIG_up_priv = 0.1d0
  SIG_theta   = 0.03d0
  SIG_leg_0   = 0.05d0
  SIG_leg_1   = 0.2d0
  SIG_up_leg_0= 0.05d0
  SIG_up_leg_1= 0.2d0
  
  dPSI_open    = 0.11
  dPSI_outer   = 0.11
  dPSI_inner   = 0.11
  dPSI_private = 0.03
  dPSI_up_priv = 0.03
  
  SDN_threshold = 1.d-4
  
  R_geo     = 10.d0
  Z_geo     = 0.d0
  amin      = 1.d0

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

  bootstrap = .false.

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

  xr1  = 9999.d0
  sig1 = 9999.d0
  xr2  = 99999.d0
  sig2 = 99999.d0

  R_begin = -0.1d0
  R_end   =  0.1d0
  Z_begin = -0.1d0
  Z_end   = 0.1d0
  
  ZK_i_perp(:) = 0.d0
  ZK_i_perp(1) = 1.d-5; ZK_i_perp(2) = 0.d0; ZK_i_perp(3)= 0.d0; ZK_i_perp(4)= 99.d0; ZK_i_perp(5) = 99.d0
  ZK_e_perp(:) = 0.d0
  ZK_e_perp(1) = 1.d-5; ZK_e_perp(2) = 0.d0; ZK_e_perp(3)= 0.d0; ZK_e_perp(4)= 99.d0; ZK_e_perp(5) = 99.d0
  ZK_par     = 1.d0
  ZK_par_max = 1.d20
  D_par      = 0.d0
  ZK_perp(1:5) = (/ 1.d-5, 0.d0, 0.d0, 99.d0, 99.d0 /)
  ZK_par       = 1.d0
  ZK_par_max   = 1.d20
  D_perp(1:5)  = (/ 1.d-5, 0.d0, 0.d0, 99.d0, 99.d0 /)
  D_par        = 0.d0
  D_perp_imp(1:5)  = (/ 1.d-5, 0.d0, 0.d0, 99.d0, 99.d0 /)
  D_par_imp        = 0.d0

  D_prof_neg         = 1.d-5
  D_prof_neg_thresh  = 0.d0 ! default is zero for keeping the old behavior
  ZK_prof_neg        = 1.d-5
  ZK_par_neg         = 1.d-3
  ZK_prof_neg_thresh = 0.d0 ! default is zero for keeping the old behavior
  ZK_par_neg_thresh  = 0.d0

  ne_SI_min          = 1.d18
  Te_eV_min          = 5.
  rn0_min            = 1.d-8

  T_min              =-1.0d20
  rho_min            =-1.0d20
  
  corr_neg_temp_coef(:) = (/ 0.5, 0.5 /)
  corr_neg_dens_coef(:) = (/ 0.5, 0.5 /)

  eta_num       = 0.d0
  visco_num     = 0.d0
  visco_par_num = 0.d0
  D_perp_num    = 0.d0
  ZK_perp_num   = 0.d0
  Dn_perp_num   = 0.d0

  heatsource          = 1.e-7
  heatsource_i        = 5.d-8
  heatsource_e        = 5.d-8
  heatsource_psin     = 1.0d0
  heatsource_sig      = 0.1d0
  particlesource      = 1.e-5
  particlesource_psin = 1.0d0
  particlesource_sig  = 0.1d0
  edgeparticlesource      = 0.d0
  edgeparticlesource_psin = 0.98
  edgeparticlesource_sig  = 0.01
  heatsource_gauss          = 0.d0
  heatsource_gauss_psin     = 0.9d0
  heatsource_gauss_sig      = 0.1d0
  particlesource_gauss      = 0.d0
  particlesource_gauss_psin = 0.9d0
  particlesource_gauss_sig  = 0.1d0
  neutral_line_source       = 0.d0
  neutral_line_R_start      = 1.d20
  neutral_line_Z_start      = 1.d20
  neutral_line_R_end        = 2.d20
  neutral_line_Z_end        = 2.d20

  
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
  
  zj_coef     = 0.d0;  zj_coef(1)  = -1.d0
  T_coef      = 0.d0;  T_coef(1)   = -1.d0
  Te_coef     = 0.d0;  Te_coef(1)  = -1.d0
  Ti_coef     = 0.d0;  Ti_coef(1)  = -1.d0
  rho_coef    = 0.d0;  rho_coef(1) =  0.d0
  FF_coef     = 0.d0;  FF_coef(1)  = -1.d0

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
  pellet_density    = 3.d8       ! pellet density (in units 10^20 m^-3)
  pellet_density_bg = 3.d8
  use_pellet        = .false.
  
  t_now       = 0.d0
  t_start     = 0.d0
  index_start = 0

  nout = 9999999

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
  Fprofile_file      = 'none'
  ffprime_file       = 'none'
  d_perp_file        = 'none'
  d_perp_imp_file    = 'none'
  zk_perp_file       = 'none'
  zk_e_perp_file     = 'none'
  zk_i_perp_file     = 'none'
  R_Z_psi_bnd_file   = 'none'
  wall_file          = 'none'
  rot_file           = 'none'
  normalized_velocity_profile = .true.

  n_Fprofile_internal = 300 ! model710 only: size of internal numerical F-profile
  Fprofile_psi_max    = 1.5 ! model710 only: max-psi_norm of internal numerical F-profile
  Fprofile_tolerance  = 1.0 ! model710 only: tolerance of average different between final FFprime and requested FFprime

  produce_live_data  = .true.
  
  keep_n0_const      = .false.
  linear_run         = .false.
  
  export_for_nemec   = .false.
  
  gmres              = .true.               ! Use iterative solver
  gmres_max_iter     = 200                  ! Max number of GMRES iterations
  gmres_tol          = 1.d-8                ! converge tolerance GMRES
  gmres_4            = 1.d3                 ! error estimate GMRES (ratio preconditioned versus non-preconditioned error
  gmres_m            = 20                   ! gmres restart parameter
  iter_precon        = 10                   ! redo preconditioner when gmres iterations > iter_precon
  max_steps_noUpdate = 10000000             ! redo preconditioner when steps without preconditioning matrix update > max_steps_noUpdate
  centralize_harm_mat= .false.              ! centralize harmonic matrices on toroidal master rank 
  
  ! --- deprecated, code will stop if these parameters are set to .true. ---
  use_murge          = .false.
  use_murge_element  = .false.
  ! ------------------------------------------------------------------------

  tgnum              = 0.d0                 ! Taylor-Galerkin Stabilisation coefficients (0.d0 == TG not used)

  keep_current_prof  = .true.               ! Keep the current_source term
  
  use_mumps          = .false.              ! Use MUMPS solver
  use_pastix         = .true.               ! Use PASTIX solver
  use_strumpack      = .false.              ! Use STRUMPACK solver  
  use_wsmp           = .false.              ! Use WSMP solver (use with care, still in development!)
  
  use_mumps_eq       = .false.              ! Use MUMPS equilibrium solver
  use_pastix_eq      = .false.              ! Use PASTIX equilibrium solver
  use_strumpack_eq   = .false.              ! Use STRUMPACK equilibrium olver  
  
  refinement         = .false.              ! enable mesh refinement
  force_central_node = .true.               ! force all nodes in the grid center to have the same values in flux surface aligned grids
  fix_axis_nodes     = .false.              ! Fix t-derivative and cross st-derivative on axis to avoid noise
  
  grid_to_wall       = .false.              ! extend the grid to a physical wall
  RZ_grid_inside_wall= .false.              ! build the rectangular grid inside first wall
  
  ! --- Option to manipulate psi_boundary, switched off by default
  manipulate_psi_map(:,1) = 0.
  manipulate_psi_map(:,2) = 99.
  manipulate_psi_map(:,3) = 99.
  manipulate_psi_map(:,4) = 0.1
  manipulate_psi_map(:,5) = 0.1
  
  adaptive_time      = .false.              ! requires no_mpi for Pastix library
  
  equil              = .true.               ! compute equilibrium
  
  Mach1_openBC       = .true.               ! Full-MHD: Apply Mach-1 BCs inside mod_boundary_matrix_open.f90 (or mod_boundary_conditions.f90)

  eta_ARAZ_on        = .true.               !< Full-MHD: to switch on/off resistive   terms for AR and AZ equations
  tauIC_ARAZ_on      = .true.               !< Full-MHD: to switch on/off diamagnetic terms for AR and AZ equations

  fix_axis_nodes     = .false.              !< Fix t-derivative on axis to avoid noise)
  
  bench_without_plot = .false.              ! .true. for benchmark (mesuring elapsed time without plot phases) 
  no_zeros_pastix    = .false.              ! .true. to remove nonzeros in the preconditioning matrix with MUMPS
  no_zeros_mumps     = .false.              ! .true. to remove nonzeros in the preconditioning matrix with PaStiX

  mumps_ordering     = 7                    ! MUMPS ordering option (7:automatic, 3:Scotch, 4:PORD, 5:METIS)
  use_BLR_compression = .false.             ! Use MUMPS / PaStiX 6 solver with Block-low-rank (BLR) compression
  pastix_blr_abs_tol = .true.               ! Use absolute tolerance
  epsilon_BLR        = 0.                   ! Accuracy of BLR compression (0. = lossless)
  just_in_time_BLR   = .true.               ! Use Just-in-time strategy for BLR compression (.false. = memory-optimal)

  
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
  
  n_wall_blocks        = 0
  n_ext_block          = 0
  n_block_points_left  = 0
  R_block_points_left  = 0.d0
  Z_block_points_left  = 0.d0
  n_block_points_right = 0
  R_block_points_right = 0.d0
  Z_block_points_right = 0.d0
 
 !======================MB rotation profile
  V_0 = 0.d0
  V_1 = 0.d0
  V_coef(1:5) = (/ 0.d0, 0.d0, 0.d0, 0.1d0, 1.0d0 /)
!======================MB

  JET_MGI = .false.
  ASDEX_MGI = .false.
  ns_amplitude = 0.d0
  ns_R      = 3.2d0
  ns_Z      =  1.5d0
  ns_phi    = 1.57d0
  ns_radius =   0.08d0
  ns_sig    =  0.05
  ns_deltaphi =  0.5
  ns_tor_norm = 1.
  ksi_ion = 1.84d-24
  D_neutral_x = 1.d-5
  D_neutral_y = 1.d-5
  D_neutral_p = 1.d-5
  delta_n_convection = 0
  nimp_bg = 0.
  !====== JET DMV-2 parameters
  L_tube = 2.4d0
  K_Dmv = 4.d-2
  A_Dmv = 1.77d-2
  V_Dmv = 9.75d-4
  t_ns  = 2.d3
  !======= Additional parameters for SPI =======
  spi_Vel_Rref    = 0.0d0
  spi_Vel_Zref    = 0.0d0
  spi_Vel_RxZref  = 0.0d0
  spi_quantity    = 0.0
  spi_quantity_bg = 0.0
  ng_radius_ratio = 1.4d0
  ng_radius_min   = 8.d-2
  spi_Vel_diff    = 0.0
  spi_angle       = 0.0
  spi_L_inj       = 0.25
  ns_phi_rotate   = 0.0
  tor_frequency   = 0.0
  n_spi           = 0
  n_spi_tot       = 0
  n_inj           = 1
  spi_rnd_seed    = 0
  spi_abl_model   = 0
  spi_shard_file  = 'none'
  spi_tor_rot     = .false.
  using_spi       = .false.

  output_rad_phi  = .false.
  n_adas          = 0
  adas_dir        = ''
  gas_type        = ''

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

end subroutine preset_parameters
