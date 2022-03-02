!> Input parameters and physical variables.
module phys_module
  
  use mod_parameters
  use constants
  use data_structure              !< Added in order to dynamically allocate pellets
  use mod_openadas
  use mod_coronal

  implicit none
  
  !> @name Various parameters
  real*8  :: eta                  !< Resistivity at plasma cener (normalized)
  real*8  :: eta_T_0              !< Initial resistivity
  real*8  :: eta_ohmic            !< Resistivity at core for the Ohmic heating term
  logical :: eta_T_dependent      !< Resistivity dependent on temperature? Otherwise constant
  real*8  :: T_max_eta            !< Temperature above which the resistivity is truncated (use with care; only for numerical reasons)
  real*8  :: T_max_eta_ohm        !< Temperature above which the resistivity used in the Ohmic heating term is truncated (use with care; only for numerical reasons)
  real*8  :: visco                !< Viscosity at plasma center (normalized)
  real*8  :: visco_rst            !< visco value from restart file
  real*8  :: visco_par_rst        !< visco_par value from restart file
  real*8  :: eta_rst              !< eta value from restart file
  logical :: visco_T_dependent    !< Viscosity dependent on temperature? Otherwise constant.
  real*8  :: visco_par            !< Parallel viscosity (normalized)
  real*8  :: F0                   !< Determines fixed toroidal magnetic field: \f$ B_\phi = F_0/R \f$
  real*8  :: central_density      !< particle density at the magnetic axis (in units of \f$10^{20} m^{-3}\f$)
  real*8  :: central_mass         !< average ion mass in atomic mass units (constant in time and space)
  real*8  :: sqrt_mu0_rho0        !< Normalization factor \f$\sqrt(\mu_0 \rho_0)\f$ calculated from input
  real*8  :: sqrt_mu0_over_rho0   !< Normalization factor \f$\sqrt(\mu_0/\rho_0)\f$ calculated from input
  real*8  :: gamma                !< ratio of specific heat (typically 5/3)
  real*8  :: Q_bar                !< (model400)
  real*8  :: sigma                !< (model400)
  real*8  :: tauIC                !< Scaling factor for diamagnetic terms (see [[diamag|diamagnetic]])
  real*8  :: tauIC_nominal        !< Nominal scaling factor (considering Ti=Te) for diamagnetic terms (see [[diamag|diamagnetic]])
  real*8  :: eta_spitzer          !< Spitzer resistivity in the core (considering main ion charge Z=1, effective ion charge Zeff=1)
  logical :: Wdia                 !< Include diamagnetic flows in viscosity terms? (see [[wdia|here]])
  logical :: U_sheath             !< Use Stangeby BCs for electric potential
  logical :: renormalise          !< Set true to give all input MHD parameters in S.I. units (ie. renormalise them before equations)
  real*8  :: gamma_sheath         !< sheath boundary condition on open fieldlines (JOREK units); you can also provide gamma_stangeby in normal units instead!
  real*8  :: gamma_stangeby       !< Sheath tranmission coefficient given by P. Stangeby in (The plasma boundary of magnetic fusion devices)
  real*8  :: gamma_sheath_e       !< sheath boundary condition on open fieldlines (JOREK units); you can also provide gamma_stangeby in normal units instead!
  real*8  :: gamma_e_stangeby     !< Sheath tranmission coefficient given by P. Stangeby in (The plasma boundary of magnetic fusion devices)
  real*8  :: gamma_sheath_i       !< sheath boundary condition on open fieldlines (JOREK units); you can also provide gamma_stangeby in normal units instead!
  real*8  :: gamma_i_stangeby     !< Sheath tranmission coefficient given by P. Stangeby in (The plasma boundary of magnetic fusion devices)
  real*8  :: density_reflection   !< density reflection coeefficient on open fieldlines
  real*8  :: neutral_reflection   !< reflection coefficient of ions into neutrals (model500)
  real*8  :: imp_reflection       !< impurity reflection coefficient on open fieldlines
  logical :: old_deuterium_atomic !< use old fit to calculate atomic coefficients for D (ionization, recombination, radiation), otherwise a better fit is used
  logical :: deuterium_adas       !< use OPEN ADAS to calculate ionization, recombination and radiation coeffients for deuterium                        
  logical :: deuterium_adas_1e20  !< use OPEN ADAS with fixed density=1e20 to calculate ionization, recombination and radiation coeffients for deuterium
  logical :: mach_one_bnd_integral!< use a boundary integral (boundary_matrix_open) to implement Mach=one boundary condition
  logical :: vpar_smoothing       !< apply a smoothing function to smooth jumps in Vpar at B.n=0
  real*8  :: vpar_smoothing_coef(3) !< coefficients for the smoothing profile of the parallel velocity
  real*8  :: min_sheath_angle     !< For sheath boundary conditions: Minimum incident angle for heat and particle fluxes (in degrees)
  integer :: mode(n_tor)          !< Toroidal mode number corresponding to the JOREK modes, e.g., for n_period=8 and n_tor=3, mode(:)=0,8,8
  integer :: nout                 !< Output a restart file every nout timesteps
  integer :: xcase                !< 1->LowerXpoint. 2->UpperXpoint. 3->doubleNull
  real*8  :: SDN_threshold        !< threshold, in absolute psi, for a symmetric-double-null grid construction
  integer :: rst_format           !< 0 == old format, 1 == new format for restart file
  logical :: restart              !< Restart a code run from the restart file jorek_restart.h5?
  logical :: regrid               !< Re-generate the flux-aligned grid (does not work currently)?
  logical :: import_equil         !< (presently unused)
  logical :: xpoint               !< X-point plasma or not? see also xcase
  logical :: bootstrap            !< Evolve the Bootstrap current consistently with time?
  real*8  :: minRad               !< Approximation of minor radius for bootstrap current calculation
  logical :: refinement           !< Use mesh refinement? (not presently available)
  logical :: force_central_node   !< Force all nodes in the center to have the same values in flux aligned grids or independent values?
  logical :: fix_axis_nodes       !< Fix t-derivative and cross st-derivative on axis to avoid noise
  logical :: bc_natural_flux      !< boundary conditions for flux surface boundaries (2 and 3)
  logical :: bc_natural_open      !< use natural boundary conditions on the open fieldlines
  logical :: produce_live_data    !< Write data 'macroscopic_vars.dat' during the code run allowing to use plot_live_data.sh?
  logical :: grid_to_wall         !< extend the grid to a physical wall
  logical :: RZ_grid_inside_wall  !< build the rectangular grid inside first wall
  real*8  :: RZ_grid_jump_thres   !< threshold to change R-resolution as RZ-grid gets sqeezed by limiter contour
  real*8  :: manipulate_psi_map(5,5) !< Option to manipulate Psi_boundary for the initial grid
  logical :: adaptive_time        !< (presently not useful)
  logical :: equil                !< compute equilibrium
  logical :: no_mach1_bc          !< Never apply Mach-1 BCs
  logical :: Mach1_openBC         !< Full-MHD: Apply Mach-1 BCs inside mod_boundary_matrix_open.f90 (or mod_boundary_conditions.f90)

  ! --- RESISTIVITY SWITCHES FOR AR AND AZ EQUATIONS
  ! --- 1.
  ! --- Default set-up is eta_ARAZ_on = .true.
  ! --- In this case, the same resistivity is used for all AR,AZ,A3 components
  ! --- 2.
  ! --- If (eta_ARAZ_on = .true.) and (eta_ARAZ_simple = .true.), then the Fprof component of Bphi is removed from the resistive term.
  ! --- Note that in a stationary equilibrium, the Fprof component of Bphi from the current and the current source should exactly cancel anyway.
  ! --- However, this Fprof component can lead to strong noise at high resistivity, so it is often better to remove it.
  ! --- 3.
  ! --- If (eta_ARAZ_on = .false.), then resistivity is switched off by default for the AR and AZ equations.
  ! --- However, regardless of which eta-model is used for A3, (ie. eta_T_dependent or not), you can still set eta_ARAZ_const 
  ! --- which will use a constant resistivity at the level eta_ARAZ_const.
  ! --- 4.
  ! --- Here again, setting (eta_ARAZ_simple = .true.) with eta_ARAZ_const removes the Fprof dependence of Bphi in the resistive term and the current term for AR&AZ
  real*8  :: eta_ARAZ_const       !< Use uniform resistivity for AR and AZ equations, used only if eta_ARAZ_on=.false.
  logical :: eta_ARAZ_on          !< Full-MHD: to switch on/off resistive terms for AR and AZ equations
  logical :: eta_ARAZ_simple      !< Full-MHD: remove the Fprof dependence of Bphi in the resistive terms for AR and AZ (which should be compensated by current source anyway)

  logical :: tauIC_ARAZ_on        !< Full-MHD: to switch on/off diamagnetic terms for AR and AZ equations
  logical :: bench_without_plot   !< if .true., do not produce certain output plots (e.g., for benchmarking)
  logical :: gmres                !< Use iterative GMRES solver
  integer :: gmres_max_iter       !< Maximum number of GMRES iterations
  logical :: keep_n0_const        !< Perform a linear run where the equilibrium quantities (i_tor=1) do not change with time?
  logical :: linear_run           !< Same as keep_n0_const, to be replaced soon by true linear run where modes are independent
  logical :: export_for_nemec     !< Export equilibrium information for the NEMEC code?
  logical :: use_murge            !< (Deprecated, Cannot be used any more)
  logical :: use_murge_element    !< (Deprecated, Cannot be used any more)
  logical :: use_BLR_compression  !< Use Block-Low-Rank (BLR) compression in MUMPS / PaStiX 6 solvers
  real*8  :: epsilon_BLR          !< Accuracy of BLR compression
  logical :: just_in_time_BLR     !< Use Just-in-time strategy for BLR compression (speed optimized)
  logical :: pastix_blr_abs_tol   !< Use absolute tolerance for BLR
  logical :: write_ps             !< Write postscript file at the end of the run
  logical :: use_mumps            !< Use Mumps solver
  logical :: use_pastix           !< Use Pastix solver
  logical :: use_strumpack        !< Use Strumpack solver
  logical :: use_mumps_eq         !< Use Mumps equilibrium solver
  logical :: use_pastix_eq        !< Use Pastix equilibrium solver
  logical :: use_strumpack_eq     !< Use Strumpack equilibrium solver  
  logical :: use_wsmp             !< Use WSMP solver
  logical :: centralize_harm_mat  !< Centralize harmonic matrices on toridal master ranks; switch for STRUMPACK solver
  real*8  :: prev_FB_fact = 1.d0  !< FB_factor that had been applied when importing the restart file

  ! ------------------------------------------------
  ! --- Structures to implement BCs in model600
  ! ------------------------------------------------
  ! --- For more info see  https://www.jorek.eu/wiki/doku.php?id=choose_boundary_conditions
  integer, parameter :: max_bnd_types=30

  type type_dirichlet_bc                           
    logical :: psi  
    logical :: u    
    logical :: zj   
    logical :: w    
    logical :: rho  
    logical :: T    
    logical :: Ti   
    logical :: Te   
    logical :: Vpar 
    logical :: rhon 
    logical :: rho_imp 
    logical :: nre  
    logical :: AR   
    logical :: AZ   
    logical :: A3  
  end type type_dirichlet_bc

  type type_natural_bc                           
    logical :: rho  
    logical :: T    
    logical :: Ti   
    logical :: Te   
    logical :: Vpar 
    logical :: rhon 
    logical :: rho_imp 
    logical :: nre  
  end type type_natural_bc

  type type_bcs                           
    type (type_dirichlet_bc) :: dirichlet
    type (type_natural_bc)   :: natural
    logical                  :: mach1 
  end type type_bcs

  type (type_bcs), dimension(max_bnd_types) :: bcs   
  ! ------------------------------------------------


  character(20)       :: numfmt     = "'_d',i5.5"
  character(20)       :: numfmt_rst = "'_r',i3.3"
  ! Identity of the processor
  integer  :: pglobal_id
  
  real*8, allocatable :: energies(:,:,:)   !< Magnetic and kinetic mode energies at timesteps.
  real*8, allocatable :: energies2(:,:,:)  !< global density and temperature at timesteps.
  real*8, allocatable :: energies3(:,:,:)  !< global currents (general and total eccd) at timesteps.
  real*8, allocatable :: energies4(:,:,:)  !< global applied eccd currents j1 and j2 at timesteps.

  character(len=3)    :: mode_type(n_tor) !< 'cos' or 'sin'
  
  !> Points used as limiters (see routine find_limiter)
  integer, parameter :: max_limiter = 1000 !< Maximum number of limiter points
  integer :: n_limiter                     !< Number of limiter points
  real*8  :: R_limiter(max_limiter)        !< R-positions of the limiter points
  real*8  :: Z_limiter(max_limiter)        !< Z-positions of the limiter points
  integer :: first_target_point		   !< index of the first target point on the limiter (for xpoint_grid_wall)
  integer :: last_target_point		   !< index of the last  target point on the limiter (does NOT need to be > first_target_point)
  
  !> Points used as blocks to extend grid into complex wall structures, see https://www.jorek.eu/wiki/doku.php?id=wallgrid_tutorial
  real*8  :: eqdsk_psi_fact                                                     !< multiply eqdsk psi by factor for grid_inside_wall
  logical :: extend_existing_grid                                               !< Add patches to existing grid from restart file
  integer, parameter :: n_wall_blocks_max = 30                                  !< Maximum number of blocks (30 should be enough)
  integer :: n_wall_blocks                                                      !< Number of blocks
  integer, parameter :: n_wall_block_points_max = 20                            !< Max number of blocks points
  integer :: corner_block(n_wall_blocks_max)                                    !< =1 for a corner block ("left" side will also be wall-aligned)
  integer :: n_ext_block(n_wall_blocks_max)                                     !< Number of 'radial' grid points from the outermost flux surface to wall)
  logical :: n_ext_equidistant(n_wall_blocks_max)                               !< if true, radial spacing of grid points will be equidistant (not adapted)
  integer :: n_block_points_left (n_wall_blocks_max)                            !< Number of points on left side of block
  real*8  :: R_block_points_left (n_wall_blocks_max,n_wall_block_points_max)    !< R-positions of points on left side of block
  real*8  :: Z_block_points_left (n_wall_blocks_max,n_wall_block_points_max)    !< Z-positions of points on left side of block
  integer :: n_block_points_right(n_wall_blocks_max)                            !< Number of points on left side of block
  real*8  :: R_block_points_right(n_wall_blocks_max,n_wall_block_points_max)    !< R-positions of points on left side of block
  real*8  :: Z_block_points_right(n_wall_blocks_max,n_wall_block_points_max)    !< Z-positions of points on left side of block
  logical :: use_simple_bnd_types                                               !< convert Stan's bnd_types to Guido's bnd_types
  
  !> @name Define X-point geometry by geometrical properties
  !!
  !! \f[
  !! \Psi(\theta) =
  !!        -x_{shift}\sin(\theta)
  !!        +x_{left}\cos(\theta)
  !!        +x_{ampl}\left[
  !!            \left(\frac{x_{width}\cdot(\theta-x_{theta})}{x_{sig}}\right)^2-1
  !!          \right]exp\left[-\left(\frac{\theta-x_{theta}}{x_{sig}}\right)^2\right]
  !! \f]
  real*8  :: xampl    !< Allows to construct simple X-point cases by coefficients (modifies Psi boundary condition)
  real*8  :: xwidth   !< Allows to construct simple X-point cases by coefficients (modifies Psi boundary condition)
  real*8  :: xsig     !< Allows to construct simple X-point cases by coefficients (modifies Psi boundary condition)
  real*8  :: xtheta   !< Allows to construct simple X-point cases by coefficients (modifies Psi boundary condition)
  real*8  :: xshift   !< Allows to construct simple X-point cases by coefficients (modifies Psi boundary condition)
  real*8  :: xleft    !< Allows to construct simple X-point cases by coefficients (modifies Psi boundary condition)
  
  !> @name Heat and particle sources
  !!
  !! \f[
  !! S(\Psi_N) = S_0 \cdot \left[0.5 - 0.5 \tanh\left(\frac{\Psi_N - \Psi_{N,0}}{\sigma}\right) \right]
  !! \f]
  !!
  !! The following parameters can be set via the namelist input file:
  !! - \f$ S_0 \f$ denotes the source strenght (e.g., heatsource)
  !! - \f$ \Psi_{N,0} \f$ denotes the position around which the source is ramped down (e.g., heatsource_psin)
  !! - \f$ \sigma \f$ denotes the width over which the source is ramped down (e.g., heatsource_sig)
  !!
  real*8  :: particlesource            !< Particle source amplitude
  real*8  :: particlesource_psin       !< Position around which the source is ramped down
  real*8  :: particlesource_sig        !< Width over which the source is ramped down
  real*8  :: particlesource_gauss      !< Additional Gaussian particle source amplitude
  real*8  :: particlesource_gauss_psin !< Position around which Gaussian source is set
  real*8  :: particlesource_gauss_sig  !< Width over which Gaussian source is set
  real*8  :: edgeparticlesource        !< Edge particle source amplitude
  real*8  :: edgeparticlesource_psin   !< Position around which the edge particle source is located
  real*8  :: edgeparticlesource_sig    !< Width over which edge particle source extends
  real*8  :: neutral_line_source(10)   !< neutral inflow source
  real*8  :: neutral_line_R_start(10)  !< neutral inflow source (starting point of line source)
  real*8  :: neutral_line_Z_start(10)  !< neutral inflow source
  real*8  :: neutral_line_R_end(10)    !< neutral inflow source (end point of line source)
  real*8  :: neutral_line_Z_end(10)    !< neutral inflow source
  real*8  :: heatsource                !< Heat source amplitude
  real*8  :: heatsource_e              !< Electron heat source amplitude
  real*8  :: heatsource_i              !< Ion heat source amplitude
  real*8  :: heatsource_psin           !< Position around which the source is ramped down
  real*8  :: heatsource_sig            !< Width over which the source is ramped down
  real*8  :: heatsource_e_psin         !< Position around which the electron source is ramped down
  real*8  :: heatsource_e_sig          !< Width over which the electron source is ramped down
  real*8  :: heatsource_i_psin         !< Position around which the ion source is ramped down
  real*8  :: heatsource_i_sig          !< Width over which the ion source is ramped down
  real*8  :: heatsource_gauss          !< Additional Gaussian heat source amplitude
  real*8  :: heatsource_gauss_psin     !< Position around which Gaussian source is located
  real*8  :: heatsource_gauss_sig      !< Width over which Gaussian source extends
  real*8  :: heatsource_gauss_e        !< Gaussian heat source for electrons
  real*8  :: heatsource_gauss_i        !< Gaussian heat source for ions
  real*8  :: heatsource_gauss_e_psin   !< Position around which electrons Gaussian source is located
  real*8  :: heatsource_gauss_e_sig    !< Width over which electrons Gaussian source extends
  real*8  :: heatsource_gauss_i_psin   !< Position around which ions Gaussian source is located
  real*8  :: heatsource_gauss_i_sig    !< Width over which ions Gaussian source extends
  
  !> @name Hyper-resistivity, -viscosity and -diffusivities
  real*8  :: eta_num, visco_num, visco_par_num, D_perp_num, Zk_perp_num, Dn_perp_num, Zk_i_perp_num, Zk_e_perp_num
  !> @name Shock-capturing terms
  logical :: use_sc  !< Use shock-capturing stabilization
  real*8  :: D_perp_sc_num, D_par_sc_num, Dn_pol_sc_num, Dn_p_sc_num
  real*8  :: ZK_perp_sc_num, ZK_par_sc_num, ZK_i_perp_sc_num, ZK_i_par_sc_num, ZK_e_perp_sc_num, ZK_e_par_sc_num
  real*8  :: visco_sc_num, visco_par_sc_num
  logical :: eta_num_T_dependent     !< Hyper-resistivity dependent on temperature? Otherwise constant.
  logical :: eta_num_psin_dependent  !< Give profile for Hyper-resistivity as function of \psi_N? Useful for 2D current flattening
  real*8  :: eta_num_prof(10)        !< Coefficients to specify \psi_N profile for hyper-resistivity
  logical :: visco_num_T_dependent!< Hyper-visocsity dependent on temperature? Otherwise constant.
  logical :: add_sources_in_sc    !< Whether to add effect of sources in shock-capturing stabilization or not

  !> @name Timestepping parameters
  real*8  :: tstep             		!< Size of the timesteps (\f$ \Delta t \f$)
  real*8  :: tstep_prev                 !< Previous time-step if using variable dt Gears
  real*8  :: tstep_n(10)       		!< Alternative to tstep: Up to ten values may be given
  integer :: nstep             		!< Number of timesteps to perform
  integer :: nstep_n(10)       		!< Alternative to nstep: Up to ten values may be given
  real*8  :: t_start           		!< Time value at the start of the code run (zero or from restart file)
  real*8  :: t_now             		!< Current time value in the simulation
  integer :: index_start       		!< Time step index at the beginning of the code run (zero or from restart file)
  integer :: index_now         		!< Current time step index
  real*8, allocatable :: xtime(:) 	!< Time values corresponding to the timesteps.
  character(len=80) :: time_evol_scheme !< Time evolution scheme to use (see [[time-integration|time_integration]])
  real*8  :: time_evol_theta   		!< Time evolution parameter theta (see [[time-integration|time_integration]])
  real*8  :: time_evol_zeta    		!< Time evolution parameter zeta (see [[time-integration|time_integration]])

  integer :: rst_hdf5                   !< Write hdf5 restart files if set to 1
  integer :: rst_hdf5_version           !< Write which version of hdf5 files?
  integer, parameter :: rst_hdf5_version_supported = 2 !< What is the highest version number supported?
  
  !> @name Machine name
  character(len=512) :: tokamak_device 	!< Name of the tokamak device we are simulating

  !> @name Analytical boundary of initial grid
  !!
  !! Analytical definition of the boundary of the non flux-aligned initial polar grid.
  !!
  !! - \f$ Z=Z_{geo} + a_{min} \epsilon \sin(\theta) \f$
  !!
  !! - for \f$ \theta < \pi \f$:
  !!   \f$ R=R_{geo} + a_{min} \cos\left[\theta+T_u\sin(\theta)+Q_u\sin(2\theta)\right] \f$
  !!
  !! - for \f$ \theta \ge \pi \f$:
  !!   \f$ R=R_{geo} + a_{min} \cos\left[\theta+T_l\sin(\theta)+Q_l\sin(2\theta)\right] \f$
  !!
  real*8  :: amin              !< Minor radius for polar grid construction, set to 1 if boundary is specified with R,Z points
  real*8  :: ellip             !< Ellipticity of polar grid (see analytical definition in phys_module.f90)
  real*8  :: tria_u            !< Upper triangularity of polar grid (see analytical definition in phys_module.f90)
  real*8  :: tria_l            !< Lower triangularity of polar grid (see analytical definition in phys_module.f90)
  real*8  :: quad_u            !< Upper quadrangularity of polar grid (see analytical definition in phys_module.f90)
  real*8  :: quad_l            !< Lower quadrangularity of polar grid (see analytical definition in phys_module.f90)
  
  !> @name Fourier expanded boundary of initial grid
  !! Boundary of the non flux-aligned initial polar grid given as Fourier series
  integer, parameter :: n_bnd_max = 3000 	!< Max number of entries in boundary points
  integer :: mf              		 	!< Number of entries in fbnd and fpsi
  real*8  :: fbnd(n_bnd_max)        		!< Fourier expansion of boundary
  real*8  :: fpsi(n_bnd_max)        		!< Fourier expansion of the poloidal flux at the boundary
  
  !> @name Numerical boundary of initial grid
  !! Numerical definition of the boundary of the non flux-aligned initial polar grid.
  integer :: n_boundary       			!< Number of points in R_boundary, Z_boundary, psi_boundary.
  real*8  :: R_boundary  (n_bnd_max)		!< Numerical R values defining the boundary
  real*8  :: Z_boundary  (n_bnd_max)		!< Numerical Z values defining the boundary
  real*8  :: psi_boundary(n_bnd_max)		!< Numerical values giving the poloidal flux at the boundary
  
  !> @name PF coils definition for initial equilibrium (MAST)
  !! Numerical definition of the PF coils definition for initial equilibrium (MAST)
  integer :: n_pfc            !< Number of coils, (OLD. for MAST...) use JOREK-STARWALL for coils instead [[jorek-starwall|JOREK-STARWALL]]
  real*8  :: Rmin_pfc(40)     !< Minimum R of coil, (OLD. for MAST...) use JOREK-STARWALL for coils instead [[jorek-starwall|JOREK-STARWALL]]
  real*8  :: Rmax_pfc(40)     !< Maximum R of coil, (OLD. for MAST...) use JOREK-STARWALL for coils instead [[jorek-starwall|JOREK-STARWALL]]
  real*8  :: Zmin_pfc(40)     !< Minimum Z of coil, (OLD. for MAST...) use JOREK-STARWALL for coils instead [[jorek-starwall|JOREK-STARWALL]]
  real*8  :: Zmax_pfc(40)     !< Maximum Z of coil, (OLD. for MAST...) use JOREK-STARWALL for coils instead [[jorek-starwall|JOREK-STARWALL]]
  real*8  :: current_pfc(40)  !< Current density in the coil, (OLD. for MAST...) use JOREK-STARWALL for coils instead [[jorek-starwall|JOREK-STARWALL]]
  
  !> @name current ropes definition for initial equilibrium (eg. merging flux ropes)
  !! Numerical definition of current ropes for initial equilibrium (eg. merging flux ropes)
  integer :: n_jropes          !< Number of ropes, 
  real*8  :: R_jropes(10)      !< R centre of rope
  real*8  :: Z_jropes(10)      !< Z centre of rope
  real*8  :: w_jropes(10)      !< width of rope
  real*8  :: current_jropes(10)!< Current inside the rope
  real*8  :: rho_jropes(10)    !< Density inside the rope
  real*8  :: T_jropes(10)      !< Temperature inside the rope
  
  !> @name Pellet-related input parameters
  real*8  :: pellet_amplitude  !< amplitude of density source (when pellet modelled as density source)
  real*8  :: pellet_R          !< major radius position pellet
  real*8  :: pellet_Z          !< Z position pellet
  real*8  :: pellet_phi        !< width of the pellet cloud (density source) in toroidal angle
  real*8  :: pellet_ellipse    !< the ellipticity of the pellet source
  real*8  :: pellet_radius     !< radius of the simulation pellet
  real*8  :: pellet_sig        !< width of smoothing of density source (arctan((r-pellet_radius)/pellet_sig))
  real*8  :: pellet_length     !< width of smoothing of density source in toroidal angle
  real*8  :: pellet_theta      !< orientation of the pellet ellipse
  real*8  :: pellet_psi        !< pellet_width in poloidal flux
  real*8  :: pellet_delta_psi  !< width of smoothing in poloidal flux
  real*8  :: pellet_velocity_R !< pellet velocity component radial direction
  real*8  :: pellet_velocity_Z !< pellet velocity component Z direction
  real*8  :: pellet_density    !< pellet atom number density (in units \f$10^{20} m^{-3}\f$)
  real*8  :: pellet_density_bg !< background species pellet atom number density (in units \f$10^{20} m^{-3}\f$)
  real*8  :: pellet_particles  !< the number of particles in the pellet (in units of \f$10^{20}\f$)
  logical :: use_pellet

  !> @name shared between MGI and SPI applications
  integer, parameter :: n_inj_max = 10 ! The hard coded maximum number of injections
  integer, parameter :: n_imp_max = 5  ! The hard coded maximum number of impurity species

  real*8  :: t_ns(n_inj_max)   !< MGI onset time (JOREK units)
  real*8  :: ns_amplitude(n_inj_max)  !< Amplitude of gas source
  real*8  :: ns_R(n_inj_max)   !< R position of gas source
  real*8  :: ns_Z(n_inj_max)   !< Z position of gas source
  real*8  :: ns_phi(n_inj_max) !< Phi position of gas source
  real*8  :: ns_radius         !< Poloidal radius of gas source
  real*8  :: ns_deltaphi       !< Toroidal extension of gas source
  real*8  :: ns_deltaminrad    !< Extension of gas source in the minor radial direction (if greater than 0.)
  real*8  :: ns_tor_norm       !< Gas source normalization factor related to its toroidal shape

  character(len=80) :: imp_type(n_imp_max) !< Type of injected material or background impurity species: Argon, neon, ...
  logical :: use_imp_adas       !< Use open adas to calculate ionization, recombination and radiation coeffients for impurities

  !> @name Massive gas injection-related input parameters
  
  logical :: JET_MGI           !< Switch to use a JET-like MGI
  logical :: ASDEX_MGI         !< Switch to use an ASDEX-like MGI
  real*8  :: V_Dmv             !< Volume of the DMV reservoir
  real*8  :: P_Dmv             !< Pressure in the DMV reservoir (bar)
  real*8  :: A_Dmv             !< Cross sectional area of DMV (Disruption mitigation valve) pipe
  real*8  :: K_Dmv             !< Correction parameter describing the gas expansion near the pipe orifice
  real*8  :: L_tube            !< Pipe length
  real*8  :: ksi_ion            !< Energy cost of each ionization
  real*8  :: delta_n_convection !< Switch to activate the convection term for neutrals (at the plasma velocity)
  real*8  :: nimp_bg(n_imp_max) !< Density of background impurities (in \f$m^{-3}\f$)
  integer :: index_main_imp     !< Index of the main impurity species (in imp_type and nimp_bg) solved with continuity equation
                               
  !> @name Shattered Pellet Injection related input parameters
  ! Note that the SPI share many of the MGI parameters. The code should return to simple MGI upon using_spi = false
  ! The reference spatial coordinate for shattered pellets are calculated using ns_R etc. 
  ! More information on the wiki: https://www.jorek.eu/wiki/doku.php?id=spi_tutorial
  logical :: using_spi          !< This determines whether to use SPI or traditional MGI; see [[spi_tutorial|SPI Tutorial]]
  real*8  :: spi_Vel_Rref(n_inj_max)   !< Reference velocity of pellet center along R upon injection
  real*8  :: spi_Vel_Zref(n_inj_max)   !< Reference velocity of pellet center along Z upon injection
  real*8  :: spi_Vel_RxZref(n_inj_max) !< Reference velocity of pellet center along RxZ direction upon injection
  real*8  :: spi_quantity(n_inj_max)   !< Total injected atom number for impurity SPI
  real*8  :: spi_quantity_bg(n_inj_max)!< Total injected atom number for background species SPI
  real*8  :: ng_radius_ratio           !< We are assuming a constant ratio between the radius of NG clouds
                                       !< and that of shattered pellets

  real*8  :: spi_Vel_diff(n_inj_max)   !< The velocity difference from the reference velocity
  real*8  :: spi_angle                 !< The vertex angle of spi spreading in terms of rad
  real*8  :: spi_L_inj(n_inj_max)      !< Distance between SPI nozzle and ns_R, ns_Z, ns_phi
  real*8  :: spi_L_inj_diff(n_inj_max) !< The position difference with respect to the point (ns_R, ns_Z, ns_phi)
  real*8  :: ns_phi_rotate             !< The toroidal position of rotated injection point
  real*8  :: tor_frequency             !< The rigid body rotation frequency

  real*8  :: ng_radius_min      !< This defines the minimum radius of neutral cloud for numerical reasons (in m)

  real*8, allocatable  :: xtime_spi_ablation(:,:)         !< The time history of SPI ablation
  real*8, allocatable  :: xtime_spi_ablation_rate(:,:)    !< The time history of SPI ablation rate
  real*8, allocatable  :: xtime_spi_ablation_bg(:,:)      !< The time history of SPI ablation for background species
  real*8, allocatable  :: xtime_spi_ablation_bg_rate(:,:) ! <The time history of SPI ablation rate for bg species

  real*8, allocatable  :: xtime_radiation(:)    !< The time history of radiated energy in SI unit
  real*8, allocatable  :: xtime_rad_power(:)    !< The time history of radiated power in SI unit

  real*8, allocatable  :: xtime_E_ion(:)        !< The time history of the ionization potential energy in SI unit
  real*8, allocatable  :: xtime_E_ion_power(:)  !< Time derivative of xtime_E_ion
  real*8, allocatable  :: xtime_P_ei(:)         !< The time history of electron-ion energy exchange power

  integer :: n_spi(n_inj_max)   !< Number of shattered fragment injected for each injection
  integer :: n_spi_tot          !< Total number of shattered fragments injected
  integer :: n_inj              !< Number of injections
  integer :: spi_abl_model      !< Determine which type of ablation model is used.
                                !< 0 for constant release rate, 1 for NGS model,
                                !< 2 for Sergeev formula, 3 for Parks formula.
                                !< For details see Nucl. Fusion 61 (2021) 026015 (23pp), 
                                !< https://iopscience.iop.org/article/10.1088/1741-4326/abcbcb
  integer :: spi_rnd_seed(40)   !< Random seed array used for the generation of the SPI velocity spread

  character(len=256) :: spi_shard_file(n_inj_max)!< The name of the shard size file
  character(len=256) :: spi_plume_file(n_inj_max)!< The name of the shard information datafile (array)
  logical            :: spi_plume_hdf5           !< if 'spi_plume_file' is in HDF5format?

  integer :: n_adas             !< Number of species to be traced by ADAS

  logical :: spi_tor_rot        !< Flag to turn on a rigid body toroidal plasma rotation for SPI
  logical :: spi_num_vol        !< Flag to turn on numerical integration of the gas source volumes from SPI

  type (type_SPI), allocatable :: pellets(:) !< Each element corresponds to one injected pellet (shard)

  character(len=512)            :: adas_dir    !< The directory of ADAS data file to be read
  type (adf11_all), allocatable :: imp_adas(:) !< The ADAS data for impurities
  type (coronal), allocatable   :: imp_cor(:)  !< The coronal equilibrium distribution of impurities

  logical :: output_prad_phi    !< Output Prad(phi) into a file using integrals_3D
  
  !> @name Fix boundary equilibrium parameters
  real*8  :: amix              !< Mix Poisson solution with previous one with a given factor
  real*8  :: equil_accuracy    !< Tolerance of the convergence for the fix-boundary equilibrium
  real*8  :: axis_srch_radius  !< Magnetic axis will be searched inside a circle with this radius
  real*8  :: delta_psi_GS      !< Expected psi_bnd - psi_axis for the final equilibrium  
  logical :: newton_GS_fixbnd  !< Newton instead of Picard iterations for fixed-boundary equilibria?
  logical :: newton_GS_freebnd !< Newton instead of Picard iterations for free-boundary equilibria?
 
  !> @name Free boundary extension
  !! Input parameters related to the free boundary extension (folder vacuum/).
  logical :: freeboundary_equil      !< use a free or fixed boundary equilibrium? ([[jorek-starwall|JOREK-STARWALL]])
  logical :: freeboundary            !< use free or fixed boundary conditions in time-evolution? ([[jorek-starwall|JOREK-STARWALL]])
  logical :: resistive_wall          !< use a resistive or ideal wall? ([[jorek-starwall|JOREK-STARWALL]])
  logical :: freeb_equil_iterate_area !< iterate to a target area during freeboundary equilibrium limiter cases [[jorek-starwall-faqs|jorek_starwall]]
  real*8  :: amix_freeb              !< choose amix for freeboundary equilibrium
  real*8  :: equil_accuracy_freeb    !< Tolerance of the convergence for the freeboundary equilibrium
  logical :: freeb_change_indices    !< Exchange grid node indices to parallelize boundary integral
  
  !> @name Rectangular Grid
  !! Parameters defining a rectangular grid in R- and Z-directions in the poloidal plane.
  integer :: n_R               !< Number of grid points in R-direction (for rectangular grid) (see also [[grids#tutorials|here]])
  integer :: n_Z               !< Number of grid points in Z-direction (for rectangular grid)
  real*8  :: R_begin           !< Left boundary of grid in R-direction (for rectangular grid)
  real*8  :: R_end             !< Right boundary of grid in R-direction (for rectangular grid)
  real*8  :: Z_begin           !< Lower boundary of grid in Z-direction (for rectangular grid)
  real*8  :: Z_end             !< Upper boundary of grid in Z-direction (for rectangular grid)
  real*8  :: rect_grid_vac_psi !< Use a vacuum psi-bnd condition for squared-grid, ie. (rect_grid_vac_psi * R**2)

  
  !> @name Polar Grid
  !! Parameters defining a non flux-aligned polar grid in the poloidal plane.
  logical :: force_horizontal_Xline !< Force the grid line through Xpoint to be horizontal (instead of perp. to line between Xpoint and axis)
  integer :: n_radial          	    !< Number of radial grid points (for polar grid) (see also [[grids|here]])
  integer :: n_pol             	    !< Number of poloidal grid points (for polar grid)
  real*8  :: R_geo             	    !< Center of the grid (for polar grid)
  real*8  :: Z_geo             	    !< Center of the grid (for polar grid)
  real*8  :: psi_axis_init     	    !< Initial guess for Psi at the magnetic axis (for polar grid)
  real*8  :: XR_r(2)           	    !< Psi_N position of radial grid accumulation (two positions) (for polar grid) (also used for R-position in square-grid)
  real*8  :: SIG_r(2)          	    !< Width of grid accumulation (two positions) (for polar grid) (also used for R-width in square-grid)
  real*8  :: XR_tht(2)         	    !< Position of poloidal grid accumulation (0...1, two positions) (for polar grid)
  real*8  :: SIG_tht(2)        	    !< Width of grid accumulation (two positions) (for polar grid)
  real*8  :: XR_z(2)           	    !< Z-position of square grid accumulation (two positions) (for square grid)
  real*8  :: SIG_z(2)          	    !< Z-Width of grid accumulation (two positions) (for square grid)
  real*8  :: bgf_r, bgf_z           !< Background for meshac distribution or R and Z accumulation (only for square grid!)
  
  !> @name Flux surface grid
  !! Parameters defining a flux-aligned grid without X-point in the poloidal plane.
  integer :: n_flux            !< Number of radial grid points (for flux-aligned grid) (see also [[grids#tutorials|here]])
  integer :: n_tht             !< Number of poloidal grid points (for flux-aligned grid)
  real*8  :: xr1               !< Grid accumulation parameter (for flux-aligned grid)
  real*8  :: xr2               !< Grid accumulation parameter (for flux-aligned grid)
  real*8  :: sig1              !< Grid accumulation parameter (for flux-aligned grid)
  real*8  :: sig2              !< Grid accumulation parameter (for flux-aligned grid)
  
  !> @name Flux surface grid with X-point
  !! Parameters defining a flux-aligned grid with X-point in the poloidal plane.
  integer :: n_open            !< Number of 'radial' grid points in the open flux region - between the two separatrices if double-null
  integer :: n_outer           !< Number of 'radial' grid points in the open flux region on the outer side (LFS) if double-null
  integer :: n_inner           !< Number of 'radial' grid points in the open flux region on the inner side (HFS) if double-null
  integer :: n_private         !< Number of 'radial' grid points in the private flux region at the bottom
  integer :: n_leg             !< Number of 'poloidal' grid points along the divertor legs at the bottom
  integer :: n_leg_out         !< Number of 'poloidal' grid points along the divertor legs at the bottom on the LFS
  integer :: n_up_priv         !< Number of 'radial' grid points in the private flux region at the top (upper Xpoint or double-null)
  integer :: n_up_leg          !< Number of 'poloidal' grid points along the divertor legs at the top (upper Xpoint or double-null)
  integer :: n_up_leg_out      !< Number of 'poloidal' grid points along the divertor legs on the top on the LFS (upper Xpoint or double-null)
  integer :: n_ext             !< Number of 'radial' grid points from the outermost flux surface to wall)
  logical :: n_tht_equidistant !< switch on to get an equidistant poloidal distribution of elements in the core of the grid (psi<0.5)
  real*8  :: SIG_closed        !< Width with grid accumulation (for flux-aligned grid)
  real*8  :: SIG_open          !< Width with grid accumulation (for flux-aligned grid)
  real*8  :: SIG_outer         !< Width with grid accumulation (for flux-aligned grid)
  real*8  :: SIG_inner         !< Width with grid accumulation (for flux-aligned grid)
  real*8  :: SIG_private       !< Width with grid accumulation (for flux-aligned grid)
  real*8  :: SIG_up_priv       !< Width with grid accumulation (for flux-aligned grid)
  real*8  :: SIG_theta         !< Width with grid accumulation (for flux-aligned grid)
  real*8  :: SIG_leg_0         !< Width with grid accumulation (for flux-aligned grid)
  real*8  :: SIG_leg_1         !< Width with grid accumulation (for flux-aligned grid)
  real*8  :: SIG_up_leg_0      !< Width with grid accumulation (for flux-aligned grid)
  real*8  :: SIG_up_leg_1      !< Width with grid accumulation (for flux-aligned grid)
  real*8  :: dPSI_open         !< Delta Psi grid extends into the open flux region (for flux-aligned grid)
  real*8  :: dPSI_outer        !< Delta Psi grid extends into the open flux region (for flux-aligned grid)
  real*8  :: dPSI_inner        !< Delta Psi grid extends into the open flux region (for flux-aligned grid)
  real*8  :: dPSI_private      !< Delta Psi grid extends into the private flux region (for flux-aligned grid)
  real*8  :: dPSI_up_priv      !< Delta Psi grid extends into the private flux region (for flux-aligned grid)
  
  !> @name Analytical heat, particle and neutral particles diffusivity parameters
  real*8  :: D_perp(10)    = 0.d0 !< Coefficients for perpendicular particle diffusion profile
  real*8  :: D_par                !< Parallel particle diffusion (usually not useful)
  real*8  :: D_perp_imp(10)= 0.d0 !< Coefficients for perpendicular imp particle diffusion profile
  real*8  :: D_par_imp            !< Parallel impurity particle diffusion (usually not useful)
  real*8  :: ZK_perp(10)   = 0.d0 !< Coefficients for perpendicular heat diffusion profile
  real*8  :: ZK_par               !< Parallel heat diffusion value in the plasma center
  real*8  :: ZK_par_max           !< Do not use larger parallel heat diffusion values for numerical reasons
  real*8  :: ZK_par_SpitzerHaerm  !< Spitzer-Haerm parallel heat diffusion value in the plasma center (assuming a Z=1 plasma with Te=Ti)
  real*8  :: ZK_i_perp(10) = 0.d0 !< Coefficients for perpendicular ion heat diffusion profile
  real*8  :: ZK_e_perp(10) = 0.d0 !< Coefficients for perpendicular electron heat diffusion profile
  real*8  :: ZK_i_par             !< Ion parallel heat diffusion coefficient in the plasma center
  real*8  :: ZK_e_par             !< Electron parallel heat diffusion coefficient in the plasma center
  real*8  :: ZK_i_par_SpitzerHaerm!< Spitzer-Haerm ion parallel heat diffusion value in the plasma center (assuming a Z=1 plasma)
  real*8  :: ZK_e_par_SpitzerHaerm!< Spitzer-Haerm electron parallel heat diffusion value in the plasma center (assuming a Z=1 plasma)
  real*8  :: D_neutral_x          !< Neutral particle diffusivity in R-direction
  real*8  :: D_neutral_y          !< Neutral particle diffusivity in Z-direction
  real*8  :: D_neutral_p          !< Neutral particle diffusivity in phi-direction
  logical :: ZKpar_T_dependent    !< Use a temperature dependent parallel heat diffusivity

  !> @name Numerical heat and particle diffusivity profiles
  character(len=512)  :: d_perp_file        !< ASCII file with perpendicular particle diffusion profile
  character(len=512)  :: d_perp_imp_file    !< ASCII file with perpendicular particle diffusion profile
  character(len=512)  :: zk_perp_file       !< ASCII file with perpendicular heat diffusion profile
  character(len=512)  :: zk_e_perp_file     !< ASCII file with perpendicular electron heat diffusion profile
  character(len=512)  :: zk_i_perp_file     !< ASCII file wtih perpendicular ion heat diffusion profile
  logical             :: num_d_perp         !< automatically set true if d_perp_file /= 'none'
  logical             :: num_d_perp_imp     !< automatically set true if d_perp_file /= 'none'
  logical             :: num_zk_perp        !< automatically set true if zk_perp_file /= 'none'
  logical             :: num_zk_e_perp      !< automatically set true if zk_e_perp_file /= 'none'
  logical             :: num_zk_i_perp      !< automatically set true if zk_i_perp_file /= 'none'
  integer             :: num_d_perp_len     !< Number of datapoints in d_perp profile
  integer             :: num_d_perp_len_imp !< Number of datapoints in d_perp profile for impurity
  integer             :: num_zk_perp_len    !< Number of datapoints in zk_perp profile
  integer             :: num_zk_e_perp_len  !< Number of datapoints in zk_e_perp profile
  integer             :: num_zk_i_perp_len  !< Number of datapoints in zk_i_perp profile
  real*8, allocatable :: num_d_perp_x(:)    !< Psi_N values of d_perp  profile
  real*8, allocatable :: num_d_perp_y(:)    !< D_perp values of d_perp profile
  real*8, allocatable :: num_d_perp_x_imp(:)!< Psi_N values of d_perp  profile for impurity
  real*8, allocatable :: num_d_perp_y_imp(:)!< D_perp values of d_perp profile for impurity
  real*8, allocatable :: num_zk_perp_x(:)   !< Psi_N values of zk_perp profile
  real*8, allocatable :: num_zk_perp_y(:)   !< ZK_perp values of zk_perp profile
  real*8, allocatable :: num_zk_e_perp_x(:) !< Psi_N values of zk_e_perp profile
  real*8, allocatable :: num_zk_e_perp_y(:) !< ZK_perp values of zk_e_perp profile
  real*8, allocatable :: num_zk_i_perp_x(:) !< Psi_N values of zk_i_perp profile
  real*8, allocatable :: num_zk_i_perp_y(:) !< ZK_perp values of zk_i_perp profile
  
  !> @name Analytical input profile for the density
  real*8  :: rho_0             !< Central normalized density (usually 1)
  real*8  :: rho_1             !< SOL normalized density
  real*8  :: rho_coef(10)      !< Density profile coefficients
  
  !> @name Numerical input profile for the density
  character(len=512)  :: rho_file        !< ASCII file the density profile is read from.
  logical             :: num_rho         !< automatically set true if rho_file /= 'none'
  integer             :: num_rho_len     !< Number of points in rho profile
  real*8, allocatable :: num_rho_x(:)    !< Psi_N values of rho profile points
  real*8, allocatable :: num_rho_y0(:)   !< Density values of rho profile
  real*8, allocatable :: num_rho_y1(:)   !< First derivatives of density profile (\f$ d\rho/d\Psi_N \f$)
  real*8, allocatable :: num_rho_y2(:)   !< Second derivatives of density profile (\f$ d^2\rho/d\Psi_N^2 \f$)
  real*8, allocatable :: num_rho_y3(:)   !< Third derivatives of density profile (\f$ d^3\rho/d\Psi_N^3 \f$)

  !> @name Analytical input profile for the temperature
  real*8  :: T_0            !< Central normalized temperature
  real*8  :: T_1            !< SOL normalized temperature
  real*8  :: T_coef(10)     !< Temperature profile coefficients
  real*8  :: Ti_0           !< Central ion normalized temperature
  real*8  :: Ti_1           !< SOL ion normalized temperature
  real*8  :: Ti_coef(10)    !< Ion temperature profile coefficients
  real*8  :: Te_0           !< Central ion normalized temperature
  real*8  :: Te_1           !< SOL ion normalized temperature
  real*8  :: Te_coef(10)    !< Ion temperature profile coefficients
  
  !> @name Numerical input profile for the temperature
  character(len=512)  :: T_file          !< ASCII file the temperature profile is read from.
  logical             :: num_T           !< automatically set true if T_file /= 'none'
  integer             :: num_T_len       !< Number of points in T profile
  real*8, allocatable :: num_T_x(:)      !< PsiN values of T profile points (PsiN values)
  real*8, allocatable :: num_T_y0(:)     !< Temperature values of T profile
  real*8, allocatable :: num_T_y1(:)     !< First derivatives of temperature profile (\f$ dT/d\Psi_N \f$)
  real*8, allocatable :: num_T_y2(:)     !< Second derivatives of temperature profile (\f$ d^2T/d\Psi_N^2 \f$)
  real*8, allocatable :: num_T_y3(:)     !< Third derivatives of temperature profile (\f$ d^3T/d\Psi_N^3 \f$)
  
  !> @name Numerical input profile for the ion temperature (model400)
  character(len=512)  :: Ti_file         !< ASCII file the ion temperature profile is read from.
  logical             :: num_Ti          !< is set true if T_file /= 'none'
  integer             :: num_Ti_len      !< Number of points in profile
  real*8, allocatable :: num_Ti_x(:)     !< Radial positions of profile points (PsiN values)
  real*8, allocatable :: num_Ti_y0(:)    !< Values of temperature profile
  real*8, allocatable :: num_Ti_y1(:)    !< First derivatives of temperature profile (\f$ dT/d\Psi_N \f$)
  real*8, allocatable :: num_Ti_y2(:)    !< Second derivatives of temperature profile (\f$ d^2T/d\Psi_N^2 \f$)
  real*8, allocatable :: num_Ti_y3(:)    !< Third derivatives of temperature profile (\f$ d^3T/d\Psi_N^3 \f$)
  
  !> @name Numerical input profile for the electron temperature (model400)
  character(len=512)  :: Te_file         !< ASCII file the electron temperature profile is read from.
  logical             :: num_Te          !< is set true if T_file /= 'none'
  integer             :: num_Te_len      !< Number of points in profile
  real*8, allocatable :: num_Te_x(:)     !< Radial positions of profile points (PsiN values)
  real*8, allocatable :: num_Te_y0(:)    !< Values of temperature profile
  real*8, allocatable :: num_Te_y1(:)    !< First derivatives of temperature profile (\f$ dT/d\Psi_N \f$)
  real*8, allocatable :: num_Te_y2(:)    !< Second derivatives of temperature profile (\f$ d^2T/d\Psi_N^2 \f$)
  real*8, allocatable :: num_Te_y3(:)    !< Third derivatives of temperature profile (\f$ d^3T/d\Psi_N^3 \f$)  
  
  !> @name Analytical input profile for the neutral density (model 500)
  real*8  :: rhon_0           !< Central value for the initial normalized neutral density
  real*8  :: rhon_1           !< SOL value for the initial normalized neutral density
  real*8  :: rhon_coef(10)    !< Coefficients for the intitial neutral density profile
  
  !> @name Numerical input profile for the neutral density (model 500)
  character(len=512)  :: rhon_file        !< ASCII file the neutral density profile is read from.
  logical             :: num_rhon         !< is set true if rho_file /= 'none'
  integer             :: num_rhon_len     !< Number of points in profile
  real*8, allocatable :: num_rhon_x(:)    !< Radial positions of profile points (PsiN values)
  real*8, allocatable :: num_rhon_y0(:)   !< Values of neutral density profile
  real*8, allocatable :: num_rhon_y1(:)   !< First derivatives of neutral density profile (\f$ d\rhon/d\Psi_N \f$)
  real*8, allocatable :: num_rhon_y2(:)   !< Second derivatives of neutral density profile (\f$ d^2\rhon/d\Psi_N^2 \f$)
  real*8, allocatable :: num_rhon_y3(:)   !< Third derivatives of neutral density profile (\f$ d^3\rhon/d\Psi_N^3 \f$)
  
  !> @name Numerical input profile for Fprofile
  character(len=512)  :: Fprofile_file      !< ASCII file the Fprofile is read from.
  logical             :: num_Fprofile       !< is set true if Fprofile_file /= 'none'
  integer             :: num_Fprofile_len   !< Number of points in profile
  real*8, allocatable :: num_Fprofile_x(:)  !< Radial positions of profile points (PsiN values)
  real*8, allocatable :: num_Fprofile_y0(:) !< Values of FFprime profile
  real*8, allocatable :: num_Fprofile_y1(:) !< First derivatives of Fprofile profile (\f$ dF/d\Psi_N \f$)
  real*8, allocatable :: num_Fprofile_y2(:) !< Second derivatives of Fprofile profile (\f$ d^2F/d\Psi_N^2 \f$)
  real*8, allocatable :: num_Fprofile_y3(:) !< Third derivatives of Fprofile profile (\f$ d^2F/d\Psi_N^2 \f$)

  !> @name Numerical input profile for Fprofile
  integer, parameter  :: n_Fprofile_internal_max = 300                !< INTERNAL Max Size of F-profile
  integer             :: n_Fprofile_internal                          !< INTERNAL Size of F-profile
  real*8              :: Fprofile_internal   (n_Fprofile_internal_max)!< INTERNAL F-profile, from  FFprime integration
  real*8              :: Fprofile_internal_d1(n_Fprofile_internal_max)!< INTERNAL F-profile, from  FFprime integration (first derivative)
  real*8              :: Fprofile_internal_d2(n_Fprofile_internal_max)!< INTERNAL F-profile, from  FFprime integration (second derivative)
  real*8              :: Fprofile_internal_d3(n_Fprofile_internal_max)!< INTERNAL F-profile, from  FFprime integration (third derivative)
  real*8              :: Fprofile_psi_max                             !< INTERNAL max psi_norm of F-profile
  real*8              :: Fprofile_tolerance                           !< INTERNAL tolerance (in %) for accuracy of F-profile compared to input FFprime

  !> @name Analytical input profile for FFprime
  real*8  :: FF_0              !< FF' value in the plasma center
  real*8  :: FF_1              !< FF' value in the SOL
  real*8  :: FF_coef(10)       !< Coefficients for FF' profile
  
  !> @name Numerical input profile for FFprime
  character(len=512)  :: ffprime_file      !< ASCII file the FF' profile is read from.
  logical             :: num_ffprime       !< is set true if ffprime_file /= 'none'
  integer             :: num_ffprime_len   !< Number of points in profile
  real*8, allocatable :: num_ffprime_x(:)  !< Radial positions of profile points (PsiN values)
  real*8, allocatable :: num_ffprime_y0(:) !< Values of FFprime profile
  real*8, allocatable :: num_ffprime_y1(:) !< First derivatives of FFprime profile (\f$ dFF'/d\Psi_N \f$)
  real*8, allocatable :: num_ffprime_y2(:) !< Second derivatives of FFprime profile (\f$ d^2FF'/d\Psi_N^2 \f$)

  !> --- Numerical input profiles for neoclassical coefficients
  logical             :: NEO              !< If .true. neoclassical effects are considered, (see [[neo|here]])
  character(len=512)  :: neo_file         !< ASCII file the aki and amu profiles is read from.
  logical             :: num_neo_file     !< automatically set true if neo_file /= 'none'
  integer             :: num_neo_len      !< Number of points in aki_neo, mu_neo profiles
  real*8, allocatable :: num_neo_psi(:)   !< Radial positions of profile points (PsiN values)
  real*8, allocatable :: num_aki_value(:) !< numerical aki profile (PsiN values)
  real*8, allocatable :: num_amu_value(:) !< numerical amu profile (PsiN values)
  real*8              :: aki_neo_const    !< if ( (NEO) .and. (neo_file=='none')), this constant value is used for aki_neo
  real*8              :: amu_neo_const    !< if ( (NEO) .and. (neo_file=='none')), this constant value is used for amu_neo

  !> @name RMP profiles
  logical :: output_bnd_elements !< If .true., writes bnd nodes and bnd elements in files 'boundary_nodes.dat' and 'boundary_elements.dat'
  logical :: RMP_on              !< Activates RMPs on boundary if .true. (the old version without STARWALL)
  character(len=512)  :: RMP_psi_cos_file  !< ASCII file the profiles of psi_RMP_cos and derivatives are read from
  character(len=512)  :: RMP_psi_sin_file  !< ASCII file the profiles of psi_RMP_sin and derivatives are read from
  real*8  :: RMP_growth_rate, RMP_ramp_up_time  !< parameters for time dependence of psi_RMP: Sigmoid f(t)= 1/ (1 + exp(-RMP_growth_rate*(t-RMP_ramp_up_time/2)))
  real*8  :: RMP_start_time    !< time when RMP coils are activated (RMP_on = .t.)
  real*8, allocatable :: psi_RMP_cos(:)
  real*8, allocatable :: dpsi_RMP_cos_dR(:)
  real*8, allocatable :: dpsi_RMP_cos_dZ(:)
  real*8, allocatable :: psi_RMP_sin(:)
  real*8, allocatable :: dpsi_RMP_sin_dR(:)
  real*8, allocatable :: dpsi_RMP_sin_dZ(:)
  integer             :: RMP_har_cos,RMP_har_sin ! Harmonics numbers for RMP-cos and RMP-sin(for ex. ntor=3, nperiod=2,RMP_har_cos=2, RMP_har_sin=3)
  integer, parameter  :: N_RMP_max = 10                  ! Maximum of RMP harmonics to take into account
  integer             :: Number_RMP_harmonics            ! Number_RMP_harmonics < N_RMP_max. If only one harmonic,  Number_RMP_harmonics=1, by default it's =1 in models/preset_parameters.f90 
  integer             :: RMP_har_cos_spectrum(N_RMP_max) = 0 ! If only one harmonic,by default RMP_har_cos_spectrum(1)=RMP_har_cos; 
  integer             :: RMP_har_sin_spectrum(N_RMP_max) = 0 ! If only one harmonic,by default RMP_har_sin_spectrum(1)=RMP_har_sin;


  !> @name toroidal rotation profile
  real*8              :: V_0               !< analytical parallel rotation profile -- central value
  real*8              :: V_1               !< analytical parallel rotation profile -- SOL value
  real*8              :: V_coef(10) = 0.d0 !< analytical parallel rotation profile -- coefficients
  character(len=512)  :: R_Z_psi_bnd_file  !< ASCII file for R_boundary,Z_boundary, psi_boundary, with n_boundary size.
  character(len=512)  :: wall_file         !< ASCII file for external wall geometry, if n_ext is greater than zero.
  
  !> @name Numerical input profile for the toroidal rotation
  character(len=512)  :: rot_file        !< ASCII file the parallel rotation profile is read from (see normalized_velocity_profile)
  logical             :: num_rot         !< automatically set true if rot_file /= 'none'
  integer             :: num_rot_len     !< Number of points in rotation profile
  real*8, allocatable :: num_rot_x(:)    !< Radial positions of profile points (PsiN values)
  real*8, allocatable :: num_rot_y0(:)   !< Values of toroidal rotation profile
  real*8, allocatable :: num_rot_y1(:)   !< First derivatives of toroidal rotation profile with respect to $\Psi_{N}$
  real*8, allocatable :: num_rot_y2(:)   !< Second derivatives of toroidal rotation profile with respect to $\Psi_{N}$
  real*8, allocatable :: num_rot_y3(:)   !< Third derivatives of toroidal rotation profile with respect to $\Psi_{N}$
  logical             :: normalized_velocity_profile !< if true, reads the normalized velocity profile as flux function, else Omega_tor is read as flux function. 
  
  !> @name Global quantities determined in each time step
  real*8, allocatable :: R_axis_t(:), Z_axis_t(:), psi_axis_t(:), R_xpoint_t(:,:), Z_xpoint_t(:,:),           &
    psi_xpoint_t(:,:), R_bnd_t(:), Z_bnd_t(:), psi_bnd_t(:),                                                  &
    current_t(:), beta_p_t(:), beta_t_t(:), beta_n_t(:), density_in_t(:), density_out_t(:), pressure_in_t(:), &
    pressure_out_t(:), heat_src_in_t(:), heat_src_out_t(:), part_src_in_t(:), part_src_out_t(:),   &
    E_tot_t(:), Helicity_tot_t(:), Kin_perp_tot_t(:), thermal_tot_t(:), kin_par_tot_t(:), ohmic_tot_t(:),      &
    Wmag_tot_t(:), Ip_tot_t(:), flux_Pvn_t(:), flux_qpar_t(:), dE_tot_dt(:), flux_qperp_t(:), flux_kinpar_t(:), &
    dWmag_tot_dt(:), dthermal_tot_dt(:), dkinpar_tot_dt(:), dkinperp_tot_dt(:), friction_dissip_tot_t(:), &
    Magwork_tot_t(:), thmwork_tot_t(:), viscopar_dissip_tot_t(:), viscopar_flux_t(:), li3_t(:),      &
    li3_tot_t(:), part_src_tot_t(:), heat_src_tot_t(:), volume_t(:), area_t(:), mag_ener_src_tot(:), &
    dpart_tot_dt(:), part_flux_Dpar_t(:), part_flux_Dperp_t(:), part_flux_vpar_t(:), part_flux_vperp_t(:), & 
    dnpart_tot_dt(:), npart_tot_t(:), npart_flux_t(:), density_tot_t(:), flux_poynting_t(:),         &
    thermal_e_tot_t(:), thermal_i_tot_t(:)

  !> @name gmres parameters
  integer             :: iter_precon        !< whenever the number of gmres iterations exceeds iter_precon, the preconditioning matrix is updated
  integer             :: max_steps_noUpdate !< whenever the steps without preconditioning matrix update exceeds max_steps_noUpdate, the preconditioning matrix is updated
  integer             :: gmres_m            !< gmres restart parameter (dimension)
  real*8              :: gmres_4            !< see gmres manual (error ratio between preconditioned and non-preconditioned error)
  real*8              :: gmres_tol          !< the tolerance for the gmres iterations to be seen as converged

  !> @name Taylor-Galerkin Stabilisation coefficients
  real*8              :: tgnum(n_var)   !< Coefficients for Taylor Galerkin stabilization for each equation separately
  real*8              :: tgnum_psi      !< Same as previous line, but avoiding equation indexing for model families 
  real*8              :: tgnum_u      
  real*8              :: tgnum_zj     
  real*8              :: tgnum_w      
  real*8              :: tgnum_rho    
  real*8              :: tgnum_T      
  real*8              :: tgnum_Ti     
  real*8              :: tgnum_Te     
  real*8              :: tgnum_vpar   
  real*8              :: tgnum_rhon   
  real*8              :: tgnum_nre    
  real*8              :: tgnum_AR     
  real*8              :: tgnum_AZ    
  real*8              :: tgnum_A3    

  !> @name Flag to determine whether or not we keep current source term  
  logical             :: keep_current_prof !< Artificial current source to approximately keep the initial current profile, i.e., \f$\eta(j-j0)\f$?
  
  !> @name Numerical parameters
  real*8              :: D_prof_neg         !< Particle diffusion coefficient in regions with negative density
  real*8              :: D_prof_neg_thresh  !< D_prof_neg becomes effective if rho < D_prof_neg_thresh
  real*8              :: ZK_prof_neg        !< Perp. heat diffusion coefficient in regions with negative temperature
  real*8              :: ZK_par_neg         !< Parallel diffusion coefficient in regions with negative temperature
  real*8              :: ZK_prof_neg_thresh !< ZK_prof_neg becomes effective if T < ZK_prof_neg_thresh
  real*8              :: ZK_par_neg_thresh  !< ZK_par_neg becomes effective if T < ZK_par_neg_thresh
  real*8              :: ZK_e_prof_neg        !< Perp. heat diffusion coefficient in regions with negative temperature
  real*8              :: ZK_e_par_neg         !< Parallel diffusion coefficient in regions with negative temperature
  real*8              :: ZK_e_prof_neg_thresh !< ZK_e_prof_neg becomes effective if T < ZK_e_prof_neg_thresh
  real*8              :: ZK_e_par_neg_thresh  !< ZK_e_par_neg becomes effective if T < ZK_e_par_neg_thresh
  real*8              :: ZK_i_prof_neg        !< Perp. heat diffusion coefficient in regions with negative temperature
  real*8              :: ZK_i_par_neg         !< Parallel diffusion coefficient in regions with negative temperature
  real*8              :: ZK_i_prof_neg_thresh !< ZK_i_prof_neg becomes effective if T < ZK_i_prof_neg_thresh
  real*8              :: ZK_i_par_neg_thresh  !< ZK_i_par_neg becomes effective if T < ZK_i_par_neg_thresh
  real*8              :: D_imp_extra_R           !< Additional impurity diffusivity in R-direction
  real*8              :: D_imp_extra_Z           !< Additional impurity diffusivity in Z-direction
  real*8              :: D_imp_extra_p           !< Additional impurity diffusivity in phi-direction
  real*8              :: D_imp_extra_neg         !< Additional impurity diffusion coefficient in regions with negative impurity density
  real*8              :: D_imp_extra_neg_thresh  !< D_imp_extra_neg becomes effective if rho_imp < D_imp_extra_neg_thresh
  real*8              :: T_min              !< minimum temperature (limits on the temperature dependence of resistivity etc.)
  real*8              :: rho_min            !< minimum density

  real*8              :: ne_SI_min          !< minimum e density (in SI unit) below which we cut-off the radiation loss
  real*8              :: Te_eV_min          !< minimum temperature (in eV) below which we cut-off the radiation loss
  real*8              :: rn0_min            !< minimum impurity density (in JU) for radiation loss cut-off

  integer             :: n_tor_fft_thresh   !< If n_tor >= n_tor_fft_thresh, element_matrix_fft will be used
  integer*8           :: fftw_plan          !< Required for FFTW library
  real*8              :: corr_neg_temp_coef(2) !< Parameters used in models/corr_neg.f90
  real*8              :: corr_neg_dens_coef(2) !< Parameters used in models/corr_neg.f90

  !> @name ECCD current sources
  real*8  :: jecamp             ! parameter, not to be confused with jec_source in element_matrix.f90
  real*8  :: jec_pos1, jec_pos2, jec_pos3, jec_pos4
  real*8  :: jec_width, jec_width2
  real*8  :: nu_jec_fast         ! 1/collision frequency
  real*8  :: nu_jec1_fast,nu_jec2_fast         ! 1/collision frequency
  real*8  :: mod_jec            ! extra parameters for ECCD
  real*8  :: JJ_par             ! velocity of resonent electrons
  real*8  :: jw1,jw2,jw3        ! parameters to determine current source

  !> @name Flag for thermalization term
  logical             :: thermalization ! If true turns on the ion-electron thermalization term

  !> @name (Currently unused)
  real*8  :: zjz_0, zjz_1,  zj_coef(10)
  real*8  :: D_neutral

  !> @name Particles-related input parameters
  logical :: restart_particles
  logical :: use_ncs          ! use neutral particles
  logical :: use_ccs          ! use current coupling scheme for fast particles
  logical :: use_pcs          ! use pressure coupling scheme for fast particles
  logical :: use_pcs_full     ! use full tensor pressure coupling scheme for fast particles
  logical :: use_cx           ! switch on sputtering         (in particle module)
  logical :: use_sputtering   ! switch on charge-exchange    (in particle module)
  logical :: use_ionisation   ! switch on ionisation         (in particle module)
  real*8  :: n_particles      ! the number of particles (real on purpose)
  real*8  :: tstep_particles  ! the time step for the particles
  integer :: nstep_particles  ! the number of particle time steps
  integer :: nsubstep_particles ! the number of particles substeps (without projection)
  real*8  :: filter_perp      ! particle projection smoothing parameter, poloidal plane
  real*8  :: filter_hyper     ! particle projection smoothing parameter, poloidal plane
  real*8  :: filter_par       ! particle projection smoothing parameter, parallel direction
  real*8  :: filter_perp_n0   ! particle projection smoothing parameter, poloidal plane (n=0)
  real*8  :: filter_hyper_n0  ! particle projection smoothing parameter, poloidal plane (n=0)
  real*8  :: filter_par_n0    ! particle projection smoothing parameter, parallel direction (n=0)
  
  !> @name Mode families preconditioner parameters
  integer, parameter :: n_fam_max = 100               !< maximum number of families
  integer :: n_mode_families                          !< number of families
  logical :: autodistribute_modes                     !< use automatic or manual mode distribution
  integer :: modes_per_family(n_fam_max)              !< Number of modes in families
  integer :: mode_families_modes(n_fam_max,n_fam_max) !< Mode numbers (i_tor) belonging to each family; first index: family number
  real*8  :: weights_per_family(n_fam_max)            !< Multiplication factor of family's contribution to the full solution
  logical :: autodistribute_ranks                     !< use automatic or manual rank distribution
  integer :: ranks_per_family(n_fam_max)              !< Number of MPI ranks per mode families
  
  contains
  
end module phys_module
