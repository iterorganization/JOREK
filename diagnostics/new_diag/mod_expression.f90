!> This module calculates "arbitrary" expressions at arbitrary positions for diagnostic purposes.
!!
!! Structure of the result(:,:,:,:) array:
!! - 1st index corresponds to toroidal positions
!! - 2nd and 3rd index correspond to positions in poloidal plane (2nd poloidal, 3rd radial)
!! - 4th index corresponds to the expressions
!!
!! To obtain good performance, evaluate several expressions at several positions in a single call.
!!
!! Note: Add a new expressions to init_expr() and eval_expr() consistently.
module mod_expression
  
  
  
  
  
  use constants
  use mod_parameters
  use mod_position
  use phys_module
  use diffusivities
  use corr_neg
  use mod_basisfunctions
  use mod_bootstrap_functions
  use mod_poloidal_currents
  use mod_impurity, only: radiation_function, radiation_function_linear
  use mod_atomic_coeff_deuterium, only : atomic_coeff_deuterium
  
  implicit none
  
  
  
  
  
  public
  private add
  
  
  
  
  
  ! --- Constants
  character(len=14), parameter, private :: THIS_MOD_NAME = 'mod_expression'
  integer,           parameter, private :: LEN_NAME      = 12
  integer,           parameter, private :: LEN_DESCR     = 54
  integer,           parameter, private :: LEN_DOMAIN    = 12
  integer,           parameter, private :: N_EXPR_MAX    = 2000
  !   --- Constants selecting a unit system for expression output
  integer,           parameter          :: JOREK_UNITS   = 0 !< Output expressions in JOREK units
  integer,           parameter          :: SI_UNITS      = 1 !< Output expressions in SI units
  
  
  
  
  
  ! --- Externally visible variables.
  !> Expressions evaluate to this value for positions outside the JOREK domain.
  !! This value may be changed by programs using the new_diag package!
  real*8, save :: expr_outside_value = 0.d0
  !> Coordinate expressions evaluate to this value for positions outside the JOREK domain.
  !! This value may be changed by programs using the new_diag package!
  real*8, save :: expr_outside_coord = 0.d0
  
  
  
  
  
  !> Datatype containing information on a single expression.
  type :: t_expr
    character(len=LEN_NAME)   :: name  !< Short name for the expression
    character(len=LEN_DESCR)  :: descr !< Brief explanation of the expression
    character(len=LEN_DOMAIN) :: domain!< Domain where this expression is defined
  end type t_expr
  
  ! > List of expressions.
  type :: t_expr_list
    integer                   :: n_expr  = 0      !< Number of expressions in this list
    integer                   :: n_coord = 0      !< Treat the first n_coord expr.s as coordinates
    type(t_expr)              :: expr(N_EXPR_MAX) !< The expressions
  end type t_expr_list
  
  
  
  
  
  ! --- Standard lists of expressions (require init_expr call first!).
  type(t_expr_list), save :: exprs_all, exprs_all_int, exprs_all_four     !< All available expressions.
  
  
  
  
  
  contains
  
  
  
  
  
  !> Initialization for the module. Should be called before any other 
  subroutine init_expr()
    
    exprs_all%n_expr     = 0
    exprs_all_int%n_expr = 0
    exprs_all_four%n_expr = 0
    call add(exprs_all, 'index_now   ', 'Restart file index (or number of run tsteps)          ')
    call add(exprs_all, 'R           ', 'Cylindrical Coordinate R (== Major Radius)            ')
    call add(exprs_all, 'Z           ', 'Cylindrical Coordinate Z                              ')
    call add(exprs_all, 'phi         ', 'Cylindrical Coordinate phi                            ')
    call add(exprs_all, 'theta       ', 'Poloidal Angle With Respect to Magnetic Axis          ')
    call add(exprs_all, 'theta_star  ', 'Poloidal Straight Field Line Angle (for flux surfaces)')
    call add(exprs_all, 'length      ', 'Length Along Poloidal Line (for poloidal lines)       ')
    call add(exprs_all, 'r_minor     ', 'Minor Radius From A = r_minor^2 pi (for flux surfaces)')
    call add(exprs_all, 'x           ', 'Cartesian Coordinate x                                ')
    call add(exprs_all, 'y           ', 'Cartesian Coordinate y                                ')
    call add(exprs_all, 'z           ', 'Cartesian Coordinate z (== Cylindrical Z)             ')
    call add(exprs_all, 'Psi_N       ', 'Normalized Poloidal Magnetic Flux                     ')
    call add(exprs_all, 'xjac        ', '2D Jacobian in the Poloidal Plane                     ')
    call add(exprs_all, 't           ', 'Simulation time                                       ')
    call add(exprs_all, 'Psi         ', 'Poloidal Magnetic Flux                                ')
    call add(exprs_all, 'dPsi_dt     ', 'Time derivative of poloidal magnetic flux             ')
    call add(exprs_all, 'u           ', 'Velocity Stream Function                              ')
    call add(exprs_all, 'Phi         ', 'Electric Potential Phi                                ')
    call add(exprs_all, 'zj          ', 'Toroidal Current Density Multiplied by 1/R            ')
    call add(exprs_all, 'currdens    ', 'Physical Toroidal Current Density (== zj/R)           ')
    call add(exprs_all, 'JR          ', 'Physical current density (R component)                ')
    call add(exprs_all, 'JZ          ', 'Physical current density (Z component)                ')
    call add(exprs_all, 'Jtor        ', 'Physical current density (phi component)              ')
    call add(exprs_all, 'Jpol        ', 'Poloidal current value in the poloidal field direction')
    call add(exprs_all, 'FFprime_loc ', 'Local FFprime value, calculated from 3D JxB=\grad p   ')
    call add(exprs_all, 'JxB_R       ', 'JxB force (R component)                               ')
    call add(exprs_all, 'JxB_Z       ', 'JxB force (Z component)                               ')
    call add(exprs_all, 'JxB_phi     ', 'JxB force (phi component)                             ')
    call add(exprs_all, 'gradP_R     ', 'Pressure gradient force (R component)                 ')
    call add(exprs_all, 'gradP_Z     ', 'Pressure gradient force (Z component)                 ')
    call add(exprs_all, 'gradP_phi   ', 'Pressure gradient force (phi component)               ')
    call add(exprs_all, 'gradP_B     ', 'Parallel pressure gradient (along B)                  ')
    call add(exprs_all, 'gradPdotCurv', 'grad p dot curvature                                  ')
    call add(exprs_all, 'curvat_R    ', 'curvature (= b . grad ( b )) in the R direction       ')
    call add(exprs_all, 'curvat_Z    ', 'curvature (= b . grad ( b )) in the Z direction       ')
    call add(exprs_all, 'curvat_phi  ', 'curvature (= b . grad ( b )) in the phi direction     ')
    call add(exprs_all, 'omega       ', 'Toroidal Vorticity Component                          ')
    call add(exprs_all, 'rho         ', 'Mass Density                                          ')
    call add(exprs_all, 'ne          ', 'Electron Density                                      ')
#ifdef WITH_Impurities
    call add(exprs_all, 'nimp        ', 'Impurity Density                                      ')
#endif
    call add(exprs_all, 'T           ', 'Temperature (Electrons plus Ions)                     ')
    call add(exprs_all, 'Te          ', 'Electron temperature (assuming Ti=Te)                 ')
    call add(exprs_all, 'vpar        ', 'Parallel Velocity (along magnetic field lines)        ')
    call add(exprs_all, 'eta_T       ', 'Temperature Dependent Resistivity                     ')
    call add(exprs_all, 'visco_T     ', 'Temperature Dependent Viscosity                       ')
    call add(exprs_all, 'zkpar_T     ', 'Temperature Dependent Parallel Heat Diffusivity       ')
    call add(exprs_all, 'dprof       ', 'Particle Diffusivity                                  ')
    call add(exprs_all, 'zkprof      ', 'Perpendicular Heat Diffusivity                        ')
    call add(exprs_all, 'pres        ', 'Total Pressure                                        ')
    call add(exprs_all, 'B_abs       ', 'Norm of the Magnetic Field Vector                     ')
    call add(exprs_all, 'Btor        ', 'Toroidal Magnetic Field Component                     ')
    call add(exprs_all, 'BR          ', 'Magnetic Field Component Along R                      ')
    call add(exprs_all, 'BZ          ', 'Vertical Magnetic Field Component                     ')
    call add(exprs_all, 'Btheta      ', 'Poloidal Magnetic Field Component                     ')
    call add(exprs_all, 'A_R         ', 'Radial Component of Magnetic Vector Potential         ')
    call add(exprs_all, 'A_Z         ', 'Vertical Component of Magnetic Vector Potential       ')
    call add(exprs_all, 'A_3         ', 'Tor. Comp. of Magn. Vect. Pot. Multipled by R (==Psi) ')
    call add(exprs_all, 'Er          ', 'Radial Electric Field                                 ')
    call add(exprs_all, 'Vtheta_i    ', 'Ion Poloidal Velocity                                 ')
    call add(exprs_all, 'Mach_par    ', 'Parallel Mach Number                                  ')
    call add(exprs_all, 'Mach_pol    ', 'Poloidal Mach Number                                  ')
    call add(exprs_all, 'V_sound     ', 'Sound Speed                                           ')
    call add(exprs_all, 'V_neo       ', 'Neoclassical Velocity                                 ')
    call add(exprs_all, 'Vperp_e     ', 'Electron Perpendicular Velocity                       ')
    call add(exprs_all, 'Vperp_i     ', 'Ion Perpendicular Velocity                            ')
    call add(exprs_all, 'V_ExB       ', 'ExB Velocity                                          ')
    call add(exprs_all, 'Vstar_e     ', 'Electron Diamagnetic Velocity                         ')
    call add(exprs_all, 'Vstar_i     ', 'Ion Diamagnetic Velocity                              ')
    call add(exprs_all, 'ki_neo      ', 'Neoclassical Heat Diffusivity                         ')
    call add(exprs_all, 'mu_neo      ', 'Neoclassical Friction Coefficient                     ')
    call add(exprs_all, 'T_e         ', 'Electron temperature                                  ')
    call add(exprs_all, 'T_i         ', 'Ion temperature                                       ')
    call add(exprs_all, 'E_||        ', 'E_|| for RE acceleration                              ')
    call add(exprs_all, 'E_crit      ', 'E_crit for RE avalanching (Connor-Hastie)             ')
    call add(exprs_all, 'E_dreicer   ', 'Electrical field for Dreicer RE primary source        ')
    call add(exprs_all, 'theta_geo   ', 'Polar angle with respect to Rgeo, Zgeo                ')
    call add(exprs_all, 'bnd_normal_R', 'R component of unit vector pointing outside JOREKs bnd', 'boundary    ')
    call add(exprs_all, 'bnd_normal_Z', 'Z component of unit vector pointing outside JOREKs bnd', 'boundary    ')
    call add(exprs_all, 'Bnorm       ', 'Normal     magnetic field to the JOREKs boundary      ', 'boundary    ')
    call add(exprs_all, 'Btan        ', 'Tangential magnetic field to the JOREKs boundary      ', 'boundary    ')
    call add(exprs_all, 'Jnorm       ', 'Normal current density to the JOREKs boundary         ', 'boundary    ')
    call add(exprs_all, 'Jpar        ', 'Parallel current density to the magnetic field        ', 'boundary    ')
    call add(exprs_all, 'Jpar_ionsat ', 'Parallel ion saturation current density               ', 'boundary    ')
    call add(exprs_all, 'vpar_norm   ', 'Perpendicular velocity to the boundary (vpar contrib) ', 'boundary    ')
    call add(exprs_all, 'vu_norm     ', 'Perpendicular velocity to the boundary (u contrib)    ', 'boundary    ')
    call add(exprs_all, 'vtot_norm   ', 'Total perpendicular velocity to the JOREKs boundary   ', 'boundary    ')
    call add(exprs_all, 'heatF_sheath', 'Sheath theory heatflux (gamma_sh nT vpar\cdot n)      ', 'boundary    ')
    call add(exprs_all, 'heatF_par_cd', 'Conductive parallel heat flux (normal to the boundary)', 'boundary    ')
    call add(exprs_all, 'heatF_prp_cd', 'Conductive perpend  heat flux (normal to the boundary)', 'boundary    ')
    call add(exprs_all, 'heatF_tot_cd', 'Conductive total    heat flux (normal to the boundary)', 'boundary    ')
    call add(exprs_all, 'heatF_par_cv', 'Convective parallel heat flux (normal to the boundary)', 'boundary    ')
    call add(exprs_all, 'heatF_prp_cv', 'Convective perpend  heat flux (normal to the boundary)', 'boundary    ')
    call add(exprs_all, 'heatF_tot_cv', 'Convective total    heat flux (normal to the boundary)', 'boundary    ')
    call add(exprs_all, 'heatF_total ', 'Total heat flux (normal to the boundary)              ', 'boundary    ')
    call add(exprs_all, 'kinEn_F_perp', 'Perpend kinetic energy flux (normal to the boundary)  ', 'boundary    ')
    call add(exprs_all, 'kinEn_F_par ', 'Parall  kinetic energy flux (normal to the boundary)  ', 'boundary    ')
    call add(exprs_all, 'kinEn_F_tot ', 'Total   kinetic energy flux (normal to the boundary)  ', 'boundary    ')
    call add(exprs_all, 'partF_par_cd', 'Conductive parallel particle flux (normal to the bnd) ', 'boundary    ')
    call add(exprs_all, 'partF_prp_cd', 'Conductive perpend  particle flux (normal to the bnd) ', 'boundary    ')
    call add(exprs_all, 'partF_par_cv', 'Convective parallel particle flux (normal to the bnd) ', 'boundary    ')
    call add(exprs_all, 'partF_prp_cv', 'Convective perpend  particle flux (normal to the bnd) ', 'boundary    ')
    call add(exprs_all, 'partF_total ', 'Total particle flux (normal to the boundary)          ', 'boundary    ')
    call add(exprs_all, 'npartF_total', 'Total neutral particle flux (normal to the boundary)  ', 'boundary    ')
    call add(exprs_all, 'ExB_norm    ', 'EM energy flux, Poynting vector (normal to boundary)  ', 'boundary    ')
    call add(exprs_all, 'gradP_norm  ', 'Total pressure gradient normal to the boundary        ', 'boundary    ')
#if JOREK_MODEL >= 303
    call add(exprs_all, 'J_bootstrap ', 'Bootstrap Current                                     ')
#endif
#if (defined WITH_Neutrals) || (defined WITH_Impurities)
    call add(exprs_all, 'radiation   ', 'Radiation terms for bolometry diagnostic              ')
#endif
#if (defined WITH_Neutrals) && (!defined WITH_Impurities)
    call add(exprs_all, 'brem        ', 'Brem terms for bolometry diagnostic                   ')
#endif
    ! --- List of volume and boundary integrals
    call add(exprs_all_int, 'index_now   ', 'Restart file index (or number of run tsteps)          ')
    call add(exprs_all_int, 'psi_axis    ', 'psi at magnetic axis                                  ')
    call add(exprs_all_int, 'R_axis      ', 'R of magnetic axis                                    ')
    call add(exprs_all_int, 'Z_axis      ', 'Z of magnetic axis                                    ')
    call add(exprs_all_int, 'psi_bnd     ', 'psi at boundary point (defining LCFS)                 ')
    call add(exprs_all_int, 'R_bnd       ', 'R of boundary point                                   ')
    call add(exprs_all_int, 'Z_bnd       ', 'Z of boundary point                                   ')
    call add(exprs_all_int, 'E_tot       ', 'Total energy                                          ')
    call add(exprs_all_int, 'E_in        ', 'Energy (inside  LCFS)                                 ')
    call add(exprs_all_int, 'E_out       ', 'Energy (outside LCFS)                                 ')
    call add(exprs_all_int, 'Wmag_tot    ', 'Total magnetic energy                                 ')
    call add(exprs_all_int, 'Wmag_in     ', 'Magnetic energy (inside  LCFS)                        ')
    call add(exprs_all_int, 'Wmag_out    ', 'Magnetic energy (outside LCFS)                        ')
    call add(exprs_all_int, 'Thermal_tot ', 'Total thermal energy                                  ')
    call add(exprs_all_int, 'Thermal_in  ', 'Thermal energy (inside  LCFS)                         ')
    call add(exprs_all_int, 'Thermal_out ', 'Thermal energy (outside LCFS)                         ')
    call add(exprs_all_int, 'Thermal_e_tot','Total electron thermal energy                                  ')
    call add(exprs_all_int, 'Thermal_e_in ','Thermal electron energy (inside  LCFS)                         ')
    call add(exprs_all_int, 'Thermal_e_out','Thermal electron energy (outside LCFS)                         ')
    call add(exprs_all_int, 'Thermal_i_tot','Total ion thermal energy                                  ')
    call add(exprs_all_int, 'Thermal_i_in ','Thermal ion energy (inside  LCFS)                         ')
    call add(exprs_all_int, 'Thermal_i_out','Thermal ion energy (outside LCFS)                         ')
    call add(exprs_all_int, 'Kin_par_tot ', 'Total parallel kinetic energy                         ')
    call add(exprs_all_int, 'Kin_par_in  ', 'Parallel kinetic energy (inside  LCFS)                ')
    call add(exprs_all_int, 'Kin_par_out ', 'Parallel kinetic energy (outside LCFS)                ')
    call add(exprs_all_int, 'Kin_perp_tot', 'Total perpendicular kinetic energy                    ')
    call add(exprs_all_int, 'Kin_perp_in ', 'Perpendicular kinetic energy (inside  LCFS)           ')
    call add(exprs_all_int, 'Kin_perp_out', 'Perpendicular kinetic energy (outside LCFS)           ')
    call add(exprs_all_int, 'Part_tot    ', 'Total number of ions                                  ')
    call add(exprs_all_int, 'Part_in     ', 'Number of ions  (inside  LCFS)                        ')
    call add(exprs_all_int, 'Part_out    ', 'Number of ions  (outside LCFS)                        ')
    call add(exprs_all_int, 'NPart_tot   ', 'Total number of neutral particles                     ')
    call add(exprs_all_int, 'Helicity_tot', 'Total magnetic helicity                               ')
    call add(exprs_all_int, 'Mag_work_tot', 'Total magnetic work = -\int v\cdot(JxB) dV            ')
    call add(exprs_all_int, 'Thm_work_tot', 'Total thermal work  = \int vpar\cdot\nabla p dV       ')
    call add(exprs_all_int, 'Part_src_tot', 'Total particle source                                 ')
    call add(exprs_all_int, 'Part_src_in ', 'Particle source (inside  LCFS)                        ')
    call add(exprs_all_int, 'Part_src_out', 'Particle source (outside LCFS)                        ')
    call add(exprs_all_int, 'Heat_src_tot', 'Total heat source                                     ')
    call add(exprs_all_int, 'Heat_src_in ', 'Heat source (inside  LCFS)                            ')
    call add(exprs_all_int, 'Heat_src_out', 'Heat source (outside LCFS)                            ')
    call add(exprs_all_int, 'Viscpar_diss', 'Total parallel viscosity dissipation                  ')
    call add(exprs_all_int, 'Friction_diss','Total frictional dissipation                          ')
    call add(exprs_all_int, 'Wmag_src_tot', 'Total magnetic energy source (from current source)    ')
    call add(exprs_all_int, 'Ohmic_tot   ', 'Total ohmic heating                                   ')
    call add(exprs_all_int, 'Ohmic_in    ', 'Ohmic heating (inside  LCFS)                          ')
    call add(exprs_all_int, 'Ohmic_out   ', 'Ohmic heating (outside LCFS)                          ')
    call add(exprs_all_int, 'P_vn        ', 'Boundary flux of outgoing pressure                    ')
    call add(exprs_all_int, 'qn_par      ', 'Boundary flux of the parallel thermal conduction      ')
    call add(exprs_all_int, 'qn_perp     ', 'Boundary flux of the perpendicular thermal conduction ')
    call add(exprs_all_int, 'sheath_heat ', 'Total heat flux (imposed by sheath BCs)               ')
    call add(exprs_all_int, 'kinpar_flux ', 'Boundary flux of parallel kinetic energy              ')
    call add(exprs_all_int, 'Poynting_flx', 'Boundary flux of magnetic energy (Poynting flux)      ')
    call add(exprs_all_int, 'vispar_flux ', 'Boundary flux of parallel viscous energy              ')
    call add(exprs_all_int, 'Dpar_pt_flx ', 'Boundary flux of parallel particle diffusion          ')
    call add(exprs_all_int, 'Dperp_pt_flx', 'Boundary flux of perpendicular particle diffusion     ')
    call add(exprs_all_int, 'vpar_pt_flx ', 'Boundary flux of parallel particle convection         ')
    call add(exprs_all_int, 'vperp_pt_flx', 'Boundary flux of perpendicular particle convection    ')
    call add(exprs_all_int, 'neut_pt_flx ', 'Boundary flux of neutral particles                    ')
    call add(exprs_all_int, 'Ip_tot      ', 'Total toroidal plasma current                         ')
    call add(exprs_all_int, 'Ip_in       ', 'Toroidal plasma current (inside  LCFS)                ')
    call add(exprs_all_int, 'Ip_out      ', 'Toroidal plasma current (outside LCFS)                ')
    call add(exprs_all_int, 'li3         ', 'Internal inductance inside LCFS, li(3)                ')
    call add(exprs_all_int, 'li3_tot     ', 'Internal inductance inside grid, li(3)                ')
    call add(exprs_all_int, 'beta_p      ', 'Poloidal beta, of the plasma inside LCFS              ')
    call add(exprs_all_int, 'beta_t      ', 'Toroidal beta, of the plasma inside LCFS              ')
    call add(exprs_all_int, 'beta_n      ', 'Normalized beta, of the plasma inside LCFS            ')
    call add(exprs_all_int, 'area        ', 'Poloidal cross section area inside LCFS               ')
    call add(exprs_all_int, 'volume      ', 'Plasma volume, inside LCFS                            ')
    call add(exprs_all_int, 'q02         ', 'Safety factor at psin=0.02                            ')
    call add(exprs_all_int, 'q95         ', 'Safety factor at psin=0.95                            ')
    call add(exprs_all_int, 'q99         ', 'Safety factor at psin=0.99                            ')
    call add(exprs_all_int, 'I_halo      ', 'Total poloidal halo currents                          ')
    call add(exprs_all_int, 'TPF_halo    ', 'Toroidal peaking factor of the poloidal halos         ')

    call add(exprs_all_four, 'absolute    ', 'Absolute value of 2D Fourier analysis                 ')
    call add(exprs_all_four, 'real        ', 'Real part      of 2D Fourier analysis                 ')
    call add(exprs_all_four, 'imaginary   ', 'Imaginary part of 2D Fourier analysis                 ')
    call add(exprs_all_four, 'phase       ', 'Complex phase  of 2D Fourier analysis                 ')

  end subroutine init_expr
  
  
  
  
  
  !> [Private] Auxilliary routine for init_expr.
  subroutine add(expr_list, name, descr, domain)
    
    ! --- Routine parameters
    type(t_expr_list),                   intent(inout) :: expr_list
    character(len=LEN_NAME),             intent(in)    :: name
    character(len=LEN_DESCR),            intent(in)    :: descr
    character(len=LEN_DOMAIN), optional, intent(in)    :: domain
    
    expr_list%n_expr = expr_list%n_expr + 1
    expr_list%expr(expr_list%n_expr)%name     = name
    expr_list%expr(expr_list%n_expr)%descr    = descr
    expr_list%expr(expr_list%n_expr)%domain   = 'all'
    if ( present(domain) ) &
      expr_list%expr(expr_list%n_expr)%domain = domain
    
  end subroutine add
  
  
  !> Creates a subset of all available expressions.
  function exprs(name, n_expr, n_coord, exprs_all_local) result(expr_list)
    type(t_expr_list) :: expr_list
    
    character(len=64), parameter :: THIS_ROUTINE_NAME = trim(THIS_MOD_NAME) // ':exprs'
    
    ! --- Routine parameters
    character(len=*),            intent(in) :: name(n_expr)
    integer,                     intent(in) :: n_expr
    integer, optional,           intent(in) :: n_coord
    type(t_expr_list), optional, intent(in) :: exprs_all_local
    type(t_expr_list)                       :: exprs_all_local0
    
    ! --- Local variables
    integer :: i, j, k
    
    exprs_all_local0 = exprs_all
    if ( present(exprs_all_local) ) exprs_all_local0 = exprs_all_local

    k = 0
    do i = 1, n_expr
      j = get_expr_num(exprs_all_local0, trim(name(i)))
      if ( j < 1 ) then
        write(*,*) 'WARNING in '//trim(THIS_ROUTINE_NAME)//': Unknown expression "'//trim(name(i)) &
          //'" ignored.'
        cycle
      end if
      k = k + 1
      expr_list%expr(k) = exprs_all_local0%expr(j)
    end do
    expr_list%n_expr = k
    
    expr_list%n_coord = 0
    if ( present(n_coord) ) expr_list%n_coord = n_coord
    
  end function exprs

  
  
  !> Merge several expression lists.    ### NOT USED AT PRESENT
  recursive function join_exprs(list1, list2, list3, list4, list5, list6, list7, list8, list9)     &
    result(expr_list)
    type(t_expr_list) :: expr_list
    
    ! --- Routine parameters
    type(t_expr_list),           intent(in) :: list1, list2
    type(t_expr_list), optional, intent(in) :: list3, list4, list5, list6, list7, list8, list9
    
    integer :: i, j
    logical :: found
    
    expr_list = list1
    
    do i = 1, list2%n_expr
      
      ! --- Already in the list? (avoid duplicates)
      found = .false.
      do j = 1, expr_list%n_expr
        if ( expr_list%expr(j)%name == list2%expr(i)%name ) then
          found = .true.
          exit
        end if
      end do
      
      ! --- Add to expr_list
      if ( .not. found ) then
        expr_list%n_expr = expr_list%n_expr + 1
        expr_list%expr(expr_list%n_expr) = list2%expr(i)
      end if
      
    end do
    
    ! --- Recursive joins of the remaining lists
    if ( present(list3) ) expr_list = join_exprs(expr_list, list3)
    if ( present(list4) ) expr_list = join_exprs(expr_list, list4)
    if ( present(list5) ) expr_list = join_exprs(expr_list, list5)
    if ( present(list6) ) expr_list = join_exprs(expr_list, list6)
    if ( present(list7) ) expr_list = join_exprs(expr_list, list7)
    if ( present(list8) ) expr_list = join_exprs(expr_list, list8)
    if ( present(list9) ) expr_list = join_exprs(expr_list, list9)
    
  end function join_exprs
  
  
  
  
  
  !> Prints an expression list in a readable table format.
  subroutine print_exprs(expr_list, short)
    
    ! --- Routine parameters
    type(t_expr_list),  intent(in) :: expr_list
    logical, optional,  intent(in) :: short
    
    ! --- Local variables
    integer :: i
    
    if ( present(short) .and. (short) ) then
      
      905 format(1x,a,',')
      do i = 1, expr_list%n_expr
        write(*,905,advance='no') trim(expr_list%expr(i)%name)
      end do
      write(*,*)
      
    else
      
      900 format(1x,a)
      901 format(1x,i3.3,' | ',a,' | ',a,' | ',a)
      902 format(1x,85('-'))
      
      write(*,*)
      write(*,*) 'List of Diagnostic Expressions:'
      write(*,*)
      
      write(*,902)
      write(*,900) 'Num | Name         | Description                                            | Domain'
      write(*,902)
      do i = 1, expr_list%n_expr
        write(*,901) i, expr_list%expr(i)%name, expr_list%expr(i)%descr, trim(expr_list%expr(i)%domain)
      end do
      write(*,902)
      write(*,*)
      
    end if
    
  end subroutine print_exprs
  
  
  
  
  !> Find out expression number in an expression list.
  integer function get_expr_num(expr_list, name) result(num)
    
    ! --- Routine parameters
    type(t_expr_list),       intent(in) :: expr_list
    character(len=*), intent(in) :: name
    
    ! --- Local variables
    integer :: i
    
    num = -99
    do i = 1, exprs_all%n_expr
      if ( trim(expr_list%expr(i)%name) == trim(name) ) then
        num = i
        exit
      end if
    end do
    
  end function get_expr_num
  
  
  
  
  
  !> Check that an expression exists.
  logical function expr_exists(expr_list, name) result(exists)
    
    ! --- Routine parameters
    type(t_expr_list),       intent(in) :: expr_list
    character(len=LEN_NAME), intent(in) :: name
    
    exists = ( get_expr_num(expr_list, name) > 0 )
    
  end function expr_exists
  
  
  
  
  
  !> Get the description of an expression.
  character(len=LEN_DESCR) function expr_descr(name) result(descr)
    
    ! --- Routine parameters
    character(len=LEN_NAME) :: name
    
    integer :: num
    
    num = get_expr_num(exprs_all, name)
    
    if ( num > 0 ) then
      descr = exprs_all%expr(num)%descr
    else
      descr = '<UNKNOWN EXPRESSION>'
    end if
    
  end function expr_descr
  
  
  
  
  
  !> Evaluate one/several expressions at one/several poloidal and one/several toroidal positions.
  subroutine eval_expr(eq, units, expr_list, pol_pos_list, tor_pos_list, result, ierr)

    character(len=64), parameter :: THIS_ROUTINE_NAME = trim(THIS_MOD_NAME) // ':eval_expr'
    
    ! --- Routine parameters
    type(t_equil_state),          intent(in)    :: eq
    integer,                      intent(in)    :: units !< Output in which units?
    type(t_expr_list),            intent(in)    :: expr_list
    type(t_pol_pos_list), target, intent(in)    :: pol_pos_list
    type(t_tor_pos_list), target, intent(in)    :: tor_pos_list
    real*8, allocatable,          intent(inout) :: result(:,:,:,:)
    integer,                      intent(out)   :: ierr
    
    ! --- Local variables
    type(t_pol_pos), pointer :: pol_pos
    type(t_tor_pos), pointer :: tor_pos
    type(type_element)       :: element
    type(type_node)          :: nodes(n_vertex_max)
    integer :: ipolpos, jpolpos, itorpos, iexpr, ielm, i, j, k, i_tor
    real*8  :: xjac, xjac_R, xjac_Z, R, R_s, R_t, R_st, R_ss, R_tt, Z, Z_s, Z_t, Z_st, Z_ss, Z_tt, &
      s, t, H(n_vertex_max,n_order+1), H_s(n_vertex_max,n_order+1), H_t(n_vertex_max,n_order+1),   &
      H_st(n_vertex_max,n_order+1), H_ss(n_vertex_max,n_order+1), H_tt(n_vertex_max,n_order+1),    &
      HZ(n_tor), HZ_p(n_tor), HZ_pp(n_tor), phi, res, BigR, BigR_R, x_cart, y_cart, theta
    real*8  :: ps0, ps0_s, ps0_t, ps0_ss, ps0_tt, ps0_st, ps0_p, ps0_pp, u0, u0_s, u0_t, u0_ss,    &
      u0_tt, u0_st, u0_p, u0_pp, zj0, zj0_s, zj0_t, zj0_ss, zj0_tt, zj0_st, zj0_p, zj0_pp, w0,     &
      w0_s, w0_t, w0_ss, w0_tt, w0_st, w0_p, w0_pp, r0, r0_s, r0_t, r0_ss, r0_tt, r0_st, r0_p,     &
      r0_pp, T0, T0_s, T0_t, T0_ss, T0_tt, T0_st, T0_p, T0_pp, Vpar0, Vpar0_s, Vpar0_t, Vpar0_ss,  &
      Vpar0_tt, Vpar0_st, Vpar0_p, Vpar0_pp, psi_norm, AR0, AR0_p, AR0_s, AR0_t, AR0_sp, AR0_tp,   &
      AR0_ss, AR0_tt, AR0_st, AR0_pp, AZ0, AZ0_p, AZ0_s, AZ0_t, AZ0_sp, AZ0_tp, AZ0_ss, AZ0_tt,    &
      AZ0_st, AZ0_pp, A30, A30_p, A30_s, A30_t, A30_ss, A30_tt, A30_st, A30_sp, A30_tp, A30_pp,    &
      Fprofile, Fprofile_s, Fprofile_t, Fprofile_R, Fprofile_Z
    real*8  :: ps0_R, ps0_Z, ps0_RR, ps0_ZZ, ps0_RZ, u0_R, u0_Z, u0_RR, u0_ZZ, u0_RZ, vv2, zj0_R,  &
      zj0_Z, zj0_RR, zj0_ZZ, zj0_RZ, w0_R, w0_Z, w0_RR, w0_ZZ, w0_RZ, r0_R, r0_Z, r0_RR, r0_ZZ,    &
      r0_RZ, r0_hat, r0_R_hat, r0_Z_hat, T0_R, T0_Z, T0_RR, T0_ZZ, T0_RZ, T0_ps0_R, T0_ps0_Z,      &
      Vpar0_R, Vpar0_Z, Vpar0_RR, Vpar0_ZZ, Vpar0_RZ, P0, P0_R, P0_Z, P0_s, P0_t, P0_p, P0_pp,     &
      P0_RR, P0_ZZ, P0_RZ, BB2, Btor, BR, BZ, BR_Z, BZ_R, Btheta, psi_abs, E_par, E_crit,          &
      E_dreicer, AR0_R, AR0_Z, AZ0_R, AZ0_Z, A30_R, A30_Z, AR0_Rp, AZ0_Zp, A30_RR, A30_ZZ, AR0_ZZ, &
      AR0_RR, AZ0_RR, AZ0_ZZ, A30_Rp, A30_Zp, AR0_RZ, AZ0_RZ, AR0_Zp, AZ0_Rp, A30_RZ, BR_p, BZ_p,  &
      BP_Z, BP_R, BR_R, BZ_Z, B_R, B_Z, Kappa_R, Kappa_Z, Kappa_phi
    real*8  :: eta_T, deta_dT, d2eta_d2T, visco_T, dvisco_dT, ZKpar_T, dZKpar_dT, D_prof, ZK_prof
    real*8 :: Ti0, Ti0_s, Ti0_t, Ti0_st, Ti0_ss, Ti0_tt, Ti0_p, Ti0_pp, Te0, Te0_s, Te0_t, Te0_st, &
      Te0_ss, Te0_tt, Te0_p, Te0_pp, Ti0_R, Ti0_Z, Te0_R, Te0_Z, Er, Vtheta, Mach_par, Mach_pol,   &
      Vsound, Vneo, Vperp_e, Vperp_i, V_ExB, Vstar_e, Vstar_i, mu_neo, ki_neo, J_boot, Te0_eV,     &
      ne0_20, ln_Lambda, ln_Lambda0, dpsi_dt 
    real*8 :: FFprime_loc, Jpol, JpolR, JpolZ, Btot, Jpar, Jpar_ionsat, fact_jsat, Bnorm, Btan, Jtor
    real*8 :: nmlR, nmlZ, theta_geo, VR, VZ, V_phi, Vpar_tot, VperpR, VperpZ
    real*8 :: hh, hh_s, hh_t, hh_ss, hh_tt, hh_st, hhz, hhz_p, hhz_pp, sz, vv(0:n_var)
    real*8 :: delta_g(n_var), delta_s(n_var), delta_t(n_var)
    ! --- Fluxes
    real*8  ::  ZKpar_flux, ZKperp_flux, Dpar_flux, Dperp_flux, partF_cnv_par, partF_cnv_tot
    real*8  ::  pres_flux_par, pres_flux_tot, kin_flux_par, kin_flux_tot, neut_part_flux, ExB_norm 
    ! --- Normalization factors
    real*8  :: rho_norm, fact_time, fact_mu_zero, fact_ne, fact_rho, fact_T, fact_vpar,            &
      fact_resistiv, fact_Er, fact_flux, fact_rad
    real*8  :: rn0, rn0_s, rn0_t, rn0_ss, rn0_tt, rn0_st, rn0_p, rn0_pp, rn0_R, rn0_Z

#if (defined WITH_Neutrals) || (defined WITH_Impurities)
    real*8  :: Te_corr_eV
    real*8  :: LradDrays_T, LradDcont_T, Sion_T, Srec_T
    real*8  :: dLradDrays_dT, dLradDcont_dT, dSion_dT, dSrec_dT
    real*8  :: ne_SI                              ! Electron density used in radiation rate
#endif
#if (defined WITH_Neutrals) && (!defined WITH_Impurities)
    real*8  :: Arad_bg, Brad_bg, Crad_bg, frad_bg
    real*8  :: Lrad_imp, r_imp, coef_rad_imp
#endif
#ifdef WITH_Impurities
    ! See https://www.jorek.eu/wiki/doku.php?id=model500_501_555 for details
    real*8  :: coef_rad_1, Te_eV
    real*8  :: T0_corr, r0_corr, rn0_corr
    ! Atomic physics coefficients:
    !   -Mass ratio between main ions and impurites (m_i/m_imp)
    real*8  :: m_i_over_m_imp
    !   -Mean impurity ionization state
    real*8  :: Z_imp, T0_Zimp, alpha_Zimp
    !   -Coefficients related to Z_imp
    real*8  :: alpha_imp
    real*8  :: beta_imp
    !   -Radiation from injected impurities
    real*8  :: Lrad                                ! Radiation rate
    real*8  :: A0_rad, A1_rad, T1_rad, sig1_rad    ! Radiation rate parameters
    real*8  :: A2_rad, T2_rad, sig2_rad
    !   -Temporary variable for charge state distribution
    real*8, allocatable :: P_imp(:)
    real*8  :: E_ion
    integer*8  :: ion_i, ion_k
#endif
    
    ierr = 0
    
    ! --- Some sanity checks
    if ( expr_list%n_expr < 1 ) then
      write(*,*) 'ERROR in '//trim(THIS_ROUTINE_NAME)//': expr_list%n_expr < 1 encountered.'
      ierr = -101
      return
    else if ( minval(pol_pos_list%n_pos) < 1 ) then
      write(*,*) 'ERROR in '//trim(THIS_ROUTINE_NAME)//': pol_pos_list%n_pos < 1 encountered.'
      ierr = -102
      return
    else if ( .not. allocated(pol_pos_list%pos) ) then
      write(*,*) 'ERROR in '//trim(THIS_ROUTINE_NAME)//': pol_pos_list%pos not allocated.'
      ierr = -103
      return
    else if ( tor_pos_list%n_pos < 1 ) then
      write(*,*) 'ERROR in '//trim(THIS_ROUTINE_NAME)//': tor_pos_list%n_pos < 1 encountered.'
      ierr = -102
      return
    else if ( .not. allocated(tor_pos_list%pos) ) then
      write(*,*) 'ERROR in '//trim(THIS_ROUTINE_NAME)//': tor_pos_list%pos not allocated.'
      ierr = -103
      return
    end if

#ifdef WITH_Impurities
     select case ( trim(imp_type) )
       case('D2')
         m_i_over_m_imp = central_mass/2.  ! Deuterium mass = 2 u
       case('Ar')
         m_i_over_m_imp = central_mass/40. ! Argon mass = 40 u 
       case('Ne')
         m_i_over_m_imp = central_mass/20. ! Neon mass = 20 u 
       case default
         write(*,*) 'ERROR: Unknown imp_type.'
         stop
     end select
#endif
    
    if ( allocated(result) ) deallocate(result)
    allocate( result(tor_pos_list%n_pos, pol_pos_list%n_pos(1), pol_pos_list%n_pos(2),             &
      expr_list%n_expr) )
    
    ! --- Loop over positions in the poloidal plane
    loop_pol1: do ipolpos = 1, pol_pos_list%n_pos(1)
      loop_pol2: do jpolpos = 1, pol_pos_list%n_pos(2)
        pol_pos => pol_pos_list%pos(ipolpos,jpolpos)
        
        if ( pol_pos%outside ) then
          result(:, ipolpos, jpolpos, :expr_list%n_coord  ) = expr_outside_coord
          result(:, ipolpos, jpolpos, expr_list%n_coord+1:) = expr_outside_value
          cycle
        end if
        
        ! --- Poloidal Coordinates and Basis Functions
        R    = pol_pos%R
        R_s  = pol_pos%R_s
        R_t  = pol_pos%R_t
        R_st = pol_pos%R_st
        R_ss = pol_pos%R_ss
        R_tt = pol_pos%R_tt
        Z    = pol_pos%Z
        Z_s  = pol_pos%Z_s
        Z_t  = pol_pos%Z_t
        Z_st = pol_pos%Z_st
        Z_ss = pol_pos%Z_ss
        Z_tt = pol_pos%Z_tt
        s    = pol_pos%s
        t    = pol_pos%t
        ielm = pol_pos%ielm
        nmlR = pol_pos%bnd_normal(1)
        nmlZ = pol_pos%bnd_normal(2)
        BigR   = R    ! Just two different names for R
        BigR_R = 1.d0 ! Trivial derivative
        call basisfunctions(s, t, H, H_s, H_t, H_st, H_ss, H_tt)
        
        ! --- Poloidal angle theta
        theta = atan2(Z-eq%Z_axis, R-eq%R_axis)
        if ( theta < 0.d0 ) theta = theta + 2.d0*PI

        ! --- Geometrical poloidal angle
        theta_geo = atan2(Z-Z_geo, R-R_geo)
        if ( theta_geo < 0.d0 ) theta_geo = theta_geo + 2.d0*PI
        
        ! --- 2D Jacobian
        xjac   = R_s * Z_t - R_t * Z_s
        xjac_R = ( R_ss * Z_t**2 - 2*R_st * Z_s*Z_t + R_tt * Z_s**2 + R_s * (Z_st*Z_t - Z_tt*Z_s )   &
                 + R_t * (Z_st*Z_s - Z_ss*Z_t ) ) / xjac
        xjac_Z = ( Z_ss * R_t**2 - 2*Z_st * R_s*R_t + Z_tt * R_s**2 + Z_s * (R_st*R_t - R_tt*R_s )   &
                 + Z_t * (R_st*R_s - R_ss*R_t ) ) / xjac
        
        ! --- Elements and nodes
        element  = pol_pos%element
        nodes(:) = pol_pos%nodes(:)
        
        ! --- Loop over toroidal positions
        loop_tor: do itorpos = 1, tor_pos_list%n_pos
          tor_pos => tor_pos_list%pos(itorpos)
          
          ! --- Toroidal Coordinate and Basis Functions
          phi     = tor_pos%phi
          HZ(1)   = 1.d0
          HZ_p(1) = 0.d0
          do i = 1, (n_tor-1) / 2
            HZ(2*i)      =                           cos(mode(2*i)  *phi)
            HZ_p(2*i)    = - float(mode(2*i))      * sin(mode(2*i)  *phi)
            HZ_pp(2*i)   = - float(mode(2*i))**2   * cos(mode(2*i)  *phi)
            HZ(2*i+1)    =                           sin(mode(2*i+1)*phi)
            HZ_p(2*i+1)  = + float(mode(2*i+1))    * cos(mode(2*i+1)*phi)
            HZ_pp(2*i+1) = - float(mode(2*i+1))**2 * sin(mode(2*i+1)*phi)
          end do
          
          ! --- Cartesian Coordinates
          x_cart = + R*cos(phi)
          y_cart = - R*sin(phi)
          
          ps0   = 0.d0; ps0_s   = 0.d0; ps0_t   = 0.d0; ps0_ss   = 0.d0; ps0_tt   = 0.d0; ps0_st   = 0.d0; ps0_p   = 0.d0; ps0_pp   = 0.d0
          u0    = 0.d0; u0_s    = 0.d0; u0_t    = 0.d0; u0_ss    = 0.d0; u0_tt    = 0.d0; u0_st    = 0.d0; u0_p    = 0.d0; u0_pp    = 0.d0
          zj0   = 0.d0; zj0_s   = 0.d0; zj0_t   = 0.d0; zj0_ss   = 0.d0; zj0_tt   = 0.d0; zj0_st   = 0.d0; zj0_p   = 0.d0; zj0_pp   = 0.d0
          w0    = 0.d0; w0_s    = 0.d0; w0_t    = 0.d0; w0_ss    = 0.d0; w0_tt    = 0.d0; w0_st    = 0.d0; w0_p    = 0.d0; w0_pp    = 0.d0
          r0    = 0.d0; r0_s    = 0.d0; r0_t    = 0.d0; r0_ss    = 0.d0; r0_tt    = 0.d0; r0_st    = 0.d0; r0_p    = 0.d0; r0_pp    = 0.d0
          T0    = 0.d0; T0_s    = 0.d0; T0_t    = 0.d0; T0_ss    = 0.d0; T0_tt    = 0.d0; T0_st    = 0.d0; T0_p    = 0.d0; T0_pp    = 0.d0
          Ti0   = 0.d0; Ti0_s   = 0.d0; Ti0_t   = 0.d0; Ti0_ss   = 0.d0; Ti0_tt   = 0.d0; Ti0_st   = 0.d0; Ti0_p   = 0.d0; Ti0_pp   = 0.d0
          Te0   = 0.d0; Te0_s   = 0.d0; Te0_t   = 0.d0; Te0_ss   = 0.d0; Te0_tt   = 0.d0; Te0_st   = 0.d0; Te0_p   = 0.d0; Te0_pp   = 0.d0
          Vpar0 = 0.d0; Vpar0_s = 0.d0; Vpar0_t = 0.d0; Vpar0_ss = 0.d0; Vpar0_tt = 0.d0; Vpar0_st = 0.d0; Vpar0_p = 0.d0; Vpar0_pp = 0.d0
          AR0   = 0.d0; AR0_s   = 0.d0; AR0_t   = 0.d0; AR0_ss   = 0.d0; AR0_tt   = 0.d0; AR0_st   = 0.d0; AR0_p   = 0.d0; AR0_pp   = 0.d0               
          AZ0   = 0.d0; AZ0_s   = 0.d0; AZ0_t   = 0.d0; AZ0_ss   = 0.d0; AZ0_tt   = 0.d0; AZ0_st   = 0.d0; AZ0_p   = 0.d0; AZ0_pp   = 0.d0
          A30   = 0.d0; A30_s   = 0.d0; A30_t   = 0.d0; A30_ss   = 0.d0; A30_tt   = 0.d0; A30_st   = 0.d0; A30_p   = 0.d0; A30_pp   = 0.d0
          rn0   = 0.d0; rn0_s   = 0.d0; rn0_t   = 0.d0; rn0_ss   = 0.d0; rn0_tt   = 0.d0; rn0_st   = 0.d0; rn0_p   = 0.d0; rn0_pp   = 0.d0

          ! Extra derivatives for current density calculation
          AR0_sp   = 0.d0; AR0_tp   = 0.d0
          AZ0_sp   = 0.d0; AZ0_tp   = 0.d0
          A30_sp   = 0.d0; A30_tp   = 0.d0

          delta_g(:) = 0.d0; delta_s(:) = 0.d0; delta_t(:) = 0.d0
          Fprofile = 0.d0;  Fprofile_s = 0.d0;  Fprofile_t = 0.d0
          
          ! --- Reconstruct variables
          do i = 1, n_vertex_max
            do j = 1, n_order+1
              
              sz    = element%size(i,j)
              hh    = H   (i,j)
              hh_s  = H_s (i,j)
              hh_t  = H_t (i,j)
              hh_ss = H_ss(i,j)
              hh_tt = H_tt(i,j)
              hh_st = H_st(i,j)
	      
#ifdef fullmhd
              Fprofile   = Fprofile   + nodes(i)%Fprof_eq(j) * sz * hh  
              Fprofile_s = Fprofile_s + nodes(i)%Fprof_eq(j) * sz * hh_s  
              Fprofile_t = Fprofile_t + nodes(i)%Fprof_eq(j) * sz * hh_t  
#endif
              
              do i_tor = 1, n_tor
                
                hhz    = HZ   (i_tor)
                hhz_p  = HZ_p (i_tor)
                hhz_pp = HZ_pp(i_tor)
                vv(:)  = 0.d0
                vv(1:n_var)  = nodes(i)%values(i_tor,j,:)
                
                ! --- Poloidal Flux
                ps0      = ps0      + vv(var_psi) * sz * hh    * hhz
                ps0_s    = ps0_s    + vv(var_psi) * sz * hh_s  * hhz
                ps0_t    = ps0_t    + vv(var_psi) * sz * hh_t  * hhz
                ps0_ss   = ps0_ss   + vv(var_psi) * sz * hh_ss * hhz
                ps0_tt   = ps0_tt   + vv(var_psi) * sz * hh_tt * hhz
                ps0_st   = ps0_st   + vv(var_psi) * sz * hh_st * hhz
                ps0_p    = ps0_p    + vv(var_psi) * sz * hh    * hhz_p
                ps0_pp   = ps0_pp   + vv(var_psi) * sz * hh    * hhz_pp
                
                ! --- Stream Function
                u0       = u0       + vv(var_u) * sz * hh    * hhz
                u0_s     = u0_s     + vv(var_u) * sz * hh_s  * hhz
                u0_t     = u0_t     + vv(var_u) * sz * hh_t  * hhz
                u0_ss    = u0_ss    + vv(var_u) * sz * hh_ss * hhz
                u0_tt    = u0_tt    + vv(var_u) * sz * hh_tt * hhz
                u0_st    = u0_st    + vv(var_u) * sz * hh_st * hhz
                u0_p     = u0_p     + vv(var_u) * sz * hh    * hhz_p
                u0_pp    = u0_pp    + vv(var_u) * sz * hh    * hhz_pp
                
                ! --- Current
                zj0      = zj0      + vv(var_zj) * sz * hh    * hhz
                zj0_s    = zj0_s    + vv(var_zj) * sz * hh_s  * hhz
                zj0_t    = zj0_t    + vv(var_zj) * sz * hh_t  * hhz
                zj0_ss   = zj0_ss   + vv(var_zj) * sz * hh_ss * hhz
                zj0_tt   = zj0_tt   + vv(var_zj) * sz * hh_tt * hhz
                zj0_st   = zj0_st   + vv(var_zj) * sz * hh_st * hhz
                zj0_p    = zj0_p    + vv(var_zj) * sz * hh    * hhz_p
                zj0_pp   = zj0_pp   + vv(var_zj) * sz * hh    * hhz_pp
                
                ! --- Vorticity
                w0       = w0       + vv(var_w) * sz * hh    * hhz
                w0_s     = w0_s     + vv(var_w) * sz * hh_s  * hhz
                w0_t     = w0_t     + vv(var_w) * sz * hh_t  * hhz
                w0_ss    = w0_ss    + vv(var_w) * sz * hh_ss * hhz
                w0_tt    = w0_tt    + vv(var_w) * sz * hh_tt * hhz
                w0_st    = w0_st    + vv(var_w) * sz * hh_st * hhz
                w0_p     = w0_p     + vv(var_w) * sz * hh    * hhz_p
                w0_pp    = w0_pp    + vv(var_w) * sz * hh    * hhz_pp
                
                ! --- Density
                r0       = r0       + vv(var_rho) * sz * hh    * hhz
                r0_s     = r0_s     + vv(var_rho) * sz * hh_s  * hhz
                r0_t     = r0_t     + vv(var_rho) * sz * hh_t  * hhz
                r0_ss    = r0_ss    + vv(var_rho) * sz * hh_ss * hhz
                r0_tt    = r0_tt    + vv(var_rho) * sz * hh_tt * hhz
                r0_st    = r0_st    + vv(var_rho) * sz * hh_st * hhz
                r0_p     = r0_p     + vv(var_rho) * sz * hh    * hhz_p
                r0_pp    = r0_pp    + vv(var_rho) * sz * hh    * hhz_pp
                
                ! --- Ion temperature
                Ti0       = Ti0       + ( vv(var_Ti) + vv(var_T)/2.d0 ) * sz * hh    * hhz
                Ti0_s     = Ti0_s     + ( vv(var_Ti) + vv(var_T)/2.d0 ) * sz * hh_s  * hhz
                Ti0_t     = Ti0_t     + ( vv(var_Ti) + vv(var_T)/2.d0 ) * sz * hh_t  * hhz
                Ti0_ss    = Ti0_ss    + ( vv(var_Ti) + vv(var_T)/2.d0 ) * sz * hh_ss * hhz
                Ti0_tt    = Ti0_tt    + ( vv(var_Ti) + vv(var_T)/2.d0 ) * sz * hh_tt * hhz
                Ti0_st    = Ti0_st    + ( vv(var_Ti) + vv(var_T)/2.d0 ) * sz * hh_st * hhz
                Ti0_p     = Ti0_p     + ( vv(var_Ti) + vv(var_T)/2.d0 ) * sz * hh    * hhz_p
                Ti0_pp    = Ti0_pp    + ( vv(var_Ti) + vv(var_T)/2.d0 ) * sz * hh    * hhz_pp
                
                ! --- Electron temperature
                Te0       = Te0       + ( vv(var_Te) + vv(var_T)/2.d0 ) * sz * hh    * hhz
                Te0_s     = Te0_s     + ( vv(var_Te) + vv(var_T)/2.d0 ) * sz * hh_s  * hhz
                Te0_t     = Te0_t     + ( vv(var_Te) + vv(var_T)/2.d0 ) * sz * hh_t  * hhz
                Te0_ss    = Te0_ss    + ( vv(var_Te) + vv(var_T)/2.d0 ) * sz * hh_ss * hhz
                Te0_tt    = Te0_tt    + ( vv(var_Te) + vv(var_T)/2.d0 ) * sz * hh_tt * hhz
                Te0_st    = Te0_st    + ( vv(var_Te) + vv(var_T)/2.d0 ) * sz * hh_st * hhz
                Te0_p     = Te0_p     + ( vv(var_Te) + vv(var_T)/2.d0 ) * sz * hh    * hhz_p
                Te0_pp    = Te0_pp    + ( vv(var_Te) + vv(var_T)/2.d0 ) * sz * hh    * hhz_pp

                ! --- Temperature (ion + electron) in models .ne. 400 & 502
                T0       = T0       + ( vv(var_T) + vv(var_Ti) + vv(var_Te) ) * sz * hh    * hhz
                T0_s     = T0_s     + ( vv(var_T) + vv(var_Ti) + vv(var_Te) ) * sz * hh_s  * hhz
                T0_t     = T0_t     + ( vv(var_T) + vv(var_Ti) + vv(var_Te) ) * sz * hh_t  * hhz
                T0_ss    = T0_ss    + ( vv(var_T) + vv(var_Ti) + vv(var_Te) ) * sz * hh_ss * hhz
                T0_tt    = T0_tt    + ( vv(var_T) + vv(var_Ti) + vv(var_Te) ) * sz * hh_tt * hhz
                T0_st    = T0_st    + ( vv(var_T) + vv(var_Ti) + vv(var_Te) ) * sz * hh_st * hhz
                T0_p     = T0_p     + ( vv(var_T) + vv(var_Ti) + vv(var_Te) ) * sz * hh    * hhz_p
                T0_pp    = T0_pp    + ( vv(var_T) + vv(var_Ti) + vv(var_Te) ) * sz * hh    * hhz_pp

                ! --- Parallel Velocity

                Vpar0    = Vpar0    + vv(var_Vpar) * sz * hh    * hhz
                Vpar0_s  = Vpar0_s  + vv(var_Vpar) * sz * hh_s  * hhz
                Vpar0_t  = Vpar0_t  + vv(var_Vpar) * sz * hh_t  * hhz
                Vpar0_ss = Vpar0_ss + vv(var_Vpar) * sz * hh_ss * hhz
                Vpar0_tt = Vpar0_tt + vv(var_Vpar) * sz * hh_tt * hhz
                Vpar0_st = Vpar0_st + vv(var_Vpar) * sz * hh_st * hhz
                Vpar0_p  = Vpar0_p  + vv(var_Vpar) * sz * hh    * hhz_p
                Vpar0_pp = Vpar0_pp + vv(var_Vpar) * sz * hh    * hhz_pp

                rn0       = rn0       + vv(var_rhon) * sz * hh    * hhz
                rn0_s     = rn0_s     + vv(var_rhon) * sz * hh_s  * hhz
                rn0_t     = rn0_t     + vv(var_rhon) * sz * hh_t  * hhz
                rn0_ss    = rn0_ss    + vv(var_rhon) * sz * hh_ss * hhz
                rn0_tt    = rn0_tt    + vv(var_rhon) * sz * hh_tt * hhz
                rn0_st    = rn0_st    + vv(var_rhon) * sz * hh_st * hhz
                rn0_p     = rn0_p     + vv(var_rhon) * sz * hh    * hhz_p
                rn0_pp    = rn0_pp    + vv(var_rhon) * sz * hh    * hhz_pp

                ! --- AR
                AR0      = AR0      + vv(var_AR) * sz * hh    * hhz
                AR0_p    = AR0_p    + vv(var_AR) * sz * hh    * hhz_p
                AR0_s    = AR0_s    + vv(var_AR) * sz * hh_s  * hhz
                AR0_t    = AR0_t    + vv(var_AR) * sz * hh_t  * hhz
                AR0_ss   = AR0_ss   + vv(var_AR) * sz * hh_ss * hhz
                AR0_tt   = AR0_tt   + vv(var_AR) * sz * hh_tt * hhz
                AR0_st   = AR0_st   + vv(var_AR) * sz * hh_st * hhz
                AR0_sp   = AR0_sp   + vv(var_AR) * sz * hh_s  * hhz_p
                AR0_tp   = AR0_tp   + vv(var_AR) * sz * hh_t  * hhz_p
                AR0_pp   = AR0_pp   + vv(var_AR) * sz * hh    * hhz_pp

                ! --- AZ
                AZ0      = AZ0      + vv(var_AZ) * sz * hh    * hhz
                AZ0_p    = AZ0_p    + vv(var_AZ) * sz * hh    * hhz_p
                AZ0_s    = AZ0_s    + vv(var_AZ) * sz * hh_s  * hhz
                AZ0_t    = AZ0_t    + vv(var_AZ) * sz * hh_t  * hhz
                AZ0_ss   = AZ0_ss   + vv(var_AZ) * sz * hh_ss * hhz
                AZ0_tt   = AZ0_tt   + vv(var_AZ) * sz * hh_tt * hhz
                AZ0_st   = AZ0_st   + vv(var_AZ) * sz * hh_st * hhz
                AZ0_sp   = AZ0_sp   + vv(var_AZ) * sz * hh_s  * hhz_p
                AZ0_tp   = AZ0_tp   + vv(var_AZ) * sz * hh_t  * hhz_p
                AZ0_pp   = AZ0_pp   + vv(var_AZ) * sz * hh    * hhz_pp

                ! --- A3
                A30      = A30      + vv(var_A3) * sz * hh    * hhz
                A30_p    = A30_p    + vv(var_A3) * sz * hh    * hhz_p
                A30_s    = A30_s    + vv(var_A3) * sz * hh_s  * hhz
                A30_t    = A30_t    + vv(var_A3) * sz * hh_t  * hhz
                A30_ss   = A30_ss   + vv(var_A3) * sz * hh_ss * hhz
                A30_tt   = A30_tt   + vv(var_A3) * sz * hh_tt * hhz
                A30_st   = A30_st   + vv(var_A3) * sz * hh_st * hhz
                A30_pp   = A30_pp   + vv(var_A3) * sz * hh    * hhz_pp
                A30_sp   = A30_sp   + vv(var_A3) * sz * hh_s  * hhz_p
                A30_tp   = A30_tp   + vv(var_A3) * sz * hh_t  * hhz_p


                ! --- Deltas
                do k = 1, n_var
                  delta_g(k) = delta_g(k) + nodes(i)%deltas(i_tor,j,k) * sz * hh    * hhz
                  delta_s(k) = delta_s(k) + nodes(i)%deltas(i_tor,j,k) * sz * hh_s  * hhz
                  delta_t(k) = delta_t(k) + nodes(i)%deltas(i_tor,j,k) * sz * hh_t  * hhz
                end do
                
              end do
            end do
          end do
          
          ! --- Construct Cartesian Derivatives of Variables.
          ps0_R    = (   Z_t * ps0_s - Z_s * ps0_t ) / xjac
          ps0_Z    = ( - R_t * ps0_s + R_s * ps0_t ) / xjac
          ps0_RR   = (ps0_ss * Z_t**2 - 2.d0*ps0_st * Z_s*Z_t + ps0_tt * Z_s**2  &
                      + ps0_s * (Z_st*Z_t - Z_tt*Z_s )                           &
                      + ps0_t * (Z_st*Z_s - Z_ss*Z_t ) )       / xjac**2         &
                    - xjac_R * (ps0_s * Z_t - ps0_t * Z_s)     / xjac**2
          ps0_ZZ   = (ps0_ss * R_t**2 - 2.d0*ps0_st * R_s*R_t + ps0_tt * R_s**2  &
                      + ps0_s * (R_st*R_t - R_tt*R_s )                           &
                      + ps0_t * (R_st*R_s - R_ss*R_t ) )       / xjac**2         &
                    - xjac_Z * (- ps0_s * R_t + ps0_t * R_s )  / xjac**2
          ps0_RZ   = (- ps0_ss * Z_t*R_t - ps0_tt * R_s*Z_s                      &
                      + ps0_st * (Z_s*R_t  + Z_t*R_s  )                          &
                      - ps0_s  * (R_st*Z_t - R_tt*Z_s )                          &
                      - ps0_t  * (R_st*Z_s - R_ss*Z_t )  )     / xjac**2         &
                    - xjac_R * (- ps0_s * R_t + ps0_t * R_s )  / xjac**2
          
          u0_R     = (   Z_t * u0_s - Z_s * u0_t ) / xjac
          u0_Z     = ( - R_t * u0_s + R_s * u0_t ) / xjac
          u0_RR    = (u0_ss * Z_t**2 - 2.d0*u0_st * Z_s*Z_t + u0_tt * Z_s**2  &
                      + u0_s * (Z_st*Z_t - Z_tt*Z_s )                         &
                      + u0_t * (Z_st*Z_s - Z_ss*Z_t ) )      / xjac**2        &
                    - xjac_R * (u0_s * Z_t - u0_t * Z_s)     / xjac**2
          u0_ZZ    = (u0_ss * R_t**2 - 2.d0*u0_st * R_s*R_t + u0_tt * R_s**2  &
                      + u0_s * (R_st*R_t - R_tt*R_s )                         &
                      + u0_t * (R_st*R_s - R_ss*R_t ) )      / xjac**2        &
                    - xjac_Z * (- u0_s * R_t + u0_t * R_s )  / xjac**2
          u0_RZ    = (- u0_ss * Z_t*R_t - u0_tt * R_s*Z_s                     &
                      + u0_st * (Z_s*R_t  + Z_t*R_s  )                        &
                      - u0_s  * (R_st*Z_t - R_tt*Z_s )                        &
                      - u0_t  * (R_st*Z_s - R_ss*Z_t )  )    / xjac**2        &
                    - xjac_R * (- u0_s * R_t + u0_t * R_s )  / xjac**2
          vv2      = BigR**2 *  ( u0_R * u0_R + u0_Z *u0_Z  )
          !----------------- simplified version of 2nd derivatives (for some unknown reason this is more stable!)
          !u0_RR = (  u0_ss * Z_t**2  + u0_tt * Z_s**2  - 2.d0*u0_st * Z_s*Z_t                  ) / xjac**2
          !u0_ZZ = (  u0_ss * R_t**2  + u0_tt * R_s**2  - 2.d0*u0_st * R_s*R_t                  ) / xjac**2
          !u0_RZ = (- u0_ss * Z_t*R_t - u0_tt * R_s*Z_s +      u0_st * (Z_s*R_t + Z_t*R_s) ) / xjac**2
          
          zj0_R    = (   Z_t * zj0_s - Z_s * zj0_t ) / xjac
          zj0_Z    = ( - R_t * zj0_s + R_s * zj0_t ) / xjac
          zj0_RR   = (zj0_ss * Z_t**2 - 2.d0*zj0_st * Z_s*Z_t + zj0_tt * Z_s**2  &
                      + zj0_s * (Z_st*Z_t - Z_tt*Z_s )                           &
                      + zj0_t * (Z_st*Z_s - Z_ss*Z_t ) )       / xjac**2         &
                    - xjac_R * (zj0_s * Z_t - zj0_t * Z_s)     / xjac**2
          zj0_ZZ   = (zj0_ss * R_t**2 - 2.d0*zj0_st * R_s*R_t + zj0_tt * R_s**2  &
                      + zj0_s * (R_st*R_t - R_tt*R_s )                           &
                      + zj0_t * (R_st*R_s - R_ss*R_t ) )       / xjac**2         &
                    - xjac_Z * (- zj0_s * R_t + zj0_t * R_s )  / xjac**2
          zj0_RZ   = (- zj0_ss * Z_t*R_t - zj0_tt * R_s*Z_s                      &
                      + zj0_st * (Z_s*R_t  + Z_t*R_s  )                          &
                      - zj0_s  * (R_st*Z_t - R_tt*Z_s )                          &
                      - zj0_t  * (R_st*Z_s - R_ss*Z_t )  )     / xjac**2         &
                    - xjac_R * (- zj0_s * R_t + zj0_t * R_s )  / xjac**2
          
          w0_R     = (   Z_t * w0_s - Z_s * w0_t ) / xjac
          w0_Z     = ( - R_t * w0_s + R_s * w0_t ) / xjac
          w0_RR    = (w0_ss * Z_t**2 - 2.d0*w0_st * Z_s*Z_t + w0_tt * Z_s**2  &
                      + w0_s * (Z_st*Z_t - Z_tt*Z_s )                         &
                      + w0_t * (Z_st*Z_s - Z_ss*Z_t ) )      / xjac**2        &
                    - xjac_R * (w0_s * Z_t - w0_t * Z_s)     / xjac**2
          w0_ZZ    = (w0_ss * R_t**2 - 2.d0*w0_st * R_s*R_t + w0_tt * R_s**2  &
                      + w0_s * (R_st*R_t - R_tt*R_s )                         &
                      + w0_t * (R_st*R_s - R_ss*R_t ) )      / xjac**2        &
                    - xjac_Z * (- w0_s * R_t + w0_t * R_s )  / xjac**2
          w0_RZ    = (- w0_ss * Z_t*R_t - w0_tt * R_s*Z_s                     &
                      + w0_st * (Z_s*R_t  + Z_t*R_s  )                        &
                      - w0_s  * (R_st*Z_t - R_tt*Z_s )                        &
                      - w0_t  * (R_st*Z_s - R_ss*Z_t )  )    / xjac**2        &
                    - xjac_R * (- w0_s * R_t + w0_t * R_s )  / xjac**2
          !----------------- simplified version of 2nd derivatives (for some unknown reason this is more stable!)
          !w0_RR = (  w0_ss * Z_t**2  + w0_tt * Z_s**2  - 2.d0*w0_st * Z_s*Z_t                  ) / xjac**2
          !w0_ZZ = (  w0_ss * R_t**2  + w0_tt * R_s**2  - 2.d0*w0_st * R_s*R_t                  ) / xjac**2
          !w0_RZ = (- w0_ss * Z_t*R_t - w0_tt * R_s*Z_s +      w0_st * (Z_s*R_t + Z_t*R_s) ) / xjac**2
          
          r0_R     = (   Z_t * r0_s - Z_s * r0_t ) / xjac
          r0_Z     = ( - R_t * r0_s + R_s * r0_t ) / xjac
          r0_RR    = (r0_ss * Z_t**2 - 2.d0*r0_st * Z_s*Z_t + r0_tt * Z_s**2  &
                      + r0_s * (Z_st*Z_t - Z_tt*Z_s )                         &
                      + r0_t * (Z_st*Z_s - Z_ss*Z_t ) )      / xjac**2        &
                    - xjac_R * (r0_s * Z_t - r0_t * Z_s)     / xjac**2
          r0_ZZ    = (r0_ss * R_t**2 - 2.d0*r0_st * R_s*R_t + r0_tt * R_s**2  &
                      + r0_s * (R_st*R_t - R_tt*R_s )                         &
                      + r0_t * (R_st*R_s - R_ss*R_t ) )      / xjac**2        &
                    - xjac_Z * (- r0_s * R_t + r0_t * R_s )  / xjac**2
          r0_RZ    = (- r0_ss * Z_t*R_t - r0_tt * R_s*Z_s                     &
                      + r0_st * (Z_s*R_t  + Z_t*R_s  )                        &
                      - r0_s  * (R_st*Z_t - R_tt*Z_s )                        &
                      - r0_t  * (R_st*Z_s - R_ss*Z_t )  )    / xjac**2        &
                    - xjac_R * (- r0_s * R_t + r0_t * R_s )  / xjac**2
          r0_hat   = BigR**2 * r0
          r0_R_hat = 2.d0 * BigR * BigR_R  * r0 + BigR**2 * r0_R
          r0_Z_hat = BigR**2 * r0_Z
          
          T0_R     = (   Z_t * T0_s  - Z_s * T0_t ) / xjac
          T0_Z     = ( - R_t * T0_s  + R_s * T0_t ) / xjac

          Ti0_R    = (   Z_t * Ti0_s  - Z_s * Ti0_t ) / xjac
          Ti0_Z    = ( - R_t * Ti0_s  + R_s * Ti0_t ) / xjac
          Te0_R    = (   Z_t * Te0_s  - Z_s * Te0_t ) / xjac
          Te0_Z    = ( - R_t * Te0_s  + R_s * Te0_t ) / xjac

          T0_RR    = (T0_ss * Z_t**2 - 2.d0*T0_st * Z_s*Z_t + T0_tt * Z_s**2 &
                      + T0_s * (Z_st*Z_t - Z_tt*Z_s )                                 &
                      + T0_t * (Z_st*Z_s - Z_ss*Z_t ) )        / xjac**2        &
                     - xjac_R * (T0_s * Z_t - T0_t * Z_s)     / xjac**2
          T0_ZZ    = (T0_ss * R_t**2 - 2.d0*T0_st * R_s*R_t + T0_tt * R_s**2 &
                      + T0_s * (R_st*R_t - R_tt*R_s )                                 &
                      + T0_t * (R_st*R_s - R_ss*R_t ) )        / xjac**2        &
                     - xjac_Z * (- T0_s * R_t + T0_t * R_s )  / xjac**2
          T0_RZ    = (- T0_ss * Z_t*R_t - T0_tt * R_s*Z_s                         &
                       + T0_st * (Z_s*R_t  + Z_t*R_s  )                           &
                       - T0_s  * (R_st*Z_t - R_tt*Z_s )                           &
                       - T0_t  * (R_st*Z_s - R_ss*Z_t )  )       / xjac**2      &
                       - xjac_R * (- T0_s * R_t + T0_t * R_s )  / xjac**2
          T0_ps0_R = T0_RR * ps0_Z - T0_RZ * ps0_R + T0_R * ps0_RZ - T0_Z * ps0_RR
          T0_ps0_Z = T0_RZ * ps0_Z - T0_ZZ * ps0_R + T0_R * ps0_ZZ - T0_Z * ps0_RZ
          !----------------- simplified version of 2nd derivatives (for some unknown reason this is more stable!)
          !T0_RR = (  T0_ss * Z_t**2  + T0_tt * Z_s**2  - 2.d0*T0_st * Z_s*Z_t             ) / xjac**2
          !T0_ZZ = (  T0_ss * R_t**2  + T0_tt * R_s**2  - 2.d0*T0_st * R_s*R_t             ) / xjac**2
          !T0_RZ = (- T0_ss * Z_t*R_t - T0_tt * R_s*Z_s +      T0_st * (Z_s*R_t + Z_t*R_s) ) / xjac**2
          
          Vpar0_R  = (   Z_t * Vpar0_s - Z_s * Vpar0_t ) / xjac
          Vpar0_Z  = ( - R_t * Vpar0_s + R_s * Vpar0_t ) / xjac
          Vpar0_RR = (Vpar0_ss * Z_t**2 - 2.d0*Vpar0_st * Z_s*Z_t + Vpar0_tt * Z_s**2  &
                      + Vpar0_s * (Z_st*Z_t - Z_tt*Z_s )                               &
                      + Vpar0_t * (Z_st*Z_s - Z_ss*Z_t ) )         / xjac**2           &
                    - xjac_R * (Vpar0_s * Z_t - Vpar0_t * Z_s)     / xjac**2
          Vpar0_ZZ = (Vpar0_ss * R_t**2 - 2.d0*Vpar0_st * R_s*R_t + Vpar0_tt * R_s**2  &
                      + Vpar0_s * (R_st*R_t - R_tt*R_s )                               &
                      + Vpar0_t * (R_st*R_s - R_ss*R_t ) )         / xjac**2           &
                    - xjac_Z * (- Vpar0_s * R_t + Vpar0_t * R_s )  / xjac**2
          Vpar0_RZ = (- Vpar0_ss * Z_t*R_t - Vpar0_tt * R_s*Z_s                        &
                      + Vpar0_st * (Z_s*R_t  + Z_t*R_s  )                                 &
                      - Vpar0_s  * (R_st*Z_t - R_tt*Z_s )                                 &
                      - Vpar0_t  * (R_st*Z_s - R_ss*Z_t )  )       / xjac**2           &
                    - xjac_R * (- Vpar0_s * R_t + Vpar0_t * R_s )  / xjac**2
          
          !delta_u_R  = (   Z_t * delta_s(2) - Z_s * delta_t(2) ) / xjac
          !delta_u_Z  = ( - R_t * delta_s(2) + R_s * delta_t(2) ) / xjac
          !delta_ps_R = (   Z_t * delta_s(1) - Z_s * delta_t(1) ) / xjac
          !delta_ps_Z = ( - R_t * delta_s(1) + R_s * delta_t(1) ) / xjac
          
          AR0_R = (   Z_t * AR0_s  - Z_s * AR0_t ) / xjac
          AR0_Z = ( - R_t * AR0_s  + R_s * AR0_t ) / xjac
          AZ0_R = (   Z_t * AZ0_s  - Z_s * AZ0_t ) / xjac
          AZ0_Z = ( - R_t * AZ0_s  + R_s * AZ0_t ) / xjac
          A30_R = (   Z_t * A30_s  - Z_s * A30_t ) / xjac
          A30_Z = ( - R_t * A30_s  + R_s * A30_t ) / xjac

          AR0_Rp = (   Z_t * AR0_sp  - Z_s * AR0_tp ) / xjac
          AR0_Zp = ( - R_t * AR0_sp  + R_s * AR0_tp ) / xjac
          AZ0_Rp = (   Z_t * AZ0_sp  - Z_s * AZ0_tp ) / xjac
          AZ0_Zp = ( - R_t * AZ0_sp  + R_s * AZ0_tp ) / xjac 
          A30_Rp = (   Z_t * A30_sp  - Z_s * A30_tp ) / xjac
          A30_Zp = ( - R_t * A30_sp  + R_s * A30_tp ) / xjac 

          AR0_RR   = (AR0_ss * Z_t**2 - 2.d0*AR0_st * Z_s*Z_t + AR0_tt * Z_s**2    &
                      + AR0_s * (Z_st*Z_t - Z_tt*Z_s )                             &
                      + AR0_t * (Z_st*Z_s - Z_ss*Z_t ) )       / xjac**2           &
                      - xjac_R * (AR0_s * Z_t - AR0_t * Z_s)     / xjac**2
          AR0_ZZ   = (AR0_ss * R_t**2 - 2.d0*AR0_st * R_s*R_t + AR0_tt * R_s**2    &
                      + AR0_s * (R_st*R_t - R_tt*R_s )                             &
                      + AR0_t * (R_st*R_s - R_ss*R_t ) )       / xjac**2           &
                      - xjac_Z * (- AR0_s * R_t + AR0_t * R_s )  / xjac**2
          AR0_RZ   = (- AR0_ss * Z_t*R_t - AR0_tt * R_s*Z_s                        & 
                      + AR0_st * (Z_s*R_t  + Z_t*R_s  )                            &
                      - AR0_s  * (R_st*Z_t - R_tt*Z_s )                            &
                      - AR0_t * (R_st*Z_s  - R_ss*Z_t ) )  / xjac**2               &
                      - xjac_R * (- AR0_s * R_t + AR0_t * R_s )   / xjac**2        

          AZ0_RR   = (AZ0_ss * Z_t**2 - 2.d0*AZ0_st * Z_s*Z_t + AZ0_tt * Z_s**2    &
                      + AZ0_s * (Z_st*Z_t - Z_tt*Z_s )                             &
                      + AZ0_t * (Z_st*Z_s - Z_ss*Z_t ) )       / xjac**2           &
                      - xjac_R * (AZ0_s * Z_t - AZ0_t * Z_s)     / xjac**2
          AZ0_ZZ   = (AZ0_ss * R_t**2 - 2.d0*AZ0_st * R_s*R_t + AZ0_tt * R_s**2    &
                      + AZ0_s * (R_st*R_t - R_tt*R_s )                             &
                      + AZ0_t * (R_st*R_s - R_ss*R_t ) )       / xjac**2           &
                      - xjac_Z * (- AZ0_s * R_t + AZ0_t * R_s )  / xjac**2
          AZ0_RZ   = (- AZ0_ss * Z_t*R_t - AZ0_tt * R_s*Z_s                        &
                      + AZ0_st * (Z_s*R_t  + Z_t*R_s  )                            &
                      - AZ0_s  * (R_st*Z_t - R_tt*Z_s )                            &
                      - AZ0_t * (R_st*Z_s  - R_ss*Z_t ) )  / xjac**2               &
                      - xjac_R * (- AZ0_s * R_t + AZ0_t * R_s )   / xjac**2        

          A30_RR   = (A30_ss * Z_t**2 - 2.d0*A30_st * Z_s*Z_t + A30_tt * Z_s**2    &
                      + A30_s * (Z_st*Z_t - Z_tt*Z_s )                             &
                      + A30_t * (Z_st*Z_s - Z_ss*Z_t ) )       / xjac**2           &
                      - xjac_R * (A30_s * Z_t - A30_t * Z_s)     / xjac**2
          A30_ZZ   = (A30_ss * R_t**2 - 2.d0*A30_st * R_s*R_t + A30_tt * R_s**2    &
                      + A30_s * (R_st*R_t - R_tt*R_s )                             &
                      + A30_t * (R_st*R_s - R_ss*R_t ) )       / xjac**2           &
                      - xjac_Z * (- A30_s * R_t + A30_t * R_s )  / xjac**2
          A30_RZ   = (- A30_ss * Z_t*R_t - A30_tt * R_s*Z_s                        &
                      + A30_st * (Z_s*R_t  + Z_t*R_s  )                            &
                      - A30_s  * (R_st*Z_t - R_tt*Z_s )                            &
                      - A30_t * (R_st*Z_s  - R_ss*Z_t ) )  / xjac**2               &
                      - xjac_R * (- A30_s * R_t + A30_t * R_s )   / xjac**2

          Fprofile_R = (   Z_t * Fprofile_s  - Z_s * Fprofile_t ) / xjac
          Fprofile_Z = ( - R_t * Fprofile_s  + R_s * Fprofile_t ) / xjac
 
          ! --- Pressure
          P0       = r0    * T0
          P0_R     = r0_R  * T0 + r0 * T0_R
          P0_Z     = r0_Z  * T0 + r0 * T0_Z
          P0_s     = r0_s  * T0 + r0 * T0_s
          P0_t     = r0_t  * T0 + r0 * T0_t
          P0_p     = r0_p  * T0 + r0 * T0_p
          P0_pp    = r0_pp * T0 + r0 * T0_pp + 2.d0 * r0_p * T0_p
          P0_RR    = r0_RR * T0 + r0 * T0_RR + 2.d0 * r0_R * T0_R
          P0_ZZ    = r0_ZZ * T0 + r0 * T0_ZZ + 2.d0 * r0_Z * T0_Z
          P0_RZ    = r0_RZ * T0 + r0 * T0_RZ + r0_R * T0_Z + r0_Z * T0_R
 
          rn0_R    = (   Z_t * rn0_s - Z_s * rn0_t ) / xjac
          rn0_Z    = ( - R_t * rn0_s + R_s * rn0_t ) / xjac

          ! --- Some things related to the magnetic field
#ifdef fullmhd
          BR   = ( A30_Z - AZ0_p )/ BigR
          BZ   = ( AR0_p - A30_R )/ BigR
          Btor = ( AZ0_R - AR0_Z )    + Fprofile / BigR
          BR_R = -1/BigR**2 * ( A30_Z - AZ0_p ) + ( A30_RZ - AZ0_Rp )/ BigR
          BR_Z = ( A30_ZZ - AZ0_Zp )/ BigR
          BR_p = ( A30_Zp - AZ0_pp )/ BigR
          BZ_R = -1/BigR**2 * ( AR0_p - A30_R ) + ( AR0_Rp - A30_RR )/ BigR
          BZ_Z = ( AR0_Zp - A30_RZ )/ BigR
          BZ_p = ( AR0_pp - A30_Rp )/BigR
          BP_R = ( AZ0_RR - AR0_RZ ) + Fprofile_R/BigR - Fprofile/BigR**2
          BP_Z = ( AZ0_RZ - AR0_ZZ ) + Fprofile_Z/BigR
          B_R  = ( BR_R + BZ_R + Bp_R ) / Btot
          B_Z  = ( BR_Z + BZ_Z + Bp_Z ) / Btot

          Btheta   = sqrt( BR*BR + BZ*BZ )
          BB2      = Btor**2 + BR**2 + BZ**2
          Btot     = sqrt(BB2)
          zj0      = - (BZ_R - BR_Z) * BigR
          JpolR    = (R*BP_Z - BZ_p) / R
          JpolZ    = (BR_p - R*BP_R - Btor) / R

          psi_norm = get_psi_n(A30, Z)
          psi_abs  = sqrt(A30_R*A30_R + A30_Z * A30_Z)
          if (psi_abs > 1.d-6) then
            FFprime_loc = zj0 + (R**2.d0) * (A30_R*P0_R + A30_Z*P0_Z)/(psi_abs**2.d0)
          else
            FFprime_loc = zj0 !--- not fully correct, but better than to put 0...
          endif

#else
          BB2      = (F0*F0 + ps0_R * ps0_R + ps0_Z * ps0_Z ) / BigR**2
          Btot     = sqrt(BB2)
          BR       = + ps0_Z / BigR
          BZ       = - ps0_R / BigR
          Btor     = + F0    / BigR
          Bnorm    = BR*nmlR + BZ*nmlZ
          Btan     = BR*nmlZ - BZ*nmlR
          Btheta   = sqrt(ps0_R*ps0_R + ps0_Z * ps0_Z) / BigR
          JpolR    = ( -zj0 * BR - R * P0_Z ) / F0
          JpolZ    = ( -zj0 * BZ + R * P0_R ) / F0
          psi_norm = get_psi_n(ps0, Z)
          psi_abs  = sqrt(ps0_R*ps0_R + ps0_Z * ps0_Z)

          BR_R     = + ps0_RZ / BigR - ps0_Z / BigR**2
          BR_Z     = + ps0_ZZ / BigR
          BZ_R     = - ps0_RR / BigR + ps0_R / BigR**2
          BZ_Z     = - ps0_RZ / BigR
          Bp_R     = - F0     / BigR**2
          Bp_Z     =   0.
          B_R      = ( BR_R + BZ_R + Bp_R ) / Btot
          B_Z      = ( BR_Z + BZ_Z + Bp_Z ) / Btot

          if (psi_abs > 1.d-6) then
            FFprime_loc = zj0 + (R**2.d0) * (ps0_R*P0_R + ps0_Z*P0_Z)/(psi_abs**2.d0)
          else
            FFprime_loc = zj0 !--- not fully correct, but better than to put 0...
          endif
#endif

          Kappa_R    = ( Btot*BR*BR_R - BR*BR*B_R   - BZ*BR*B_Z   + Btot*BZ*BR_Z - Btot*Btor**2/BigR ) / Btot**3.
          Kappa_Z    = ( Btot*BR*BZ_R - BR*BZ*B_R   - BZ*BZ*B_Z   + Btot*BZ*BZ_Z                     ) / Btot**3.
          Kappa_phi  = ( Btot*BR*Bp_R - BR*Btor*B_R - BZ*Btor*B_Z + Btot*BZ*Bp_Z + Btot*Btor*BR/BigR ) / Btot**3.

          Jtor        = -zj0/BigR
          Jpol        = FFprime_loc * Btheta     / F0     !Jpol = F' Bpol
          Jpar        = (JpolR*BR + JpolZ*BZ + Jtor*Btor) / Btot
          Jpar_ionsat = r0 * vpar0 * Btot 

          ! --- Velocity
          VR       = -R*u0_Z + vpar0*ps0_Z/R 
          VZ       =  R*u0_R - vpar0*ps0_R/R
          V_phi    =  F0 * vpar0/R
          Vpar_tot =  (VR*BR + VZ*BZ + V_phi*Btor) / Btot
          VperpR   =  VR - Vpar_tot * BR / Btot
          VperpZ   =  VZ - Vpar_tot * BZ / Btot

          ! --- Some input profiles
          if ( eta_T_dependent ) then
            eta_T     = eta   * (corr_neg_temp(T0)/T_0)**(-1.5d0)
            deta_dT   = - eta   * (1.5d0)  * corr_neg_temp(T0)**(-2.5d0) * T_0**(1.5d0)
            d2eta_d2T =   eta   * (3.75d0) * corr_neg_temp(T0)**(-3.5d0) * T_0**(1.5d0)
          else
            eta_T     = eta
            deta_dT   = 0.d0
            d2eta_d2T = 0.d0
          end if
          
          if ( visco_T_dependent ) then
            visco_T   = visco * (corr_neg_temp(T0)/T_0)**(-1.5d0)
            dvisco_dT = - visco * (1.5d0)  * corr_neg_temp(T0)**(-2.5d0) * T_0**(1.5d0)
          else
            visco_T   = visco
            dvisco_dT = 0.d0
          end if
          
          if ( ZKpar_T_dependent ) then
            ZKpar_T   = ZK_par * (corr_neg_temp(T0)/T_0)**(+2.5d0)
            dZKpar_dT = ZK_par * (2.5d0)  * corr_neg_temp(T0)**(+1.5d0) * T_0**(-2.5d0)
            if (ZKpar_T .gt. ZK_par_max) then
              ZKpar_T   = Zk_par_max
              dZKpar_dT = 0.d0
            end if
          else
            ZKpar_T   = ZK_par
            dZKpar_dT = 0.d0
          end if
          
          D_prof  = get_dperp (psi_norm)
          ZK_prof = get_zkperp(psi_norm)

          ! --- Fluxes 
          pres_flux_par =  gamma/(gamma-1.d0) * r0 * T0 * Vpar_tot * Bnorm / Btot          !  p v_par·n
          pres_flux_tot =  gamma/(gamma-1.d0) * r0 * T0 * (VR*nmlR + VZ*nmlZ)              !  p v·n

          kin_flux_par  = 0.5d0*r0* (VR*VR + VZ*VZ + V_phi*V_phi)* Vpar_tot * Bnorm / Btot ! 0.5 nv^2 v_par·n
          kin_flux_tot  = 0.5d0*r0* (VR*VR + VZ*VZ + V_phi*V_phi)* (VR*nmlR + VZ*nmlZ)     ! 0.5 nv^2 v·n 

          ZKpar_flux    = - ZKpar_T *(BR*T0_R + BZ*T0_Z + Btor*T0_p/R) * Bnorm / BB2 / (gamma-1.d0) ! q_par·n 
          ZKperp_flux   = - ZK_prof *( T0_R*nmlR + T0_Z*nmlZ)        / (gamma-1.d0) &                ! q_perp·n
                          + ZK_prof *(BR*T0_R + BZ*T0_Z + Btor*T0_p/R) * Bnorm / BB2 / (gamma-1.d0) 
    
          Dpar_flux     = - D_par  * (BR*r0_R + BZ*T0_Z + Btor*T0_p/R) * Bnorm / BB2
          Dperp_flux    = - D_prof * ( r0_R*nmlR + T0_Z*nmlZ)                       &                              
                          + D_prof * (BR*r0_R + BZ*T0_Z + Btor*T0_p/R) * Bnorm / BB2 
    
          partF_cnv_par =   r0 * Vpar_tot * Bnorm / Btot                           !  p v_par·n
          partF_cnv_tot =   r0 * ( VR * nmlR + VZ * nmlZ )                         !  n v·n
    
#if (defined WITH_Neutrals) && (!defined WITH_Impurities)
          neut_part_flux= -D_neutral_x*rn0_R * nmlR - D_neutral_y * rn0_Z * nmlZ
#else
          neut_part_flux= 0.d0
#endif    
          dpsi_dt   = BigR*(ps0_s*u0_t - ps0_t*u0_s)/xjac + eta_T*zj0 - F0*u0_p 
          ExB_norm  = -dpsi_dt * (ps0_R*nmlR + ps0_Z*nmlZ) / (BigR**2.d0) 
         
          ! --- Other parameters (combination of the main variables)
          Er       = 0.d0
          Vtheta   = 0.d0
	  mach_par = 0.d0
          mach_pol = 0.d0
          vsound   = 0.d0
          Vneo     = 0.d0
          ki_neo   = 0.d0
          mu_neo   = 0.d0
          Vperp_e  = 0.d0
          Vperp_i  = 0.d0
          V_ExB    = 0.d0 
          Vstar_e  = 0.d0
          Vstar_i  = 0.d0
          
          if ( (psi_abs > 1.d-6) .and. (r0 > 1.d-6) .and. (abs(Btheta) > 1.d-6) ) then
            
            Er       = -(u0_R * ps0_R + u0_Z * ps0_Z) / psi_abs   ! radial electric field
            
            Vsound   = sqrt(GAMMA*T0) / sqrt(BB2)                 ! sound speed
            Mach_par = Vpar0 / Vsound                             ! parallel Mach number
            Mach_pol = Vtheta / Vsound                            ! poloidal Mach number
            
            Vtheta   = -1./Btheta * (  ( u0_R + tauIC/r0 * (T0_R*r0 + r0_R*T0) ) * ps0_R  +        &
              ( u0_Z + tauIC/r0 * (T0_Z*r0 + r0_Z*T0) ) * ps0_Z) + Vpar0 * Btheta
            
            Vperp_i  = -1./Btheta * (  ( u0_R + tauIC/r0 * (T0_R*r0 + r0_R*T0) ) * ps0_R  +        &
              ( u0_Z + tauIC/r0 * (T0_Z*r0 + r0_Z*T0) ) * ps0_Z )
            
            Vperp_e  = -1./Btheta * (  ( u0_R - tauIC/r0 * (T0_R*r0 + r0_R*T0) ) * ps0_R  +        &
              ( u0_Z - tauIC/r0 * (T0_Z*r0 + r0_Z*T0) ) * ps0_Z )
            
            V_ExB    = -1./Btheta* ( u0_R*ps0_R + u0_Z*ps0_Z )
            
            Vstar_i  = -1./Btheta * (  tauIC/r0 * (T0_R*r0 + r0_R*T0) * ps0_R  +                   &
              tauIC/r0 * (T0_Z*r0 + r0_Z*T0) * ps0_Z )
            
            Vstar_e  = +1./Btheta * (  tauIC/r0 * (T0_R*r0 + r0_R*T0) * ps0_R  +                   &
              tauIC/r0 * (T0_Z*r0 + r0_Z*T0) * ps0_Z )
            ! ### Warning : in jorek_model=400, Vstar_i .ne. -Vstar_e since T_i .ne. T_e
          end if
          
          if (NEO) then
            if (num_neo_file) then ! (read neoclassical profiles from ascii file)
              call neo_coef( eq%xpoint, eq%xcase, Z, eq%Z_xpoint, Ps0 ,eq%psi_axis, eq%psi_bnd,    &
                mu_neo, ki_neo)
              Vneo = ki_neo / Btheta * tauIC  * ( ps0_R*T0_R + ps0_Z*T0_Z )
            else ! (use constant neoclassical coefficients from namelist input file)
              mu_neo = amu_neo_const
              ki_neo = aki_neo_const
              Vneo   = aki_neo_const / Btheta * tauIC * ( ps0_R*T0_R + ps0_Z*T0_Z )
            end if
          end if
          
          ! --- Coulomb logarithms calculated according to Ref. [L. Hesselow et al, J Plasma Phys 84,
          !     p. 905840605 (2018); doi:10.1017/S0022377818001113] Eq. (2.7) and (2.9):
          Te0_eV     = Te0 / ( EL_CHG * MU_ZERO * central_density * 1.d20 )
          ne0_20     = max(1.d-8, r0) * central_density
          ln_Lambda0 = 14.9 - 0.5 * log( ne0_20 ) + log( Te0_eV / 1000.d0 ) ! Eq. (2.7) at thermal speeds
          ln_Lambda  = 14.6 + 0.5 * log( Te0_eV / ne0_20 )                  ! Eq. (2.9) at relativistic energies
          
          E_par = - R * ( eta_T * zj0 / R**2                                                       &
                        + tauIC / r0 * ( (P0_R * Ps0_Z - P0_Z * Ps0_R) / R + F0 * P0_p / R**2 ) )
          
          E_crit = C_LIGHT**2 * EL_CHG**3 * ln_Lambda * MU_ZERO**2.5 * (central_density*1.d20*central_mass*MASS_PROTON)**1.5 * r0 / ( 4 * PI * MASS_ELECTRON * MASS_PROTON * central_mass )
          
          E_dreicer = EL_CHG**3 * ln_Lambda0 * MU_ZERO**1.5 * (central_density*1.d20*central_mass*MASS_PROTON)**2.5 * r0 / ( 2.d0 * PI * EPS_ZERO**2 * (MASS_PROTON*central_mass)**2 * T0 )
          
#if JOREK_MODEL >= 303
          if (bootstrap) then
            call bootstrap_current(R, Z, eq%R_axis, eq%Z_axis, eq%psi_axis, eq%R_xpoint, eq%Z_xpoint, eq%psi_bnd, psi_norm, ps0, ps0_R,    &
              ps0_Z, r0,  r0_R, r0_Z, Ti0, Ti0_R, Ti0_Z, Te0, Te0_R, Te0_Z, J_boot)
          endif
#else
          J_boot = 0.d0
#endif

#if (defined WITH_Neutrals) && (!defined WITH_Impurities)

   Te_corr_eV = corr_neg_temp(T0)/(2.d0*EL_CHG*MU_ZERO*central_density * 1.d20)

   if (use_imp_adas) then
     call atomic_coeff_deuterium(Te0, Sion_T, dSion_dT, Srec_T, dSrec_dT,        &
                                LradDcont_T, dLradDcont_dT, LradDrays_T, dLradDrays_dT, r0 ) 
     ! Note the input Te0 for atomic_coeff_deuterium should be in JOREK units!!!

    !--------------------------------------------------------
    ! --- Radiation from background impurity
    !--------------------------------------------------------
      ne_SI = corr_neg_dens(r0) * 1.d20 * central_density !electron density (SI)
      r_imp = nimp_bg / (1.d20 * central_density)  ! Background impurity density in JU     

      if (ne_SI > ne_SI_min .and. Te_corr_eV > Te_eV_min .and. nimp_bg > 0) then
        ! Normalization coefficient for radiation rate from SI units (W.m^3) to JOREK units:
        coef_rad_imp = 2.d0/3.d0*MU_ZERO**1.5d0*(central_mass*MASS_PROTON)**0.5d0&
                     *(central_density*1.d20)**2.5d0
        Lrad_imp = 0.0
        call radiation_function_linear(imp_adas(1),imp_cor(1),log10(ne_SI),log10(Te_corr_eV*EL_CHG/K_BOLTZ),Lrad_imp)         
        if (Lrad_imp < 0.) then
          Lrad_imp = 0.
        end if
        if ( units == SI_UNITS ) then
          frad_bg = nimp_bg * Lrad_imp
        else if ( units == JOREK_UNITS ) then
          frad_bg = r_imp * Lrad_imp * coef_rad_imp 
        endif
      else     
        Lrad_imp = 0.
        frad_bg = 0.
      end if   
    else
      if ( trim(imp_type) == 'Ar') then ! Hard-coded fitting exists for argon
        Arad_bg = 2.4d-31
        Brad_bg = 20.
        Crad_bg = 0.8
        if ( units == SI_UNITS ) then
          frad_bg = nimp_bg*Arad_bg*exp(-((log(Te_corr_eV)-log(Brad_bg))**2.)/Crad_bg**2.)
        else if ( units == JOREK_UNITS ) then
          frad_bg = (2./3.)*(1./(central_mass*MASS_PROTON))*((MU_ZERO*central_mass*MASS_PROTON*central_density*1.d20)**(1.5d0))*nimp_bg*Arad_bg*exp(-((log(Te_corr_eV)-log(Brad_bg))**2.)/Crad_bg**2.)
        end if
      else
        write(*,*) "WARNING: hard-coded fitting doesn't exist for  ", trim(imp_type), ",use open adas instead!"
        stop
      end if      
    end if  

#endif

#ifdef WITH_Impurities

          T0_corr = corr_neg_temp(T0,(/5.d-1,5.d-1/))
          Te_corr_eV   = T0_corr/(2.d0*EL_CHG*MU_ZERO*central_density*1.d20)  ! Te in eV
          Te_eV = T0/(2.d0*EL_CHG*MU_ZERO*central_density * 1.d20)
  
          call imp_cor(1)%interp_linear(density=20.,temperature=log10(Te_corr_eV*EL_CHG/K_BOLTZ),z_avg=Z_imp)
	  
          alpha_imp = 0.5*m_i_over_m_imp*(Z_imp+1.) - 1.
          beta_imp  = m_i_over_m_imp*Z_imp - 1.
                  
          r0_corr  = corr_neg_dens(r0,(/1.d-9,1.d-5/),1.d-3)
          rn0_corr = corr_neg_dens(rn0,(/1.d-9,1.d-5 /),1.d-3)
          ne_SI   = (r0_corr + beta_imp * rn0_corr) * 1.d20 * central_density ! electron density (SI)
  
          ! Normalization coefficient for radiation rate from SI units (W.m^3) to JOREK units:
          coef_rad_1 = 2.d0/3.d0*MU_ZERO**1.5d0*(central_mass*MASS_PROTON)**0.5d0*(central_density*1.d20)**2.5d0*m_i_over_m_imp
  
          if (ne_SI > ne_SI_min .and. Te_eV > Te_eV_min .and. rn0 > rn0_min) then
            Lrad = 0.
            call radiation_function(imp_adas(1),imp_cor(1),log10(ne_SI),log10(Te_corr_eV*EL_CHG/K_BOLTZ),Lrad)
            Lrad = Lrad * coef_rad_1
          else
            Lrad = 0.
          end if
  
          ne_SI = ne_SI / 1.d20 / central_density ! Put ne_SI back to JOREK units to have consistent fact_ne factor with other models (see below)
  
#endif

          ! --- Factors for switching between JOREK normalized and SI units.
          if ( units == SI_UNITS ) then
             rho_norm      = central_density *1.d20 * central_mass * mass_proton   ! rho_0 = central mass density
             fact_time     = sqrt(MU_zero*rho_norm)                                ! time factor
             fact_mu_zero  = MU_zero                                               ! division by mu_zero for P and J
             fact_ne       = central_density * 1.d20                               ! factor for n_e
             fact_rho      = central_density * 1.d20 * central_mass*MASS_PROTON    ! factor for rho
             fact_T        = 1.d0 / ( MU_zero * central_density * 1.d20 * EL_CHG ) ! factor for T
             fact_vpar     = sqrt(BB2) / fact_time                                 ! factor for Vpar
             fact_resistiv = sqrt ( MU_zero / rho_norm )                           ! factor for eta == 1 / (factor for visco)
             fact_Er       = F0 / fact_time
             fact_rad      = 1.d0/(2.d0/3.d0*MU_ZERO**1.5d0*(central_mass*MASS_PROTON*central_density*1.d20)**0.5d0)
             fact_flux     = 1.d0/(mu_zero*fact_time)  
          else if ( units == JOREK_UNITS ) then
             rho_norm      = 1.d0
             fact_time     = 1.d0
             fact_mu_zero  = 1.d0
             fact_ne       = 1.d0
             fact_rho      = 1.d0
             fact_T        = 1.d0
             fact_vpar     = 1.d0
             fact_resistiv = 1.d0
             fact_Er       = 1.d0
             fact_rad      = 1.d0
             fact_flux     = 1.d0 
          end if
          
          ! --- factor to calculate ion saturation current in JOREK units
          fact_jsat = EL_CHG * 1.d20 * central_density * sqrt(MU_ZERO/rho_norm) 

          ! --- Now that everything is prepared, evaluate the requested expressions.
          loop_expr: do iexpr = 1, expr_list%n_expr
            
            select case ( trim(expr_list%expr(iexpr)%name) )

              case ( 'index_now' )
                res = real(index_now)

              case ( 'R' )
                res = R
                
              case ( 'Z', 'z' )
                res = Z
                
              case ( 'phi' )
                res = phi
                
              case ( 'theta' )
                res = theta
                
              case ( 'theta_star' )
                res = pol_pos%theta_star
                
              case ( 'length' )
                res = pol_pos%length
                
              case ( 'r_minor' )
                res = pol_pos%r_minor
                
              case ( 'xjac' )
                res = xjac
                
              case ( 'x' )
                res = x_cart
                
              case ( 'y' )
                res = y_cart
                
              case ( 't' )
                res = t_now * fact_time
                
              case ( 'Psi' )
                res = ps0
 
              case ( 'dPsi_dt' )
                res = dpsi_dt
                
              case ( 'Psi_N' )
                res = psi_norm
                
              case ( 'u' )
                res = u0
                
              case ( 'Phi' )
                res = u0 * F0 !### sign?
                
              case ( 'zj' )
                res = zj0 / fact_mu_zero
                
              case ( 'omega' )
                res = w0
                
              case ( 'rho' )
                res = r0 * fact_rho
                
              case ( 'ne' )
#ifdef WITH_Impurities
                res = ne_SI * fact_ne 
#else
                res = r0 * fact_ne
#endif

#ifdef WITH_Impurities
              case ( 'nimp' )
                res = rn0 * fact_ne * m_i_over_m_imp
#endif

              case ( 'T' )
                res = T0 * fact_T
              
              case ( 'Te' )
                res = T0 * fact_T / 2.d0
              
              case ( 'vpar' )
                res = Vpar0 * fact_vpar
                
              case ( 'eta_T' )
                res = eta_t * fact_resistiv
                
              case ( 'visco_T' )
                res = visco_t / fact_resistiv
                
              case ( 'zkpar_T' )
                res = zkpar_t / fact_time
                
              case ( 'dprof' ) 
                res = d_prof / fact_time
                
              case ( 'zkprof' )
                res = zk_prof / fact_time
                
              case ( 'pres' )
                res = P0 / fact_mu_zero
                
              case ( 'B_abs' )
                res = sqrt(BB2)
                
              case ( 'Btor' )
                res = Btor
                
              case ( 'BR' )
                res = BR
                
              case ( 'BZ' )
                res = BZ
                
              case ( 'Btheta' )
                res = Btheta
	        
              case ( 'A_R' )
                res = AR0
                
              case ( 'A_Z' )
                res = AZ0  
                
              case ( 'A_3' )
                res = A30
                
              case ( 'currdens' )
                res = -zj0 / R / fact_mu_zero

              case ( 'JR' )
                res = JpolR / fact_mu_zero

              case ( 'JZ' )
                res = JpolZ / fact_mu_zero

              case ( 'Jtor' )
                res = Jtor / fact_mu_zero

              case ( 'FFprime_loc' )
                res = FFprime_loc

              case ( 'Jpol' )
                res = Jpol / fact_mu_zero
                
              case ( 'JxB_R' )
                res = (JpolZ*Btor-Jtor*BZ) / fact_mu_zero

              case ( 'JxB_Z' )
                res = (-JpolR*Btor+Jtor*BR) / fact_mu_zero

              case ( 'JxB_phi' )
                res = (JpolR*BZ-JpolZ*BR) / fact_mu_zero

              case ( 'gradP_R' )
                res = P0_R / fact_mu_zero
                
              case ( 'gradP_Z' )
                res = P0_Z / fact_mu_zero
 
              case ( 'gradP_phi' )
                res = P0_p / BigR / fact_mu_zero

              case ( 'gradP_B' )
                res = (BR*P0_R+BZ*P0_Z+Btor*P0_p/BigR) / Btot / fact_mu_zero
 
              case ( 'gradPdotCurv' )
                res = ( P0_R*Kappa_R + P0_Z*Kappa_Z + P0_p / BigR * Kappa_phi ) / fact_mu_zero
      
              case ( 'curvat_R' )
                res = Kappa_R
   
              case ( 'curvat_Z' )
                res = Kappa_Z

              case ( 'curvat_phi' )
                res = Kappa_phi

              case ( 'Er' )
                res = Er * fact_Er
                
              case ( 'Vtheta_i' )
                res = Vtheta / fact_time
                
              case ( 'Mach_par' )
                res = Mach_par
                
              case ( 'Mach_pol' )
                res = Mach_pol
                
              case ( 'V_sound' )
                res = Vsound / fact_time
                
              case ( 'V_neo' )
                res = Vneo / fact_time
                
              case ( 'Vperp_e' )
                res = Vperp_e / fact_time
                
              case ( 'Vperp_i' )
                res = Vperp_i / fact_time
                
              case ( 'V_ExB' )
                res = V_ExB / fact_time
                
              case ( 'Vstar_e' )
                res = Vstar_e / fact_time
                
              case ( 'Vstar_i' )
                res = Vstar_i / fact_time
                
              case ( 'ki_neo' )
                res = ki_neo
                
              case ( 'mu_neo' )
                res = mu_neo / fact_time
                
              case ( 'T_e' )
                res = Te0 * fact_T
                
              case ( 'T_i' )
                res = Ti0 * fact_T
                
              case ( 'E_||' )
                res = E_par / fact_time
                
              case ( 'E_crit' )
                res = E_crit / fact_time
                
              case ( 'E_dreicer' )
                res = E_dreicer / fact_time

              case ( 'theta_geo'    )
                res = theta_geo

              case ( 'bnd_normal_R' )
                res = nmlR

              case ( 'bnd_normal_Z' )
                res = nmlZ

              case ( 'Bnorm'        )
                res = Bnorm 

              case ( 'Btan'         )
                res = Btan

              case ( 'Jnorm'        )
                res = (JpolR*nmlR + JpolZ*nmlZ) / fact_mu_zero

              case ( 'gradP_norm' )
                res = (P0_R *nmlR +  P0_Z*nmlZ) / fact_mu_zero
 
              case ( 'Jpar'         )
                res = Jpar/fact_mu_zero

              case ( 'Jpar_ionsat'  )
                res = Jpar_ionsat * fact_jsat / fact_mu_zero

              case ( 'vpar_norm'    )
                res = vpar0 * Bnorm / fact_time 

              case ( 'vu_norm'   )
                res = (VR*nmlR + VZ*nmlZ - vpar0*Bnorm) / fact_time

              case ( 'vtot_norm'   )
                res = (VR*nmlR + VZ*nmlZ) / fact_time

              case ( 'heatF_sheath' )
                res = gamma_stangeby*r0*T0/2.d0*vpar0*Bnorm*fact_flux

              case ( 'heatF_par_cd' )
                res = ZKpar_flux * fact_flux

              case ( 'heatF_prp_cd' )
                res = ZKperp_flux * fact_flux

              case ( 'heatF_tot_cd' )
                res = (ZKperp_flux + ZKpar_flux) * fact_flux

              case ( 'heatF_par_cv' )
                res = pres_flux_par * fact_flux

              case ( 'heatF_prp_cv' )
                res = (pres_flux_tot-pres_flux_par) * fact_flux

              case ( 'heatF_tot_cv' )
                res = pres_flux_tot * fact_flux

              case ( 'heatF_total'  )
                res = (pres_flux_tot + ZKperp_flux + ZKpar_flux) * fact_flux

              case ( 'kinEn_F_par' )
                res = kin_flux_par * fact_flux

              case ( 'kinEn_F_perp ' )
                res = (kin_flux_tot-kin_flux_par) * fact_flux

              case ( 'kinEn_F_tot ' )
                res = kin_flux_tot * fact_flux

              case ( 'partF_par_cd' )
                res = Dpar_flux * fact_ne / fact_time

              case ( 'partF_prp_cd' )
                res = Dperp_flux * fact_ne / fact_time

              case ( 'partF_par_cv' )
                res = partF_cnv_par * fact_ne / fact_time

              case ( 'partF_prp_cv' )
                res = (partF_cnv_tot - partF_cnv_par) * fact_ne / fact_time

              case ( 'partF_total'  )
                res = partF_cnv_tot * fact_ne / fact_time

              case ( 'npartF_total'  )
                res = neut_part_flux * fact_ne / fact_time

              case ( 'ExB_norm'  )
                res = ExB_norm * fact_flux
  
              case ( 'J_bootstrap' )
                res = J_boot / R / fact_mu_zero

#if (defined WITH_Neutrals) && (!defined WITH_Impurities)
              case ( 'radiation' )

                if (rn0 .lt. 0.d0) then
                  res = r0 * fact_ne * r0 * fact_ne * LradDcont_T &
                       + r0 * fact_ne * frad_bg
                else
                  res = r0 * fact_ne * rn0 * fact_ne * LradDrays_T &
                       + r0 * fact_ne * r0 * fact_ne * LradDcont_T &
                       + r0 * fact_ne * frad_bg
                endif

              case ( 'brem' )
                res = r0 * fact_ne * r0 * fact_ne * LradDcont_T
#endif
#ifdef WITH_Impurities
              case ( 'radiation' )
                res = (r0_corr + beta_imp*rn0_corr) * rn0_corr * Lrad * fact_rad
#endif

              case default
                write(*,*) 'ERROR in '//trim(THIS_ROUTINE_NAME)//': Illegal expression ("' //      &
                  trim(expr_list%expr(iexpr)%name) // '")'
                ierr = 100
                return
                
            end select
            
            result(itorpos, ipolpos, jpolpos, iexpr) = res
            
          end do loop_expr
        end do loop_tor
      end do loop_pol2
    end do loop_pol1
    
  end subroutine eval_expr
  
  
  
  
  
end module mod_expression
