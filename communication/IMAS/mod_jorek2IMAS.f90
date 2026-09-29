!>  Module that contais functions to fill IDSs from JOREK data
module mod_jorek2IMAS

#ifdef USE_IMAS
  use ids_schemas !, only: ids_equilibrium
  use ids_routines

  use mod_parameters 
  use mod_new_diag
  use data_structure
  use nodes_elements
  use constants
  use mod_expression, only: exprs, exprs_all_int
  use exec_commands,  only: average, expr_list, clean_up, step_imported, qprofile, &
                            zeroD_quantities, separatrix, rectangle, boundary_quantities
  use parse_commands, only: type_command
  use settings,       only: set_setting
  
  implicit none
  
 
  public
 
  ! Transform to COCOS convention
  real*8 :: fact_psi = -2.d0 * PI ! -2pi for COCOS 17 in DD4, 2pi for COCOS 11 in DD3
  real*8 :: fact_Ip  = -1.d0
  real*8 :: fact_phi_dir = -1.d0  ! -1 to go from 8-->11 or 8-->17

  ! Structure for rectangular grid parameters
  type :: t_rect_grid_params
    integer :: nR = 70         ! Number of points in the major radius direction
    integer :: nZ = 120        ! Number of points in the vertical direction
    real*8  :: R_min =  3.d0   ! Minimum R
    real*8  :: R_max = 10.d0   ! Maximum R
    real*8  :: Z_min = -6.d0   ! Minimum Z
    real*8  :: Z_max =  6.d0   ! Maximum Z
  end type t_rect_grid_params


  ! Strucure for PFC thin wall grid, discretized with triangles
  type :: t_PFC_triang_grid
    integer :: n_nodes
    integer :: n_triangles
    integer :: n_pol
    integer :: n_phi
    real*8, allocatable  :: tria_xyz(:,:)
    integer, allocatable :: tria_connectivity(:,:)
  end type t_PFC_triang_grid

  type(t_expr_list), save :: expr_list_bnd

  
  ! *******************************************************************************************************
  ! *************** Data structures needed to read geometry of the STARWALL coils *************************
  ! *******************************************************************************************************
  ! --- Constants for STARWALL coils 
  character(len=12), parameter    :: AXISYM_THICK = 'axisym_thick' !< Axisymmetric coils with dR,dZ extent
  character(len=12), parameter    :: AXISYM_FILA  = 'axisym_fila ' !< Axisymmetric coils specified by filaments
  character(len=12), parameter    :: GENERAL_THIN = 'general_thin' !< 3D coils specified by a list of points
  integer, parameter              :: N_MAX_FILA  = 2001       !< Maximum number of filaments per c
  integer, parameter              :: N_MAX_PARTS = 2001       !< Maximum number of parts per coil 
  integer, parameter              :: N_MAX_PTS   = 2001       !< Maximum number of points per coil

  ! --- Data structure for a single STARWALL coil of arbitrary type
  type :: t_coil_starwall
    ! --- For all coil types
    character(len=128)   :: name     = 'NONE'        !< Name of the coil
    character(len=32)    :: coil_type                !< Coil type according to constants specified above
    real*8               :: resist   = -999.         !< Total coil resistivity in Ohm
    real*8               :: nturns   = -1d4          !< Total number of turns in the coil
    integer              :: ntorpts  = 0             !< Number of toroidal points (axisym coils)
    ! --- For axisymmetric coils with dR and dZ extent
    real*8               :: R(N_MAX_PARTS)        = -999.      !< Position in R (center of coil)
    real*8               :: Z(N_MAX_PARTS)        = -999.      !< Position in Z (center of coil)
    real*8               :: dR(N_MAX_PARTS)       = -999.      !< Full width in R
    real*8               :: dZ(N_MAX_PARTS)       = -999.      !< Full width in Z
    integer              :: nbands_R = 0                       !< #radial bands for each coil part
    integer              :: nbands_Z = 0                       !< #toroidal bands for each coil part
    real*8               :: n_thick_turns(N_MAX_PARTS) = -1.d4 !< Number of turns in each coil part
    integer              :: nparts_coil  = 1                  !< Number of bands
    ! --- For coils specified by filaments
    integer              :: n_fila   = 0                  !< Number of filaments for the coil
    real*8               :: n_fila_turns(N_MAX_FILA) = 1. !< Number of turns in each filament
    real*8               :: R_fila(N_MAX_FILA) = -999.    !< R-position of filaments
    real*8               :: Z_fila(N_MAX_FILA) = -999.    !< Z-position of filaments
    ! --- For 3D coil specified by a list of points
    integer              :: n_pts    = 0             !< Number of points along the coil
    real*8               :: width    = -999.         !< Width of coil
    real*8               :: xpts(N_MAX_PTS) = -999.  !< x-Position of points
    real*8               :: ypts(N_MAX_PTS) = -999.  !< y-Position of points
    real*8               :: zpts(N_MAX_PTS) = -999.  !< z-Position of points
    real*8               :: dir3d(3*N_MAX_PTS) = -999. !< optional, band direction for each point
  end type t_coil_starwall

  ! --- Data structure for a coil set
  type :: t_coil_set_starwall
    character(len=32)    :: name         = 'NONE'    !< Name of coil set (should not be changed by user)
    character(len=256)   :: description  = 'NONE'    !< Human readable description
    integer              :: ncoil        = 0         !< Number of coils in the set
    integer              :: index_start  = -999      !< Index of first coil of this set in all_coil_set
    type(t_coil_starwall), allocatable :: coil(:)             !< Coils in this coil set
  end type t_coil_set_starwall
  ! *******************************************************************************************************


  contains


  subroutine fill_profiles_w_JOREK_var(first_step, time_SI, plasma_profiles_ids)  

    use phys_module, only : F0, central_density, sqrt_mu0_rho0, &
                           sqrt_mu0_over_rho0, central_mass

    implicit none

    ! --- External parameters
    logical,   intent(in) :: first_step   ! is this the first step?
    real*8,    intent(in) :: time_SI
    
    type(ids_plasma_profiles), target,     intent(inout) :: plasma_profiles_ids

   
    ! --- Local parameters 
    integer    :: i, j, k, m, etype, irst, int, i_var, i_tor, index, index_node, my_id, ierr
    real*8     :: fact_T, fact_v, fact_zj, rho0, fact_phi, fact_rho, fact_w, fact_nimp
    
    
    ! **********************************************************************************
    ! ******************************* IMAS **********************************************
    ! **********************************************************************************
    type(ids_generic_grid_scalar),            pointer :: ggd_scalar
    type(ids_generic_grid_vector_components), pointer :: ggd_vector
    type(ids_generic_grid_aos3_root),         pointer :: grid
    
    integer:: num_nodes
    
    integer :: n_slice, i_slice, grid_ind, grid_sub_ind, n_grid_sub, n_grid, i_ion, a_elm, z_elm, i_neut
    ! **********************************************************************************
  
    ! --- Number of grids and grid subsets
    n_grid       = 1
    n_grid_sub   = 1
    grid_ind     = 1  ! Index
    grid_sub_ind = 1  ! Index
  
    if (first_step) then
      ! --- Put the grid in GGD
      allocate( plasma_profiles_ids%grid_ggd(n_grid) )
      grid => plasma_profiles_ids%grid_ggd(grid_ind)
      grid%time = time_SI
      call grid2ggd( grid, node_list, element_list, bnd_node_list, bnd_elm_list )
    else 
      if ( associated(plasma_profiles_ids%grid_ggd) ) then
        call ids_deallocate_struct(plasma_profiles_ids%grid_ggd(grid_ind), .false.)
        deallocate(plasma_profiles_ids%grid_ggd)
      endif
    endif

    ! --- Normalization factors for IMAS
    rho0               = central_density * 1.d20 * central_mass * ATOMIC_MASS_UNIT
    sqrt_mu0_rho0      = sqrt( mu_zero * rho0 )
    sqrt_mu0_over_rho0 = sqrt( mu_zero / rho0 )
  
    fact_v    =  1.d0 /  sqrt_mu0_rho0 
    fact_w    = -1.d0 /  sqrt_mu0_rho0      ! Transform for COCOS convention of toroidal direction (anti-clockwise) 
    fact_phi  = -1.d0 /  sqrt_mu0_rho0 * F0 ! COCOS convection: F0 depends on phi direction
    fact_zj   = -1.d0 / mu_zero * fact_Ip   ! Last sign due to COCOS transformation
    fact_rho  =  rho0 
    fact_T    =  1.d0 / ( EL_CHG * mu_zero * central_density * 1.d20 )  
  
    ! --- Set times
    n_slice = 1  
    i_slice = 1
    allocate(  plasma_profiles_ids%time(n_slice) )
    allocate(  plasma_profiles_ids%ggd(n_slice ) )

    plasma_profiles_ids%ids_properties%homogeneous_time = IDS_TIME_MODE_HETEROGENEOUS
  
    plasma_profiles_ids%time(i_slice)     = time_SI 
    plasma_profiles_ids%ggd(i_slice)%time = time_SI

    ! --- Fill MHD data
    do i=1, n_var 
  
      ! --- Poloidal magnetic flux
      if (variable_names(i) == 'Psi') then      
        allocate( plasma_profiles_ids%ggd(i_slice)%psi(n_grid_sub))
        ggd_scalar => plasma_profiles_ids%ggd(i_slice)%psi(grid_sub_ind)
        call fill_Bezier_coefficients( ggd_scalar, node_list, var_psi, grid_ind, grid_sub_ind, fact_psi )
      endif
  
      ! --- Electrostatic potential 
      if (variable_names(i) == 'u') then      
        allocate( plasma_profiles_ids%ggd(i_slice)%phi_potential(n_grid_sub))
        ggd_scalar => plasma_profiles_ids%ggd(i_slice)%phi_potential(grid_sub_ind)
        call fill_Bezier_coefficients( ggd_scalar, node_list, var_u, grid_ind, grid_sub_ind, fact_phi )
      endif
  
      ! --- Toroidal current density * R
      if (variable_names(i) == 'zj') then      
        allocate( plasma_profiles_ids%ggd(i_slice)%r_j_total_phi(n_grid_sub))
        ggd_scalar => plasma_profiles_ids%ggd(i_slice)%r_j_total_phi(grid_sub_ind)
        call fill_Bezier_coefficients( ggd_scalar, node_list, var_zj, grid_ind, grid_sub_ind, fact_zj )
      endif
  
      ! --- Toroidal vorticity / R 
      if (variable_names(i) == 'omega') then      
        allocate( plasma_profiles_ids%ggd(i_slice)%vorticity_over_r(n_grid_sub))
        ggd_vector => plasma_profiles_ids%ggd(i_slice)%vorticity_over_r(grid_sub_ind)
        call fill_Bezier_vector_coefficients( ggd_vector, node_list, var_w, grid_ind, grid_sub_ind, fact_w, 'phi')
      endif
  
      ! --- Mass density
      if ((variable_names(i) == 'rho') .and. with_rho) then      
        allocate( plasma_profiles_ids%ggd(i_slice)%mass_density(n_grid_sub))
        ggd_scalar => plasma_profiles_ids%ggd(i_slice)%mass_density(grid_sub_ind)
        call fill_Bezier_coefficients( ggd_scalar, node_list, var_rho, grid_ind, grid_sub_ind, fact_rho )
      endif
  
      ! --- Total temperature 
      if (variable_names(i) == 'T') then      
        ! --- Te
        allocate( plasma_profiles_ids%ggd(i_slice)%electrons%temperature(n_grid_sub))
        ggd_scalar => plasma_profiles_ids%ggd(i_slice)%electrons%temperature(grid_sub_ind)
        call fill_Bezier_coefficients( ggd_scalar, node_list, var_T, grid_ind, grid_sub_ind, fact_T*0.5d0 )
  
        ! --- Ti
        allocate( plasma_profiles_ids%ggd(i_slice)%t_i_average(n_grid_sub))
        ggd_scalar => plasma_profiles_ids%ggd(i_slice)%t_i_average(grid_sub_ind)
        call fill_Bezier_coefficients( ggd_scalar, node_list, var_T, grid_ind, grid_sub_ind, fact_T*0.5d0 )
      endif
  
      ! --- Ion temperature
      if (variable_names(i) == 'T_i') then      
        allocate( plasma_profiles_ids%ggd(i_slice)%t_i_average(n_grid_sub))
        ggd_scalar => plasma_profiles_ids%ggd(i_slice)%t_i_average(grid_sub_ind)
        call fill_Bezier_coefficients( ggd_scalar, node_list, var_Ti, grid_ind, grid_sub_ind, fact_T )
      endif
  
      ! --- Electron temperature
      if (variable_names(i) == 'T_e') then      
        allocate( plasma_profiles_ids%ggd(i_slice)%electrons%temperature(n_grid_sub))
        ggd_scalar => plasma_profiles_ids%ggd(i_slice)%electrons%temperature(grid_sub_ind)
        call fill_Bezier_coefficients( ggd_scalar, node_list, var_Te, grid_ind, grid_sub_ind, fact_T )
      endif
  
      ! --- Parallel velocity
      if (variable_names(i) == 'v_par') then      
        allocate( plasma_profiles_ids%ggd(i_slice)%velocity_parallel_over_b_field(n_grid_sub))
        ggd_scalar => plasma_profiles_ids%ggd(i_slice)%velocity_parallel_over_b_field(grid_sub_ind)
        call fill_Bezier_coefficients( ggd_scalar, node_list, var_vpar, grid_ind, grid_sub_ind, fact_v )
      endif

      ! --- Impurity density
      if (variable_names(i) == 'rho_imp') then      
        i_ion = 1
        allocate(plasma_profiles_ids%ggd(i_slice)%ion(i_ion))
        allocate(plasma_profiles_ids%ggd(i_slice)%ion(i_ion)%element(1))

        ! --- Element identifier
        call get_element_atomic_numbers(imp_type(index_main_imp), a_elm, z_elm )
        plasma_profiles_ids%ggd(i_slice)%ion(i_ion)%element(1)%a   = a_elm
        plasma_profiles_ids%ggd(i_slice)%ion(i_ion)%element(1)%z_n = z_elm

        fact_nimp =  1.d20 * central_density * central_mass / a_elm

        ! --- Density in GGD
        allocate(plasma_profiles_ids%ggd(i_slice)%ion(i_ion)%density(n_grid_sub))
        ggd_scalar => plasma_profiles_ids%ggd(i_slice)%ion(i_ion)%density(grid_sub_ind)
        call fill_Bezier_coefficients( ggd_scalar, node_list, var_rhoimp, grid_ind, grid_sub_ind, fact_nimp )
      endif

      ! --- Main ion neutrals
      if (variable_names(i) == 'rho_n') then      
        i_neut = 1
        allocate(plasma_profiles_ids%ggd(i_slice)%neutral(i_neut))
        allocate(plasma_profiles_ids%ggd(i_slice)%neutral(i_neut)%element(1))

        plasma_profiles_ids%ggd(i_slice)%neutral(i_neut)%element(1)%a   = nint(central_mass)
        plasma_profiles_ids%ggd(i_slice)%neutral(i_neut)%element(1)%z_n = 1

        ! --- Density in GGD
        allocate(plasma_profiles_ids%ggd(i_slice)%neutral(i_neut)%density(n_grid_sub))
        ggd_scalar => plasma_profiles_ids%ggd(i_slice)%neutral(i_neut)%density(grid_sub_ind)
        call fill_Bezier_coefficients( ggd_scalar, node_list, var_rhon, grid_ind, grid_sub_ind, central_density*1d20 )
      endif
  
    enddo

  end subroutine fill_profiles_w_JOREK_var




  ! --- Fill a wall IDS from STARWALL data, needs to be adapted to CARIDDI!
  subroutine fill_wall_IDS(first_step, time_SI, wall_thickness, wall_ids, res_bnd, PFC_triang_grid)  

    use phys_module, only : F0, central_density, sqrt_mu0_rho0, &
                           sqrt_mu0_over_rho0, central_mass, imp_type, &
                           gamma, index_main_imp
    use vacuum
    use vacuum_response,  only: reconstruct_triangle_potentials

    implicit none

    ! --- External parameters
    logical,                 intent(in)    :: first_step      ! is this the first step?
    real*8,                  intent(in)    :: time_SI
    real*8,                  intent(in)    :: wall_thickness  !< effective thin wall thickness
    type(ids_wall),  target, intent(inout) :: wall_ids
    real*8,                  intent(in)    :: res_bnd(:,:,:)   !< data for PFC grid
    type(t_PFC_triang_grid), intent(in)    :: PFC_triang_grid  !< grid for PFC wall, discretized with triangles
   
    ! --- Local parameters 
    integer    :: i, j, k, m, var_rad, i_var, i_tor, index, index_node, my_id=0, ierr, i_tri
    integer    :: i_phi, i_pol, n_phi, n_pol, ipol1, ipol2, ipol3, ipol4, itor1, itor2, itor3, itor4
    integer    :: i1, i2, i3, i4, i_exp
    real*8     :: rho0, fact_rad, co, si, Rtri
    real*8     :: phi1, phi2, phi3, r1(3), r2(3), r3(3), r21(3), r32(3), r13(3), r21_cross_r32(3)
    real*8     :: j_lin(3), r_mid(3), Iw_net_tor
    real*8, allocatable :: tripot_w(:)
    
    ! **********************************************************************************
    ! ******************************* IMAS **********************************************
    ! **********************************************************************************
    type(ids_generic_grid_scalar),      pointer :: ggd_scalar
    type(ids_generic_grid_aos3_root),   pointer :: grid
    
    integer:: n_wall_nodes, n_wall_triangles, n_walls, i_vv, i_fw
    
    integer :: n_slice, i_slice, grid_ind, grid_sub_ind, n_grid_sub, n_grid
    ! **********************************************************************************
  
    ! --- Number of grids and grid subsets
    n_grid       = 1
    n_grid_sub   = 1
    grid_ind     = 1  ! Index
    grid_sub_ind = 1  ! Index    
    n_walls      = 2  ! 1 for vessel current, 2 for first wall fluxes
    i_vv         = 1
    i_fw         = 2

    wall_ids%ids_properties%homogeneous_time = IDS_TIME_MODE_HETEROGENEOUS

    ! --- Set times
    n_slice = 1  
    i_slice = 1
    allocate(  wall_ids%time(n_slice) )

    n_wall_triangles = sr%ntri_w

    if (first_step) then
      allocate( wall_ids%description_ggd(n_walls))
    endif

    ! *******************************************************************************
    ! ************* Export STARWALL wall currents to wall_IDS ***********************
    ! *******************************************************************************
    if (first_step) then

      ! --- Put the wall grid in GGD
      allocate( wall_ids%description_ggd(i_vv)%grid_ggd(n_grid) )

      wall_ids%description_ggd(i_vv)%type%index  = 2  ! For thin wall description
      
      grid => wall_ids%description_ggd(i_vv)%grid_ggd(grid_ind)
      grid%time = time_SI
      
      call triang_thin_wall_2ggd( grid, sr%xyzpot_w, sr%jpot_w, subset_choice="triangles")

      ! --- Save wall thickness
      allocate(wall_ids%description_ggd(i_vv)%thickness(n_grid))
      allocate(wall_ids%description_ggd(i_vv)%thickness(1)%grid_subset(n_grid_sub))
      allocate(wall_ids%description_ggd(i_vv)%thickness(1)%grid_subset(1)%values(n_wall_triangles))

      wall_ids%description_ggd(i_vv)%thickness(1)%grid_subset(1)%grid_index        = grid_ind
      wall_ids%description_ggd(i_vv)%thickness(1)%grid_subset(1)%grid_subset_index = grid_sub_ind

      wall_ids%description_ggd(i_vv)%thickness(1)%grid_subset(1)%values(:) = wall_thickness
      wall_ids%description_ggd(i_vv)%thickness(1)%time = time_SI

      !--- Information about the wall component
      allocate(wall_ids%description_ggd(i_vv)%component(n_grid))
      allocate(wall_ids%description_ggd(i_vv)%component(1)%type(1))
      allocate(wall_ids%description_ggd(i_vv)%component(1)%identifiers(1))
      wall_ids%description_ggd(i_vv)%component(1)%identifiers = ''
      wall_ids%description_ggd(i_vv)%component(1)%type(1)%grid_index        = grid_ind
      wall_ids%description_ggd(i_vv)%component(1)%type(1)%grid_subset_index = grid_sub_ind
      wall_ids%description_ggd(i_vv)%component(1)%type(1)%identifier%index  = 0  ! --- 0 means not specified yet
      allocate(wall_ids%description_ggd(i_vv)%component(1)%type(1)%identifier%name(1))
      wall_ids%description_ggd(i_vv)%component(1)%type(1)%identifier%name = "Vacuum Vessel"
      allocate(wall_ids%description_ggd(i_vv)%component(1)%type(1)%identifier%description(1))
      wall_ids%description_ggd(i_vv)%component(1)%type(1)%identifier%description = "A single layer of the &
          vacuum vessel discretized with linear thin triangles"
      wall_ids%description_ggd(i_vv)%component(1)%time = time_SI
    else
      if ( associated(wall_ids%description_ggd(i_vv)%grid_ggd)) then
        call ids_deallocate_struct(wall_ids%description_ggd(i_vv)%grid_ggd(grid_ind), .false.)
        deallocate( wall_ids%description_ggd(i_vv)%thickness)
        deallocate( wall_ids%description_ggd(i_vv)%component)
        deallocate(wall_ids%description_ggd(i_vv)%grid_ggd)
      end if
    endif

    ! --- Set times
    allocate( wall_ids%description_ggd(i_vv)%ggd(n_slice) )
  
    wall_ids%time(i_slice) = time_SI 
    wall_ids%description_ggd(i_vv)%ggd(i_slice)%time = time_SI
    
    allocate( wall_ids%description_ggd(i_vv)%ggd(i_slice)%j_total(n_grid_sub) )
    allocate( wall_ids%description_ggd(i_vv)%ggd(i_slice)%j_total(1)%r(n_wall_triangles) ) ! --- one value per triangle
    allocate( wall_ids%description_ggd(i_vv)%ggd(i_slice)%j_total(1)%z(n_wall_triangles) ) ! --- one value per triangle
    allocate( wall_ids%description_ggd(i_vv)%ggd(i_slice)%j_total(1)%phi(n_wall_triangles) ) ! --- one value per triangle

    wall_ids%description_ggd(i_vv)%ggd(i_slice)%j_total(1)%grid_index        = 1
    wall_ids%description_ggd(i_vv)%ggd(i_slice)%j_total(1)%grid_subset_index = 1

    ! --- Get the current density per triangle
    ! --- The triangle current density
    call reconstruct_triangle_potentials(tripot_w, wall_curr, my_id, Iw_net_tor)

    do i = 1, n_wall_triangles

      ! --- Wall potential at triangle nodes
      phi1   = tripot_w(sr%jpot_w(i,1)) + Iw_net_tor*sr%phi0_w(i,1) 
      phi2   = tripot_w(sr%jpot_w(i,2)) + Iw_net_tor*sr%phi0_w(i,2) 
      phi3   = tripot_w(sr%jpot_w(i,3)) + Iw_net_tor*sr%phi0_w(i,3) 

      ! --- Position of triangle nodes
      r1(:)  = sr%xyzpot_w(sr%jpot_w(i,1),:)
      r2(:)  = sr%xyzpot_w(sr%jpot_w(i,2),:)
      r3(:)  = sr%xyzpot_w(sr%jpot_w(i,3),:)
      r21(:) = r1(:)-r2(:)
      r32(:) = r2(:)-r3(:)
      r21_cross_r32(:) = (/ r21(2)*r32(3) - r21(3)*r32(2), r21(3)*r32(1) - r21(1)*r32(3),          &
        r21(1)*r32(2) - r21(2)*r32(1) /)

      j_lin(:) = ( phi1*(r3-r2)+phi2*(r1-r3)+phi3*(r2-r1) ) / sqrt(sum(r21_cross_r32**2) ) / mu_zero 

      r_mid(:) = (r1 + r2 + r3) / 3.d0   ! Middle point of the triangle

      Rtri = sqrt( r_mid(1)**2.d0 + r_mid(2)**2.d0 )
      co   =  r_mid(1) / Rtri
      si   = -r_mid(2) / Rtri

      wall_ids%description_ggd(i_vv)%ggd(i_slice)%j_total(1)%r(i)        =  ( j_lin(1)*co - j_lin(2)*si ) / wall_thickness
      wall_ids%description_ggd(i_vv)%ggd(i_slice)%j_total(1)%z(i)        =    j_lin(3)                    / wall_thickness
      wall_ids%description_ggd(i_vv)%ggd(i_slice)%j_total(1)%phi(i)      =  (-j_lin(1)*si - j_lin(2)*co ) / wall_thickness * fact_Ip   
    enddo

    if (first_step) then
      allocate( wall_ids%description_ggd(i_vv)%ggd(i_slice)%resistivity(n_grid_sub) )
      allocate( wall_ids%description_ggd(i_vv)%ggd(i_slice)%resistivity(1)%values(n_wall_triangles) ) 
      wall_ids%description_ggd(i_vv)%ggd(i_slice)%resistivity(1)%values(:) = sr%eta_thin_w * wall_thickness * wall_resistivity_fact

      wall_ids%description_ggd(i_vv)%ggd(i_slice)%resistivity(1)%grid_index        = 1
      wall_ids%description_ggd(i_vv)%ggd(i_slice)%resistivity(1)%grid_subset_index = 1
    endif
    ! *******************************************************************************
   

    ! *******************************************************************************
    ! ************* Export FW and divertor fluxes and currents **********************
    ! *******************************************************************************
    ! --- Export grid boundary geometry as thin triangles (FW + divertor)
    if (first_step) then

      ! --- Put the wall grid in GGD
      allocate( wall_ids%description_ggd(i_fw)%grid_ggd(n_grid) )

      wall_ids%description_ggd(i_fw)%type%index  = 2  ! For thin wall description
      
      grid => wall_ids%description_ggd(i_fw)%grid_ggd(grid_ind)
      grid%time = time_SI
      
      call triang_thin_wall_2ggd( grid, PFC_triang_grid%tria_xyz, PFC_triang_grid%tria_connectivity, subset_choice="nodes")

      !--- Information about the wall component
      allocate(wall_ids%description_ggd(i_fw)%component(n_grid))
      allocate(wall_ids%description_ggd(i_fw)%component(1)%type(1))
      allocate(wall_ids%description_ggd(i_fw)%component(1)%identifiers(1))
      wall_ids%description_ggd(i_fw)%component(1)%identifiers= ''
      wall_ids%description_ggd(i_fw)%component(1)%type(1)%grid_index        = grid_ind
      wall_ids%description_ggd(i_fw)%component(1)%type(1)%grid_subset_index = grid_sub_ind
      wall_ids%description_ggd(i_fw)%component(1)%type(1)%identifier%index  = 0  ! --- 0 means not specified yet
      allocate(wall_ids%description_ggd(i_fw)%component(1)%type(1)%identifier%name(1))
      wall_ids%description_ggd(i_fw)%component(1)%type(1)%identifier%name = "FW + divertor"
      allocate(wall_ids%description_ggd(i_fw)%component(1)%type(1)%identifier%description(1))
      wall_ids%description_ggd(i_fw)%component(1)%type(1)%identifier%description = "A single layer representing the first wall &
          + divertor surfaces discretized with linear thin triangles"
      wall_ids%description_ggd(i_fw)%component(1)%time = time_SI
    else   
      if ( associated(wall_ids%description_ggd(i_fw)%grid_ggd)) then
        call ids_deallocate_struct(wall_ids%description_ggd(i_fw)%grid_ggd(grid_ind), .false.)
        deallocate( wall_ids%description_ggd(i_fw)%component)
        deallocate(wall_ids%description_ggd(i_fw)%grid_ggd)    
      end if                                                     
    endif

    allocate( wall_ids%description_ggd(i_fw)%ggd(n_slice) )
  
    wall_ids%time(i_slice) = time_SI 
    wall_ids%description_ggd(i_fw)%ggd(i_slice)%time = time_SI

    ! --- Fill expressions in the wall nodes
    allocate( wall_ids%description_ggd(i_fw)%ggd(i_slice)%j_total(n_grid_sub) )
    wall_ids%description_ggd(i_fw)%ggd(i_slice)%j_total(1)%grid_index        = 1
    wall_ids%description_ggd(i_fw)%ggd(i_slice)%j_total(1)%grid_subset_index = 1

    n_pol = PFC_triang_grid%n_pol
    n_phi = PFC_triang_grid%n_phi
    n_wall_nodes = PFC_triang_grid%n_nodes

    do i_exp=1, expr_list_bnd%n_expr
      
      ! --- Total heatflux
      if (expr_list_bnd%expr(i_exp)%name=='heatF_total') then
        allocate( wall_ids%description_ggd(i_fw)%ggd(i_slice)%power_density(n_grid_sub) )
        allocate( wall_ids%description_ggd(i_fw)%ggd(i_slice)%power_density(1)%values(n_wall_nodes) ) ! --- one value per node

        wall_ids%description_ggd(i_fw)%ggd(i_slice)%power_density(1)%grid_index        = 1
        wall_ids%description_ggd(i_fw)%ggd(i_slice)%power_density(1)%grid_subset_index = 1

        do i_phi=1, n_phi
          do i_pol=1, n_pol 
            i = i_pol + (i_phi-1)*n_pol  !< Global index of refence node
            wall_ids%description_ggd(i_fw)%ggd(i_slice)%power_density(1)%values(i) = res_bnd(i_phi,i_pol,i_exp)
          enddo
        enddo
      endif

      ! --- Psi
      if (expr_list_bnd%expr(i_exp)%name=='Psi') then
        allocate( wall_ids%description_ggd(i_fw)%ggd(i_slice)%psi(n_grid_sub) )
        allocate( wall_ids%description_ggd(i_fw)%ggd(i_slice)%psi(1)%values(n_wall_nodes) ) ! --- one value per node

        wall_ids%description_ggd(i_fw)%ggd(i_slice)%psi(1)%grid_index        = 1
        wall_ids%description_ggd(i_fw)%ggd(i_slice)%psi(1)%grid_subset_index = 1

        do i_phi=1, n_phi
          do i_pol=1, n_pol 
            i = i_pol + (i_phi-1)*n_pol  !< Global index of refence node
            wall_ids%description_ggd(i_fw)%ggd(i_slice)%psi(1)%values(i) = res_bnd(i_phi,i_pol,i_exp)
          enddo
        enddo
      endif

      ! --- J_phi
      if (expr_list_bnd%expr(i_exp)%name=='Jtor') then
        allocate( wall_ids%description_ggd(i_fw)%ggd(i_slice)%j_total(1)%phi(n_wall_nodes) ) ! --- one value per node

        do i_phi=1, n_phi
          do i_pol=1, n_pol 
            i = i_pol + (i_phi-1)*n_pol  !< Global index of refence node
            wall_ids%description_ggd(i_fw)%ggd(i_slice)%j_total(1)%phi(i) = res_bnd(i_phi,i_pol,i_exp) * fact_Ip
          enddo
        enddo
      endif

      ! --- J_R
      if (expr_list_bnd%expr(i_exp)%name=='JR') then
        allocate( wall_ids%description_ggd(i_fw)%ggd(i_slice)%j_total(1)%r(n_wall_nodes) ) ! --- one value per node

        do i_phi=1, n_phi
          do i_pol=1, n_pol 
            i = i_pol + (i_phi-1)*n_pol  !< Global index of refence node
            wall_ids%description_ggd(i_fw)%ggd(i_slice)%j_total(1)%r(i) = res_bnd(i_phi,i_pol,i_exp) 
          enddo
        enddo
      endif

      ! --- J_Z
      if (expr_list_bnd%expr(i_exp)%name=='JZ') then
        allocate( wall_ids%description_ggd(i_fw)%ggd(i_slice)%j_total(1)%z(n_wall_nodes) ) ! --- one value per node

        do i_phi=1, n_phi
          do i_pol=1, n_pol 
            i = i_pol + (i_phi-1)*n_pol  !< Global index of refence node
            wall_ids%description_ggd(i_fw)%ggd(i_slice)%j_total(1)%z(i) = res_bnd(i_phi,i_pol,i_exp) 
          enddo
        enddo
      endif
    
    enddo
    ! *******************************************************************************
  
  end subroutine fill_wall_IDS
 




  


  ! --- Fill a plasma transport IDS with boundary heat fluxes
  subroutine fill_transport_IDS(first_step, time_SI, transport_ids, res_bnd, PFC_triang_grid)

    use phys_module, only: central_mass

    implicit none

    ! --- External parameters
    logical,                 intent(in)    :: first_step      ! is this the first step?
    real*8,                  intent(in)    :: time_SI
    real*8,                  intent(in)    :: res_bnd(:,:,:)   !< data for PFC grid
    type(t_PFC_triang_grid), intent(in)    :: PFC_triang_grid  !< grid for PFC wall, discretized with triangles
    type(ids_plasma_transport),  target, intent(inout) :: transport_ids
   
    ! --- Local parameters 
    integer    :: i, j, k, m, i_tor, index, index_node, my_id=0, ierr, i_tri
    integer    :: i_phi, i_pol, n_phi, n_pol, n_wall_nodes
    integer    :: i_exp, i_qpar_e, i_qpar_i, i_partf
    real*8     :: qpar_e, qpar_i
    real*8     :: partf

    ! **********************************************************************************
    ! ******************************* IMAS **********************************************
    ! **********************************************************************************
    type(ids_generic_grid_scalar),      pointer :: ggd_scalar
    type(ids_generic_grid_aos3_root),   pointer :: grid
    
    integer :: n_slice, i_slice, grid_ind, grid_sub_ind, n_grid_sub, n_grid, i_model
    ! **********************************************************************************
  
    ! --- Number of grids and grid subsets
    n_grid       = 1
    n_grid_sub   = 1
    grid_ind     = 1  ! Index
    grid_sub_ind = 1  ! Index
    i_model      = 1
    n_slice      = 1
    i_slice      = 1

    transport_ids%ids_properties%homogeneous_time = IDS_TIME_MODE_HETEROGENEOUS
    allocate(transport_ids%time(n_slice))

    ! *******************************************************************************
    ! ************* Export FW and divertor fluxes and currents **********************
    ! *******************************************************************************
    ! --- Export grid boundary geometry as thin triangles (FW + divertor)
    if (first_step) then
      ! --- Put the wall grid in GGD
      allocate( transport_ids%grid_ggd(n_grid) )
      grid => transport_ids%grid_ggd(grid_ind)
      grid%time = time_SI
    
      call triang_thin_wall_2ggd( grid, PFC_triang_grid%tria_xyz, PFC_triang_grid%tria_connectivity, subset_choice="nodes")
    else
      if ( associated(transport_ids%grid_ggd)) then
        call ids_deallocate_struct(transport_ids%grid_ggd(grid_ind), .false.)
        deallocate(transport_ids%grid_ggd)
      end if
    endif

    allocate( transport_ids%model(1) )
    allocate( transport_ids%model(i_model)%ggd(n_slice) )
  
    transport_ids%time(i_slice) = time_SI
    transport_ids%model(i_model)%ggd(i_slice)%time = time_SI

    ! --- Fill expressions in the wall nodes
    n_pol = PFC_triang_grid%n_pol
    n_phi = PFC_triang_grid%n_phi
    n_wall_nodes = PFC_triang_grid%n_nodes

    i_qpar_e   = 0
    i_qpar_i   = 0
    i_partf    = 0
    do i_exp=1, expr_list_bnd%n_expr
      select case (trim(expr_list_bnd%expr(i_exp)%name))
      case ('qpar_e')
        i_qpar_e = i_exp
      case ('qpar_i')
        i_qpar_i = i_exp
      case ('partF_par')
        i_partf = i_exp
      end select
    enddo

    if (i_qpar_e == 0 .or. i_qpar_i == 0 .or. i_partf == 0) then
      write(*,*) 'ERROR: heat and particle flux expressions are missing for the plasma_transport IDS'
      stop
    endif

    allocate(transport_ids%model(i_model)%ggd(i_slice)%electrons%energy%flux_parallel(n_grid_sub))
    allocate(transport_ids%model(i_model)%ggd(i_slice)%electrons%energy%flux_parallel(1)%values(n_wall_nodes))
    transport_ids%model(i_model)%ggd(i_slice)%electrons%energy%flux_parallel(1)%grid_index = grid_ind
    transport_ids%model(i_model)%ggd(i_slice)%electrons%energy%flux_parallel(1)%grid_subset_index = grid_sub_ind

    allocate(transport_ids%model(i_model)%ggd(i_slice)%total_ion_energy%flux_parallel(n_grid_sub))
    allocate(transport_ids%model(i_model)%ggd(i_slice)%total_ion_energy%flux_parallel(1)%values(n_wall_nodes))
    transport_ids%model(i_model)%ggd(i_slice)%total_ion_energy%flux_parallel(1)%grid_index = grid_ind
    transport_ids%model(i_model)%ggd(i_slice)%total_ion_energy%flux_parallel(1)%grid_subset_index = grid_sub_ind

    allocate(transport_ids%model(i_model)%ggd(i_slice)%ion(1))
    allocate(transport_ids%model(i_model)%ggd(i_slice)%ion(1)%element(1))
    transport_ids%model(i_model)%ggd(i_slice)%ion(1)%element(1)%a = central_mass
    transport_ids%model(i_model)%ggd(i_slice)%ion(1)%element(1)%z_n = 1
    allocate(transport_ids%model(i_model)%ggd(i_slice)%ion(1)%particles%flux_parallel(n_grid_sub))
    allocate(transport_ids%model(i_model)%ggd(i_slice)%ion(1)%particles%flux_parallel(1)%values(n_wall_nodes))
    transport_ids%model(i_model)%ggd(i_slice)%ion(1)%particles%flux_parallel(1)%grid_index = grid_ind
    transport_ids%model(i_model)%ggd(i_slice)%ion(1)%particles%flux_parallel(1)%grid_subset_index = grid_sub_ind

    do i_phi=1, n_phi
      do i_pol=1, n_pol
        i = i_pol + (i_phi-1)*n_pol  !< Global index of reference node
        qpar_e   = res_bnd(i_phi,i_pol,i_qpar_e)
        qpar_i   = res_bnd(i_phi,i_pol,i_qpar_i)
        partf = res_bnd(i_phi,i_pol,i_partf)

        transport_ids%model(i_model)%ggd(i_slice)%electrons%energy%flux_parallel(1)%values(i) = qpar_e
        transport_ids%model(i_model)%ggd(i_slice)%total_ion_energy%flux_parallel(1)%values(i) = qpar_i
        transport_ids%model(i_model)%ggd(i_slice)%ion(1)%particles%flux_parallel(1)%values(i) = partf
      enddo
    enddo

    ! *******************************************************************************
  
  end subroutine fill_transport_IDS
 




  ! --- Fill a pf_passive IDS from STARWALL data, needs to be adapted to CARIDDI!
  subroutine fill_pf_passive_IDS(first_step, time_SI, pf_passive, passive_coil_geo_file)  

    use phys_module, only : F0, central_density, sqrt_mu0_rho0, &
                           sqrt_mu0_over_rho0, central_mass, imp_type, &
                           gamma, index_main_imp
    use vacuum
    use vacuum_response,  only: reconstruct_coil_potentials

    implicit none

    ! --- External parameters
    logical,                 intent(in) :: first_step   ! is this the first step?
    real*8,                  intent(in) :: time_SI
    type(ids_pf_passive), intent(inout) :: pf_passive
    character(len=*),        intent(in) :: passive_coil_geo_file
   
    ! --- Local parameters 
    integer :: i, j, k, m, my_id=0, ierr, i_coil
    integer :: n_slice=1, i_slice=1
    logical :: found_coil
    real*8, allocatable :: pot_c(:) 
    type(t_coil_set_starwall) :: coil_set
   
    pf_passive%ids_properties%homogeneous_time = IDS_TIME_MODE_HETEROGENEOUS

    if (sr%n_diag_coils == 0) then
      write(*,*) '  No diagnostic coils in the STARWALL response file'
      write(*,*) '  needed for pf_passive IDS'
      stop
    endif

    ! --- Read STARWALL coil geometry input file for passive conductors
    call read_coil_set_starwall(passive_coil_geo_file, coil_set)

    ! --- Set times
    allocate( pf_passive%time(n_slice) )
    pf_passive%time(i_slice) = time_SI 

    if (first_step) then
      allocate(pf_passive%loop(coil_set%ncoil))
    endif

    ! --- Reconstruct currents from wall_curr
    call reconstruct_coil_potentials(pot_c, wall_curr, my_id)

    do i_coil=1, coil_set%ncoil
       ! --- Check if the coil given in the input file exists in the JOREK restart
      found_coil = .false.
      do i=1, sr%ncoil
        if (INDEX(trim(sr%coil_name(i)),trim(coil_set%coil(i_coil)%name)) > 0) then
          found_coil = .true.
          exit
        endif
      enddo
      if ( .not. found_coil) then
        write(*,*) coil_set%coil(i_coil)%name, " is not part of the starwall_response.dat file data!"
        stop
      endif

      if (coil_set%coil(i_coil)%coil_type /= AXISYM_THICK) then
        write(*,*) coil_set%coil(i_coil)%coil_type," coil type in pf_passive not supported yet in jorek2_IDS"
        stop
      endif

      ! --- Fill coil geometry for first time step
      if (first_step) then
        allocate( pf_passive%loop(i_coil)%name(1) )
        pf_passive%loop(i_coil)%name       = trim(coil_set%coil(i_coil)%name)   
        pf_passive%loop(i_coil)%resistance = coil_set%coil(i_coil)%resist * wall_resistivity_fact
        allocate( pf_passive%loop(i_coil)%element(coil_set%coil(i_coil)%nparts_coil) )
        do k=1, coil_set%coil(i_coil)%nparts_coil
          pf_passive%loop(i_coil)%element(k)%turns_with_sign           = coil_set%coil(i_coil)%n_thick_turns(k)
          pf_passive%loop(i_coil)%element(k)%geometry%geometry_type    = 2 ! Rectangle type
          pf_passive%loop(i_coil)%element(k)%geometry%rectangle%r      = coil_set%coil(i_coil)%R(k)
          pf_passive%loop(i_coil)%element(k)%geometry%rectangle%z      = coil_set%coil(i_coil)%Z(k)
          pf_passive%loop(i_coil)%element(k)%geometry%rectangle%width  = coil_set%coil(i_coil)%dR(k)
          pf_passive%loop(i_coil)%element(k)%geometry%rectangle%height = coil_set%coil(i_coil)%dZ(k)
        enddo
      endif

      allocate( pf_passive%loop(i_coil)%current(n_slice) )
      allocate( pf_passive%loop(i_coil)%time(n_slice) )

      pf_passive%loop(i_coil)%current(i_slice) = pot_c(i) / MU_ZERO * fact_Ip
      pf_passive%loop(i_coil)%time(i_slice)    = time_SI
    enddo

  end subroutine fill_pf_passive_IDS





    ! --- Fill an SPI IDS
  subroutine fill_spi_IDS(first_step, time_SI, spi)  

    use phys_module, only : F0, central_density, sqrt_mu0_rho0, &
                           sqrt_mu0_over_rho0, central_mass, imp_type, &
                           gamma, index_main_imp, pellets, spi_vel_Rref, spi_Vel_Zref, &
                           spi_Vel_RxZref

    implicit none

    ! --- External parameters
    logical,                 intent(in) :: first_step   ! is this the first step?
    real*8,                  intent(in) :: time_SI
    type(ids_spi),        intent(inout) :: spi
   
    ! --- Local parameters 
    integer :: i_inj, i_frag, i_frag_glob, j, k, m, ierr, n_spi_begin
    integer :: n_slice=1, i_slice=1, a_imp, z_imp
    real*8  :: pellet_vol_imp, pellet_vol_bg, pellet_vol, rho0 
   
    if (.not. using_spi) then
      write(*,*) 'SPI IDS could not be exported! using_spi must be true'
      return
    endif

    spi%ids_properties%homogeneous_time = IDS_TIME_MODE_HOMOGENEOUS

    ! --- Set times
    allocate( spi%time(n_slice) )
    spi%time(i_slice) = time_SI 
    rho0              = central_density * 1.d20 * central_mass * ATOMIC_MASS_UNIT
    sqrt_mu0_rho0     = sqrt( mu_zero * rho0 )

    ! --- Loop over injectors
    allocate( spi%injector(n_inj) )
    n_spi_begin = 1
    do i_inj = 1, n_inj 

      ! --- Pellet properties
      ! --- Total atoms in the unshattered pellet
      spi%injector(i_inj)%pellet%core%atoms_n = spi_quantity(i_inj) + spi_quantity_bg(i_inj) 

      ! --- Info about species
      pellet_vol_imp = spi_quantity(i_inj)    / (pellet_density*1.d20)     ! total volume occupied by impurities in the pellet
      pellet_vol_bg  = spi_quantity_bg(i_inj) / (pellet_density_bg*1.d20)  ! total volume occupied by backg. species in the pellet
      pellet_vol     = pellet_vol_imp + pellet_vol_bg
      
      allocate( spi%injector(i_inj)%pellet%core%species(2) )  ! 1 impurity and 1 background species

      ! --- Impurity species
      call get_element_atomic_numbers(imp_type(index_main_imp), a_imp, z_imp )

      spi%injector(i_inj)%pellet%core%species(1)%a       = a_imp
      spi%injector(i_inj)%pellet%core%species(1)%z_n     = z_imp
      spi%injector(i_inj)%pellet%core%species(1)%density = spi_quantity(i_inj) / pellet_vol

      ! --- Background species
      spi%injector(i_inj)%pellet%core%species(2)%a       = nint(central_mass)
      spi%injector(i_inj)%pellet%core%species(2)%z_n     = 1
      spi%injector(i_inj)%pellet%core%species(2)%density = spi_quantity_bg(i_inj) / pellet_vol

      ! --- Velocity of centre of mass at shattering location
      spi%injector(i_inj)%velocity_mass_centre_fragments_r   = spi_vel_Rref(i_inj)
      spi%injector(i_inj)%velocity_mass_centre_fragments_z   = spi_vel_Zref(i_inj)
      spi%injector(i_inj)%velocity_mass_centre_fragments_phi = spi_vel_RxZref(i_inj) * fact_phi_dir

      ! --- Fragment properties
      allocate( spi%injector(i_inj)%fragment( n_spi(i_inj) ) )
      
      do i_frag=1, n_spi(i_inj)
        allocate( spi%injector(i_inj)%fragment(i_frag)%position%r(n_slice)   )
        allocate( spi%injector(i_inj)%fragment(i_frag)%position%z(n_slice)   )
        allocate( spi%injector(i_inj)%fragment(i_frag)%position%phi(n_slice) )
        allocate( spi%injector(i_inj)%fragment(i_frag)%velocity_r(n_slice)   )
        allocate( spi%injector(i_inj)%fragment(i_frag)%velocity_z(n_slice)   )
        allocate( spi%injector(i_inj)%fragment(i_frag)%velocity_phi(n_slice) )
        allocate( spi%injector(i_inj)%fragment(i_frag)%volume(n_slice)       ) 
        
        i_frag_glob = i_frag - 1 + n_spi_begin

        spi%injector(i_inj)%fragment(i_frag)%position%r(i_slice)   = pellets(i_frag_glob)%spi_r
        spi%injector(i_inj)%fragment(i_frag)%position%z(i_slice)   = pellets(i_frag_glob)%spi_z
        spi%injector(i_inj)%fragment(i_frag)%position%phi(i_slice) = pellets(i_frag_glob)%spi_phi * fact_phi_dir

        spi%injector(i_inj)%fragment(i_frag)%velocity_r(i_slice)   = pellets(i_frag_glob)%spi_vel_r
        spi%injector(i_inj)%fragment(i_frag)%velocity_z(i_slice)   = pellets(i_frag_glob)%spi_vel_z
        spi%injector(i_inj)%fragment(i_frag)%velocity_phi(i_slice) = pellets(i_frag_glob)%spi_vel_rxz * fact_phi_dir

        spi%injector(i_inj)%fragment(i_frag)%volume(i_slice) = 4.d0/3.d0*PI*pellets(i_frag_glob)%spi_radius**3.d0 
      enddo ! --- fragments

      n_spi_begin = n_spi_begin + n_spi(i_inj)
    end do  ! --- injectors

  end subroutine fill_spi_IDS




   ! --- Fill a pf_active IDS from STARWALL data, needs to be adapted to CARIDDI!
  subroutine fill_pf_active_IDS(first_step, time_SI, pf_active, active_coil_geo_file)  

    use phys_module, only : F0, central_density, sqrt_mu0_rho0, &
                           sqrt_mu0_over_rho0, central_mass, imp_type, &
                           gamma, index_main_imp
    use vacuum
    use vacuum_response,  only: reconstruct_coil_potentials

    implicit none

    ! --- External parameters
    logical,                 intent(in) :: first_step   ! is this the first step?
    real*8,                  intent(in) :: time_SI
    type(ids_pf_active),  intent(inout) :: pf_active
    character(len=*),        intent(in) :: active_coil_geo_file
   
    ! --- Local parameters 
    integer :: i, j, k, m, my_id=0, ierr, i_coil
    integer :: n_slice=1, i_slice=1
    logical :: found_coil
    real*8, allocatable :: pot_c(:) 
    type(t_coil_set_starwall) :: coil_set
   
    pf_active%ids_properties%homogeneous_time = IDS_TIME_MODE_HETEROGENEOUS

    if (sr%n_pol_coils == 0) then
      write(*,*) '  No poloidal active coils in the STARWALL response file'
      write(*,*) '  needed for pf_active IDS'
      stop
    endif

    ! --- Read STARWALL coil geometry input file for active conductors
    call read_coil_set_starwall(active_coil_geo_file, coil_set)

    ! --- Set times
    allocate( pf_active%time(n_slice) )
    pf_active%time(i_slice) = time_SI 

    if (first_step) then
      allocate(pf_active%coil(coil_set%ncoil))
    endif

    ! --- Reconstruct currents from wall_curr
    call reconstruct_coil_potentials(pot_c, wall_curr, my_id)

    do i_coil=1, coil_set%ncoil
       ! --- Check if the coil given in the input file exists in the JOREK restart
      found_coil = .false.
      do i=1, sr%ncoil
        if (INDEX(trim(sr%coil_name(i)),trim(coil_set%coil(i_coil)%name)) > 0) then
          found_coil = .true.
          exit
        endif
      enddo
      if ( .not. found_coil) then
        write(*,*) coil_set%coil(i_coil)%name, " is not part of the starwall_response.dat file data!"
        stop
      endif

      if (coil_set%coil(i_coil)%coil_type /= AXISYM_THICK) then
        write(*,*) coil_set%coil(i_coil)%coil_type," coil type in pf_active not supported yet in jorek2_IDS"
        stop
      endif

      ! --- Fill coil geometry for first time step
      if (first_step) then
        allocate( pf_active%coil(i_coil)%name(1) )
        pf_active%coil(i_coil)%name       = trim(coil_set%coil(i_coil)%name)   
        pf_active%coil(i_coil)%resistance = coil_set%coil(i_coil)%resist * wall_resistivity_fact
        allocate( pf_active%coil(i_coil)%element(coil_set%coil(i_coil)%nparts_coil) )
        do k=1, coil_set%coil(i_coil)%nparts_coil
          pf_active%coil(i_coil)%element(k)%turns_with_sign           = coil_set%coil(i_coil)%n_thick_turns(k)
          pf_active%coil(i_coil)%element(k)%geometry%geometry_type    = 2 ! Rectangle type
          pf_active%coil(i_coil)%element(k)%geometry%rectangle%r      = coil_set%coil(i_coil)%R(k)
          pf_active%coil(i_coil)%element(k)%geometry%rectangle%z      = coil_set%coil(i_coil)%Z(k)
          pf_active%coil(i_coil)%element(k)%geometry%rectangle%width  = coil_set%coil(i_coil)%dR(k)
          pf_active%coil(i_coil)%element(k)%geometry%rectangle%height = coil_set%coil(i_coil)%dZ(k)
        enddo
      endif

      allocate( pf_active%coil(i_coil)%current%data(n_slice) )
      allocate( pf_active%coil(i_coil)%current%time(n_slice) )

      pf_active%coil(i_coil)%current%data(i_slice) = pot_c(i) / MU_ZERO * fact_Ip
      pf_active%coil(i_coil)%current%time(i_slice) = time_SI

    enddo

  end subroutine fill_pf_active_IDS






  subroutine fill_radiation_IDS(first_step, time_SI, expr_avg_list, avg, n_grid_1d, radiation_ids)  

    use phys_module, only : F0, central_density, sqrt_mu0_rho0, &
                           sqrt_mu0_over_rho0, central_mass, imp_type, &
                           gamma, index_main_imp
    implicit none

    ! --- External parameters
    logical,                 intent(in) :: first_step   ! is this the first step?
    real*8,                  intent(in) :: time_SI
    type(t_expr_list),       intent(in) :: expr_avg_list ! list of expressions for average
    real*8,                  intent(in) :: avg(:,:)
    integer,                 intent(in) :: n_grid_1d ! (Number of points for 1D profile)
    type(ids_radiation),  intent(inout) :: radiation_ids
   
    ! --- Local parameters 
    integer    :: i, j, k, m, i_exp, var_rad, my_id
    real*8     :: rho0, fact_rad

    ! **********************************************************************************
    ! ******************************* IMAS **********************************************
    ! **********************************************************************************
    type(ids_generic_grid_scalar),      pointer :: ggd_scalar
    type(ids_generic_grid_aos3_root),   pointer :: grid
    
    integer:: num_nodes
    
    integer :: n_slice, i_slice, grid_ind, grid_sub_ind, n_grid_sub, n_grid, i_psi
    ! **********************************************************************************
  
    ! --- Number of grids and grid subsets
    n_grid       = 1
    n_grid_sub   = 1
    grid_ind     = 1  ! Index
    grid_sub_ind = 1  ! Index

    if (use_marker) then
      if (first_step) then
        ! --- Put the grid in GGD
        allocate( radiation_ids%grid_ggd(n_grid) )
        grid => radiation_ids%grid_ggd(grid_ind)
        grid%time = time_SI
        call grid2ggd( grid, node_list, element_list, bnd_node_list, bnd_elm_list )
      else
        if ( associated(radiation_ids%grid_ggd)) then
          call ids_deallocate_struct(radiation_ids%grid_ggd(grid_ind), .false.)     
          deallocate(radiation_ids%grid_ggd)
        endif
      endif
    end if
    
    ! --- Normalization factors for IMAS
    rho0               = central_density * 1.d20 * central_mass * ATOMIC_MASS_UNIT
    sqrt_mu0_rho0      = sqrt( mu_zero * rho0 )
    sqrt_mu0_over_rho0 = sqrt( mu_zero / rho0 )

    fact_rad = 1.d0 / ( (gamma-1.d0) * MU_ZERO * sqrt_mu0_rho0 )

    ! --- Set times
    n_slice = 1  
    i_slice = 1
    allocate(  radiation_ids%time(n_slice) )
    radiation_ids%time(i_slice)                = time_SI 

    radiation_ids%ids_properties%homogeneous_time = IDS_TIME_MODE_HETEROGENEOUS
    allocate( radiation_ids%process(1))   ! --- 1 type of radiation
    allocate( radiation_ids%process(1)%identifier%name(1) )
    allocate( radiation_ids%process(1)%identifier%description(1) )

    radiation_ids%process(1)%identifier%name         = "Line radiation"
    radiation_ids%process(1)%identifier%description  = "Total line radiation"
    radiation_ids%process(1)%identifier%index        = 10

    
    if (use_marker) then
      allocate( radiation_ids%process(1)%ggd(n_slice) )
      radiation_ids%process(1)%ggd(i_slice)%time = time_SI

      ! --- Fill radiation data 
      var_rad = 2

      allocate( radiation_ids%process(1)%ggd(i_slice)%ion(1))
      allocate( radiation_ids%process(1)%ggd(i_slice)%ion(1)%emissivity(n_grid_sub))
      allocate( radiation_ids%process(1)%ggd(i_slice)%ion(1)%name(1) )  

      radiation_ids%process(1)%ggd(i_slice)%ion(1)%name = imp_type(index_main_imp) 

      ggd_scalar => radiation_ids%process(1)%ggd(i_slice)%ion(1)%emissivity(grid_sub_ind)
      call fill_Bezier_coefficients( ggd_scalar, aux_node_list, var_rad, grid_ind, grid_sub_ind, fact_rad )
    else
      allocate( radiation_ids%process(1)%profiles_1d(1))
      radiation_ids%process(1)%profiles_1d(1)%time = time_SI

      ! --- Fill expressions in IDSs
      do i_exp=1, expr_avg_list%n_expr
        
        ! --- psi
        if (expr_avg_list%expr(i_exp)%name=='Psi_N') then
          ! --- Psi_N
          allocate( radiation_ids%process(1)%profiles_1d(i_slice)%grid%rho_pol_norm(n_grid_1d) )
          radiation_ids%process(1)%profiles_1d(i_slice)%grid%psi_magnetic_axis = ES%Psi_axis * fact_psi
          radiation_ids%process(1)%profiles_1d(i_slice)%grid%psi_boundary      = ES%Psi_bnd  * fact_psi
          radiation_ids%process(1)%profiles_1d(i_slice)%grid%rho_pol_norm(:)   = sqrt(avg(:,i_exp))
        end if

        if (expr_avg_list%expr(i_exp)%name=='Psi') then
          ! --- Psi
          allocate( radiation_ids%process(1)%profiles_1d(i_slice)%grid%psi(n_grid_1d) )
          radiation_ids%process(1)%profiles_1d(i_slice)%grid%psi(:)   = avg(:,i_exp) * fact_psi
        end if
        
        if (expr_avg_list%expr(i_exp)%name=='radiation') then
          ! --- radiation
          allocate( radiation_ids%process(1)%profiles_1d(1)%emissivity_ion_total(n_grid_1d))
          radiation_ids%process(1)%profiles_1d(1)%emissivity_ion_total(:) = avg(:,i_exp) 
        end if
        
      end do
      
      allocate( radiation_ids%process(1)%profiles_1d(i_slice)%grid%rho_tor_norm(n_grid_1d) )
      i = expr_avg_list%n_expr + 2
      radiation_ids%process(1)%profiles_1d(i_slice)%grid%rho_tor_norm(:) = sqrt( avg(:,i)/avg(n_grid_1d,i) )

    end if
    
  end subroutine fill_radiation_IDS






  subroutine fill_plasma_profiles_IDS(first_step, time_SI, expr_avg_list, avg, plasma_profiles_ids, n_grid)  

    use phys_module, only : F0, central_density, sqrt_mu0_rho0, &
                           sqrt_mu0_over_rho0, central_mass, imp_type, &
                           gamma, index_main_imp
    implicit none

    ! --- External parameters
    logical,      intent(in) :: first_step   ! is this the first step?
    real*8,       intent(in) :: time_SI
    type(t_expr_list) , intent(inout) :: expr_avg_list ! list of expressions for average
    real*8,             intent(in)    :: avg(:,:)
    integer,            intent(in)    :: n_grid       ! Number of flux surfaces to compute average
    type(ids_plasma_profiles),   intent(inout) :: plasma_profiles_ids
   
    ! --- Local parameters 
    integer    :: i, j, k, m, var_rad, i_var, i_tor, index, index_node, my_id, ierr
    integer    :: a_imp, z_imp, i_ion_main, i_ion_imp
    real*8     :: rho0
    
    ! **********************************************************************************
    ! ******************************* IMAS **********************************************
    ! **********************************************************************************
    integer :: n_slice, i_slice, i_exp, i_psi
    ! **********************************************************************************

    ! --- Set times
    n_slice = 1;   i_slice = 1

    allocate( plasma_profiles_ids%profiles_1d(n_slice) )
    allocate( plasma_profiles_ids%time(n_slice) )

    ! --- Normalization factors for IMAS
    rho0               = central_density * 1.d20 * central_mass * ATOMIC_MASS_UNIT
    sqrt_mu0_rho0      = sqrt( mu_zero * rho0 )
    
    plasma_profiles_ids%ids_properties%homogeneous_time = IDS_TIME_MODE_HOMOGENEOUS     
    plasma_profiles_ids%time(i_slice) = time_SI 

    ! --- Some allocations
    allocate( plasma_profiles_ids%profiles_1d(i_slice)%ion(1+n_adas) ) ! First index is for main ions
    allocate( plasma_profiles_ids%profiles_1d(i_slice)%neutral(1+n_adas) )
    i_ion_main = 1;  i_ion_imp = 2;

    ! --- Fill expressions in IDSs
    do i_exp=1, expr_avg_list%n_expr

      ! --- Psi_N
      if (expr_avg_list%expr(i_exp)%name=='Psi_N') then
        allocate( plasma_profiles_ids%profiles_1d(i_slice)%grid%rho_pol_norm(n_grid) )
        plasma_profiles_ids%profiles_1d(i_slice)%grid%psi_magnetic_axis = ES%Psi_axis * fact_psi
        plasma_profiles_ids%profiles_1d(i_slice)%grid%psi_boundary      = ES%Psi_bnd  * fact_psi
        plasma_profiles_ids%profiles_1d(i_slice)%grid%rho_pol_norm(:)   = sqrt(avg(:,i_exp))
      endif

      ! --- Psi
      if (expr_avg_list%expr(i_exp)%name=='Psi') then
        allocate( plasma_profiles_ids%profiles_1d(i_slice)%grid%psi(n_grid) )
        plasma_profiles_ids%profiles_1d(i_slice)%grid%psi(:)   = avg(:,i_exp) * fact_psi
      endif

      ! --- Ion temperature
      if (expr_avg_list%expr(i_exp)%name=='T_i') then
        allocate( plasma_profiles_ids%profiles_1d(i_slice)%t_i_average(n_grid) )
        plasma_profiles_ids%profiles_1d(i_slice)%t_i_average(:) = avg(:,i_exp)
      endif

      ! --- Electron temperature
      if (expr_avg_list%expr(i_exp)%name=='T_e') then
        allocate( plasma_profiles_ids%profiles_1d(i_slice)%electrons%temperature(n_grid) )
        plasma_profiles_ids%profiles_1d(i_slice)%electrons%temperature(:) = avg(:,i_exp)
      endif

      ! --- Electron density
      if (expr_avg_list%expr(i_exp)%name=='ne') then
        allocate( plasma_profiles_ids%profiles_1d(i_slice)%electrons%density(n_grid) )
        plasma_profiles_ids%profiles_1d(i_slice)%electrons%density(:) = avg(:,i_exp)
      endif

      ! --- Total pressure
      if (expr_avg_list%expr(i_exp)%name=='pres') then
        allocate( plasma_profiles_ids%profiles_1d(i_slice)%pressure_thermal(n_grid) )
        plasma_profiles_ids%profiles_1d(i_slice)%pressure_thermal(:) = avg(:,i_exp)
      endif

      ! --- Electrostatic potential
      if (expr_avg_list%expr(i_exp)%name=='Phi') then
        allocate( plasma_profiles_ids%profiles_1d(i_slice)%phi_potential(n_grid) )
        plasma_profiles_ids%profiles_1d(i_slice)%phi_potential(:) = avg(:,i_exp)
      endif

      ! --- Parallel conductivity
      if (expr_avg_list%expr(i_exp)%name=='eta_T') then
        allocate( plasma_profiles_ids%profiles_1d(i_slice)%conductivity_parallel(n_grid) )
        plasma_profiles_ids%profiles_1d(i_slice)%conductivity_parallel(:) = 1.d0 / avg(:,i_exp)
      endif

      ! --- Parallel current density
      if (expr_avg_list%expr(i_exp)%name=='Jpar') then
        allocate( plasma_profiles_ids%profiles_1d(i_slice)%j_total(n_grid) )
        plasma_profiles_ids%profiles_1d(i_slice)%j_total(:) = avg(:,i_exp)
      endif

      ! --- Parallel electric field
      if (expr_avg_list%expr(i_exp)%name=='E_||') then
        allocate( plasma_profiles_ids%profiles_1d(i_slice)%e_field%parallel(n_grid) )
        plasma_profiles_ids%profiles_1d(i_slice)%e_field%parallel(:) = avg(:,i_exp)
      endif

      ! --- Radial electric field
      if (expr_avg_list%expr(i_exp)%name=='Er') then
        allocate( plasma_profiles_ids%profiles_1d(i_slice)%e_field%radial(n_grid) )
        plasma_profiles_ids%profiles_1d(i_slice)%e_field%radial(:) = avg(:,i_exp)
      endif

      ! --- Parallel velocity
      if (expr_avg_list%expr(i_exp)%name=='vpar') then
        allocate( plasma_profiles_ids%profiles_1d(i_slice)%ion(i_ion_main)%velocity%parallel(n_grid) )
        plasma_profiles_ids%profiles_1d(i_slice)%ion(i_ion_main)%velocity%parallel(:) = avg(:,i_exp)
      endif

      ! --- Poloidal velocity
      if (expr_avg_list%expr(i_exp)%name=='Vtheta_i') then
        allocate( plasma_profiles_ids%profiles_1d(i_slice)%ion(i_ion_main)%velocity%poloidal(n_grid) )
        plasma_profiles_ids%profiles_1d(i_slice)%ion(i_ion_main)%velocity%poloidal(:) = avg(:,i_exp)
      endif

      ! --- Diamagnetic velocity
      if (expr_avg_list%expr(i_exp)%name=='Vstar_i') then
        allocate( plasma_profiles_ids%profiles_1d(i_slice)%ion(i_ion_main)%velocity%diamagnetic(n_grid) )
        plasma_profiles_ids%profiles_1d(i_slice)%ion(i_ion_main)%velocity%diamagnetic(:) = avg(:,i_exp)
      endif

      ! --- Z_eff
      if (expr_avg_list%expr(i_exp)%name=='Z_eff') then
        allocate( plasma_profiles_ids%profiles_1d(i_slice)%zeff(n_grid) )
        plasma_profiles_ids%profiles_1d(i_slice)%zeff(:) = avg(:,i_exp)
      endif

      ! --- Ion density
      if (expr_avg_list%expr(i_exp)%name=='ni_main') then
        allocate( plasma_profiles_ids%profiles_1d(i_slice)%ion(i_ion_main)%density(n_grid) )
        plasma_profiles_ids%profiles_1d(i_slice)%ion(i_ion_main)%density(:) = avg(:,i_exp) 
        allocate( plasma_profiles_ids%profiles_1d(i_slice)%ion(i_ion_main)%element(1) )
        plasma_profiles_ids%profiles_1d(i_slice)%ion(i_ion_main)%element(1)%a   = central_mass
        plasma_profiles_ids%profiles_1d(i_slice)%ion(i_ion_main)%element(1)%z_n = 1   ! Main ions have Z=1
      endif

      ! --- Neutral density (of main ions)
      if (expr_avg_list%expr(i_exp)%name=='nn_main') then
        allocate( plasma_profiles_ids%profiles_1d(i_slice)%neutral(i_ion_main)%density(n_grid) )
        plasma_profiles_ids%profiles_1d(i_slice)%neutral(i_ion_main)%density(:) = avg(:,i_exp) 
        allocate( plasma_profiles_ids%profiles_1d(i_slice)%neutral(i_ion_main)%element(1) )
        plasma_profiles_ids%profiles_1d(i_slice)%neutral(i_ion_main)%element(1)%a   = central_mass
        plasma_profiles_ids%profiles_1d(i_slice)%neutral(i_ion_main)%element(1)%z_n = 1   ! Main ions have Z=1
      endif

      ! --- Main impurity density
      if (expr_avg_list%expr(i_exp)%name=='nimp') then   ! ion index 2 is for the main impurity species
        allocate( plasma_profiles_ids%profiles_1d(i_slice)%ion(i_ion_imp)%density(n_grid) )
        plasma_profiles_ids%profiles_1d(i_slice)%ion(i_ion_imp)%density(:) = avg(:,i_exp)
        allocate( plasma_profiles_ids%profiles_1d(i_slice)%ion(i_ion_imp)%element(1) )
        call get_element_atomic_numbers(imp_type(index_main_imp), a_imp, z_imp )

        plasma_profiles_ids%profiles_1d(i_slice)%ion(i_ion_imp)%element(1)%a   = a_imp
        plasma_profiles_ids%profiles_1d(i_slice)%ion(i_ion_imp)%element(1)%z_n = z_imp
      endif

   end do
   

   index = 1 ! 1 is main impurity
   do i = 1, n_adas
     index = 1 + index
     if (i == index_main_imp) cycle
     if (nimp_bg(i)>1.d-16) then
       allocate( plasma_profiles_ids%profiles_1d(i_slice)%ion(index)%density(n_grid) )
       plasma_profiles_ids%profiles_1d(i_slice)%ion(index)%density(:) = nimp_bg(i) ! nimp_bg is given as 1/m^3
       allocate( plasma_profiles_ids%profiles_1d(i_slice)%ion(index)%element(1) )
       call get_element_atomic_numbers(imp_type(i), a_imp, z_imp )
       plasma_profiles_ids%profiles_1d(i_slice)%ion(index)%element(1)%a   = a_imp
       plasma_profiles_ids%profiles_1d(i_slice)%ion(index)%element(1)%z_n = z_imp
     end if
   end do

    ! --- q-profile
   allocate( plasma_profiles_ids%profiles_1d(i_slice)%q(n_grid) )
   index = expr_avg_list%n_expr + 1
   plasma_profiles_ids%profiles_1d(i_slice)%q(:) = avg(:,index)
    
   ! --- Get rho_norm_tor from psi_N and q_profile
   allocate( plasma_profiles_ids%profiles_1d(i_slice)%grid%rho_tor_norm(n_grid) )
   index = expr_avg_list%n_expr + 2
   plasma_profiles_ids%profiles_1d(i_slice)%grid%rho_tor_norm(:) = sqrt( avg(:,index)/avg(n_grid,index) )

  end subroutine fill_plasma_profiles_IDS






  ! --- This subroutine calls the 0D data subroutine
  ! --- by calling mod_integrals3D only once and avoid duplications
  subroutine get_0D_data(first_step, res0D)

    implicit none

    ! --- External parameters
    logical, intent(in)                :: first_step
    real*8, allocatable, intent(inout) :: res0D(:)  ! result of 0D data
    ! --- Local parameters
    integer             :: ierr
    type(type_command)  :: command_tmp

    step_imported = .true.

    command_tmp%n_args = 0
    if (allocated(res0D)) deallocate(res0D)
    call clean_up()
    call zeroD_quantities(command_tmp, first_step, ierr, res0D)
    call clean_up()
    
  end subroutine get_0D_data


  subroutine get_average_data(first_step, n_grid, expr_avg_list, avg)
    
    implicit none
    
    ! --- Input Variables
    logical, intent(in)                :: first_step
    integer, intent(in)                :: n_grid        ! Number of surfaces
    type(t_expr_list) , intent(inout)  :: expr_avg_list ! list of expressions for average
    real*8, allocatable, intent(inout) :: avg(:,:)  ! array: average + qprof + rho_tor
    ! --- Local Variables
    type(type_command)  :: command_tmp
    real*8, allocatable :: res(:,:)
    real*8, allocatable :: q_prof(:), rho_tor(:)
    integer :: i_psi, ierr

    ! Keep each expression only once: the IDS fill routines allocate one
    ! destination per expression name.
    expr_list = exprs((/'Psi_N', 'Psi', 'pres', 'FFprime_loc', 'p_prime_loc', &
        'Jpar', 'T_i', 'T_e', 'ne', 'Phi', 'eta_T', 'E_||', 'Er', 'vpar', &
        'Vtheta_i', 'Vstar_i', 'rho', 'Z_eff', 'nimp', 'ni_main', 'nn_main', &
        'radiation' /), 22)
    
    command_tmp%n_args = 0
    step_imported = .true.
    
    call clean_up()
    if ( .not. ES%LCFS_is_lost ) then
      if (allocated(res)) deallocate(res)
      call average(command_tmp, first_step, ierr, res, .true.)
      call clean_up()
      call qprofile(command_tmp, first_step, ierr, q_prof)
    else 
      if(allocated(res)) deallocate(res)
      if(allocated(q_prof)) deallocate(q_prof)
      allocate(res(n_grid, expr_list%n_expr))
      allocate(q_prof(n_grid))
      res = -1.d99;   q_prof = -1.d99;
      write(*,*) '  profiles cannot be produced without closed flux surfaces'
    endif
    ! --- Correct first and last points
    q_prof(1)      = q_prof(2)        + (q_prof(2)-q_prof(3))
    q_prof(n_grid) = q_prof(n_grid-1) + (q_prof(n_grid-1)-q_prof(n_grid-2))
    allocate(rho_tor(n_grid))
    rho_tor(:) = 0.d0
    do i_psi=2, n_grid
      rho_tor(i_psi) = sum(q_prof(1:i_psi)) ! Assuming equidistant psi grid!!
    end do

    expr_avg_list = expr_list
    if (allocated(avg)) deallocate(avg)
    allocate(avg(n_grid, expr_list%n_expr+2))
    avg(:,1:expr_list%n_expr) = res(:,:)
    avg(:,expr_list%n_expr + 1) = q_prof(:)
    avg(:,expr_list%n_expr + 2) = rho_tor(:)

    deallocate(res, q_prof, rho_tor)
    
  end subroutine get_average_data






  ! --- Get values at JOREK boundary
  subroutine get_boundary_data(n_phi, result)

    implicit none

    real*8, allocatable, intent(inout)  :: result(:,:,:)
    integer, intent(in)                 :: n_phi   ! Number of points in toroidal direction for boundary data
    character(10)        :: str
    integer              :: ierr, n_phi_points
    type(type_command)   :: command_tmp

    ! --- Arguments for boundary quantities
    step_imported      = .true.
    command_tmp%n_args = 3
    command_tmp%args(1) = '0'              ! phimin
    command_tmp%args(2) = '6.28318530718'  ! phimax
    n_phi_points = max(n_phi, n_tor* 3)
    write(str, '(I0)') n_phi_points
    command_tmp%args(3) = str  
    call clean_up()
    expr_list = exprs((/'x', 'y', 'z', 'Psi', 'A_R', 'A_Z', &
                        'BR', 'BZ', 'Btor', 'JR', 'JZ', 'Jtor',  &
                        'heatF_total', 'partF_total', 'npartF_total', &
                        'qpar_e', 'qpar_i', 'partF_par'/), 18)

    expr_list_bnd = expr_list

    call boundary_quantities(command_tmp, .true., ierr, result)
    call clean_up()
  
  end subroutine get_boundary_data


  


  
  subroutine fill_equilibrium_IDS(first_step, time_SI, n_grid, res0D, expr_avg_list, avg, equilibrium_ids, rect_grid_params)  

    use phys_module, only : F0, central_density, sqrt_mu0_rho0, &
                           sqrt_mu0_over_rho0, central_mass, imp_type, &
                           gamma, index_main_imp
    implicit none

    ! --- External parameters
    logical,             intent(in) :: first_step   ! is this the first step?
    real*8,              intent(in) :: time_SI
    integer,             intent(in) :: n_grid       ! Number of flux surfaces to compute average
    real*8, allocatable, intent(in) :: res0D(:)     ! List of 0D quantities defined in exprs_all_int (mod_expressions.f90)
    type(t_expr_list) , intent(inout) :: expr_avg_list ! list of expressions for average
    real*8,             intent(in)    :: avg(:,:)
    type(t_rect_grid_params), intent(in) :: rect_grid_params ! Parameters that define grid for profiles_2d
    
    type(ids_equilibrium),  intent(inout)  :: equilibrium_ids
   
    ! --- Local parameters 
    integer    :: i, j, k, m, my_id, ierr, ixp1, ixp2, nR, nZ
    real*8     :: rho0, fact_rad, R_min, Z_min, R_max, Z_max, R_node, Z_node
    real*8, allocatable :: R_sep(:), Z_sep(:)
    real*8, allocatable :: result2D(:,:,:), R_vec(:), Z_vec(:)
    character(30)       :: str
    type(type_command)  :: command_tmp
    
    ! **********************************************************************************
    ! ******************************* IMAS **********************************************
    ! **********************************************************************************
    integer :: n_slice, i_slice, i_exp, i_psi
    ! **********************************************************************************

    ! --- Set times
    n_slice = 1;   i_slice = 1

    allocate( equilibrium_ids%time_slice(n_slice) )
    allocate( equilibrium_ids%time(n_slice) )

    ! --- Normalization factors for IMAS
    rho0               = central_density * 1.d20 * central_mass * ATOMIC_MASS_UNIT
    sqrt_mu0_rho0      = sqrt( mu_zero * rho0 )
    
    equilibrium_ids%ids_properties%homogeneous_time = IDS_TIME_MODE_HOMOGENEOUS    
    equilibrium_ids%time(i_slice) = time_SI 

    ! --- Fill profiles
    do i_exp=1, expr_avg_list%n_expr

      ! --- psi
      if (expr_avg_list%expr(i_exp)%name=='Psi') then
        allocate( equilibrium_ids%time_slice(i_slice)%profiles_1d%psi(n_grid) )
        equilibrium_ids%time_slice(i_slice)%profiles_1d%psi(:)   = avg(:,i_exp) * fact_psi
      endif

      ! --- pressure
      if (expr_avg_list%expr(i_exp)%name=='pres') then
        allocate( equilibrium_ids%time_slice(i_slice)%profiles_1d%pressure(n_grid) )
        equilibrium_ids%time_slice(i_slice)%profiles_1d%pressure(:)   = avg(:,i_exp) 
      endif

      ! --- p'
      if (expr_avg_list%expr(i_exp)%name=='p_prime_loc') then
        allocate( equilibrium_ids%time_slice(i_slice)%profiles_1d%dpressure_dpsi(n_grid) )
        equilibrium_ids%time_slice(i_slice)%profiles_1d%dpressure_dpsi(:)   = avg(:,i_exp) / fact_psi
      endif

      ! --- FF'
      if (expr_avg_list%expr(i_exp)%name=='FFprime_loc') then
        allocate( equilibrium_ids%time_slice(i_slice)%profiles_1d%f_df_dpsi(n_grid) )
        equilibrium_ids%time_slice(i_slice)%profiles_1d%f_df_dpsi(:) = avg(:,i_exp) / fact_psi
      endif

      ! --- Parallel current
      if (expr_avg_list%expr(i_exp)%name=='Jpar') then
        allocate( equilibrium_ids%time_slice(i_slice)%profiles_1d%j_parallel(n_grid) )
        equilibrium_ids%time_slice(i_slice)%profiles_1d%j_parallel(:) = avg(:,i_exp) 
      endif

    end do
    
    ! --- q-profile
    allocate( equilibrium_ids%time_slice(i_slice)%profiles_1d%q(n_grid) )
    k = expr_avg_list%n_expr + 1
    equilibrium_ids%time_slice(i_slice)%profiles_1d%q(:) = avg(:,k)

    ! --- Get rho_norm_tor from psi_N and q_profile
    allocate( equilibrium_ids%time_slice(i_slice)%profiles_1d%rho_tor_norm(n_grid) )
    k = expr_avg_list%n_expr + 2
    equilibrium_ids%time_slice(i_slice)%profiles_1d%rho_tor_norm(:) = sqrt( avg(:,k)/avg(n_grid,k) )

    ! --- Information about the toroidal field
    equilibrium_ids%vacuum_toroidal_field%r0 = R_geo
    allocate(equilibrium_ids%vacuum_toroidal_field%b0(n_slice))
    equilibrium_ids%vacuum_toroidal_field%b0(i_slice) = F0/R_geo * fact_Ip
    
    ! --- Fill global quantities (call mod_integrals3D)
    equilibrium_ids%time_slice(i_slice)%global_quantities%psi_magnetic_axis = ES%Psi_axis * fact_psi
    equilibrium_ids%time_slice(i_slice)%global_quantities%magnetic_axis%r = ES%R_axis
    equilibrium_ids%time_slice(i_slice)%global_quantities%magnetic_axis%z = ES%Z_axis
    
    do i_exp=1, exprs_all_int%n_expr

      ! --- Beta poloidal
      if (exprs_all_int%expr(i_exp)%name=='beta_p') then
        equilibrium_ids%time_slice(i_slice)%global_quantities%beta_pol   = res0D(i_exp)
      endif

      ! --- Beta poloidal
      if (exprs_all_int%expr(i_exp)%name=='beta_t') then
        equilibrium_ids%time_slice(i_slice)%global_quantities%beta_tor   = res0D(i_exp)
      endif

      ! --- Normalized beta
      if (exprs_all_int%expr(i_exp)%name=='beta_n') then
        equilibrium_ids%time_slice(i_slice)%global_quantities%beta_tor_norm = abs(res0D(i_exp))
      endif

      ! --- Total current
      if (exprs_all_int%expr(i_exp)%name=='Ip_tot') then
        equilibrium_ids%time_slice(i_slice)%global_quantities%ip = res0D(i_exp) * fact_Ip
      endif

      ! --- li(3)
      if (exprs_all_int%expr(i_exp)%name=='li3') then
        equilibrium_ids%time_slice(i_slice)%global_quantities%li_3 = res0D(i_exp)
      endif

      ! --- Volume
      if (exprs_all_int%expr(i_exp)%name=='volume') then
        equilibrium_ids%time_slice(i_slice)%global_quantities%volume = res0D(i_exp)
      endif

      ! --- Area of poloidal cross section inside LCFS
      if (exprs_all_int%expr(i_exp)%name=='area') then
        equilibrium_ids%time_slice(i_slice)%global_quantities%area = res0D(i_exp)
      endif

      ! --- Current centre - R
      if (exprs_all_int%expr(i_exp)%name=='R_curr_cent') then
        equilibrium_ids%time_slice(i_slice)%global_quantities%current_centre%r = res0D(i_exp)
      endif

      ! --- Current centre - Z
      if (exprs_all_int%expr(i_exp)%name=='Z_curr_cent') then
        equilibrium_ids%time_slice(i_slice)%global_quantities%current_centre%z = res0D(i_exp)
      endif

      ! --- q_axis
      if (exprs_all_int%expr(i_exp)%name=='q02') then
        equilibrium_ids%time_slice(i_slice)%global_quantities%q_axis = res0D(i_exp)
      endif

      ! --- q_95
      if (exprs_all_int%expr(i_exp)%name=='q95') then
        equilibrium_ids%time_slice(i_slice)%global_quantities%q_95 = res0D(i_exp)
      endif

      ! --- Thermal energy
      if (exprs_all_int%expr(i_exp)%name=='Thermal_tot') then
        equilibrium_ids%time_slice(i_slice)%global_quantities%energy_mhd = res0D(i_exp)
      endif

    end do

    ! --- Shaping parameters, T. Luce, PPCF 55 (2013) 095009, equations (1-6)
    if (ES%limiter_plasma) then
      equilibrium_ids%time_slice(i_slice)%boundary%type = 0
    else
      equilibrium_ids%time_slice(i_slice)%boundary%type = 1
    endif
    equilibrium_ids%time_slice(i_slice)%boundary%psi                    = ES%Psi_bnd * fact_psi
    equilibrium_ids%time_slice(i_slice)%boundary%minor_radius           = ES%LCFS_a
    equilibrium_ids%time_slice(i_slice)%boundary%elongation             = ES%LCFS_kappa
    equilibrium_ids%time_slice(i_slice)%boundary%triangularity_upper    = ES%LCFS_deltaU
    equilibrium_ids%time_slice(i_slice)%boundary%triangularity_lower    = ES%LCFS_deltaL
    equilibrium_ids%time_slice(i_slice)%boundary%geometric_axis%r       = ES%LCFS_Rgeo
    equilibrium_ids%time_slice(i_slice)%boundary%geometric_axis%z       = ES%LCFS_Zgeo
    equilibrium_ids%time_slice(i_slice)%boundary%closest_wall_point%r   = ES%R_lim 
    equilibrium_ids%time_slice(i_slice)%boundary%closest_wall_point%z   = ES%Z_lim 
    equilibrium_ids%time_slice(i_slice)%boundary%closest_wall_point%distance = 0.d0

    ! -- Save two special points (axis and X-points)
    allocate(equilibrium_ids%time_slice(i_slice)%contour_tree%node(3))
    !--- Save axis
    if (ES%axis_is_psi_minimum) then
       equilibrium_ids%time_slice(i_slice)%contour_tree%node(1)%critical_type = 0 ! minimum
    else
       equilibrium_ids%time_slice(i_slice)%contour_tree%node(1)%critical_type = 2 ! maximum
    end if
    equilibrium_ids%time_slice(i_slice)%contour_tree%node(1)%r   = ES%R_axis
    equilibrium_ids%time_slice(i_slice)%contour_tree%node(1)%z   = ES%Z_axis
    equilibrium_ids%time_slice(i_slice)%contour_tree%node(1)%psi = ES%psi_axis

    !--- Save X-points
    ixp1 = 1;  ixp2 = 2;
    if (ES%active_xpoint == 2) then
      ixp1 = 2;  ixp2 = 1;
    endif
    equilibrium_ids%time_slice(i_slice)%contour_tree%node(2)%critical_type = 1 ! saddle
    equilibrium_ids%time_slice(i_slice)%contour_tree%node(2)%r   = ES%R_xpoint(ixp1)
    equilibrium_ids%time_slice(i_slice)%contour_tree%node(2)%z   = ES%Z_xpoint(ixp1)
    equilibrium_ids%time_slice(i_slice)%contour_tree%node(2)%psi = ES%Psi_xpoint(ixp1)

    equilibrium_ids%time_slice(i_slice)%contour_tree%node(3)%critical_type = 1 ! saddle
    equilibrium_ids%time_slice(i_slice)%contour_tree%node(3)%r   = ES%R_xpoint(ixp2)
    equilibrium_ids%time_slice(i_slice)%contour_tree%node(3)%z   = ES%Z_xpoint(ixp2)
    equilibrium_ids%time_slice(i_slice)%contour_tree%node(3)%psi = ES%Psi_xpoint(ixp2)

    ! --- Export separatrix
    command_tmp%n_args = 0
    call separatrix(command_tmp, ierr, R_sep, Z_sep)
    allocate(equilibrium_ids%time_slice(i_slice)%boundary%outline%r(size(R_sep)))
    allocate(equilibrium_ids%time_slice(i_slice)%boundary%outline%z(size(R_sep)))
    equilibrium_ids%time_slice(i_slice)%boundary%outline%r(:) = R_sep(:)
    equilibrium_ids%time_slice(i_slice)%boundary%outline%z(:) = Z_sep(:)

    ! --- Export 2D quantities on a rectangular grid
    call clean_up()

    nR    = rect_grid_params%nR
    nZ    = rect_grid_params%nZ
    R_min = rect_grid_params%R_min
    R_max = rect_grid_params%R_max
    Z_min = rect_grid_params%Z_min
    Z_max = rect_grid_params%Z_max

    command_tmp%n_args = 7
    write(str, '(F16.12)') R_min
    command_tmp%args(1) = str  ! Rmin
    write(str, '(F16.12)') R_max
    command_tmp%args(2) = str  ! Rmax
    write(str, '(I0)') nR
    command_tmp%args(3) = str  ! nR
    write(str, '(F16.12)') Z_min
    command_tmp%args(4) = str  ! Zmin
    write(str, '(F16.12)') Z_max
    command_tmp%args(5) = str  ! Zmax
    write(str, '(I0)') nZ
    command_tmp%args(6) = str  ! nZ
    command_tmp%args(7) = '0'  ! phi

    allocate(R_vec(nR), Z_vec(nZ))
    R_vec = [(R_min + float((i-1)) * (R_max-R_min) / float((nR - 1)), i = 1, nR)]
    Z_vec = [(Z_min + float((i-1)) * (Z_max-Z_min) / float((nZ - 1)), i = 1, nZ)]
    
    expr_list = exprs((/'Psi', 'Jtor', 'BR', 'BZ', 'Btor'/), 5)
    call rectangle(command_tmp, first_step, ierr, only_n0=.true., res2D_out=result2D)
    
    allocate(equilibrium_ids%time_slice(i_slice)%profiles_2d(1))
    equilibrium_ids%time_slice(i_slice)%profiles_2d(1)%type%index      = 0
    equilibrium_ids%time_slice(i_slice)%profiles_2d(1)%grid_type%index = 1 ! --- Rectangular
    
    allocate(equilibrium_ids%time_slice(i_slice)%profiles_2d(1)%grid%dim1(nR))
    allocate(equilibrium_ids%time_slice(i_slice)%profiles_2d(1)%grid%dim2(nZ))
    
    equilibrium_ids%time_slice(i_slice)%profiles_2d(1)%grid%dim1(:) = R_vec(:)
    equilibrium_ids%time_slice(i_slice)%profiles_2d(1)%grid%dim2(:) = Z_vec(:)
    
    allocate( equilibrium_ids%time_slice(i_slice)%profiles_2d(1)%r(nR, nZ) )
    allocate( equilibrium_ids%time_slice(i_slice)%profiles_2d(1)%z(nR, nZ) )
    do i=1, nR
      do j=1, nZ
        equilibrium_ids%time_slice(i_slice)%profiles_2d(1)%r(i,j) = R_vec(i)  ! --- Som plotting tools use this field as well
        equilibrium_ids%time_slice(i_slice)%profiles_2d(1)%z(i,j) = Z_vec(j)
      enddo
    enddo

    ! --- Fill profiles
    do i_exp=1, expr_list%n_expr

      ! --- psi
      if (expr_list%expr(i_exp)%name=='Psi') then
        allocate( equilibrium_ids%time_slice(i_slice)%profiles_2d(1)%psi(nR, nZ) )
        equilibrium_ids%time_slice(i_slice)%profiles_2d(1)%psi(:,:) = result2D(:,:,i_exp) * fact_psi
      endif

      ! --- Jtor
      if (expr_list%expr(i_exp)%name=='Jtor') then
        allocate( equilibrium_ids%time_slice(i_slice)%profiles_2d(1)%j_phi(nR, nZ) )
        equilibrium_ids%time_slice(i_slice)%profiles_2d(1)%j_phi(:,:) = result2D(:,:,i_exp) * fact_Ip
      endif

      ! --- B_R
      if (expr_list%expr(i_exp)%name=='BR') then
        allocate( equilibrium_ids%time_slice(i_slice)%profiles_2d(1)%b_field_r(nR, nZ) )
        equilibrium_ids%time_slice(i_slice)%profiles_2d(1)%b_field_r(:,:) = result2D(:,:,i_exp)
      endif

      ! --- B_Z
      if (expr_list%expr(i_exp)%name=='BZ') then
        allocate( equilibrium_ids%time_slice(i_slice)%profiles_2d(1)%b_field_z(nR, nZ) )
        equilibrium_ids%time_slice(i_slice)%profiles_2d(1)%b_field_z(:,:) = result2D(:,:,i_exp)
      endif

      ! --- B_tor
      if (expr_list%expr(i_exp)%name=='Btor') then
        allocate( equilibrium_ids%time_slice(i_slice)%profiles_2d(1)%b_field_phi(nR, nZ) )
        equilibrium_ids%time_slice(i_slice)%profiles_2d(1)%b_field_phi(:,:) = result2D(:,:,i_exp) * fact_Ip
      endif

    enddo

  end subroutine fill_equilibrium_IDS






  subroutine fill_summary_IDS(first_step, time_SI, res0D, summary_ids, simulation_description, &
                              summary_disruption_type, summary_disruption_vertical_displacement)

    implicit none

    ! --- External parameters
    logical,             intent(in) :: first_step   ! is this the first step?
    real*8,              intent(in) :: time_SI
    real*8, allocatable, intent(in) :: res0D(:)     ! List of 0D quantities defined in exprs_all_int (mod_expressions.f90)
    type(ids_summary),  intent(inout)  :: summary_ids
    character(len=1000)             :: simulation_description
    integer, intent(in)             :: summary_disruption_type
    integer, intent(in)             :: summary_disruption_vertical_displacement
   
    ! --- Local parameters 
    integer    :: i, j, k, m, var_rad, i_var, i_tor, index, index_node, my_id, ierr
    
    ! **********************************************************************************
    ! ******************************* IMAS **********************************************
    ! **********************************************************************************
    integer :: n_slice, i_slice, i_exp, i_psi
    ! **********************************************************************************

    ! --- Set times
    n_slice = 1;   i_slice = 1
    allocate( summary_ids%time(n_slice) )
    summary_ids%ids_properties%homogeneous_time = IDS_TIME_MODE_HOMOGENEOUS   
    if (first_step) then 
      allocate(summary_ids%ids_properties%comment(1))
      summary_ids%ids_properties%comment = simulation_description
    endif
    summary_ids%time(i_slice) = time_SI 

    ! --- Information about the toroidal field
    summary_ids%global_quantities%r0%value = R_geo
    allocate(summary_ids%global_quantities%b0%value(n_slice))
    summary_ids%global_quantities%b0%value(i_slice) = F0/R_geo * fact_Ip
    
    do i_exp=1, exprs_all_int%n_expr

      ! --- Beta poloidal
      if (exprs_all_int%expr(i_exp)%name=='beta_p') then
        allocate(summary_ids%global_quantities%beta_pol%value(n_slice))
        summary_ids%global_quantities%beta_pol%value(i_slice) = res0D(i_exp)
      endif

      ! --- Beta toroidal
      if (exprs_all_int%expr(i_exp)%name=='beta_t') then
        allocate(summary_ids%global_quantities%beta_tor%value(n_slice))
        summary_ids%global_quantities%beta_tor%value(i_slice) = res0D(i_exp)
      endif

      ! --- Normalized beta
      if (exprs_all_int%expr(i_exp)%name=='beta_n') then
        allocate(summary_ids%global_quantities%beta_tor_norm%value(n_slice))
        summary_ids%global_quantities%beta_tor_norm%value(i_slice) = res0D(i_exp)
      endif

      ! --- Total current
      if (exprs_all_int%expr(i_exp)%name=='Ip_tot') then
        allocate(summary_ids%global_quantities%ip%value(n_slice))
        summary_ids%global_quantities%ip%value(i_slice) = res0D(i_exp) * fact_Ip
      endif

      ! --- li(3)
      if (exprs_all_int%expr(i_exp)%name=='li3') then
        allocate(summary_ids%global_quantities%li_3%value(n_slice))
        summary_ids%global_quantities%li_3%value(i_slice) = res0D(i_exp)
      endif

      ! --- Volume
      if (exprs_all_int%expr(i_exp)%name=='volume') then
        allocate(summary_ids%global_quantities%volume%value(n_slice))
        summary_ids%global_quantities%volume%value(i_slice) = res0D(i_exp)
      endif

      ! --- q_95
      if (exprs_all_int%expr(i_exp)%name=='q95') then
        allocate(summary_ids%global_quantities%q_95%value(n_slice))
        summary_ids%global_quantities%q_95%value(i_slice) = res0D(i_exp)
      endif

      ! --- Thermal energy
      if (exprs_all_int%expr(i_exp)%name=='Thermal_in') then
        allocate(summary_ids%global_quantities%energy_thermal%value(n_slice))
        summary_ids%global_quantities%energy_thermal%value(i_slice) = res0D(i_exp)
      endif

      ! --- Magnetic energy
      if (exprs_all_int%expr(i_exp)%name=='Wmag_in') then
        allocate(summary_ids%global_quantities%energy_b_field_pol%value(n_slice))
        summary_ids%global_quantities%energy_b_field_pol%value(i_slice) = res0D(i_exp)
      endif

      ! --- Ohmic power
      if (exprs_all_int%expr(i_exp)%name=='Ohmic_tot') then
        allocate(summary_ids%global_quantities%power_ohm%value(n_slice))
        summary_ids%global_quantities%power_ohm%value(i_slice) = res0D(i_exp)
      endif

      ! --- Radiated power
      if (exprs_all_int%expr(i_exp)%name=='Rad_tot') then
        allocate(summary_ids%global_quantities%power_radiated%value(n_slice))
        summary_ids%global_quantities%power_radiated%value(i_slice) = res0D(i_exp)
      endif

      ! --- Heating power
      if (exprs_all_int%expr(i_exp)%name=='Heat_src_tot') then
        allocate(summary_ids%heating_current_drive%power_additional%value(n_slice))
        summary_ids%heating_current_drive%power_additional%value(i_slice) = res0D(i_exp)
      endif

      ! --- R_axis
      if (exprs_all_int%expr(i_exp)%name=='R_axis') then
        allocate(summary_ids%local%magnetic_axis%position%r(n_slice))
        summary_ids%local%magnetic_axis%position%r(i_slice) = res0D(i_exp)
      endif

      ! --- Z_axis
      if (exprs_all_int%expr(i_exp)%name=='Z_axis') then
        allocate(summary_ids%local%magnetic_axis%position%z(n_slice))
        summary_ids%local%magnetic_axis%position%z(i_slice) = res0D(i_exp)
      endif

      ! --- psi_axis
      if (exprs_all_int%expr(i_exp)%name=='psi_axis') then
        allocate(summary_ids%local%magnetic_axis%position%psi(n_slice))
        summary_ids%local%magnetic_axis%position%psi(i_slice) = res0D(i_exp) * fact_psi
      endif

    end do

    ! --- Shaping parameters, T. Luce, PPCF 55 (2013) 095009, equations (1-6)
    allocate(summary_ids%boundary%type%value(n_slice))
    if (ES%limiter_plasma) then
      summary_ids%boundary%type%value(i_slice)  = 0
    else
      summary_ids%boundary%type%value(i_slice)  = 1
    endif
    allocate(summary_ids%boundary%minor_radius%value(n_slice))
    summary_ids%boundary%minor_radius%value(i_slice)           = ES%LCFS_a
    allocate(summary_ids%boundary%elongation%value(n_slice))
    summary_ids%boundary%elongation%value(i_slice)             = ES%LCFS_kappa
    allocate(summary_ids%boundary%triangularity_upper%value(n_slice))
    summary_ids%boundary%triangularity_upper%value(i_slice)    = ES%LCFS_deltaU
    allocate(summary_ids%boundary%triangularity_lower%value(n_slice))
    summary_ids%boundary%triangularity_lower%value(i_slice)    = ES%LCFS_deltaL
    allocate(summary_ids%boundary%geometric_axis_r%value(n_slice))
    summary_ids%boundary%geometric_axis_r%value(i_slice)       = ES%LCFS_Rgeo
    allocate(summary_ids%boundary%geometric_axis_z%value(n_slice))
    summary_ids%boundary%geometric_axis_z%value(i_slice)       = ES%LCFS_Zgeo

    ! --- Optional disruption metadata. Leave these nodes empty unless the
    !     corresponding values were explicitly provided in imas.nml.
    if (summary_disruption_type /= -999) then
      summary_ids%disruption%type%index = summary_disruption_type
    endif
    if (summary_disruption_vertical_displacement /= -999) then
      summary_ids%disruption%vertical_displacement%value = summary_disruption_vertical_displacement
    endif

  end subroutine fill_summary_IDS






  subroutine fill_disruption_IDS(first_step, time_SI, res0D, disruption_ids)  

    implicit none

    ! --- External parameters
    logical,             intent(in) :: first_step   ! is this the first step?
    real*8,              intent(in) :: time_SI
    real*8, allocatable, intent(in) :: res0D(:)     ! List of 0D quantities defined in exprs_all_int (mod_expressions.f90)
    type(ids_disruption),  intent(inout)  :: disruption_ids
   
    ! --- Local parameters 
    integer    :: i, j, k, m, var_rad, i_var, i_tor, index, index_node, my_id, ierr
    real*8     :: vpar_power, kinpar_power
    
    ! **********************************************************************************
    ! ******************************* IMAS **********************************************
    ! **********************************************************************************
    integer :: n_slice, i_slice, i_exp, i_psi
    ! **********************************************************************************

    ! --- Set times
    n_slice = 1;   i_slice = 1
    allocate( disruption_ids%time(n_slice) )
    disruption_ids%ids_properties%homogeneous_time = IDS_TIME_MODE_HOMOGENEOUS    
    disruption_ids%time(i_slice) = time_SI 

    do i_exp=1, exprs_all_int%n_expr

      ! --- Poloidal halos (defined as positive)
      if (exprs_all_int%expr(i_exp)%name=='I_halo') then
        allocate(disruption_ids%global_quantities%current_halo_pol(n_slice))
        disruption_ids%global_quantities%current_halo_pol(i_slice) = res0D(i_exp) * 1.0d6  ! This diagnostic is in MA
      endif

      ! --- Toroidal halos
      if (exprs_all_int%expr(i_exp)%name=='Ip_out') then
        allocate(disruption_ids%global_quantities%current_halo_phi(n_slice))
        disruption_ids%global_quantities%current_halo_phi(i_slice) = res0D(i_exp) * fact_Ip
      endif

      ! --- Total ohmic power
      if (exprs_all_int%expr(i_exp)%name=='Ohmic_tot') then
        allocate(disruption_ids%global_quantities%power_ohm(n_slice))
        disruption_ids%global_quantities%power_ohm(i_slice) = res0D(i_exp) 
      endif

      ! --- Halo ohmic power
      if (exprs_all_int%expr(i_exp)%name=='Ohmic_out') then
        allocate(disruption_ids%global_quantities%power_ohm_halo(n_slice))
        disruption_ids%global_quantities%power_ohm_halo(i_slice) = res0D(i_exp) 
      endif

      ! --- Total thermal and kinetic power flowing into the wall 
      if (exprs_all_int%expr(i_exp)%name=='qn_par') then
        allocate(disruption_ids%global_quantities%power_parallel_halo(n_slice))
        vpar_power   = 0.d0
        kinpar_power = 0.d0
        ! --- Get other contributions to this power
        do j=1, exprs_all_int%n_expr
          if (exprs_all_int%expr(j)%name=='P_vn') then
            vpar_power = res0D(j) 
          endif
          if (exprs_all_int%expr(j)%name=='kinpar_flux') then
            kinpar_power = res0D(j) 
          endif
        enddo
        disruption_ids%global_quantities%power_parallel_halo(i_slice) = res0D(i_exp) + vpar_power + kinpar_power 
      endif

      ! --- Total radiated power
      if (exprs_all_int%expr(i_exp)%name=='Rad_tot') then
        allocate(disruption_ids%global_quantities%power_radiated_electrons_impurities(n_slice))
        disruption_ids%global_quantities%power_radiated_electrons_impurities(i_slice) = res0D(i_exp) 
      endif

    end do

  end subroutine fill_disruption_IDS






  ! --- Fill an IDS with a rectangular grid containing BR and BZ, extending to the vacuum
  ! --- It needs full free-boundary
  subroutine fill_fields_vacuum_extension(first_step, time_SI, plasma_profiles_ids, rect_grid_params, equilibrium_ids)  

    use mod_vacuum_fields, only: mag_field_including_vacuum

    implicit none

    ! --- External parameters
    logical,   intent(in) :: first_step   !< Is this the first step?
    real*8,    intent(in) :: time_SI      !< Time in SI units
    
    type(ids_plasma_profiles), target,     intent(inout) :: plasma_profiles_ids
    type(t_rect_grid_params),                 intent(in) :: rect_grid_params
    type(ids_equilibrium), optional,      intent(inout)  :: equilibrium_ids
    
    ! --- Local parameters
    type(ids_generic_grid_scalar),            pointer :: ggd_scalar
    type(ids_generic_grid_vector_components), pointer :: ggd_vector
    type(ids_generic_grid_aos3_root),         pointer :: grid
    
    integer :: n_slice, i_slice, grid_ind, grid_sub_ind, n_grid_sub, n_grid, num_nodes
    integer :: iR, iZ, i_pol, i_element, n_element, nR, nZ
    integer, allocatable :: rect_elms_vertices_tmp(:,:), rect_elms_vertices(:,:)
    real*8, allocatable  :: B_tot(:,:,:), RZ(:,:), psi_tot(:,:)
    real*8               :: Rmin, Rmax, Zmin, Zmax

    nR   = rect_grid_params%nR
    nZ   = rect_grid_params%nZ
    Rmin = rect_grid_params%R_min
    Rmax = rect_grid_params%R_max
    Zmin = rect_grid_params%Z_min
    Zmax = rect_grid_params%Z_max

    ! --- Construct rectangular RZ grid in the poloidal plane
    allocate(RZ(nR*nZ,2))
    do iR=1, nR
      do iZ=1, nZ
        i_pol = iR + (iZ-1)*nR
        RZ(i_pol,1) = Rmin + float(iR-1)*(Rmax-Rmin) / float(nR-1)
        RZ(i_pol,2) = Zmin + float(iZ-1)*(Zmax-Zmin) / float(nZ-1)
      enddo
    enddo

    ! --- Connectivity matrix for quadrilateral elements
    allocate(rect_elms_vertices_tmp(nR*nZ,4))
    i_element  = 0
    do iZ=1,nZ-1
      do iR=1,nR-1               
        i_element = i_element + 1
        rect_elms_vertices_tmp(i_element,1) = (iZ-1)*nR + iR
        rect_elms_vertices_tmp(i_element,2) = (iZ-1)*nR + iR + 1
        rect_elms_vertices_tmp(i_element,3) = (iZ  )*nR + iR + 1
        rect_elms_vertices_tmp(i_element,4) = (iZ  )*nR + iR
      enddo
    enddo
    n_element = i_element
    allocate(rect_elms_vertices(n_element,4))
    rect_elms_vertices(:,:) = rect_elms_vertices_tmp(1:n_element,:) 
    deallocate(rect_elms_vertices_tmp)
    
    ! --- Get magnetic field in the given grid for different harmonics
    call mag_field_including_vacuum(RZ, B_tot, psi_tot)

    ! --- Fill IDS
    ! --- Number of grids and grid subsets
    n_grid       = 1
    n_grid_sub   = 1
    grid_ind     = 1  ! Index
    grid_sub_ind = 1  ! Index
  
    if (first_step) then
      ! --- Put the grid in GGD
      allocate( plasma_profiles_ids%grid_ggd(n_grid) )
      grid => plasma_profiles_ids%grid_ggd(grid_ind)
      grid%time = time_SI
      call rect_grid2ggd( grid, rect_elms_vertices, RZ )
    else
      if ( associated(plasma_profiles_ids%grid_ggd)) then
        call ids_deallocate_struct(plasma_profiles_ids%grid_ggd(grid_ind), .false.)  
        deallocate(plasma_profiles_ids%grid_ggd)
      endif
    endif

    ! --- Set times
    n_slice = 1  
    i_slice = 1
    allocate(  plasma_profiles_ids%time(n_slice) )
    allocate(  plasma_profiles_ids%ggd(n_slice ) )

    plasma_profiles_ids%ids_properties%homogeneous_time = IDS_TIME_MODE_HETEROGENEOUS 

    allocate(plasma_profiles_ids%ids_properties%comment(1))
  
    plasma_profiles_ids%time(i_slice)     = time_SI 
    plasma_profiles_ids%ggd(i_slice)%time = time_SI

    plasma_profiles_ids%ids_properties%comment = "The magnetic field exported in a rectangular grid for each &
                                                  Fourier harmonic. It can be used for field line tracing or &
                                                  to produce synthetic magnetic diagnostics (for example). "

    
    allocate( plasma_profiles_ids%ggd(i_slice)%b_field(n_grid_sub))
    ggd_vector => plasma_profiles_ids%ggd(i_slice)%b_field(grid_sub_ind)

    ! --- Fill BR
    call fill_values_vector_with_harmonics( ggd_vector, B_tot(:,:,1), grid_ind, grid_sub_ind, 1.d0, 'r')

    ! --- Fill BZ
    call fill_values_vector_with_harmonics( ggd_vector, B_tot(:,:,2), grid_ind, grid_sub_ind, 1.d0, 'z')

    ! --- Fill in n=0 for equilibrium IDS if provided
    if (present(equilibrium_ids)) then
      do iR=1, nR
        do iZ=1, nZ
          i_pol = iR + (iZ-1)*nR
          equilibrium_ids%time_slice(i_slice)%profiles_2d(1)%psi(iR,iZ)       =  psi_tot(i_pol,1) * fact_psi
          equilibrium_ids%time_slice(i_slice)%profiles_2d(1)%b_field_r(iR,iZ) =  B_tot(i_pol,1,1)
          equilibrium_ids%time_slice(i_slice)%profiles_2d(1)%b_field_z(iR,iZ) =  B_tot(i_pol,1,2)
        enddo
      enddo
    endif

  end subroutine fill_fields_vacuum_extension





  ! --- Fills Bezier coefficients and values in GGD
  subroutine fill_values_coefficients( values, coefficients, node_list, var_index, res_fact )
  
    implicit none
  
    ! --- External parameters
    real*8, pointer,             intent(inout) :: values(:), coefficients(:,:)
    type (type_node_list), intent(in)  :: node_list
    integer,               intent(in)  :: var_index
    real*8,                intent(in)  :: res_fact
  
    ! --- Local parameters
    integer :: num_nodes, inode, inode_glob, itor
    
    num_nodes = node_list%n_nodes

    if ( associated(coefficients) ) then
      deallocate(coefficients )
      allocate(coefficients(num_nodes*n_tor, n_degrees) ) 
    else
      allocate(coefficients(num_nodes*n_tor, n_degrees) )
    endif

    if ( associated(values) ) then
      deallocate( values )
      allocate( values(num_nodes*n_tor) ) 
    else
      allocate( values(num_nodes*n_tor) )
    endif
  
    do inode=1, num_nodes    
      do itor=1, n_tor
        
         inode_glob = inode + (itor-1)*num_nodes
         
        if ( (itor .eq. 1) .or. mod(itor, 2) .eq. 0 ) then ! reverse sign of sine modes to transform COCOS convention
           coefficients(inode_glob,:) =   node_list%node(inode)%values(itor,:,var_index) * res_fact
           values(inode_glob)         =   node_list%node(inode)%values(itor,1,var_index) * res_fact
        else
           coefficients(inode_glob,:) = fact_phi_dir *  node_list%node(inode)%values(itor,:,var_index) * res_fact
           values(inode_glob)         = fact_phi_dir *  node_list%node(inode)%values(itor,1,var_index) * res_fact
        end if

      enddo
    enddo

  end subroutine fill_values_coefficients






  ! --- Fills Bezier coefficients and values in GGD
  subroutine fill_Bezier_coefficients( ggd_scalar, node_list, var_index, grid_ind, grid_sub_ind, res_fact )
  
    implicit none
  
    ! --- External parameters
    type(ids_generic_grid_scalar),  intent(inout) ::  ggd_scalar
    type (type_node_list), intent(in)  :: node_list
    integer,               intent(in)  :: var_index, grid_ind, grid_sub_ind
    real*8,                intent(in)  :: res_fact
  
    ggd_scalar%grid_index        = grid_ind
    ggd_scalar%grid_subset_index = grid_sub_ind
 
    call fill_values_coefficients( ggd_scalar%values, ggd_scalar%coefficients, node_list, var_index, res_fact )
  
  end subroutine fill_Bezier_coefficients
  




  ! --- Fills Bezier coefficients in GGD for vectors
  subroutine fill_Bezier_vector_coefficients( ggd_vector, node_list, var_index, grid_ind, grid_sub_ind, res_fact, component)
  
    implicit none
  
    ! --- External parameters
    type(ids_generic_grid_vector_components),  intent(inout) ::  ggd_vector
    type (type_node_list), intent(in)  :: node_list
    integer,               intent(in)  :: var_index, grid_ind, grid_sub_ind
    real*8,                intent(in)  :: res_fact
    character(len=*),      intent(in)  :: component

    ggd_vector%grid_index        = grid_ind
    ggd_vector%grid_subset_index = grid_sub_ind

    select case (trim(component))
      case('radial')
        call fill_values_coefficients( ggd_vector%radial,      ggd_vector%radial_coefficients,      node_list, var_index, res_fact )
      case('diamagnetic')
        call fill_values_coefficients( ggd_vector%diamagnetic, ggd_vector%diamagnetic_coefficients, node_list, var_index, res_fact )
      case('parallel')
        call fill_values_coefficients( ggd_vector%parallel,    ggd_vector%parallel_coefficients,    node_list, var_index, res_fact )
      case('poloidal')
        call fill_values_coefficients( ggd_vector%poloidal,    ggd_vector%poloidal_coefficients,    node_list, var_index, res_fact )
      case('r')
        call fill_values_coefficients( ggd_vector%r,           ggd_vector%r_coefficients,           node_list, var_index, res_fact )
      case('phi')
        call fill_values_coefficients( ggd_vector%phi,         ggd_vector%phi_coefficients,         node_list, var_index, res_fact )
      case('z')
        call fill_values_coefficients( ggd_vector%z,           ggd_vector%z_coefficients,           node_list, var_index, res_fact )
    end select
  
  end subroutine fill_Bezier_vector_coefficients
  




    ! --- Fills simple values in GGD with harmonics
  subroutine fill_values_harmonics_simple( values, node_values, res_fact )
  
    implicit none
  
    ! --- External parameters
    real*8, pointer,     intent(inout) :: values(:)
    real*8, dimension(:,:), intent(in) :: node_values  !< values at nodes on the poloidal plane (i_pol, i_harmonic)
    real*8,                intent(in)  :: res_fact
  
    ! --- Local parameters
    integer :: num_nodes, inode, inode_glob, itor
    
    num_nodes = size(node_values(:,1),1)

    if ( associated(values) ) then
      deallocate( values )
      allocate( values(num_nodes*n_tor) ) 
    else
      allocate( values(num_nodes*n_tor) )
    endif
  
    do inode=1, num_nodes    
      do itor=1, n_tor
        inode_glob = inode + (itor-1)*num_nodes 
        values(inode_glob)= node_values(inode,itor) * res_fact
      enddo
    enddo

  end subroutine fill_values_harmonics_simple
  
  



  ! --- Fills values in GGD for each toroidal harmonic
  subroutine fill_node_values_with_harmonics( ggd_scalar, node_values, grid_ind, grid_sub_ind, res_fact )
  
    implicit none
  
    ! --- External parameters
    type(ids_generic_grid_scalar),  intent(inout) ::  ggd_scalar
    real*8, dimension(:,:), intent(in)  :: node_values  !< values at nodes on the poloidal plane (i_pol, i_harmonic)
    integer,                intent(in)  :: grid_ind, grid_sub_ind
    real*8,                 intent(in)  :: res_fact
  
    ggd_scalar%grid_index = grid_ind
    ggd_scalar%grid_subset_index = grid_sub_ind

    call fill_values_harmonics_simple( ggd_scalar%values, node_values, res_fact )
  
  end subroutine fill_node_values_with_harmonics
  
  




  ! --- Fills simple values in GGD for vectors with harmonics
  subroutine fill_values_vector_with_harmonics( ggd_vector, node_values, grid_ind, grid_sub_ind, res_fact, component)
  
    implicit none
  
    ! --- External parameters
    type(ids_generic_grid_vector_components),  intent(inout) ::  ggd_vector
    real*8, dimension(:,:), intent(in) :: node_values  !< values at nodes on the poloidal plane (i_pol, i_harmonic)
    integer,               intent(in)  :: grid_ind, grid_sub_ind
    real*8,                intent(in)  :: res_fact
    character(len=*),      intent(in)  :: component

    ggd_vector%grid_index        = grid_ind
    ggd_vector%grid_subset_index = grid_sub_ind

    select case (trim(component))
      case('radial')
        call fill_values_harmonics_simple( ggd_vector%radial,      node_values, res_fact )
      case('diamagnetic')
        call fill_values_harmonics_simple( ggd_vector%diamagnetic, node_values, res_fact )
      case('parallel')
        call fill_values_harmonics_simple( ggd_vector%parallel,    node_values, res_fact )
      case('poloidal')
        call fill_values_harmonics_simple( ggd_vector%poloidal,    node_values, res_fact )
      case('r')
        call fill_values_harmonics_simple( ggd_vector%r,           node_values, res_fact )
      case('phi')
        call fill_values_harmonics_simple( ggd_vector%phi,         node_values, res_fact )
      case('z')
        call fill_values_harmonics_simple( ggd_vector%z,           node_values, res_fact )
    end select
  
  end subroutine fill_values_vector_with_harmonics






  
  !< Fills JOREK grid into GGD
  subroutine grid2ggd( grid, node_list, element_list, bnd_node_list, bnd_elm_list )
  
    implicit none
  
    ! --- External parameters
    type(ids_generic_grid_aos3_root),     pointer   :: grid
    type (type_node_list),    intent(in)            :: node_list
    type (type_element_list), intent(in)            :: element_list
    type (type_bnd_node_list),    intent(in)        :: bnd_node_list
    type (type_bnd_element_list), intent(in)        :: bnd_elm_list
  
  
    ! --- Local parameters
    type(ids_generic_grid_dynamic_space), pointer   ::  space_RZ
    type(ids_generic_grid_dynamic_space), pointer   ::  space_fourier
    type(ids_generic_grid_dynamic_space_dimension),     pointer :: ids_cells
    type(ids_generic_grid_dynamic_grid_subset_element), pointer :: sub_elm
    type(type_bnd_element)                             :: bnd_elm


    
    integer:: idx, shot_number, run_number, num_nodes, num_cells   
    integer :: gs_index, i, j
    integer, allocatable :: vertex_elm_array(:,:)
    real*8,  allocatable :: RZ(:,:)
    
    integer :: itor, idof, n_slice, i_slice, grid_ind, grid_sub_ind, n_grid_sub
  
    ! --- create vertex - elements array
    allocate(  vertex_elm_array(n_vertex_max, element_list%n_elements)  )
    do i=1, element_list%n_elements
      vertex_elm_array(:,i) = element_list%element(i)%vertex
    enddo 
  
    ! Get values of R and Z at the nodes
    allocate( RZ(node_list%n_nodes, 2) )
    do i=1, node_list%n_nodes
      RZ(i,:) = node_list%node(i)%x(1,1,:)  
    enddo
  
    ! Write grid geometry
    allocate(  grid%space(2)                           )
    allocate(  grid%space(1)%objects_per_dimension(3)  )
    allocate(  grid%space(1)%coordinates_type(2)       )
  
    ! Set coordinates type to [R, Z]
    grid%space(1)%coordinates_type(1)%index = 4
    grid%space(1)%coordinates_type(2)%index = 3
  
    allocate(grid%identifier%description(1))
    allocate(grid%identifier%name(1))
    grid%identifier%description(1) = "Mesh coming from the JOREK code: combined 2D finite elements space in the poloidal plane &
                                      with Fourier space for the toroidal angle dependence"
    grid%identifier%name  = "JOREK mesh"
    grid%identifier%index = 0   ! Unspecified
  
    num_nodes = size(RZ,1)
    space_RZ  => grid%space(1)
  
    ! Fill simplified grid nodes (uses geometry instead of geometry_2D and misses Bezier representation) 
    allocate( space_RZ%objects_per_dimension(1)%object(num_nodes) )  ! Allocate to number of nodes
    do i=1, num_nodes 
      allocate( space_RZ%objects_per_dimension(1)%object(i)%geometry(2) ) ! Allocate dimensions per each node
      space_RZ%objects_per_dimension(1)%object(i)%geometry(:) = RZ(i,:)
    enddo
  
    ! Fill dummy variables for 1D elements (edges)
    allocate( space_RZ%objects_per_dimension(2)%object(1) )          ! Allocate just one edge    
    allocate( space_RZ%objects_per_dimension(2)%object(1)%nodes(1))  
    space_RZ%objects_per_dimension(2)%object(1)%nodes(1) = 0
  
    ! Fill JOREK 2D elements (or cells)
    num_cells = element_list%n_elements
    ids_cells => space_RZ%objects_per_dimension(3)
    allocate(    ids_cells%object(num_cells) )
    do i=1, num_cells
      allocate(  ids_cells%object(i)%nodes(n_vertex_max)  )
      ids_cells%object(i)%nodes(:) = vertex_elm_array(:,i) 
    enddo
  
    ! Writing grid_subsets
    allocate(grid%grid_subset(2))  
  
    ! Subset for nodes in the combined space
    gs_index = 1
  
    allocate( grid%grid_subset(gs_index)%identifier%name(1)         )
    allocate( grid%grid_subset(gs_index)%identifier%description(1)  )
    grid%grid_subset(gs_index)%identifier%name(1)        = "nodes"
    grid%grid_subset(gs_index)%identifier%index          = 1 
    grid%grid_subset(gs_index)%identifier%description(1) = "The elements of the grid subset are the 0D nodes &
                                                         of the combined RZ x Fourier space (number of nodes &
                                                         is N_poloidal_nodes x N_fourier). "
    grid%grid_subset(gs_index)%dimension                 = 1    ! 1 is the convention for 0D nodes


    
    gs_index = 2  ! for boundary 
  
    allocate( grid%grid_subset(gs_index)%identifier%name(1)         )
    allocate( grid%grid_subset(gs_index)%identifier%description(1)  )
    grid%grid_subset(gs_index)%identifier%name(1)        = "bnd elements"
    grid%grid_subset(gs_index)%identifier%index          = 44 
    grid%grid_subset(gs_index)%identifier%description(1) = "All elements belonging to the boundary &
                                                            and index to the nodes belonging to it."
    grid%grid_subset(gs_index)%dimension                 = 2    ! 2 is the convention for edge
    allocate(grid%grid_subset(gs_index)%element(bnd_elm_list%n_bnd_elements))
    do i =1, bnd_elm_list%n_bnd_elements
       bnd_elm =  bnd_elm_list%bnd_element(i)
       sub_elm => grid%grid_subset(gs_index)%element(i)
       allocate(sub_elm%object(2))
       do j = 1,2
          sub_elm%object(j)%dimension = 1   ! 0D nodes
          sub_elm%object(j)%space     = 1   ! from poloidal space
          sub_elm%object(j)%index = bnd_elm%vertex(j)
       end do
    end do

    ! Fill toroidal space 
    space_fourier  => grid%space(2)
    allocate(    space_fourier%coordinates_type(1)    )
    allocate(    space_fourier%identifier%description(1)  )
    allocate(    space_fourier%identifier%name(1)  )
    space_fourier%coordinates_type(1)%index = 5          ! The coordinate type is 5, phi angle
    space_fourier%geometry_type%index    = n_period   ! Fourier periodicity
    space_fourier%identifier%name        = "Toroidal Fourier space"             
    space_fourier%identifier%description = "Description of the toroidal Fourier harmonics series"  
    space_fourier%identifier%index       = 3          ! 3= Secondary space extending dimensions           
  
    allocate(  space_fourier%objects_per_dimension(1)               )  ! We have only one dimension of
    allocate(  space_fourier%objects_per_dimension(1)%object(n_tor) )  ! toroidal harmonics
  
    do i=1, n_tor
      allocate(  space_fourier%objects_per_dimension(1)%object(i)%geometry(1) )  ! toroidal harmonics
      space_fourier%objects_per_dimension(1)%object(i)%geometry(1) = i
    enddo
  
    ! Fill in grid Bezier coefficients
    ! Needs generalization for JOREK 3D STELLERATOR EXTENSION!!
    space_RZ%geometry_type%index = 0  ! Standard geometry (non Fourier)
    do i=1, num_nodes
      allocate( space_RZ%objects_per_dimension(1)%object(i)%geometry_2d(2,n_degrees) )
      space_RZ%objects_per_dimension(1)%object(i)%geometry_2d(1,:) = node_list%node(i)%x(1,:,1 )   ! R dofs
      space_RZ%objects_per_dimension(1)%object(i)%geometry_2d(2,:) = node_list%node(i)%x(1,:,2 )   ! Z dofs
    enddo
  
    ! JOREK element sizes
    do i=1, num_cells
      allocate( space_RZ%objects_per_dimension(3)%object(i)%geometry_2d(n_degrees, n_vertex_max) )
      do j=1, n_vertex_max
        space_RZ%objects_per_dimension(3)%object(i)%geometry_2d(:,j) = element_list%element(i)%size(j,:)   
      enddo
    enddo
  
  end subroutine grid2ggd






   !< Fills simple rectangular RZ grid including Fourier space into GGD
  subroutine rect_grid2ggd( grid, vertex_elm_array, RZ )
  
    implicit none
  
    ! --- External parameters
    type(ids_generic_grid_aos3_root),     pointer   :: grid
    real*8,  intent(in) :: RZ(:,:)
    integer, intent(in) :: vertex_elm_array(:,:)
  
    ! --- Local parameters
    type(ids_generic_grid_dynamic_space), pointer   ::  space_RZ
    type(ids_generic_grid_dynamic_space), pointer   ::  space_fourier
    type(ids_generic_grid_dynamic_space_dimension), pointer :: ids_cells
    
    integer:: idx, shot_number, run_number, num_nodes, num_cells   
    integer :: gs_index, i, j
    
    integer :: itor, idof, n_slice, i_slice, grid_ind, grid_sub_ind, n_grid_sub

    num_nodes = size(RZ(:,1),1)
    num_cells = size(vertex_elm_array(:,1),1)

    ! Write grid geometry
    allocate(  grid%space(2)                           )
    allocate(  grid%space(1)%objects_per_dimension(3)  )
    allocate(  grid%space(1)%coordinates_type(2)       )
  
    ! Set coordinates type to [R, Z]
    grid%space(1)%coordinates_type(1)%index = 4
    grid%space(1)%coordinates_type(2)%index = 3
  
    allocate(grid%identifier%description(1))
    allocate(grid%identifier%name(1))
    grid%identifier%description = "Mesh coming from the JOREK code: combined 2D space in the poloidal plane &
                                   with Fourier space for the toroidal angle dependence"
    grid%identifier%name  = "JOREK simple rectangular mesh"
    grid%identifier%index = 0   ! Unspecified
  
    space_RZ  => grid%space(1)
  
    ! Fill simplified grid nodes (uses geometry instead of geometry_2D and misses Bezier representation) 
    allocate( space_RZ%objects_per_dimension(1)%object(num_nodes) )  ! Allocate to number of nodes
    do i=1, num_nodes 
      allocate( space_RZ%objects_per_dimension(1)%object(i)%geometry(2) ) ! Allocate dimensions per each node
      space_RZ%objects_per_dimension(1)%object(i)%geometry(:) = RZ(i,:)
    enddo
  
    ! Fill dummy variables for 1D elements (edges)
    allocate( space_RZ%objects_per_dimension(2)%object(1) )          ! Allocate just one edge    
    allocate( space_RZ%objects_per_dimension(2)%object(1)%nodes(1))  
    space_RZ%objects_per_dimension(2)%object(1)%nodes(1) = 0

    space_RZ%geometry_type%index = 0  ! Standard geometry (non Fourier)
  
    ! Fill JOREK 2D elements (or cells)
    ids_cells => space_RZ%objects_per_dimension(3)
    allocate( ids_cells%object(num_cells) )
    do i=1, num_cells
      allocate(  ids_cells%object(i)%nodes(4)  )
      ids_cells%object(i)%nodes(:) = vertex_elm_array(i,:) 
    enddo
  
    ! Writing grid_subsets
    allocate(grid%grid_subset(1))  
  
    ! Subset for nodes in the combined space
    gs_index = 1
  
    allocate( grid%grid_subset(gs_index)%identifier%name(1)         )
    allocate( grid%grid_subset(gs_index)%identifier%description(1)  )
    grid%grid_subset(gs_index)%identifier%name        = "nodes"
    grid%grid_subset(gs_index)%identifier%index       = 1 
    grid%grid_subset(gs_index)%identifier%description = "The elements of the grid subset are the 0D nodes &
                                                      of the combined RZ x Fourier space (number of nodes &
                                                      is N_poloidal_nodes x N_fourier). "
    grid%grid_subset(gs_index)%dimension              = 1    ! 1 is the convention for 0D nodes
  
    ! Fill toroidal space 
    space_fourier  => grid%space(2)
    allocate(    space_fourier%coordinates_type(1)    )
    allocate(    space_fourier%identifier%description(1)  )
    allocate(    space_fourier%identifier%name(1)  )
    space_fourier%coordinates_type(1)%index = 5          ! The coordinate type is 5, phi angle
    space_fourier%identifier%description    = "Description of the toroidal Fourier harmonics series"        
    space_fourier%identifier%name           = "Toroidal Fourier space"    
    space_fourier%identifier%index          = 3          ! 3= Secondary space extending dimensions
    space_fourier%geometry_type%index       = n_period   ! Fourier periodicity
  
    allocate(  space_fourier%objects_per_dimension(1)               )  ! We have only one dimension of
    allocate(  space_fourier%objects_per_dimension(1)%object(n_tor) )  ! toroidal harmonics
  
    do i=1, n_tor
      allocate(  space_fourier%objects_per_dimension(1)%object(i)%geometry(1) )  ! toroidal harmonics
      space_fourier%objects_per_dimension(1)%object(i)%geometry(1) = i
    enddo


  end subroutine rect_grid2ggd





  !< Fills a triangular thin wall grid into GGD. The wall is described with linear triangles, so only node coordinates and triangle connectivity are needed.
  subroutine triang_thin_wall_2ggd( grid, nodes_xyz, tria_connectivity, subset_choice)
  
    implicit none
  
    ! --- External parameters
    type(ids_generic_grid_aos3_root),  pointer, intent(inout)   :: grid
    real*8,           intent(in) :: nodes_xyz(:,:)           !< First index is node number, second index is coordinate (1,2,3 for x,y,z)
    integer,          intent(in) :: tria_connectivity(:,:)   !< First index is triangle number, second index is global node index of each triangle node (3 nodes per triangle)
    character(len=*), intent(in) :: subset_choice            !< Choice of grid subset to create: "nodes" or "triangles"

    ! --- Local parameters
    integer :: n_wall_nodes, n_wall_triangles, i

    n_wall_nodes    = size(nodes_xyz,1)
    n_wall_triangles = size(tria_connectivity,1)

    grid%identifier%index = 0   ! Unspecified
    allocate( grid%identifier%description(1))
    allocate( grid%identifier%name(1))
    grid%identifier%description = "Thin wall described with linear triangles"
    grid%identifier%name        = "Thin triangular wall"

    allocate(grid%space(1))

    ! --- Identifier
    allocate( grid%space(1)%identifier%description(1))
    allocate( grid%space(1)%identifier%name(1))
    grid%space(1)%identifier%index       = 1  ! Primary space
    grid%space(1)%identifier%description = "This is just a 3D cartesian space"
    grid%space(1)%identifier%name        = "Primary space"

    grid%space(1)%geometry_type%index = 0  ! Standard (not Fourier)
    allocate(grid%space(1)%coordinates_type(3))

    ! --- Identifiers for (x,y,z) type coordinates (1,2,3)
    grid%space(1)%coordinates_type(1)%index = 1
    grid%space(1)%coordinates_type(2)%index = 2
    grid%space(1)%coordinates_type(3)%index = 3

    allocate(grid%space(1)%objects_per_dimension(3))

    ! --- Save wall grid nodes
    allocate(grid%space(1)%objects_per_dimension(1)%object(n_wall_nodes))
    grid%space(1)%objects_per_dimension(1)%geometry_content%index = 1  ! node coordinates
    do i=1, n_wall_nodes
      allocate( grid%space(1)%objects_per_dimension(1)%object(i)%geometry(3) ) ! Allocate dimensions for each node
      grid%space(1)%objects_per_dimension(1)%object(i)%geometry(:) = nodes_xyz(i, :)  
    enddo

    ! --- Save thin wall triangles 
    allocate(grid%space(1)%objects_per_dimension(3)%object(n_wall_triangles))  ! Index 3 for 2D objects (faces)
    do i=1, n_wall_triangles
      allocate( grid%space(1)%objects_per_dimension(3)%object(i)%nodes(3) ) ! 3 nodes per triangle
      grid%space(1)%objects_per_dimension(3)%object(i)%nodes(:) = tria_connectivity(i, :)
    enddo

    ! --- Create the grid substet for either nodes or triangles 
    ! ---   (values are asigned either at triangle centers or at nodes, so only one of the two subsets is needed) 
    allocate(grid%grid_subset(1))
    allocate(grid%grid_subset(1)%identifier%name(1))
    allocate(grid%grid_subset(1)%identifier%description(1))
    
    if (trim(subset_choice) == "nodes") then
      grid%grid_subset(1)%identifier%index = 1  ! Identifier index for 0D nodes in the dictionary
      grid%grid_subset(1)%identifier%name  = "0D nodes"
      grid%grid_subset(1)%identifier%description = "Triangle nodes of the grid"
      grid%grid_subset(1)%dimension        = 1  ! Index 1 means 0 dimensions in the dictionary (for 0D nodes)
    else if (trim(subset_choice) == "triangles") then
      grid%grid_subset(1)%identifier%index = 5  ! Identifier index for 2D cells in the dictionary
      grid%grid_subset(1)%identifier%name = "2D triangles"
      grid%grid_subset(1)%identifier%description = "2D cells representing the linear thin triangles"
      grid%grid_subset(1)%dimension        = 3  ! Index 3 means 2 dimensions in the dictionary
    else
      write(*,*) "Error: subset_choice should be either 'nodes' or 'triangles'"
      stop
    end if
  
  end subroutine triang_thin_wall_2ggd





  ! --- Create a triangular thing wall from the JOREK boundary
  subroutine triangulate_thin_wall_from_bnd(bnd_points, PFC_wall)

    implicit none

    real*8,  intent(in) :: bnd_points(:,:,:)  !< indexing: (i_phi, i_pol, i_quantity), i_quantity=1,2,3 for x,y,z coordinates
    type(t_PFC_triang_grid), intent(inout) :: PFC_wall

    integer :: i_phi, i_pol, ipol1, ipol2, ipol3, ipol4, itor1, itor2, itor3, itor4, n_phi, n_pol
    integer :: i, i1, i2, i3, i4, i_tri, n_wall_triangles, n_wall_nodes

    n_phi          = size(bnd_points, 1)
    n_pol          = size(bnd_points, 2)
    PFC_wall%n_pol = n_pol
    PFC_wall%n_phi = n_phi

    n_wall_triangles      = 2 * n_phi * n_pol
    n_wall_nodes          = n_phi * n_pol
    PFC_wall%n_triangles  = n_wall_triangles
    PFC_wall%n_nodes      = n_wall_nodes

    if (allocated(PFC_wall%tria_xyz)) deallocate(PFC_wall%tria_xyz)
    if (allocated(PFC_wall%tria_connectivity)) deallocate(PFC_wall%tria_connectivity)
    allocate(PFC_wall%tria_xyz(n_wall_nodes,3))
    allocate(PFC_wall%tria_connectivity(n_wall_triangles,3))

    ! --- Construct triangles and connectivity
    i_tri = 0
    do i_phi=1, n_phi
      do i_pol=1, n_pol
        
        i = i_pol + (i_phi-1)*n_pol  !< Global index of refence node

        ! --- Get global indices of the 4 nodes forming a quadrilateral, and make two triangles out of it
        ipol1 = i_pol;                  itor1 = i_phi;
        ipol2 = mod(i_pol,n_pol) + 1;   itor2 = i_phi;
        ipol3 = mod(i_pol,n_pol) + 1;   itor3 = mod(i_phi,n_phi) + 1;
        ipol4 = i_pol;                  itor4 = mod(i_phi,n_phi) + 1;

        i1 = ipol1 + (itor1-1)*n_pol
        i2 = ipol2 + (itor2-1)*n_pol
        i3 = ipol3 + (itor3-1)*n_pol
        i4 = ipol4 + (itor4-1)*n_pol
        
        ! --- Fill in reference node coordinates
        PFC_wall%tria_xyz(i,:) =  (/ bnd_points(itor1,ipol1,1), bnd_points(itor1,ipol1,2), bnd_points(itor1,ipol1,3) /) ! x, y, z coordinates

        ! --- Fill in two triangles
        ! --- Triangle 1
        i_tri = i_tri + 1
        PFC_wall%tria_connectivity(i_tri,:) = (/ i1, i2, i3 /)  ! The node indices of this triangle

        ! --- Triangle 2
        i_tri = i_tri + 1
        PFC_wall%tria_connectivity(i_tri,:) = (/ i1, i3, i4 /)  ! The node indices of this triangle
      enddo
    enddo

  end subroutine triangulate_thin_wall_from_bnd







  !> Read one coils set from a STARWALL coil file
  subroutine read_coil_set_starwall(filename, coil_set)
    
    ! --- Routine parameters
    character(len=*), intent(in)              :: filename
    type(t_coil_set_starwall), intent(inout)  :: coil_set
    
    ! --- Local variables
    character(len=256)        :: description
    integer                   :: i, ncoil, err
    integer, parameter        :: IOCH = 167
    type(t_coil_starwall), allocatable :: coil(:)
    
    ! --- Namelists
    namelist /coil_set_nml/ description, ncoil
    namelist /coils_nml   / coil
    
    write(*,*) '  reading coil file =  ', trim(filename)
      
    open(IOCH, file=trim(filename), status='old', action='read', iostat=err)
    if ( err /= 0 ) then
      write(*,*) '  ERROR opening file', trim(filename), '!'
      stop
    end if
      
    ! --- First namelist
    read(IOCH, coil_set_nml)
    coil_set%description = description
    coil_set%ncoil       = ncoil
      
    ! --- Second namelist
    allocate( coil_set%coil(coil_set%ncoil), coil(coil_set%ncoil) )
    read(IOCH, coils_nml)

    ! --- Calculate turns for AXISYM_FILA coils from n_fila_turns
    do i = 1, coil_set%ncoil
      if ( coil(i)%coil_type == AXISYM_FILA ) then
        coil(i)%nturns = sum( abs( coil(i)%n_fila_turns(1:coil(i)%n_fila) ) )
      end if
      if (  coil(i)%coil_type == AXISYM_THICK) then
        if (minval(coil(i)%n_thick_turns(1:coil(i)%nparts_coil)) > -1.d4) then
          coil(i)%nturns = sum( abs(coil(i)%n_thick_turns(1:coil(i)%nparts_coil)) )
          !write(*,*) 'Warning: As n_thick_turns was specified nturns will automatically be calculated from it.'
        else if (coil(i)%nparts_coil .eq. 1) then
          coil(i)%n_thick_turns(1) = coil(i)%nturns
          !write(*,*) 'Warning: The number of turns is now specified by n_thick_turns'
          !write(*,*) 'To model the old behavior n_thick_turns is set to nturns'
        else
          write(*,*) 'Error: When nparts_coil is bigger than 1'
          write(*,*) 'the value the turns has to be given in n_thick_turns'
          stop
        end if
        !write(*,*) 'Coil ',i,' has ', coil(i)%nturns, ' turns'
      end if
    end do
    coil_set%coil(:) = coil(:)
    deallocate(coil)
    
    close(IOCH)
      
    
  end subroutine read_coil_set_starwall






  ! --- Get atomic element identifiers for impurities
  subroutine get_element_atomic_numbers(impurity_name, a_imp, z_imp)

    implicit none

    character(len=*),    intent(in)    :: impurity_name
    integer,             intent(inout) :: a_imp, z_imp

    select case ( trim(impurity_name) )
    case('D2')
      a_imp  = 2
      z_imp  = 1
    case('Ar')
      a_imp  = 40
      z_imp  = 18
    case('Ne')
      a_imp  = 20
      z_imp  = 10
    case('O')
      a_imp  = 16
      z_imp  = 8
    case('Be')
        a_imp  = 9
        z_imp  = 4
    case('W')
        a_imp  = 184
        z_imp  = 74
    case default
        write(*,*) '!! Impurity type "', trim(imp_type(index_main_imp)), '" unknown !!'
        write(*,*) '=> We assume D2.'
        a_imp  = 2
        z_imp  = 1
    end select

  end subroutine get_element_atomic_numbers






  subroutine initialise_postproc_settings(first_step, n_grid)
    
    implicit none

    logical, intent(in) :: first_step
    integer, intent(in) :: n_grid
    
    character(30)       :: str
    integer             :: ierr

    if (first_step) then 
      call init_new_diag(.false.)
      write(str, '(I0)') n_grid
      call set_setting('units',           '1',     ierr, 'Calculate quantities in which units (0=JOREK, 1=SI)')
      call set_setting('loop_units',      '1',     ierr, 'Use which units for time-loops (0=JOREK, 1=SI)'     )
      call set_setting('linepoints',      '200',   ierr, 'Number of points along a line e.g. for pol_line'    )
      call set_setting('tor_points',      '200',   ierr, 'Number of toroidal points e.g. for tor_line'        )
      call set_setting('surfaces',         str,    ierr, 'number for flux surfaces e.g. for qprofile'         )
      call set_setting('nsmallsteps',     '5',     ierr, 'numerical parameter for field line tracing'         )
      call set_setting('nmaxsteps',       '2500',  ierr, 'numerical parameter for field line tracing'         )
      call set_setting('deltaphi',        '0.1',   ierr, 'numerical parameter for field line tracing'         )
      call set_setting('rad_range_min',   '0.01',  ierr, 'numerical parameter for field line tracing'         )
      call set_setting('rad_range_max',   '0.995', ierr, 'numerical parameter for field line tracing'         )
      call set_setting('nTht',            '32',    ierr, 'numerical parameter for field line tracing'         )
      call set_setting('nsub_bnd',         '4',    ierr, 'numerical parameter for field line tracing'         )
    endif
  
  end subroutine initialise_postproc_settings

#endif

end module mod_jorek2IMAS
