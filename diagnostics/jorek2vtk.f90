!> Program to convert a JOREK2 restart file into binary VTK format
program jorek2vtk

use mod_parameters, only: n_var, variable_names
use data_structure
use phys_module
use basis_at_gaussian
use diffusivities, only: get_dperp, get_zkperp
use pellet_module
use mpi_mod
use mod_bootstrap_functions
use corr_neg
use mod_import_restart
use equil_info
use mod_boundary
use mod_vtk
use mod_interp
use mod_poloidal_currents
use mod_impurity, only: init_imp_adas, radiation_function, radiation_function_linear
use mod_atomic_coeff_deuterium, only : atomic_coeff_deuterium
use mod_openadas , only : read_adf11
use mod_atomic_coeff_deuterium, only : ad_deuterium , atomic_coeff_deuterium

implicit none

type (type_node_list)   ,     pointer :: node_list
type (type_node_list)   ,     pointer :: aux_node_list
type (type_element_list),     pointer :: element_list
type (type_bnd_element_list), pointer :: bnd_elm_list    
type (type_bnd_node_list),    pointer :: bnd_node_list 

integer               :: nnoel, nnos, nel, nsub, inode, ielm, n_scalars, n_vectors
real*4,allocatable    :: currdens(:), xyz (:,:), scalars(:,:), vectors(:,:,:)
integer,allocatable   :: ien (:,:)
integer, parameter    :: ivtk = 22 ! an arbitrary unit number for the VTK output file
integer               :: i, j, k, m, etype, irst, int, i_var, i_tor, i_tor_old, i_plane, index, index_node, my_id
character             :: buffer*80, lf*1, str1*12, str2*12
character*36, allocatable :: scalar_names(:), vector_names(:), variable_names_si(:)
real*8                :: s, t
real*8                :: P,P_s,P_t,P_st,P_ss,P_tt
real*8                :: R,R_s,R_t,R_st,R_ss,R_tt
real*8                :: Z,Z_s,Z_t,Z_st,Z_ss,Z_tt
real*8                :: Ps0, Ps0_s, Ps0_t, Ps0_st, Ps0_ss, Ps0_tt, Psi, Ps_s, Ps_t, Ps_st, Ps_ss, Ps_tt
real*8                :: ZJ0, ZJ0_s, ZJ0_t, ZJ0_st, ZJ0_ss, ZJ0_tt, ZJ,  ZJ_s, ZJ_t, ZJ_st, ZJ_ss, ZJ_tt
real*8                :: U0,  U0_s,  U0_t,  U0_st,  U0_ss,  U0_tt,  U,   U_s,  U_t,  U_st,  U_ss,  U_tt
real*8                :: W0,  W0_s,  W0_t,  W0_st,  W0_ss,  W0_tt,  W,   W_s,  W_t,  W_st,  W_ss,  W_tt
real*8                :: ZN0, ZN0_s, ZN0_t, ZN0_st, ZN0_ss, ZN0_tt, RHO, RHO_s,RHO_t,RHO_st,RHO_ss,RHO_tt
real*8                :: T0,  T0_s,  T0_t,  T0_st,  T0_ss,  T0_tt,  TT,  TT_s, TT_t, TT_st, TT_ss, TT_tt
real*8                :: Ti0, Ti0_s, Ti0_t, Ti0_st, Ti0_ss, Ti0_tt, Ti,  Ti_s, Ti_t, Ti_st, Ti_ss, Ti_tt
real*8                :: Te0, Te0_s, Te0_t, Te0_st, Te0_ss, Te0_tt, Te, Te_s,  Te_t, Te_st, Te_ss, Te_tt
real*8                :: V0,  V0_s,  V0_t,  V0_st,  V0_ss,  V0_tt,  V, V_s, V_t, V_st, V_ss, V_tt, V_sum
real*8                :: dPsi, dPs_s, dPs_t, dPs_st, dPs_ss, dPs_tt
real*8                :: dU,    dU_s,  dU_t,  dU_st,  dU_ss,  dU_tt
real*8                :: ps0_x, ps0_y, psi_sum, ps_x, ps_y, ps_p
real*8                :: u0_x,  u0_y,  u_sum,   u_x,  u_y,  u_p
real*8                :: zj0_x, zj0_y, zj_sum,  zj_x, zj_y, zj_p
real*8                :: w0_x,  w0_y,  w_sum,   w0_xx, w0_yy, w_x, w_y, w_p, w_xx, w_yy
real*8                :: zn0_x, zn0_y, zn_sum,  zn_x, zn_y, zn_p, rho_x, rho_y, rho_p
real*8                :: T0_x,  T0_y,  T_sum,   TT_x, TT_y, TT_p
real*8                :: Ti0_x, Ti0_y, Ti_sum,  Ti_x, Ti_y, Ti_p
real*8                :: Te0_x, Te0_y, Te_sum,  Te_x, Te_y, Te_p
real*8                :: AR0, AR0_s, AR0_t, AR0_st, AR0_ss, AR0_tt
real*8                :: AZ0, AZ0_s, AZ0_t, AZ0_st, AZ0_ss, AZ0_tt
real*8                :: A30, A30_s, A30_t, A30_st, A30_ss, A30_tt
real*8                :: AR,  AR_s,  AR_t,  AR_st,  AR_ss,  AR_tt, AR_R, AR_Z, AR_p, AR_RR, AR_ZZ, AR_RZ, AR_Rp, AR_Zp, AR_pp
real*8                :: AZ,  AZ_s,  AZ_t,  AZ_st,  AZ_ss,  AZ_tt, AZ_R, AZ_Z, AZ_p, AZ_RR, AZ_ZZ, AZ_RZ, AZ_Rp, AZ_Zp, AZ_pp
real*8                :: A3,  A3_s,  A3_t,  A3_st,  A3_ss,  A3_tt, A3_R, A3_Z, A3_p, A3_RR, A3_ZZ, A3_RZ, A3_Rp, A3_Zp, A3_pp
real*8                :: VR0, VR0_s, VR0_t, VR0_st, VR0_ss, VR0_tt
real*8                :: VZ0, VZ0_s, VZ0_t, VZ0_st, VZ0_ss, VZ0_tt
real*8                :: VP0, VP0_s, VP0_t, VP0_st, VP0_ss, VP0_tt
real*8                :: VR,  VR_s,  VR_t,  VR_st,  VR_ss,  VR_tt
real*8                :: VZ,  VZ_s,  VZ_t,  VZ_st,  VZ_ss,  VZ_tt
real*8                :: VP,  VP_s,  VP_t,  VP_st,  VP_ss,  VP_tt
real*8                :: Fprof
real*8                :: BR, BR_R, BR_Z, BR_p
real*8                :: BZ, BZ_R, BZ_Z, BZ_p
real*8                :: BP, BP_R, BP_Z, BP_p
real*8                :: JxB_R,   JxB_Z,   JxB_p,   JxB_pol
real*8                :: GradP_R, GradP_Z, GradP_p, GradP_pol
real*8                :: psi_axis,      R_axis,      Z_axis,      s_axis,      t_axis
real*8                :: psi_xpoint(2), R_xpoint(2), Z_xpoint(2), s_xpoint(2), t_xpoint(2)
real*8                :: psi_norm, psi_bnd, grad_psi
real*8                :: J_phi, J_R, J_Z, eta_T
real*8                :: E_phi, E_R, E_Z, dU_x, dU_y, Jpol_R, Jpol_Z, FFp
real*8                :: xjac, xjac_x, xjac_y, v_perp, Psi_J, R_p, error, Btot, BigR
real*8                :: particle_source, D_prof, ZK_prof, source_pellet, ZKpar_T
real*8                :: Jb,rho_norm,t_norm
integer               :: i_elm_axis, i_elm_xpoint(2), k_tor, ifail, ierr
logical               :: without_n0_mode, SI_units
logical               :: include_fluxes, include_neo, include_magnetic_field, include_velocity_field
logical               :: include_bootstrap, include_psi_norm, include_electric_field, include_Jpol, RphiZ_coords
logical               :: include_projections
character*80          :: proj_basename, filename_proj
real*8                :: toroidal_angle

real*8                :: Er, psi_abs, Vtheta, Btheta, Mach_par,Mach_pol,Vsound, Vneo
real*8                :: amu_neo_node, aki_neo_node
real*8                :: Vperp_e, Psi_tot

real*8                :: angle, source_volume, local_density, local_temperature, local_pressure, local_psi, local_source

logical               :: include_radiation
real*8                :: Arad_bg, Brad_bg, Crad_bg, frad_bg
real*8                :: Te_eV, ne_SI, Lrad_imp, r_imp
real*8                :: Te_corr_eV, coef_rad_1, Sion_T, eta_Sp, ksiion, LradDcont_T
real*8                :: LradDrays_T, coef_ion_1, coef_ion_2, coef_ion_3, S_ion_puiss
real*8                :: r0_real8, rn0_real8
real*8                :: T0_corr, r0_corr, rn0_corr
integer               :: i_imp, offset_bgimp, i_bg     ! Loop for more than one background impurity
integer               :: i_proj
integer               :: i_psin, i_test, iimp(6), ineu(7), ibg_tot, i_pellet(2), i_flux(8), i_neo(10), i_boot(2)
integer               :: i_full(11), i_vec_B, i_vec_V, i_vec_E, i_vec_Jpol
integer, allocatable  :: iibg(:), iproj(:)
character*36          :: imp_label, proj_label


#ifdef WITH_Impurities
! See https://www.jorek.eu/wiki/doku.php?id=model500_501_555 for details
! Atomic physics coefficients:
!   -Mass ratio between main ions and impurites (m_i/m_imp)
real*8     :: m_i_over_m_imp
!   -Mean impurity ionization state
real*8     :: Z_imp, T0_Zimp, alpha_Zimp, n_imp
!   -Coefficients related to Z_imp
real*8     :: alpha_imp
real*8     :: beta_imp
!   - Effective charge state
real*8     :: Z_eff, eta_coef
!   -Radiation from injected impurities
real*8     :: Lrad
real*8     :: A0_rad, A1_rad, T1_rad, sig1_rad    ! Radiation rate parameters
real*8     :: A2_rad, T2_rad, sig2_rad
!   -Temporary variable for charge state distribution
real*8, allocatable :: P_imp(:)
real*8     :: E_ion
integer*8  :: ion_i, ion_k
#endif
real*8                :: T_real8, dLradDrays_dT, dLradDcont_dT, dSion_dT, dSrec_dT

real*8                :: psi_equi, psi_equi_s, psi_equi_t, psi_equi_R, psi_equi_Z
real*8                :: F_prof,   F_prof_s,   F_prof_t,   dF_dR,      dF_dZ,     dF_dpsi


!====================== --- Variables related to neutral density evolution (model 500 or 555)
logical               :: include_neutral_dens
real*8                :: IonN, RecN, AblN, coef_rec_1, Srec_T

integer, parameter :: nplot = 200
integer :: iplot, i_elm
real*8  :: stmp(200)
real*8  :: Rp_start, Zp_start, Rp_end, Zp_end
real*8  :: Rp, Zp, Rmin, Rmax, Zmin, Zmax, s_out, t_out, R_out, Z_out

namelist /vtk_params/ nsub, i_tor, i_plane, without_n0_mode, SI_units, &
                      include_fluxes, include_neo, include_magnetic_field, include_velocity_field,&
                      include_bootstrap, include_psi_norm, include_electric_field, include_Jpol, RphiZ_coords,&
                      include_projections, proj_basename


write(*,*) '***************************************'
write(*,*) '*       jorek2vtk                     *'
write(*,*) '***************************************'
write(*,*) ' if your VTK is smaller than expected,'
write(*,*) ' please consider the new parameters:'
write(*,*) '   -include_fluxes'
write(*,*) '   -include_neo'
write(*,*) '   -include_magnetic_field'
write(*,*) '   -include_velocity_field'
write(*,*) '   -include_electric_field'
write(*,*) '   -include_Jpol'
write(*,*) '   -include_bootstrap'
write(*,*) '   -include_psi_norm'
write(*,*) '   -include_projections'
write(*,*) '***************************************'

call flush_it(6)

allocate(node_list)
allocate(aux_node_list)
allocate(element_list)
allocate(bnd_elm_list)
allocate(bnd_node_list)

! --- Initialise input parameters and read the input namelist.
my_id     = 0
call initialise_parameters(my_id, "__NO_FILENAME__")

! --- Preset parameters
nsub                   = 5       ! Number of subdivisions of the cubic finite elements into linear pieces
i_tor                  = -1      ! If i_tor > 0, only this mode will be included in the vtk file...
i_plane                = 1       ! ... otherwise, all modes will be summed up at the toroidal plane i_plane
without_n0_mode        = .false. ! If true, do not include the n=0 mode (i_tor=1)
SI_units               = .false. ! when true, write variables in SI units
include_fluxes         = .false. ! include energy and density fluxes (or not)
include_neo            = .false. ! include neoclassical and more terms (or not)
include_magnetic_field = .false. ! include vector of magnetic field (or not)
include_velocity_field = .false. ! include vector of velocity field (or not)
include_electric_field = .false. ! include vector of E-field (or not), evaluated at t-dt/2 
include_Jpol           = .false. ! include poloidal current vector (J_phi=0 for visualization)
include_bootstrap      = .false. ! include bootstrap current and averaged current
include_psi_norm       = .true.  ! include normalized flux
include_projections    = .false. ! include projections from particles
proj_basename          = 'projections' ! basename for particle projection output files
RphiZ_coords           = .false. ! use xyz transformation (R,0,Z) instead of (R,Z,0)

include_radiation = .false. 
#if (defined WITH_Neutrals) || (defined WITH_Impurities)
include_radiation = .true.
include_neutral_dens = .true.
! --- Read ADAS data and generate coronal equilibrium if needed
call init_imp_adas(my_id)
#else
if (use_imp_adas .and. (nimp_bg(1) > 0.d0)) then
  call init_imp_adas(my_id)
  include_radiation = .true.
endif
#endif

! --- Read parameters from namelist file 'vtk.nml' if it exists
open(42, file='vtk.nml', action='read', status='old', iostat=ierr)
if ( ierr == 0 ) then
  write(*,*) 'Reading parameters from vtk.nml namelist.'
  read(42,vtk_params)
  close(42)
end if

write(*,*)
write(*,*) 'Parameters:'
write(*,*) '-----------'
write(*,*) 'nsub            =', nsub
write(*,*) 'i_tor           =', i_tor
write(*,*) 'i_plane         =', i_plane
write(*,*) 'without_n0_mode =', without_n0_mode
write(*,*) 'si_units        =', si_units
write(*,*) 'include_fluxes  =', include_fluxes
write(*,*) 'include_neo     =', include_neo
write(*,*) 'include_magnetic_field =',include_magnetic_field
write(*,*) 'include_velocity_field =',include_velocity_field
write(*,*) 'include_electric_field =',include_electric_field
write(*,*) 'include_Jpol      =', include_Jpol
write(*,*) 'include_bootstrap =', include_bootstrap
write(*,*) 'include_psi_norm  =', include_psi_norm
write(*,*) 'include_projections =', include_projections

if (include_projections) then
  write(*,*) ' -proj_basename =', trim(proj_basename)
end if

write(*,*) '-----------'
write(*,*) 'n_tor           =', n_tor
write(*,*) 'n_period        =', n_period
write(*,*) 'F0              =', F0
write(*,*) 'R_geo,Z_geo     =', R_geo, Z_geo
write(*,*)
call flush_it(6)

! --- Assign variable names in SI units
allocate(variable_names_si(0:n_var))
do i=1, n_var
  variable_names_si(i) = variable_names(i)
enddo

if ( SI_units ) then
#ifdef fullmhd
   variable_names_si(var_rho)='ne20_m-3    '
   variable_names_si(var_T  )='Te_keV      '
   variable_names_si(var_UR )='VR_km/s     '
   variable_names_si(var_UZ )='VZ_km/s     '
   variable_names_si(var_Up )='Vp_km/s     '
#else
   variable_names_si(var_u)='u_m/s       '
   variable_names_si(var_zj)='j_MA/m2     '
   variable_names_si(var_rho)='ne20_m-3    '
   if (with_TiTe) then
      variable_names_si(var_Ti)='Ti_keV      '
      variable_names_si(var_Te)='Te_keV      '
   else
      variable_names_si(var_T)='Te_keV      '
   endif
   if (with_Vpar) then
      variable_names_si(var_Vpar)='Vpar_km/s   '
   endif
#endif

#if (defined WITH_Neutrals) || (defined WITH_Impurities)
   variable_names_si(var_rhon)='N_dens_1d20  '
#endif

endif

! --- Add scalar quantities
n_scalars = 0
do i=1, n_var
  call add_vtk_entry(variable_names(i), variable_names_si(i), i_test, n_scalars, si_units, scalar_names) 
enddo

if (include_bootstrap) then
  if (.not. bootstrap) write(*,*)'VTK WARNING: if you want the bootstrap, please set bootstrap=.t. in your input file!'
  call add_vtk_entry('j_bootstrap ', 'j_b_MA/m2   ',   i_boot( 1), n_scalars, si_units, scalar_names) 
  call add_vtk_entry('j_averaged  ', 'j_av_MA/m2  ',   i_boot( 2), n_scalars, si_units, scalar_names) 
endif

if (include_neo) then
  call add_vtk_entry('Er          ', 'Er_kV/m     ',    i_neo( 1), n_scalars, si_units, scalar_names) 
  call add_vtk_entry('Vtheta      ', 'Vtheta_km/s ',    i_neo( 2), n_scalars, si_units, scalar_names) 
  call add_vtk_entry('Mach_par    ', 'Mach_par    ',    i_neo( 3), n_scalars, si_units, scalar_names) 
  call add_vtk_entry('Mach_pol    ', 'Mach_pol    ',    i_neo( 4), n_scalars, si_units, scalar_names) 
  call add_vtk_entry('Vsound      ', 'Vsound_km/s ',    i_neo( 5), n_scalars, si_units, scalar_names) 
  call add_vtk_entry('Btot        ', 'Btot_T      ',    i_neo( 6), n_scalars, si_units, scalar_names) 
  call add_vtk_entry('Vneo        ', 'Vneo_km/s   ',    i_neo( 7), n_scalars, si_units, scalar_names) 
  call add_vtk_entry('Vperp_e     ', 'Vperp_e_km/s',    i_neo( 8), n_scalars, si_units, scalar_names) 
  call add_vtk_entry('ki_neo      ', 'ki_neo      ',    i_neo( 9), n_scalars, si_units, scalar_names) 
  call add_vtk_entry('mu_neo      ', 'mu_neo      ',    i_neo(10), n_scalars, si_units, scalar_names) 
endif

if (include_fluxes) then
  call add_vtk_entry('pressure    ', 'pressure    ',    i_flux(1), n_scalars, si_units, scalar_names) 
  call add_vtk_entry('E_flux_Kpar ', 'E_flux_Kpar ',    i_flux(2), n_scalars, si_units, scalar_names) 
  call add_vtk_entry('E_flux_kperp', 'E_flux_kperp',    i_flux(3), n_scalars, si_units, scalar_names) 
  call add_vtk_entry('E_flux_Vpar ', 'E_flux_Vpar ',    i_flux(4), n_scalars, si_units, scalar_names) 
  call add_vtk_entry('E_flux_Vperp', 'E_flux_Vperp',    i_flux(5), n_scalars, si_units, scalar_names) 
  call add_vtk_entry('D_flux_Dperp', 'D_flux_Dperp',    i_flux(6), n_scalars, si_units, scalar_names) 
  call add_vtk_entry('D_flux_Vpar ', 'D_flux_Vpar ',    i_flux(7), n_scalars, si_units, scalar_names) 
  call add_vtk_entry('D_flux_Vperp', 'D_flux_Vperp',    i_flux(8), n_scalars, si_units, scalar_names) 
endif

if (use_pellet) then
  call add_vtk_entry('Pressure    ', 'Pressure    ',  i_pellet(1), n_scalars, si_units, scalar_names) 
  call add_vtk_entry('Pellet      ', 'Pellet      ',  i_pellet(2), n_scalars, si_units, scalar_names) 
endif

if (include_psi_norm) then
  call add_vtk_entry('psi_norm    ', 'psi_norm    ',    i_psin, n_scalars, si_units, scalar_names) 
endif

allocate(iibg(n_adas),iproj(n_var))

#if (defined WITH_Neutrals) && (!defined WITH_Impurities)
  if (include_neutral_dens) then
    call add_vtk_entry('IonN_s-1    ', 'IonN_s-1    ',  ineu(6), n_scalars, si_units, scalar_names) 
    call add_vtk_entry('RecN_s-1    ', 'RecN_s-1    ',  ineu(7), n_scalars, si_units, scalar_names) 
  endif
#endif

if (include_radiation) then

#if (defined WITH_Neutrals) && (!defined WITH_Impurities)
  call add_vtk_entry('Neut_Ionis  ', 'NeutIonis_Wm-3',  ineu(1), n_scalars, si_units, scalar_names) 
  call add_vtk_entry('Lin_rad     ', 'Lin_radWm-3 ',    ineu(2), n_scalars, si_units, scalar_names) 
  call add_vtk_entry('Brems       ', 'Brems_Wm-3  ',    ineu(3), n_scalars, si_units, scalar_names) 
  call add_vtk_entry('Joule       ', 'Joule_Wm-3  ',    ineu(4), n_scalars, si_units, scalar_names) 
#endif

if (include_projections) then
  do i = 1, n_var
    write(proj_label, '(a4,i2.2)') 'aux_', i
    call add_vtk_entry(proj_label, proj_label, iproj(i), n_scalars, si_units, scalar_names) 
  end do
end if

#ifdef WITH_Impurities
  call add_vtk_entry('Ionis       ', 'Ionis_Jm-3  ',    iimp(1), n_scalars, si_units, scalar_names) 
  call add_vtk_entry('Coronal_rad ', 'Coronal_radWm-3', iimp(2), n_scalars, si_units, scalar_names) 
  call add_vtk_entry('Joule       ', 'Joule_Wm-3  ',    iimp(3), n_scalars, si_units, scalar_names) 
  call add_vtk_entry('Z_imp       ', 'Z_imp       ',    iimp(4), n_scalars, si_units, scalar_names) 
  call add_vtk_entry('Z_eff       ', 'Z_eff       ',    iimp(5), n_scalars, si_units, scalar_names) 
#endif

  ! --- Background radiation
  do i = 1,n_adas
    imp_label = 'Imp_bg_'//trim(imp_type(i))//'_Wm-3'
    call add_vtk_entry(imp_label, imp_label,  iibg(i), n_scalars, si_units, scalar_names) 
  end do

  call add_vtk_entry('Imp_bg      ', 'Imp_bg_Wm-3 ',  ibg_tot, n_scalars, si_units, scalar_names)  ! total bg radiation 

endif

#ifdef fullmhd
  call add_vtk_entry('B_R         ', 'B_R         ',    i_full( 1), n_scalars, si_units, scalar_names) 
  call add_vtk_entry('B_Z         ', 'B_Z         ',    i_full( 2), n_scalars, si_units, scalar_names) 
  call add_vtk_entry('B_phi       ', 'B_phi       ',    i_full( 3), n_scalars, si_units, scalar_names) 
  call add_vtk_entry('J_R         ', 'J_R         ',    i_full( 4), n_scalars, si_units, scalar_names) 
  call add_vtk_entry('J_Z         ', 'J_Z         ',    i_full( 5), n_scalars, si_units, scalar_names) 
  call add_vtk_entry('J_phi       ', 'J_phi       ',    i_full( 6), n_scalars, si_units, scalar_names) 
  call add_vtk_entry('FFprime_equi', 'FFprime_equi',    i_full( 7), n_scalars, si_units, scalar_names) 
  call add_vtk_entry('F_prof_equi ', 'F_prof_equi ',    i_full( 8), n_scalars, si_units, scalar_names) 
  call add_vtk_entry('Grad_P      ', 'Grad_P      ',    i_full( 9), n_scalars, si_units, scalar_names) 
  call add_vtk_entry('JxB         ', 'JxB         ',    i_full(10), n_scalars, si_units, scalar_names) 
  call add_vtk_entry('V_parallel  ', 'V_par_km/s  ',    i_full(11), n_scalars, si_units, scalar_names) 
#endif /*fullmhd*/
! --- end adding scalars


! --- Add vectors
n_vectors   = 0

if (include_magnetic_field) then
  call add_vtk_entry('B_field     ', 'B_field     ',    i_vec_B, n_vectors, si_units, vector_names) 
endif

if (include_velocity_field) then
  call add_vtk_entry('v_field     ', 'v_field     ',    i_vec_V, n_vectors, si_units, vector_names) 
endif

if (include_electric_field) then
  call add_vtk_entry('E_field_tmid', 'E_field_tmid',    i_vec_E, n_vectors, si_units, vector_names) 
endif

if (include_Jpol) then
  call add_vtk_entry('Jpol        ', 'Jpol (MA/m2)', i_vec_Jpol, n_vectors, si_units, vector_names) 
endif
! --- end adding vectors



do k_tor=1, n_tor
  mode(k_tor) = + int(k_tor / 2) * n_period
enddo

if (include_projections) then
  filename_proj = trim(proj_basename)//'_restart.h5' ! only hdf5 format supported for particle projections
  call import_hdf5_restart_aux(aux_node_list, filename_proj, rst_format, ierr)
  if (ierr .ne. 0) then
    write(*,*) 'ERROR: Cannot find projection restart file. Check if proj_basename is set properly.'
    stop
  end if
end if

call import_restart(node_list, element_list, 'jorek_restart', rst_format, ierr, .true.)

call initialise_basis                              ! define the basis functions at the Gaussian points

nnos = nsub*nsub*element_list%n_elements
allocate(currdens(nnos),xyz(3,nnos),scalars(nnos,1:n_scalars),vectors(nnos,3,1:n_vectors))
currdens = 0.

nnoel = 4
nel   = (nsub-1)*(nsub-1)*element_list%n_elements
allocate(ien(nnoel,nel))

inode   = 0
ielm    = 0
scalars = 0.d0
vectors = 0.d0
xyz     = 0
ien     = 0

call boundary_from_grid(node_list, element_list, bnd_node_list, bnd_elm_list, .false.)

minRad = 0.0
if (bootstrap) then
  call bootstrap_find_minRad(node_list, element_list, ES%R_axis, ES%Z_axis, ES%psi_axis, ES%psi_bnd)
  call bootstrap_get_q_and_ft_splines(node_list, element_list, ES%psi_axis, ES%psi_xpoint, ES%R_xpoint, ES%Z_xpoint)
  call bootstrap_get_averaged_j_spline(node_list, element_list, ES%psi_axis, ES%psi_xpoint, ES%R_xpoint, ES%Z_xpoint)
endif

grad_psi = 0.d0

! --- You may choose to print your poloidal snapshot at a different toroidal angle
toroidal_angle = 0.d0 ! 2*PI / 6
if (toroidal_angle .ne. 0.d0) then
  do k_tor=1, n_tor
    mode(k_tor) = + int(k_tor / 2) * n_period
  enddo
  HZ(1,i_plane)   = 1.d0
  do i=1,(n_tor-1)/2
    HZ(2*i,i_plane)      = cos(mode(2*i)   * toroidal_angle )
    HZ(2*i+1,i_plane)    = sin(mode(2*i+1) * toroidal_angle )
  enddo
endif

do i=1,element_list%n_elements

   ! if(element_list%element(i)%n_sons.eq.0) then

  do j=1,nsub

    s = float(j-1)/float(nsub-1)

    do k=1,nsub

      t = float(k-1)/float(nsub-1)

      call interp_RZ(node_list,element_list,i,s,t,R,R_s,R_t,R_st,R_ss,R_tt,Z,Z_s,Z_t,Z_st,Z_ss,Z_tt)

      xjac  = R_s * Z_t - R_t * Z_s

      if ( xjac == 0.d0 ) xjac = 1.d-8

      BigR  = R

      xjac_x  = (R_ss*Z_t**2 - Z_ss*R_t*Z_t - 2.d0*R_st*Z_s*Z_t   &
              + Z_st*(R_s*Z_t + R_t*Z_s) + R_tt*Z_s**2 - Z_tt*R_s*Z_s) / xjac

      xjac_y  = (Z_tt*R_s**2 - R_tt*Z_s*R_s - 2.d0*Z_st*R_t*R_s   &
              + R_st*(Z_t*R_s + Z_s*R_t) + Z_ss*R_t**2 - R_ss*Z_t*R_t) / xjac

      inode = inode+1
      
      if (RphiZ_coords) then
        xyz(1:3,inode) = (/ R, 0.d0,    Z /)
      else
        xyz(1:3,inode) = (/ R,    Z, 0.d0 /)
      endif
      !====================== --- specific for axisymmetric quantities
      ! Put here all quantities that are axisymmetric (n=0 mode only) and should not be summed
      ! over all harmonics: for instance, to compute Vtheta, Er, Vneo, etc.
      ! ===> this corresponds to forcing i_tor = 1 (thus n=0 only)

      ! save old values
      i_tor_old = i_tor
      i_tor     = 1
      ! compute all derivatives, as in loop below
      if ( (xjac .gt. 1.d-6) .and. (jorek_model .ge. 100) ) then

#ifndef fullmhd
        call interp(node_list,element_list,i,var_psi,i_tor,s,t,Ps0,Ps0_s,Ps0_t,Ps0_st,Ps0_ss,Ps0_tt)
        call interp(node_list,element_list,i,var_u,  i_tor,s,t,U0, U0_s, U0_t, U0_st, U0_ss, U0_tt)
        call interp(node_list,element_list,i,var_zj, i_tor,s,t,ZJ0,ZJ0_s,ZJ0_t,ZJ0_st,ZJ0_ss,ZJ0_tt)
        call interp(node_list,element_list,i,var_w,  i_tor,s,t,W0, W0_s, W0_t, W0_st, W0_ss, W0_tt)
        call interp(node_list,element_list,i,var_rho,i_tor,s,t,ZN0,ZN0_s,ZN0_t,ZN0_st,ZN0_ss,ZN0_tt)

        if (with_Vpar) then
          call interp(node_list,element_list,i,var_Vpar,i_tor,s,t,V0,V0_s,V0_t,V0_st,V0_ss,V0_tt)
        else
          V0=0; V0_s=0; V0_t=0; V0_st=0; V0_ss=0; V0_tt=0
        end if

        if (with_TiTe) then
           call interp(node_list,element_list,i,var_Te,  i_tor,s,t,Te0, Te0_s, Te0_t, Te0_st, Te0_ss, Te0_tt)
           call interp(node_list,element_list,i,var_Ti,  i_tor,s,t,Ti0, Ti0_s, Ti0_t, Ti0_st, Ti0_ss, Ti0_tt)
           T0   = Ti0   + Te0
           T0_s = Ti0_s + Te0_s
           T0_t = Ti0_t + Te0_t
        else
          call interp(node_list,element_list,i,var_T,  i_tor,s,t,T0, T0_s, T0_t, T0_st, T0_ss, T0_tt)
          Te0    = T0  /2.d0;     Ti0    = T0  /2.d0
          Te0_s  = T0_s/2.d0;     Ti0_s  = T0_s/2.d0
          Te0_t  = T0_t/2.d0;     Ti0_t  = T0_t/2.d0
        endif

        u0_x   = (   Z_t * U0_s  - Z_s * U0_t ) / xjac
        u0_y   = ( - R_t * U0_s  + R_s * U0_t ) / xjac

        ps0_x  = (   Z_t * Ps0_s - Z_s * Ps0_t ) / xjac
        ps0_y  = ( - R_t * Ps0_s + R_s * Ps0_t ) / xjac

        T0_x   = (   Z_t * T0_s  - Z_s * T0_t ) / xjac
        T0_y   = ( - R_t * T0_s  + R_s * T0_t ) / xjac

        zj0_x  = (   Z_t * ZJ0_s - Z_s * ZJ0_t ) / xjac
        zj0_y  = ( - R_t * ZJ0_s + R_s * ZJ0_t ) / xjac

        zn0_x  = (   Z_t * zn0_s - Z_s * zn0_t ) / xjac
        zn0_y  = ( - R_t * zn0_s + R_s * zn0_t ) / xjac
#endif
        if (include_neo) then

#ifdef fullmhd 
          ! not yet implemented in model710
          scalars(inode,i_neo(1))  = Er
          scalars(inode,i_neo(2))  = Vtheta
          scalars(inode,i_neo(3))  = Mach_par
          scalars(inode,i_neo(4))  = Mach_pol
          scalars(inode,i_neo(5))  = Vsound
          scalars(inode,i_neo(6))  = Btot
          scalars(inode,i_neo(7))  = Vneo
          scalars(inode,i_neo(8))  = Vperp_e
          scalars(inode,i_neo(9))  = aki_neo_node
          scalars(inode,i_neo(10)) = amu_neo_node
#else /* not full-MHD */

            !*** compute diagnostics ***
          psi_abs = sqrt(ps0_x*ps0_x + ps0_y * ps0_y)
          Btheta  = (psi_abs/R)
          Vtheta  = 0.d0
          Vperp_e = 0.0
          Vneo    = 0.d0
          Er      = 0.d0
          mach_par= 0.d0
          mach_pol= 0.d0
          vsound  = 0.d0
          amu_neo_node = 0.d0
          aki_neo_node = 0.d0

          if ((psi_abs .gt. 1.d-6) .and. (ZN0.gt.1.d-6) .and. (abs(Btheta).gt.1.d-6)) then

            Vtheta  = -1./Btheta*((u0_x + tauIC/ZN0*(T0_x*ZN0 + ZN0_x*T0))*ps0_x  + &
                                  (u0_y + tauIC/ZN0*(T0_y*ZN0 + ZN0_y*T0))*ps0_y) + V0*Btheta

            Vperp_e = -1./Btheta*((u0_x - tauIC/ZN0*(T0_x*ZN0 + ZN0_x*T0))*ps0_x  + &
                                  (u0_y - tauIC/ZN0*(T0_y*ZN0 + ZN0_y*T0))*ps0_y)

            if (NEO) then
!                num_neo_file= ( neo_file /= 'none')
              if (num_neo_file) then
!                  write(*,*) 'neo_file=',neo_file
!                  write(*,*) 'using ki and mui profiles from file "'//trim(neo_file)//'"'
                call neo_coef( xpoint, xcase, Z, ES%Z_xpoint, Ps0 ,ES%psi_axis, ES%psi_bnd, &
                               amu_neo_node, aki_neo_node)
                Vneo   = aki_neo_node / Btheta*tauIC  * (ps0_x*T0_x + ps0_y*T0_y)
              else
                Vneo   = aki_neo_const / Btheta*tauIC * (ps0_x*T0_x + ps0_y*T0_y)
              endif

            endif  ! NEO

            Er       = -(U0_x * ps0_x + U0_y * ps0_y)/psi_abs         ! radial electric field
            Btot     = sqrt(F0**2 + ps0_x**2 + ps0_y**2) / BigR       ! total magnetic field (equilibrium)
            Vsound   = sqrt(GAMMA*T0)/Btot                            ! sound speed
            Mach_par = V0/Vsound                                      ! parallel Mach number
            Mach_pol = Vtheta/Vsound                                  ! poloidal Mach number

          endif !psi_abs

          ! save those specific values of axisymmetric parameters
          if (grad_psi .ne. 0.d0) then
            scalars(inode,i_neo(1)) = Er
            scalars(inode,i_neo(2)) = Vtheta
            scalars(inode,i_neo(3)) = Mach_par
            scalars(inode,i_neo(4)) = Mach_pol
            scalars(inode,i_neo(5)) = Vsound
            scalars(inode,i_neo(6)) = Btot
            scalars(inode,i_neo(7)) = Vneo
            scalars(inode,i_neo(8)) = Vperp_e
            if (NEO) then
               if (num_neo_file) then
                  scalars(inode,i_neo(9))  = aki_neo_node
                  scalars(inode,i_neo(10)) = amu_neo_node
               else
                  scalars(inode,i_neo(9))  = aki_neo_const
                  scalars(inode,i_neo(10)) = amu_neo_const
               endif
            endif   ! NEO

          endif     ! grad_psi
#endif /* non-full-MHD part */

        endif       ! include_neo

      endif         ! xjac

      ! old values back to normal
      i_tor = i_tor_old

      !====================== --- specific for NON-axisymmetric quantities
      ! 2 cases, depending on the value of i_tor chosen
      if ((i_tor .ge. 1) .and. (i_tor .le. n_tor)) then

        do m=1,n_var
          call interp(node_list,element_list,i,m,i_tor,s,t,P,P_s,P_t,P_st,P_ss,P_tt)
          scalars(inode,m) = P * HZ(i_tor,i_plane)
        enddo

        if (include_projections) then
          do i_proj=1,n_var
            call interp(aux_node_list,element_list,i,i_proj,i_tor,s,t,P,P_s,P_t,P_st,P_ss,P_tt)
            scalars(inode,iproj(i_proj)) = P * HZ(i_tor,i_plane)
          end do
        end if

        if (jorek_model .lt. 100) cycle
        
        ! The real current density
        currdens(inode) = -scalars(inode,3)/BigR

        if ((xjac .gt. 1.d-6)) then
#ifndef fullmhd
          call interp(node_list,element_list,i,var_psi,i_tor,s,t,Psi,Ps_s,Ps_t,Ps_st,Ps_ss,Ps_tt)
          call interp(node_list,element_list,i,var_u,  i_tor,s,t,U,U_s,U_t,U_st,U_ss,U_tt)
          call interp(node_list,element_list,i,var_zj, i_tor,s,t,ZJ,ZJ_s,ZJ_t,ZJ_st,ZJ_ss,ZJ_tt)
          call interp(node_list,element_list,i,var_w,  i_tor,s,t,W,W_s,W_t,W_st,W_ss,W_tt)
          call interp(node_list,element_list,i,var_rho,i_tor,s,t,RHO,RHO_s,RHO_t,RHO_st,RHO_ss,RHO_tt)
          if (with_TiTe) then
             call interp(node_list,element_list,i,var_Ti,i_tor,s,t,Ti,Ti_s,Ti_t,Ti_st,Ti_ss,Ti_tt)
             call interp(node_list,element_list,i,var_Te,i_tor,s,t,Te,Te_s,Te_t,Te_st,Te_ss,Te_tt)
             TT   = Ti   + Te
             TT_s = Ti_s + Te_s
             TT_t = Ti_t + Te_t
          else
             call interp(node_list,element_list,i,var_T,  i_tor,s,t,TT ,TT_s, TT_t, TT_st, TT_ss, TT_tt)
          endif
 
          if (with_Vpar) then
            call interp(node_list,element_list,i,var_Vpar,i_tor,s,t,V,V_s,V_t,V_st,V_ss,V_tt)
          else
            V=0; V_s=0; V_t=0; V_st=0; V_ss=0; V_tt=0
          end if

          u_x   = (   Z_t * U_s  - Z_s * U_t ) / xjac
          u_y   = ( - R_t * U_s  + R_s * U_t ) / xjac

          ps_x  = (   Z_t * PS_s - Z_s * PS_t ) / xjac
          ps_y  = ( - R_t * PS_s + R_s * PS_t ) / xjac

          TT_x  = (   Z_t * TT_s - Z_s * TT_t ) / xjac
          TT_y  = ( - R_t * TT_s + R_s * TT_t ) / xjac

          zj_x  = (   Z_t * ZJ_s - Z_s * ZJ_t ) / xjac
          zj_y  = ( - R_t * ZJ_s + R_s * ZJ_t ) / xjac

           !*** compute diagnostics ***
          v_perp  = R * sqrt(u_x*u_x + u_y*u_y)

          psi_J = (Ps_s * ZJ_t - PS_t * ZJ_s ) / xjac
          R_p   = (2.d0 * R * (R_s * (RHO_t * TT + RHO * TT_t) - R_t * (RHO_s * TT + RHO * TT_s) )) / xjac
          error = psi_J - R_p  ! "error" in Grad_Shafranov equilibrium force balance
#endif
        endif  ! xjac check

      else  ! i_tor

#ifdef fullmhd

        scalars(inode,i_full(1):i_full(size(i_full))) = 0.d0
        
        A3    = 0.d0
        A3_R  = 0.d0  ;  AR_R  = 0.d0  ;  AZ_R  = 0.d0  ;  rho   = 0.d0  ;  TT    = 0.d0
        A3_Z  = 0.d0  ;  AR_Z  = 0.d0  ;  AZ_Z  = 0.d0  ;  rho_x = 0.d0  ;  TT_x  = 0.d0
        A3_RR = 0.d0  ;  AR_RR = 0.d0  ;  AZ_RR = 0.d0  ;  rho_y = 0.d0  ;  TT_y  = 0.d0
        A3_ZZ = 0.d0  ;  AR_ZZ = 0.d0  ;  AZ_ZZ = 0.d0  ;  rho_p = 0.d0  ;  TT_p  = 0.d0
        A3_RZ = 0.d0  ;  AR_RZ = 0.d0  ;  AZ_RZ = 0.d0
        A3_p  = 0.d0  ;  AR_p  = 0.d0  ;  AZ_p  = 0.d0
        A3_Rp = 0.d0  ;  AR_Rp = 0.d0  ;  AZ_Rp = 0.d0
        A3_Zp = 0.d0  ;  AR_Zp = 0.d0  ;  AZ_Zp = 0.d0
        A3_pp = 0.d0  ;  AR_pp = 0.d0  ;  AZ_pp = 0.d0

        VR    = 0.d0  ;  VZ    = 0.d0  ;  Vp    = 0.d0

        do i_tor = 1, n_tor

          if ( ( i_tor == 1 ) .and. ( without_n0_mode ) ) cycle ! Do not include the n=0 mode

          do m=1,n_var
             call interp(node_list,element_list,i,m,i_tor,s,t,P,P_s,P_t,P_st,P_ss,P_tt)
             scalars(inode,m) = scalars(inode,m) + P * HZ(i_tor,i_plane)
          enddo
          
          call interp(node_list,element_list,i,var_AR, i_tor,s,t,AR0,AR0_s,AR0_t,AR0_st,AR0_ss,AR0_tt)
          call interp(node_list,element_list,i,var_AZ, i_tor,s,t,AZ0,AZ0_s,AZ0_t,AZ0_st,AZ0_ss,AZ0_tt)
          call interp(node_list,element_list,i,var_A3, i_tor,s,t,A30,A30_s,A30_t,A30_st,A30_ss,A30_tt)
          call interp(node_list,element_list,i,var_uR, i_tor,s,t,VR0,VR0_s,VR0_t,VR0_st,VR0_ss,VR0_tt)
          call interp(node_list,element_list,i,var_uZ, i_tor,s,t,VZ0,VZ0_s,VZ0_t,VZ0_st,VZ0_ss,VZ0_tt)
          call interp(node_list,element_list,i,var_uP, i_tor,s,t,VP0,VP0_s,VP0_t,VP0_st,VP0_ss,VP0_tt)
          call interp(node_list,element_list,i,var_rho,i_tor,s,t,ZN0,ZN0_s,ZN0_t,ZN0_st,ZN0_ss,ZN0_tt)

          if (with_TiTe) then
             call interp(node_list,element_list,i,var_Te,  i_tor,s,t,Te0, Te0_s, Te0_t, Te0_st, Te0_ss, Te0_tt)
             call interp(node_list,element_list,i,var_Ti,  i_tor,s,t,Ti0, Ti0_s, Ti0_t, Ti0_st, Ti0_ss, Ti0_tt)
             T0   = Ti0   + Te0
             T0_s = Ti0_s + Te0_s
             T0_t = Ti0_t + Te0_t
          else
            call interp(node_list,element_list,i,var_T,  i_tor,s,t,T0, T0_s, T0_t, T0_st, T0_ss, T0_tt)
            Te0    = T0  /2.d0;     Ti0    = T0  /2.d0
            Te0_s  = T0_s/2.d0;     Ti0_s  = T0_s/2.d0
            Te0_t  = T0_t/2.d0;     Ti0_t  = T0_t/2.d0
          endif


          if (i_tor == 1) then
            call interp(node_list,element_list,i,710,i_tor,s,t,F_prof  ,F_prof_s  ,F_prof_t  ,W_st,W_ss,W_tt)
            call interp(node_list,element_list,i,711,i_tor,s,t,psi_equi,psi_equi_s,psi_equi_t,W_st,W_ss,W_tt)
          endif

          AR_p  = AR_p  + AR0 * HZ_p(i_tor,i_plane)
          AR_pp = AR_pp + AR0 * HZ_pp(i_tor,i_plane)

          AZ_p  = AZ_p  + AZ0 * HZ_p(i_tor,i_plane)
          AZ_pp = AZ_pp + AZ0 * HZ_pp(i_tor,i_plane)

          A3    = A3    + A30 * HZ(i_tor,i_plane)
          A3_p  = A3_p  + A30 * HZ_p(i_tor,i_plane)
          A3_pp = A3_pp + A30 * HZ_pp(i_tor,i_plane)

          VR    = VR   + VR0 * HZ(i_tor,i_plane)
          VZ    = VZ   + VZ0 * HZ(i_tor,i_plane)
          Vp    = Vp   + Vp0 * HZ(i_tor,i_plane)

          TT    = TT   + T0 * HZ(i_tor,i_plane)
          TT_p  = TT_p + T0 * HZ_p(i_tor,i_plane)

          rho   = rho   + zn0* HZ(i_tor,i_plane)
          rho_p = rho_p + zn0* HZ_p(i_tor,i_plane)

          if ((xjac .gt. 1.d-6)) then  ! avoid the axis

            AR_R  = AR_R  + (   Z_t * AR0_s - Z_s * AR0_t ) / xjac * HZ(i_tor,i_plane)
            AR_Z  = AR_Z  + ( - R_t * AR0_s + R_s * AR0_t ) / xjac * HZ(i_tor,i_plane)
            AR_RR = AR_RR + ( (AR0_ss * Z_t**2 - 2.d0*AR0_st * Z_s*Z_t + AR0_tt * Z_s**2  &
                             + AR0_s * (Z_st*Z_t - Z_tt*Z_s )                             &
                             + AR0_t * (Z_st*Z_s - Z_ss*Z_t ) )    / xjac**2              &
                             - xjac_x * (AR0_s* Z_t - AR0_t * Z_s)  / xjac**2             ) * HZ(i_tor,i_plane)
            AR_ZZ = AR_ZZ + ( (AR0_ss * R_t**2 - 2.d0*AR0_st * R_s*R_t + AR0_tt * R_s**2  &
                             + AR0_s * (R_st*R_t - R_tt*R_s )                             &
                             + AR0_t * (R_st*R_s - R_ss*R_t ) )    / xjac**2              &
                             - xjac_y * (- AR0_s * R_t + AR0_t * R_s )  / xjac**2         ) * HZ(i_tor,i_plane)
            AR_RZ = AR_RZ + ( (- AR0_ss * Z_t*R_t - AR0_tt * R_s*Z_s                        &
                               + AR0_st * (Z_s*R_t  + Z_t*R_s  )                            &
                               - AR0_s  * (R_st*Z_t - R_tt*Z_s )                            &
                               - AR0_t * (R_st*Z_s  - R_ss*Z_t ) )  / xjac**2               &
                               - xjac_x * (- AR0_s * R_t + AR0_t * R_s )   / xjac**2        ) * HZ(i_tor,i_plane)
            AR_Rp = AR_Rp + (   Z_t * AR0_s - Z_s * AR0_t ) / xjac * HZ_p(i_tor,i_plane)
            AR_Zp = AR_Zp + ( - R_t * AR0_s + R_s * AR0_t ) / xjac * HZ_p(i_tor,i_plane)

            AZ_R  = AZ_R  + (   Z_t * AZ0_s - Z_s * AZ0_t ) / xjac * HZ(i_tor,i_plane)
            AZ_Z  = AZ_Z  + ( - R_t * AZ0_s + R_s * AZ0_t ) / xjac * HZ(i_tor,i_plane)
            AZ_RR = AZ_RR + ( (AZ0_ss * Z_t**2 - 2.d0*AZ0_st * Z_s*Z_t + AZ0_tt * Z_s**2  &
                             + AZ0_s * (Z_st*Z_t - Z_tt*Z_s )                             &
                             + AZ0_t * (Z_st*Z_s - Z_ss*Z_t ) )    / xjac**2              &
                             - xjac_x * (AZ0_s* Z_t - AZ0_t * Z_s)  / xjac**2             ) * HZ(i_tor,i_plane)
            AZ_ZZ = AZ_ZZ + ( (AZ0_ss * R_t**2 - 2.d0*AZ0_st * R_s*R_t + AZ0_tt * R_s**2  &
                             + AZ0_s * (R_st*R_t - R_tt*R_s )                             &
                             + AZ0_t * (R_st*R_s - R_ss*R_t ) )    / xjac**2              &
                             - xjac_y * (- AZ0_s * R_t + AZ0_t * R_s )  / xjac**2         ) * HZ(i_tor,i_plane)
            AZ_RZ = AZ_RZ + ( (- AZ0_ss * Z_t*R_t - AZ0_tt * R_s*Z_s                        &
                               + AZ0_st * (Z_s*R_t  + Z_t*R_s  )                            &
                               - AZ0_s  * (R_st*Z_t - R_tt*Z_s )                            &
                               - AZ0_t * (R_st*Z_s  - R_ss*Z_t ) )  / xjac**2               &
                               - xjac_x * (- AZ0_s * R_t + AZ0_t * R_s )   / xjac**2        ) * HZ(i_tor,i_plane)
            AZ_Rp = AZ_Rp + (   Z_t * AZ0_s - Z_s * AZ0_t ) / xjac * HZ_p(i_tor,i_plane)
            AZ_Zp = AZ_Zp + ( - R_t * AZ0_s + R_s * AZ0_t ) / xjac * HZ_p(i_tor,i_plane)

            A3_R  = A3_R  + (   Z_t * A30_s - Z_s * A30_t ) / xjac * HZ(i_tor,i_plane)
            A3_Z  = A3_Z  + ( - R_t * A30_s + R_s * A30_t ) / xjac * HZ(i_tor,i_plane)
            A3_RR = A3_RR + ( (A30_ss * Z_t**2 - 2.d0*A30_st * Z_s*Z_t + A30_tt * Z_s**2  &
                             + A30_s * (Z_st*Z_t - Z_tt*Z_s )                             &
                             + A30_t * (Z_st*Z_s - Z_ss*Z_t ) )    / xjac**2              &
                             - xjac_x * (A30_s* Z_t - A30_t * Z_s)  / xjac**2             ) * HZ(i_tor,i_plane)
            A3_ZZ = A3_ZZ + ( (A30_ss * R_t**2 - 2.d0*A30_st * R_s*R_t + A30_tt * R_s**2  &
                             + A30_s * (R_st*R_t - R_tt*R_s )                             &
                             + A30_t * (R_st*R_s - R_ss*R_t ) )    / xjac**2              &
                             - xjac_y * (- A30_s * R_t + A30_t * R_s )  / xjac**2         ) * HZ(i_tor,i_plane)
            A3_RZ = A3_RZ + ( (- A30_ss * Z_t*R_t - A30_tt * R_s*Z_s                        &
                               + A30_st * (Z_s*R_t  + Z_t*R_s  )                            &
                               - A30_s  * (R_st*Z_t - R_tt*Z_s )                            &
                               - A30_t * (R_st*Z_s  - R_ss*Z_t ) )  / xjac**2               &
                               - xjac_x * (- A30_s * R_t + A30_t * R_s )   / xjac**2        ) * HZ(i_tor,i_plane)
            A3_Rp = A3_Rp + (   Z_t * A30_s - Z_s * A30_t ) / xjac * HZ_p(i_tor,i_plane)
            A3_Zp = A3_Zp + ( - R_t * A30_s + R_s * A30_t ) / xjac * HZ_p(i_tor,i_plane)

            TT_x  = TT_x + (   Z_t * T0_s - Z_s * T0_t )   / xjac * HZ(i_tor,i_plane)
            TT_y  = TT_y + ( - R_t * T0_s + R_s * T0_t )   / xjac * HZ(i_tor,i_plane)

            rho_x = rho_x  + (   Z_t * ZN0_s - Z_s * ZN0_t ) / xjac * HZ(i_tor,i_plane)
            rho_y = rho_y  + ( - R_t * ZN0_s + R_s * ZN0_t ) / xjac * HZ(i_tor,i_plane)
          endif

        enddo  ! end loop toroidal harmonics

        ! --- Derivatives of F-profile
        dF_dR      = (   Z_t * F_prof_s   - Z_s * F_prof_t   ) / xjac
        dF_dZ      = ( - R_t * F_prof_s   + R_s * F_prof_t   ) / xjac
        psi_equi_R = (   Z_t * psi_equi_s - Z_s * psi_equi_t ) / xjac
        psi_equi_Z = ( - R_t * psi_equi_s + R_s * psi_equi_t ) / xjac
        dF_dpsi = (dF_dR*psi_equi_R + dF_dZ*psi_equi_Z) / max(1.d-13,(psi_equi_R**2 + psi_equi_Z**2))

        BR = ( A3_Z - AZ_p )/ R
        BZ = ( AR_p - A3_R )/ R 
        BP = ( AZ_R - AR_Z )    + F_prof/ R

        ZKpar_T = ZK_par * ((max(TT, T_min_ZKpar ))/T_0)**2.5

        BR_R = -1/R**2 * ( A3_Z - AZ_p ) + ( A3_RZ - AZ_Rp )/R
        BR_Z = ( A3_ZZ - AZ_Zp )/R
        BR_p = ( A3_Zp - AZ_pp )/R
        BZ_R = -1/R**2 * ( AR_p - A3_R ) + ( AR_Rp - A3_RR )/R
        BZ_Z = ( AR_Zp - A3_RZ )/R
        BZ_p = ( AR_pp - A3_Rp )/R
        BP_R = ( AZ_RR - AR_RZ ) + dF_dpsi*A3_R/R - F_prof/R**2
        BP_Z = ( AZ_RZ - AR_ZZ ) + dF_dpsi*A3_Z/R
        BP_p = ( AZ_Rp - AR_Zp )

        J_R   = (R*BP_Z - BZ_p) / R
        J_Z   = (BR_p - R*BP_R - BP) / R
        J_phi = (BZ_R - BR_Z) ! defined as the physical component J_phi * e_phi, like Bp and Vp

        JxB_R      = J_Z  *BP - J_phi*BZ
        JxB_Z      = J_phi*BR - J_R  *BP
        JxB_p      = J_R  *BZ - J_Z  *BR
        JxB_pol    = (JxB_R**2 + JxB_Z**2)**0.5
        GradP_R    = rho_x*TT + rho*TT_x
        GradP_Z    = rho_y*TT + rho*TT_y
        GradP_p    = rho_p*TT + rho*TT_p
        GradP_pol  = (GradP_R**2 + GradP_Z**2)**0.5
        
        scalars(inode,i_full(1)) = BR
        scalars(inode,i_full(2)) = BZ
        scalars(inode,i_full(3)) = BP
        scalars(inode,i_full(4)) = J_R
        scalars(inode,i_full(5)) = J_Z
        scalars(inode,i_full(6)) = J_phi
        scalars(inode,i_full(7)) = F_prof*dF_dpsi
        scalars(inode,i_full(8)) = F_prof 
        
        ! --- Choose which direction
        scalars(inode,i_full(9)) = GradP_pol ! GradP_p ! GradP_Z  ! GradP_R !
        scalars(inode,i_full(10))= JxB_pol   ! JxB_p   ! JxB_Z    ! JxB_R   !

        ! --- V_parallel
        scalars(inode,i_full(11))= (VR*BR + VZ*BZ + Vp*Bp)  / sqrt(BR**2 + BZ**2 + Bp**2)

        psi_norm = get_psi_n(A3, Z)

        grad_psi = sqrt(A3_R**2 + A3_Z**2)

        D_prof  = get_dperp (psi_norm)
        ZK_prof = get_zkperp(psi_norm)

        if ( eta_T_dependent ) then
          eta_T = eta * (max(TT,0.d0)/T_0)**(-1.5d0)
        else
          eta_T = eta
        end if

        if (include_bootstrap) then
          call bootstrap_current(R, Z, ES%R_axis, ES%Z_axis, ES%psi_axis, ES%R_xpoint, ES%Z_xpoint, ES%psi_bnd, psi_norm,&
                                 A3,     A3_R, A3_Z, rho, rho_x, rho_y,      &
                                 TT/2.0, TT_x/2.0, TT_y/2.0, TT/2.0, TT_x/2.0, TT_y/2.0, Jb   )
          scalars(inode,i_boot(1)) = Jb
          scalars(inode,i_boot(2)) = J_phi * R ! to compare against the RMHD-303 current that has a factor R in it
        else
          Jb = 0.d0
        endif

        if (include_fluxes) then

          scalars(inode,i_flux(1))   = scalars(inode,var_rho) * scalars(inode,var_T)

          if (grad_psi .ne. 0.d0) then

            scalars(inode,i_flux(2))  = ZKpar_T * ( BR * TT_x + BZ * TT_y + BP * TT_P / R) / sqrt(BR**2 + BZ**2 + Bp**2)

            scalars(inode,i_flux(3))  = ZK_prof * (TT_x * A3_R + TT_y * A3_Z) / grad_psi

            scalars(inode,i_flux(4))  = rho * TT * (VR * BR + VZ * BZ + VP * BP) / sqrt(BR**2 + BZ**2 + Bp**2)

            scalars(inode,i_flux(5))  = ( (VZ*Bp-Vp*BZ)**2 + (Vp*BR-VR*Bp)**2 + (VR*BZ-VZ*BR)**2 ) **0.5
            scalars(inode,i_flux(5))  = rho * TT * scalars(inode,i_flux(5)) / sqrt(BR**2 + BZ**2 + Bp**2)

            scalars(inode,i_flux(6))  = D_prof * (rho_x * A3_R + rho_y * A3_Z) / grad_psi

            scalars(inode,i_flux(7))  = rho * (VR * BR + VZ * BZ + VP * BP) / sqrt(BR**2 + BZ**2 + Bp**2)

            scalars(inode,i_flux(8))  = ( (VZ*Bp-Vp*BZ)**2 + (Vp*BR-VR*Bp)**2 + (VR*BZ-VZ*BR)**2 ) **0.5
            scalars(inode,i_flux(8))  = rho * scalars(inode,i_flux(8)) / sqrt(BR**2 + BZ**2 + Bp**2)

          endif ! grad_psi
        endif ! include_fluxes

        if (include_psi_norm) then
          scalars(inode,i_psin) = psi_norm
        endif

        if (include_magnetic_field) then
          vectors(inode,:,i_vec_B) = (/ BR, BZ, Bp /)
        endif

        if (include_velocity_field) then
          vectors(inode,:,i_vec_V) = (/ VR, VZ, Vp /)
        endif

        if (include_electric_field) then
          E_R   = - (VZ*BP - VP*BZ) + eta_T * J_R
          E_Z   = - (VP*BR - VR*BP) + eta_T * J_Z
          E_phi = - (VR*BZ - VZ*BR) + eta_T * J_phi
          vectors(inode,:,i_vec_E) =  (/ E_R, E_Z, E_phi /)
        endif

        if (include_Jpol) then
          vectors(inode,:, i_vec_Jpol) =  (/ J_R, J_Z, J_phi /)
        endif

#else /* not fullmhd */ 

        u_sum   = 0.d0; u_x  = 0.d0; u_y  = 0.d0; u_p  = 0.d0
        psi_sum = 0.d0; ps_x = 0.d0; ps_y = 0.d0; ps_p = 0.d0
        zj_sum  = 0.d0; zj_x = 0.d0; zj_y = 0.d0; zj_p = 0.d0
        T_sum   = 0.d0; TT_x = 0.d0; TT_y = 0.d0; TT_p = 0.d0
        zn_sum  = 0.d0; zn_x = 0.d0; zn_y = 0.d0; zn_p = 0.d0;
        Ti_sum  = 0.d0; Ti_x = 0.d0; Ti_y = 0.d0; Ti_p = 0.d0
        Te_sum  = 0.d0; Te_x = 0.d0; Te_y = 0.d0; Te_p = 0.d0
        w_sum   = 0.d0; w_x  = 0.d0; w_y  = 0.d0; w_p  = 0.d0; w_xx = 0.d0; w_yy = 0.d0
        E_R     = 0.d0; E_z  = 0.d0; E_phi= 0.d0
        dU_x    = 0.d0; dU_y = 0.d0
        V_sum = 0.d0;

        do i_tor = 1, n_tor

          if ( ( i_tor == 1 ) .and. ( without_n0_mode ) ) cycle ! Do not include the n=0 mode

          do m=1,n_var
             call interp(node_list,element_list,i,m,i_tor,s,t,P,P_s,P_t,P_st,P_ss,P_tt)
             scalars(inode,m) = scalars(inode,m) + P * HZ(i_tor,i_plane)
          enddo

          if (include_projections) then
            do i_proj=1,n_var
              call interp(aux_node_list,element_list,i,i_proj,i_tor,s,t,P,P_s,P_t,P_st,P_ss,P_tt)
              scalars(inode,iproj(i_proj)) = scalars(inode,iproj(i_proj)) + P * HZ(i_tor,i_plane)
            end do
          end if

          if (jorek_model .lt. 100) cycle
          
          call interp_delta(node_list,element_list,i,var_psi,i_tor,s,t,dpsi,dPs_s, dPs_t, dPs_st, dPs_ss, dPs_tt)
          call interp_delta(node_list,element_list,i,var_u,  i_tor,s,t,dU,dU_s, dU_t, dU_st, dU_ss, dU_tt)         

          call interp(node_list,element_list,i,var_psi,i_tor,s,t,Psi,Ps_s, Ps_t, Ps_st, Ps_ss, Ps_tt)
          call interp(node_list,element_list,i,var_u,  i_tor,s,t,U  ,U_s,  U_t,  U_st,  U_ss,  U_tt)
          call interp(node_list,element_list,i,var_zj, i_tor,s,t,ZJ ,ZJ_s, ZJ_t, ZJ_st, ZJ_ss, ZJ_tt)
          call interp(node_list,element_list,i,var_w,  i_tor,s,t,W  ,W_s,  W_t,  W_st,  W_ss,  W_tt)
          call interp(node_list,element_list,i,var_rho,i_tor,s,t,RHO,RHO_s,RHO_t,RHO_st,RHO_ss,RHO_tt)
          if (with_TiTe) then
             call interp(node_list,element_list,i,var_Ti,i_tor,s,t,Ti,Ti_s,Ti_t,Ti_st,Ti_ss,Ti_tt)
             call interp(node_list,element_list,i,var_Te,i_tor,s,t,Te,Te_s,Te_t,Te_st,Te_ss,Te_tt)
             TT   = Ti   + Te
             TT_s = Ti_s + Te_s
             TT_t = Ti_t + Te_t
          else
             call interp(node_list,element_list,i,var_T,  i_tor,s,t,TT ,TT_s, TT_t, TT_st, TT_ss, TT_tt)
          endif
         
          if (with_Vpar) then
             call interp(node_list,element_list,i,var_Vpar,i_tor,s,t,V,V_s,V_t,V_st,V_ss,V_tt)
          else
             V=0; V_s=0; V_t=0; V_st=0; V_ss=0; V_tt=0
          endif

          V_sum   = V_sum   + V   * HZ(i_tor,i_plane)
          psi_sum = psi_sum + psi * HZ(i_tor,i_plane)
          zj_sum  = zj_sum  + zj  * HZ(i_tor,i_plane)
          u_sum   = u_sum   + U   * HZ(i_tor,i_plane)
          w_sum   = w_sum   + w   * HZ(i_tor,i_plane)
          zn_sum  = zn_sum  + RHO * HZ(i_tor,i_plane)
          if (with_TiTe) then
            Ti_sum  = Ti_sum + Ti * HZ(i_tor,i_plane)
            Te_sum  = Te_sum + Te * HZ(i_tor,i_plane)
          else
            Ti_sum  = Ti_sum + 0.5d0*TT * HZ(i_tor,i_plane)
            Te_sum  = Te_sum + 0.5d0*TT * HZ(i_tor,i_plane)
          endif

          if ((xjac .gt. 1.d-6)) then  ! avoid the axis

             u_x  = u_x   + (   Z_t * U_s - Z_s * U_t )     / xjac * HZ(i_tor,i_plane)
             u_y  = u_y   + ( - R_t * U_s + R_s * U_t )     / xjac * HZ(i_tor,i_plane)

             du_x = du_x  + (   Z_t * dU_s - Z_s * dU_t )   / xjac * HZ(i_tor,i_plane)
             du_y = du_y  + ( - R_t * dU_s + R_s * dU_t )   / xjac * HZ(i_tor,i_plane)

             ps_x  = ps_x + (   Z_t * PS_s - Z_s * PS_t )   / xjac * HZ(i_tor,i_plane)
             ps_y  = ps_y + ( - R_t * PS_s + R_s * PS_t )   / xjac * HZ(i_tor,i_plane)

             zj_x  = zj_x + (   Z_t * ZJ_s - Z_s * ZJ_t )   / xjac * HZ(i_tor,i_plane)
             zj_y  = zj_y + ( - R_t * ZJ_s + R_s * ZJ_t )   / xjac * HZ(i_tor,i_plane)

             TT_x  = TT_x + (   Z_t * TT_s - Z_s * TT_t )   / xjac * HZ(i_tor,i_plane)
             TT_y  = TT_y + ( - R_t * TT_s + R_s * TT_t )   / xjac * HZ(i_tor,i_plane)
             TT_p  = TT_p + TT * HZ_p(i_tor,i_plane)

             if (with_TiTe) then
               Ti_x  = Ti_x + (   Z_t * Ti_s - Z_s * Ti_t )   / xjac * HZ(i_tor,i_plane)
               Ti_y  = Ti_y + ( - R_t * Ti_s + R_s * Ti_t )   / xjac * HZ(i_tor,i_plane)
               Ti_p  = Ti_p + Ti * HZ_p(i_tor,i_plane)
               Te_x  = Te_x + (   Z_t * Te_s - Z_s * Te_t )   / xjac * HZ(i_tor,i_plane)
               Te_y  = Te_y + ( - R_t * Te_s + R_s * Te_t )   / xjac * HZ(i_tor,i_plane)
               Te_p  = Te_p + Te * HZ_p(i_tor,i_plane)
             else
               Ti_x  = TT_x / 2.0
               Ti_y  = TT_y / 2.0
               Ti_p  = TT_p / 2.0
               Te_x  = TT_x / 2.0
               Te_y  = TT_y / 2.0
               Te_p  = TT_p / 2.0
             end if

             zn_x = zn_x  + (   Z_t * RHO_s - Z_s * RHO_t ) / xjac * HZ(i_tor,i_plane)
             zn_y = zn_y  + ( - R_t * RHO_s + R_s * RHO_t ) / xjac * HZ(i_tor,i_plane)
             zn_p = zn_p  + RHO * HZ_p(i_tor,i_plane)

             w_x  = w_x   + (   Z_t * w_s - Z_s * w_t )     / xjac * HZ(i_tor,i_plane)
             w_y  = w_y   + ( - R_t * w_s + R_s * w_t )     / xjac * HZ(i_tor,i_plane)

             w_xx = w_xx  + (w_ss * Z_t**2 - 2.d0*w_st * Z_s*Z_t + w_tt * Z_s**2        &
                  + w_s * (Z_st*Z_t - Z_tt*Z_s )                                        &
                  + w_t * (Z_st*Z_s - Z_ss*Z_t ) )     / xjac**2                        &
                  - xjac_x * (w_s* Z_t - w_t * Z_s)  / xjac**2

             w_yy = w_yy  + (w_ss * R_t**2 - 2.d0*w_st * R_s*R_t + w_tt * R_s**2        &
                  + w_s * (R_st*R_t - R_tt*R_s )                                        &
                  + w_t * (R_st*R_s - R_ss*R_t ) )       / xjac**2                      &
                  - xjac_y * (- w_s * R_t + w_t * R_s )  / xjac**2

            ! --- Full toroidal electric field evaluated at t_now - dt/2
            E_R   = E_R   - F0*(U_x-0.5d0*dU_x)*HZ(i_tor,i_plane)
            E_Z   = E_Z   - F0*(U_y-0.5d0*dU_y)*HZ(i_tor,i_plane) 
            E_phi = E_phi - dpsi/tstep * HZ(i_tor,i_plane)/BigR - F0*(U-0.5d0*dU)*HZ_p(i_tor,i_plane)/BigR 

          endif ! xjac

        enddo  ! end loop toroidal harmonics
        if (jorek_model .lt. 100) cycle

        Psi_tot = 0.d0
        do i_tor =1, n_tor
           call interp(node_list,element_list,i,1,i_tor,s,t,P,P_s,P_t,P_st,P_ss,P_tt)
           Psi_tot = Psi_tot + P * HZ(i_tor,i_plane)
        enddo
 
        psi_norm = get_psi_n(Psi_tot, Z)

        if (include_bootstrap) then
          call bootstrap_current(R, Z, ES%R_axis, ES%Z_axis, ES%psi_axis, ES%R_xpoint, ES%Z_xpoint, ES%psi_bnd, psi_norm,&
                                 psi_sum, ps_x, ps_y, zn_sum,  zn_x, zn_y,      &
                                 Ti_sum,  Ti_x, Ti_y, Te_sum,  Te_x, Te_y, Jb   )
          scalars(inode,i_boot(1)) = Jb
          scalars(inode,i_boot(2)) = bootstrap_spline3_eval(n_spline_vtk-1,psi_knots_vtk,j_knots_vtk,j_spline_vtk,psi_norm)
          ! --- The JOREK bootstrap is not constant on flux surface
          ! --- Because J_jorek = R*J_physical, and the physical bootstrap is constant on a surface
          ! --- Hence, it makes more sense to look at R*j_average to compare with the bootstrap...
          scalars(inode,i_boot(2)) = R* scalars(inode, i_boot(2)) / ES%R_axis
        else
          Jb = 0.d0
        endif

        v_perp  = R * sqrt(u_x*u_x + u_y * u_y)
        Btot    = sqrt(F0**2 + ps_x**2 + ps_y**2) / BigR
        D_prof  = get_dperp (psi_norm)
        ZK_prof = get_zkperp(psi_norm)

        ZKpar_T = ZK_par * ((max( scalars(inode,6), T_min ))/T_0)**2.5

        grad_psi = sqrt(ps_x*ps_x + ps_y*ps_y)

        if (SI_units .and. with_Vpar) then
          scalars(inode,var_Vpar) = scalars(inode,var_Vpar) * sign(Btot,F0)   ! with si-units= .f. gives jorek variable, otherwise physical v_par
        endif

        !   'E_flux_Kpar ','E_flux_kperp','E_flux_Vpar ','E_flux_Vperp','D_flux_Dperp','D_flux_Vpar ','D_flux_Vperp'/)

        if (include_fluxes) then

          scalars(inode,i_flux(1))  = scalars(inode,5) * scalars(inode,6)

          if (grad_psi .ne. 0.d0) then

            scalars(inode,i_flux(2))  = ZKpar_T * ( F0 * TT_p / BigR**2  + (TT_x * ps_y - TT_y * ps_x) / BigR ) / Btot

            scalars(inode,i_flux(3))  = ZK_prof * (TT_x * ps_x + TT_y * ps_y) / grad_psi

            scalars(inode,i_flux(4))  = scalars(inode,5) * scalars(inode,6) * scalars(inode,7) * Btot

            scalars(inode,i_flux(5))  = BigR   * (u_x * ps_y - u_y * ps_x) / sqrt(ps_x*ps_x + ps_y*ps_y) * scalars(inode,5) * scalars(inode,6)

            scalars(inode,i_flux(6))  = D_prof * (zn_x * ps_x + zn_y * ps_y) / sqrt(ps_x*ps_x + ps_y*ps_y)

            scalars(inode,i_flux(7))  = scalars(inode,5) * scalars(inode,7) * Btot

            scalars(inode,i_flux(8))  = BigR   * (u_x * ps_y - u_y * ps_x) / sqrt(ps_x*ps_x + ps_y*ps_y) * scalars(inode,5)

          endif ! grad_psi
        endif ! include_fluxes

        if (include_magnetic_field) then
          vectors(inode,:,i_vec_B) = (/ ps_y/BigR, -ps_x/BigR, F0/BigR /)          
        endif

        if (include_velocity_field) then
          vectors(inode,:,i_vec_V) = (/ -BigR*u_y + V_sum/BigR*ps_y, BigR*u_x - V_sum/BigR*ps_x, V_sum*F0/BigR /)          
        endif

        if (include_electric_field) then
          vectors(inode,:,i_vec_E) =  (/ E_R, E_Z, E_phi /)
        endif

        if (include_Jpol) then
          call J_pol(node_list, element_list, i, s, t, i_plane,.false., Jpol_R, Jpol_Z, FFp)
          vectors(inode,:, i_vec_Jpol) =  (/ Jpol_R, Jpol_Z, 0.d0 /)
        endif
        
        if (include_psi_norm) then
           scalars(inode,i_psin) = psi_norm
        endif

        if (use_pellet) then

           local_density     = scalars(inode,5)
           local_temperature = scalars(inode,6)/2.
           local_psi         = scalars(inode,1)
           local_pressure    = local_density * local_temperature
           scalars(inode,i_pellet(1))  = local_pressure

           angle = 0.0
           call pellet_source2(pellet_amplitude,pellet_R,pellet_Z,pellet_psi,pellet_phi, &
                pellet_radius, pellet_delta_psi, pellet_sig, pellet_length, pellet_ellipse, &
                pellet_theta, R, Z, local_psi, angle, &
                local_density,local_temperature, &
                central_density, pellet_particles, pellet_density, &
                total_pellet_volume, local_source, source_volume)

           scalars(inode,i_pellet(2)) = local_source
        endif ! use_pellet

        ! vectors(inode,:,1) = (/ - R * u0_y ,   + R * u0_x ,   0.d0 /)
        ! vectors(inode,:,2) = (/ + ps_y /R * scalars(inode,7), - ps_x /R * scalars(inode,7), 0.d0 /) * Btot
        ! vectors(inode,:,3) = (/ - R * u0_y + ps_y /R * scalars(inode,7) * Btot, + R * u0_x - ps_x /R * scalars(inode,7) * Btot, 0.d0 /)

#endif /* end of non-full-MHD part */

        currdens(inode) = -scalars(inode,3)/BigR

      endif ! i_tor from 1 to n_tor

    enddo  ! nsub
  enddo     ! nsub

  do j=1,nsub-1
    do k=1,nsub-1
      ielm	  = ielm+1
      ien(1,ielm) = inode - nsub*nsub + nsub*(j-1) + k-1       ! 0 based indices for VTK
      ien(2,ielm) = inode - nsub*nsub + nsub*(j  ) + k-1
      ien(3,ielm) = inode - nsub*nsub + nsub*(j  ) + k
      ien(4,ielm) = inode - nsub*nsub + nsub*(j-1) + k
    enddo
  enddo

enddo  ! n_elements

#if (!defined WITH_Impurities)
  if (deuterium_adas)  ad_deuterium =  read_adf11(0,'96_h') !< for both include_radiation and include_neutral_dens
  if (include_radiation) then
    do i=1,nnos
      r0_real8  = scalars(i,var_rho)
      if ( with_TiTe ) then
        T_real8 = scalars(i,var_Te)
        Te_corr_eV = corr_neg_temp(T_real8*2.d0)/(2.d0*EL_CHG*MU_ZERO*central_density*1.d20)
        Te_eV = T_real8/(EL_CHG*MU_ZERO*central_density*1.d20)
      else
        T_real8 = scalars(i,var_T)
        Te_corr_eV = corr_neg_temp(T_real8)/(2.d0*EL_CHG*MU_ZERO*central_density*1.d20)
        Te_eV = T_real8/(2.d0*EL_CHG*MU_ZERO*central_density*1.d20)
      endif

#if (defined WITH_Neutrals) 
      coef_ion_3 = 27.2d0*EL_CHG*MU_ZERO*central_density*1.d20
      coef_ion_2 = 0.232d0
      coef_ion_1 = (MU_ZERO*central_mass*MASS_PROTON)**(0.5d0)*0.2917d-13*(central_density*1.d20)**(1.5d0)
      S_ion_puiss = 3.9d-1

      ksiion = ksi_ion * central_density * 1.d20
      rn0_real8 = scalars(i,var_rhon)

      if ( with_TiTe ) then
        call atomic_coeff_deuterium(T_real8, Sion_T, dSion_dT, Srec_T, dSrec_dT,        &
                                  LradDcont_T, dLradDcont_dT, LradDrays_T, dLradDrays_dT,r0_real8,rn0_real8,.true. ) !< add scalars(i,var_rho) as last optional parameter for density dependence
      else
        call atomic_coeff_deuterium(0.5d0*T_real8, Sion_T, dSion_dT, Srec_T, dSrec_dT,        &
                                  LradDcont_T, dLradDcont_dT, LradDrays_T, dLradDrays_dT,r0_real8,rn0_real8,.true. ) 
      endif


      eta_Sp = 1.65d-9*17*(1.d-3*Te_corr_eV)**(-1.5d0) &
                              *(central_mass*MASS_PROTON*central_density * 1.d20/MU_ZERO)**(0.5d0)

      scalars(i,ineu(1)) = ksiion * scalars(i,var_rho) * scalars(i,var_rhon) * Sion_T
      scalars(i,ineu(2)) = scalars(i,var_rho) * scalars(i,var_rhon) * LradDrays_T
      scalars(i,ineu(3)) = LradDcont_T * scalars(i,var_rho)**2.d0
#ifdef fullmhd
      scalars(i,ineu(4)) = 0.d0   ! NEEDS BE CALCULATED FOR FULL MHD ELESEWHERE! 
#else /* not fullmhd */ 
      scalars(i,ineu(4)) = (2/(3 * BigR**2)) * eta_Sp * scalars(i,var_zj)**2.d0
#endif /* with neutrals */
#endif /* with neutrals */
      
      !--------------------------------------------------------
      ! --- Radiation from background impurity
      !--------------------------------------------------------


      r0_corr   = corr_neg_dens(r0_real8)
      ne_SI = r0_corr * 1.d20 * central_density !electron density (SI)
   
      if (use_imp_adas) then  ! use open adas by default
        frad_bg = 0. 
        do i_imp =1, n_adas
          r_imp = nimp_bg(i_imp) / (1.d20 * central_density)  ! Background impurity density in JU 
          if (ne_SI > ne_SI_min .and. Te_eV > Te_eV_min .and. r_imp > 0) then
            Lrad_imp = 0.0
            call radiation_function_linear(imp_adas(i_imp),imp_cor(i_imp),log10(ne_SI),    & 
                                           log10(Te_corr_eV*EL_CHG/K_BOLTZ),.true.,Lrad_imp)    
          else     
            Lrad_imp = 0.
          end if
          frad_bg = frad_bg + r_imp * Lrad_imp
          scalars(i,iibg(i_imp)) = r_imp * Lrad_imp * scalars(i,var_rho)
        end do
      else
        if ( trim(imp_type(1)) == 'Ar') then ! Hard-coded fitting exists for argon
          Arad_bg = 2.4d-31
          Brad_bg = 20.
          Crad_bg = 0.8
          frad_bg = (2./3.)*(1./(central_mass*MASS_PROTON))                               &
                     *((MU_ZERO*central_mass*MASS_PROTON*central_density*1.d20)**(1.5d0)) &
                     *nimp_bg(1)* Arad_bg*exp(-((log(Te_corr_eV)-log(Brad_bg))**2.)/Crad_bg**2.)
        else
          write(*,*) "WARNING: hard-coded fitting doesn't exist for  ", trim(imp_type(1)), ", use open adas instead!"
          stop
        end if
      end if   


      scalars(i,ibg_tot) = scalars(i,var_rho) * frad_bg
   
    enddo
  endif
#endif /* (.not. with_impurities) */

#ifdef WITH_Impurities

 if (include_radiation) then

  !-------------------------------------------
  ! Atomic physics parameters for Impurities
  !-------------------------------------------

   select case ( trim(imp_type(index_main_imp)) )
     case('D2')
       m_i_over_m_imp = central_mass/2.  ! Deuterium mass = 2 u
     case('Ar')
       m_i_over_m_imp = central_mass/40. ! Argon mass = 40 u
     case('Ne')
       m_i_over_m_imp = central_mass/20. ! Neon mass = 20 u
     case default
       write(*,*) '!! Gas type "', trim(imp_type(index_main_imp)), '" unknown (in mod_injection_source.f90) !!'
       write(*,*) '=> We assume the gas is D2.'
       m_i_over_m_imp = central_mass/2.
   end select

   do i=1,nnos
     if ( with_TiTe ) then
       T_real8 = scalars(i,var_Te)
       Te_corr_eV = corr_neg_temp(T_real8*2.d0)/(2.d0*EL_CHG*MU_ZERO*central_density*1.d20)
       Te_eV = T_real8/(EL_CHG*MU_ZERO*central_density*1.d20)
     else
       T_real8 = scalars(i,var_T)
       Te_corr_eV = corr_neg_temp(T_real8)/(2.d0*EL_CHG*MU_ZERO*central_density*1.d20)
       Te_eV = T_real8/(2.d0*EL_CHG*MU_ZERO*central_density*1.d20)
     endif
     
     eta_Sp = 1.65d-9*17*(1.d-3*Te_corr_eV)**(-1.5d0) &
                        *(central_mass*MASS_PROTON*central_density * 1.d20/MU_ZERO)**(0.5d0)

     r0_real8 = scalars(i,var_rho)
     rn0_real8 = scalars(i,var_rhon)

     r0_corr = corr_neg_dens(r0_real8,(/1.d-9,1.d-5/),1.d-3)
     rn0_corr = corr_neg_dens(rn0_real8,(/1.d-9,1.d-5/),1.d-3)

     ! We estimate the effective charge by a test density 10^20/m^3
     ! Later maybe we should implement a iterative method

     if (allocated(imp_adas(index_main_imp)%ionisation_energy) .and. (.not. use_marker)) then
       if (allocated(P_imp)) deallocate(P_imp)
       allocate(P_imp(0:imp_adas(index_main_imp)%n_Z))

       call imp_cor(index_main_imp)%interp_linear(density=20.,temperature=log10(Te_corr_eV*EL_CHG/K_BOLTZ),&
                                     p_out=P_imp,z_avg=Z_imp)

       ! Calculate the ionization potential energy and derivative wrt. temperature
       E_ion     = 0.

       do ion_i=1, imp_adas(index_main_imp)%n_Z
         do ion_k=1, ion_i
           E_ion     = E_ion + P_imp(ion_i)*imp_adas(index_main_imp)%ionisation_energy(ion_k)
         end do
       end do
     ! Convert from eV to JOREK unit
       E_ion     = E_ion * EL_CHG*MU_ZERO*central_density*1.d20*m_i_over_m_imp
     else if (.not. use_marker) then
       call imp_cor(index_main_imp)%interp_linear(density=20.,temperature=log10(Te_corr_eV*EL_CHG/K_BOLTZ),z_avg=Z_imp)
       E_ion     = 0.
       write(*,*) 'The ionization energy file of the main impuritity is required'
       stop
     else ! use_marker
       E_ion     = 0.
       Z_imp     = 0.
       Z_eff     = 0.

       Z_eff     = max(scalars(i,iproj(3)),0.0)
       Z_imp     = max(scalars(i,iproj(4)),0.0)
       n_imp     = max(scalars(i,iproj(5)),1.d-2)
       Z_eff     = Z_eff / n_imp
       Z_imp     = Z_imp / n_imp
     end if

     !alpha_imp    = 0.5*m_i_over_m_imp*(Z_imp+1.) - 1.
     beta_imp     = m_i_over_m_imp*Z_imp - 1.

     ne_SI       = (r0_corr + beta_imp * rn0_corr) * 1.d20 * central_density ! electron density (SI)
     scalars(i,var_rho) = (r0_corr + beta_imp * rn0_corr)                           ! electron density (JOREK units)

     !Calculate the Z_eff, as it is done in mod_elt_matrix
     if (.not. use_marker) then
       Z_eff = r0_corr - rn0_corr
       do ion_i=1, imp_adas(index_main_imp)%n_Z
         Z_eff = Z_eff + m_i_over_m_imp * rn0_corr * P_imp(ion_i) * real(ion_i,8)**2
       end do
       Z_eff = Z_eff / scalars(i,var_rho)
     else
       Z_eff = Z_eff * n_imp
       Z_eff = Z_eff + max((r0_corr - rn0_corr),0.)
       Z_eff = max(Z_eff / scalars(i,var_rho), 1.)
     end if

     scalars(inode,iimp(5)) = Z_eff

     ! This is to represent the dependence on Z_eff in resistivity
     eta_coef = Z_eff*(1.+1.198*Z_eff+0.222*Z_eff**2)/(1.+2.966*Z_eff+0.753*Z_eff**2) 
     eta_coef = eta_coef / ((1.+1.198+0.222)/(1.+2.966+0.753))

     eta_Sp = eta_Sp * eta_coef

  !-------------------------------------------
  ! --- Radiative function, using interpolation
  ! ------------------------------------------
     if (ne_SI > ne_SI_min .and. Te_eV > Te_eV_min .and. rn0_real8 > rn0_min) then

       Lrad = 0.0
       
       ! Here we are temperarily only considering one impurity species, in the
       ! future maybe a do loop will be needed
       call radiation_function_linear(imp_adas(index_main_imp),imp_cor(index_main_imp),log10(ne_SI),log10(Te_corr_eV*EL_CHG/K_BOLTZ),.true.,Lrad)
       Lrad = Lrad * m_i_over_m_imp

     else
       Lrad = 0.
       E_ion = 0.
     end if

     frad_bg = 0.
     do i_imp =1, n_adas
       r_imp = nimp_bg(i_imp) / (1.d20 * central_density)  ! Background impurity density in JU     
       if (ne_SI > ne_SI_min .and. Te_eV > Te_eV_min .and. r_imp > 0) then
         Lrad_imp = 0.0
         call radiation_function_linear(imp_adas(i_imp),imp_cor(i_imp),log10(ne_SI),    & 
                                        log10(Te_eV*EL_CHG/K_BOLTZ),.true.,Lrad_imp)           
       else     
         Lrad_imp = 0.
       end if
       frad_bg = frad_bg + r_imp * Lrad_imp
       scalars(i,iibg(i_imp)) = scalars(i,var_rho) * r_imp * Lrad_imp
     end do
     scalars(i,iimp(1)) = (2./3.) * scalars(i,var_rhon) * E_ion
     scalars(i,iimp(2)) = (r0_corr+beta_imp*rn0_corr) * rn0_corr * Lrad
     scalars(i,iimp(3)) = (2./(3. * BigR**2)) * eta_Sp * scalars(i,var_zj)**2.d0
     scalars(i,iimp(4)) = Z_imp
     scalars(i,iimp(5)) = Z_eff
     scalars(i,ibg_tot) = scalars(i,var_rho) * frad_bg

   end do
 endif
#endif /*WITH_Impurities*/

#if (defined WITH_Neutrals) && (!defined WITH_Impurities)
  if (include_neutral_dens) then

    do i=1,nnos

      r0_real8  = scalars(i,var_rho)
      rn0_real8 = scalars(i,var_rhon) 

      if ( with_TiTe ) then
        T_real8 = scalars(i,var_Te)
        call atomic_coeff_deuterium(T_real8, Sion_T, dSion_dT, Srec_T, dSrec_dT,        &
                                  LradDcont_T, dLradDcont_dT, LradDrays_T, dLradDrays_dT, r0_real8, rn0_real8, .true. )
      else
        T_real8 = scalars(i,var_T)
        call atomic_coeff_deuterium(0.5d0*T_real8, Sion_T, dSion_dT, Srec_T, dSrec_dT,        &
                                  LradDcont_T, dLradDcont_dT, LradDrays_T, dLradDrays_dT, r0_real8, rn0_real8, .true. )
      endif

      r0_corr   = corr_neg_dens(r0_real8)
      rn0_corr  = corr_neg_dens(rn0_real8, (/ 0.d-5, 1.d-5 /))

      IonN      = -(r0_corr) * (rn0_corr) * Sion_T
      RecN      = (r0_corr)**2 * Srec_T 

      scalars(i,ineu(6)) = IonN
      scalars(i,ineu(7)) = RecN

    end do
  end if
#endif /* WITH_Neutrals but not WITH_Impurities */


if (SI_units) then

  !===========================================================real values=============
  rho_norm = central_density*1.d20 * central_mass * mass_proton
  t_norm   = sqrt(MU_zero*rho_norm)

  !=================================================real values============
  do i=1,nnos

#ifdef fullmhd
    !===========================================density in 1e20m-3
    scalars(i,var_rho) = scalars(i,var_rho) * central_density
    !===========================================electron temperature in keV
    if (with_TiTe) then
      scalars(i,var_Ti ) = scalars(i,var_Ti)  / MU_zero / (central_density * 1d20) / EL_CHG /1.e3 
      scalars(i,var_Te ) = scalars(i,var_Te)  / MU_zero / (central_density * 1d20) / EL_CHG /1.e3 
    else
      scalars(i,var_T  ) = scalars(i,var_T)   / MU_zero / (central_density * 1d20) / EL_CHG /2./1.e3 !(assumes Te=Ti=T/2)
    endif
    !===========================================Velocity
    scalars(i,var_UR) = scalars(i,var_UR) /t_norm/1.e3
    scalars(i,var_UZ) = scalars(i,var_UZ) /t_norm/1.e3
    scalars(i,var_Up) = scalars(i,var_Up) /t_norm/1.e3
    scalars(i,i_full(11)) = scalars(i,i_full(11)) /t_norm/1.e3 ! V_parallel
    !=====================Pressure in kPa
    if (include_fluxes) scalars(i,i_flux(1)) = scalars(i,i_flux(1)) / MU_zero/1.e3
    ! not yet implemented
    !if (include_neo) then
    !  !============================Er in kV/m
    !  scalars(i,s_neo+ 1) = F0*scalars(i,s_neo+ 1) / t_norm/1.e3
    !  !====================================Vtheta km/s
    !  scalars(i,s_neo+ 2) =    scalars(i,s_neo+ 2) / t_norm/1.e3
    !  !===================================Vsound in km/s
    !  scalars(i,s_neo+ 5) =    scalars(i,s_neo+ 5) / t_norm/1.e3
    !  !===================================Vneo in km/s
    !  scalars(i,s_neo+ 7) =    scalars(i,s_neo+ 7) / t_norm/1.e3
    !  !===================================Vperp_e in km/s
    !  scalars(i,s_neo+ 8) =    scalars(i,s_neo+ 8) / t_norm/1.e3
    !  ! ===================================mu_neo in SI units
    !  scalars(i,s_neo+10) =    scalars(i,s_neo+10) / sqrt(rho_norm*MU_zero)
    !endif
    !============================================j_bootstrap, javeraged in MA/m2
    if (include_bootstrap) then
    scalars(i,i_boot(1))=scalars(i,i_boot(1))/MU_zero*1.e-6
    scalars(i,i_boot(2))=scalars(i,i_boot(2))/MU_zero*1.e-6
    endif
    if (include_velocity_field) then 
      vectors(i,:,i_vec_V) = vectors(i,:,i_vec_V)/t_norm
    endif
    if (include_electric_field) then 
      vectors(i,:,i_vec_E) = vectors(i,:,i_vec_E)/t_norm
    endif
    if (include_Jpol) then
      vectors(i,:, i_vec_Jpol) = vectors(i,:,i_vec_Jpol)/MU_zero*1e-6
    endif

#else /* not full-MHD */

    !============================================u in m/s
    scalars(i,var_u) = scalars(i,var_u)/t_norm
    !============================================j_phi in MA/m2
    scalars(i,var_zj) = currdens(i) / MU_zero * 1.e-6
    !============================================density in 1e20m-3
    scalars(i,var_rho) = scalars(i,var_rho) * central_density
    if (with_TiTe) then
      !===========================================ion and electron temperatures in keV
      scalars(i,var_Ti) = scalars(i,var_Ti) / MU_zero / (central_density * 1d20) / EL_CHG /1.e3 !
      scalars(i,var_Te) = scalars(i,var_Te) / MU_zero / (central_density * 1d20) / EL_CHG /1.e3 !
    else
    !===========================================electron temperature in keV
      scalars(i,var_T) = scalars(i,var_T) / MU_zero / (central_density * 1d20) / EL_CHG /2./1.e3 !(assumes Te=Ti=T/2)
    endif
    !=====================================Vparal in km/s *Btot!!!
    if (with_Vpar) then
       scalars(i,var_Vpar) = scalars(i,var_Vpar) /t_norm/1.e3
    endif
#if (defined WITH_Neutrals) && (!defined WITH_Impurities)
    !===================================== Neutral density in 1e20m-3
    scalars(i,var_rhon) = scalars(i,var_rhon) * central_density
#endif

#ifdef WITH_Impurities
    !===================================== Impurity density in 1e20m-3
    scalars(i,var_rhon) = scalars(i,var_rhon) * central_density * m_i_over_m_imp
#endif
    !=====================Pressure in kPa
    if (include_fluxes) scalars(i,i_flux(1)) = scalars(i,i_flux(1)) / MU_zero/1.e3
    if (include_neo) then
      !============================Er in kV/m
      scalars(i,i_neo( 1)) = F0*scalars(i,i_neo( 1)) / t_norm/1.e3
      !====================================Vtheta km/s
      scalars(i,i_neo( 2)) =    scalars(i,i_neo( 2)) / t_norm/1.e3
      !===================================Vsound in km/s
      scalars(i,i_neo( 5)) =    scalars(i,i_neo( 5)) / t_norm/1.e3
      !===================================Vneo in km/s
      scalars(i,i_neo( 7)) =    scalars(i,i_neo( 7)) / t_norm/1.e3
      !===================================Vperp_e in km/s
      scalars(i,i_neo( 8)) =    scalars(i,i_neo( 8)) / t_norm/1.e3
      ! ===================================mu_neo in SI units
      scalars(i,i_neo(10)) =    scalars(i,i_neo(10)) / sqrt(rho_norm*MU_zero)
    endif
    !============================================j_bootstrap, javeraged in MA/m2

    if (include_bootstrap) then
    scalars(i,i_boot(1))=scalars(i,i_boot(1))/MU_zero*1.e-6
    scalars(i,i_boot(2))=scalars(i,i_boot(2))/MU_zero*1.e-6
    endif
    if (include_velocity_field) then 
      vectors(i,:,i_vec_V) = vectors(i,:,i_vec_V)/t_norm
    endif
    if (include_electric_field) then 
      vectors(i,:,i_vec_E) = vectors(i,:,i_vec_E)/t_norm
    endif
    if (include_Jpol) then
      vectors(i,:, i_vec_Jpol) = vectors(i,:,i_vec_Jpol)/MU_zero*1e-6
    endif
 
    !========================================================

#if (!defined WITH_Impurities)
    if (include_radiation) then
      r0_real8  = scalars(i,var_rho)/central_density ! Back to JOREK unit for calling atomic_coeff_deuterium
      if ( with_TiTe ) then
        T_real8 = scalars(i,var_Te)*1.e3*EL_CHG*MU_zero*(central_density * 1.d20) ! T_real8 back to JOREK units
        Te_corr_eV = corr_neg_temp(T_real8*2.d0)/(2.d0*EL_CHG*MU_ZERO*central_density*1.d20)
        Te_eV = T_real8/(EL_CHG*MU_ZERO*central_density*1.d20)
      else
        T_real8 = scalars(i,var_T)*1.e3*2.*EL_CHG*MU_zero*(central_density * 1.d20)
        Te_corr_eV = corr_neg_temp(T_real8)/(2.d0*EL_CHG*MU_ZERO*central_density*1.d20)
        Te_eV = T_real8/(2.d0*EL_CHG*MU_ZERO*central_density*1.d20)
      endif

#if (defined WITH_Neutrals) 
      coef_ion_1 = (MU_ZERO*central_mass*MASS_PROTON)**(0.5d0)*(central_density*1.d20)**(1.5d0)
      coef_rad_1 = (gamma-1.d0)*MU_ZERO**1.5d0*(central_mass*MASS_PROTON)**0.5d0*(central_density*1.d20)**2.5d0

      ksiion = ksi_ion * central_density * 1.d20
      rn0_real8 = scalars(i,8)/central_density

      if ( with_TiTe ) then
        call atomic_coeff_deuterium(T_real8, Sion_T, dSion_dT, Srec_T, dSrec_dT,        &
                                  LradDcont_T, dLradDcont_dT, LradDrays_T, dLradDrays_dT, r0_real8,rn0_real8,.true. )  ! T, rho and rhon should be in JOREK unit here
      else
        call atomic_coeff_deuterium(0.5d0*T_real8, Sion_T, dSion_dT, Srec_T, dSrec_dT,        &
                                  LradDcont_T, dLradDcont_dT, LradDrays_T, dLradDrays_dT, r0_real8,rn0_real8,.true. ) 
      endif

      eta_Sp = 1.65d-9*17*(1.d-3*Te_corr_eV)**(-1.5d0)
  
      scalars(i,ineu(1)) = ksiion* (1.5d0)/(MU_zero*central_density*1.d20)      &
                                          * scalars(i,var_rho) * 1.d20 * scalars(i,var_rhon) * 1.d20 * Sion_T / coef_ion_1

      scalars(i,ineu(2)) = scalars(i,var_rho)* 1.d20 * scalars(i,var_rhon) * 1.d20 * LradDrays_T/ coef_rad_1

      scalars(i,ineu(3)) = LradDcont_T * (scalars(i,var_rho)*1.d20)**2.d0 / coef_rad_1

      scalars(i,ineu(4)) = eta_Sp * (1.d6*scalars(i,var_zj))**2.d0
#endif /* WITH_Neutrals but not WITH_Impurities */
      !--------------------------------------------------------
      ! --- Radiation from background impurity
      !--------------------------------------------------------
      r0_corr   = corr_neg_dens(r0_real8)
      ne_SI = r0_corr * 1.d20 * central_density !electron density (SI)

      if (use_imp_adas) then  ! use open adas by default
        ! Use radiation coefficients from ADAS
        frad_bg = 0. 
        do i_imp =1, n_adas  
          if (ne_SI > ne_SI_min .and. Te_eV > Te_eV_min .and. nimp_bg(i_imp) > 0) then
            Lrad_imp = 0.0
            call radiation_function_linear(imp_adas(i_imp),imp_cor(i_imp),log10(ne_SI),    & 
                                           log10(Te_corr_eV*EL_CHG/K_BOLTZ),.false.,Lrad_imp) 
          else     
            Lrad_imp = 0.
          end if
          scalars(i,iibg(i_imp)) = 1.d20 * scalars(i,var_rho) * nimp_bg(i_imp) * Lrad_imp 
          frad_bg = frad_bg + nimp_bg(i_imp) * Lrad_imp
        end do
      else
        if ( trim(imp_type(1)) == 'Ar' ) then ! Hard-coded fitting exists for argon
          Arad_bg = 2.4d-31
          Brad_bg = 20.
          Crad_bg = 0.8
          frad_bg = nimp_bg(1) * Arad_bg * exp(-((log(Te_corr_eV)-log(Brad_bg))**2.)/Crad_bg**2.)
        else 
          write(*,*) "WARNING: hard-coded fitting doesn't exist for  ", trim(imp_type(1)), ", use open adas instead!"
          stop
        end if
      end if
      scalars(i,ibg_tot) = scalars(i,var_rho)*1.d20 * frad_bg
    endif
#endif /*(.not. with_Impurities)*/

#ifdef WITH_Impurities
  if (include_radiation) then
   scalars(i,iimp(1)) = scalars(i,iimp(1))/(K_BOLTZ*MU_ZERO)
   scalars(i,iimp(2)) = scalars(i,iimp(2))/(2.d0/3.d0*MU_ZERO**1.5d0*(central_mass*MASS_PROTON*central_density*1.d20)**0.5d0)
   scalars(i,iimp(3)) = scalars(i,iimp(3))/(2.d0/3.d0*((central_mass*MASS_PROTON*central_density*1.d20)**0.5)*(MU_ZERO**1.5)) 
   scalars(i,iimp(4)) = scalars(i,iimp(4))
   scalars(i,ibg_tot) = scalars(i,ibg_tot) &
       /((GAMMA-1.d0)*MU_ZERO**1.5d0*(central_mass*MASS_PROTON*central_density*1.d20)**0.5d0)
   do i_imp=1,n_adas
     scalars(i,iibg(i_imp)) = scalars(i,iibg(i_imp)) / ((GAMMA-1.d0)*MU_ZERO**1.5d0*(central_mass*MASS_PROTON*central_density*1.d20)**0.5d0)
   end do
  end if
#endif /* WITH_Impurities */
#endif /* end of non-full-MHD part*/

  enddo  ! nnos

endif ! SI_UNITS

!--------------------------------------------------- write the binary VTK file
etype = 9  ! for vtk_quad

call write_vtk('jorek_tmp.vtk',xyz,ien,etype,scalar_names,scalars,vector_names,vectors)

write(*,*) 'done.'



contains



  subroutine add_vtk_entry(name, name_si, index, counter, si_units, names_list)
  
    implicit none
  
    !--- Routine parameters
    integer,                   intent(inout) :: index
    integer,                   intent(inout) :: counter
    character*36, allocatable, intent(inout) :: names_list(:)
    character*12,              intent(in)    :: name, name_si           
    logical,                   intent(in)    :: si_units
    
    !--- Local parameters  
    character*36, allocatable :: names_tmp(:)
    character*36              :: final_name
    integer                   :: n_names_old, n_names
  
    ! --- Choose which name must be added depending on units
    if (si_units) then
      final_name = name_si
    else
      final_name = name
    endif
  
    ! --- First added value or append
    if (.not. allocated(names_list) ) then
      allocate(names_list(1))
  
      names_list(1) = final_name 
    else
  
      n_names_old = size(names_list,1) 
      n_names = n_names_old + 1 
  
      allocate(names_tmp(n_names_old))
      names_tmp = names_list
  
      !append to the old vector
      deallocate(names_list)
      allocate(names_list(n_names))
      names_list(1:n_names_old) = names_tmp
      names_list(n_names)       = final_name
      
      !clean
      deallocate(names_tmp)
  
    endif
  
    counter = counter + 1 
    index   = counter
  
  end subroutine add_vtk_entry


end program jorek2vtk



