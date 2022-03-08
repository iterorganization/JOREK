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
type (type_element_list),     pointer :: element_list
type (type_bnd_element_list), pointer :: bnd_elm_list    
type (type_bnd_node_list),    pointer :: bnd_node_list 

integer               :: nnoel, nnos, nel, nsub, inode, ielm, n_scalars, n_vectors
real*4,allocatable    :: currdens(:), xyz (:,:), scalars(:,:), vectors(:,:,:)
integer,allocatable   :: ien (:,:)
integer, parameter    :: ivtk = 22 ! an arbitrary unit number for the VTK output file
integer               :: i, j, k, m, etype, irst, int, i_var, i_tor, i_tor_old, i_plane, index, index_node, my_id
character             :: buffer*80, lf*1, str1*12, str2*12
character*36, allocatable :: scalar_names(:), vector_names(:)
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
real*8                :: V0,  V0_s,  V0_t,  V0_st,  V0_ss,  V0_tt,  V, V_s, V_t, V_st, V_ss, V_tt
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
integer               :: n_fluxes, n_neo, n_bfield, n_vfield,n_pellet,n_bootstrap, n_psi_norm, n_Efield
integer               :: n_Jpol
integer               :: s_fluxes, s_neo, s_bfield, s_vfield,s_pellet,s_bootstrap, s_psi_norm, s_Efield
integer               :: s_Jpol
real*8                :: Jb,rho_norm,t_norm
integer               :: i_elm_axis, i_elm_xpoint(2), k_tor, ifail, ierr
logical               :: without_n0_mode, SI_units
logical               :: include_fluxes, include_neo, include_magnetic_field, include_velocity_field
logical               :: include_bootstrap, include_psi_norm, include_electric_field, include_Jpol, RphiZ_coords
real*8                :: toroidal_angle

real*8                :: Er, psi_abs, Vtheta, Btheta, Mach_par,Mach_pol,Vsound, Vneo
real*8                :: amu_neo_node, aki_neo_node
real*8                :: Vperp_e, Psi_tot

real*8                :: angle, source_volume, local_density, local_temperature, local_pressure, local_psi, local_source

logical               :: include_radiation
integer               :: n_radiation,s_radiation
real*8                :: Arad_bg, Brad_bg, Crad_bg, frad_bg
real*8                :: Te_eV, ne_SI, Lrad_imp, r_imp
real*8                :: Te_corr_eV, coef_rad_1, Sion_T, eta_Sp, ksiion, LradDcont_T
real*8                :: LradDrays_T, coef_ion_1, coef_ion_2, coef_ion_3, S_ion_puiss
real*8                :: r0_real8, rn0_real8
real*8                :: T0_corr, r0_corr, rn0_corr
integer               :: i_imp, offset_bgimp, i_bg     ! Loop for more than one background impurity


#ifdef WITH_Impurities
! See https://www.jorek.eu/wiki/doku.php?id=model500_501_555 for details
! Atomic physics coefficients:
!   -Mass ratio between main ions and impurites (m_i/m_imp)
real*8     :: m_i_over_m_imp
!   -Mean impurity ionization state
real*8     :: Z_imp, T0_Zimp, alpha_Zimp
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
integer               :: n_rn0, s_rn0
real*8                :: IonN, RecN, AblN, coef_rec_1, Srec_T

#ifdef fullmhd
!====================== --- Variables related to full mhd 
integer               :: n_fullmhd,s_fullmhd
#endif /*fullmhd*/

integer, parameter :: nplot = 200
integer :: iplot, i_elm
real*8  :: stmp(200)
real*8  :: Rp_start, Zp_start, Rp_end, Zp_end
real*8  :: Rp, Zp, Rmin, Rmax, Zmin, Zmax, s_out, t_out, R_out, Z_out

namelist /vtk_params/ nsub, i_tor, i_plane, without_n0_mode, SI_units, &
                      include_fluxes, include_neo, include_magnetic_field, include_velocity_field,&
                      include_bootstrap, include_psi_norm, include_electric_field, include_Jpol, RphiZ_coords


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
write(*,*) '***************************************'

call flush_it(6)

allocate(node_list)
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
RphiZ_coords           = .false. ! use xyz transformation (R,0,Z) instead of (R,Z,0)

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


write(*,*) '-----------'
write(*,*) 'n_tor           =', n_tor
write(*,*) 'n_period        =', n_period
write(*,*) 'F0              =', F0
write(*,*) 'R_geo,Z_geo     =', R_geo, Z_geo
write(*,*)
call flush_it(6)

! --- Number of scalars to write to the VTK output file
n_scalars   = n_var
n_vectors   = 0
n_fluxes    = 0
n_neo       = 0
n_bfield    = 0
n_vfield    = 0
n_Efield    = 0
n_Jpol      = 0
n_pellet    = 0
n_bootstrap = 0
n_psi_norm  = 0

if (include_fluxes) then
  n_fluxes  = 8
  s_fluxes  = n_scalars
  n_scalars = n_scalars + n_fluxes
endif
if (include_neo) then
  n_neo     = 10
  s_neo     = n_scalars
  n_scalars = n_scalars + n_neo
endif
if (include_magnetic_field) then
  n_bfield  = 1
  s_bfield  = n_vectors
  n_vectors = n_vectors + n_bfield
endif
if (include_velocity_field) then
  n_vfield  = 1
  s_vfield  = n_vectors
  n_vectors = n_vectors + n_vfield
endif
if (include_electric_field) then
  n_Efield  = 1
  s_Efield  = n_vectors
  n_vectors = n_vectors + n_Efield
endif
if (include_Jpol) then
  n_Jpol    = 1
  s_Jpol    = n_vectors
  n_vectors = n_vectors + n_Jpol
endif
if (use_pellet) then
  n_pellet  = 2  ! pellet and pressuren
  s_pellet  = n_scalars
  n_scalars = n_scalars + n_pellet
endif
if (include_bootstrap) then
  n_bootstrap = 2
  s_bootstrap = n_scalars
  n_scalars   = n_scalars + n_bootstrap
endif
if (include_psi_norm) then
   n_psi_norm = 1
   s_psi_norm = n_scalars
   n_scalars  = n_scalars + n_psi_norm
endif

#if (defined WITH_Neutrals) && (!defined WITH_Impurities)
    n_radiation = 0
 if (include_radiation) then
    n_radiation = 5 + n_adas
    s_radiation = n_scalars
    n_scalars   = n_scalars + n_radiation
    offset_bgimp = s_radiation+5 ! offset for plotting single bgimp species
  endif
  
    n_rn0 = 0
 if (include_neutral_dens) then
    n_rn0       = 2
    s_rn0       = n_scalars
    n_scalars   = n_scalars + n_rn0
  endif
#elif (!defined WITH_Impurities)
  if (include_radiation) then
    n_radiation = 1 + n_adas
    s_radiation = n_scalars
    n_scalars   = n_scalars + n_radiation
    offset_bgimp = s_radiation+1 ! offset for plotting single bgimp species
  end if
#endif

#ifdef WITH_Impurities
  n_radiation = 0
 if (include_radiation) then
    n_radiation  = 6 + n_adas-1 !(main impurity is not in background)
    s_radiation  = n_scalars
    n_scalars    = n_scalars + n_radiation
    offset_bgimp = s_radiation+6
  endif
#endif

#if fullmhd
 n_fullmhd = 11
 s_fullmhd = n_scalars
 n_scalars = n_scalars + n_fullmhd
#endif /*fullmhd*/

allocate(scalar_names(n_scalars), vector_names(n_vectors))

grad_psi = 0.d0

scalar_names(1:n_var) = variable_names((/(i, i=1,n_var)/))
if ( SI_units ) then
#ifdef fullmhd
   scalar_names(var_rho)='n_e20m-3    '
   scalar_names(var_T  )='Te_keV      '
   scalar_names(var_UR )='VR_km/s     '
   scalar_names(var_UZ )='VZ_km/s     '
   scalar_names(var_Up )='Vp_km/s     '
#else
   scalar_names(var_u)='u_m/s       '
   scalar_names(var_zj)='j_MA/m2     '
   scalar_names(var_rho)='n_e20m-3    '
   if (with_TiTe) then
      scalar_names(var_Ti)='Ti_keV      '
      scalar_names(var_Te)='Te_keV      '
   else
      scalar_names(var_T)='Te_keV      '
   endif
   if (with_Vpar) then
      scalar_names(var_Vpar)='Vpar_km/s   '
   endif
#endif

#if (defined WITH_Neutrals) || (defined WITH_Impurities)
   scalar_names(var_rhon)='N_dens_1d20  '
#endif

endif

if ( SI_units ) then

  if (include_fluxes) then
     scalar_names(s_fluxes+1:s_fluxes+n_fluxes) = (/ &
      'P_kPa       ', 'E_flux_Kpar ', 'E_flux_kperp', 'E_flux_Vpar ', &
      'E_flux_Vperp', 'D_flux_Dperp', 'D_flux_Vpar ', 'D_flux_Vperp'/)
  endif
  if (include_neo) then
     scalar_names(s_neo+1:s_neo+n_neo) = (/ &
      'Er_kV/m     ', 'Vtheta_km/s ', 'Mach_par    ', 'Mach_pol    ', &
      'Vsound_km/s ', 'Btot_T      ', 'Vneo_km/s   ', 'Vperp_e_km/s', &
      'ki_neo      ', 'mu_neo      '/)
  endif
  if (include_bootstrap)then
     scalar_names(s_bootstrap+1:s_bootstrap+n_bootstrap) = (/'j_b_MA/m2   ', 'j_av_MA/m2  '/)
  endif

  if (include_Jpol)  vector_names(s_Jpol  +1:s_Jpol  +n_Jpol  ) = 'Jpol (MA/m2)'

else

  if (include_fluxes) then
    scalar_names(s_fluxes+1:s_fluxes+n_fluxes) = (/ &
      'pressure    ', 'E_flux_Kpar ', 'E_flux_kperp', 'E_flux_Vpar ',&
      'E_flux_Vperp', 'D_flux_Dperp', 'D_flux_Vpar ', 'D_flux_Vperp'/)
  endif
  if (include_neo) then
    scalar_names(s_neo+1:s_neo+n_neo) = (/ &
      'Er          ', 'Vtheta      ', 'Mach_par    ', 'Mach_pol    ', &
      'Vsound      ', 'Btot        ', 'Vneo        ', 'Vperp_e     ', &
      'ki_neo      ', 'mu_neo      '/)
   endif
   if (include_bootstrap) then
      if (.not. bootstrap) write(*,*)'VTK WARNING: if you want the bootstrap, please set bootstrap=.t. in your input file!'
      scalar_names(s_bootstrap+1:s_bootstrap+n_bootstrap) = (/ 'j_bootstrap ', 'j_averaged  ' /)
   endif

   if (include_Jpol)  vector_names(s_Jpol  +1:s_Jpol  +n_Jpol  ) = 'Jpol'

!======================end SI units
endif

if (use_pellet) then
   scalar_names(s_pellet+1:s_pellet+n_pellet) =          (/ 'Pressure    ', 'Pellet      ' /)
endif
if (include_psi_norm) then
   scalar_names(s_psi_norm+1:s_psi_norm+n_psi_norm) = ('psi_norm    ')
endif

#if (defined WITH_Neutrals) && (!defined WITH_Impurities)
 if (include_radiation) then
  
   scalar_names(s_radiation+1:s_radiation+n_radiation-n_adas)                                   &
       = (/ 'Ionis_Wm-3  ', 'Lin_radWm-3 ', 'Brems_Wm-3  ', 'Joule_Wm-3  ', 'Imp_bg_Wm-3 '/)
   do i = 1,n_adas
     scalar_names(offset_bgimp + i) = 'Imp_bg_'//trim(imp_type(i))//'_Wm-3'
   end do
 endif

 if (include_neutral_dens) then
   scalar_names(s_rn0+1:s_rn0+n_rn0)&
                  = (/ 'IonN_s-1     ', 'RecN_s-1     '/)

 endif
#elif (!defined WITH_Impurities)
 if (include_radiation) then
   scalar_names(s_radiation+1:s_radiation+n_radiation)                                   &
       = (/  'Imp_bg_Wm-3 '/)
   do i = 1,n_adas
     scalar_names(offset_bgimp + i) = 'Imp_bg_'//trim(imp_type(i))//'_Wm-3'
   end do

 endif
#endif

#ifdef WITH_Impurities
 if (include_radiation) then
     scalar_names(s_radiation+1:s_radiation+n_radiation) &
                  = (/ 'Ionis_Jm-3  ', 'Coronal_radWm-3 ', 'Joule_Wm-3  ', 'Z_imp   ', 'Z_eff       ','Imp_bg_Wm-3 '/)
   i_bg=1
   do i = 1,n_adas
     if (i == index_main_imp) cycle
     scalar_names(offset_bgimp + i_bg) = 'Imp_bg_'//trim(imp_type(i))//'_Wm-3'
     i_bg=i_bg+1
   end do
 endif
#endif

#ifdef fullmhd
scalar_names(s_fullmhd+1:s_fullmhd+n_fullmhd) = (/  'B_R         ', 'B_Z         ', 'B_phi       ', &
                                                    'J_R         ', 'J_Z         ', 'J_phi       ', 'FFprime_equi', &
                                                    'F_prof_equi ', 'Grad_P      ', 'JxB         ', 'V_parallel  '/)
if ( SI_units ) then
scalar_names(s_fullmhd+1:s_fullmhd+11) = 'V_par_km/s  '
endif
#endif /*fullmhd*/

if (include_magnetic_field)  vector_names(s_bfield+1:s_bfield+n_bfield) = 'B_field' 
if (include_velocity_field)  vector_names(s_vfield+1:s_vfield+n_vfield) = 'v_field'
if (include_electric_field)  vector_names(s_Efield+1:s_Efield+n_Efield) = 'E_field_tmid'

do k_tor=1, n_tor
  mode(k_tor) = + int(k_tor / 2) * n_period
enddo


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
        call interp(node_list,element_list,i,var_T,  i_tor,s,t,T0, T0_s, T0_t, T0_st, T0_ss, T0_tt)

        if (with_Vpar) then
          call interp(node_list,element_list,i,var_Vpar,i_tor,s,t,V0,V0_s,V0_t,V0_st,V0_ss,V0_tt)
        else
          V0=0; V0_s=0; V0_t=0; V0_st=0; V0_ss=0; V0_tt=0
        end if

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
          scalars(inode,s_neo+1) = Er
          scalars(inode,s_neo+2) = Vtheta
          scalars(inode,s_neo+3) = Mach_par
          scalars(inode,s_neo+4) = Mach_pol
          scalars(inode,s_neo+5) = Vsound
          scalars(inode,s_neo+6) = Btot
          scalars(inode,s_neo+7) = Vneo
          scalars(inode,s_neo+8) = Vperp_e
          scalars(inode,s_neo+9) = aki_neo_node
          scalars(inode,s_neo+10) = amu_neo_node
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
            scalars(inode,s_neo+1) = Er
            scalars(inode,s_neo+2) = Vtheta
            scalars(inode,s_neo+3) = Mach_par
            scalars(inode,s_neo+4) = Mach_pol
            scalars(inode,s_neo+5) = Vsound
            scalars(inode,s_neo+6) = Btot
            scalars(inode,s_neo+7) = Vneo
            scalars(inode,s_neo+8) = Vperp_e
            if (NEO) then
               if (num_neo_file) then
                  scalars(inode,s_neo+9) = aki_neo_node
                  scalars(inode,s_neo+10) = amu_neo_node
               else
                  scalars(inode,s_neo+9) = aki_neo_const
                  scalars(inode,s_neo+10) = amu_neo_const
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
          call interp(node_list,element_list,i,var_T,  i_tor,s,t,TT,TT_s,TT_t,TT_st,TT_ss,TT_tt)
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

        scalars(inode,s_fullmhd+1:s_fullmhd+n_fullmhd) = 0.d0
        
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
          call interp(node_list,element_list,i,var_T,  i_tor,s,t,T0 ,T0_s, T0_t, T0_st, T0_ss, T0_tt)
          call interp(node_list,element_list,i,var_rho,i_tor,s,t,ZN0,ZN0_s,ZN0_t,ZN0_st,ZN0_ss,ZN0_tt)

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

        ZKpar_T = ZK_par * ((max(TT, T_min ))/T_0)**2.5

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
        
        scalars(inode,s_fullmhd+1) = BR
        scalars(inode,s_fullmhd+2) = BZ
        scalars(inode,s_fullmhd+3) = BP
        scalars(inode,s_fullmhd+4) = J_R
        scalars(inode,s_fullmhd+5) = J_Z
        scalars(inode,s_fullmhd+6) = J_phi
        scalars(inode,s_fullmhd+7) = F_prof*dF_dpsi
        scalars(inode,s_fullmhd+8) = F_prof 
        
        ! --- Choose which direction
        scalars(inode,s_fullmhd+9)  = GradP_pol ! GradP_p ! GradP_Z  ! GradP_R !
        scalars(inode,s_fullmhd+10) = JxB_pol   ! JxB_p   ! JxB_Z    ! JxB_R   !

        ! --- V_parallel
        scalars(inode,s_fullmhd+11)= (VR*BR + VZ*BZ + Vp*Bp)  / sqrt(BR**2 + BZ**2 + Bp**2)

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
          scalars(inode,s_bootstrap+1) = Jb
          scalars(inode,s_bootstrap+2) = J_phi * R ! to compare against the RMHD-303 current that has a factor R in it
        else
          Jb = 0.d0
        endif

        if (include_fluxes) then

          scalars(inode,s_fluxes+1)   = scalars(inode,var_rho) * scalars(inode,var_T)

          if (grad_psi .ne. 0.d0) then

            scalars(inode,s_fluxes+2)  = ZKpar_T * ( BR * TT_x + BZ * TT_y + BP * TT_P / R) / sqrt(BR**2 + BZ**2 + Bp**2)

            scalars(inode,s_fluxes+3)  = ZK_prof * (TT_x * A3_R + TT_y * A3_Z) / grad_psi

            scalars(inode,s_fluxes+4)  = rho * TT * (VR * BR + VZ * BZ + VP * BP) / sqrt(BR**2 + BZ**2 + Bp**2)

            scalars(inode,s_fluxes+5)  = ( (VZ*Bp-Vp*BZ)**2 + (Vp*BR-VR*Bp)**2 + (VR*BZ-VZ*BR)**2 ) **0.5
            scalars(inode,s_fluxes+5)  = rho * TT * scalars(inode,s_fluxes+5) / sqrt(BR**2 + BZ**2 + Bp**2)

            scalars(inode,s_fluxes+6)  = D_prof * (rho_x * A3_R + rho_y * A3_Z) / grad_psi

            scalars(inode,s_fluxes+7)  = rho * (VR * BR + VZ * BZ + VP * BP) / sqrt(BR**2 + BZ**2 + Bp**2)

            scalars(inode,s_fluxes+8)  = ( (VZ*Bp-Vp*BZ)**2 + (Vp*BR-VR*Bp)**2 + (VR*BZ-VZ*BR)**2 ) **0.5
            scalars(inode,s_fluxes+8)  = rho * scalars(inode,s_fluxes+8) / sqrt(BR**2 + BZ**2 + Bp**2)

          endif ! grad_psi
        endif ! include_fluxes

        if (include_psi_norm) then
          scalars(inode,s_psi_norm+1) = psi_norm
        endif

        if (include_magnetic_field) then
          vectors(inode,:,s_Bfield + 1) = (/ BR, BZ, Bp /)
        endif

        if (include_velocity_field) then
          vectors(inode,:,s_vfield + 1) = (/ VR, VZ, Vp /)
        endif

        if (include_electric_field) then
          E_R   = - (VZ*BP - VP*BZ) + eta_T * J_R
          E_Z   = - (VP*BR - VR*BP) + eta_T * J_Z
          E_phi = - (VR*BZ - VZ*BR) + eta_T * J_phi
          vectors(inode,:,s_Efield + 1) =  (/ E_R, E_Z, E_phi /)
        endif

        if (include_Jpol) then
          vectors(inode,:, s_Jpol  + 1) =  (/ J_R, J_Z, J_phi /)
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

        do i_tor = 1, n_tor

          if ( ( i_tor == 1 ) .and. ( without_n0_mode ) ) cycle ! Do not include the n=0 mode

          do m=1,n_var
             call interp(node_list,element_list,i,m,i_tor,s,t,P,P_s,P_t,P_st,P_ss,P_tt)
             scalars(inode,m) = scalars(inode,m) + P * HZ(i_tor,i_plane)
          enddo
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
          else
             call interp(node_list,element_list,i,var_T,  i_tor,s,t,TT ,TT_s, TT_t, TT_st, TT_ss, TT_tt)
          endif
         
          if (with_Vpar) then
             call interp(node_list,element_list,i,var_Vpar,i_tor,s,t,V,V_s,V_t,V_st,V_ss,V_tt)
          else
             V=0; V_s=0; V_t=0; V_st=0; V_ss=0; V_tt=0
          endif

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
          scalars(inode,s_bootstrap+1) = Jb
          scalars(inode,s_bootstrap+2) = bootstrap_spline3_eval(n_spline_vtk-1,psi_knots_vtk,j_knots_vtk,j_spline_vtk,psi_norm)
          ! --- The JOREK bootstrap is not constant on flux surface
          ! --- Because J_jorek = R*J_physical, and the physical bootstrap is constant on a surface
          ! --- Hence, it makes more sense to look at R*j_average to compare with the bootstrap...
          scalars(inode,s_bootstrap+2) = R* scalars(inode,n_var+2+n_fluxes+n_neo+n_pellet) / ES%R_axis
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

          scalars(inode,s_fluxes+1)   = scalars(inode,5) * scalars(inode,6)

          if (grad_psi .ne. 0.d0) then

            scalars(inode,s_fluxes+2)  = ZKpar_T * ( F0 * TT_p / BigR**2  + (TT_x * ps_y - TT_y * ps_x) / BigR ) / Btot

            scalars(inode,s_fluxes+3)  = ZK_prof * (TT_x * ps_x + TT_y * ps_y) / grad_psi

            scalars(inode,s_fluxes+4)  = scalars(inode,5) * scalars(inode,6) * scalars(inode,7) * Btot

            scalars(inode,s_fluxes+5)  = BigR   * (u_x * ps_y - u_y * ps_x) / sqrt(ps_x*ps_x + ps_y*ps_y) * scalars(inode,5) * scalars(inode,6)

            scalars(inode,s_fluxes+6)  = D_prof * (zn_x * ps_x + zn_y * ps_y) / sqrt(ps_x*ps_x + ps_y*ps_y)

            scalars(inode,s_fluxes+7)  = scalars(inode,5) * scalars(inode,7) * Btot

            scalars(inode,s_fluxes+8)  = BigR   * (u_x * ps_y - u_y * ps_x) / sqrt(ps_x*ps_x + ps_y*ps_y) * scalars(inode,5)

          endif ! grad_psi
        endif ! include_fluxes

        if (include_magnetic_field) then
          vectors(inode,:,s_Bfield + 1) = (/ ps_y/BigR, -ps_x/BigR, F0/BigR /)          
        endif

        if (include_velocity_field) then
          vectors(inode,:,s_vfield + 1) = (/ -BigR*u_y + V/BigR*ps_y, BigR*u_x - V/BigR*ps_x, V*F0/BigR /)          
        endif

        if (include_electric_field) then
          vectors(inode,:,s_Efield + 1) =  (/ E_R, E_Z, E_phi /)
        endif

        if (include_Jpol) then
          call J_pol(node_list, element_list, i, s, t, i_plane,.false., Jpol_R, Jpol_Z, FFp)
          vectors(inode,:, s_Jpol  + 1) =  (/ Jpol_R, Jpol_Z, 0.d0 /)
        endif
        
        if (include_psi_norm) then
           scalars(inode,s_psi_norm+1) = psi_norm
        endif

        if (use_pellet) then

           local_density     = scalars(inode,5)
           local_temperature = scalars(inode,6)/2.
           local_psi         = scalars(inode,1)
           local_pressure    = local_density * local_temperature
           scalars(inode,s_pellet+1)  = local_pressure

           angle = 0.0
           call pellet_source2(pellet_amplitude,pellet_R,pellet_Z,pellet_psi,pellet_phi, &
                pellet_radius, pellet_delta_psi, pellet_sig, pellet_length, pellet_ellipse, &
                pellet_theta, R, Z, local_psi, angle, &
                local_density,local_temperature, &
                central_density, pellet_particles, pellet_density, &
                total_pellet_volume, local_source, source_volume)

           scalars(inode,n_var+n_fluxes+n_neo+n_pellet) = local_source
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
        Te_corr_eV = corr_neg_temp(T_real8)/(EL_CHG*MU_ZERO*central_density*1.d20)
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

      scalars(i,s_radiation+1) = ksiion * scalars(i,var_rho) * scalars(i,var_rhon) * Sion_T
      scalars(i,s_radiation+2) = scalars(i,var_rho) * scalars(i,var_rhon) * LradDrays_T
      scalars(i,s_radiation+3) = LradDcont_T * scalars(i,var_rho)**2.d0
      scalars(i,s_radiation+4) = (2/(3 * BigR**2)) * eta_Sp * scalars(i,var_zj)**2.d0
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
          scalars(i,offset_bgimp + i_imp) = r_imp * Lrad_imp * scalars(i,var_rho)
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


      scalars(i,s_radiation+n_radiation-n_adas) = scalars(i,var_rho) * frad_bg
   
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
       Te_corr_eV = corr_neg_temp(T_real8,(/5.d-1,5.d-1/),max(T_min,Te_1))/(EL_CHG*MU_ZERO*central_density*1.d20)
       Te_eV = T_real8/(EL_CHG*MU_ZERO*central_density*1.d20)
     else
       T_real8 = scalars(i,var_T)
       Te_corr_eV = corr_neg_temp(T_real8,(/5.d-1,5.d-1/),max(2.*T_min,2.*T_1))/(2.d0*EL_CHG*MU_ZERO*central_density*1.d20)
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

     if (allocated(imp_adas(index_main_imp)%ionisation_energy)) then
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
       E_ion     = E_ion * EL_CHG*MU_ZERO*central_density*1.d20
     else
       call imp_cor(index_main_imp)%interp_linear(density=20.,temperature=log10(Te_corr_eV*EL_CHG/K_BOLTZ),z_avg=Z_imp)
       E_ion     = 0.
     end if

     alpha_imp    = 0.5*m_i_over_m_imp*(Z_imp+1.) - 1.
     beta_imp     = m_i_over_m_imp*Z_imp - 1.

     ne_SI       = (r0_corr + beta_imp * rn0_corr) * 1.d20 * central_density ! electron density (SI)
     scalars(i,var_rho) = (r0_corr + beta_imp * rn0_corr)                           ! electron density (JOREK units)

     !Calculate the Z_eff, as it is done in mod_elt_matrix
     Z_eff = r0_corr - rn0_corr
     do ion_i=1, imp_adas(index_main_imp)%n_Z
       Z_eff = Z_eff + m_i_over_m_imp * rn0_corr * P_imp(ion_i) * real(ion_i,8)**2
     end do
     Z_eff = Z_eff / scalars(i,var_rho)  
     scalars(inode,s_radiation+5) = Z_eff

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
     i_bg = 1
     do i_imp =1, n_adas
       if (i_imp == index_main_imp) cycle

       r_imp = nimp_bg(i_imp) / (1.d20 * central_density)  ! Background impurity density in JU     
       if (ne_SI > ne_SI_min .and. Te_eV > Te_eV_min .and. r_imp > 0) then
         Lrad_imp = 0.0
         call radiation_function_linear(imp_adas(i_imp),imp_cor(i_imp),log10(ne_SI),    & 
                                        log10(Te_eV*EL_CHG/K_BOLTZ),.true.,Lrad_imp)           
       else     
         Lrad_imp = 0.
       end if
       frad_bg = frad_bg + r_imp * Lrad_imp
       scalars(i,offset_bgimp+i_bg) = scalars(i,var_rho) * r_imp * Lrad_imp
       i_bg = i_bg+1
     end do
     scalars(i,s_radiation+1) = (2./3.) * scalars(i,var_rhon) * E_ion
     scalars(i,s_radiation+2) = (r0_corr+beta_imp*rn0_corr) * rn0_corr * Lrad
     scalars(i,s_radiation+3) = (2./(3. * BigR**2)) * eta_Sp * scalars(i,var_zj)**2.d0
     scalars(i,s_radiation+4) = Z_imp
     scalars(i,s_radiation+5) = Z_eff
     scalars(i,s_radiation+6) = scalars(i,var_rho) * frad_bg

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

      scalars(i,s_rn0+1) = IonN
      scalars(i,s_rn0+2) = RecN

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
    scalars(i,var_T  ) = scalars(i,var_T)   / MU_zero / (central_density * 1d20) / EL_CHG /2./1.e3 !(assumes Te=Ti=T/2)
    !===========================================Velocity
    scalars(i,var_UR) = scalars(i,var_UR) /t_norm/1.e3
    scalars(i,var_UZ) = scalars(i,var_UZ) /t_norm/1.e3
    scalars(i,var_Up) = scalars(i,var_Up) /t_norm/1.e3
    scalars(i,s_fullmhd+11) = scalars(i,s_fullmhd+11) /t_norm/1.e3 ! V_parallel
    !=====================Pressure in kPa
    if (include_fluxes) scalars(i,s_fluxes+1) = scalars(i,s_fluxes+1) / MU_zero/1.e3
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
    scalars(i,s_bootstrap+1)=scalars(i,s_bootstrap+1)/MU_zero*1.e-6
    scalars(i,s_bootstrap+2)=scalars(i,s_bootstrap+2)/MU_zero*1.e-6
    endif
    if (include_velocity_field) then 
      vectors(i,:,s_vfield + 1) = vectors(i,:,s_vfield + 1)/t_norm
    endif
    if (include_electric_field) then 
      vectors(i,:,s_Efield + 1) = vectors(i,:,s_Efield + 1)/t_norm
    endif
    if (include_Jpol) then
      vectors(i,:, s_Jpol  + 1) = vectors(i,:,s_Jpol   + 1)/MU_zero*1e-6
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
    if (include_fluxes) scalars(i,s_fluxes+1) = scalars(i,s_fluxes+1) / MU_zero/1.e3
    if (include_neo) then
      !============================Er in kV/m
      scalars(i,s_neo+ 1) = F0*scalars(i,s_neo+ 1) / t_norm/1.e3
      !====================================Vtheta km/s
      scalars(i,s_neo+ 2) =    scalars(i,s_neo+ 2) / t_norm/1.e3
      !===================================Vsound in km/s
      scalars(i,s_neo+ 5) =    scalars(i,s_neo+ 5) / t_norm/1.e3
      !===================================Vneo in km/s
      scalars(i,s_neo+ 7) =    scalars(i,s_neo+ 7) / t_norm/1.e3
      !===================================Vperp_e in km/s
      scalars(i,s_neo+ 8) =    scalars(i,s_neo+ 8) / t_norm/1.e3
      ! ===================================mu_neo in SI units
      scalars(i,s_neo+10) =    scalars(i,s_neo+10) / sqrt(rho_norm*MU_zero)
    endif
    !============================================j_bootstrap, javeraged in MA/m2

    if (include_bootstrap) then
    scalars(i,s_bootstrap+1)=scalars(i,s_bootstrap+1)/MU_zero*1.e-6
    scalars(i,s_bootstrap+2)=scalars(i,s_bootstrap+2)/MU_zero*1.e-6
    endif
    if (include_velocity_field) then 
      vectors(i,:,s_vfield + 1) = vectors(i,:,s_vfield + 1)/t_norm
    endif
    if (include_electric_field) then 
      vectors(i,:,s_Efield + 1) = vectors(i,:,s_Efield + 1)/t_norm
    endif
    if (include_Jpol) then
      vectors(i,:, s_Jpol  + 1) = vectors(i,:,s_Jpol   + 1)/MU_zero*1e-6
    endif
 
    !========================================================

#if (!defined WITH_Impurities)
    if (include_radiation) then
      r0_real8  = scalars(i,var_rho)/central_density ! Back to JOREK unit for calling atomic_coeff_deuterium
      if ( with_TiTe ) then
        T_real8 = scalars(i,var_Te)*1.e3*EL_CHG*MU_zero*(central_density * 1.d20) ! T_real8 back to JOREK units
        Te_corr_eV = corr_neg_temp(T_real8)/(EL_CHG*MU_ZERO*central_density*1.d20)
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
  
      scalars(i,s_radiation+1) = ksiion* (1.5d0)/(MU_zero*central_density*1.d20)      &
                                          * scalars(i,var_rho) * 1.d20 * scalars(i,var_rhon) * 1.d20 * Sion_T / coef_ion_1

      scalars(i,s_radiation+2) = scalars(i,var_rho)* 1.d20 * scalars(i,var_rhon) * 1.d20 * LradDrays_T/ coef_rad_1

      scalars(i,s_radiation+3) = LradDcont_T * (scalars(i,var_rho)*1.d20)**2.d0 / coef_rad_1

      scalars(i,s_radiation+4) = eta_Sp * (1.d6*scalars(i,var_zj))**2.d0
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
          scalars(i,offset_bgimp + i_imp) = 1.d20 * scalars(i,var_rho) * nimp_bg(i_imp) * Lrad_imp 
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
      scalars(i,s_radiation+n_radiation-n_adas) = scalars(i,var_rho)*1.d20 * frad_bg
    endif
#endif /*(.not. with_Impurities)*/

#ifdef WITH_Impurities
  if (include_radiation) then
   scalars(i,s_radiation+1) = scalars(i,s_radiation+1)/(K_BOLTZ*MU_ZERO)
   scalars(i,s_radiation+2) = scalars(i,s_radiation+2)/(2.d0/3.d0*MU_ZERO**1.5d0*(central_mass*MASS_PROTON*central_density*1.d20)**0.5d0)
   scalars(i,s_radiation+3) = scalars(i,s_radiation+3)/(2.d0/3.d0*((central_mass*MASS_PROTON*central_density*1.d20)**0.5)*(MU_ZERO**1.5)) 
   scalars(i,s_radiation+4) = scalars(i,s_radiation+4)
   scalars(i,s_radiation+6) = scalars(i,s_radiation+6) &
       /((GAMMA-1.d0)*MU_ZERO**1.5d0*(central_mass*MASS_PROTON*central_density*1.d20)**0.5d0)
   do i_imp=1,n_adas
     scalars(i,offset_bgimp+i_imp) = scalars(i,offset_bgimp+i_imp) / ((GAMMA-1.d0)*MU_ZERO**1.5d0*(central_mass*MASS_PROTON*central_density*1.d20)**0.5d0)
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

end program jorek2vtk
