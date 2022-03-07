!> Calculates 3D integrals and boundary fluxes.
!! The NOMPIVERSION preprocessor flag is needed for mod_integrals3D_nompi, a version of this subroutine not requiring MPI.
#ifndef NOMPIVERSION
module mod_integrals3D
#endif

  use constants
  use mod_parameters
  use data_structure
  use gauss
  use basis_at_gaussian
  use tr_module
  use phys_module
  use mod_interp
  use convert_character
  use mpi_mod
  use mod_expression
  use mod_resistivity
  use mod_poloidal_currents, only : integrated_normal_bnd_curr 
  use corr_neg
  use pellet_module
#if (defined WITH_Neutrals) && (!defined WITH_Impurities)
  use mod_neutral_source, only: total_neutral_source, total_n_particles, total_n_particles_inj, total_n_particles_inj_all
  use mod_source_shape, only: source_shape
#endif
#ifdef WITH_Impurities
  use mod_injection_source, only: total_imp_source, total_n_particles, total_n_particles_inj, total_n_particles_inj_all
  use mod_source_shape, only: source_shape
#endif
  use mod_impurity, only: radiation_function, radiation_function_linear
  use equil_info, only : get_psi_n, ES
  use mod_atomic_coeff_deuterium, only : atomic_coeff_deuterium
  use mod_sources

  implicit none
  
  private
  
  public :: int3d_new 
  
  contains


subroutine int3d_new(my_id, node_list, element_list, bnd_node_list, bnd_elm_list, expr_list, res, units)

!$ use omp_lib
 
implicit none

type (type_node_list),        intent(in)    :: node_list
type (type_element_list),     intent(in)    :: element_list   
type (type_bnd_node_list),    intent(in)    :: bnd_node_list
type (type_bnd_element_list), intent(in)    :: bnd_elm_list   
type (t_expr_list),           intent(in)    :: expr_list
real*8,                    intent(inout)    :: res(:)
integer,                      intent(in)    :: units

! --- Local variables
type (type_element)      :: element, elm_k
type (type_node)         :: nodes(n_vertex_max), node_k
type (type_bnd_element)  :: bndelem
type (type_surface_list) :: surface_list

real*8  :: psi_axis, psi_bnd
real*8  :: x_g(n_gauss,n_gauss),        x_s(n_gauss,n_gauss),        x_t(n_gauss,n_gauss),        x_ss(n_gauss,n_gauss),        x_tt(n_gauss,n_gauss),        x_st(n_gauss,n_gauss)
real*8  :: y_g(n_gauss,n_gauss),        y_s(n_gauss,n_gauss),        y_t(n_gauss,n_gauss),        y_ss(n_gauss,n_gauss),        y_tt(n_gauss,n_gauss),        y_st(n_gauss,n_gauss)
real*8  :: eq_g(n_plane,0:n_var,n_gauss,n_gauss), eq_s(n_plane,0:n_var,n_gauss,n_gauss)
real*8  :: eq_t(n_plane,0:n_var,n_gauss,n_gauss), eq_p(n_plane,0:n_var,n_gauss,n_gauss)
real*8  :: eq_ss(n_plane,0:n_var,n_gauss,n_gauss), eq_tt(n_plane,0:n_var,n_gauss,n_gauss), eq_st(n_plane,0:n_var,n_gauss,n_gauss)
real*8  :: eq_sp(n_plane,0:n_var,n_gauss,n_gauss), eq_tp(n_plane,0:n_var,n_gauss,n_gauss)
real*8  :: wgauss_copy(n_gauss)
real*8  :: psi_axisym(n_gauss,n_gauss)

real*8  :: x_g_1D(n_gauss),  x_s_1D(n_gauss),   x_t_1D(n_gauss)
real*8  :: y_g_1D(n_gauss),  y_s_1D(n_gauss),   y_t_1D(n_gauss)
real*8  :: eq_g_1D(n_plane,0:n_var,n_gauss), eq_s_1D(n_plane,0:n_var,n_gauss)
real*8  :: delta_g_1D(n_plane,0:n_var,n_gauss)
real*8  :: eq_t_1D(n_plane,0:n_var,n_gauss), eq_p_1D(n_plane,0:n_var,n_gauss)

real*8  :: current_source, particle_source, heat_source, heat_source_i, heat_source_e, rotation_source
real*8  :: xt, t_norm, rho_norm, t_norm2
real*8  :: eq_zne(n_gauss,n_gauss), eq_zTe(n_gauss,n_gauss)
real*8  :: dn_dpsi,dn_dz,dn_dpsi2,dn_dz2,dn_dpsi_dz,dn_dpsi3,dn_dpsi_dz2, dn_dpsi2_dz
real*8  :: dT_dpsi,dT_dz,dT_dpsi2,dT_dz2,dT_dpsi_dz,dT_dpsi3,dT_dpsi_dz2, dT_dpsi2_dz

integer :: i, j, k, in, ms, mt, mp, iv, inode, ife, n_elements, i_elm_axis, i_elm_xpoint(2), ifail
integer :: ierr, n_cpu, my_id, ife_delta, ife_min, ife_max, omp_nthreads, omp_tid
integer :: k_vertex, k_dof, k_node, k_dir, k_dir_perp, m_bndelem, dir_perp(2), mv1, m_elm
integer :: iexpr
real*8  :: R_c, Z_c, vec_inside(2), grad_t(2)
real*8  :: k_size, k_size_perp
real*8  :: G(4,4), sign_out, psi_n, ps0_sbnd, u0_sbnd
real*8  :: dt_back, dt_now, r_dt, r_dt2
real*8  :: I_halo, TPF, q02, q95, q99
real*8, allocatable :: qval(:), radav(:)

real*8  :: R_axis,Z_axis,s_axis,t_axis
real*8  :: current_tot, beta_p, beta_n, beta_t, aminor, current_MA
real*8  :: xjac, xjac_R, xjac_Z, BigR, wst, P_int, P_e_int, P_i_int, C_intern, zj0, ps0, r0, T0, T0e, T0i
real*8  :: Vol, Volume, Area, Bgeo, area1 
real*8  :: psi_as_coord
real*8  :: AR0, AR0_p, AR0_s, AR0_t, AR0_sp, AR0_tp, AR0_Rp, AZ0, AZ0_p, AZ0_s, AZ0_t, AZ0_sp, AZ0_tp, AZ0_Zp, A30
real*8  :: A30_p, A30_s, A30_t, A30_ss, A30_tt, A30_st, A30_R, A30_RR, A30_ZZ
real*8  :: BR_Z, BZ_R
real*8  :: r0_corr, T0_corr, T0e_corr, T0i_corr
real*8  :: density_tot, density_in, density_out,  pressure, pressure_in, pressure_out
real*8  :: pressure_e, pressure_e_in, pressure_e_out, pressure_i, pressure_i_in, pressure_i_out
real*8  :: current_in, current_out, D_int, D_ext, P_ext, C_ext, delta_phi, phi, P_tot, D_tot
real*8  :: P_e_ext, P_i_ext, P_e_tot, P_i_tot
real*8  :: VP_int, VP_ext, VK_int, VK_ext, vpar0, BB2, VP_tot, VK_tot
real*8  :: kin_par_in, kin_par_out, kin_par_tot, kin_perp_in, kin_perp_out, kin_perp_tot
real*8  :: VM_int, VM_ext, VM_tot, mag_in, mag_out, mag_tot, J2_int, J2_ext, J2_tot, ohm_in, ohm_tot, ohm_out
real*8  :: heli_tot, thm_wk, thm_wk_tot, mag_wk, mag_wk_tot, thermal_work_tot
real*8  :: vpar_disp_tot, vpar_disp, viscopar_dissip_tot, source_tot, heating_tot
real*8  :: fric_disp_tot, fric_disp, friction_dissip_tot
real*8  :: H_int, H_ext, S_int, S_ext, heating_in, heating_out, source_in, source_out
real*8  :: psi_xpoint(2),R_xpoint(2),Z_xpoint(2),s_xpoint(2),t_xpoint(2)
real*8  :: dTdx, dTdy, drhodx, drhody, dPdx, dPdy, dpsidx, dpsidy, dudx, dudy, drhondx, drhondy
real*8  :: dTedx, dTedy, dTidx, dTidy, dPedx, dPedy, dPidx, dPidy
real*8  :: source_volume, source_pellet, eta_T, eta_T_ohm
real*8  :: local_pellet_particles, local_plasma_particles, local_pellet_volume
real*8  :: local_n_particles_inj, local_n_particles, rn0, rn0_corr, neut_particles_tot
real*8  :: E_tot, E_in, E_out, Zkpar_T, D_prof, ZK_prof, sheath_heatflux
real*8  :: ZK_e_prof, ZK_i_prof, ZK_e_par_T, ZK_i_par_T
real*8  :: fact_mu0, fact_flux, fact_part
real*8  :: hel1, heli, helicity_tot, psi_off, curr, Ip, vn_p0, qn, pflow, kinflow, cond_par, cond_perp
real*8  :: kinpar_flux, qn_par, qn_perp, mag_work_tot, mag_src_tot, mag_source_tot
real*8  :: vpar_part_flux, vperp_part_flux, Dperp_part_flux, Dpar_part_flux, neut_part_flux
real*8  :: vpar_part_flow, vperp_part_flow, Dperp_part_flow, Dpar_part_flow, neut_part_flow
real*8  :: poynting_flux, poynting_tmp, dpsi_dt
real*8  :: s_or_t,sg,tg,R,R_s,R_t,R_st,R_ss,R_tt,Z,Z_s,Z_t,Z_st,Z_ss,Z_tt
real*8  :: RH,RH_s,RH_t,RH_st,RH_ss,RH_tt
real*8  :: TT,TT_s,TT_t,TT_st,TT_ss,TT_tt 
real*8  :: UU,UU_s,UU_t,UU_st,UU_ss,UU_tt 
real*8  :: PS,PS_s,PS_t,PS_st,PS_ss,PS_tt 
real*8  :: vp,vp_s,vp_t,vp_st,vp_ss,vp_tt 
real*8  :: rn,rn_s,rn_t,rn_st,rn_ss,rn_tt 
real*8  :: psi_s, psi_t, rho_s, rho_t, T_s, T_t, Ti, Ti_s, Ti_t, Te, Te_s, Te_t, p0_s, p0_t, u0_s, u0_t, ps0_s, ps0_t, p0_p, rhon_s, rhon_t
real*8  :: u0_p, u_s, u_t, u_p
real*8  :: u0_x, u0_y
real*8  :: viscopar_flux, viscopar_f, vpar_s, vpar_t, vpar_x, vpar_y, li3_tot, li3
real*8  :: varmin(n_var), varmax(n_var), V_min(n_var), V_max(n_var)
real*8  :: R_curr_cent, Z_curr_cent, Zcurr_tmp, R2curr_tmp, R2curr


#if (defined WITH_Neutrals) && (!defined WITH_Impurities)
real*8  :: source_neutral
#endif
#ifdef WITH_Impurities
real*8  :: source_bg, source_imp
#endif
#if (defined WITH_Neutrals) || (defined WITH_Impurities)
real*8  :: local_radiation, local_E_ion, total_radiation, total_E_ion, local_P_ei, total_P_ei
real*8  :: local_P_ion, total_P_ion
real*8  :: local_radiation_phi(n_plane), total_radiation_phi(n_plane)
real*8  :: ne_SI, Te_eV, Te_corr_eV, Ti_eV

! SPI-related variables
integer    :: spi_i
integer    :: i_inj,  n_spi_tmp
real*8     :: spi_R_tmp
real*8     :: spi_Z_tmp
real*8     :: spi_phi_tmp
real*8     :: spi_psi_tmp
real*8     :: spi_grad_psi_tmp
real*8     :: ng_radius_tmp !< Radius of neutral gas cloud as a result of the ablation
real*8     :: source_tmp
real*8     :: ns_shape ! variable for numerical integration of source volume
real*8     :: V_ns, V_ns_drift
real*8, allocatable :: local_source_volume(:), local_source_volume_drift(:)

#endif

! Additional diagnostic variables for impurity model
! See https://www.jorek.eu/wiki/doku.php?id=model500_501_555 for details
#ifdef WITH_Impurities
! Atomic physics coefficients:
!   -Mass ratio between main ions and impurites (m_i/m_imp)
real*8  :: m_i_over_m_imp, m_imp
!   -Mean impurity ionization state
real*8  :: Z_imp, dZ_imp_dT, Z_eff, eta_coef, ne_JOREK, dne_JOREK_dx, dne_JOREK_dy
real*8  :: Z_eff_imp
!   -Corrected plasma temperature and density for radiation calculation
real*8  :: dT0e_corr_dT, Ti_corr_eV
!   -Temporary variable for charge state distribution
real*8, allocatable :: P_imp(:)
real*8     :: E_ion, Lrad, E_ion_bg
real*8     :: Lrad_imp, frad_bg
integer    :: ion_i, ion_k, i_phi, i_imp
#endif
#ifdef WITH_Impurities
#ifdef WITH_TiTe
!   -Coefficients related to Z_imp
real*8  :: alpha_i, alpha_e, dalpha_e_dT

!   -Ion-electron energy transfer
real*8     :: nu_e_imp, nu_e_bg, lambda_e_imp, lambda_e_bg, dTi_e, dTe_i
#else /* WITH_TiTe */
!   -Coefficients related to Z_imp
real*8  :: alpha_imp, dalpha_imp_dT, beta_imp, dbeta_imp_dT
#endif /* WITH_TiTe */
#endif /* WITH_Impurities */

! Additional variables related to the radiated power
#if (defined WITH_Neutrals) && (!defined WITH_Impurities)
! Atomic physics coefficients:
!   -Ionization
real*8     :: Sion_T, dSion_dT                                ! Ionization rate and its derivative wrt. temperature
real*8     :: ksiion                                          ! Ionization energy
!   -Recombination
real*8     :: Srec_T, dSrec_dT                                ! Recombination rate and its derivative wrt. temperature
!   -Radiation from injected gas/impurities
real*8     :: LradDrays_T, dLradDrays_dT                      ! Line (/rays) radiation rate and its derivative wrt. temperature
real*8     :: LradDcont_T, dLradDcont_dT                      ! Continuum (Brem.) radiation rate and its derivative wrt. T
!   -Radiation from background impurities
real*8     :: Arad_bg, Brad_bg, Crad_bg, frad_bg              ! Retain hard-coded fitting for argon
real*8     :: Lrad_imp
integer*8  :: i_phi
real*8     :: coef_prad_si                                    ! Prad,SI = coef_prad_si * Prad,jorek
integer    :: i_imp                                           ! Loop for more than one background impurity
#endif

#ifndef NOMPIVERSION
call MPI_COMM_SIZE(MPI_COMM_WORLD, n_cpu, ierr) ! number of MPI procs
n_cpu = max(n_cpu,1)
#else
n_cpu = 1
#endif

if (my_id .eq. 0) then
  write(*,*) '***************************************'
  write(*,*) '* Integrals  (3D)                     *'
  write(*,*) '***************************************'
  !write(*,*) ' n_plane : ',n_plane
  !write(*,*) ' n_cpu   : ',n_cpu
endif

density_tot  = 0.d0
pressure = 0.d0
pressure_i = 0.d0
pressure_e = 0.d0
D_int    = 0.d0
P_int    = 0.d0
P_e_int  = 0.d0
P_i_int  = 0.d0
C_intern = 0.d0
H_int    = 0.d0
S_int    = 0.d0
VP_int   = 0.d0
VK_int   = 0.d0
VM_int   = 0.d0
J2_int   = 0.d0
D_ext    = 0.d0
P_ext    = 0.d0
P_e_ext  = 0.d0
P_i_ext  = 0.d0
C_ext    = 0.d0
H_ext    = 0.d0
S_ext    = 0.d0
VP_ext   = 0.d0
VK_ext   = 0.d0
VM_ext   = 0.d0
J2_ext   = 0.d0
Vol      = 0.d0
area1    = 0.d0
P_tot    = 0.d0
P_e_tot  = 0.d0
P_i_tot  = 0.d0
D_tot    = 0.d0
wgauss_copy = wgauss
VP_tot   = 0.d0
VK_tot   = 0.d0
VM_tot   = 0.d0
J2_tot   = 0.d0
hel1     = 0.d0
heli_tot = 0.d0
vn_p0    = 0.d0
qn_par   = 0.d0
qn_perp  = 0.d0
kinpar_flux  = 0.d0
poynting_flux= 0.d0
viscopar_flux= 0.d0
vpar_disp_tot= 0.d0
fric_disp_tot= 0.d0
psi_off      = 0.d0
thm_wk_tot   = 0.d0
mag_wk_tot   = 0.d0
mag_src_tot  = 0.d0
varmin   = +1.d99
varmax   = -1.d99
R2curr_tmp = 0.d0
Zcurr_tmp  = 0.d0

Dpar_part_flux   = 0.d0 
Dperp_part_flux  = 0.d0
vpar_part_flux   = 0.d0
vperp_part_flux  = 0.d0
neut_part_flux   = 0.d0

local_pellet_particles = 0.d0
local_plasma_particles = 0.d0
local_pellet_volume    = 0.d0

local_n_particles_inj = 0.d0
local_n_particles     = 0.d0

#if (defined WITH_Neutrals) || (defined WITH_Impurities)
local_radiation       = 0.d0
local_radiation_phi   = 0.d0
local_E_ion           = 0.d0
local_P_ei            = 0.d0
local_P_ion           = 0.d0

if (using_spi) then
   if (allocated(local_source_volume)) then
      deallocate(local_source_volume)
   end if
   if (allocated(local_source_volume_drift)) then
      deallocate(local_source_volume_drift)
   end if

   allocate (local_source_volume(n_spi_tot))
   allocate (local_source_volume_drift(n_spi_tot))

   do spi_i=1, n_spi_tot
      local_source_volume(spi_i)            = 0.d0
      local_source_volume_drift(spi_i)      = 0.d0
   end do
end if
if (.not. allocated(local_source_volume)) allocate (local_source_volume(1)) ! Allocate a dummy array for omp
if (.not. allocated(local_source_volume_drift)) allocate (local_source_volume_drift(1)) 
#endif

delta_phi     = 2.d0 * PI / float(n_plane) / float(n_period)

psi_axis   = ES%psi_axis;        R_axis = ES%R_axis;        Z_axis = ES%Z_axis
psi_xpoint = ES%psi_xpoint;    R_xpoint = ES%R_xpoint;    Z_xpoint = ES%Z_xpoint 
psi_bnd    = ES%psi_bnd

ife_delta = ceiling(float(element_list%n_elements) / n_cpu)
ife_min   =      my_id     * ife_delta + 1
ife_max   = min((my_id +1) * ife_delta, element_list%n_elements)

!$omp parallel default(none)                                                                   &
!$omp   shared(element_list,node_list, H, H_s, H_t, HZ, HZ_p, ife_min, ife_max, xpoint, xcase, &
!$omp          H_ss, H_tt, H_st,                                                               &
!$omp          R_xpoint, Z_xpoint, my_id, use_pellet, delta_phi, R_axis, Z_axis, psi_axis, psi_bnd, &
!$omp          D_tot, D_int, D_Ext, P_tot, P_int, P_ext, Vol, C_intern, C_ext, VP_ext, VP_int, &
!$omp          VK_ext, VK_int, VK_tot, VM_ext, VM_int, VM_tot, J2_tot, J2_ext, J2_int,         &
!$omp          H_int, H_ext, S_int, S_ext,psi_xpoint,  F0, VP_tot,eta, T_0, Te_0, T_min,       &
!$omp          ne_SI_min, Te_eV_min, rn0_min, P_e_tot, P_i_tot, P_e_int, P_i_int, P_e_ext, P_i_ext, &
!$omp          pellet_amplitude,pellet_R,pellet_Z,pellet_psi,pellet_phi,                       &
!$omp          pellet_radius, pellet_delta_psi, pellet_sig, pellet_length, pellet_ellipse, pellet_theta,  &
!$omp          central_density, pellet_particles,pellet_density, pellet_volume,                &
!$omp          local_pellet_particles, local_plasma_particles, local_pellet_volume,            &
!$omp          heli_tot,  keep_current_prof, psi_off, visco_par, thm_wk_tot,                   &
!$omp          mag_wk_tot, vpar_disp_tot, fric_disp_tot, area1, mag_src_tot,                   &
!$omp          eta_ohmic, central_mass, R2curr_tmp, Zcurr_tmp,                                 &
#if (defined WITH_Neutrals) || (defined WITH_Impurities)
!$omp          spi_num_vol, local_source_volume, local_source_volume_drift, drift_distance,    &
!$omp          using_spi, n_spi_tot, n_inj, n_spi,                                             &
!$omp          pellets, ng_radius_ratio, ng_radius_min,                                        &
!$omp          local_n_particles_inj, local_n_particles, ns_amplitude, ns_R, ns_Z,             &
!$omp          ns_phi, ns_radius, ns_deltaphi, ns_deltaminrad, ns_tor_norm, spi_tor_rot, local_E_ion,          &
!$omp          t_now, A_Dmv, K_Dmv, V_Dmv, P_Dmv, t_ns, L_tube, JET_MGI,ASDEX_MGI, local_P_ion,&
!$omp          local_radiation, local_radiation_phi, imp_cor, imp_adas, imp_type, local_P_ei,  &
!$omp          n_adas, nimp_bg,                                                                &
#endif
#if (defined WITH_Neutrals) && (!defined WITH_Impurities)
!$omp          ksi_ion, GAMMA, use_imp_adas,                                                   &
#endif
#if (defined WITH_Impurities)
!$omp          index_main_imp,                                                                 &
#endif
!$omp          T_1, T_max_eta, T_max_eta_ohm, eta_T_dependent,                                 &
!$omp          wgauss_copy, varmin, varmax)                                                    &
!$omp   private(ife,iv,inode,element,nodes,i,j, k,in, mp, ms, mt,                              &
!$omp           x_g, y_g, x_s, y_s, x_t, y_t, xjac, xjac_R, xjac_Z, eq_g, eq_s, eq_t, eq_p,    &
!$omp           x_ss, x_tt, x_st, y_ss, y_tt, y_st, eq_ss, eq_tt, eq_st, eq_sp, eq_tp,         &
!$omp           psi_axisym,                                                                    &
!$omp           wst, BigR, r0, T0, T0e, zj0, ps0, dTdx, dTdy, drhodx, drhody, dpsidx, dpsidy, dudx, dudy,  &
!$omp           dpdx, dpdy, phi, T0i, psi_as_coord,                                            &
!$omp           source_pellet, source_volume, eq_zne, eq_zTe, vpar0, BB2,                      &
!$omp           heat_source, heat_source_i, heat_source_e, particle_source, current_source, rotation_source, &
!$omp           dn_dpsi,dn_dz,dn_dpsi2,dn_dz2,dn_dpsi_dz,dn_dpsi3,dn_dpsi_dz2, dn_dpsi2_dz,    &
!$omp           dT_dpsi,dT_dz,dT_dpsi2,dT_dz2,dT_dpsi_dz,dT_dpsi3,dT_dpsi_dz2, dT_dpsi2_dz,    &
!$omp           hel1, vpar_x, vpar_y, ps0_s, ps0_t, u0_s, u0_t, p0_s, p0_t, vpar_s, vpar_t,    &
!$omp           u0_x, u0_y, T0e_corr, T0i_corr,                                                &
!$omp           thm_wk, mag_wk, eta_T, vpar_disp, fric_disp, p0_p, T0_corr, r0_corr, u0_p,     &
!$omp           AR0, AR0_p, AR0_s, AR0_t, AR0_sp, AR0_tp, AR0_Rp, AZ0, AZ0_p, AZ0_s, AZ0_t, AZ0_sp, AZ0_tp, AZ0_Zp, A30, &
!$omp           A30_p, A30_s, A30_t, A30_ss, A30_tt, A30_st, A30_R, A30_RR, A30_ZZ, BR_Z, BZ_R,&
!$omp           eta_T_ohm,                                                              &

#if (defined WITH_Neutrals) || (defined WITH_Impurities)
!$omp           rn0, rn0_corr, i_imp, frad_bg, Lrad_imp, Te_corr_eV, Te_eV, ne_SI, Ti_eV,      &
!$omp           spi_R_tmp, spi_Z_tmp, spi_phi_tmp, ng_radius_tmp,                              &
!$omp           spi_psi_tmp, spi_grad_psi_tmp,                                                 &
!$omp           n_spi_tmp, source_tmp, ns_shape,                                               &
#endif
#ifdef WITH_Impurities
!$omp           source_bg, source_imp,                                                         &
!$omp           m_i_over_m_imp, m_imp, Z_imp, dZ_imp_dT,                                       &
!$omp           ne_JOREK, P_imp, Lrad, E_ion, E_ion_bg, ion_i,                                 &
!$omp           ion_k, Z_eff, Z_eff_imp, eta_coef, Ti_corr_eV,                                 &
#endif
#if (defined WITH_Impurities) && (defined WITH_TiTe)
!$omp           alpha_i, alpha_e, nu_e_imp, nu_e_bg, lambda_e_imp, lambda_e_bg, dTi_e, dTe_i,  &
!$omp           dalpha_e_dT,                                                                   &
#endif
#if (defined WITH_Impurities) && (!defined WITH_TiTe) 
!$omp           alpha_imp, beta_imp, dalpha_imp_dT,                                            &
#endif
#if (defined WITH_Neutrals) && (!defined WITH_Impurities)
!$omp           Sion_T, dSion_dT, Srec_T, dSrec_dT, ksiion, source_neutral,                    &
!$omp           LradDrays_T, LradDcont_T, dLradDrays_dT, dLradDcont_dT,          &
!$omp           Arad_bg, Brad_bg, Crad_bg,                                                     &
!$omp           coef_prad_si,                                                                  &
#endif
!$omp           omp_nthreads,omp_tid)


#ifdef OPENMP
omp_nthreads = omp_get_num_threads()
omp_tid      = omp_get_thread_num()
#else
omp_nthreads = 1
omp_tid      = 0
#endif

!$omp do reduction(+:local_pellet_particles, local_plasma_particles, local_pellet_volume,     &
#if (defined WITH_Neutrals) || (defined WITH_Impurities)
!$omp                local_n_particles_inj,  local_n_particles,                               &
!$omp                local_radiation, local_radiation_phi, local_E_ion, local_P_ei, local_P_ion, &
!$omp                local_source_volume, local_source_volume_drift,                          &
#endif
!$omp                D_int, D_ext, P_int, H_int, S_int, H_ext, S_ext, P_ext, C_intern, C_ext, &
!$omp                P_e_int, P_i_int, P_e_ext, P_i_ext, P_e_tot, P_i_tot,                    &
!$omp                VP_int, VP_ext, VP_tot, VK_tot, VK_int, VK_ext, VM_ext,                  &
!$omp                VM_int, VM_tot, Vol, P_tot, D_tot,J2_tot, J2_int, J2_ext,                &
!$omp                heli_tot, mag_wk_tot, vpar_disp_tot, thm_wk_tot, area1, mag_src_tot,     &
!$omp                fric_disp_tot, R2curr_tmp, Zcurr_tmp)

do ife = ife_min, ife_max

  element = element_list%element(ife)

  do iv = 1, n_vertex_max
    inode     = element%vertex(iv)
    nodes(iv) = node_list%node(inode)
  enddo

  x_g(:,:)    = 0.d0; x_s(:,:)    = 0.d0; x_t(:,:)    = 0.d0; x_ss(:,:)    = 0.d0; x_tt(:,:)    = 0.d0; x_st(:,:)    = 0.d0;
  y_g(:,:)    = 0.d0; y_s(:,:)    = 0.d0; y_t(:,:)    = 0.d0; y_ss(:,:)    = 0.d0; y_tt(:,:)    = 0.d0; y_st(:,:)    = 0.d0;
  psi_axisym(:,:) = 0.d0

  do i=1,n_vertex_max
    do j=1,n_order+1

      do ms=1, n_gauss
        do mt=1, n_gauss

          x_g(ms,mt) = x_g(ms,mt) + nodes(i)%x(1,j,1) * element%size(i,j) * H(i,j,ms,mt)
          y_g(ms,mt) = y_g(ms,mt) + nodes(i)%x(1,j,2) * element%size(i,j) * H(i,j,ms,mt)

          x_s(ms,mt) = x_s(ms,mt) + nodes(i)%x(1,j,1) * element%size(i,j) * H_s(i,j,ms,mt)
          x_t(ms,mt) = x_t(ms,mt) + nodes(i)%x(1,j,1) * element%size(i,j) * H_t(i,j,ms,mt)
          x_ss(ms,mt) = x_ss(ms,mt) + nodes(i)%x(1,j,1) * element%size(i,j) * H_ss(i,j,ms,mt)
          x_tt(ms,mt) = x_tt(ms,mt) + nodes(i)%x(1,j,1) * element%size(i,j) * H_tt(i,j,ms,mt)
          x_st(ms,mt) = x_st(ms,mt) + nodes(i)%x(1,j,1) * element%size(i,j) * H_st(i,j,ms,mt)
          
          y_s(ms,mt) = y_s(ms,mt) + nodes(i)%x(1,j,2) * element%size(i,j) * H_s(i,j,ms,mt)
          y_t(ms,mt) = y_t(ms,mt) + nodes(i)%x(1,j,2) * element%size(i,j) * H_t(i,j,ms,mt)
          y_ss(ms,mt) = y_ss(ms,mt) + nodes(i)%x(1,j,2) * element%size(i,j) * H_ss(i,j,ms,mt)
          y_tt(ms,mt) = y_tt(ms,mt) + nodes(i)%x(1,j,2) * element%size(i,j) * H_tt(i,j,ms,mt)
          y_st(ms,mt) = y_st(ms,mt) + nodes(i)%x(1,j,2) * element%size(i,j) * H_st(i,j,ms,mt)

#ifdef fullmhd
          ! --- Equilibrium psi (n=0 only)
          psi_axisym(ms,mt) = psi_axisym(ms,mt) + nodes(i)%values(1,j,var_A3) * element%size(i,j) * H(i,j,ms,mt)
#endif

        enddo
      enddo
    enddo
  enddo

  eq_g(:,:,:,:) = 0.d0; eq_s(:,:,:,:) = 0.d0; eq_t(:,:,:,:) = 0.d0; eq_p(:,:,:,:) = 0.d0; eq_ss(:,:,:,:) = 0.d0; eq_tt(:,:,:,:) = 0.d0; eq_st(:,:,:,:) = 0.d0; 
  eq_sp(:,:,:,:) = 0.d0; eq_tp(:,:,:,:) = 0.d0;

  do i=1,n_vertex_max
    do j=1,n_order+1

      do mp=1,n_plane
        do ms=1, n_gauss
          do mt=1, n_gauss

            do k=1,n_var
              do in=1,n_tor
                eq_g(mp,k,ms,mt) = eq_g(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H(i,j,ms,mt)  * HZ(in,mp)
                eq_s(mp,k,ms,mt) = eq_s(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_s(i,j,ms,mt)* HZ(in,mp)
                eq_t(mp,k,ms,mt) = eq_t(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_t(i,j,ms,mt)* HZ(in,mp)
                eq_p(mp,k,ms,mt) = eq_p(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H(i,j,ms,mt)  * HZ_p(in,mp)
                eq_sp(mp,k,ms,mt) = eq_sp(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_s(i,j,ms,mt)* HZ_p(in,mp)
                eq_tp(mp,k,ms,mt) = eq_tp(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_t(i,j,ms,mt)* HZ_p(in,mp)
                eq_ss(mp,k,ms,mt) = eq_ss(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_ss(i,j,ms,mt)* HZ(in,mp)
                eq_tt(mp,k,ms,mt) = eq_tt(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_tt(i,j,ms,mt)* HZ(in,mp)
                eq_st(mp,k,ms,mt) = eq_st(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_st(i,j,ms,mt)* HZ(in,mp)
              enddo
            enddo

	  enddo
        enddo
      enddo

    enddo
  enddo
  
  ! --- Determine smallest and largest values of the variables in the whole domain (on Gauss points and toroidal integration surfaces)
  !$omp critical
  do k = 1, n_var
    varmin(k) = min( varmin(k), minval(eq_g(:,k,:,:)) )
    varmax(k) = max( varmax(k), maxval(eq_g(:,k,:,:)) )
  end do
  !$omp end critical
  
  do ms=1, n_gauss
    do mt=1, n_gauss
      call density(xpoint, xcase, y_g(ms,mt), Z_xpoint, eq_g(1,var_psi,ms,mt),psi_axis,psi_bnd,eq_zne(ms,mt), &
                   dn_dpsi,dn_dz,dn_dpsi2,dn_dz2,dn_dpsi_dz,dn_dpsi3,dn_dpsi_dz2, dn_dpsi2_dz)

#ifdef WITH_TiTe
      call temperature_e(xpoint, xcase, y_g(ms,mt), Z_xpoint, eq_g(1,1,ms,mt),psi_axis,psi_bnd,eq_zTe(ms,mt), &
                       dT_dpsi,dT_dz,dT_dpsi2,dT_dz2,dT_dpsi_dz,dT_dpsi3,dT_dpsi_dz2, dT_dpsi2_dz)
#else
      call temperature(xpoint, xcase, y_g(ms,mt), Z_xpoint, eq_g(1,1,ms,mt),psi_axis,psi_bnd,eq_zTe(ms,mt),   &
                       dT_dpsi,dT_dz,dT_dpsi2,dT_dz2,dT_dpsi_dz,dT_dpsi3,dT_dpsi_dz2, dT_dpsi2_dz)
      eq_zTe(ms,mt) = eq_zTe(ms,mt) / 2.d0	! electron temperature
#endif
    enddo
  enddo

  !--------------------------------------------------- sum over the Gaussian integration points
  do mp=1,n_plane

    phi       = 2.d0*PI*float(mp-1)/float(n_plane) / float(n_period)

    do ms=1, n_gauss
      do mt=1, n_gauss

        wst  = wgauss_copy(ms)*wgauss_copy(mt)

        xjac = x_s(ms,mt)*y_t(ms,mt) - x_t(ms,mt)*y_s(ms,mt)
        xjac_R = ( x_ss(ms,mt) * y_t(ms,mt)**2 - 2*x_st(ms,mt) * y_s(ms,mt)*y_t(ms,mt) + x_tt(ms,mt) * y_s(ms,mt)**2 + x_s(ms,mt) * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )   &
                 + x_t(ms,mt) * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) ) ) / xjac
        xjac_Z = ( y_ss(ms,mt) * x_t(ms,mt)**2 - 2*y_st(ms,mt) * x_s(ms,mt)*x_t(ms,mt) + y_tt(ms,mt) * x_s(ms,mt)**2 + y_s(ms,mt) * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )   &
                 + y_t(ms,mt) * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) ) ) / xjac
        BigR = x_g(ms,mt)

        r0     = eq_g(mp,var_rho,ms,mt)
        r0_corr = corr_neg_dens1(r0)
#ifdef WITH_TiTe
        T0i    = eq_g(mp,var_Ti,ms,mt)
        T0e    = eq_g(mp,var_Te,ms,mt)
        T0e_corr = corr_neg_temp1(T0e)
        T0i_corr = corr_neg_temp1(T0i)
#else
        T0     = eq_g(mp,var_T,ms,mt)
        T0i    = eq_g(mp,var_T,ms,mt) /2.d0
        T0e    = eq_g(mp,var_T,ms,mt) /2.d0
        T0_corr = corr_neg_temp1(T0)
        T0e_corr = corr_neg_temp1(T0e)
        T0i_corr = corr_neg_temp1(T0i)
#endif
        zj0    = eq_g(mp,var_zj,ms,mt)
        ps0    = eq_g(mp,var_psi,ms,mt)
        ps0_s  = eq_s(mp,var_psi,ms,mt) 
        ps0_t  = eq_t(mp,var_psi,ms,mt)
        u0_s   = eq_s(mp,var_u,ms,mt) 
        u0_t   = eq_t(mp,var_u,ms,mt)
        u0_p   = eq_p(mp,var_u,ms,mt)
        u0_x   = (   y_t(ms,mt) * eq_s(mp,var_u,ms,mt) - y_s(ms,mt) * eq_t(mp,var_u,ms,mt) ) / xjac
        u0_y   = ( - x_t(ms,mt) * eq_s(mp,var_u,ms,mt) + x_s(ms,mt) * eq_t(mp,var_u,ms,mt) ) / xjac

        vpar0   = eq_g(mp,var_Vpar,ms,mt)
        vpar_s  = eq_s(mp,var_Vpar,ms,mt)
        vpar_t  = eq_t(mp,var_Vpar,ms,mt)

#if (defined WITH_Neutrals) || (defined WITH_Impurities)
        rn0    = eq_g(mp,var_rhon,ms,mt)
        rn0_corr = corr_neg_dens(rn0, (/ 0.d-5, 1.d-5 /)) ! Correction for negative rn0
#endif
      
#ifdef fullmhd
        AR0      = eq_g(mp,var_AR,ms,mt)
        AR0_p    = eq_p(mp,var_AR,ms,mt)
        AR0_s    = eq_s(mp,var_AR,ms,mt)
        AR0_t    = eq_t(mp,var_AR,ms,mt)
        AR0_sp   = eq_sp(mp,var_AR,ms,mt)
        AR0_tp   = eq_tp(mp,var_AR,ms,mt)
        AR0_Rp   = (   y_t(ms,mt) * AR0_sp  - y_s(ms,mt) * AR0_tp ) / xjac
        AZ0      = eq_g(mp,var_AZ,ms,mt)
        AZ0_p    = eq_p(mp,var_AZ,ms,mt)
        AZ0_s    = eq_s(mp,var_AZ,ms,mt)
        AZ0_t    = eq_t(mp,var_AZ,ms,mt)
        AZ0_sp   = eq_sp(mp,var_AZ,ms,mt)
        AZ0_tp   = eq_tp(mp,var_AZ,ms,mt)
        AZ0_Zp   = ( - x_t(ms,mt) * AZ0_sp  + x_s(ms,mt) * AZ0_tp ) / xjac
        A30      = eq_g(mp,var_A3,ms,mt)
        A30_p    = eq_p(mp,var_A3,ms,mt)
        A30_s    = eq_s(mp,var_A3,ms,mt)
        A30_t    = eq_t(mp,var_A3,ms,mt)
        A30_ss   = eq_ss(mp,var_A3,ms,mt)
        A30_tt   = eq_tt(mp,var_A3,ms,mt)
        A30_st   = eq_st(mp,var_A3,ms,mt)
        A30_R = (   y_t(ms,mt) * A30_s  - y_s(ms,mt) * A30_t ) / xjac
        A30_RR   = (A30_ss * y_t(ms,mt)**2 - 2.d0*A30_st * y_s(ms,mt)*y_t(ms,mt) + A30_tt * y_s(ms,mt)**2  &
                    + A30_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                           &
                    + A30_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) ) )       / xjac**2         &
                    - xjac_R * (A30_s * y_t(ms,mt) - A30_t * y_s(ms,mt))     / xjac**2
        A30_ZZ   = (A30_ss * x_t(ms,mt)**2 - 2.d0*A30_st * x_s(ms,mt)*x_t(ms,mt) + A30_tt * x_s(ms,mt)**2  &
                    + A30_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                           &
                    + A30_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) ) )       / xjac**2         &
                    - xjac_Z * (- A30_s * x_t(ms,mt) + A30_t * x_s(ms,mt) )  / xjac**2
        BR_Z = ( A30_ZZ - AZ0_Zp )/ BigR
        BZ_R = -1/BigR**2 * ( AR0_p - A30_R ) + ( AR0_Rp - A30_RR )/ BigR
        zj0 = - (BZ_R - BR_Z) * BigR
        psi_as_coord = psi_axisym(ms,mt)
#else
        psi_as_coord = ps0
#endif
#ifdef WITH_TiTe
        eta_T         = resistivity(eta, T0e_corr, T_max_eta, Te_0)  
        eta_T_ohm     = resistivity(eta_ohmic, T0e_corr, T_max_eta_ohm, Te_0)
#else
        eta_T         = resistivity(eta, T0_corr, T_max_eta, T_0)  
        eta_T_ohm     = resistivity(eta_ohmic, T0_corr, T_max_eta_ohm, T_0)
#endif
        ! This is currently broken for two temperature models !
        ! Some of these do not seem to be doing anything so I'm commenting them off !
        !dTdx   = (   y_t(ms,mt) * eq_s(mp,var_T,ms,mt) - y_s(ms,mt) * eq_t(mp,var_T,ms,mt) ) / xjac
        !dTdy   = ( - x_t(ms,mt) * eq_s(mp,var_T,ms,mt) + x_s(ms,mt) * eq_t(mp,var_T,ms,mt) ) / xjac
        !drhodx = (   y_t(ms,mt) * eq_s(mp,var_rho,ms,mt) - y_s(ms,mt) * eq_t(mp,var_rho,ms,mt) ) / xjac
        !drhody = ( - x_t(ms,mt) * eq_s(mp,var_rho,ms,mt) + x_s(ms,mt) * eq_t(mp,var_rho,ms,mt) ) / xjac

        dpsidx = (   y_t(ms,mt) * eq_s(mp,var_psi,ms,mt) - y_s(ms,mt) * eq_t(mp,var_psi,ms,mt) ) / xjac
        dpsidy = ( - x_t(ms,mt) * eq_s(mp,var_psi,ms,mt) + x_s(ms,mt) * eq_t(mp,var_psi,ms,mt) ) / xjac

        dudx = (   y_t(ms,mt) * eq_s(mp,var_u,ms,mt) - y_s(ms,mt) * eq_t(mp,var_u,ms,mt) ) / xjac
        dudy = ( - x_t(ms,mt) * eq_s(mp,var_u,ms,mt) + x_s(ms,mt) * eq_t(mp,var_u,ms,mt) ) / xjac

        vpar_x = (   y_t(ms,mt) * vpar_s - y_s(ms,mt) * vpar_t ) / xjac
        vpar_y = ( - x_t(ms,mt) * vpar_s + x_s(ms,mt) * vpar_t ) / xjac

        BB2 = (F0*F0 + dpsidx*dpsidx + dpsidy*dpsidy) / BigR**2

        !dPdx = r0 * dTdx + T0 * drhodx
        !dPdy = r0 * dTdy + T0 * drhody

#ifdef WITH_TiTe
        call sources(xpoint, xcase, y_g(ms,mt), Z_xpoint, ps0, psi_axis, psi_bnd, &
                     particle_source,heat_source_i,heat_source_e)
        heat_source = heat_source_i + heat_source_e
#else
        call sources(xpoint, xcase, y_g(ms,mt), Z_xpoint, psi_as_coord, psi_axis, psi_bnd, &
                     particle_source,heat_source)
#endif
        if (keep_current_prof) then
          call current(xpoint, xcase, x_g(ms,mt),y_g(ms,mt), Z_xpoint, psi_as_coord,&
                       psi_axis,psi_bnd,current_source)
        else
          current_source = 0.d0
        endif

!-------------------------------------------
! --- Radiation and ionization power
! ------------------------------------------
#if ( (defined WITH_Neutrals) && (! defined WITH_Impurities) )
        ! --- Get ionization, recombination and radiation coefficients for Deuterium 
#ifdef WITH_TiTe
        call atomic_coeff_deuterium  (   T0e, Sion_T, dSion_dT, Srec_T, dSrec_dT,        &
                                            LradDcont_T, dLradDcont_dT, LradDrays_T, dLradDrays_dT, r0, rn0, .true. ) 
#else
        call atomic_coeff_deuterium(0.5d0*T0, Sion_T, dSion_dT, Srec_T, dSrec_dT,        &
                                            LradDcont_T, dLradDcont_dT, LradDrays_T, dLradDrays_dT, r0, rn0, .true. ) 
#endif
     
        ! Get coefficient:  Prad,SI = coef_prad_si * Prad,jorek
        coef_prad_si = 1./((GAMMA-1)*MU_ZERO*(MU_ZERO*central_mass*MASS_PROTON*central_density*1.d20)**0.5) 
      
        ksiion = central_density * 1.d20 * ksi_ion   ! Normalisation of the ionization energy cost for Deuterium
      
        ! --- Radiation from background impurity
        ne_SI = r0_corr * 1.d20 * central_density !electron density (SI)
        Te_corr_eV = T0e_corr/(EL_CHG*MU_ZERO*central_density*1.d20) ! Te in eV
        Te_eV = T0e/(EL_CHG*MU_ZERO*central_density*1.d20) ! Te in eV, uncorrected

        if (use_imp_adas) then  ! use open adas by default  
          ! Use radiation coefficients from ADAS
          frad_bg = 0. 
          do i_imp = 1, n_adas     
            if (ne_SI > ne_SI_min .and. Te_eV > Te_eV_min .and. nimp_bg(i_imp) > 0) then
              Lrad_imp = 0.0
              call radiation_function_linear(imp_adas(i_imp),imp_cor(i_imp),log10(ne_SI),    & 
                                           log10(Te_corr_eV*EL_CHG/K_BOLTZ),.false.,Lrad_imp)           
            else     
              Lrad_imp = 0.
            end if
            frad_bg = frad_bg + nimp_bg(i_imp) * Lrad_imp
          end do

          local_radiation_phi(mp) = local_radiation_phi(mp) + ( (r0_corr * rn0_corr  * LradDrays_T    &
                                     + r0_corr ** 2 * LradDcont_T) * coef_prad_si                     & 
                                     + ne_SI * frad_bg) * bigR * xjac * wst * delta_phi  
          local_radiation         = local_radiation + ( (r0_corr * rn0_corr  * LradDrays_T            &
                                     + r0_corr ** 2 * LradDcont_T) * coef_prad_si                     & 
                                     + ne_SI * frad_bg) * bigR * xjac * wst * delta_phi 
          local_P_ion             = local_P_ion + ksiion * r0_corr * rn0_corr * Sion_T * coef_prad_si &
                                   * bigR * xjac * wst * delta_phi
        else
          if ( trim(imp_type(1)) == 'Ar') then ! Hard-coded fitting exists for argon
            Arad_bg = 2.4d-31 
            Brad_bg = 20.
            Crad_bg = 0.8
            frad_bg = (2./3.)*(1./(central_mass*MASS_PROTON))*((MU_ZERO*central_mass*MASS_PROTON*central_density*1.d20)**(1.5d0))                &
                            *nimp_bg(1)*Arad_bg*exp(-((log(Te_corr_eV)-log(Brad_bg))**2.)/Crad_bg**2.)
                    
            local_radiation_phi(mp) = local_radiation_phi(mp) + (r0_corr * rn0_corr  * LradDrays_T &
                                       + r0_corr ** 2 * LradDcont_T + r0_corr * frad_bg) * coef_prad_si & 
                                       * bigR * xjac * wst * delta_phi  
            local_radiation         = local_radiation + (r0_corr * rn0_corr  * LradDrays_T &
                                       + r0_corr ** 2 * LradDcont_T + r0_corr * frad_bg) * coef_prad_si & 
                                       * bigR * xjac * wst * delta_phi 
            local_P_ion             = local_P_ion + ksiion * r0_corr * rn0_corr * Sion_T * coef_prad_si &
                                       * bigR * xjac * wst * delta_phi
          else
            write(*,*) "WARNING: hard-coded fitting doesn't exist for  ", trim(imp_type(1)), ", use open adas instead!"
            stop
          end if
        end if

#endif

#ifdef WITH_Impurities
        !-------------------------------------------
        ! Atomic physics parameters for Impurities
        !-------------------------------------------

        select case ( trim(imp_type(index_main_imp)) )
          case('D2')
            m_i_over_m_imp = central_mass/2.  ! Deuterium mass = 2 u
            m_imp          = 2.
          case('Ar')
            m_i_over_m_imp = central_mass/40. ! Argon mass = 40 u
            m_imp          = 40.
          case('Ne')
            m_i_over_m_imp = central_mass/20. ! Neon mass = 20 u
            m_imp          = 20.
          case default
            write(*,*) '!! Gas type "', trim(imp_type(index_main_imp)), '" unknown (in mod_injection_source.f90) !!'
            write(*,*) '=> We assume the gas is D2.'
            m_i_over_m_imp = central_mass/2.
            m_imp          = 2.
        end select

        ! Te in eV:
        Te_corr_eV = T0e_corr/(EL_CHG*MU_ZERO*central_density*1.d20)
        Te_eV = T0e/(EL_CHG*MU_ZERO*central_density*1.d20)
   
        ! Ti in eV:
        Ti_corr_eV = T0i_corr/(EL_CHG*MU_ZERO*central_density*1.d20)
        Ti_eV = T0i/(EL_CHG*MU_ZERO*central_density*1.d20)

        if (allocated(P_imp)) deallocate(P_imp)
        allocate(P_imp(0:imp_adas(index_main_imp)%n_Z))
        call imp_cor(index_main_imp)%interp_linear(density=20.,temperature=log10(Te_corr_eV*EL_CHG/K_BOLTZ),&
                                      p_out=P_imp,z_avg=Z_imp, z_avg_Te=dZ_imp_dT)

        if (allocated(imp_adas(index_main_imp)%ionisation_energy)) then
   
          ! Calculate the ionization potential energy and its derivative wrt. temperature
          E_ion     = 0.
          E_ion_bg  = 13.6
          do ion_i=1, imp_adas(index_main_imp)%n_Z
            do ion_k=1, ion_i
              E_ion     = E_ion + P_imp(ion_i)*imp_adas(index_main_imp)%ionisation_energy(ion_k)
            end do
          end do
          ! Convert from eV to SI unit
          E_ion     = E_ion * EL_CHG
          E_ion_bg  = E_ion_bg * EL_CHG
        else
          E_ion     = 0.
          E_ion_bg  = 0.
        end if

#ifdef WITH_TiTe
        alpha_i       = m_i_over_m_imp - 1.
        alpha_e       = m_i_over_m_imp*Z_imp - 1.
        dalpha_e_dT   = m_i_over_m_imp*dZ_imp_dT
  
        ne_SI        = (r0_corr + alpha_e * rn0_corr) * 1.d20 * central_density ! electron density (SI)
        ne_JOREK     = r0_corr + alpha_e * rn0_corr ! Electron density in JOREK unit
        ne_JOREK     = corr_neg_dens(ne_JOREK,(/1.d-1,1.d-1/),1.d-3) ! Correction for negative electron density
                                                            ! Too small rho_1 will cause a problem
#else /* WITH_TiTe */
        alpha_imp    = 0.5*m_i_over_m_imp*(Z_imp+1.) - 1.
        dalpha_imp_dT= 0.5*m_i_over_m_imp*dZ_imp_dT
        beta_imp     = m_i_over_m_imp*Z_imp - 1.
        ne_SI        = (r0_corr + beta_imp * rn0_corr) * 1.d20 * central_density !electron density (SI)
        ne_JOREK     = r0_corr + beta_imp * rn0_corr ! Electron density in JOREK unit
        ne_JOREK     = corr_neg_dens(ne_JOREK,(/1.d-1,1.d-1/),1.d-3) ! Correction for negative electron density
                                                               ! Too small rho_1 will cause a problem
#endif /* WITH_TiTe */

        ! Calculate the effective charge of all species
        Z_eff        = 0.
        Z_eff_imp    = 0.
   
        ! First get the value of Z_eff
        Z_eff        = r0_corr - rn0_corr
        do ion_i=1, imp_adas(index_main_imp)%n_Z
          Z_eff      = Z_eff + m_i_over_m_imp * rn0_corr * P_imp(ion_i) * real(ion_i,8)**2
          Z_eff_imp  = Z_eff_imp + P_imp(ion_i) * real(ion_i,8)**2 ! The summation of normalized nZ**2 for impurity
        end do
        Z_eff        = Z_eff / ne_JOREK
  
        if (Z_eff < 1) Z_eff = 1.
   
        ! This is to represent the dependence on Z_eff in resistivity
        if ( eta_T_dependent ) then
          eta_coef     = Z_eff*(1.+1.198*Z_eff+0.222*Z_eff**2)/(1.+2.966*Z_eff+0.753*Z_eff**2)
          eta_coef     = eta_coef / ((1.+1.198+0.222)/(1.+2.966+0.753))
          eta_T        = eta_T * eta_coef
          eta_T_ohm    = eta_T_ohm * eta_coef
        endif

        if (ne_SI > ne_SI_min .and. Te_eV > Te_eV_min .and. rn0 > rn0_min) then
          Lrad = 0.0
          !call radiation_function(imp_adas(index_main_imp),imp_cor(index_main_imp),log10(ne_SI),log10(Te_corr_eV*EL_CHG/K_BOLTZ),Lrad)
          call radiation_function_linear(imp_adas(index_main_imp),imp_cor(index_main_imp),log10(ne_SI),log10(Te_corr_eV*EL_CHG/K_BOLTZ),.false.,Lrad)
        else
          Lrad = 0.
        end if

        Lrad = Lrad * m_i_over_m_imp
        E_ion = E_ion * m_i_over_m_imp

        frad_bg = 0. 
        do i_imp = 1, n_adas     
          if (i_imp == index_main_imp) cycle
          if (ne_SI > ne_SI_min .and. Te_eV > Te_eV_min .and. nimp_bg(i_imp) > 0) then
            Lrad_imp = 0.0
            call radiation_function_linear(imp_adas(i_imp),imp_cor(i_imp),log10(ne_SI),    & 
                                         log10(Te_eV*EL_CHG/K_BOLTZ),.false.,Lrad_imp)           
          else     
            Lrad_imp = 0.
          end if
          frad_bg = frad_bg + nimp_bg(i_imp) * Lrad_imp
        end do

        local_radiation_phi(mp) = local_radiation_phi(mp) &
                                  + ne_SI * (rn0_corr * central_density * 1.d20 * Lrad + frad_bg)&
                                  * bigR * xjac * wst * delta_phi        
        local_radiation = local_radiation &
                          + ne_SI * (rn0_corr * central_density * 1.d20 * Lrad + frad_bg)&
                          * bigR * xjac * wst * delta_phi 
        local_E_ion     = local_E_ion + rn0 * central_density * 1.d20 * E_ion             &
                          * bigR * xjac * wst * delta_phi
        local_E_ion     = local_E_ion + (r0 - rn0) * central_density * 1.d20 * E_ion_bg   &
                          * bigR * xjac * wst * delta_phi
#ifdef WITH_TiTe
       !--------------------------------------------------------
       ! --- Ion-electron energy transfer
       !--------------------------------------------------------
        if (Te_corr_eV < 10.*Z_imp**2) then
          lambda_e_imp = 23. - log((ne_SI*1.d-6)**0.5*Z_imp*Te_corr_eV**(-1.5))
        else
          lambda_e_imp = 24. - log((ne_SI*1.d-6)**0.5*Te_corr_eV**(-1.0))
        endif
        if (Te_corr_eV < 10.) then
          lambda_e_bg  = 23. - log((ne_SI*1.d-6)**0.5*Te_corr_eV**(-1.5)) ! Assuming bg_charge is 1! 
        else
          lambda_e_bg  = 24. - log((ne_SI*1.d-6)**0.5*Te_corr_eV**(-1.0))
        endif
        nu_e_imp     = 1.8d-19*(1.d6*MASS_ELECTRON*MASS_PROTON*m_imp) ** 0.5&
                       * Z_eff_imp * (1.d14*central_density*rn0_corr*m_i_over_m_imp) * lambda_e_imp &
                       / (1.d3*(MASS_ELECTRON*T0i_corr+T0e_corr*MASS_PROTON*m_imp)&
                       / (EL_CHG * MU_ZERO * central_density * 1.d20)) ** 1.5
        nu_e_bg      = 1.8d-19*(1.d6*MASS_ELECTRON*MASS_PROTON*central_mass) ** 0.5&
                       * (1.d14*central_density*(r0_corr-rn0_corr)) * lambda_e_bg &
                       / (1.d3*(MASS_ELECTRON*T0i_corr+T0e_corr*MASS_PROTON*central_mass)&
                       / (EL_CHG * MU_ZERO * central_density * 1.d20)) ** 1.5 ! Assuming bg_charge is 1!
    
        if (nu_e_imp < 0.) nu_e_imp = 0.
        if (nu_e_bg < 0.)  nu_e_bg  = 0.
    
        dTe_i    = (nu_e_imp + nu_e_bg) * (Ti_corr_eV - Te_corr_eV) * EL_CHG ! To SI unit
        local_P_ei = local_P_ei + ne_SI * dTe_i * bigR * xjac * wst * delta_phi 
#endif /* WITH_TiTe */
#endif /* WITH_Impurities */
#ifdef WITH_Impurities
        D_tot  = D_tot  + (r0-rn0) * xjac * BigR * wst * delta_phi 
#ifdef WITH_TiTe
        P_e_tot = P_e_tot + (r0+alpha_e*rn0) * T0e * xjac * BigR * wst * delta_phi
        P_i_tot = P_i_tot + (r0+alpha_i*rn0) * T0i * xjac * BigR * wst * delta_phi
        P_tot   = P_e_tot + P_i_tot

        p0_s   = (r0+alpha_i*rn0)*eq_s(mp,var_Ti,ms,mt) &
                 + T0i * (eq_s(mp,var_rho,ms,mt)+alpha_i*eq_s(mp,var_rhon,ms,mt))&
                 + (r0+alpha_e*rn0+dalpha_e_dT*rn0*T0e)*eq_s(mp,var_Te,ms,mt)&
                 + T0e * (eq_s(mp,var_rho,ms,mt)+alpha_e*eq_s(mp,var_rhon,ms,mt))
        p0_t   = (r0+alpha_i*rn0)*eq_t(mp,var_Ti,ms,mt) &
                 + T0i * (eq_t(mp,var_rho,ms,mt)+alpha_i*eq_t(mp,var_rhon,ms,mt))&
                 + (r0+alpha_e*rn0+dalpha_e_dT*rn0*T0e)*eq_t(mp,var_Te,ms,mt)&
                 + T0e * (eq_t(mp,var_rho,ms,mt)+alpha_e*eq_t(mp,var_rhon,ms,mt))
        p0_p   = (r0+alpha_i*rn0)*eq_p(mp,var_Ti,ms,mt) &
                 + T0i * (eq_p(mp,var_rho,ms,mt)+alpha_i*eq_p(mp,var_rhon,ms,mt))&
                 + (r0+alpha_e*rn0+dalpha_e_dT*rn0*T0e)*eq_p(mp,var_Te,ms,mt)&
                 + T0e * (eq_p(mp,var_rho,ms,mt)+alpha_e*eq_p(mp,var_rhon,ms,mt))
#else /* WITH_TiTe */
        P_tot  = P_tot  + (r0+alpha_imp*rn0) * T0 * xjac * BigR * wst * delta_phi
        P_e_tot = P_tot / 2.
        P_i_tot = P_e_tot

        p0_s   = (r0+alpha_imp*rn0+dalpha_imp_dT*rn0*T0)*eq_s(mp,var_T,ms,mt) &
                 + T0 * (eq_s(mp,var_rho,ms,mt)+alpha_imp*eq_s(mp,var_rhon,ms,mt))
        p0_t   = (r0+alpha_imp*rn0+dalpha_imp_dT*rn0*T0)*eq_t(mp,var_T,ms,mt) &
                 + T0 * (eq_t(mp,var_rho,ms,mt)+alpha_imp*eq_t(mp,var_rhon,ms,mt))
        p0_p   = (r0+alpha_imp*rn0+dalpha_imp_dT*rn0*T0)*eq_p(mp,var_T,ms,mt) &
                 + T0 * (eq_p(mp,var_rho,ms,mt)+alpha_imp*eq_p(mp,var_rhon,ms,mt))
#endif /* WITH_TiTe */
#else /* WITH_Impurities */
        D_tot  = D_tot  + r0       * xjac * BigR * wst * delta_phi
#ifdef WITH_TiTe
        P_e_tot = P_e_tot + r0 * T0e * xjac * BigR * wst * delta_phi
        P_i_tot = P_i_tot + r0 * T0i * xjac * BigR * wst * delta_phi
        P_tot   = P_e_tot + P_i_tot

        p0_s   = r0*eq_s(mp,var_Te,ms,mt) + T0e * eq_s(mp,var_rho,ms,mt) &
                 +r0*eq_s(mp,var_Ti,ms,mt) + T0i * eq_s(mp,var_rho,ms,mt)
        p0_t   = r0*eq_t(mp,var_Te,ms,mt) + T0e * eq_t(mp,var_rho,ms,mt) &
                 +r0*eq_t(mp,var_Ti,ms,mt) + T0i * eq_t(mp,var_rho,ms,mt)
        p0_p   = r0*eq_p(mp,var_Te,ms,mt) + T0e * eq_p(mp,var_rho,ms,mt) &
                 +r0*eq_p(mp,var_Ti,ms,mt) + T0i * eq_p(mp,var_rho,ms,mt)
#else /* WITH_TiTe */
        P_tot  = P_tot  + r0 * T0 * xjac * BigR * wst * delta_phi
        P_e_tot = P_tot / 2.
        P_i_tot = P_e_tot

        p0_s   = r0*eq_s(mp,var_T,ms,mt) + T0 * eq_s(mp,var_rho,ms,mt) 
        p0_t   = r0*eq_t(mp,var_T,ms,mt) + T0 * eq_t(mp,var_rho,ms,mt) 
        p0_p   = r0*eq_p(mp,var_T,ms,mt) + T0 * eq_p(mp,var_rho,ms,mt) 
#endif /* WITH_TiTe */
#endif /* WITH_Impurities */

        thm_wk     = vpar0 * (p0_s*ps0_t - p0_t*ps0_s) + vpar0 * F0/BigR*p0_p*xjac 
        hel1       = F0* ( (ps0 - psi_off) - y_g(ms,mt)*dpsidy) / (BigR**2.d0)
        mag_wk     = - (ps0_s*u0_t - ps0_t*u0_s) / xjac * zj0 / BigR  &
                     + F0 * zj0 * u0_p / (BigR**2.d0)
        vpar_disp  = visco_par * (F0/BigR)**2.d0 * (vpar_x**2.d0+vpar_y**2.d0 ) 

        VP_tot = VP_tot + r0 * vpar0**2 * BB2 * xjac * BigR * wst * delta_phi
        VK_tot = VK_tot + r0 * (dudx**2 + dudy**2) * BigR**2 * xjac * BigR * wst * delta_phi
        VM_tot = VM_tot + (dpsidx**2+dpsidy**2)/BigR**2 * xjac * BigR * wst * delta_phi
        J2_tot = J2_tot + eta_T_ohm *(ZJ0/BigR)**2.d0 * xjac * BigR * wst * delta_phi

        mag_src_tot   = mag_src_tot + eta_T*ZJ0*current_source/(BigR**2) * xjac * BigR * wst * delta_phi

        heli_tot      = heli_tot   + hel1         * BigR * xjac * wst * delta_phi
        mag_wk_tot    = mag_wk_tot + mag_wk       * BigR * xjac * wst * delta_phi
        thm_wk_tot    = thm_wk_tot + thm_wk                     * wst * delta_phi
        vpar_disp_tot = vpar_disp_tot + vpar_disp * BigR * xjac * wst * delta_phi

        R2curr_tmp    = R2curr_tmp - x_g(ms,mt)**2.0 * zj0 /BigR * xjac * wst * delta_phi    
        Zcurr_tmp     = Zcurr_tmp  - y_g(ms,mt)      * zj0 /BigR * xjac * wst * delta_phi    

        if (use_pellet) then
          call pellet_source2(pellet_amplitude,pellet_R,pellet_Z,pellet_psi,pellet_phi, &
                              pellet_radius, pellet_delta_psi, pellet_sig, pellet_length, pellet_ellipse, pellet_theta, &
                              x_g(ms,mt),y_g(ms,mt), psi_as_coord, phi, r0_corr, T0_corr/2.d0, central_density, &
                              pellet_particles, pellet_density, pellet_volume, source_pellet, source_volume)

          local_pellet_particles = local_pellet_particles + source_pellet * bigR * xjac * wst * delta_phi
          local_plasma_particles = local_plasma_particles + r0            * bigR * xjac * wst * delta_phi
          local_pellet_volume    = local_pellet_volume    + source_volume * bigR * xjac * wst * delta_phi
        endif

#if (defined WITH_Neutrals) || (defined WITH_Impurities)
        if (using_spi) then

           if (JET_MGI .or. ASDEX_MGI) then
              write(*,*) "WARNING: Using SPI, disabling MGI settings"
              JET_MGI = .false.
              ASDEX_MGI = .false.
           end if

           do spi_i=1, n_spi_tot

              n_spi_tmp = 0
              do i_inj = 1, n_inj
                 n_spi_tmp = n_spi_tmp + n_spi(i_inj)
                 if (spi_i <= n_spi_tmp)  exit !< Determine the injection location index of the fragment
              end do

              if (t_now >= t_ns(i_inj)) then

                 source_tmp = 0.d0
                 ns_shape = 0.d0

                 spi_R_tmp   = pellets(spi_i)%spi_R
                 spi_Z_tmp   = pellets(spi_i)%spi_Z
                 spi_phi_tmp = pellets(spi_i)%spi_phi

                 spi_psi_tmp = pellets(spi_i)%spi_psi
                 spi_grad_psi_tmp = pellets(spi_i)%spi_grad_psi
                 
                 ng_radius_tmp   = pellets(spi_i)%spi_radius * ng_radius_ratio

                 if (ng_radius_tmp < ng_radius_min) then
                    ng_radius_tmp = ng_radius_min
                 end if

                 ! Compute the source shape
                 ns_shape = source_shape(x_g(ms,mt),y_g(ms,mt),phi, &
                      spi_R_tmp,spi_Z_tmp,spi_phi_tmp,              &
                      ng_radius_tmp,ns_deltaphi,                    &
                      ps0,spi_psi_tmp,spi_grad_psi_tmp,ns_deltaminrad)

                 local_source_volume(spi_i) = local_source_volume(spi_i) &
                      + ns_shape * bigR * xjac * wst * delta_phi

                 if (drift_distance /= 0) then ! Get the volume at the post-drift location (for normalization)
                   ns_shape = source_shape(x_g(ms,mt),y_g(ms,mt),phi,     &
                        spi_R_tmp+drift_distance,spi_Z_tmp,spi_phi_tmp,   &
                        ng_radius_tmp,ns_deltaphi,                        &
                        ps0,pellets(spi_i)%spi_psi_drift,                 &
                        pellets(spi_i)%spi_grad_psi_drift,                &
                        ns_deltaminrad)

                   local_source_volume_drift(spi_i) = local_source_volume_drift(spi_i) &
                        + ns_shape * bigR * xjac * wst * delta_phi
                 end if

              end if

           end do
        end if
#endif

#if ( (defined WITH_Neutrals) && (! defined WITH_Impurities) )
        !--- We calculate here the number of neutrals particles injected per second with n_particles_inj and the number of neutrals in the plasma

        source_neutral = 0.d0

        call total_neutral_source(x_g(ms,mt),y_g(ms,mt),phi,ps0,source_neutral)

        ! Neutral injection rate in particles/s
        local_n_particles_inj = local_n_particles_inj + 0.5d0 * central_density * 1.d20 * source_neutral * bigR *&
                                 xjac * wst * delta_phi / sqrt(MU_ZERO*central_mass*MASS_PROTON*central_density*1.d20)
        ! Total neutrals in particles
        local_n_particles     = local_n_particles     +  rn0 * central_density * 1.d20 * bigR * xjac * wst * delta_phi
        ! Frictional heat source
        fric_disp     =   0.5 * BigR**2 * (u0_x**2.0 + u0_y**2.0) * (r0_corr * rn0 * Sion_T)&
                        + 0.5 * vpar0**2 * BB2 * (r0_corr * rn0 * Sion_T)
        fric_disp_tot = fric_disp_tot + fric_disp * BigR * xjac * wst * delta_phi 
#endif
#ifdef WITH_Impurities
        !--- Calculate the neutral injection rate and the number of neutrals in the plasma

        source_imp = 0.d0
        source_bg  = 0.d0

        call total_imp_source(x_g(ms,mt),y_g(ms,mt),phi,ps0,source_bg,source_imp,m_i_over_m_imp,index_main_imp)

        ! Frictional heat source
        fric_disp     =   0.5 * BigR**2 * (u0_x**2.0 + u0_y**2.0) * (source_bg + source_imp)&
                        + 0.5 * vpar0**2 * BB2 * (source_bg + source_imp)
        fric_disp_tot = fric_disp_tot + fric_disp * BigR * xjac * wst * delta_phi 

        ! Neutral injection rate in particles/s
        local_n_particles_inj = local_n_particles_inj + 0.5d0 * central_density * 1.d20 * source_imp * m_i_over_m_imp * bigR &
	                         * xjac * wst * delta_phi / sqrt(MU_ZERO*central_mass*MASS_PROTON*central_density*1.d20)
        ! Total neutrals in particles
        local_n_particles     = local_n_particles + central_density * 1.d20 * rn0 * m_i_over_m_imp * bigR * xjac * wst * delta_phi
#endif
        if ( get_psi_n(psi_as_coord, y_g(ms,mt)) <= 1.d0 ) then   !inside LCFS
#ifdef WITH_Impurities
          D_int = D_int + (r0-rn0) * xjac * BigR * wst * delta_phi
#ifdef WITH_TiTe
          P_e_int = P_e_int + (r0+alpha_e*rn0) * T0e * xjac * BigR * wst * delta_phi
          P_i_int = P_i_int + (r0+alpha_i*rn0) * T0i * xjac * BigR * wst * delta_phi
          P_int   = P_e_int + P_i_int
#else /* WITH_TiTe */
          P_int = P_int + (r0+alpha_imp*rn0) * T0   * xjac * BigR * wst * delta_phi
          P_e_int = P_int / 2.
          P_i_int = P_e_int
#endif /* WITH_TiTe */
#else /* WITH_Impurities */
          D_int = D_int + r0        * xjac * BigR * wst * delta_phi
#ifdef WITH_TiTe
          P_e_int = P_e_int + r0 * T0e * xjac * BigR * wst * delta_phi
          P_i_int = P_i_int + r0 * T0i * xjac * BigR * wst * delta_phi
          P_int   = P_e_int + P_i_int
#else /* WITH_TiTe */
          P_int = P_int + r0 * T0   * xjac * BigR * wst * delta_phi
          P_e_int = P_int / 2.
          P_i_int = P_e_int
#endif /* WITH_TiTe */
#endif /* WITH_Impurities */
          C_intern = C_intern - zj0 /BigR * xjac *        wst * delta_phi    ! 2D integral
          area1    = area1    +  xjac * wst * delta_phi         
          Vol   = Vol   +             xjac * BigR * wst * delta_phi
          H_int = H_int + heat_source     * xjac * BigR * wst * delta_phi
          S_int = S_int + particle_source * xjac * BigR * wst * delta_phi
          VP_int = VP_int + r0 * vpar0**2 * BB2 * xjac * BigR * wst * delta_phi
          VK_int = VK_int + r0 * (dudx**2 + dudy**2) * BigR**2 * xjac * BigR * wst * delta_phi
          VM_int = VM_int + (dpsidx**2+dpsidy**2)/BigR**2 * xjac * BigR * wst * delta_phi
          J2_int = J2_int + eta_T_ohm * (ZJ0/BigR)**2.d0 * xjac * BigR * wst * delta_phi
        else
#ifdef WITH_Impurities
          D_ext = D_ext + (r0-rn0) * xjac * BigR * wst * delta_phi
#ifdef WITH_TiTe
          P_e_ext = P_e_ext + (r0+alpha_e*rn0) * T0e * xjac * BigR * wst * delta_phi
          P_i_ext = P_i_ext + (r0+alpha_i*rn0) * T0i * xjac * BigR * wst * delta_phi
          P_ext   = P_e_ext + P_i_ext
#else /* WITH_TiTe */
          P_ext = P_ext + (r0+alpha_imp*rn0) * T0   * xjac * BigR * wst * delta_phi
          P_e_ext = P_ext / 2.
          P_i_ext = P_e_ext
#endif /* WITH_TiTe */
#else /* WITH_Impurities */
          D_ext = D_ext + r0 * xjac * BigR * wst * delta_phi
#ifdef WITH_TiTe
          P_e_ext = P_e_ext + r0 * T0e * xjac * BigR * wst * delta_phi
          P_i_ext = P_i_ext + r0 * T0i * xjac * BigR * wst * delta_phi
          P_ext   = P_e_ext + P_i_ext
#else /* WITH_TiTe */
          P_ext = P_ext + r0 * T0   * xjac * BigR * wst * delta_phi
          P_e_ext = P_ext / 2.
          P_i_ext = P_e_ext
#endif /* WITH_TiTe */
#endif /* WITH_Impurities */
          C_ext = C_ext - zj0 / BigR * xjac *        wst * delta_phi  ! 2D integral
          H_ext = H_ext + heat_source     * xjac * BigR * wst * delta_phi
          S_ext = S_ext + particle_source * xjac * BigR * wst * delta_phi
          VP_ext = VP_ext + r0 * vpar0**2 * BB2 * xjac * BigR * wst * delta_phi
          VK_ext = VK_ext + r0 * (dudx**2 + dudy**2) * BigR**2 * xjac * BigR * wst * delta_phi
          VM_ext = VM_ext + (dpsidx**2+dpsidy**2)/BigR**2 * xjac * BigR * wst * delta_phi
          J2_ext = J2_ext + eta_T_ohm * (ZJ0/BigR)**2.d0 * xjac * BigR * wst * delta_phi
        endif

      enddo
    enddo
  enddo

enddo
!$omp end do
!$omp end parallel

!------ Calculate boundary fluxes --------------------------------------------------------
!--- go through the boundary elements
do m_bndelem = 1, bnd_elm_list%n_bnd_elements

  bndelem = bnd_elm_list%bnd_element(m_bndelem)
  elm_k   = element_list%element(bndelem%element)
  mv1     = bnd_elm_list%bnd_element(m_bndelem)%side
  m_elm   = bnd_elm_list%bnd_element(m_bndelem)%element

  !--- calculate values at gaussian points on the element
  x_g_1D(:)  = 0.d0; x_s_1D(:)  = 0.d0;  x_t_1D(:)    = 0.d0;
  y_g_1D(:)  = 0.d0; y_s_1D(:)  = 0.d0;  y_t_1D(:)    = 0.d0;

  eq_g_1D(:,:,:) = 0.d0; eq_s_1D(:,:,:) = 0.d0; delta_g_1D(:,:,:) = 0.d0;

  do k_vertex = 1, 2
    do k_dof = 1, 2
      k_node      = bndelem%vertex(k_vertex)
      k_dir       = bndelem%direction(k_vertex,k_dof)
      k_size      = bndelem%size(k_vertex,k_dof)
      node_k      = node_list%node(k_node)
  
      x_g_1D(:)   = x_g_1D(:)  + node_k%x(1,k_dir,1) * k_size * H1  (k_vertex,k_dof,:)
      y_g_1D(:)   = y_g_1D(:)  + node_k%x(1,k_dir,2) * k_size * H1  (k_vertex,k_dof,:)
      x_s_1D(:)   = x_s_1D(:)  + node_k%x(1,k_dir,1) * k_size * H1_s(k_vertex,k_dof,:)
      y_s_1D(:)   = y_s_1D(:)  + node_k%x(1,k_dir,2) * k_size * H1_s(k_vertex,k_dof,:)

      do k=1,n_var
        do mp=1, n_plane 
          do in=1,n_tor
            eq_g_1D(mp,k,:) = eq_g_1D(mp,k,:) + node_k%values(in,k_dir,k) * k_size * H1(k_vertex,k_dof,:)   * HZ(in,mp)
            eq_s_1D(mp,k,:) = eq_s_1D(mp,k,:) + node_k%values(in,k_dir,k) * k_size * H1_s(k_vertex,k_dof,:) * HZ(in,mp)

            delta_g_1D(mp,k,:) = delta_g_1D(mp,k,:) + node_k%deltas(in,k_dir,k) * k_size * H1(k_vertex,k_dof,:)   * HZ(in,mp)
          enddo
        enddo
      enddo

    end do
  end do


  !--- Find out correct sign of the normal (it has to point outwards the domain)
  !---------------------------------------------------------------------------------- 
  ! --- Calculate an inside point on the element to calculate the
  ! direction of bnd normals
  call basisfunctions(xgauss(2),xgauss(2), G)  
  R_c = 0.d0 ;  Z_c = 0.d0 
  do i = 1, n_vertex_max
    do j = 1, n_order+1
      node_k = node_list%node(elm_k%vertex(i)) 
      R_c    = R_c + node_k%x(1,j,1) * elm_k%size(i,j) * G(i,j)
      Z_c    = Z_c + node_k%x(1,j,2) * elm_k%size(i,j) * G(i,j)
    enddo
  enddo  
  vec_inside = (/ R_c - x_g_1D(2), Z_c - y_g_1D(2) /)       ! vector pointing towards the domain
  grad_t     = (/ -y_s_1D(2) , x_s_1D(2) /)     ! gradient of the coordinate t (normal to the boundary here)
  sign_out   = -1.d0 * sign( 1.d0, ( vec_inside(1)*grad_t(1) + vec_inside(2)*grad_t(2) ) )  
  !--------------------------------------------------------------------------------

  !--- Integrate quantity on the element surface
  do ms=1, n_gauss

    s_or_t = xgauss(ms) 

    ! --- Which s and t values correspond to the current point and is the
    !     boundary element an s=const or t=const side of the 2D element?
    select case (mv1)
    case (1)
      sg = s_or_t;  tg = 0.d0;   
    case (2)
      sg = 1.d0;    tg = s_or_t; 
    case (3)
      sg = s_or_t;  tg = 1.d0;  
    case (4)
      sg = 0.d0;    tg = s_or_t; 
    end select

    call interp_RZ(node_list,element_list,m_elm,sg,tg,R,R_s,R_t,R_st,R_ss,R_tt,Z,Z_s,Z_t,Z_st,Z_ss,Z_tt)

    BigR   = R
    xjac   = R_s * Z_t - R_t * Z_s
    grad_t = (/ -y_s_1D(ms) , x_s_1D(ms) /)   ! --- normal vector to the boundary 

    do mp=1, n_plane

      ps0      = eq_g_1D(mp,var_psi ,ms)  !--- here sbnd is the direction along the boundary!!
      ps0_sbnd = eq_s_1D(mp,var_psi ,ms)
      u0_sbnd  = eq_s_1D(mp,var_u   ,ms)
      zj0      = eq_g_1D(mp,var_Zj  ,ms) 
      r0       = eq_g_1D(mp,var_rho ,ms) 
      r0_corr  = corr_neg_dens1(r0)
      T0       = eq_g_1D(mp,var_T   ,ms) 
#ifdef WITH_TiTe
      T0i    = eq_g_1D(mp,var_Ti,ms)
      T0e    = eq_g_1D(mp,var_Te,ms)
#else
      T0i    = eq_g_1D(mp,var_T,ms) /2.d0
      T0e    = eq_g_1D(mp,var_T,ms) /2.d0
#endif

#ifdef WITH_Vpar
      vpar0    = eq_g_1D(mp,var_vpar,ms)
#else
      vpar0    = 0.d0
#endif

#if (defined WITH_Neutrals) || (defined WITH_Impurities)
      rn0      = eq_g_1D(mp,var_rhon,ms)
      rn0_corr = corr_neg_dens1(rn0) ! Correction for negative rn0
#else
      rn0      = 0.d0
      rn0_corr = 0.d0 
#endif

      if (keep_current_prof) then
        call current(xpoint, xcase, R, Z, Z_xpoint, ps0, psi_axis, psi_bnd, current_source)
      else
        current_source = 0.d0
      endif
 
      !--- calculate derivates in real s, t (s_1D is a coordinate that can be s or t)
      psi_s  = 0.d0; psi_t  = 0.d0;
      u_s    = 0.d0; u_t    = 0.d0; u_p = 0.d0;
      rho_s  = 0.d0; rho_t  = 0.d0;
      rhon_s = 0.d0; rhon_t = 0.d0;
      T_s    = 0.d0; T_t    = 0.d0;
      Ti_s   = 0.d0; Ti_t   = 0.d0;
      Te_s   = 0.d0; Te_t   = 0.d0;
      vpar_s = 0.d0; vpar_t = 0.d0; 

      do in = 1,n_tor
        call interp(node_list,element_list,m_elm,var_psi,in,sg,tg,PS,PS_s,PS_t,PS_st,PS_ss,PS_tt)
        psi_s = psi_s + PS_s * HZ(in,mp)
        psi_t = psi_t + PS_t * HZ(in,mp)
        
        if ( var_u /= 0 ) then
          call interp(node_list,element_list,m_elm,var_u,in,sg,tg,UU,UU_s,UU_t,UU_st,UU_ss,UU_tt)
          u_s   = u_s   + UU_s * HZ(in,mp)
          u_t   = u_t   + UU_t * HZ(in,mp)
          u_p   = u_p   + UU   * HZ_p(in,mp)
        else
          u_s = 0.d0
          u_t = 0.d0
          u_p = 0.d0
        end if

        call interp(node_list,element_list,m_elm,var_rho,in,sg,tg,RH,RH_s,RH_t,RH_st,RH_ss,RH_tt)
        rho_s = rho_s + RH_s * HZ(in,mp)
        rho_t = rho_t + RH_t * HZ(in,mp)

#ifdef WITH_TiTe
        call interp(node_list,element_list,m_elm,var_Ti,in,sg,tg,TT,TT_s,TT_t,TT_st,TT_ss,TT_tt)
        Ti_s = Ti_s + TT_s * HZ(in,mp)
        Ti_t = Ti_t + TT_t * HZ(in,mp)
        call interp(node_list,element_list,m_elm,var_Te,in,sg,tg,TT,TT_s,TT_t,TT_st,TT_ss,TT_tt)
        Te_s = Te_s + TT_s * HZ(in,mp)
        Te_t = Te_t + TT_t * HZ(in,mp)
        T_s = Te_s + Ti_s
        T_t = Te_t + Ti_t
#else
        call interp(node_list,element_list,m_elm,var_T,in,sg,tg,TT,TT_s,TT_t,TT_st,TT_ss,TT_tt)
        T_s = T_s + TT_s * HZ(in,mp)
        T_t = T_t + TT_t * HZ(in,mp)
        Te_s = T_s / 2.0
        Te_t = T_t / 2.0
        Ti_s = T_s / 2.0
        Ti_t = T_t / 2.0
#endif

#ifdef WITH_Vpar
        call interp(node_list,element_list,m_elm,var_Vpar,in,sg,tg,vp,vp_s,vp_t,vp_st,vp_ss,vp_tt)
        vpar_s = vpar_s + vp_s * HZ(in,mp)
        vpar_t = vpar_t + vp_t * HZ(in,mp)
#else
        vpar_s = 0.d0
        vpar_t = 0.d0
#endif

#if (defined WITH_Neutrals) || (defined WITH_Impurities)
        call interp(node_list,element_list,m_elm,var_rhon,in,sg,tg,rn,rn_s,rn_t,rn_st,rn_ss,rn_tt)
        rhon_s = rhon_s + rn_s * HZ(in,mp)
        rhon_t = rhon_t + rn_t * HZ(in,mp)
#else
        rhon_s = 0.d0
        rhon_t = 0.d0
#endif

      enddo

      dTedx  = (   Z_t * Te_s  - Z_s * Te_t  ) / xjac
      dTedy  = ( - R_t * Te_s  + R_s * Te_t  ) / xjac
      dTidx  = (   Z_t * Ti_s  - Z_s * Ti_t  ) / xjac
      dTidy  = ( - R_t * Ti_s  + R_s * Ti_t  ) / xjac
      dTdx   = (   Z_t * T_s   - Z_s * T_t   ) / xjac
      dTdy   = ( - R_t * T_s   + R_s * T_t   ) / xjac
      drhodx = (   Z_t * rho_s - Z_s * rho_t ) / xjac
      drhody = ( - R_t * rho_s + R_s * rho_t ) / xjac

      dpsidx = (   Z_t * psi_s - Z_s * psi_t ) / xjac
      dpsidy = ( - R_t * psi_s + R_s * psi_t ) / xjac

      vpar_x = (   Z_t * vpar_s - Z_s * vpar_t ) / xjac
      vpar_y = ( - R_t * vpar_s + R_s * vpar_t ) / xjac

      drhondx = (   Z_t * rhon_s - Z_s * rhon_t ) / xjac
      drhondy = ( - R_t * rhon_s + R_s * rhon_t ) / xjac

      BB2    = (F0*F0 + dpsidx*dpsidx + dpsidy*dpsidy) / BigR**2

      ! --- get normalized flux 
      psi_n = get_psi_n(ps0,Z)
 
#ifdef WITH_TiTe
      T0e_corr      = corr_neg_temp1(T0e)
      eta_T         = resistivity(eta, T0e_corr, T_max_eta, Te_0)  
      eta_T_ohm     = resistivity(eta_ohmic, T0e_corr, T_max_eta_ohm, Te_0)

      ZK_e_prof     = get_zk_eperp(psi_n)
      ZK_i_prof     = get_zk_iperp(psi_n)
 
      if (ZKpar_T_dependent) then
        ZK_e_par_T = ZK_e_par * (max(T0e,T_min)/Te_0)**( 2.5d0)
        ZK_i_par_T = ZK_i_par * (max(T0i,T_min)/Ti_0)**( 2.5d0)
      else
        ZK_e_par_T = ZK_e_par 
        ZK_i_par_T = ZK_i_par 
      endif
#else
      T0_corr       = corr_neg_temp1(T0)
      eta_T         = resistivity(eta, T0_corr, T_max_eta, T_0)  
      eta_T_ohm     = resistivity(eta_ohmic, T0_corr, T_max_eta_ohm, T_0)

      ZK_prof = get_zkperp(psi_n)
 
      if (ZKpar_T_dependent) then
        ZKpar_T = ZK_par * (max(T0,T_min)/T_0)**( 2.5d0)
      else
        ZKpar_T = ZK_par
      endif
#endif

#ifdef WITH_Impurities
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

      T0e_corr = corr_neg_temp1(T0e)
      dT0e_corr_dT = dcorr_neg_temp_dT(T0e,(/5.d-1,5.d-1/),T_min)
      Te_corr_eV = T0e_corr/(EL_CHG*MU_ZERO*central_density*1.d20)
      Te_eV = T0e/(EL_CHG*MU_ZERO*central_density*1.d20)
   
      if (allocated(P_imp)) deallocate(P_imp)
      allocate(P_imp(0:imp_adas(index_main_imp)%n_Z))
      call imp_cor(index_main_imp)%interp_linear(density=20.,temperature=log10(Te_corr_eV*EL_CHG/K_BOLTZ),&
                                    p_out=P_imp,z_avg=Z_imp,z_avg_Te=dZ_imp_dT)
      ! Convert gradient in T(K) in to gradient in T (eV)
      dZ_imp_dT = dZ_imp_dT *EL_CHG / K_BOLTZ
      ! Derivative wrt to T, with T in JOREK units
      dZ_imp_dT = dZ_imp_dT / (EL_CHG*MU_ZERO*central_density*1.d20)
      dZ_imp_dT = dZ_imp_dT * dT0e_corr_dT
#ifdef WITH_TiTe
      alpha_i       = m_i_over_m_imp - 1.
      alpha_e       = m_i_over_m_imp*Z_imp - 1.
      dalpha_e_dT   = m_i_over_m_imp*dZ_imp_dT

      ne_SI        = (r0_corr + alpha_e * rn0_corr) * 1.d20 * central_density ! electron density (SI)
      ne_JOREK     = r0_corr + alpha_e * rn0_corr ! Electron density in JOREK unit
      ne_JOREK     = corr_neg_dens(ne_JOREK,(/1.d-1,1.d-1/),1.d-3) ! Correction for negative electron density
                                                          ! Too small rho_1 will cause a problem
      dne_JOREK_dx = drhodx + alpha_e * drhondx + rn0_corr * dalpha_e_dT * dTedx
      dne_JOREK_dy = drhody + alpha_e * drhondy + rn0_corr * dalpha_e_dT * dTedy
#else /* WITH_TiTe */
      alpha_imp    = 0.5*m_i_over_m_imp*(Z_imp+1.) - 1.
      beta_imp     = m_i_over_m_imp*Z_imp - 1.
      dbeta_imp_dT = m_i_over_m_imp*dZ_imp_dT
      ne_SI        = (r0_corr + beta_imp * rn0_corr) * 1.d20 * central_density !electron density (SI)
      ne_JOREK     = r0_corr + beta_imp * rn0_corr ! Electron density in JOREK unit
      ne_JOREK     = corr_neg_dens(ne_JOREK,(/1.d-1,1.d-1/),1.d-3) ! Correction for negative electron density
                                                             ! Too small rho_1 will cause a problem
      dne_JOREK_dx = drhodx + alpha_imp * drhondx + rn0_corr * dbeta_imp_dT * dTedx
      dne_JOREK_dy = drhody + alpha_imp * drhondy + rn0_corr * dbeta_imp_dT * dTedy
#endif /* WITH_TiTe */

      ! Calculate the effective charge of all species
      Z_eff        = 0.

      ! First get the value of Z_eff
      Z_eff        = r0_corr - rn0_corr
      do ion_i=1, imp_adas(index_main_imp)%n_Z
        Z_eff      = Z_eff + m_i_over_m_imp * rn0_corr * P_imp(ion_i) * real(ion_i,8)**2
      end do
      Z_eff        = Z_eff / ne_JOREK

      if (Z_eff < 1) Z_eff = 1.

      ! This is to represent the dependence on Z_eff in resistivity
      if ( eta_T_dependent ) then
        eta_coef     = Z_eff*(1.+1.198*Z_eff+0.222*Z_eff**2)/(1.+2.966*Z_eff+0.753*Z_eff**2)
        eta_coef     = eta_coef / ((1.+1.198+0.222)/(1.+2.966+0.753))
        eta_T        = eta_T * eta_coef
        eta_T_ohm    = eta_T_ohm * eta_coef
      endif

      dPedx  = ne_JOREK * dTedx + T0e * dne_JOREK_dx
      dPedy  = ne_JOREK * dTedy + T0e * dne_JOREK_dy
      dPidx  = (r0 + rn0) * dTidx + T0i * (drhodx + drhondx)
      dPidy  = (r0 + rn0) * dTidy + T0i * (drhody + drhondy)
      dPdx   = dPedx + dPidx
      dPdy   = dPedy + dPidy
#else /* WITH_Impurities */
      dPdx   = r0 * dTdx + T0 * drhodx
      dPdy   = r0 * dTdy + T0 * drhody
#endif /* WITH_Impurities */

      D_prof  = get_dperp (psi_n)

      pflow       = - gamma/(gamma-1.d0) * r0 * T0 * vpar0 * ps0_sbnd * sign_out 
      kinflow     = - 0.5d0*r0*vpar0**3.d0*BB2* ps0_sbnd * sign_out 

#ifdef WITH_TiTe
      cond_par    =   ZK_e_par_T *(dTedx * dpsidy - dTedy * dpsidx) /BigR/BB2 * ps0_sbnd * sign_out / (gamma-1.d0)&
                    + ZK_i_par_T *(dTidx * dpsidy - dTidy * dpsidx) /BigR/BB2 * ps0_sbnd * sign_out / (gamma-1.d0) 
      cond_perp   = - ZK_e_prof *(dTedx*grad_t(1) + dTedy*grad_t(2)) * BigR              * sign_out / (gamma-1.d0)&
                    - ZK_e_prof *(dTedx * dpsidy  - dTedy * dpsidx) /BigR/BB2 * ps0_sbnd * sign_out / (gamma-1.d0)&
                    - ZK_i_prof *(dTidx*grad_t(1) + dTidy*grad_t(2)) * BigR              * sign_out / (gamma-1.d0)&
                    - ZK_i_prof *(dTidx * dpsidy  - dTidy * dpsidx) /BigR/BB2 * ps0_sbnd * sign_out / (gamma-1.d0)
#else
      cond_par    =   ZKpar_T *( dTdx * dpsidy  - dTdy * dpsidx ) /BigR/BB2 * ps0_sbnd * sign_out / (gamma-1.d0) 
      cond_perp   = - ZK_prof *( dTdx*grad_t(1) + dTdy*grad_t(2)) * BigR               * sign_out / (gamma-1.d0) &
                    - ZK_prof *( dTdx * dpsidy  - dTdy * dpsidx ) /BigR/BB2 * ps0_sbnd * sign_out / (gamma-1.d0) 
#endif

      Dpar_part_flow    =   D_par  * (drhodx*dpsidy    - drhody*dpsidx    )/BigR/BB2 * ps0_sbnd * sign_out 
      Dperp_part_flow   = - D_prof * (drhodx*grad_t(1) + drhody*grad_t(2) ) * BigR              * sign_out &                               
                          - D_prof * (drhodx*dpsidy    - drhody*dpsidx    )/BigR/BB2 * ps0_sbnd * sign_out

      vpar_part_flow    =  - r0 *  vpar0     * ps0_sbnd * sign_out 
      vperp_part_flow   =    r0 * BigR**2.d0 * u0_sbnd  * sign_out 

      neut_part_flow    = - (D_neutral_x*drhondx*grad_t(1) + D_neutral_y*drhondy*grad_t(2)) * BigR  * sign_out

      viscopar_f  = visco_par * (F0/BigR)**2.d0 *  (vpar_x*grad_t(1) + vpar_y *grad_t(2) ) &
                  * sign_out  * BigR * vpar0

!      dpsi_dt     = BigR*(psi_s*u_t - psi_t*u_s)/xjac + eta_T*(zj0-current_source) - F0*u_p 
      dpsi_dt      = delta_g_1D(mp, var_psi, ms) / tstep  
      poynting_tmp = dpsi_dt * (dpsidx*grad_t(1) + dpsidy*grad_t(2)) * sign_out / BigR 


      vn_p0         = vn_p0          +   pflow      * wgauss(ms) * delta_phi 
      kinpar_flux   = kinpar_flux    + kinflow      * wgauss(ms) * delta_phi 
      qn_par        = qn_par         + cond_par     * wgauss(ms) * delta_phi 
      qn_perp       = qn_perp        + cond_perp    * wgauss(ms) * delta_phi 

      Dpar_part_flux   =  Dpar_part_flux     + Dpar_part_flow    * wgauss(ms) * delta_phi 
      Dperp_part_flux  =  Dperp_part_flux    + Dperp_part_flow   * wgauss(ms) * delta_phi 
      vpar_part_flux   =  vpar_part_flux     + vpar_part_flow    * wgauss(ms) * delta_phi 
      vperp_part_flux  =  vperp_part_flux    + vperp_part_flow   * wgauss(ms) * delta_phi  
      neut_part_flux   =  neut_part_flux     + neut_part_flow    * wgauss(ms) * delta_phi 

      viscopar_flux = viscopar_flux  +  viscopar_f  * wgauss(ms) * delta_phi
      poynting_flux = poynting_flux  + poynting_tmp * wgauss(ms) * delta_phi

    enddo
  enddo

enddo !--- bnd elements, end of calculation of boundary fluxes

! --- gather contribution from all MPI processes
#ifndef NOMPIVERSION
call MPI_AllReduce(D_int,density_in,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(D_ext,density_out,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(P_int,pressure_in,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(P_ext,pressure_out,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(P_e_int,pressure_e_in,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(P_e_ext,pressure_e_out,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(P_i_int,pressure_i_in,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(P_i_ext,pressure_i_out,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(C_intern,current_in,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(C_ext,current_out,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(R2curr_tmp,      R2curr,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(Zcurr_tmp , Z_curr_cent,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(Vol,Volume,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(area1,area,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(D_tot,density_tot,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(P_tot,pressure,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(P_e_tot,pressure_e,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(P_i_tot,pressure_i,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(H_ext,heating_out,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(H_int,heating_in,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(S_ext,source_out,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(S_int,source_in,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(VP_int,kin_par_in,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(VP_ext,kin_par_out,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(VP_tot,kin_par_tot,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(VK_int,kin_perp_in,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(VK_ext,kin_perp_out,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(VK_tot,kin_perp_tot,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(VM_int,mag_in,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(VM_ext,mag_out,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(VM_tot,mag_tot,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(J2_int,ohm_in,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(J2_ext,ohm_out,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(J2_tot,ohm_tot,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(heli_tot, helicity_tot,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(thm_wk_tot, thermal_work_tot,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(mag_wk_tot, mag_work_tot,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(vpar_disp_tot, viscopar_dissip_tot,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(fric_disp_tot, friction_dissip_tot,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(mag_src_tot, mag_source_tot,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(varmin,V_min,n_var,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(varmax,V_max,n_var,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
#if (defined WITH_Neutrals) || (defined WITH_Impurities)
call MPI_AllReduce(local_radiation, total_radiation,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(local_E_ion, total_E_ion,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(local_P_ei, total_P_ei,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(local_P_ion, total_P_ion,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(local_radiation_phi, total_radiation_phi,n_plane,&
                   MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
#endif /* (defined WITH_Neutrals) || (defined WITH_Impurities) */
#else /* NOMPIVERSION */
density_in           = D_int
density_out          = D_ext
pressure_in          = P_int
pressure_out         = P_ext
pressure_e_in        = P_e_int
pressure_e_out       = P_e_ext
pressure_i_in        = P_i_int
pressure_i_out       = P_i_ext
current_in           = C_intern
current_out          = C_ext
R2curr               = R2curr_tmp
Z_curr_cent          = Zcurr_tmp
Volume               = Vol
area                 = area1
density_tot          = D_tot
pressure             = P_tot
pressure_e           = P_e_tot
pressure_i           = P_i_tot
heating_out          = H_ext
heating_in           = H_int
source_out           = S_ext
source_in            = S_int
kin_par_in           = VP_int
kin_par_out          = VP_ext
kin_par_tot          = VP_tot
kin_perp_in          = VK_int
kin_perp_out         = VK_ext
kin_perp_tot         = VK_tot
mag_in               = VM_int
mag_out              = VM_ext
mag_tot              = VM_tot
ohm_in               = J2_int
ohm_out              = J2_ext
ohm_tot              = J2_tot
helicity_tot         = heli_tot
thermal_work_tot     = thm_wk_tot
mag_work_tot         = mag_wk_tot
viscopar_dissip_tot  = vpar_disp_tot
friction_dissip_tot  = fric_disp_tot
mag_source_tot       = mag_src_tot
V_min                = varmin
V_max                = varmax

#if (defined WITH_Neutrals) || (defined WITH_Impurities)
total_radiation      = local_radiation
total_E_ion          = local_E_ion
total_P_ei           = local_P_ei
total_P_ion          = local_P_ion
total_radiation_phi  = local_radiation_phi
#endif /* (defined WITH_Neutrals) || (defined WITH_Impurities) */
#endif /* NOMPIVERSION */

if (use_pellet) then
#ifndef NOMPIVERSION
  call MPI_AllReduce(local_pellet_particles,total_pellet_particles,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_AllReduce(local_plasma_particles,total_plasma_particles,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_AllReduce(local_pellet_volume,total_pellet_volume,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
#else /* NOMPIVERSION */
  total_pellet_particles = local_pellet_particles
  total_plasma_particles = local_plasma_particles
  total_pellet_volume    = local_pellet_volume
#endif /* NOMPIVERSION */
endif

#if (defined WITH_Neutrals) || (defined WITH_Impurities)
if (using_spi) then   
   do spi_i=1, n_spi_tot
#ifndef NOMPIVERSION
      call MPI_AllReduce(local_source_volume(spi_i),pellets(spi_i)%spi_vol,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_AllReduce(local_source_volume_drift(spi_i),pellets(spi_i)%spi_vol_drift,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
#else /* NOMPIVERSION */
      pellets(spi_i)%spi_vol = local_source_volume(spi_i)
      pellets(spi_i)%spi_vol_drift = local_source_volume_drift(spi_i)
#endif /* NOMPIVERSION */
   end do
   deallocate(local_source_volume)
   deallocate(local_source_volume_drift)
end if
if (allocated(local_source_volume)) deallocate(local_source_volume) !In case of dummy array
if (allocated(local_source_volume_drift)) deallocate(local_source_volume_drift)
#endif

#if (defined WITH_Neutrals) || (defined WITH_Impurities)
#ifndef NOMPIVERSION
  call MPI_AllReduce(local_n_particles_inj, total_n_particles_inj,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_AllReduce(local_n_particles, total_n_particles,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
#else /* NOMPIVERSION */
  total_n_particles_inj = local_n_particles_inj
  total_n_particles     = local_n_particles
#endif /* NOMPIVERSION */
  neut_particles_tot    = total_n_particles / (central_density * 1.d20)
#else
  neut_particles_tot = 0.d0
#endif

! --- Normalization factors
rho_norm = central_density*1.d20 * central_mass * MASS_PROTON 
t_norm   = sqrt(MU_zero*rho_norm)

if (units == SI_UNITS) then
  fact_mu0  = 1.d0/mu_zero
  fact_flux = 1.d0/(mu_zero*t_norm)  
  t_norm2   = t_norm
  fact_part = central_density * 1.d20
else
  fact_mu0  = 1.d0
  fact_flux = 1.d0
  t_norm2   = 1.d0
  fact_part = 1.d0
endif

! --- Volume integrals
current_in           = n_period * current_in  * fact_mu0  / (2.d0 * PI)
current_out          = n_period * current_out * fact_mu0  / (2.d0 * PI)
R2curr               = n_period * R2curr      * fact_mu0  / (2.d0 * PI)
Z_curr_cent          = n_period * Z_curr_cent * fact_mu0  / (2.d0 * PI)
pressure             = n_period * pressure    * fact_mu0  / (GAMMA-1.d0)
pressure_in          = n_period * pressure_in * fact_mu0  / (GAMMA-1.d0)
pressure_out         = n_period * pressure_out* fact_mu0  / (GAMMA-1.d0)
pressure_e           = n_period * pressure_e  * fact_mu0  / (GAMMA-1.d0)
pressure_e_in        = n_period * pressure_e_in * fact_mu0  / (GAMMA-1.d0)
pressure_e_out       = n_period * pressure_e_out* fact_mu0  / (GAMMA-1.d0)
pressure_i           = n_period * pressure_i    * fact_mu0  / (GAMMA-1.d0)
pressure_i_in        = n_period * pressure_i_in * fact_mu0  / (GAMMA-1.d0)
pressure_i_out       = n_period * pressure_i_out* fact_mu0  / (GAMMA-1.d0)
kin_par_tot          = n_period * kin_par_tot * fact_mu0  * 0.5d0
kin_par_in           = n_period * kin_par_in  * fact_mu0  * 0.5d0
kin_par_out          = n_period * kin_par_out * fact_mu0  * 0.5d0
kin_perp_tot         = n_period * kin_perp_tot* fact_mu0  * 0.5d0
kin_perp_in          = n_period * kin_perp_in * fact_mu0  * 0.5d0
kin_perp_out         = n_period * kin_perp_out* fact_mu0  * 0.5d0
mag_tot              = n_period * mag_tot     * fact_mu0  * 0.5d0
mag_in               = n_period * mag_in      * fact_mu0  * 0.5d0
mag_out              = n_period * mag_out     * fact_mu0  * 0.5d0
ohm_tot              = n_period * ohm_tot     * fact_flux
ohm_in               = n_period * ohm_in      * fact_flux
ohm_out              = n_period * ohm_out     * fact_flux
heating_out          = n_period * heating_out * fact_flux / (GAMMA-1.d0)
heating_in           = n_period * heating_in  * fact_flux / (GAMMA-1.d0)
source_out           = n_period * source_out  * fact_part / t_norm2
source_in            = n_period * source_in   * fact_part / t_norm2
density_tot          = n_period * density_tot * fact_part 
density_in           = n_period * density_in  * fact_part
density_out          = n_period * density_out * fact_part
neut_particles_tot   = n_period * neut_particles_tot * fact_part
helicity_tot         = n_period * helicity_tot
thermal_work_tot     = n_period * thermal_work_tot    * fact_flux 
mag_work_tot         = n_period * mag_work_tot        * fact_flux
viscopar_dissip_tot  = n_period * viscopar_dissip_tot * fact_flux
friction_dissip_tot  = n_period * friction_dissip_tot * fact_flux
mag_source_tot       = n_period * mag_source_tot      * fact_flux
volume               = n_period * volume
area                 = n_period * area / (2.d0 * PI)

#if (defined WITH_Neutrals) || (defined WITH_Impurities)
total_radiation     = n_period * total_radiation
total_radiation_phi = n_period * total_radiation_phi
total_E_ion         = n_period * total_E_ion
total_P_ei          = n_period * total_P_ei
total_P_ion         = n_period * total_P_ion
#endif

! --- Boundary integrals
vn_p0                =  n_period * vn_p0          * fact_flux 
qn_par               =  n_period * qn_par         * fact_flux  
qn_perp              =  n_period * qn_perp        * fact_flux 
kinpar_flux          =  n_period * kinpar_flux    * fact_flux  
viscopar_flux        =  n_period * viscopar_flux  * fact_flux
Dpar_part_flux       =  n_period * Dpar_part_flux * fact_part / t_norm2  
Dperp_part_flux      =  n_period * Dperp_part_flux* fact_part / t_norm2
vpar_part_flux       =  n_period * vpar_part_flux * fact_part / t_norm2
vperp_part_flux      =  n_period * vperp_part_flux* fact_part / t_norm2
neut_part_flux       =  n_period * neut_part_flux * fact_part / t_norm2
poynting_flux        =  n_period * poynting_flux  * fact_flux

! --- Derived quantities
E_tot        = mag_tot + pressure     + kin_par_tot + kin_perp_tot 
E_in         = mag_in  + pressure_in  + kin_par_in  + kin_perp_in 
E_out        = mag_out + pressure_out + kin_par_out + kin_perp_out 
current_tot  = current_in + current_out
heating_tot  = heating_in + heating_out
source_tot   = source_in  + source_out
Bgeo         = F0 / R_geo
current_MA   = current_in * 1.d-6 * (1.d0/fact_mu0) * (1/mu_zero)
beta_p       = 4.d0 * pressure_in/(R_geo * current_in**2 )     * (GAMMA-1)*fact_mu0
beta_t       = 2.d0 * pressure_in / volume / Bgeo**2           * (GAMMA-1)/fact_mu0
beta_n       = 100.d0 * beta_t * Bgeo/current_MA * sqrt(area/PI)
li3          = 2.d0 * mag_in /0.5  /( current_in**2 * R_geo ) * fact_mu0
li3_tot      = 2.d0 * mag_tot/0.5  /(current_tot**2 * R_geo ) * fact_mu0
sheath_heatflux =  gamma_stangeby * (gamma-1)/(2.d0*gamma) * vn_p0 ! the factor comes to obtain n T_e v from vn_p0
R_curr_cent  = sqrt(R2curr / current_tot) 
Z_curr_cent  = Z_curr_cent / current_tot 


! --- Externally calculated quantities
! --- Halo currents
call integrated_normal_bnd_curr(node_list, bnd_node_list, bnd_elm_list, I_halo, TPF, .false.)

! --- Safety factor at important locations
surface_list%n_psi = 4 
allocate( surface_list%psi_values(surface_list%n_psi) )
allocate( qval(surface_list%n_psi), radav(surface_list%n_psi) )
surface_list%psi_values(1) = ES%psi_axis + (ES%psi_bnd - ES%psi_axis) * 0.01d0
surface_list%psi_values(2) = ES%psi_axis + (ES%psi_bnd - ES%psi_axis) * 0.02d0
surface_list%psi_values(3) = ES%psi_axis + (ES%psi_bnd - ES%psi_axis) * 0.95d0 
surface_list%psi_values(4) = ES%psi_axis + (ES%psi_bnd - ES%psi_axis) * 0.99d0

call find_flux_surfaces(my_id,xpoint, xcase, node_list, element_list, surface_list)
call determine_q_profile(node_list, element_list, surface_list, ES%psi_axis, ES%psi_xpoint,    &
     ES%Z_xpoint, qval, radav)

q02 = qval(2)
q95 = qval(3)
q99 = qval(4) 


if (my_id .eq. 0) then 

  ! --- Get current time
  if (index_now > 0) then
    xt = xtime(index_now)
  else
    xt = 0.d0
  endif

  ! --- Export or save quantities
  res(1)  =  xt*t_norm2
  loop_expr: do iexpr = 1, expr_list%n_expr
            
    select case ( trim(expr_list%expr(iexpr)%name) )

      case ( 'index_now' )
        res(iexpr+1) = real(index_now) 

      case ( 'psi_axis' )
        res(iexpr+1) = ES%psi_axis 

      case ( 'R_axis' )
        res(iexpr+1) = ES%R_axis 

      case ( 'Z_axis' )
        res(iexpr+1) = ES%Z_axis

      case ( 'R_curr_cent' )
        res(iexpr+1) = R_curr_cent 

      case ( 'Z_curr_cent' )
        res(iexpr+1) = Z_curr_cent

      case ( 'psi_bnd' )
        res(iexpr+1) = ES%psi_bnd 

      case ( 'R_bnd' )
        res(iexpr+1) = ES%R_bnd 

      case ( 'Z_bnd' )
        res(iexpr+1) = ES%Z_bnd

      case ( 'E_tot' )
        res(iexpr+1) = E_tot 

      case ( 'E_in' )
        res(iexpr+1) = E_in 

      case ( 'E_out' )
        res(iexpr+1) = E_out

      case ( 'Wmag_tot' )
        res(iexpr+1) = mag_tot 

      case ( 'Wmag_in' )
        res(iexpr+1) = mag_in 

      case ( 'Wmag_out' )
        res(iexpr+1) = mag_out 
     
      case ( 'Thermal_tot' )
        res(iexpr+1) = pressure 

      case ( 'Thermal_in' )
        res(iexpr+1) = pressure_in 

      case ( 'Thermal_out' )
        res(iexpr+1) = pressure_out 

      case ( 'Thermal_e_tot' )
        res(iexpr+1) = pressure_e

      case ( 'Thermal_e_in' )
        res(iexpr+1) = pressure_e_in 

      case ( 'Thermal_e_out' )
        res(iexpr+1) = pressure_e_out 

      case ( 'Thermal_i_tot' )
        res(iexpr+1) = pressure_i

      case ( 'Thermal_i_in' )
        res(iexpr+1) = pressure_i_in 

      case ( 'Thermal_i_out' )
        res(iexpr+1) = pressure_i_out 

      case ( 'Kin_par_tot' )
        res(iexpr+1) = kin_par_tot 

      case ( 'Kin_par_in' )
        res(iexpr+1) = kin_par_in 

      case ( 'Kin_par_out' )
        res(iexpr+1) = kin_par_out 

      case ( 'Kin_perp_tot' )
        res(iexpr+1) = kin_perp_tot 

      case ( 'Kin_perp_in' )
        res(iexpr+1) = kin_perp_in 

      case ( 'Kin_perp_out' )
        res(iexpr+1) = kin_perp_out 

      case ( 'Part_tot' ) 
        res(iexpr+1) = density_tot 

      case ( 'Part_in' ) 
        res(iexpr+1) = density_in 

      case ( 'Part_out' ) 
        res(iexpr+1) = density_out 

      case ( 'NPart_tot' ) 
        res(iexpr+1) = neut_particles_tot  

      case ( 'Helicity_tot' )
        res(iexpr+1) = helicity_tot 

      case ( 'Mag_work_tot' )
        res(iexpr+1) = mag_work_tot

      case ( 'Thm_work_tot' )
        res(iexpr+1) = thermal_work_tot

      case ( 'Part_src_tot' )
        res(iexpr+1) = source_tot

      case ( 'Part_src_in' )
        res(iexpr+1) = source_in

      case ( 'Part_src_out' )
        res(iexpr+1) = source_out

      case ( 'Heat_src_tot' )
        res(iexpr+1) = heating_tot

      case ( 'Heat_src_in' )
        res(iexpr+1) = heating_in

      case ( 'Heat_src_out' )
        res(iexpr+1) = heating_out

      case ( 'Viscpar_diss' )
        res(iexpr+1) = viscopar_dissip_tot

      case ( 'Friction_diss' )
        res(iexpr+1) = friction_dissip_tot

      case ( 'Wmag_src_tot' )
        res(iexpr+1) = mag_source_tot 

      case ( 'Ohmic_tot' )
        res(iexpr+1) = ohm_tot 

      case ( 'Ohmic_in' )
        res(iexpr+1) = ohm_in 

      case ( 'Ohmic_out' )
        res(iexpr+1) = ohm_out 

      case ( 'P_vn' )
        res(iexpr+1) = vn_p0 

      case ( 'qn_par' )
        res(iexpr+1) = qn_par 

      case ( 'qn_perp' )
        res(iexpr+1) = qn_perp

      case ( 'sheath_heat' )
        res(iexpr+1) = sheath_heatflux 

      case ( 'kinpar_flux' )
        res(iexpr+1) = kinpar_flux

      case ( 'Poynting_flx' )
        res(iexpr+1) = poynting_flux

      case ( 'vispar_flux' )
        res(iexpr+1) = viscopar_flux

      case ( 'Dpar_pt_flx' )
        res(iexpr+1) = Dpar_part_flux

      case ( 'Dperp_pt_flx' )
        res(iexpr+1) = Dperp_part_flux

      case ( 'vpar_pt_flx' )
        res(iexpr+1) = vpar_part_flux

      case ( 'vperp_pt_flx' )
        res(iexpr+1) = vperp_part_flux

      case ( 'neut_pt_flx' )
        res(iexpr+1) = neut_part_flux

      case ( 'Ip_tot' )
        res(iexpr+1) = current_tot 

      case ( 'Ip_in' )
        res(iexpr+1) = current_in 

      case ( 'Ip_out' )
        res(iexpr+1) = current_out 

      case ( 'li3' )
        res(iexpr+1) = li3

      case ( 'li3_tot' )
        res(iexpr+1) = li3_tot

      case ( 'beta_p' )
        res(iexpr+1) = beta_p 

      case ( 'beta_t' )
        res(iexpr+1) = beta_t 

      case ( 'beta_n' )
        res(iexpr+1) = beta_n 

      case ( 'area' )
        res(iexpr+1) = area 

      case ( 'volume' )
        res(iexpr+1) = volume

      case ( 'q02' )
        res(iexpr+1) = q02 

      case ( 'q95' )
        res(iexpr+1) = q95 

      case ( 'q99' )
        res(iexpr+1) = q99 

      case ( 'I_halo' )
        res(iexpr+1) = I_halo 

      case ( 'TPF_halo' )
        res(iexpr+1) = TPF 

      case ( 'LCFS_Rgeo' )
        res(iexpr+1) = ES%LCFS_Rgeo 

      case ( 'LCFS_a' )
        res(iexpr+1) = ES%LCFS_a

      case ( 'LCFS_epsilon' )
        res(iexpr+1) = ES%LCFS_epsilon

      case ( 'LCFS_kappa' )
        res(iexpr+1) = ES%LCFS_kappa

      case ( 'LCFS_deltaU' )
        res(iexpr+1) = ES%LCFS_deltaU 

      case ( 'LCFS_deltaL' )
        res(iexpr+1) = ES%LCFS_deltaL

    end select
            
  end do loop_expr

  ! ---- Print out some data 
  write(*,'(A,3e14.6,A)') ' Time : ',xt,xt*t_norm,t_norm, ' [s]'
  if (use_pellet) then 
    write(*,'(A,4e14.6)')   ' Integrals_3D, PELLET            : ',pellet_volume, total_pellet_volume, total_pellet_particles, total_plasma_particles
  endif
#if (defined WITH_Neutrals) || (defined WITH_Impurities)
  if (using_spi) then
    write(*,'(A)')   ' Integrals_3D, SPI               : '
    do i = 1, n_spi_tot
       if (pellets(i)%spi_radius > 0. .and. pellets(i)%spi_abl > 0.) then
          write(*,'(A,i14)')    "Pellet number                = ", i
          write(*,'(A,3f14.6)') "Pellet coordinates (R,Z,phi) = ", pellets(i)%spi_R, pellets(i)%spi_Z, pellets(i)%spi_phi
          write(*,'(A,3f14.6)') "Pellet velocity    (R,Z,phi) = ", pellets(i)%spi_Vel_R, pellets(i)%spi_Vel_Z, &
               pellets(i)%spi_Vel_RxZ
          write(*,'(A,3es14.6)')"Pellet ablation (radius,abl) = ", pellets(i)%spi_radius, pellets(i)%spi_abl
          write(*,'(A,f14.6)')  "Pellet species               = ", pellets(i)%spi_species

          ng_radius_tmp   = pellets(i)%spi_radius * ng_radius_ratio

          if (ng_radius_tmp < ng_radius_min) then
             ng_radius_tmp = ng_radius_min
          end if

          if (ns_deltaminrad .gt. 0.) then
             ! i.e., with poloidally elongated ablation cloud
             ! in this case the analytical formula below is approximate (usually it agrees with the numerical integral within a few percents)
             V_ns  = PI * pellets(i)%spi_R * ns_tor_norm * ng_radius_tmp * min(ns_deltaminrad,ng_radius_tmp)
             if (drift_distance /= 0.d0) then
               V_ns_drift  = PI * (pellets(i)%spi_R + drift_distance) * ns_tor_norm * ng_radius_tmp * min(ns_deltaminrad,ng_radius_tmp)
             end if
          else
             ! i.e., standard case with circular ablation cloud in the poloidal plane
             ! in this case the ablation source volume is given by the exact analytical formula as derived by E. Nardon
             V_ns  = PI * pellets(i)%spi_R * ns_tor_norm * ng_radius_tmp**2.d0
             if (drift_distance /= 0.d0) then
               V_ns_drift  = PI * (pellets(i)%spi_R + drift_distance) * ns_tor_norm * ng_radius_tmp**2.d0
             end if
          endif
          
          write(*,'(A,2es14.6,f14.6)') "Source vol (num,an,diff %)   = ", pellets(i)%spi_vol, V_ns, 1d2*(pellets(i)%spi_vol - V_ns)/V_ns
          if (abs((pellets(i)%spi_vol - V_ns)/V_ns) .gt. 0.1d0) write(*,*) "WARNING: Difference larger than 10% "

          if (drift_distance /= 0.d0) then 
            write(*,'(A,2es14.6,f14.6)') "Drifted source vol (num,an,diff %)   = ", pellets(i)%spi_vol_drift, V_ns_drift, 1d2*(pellets(i)%spi_vol_drift - V_ns_drift)/V_ns_drift
            if (abs((pellets(i)%spi_vol_drift - V_ns_drift)/V_ns_drift) .gt. 0.1d0) write(*,*) "WARNING: Difference larger than 10% "
          end if

          ! recommended ablation source radius in the poloidal direction from ng_radius / (R*ns_deltaphi) = B_pol/B_tor
          write(*,'(A,2f14.6)') "Source pol rad (actual,recom)= ", ng_radius_tmp, pellets(i)%spi_R * ns_deltaphi * pellets(i)%spi_grad_psi / abs(F0)
          if (drift_distance /= 0.d0) then 
            write(*,'(A,2f14.6)') "Drifted source pol rad (actual,recom)= ", ng_radius_tmp, (pellets(i)%spi_R + drift_distance) * ns_deltaphi * pellets(i)%spi_grad_psi_drift / abs(F0)
          end if
       end if
    end do
  endif
#endif
  write(*,'(A,2es14.6,A)') ' Volume                          : ',xt,volume,' [m^3]'
  write(*,'(A,4es14.6,A)') ' density  (total/in/out)         : ',xt,density_tot,  density_in,  density_out,'[ 10^20/m^3]'
  write(*,'(A,4es14.6,A)') ' pressure (total/in/out)         : ',xt,pressure/1.d6, pressure_in/1.d6, pressure_out/1.d6,' [MJ]'
  write(*,'(A,4es14.6,A)') ' kinetic parallel (total/in/out) : ',xt,kin_par_tot/1.d6, kin_par_in/1.d6, kin_par_out/1.d6,' [MJ]'
  write(*,'(A,4es14.6,A)') ' kinetic perp (total/in/out)     : ',xt,kin_perp_tot/1.d6, kin_perp_in/1.d6, kin_perp_out/1.d6,' [MJ]'
  write(*,'(A,4es14.6,A)') ' magnetic (total/in/out)         : ',xt,mag_tot/1.d6, mag_in/1.d6, mag_out/1.d6,' [MJ]'
  write(*,'(A,3es14.6,A)') ' current  (in/out)               : ',xt,current_in/1.d6, current_out/1.d6, ' [MA]'
  write(*,'(A,3es14.6,A)') ' heating  (in/out)               : ',xt,heating_in/1d6, heating_out/1.d6 ,' [MW]'
  write(*,'(A,3es14.6,A)') ' source   (in/out)               : ',xt,source_in, source_out,' [10^20/m^3/s]'
  write(*,'(A,4es14.6,A)') ' Ohmic    (in/out)               : ',xt,Ohm_tot/1.d6, Ohm_in/1.d6, Ohm_out/1.d6,' [MW]'

  write(*,'(A,2es14.6)')   ' li(3)                           : ',xt, li3 
  write(*,'(A,2es14.6)')   ' betap(1)                        : ',xt, beta_p

  write(*,'(A)')           ' sum ,time ,density_tot, pressure, Wkin_par, Wkin_perp, Wmag, Ohm, heating, source'

  write(*,'(A,20es14.6)')  ' sum ',xt,density_tot,pressure/1.d6,kin_par_tot/1.d6,kin_perp_tot/1.d6,mag_tot/1.d6, &
                                 Ohm_tot/1.d6,heating_in/1d6+heating_out/1.d6 ,source_in+source_out

#if (defined WITH_Neutrals) || (defined WITH_Impurities)
  write(*,'(A,4es14.6)')   ' Integrals_3D, MGI              : ', total_n_particles_inj, total_n_particles
  write(*,'(A,1e14.6,A)')  ' Radiation power                : ', total_radiation/1.d6, ' [MW]'
  write(*,'(A,1e14.6,A)')  ' Radiation power SANITY         : ', sum(total_radiation_phi)/1.d6, ' [MW]'
  if (with_neutrals) then
    write(*,'(A,1e14.6,A)') ' Ionization power              : ', total_P_ion/1.d6, ' [MW]'
  else if (with_impurities) then ! With CE assumption, it's easier to obtain the total ionization energy then get the ionization power by finite difference
    write(*,'(A,1e14.6,A)') ' Ionization energy              : ', total_E_ion/1.d6, ' [MJ]'
  endif
  if (with_TiTe)  write(*,'(A,1e14.6,A)') ' Electron-ion energy exchange    : ', total_P_ei/1.d6, ' [MW]'

  if (index_now > 1) then
    xtime_radiation(index_now) = xtime_radiation(index_now-1) + t_norm * tstep * total_radiation
  else if (index_now == 1) then
    xtime_radiation(index_now) = t_norm * tstep * total_radiation
  end if
  if (index_now > 0) then
  xtime_rad_power(index_now) = total_radiation
  end if

  if (output_prad_phi) then
    open(20,file="total_radiation_phi.dat",action="write",position="append")
    do i_phi = 1, n_plane
      write(20,'(1e14.6)') total_radiation_phi(i_phi)/1.d6
    end do
    close (20)
  end if

  if (with_neutrals) then
    if (index_now > 1) then
      xtime_E_ion(index_now) = xtime_E_ion(index_now-1) + t_norm * tstep * total_P_ion
    else if (index_now == 1) then
      xtime_E_ion(index_now) = t_norm * tstep * total_P_ion
    end if
    if (index_now > 0) then
      xtime_E_ion_power(index_now) = total_P_ion
    end if
  else if (with_impurities) then ! For CE assumption, we directly give the total ionization energy
    if (index_now > 0) then
      xtime_E_ion(index_now) = total_E_ion
    end if
  endif
  if (with_TiTe) then
    if (index_now > 0) xtime_P_ei(index_now) = total_P_ei
  endif
#endif

  do k = 1, n_var
    write(*,'(A,i3,A20,2es14.6)') ' min/max', k, trim(variable_names(k)), V_min(k), V_max(k)
  end do


  if (index_now > 0 ) then

    E_tot_t(index_now)               = E_tot 
    Wmag_tot_t(index_now)            = mag_tot 
    Ohmic_tot_t(index_now)           = ohm_tot
    viscopar_dissip_tot_t(index_now) = viscopar_dissip_tot
    friction_dissip_tot_t(index_now) = friction_dissip_tot
    thmwork_tot_t(index_now)         = thermal_work_tot 
    magwork_tot_t(index_now)         = mag_work_tot 
    Thermal_tot_t(index_now)         = pressure 
    Thermal_e_tot_t(index_now)       = pressure_e
    Thermal_i_tot_t(index_now)       = pressure_i
    Helicity_tot_t(index_now)        = helicity_tot
    Ip_tot_t(index_now)              = current_tot 
    current_t(index_now)             = current_in 
    Kin_par_tot_t(index_now)         = kin_par_tot
    Kin_perp_tot_t(index_now)        = kin_perp_tot 
    flux_Pvn_t(index_now)            = vn_p0
    flux_qpar_t(index_now)           = qn_par 
    flux_qperp_t(index_now)          = qn_perp 
    flux_kinpar_t(index_now)         = kinpar_flux 
    flux_poynting_t(index_now)       = poynting_flux 
    li3_t(index_now)                 = li3
    li3_tot_t(index_now)             = li3_tot
    viscopar_flux_t(index_now)       = viscopar_flux
    heat_src_tot_t(index_now)        = heating_tot
    heat_src_in_t(index_now)         = heating_in 
    heat_src_out_t(index_now)        = heating_out           
    part_src_tot_t(index_now)        = source_tot
    part_src_in_t(index_now)         = source_in
    part_src_out_t(index_now)        = source_out
    area_t(index_now)                = area
    volume_t(index_now)              = volume
    mag_ener_src_tot(index_now)      = mag_source_tot
    beta_n_t(index_now)              = beta_n
    beta_t_t(index_now)              = beta_t
    beta_p_t(index_now)              = beta_p
    npart_tot_t(index_now)           = neut_particles_tot 
    density_tot_t(index_now)         = density_tot
    density_in_t(index_now)          = density_in
    density_out_t(index_now)         = density_out
    pressure_in_t(index_now)         = pressure_in
    pressure_out_t(index_now)        = pressure_out     
    part_flux_Dpar_t(index_now)      = Dpar_part_flux
    part_flux_Dperp_t(index_now)     = Dperp_part_flux
    part_flux_vpar_t(index_now)      = vpar_part_flux
    part_flux_vperp_t(index_now)     = vperp_part_flux
    npart_flux_t(index_now)          = neut_part_flux 
 
    !--- Calculate time derivatives at previous step (second order accuracy)
    if (index_now > 2) then
      dt_back   = xtime(index_now-1) - xtime(index_now - 2)
      dt_now    = xtime(index_now)   - xtime(index_now - 1)
      r_dt2     = (dt_now/dt_back)**2.d0

      dE_tot_dt(index_now-1) = (E_tot_t(index_now) - r_dt2*E_tot_t(index_now-2) &
        -(1.d0-r_dt2)*E_tot_t(index_now-1))  / (dt_now + dt_back*r_dt2) / t_norm

      dWmag_tot_dt(index_now-1) = (Wmag_tot_t(index_now) - r_dt2*Wmag_tot_t(index_now-2) &
        -(1.d0-r_dt2)*Wmag_tot_t(index_now-1))  / (dt_now + dt_back*r_dt2) / t_norm

      dthermal_tot_dt(index_now-1) = (thermal_tot_t(index_now) - r_dt2*thermal_tot_t(index_now-2) &
        -(1.d0-r_dt2)*thermal_tot_t(index_now-1))  / (dt_now + dt_back*r_dt2) / t_norm

      dkinperp_tot_dt(index_now-1) = (kin_perp_tot_t(index_now) - r_dt2*kin_perp_tot_t(index_now-2) &
        -(1.d0-r_dt2)*kin_perp_tot_t(index_now-1))  / (dt_now + dt_back*r_dt2) / t_norm

      dkinpar_tot_dt(index_now-1) = (kin_par_tot_t(index_now) - r_dt2*kin_par_tot_t(index_now-2) &
        -(1.d0-r_dt2)*kin_par_tot_t(index_now-1))  / (dt_now + dt_back*r_dt2) / t_norm

      dpart_tot_dt(index_now-1) = (density_tot_t(index_now) - r_dt2*density_tot_t(index_now-2) &
        -(1.d0-r_dt2)*density_tot_t(index_now-1))  / (dt_now + dt_back*r_dt2) / t_norm

      dnpart_tot_dt(index_now-1) = (npart_tot_t(index_now) - r_dt2*npart_tot_t(index_now-2) &
        -(1.d0-r_dt2)*npart_tot_t(index_now-1))  / (dt_now + dt_back*r_dt2) / t_norm

      if (with_impurities) then ! Calculate the ionization energy change
        xtime_E_ion_power(index_now-1) = (xtime_E_ion(index_now) - r_dt2*xtime_E_ion(index_now-2) &
        -(1.d0-r_dt2)*xtime_E_ion(index_now-1))  / (dt_now + dt_back*r_dt2) / t_norm
      endif

    endif

    !--- Estimate time derivatives for 1st tstep (1st order accuracy)
    if (index_now == 2) then
      dt_now    = xtime(index_now) - xtime(index_now - 1)

      dE_tot_dt(index_now-1)       = (E_tot_t(index_now)-E_tot_t(index_now-1))               / dt_now / t_norm

      dWmag_tot_dt(index_now-1)    = (Wmag_tot_t(index_now)-Wmag_tot_t(index_now-1))         / dt_now / t_norm

      dthermal_tot_dt(index_now-1) = (thermal_tot_t(index_now)-thermal_tot_t(index_now-1))   / dt_now / t_norm

      dkinperp_tot_dt(index_now-1) = (kin_perp_tot_t(index_now)-kin_perp_tot_t(index_now-1)) / dt_now / t_norm

      dkinpar_tot_dt(index_now-1)  = (kin_par_tot_t(index_now) - kin_par_tot_t(index_now-1)) / dt_now / t_norm

      dpart_tot_dt(index_now-1)    = (density_tot_t(index_now)-density_tot_t(index_now-1))   / dt_now / t_norm

      dnpart_tot_dt(index_now-1)   = (npart_tot_t(index_now)-npart_tot_t(index_now-1))       / dt_now / t_norm

      if (with_impurities) then ! Calculate the ionization energy change
        xtime_E_ion_power(index_now-1) = (xtime_E_ion(index_now)-xtime_E_ion(index_now-1))   / dt_now / t_norm
      endif

    endif


  endif

endif !--- my_id
end subroutine int3d_new 

#ifndef NOMPIVERSION
end module mod_integrals3D
#endif
