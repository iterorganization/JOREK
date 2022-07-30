module mod_elt_matrix_fft

  implicit none

contains

subroutine element_matrix_fft(element, nodes, xpoint2, xcase2, R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint, Z_xpoint, ELM, RHS, tid, &
                              ELM_p, ELM_n, ELM_k, ELM_kn, RHS_p, RHS_k,  eq_g, eq_s, eq_t, eq_p, eq_ss, eq_st, eq_tt,               &
                              delta_g, delta_s, delta_t, i_tor_min, i_tor_max, aux_nodes, ELM_pnn)
!---------------------------------------------------------------
! calculates the matrix contribution of one element
!---------------------------------------------------------------
use constants
use mod_parameters
use data_structure
use gauss
use basis_at_gaussian
use phys_module
use pellet_module
use diffusivities, only: get_dperp, get_zk_eperp, get_zk_iperp
use corr_neg
use mod_injection_source
use mod_impurity
use mod_coronal
use mod_bootstrap_functions
use equil_info, only : get_psi_n
use mod_sources

implicit none

type (type_element)        :: element
type (type_node)           :: nodes(n_vertex_max)
type (type_node), optional :: aux_nodes(n_vertex_max)

#define DIM0 n_tor*n_vertex_max*n_degrees*n_var

real*8, dimension (DIM0,DIM0)  :: ELM
real*8, dimension (DIM0)       :: RHS
integer, intent(in)            :: tid
integer, intent(in)            :: i_tor_min, i_tor_max
integer    :: i, j, ms, mt, mp, k, l, index_ij, index_kl, index, index_k, index_m, m, ik, xcase2, i_inj, n_spi_tmp
integer    :: in, im, ij1, ij2, ij3, ij4, ij5, ij6, ij7, kl1, kl2, kl3, kl4, kl5, kl6, kl7
real*8     :: wst, xjac, xjac_s, xjac_t, xjac_x, xjac_y, BigR, r2, phi, delta_phi, eps_cyl
real*8     :: current_source(n_gauss,n_gauss), particle_source(n_gauss,n_gauss)
real*8     :: heat_source_i(n_gauss,n_gauss), heat_source_e(n_gauss,n_gauss)
real*8     :: R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint(2), Z_xpoint(2), dj_dpsi, dj_dz, source_pellet, source_volume
real*8     :: Bgrad_rho_star,     Bgrad_rho,     Bgrad_rhon,    Bgrad_T_star,  Bgrad_Ti, Bgrad_Te, BB2
real*8     :: Bgrad_rho_star_psi, Bgrad_rho_psi, Bgrad_rhon_psi, Bgrad_rho_rho, Bgrad_rho_rhon
real*8     :: Bgrad_T_star_psi, Bgrad_Ti_psi, Bgrad_Ti_T, Bgrad_Te_psi, Bgrad_Te_T, BB2_psi
real*8     :: Bgrad_rho_rho_n, Bgrad_rho_rhon_n, Bgrad_Te_T_n, Bgrad_rho_k_star, Bgrad_T_k_star
real*8     :: ZK_i_par_T, dZK_i_par_dT, ZK_e_par_T, dZK_e_par_dT, Bgrad_Ti_T_n
real*8     :: D_prof, ZK_i_prof, ZK_e_prof, psi_norm, theta, zeta, delta_u_x, delta_u_y, delta_ps_x, delta_ps_y
real*8     :: rhs_ij_1,   rhs_ij_2,   rhs_ij_3,   rhs_ij_4,   rhs_ij_5,   rhs_ij_6, rhs_ij_7
real*8     :: rhs_ij_5_k, rhs_ij_6_k, rhs_ij_7_k
real*8     :: rhs_stab_1, rhs_stab_2, rhs_stab_3, rhs_stab_4, rhs_stab_5, rhs_stab_6
real*8     :: D_prof_imp
real*8     :: v, v_x, v_y, v_s, v_t, v_p, v_ss, v_st, v_tt, v_xx, v_xy, v_yy
real*8     :: ps0, ps0_x, ps0_y, ps0_p, ps0_s, ps0_t, ps0_ss, ps0_tt, ps0_st, ps0_xx, ps0_yy, ps0_xy
real*8     :: zj0, zj0_x, zj0_y, zj0_p, zj0_s, zj0_t
real*8     :: u0, u0_x, u0_y, u0_p, u0_s, u0_t, u0_ss, u0_st, u0_tt, u0_xx, u0_xy, u0_yy 
real*8     :: w0, w0_x, w0_y, w0_p, w0_s, w0_t, w0_ss, w0_st, w0_tt, w0_xx, w0_xy, w0_yy
real*8     :: r0, r0_x, r0_y, r0_p, r0_s, r0_t, r0_ss, r0_st, r0_tt, r0_xx, r0_xy, r0_yy
real*8     :: r0_corr, rn0_corr, r0_hat, r0_x_hat, r0_y_hat
real*8     :: Ti0_corr, dTi0_corr_dT, d2Ti0_corr_dT2, Te0_corr, dTe0_corr_dT, d2Te0_corr_dT2
real*8     :: dr0_corr_dn, drn0_corr_dn
real*8     :: Ti0, Ti0_x, Ti0_y, Ti0_p, Ti0_s, Ti0_t, Ti0_ss, Ti0_st, Ti0_tt, Ti0_xx, Ti0_xy, Ti0_yy
real*8     :: Te0, Te0_x, Te0_y, Te0_p, Te0_s, Te0_t, Te0_ss, Te0_st, Te0_tt, Te0_xx, Te0_xy, Te0_yy
real*8     :: psi, psi_x, psi_y, psi_p, psi_s, psi_t, psi_ss, psi_st, psi_tt, psi_xx, psi_yy, psi_xy
real*8     :: zj, zj_x, zj_y, zj_p, zj_s, zj_t, zj_ss, zj_st, zj_tt
real*8     :: u, u_x, u_y, u_p, u_s, u_t, u_ss, u_st, u_tt, u_xx, u_xy, u_yy
real*8     :: w, w_x, w_y, w_p, w_s, w_t, w_ss, w_st, w_tt, w_xx, w_xy, w_yy
real*8     :: rho, rho_x, rho_y, rho_s, rho_t, rho_p, rho_hat, rho_x_hat, rho_y_hat, rho_ss, rho_st, rho_tt, rho_xx, rho_xy, rho_yy
real*8     :: Ti, Ti_x, Ti_y, Ti_s, Ti_t, Ti_p, Ti_ss, Ti_st, Ti_tt, Ti_xx, Ti_xy, Ti_yy
real*8     :: Te, Te_x, Te_y, Te_s, Te_t, Te_p, Te_ss, Te_st, Te_tt, Te_xx, Te_xy, Te_yy
real*8     :: zTi, zTi_x, zTi_y, zTe, zTe_x, zTe_y, zn_x, zn_y
real*8     :: Jb_0 , Jb
real*8     :: Vpar, Vpar_x, Vpar_y, Vpar_p, Vpar_s, Vpar_t, Vpar_ss, Vpar_st, Vpar_tt, Vpar_xx, Vpar_yy, Vpar_xy
real*8     :: P0, P0_s, P0_t, P0_x, P0_y, P0_p, P0_ss, P0_st, P0_tt, P0_xx, P0_xy, P0_yy
real*8     :: Pi0, Pi0_x, Pi0_y, Pi0_s, Pi0_t, Pi0_ss, Pi0_st, Pi0_tt, Pi0_p, Pi0_xx, Pi0_xy, Pi0_yy
real*8     :: Pe0, Pe0_x, Pe0_y, Pe0_s, Pe0_t, Pe0_ss, Pe0_st, Pe0_tt, Pe0_p, Pe0_xx, Pe0_xy, Pe0_yy
real*8     :: Vpar0, Vpar0_s, Vpar0_t, Vpar0_p, Vpar0_x, Vpar0_y, Vpar0_ss, Vpar0_st, Vpar0_tt, Vpar0_xx, Vpar0_yy,Vpar0_xy
real*8     :: BigR_x, vv2, eta_T, visco_T, deta_dT, d2eta_d2T, dvisco_dT, visco_num_T, eta_num_T
real*8     :: deta_dr0, deta_drn0, deta_num_dT, dvisco_num_dT
real*8     :: eta_T_ohm, deta_dT_ohm, deta_dr0_ohm, deta_drn0_ohm
real*8     :: ZK_par_num, Ti0_ps0_x, Ti_ps0_x, Ti0_psi_x, Ti0_ps0_y, Ti_ps0_y, Ti0_psi_y, v_ps0_x, v_psi_x, v_ps0_y, v_psi_y
real*8     :: Te0_ps0_x, Te_ps0_x, Te0_psi_x, Te0_ps0_y, Te_ps0_y, Te0_psi_y
real*8     :: amat_11, amat_12, amat_21, amat_22, amat_23, amat_24, amat_25, amat_26, amat_33, amat_31, amat_44, amat_42
real*8     :: amat_51, amat_51_n, amat_52, amat_55, amat_56, amat_57, amat_57_k, amat_58_k, amat_58_n
real*8     :: amat_61, amat_62, amat_63, amat_65, amat_66, amat_67, amat_65_k, amat_65_kn, amat_67_k, amat_16, amat_13
real*8     :: amat_71, amat_72, amat_75, amat_76, amat_75_n, amat_76_n, amat_15, amat_15_n, amat_16_n, amat_18
real*8     :: amat_19, amat_29
real*8     :: amat_12_n, amat_23_n, amat_51_k, amat_55_kn, amat_55_k, amat_55_n, amat_57_n, amat_57_kn, amat_58_kn
real*8     :: amat_75_k, amat_75_kn, amat_77, amat_77_k, amat_77_n, amat_77_kn
real*8     :: amat_61_k, amat_65_n, amat_66_kn, amat_66_k, amat_66_n, amat_67_n, amat_68_n, amat_68_k, amat_68_kn
real*8     :: TG_num1, TG_num2, TG_num5, TG_num6, TG_num7, TG_num8, TG_num9
logical    :: xpoint2

!==================MB: velocity profile is kept by a source which compensating diffusion
real*8     :: Vt0,Vt0_x,Vt0_y
real*8     :: V_source(n_gauss,n_gauss)
real*8     :: dV_dpsi_source(n_gauss,n_gauss),dV_dz_source(n_gauss,n_gauss)
real*8     :: dV_dpsi2,dV_dz2,dV_dpsi_dz,dV_dpsi3,dV_dpsi_dz2,dV_dpsi2_dz
!=======================================
real*8     :: eq_zne(n_gauss,n_gauss), eq_zTe(n_gauss,n_gauss), eq_zTi(n_gauss,n_gauss)
real*8     :: dn_dpsi(n_gauss,n_gauss),dn_dz,dn_dpsi2,dn_dz2,dn_dpsi_dz,dn_dpsi3,dn_dpsi_dz2,dn_dpsi2_dz
real*8     :: dTi_dpsi(n_gauss,n_gauss),dTi_dz,dTi_dpsi2,dTi_dz2,dTi_dpsi_dz,dTi_dpsi3,dTi_dpsi_dz2,dTi_dpsi2_dz
real*8     :: dTe_dpsi(n_gauss,n_gauss),dTe_dz,dTe_dpsi2,dTe_dz2,dTe_dpsi_dz,dTe_dpsi3,dTe_dpsi_dz2,dTe_dpsi2_dz
real*8     :: w00_xx, w00_yy  
!======================================= NEO
real*8     :: amat_27, Btheta2
real*8     :: epsil, Btheta2_psi
real*8, dimension(n_gauss,n_gauss)    :: amu_neo_prof, aki_neo_prof
!======================================= NEO

!================== Parameters specific to model5XX
! Matrix, RHS and neutrals-related variables
real*8     :: amat_28, amat_58, amat_68, amat_78, amat_78_n, amat_69, amat_79, amat_79_n
real*8     :: amat_81, amat_81_k, amat_82, amat_85, amat_86, amat_87, amat_87_k, amat_87_n 
real*8     :: amat_88, amat_88_k, amat_88_n, amat_88_kn
real*8     :: amat_91, amat_92, amat_93, amat_95, amat_96, amat_97, amat_98, amat_99
real*8     :: amat_91_k, amat_95_n, amat_95_k, amat_95_kn, amat_97_n, amat_97_k, amat_98_n, amat_98_k, amat_98_kn
real*8     :: amat_99_n, amat_99_k, amat_99_kn 
real*8     :: rhs_ij_8, rhs_ij_8_k
real*8     :: rhs_ij_9, rhs_ij_9_k
real*8     :: ij8, kl8, ij9, kl9
real*8     :: rn0, rn0_x, rn0_y, rn0_p, rn0_s, rn0_t, rn0_ss, rn0_st, rn0_tt, rn0_hat, rn0_x_hat, rn0_y_hat
real*8     :: rhon, rhon_x, rhon_y, rhon_s, rhon_t, rhon_p, rhon_ss, rhon_st, rhon_tt, rhon_hat, rhon_x_hat, rhon_y_hat  

real*8     :: rn0_xx, rn0_yy, rn0_xy, rhon_xx, rhon_yy

! New momentum equation related rhs and amat
real*8     :: amat_25_n, amat_27_n

! Neutral source
real*8     :: source_imp, source_bg

! time normalization
real*8     :: t_norm

real*8     :: Dn0x, Dn0y, Dn0p

! Atomic physics coefficients:
integer    :: i_imp
!   -Mass ratio between main ions and impurites (m_i/m_imp)
real*8     :: m_i_over_m_imp, m_imp
!   -Mean impurity ionization state
real*8     :: Z_imp, dZ_imp_dT, d2Z_imp_dT2, T0_Zimp, alpha_Zimp, Z_eff, dZ_eff_dT, eta_coef, deta_coef_dZeff
real*8     :: dZ_eff_dr0, dZ_eff_drn0, Z_eff_imp, dZ_eff_imp_dT
!   -Coefficients related to Z_imp
real*8     :: alpha_i, dalpha_i_dT, d2alpha_i_dT2
real*8     :: alpha_e, dalpha_e_dT, d2alpha_e_dT2, alpha_e_bis, alpha_e_tri
!   -Radiation from injected impurities
real*8     :: Lrad, dLrad_dT                                  ! Radiation rate and its derivative wrt. temperature
real*8     :: Lrad_imp_bg, dLrad_imp_bg_dT                    ! Radiation rate and its derivative wrt. temperature
real*8     :: r_imp                                           ! Background impurity density in JOREK unit
real*8     :: Te_corr_eV, dTe_corr_eV_dT                      ! Temperature used in radiation rate
real*8     :: Te_eV                                           ! Uncorrected temperature
real*8     :: ne_SI                                           ! Electron density used in radiation rate
real*8     :: ne_JOREK                                        ! Electron density in JOREK unit 
real*8     :: A0_rad, A1_rad, T1_rad, sig1_rad                ! Radiation rate parameters
real*8     :: A2_rad, T2_rad, sig2_rad
!   -Radiation from background impurities
real*8     :: Arad_bg, Brad_bg, Crad_bg, frad_bg, dfrad_bg_dT

!   -Ion-electron energy transfer
real*8     :: nu_e_imp, nu_e_bg, lambda_e_imp, lambda_e_bg, dTi_e, dTe_i
real*8     :: dnu_e_imp_dTi, dnu_e_imp_dTe, dnu_e_bg_dTi, dnu_e_bg_dTe
real*8     :: dnu_e_imp_drho, dnu_e_imp_drhon, dnu_e_bg_drho, dnu_e_bg_drhon
real*8     :: ddTi_e_dTi, ddTi_e_dTe, ddTi_e_drho, ddTi_e_drhon
real*8     :: ddTe_i_dTi, ddTe_i_dTe, ddTe_i_drho, ddTe_i_drhon

!   -Temporary variable for charge state distribution
real*8, allocatable :: dP_imp_dT(:), P_imp(:)
real*8     :: E_ion, dE_ion_dT, E_ion_bg
integer*8  :: ion_i, ion_k

! Parameters related to negative temperature handling
real*8     :: T_neg, delta_neg
!===============================

real*8     :: in_fft(1:n_plane)
complex*16 :: out_fft(1:n_plane)
integer*8  :: plan


#define DIM1 n_plane
#define DIM2 1:n_vertex_max*n_var*n_degrees

real*8, dimension(DIM1, DIM2, DIM2) :: ELM_p
real*8, dimension(DIM1, DIM2, DIM2) :: ELM_n
real*8, dimension(DIM1, DIM2, DIM2) :: ELM_k
real*8, dimension(DIM1, DIM2, DIM2) :: ELM_kn
real*8, dimension(DIM1, DIM2)       :: RHS_p
real*8, dimension(DIM1, DIM2)       :: RHS_k
real*8, dimension(DIM1, DIM2, DIM2) :: ELM_pnn

real*8, dimension(n_gauss,n_gauss)    :: x_g, x_s, x_t
real*8, dimension(n_gauss,n_gauss)    :: x_ss, x_st, x_tt
real*8, dimension(n_gauss,n_gauss)    :: y_g, y_s, y_t
real*8, dimension(n_gauss,n_gauss)    :: y_ss, y_st, y_tt

real*8, dimension(n_plane,n_var,n_gauss,n_gauss) :: eq_g, eq_s, eq_t
real*8, dimension(n_plane,n_var,n_gauss,n_gauss) :: eq_p
real*8, dimension(n_plane,n_var,n_gauss,n_gauss) :: eq_ss, eq_st, eq_tt
real*8, dimension(n_plane,n_var,n_gauss,n_gauss) :: delta_g, delta_s, delta_t

ELM_p = 0.d0
ELM_n = 0.d0
ELM_k = 0.d0
ELM_kn = 0.d0
RHS_p = 0.d0
RHS_k = 0.d0
ELM   = 0.d0
RHS   = 0.d0
!======================================= NEO
epsil=1.d-3
!======================================= NEO

zk_par_num = 0.d0
amat_57_kn = 0.d0

! --- Taylor-Galerkin Stabilisation coefficients
TG_num1    = TGNUM(1); TG_num2    = TGNUM(2); TG_num5    = TGNUM(5); TG_num6    = TGNUM(6); TG_num7    = TGNUM(7);

TG_num8    = TGNUM(8); TG_num9    = TGNUM(9)

! --- Take time evolution parameters from phys_module
theta = time_evol_theta
!zeta  = time_evol_zeta
! change zeta for variable dt
zeta  = time_evol_zeta * 2.0d0 * tstep / (tstep + tstep_prev)

!---------------------------------------------------- value of (x,y) and derivatives on Gaussian points
x_g  = 0.d0; x_s  = 0.d0; x_t  = 0.d0; x_st  = 0.d0; x_ss  = 0.d0; x_tt  = 0.d0;
y_g  = 0.d0; y_s  = 0.d0; y_t  = 0.d0; y_st  = 0.d0; y_ss  = 0.d0; y_tt  = 0.d0;
eq_g = 0.d0; eq_s = 0.d0; eq_t = 0.d0; eq_st = 0.d0; eq_ss = 0.d0; eq_tt = 0.d0; eq_p = 0.d0;

delta_g = 0.d0; delta_s = 0.d0; delta_t = 0.d0

current_source  = 0.d0
particle_source = 0.d0
heat_source_e   = 0.d0
heat_source_i   = 0.d0
V_source=0.d0
dV_dpsi_source=0.d0
dV_dz_source=0.d0
eq_zne          = 0.d0
eq_zTe          = 0.d0         
eq_zTi          = 0.d0
!======================================= NEO
if ( NEO ) then 
   amu_neo_prof   = 0.d0
   aki_neo_prof   = 0.d0
endif
!======================================= NEO

do i=1,n_vertex_max
 do j=1,n_degrees

   do ms=1, n_gauss
     do mt=1, n_gauss

       x_g(ms,mt)  = x_g(ms,mt)  + nodes(i)%x(1,j,1) * element%size(i,j) * H(i,j,ms,mt)
       x_s(ms,mt)  = x_s(ms,mt)  + nodes(i)%x(1,j,1) * element%size(i,j) * H_s(i,j,ms,mt)
       x_t(ms,mt)  = x_t(ms,mt)  + nodes(i)%x(1,j,1) * element%size(i,j) * H_t(i,j,ms,mt)

       x_ss(ms,mt) = x_ss(ms,mt) + nodes(i)%x(1,j,1) * element%size(i,j) * H_ss(i,j,ms,mt)
       x_st(ms,mt) = x_st(ms,mt) + nodes(i)%x(1,j,1) * element%size(i,j) * H_st(i,j,ms,mt)
       x_tt(ms,mt) = x_tt(ms,mt) + nodes(i)%x(1,j,1) * element%size(i,j) * H_tt(i,j,ms,mt)

       y_g(ms,mt)  = y_g(ms,mt)  + nodes(i)%x(1,j,2) * element%size(i,j) * H(i,j,ms,mt)
       y_s(ms,mt)  = y_s(ms,mt)  + nodes(i)%x(1,j,2) * element%size(i,j) * H_s(i,j,ms,mt)
       y_t(ms,mt)  = y_t(ms,mt)  + nodes(i)%x(1,j,2) * element%size(i,j) * H_t(i,j,ms,mt)

       y_ss(ms,mt) = y_ss(ms,mt) + nodes(i)%x(1,j,2) * element%size(i,j) * H_ss(i,j,ms,mt)
       y_st(ms,mt) = y_st(ms,mt) + nodes(i)%x(1,j,2) * element%size(i,j) * H_st(i,j,ms,mt)
       y_tt(ms,mt) = y_tt(ms,mt) + nodes(i)%x(1,j,2) * element%size(i,j) * H_tt(i,j,ms,mt)

       do mp=1,n_plane

         do k=1,n_var

           do in=1,n_tor

             eq_g(mp,k,ms,mt) = eq_g(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H(i,j,ms,mt)  * HZ(in,mp)
             eq_s(mp,k,ms,mt) = eq_s(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_s(i,j,ms,mt)* HZ(in,mp)
             eq_t(mp,k,ms,mt) = eq_t(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_t(i,j,ms,mt)* HZ(in,mp)
             eq_p(mp,k,ms,mt) = eq_p(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H(i,j,ms,mt)  * HZ_p(in,mp)

             eq_ss(mp,k,ms,mt) = eq_ss(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_ss(i,j,ms,mt)* HZ(in,mp)
             eq_st(mp,k,ms,mt) = eq_st(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_st(i,j,ms,mt)* HZ(in,mp)
             eq_tt(mp,k,ms,mt) = eq_tt(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_tt(i,j,ms,mt)* HZ(in,mp)

             delta_g(mp,k,ms,mt) = delta_g(mp,k,ms,mt) + nodes(i)%deltas(in,j,k) * element%size(i,j) * H(i,j,ms,mt)   * HZ(in,mp)
             delta_s(mp,k,ms,mt) = delta_s(mp,k,ms,mt) + nodes(i)%deltas(in,j,k) * element%size(i,j) * H_s(i,j,ms,mt) * HZ(in,mp)
             delta_t(mp,k,ms,mt) = delta_t(mp,k,ms,mt) + nodes(i)%deltas(in,j,k) * element%size(i,j) * H_t(i,j,ms,mt) * HZ(in,mp)

           enddo

         enddo

       enddo
     enddo
   enddo
 enddo
enddo

! changes deltas for variable time steps
delta_g = delta_g * tstep / tstep_prev
delta_s = delta_s * tstep / tstep_prev
delta_t = delta_t * tstep / tstep_prev

do ms=1, n_gauss
  do mt=1, n_gauss

    if (keep_current_prof) &
      call current(xpoint2, xcase2, x_g(ms,mt),y_g(ms,mt), Z_xpoint, eq_g(1,1,ms,mt),psi_axis,psi_bnd,current_source(ms,mt))

    call sources(xpoint2, xcase2, y_g(ms,mt), Z_xpoint, eq_g(1,1,ms,mt),psi_axis,psi_bnd,particle_source(ms,mt),heat_source_i(ms,mt),heat_source_e(ms,mt))
    
    call density(xpoint2, xcase2, y_g(ms,mt), Z_xpoint, eq_g(1,1,ms,mt),psi_axis,psi_bnd,eq_zne(ms,mt), &
                 dn_dpsi(ms,mt),dn_dz,dn_dpsi2,dn_dz2,dn_dpsi_dz,dn_dpsi3,dn_dpsi_dz2, dn_dpsi2_dz)

    call temperature_i(xpoint2, xcase2, y_g(ms,mt), Z_xpoint, eq_g(1,1,ms,mt),psi_axis,psi_bnd,eq_zTi(ms,mt), &
                          dTi_dpsi(ms,mt),dTi_dz,dTi_dpsi2,dTi_dz2,dTi_dpsi_dz,dTi_dpsi3,dTi_dpsi_dz2, dTi_dpsi2_dz)
    call temperature_e(xpoint2, xcase2, y_g(ms,mt), Z_xpoint, eq_g(1,1,ms,mt),psi_axis,psi_bnd,eq_zTe(ms,mt), &
                          dTe_dpsi(ms,mt),dTe_dz,dTe_dpsi2,dTe_dz2,dTe_dpsi_dz,dTe_dpsi3,dTe_dpsi_dz2, dTe_dpsi2_dz)

  enddo
enddo

!--------------------------------------------------- sum over the Gaussian integration points
do ms=1, n_gauss

 do mt=1, n_gauss

   wst = wgauss(ms)*wgauss(mt)

   xjac    = x_s(ms,mt)*y_t(ms,mt) - x_t(ms,mt)*y_s(ms,mt)

   xjac_x  = (x_ss(ms,mt)*y_t(ms,mt)**2 - y_ss(ms,mt)*x_t(ms,mt)*y_t(ms,mt) - 2.d0*x_st(ms,mt)*y_s(ms,mt)*y_t(ms,mt)   &
	   + y_st(ms,mt)*(x_s(ms,mt)*y_t(ms,mt) + x_t(ms,mt)*y_s(ms,mt))                                               &
       + x_tt(ms,mt)*y_s(ms,mt)**2 - y_tt(ms,mt)*x_s(ms,mt)*y_s(ms,mt)) / xjac

   xjac_y  = (y_tt(ms,mt)*x_s(ms,mt)**2 - x_tt(ms,mt)*y_s(ms,mt)*x_s(ms,mt) - 2.d0*y_st(ms,mt)*x_t(ms,mt)*x_s(ms,mt)   &
	   + x_st(ms,mt)*(y_t(ms,mt)*x_s(ms,mt) + y_s(ms,mt)*x_t(ms,mt))                                               &
	   + y_ss(ms,mt)*x_t(ms,mt)**2 - x_ss(ms,mt)*y_t(ms,mt)*x_t(ms,mt)) / xjac

   BigR    = x_g(ms,mt)
   BigR_x  = 1.d0

   eps_cyl = 1.d0          ! for cylinder geometry : epscyl = eps

   do mp = 1, n_plane

     ps0    = eq_g(mp,1,ms,mt)
     ps0_x  = (   y_t(ms,mt) * eq_s(mp,1,ms,mt) - y_s(ms,mt) * eq_t(mp,1,ms,mt) ) / xjac
     ps0_y  = ( - x_t(ms,mt) * eq_s(mp,1,ms,mt) + x_s(ms,mt) * eq_t(mp,1,ms,mt) ) / xjac
     ps0_p  = eq_p(mp,1,ms,mt)
     ps0_s  = eq_s(mp,1,ms,mt)
     ps0_t  = eq_t(mp,1,ms,mt)
     ps0_ss = eq_ss(mp,1,ms,mt)
     ps0_tt = eq_tt(mp,1,ms,mt)
     ps0_st = eq_st(mp,1,ms,mt)

     u0    = eq_g(mp,2,ms,mt)
     u0_x  = (   y_t(ms,mt) * eq_s(mp,2,ms,mt) - y_s(ms,mt) * eq_t(mp,2,ms,mt) ) / xjac
     u0_y  = ( - x_t(ms,mt) * eq_s(mp,2,ms,mt) + x_s(ms,mt) * eq_t(mp,2,ms,mt) ) / xjac
     u0_p  = eq_p(mp,2,ms,mt)
     u0_s  = eq_s(mp,2,ms,mt)
     u0_t  = eq_t(mp,2,ms,mt)
     u0_ss = eq_ss(mp,2,ms,mt)
     u0_tt = eq_tt(mp,2,ms,mt)
     u0_st = eq_st(mp,2,ms,mt)

     vv2   = BigR**2 *  ( u0_x * u0_x + u0_y *u0_y  )

     zj0   = eq_g(mp,3,ms,mt)
     zj0_x = (   y_t(ms,mt) * eq_s(mp,3,ms,mt) - y_s(ms,mt) * eq_t(mp,3,ms,mt) ) / xjac
     zj0_y = ( - x_t(ms,mt) * eq_s(mp,3,ms,mt) + x_s(ms,mt) * eq_t(mp,3,ms,mt) ) / xjac
     zj0_p = eq_p(mp,3,ms,mt)
     zj0_s = eq_s(mp,3,ms,mt)
     zj0_t = eq_t(mp,3,ms,mt)

     w0    = eq_g(mp,4,ms,mt)
     w0_x  = (   y_t(ms,mt) * eq_s(mp,4,ms,mt) - y_s(ms,mt) * eq_t(mp,4,ms,mt) ) / xjac
     w0_y  = ( - x_t(ms,mt) * eq_s(mp,4,ms,mt) + x_s(ms,mt) * eq_t(mp,4,ms,mt) ) / xjac
     w0_p  = eq_p(mp,4,ms,mt)
     w0_s  = eq_s(mp,4,ms,mt)
     w0_t  = eq_t(mp,4,ms,mt)
     w0_ss = eq_ss(mp,4,ms,mt)
     w0_tt = eq_tt(mp,4,ms,mt)
     w0_st = eq_st(mp,4,ms,mt)

     r0    = eq_g(mp,5,ms,mt)
     r0_x  = (   y_t(ms,mt) * eq_s(mp,5,ms,mt) - y_s(ms,mt) * eq_t(mp,5,ms,mt) ) / xjac
     r0_y  = ( - x_t(ms,mt) * eq_s(mp,5,ms,mt) + x_s(ms,mt) * eq_t(mp,5,ms,mt) ) / xjac
     r0_p  = eq_p(mp,5,ms,mt)
     r0_s  = eq_s(mp,5,ms,mt)
     r0_t  = eq_t(mp,5,ms,mt)
     r0_ss = eq_ss(mp,5,ms,mt)
     r0_st = eq_st(mp,5,ms,mt)
     r0_tt = eq_tt(mp,5,ms,mt)

     r0_corr = corr_neg_dens(r0,(/1.d-9,1.d-5/),1.d-3) ! Correction for negative r0 ...
     dr0_corr_dn = dcorr_neg_dens_drho(r0,(/1.d-9,1.d-5/),1.d-3)

     r0_hat   = BigR**2 * r0
     r0_x_hat = 2.d0 * BigR * BigR_x  * r0 + BigR**2 * r0_x
     r0_y_hat = BigR**2 * r0_y

     rn0    = eq_g(mp,8,ms,mt)
     rn0_x  = (   y_t(ms,mt) * eq_s(mp,8,ms,mt) - y_s(ms,mt) * eq_t(mp,8,ms,mt) ) / xjac    
     rn0_y  = ( - x_t(ms,mt) * eq_s(mp,8,ms,mt) + x_s(ms,mt) * eq_t(mp,8,ms,mt) ) / xjac   
     rn0_p  = eq_p(mp,8,ms,mt)                                                             
     rn0_s  = eq_s(mp,8,ms,mt)                                                             
     rn0_t  = eq_t(mp,8,ms,mt)                                                             
     rn0_ss = eq_ss(mp,8,ms,mt)                                                            
     rn0_st = eq_st(mp,8,ms,mt)                                                            
     rn0_tt = eq_tt(mp,8,ms,mt)                                                            

     rn0_corr = corr_neg_dens(rn0,(/1.d-9,1.d-5/),1.d-3) ! Correction for negative rn0 ...
     drn0_corr_dn = dcorr_neg_dens_drho(rn0, (/ 1.d-9, 1.d-5 /),1.d-3)

     rn0_xx = (rn0_ss * y_t(ms,mt)**2 - 2.d0*rn0_st * y_s(ms,mt)*y_t(ms,mt) + rn0_tt * y_s(ms,mt)**2     &
            + rn0_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                                 &
            + rn0_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) ) )    /    xjac**2               &
            - xjac_x * (rn0_s* y_t(ms,mt) - rn0_t * y_s(ms,mt))  / xjac**2

     rn0_yy = (rn0_ss * x_t(ms,mt)**2 - 2.d0*rn0_st * x_s(ms,mt)*x_t(ms,mt) + rn0_tt * x_s(ms,mt)**2     &
            + rn0_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                                 &
            + rn0_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) ) )    /    xjac**2               &
            - xjac_y * (- rn0_s * x_t(ms,mt) + rn0_t * x_s(ms,mt) )  / xjac**2

     rn0_xy = (- rn0_ss * y_t(ms,mt)*x_t(ms,mt) - rn0_tt * x_s(ms,mt)*y_s(ms,mt) &
              + rn0_st * (y_s(ms,mt)*x_t(ms,mt)  + y_t(ms,mt)*x_s(ms,mt)  ) &
              - rn0_s  * (x_st(ms,mt)*y_t(ms,mt) - x_tt(ms,mt)*y_s(ms,mt) ) &
              - rn0_t * (x_st(ms,mt)*y_s(ms,mt)  - x_ss(ms,mt)*y_t(ms,mt) ) )  / xjac**2              &
              - xjac_x * (- rn0_s * x_t(ms,mt) + rn0_t * x_s(ms,mt) )   / xjac**2



     rn0_hat   = BigR**2 * rn0                                                        
     rn0_x_hat = 2.d0 * BigR * BigR_x  * rn0 + BigR**2 * rn0_x                             
     rn0_y_hat = BigR**2 * rn0_y                                                            

     Ti0    = eq_g(mp,6,ms,mt)
     Ti0_x  = (   y_t(ms,mt) * eq_s(mp,6,ms,mt) - y_s(ms,mt) * eq_t(mp,6,ms,mt) ) / xjac
     Ti0_y  = ( - x_t(ms,mt) * eq_s(mp,6,ms,mt) + x_s(ms,mt) * eq_t(mp,6,ms,mt) ) / xjac
     Ti0_p  = eq_p(mp,6,ms,mt)
     Ti0_s  = eq_s(mp,6,ms,mt)
     Ti0_t  = eq_t(mp,6,ms,mt)
     Ti0_ss = eq_ss(mp,6,ms,mt)
     Ti0_tt = eq_tt(mp,6,ms,mt)
     Ti0_st = eq_st(mp,6,ms,mt)

     Te0    = eq_g(mp,9,ms,mt)
     Te0_x  = (   y_t(ms,mt) * eq_s(mp,9,ms,mt) - y_s(ms,mt) * eq_t(mp,9,ms,mt) ) / xjac
     Te0_y  = ( - x_t(ms,mt) * eq_s(mp,9,ms,mt) + x_s(ms,mt) * eq_t(mp,9,ms,mt) ) / xjac
     Te0_p  = eq_p(mp,9,ms,mt)
     Te0_s  = eq_s(mp,9,ms,mt)
     Te0_t  = eq_t(mp,9,ms,mt)
     Te0_ss = eq_ss(mp,9,ms,mt)
     Te0_tt = eq_tt(mp,9,ms,mt)
     Te0_st = eq_st(mp,9,ms,mt)

     ! The corrected temperature should not go below 1eV otherwise the spline
     ! fit of effective charge becomes inaccurate
     if (T_min > Ti_1) then
       Ti0_corr = corr_neg_temp(Ti0,(/5.d-1,5.d-1/),T_min) ! For use in eta(T), visco(T), ...
       dTi0_corr_dT = dcorr_neg_temp_dT(Ti0,(/5.d-1,5.d-1/),T_min) ! Improve the correction
       d2Ti0_corr_dT2 = d2corr_neg_temp_dT2(Ti0,(/5.d-1,5.d-1/),T_min)
     else
       Ti0_corr = corr_neg_temp(Ti0,(/5.d-1,5.d-1/),Ti_1) ! For use in eta(T), visco(T), ...
       dTi0_corr_dT = dcorr_neg_temp_dT(Ti0,(/5.d-1,5.d-1/),Ti_1) ! Improve the correction
       d2Ti0_corr_dT2 = d2corr_neg_temp_dT2(Ti0,(/5.d-1,5.d-1/),Ti_1)
     end if

     if (T_min > Te_1) then
       Te0_corr = corr_neg_temp(Te0,(/5.d-1,5.d-1/),T_min) ! For use in eta(T), visco(T), ...
       dTe0_corr_dT = dcorr_neg_temp_dT(Te0,(/5.d-1,5.d-1/),T_min) ! Improve the correction
       d2Te0_corr_dT2 = d2corr_neg_temp_dT2(Te0,(/5.d-1,5.d-1/),T_min)
     else
       Te0_corr = corr_neg_temp(Te0,(/5.d-1,5.d-1/),Te_1) ! For use in eta(T), visco(T), ...
       dTe0_corr_dT = dcorr_neg_temp_dT(Te0,(/5.d-1,5.d-1/),Te_1) ! Improve the correction
       d2Te0_corr_dT2 = d2corr_neg_temp_dT2(Te0,(/5.d-1,5.d-1/),Te_1)
     end if


     Vpar0    = eq_g(mp,7,ms,mt)
     Vpar0_x  = (   y_t(ms,mt) * eq_s(mp,7,ms,mt) - y_s(ms,mt) * eq_t(mp,7,ms,mt) ) / xjac
     Vpar0_y  = ( - x_t(ms,mt) * eq_s(mp,7,ms,mt) + x_s(ms,mt) * eq_t(mp,7,ms,mt) ) / xjac
     Vpar0_p  = eq_p(mp,7,ms,mt)
     Vpar0_s  = eq_s(mp,7,ms,mt)
     Vpar0_t  = eq_t(mp,7,ms,mt)
     Vpar0_ss = eq_ss(mp,7,ms,mt)
     Vpar0_st = eq_st(mp,7,ms,mt)
     Vpar0_tt = eq_tt(mp,7,ms,mt)

     ps0_xx = (ps0_ss * y_t(ms,mt)**2 - 2.d0*ps0_st * y_s(ms,mt)*y_t(ms,mt) + ps0_tt * y_s(ms,mt)**2 &
             + ps0_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                            &
             + ps0_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) ) )    / xjac**2             &
             - xjac_x * (ps0_s* y_t(ms,mt) - ps0_t * y_s(ms,mt))  / xjac**2

     ps0_yy = (ps0_ss * x_t(ms,mt)**2 - 2.d0*ps0_st * x_s(ms,mt)*x_t(ms,mt) + ps0_tt * x_s(ms,mt)**2 &
             + ps0_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                            &
             + ps0_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) ) )    / xjac**2             &
             - xjac_y * (- ps0_s * x_t(ms,mt) + ps0_t * x_s(ms,mt) )  / xjac**2

     ps0_xy = (- ps0_ss * y_t(ms,mt)*x_t(ms,mt) - ps0_tt * x_s(ms,mt)*y_s(ms,mt)                     &
              + ps0_st * (y_s(ms,mt)*x_t(ms,mt)  + y_t(ms,mt)*x_s(ms,mt)  )                          &
              - ps0_s  * (x_st(ms,mt)*y_t(ms,mt) - x_tt(ms,mt)*y_s(ms,mt) )                          &
              - ps0_t * (x_st(ms,mt)*y_s(ms,mt)  - x_ss(ms,mt)*y_t(ms,mt) ) )  / xjac**2             &
              - xjac_x * (- ps0_s * x_t(ms,mt) + ps0_t * x_s(ms,mt) )   / xjac**2

     w0_xy = (- w0_ss * y_t(ms,mt)*x_t(ms,mt) - w0_tt * x_s(ms,mt)*y_s(ms,mt)                        &
              + w0_st * (y_s(ms,mt)*x_t(ms,mt)  + y_t(ms,mt)*x_s(ms,mt)  )                           &
              - w0_s  * (x_st(ms,mt)*y_t(ms,mt) - x_tt(ms,mt)*y_s(ms,mt) )                           &
              - w0_t  * (x_st(ms,mt)*y_s(ms,mt) - x_ss(ms,mt)*y_t(ms,mt) ) )  / xjac**2              &
              - xjac_x * (- w0_s * x_t(ms,mt) + w0_t * x_s(ms,mt) )   / xjac**2

!----------------- simplified version of 2nd derivatives (for some unknown reason this is more stable!)
     w0_xx = (w0_ss * y_t(ms,mt)**2 - 2.d0*w0_st * y_s(ms,mt)*y_t(ms,mt) + w0_tt * y_s(ms,mt)**2 )     / xjac**2
     w0_yy = (w0_ss * x_t(ms,mt)**2 - 2.d0*w0_st * x_s(ms,mt)*x_t(ms,mt) + w0_tt * x_s(ms,mt)**2 )     / xjac**2
     
     r0_xx = (r0_ss * y_t(ms,mt)**2 - 2.d0*r0_st * y_s(ms,mt)*y_t(ms,mt) + r0_tt * y_s(ms,mt)**2  &
	   + r0_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                            &
	   + r0_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) ) )    / xjac**2             & 
           - xjac_x * (r0_s* y_t(ms,mt) - r0_t * y_s(ms,mt))  / xjac**2

     r0_yy = (r0_ss * x_t(ms,mt)**2 - 2.d0*r0_st * x_s(ms,mt)*x_t(ms,mt) + r0_tt * x_s(ms,mt)**2  &
           + r0_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                            &
	   + r0_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) ) )    / xjac**2             & 
           - xjac_y * (- r0_s * x_t(ms,mt) + r0_t * x_s(ms,mt) )  / xjac**2

     r0_xy = (- r0_ss * y_t(ms,mt)*x_t(ms,mt) - r0_tt * x_s(ms,mt)*y_s(ms,mt)                     &
     	      + r0_st * (y_s(ms,mt)*x_t(ms,mt)  + y_t(ms,mt)*x_s(ms,mt)  )                        &        
              - r0_s  * (x_st(ms,mt)*y_t(ms,mt) - x_tt(ms,mt)*y_s(ms,mt) )                        &	   
	      - r0_t * (x_st(ms,mt)*y_s(ms,mt)  - x_ss(ms,mt)*y_t(ms,mt) ) )  / xjac**2           &		
           - xjac_x * (- r0_s * x_t(ms,mt) + r0_t * x_s(ms,mt) )   / xjac**2

     Ti0_xx = (Ti0_ss * y_t(ms,mt)**2 - 2.d0*Ti0_st * y_s(ms,mt)*y_t(ms,mt) + Ti0_tt * y_s(ms,mt)**2     &
              + Ti0_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                              &
              + Ti0_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) ) )    / xjac**2               &	
              - xjac_x * (Ti0_s * y_t(ms,mt) - Ti0_t * y_s(ms,mt))  / xjac**2

     Ti0_yy = (Ti0_ss * x_t(ms,mt)**2 - 2.d0*Ti0_st * x_s(ms,mt)*x_t(ms,mt) + Ti0_tt * x_s(ms,mt)**2     &
              + Ti0_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                              &
              + Ti0_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) ) )       / xjac**2            &	
              - xjac_y * (- Ti0_s * x_t(ms,mt) + Ti0_t * x_s(ms,mt) )  / xjac**2

     Ti0_xy = (- Ti0_ss * y_t(ms,mt)*x_t(ms,mt) - Ti0_tt * x_s(ms,mt)*y_s(ms,mt)                        &
              + Ti0_st * (y_s(ms,mt)*x_t(ms,mt)  + y_t(ms,mt)*x_s(ms,mt)  )                           &
              - Ti0_s  * (x_st(ms,mt)*y_t(ms,mt) - x_tt(ms,mt)*y_s(ms,mt) )                           &
              - Ti0_t  * (x_st(ms,mt)*y_s(ms,mt) - x_ss(ms,mt)*y_t(ms,mt) )  )  / xjac**2             &
              - xjac_x * (- Ti0_s * x_t(ms,mt) + Ti0_t * x_s(ms,mt) )   / xjac**2

     Te0_xx = (Te0_ss * y_t(ms,mt)**2 - 2.d0*Te0_st * y_s(ms,mt)*y_t(ms,mt) + Te0_tt * y_s(ms,mt)**2     &
              + Te0_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                              &
              + Te0_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) ) )    / xjac**2               &
              - xjac_x * (Te0_s * y_t(ms,mt) - Te0_t * y_s(ms,mt))  / xjac**2

     Te0_yy = (Te0_ss * x_t(ms,mt)**2 - 2.d0*Te0_st * x_s(ms,mt)*x_t(ms,mt) + Te0_tt * x_s(ms,mt)**2     &
              + Te0_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                              &
              + Te0_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) ) )       / xjac**2            &
              - xjac_y * (- Te0_s * x_t(ms,mt) + Te0_t * x_s(ms,mt) )  / xjac**2

     Te0_xy = (- Te0_ss * y_t(ms,mt)*x_t(ms,mt) - Te0_tt * x_s(ms,mt)*y_s(ms,mt)                        &
              + Te0_st * (y_s(ms,mt)*x_t(ms,mt)  + y_t(ms,mt)*x_s(ms,mt)  )                           &
              - Te0_s  * (x_st(ms,mt)*y_t(ms,mt) - x_tt(ms,mt)*y_s(ms,mt) )                           &
              - Te0_t  * (x_st(ms,mt)*y_s(ms,mt) - x_ss(ms,mt)*y_t(ms,mt) )  )  / xjac**2             &
              - xjac_x * (- Te0_s * x_t(ms,mt) + Te0_t * x_s(ms,mt) )   / xjac**2

     vpar0_xx = (vpar0_ss * y_t(ms,mt)**2 - 2.d0*vpar0_st * y_s(ms,mt)*y_t(ms,mt) + vpar0_tt * y_s(ms,mt)**2 &
               + vpar0_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                                &
               + vpar0_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) ) )   / xjac**2                  &
               - xjac_x * (vpar0_s * y_t(ms,mt) - vpar0_t * y_s(ms,mt)) / xjac**2

     vpar0_yy = (vpar0_ss * x_t(ms,mt)**2 - 2.d0*vpar0_st * x_s(ms,mt)*x_t(ms,mt) + vpar0_tt * x_s(ms,mt)**2 &
               + vpar0_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                                &
               + vpar0_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) ) )      / xjac**2               &
               - xjac_y * (- vpar0_s * x_t(ms,mt) + vpar0_t * x_s(ms,mt) ) / xjac**2

     vpar0_xy = (- vpar0_ss * y_t(ms,mt)*x_t(ms,mt) - vpar0_tt * x_s(ms,mt)*y_s(ms,mt)                       &
                 + vpar0_st * (y_s(ms,mt)*x_t(ms,mt)  + y_t(ms,mt)*x_s(ms,mt)  )                             &
                 - vpar0_s  * (x_st(ms,mt)*y_t(ms,mt) - x_tt(ms,mt)*y_s(ms,mt) )                             &
                 - vpar0_t * (x_st(ms,mt)*y_s(ms,mt)  - x_ss(ms,mt)*y_t(ms,mt) )  )  / xjac**2               &
                 - xjac_x * (- vpar0_s * x_t(ms,mt) + vpar0_t * x_s(ms,mt) )   / xjac**2
 
     Ti0_ps0_x = Ti0_xx * ps0_y - Ti0_xy * ps0_x + Ti0_x * ps0_xy - Ti0_y * ps0_xx
     Ti0_ps0_y = Ti0_xy * ps0_y - Ti0_yy * ps0_x + Ti0_x * ps0_yy - Ti0_y * ps0_xy

     Te0_ps0_x = Te0_xx * ps0_y - Te0_xy * ps0_x + Te0_x * ps0_xy - Te0_y * ps0_xx
     Te0_ps0_y = Te0_xy * ps0_y - Te0_yy * ps0_x + Te0_x * ps0_yy - Te0_y * ps0_xy

     u0_xx = (u0_ss * y_t(ms,mt)**2 - 2.d0*u0_st * y_s(ms,mt)*y_t(ms,mt) + u0_tt * y_s(ms,mt)**2  &
           + u0_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                               &
           + u0_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) ) )      / xjac**2              &
           - xjac_x * (u0_s * y_t(ms,mt) - u0_t * y_s(ms,mt)) / xjac**2

     u0_yy = (u0_ss * x_t(ms,mt)**2 - 2.d0*u0_st * x_s(ms,mt)*x_t(ms,mt) + u0_tt * x_s(ms,mt)**2     &
           + u0_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                               &
           + u0_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) ) )      / xjac**2              &
           - xjac_y * (- u0_s * x_t(ms,mt) + u0_t * x_s(ms,mt) ) / xjac**2

     u0_xy = (- u0_ss * y_t(ms,mt)*x_t(ms,mt) - u0_tt * x_s(ms,mt)*y_s(ms,mt)                        &
              + u0_st * (y_s(ms,mt)*x_t(ms,mt)  + y_t(ms,mt)*x_s(ms,mt)  )                           &
              - u0_s  * (x_st(ms,mt)*y_t(ms,mt) - x_tt(ms,mt)*y_s(ms,mt) )                           &
              - u0_t * (x_st(ms,mt)*y_s(ms,mt)  - x_ss(ms,mt)*y_t(ms,mt) )  )  / xjac**2             &
              - xjac_x * (- u0_s * x_t(ms,mt) + u0_t * x_s(ms,mt) )   / xjac**2

     delta_u_x = (   y_t(ms,mt) * delta_s(mp,2,ms,mt) - y_s(ms,mt) * delta_t(mp,2,ms,mt) ) / xjac
     delta_u_y = ( - x_t(ms,mt) * delta_s(mp,2,ms,mt) + x_s(ms,mt) * delta_t(mp,2,ms,mt) ) / xjac

     delta_ps_x = (   y_t(ms,mt) * delta_s(mp,1,ms,mt) - y_s(ms,mt) * delta_t(mp,1,ms,mt) ) / xjac
     delta_ps_y = ( - x_t(ms,mt) * delta_s(mp,1,ms,mt) + x_s(ms,mt) * delta_t(mp,1,ms,mt) ) / xjac
     
     ! --- Temperature dependent resistivity
     if ( eta_T_dependent .and. Te0_corr <= T_max_eta) then
       eta_T     = eta   * (Te0_corr/Te_0)**(-1.5d0)
       deta_dT   = ( - eta   * (1.5d0)  * Te0_corr**(-2.5d0) * Te_0**(1.5d0) ) * dTe0_corr_dT
       d2eta_d2T = (  eta   * (3.75d0) * Te0_corr**(-3.5d0) * Te_0**(1.5d0) ) * d2Te0_corr_dT2
       deta_dr0  = 0.
       deta_drn0 = 0.
     else if (eta_T_dependent .and. Te0_corr > T_max_eta) then
       eta_T     = eta   * (T_max_eta/Te_0)**(-1.5d0)
       deta_dT   = 0.
       d2eta_d2T = 0.
       deta_dr0  = 0.
       deta_drn0 = 0.
     else
       eta_T     = eta
       deta_dT   = 0.d0
       d2eta_d2T = 0.d0
       deta_dr0  = 0.
       deta_drn0 = 0.
     end if

     if ( eta_T_dependent .and. Te0_corr <= T_max_eta_ohm ) then
       eta_T_ohm     = eta_ohmic   * (Te0_corr/Te_0)**(-1.5d0)
       deta_dT_ohm   = ( - eta_ohmic   * (1.5d0)  * Te0_corr**(-2.5d0) * Te_0**(1.5d0) ) * dTe0_corr_dT
       deta_dr0_ohm  = 0.
       deta_drn0_ohm = 0.
     else if (eta_T_dependent .and. Te0_corr > T_max_eta_ohm) then
       eta_T_ohm     = eta_ohmic   * (T_max_eta_ohm/Te_0)**(-1.5d0)
       deta_dT_ohm   = 0.
       deta_dr0_ohm  = 0.
       deta_drn0_ohm = 0.
     else
       eta_T_ohm      = eta_ohmic
       deta_dT_ohm    = 0.d0
       deta_dr0_ohm   = 0.
       deta_drn0_ohm  = 0.
     end if

     ! --- Temperature dependent viscosity
     ! --- Note: No good physics basis, simply for keeping the magnetic Prandtl number constant
     if ( visco_T_dependent .and. Te0_corr <= T_max_eta ) then       
       visco_T   = visco * (Te0_corr/Te_0)**(-1.5d0)
       dvisco_dT = - visco * (1.5d0)  * Te0_corr**(-2.5d0) * Te_0**(1.5d0) * dTe0_corr_dT
     else if (visco_T_dependent .and. Te0_corr > T_max_eta) then
       visco_T   = visco * (T_max_eta/Te_0)**(-1.5d0)
       dvisco_dT = 0.
     else
       visco_T   = visco
       dvisco_dT = 0.d0
     end if
     
     ! --- Temperature dependent parallel heat diffusivity
     if ( ZKpar_T_dependent ) then
       ZK_i_par_T   = ZK_i_par * (Ti0_corr/Ti_0)**(+2.5d0)              ! temperature dependent parallel conductivity
       dZK_i_par_dT = ZK_i_par * (2.5d0)  * Ti0_corr**(+1.5d0) * Ti_0**(-2.5d0) * dTi0_corr_dT
       if (ZK_i_par_T .gt. ZK_par_max) then
         ZK_i_par_T   = Zk_par_max
         dZK_i_par_dT = 0.d0
       endif

       ZK_e_par_T   = ZK_e_par * (Te0_corr/Te_0)**(+2.5d0)              ! temperature dependent parallel conductivity
       dZK_e_par_dT = ZK_e_par * (2.5d0)  * Te0_corr**(+1.5d0) * Te_0**(-2.5d0) * dTe0_corr_dT
       if (ZK_e_par_T .gt. ZK_par_max) then
         ZK_e_par_T   = Zk_par_max
         dZK_e_par_dT = 0.d0
       endif
     else
       ZK_i_par_T   = ZK_i_par                                            ! parallel conductivity
       dZK_i_par_dT = 0.d0
       ZK_e_par_T   = ZK_e_par                                            ! parallel conductivity
       dZK_e_par_dT = 0.d0
     endif


     ! --- Temperature dependent hyper-resistivity resistivity, there is no
     ! physical reason for this dependence whatsoever, just to keep a constant
     ! ratio between the resistivity and hyper-resistivity
     if ( eta_num_T_dependent .and. Te0_corr <= T_max_eta) then
       eta_num_T = eta_num * (Te0_corr/Te_0)**(-1.5d0)
       deta_num_dT = ( - eta_num * (1.5d0)  * Te0_corr**(-2.5d0) * Te_0**(1.5d0) ) * dTe0_corr_dT
     else if (eta_num_T_dependent .and. Te0_corr > T_max_eta) then
       eta_num_T = eta_num * (T_max_eta/Te_0)**(-1.5d0)
       deta_num_dT = 0.
     else
       eta_num_T = eta_num
       deta_num_dT = 0.
     end if

     if ( visco_num_T_dependent .and. Te0_corr <= T_max_eta) then
       visco_num_T = visco_num * (Te0_corr/Te_0)**(-1.5d0)
       dvisco_num_dT = ( - visco_num * (1.5d0)  * Te0_corr**(-2.5d0) * Te_0**(1.5d0) ) * dTe0_corr_dT
     else if (visco_num_T_dependent .and. Te0_corr > T_max_eta) then
       visco_num_T = visco_num * (T_max_eta/Te_0)**(-1.5d0)
       dvisco_num_dT = 0.
     else
       visco_num_T = visco_num
       dvisco_num_dT = 0.
     end if


     !eta_num_T   = eta_num                         ! hyperresistivity
     !visco_num_T = visco_num                       ! hyperviscosity

     psi_norm = get_psi_n(ps0, y_g(ms,mt))

     ! --- Bootstrap current 
     if (bootstrap) then
       call bootstrap_current(bigR, y_g(ms,mt),                     &
                              R_axis,   Z_axis,   psi_axis,         &
                              R_xpoint, Z_xpoint, psi_bnd, psi_norm,&
                              ps0, ps0_x, ps0_y,                    &
                              r0,  r0_x,  r0_y,                     &
                              Ti0, Ti0_x, Ti0_y,                    &
                              Te0, Te0_x, Te0_y,                  Jb)
       
       
       ! --- Full Sauter formula for initial profiles
       
       zTi   = eq_zTi(ms,mt)              
       zTi_x = dTi_dpsi(ms,mt) * ps0_x
       zTi_y = dTi_dpsi(ms,mt) * ps0_y
       zTe   = eq_zTe(ms,mt)
       zTe_x = dTe_dpsi(ms,mt) * ps0_x
       zTe_y = dTe_dpsi(ms,mt) * ps0_x
       zn_x  = dn_dpsi(ms,mt) * ps0_x
       zn_y  = dn_dpsi(ms,mt) * ps0_y
       call bootstrap_current(bigR, y_g(ms,mt),                       &
                              R_axis,   Z_axis,   psi_axis,           &
                              R_xpoint, Z_xpoint, psi_bnd, psi_norm,  &
                              ps0, ps0_x, ps0_y,                      &
                              eq_zne(ms,mt),  zn_x,  zn_y,            &
                              zTi, zTi_x, zTi_y,                      &
                              zTe, zTe_x, zTe_y,                    Jb_0)
       ! --- Subtract the initial equilibrium part
       Jb = Jb - Jb_0
     else
       Jb = 0.d0
     endif

     D_prof     = get_dperp(psi_norm)
     if (num_d_perp_imp) then
       D_prof_imp = get_dperp(psi_norm,num_d_prof_x=num_d_perp_x_imp,&
                              num_d_prof_y=num_d_perp_y_imp,num_d_prof_len=num_d_perp_len_imp)
     else
       D_prof_imp = get_dperp(psi_norm,D_perp_sp=D_perp_imp)
     end if
     ZK_e_prof = get_zk_eperp(psi_norm)
     ZK_i_prof = get_zk_iperp(psi_norm)

     ! --- Increase diffusivity if very small density/temperature
     if (xpoint2) then
       if (r0 .lt. D_prof_neg_thresh)  then
         D_prof  = D_prof_neg
         D_prof_imp = D_prof_neg
         D_par   = D_prof_neg
         D_par_imp = D_prof_neg
       endif
       if (Ti0 .lt. ZK_i_prof_neg_thresh) then
         ZK_i_prof = ZK_i_prof_neg
       end if
       if (Ti0 .lt. ZK_i_par_neg_thresh) then
         ZK_i_par_T = ZK_i_par_neg
       endif
       if (Te0 .lt. ZK_e_prof_neg_thresh) then
         ZK_e_prof = ZK_e_prof_neg
       end if
       if (Te0 .lt. ZK_e_par_neg_thresh) then
         ZK_e_par_T = ZK_e_par_neg
       endif
     endif

     phi       = 2.d0*PI*float(mp-1)/float(n_plane) / float(n_period)
     delta_phi = 2.d0*PI/float(n_plane) / float(n_period)

     source_pellet = 0.d0
     source_volume = 0.d0

     Dn0x = D_imp_extra_R      
     Dn0y = D_imp_extra_Z      
     Dn0p = D_imp_extra_p      

     if (xpoint2) then
       if (rn0 .lt. D_imp_extra_neg_thresh)  then
        Dn0x = D_imp_extra_neg
        Dn0y = D_imp_extra_neg
        Dn0p = D_imp_extra_neg
       endif
     endif
   
  !-------------------------------------------
  ! Atomic physics parameters for Argon
  !-------------------------------------------

     select case ( trim(imp_type(index_main_imp)) )
       case('D2')
         m_i_over_m_imp = central_mass/2.
         m_imp          = 2.
       case('Ar')
         m_i_over_m_imp = central_mass/40. ! Argon mass = 40 u and main ion (D) mass = 2 u
         m_imp          = 40.
       case('Ne')
         m_i_over_m_imp = central_mass/20. ! Neon mass = 20 u and main ion (D) mass = 2 u
         m_imp          = 20.
       case default
         write(*,*) '!! Gas type "', trim(imp_type(index_main_imp)), '" unknown (in inj_source.f90) !!'
         write(*,*) '=> We assume the gas is D2.'
         m_i_over_m_imp = central_mass/2.
         m_imp          = 2.
     end select

     Z_imp = 0.
     dZ_imp_dT = 0.
     d2Z_imp_dT2 = 0.

     ! Te in eV:
     Te_corr_eV = Te0_corr/(EL_CHG*MU_ZERO*central_density*1.d20)
     dTe_corr_eV_dT = dTe0_corr_dT/(EL_CHG*MU_ZERO*central_density*1.d20)
     Te_eV = Te0/(EL_CHG*MU_ZERO*central_density*1.d20)

     ! We estimate the effective charge by a test density 10^20/m^3
     ! Later maybe we should implement a iterative method

     if (allocated(imp_adas(index_main_imp)%ionisation_energy)) then

       if (allocated(P_imp)) deallocate(P_imp)
       if (allocated(dP_imp_dT)) deallocate(dP_imp_dT)

       allocate(P_imp(0:imp_adas(index_main_imp)%n_Z))
       allocate(dP_imp_dT(0:imp_adas(index_main_imp)%n_Z))

!       call imp_cor(index_main_imp)%interp(density=20.,temperature=log10(Te_corr_eV*EL_CHG/K_BOLTZ),&
!                              p_out=P_imp,p_Te_out=dP_imp_dT,z_out=Z_imp,z_Te_out=dZ_imp_dT,&
!                              z_TeTe_out=d2Z_imp_dT2)
       call imp_cor(index_main_imp)%interp_linear(density=20.,temperature=log10(Te_corr_eV*EL_CHG/K_BOLTZ),&
                              p_out=P_imp,p_Te_out=dP_imp_dT,z_avg=Z_imp,z_avg_Te=dZ_imp_dT,&
                              z_avg_TeTe=d2Z_imp_dT2)


       ! Calculate the ionization potential energy and its derivative wrt. temperature
       E_ion     = 0.
       dE_ion_dT = 0.
       E_ion_bg  = 13.6 ! Hydrogen and deterium seem to have different ionization energy, 
                        ! but the difference is of the next order.

       do ion_i=1, imp_adas(index_main_imp)%n_Z
         do ion_k=1, ion_i
           E_ion     = E_ion + P_imp(ion_i)*imp_adas(index_main_imp)%ionisation_energy(ion_k)
           dE_ion_dT = dE_ion_dT + dP_imp_dT(ion_i)*imp_adas(index_main_imp)%ionisation_energy(ion_k)
         end do
       end do
       ! Convert from eV to JOREK unit
       E_ion     = E_ion * EL_CHG*MU_ZERO*central_density*1.d20*m_i_over_m_imp
       dE_ion_dT = dE_ion_dT * EL_CHG*MU_ZERO*central_density*1.d20*m_i_over_m_imp
       E_ion_bg  = E_ion_bg * EL_CHG*MU_ZERO*central_density*1.d20
       ! Convert the gradient in K to gradient in JOREK unit
       dE_ion_dT = dE_ion_dT * dTe_corr_eV_dT * EL_CHG / K_BOLTZ

     else

       if (allocated(P_imp)) deallocate(P_imp)
       if (allocated(dP_imp_dT)) deallocate(dP_imp_dT)

       allocate(P_imp(0:imp_adas(index_main_imp)%n_Z))
       allocate(dP_imp_dT(0:imp_adas(index_main_imp)%n_Z))

!       call imp_cor(index_main_imp)%interp(density=20.,temperature=log10(Te_corr_eV*EL_CHG/K_BOLTZ),&
!                                          z_out=Z_imp,z_Te_out=dZ_imp_dT,z_TeTe_out=d2Z_imp_dT2)
       call imp_cor(index_main_imp)%interp_linear(density=20.,temperature=log10(Te_corr_eV*EL_CHG/K_BOLTZ),&
                                     p_out=P_imp,p_Te_out=dP_imp_dT,                          &
                                     z_avg=Z_imp,z_avg_Te=dZ_imp_dT,z_avg_TeTe=d2Z_imp_dT2)

       E_ion     = 0.
       dE_ion_dT = 0.
       E_ion_bg  = 0.
     end if

     ! Convert gradient in T(K) in to gradient in T (eV)
     dZ_imp_dT = dZ_imp_dT *EL_CHG / K_BOLTZ
     ! Derivative wrt to T, with T in JOREK units
     dZ_imp_dT = dZ_imp_dT / (EL_CHG*MU_ZERO*central_density*1.d20)
     dZ_imp_dT = dZ_imp_dT * dTe0_corr_dT

     if (Te_corr_eV < 0.1) then
       Z_imp = 0.
       dZ_imp_dT = 0.
       d2Z_imp_dT2 = 0.
     endif

     if (Z_imp /= Z_imp .or. dZ_imp_dT /= dZ_imp_dT) then
       write(*,*) "WARNING!!! Z_imp:", Z_imp, dZ_imp_dT
       write(*,*) "Te_corr_eV =", Te_corr_eV
       stop
     end if

     if (dZ_imp_dT < 0) then
       write(*,*) "WARNING, ERROR with dZ_imp_dT = ", dZ_imp_dT
       write(*,*) "Z_imp, T_e", Z_imp, Te_corr_eV, Te0
       stop
     end if

     alpha_i       = m_i_over_m_imp - 1.
     dalpha_i_dT   = 0.
     d2alpha_i_dT2 = 0.

     alpha_e       = m_i_over_m_imp*Z_imp - 1.
     dalpha_e_dT   = m_i_over_m_imp*dZ_imp_dT
     d2alpha_e_dT2 = m_i_over_m_imp*d2Z_imp_dT2
     alpha_e_bis   = alpha_e + dalpha_e_dT*Te0
     alpha_e_tri   = 2. * dalpha_e_dT + d2alpha_e_dT2 * Te0

     ne_SI       = (r0_corr + alpha_e * rn0_corr) * 1.d20 * central_density ! electron density (SI)
     ne_JOREK     = r0_corr + alpha_e * rn0_corr ! Electron density in JOREK unit
     ne_JOREK     = corr_neg_dens(ne_JOREK,(/1.d-1,1.d-1/),1.d-3) ! Correction for negative electron density
                                                            ! Too small rho_1 will cause a problem
     if (ne_SI < 1.d16) ne_SI = 1.d16

     ! Calculate the effective charge of all species
     Z_eff        = 0.
     dZ_eff_dT    = 0.
     dZ_eff_dr0   = 0.
     dZ_eff_drn0  = 0.

     Z_eff_imp    = 0.
     dZ_eff_imp_dT= 0.

     ! First get the value of Z_eff
     Z_eff        = r0_corr - rn0_corr
     do ion_i=1, imp_adas(index_main_imp)%n_Z
       Z_eff      = Z_eff + m_i_over_m_imp * rn0_corr * P_imp(ion_i) * real(ion_i,8)**2
       Z_eff_imp  = Z_eff_imp + P_imp(ion_i) * real(ion_i,8)**2 ! The summation of normalized nZ**2 for impurity
       dZ_eff_imp_dT = dZ_eff_imp_dT + dP_imp_dT(ion_i) * real(ion_i,8)**2 ! Its temperature gradient
     end do
     dZ_eff_imp_dT = dZ_eff_imp_dT * dTe_corr_eV_dT * EL_CHG / K_BOLTZ ! convert from K to JOREK unit
     Z_eff         = Z_eff / ne_JOREK
     if (Z_eff < 1.) Z_eff = 1.
     if (Z_eff > (imp_adas(index_main_imp)%n_Z)**2) Z_eff = (imp_adas(index_main_imp)%n_Z)**2

     ! Then three(!) gradients
     if (Z_eff >= 1.) then
       do ion_i=1, imp_adas(index_main_imp)%n_Z
         dZ_eff_dT  = dZ_eff_dT + m_i_over_m_imp * rn0_corr * dP_imp_dT(ion_i) * real(ion_i,8)**2
       end do
       dZ_eff_dT    = dZ_eff_dT / ne_JOREK
       dZ_eff_dT    = dZ_eff_dT * dTe_corr_eV_dT * EL_CHG / K_BOLTZ ! convert from K to JOREK unit
       dZ_eff_dT    = dZ_eff_dT - Z_eff * dalpha_e_dT * rn0_corr / ne_JOREK
  
       dZ_eff_dr0   = (1. - Z_eff)/ne_JOREK
  
       dZ_eff_drn0  = dZ_eff_drn0 - 1.
       do ion_i=1, imp_adas(index_main_imp)%n_Z
         dZ_eff_drn0= dZ_eff_drn0 + m_i_over_m_imp * P_imp(ion_i) * real(ion_i,8)**2
       end do
       dZ_eff_drn0  = dZ_eff_drn0 / ne_JOREK
       dZ_eff_drn0  = dZ_eff_drn0 - Z_eff * alpha_e / ne_JOREK
     else
       Z_eff        = 1.
     end if

     ! This is to represent the dependence on Z_eff in resistivity
     eta_coef     = Z_eff*(1.+1.198*Z_eff+0.222*Z_eff**2)/(1.+2.966*Z_eff+0.753*Z_eff**2)
     eta_coef     = eta_coef / ((1.+1.198+0.222)/(1.+2.966+0.753))

     deta_coef_dZeff = (1.+1.198*Z_eff+0.222*Z_eff**2)/(1.+2.966*Z_eff+0.753*Z_eff**2)
     deta_coef_dZeff = deta_coef_dZeff + Z_eff*(1.198+2.*0.222*Z_eff)/(1.+2.966*Z_eff+0.753*Z_eff**2)
     deta_coef_dZeff = deta_coef_dZeff - Z_eff*(1.+1.198*Z_eff+0.222*Z_eff**2)*(2.966+2.*0.753*Z_eff)/((1.+2.966*Z_eff+0.753*Z_eff**2)**2)
     deta_coef_dZeff = deta_coef_dZeff / ((1.+1.198+0.222)/(1.+2.966+0.753))

     if ( eta_T_dependent ) then
       deta_dr0  = eta_T * deta_coef_dZeff * dZ_eff_dr0 * dr0_corr_dn
       deta_drn0 = eta_T * deta_coef_dZeff * dZ_eff_drn0 * drn0_corr_dn
       deta_dT   = deta_dT * eta_coef + eta_T * deta_coef_dZeff * dZ_eff_dT
       eta_T     = eta_T * eta_coef

       deta_dr0_ohm  = eta_T_ohm * deta_coef_dZeff * dZ_eff_dr0 * dr0_corr_dn
       deta_drn0_ohm = eta_T_ohm * deta_coef_dZeff * dZ_eff_drn0 * drn0_corr_dn
       deta_dT_ohm   = deta_dT_ohm * eta_coef + eta_T_ohm * deta_coef_dZeff * dZ_eff_dT
       eta_T_ohm = eta_T_ohm * eta_coef
     end if


  !-------------------------------------------
  ! --- Radiative function, using interpolation
  ! ------------------------------------------

     if (ne_SI > ne_SI_min .and. Te_eV > Te_eV_min .and. rn0 > rn0_min) then

       Lrad = 0.0
       dLrad_dT = 0.0

       ! Here we are temperarily only considering one impurity species, in the
       ! future maybe a do loop will is needed
!       call radiation_function(imp_adas(index_main_imp),imp_cor(index_main_imp),log10(ne_SI),log10(Te_corr_eV*EL_CHG/K_BOLTZ),Lrad,dLrad_dT)
       call radiation_function_linear(imp_adas(index_main_imp),imp_cor(index_main_imp),log10(ne_SI),log10(Te_corr_eV*EL_CHG/K_BOLTZ),.true.,Lrad,dLrad_dT)
       Lrad = Lrad * m_i_over_m_imp
       dLrad_dT = dLrad_dT * m_i_over_m_imp * dTe0_corr_dT            

     else

       Lrad = 0.
       dLrad_dT = 0.

       !E_ion = 0.
       !dE_ion_dT = 0.
     end if
   
     ! This is to detect N/A
     if (Lrad/=Lrad .or. dLrad_dT/=dLrad_dT .or. E_ion/=E_ion .or. dE_ion_dT/=dE_ion_dT) then
       write(*,*) "WARNING: Lrad, dLrad_dT, E_ion/=E_ion, dE_ion_dT/=dE_ion_dT = ",&
                            Lrad, dLrad_dT, E_ion, dE_ion_dT
       stop
     end if


   !--------------------------------------------------------
   ! --- Source of neutrals from Massive Gas Injection (MGI)
   !--------------------------------------------------------

     source_imp = 0.d0
     source_bg  = 0.d0

!============================================================!
! Important note: in order to implementing more complicated  !
!    model, we should add more arguments to inj_source       !
!============================================================!

     call total_imp_source(x_g(ms,mt),y_g(ms,mt),phi,ps0,source_bg,source_imp,m_i_over_m_imp,index_main_imp)

     ! This is to detect N/A
     if (source_imp /= source_imp .or. source_bg /= source_bg) then
       write(*,*) "WARNING: source_imp = ", source_imp
       write(*,*) "WARNING: source_bg = ", source_bg
       stop
     end if
     
     if (source_imp .lt. 0.d0) then
      source_imp = 0.d0
     endif
     if (source_bg .lt. 0.d0) then
      source_bg = 0.d0
     endif

   !--------------------------------------------------------
   ! --- Radiation from background impurity
   !--------------------------------------------------------

    frad_bg = 0.
    dfrad_bg_dT = 0.
    do i_imp =1, n_adas
      if (i_imp == index_main_imp) cycle
      r_imp = nimp_bg(i_imp) / (1.d20 * central_density)  ! Background impurity density in JU
      if (ne_SI > ne_SI_min .and. Te_eV > Te_eV_min .and. r_imp > 0) then
        Lrad_imp_bg = 0.0
        dLrad_imp_bg_dT = 0.0
        call radiation_function_linear(imp_adas(i_imp),imp_cor(i_imp),log10(ne_SI),    &
                                       log10(Te_corr_eV*EL_CHG/K_BOLTZ),.true.,Lrad_imp_bg,dLrad_imp_bg_dT)
        dLrad_imp_bg_dT = dLrad_imp_bg_dT * dTe0_corr_dT
      else
        Lrad_imp_bg = 0.
        dLrad_imp_bg_dT = 0.
      end if
      if (dLrad_imp_bg_dT/=dLrad_imp_bg_dT) then
        write(*,*) "WARNING: dLrad_imp_bg_dT ", dLrad_imp_bg_dT
        stop
      end if

      frad_bg = frad_bg + r_imp * Lrad_imp_bg
      dfrad_bg_dT =  dfrad_bg_dT + r_imp * dLrad_imp_bg_dT

    end do

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
                   / (1.d3*(MASS_ELECTRON*Ti0_corr+Te0_corr*MASS_PROTON*m_imp)&
                   / (EL_CHG * MU_ZERO * central_density * 1.d20)) ** 1.5
    nu_e_bg      = 1.8d-19*(1.d6*MASS_ELECTRON*MASS_PROTON*central_mass) ** 0.5&
                   * (1.d14*central_density*(r0_corr-rn0_corr)) * lambda_e_bg &
                   / (1.d3*(MASS_ELECTRON*Ti0_corr+Te0_corr*MASS_PROTON*central_mass)&
                   / (EL_CHG * MU_ZERO * central_density * 1.d20)) ** 1.5 ! Assuming bg_charge is 1!

    if (nu_e_imp < 0.) nu_e_imp = 0.
    if (nu_e_bg < 0.)  nu_e_bg  = 0.

    !Converting the energy transfer rate from s^-1 to JOREK unit
    t_norm   = sqrt(MU_ZERO * central_mass * MASS_PROTON * central_density * 1.d20)
    nu_e_imp = nu_e_imp * t_norm
    nu_e_bg  = nu_e_bg * t_norm

    dTe_i    = (nu_e_imp + nu_e_bg) * (Ti0_corr - Te0_corr) * (r0_corr + alpha_e*rn0_corr)

    !Calculating the density and temperature derivative for amats
    !We negelect the coulomb log's dericatives due to their smallness
    dnu_e_imp_dTi   = -1.5*MASS_ELECTRON*nu_e_imp*dTi0_corr_dT / (MASS_ELECTRON*Ti0_corr + MASS_PROTON*m_imp*Te0_corr)
    dnu_e_imp_dTe   = -1.5*MASS_PROTON*m_imp*nu_e_imp*dTe0_corr_dT / (MASS_ELECTRON*Ti0_corr + MASS_PROTON*m_imp*Te0_corr) &
                      + nu_e_imp * dZ_eff_imp_dT / Z_eff_imp

    dnu_e_imp_drhon = nu_e_imp * drn0_corr_dn / rn0_corr
    dnu_e_imp_drho  = 0.

    dnu_e_bg_dTi    = -1.5*MASS_ELECTRON*nu_e_bg*dTi0_corr_dT / (MASS_ELECTRON*Ti0_corr + MASS_PROTON*central_mass*Te0_corr)
    dnu_e_bg_dTe    = -1.5*MASS_PROTON*central_mass*nu_e_bg*dTe0_corr_dT &
                      / (MASS_ELECTRON*Ti0_corr + MASS_PROTON*central_mass*Te0_corr)
    if (r0_corr-rn0_corr <= 0.) then
      dnu_e_bg_drhon = 0.
      dnu_e_bg_drho  = 0.
    else
      dnu_e_bg_drhon = -nu_e_bg * drn0_corr_dn / (r0_corr-rn0_corr)
      dnu_e_bg_drho  = nu_e_bg * dr0_corr_dn / (r0_corr-rn0_corr)
    end if

    ddTe_i_dTi      = (dnu_e_imp_dTi + dnu_e_bg_dTi) * (Ti0_corr - Te0_corr) * (r0_corr + alpha_e*rn0_corr)&
                      + (nu_e_imp + nu_e_bg) * dTi0_corr_dT * (r0_corr + alpha_e*rn0_corr)
    ddTe_i_dTe      = (dnu_e_imp_dTe + dnu_e_bg_dTe) * (Ti0_corr - Te0_corr) * (r0_corr + alpha_e*rn0_corr)&
                      - (nu_e_imp + nu_e_bg) * dTe0_corr_dT * (r0_corr + alpha_e*rn0_corr)                 &
                      + (nu_e_imp + nu_e_bg) * (Ti0_corr - Te0_corr) * dalpha_e_dT * rn0_corr
    ddTe_i_drhon    = (dnu_e_imp_drhon + dnu_e_bg_drhon) * (Ti0_corr - Te0_corr) * (r0_corr + alpha_e*rn0_corr)&
                      +(nu_e_imp + nu_e_bg) * (Ti0_corr - Te0_corr) * alpha_e * drn0_corr_dn
    ddTe_i_drho     = (dnu_e_imp_drho + dnu_e_bg_drho) * (Ti0_corr - Te0_corr) * (r0_corr + alpha_e*rn0_corr)&
                      +(nu_e_imp + nu_e_bg) * (Ti0_corr - Te0_corr) * dr0_corr_dn

    if (r0_corr+alpha_e*rn0_corr < 0.) then
      dTe_i         = 0.
      ddTe_i_dTi    = 0.
      ddTe_i_dTe    = 0.
      ddTe_i_drhon  = 0.
      ddTe_i_drho   = 0.
    end if

    dTi_e           = -dTe_i
    ddTi_e_dTi      = -ddTe_i_dTi
    ddTi_e_dTe      = -ddTe_i_dTe
    ddTi_e_drhon    = -ddTe_i_drhon
    ddTi_e_drho     = -ddTe_i_drho

!--------------------------------------------------------

     Pi0    = (r0+rn0*alpha_i) * Ti0
     Pi0_x  = (r0_x+rn0_x*alpha_i) * Ti0 + (r0+rn0*alpha_i) * Ti0_x
     Pi0_y  = (r0_y+rn0_y*alpha_i) * Ti0 + (r0+rn0*alpha_i) * Ti0_y
     Pi0_s  = (r0_s+rn0_s*alpha_i) * Ti0 + (r0+rn0*alpha_i) * Ti0_s
     Pi0_t  = (r0_t+rn0_t*alpha_i) * Ti0 + (r0+rn0*alpha_i) * Ti0_t
     Pi0_p  = (r0_p+rn0_p*alpha_i) * Ti0 + (r0+rn0*alpha_i) * Ti0_p
     Pi0_ss = (r0_ss+rn0_ss*alpha_i) * Ti0 + 2.d0 * (r0_s+rn0_s*alpha_i) * Ti0_s + (r0+rn0*alpha_i) * Ti0_ss
     Pi0_tt = (r0_tt+rn0_tt*alpha_i) * Ti0 + 2.d0 * (r0_t+rn0_t*alpha_i) * Ti0_t + (r0+rn0*alpha_i) * Ti0_tt
     Pi0_st = (r0_st+rn0_st*alpha_i) * Ti0 + (r0_t+rn0_t*alpha_i) * Ti0_s + (r0_s+rn0_s*alpha_i) * Ti0_t               &
              + (r0+rn0*alpha_i) * Ti0_st
     Pi0_xx = (r0_xx+rn0_xx*alpha_i) * Ti0 + 2.d0 * (r0_x+rn0_x*alpha_i) * Ti0_x + (r0+rn0*alpha_i) * Ti0_xx
     Pi0_yy = (r0_yy+rn0_yy*alpha_i) * Ti0 + 2.d0 * (r0_y+rn0_y*alpha_i) * Ti0_y + (r0+rn0*alpha_i) * Ti0_yy
     Pi0_xy = (r0_xy+rn0_xy*alpha_i) * Ti0 + (r0_y+rn0_y*alpha_i) * Ti0_x + (r0_x+rn0_x*alpha_i) * Ti0_y               &
              + (r0+rn0*alpha_i) * Ti0_xy

     Pe0    = (r0+rn0*alpha_e) * Te0
     Pe0_x  = (r0_x+rn0_x*alpha_e) * Te0 + (r0+rn0*alpha_e_bis) * Te0_x
     Pe0_y  = (r0_y+rn0_y*alpha_e) * Te0 + (r0+rn0*alpha_e_bis) * Te0_y
     Pe0_s  = (r0_s+rn0_s*alpha_e) * Te0 + (r0+rn0*alpha_e_bis) * Te0_s
     Pe0_t  = (r0_t+rn0_t*alpha_e) * Te0 + (r0+rn0*alpha_e_bis) * Te0_t
     Pe0_p  = (r0_p+rn0_p*alpha_e) * Te0 + (r0+rn0*alpha_e_bis) * Te0_p
     Pe0_ss = (r0_ss+rn0_ss*alpha_e) * Te0 + 2.d0 * (r0_s+rn0_s*alpha_e_bis) * Te0_s + (r0+rn0*alpha_e_bis) * Te0_ss &
              + rn0 * (2.d0*dalpha_e_dT + d2alpha_e_dT2*Te0) * (Te0_s)**2.d0
     Pe0_tt = (r0_tt+rn0_tt*alpha_e) * Te0 + 2.d0 * (r0_t+rn0_t*alpha_e_bis) * Te0_t + (r0+rn0*alpha_e_bis) * Te0_tt &
              + rn0 * (2.d0*dalpha_e_dT + d2alpha_e_dT2*Te0) * (Te0_t)**2.d0
     Pe0_st = (r0_st+rn0_st*alpha_e) * Te0 + (r0_t+rn0_t*alpha_e_bis) * Te0_s + (r0_s+rn0_s*alpha_e_bis) * Te0_t     &
              + rn0 * (2.d0*dalpha_e_dT + d2alpha_e_dT2*Te0) * Te0_s * Te0_t + (r0+rn0*alpha_e_bis) * Te0_st
     Pe0_xx = (r0_xx+rn0_xx*alpha_e) * Te0 + 2.d0 * (r0_x+rn0_x*alpha_e_bis) * Te0_x + (r0+rn0*alpha_e_bis) * Te0_xx &
              + rn0 * (2.d0*dalpha_e_dT + d2alpha_e_dT2*Te0) * (Te0_x)**2.d0
     Pe0_yy = (r0_yy+rn0_yy*alpha_e) * Te0 + 2.d0 * (r0_y+rn0_y*alpha_e_bis) * Te0_y + (r0+rn0*alpha_e_bis) * Te0_yy &
              + rn0 * (2.d0*dalpha_e_dT + d2alpha_e_dT2*Te0) * (Te0_y)**2.d0
     Pe0_xy = (r0_xy+rn0_xy*alpha_e) * Te0 + (r0_y+rn0_y*alpha_e_bis) * Te0_x + (r0_x+rn0_x*alpha_e_bis) * Te0_y     &
              + rn0 * (2.d0*dalpha_e_dT + d2alpha_e_dT2*Te0) * Te0_x * Te0_y + (r0+rn0*alpha_e_bis) * Te0_xy

     P0     = Pi0 + Pe0
     P0_x   = Pi0_x + Pe0_x
     P0_y   = Pi0_y + Pe0_y
     P0_s   = Pi0_s + Pe0_s
     P0_t   = Pi0_t + Pe0_t
     P0_p   = Pi0_p + Pe0_p
     P0_ss  = Pi0_ss + Pe0_ss
     P0_tt  = Pi0_tt + Pe0_tt
     P0_st  = Pi0_st + Pe0_st
     P0_xx  = Pi0_xx + Pe0_xx
     P0_yy  = Pi0_yy + Pe0_yy
     P0_xy  = Pi0_xy + Pe0_xy

!-------------------------------------------------------- 


     do i=1,n_vertex_max

       do j=1,n_degrees

         index_ij = n_var*n_degrees*(i-1) + n_var * (j-1) + 1   ! index in the ELM matrix

         v   =  H(i,j,ms,mt) * element%size(i,j)
         v_x = (  y_t(ms,mt) * h_s(i,j,ms,mt) - y_s(ms,mt) * h_t(i,j,ms,mt) ) * element%size(i,j) / xjac
         v_y = (- x_t(ms,mt) * h_s(i,j,ms,mt) + x_s(ms,mt) * h_t(i,j,ms,mt) ) * element%size(i,j) / xjac

         v_s = h_s(i,j,ms,mt) * element%size(i,j)
         v_t = h_t(i,j,ms,mt) * element%size(i,j)
         v_p = H(i,j,ms,mt)   * element%size(i,j)

         v_ss = h_ss(i,j,ms,mt) * element%size(i,j)
         v_tt = h_tt(i,j,ms,mt) * element%size(i,j)
         v_st = h_st(i,j,ms,mt) * element%size(i,j)

         v_xx = (v_ss * y_t(ms,mt)**2 - 2.d0*v_st * y_s(ms,mt)*y_t(ms,mt) + v_tt * y_s(ms,mt)**2    &	        
	        + v_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                          &	   
	        + v_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) ) )  / xjac**2             &		
		- xjac_x * (v_s * y_t(ms,mt) - v_t * y_s(ms,mt)) / xjac**2

         v_yy = (v_ss * x_t(ms,mt)**2 - 2.d0*v_st * x_s(ms,mt)*x_t(ms,mt) + v_tt * x_s(ms,mt)**2    &	        
		+ v_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                          &	   
	        + v_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) ) )     / xjac**2          &		
		- xjac_y * (- v_s * x_t(ms,mt) + v_t * x_s(ms,mt) ) / xjac**2

         v_xy = (- v_ss * y_t(ms,mt)*x_t(ms,mt) - v_tt * x_s(ms,mt)*y_s(ms,mt)                      &
     	        + v_st * (y_s(ms,mt)*x_t(ms,mt)  + y_t(ms,mt)*x_s(ms,mt)  )                         &        
                - v_s  * (x_st(ms,mt)*y_t(ms,mt) - x_tt(ms,mt)*y_s(ms,mt) )                         &	   
	        - v_t  * (x_st(ms,mt)*y_s(ms,mt) - x_ss(ms,mt)*y_t(ms,mt) )  )  / xjac**2           &		
                - xjac_x * (- v_s * x_t(ms,mt) + v_t * x_s(ms,mt) )   / xjac**2

         Bgrad_rho_star   = ( v_x  * ps0_y - v_y  * ps0_x ) / BigR                         
         Bgrad_rho_k_star = ( F0 / BigR * v_p )           / BigR                           
         Bgrad_rho        = ( F0 / BigR * r0_p +  r0_x * ps0_y - r0_y * ps0_x ) / BigR
         Bgrad_rhon       = ( F0 / BigR * rn0_p + rn0_x * ps0_y - rn0_y * ps0_x ) / BigR

         Bgrad_T_star     = ( v_x  * ps0_y - v_y  * ps0_x ) / BigR                         
         Bgrad_T_k_star   = ( F0 / BigR * v_p           ) / BigR                           

         Bgrad_Ti       = ( F0 / BigR * Ti0_p +  Ti0_x * ps0_y - Ti0_y * ps0_x ) / BigR 
         Bgrad_Te       = ( F0 / BigR * Te0_p +  Te0_x * ps0_y - Te0_y * ps0_x ) / BigR 

         BB2            = (F0*F0 + ps0_x * ps0_x + ps0_y * ps0_y )/BigR**2
         Btheta2        = (ps0_x * ps0_x + ps0_y * ps0_y )/BigR**2

         v_ps0_x  = v_xx  * ps0_y - v_xy  * ps0_x + v_x  * ps0_xy - v_y * ps0_xx
         v_ps0_y  = v_xy  * ps0_y - v_yy  * ps0_x + v_x  * ps0_yy - v_y * ps0_xy

!###################################################################################################
!#  equation 1   (induction equation)                                                              #
!###################################################################################################


           rhs_ij_1 =   v * (eta_T  * (zj0-current_source(ms,mt)-Jb))/ BigR  * xjac * tstep &
                      + v * (ps0_s * u0_t - ps0_t * u0_s)                        * tstep &
                      - v * eps_cyl * F0 / BigR  * u0_p                   * xjac * tstep &
                      + eta_num_T * (v_x * zj0_x + v_y * zj0_y)           * xjac * tstep &

                      - v * tauIC/(r0_corr*BB2) * F0**2/BigR**2 * (ps0_s * pi0_t - ps0_t * pi0_s) * tstep &
                      + v * tauIC/(r0_corr*BB2) * F0**3/BigR**3 * eps_cyl * pi0_p * xjac * tstep &

                      + zeta * v * delta_g(mp,1,ms,mt) / BigR             * xjac


!###################################################################################################
!#  equation 2   (perpendicular momentum equation)                                                 #
!###################################################################################################

         rhs_ij_2 = - 0.5d0 * vv2 * (v_x * r0_y_hat - v_y * r0_x_hat)     * xjac * tstep &
                      - r0_hat * BigR**2 * w0 * (v_s * u0_t - v_t * u0_s)        * tstep &
                      + v * (ps0_s * zj0_t - ps0_t * zj0_s )                     * tstep &
                      - visco_T * BigR * (v_x * w0_x + v_y * w0_y)        * xjac * tstep &
                      - v * eps_cyl * F0 / BigR * zj0_p                   * xjac * tstep &
                      + BigR**2 * (v_s * p0_t - v_t * p0_s)                      * tstep &

                      - visco_num_T * (v_xx + v_x/Bigr + v_yy)*(w0_xx + w0_x/Bigr + w0_yy) * xjac * tstep &

                      - TG_num2 * 0.25d0 * r0_hat * BigR**3 * (w0_x * u0_y - w0_y * u0_x) &
                                * ( v_x * u0_y - v_y * u0_x) * xjac * tstep * tstep &

!====================================New TG_num terms=================================
                      - TG_num2 * 0.25d0 * w0 * BigR**3 * (r0_x_hat * u0_y - r0_y_hat * u0_x) &
                                * ( v_x * u0_y - v_y * u0_x) * xjac * tstep * tstep &
!===============================End of NewTG_num terms==============================

                      - v * tauIC * BigR**4 * (pi0_s * w0_t - pi0_t * w0_s)        * tstep &

                      - tauIC * BigR**3 * pi0_y * (v_x* u0_x + v_y * u0_y) * xjac * tstep &

                      - v * tauIC * BigR**4 * (u0_xy * (pi0_xx - pi0_yy) - pi0_xy * (u0_xx - u0_yy) ) * xjac * tstep &

                      - zeta * BigR * r0_hat * (v_x * delta_u_x + v_y * delta_u_y) * xjac &

                      ! New terms coming from -(\partial_t \rho + \nabla \cdot (\rho \mathbf{v})) \mathbf{v} in RHS of momentum equation
                      ! (see wiki: https://www.jorek.eu/wiki/doku.php?id=model500_501_555#equations):		      
                      - zeta * BigR * BigR**2 * delta_g(mp,5,ms,mt) * (v_x * u0_x + v_y * u0_y) * xjac &		      
                      - BigR**2 * (r0_x_hat * u0_y - r0_y_hat * u0_x) * (v_x * u0_x + v_y * u0_y)      * xjac * tstep &
                      + BigR * F0 * (r0 * vpar0_p + vpar0 * r0_p) * (v_x * u0_x + v_y * u0_y)          * xjac * tstep &		      		      
                      + BigR**2 * r0 * (vpar0_x * ps0_y - vpar0_y * ps0_x) * (v_x * u0_x + v_y * u0_y) * xjac * tstep  &
                      + BigR**2 * vpar0 * (r0_x * ps0_y - r0_y * ps0_x)    * (v_x * u0_x + v_y * u0_y) * xjac * tstep

                      ! Old term (not to be included anymore due to implementations of terms above):
                      ! + BigR**3 * (particle_source(ms,mt) + source_pellet) * (v_x * u0_x + v_y * u0_y) * xjac * tstep

!###################################################################################################
!#  equation 3   (current definition)                                                              #
!###################################################################################################

         rhs_ij_3 = - ( v_x * ps0_x  + v_y * ps0_y + v*zj0) / BigR * xjac

!###################################################################################################
!#  equation 4   (vorticity definition)                                                            #
!###################################################################################################

         rhs_ij_4 = - ( v_x * u0_x   + v_y * u0_y  + v*w0)  * BigR * xjac 

!###################################################################################################
!#  equation 5 (total density equation)                                                                  #
!###################################################################################################

         rhs_ij_5   = v * BigR * (particle_source(ms,mt) + source_pellet+source_bg+source_imp) * xjac * tstep &
                    + v * BigR**2 * ( r0_s * u0_t - r0_t * u0_s)                                      * tstep &
                    + v * 2.d0 * BigR * r0 * u0_y                                              * xjac * tstep &
                    - (D_par-D_prof) * BigR / BB2 * Bgrad_rho_star * (Bgrad_rho-Bgrad_rhon)    * xjac * tstep &
                    - (D_par_imp-D_prof_imp) * BigR / BB2 * Bgrad_rho_star * (Bgrad_rhon)      * xjac * tstep &
                    - D_prof * BigR  * (v_x*(r0_x-rn0_x) + v_y*(r0_y-rn0_y)                  ) * xjac * tstep &
                    - D_prof_imp * BigR  * (v_x*(rn0_x) + v_y*(rn0_y)                        ) * xjac * tstep &
                    ! The old diffusion scheme for the impurities
                    !- (D_par-D_prof) * BigR / BB2 * Bgrad_rho_star * (Bgrad_rho)               * xjac * tstep &
                    ! The old diffusion scheme for the impurities
                    !- D_prof * BigR  * (v_x*(r0_x) + v_y*(r0_y)             )                  * xjac * tstep &
                    + BigR* (- Dn0x * rn0_x * v_x - Dn0y * rn0_y * v_y)                        * xjac * tstep &  
                    - v * F0 / BigR * Vpar0 * r0_p                                             * xjac * tstep &
                    - v * Vpar0 * (r0_s * ps0_t - r0_t * ps0_s)                                       * tstep &
                    - v * F0 / BigR * r0 * vpar0_p                                             * xjac * tstep &
                    - v * r0 * (vpar0_s * ps0_t - vpar0_t * ps0_s)                                    * tstep &

                    + v * 2.d0 * tauIC * pi0_y * BigR                                          * xjac * tstep &

                    + zeta * v * delta_g(mp,5,ms,mt) * BigR                                    * xjac  &

                    - D_perp_num * (v_xx + v_x/Bigr + v_yy)*(r0_xx + r0_x/Bigr + r0_yy) * BigR * xjac * tstep &

                    - TG_num5 * 0.25d0 * BigR**3 * (r0_x * u0_y - r0_y * u0_x)                                &
                                                 * ( v_x * u0_y - v_y * u0_x) * xjac * tstep * tstep          &
                    - TG_num5 * 0.25d0 / BigR * vpar0**2                                                      &
                              * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                              &
                              * ( v_x * ps0_y -  v_y * ps0_x                   ) * xjac * tstep * tstep       
		    
 

         rhs_ij_5_k =  &
                      - (D_par-D_prof) * BigR / BB2 * Bgrad_rho_k_star * (Bgrad_rho-Bgrad_rhon)   * xjac * tstep &
                      - (D_par_imp-D_prof_imp) * BigR / BB2 * Bgrad_rho_k_star * (Bgrad_rhon)     * xjac * tstep &
                      - D_prof * BigR  * (          v_p*(r0_p-rn0_p) * eps_cyl**2 /BigR**2 )      * xjac * tstep &
                      - D_prof_imp * BigR  * (           v_p*(rn0_p) * eps_cyl**2 /BigR**2 )      * xjac * tstep &
                       !Old diffusion scheme for impurities
                       !- (D_par-D_prof) * BigR / BB2 * Bgrad_rho_k_star * (Bgrad_rho)             * xjac * tstep &
                       !- D_prof * BigR  * (          v_p*(r0_p) * eps_cyl**2 /BigR**2 )           * xjac * tstep &
                       + BigR* ( - Dn0p * rn0_p * v_p*eps_cyl**2/BigR**2)                         * xjac * tstep & 
                    - TG_num5 * 0.25d0 / BigR * vpar0**2 &
                              * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                              &
                              * (                            + F0 / BigR * v_p) * xjac * tstep * tstep

!###################################################################################################
!#  equation 6 (ion energy  equation)                                                              #
!###################################################################################################

         rhs_ij_6 =   v * BigR * heat_source_i(ms,mt)                                  * xjac * tstep &

                    + v * (r0 + rn0*alpha_i) * BigR**2 * ( Ti0_s * u0_t - Ti0_t * u0_s)       * tstep &
                    + v * Ti0 * BigR**2 * (r0_s * u0_t - r0_t * u0_s)                         * tstep &
                    + v * alpha_i * Ti0 * BigR**2 * (rn0_s * u0_t - rn0_t * u0_s)             * tstep &

                    + v * (r0 + rn0*alpha_i) * Ti0 * 2.d0* GAMMA * BigR * u0_y         * xjac * tstep &

                    - v * (r0 + rn0*alpha_i) * F0 / BigR * Vpar0 * Ti0_p               * xjac * tstep &
                    - v * Ti0 * F0 / BigR * Vpar0 * (r0_p + alpha_i * rn0_p)           * xjac * tstep &

                    - v * (r0 + rn0*alpha_i) * Vpar0 * (Ti0_s * ps0_t - Ti0_t * ps0_s)        * tstep &
                    - v * Ti0 * Vpar0 * (r0_s * ps0_t - r0_t * ps0_s)                         * tstep &
                    - v * Ti0 * Vpar0 * alpha_i * (rn0_s * ps0_t - rn0_t * ps0_s)             * tstep &

                    - v * (r0 + rn0*alpha_i) * Ti0 * GAMMA * (vpar0_s * ps0_t - vpar0_t * ps0_s)        * tstep &
                    - v * (r0 + rn0*alpha_i) * Ti0 * GAMMA * F0 / BigR * vpar0_p                 * xjac * tstep &

                    - (ZK_i_par_T-ZK_i_prof) * BigR / BB2 * Bgrad_T_star * Bgrad_Ti    * xjac * tstep &
                    - ZK_i_prof * BigR * (v_x*Ti0_x + v_y*Ti0_y                      ) * xjac * tstep &

                    - ZK_perp_num  *  (v_xx + v_x/Bigr + v_yy)*(Ti0_xx + Ti0_x/Bigr + Ti0_yy) * BigR * xjac * tstep &

                    - ZK_par_num * (v_ps0_x  * ps0_y - v_ps0_y  * ps0_x)                              &
                                 * (Ti0_ps0_x * ps0_y - Ti0_ps0_y * ps0_x)             * xjac * tstep &

                    - TG_num6 * 0.25d0 * BigR**3 * Ti0 * ((r0_x+alpha_i*rn0_x) * u0_y - (r0_y+alpha_i*rn0_y) * u0_x) &
                                       * ( v_x * u0_y - v_y * u0_x) * xjac * tstep * tstep  &
                    - TG_num6 * 0.25d0 * BigR**3 * (r0+alpha_i*rn0) * (Ti0_x * u0_y - Ti0_y * u0_x) &
                                       * ( v_x * u0_y - v_y * u0_x) * xjac * tstep * tstep  &

                    - TG_num6 * 0.25d0 / BigR * vpar0**2 &
                              * Ti0 * ((r0_x+alpha_i*rn0_x) * ps0_y - (r0_y+alpha_i*rn0_y) * ps0_x &
                                        + F0 / BigR * (r0_p+alpha_i*rn0_p))                        &
                              * ( v_x * ps0_y -  v_y * ps0_x                  ) * xjac * tstep * tstep        &

                    - TG_num6 * 0.25d0 / BigR * vpar0**2 &
                              * (r0+alpha_i*rn0) * (Ti0_x * ps0_y - Ti0_y * ps0_x + F0 / BigR * Ti0_p)        &
                              * ( v_x * ps0_y -  v_y * ps0_x                  ) * xjac * tstep * tstep        &

                    + zeta * v * (r0 + rn0 * alpha_i) * delta_g(mp,6,ms,mt) * BigR                     * xjac &
                    + zeta * v * Ti0 * delta_g(mp,5,ms,mt) * BigR                                      * xjac &
                    + zeta * v * alpha_i * Ti0 * delta_g(mp,8,ms,mt) * BigR                            * xjac &
!============================Behold, the parallel viscous heating terms!=============
                    + (GAMMA - 1.) * v * BigR * visco_par_heating * (vpar0_x * vpar0_x + vpar0_y * vpar0_y) * xjac * tstep &
                    + (GAMMA - 1.) * vpar0 * BigR * visco_par_heating * (v_x * vpar0_x     + v_y * vpar0_y) * xjac * tstep &
!==========================End of viscous heating terms==============================
!===================== Additional terms from friction terms============
                    + v * BigR * ((GAMMA - 1.)/2.) * vpar0**2 * BB2 * (source_bg + source_imp) * xjac * tstep &
                    + v * BigR * ((GAMMA - 1.)/2.) * vv2 * (source_bg + source_imp)            * xjac * tstep &
!==============================End of friction terms=================
                    ! Energy exchange term
                    + v * BigR * dTi_e                                                         * xjac * tstep


         rhs_ij_6_k =  - (ZK_i_par_T-ZK_i_prof) * BigR / BB2 * Bgrad_T_k_star * Bgrad_Ti * xjac * tstep  &
                       - ZK_i_prof * BigR * (                + v_p*Ti0_p /BigR**2 )      * xjac * tstep  &
                       - TG_num6 * 0.25d0 / BigR * vpar0**2 &
                               * Ti0 * ((r0_x+alpha_i*rn0_x) * ps0_y - (r0_y+alpha_i*rn0_y) * ps0_x &
                                         + F0 / BigR * (r0_p+alpha_i*rn0_p))                        &
                               * (                                 + F0 / BigR * v_p) * xjac * tstep * tstep   &
   
                       - TG_num6 * 0.25d0 / BigR * vpar0**2 &
                               * (r0+alpha_i*rn0) * (Ti0_x * ps0_y - Ti0_y * ps0_x + F0 / BigR * Ti0_p)        &
                               * (                                 + F0 / BigR * v_p) * xjac * tstep * tstep 


!###################################################################################################
!#  equation 7   (parallel velocity equation)                                                      #
!###################################################################################################

         rhs_ij_7 = - v * F0 / BigR * P0_p                                              * xjac * tstep &
                    - v * (P0_s * ps0_t - P0_t * ps0_s)                                        * tstep &

                      - visco_par * (v_x * vpar0_x + v_y * vpar0_y) * BigR                * xjac * tstep &

                    - 0.5d0 * r0 * vpar0**2 * BB2 * (ps0_s * v_t - ps0_t * v_s)                * tstep &

                    - 0.5d0 * v  * vpar0**2 * BB2 * (ps0_s * r0_t - ps0_t * r0_s)              * tstep &         
                    + 0.5d0 * v  * vpar0**2 * BB2 * F0 / BigR * r0_p                    * xjac * tstep &

                    - visco_par_num * (v_xx + v_x/Bigr + v_yy)*(vpar0_xx + vpar0_x/Bigr + vpar0_yy) * BigR * xjac * tstep &

                    + zeta * v * delta_g(mp,7,ms,mt) * R0 * F0**2 / BigR                        * xjac &
                    + zeta * v * r0 * vpar0 * (ps0_x * delta_ps_x + ps0_y * delta_ps_y) / BigR  * xjac &

                    ! New terms coming from -(\partial_t \rho + \nabla \cdot (\rho \mathbf{v})) \mathbf{v} in RHS of momentum equation
                    ! (see wiki: https://www.jorek.eu/wiki/doku.php?id=model500_501_555#equations):
                    + zeta * v * delta_g(mp,5,ms,mt) * vpar0 * F0**2 / BigR              * xjac  &
                    + v * (r0_x_hat * u0_y - r0_y_hat * u0_x)       * vpar0 * BB2 * xjac * tstep &
                    - v * F0 / BigR * (r0 * vpar0_p + r0_p * vpar0) * vpar0 * BB2 * xjac * tstep &
                    - v * r0 * (vpar0_x * ps0_y - vpar0_y * ps0_x)  * vpar0 * BB2 * xjac * tstep &
                    - v * vpar0 * (r0_x * ps0_y - r0_y * ps0_x)     * vpar0 * BB2 * xjac * tstep &

                    ! Old term (not to be included anymore due to implementations of terms above):
                    ! - v*(particle_source(ms,mt) + source_pellet) * vpar0 * BB2   * BigR * xjac * tstep &

            - TG_NUM7 * 0.25d0 * r0 * Vpar0**2 * BB2 &
                      * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac + F0 / BigR * vpar0_p) / BigR  &
                      * (-(ps0_s * v_t     - ps0_t * v_s)    /xjac                      )  * xjac * tstep * tstep &
!            - TG_NUM7 * 0.25d0 * v  * Vpar0**2 * BB2 &
!                      * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac + F0 / BigR * vpar0_p) / BigR  &
!                      * (-(ps0_s * r0_t    - ps0_t * r0_s)   /xjac + F0 / BigR * r0_p)  * xjac * tstep * tstep &

!=============================== New TG_num terms==================================
            - TG_NUM7 * 0.25d0 * vpar0 * Vpar0**2 * BB2 &
                      * (-(ps0_s * r0_t - ps0_t * r0_s)/xjac + F0 / BigR * r0_p) / BigR  &
                      * (-(ps0_s * v_t     - ps0_t * v_s)    /xjac)  * xjac * tstep * tstep
!===============================End of new TG_num terms============================


!=============================== New TG_num terms==================================
!
!            - TG_NUM7 * 0.25d0 * vpar0 * Vpar0**2 * BB2 &
!                      * (-(ps0_s * r0_t - ps0_t * r0_s)/xjac + F0 / BigR * r0_p) / BigR  &
!                      * (-(ps0_s * v_t  - ps0_t * v_s) /xjac   )  * xjac * tstep * tstep &
!            - TG_NUM7 * 0.25d0 * v  * Vpar0**2 * BB2 &
!                      * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac + F0 / BigR * vpar0_p) / BigR  &
!                      * (-(ps0_s * r0_t    - ps0_t * r0_s)   /xjac + F0 / BigR * r0_p)  * xjac * tstep * tstep
!
!===============================End of new TG_num terms============================

    
             rhs_ij_7_k = + 0.5d0 * r0 * vpar0**2 * BB2 * F0 / BigR * v_p                     * xjac * tstep &

            - TG_NUM7 * 0.25d0 * r0 * Vpar0**2 * BB2 &
                      * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac + F0 / BigR * vpar0_p) / BigR  &
                      * (                                          + F0 / BigR * v_p)  * xjac * tstep * tstep &
 
!=============================== New TG_num terms==================================

            - TG_NUM7 * 0.25d0 * vpar0 * Vpar0**2 * BB2 &
                      * (-(ps0_s * r0_t - ps0_t * r0_s)/xjac + F0 / BigR * r0_p) / BigR  &
                      * (                                          + F0 / BigR * v_p)  * xjac * tstep * tstep

!===============================End of new TG_num terms============================


!################################################################################################### 
!#  equation 8 (impurity density equation)                                                          # 
!################################################################################################### 


	   rhs_ij_8 =   BigR* (- Dn0x * rn0_x * v_x - Dn0y * rn0_y * v_y)                                            * xjac * tstep &       
                    ! The new diffusion scheme for the impurities
                      - (D_par_imp-D_prof_imp) * BigR / BB2 * Bgrad_rho_star * (Bgrad_rhon) * xjac * tstep &
                    ! The new diffusion scheme for the impurities
                      - D_prof_imp * BigR  * (v_x*(rn0_x) + v_y*(rn0_y)              )  * xjac * tstep &


                      + v * BigR**2 * ( rn0_s * u0_t - rn0_t * u0_s)                                                    * tstep &
                      + v * 2.d0 * BigR * rn0 * u0_y                                                             * xjac * tstep &
                      - v * F0 / BigR * Vpar0 * rn0_p                                                            * xjac * tstep &
                      - v * Vpar0 * (rn0_s * ps0_t - rn0_t * ps0_s)                                                     * tstep &
                      - v * F0 / BigR * rn0 * vpar0_p                                                            * xjac * tstep &
                      - v * rn0 * (vpar0_s * ps0_t - vpar0_t * ps0_s)                                                   * tstep &

                      - TG_num8 * 0.25d0 * BigR**3 * (rn0_x * u0_y - rn0_y * u0_x)                                              &
                                                 * ( v_x * u0_y - v_y * u0_x) * xjac * tstep * tstep                            &
                      - TG_num8 * 0.25d0 / BigR * vpar0**2                                                                      &
                              * (rn0_x * ps0_y - rn0_y * ps0_x + F0 / BigR * rn0_p)                                             &
                              * ( v_x * ps0_y -  v_y * ps0_x                   ) * xjac * tstep * tstep                         &

         	      + BigR * v * source_imp                                    * xjac * tstep                                 &

                      + v * delta_g(mp,8,ms,mt) * BigR * xjac * zeta                          &
                      - Dn_perp_num * (v_xx + v_x/Bigr + v_yy)*(rn0_xx + rn0_x/Bigr + rn0_yy) &
                        * BigR * xjac * tstep


           rhs_ij_8_k =  BigR* ( - Dn0p * rn0_p * v_p*eps_cyl**2/BigR**2)                              * xjac * tstep            & 
                       ! The new diffusion scheme for the impurities
                       - (D_par_imp-D_prof_imp) * BigR / BB2 * Bgrad_rho_k_star * (Bgrad_rhon)         * xjac * tstep &
                       - D_prof_imp * BigR  * (          v_p*(rn0_p) * eps_cyl**2  /BigR**2)           * xjac * tstep &

                        - TG_num8 * 0.25d0 / BigR * vpar0**2                                                                    &
                              * (rn0_x * ps0_y - rn0_y * ps0_x + F0 / BigR * rn0_p)                                             &
                              * (                            + F0 / BigR * v_p) * xjac * tstep * tstep 

!###################################################################################################
!#  equation 9 (electron energy  equation)                                                         #
!###################################################################################################

         rhs_ij_9 =   v * BigR * heat_source_e(ms,mt)                                  * xjac * tstep &
 
                    + v * (r0 + rn0*alpha_e_bis) * BigR**2 * ( Te0_s * u0_t - Te0_t * u0_s)   * tstep &
                    + v * Te0 * BigR**2 * ( r0_s * u0_t - r0_t * u0_s)                        * tstep &
                    + v * alpha_e * Te0 * BigR**2 * (rn0_s * u0_t - rn0_t * u0_s)             * tstep &

                    + v * (r0 + rn0*alpha_e) * Te0 * 2.d0* GAMMA * BigR * u0_y         * xjac * tstep &

                    - v * (r0 + rn0*alpha_e_bis) * F0 / BigR * Vpar0 * Te0_p           * xjac * tstep &
                    - v * Te0 * F0 / BigR * Vpar0 * (r0_p + alpha_e * rn0_p)           * xjac * tstep &

                    - v * (r0 + rn0*alpha_e_bis) * Vpar0 * (Te0_s * ps0_t - Te0_t * ps0_s)    * tstep &
                    - v * Te0 * Vpar0 * (r0_s * ps0_t - r0_t * ps0_s)                         * tstep &
                    - v * Te0 * Vpar0 * alpha_e * (rn0_s * ps0_t - rn0_t * ps0_s)             * tstep &

                    - v * (r0 + rn0*alpha_e) * Te0 * GAMMA * (vpar0_s * ps0_t - vpar0_t * ps0_s)        * tstep &
                    - v * (r0 + rn0*alpha_e) * Te0 * GAMMA * F0 / BigR * vpar0_p                 * xjac * tstep &

                    - (ZK_e_par_T-ZK_e_prof) * BigR / BB2 * Bgrad_T_star * Bgrad_Te    * xjac * tstep &
                    - ZK_e_prof * BigR * (v_x*Te0_x + v_y*Te0_y                      ) * xjac * tstep &

                    - ZK_perp_num  *  (v_xx + v_x/Bigr + v_yy)*(Te0_xx + Te0_x/Bigr + Te0_yy) * BigR * xjac * tstep &
                     
                    - ZK_par_num * (v_ps0_x  * ps0_y - v_ps0_y  * ps0_x)                              &
                                 * (Te0_ps0_x * ps0_y - Te0_ps0_y * ps0_x)             * xjac * tstep &

                    - TG_num9 * 0.25d0 * BigR**3 * Te0 * ((r0_x+alpha_e*rn0_x) * u0_y - (r0_y+alpha_e*rn0_y) * u0_x) &
                                       * ( v_x * u0_y - v_y * u0_x) * xjac * tstep * tstep  &
                    - TG_num9 * 0.25d0 * BigR**3 * (r0+alpha_e_bis*rn0) * (Te0_x * u0_y - Te0_y * u0_x) &
                                       * ( v_x * u0_y - v_y * u0_x) * xjac * tstep * tstep  &

                    - TG_num9 * 0.25d0 / BigR * vpar0**2 &
                              * Te0 * ((r0_x+alpha_e*rn0_x) * ps0_y - (r0_y+alpha_e*rn0_y) * ps0_x &
                                        + F0 / BigR * (r0_p+alpha_e*rn0_p))                        &
                              * ( v_x * ps0_y -  v_y * ps0_x                  ) * xjac * tstep * tstep        &

                    - TG_num9 * 0.25d0 / BigR * vpar0**2 &
                              * (r0+alpha_e_bis*rn0) * (Te0_x * ps0_y - Te0_y * ps0_x + F0 / BigR * Te0_p)    &
                              * ( v_x * ps0_y -  v_y * ps0_x                  ) * xjac * tstep * tstep        &

                    + zeta * v * (r0 + rn0 * alpha_e_bis) * delta_g(mp,9,ms,mt) * BigR                 * xjac &
                    + zeta * v * Te0 * delta_g(mp,5,ms,mt) * BigR                                      * xjac &
                    + zeta * v * alpha_e * Te0 * delta_g(mp,8,ms,mt) * BigR                            * xjac &   

!===================== Additional terms from ionization energy terms============
                    + (GAMMA-1.) * zeta * v * E_ion * delta_g(mp,8,ms,mt) *BigR                           * xjac &
                    + (GAMMA-1.) * zeta * v * dE_ion_dT * rn0 * delta_g(mp,9,ms,mt) *BigR                 * xjac &
                    + (GAMMA-1.) * zeta * v * E_ion_bg * (delta_g(mp,5,ms,mt) - delta_g(mp,8,ms,mt))*BigR * xjac &

                    + (GAMMA-1.) * v * rn0 * dE_ion_dT * BigR**2 * ( Te0_s * u0_t - Te0_t * u0_s)   * tstep &
                    + (GAMMA-1.) * v * E_ion * BigR**2 * (rn0_s * u0_t - rn0_t * u0_s)              * tstep &
                    + (GAMMA-1.) * v * E_ion_bg * BigR**2*((r0_s-rn0_s)*u0_t - (r0_t-rn0_t)*u0_s)   * tstep &

                    - (GAMMA-1.) * v * rn0 * dE_ion_dT * F0 / BigR * Vpar0 * Te0_p           * xjac * tstep &
                    - (GAMMA-1.) * v * E_ion * F0 / BigR * Vpar0 * rn0_p                     * xjac * tstep &
                    - (GAMMA-1.) * v * E_ion_bg * F0 / BigR * Vpar0 * (r0_p - rn0_p)         * xjac * tstep &

                    - (GAMMA-1.) * v * rn0 * dE_ion_dT * Vpar0 * (Te0_s * ps0_t - Te0_t * ps0_s)    * tstep &
                    - (GAMMA-1.) * v * E_ion * Vpar0 * (rn0_s * ps0_t - rn0_t * ps0_s)              * tstep &
                    - (GAMMA-1.) * v * E_ion_bg * Vpar0*((r0_s-rn0_s)*ps0_t - (r0_t-rn0_t)*ps0_s)   * tstep &

                    + (GAMMA-1.) * v * E_ion * rn0 * 2.d0 * BigR * u0_y                      * xjac * tstep &
                    - (GAMMA-1.) * v * E_ion * rn0 * (vpar0_s * ps0_t - vpar0_t * ps0_s)            * tstep &
                    - (GAMMA-1.) * v * E_ion * rn0 * F0 / BigR * vpar0_p                     * xjac * tstep &

                    + (GAMMA-1.) * v * E_ion_bg * (r0-rn0) * 2.d0 * BigR * u0_y              * xjac * tstep &
                    - (GAMMA-1.) * v * E_ion_bg * (r0-rn0) * (vpar0_s * ps0_t - vpar0_t * ps0_s)    * tstep &
                    - (GAMMA-1.) * v * E_ion_bg * (r0-rn0) * F0 / BigR * vpar0_p             * xjac * tstep &

                      ! New diffusive flux of the ionization potential energy for impurities
                    - (GAMMA - 1.) * E_ion * (D_par_imp-D_prof_imp) * BigR / BB2 * Bgrad_rho_star * (Bgrad_rhon) * xjac * tstep &
                    - (GAMMA - 1.) * E_ion * D_prof_imp * BigR  * (v_x*(rn0_x) + v_y*(rn0_y)                                &
                                                                                                 )       * xjac * tstep &
                    - (GAMMA - 1.) * E_ion_bg * (D_par-D_prof) * BigR / BB2 * Bgrad_rho_star * (Bgrad_rho-Bgrad_rhon)   &
                                                                                                         * xjac * tstep &
                    - (GAMMA - 1.) * E_ion_bg * D_prof * BigR  * (v_x*(r0_x-rn0_x) + v_y*(r0_y-rn0_y)                   &
                                                                                                      )  * xjac * tstep &

!==============================End of ionization energy terms=================
                    + v * BigR * ((GAMMA-1.)/(BigR**2)) * eta_T_ohm * zj0**2                 * xjac * tstep &
                    - v * BigR * (r0_corr+alpha_e*rn0_corr) * rn0_corr * Lrad                * xjac * tstep &
                    - v * BigR * (r0_corr+alpha_e*rn0_corr) * frad_bg                        * xjac * tstep &
                    ! Energy exchange term
                    + v * BigR * dTe_i                                                       * xjac * tstep

         rhs_ij_9_k =  - (ZK_e_par_T-ZK_e_prof) * BigR / BB2 * Bgrad_T_k_star * Bgrad_Te * xjac * tstep  &
                       - ZK_e_prof * BigR * (                + v_p*Te0_p /BigR**2 )      * xjac * tstep  &
!===================== Additional terms from ionization energy terms============
                    - (GAMMA - 1.) * E_ion * (D_par_imp-D_prof_imp) * BigR / BB2 * Bgrad_rho_k_star * (Bgrad_rhon) * xjac * tstep &
                    - (GAMMA - 1.) * E_ion * D_prof_imp * BigR * (                                                            &
                                                              + v_p*(rn0_p) * eps_cyl**2 /BigR**2 )        * xjac * tstep &
                    - (GAMMA - 1.) * E_ion_bg * (D_par-D_prof) * BigR / BB2 * Bgrad_rho_k_star * (Bgrad_rho-Bgrad_rhon)   &
                                                                                                           * xjac * tstep &
                    - (GAMMA - 1.) * E_ion_bg * D_prof * BigR * (                                                         &
                                                              + v_p*(r0_p-rn0_p) * eps_cyl**2 /BigR**2 )   * xjac * tstep &

!==============================End of ionization energy terms=================
                       - TG_num9 * 0.25d0 / BigR * vpar0**2 &
                               * Te0 * ((r0_x+alpha_e*rn0_x) * ps0_y - (r0_y+alpha_e*rn0_y) * ps0_x &
                                         + F0 / BigR * (r0_p+alpha_e*rn0_p))                        &
                               * (                                 + F0 / BigR * v_p) * xjac * tstep * tstep   &

                       - TG_num9 * 0.25d0 / BigR * vpar0**2 &
                               * (r0+alpha_e_bis*rn0) * (Te0_x * ps0_y - Te0_y * ps0_x + F0 / BigR * Te0_p)    &
                               * (                                 + F0 / BigR * v_p) * xjac * tstep * tstep
 
!###################################################################################################
!#  RHS equations end                                                                              #
!###################################################################################################

         ij1 = index_ij
         ij2 = index_ij + 1
         ij3 = index_ij + 2
         ij4 = index_ij + 3
         ij5 = index_ij + 4
         ij6 = index_ij + 5
         ij7 = index_ij + 6
         ij8 = index_ij + 7
         ij9 = index_ij + 8

         RHS_p(mp,ij1) = RHS_p(mp,ij1) + rhs_ij_1 * wst
         RHS_p(mp,ij2) = RHS_p(mp,ij2) + rhs_ij_2 * wst
         RHS_p(mp,ij3) = RHS_p(mp,ij3) + rhs_ij_3 * wst
         RHS_p(mp,ij4) = RHS_p(mp,ij4) + rhs_ij_4 * wst
         RHS_p(mp,ij5) = RHS_p(mp,ij5) + rhs_ij_5 * wst
         RHS_p(mp,ij6) = RHS_p(mp,ij6) + rhs_ij_6 * wst
         RHS_p(mp,ij7) = RHS_p(mp,ij7) + rhs_ij_7 * wst
         RHS_p(mp,ij8) = RHS_p(mp,ij8) + rhs_ij_8 * wst
         RHS_p(mp,ij9) = RHS_p(mp,ij9) + rhs_ij_9 * wst

         RHS_k(mp,ij5) = RHS_k(mp,ij5) + rhs_ij_5_k * wst
         RHS_k(mp,ij6) = RHS_k(mp,ij6) + rhs_ij_6_k * wst
         RHS_k(mp,ij7) = RHS_k(mp,ij7) + rhs_ij_7_k * wst
         RHS_k(mp,ij8) = RHS_k(mp,ij8) + rhs_ij_8_k * wst
         RHS_k(mp,ij9) = RHS_k(mp,ij9) + rhs_ij_9_k * wst

         do k=1,n_vertex_max

           do l=1,n_degrees

             psi   = H(k,l,ms,mt) * element%size(k,l)

             psi_x = (   y_t(ms,mt) * h_s(k,l,ms,mt) - y_s(ms,mt) * h_t(k,l,ms,mt) ) / xjac * element%size(k,l)
             psi_y = ( - x_t(ms,mt) * h_s(k,l,ms,mt) + x_s(ms,mt) * h_t(k,l,ms,mt) ) / xjac * element%size(k,l)

             psi_p  = H(k,l,ms,mt)   * element%size(k,l)
             psi_s  = h_s(k,l,ms,mt) * element%size(k,l)
             psi_t  = h_t(k,l,ms,mt) * element%size(k,l)
             psi_ss = h_ss(k,l,ms,mt) * element%size(k,l)
             psi_tt = h_tt(k,l,ms,mt) * element%size(k,l)
             psi_st = h_st(k,l,ms,mt) * element%size(k,l)

             psi_xx = (psi_ss * y_t(ms,mt)**2 - 2.d0*psi_st * y_s(ms,mt)*y_t(ms,mt) + psi_tt * y_s(ms,mt)**2  &
                    + psi_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                              &
	            + psi_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) ) )    / xjac**2               & 
		    - xjac_x * (psi_s * y_t(ms,mt) - psi_t * y_s(ms,mt)) / xjac**2

             psi_yy = (psi_ss * x_t(ms,mt)**2 - 2.d0*psi_st * x_s(ms,mt)*x_t(ms,mt) + psi_tt * x_s(ms,mt)**2  &
	            + psi_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                              &
	            + psi_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) ) )       / xjac**2            &
		    - xjac_y * (- psi_s * x_t(ms,mt) + psi_t * x_s(ms,mt) ) / xjac**2
           
	     psi_xy = (- psi_ss * y_t(ms,mt)*x_t(ms,mt) - psi_tt * x_s(ms,mt)*y_s(ms,mt)                      &
     	            + psi_st * (y_s(ms,mt)*x_t(ms,mt)  + y_t(ms,mt)*x_s(ms,mt)  )                             &        
                    - psi_s  * (x_st(ms,mt)*y_t(ms,mt) - x_tt(ms,mt)*y_s(ms,mt) )                             &	   
	            - psi_t  * (x_st(ms,mt)*y_s(ms,mt) - x_ss(ms,mt)*y_t(ms,mt) )  )  / xjac**2               &		
                    - xjac_x * (- psi_s * x_t(ms,mt) + psi_t * x_s(ms,mt) )   / xjac**2

             u    = psi    ;  zj    = psi    ;  w    = psi    ; rho    = psi    ;  Ti    = psi    ;  Te    = psi    ; vpar    = psi    ;    rhon    = psi    ;
             u_x  = psi_x  ;  zj_x  = psi_x  ;  w_x  = psi_x  ; rho_x  = psi_x  ;  Ti_x  = psi_x  ;  Te_x  = psi_x  ; vpar_x  = psi_x  ;    rhon_x  = psi_x  ;
             u_y  = psi_y  ;  zj_y  = psi_y  ;  w_y  = psi_y  ; rho_y  = psi_y  ;  Ti_y  = psi_y  ;  Te_y  = psi_y  ; vpar_y  = psi_y  ;    rhon_y  = psi_y  ;
             u_p  = psi_p  ;  zj_p  = psi_p  ;  w_p  = psi_p  ; rho_p  = psi_p  ;  Ti_p  = psi_p  ;  Te_p  = psi_p  ; vpar_p  = psi_p  ;    rhon_p  = psi_p  ;
             u_s  = psi_s  ;  zj_s  = psi_s  ;  w_s  = psi_s  ; rho_s  = psi_s  ;  Ti_s  = psi_s  ;  Te_s  = psi_s  ; vpar_s  = psi_s  ;    rhon_s  = psi_s  ;
             u_t  = psi_t  ;  zj_t  = psi_t  ;  w_t  = psi_t  ; rho_t  = psi_t  ;  Ti_t  = psi_t  ;  Te_t  = psi_t  ; vpar_t  = psi_t  ;    rhon_t  = psi_t  ;
             u_ss = psi_ss ;  zj_ss = psi_ss ;  w_ss = psi_ss ; rho_ss = psi_ss ;  Ti_ss = psi_ss ;  Te_ss = psi_ss ; vpar_ss = psi_ss ;    rhon_ss = psi_ss ;
             u_tt = psi_tt ;  zj_tt = psi_tt ;  w_tt = psi_tt ; rho_tt = psi_tt ;  Ti_tt = psi_tt ;  Te_tt = psi_tt ; vpar_tt = psi_tt ;    rhon_tt = psi_tt ;
             u_st = psi_st ;  zj_st = psi_st ;  w_st = psi_st ; rho_st = psi_st ;  Ti_st = psi_st ;  Te_st = psi_st ; vpar_st = psi_st ;    rhon_st = psi_st ;

             u_xx = psi_xx ;                                    rho_xx = psi_xx ;  Ti_xx = psi_xx ;  Te_xx = psi_xx ; vpar_xx = psi_xx ;     rhon_xx = psi_xx ;
             u_yy = psi_yy ;                                    rho_yy = psi_yy ;  Ti_yy = psi_yy ;  Te_yy = psi_yy ; vpar_yy = psi_yy ;     rhon_yy = psi_yy ;
             u_xy = psi_xy ;                                    rho_xy = psi_xy ;  Ti_xy = psi_xy ;  Te_xy = psi_xy ; vpar_xy = psi_xy
             
             w_xx = (psi_ss * y_t(ms,mt)**2 - 2.d0*psi_st * y_s(ms,mt)*y_t(ms,mt) + psi_tt * y_s(ms,mt)**2 ) / xjac**2
             w_yy = (psi_ss * x_t(ms,mt)**2 - 2.d0*psi_st * x_s(ms,mt)*x_t(ms,mt) + psi_tt * x_s(ms,mt)**2 ) / xjac**2
             w_xy = psi_xy

             rho_hat   = BigR**2 * rho
             rho_x_hat = 2.d0 * BigR * BigR_x  * rho + BigR**2 * rho_x
             rho_y_hat = BigR**2 * rho_y

             Bgrad_rho_star_psi = ( v_x  * psi_y - v_y  * psi_x ) / BigR
             Bgrad_rho_psi      = ( r0_x * psi_y - r0_y * psi_x ) / BigR
             Bgrad_rhon_psi     = ( rn0_x * psi_y - rn0_y * psi_x ) / BigR
             Bgrad_rho_rho      = ( rho_x * ps0_y - rho_y * ps0_x ) / BigR
             Bgrad_rho_rhon     = ( rhon_x * ps0_y - rhon_y * ps0_x ) / BigR
             Bgrad_rho_rho_n    = ( F0 / BigR * rho_p ) / BigR
             Bgrad_rho_rhon_n    = ( F0 / BigR * rhon_p ) / BigR
             BB2_psi            = 2.d0 * (psi_x * ps0_x + psi_y * ps0_y ) /BigR**2

             Btheta2_psi  = 2.d0 * (psi_x * ps0_x + psi_y * ps0_y ) /BigR**2
	         rhon_hat   = BigR**2 * rhon                                     
             rhon_x_hat = 2.d0 * BigR * BigR_x  * rhon + BigR**2 * rhon_x    
             rhon_y_hat = BigR**2 * rhon_y                                   

             index_kl = n_var*n_degrees*(k-1) + n_var * (l-1) + 1   ! index in the ELM matrix

!###################################################################################################
!#  equation 1   (induction equation)                                                              #
!###################################################################################################

             amat_11 = v * psi / BigR * xjac * (1.d0 + zeta)                                              &
                     - v * (psi_s * u0_t - psi_t * u0_s)                                  * theta * tstep &
                     + v * tauIC/(r0*BB2) * F0**2/BigR**2 * (psi_s * pi0_t - psi_t * pi0_s) * theta * tstep & 
                     - v * tauIC/(r0_corr*BB2**2) * BB2_psi * F0**2/BigR**2 * (ps0_x*p0_y - ps0_y*p0_x) * xjac * theta * tstep &
                     + v * tauIC/(r0_corr*BB2**2) * BB2_psi * F0**3/BigR**3 * p0_p           * xjac * theta * tstep 

             amat_12 = -  v * (ps0_s * u_t - ps0_t * u_s)                             * theta * tstep

             amat_12_n = +  eps_cyl * F0 / BigR * v * u_p * xjac                      * theta * tstep

             amat_13 = - eta_num_T * (v_x * zj_x + v_y * zj_y)                 * xjac * theta * tstep  &
                       - eta_T * v * zj / BigR                                 * xjac * theta * tstep

	     amat_15 = + v * tauIC/(r0_corr*BB2) * F0**2/BigR**2 * Ti0  * (ps0_s * rho_t - ps0_t * rho_s) * theta * tstep &
		       + v * tauIC/(r0_corr*BB2) * F0**2/BigR**2 * rho * (ps0_s * Ti0_t  - ps0_t * Ti0_s)  * theta * tstep &
		       - v * tauIC/(r0_corr*BB2) * F0**3/BigR**3 * eps_cyl * rho * Ti0_p  * xjac          * theta * tstep &

                      - v * tauIC * rho /(r0_corr**2 * BB2) * F0**2/BigR**2 * (ps0_s * pi0_t - ps0_t * pi0_s) * theta * tstep &
                      + v * tauIC * rho /(r0_corr**2 * BB2) * F0**3/BigR**3 * eps_cyl * pi0_p * xjac         * theta * tstep &
                      ! The density gradient term from Z_eff
                      - deta_dr0 * v * rho * (zj0-current_source(ms,mt)) / BigR * xjac * theta * tstep

             amat_15_n = - v * tauIC/(r0_corr*BB2) * F0**3/BigR**3 * eps_cyl * Ti0  * rho_p * xjac * theta * tstep 

             amat_16 = + v * tauIC/(r0_corr*BB2) * F0**2/BigR**2 * r0 * (ps0_s * Ti_t  - ps0_t * Ti_s) * theta * tstep &
                       + v * tauIC/(r0_corr*BB2) * F0**2/BigR**2 * Ti  * (ps0_s * r0_t - ps0_t * r0_s) * theta * tstep &
                       - v * tauIC/(r0_corr*BB2) * F0**3/BigR**3 * eps_cyl * Ti  * r0_p         * xjac * theta * tstep 

             amat_16_n = - v * tauIC/(r0_corr*BB2) * F0**3/BigR**3 * eps_cyl * r0 * Ti_p        * xjac * theta * tstep 

             ! The density gradient term from Z_eff
             amat_18 = - deta_drn0 * v * rhon * (zj0-current_source(ms,mt)) / BigR * xjac * theta * tstep

             amat_19 = - deta_dT * v * Te * (zj0-current_source(ms,mt)-Jb) / BigR  * xjac * theta * tstep &
                       ! Temperature dependent hyper-resistivity
                       - deta_num_dT * Te * (v_x * zj0_x + v_y * zj0_y)            * xjac * theta * tstep

!###################################################################################################
!#  equation 2   (perpendicular momentum equation)                                                 #
!###################################################################################################

             amat_21 = - v * (psi_s * zj0_t - psi_t * zj0_s)                          * theta * tstep &

                       ! New terms coming from -(\partial_t \rho + \nabla \cdot (\rho \mathbf{v})) \mathbf{v} in RHS of momentum equation
                       ! (see wiki: https://www.jorek.eu/wiki/doku.php?id=model500_501_555#equations):
                       - BigR**2 * r0 * (vpar0_x * psi_y - vpar0_y * psi_x) * (v_x * u0_x + v_y * u0_y) * xjac * theta * tstep &
                       - BigR**2 * vpar0 * (r0_x * psi_y - r0_y * psi_x)    * (v_x * u0_x + v_y * u0_y) * xjac * theta * tstep


             amat_22 = - BigR * r0_hat * (v_x * u_x + v_y * u_y) * xjac  * (1.d0 + zeta)           &
                       + r0_hat * BigR**2 * w0 * (v_s * u_t  - v_t  * u_s)                              * theta * tstep &
                       + BigR**2 * (u_x * u0_x + u_y * u0_y) * (v_x * r0_y_hat - v_y * r0_x_hat) * xjac * theta * tstep &

                       + tauIC * BigR**3 * pi0_y * (v_x* u_x + v_y * u_y)                         * xjac * theta * tstep &
                       + v * tauIC * BigR**4 * (u_xy * (pi0_xx - pi0_yy) - pi0_xy * (u_xx - u_yy))  * xjac * theta * tstep &

                       ! New terms coming from -(\partial_t \rho + \nabla \cdot (\rho \mathbf{v})) \mathbf{v} in RHS of momentum equation
                       ! (see wiki: https://www.jorek.eu/wiki/doku.php?id=model500_501_555#equations):
                       + BigR**2 * (r0_x_hat * u0_y - r0_y_hat * u0_x)      * (v_x * u_x  + v_y * u_y)  * xjac * theta * tstep &
                       + BigR**2 * (r0_x_hat * u_y  - r0_y_hat * u_x)       * (v_x * u0_x + v_y * u0_y) * xjac * theta * tstep &
                       - BigR * F0 * (r0 * vpar0_p + vpar0 * r0_p)          * (v_x * u_x  + v_y * u_y)  * xjac * theta * tstep &
                       - BigR**2 * r0 * (vpar0_x * ps0_y - vpar0_y * ps0_x) * (v_x * u_x  + v_y * u_y)  * xjac * theta * tstep &
                       - BigR**2 * vpar0 * (r0_x * ps0_y - r0_y * ps0_x)    * (v_x * u_x + v_y * u_y)   * xjac * theta * tstep &

                       ! Old term (not to be included anymore due to implementations of terms above):
                       !- BigR**3 * (particle_source(ms,mt)+source_pellet) * (v_x * u_x + v_y * u_y) * xjac * theta * tstep &

                       + TG_num2 * 0.25d0 * r0_hat * BigR**3 * (w0_x * u_y - w0_y * u_x)       &
                                 * ( v_x * u0_y - v_y * u0_x) * xjac * theta * tstep * tstep   &
                       + TG_num2 * 0.25d0 * r0_hat * BigR**3 * (w0_x * u0_y - w0_y * u0_x)     &
                                 * ( v_x * u_y - v_y * u_x)   * xjac * theta * tstep * tstep   &

!====================================New TG_num terms=================================
                      + TG_num2 * 0.25d0 * w0 * BigR**3 * (r0_x_hat * u_y - r0_y_hat * u_x) &
                                * ( v_x * u0_y - v_y * u0_x) * theta * xjac * tstep * tstep &

                      + TG_num2 * 0.25d0 * w0 * BigR**3 * (r0_x_hat * u0_y - r0_y_hat * u0_x) &
                                * ( v_x * u_y - v_y * u_x) * theta * xjac * tstep * tstep 

!===============================End of NewTG_num terms==============================


             amat_23 = - v * (ps0_s * zj_t  - ps0_t * zj_s)                * theta * tstep

             amat_23_n = + eps_cyl * F0 / BigR * v * zj_p  * xjac          * theta * tstep

             amat_24 = r0_hat * BigR**2 * w  * ( v_s * u0_t - v_t * u0_s)  * theta * tstep  &
                     + BigR * ( v_x * w_x + v_y * w_y) * visco_T  * xjac   * theta * tstep  &
                     + v * tauIC * BigR**4 * (pi0_s * w_t - pi0_t * w_s)   * theta * tstep  & 

                     + visco_num_T * (v_xx + v_x/BigR + v_yy)*(w_xx + w_x/BigR + w_yy) * xjac * theta * tstep    &

                     + TG_num2 * 0.25d0 * r0_hat * BigR**3 * (w_x * u0_y - w_y * u0_x)     &
                               * ( v_x * u0_y - v_y * u0_x) * xjac * theta * tstep * tstep &

!====================================New TG_num terms=================================
                      + TG_num2 * 0.25d0 * w * BigR**3 * (r0_x_hat * u0_y - r0_y_hat * u0_x) &
                                * ( v_x * u0_y - v_y * u0_x) * theta * xjac * tstep * tstep 

!===============================End of NewTG_num terms==============================


             amat_25 = + 0.5d0 * vv2 * (v_x * rho_y_hat - v_y * rho_x_hat)   * xjac * theta * tstep &
                       + rho_hat * BigR**2 * w0 * (v_s * u0_t - v_t * u0_s)         * theta * tstep &
                       - BigR**2 * (v_s*rho_t*(Ti0+Te0) - v_t*rho_s*(Ti0+Te0))      * theta * tstep &
                       - BigR**2 * (v_s*rho*(Ti0_t+Te0_t) - v_t*rho*(Ti0_s+Te0_s))  * theta * tstep &
 
                       ! New terms coming from -(\partial_t \rho + \nabla \cdot (\rho \mathbf{v})) \mathbf{v} in RHS of momentum equation
                       ! (see wiki: https://www.jorek.eu/wiki/doku.php?id=model500_501_555#equations):	     	    	     
	               - BigR**3 * rho * (v_x * u0_x + v_y * u0_y) * xjac  * (1.d0 + zeta)          &
                       + BigR**2 * (rho_x_hat * u0_y - rho_y_hat * u0_x)     * (v_x * u0_x + v_y * u0_y) * xjac * theta * tstep &
                       - BigR * rho * F0 * vpar0_p                           * (v_x * u0_x + v_y * u0_y) * xjac * theta * tstep &
                       - BigR**2 * rho * (vpar0_x * ps0_y - vpar0_y * ps0_x) * (v_x * u0_x + v_y * u0_y) * xjac * theta * tstep &
                       - BigR**2 * vpar0 * (rho_x * ps0_y - rho_y * ps0_x)   * (v_x * u0_x + v_y * u0_y) * xjac * theta * tstep &

                       + v * tauIC * BigR**4 * Ti0  * (rho_s * w0_t - rho_t * w0_s)  * theta * tstep &
                       + v * tauIC * BigR**4 * rho * (Ti0_s  * w0_t - Ti0_t  * w0_s) * theta * tstep &
                       + tauIC * BigR**3 * (Ti0_y * rho + Ti0 * rho_y) * (v_x* u0_x + v_y * u0_y) * xjac * theta * tstep   &
                       + v * tauIC * BigR**4 * ( (u0_xy * (rho_xx*Ti0 + 2.d0*rho_x*Ti0_x + rho*Ti0_xx                      &
                             -  rho_yy*Ti0 - 2.d0*rho_y*Ti0_y - rho*Ti0_yy))                        &                 
                             - (rho_xy * Ti0 + rho_x*Ti0_y + rho_y*Ti0_x + rho*Ti0_xy) * (u0_xx - u0_yy)  )   &
                                                                                                  * xjac * theta * tstep   &
                       + TG_num2 * 0.25d0 * rho_hat * BigR**3 * (w0_x * u0_y - w0_y * u0_x)         &
                                 * ( v_x * u0_y - v_y * u0_x) * xjac * theta * tstep * tstep        &

!====================================New TG_num terms=================================
                      + TG_num2 * 0.25d0 * w0 * BigR**3 * (rho_x_hat * u0_y - rho_y_hat * u0_x) &
                                * ( v_x * u0_y - v_y * u0_x) * theta * xjac * tstep * tstep

!===============================End of NewTG_num terms==============================

             ! New term coming from -(\partial_t \rho + \nabla \cdot (\rho \mathbf{v})) \mathbf{v} in RHS of momentum equation
             ! (see wiki: https://www.jorek.eu/wiki/doku.php?id=model500_501_555#equations):	 
             amat_25_n = - BigR * vpar0 * F0 * rho_p * (v_x * u0_x + v_y * u0_y) * xjac * theta * tstep

             amat_26 = - BigR**2 * (v_s * r0_t * Ti   - v_t * r0_s * Ti)      * theta * tstep  &
                       - BigR**2 * (v_s * r0   * Ti_t - v_t * r0   * Ti_s)    * theta * tstep  &
                       - BigR**2 * (v_s * rn0_t * alpha_i * Ti - v_t * rn0_s * alpha_i * Ti) * theta * tstep  &
                       - BigR**2 * (v_s * rn0 * alpha_i * Ti_t - v_t * rn0 * alpha_i * Ti_s) * theta * tstep  &

                       + v * tauIC * BigR**4 * r0 * (Ti_s  * w0_t - Ti_t  * w0_s)  * theta * tstep &
                       + v * tauIC * BigR**4 * Ti   * (r0_s * w0_t - r0_t * w0_s)  * theta * tstep &

                       + tauIC * BigR**3 * (r0_y * Ti + r0 * Ti_y) * (v_x * u0_x + v_y * u0_y) * xjac * theta * tstep &
                       + v * tauIC * BigR**4 * ( (u0_xy * (Ti_xx * r0 + 2.d0 * Ti_x * r0_x + Ti * r0_xx         &
                                                           - Ti_yy*r0 - 2.d0 * Ti_y * r0_y - Ti * r0_yy))       &
                                                - (Ti_xy * r0 + Ti_x*r0_y + Ti_y*r0_x + Ti*r0_xy) * (u0_xx - u0_yy)  ) &
                                                                                             * xjac * theta * tstep

             amat_27 = &
                       ! New terms coming from -(\partial_t \rho + \nabla \cdot (\rho \mathbf{v})) \mathbf{v} in RHS of momentum equation
                       ! (see wiki: https://www.jorek.eu/wiki/doku.php?id=model500_501_555#equations):	     	    	     	  
                       - BigR * vpar * F0 * r0_p                          * (v_x * u0_x + v_y * u0_y) * xjac * theta * tstep &
                       - BigR**2 * r0 * (vpar_x * ps0_y - vpar_y * ps0_x) * (v_x * u0_x + v_y * u0_y) * xjac * theta * tstep &
                       - BigR**2 * vpar * (r0_x * ps0_y - r0_y * ps0_x)   * (v_x * u0_x + v_y * u0_y) * xjac * theta * tstep

             ! New term coming from -(\partial_t \rho + \nabla \cdot (\rho \mathbf{v})) \mathbf{v} in RHS of momentum equation
             ! (see wiki: https://www.jorek.eu/wiki/doku.php?id=model500_501_555#equations):	 
             amat_27_n = - BigR * r0 * F0 * vpar_p * (v_x * u0_x + v_y * u0_y) * xjac * theta * tstep

             amat_28 = - BigR**2 * (v_s * rhon_t * alpha_i * Ti0 - v_t * rhon_s * alpha_i * Ti0)         * theta * tstep &
                       - BigR**2 * (v_s * rhon * alpha_i * Ti0_t - v_t * rhon * alpha_i * Ti0_s)         * theta * tstep &
                       - BigR**2 * (v_s * rhon_t * alpha_e * Te0     - v_t * rhon_s * alpha_e * Te0)     * theta * tstep &
                       - BigR**2 * (v_s * rhon * alpha_e_bis * Te0_t - v_t * rhon * alpha_e_bis * Te0_s) * theta * tstep

             amat_29 = - BigR**2 * (v_s * r0_t * Te   - v_t * r0_s * Te)      * theta * tstep  &
                       - BigR**2 * (v_s * r0   * Te_t - v_t * r0   * Te_s)    * theta * tstep  &
                       - BigR**2 * (v_s * rn0_t * alpha_e_bis * Te - v_t * rn0_s * alpha_e_bis * Te) * theta * tstep  &
                       - BigR**2 * (v_s * rn0 * alpha_e_bis * Te_t - v_t * rn0 * alpha_e_bis * Te_s) * theta * tstep  &
                       - BigR**2 * (v_s * rn0 * alpha_e_tri * Te * Te0_t &
                                    - v_t * rn0 * alpha_e_tri * Te * Te0_s) * theta * tstep  &
                       + dvisco_dT * Te * ( v_x * w0_x + v_y * w0_y ) * BigR * xjac  * theta * tstep &
                       ! Hyper-viscosity terms
                       + dvisco_num_dT * Te * (v_xx + v_x/Bigr + v_yy)*(w0_xx + w0_x/Bigr + w0_yy)*xjac * theta * tstep 

!###################################################################################################
!#  equation 3   (current definition)                                                              #
!###################################################################################################

             amat_33 = v * zj / BigR * xjac                               
             amat_31 = (v_x * psi_x + v_y * psi_y ) / BigR * xjac          

!###################################################################################################
!#  equation 4   (vorticity definition)                                                            #
!###################################################################################################

             amat_44 =  v * w * BigR * xjac                               
             amat_42 = (v_x * u_x + v_y * u_y) * BigR * xjac               

!###################################################################################################
!#  equation 5   (density equation)                                                                #
!###################################################################################################

             amat_51 = &
                       - (D_par-D_prof) * BigR * BB2_psi/ BB2**2 * Bgrad_rho_star     * (Bgrad_rho-Bgrad_rhon) * xjac * theta * tstep &
                       + (D_par-D_prof) * BigR / BB2             * Bgrad_rho_star_psi * (Bgrad_rho-Bgrad_rhon) * xjac * theta * tstep &
                       + (D_par-D_prof) * BigR / BB2             * Bgrad_rho_star     * (Bgrad_rho_psi-Bgrad_rhon_psi) * xjac * theta * tstep &
                       - (D_par_imp-D_prof_imp) * BigR * BB2_psi/ BB2**2 * Bgrad_rho_star     * (Bgrad_rhon)   * xjac * theta * tstep &
                       + (D_par_imp-D_prof_imp) * BigR / BB2             * Bgrad_rho_star_psi * (Bgrad_rhon)   * xjac * theta * tstep &
                       + (D_par_imp-D_prof_imp) * BigR / BB2             * Bgrad_rho_star     * (Bgrad_rhon_psi)       * xjac * theta * tstep &  
                       !Old diffusion scheme for impurities
                       !- (D_par-D_prof) * BigR * BB2_psi/ BB2**2 *  Bgrad_rho_star     * (Bgrad_rho) * xjac * theta * tstep &
                       !+ (D_par-D_prof) * BigR / BB2             * Bgrad_rho_star_psi * (Bgrad_rho) * xjac * theta * tstep &
                       !+ (D_par-D_prof) * BigR / BB2             * Bgrad_rho_star     * (Bgrad_rho_psi) * xjac * theta * tstep &
                       + v * Vpar0 * (r0_s * psi_t - r0_t * psi_s)                                           * theta * tstep &
                       + v * r0 * (vpar0_s * psi_t - vpar0_t * psi_s)                                        * theta * tstep &
 
                       + TG_num5 * 0.25d0 / BigR * vpar0**2                                                              &
                                 * (r0_x * psi_y - r0_y * psi_x)                                                         &
                                 * ( v_x * ps0_y -  v_y * ps0_x                   ) * xjac * theta * tstep * tstep       &
                       + TG_num5 * 0.25d0 / BigR * vpar0**2                                                              &
                                 * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                                      &
                                 * ( v_x * psi_y -  v_y * psi_x                   ) * xjac * theta * tstep * tstep


             amat_51_k = &
                         - (D_par-D_prof) * BigR * BB2_psi/ BB2**2 * Bgrad_rho_k_star * (Bgrad_rho-Bgrad_rhon)         * xjac * theta * tstep &
                         + (D_par-D_prof) * BigR / BB2             * Bgrad_rho_k_star * (Bgrad_rho_psi-Bgrad_rhon_psi) * xjac * theta * tstep &
                         - (D_par_imp-D_prof_imp) * BigR * BB2_psi/ BB2**2 * Bgrad_rho_k_star * (Bgrad_rhon)           * xjac * theta * tstep &
                         + (D_par_imp-D_prof_imp) * BigR / BB2             * Bgrad_rho_k_star * (Bgrad_rhon_psi)       * xjac * theta * tstep &
                         !- (D_par-D_prof) * BigR * BB2_psi/ BB2**2 * Bgrad_rho_k_star * (Bgrad_rho)         * xjac * theta * tstep &
                         !+ (D_par-D_prof) * BigR / BB2             * Bgrad_rho_k_star * (Bgrad_rho_psi) * xjac * theta * tstep &

                         + TG_num5 * 0.25d0 / BigR * vpar0**2                                                            &
                                 * (r0_x * psi_y - r0_y * psi_x)                                                         &
                                 * (                            + F0 / BigR * v_p) * xjac * theta * tstep * tstep


             amat_51_n =  0.d0

             amat_52 =  - v * BigR**2 * ( r0_s * u_t - r0_t * u_s)                                        * theta * tstep &
                        - v * 2.d0 * BigR * r0 * u_y                                               * xjac * theta * tstep &

                        + TG_num5 * 0.25d0 * BigR**3 * (r0_x * u_y  - r0_y * u_x)                                      &
                                                     * ( v_x * u0_y - v_y  * u0_x) * xjac * theta * tstep * tstep      &
                        + TG_num5 * 0.25d0 * BigR**3 * (r0_x * u0_y - r0_y * u0_x)                                     &
                                                     * ( v_x * u_y  - v_y  * u_x)  * xjac * theta * tstep * tstep 

             amat_55 = v * rho * BigR * (1.d0 + zeta)                                              * xjac   &
                     - v * BigR**2 * ( rho_s * u0_t - rho_t * u0_s)                                       * theta * tstep &
                     - v * 2.d0 * BigR * rho * u0_y                                                * xjac * theta * tstep &
                     + (D_par-D_prof) * BigR / BB2 * Bgrad_rho_star * Bgrad_rho_rho                * xjac * theta * tstep &
                     + D_prof * BigR  * (v_x*rho_x + v_y*rho_y )                                   * xjac * theta * tstep &
                     + v * Vpar0 * (rho_s * ps0_t - rho_t * ps0_s)                                        * theta * tstep &
                     + v * rho * (vpar0_s * ps0_t - vpar0_t * ps0_s)                                      * theta * tstep &
                     + v * rho * F0 / BigR * vpar0_p                                               * xjac * theta * tstep &

                     - v * 2.d0 * tauIC * (rho_y * Ti0 + rho*Ti0_y) * BigR                         * xjac * theta * tstep &

                     + D_perp_num * (v_xx + v_x/BigR + v_yy)*(rho_xx + rho_x/BigR + rho_yy)   * BigR * xjac * theta * tstep &

                     + TG_num5 * 0.25d0 * BigR**3 * (rho_x * u0_y - rho_y * u0_x)                                &
                                                  * ( v_x  * u0_y - v_y   * u0_x) * xjac * theta * tstep * tstep &

                     + TG_num5 * 0.25d0 / BigR * vpar0**2                                                        &
                               * (rho_x * ps0_y - rho_y * ps0_x )                             &
                               * ( v_x * ps0_y -  v_y * ps0_x   ) * xjac * theta * tstep * tstep 


             amat_55_k = + (D_par-D_prof) * BigR / BB2 * Bgrad_rho_k_star * Bgrad_rho_rho          * xjac * theta * tstep &
                        
                    + TG_num5 * 0.25d0 / BigR * vpar0**2                                                        &
                               * (rho_x * ps0_y - rho_y * ps0_x                  )                              &
                               * (                              + F0 / BigR * v_p) * xjac * theta * tstep * tstep

             amat_55_n = + (D_par-D_prof) * BigR / BB2 * Bgrad_rho_star   * Bgrad_rho_rho_n        * xjac * theta * tstep &
                         + v * F0 / BigR * Vpar0 * rho_p                                           * xjac * theta * tstep &

                     + TG_num5 * 0.25d0 / BigR * vpar0**2                                                        &
                               * (                              + F0 / BigR * rho_p)                             &
                               * ( v_x * ps0_y -  v_y * ps0_x                      ) * xjac * theta * tstep * tstep

             amat_55_kn = + (D_par-D_prof) * BigR / BB2 * Bgrad_rho_k_star * Bgrad_rho_rho_n       * xjac * theta * tstep &
                          + D_prof * BigR  * ( v_p*rho_p * eps_cyl**2 /BigR**2 )                   * xjac * theta * tstep &

                     + TG_num5 * 0.25d0 / BigR * vpar0**2                                                        &
                               * ( + F0 / BigR * rho_p)                                                          &
                               * ( + F0 / BigR * v_p) * xjac * theta * tstep * tstep

             amat_56   = - v * 2.d0 * tauIC * (Ti_y * r0 + Ti*r0_y) * BigR                         * xjac * theta * tstep 

             amat_57   = + v * F0 / BigR * Vpar * r0_p                                             * xjac * theta * tstep &
                         + v * Vpar * (r0_s * ps0_t - r0_t * ps0_s)                                       * theta * tstep &
                         + v * r0 * (vpar_s * ps0_t - vpar_t * ps0_s)                                     * theta * tstep &

                         + TG_num5 * 0.25d0 / BigR * 2.d0*vpar0*vpar                                   &
                              * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                       &
                              * ( v_x * ps0_y -  v_y * ps0_x                   ) * xjac * theta * tstep * tstep 

             amat_57_k = + TG_num5 * 0.25d0 / BigR * 2.d0*vpar0*vpar                                   &
                              * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                       &
                              * (                            + F0 / BigR * v_p) * xjac * theta * tstep * tstep 

             amat_57_n = + v * r0 * F0 / BigR * vpar_p                                             * xjac * theta * tstep
	     
             amat_58   = &
                         + BigR * (Dn0x * rhon_x * v_x + Dn0y * rhon_y * v_y)                      * xjac * theta * tstep &
                         - (D_par-D_prof) * BigR / BB2 * Bgrad_rho_star * Bgrad_rho_rhon           * xjac * theta * tstep &
                         + (D_par_imp-D_prof_imp) * BigR / BB2 * Bgrad_rho_star * Bgrad_rho_rhon   * xjac * theta * tstep &
                         - D_prof * BigR  * (v_x*rhon_x + v_y*rhon_y )                             * xjac * theta * tstep &
                         + D_prof_imp * BigR  * (v_x*rhon_x + v_y*rhon_y )                         * xjac * theta * tstep

             amat_58_k  = - (D_par-D_prof) * BigR / BB2 * Bgrad_rho_k_star * Bgrad_rho_rhon       * xjac * theta * tstep &
                          + (D_par_imp-D_prof_imp) * BigR / BB2 * Bgrad_rho_k_star * Bgrad_rho_rhon * xjac * theta * tstep

             amat_58_n  = - (D_par-D_prof) * BigR / BB2 * Bgrad_rho_star   * Bgrad_rho_rhon_n     * xjac * theta * tstep &
                          + (D_par_imp-D_prof_imp) * BigR / BB2 * Bgrad_rho_star * Bgrad_rho_rhon_n * xjac * theta * tstep

             amat_58_kn = &
                          + Dn0p * rhon_p * v_p*eps_cyl**2/BigR * xjac * theta * tstep &
                          - (D_par-D_prof) * BigR / BB2 * Bgrad_rho_k_star * Bgrad_rho_rhon_n       * xjac * theta * tstep &
                          + (D_par_imp-D_prof_imp) * BigR / BB2 * Bgrad_rho_k_star * Bgrad_rho_rhon_n * xjac * theta * tstep &
                          - D_prof * BigR  * ( v_p*rhon_p * eps_cyl**2 /BigR**2 )                   * xjac * theta * tstep &
                          + D_prof_imp * BigR  * ( v_p*rhon_p * eps_cyl**2 /BigR**2 )               * xjac * theta * tstep

!###################################################################################################
!#  equation 6   (ion energy equation)                                                             #
!###################################################################################################
             Bgrad_T_star_psi = ( v_x  * psi_y - v_y  * psi_x  ) / BigR
             Bgrad_Ti_psi     = ( Ti0_x * psi_y - Ti0_y * psi_x )  / BigR
             Bgrad_Ti_T       = ( Ti_x * ps0_y - Ti_y * ps0_x ) / BigR 
             Bgrad_Ti_T_n     = ( F0 / BigR * Ti_p) / BigR

             Ti_ps0_x = Ti_xx * ps0_y - Ti_xy * ps0_x + Ti_x * ps0_xy - Ti_y * ps0_xx
             Ti_ps0_y = Ti_xy * ps0_y - Ti_yy * ps0_x + Ti_x * ps0_yy - Ti_y * ps0_xy

             Ti0_psi_x = Ti0_xx * psi_y - Ti0_xy * psi_x + Ti0_x * psi_xy - Ti0_y * psi_xx
             Ti0_psi_y = Ti0_xy * psi_y - Ti0_yy * psi_x + Ti0_x * psi_yy - Ti0_y * psi_xy

             v_psi_x = v_xx * psi_y - v_xy * psi_x + v_x * psi_xy - v_y * psi_xx
             v_psi_y = v_xy * psi_y - v_yy * psi_x + v_x * psi_yy - v_y * psi_xy
             
             amat_61 = - (ZK_i_par_T-ZK_i_prof)*BigR * BB2_psi / BB2**2 * Bgrad_T_star * Bgrad_Ti * xjac * theta * tstep &
                       + (ZK_i_par_T-ZK_i_prof) * BigR / BB2     * Bgrad_T_star_psi    * Bgrad_Ti * xjac * theta * tstep &
                       + (ZK_i_par_T-ZK_i_prof) * BigR / BB2     * Bgrad_T_star    * Bgrad_Ti_psi * xjac * theta * tstep &
                       + v * (r0 + rn0 * alpha_i) * Vpar0 * (Ti0_s * psi_t - Ti0_t * psi_s)              * theta * tstep &
                       + v * Ti0 * Vpar0 * ((r0_s+rn0_s*alpha_i)*psi_t - (r0_t+rn0_t*alpha_i)*psi_s)     * theta * tstep &
                       + v * (r0 + rn0 * alpha_i) * GAMMA * Ti0 * (vpar0_s * psi_t - vpar0_t * psi_s)    * theta * tstep &

!===================== Additional terms from friction terms============
                       - v * ((GAMMA - 1.) / BigR) * vpar0**2 * (psi_x * ps0_x + psi_y * ps0_y)                          &
                           * (source_bg + source_imp)                                             * xjac * theta * tstep &
!==============================End of friction terms=================

                       + ZK_par_num * (v_psi_x  * ps0_y - v_psi_y  * ps0_x + v_ps0_x * psi_y - v_ps0_y * psi_x)          &
                                    * (Ti0_ps0_x * ps0_y - Ti0_ps0_y * ps0_x)                     * xjac * theta * tstep &
                       + ZK_par_num * (Ti0_psi_x * ps0_y - Ti0_psi_y * ps0_x + Ti0_ps0_x * psi_y - Ti0_ps0_y * psi_x)    &
                                    * (v_ps0_x  * ps0_y - v_ps0_y  * ps0_x)                       * xjac * theta * tstep &

                       + TG_num6 * 0.25d0 / BigR * vpar0**2                                                       &
                                 * Ti0 * ((r0_x+alpha_i*rn0_x) * psi_y - (r0_y+alpha_i*rn0_y) * psi_x)            &
                                 * ( v_x * ps0_y -  v_y * ps0_x                  ) * xjac * theta * tstep * tstep &

                       + TG_num6 * 0.25d0 / BigR * vpar0**2                                                       &
                                 * (r0+alpha_i*rn0) * (Ti0_x * psi_y - Ti0_y * psi_x)                             &
                                 * ( v_x * ps0_y -  v_y * ps0_x                  ) * xjac * theta * tstep * tstep &

                       + TG_num6 * 0.25d0 / BigR * vpar0**2                                                       &
                                 * Ti0 * ((r0_x+alpha_i*rn0_x) * ps0_y - (r0_y+alpha_i*rn0_y) * ps0_x             &
                                         + F0 / BigR * (r0_p+alpha_i*rn0_p))                                      &
                                 * ( v_x * psi_y -  v_y * psi_x ) * xjac * theta * tstep * tstep                  &

                       + TG_num6 * 0.25d0 / BigR * vpar0**2                                                       &
                                 * (r0+alpha_i*rn0) * (Ti0_x * ps0_y - Ti0_y * ps0_x + F0 / BigR * Ti0_p)         &
                                 * ( v_x * psi_y -  v_y * psi_x ) * xjac * theta * tstep * tstep

             amat_61_k = - (ZK_i_par_T-ZK_i_prof)*BigR*BB2_psi / BB2**2 * Bgrad_T_k_star * Bgrad_Ti * xjac * theta * tstep &
                         + (ZK_i_par_T-ZK_i_prof) * BigR / BB2      * Bgrad_T_k_star * Bgrad_Ti_psi * xjac * theta * tstep &

                         + TG_num6 * 0.25d0 / BigR * vpar0**2                                                       &
                                   * Ti0 * ((r0_x+alpha_i*rn0_x) * psi_y - (r0_y+alpha_i*rn0_y) * psi_x)            &
                                   * (                            + F0 / BigR * v_p) * xjac * theta * tstep * tstep &

                         + TG_num6 * 0.25d0 / BigR * vpar0**2                                                       &
                                   * (r0+alpha_i*rn0) * (Ti0_x * psi_y - Ti0_y * psi_x)                             &
                                   * (                            + F0 / BigR * v_p) * xjac * theta * tstep * tstep 


             amat_62 = - v * (r0 + rn0 * alpha_i) * BigR**2 * ( Ti0_s * u_t - Ti0_t * u_s)             * theta * tstep &
                       - v * Ti0 * BigR**2 * ((r0_s+rn0_s*alpha_i)*u_t - (r0_t+rn0_t*alpha_i)*u_s)     * theta * tstep &
                       - v * (r0 + rn0 * alpha_i) * 2.d0* GAMMA * BigR * Ti0 * u_y              * xjac * theta * tstep &

!===================== Additional terms from friction terms============
                       - v * BigR**3 * (GAMMA - 1.) * (u_x * u0_x + u_y * u0_y)                                        &
                           * (source_bg + source_imp)                                           * xjac * theta * tstep &
!==============================End of friction terms=================

                       + TG_num6 * 0.25d0 * BigR**2 * Ti0* ((r0_x+alpha_i*rn0_x) * u_y - (r0_y+alpha_i*rn0_y) * u_x) &
                                 * ( v_x * u0_y - v_y * u0_x) * xjac * theta*tstep*tstep  &

                       + TG_num6 * 0.25d0 * BigR**2 * (r0+alpha_i*rn0)* (Ti0_x * u_y - Ti0_y * u_x)                  &
                                 * ( v_x * u0_y - v_y * u0_x) * xjac * theta*tstep*tstep  &

                       + TG_num6 * 0.25d0 * BigR**2 * Ti0* ((r0_x+alpha_i*rn0_x)*u0_y - (r0_y+alpha_i*rn0_y)*u0_x)   &
                                 * ( v_x * u_y - v_y * u_x) * xjac * theta*tstep*tstep    &

                       + TG_num6 * 0.25d0 * BigR**2 * (r0+alpha_i*rn0)* (Ti0_x * u0_y - Ti0_y * u0_x)                &
                                 * ( v_x * u_y - v_y * u_x) * xjac * theta*tstep*tstep


             amat_63 = 0.


             amat_65 =  + v * rho * Ti0 * BigR * xjac * (1.d0 + zeta)    &

                        - v * rho * BigR**2 * ( Ti0_s  * u0_t - Ti0_t  * u0_s)                      * theta * tstep &
                        - v * Ti0  * BigR**2 * ( rho_s * u0_t - rho_t * u0_s)                       * theta * tstep &

                        - v * rho * 2.d0* GAMMA * BigR * Ti0 * u0_y                          * xjac * theta * tstep &

                        + v * rho * F0 / BigR * Vpar0 * Ti0_p                                * xjac * theta * tstep &

                        + v * rho * Vpar0 * (Ti0_s  * ps0_t - Ti0_t  * ps0_s)                       * theta * tstep &
                        + v * Ti0  * Vpar0 * (rho_s * ps0_t - rho_t * ps0_s)                        * theta * tstep &

                        + v * rho * GAMMA * Ti0 * (vpar0_s * ps0_t - vpar0_t * ps0_s)               * theta * tstep &
                        + v * rho * GAMMA * Ti0 * F0 / BigR * vpar0_p                        * xjac * theta * tstep &

                        ! Energy exchange term
                        - v * BigR * ddTi_e_drho * rho                                       * xjac * theta * tstep &

                        + TG_num6 * 0.25d0 * BigR**2 * Ti0* (rho_x * u0_y - rho_y * u0_x)     &
                                  * ( v_x * u0_y - v_y * u0_x) * xjac * theta*tstep*tstep     &

                        + TG_num6 * 0.25d0 * BigR**2 * rho * (Ti0_x * u0_y - Ti0_y * u0_x)    &
                                  * ( v_x * u0_y - v_y * u0_x) * xjac* theta*tstep*tstep      &

                        + TG_num6 * 0.25d0 / BigR * vpar0**2 &
                                  * Ti0 * (rho_x * ps0_y - rho_y * ps0_x                    )                     &
                                  * ( v_x * ps0_y -  v_y * ps0_x                  ) * xjac * theta * tstep * tstep&

                        + TG_num6 * 0.25d0 / BigR * vpar0**2 &
                                  * rho * (Ti0_x * ps0_y - Ti0_y * ps0_x + F0 / BigR * Ti0_p)                      &
                                  * ( v_x * ps0_y -  v_y * ps0_x                  ) * xjac * theta * tstep * tstep

             amat_65_n = + v * Ti0  * F0 / BigR * Vpar0 * rho_p                             * xjac * theta * tstep &

                    + TG_num6 * 0.25d0 / BigR * vpar0**2 &
                              * Ti0 * (                             + F0 / BigR * rho_p)                      &
                              * ( v_x * ps0_y -  v_y * ps0_x                  ) * xjac * theta * tstep * tstep

             amat_65_k = + TG_num6 * 0.25d0 / BigR * vpar0**2 &
                              * Ti0 * (rho_x * ps0_y - rho_y * ps0_x                   )                      &
                              * (                            + F0 / BigR * v_p) * xjac * theta * tstep * tstep&
                         + TG_num6 * 0.25d0 / BigR * vpar0**2 &
                              * rho * (Ti0_x * ps0_y - Ti0_y * ps0_x + F0 / BigR * Ti0_p)                     &
                              * (                            + F0 / BigR * v_p) * xjac * theta * tstep * tstep

             amat_65_kn = + TG_num6 * 0.25d0 / BigR * vpar0**2 &
                              * Ti0 * (+ F0 / BigR * rho_p)                      &
                              * (      + F0 / BigR * v_p) * xjac * theta * tstep * tstep

             amat_66 =   v * (r0 + rn0 * alpha_i) * Ti * BigR * xjac * (1.d0 + zeta)                           &
                       - v * (r0 + rn0 * alpha_i) * BigR**2 * ( Ti_s  * u0_t - Ti_t  * u0_s)   * theta * tstep &
                       - v * Ti  * BigR**2 * ( r0_s * u0_t - r0_t * u0_s)                      * theta * tstep &
                       - v * alpha_i * Ti * BigR**2 * (rn0_s * u0_t - rn0_t * u0_s)            * theta * tstep &

                       - v * (r0 + rn0 * alpha_i) * 2.d0* GAMMA * BigR * Ti * u0_y      * xjac * theta * tstep &

                       + v * Ti * F0  / BigR * Vpar0 * (r0_p + rn0_p * alpha_i)         * xjac * theta * tstep &

                       + v * (r0 + rn0 * alpha_i) * Vpar0 * (Ti_s  * ps0_t - Ti_t  * ps0_s)    * theta * tstep &
                       + v * Ti  * Vpar0 * (r0_s * ps0_t - r0_t * ps0_s)                       * theta * tstep &
                       + v * alpha_i   * Ti * Vpar0 * (rn0_s * ps0_t - rn0_t * ps0_s)          * theta * tstep &

                       + v * (r0 + rn0 * alpha_i) * GAMMA * Ti * (vpar0_s * ps0_t - vpar0_t * ps0_s) * theta * tstep &
                       + v * (r0 + rn0 * alpha_i) * GAMMA * Ti * F0 / BigR * vpar0_p          * xjac * theta * tstep &

                       ! Energy exchange term
                       - v * BigR * ddTi_e_dTi * Ti                                           * xjac * theta * tstep &

                       + (ZK_i_par_T-ZK_i_prof) * BigR / BB2 * Bgrad_T_star * Bgrad_Ti_T      * xjac * theta * tstep &
                       + ZK_i_prof * BigR * (v_x*Ti_x + v_y*Ti_y                     )        * xjac * theta * tstep &

                       + dZK_i_par_dT * Ti * BigR / BB2 * Bgrad_T_star * Bgrad_Ti             * xjac * theta * tstep &

                       + ZK_perp_num * (v_xx + v_x/BigR + v_yy)*(Ti_xx+Ti_x/BigR+Ti_yy)*BigR  * xjac * theta * tstep &

                       + TG_num6 * 0.25d0 * BigR**2 * Ti * ((r0_x+alpha_i*rn0_x)*u0_y - (r0_y+alpha_i*rn0_y)*u0_x)   &
                                 * ( v_x * u0_y - v_y * u0_x) * xjac * theta * tstep * tstep &

                       + TG_num6 * 0.25d0 * BigR**2 * (r0+alpha_i*rn0)* (Ti_x * u0_y - Ti_y * u0_x)                  &
                                 * ( v_x * u0_y - v_y * u0_x) * xjac * theta * tstep * tstep &

                       + TG_num6 * 0.25d0 / BigR * vpar0**2 &
                                 * Ti * ((r0_x+alpha_i*rn0_x)*ps0_y - (r0_y+alpha_i*rn0_y)*ps0_x                 &
                                         + F0/BigR*(r0_p+alpha_i*rn0_p))                                          &
                                 * ( v_x * ps0_y -  v_y * ps0_x                  ) * xjac * theta * tstep * tstep&

                       + TG_num6 * 0.25d0 / BigR * vpar0**2 &
                                 * (r0+alpha_i*rn0) * (Ti_x * ps0_y - Ti_y * ps0_x                  )            &
                                 * ( v_x * ps0_y -  v_y * ps0_x                  ) * xjac * theta * tstep * tstep

             amat_66_k = + (ZK_i_par_T-ZK_i_prof) * BigR / BB2 * Bgrad_T_k_star * Bgrad_Ti_T * xjac * theta * tstep &
                         + dZK_i_par_dT * Ti * BigR / BB2 * Bgrad_T_k_star * Bgrad_Ti        * xjac * theta * tstep &
                         + TG_num6 * 0.25d0 / BigR * vpar0**2 &
                                   * Ti * ((r0_x+alpha_i*rn0_x)*ps0_y - (r0_y+alpha_i*rn0_y)*ps0_x                 &
                                           + F0/BigR*(r0_p+alpha_i*rn0_p))                                          &
                                   * (                          + F0 / BigR * v_p) * xjac * theta * tstep * tstep&

                         + TG_num6 * 0.25d0 / BigR * vpar0**2 &
                                   * (r0+alpha_i*rn0) * (Ti_x * ps0_y - Ti_y * ps0_x                )            &
                                   * (                          + F0 / BigR * v_p) * xjac * theta * tstep * tstep


             amat_66_n = + (ZK_i_par_T-ZK_i_prof) * BigR / BB2 * Bgrad_T_star * Bgrad_Ti_T_n * xjac * theta * tstep &

                         + v * (r0 + rn0 * alpha_i) * F0 / BigR * Vpar0 * Ti_p            * xjac * theta * tstep &

                         + TG_num6 * 0.25d0 / BigR * vpar0**2 &
                                   * (r0+alpha_i*rn0) * (                         + F0 / BigR * Ti_p)            &
                                   * ( v_x * ps0_y -  v_y * ps0_x                ) * xjac * theta * tstep * tstep

             amat_66_kn = + (ZK_i_par_T-ZK_i_prof) * BigR / BB2 * Bgrad_T_k_star * Bgrad_Ti_T_n * xjac * theta * tstep &
                          + ZK_i_prof * BigR * (                    + v_p*Ti_p /BigR**2 )       * xjac * theta * tstep &

                          + TG_num6 * 0.25d0 / BigR * vpar0**2 &
                                    * (r0+alpha_i*rn0) * (                        + F0 / BigR * Ti_p)            &
                                    * (                         + F0 / BigR * v_p) * xjac * theta * tstep * tstep

             amat_67 = + v * (r0 + rn0 * alpha_i) * F0 / BigR * Vpar * Ti0_p            * xjac * theta * tstep &
                       + v * Ti0 * F0 / BigR * Vpar * (r0_p + rn0_p * alpha_i)          * xjac * theta * tstep &

                       + v * (r0 + rn0 * alpha_i) * Vpar * (Ti0_s * ps0_t - Ti0_t * ps0_s)            * theta * tstep &
                       + v * Ti0 * Vpar * ((r0_s+rn0_s*alpha_i)*ps0_t - (r0_t+rn0_t*alpha_i)*ps0_s)   * theta * tstep &

                       + v * (r0 + rn0 * alpha_i) * GAMMA * Ti0 * (vpar_s * ps0_t - vpar_t * ps0_s)   * theta * tstep &

!===================== Additional terms from friction terms============
                       - v * BigR *(GAMMA - 1.) * vpar0 * Vpar * BB2 * (source_bg + source_imp) * xjac * theta * tstep &
!==============================End of friction terms=================
!============================Behold, the parallel viscous heating terms!=============
                       - (GAMMA - 1.) * v * BigR * visco_par_heating * 2.d0 * (vpar_x*vpar0_x + vpar_y*vpar0_y) * xjac * theta * tstep  &
                       - (GAMMA - 1.) * vpar0 * BigR * visco_par_heating    * (vpar_x*v_x     + vpar_y*v_y)     * xjac * theta * tstep  &
                       - (GAMMA - 1.) * vpar * BigR * visco_par_heating    * (vpar0_x*v_x     + vpar0_y*v_y)    * xjac * theta * tstep  &
!==========================End of viscous heating terms==============================

                       + TG_num6 * 0.25d0 / BigR * 2.d0 * vpar0*vpar &
                          * Ti0 * ((r0_x+alpha_i*rn0_x)*ps0_y - (r0_y+alpha_i*rn0_y)*ps0_x                 &
                                  + F0 / BigR * (r0_p+alpha_i*rn0_p))                                      &
                          * ( v_x * ps0_y -  v_y * ps0_x                  ) * xjac * theta * tstep * tstep &

                       + TG_num6 * 0.25d0 / BigR * 2.d0 * vpar0*vpar &
                          * (r0+alpha_i*rn0) * (Ti0_x * ps0_y - Ti0_y * ps0_x + F0 / BigR * Ti0_p)         &
                          * ( v_x * ps0_y -  v_y * ps0_x                  ) * xjac * theta * tstep * tstep


             amat_67_k =  &
                       + TG_num6 * 0.25d0 / BigR * 2.d0 * vpar0*vpar &
                          * Ti0 * ((r0_x+alpha_i*rn0_x)*ps0_y - (r0_y+alpha_i*rn0_y)*ps0_x                 &
                                  + F0 / BigR * (r0_p+alpha_i*rn0_p))                                      &
                          * (                            + F0 / BigR * v_p) * xjac * theta * tstep * tstep &

                       + TG_num6 * 0.25d0 / BigR * 2.d0 * vpar0*vpar &
                          * (r0+alpha_i*rn0) * (Ti0_x * ps0_y - Ti0_y * ps0_x + F0 / BigR * Ti0_p)         &
                          * (                              F0 / BigR * v_p) * xjac * theta * tstep * tstep

             amat_67_n = + v * (r0 + rn0 * alpha_i) * GAMMA * Ti0 * F0 / BigR * vpar_p * xjac * theta * tstep 


             amat_68 =   v * rhon * alpha_i * Ti0 * BigR * xjac * (1.d0 + zeta)                                  &
!=========================New TG_num terms====================================
                       + TG_num6 * 0.25d0 * BigR**2 * Ti0 * alpha_i * (rhon_x * u0_y - rhon_y * u0_x)        &
                                 * ( v_x * u0_y - v_y * u0_x) * xjac * theta*tstep*tstep     &

                       + TG_num6 * 0.25d0 * BigR**2 * alpha_i * rhon * (Ti0_x * u0_y - Ti0_y * u0_x)         &
                                 * ( v_x * u0_y - v_y * u0_x) * xjac* theta*tstep*tstep      &

                       + TG_num6 * 0.25d0 / BigR * vpar0**2 &
                          * Ti0 * alpha_i * (rhon_x * ps0_y - rhon_y * ps0_x                     )           &
                          * ( v_x * ps0_y -  v_y * ps0_x                  ) * xjac * theta * tstep * tstep   &

                       + TG_num6 * 0.25d0 / BigR * vpar0**2 &
                          * alpha_i * rhon * (Ti0_x * ps0_y - Ti0_y * ps0_x + F0 / BigR * Ti0_p)           &
                          * ( v_x * ps0_y -  v_y * ps0_x                  ) * xjac * theta * tstep * tstep &
!===========================End of new TG_num terms===========================

                       - v * rhon * BigR**2 * alpha_i * (Ti0_s * u0_t - Ti0_t * u0_s)         * theta * tstep &
                       - v * alpha_i * Ti0 * BigR**2 * (rhon_s * u0_t - rhon_t * u0_s)        * theta * tstep &
                       + v * rhon * F0 / BigR * Vpar0 * alpha_i * Ti0_p                * xjac * theta * tstep &
                       + v * rhon * Vpar0 * alpha_i * (Ti0_s * ps0_t - Ti0_t * ps0_s)         * theta * tstep &
                       + v * alpha_i * Ti0 * Vpar0 * (rhon_s * ps0_t - rhon_t * ps0_s)        * theta * tstep &

                       - v * alpha_i * rhon * 2.d0* GAMMA * BigR * Ti0 * u0_y                  * xjac * theta * tstep &
                       + v * alpha_i * rhon * GAMMA * Ti0 * (vpar0_s * ps0_t - vpar0_t * ps0_s)       * theta * tstep &
                       + v * alpha_i * rhon * GAMMA * Ti0 * F0 / BigR * vpar0_p                * xjac * theta * tstep &
                       ! Energy exchange term
                       - v * BigR * ddTi_e_drhon * rhon                                        * xjac * theta * tstep

             amat_68_k = &
                       + TG_num6 * 0.25d0 / BigR * vpar0**2 &
                          * Ti0 * alpha_i * (rhon_x * ps0_y - rhon_y * ps0_x                     )           &
                          * (                            + F0 / BigR * v_p) * xjac * theta * tstep * tstep   &

                       + TG_num6 * 0.25d0 / BigR * vpar0**2 &
                          * alpha_i * rhon * (Ti0_x * ps0_y - Ti0_y * ps0_x + F0 / BigR * Ti0_p)           &
                          * (                            + F0 / BigR * v_p) * xjac * theta * tstep * tstep


             amat_68_n = + v * alpha_i * Ti0 * F0 / BigR * Vpar0 * rhon_p              * xjac * theta * tstep &
                       + TG_num6 * 0.25d0 / BigR * vpar0**2 &
                          * Ti0 * alpha_i * (                                + F0 / BigR * rhon_p)           &
                          * ( v_x * ps0_y -  v_y * ps0_x                  ) * xjac * theta * tstep * tstep

             amat_68_kn = &
                       + TG_num6 * 0.25d0 / BigR * vpar0**2 &
                          * Ti0 * alpha_i * (                                + F0 / BigR * rhon_p)           &
                          * (                            + F0 / BigR * v_p) * xjac * theta * tstep * tstep

             amat_69 = - v * BigR * ddTi_e_dTe * Te                                            * xjac * theta * tstep                      

!###################################################################################################
!#  equation 7   (parallel velocity equation)                                                      #
!###################################################################################################

             amat_71 =   v * r0 * vpar0 / BigR * (ps0_x * psi_x + ps0_y * psi_y) * xjac * (1.d0 + zeta) &

                       + v * (P0_s * psi_t - P0_t * psi_s)                                       * theta * tstep &
                       + 0.5d0 * r0 * vpar0**2 * BB2     * (psi_x * v_y  - psi_y * v_x)   * xjac * theta * tstep &
                       + 0.5d0 * r0 * vpar0**2 * BB2_psi * (ps0_x * v_y  - ps0_y * v_x)   * xjac * theta * tstep &
                       + 0.5d0 * v  * vpar0**2 * BB2     * (psi_x * r0_y - psi_y * r0_x)  * xjac * theta * tstep &
                       + 0.5d0 * v  * vpar0**2 * BB2_psi * (ps0_x * r0_y - ps0_y * r0_x)  * xjac * theta * tstep &
                       - 0.5d0 * v  * vpar0**2 * BB2_psi * F0 / BigR * r0_p               * xjac * theta * tstep &

                       ! New terms coming from -(\partial_t \rho + \nabla \cdot (\rho \mathbf{v})) \mathbf{v} in RHS of momentum equation
                       ! (see wiki: https://www.jorek.eu/wiki/doku.php?id=model500_501_555#equations):
                       + v * r0 * (vpar0_x * psi_y - vpar0_y * psi_x) * vpar0 * BB2 * xjac * theta * tstep  &
                       + v * vpar0 * (r0_x * psi_y - r0_y * psi_x)    * vpar0 * BB2 * xjac * theta * tstep  &
                       - v * (r0_x_hat * u0_y - r0_y_hat * u0_x)      * vpar0 * BB2_psi * xjac * theta * tstep &
                       + v * F0 / BigR * (r0 * vpar0_p + r0_p * vpar0)* vpar0 * BB2_psi * xjac * theta * tstep &
                       + v * r0 * (vpar0_x * ps0_y - vpar0_y * ps0_x) * vpar0 * BB2_psi * xjac * theta * tstep &
                       + v * vpar0 * (r0_x * ps0_y - r0_y * ps0_x)    * vpar0 * BB2_psi * xjac * theta * tstep &

                      ! Old term (not to be included anymore due to implementations of terms above):
                      !+ v * (particle_source(ms,mt) + source_pellet) * vpar0 * BB2_psi * BigR * xjac * theta * tstep &
 
                       + TG_NUM7 * 0.25d0 * r0 * Vpar0**2 * BB2 &
                                 * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac) / BigR  &
                                 * (-(psi_s * v_t     - psi_t * v_s)    /xjac)  * xjac * theta * tstep*tstep &

                       + TG_NUM7 * 0.25d0 * r0 * Vpar0**2 * BB2 &
                                 * (-(psi_s * vpar0_t - psi_t * vpar0_s)/xjac) / BigR  &
                                 * (-(ps0_s * v_t     - ps0_t * v_s)    /xjac)  * xjac * theta * tstep*tstep &

!                       + TG_NUM7 * 0.25d0 * v  * Vpar0**2 * BB2 &
!                                 * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac) / BigR  &
!                                 * (-(psi_s * r0_t    - psi_t * r0_s)   /xjac)  * xjac * theta * tstep*tstep &

!                       + TG_NUM7 * 0.25d0 * v  * Vpar0**2 * BB2 &
!                                 * (-(psi_s * vpar0_t - psi_t * vpar0_s)/xjac) / BigR  &
!                                 * (-(ps0_s * r0_t    - ps0_t * r0_s)   /xjac)  * xjac * theta * tstep*tstep &


!=============================== New TG_num terms==================================
                       + TG_NUM7 * 0.25d0 * vpar0 * Vpar0**2 * BB2 &
                                 * (-(ps0_s * r0_t - ps0_t * r0_s)/xjac) / BigR  &
                                 * (-(psi_s * v_t     - psi_t * v_s)    /xjac) * xjac * theta * tstep*tstep &

                       + TG_NUM7 * 0.25d0 * vpar0 * Vpar0**2 * BB2 &
                                 * (-(psi_s * r0_t - psi_t * r0_s)/xjac) / BigR  &
                                 * (-(ps0_s * v_t     - ps0_t * v_s)    /xjac) * xjac * theta * tstep*tstep
!===============================End of new TG_num terms============================

             ! New term coming from -(\partial_t \rho + \nabla \cdot (\rho \mathbf{v})) \mathbf{v} in RHS of momentum equation
             ! (see wiki: https://www.jorek.eu/wiki/doku.php?id=model500_501_555#equations):
             amat_72 = - v * (r0_x_hat * u_y - r0_y_hat * u_x) * vpar0 * BB2 * theta * xjac * tstep
	     

             amat_75 = + v * (rho_s * (Ti0+Te0) * ps0_t - rho_t * (Ti0+Te0) * ps0_s)           * theta * tstep &
                       + v * (rho * (Ti0_s+Te0_s) * ps0_t - rho * (Ti0_t+Te0_t) * ps0_s)       * theta * tstep &
                       + v * F0 / BigR * (                  + rho * (Ti0_p+Te0_p))      * xjac * theta * tstep &

                       + 0.5d0 * rho * vpar0**2 * BB2 * (ps0_s * v_t - ps0_t * v_s)    * theta * tstep &
                       + 0.5d0 * v   * vpar0**2 * BB2 * (ps0_s * rho_t - ps0_t * rho_s)* theta * tstep & 

                       ! New terms coming from -(\partial_t \rho + \nabla \cdot (\rho \mathbf{v})) \mathbf{v} in RHS of momentum equation
                       ! (see wiki: https://www.jorek.eu/wiki/doku.php?id=model500_501_555#equations):
                       + v * rho * vpar0 * F0**2 / BigR * xjac * (1.d0 + zeta)  &
                       - v * (rho_x_hat * u0_y - rho_y_hat * u0_x)       * vpar0 * BB2 * theta * xjac * tstep &   
                       + v * F0 / BigR * rho * vpar0_p                   * vpar0 * BB2 * theta * xjac * tstep &	                            		
                       + v * rho * (vpar0_x * ps0_y - vpar0_y * ps0_x)   * vpar0 * BB2 * theta * xjac * tstep &
                       + v * vpar0 * (rho_x * ps0_y - rho_y * ps0_x)     * vpar0 * BB2 * theta * xjac * tstep &

                       + TG_NUM7 * 0.25d0 * rho * Vpar0**2 * BB2 &
                                 * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac + F0 / BigR * vpar0_p) / BigR  &
                                 * (-(ps0_s * v_t     - ps0_t * v_s)    /xjac          )  * xjac * theta * tstep*tstep &

!                       + TG_NUM7 * 0.25d0 * v * Vpar0**2 * BB2 &
!                                 * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac + F0 / BigR * vpar0_p) / BigR  &
!                                 * (-(ps0_s * rho_t   - ps0_t * rho_s)  /xjac           ) * xjac * theta * tstep*tstep 

!=============================== New TG_num terms==================================

                       + TG_NUM7 * 0.25d0 * vpar0 * Vpar0**2 * BB2 &
                                 * (-(ps0_s * rho_t - ps0_t * rho_s)/xjac                          ) / BigR  &
                                 * (-(ps0_s * v_t     - ps0_t * v_s)    /xjac)  * xjac * theta * tstep*tstep

!===============================End of new TG_num terms============================


             amat_75_k = - 0.5d0 * rho * vpar0**2 * BB2 * F0 / BigR * v_p       * xjac * theta * tstep &

                         + TG_NUM7 * 0.25d0 * rho * Vpar0**2 * BB2 &
                                   * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac + F0 / BigR * vpar0_p) / BigR  &
                                   * (                                          + F0 / BigR * v_p)  * xjac * theta * tstep*tstep &

!=============================== New TG_num terms==================================

                       + TG_NUM7 * 0.25d0 * vpar0 * Vpar0**2 * BB2 &
                                 * (-(ps0_s * rho_t - ps0_t * rho_s)/xjac) / BigR  &
                                 * (+ F0 / BigR * v_p) * xjac * theta * tstep*tstep


!===============================End of new TG_num terms============================


             amat_75_n = + v * F0 / BigR * (rho_p * (Ti0+Te0)              )      * xjac * theta * tstep &
                         - 0.5d0 * v   * vpar0**2 * BB2 * F0 / BigR * rho_p       * xjac * theta * tstep &

                         ! New term coming from -(\partial_t \rho + \nabla \cdot (\rho \mathbf{v})) \mathbf{v} in RHS of momentum equation
                         ! (see wiki: https://www.jorek.eu/wiki/doku.php?id=model500_501_555#equations):	 
                         + v * vpar0 * F0 / BigR * rho_p * vpar0 * BB2          * theta * xjac * tstep & 

!                         + TG_NUM7 * 0.25d0 * v * Vpar0**2 * BB2 &
!                                   * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac + F0 / BigR * vpar0_p) / BigR  &
!                                   * (                                          + F0 / BigR * rho_p)* xjac * theta * tstep*tstep 

!=============================== New TG_num terms==================================

                       + TG_NUM7 * 0.25d0 * vpar0 * Vpar0**2 * BB2 &
                                 * (+ F0 / BigR * rho_p) / BigR  &
                                 * (-(ps0_s * v_t     - ps0_t * v_s)    /xjac) * xjac * theta * tstep*tstep

!===============================End of new TG_num terms============================


!=============================== New TG_num terms==================================

             amat_75_kn = + TG_NUM7 * 0.25d0 * vpar0 * Vpar0**2 * BB2 &
                                 * (+ F0 / BigR * rho_p) / BigR  &
                                 * (+ F0 / BigR * v_p) * xjac * theta * tstep*tstep

!===============================End of new TG_num terms============================

             amat_76 = + v * (Ti_s * r0 * ps0_t - Ti_t * r0 * ps0_s)                          * theta * tstep &
                       + v * (Ti * r0_s * ps0_t - Ti * r0_t * ps0_s)                          * theta * tstep &
                       + v * F0 / BigR * (          + Ti * r0_p)                       * xjac * theta * tstep &
                       + v * (Ti_s * rn0 * alpha_i * ps0_t - Ti_t * rn0 * alpha_i * ps0_s)    * theta * tstep &
                       + v * (Ti * rn0_s * alpha_i * ps0_t - Ti * rn0_t * alpha_i * ps0_s)    * theta * tstep &
                       + v * F0 / BigR * (                     + Ti * rn0_p * alpha_i) * xjac * theta * tstep


             amat_76_n = + v * F0 / BigR * (Ti_p * r0           )                     * xjac * theta * tstep &
                         + v * F0 / BigR * (Ti_p * rn0 * alpha_i                    ) * xjac * theta * tstep


             amat_77 = v * Vpar * r0 * F0**2 / BigR * xjac * (1.d0 + zeta) &
                     + visco_par * (v_x * Vpar_x + v_y * Vpar_y) * BigR           * xjac * theta * tstep &

                       ! New terms coming from -(\partial_t \rho + \nabla \cdot (\rho \mathbf{v})) \mathbf{v} in RHS of momentum equation
                       ! (see wiki: https://www.jorek.eu/wiki/doku.php?id=model500_501_555#equations):
                       - v * (r0_x_hat * u0_y - r0_y_hat * u0_x) * vpar         * BB2 * xjac * theta * tstep &
                       + 2.0 * v * vpar0 * F0 / BigR * r0_p      * vpar         * BB2 * xjac * theta * tstep &
                       + v * r0 * F0 /BigR * vpar0_p             * vpar         * BB2 * xjac * theta * tstep &		      
                       + v * r0 * (vpar_x * ps0_y - vpar_y * ps0_x) * vpar0     * BB2 * xjac * theta * tstep &		      
                       + v * r0 * (vpar0_x * ps0_y - vpar0_y * ps0_x) * vpar    * BB2 * xjac * theta * tstep &		      
                       + 2.0 * v * vpar0 * (r0_x * ps0_y - r0_y * ps0_x) * vpar * BB2 * xjac * theta * tstep  &

                       ! Old term (not to be included anymore due to implementations of terms above):
                       ! + v * (particle_source(ms,mt) + source_pellet)*vpar*BB2 * BigR  * xjac * theta * tstep &

                     + r0 * vpar0 * vpar * BB2 * (ps0_s * v_t - ps0_t * v_s)             * theta * tstep &
                     + v  * vpar0 * vpar * BB2 * (ps0_s * r0_t - ps0_t * r0_s)           * theta * tstep &
                     - v  * vpar0 * vpar * BB2 * F0 / BigR * r0_p                 * xjac * theta * tstep &

                     + visco_par_num * (v_xx + v_x/BigR + v_yy)*(vpar_xx + vpar_x/BigR + vpar_yy) * BigR * xjac * theta * tstep&

            + TG_NUM7 * 0.5d0 * r0 * Vpar * Vpar0 * BB2 &
                      * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac + F0 / BigR * vpar0_p) / BigR                      &
                      * (-(ps0_s * v_t     - ps0_t * v_s)    /xjac                      )  * xjac * theta * tstep*tstep  &
!            + TG_NUM7 * 0.5d0 * v * Vpar * Vpar0 * BB2 &
!                      * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac + F0 / BigR * vpar0_p) / BigR                      &
!                      * (-(ps0_s * r0_t    - ps0_t * r0_s)   /xjac + F0 / BigR * r0_p)  * xjac * theta * tstep*tstep  &
            + TG_NUM7 * 0.25d0 * r0 * Vpar0**2 * BB2 &
                      * (-(ps0_s * vpar_t - ps0_t * vpar_s)/xjac                   ) / BigR                           &
                      * (-(ps0_s * v_t    - ps0_t * v_s)   /xjac                   )  * xjac * theta * tstep*tstep    &
!            + TG_NUM7 * 0.25d0 * v * Vpar0**2 * BB2 &
!                      * (-(ps0_s * vpar_t - ps0_t * vpar_s)/xjac                   ) / BigR                           &
!                      * (-(ps0_s * r0_t   - ps0_t * r0_s)  /xjac + F0 / BigR * r0_p)  * xjac * theta * tstep*tstep    
!=============================== New TG_num terms==================================

            + TG_NUM7 * 0.75d0 * Vpar * Vpar0**2 * BB2 &
                      * (-(ps0_s * r0_t - ps0_t * r0_s)/xjac + F0 / BigR * r0_p) / BigR             &
                      * (-(ps0_s * v_t     - ps0_t * v_s)    /xjac )  * xjac * theta * tstep*tstep  

!===============================End of new TG_num terms============================
 

            amat_77_k = - r0 * vpar0 * vpar * BB2 * F0 / BigR * v_p                 * xjac * theta * tstep               &

            + TG_NUM7 * 0.5d0 * r0 * Vpar * Vpar0 * BB2 &
                      * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac + F0 / BigR * vpar0_p) / BigR                     &
                      * (                                          + F0 / BigR * v_p)  * xjac * theta * tstep*tstep  &
            + TG_NUM7 * 0.25d0 * r0 * Vpar0**2 * BB2 &
                      * (-(ps0_s * vpar_t - ps0_t * vpar_s)/xjac                  ) / BigR                           &
                      * (                                        + F0 / BigR * v_p)  * xjac * theta * tstep*tstep&

!=============================== New TG_num terms==================================

            + TG_NUM7 * 0.75d0 * Vpar * Vpar0**2 * BB2 &
                      * (-(ps0_s * r0_t - ps0_t * r0_s)/xjac + F0 / BigR * r0_p) / BigR             &
                      * (+ F0 / BigR * v_p)  * xjac * theta * tstep*tstep  

!===============================End of new TG_num terms============================

            amat_77_n = &
	    
                        ! New term coming from -(\partial_t \rho + \nabla \cdot (\rho \mathbf{v})) \mathbf{v} in RHS of momentum equation
                        ! (see wiki: https://www.jorek.eu/wiki/doku.php?id=model500_501_555#equations):	    
	                + v * r0 * F0 /BigR * vpar_p * vpar0 * BB2 * xjac * theta * tstep &


                        + TG_NUM7 * 0.25d0 * r0 * Vpar0**2 * BB2 &
                          * F0 / BigR * vpar_p / BigR * (-(ps0_s * v_t - ps0_t * v_s) /xjac) * xjac * theta * tstep * tstep
			  
!            + TG_NUM7 * 0.25d0 * v * Vpar0**2 * BB2 &
!                      * (                                        + F0 / BigR * vpar_p) / BigR                        &
!                      * (-(ps0_s * r0_t   - ps0_t * r0_s)  /xjac + F0 / BigR * r0_p)  * xjac * theta * tstep*tstep 

            amat_77_kn = &

            + TG_NUM7 * 0.25d0 * r0 * Vpar0**2 * BB2 &
                      * (                                        + F0 / BigR * vpar_p) / BigR                        &
                      * (                                        + F0 / BigR * v_p)  * xjac * theta * tstep*tstep

            amat_78 = + v * (rhon_s * alpha_i * Ti0 * ps0_t - rhon_t * alpha_i * Ti0 * ps0_s)         * theta * tstep &
                      + v * (rhon * alpha_i * Ti0_s * ps0_t - rhon * alpha_i * Ti0_t * ps0_s)         * theta * tstep &
                      + v * F0 / BigR * (                       + rhon * alpha_i * Ti0_p)      * xjac * theta * tstep &
                      + v * (rhon_s * alpha_e * Te0 * ps0_t - rhon_t * alpha_e * Te0 * ps0_s)         * theta * tstep &
                      + v * (rhon * alpha_e_bis * Te0_s * ps0_t - rhon * alpha_e_bis * Te0_t * ps0_s) * theta * tstep &
                      + v * F0 / BigR * (                       + rhon * alpha_e_bis * Te0_p)  * xjac * theta * tstep

            amat_78_n = &
                      + v * F0 / BigR * (rhon_p * alpha_i * Ti0                         )      * xjac * theta * tstep &
                      + v * F0 / BigR * (rhon_p * alpha_e * Te0                             )  * xjac * theta * tstep

            amat_79 = + v * (Te_s * r0 * ps0_t - Te_t * r0 * ps0_s)                               * theta * tstep &
                      + v * (Te * r0_s * ps0_t - Te * r0_t * ps0_s)                               * theta * tstep &
                      + v * F0 / BigR * (          + Te * r0_p)                            * xjac * theta * tstep &
                      + v * (Te_s * rn0 * alpha_e_bis * ps0_t  - Te_t * rn0 * alpha_e_bis *  ps0_s) * theta * tstep &
                      + v * (Te0_s* rn0 * alpha_e_tri*Te*ps0_t - Te0_t* rn0 * alpha_e_tri*Te*ps0_s) * theta * tstep &
                      + v * (Te * rn0_s * alpha_e_bis * ps0_t - Te * rn0_t * alpha_e_bis * ps0_s) * theta * tstep &
                      + v * F0 / BigR * (                     + Te * rn0_p * alpha_e_bis)  * xjac * theta * tstep &
                      + v * F0 / BigR * Te0_p * rn0 * Te                   * alpha_e_tri   * xjac * theta * tstep

            amat_79_n = &
                      + v * F0 / BigR * (Te_p * r0            )                            * xjac * theta * tstep &
                      + v * F0 / BigR * (Te_p * rn0 * alpha_e_bis                       )  * xjac * theta * tstep    

!################################################################################################### 
!#  equation 8   neutral density equation                                                          # 
!################################################################################################### 

         amat_81 = &
                       !New diffusion scheme for impurities
                   - (D_par_imp-D_prof_imp) * BigR * BB2_psi/ BB2**2 * Bgrad_rho_star * (Bgrad_rhon)   * xjac * theta * tstep &
                   + (D_par_imp-D_prof_imp) * BigR / BB2             * Bgrad_rho_star_psi * (Bgrad_rhon) * xjac * theta * tstep &
                   + (D_par_imp-D_prof_imp) * BigR / BB2             * Bgrad_rho_star     * (Bgrad_rhon_psi) * xjac * theta * tstep &


                   + v * Vpar0 * (rn0_s * psi_t - rn0_t * psi_s)                                      * theta * tstep &
                   + v * rn0 * (vpar0_s * psi_t - vpar0_t * psi_s)                                    * theta * tstep &

                   + TG_num8 * 0.25d0 / BigR * vpar0**2                                                           &
                             * (rn0_x * psi_y - rn0_y * psi_x)                                                    &
                             * ( v_x * ps0_y -  v_y * ps0_x                   ) * xjac * theta * tstep * tstep    &
                   + TG_num8 * 0.25d0 / BigR * vpar0**2                                                           &
                             * (rn0_x * ps0_y - rn0_y * ps0_x + F0 / BigR * rn0_p)                                &
                             * ( v_x * psi_y -  v_y * psi_x                   ) * xjac * theta * tstep * tstep

         amat_81_k =  &
                   - (D_par_imp-D_prof_imp) * BigR * BB2_psi/ BB2**2 * Bgrad_rho_k_star * (Bgrad_rhon) * xjac * theta * tstep &
                   + (D_par_imp-D_prof_imp) * BigR / BB2             * Bgrad_rho_k_star * (Bgrad_rhon_psi) * xjac * theta * tstep &

                   + TG_num8 * 0.25d0 / BigR * vpar0**2                                                               &
                                 * (rn0_x * psi_y - rn0_y * psi_x)                                                    &
                                 * (                            + F0 / BigR * v_p) * xjac * theta * tstep * tstep

         amat_82 = - v * BigR**2 * ( rn0_s * u_t - rn0_t * u_s)                                       * theta * tstep &
                   - v * 2.d0 * BigR * rn0 * u_y                                               * xjac * theta * tstep &
                   + TG_num8 * 0.25d0 * BigR**3 * (rn0_x * u_y  - rn0_y * u_x)                                        &
                                                     * ( v_x * u0_y - v_y  * u0_x) * xjac * theta * tstep * tstep     &
                   + TG_num8 * 0.25d0 * BigR**3 * (rn0_x * u0_y - rn0_y * u0_x)                                       &
                                                     * ( v_x * u_y  - v_y  * u_x)  * xjac * theta * tstep * tstep 

         amat_85 = 0 ! Place holder
         amat_86 = 0 ! Place holder

      
         amat_87 = + v * F0 / BigR * Vpar * rn0_p                                             *  xjac * theta * tstep &
                   + v * Vpar * (rn0_s * ps0_t - rn0_t * ps0_s)                                       * theta * tstep &
                   + v * rn0 * (vpar_s * ps0_t - vpar_t * ps0_s)                                      * theta * tstep &

                   + TG_num8 * 0.25d0 / BigR * 2.d0*vpar0*vpar                                                        &
                              * (rn0_x * ps0_y - rn0_y * ps0_x + F0 / BigR * rn0_p)                                   &
                              * ( v_x * ps0_y -  v_y * ps0_x                   ) * xjac * theta * tstep * tstep 

         amat_87_k = + TG_num8 * 0.25d0 / BigR * 2.d0*vpar0*vpar                                                      &
                              * (rn0_x * ps0_y - rn0_y * ps0_x + F0 / BigR * rn0_p)                                   &
                              * (                            + F0 / BigR * v_p) * xjac * theta * tstep * tstep 
                            
         amat_87_n = + v * rn0 * F0 / BigR * vpar_p                         * xjac * theta * tstep 
          
	 amat_88  = + v * rhon * BigR * xjac * (1.d0 + zeta)                                                          &
                    - v * BigR**2 * ( rhon_s * u0_t - rhon_t * u0_s)                                  * theta * tstep &
                    - v * 2.d0 * BigR * rhon * u0_y                                            * xjac * theta * tstep &
                    + v * Vpar0 * (rhon_s * ps0_t - rhon_t * ps0_s)                                   * theta * tstep &
                    + v * rhon * (vpar0_s * ps0_t - vpar0_t * ps0_s)                                  * theta * tstep &
                    + v * F0 / BigR * rhon * vpar0_p                                           * xjac * theta * tstep &
                     ! New diffusion scheme for impurities
                    + (D_par_imp-D_prof_imp) * BigR / BB2 * Bgrad_rho_star * Bgrad_rho_rhon               * xjac * theta * tstep &
                    + D_prof_imp * BigR  * (v_x*rhon_x + v_y*rhon_y ) * xjac * theta * tstep &
                    + TG_num8 * 0.25d0 * BigR**3 * (rhon_x * u0_y - rhon_y * u0_x)                                    &
                                                  * ( v_x  * u0_y - v_y   * u0_x) * xjac * theta * tstep * tstep      &

                    + TG_num8 * 0.25d0 / BigR * vpar0**2                                                              &
                               * (rhon_x * ps0_y - rhon_y * ps0_x )                                                   &
                               * ( v_x * ps0_y -  v_y * ps0_x   ) * xjac * theta * tstep * tstep                      &
                   + BigR * (Dn0x * rhon_x * v_x + Dn0y * rhon_y * v_y)                        * xjac * theta * tstep  &
                   + Dn_perp_num * (v_xx + v_x/BigR + v_yy)*(rhon_xx + rhon_x/BigR + rhon_yy)  * BigR * xjac * theta * tstep 
          

         amat_88_k = &
                     ! New diffusion scheme for impurities
                     + (D_par_imp-D_prof_imp) * BigR / BB2 * Bgrad_rho_k_star * Bgrad_rho_rhon         * xjac * theta * tstep &

                     + TG_num8 * 0.25d0 / BigR * vpar0**2                                                             &
                               * (rhon_x * ps0_y - rhon_y * ps0_x                  )                                  &
                               * (                              + F0 / BigR * v_p) * xjac * theta * tstep * tstep

         amat_88_n = + v * F0 / BigR * Vpar0 * rhon_p                      * xjac * theta * tstep                     &
                     ! New diffusion scheme for impurities
                     + (D_par_imp-D_prof_imp) * BigR / BB2 * Bgrad_rho_star   * Bgrad_rho_rhon_n       * xjac * theta * tstep &
                     + TG_num8 * 0.25d0 / BigR * vpar0**2                                                             &
                               * (                              + F0 / BigR * rhon_p)                                 &
                               * ( v_x * ps0_y -  v_y * ps0_x                      ) * xjac * theta * tstep * tstep

	          
         amat_88_kn = + Dn0p * rhon_p * v_p*eps_cyl**2/BigR * xjac * theta * tstep                                    &
                     ! New diffusion scheme for impurities
                      + (D_par_imp-D_prof_imp) * BigR / BB2 * Bgrad_rho_k_star * Bgrad_rho_rhon_n    * xjac * theta * tstep &
                      + D_prof_imp * BigR  * ( v_p*rhon_p * eps_cyl**2 /BigR**2 )                * xjac * theta * tstep &
                      + TG_num8 * 0.25d0 / BigR * vpar0**2                                                            &
                               * ( + F0 / BigR * rhon_p)                                                              &
                               * ( + F0 / BigR * v_p) * xjac * theta * tstep * tstep
		
		


!###################################################################################################
!#  equation 9   (electron energy equation)                                                        #
!###################################################################################################

             v_psi_x = v_xx * psi_y - v_xy * psi_x + v_x * psi_xy - v_y * psi_xx
             v_psi_y = v_xy * psi_y - v_yy * psi_x + v_x * psi_yy - v_y * psi_xy

             Bgrad_Te_psi     = ( Te0_x * psi_y - Te0_y * psi_x )  / BigR
             Bgrad_Te_T       = ( Te_x * ps0_y - Te_y * ps0_x ) / BigR 
             Bgrad_Te_T_n     = ( F0 / BigR * Te_p) / BigR

             Te_ps0_x = Te_xx * ps0_y - Te_xy * ps0_x + Te_x * ps0_xy - Te_y * ps0_xx
             Te_ps0_y = Te_xy * ps0_y - Te_yy * ps0_x + Te_x * ps0_yy - Te_y * ps0_xy

             Te0_psi_x = Te0_xx * psi_y - Te0_xy * psi_x + Te0_x * psi_xy - Te0_y * psi_xx
             Te0_psi_y = Te0_xy * psi_y - Te0_yy * psi_x + Te0_x * psi_yy - Te0_y * psi_xy

             amat_91 = - (ZK_e_par_T-ZK_e_prof)*BigR * BB2_psi / BB2**2 * Bgrad_T_star * Bgrad_Te * xjac * theta * tstep &
                       + (ZK_e_par_T-ZK_e_prof) * BigR / BB2 * Bgrad_T_star_psi        * Bgrad_Te * xjac * theta * tstep &
                       + (ZK_e_par_T-ZK_e_prof) * BigR / BB2 * Bgrad_T_star        * Bgrad_Te_psi * xjac * theta * tstep &
                       + v * (r0 + rn0 * alpha_e_bis) * Vpar0 * (Te0_s * psi_t - Te0_t * psi_s)          * theta * tstep &
                       + v * Te0 * Vpar0 * ((r0_s+rn0_s*alpha_e)*psi_t - (r0_t+rn0_t*alpha_e)*psi_s)     * theta * tstep &
                       + v * (r0 + rn0 * alpha_e) * GAMMA * Te0 * (vpar0_s * psi_t - vpar0_t * psi_s)    * theta * tstep &
!=============== The ionization potential energy term=========================
                       + (GAMMA-1.) * v * rn0 * dE_ion_dT * Vpar0 * (Te0_s * psi_t - Te0_t * psi_s)         * theta * tstep &
                       + (GAMMA-1.) * v * E_ion * Vpar0 * (rn0_s * psi_t - rn0_t * psi_s)                   * theta * tstep &
                       + (GAMMA-1.) * v * E_ion_bg *Vpar0*((r0_s-rn0_s)*psi_t - (r0_t-rn0_t)*psi_s)         * theta * tstep &
                       + (GAMMA-1.) * v * E_ion * rn0 * (vpar0_s * psi_t - vpar0_t * psi_s)                 * theta * tstep &
                       + (GAMMA-1.) * v * E_ion_bg * (r0-rn0) * (vpar0_s * psi_t - vpar0_t * psi_s)         * theta * tstep &

                          ! New diffusive flux of the ionization potential energy for impurities
                       - (GAMMA - 1.) * E_ion * (D_par_imp-D_prof_imp) * BigR * BB2_psi / BB2**2 * Bgrad_rho_star * (Bgrad_rhon) * xjac * theta * tstep &
                       + (GAMMA - 1.) * E_ion * (D_par_imp-D_prof_imp) * BigR / BB2          * Bgrad_rho_star_psi * (Bgrad_rhon) * xjac * theta * tstep &
                       + (GAMMA - 1.) * E_ion * (D_par_imp-D_prof_imp) * BigR / BB2          * Bgrad_rho_star * (Bgrad_rhon_psi) * xjac * theta * tstep &
                       - (GAMMA - 1.) * E_ion_bg * (D_par-D_prof) * BigR * BB2_psi / BB2**2 * Bgrad_rho_star * (Bgrad_rho-Bgrad_rhon)         * xjac * theta * tstep &
                       + (GAMMA - 1.) * E_ion_bg * (D_par-D_prof) * BigR / BB2              * Bgrad_rho_star_psi * (Bgrad_rho-Bgrad_rhon)     * xjac * theta * tstep &
                       + (GAMMA - 1.) * E_ion_bg * (D_par-D_prof) * BigR / BB2              * Bgrad_rho_star * (Bgrad_rho_psi-Bgrad_rhon_psi) * xjac * theta * tstep &
!================= End ionization potential energy ===========================
                       + ZK_par_num * (v_psi_x  * ps0_y - v_psi_y  * ps0_x + v_ps0_x * psi_y - v_ps0_y * psi_x)          &
                                    * (Te0_ps0_x * ps0_y - Te0_ps0_y * ps0_x)                     * xjac * theta * tstep &
                       + ZK_par_num * (Te0_psi_x * ps0_y - Te0_psi_y * ps0_x + Te0_ps0_x * psi_y - Te0_ps0_y * psi_x)    &
                                    * (v_ps0_x  * ps0_y - v_ps0_y  * ps0_x)                       * xjac * theta * tstep &

                       + TG_num9 * 0.25d0 / BigR * vpar0**2                                                       &
                                 * Te0 * ((r0_x+alpha_e*rn0_x) * psi_y - (r0_y+alpha_e*rn0_y) * psi_x)            &
                                 * ( v_x * ps0_y -  v_y * ps0_x                  ) * xjac * theta * tstep * tstep &

                       + TG_num9 * 0.25d0 / BigR * vpar0**2                                                       &
                                 * (r0+alpha_e_bis*rn0) * (Te0_x * psi_y - Te0_y * psi_x)                         &
                                 * ( v_x * ps0_y -  v_y * ps0_x                  ) * xjac * theta * tstep * tstep &

                       + TG_num9 * 0.25d0 / BigR * vpar0**2                                                       &
                                 * Te0 * ((r0_x+alpha_e*rn0_x) * ps0_y - (r0_y+alpha_e*rn0_y) * ps0_x             &
                                         + F0 / BigR * (r0_p+alpha_e*rn0_p))                                      &
                                 * ( v_x * psi_y -  v_y * psi_x ) * xjac * theta * tstep * tstep                  &

                       + TG_num9 * 0.25d0 / BigR * vpar0**2                                                       &
                                 * (r0+alpha_e_bis*rn0) * (Te0_x * ps0_y - Te0_y * ps0_x + F0 / BigR * Te0_p)     &
                                 * ( v_x * psi_y -  v_y * psi_x ) * xjac * theta * tstep * tstep

             amat_91_k = - (ZK_e_par_T-ZK_e_prof)*BigR*BB2_psi / BB2**2 * Bgrad_T_k_star * Bgrad_Te * xjac * theta * tstep &
                         + (ZK_e_par_T-ZK_e_prof)*BigR / BB2 * Bgrad_T_k_star        * Bgrad_Te_psi * xjac * theta * tstep &
!=============== The ionization potential energy term=========================
                       ! New diffusive flux of the ionization potential energy for impurities
                         - (GAMMA - 1.) * E_ion * (D_par_imp-D_prof_imp) * BigR * BB2_psi / BB2**2 * Bgrad_rho_k_star * (Bgrad_rhon) * xjac * theta * tstep &
                         + (GAMMA - 1.) * E_ion * (D_par_imp-D_prof_imp) * BigR / BB2          * Bgrad_rho_k_star * (Bgrad_rhon_psi) * xjac * theta * tstep &
                         - (GAMMA - 1.) * E_ion_bg * (D_par-D_prof) * BigR * BB2_psi / BB2**2 * Bgrad_rho_k_star * (Bgrad_rho-Bgrad_rhon)         * xjac * theta * tstep &
                         + (GAMMA - 1.) * E_ion_bg * (D_par-D_prof) * BigR / BB2              * Bgrad_rho_k_star * (Bgrad_rho_psi-Bgrad_rhon_psi) * xjac * theta * tstep &
!================= End ionization potential energy ===========================
                         + TG_num9 * 0.25d0 / BigR * vpar0**2                                                       &
                                   * Te0 * ((r0_x+alpha_e*rn0_x) * psi_y - (r0_y+alpha_e*rn0_y) * psi_x)            &
                                   * (                            + F0 / BigR * v_p) * xjac * theta * tstep * tstep &

                         + TG_num9 * 0.25d0 / BigR * vpar0**2                                                       &
                                   * (r0+alpha_e_bis*rn0) * (Te0_x * psi_y - Te0_y * psi_x)                         &
                                   * (                            + F0 / BigR * v_p) * xjac * theta * tstep * tstep 

             amat_92 = - v * (r0 + rn0 * alpha_e_bis) * BigR**2 * ( Te0_s * u_t - Te0_t * u_s)       * theta * tstep &
                       - v * Te0 * BigR**2 * ((r0_s+rn0_s*alpha_e)*u_t - (r0_t+rn0_t*alpha_e)*u_s)   * theta * tstep &
                       - v * (r0 + rn0 * alpha_e) * 2.d0* GAMMA * BigR * Te0 * u_y            * xjac * theta * tstep &
!=============== The ionization potential energy term=========================
                       - (GAMMA-1.) * v * rn0 * dE_ion_dT * BigR**2 * ( Te0_s * u_t - Te0_t * u_s)        * theta * tstep &
                       - (GAMMA-1.) * v * E_ion * BigR**2 * (rn0_s * u_t - rn0_t * u_s)                   * theta * tstep &
                       - (GAMMA-1.) * v * E_ion_bg * BigR**2 *((r0_s-rn0_s)*u_t - (r0_t-rn0_t)*u_s)       * theta * tstep &
                       - (GAMMA-1.) * v * E_ion * rn0 * 2.d0 * BigR * u_y                          * xjac * theta * tstep &
                       - (GAMMA-1.) * v * E_ion_bg * (r0-rn0) * 2.d0 * BigR * u_y                  * xjac * theta * tstep &
!================= End ionization potential energy ===========================

                       + TG_num9 * 0.25d0 * BigR**2 * Te0* ((r0_x+alpha_e*rn0_x) * u_y - (r0_y+alpha_e*rn0_y) * u_x) &
                                 * ( v_x * u0_y - v_y * u0_x) * xjac * theta*tstep*tstep  &

                       + TG_num9 * 0.25d0 * BigR**2 * (r0+alpha_e_bis*rn0) * (Te0_x * u_y - Te0_y * u_x)             &
                                 * ( v_x * u0_y - v_y * u0_x) * xjac * theta*tstep*tstep  &

                       + TG_num9 * 0.25d0 * BigR**2 * Te0* ((r0_x+alpha_e*rn0_x)*u0_y - (r0_y+alpha_e*rn0_y)*u0_x)   &
                                 * ( v_x * u_y - v_y * u_x) * xjac * theta*tstep*tstep    &

                       + TG_num9 * 0.25d0 * BigR**2 * (r0+alpha_e_bis*rn0)* (Te0_x * u0_y - Te0_y * u0_x)            &
                                 * ( v_x * u_y - v_y * u_x) * xjac * theta*tstep*tstep

             amat_93 = - v * BigR * zj * 2. * ((GAMMA-1.)/BigR**2) * eta_T_ohm * zj0        * xjac * theta * tstep

             amat_95 =   v * rho * Te0   * BigR * xjac * (1.d0 + zeta)    &

                       - v * rho * BigR**2 * ( Te0_s  * u0_t - Te0_t  * u0_s)                      * theta * tstep &
                       - v * Te0  * BigR**2 * ( rho_s * u0_t - rho_t * u0_s)                       * theta * tstep &

                       - v * rho * 2.d0* GAMMA * BigR * Te0 * u0_y                          * xjac * theta * tstep &

                       + v * rho * F0 / BigR * Vpar0 * Te0_p                                * xjac * theta * tstep &

                       + v * rho * Vpar0 * (Te0_s * ps0_t - Te0_t * ps0_s)                         * theta * tstep &
                       + v * Te0 * Vpar0 * (rho_s * ps0_t - rho_t * ps0_s)                         * theta * tstep &

                       + v * rho * GAMMA * Te0 * (vpar0_s * ps0_t - vpar0_t * ps0_s)               * theta * tstep &
                       + v * rho * GAMMA * Te0 * F0 / BigR * vpar0_p                        * xjac * theta * tstep &

                       ! Energy exchange term
                       - v * BigR * ddTe_i_drho * rho                                       * xjac * theta * tstep &

                       + TG_num9 * 0.25d0 * BigR**2 * Te0* (rho_x * u0_y - rho_y * u0_x)     &
                                 * ( v_x * u0_y - v_y * u0_x) * xjac * theta*tstep*tstep     &

                       + TG_num9 * 0.25d0 * BigR**2 * rho * (Te0_x * u0_y - Te0_y * u0_x)    &
                                 * ( v_x * u0_y - v_y * u0_x) * xjac* theta*tstep*tstep      &

                       + TG_num9 * 0.25d0 / BigR * vpar0**2 &
                          * Te0 * (rho_x * ps0_y - rho_y * ps0_x                    )                     &
                          * ( v_x * ps0_y -  v_y * ps0_x                  ) * xjac * theta * tstep * tstep&

                       + TG_num9 * 0.25d0 / BigR * vpar0**2 &
                          * rho * (Te0_x * ps0_y - Te0_y * ps0_x + F0 / BigR * Te0_p)                      &
                          * ( v_x * ps0_y -  v_y * ps0_x                  ) * xjac * theta * tstep * tstep &

                       + v * BigR * rho * dr0_corr_dn * rn0_corr * Lrad                     * xjac * theta * tstep &
                       + v * BigR * rho * dr0_corr_dn * frad_bg                             * xjac * theta * tstep &
                      ! New term from Z_eff
                       - v * BigR * rho * ((GAMMA-1.)/BigR**2) * deta_dr0_ohm * zj0**2      * xjac * theta * tstep &
!=============== The ionization potential energy term=========================
                       + (GAMMA-1.) * v * rho * E_ion_bg * BigR * xjac * (1.d0 + zeta)  &
                       - (GAMMA-1.) * v * E_ion_bg * BigR**2 * (rho_s * u0_t - rho_t * u0_s) * theta * tstep &

                       + (GAMMA-1.) * v * E_ion_bg * Vpar0 * (rho_s * ps0_t - rho_t * ps0_s) * theta * tstep &

                       - (GAMMA-1.) * v * E_ion_bg * rho * 2.d0 * BigR * u0_y         * xjac * theta * tstep &
                       + (GAMMA-1.) * v * E_ion_bg * rho * (vpar0_s*ps0_t - vpar0_t*ps0_s)   * theta * tstep &
                       + (GAMMA-1.) * v * E_ion_bg * rho * F0 / BigR * vpar0_p        * xjac * theta * tstep &

                    ! New diffusive ionization energy flux term
                       + (GAMMA - 1.) * E_ion_bg * (D_par-D_prof) * BigR / BB2 * Bgrad_rho_star * Bgrad_rho_rho                * xjac * theta * tstep &
                       + (GAMMA - 1.) * E_ion_bg * D_prof * BigR  * (v_x*rho_x + v_y*rho_y                                   ) * xjac * theta * tstep
!================= End ionization potential energy ===========================

             amat_95_n = &
!=============== The ionization potential energy term=========================
                       + (GAMMA - 1.) * v * E_ion_bg * F0 / BigR * Vpar0 * rho_p                                  * xjac * theta * tstep &
                    ! New diffusive ionization energy flux term
                       + (GAMMA - 1.) * E_ion_bg * (D_par-D_prof) * BigR / BB2 * Bgrad_rho_star * Bgrad_rho_rho_n * xjac * theta * tstep &
!================= End ionization potential energy =========================== 
                       + v * Te0 * F0 / BigR * Vpar0 * rho_p                                * xjac * theta * tstep &
                       + TG_num9 * 0.25d0 / BigR * vpar0**2 &
                           * Te0 * (                                + F0 / BigR * rho_p)                     &
                           * ( v_x * ps0_y -  v_y * ps0_x                     ) * xjac * theta * tstep * tstep 

             amat_95_k = &
!=============== The ionization potential energy term=========================
                    ! New diffusive ionization energy flux term
                       + (GAMMA - 1.) * E_ion_bg * (D_par-D_prof) * BigR / BB2 * Bgrad_rho_k_star * Bgrad_rho_rho * xjac * theta * tstep &

!================= End ionization potential energy ===========================

                       + TG_num9 * 0.25d0 / BigR * vpar0**2 &
                          * Te0 * (rho_x * ps0_y - rho_y * ps0_x                    )                     &
                          * (                            + F0 / BigR * v_p) * xjac * theta * tstep * tstep&

                       + TG_num9 * 0.25d0 / BigR * vpar0**2 &
                          * rho * (Te0_x * ps0_y - Te0_y * ps0_x + F0 / BigR * Te0_p)                      &
                          * (                            + F0 / BigR * v_p) * xjac * theta * tstep * tstep

             amat_95_kn = &
!=============== The ionization potential energy term=========================
                    ! New diffusive ionization energy flux term
                    + (GAMMA - 1.) * E_ion_bg * (D_par-D_prof) * BigR / BB2 * Bgrad_rho_k_star * Bgrad_rho_rho_n            * xjac * theta * tstep &
                    + (GAMMA - 1.) * E_ion_bg * D_prof * BigR  * (                      + v_p*rho_p * eps_cyl**2 /BigR**2 ) * xjac * theta * tstep &

!================= End ionization potential energy ===========================
                       + TG_num9 * 0.25d0 / BigR * vpar0**2 &
                          * Te0 * (                              + F0 / BigR * rho_p)                     &
                          * (                            + F0 / BigR * v_p) * xjac * theta * tstep * tstep

             amat_96 = - v * BigR * ddTe_i_dTi * Ti                                     * xjac * theta * tstep

             amat_97 = + v * (r0 + rn0 * alpha_e_bis) * F0 / BigR * Vpar * Te0_p        * xjac * theta * tstep &
                       + v * Te0 * F0 / BigR * Vpar * (r0_p + rn0_p * alpha_e)          * xjac * theta * tstep &

                       + v * (r0 + rn0 * alpha_e_bis) * Vpar * (Te0_s * ps0_t - Te0_t * ps0_s)         * theta * tstep &
                       + v * Te0 * Vpar * ((r0_s+rn0_s*alpha_e)*ps0_t - (r0_t+rn0_t*alpha_e)*ps0_s)    * theta * tstep &

                       + v * (r0 + rn0 * alpha_e) * GAMMA * Te0 * (vpar_s * ps0_t - vpar_t * ps0_s)    * theta * tstep &
!=============== The ionization potential energy term=========================
                       + (GAMMA-1.) * v * rn0 * dE_ion_dT * F0 / BigR * Vpar * Te0_p       * xjac * theta * tstep  &
                       + (GAMMA-1.) * v * E_ion * F0 / BigR * Vpar * rn0_p                 * xjac * theta * tstep  &
                       + (GAMMA-1.) * v * E_ion_bg * F0 / BigR * Vpar * (r0_p-rn0_p)       * xjac * theta * tstep  &

                       + (GAMMA-1.) * v * rn0 * dE_ion_dT * Vpar * (Te0_s * ps0_t - Te0_t * ps0_s)* theta * tstep  &
                       + (GAMMA-1.) * v * E_ion * Vpar * (rn0_s * ps0_t - rn0_t * ps0_s)          * theta * tstep  &
                       + (GAMMA-1.) * v * E_ion_bg*Vpar*((r0_s-rn0_s)*ps0_t - (r0_t-rn0_t)*ps0_s) * theta * tstep  &

                       + (GAMMA-1.) * v * E_ion * rn0 * (vpar_s * ps0_t - vpar_t * ps0_s)         * theta * tstep  &

                       + (GAMMA-1.) * v * E_ion_bg * (r0-rn0) * (vpar_s * ps0_t - vpar_t * ps0_s) * theta * tstep  &
!================= End ionization potential energy ===========================

                           + TG_num9 * 0.25d0 / BigR * 2.d0 * vpar0*vpar &
                              * Te0 * ((r0_x+alpha_e*rn0_x)*ps0_y - (r0_y+alpha_e*rn0_y)*ps0_x                 &
                                      + F0 / BigR * (r0_p+alpha_e*rn0_p))                                      &
                              * ( v_x * ps0_y -  v_y * ps0_x                  ) * xjac * theta * tstep * tstep &

                           + TG_num9 * 0.25d0 / BigR * 2.d0 * vpar0*vpar &
                              * (r0+alpha_e_bis*rn0) * (Te0_x * ps0_y - Te0_y * ps0_x + F0 / BigR * Te0_p)     &
                              * ( v_x * ps0_y -  v_y * ps0_x                  ) * xjac * theta * tstep * tstep

             amat_97_k =  &
                           + TG_num9 * 0.25d0 / BigR * 2.d0 * vpar0*vpar &
                              * Te0 * ((r0_x+alpha_e*rn0_x)*ps0_y - (r0_y+alpha_e*rn0_y)*ps0_x                 &
                                      + F0 / BigR * (r0_p+alpha_e*rn0_p))                                      &
                              * (                            + F0 / BigR * v_p) * xjac * theta * tstep * tstep &

                           + TG_num9 * 0.25d0 / BigR * 2.d0 * vpar0*vpar &
                              * (r0+alpha_e_bis*rn0) * (Te0_x * ps0_y - Te0_y * ps0_x + F0 / BigR * Te0_p)     &
                              * (                            + F0 / BigR * v_p) * xjac * theta * tstep * tstep

             amat_97_n = &
                       + v * (r0 + rn0 * alpha_e) * GAMMA * Te0 * F0 / BigR * vpar_p            * xjac * theta * tstep &
                       + (GAMMA-1.) * v * E_ion * rn0 * F0 / BigR * vpar_p                      * xjac * theta * tstep &
                       + (GAMMA-1.) * v * E_ion_bg * (r0-rn0) * F0 / BigR * vpar_p              * xjac * theta * tstep 

             amat_98 =   v * rhon * alpha_e * Te0 * BigR * xjac * (1.d0 + zeta)                              &
!=============== The ionization potential energy term=========================
                       + (GAMMA-1.) * v * rhon * (E_ion - E_ion_bg) * BigR * xjac * (1.d0 + zeta)                &

                       - (GAMMA-1.) * v * rhon * dE_ion_dT * BigR**2 * (Te0_s*u0_t - Te0_t*u0_s)  * theta * tstep&
                       - (GAMMA-1.) * v * (E_ion-E_ion_bg) * BigR**2 * (rhon_s*u0_t - rhon_t*u0_s)* theta * tstep&

                       + (GAMMA-1.) * v * rhon * dE_ion_dT * F0 / BigR * Vpar0 * Te0_p     * xjac * theta * tstep&

                       + (GAMMA-1.) * v * rhon * dE_ion_dT * Vpar0 * (Te0_s*ps0_t - Te0_t*ps0_s)  * theta * tstep&
                       + (GAMMA-1.) * v * (E_ion-E_ion_bg) * Vpar0 * (rhon_s*ps0_t - rhon_t*ps0_s)* theta * tstep&

                       - (GAMMA-1.) * v * (E_ion-E_ion_bg) * rhon * 2.d0 * BigR * u0_y     * xjac * theta * tstep&
                       + (GAMMA-1.) * v * (E_ion-E_ion_bg) * rhon*(vpar0_s*ps0_t - vpar0_t*ps0_s) * theta * tstep&
                       + (GAMMA-1.) * v * (E_ion-E_ion_bg) * rhon * F0 / BigR * vpar0_p    * xjac * theta * tstep&

                       ! New diffusive ionization energy flux term
                       + (GAMMA - 1.) * E_ion * (D_par_imp-D_prof_imp) * BigR / BB2 * Bgrad_rho_star * Bgrad_rho_rhon                  * xjac * theta * tstep &
                       + (GAMMA - 1.) * E_ion * D_prof_imp * BigR  * (v_x*rhon_x + v_y*rhon_y                                    ) * xjac * theta * tstep &
                       - (GAMMA - 1.) * E_ion_bg * (D_par-D_prof) * BigR / BB2 * Bgrad_rho_star * Bgrad_rho_rhon                  * xjac * theta * tstep &
                       - (GAMMA - 1.) * E_ion_bg * D_prof * BigR  * (v_x*rhon_x + v_y*rhon_y                                    ) * xjac * theta * tstep &
!================= End ionization potential energy ===========================
!=========================New TG_num terms====================================
                       + TG_num9 * 0.25d0 * BigR**2 * Te0 * alpha_e * (rhon_x * u0_y - rhon_y * u0_x)        &
                                 * ( v_x * u0_y - v_y * u0_x) * xjac * theta*tstep*tstep     &

                       + TG_num9 * 0.25d0 * BigR**2 * alpha_e_bis * rhon * (Te0_x * u0_y - Te0_y * u0_x)     &
                                 * ( v_x * u0_y - v_y * u0_x) * xjac* theta*tstep*tstep      &

                       + TG_num9 * 0.25d0 / BigR * vpar0**2 &
                          * Te0 * alpha_e * (rhon_x * ps0_y - rhon_y * ps0_x                     )           &
                          * ( v_x * ps0_y -  v_y * ps0_x                  ) * xjac * theta * tstep * tstep   &

                       + TG_num9 * 0.25d0 / BigR * vpar0**2 &
                          * alpha_e_bis * rhon * (Te0_x * ps0_y - Te0_y * ps0_x + F0 / BigR * Te0_p)         &
                          * ( v_x * ps0_y -  v_y * ps0_x                    ) * xjac * theta * tstep * tstep &
!===========================End of new TG_num terms===========================
                       ! New term from Z_eff
                       - v * BigR * rhon*((GAMMA-1.)/BigR**2) * deta_drn0_ohm * zj0**2 * xjac * theta * tstep &
                       - v * rhon * BigR**2 * alpha_e_bis * (Te0_s * u0_t - Te0_t * u0_s)     * theta * tstep &
                       - v * alpha_e * Te0 * BigR**2 * (rhon_s * u0_t - rhon_t * u0_s)        * theta * tstep &
                       + v * rhon * F0 / BigR * Vpar0 * alpha_e_bis * Te0_p            * xjac * theta * tstep &
                       + v * rhon * Vpar0 * alpha_e_bis * (Te0_s * ps0_t - Te0_t * ps0_s)     * theta * tstep &
                       + v * alpha_e * Te0 * Vpar0 * (rhon_s * ps0_t - rhon_t * ps0_s)        * theta * tstep &

                       - v * alpha_e * rhon * 2.d0 * GAMMA * BigR * Te0 * u0_y                  * xjac * theta * tstep &
                       + v * alpha_e * rhon * GAMMA * Te0 * (vpar0_s * ps0_t - vpar0_t * ps0_s)        * theta * tstep &
                       + v * alpha_e * rhon * GAMMA * Te0 * F0 / BigR * vpar0_p                 * xjac * theta * tstep &

                       ! Energy exchange term
                       - v * BigR * ddTe_i_drhon * rhon                                         * xjac * theta * tstep &

                       + v * BigR * rhon * drn0_corr_dn * (r0_corr + 2.*alpha_e*rn0_corr) * Lrad* xjac * theta * tstep &
                       + v * BigR * rhon * drn0_corr_dn * alpha_e * frad_bg                     * xjac * theta * tstep


             amat_98_k = &
!=============== The ionization potential energy term=========================
                       ! New diffusive ionization energy flux term
                       + (GAMMA - 1.) * E_ion * (D_par_imp-D_prof_imp) * BigR / BB2 * Bgrad_rho_k_star * Bgrad_rho_rhon                * xjac * theta * tstep &
                       - (GAMMA - 1.) * E_ion_bg * (D_par-D_prof) * BigR / BB2 * Bgrad_rho_k_star * Bgrad_rho_rhon             * xjac * theta * tstep &
!================= End ionization potential energy ===========================
                       + TG_num9 * 0.25d0 / BigR * vpar0**2 &
                          * Te0 * alpha_e * (rhon_x * ps0_y - rhon_y * ps0_x                     )           &
                          * (                            + F0 / BigR * v_p) * xjac * theta * tstep * tstep   &

                       + TG_num9 * 0.25d0 / BigR * vpar0**2 &
                          * alpha_e_bis * rhon * (Te0_x * ps0_y - Te0_y * ps0_x + F0 / BigR * Te0_p)         &
                          * (                            + F0 / BigR * v_p) * xjac * theta * tstep * tstep 

             amat_98_n = &
!=============== The ionization potential energy term=========================
                       + (GAMMA - 1.) * v * (E_ion-E_ion_bg) * F0 / BigR * Vpar0 * rhon_p   * xjac * theta * tstep &
                       ! New diffusive ionization energy flux term
                       + (GAMMA - 1.) * E_ion * (D_par_imp-D_prof_imp) * BigR / BB2 * Bgrad_rho_star * Bgrad_rho_rhon_n                * xjac * theta * tstep &
                       - (GAMMA - 1.) * E_ion_bg * (D_par-D_prof) * BigR / BB2 * Bgrad_rho_star * Bgrad_rho_rhon_n             * xjac * theta * tstep &
!================= End ionization potential energy ===========================
                       + v * alpha_e * Te0 * F0 / BigR * Vpar0 * rhon_p                    * xjac * theta * tstep &
                       + TG_num9 * 0.25d0 / BigR * vpar0**2 &
                          * Te0 * alpha_e * (                                + F0 / BigR * rhon_p)           &
                          * ( v_x * ps0_y -  v_y * ps0_x                  ) * xjac * theta * tstep * tstep

             amat_98_kn = &
!=============== The ionization potential energy term=========================
                       ! New diffusive ionization energy flux term
                       + (GAMMA - 1.) * E_ion * (D_par_imp-D_prof_imp) * BigR / BB2 * Bgrad_rho_k_star * Bgrad_rho_rhon_n              * xjac * theta * tstep &
                       + (GAMMA - 1.) * E_ion * D_prof_imp * BigR  * (                        + v_p*rhon_p * eps_cyl**2 /BigR**2 ) * xjac * theta * tstep &
                       - (GAMMA - 1.) * E_ion_bg * (D_par-D_prof) * BigR / BB2 * Bgrad_rho_k_star * Bgrad_rho_rhon_n              * xjac * theta * tstep &
                       - (GAMMA - 1.) * E_ion_bg * D_prof * BigR  * (                        + v_p*rhon_p * eps_cyl**2 /BigR**2 ) * xjac * theta * tstep &
!================= End ionization potential energy ===========================
                       + TG_num9 * 0.25d0 / BigR * vpar0**2 &
                          * Te0 * alpha_e * (                                + F0 / BigR * rhon_p)           &
                          * (                            + F0 / BigR * v_p) * xjac * theta * tstep * tstep

             amat_99 =   v * (r0 + rn0 * alpha_e_bis) * Te * BigR * xjac * (1.d0 + zeta)                     &
!=============== The ionization potential energy term=========================
                       + (GAMMA-1.) * v * rn0 * dE_ion_dT  * Te * BigR * xjac * (1.d0 + zeta)                &
                       - (GAMMA-1.) * v * rn0 * dE_ion_dT * BigR**2 * (Te_s*u0_t - Te_t*u0_s)  * theta * tstep &

                       + (GAMMA-1.) * v * rn0 * dE_ion_dT * Vpar0 * (Te_s*ps0_t - Te_t*ps0_s)  * theta * tstep &

                       ! New diffusive ionization energy flux term
                       + (GAMMA - 1.) * dE_ion_dT * Te * (D_par_imp-D_prof_imp) * BigR / BB2 * Bgrad_rho_star * (Bgrad_rhon)                       * xjac * tstep &
                       + (GAMMA - 1.) * dE_ion_dT * Te * D_prof_imp * BigR  * (v_x*(rn0_x) + v_y*(rn0_y)                                     ) * xjac * tstep &
!================= End ionization potential energy ===========================
                       - v * (r0 + rn0 * alpha_e_bis) * BigR**2 * (Te_s * u0_t - Te_t  * u0_s) * theta * tstep &
                       - v * (rn0 * alpha_e_tri) * Te * BigR**2 * (Te0_s* u0_t - Te0_t * u0_s) * theta * tstep &
                       - v * Te  * BigR**2 * (r0_s * u0_t - r0_t * u0_s)                       * theta * tstep &
                       - v * alpha_e_bis * Te * BigR**2 * (rn0_s * u0_t - rn0_t * u0_s)        * theta * tstep &

                       - v * (r0 + rn0 * alpha_e_bis) * 2.d0* GAMMA * BigR * Te * u0_y  * xjac * theta * tstep &

                       + v * Te * F0  / BigR * Vpar0 * (r0_p + rn0_p * alpha_e_bis)     * xjac * theta * tstep &
                       + v * (rn0 * alpha_e_tri) * Te * F0 / BigR * Vpar0 * Te0_p       * xjac * theta * tstep &

                       + v * (r0 + rn0 * alpha_e_bis) * Vpar0 * (Te_s  * ps0_t - Te_t  * ps0_s)* theta * tstep &
                       + v * (rn0 * alpha_e_tri) * Te * Vpar0 * (Te0_s * ps0_t - Te0_t * ps0_s)* theta * tstep &
                       + v * Te  * Vpar0 * (r0_s * ps0_t - r0_t * ps0_s)                       * theta * tstep &
                       + v * alpha_e_bis * Te * Vpar0 * (rn0_s * ps0_t - rn0_t * ps0_s)        * theta * tstep &

                       + v * (r0 + rn0 * alpha_e_bis) * GAMMA * Te * (vpar0_s * ps0_t - vpar0_t * ps0_s)* theta * tstep &
                       + v * (r0 + rn0 * alpha_e_bis) * GAMMA * Te * F0 / BigR * vpar0_p         * xjac * theta * tstep &

                       ! Energy exchange term
                       - v * BigR * ddTe_i_dTe * Te                                              * xjac * theta * tstep &

                       + (ZK_e_par_T-ZK_e_prof) * BigR / BB2 * Bgrad_T_star * Bgrad_Te_T         * xjac * theta * tstep &
                       + ZK_e_prof * BigR * (v_x*Te_x + v_y*Te_y                     )           * xjac * theta * tstep &

                       + dZK_e_par_dT * Te * BigR / BB2 * Bgrad_T_star * Bgrad_Te                * xjac * theta * tstep &

                       + ZK_perp_num*(v_xx + v_x/BigR + v_yy)*(Te_xx + Te_x/BigR + Te_yy) * BigR * xjac * theta * tstep &

                       + TG_num9 * 0.25d0 * BigR**2 * Te * ((r0_x+alpha_e_bis*rn0_x)*u0_y &
                                                            - (r0_y+alpha_e_bis*rn0_y)*u0_x)   &
                                 * ( v_x * u0_y - v_y * u0_x) * xjac * theta * tstep * tstep &

                       + TG_num9 * 0.25d0 * BigR**2 * (r0+alpha_e_bis*rn0)* (Te_x * u0_y - Te_y * u0_x)  &
                                 * ( v_x * u0_y - v_y * u0_x) * xjac * theta * tstep * tstep &

                       + TG_num9 * 0.25d0 * BigR**2 * (alpha_e_tri*rn0)*Te* (Te0_x* u0_y - Te0_y* u0_x)         &
                             * ( v_x * u0_y - v_y * u0_x) * xjac * theta * tstep * tstep &

                       + TG_num9 * 0.25d0 / BigR * vpar0**2 &
                          * Te * ((r0_x+alpha_e_bis*rn0_x)*ps0_y - (r0_y+alpha_e_bis*rn0_y)*ps0_x         &
                                 + F0/BigR*(r0_p+alpha_e_bis*rn0_p))                                      &
                          * ( v_x * ps0_y -  v_y * ps0_x                  ) * xjac * theta * tstep * tstep&

                       + TG_num9 * 0.25d0 / BigR * vpar0**2 &
                          * (r0+alpha_e_bis*rn0) * (Te_x * ps0_y - Te_y * ps0_x                   )        &
                          * ( v_x * ps0_y -  v_y * ps0_x                  ) * xjac * theta * tstep * tstep &

                       + TG_num9 * 0.25d0 / BigR * vpar0**2 &
                                 * (alpha_e_tri*rn0)* Te * (Te0_x * ps0_y - Te0_y * ps0_x + F0 / BigR * Te0_p)      &
                                 * ( v_x * ps0_y -  v_y * ps0_x                  ) * xjac * theta * tstep * tstep &

                       - v * BigR * Te * ((GAMMA-1.)/BigR**2) * deta_dT_ohm * zj0**2        * xjac * theta * tstep  &
                       + v * BigR * Te * (r0_corr + alpha_e*rn0_corr) * rn0_corr * dLrad_dT * xjac * theta * tstep  &
                       + v * BigR * Te * dalpha_e_dT * rn0_corr**2 * Lrad                   * xjac * theta * tstep  &
                       + v * BigR * Te * (r0_corr + alpha_e*rn0_corr) * dfrad_bg_dT         * xjac * theta * tstep  &
                       + v * BigR * Te * dalpha_e_dT * rn0_corr * frad_bg                   * xjac * theta * tstep


             amat_99_k = &
                       + (ZK_e_par_T-ZK_e_prof) * BigR / BB2 * Bgrad_T_k_star * Bgrad_Te_T   * xjac * theta * tstep &
                       + dZK_e_par_dT * Te * BigR / BB2 * Bgrad_T_k_star * Bgrad_Te          * xjac * theta * tstep &
!=============== The ionization potential energy term=========================
                       ! New diffusive ionization energy flux term
                       + (GAMMA - 1.) * dE_ion_dT * Te * (D_par_imp-D_prof_imp) * BigR / BB2 * Bgrad_rho_k_star * (Bgrad_rhon)                     * xjac * tstep &
                       + (GAMMA - 1.) * dE_ion_dT * Te * D_prof_imp * BigR  * (                          + v_p*(rn0_p) * eps_cyl**2 /BigR**2 ) * xjac * tstep &

!================= End ionization potential energy ===========================
                       + TG_num9 * 0.25d0 / BigR * vpar0**2 &
                          * Te * ((r0_x+alpha_e_bis*rn0_x)*ps0_y - (r0_y+alpha_e_bis*rn0_y)*ps0_x         &
                                 + F0/BigR*(r0_p+alpha_e_bis*rn0_p))                                      &
                          * (                            + F0 / BigR * v_p) * xjac * theta * tstep * tstep&

                       + TG_num9 * 0.25d0 / BigR * vpar0**2 &
                                 * (alpha_e_tri*rn0)* Te * (Te0_x * ps0_y - Te0_y * ps0_x + F0 / BigR * Te0_p)    &
                                 * (                            + F0 / BigR * v_p) * xjac * theta * tstep * tstep &

                       + TG_num9 * 0.25d0 / BigR * vpar0**2 &
                          * (r0+alpha_e_bis*rn0) * (Te_x * ps0_y - Te_y * ps0_x                  )        &
                          * (                            + F0 / BigR * v_p) * xjac * theta * tstep * tstep 

             amat_99_n = &
                       + (ZK_e_par_T-ZK_e_prof) * BigR / BB2 * Bgrad_T_star * Bgrad_Te_T_n   * xjac * theta * tstep &
                       + (GAMMA-1.) * v * rn0 * dE_ion_dT * F0 / BigR * Vpar0 * Te_p         * xjac * theta * tstep &
                       + v * (r0 + rn0 * alpha_e_bis) * F0 / BigR * Vpar0 * Te_p             * xjac * theta * tstep &
                       + TG_num9 * 0.25d0 / BigR * vpar0**2 &
                          * (r0+alpha_e_bis*rn0) * (                                 + F0 / BigR * Te_p)        &
                          * ( v_x * ps0_y -  v_y * ps0_x                           ) * xjac * theta * tstep * tstep

             amat_99_kn = &
                       + (ZK_e_par_T-ZK_e_prof) * BigR / BB2 * Bgrad_T_k_star * Bgrad_Te_T_n * xjac * theta * tstep &
                       + ZK_e_prof * BigR * (                    + v_p*Te_p /BigR**2 )       * xjac * theta * tstep &
                       + TG_num9 * 0.25d0 / BigR * vpar0**2 &
                          * (r0+alpha_e_bis*rn0) * (                                  + F0 / BigR * Te_p)        &
                          * (                                      + F0 / BigR * v_p) * xjac * theta * tstep * tstep


!###################################################################################################
!# end equation 9                                                                                  #
!###################################################################################################

             kl1 = index_kl
             kl2 = index_kl + 1
             kl3 = index_kl + 2
             kl4 = index_kl + 3
             kl5 = index_kl + 4
             kl6 = index_kl + 5
             kl7 = index_kl + 6
     	     kl8 = index_kl + 7
             kl9 = index_kl + 8


             ELM_p(mp,ij1,kl1)  =  ELM_p(mp,ij1,kl1) + wst * amat_11
             ELM_n(mp,ij1,kl2)  =  ELM_n(mp,ij1,kl2) + wst * amat_12_n

             ELM_p(mp,ij1,kl2)  =  ELM_p(mp,ij1,kl2) + wst * amat_12
             ELM_p(mp,ij1,kl3)  =  ELM_p(mp,ij1,kl3) + wst * amat_13
             ELM_p(mp,ij1,kl5)  =  ELM_p(mp,ij1,kl5) + wst * amat_15
             ELM_n(mp,ij1,kl5)  =  ELM_n(mp,ij1,kl5) + wst * amat_15_n
             ELM_p(mp,ij1,kl6)  =  ELM_p(mp,ij1,kl6) + wst * amat_16
             ELM_n(mp,ij1,kl6)  =  ELM_n(mp,ij1,kl6) + wst * amat_16_n
             ELM_p(mp,ij1,kl8)  =  ELM_p(mp,ij1,kl8) + wst * amat_18 ! New Z_eff term
             ELM_p(mp,ij1,kl9)  =  ELM_p(mp,ij1,kl9) + wst * amat_19

             ELM_p(mp,ij2,kl1)  =  ELM_p(mp,ij2,kl1) + wst * amat_21
             ELM_p(mp,ij2,kl2)  =  ELM_p(mp,ij2,kl2) + wst * amat_22
             ELM_p(mp,ij2,kl3)  =  ELM_p(mp,ij2,kl3) + wst * amat_23
             ELM_n(mp,ij2,kl3)  =  ELM_n(mp,ij2,kl3) + wst * amat_23_n

             ELM_p(mp,ij2,kl4)  =  ELM_p(mp,ij2,kl4) + wst * amat_24
             ELM_p(mp,ij2,kl5)  =  ELM_p(mp,ij2,kl5) + wst * amat_25
             ELM_n(mp,ij2,kl5)  =  ELM_n(mp,ij2,kl5) + wst * amat_25_n
             ELM_p(mp,ij2,kl6)  =  ELM_p(mp,ij2,kl6) + wst * amat_26
             ELM_p(mp,ij2,kl7)  =  ELM_p(mp,ij2,kl7) + wst * amat_27  
             ELM_n(mp,ij2,kl7)  =  ELM_n(mp,ij2,kl7) + wst * amat_27_n
     	     ELM_p(mp,ij2,kl8)  =  ELM_p(mp,ij2,kl8) + wst * amat_28
             ELM_p(mp,ij2,kl9)  =  ELM_p(mp,ij2,kl9) + wst * amat_29

             ELM_p(mp,ij3,kl1)  =  ELM_p(mp,ij3,kl1) + wst * amat_31
             ELM_p(mp,ij3,kl3)  =  ELM_p(mp,ij3,kl3) + wst * amat_33

             ELM_p(mp,ij4,kl2)  =  ELM_p(mp,ij4,kl2) + wst * amat_42
             ELM_p(mp,ij4,kl4)  =  ELM_p(mp,ij4,kl4) + wst * amat_44


             ELM_p(mp,ij5,kl1)  =  ELM_p(mp,ij5,kl1)  + wst * amat_51
             ELM_k(mp,ij5,kl1)  =  ELM_k(mp,ij5,kl1)  + wst * amat_51_k

             ELM_p(mp,ij5,kl2)  =  ELM_p(mp,ij5,kl2)  + wst * amat_52

             ELM_p(mp,ij5,kl5)  =  ELM_p(mp,ij5,kl5)  + wst * amat_55
             ELM_k(mp,ij5,kl5)  =  ELM_k(mp,ij5,kl5)  + wst * amat_55_k
             ELM_n(mp,ij5,kl5)  =  ELM_n(mp,ij5,kl5)  + wst * amat_55_n
             ELM_kn(mp,ij5,kl5) =  ELM_kn(mp,ij5,kl5) + wst * amat_55_kn

             ELM_p(mp,ij5,kl6)  =  ELM_p(mp,ij5,kl6)  + wst * amat_56

             ELM_p(mp,ij5,kl7)  =  ELM_p(mp,ij5,kl7)  + wst * amat_57
             ELM_n(mp,ij5,kl7)  =  ELM_n(mp,ij5,kl7)  + wst * amat_57_n
             ELM_k(mp,ij5,kl7)  =  ELM_k(mp,ij5,kl7)  + wst * amat_57_k
             ELM_kn(mp,ij5,kl7) =  ELM_kn(mp,ij5,kl7) + wst * amat_57_kn
	     
             ELM_p(mp,ij5,kl8)  =  ELM_p(mp,ij5,kl8)  + wst * amat_58
             ELM_kn(mp,ij5,kl8) =  ELM_kn(mp,ij5,kl8) + wst * amat_58_kn !This term was not implemented previously...

             ELM_p(mp,ij6,kl1)  =  ELM_p(mp,ij6,kl1)  + wst * amat_61
             ELM_k(mp,ij6,kl1)  =  ELM_k(mp,ij6,kl1)  + wst * amat_61_k

             ELM_p(mp,ij6,kl2)  =  ELM_p(mp,ij6,kl2)  + wst * amat_62
             ELM_p(mp,ij6,kl3)  =  ELM_p(mp,ij6,kl3)  + wst * amat_63
             ELM_p(mp,ij6,kl5)  =  ELM_p(mp,ij6,kl5)  + wst * amat_65
             ELM_n(mp,ij6,kl5)  =  ELM_n(mp,ij6,kl5)  + wst * amat_65_n
             ELM_k(mp,ij6,kl5)  =  ELM_k(mp,ij6,kl5)  + wst * amat_65_k
             ELM_kn(mp,ij6,kl5) =  ELM_kn(mp,ij6,kl5) + wst * amat_65_kn

             ELM_p(mp,ij6,kl6)  =  ELM_p(mp,ij6,kl6)  + wst * amat_66
             ELM_k(mp,ij6,kl6)  =  ELM_k(mp,ij6,kl6)  + wst * amat_66_k
             ELM_n(mp,ij6,kl6)  =  ELM_n(mp,ij6,kl6)  + wst * amat_66_n
             ELM_kn(mp,ij6,kl6) =  ELM_kn(mp,ij6,kl6) + wst * amat_66_kn

             ELM_p(mp,ij6,kl7)  =  ELM_p(mp,ij6,kl7)  + wst * amat_67
             ELM_k(mp,ij6,kl7)  =  ELM_k(mp,ij6,kl7)  + wst * amat_67_k
             ELM_n(mp,ij6,kl7)  =  ELM_n(mp,ij6,kl7)  + wst * amat_67_n

             ELM_p(mp,ij6,kl8)  =  ELM_p(mp,ij6,kl8)  + wst * amat_68
             ELM_n(mp,ij6,kl8)  =  ELM_n(mp,ij6,kl8)  + wst * amat_68_n ! Not implemented previously
             ELM_k(mp,ij6,kl8)  =  ELM_k(mp,ij6,kl8)  + wst * amat_68_k ! Not implemented previously
             ELM_kn(mp,ij6,kl8) =  ELM_kn(mp,ij6,kl8) + wst * amat_68_kn ! Not implemented previously

             ELM_p(mp,ij6,kl9)  =  ELM_p(mp,ij6,kl9)  + wst * amat_69

             ELM_p(mp,ij7,kl1)  =  ELM_p(mp,ij7,kl1)  + wst * amat_71

             ELM_p(mp,ij7,kl2)  =  ELM_p(mp,ij7,kl2)  + wst * amat_72  ! NEW

             ELM_p(mp,ij7,kl5)  =  ELM_p(mp,ij7,kl5)  + wst * amat_75
             ELM_n(mp,ij7,kl5)  =  ELM_n(mp,ij7,kl5)  + wst * amat_75_n
             ELM_k(mp,ij7,kl5)  =  ELM_k(mp,ij7,kl5)  + wst * amat_75_k
             ELM_kn(mp,ij7,kl5) =  ELM_kn(mp,ij7,kl5) + wst * amat_75_kn ! New TG term

             ELM_p(mp,ij7,kl6)  =  ELM_p(mp,ij7,kl6)  + wst * amat_76
             ELM_n(mp,ij7,kl6)  =  ELM_n(mp,ij7,kl6)  + wst * amat_76_n

             ELM_p(mp,ij7,kl7)  =  ELM_p(mp,ij7,kl7)  + wst * amat_77
             ELM_k(mp,ij7,kl7)  =  ELM_k(mp,ij7,kl7)  + wst * amat_77_k
             ELM_n(mp,ij7,kl7)  =  ELM_n(mp,ij7,kl7)  + wst * amat_77_n ! New term due to the new momentum eq.
             ELM_kn(mp,ij7,kl7) =  ELM_kn(mp,ij7,kl7) + wst * amat_77_kn
     	     ELM_p(mp,ij7,kl8)  =  ELM_p(mp,ij7,kl8)  + wst * amat_78
             ELM_n(mp,ij7,kl8)  =  ELM_n(mp,ij7,kl8)  + wst * amat_78_n
             ELM_p(mp,ij7,kl9)  =  ELM_p(mp,ij7,kl9)  + wst * amat_79
             ELM_n(mp,ij7,kl9)  =  ELM_n(mp,ij7,kl9)  + wst * amat_79_n

             ELM_p(mp,ij8,kl1)  =  ELM_p(mp,ij8,kl1)  + wst * amat_81
             ELM_k(mp,ij8,kl1)  =  ELM_k(mp,ij8,kl1)  + wst * amat_81_k ! Not implemented previously
             ELM_p(mp,ij8,kl2)  =  ELM_p(mp,ij8,kl2)  + wst * amat_82
             ELM_p(mp,ij8,kl5)  =  ELM_p(mp,ij8,kl5)  + wst * amat_85        
             ELM_p(mp,ij8,kl6)  =  ELM_p(mp,ij8,kl6)  + wst * amat_86     
             ELM_p(mp,ij8,kl7)  =  ELM_p(mp,ij8,kl7)  + wst * amat_87
             ELM_k(mp,ij8,kl7)  =  ELM_k(mp,ij8,kl7)  + wst * amat_87_k ! Not implemented previously
             ELM_n(mp,ij8,kl7)  =  ELM_n(mp,ij8,kl7)  + wst * amat_87_n   
             ELM_p(mp,ij8,kl8)  =  ELM_p(mp,ij8,kl8)  + wst * amat_88  
             ELM_k(mp,ij8,kl8)  =  ELM_k(mp,ij8,kl8)  + wst * amat_88_k ! Not implemented previously
             ELM_n(mp,ij8,kl8)  =  ELM_n(mp,ij8,kl8)  + wst * amat_88_n      
             ELM_kn(mp,ij8,kl8) =  ELM_kn(mp,ij8,kl8) + wst * amat_88_kn     

             ELM_p(mp,ij9,kl1)  =  ELM_p(mp,ij9,kl1)  + wst * amat_91
             ELM_k(mp,ij9,kl1)  =  ELM_k(mp,ij9,kl1)  + wst * amat_91_k

             ELM_p(mp,ij9,kl2)  =  ELM_p(mp,ij9,kl2)  + wst * amat_92
             ELM_p(mp,ij9,kl3)  =  ELM_p(mp,ij9,kl3)  + wst * amat_93
             ELM_p(mp,ij9,kl5)  =  ELM_p(mp,ij9,kl5)  + wst * amat_95
             ELM_n(mp,ij9,kl5)  =  ELM_n(mp,ij9,kl5)  + wst * amat_95_n
             ELM_k(mp,ij9,kl5)  =  ELM_k(mp,ij9,kl5)  + wst * amat_95_k
             ELM_kn(mp,ij9,kl5) =  ELM_kn(mp,ij9,kl5) + wst * amat_95_kn

             ELM_p(mp,ij9,kl6)  =  ELM_p(mp,ij9,kl6)  + wst * amat_96

             ELM_p(mp,ij9,kl7)  =  ELM_p(mp,ij9,kl7)  + wst * amat_97
             ELM_k(mp,ij9,kl7)  =  ELM_k(mp,ij9,kl7)  + wst * amat_97_k
             ELM_n(mp,ij9,kl7)  =  ELM_n(mp,ij9,kl7)  + wst * amat_97_n

             ELM_p(mp,ij9,kl8)  =  ELM_p(mp,ij9,kl8)  + wst * amat_98
             ELM_n(mp,ij9,kl8)  =  ELM_n(mp,ij9,kl8)  + wst * amat_98_n ! Not implemented previously
             ELM_k(mp,ij9,kl8)  =  ELM_k(mp,ij9,kl8)  + wst * amat_98_k ! Not implemented previously
             ELM_kn(mp,ij9,kl8) =  ELM_kn(mp,ij9,kl8) + wst * amat_98_kn ! Not implemented previously

             ELM_p(mp,ij9,kl9)  =  ELM_p(mp,ij9,kl9)  + wst * amat_99
             ELM_k(mp,ij9,kl9)  =  ELM_k(mp,ij9,kl9)  + wst * amat_99_k
             ELM_n(mp,ij9,kl9)  =  ELM_n(mp,ij9,kl9)  + wst * amat_99_n
             ELM_kn(mp,ij9,kl9) =  ELM_kn(mp,ij9,kl9) + wst * amat_99_kn

           enddo

         enddo
       enddo

     enddo
   enddo

 enddo
enddo

do i=1,n_vertex_max*n_var*n_degrees

  do j=1, n_vertex_max*n_var*n_degrees

    in_fft =  ELM_p(1:n_plane,i,j)
#ifdef USE_FFTW
    call dfftw_execute_dft_r2c(fftw_plan, in_fft, out_fft)
#else
    call my_fft(in_fft, out_fft, n_plane)
#endif

    do k=1,(n_tor+1)/2

      index_k = n_tor*(i-1) + max(2*(k-1),1)

      do m=1,(n_tor+1)/2

        index_m = n_tor*(j-1) + max(2*(m-1),1)

        l = (k-1) + (m-1)

        if ( (l .ge. 0) .and. (l .le. n_plane/2) ) then

          ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   + real(out_fft(l+1))
          ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   - imag(out_fft(l+1))
          ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) - imag(out_fft(l+1))
          ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) - real(out_fft(l+1))

        elseif ( (l .lt. 0) .and. (abs(l) .le. n_plane/2) ) then

          ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   + real(out_fft(abs(l)+1))
          ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   + imag(out_fft(abs(l)+1))
          ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) + imag(out_fft(abs(l)+1))
          ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) - real(out_fft(abs(l)+1))

        endif

        l = (k-1) - (m-1)

        if ( (l .ge. 0) .and. (l .le. n_plane/2) ) then

          ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   + real(out_fft(l+1))
          ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   - imag(out_fft(l+1))
          ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) + imag(out_fft(l+1))
          ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) + real(out_fft(l+1))

        elseif ( (l .lt. 0) .and. (abs(l) .le. n_plane/2) ) then

          ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   + real(out_fft(abs(l)+1))
          ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   + imag(out_fft(abs(l)+1))
          ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) - imag(out_fft(abs(l)+1))
          ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) + real(out_fft(abs(l)+1))

        endif

      enddo

    enddo

  enddo

enddo

do i=1,n_vertex_max*n_var*n_degrees

  do j=1, n_vertex_max*n_var*n_degrees

  if (maxval(abs(ELM_n(1:n_plane,i,j))) .ne. 0.d0) then

    in_fft =  ELM_n(1:n_plane,i,j)
#ifdef USE_FFTW
    call dfftw_execute_dft_r2c(fftw_plan, in_fft, out_fft)
#else
    call my_fft(in_fft, out_fft, n_plane)
#endif

    
    do k=1,(n_tor+1)/2

      index_k = n_tor*(i-1) + max(2*(k-1),1)

      do m=1,(n_tor+1)/2

        im = max(2*(m-1),1)
        index_m = n_tor*(j-1) + max(2*(m-1),1)

        l = (k-1) + (m-1)

        if ( (l .ge. 0) .and. (l .le. n_plane/2) ) then

          ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   + imag(out_fft(l+1)) * float(mode(im))
          ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   + real(out_fft(l+1)) * float(mode(im))
          ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) + real(out_fft(l+1)) * float(mode(im))
          ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) - imag(out_fft(l+1)) * float(mode(im))

        elseif ( (l .lt. 0) .and. (abs(l) .le. n_plane/2) ) then

          ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   - imag(out_fft(abs(l)+1)) * float(mode(im))
          ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   + real(out_fft(abs(l)+1)) * float(mode(im))
          ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) + real(out_fft(abs(l)+1)) * float(mode(im))
          ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) + imag(out_fft(abs(l)+1)) * float(mode(im))

        endif

        l = (k-1) - (m-1)

        if ( (l .ge. 0) .and. (l .le. n_plane/2) ) then

          ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   - imag(out_fft(l+1)) * float(mode(im))
          ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   - real(out_fft(l+1)) * float(mode(im))
          ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) + real(out_fft(l+1)) * float(mode(im))
          ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) - imag(out_fft(l+1)) * float(mode(im))

        elseif ( (l .lt. 0) .and. (abs(l) .le. n_plane/2) ) then

          ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   + imag(out_fft(abs(l)+1)) * float(mode(im))
          ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   - real(out_fft(abs(l)+1)) * float(mode(im))
          ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) + real(out_fft(abs(l)+1)) * float(mode(im))
          ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) + imag(out_fft(abs(l)+1)) * float(mode(im))

        endif

      enddo

    enddo

  endif
  enddo

enddo

do i=1,n_vertex_max*n_var*n_degrees

  do j=1, n_vertex_max*n_var*n_degrees

  if (maxval(abs(ELM_k(1:n_plane,i,j))) .ne. 0.d0) then

    in_fft =  ELM_k(1:n_plane,i,j)
#ifdef USE_FFTW
    call dfftw_execute_dft_r2c(fftw_plan, in_fft, out_fft)
#else
    call my_fft(in_fft, out_fft, n_plane)
#endif
    
    do k=1,(n_tor+1)/2

      ik      = max(2*(k-1),1)
      index_k = n_tor*(i-1) + max(2*(k-1),1)

      do m=1,(n_tor+1)/2

        index_m = n_tor*(j-1) + max(2*(m-1),1)

        l = (k-1) + (m-1)

        if ( (l .ge. 0) .and. (l .le. n_plane/2) ) then

          ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   + imag(out_fft(l+1)) * float(mode(ik))
          ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   + real(out_fft(l+1)) * float(mode(ik))
          ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) + real(out_fft(l+1)) * float(mode(ik))
          ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) - imag(out_fft(l+1)) * float(mode(ik))

        elseif ( (l .lt. 0) .and. (abs(l) .le. n_plane/2) ) then

          ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   - imag(out_fft(abs(l)+1)) * float(mode(ik))
          ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   + real(out_fft(abs(l)+1)) * float(mode(ik))
          ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) + real(out_fft(abs(l)+1)) * float(mode(ik))
          ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) + imag(out_fft(abs(l)+1)) * float(mode(ik))

        endif

        l = (k-1) - (m-1)

        if ( (l .ge. 0) .and. (l .le. n_plane/2) ) then

          ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   + imag(out_fft(l+1)) * float(mode(ik))
          ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   + real(out_fft(l+1)) * float(mode(ik))
          ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) - real(out_fft(l+1)) * float(mode(ik))
          ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) + imag(out_fft(l+1)) * float(mode(ik))

        elseif ( (l .lt. 0) .and. (abs(l) .le. n_plane/2) ) then

          ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   - imag(out_fft(abs(l)+1)) * float(mode(ik))
          ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   + real(out_fft(abs(l)+1)) * float(mode(ik))
          ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) - real(out_fft(abs(l)+1)) * float(mode(ik))
          ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) - imag(out_fft(abs(l)+1)) * float(mode(ik))

        endif

      enddo

    enddo

  endif
  enddo

enddo

do i=1,n_vertex_max*n_var*n_degrees

  do j=1, n_vertex_max*n_var*n_degrees

  if (maxval(abs(ELM_kn(1:n_plane,i,j))) .ne. 0.d0) then

    in_fft =  ELM_kn(1:n_plane,i,j)
#ifdef USE_FFTW
    call dfftw_execute_dft_r2c(fftw_plan, in_fft, out_fft)
#else
    call my_fft(in_fft, out_fft, n_plane)
#endif
    
    do k=1,(n_tor+1)/2

      ik      = max(2*(k-1),1)
      index_k = n_tor*(i-1) + max(2*(k-1),1)

      do m=1,(n_tor+1)/2

        im      = max(2*(m-1),1)
        index_m = n_tor*(j-1) + max(2*(m-1),1)

        l = (k-1) + (m-1)

        if ( (l .ge. 0) .and. (l .le. n_plane/2) ) then

           ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   - real(out_fft(l+1)) * float(mode(im)) * float(mode(ik))
           ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   + imag(out_fft(l+1)) * float(mode(im)) * float(mode(ik))
           ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) + imag(out_fft(l+1)) * float(mode(im)) * float(mode(ik))
           ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) + real(out_fft(l+1)) * float(mode(im)) * float(mode(ik))

        elseif ( (l .lt. 0) .and. (abs(l) .le. n_plane/2) ) then

           ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   - real(out_fft(abs(l)+1)) * float(mode(im)) * float(mode(ik))
           ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   - imag(out_fft(abs(l)+1)) * float(mode(im)) * float(mode(ik))
           ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   - imag(out_fft(abs(l)+1)) * float(mode(im)) * float(mode(ik))
           ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) - imag(out_fft(abs(l)+1)) * float(mode(im)) * float(mode(ik))
           ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) + real(out_fft(abs(l)+1)) * float(mode(im)) * float(mode(ik))

        endif

        l = (k-1) - (m-1)

        if ( (l .ge. 0) .and. (l .le. n_plane/2) ) then

          ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   + real(out_fft(l+1)) * float(mode(im)) * float(mode(ik))
          ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   - imag(out_fft(l+1)) * float(mode(im)) * float(mode(ik))
          ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) + imag(out_fft(l+1)) * float(mode(im)) * float(mode(ik))
          ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) + real(out_fft(l+1)) * float(mode(im)) * float(mode(ik))

        elseif ( (l .lt. 0) .and. (abs(l) .le. n_plane/2) ) then

          ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   + real(out_fft(abs(l)+1)) * float(mode(im)) * float(mode(ik))
          ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   + imag(out_fft(abs(l)+1)) * float(mode(im)) * float(mode(ik))
          ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) - imag(out_fft(abs(l)+1)) * float(mode(im)) * float(mode(ik))
          ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) + real(out_fft(abs(l)+1)) * float(mode(im)) * float(mode(ik))

        endif

      enddo

    enddo

  endif
  enddo

enddo

ELM = 0.5d0 * ELM

do j=1, n_vertex_max*n_var*n_degrees

  in_fft = RHS_p(1:n_plane,j)
#ifdef USE_FFTW
  call dfftw_execute_dft_r2c(fftw_plan, in_fft, out_fft)
#else
  call my_fft(in_fft, out_fft, n_plane)
#endif
    
  index = n_tor*(j-1) + 1

  RHS(index) = real(out_fft(1))

  do k=2,(n_tor+1)/2

    index = n_tor*(j-1) + 2*(k-1)

    RHS(index)   =   real(out_fft(k))
    RHS(index+1) = - imag(out_fft(k))

  enddo

enddo

do j=1, n_vertex_max*n_var*n_degrees

  in_fft = RHS_k(1:n_plane,j)
#ifdef USE_FFTW
  call dfftw_execute_dft_r2c(fftw_plan, in_fft, out_fft)
#else
  call my_fft(in_fft, out_fft, n_plane)
#endif
  
  index = n_tor*(j-1) + 1
  ik    = 1

  RHS(index) = RHS(index) + imag(out_fft(1)) * float(mode(ik))

  do k=2,(n_tor+1)/2

    ik    = max(2*(k-1),1)
    index = n_tor*(j-1) + 2*(k-1)

    RHS(index)   = RHS(index)   + imag(out_fft(k)) * float(mode(ik))
    RHS(index+1) = RHS(index+1) + real(out_fft(k)) * float(mode(ik))

  enddo

enddo

return
end subroutine element_matrix_fft

subroutine my_fft(in_fft,out_fft,n)

real*8     :: in_fft(*)
complex*16 :: out_fft(*)

integer    :: i, n
real*8     :: tmp_fft(2*n+2)

tmp_fft(1:n) = in_fft(1:n)

call RFT2(tmp_fft,n,1)

do i=1,n
  out_fft(i) = cmplx(tmp_fft(2*i-1),tmp_fft(2*i))
enddo

return
end subroutine my_fft
end module mod_elt_matrix_fft
