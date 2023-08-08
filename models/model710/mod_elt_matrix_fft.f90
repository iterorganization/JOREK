module mod_elt_matrix_fft
  implicit none
contains

#include "corr_neg_include.f90"

subroutine element_matrix_fft(element, nodes, xpoint2, xcase2, R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint, Z_xpoint, &
                              ELM, RHS, tid, ELM_p, ELM_n, ELM_k, ELM_kn, RHS_p, RHS_k,                               &
                              eq_g, eq_s, eq_t, eq_p, eq_ss, eq_st, eq_tt, delta_g, delta_s, delta_t,                 & 
                              i_tor_min, i_tor_max, aux_nodes, ELM_pnn)

use mod_parameters
use data_structure
use gauss
use basis_at_gaussian
use phys_module
use diffusivities, only: get_dperp, get_zkperp
use equil_info, only : get_psi_n, ES
use mod_F_profile
use mod_bootstrap_functions
use mod_sources

implicit none

! --- Input Variables
type (type_element)        :: element
type (type_node)           :: nodes(n_vertex_max)     ! fluid variables
type (type_node), optional :: aux_nodes(n_vertex_max) ! particle moments

logical, intent(in)    :: xpoint2
integer, intent(in)    :: xcase2
real*8,  intent(in)    :: R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint(2), Z_xpoint(2)

#define DIM0 n_tor*n_vertex_max*n_degrees*n_var
#define DIM1 n_plane
#define DIM2 1:n_vertex_max*n_var*n_degrees

real*8, dimension (DIM0,DIM0)  :: ELM
real*8, dimension (DIM0)       :: RHS
integer, intent(in)            :: tid, i_tor_min, i_tor_max
real*8, dimension(DIM1, DIM2, DIM2), intent(inout) :: ELM_p
real*8, dimension(DIM1, DIM2, DIM2), intent(inout) :: ELM_n
real*8, dimension(DIM1, DIM2, DIM2), intent(inout) :: ELM_k
real*8, dimension(DIM1, DIM2, DIM2), intent(inout) :: ELM_kn
real*8, dimension(DIM1, DIM2),       intent(inout) :: RHS_p
real*8, dimension(DIM1, DIM2),       intent(inout) :: RHS_k

real*8, dimension(n_plane,n_var,n_gauss,n_gauss), intent(inout) :: eq_g, eq_s, eq_t
real*8, dimension(n_plane,n_var,n_gauss,n_gauss), intent(inout) :: eq_p
real*8, dimension(n_plane,n_var,n_gauss,n_gauss), intent(inout) :: eq_ss, eq_st, eq_tt
real*8, dimension(n_plane,n_var,n_gauss,n_gauss), intent(inout) :: delta_g, delta_s, delta_t
real*8, dimension(n_plane,n_var,n_gauss,n_gauss) :: eq_aux_g, eq_aux_s, eq_aux_t, eq_aux_p

! --- Variables outside the OMP loop
integer    :: n_tor_start, n_tor_end, n_tor_local
integer    :: i, j, index, index_k, index_m, i_v, j_loc, i_loc, m, ik
integer    :: in, im, ivar, kvar, ms, mt, mp
real*8     :: wst, xjac, xjac_R, xjac_Z, R, Z, theta, zeta
real*8     :: current_source_JR(n_gauss,n_gauss), current_source_JZ(n_gauss,n_gauss), current_source_Jp(n_gauss,n_gauss)
real*8     :: particle_source(n_gauss,n_gauss),heat_source(n_gauss,n_gauss)
real*8     :: in_fft(1:n_plane)
complex*16 :: out_fft(1:n_plane)
integer    :: VmsType=0, ViscType=0
real*8     :: TG_NUM_Eq, CoefAdv=0.0
real*8     :: Coef_DivV
real*8     :: psi_axisym(n_gauss,n_gauss), psi_axisym_s(n_gauss,n_gauss), psi_axisym_t(n_gauss,n_gauss)
real*8     ::                              psi_axisym_R(n_gauss,n_gauss), psi_axisym_Z(n_gauss,n_gauss)
real*8     :: Fprof_time_dep,dF_dpsi(n_gauss,n_gauss)      ,dF_dz      ,dF_dpsi2      ,dF_dz2      ,dF_dpsi_dz
real*8     :: zFFprime      ,dFFprime_dpsi,dFFprime_dz,dFFprime_dpsi2,dFFprime_dz2,dFFprime_dpsi_dz
real*8     :: rho_initial(n_gauss,n_gauss),dn_dpsi,dn_dz,dn_dpsi2,dn_dz2,dn_dpsi_dz,dn_dpsi3,dn_dpsi_dz2, dn_dpsi2_dz
real*8     :: T_initial  (n_gauss,n_gauss),dT_dpsi,dT_dz,dT_dpsi2,dT_dz2,dT_dpsi_dz,dT_dpsi3,dT_dpsi_dz2, dT_dpsi2_dz
real*8     :: Jb, Jb_0
real*8     :: Ti0, Ti0_R, Ti0_Z
real*8     :: Te0, Te0_R, Te0_Z
real*8     :: amu_neo_prof(n_gauss,n_gauss), aki_neo_prof(n_gauss,n_gauss)
real*8     :: V_source(n_gauss,n_gauss), dV_dpsi_source(n_gauss,n_gauss),dV_dz_source(n_gauss,n_gauss)
real*8     :: dV_dpsi2,dV_dz2,dV_dpsi_dz,dV_dpsi3,dV_dpsi_dz2, dV_dpsi2_dz
real*8     :: eta_ARAZ, tauIC_ARAZ
logical    :: use_fft


real*8, dimension(n_gauss,n_gauss)    :: x_g, x_s, x_t, x_ss, x_st, x_tt
real*8, dimension(n_gauss,n_gauss)    :: y_g, y_s, y_t, y_ss, y_st, y_tt
real*8, dimension(n_gauss,n_gauss)    :: Fprofile
real*8, dimension(n_var)              :: TG_NUM
real*8, dimension(n_tor,n_plane) :: HHZ, HHZ_p, HHZ_pp

! --- Variables inside the OMP loop
integer    :: k, l, index_ij, index_kl, ij, kl

real*8     :: AR0,  AR0_R,  AR0_Z,  AR0_p,  AR0_s,  AR0_t, AR0_ss, AR0_tt, AR0_st, AR0_RR, AR0_ZZ, AR0_RZ, AR0_pp
real*8     :: AZ0,  AZ0_R,  AZ0_Z,  AZ0_p,  AZ0_s,  AZ0_t, AZ0_ss, AZ0_tt, AZ0_st, AZ0_RR, AZ0_ZZ, AZ0_RZ, AZ0_pp
real*8     :: A30,  A30_R,  A30_Z,  A30_p,  A30_s,  A30_t, A30_ss, A30_tt, A30_st, A30_RR, A30_ZZ, A30_RZ, A30_pp
real*8     :: UR0,  UR0_R,  UR0_Z,  UR0_p,  UR0_s,  UR0_t, UR0_ss, UR0_st, UR0_tt, UR0_RR, UR0_ZZ, UR0_RZ, UR0_pp
real*8     :: UZ0,  UZ0_R,  UZ0_Z,  UZ0_p,  UZ0_s,  UZ0_t, UZ0_ss, UZ0_st, UZ0_tt, UZ0_RR, UZ0_ZZ, UZ0_RZ, UZ0_pp
real*8     :: Up0,  Up0_R,  Up0_Z,  Up0_p,  Up0_s,  Up0_t, Up0_ss, Up0_st, Up0_tt, Up0_RR, Up0_ZZ, Up0_RZ, Up0_pp
real*8     :: rho0, rho0_R, rho0_Z, rho0_p, rho0_s, rho0_t,rho0_ss,rho0_st,rho0_tt,rho0_RR,rho0_ZZ,rho0_RZ,rho0_pp, rho0_corr
real*8     :: T0,   T0_R,   T0_Z,   T0_p,   T0_s,   T0_t,  T0_ss,  T0_st,  T0_tt,  T0_RR,  T0_ZZ,  T0_RZ,  T0_pp,   T0_corr
real*8     :: p0,   p0_R,   p0_Z,   p0_p,   p0_s,   p0_t,  p0_corr

real*8     :: AR,  AR_R,  AR_Z,  AR_p,  AR_s,  AR_t
real*8     :: AZ,  AZ_R,  AZ_Z,  AZ_p,  AZ_s,  AZ_t
real*8     :: A3,  A3_R,  A3_Z,  A3_p,  A3_s,  A3_t
real*8     :: UR,  UR_R,  UR_Z,  UR_p,  UR_s,  UR_t
real*8     :: UZ,  UZ_R,  UZ_Z,  UZ_p,  UZ_s,  UZ_t
real*8     :: Up,  Up_R,  Up_Z,  Up_p,  Up_s,  Up_t
real*8     :: T,   T_R,   T_Z,   T_p,   T_s,   T_t
real*8     :: rho, rho_R, rho_Z, rho_p, rho_s, rho_t

real*8     :: v,  v_R,  v_Z,  v_s,  v_t,  v_p,  v_ss,  v_st,  v_tt,  v_RR,  v_ZZ
real*8     :: bf, bf_R, bf_Z, bf_s, bf_t, bf_p, bf_ss, bf_st, bf_tt, bf_RR, bf_ZZ

real*8     :: Fprof
real*8     :: BR0, BR0_AR,    BR0_AZ__n, BR0_A3
real*8     :: BZ0, BZ0_AR__n, BZ0_AZ,    BZ0_A3
real*8     :: Bp0, Bp0_AR,    Bp0_AZ,    Bp0_A3
real*8     :: BB2, BB2_AR__p, BB2_AR__n, BB2_AZ__p, BB2_AZ__n, BB2_A3
real*8     :: Bp00

real*8     :: BgradT, BgradT_AR__p, BgradT_AR__n, BgradT_AZ__p, BgradT_AZ__n, BgradT_A3, BgradT_T__p, BgradT_T__n

real*8     :: BgradRho, BgradRho_AR__p, BgradRho_AR__n, BgradRho_AZ__p, BgradRho_AZ__n, BgradRho_A3, BgradRho_rho__p, BgradRho_rho__n

real*8     :: BgradP, BgradP_AR__p, BgradP_AR__n, BgradP_AZ__p, BgradP_AZ__n, BgradP_A3, BgradP_T__p, BgradP_T__n, BgradP_rho__p, BgradP_rho__n

real*8     :: BgradVstar__p, BgradVstar__k
real*8     :: BgradVstar_AR__p, BgradVstar_AR__k, BgradVstar_AR__n
real*8     :: BgradVstar_AZ__p, BgradVstar_AZ__k, BgradVstar_AZ__n
real*8     :: BgradVstar_A3__p, BgradVstar_A3__k

real*8     :: UgradRho, UgradRho_UR, UgradRho_UZ, UgradRho_Up, UgradRho_rho__p, UgradRho_rho__n
real*8     :: UgradT,   UgradT_UR,   UgradT_UZ,   UgradT_Up,   UgradT_T__p,     UgradT_T__n

real*8     :: UgradVstar__p, UgradVstar__k, UgradVstar_UR, UgradVstar_UZ, UgradVstar_Up__k

real*8     :: gradRho_gradVstar__p, gradRho_gradVstar__k, gradRho_gradVstar_rho__p, gradRho_gradVstar_rho__kn
real*8     :: gradT_gradVstar__p,   gradT_gradVstar__k,   gradT_gradVstar_T__p,   gradT_gradVstar_T__kn

real*8     :: gradBF_gradVstar__p, gradBF_gradVstar__kn, UgradBF__p, UgradBF__n, BgradBF__p, BgradBF__n

real*8     :: UgradUR, UgradUR_UR__p, UgradUR_UR__n, UgradUR_UZ,    UgradUR_Up
real*8     :: UgradUZ, UgradUZ_UR,    UgradUZ_UZ__p, UgradUZ_UZ__n, UgradUZ_Up
real*8     :: UgradUp, UgradUp_UR,    UgradUp_UZ   , UgradUp_Up__p, UgradUp_Up__n

real*8     :: divU, divU_UR, divU_UZ, divU_Up__n

real*8     :: divRhoU, divRhoU_UR, divRhoU_UZ, divRhoU_Up__p, divRhoU_Up__n, divRhoU_rho__p, divRhoU_rho__n

real*8     :: Vt0, Vt0_R, Vt0_Z

real*8     :: tau_IC, distance_bnd

real*8     :: VdiaR0
real*8     :: VdiaR0_AR__p,  VdiaR0_AR__n
real*8     :: VdiaR0_AZ__p,  VdiaR0_AZ__n
real*8     :: VdiaR0_A3__p,  VdiaR0_A3__n
real*8     :: VdiaR0_rho__p, VdiaR0_rho__n
real*8     :: VdiaR0_T__p,   VdiaR0_T__n

real*8     :: VdiaZ0
real*8     :: VdiaZ0_AR__p,  VdiaZ0_AR__n
real*8     :: VdiaZ0_AZ__p,  VdiaZ0_AZ__n
real*8     :: VdiaZ0_A3__p,  VdiaZ0_A3__n
real*8     :: VdiaZ0_rho__p, VdiaZ0_rho__n
real*8     :: VdiaZ0_T__p,   VdiaZ0_T__n

real*8     :: VdiaP0
real*8     :: VdiaP0_AR__p,  VdiaP0_AR__n
real*8     :: VdiaP0_AZ__p,  VdiaP0_AZ__n
real*8     :: VdiaP0_A3__p,  VdiaP0_A3__n
real*8     :: VdiaP0_rho__p, VdiaP0_rho__n
real*8     :: VdiaP0_T__p,   VdiaP0_T__n

real*8     :: VdiaGradUR
real*8     :: VdiaGradUR_AR__p,  VdiaGradUR_AR__n
real*8     :: VdiaGradUR_AZ__p,  VdiaGradUR_AZ__n
real*8     :: VdiaGradUR_A3__p,  VdiaGradUR_A3__n
real*8     :: VdiaGradUR_rho__p, VdiaGradUR_rho__n
real*8     :: VdiaGradUR_T__p,   VdiaGradUR_T__n
real*8     :: VdiaGradUR_UR__p,  VdiaGradUR_UR__n

real*8     :: VdiaGradUZ
real*8     :: VdiaGradUZ_AR__p,  VdiaGradUZ_AR__n
real*8     :: VdiaGradUZ_AZ__p,  VdiaGradUZ_AZ__n
real*8     :: VdiaGradUZ_A3__p,  VdiaGradUZ_A3__n
real*8     :: VdiaGradUZ_rho__p, VdiaGradUZ_rho__n
real*8     :: VdiaGradUZ_T__p,   VdiaGradUZ_T__n
real*8     :: VdiaGradUZ_UZ__p,  VdiaGradUZ_UZ__n

real*8     :: VdiaGradUp
real*8     :: VdiaGradUp_AR__p,  VdiaGradUp_AR__n
real*8     :: VdiaGradUp_AZ__p,  VdiaGradUp_AZ__n
real*8     :: VdiaGradUp_A3__p,  VdiaGradUp_A3__n
real*8     :: VdiaGradUp_rho__p, VdiaGradUp_rho__n
real*8     :: VdiaGradUp_T__p,   VdiaGradUp_T__n
real*8     :: VdiaGradUp_Up__p,  VdiaGradUp_Up__n

real*8     :: VdiaGradVstar__p, VdiaGradVstar__k
real*8     :: VdiaGradVstar_AR__p,  VdiaGradVstar_AR__n,  VdiaGradVstar_AR__k,  VdiaGradVstar_AR__kn
real*8     :: VdiaGradVstar_AZ__p,  VdiaGradVstar_AZ__n,  VdiaGradVstar_AZ__k,  VdiaGradVstar_AZ__kn
real*8     :: VdiaGradVstar_A3__p,  VdiaGradVstar_A3__n,  VdiaGradVstar_A3__k,  VdiaGradVstar_A3__kn
real*8     :: VdiaGradVstar_rho__p, VdiaGradVstar_rho__n, VdiaGradVstar_rho__k, VdiaGradVstar_rho__kn
real*8     :: VdiaGradVstar_T__p,   VdiaGradVstar_T__n,   VdiaGradVstar_T__k,   VdiaGradVstar_T__kn

real*8     :: Btht, Btht_AR__p, Btht_AR__n, Btht_AZ__p, Btht_AZ__n, Btht_A3
real*8     :: Vtht
real*8     :: Vtht_AR__p,  Vtht_AR__n
real*8     :: Vtht_AZ__p,  Vtht_AZ__n
real*8     :: Vtht_A3__p,  Vtht_A3__n
real*8     :: Vtht_rho__p, Vtht_rho__n
real*8     :: Vtht_T__p,   Vtht_T__n
real*8     :: Vtht_UR, Vtht_UZ
real*8     :: Vneo
real*8     :: Vneo_AR__p,  Vneo_AR__n
real*8     :: Vneo_AZ__p,  Vneo_AZ__n
real*8     :: Vneo_A3__p,  Vneo_A3__n
real*8     :: Vneo_T__p,   Vneo_T__n
real*8     :: PneoR
real*8     :: PneoR_AR__p,  PneoR_AR__n
real*8     :: PneoR_AZ__p,  PneoR_AZ__n
real*8     :: PneoR_A3__p,  PneoR_A3__n
real*8     :: PneoR_rho__p, PneoR_rho__n
real*8     :: PneoR_T__p,   PneoR_T__n
real*8     :: PneoR_UR,     PneoR_UZ
real*8     :: PneoZ
real*8     :: PneoZ_AR__p,  PneoZ_AR__n
real*8     :: PneoZ_AZ__p,  PneoZ_AZ__n
real*8     :: PneoZ_A3__p,  PneoZ_A3__n
real*8     :: PneoZ_rho__p, PneoZ_rho__n
real*8     :: PneoZ_T__p,   PneoZ_T__n
real*8     :: PneoZ_UR,     PneoZ_UZ

real*8     :: ZK_prof, D_prof, psi_norm

real*8     :: eta_T, visco_T, deta_dT, d2eta_d2T, dvisco_dT, visco_num_T, visco_divV, dvisco_divV_dT
real*8     :: eta_num_T, eta_R, eta_Z, eta_p, Zkpar_T, dZKpar_dt
real*8     :: eta_T_T, eta_R_T, eta_Z_T, eta_p_T__p, eta_p_T__n

real*8     :: Qconv_UR
real*8     :: Qconv_UR_AR__p,  Qconv_UR_AR__n
real*8     :: Qconv_UR_AZ__p,  Qconv_UR_AZ__n
real*8     :: Qconv_UR_A3__p,  Qconv_UR_A3__n
real*8     :: Qconv_UR_UR__p,  Qconv_UR_UR__n
real*8     :: Qconv_UR_UZ__p,  Qconv_UR_UZ__n
real*8     :: Qconv_UR_Up__p,  Qconv_UR_Up__n
real*8     :: Qconv_UR_rho__p, Qconv_UR_rho__n
real*8     :: Qconv_UR_T__p,   Qconv_UR_T__n

real*8     :: Qconv_UZ
real*8     :: Qconv_UZ_AR__p,  Qconv_UZ_AR__n
real*8     :: Qconv_UZ_AZ__p,  Qconv_UZ_AZ__n
real*8     :: Qconv_UZ_A3__p,  Qconv_UZ_A3__n
real*8     :: Qconv_UZ_UR__p,  Qconv_UZ_UR__n
real*8     :: Qconv_UZ_UZ__p,  Qconv_UZ_UZ__n
real*8     :: Qconv_UZ_Up__p,  Qconv_UZ_Up__n
real*8     :: Qconv_UZ_rho__p, Qconv_UZ_rho__n
real*8     :: Qconv_UZ_T__p,   Qconv_UZ_T__n

real*8     :: Qconv_Up
real*8     :: Qconv_Up_AR__p,  Qconv_Up_AR__n
real*8     :: Qconv_Up_AZ__p,  Qconv_Up_AZ__n
real*8     :: Qconv_Up_A3__p,  Qconv_Up_A3__n
real*8     :: Qconv_Up_UR__p,  Qconv_Up_UR__n
real*8     :: Qconv_Up_UZ__p,  Qconv_Up_UZ__n
real*8     :: Qconv_Up_Up__p,  Qconv_Up_Up__n
real*8     :: Qconv_Up_rho__p, Qconv_Up_rho__n
real*8     :: Qconv_Up_T__p,   Qconv_Up_T__n

real*8     :: JxB_UR__p, JxB_UR__k
real*8     :: JxB_UZ__p, JxB_UZ__k
real*8     :: JxB_Up__p, JxB_Up__k
real*8     :: JxB_UR_AR__p, JxB_UR_AR__k, JxB_UR_AR__n, JxB_UR_AR__kn
real*8     :: JxB_UZ_AR__p, JxB_UZ_AR__k, JxB_UZ_AR__n, JxB_UZ_AR__kn
real*8     :: JxB_Up_AR__p, JxB_Up_AR__k, JxB_Up_AR__n, JxB_Up_AR__kn
real*8     :: JxB_UR_AZ__p, JxB_UR_AZ__k, JxB_UR_AZ__n, JxB_UR_AZ__kn
real*8     :: JxB_UZ_AZ__p, JxB_UZ_AZ__k, JxB_UZ_AZ__n, JxB_UZ_AZ__kn
real*8     :: JxB_Up_AZ__p, JxB_Up_AZ__k, JxB_Up_AZ__n, JxB_Up_AZ__kn
real*8     :: JxB_UR_A3__p, JxB_UR_A3__k, JxB_UR_A3__n, JxB_UR_A3__kn
real*8     :: JxB_UZ_A3__p, JxB_UZ_A3__k, JxB_UZ_A3__n, JxB_UZ_A3__kn
real*8     :: JxB_Up_A3__p, JxB_Up_A3__k, JxB_Up_A3__n, JxB_Up_A3__kn

real*8     :: Qvisc_UR__p, Qvisc_UR__k
real*8     :: Qvisc_UR_UR__p, Qvisc_UR_UR__k, Qvisc_UR_UR__n, Qvisc_UR_UR__kn
real*8     :: Qvisc_UR_UZ__p, Qvisc_UR_UZ__k, Qvisc_UR_UZ__n, Qvisc_UR_UZ__kn
real*8     :: Qvisc_UR_Up__p, Qvisc_UR_Up__k, Qvisc_UR_Up__n, Qvisc_UR_Up__kn

real*8     :: Qvisc_UZ__p, Qvisc_UZ__k
real*8     :: Qvisc_UZ_UR__p, Qvisc_UZ_UR__k, Qvisc_UZ_UR__n, Qvisc_UZ_UR__kn
real*8     :: Qvisc_UZ_UZ__p, Qvisc_UZ_UZ__k, Qvisc_UZ_UZ__n, Qvisc_UZ_UZ__kn
real*8     :: Qvisc_UZ_Up__p, Qvisc_UZ_Up__k, Qvisc_UZ_Up__n, Qvisc_UZ_Up__kn

real*8     :: Qvisc_Up__p, Qvisc_Up__k
real*8     :: Qvisc_Up_UR__p, Qvisc_Up_UR__k, Qvisc_Up_UR__n, Qvisc_Up_UR__kn
real*8     :: Qvisc_Up_UZ__p, Qvisc_Up_UZ__k, Qvisc_Up_UZ__n, Qvisc_Up_UZ__kn
real*8     :: Qvisc_Up_Up__p, Qvisc_Up_Up__k, Qvisc_Up_Up__n, Qvisc_Up_Up__kn

real*8     :: QviscT0
real*8     :: QviscT0_UR__p, QviscT0_UR__n
real*8     :: QviscT0_UZ__p, QviscT0_UZ__n
real*8     :: QviscT0_Up__p, QviscT0_Up__n
real*8     :: QviscT0_T__p,  QviscT0_T__n
real*8     :: Qvisc_T
real*8     :: Qvisc_T_UR__p, Qvisc_T_UR__n
real*8     :: Qvisc_T_UZ__p, Qvisc_T_UZ__n
real*8     :: Qvisc_T_Up__p, Qvisc_T_Up__n
real*8     :: Qvisc_T_T__p,  Qvisc_T_T__n

! --- Auxiliary values of pressure coupling
real*8     :: aux_PIRR, aux_PIRR_s, aux_PIRR_t, aux_PIRR_R, aux_PIRR_Z,aux_PIRR_p
real*8     :: aux_PIZZ, aux_PIZZ_s, aux_PIZZ_t, aux_PIZZ_R, aux_PIZZ_Z,aux_PIZZ_p
real*8     :: aux_PIPP, aux_PIPP_s, aux_PIPP_t, aux_PIPP_R, aux_PIPP_Z,aux_PIPP_p
real*8     :: aux_PIZR, aux_PIZR_s, aux_PIZR_t, aux_PIZR_R, aux_PIZR_Z,aux_PIZR_p
real*8     :: aux_PIRP, aux_PIRP_s, aux_PIRP_t, aux_PIRP_R, aux_PIRP_Z,aux_PIRP_p
real*8     :: aux_PIZP, aux_PIZP_s, aux_PIZP_t, aux_PIZP_R, aux_PIZP_Z,aux_PIZP_p
real*8     :: aux_divPIR_perp, aux_divPIZ_perp,     aux_divPIp_perp
real*8     :: aux_divPIR,      aux_divPIZ,          aux_divPIp
real*8     :: aux_BdivPI



! --- fourth order stabilization (numerical diffusion terms)
real*8     :: lap_Vstar, lap_bf
real*8     :: lap_AR, lap_AZ, lap_A3
real*8     :: lap_UR, lap_UZ, lap_Up
real*8     :: lap_rho, lap_T

! --- VMS
real*8     :: VdotB
real*8     :: CvR0, CvZ0, Cvp0, CvGradAR0, CvGradAZ0, CvGradA30, CvGradr0, CvGradT0
real*8     :: VbR0, VbZ0, Vbp0, VbGradAR0, VbGradAZ0, VbGradA30, VbGradr0, VbGradT0
real*8     :: CvGradUR0, CvGradUZ0, CvGradUp0, VbGradUR0, VbGradUZ0, VbGradUp0

real*8     :: CvGradVi__p, VbGradVi__p, CvGradVi__k, VbGradVi__k, DiveRMVi, DiveZMVi, DivePMVi__k
real*8     :: CvGradVj__p, CvGradVj__n, VbGradVj__p, VbGradVj__n, DiveRMVj, DiveZMVj, DivePMVj__n

real*8     :: VmsCoefF, VmsCoefF_T

! --- Matrix
real*8, dimension(n_var,n_var)   :: QvmsAd_p, QvmsAd_n, QvmsAd_k, QvmsAd_kn, QvmsF_p, QvmsF_n, QvmsF_k, QvmsF_kn
real*8, dimension(n_var      )   :: rhs_p_ij, rhs_k_ij, Pvec_prev, Qvec_p, Qvec_k, VMS__p, VMS__k
real*8, dimension(n_var,n_var)   :: amat, Pjac, Qjac_p, Qjac_k, Qjac_n, Qjac_kn

! --- Ohmic heating, for details please see:
! https://www.jorek.eu/wiki/doku.php?id=ohmic_heating
real*8, dimension(DIM1, DIM2, DIM2), intent(inout) :: ELM_pnn
real*8, dimension(n_plane,n_var,n_gauss,n_gauss) :: eq_pp, eq_sp, eq_tp
real*8, dimension(n_var,n_var)   :: Qjac_pnn
real*8     :: AR0_Rp, AR0_Zp, AR0_sp, AR0_tp
real*8     :: AZ0_Rp, AZ0_Zp, AZ0_sp, AZ0_tp
real*8     :: A30_Rp, A30_Zp, A30_sp, A30_tp
real*8     :: bf_RZ, bf_pp, bf_sp, bf_tp, bf_Rp, bf_Zp
real*8     :: AR_RR, AR_RZ, AR_ZZ, AR_pp, AR_Rp, AR_Zp
real*8     :: AZ_RR, AZ_RZ, AZ_ZZ, AZ_pp, AZ_Rp, AZ_Zp
real*8     :: A3_RR, A3_RZ, A3_ZZ, A3_pp, A3_Rp, A3_Zp
real*8     :: eta_T_ohm, deta_dT_ohm
real*8     :: psieq_R,  psieq_Z
real*8     :: BR0_R, BR0_Z, BR0_p, BZ0_R, BZ0_Z, BZ0_p, Bp0_R, Bp0_Z, Bp0_p
real*8     :: JR0 , JR0_AR__p, JR0_AR__n, JR0_AR__nn, JR0_AZ__p, JR0_AZ__n, JR0_AZ__nn, JR0_A3__p, JR0_A3__n, JR0_A3__nn
real*8     :: JZ0 , JZ0_AR__p, JZ0_AR__n, JZ0_AR__nn, JZ0_AZ__p, JZ0_AZ__n, JZ0_AZ__nn, JZ0_A3__p, JZ0_A3__n, JZ0_A3__nn
real*8     :: Jp0 , Jp0_AR__p, Jp0_AR__n, Jp0_AR__nn, Jp0_AZ__p, Jp0_AZ__n, Jp0_AZ__nn, Jp0_A3__p, Jp0_A3__n, Jp0_A3__nn
real*8     :: JJ2,  JJ2_AR__p, JJ2_AR__n, JJ2_AR__nn, JJ2_AZ__p, JJ2_AZ__n, JJ2_AZ__nn, JJ2_A3__p, JJ2_A3__n, JJ2_A3__nn

rho_min = 0.005 ! should be moved to namelist input

! --- Main switches
Coef_DivV = 0.0d0 ! this is a stabilisation term !
ViscType  = 22
VmsType   = -1    ! this is a stabilisation term if >=0 !
CoefAdv   = 0.d0  ! this is a stabilisation term (part of VMS) !

! --- Switches for numerical stability of resistive and diamagnetic terms in AR and AZ equations
eta_ARAZ  = 0.d0  ! =0.0 to switch off resistive   terms for AR and AZ equations
tauIC_ARAZ= 0.d0  ! =0.0 to switch off diamagnetic terms for AR and AZ equations
if (eta_ARAZ_on  ) eta_ARAZ   = 1.d0 ! switched on by default
if (tauIC_ARAZ_on) tauIC_ARAZ = 1.d0 ! switched on by default
if (.not. present(aux_nodes)) then
  do i=1,n_vertex_max
    aux_nodes(i)%values(:,:,:) = 0.d0
  enddo
endif

! --- Info about ViscType:
! --- =0  : Full-MHD resistivity
! --- =10 : Generic form of the full viscous tensor (from B.Nkonga and A.Bhole)
! --- =15 : Mock-up of reduced-MHD resistivity (from B.Nkonga)
! --- =20 : Viscous terms as implemented by W.Haverkort
! --- =22 : The correct viscous tensor (there is only one.)


! --- The settings that will reproduce the old 710 set-up (assuming you are also
! --- changing the F-profile to be wrong, which I have not done here obviously)
!Coef_DivV           = 0.0d0 ! =0 means visco_divV=visco_T, otherwise visco_divV=visco_T+Coef_DivV
!ViscType            = 0
!VmsType             = 0
!CoefAdv             = 0.0


! --- Initialise
ELM_p  = 0.d0
ELM_n  = 0.d0
ELM_k  = 0.d0
ELM_kn = 0.d0
RHS_p  = 0.d0
RHS_k  = 0.d0
ELM    = 0.d0
RHS    = 0.d0
ELM_pnn=0.d0

rhs_p_ij  = 0.d0
rhs_k_ij  = 0.d0
amat      = 0.d0
Pjac      = 0.d0
Qjac_p    = 0.d0
Qjac_k    = 0.d0
Qjac_n    = 0.d0
Qjac_kn   = 0.d0
Qjac_pnn  = 0.d0
QvmsF_p   = 0.d0
QvmsF_n   = 0.d0
QvmsF_k   = 0.d0
QvmsF_kn  = 0.d0
QvmsAd_p  = 0.d0
QvmsAd_n  = 0.d0
QvmsAd_k  = 0.d0
QvmsAd_kn = 0.d0
Pvec_prev = 0.d0
Qvec_p    = 0.d0
Qvec_k    = 0.d0

!Pressure coupling should be 0
aux_divPIR_perp = 0.d0 
aux_divPIZ_perp = 0.d0     
aux_divPIp_perp = 0.d0
! --- Implicit scheme
theta = time_evol_theta
!zeta  = time_evol_zeta
! change zeta for variable dt
zeta  = time_evol_zeta * 2.0d0 * tstep / (tstep + tstep_prev)

! --- Do we need to use the FFT or non-FFT version?
if ( (i_tor_min == 1) .and. (i_tor_max == n_tor) ) then
  ! In case of global matrix construction:
  use_fft = n_tor > n_tor_fft_thresh
else
  ! In case of "direct construction" of harmonic matrix never FFT:
  use_fft = .false.
end if

if ( use_fft ) then
  ! In case of FFT, don't loop over toroidal harmonics:
  n_tor_start = 1
  n_tor_end   = 1
else
  n_tor_start = i_tor_min
  n_tor_end   = i_tor_max
end if

n_tor_local = n_tor_end - n_tor_start +1

! --- Toroidal basis functions
if (use_fft) then
  ! --- Not needed in case of FFT
  HHZ    = 1.d0
  HHZ_p  = 1.d0
  HHZ_pp = 1.d0
else
  ! --- Needed in case of non-FFT
  do in = 1,n_tor
    do mp=1,n_plane
      HHZ   (in,mp) = HZ   (in,mp)
      HHZ_p (in,mp) = HZ_p (in,mp)
      HHZ_pp(in,mp) = HZ_pp(in,mp)
    enddo
  enddo
endif

! --- TG-stabilisation
TG_NUM = 0.25*tstep
TG_NUM(var_AR) = 0.d0
TG_NUM(var_AZ) = 0.d0
TG_NUM(var_A3) = 0.d0

! jorek00$i.rst jorek_restart.rst
TG_NUM            = 0.0
TG_NUM(var_UR)    = 0.25*tstep
TG_NUM(var_UZ)    = 0.25*tstep
TG_NUM(var_Up)    = 0.25*tstep
TG_NUM(var_UR)    = 0.25*MAX(0.1, tstep)
TG_NUM(var_UZ)    = 0.25*MAX(0.1, tstep)
TG_NUM(var_Up)    = 0.25*MAX(0.1, tstep)
TG_NUM_Eq = 0.0
TG_NUM    = 0.0

! --- RZ variables, Equations variables, and GS-Equilibrium variables
x_g  = 0.d0; x_s  = 0.d0; x_t  = 0.d0; x_ss = 0.d0; x_st = 0.d0; x_tt = 0.d0 
y_g  = 0.d0; y_s  = 0.d0; y_t  = 0.d0; y_ss = 0.d0; y_st = 0.d0; y_tt = 0.d0


eq_g = 0.d0; eq_s = 0.d0; eq_t = 0.d0; eq_p = 0.d0; eq_ss = 0.d0; eq_st = 0.d0; eq_tt = 0.d0
eq_aux_g = 0.d0; eq_aux_s = 0.d0; eq_aux_t = 0.d0; eq_aux_p = 0.d0;
eq_pp    = 0.d0; eq_sp    = 0.d0; eq_tp=0.d0

psi_axisym = 0.d0 ; psi_axisym_s = 0.d0 ; psi_axisym_t = 0.d0
delta_g = 0.d0; delta_s = 0.d0; delta_t = 0.d0
Fprofile= 0.d0

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

        Fprofile  (ms,mt) = Fprofile  (ms,mt) + nodes(i)%Fprof_eq(j) * element%size(i,j) * H  (i,j,ms,mt)

        ! --- Equilibrium psi (n=0 only) for sources
        psi_axisym(ms,mt) = psi_axisym(ms,mt) + nodes(i)%values(1,j,var_A3) * element%size(i,j) * H(i,j,ms,mt)  * HZ(1,1)
        psi_axisym_s(ms,mt) = psi_axisym_s(ms,mt) + nodes(i)%values(1,j,var_A3) * element%size(i,j) * H_s(i,j,ms,mt)
        psi_axisym_t(ms,mt) = psi_axisym_t(ms,mt) + nodes(i)%values(1,j,var_A3) * element%size(i,j) * H_t(i,j,ms,mt)
      end do
    end do

    do ms=1, n_gauss
      do mt=1, n_gauss

        do k=1,n_var

          do in=1,n_tor
            do mp=1,n_plane
              ! --- store variables
              eq_g(mp,k,ms,mt) = eq_g(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H(i,j,ms,mt)  * HZ(in,mp)
              eq_s(mp,k,ms,mt) = eq_s(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_s(i,j,ms,mt)* HZ(in,mp)
              eq_t(mp,k,ms,mt) = eq_t(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_t(i,j,ms,mt)* HZ(in,mp)
              eq_p(mp,k,ms,mt) = eq_p(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H(i,j,ms,mt)  * HZ_p(in,mp)

              eq_ss(mp,k,ms,mt) = eq_ss(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_ss(i,j,ms,mt)* HZ(in,mp)
              eq_st(mp,k,ms,mt) = eq_st(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_st(i,j,ms,mt)* HZ(in,mp)
              eq_tt(mp,k,ms,mt) = eq_tt(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_tt(i,j,ms,mt)* HZ(in,mp)
              
              eq_pp(mp,k,ms,mt) = eq_pp(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H(i,j,ms,mt)  * HZ_pp(in,mp)
              eq_sp(mp,k,ms,mt) = eq_sp(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_s(i,j,ms,mt)* HZ_p(in,mp)
              eq_tp(mp,k,ms,mt) = eq_tp(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_t(i,j,ms,mt)* HZ_p(in,mp)

              if (present(aux_nodes)) then
                eq_aux_g(mp,k,ms,mt) =  eq_aux_g(mp,k,ms,mt) + aux_nodes(i)%values(in,j,k) * element%size(i,j) * H(i,j,ms,mt)   * HZ(in,mp)
                eq_aux_s(mp,k,ms,mt) =  eq_aux_s(mp,k,ms,mt) + aux_nodes(i)%values(in,j,k) * element%size(i,j) * H_s(i,j,ms,mt) * HZ(in,mp)
                eq_aux_t(mp,k,ms,mt) =  eq_aux_t(mp,k,ms,mt) + aux_nodes(i)%values(in,j,k) * element%size(i,j) * H_t(i,j,ms,mt) * HZ(in,mp)
                eq_aux_p(mp,k,ms,mt) =  eq_aux_p(mp,k,ms,mt) + aux_nodes(i)%values(in,j,k) * element%size(i,j) * H(i,j,ms,mt)   * HZ_p(in,mp)
              endif

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

! --- Sources
! --- Note about the current sources:
! --- The toroidal current source can be taken from the routine current.f90, as usual.
! --- The JR and JZ current sources, however, need to be calculated using the initial Grad-Shafranov equilibrium
! --- At t=0, with GS equilibrium, we have:
! --- JR = + psi_Z * dF_dpsi / R = + (psi_bnd_init - psi_axis_init) * psi_norm_Z * dF_dpsi / R
! --- JZ = - psi_R * dF_dpsi / R = - (psi_bnd_init - psi_axis_init) * psi_norm_R * dF_dpsi / R
! --- Thus, these are our current sources as a funcion of psi_norm. Now, at any time, the psi-map will change
! --- And so denormalising psi_norm again, we get:
! --- JR = + (psi_bnd_init - psi_axis_init) / (psi_bnd - psi_axis) * psi_Z * dF_dpsi / R
! --- JZ = - (psi_bnd_init - psi_axis_init) / (psi_bnd - psi_axis) * psi_R * dF_dpsi / R
current_source_JR = 0.d0
current_source_JZ = 0.d0
current_source_Jp = 0.d0
particle_source   = 0.d0
heat_source       = 0.d0
V_source          = 0.d0
dV_dpsi_source    = 0.d0
dV_dz_source      = 0.d0
dF_dpsi           = 0.d0
do ms=1, n_gauss
  do mt=1, n_gauss
    ! --- These are just temporary, only for the sources (they are recalculated afterwards at each n_plane)
    R     = x_g(ms,mt)
    Z     = y_g(ms,mt)
    xjac  = x_s(ms,mt)*y_t(ms,mt) - x_t(ms,mt)*y_s(ms,mt)
    ! --- This will be used in the main loop!
    psi_axisym_R(ms,mt) = (   y_t(ms,mt) * psi_axisym_s(ms,mt)  - y_s(ms,mt) * psi_axisym_t(ms,mt) ) / xjac
    psi_axisym_Z(ms,mt) = ( - x_t(ms,mt) * psi_axisym_s(ms,mt)  + x_s(ms,mt) * psi_axisym_t(ms,mt) ) / xjac
    ! --- These are just temporary, only for the sources (they are recalculated afterwards at each n_plane)
    A30   = eq_g(1,var_A3,ms,mt)
    A30_s = eq_s(1,var_A3,ms,mt)
    A30_t = eq_t(1,var_A3,ms,mt)
    A30_R = (   y_t(ms,mt) * A30_s  - y_s(ms,mt) * A30_t ) / xjac
    A30_Z = ( - x_t(ms,mt) * A30_s  + x_s(ms,mt) * A30_t ) / xjac
    ! --- The dF_dpsi function calculated on time-dependent psi_norm
    ! --- Note: Fprof is be taken from the node values (cleaner)
    call F_profile(xpoint2, xcase2, Z, Z_xpoint, psi_axisym(ms,mt), psi_axis, psi_bnd, &
                   Fprof_time_dep,dF_dpsi(ms,mt) ,dF_dz      ,dF_dpsi2, dF_dz2      ,dF_dpsi_dz , &
                   zFFprime,dFFprime_dpsi,dFFprime_dz,dFFprime_dpsi2,dFFprime_dz2,dFFprime_dpsi_dz)
    if (keep_current_prof) then
      ! --- Toroidal current source. Historically JOREK uses a negative current, so we need to reverse it.
      call current(xpoint2, xcase2, R,Z, Z_xpoint, psi_axisym(ms,mt),psi_axis,psi_bnd,current_source_Jp(ms,mt))
      current_source_Jp(ms,mt) = - current_source_Jp(ms,mt)
      ! --- Poloidal current sources
      current_source_JR(ms,mt) = + (ES%psi_bnd_init - ES%psi_axis_init) / (psi_bnd - psi_axis) * psi_axisym_Z(ms,mt) * dF_dpsi(ms,mt) / R
      current_source_JZ(ms,mt) = - (ES%psi_bnd_init - ES%psi_axis_init) / (psi_bnd - psi_axis) * psi_axisym_R(ms,mt) * dF_dpsi(ms,mt) / R
      if (eta_ARAZ_simple) then
        current_source_JR(ms,mt) = 0.d0
        current_source_JZ(ms,mt) = 0.d0
      endif
    endif
    call sources(xpoint2, xcase2, Z, Z_xpoint, psi_axisym(ms,mt),psi_axis,psi_bnd,particle_source(ms,mt),heat_source(ms,mt))
    ! --- Bootstrap current 
    if (bootstrap) then
      T0     = eq_g(1,var_T,ms,mt)
      T0_s   = eq_s(1,var_T,ms,mt)
      T0_t   = eq_t(1,var_T,ms,mt)
      T0_R   = (   y_t(ms,mt) * T0_s  - y_s(ms,mt) * T0_t ) / xjac
      T0_Z   = ( - x_t(ms,mt) * T0_s  + x_s(ms,mt) * T0_t ) / xjac
      rho0   = eq_g(1,var_rho,ms,mt)
      rho0_s = eq_s(1,var_rho,ms,mt)
      rho0_t = eq_t(1,var_rho,ms,mt)
      rho0_R = (   y_t(ms,mt) * rho0_s  - y_s(ms,mt) * rho0_t ) / xjac
      rho0_Z = ( - x_t(ms,mt) * rho0_s  + x_s(ms,mt) * rho0_t ) / xjac
      ! --- Full Sauter formula
      Ti0   = T0   / 2.d0 ; Te0   = T0   / 2.d0
      Ti0_R = T0_R / 2.d0 ; Te0_R = T0_R / 2.d0
      Ti0_Z = T0_Z / 2.d0 ; Te0_Z = T0_Z / 2.d0
      call bootstrap_current(R, Z,                                                        &
                             R_axis,   Z_axis,   psi_axis,                                &
                             R_xpoint, Z_xpoint, psi_bnd, psi_norm,                       &
                             psi_axisym(ms,mt), psi_axisym_R(ms,mt), psi_axisym_Z(ms,mt), &
                             rho0, rho0_R, rho0_Z,                                        &
                             Ti0,  Ti0_R,  Ti0_Z,                                         &
                             Te0,  Te0_R,  Te0_Z,                                         &
                             Jb)
      ! --- Full Sauter formula for initial profiles
      call density    (xpoint2, xcase2, Z, Z_xpoint, psi_axisym(ms,mt),psi_axis,psi_bnd, &
                       rho_initial(ms,mt),dn_dpsi,dn_dz,dn_dpsi2,dn_dz2,dn_dpsi_dz,dn_dpsi3,dn_dpsi_dz2, dn_dpsi2_dz)
      call temperature(xpoint2, xcase2, y_g(ms,mt), Z_xpoint, eq_g(1,1,ms,mt),psi_axis,psi_bnd, &
                       T_initial  (ms,mt),dT_dpsi,dT_dz,dT_dpsi2,dT_dz2,dT_dpsi_dz,dT_dpsi3,dT_dpsi_dz2, dT_dpsi2_dz)
      Ti0   = T_initial(ms,mt)       / 2.d0
      Ti0_R = dT_dpsi * psi_axisym_R(ms,mt) / 2.d0
      Ti0_Z = dT_dpsi * psi_axisym_Z(ms,mt) / 2.d0
      Te0   = Ti0  
      Te0_R = Ti0_R
      Te0_Z = Ti0_Z
      rho0   = rho_initial(ms,mt)
      rho0_R = dn_dpsi * psi_axisym_R(ms,mt)
      rho0_Z = dn_dpsi * psi_axisym_Z(ms,mt)
      call bootstrap_current(R, Z,                                                        &
                             R_axis,   Z_axis,   psi_axis,                                &
                             R_xpoint, Z_xpoint, psi_bnd, psi_norm,                       &
                             psi_axisym(ms,mt), psi_axisym_R(ms,mt), psi_axisym_Z(ms,mt), &
                             rho0, rho0_R, rho0_Z,                                        &
                             Ti0,  Ti0_R,  Ti0_Z,                                         &
                             Te0,  Te0_R,  Te0_Z,                                         &
                             Jb_0)
      ! --- Subtract the initial equilibrium part
      Jb = Jb - Jb_0
      ! --- Historically JOREK uses a negative current, so we need to reverse it.
      Jb = -Jb
      ! --- Add bootstrap to main current source
      current_source_Jp(ms,mt) = current_source_Jp(ms,mt) + Jb
    endif
    ! --- Neoclassical friction
    if ( NEO ) then 
      if (num_neo_file) then
        call neo_coef(xpoint2, xcase2, Z, Z_xpoint, psi_axisym(ms,mt),psi_axis,psi_bnd, amu_neo_prof(ms,mt), aki_neo_prof(ms,mt))
      else
        amu_neo_prof(ms,mt) = amu_neo_const
        aki_neo_prof(ms,mt) = aki_neo_const
      endif
      psi_norm = get_psi_n(psi_axisym(ms,mt), Z)
      if (psi_norm .gt. 0.995) amu_neo_prof(ms,mt) = 0.d0
      if (psi_norm .lt. 0.100) amu_neo_prof(ms,mt) = 0.d0
    endif
    ! --- Source of toroidal velocity
    if ( ( abs(V_0) .ge. 1.e-12 ) .or. ( num_rot ) ) then
      call velocity(xpoint2, xcase2, Z, Z_xpoint, psi_axisym(ms,mt), psi_axis, psi_bnd, &
                    V_source(ms,mt), dV_dpsi_source(ms,mt),dV_dz_source(ms,mt),         &
                    dV_dpsi2,dV_dz2,dV_dpsi_dz,dV_dpsi3,dV_dpsi_dz2, dV_dpsi2_dz)
    else
      V_source(ms,mt)       = 0.d0
      dV_dpsi_source(ms,mt) = 0.d0
      dV_dz_source(ms,mt)   = 0.d0
    endif

  enddo
enddo


! --- Main loops
do i=1,n_vertex_max
  do j=1,n_degrees
    ELM_p(:,:,1:n_var)  = 0
    ELM_n(:,:,1:n_var)  = 0
    ELM_k(:,:,1:n_var)  = 0
    ELM_kn(:,:,1:n_var) = 0

    do ms=1, n_gauss
      do mt=1, n_gauss

        wst = wgauss(ms)*wgauss(mt)

        xjac    = x_s(ms,mt)*y_t(ms,mt) - x_t(ms,mt)*y_s(ms,mt)
        xjac_R  = (x_ss(ms,mt)*y_t(ms,mt)**2 - y_ss(ms,mt)*x_t(ms,mt)*y_t(ms,mt) - 2.d0*x_st(ms,mt)*y_s(ms,mt)*y_t(ms,mt)   &
                + y_st(ms,mt)*(x_s(ms,mt)*y_t(ms,mt) + x_t(ms,mt)*y_s(ms,mt))                                               &
                + x_tt(ms,mt)*y_s(ms,mt)**2 - y_tt(ms,mt)*x_s(ms,mt)*y_s(ms,mt)) / xjac
        xjac_Z  = (y_tt(ms,mt)*x_s(ms,mt)**2 - x_tt(ms,mt)*y_s(ms,mt)*x_s(ms,mt) - 2.d0*y_st(ms,mt)*x_t(ms,mt)*x_s(ms,mt)   &
                + x_st(ms,mt)*(y_t(ms,mt)*x_s(ms,mt) + y_s(ms,mt)*x_t(ms,mt))                                               &
                + y_ss(ms,mt)*x_t(ms,mt)**2 - x_ss(ms,mt)*y_t(ms,mt)*x_t(ms,mt)) / xjac

        R = x_g(ms,mt)
        Z = y_g(ms,mt)

        do mp = 1, n_plane

          ! --- AR
          AR0   = eq_g(mp,var_AR,ms,mt)
          AR0_p = eq_p(mp,var_AR,ms,mt)
          AR0_s = eq_s(mp,var_AR,ms,mt)
          AR0_t = eq_t(mp,var_AR,ms,mt)
          AR0_R = (   y_t(ms,mt) * AR0_s  - y_s(ms,mt) * AR0_t ) / xjac
          AR0_Z = ( - x_t(ms,mt) * AR0_s  + x_s(ms,mt) * AR0_t ) / xjac
          AR0_ss = eq_ss(mp,var_AR,ms,mt)
          AR0_st = eq_st(mp,var_AR,ms,mt)
          AR0_tt = eq_tt(mp,var_AR,ms,mt)
          AR0_RR = (AR0_ss * y_t(ms,mt)**2 - 2.d0*AR0_st * y_s(ms,mt)*y_t(ms,mt) + AR0_tt * y_s(ms,mt)**2 &
                 + AR0_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) ) &
                 + AR0_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) ) ) / xjac**2             &
                 - xjac_R * (AR0_s* y_t(ms,mt) - AR0_t * y_s(ms,mt))  / xjac**2
          AR0_ZZ = (AR0_ss * x_t(ms,mt)**2 - 2.d0*AR0_st * x_s(ms,mt)*x_t(ms,mt) + AR0_tt * x_s(ms,mt)**2 &
                 + AR0_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) ) &
                 + AR0_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) ) ) / xjac**2             &
                 - xjac_Z * (- AR0_s * x_t(ms,mt) + AR0_t * x_s(ms,mt) )  / xjac**2

          AR0_sp = eq_sp(mp,var_AR,ms,mt)
          AR0_tp = eq_tp(mp,var_AR,ms,mt)
          AR0_pp = eq_pp(mp,var_AR,ms,mt)
          AR0_Rp = (   y_t(ms,mt) * AR0_sp  - y_s(ms,mt) * AR0_tp ) / xjac
          AR0_Zp = ( - x_t(ms,mt) * AR0_sp  + x_s(ms,mt) * AR0_tp ) / xjac
          AR0_RZ = (- AR0_ss * y_t(ms,mt)*x_t(ms,mt) - AR0_tt * x_s(ms,mt)*y_s(ms,mt) &
                 + AR0_st * (y_s(ms,mt)*x_t(ms,mt)  + y_t(ms,mt)*x_s(ms,mt)  ) &
                 - AR0_s  * (x_st(ms,mt)*y_t(ms,mt) - x_tt(ms,mt)*y_s(ms,mt) ) &
                 - AR0_t * (x_st(ms,mt)*y_s(ms,mt)  - x_ss(ms,mt)*y_t(ms,mt) ) ) / xjac**2             &
                 - xjac_R * (- AR0_s * x_t(ms,mt) + AR0_t * x_s(ms,mt) )   / xjac**2

          ! --- AZ
          AZ0   = eq_g(mp,var_AZ,ms,mt)
          AZ0_p = eq_p(mp,var_AZ,ms,mt)
          AZ0_s = eq_s(mp,var_AZ,ms,mt)
          AZ0_t = eq_t(mp,var_AZ,ms,mt)
          AZ0_R = (   y_t(ms,mt) * AZ0_s  - y_s(ms,mt) * AZ0_t ) / xjac
          AZ0_Z = ( - x_t(ms,mt) * AZ0_s  + x_s(ms,mt) * AZ0_t ) / xjac
          AZ0_ss = eq_ss(mp,var_AZ,ms,mt)
          AZ0_st = eq_st(mp,var_AZ,ms,mt)
          AZ0_tt = eq_tt(mp,var_AZ,ms,mt)
          AZ0_RR = (AZ0_ss * y_t(ms,mt)**2 - 2.d0*AZ0_st * y_s(ms,mt)*y_t(ms,mt) + AZ0_tt * y_s(ms,mt)**2 &
                 + AZ0_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) ) &
                 + AZ0_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) ) ) / xjac**2             &
                 - xjac_R * (AZ0_s* y_t(ms,mt) - AZ0_t * y_s(ms,mt))  / xjac**2
          AZ0_ZZ = (AZ0_ss * x_t(ms,mt)**2 - 2.d0*AZ0_st * x_s(ms,mt)*x_t(ms,mt) + AZ0_tt * x_s(ms,mt)**2 &
                 + AZ0_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) ) &
                 + AZ0_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) ) ) / xjac**2             &
                 - xjac_Z * (- AZ0_s * x_t(ms,mt) + AZ0_t * x_s(ms,mt) )  / xjac**2

          AZ0_sp = eq_sp(mp,var_AZ,ms,mt)
          AZ0_tp = eq_tp(mp,var_AZ,ms,mt)
          AZ0_pp = eq_pp(mp,var_AZ,ms,mt)
          AZ0_Rp = (   y_t(ms,mt) * AZ0_sp  - y_s(ms,mt) * AZ0_tp ) / xjac
          AZ0_Zp = ( - x_t(ms,mt) * AZ0_sp  + x_s(ms,mt) * AZ0_tp ) / xjac
          AZ0_RZ = (- AZ0_ss * y_t(ms,mt)*x_t(ms,mt) - AZ0_tt * x_s(ms,mt)*y_s(ms,mt) &
                 + AZ0_st * (y_s(ms,mt)*x_t(ms,mt)  + y_t(ms,mt)*x_s(ms,mt)  ) &
                 - AZ0_s  * (x_st(ms,mt)*y_t(ms,mt) - x_tt(ms,mt)*y_s(ms,mt) ) &
                 - AZ0_t * (x_st(ms,mt)*y_s(ms,mt)  - x_ss(ms,mt)*y_t(ms,mt) ) ) / xjac**2             &
                 - xjac_R * (- AZ0_s * x_t(ms,mt) + AZ0_t * x_s(ms,mt) )   / xjac**2

          ! --- A3:=psi is defined as: A = ... + A03 * grad(phi)
          ! --- as opposed to the magnetic and velocity fieds that are defined as
          ! --- B = ... + Bp * e_phi             and     V = ... + Vp * e_phi       ie.
          ! --- B = ... + Bp * R * grad(phi)     and     V = ... + Vp * R * grad(phi)
          A30   = eq_g(mp,var_A3,ms,mt)
          A30_p = eq_p(mp,var_A3,ms,mt)
          A30_s = eq_s(mp,var_A3,ms,mt)
          A30_t = eq_t(mp,var_A3,ms,mt)
          A30_R = (   y_t(ms,mt) * A30_s  - y_s(ms,mt) * A30_t ) / xjac
          A30_Z = ( - x_t(ms,mt) * A30_s  + x_s(ms,mt) * A30_t ) / xjac
          A30_ss = eq_ss(mp,var_A3,ms,mt)
          A30_st = eq_st(mp,var_A3,ms,mt)
          A30_tt = eq_tt(mp,var_A3,ms,mt)
          A30_RR = (A30_ss * y_t(ms,mt)**2 - 2.d0*A30_st * y_s(ms,mt)*y_t(ms,mt) + A30_tt * y_s(ms,mt)**2 &
                 + A30_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) ) &
                 + A30_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) ) ) / xjac**2             &
                 - xjac_R * (A30_s* y_t(ms,mt) - A30_t * y_s(ms,mt))  / xjac**2
          A30_ZZ = (A30_ss * x_t(ms,mt)**2 - 2.d0*A30_st * x_s(ms,mt)*x_t(ms,mt) + A30_tt * x_s(ms,mt)**2 &
                 + A30_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) ) &
                 + A30_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) ) ) / xjac**2             &
                 - xjac_Z * (- A30_s * x_t(ms,mt) + A30_t * x_s(ms,mt) )  / xjac**2

          A30_sp = eq_sp(mp,var_A3,ms,mt)
          A30_tp = eq_tp(mp,var_A3,ms,mt)
          A30_pp = eq_pp(mp,var_A3,ms,mt)
          A30_Rp = (   y_t(ms,mt) * A30_sp  - y_s(ms,mt) * A30_tp ) / xjac
          A30_Zp = ( - x_t(ms,mt) * A30_sp  + x_s(ms,mt) * A30_tp ) / xjac
          A30_RZ = (- A30_ss * y_t(ms,mt)*x_t(ms,mt) - A30_tt * x_s(ms,mt)*y_s(ms,mt)                     &
                 + A30_st * (y_s(ms,mt)*x_t(ms,mt)  + y_t(ms,mt)*x_s(ms,mt)  ) &
                 - A30_s  * (x_st(ms,mt)*y_t(ms,mt) - x_tt(ms,mt)*y_s(ms,mt) ) &
                 - A30_t * (x_st(ms,mt)*y_s(ms,mt)  - x_ss(ms,mt)*y_t(ms,mt) ) ) / xjac**2             &
                 - xjac_R * (- A30_s * x_t(ms,mt) + A30_t * x_s(ms,mt) )   / xjac**2

          ! --- UR
          UR0   = eq_g(mp,var_UR,ms,mt)
          UR0_p = eq_p(mp,var_UR,ms,mt)
          UR0_s = eq_s(mp,var_UR,ms,mt)
          UR0_t = eq_t(mp,var_UR,ms,mt)
          UR0_R = (   y_t(ms,mt) * UR0_s  - y_s(ms,mt) * UR0_t ) / xjac
          UR0_Z = ( - x_t(ms,mt) * UR0_s  + x_s(ms,mt) * UR0_t ) / xjac
          UR0_ss = eq_ss(mp,var_UR,ms,mt)
          UR0_st = eq_st(mp,var_UR,ms,mt)
          UR0_tt = eq_tt(mp,var_UR,ms,mt)
          UR0_RR =  (   UR0_ss * y_t(ms,mt)**2 - 2.d0*UR0_st * y_s(ms,mt)*y_t(ms,mt)              &
                      + UR0_tt * y_s(ms,mt)**2                                                    &
                      + UR0_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                &
                      + UR0_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) )   )/ xjac**2   &
                      - UR0_R * xjac_R / xjac
          UR0_ZZ =  (   UR0_ss * x_t(ms,mt)**2 - 2.d0*UR0_st * x_s(ms,mt)*x_t(ms,mt)              &
                      + UR0_tt * x_s(ms,mt)**2                                                    &
                      + UR0_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                &
                      + UR0_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) )   )/ xjac**2   &
                      - UR0_Z * xjac_Z / xjac

          ! --- UZ
          UZ0   = eq_g(mp,var_UZ,ms,mt)
          UZ0_p = eq_p(mp,var_UZ,ms,mt)
          UZ0_s = eq_s(mp,var_UZ,ms,mt)
          UZ0_t = eq_t(mp,var_UZ,ms,mt)
          UZ0_R = (   y_t(ms,mt) * UZ0_s  - y_s(ms,mt) * UZ0_t ) / xjac
          UZ0_Z = ( - x_t(ms,mt) * UZ0_s  + x_s(ms,mt) * UZ0_t ) / xjac
          UZ0_ss = eq_ss(mp,var_UZ,ms,mt)
          UZ0_st = eq_st(mp,var_UZ,ms,mt)
          UZ0_tt = eq_tt(mp,var_UZ,ms,mt)
          UZ0_RR =  (   UZ0_ss * y_t(ms,mt)**2 - 2.d0*UZ0_st * y_s(ms,mt)*y_t(ms,mt)              &
                      + UZ0_tt * y_s(ms,mt)**2                                                    &
                      + UZ0_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                &
                      + UZ0_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) )   )/ xjac**2   &
                      - UZ0_R * xjac_R / xjac
          UZ0_ZZ =  (   UZ0_ss * x_t(ms,mt)**2 - 2.d0*UZ0_st * x_s(ms,mt)*x_t(ms,mt)              &
                      + UZ0_tt * x_s(ms,mt)**2                                                    &
                      + UZ0_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                &
                      + UZ0_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) )   )/ xjac**2   &
                      - UZ0_Z * xjac_Z / xjac

          ! --- Up is defined u0_phi : V = .. + Up0 * e_phi (physical component)
          Up0   = eq_g(mp,var_Up,ms,mt)
          Up0_p = eq_p(mp,var_Up,ms,mt)
          Up0_s = eq_s(mp,var_Up,ms,mt)
          Up0_t = eq_t(mp,var_Up,ms,mt)
          Up0_R = (   y_t(ms,mt) * Up0_s  - y_s(ms,mt) * Up0_t ) / xjac
          Up0_Z = ( - x_t(ms,mt) * Up0_s  + x_s(ms,mt) * Up0_t ) / xjac
          Up0_ss = eq_ss(mp,var_Up,ms,mt)
          Up0_st = eq_st(mp,var_Up,ms,mt)
          Up0_tt = eq_tt(mp,var_Up,ms,mt)
          Up0_RR =  (   Up0_ss * y_t(ms,mt)**2 - 2.d0*Up0_st * y_s(ms,mt)*y_t(ms,mt)              &
                      + Up0_tt * y_s(ms,mt)**2                                                    &
                      + Up0_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                &
                      + Up0_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) )   )/ xjac**2   &
                      - Up0_R * xjac_R / xjac
          Up0_ZZ =  (   Up0_ss * x_t(ms,mt)**2 - 2.d0*Up0_st * x_s(ms,mt)*x_t(ms,mt)              &
                      + Up0_tt * x_s(ms,mt)**2                                                    &
                      + Up0_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                &
                      + Up0_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) )   )/ xjac**2   &
                      - Up0_Z * xjac_Z / xjac

          ! --- rho
          rho0      = eq_g(mp,var_rho,ms,mt)
          rho0_corr = max(rho0,1.d-12)!corr_neg_dens1(rho0) ! CAREFUL! FULL-MHD DOESN'T LIKE THE CORR FUNCTIONS AT ALL
          rho0_p    = eq_p(mp,var_rho,ms,mt)
          rho0_s    = eq_s(mp,var_rho,ms,mt)
          rho0_t    = eq_t(mp,var_rho,ms,mt)
          rho0_R    = (   y_t(ms,mt) * rho0_s  - y_s(ms,mt) * rho0_t ) / xjac
          rho0_Z    = ( - x_t(ms,mt) * rho0_s  + x_s(ms,mt) * rho0_t ) / xjac
          rho0_ss   = eq_ss(mp,var_rho,ms,mt)
          rho0_st   = eq_st(mp,var_rho,ms,mt)
          rho0_tt   = eq_tt(mp,var_rho,ms,mt)
          rho0_RR   = ( rho0_ss * y_t(ms,mt)**2 - 2.d0*rho0_st * y_s(ms,mt)*y_t(ms,mt)              &
                      + rho0_tt * y_s(ms,mt)**2 &
                      + rho0_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                &
                      + rho0_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) )   )/ xjac**2   &
                      - rho0_R * xjac_R / xjac
          rho0_ZZ   = ( rho0_ss * x_t(ms,mt)**2 - 2.d0*rho0_st * x_s(ms,mt)*x_t(ms,mt)              &
                      + rho0_tt * x_s(ms,mt)**2 &
                      + rho0_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                &
                      + rho0_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) )   )/ xjac**2   &
                      - rho0_Z * xjac_Z / xjac

          ! --- T
          T0      = eq_g(mp,var_T,ms,mt)
          T0_corr = max(T0,1.d-12)!corr_neg_temp1(T0) ! CAREFUL! FULL-MHD DOESN'T LIKE THE CORR FUNCTIONS AT ALL
          T0_p    = eq_p(mp,var_T,ms,mt)
          T0_s    = eq_s(mp,var_T,ms,mt)
          T0_t    = eq_t(mp,var_T,ms,mt)
          T0_R    = (   y_t(ms,mt) * T0_s  - y_s(ms,mt) * T0_t ) / xjac
          T0_Z    = ( - x_t(ms,mt) * T0_s  + x_s(ms,mt) * T0_t ) / xjac
          T0_ss   = eq_ss(mp,var_T,ms,mt)
          T0_st   = eq_st(mp,var_T,ms,mt)
          T0_tt   = eq_tt(mp,var_T,ms,mt)
          T0_RR   = ( T0_ss * y_t(ms,mt)**2 - 2.d0*T0_st * y_s(ms,mt)*y_t(ms,mt)              &
                      + T0_tt * y_s(ms,mt)**2 &
                      + T0_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                &
                      + T0_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) )   )/ xjac**2   &
                      - T0_R * xjac_R / xjac
          T0_ZZ   = ( T0_ss * x_t(ms,mt)**2 - 2.d0*T0_st * x_s(ms,mt)*x_t(ms,mt)              &
                      + T0_tt * x_s(ms,mt)**2 &
                      + T0_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                &
                      + T0_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) )   )/ xjac**2   &
                      - T0_Z * xjac_Z / xjac

          ! --- P
          p0      = rho0 * T0
          p0_corr = rho0_corr * T0_corr
          p0_R    = rho0_R * T0 + rho0 * T0_R
          p0_Z    = rho0_Z * T0 + rho0 * T0_Z
          p0_s    = rho0_s * T0 + rho0 * T0_s
          p0_t    = rho0_t * T0 + rho0 * T0_t
          p0_p    = rho0_p * T0 + rho0 * T0_p

          ! --- psi_norm
          psi_norm = get_psi_n(psi_axisym(ms,mt), y_g(ms,mt))

          ! --- Diffusions
          D_prof  = get_dperp (psi_norm)
          ZK_prof = get_zkperp(psi_norm)

          ! --- Resistivity
          if ( eta_T_dependent .and. T0_corr <= T_max_eta) then
            eta_T     =   eta   * (T0_corr/T_0)**(-1.5d0)
            deta_dT   = - eta   * 1.5d0  * T0_corr**(-2.5d0) * T_0**(1.5d0)
            d2eta_d2T =   eta   * 3.75d0 * T0_corr**(-3.5d0) * T_0**(1.5d0)
          else if ( eta_T_dependent .and. T0_corr > T_max_eta) then
            eta_T     = eta * (T_max_eta/T_0)**(-1.5d0)
            deta_dT   = 0.d0
            d2eta_d2T = 0.d0
          else
            eta_T     = eta
            deta_dT   = 0.d0
            d2eta_d2T = 0.d0
          end if
          if ( eta_T_dependent .and.  xpoint2 .and. (T0 .lt. T_min) ) then
              eta_T     = eta    * (T_min/T_0)**(-1.5d0)
              deta_dT   = 0.d0
              d2eta_d2T = 0.d0
          end if
          eta_R = deta_dT * T0_R
          eta_Z = deta_dT * T0_Z
          eta_p = deta_dT * T0_p

          ! --- Resistivity for Ohmic heating
          if ( eta_T_dependent .and. T0_corr <= T_max_eta) then
            eta_T_ohm   =   eta_ohmic   * (T0_corr/T_0)**(-1.5d0)
            deta_dT_ohm = - eta_ohmic   * 1.5d0  * T0_corr**(-2.5d0) * T_0**(1.5d0)
          else if ( eta_T_dependent .and. T0_corr > T_max_eta) then
            eta_T_ohm   = eta_ohmic * (T_max_eta/T_0)**(-1.5d0)
            deta_dT_ohm = 0.d0
          else
            eta_T_ohm     = eta_ohmic
            deta_dT_ohm   = 0.d0
          end if
          if ( eta_T_dependent .and.  xpoint2 .and. (T0 .lt. T_min) ) then
              eta_T_ohm     = eta_ohmic    * (T_min/T_0)**(-1.5d0)
              deta_dT_ohm   = 0.d0
          end if

          ! --- Viscosity
          if ( visco_T_dependent ) then
            visco_T   = visco * (T0_corr/T_0)**(-1.5d0)
            dvisco_dT = - visco * (1.5d0)  * T0_corr**(-2.5d0) * T_0**(1.5d0)
            if ( xpoint2 .and. (T0 .lt. T_min) ) then
              visco_T     = visco  * (T_min/T_0)**(-1.5d0)
              dvisco_dT   = 0.d0
            endif
          else
            visco_T   = visco
            dvisco_dT = 0.d0
          end if

          ! --- Kpar
          if ( ZKpar_T_dependent ) then
            ZKpar_T   = ZK_par * (T0_corr/T_0)**(+2.5d0)
            dZKpar_dT = ZK_par * (2.5d0)  * T0_corr**(+1.5d0) * T_0**(-2.5d0)
            if (ZKpar_T .gt. ZK_par_max) then
              ZKpar_T   = Zk_par_max
              dZKpar_dT = 0.d0
            endif
            if ( xpoint2 .and. (T0 .lt. T_min_ZKpar) ) then
              ZKpar_T   = ZK_par * (T_min_ZKpar/T_0)**(+2.5d0)
              dZKpar_dT = 0.d0
            endif
          else
            ZKpar_T   = ZK_par
            dZKpar_dT = 0.d0
          endif

          ! --- Magnetic field
          Fprof = Fprofile(ms,mt)
          BR0 = ( A30_Z - AZ0_p )/ R
          BZ0 = ( AR0_p - A30_R )/ R
          Bp0 = ( AZ0_R - AR0_Z )    + Fprof / R
          BB2 = Bp0**2 + BR0**2 + BZ0**2

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

          ! --- Toroidal magnetic field without F contribution (for eta_ARAZ_simple)
          Bp00 = Bp0
          if (eta_ARAZ_simple) Bp00 = ( AZ0_R - AR0_Z )

          ! --- B.grad and V.grad
          BgradT   = BR0 * T0_R   + BZ0 * T0_Z   + Bp0 * T0_p   / R
          UgradT   = UR0 * T0_R   + UZ0 * T0_Z   + Up0 * T0_p   / R
          BgradRho = BR0 * rho0_R + BZ0 * rho0_Z + Bp0 * rho0_p / R
          UgradRho = UR0 * rho0_R + UZ0 * rho0_Z + Up0 * rho0_p / R
          BgradP   = rho0*BgradT + T0*BgradRho

          ! --- V.grad.V
          UgradUR  = UR0 * UR0_R + UZ0 * UR0_Z + Up0 * UR0_p / R
          UgradUZ  = UR0 * UZ0_R + UZ0 * UZ0_Z + Up0 * UZ0_p / R
          UgradUp  = UR0 * Up0_R + UZ0 * Up0_Z + Up0 * Up0_p / R

          ! --- div.V
          divU     = UR0_R + UR0/R + UZ0_Z + Up0_p/R
          divRhoU  = rho0 * divU + UgradRho

          ! --- Anisotropic pressure coupling scheme
          if (use_pcs) then
             aux_PIRR    = eq_aux_g(mp,1,ms,mt);
             aux_PIRR_s  = eq_aux_s(mp,1,ms,mt);
             aux_PIRR_t  = eq_aux_t(mp,1,ms,mt);
             aux_PIRR_R  = (   y_t(ms,mt) * aux_PIRR_s - y_s(ms,mt) * aux_PIRR_t ) / xjac
             aux_PIRR_Z  = ( - x_t(ms,mt) * aux_PIRR_s + x_s(ms,mt) * aux_PIRR_t ) / xjac
             aux_PIRR_p  = eq_aux_p(mp,1,ms,mt);

             aux_PIZZ    = eq_aux_g(mp,2,ms,mt);
             aux_PIZZ_s  = eq_aux_s(mp,2,ms,mt);
             aux_PIZZ_t  = eq_aux_t(mp,2,ms,mt);
             aux_PIZZ_R  = (   y_t(ms,mt) * aux_PIZZ_s - y_s(ms,mt) * aux_PIZZ_t ) / xjac
             aux_PIZZ_Z  = ( - x_t(ms,mt) * aux_PIZZ_s + x_s(ms,mt) * aux_PIZZ_t ) / xjac
             aux_PIZZ_p  = eq_aux_p(mp,2,ms,mt);
             
             aux_PIPP    = eq_aux_g(mp,3,ms,mt);
             aux_PIPP_s  = eq_aux_s(mp,3,ms,mt);
             aux_PIPP_t  = eq_aux_t(mp,3,ms,mt);
             aux_PIPP_R  = (   y_t(ms,mt) * aux_PIPP_s - y_s(ms,mt) * aux_PIPP_t ) / xjac
             aux_PIPP_Z  = ( - x_t(ms,mt) * aux_PIPP_s + x_s(ms,mt) * aux_PIPP_t ) / xjac
             aux_PIPP_p  = eq_aux_p(mp,3,ms,mt);             

             aux_PIZR    = eq_aux_g(mp,4,ms,mt);
             aux_PIZR_s  = eq_aux_s(mp,4,ms,mt);
             aux_PIZR_t  = eq_aux_t(mp,4,ms,mt);
             aux_PIZR_R  = (   y_t(ms,mt) * aux_PIZR_s - y_s(ms,mt) * aux_PIZR_t ) / xjac
             aux_PIZR_Z  = ( - x_t(ms,mt) * aux_PIZR_s + x_s(ms,mt) * aux_PIZR_t ) / xjac
             aux_PIZR_p  = eq_aux_p(mp,4,ms,mt);             

             aux_PIRP    = eq_aux_g(mp,5,ms,mt);
             aux_PIRP_s  = eq_aux_s(mp,5,ms,mt);
             aux_PIRP_t  = eq_aux_t(mp,5,ms,mt);
             aux_PIRP_R  = (   y_t(ms,mt) * aux_PIRP_s - y_s(ms,mt) * aux_PIRP_t ) / xjac
             aux_PIRP_Z  = ( - x_t(ms,mt) * aux_PIRP_s + x_s(ms,mt) * aux_PIRP_t ) / xjac
             aux_PIRP_p  = eq_aux_p(mp,5,ms,mt);             

             aux_PIZP    = eq_aux_g(mp,6,ms,mt);
             aux_PIZP_s  = eq_aux_s(mp,6,ms,mt);
             aux_PIZP_t  = eq_aux_t(mp,6,ms,mt);
             aux_PIZP_R  = (   y_t(ms,mt) * aux_PIZP_s - y_s(ms,mt) * aux_PIZP_t ) / xjac
             aux_PIZP_Z  = ( - x_t(ms,mt) * aux_PIZP_s + x_s(ms,mt) * aux_PIZP_t ) / xjac
             aux_PIZP_p  = eq_aux_p(mp,6,ms,mt);

             !See https://www.jorek.eu/wiki/doku.php?id=coordinates, div Pi (tensor)
             !such that R component=alpha, Z component=beta, phi component=gamma, in physical components
             !These calculations are the same in reduced MHD. Only the way in which these terms enter
             !the projected equations differs.
             aux_divPIR  = aux_PIRR_R + aux_PIZR_Z + (aux_PIRP_p + aux_PIRR - aux_PIPP) / R
             aux_divPIZ  = aux_PIZZ_Z + aux_PIZR_R + (aux_PIZP_p + aux_PIZR) / R
             aux_divPIp  = aux_PIPP_p / R + aux_PIRP_R + aux_PIZP_Z + 2.d0 / R * aux_PIRP
  
             ! B. div P
             aux_BdivPI      = BR0 * aux_divPIR + BZ0 * aux_divPIZ + Bp0 * aux_divPIp
             
             aux_divPIR_perp = aux_divPIR - BR0 / BB2 * aux_BdivPI
             aux_divPIZ_perp = aux_divPIZ - BZ0 / BB2 * aux_BdivPI
             aux_divPIp_perp = aux_divPIp - Bp0 / BB2 * aux_BdivPI
     
             
          endif
          
          ! --- Diamagnetic velocity Vdia = tau / (rho*BB2) * B x grad(p)
          ! --- Note-1: Phi component defined as physical component VdiaP*e_phi, like V and B
          ! --- Note-2: Factor of F0 is here so that we have the same definition of tau_IC in RMHD and FMHD
          tau_IC = tauIC
          ! --- Switch off at targets?
          !distance_bnd = 1.d10
          !do im=1,n_vertex_max
          !  if (nodes(im)%boundary .eq. 1) then
          !    distance_bnd = min(distance_bnd, sqrt((R-nodes(im)%x(1,1,1))**2 + (Z-nodes(im)%x(1,1,2))**2) )
          !  endif
          !enddo
          !tau_IC = tauIC * (0.5d0 - 0.5d0 * tanh(-(distance_bnd - 0.02)/0.01) )
          VdiaR0 = tau_IC*F0 / (R * rho0_corr * BB2) * (  BZ0*p0_p - R*Bp0*p0_Z)
          VdiaZ0 = tau_IC*F0 / (R * rho0_corr * BB2) * (R*BP0*p0_R -   BR0*p0_p)
          VdiaP0 = tau_IC*F0 / (    rho0_corr * BB2) * (  BR0*p0_Z -   BZ0*p0_R)

          ! --- Vdia.grad.V
          VdiaGradUR  = VdiaR0 * UR0_R + VdiaZ0 * UR0_Z + VdiaP0 * UR0_p / R
          VdiaGradUZ  = VdiaR0 * UZ0_R + VdiaZ0 * UZ0_Z + VdiaP0 * UZ0_p / R
          VdiaGradUp  = VdiaR0 * Up0_R + VdiaZ0 * Up0_Z + VdiaP0 * Up0_p / R
          
          ! --- Toroidal velocity source
          Vt0   = V_source(ms,mt)
          Vt0_R = dV_dpsi_source(ms,mt) * psi_axisym_R(ms,mt)
          Vt0_Z = dV_dz_source(ms,mt) + dV_dpsi_source(ms,mt) * psi_axisym_Z(ms,mt)

          ! --- Neoclassical friction, we assume that the magnetic field is from the axisymmetric component only
          if (NEO) then
            Btht  = sqrt(BR0**2 + BZ0**2)
            Btht  = max(1.d-10,Btht)
            Vtht  = ( BR0*(UR0+VdiaR0) + BZ0*(UZ0+VdiaZ0) ) / Btht
            Vneo  = +aki_neo_prof(ms,mt) * tau_IC*F0 / BB2 / Btht * ( BR0*Bp0*T0_Z - BZ0*Bp0*T0_R)
            PneoR = +amu_neo_prof(ms,mt) * rho0 * BB2/Btht**3 * BR0 * (Vtht - Vneo)
            PneoZ = +amu_neo_prof(ms,mt) * rho0 * BB2/Btht**3 * BZ0 * (Vtht - Vneo)
          else
            PneoR = 0.d0
            PneoZ = 0.d0
          endif
          
          ! --- VMS Variables
          VdotB  = UR0*BR0 + Up0*Bp0 + UZ0*BZ0

          VbR0   = VdotB*BR0/BB2
          VbZ0   = VdotB*BZ0/BB2
          Vbp0   = VdotB*Bp0/BB2

          CvR0   = UR0 - VbR0
          CvZ0   = UZ0 - VbZ0
          Cvp0   = Up0 - Vbp0

          CvGradAR0       = CvR0 * AR0_R + CvZ0 * AR0_Z + Cvp0 * AR0_p / R
          CvGradAZ0       = CvR0 * AZ0_R + CvZ0 * AZ0_Z + Cvp0 * AZ0_p / R
          CvGradA30       = CvR0 * A30_R + CvZ0 * A30_Z + Cvp0 * A30_p / R
          
          CvGradUR0       = CvR0 * UR0_R + CvZ0 * UR0_Z + Cvp0 * UR0_p / R
          CvGradUZ0       = CvR0 * UZ0_R + CvZ0 * UZ0_Z + Cvp0 * UZ0_p / R
          CvGradUp0       = CvR0 * Up0_R + CvZ0 * Up0_Z + Cvp0 * Up0_p / R

          CvGradr0        = CvR0 * rho0_R  + CvZ0 * rho0_Z  + Cvp0 * rho0_p / R
          CvGradT0        = CvR0 * T0_R    + CvZ0 * T0_Z    + Cvp0 * T0_p   / R

          VbGradUR0       = VbR0 * UR0_R + VbZ0 * UR0_Z + Vbp0 * UR0_p / R
          VbGradUZ0       = VbR0 * UZ0_R + VbZ0 * UZ0_Z + Vbp0 * UZ0_p / R
          VbGradUp0       = VbR0 * Up0_R + VbZ0 * Up0_Z + Vbp0 * Up0_p / R
          
          VbGradr0        = VbR0 * rho0_R  + VbZ0 * rho0_Z  + Vbp0 * rho0_p / R
          VbGradT0        = VbR0 * T0_R    + VbZ0 * T0_Z    + Vbp0 * T0_p   / R

          ! --- Coeficients for Acoustic waves stabilization (set to zero to remove)
          VmsCoefF        = (gamma * T0 + BB2/MAX(rho0, Rho_min) )
          VmsCoefF_T      = gamma * T0 
          VmsCoefF_T      = 0.0d0

          ! --- Force Div(V)=0
          visco_divV      = 0.d0
          dvisco_divV_dT  = 0.d0
          select case(ViscType)
          case(0)
            visco_divV      =  visco_T + Coef_DivV
            dvisco_divV_dT  =  dvisco_dT
          case(10:15)
            visco_divV      =  Coef_DivV
            dvisco_divV_dT  =  0.0
          case(20)
            visco_divV      =  0.0
            dvisco_divV_dT  =  0.0
          case(30)
            visco_divV      =  visco_T
            dvisco_divV_dT  =  dvisco_dT
          case default
            visco_divV      =  0.0
            dvisco_divV_dT  =  0.0
          end select

          psieq_R = (   y_t(ms,mt) * psi_axisym_s(ms,mt)  - y_s(ms,mt) * psi_axisym_t(ms,mt) ) / xjac
          psieq_Z = ( - x_t(ms,mt) * psi_axisym_s(ms,mt)  + x_s(ms,mt) * psi_axisym_t(ms,mt) ) / xjac

          BR0_R = -1.0d0/R**2 * ( A30_Z - AZ0_p ) + ( A30_RZ - AZ0_Rp )/R
          BR0_Z = ( A30_ZZ - AZ0_Zp ) / R
          BR0_p = ( A30_Zp - AZ0_pp ) / R
          BZ0_R = -1.0d0/R**2 * ( AR0_p - A30_R ) + ( AR0_Rp - A30_RR )/R
          BZ0_Z = ( AR0_Zp - A30_RZ ) / R
          BZ0_p = ( AR0_pp - A30_Rp ) / R
          Bp0_R = ( AZ0_RR - AR0_RZ )  &
                + (psi_bnd - psi_axis) / (ES%psi_bnd_init - ES%psi_axis_init) * psieq_R * dF_dpsi(ms,mt) / R  &
                - Fprofile(ms,mt)/R**2
          Bp0_Z = ( AZ0_RZ - AR0_ZZ ) &
                + (psi_bnd - psi_axis) / (ES%psi_bnd_init - ES%psi_axis_init) * psieq_Z * dF_dpsi(ms,mt) / R
          Bp0_p = ( AZ0_Rp - AR0_Zp )

          JR0 = Bp0_Z - BZ0_p / R
          JZ0 = (BR0_p - R*Bp0_R - Bp0) / R
          Jp0 = BZ0_R - BR0_Z

          JJ2 = JR0*JR0 + JZ0*JZ0 + Jp0*Jp0

          do im=n_tor_start, n_tor_end

            ! --- test functions (V*)
            v   = H(i,j,ms,mt)   * element%size(i,j) * HHZ(im,mp)
            v_s = H_s(i,j,ms,mt) * element%size(i,j) * HHZ(im,mp)
            v_t = H_t(i,j,ms,mt) * element%size(i,j) * HHZ(im,mp)
            v_p = H(i,j,ms,mt)   * element%size(i,j) * HHZ_p(im,mp)
            v_R = (  y_t(ms,mt) * v_s - y_s(ms,mt) * v_t ) / xjac
            v_Z = (- x_t(ms,mt) * v_s + x_s(ms,mt) * v_t ) / xjac 
            v_ss= H_ss(i,j,ms,mt) * element%size(i,j) * HHZ(im,mp)
            v_st= H_st(i,j,ms,mt) * element%size(i,j) * HHZ(im,mp)
            v_tt= H_tt(i,j,ms,mt) * element%size(i,j) * HHZ(im,mp)
            v_RR= ( v_ss * y_t(ms,mt)**2 - 2.d0*v_st * y_s(ms,mt)*y_t(ms,mt) &
                  + v_tt * y_s(ms,mt)**2 &
                  + v_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                &
                  + v_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) )   )/ xjac**2   &
                  - v_R * xjac_R / xjac
            v_ZZ= ( v_ss * x_t(ms,mt)**2 - 2.d0*v_st * x_s(ms,mt)*x_t(ms,mt) &
                  + v_tt * x_s(ms,mt)**2 &
                  + v_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                &
                  + v_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) )   )/ xjac**2   &
                  - v_Z * xjac_Z / xjac

            ! --- Terms for Integration-by-parts
            BgradVstar__p = BR0 * v_R + BZ0 * v_Z
            BgradVstar__k = Bp0 * v_p / R

            UgradVstar__p = UR0 * v_R + UZ0 * v_Z
            UgradVstar__k = Up0 * v_p / R

            gradRho_gradVstar__p = rho0_R * v_R  + rho0_Z * v_Z
            gradRho_gradVstar__k = (rho0_p / R) * (v_p  / R)

            gradT_gradVstar__p = T0_R * v_R  + T0_Z * v_Z
            gradT_gradVstar__k = (T0_p / R) * (v_p  / R)

            ! --- (V.grad)V        normal part                       diamagnetic part
            Qconv_UR = - v * rho0 * ( UgradUR - Up0**2  / R )    - v * rho0 * ( VdiaGradUR - VdiaP0*Up0 / R )
            Qconv_UZ = - v * rho0 * UgradUZ                      - v * rho0 * VdiaGradUZ
            Qconv_Up = - v * rho0 * ( UgradUp + UR0*Up0 / R )    - v * rho0 * ( VdiaGradUp + VdiaP0*UR0 / R )

            JxB_UR__p = - BR0 * BgradVstar__p - Bp0**2 * v / R + ( v_R + v / R ) * ( 0.5d0 * BB2 )
            JxB_UR__k = - BR0 * BgradVstar__k
            JxB_UZ__p = - BZ0 * BgradVstar__p + v_Z * ( 0.5d0 * BB2)
            JxB_UZ__k = - BZ0 * BgradVstar__k
            JxB_Up__p = - Bp0 * BgradVstar__p + Bp0*BR0* v / R
            JxB_Up__k = - Bp0 * BgradVstar__k + ( v_p/ R ) * ( 0.5d0 * BB2 )

            ! --- Diamagnetic terms
            VdiaGradVstar__p = VdiaR0 * v_R + VdiaZ0 * v_Z
            VdiaGradVstar__k = VdiaP0 * v_p / R

            ! --- VMS variables
            CvGradVi__p     = CvR0 * v_R   + CvZ0 * v_Z
            CvGradVi__k     = Cvp0 * v_p / R
            VbGradVi__p     = VbR0 * v_R   + VbZ0 * v_Z
            VbGradVi__k     = Vbp0 * v_p / R

            DiveRMVi    = v_R + v / R
            DiveZMVi    = v_z
            DivePMVi__k = v_p / R

            ! --- Viscous terms
            Qvisc_UR__p = 0.d0 ; Qvisc_UR__k = 0.d0
            Qvisc_UZ__p = 0.d0 ; Qvisc_UZ__k = 0.d0
            Qvisc_Up__p = 0.d0 ; Qvisc_Up__k = 0.d0
            Qvisc_T     = 0.d0
            select case(ViscType)
            case(0) ! --- Full-MHD viscosity          
              Qvisc_UR__p = - (v_R * UR0_R + v_Z * UR0_Z) - v * ( UR0 / R**2 + 2.d0 * Up0_p / R**2 ) 
              Qvisc_UR__k = - (v_p * UR0_p / R**2)
              Qvisc_UZ__p = - (v_R * UZ0_R + v_Z * UZ0_Z)
              Qvisc_UZ__k = - (v_p * UZ0_p / R**2)
              Qvisc_Up__p = - (v_R * UP0_R + v_Z * UP0_Z) - v * ( Up0 / R**2 - 2.d0 * UR0_p / R**2 )
              Qvisc_Up__k = - (v_p * uP0_p / R**2)
            case(15) ! --- Mock-up of reduced-MHD viscosity
              Qvisc_UR__p = + ( UZ0_R - UR0_Z ) * v_Z
              Qvisc_UZ__p = - ( UZ0_R - UR0_Z ) * v_R
              Qvisc_Up__p = - ( Up0_R + Up0/R ) * ( v_R + v/R ) - Up0_z * v_z 
            case(10) ! --- Generic form of full viscosity tensor (from B.Nkonga and A.Bhole)
              Qvisc_UR__p = - ( v_R * UR0_R + v_Z * UR0_Z )   &
                            - ( v_R * UR0_R + v_Z * UZ0_R )   &
                            - ( v * UR0  + v *  Up0_p )/ R**2 &
                            - v * ( UR0 +  Up0_R )/R**2
              Qvisc_UR__k = - ( v_p * UR0_p / R**2 ) &
                            - ( v_p * Up0_R / R    ) &
                            - ( - v_p * Up0 )/ R**2
              
              Qvisc_UZ__p = - ( v_R * UZ0_R + v_Z * UZ0_Z ) &
                            - ( v_R * UR0_Z + v_Z * UZ0_Z )
              Qvisc_UZ__k = - ( v_p * UZ0_p / R**2 ) &
                            - ( v_p * Up0_Z / R    )
              
              Qvisc_Up__p = - ( v_R * UP0_R + v_Z * UP0_Z )   &
                            - ( v_R * UR0_p + v_Z * UZ0_p )/R &
                            - ( - v * UR0_p + v * Up0 )/ R**2 &
                            + ( v_R * Up0 + v * Up0_R )/ R
              Qvisc_Up__k = - ( v_p * uP0_p / R**2) &
                            - ( v_p * Up0_p / R )/R &
                            - ( v_p * UR0 )/ R**2   &
                            - ( v_p * UR0 )/ R**2

              QviscT0     =  - ( UR0_R * UR0_R + UR0_Z * UR0_Z + UR0_p * UR0_p / R**2)      &
                             - ( UR0_R * UR0_R + UR0_Z * UZ0_R + UR0_p * Up0_R / R )        &
                             - ( UR0 * UR0  + UR0 *  Up0_p - UR0_p * Up0 )/ R**2            &
                             - UR0 * ( UR0 +  Up0_R )/R**2                                  &
                             - ( UZ0_R * UZ0_R + UZ0_Z * UZ0_Z + UZ0_p * UZ0_p / R**2)      &
                             - ( UZ0_R * UR0_Z + UZ0_Z * UZ0_Z + UZ0_p * Up0_Z / R )        &
                             - ( Up0_R * UP0_R + Up0_Z * UP0_Z + Up0_p * uP0_p / R**2)      &
                             - ( Up0_R * UR0_p + Up0_Z * UZ0_p + Up0_p * Up0_p / R )/R      &
                             - ( Up0_p * UR0   -   Up0 * UR0_p +   Up0 * Up0 )/ R**2        &
                             - ( Up0_p * UR0 )/ R**2                                        &
                             + ( Up0_R * Up0 +  Up0 * Up0_R )/ R

              Qvisc_T     = visco_divV * divU**2 - visco_T * QviscT0
            case(20) ! --- Viscosity as implemented by W.Haverkort
              Qvisc_UR__p = - (2.d0 * Up0_p * v / R**2 + 2.d0 * UR0_R * v_R + (UR0_Z + UZ0_R ) * v_Z )
              Qvisc_UR__k = - ((UR0_p + R * Up0_R - Up0) * v_p / R**2)

              Qvisc_UZ__p = - (2.d0 * UZ0_Z * v_Z + (UZ0_R + UR0_Z ) * v_R)
              Qvisc_UZ__k = - ((UZ0_p + R * Up0_Z ) * v_p / R**2)

              Qvisc_Up__p = - (+ ( 3.d0 * Up0 / R -Up0_R - UR0_p / R ) * v / R &
                               + ( -Up0 / R + Up0_R + UR0_p / R ) * v_R           &
                               + (R * Up0_Z + UZ0_p) * v_Z / R                    )
              Qvisc_Up__k = - (2.d0 * ( UR0 + Up0_p ) * v_p / R**2)
            case(22) ! --- The correct viscosity (there is only one...)
              Qvisc_UR__p = - (v_R + v/R) * divU + v_Z * (UZ0_R - UR0_Z)
              Qvisc_UR__k = - v_p/R**2 * (UR0_p - Up0 - R*Up0_R)

              Qvisc_UZ__p = - v_Z * divU - v_R * (UZ0_R - UR0_Z)
              Qvisc_UZ__k = + v_p/R**2 * (R*Up0_Z - UZ0_p)

              Qvisc_Up__p = - v_Z/R * (R*Up0_Z - UZ0_p) + (v_R + v/R)/R * (UR0_p - Up0 - R*Up0_R)
              Qvisc_Up__k = - v_p/R * divU
              
              ! --- The toroidal momentum source (same as above replacing all Up0 by -Vt0)
              Qvisc_UR__k = Qvisc_UR__k + v_p/R**2 * (- Vt0 - R*Vt0_R)
              Qvisc_UZ__k = Qvisc_UZ__k - v_p/R**2 * (R*Vt0_Z)
              Qvisc_Up__p = Qvisc_Up__p + v_Z/R * (R*Vt0_Z) - (v_R + v/R)/R * (- Vt0 - R*Vt0_R)
            case(30)
            case default
            end select

            ! --- fourth order diffusion terms
            lap_Vstar = v_R / R + v_RR + v_ZZ

            lap_AR  = AR0_R  / R  + AR0_RR   + AR0_ZZ
            lap_AZ  = AZ0_R  / R  + AZ0_RR   + AZ0_ZZ 
            lap_A3  = A30_R  / R  + A30_RR   + A30_ZZ 
            lap_UR  = UR0_R  / R  + UR0_RR   + UR0_ZZ 
            lap_UZ  = UZ0_R  / R  + UZ0_RR   + UZ0_ZZ 
            lap_Up  = Up0_R  / R  + Up0_RR   + Up0_ZZ 
            lap_rho = rho0_R / R  + rho0_RR  + rho0_ZZ
            lap_T   = T0_R   / R   + T0_RR   + T0_ZZ  

            rhs_p_ij  = 0.d0
            rhs_k_ij  = 0.d0
            VMS__p    = 0.d0
            VMS__k    = 0.d0
            Pvec_prev = 0.d0 ! The time derivative part
            Qvec_p    = 0.d0 ! The rest of the RHS (poloidal part)
            Qvec_k    = 0.d0 ! The rest of the RHS (toroidal part that has phi-derivatives of the test-function)

            !###################################################################################################
            !#  equation 1 (R component induction equation)                                                    #
            !###################################################################################################
            Pvec_prev(var_AR) =   v * delta_g(mp,var_AR,ms,mt)

            Qvec_p(var_AR) = + v * (UZ0 * Bp0 - Up0 * BZ0)                             &
                             + tauIC_ARAZ * v * tau_IC*F0/rho0_corr/BB2 * BR0 * BgradP &
                             + eta_ARAZ * v * (eta_Z * Bp00 - eta_p * BZ0 / R )        &
                             - eta_ARAZ * eta_T * ( - v_Z * Bp00 )                     &
                             - eta_ARAZ_const   * ( - v_Z * Bp00)                      &
                             + eta_ARAZ * eta_T * v * current_source_JR(ms,mt)         &
                             - eta_ARAZ * eta_num * lap_Vstar * lap_AR
            Qvec_k(var_AR) = - eta_ARAZ * eta_T * ( + v_p * BZ0 / R)

            !###################################################################################################
            !#  equation 2 (Z component induction equation)                                                    #
            !###################################################################################################
            Pvec_prev(var_AZ) =   v * delta_g(mp,var_AZ,ms,mt) 

            Qvec_p(var_AZ) = + v * (Up0 * BR0 - UR0 * Bp0)                             &
                             + tauIC_ARAZ * v * tau_IC*F0/rho0_corr/BB2 * BZ0 * BgradP &
                             + eta_ARAZ * v * (eta_p / R * BR0 - eta_R * Bp00)         &
                             - eta_ARAZ * eta_T * ( + v_R * Bp00)                      &
                             - eta_ARAZ_const   * ( + v_R * Bp00)                      &
                             + eta_ARAZ * eta_T * v * current_source_JZ(ms,mt)         &
                             - eta_ARAZ * eta_num * lap_Vstar * lap_AZ
            Qvec_k(var_AZ) = - eta_ARAZ * eta_T * ( - v_p * BR0 / R)

            !###################################################################################################
            !#  equation 3 (PHI component induction equation)                                                  #
            !###################################################################################################
            Pvec_prev(var_A3) =   v * delta_g(mp,var_A3,ms,mt)

            Qvec_p(var_A3) = + R * v * (UR0 * BZ0 - UZ0 * BR0)              &
                             + v * tau_IC*F0/rho0_corr/BB2 * R*Bp0 * BgradP &
                             - eta_T * v_Z * R                * BR0         &
                             + eta_T * ( 2.d0 * v + R * v_R ) * BZ0         &
                             + R * v * (eta_R * BZ0 - eta_Z * BR0)          &
                             + eta_T * v * current_source_Jp(ms,mt)         &
                             - eta_num * lap_Vstar * lap_A3

            !###################################################################################################
            !#  equation 4   (R component momentum equation)                                                   #
            !###################################################################################################
            Pvec_prev(var_UR) =   v * rho0_corr * delta_g(mp,var_UR, ms,mt)
                      
            Qvec_p(var_UR)    = + Qconv_UR              &
                                + (v_R+v/R) * p0        &
                                + JxB_UR__p             &
                                + visco_T * Qvisc_UR__p &
                                - v * PneoR             &
                                - v * particle_source(ms,mt) * UR0 &
                                - visco_num * lap_Vstar * lap_UR &
                                - v * aux_divPIR_perp                                
            Qvec_k(var_UR)    = + JxB_UR__k             &
                                + visco_T * Qvisc_UR__k

            !###################################################################################################
            !#  equation 5   (Z component momentum equation)                                                   #
            !###################################################################################################
            Pvec_prev(var_UZ) =   v * rho0_corr * delta_g(mp,var_UZ, ms,mt)

            Qvec_p(var_UZ)  = + Qconv_UZ              &
                              + v_Z * p0              &
                              + JxB_UZ__p             &
                              + visco_T * Qvisc_UZ__p &
                              - v * PneoZ             &
                              - v * particle_source(ms,mt) * UZ0 &
                              - visco_num * lap_Vstar * lap_UZ &                        
                              - v * aux_divPIZ_perp
            Qvec_k(var_UZ)  = + JxB_UZ__k             &
                              + visco_T * Qvisc_UZ__k
                           
            !###################################################################################################
            !#  equation 6 (Phi component momentum equation)                                                   #
            !###################################################################################################
            Pvec_prev(var_Up) = + v * rho0_corr * BR0 * delta_g(mp,var_UR, ms,mt) &
                                + v * rho0_corr * BZ0 * delta_g(mp,var_UZ, ms,mt) &
                                + v * rho0_corr * BP0 * delta_g(mp,var_Up, ms,mt)

            Qvec_p(var_Up)    = + BR0 * Qconv_UR              &
                                + BZ0 * Qconv_UZ              &
                                + Bp0 * Qconv_Up              &
                                + BR0 * visco_T * Qvisc_UR__p &
                                + BZ0 * visco_T * Qvisc_UZ__p &
                                + Bp0 * visco_T * Qvisc_Up__p &
                                + p0 * BgradVstar__p          &
                                - v * (BR0*PneoR + BZ0*PneoZ) &
                                - v * particle_source(ms,mt) * (BR0*UR0 + BZ0*UZ0 + Bp0*Up0) &
                                - visco_num * lap_Vstar * (BR0 * lap_UR + BZ0 * lap_UZ + Bp0 * lap_Up)
            Qvec_k(var_Up)    = + BR0 * visco_T * Qvisc_UR__k &
                                + BZ0 * visco_T * Qvisc_UZ__k &
                                + Bp0 * visco_T * Qvisc_Up__k &
                                + p0 * BgradVstar__k
            
            !###################################################################################################
            !#  equation 7 (Density equation)                                                                  #
            !###################################################################################################
            Pvec_prev(var_rho) =   v * delta_g(mp,var_rho,ms,mt)

            Qvec_p(var_rho) = - v * ( rho0 * divU + UgradRho )                  &
                              + rho0 * VdiaGradVstar__p                         &
                              - D_prof * gradRho_gradVstar__p                   &
                              - (D_par-D_prof) * BgradVstar__p * BgradRho / BB2 &
                              + v * particle_source(ms,mt)                      &
                              - D_perp_num * lap_Vstar * lap_Rho
            Qvec_k(var_rho) = + rho0 * VdiaGradVstar__k                         &
                              - D_prof * gradRho_gradVstar__k                   &
                              - (D_par-D_prof) * BgradVstar__k * BgradRho / BB2 
            !###################################################################################################
            !#  equation 8 (Pressure equation)                                                                 #
            !###################################################################################################
            Pvec_prev(var_T) =   v * rho0_corr * delta_g(mp,var_T,  ms,mt) &
                               + v * T0_corr   * delta_g(mp,var_rho,ms,mt)

            Qvec_p(var_T) = + v * ( - rho0 * UgradT  -  T0 * UgradRho  -  gamma * p0 * divU ) &
                            + v * heat_source(ms,mt)                                          &
                            + v * (gamma-1.d0) * Qvisc_T                                      &
                            - ZK_prof * gradT_gradVstar__p                                    &
                            - (ZKpar_T-ZK_prof) * BgradVstar__p * BgradT / BB2                &
                            + v * (gamma-1.0d0) * eta_T_ohm * JJ2                             &
                            - ZK_perp_num * lap_Vstar * lap_T
            Qvec_k(var_T) = - ZK_prof * gradT_gradVstar__k                                    &
                            - (ZKpar_T-ZK_prof) * BgradVstar__k * BgradT / BB2

            !###################################################################################################
            !#  VMS STABILISATION                                                                              #
            !###################################################################################################
            SELECT CASE(VmsType)
            CASE(0:10)
              VMS__p(var_AR)   = CvGradAr0 * CvGradVi__p
              VMS__p(var_AZ)   = CvGradAZ0 * CvGradVi__p
              VMS__p(var_A3)   = CvGradA30 * CvGradVi__p
              Qvec_p(var_AR)   = Qvec_p(var_AR) - TG_NUM_Eq*TG_NUM(var_AR) * CoefAdv * VMS__p(var_AR)
              Qvec_p(var_AZ)   = Qvec_p(var_AZ) - TG_NUM_Eq*TG_NUM(var_AZ) * CoefAdv * VMS__p(var_AZ)
              Qvec_p(var_A3)   = Qvec_p(var_A3) - TG_NUM_Eq*TG_NUM(var_A3) * CoefAdv * VMS__p(var_A3)
             
              VMS__k(var_AR)   = CvGradAr0 * CvGradVi__k
              VMS__k(var_AZ)   = CvGradAZ0 * CvGradVi__k
              VMS__k(var_A3)   = CvGradA30 * CvGradVi__k
              Qvec_k(var_AR)   = Qvec_k(var_AR) - TG_NUM_Eq*TG_NUM(var_AR) * CoefAdv * VMS__k(var_AR)
              Qvec_k(var_AZ)   = Qvec_k(var_AZ) - TG_NUM_Eq*TG_NUM(var_AZ) * CoefAdv * VMS__k(var_AZ)
              Qvec_k(var_A3)   = Qvec_k(var_A3) - TG_NUM_Eq*TG_NUM(var_A3) * CoefAdv * VMS__k(var_A3)
             
              VMS__p(var_UR)   = rho0 * ( CvGradUR0 * CvGradVi__p + VbGradUR0 * VbGradVi__p )
              VMS__p(var_UZ)   = rho0 * ( CvGradUZ0 * CvGradVi__p + VbGradUZ0 * VbGradVi__p )
              VMS__p(var_Up)   = rho0 * ( CvGradUp0 * CvGradVi__p + VbGradUp0 * VbGradVi__p )
              Qvec_p(var_UR)   = Qvec_p(var_UR) - TG_NUM_Eq*TG_NUM(var_UR) * CoefAdv * VMS__p(var_UR)
              Qvec_p(var_UZ)   = Qvec_p(var_UZ) - TG_NUM_Eq*TG_NUM(var_UZ) * CoefAdv * VMS__p(var_UZ)
              Qvec_p(var_Up) = Qvec_p(var_Up) - TG_NUM_Eq * TG_NUM(var_Up) * CoefAdv * (  BR0*VMS__p(var_UR) &
                                                                                        + BZ0*VMS__p(var_UZ) &
                                                                                        + Bp0*VMS__p(var_Up) )
             
              VMS__k(var_UR)   = rho0 * ( CvGradUR0 * CvGradVi__k + VbGradUR0 * VbGradVi__k )
              VMS__k(var_UZ)   = rho0 * ( CvGradUZ0 * CvGradVi__k + VbGradUZ0 * VbGradVi__k )
              VMS__k(var_Up)   = rho0 * ( CvGradUp0 * CvGradVi__k + VbGradUp0 * VbGradVi__k )
              Qvec_k(var_UR)   = Qvec_k(var_UR) - TG_NUM_Eq*TG_NUM(var_UR) * CoefAdv * VMS__k(var_UR)
              Qvec_k(var_UZ)   = Qvec_k(var_UZ) - TG_NUM_Eq*TG_NUM(var_UZ) * CoefAdv * VMS__k(var_UZ)
              Qvec_k(var_Up) = Qvec_k(var_Up) - TG_NUM_Eq * TG_NUM(var_Up) * CoefAdv * (  BR0*VMS__k(var_UR) &
                                                                                        + BZ0*VMS__k(var_UZ) &
                                                                                        + Bp0*VMS__k(var_Up) )
             
              VMS__p(var_rho) =         CvGradr0 * CvGradVi__p + VbGradr0 * VbGradVi__p
              VMS__p(var_T)   = rho0 *( CvGradT0 * CvGradVi__p + VbGradT0 * VbGradVi__p )
              Qvec_p(var_rho) = Qvec_p(var_rho) - TG_NUM_Eq * TG_NUM(var_rho) * CoefAdv * VMS__p(var_rho)
              Qvec_p(var_T)   = Qvec_p(var_T)   - TG_NUM_Eq * TG_NUM(var_T)   * CoefAdv * VMS__p(var_T)
             
              VMS__k(var_rho) =         CvGradr0 * CvGradVi__p + VbGradr0 * VbGradVi__k
              VMS__k(var_T)   = rho0 *( CvGradT0 * CvGradVi__p + VbGradT0 * VbGradVi__k )
              Qvec_k(var_rho) = Qvec_k(var_rho) - TG_NUM_Eq * TG_NUM(var_rho) * CoefAdv * VMS__k(var_rho)
              Qvec_k(var_T)   = Qvec_k(var_T)   - TG_NUM_Eq * TG_NUM(var_T)   * CoefAdv * VMS__k(var_T)
            CASE(-1)
              ! --- No VMS
            END SELECT


            !###################################################################################################
            !#  Div(V) STABILISATION                                                                           #
            !###################################################################################################
            Qvec_p(var_UR) = Qvec_p(var_UR) - visco_divV * divU * ( v_R + v / R )
            Qvec_p(var_UZ) = Qvec_p(var_UZ) - visco_divV * divU * v_Z 
            Qvec_p(var_Up) = Qvec_p(var_Up) - visco_divV * divU * BgradVstar__p
            Qvec_k(var_Up) = Qvec_k(var_Up) - visco_divV * divU * BgradVstar__k




            ! --- Fill Up the RHS
            if (use_fft) then
              index_ij =       n_var*n_degrees*(i-1) +       n_var*(j-1) + 1
              do ivar= 1,n_var
                RHS_p_ij(ivar) = tstep * Qvec_p(ivar) + zeta * Pvec_prev(ivar) * tstep / tstep_prev
                RHS_k_ij(ivar) = tstep * Qvec_k(ivar)

                ij = index_ij + (ivar-1)
                RHS_p(mp,ij)   =  RHS_p(mp,ij) + RHS_p_ij(ivar) * wst * R * xjac
                RHS_k(mp,ij)   =  RHS_k(mp,ij) + RHS_k_ij(ivar) * wst * R * xjac
              enddo
            else
              index_ij = n_tor_local*n_var*n_degrees*(i-1) + n_tor_local*n_var*(j-1) + im - n_tor_start +1
              do ivar= 1,n_var
                RHS_p_ij(ivar) = tstep * Qvec_p(ivar) + zeta * Pvec_prev(ivar) * tstep / tstep_prev
                RHS_k_ij(ivar) = tstep * Qvec_k(ivar)

                ij = index_ij + (ivar-1)*n_tor_local
                RHS(ij)        =  RHS(ij) + (RHS_p_ij(ivar) + RHS_k_ij(ivar)) * wst * R * xjac
              enddo
            endif

            do k=1,n_vertex_max

              do l=1,n_degrees

                do in =  n_tor_start, n_tor_end

                  ! --- Basis functions
                  bf    = H(k,l,ms,mt)   * element%size(k,l) * HHZ(in,mp)
                  bf_p  = H(k,l,ms,mt)   * element%size(k,l) * HHZ_p(in,mp)
                  bf_s  = H_s(k,l,ms,mt) * element%size(k,l) * HHZ(in,mp)
                  bf_t  = H_t(k,l,ms,mt) * element%size(k,l) * HHZ(in,mp)
                  bf_R  = (   y_t(ms,mt) * bf_s - y_s(ms,mt) * bf_t ) / xjac
                  bf_Z  = ( - x_t(ms,mt) * bf_s + x_s(ms,mt) * bf_t ) / xjac
                  bf_ss = H_ss(k,l,ms,mt) * element%size(k,l) * HHZ(in,mp)
                  bf_tt = H_tt(k,l,ms,mt) * element%size(k,l) * HHZ(in,mp)
                  bf_st = H_st(k,l,ms,mt) * element%size(k,l) * HHZ(in,mp)
                  bf_RR = (bf_ss * y_t(ms,mt)**2 - 2.d0*bf_st * y_s(ms,mt)*y_t(ms,mt) + bf_tt * y_s(ms,mt)**2  &
                         + bf_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                              &
                         + bf_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) ) )    / xjac**2               & 
                         - xjac_R * (bf_s * y_t(ms,mt) - bf_t * y_s(ms,mt)) / xjac**2
                  bf_ZZ = (bf_ss * x_t(ms,mt)**2 - 2.d0*bf_st * x_s(ms,mt)*x_t(ms,mt) + bf_tt * x_s(ms,mt)**2  &
                         + bf_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                              &
                         + bf_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) ) )       / xjac**2            & 
                         - xjac_Z * (- bf_s * x_t(ms,mt) + bf_t * x_s(ms,mt) ) / xjac**2

                  bf_pp  = H(k,l,ms,mt)   * element%size(k,l) * HHZ_pp(in,mp)
                  bf_sp  = H_s(k,l,ms,mt) * element%size(k,l) * HHZ_p(in,mp)
                  bf_tp  = H_t(k,l,ms,mt) * element%size(k,l) * HHZ_p(in,mp)
                  bf_Rp = (   y_t(ms,mt) * bf_sp - y_s(ms,mt) * bf_tp ) / xjac
                  bf_Zp = ( - x_t(ms,mt) * bf_sp + x_s(ms,mt) * bf_tp ) / xjac
                  bf_RZ = (- bf_ss * y_t(ms,mt)*x_t(ms,mt) - bf_tt * x_s(ms,mt)*y_s(ms,mt) &
                         + bf_st * (y_s(ms,mt)*x_t(ms,mt)  + y_t(ms,mt)*x_s(ms,mt)  ) &
                         - bf_s  * (x_st(ms,mt)*y_t(ms,mt) - x_tt(ms,mt)*y_s(ms,mt) ) &
                         - bf_t * (x_st(ms,mt)*y_s(ms,mt)  - x_ss(ms,mt)*y_t(ms,mt) ) ) / xjac**2             &
                         - xjac_R * (- bf_s * x_t(ms,mt) + bf_t * x_s(ms,mt) )   / xjac**2

                  UR    = bf    ;  UZ    = bf    ;  Up    = bf
                  UR_R  = bf_R  ;  UZ_R  = bf_R  ;  Up_R  = bf_R
                  UR_Z  = bf_Z  ;  UZ_Z  = bf_Z  ;  Up_Z  = bf_Z
                  UR_p  = bf_p  ;  UZ_p  = bf_p  ;  Up_p  = bf_p
                  UR_s  = bf_s  ;  UZ_s  = bf_s  ;  Up_s  = bf_s
                  UR_t  = bf_t  ;  UZ_t  = bf_t  ;  Up_t  = bf_t

                  AR    = bf    ;  AZ    = bf    ;  A3    = bf    ; T    = bf    ; rho    = bf
                  AR_R  = bf_R  ;  AZ_R  = bf_R  ;  A3_R  = bf_R  ; T_R  = bf_R  ; rho_R  = bf_R 
                  AR_Z  = bf_Z  ;  AZ_Z  = bf_Z  ;  A3_Z  = bf_Z  ; T_Z  = bf_Z  ; rho_Z  = bf_Z
                  AR_p  = bf_p  ;  AZ_p  = bf_p  ;  A3_p  = bf_p  ; T_p  = bf_p  ; rho_p  = bf_p
                  AR_s  = bf_s  ;  AZ_s  = bf_s  ;  A3_s  = bf_s  ; T_s  = bf_s  ; rho_s  = bf_s 
                  AR_t  = bf_t  ;  AZ_t  = bf_t  ;  A3_t  = bf_t  ; T_t  = bf_T  ; rho_t  = bf_T 

                  AR_RR = bf_RR ; AR_RZ = bf_RZ ; AR_ZZ = bf_ZZ
                  AZ_RR = bf_RR ; AZ_RZ = bf_RZ ; AZ_ZZ = bf_ZZ
                  A3_RR = bf_RR ; A3_RZ = bf_RZ ; A3_ZZ = bf_ZZ

                  AR_Rp = bf_Rp ; AR_Zp  = bf_Zp ; AR_pp = bf_pp
                  AZ_Rp = bf_Rp ; AZ_Zp  = bf_Zp ; AZ_pp = bf_pp
                  A3_Rp = bf_Rp ; A3_Zp  = bf_Zp ; A3_pp = bf_pp

                  ! --- Linearised quantities
                  BR0_AR    =   0.d0     ; BR0_AZ__n = - AZ_p / R ; BR0_A3 =   A3_Z / R
                  BZ0_AR__n =   AR_p / R ; BZ0_AZ    =   0.d0     ; BZ0_A3 = - A3_R / R
                  Bp0_AR    = - AR_Z     ; Bp0_AZ    =   AZ_R     ; Bp0_A3 =   0.d0

                  BB2_AR__p = 2.d0*(BR0_AR    * BR0                   + Bp0_AR * Bp0 )
                  BB2_AR__n = 2.d0*(                + BZ0_AR__n * BZ0                )
                  BB2_AZ__p = 2.d0*(                + BZ0_AZ    * BZ0 + Bp0_AZ * Bp0 )
                  BB2_AZ__n = 2.d0*(BR0_AZ__n * BR0                                  )
                  BB2_A3    = 2.d0*(BR0_A3    * BR0 + BZ0_A3    * BZ0 + Bp0_A3 * Bp0 )

                  BgradT_AR__p = BR0_AR    * T0_R                    + Bp0_AR * T0_p / R
                  BgradT_AR__n =                  + BZ0_AR__n * T0_Z
                  BgradT_AZ__p =                  + BZ0_AZ    * T0_Z + Bp0_AZ * T0_p / R
                  BgradT_AZ__n = BR0_AZ__n * T0_R
                  BgradT_A3    = BR0_A3    * T0_R + BZ0_A3    * T0_Z + Bp0_A3 * T0_p / R
                  BgradT_T__p  = BR0 * T_R + BZ0 * T_Z
                  BgradT_T__n  = Bp0 * T_p / R

                  BgradRho_AR__p  = BR0_AR    * rho0_R                      + Bp0_AR * rho0_p / R
                  BgradRho_AR__n  =                    + BZ0_AR__n * rho0_Z
                  BgradRho_AZ__p  =                    + BZ0_AZ    * rho0_Z + Bp0_AZ * rho0_p / R
                  BgradRho_AZ__n  = BR0_AZ__n * rho0_R
                  BgradRho_A3     = BR0_A3    * rho0_R + BZ0_A3    * rho0_Z + Bp0_A3 * rho0_p / R
                  BgradRho_rho__p = BR0 * rho_R + BZ0 * rho_Z
                  BgradRho_rho__n = Bp0 * rho_p / R

                  BgradP_AR__p  = rho0*BgradT_AR__p + T0*BgradRho_AR__p
                  BgradP_AR__n  = rho0*BgradT_AR__n + T0*BgradRho_AR__n
                  BgradP_AZ__p  = rho0*BgradT_AZ__p + T0*BgradRho_AZ__p
                  BgradP_AZ__n  = rho0*BgradT_AZ__n + T0*BgradRho_AZ__n
                  BgradP_A3     = rho0*BgradT_A3    + T0*BgradRho_A3   
                  BgradP_rho__p = rho*BgradT        + T0*BgradRho_rho__p
                  BgradP_rho__n =                   + T0*BgradRho_rho__n
                  BgradP_T__p   = rho0*BgradT_T__p  + T *BgradRho
                  BgradP_T__n   = rho0*BgradT_T__n

                  gradBF_gradVstar__p  = bf_R * v_R + bf_Z * v_Z
                  gradBF_gradVstar__kn = (bf_p/R) * (v_p/R)
                  BgradBF__p = BR0 * bf_R  + BZ0 * bf_Z
                  BgradBF__n = Bp0 * bf_p  / R
                  UgradBF__p = UR0 * bf_R  + UZ0 * bf_Z
                  UgradBF__n = Up0 * bf_p  / R

                  UgradUR_UR__p = UR  * UR0_R + UR0 * UR_R + UZ0 * UR_Z
                  UgradUR_UR__n = Up0 * UR_p  / R
                  UgradUR_UZ    = UZ  * UR0_Z
                  UgradUR_Up    = Up  * UR0_p / R

                  UgradUZ_UR    = UR  * UZ0_R
                  UgradUZ_UZ__p = UR0 * UZ_R + UZ * UZ0_Z + UZ0 * UZ_Z
                  UgradUZ_UZ__n = Up0 * UZ_p  / R
                  UgradUZ_Up    = Up  * UZ0_p / R

                  UgradUp_UR    = UR  * Up0_R
                  UgradUp_UZ    = UZ  * Up0_Z
                  UgradUp_Up__p = UR0 * Up_R + UZ0 * Up_Z + Up  * Up0_p / R
                  UgradUp_Up__n = Up0 * Up_p  / R

                  UgradT_UR   = UR * T0_R
                  UgradT_UZ   = UZ * T0_Z
                  UgradT_Up   = Up * T0_p / R
                  UgradT_T__p = UR0 * T_R + UZ0 * T_Z
                  UgradT_T__n = Up0 * T_p / R

                  UgradRho_UR     = UR * rho0_R
                  UgradRho_UZ     = UZ * rho0_Z
                  UgradRho_Up     = Up * rho0_p / R
                  UgradRho_rho__p = UR0 * rho_R + UZ0 * rho_Z
                  UgradRho_rho__n = Up0 * rho_p / R

                  divU_UR    = UR_R + UR / R
                  divU_UZ    = UZ_Z
                  divU_Up__n = Up_p / R

                  divRhoU_UR    = rho0 * divU_UR + UgradRho_UR
                  divRhoU_UZ    = rho0 * divU_UZ + UgradRho_UZ
                  divRhoU_Up__p = UgradRho_Up
                  divRhoU_Up__n = rho0 * divU_Up__n
                  divRhoU_rho__p = rho * divU + UgradRho_rho__p
                  divRhoU_rho__n = UgradRho_rho__n

                  VdiaR0_AR__p  = + tau_IC*F0 / (R * rho0_corr * BB2   ) * (                 - R*Bp0_AR*p0_Z) &
                                  - tau_IC*F0 / (R * rho0_corr * BB2**2) * (  BZ0      *p0_p - R*Bp0   *p0_Z) * BB2_AR__p
                  VdiaR0_AR__n  = + tau_IC*F0 / (R * rho0_corr * BB2   ) * (  BZ0_AR__n*p0_p                ) &
                                  - tau_IC*F0 / (R * rho0_corr * BB2**2) * (  BZ0      *p0_p - R*Bp0   *p0_Z) * BB2_AR__n
                  VdiaR0_AZ__p  = + tau_IC*F0 / (R * rho0_corr * BB2   ) * (  BZ0_AZ   *p0_p - R*Bp0_AZ*p0_Z) &
                                  - tau_IC*F0 / (R * rho0_corr * BB2**2) * (  BZ0      *p0_p - R*Bp0   *p0_Z) * BB2_AZ__p
                  VdiaR0_AZ__n  = - tau_IC*F0 / (R * rho0_corr * BB2**2) * (  BZ0      *p0_p - R*Bp0   *p0_Z) * BB2_AZ__n
                  VdiaR0_A3__p  = + tau_IC*F0 / (R * rho0_corr * BB2   ) * (  BZ0_A3   *p0_p - R*Bp0_A3*p0_Z) &
                                  - tau_IC*F0 / (R * rho0_corr * BB2**2) * (  BZ0      *p0_p - R*Bp0   *p0_Z) * BB2_A3
                  VdiaR0_A3__n  = 0.d0
                  VdiaR0_rho__p = + tau_IC*F0 / (R * rho0_corr * BB2   ) * (  BZ0*(rho*T0_p) - R*Bp0*(rho*T0_Z+rho_Z*T0)) &
                                  - tau_IC*F0 / (R * rho0_corr**2 * BB2) * (  BZ0*p0_p       - R*Bp0*p0_Z               ) * rho
                  VdiaR0_rho__n = + tau_IC*F0 / (R * rho0_corr * BB2   ) * (  BZ0*(rho_p*T0)                            )
                  VdiaR0_T__p   = + tau_IC*F0 / (R * rho0_corr * BB2   ) * (  BZ0*(rho0_p*T) - R*Bp0*(rho0*T_Z+rho0_Z*T))
                  VdiaR0_T__n   = + tau_IC*F0 / (R * rho0_corr * BB2   ) * (  BZ0*(rho0*T_p)                            )

                  VdiaZ0_AR__p  = + tau_IC*F0 / (R * rho0_corr * BB2   ) * (R*BP0_AR*p0_R -   BR0_AR   *p0_p) &
                                  - tau_IC*F0 / (R * rho0_corr * BB2**2) * (R*BP0   *p0_R -   BR0      *p0_p) * BB2_AR__p
                  VdiaZ0_AR__n  = - tau_IC*F0 / (R * rho0_corr * BB2**2) * (R*BP0   *p0_R -   BR0      *p0_p) * BB2_AR__n
                  VdiaZ0_AZ__p  = + tau_IC*F0 / (R * rho0_corr * BB2   ) * (R*BP0_AZ*p0_R                   ) &
                                  - tau_IC*F0 / (R * rho0_corr * BB2**2) * (R*BP0   *p0_R -   BR0      *p0_p) * BB2_AZ__p
                  VdiaZ0_AZ__n  = + tau_IC*F0 / (R * rho0_corr * BB2   ) * (              -   BR0_AZ__n*p0_p) &
                                  - tau_IC*F0 / (R * rho0_corr * BB2**2) * (R*BP0   *p0_R -   BR0      *p0_p) * BB2_AZ__n
                  VdiaZ0_A3__p  = + tau_IC*F0 / (R * rho0_corr * BB2   ) * (R*BP0_A3*p0_R -   BR0_A3   *p0_p) &
                                  - tau_IC*F0 / (R * rho0_corr * BB2**2) * (R*BP0   *p0_R -   BR0      *p0_p) * BB2_A3
                  VdiaZ0_A3__n  = 0.d0
                  VdiaZ0_rho__p = + tau_IC*F0 / (R * rho0_corr * BB2   ) * (R*BP0*(rho*T0_R+rho_R*T0) -   BR0*(rho*T0_p)) &
                                  - tau_IC*F0 / (R * rho0_corr**2 * BB2) * (R*BP0*p0_R                -   BR0*p0_p      ) * rho
                  VdiaZ0_rho__n = + tau_IC*F0 / (R * rho0_corr * BB2   ) * (                          -   BR0*(rho_p*T0))
                  VdiaZ0_T__p   = + tau_IC*F0 / (R * rho0_corr * BB2   ) * (R*BP0*(T*rho0_R+T_R*rho0) -   BR0*(T*rho0_p))
                  VdiaZ0_T__n   = + tau_IC*F0 / (R * rho0_corr * BB2   ) * (                          -   BR0*(T_p*rho0))

                  VdiaP0_AR__p  = + tau_IC*F0 / (    rho0_corr * BB2   ) * (  BR0_AR   *p0_Z                   ) &
                                  - tau_IC*F0 / (    rho0_corr * BB2**2) * (  BR0      *p0_Z -   BZ0      *p0_R) * BB2_AR__p
                  VdiaP0_AR__n  = + tau_IC*F0 / (    rho0_corr * BB2   ) * (                 -   BZ0_AR__n*p0_R) &
                                  - tau_IC*F0 / (    rho0_corr * BB2**2) * (  BR0      *p0_Z -   BZ0      *p0_R) * BB2_AR__n
                  VdiaP0_AZ__p  = + tau_IC*F0 / (    rho0_corr * BB2   ) * (                 -   BZ0_AZ   *p0_R) &
                                  - tau_IC*F0 / (    rho0_corr * BB2**2) * (  BR0      *p0_Z -   BZ0      *p0_R) * BB2_AZ__p
                  VdiaP0_AZ__n  = + tau_IC*F0 / (    rho0_corr * BB2   ) * (  BR0_AZ__n*p0_Z                   ) &
                                  - tau_IC*F0 / (    rho0_corr * BB2**2) * (  BR0      *p0_Z -   BZ0      *p0_R) * BB2_AZ__n
                  VdiaP0_A3__p  = + tau_IC*F0 / (    rho0_corr * BB2   ) * (  BR0_A3   *p0_Z -   BZ0_A3   *p0_R) &
                                  - tau_IC*F0 / (    rho0_corr * BB2**2) * (  BR0      *p0_Z -   BZ0      *p0_R) * BB2_A3
                  VdiaP0_A3__n  = 0.d0
                  VdiaP0_rho__p = + tau_IC*F0 / (    rho0_corr * BB2   ) * (  BR0*(rho*T0_Z+rho_Z*T0) -   BZ0*(rho*T0_R+rho_R*T0)) &
                                  - tau_IC*F0 / (    rho0_corr**2 * BB2) * (  BR0*p0_Z                -   BZ0*p0_R               ) * rho
                  VdiaP0_rho__n = 0.d0
                  VdiaP0_T__p   = + tau_IC*F0 / (    rho0_corr * BB2   ) * (  BR0*(rho0*T_Z+rho0_Z*T) -   BZ0*(rho0*T_R+rho0_R*T))
                  VdiaP0_T__n   = 0.d0

                  VdiaGradUR_AR__p  = VdiaR0_AR__p  * UR0_R + VdiaZ0_AR__p  * UR0_Z + VdiaP0_AR__p  * UR0_p / R
                  VdiaGradUR_AR__n  = VdiaR0_AR__n  * UR0_R + VdiaZ0_AR__n  * UR0_Z + VdiaP0_AR__n  * UR0_p / R
                  VdiaGradUR_AZ__p  = VdiaR0_AZ__p  * UR0_R + VdiaZ0_AZ__p  * UR0_Z + VdiaP0_AZ__p  * UR0_p / R
                  VdiaGradUR_AZ__n  = VdiaR0_AZ__n  * UR0_R + VdiaZ0_AZ__n  * UR0_Z + VdiaP0_AZ__n  * UR0_p / R
                  VdiaGradUR_A3__p  = VdiaR0_A3__p  * UR0_R + VdiaZ0_A3__p  * UR0_Z + VdiaP0_A3__p  * UR0_p / R
                  VdiaGradUR_A3__n  = VdiaR0_A3__n  * UR0_R + VdiaZ0_A3__n  * UR0_Z + VdiaP0_A3__n  * UR0_p / R
                  VdiaGradUR_rho__p = VdiaR0_rho__p * UR0_R + VdiaZ0_rho__p * UR0_Z + VdiaP0_rho__p * UR0_p / R
                  VdiaGradUR_rho__n = VdiaR0_rho__n * UR0_R + VdiaZ0_rho__n * UR0_Z + VdiaP0_rho__n * UR0_p / R
                  VdiaGradUR_T__p   = VdiaR0_T__p   * UR0_R + VdiaZ0_T__p   * UR0_Z + VdiaP0_T__p   * UR0_p / R
                  VdiaGradUR_T__n   = VdiaR0_T__n   * UR0_R + VdiaZ0_T__n   * UR0_Z + VdiaP0_T__n   * UR0_p / R
                  VdiaGradUR_UR__p  = VdiaR0        * UR_R  + VdiaZ0        * UR_Z
                  VdiaGradUR_UR__n  =                                               + VdiaP0        * UR_p  / R

                  VdiaGradUZ_AR__p  = VdiaR0_AR__p  * UZ0_R + VdiaZ0_AR__p  * UZ0_Z + VdiaP0_AR__p  * UZ0_p / R
                  VdiaGradUZ_AR__n  = VdiaR0_AR__n  * UZ0_R + VdiaZ0_AR__n  * UZ0_Z + VdiaP0_AR__n  * UZ0_p / R
                  VdiaGradUZ_AZ__p  = VdiaR0_AZ__p  * UZ0_R + VdiaZ0_AZ__p  * UZ0_Z + VdiaP0_AZ__p  * UZ0_p / R
                  VdiaGradUZ_AZ__n  = VdiaR0_AZ__n  * UZ0_R + VdiaZ0_AZ__n  * UZ0_Z + VdiaP0_AZ__n  * UZ0_p / R
                  VdiaGradUZ_A3__p  = VdiaR0_A3__p  * UZ0_R + VdiaZ0_A3__p  * UZ0_Z + VdiaP0_A3__p  * UZ0_p / R
                  VdiaGradUZ_A3__n  = VdiaR0_A3__n  * UZ0_R + VdiaZ0_A3__n  * UZ0_Z + VdiaP0_A3__n  * UZ0_p / R
                  VdiaGradUZ_rho__p = VdiaR0_rho__p * UZ0_R + VdiaZ0_rho__p * UZ0_Z + VdiaP0_rho__p * UZ0_p / R
                  VdiaGradUZ_rho__n = VdiaR0_rho__n * UZ0_R + VdiaZ0_rho__n * UZ0_Z + VdiaP0_rho__n * UZ0_p / R
                  VdiaGradUZ_T__p   = VdiaR0_T__p   * UZ0_R + VdiaZ0_T__p   * UZ0_Z + VdiaP0_T__p   * UZ0_p / R
                  VdiaGradUZ_T__n   = VdiaR0_T__n   * UZ0_R + VdiaZ0_T__n   * UZ0_Z + VdiaP0_T__n   * UZ0_p / R
                  VdiaGradUZ_UZ__p  = VdiaR0        * UZ_R  + VdiaZ0        * UZ_Z
                  VdiaGradUZ_UZ__n  =                                               + VdiaP0        * UZ_p  / R

                  VdiaGradUp_AR__p  = VdiaR0_AR__p  * Up0_R + VdiaZ0_AR__p  * Up0_Z + VdiaP0_AR__p  * Up0_p / R
                  VdiaGradUp_AR__n  = VdiaR0_AR__n  * Up0_R + VdiaZ0_AR__n  * Up0_Z + VdiaP0_AR__n  * Up0_p / R
                  VdiaGradUp_AZ__p  = VdiaR0_AZ__p  * Up0_R + VdiaZ0_AZ__p  * Up0_Z + VdiaP0_AZ__p  * Up0_p / R
                  VdiaGradUp_AZ__n  = VdiaR0_AZ__n  * Up0_R + VdiaZ0_AZ__n  * Up0_Z + VdiaP0_AZ__n  * Up0_p / R
                  VdiaGradUp_A3__p  = VdiaR0_A3__p  * Up0_R + VdiaZ0_A3__p  * Up0_Z + VdiaP0_A3__p  * Up0_p / R
                  VdiaGradUp_A3__n  = VdiaR0_A3__n  * Up0_R + VdiaZ0_A3__n  * Up0_Z + VdiaP0_A3__n  * Up0_p / R
                  VdiaGradUp_rho__p = VdiaR0_rho__p * Up0_R + VdiaZ0_rho__p * Up0_Z + VdiaP0_rho__p * Up0_p / R
                  VdiaGradUp_rho__n = VdiaR0_rho__n * Up0_R + VdiaZ0_rho__n * Up0_Z + VdiaP0_rho__n * Up0_p / R
                  VdiaGradUp_T__p   = VdiaR0_T__p   * Up0_R + VdiaZ0_T__p   * Up0_Z + VdiaP0_T__p   * Up0_p / R
                  VdiaGradUp_T__n   = VdiaR0_T__n   * Up0_R + VdiaZ0_T__n   * Up0_Z + VdiaP0_T__n   * Up0_p / R
                  VdiaGradUp_Up__p  = VdiaR0        * Up_R  + VdiaZ0        * Up_Z
                  VdiaGradUp_Up__n  =                                               + VdiaP0        * Up_p  / R

                  !if (NEO) then ! can't remember, I think only the RHS is needed because it's small. maybe all this can be removed...
                  if (.false.) then
                    if (Btht.le. 1.d-10) then
                      Btht_AR__p = 0.d0
                      Btht_AR__n = 0.d0
                      Btht_AZ__p = 0.d0
                      Btht_AZ__n = 0.d0
                      Btht_A3    = 0.d0
                    else
                      Btht_AR__p = 0.d0
                      Btht_AR__n = 0.5*(BR0**2 + BZ0**2)**(-0.5) * 2.0*(BZ0*BZ0_AR__n)
                      Btht_AZ__p = 0.d0
                      Btht_AZ__n = 0.5*(BR0**2 + BZ0**2)**(-0.5) * 2.0*(BR0*BR0_AZ__n)
                      Btht_A3    = 0.5*(BR0**2 + BZ0**2)**(-0.5) * 2.0*(BR0*BR0_A3 + BZ0*BZ0_A3)
                    endif
                    
                    Vtht_AR__p  = + ( BR0      *(VdiaR0_AR__p) + BZ0      *(VdiaZ0_AR__p) ) / Btht &
                                  - ( BR0      *(UR0+VdiaR0)   + BZ0      *(UZ0+VdiaZ0)   ) / Btht**2 * Btht_AR__p
                    Vtht_AR__n  = + ( BR0      *(VdiaR0_AR__n) + BZ0      *(VdiaZ0_AR__n) ) / Btht &
                                  + (                          + BZ0_AR__n*(UZ0+VdiaZ0)   ) / Btht &
                                  - ( BR0      *(UR0+VdiaR0)   + BZ0      *(UZ0+VdiaZ0)   ) / Btht**2 * Btht_AR__n
                    Vtht_AZ__p  = + ( BR0      *(VdiaR0_AZ__p) + BZ0      *(VdiaZ0_AZ__p) ) / Btht &
                                  - ( BR0      *(UR0+VdiaR0)   + BZ0      *(UZ0+VdiaZ0)   ) / Btht**2 * Btht_AZ__p
                    Vtht_AZ__n  = + ( BR0      *(VdiaR0_AZ__n) + BZ0      *(VdiaZ0_AZ__n) ) / Btht &
                                  + ( BR0_AZ__n*(UR0+VdiaR0)                              ) / Btht &
                                  - ( BR0      *(UR0+VdiaR0)   + BZ0      *(UZ0+VdiaZ0)   ) / Btht**2 * Btht_AZ__n
                    Vtht_A3__p  = + ( BR0      *(VdiaR0_A3__p) + BZ0      *(VdiaZ0_A3__p) ) / Btht &
                                  + ( BR0_A3   *(UR0+VdiaR0)   + BZ0_A3   *(UZ0+VdiaZ0)   ) / Btht &
                                  - ( BR0      *(UR0+VdiaR0)   + BZ0      *(UZ0+VdiaZ0)   ) / Btht**2 * Btht_A3
                    Vtht_A3__n  = + ( BR0      *(VdiaR0_A3__n) + BZ0      *(VdiaZ0_A3__n) ) / Btht
                    Vtht_rho__p = + ( BR0*(UR0+VdiaR0_rho__p) + BZ0*(UZ0+VdiaZ0_rho__p) ) / Btht
                    Vtht_rho__n = + ( BR0*(UR0+VdiaR0_rho__n) + BZ0*(UZ0+VdiaZ0_rho__n) ) / Btht
                    Vtht_T__p   = + ( BR0*(UR0+VdiaR0_T__p  ) + BZ0*(UZ0+VdiaZ0_T__p  ) ) / Btht
                    Vtht_T__n   = + ( BR0*(UR0+VdiaR0_T__n  ) + BZ0*(UZ0+VdiaZ0_T__n  ) ) / Btht
                    Vtht_UR     = + ( BR0*UR ) / Btht
                    Vtht_UZ     = + ( BZ0*UZ ) / Btht
                    
                    Vneo_AR__p = -aki_neo_prof(ms,mt) * tau_IC*F0 / BB2    / Btht    * ( BR0*Bp0_AR*T0_Z    - BZ0*Bp0_AR*T0_R   ) &
                                 +aki_neo_prof(ms,mt) * tau_IC*F0 / BB2**2 / Btht    * ( BR0*Bp0*T0_Z       - BZ0*Bp0*T0_R      ) * BB2_AR__p &
                                 +aki_neo_prof(ms,mt) * tau_IC*F0 / BB2    / Btht**2 * ( BR0*Bp0*T0_Z       - BZ0*Bp0*T0_R      ) * Btht_AR__p
                    Vneo_AR__n = -aki_neo_prof(ms,mt) * tau_IC*F0 / BB2    / Btht    * (                    - BZ0_AR__n*Bp0*T0_R) &
                                 +aki_neo_prof(ms,mt) * tau_IC*F0 / BB2**2 / Btht    * ( BR0*Bp0*T0_Z       - BZ0*Bp0*T0_R      ) * BB2_AR__n &
                                 +aki_neo_prof(ms,mt) * tau_IC*F0 / BB2    / Btht**2 * ( BR0*Bp0*T0_Z       - BZ0*Bp0*T0_R      ) * Btht_AR__n
                    Vneo_AZ__p = -aki_neo_prof(ms,mt) * tau_IC*F0 / BB2    / Btht    * ( BR0*Bp0_AZ*T0_Z    - BZ0*Bp0_AZ*T0_R   ) &
                                 +aki_neo_prof(ms,mt) * tau_IC*F0 / BB2**2 / Btht    * ( BR0*Bp0*T0_Z       - BZ0*Bp0*T0_R      ) * BB2_AZ__p &
                                 +aki_neo_prof(ms,mt) * tau_IC*F0 / BB2    / Btht**2 * ( BR0*Bp0*T0_Z       - BZ0*Bp0*T0_R      ) * Btht_AZ__p
                    Vneo_AZ__n = -aki_neo_prof(ms,mt) * tau_IC*F0 / BB2    / Btht    * ( BR0_AZ__n*Bp0*T0_Z                     ) &
                                 +aki_neo_prof(ms,mt) * tau_IC*F0 / BB2**2 / Btht    * ( BR0*Bp0*T0_Z       - BZ0*Bp0*T0_R      ) * BB2_AZ__n &
                                 +aki_neo_prof(ms,mt) * tau_IC*F0 / BB2    / Btht**2 * ( BR0*Bp0*T0_Z       - BZ0*Bp0*T0_R      ) * Btht_AZ__n
                    Vneo_A3__p = -aki_neo_prof(ms,mt) * tau_IC*F0 / BB2    / Btht    * ( BR0_A3*Bp0*T0_Z    - BZ0_A3*Bp0*T0_R   ) &
                                 +aki_neo_prof(ms,mt) * tau_IC*F0 / BB2**2 / Btht    * ( BR0*Bp0*T0_Z       - BZ0*Bp0*T0_R      ) * BB2_A3 &
                                 +aki_neo_prof(ms,mt) * tau_IC*F0 / BB2    / Btht**2 * ( BR0*Bp0*T0_Z       - BZ0*Bp0*T0_R      ) * Btht_A3
                    Vneo_A3__n = 0.d0
                    Vneo_T__p  = -aki_neo_prof(ms,mt) * tau_IC*F0 / BB2    / Btht    * ( BR0*Bp0*T_Z        - BZ0*Bp0*T_R       )
                    Vneo_T__n  = 0.d0
                    
                    PneoR_AR__p  = +amu_neo_prof(ms,mt) * rho0 * BB2_AR__p/Btht**3 * BR0       * (Vtht        - Vneo      ) &
                                   +amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**3 * BR0       * (Vtht_AR__p  - Vneo_AR__p) &
                                   -amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**4 * BR0       * (Vtht        - Vneo      ) * 3.0 * Btht_AR__p
                    PneoR_AR__n  = +amu_neo_prof(ms,mt) * rho0 * BB2_AR__n/Btht**3 * BR0       * (Vtht        - Vneo      ) &
                                   +amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**3 * BR0       * (Vtht_AR__n  - Vneo_AR__n) &
                                   -amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**4 * BR0       * (Vtht        - Vneo      ) * 3.0 * Btht_AR__n
                    PneoR_AZ__p  = +amu_neo_prof(ms,mt) * rho0 * BB2_AZ__p/Btht**3 * BR0       * (Vtht        - Vneo      ) &
                                   +amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**3 * BR0       * (Vtht_AZ__p  - Vneo_AZ__p) &
                                   -amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**4 * BR0       * (Vtht        - Vneo      ) * 3.0 * Btht_AZ__p
                    PneoR_AZ__n  = +amu_neo_prof(ms,mt) * rho0 * BB2_AZ__n/Btht**3 * BR0       * (Vtht        - Vneo      ) &
                                   +amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**3 * BR0_AZ__n * (Vtht        - Vneo      ) &
                                   +amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**3 * BR0       * (Vtht_AZ__n  - Vneo_AZ__n) &
                                   -amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**4 * BR0       * (Vtht        - Vneo      ) * 3.0 * Btht_AZ__n
                    PneoR_A3__p  = +amu_neo_prof(ms,mt) * rho0 * BB2_A3   /Btht**3 * BR0       * (Vtht        - Vneo      ) &
                                   +amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**3 * BR0_A3    * (Vtht        - Vneo      ) &
                                   +amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**3 * BR0       * (Vtht_A3__p  - Vneo_A3__p) &
                                   -amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**4 * BR0       * (Vtht        - Vneo      ) * 3.0 * Btht_A3
                    PneoR_A3__n  = +amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**3 * BR0       * (Vtht_A3__n  - Vneo_A3__n)
                    PneoR_rho__p = +amu_neo_prof(ms,mt) * rho  * BB2      /Btht**3 * BR0       * (Vtht        - Vneo      ) &
                                   +amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**3 * BR0       * (Vtht_rho__p             )
                    PneoR_rho__n = +amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**3 * BR0       * (Vtht_rho__n             )
                    PneoR_T__p   = +amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**3 * BR0       * (Vtht_T__p   - Vneo_T__p )
                    PneoR_T__n   = +amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**3 * BR0       * (Vtht_T__n   - Vneo_T__n )
                    PneoR_UR     = +amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**3 * BR0       * (Vtht_UR                 )
                    PneoR_UZ     = +amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**3 * BR0       * (Vtht_UZ                 )
                    
                    PneoZ_AR__p  = +amu_neo_prof(ms,mt) * rho0 * BB2_AR__p/Btht**3 * BZ0       * (Vtht        - Vneo      ) &
                                   +amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**3 * BZ0       * (Vtht_AR__p  - Vneo_AR__p) &
                                   -amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**4 * BZ0       * (Vtht        - Vneo      ) * 3.0 * Btht_AR__p
                    PneoZ_AR__n  = +amu_neo_prof(ms,mt) * rho0 * BB2_AR__n/Btht**3 * BZ0       * (Vtht        - Vneo      ) &
                                   +amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**3 * BZ0_AR__n * (Vtht        - Vneo      ) &
                                   +amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**3 * BZ0       * (Vtht_AR__n  - Vneo_AR__n) &
                                   -amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**4 * BZ0       * (Vtht        - Vneo      ) * 3.0 * Btht_AR__n
                    PneoZ_AZ__p  = +amu_neo_prof(ms,mt) * rho0 * BB2_AZ__p/Btht**3 * BZ0       * (Vtht        - Vneo      ) &
                                   +amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**3 * BZ0       * (Vtht_AZ__p  - Vneo_AZ__p) &
                                   -amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**4 * BZ0       * (Vtht        - Vneo      ) * 3.0 * Btht_AZ__p
                    PneoZ_AZ__n  = +amu_neo_prof(ms,mt) * rho0 * BB2_AZ__n/Btht**3 * BZ0       * (Vtht        - Vneo      ) &
                                   +amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**3 * BZ0       * (Vtht_AZ__n  - Vneo_AZ__n) &
                                   -amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**4 * BZ0       * (Vtht        - Vneo      ) * 3.0 * Btht_AZ__n
                    PneoZ_A3__p  = +amu_neo_prof(ms,mt) * rho0 * BB2_A3   /Btht**3 * BZ0       * (Vtht        - Vneo      ) &
                                   +amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**3 * BZ0_A3    * (Vtht        - Vneo      ) &
                                   +amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**3 * BZ0       * (Vtht_A3__p  - Vneo_A3__p) &
                                   -amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**4 * BZ0       * (Vtht        - Vneo      ) * 3.0 * Btht_A3
                    PneoZ_A3__n  = +amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**3 * BZ0       * (Vtht_A3__n  - Vneo_A3__n)
                    PneoZ_rho__p = +amu_neo_prof(ms,mt) * rho  * BB2      /Btht**3 * BZ0       * (Vtht        - Vneo      ) &
                                   +amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**3 * BZ0       * (Vtht_rho__p             )
                    PneoZ_rho__n = +amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**3 * BZ0       * (Vtht_rho__n             )
                    PneoZ_T__p   = +amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**3 * BZ0       * (Vtht_T__p   - Vneo_T__p )
                    PneoZ_T__n   = +amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**3 * BZ0       * (Vtht_T__n   - Vneo_T__n )
                    PneoZ_UR     = +amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**3 * BZ0       * (Vtht_UR                 )
                    PneoZ_UZ     = +amu_neo_prof(ms,mt) * rho0 * BB2      /Btht**3 * BZ0       * (Vtht_UZ                 )
                  else
                    PneoR_AR__p  = 0.d0 ; PneoZ_AR__p  = 0.d0
                    PneoR_AR__n  = 0.d0 ; PneoZ_AR__n  = 0.d0
                    PneoR_AZ__p  = 0.d0 ; PneoZ_AZ__p  = 0.d0
                    PneoR_AZ__n  = 0.d0 ; PneoZ_AZ__n  = 0.d0
                    PneoR_A3__p  = 0.d0 ; PneoZ_A3__p  = 0.d0
                    PneoR_A3__n  = 0.d0 ; PneoZ_A3__n  = 0.d0
                    PneoR_rho__p = 0.d0 ; PneoZ_rho__p = 0.d0
                    PneoR_rho__n = 0.d0 ; PneoZ_rho__n = 0.d0
                    PneoR_T__p   = 0.d0 ; PneoZ_T__p   = 0.d0
                    PneoR_T__n   = 0.d0 ; PneoZ_T__n   = 0.d0
                    PneoR_UR     = 0.d0 ; PneoZ_UR     = 0.d0
                    PneoR_UZ     = 0.d0 ; PneoZ_UZ     = 0.d0
                  endif

                  Qconv_UR_AR__p  =                                           - v * rho0 * ( VdiaGradUR_AR__p  - VdiaP0_AR__p *Up0 / R )
                  Qconv_UR_AR__n  =                                           - v * rho0 * ( VdiaGradUR_AR__n  - VdiaP0_AR__n *Up0 / R )
                  Qconv_UR_AZ__p  =                                           - v * rho0 * ( VdiaGradUR_AZ__p  - VdiaP0_AZ__p *Up0 / R )
                  Qconv_UR_AZ__n  =                                           - v * rho0 * ( VdiaGradUR_AZ__n  - VdiaP0_AZ__n *Up0 / R )
                  Qconv_UR_A3__p  =                                           - v * rho0 * ( VdiaGradUR_A3__p  - VdiaP0_A3__p *Up0 / R )
                  Qconv_UR_A3__n  =                                           - v * rho0 * ( VdiaGradUR_A3__n  - VdiaP0_A3__n *Up0 / R )
                  Qconv_UR_rho__p = - v * rho  * ( UgradUR - Up0**2  / R )    - v * rho0 * ( VdiaGradUR_rho__p - VdiaP0_rho__p*Up0 / R ) &
                                                                              - v * rho  * ( VdiaGradUR        - VdiaP0       *Up0 / R )
                  Qconv_UR_rho__n =                                           - v * rho0 * ( VdiaGradUR_rho__n - VdiaP0_rho__n*Up0 / R )
                  Qconv_UR_T__p   =                                           - v * rho0 * ( VdiaGradUR_T__p   - VdiaP0_T__p  *Up0 / R )
                  Qconv_UR_T__n   =                                           - v * rho0 * ( VdiaGradUR_T__n   - VdiaP0_T__n  *Up0 / R )
                  Qconv_UR_UR__p  = - v * rho0 * ( UgradUR_UR__p         )    - v * rho0 * ( VdiaGradUR_UR__p                          )
                  Qconv_UR_UR__n  = - v * rho0 * ( UgradUR_UR__n         )    - v * rho0 * ( VdiaGradUR_UR__n                          )
                  Qconv_UR_UZ__p  = - v * rho0 * ( UgradUR_UZ            )
                  Qconv_UR_UZ__n  = 0.d0
                  Qconv_UR_Up__p  = - v * rho0 * ( UgradUR_Up-2.0*Up*Up0**2/R)- v * rho0 * (                   - VdiaP0       *Up  / R )
                  Qconv_UR_Up__n  = 0.d0

                  Qconv_UZ_AR__p  =                                           - v * rho0 * VdiaGradUZ_AR__p 
                  Qconv_UZ_AR__n  =                                           - v * rho0 * VdiaGradUZ_AR__n 
                  Qconv_UZ_AZ__p  =                                           - v * rho0 * VdiaGradUZ_AZ__p 
                  Qconv_UZ_AZ__n  =                                           - v * rho0 * VdiaGradUZ_AZ__n 
                  Qconv_UZ_A3__p  =                                           - v * rho0 * VdiaGradUZ_A3__p 
                  Qconv_UZ_A3__n  =                                           - v * rho0 * VdiaGradUZ_A3__n 
                  Qconv_UZ_rho__p = - v * rho  * UgradUZ                      - v * rho0 * VdiaGradUZ_rho__p - v * rho  * VdiaGradUZ
                  Qconv_UZ_rho__n =                                           - v * rho0 * VdiaGradUZ_rho__n
                  Qconv_UZ_T__p   =                                           - v * rho0 * VdiaGradUZ_T__p
                  Qconv_UZ_T__n   =                                           - v * rho0 * VdiaGradUZ_T__n
                  Qconv_UZ_UR__p  = - v * rho0 * UgradUZ_UR
                  Qconv_UZ_UR__n  = 0.d0
                  Qconv_UZ_UZ__p  = - v * rho0 * UgradUZ_UZ__p                - v * rho0 * VdiaGradUZ_UZ__p
                  Qconv_UZ_UZ__n  = - v * rho0 * UgradUZ_UZ__n                - v * rho0 * VdiaGradUZ_UZ__n
                  Qconv_UZ_Up__p  = - v * rho0 * UgradUZ_Up
                  Qconv_UZ_Up__n  = 0.d0

                  Qconv_Up_AR__p  =                                           - v * rho0 * ( VdiaGradUp_AR__p  + VdiaP0_AR__p *UR0 / R )
                  Qconv_Up_AR__n  =                                           - v * rho0 * ( VdiaGradUp_AR__n  + VdiaP0_AR__n *UR0 / R )
                  Qconv_Up_AZ__p  =                                           - v * rho0 * ( VdiaGradUp_AZ__p  + VdiaP0_AZ__p *UR0 / R )
                  Qconv_Up_AZ__n  =                                           - v * rho0 * ( VdiaGradUp_AZ__n  + VdiaP0_AZ__n *UR0 / R )
                  Qconv_Up_A3__p  =                                           - v * rho0 * ( VdiaGradUp_A3__p  + VdiaP0_A3__p *UR0 / R )
                  Qconv_Up_A3__n  =                                           - v * rho0 * ( VdiaGradUp_A3__n  + VdiaP0_A3__n *UR0 / R )
                  Qconv_Up_rho__p = - v * rho  * ( UgradUp + UR0*Up0 / R )    - v * rho0 * ( VdiaGradUp_rho__p + VdiaP0_rho__p*UR0 / R ) &
                                                                              - v * rho  * ( VdiaGradUp        + VdiaP0       *UR0 / R )
                  Qconv_Up_rho__n =                                           - v * rho0 * ( VdiaGradUp_rho__n + VdiaP0_rho__n*UR0 / R )
                  Qconv_Up_T__p   =                                           - v * rho0 * ( VdiaGradUp_T__p   + VdiaP0_T__p  *UR0 / R )
                  Qconv_Up_T__n   =                                           - v * rho0 * ( VdiaGradUp_T__n   + VdiaP0_T__n  *UR0 / R )
                  Qconv_Up_UR__p  = - v * rho0 * ( UgradUp_UR    + UR*Up0 / R)- v * rho0 * (                   + VdiaP0       *UR  / R )
                  Qconv_Up_UR__n  = 0.d0
                  Qconv_Up_UZ__p  = - v * rho0 * ( UgradUp_UZ                )
                  Qconv_Up_UZ__n  = 0.d0
                  Qconv_Up_Up__p  = - v * rho0 * ( UgradUp_Up__p + UR0*Up / R)- v * rho0 * ( VdiaGradUp_Up__p                          )
                  Qconv_Up_Up__n  = - v * rho0 * ( UgradUp_Up__n             )- v * rho0 * ( VdiaGradUp_Up__n                          )

                  BgradVstar_AR__p = BR0_AR    * v_R
                  BgradVstar_AR__n = BZ0_AR__n * v_Z
                  BgradVstar_AR__k = Bp0_AR    * v_p / R

                  BgradVstar_AZ__n = BR0_AZ__n * v_R
                  BgradVstar_AZ__p = BZ0_AZ    * v_Z
                  BgradVstar_AZ__k = Bp0_AZ    * v_p / R

                  BgradVstar_A3__p = BR0_A3 * v_R  + BZ0_A3 * v_Z
                  BgradVstar_A3__k = Bp0_A3 * v_p / R

                  VdiaGradVstar_AR__p   = VdiaR0_AR__p * v_R + VdiaZ0_AR__p * v_Z
                  VdiaGradVstar_AR__n   = VdiaR0_AR__n * v_R + VdiaZ0_AR__n * v_Z
                  VdiaGradVstar_AR__k   = VdiaP0_AR__p * v_p / R
                  VdiaGradVstar_AR__kn  = VdiaP0_AR__n * v_p / R

                  VdiaGradVstar_AZ__p   = VdiaR0_AZ__p * v_R + VdiaZ0_AZ__p * v_Z
                  VdiaGradVstar_AZ__n   = VdiaR0_AZ__n * v_R + VdiaZ0_AZ__n * v_Z
                  VdiaGradVstar_AZ__k   = VdiaP0_AZ__p * v_p / R
                  VdiaGradVstar_AZ__kn  = VdiaP0_AZ__n * v_p / R

                  VdiaGradVstar_A3__p   = VdiaR0_A3__p * v_R + VdiaZ0_A3__p * v_Z
                  VdiaGradVstar_A3__n   = VdiaR0_A3__n * v_R + VdiaZ0_A3__n * v_Z
                  VdiaGradVstar_A3__k   = VdiaP0_A3__p * v_p / R
                  VdiaGradVstar_A3__kn  = VdiaP0_A3__n * v_p / R

                  VdiaGradVstar_rho__p  = VdiaR0_rho__p * v_R + VdiaZ0_rho__p * v_Z
                  VdiaGradVstar_rho__n  = VdiaR0_rho__n * v_R + VdiaZ0_rho__n * v_Z
                  VdiaGradVstar_rho__k  = VdiaP0_rho__p * v_p / R
                  VdiaGradVstar_rho__kn = VdiaP0_rho__n * v_p / R

                  VdiaGradVstar_T__p    = VdiaR0_T__p * v_R + VdiaZ0_T__p * v_Z
                  VdiaGradVstar_T__n    = VdiaR0_T__n * v_R + VdiaZ0_T__n * v_Z
                  VdiaGradVstar_T__k    = VdiaP0_T__p * v_p / R
                  VdiaGradVstar_T__kn   = VdiaP0_T__n * v_p / R

                  JxB_UR_AR__p  = - BR0_AR * BgradVstar__p - BR0 * BgradVstar_AR__p &
                                  - 2.0*Bp0*Bp0_AR * v / R + ( v_R + v / R ) * ( 0.5d0 * BB2_AR__p )
                  JxB_UR_AR__n  = - BR0 * BgradVstar_AR__n + ( v_R + v / R ) * ( 0.5d0 * BB2_AR__n )
                  JxB_UR_AR__k  = - BR0_AR * BgradVstar__k - BR0 * BgradVstar_AR__k
                  JxB_UR_AR__kn = 0.d0

                  JxB_UR_AZ__p  = - BR0       * BgradVstar_AZ__p - 2.0*Bp0*Bp0_AZ * v / R + ( v_R + v / R ) * ( 0.5d0 * BB2_AZ__p )
                  JxB_UR_AZ__n  = - BR0_AZ__n * BgradVstar__p    - BR0 * BgradVstar_AZ__n + ( v_R + v / R ) * ( 0.5d0 * BB2_AZ__n )
                  JxB_UR_AZ__k  = - BR0       * BgradVstar_AZ__k
                  JxB_UR_AZ__kn = - BR0_AZ__n * BgradVstar__k

                  JxB_UR_A3__p  = - BR0_A3 * BgradVstar__p - BR0 * BgradVstar_A3__p &
                                  - 2.0*Bp0*Bp0_A3 * v / R + ( v_R + v / R ) * ( 0.5d0 * BB2_A3 )
                  JxB_UR_A3__n  = 0.d0
                  JxB_UR_A3__k  = - BR0_A3 * BgradVstar__k - BR0 * BgradVstar_A3__k
                  JxB_UR_A3__kn = 0.d0

                  JxB_UZ_AR__p  = - BZ0       * BgradVstar_AR__p                          + v_Z * ( 0.5d0 * BB2_AR__p)
                  JxB_UZ_AR__n  = - BZ0_AR__n * BgradVstar__p - BZ0 * BgradVstar_AR__n    + v_Z * ( 0.5d0 * BB2_AR__n)
                  JxB_UZ_AR__k  = - BZ0       * BgradVstar_AR__k
                  JxB_UZ_AR__kn = - BZ0_AR__n * BgradVstar__k

                  JxB_UZ_AZ__p  = - BZ0_AZ * BgradVstar__p - BZ0 * BgradVstar_AZ__p + v_Z * ( 0.5d0 * BB2_AZ__p)
                  JxB_UZ_AZ__n  = - BZ0    * BgradVstar_AZ__n                       + v_Z * ( 0.5d0 * BB2_AZ__n)
                  JxB_UZ_AZ__k  = - BZ0_AZ * BgradVstar__k - BZ0 * BgradVstar_AZ__k
                  JxB_UZ_AZ__kn = 0.d0

                  JxB_UZ_A3__p  = - BZ0_A3 * BgradVstar__p - BZ0 * BgradVstar_A3__p + v_Z * ( 0.5d0 * BB2_A3)
                  JxB_UZ_A3__n  = 0.d0
                  JxB_UZ_A3__k  = - BZ0_A3 * BgradVstar__k - BZ0 * BgradVstar_A3__k
                  JxB_UZ_A3__kn = 0.d0

                  JxB_Up_AR__p  = - Bp0_AR * BgradVstar__p - Bp0 * BgradVstar_AR__p + Bp0_AR*BR0* v / R + Bp0*BR0_AR* v / R
                  JxB_Up_AR__n  = - Bp0    * BgradVstar_AR__n
                  JxB_Up_AR__k  = - Bp0_AR * BgradVstar__k - Bp0 * BgradVstar_AR__k + ( v_p/ R ) * ( 0.5d0 * BB2_AR__p )
                  JxB_Up_AR__kn = + ( v_p/ R ) * ( 0.5d0 * BB2_AR__n )

                  JxB_Up_AZ__p  = - Bp0_AZ * BgradVstar__p - Bp0 * BgradVstar_AZ__p + Bp0_AZ*BR0* v / R
                  JxB_Up_AZ__n  = - Bp0    * BgradVstar_AZ__n                       + Bp0*BR0_AZ__n* v / R
                  JxB_Up_AZ__k  = - Bp0_AZ * BgradVstar__k - Bp0 * BgradVstar_AZ__k + ( v_p/ R ) * ( 0.5d0 * BB2_AZ__p )
                  JxB_Up_AZ__kn = + ( v_p/ R ) * ( 0.5d0 * BB2_AZ__n )

                  JxB_Up_A3__p  = - Bp0_A3 * BgradVstar__p - Bp0 * BgradVstar_A3__p + Bp0_A3*BR0* v / R + Bp0*BR0_A3* v / R
                  JxB_Up_A3__n  = 0.d0
                  JxB_Up_A3__k  = - Bp0_A3 * BgradVstar__k - Bp0 * BgradVstar_A3__k + ( v_p/ R ) * ( 0.5d0 * BB2_A3 )
                  JxB_Up_A3__kn = 0.d0

                  gradRho_gradVstar_rho__p  = rho_R * v_R  + rho_Z * v_Z
                  gradRho_gradVstar_rho__kn = (rho_p / R) * (v_p  / R)

                  gradT_gradVstar_T__p  = T_R * v_R  + T_Z * v_Z
                  gradT_gradVstar_T__kn = (T_p / R) * (v_p  / R)

                  eta_T_T    = deta_dT * T
                  eta_R_T    = d2eta_d2T * T * T0_R + deta_dT * T_R
                  eta_Z_T    = d2eta_d2T * T * T0_Z + deta_dT * T_Z
                  eta_p_T__p = d2eta_d2T * T * T0_p
                  eta_p_T__n = deta_dT * T_p

                  JR0_AR__p  = 0.d0  ; JR0_AR__n  = 0.d0  ; JR0_AR__nn = 0.d0
                  JR0_AZ__p  = 0.d0  ; JR0_AZ__n  = 0.d0  ; JR0_AZ__nn = 0.d0
                  JR0_A3__p  = 0.d0  ; JR0_A3__n  = 0.d0  ; JR0_A3__nn = 0.d0

                  JZ0_AR__p  = 0.d0  ; JZ0_AR__n  = 0.d0  ; JZ0_AR__nn = 0.d0
                  JZ0_AZ__p  = 0.d0  ; JZ0_AZ__n  = 0.d0  ; JZ0_AZ__nn = 0.d0
                  JZ0_A3__p  = 0.d0  ; JZ0_A3__n  = 0.d0  ; JZ0_A3__nn = 0.d0

                  Jp0_AR__p  = 0.d0  ; Jp0_AR__n  = 0.d0  ; Jp0_AR__nn = 0.d0
                  Jp0_AZ__p  = 0.d0  ; Jp0_AZ__n  = 0.d0  ; Jp0_AZ__nn = 0.d0
                  Jp0_A3__p  = 0.d0  ; Jp0_A3__n  = 0.d0  ; Jp0_A3__nn = 0.d0

                  JJ2_AR__p  = 0.d0  ; JJ2_AR__n  = 0.d0  ; JJ2_AR__nn = 0.d0
                  JJ2_AZ__p  = 0.d0  ; JJ2_AZ__n  = 0.d0  ; JJ2_AZ__nn = 0.d0
                  JJ2_A3__p  = 0.d0  ; JJ2_A3__n  = 0.d0  ; JJ2_A3__nn = 0.d0

                  if(eta_ohmic .gt. 1.d-12)then
                    JR0_AR__p  = - AR_ZZ
                    JR0_AR__n  = 0.d0
                    JR0_AR__nn = - AR_pp / R**2
                    JR0_AZ__p  = AZ_RZ
                    JR0_AZ__n  = 0.0d0
                    JR0_AZ__nn = 0.0d0
                    JR0_A3__p  = dF_dpsi(ms,mt) * A3_Z / R
                    JR0_A3__n  = A3_Rp / R**2
                    JR0_A3__nn = 0.d0

                    JZ0_AR__p  = AR_RZ / R**2 + AR_Z / R
                    JZ0_AR__n  = 0.0d0
                    JZ0_AR__nn = 0.0d0
                    JZ0_AZ__p  = - AZ_RR - AZ_R / R
                    JZ0_AZ__n  = 0.d0
                    JZ0_AZ__nn = - AZ_pp / R**2
                    JZ0_A3__p  = - (dF_dpsi(ms,mt)*A3_R) / R ! + Fprofile(ms,mt)/R**2 - Fprofile(ms,mt)/R**2
                    JZ0_A3__n  = A3_Zp / R**2
                    JZ0_A3__nn = 0.d0

                    Jp0_AR__p  = 0.0d0
                    Jp0_AR__n  = AR_Rp / R - AR_p / R**2
                    Jp0_AR__nn = 0.d0
                    Jp0_AZ__p  = 0.0d0
                    Jp0_AZ__n  = AZ_Zp / R
                    Jp0_AZ__nn = 0.d0
                    Jp0_A3__p  = - A3_RR / R + A3_R / R**2 - A3_ZZ / R
                    Jp0_A3__n  = 0.0d0
                    Jp0_A3__nn = 0.0d0

                    JJ2_AR__p  = 2.0d0 * (JR0*JR0_AR__p  + JZ0*JZ0_AR__p  + Jp0*Jp0_AR__p) 
                    JJ2_AR__n  = 2.0d0 * (JR0*JR0_AR__n  + JZ0*JZ0_AR__n  + Jp0*Jp0_AR__n)
                    JJ2_AR__nn = 2.0d0 * (JR0*JR0_AR__nn + JZ0*JZ0_AR__nn + Jp0*Jp0_AR__nn)

                    JJ2_AZ__p  = 2.0d0 * (JR0*JR0_AZ__p  + JZ0*JZ0_AZ__p  + Jp0*Jp0_AZ__p)
                    JJ2_AZ__n  = 2.0d0 * (JR0*JR0_AZ__n  + JZ0*JZ0_AZ__n  + Jp0*Jp0_AZ__n)
                    JJ2_AZ__nn = 2.0d0 * (JR0*JR0_AZ__nn + JZ0*JZ0_AZ__nn + Jp0*Jp0_AZ__nn)

                    JJ2_A3__p  = 2.0d0 * (JR0*JR0_A3__p  + JZ0*JZ0_A3__p  + Jp0*Jp0_A3__p)
                    JJ2_A3__n  = 2.0d0 * (JR0*JR0_A3__n  + JZ0*JZ0_A3__n  + Jp0*Jp0_A3__n)
                    JJ2_A3__nn = 2.0d0 * (JR0*JR0_A3__nn + JZ0*JZ0_A3__nn + Jp0*Jp0_A3__nn)

                  endif

                  ! --- VMS variables
                  CvGradVj__p = CvR0 * bf_R + CvZ0 * bf_Z
                  CvGradVj__n = Cvp0 * bf_p / R
                  VbGradVj__p = VbR0 * bf_R + VbZ0 * bf_Z
                  VbGradVj__n = Vbp0 * bf_p / R

                  DiveRMVj    = bf_R + bf / R
                  DiveZMVj    = bf_z
                  DivePMVj__n = bf_p / R

                  ! --- Viscous terms
                  Qvisc_UR_UR__p = 0.d0 ; Qvisc_UR_UR__n = 0.d0 ; Qvisc_UR_UR__k = 0.d0 ; Qvisc_UR_UR__kn = 0.d0
                  Qvisc_UR_UZ__p = 0.d0 ; Qvisc_UR_UZ__n = 0.d0 ; Qvisc_UR_UZ__k = 0.d0 ; Qvisc_UR_UZ__kn = 0.d0
                  Qvisc_UR_Up__p = 0.d0 ; Qvisc_UR_Up__n = 0.d0 ; Qvisc_UR_Up__k = 0.d0 ; Qvisc_UR_Up__kn = 0.d0
                  Qvisc_UZ_UR__p = 0.d0 ; Qvisc_UZ_UR__n = 0.d0 ; Qvisc_UZ_UR__k = 0.d0 ; Qvisc_UZ_UR__kn = 0.d0
                  Qvisc_UZ_UZ__p = 0.d0 ; Qvisc_UZ_UZ__n = 0.d0 ; Qvisc_UZ_UZ__k = 0.d0 ; Qvisc_UZ_UZ__kn = 0.d0
                  Qvisc_UZ_Up__p = 0.d0 ; Qvisc_UZ_Up__n = 0.d0 ; Qvisc_UZ_Up__k = 0.d0 ; Qvisc_UZ_Up__kn = 0.d0
                  Qvisc_Up_UR__p = 0.d0 ; Qvisc_Up_UR__n = 0.d0 ; Qvisc_Up_UR__k = 0.d0 ; Qvisc_Up_UR__kn = 0.d0
                  Qvisc_Up_UZ__p = 0.d0 ; Qvisc_Up_UZ__n = 0.d0 ; Qvisc_Up_UZ__k = 0.d0 ; Qvisc_Up_UZ__kn = 0.d0
                  Qvisc_Up_Up__p = 0.d0 ; Qvisc_Up_Up__n = 0.d0 ; Qvisc_Up_Up__k = 0.d0 ; Qvisc_Up_Up__kn = 0.d0
                  Qvisc_T_UR__p  = 0.d0 ; Qvisc_T_UR__n  = 0.d0
                  Qvisc_T_UZ__p  = 0.d0 ; Qvisc_T_UZ__n  = 0.d0
                  Qvisc_T_Up__p  = 0.d0 ; Qvisc_T_Up__n  = 0.d0
                  Qvisc_T_T__p   = 0.d0 ; Qvisc_T_T__n   = 0.d0
                  select case(ViscType)
                  case(0) ! --- Full-MHD viscosity
                    Qvisc_UR_UR__p  = - (v_R * UR_R + v_Z * UR_Z) - v * ( UR / R**2 ) 
                    Qvisc_UR_Up__n  = - v * (+ 2.d0 * Up_p / R**2 ) 
                    Qvisc_UR_UR__kn = - (v_p * UR_p / R**2)
                    Qvisc_UZ_UZ__p  = - (v_R * UZ_R + v_Z * UZ_Z)
                    Qvisc_UZ_UZ__kn = - (v_p * UZ_p / R**2)
                    Qvisc_Up_Up__p  = - (v_R * UP_R + v_Z * UP_Z) - v * ( Up / R**2 )
                    Qvisc_Up_UR__n  = - v * ( - 2.d0 * UR_p / R**2 )
                    Qvisc_Up_Up__kn = - (v_p * uP_p / R**2)
                  case(15) ! --- Mock-up of reduced-MHD viscosity
                    Qvisc_UR_UR__p = + ( - UR_Z ) * v_Z
                    Qvisc_UR_UZ__p = + (   UZ_R ) * v_Z
                    Qvisc_UZ_UR__p = - ( - UR_Z ) * v_R
                    Qvisc_UZ_UZ__p = - (   UZ_R ) * v_R
                    Qvisc_Up_Up__p = - ( Up_R + Up/R ) * ( v_R + v/R ) - Up_z * v_z 
                  case(10) ! --- Generic form of full viscosity tensor (from B.Nkonga and A.Bhole)
                    Qvisc_UR_UR__p  = - ( v_R * UR_R + v_Z * UR_Z )      &
                                      - ( v_R * UR_R              )      &
                                      - v * UR / R**2                    &
                                      - v * UR / R**2
                    Qvisc_UR_UZ__p  = - (            + v_Z * UZ_R )
                    Qvisc_UR_Up__p  = - v * Up_R /R**2
                    Qvisc_UR_Up__n  = - (            + v *   Up_p )/ R**2

                    Qvisc_UR_UR__kn = - ( v_p * UR_p / R**2 )
                    Qvisc_UR_Up__k  = - ( v_p * Up_R / R    ) &
                                      - ( - v_p * Up )/ R**2

                    Qvisc_UZ_UR__p  = - ( v_R * UR_Z )
                    Qvisc_UZ_UZ__p  = - ( v_R * UZ_R + v_Z * UZ_Z ) &
                                      - (            + v_Z * UZ_Z )

                    Qvisc_UZ_UZ__kn = - ( v_p * UZ_p / R**2 )
                    Qvisc_UZ_Up__k  = - ( v_p * Up_Z / R    )

                    Qvisc_Up_UR__n  = - ( v_R * UR_p )/R    &
                                      - ( - v * UR_p )/R**2 
                    Qvisc_Up_UZ__n  = - ( v_Z * UZ_p )/R
                    Qvisc_Up_Up__p  = - ( v_R * Up_R + v_Z * Up_Z )       &
                                      - (            + v   * Up   )/ R**2 &
                                      + ( v_R * Up   + v   * Up_R )/ R

                    Qvisc_Up_UR__k  = - ( v_p * UR )/ R**2        &
                                      - ( v_p * UR )/ R**2
                    Qvisc_Up_Up__kn = - ( v_p * Up_p / R**2)      &
                                      - ( v_p * Up_p / R   )/R

                    QviscT0_UR__p   = - ( 2.0 * UR0_R * UR_R + 2.0 * UR0_Z * UR_Z ) &
                                      - ( 2.0 * UR0_R * UR_R + UR_Z * UZ0_R )       &
                                      - ( 2.0 * UR0 * UR  + UR * Up0_p )/ R**2      &
                                      - UR  * ( UR0 +  Up0_R )/R**2                 &
                                      - UR0 * ( UR           )/R**2                 &
                                      - ( UZ0_R * UR_Z )                            &
                                      - ( Up0_p * UR )/ R**2                        &
                                      - ( Up0_p * UR )/ R**2
                    QviscT0_UR__n   = - ( 2.0 * UR0_p * UR_p / R**2) &
                                      - (       UR_p  * Up0_R / R  ) &
                                      - ( - UR_p * Up0 )/ R**2       &
                                      - ( Up0_R * UR_p )/R           &
                                      - ( - Up0 * UR_p )/ R**2
                    QviscT0_UZ__p   = - ( UR0_Z * UZ_R )                            &
                                      - ( 2.0 * UZ0_R * UZ_R + 2.0 * UZ0_Z * UZ_Z ) &
                                      - ( UZ_R * UR0_Z + 2.0 * UZ0_Z * UZ_Z )
                    QviscT0_UZ__n   = - ( 2.0 * UZ0_p * UZ_p / R**2) &
                                      - ( UZ_p * Up0_Z / R )         &
                                      - ( Up0_Z * UZ_p )/R
                    QviscT0_Up__p   = - ( UR0_p * Up_R / R )                        &
                                      - ( - UR0_p * Up )/ R**2                      &
                                      - UR0 * ( Up_R )/R**2                         &
                                      - ( UZ0_p * Up_Z / R )                        &
                                      - ( 2.0 * Up0_R * Up_R + 2.0 * Up0_Z * Up_Z ) &
                                      - ( Up_R * UR0_p + Up_Z * UZ0_p )/R           &
                                      - ( - Up * UR0_p + 2.0 * Up0 * Up )/ R**2     &
                                      + ( Up_R * Up0 + Up0_R * Up +  Up * Up0_R +  Up0 * Up_R )/ R
                    QviscT0_Up__n   = - ( UR0 * Up_p )/ R**2            &
                                      - ( 2.0 * Up0_p * Up_p / R**2)    &
                                      - ( 2.0 * Up0_p * Up_p / R )/R    &
                                      - ( Up_p * UR0 )/ R**2            &
                                      - ( Up_p * UR0 )/ R**2

                    Qvisc_T_UR__p   = visco_divV * 2.0*divU*divU_UR - visco_T * QviscT0_UR__p
                    Qvisc_T_UR__n   =                               - visco_T * QviscT0_UR__n

                    Qvisc_T_UZ__p   = visco_divV * 2.0*divU*divU_UZ - visco_T * QviscT0_UZ__p
                    Qvisc_T_UZ__n   =                               - visco_T * QviscT0_UZ__n

                    Qvisc_T_Up__p   =                                  - visco_T * QviscT0_Up__p
                    Qvisc_T_Up__n   = visco_divV * 2.0*divU*divU_Up__n - visco_T * QviscT0_Up__n
                  case (20) ! --- Viscosity as implemented by W.Haverkort
                    Qvisc_UR_UR__p  = - (2.d0 * UR_R * v_R + (UR_Z ) * v_Z )
                    Qvisc_UR_UZ__p  = - ( (UZ_R ) * v_Z )
                    Qvisc_UR_Up__n  = - (2.d0 * Up_p * v / R**2 )

                    Qvisc_UR_UR__kn = - ((UR_p) * v_p / R**2)
                    Qvisc_UR_Up__k  = - ((+ R * Up_R - Up) * v_p / R**2)

                    Qvisc_UZ_UR__p  = - ( (UR_Z ) * v_R)
                    Qvisc_UZ_UZ__p  = - (2.d0 * UZ_Z * v_Z + (UZ_R) * v_R)

                    Qvisc_UZ_UZ__kn = - ((UZ_p) * v_p / R**2)
                    Qvisc_UZ_Up__k  = - (( R * Up_Z ) * v_p / R**2)

                    Qvisc_Up_UR__n  = - (+ (- UR_p / R ) * v / R &
                                         + (+ UR_p / R ) * v_R )
                    Qvisc_Up_UZ__n  = - (+ (+ UZ_p) * v_Z / R )
                    Qvisc_Up_Up__p  = - (+ ( 3.d0 * Up / R - Up_R  ) * v / R &
                                         + ( - Up / R + Up_R ) * v_R            &
                                         + (R * Up_Z) * v_Z / R )

                    Qvisc_Up_UR__k  = - (2.d0 * ( UR ) * v_p / R**2)
                    Qvisc_Up_Up__kn = - (2.d0 * ( Up_p ) * v_p / R**2)
                  case (22) ! --- The correct viscosity (there is only one...)
                    Qvisc_UR_UR__p  = - (v_R + v/R) * divU_UR + v_Z * (- UR_Z)
                    Qvisc_UR_UR__k  = 0.d0
                    Qvisc_UR_UR__n  = 0.d0
                    Qvisc_UR_UR__kn = - v_p/R**2 * (UR_p)

                    Qvisc_UR_UZ__p  = - (v_R + v/R) * divU_UZ + v_Z * (UZ_R)
                    Qvisc_UR_UZ__k  = 0.d0
                    Qvisc_UR_UZ__n  = 0.d0
                    Qvisc_UR_UZ__kn = 0.d0

                    Qvisc_UR_Up__p  = 0.d0
                    Qvisc_UR_Up__k  = - v_p/R**2 * (- Up - R*Up_R)
                    Qvisc_UR_Up__n  = - (v_R + v/R) * divU_Up__n
                    Qvisc_UR_Up__kn = 0.d0

                    Qvisc_UZ_UR__p  = - v_Z * divU_UR - v_R * (- UR_Z)
                    Qvisc_UZ_UR__k  = 0.d0
                    Qvisc_UZ_UR__n  = 0.d0
                    Qvisc_UZ_UR__kn = 0.d0

                    Qvisc_UZ_UZ__p  = - v_Z * divU_UZ - v_R * (UZ_R)
                    Qvisc_UZ_UZ__k  = 0.d0
                    Qvisc_UZ_UZ__n  = 0.d0
                    Qvisc_UZ_UZ__kn = + v_p/R**2 * (- UZ_p)

                    Qvisc_UZ_Up__p  = 0.d0
                    Qvisc_UZ_Up__k  = + v_p/R**2 * (R*Up_Z)
                    Qvisc_UZ_Up__n  = - v_Z * divU_Up__n
                    Qvisc_UZ_Up__kn = 0.d0

                    Qvisc_Up_UR__p  = 0.d0
                    Qvisc_Up_UR__k  = - v_p/R * divU_UR
                    Qvisc_Up_UR__n  = + (v_R + v/R)/R * (UR_p)
                    Qvisc_Up_UR__kn = 0.d0

                    Qvisc_Up_UZ__p  = 0.d0
                    Qvisc_Up_UZ__k  = - v_p/R * divU_UZ
                    Qvisc_Up_UZ__n  = - v_Z/R * (- UZ_p)
                    Qvisc_Up_UZ__kn = 0.d0

                    Qvisc_Up_Up__p  = - v_Z/R * (R*Up_Z) + (v_R + v/R)/R * (- Up - R*Up_R)
                    Qvisc_Up_Up__k  = 0.d0
                    Qvisc_Up_Up__n  = 0.d0
                    Qvisc_Up_Up__kn = - v_p/R * divU_Up__n
                  case (30)
                  case default
                  end select

                  lap_bf = bf_R / R + bf_RR + bf_ZZ

                  amat      = 0.d0
                  Pjac      = 0.d0 ! time derivative part
                  Qjac_p    = 0.d0 ! rest of the LHS (poloidal part)
                  Qjac_k    = 0.d0 ! rest of the LHS (toroidal part with phi-derivatives of the test-function)
                  Qjac_n    = 0.d0 ! rest of the LHS (toroidal part with phi-derivatives of the basis-functions)
                  Qjac_kn   = 0.d0 ! rest of the LHS (toroidal part with phi-derivatives of the basis-functions and the test-functions)
                  Qjac_pnn = 0.d0

                  QvmsF_p   = 0.d0 ! VMS terms
                  QvmsF_n   = 0.d0
                  QvmsF_k   = 0.d0
                  QvmsF_kn  = 0.d0
                  QvmsAd_p  = 0.d0
                  QvmsAd_n  = 0.d0
                  QvmsAd_k  = 0.d0
                  QvmsAd_kn = 0.d0

                  !###################################################################################################
                  !#  equation 1   (R component induction equation)                                                  #
                  !###################################################################################################
                  Pjac   (var_AR,var_AR) =   v * AR 

                  Qjac_p (var_AR,var_AR) = + v * (UZ0   * Bp0_AR ) &
                                           + tauIC_ARAZ * v * tau_IC*F0/rho0_corr/BB2    * BR0 * BgradP_AR__p &
                                           - tauIC_ARAZ * v * tau_IC*F0/rho0_corr/BB2**2 * BR0 * BgradP * BB2_AR__p &
                                           + eta_ARAZ * v * (eta_Z * Bp0_AR ) &
                                           - eta_ARAZ * eta_T * ( - v_Z * Bp0_AR ) &
                                           - eta_ARAZ_const   * ( - v_Z * Bp0_AR ) &
                                           - eta_ARAZ * eta_num * lap_Vstar * lap_bf
                  Qjac_n (var_AR,var_AR) = + v * (- Up0   * BZ0_AR__n     ) &
                                           + tauIC_ARAZ * v * tau_IC*F0/rho0_corr/BB2    * BR0 * BgradP_AR__n &
                                           - tauIC_ARAZ * v * tau_IC*F0/rho0_corr/BB2**2 * BR0 * BgradP * BB2_AR__n &
                                           + eta_ARAZ * v * (- eta_p * BZ0_AR__n / R )
                  Qjac_kn(var_AR,var_AR) = - eta_ARAZ * eta_T * ( + v_p * BZ0_AR__n / R)

                  Qjac_p (var_AR,var_AZ) = + v * (UZ0 * Bp0_AZ - Up0 * BZ0_AZ)          &
                                           + tauIC_ARAZ * v * tau_IC*F0/rho0_corr/BB2    * BR0 * BgradP_AZ__p &
                                           - tauIC_ARAZ * v * tau_IC*F0/rho0_corr/BB2**2 * BR0 * BgradP * BB2_AZ__p &
                                           + eta_ARAZ * v * (eta_Z * Bp0_AZ - eta_p * BZ0_AZ / R ) &
                                           - eta_ARAZ * eta_T * ( - v_Z * Bp0_AZ ) &
                                           - eta_ARAZ_const   * ( - v_Z * Bp0_AZ )
                  Qjac_n (var_AR,var_AZ) = + tauIC_ARAZ * v * tau_IC*F0/rho0_corr/BB2    * BR0_AZ__n * BgradP &
                                           + tauIC_ARAZ * v * tau_IC*F0/rho0_corr/BB2    * BR0       * BgradP_AZ__n &
                                           - tauIC_ARAZ * v * tau_IC*F0/rho0_corr/BB2**2 * BR0       * BgradP * BB2_AZ__n
                  Qjac_k (var_AR,var_AZ) = - eta_ARAZ * eta_T * ( + v_p * BZ0_AZ / R)

                  Qjac_p (var_AR,var_A3) = + v * (UZ0 * Bp0_A3 - Up0 * BZ0_A3)          &
                                           + tauIC_ARAZ * v * tau_IC*F0/rho0_corr/BB2    * BR0_A3 * BgradP &
                                           + tauIC_ARAZ * v * tau_IC*F0/rho0_corr/BB2    * BR0    * BgradP_A3 &
                                           - tauIC_ARAZ * v * tau_IC*F0/rho0_corr/BB2**2 * BR0    * BgradP * BB2_A3 &
                                           + eta_ARAZ * v * (eta_Z * Bp0_A3 - eta_p * BZ0_A3 / R ) &
                                           - eta_ARAZ * eta_T * ( - v_Z * Bp0_A3 ) &
                                           - eta_ARAZ_const   * ( - v_Z * Bp0_A3 )
                  Qjac_k (var_AR,var_A3) = - eta_ARAZ * eta_T * ( + v_p * BZ0_A3 / R)

                  Qjac_p (var_AR,var_UZ) = + v * (  UZ * Bp0)
                  Qjac_p (var_AR,var_Up) = + v * (- Up * BZ0)

                  Qjac_p (var_AR,var_rho)= + tauIC_ARAZ * v * tau_IC*F0/rho0_corr   /BB2 * BR0 * BgradP_rho__p &
                                           - tauIC_ARAZ * v * tau_IC*F0/rho0_corr**2/BB2 * BR0 * BgradP * rho
                  Qjac_n (var_AR,var_rho)= - tauIC_ARAZ * v * tau_IC*F0/rho0_corr   /BB2 * BR0 * BgradP_rho__n

                  Qjac_p (var_AR,var_T ) = + tauIC_ARAZ * v * tau_IC*F0/rho0_corr/BB2 * BR0 * BgradP_T__p &
                                           + eta_ARAZ * v * (eta_Z_T * Bp00 - eta_p_T__p * BZ0 / R ) &
                                           - eta_ARAZ * eta_T_T * ( - v_Z * Bp00 )                   &
                                           + eta_ARAZ * eta_T_T * v * current_source_JR(ms,mt)
                  Qjac_n (var_AR,var_T ) = + tauIC_ARAZ * v * tau_IC*F0/rho0_corr/BB2 * BR0 * BgradP_T__n &
                                           + eta_ARAZ * v * (              - eta_p_T__n * BZ0 / R )
                  Qjac_k (var_AR,var_T ) = - eta_ARAZ * eta_T_T * ( + v_p * BZ0 / R)

                  !###################################################################################################
                  !#  equation 2   (Z component induction equation)                                                  #
                  !###################################################################################################
                  Pjac   (var_AZ,var_AZ) =   v * AZ 

                  Qjac_p (var_AZ,var_AR) = + v * (Up0 * BR0_AR - UR0 * Bp0_AR)         &
                                           + tauIC_ARAZ * v * tau_IC*F0/rho0_corr/BB2    * BZ0 * BgradP_AR__p &
                                           - tauIC_ARAZ * v * tau_IC*F0/rho0_corr/BB2**2 * BZ0 * BgradP * BB2_AR__p &
                                           + eta_ARAZ * v * (eta_p / R * BR0_AR - eta_R * Bp0_AR) &
                                           - eta_ARAZ * eta_T * ( + v_R * Bp0_AR ) &
                                           - eta_ARAZ_const   * ( + v_R * Bp0_AR )
                  Qjac_n (var_AZ,var_AR) = + tauIC_ARAZ * v * tau_IC*F0/rho0_corr/BB2    * BZ0_AR__n * BgradP &
                                           + tauIC_ARAZ * v * tau_IC*F0/rho0_corr/BB2    * BZ0       * BgradP_AR__n &
                                           - tauIC_ARAZ * v * tau_IC*F0/rho0_corr/BB2**2 * BZ0       * BgradP * BB2_AR__n
                  Qjac_k (var_AZ,var_AR) = - eta_ARAZ * eta_T * ( - v_p * BR0_AR / R)

                  Qjac_p (var_AZ,var_AZ) = + v * (- UR0   * Bp0_AZ) &
                                           + tauIC_ARAZ * v * tau_IC*F0/rho0_corr/BB2    * BZ0 * BgradP_AZ__p &
                                           - tauIC_ARAZ * v * tau_IC*F0/rho0_corr/BB2**2 * BZ0 * BgradP * BB2_AZ__p &
                                           + eta_ARAZ * v * (- eta_R * Bp0_AZ) &
                                           - eta_ARAZ * eta_T * ( + v_R * Bp0_AZ ) &
                                           - eta_ARAZ_const   * ( + v_R * Bp0_AZ ) &
                                           - eta_ARAZ * eta_num * lap_Vstar * lap_bf
                  Qjac_n (var_AZ,var_AZ) = + v * (Up0 * BR0_AZ__n)       &
                                           + tauIC_ARAZ * v * tau_IC*F0/rho0_corr/BB2    * BZ0 * BgradP_AZ__n &
                                           - tauIC_ARAZ * v * tau_IC*F0/rho0_corr/BB2**2 * BZ0 * BgradP * BB2_AZ__n &
                                           + eta_ARAZ * v * (eta_p / R * BR0_AZ__n)
                  Qjac_kn(var_AZ,var_AZ) = - eta_ARAZ * eta_T * ( - v_p * BR0_AZ__n / R)

                  Qjac_p (var_AZ,var_A3) = + v * (Up0 * BR0_A3 - UR0 * Bp0_A3)         &
                                           + tauIC_ARAZ * v * tau_IC*F0/rho0_corr/BB2    * BZ0_A3 * BgradP &
                                           + tauIC_ARAZ * v * tau_IC*F0/rho0_corr/BB2    * BZ0    * BgradP_A3 &
                                           - tauIC_ARAZ * v * tau_IC*F0/rho0_corr/BB2**2 * BZ0    * BgradP * BB2_A3 &
                                           + eta_ARAZ * v * (eta_p / R * BR0_A3 - eta_R * Bp0_A3) &
                                           - eta_ARAZ * eta_T * ( + v_R * Bp0_A3 ) &
                                           - eta_ARAZ_const   * ( + v_R * Bp0_A3 )
                  Qjac_k (var_AZ,var_A3) = - eta_ARAZ * eta_T * ( - v_p * BR0_A3 / R)

                  Qjac_p (var_AZ,var_UR) = + v * (- UR * Bp0)
                  Qjac_p (var_AZ,var_Up) = + v * (  Up * BR0)

                  Qjac_p (var_AZ,var_rho)= + tauIC_ARAZ * v * tau_IC*F0/rho0_corr   /BB2 * BZ0 * BgradP_rho__p &
                                           - tauIC_ARAZ * v * tau_IC*F0/rho0_corr**2/BB2 * BZ0 * BgradP * rho
                  Qjac_n (var_AZ,var_rho)= + tauIC_ARAZ * v * tau_IC*F0/rho0_corr   /BB2 * BZ0 * BgradP_rho__n

                  Qjac_p (var_AZ,var_T ) = + tauIC_ARAZ * v * tau_IC*F0/rho0_corr/BB2 * BZ0 * BgradP_T__p &
                                           + eta_ARAZ * v * (eta_p_T__p / R * BR0 - eta_R_T * Bp00) &
                                           - eta_ARAZ * eta_T_T * ( + v_R * Bp00 )                  &
                                           + eta_ARAZ * eta_T_T * v * current_source_JZ(ms,mt)
                  Qjac_n (var_AZ,var_T ) = + tauIC_ARAZ * v * tau_IC*F0/rho0_corr/BB2 * BZ0 * BgradP_T__n &
                                           + eta_ARAZ * v * (eta_p_T__n / R * BR0 )
                  Qjac_k (var_AZ,var_T ) = - eta_ARAZ * eta_T_T * ( - v_p * BR0 / R)

                  !###################################################################################################
                  !#  equation 3   (Phi component induction equation)                                                #
                  !###################################################################################################
                  Pjac   (var_A3,var_A3) =   v * A3

                  Qjac_p (var_A3,var_AR) = + R * v * (- UZ0 * BR0_AR)                &
                                           + v * tau_IC*F0/rho0_corr/BB2    * R*Bp0_AR * BgradP &
                                           + v * tau_IC*F0/rho0_corr/BB2    * R*Bp0    * BgradP_AR__p &
                                           - v * tau_IC*F0/rho0_corr/BB2**2 * R*Bp0    * BgradP * BB2_AR__p &
                                           - eta_T * v_Z * R                * BR0_AR &
                                           + R * v * (- eta_Z * BR0_AR)
                  Qjac_n (var_A3,var_AR) = + R * v * (UR0 * BZ0_AR__n)                  &
                                           + v * tau_IC*F0/rho0_corr/BB2    * R*Bp0    * BgradP_AR__n &
                                           - v * tau_IC*F0/rho0_corr/BB2**2 * R*Bp0    * BgradP * BB2_AR__n &
                                           + eta_T * ( 2.d0 * v + R * v_R ) * BZ0_AR__n &
                                           + R * v * (eta_R * BZ0_AR__n)

                  Qjac_p (var_A3,var_AZ) = + R * v * (UR0 * BZ0_AZ)                  &
                                           + v * tau_IC*F0/rho0_corr/BB2    * R*Bp0_AZ * BgradP &
                                           + v * tau_IC*F0/rho0_corr/BB2    * R*Bp0    * BgradP_AZ__p &
                                           - v * tau_IC*F0/rho0_corr/BB2**2 * R*Bp0    * BgradP * BB2_AZ__p &
                                           + eta_T * ( 2.d0 * v + R * v_R ) * BZ0_AZ &
                                            + R * v * (eta_R * BZ0_AZ)
                  Qjac_n (var_A3,var_AZ) = + R * v * (- UZ0 * BR0_AZ__n)                &
                                           + v * tau_IC*F0/rho0_corr/BB2    * R*Bp0    * BgradP_AZ__n &
                                           - v * tau_IC*F0/rho0_corr/BB2**2 * R*Bp0    * BgradP * BB2_AZ__n &
                                           - eta_T * v_Z * R                * BR0_AZ__n &
                                           + R * v * (- eta_Z * BR0_AZ__n)

                  Qjac_p (var_A3,var_A3) = + R * v * (UR0 * BZ0_A3 - UZ0 * BR0_A3)   &
                                           + v * tau_IC*F0/rho0_corr/BB2    * R*Bp0    * BgradP_A3 &
                                           - v * tau_IC*F0/rho0_corr/BB2**2 * R*Bp0    * BgradP * BB2_A3 &
                                           - eta_T * v_Z * R                * BR0_A3 &
                                           + eta_T * ( 2.d0 * v + R * v_R ) * BZ0_A3 &
                                           + R * v * (eta_R * BZ0_A3 - eta_Z * BR0_A3) &
                                           - eta_num * lap_Vstar * lap_bf

                  Qjac_p (var_A3,var_UR) = + R * v * (  UR * BZ0)
                  Qjac_p (var_A3,var_UZ) = + R * v * (- UZ * BR0)

                  Qjac_p (var_A3,var_rho)= + v * tau_IC*F0/rho0_corr   /BB2 * R*Bp0 * BgradP_rho__p &
                                           - v * tau_IC*F0/rho0_corr**2/BB2 * R*Bp0 * BgradP * rho
                  Qjac_n (var_A3,var_rho)= + v * tau_IC*F0/rho0_corr/BB2 * R*Bp0 * BgradP_rho__n

                  Qjac_p (var_A3,var_T ) = + v * tau_IC*F0/rho0_corr/BB2 * R*Bp0 * BgradP_T__p &
                                           - eta_T_T * v_Z * R                * BR0 &
                                           + eta_T_T * ( 2.d0 * v + R * v_R ) * BZ0 &
                                           + eta_T_T * v * current_source_Jp(ms,mt) &
                                           + R * v * (eta_R_T * BZ0 - eta_Z_T * BR0)
                  Qjac_n (var_A3,var_T ) = + v * tau_IC*F0/rho0_corr/BB2 * R*Bp0 * BgradP_T__n

                  !###################################################################################################
                  !#  equation 4   (R component momentum equation)                                                   #
                  !###################################################################################################
                  Pjac(var_UR,var_UR )    =   v * rho0_corr * UR

                  Qjac_p (var_UR,var_AR ) =  Qconv_UR_AR__p + JxB_UR_AR__p  - v * PneoR_AR__p
                  Qjac_n (var_UR,var_AR ) =  Qconv_UR_AR__n + JxB_UR_AR__n  - v * PneoR_AR__n
                  Qjac_k (var_UR,var_AR ) =                 + JxB_UR_AR__k
                  Qjac_kn(var_UR,var_AR ) =                 + JxB_UR_AR__kn

                  Qjac_p (var_UR,var_AZ ) =  Qconv_UR_AZ__p + JxB_UR_AZ__p  - v * PneoR_AZ__p
                  Qjac_n (var_UR,var_AZ ) =  Qconv_UR_AZ__n + JxB_UR_AZ__n  - v * PneoR_AZ__n
                  Qjac_k (var_UR,var_AZ ) =                 + JxB_UR_AZ__k
                  Qjac_kn(var_UR,var_AZ ) =                 + JxB_UR_AZ__kn

                  Qjac_p (var_UR,var_A3 ) =  Qconv_UR_A3__p + JxB_UR_A3__p  - v * PneoR_A3__p
                  Qjac_n (var_UR,var_A3 ) =  Qconv_UR_A3__n + JxB_UR_A3__n  - v * PneoR_A3__n
                  Qjac_k (var_UR,var_A3 ) =                 + JxB_UR_A3__k
                  Qjac_kn(var_UR,var_A3 ) =                 + JxB_UR_A3__kn

                  Qjac_p (var_UR,var_UR ) =  Qconv_UR_UR__p + visco_T * Qvisc_UR_UR__p  - v * PneoR_UR - v * particle_source(ms,mt) * UR &
                                            - visco_num * lap_Vstar * lap_bf
                  Qjac_n (var_UR,var_UR ) =  Qconv_UR_UR__n + visco_T * Qvisc_UR_UR__n
                  Qjac_k (var_UR,var_UR ) =                 + visco_T * Qvisc_UR_UR__k
                  Qjac_kn(var_UR,var_UR ) =                 + visco_T * Qvisc_UR_UR__kn

                  Qjac_p (var_UR,var_UZ ) =  Qconv_UR_UZ__p + visco_T * Qvisc_UR_UZ__p  - v * PneoR_UZ
                  Qjac_n (var_UR,var_UZ ) =  Qconv_UR_UZ__n + visco_T * Qvisc_UR_UZ__n
                  Qjac_k (var_UR,var_UZ ) =                 + visco_T * Qvisc_UR_UZ__k
                  Qjac_kn(var_UR,var_UZ ) =                 + visco_T * Qvisc_UR_UZ__kn

                  Qjac_p (var_UR,var_Up ) =  Qconv_UR_Up__p + visco_T * Qvisc_UR_Up__p
                  Qjac_n (var_UR,var_Up ) =  Qconv_UR_Up__n + visco_T * Qvisc_UR_Up__n
                  Qjac_k (var_UR,var_Up ) =                 + visco_T * Qvisc_UR_Up__k
                  Qjac_kn(var_UR,var_Up ) =                 + visco_T * Qvisc_UR_Up__kn

                  Qjac_p (var_UR,var_rho) =  Qconv_UR_rho__p + ( v_R + v / R ) * (rho*T0) - v * PneoR_rho__p
                  Qjac_n (var_UR,var_rho) =  Qconv_UR_rho__n                              - v * PneoR_rho__n

                  Qjac_p (var_UR,var_T  ) =  Qconv_UR_T__p   + ( v_R + v / R ) * (rho0*T) + dvisco_dT * T * Qvisc_UR__p - v * PneoR_T__p
                  Qjac_n (var_UR,var_T  ) =  Qconv_UR_T__n                                                              - v * PneoR_T__n
                  Qjac_k (var_UR,var_T  ) =                                               + dvisco_dT * T * Qvisc_UR__k

                  !###################################################################################################
                  !#  equation 5   (Z component momentum equation)                                                   #
                  !###################################################################################################
                  Pjac(var_UZ,var_UZ)      =   v * rho0_corr * UZ

                  Qjac_p (var_UZ,var_AR )  =  Qconv_UZ_AR__p + JxB_UZ_AR__p  - v * PneoZ_AR__p
                  Qjac_n (var_UZ,var_AR )  =  Qconv_UZ_AR__n + JxB_UZ_AR__n  - v * PneoZ_AR__n
                  Qjac_k (var_UZ,var_AR )  =                 + JxB_UZ_AR__k
                  Qjac_kn(var_UZ,var_AR )  =                 + JxB_UZ_AR__kn

                  Qjac_p (var_UZ,var_AZ )  =  Qconv_UZ_AZ__p + JxB_UZ_AZ__p  - v * PneoZ_AZ__p
                  Qjac_n (var_UZ,var_AZ )  =  Qconv_UZ_AZ__n + JxB_UZ_AZ__n  - v * PneoZ_AZ__n
                  Qjac_k (var_UZ,var_AZ )  =                 + JxB_UZ_AZ__k
                  Qjac_kn(var_UZ,var_AZ )  =                 + JxB_UZ_AZ__kn

                  Qjac_p (var_UZ,var_A3 )  =  Qconv_UZ_A3__p + JxB_UZ_A3__p  - v * PneoZ_A3__p
                  Qjac_n (var_UZ,var_A3 )  =  Qconv_UZ_A3__n + JxB_UZ_A3__n  - v * PneoZ_A3__n
                  Qjac_k (var_UZ,var_A3 )  =                 + JxB_UZ_A3__k
                  Qjac_kn(var_UZ,var_A3 )  =                 + JxB_UZ_A3__kn

                  Qjac_p (var_UZ,var_UR )  =  Qconv_UZ_UR__p + visco_T * Qvisc_UZ_UR__p - v * PneoZ_UR
                  Qjac_n (var_UZ,var_UR )  =  Qconv_UZ_UR__n + visco_T * Qvisc_UZ_UR__n
                  Qjac_k (var_UZ,var_UR )  =                 + visco_T * Qvisc_UZ_UR__k
                  Qjac_kn(var_UZ,var_UR )  =                 + visco_T * Qvisc_UZ_UR__kn

                  Qjac_p (var_UZ,var_UZ )  =  Qconv_UZ_UZ__p + visco_T * Qvisc_UZ_UZ__p - v * PneoZ_UZ - v * particle_source(ms,mt) * UZ &
                                            - visco_num * lap_Vstar * lap_bf
                  Qjac_n (var_UZ,var_UZ )  =  Qconv_UZ_UZ__n + visco_T * Qvisc_UZ_UZ__n
                  Qjac_k (var_UZ,var_UZ )  =                 + visco_T * Qvisc_UZ_UZ__k
                  Qjac_kn(var_UZ,var_UZ )  =                 + visco_T * Qvisc_UZ_UZ__kn

                  Qjac_p (var_UZ,var_Up )  =  Qconv_UZ_Up__p + visco_T * Qvisc_UZ_Up__p
                  Qjac_n (var_UZ,var_Up )  =  Qconv_UZ_Up__n + visco_T * Qvisc_UZ_Up__n
                  Qjac_k (var_UZ,var_Up )  =                 + visco_T * Qvisc_UZ_Up__k
                  Qjac_kn(var_UZ,var_Up )  =                 + visco_T * Qvisc_UZ_Up__kn

                  Qjac_p (var_UZ,var_rho)  =  Qconv_UZ_rho__p + v_Z * (rho*T0) - v * PneoZ_rho__p
                  Qjac_n (var_UZ,var_rho)  =  Qconv_UZ_rho__n                  - v * PneoZ_rho__n

                  Qjac_p (var_UZ,var_T  )  =  Qconv_UZ_T__p   + v_Z * (rho0*T) + dvisco_dT * T * Qvisc_UZ__p - v * PneoZ_T__p
                  Qjac_n (var_UZ,var_T  )  =  Qconv_UZ_T__n                                                  - v * PneoZ_T__n
                  Qjac_k (var_UZ,var_T  )  =                                   + dvisco_dT * T * Qvisc_UZ__k

                  !###################################################################################################
                  !#  equation 6   (Phi component momentum equation)                                                 #
                  !###################################################################################################
                  Pjac(var_Up,var_UR)    = v * BR0 * rho0_corr * UR
                  Pjac(var_Up,var_UZ)    = v * BZ0 * rho0_corr * UZ
                  Pjac(var_Up,var_Up)    = v * Bp0 * rho0_corr * Up 

                  Qjac_p (var_Up,var_AR)  = + BR0_AR * Qconv_UR              &
                                            + Bp0_AR * Qconv_Up              &
                                            + BR0    * Qconv_UR_AR__p        &
                                            + BZ0    * Qconv_UZ_AR__p        &
                                            + Bp0    * Qconv_Up_AR__p        &
                                            + BR0_AR * visco_T * Qvisc_UR__p &
                                            + Bp0_AR * visco_T * Qvisc_Up__p &
                                            + p0 * BgradVstar_AR__p          &
                                            - v * (BR0*PneoR_AR__p + BZ0*PneoZ_AR__p) &
                                            - v * particle_source(ms,mt) * (Bp0_AR*Up0) &
                                            - visco_num * lap_Vstar * (BR0_AR * lap_UR               + Bp0_AR * lap_Up)
                  Qjac_n (var_Up,var_AR)  = + BZ0_AR__n * Qconv_UZ              &
                                            + BR0       * Qconv_UR_AR__n        &
                                            + BZ0       * Qconv_UZ_AR__n        &
                                            + Bp0       * Qconv_Up_AR__n        &
                                            + BZ0_AR__n * visco_T * Qvisc_UZ__p &
                                            + p0 * BgradVstar_AR__n             &
                                            - v * (BR0*PneoR_AR__n + BZ0*PneoZ_AR__n) &
                                            - v * (                + BZ0_AR__n*PneoZ) &
                                            - v * particle_source(ms,mt) * (BZ0_AR__n*UZ0) &
                                            - visco_num * lap_Vstar * (               + BZ0_AR__n * lap_UZ)
                  Qjac_k (var_Up,var_AR)  = + BR0_AR * visco_T * Qvisc_UR__k &
                                            + Bp0_AR * visco_T * Qvisc_Up__k &
                                            + p0 * BgradVstar_AR__k
                  Qjac_kn(var_Up,var_AR)  = + BZ0_AR__n * visco_T * Qvisc_UZ__k

                  Qjac_p (var_Up,var_AZ)  = + BZ0_AZ * Qconv_UZ              &
                                            + Bp0_AZ * Qconv_Up              &
                                            + BR0    * Qconv_UR_AZ__p        &
                                            + BZ0    * Qconv_UZ_AZ__p        &
                                            + Bp0    * Qconv_Up_AZ__p        &
                                            + BZ0_AZ * visco_T * Qvisc_UZ__p &
                                            + Bp0_AZ * visco_T * Qvisc_Up__p &
                                            + p0 * BgradVstar_AZ__p          &
                                            - v * (BR0*PneoR_AZ__p + BZ0*PneoZ_AZ__p) &
                                            - v * particle_source(ms,mt) * (Bp0_AZ*Up0) &
                                            - visco_num * lap_Vstar * (              + BZ0_AZ * lap_UZ + Bp0_AZ * lap_Up)
                  Qjac_n (var_Up,var_AZ)  = + BR0_AZ__n * Qconv_UR              &
                                            + BR0       * Qconv_UR_AZ__n        &
                                            + BZ0       * Qconv_UZ_AZ__n        &
                                            + Bp0       * Qconv_Up_AZ__n        &
                                            + BR0_AZ__n * visco_T * Qvisc_UR__p &
                                            + p0 * BgradVstar_AZ__n             &
                                            - v * (BR0*PneoR_AZ__n + BZ0*PneoZ_AZ__n) &
                                            - v * (BR0_AZ__n*PneoR                  ) &
                                            - v * particle_source(ms,mt) * (BR0_AZ__n*UR0) &
                                            - visco_num * lap_Vstar * (BR0_AZ__n * lap_UR )
                  Qjac_k (var_Up,var_AZ)  = + BZ0_AZ * visco_T * Qvisc_UZ__k &
                                            + Bp0_AZ * visco_T * Qvisc_Up__k &
                                            + p0 * BgradVstar_AZ__k
                  Qjac_kn(var_Up,var_AZ)  = + BR0_AZ__n * visco_T * Qvisc_UR__k

                  Qjac_p (var_Up,var_A3)  = + BR0_A3 * Qconv_UR              &
                                            + BZ0_A3 * Qconv_UZ              &
                                            + Bp0_A3 * Qconv_Up              &
                                            + BR0    * Qconv_UR_A3__p        &
                                            + BZ0    * Qconv_UZ_A3__p        &
                                            + Bp0    * Qconv_Up_A3__p        &
                                            + BR0_A3 * visco_T * Qvisc_UR__p &
                                            + BZ0_A3 * visco_T * Qvisc_UZ__p &
                                            + Bp0_A3 * visco_T * Qvisc_Up__p &
                                            + p0 * BgradVstar_A3__p          &
                                            - v * (BR0*PneoR_A3__p + BZ0*PneoZ_A3__p) &
                                            - v * (BR0_A3*PneoR    + BZ0_A3*PneoZ   ) &
                                            - v * particle_source(ms,mt) * (BR0_A3*UR0 + BZ0_A3*UZ0) &
                                            - visco_num * lap_Vstar * (BR0_A3 * lap_UR + BZ0_A3 * lap_UZ + Bp0_A3 * lap_Up)
                  Qjac_n (var_Up,var_A3)  = + BR0    * Qconv_UR_A3__n        &
                                            + BZ0    * Qconv_UZ_A3__n        &
                                            + Bp0    * Qconv_Up_A3__n        &
                                            - v * (BR0*PneoR_A3__n + BZ0*PneoZ_A3__n)
                  Qjac_k (var_Up,var_A3)  = + BR0_A3 * visco_T * Qvisc_UR__k &
                                            + BZ0_A3 * visco_T * Qvisc_UZ__k &
                                            + Bp0_A3 * visco_T * Qvisc_Up__k &
                                            + p0 * BgradVstar_A3__k

                  Qjac_p (var_Up,var_UR)  = + BR0 * Qconv_UR_UR__p           &
                                            + BZ0 * Qconv_UZ_UR__p           &
                                            + Bp0 * Qconv_Up_UR__p           &
                                            + BR0 * visco_T * Qvisc_UR_UR__p &
                                            + BZ0 * visco_T * Qvisc_UZ_UR__p &
                                            + Bp0 * visco_T * Qvisc_Up_UR__p &
                                            - v * (BR0*PneoR_UR + BZ0*PneoZ_UR) &
                                            - v * particle_source(ms,mt) * (BR0*UR) &
                                            - visco_num * lap_Vstar * BR0 * lap_bf
                  Qjac_n (var_Up,var_UR)  = + BR0 * Qconv_UR_UR__n           &
                                            + BZ0 * Qconv_UZ_UR__n           &
                                            + Bp0 * Qconv_Up_UR__n           &
                                            + BR0 * visco_T * Qvisc_UR_UR__n &
                                            + BZ0 * visco_T * Qvisc_UZ_UR__n &
                                            + Bp0 * visco_T * Qvisc_Up_UR__n
                  Qjac_k (var_Up,var_UR)  = + BR0 * visco_T * Qvisc_UR_UR__k &
                                            + BZ0 * visco_T * Qvisc_UZ_UR__k &
                                            + Bp0 * visco_T * Qvisc_Up_UR__k
                  Qjac_kn(var_Up,var_UR)  = + BR0 * visco_T * Qvisc_UR_UR__kn &
                                            + BZ0 * visco_T * Qvisc_UZ_UR__kn &
                                            + Bp0 * visco_T * Qvisc_Up_UR__kn

                  Qjac_p (var_Up,var_UZ)  = + BR0 * Qconv_UR_UZ__p           &
                                            + BZ0 * Qconv_UZ_UZ__p           &
                                            + Bp0 * Qconv_Up_UZ__p           &
                                            + BR0 * visco_T * Qvisc_UR_UZ__p &
                                            + BZ0 * visco_T * Qvisc_UZ_UZ__p &
                                            + Bp0 * visco_T * Qvisc_Up_UZ__p &
                                            - v * (BR0*PneoR_UZ + BZ0*PneoZ_UZ) &
                                            - v * particle_source(ms,mt) * (BZ0*UZ) &
                                            - visco_num * lap_Vstar * BZ0 * lap_bf
                  Qjac_n (var_Up,var_UZ)  = + BR0 * Qconv_UR_UZ__n           &
                                            + BZ0 * Qconv_UZ_UZ__n           &
                                            + Bp0 * Qconv_Up_UZ__n           &
                                            + BR0 * visco_T * Qvisc_UR_UZ__n &
                                            + BZ0 * visco_T * Qvisc_UZ_UZ__n &
                                            + Bp0 * visco_T * Qvisc_Up_UZ__n
                  Qjac_k (var_Up,var_UZ)  = + BR0 * visco_T * Qvisc_UR_UZ__k &
                                            + BZ0 * visco_T * Qvisc_UZ_UZ__k &
                                            + Bp0 * visco_T * Qvisc_Up_UZ__k
                  Qjac_kn(var_Up,var_UZ)  = + BR0 * visco_T * Qvisc_UR_UZ__kn &
                                            + BZ0 * visco_T * Qvisc_UZ_UZ__kn &
                                            + Bp0 * visco_T * Qvisc_Up_UZ__kn

                  Qjac_p (var_Up,var_Up)  = + BR0 * Qconv_UR_Up__p           &
                                            + BZ0 * Qconv_UZ_Up__p           &
                                            + Bp0 * Qconv_Up_Up__p           &
                                            + BR0 * visco_T * Qvisc_UR_Up__p &
                                            + BZ0 * visco_T * Qvisc_UZ_Up__p &
                                            + Bp0 * visco_T * Qvisc_Up_Up__p &
                                            - v * particle_source(ms,mt) * (Bp0*Up) &
                                            - visco_num * lap_Vstar * Bp0 * lap_bf
                  Qjac_n (var_Up,var_Up)  = + BR0 * Qconv_UR_Up__n              &
                                            + BZ0 * Qconv_UZ_Up__n              &
                                            + Bp0 * Qconv_Up_Up__n              &
                                            + BR0 * visco_T * Qvisc_UR_Up__n &
                                            + BZ0 * visco_T * Qvisc_UZ_Up__n &
                                            + Bp0 * visco_T * Qvisc_Up_Up__n
                  Qjac_k (var_Up,var_Up)  = + BR0 * visco_T * Qvisc_UR_Up__k &
                                            + BZ0 * visco_T * Qvisc_UZ_Up__k &
                                            + Bp0 * visco_T * Qvisc_Up_Up__k
                  Qjac_kn(var_Up,var_Up)  = + BR0 * visco_T * Qvisc_UR_Up__kn &
                                            + BZ0 * visco_T * Qvisc_UZ_Up__kn &
                                            + Bp0 * visco_T * Qvisc_Up_Up__kn

                  Qjac_p (var_Up,var_rho) = + BR0 * Qconv_UR_rho__p              &
                                            + BZ0 * Qconv_UZ_rho__p              &
                                            + Bp0 * Qconv_Up_rho__p              &
                                            + (rho*T0) * BgradVstar__p           &
                                            - v * (BR0*PneoR_rho__p + BZ0*PneoZ_rho__p)
                  Qjac_n (var_Up,var_rho) = + BR0 * Qconv_UR_rho__n              &
                                            + BZ0 * Qconv_UZ_rho__n              &
                                            + Bp0 * Qconv_Up_rho__n              &
                                            - v * (BR0*PneoR_rho__n + BZ0*PneoZ_rho__n)
                  Qjac_k (var_Up,var_rho) = + (rho*T0) * BgradVstar__k

                  Qjac_p (var_Up,var_T )  = + BR0 * Qconv_UR_T__p               &
                                            + BZ0 * Qconv_UZ_T__p               &
                                            + Bp0 * Qconv_Up_T__p               &
                                            + BR0 * dvisco_dT * T * Qvisc_UR__p &
                                            + BZ0 * dvisco_dT * T * Qvisc_UZ__p &
                                            + Bp0 * dvisco_dT * T * Qvisc_Up__p &
                                            + (rho0*T) * BgradVstar__p          &
                                            - v * (BR0*PneoR_T__p + BZ0*PneoZ_T__p)
                  Qjac_n (var_Up,var_T )  = + BR0 * Qconv_UR_T__n               &
                                            + BZ0 * Qconv_UZ_T__n               &
                                            + Bp0 * Qconv_Up_T__n               &
                                            - v * (BR0*PneoR_T__n + BZ0*PneoZ_T__n)
                  Qjac_k (var_Up,var_T )  = + BR0 * dvisco_dT * T * Qvisc_UR__k &
                                            + BZ0 * dvisco_dT * T * Qvisc_UZ__k &
                                            + Bp0 * dvisco_dT * T * Qvisc_Up__k &
                                            + (rho0*T) * BgradVstar__k

                  !###################################################################################################
                  !#  equation 7   (Density equation)                                                                #
                  !###################################################################################################
                  Pjac   (var_rho,var_rho) =   v * rho

                  Qjac_p (var_rho,var_AR)  = + rho0 * VdiaGradVstar_AR__p &
                                             - (D_par-D_prof) * BgradVstar_AR__p * BgradRho       / BB2 &
                                             - (D_par-D_prof) * BgradVstar__p    * BgradRho_AR__p / BB2 &
                                             + (D_par-D_prof) * BgradVstar__p    * BgradRho       / BB2**2 * BB2_AR__p
                  Qjac_n (var_rho,var_AR)  = + rho0 * VdiaGradVstar_AR__n &
                                             - (D_par-D_prof) * BgradVstar_AR__n * BgradRho       / BB2 &
                                             - (D_par-D_prof) * BgradVstar__p    * BgradRho_AR__n / BB2 &
                                             + (D_par-D_prof) * BgradVstar__p    * BgradRho       / BB2**2 * BB2_AR__n
                  Qjac_k (var_rho,var_AR)  = + rho0 * VdiaGradVstar_AR__k &
                                             - (D_par-D_prof) * BgradVstar_AR__k * BgradRho       / BB2 &
                                             - (D_par-D_prof) * BgradVstar__k    * BgradRho_AR__p / BB2 &
                                             + (D_par-D_prof) * BgradVstar__k    * BgradRho       / BB2**2 * BB2_AR__p
                  Qjac_kn(var_rho,var_AR)  = + rho0 * VdiaGradVstar_AR__kn &
                                             - (D_par-D_prof) * BgradVstar__k    * BgradRho_AR__n / BB2 &
                                             + (D_par-D_prof) * BgradVstar__k    * BgradRho       / BB2**2 * BB2_AR__n

                  Qjac_p (var_rho,var_AZ)  = + rho0 * VdiaGradVstar_AZ__p &
                                             - (D_par-D_prof) * BgradVstar_AZ__p * BgradRho       / BB2 &
                                             - (D_par-D_prof) * BgradVstar__p    * BgradRho_AZ__p / BB2 &
                                             + (D_par-D_prof) * BgradVstar__p    * BgradRho       / BB2**2 * BB2_AZ__p
                  Qjac_n (var_rho,var_AZ)  = + rho0 * VdiaGradVstar_AZ__n &
                                             - (D_par-D_prof) * BgradVstar_AZ__n * BgradRho       / BB2 &
                                             - (D_par-D_prof) * BgradVstar__p    * BgradRho_AZ__n / BB2 &
                                             + (D_par-D_prof) * BgradVstar__p    * BgradRho       / BB2**2 * BB2_AZ__n
                  Qjac_k (var_rho,var_AZ)  = + rho0 * VdiaGradVstar_AZ__k &
                                             - (D_par-D_prof) * BgradVstar_AZ__k * BgradRho       / BB2 &
                                             - (D_par-D_prof) * BgradVstar__k    * BgradRho_AZ__p / BB2 &
                                             + (D_par-D_prof) * BgradVstar__k    * BgradRho       / BB2**2 * BB2_AZ__p
                  Qjac_kn(var_rho,var_AZ)  = + rho0 * VdiaGradVstar_AZ__kn &
                                             - (D_par-D_prof) * BgradVstar__k    * BgradRho_AZ__n / BB2 &
                                             + (D_par-D_prof) * BgradVstar__k    * BgradRho       / BB2**2 * BB2_AZ__n

                  Qjac_p (var_rho,var_A3)  = + rho0 * VdiaGradVstar_A3__p &
                                             - (D_par-D_prof) * BgradVstar_A3__p * BgradRho    / BB2 &
                                             - (D_par-D_prof) * BgradVstar__p    * BgradRho_A3 / BB2 &
                                             + (D_par-D_prof) * BgradVstar__p    * BgradRho    / BB2**2 * BB2_A3
                  Qjac_n (var_rho,var_A3)  = + rho0 * VdiaGradVstar_A3__n
                  Qjac_k (var_rho,var_A3)  = + rho0 * VdiaGradVstar_A3__k &
                                             - (D_par-D_prof) * BgradVstar_A3__k * BgradRho    / BB2 &
                                             - (D_par-D_prof) * BgradVstar__k    * BgradRho_A3 / BB2 &
                                             + (D_par-D_prof) * BgradVstar__k    * BgradRho    / BB2**2 * BB2_A3
                  Qjac_kn(var_rho,var_A3)  = + rho0 * VdiaGradVstar_A3__kn

                  Qjac_p (var_rho,var_UR)  = - v * ( rho0 * divU_UR + UgradRho_UR )

                  Qjac_p (var_rho,var_UZ)  = - v * ( rho0 * divU_UZ + UgradRho_UZ )

                  Qjac_p (var_rho,var_Up)  = - v * (                + UgradRho_Up )
                  Qjac_n (var_rho,var_Up)  = - v * ( rho0 * divU_Up__n )

                  Qjac_p (var_rho,var_rho) = - v * ( rho * divU + UgradRho_rho__p )                   &
                                             + rho  * VdiaGradVstar__p &
                                             + rho0 * VdiaGradVstar_rho__p &
                                             - D_prof * gradRho_gradVstar_rho__p                      &
                                             - (D_par-D_prof) * BgradVstar__p * BgradRho_rho__p / BB2 &
                                             - D_perp_num * lap_Vstar * lap_bf
                  Qjac_n (var_rho,var_rho) = - v * (UgradRho_rho__n )                                 &
                                             + rho0 * VdiaGradVstar_rho__n &
                                             - (D_par-D_prof) * BgradVstar__p * BgradRho_rho__n / BB2
                  Qjac_k (var_rho,var_rho) = + rho0 * VdiaGradVstar_rho__k &
                                             - (D_par-D_prof) * BgradVstar__k * BgradRho_rho__p / BB2
                  Qjac_kn(var_rho,var_rho) = + rho0 * VdiaGradVstar_rho__kn &
                                             - D_prof * gradRho_gradVstar_rho__kn                     &
                                             - (D_par-D_prof) * BgradVstar__k * BgradRho_rho__n / BB2

                  Qjac_p (var_rho,var_T  ) = + rho0 * VdiaGradVstar_T__p
                  Qjac_n (var_rho,var_T  ) = + rho0 * VdiaGradVstar_T__n
                  Qjac_k (var_rho,var_T  ) = + rho0 * VdiaGradVstar_T__k
                  Qjac_kn(var_rho,var_T  ) = + rho0 * VdiaGradVstar_T__kn

                  !###################################################################################################
                  !#  equation 8   (Temperature  equation)                                                           #
                  !###################################################################################################
                  Pjac   (var_T,var_T)   =   v * rho0_corr * T
                  Pjac   (var_T,var_rho) =   v * rho       * T0_corr

                  Qjac_p (var_T,var_AR)  = - (ZKpar_T-ZK_prof) * BgradVstar_AR__p * BgradT       / BB2                &
                                           - (ZKpar_T-ZK_prof) * BgradVstar__p    * BgradT_AR__p / BB2                &
                                           + (ZKpar_T-ZK_prof) * BgradVstar__p    * BgradT       / BB2**2 * BB2_AR__p &
                                           + v * (gamma-1.0d0) * eta_T_ohm * JJ2_AR__p
                  Qjac_n (var_T,var_AR)  = - (ZKpar_T-ZK_prof) * BgradVstar_AR__n * BgradT       / BB2                &
                                           - (ZKpar_T-ZK_prof) * BgradVstar__p    * BgradT_AR__n / BB2                &
                                           + (ZKpar_T-ZK_prof) * BgradVstar__p    * BgradT       / BB2**2 * BB2_AR__n &
                                           + v * (gamma-1.0d0) * eta_T_ohm * JJ2_AR__n
                  Qjac_k (var_T,var_AR)  = - (ZKpar_T-ZK_prof) * BgradVstar_AR__k * BgradT       / BB2                &
                                           - (ZKpar_T-ZK_prof) * BgradVstar__k    * BgradT_AR__p / BB2                &
                                           + (ZKpar_T-ZK_prof) * BgradVstar__k    * BgradT       / BB2**2 * BB2_AR__p
                  Qjac_kn(var_T,var_AR)  = - (ZKpar_T-ZK_prof) * BgradVstar__k    * BgradT_AR__n / BB2                &
                                           + (ZKpar_T-ZK_prof) * BgradVstar__k    * BgradT       / BB2**2 * BB2_AR__n
                  Qjac_pnn(var_T,var_AR)  = v * (gamma-1.0d0) * eta_T_ohm * JJ2_AR__nn

                  Qjac_p (var_T,var_AZ)  = - (ZKpar_T-ZK_prof) * BgradVstar_AZ__p * BgradT       / BB2                &
                                           - (ZKpar_T-ZK_prof) * BgradVstar__p    * BgradT_AZ__p / BB2                &
                                           + (ZKpar_T-ZK_prof) * BgradVstar__p    * BgradT       / BB2**2 * BB2_AZ__p &
                                           + v * (gamma-1.0d0) * eta_T_ohm * JJ2_AZ__p
                  Qjac_n (var_T,var_AZ)  = - (ZKpar_T-ZK_prof) * BgradVstar_AZ__n * BgradT       / BB2                &
                                           - (ZKpar_T-ZK_prof) * BgradVstar__p    * BgradT_AZ__n / BB2                &
                                           + (ZKpar_T-ZK_prof) * BgradVstar__p    * BgradT       / BB2**2 * BB2_AZ__n &
                                           + v * (gamma-1.0d0) * eta_T_ohm * JJ2_AZ__n
                  Qjac_k (var_T,var_AZ)  = - (ZKpar_T-ZK_prof) * BgradVstar_AZ__k * BgradT       / BB2                &
                                           - (ZKpar_T-ZK_prof) * BgradVstar__k    * BgradT_AZ__p / BB2                &
                                           + (ZKpar_T-ZK_prof) * BgradVstar__k    * BgradT       / BB2**2 * BB2_AZ__p
                  Qjac_kn(var_T,var_AZ)  = - (ZKpar_T-ZK_prof) * BgradVstar__k    * BgradT_AZ__n / BB2                &
                                           + (ZKpar_T-ZK_prof) * BgradVstar__k    * BgradT       / BB2**2 * BB2_AZ__n
                  Qjac_pnn(var_T,var_AZ)  = v * (gamma-1.0d0) * eta_T_ohm * JJ2_AZ__nn

                  Qjac_p (var_T,var_A3)  = - (ZKpar_T-ZK_prof) * BgradVstar_A3__p * BgradT    / BB2             &
                                           - (ZKpar_T-ZK_prof) * BgradVstar__p    * BgradT_A3 / BB2             &
                                           + (ZKpar_T-ZK_prof) * BgradVstar__p    * BgradT    / BB2**2 * BB2_A3 &
                                           + v * (gamma-1.0d0) * eta_T_ohm * JJ2_A3__p
                  Qjac_n (var_T,var_A3)  = + v * (gamma-1.0d0) * eta_T_ohm * JJ2_A3__n
                  Qjac_k (var_T,var_A3)  = - (ZKpar_T-ZK_prof) * BgradVstar_A3__k * BgradT    / BB2             &
                                           - (ZKpar_T-ZK_prof) * BgradVstar__k    * BgradT_A3 / BB2             &
                                           + (ZKpar_T-ZK_prof) * BgradVstar__k    * BgradT    / BB2**2 * BB2_A3
                  Qjac_pnn(var_T,var_A3)  = v * (gamma-1.0d0) * eta_T_ohm * JJ2_A3__nn

                  Qjac_p (var_T,var_UR)  = + v * ( - rho0 * UgradT_UR  -  T0 * UgradRho_UR  -  gamma * p0 * divU_UR ) &
                                           + v * (gamma-1.d0) * Qvisc_T_UR__p
                  Qjac_n (var_T,var_UR)  = + v * (gamma-1.d0) * Qvisc_T_UR__n

                  Qjac_p (var_T,var_UZ)  = + v * ( - rho0 * UgradT_UZ  -  T0 * UgradRho_UZ  -  gamma * p0 * divU_UZ ) &
                                           + v * (gamma-1.d0) * Qvisc_T_UZ__p
                  Qjac_n (var_T,var_UZ)  = + v * (gamma-1.d0) * Qvisc_T_UZ__n

                  Qjac_p (var_T,var_Up)  = + v * ( - rho0 * UgradT_Up  -  T0 * UgradRho_Up                             ) &
                                           + v * (gamma-1.d0) * Qvisc_T_Up__p
                  Qjac_n (var_T,var_Up)  = + v * (                                          -  gamma * p0 * divU_Up__n ) &
                                           + v * (gamma-1.d0) * Qvisc_T_Up__n

                  Qjac_p (var_T,var_rho) = + v * ( - rho * UgradT  -  T0 * UgradRho_rho__p  -  gamma * (rho*T0) * divU )
                  Qjac_n (var_T,var_rho) = + v * (                 -  T0 * UgradRho_rho__n                             )

                  Qjac_p (var_T,var_T )  = + v * ( - rho0 * UgradT_T__p  -  T * UgradRho  -  gamma * (rho0*T) * divU ) &
                                           + v * (gamma-1.d0) * Qvisc_T_T__p                                           &
                                           - ZK_prof * gradT_gradVstar_T__p                                            &
                                           - (ZKpar_T-ZK_prof) * BgradVstar__p * BgradT_T__p / BB2                     &
                                           - (dZKpar_dT*T    ) * BgradVstar__p * BgradT      / BB2                     &
                                           - ZK_perp_num * lap_Vstar * lap_bf                                          &
                                           + v * (gamma-1.0d0) * deta_dT_ohm * T * JJ2 
                  Qjac_n (var_T,var_T )  = + v * ( - rho0 * UgradT_T__n                                              ) &
                                           + v * (gamma-1.d0) * Qvisc_T_T__n                                           &
                                           - (ZKpar_T-ZK_prof) * BgradVstar__p * BgradT_T__n / BB2
                  Qjac_k (var_T,var_T )  = - (ZKpar_T-ZK_prof) * BgradVstar__k * BgradT_T__p / BB2                     &
                                           - (dZKpar_dT*T    ) * BgradVstar__k * BgradT      / BB2
                  Qjac_kn(var_T,var_T )  = - ZK_prof * gradT_gradVstar_T__kn                                           &
                                           - (ZKpar_T-ZK_prof) * BgradVstar__k * BgradT_T__n / BB2

                  !###################################################################################################
                  !#  VMS STABILISATION                                                                              #
                  !###################################################################################################

                  SELECT CASE(VmsType)
                  CASE(10)
                     
                     !#  equations 1,2,3                                                                                #
                     QvmsAd_p (var_AR,var_AR) = CvGradVj__p * CvGradVi__p + Cvp0 * AR  * Cvp0 * v / R**2
                     QvmsAd_n (var_AR,var_AR) = CvGradVj__n * CvGradVi__p
                     QvmsAd_k (var_AR,var_AR) = CvGradVj__p * CvGradVi__k
                     QvmsAd_kn(var_AR,var_AR) = CvGradVj__n * CvGradVi__k

                     QvmsAd_p (var_AR,var_A3) = Cvp0 * A3 / R**2  * CvGradVi__p + ( - CvGradVj__p + CvR0 * A3 / R ) * Cvp0 * v / R**2
                     QvmsAd_n (var_AR,var_A3) =                                 + ( - CvGradVj__n + CvR0 * A3 / R ) * Cvp0 * v / R**2
                     QvmsAd_k (var_AR,var_A3) = Cvp0 * A3 / R**2  * CvGradVi__k
                     QvmsAd_kn(var_AR,var_A3) = 0.d0

                     QvmsAd_p (var_AZ,var_AZ) = CvGradVj__p * CvGradVi__p
                     QvmsAd_n (var_AZ,var_AZ) = CvGradVj__n * CvGradVi__p
                     QvmsAd_k (var_AZ,var_AZ) = CvGradVj__p * CvGradVi__k
                     QvmsAd_kn(var_AZ,var_AZ) = CvGradVj__n * CvGradVi__k

                     QvmsAd_p (var_A3,var_A3) = (Cvp0 * A3 / R**2 ) * Cvp0 * v                                    &
                                                - ( - CvGradVj__p + CvR0 * A3 / R ) * (CvGradVi__p + CvR0 * v / R)
                     QvmsAd_n (var_A3,var_A3) = - ( - CvGradVj__n                 ) * (CvGradVi__p + CvR0 * v / R)
                     QvmsAd_k (var_A3,var_A3) = - ( - CvGradVj__p + CvR0 * A3 / R ) * (CvGradVi__k                  )
                     QvmsAd_kn(var_A3,var_A3) = - ( - CvGradVj__n                 ) * (CvGradVi__k                  )

                     QvmsAd_p (var_A3,var_AR) =    CvGradVj__p * Cvp0 * v                       &
                                                - (  Cvp0 * AR ) * (CvGradVi__p + CvR0 * v / R)
                     QvmsAd_n (var_A3,var_AR) =    CvGradVj__n * Cvp0 * v
                     QvmsAd_k (var_A3,var_AR) = - (  Cvp0 * AR ) * (CvGradVi__k               )
                     
                     !#  equation 4   (R component momentum equation)                                                   #
                     QvmsAd_p (var_UR,var_rho ) =  0.0  ! to be completed  !$BNK
                     QvmsAd_n (var_UR,var_rho ) =  0.0  ! to be completed  !$BNK
                     QvmsAd_k (var_UR,var_rho ) =  0.0  ! to be completed  !$BNK
                     QvmsAd_kn(var_UR,var_rho ) =  0.0  ! to be completed  !$BNK

                     QvmsAd_p (var_UR,var_UR) =   ( rho0 * CvGradVj__p + CvGradr0 * UR  ) *   CvGradVi__p      &
                                                + ( rho0 * Cvp0 * UR / R                ) * ( Cvp0 * v /R )    &
                                                + ( rho0 * CvGradVj__p + CvGradr0 * UR  ) *   CvGradVi__p      &
                                                + ( rho0 * Cvp0 * UR / R                ) * ( Cvp0 * v /R )
                     QvmsAd_n (var_UR,var_UR) =   ( rho0 * CvGradVj__n                  ) *   CvGradVi__p      &
                                                + ( rho0 * CvGradVj__n                  ) *   CvGradVi__p
                     QvmsAd_k (var_UR,var_UR) =   ( rho0 * CvGradVj__p + CvGradr0 * UR  ) *   CvGradVi__k      &
                                                + ( rho0 * CvGradVj__p + CvGradr0 * UR  ) *   CvGradVi__k
                     QvmsAd_kn(var_UR,var_UR) =   ( rho0 * CvGradVj__n                  ) *   CvGradVi__k      &
                                                + ( rho0 * CvGradVj__n                  ) *   CvGradVi__k

                     QvmsAd_p (var_UR,var_Up) =   ( - rho0 * Cvp0 * Up / R                ) *   CvGradVi__p   &
                                                + (   rho0 * CvGradVj__p + CvGradr0 * Up  ) * ( Cvp0 * v /R ) &
                                                + ( - rho0 * Vbp0 * Up / R                ) *   VbGradVi__p   &
                                                + (   rho0 * VbGradVj__p + VbGradr0 * Up  ) * ( Vbp0 * v /R )
                     QvmsAd_n (var_UR,var_Up) = + (   rho0 * CvGradVj__n                  ) * ( Cvp0 * v /R ) &
                                                + (   rho0 * VbGradVj__n                  ) * ( Vbp0 * v /R )
                     QvmsAd_k (var_UR,var_Up) =   ( - rho0 * Cvp0 * Up / R                ) *   CvGradVi__k &
                                                + ( - rho0 * Vbp0 * Up / R                ) *   VbGradVi__k

                     QvmsF_p (var_UR,var_rho)  =  VmsCoefF * ( UgradBF__p     + rho * divU           ) * DiveRMVi
                     QvmsF_n (var_UR,var_rho)  =  VmsCoefF * ( UgradBF__n                            ) * DiveRMVi

                     QvmsF_p (var_UR,var_UR )  =  VmsCoefF * ( rho0 * DiveRMVj + UR * rho0_R         ) * DiveRMVi

                     QvmsF_p (var_UR,var_UZ )  =  VmsCoefF * ( rho0 * DiveZMVj + UZ * rho0_Z         ) * DiveRMVi

                     QvmsF_p (var_UR,var_Up )  =  VmsCoefF * (                    + Up * rho0_p/R    ) * DiveRMVi
                     QvmsF_n (var_UR,var_Up )  =  VmsCoefF * ( rho0 * DivePMVj__n                    ) * DiveRMVi

                     !#  equation 5   (Z component momentum equation)                                                   #
                     QvmsAd_p (var_UZ,var_rho ) =  0.0 ! to be completed  !$BNK
                     QvmsAd_n (var_UZ,var_rho ) =  0.0 ! to be completed  !$BNK
                     QvmsAd_k (var_UZ,var_rho ) =  0.0 ! to be completed  !$BNK
                     QvmsAd_kn(var_UZ,var_rho ) =  0.0 ! to be completed  !$BNK

                     QvmsAd_p (var_UZ,var_UZ) =   ( rho0 * CvGradVj__p + CvGradr0 * UZ )  * CvGradVi__p  &
                                                + ( rho0 * VbGradVj__p + VbGradr0 * UZ )  * VbGradVi__p
                     QvmsAd_n (var_UZ,var_UZ) =   ( rho0 * CvGradVj__n                 )  * CvGradVi__p  &
                                                + ( rho0 * VbGradVj__n                 )  * VbGradVi__p
                     QvmsAd_k (var_UZ,var_UZ) =   ( rho0 * CvGradVj__p + CvGradr0 * UZ )  * CvGradVi__k  &
                                                + ( rho0 * VbGradVj__p + VbGradr0 * UZ )  * VbGradVi__k
                     QvmsAd_kn(var_UZ,var_UZ) =   ( rho0 * CvGradVj__n                 )  * CvGradVi__k  &
                                                + ( rho0 * VbGradVj__n                 )  * VbGradVi__k

                     QvmsF_p (var_UZ,var_rho)  =  VmsCoefF * ( UgradBF__p     + rho * divU          ) * DiveZMVi
                     QvmsF_n (var_UZ,var_rho)  =  VmsCoefF * ( UgradBF__n                           ) * DiveZMVi

                     QvmsF_p (var_UZ,var_UR )  =  VmsCoefF * ( rho0 * DiveRMVj    + UR * rho0_R     ) * DiveZMVi

                     QvmsF_p (var_UZ,var_UZ )  =  VmsCoefF * ( rho0 * DiveZMVj    + UZ * rho0_Z     ) * DiveZMVi

                     QvmsF_p (var_UZ,var_Up )  =  VmsCoefF * (                    + Up * rho0_p/R   ) * DiveZMVi
                     QvmsF_n (var_UZ,var_Up )  =  VmsCoefF * ( rho0 * DivePMVj__n                   ) * DiveZMVi

                     !#  equation 6   (Phi component momentum equation)                                                 #
                     QvmsAd_p (var_Up,var_rho ) = 0.0  ! to be completed  !$BNK
                     QvmsAd_n (var_Up,var_rho ) = 0.0  ! to be completed  !$BNK
                     QvmsAd_k (var_Up,var_rho ) = 0.0  ! to be completed  !$BNK
                     QvmsAd_kn(var_Up,var_rho ) = 0.0  ! to be completed  !$BNK

                     QvmsAd_p (var_Up,var_UR) =  -  ( rho0 * CvGradVj__p + CvGradr0 * UR ) *   Cvp0 * v / R                 &
                                                 +  ( rho0 * Cvp0 * UR / R               ) * ( CvGradVi__p + CvR0 * v /R )  &
                                                 -  ( rho0 * VbGradVj__p + VbGradr0 * UR ) *   Vbp0 * v / R                 &
                                                 +  ( rho0 * Vbp0 * UR / R               ) * ( VbGradVi__p + VbR0 * v /R )  
                     QvmsAd_n (var_Up,var_UR) =  -  ( rho0 * CvGradVj__n                 ) *   Cvp0 * v / R                 &
                                                 -  ( rho0 * VbGradVj__n + VbGradr0 * UR ) *   Vbp0 * v / R
                     QvmsAd_k (var_Up,var_UR) =  +  ( rho0 * Cvp0 * UR / R               ) * ( CvGradVi__k + CvR0 * v /R )  &
                                                 +  ( rho0 * Vbp0 * UR / R               ) * ( VbGradVi__k + VbR0 * v /R )  

                     QvmsAd_p (var_Up,var_Up) =  - ( - rho0 * Cvp0 * Up / R                ) *   Cvp0 * v / R                    &
                                                 + (   rho0 * CvGradVj__p + Up * CvGradr0  ) * ( CvGradVi__p + CvR0 * v /R    )  &
                                                 - ( - rho0 * Vbp0 * Up / R                ) *   Vbp0 * v / R                    &
                                                 + (   rho0 * VbGradVj__p + Up * VbGradr0  ) * ( VbGradVi__p + VbR0 * v /R    )
                     QvmsAd_n (var_Up,var_Up) =  + (   rho0 * CvGradVj__n                  ) * ( CvGradVi__p + CvR0 * v /R    )  &
                                                 + (   rho0 * VbGradVj__n                  ) * ( VbGradVi__p + VbR0 * v /R    )
                     QvmsAd_k (var_Up,var_Up) =  + (   rho0 * CvGradVj__p + Up * CvGradr0  ) * ( CvGradVi__k                  )  &
                                                 + (   rho0 * VbGradVj__p + Up * VbGradr0  ) * ( VbGradVi__k                  )
                     QvmsAd_kn(var_Up,var_Up) =  + (   rho0 * CvGradVj__n                  ) * ( CvGradVi__k                  )  &
                                                 + (   rho0 * VbGradVj__n                  ) * ( VbGradVi__k                  )

                     QvmsF_k (var_Up,var_rho)  =  VmsCoefF * ( UgradBF__p     + rho * divU           ) * DivePMVi__k
                     QvmsF_kn(var_Up,var_rho)  =  VmsCoefF * ( UgradBF__n                            ) * DivePMVi__k

                     QvmsF_k (var_Up,var_UR )  =  VmsCoefF * ( rho0 * DiveRMVj    + UR * rho0_R      ) * DivePMVi__k

                     QvmsF_k (var_Up,var_UZ )  =  VmsCoefF * ( rho0 * DiveZMVj    + UZ * rho0_Z      ) * DivePMVi__k

                     QvmsF_k (var_Up,var_Up )  =  VmsCoefF * (                    + Up * rho0_p/R    ) * DivePMVi__k
                     QvmsF_kn(var_Up,var_Up )  =  VmsCoefF * ( rho0 * DivePMVj__n                    ) * DivePMVi__k

                     !#  equation 7   (Density equation)                                                                #
                     QvmsAd_p (var_rho,var_rho) =  CvGradVj__p * CvGradVi__p  + VbGradVi__p * VbGradVj__p
                     QvmsAd_n (var_rho,var_rho) =  CvGradVj__n * CvGradVi__p  + VbGradVi__p * VbGradVj__n
                     QvmsAd_k (var_rho,var_rho) =  CvGradVj__p * CvGradVi__k  + VbGradVi__k * VbGradVj__p
                     QvmsAd_kn(var_rho,var_rho) =  CvGradVj__n * CvGradVi__k  + VbGradVi__k * VbGradVj__n

                     !#  equation 8   (Temperature  equation)                                                           #
                     QvmsAd_p (var_T,var_rho) =   ( T0 * CvGradVj__p  + rho * CvGradT0 ) * CvGradVi__p   &
                                                + ( T0 * VbGradVj__p  + rho * VbGradT0 ) * VbGradVi__p   
                     QvmsAd_n (var_T,var_rho) =   ( T0 * CvGradVj__n                   ) * CvGradVi__p   &
                                                + ( T0 * VbGradVj__n                   ) * VbGradVi__p   
                     QvmsAd_k (var_T,var_rho) =   ( T0 * CvGradVj__p  + rho * CvGradT0 ) * CvGradVi__k   &
                                                + ( T0 * VbGradVj__p  + rho * VbGradT0 ) * VbGradVi__k   
                     QvmsAd_kn(var_T,var_rho) =   ( T0 * CvGradVj__n                   ) * CvGradVi__k   &
                                                + ( T0 * VbGradVj__n                   ) * VbGradVi__k   

                     QvmsAd_p (var_T,var_T) =   ( rho0 * CvGradVj__p  + T * CvGradr0 ) * CvGradVi__p   &
                                              + ( rho0 * VbGradVj__p  + T * VbGradr0 ) * VbGradVi__p   
                     QvmsAd_n (var_T,var_T) =   ( rho0 * CvGradVj__n                 ) * CvGradVi__p   &
                                              + ( rho0 * VbGradVj__n                 ) * VbGradVi__p   
                     QvmsAd_k (var_T,var_T) =   ( rho0 * CvGradVj__p  + T * CvGradr0 ) * CvGradVi__k   &
                                              + ( rho0 * VbGradVj__p  + T * VbGradr0 ) * VbGradVi__k   
                     QvmsAd_kn(var_T,var_T) =   ( rho0 * CvGradVj__n                 ) * CvGradVi__k   &
                                              + ( rho0 * VbGradVj__n                 ) * VbGradVi__k   

                     QvmsF_p (var_T,var_rho)  =  VmsCoefF_T * ( T0*gradBF_gradVstar__p  + gradT_gradVstar__p * rho )
                     QvmsF_k (var_T,var_rho)  =  VmsCoefF_T * (                         + gradT_gradVstar__k * rho )
                     QvmsF_kn(var_T,var_rho)  =  VmsCoefF_T * ( T0*gradBF_gradVstar__kn                            )

                     QvmsF_p (var_T,var_T  )  =  VmsCoefF_T * ( rho0*gradBF_gradVstar__p  + gradRho_gradVstar__p * T )
                     QvmsF_k (var_T,var_T  )  =  VmsCoefF_T * (                           + gradRho_gradVstar__k * T )
                     QvmsF_kn(var_T,var_T  )  =  VmsCoefF_T * ( rho0*gradBF_gradVstar__kn                            )
                     
                     ! --- Add VMS to LHS
                     Qjac_p (var_AR,:) = Qjac_p (var_AR,:) - TG_NUM(var_AR) *  CoefAdv * QvmsAd_p (var_AR,:)
                     Qjac_n (var_AR,:) = Qjac_n (var_AR,:) - TG_NUM(var_AR) *  CoefAdv * QvmsAd_n (var_AR,:)
                     Qjac_k (var_AR,:) = Qjac_k (var_AR,:) - TG_NUM(var_AR) *  CoefAdv * QvmsAd_k (var_AR,:)
                     Qjac_kn(var_AR,:) = Qjac_kn(var_AR,:) - TG_NUM(var_AR) *  CoefAdv * QvmsAd_kn(var_AR,:)

                     Qjac_p (var_AZ,:) = Qjac_p (var_AZ,:) - TG_NUM(var_AZ) *  CoefAdv * QvmsAd_p (var_AZ,:)
                     Qjac_n (var_AZ,:) = Qjac_n (var_AZ,:) - TG_NUM(var_AZ) *  CoefAdv * QvmsAd_n (var_AZ,:)
                     Qjac_k (var_AZ,:) = Qjac_k (var_AZ,:) - TG_NUM(var_AZ) *  CoefAdv * QvmsAd_k (var_AZ,:)
                     Qjac_kn(var_AZ,:) = Qjac_kn(var_AZ,:) - TG_NUM(var_AZ) *  CoefAdv * QvmsAd_kn(var_AZ,:)

                     Qjac_p (var_A3,:) = Qjac_p (var_A3,:) - TG_NUM(var_A3) *  CoefAdv * QvmsAd_p (var_A3,:)
                     Qjac_n (var_A3,:) = Qjac_n (var_A3,:) - TG_NUM(var_A3) *  CoefAdv * QvmsAd_n (var_A3,:)
                     Qjac_k (var_A3,:) = Qjac_k (var_A3,:) - TG_NUM(var_A3) *  CoefAdv * QvmsAd_k (var_A3,:)
                     Qjac_kn(var_A3,:) = Qjac_kn(var_A3,:) - TG_NUM(var_A3) *  CoefAdv * QvmsAd_kn(var_A3,:)

                     Qjac_p (var_UR,:) = Qjac_p (var_UR,:) - TG_NUM(var_UR) * (CoefAdv * QvmsAd_p (var_UR,:) + QvmsF_p (var_UR,:) )
                     Qjac_n (var_UR,:) = Qjac_n (var_UR,:) - TG_NUM(var_UR) * (CoefAdv * QvmsAd_n (var_UR,:) + QvmsF_n (var_UR,:) )
                     Qjac_k (var_UR,:) = Qjac_k (var_UR,:) - TG_NUM(var_UR) * (CoefAdv * QvmsAd_k (var_UR,:) + QvmsF_k (var_UR,:) )
                     Qjac_kn(var_UR,:) = Qjac_kn(var_UR,:) - TG_NUM(var_UR) * (CoefAdv * QvmsAd_kn(var_UR,:) + QvmsF_kn(var_UR,:) )

                     Qjac_p (var_UZ,:) = Qjac_p (var_UZ,:) - TG_NUM(var_UZ) * (CoefAdv * QvmsAd_p (var_UZ,:) + QvmsF_p (var_UZ,:) )
                     Qjac_n (var_UZ,:) = Qjac_n (var_UZ,:) - TG_NUM(var_UZ) * (CoefAdv * QvmsAd_n (var_UZ,:) + QvmsF_n (var_UZ,:) )
                     Qjac_k (var_UZ,:) = Qjac_k (var_UZ,:) - TG_NUM(var_UZ) * (CoefAdv * QvmsAd_k (var_UZ,:) + QvmsF_k (var_UZ,:) )
                     Qjac_kn(var_UZ,:) = Qjac_kn(var_UZ,:) - TG_NUM(var_UZ) * (CoefAdv * QvmsAd_kn(var_UZ,:) + QvmsF_kn(var_UZ,:) )

                     Qjac_p (var_Up, :) = Qjac_p (var_Up, :) - TG_NUM(var_Up) * (  BR0 *( CoefAdv * QvmsAd_p (var_UR,:) + QvmsF_p (var_UR, :) ) &
                                                                                 + BZ0 *( CoefAdv * QvmsAd_p (var_UZ,:) + QvmsF_p (var_UZ, :) ) &
                                                                                 + Bp0 *( CoefAdv * QvmsAd_p (var_Up,:) + QvmsF_p (var_Up, :) ) )
                     Qjac_n (var_Up, :) = Qjac_n (var_Up, :) - TG_NUM(var_Up) * (  BR0 *( CoefAdv * QvmsAd_n (var_UR,:) + QvmsF_n (var_UR, :) ) &
                                                                                 + BZ0 *( CoefAdv * QvmsAd_n (var_UZ,:) + QvmsF_n (var_UZ, :) ) &
                                                                                 + Bp0 *( CoefAdv * QvmsAd_n (var_Up,:) + QvmsF_n (var_Up, :) ) )
                     Qjac_k (var_Up, :) = Qjac_k (var_Up, :) - TG_NUM(var_Up) * (  BR0 *( CoefAdv * QvmsAd_k (var_UR,:) + QvmsF_k (var_UR, :) ) &
                                                                                 + BZ0 *( CoefAdv * QvmsAd_k (var_UZ,:) + QvmsF_k (var_UZ, :) ) &
                                                                                 + Bp0 *( CoefAdv * QvmsAd_k (var_Up,:) + QvmsF_k (var_Up, :) ) )
                     Qjac_kn(var_Up, :) = Qjac_kn(var_Up, :) - TG_NUM(var_Up) * (  BR0 *( CoefAdv * QvmsAd_kn(var_UR,:) + QvmsF_kn(var_UR, :) ) &
                                                                                 + BZ0 *( CoefAdv * QvmsAd_kn(var_UZ,:) + QvmsF_kn(var_UZ, :) ) &
                                                                                 + Bp0 *( CoefAdv * QvmsAd_kn(var_Up,:) + QvmsF_kn(var_Up, :) ) )

                     Qjac_p (var_rho, :) = Qjac_p (var_rho, :) - TG_NUM(var_rho ) *  CoefAdv * QvmsAd_p (var_rho, :) 
                     Qjac_n (var_rho, :) = Qjac_n (var_rho, :) - TG_NUM(var_rho ) *  CoefAdv * QvmsAd_n (var_rho, :) 
                     Qjac_k (var_rho, :) = Qjac_k (var_rho, :) - TG_NUM(var_rho ) *  CoefAdv * QvmsAd_k (var_rho, :) 
                     Qjac_kn(var_rho, :) = Qjac_kn(var_rho, :) - TG_NUM(var_rho ) *  CoefAdv * QvmsAd_kn(var_rho, :) 

                     Qjac_p (var_T, :) = Qjac_p (var_T, :) - TG_NUM(var_T ) * (CoefAdv * QvmsAd_p (var_T, :) + QvmsF_p (var_T, :) )
                     Qjac_n (var_T, :) = Qjac_n (var_T, :) - TG_NUM(var_T ) * (CoefAdv * QvmsAd_n (var_T, :) + QvmsF_n (var_T, :) )
                     Qjac_k (var_T, :) = Qjac_k (var_T, :) - TG_NUM(var_T ) * (CoefAdv * QvmsAd_k (var_T, :) + QvmsF_k (var_T, :) )
                     Qjac_kn(var_T, :) = Qjac_kn(var_T, :) - TG_NUM(var_T ) * (CoefAdv * QvmsAd_kn(var_T, :) + QvmsF_kn(var_T, :) )

                  CASE(0)
                     QvmsF_p (var_UR,var_UR )  =  VmsCoefF * ( rho0 * DiveRMVj     ) * DiveRMVi
                     QvmsF_p (var_UR,var_UZ )  =  VmsCoefF * ( rho0 * DiveZMVj     ) * DiveRMVi
                     QvmsF_n (var_UR,var_Up )  =  VmsCoefF * ( rho0 * DivePMVj__n  ) * DiveRMVi
                    
                     QvmsF_p (var_UZ,var_UR )  =  VmsCoefF * ( rho0 * DiveRMVj     ) * DiveZMVi
                     QvmsF_p (var_UZ,var_UZ )  =  VmsCoefF * ( rho0 * DiveZMVj     ) * DiveZMVi
                     QvmsF_n (var_UZ,var_Up )  =  VmsCoefF * ( rho0 * DivePMVj__n  ) * DiveZMVi
                    
                     QvmsF_k (var_Up,var_UR )  =  VmsCoefF * ( rho0 * DiveRMVj     ) * DivePMVi__k
                     QvmsF_k (var_Up,var_UZ )  =  VmsCoefF * ( rho0 * DiveZMVj     ) * DivePMVi__k
                     QvmsF_kn(var_Up,var_Up )  =  VmsCoefF * ( rho0 * DivePMVj__n  ) * DivePMVi__k
                     
                     ! --- Add VMS to LHS
                     Qjac_p (var_UR, :)   = Qjac_p (var_UR, :) - TG_NUM(var_UR) * QvmsF_p (var_UR, :)
                     Qjac_n (var_UR, :)   = Qjac_n (var_UR, :) - TG_NUM(var_UR) * QvmsF_n (var_UR, :)
                     Qjac_k (var_UR, :)   = Qjac_k (var_UR, :) - TG_NUM(var_UR) * QvmsF_k (var_UR, :)
                     Qjac_kn(var_UR, :)   = Qjac_kn(var_UR, :) - TG_NUM(var_UR) * QvmsF_kn(var_UR, :)

                     Qjac_p (var_UZ, :)   = Qjac_p (var_UZ, :) - TG_NUM(var_UZ) * QvmsF_p (var_UZ, :)
                     Qjac_n (var_UZ, :)   = Qjac_n (var_UZ, :) - TG_NUM(var_UZ) * QvmsF_n (var_UZ, :)
                     Qjac_k (var_UZ, :)   = Qjac_k (var_UZ, :) - TG_NUM(var_UZ) * QvmsF_k (var_UZ, :)
                     Qjac_kn(var_UZ, :)   = Qjac_kn(var_UZ, :) - TG_NUM(var_UZ) * QvmsF_kn(var_UZ, :)

                     Qjac_p (var_Up, :) = Qjac_p (var_Up, :) - TG_NUM(var_Up) *  (  BR0 *( QvmsF_p (var_UR, :) ) &
                                                                                  + BZ0 *( QvmsF_p (var_UZ, :) ) &
                                                                                  + Bp0 *( QvmsF_p (var_Up, :) ) )
                     Qjac_n (var_Up, :) = Qjac_n (var_Up, :) - TG_NUM(var_Up) *  (  BR0 *( QvmsF_n (var_UR, :) ) &
                                                                                  + BZ0 *( QvmsF_n (var_UZ, :) ) &
                                                                                  + Bp0 *( QvmsF_n (var_Up, :) ) )
                     Qjac_k (var_Up, :) = Qjac_k (var_Up, :) - TG_NUM(var_Up) *  (  BR0 *( QvmsF_k (var_UR, :) ) &
                                                                                  + BZ0 *( QvmsF_k (var_UZ, :) ) &
                                                                                  + Bp0 *( QvmsF_k (var_Up, :) ) )
                     Qjac_kn(var_Up, :) = Qjac_kn(var_Up, :) - TG_NUM(var_Up) *  (  BR0 *( QvmsF_kn(var_UR, :) ) &
                                                                                  + BZ0 *( QvmsF_kn(var_UZ, :) ) &
                                                                                  + Bp0 *( QvmsF_kn(var_Up, :) ) )
                  CASE(-1)
                     ! --- no VMS
                  END SELECT


                  !###################################################################################################
                  !#  Div(V) STABILISATION                                                                           #
                  !###################################################################################################
                  Qjac_p (var_UR,var_UR) = Qjac_p (var_UR,var_UR) -  visco_divV * divU_UR     * ( v_R + v / R )
                  Qjac_p (var_UR,var_UZ) = Qjac_p (var_UR,var_UZ) -  visco_divV * divU_UZ     * ( v_R + v / R )
                  Qjac_n (var_UR,var_Up) = Qjac_n (var_UR,var_Up) -  visco_divV * divU_Up__n  * ( v_R + v / R )
                  Qjac_p (var_UR,var_T)  = Qjac_p (var_UR,var_T)  - dvisco_divV_dT * T * divU * ( v_R + v / R ) 

                  Qjac_p (var_UZ,var_UR) = Qjac_p (var_UZ,var_UR) -  visco_divV * divU_UR     * v_Z 
                  Qjac_p (var_UZ,var_UZ) = Qjac_p (var_UZ,var_UZ) -  visco_divV * divU_UZ     * v_Z 
                  Qjac_n (var_UZ,var_Up) = Qjac_n (var_UZ,var_Up) -  visco_divV * divU_Up__n  * v_Z 
                  Qjac_p (var_UZ,var_T ) = Qjac_p (var_UZ,var_T ) - dvisco_divV_dT * T * divU * v_Z

                  Qjac_p (var_Up,var_AR) = Qjac_p (var_Up,var_AR) - visco_divV * divU    * BgradVstar_AR__p
                  Qjac_n (var_Up,var_AR) = Qjac_n (var_Up,var_AR) - visco_divV * divU    * BgradVstar_AR__n
                  Qjac_k (var_Up,var_AR) = Qjac_k (var_Up,var_AR) - visco_divV * divU    * BgradVstar_AR__k

                  Qjac_p (var_Up,var_AZ) = Qjac_p (var_Up,var_AZ) - visco_divV * divU    * BgradVstar_AZ__p
                  Qjac_n (var_Up,var_AZ) = Qjac_n (var_Up,var_AZ) - visco_divV * divU    * BgradVstar_AZ__n
                  Qjac_k (var_Up,var_AZ) = Qjac_k (var_Up,var_AZ) - visco_divV * divU    * BgradVstar_AZ__k

                  Qjac_p (var_Up,var_A3) = Qjac_p (var_Up,var_A3) - visco_divV * divU    * BgradVstar_A3__p
                  Qjac_k (var_Up,var_A3) = Qjac_k (var_Up,var_A3) - visco_divV * divU    * BgradVstar_A3__k

                  Qjac_p (var_Up,var_UR) = Qjac_p (var_Up,var_UR) -  visco_divV * divU_UR     * BgradVstar__p
                  Qjac_k (var_Up,var_UR) = Qjac_k (var_Up,var_UR) -  visco_divV * divU_UR     * BgradVstar__k
                  Qjac_p (var_Up,var_UZ) = Qjac_p (var_Up,var_UZ) -  visco_divV * divU_UZ     * BgradVstar__p
                  Qjac_k (var_Up,var_UZ) = Qjac_k (var_Up,var_UZ) -  visco_divV * divU_UZ     * BgradVstar__k
                  Qjac_n (var_Up,var_Up) = Qjac_n (var_Up,var_Up) -  visco_divV * divU_Up__n  * BgradVstar__p
                  Qjac_kn(var_Up,var_Up) = Qjac_kn(var_Up,var_Up) -  visco_divV * divU_Up__n  * BgradVstar__k
                  Qjac_p (var_Up,var_T ) = Qjac_p (var_Up,var_T ) - dvisco_divV_dT * T * divU * BgradVstar__p
                  Qjac_k (var_Up,var_T ) = Qjac_k (var_Up,var_T ) - dvisco_divV_dT * T * divU * BgradVstar__k

                  if (use_fft) then
                    index_kl =       n_var*n_degrees*(k-1) +       n_var*(l-1) + 1
                    do ivar= 1,n_var
                      do kvar= 1,n_var
                        ij = ivar
                        kl = index_kl + (kvar-1)

                        amat(ivar,kvar)  = (1.d0+zeta)*Pjac(ivar,kvar) - tstep * theta * Qjac_p(ivar,kvar)
                        ELM_p(mp,kl,ij)  = ELM_p(mp,kl,ij)  +  wst * amat(ivar,kvar) * R * xjac

                        amat(ivar,kvar)  = - tstep * theta * Qjac_k(ivar,kvar)
                        ELM_k(mp,kl,ij)  = ELM_k(mp,kl,ij)  +  wst * amat(ivar,kvar) * R * xjac

                        amat(ivar,kvar)  = - tstep * theta * Qjac_n(ivar,kvar)
                        ELM_n(mp,kl,ij)  = ELM_n(mp,kl,ij)  +  wst * amat(ivar,kvar) * R * xjac

                        amat(ivar,kvar)  = - tstep * theta * Qjac_kn(ivar,kvar)
                        ELM_kn(mp,kl,ij) = ELM_kn(mp,kl,ij) +  wst * amat(ivar,kvar) * R * xjac

                        amat(ivar,kvar)  = - tstep * theta * Qjac_pnn(ivar,kvar)
                        ELM_pnn(mp,kl,ij) = ELM_pnn(mp,kl,ij) +  wst * amat(ivar,kvar) * R * xjac
                      enddo
                    enddo
                  else
                    index_kl = n_tor_local*n_var*n_degrees*(k-1) + n_tor_local*n_var*(l-1) + in - n_tor_start +1
                    do ivar= 1,n_var
                      do kvar= 1,n_var
                        ij = index_ij + (ivar-1) * n_tor_local
                        kl = index_kl + (kvar-1) * n_tor_local

                        amat(ivar,kvar) = (1.d0+zeta)*Pjac(ivar,kvar) - tstep * theta * Qjac_p(ivar,kvar)
                        ELM(ij,kl)      = ELM(ij,kl) + wst * amat(ivar,kvar) * R * xjac

                        amat(ivar,kvar) = - tstep * theta * Qjac_k(ivar,kvar)
                        ELM(ij,kl)      = ELM(ij,kl) + wst * amat(ivar,kvar) * R * xjac

                        amat(ivar,kvar) = - tstep * theta * Qjac_n(ivar,kvar)
                        ELM(ij,kl)      = ELM(ij,kl) + wst * amat(ivar,kvar) * R * xjac

                        amat(ivar,kvar) = - tstep * theta * Qjac_kn(ivar,kvar)
                        ELM(ij,kl)      = ELM(ij,kl) + wst * amat(ivar,kvar) * R * xjac

                        amat(ivar,kvar) = - tstep * theta * Qjac_pnn(ivar,kvar)
                        ELM(ij,kl)      = ELM(ij,kl) + wst * amat(ivar,kvar) * R * xjac
                      enddo
                    enddo
                  endif

                enddo ! in loop (n_tor, or not...)

              enddo ! l loop n_degrees
            enddo ! k loop (n_vertex)

          enddo ! im loop (n_tor, or not...)

        enddo ! mp loop (n_plane)

      enddo ! ms loop
    enddo ! mt loop


    if (use_fft) then
      do i_v = 1, n_var
        do j_loc=1, n_vertex_max*n_var*n_degrees
          i_loc = n_var*n_degrees*(i-1) + n_var * (j-1) + i_v 
          in_fft =  ELM_p(1:n_plane,j_loc,i_v)
#ifdef USE_FFTW
          call dfftw_execute_dft_r2c(fftw_plan, in_fft, out_fft)
#else
          call my_fft(in_fft, out_fft, n_plane)
#endif
          do m=1,(n_tor+1)/2

            index_m = n_tor*(j_loc-1) + max(2*(m-1),1)

            do k=1,(n_tor+1)/2

              index_k = n_tor*(i_loc-1) + max(2*(k-1),1)

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

          if (maxval(abs(ELM_n(1:n_plane,j_loc, i_v))) .ne. 0.d0) then

            in_fft =  ELM_n(1:n_plane,j_loc, i_v)
#ifdef USE_FFTW
            call dfftw_execute_dft_r2c(fftw_plan, in_fft, out_fft)
#else
            call my_fft(in_fft, out_fft, n_plane)
#endif

            do m=1,(n_tor+1)/2
              im = max(2*(m-1),1)
              index_m = n_tor*(j_loc-1) + max(2*(m-1),1)

              do k=1,(n_tor+1)/2

                index_k = n_tor*(i_loc-1) + max(2*(k-1),1)

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

          if (maxval(abs(ELM_k(1:n_plane,j_loc, i_v))) .ne. 0.d0) then

            in_fft =  ELM_k(1:n_plane,j_loc, i_v)

#ifdef USE_FFTW
            call dfftw_execute_dft_r2c(fftw_plan, in_fft, out_fft)
#else
            call my_fft(in_fft, out_fft, n_plane)
#endif

            do m=1,(n_tor+1)/2

              index_m = n_tor*(j_loc-1) + max(2*(m-1),1)

              do k=1,(n_tor+1)/2

                ik      = max(2*(k-1),1)
                index_k = n_tor*(i_loc-1) + max(2*(k-1),1)

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

          if (maxval(abs(ELM_kn(1:n_plane,j_loc, i_v))) .ne. 0.d0) then

            in_fft =  ELM_kn(1:n_plane,j_loc, i_v)

#ifdef USE_FFTW
            call dfftw_execute_dft_r2c(fftw_plan, in_fft, out_fft)
#else
            call my_fft(in_fft, out_fft, n_plane)
#endif

            do m=1,(n_tor+1)/2

              im      = max(2*(m-1),1)
              index_m = n_tor*(j_loc-1) + max(2*(m-1),1)

              do k=1,(n_tor+1)/2

                ik      = max(2*(k-1),1)
                index_k = n_tor*(i_loc-1) + max(2*(k-1),1)

                l = (k-1) + (m-1)

                if ( (l .ge. 0) .and. (l .le. n_plane/2) ) then
                  ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   - real(out_fft(l+1)) * float(mode(im)) * float(mode(ik))
                  ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   + imag(out_fft(l+1)) * float(mode(im)) * float(mode(ik))
                  ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) + imag(out_fft(l+1)) * float(mode(im)) * float(mode(ik))
                  ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) + real(out_fft(l+1)) * float(mode(im)) * float(mode(ik))
                elseif ( (l .lt. 0) .and. (abs(l) .le. n_plane/2) ) then
                  ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   - real(out_fft(abs(l)+1)) * float(mode(im)) * float(mode(ik))
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

          if (maxval(abs(ELM_pnn(1:n_plane,j_loc, i_v))) .ne. 0.d0) then

            in_fft =  ELM_pnn(1:n_plane,j_loc, i_v)

#ifdef USE_FFTW
            call dfftw_execute_dft_r2c(fftw_plan, in_fft, out_fft)
#else
            call my_fft(in_fft, out_fft, n_plane)
#endif

            do m=1,(n_tor+1)/2

              im      = max(2*(m-1),1)
              index_m = n_tor*(j_loc-1) + max(2*(m-1),1)

              do k=1,(n_tor+1)/2

                ik      = max(2*(k-1),1)
                index_k = n_tor*(i_loc-1) + max(2*(k-1),1)

                l = (k-1) + (m-1)

                if ( (l .ge. 0) .and. (l .le. n_plane/2) ) then
                  ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   - real(out_fft(l+1)) * float(mode(im)**2)
                  ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   + imag(out_fft(l+1)) * float(mode(im)**2)
                  ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) + imag(out_fft(l+1)) * float(mode(im)**2)
                  ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) + real(out_fft(l+1)) * float(mode(im)**2)
                elseif ( (l .lt. 0) .and. (abs(l) .le. n_plane/2) ) then
                  ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   + real(out_fft(abs(l)+1)) * float(mode(im)**2)
                  ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   - imag(out_fft(abs(l)+1)) * float(mode(im)**2)
                  ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) - imag(out_fft(abs(l)+1)) * float(mode(im)**2)
                  ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) - real(out_fft(abs(l)+1)) * float(mode(im)**2)
                endif

                l = (k-1) - (m-1)

                if ( (l .ge. 0) .and. (l .le. n_plane/2) ) then
                  ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   - real(out_fft(l+1)) * float(mode(im)**2)
                  ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   + imag(out_fft(l+1)) * float(mode(im)**2)
                  ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) - imag(out_fft(l+1)) * float(mode(im)**2)
                  ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) - real(out_fft(l+1)) * float(mode(im)**2)
                elseif ( (l .lt. 0) .and. (abs(l) .le. n_plane/2) ) then
                  ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   - real(out_fft(abs(l)+1)) * float(mode(im)**2)
                  ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   - imag(out_fft(abs(l)+1)) * float(mode(im)**2)
                  ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) + imag(out_fft(abs(l)+1)) * float(mode(im)**2)
                  ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) - real(out_fft(abs(l)+1)) * float(mode(im)**2)
                endif

              enddo

            enddo

          endif
 
        enddo

      enddo

    endif ! apply fft (or not)

  enddo ! j loop n_degrees
enddo ! i loop (n_vertex)

if (n_tor .le. n_tor_fft_thresh) return

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


CONTAINS


subroutine my_fft(in_fft,out_fft,n)
  implicit none
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


end subroutine element_matrix_fft

end module mod_elt_matrix_fft
