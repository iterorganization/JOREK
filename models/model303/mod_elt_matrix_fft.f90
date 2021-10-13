module mod_elt_matrix_fft

  implicit none

contains

#include "corr_neg_include.f90"

subroutine element_matrix_fft(element, nodes, xpoint2, xcase2, R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint, Z_xpoint, &
                              ELM, RHS, tid, ELM_p, ELM_n, ELM_k, ELM_kn, RHS_p, RHS_k,                               &
                              eq_g, eq_s, eq_t, eq_p, eq_ss, eq_st, eq_tt, delta_g, delta_s, delta_t,                 &
                              i_tor_min, i_tor_max, aux_nodes)
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
use diffusivities, only: get_dperp, get_zkperp
use equil_info, only : get_psi_n
use mod_bootstrap_functions
use mod_sources

implicit none

type (type_element)        :: element                 ! finite element
type (type_node)           :: nodes(n_vertex_max)     ! fluid variables
type (type_node), optional :: aux_nodes(n_vertex_max) ! moments of particles


#define DIM0 n_tor*n_vertex_max*(n_order+1)*n_var

real*8, dimension (DIM0,DIM0)  :: ELM
real*8, dimension (DIM0)       :: RHS
integer, intent(in)            :: tid
integer, intent(in)            :: i_tor_min, i_tor_max

integer    :: i, j, ms, mt, mp, k, l, index_ij, index_kl, index, index_k, index_m, m, ik, xcase2
integer    :: n_tor_start, n_tor_end, n_tor_local
integer    :: in, im, ij1, ij2, ij3, ij4, ij5, ij6, ij7, kl1, kl2, kl3, kl4, kl5, kl6, kl7
real*8     :: wst, xjac, xjac_s, xjac_t, xjac_x, xjac_y, BigR, r2, phi, delta_phi, eps_cyl
real*8     :: current_source(n_gauss,n_gauss),particle_source(n_gauss,n_gauss),heat_source(n_gauss,n_gauss)
real*8     :: R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint(2), Z_xpoint(2), dj_dpsi, dj_dz, source_pellet, source_volume
real*8     :: Bgrad_rho_star,     Bgrad_rho,     Bgrad_T_star,  Bgrad_T, BB2
real*8     :: Bgrad_rho_star_psi, Bgrad_rho_psi, Bgrad_rho_rho, Bgrad_T_star_psi, Bgrad_T_psi, Bgrad_T_T, BB2_psi
real*8     :: Bgrad_rho_rho_n, Bgrad_T_T_n, Bgrad_rho_k_star, Bgrad_T_k_star, ZKpar_T, dZKpar_dT
real*8     :: D_prof, ZK_prof, psi_norm, theta, zeta, delta_u_x, delta_u_y, delta_ps_x, delta_ps_y
real*8     :: rhs_ij_1,   rhs_ij_2,   rhs_ij_3,   rhs_ij_4,   rhs_ij_5,   rhs_ij_6,   rhs_ij_7
real*8     :: rhs_ij_1_k, rhs_ij_2_k, rhs_ij_3_k, rhs_ij_4_k, rhs_ij_5_k, rhs_ij_6_k, rhs_ij_7_k
real*8     :: rhs_stab_1, rhs_stab_2, rhs_stab_3, rhs_stab_4, rhs_stab_5, rhs_stab_6, rhs_stab_7

real*8     :: v, v_x, v_y, v_s, v_t, v_p, v_ss, v_st, v_tt, v_xx, v_xy, v_yy
real*8     :: ps0, ps0_x, ps0_y, ps0_p, ps0_s, ps0_t, ps0_ss, ps0_tt, ps0_st, ps0_xx, ps0_yy, ps0_xy
real*8     :: zj0, zj0_x, zj0_y, zj0_p, zj0_s, zj0_t
real*8     :: u0, u0_x, u0_y, u0_p, u0_s, u0_t, u0_ss, u0_tt, u0_st, u0_xx, u0_xy, u0_yy
real*8     :: w0, w0_x, w0_y, w0_p, w0_s, w0_t, w0_ss, w0_st, w0_tt, w0_xx, w0_xy, w0_yy
real*8     :: r0, r0_x, r0_y, r0_p, r0_s, r0_t, r0_ss, r0_st, r0_tt, r0_xx, r0_xy, r0_yy, r0_hat, r0_x_hat, r0_y_hat, r0_corr
real*8     :: T0, T0_x, T0_y, T0_p, T0_s, T0_t, T0_ss, T0_st, T0_tt, T0_xx, T0_xy, T0_yy, T0_corr
real*8     :: psi, psi_x, psi_y, psi_p, psi_s, psi_t, psi_ss, psi_st, psi_tt, psi_xx, psi_xy, psi_yy
real*8     :: zj, zj_x, zj_y, zj_p, zj_s, zj_t, zj_ss, zj_st, zj_tt
real*8     :: u, u_x, u_y, u_p, u_s, u_t, u_ss, u_st, u_tt, u_xx, u_xy, u_yy
real*8     :: w, w_x, w_y, w_p, w_s, w_t, w_ss, w_st, w_tt, w_xx, w_xy, w_yy
real*8     :: rho, rho_x, rho_y, rho_s, rho_t, rho_p, rho_hat, rho_x_hat, rho_y_hat, rho_ss, rho_st, rho_tt, rho_xx, rho_xy, rho_yy
real*8     :: T, T_x, T_y, T_s, T_t, T_p, T_ss, T_st, T_tt, T_xx, T_xy, T_yy
real*8     :: Ti0, Ti0_x, Ti0_y, Te0, Te0_x, Te0_y
real*8     :: zTi, zTi_x, zTi_y, zTe, zTe_x, zTe_y, zn_x, zn_y
real*8	   :: Jb_0 , Jb
real*8     :: Vpar, Vpar_x, Vpar_y, Vpar_p, Vpar_s, Vpar_t, Vpar_ss, Vpar_st, Vpar_tt, Vpar_xx, Vpar_yy, Vpar_xy
real*8     :: P0, P0_s, P0_t, P0_x, P0_y, P0_p, P0_ss, P0_st, P0_tt, P0_xx, P0_xy, P0_yy
real*8     :: P0_x_rho, P0_xx_rho, P0_y_rho, P0_yy_rho, P0_xy_rho
real*8     :: P0_x_T,   P0_xx_T,   P0_y_T,   P0_yy_T,   P0_xy_T
real*8     :: Vpar0, Vpar0_s, Vpar0_t, Vpar0_p, Vpar0_x, Vpar0_y, Vpar0_ss, Vpar0_st, Vpar0_tt, Vpar0_xx, Vpar0_yy,Vpar0_xy
real*8     :: BigR_x, vv2, eta_T, visco_T, deta_dT, d2eta_d2T, dvisco_dT, d2visco_dT2, visco_num_T, dvisco_num_dT, eta_num_T, deta_num_dT, W_dia, W_dia_rho, W_dia_T
real*8     :: eta_T_ohm, deta_dT_ohm
real*8     :: ZK_par_num, T0_ps0_x, T_ps0_x, T0_psi_x, T0_ps0_y, T_ps0_y, T0_psi_y, v_ps0_x, v_psi_x, v_ps0_y, v_psi_y
real*8     :: amat_11, amat_12, amat_13, amat_14, amat_15, amat_16, amat_17
real*8     :: amat_21, amat_22, amat_23, amat_24, amat_25, amat_26, amat_27
real*8     :: amat_31, amat_32, amat_33, amat_34, amat_35, amat_36, amat_37
real*8     :: amat_41, amat_42, amat_43, amat_44, amat_45, amat_46, amat_47
real*8     :: amat_51, amat_52, amat_53, amat_54, amat_55, amat_56, amat_57
real*8     :: amat_61, amat_62, amat_63, amat_64, amat_65, amat_66, amat_67
real*8     :: amat_71, amat_72, amat_73, amat_74, amat_75, amat_76, amat_77
real*8     :: amat_11_k, amat_12_k, amat_13_k, amat_14_k, amat_15_k, amat_16_k, amat_17_k
real*8     :: amat_21_k, amat_22_k, amat_23_k, amat_24_k, amat_25_k, amat_26_k, amat_27_k
real*8     :: amat_31_k, amat_32_k, amat_33_k, amat_34_k, amat_35_k, amat_36_k, amat_37_k
real*8     :: amat_41_k, amat_42_k, amat_43_k, amat_44_k, amat_45_k, amat_46_k, amat_47_k
real*8     :: amat_51_k, amat_52_k, amat_53_k, amat_54_k, amat_55_k, amat_56_k, amat_57_k
real*8     :: amat_61_k, amat_62_k, amat_63_k, amat_64_k, amat_65_k, amat_66_k, amat_67_k
real*8     :: amat_71_k, amat_72_k, amat_73_k, amat_74_k, amat_75_k, amat_76_k, amat_77_k
real*8     :: amat_11_n, amat_12_n, amat_13_n, amat_14_n, amat_15_n, amat_16_n, amat_17_n
real*8     :: amat_21_n, amat_22_n, amat_23_n, amat_24_n, amat_25_n, amat_26_n, amat_27_n
real*8     :: amat_31_n, amat_32_n, amat_33_n, amat_34_n, amat_35_n, amat_36_n, amat_37_n
real*8     :: amat_41_n, amat_42_n, amat_43_n, amat_44_n, amat_45_n, amat_46_n, amat_47_n
real*8     :: amat_51_n, amat_52_n, amat_53_n, amat_54_n, amat_55_n, amat_56_n, amat_57_n
real*8     :: amat_61_n, amat_62_n, amat_63_n, amat_64_n, amat_65_n, amat_66_n, amat_67_n
real*8     :: amat_71_n, amat_72_n, amat_73_n, amat_74_n, amat_75_n, amat_76_n, amat_77_n
real*8     :: amat_11_kn, amat_12_kn, amat_13_kn, amat_14_kn, amat_15_kn, amat_16_kn, amat_17_kn
real*8     :: amat_21_kn, amat_22_kn, amat_23_kn, amat_24_kn, amat_25_kn, amat_26_kn, amat_27_kn
real*8     :: amat_31_kn, amat_32_kn, amat_33_kn, amat_34_kn, amat_35_kn, amat_36_kn, amat_37_kn
real*8     :: amat_41_kn, amat_42_kn, amat_43_kn, amat_44_kn, amat_45_kn, amat_46_kn, amat_47_kn
real*8     :: amat_51_kn, amat_52_kn, amat_53_kn, amat_54_kn, amat_55_kn, amat_56_kn, amat_57_kn
real*8     :: amat_61_kn, amat_62_kn, amat_63_kn, amat_64_kn, amat_65_kn, amat_66_kn, amat_67_kn
real*8     :: amat_71_kn, amat_72_kn, amat_73_kn, amat_74_kn, amat_75_kn, amat_76_kn, amat_77_kn
real*8     :: TG_num1, TG_num2, TG_num5, TG_num6, TG_num7
!==================MB: velocity profile is kept by a source which compensating diffusion
real*8     :: Vt0,Omega_tor0_x,Omega_tor0_y,Vt0_x,Vt0_y
real*8     :: V_source(n_gauss,n_gauss), Vt_x_psi, Vt_y_psi, Omega_tor_x_psi, Omega_tor_y_psi
real*8     :: dV_dpsi_source(n_gauss,n_gauss),dV_dz_source(n_gauss,n_gauss)
real*8     :: dV_dpsi2,dV_dz2,dV_dpsi_dz,dV_dpsi3,dV_dpsi_dz2,dV_dpsi2_dz
!=======================================
real*8     :: eq_zne(n_gauss,n_gauss), eq_zTe(n_gauss,n_gauss)
real*8     :: dn_dpsi(n_gauss,n_gauss),dn_dz,dn_dpsi2,dn_dz2,dn_dpsi_dz,dn_dpsi3,dn_dpsi_dz2,dn_dpsi2_dz
real*8     :: dT_dpsi(n_gauss,n_gauss),dT_dz,dT_dpsi2,dT_dz2,dT_dpsi_dz,dT_dpsi3,dT_dpsi_dz2,dT_dpsi2_dz
logical    :: xpoint2, use_fft
real*8     :: w00_xx, w00_yy                                                                                                                                                                
!======================================= NEO
real*8     :: Btheta2
real*8     :: epsil, Btheta2_psi
real*8, dimension(n_gauss,n_gauss)    :: amu_neo_prof, aki_neo_prof
!======================================= NEO

real*8     :: in_fft(1:n_plane)
complex*16 :: out_fft(1:n_plane)
integer*8  :: plan

integer    :: i_v, i_loc, j_loc
#define DIM1 n_plane
#define DIM2 1:n_vertex_max*n_var*(n_order+1)

real*8, dimension(DIM1, DIM2, DIM2) :: ELM_p
real*8, dimension(DIM1, DIM2, DIM2) :: ELM_n
real*8, dimension(DIM1, DIM2, DIM2) :: ELM_k
real*8, dimension(DIM1, DIM2, DIM2) :: ELM_kn
real*8, dimension(DIM1, DIM2)       :: RHS_p
real*8, dimension(DIM1, DIM2)       :: RHS_k

real*8, dimension(n_gauss,n_gauss)    :: x_g, x_s, x_t
real*8, dimension(n_gauss,n_gauss)    :: x_ss, x_st, x_tt
real*8, dimension(n_gauss,n_gauss)    :: y_g, y_s, y_t
real*8, dimension(n_gauss,n_gauss)    :: y_ss, y_st, y_tt

real*8, dimension(n_plane,n_var,n_gauss,n_gauss) :: eq_g, eq_s, eq_t
real*8, dimension(n_plane,n_var,n_gauss,n_gauss) :: eq_p
real*8, dimension(n_plane,n_var,n_gauss,n_gauss) :: eq_ss, eq_st, eq_tt
real*8, dimension(n_plane,n_var,n_gauss,n_gauss) :: delta_g, delta_s, delta_t

real*8, dimension(n_tor,n_plane) :: HHZ, HHZ_p, HHZ_pp

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

! --- Take time evolution parameters from phys_module
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

n_tor_local = n_tor_end - n_tor_start + 1
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
  
rhs_ij_1   = 0.d0 ; rhs_ij_2   = 0.d0 ; rhs_ij_3   = 0.d0 ; rhs_ij_4   = 0.d0 ; rhs_ij_5   = 0.d0 ; rhs_ij_6   = 0.d0 ; rhs_ij_7   = 0.d0
rhs_ij_1_k = 0.d0 ; rhs_ij_2_k = 0.d0 ; rhs_ij_3_k = 0.d0 ; rhs_ij_4_k = 0.d0 ; rhs_ij_5_k = 0.d0 ; rhs_ij_6_k = 0.d0 ; rhs_ij_7_k = 0.d0
rhs_stab_1 = 0.d0 ; rhs_stab_2 = 0.d0 ; rhs_stab_3 = 0.d0 ; rhs_stab_4 = 0.d0 ; rhs_stab_5 = 0.d0 ; rhs_stab_6 = 0.d0 ; rhs_stab_7 = 0.d0

amat_11 = 0.d0 ; amat_12 = 0.d0 ; amat_13 = 0.d0 ; amat_14 = 0.d0 ; amat_15 = 0.d0 ; amat_16 = 0.d0 ; amat_17 = 0.d0
amat_21 = 0.d0 ; amat_22 = 0.d0 ; amat_23 = 0.d0 ; amat_24 = 0.d0 ; amat_25 = 0.d0 ; amat_26 = 0.d0 ; amat_27 = 0.d0
amat_31 = 0.d0 ; amat_32 = 0.d0 ; amat_33 = 0.d0 ; amat_34 = 0.d0 ; amat_35 = 0.d0 ; amat_36 = 0.d0 ; amat_37 = 0.d0
amat_41 = 0.d0 ; amat_42 = 0.d0 ; amat_43 = 0.d0 ; amat_44 = 0.d0 ; amat_45 = 0.d0 ; amat_46 = 0.d0 ; amat_47 = 0.d0
amat_51 = 0.d0 ; amat_52 = 0.d0 ; amat_53 = 0.d0 ; amat_54 = 0.d0 ; amat_55 = 0.d0 ; amat_56 = 0.d0 ; amat_57 = 0.d0
amat_61 = 0.d0 ; amat_62 = 0.d0 ; amat_63 = 0.d0 ; amat_64 = 0.d0 ; amat_65 = 0.d0 ; amat_66 = 0.d0 ; amat_67 = 0.d0
amat_71 = 0.d0 ; amat_72 = 0.d0 ; amat_73 = 0.d0 ; amat_74 = 0.d0 ; amat_75 = 0.d0 ; amat_76 = 0.d0 ; amat_77 = 0.d0
amat_11_k = 0.d0 ; amat_12_k = 0.d0 ; amat_13_k = 0.d0 ; amat_14_k = 0.d0 ; amat_15_k = 0.d0 ; amat_16_k = 0.d0 ; amat_17_k = 0.d0
amat_21_k = 0.d0 ; amat_22_k = 0.d0 ; amat_23_k = 0.d0 ; amat_24_k = 0.d0 ; amat_25_k = 0.d0 ; amat_26_k = 0.d0 ; amat_27_k = 0.d0
amat_31_k = 0.d0 ; amat_32_k = 0.d0 ; amat_33_k = 0.d0 ; amat_34_k = 0.d0 ; amat_35_k = 0.d0 ; amat_36_k = 0.d0 ; amat_37_k = 0.d0
amat_41_k = 0.d0 ; amat_42_k = 0.d0 ; amat_43_k = 0.d0 ; amat_44_k = 0.d0 ; amat_45_k = 0.d0 ; amat_46_k = 0.d0 ; amat_47_k = 0.d0
amat_51_k = 0.d0 ; amat_52_k = 0.d0 ; amat_53_k = 0.d0 ; amat_54_k = 0.d0 ; amat_55_k = 0.d0 ; amat_56_k = 0.d0 ; amat_57_k = 0.d0
amat_61_k = 0.d0 ; amat_62_k = 0.d0 ; amat_63_k = 0.d0 ; amat_64_k = 0.d0 ; amat_65_k = 0.d0 ; amat_66_k = 0.d0 ; amat_67_k = 0.d0
amat_71_k = 0.d0 ; amat_72_k = 0.d0 ; amat_73_k = 0.d0 ; amat_74_k = 0.d0 ; amat_75_k = 0.d0 ; amat_76_k = 0.d0 ; amat_77_k = 0.d0
amat_11_n = 0.d0 ; amat_12_n = 0.d0 ; amat_13_n = 0.d0 ; amat_14_n = 0.d0 ; amat_15_n = 0.d0 ; amat_16_n = 0.d0 ; amat_17_n = 0.d0
amat_21_n = 0.d0 ; amat_22_n = 0.d0 ; amat_23_n = 0.d0 ; amat_24_n = 0.d0 ; amat_25_n = 0.d0 ; amat_26_n = 0.d0 ; amat_27_n = 0.d0
amat_31_n = 0.d0 ; amat_32_n = 0.d0 ; amat_33_n = 0.d0 ; amat_34_n = 0.d0 ; amat_35_n = 0.d0 ; amat_36_n = 0.d0 ; amat_37_n = 0.d0
amat_41_n = 0.d0 ; amat_42_n = 0.d0 ; amat_43_n = 0.d0 ; amat_44_n = 0.d0 ; amat_45_n = 0.d0 ; amat_46_n = 0.d0 ; amat_47_n = 0.d0
amat_51_n = 0.d0 ; amat_52_n = 0.d0 ; amat_53_n = 0.d0 ; amat_54_n = 0.d0 ; amat_55_n = 0.d0 ; amat_56_n = 0.d0 ; amat_57_n = 0.d0
amat_61_n = 0.d0 ; amat_62_n = 0.d0 ; amat_63_n = 0.d0 ; amat_64_n = 0.d0 ; amat_65_n = 0.d0 ; amat_66_n = 0.d0 ; amat_67_n = 0.d0
amat_71_n = 0.d0 ; amat_72_n = 0.d0 ; amat_73_n = 0.d0 ; amat_74_n = 0.d0 ; amat_75_n = 0.d0 ; amat_76_n = 0.d0 ; amat_77_n = 0.d0
amat_11_kn = 0.d0 ; amat_12_kn = 0.d0 ; amat_13_kn = 0.d0 ; amat_14_kn = 0.d0 ; amat_15_kn = 0.d0 ; amat_16_kn = 0.d0 ; amat_17_kn = 0.d0
amat_21_kn = 0.d0 ; amat_22_kn = 0.d0 ; amat_23_kn = 0.d0 ; amat_24_kn = 0.d0 ; amat_25_kn = 0.d0 ; amat_26_kn = 0.d0 ; amat_27_kn = 0.d0
amat_31_kn = 0.d0 ; amat_32_kn = 0.d0 ; amat_33_kn = 0.d0 ; amat_34_kn = 0.d0 ; amat_35_kn = 0.d0 ; amat_36_kn = 0.d0 ; amat_37_kn = 0.d0
amat_41_kn = 0.d0 ; amat_42_kn = 0.d0 ; amat_43_kn = 0.d0 ; amat_44_kn = 0.d0 ; amat_45_kn = 0.d0 ; amat_46_kn = 0.d0 ; amat_47_kn = 0.d0
amat_51_kn = 0.d0 ; amat_52_kn = 0.d0 ; amat_53_kn = 0.d0 ; amat_54_kn = 0.d0 ; amat_55_kn = 0.d0 ; amat_56_kn = 0.d0 ; amat_57_kn = 0.d0
amat_61_kn = 0.d0 ; amat_62_kn = 0.d0 ; amat_63_kn = 0.d0 ; amat_64_kn = 0.d0 ; amat_65_kn = 0.d0 ; amat_66_kn = 0.d0 ; amat_67_kn = 0.d0
amat_71_kn = 0.d0 ; amat_72_kn = 0.d0 ; amat_73_kn = 0.d0 ; amat_74_kn = 0.d0 ; amat_75_kn = 0.d0 ; amat_76_kn = 0.d0 ; amat_77_kn = 0.d0

!---------------------------------------------------- value of (x,y) and derivatives on Gaussian points
x_g  = 0.d0; x_s  = 0.d0; x_t  = 0.d0; x_st  = 0.d0; x_ss  = 0.d0; x_tt  = 0.d0;
y_g  = 0.d0; y_s  = 0.d0; y_t  = 0.d0; y_st  = 0.d0; y_ss  = 0.d0; y_tt  = 0.d0;
eq_g = 0.d0; eq_s = 0.d0; eq_t = 0.d0; eq_st = 0.d0; eq_ss = 0.d0; eq_tt = 0.d0; eq_p = 0.d0;

delta_g = 0.d0; delta_s = 0.d0; delta_t = 0.d0

current_source  = 0.d0
particle_source = 0.d0
heat_source     = 0.d0
V_source=0.d0
dV_dpsi_source=0.d0
dV_dz_source=0.d0
eq_zne          = 0.d0
eq_zTe          = 0.d0         
!======================================= NEO
if ( NEO ) then 
   amu_neo_prof   = 0.d0
   aki_neo_prof   = 0.d0
endif
!======================================= NEO

do i=1,n_vertex_max
  do j=1,n_order+1
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

      end do
    end do

    do ms=1, n_gauss
      do mt=1, n_gauss
        do k=1,n_var
          do in=1,n_tor
            do mp=1,n_plane
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
    call sources(xpoint2, xcase2, y_g(ms,mt), Z_xpoint, eq_g(1,1,ms,mt),psi_axis,psi_bnd,particle_source(ms,mt),heat_source(ms,mt))

    ! Source of parallel velocity
    if ( ( abs(V_0) .ge. 1.e-12 ) .or. ( num_rot ) ) then
      call velocity(xpoint2, xcase2, y_g(ms,mt), Z_xpoint, eq_g(1,1,ms,mt), psi_axis, psi_bnd, V_source(ms,mt), &
           dV_dpsi_source(ms,mt),dV_dz_source(ms,mt),dV_dpsi2,dV_dz2,dV_dpsi_dz,dV_dpsi3,dV_dpsi_dz2, dV_dpsi2_dz)
    endif

    call density(xpoint2, xcase2, y_g(ms,mt), Z_xpoint, eq_g(1,1,ms,mt),psi_axis,psi_bnd,eq_zne(ms,mt), &
                 dn_dpsi(ms,mt),dn_dz,dn_dpsi2,dn_dz2,dn_dpsi_dz,dn_dpsi3,dn_dpsi_dz2, dn_dpsi2_dz)

    call temperature(xpoint2, xcase2, y_g(ms,mt), Z_xpoint, eq_g(1,1,ms,mt),psi_axis,psi_bnd,eq_zTe(ms,mt), &
                     dT_dpsi(ms,mt),dT_dz,dT_dpsi2,dT_dz2,dT_dpsi_dz,dT_dpsi3,dT_dpsi_dz2, dT_dpsi2_dz)

    !======================================= NEO
    if ( NEO ) then 
      if (num_neo_file) then
        call neo_coef( xpoint2, xcase2,  &
             y_g(ms,mt), Z_xpoint, eq_g(1,1,ms,mt),psi_axis,psi_bnd, amu_neo_prof(ms,mt), aki_neo_prof(ms,mt))
      else
        amu_neo_prof(ms,mt) = amu_neo_const
        aki_neo_prof(ms,mt) = aki_neo_const
      endif
    endif
    !======================================= NEO

  enddo
enddo

eq_zTe = eq_zTe / 2.d0  ! electron temperature  

!--------------------------------------------------- sum over the Gaussian integration points
do i=1,n_vertex_max
  do j=1,n_order+1
    ELM_p(:,:,1:n_var)  = 0
    ELM_n(:,:,1:n_var)  = 0
    ELM_k(:,:,1:n_var)  = 0
    ELM_kn(:,:,1:n_var) = 0

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

          !=======================================:parallel velocity 
          Vt0   = V_source(ms,mt)
          if (normalized_velocity_profile) then
            Vt0_x = dV_dpsi_source(ms,mt)*ps0_x
            Vt0_y = dV_dz_source(ms,mt)+dV_dpsi_source(ms,mt)*ps0_y
          else
            Omega_tor0_x = dV_dpsi_source(ms,mt)*ps0_x
            Omega_tor0_y = dV_dz_source(ms,mt)+dV_dpsi_source(ms,mt)*ps0_y
          endif
          !=======================================MB

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
          r0_corr = corr_neg_dens1(r0)
          r0_x  = (   y_t(ms,mt) * eq_s(mp,5,ms,mt) - y_s(ms,mt) * eq_t(mp,5,ms,mt) ) / xjac
          r0_y  = ( - x_t(ms,mt) * eq_s(mp,5,ms,mt) + x_s(ms,mt) * eq_t(mp,5,ms,mt) ) / xjac
          r0_p  = eq_p(mp,5,ms,mt)
          r0_s  = eq_s(mp,5,ms,mt)
          r0_t  = eq_t(mp,5,ms,mt)
          r0_ss = eq_ss(mp,5,ms,mt)
          r0_st = eq_st(mp,5,ms,mt)
          r0_tt = eq_tt(mp,5,ms,mt)

          r0_hat   = BigR**2 * r0
          r0_x_hat = 2.d0 * BigR * BigR_x  * r0 + BigR**2 * r0_x
          r0_y_hat = BigR**2 * r0_y

          T0    = eq_g(mp,6,ms,mt)
          T0_corr = corr_neg_temp1(T0)
          T0_x  = (   y_t(ms,mt) * eq_s(mp,6,ms,mt) - y_s(ms,mt) * eq_t(mp,6,ms,mt) ) / xjac
          T0_y  = ( - x_t(ms,mt) * eq_s(mp,6,ms,mt) + x_s(ms,mt) * eq_t(mp,6,ms,mt) ) / xjac
          T0_p  = eq_p(mp,6,ms,mt)
          T0_s  = eq_s(mp,6,ms,mt)
          T0_t  = eq_t(mp,6,ms,mt)
          T0_ss = eq_ss(mp,6,ms,mt)
          T0_tt = eq_tt(mp,6,ms,mt)
          T0_st = eq_st(mp,6,ms,mt)

          Vpar0    = eq_g(mp,7,ms,mt)
          Vpar0_x  = (   y_t(ms,mt) * eq_s(mp,7,ms,mt) - y_s(ms,mt) * eq_t(mp,7,ms,mt) ) / xjac
          Vpar0_y  = ( - x_t(ms,mt) * eq_s(mp,7,ms,mt) + x_s(ms,mt) * eq_t(mp,7,ms,mt) ) / xjac
          Vpar0_p  = eq_p(mp,7,ms,mt)
          Vpar0_s  = eq_s(mp,7,ms,mt)
          Vpar0_t  = eq_t(mp,7,ms,mt)
          Vpar0_ss = eq_ss(mp,7,ms,mt)
          Vpar0_st = eq_st(mp,7,ms,mt)
          Vpar0_tt = eq_tt(mp,7,ms,mt)

          P0    = r0 * T0
          P0_x  = r0_x * T0 + r0 * T0_x
          P0_y  = r0_y * T0 + r0 * T0_y
          P0_s  = r0_s * T0 + r0 * T0_s
          P0_t  = r0_t * T0 + r0 * T0_t
          P0_p  = r0_p * T0 + r0 * T0_p
          P0_ss = r0_ss * T0 + 2.d0 * r0_s * T0_s + r0 * T0_ss
          P0_tt = r0_tt * T0 + 2.d0 * r0_t * T0_t + r0 * T0_tt
          P0_st = r0_st * T0 + r0_s * T0_t + r0_t * T0_s + r0 * T0_st

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

          w0_xx = (w0_ss * y_t(ms,mt)**2 - 2.d0*w0_st * y_s(ms,mt)*y_t(ms,mt) + w0_tt * y_s(ms,mt)**2     &
                 + w0_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                              &
                 + w0_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) ) )     / xjac**2              &
                 - xjac_x * (w0_s* y_t(ms,mt) - w0_t * y_s(ms,mt))  / xjac**2

          w0_yy = (w0_ss * x_t(ms,mt)**2 - 2.d0*w0_st * x_s(ms,mt)*x_t(ms,mt) + w0_tt * x_s(ms,mt)**2     &
                 + w0_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                              &
                 + w0_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) ) )     / xjac**2              &
                 - xjac_y * (- w0_s * x_t(ms,mt) + w0_t * x_s(ms,mt) )  / xjac**2

          w0_xy = (- w0_ss * y_t(ms,mt)*x_t(ms,mt) - w0_tt * x_s(ms,mt)*y_s(ms,mt)                        &
                   + w0_st * (y_s(ms,mt)*x_t(ms,mt)  + y_t(ms,mt)*x_s(ms,mt)  )                           &
                   - w0_s  * (x_st(ms,mt)*y_t(ms,mt) - x_tt(ms,mt)*y_s(ms,mt) )                           &
                   - w0_t  * (x_st(ms,mt)*y_s(ms,mt) - x_ss(ms,mt)*y_t(ms,mt) ) )  / xjac**2              &
                   - xjac_x * (- w0_s * x_t(ms,mt) + w0_t * x_s(ms,mt) )   / xjac**2

          !----------------- simplified version of 2nd derivatives (for some unknown reason this is more stable!)
          !w0_xx = (w0_ss * y_t(ms,mt)**2 - 2.d0*w0_st * y_s(ms,mt)*y_t(ms,mt) + w0_tt * y_s(ms,mt)**2 )     / xjac**2
          !w0_yy = (w0_ss * x_t(ms,mt)**2 - 2.d0*w0_st * x_s(ms,mt)*x_t(ms,mt) + w0_tt * x_s(ms,mt)**2 )     / xjac**2

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

          T0_xx = (T0_ss * y_t(ms,mt)**2 - 2.d0*T0_st * y_s(ms,mt)*y_t(ms,mt) + T0_tt * y_s(ms,mt)**2  &
                + T0_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                            &
                + T0_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) ) )    / xjac**2             &     
                - xjac_x * (T0_s * y_t(ms,mt) - T0_t * y_s(ms,mt))  / xjac**2

          T0_yy = (T0_ss * x_t(ms,mt)**2 - 2.d0*T0_st * x_s(ms,mt)*x_t(ms,mt) + T0_tt * x_s(ms,mt)**2  &
                + T0_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                            &
                + T0_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) ) )    / xjac**2             &
                - xjac_y * (- T0_s * x_t(ms,mt) + T0_t * x_s(ms,mt) )  / xjac**2

          T0_xy = (- T0_ss * y_t(ms,mt)*x_t(ms,mt) - T0_tt * x_s(ms,mt)*y_s(ms,mt)                     &
                   + T0_st * (y_s(ms,mt)*x_t(ms,mt)  + y_t(ms,mt)*x_s(ms,mt)  )                        &
                   - T0_s  * (x_st(ms,mt)*y_t(ms,mt) - x_tt(ms,mt)*y_s(ms,mt) )                        &
                   - T0_t * (x_st(ms,mt)*y_s(ms,mt)  - x_ss(ms,mt)*y_t(ms,mt) )  )  / xjac**2          &
                - xjac_x * (- T0_s * x_t(ms,mt) + T0_t * x_s(ms,mt) )   / xjac**2

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


          P0_xx = r0_xx * T0 + 2.d0 * r0_x * T0_x + r0 * T0_xx
          P0_yy = r0_yy * T0 + 2.d0 * r0_y * T0_y + r0 * T0_yy
          P0_xy = r0_xy * T0 + r0_x * T0_y + r0_y * T0_x + r0 * T0_xy

          T0_ps0_x = T0_xx * ps0_y - T0_xy * ps0_x + T0_x * ps0_xy - T0_y * ps0_xx
          T0_ps0_y = T0_xy * ps0_y - T0_yy * ps0_x + T0_x * ps0_yy - T0_y * ps0_xy

          u0_xx = (u0_ss * y_t(ms,mt)**2 - 2.d0*u0_st * y_s(ms,mt)*y_t(ms,mt) + u0_tt * y_s(ms,mt)**2  &
                + u0_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                            &
                + u0_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) ) )   / xjac**2              &
                - xjac_x * (u0_s * y_t(ms,mt) - u0_t * y_s(ms,mt)) / xjac**2

          u0_yy = (u0_ss * x_t(ms,mt)**2 - 2.d0*u0_st * x_s(ms,mt)*x_t(ms,mt) + u0_tt * x_s(ms,mt)**2  &
                + u0_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                            &
                + u0_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) ) )      / xjac**2           &
                - xjac_y * (- u0_s * x_t(ms,mt) + u0_t * x_s(ms,mt) ) / xjac**2

          u0_xy = (- u0_ss * y_t(ms,mt)*x_t(ms,mt) - u0_tt * x_s(ms,mt)*y_s(ms,mt)                     &
                   + u0_st * (y_s(ms,mt)*x_t(ms,mt)  + y_t(ms,mt)*x_s(ms,mt)  )                        &
                   - u0_s  * (x_st(ms,mt)*y_t(ms,mt) - x_tt(ms,mt)*y_s(ms,mt) )                        &
                   - u0_t * (x_st(ms,mt)*y_s(ms,mt)  - x_ss(ms,mt)*y_t(ms,mt) )  )  / xjac**2          &
                - xjac_x * (- u0_s * x_t(ms,mt) + u0_t * x_s(ms,mt) )   / xjac**2


          delta_u_x = (   y_t(ms,mt) * delta_s(mp,2,ms,mt) - y_s(ms,mt) * delta_t(mp,2,ms,mt) ) / xjac
          delta_u_y = ( - x_t(ms,mt) * delta_s(mp,2,ms,mt) + x_s(ms,mt) * delta_t(mp,2,ms,mt) ) / xjac

          delta_ps_x = (   y_t(ms,mt) * delta_s(mp,1,ms,mt) - y_s(ms,mt) * delta_t(mp,1,ms,mt) ) / xjac
          delta_ps_y = ( - x_t(ms,mt) * delta_s(mp,1,ms,mt) + x_s(ms,mt) * delta_t(mp,1,ms,mt) ) / xjac

          ! --- Temperature dependent resistivity
          if ( eta_T_dependent .and. corr_neg_temp1(T0) <= T_max_eta) then
            eta_T     = eta   * (corr_neg_temp1(T0)/T_0)**(-1.5d0)
            deta_dT   = - eta   * (1.5d0)  * corr_neg_temp1(T0)**(-2.5d0) * T_0**(1.5d0)
            d2eta_d2T =   eta   * (3.75d0) * corr_neg_temp1(T0)**(-3.5d0) * T_0**(1.5d0)
          else if ( eta_T_dependent .and. corr_neg_temp1(T0) > T_max_eta) then
            eta_T     = eta   * (T_max_eta/T_0)**(-1.5d0)
            deta_dT   = 0.
            d2eta_d2T = 0.     
          else
            eta_T     = eta
            deta_dT   = 0.d0
            d2eta_d2T = 0.d0
          end if

          if ( eta_T_dependent .and.  xpoint2 .and. (T0 .lt. T_min) ) then
              eta_T     = eta    * (max(T0,T_min)/T_0)**(-1.5d0)
              deta_dT   = 0.d0
              d2eta_d2T = 0.d0
          end if

          ! --- Eta for ohmic heating
          if ( eta_T_dependent .and. corr_neg_temp1(T0) <= T_max_eta_ohm) then
            eta_T_ohm     = eta_ohmic   * (corr_neg_temp1(T0)/T_0)**(-1.5d0)
            deta_dT_ohm   = - eta_ohmic   * (1.5d0)  * corr_neg_temp1(T0)**(-2.5d0) * T_0**(1.5d0)
          else if ( eta_T_dependent .and. corr_neg_temp1(T0) > T_max_eta_ohm) then
            eta_T_ohm     = eta_ohmic   * (T_max_eta_ohm/T_0)**(-1.5d0)
            deta_dT_ohm   = 0.    
          else
            eta_T_ohm     = eta_ohmic
            deta_dT_ohm   = 0.d0
          end if

          ! --- Temperature dependent viscosity
          if ( visco_T_dependent ) then
            visco_T     =   visco * (T0_corr/T_0)**(-1.5d0)
            dvisco_dT   = - visco * (1.5d0)  * T0_corr**(-2.5d0) * T_0**(1.5d0)
            d2visco_dT2 =   visco * (3.75d0) * T0_corr**(-3.5d0) * T_0**(1.5d0)
            if ( xpoint2 .and. (T0 .lt. T_min) ) then
              visco_T     = visco  * (max(T0,T_min)/T_0)**(-1.5d0)
              dvisco_dT   = 0.d0
              d2visco_dT2 = 0.d0
            endif
          else
            visco_T     = visco
            dvisco_dT   = 0.d0
            d2visco_dT2 = 0.d0
          end if

          ! --- Temperature dependent parallel heat diffusivity
          if ( ZKpar_T_dependent ) then
            ZKpar_T   = ZK_par * (T0_corr/T_0)**(+2.5d0)              ! temperature dependent parallel conductivity
            dZKpar_dT = ZK_par * (2.5d0)  * T0_corr**(+1.5d0) * T_0**(-2.5d0)
            if (ZKpar_T .gt. ZK_par_max) then
              ZKpar_T   = Zk_par_max
              dZKpar_dT = 0.d0
            endif
            if ( xpoint2 .and. (T0 .lt. T_min) ) then
              ZKpar_T   = ZK_par * (max(T0,T_min)/T_0)**(+2.5d0)
              dZKpar_dT = 0.d0
            endif
          else
            ZKpar_T   = ZK_par                                            ! parallel conductivity
            dZKpar_dT = 0.d0
          endif

          ! --- Diamagnetic viscosity
          if (Wdia) then
            W_dia = + tauIC /r0_corr    * (p0_xx + p0_x/bigR + p0_yy) &
                    - tauIC /r0_corr**2 * (r0_x*p0_x + r0_y*p0_y)
          else
            W_dia = 0.d0
          endif
          
          ! --- Temperature dependent hyper-resistivity. There is no physical
          ! reason for this dependence whatsoever, this is just to keep a constant
          ! ratio between the resistivity and hyper-resistivity.
          if ( eta_num_T_dependent .and. T0_corr <= T_max_eta) then
             eta_num_T = eta_num * (T0_corr/T_0)**(-1.5d0)
             deta_num_dT =  - eta_num * (1.5d0)  * T0_corr**(-2.5d0) * T_0**(1.5d0)
          else if (eta_num_T_dependent .and. T0_corr > T_max_eta) then
             eta_num_T = eta_num * (T_max_eta/T_0)**(-1.5d0)
             deta_num_dT = 0.
          else
             eta_num_T = eta_num
             deta_num_dT = 0.
          end if

          ! --- Same for the hyper-viscosity.
          if ( visco_num_T_dependent .and. T0_corr <= T_max_eta) then
             visco_num_T = visco_num * (T0_corr/T_0)**(-1.5d0)
             dvisco_num_dT =  - visco_num * (1.5d0)  * T0_corr**(-2.5d0) * T_0**(1.5d0) 
          else if (visco_num_T_dependent .and. T0_corr > T_max_eta) then
             visco_num_T = visco_num * (T_max_eta/T_0)**(-1.5d0)
             dvisco_num_dT = 0.
          else
             visco_num_T = visco_num
             dvisco_num_dT = 0.
          end if

          psi_norm = get_psi_n( ps0, y_g(ms,mt))

          ! --- Bootstrap current 
          if (bootstrap) then
            ! --- Full Sauter formula
            Ti0   = T0   / 2.d0 ; Te0        = T0   / 2.d0
            Ti0_x = T0_x / 2.d0 ; Te0_x = T0_x / 2.d0
            Ti0_y = T0_y / 2.d0 ; Te0_y = T0_y / 2.d0
            call bootstrap_current(bigR, y_g(ms,mt),                     &
                                   R_axis,   Z_axis,   psi_axis,         &
                                   R_xpoint, Z_xpoint, psi_bnd, psi_norm,&
                                   ps0, ps0_x, ps0_y,                    &
                                   r0,  r0_x,  r0_y,                     &
                                   Ti0, Ti0_x, Ti0_y,                    &
                                   Te0, Te0_x, Te0_y,                  Jb)
            
            
            ! --- Full Sauter formula for initial profiles
            
            zTi   = eq_zTe(ms,mt)  ! Dividing by 2.0 not necessary because it's been done above already            
            zTi_x = dT_dpsi(ms,mt) * ps0_x / 2.d0
            zTi_y = dT_dpsi(ms,mt) * ps0_y / 2.d0
            zTe   = zTi  
            zTe_x = zTi_x
            zTe_y = zTi_y
            zn_x  = dn_dpsi(ms,mt) * ps0_x
            zn_y  = dn_dpsi(ms,mt) * ps0_y
            call bootstrap_current(bigR, y_g(ms,mt) ,                      &
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



          D_prof  = get_dperp (psi_norm)
          ZK_prof = get_zkperp(psi_norm)

          ! --- Increase diffusivity if very small density/temperature
          if (xpoint2) then
            if (r0 .lt. D_prof_neg_thresh)  then
              D_prof  = D_prof_neg
            endif
            if (T0 .lt. ZK_prof_neg_thresh) then
              ZK_prof = ZK_prof_neg
            endif
          endif

          phi       = 2.d0*PI*float(mp-1)/float(n_plane) / float(n_period)
          delta_phi = 2.d0*PI/float(n_plane) / float(n_period)

          source_pellet = 0.d0
          source_volume = 0.d0

          if (use_pellet) then

            call pellet_source2(pellet_amplitude,pellet_R,pellet_Z,pellet_psi,pellet_phi, &
                                pellet_radius, pellet_delta_psi, pellet_sig, pellet_length, pellet_ellipse, pellet_theta, &
                                x_g(ms,mt),y_g(ms,mt), ps0, phi, r0_corr, T0_corr/2.d0, &
                                central_density, pellet_particles, pellet_density, total_pellet_volume, &
                                source_pellet, source_volume)
          endif

          do im=n_tor_start, n_tor_end

            v   =  H(i,j,ms,mt) * element%size(i,j) * HHZ(im,mp)
            v_x = (  y_t(ms,mt) * h_s(i,j,ms,mt) - y_s(ms,mt) * h_t(i,j,ms,mt) ) * element%size(i,j) / xjac * HHZ(im,mp)
            v_y = (- x_t(ms,mt) * h_s(i,j,ms,mt) + x_s(ms,mt) * h_t(i,j,ms,mt) ) * element%size(i,j) / xjac * HHZ(im,mp)

            v_s = h_s(i,j,ms,mt) * element%size(i,j) * HHZ(im,mp)
            v_t = h_t(i,j,ms,mt) * element%size(i,j) * HHZ(im,mp)
            v_p = H(i,j,ms,mt)   * element%size(i,j) * HHZ_p(im,mp)

            v_ss = h_ss(i,j,ms,mt) * element%size(i,j) * HHZ(im,mp)
            v_tt = h_tt(i,j,ms,mt) * element%size(i,j) * HHZ(im,mp)
            v_st = h_st(i,j,ms,mt) * element%size(i,j) * HHZ(im,mp)

            v_xx = (v_ss * y_t(ms,mt)**2 - 2.d0*v_st * y_s(ms,mt)*y_t(ms,mt) + v_tt * y_s(ms,mt)**2    &
                   + v_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                          &
                   + v_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) ) )  / xjac**2             &   
                   - xjac_x * (v_s * y_t(ms,mt) - v_t * y_s(ms,mt)) / xjac**2

            v_yy = (v_ss * x_t(ms,mt)**2 - 2.d0*v_st * x_s(ms,mt)*x_t(ms,mt) + v_tt * x_s(ms,mt)**2    &
                   + v_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                          &
                   + v_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) ) )  / xjac**2             &   
                   - xjac_y * (- v_s * x_t(ms,mt) + v_t * x_s(ms,mt) ) / xjac**2

            v_xy = (- v_ss * y_t(ms,mt)*x_t(ms,mt) - v_tt * x_s(ms,mt)*y_s(ms,mt)                      &
                   + v_st * (y_s(ms,mt)*x_t(ms,mt)  + y_t(ms,mt)*x_s(ms,mt)  )                         &
                   - v_s  * (x_st(ms,mt)*y_t(ms,mt) - x_tt(ms,mt)*y_s(ms,mt) )                         &
                   - v_t  * (x_st(ms,mt)*y_s(ms,mt) - x_ss(ms,mt)*y_t(ms,mt) )  )  / xjac**2           &
                   - xjac_x * (- v_s * x_t(ms,mt) + v_t * x_s(ms,mt) )   / xjac**2

            Bgrad_rho_star   = ( v_x  * ps0_y - v_y  * ps0_x ) / BigR                         
            Bgrad_rho_k_star = ( F0 / BigR * v_p )           / BigR                           
            Bgrad_rho        = ( F0 / BigR * r0_p +  r0_x * ps0_y - r0_y * ps0_x ) / BigR

            Bgrad_T_star     = ( v_x  * ps0_y - v_y  * ps0_x ) / BigR                         
            Bgrad_T_k_star   = ( F0 / BigR * v_p           ) / BigR                           

            Bgrad_T          = ( F0 / BigR * T0_p +  T0_x * ps0_y - T0_y * ps0_x ) / BigR

            BB2              = (F0*F0 + ps0_x * ps0_x + ps0_y * ps0_y )/BigR**2
            Btheta2          = (ps0_x * ps0_x + ps0_y * ps0_y )/BigR**2

            v_ps0_x  = v_xx  * ps0_y - v_xy  * ps0_x + v_x  * ps0_xy - v_y * ps0_xx
            v_ps0_y  = v_xy  * ps0_y - v_yy  * ps0_x + v_x  * ps0_yy - v_y * ps0_xy

            !###################################################################################################
            !#  equation 1   (induction equation)                                                              #
            !###################################################################################################

            rhs_ij_1 =   v * eta_T  * (zj0 - current_source(ms,mt) - Jb)/ BigR    * xjac * tstep &
                         + v * (ps0_s * u0_t - ps0_t * u0_s)                        * tstep &
                         - v * eps_cyl * F0 / BigR  * u0_p                   * xjac * tstep &
                         + eta_num_T * (v_x * zj0_x + v_y * zj0_y)           * xjac * tstep &

                         - v * tauIC/(r0_corr*BB2) * F0**2/BigR**2 * (ps0_s * p0_t - ps0_t * p0_s) * tstep &
                         + v * tauIC/(r0_corr*BB2) * F0**3/BigR**3 * eps_cyl * p0_p * xjac * tstep &

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

                         - v * tauIC * BigR**4 * (p0_s * w0_t - p0_t * w0_s)        * tstep &

                         - tauIC * BigR**3 * p0_y * (v_x* u0_x + v_y * u0_y) * xjac * tstep &

                         - v * tauIC * BigR**4 * (u0_xy * (p0_xx - p0_yy) - p0_xy * (u0_xx - u0_yy) ) * xjac * tstep &

                         ! --- Diamagnetic viscosity
                         + dvisco_dT * bigR * W_dia * (v_x*T0_x + v_y*T0_y)    * xjac * tstep  &
                         + visco_T   * bigR * W_dia * (v_xx + v_x/bigR + v_yy) * xjac * tstep  &

                         + BigR**3 * (particle_source(ms,mt) + source_pellet) * (v_x * u0_x + v_y * u0_y) * xjac* tstep &

                         - zeta * BigR * r0_hat * (v_x * delta_u_x + v_y * delta_u_y) * xjac
            
            !------------------------------------------------------------------------ NEO
            if (NEO) then
              rhs_ij_2 =  rhs_ij_2  + amu_neo_prof(ms,mt)*BB2/((Btheta2+epsil)**2)*(ps0_x*v_x+ps0_y*v_y)* &
                   (r0*(ps0_x*u0_x+ps0_y*u0_y)+&
                   tauIC*(ps0_x*P0_x+ps0_y*P0_y) &
                   +aki_neo_prof(ms,mt)*tauIC*r0*(ps0_x*T0_x+ps0_y*T0_y)&
                   -r0*Vpar0*Btheta2)*xjac*tstep*BigR 
            endif
            !------------------------------------------------------------------------ NEO

            !###################################################################################################
            !#  equation 3   (current definition)                                                              #
            !###################################################################################################

            rhs_ij_3 = - ( v_x * ps0_x  + v_y * ps0_y + v*zj0) / BigR * xjac

            !###################################################################################################
            !#  equation 4   (vorticity definition)                                                            #
            !###################################################################################################

            rhs_ij_4 = - ( v_x * u0_x   + v_y * u0_y  + v*w0)  * BigR * xjac 

            !###################################################################################################
            !#  equation 5   (density equation)                                                                #
            !###################################################################################################

            rhs_ij_5   = v * BigR * (particle_source(ms,mt) + source_pellet)                      * xjac * tstep &
                       + v * BigR**2 * ( r0_s * u0_t - r0_t * u0_s)                                      * tstep &
                       + v * 2.d0 * BigR * r0 * u0_y                                              * xjac * tstep &
                       - (D_par-D_prof) * BigR / BB2 * Bgrad_rho_star * Bgrad_rho                 * xjac * tstep &
                       - D_prof * BigR  * (v_x*r0_x + v_y*r0_y                                  ) * xjac * tstep &
                       - v * F0 / BigR * Vpar0 * r0_p                                             * xjac * tstep &
                       - v * Vpar0 * (r0_s * ps0_t - r0_t * ps0_s)                                       * tstep &
                       - v * F0 / BigR * r0 * vpar0_p                                             * xjac * tstep &
                       - v * r0 * (vpar0_s * ps0_t - vpar0_t * ps0_s)                                    * tstep &

                       + v * 2.d0 * tauIC * p0_y * BigR                                           * xjac * tstep &

                       + zeta * v * delta_g(mp,5,ms,mt) * BigR                                    * xjac  &

                       - D_perp_num * (v_xx + v_x/Bigr + v_yy)*(r0_xx + r0_x/Bigr + r0_yy) * BigR * xjac * tstep &

                       - TG_num5 * 0.25d0 * BigR**3 * (r0_x * u0_y - r0_y * u0_x)                                &
                                                    * ( v_x * u0_y - v_y * u0_x) * xjac * tstep * tstep          &
                       - TG_num5 * 0.25d0 / BigR * vpar0**2                                                      &
                                 * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                              &
                                 * ( v_x * ps0_y -  v_y * ps0_x                   ) * xjac * tstep * tstep


            rhs_ij_5_k =  - (D_par-D_prof) * BigR / BB2 * Bgrad_rho_k_star * Bgrad_rho            * xjac * tstep &
                          - D_prof * BigR  * (                  v_p*r0_p * eps_cyl**2 /BigR**2 )  * xjac * tstep &

                       - TG_num5 * 0.25d0 / BigR * vpar0**2 &
                                 * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                              &
                                 * (                            + F0 / BigR * v_p) * xjac * tstep * tstep

            !###################################################################################################
            !#  equation 6   (energy equation)                                                                 #
            !###################################################################################################

            rhs_ij_6 =   v * BigR * heat_source(ms,mt)                                    * xjac * tstep &
            
                       + v * r0 * BigR**2 * ( T0_s * u0_t - T0_t * u0_s)                         * tstep &
                       + v * T0 * BigR**2 * ( r0_s * u0_t - r0_t * u0_s)                         * tstep &

                       + v * r0 * T0 * 2.d0* GAMMA * BigR * u0_y                          * xjac * tstep &

                       - v * r0 * F0 / BigR * Vpar0 * T0_p                                * xjac * tstep &
                       - v * T0 * F0 / BigR * Vpar0 * r0_p                                * xjac * tstep &

                       - v * r0 * Vpar0 * (T0_s * ps0_t - T0_t * ps0_s)                          * tstep &
                       - v * T0 * Vpar0 * (r0_s * ps0_t - r0_t * ps0_s)                          * tstep &

                       - v * r0 * T0 * GAMMA * (vpar0_s * ps0_t - vpar0_t * ps0_s)               * tstep &
                       - v * r0 * T0 * GAMMA * F0 / BigR * vpar0_p                        * xjac * tstep &

                       - (ZKpar_T-ZK_prof) * BigR / BB2 * Bgrad_T_star * Bgrad_T     * xjac * tstep &
                       - ZK_prof * BigR * (v_x*T0_x + v_y*T0_y                     ) * xjac * tstep &
 
                       - ZK_perp_num  *  (v_xx + v_x/Bigr + v_yy)*(T0_xx + T0_x/Bigr + T0_yy) * BigR * xjac * tstep &

                       - TG_num6 * 0.25d0 * BigR**3 * T0 * (r0_x * u0_y - r0_y * u0_x) &
                                          * ( v_x * u0_y - v_y * u0_x) * xjac * tstep * tstep  &
                       - TG_num6 * 0.25d0 * BigR**3 * r0 * (T0_x * u0_y - T0_y * u0_x) &
                                          * ( v_x * u0_y - v_y * u0_x) * xjac * tstep * tstep  &

                       - TG_num6 * 0.25d0 / BigR * vpar0**2 &
                                 * T0 * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                         &
                                 * ( v_x * ps0_y -  v_y * ps0_x                        ) * xjac * tstep * tstep  &
                       - TG_num6 * 0.25d0 / BigR * vpar0**2 &
                                 * r0 * (T0_x * ps0_y - T0_y * ps0_x + F0 / BigR * T0_p)                         &
                                 * ( v_x * ps0_y -  v_y * ps0_x                        ) * xjac * tstep * tstep  &

                       + zeta * v * r0 * delta_g(mp,6,ms,mt) * BigR                       * xjac &
                       + zeta * v * T0 * delta_g(mp,5,ms,mt) * BigR                       * xjac &

                       + v * (gamma-1.d0) * eta_T_ohm * (zj0 / BigR)**2.d0         * BigR * xjac  * tstep


            rhs_ij_6_k =  - (ZKpar_T-ZK_prof) * BigR / BB2 * Bgrad_T_k_star * Bgrad_T  * xjac * tstep &
                          - ZK_prof * BigR * (                + v_p*T0_p /BigR**2 )   * xjac * tstep  &

                         - TG_num6 * 0.25d0 / BigR * vpar0**2 &
                                 * T0 * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                         &
                                 * (                                 + F0 / BigR * v_p) * xjac * tstep * tstep   &
                         - TG_num6 * 0.25d0 / BigR * vpar0**2 &
                                 * r0 * (T0_x * ps0_y - T0_y * ps0_x + F0 / BigR * T0_p)                         &
                                 * (                                 + F0 / BigR * v_p) * xjac * tstep * tstep

            !###################################################################################################
            !#  equation 7   (parallel velocity equation)                                                      #
            !###################################################################################################

            rhs_ij_7 = - v * F0 / BigR * P0_p                                              * xjac * tstep &
                       - v * (P0_s * ps0_t - P0_t * ps0_s)                                        * tstep &

                     - v*(particle_source(ms,mt) + source_pellet) * vpar0 * BB2   * BigR * xjac * tstep &

                     - 0.5d0 * r0 * vpar0**2 * BB2 * (ps0_s * v_t - ps0_t * v_s)                * tstep &
                     - 0.5d0 * v  * vpar0**2 * BB2 * (ps0_s * r0_t - ps0_t * r0_s)              * tstep &
                     + 0.5d0 * v  * vpar0**2 * BB2 * F0 / BigR * r0_p                    * xjac * tstep &

                     - visco_par_num * (v_xx + v_x/Bigr + v_yy)*(vpar0_xx + vpar0_x/Bigr + vpar0_yy) * BigR * xjac * tstep &

                     + zeta * v * delta_g(mp,7,ms,mt) * R0 * F0**2 / BigR                        * xjac &
                     + zeta * v * r0 * vpar0 * (ps0_x * delta_ps_x + ps0_y * delta_ps_y) / BigR  * xjac &

             - TG_NUM7 * 0.25d0 * r0 * Vpar0**2 * BB2 &
                       * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac + F0 / BigR * vpar0_p) / BigR  &
                       * (-(ps0_s * v_t     - ps0_t * v_s)    /xjac                      )  * xjac * tstep * tstep &
             - TG_NUM7 * 0.25d0 * v  * Vpar0**2 * BB2 &
                       * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac + F0 / BigR * vpar0_p) / BigR  &
                       * (-(ps0_s * r0_t    - ps0_t * r0_s)   /xjac + F0 / BigR * r0_p)  * xjac * tstep * tstep 

            !====================================================parallel velocity source:
            if (normalized_velocity_profile) then
              rhs_ij_7 = rhs_ij_7  - visco_par * (v_x * (vpar0_x-Vt0_x) + v_y * (vpar0_y-Vt0_y)) * BigR* xjac * tstep 
            else
              rhs_ij_7 = rhs_ij_7 - visco_par * (  v_x * (vpar0_x * F0**2 / BigR**2 - 2 * vpar0 * F0**2 / BigR**3  - 2 * PI * F0 * Omega_tor0_x ) & 
                                                 + v_y * (vpar0_y * F0**2 / BigR**2 -2 * PI * F0 * Omega_tor0_y) ) * BigR* xjac * tstep 
            endif
            !=============================================================================

            ! ------------------------------------------------------ NEO
            if (NEO) then
              rhs_ij_7 =  rhs_ij_7  + amu_neo_prof(ms,mt)*BB2/(Btheta2+epsil)*v*(r0*(ps0_x*u0_x+ps0_y*u0_y)+&
                   tauIC*(ps0_x*P0_x+ps0_y*P0_y) &
                   +aki_neo_prof(ms,mt)*tauIC*r0*(ps0_x*T0_x+ps0_y*T0_y)&
                   -r0*Vpar0*Btheta2)*xjac*tstep*BigR
            endif
            ! ------------------------------------------------------ NEO


            rhs_ij_7_k = + 0.5d0 * r0 * vpar0**2 * BB2 * F0 / BigR * v_p                     * xjac * tstep &

               - TG_NUM7 * 0.25d0 * r0 * Vpar0**2 * BB2 &
                         * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac + F0 / BigR * vpar0_p) / BigR  &
                         * (                                          + F0 / BigR * v_p)  * xjac * tstep * tstep 

            !###################################################################################################
            !#  RHS equations end                                                                                  #
            !###################################################################################################

            if (use_fft) then

              index_ij =       n_var*(n_order+1)*(i-1) +       n_var*(j-1) + 1
            else
              index_ij = n_tor_local*n_var*(n_order+1)*(i-1) + n_tor_local * n_var * (j-1) + im - n_tor_start +1 
            endif


            ! --- Fill up the matrix
            if (use_fft) then

              ij1 = index_ij
              ij2 = index_ij + 1
              ij3 = index_ij + 2
              ij4 = index_ij + 3
              ij5 = index_ij + 4
              ij6 = index_ij + 5
              ij7 = index_ij + 6
              RHS_p(mp,ij1) = RHS_p(mp,ij1) + rhs_ij_1 * wst
              RHS_p(mp,ij2) = RHS_p(mp,ij2) + rhs_ij_2 * wst
              RHS_p(mp,ij3) = RHS_p(mp,ij3) + rhs_ij_3 * wst
              RHS_p(mp,ij4) = RHS_p(mp,ij4) + rhs_ij_4 * wst
              RHS_p(mp,ij5) = RHS_p(mp,ij5) + rhs_ij_5 * wst
              RHS_p(mp,ij6) = RHS_p(mp,ij6) + rhs_ij_6 * wst
              RHS_p(mp,ij7) = RHS_p(mp,ij7) + rhs_ij_7 * wst

              RHS_k(mp,ij5) = RHS_k(mp,ij5) + rhs_ij_5_k * wst
              RHS_k(mp,ij6) = RHS_k(mp,ij6) + rhs_ij_6_k * wst
              RHS_k(mp,ij7) = RHS_k(mp,ij7) + rhs_ij_7_k * wst
            else
              ij1 = index_ij
              ij2 = index_ij + 1*(n_tor_end - n_tor_start + 1)
              ij3 = index_ij + 2*(n_tor_end - n_tor_start + 1)
              ij4 = index_ij + 3*(n_tor_end - n_tor_start + 1)
              ij5 = index_ij + 4*(n_tor_end - n_tor_start + 1)
              ij6 = index_ij + 5*(n_tor_end - n_tor_start + 1)
              ij7 = index_ij + 6*(n_tor_end - n_tor_start + 1)
              RHS(ij1) = RHS(ij1) + (rhs_ij_1 + rhs_ij_1_k) * wst
              RHS(ij2) = RHS(ij2) + (rhs_ij_2 + rhs_ij_2_k) * wst
              RHS(ij3) = RHS(ij3) + (rhs_ij_3 + rhs_ij_3_k) * wst
              RHS(ij4) = RHS(ij4) + (rhs_ij_4 + rhs_ij_4_k) * wst
              RHS(ij5) = RHS(ij5) + (rhs_ij_5 + rhs_ij_5_k) * wst
              RHS(ij6) = RHS(ij6) + (rhs_ij_6 + rhs_ij_6_k) * wst
              RHS(ij7) = RHS(ij7) + (rhs_ij_7 + rhs_ij_7_k) * wst
            endif

            do k=1,n_vertex_max

              do l=1,n_order+1

                do in = n_tor_start, n_tor_end

                  psi   = H(k,l,ms,mt) * element%size(k,l) * HHZ(in,mp)

                  psi_x = (   y_t(ms,mt) * h_s(k,l,ms,mt) - y_s(ms,mt) * h_t(k,l,ms,mt) ) / xjac * element%size(k,l) * HHZ(in,mp)
                  psi_y = ( - x_t(ms,mt) * h_s(k,l,ms,mt) + x_s(ms,mt) * h_t(k,l,ms,mt) ) / xjac * element%size(k,l) * HHZ(in,mp)

                  psi_p  = H(k,l,ms,mt)   * element%size(k,l) * HHZ_p(in,mp)
                  psi_s  = h_s(k,l,ms,mt) * element%size(k,l) * HHZ(in,mp)
                  psi_t  = h_t(k,l,ms,mt) * element%size(k,l) * HHZ(in,mp)
                  psi_ss = h_ss(k,l,ms,mt) * element%size(k,l) * HHZ(in,mp)
                  psi_tt = h_tt(k,l,ms,mt) * element%size(k,l) * HHZ(in,mp)
                  psi_st = h_st(k,l,ms,mt) * element%size(k,l) * HHZ(in,mp)

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

                  u    = psi    ;  zj    = psi    ;  w    = psi    ; rho    = psi    ;  T    = psi    ; vpar    = psi
                  u_x  = psi_x  ;  zj_x  = psi_x  ;  w_x  = psi_x  ; rho_x  = psi_x  ;  T_x  = psi_x  ; vpar_x  = psi_x
                  u_y  = psi_y  ;  zj_y  = psi_y  ;  w_y  = psi_y  ; rho_y  = psi_y  ;  T_y  = psi_y  ; vpar_y  = psi_y
                  u_p  = psi_p  ;  zj_p  = psi_p  ;  w_p  = psi_p  ; rho_p  = psi_p  ;  T_p  = psi_p  ; vpar_p  = psi_p
                  u_s  = psi_s  ;  zj_s  = psi_s  ;  w_s  = psi_s  ; rho_s  = psi_s  ;  T_s  = psi_s  ; vpar_s  = psi_s
                  u_t  = psi_t  ;  zj_t  = psi_t  ;  w_t  = psi_t  ; rho_t  = psi_t  ;  T_t  = psi_t  ; vpar_t  = psi_t
                  u_ss = psi_ss ;  zj_ss = psi_ss ;  w_ss = psi_ss ; rho_ss = psi_ss ;  T_ss = psi_ss ; vpar_ss = psi_ss
                  u_tt = psi_tt ;  zj_tt = psi_tt ;  w_tt = psi_tt ; rho_tt = psi_tt ;  T_tt = psi_tt ; vpar_tt = psi_tt
                  u_st = psi_st ;  zj_st = psi_st ;  w_st = psi_st ; rho_st = psi_st ;  T_st = psi_st ; vpar_st = psi_st

                  u_xx = psi_xx ;                    w_xx = psi_xx ; rho_xx = psi_xx ;  T_xx = psi_xx ; vpar_xx = psi_xx
                  u_yy = psi_yy ;                    w_yy = psi_yy ; rho_yy = psi_yy ;  T_yy = psi_yy ; vpar_yy = psi_yy
                  u_xy = psi_xy ;                    w_xy = psi_xy ; rho_xy = psi_xy ;  T_xy = psi_xy ; vpar_xy = psi_xy

                  !w_xx = (psi_ss * y_t(ms,mt)**2 - 2.d0*psi_st * y_s(ms,mt)*y_t(ms,mt) + psi_tt * y_s(ms,mt)**2 ) / xjac**2
                  !w_yy = (psi_ss * x_t(ms,mt)**2 - 2.d0*psi_st * x_s(ms,mt)*x_t(ms,mt) + psi_tt * x_s(ms,mt)**2 ) / xjac**2
                  !w_xy = psi_xy

                  rho_hat   = BigR**2 * rho
                  rho_x_hat = 2.d0 * BigR * BigR_x  * rho + BigR**2 * rho_x
                  rho_y_hat = BigR**2 * rho_y

                  Btheta2_psi  = 2.d0 * (psi_x * ps0_x + psi_y * ps0_y ) /BigR**2
                  Bgrad_rho_star_psi = ( v_x  * psi_y - v_y  * psi_x ) / BigR
                  Bgrad_rho_psi      = ( r0_x * psi_y - r0_y * psi_x ) / BigR
                  Bgrad_rho_rho      = ( rho_x * ps0_y - rho_y * ps0_x ) / BigR
                  Bgrad_rho_rho_n    = ( F0 / BigR * rho_p ) / BigR
                  BB2_psi            = 2.d0 * (psi_x * ps0_x + psi_y * ps0_y ) /BigR**2

                  P0_x_rho  = rho_x*T0 + rho*T0_x
                  P0_xx_rho = rho_xx*T0 + 2.0*rho_x*T0_x + rho*T0_xx
                  P0_y_rho  = rho_y*T0 + rho*T0_y
                  P0_yy_rho = rho_yy*T0 + 2.0*rho_y*T0_y + rho*T0_yy
                  P0_xy_rho = rho_xy*T0 + rho*T0_xy + rho_x*T0_y + rho_y*T0_x
                 
                  P0_x_T  = r0_x*T + r0*T_x
                  P0_xx_T = r0_xx*T + 2.0*r0_x*T_x + r0*T_xx
                  P0_y_T  = r0_y*T + r0*T_y
                  P0_yy_T = r0_yy*T + 2.0*r0_y*T_y + r0*T_yy
                  P0_xy_T = r0_xy*T + r0*T_xy + r0_x*T_y + r0_y*T_x
  
                  if (Wdia) then
                    W_dia_rho = - tauIC *     rho/r0_corr**2 * (p0_xx     + p0_x    /bigR + p0_yy    ) &
                                + tauIC          /r0_corr    * (p0_xx_rho + p0_x_rho/bigR + p0_yy_rho) &
                                + tauIC * 2.0*rho/r0_corr**3 * (r0_x *p0_x     + r0_y *p0_y    )    &
                                - tauIC          /r0_corr**2 * (rho_x*p0_x     + rho_y*p0_y    )    &
                                - tauIC          /r0_corr**2 * (r0_x *p0_x_rho + r0_y *p0_y_rho)
                    W_dia_T   = + tauIC          /r0_corr    * (p0_xx_T + p0_x_T/bigR + p0_yy_T) &
                                - tauIC          /r0_corr**2 * (r0_x  *p0_x_T + r0_y  *p0_y_T)
                  else
                    W_dia_rho = 0.d0
                    W_dia_T   = 0.d0
                  endif
  
                  if (normalized_velocity_profile) then
                     Vt_x_psi    = dV_dpsi_source(ms,mt) * psi_x
                     Vt_y_psi    = dV_dpsi_source(ms,mt) * psi_y
                  else
                     Omega_tor_x_psi = dV_dpsi_source(ms,mt)*psi_x
                     Omega_tor_y_psi = dV_dpsi_source(ms,mt)*psi_y
                  endif

                  !###################################################################################################
                  !#  equation 1   (induction equation)                                                              #
                  !###################################################################################################

                  amat_11 = v * psi / BigR * xjac * (1.d0 + zeta)                                              &
                          - v * (psi_s * u0_t - psi_t * u0_s)                                  * theta * tstep &
                          + v * tauIC/(r0_corr*BB2) * F0**2/BigR**2 * (psi_s * p0_t - psi_t * p0_s) * theta * tstep &
                     - v * tauIC/(r0_corr*BB2**2) * BB2_psi * F0**2/BigR**2 * (ps0_x*p0_y - ps0_y*p0_x)              * xjac * theta * tstep &
                     + v * tauIC/(r0_corr*BB2**2) * BB2_psi * F0**3/BigR**3 * eps_cyl * p0_p                 * xjac * theta * tstep

                  amat_12 = -  v * (ps0_s * u_t - ps0_t * u_s)                             * theta * tstep

                  amat_12_n = +  eps_cyl * F0 / BigR * v * u_p * xjac                      * theta * tstep

                  amat_13 = - eta_num_T * (v_x * zj_x + v_y * zj_y)                 * xjac * theta * tstep  &
                            - eta_T * v * zj / BigR                                 * xjac * theta * tstep

                  amat_15 = + v * tauIC/(r0_corr*BB2) * F0**2/BigR**2 * T0  * (ps0_s * rho_t - ps0_t * rho_s) * theta * tstep &
                            + v * tauIC/(r0_corr*BB2) * F0**2/BigR**2 * rho * (ps0_s * T0_t  - ps0_t * T0_s)  * theta * tstep &
                            - v * tauIC/(r0_corr*BB2) * F0**3/BigR**3 * eps_cyl * rho * T0_p  * xjac          * theta * tstep &

                            - v * tauIC * rho /(r0_corr**2 * BB2) * F0**2/BigR**2 * (ps0_s * p0_t - ps0_t * p0_s) * theta * tstep &
                            + v * tauIC * rho /(r0_corr**2 * BB2) * F0**3/BigR**3 * eps_cyl * p0_p * xjac         * theta * tstep 

                  amat_15_n = - v * tauIC/(r0_corr*BB2) * F0**3/BigR**3 * eps_cyl * T0  * rho_p * xjac * theta * tstep 

                  amat_16 = - deta_dT * v * T * (zj0 - current_source(ms,mt) - Jb)/ BigR * xjac         * theta * tstep &
                          + v * tauIC/(r0_corr*BB2) * F0**2/BigR**2 * r0 * (ps0_s * T_t  - ps0_t * T_s) * theta * tstep &
                          + v * tauIC/(r0_corr*BB2) * F0**2/BigR**2 * T  * (ps0_s * r0_t - ps0_t * r0_s)* theta * tstep &
                          - v * tauIC/(r0_corr*BB2) * F0**3/BigR**3 * eps_cyl * T  * r0_p * xjac        * theta * tstep 

                  amat_16_n = - v * tauIC/(r0_corr*BB2) * F0**3/BigR**3 * eps_cyl * r0 * T_p  * xjac * theta * tstep 

                  !###################################################################################################
                  !#  equation 2   (perpendicular momentum equation)                                                 #
                  !###################################################################################################

                  amat_21 = - v * (psi_s * zj0_t - psi_t * zj0_s)                          * theta * tstep

                  ! ------------------------------------------------------ NEO
                  if (NEO) then
                    amat_21 = amat_21 &
                         -amu_neo_prof(ms,mt)*BB2/((Btheta2+epsil)**2)*(psi_x*v_x+psi_y*v_y)*&
                         (r0*(ps0_x*u0_x+ps0_y*u0_y)+&
                         tauIC*(ps0_x*P0_x+ps0_y*P0_y) &
                         +aki_neo_prof(ms,mt)*tauIC*r0*(ps0_x*T0_x+ps0_y*T0_y)&
                         -r0*Vpar0*Btheta2)*BigR*xjac*theta*tstep &
                         
                         -amu_neo_prof(ms,mt)*BB2/((Btheta2+epsil)**2)*(ps0_x*v_x+ps0_y*v_y)*(r0*(psi_x*u0_x+psi_y*u0_y)+&
                         tauIC*(psi_x*P0_x+psi_y*P0_y) &
                         +aki_neo_prof(ms,mt)*tauIC*r0*(psi_x*T0_x+psi_y*T0_y))*BigR*xjac*theta*tstep &

                         ! ========= linearization of 1/(Btheta2**i) , i=2 or 1
                         -amu_neo_prof(ms,mt)*BB2*(-2*Btheta2_psi)/((Btheta2+epsil)**3)*(ps0_x*v_x+ps0_y*v_y)*&
                         (r0*(ps0_x*u0_x+ps0_y*u0_y)+&
                         tauIC*(ps0_x*P0_x+ps0_y*P0_y) &
                         +aki_neo_prof(ms,mt)*tauIC*r0*(ps0_x*T0_x+ps0_y*T0_y))&
                         *BigR*xjac*theta*tstep &

                         +amu_neo_prof(ms,mt)*BB2*(-Btheta2_psi)/((Btheta2+epsil)**2)*r0*vpar0*(ps0_x*v_x+ps0_y*v_y)&
                         *BigR*xjac*tstep*theta 
                  endif
                  ! ------------------------------------------------------ NEO

                  amat_22 = - BigR**3 * r0_corr * (v_x * u_x + v_y * u_y) * xjac * (1.d0 + zeta)                                 &
                            + r0_hat * BigR**2 * w0 * (v_s * u_t  - v_t  * u_s)                              * theta * tstep &
                            + BigR**2 * (u_x * u0_x + u_y * u0_y) * (v_x * r0_y_hat - v_y * r0_x_hat) * xjac * theta * tstep &
                            
                            + tauIC * BigR**3 * p0_y * (v_x* u_x + v_y * u_y)                         * xjac * theta * tstep &
                            
                            + v * tauIC * BigR**4 * (u_xy * (p0_xx - p0_yy) - p0_xy * (u_xx - u_yy))  * xjac * theta * tstep &
                            
                            - BigR**3 * (particle_source(ms,mt)+source_pellet) * (v_x * u_x + v_y * u_y) * xjac * theta * tstep &
                            
                            + TG_num2 * 0.25d0 * r0_hat * BigR**3 * (w0_x * u_y - w0_y * u_x)       &
                                      * ( v_x * u0_y - v_y * u0_x) * xjac * theta * tstep * tstep   &
                                      
                            + TG_num2 * 0.25d0 * r0_hat * BigR**3 * (w0_x * u0_y - w0_y * u0_x)     &
                                      * ( v_x * u_y - v_y * u_x)   * xjac * theta * tstep * tstep
                  

                  !---------------------------------------- NEO
                  if ( NEO ) then
                    amat_22 =  amat_22 &
                         -amu_neo_prof(ms,mt)*BB2/((Btheta2+epsil)**2)*(ps0_x*v_x+ps0_y*v_y)*r0*(ps0_x*u_x+ps0_y*u_y)&
                         *BigR*xjac * theta * tstep 
                  endif
                  !---------------------------------------- NEO

                  amat_23 = - v * (ps0_s * zj_t  - ps0_t * zj_s)                * theta * tstep

                  amat_23_n = + eps_cyl * F0 / BigR * v * zj_p  * xjac          * theta * tstep

                  amat_24 = r0_hat * BigR**2 * w  * ( v_s * u0_t - v_t * u0_s)  * theta * tstep  &
                          + BigR * ( v_x * w_x + v_y * w_y) * visco_T  * xjac   * theta * tstep  &
                          + v * tauIC * BigR**4 * (p0_s * w_t - p0_t * w_s)     * theta * tstep  & 

                          + visco_num_T * (v_xx + v_x/BigR + v_yy)*(w_xx + w_x/BigR + w_yy) * xjac * theta * tstep    &

                          + TG_num2 * 0.25d0 * r0_hat * BigR**3 * (w_x * u0_y - w_y * u0_x)     &
                                    * ( v_x * u0_y - v_y * u0_x) * xjac * theta * tstep * tstep

                  amat_25 = + 0.5d0 * vv2 * (v_x * rho_y_hat - v_y * rho_x_hat)   * xjac * theta * tstep &
                            + rho_hat * BigR**2 * w0 * (v_s * u0_t - v_t * u0_s)         * theta * tstep &
                            - BigR**2 * (v_s * rho_t * T0   - v_t * rho_s * T0  )        * theta * tstep &
                            - BigR**2 * (v_s * rho   * T0_t - v_t * rho   * T0_s)        * theta * tstep &

                            + v * tauIC * BigR**4 * T0  * (rho_s * w0_t - rho_t * w0_s)  * theta * tstep &
                            + v * tauIC * BigR**4 * rho * (T0_s  * w0_t - T0_t  * w0_s)  * theta * tstep &
                            + tauIC * BigR**3 * (T0_y * rho + T0 * rho_y) * (v_x* u0_x + v_y * u0_y) * xjac * theta * tstep      &
                            + v * tauIC * BigR**4 * ( (u0_xy * (rho_xx*T0 + 2.d0*rho_x*T0_x + rho*T0_xx                          &
                                  -  rho_yy*T0 - 2.d0*rho_y*T0_y - rho*T0_yy))                           &
                                  - (rho_xy * T0 + rho_x*T0_y + rho_y*T0_x + rho*T0_xy) * (u0_xx - u0_yy)  )   &
                                                      * xjac * theta * tstep                             &
                            ! --- Diamagnetic viscosity
                            - dvisco_dT * bigR * W_dia_rho * (v_x*T0_x + v_y*T0_y)    * xjac * theta * tstep  &
                            - visco_T   * bigR * W_dia_rho * (v_xx + v_x/bigR + v_yy) * xjac * theta * tstep  &
                                
                            + TG_num2 * 0.25d0 * rho_hat * BigR**3 * (w0_x * u0_y - w0_y * u0_x)         &
                                      * ( v_x * u0_y - v_y * u0_x) * xjac * theta * tstep * tstep

                  !---------------------------------------- NEO
                  if ( NEO ) then
                    amat_25 = amat_25 &
                         -amu_neo_prof(ms,mt)*BB2/((Btheta2+epsil)**2)*(ps0_x*v_x+ps0_y*v_y)*&
                         (rho*(ps0_x*u0_x+ps0_y*u0_y)+&
                         tauIC*(ps0_x*(rho_x*T0+rho*T0_x)+ps0_y*(rho_y*T0+rho*T0_y)) &
                         +aki_neo_prof(ms,mt)*tauIC*rho*(ps0_x*T0_x+ps0_y*T0_y)&
                         -rho*Vpar0*Btheta2)*BigR*xjac*tstep*theta
                  endif
                  !---------------------------------------- NEO


                  amat_26 = - BigR**2 * (v_s * r0_t * T   - v_t * r0_s * T)      * theta * tstep  &
                            - BigR**2 * (v_s * r0   * T_t - v_t * r0   * T_s)    * theta * tstep  &
                            + dvisco_dT * T * ( v_x * w0_x + v_y * w0_y ) * BigR * xjac * theta * tstep &

                            + v * tauIC * BigR**4 * r0 * (T_s * w0_t - T_t * w0_s)    * theta * tstep &
                            + v * tauIC * BigR**4 * T  * (r0_s * w0_t - r0_t * w0_s)  * theta * tstep &
                            + tauIC * BigR**3 * (r0_y * T + r0 * T_y) * (v_x* u0_x + v_y * u0_y) * xjac * theta * tstep &
                            + v * tauIC * BigR**4 * ( (u0_xy * (T_xx*r0 + 2.d0*T_x*r0_x + T*r0_xx                       &
                                                              - T_yy*r0 - 2.d0*T_y*r0_y - T*r0_yy))                     &
                                                    - (T_xy * r0 + T_x*r0_y + T_y*r0_x + T*r0_xy) * (u0_xx - u0_yy)  )  &
                                                  * xjac * theta * tstep                                                &
                            ! --- Diamagnetic viscosity
                            - d2visco_dT2*T * bigR * W_dia   * (v_x*T0_x + v_y*T0_y)    * xjac * theta * tstep  &
                            - dvisco_dT*T   * bigR * W_dia   * (v_xx + v_x/bigR + v_yy) * xjac * theta * tstep  &
                            - dvisco_dT     * bigR * W_dia_T * (v_x*T0_x + v_y*T0_y)    * xjac * theta * tstep  &
                            - visco_T       * bigR * W_dia_T * (v_xx + v_x/bigR + v_yy) * xjac * theta * tstep  &
                            - dvisco_dT     * bigR * W_dia   * (v_x*T_x  + v_y*T_y )    * xjac * theta * tstep  
                            

                  !---------------------------------------- NEO
                  if ( NEO ) then
                    amat_26 = amat_26 -amu_neo_prof(ms,mt)*BB2/((Btheta2+epsil)**2)*(ps0_x*v_x+ps0_y*v_y)*&
                         (tauIC*(ps0_x*(r0_x*T+r0*T_x)+ps0_y*(r0_y*T+r0*T_y)) &
                         +aki_neo_prof(ms,mt)*tauIC*r0*(ps0_x*T_x+ps0_y*T_y))*BigR*xjac*tstep*theta
                    
                    amat_27= amu_neo_prof(ms,mt)*BB2/(Btheta2+epsil)*r0*vpar*(ps0_x*v_x+ps0_y*v_y)&
                         *BigR*xjac*tstep*theta 
                  else
                    amat_27 = 0
                  endif
                  !---------------------------------------- NEO

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


                  amat_51 = - (D_par-D_prof) * BigR * BB2_psi/ BB2**2 * Bgrad_rho_star     * Bgrad_rho     * xjac * theta * tstep &
                            + (D_par-D_prof) * BigR / BB2             * Bgrad_rho_star_psi * Bgrad_rho     * xjac * theta * tstep &
                            + (D_par-D_prof) * BigR / BB2             * Bgrad_rho_star     * Bgrad_rho_psi * xjac * theta * tstep &
                            + v * Vpar0 * (r0_s * psi_t - r0_t * psi_s)                                           * theta * tstep &
                            + v * r0 * (vpar0_s * psi_t - vpar0_t * psi_s)                                        * theta * tstep &
 
                            + TG_num5 * 0.25d0 / BigR * vpar0**2                                                              &
                                      * (r0_x * psi_y - r0_y * psi_x)                                                         &
                                      * ( v_x * ps0_y -  v_y * ps0_x                   ) * xjac * theta * tstep * tstep       &
                            + TG_num5 * 0.25d0 / BigR * vpar0**2                                                              &
                                      * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                                      &
                                      * ( v_x * psi_y -  v_y * psi_x                   ) * xjac * theta * tstep * tstep


                  amat_51_k = - (D_par-D_prof) * BigR * BB2_psi/ BB2**2 * Bgrad_rho_k_star * Bgrad_rho     * xjac * theta * tstep &
                              + (D_par-D_prof) * BigR / BB2             * Bgrad_rho_k_star * Bgrad_rho_psi * xjac * theta * tstep &

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

                          - v * 2.d0 * tauIC * (rho_y * T0 + rho*T0_y) * BigR                           * xjac * theta * tstep &

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

                  amat_56   = - v * 2.d0 * tauIC * (T_y * r0 + T*r0_y) * BigR                           * xjac * theta * tstep 

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

                  !###################################################################################################
                  !#  equation 6   (energy equation)                                                                 #
                  !###################################################################################################
                  Bgrad_T_star_psi = ( v_x  * psi_y - v_y  * psi_x  ) / BigR
                  Bgrad_T_psi      = ( T0_x * psi_y - T0_y * psi_x )  / BigR
                  Bgrad_T_T        = ( T_x * ps0_y - T_y * ps0_x ) / BigR
                  Bgrad_T_T_n      = ( F0 / BigR * T_p) / BigR

                  T_ps0_x = T_xx * ps0_y - T_xy * ps0_x + T_x * ps0_xy - T_y * ps0_xx
                  T_ps0_y = T_xy * ps0_y - T_yy * ps0_x + T_x * ps0_yy - T_y * ps0_xy
  
                  T0_psi_x = T0_xx * psi_y - T0_xy * psi_x + T0_x * psi_xy - T0_y * psi_xx
                  T0_psi_y = T0_xy * psi_y - T0_yy * psi_x + T0_x * psi_yy - T0_y * psi_xy
                  
                  v_psi_x = v_xx * psi_y - v_xy * psi_x + v_x * psi_xy - v_y * psi_xx
                  v_psi_y = v_xy * psi_y - v_yy * psi_x + v_x * psi_yy - v_y * psi_xy
                  
                  amat_61 = - (ZKpar_T-ZK_prof) * BigR * BB2_psi / BB2**2 * Bgrad_T_star     * Bgrad_T     * xjac * theta * tstep &
                            + (ZKpar_T-ZK_prof) * BigR / BB2              * Bgrad_T_star_psi * Bgrad_T     * xjac * theta * tstep &
                            + (ZKpar_T-ZK_prof) * BigR / BB2              * Bgrad_T_star     * Bgrad_T_psi * xjac * theta * tstep &

                            + v * r0 * Vpar0 * (T0_s * psi_t - T0_t * psi_s)                                     * theta * tstep &
                            + v * T0 * Vpar0 * (r0_s * psi_t - r0_t * psi_s)                                     * theta * tstep &
                            + v * r0 * GAMMA * T0 * (vpar0_s * psi_t - vpar0_t * psi_s)                          * theta * tstep &

                            + ZK_par_num * (v_psi_x  * ps0_y - v_psi_y  * ps0_x + v_ps0_x * psi_y - v_ps0_y * psi_x)          &
                                         * (T0_ps0_x * ps0_y - T0_ps0_y * ps0_x)                       * xjac * theta * tstep &
                            + ZK_par_num * (T0_psi_x * ps0_y - T0_psi_y * ps0_x + T0_ps0_x * psi_y - T0_ps0_y * psi_x)        &
                                         * (v_ps0_x  * ps0_y - v_ps0_y  * ps0_x)                       * xjac * theta * tstep &

                         + TG_num6 * 0.25d0 / BigR * vpar0**2                                                       &
                                   * T0 * (r0_x * psi_y - r0_y * psi_x)                                             &
                                   * ( v_x * ps0_y -  v_y * ps0_x                  ) * xjac * theta * tstep * tstep &
                         + TG_num6 * 0.25d0 / BigR * vpar0**2                                                       &
                                   * r0 * (T0_x * psi_y - T0_y * psi_x)                                             &
                                   * ( v_x * ps0_y -  v_y * ps0_x                  ) * xjac * theta * tstep * tstep &
                         + TG_num6 * 0.25d0 / BigR * vpar0**2                                                       &
                                   * T0 * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                          &
                                   * ( v_x * psi_y -  v_y * psi_x ) * xjac * theta * tstep * tstep                  &
                         + TG_num6 * 0.25d0 / BigR * vpar0**2                                                       &
                                   * r0 * (T0_x * ps0_y - T0_y * ps0_x + F0 / BigR * T0_p)                          &
                                   * ( v_x * psi_y -  v_y * psi_x ) * xjac * theta * tstep * tstep

                  amat_61_k = - (ZKpar_T-ZK_prof) * BigR * BB2_psi / BB2**2 * Bgrad_T_k_star * Bgrad_T     * xjac * theta * tstep &
                              + (ZKpar_T-ZK_prof) * BigR / BB2              * Bgrad_T_k_star * Bgrad_T_psi * xjac * theta * tstep &
  
                        + TG_num6 * 0.25d0 / BigR * vpar0**2                                                       &
                                  * T0 * (r0_x * psi_y - r0_y * psi_x)                                             &
                                  * (                            + F0 / BigR * v_p) * xjac * theta * tstep * tstep &
                        + TG_num6 * 0.25d0 / BigR * vpar0**2                                                       &
                                  * r0 * (T0_x * psi_y - T0_y * psi_x)                                             &
                                  * (                            + F0 / BigR * v_p) * xjac * theta * tstep * tstep


                  amat_62 = - v * r0 * BigR**2 * ( T0_x * u_y - T0_y * u_x)             * xjac * theta * tstep &
                            - v * T0 * BigR**2 * ( r0_x * u_y - r0_y * u_x)             * xjac * theta * tstep &
                            - v * r0 * 2.d0* GAMMA * BigR * T0 * u_y                    * xjac * theta * tstep &

                         + TG_num6 * 0.25d0 * BigR**2 * T0* (r0_x * u_y - r0_y * u_x)                &
                                            * ( v_x * u0_y - v_y * u0_x) * xjac * theta*tstep*tstep  &
                         + TG_num6 * 0.25d0 * BigR**2 * r0* (T0_x * u_y - T0_y * u_x)                &
                                            * ( v_x * u0_y - v_y * u0_x) * xjac * theta*tstep*tstep  &
                         + TG_num6 * 0.25d0 * BigR**2 * T0* (r0_x * u0_y - r0_y * u0_x)              &
                                            * ( v_x * u_y - v_y * u_x) * xjac * theta*tstep*tstep    &
                         + TG_num6 * 0.25d0 * BigR**2 * r0* (T0_x * u0_y - T0_y * u0_x)              &
                                            * ( v_x * u_y - v_y * u_x) * xjac * theta*tstep*tstep 

                  amat_63 = - v * (gamma-1.d0) * eta_T_ohm * 2.d0 * zj * zj0/(BigR**2.d0) * BigR * xjac * theta * tstep

                  amat_65 =   v * rho * T0   * BigR * xjac * (1.d0 + zeta)     &
                            - v * rho * BigR**2 * ( T0_s * u0_t - T0_t * u0_s)                          * theta * tstep &
                            - v * T0  * BigR**2 * ( rho_s * u0_t - rho_t * u0_s)                        * theta * tstep &
                            - v * rho * 2.d0* GAMMA * BigR * T0 * u0_y                           * xjac * theta * tstep &
                            + v * rho * F0 / BigR * Vpar0 * T0_p                                 * xjac * theta * tstep &

                            + v * rho * Vpar0 * (T0_s  * ps0_t - T0_t * ps0_s)                          * theta * tstep &
                            + v * T0  * Vpar0 * (rho_s * ps0_t - rho_t * ps0_s)                         * theta * tstep & 

                            + v * rho * GAMMA * T0 * (vpar0_s * ps0_t - vpar0_t * ps0_s)        * theta * tstep &
                            + v * rho * GAMMA * T0 * F0 / BigR * vpar0_p                 * xjac * theta * tstep &

                         + TG_num6 * 0.25d0 * BigR**2 * T0* (rho_x * u0_y - rho_y * u0_x)      &
                                   * ( v_x * u0_y - v_y * u0_x) * xjac * theta*tstep*tstep     &
                         + TG_num6 * 0.25d0 * BigR**2 * rho * (T0_x * u0_y - T0_y * u0_x)      &
                                   * ( v_x * u0_y - v_y * u0_x) * xjac* theta*tstep*tstep      &
                         + TG_num6 * 0.25d0 / BigR * vpar0**2 &
                                   * T0 * (rho_x * ps0_y - rho_y * ps0_x )                      &
                                   * ( v_x * ps0_y -  v_y * ps0_x ) * xjac * theta * tstep * tstep&
                         + TG_num6 * 0.25d0 / BigR * vpar0**2 &
                                   * rho * (T0_x * ps0_y - T0_y * ps0_x + F0 / BigR * T0_p)                        &
                                   * ( v_x * ps0_y -  v_y * ps0_x ) * xjac * theta * tstep * tstep

                  amat_65_n = + v * T0  * F0 / BigR * Vpar0 * rho_p                      * xjac * theta * tstep    &

                         + TG_num6 * 0.25d0 / BigR * vpar0**2 &
                                   * T0 * (                              + F0 / BigR * rho_p)                      &
                                   * ( v_x * ps0_y -  v_y * ps0_x                  ) * xjac * theta * tstep * tstep

                  amat_65_k = + TG_num6 * 0.25d0 / BigR * vpar0**2 &
                                   * T0 * (rho_x * ps0_y - rho_y * ps0_x                    )                      &
                                   * (                            + F0 / BigR * v_p) * xjac * theta * tstep * tstep&
                              + TG_num6 * 0.25d0 / BigR * vpar0**2 &
                                   * rho * (T0_x * ps0_y - T0_y * ps0_x + F0 / BigR * T0_p)                        &
                                   * (                            + F0 / BigR * v_p) * xjac * theta * tstep * tstep

                  amat_65_kn = + TG_num6 * 0.25d0 / BigR * vpar0**2 &
                                   * T0 * (+ F0 / BigR * rho_p)                      &
                                   * (     + F0 / BigR * v_p) * xjac * theta * tstep * tstep


                  amat_66 =   v * r0_corr * T   * BigR * xjac * (1.d0 + zeta)     &
                            - v * r0 * BigR**2 * ( T_s  * u0_t - T_t  * u0_s)                        * theta * tstep &
                            - v * T  * BigR**2 * ( r0_s * u0_t - r0_t * u0_s)                        * theta * tstep &

                            - v * r0 * 2.d0* GAMMA * BigR * T * u0_y                 * xjac * theta * tstep &

                            + v * r0 * Vpar0 * (T_s  * ps0_t - T_t  * ps0_s)                         * theta * tstep &
                            + v * T  * Vpar0 * (r0_s * ps0_t - r0_t * ps0_s)                         * theta * tstep & 

                            + v * r0 * GAMMA * T * (vpar0_s * ps0_t - vpar0_t * ps0_s)      * theta * tstep &
                            + v * r0 * GAMMA * T * F0 / BigR * vpar0_p               * xjac * theta * tstep &

                            + v * T * F0 / BigR * Vpar0 * r0_p                              * xjac * theta * tstep &

                            + (ZKpar_T-ZK_prof) * BigR / BB2 * Bgrad_T_star * Bgrad_T_T * xjac * theta * tstep &
                            + ZK_prof * BigR * ( v_x*T_x + v_y*T_y )                    * xjac * theta * tstep &

                            + dZKpar_dT * T * BigR / BB2 * Bgrad_T_star * Bgrad_T           * xjac * theta * tstep &
  
                            + ZK_perp_num * (v_xx + v_x/BigR + v_yy)*(T_xx + T_x/BigR + T_yy) * BigR * xjac * theta * tstep &

                        + TG_num6 * 0.25d0 * BigR**2 * T* (r0_x * u0_y - r0_y * u0_x)         &
                                  * ( v_x * u0_y - v_y * u0_x) * xjac * theta * tstep * tstep &
                        + TG_num6 * 0.25d0 * BigR**2 * r0* (T_x * u0_y - T_y * u0_x)          &
                                  * ( v_x * u0_y - v_y * u0_x) * xjac * theta * tstep * tstep &
                        + TG_num6 * 0.25d0 / BigR * vpar0**2 &
                                  * T * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                          &
                                  * ( v_x * ps0_y -  v_y * ps0_x                  ) * xjac * theta * tstep * tstep&
                        + TG_num6 * 0.25d0 / BigR * vpar0**2 &
                                  * r0 * (T_x * ps0_y - T_y * ps0_x               )                            &
                                  * ( v_x * ps0_y -  v_y * ps0_x                  ) * xjac * theta * tstep * tstep &

                        -v * T * (gamma-1.d0) * deta_dT_ohm * (zj0 / BigR)**2.d0        * BigR * xjac * theta * tstep



                  amat_66_k = + (ZKpar_T-ZK_prof) * BigR / BB2 * Bgrad_T_k_star * Bgrad_T_T    * xjac * theta * tstep &

                              + dZKpar_dT * T * BigR / BB2 * Bgrad_T_k_star * Bgrad_T          * xjac * theta * tstep &

                        + TG_num6 * 0.25d0 / BigR * vpar0**2 &
                                  * T * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                              &
                                  * (                                + F0 / BigR * v_p) * xjac * theta * tstep * tstep&
                        + TG_num6 * 0.25d0 / BigR * vpar0**2 &
                                  * r0 * (T_x * ps0_y - T_y * ps0_x                  )                                &
                                  * (                               + F0 / BigR * v_p) * xjac * theta * tstep * tstep

                  amat_66_n = + (ZKpar_T-ZK_prof) * BigR / BB2 * Bgrad_T_star   * Bgrad_T_T_n  * xjac * theta * tstep &

                              + v * r0 * F0 / BigR * Vpar0 * T_p                               * xjac * theta * tstep &
  
                        + TG_num6 * 0.25d0 / BigR * vpar0**2 &
                                  * r0 * ( + F0 / BigR * T_p)                            &
                                  * ( v_x * ps0_y -  v_y * ps0_x ) * xjac * theta * tstep * tstep

                  amat_66_kn = + (ZKpar_T-ZK_prof) * BigR / BB2 * Bgrad_T_k_star * Bgrad_T_T_n * xjac * theta * tstep &
                               + ZK_prof * BigR   * (v_p*T_p /BigR**2 )                       * xjac * theta * tstep  &

                        + TG_num6 * 0.25d0 / BigR * vpar0**2 &
                                  * r0 * ( + F0 / BigR * T_p)                            &
                                  * (      + F0 / BigR * v_p) * xjac * theta * tstep * tstep

                  amat_67 = + v * r0 * F0 / BigR * Vpar * T0_p                              * xjac * theta * tstep &
                            + v * T0 * F0 / BigR * Vpar * r0_p                              * xjac * theta * tstep &

                            + v * r0 * Vpar * (T0_s * ps0_t - T0_t * ps0_s)                        * theta * tstep &
                            + v * T0 * Vpar * (r0_s * ps0_t - r0_t * ps0_s)                        * theta * tstep & 

                            + v * r0 * GAMMA * T0 * (vpar_s * ps0_t - vpar_t * ps0_s)       * theta * tstep        &
  
                        + TG_num6 * 0.25d0 / BigR * 2.d0 * vpar0*vpar &
                                  * T0 * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                          &
                                  * ( v_x * ps0_y -  v_y * ps0_x                        ) * xjac * theta * tstep * tstep &
                        + TG_num6 * 0.25d0 / BigR * 2.d0 * vpar0*vpar &
                                  * r0 * (T0_x * ps0_y - T0_y * ps0_x + F0 / BigR * T0_p)                          &
                                  * ( v_x * ps0_y -  v_y * ps0_x                        ) * xjac * theta * tstep * tstep 

                  amat_67_k =  &
                        + TG_num6 * 0.25d0 / BigR * 2.d0 * vpar0*vpar &
                                  * T0 * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                          &
                                  * (                            + F0 / BigR * v_p) * xjac * theta * tstep * tstep &
                        + TG_num6 * 0.25d0 / BigR * 2.d0 * vpar0*vpar &
                                  * r0 * (T0_x * ps0_y - T0_y * ps0_x + F0 / BigR * T0_p)                          &
                                  * (                            + F0 / BigR * v_p) * xjac * theta * tstep * tstep 


                  amat_67_n = + v * r0 * GAMMA * T0 * F0 / BigR * vpar_p             * xjac * theta * tstep


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

                            + v*(particle_source(ms,mt) + source_pellet)*vpar0* BB2_psi * BigR * xjac * theta * tstep &
  
                            + TG_NUM7 * 0.25d0 * r0 * Vpar0**2 * BB2 &
                                      * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac) / BigR  &
                                      * (-(psi_s * v_t     - psi_t * v_s)    /xjac)  * xjac * theta * tstep*tstep &

                            + TG_NUM7 * 0.25d0 * r0 * Vpar0**2 * BB2 &
                                      * (-(psi_s * vpar0_t - psi_t * vpar0_s)/xjac) / BigR  &
                                      * (-(ps0_s * v_t     - ps0_t * v_s)    /xjac)  * xjac * theta * tstep*tstep &

                            + TG_NUM7 * 0.25d0 * v  * Vpar0**2 * BB2 &
                                      * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac) / BigR  &
                                      * (-(psi_s * r0_t    - psi_t * r0_s)   /xjac)  * xjac * theta * tstep*tstep &

                            + TG_NUM7 * 0.25d0 * v  * Vpar0**2 * BB2 &
                                      * (-(psi_s * vpar0_t - psi_t * vpar0_s)/xjac) / BigR  &
                                      * (-(ps0_s * r0_t    - ps0_t * r0_s)   /xjac)  * xjac * theta * tstep*tstep 

                  if (normalized_velocity_profile) then
                    amat_71 = amat_71  - visco_par * (v_x * Vt_x_psi   + v_y * Vt_y_psi) * BigR * xjac * theta * tstep      
                  else
                    amat_71 = amat_71  - visco_par * 2.d0 * PI * F0 * (v_x * Omega_tor_x_psi + v_y * Omega_tor_y_psi) * BigR * xjac * theta * tstep 
                  endif

                  amat_71_k = - 0.5d0 * r0 * vpar0**2 * BB2_psi * F0 / BigR * v_p    * xjac * theta * tstep 

                  amat_72 = 0.d0 

                  !---------------------------------------- NEO
                  if ( NEO ) then
                    amat_71 = amat_71 &
                         -v*amu_neo_prof(ms,mt)*BB2/(Btheta2+epsil)*(r0*(psi_x*u0_x+psi_y*u0_y)+&
                         tauIC*(psi_x*P0_x+psi_y*P0_y) &
                         +aki_neo_prof(ms,mt)*tauIC*r0*(psi_x*T0_x+psi_y*T0_y))*BigR*xjac*theta*tstep &

                         -v*amu_neo_prof(ms,mt)*(-Btheta2_psi)*BB2/(Btheta2**2)*&
                         (r0*(ps0_x*u0_x+ps0_y*u0_y)+&
                         tauIC*(ps0_x*P0_x+ps0_y*P0_y) &                     
                         +aki_neo_prof(ms,mt)*tauIC*r0*(ps0_x*T0_x+ps0_y*T0_y))*BigR*xjac*theta*tstep

                    amat_72 =  amat_72 &
                         -v*amu_neo_prof(ms,mt)*BB2/(Btheta2+epsil)*r0*(ps0_x*u_x+ps0_y*u_y)&
                         *BigR*xjac * theta * tstep 
                  endif
                  !---------------------------------------- NEO

                  amat_75 = + v * (rho_s * T0 * ps0_t - rho_t * T0 * ps0_s)                 * theta * tstep &
                            + v * (rho * T0_s * ps0_t - rho * T0_t * ps0_s)                 * theta * tstep &
                            + v * F0 / BigR * rho * T0_p                             * xjac * theta * tstep &

                            + 0.5d0 * rho * vpar0**2 * BB2 * (ps0_s * v_t - ps0_t * v_s)    * theta * tstep &
                            + 0.5d0 * v   * vpar0**2 * BB2 * (ps0_s * rho_t - ps0_t * rho_s)* theta * tstep &

                            + TG_NUM7 * 0.25d0 * rho * Vpar0**2 * BB2 &
                                      * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac + F0 / BigR * vpar0_p) / BigR  &
                                      * (-(ps0_s * v_t     - ps0_t * v_s)    /xjac          )  * xjac * theta * tstep*tstep &

                            + TG_NUM7 * 0.25d0 * v * Vpar0**2 * BB2 &
                                      * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac + F0 / BigR * vpar0_p) / BigR  &
                                      * (-(ps0_s * rho_t   - ps0_t * rho_s)  /xjac           ) * xjac * theta * tstep*tstep

                  amat_75_k = - 0.5d0 * rho * vpar0**2 * BB2 * F0 / BigR * v_p       * xjac * theta * tstep &

                              + TG_NUM7 * 0.25d0 * rho * Vpar0**2 * BB2 &
                                        * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac + F0 / BigR * vpar0_p) / BigR  &
                                        * (                                          + F0 / BigR * v_p)  * xjac * theta * tstep*tstep

                  amat_75_n = + v * F0 / BigR * rho_p * T0                           * xjac * theta * tstep &
                              - 0.5d0 * v   * vpar0**2 * BB2 * F0 / BigR * rho_p       * xjac * theta * tstep &

                              + TG_NUM7 * 0.25d0 * v * Vpar0**2 * BB2 &
                                        * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac + F0 / BigR * vpar0_p) / BigR  &
                                        * (                                          + F0 / BigR * rho_p)* xjac * theta * tstep*tstep
                  
                  amat_76 = + v * (T_s * R0 * ps0_t - T_t * R0 * ps0_s)                     * theta * tstep &
                            + v * (T * R0_s * ps0_t - T * R0_t * ps0_s)                     * theta * tstep &
                            + v * F0 / BigR * T * R0_p                               * xjac * theta * tstep

                  amat_76_n= + v * F0 / BigR * T_p * R0                              * xjac * theta * tstep

                  amat_77 = v * Vpar * r0_corr * F0**2 / BigR * xjac * (1.d0 + zeta) &

                          + v*(particle_source(ms,mt) + source_pellet)*vpar*BB2 * BigR * xjac * theta * tstep &

                          + r0 * vpar0 * vpar * BB2 * (ps0_s * v_t - ps0_t * v_s)             * theta * tstep &
                          + v  * vpar0 * vpar * BB2 * (ps0_s * r0_t - ps0_t * r0_s)           * theta * tstep &
                          - v  * vpar0 * vpar * BB2 * F0 / BigR * r0_p                 * xjac * theta * tstep &

                          + visco_par_num * (v_xx + v_x/BigR + v_yy)*(vpar_xx + vpar_x/BigR + vpar_yy) * BigR * xjac * theta * tstep&

                          + TG_NUM7 * 0.5d0 * r0 * Vpar * Vpar0 * BB2 &
                                    * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac + F0 / BigR * vpar0_p) / BigR                      &
                                    * (-(ps0_s * v_t     - ps0_t * v_s)    /xjac                      )  * xjac * theta * tstep*tstep  &
                          + TG_NUM7 * 0.5d0 * v * Vpar * Vpar0 * BB2 &
                                    * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac + F0 / BigR * vpar0_p) / BigR                      &
                                    * (-(ps0_s * r0_t    - ps0_t * r0_s)   /xjac + F0 / BigR * r0_p)  * xjac * theta * tstep*tstep  &
                          + TG_NUM7 * 0.25d0 * r0 * Vpar0**2 * BB2 &
                                    * (-(ps0_s * vpar_t - ps0_t * vpar_s)/xjac                   ) / BigR                           &
                                    * (-(ps0_s * v_t    - ps0_t * v_s)   /xjac                   )  * xjac * theta * tstep*tstep    &
                          + TG_NUM7 * 0.25d0 * v * Vpar0**2 * BB2 &
                                    * (-(ps0_s * vpar_t - ps0_t * vpar_s)/xjac                   ) / BigR                           &
                                    * (-(ps0_s * r0_t   - ps0_t * r0_s)  /xjac + F0 / BigR * r0_p)  * xjac * theta * tstep*tstep 

                  if (normalized_velocity_profile) then
                    amat_77 = amat_77 + visco_par * (v_x * Vpar_x + v_y * Vpar_y) * BigR        * xjac  * theta * tstep 
                  else
                    amat_77 = amat_77 + visco_par * F0**2 / BigR**2 * (v_x * (Vpar_x - 2*vpar/BigR) + v_y * Vpar_y) * BigR * xjac  * theta * tstep 
                  endif

                  amat_77_k = - r0 * vpar0 * vpar * BB2 * F0 / BigR * v_p               * xjac * theta * tstep             &
     
                            + TG_NUM7 * 0.5d0 * r0 * Vpar * Vpar0 * BB2 &
                                      * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac + F0 / BigR * vpar0_p) / BigR                     &
                                      * (                                          + F0 / BigR * v_p)  * xjac * theta * tstep*tstep  &
                            + TG_NUM7 * 0.25d0 * r0 * Vpar0**2 * BB2 &
                                      * (-(ps0_s * vpar_t - ps0_t * vpar_s)/xjac                  ) / BigR                           &
                                      * (                                        + F0 / BigR * v_p)  * xjac * theta * tstep*tstep

                  amat_77_n = &

                            + TG_NUM7 * 0.25d0 * r0 * Vpar0**2 * BB2 &
                                      * (                                        + F0 / BigR * vpar_p) / BigR                        &
                                      * (-(ps0_s * v_t    - ps0_t * v_s)   /xjac                     )  * xjac * theta * tstep*tstep &
                            + TG_NUM7 * 0.25d0 * v * Vpar0**2 * BB2 &
                                      * (                                        + F0 / BigR * vpar_p) / BigR                        &
                                      * (-(ps0_s * r0_t   - ps0_t * r0_s)  /xjac + F0 / BigR * r0_p)  * xjac * theta * tstep*tstep 

                  amat_77_kn = &

                             + TG_NUM7 * 0.25d0 * r0 * Vpar0**2 * BB2 &
                                       * (                                        + F0 / BigR * vpar_p) / BigR                        &
                                       * (                                        + F0 / BigR * v_p)  * xjac * theta * tstep*tstep

                  !---------------------------------------- NEO
                  if ( NEO ) then
                    amat_75 = amat_75 &
                         -v*amu_neo_prof(ms,mt)*BB2/(Btheta2+epsil)* &
                         (rho*(ps0_x*u0_x+ps0_y*u0_y)+&
                         tauIC*(ps0_x*(rho_x*T0+rho*T0_x)+ps0_y*(rho_y*T0+rho*T0_y)) &
                         +aki_neo_prof(ms,mt)*tauIC*rho*(ps0_x*T0_x+ps0_y*T0_y)&
                         -rho*Vpar0*Btheta2)*BigR*xjac*tstep*theta
                    
                    amat_76 = amat_76 -v*amu_neo_prof(ms,mt)*BB2/(Btheta2+epsil)* &
                         (tauIC*(ps0_x*(r0_x*T+r0*T_x)+ps0_y*(r0_y*T+r0*T_y)) &
                         +aki_neo_prof(ms,mt)*tauIC*r0*(ps0_x*T_x+ps0_y*T_y))*BigR*xjac*tstep*theta
                  
                    amat_77= amat_77+ v*amu_neo_prof(ms,mt)*BB2*r0*vpar &
                         *BigR*xjac*tstep*theta 
                  endif
                  !---------------------------------------- NEO


                  !###################################################################################################
                  !# end equation 7                                                                                  #
                  !###################################################################################################

                  if (use_fft) then

                    index_kl =       n_var*(n_order+1)*(k-1) +       n_var*(l-1) + 1
                  else
                    index_kl = n_tor_local*n_var*(n_order+1)*(k-1) + n_tor_local*n_var*(l-1) + in - n_tor_start +1
                  endif

                  ! --- Fill up the matrix
                  if (use_fft) then

                    kl1 = index_kl
                    kl2 = index_kl + 1
                    kl3 = index_kl + 2
                    kl4 = index_kl + 3
                    kl5 = index_kl + 4
                    kl6 = index_kl + 5
                    kl7 = index_kl + 6
                    ij1 = 1
                    ij2 = 2
                    ij3 = 3
                    ij4 = 4
                    ij5 = 5
                    ij6 = 6
                    ij7 = 7

                    ELM_p(mp,kl1,ij1)  =  ELM_p(mp,kl1,ij1) + wst * amat_11
                    ELM_n(mp,kl2,ij1)  =  ELM_n(mp,kl2,ij1) + wst * amat_12_n

                    ELM_p(mp,kl2,ij1)  =  ELM_p(mp,kl2,ij1) + wst * amat_12
                    ELM_p(mp,kl3,ij1)  =  ELM_p(mp,kl3,ij1) + wst * amat_13
                    ELM_p(mp,kl5,ij1)  =  ELM_p(mp,kl5,ij1) + wst * amat_15
                    ELM_n(mp,kl5,ij1)  =  ELM_n(mp,kl5,ij1) + wst * amat_15_n
                    ELM_p(mp,kl6,ij1)  =  ELM_p(mp,kl6,ij1) + wst * amat_16
                    ELM_n(mp,kl6,ij1)  =  ELM_n(mp,kl6,ij1) + wst * amat_16_n

                    ELM_p(mp,kl1,ij2)  =  ELM_p(mp,kl1,ij2) + wst * amat_21
                    ELM_p(mp,kl2,ij2)  =  ELM_p(mp,kl2,ij2) + wst * amat_22
                    ELM_p(mp,kl3,ij2)  =  ELM_p(mp,kl3,ij2) + wst * amat_23
                    ELM_n(mp,kl3,ij2)  =  ELM_n(mp,kl3,ij2) + wst * amat_23_n

                    ELM_p(mp,kl4,ij2)  =  ELM_p(mp,kl4,ij2) + wst * amat_24
                    ELM_p(mp,kl5,ij2)  =  ELM_p(mp,kl5,ij2) + wst * amat_25
                    ELM_p(mp,kl6,ij2)  =  ELM_p(mp,kl6,ij2) + wst * amat_26
                    !---------------------------------------- NEO
                    ELM_p(mp,kl7,ij2)  =  ELM_p(mp,kl7,ij2) + wst * amat_27
                    !---------------------------------------- NEO


                    ELM_p(mp,kl1,ij3)  =  ELM_p(mp,kl1,ij3) + wst * amat_31
                    ELM_p(mp,kl3,ij3)  =  ELM_p(mp,kl3,ij3) + wst * amat_33

                    ELM_p(mp,kl2,ij4)  =  ELM_p(mp,kl2,ij4) + wst * amat_42
                    ELM_p(mp,kl4,ij4)  =  ELM_p(mp,kl4,ij4) + wst * amat_44


                    ELM_p(mp,kl1,ij5)  =  ELM_p(mp,kl1,ij5)  + wst * amat_51
                    ELM_k(mp,kl1,ij5)  =  ELM_k(mp,kl1,ij5)  + wst * amat_51_k

                    ELM_p(mp,kl2,ij5)  =  ELM_p(mp,kl2,ij5)  + wst * amat_52

                    ELM_p(mp,kl5,ij5)  =  ELM_p(mp,kl5,ij5)  + wst * amat_55
                    ELM_k(mp,kl5,ij5)  =  ELM_k(mp,kl5,ij5)  + wst * amat_55_k
                    ELM_n(mp,kl5,ij5)  =  ELM_n(mp,kl5,ij5)  + wst * amat_55_n
                    ELM_kn(mp,kl5,ij5) =  ELM_kn(mp,kl5,ij5) + wst * amat_55_kn

                    ELM_p(mp,kl6,ij5)  =  ELM_p(mp,kl6,ij5)  + wst * amat_56

                    ELM_p(mp,kl7,ij5)  =  ELM_p(mp,kl7,ij5)  + wst * amat_57
                    ELM_n(mp,kl7,ij5)  =  ELM_n(mp,kl7,ij5)  + wst * amat_57_n
                    ELM_k(mp,kl7,ij5)  =  ELM_k(mp,kl7,ij5)  + wst * amat_57_k
                    ELM_kn(mp,kl7,ij5) =  ELM_kn(mp,kl7,ij5) + wst * amat_57_kn

                    ELM_p(mp,kl1,ij6)  =  ELM_p(mp,kl1,ij6)  + wst * amat_61
                    ELM_k(mp,kl1,ij6)  =  ELM_k(mp,kl1,ij6)  + wst * amat_61_k

                    ELM_p(mp,kl2,ij6)  =  ELM_p(mp,kl2,ij6)  + wst * amat_62
                    ELM_p(mp,kl3,ij6)  =  ELM_p(mp,kl3,ij6)  + wst * amat_63
                    ELM_p(mp,kl5,ij6)  =  ELM_p(mp,kl5,ij6)  + wst * amat_65
                    ELM_n(mp,kl5,ij6)  =  ELM_n(mp,kl5,ij6)  + wst * amat_65_n
                    ELM_k(mp,kl5,ij6)  =  ELM_k(mp,kl5,ij6)  + wst * amat_65_k
                    ELM_kn(mp,kl5,ij6) =  ELM_kn(mp,kl5,ij6) + wst * amat_65_kn

                    ELM_p(mp,kl6,ij6)  =  ELM_p(mp,kl6,ij6)  + wst * amat_66
                    ELM_k(mp,kl6,ij6)  =  ELM_k(mp,kl6,ij6)  + wst * amat_66_k
                    ELM_n(mp,kl6,ij6)  =  ELM_n(mp,kl6,ij6)  + wst * amat_66_n
                    ELM_kn(mp,kl6,ij6) =  ELM_kn(mp,kl6,ij6) + wst * amat_66_kn

                    ELM_p(mp,kl7,ij6)  =  ELM_p(mp,kl7,ij6)  + wst * amat_67
                    ELM_k(mp,kl7,ij6)  =  ELM_k(mp,kl7,ij6)  + wst * amat_67_k
                    ELM_n(mp,kl7,ij6)  =  ELM_n(mp,kl7,ij6)  + wst * amat_67_n

                    ELM_p(mp,kl1,ij7)  =  ELM_p(mp,kl1,ij7)  + wst * amat_71
                    ELM_k(mp,kl1,ij7)  =  ELM_k(mp,kl1,ij7)  + wst * amat_71_k

                    ELM_p(mp,kl2,ij7)  =  ELM_p(mp,kl2,ij7)  + wst * amat_72

                    ELM_p(mp,kl5,ij7)  =  ELM_p(mp,kl5,ij7)  + wst * amat_75
                    ELM_n(mp,kl5,ij7)  =  ELM_n(mp,kl5,ij7)  + wst * amat_75_n
                    ELM_k(mp,kl5,ij7)  =  ELM_k(mp,kl5,ij7)  + wst * amat_75_k

                    ELM_p(mp,kl6,ij7)  =  ELM_p(mp,kl6,ij7)  + wst * amat_76
                    ELM_n(mp,kl6,ij7)  =  ELM_n(mp,kl6,ij7)  + wst * amat_76_n

                    ELM_p(mp,kl7,ij7)  =  ELM_p(mp,kl7,ij7)  + wst * amat_77
                    ELM_k(mp,kl7,ij7)  =  ELM_k(mp,kl7,ij7)  + wst * amat_77_k
                    ELM_n(mp,kl7,ij7)  =  ELM_n(mp,kl7,ij7)  + wst * amat_77_n
                    ELM_kn(mp,kl7,ij7) =  ELM_kn(mp,kl7,ij7) + wst * amat_77_kn
                  else
                    kl1 = index_kl
                    kl2 = index_kl + 1 * (n_tor_end - n_tor_start + 1)
                    kl3 = index_kl + 2 * (n_tor_end - n_tor_start + 1) 
                    kl4 = index_kl + 3 * (n_tor_end - n_tor_start + 1)
                    kl5 = index_kl + 4 * (n_tor_end - n_tor_start + 1)
                    kl6 = index_kl + 5 * (n_tor_end - n_tor_start + 1)
                    kl7 = index_kl + 6 * (n_tor_end - n_tor_start + 1)
                    ij1 = index_ij
                    ij2 = index_ij + 1 * (n_tor_end - n_tor_start + 1)
                    ij3 = index_ij + 2 * (n_tor_end - n_tor_start + 1)
                    ij4 = index_ij + 3 * (n_tor_end - n_tor_start + 1)
                    ij5 = index_ij + 4 * (n_tor_end - n_tor_start + 1)
                    ij6 = index_ij + 5 * (n_tor_end - n_tor_start + 1)
                    ij7 = index_ij + 6 * (n_tor_end - n_tor_start + 1)
                    
                    ELM(ij1,kl1) = ELM(ij1,kl1) + (amat_11 + amat_11_k + amat_11_n + amat_11_kn) * wst
                    ELM(ij1,kl2) = ELM(ij1,kl2) + (amat_12 + amat_12_k + amat_12_n + amat_12_kn) * wst
                    ELM(ij1,kl3) = ELM(ij1,kl3) + (amat_13 + amat_13_k + amat_13_n + amat_13_kn) * wst
                    ELM(ij1,kl4) = ELM(ij1,kl4) + (amat_14 + amat_14_k + amat_14_n + amat_14_kn) * wst
                    ELM(ij1,kl5) = ELM(ij1,kl5) + (amat_15 + amat_15_k + amat_15_n + amat_15_kn) * wst
                    ELM(ij1,kl6) = ELM(ij1,kl6) + (amat_16 + amat_16_k + amat_16_n + amat_16_kn) * wst
                    ELM(ij1,kl7) = ELM(ij1,kl7) + (amat_17 + amat_17_k + amat_17_n + amat_17_kn) * wst
                    
                    ELM(ij2,kl1) = ELM(ij2,kl1) + (amat_21 + amat_21_k + amat_21_n + amat_21_kn) * wst
                    ELM(ij2,kl2) = ELM(ij2,kl2) + (amat_22 + amat_22_k + amat_22_n + amat_22_kn) * wst
                    ELM(ij2,kl3) = ELM(ij2,kl3) + (amat_23 + amat_23_k + amat_23_n + amat_23_kn) * wst
                    ELM(ij2,kl4) = ELM(ij2,kl4) + (amat_24 + amat_24_k + amat_24_n + amat_24_kn) * wst
                    ELM(ij2,kl5) = ELM(ij2,kl5) + (amat_25 + amat_25_k + amat_25_n + amat_25_kn) * wst
                    ELM(ij2,kl6) = ELM(ij2,kl6) + (amat_26 + amat_26_k + amat_26_n + amat_26_kn) * wst
                    ELM(ij2,kl7) = ELM(ij2,kl7) + (amat_27 + amat_27_k + amat_27_n + amat_27_kn) * wst
                    
                    ELM(ij3,kl1) = ELM(ij3,kl1) + (amat_31 + amat_31_k + amat_31_n + amat_31_kn) * wst
                    ELM(ij3,kl2) = ELM(ij3,kl2) + (amat_32 + amat_32_k + amat_32_n + amat_32_kn) * wst
                    ELM(ij3,kl3) = ELM(ij3,kl3) + (amat_33 + amat_33_k + amat_33_n + amat_33_kn) * wst
                    ELM(ij3,kl4) = ELM(ij3,kl4) + (amat_34 + amat_34_k + amat_34_n + amat_34_kn) * wst
                    ELM(ij3,kl5) = ELM(ij3,kl5) + (amat_35 + amat_35_k + amat_35_n + amat_35_kn) * wst
                    ELM(ij3,kl6) = ELM(ij3,kl6) + (amat_36 + amat_36_k + amat_36_n + amat_36_kn) * wst
                    ELM(ij3,kl7) = ELM(ij3,kl7) + (amat_37 + amat_37_k + amat_37_n + amat_37_kn) * wst
                    
                    ELM(ij4,kl1) = ELM(ij4,kl1) + (amat_41 + amat_41_k + amat_41_n + amat_41_kn) * wst
                    ELM(ij4,kl2) = ELM(ij4,kl2) + (amat_42 + amat_42_k + amat_42_n + amat_42_kn) * wst
                    ELM(ij4,kl3) = ELM(ij4,kl3) + (amat_43 + amat_43_k + amat_43_n + amat_43_kn) * wst
                    ELM(ij4,kl4) = ELM(ij4,kl4) + (amat_44 + amat_44_k + amat_44_n + amat_44_kn) * wst
                    ELM(ij4,kl5) = ELM(ij4,kl5) + (amat_45 + amat_45_k + amat_45_n + amat_45_kn) * wst
                    ELM(ij4,kl6) = ELM(ij4,kl6) + (amat_46 + amat_46_k + amat_46_n + amat_46_kn) * wst
                    ELM(ij4,kl7) = ELM(ij4,kl7) + (amat_47 + amat_47_k + amat_47_n + amat_47_kn) * wst
                    
                    ELM(ij5,kl1) = ELM(ij5,kl1) + (amat_51 + amat_51_k + amat_51_n + amat_51_kn) * wst
                    ELM(ij5,kl2) = ELM(ij5,kl2) + (amat_52 + amat_52_k + amat_52_n + amat_52_kn) * wst
                    ELM(ij5,kl3) = ELM(ij5,kl3) + (amat_53 + amat_53_k + amat_53_n + amat_53_kn) * wst
                    ELM(ij5,kl4) = ELM(ij5,kl4) + (amat_54 + amat_54_k + amat_54_n + amat_54_kn) * wst
                    ELM(ij5,kl5) = ELM(ij5,kl5) + (amat_55 + amat_55_k + amat_55_n + amat_55_kn) * wst
                    ELM(ij5,kl6) = ELM(ij5,kl6) + (amat_56 + amat_56_k + amat_56_n + amat_56_kn) * wst
                    ELM(ij5,kl7) = ELM(ij5,kl7) + (amat_57 + amat_57_k + amat_57_n + amat_57_kn) * wst
                    
                    ELM(ij6,kl1) = ELM(ij6,kl1) + (amat_61 + amat_61_k + amat_61_n + amat_61_kn) * wst
                    ELM(ij6,kl2) = ELM(ij6,kl2) + (amat_62 + amat_62_k + amat_62_n + amat_62_kn) * wst
                    ELM(ij6,kl3) = ELM(ij6,kl3) + (amat_63 + amat_63_k + amat_63_n + amat_63_kn) * wst
                    ELM(ij6,kl4) = ELM(ij6,kl4) + (amat_64 + amat_64_k + amat_64_n + amat_64_kn) * wst
                    ELM(ij6,kl5) = ELM(ij6,kl5) + (amat_65 + amat_65_k + amat_65_n + amat_65_kn) * wst
                    ELM(ij6,kl6) = ELM(ij6,kl6) + (amat_66 + amat_66_k + amat_66_n + amat_66_kn) * wst
                    ELM(ij6,kl7) = ELM(ij6,kl7) + (amat_67 + amat_67_k + amat_67_n + amat_67_kn) * wst
                    
                    ELM(ij7,kl1) = ELM(ij7,kl1) + (amat_71 + amat_71_k + amat_71_n + amat_71_kn) * wst
                    ELM(ij7,kl2) = ELM(ij7,kl2) + (amat_72 + amat_72_k + amat_72_n + amat_72_kn) * wst
                    ELM(ij7,kl3) = ELM(ij7,kl3) + (amat_73 + amat_73_k + amat_73_n + amat_73_kn) * wst
                    ELM(ij7,kl4) = ELM(ij7,kl4) + (amat_74 + amat_74_k + amat_74_n + amat_74_kn) * wst
                    ELM(ij7,kl5) = ELM(ij7,kl5) + (amat_75 + amat_75_k + amat_75_n + amat_75_kn) * wst
                    ELM(ij7,kl6) = ELM(ij7,kl6) + (amat_76 + amat_76_k + amat_76_n + amat_76_kn) * wst
                    ELM(ij7,kl7) = ELM(ij7,kl7) + (amat_77 + amat_77_k + amat_77_n + amat_77_kn) * wst
                    
                  endif

                enddo ! in loop (n_tor, or not...)

              enddo ! l loop (n_order+1)
            enddo ! k loop (n_vertex)

          enddo ! im loop (n_tor, or not...)

        enddo ! mp loop (n_plane)

      enddo ! ms loop
    enddo ! mt loop


    if (use_fft) then

      do i_v = 1, n_var
        do j_loc=1, n_vertex_max*n_var*(n_order+1)
          !index_ij = n_var*(n_order+1)*(i-1) + n_var * (j-1) + 1
          !i_loc = index_ij + i_v-1
          i_loc = n_var*(n_order+1)*(i-1) + n_var * (j-1) + i_v 
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
        enddo

      enddo

    endif ! apply fft (or not)

  enddo ! j loop (n_order+1)
enddo ! i loop (n_vertex)

if (.NOT. use_fft) return

ELM = 0.5d0 * ELM

do j=1, n_vertex_max*n_var*(n_order+1)

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

do j=1, n_vertex_max*n_var*(n_order+1)

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
end module mod_elt_matrix_fft
