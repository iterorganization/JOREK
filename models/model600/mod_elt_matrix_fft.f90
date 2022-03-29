module mod_elt_matrix_fft

  implicit none

contains

#include "corr_neg_include.f90"

subroutine element_matrix_fft(element, nodes, xpoint2, xcase2, R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint, Z_xpoint, ELM, RHS, tid, &
                              ELM_p, ELM_n, ELM_k, ELM_kn, RHS_p, RHS_k,  eq_g, eq_s, eq_t, eq_p, eq_ss, eq_st, eq_tt, delta_g, delta_s, delta_t, &
                              i_tor_min, i_tor_max, aux_nodes, get_terms)
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
use diffusivities, only: get_dperp, get_zkperp, get_zk_iperp, get_zk_eperp
use equil_info, only : get_psi_n
use corr_neg
use mod_neutral_source
use mod_bootstrap_functions
use mod_atomic_coeff_deuterium, only : atomic_coeff_deuterium
use mod_impurity, only: radiation_function, radiation_function_linear
use mod_sources
use mod_model_settings

implicit none

type (type_element)       :: element
type (type_node)          :: nodes(n_vertex_max)
type (type_node),optional :: aux_nodes(n_vertex_max)

#define DIM0 n_tor*n_vertex_max*(n_order+1)*n_var

integer, intent(in)            :: tid
integer, intent(in)            :: i_tor_min, i_tor_max

real*8, dimension (DIM0,DIM0)  :: ELM
real*8, dimension (DIM0)       :: RHS

logical, intent(in), optional  :: get_terms 

integer    :: i, j, ms, mt, mp, k, l, index_ij, index_kl, index, index_k, index_m, m, ik, xcase2
integer    :: n_tor_start, n_tor_end, n_tor_local, n_tor_loop
integer    :: in, im, ij1, ij2, ij3, ij4, ij5, ij6, ij7, ij8, kl1, kl2, kl3, kl4, kl5, kl6, kl7, kl8, ij, kl
real*8     :: wst, xjac, xjac_s, xjac_t, xjac_x, xjac_y, BigR, r2, phi, delta_phi
real*8     :: current_source(n_gauss,n_gauss),particle_source(n_gauss,n_gauss),heat_source(n_gauss,n_gauss),heat_source_i(n_gauss,n_gauss),heat_source_e(n_gauss,n_gauss)
real*8     :: R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint(2), Z_xpoint(2), dj_dpsi, dj_dz, source_pellet, source_volume
real*8     :: Bgrad_rho_star,     Bgrad_rho,     Bgrad_T_star,  Bgrad_Ti, Bgrad_Te, Bgrad_T, BB2
real*8     :: Bgrad_rho_star_psi, Bgrad_rho_psi, Bgrad_rho_rho, Bgrad_T_star_psi, Bgrad_Ti_psi, Bgrad_T_psi, Bgrad_Ti_Ti, Bgrad_Te_psi, Bgrad_T_T, Bgrad_Te_Te, BB2_psi
real*8     :: Bgrad_rho_k_star, Bgrad_T_k_star, Bgrad_Ti_Ti_n, Bgrad_Te_Te_n, Bgrad_T_T_n, Bgrad_rho_rho_n
real*8     :: ZK_par_T, dZK_par_dT, ZKi_par_T, dZKi_par_dT, ZKe_par_T, dZKe_par_dT
real*8     :: D_prof, ZK_prof, ZKi_prof, ZKe_prof, psi_norm, theta, zeta, delta_u_x, delta_u_y, delta_ps_x, delta_ps_y
real*8     :: rhs_ij(n_var), rhs_ij_k(n_var)
real*8     :: amat(n_var,n_var), amat_k(n_var,n_var), amat_n(n_var,n_var), amat_kn(n_var,n_var)

real*8     :: v, v_x, v_y, v_s, v_t, v_p, v_ss, v_st, v_tt, v_xx, v_xy, v_yy
real*8     :: ps0, ps0_x, ps0_y, ps0_p, ps0_s, ps0_t, ps0_ss, ps0_tt, ps0_st, ps0_xx, ps0_yy, ps0_xy
real*8     :: zj0, zj0_x, zj0_y, zj0_p, zj0_s, zj0_t
real*8     :: u0, u0_x, u0_y, u0_p, u0_s, u0_t, u0_ss, u0_tt, u0_st, u0_xx, u0_xy, u0_yy
real*8     :: w0, w0_x, w0_y, w0_p, w0_s, w0_t, w0_ss, w0_st, w0_tt, w0_xx, w0_xy, w0_yy
real*8     :: r0, r0_x, r0_y, r0_p, r0_s, r0_t, r0_ss, r0_st, r0_tt, r0_xx, r0_xy, r0_yy, r0_hat, r0_x_hat, r0_y_hat, r0_corr
real*8     :: T0, T0_x, T0_y, T0_p, T0_s, T0_t, T0_ss, T0_st, T0_tt, T0_xx, T0_xy, T0_yy, T0_corr, dT0_corr_dT
real*8     :: Ti0, Ti0_x, Ti0_y, Ti0_p, Ti0_s, Ti0_t, Ti0_ss, Ti0_st, Ti0_tt, Ti0_xx, Ti0_xy, Ti0_yy, Ti0_corr, dTi0_corr_dT
real*8     :: Te0, Te0_x, Te0_y, Te0_p, Te0_s, Te0_t, Te0_ss, Te0_st, Te0_tt, Te0_xx, Te0_xy, Te0_yy, Te0_corr, dTe0_corr_dT
real*8     :: psi, psi_x, psi_y, psi_p, psi_s, psi_t, psi_ss, psi_st, psi_tt, psi_xx, psi_xy, psi_yy
real*8     :: zj, zj_x, zj_y, zj_p, zj_s, zj_t, zj_ss, zj_st, zj_tt
real*8     :: u, u_x, u_y, u_p, u_s, u_t, u_ss, u_st, u_tt, u_xx, u_xy, u_yy
real*8     :: w, w_x, w_y, w_p, w_s, w_t, w_ss, w_st, w_tt, w_xx, w_xy, w_yy
real*8     :: rho, rho_x, rho_y, rho_s, rho_t, rho_p, rho_hat, rho_x_hat, rho_y_hat, rho_ss, rho_st, rho_tt, rho_xx, rho_xy, rho_yy
real*8     :: T, T_x, T_y, T_s, T_t, T_p, T_ss, T_st, T_tt, T_xx, T_xy, T_yy
real*8     :: Ti, Ti_x, Ti_y, Ti_s, Ti_t, Ti_p, Ti_ss, Ti_st, Ti_tt, Ti_xx, Ti_xy, Ti_yy
real*8     :: Te, Te_x, Te_y, Te_s, Te_t, Te_p, Te_ss, Te_st, Te_tt, Te_xx, Te_xy, Te_yy
real*8	   :: zTi, zTi_x, zTi_y, zTe, zTe_x, zTe_y, zn_x, zn_y, Jb_0 , Jb
real*8     :: Vpar, Vpar_x, Vpar_y, Vpar_p, Vpar_s, Vpar_t, Vpar_ss, Vpar_st, Vpar_tt, Vpar_xx, Vpar_yy, Vpar_xy
real*8     :: rn0, rn0_s, rn0_t, rn0_x, rn0_y, rn0_p, rn0_ss, rn0_st, rn0_tt, rn0_xx, rn0_yy
real*8     :: rn0_hat, rn0_x_hat, rn0_y_hat, rn0_corr
real*8     :: rhon, rhon_s, rhon_t, rhon_p, rhon_x, rhon_y, rhon_ss, rhon_tt, rhon_st, rhon_xx, rhon_yy, rhon_xy
real*8	   :: dn0x, dn0y, dn0p
real*8     :: P0,  P0_s,  P0_t,  P0_x,  P0_y,  P0_p
real*8     :: Pi0, Pi0_s, Pi0_t, Pi0_x, Pi0_y, Pi0_p, Pi0_ss, Pi0_st, Pi0_tt, Pi0_xx, Pi0_xy, Pi0_yy
real*8     :: Pi0_x_rho, Pi0_xx_rho, Pi0_y_rho, Pi0_yy_rho, Pi0_xy_rho
real*8     :: Pi0_x_Ti,  Pi0_xx_Ti,  Pi0_y_Ti,  Pi0_yy_Ti,  Pi0_xy_Ti
real*8     :: Pe0, Pe0_s, Pe0_t, Pe0_x, Pe0_y, Pe0_p, Pe0_ss, Pe0_st, Pe0_tt, Pe0_xx, Pe0_xy, Pe0_yy
real*8     :: Vpar0, Vpar0_s, Vpar0_t, Vpar0_p, Vpar0_x, Vpar0_y, Vpar0_ss, Vpar0_st, Vpar0_tt, Vpar0_xx, Vpar0_yy,Vpar0_xy
real*8     :: BigR_x, vv2, eta_T, visco_T, deta_dT, d2eta_d2T, dvisco_dT, d2visco_dT2, visco_num_T, eta_num_T, W_dia, W_dia_rho, W_dia_Ti
real*8     :: eta_T_ohm, deta_dT_ohm,  deta_num_dT,  dvisco_num_dT
real*8     :: Ti0_ps0_x, Ti_ps0_x, Ti0_psi_x, Ti0_ps0_y, Ti_ps0_y, Ti0_psi_y, v_ps0_x, v_psi_x, v_ps0_y, v_psi_y
real*8     :: Te0_ps0_x, Te_ps0_x, Te0_psi_x, Te0_ps0_y, Te_ps0_y, Te0_psi_y

real*8     :: Vt0,Omega_tor0_x,Omega_tor0_y,Vt0_x,Vt0_y
real*8     :: V_source(n_gauss,n_gauss), Vt_x_psi, Vt_y_psi, Omega_tor_x_psi, Omega_tor_y_psi

real*8     :: dV_dpsi_source(n_gauss,n_gauss), dV_dz_source(n_gauss,n_gauss)
real*8     :: dV_dpsi2,dV_dz2,dV_dpsi_dz,dV_dpsi3,dV_dpsi_dz2,dV_dpsi2_dz
real*8     :: eq_zne(n_gauss,n_gauss),  eq_zT(n_gauss,n_gauss), eq_zTi(n_gauss,n_gauss), eq_zTe(n_gauss,n_gauss)
real*8     :: dn_dpsi(n_gauss,n_gauss), dn_dz,  dn_dpsi2,  dn_dz2,  dn_dpsi_dz,  dn_dpsi3,  dn_dpsi_dz2,  dn_dpsi2_dz
real*8     :: dT_dpsi(n_gauss,n_gauss), dT_dz,  dT_dpsi2,  dT_dz2,  dT_dpsi_dz,  dT_dpsi3,  dT_dpsi_dz2,  dT_dpsi2_dz
real*8     :: dTi_dpsi(n_gauss,n_gauss),dTi_dz, dTi_dpsi2, dTi_dz2, dTi_dpsi_dz, dTi_dpsi3, dTi_dpsi_dz2, dTi_dpsi2_dz
real*8     :: dTe_dpsi(n_gauss,n_gauss),dTe_dz, dTe_dpsi2, dTe_dz2, dTe_dpsi_dz, dTe_dpsi3, dTe_dpsi_dz2, dTe_dpsi2_dz

logical    :: xpoint2, use_fft
real*8     :: Btheta2, epsil, Btheta2_psi
real*8, dimension(n_gauss,n_gauss)    :: amu_neo_prof, aki_neo_prof

! neutral source
real*8     :: source_neutral
real*8     :: source_neutral_drift ! Neutral source deposited at R+drift_distance to impose plasmoid drift
real*8     :: power_dens_teleport_ju ! Teleported power density in JOREK unit (sink at R and source at R+drift)

! time normalisation
real*8     :: t_norm

! Atomic physics coefficients:
!   -Ionization
real*8     :: Sion_T, dSion_dT                                ! Ionization rate and its derivative wrt. temperature
real*8     :: coef_ion_1, coef_ion_2, coef_ion_3, S_ion_puiss ! Ionization rate parameters
real*8     :: ksiion                                          ! Ionization energy
!   -Recombination
real*8     :: Srec_T, dSrec_dT                                ! Recombination rate and its derivative wrt. temperature
real*8     :: coef_rec_1                                      ! Recombination rate parameters
!   -Radiation from injected gas/impurities
real*8     :: LradDrays_T, dLradDrays_dT                      ! Line (/rays) radiation rate and its derivative wrt. temperature
real*8     :: LradDcont_T, dLradDcont_dT                      ! Continuum (Brem.) radiation rate and its derivative wrt. T

!   -Radiation from background impurities
real*8     :: Arad_bg, Brad_bg, Crad_bg, frad_bg, dfrad_bg_dT ! Retain the hard-coded fitting for argon
real*8     :: Lrad_imp, dLrad_imp_dT                          ! Radiation rate and its derivative wrt. temperature
real*8     :: r_imp                                           ! Background impurity density in JOREK unit
integer    :: i_imp                                           ! Loop for more than one background impurity


real*8     :: in_fft(1:n_plane)
complex*16 :: out_fft(1:n_plane)
integer*8  :: plan

integer    :: max_terms_loop, i_term
real*8     :: factor(n_var,max_terms)

integer    :: i_v, i_loc, j_loc

!   -Ion-electron energy transfer
real*8     :: nu_e_bg, lambda_e_bg, dTi_e, dTe_i
real*8     :: dnu_e_bg_dTi, dnu_e_bg_dTe
real*8     :: dnu_e_bg_drho, dnu_e_bg_drhon
real*8     :: ddTi_e_dTi, ddTi_e_dTe, ddTi_e_drho, ddTi_e_drhon
real*8     :: ddTe_i_dTi, ddTe_i_dTe, ddTe_i_drho, ddTe_i_drhon
real*8     :: Te_eV, Te_corr_eV, dTe_corr_eV_dT              ! Electron temperature in eV
real*8     :: ne_SI                                          ! Electron density in SI unit
real*8     :: dr0_corr_dn

! --- General T for T-denpendent functions
real*8     :: T_or_Te, T_or_Te_corr, T_or_Te_0, dT_or_Te_corr_dT

! --- Factor to use the conservative form or not of the momentum equation
real*8     :: fact_conservative_u = 1.d0

#define DIM1 n_plane
#define DIM2 1:n_vertex_max*n_var*(n_order+1)

real*8, dimension(DIM1, DIM2, DIM2) :: ELM_p
real*8, dimension(DIM1, DIM2, DIM2) :: ELM_n
real*8, dimension(DIM1, DIM2, DIM2) :: ELM_k
real*8, dimension(DIM1, DIM2, DIM2) :: ELM_kn
real*8, dimension(DIM1, DIM2)       :: RHS_p
real*8, dimension(DIM1, DIM2)       :: RHS_k

real*8, dimension(n_gauss,n_gauss)    :: x_g, x_s, x_t, x_ss, x_st, x_tt
real*8, dimension(n_gauss,n_gauss)    :: y_g, y_s, y_t, y_ss, y_st, y_tt

real*8, dimension(n_plane,n_var,n_gauss,n_gauss) :: eq_g, eq_s, eq_t, eq_p, eq_ss, eq_st, eq_tt
real*8, dimension(n_plane,n_var,n_gauss,n_gauss) :: delta_g, delta_s, delta_t

real*8, dimension(n_tor,n_plane) :: HHZ, HHZ_p, HHZ_pp

!  --- For shock capturing stabilization
real*8     :: midp_edge1(1:2), midp_edge2(1:2), midp_edge3(1:2), midp_edge4(1:2)
real*8     :: len1, len2, h_e
real*8     :: Ptot, Ptot_x,  Ptot_y,  Ptot_p, Ptot_corr
real*8     :: f_p, d_p, tau_sc, R_rho, R_pi, R_pe, R_p, R_rhon, my_zero = 0.d0
real*8     :: s_p, src_rho, src_p, src_pi, src_pe, src_rhon

if (present(get_terms)) then
  max_terms_loop = max_terms
else
  max_terms_loop = 1
endif

ELM_p = 0.d0
ELM_n = 0.d0
ELM_k = 0.d0
ELM_kn = 0.d0
RHS_p = 0.d0
RHS_k = 0.d0
ELM   = 0.d0
RHS   = 0.d0

epsil=1.d-3

! --- Decide whether or not use the conservative form of the momentum equation
! --- (conservative form not ready yet with diamagnetic flows)
if (tauIC > 0.d0) fact_conservative_u = 0.d0

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

! --- Toroidal functions            
if (use_fft) then
  HHZ    = 1.d0
  HHZ_p  = 1.d0
  HHZ_pp = 1.d0
else
  do in = 1,n_tor
    do mp=1,n_plane
      HHZ   (in,mp) = HZ   (in,mp)
      HHZ_p (in,mp) = HZ_p (in,mp)
      HHZ_pp(in,mp) = HZ_pp(in,mp)
    enddo
  enddo
endif

rhs_ij  = 0.d0; rhs_ij_k  = 0.d0; 
amat    = 0.d0; amat_k    = 0.d0; amat_n = 0.d0; amat_kn = 0.d0
!---------------------------------------------------- value of (x,y) and derivatives on Gaussian points
x_g  = 0.d0; x_s  = 0.d0; x_t  = 0.d0; x_st  = 0.d0; x_ss  = 0.d0; x_tt  = 0.d0;
y_g  = 0.d0; y_s  = 0.d0; y_t  = 0.d0; y_st  = 0.d0; y_ss  = 0.d0; y_tt  = 0.d0;
eq_g = 0.d0; eq_s = 0.d0; eq_t = 0.d0; eq_st = 0.d0; eq_ss = 0.d0; eq_tt = 0.d0; eq_p = 0.d0;

delta_g = 0.d0; delta_s = 0.d0; delta_t = 0.d0

current_source  = 0.d0
particle_source = 0.d0
heat_source_i   = 0.d0
heat_source_e   = 0.d0
heat_source     = 0.d0
V_source        = 0.d0
dV_dpsi_source  = 0.d0
dV_dz_source    = 0.d0
eq_zne          = 0.d0
eq_zTi          = 0.d0         
eq_zTe          = 0.d0
eq_zT           = 0.d0

amu_neo_prof   = 0.d0
aki_neo_prof   = 0.d0

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

! approximate estimate of the element length h_e
! needed for Shock capturing stabilzation
midp_edge1(:) =  0.5d0 * ( nodes(1)%x(1,1,:) + nodes(2)%x(1,1,:) )
midp_edge2(:) =  0.5d0 * ( nodes(2)%x(1,1,:) + nodes(3)%x(1,1,:) )
midp_edge3(:) =  0.5d0 * ( nodes(3)%x(1,1,:) + nodes(4)%x(1,1,:) )
midp_edge4(:) =  0.5d0 * ( nodes(4)%x(1,1,:) + nodes(1)%x(1,1,:) )

len1 = sqrt( (midp_edge1(1)-midp_edge3(1))**2 + (midp_edge1(2)-midp_edge3(2))**2 )
len2 = sqrt( (midp_edge2(1)-midp_edge4(1))**2 + (midp_edge2(2)-midp_edge4(2))**2 )
h_e = dmin1(len1, len2)

do ms=1, n_gauss
  do mt=1, n_gauss

    if (keep_current_prof) &
      call current(xpoint2, xcase2, x_g(ms,mt),y_g(ms,mt), Z_xpoint, eq_g(1,var_psi,ms,mt),psi_axis,psi_bnd,current_source(ms,mt))
  
    if ( with_TiTe ) then
      call sources(xpoint2, xcase2, y_g(ms,mt), Z_xpoint, eq_g(1,var_psi,ms,mt),psi_axis,psi_bnd, &
                   particle_source(ms,mt),heat_source_i(ms,mt),heat_source_e(ms,mt))
    else
      call sources(xpoint2, xcase2, y_g(ms,mt), Z_xpoint, eq_g(1,var_psi,ms,mt),psi_axis,psi_bnd, &
                   particle_source(ms,mt),heat_source(ms,mt))
    end if

    ! Source of parallel velocity
    if ( ( abs(V_0) .ge. 1.e-12 ) .or. ( num_rot ) ) then
      call velocity(xpoint2, xcase2, y_g(ms,mt), Z_xpoint, eq_g(1,var_psi,ms,mt), psi_axis, psi_bnd, V_source(ms,mt), &
                    dV_dpsi_source(ms,mt),dV_dz_source(ms,mt),dV_dpsi2,dV_dz2,dV_dpsi_dz,dV_dpsi3,dV_dpsi_dz2, dV_dpsi2_dz)
    endif

    call density(xpoint2, xcase2, y_g(ms,mt), Z_xpoint, eq_g(1,var_psi,ms,mt),psi_axis,psi_bnd,eq_zne(ms,mt), &
                 dn_dpsi(ms,mt),dn_dz,dn_dpsi2,dn_dz2,dn_dpsi_dz,dn_dpsi3,dn_dpsi_dz2, dn_dpsi2_dz)

    if ( with_TiTe ) then
      call temperature_i(xpoint2, xcase2, y_g(ms,mt), Z_xpoint, eq_g(1,var_psi,ms,mt),psi_axis,psi_bnd, eq_zTi(ms,mt), &
                       dTi_dpsi(ms,mt),dTi_dz,dTi_dpsi2,dTi_dz2,dTi_dpsi_dz,dTi_dpsi3,dTi_dpsi_dz2, dTi_dpsi2_dz)
  
      call temperature_e(xpoint2, xcase2, y_g(ms,mt), Z_xpoint, eq_g(1,var_psi,ms,mt),psi_axis,psi_bnd,eq_zTe(ms,mt), &
                       dTe_dpsi(ms,mt),dTe_dz,dTe_dpsi2,dTe_dz2,dTe_dpsi_dz,dTe_dpsi3,dTe_dpsi_dz2, dTe_dpsi2_dz)
    else
      call temperature(xpoint2, xcase2, y_g(ms,mt), Z_xpoint, eq_g(1,var_psi,ms,mt),psi_axis,psi_bnd,eq_zT(ms,mt), &
                       dT_dpsi(ms,mt),dT_dz,dT_dpsi2,dT_dz2,dT_dpsi_dz,dT_dpsi3,dT_dpsi_dz2, dT_dpsi2_dz)
    end if

    if ( NEO ) then 
      if (num_neo_file) then
        call neo_coef( xpoint2, xcase2, y_g(ms,mt), Z_xpoint, eq_g(1,var_psi,ms,mt),psi_axis,psi_bnd, amu_neo_prof(ms,mt), aki_neo_prof(ms,mt))
      else
        amu_neo_prof(ms,mt) = amu_neo_const
        aki_neo_prof(ms,mt) = aki_neo_const
      endif
    endif

  enddo
enddo

!--------------------------------------------------- sum over the Gaussian integration points
do i=1,n_vertex_max
  do j=1,n_order+1

    if (.not. present(get_terms)) then
      ELM_p(:,:,1:n_var)  = 0
      ELM_n(:,:,1:n_var)  = 0
      ELM_k(:,:,1:n_var)  = 0
      ELM_kn(:,:,1:n_var) = 0
    endif

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

        do mp = 1, n_plane

          ps0    = eq_g(mp,var_psi,ms,mt)
          ps0_x  = (   y_t(ms,mt) * eq_s(mp,var_psi,ms,mt) - y_s(ms,mt) * eq_t(mp,var_psi,ms,mt) ) / xjac
          ps0_y  = ( - x_t(ms,mt) * eq_s(mp,var_psi,ms,mt) + x_s(ms,mt) * eq_t(mp,var_psi,ms,mt) ) / xjac
          ps0_p  = eq_p(mp,var_psi,ms,mt)
          ps0_s  = eq_s(mp,var_psi,ms,mt)
          ps0_t  = eq_t(mp,var_psi,ms,mt)
          ps0_ss = eq_ss(mp,var_psi,ms,mt)
          ps0_tt = eq_tt(mp,var_psi,ms,mt)
          ps0_st = eq_st(mp,var_psi,ms,mt)

          u0    = eq_g(mp,var_u,ms,mt)
          u0_x  = (   y_t(ms,mt) * eq_s(mp,var_u,ms,mt) - y_s(ms,mt) * eq_t(mp,var_u,ms,mt) ) / xjac
          u0_y  = ( - x_t(ms,mt) * eq_s(mp,var_u,ms,mt) + x_s(ms,mt) * eq_t(mp,var_u,ms,mt) ) / xjac
          u0_p  = eq_p(mp,var_u,ms,mt)
          u0_s  = eq_s(mp,var_u,ms,mt)
          u0_t  = eq_t(mp,var_u,ms,mt)
          u0_ss = eq_ss(mp,var_u,ms,mt)
          u0_tt = eq_tt(mp,var_u,ms,mt)
          u0_st = eq_st(mp,var_u,ms,mt)

          vv2   = BigR**2 *  ( u0_x * u0_x + u0_y *u0_y  )

          zj0   = eq_g(mp,var_zj,ms,mt)
          zj0_x = (   y_t(ms,mt) * eq_s(mp,var_zj,ms,mt) - y_s(ms,mt) * eq_t(mp,var_zj,ms,mt) ) / xjac
          zj0_y = ( - x_t(ms,mt) * eq_s(mp,var_zj,ms,mt) + x_s(ms,mt) * eq_t(mp,var_zj,ms,mt) ) / xjac
          zj0_p = eq_p(mp,var_zj,ms,mt)
          zj0_s = eq_s(mp,var_zj,ms,mt)
          zj0_t = eq_t(mp,var_zj,ms,mt)


          w0    = eq_g(mp,var_w,ms,mt)
          w0_x  = (   y_t(ms,mt) * eq_s(mp,var_w,ms,mt) - y_s(ms,mt) * eq_t(mp,var_w,ms,mt) ) / xjac
          w0_y  = ( - x_t(ms,mt) * eq_s(mp,var_w,ms,mt) + x_s(ms,mt) * eq_t(mp,var_w,ms,mt) ) / xjac
          w0_p  = eq_p(mp,var_w,ms,mt)
          w0_s  = eq_s(mp,var_w,ms,mt)
          w0_t  = eq_t(mp,var_w,ms,mt)
          w0_ss = eq_ss(mp,var_w,ms,mt)
          w0_tt = eq_tt(mp,var_w,ms,mt)
          w0_st = eq_st(mp,var_w,ms,mt)

          r0    = eq_g(mp,var_rho,ms,mt)
          r0_corr = corr_neg_dens(r0)
          dr0_corr_dn = dcorr_neg_dens_drho(r0)
          r0_x  = (   y_t(ms,mt) * eq_s(mp,var_rho,ms,mt) - y_s(ms,mt) * eq_t(mp,var_rho,ms,mt) ) / xjac
          r0_y  = ( - x_t(ms,mt) * eq_s(mp,var_rho,ms,mt) + x_s(ms,mt) * eq_t(mp,var_rho,ms,mt) ) / xjac
          r0_p  = eq_p(mp,var_rho,ms,mt)
          r0_s  = eq_s(mp,var_rho,ms,mt)
          r0_t  = eq_t(mp,var_rho,ms,mt)
          r0_ss = eq_ss(mp,var_rho,ms,mt)
          r0_st = eq_st(mp,var_rho,ms,mt)
          r0_tt = eq_tt(mp,var_rho,ms,mt)

          r0_hat   = BigR**2 * r0
          r0_x_hat = 2.d0 * BigR * BigR_x  * r0 + BigR**2 * r0_x
          r0_y_hat = BigR**2 * r0_y

          if ( with_TiTe ) then ! ******************************************************************
            
            T0    = 0.d0
            T0_x  = 0.d0
            T0_y  = 0.d0
            T0_p  = 0.d0
            T0_s  = 0.d0
            T0_t  = 0.d0
            T0_ss = 0.d0
            T0_tt = 0.d0
            T0_st = 0.d0
            
            T0_corr     = 0.d0
            dT0_corr_dT = 0.d0
           
            Ti0    = eq_g(mp,var_Ti,ms,mt)
            Ti0_x  = (   y_t(ms,mt) * eq_s(mp,var_Ti,ms,mt) - y_s(ms,mt) * eq_t(mp,var_Ti,ms,mt) ) / xjac
            Ti0_y  = ( - x_t(ms,mt) * eq_s(mp,var_Ti,ms,mt) + x_s(ms,mt) * eq_t(mp,var_Ti,ms,mt) ) / xjac
            Ti0_p  = eq_p(mp,var_Ti,ms,mt)
            Ti0_s  = eq_s(mp,var_Ti,ms,mt)
            Ti0_t  = eq_t(mp,var_Ti,ms,mt)
            Ti0_ss = eq_ss(mp,var_Ti,ms,mt)
            Ti0_tt = eq_tt(mp,var_Ti,ms,mt)
            Ti0_st = eq_st(mp,var_Ti,ms,mt)
           
            Ti0_corr     = corr_neg_temp(Ti0) ! For use in eta(T), visco(T), ...
            dTi0_corr_dT = dcorr_neg_temp_dT(Ti0) ! Improve the correction
           
            Te0    = eq_g(mp,var_Te,ms,mt)
            Te0_x  = (   y_t(ms,mt) * eq_s(mp,var_Te,ms,mt) - y_s(ms,mt) * eq_t(mp,var_Te,ms,mt) ) / xjac
            Te0_y  = ( - x_t(ms,mt) * eq_s(mp,var_Te,ms,mt) + x_s(ms,mt) * eq_t(mp,var_Te,ms,mt) ) / xjac
            Te0_p  = eq_p(mp,var_Te,ms,mt)
            Te0_s  = eq_s(mp,var_Te,ms,mt)
            Te0_t  = eq_t(mp,var_Te,ms,mt)
            Te0_ss = eq_ss(mp,var_Te,ms,mt)
            Te0_tt = eq_tt(mp,var_Te,ms,mt)
            Te0_st = eq_st(mp,var_Te,ms,mt)
           
            Te0_corr     = corr_neg_temp(Te0) ! For use in eta(T), visco(T), ...
            dTe0_corr_dT = dcorr_neg_temp_dT(Te0) ! Improve the correction

            zTi   =   eq_zTi(ms,mt)
            zTi_x = dTi_dpsi(ms,mt) * ps0_x
            zTi_y = dTi_dpsi(ms,mt) * ps0_y
            zTe   =   eq_zTe(ms,mt)
            zTe_x = dTe_dpsi(ms,mt) * ps0_x
            zTe_y = dTe_dpsi(ms,mt) * ps0_y
            zn_x  =  dn_dpsi(ms,mt) * ps0_x
            zn_y  =  dn_dpsi(ms,mt) * ps0_y

          else ! (with_TiTe = .f.), i.e. with single temperature *****************************************

            T0    = eq_g(mp,var_T,ms,mt)
            T0_x  = (   y_t(ms,mt) * eq_s(mp,var_T,ms,mt) - y_s(ms,mt) * eq_t(mp,var_T,ms,mt) ) / xjac
            T0_y  = ( - x_t(ms,mt) * eq_s(mp,var_T,ms,mt) + x_s(ms,mt) * eq_t(mp,var_T,ms,mt) ) / xjac
            T0_p  = eq_p(mp,var_T,ms,mt)
            T0_s  = eq_s(mp,var_T,ms,mt)
            T0_t  = eq_t(mp,var_T,ms,mt)
            T0_ss = eq_ss(mp,var_T,ms,mt)
            T0_tt = eq_tt(mp,var_T,ms,mt)
            T0_st = eq_st(mp,var_T,ms,mt)
           
            T0_corr     = corr_neg_temp(T0) ! For use in eta(T), visco(T), ...
            dT0_corr_dT = dcorr_neg_temp_dT(T0) ! Improve the correction

            Ti0    = T0    / 2.d0
            Ti0_x  = T0_x  / 2.d0
            Ti0_y  = T0_y  / 2.d0
            Ti0_p  = T0_p  / 2.d0
            Ti0_s  = T0_s  / 2.d0
            Ti0_t  = T0_t  / 2.d0
            Ti0_ss = T0_ss / 2.d0
            Ti0_tt = T0_tt / 2.d0
            Ti0_st = T0_st / 2.d0

            Ti0_corr     = corr_neg_temp(Ti0) ! For use in eta(T), visco(T), ...
            dTi0_corr_dT = dcorr_neg_temp_dT(Ti0) ! Improve the correction
           
            Te0    = Ti0
            Te0_x  = Ti0_x
            Te0_y  = Ti0_y
            Te0_p  = Ti0_p
            Te0_s  = Ti0_s
            Te0_t  = Ti0_t
            Te0_ss = Ti0_ss
            Te0_tt = Ti0_tt
            Te0_st = Ti0_st

            Te0_corr     = corr_neg_temp(Te0) ! For use in eta(T), visco(T), ...
            dTe0_corr_dT = dcorr_neg_temp_dT(Te0) ! Improve the correction

            zTi   =   eq_zT(ms,mt)         / 2.d0
            zTi_x = dT_dpsi(ms,mt) * ps0_x / 2.d0
            zTi_y = dT_dpsi(ms,mt) * ps0_y / 2.d0
            zTe   = zTi
            zTe_x = zTi_x
            zTe_y = zTi_y
            zn_x  = dn_dpsi(ms,mt) * ps0_x
            zn_y  = dn_dpsi(ms,mt) * ps0_y
          end if ! (with_TiTe) *********************************************************************

          if ( with_vpar ) then
            Vpar0    = eq_g(mp,var_Vpar,ms,mt)
            Vpar0_x  = (   y_t(ms,mt) * eq_s(mp,var_Vpar,ms,mt) - y_s(ms,mt) * eq_t(mp,var_Vpar,ms,mt) ) / xjac
            Vpar0_y  = ( - x_t(ms,mt) * eq_s(mp,var_Vpar,ms,mt) + x_s(ms,mt) * eq_t(mp,var_Vpar,ms,mt) ) / xjac
            Vpar0_p  = eq_p(mp,var_Vpar,ms,mt)
            Vpar0_s  = eq_s(mp,var_Vpar,ms,mt)
            Vpar0_t  = eq_t(mp,var_Vpar,ms,mt)
            Vpar0_ss = eq_ss(mp,var_Vpar,ms,mt)
            Vpar0_st = eq_st(mp,var_Vpar,ms,mt)
            Vpar0_tt = eq_tt(mp,var_Vpar,ms,mt)
          else
            Vpar0    = 0.d0
            Vpar0_x  = 0.d0
            Vpar0_y  = 0.d0
            Vpar0_p  = 0.d0
            Vpar0_s  = 0.d0
            Vpar0_t  = 0.d0
            Vpar0_ss = 0.d0
            Vpar0_st = 0.d0
            Vpar0_tt = 0.d0
          end if

          if (with_neutrals) then
            rn0      = eq_g(mp,var_rhon,ms,mt)
            rn0_x    = (   y_t(ms,mt) * eq_s(mp,var_rhon,ms,mt) - y_s(ms,mt) * eq_t(mp,var_rhon,ms,mt) ) / xjac    
            rn0_y    = ( - x_t(ms,mt) * eq_s(mp,var_rhon,ms,mt) + x_s(ms,mt) * eq_t(mp,var_rhon,ms,mt) ) / xjac   
            rn0_p    = eq_p(mp,var_rhon,ms,mt)                                                             
            rn0_s    = eq_s(mp,var_rhon,ms,mt)                                                             
            rn0_t    = eq_t(mp,var_rhon,ms,mt)                                                             
            rn0_ss   = eq_ss(mp,var_rhon,ms,mt)                                                            
            rn0_st   = eq_st(mp,var_rhon,ms,mt)                                                            
            rn0_tt   = eq_tt(mp,var_rhon,ms,mt)  
            rn0_corr = corr_neg_dens(rn0, (/ 0.d-5, 1.d-5 /)) ! Correction for negative rn0 ...
          else
            rn0      = 0.d0
            rn0_x    = 0.d0  
            rn0_y    = 0.d0 
            rn0_p    = 0.d0
            rn0_s    = 0.d0
            rn0_t    = 0.d0
            rn0_ss   = 0.d0
            rn0_st   = 0.d0
            rn0_tt   = 0.d0
            rn0_corr = 0.d0 
          endif
     
          rn0_xx = (rn0_ss * y_t(ms,mt)**2 - 2.d0*rn0_st * y_s(ms,mt)*y_t(ms,mt) + rn0_tt * y_s(ms,mt)**2     &
            + rn0_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                              &
            + rn0_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) ) )    / xjac**2               &
            - xjac_x * (rn0_s* y_t(ms,mt) - rn0_t * y_s(ms,mt))  / xjac**2

          rn0_yy = (rn0_ss * x_t(ms,mt)**2 - 2.d0*rn0_st * x_s(ms,mt)*x_t(ms,mt) + rn0_tt * x_s(ms,mt)**2     &
            + rn0_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                              &
            + rn0_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) ) )       / xjac**2            &
            - xjac_y * (- rn0_s * x_t(ms,mt) + rn0_t * x_s(ms,mt) )  / xjac**2

          rn0_hat   = BigR**2 * rn0                                                        
          rn0_x_hat = 2.d0 * BigR * BigR_x  * rn0 + BigR**2 * rn0_x                             
          rn0_y_hat = BigR**2 * rn0_y                                                            

          Pi0    = r0    * Ti0
          Pi0_x  = r0_x  * Ti0 + r0 * Ti0_x
          Pi0_y  = r0_y  * Ti0 + r0 * Ti0_y
          Pi0_s  = r0_s  * Ti0 + r0 * Ti0_s
          Pi0_t  = r0_t  * Ti0 + r0 * Ti0_t
          Pi0_p  = r0_p  * Ti0 + r0 * Ti0_p
          Pi0_ss = r0_ss * Ti0 + 2.d0 * r0_s * Ti0_s + r0 * Ti0_ss
          Pi0_tt = r0_tt * Ti0 + 2.d0 * r0_t * Ti0_t + r0 * Ti0_tt
          Pi0_st = r0_st * Ti0 + r0_s * Ti0_t + r0_t * Ti0_s + r0 * Ti0_st

          Pe0    = r0    * Te0
          Pe0_x  = r0_x  * Te0 + r0 * Te0_x
          Pe0_y  = r0_y  * Te0 + r0 * Te0_y
          Pe0_s  = r0_s  * Te0 + r0 * Te0_s
          Pe0_t  = r0_t  * Te0 + r0 * Te0_t
          Pe0_p  = r0_p  * Te0 + r0 * Te0_p
          Pe0_ss = r0_ss * Te0 + 2.d0 * r0_s * Te0_s + r0 * Te0_ss
          Pe0_tt = r0_tt * Te0 + 2.d0 * r0_t * Te0_t + r0 * Te0_tt
          Pe0_st = r0_st * Te0 + r0_s * Te0_t + r0_t * Te0_s + r0 * Te0_st

          P0     = Pi0   + Pe0
          P0_s   = Pi0_s + Pe0_s
          P0_t   = Pi0_t + Pe0_t
          P0_p   = Pi0_p + Pe0_p
          P0_x   = Pi0_x + Pe0_x
          P0_y   = Pi0_y + Pe0_y

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

          Ti0_xx = (Ti0_ss * y_t(ms,mt)**2 - 2.d0*Ti0_st * y_s(ms,mt)*y_t(ms,mt) + Ti0_tt * y_s(ms,mt)**2  &
                + Ti0_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                            &
                + Ti0_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) ) )    / xjac**2             &     
                - xjac_x * (Ti0_s * y_t(ms,mt) - Ti0_t * y_s(ms,mt))  / xjac**2

          Ti0_yy = (Ti0_ss * x_t(ms,mt)**2 - 2.d0*Ti0_st * x_s(ms,mt)*x_t(ms,mt) + Ti0_tt * x_s(ms,mt)**2  &
                + Ti0_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                            &
                + Ti0_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) ) )    / xjac**2             &
                - xjac_y * (- Ti0_s * x_t(ms,mt) + Ti0_t * x_s(ms,mt) )  / xjac**2

          Ti0_xy = (- Ti0_ss * y_t(ms,mt)*x_t(ms,mt) - Ti0_tt * x_s(ms,mt)*y_s(ms,mt)                     &
                   + Ti0_st * (y_s(ms,mt)*x_t(ms,mt)  + y_t(ms,mt)*x_s(ms,mt)  )                        &
                   - Ti0_s  * (x_st(ms,mt)*y_t(ms,mt) - x_tt(ms,mt)*y_s(ms,mt) )                        &
                   - Ti0_t * (x_st(ms,mt)*y_s(ms,mt)  - x_ss(ms,mt)*y_t(ms,mt) )  )  / xjac**2          &
                - xjac_x * (- Ti0_s * x_t(ms,mt) + Ti0_t * x_s(ms,mt) )   / xjac**2

          Te0_xx = (Te0_ss * y_t(ms,mt)**2 - 2.d0*Te0_st * y_s(ms,mt)*y_t(ms,mt) + Te0_tt * y_s(ms,mt)**2  &
                + Te0_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                            &
                + Te0_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) ) )    / xjac**2             &     
                - xjac_x * (Te0_s * y_t(ms,mt) - Te0_t * y_s(ms,mt))  / xjac**2

          Te0_yy = (Te0_ss * x_t(ms,mt)**2 - 2.d0*Te0_st * x_s(ms,mt)*x_t(ms,mt) + Te0_tt * x_s(ms,mt)**2  &
                + Te0_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                            &
                + Te0_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) ) )    / xjac**2             &
                - xjac_y * (- Te0_s * x_t(ms,mt) + Te0_t * x_s(ms,mt) )  / xjac**2

          Te0_xy = (- Te0_ss * y_t(ms,mt)*x_t(ms,mt) - Te0_tt * x_s(ms,mt)*y_s(ms,mt)                     &
                   + Te0_st * (y_s(ms,mt)*x_t(ms,mt)  + y_t(ms,mt)*x_s(ms,mt)  )                        &
                   - Te0_s  * (x_st(ms,mt)*y_t(ms,mt) - x_tt(ms,mt)*y_s(ms,mt) )                        &
                   - Te0_t * (x_st(ms,mt)*y_s(ms,mt)  - x_ss(ms,mt)*y_t(ms,mt) )  )  / xjac**2          &
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


          Pi0_xx = r0_xx * Ti0 + 2.d0 * r0_x * Ti0_x + r0 * Ti0_xx
          Pi0_yy = r0_yy * Ti0 + 2.d0 * r0_y * Ti0_y + r0 * Ti0_yy
          Pi0_xy = r0_xy * Ti0 + r0_x * Ti0_y + r0_y * Ti0_x + r0 * Ti0_xy

          Pe0_xx = r0_xx * Te0 + 2.d0 * r0_x * Te0_x + r0 * Te0_xx
          Pe0_yy = r0_yy * Te0 + 2.d0 * r0_y * Te0_y + r0 * Te0_yy
          Pe0_xy = r0_xy * Te0 + r0_x * Te0_y + r0_y * Te0_x + r0 * Te0_xy

          Ti0_ps0_x = Ti0_xx * ps0_y - Ti0_xy * ps0_x + Ti0_x * ps0_xy - Ti0_y * ps0_xx
          Ti0_ps0_y = Ti0_xy * ps0_y - Ti0_yy * ps0_x + Ti0_x * ps0_yy - Ti0_y * ps0_xy
          Te0_ps0_x = Te0_xx * ps0_y - Te0_xy * ps0_x + Te0_x * ps0_xy - Te0_y * ps0_xx
          Te0_ps0_y = Te0_xy * ps0_y - Te0_yy * ps0_x + Te0_x * ps0_yy - Te0_y * ps0_xy

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


          delta_u_x = (   y_t(ms,mt) * delta_s(mp,var_u,ms,mt) - y_s(ms,mt) * delta_t(mp,var_u,ms,mt) ) / xjac
          delta_u_y = ( - x_t(ms,mt) * delta_s(mp,var_u,ms,mt) + x_s(ms,mt) * delta_t(mp,var_u,ms,mt) ) / xjac

          delta_ps_x = (   y_t(ms,mt) * delta_s(mp,var_psi,ms,mt) - y_s(ms,mt) * delta_t(mp,var_psi,ms,mt) ) / xjac
          delta_ps_y = ( - x_t(ms,mt) * delta_s(mp,var_psi,ms,mt) + x_s(ms,mt) * delta_t(mp,var_psi,ms,mt) ) / xjac

          if ( with_TiTe ) then ! ******************************************************************

            ! --- Temperature dependent parallel heat diffusivity
            if ( ZKpar_T_dependent ) then
              ZKi_par_T   = ZK_i_par * (Ti0_corr/Ti_0)**(+2.5d0)              ! temperature dependent parallel conductivity
              dZKi_par_dT = ZK_i_par * (2.5d0)  * Ti0_corr**(+1.5d0) * Ti_0**(-2.5d0) * dTi0_corr_dT
              if (ZKi_par_T .gt. ZK_par_max) then
                ZKi_par_T   = Zk_par_max
                dZKi_par_dT = 0.d0
              endif
              if (Ti0 .lt. T_min) then
                ZKi_par_T   = ZK_i_par * (max(Ti0,T_min)/Ti_0)**(+2.5d0)
                dZKi_par_dT = 0.d0
              endif
              
              ZKe_par_T   = ZK_e_par * (Te0_corr/Te_0)**(+2.5d0)              ! temperature dependent parallel conductivity
              dZKe_par_dT = ZK_e_par * (2.5d0)  * Te0_corr**(+1.5d0) * Te_0**(-2.5d0) * dTe0_corr_dT
              if (ZKe_par_T .gt. ZK_par_max) then
                ZKe_par_T   = Zk_par_max
                dZKe_par_dT = 0.d0
              endif
              if (Te0 .lt. T_min) then
                ZKe_par_T   = ZK_e_par * (max(Te0,T_min)/Te_0)**(+2.5d0)
                dZKe_par_dT = 0.d0
              endif
            else
              ZKi_par_T   = ZK_i_par                                            ! parallel conductivity
              dZKi_par_dT = 0.d0
              ZKe_par_T   = ZK_e_par                                            ! parallel conductivity
              dZKe_par_dT = 0.d0
            endif
            
            ! --- Ion-electron energy transfer
            if (thermalization) then
              ! Te in eV:
              Te_corr_eV     = Te0_corr/(EL_CHG*MU_ZERO*central_density*1.d20)
              dTe_corr_eV_dT = dTe0_corr_dT/(EL_CHG*MU_ZERO*central_density*1.d20)
              
              ne_SI          = r0_corr * 1.d20 * central_density ! electron density (SI)
              if (ne_SI < 1.d16) ne_SI = 1.d16 ! To prevent absurd number in the coulomb lambda
              
              lambda_e_bg  = 23. - log((ne_SI*1.d-6)**0.5*Te_corr_eV**(-1.5)) ! Assuming bg_charge is 1! 
              nu_e_bg      = 1.8d-19*(1.d6*MASS_ELECTRON*MASS_PROTON*central_mass) ** 0.5&
                             * (1.d14*central_density*r0_corr) * lambda_e_bg &
                             / (1.d3*(MASS_ELECTRON*Ti0_corr+Te0_corr*MASS_PROTON*central_mass)&
                             / (EL_CHG * MU_ZERO * central_density * 1.d20)) ** 1.5 ! Assuming bg_charge is 1!
              
              if (nu_e_bg < 0.)  nu_e_bg  = 0.
              
              !Converting the energy transfer rate from s^-1 to JOREK unit
              t_norm   = sqrt(MU_ZERO * central_mass * MASS_PROTON * central_density * 1.d20)
              nu_e_bg  = nu_e_bg * t_norm    
              
              dTe_i    = nu_e_bg * (Ti0_corr - Te0_corr)
              dTi_e    = -dTe_i
              
              !Calculating the density and temperature derivative for amats
              !We negelect the coulomb log's dericatives due to their smallness
              
              dnu_e_bg_dTi    = -1.5*MASS_ELECTRON*nu_e_bg*dTi0_corr_dT / (MASS_ELECTRON*Ti0_corr + MASS_PROTON*central_mass*Te0_corr)
              dnu_e_bg_dTe    = -1.5*MASS_PROTON*central_mass*nu_e_bg*dTe0_corr_dT &
                                / (MASS_ELECTRON*Ti0_corr + MASS_PROTON*central_mass*Te0_corr)
             
              dnu_e_bg_drho   = nu_e_bg * dr0_corr_dn / r0_corr
              
              ddTe_i_dTi      = dnu_e_bg_dTi  * (Ti0_corr - Te0_corr) + nu_e_bg
              ddTe_i_dTe      = dnu_e_bg_dTe  * (Ti0_corr - Te0_corr) - nu_e_bg
              ddTe_i_drho     = dnu_e_bg_drho * (Ti0_corr - Te0_corr)
              
              ddTi_e_dTi      = -ddTe_i_dTi
              ddTi_e_dTe      = -ddTe_i_dTe
              ddTi_e_drho     = -ddTe_i_drho
            else
              dTe_i = 0.
              dTi_e = 0.
              ddTe_i_dTi = 0.
              ddTe_i_dTe = 0.
              ddTe_i_drho = 0.
              ddTi_e_dTi = 0.
              ddTi_e_dTe = 0.
              ddTi_e_drho = 0.
            endif

            ! ---Temperature parameters used for general T-dependent functions (eta, visco, etc)
            T_or_Te          = Te0
            T_or_Te_corr     = Te0_corr
            T_or_Te_0        = Te_0
            dT_or_Te_corr_dT = dTe0_corr_dT
            
          else ! (with_TiTe = .f.), i.e. with single temperature *****************************************

            ! --- Temperature dependent parallel heat diffusivity
            if ( ZKpar_T_dependent ) then
              ZK_par_T   = ZK_par * (T0_corr/T_0)**(+2.5d0)              ! temperature dependent parallel conductivity
              dZK_par_dT = ZK_par * (2.5d0)  * T0_corr**(+1.5d0) * T_0**(-2.5d0) * dT0_corr_dT
              if (ZK_par_T .gt. ZK_par_max) then
                ZK_par_T   = Zk_par_max
                dZK_par_dT = 0.d0
              endif
              if (T0 .lt. T_min) then
                ZK_par_T   = ZK_par * (max(T0,T_min)/T_0)**(+2.5d0)
                dZK_par_dT = 0.d0
              endif
            else
              ZK_par_T   = ZK_par                                            ! parallel conductivity
              dZK_par_dT = 0.d0
            endif

            ! --- Temperature parameters used for general T-dependent functions (eta, visco, etc)
            T_or_Te          = T0
            T_or_Te_corr     = T0_corr
            T_or_Te_0        = T_0
            dT_or_Te_corr_dT = dT0_corr_dT
            
          end if ! (with_TiTe) *********************************************************************

          ! --- Eta
          if ( eta_T_dependent .and. T_or_Te_corr <= T_max_eta) then
            eta_T     =   eta   * (T_or_Te_corr/T_or_Te_0)**(-1.5d0)
            deta_dT   = - eta   * (1.5d0)  * T_or_Te_corr**(-2.5d0) * T_or_Te_0**(1.5d0)
            d2eta_d2T =   eta   * (3.75d0) * T_or_Te_corr**(-3.5d0) * T_or_Te_0**(1.5d0)
          else if ( eta_T_dependent .and. T_or_Te_corr > T_max_eta) then
            eta_T     = eta   * (T_max_eta/T_or_Te_0)**(-1.5d0)
            deta_dT   = 0.
            d2eta_d2T = 0.     
          else
            eta_T     = eta
            deta_dT   = 0.d0
            d2eta_d2T = 0.d0
          end if
          
          ! --- Eta for ohmic heating
          if ( eta_T_dependent .and. T_or_Te_corr <= T_max_eta_ohm) then
            eta_T_ohm     =   eta_ohmic   * (T_or_Te_corr/T_or_Te_0)**(-1.5d0)
            deta_dT_ohm   = - eta_ohmic   * (1.5d0)  * T_or_Te_corr**(-2.5d0) * T_or_Te_0**(1.5d0)
          else if ( eta_T_dependent .and. T_or_Te_corr > T_max_eta_ohm) then
            eta_T_ohm     =   eta_ohmic   * (T_max_eta_ohm/T_or_Te_0)**(-1.5d0)
            deta_dT_ohm   = 0.    
          else
            eta_T_ohm     = eta_ohmic
            deta_dT_ohm   = 0.d0
          end if

          if ( eta_T_dependent .and. (T_or_Te .lt. T_min) ) then
            eta_T     = eta       * (max(T_or_Te,T_min)/T_or_Te_0)**(-1.5d0)
            deta_dT   = 0.d0
            d2eta_d2T = 0.d0

            eta_T_ohm = eta_ohmic * (max(T_or_Te,T_min)/T_or_Te_0)**(-1.5d0)
            deta_dT_ohm   = 0.
          end if

          ! --- Temperature dependent viscosity
          if ( visco_T_dependent ) then
            visco_T     =   visco * (T_or_Te_corr/T_or_Te_0)**(-1.5d0)
            dvisco_dT   = - visco * (1.5d0)  * T_or_Te_corr**(-2.5d0) * T_or_Te_0**(1.5d0)
            d2visco_dT2 =   visco * (3.75d0) * T_or_Te_corr**(-3.5d0) * T_or_Te_0**(1.5d0)
            if (T_or_Te .lt. T_min) then
              visco_T     = visco  * (max(T_or_Te,T_min)/T_or_Te_0)**(-1.5d0)
              dvisco_dT   = 0.d0
              d2visco_dT2 = 0.d0
            endif
          else
            visco_T     = visco
            dvisco_dT   = 0.d0
            d2visco_dT2 = 0.d0
          end if

          psi_norm = get_psi_n( ps0, y_g(ms,mt))
          
          ! --- Hyper-resistivity
          if ( eta_num_psin_dependent ) then
            eta_num_T   = eta_num * 0.5d0 * ( 1.d0 - tanh( (psi_norm-eta_num_prof(1))/eta_num_prof(2)) )      
            deta_num_dT = 0.d0      
          else if ( eta_num_T_dependent ) then
            eta_num_T     =   eta_num   * (T_or_Te_corr/T_or_Te_0)**(-3.d0)
            deta_num_dT   = - eta_num   * (3.d0)  * T_or_Te_corr**(-4.d0) * T_or_Te_0**(3.d0)
            if (T_or_Te .lt. T_min) then
              eta_num_T     = eta_num    * (max(T_or_Te,T_min)/T_or_Te_0)**(-3.d0)
              deta_num_dT   = 0.d0
            endif
          else
            eta_num_T     = eta_num
            deta_num_dT   = 0.d0
          end if

          
          ! --- Temperature dependent hyper-viscosity
          if ( visco_num_T_dependent ) then
            visco_num_T     =   visco_num   * (T_or_Te_corr/T_or_Te_0)**(-3.d0)
            dvisco_num_dT   = - visco_num   * (3.d0)  * T_or_Te_corr**(-4.d0) * T_or_Te_0**(3.d0)
            if (T_or_Te .lt. T_min) then
              visco_num_T     = visco_num    * (max(T_or_Te,T_min)/T_or_Te_0)**(-3.d0)
              dvisco_num_dT   = 0.d0
            endif
          else
            visco_num_T     = visco_num
            dvisco_num_dT   = 0.d0
          end if

          ! --- Diamagnetic viscosity
          if (Wdia) then
            W_dia = + tauIC*2. /r0_corr    * (Pi0_xx + Pi0_x/bigR + Pi0_yy) &
                    - tauIC*2. /r0_corr**2 * (r0_x*Pi0_x + r0_y*Pi0_y)
          else
            W_dia = 0.d0
          endif

          ! --- Bootstrap current 
          if (bootstrap) then
            ! --- Full Sauter formula
            call bootstrap_current(bigR, y_g(ms,mt),                     &
                                   R_axis,   Z_axis,   psi_axis,         &
                                   R_xpoint, Z_xpoint, psi_bnd, psi_norm,&
                                   ps0, ps0_x, ps0_y,                    &
                                   r0,  r0_x,  r0_y,                     &
                                   Ti0, Ti0_x, Ti0_y,                    &
                                   Te0, Te0_x, Te0_y,                  Jb)

            ! --- Full Sauter formula for initial profiles
            call bootstrap_current(bigR, y_g(ms,mt),                       &
                                   R_axis,   Z_axis,   psi_axis,           &
                                   R_xpoint, Z_xpoint, psi_bnd, psi_norm,  &
                                   ps0, ps0_x, ps0_y,                      &
                                   eq_zne(ms,mt),  zn_x,  zn_y,            &
                                   zTi, zTi_x, zTi_y,                      &
                                   zTe, zTe_x, zTe_y,                  Jb_0)

            ! --- Subtract the initial equilibrium part
            Jb = Jb - Jb_0
          else
            Jb = 0.d0
          endif


          D_prof   = get_dperp (psi_norm)
          if ( with_TiTe ) then
            ZKi_prof = get_zk_iperp(psi_norm)
            ZKe_prof = get_zk_eperp(psi_norm)
          else
            ZK_prof = get_zkperp(psi_norm)
          end if

          ! --- Increase diffusivity if very small density/temperature
          if (r0 .lt. D_prof_neg_thresh)  then
            D_prof  = D_prof_neg
          endif

          if ( with_TiTe ) then ! (with_TiTe) ****************************************************
            if (Ti0 .lt. ZK_i_prof_neg_thresh) then
              ZKi_prof = ZK_i_prof_neg
            end if
            if (Ti0 .lt. ZK_i_par_neg_thresh) then
              ZKi_par_T = ZK_i_par_neg
            endif
            if (Te0 .lt. ZK_e_prof_neg_thresh) then
              ZKe_prof = ZK_e_prof_neg
            end if
            if (Te0 .lt. ZK_e_par_neg_thresh) then
              ZKe_par_T = ZK_e_par_neg
            endif
          else ! (with_TiTe = .f.), i.e. with single temperature ***************************************
            if (T0 .lt. ZK_prof_neg_thresh) then
              ZK_prof = ZK_prof_neg
            end if
            if (T0 .lt. ZK_par_neg_thresh) then
              ZK_par_T = ZK_par_neg
            endif
          endif ! (with_TiTe) ********************************************************************


          phi       = 2.d0*PI*float(mp-1)/float(n_plane) / float(n_period)
          delta_phi = 2.d0*PI/float(n_plane) / float(n_period)

          source_pellet = 0.d0
          source_volume = 0.d0

          Dn0x = D_neutral_x      
          Dn0y = D_neutral_y      
          Dn0p = D_neutral_p    

          if (use_pellet) then

            call pellet_source2(pellet_amplitude,pellet_R,pellet_Z,pellet_psi,pellet_phi, &
                                pellet_radius, pellet_delta_psi, pellet_sig, pellet_length, pellet_ellipse, pellet_theta, &
                                x_g(ms,mt),y_g(ms,mt), ps0, phi, r0_corr, Te0_corr, &
                                central_density, pellet_particles, pellet_density, total_pellet_volume, &
                                source_pellet, source_volume)
          endif

          Vt0   = V_source(ms,mt)
          if (normalized_velocity_profile) then
            Vt0_x = dV_dpsi_source(ms,mt)*ps0_x
            Vt0_y = dV_dz_source(ms,mt)+dV_dpsi_source(ms,mt)*ps0_y
          else
            Omega_tor0_x = dV_dpsi_source(ms,mt)*ps0_x
            Omega_tor0_y = dV_dz_source(ms,mt)+dV_dpsi_source(ms,mt)*ps0_y
          endif

          !-------------------------------------------
          ! --- Normalisation of the ionization energy cost for Deuterium
          !------------------------------------------- 
      
          ksiion = central_density * 1.d20 * ksi_ion
 
          if (with_neutrals) then
            !-------------------------------------------------------------------------------
            !---------- Atomic data for neutrals and background constant impurity ----------
            !-------------------------------------------------------------------------------

            call atomic_coeff_deuterium(Te0, Sion_T, dSion_dT, Srec_T, dSrec_dT,        &
                                        LradDcont_T, dLradDcont_dT, LradDrays_T, dLradDrays_dT, r0, rn0, .true. )

            if (.not. with_TiTe) then
              ! --- Transform derivatives on Te to derivatives in total T
              dSion_dT      = dSion_dT      / 2.d0
              dSrec_dT      = dSrec_dT      / 2.d0
              dLradDrays_dT = dLradDrays_dT / 2.d0
              dLradDcont_dT = dLradDcont_dT / 2.d0
            endif
    
            !--------------------------------------------------------
            ! --- Source of neutrals, e.g. from MGI/SPI
            !--------------------------------------------------------
      
            source_neutral       = 0.d0    
            source_neutral_drift = 0.d0 
      
            call total_neutral_source(x_g(ms,mt),y_g(ms,mt),phi,ps0,source_neutral,source_neutral_drift)
          
            source_neutral       = max(source_neutral,0.)
            source_neutral_drift = max(source_neutral_drift,0.)

          else ! no neutrals (neutral terms are always multiplied by one of these coefficients)
            
            Sion_T        = 0.d0
            dSion_dT      = 0.d0 
            Srec_T        = 0.d0
            dSrec_dT      = 0.d0
            LradDcont_T   = 0.d0 
            dLradDcont_dT = 0.d0
            LradDrays_T   = 0.d0 
            dLradDrays_dT = 0.d0

          endif  ! with_neutrals
                 
          !------------------------------------------------------------------------------------------
          ! ---Calculate energy teleported in JOREK unit (sink at R and source at R + drift_distance)
          !------------------------------------------------------------------------------------------
          ! Input energy_teleported is in eV
          power_dens_teleport_ju = 0.d0
          if (with_neutrals .and. energy_teleported /= 0.d0) then
            power_dens_teleport_ju = (-source_neutral + source_neutral_drift)  * energy_teleported * &
                                     EL_CHG * (GAMMA-1) * MU_ZERO * 1.d20 * central_density
          end if

          !-----------------------------------------------------------------
          ! --- Radiation from background impurity, using ADAS (by default)
          !-----------------------------------------------------------------
          ne_SI = r0_corr * 1.d20 * central_density !electron density (SI)
          if (with_TiTe) then 
            Te_corr_eV =       Te0_corr/(EL_CHG*MU_ZERO*central_density*1.d20)  ! Te in eV
            Te_eV =       Te0/(EL_CHG*MU_ZERO*central_density*1.d20)  ! Te in eV, uncorrected
          else
            Te_corr_eV = 0.5d0* T0_corr/(EL_CHG*MU_ZERO*central_density*1.d20)  ! Te in eV
            Te_eV = 0.5d0* T0/(EL_CHG*MU_ZERO*central_density*1.d20)  ! Te in eV, uncorrected
          endif
          if (use_imp_adas) then  ! use open adas by default
            frad_bg = 0. 
            dfrad_bg_dT = 0.
            do i_imp =1, n_adas
              if (i_imp == index_main_imp) cycle
              r_imp = nimp_bg(i_imp) / (1.d20 * central_density)  ! Background impurity density in JU     
              if (ne_SI > ne_SI_min .and. Te_eV > Te_eV_min .and. r_imp > 0) then
                Lrad_imp = 0.0
                dLrad_imp_dT = 0.0
                call radiation_function_linear(imp_adas(i_imp),imp_cor(i_imp),log10(ne_SI),    & 
                                               log10(Te_corr_eV*EL_CHG/K_BOLTZ),.true.,Lrad_imp,dLrad_imp_dT)
                dLrad_imp_dT = dLrad_imp_dT * dT0_corr_dT            
              else     
                Lrad_imp = 0.
                dLrad_imp_dT = 0.
              end if
              if (dLrad_imp_dT/=dLrad_imp_dT) then
                write(*,*) "WARNING: dLrad_imp_dT ", dLrad_imp_dT
                stop
              end if

              frad_bg = frad_bg + r_imp * Lrad_imp
              dfrad_bg_dT =  dfrad_bg_dT + r_imp * dLrad_imp_dT 

            end do

          else 

            if ( trim(imp_type(1)) == 'Ar') then ! Hard-coded fitting exists for argon
              Arad_bg = 2.4d-31
              Brad_bg = 20.
              Crad_bg = 0.8
      
              frad_bg     = (2./3.)*(1./(central_mass*MASS_PROTON))*((MU_ZERO*central_mass*MASS_PROTON*central_density*1.d20)**(1.5d0))            &
                            *nimp_bg(1)*Arad_bg*exp(-((log(Te_corr_eV)-log(Brad_bg))**2.)/Crad_bg**2.)
      
              dfrad_bg_dT = -(1./3.)*((MU_ZERO*central_mass*MASS_PROTON*central_density*1.d20)**(0.5d0))*(1./EL_CHG)                               &
                            *2.*(nimp_bg(1)*Arad_bg/Crad_bg**2.)*(log(Te_corr_eV)-log(Brad_bg))*(1./Te_corr_eV)*exp(-((log(Te_corr_eV)-log(Brad_bg))**2.)/Crad_bg**2.)
            else
              write(*,*) "WARNING: hard-coded fitting doesn't exist for  ", trim(imp_type(1)), ", use open adas instead!"
              stop
            end if 

          end if

         ! For shock capturing stabilization
         tau_sc = 0.d0
         if (use_sc) call calculate_sc_quantities()

!--------------------------------------------------------

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

            Bgrad_Ti          = ( F0 / BigR * Ti0_p +  Ti0_x * ps0_y - Ti0_y * ps0_x ) / BigR
            Bgrad_Te          = ( F0 / BigR * Te0_p +  Te0_x * ps0_y - Te0_y * ps0_x ) / BigR
            Bgrad_T           = ( F0 / BigR * T0_p  +  T0_x  * ps0_y - T0_y  * ps0_x ) / BigR

            BB2              = (F0*F0 + ps0_x * ps0_x + ps0_y * ps0_y )/BigR**2
            Btheta2          = (ps0_x * ps0_x + ps0_y * ps0_y )/BigR**2

            v_ps0_x  = v_xx  * ps0_y - v_xy  * ps0_x + v_x  * ps0_xy - v_y * ps0_xx
            v_ps0_y  = v_xy  * ps0_y - v_yy  * ps0_x + v_x  * ps0_yy - v_y * ps0_xy

            do i_term=1, max_terms_loop

            if (present(get_terms)) then
              factor          = 0.d0
              factor(:,i_term) = 1.d0 / tstep
            else
              factor          = 1.d0
            endif

            !###################################################################################################
            !#  Induction Equation                                                                             #
            !###################################################################################################

            rhs_ij(var_psi) = v * eta_T  * (zj0 - current_source(ms,mt) - Jb)/ BigR           * xjac * tstep * factor(var_psi,1) &
                      + v * (ps0_s * u0_t - ps0_t * u0_s)                                            * tstep * factor(var_psi,2) &
                      - v * F0 / BigR  * u0_p                                                 * xjac * tstep * factor(var_psi,2) &
                      + eta_num_T * (v_x * zj0_x + v_y * zj0_y)                               * xjac * tstep * factor(var_psi,3) &

                      - v * tauIC*2./(r0_corr*BB2) * F0**2/BigR**2 * (ps0_s * Pe0_t - ps0_t * Pe0_s) * tstep * factor(var_psi,4) &
                      + v * tauIC*2./(r0_corr*BB2) * F0**3/BigR**3 * Pe0_p                    * xjac * tstep * factor(var_psi,4) &

                      + zeta * v * delta_g(mp,var_psi,ms,mt) / BigR                           * xjac * factor(var_psi,5)


            !###################################################################################################
            !#  Perpendicular Momentum Equation                                                                #
            !###################################################################################################

            rhs_ij(var_u) =  - 0.5d0 * vv2 * (v_x * r0_y_hat - v_y * r0_x_hat)                              * xjac * tstep * factor(var_u,1) &
                         - r0_hat * BigR**2 * w0 * (v_s * u0_t - v_t * u0_s)                                       * tstep * factor(var_u,1) &
                         + v * (ps0_s * zj0_t - ps0_t * zj0_s )                                                    * tstep * factor(var_u,2) &

                         - visco_T * BigR * (v_x * w0_x + v_y * w0_y)                                       * xjac * tstep * factor(var_u,3) &

                         - v * F0 / BigR * zj0_p                                                            * xjac * tstep * factor(var_u,2) &
                         + BigR**2 * (v_s * p0_t - v_t * p0_s)                                                     * tstep * factor(var_u,4) &

                         - visco_num_T * (v_xx + v_x/Bigr + v_yy)*(w0_xx + w0_x/Bigr + w0_yy)               * xjac * tstep * factor(var_u,5) &

                         - tgnum_u * 0.25d0 * r0_hat * BigR**3 * (w0_x * u0_y - w0_y * u0_x) &
                                   * ( v_x * u0_y - v_y * u0_x) * xjac * tstep * tstep * factor(var_u,6)    &
!====================================New TG_num terms=================================
                         - tgnum_u * 0.25d0 * w0 * BigR**3 * (r0_x_hat * u0_y - r0_y_hat * u0_x) &
                                   * ( v_x * u0_y - v_y * u0_x) * xjac * tstep * tstep * fact_conservative_u * factor(var_u,6) &
!===============================End of NewTG_num terms==============================

                         - v * tauIC*2. * BigR**4 * (Pi0_s * w0_t - Pi0_t * w0_s)                                  * tstep * factor(var_u,7) &

                         - tauIC*2. * BigR**3 * Pi0_y * (v_x* u0_x + v_y * u0_y)                            * xjac * tstep * factor(var_u,7) &

                         - v * tauIC*2. * BigR**4 * (u0_xy * (Pi0_xx - Pi0_yy) - Pi0_xy * (u0_xx - u0_yy) ) * xjac * tstep * factor(var_u,7) &

                         ! --- Diamagnetic viscosity
                         + dvisco_dT * bigR * W_dia * (v_x*Ti0_x + v_y*Ti0_y)                               * xjac * tstep * factor(var_u,8) &
                         + visco_T   * bigR * W_dia * (v_xx + v_x/bigR + v_yy)                              * xjac * tstep * factor(var_u,8) &

                        - zeta * BigR * r0_hat * (v_x * delta_u_x + v_y * delta_u_y) * xjac * factor(var_u,9)      &

                         ! Not to be included in conservative form
                         + BigR**3 * (particle_source(ms,mt) + source_pellet) * (v_x * u0_x + v_y * u0_y) * xjac* tstep* factor(var_u,10)  &
                                   * (1.d0 - fact_conservative_u)  &          
 
                         ! New terms coming from -(\partial_t \rho + \nabla \cdot (\rho \mathbf{v})) \mathbf{v} in RHS of momentum equation
                         ! (see wiki: https://www.jorek.eu/wiki/doku.php?id=model500_501_555#equations):		      
                         + fact_conservative_u * (                                                                           &
                             - zeta * BigR * BigR**2 * delta_g(mp,var_rho,ms,mt) * (v_x * u0_x + v_y * u0_y)  * xjac         &
                             - BigR**2 * (r0_x_hat * u0_y - r0_y_hat * u0_x) * (v_x * u0_x + v_y * u0_y)      * xjac * tstep &
                             + BigR * F0 * (r0 * vpar0_p + vpar0 * r0_p) * (v_x * u0_x + v_y * u0_y)          * xjac * tstep &
                             + BigR**2 * r0 * (vpar0_x * ps0_y - vpar0_y * ps0_x) * (v_x * u0_x + v_y * u0_y) * xjac * tstep &
                             + BigR**2 * vpar0 * (r0_x * ps0_y - r0_y * ps0_x)    * (v_x * u0_x + v_y * u0_y) * xjac * tstep &
                           ) * factor(var_u,10)   
            
            !------------------------------------------------------------------------ NEO
            if (NEO) then
              rhs_ij(var_u) = rhs_ij(var_u)  + amu_neo_prof(ms,mt)*BB2/((Btheta2+epsil)**2)     &
                        * (ps0_x * v_x + ps0_y * v_y) *                                         &  
                          ( r0 * (ps0_x * u0_x + ps0_y * u0_y)                                  &
                            + tauIC*2. * (ps0_x * Pi0_x + ps0_y * Pi0_y)                        &
                            + aki_neo_prof(ms,mt) * tauIC*2. * r0 * (ps0_x*Ti0_x + ps0_y*Ti0_y) &
                            - r0 * Vpar0 * Btheta2)     * xjac * tstep * BigR * factor(var_u,11)
            endif
            !------------------------------------------------------------------------ NEO

            !###################################################################################################
            !#  Current Definition Equation                                                                    #
            !###################################################################################################

            rhs_ij(var_zj) = - ( v_x * ps0_x  + v_y * ps0_y + v*zj0) / BigR * xjac * factor(var_zj,1)

            !###################################################################################################
            !#  Vorticity Definition Equation                                                                  #
            !###################################################################################################

            rhs_ij(var_w) = - ( v_x * u0_x   + v_y * u0_y  + v*w0)  * BigR * xjac * factor(var_w,1)

            !###################################################################################################
            !#  Density Equation                                                                               #
            !###################################################################################################

            rhs_ij(var_rho)  = v * BigR * (particle_source(ms,mt) + source_pellet)                * xjac * tstep * factor(var_rho,1) &
                       + v * BigR**2 * ( r0_s * u0_t - r0_t * u0_s)                                      * tstep * factor(var_rho,2) &
                       + v * 2.d0 * BigR * r0 * u0_y                                              * xjac * tstep * factor(var_rho,3) &

                       - ((D_par+D_par_sc_num*tau_sc) - D_prof)  * BigR / BB2 * Bgrad_rho_star * Bgrad_rho * xjac * tstep * factor(var_rho,4) &
                       - D_prof * BigR  * (v_x*r0_x + v_y*r0_y                                  ) * xjac * tstep * factor(var_rho,5) &

                       - v * F0 / BigR * Vpar0 * r0_p                                             * xjac * tstep * factor(var_rho,6) &
                       - v * Vpar0 * (r0_s * ps0_t - r0_t * ps0_s)                                       * tstep * factor(var_rho,6) &
                       - v * F0 / BigR * r0 * vpar0_p                                             * xjac * tstep * factor(var_rho,3) &
                       - v * r0 * (vpar0_s * ps0_t - vpar0_t * ps0_s)                                    * tstep * factor(var_rho,3) &

                       + v * 2.d0 * tauIC*2. * Pi0_y * BigR                                       * xjac * tstep * factor(var_rho,7) &

                       + v * r0 * rn0 * BigR * Sion_T                                             * xjac * tstep * factor(var_rho,8) &
                       - v * r0 * r0  * BigR * Srec_T                                             * xjac * tstep * factor(var_rho,9) &
                       
                       + zeta * v * delta_g(mp,var_rho,ms,mt) * BigR                              * xjac         * factor(var_rho,10)&

                       - D_perp_num * (v_xx + v_x/Bigr + v_yy)*(r0_xx + r0_x/Bigr + r0_yy) * BigR * xjac * tstep * factor(var_rho,11)&

                       - tgnum_rho * 0.25d0 * BigR**3 * (r0_x * u0_y - r0_y * u0_x)                              &
                                                    * ( v_x * u0_y - v_y * u0_x) * xjac * tstep * tstep          * factor(var_rho,12)&
                       - tgnum_rho * 0.25d0 / BigR * vpar0**2                                                    &
                                 * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                              &
                                 * ( v_x * ps0_y -  v_y * ps0_x                   ) * xjac * tstep * tstep       * factor(var_rho,12)

            rhs_ij_k(var_rho) = - ((D_par+D_par_sc_num*tau_sc)-D_prof) * BigR / BB2 * Bgrad_rho_k_star * Bgrad_rho * xjac * tstep * factor(var_rho,4) &
                          - D_prof * BigR  * (                  v_p*r0_p /BigR**2 )               * xjac * tstep * factor(var_rho,5) &

                       - tgnum_rho * 0.25d0 / BigR * vpar0**2 &
                                 * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                              &
                                 * (                            + F0 / BigR * v_p) * xjac * tstep * tstep        * factor(var_rho,12)

            !###################################################################################################
            !#  Parallel Velocity Equation                                                                     #
            !###################################################################################################

            if ( with_vpar ) then
              rhs_ij(var_vpar) = - v * F0 / BigR * P0_p                                                                 * xjac * tstep * factor(var_vpar,1) &
                                 - v * (P0_s * ps0_t - P0_t * ps0_s)                                                           * tstep * factor(var_vpar,1) &
                                                                                                                        
                                ! Not to be included in the conservative form
                                 - v*(particle_source(ms,mt) + source_pellet) * vpar0 * BB2   * BigR                    * xjac * tstep * factor(var_vpar,2) &
                                     * (1.d0 - fact_conservative_u)  &          
                                                                                                                       
                                 - 0.5d0 * r0 * vpar0**2 * BB2 * (ps0_s * v_t - ps0_t * v_s)                                   * tstep * factor(var_vpar,3) &
                                 - 0.5d0 * v  * vpar0**2 * BB2 * (ps0_s * r0_t - ps0_t * r0_s)                                 * tstep * factor(var_vpar,3) &
                                 + 0.5d0 * v  * vpar0**2 * BB2 * F0 / BigR * r0_p                                       * xjac * tstep * factor(var_vpar,3) &
                                 
                                 - visco_par_num * (v_xx + v_x/Bigr + v_yy)*(vpar0_xx + vpar0_x/Bigr + vpar0_yy) * BigR * xjac * tstep * factor(var_vpar,4) &
                                 
                                 + zeta * v * delta_g(mp,var_vpar,ms,mt) * r0_corr * F0**2 / BigR                       * xjac         * factor(var_vpar,5) &
                                 + zeta * v * r0_corr * vpar0 * (ps0_x * delta_ps_x + ps0_y * delta_ps_y) / BigR        * xjac         * factor(var_vpar,5) &

                                 ! New terms coming from -(\partial_t \rho + \nabla \cdot (\rho \mathbf{v})) \mathbf{v} in RHS of momentum equation
                                 ! (see wiki: https://www.jorek.eu/wiki/doku.php?id=model500_501_555#equations):
                                 + fact_conservative_u * (                                                 &
                                     + zeta * v * delta_g(mp,var_rho,ms,mt) * vpar0 * F0**2 / BigR * xjac  &
                                     + v * (r0_x_hat * u0_y - r0_y_hat * u0_x)       * vpar0 * BB2 * xjac * tstep  &
                                     - v * F0 / BigR * (r0 * vpar0_p + r0_p * vpar0) * vpar0 * BB2 * xjac * tstep  &
                                     - v * r0 * (vpar0_x * ps0_y - vpar0_y * ps0_x)  * vpar0 * BB2 * xjac * tstep  &
                                     - v * vpar0 * (r0_x * ps0_y - r0_y * ps0_x)     * vpar0 * BB2 * xjac * tstep  &
                                                          ) * factor(var_vpar,2)&

                           - tgnum_vpar * 0.25d0 * r0 * Vpar0**2 * BB2                                                                 &
                                     * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac + F0 / BigR * vpar0_p) / BigR                        &
                                     * (-(ps0_s * v_t     - ps0_t * v_s)    /xjac                      )        * xjac * tstep * tstep * factor(var_vpar,6)&
                           - tgnum_vpar * 0.25d0 * v  * Vpar0**2 * BB2 * (1.d0 - fact_conservative_u)                                  &
                                     * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac + F0 / BigR * vpar0_p) / BigR                        &
                                     * (-(ps0_s * r0_t    - ps0_t * r0_s)   /xjac + F0 / BigR * r0_p)           * xjac * tstep * tstep * factor(var_vpar,6)&
!=============================== New TG_num terms==================================
                           - tgnum_vpar * 0.25d0 * vpar0 * Vpar0**2 * BB2 * fact_conservative_u            &
                                     * (-(ps0_s * r0_t - ps0_t * r0_s)/xjac + F0 / BigR * r0_p) / BigR     &
                                     * (-(ps0_s * v_t     - ps0_t * v_s)    /xjac)  * xjac * tstep * tstep * factor(var_vpar,6)&
!===============================End of new TG_num terms============================
                 
                                 + (1.d0 - delta_n_convection) * (     &
                                   - v *(r0 * rn0 * Sion_T) * vpar0 * BB2 * BigR                                * xjac * tstep * factor(var_vpar,7)&
                                   + v *(r0 * r0  * Srec_T) * vpar0 * BB2 * BigR                                * xjac * tstep * factor(var_vpar,8)&
                                   ) 
                    
              if (normalized_velocity_profile) then
                rhs_ij(var_vpar) = rhs_ij(var_vpar) - (visco_par + visco_par_sc_num * tau_sc) * (v_x * (vpar0_x-Vt0_x) + v_y * (vpar0_y-Vt0_y)) * BigR* xjac * tstep * factor(var_vpar,9 ) 
              else
                rhs_ij(var_vpar) = rhs_ij(var_vpar) - (visco_par + visco_par_sc_num * tau_sc) * (v_x * (vpar0_x * F0**2 / BigR**2 - 2.d0 * vpar0 * F0**2 / BigR**3  - 2.d0 * PI * F0 * Omega_tor0_x ) & 
                                                   + v_y * (vpar0_y * F0**2 / BigR**2 - 2.d0 * PI * F0 * Omega_tor0_y) ) * BigR* xjac * tstep  * factor(var_vpar,9 )
              endif

              if (NEO) then
                rhs_ij(var_vpar) =  rhs_ij(var_vpar)  + amu_neo_prof(ms,mt)*BB2/(Btheta2+epsil)  &
                          * v * ( r0 * (ps0_x * u0_x + ps0_y * u0_y)                             &
                                + tauIC*2.   * (ps0_x * Pi0_x + ps0_y * Pi0_y)                   &
                                + aki_neo_prof(ms,mt) * tauIC*2. * r0 * (ps0_x * Ti0_x + ps0_y * Ti0_y) - r0 * Vpar0 * Btheta2) * xjac * tstep * factor(var_vpar,10) * BigR
              endif
   
              rhs_ij_k(var_vpar) = + 0.5d0 * r0 * vpar0**2 * BB2 * F0 / BigR * v_p                 * xjac * tstep * factor(var_vpar,3) &
  
                 - tgnum_vpar * 0.25d0 * r0 * Vpar0**2 * BB2 &
                           * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac + F0 / BigR * vpar0_p) / BigR  &
                           * (                                          + F0 / BigR * v_p)  * xjac * tstep * tstep * factor(var_vpar,6)& 
!=============================== New TG_num terms==================================

                 - tgnum_vpar * 0.25d0 * vpar0 * Vpar0**2 * BB2 * fact_conservative_u &
                           * (-(ps0_s * r0_t - ps0_t * r0_s)/xjac + F0 / BigR * r0_p) / BigR  &
                           * (                                          + F0 / BigR * v_p)  * xjac * tstep * tstep * factor(var_vpar,6)

!===============================End of new TG_num terms============================
            end if ! (with_vpar)
            

            if ( with_TiTe ) then ! (with_TiTe) ****************************************************
              
              !###################################################################################################
              !#  Ion Energy Equation                                                                            #
              !###################################################################################################
  
              rhs_ij(var_Ti) =  v * BigR * heat_source_i(ms,mt)                                  * xjac * tstep * factor(var_Ti,1) &
              
                              + v * r0 * BigR**2 * ( Ti0_s * u0_t - Ti0_t * u0_s)                       * tstep * factor(var_Ti,2) &
                              + v * Ti0 * BigR**2 * ( r0_s * u0_t - r0_t * u0_s)                        * tstep * factor(var_Ti,2) &
                              
                              + v * r0 * Ti0 * 2.d0* GAMMA * BigR * u0_y                         * xjac * tstep * factor(var_Ti,3) &
                              
                              - v * r0 * F0 / BigR * Vpar0 * Ti0_p                               * xjac * tstep * factor(var_Ti,4) &
                              - v * Ti0 * F0 / BigR * Vpar0 * r0_p                               * xjac * tstep * factor(var_Ti,4) &
                              
                              - v * r0 * Vpar0 * (Ti0_s * ps0_t - Ti0_t * ps0_s)                        * tstep * factor(var_Ti,4) &
                              - v * Ti0 * Vpar0 * (r0_s * ps0_t - r0_t * ps0_s)                         * tstep * factor(var_Ti,4) &
                              
                              - v * r0 * Ti0 * GAMMA * (vpar0_s * ps0_t - vpar0_t * ps0_s)              * tstep * factor(var_Ti,3) &
                              - v * r0 * Ti0 * GAMMA * F0 / BigR * vpar0_p                       * xjac * tstep * factor(var_Ti,3) &
                              
                              - (ZKi_par_T-ZKi_prof) * BigR / BB2 * Bgrad_T_star * Bgrad_Ti      * xjac * tstep * factor(var_Ti,5) &
                              - ZKi_prof * BigR * (v_x*Ti0_x + v_y*Ti0_y                   )     * xjac * tstep * factor(var_Ti,6) &
                              
                              - ZK_i_perp_num  *  (v_xx + v_x/Bigr + v_yy)*(Ti0_xx + Ti0_x/Bigr + Ti0_yy) * BigR * xjac * tstep * factor(var_Ti,7) &
  
                         - tgnum_Ti* 0.25d0 * BigR**3 * Ti0 * (r0_x * u0_y - r0_y * u0_x)                         &
                                            * ( v_x * u0_y - v_y * u0_x)                  * xjac * tstep * tstep  * factor(var_Ti,8)&
                         - tgnum_Ti* 0.25d0 * BigR**3 * r0 * (Ti0_x * u0_y - Ti0_y * u0_x)                        &
                                            * ( v_x * u0_y - v_y * u0_x)                  * xjac * tstep * tstep  * factor(var_Ti,8)&
  
                         - tgnum_Ti* 0.25d0 / BigR * vpar0**2                                                     &
                                   * Ti0 * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                       &
                                   * ( v_x * ps0_y -  v_y * ps0_x                        ) * xjac * tstep * tstep * factor(var_Ti,8)&
                         - tgnum_Ti* 0.25d0 / BigR * vpar0**2                                                     &
                                   * r0 * (Ti0_x * ps0_y - Ti0_y * ps0_x + F0 / BigR * Ti0_p)                     &
                                   * ( v_x * ps0_y -  v_y * ps0_x                        ) * xjac * tstep * tstep * factor(var_Ti,8)&

                            !===================== Additional terms from friction terms============
                            + v * BigR * ((GAMMA - 1.)/2.) * vpar0**2 * BB2 * (r0*rn0*Sion_T)      * xjac * tstep * factor(var_Ti,10) &
                            + v * BigR * ((GAMMA - 1.)/2.) * vv2 * ((r0*rn0*Sion_T))               * xjac * tstep * factor(var_Ti,10) &
                            !==============================End of friction terms=================
                            !============================Behold, the parallel viscous heating terms!=============
                            + (GAMMA - 1.) * v * BigR * visco_par_heating * (vpar0_x * vpar0_x + vpar0_y * vpar0_y) * xjac * tstep * factor(var_Ti,12) &
                            + (GAMMA - 1.) * vpar0 * BigR * visco_par_heating * (v_x * vpar0_x     + v_y * vpar0_y) * xjac * tstep * factor(var_Ti,12) &
                            !==========================End of viscous heating terms==============================

                         + zeta * v * r0_corr  * delta_g(mp,var_Ti,ms,mt) * BigR                 * xjac         * factor(var_Ti,9)   &
                         + zeta * v * Ti0_corr * delta_g(mp,var_rho,ms,mt) * BigR                * xjac         * factor(var_Ti,9)   &
                         ! Energy exchange term
                         + v * BigR * dTi_e                                                      * xjac * tstep * factor(var_Ti,11)                       
  
              rhs_ij_k(var_Ti) = - (ZKi_par_T-ZKi_prof) * BigR / BB2 * Bgrad_T_k_star * Bgrad_Ti * xjac * tstep * factor(var_Ti,5) &
                                 - ZKi_prof * BigR * (                + v_p*Ti0_p /BigR**2 )     * xjac * tstep * factor(var_Ti,6) &
  
                           - tgnum_Ti* 0.25d0 / BigR * vpar0**2                                                     &
                                   * Ti0 * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                         &
                                   * (                                  + F0 / BigR * v_p)  * xjac * tstep * tstep  * factor(var_Ti,8)&
                           - tgnum_Ti* 0.25d0 / BigR * vpar0**2                                                     &
                                   * r0 * (Ti0_x * ps0_y - Ti0_y * ps0_x + F0 / BigR * Ti0_p)                       &
                                   * (                                   + F0 / BigR * v_p) * xjac * tstep * tstep  * factor(var_Ti,8)
  
              !###################################################################################################
              !#  Electron Energy Equation                                                                       #
              !###################################################################################################
  
              rhs_ij(var_Te) =  v * BigR * heat_source_e(ms,mt)                                  * xjac * tstep * factor(var_Te,1 ) &
                                                                                                 
                              + v * r0 * BigR**2  * (Te0_s * u0_t - Te0_t * u0_s)                       * tstep * factor(var_Te,2 ) &
                              + v * Te0 * BigR**2 * ( r0_s * u0_t -  r0_t * u0_s)                       * tstep * factor(var_Te,2 ) &
                                                                                                 
                              + v * r0 * Te0 * 2.d0* GAMMA * BigR * u0_y                         * xjac * tstep * factor(var_Te,3 ) &
                                                                                                 
                              - v * r0 * F0 / BigR * Vpar0 * Te0_p                               * xjac * tstep * factor(var_Te,4 ) &
                              - v * Te0 * F0 / BigR * Vpar0 * r0_p                               * xjac * tstep * factor(var_Te,4 ) &
                                                                                                 
                              - v * r0 * Vpar0 * (Te0_s * ps0_t - Te0_t * ps0_s)                        * tstep * factor(var_Te,4 ) &
                              - v * Te0 * Vpar0 * (r0_s * ps0_t - r0_t * ps0_s)                         * tstep * factor(var_Te,4 ) &
                                                                                                 
                              - v * r0 * Te0 * GAMMA * (vpar0_s * ps0_t - vpar0_t * ps0_s)              * tstep * factor(var_Te,3 ) &
                              - v * r0 * Te0 * GAMMA * F0 / BigR * vpar0_p                       * xjac * tstep * factor(var_Te,3 ) &
                                                                                                 
                              - (ZKe_par_T-ZKe_prof) * BigR / BB2 * Bgrad_T_star * Bgrad_Te      * xjac * tstep * factor(var_Te,5 ) &
                              - ZKe_prof * BigR * (v_x*Te0_x + v_y*Te0_y                   )     * xjac * tstep * factor(var_Te,6 ) &
                              
                              - ZK_e_perp_num  *  (v_xx + v_x/Bigr + v_yy)*(Te0_xx + Te0_x/Bigr + Te0_yy) * BigR * xjac * tstep * factor(var_Te,7 ) &

                              + v * (gamma-1.d0) * eta_T_ohm * (zj0 / BigR)**2.d0         * BigR * xjac * tstep * factor(var_Te,9 ) &

                              - v * BigR * ksiion  * r0 * rn0 * Sion_T                           * xjac * tstep * factor(var_Te,12) &
                              - v * BigR * r0_corr * rn0_corr * LradDrays_T                      * xjac * tstep * factor(var_Te,13) &
                              - v * BigR * r0_corr * r0_corr  * LradDcont_T                      * xjac * tstep * factor(var_Te,14) &
                              - v * BigR * r0_corr * frad_bg                                     * xjac * tstep * factor(var_Te,15) &
  
                              ! Additional energy teleportation term for plasmoid drift
                              + v * BigR * power_dens_teleport_ju                                * xjac * tstep * factor(var_Te,16) &  
                          
                         - tgnum_Te * 0.25d0 * BigR**3 * Te0 * (r0_x * u0_y - r0_y * u0_x)                         &
                                            * ( v_x * u0_y - v_y * u0_x)                   * xjac * tstep * tstep * factor(var_Te,8 )&
                         - tgnum_Te * 0.25d0 * BigR**3 * r0 * (Te0_x * u0_y - Te0_y * u0_x)                       &
                                            * ( v_x * u0_y - v_y * u0_x)                   * xjac * tstep * tstep * factor(var_Te,8 )&
  
                         - tgnum_Te * 0.25d0 / BigR * vpar0**2                                                    &
                                   * Te0 * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                       &
                                   * ( v_x * ps0_y -  v_y * ps0_x                        ) * xjac * tstep * tstep * factor(var_Te,8 )&
                         - tgnum_Te * 0.25d0 / BigR * vpar0**2                                                    &
                                   * r0 * (Te0_x * ps0_y - Te0_y * ps0_x + F0 / BigR * Te0_p)                     &
                                   * ( v_x * ps0_y -  v_y * ps0_x                        ) * xjac * tstep * tstep * factor(var_Te,8 )&
  
                         + zeta * v * r0_corr  * delta_g(mp,var_Te,ms,mt)  * BigR                * xjac         * factor(var_Te,10)  &
                         + zeta * v * Te0_corr * delta_g(mp,var_rho,ms,mt) * BigR                * xjac         * factor(var_Te,10)  &
                         ! Energy exchange term
                         + v * BigR * dTe_i                                                      * xjac * tstep * factor(var_Te,11)
  
              rhs_ij_k(var_Te) = - (ZKe_par_T-ZKe_prof) * BigR / BB2 * Bgrad_T_k_star * Bgrad_Te * xjac * tstep * factor(var_Te,5) &
                                 - ZKe_prof * BigR * (                + v_p*Te0_p /BigR**2 )     * xjac * tstep * factor(var_Te,6) &
  
                           - tgnum_Te * 0.25d0 / BigR * vpar0**2                                                    &
                                   * Te0 * (r0_x * ps0_y - r0_y * ps0_x  + F0 / BigR * r0_p)                        &
                                   * (                                   + F0 / BigR * v_p) * xjac * tstep * tstep * factor(var_Te,8 ) &
                           - tgnum_Te * 0.25d0 / BigR * vpar0**2                                                    &
                                   * r0 * (Te0_x * ps0_y - Te0_y * ps0_x + F0 / BigR * Te0_p)                       &
                                   * (                                   + F0 / BigR * v_p) * xjac * tstep * tstep * factor(var_Te,8 )

            else ! (with_TiTe = .f.), i.e. with single temperature ***************************************
  
              !###################################################################################################
              !#  Electron + Ion Energy Equation                                                                 #
              !###################################################################################################
  
              rhs_ij(var_T) =  v * BigR * heat_source(ms,mt)                                    * xjac * tstep * factor(var_T,1 ) &
              
                             + v * r0 * BigR**2 * (T0_s  * u0_t - T0_t * u0_s)                         * tstep * factor(var_T,2 ) &
                             + v * T0 * BigR**2 * ( r0_s * u0_t - r0_t * u0_s)                         * tstep * factor(var_T,2 ) &
                                                                                            
                             + v * r0 * T0 * 2.d0* GAMMA * BigR * u0_y                          * xjac * tstep * factor(var_T,3 ) &
                            
                             - v * r0 * F0 / BigR * Vpar0 * T0_p                                * xjac * tstep * factor(var_T,4 ) &
                             - v * T0 * F0 / BigR * Vpar0 * r0_p                                * xjac * tstep * factor(var_T,4 ) &
                            
                             - v * r0 * Vpar0 * (T0_s * ps0_t - T0_t * ps0_s)                          * tstep * factor(var_T,4 ) &
                             - v * T0 * Vpar0 * (r0_s * ps0_t - r0_t * ps0_s)                          * tstep * factor(var_T,4 ) &
                                                                                             
                             - v * r0 * T0 * GAMMA * (vpar0_s * ps0_t - vpar0_t * ps0_s)               * tstep * factor(var_T,3 ) &
                             - v * r0 * T0 * GAMMA * F0 / BigR * vpar0_p                        * xjac * tstep * factor(var_T,3 ) &
                            
                             - (ZK_par_T-ZK_prof) * BigR / BB2 * Bgrad_T_star * Bgrad_T         * xjac * tstep * factor(var_T,5 ) &
                             - ZK_prof * BigR * (v_x*T0_x + v_y*T0_y                   )        * xjac * tstep * factor(var_T,6 ) &
                            
                             - ZK_perp_num  *  (v_xx + v_x/Bigr + v_yy)*(T0_xx + T0_x/Bigr + T0_yy) * BigR * xjac * tstep * factor(var_T,7 ) &
                            
                             + v * (gamma-1.d0) * eta_T_ohm * (zj0 / BigR)**2.d0         * BigR * xjac * tstep * factor(var_T,9 ) &

                             !===================== Additional terms from friction terms============
                             + v * BigR * ((GAMMA - 1.)/2.) * vpar0**2 * BB2 * (r0*rn0*Sion_T)  * xjac * tstep * factor(var_T,11) &
                             + v * BigR * ((GAMMA - 1.)/2.) * vv2 * ((r0*rn0*Sion_T))           * xjac * tstep * factor(var_T,11) &
                             !==============================End of friction terms=================
                             !============================Behold, the parallel viscous heating terms!=============
                             + (GAMMA - 1.) * v * BigR * visco_par_heating * (vpar0_x * vpar0_x + vpar0_y * vpar0_y) * xjac * tstep * factor(var_T,17) &
                             + (GAMMA - 1.) * vpar0 * BigR * visco_par_heating * (v_x * vpar0_x     + v_y * vpar0_y) * xjac * tstep * factor(var_T,17) &
                             !==========================End of viscous heating terms==============================

                             - v * BigR * ksiion  * r0 * rn0 * Sion_T                           * xjac * tstep * factor(var_T,12) &
                             - v * BigR * r0_corr * rn0_corr * LradDrays_T                      * xjac * tstep * factor(var_T,13) &
                             - v * BigR * r0_corr * r0_corr  * LradDcont_T                      * xjac * tstep * factor(var_T,14) &
                             - v * BigR * r0_corr * frad_bg                                     * xjac * tstep * factor(var_T,15) &
                             
                             ! Additional energy teleportation term for plasmoid drift
                             + v * BigR * power_dens_teleport_ju                                * xjac * tstep * factor(var_T,16) &

                            
                             - tgnum_T * 0.25d0 * BigR**3 * T0 * (r0_x * u0_y - r0_y * u0_x)                           &
                                                * ( v_x * u0_y - v_y * u0_x)                    * xjac * tstep * tstep * factor(var_T,8 )&
                             - tgnum_T * 0.25d0 * BigR**3 * r0 * (T0_x * u0_y - T0_y * u0_x)                           &
                                                * ( v_x * u0_y - v_y * u0_x)                    * xjac * tstep * tstep * factor(var_T,8 )&
                            
                             - tgnum_T * 0.25d0 / BigR * vpar0**2                                                      &
                                       * T0 * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                         &
                                       * ( v_x * ps0_y -  v_y * ps0_x                        )  * xjac * tstep * tstep * factor(var_T,8 )&
                             - tgnum_T * 0.25d0 / BigR * vpar0**2                                                      &
                                       * r0 * (T0_x * ps0_y - T0_y * ps0_x + F0 / BigR * T0_p)                         &
                                       * ( v_x * ps0_y -  v_y * ps0_x                        )  * xjac * tstep * tstep * factor(var_T,8 )&
                            
                             + zeta * v * r0_corr * delta_g(mp,var_T,ms,mt)   * BigR            * xjac * factor(var_T,10) &
                             + zeta * v * T0_corr * delta_g(mp,var_rho,ms,mt) * BigR            * xjac * factor(var_T,10)
  
              rhs_ij_k(var_T) = - (ZK_par_T-ZK_prof) * BigR / BB2 * Bgrad_T_k_star * Bgrad_T    * xjac * tstep * factor(var_T,5) &
                                - ZK_prof * BigR * (                + v_p*T0_p /BigR**2 )       * xjac * tstep * factor(var_T,6) &
  
                             - tgnum_T * 0.25d0 / BigR * vpar0**2                                                      &
                                     * T0 * (r0_x * ps0_y - r0_y * ps0_x  + F0 / BigR * r0_p)                          &
                                     * (                                   + F0 / BigR * v_p)   * xjac * tstep * tstep * factor(var_T,8 )&
                             - tgnum_T * 0.25d0 / BigR * vpar0**2                                                      &
                                     * r0 * (T0_x * ps0_y - T0_y * ps0_x + F0 / BigR * T0_p)                           &
                                     * (                                   + F0 / BigR * v_p)   * xjac * tstep * tstep * factor(var_T,8 )

            end if ! (with_TiTe) *******************************************************************

            !###################################################################################################
            !#  Neutral density equation                                                                       #
            !###################################################################################################

            if (with_neutrals) then
             
              rhs_ij(var_rhon) = BigR * (- Dn0x * rn0_x * v_x - Dn0y * rn0_y * v_y )   * xjac * tstep * factor(var_rhon,1)  &         
                      
                        + delta_n_convection*(                                                              &
                        + v * BigR**2 * ( rn0_s * u0_t - rn0_t * u0_s)                              * tstep &
                        + v * 2.d0 * BigR * rn0 * u0_y                                       * xjac * tstep &
                        - v * rn0 * (vpar0_s * ps0_t - vpar0_t * ps0_s)                      * tstep &
                        - v * Vpar0 * (rn0_s * ps0_t - rn0_t * ps0_s)                        * tstep &
                        - v * F0 / BigR * Vpar0 * rn0_p                                      * xjac * tstep &
                        - v * F0 / BigR * rn0 * vpar0_p                                      * xjac * tstep &
                        ) * factor(var_rhon,2)                                                              &

                    - BigR * v * r0 * rn0 * Sion_T                                           * xjac * tstep * factor(var_rhon,3)&  
                    + BigR * v * r0 * r0  * Srec_T                                           * xjac * tstep * factor(var_rhon,4)&
                    + BigR * v * source_neutral_drift                                        * xjac * tstep * factor(var_rhon,5)&
                    - Dn_perp_num * (v_xx + v_x/Bigr + v_yy)*(rn0_xx + rn0_x/Bigr + rn0_yy)  * BigR * xjac * tstep * factor(var_rhon,6)&

                    + v * delta_g(mp,var_rhon,ms,mt) * BigR * xjac * zeta * factor(var_rhon,7)

              rhs_ij_k(var_rhon) = BigR * ( - Dn0p * rn0_p * v_p/BigR**2)   * xjac * tstep * factor(var_rhon,1)                             

            endif ! with_neutrals

            !###################################################################################################
            !#  RHS equations end                                                                              #
            !###################################################################################################

            if (use_fft) then
              index_ij = n_var*(n_order+1)*(i-1) +       n_var*(j-1) + 1
            else
              index_ij = n_tor_local*n_var*(n_order+1)*(i-1) + n_tor_local*n_var*(j-1) + im - n_tor_start +1 
            endif


            ! --- Fill up the rhs
            if (use_fft) then

              do ij = 1, n_var
                RHS_p(mp,index_ij+ij-1) = RHS_p(mp,index_ij+ij-1) + rhs_ij(ij)   * wst
                RHS_k(mp,index_ij+ij-1) = RHS_k(mp,index_ij+ij-1) + rhs_ij_k(ij) * wst
              enddo

              !--- THIS IS ONLY USED FOR DIAGNOSTIC PURPOSES-------------------------
              !--- ELM structures are re-used to plot separate terms in vtk ---------
              if ( present(get_terms) ) then
                do ij = 1, n_var
                  ELM_p(mp,i_term,index_ij+ij-1) = ELM_p(mp,i_term,index_ij+ij-1) + rhs_ij(ij)   * wst
                  ELM_k(mp,i_term,index_ij+ij-1) = ELM_k(mp,i_term,index_ij+ij-1) + rhs_ij_k(ij) * wst
                enddo
              endif 

            else

              do ij = 1, n_var
                RHS(index_ij+(ij-1)*(n_tor_local)) = RHS(index_ij+(ij-1)*(n_tor_local)) &
                                                   + (rhs_ij(ij) + rhs_ij_k(ij)) * wst
              enddo
              !--- THIS IS ONLY USED FOR DIAGNOSTIC PURPOSES-------------------------
              !--- ELM structures are re-used to plot separate terms in vtk ---------
              if ( present(get_terms) ) then
                do ij = 1, n_var
                  ELM(i_term,index_ij+(ij-1)*(n_tor_local)) = ELM(i_term,index_ij+(ij-1)*(n_tor_local)) &
                                                            + (rhs_ij(ij) + rhs_ij_k(ij)) * wst
                enddo
              endif 

            endif ! --- FFT

            enddo ! max_terms_loop for diagnostic purposes

            if ( present(get_terms) ) cycle

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

                  u    = psi    ;  zj    = psi    ;  w    = psi    ; rho    = psi    ;  Ti    = psi    ; vpar    = psi   ; Te   = psi    ; T   = psi
                  u_x  = psi_x  ;  zj_x  = psi_x  ;  w_x  = psi_x  ; rho_x  = psi_x  ;  Ti_x  = psi_x  ; vpar_x  = psi_x ; Te_x = psi_x  ; T_x = psi_x
                  u_y  = psi_y  ;  zj_y  = psi_y  ;  w_y  = psi_y  ; rho_y  = psi_y  ;  Ti_y  = psi_y  ; vpar_y  = psi_y ; Te_y = psi_y  ; T_y = psi_y
                  u_p  = psi_p  ;  zj_p  = psi_p  ;  w_p  = psi_p  ; rho_p  = psi_p  ;  Ti_p  = psi_p  ; vpar_p  = psi_p ; Te_p = psi_p  ; T_p = psi_p
                  u_s  = psi_s  ;  zj_s  = psi_s  ;  w_s  = psi_s  ; rho_s  = psi_s  ;  Ti_s  = psi_s  ; vpar_s  = psi_s ; Te_s = psi_s  ; T_s = psi_s
                  u_t  = psi_t  ;  zj_t  = psi_t  ;  w_t  = psi_t  ; rho_t  = psi_t  ;  Ti_t  = psi_t  ; vpar_t  = psi_t ; Te_t = psi_t  ; T_t = psi_t
                  u_ss = psi_ss ;  zj_ss = psi_ss ;  w_ss = psi_ss ; rho_ss = psi_ss ;  Ti_ss = psi_ss ; vpar_ss = psi_ss; Te_ss = psi_ss; T_ss = psi_ss
                  u_tt = psi_tt ;  zj_tt = psi_tt ;  w_tt = psi_tt ; rho_tt = psi_tt ;  Ti_tt = psi_tt ; vpar_tt = psi_tt; Te_tt = psi_tt; T_tt = psi_tt
                  u_st = psi_st ;  zj_st = psi_st ;  w_st = psi_st ; rho_st = psi_st ;  Ti_st = psi_st ; vpar_st = psi_st; Te_st = psi_st; T_st = psi_st

                  u_xx = psi_xx ;                    w_xx = psi_xx ; rho_xx = psi_xx ;  Ti_xx = psi_xx ; vpar_xx = psi_xx; Te_xx = psi_xx; T_xx = psi_xx
                  u_yy = psi_yy ;                    w_yy = psi_yy ; rho_yy = psi_yy ;  Ti_yy = psi_yy ; vpar_yy = psi_yy; Te_yy = psi_yy; T_yy = psi_yy
                  u_xy = psi_xy ;                    w_xy = psi_xy ; rho_xy = psi_xy ;  Ti_xy = psi_xy ; vpar_xy = psi_xy; Te_xy = psi_xy; T_xy = psi_xy

                  rhon   = psi
                  rhon_x = psi_x
                  rhon_y = psi_y
                  rhon_p = psi_p
                  rhon_s = psi_s
                  rhon_t = psi_t
                  rhon_ss = psi_ss
                  rhon_tt = psi_tt
                  rhon_st = psi_st

                  rhon_xx = psi_xx
                  rhon_yy = psi_yy
                  rhon_xy = psi_xy


                  rho_hat   = BigR**2 * rho
                  rho_x_hat = 2.d0 * BigR * BigR_x  * rho + BigR**2 * rho_x
                  rho_y_hat = BigR**2 * rho_y

                  Btheta2_psi        = 2.d0 * (psi_x * ps0_x + psi_y * ps0_y ) /BigR**2
                  Bgrad_rho_star_psi = ( v_x   * psi_y - v_y   * psi_x ) / BigR
                  Bgrad_rho_psi      = ( r0_x  * psi_y - r0_y  * psi_x ) / BigR
                  Bgrad_rho_rho      = ( rho_x * ps0_y - rho_y * ps0_x ) / BigR
                  Bgrad_rho_rho_n    = ( F0 / BigR * rho_p ) / BigR
                  BB2_psi            = 2.d0 * (psi_x * ps0_x + psi_y * ps0_y ) /BigR**2

                  Pi0_x_rho  = rho_x  * Ti0 +       rho   * Ti0_x
                  Pi0_xx_rho = rho_xx * Ti0 + 2.0 * rho_x * Ti0_x  + rho   * Ti0_xx
                  Pi0_y_rho  = rho_y  * Ti0 +       rho   * Ti0_y
                  Pi0_yy_rho = rho_yy * Ti0 + 2.0 * rho_y * Ti0_y  + rho   * Ti0_yy
                  Pi0_xy_rho = rho_xy * Ti0 +       rho   * Ti0_xy + rho_x * Ti0_y + rho_y * Ti0_x
                 
                  Pi0_x_Ti   = r0_x  * Ti +       r0   * Ti_x
                  Pi0_xx_Ti  = r0_xx * Ti + 2.0 * r0_x * Ti_x  + r0   * Ti_xx
                  Pi0_y_Ti   = r0_y  * Ti +       r0   * Ti_y
                  Pi0_yy_Ti  = r0_yy * Ti + 2.0 * r0_y * Ti_y  + r0   * Ti_yy
                  Pi0_xy_Ti  = r0_xy * Ti +       r0   * Ti_xy + r0_x * Ti_y  + r0_y * Ti_x
  
                  if (Wdia) then
                    W_dia_rho = - tauIC*2. *     rho/r0_corr**2 * (Pi0_xx     + Pi0_x    /bigR + Pi0_yy    ) &
                                + tauIC*2.          /r0_corr    * (Pi0_xx_rho + Pi0_x_rho/bigR + Pi0_yy_rho) &
                                + tauIC*2. * 2.0*rho/r0_corr**3 * (r0_x *Pi0_x     + r0_y *Pi0_y    )        &
                                - tauIC*2.          /r0_corr**2 * (rho_x*Pi0_x     + rho_y*Pi0_y    )        &
                                - tauIC*2.          /r0_corr**2 * (r0_x *Pi0_x_rho + r0_y *Pi0_y_rho)
                    W_dia_Ti  = + tauIC*2.          /r0_corr    * (Pi0_xx_Ti + Pi0_x_Ti/bigR + Pi0_yy_Ti)    &
                                - tauIC*2.          /r0_corr**2 * (r0_x  *Pi0_x_Ti + r0_y  *Pi0_y_Ti)
                  else
                    W_dia_rho = 0.d0
                    W_dia_Ti  = 0.d0
                  endif
  
                  if (normalized_velocity_profile) then
                     Vt_x_psi    = dV_dpsi_source(ms,mt) * psi_x
                     Vt_y_psi    = dV_dpsi_source(ms,mt) * psi_y
                  else
                     Omega_tor_x_psi = dV_dpsi_source(ms,mt)*psi_x
                     Omega_tor_y_psi = dV_dpsi_source(ms,mt)*psi_y
                  endif

                  !###################################################################################################
                  !#  Induction Equation                                                                             #
                  !###################################################################################################

                  amat(var_psi,var_psi) = v * psi / BigR * xjac * (1.d0 + zeta)                                                     &
                            - v * (psi_s * u0_t - psi_t * u0_s)                                                     * theta * tstep &
                                                                                                                   
                     + v * tauIC*2./(r0_corr*BB2) * F0**2/BigR**2 * (psi_s * Pe0_t - psi_t * Pe0_s)                 * theta * tstep &
                     - v * tauIC*2./(r0_corr*BB2**2) * BB2_psi * F0**2/BigR**2 * (ps0_x*Pe0_y - ps0_y*Pe0_x) * xjac * theta * tstep &
                     + v * tauIC*2./(r0_corr*BB2**2) * BB2_psi * F0**3/BigR**3 * Pe0_p                       * xjac * theta * tstep
                                                                                                             
                  amat(var_psi,var_u) = -  v * (ps0_s * u_t - ps0_t * u_s)                                          * theta * tstep
                                                                                                                    
                  amat_n(var_psi,var_u) = + F0 / BigR * v * u_p * xjac                                              * theta * tstep
                                                                                                             
                  amat(var_psi,var_zj) = - eta_num_T * (v_x * zj_x + v_y * zj_y)                             * xjac * theta * tstep  &
                              - eta_T * v * zj / BigR                                                        * xjac * theta * tstep

                  amat(var_psi,var_rho) = + v * tauIC*2./(r0_corr*BB2) * F0**2/BigR**2 * Te0 * (ps0_s * rho_t - ps0_t * rho_s) * theta * tstep &
                              + v * tauIC*2./(r0_corr*BB2) * F0**2/BigR**2 * rho * (ps0_s * Te0_t - ps0_t * Te0_s)  * theta * tstep &
                              - v * tauIC*2./(r0_corr*BB2) * F0**3/BigR**3 * rho * Te0_p                     * xjac * theta * tstep &

                              - v * tauIC*2. * rho /(r0_corr**2 * BB2) * F0**2/BigR**2 * (ps0_s * Pe0_t - ps0_t * Pe0_s) * theta * tstep &
                              + v * tauIC*2. * rho /(r0_corr**2 * BB2) * F0**3/BigR**3 * Pe0_p                    * xjac * theta * tstep 

                  amat_n(var_psi,var_rho) = - v * tauIC*2./(r0_corr*BB2) * F0**3/BigR**3 * Te0  * rho_p           * xjac * theta * tstep 

                  if ( with_TiTe ) then ! (with_TiTe) **********************************************
                    amat(var_psi,var_Te) = - deta_dT * v * Te * (zj0 - current_source(ms,mt) - Jb)/ BigR * xjac  * theta * tstep &
                                   - deta_num_dT * Te * (v_x * zj0_x + v_y * zj0_y)              * xjac          * theta * tstep &
                              + v * tauIC*2./(r0_corr*BB2) * F0**2/BigR**2 * r0 * (ps0_s * Te_t  - ps0_t * Te_s) * theta * tstep &
                              + v * tauIC*2./(r0_corr*BB2) * F0**2/BigR**2 * Te * (ps0_s * r0_t - ps0_t * r0_s)  * theta * tstep &
                              - v * tauIC*2./(r0_corr*BB2) * F0**3/BigR**3 * Te * r0_p * xjac                    * theta * tstep

                    amat_n(var_psi,var_Te) = - v * tauIC*2./(r0_corr*BB2) * F0**3/BigR**3 * r0 * Te_p     * xjac * theta * tstep
                  else ! (with_TiTe = .f.), i.e. with single temperature *********************************
                    amat(var_psi,var_T) = - deta_dT * v * T * (zj0 - current_source(ms,mt) - Jb)/ BigR    * xjac * theta * tstep &
                                   - deta_num_dT * T * (v_x * zj0_x + v_y * zj0_y)                        * xjac * theta * tstep &
                              + v * tauIC/(r0_corr*BB2) * F0**2/BigR**2 * r0 * (ps0_s * T_t  - ps0_t * T_s)   * theta * tstep &
                              + v * tauIC/(r0_corr*BB2) * F0**2/BigR**2 * T  * (ps0_s * r0_t - ps0_t * r0_s)  * theta * tstep &
                              - v * tauIC/(r0_corr*BB2) * F0**3/BigR**3 * T  * r0_p                    * xjac * theta * tstep

                    amat_n(var_psi,var_T) = - v * tauIC/(r0_corr*BB2) * F0**3/BigR**3 * r0 * T_p       * xjac * theta * tstep
                  end if ! (with_TiTe) *************************************************************

                  !###################################################################################################
                  !#  Perpendicular Momentum Equation                                                                #
                  !###################################################################################################

                  amat(var_u,var_psi) = - v * (psi_s * zj0_t - psi_t * zj0_s)                          * theta * tstep &

                                        ! New terms coming from -(\partial_t \rho + \nabla \cdot (\rho \mathbf{v})) \mathbf{v} in RHS of momentum equation
                                        ! (see wiki: https://www.jorek.eu/wiki/doku.php?id=model500_501_555#equations):
                                        + fact_conservative_u * ( &
                                        - BigR**2 * r0 * (vpar0_x * psi_y - vpar0_y * psi_x) * (v_x * u0_x + v_y * u0_y) * xjac * theta * tstep &
                                        - BigR**2 * vpar0 * (r0_x * psi_y - r0_y * psi_x)    * (v_x * u0_x + v_y * u0_y) * xjac * theta * tstep &
                                                                 )

                  ! ------------------------------------------------------ NEO
                  if (NEO) then
                    amat(var_u,var_psi) = amat(var_u,var_psi) &
                         - amu_neo_prof(ms,mt)*BB2/((Btheta2+epsil)**2)*(psi_x*v_x+psi_y*v_y)       &
                         * (r0 * (ps0_x*u0_x + ps0_y*u0_y) + tauIC*2. * (ps0_x*Pi0_x + ps0_y*Pi0_y) &
                         + aki_neo_prof(ms,mt) * tauIC*2. * r0 * (ps0_x*Ti0_x + ps0_y*Ti0_y) -r0 * Vpar0 * Btheta2) * BigR * xjac * theta * tstep &
                         - amu_neo_prof(ms,mt) * BB2/((Btheta2+epsil)**2) * (ps0_x*v_x + ps0_y*v_y) * (r0*(psi_x*u0_x + psi_y*u0_y) &
                              + tauIC*2. * (psi_x*Pi0_x + psi_y*Pi0_y)                                     &
                         + aki_neo_prof(ms,mt) * tauIC*2. * r0 * (psi_x*Ti0_x + psi_y*Ti0_y)) * BigR * xjac * theta * tstep &

                         ! ========= linearization of 1/(Btheta2**i) , i=2 or 1
                         - amu_neo_prof(ms,mt) * BB2 * (-2.d0*Btheta2_psi)/((Btheta2+epsil)**3) * (ps0_x*v_x + ps0_y*v_y)  &
                                                     * (r0 * (ps0_x*u0_x + ps0_y*u0_y) + tauIC*2. * (ps0_x*Pi0_x + ps0_y*Pi0_y) &
                                                        + aki_neo_prof(ms,mt) * tauIC*2. * r0 * (ps0_x*Ti0_x + ps0_y*Ti0_y) &
                                                        - r0 * Vpar0 * Btheta2) * BigR * xjac * theta * tstep     &
                         + amu_neo_prof(ms,mt) * BB2 * (-Btheta2_psi)/((Btheta2+epsil)**2) * r0 * vpar0 * (ps0_x*v_x + ps0_y*v_y) &
                                               * BigR * xjac * tstep * theta
                  endif

                  amat(var_u,var_u) = - BigR**3 * r0_corr * (v_x * u_x + v_y * u_y) * xjac * (1.d0 + zeta)                                &
                                    + r0_hat * BigR**2 * w0 * (v_s * u_t  - v_t  * u_s)                                   * theta * tstep &
                                    + BigR**2 * (u_x * u0_x + u_y * u0_y) * (v_x * r0_y_hat - v_y * r0_x_hat)      * xjac * theta * tstep &

                                    + tauIC*2. * BigR**3 * Pi0_y * (v_x* u_x + v_y * u_y)                          * xjac * theta * tstep &

                                    + v * tauIC*2. * BigR**4 * (u_xy * (Pi0_xx - Pi0_yy) - Pi0_xy * (u_xx - u_yy)) * xjac * theta * tstep &

                                    ! New terms coming from -(\partial_t \rho + \nabla \cdot (\rho \mathbf{v})) \mathbf{v} in RHS of momentum equation
                                    ! (see wiki: https://www.jorek.eu/wiki/doku.php?id=model500_501_555#equations):
                                    + fact_conservative_u * ( &
                                        + BigR**2 * (r0_x_hat * u0_y - r0_y_hat * u0_x)      * (v_x * u_x  + v_y * u_y)  * xjac * theta * tstep &
                                        + BigR**2 * (r0_x_hat * u_y  - r0_y_hat * u_x)       * (v_x * u0_x + v_y * u0_y) * xjac * theta * tstep &
                                        - BigR * F0 * (r0 * vpar0_p + vpar0 * r0_p)          * (v_x * u_x  + v_y * u_y)  * xjac * theta * tstep &
                                        - BigR**2 * r0 * (vpar0_x * ps0_y - vpar0_y * ps0_x) * (v_x * u_x  + v_y * u_y)  * xjac * theta * tstep &
                                        - BigR**2 * vpar0 * (r0_x * ps0_y - r0_y * ps0_x)    * (v_x * u_x + v_y * u_y)   * xjac * theta * tstep &
                                                            ) &
                                    ! Not to be included in conservative form
                                    - BigR**3 * (particle_source(ms,mt)+source_pellet) * (v_x * u_x + v_y * u_y)   * xjac * theta * tstep &
                                              * (1.d0 - fact_conservative_u) &

                                    + tgnum_u * 0.25d0 * r0_hat * BigR**3 * (w0_x * u_y - w0_y * u_x)                                 &
                                              * ( v_x * u0_y - v_y * u0_x)                             * xjac * theta * tstep * tstep &

                                    + tgnum_u * 0.25d0 * r0_hat * BigR**3 * (w0_x * u0_y - w0_y * u0_x)                                         &
                                      * ( v_x * u_y - v_y * u_x)                                                 * xjac * theta * tstep * tstep &
!====================================New TG_num terms=================================
                                    + tgnum_u * 0.25d0 * w0 * BigR**3 * (r0_x_hat * u_y - r0_y_hat * u_x) &
                                              * ( v_x * u0_y - v_y * u0_x) * theta * xjac * tstep * tstep * fact_conservative_u &
              
                                    + tgnum_u * 0.25d0 * w0 * BigR**3 * (r0_x_hat * u0_y - r0_y_hat * u0_x) &
                                              * ( v_x * u_y - v_y * u_x) * theta * xjac * tstep * tstep   * fact_conservative_u & 

!===============================End of NewTG_num terms==============================
                  
                                    + (1.d0 - delta_n_convection) * (  &
                                    - BigR**3 * (r0*rn0*Sion_T) * (v_x * u_x + v_y * u_y)                        * xjac * theta * tstep &
                                    + BigR**3 * (r0* r0*Srec_T) * (v_x * u_x + v_y * u_y)                        * xjac * theta * tstep &
                                    )

                  if ( NEO ) then
                    amat(var_u,var_u) = amat(var_u,var_u) &
                              - amu_neo_prof(ms,mt)*BB2/((Btheta2+epsil)**2) * (ps0_x*v_x + ps0_y*v_y) * r0 *(ps0_x*u_x + ps0_y*u_y) &
                              * BigR * xjac * theta * tstep
                  endif
                  !---------------------------------------- NEO

                  amat(var_u,var_zj)   = - v * (ps0_s * zj_t  - ps0_t * zj_s)                                              * theta * tstep
                                                                                                                    
                  amat_n(var_u,var_zj) = + F0 / BigR * v * zj_p  * xjac                                                    * theta * tstep
                                                                                                                    
                  amat(var_u,var_w) = r0_hat * BigR**2 * w  * ( v_s * u0_t - v_t * u0_s)                                   * theta * tstep &
                                    + BigR * ( v_x * w_x + v_y * w_y) * visco_T  * xjac                                    * theta * tstep &
                                    + v * tauIC*2. * BigR**4 * (Pi0_s * w_t - Pi0_t * w_s)                                 * theta * tstep &
                                                                                                                    
                                    + visco_num_T * (v_xx + v_x/BigR + v_yy)*(w_xx + w_x/BigR + w_yy)               * xjac * theta * tstep &
                                                                                                                    
                                    + tgnum_u * 0.25d0 * r0_hat * BigR**3 * (w_x * u0_y - w_y * u0_x)                                      &
                                            * ( v_x * u0_y - v_y * u0_x) * xjac * theta * tstep * tstep &
!====================================New TG_num terms=================================
                                    + tgnum_u * 0.25d0 * w * BigR**3 * (r0_x_hat * u0_y - r0_y_hat * u0_x) &
                                              * ( v_x * u0_y - v_y * u0_x) * theta * xjac * tstep * tstep * fact_conservative_u 

!===============================End of NewTG_num terms==============================

                  amat(var_u,var_rho) = + 0.5d0 * vv2 * (v_x * rho_y_hat - v_y * rho_x_hat)                         * xjac * theta * tstep &
                                      + rho_hat * BigR**2 * w0 * (v_s * u0_t   - v_t * u0_s)                               * theta * tstep &
                                      - BigR**2 * (v_s * rho_t * (Ti0+Te0)     - v_t * rho_s * (Ti0+Te0))                  * theta * tstep &
                                      - BigR**2 * (v_s * rho   * (Ti0_t+Te0_t) - v_t * rho   * (Ti0_s+Te0_s))              * theta * tstep &

                                      ! New terms coming from -(\partial_t \rho + \nabla \cdot (\rho \mathbf{v})) \mathbf{v} in RHS of momentum equation
                                      ! (see wiki: https://www.jorek.eu/wiki/doku.php?id=model500_501_555#equations):	     	    	     
                                      + fact_conservative_u *( &
                                          - BigR**3 * rho * (v_x * u0_x + v_y * u0_y) * xjac  * (1.d0 + zeta)          &
                                          + BigR**2 * (rho_x_hat * u0_y - rho_y_hat * u0_x)     * (v_x * u0_x + v_y * u0_y) * xjac * theta * tstep &
                                          - BigR * rho * F0 * vpar0_p                           * (v_x * u0_x + v_y * u0_y) * xjac * theta * tstep &
                                          - BigR**2 * rho * (vpar0_x * ps0_y - vpar0_y * ps0_x) * (v_x * u0_x + v_y * u0_y) * xjac * theta * tstep &
                                          - BigR**2 * vpar0 * (rho_x * ps0_y - rho_y * ps0_x)   * (v_x * u0_x + v_y * u0_y) * xjac * theta * tstep &
                                                             ) &

                                      + v * tauIC*2. * BigR**4 * Ti0  * (rho_s * w0_t - rho_t * w0_s)                      * theta * tstep &
                                      + v * tauIC*2. * BigR**4 * rho * (Ti0_s  * w0_t - Ti0_t  * w0_s)                     * theta * tstep &
                                      + tauIC*2. * BigR**3 * (Ti0_y * rho + Ti0 * rho_y) * (v_x* u0_x + v_y * u0_y) * xjac * theta * tstep &
                                      + v * tauIC*2. * BigR**4 * ( (u0_xy * (rho_xx*Ti0 + 2.d0*rho_x*Ti0_x + rho*Ti0_xx                    &
                                            -  rho_yy*Ti0 - 2.d0*rho_y*Ti0_y - rho*Ti0_yy))                                                &
                                            - (rho_xy*Ti0 + rho_x*Ti0_y + rho_y*Ti0_x + rho*Ti0_xy) * (u0_xx - u0_yy)  ) * xjac * theta * tstep&
                                      ! --- Diamagnetic viscosity
                                      - dvisco_dT * bigR * W_dia_rho * (v_x*Ti0_x + v_y*Ti0_y)                           * xjac * theta * tstep&
                                      - visco_T   * bigR * W_dia_rho * (v_xx + v_x/bigR + v_yy)                          * xjac * theta * tstep&

                                      + tgnum_u * 0.25d0 * rho_hat * BigR**3 * (w0_x * u0_y - w0_y * u0_x)                                            &
                                      * ( v_x * u0_y - v_y * u0_x)                                           * xjac * theta * tstep * tstep     &
!====================================New TG_num terms=================================
                                      + tgnum_u * 0.25d0 * w0 * BigR**3 * (rho_x_hat * u0_y - rho_y_hat * u0_x) &
                                                * ( v_x * u0_y - v_y * u0_x) * theta * xjac * tstep * tstep * fact_conservative_u &

!===============================End of NewTG_num terms==============================

                                      + (1.d0 - delta_n_convection) * (                                                                &
                                        - BigR**3 * (rho* rn0 * Sion_T)        * (v_x * u0_x + v_y * u0_y)   * xjac * theta * tstep &
                                        + BigR**3 * (rho * 2.d0 * r0 * Srec_T) * (v_x * u0_x + v_y * u0_y)   * xjac * theta * tstep &
                                        )

                  if ( NEO ) then
                    amat(var_u,var_rho) = amat(var_u,var_rho) &
                         - amu_neo_prof(ms,mt)*BB2/((Btheta2+epsil)**2) * (ps0_x*v_x + ps0_y*v_y) * (rho*(ps0_x*u0_x+ps0_y*u0_y) &
                         + tauIC*2. * (ps0_x*(rho_x*Ti0 + rho*Ti0_x) + ps0_y*(rho_y*Ti0 + rho*Ti0_y))                                   &
                         + aki_neo_prof(ms,mt) * tauIC*2. * rho * (ps0_x*Ti0_x + ps0_y*Ti0_y) - rho * Vpar0 * Btheta2) * BigR * xjac * tstep * theta
                  endif

                  ! New term coming from -(\partial_t \rho + \nabla \cdot (\rho \mathbf{v})) \mathbf{v} in RHS of momentum equation
                  ! (see wiki: https://www.jorek.eu/wiki/doku.php?id=model500_501_555#equations):	 
                  amat_n(var_u,var_rho) = - BigR * vpar0 * F0 * rho_p * (v_x * u0_x + v_y * u0_y) * xjac * theta * tstep * fact_conservative_u

                  if (with_vpar) amat(var_u,var_vpar) = 0.d0

                  if ( with_TiTe ) then ! (with_TiTe) **********************************************
                    amat(var_u,var_Ti) = - BigR**2 * (v_s * r0_t * Ti   - v_t * r0_s * Ti)           * theta * tstep  &
                                         - BigR**2 * (v_s * r0   * Ti_t - v_t * r0   * Ti_s)         * theta * tstep  &

                                       + v * tauIC*2. * BigR**4 * r0 * (Ti_s * w0_t - Ti_t * w0_s)                      * theta * tstep &
                                       + v * tauIC*2. * BigR**4 * Ti  * (r0_s * w0_t - r0_t * w0_s)                     * theta * tstep &
                                       + tauIC*2. * BigR**3 * (r0_y * Ti + r0 * Ti_y) * (v_x* u0_x + v_y * u0_y) * xjac * theta * tstep &
                                       + v * tauIC*2. * BigR**4 * ( (u0_xy * (Ti_xx*r0 + 2.d0*Ti_x*r0_x + Ti*r0_xx                      &
                                                                         - Ti_yy*r0 - 2.d0*Ti_y*r0_y - Ti*r0_yy))                    &
                                                               - (Ti_xy * r0 + Ti_x*r0_y + Ti_y*r0_x + Ti*r0_xy) * (u0_xx - u0_yy))  &
                                                             * xjac * theta * tstep                                                  &

                                       - dvisco_dT     * bigR * W_dia_Ti * (v_x*Ti0_x + v_y*Ti0_y)  * xjac * theta * tstep  &
                                       - visco_T       * bigR * W_dia_Ti * (v_xx + v_x/bigR + v_yy) * xjac * theta * tstep  &
                                       - dvisco_dT     * bigR * W_dia   * (v_x*Ti_x  + v_y*Ti_y )   * xjac * theta * tstep

                    if ( NEO ) then
                      amat(var_u,var_Ti) = amat(var_u,var_Ti) - amu_neo_prof(ms,mt)*BB2/((Btheta2+epsil)**2)*(ps0_x*v_x + ps0_y*v_y) &
                                * (tauIC*2.*(ps0_x*(r0_x*Ti+r0*Ti_x) + ps0_y*(r0_y*Ti+r0*Ti_y))                       &
                                + aki_neo_prof(ms,mt) *tauIC*2. * r0 *(ps0_x*Ti_x + ps0_y*Ti_y)) * BigR * xjac * theta * tstep 

                      if ( with_vpar ) &
                      amat(var_u,var_vpar) =  amat(var_u,var_vpar) + amu_neo_prof(ms,mt)*BB2 * Btheta2 /((Btheta2+epsil)**2)*r0*vpar*(ps0_x*v_x+ps0_y*v_y) &
                                * BigR * xjac * tstep * theta
                    endif


                    amat(var_u,var_Te) = - BigR**2 * (v_s * r0_t * Te   - v_t * r0_s * Te)           * theta * tstep  &
                                - BigR**2 * (v_s * r0   * Te_t - v_t * r0   * Te_s)         * theta * tstep  &

                                + dvisco_dT * Te * ( v_x * w0_x + v_y * w0_y ) * BigR * xjac * theta * tstep  &

                                - d2visco_dT2*Te * bigR * W_dia   * (v_x*Ti0_x + v_y*Ti0_y)  * xjac * theta * tstep  &
                                - dvisco_dT*Te   * bigR * W_dia   * (v_xx + v_x/bigR + v_yy) * xjac * theta * tstep

                  else ! (with_TiTe = .f.), i.e. with single temperature *********************************

                    amat(var_u,var_T) = - BigR**2 * (v_s * r0_t * T   - v_t * r0_s * T)           * theta * tstep  &
                                        - BigR**2 * (v_s * r0   * T_t - v_t * r0   * T_s)         * theta * tstep  &

                              + v * tauIC * BigR**4 * r0 * (T_s  * w0_t - T_t  * w0_s)                      * theta * tstep &
                              + v * tauIC * BigR**4 * T  * (r0_s * w0_t - r0_t * w0_s)                     * theta * tstep &
                              + tauIC * BigR**3 * (r0_y * T + r0 * T_y) * (v_x* u0_x + v_y * u0_y) * xjac * theta * tstep &
                              + v * tauIC * BigR**4 * ( (u0_xy * (T_xx*r0 + 2.d0*T_x*r0_x + T*r0_xx                      &
                                                                - T_yy*r0 - 2.d0*T_y*r0_y - T*r0_yy))                    &
                                                      - (T_xy * r0 + T_x*r0_y + T_y*r0_x + T*r0_xy) * (u0_xx - u0_yy))  &
                                                    * xjac * theta * tstep                                                  &

                              + dvisco_dT * T * ( v_x * w0_x + v_y * w0_y )         * BigR * xjac * theta * tstep  &

                              ! --- Contributions of the diamagnetic viscosity 
                              - dvisco_dT     * bigR * W_dia_Ti * (v_x*Ti0_x + v_y*Ti0_y)  * xjac * theta * tstep  &
                              - visco_T       * bigR * W_dia_Ti * (v_xx + v_x/bigR + v_yy) * xjac * theta * tstep  &
                              - dvisco_dT     * bigR * W_dia    * (v_x*T_x  + v_y*T_y )    * xjac * theta * tstep  &

                              - d2visco_dT2*T * bigR * W_dia    * (v_x*Ti0_x + v_y*Ti0_y)  * xjac * theta * tstep  &
                              - dvisco_dT*T   * bigR * W_dia    * (v_xx + v_x/bigR + v_yy) * xjac * theta * tstep  &

                           + (1 - delta_n_convection) * (  &
                           - BigR**3 * (r0 * rn0 * dSion_dT * T) * (v_x * u0_x + v_y * u0_y) * xjac * theta * tstep  &
                           + BigR**3 * (r0 * r0  * dSrec_dT * T) * (v_x * u0_x + v_y * u0_y) * xjac * theta * tstep  &
                           )


                    if ( NEO ) then
                      amat(var_u,var_T) = amat(var_u,var_T) - amu_neo_prof(ms,mt)*BB2/((Btheta2+epsil)**2)*(ps0_x*v_x + ps0_y*v_y) &
                                * (tauIC*(ps0_x*(r0_x*T+r0*T_x) + ps0_y*(r0_y*T+r0*T_y))                       &
                                + aki_neo_prof(ms,mt) *tauIC * r0 *(ps0_x*T_x + ps0_y*T_y)) * BigR * xjac * theta * tstep

                      if ( with_vpar ) &
                      amat(var_u,var_vpar) = amat(var_u,var_vpar) + amu_neo_prof(ms,mt)*BB2 * Btheta2 /((Btheta2+epsil)**2)*r0*vpar*(ps0_x*v_x+ps0_y*v_y) &
                                * BigR * xjac * tstep * theta
                    endif

                  end if ! (with_TiTe) *************************************************************

                  if (with_vpar) then
                    ! New terms coming from -(\partial_t \rho + \nabla \cdot (\rho \mathbf{v})) \mathbf{v} in RHS of momentum equation
                    ! (see wiki: https://www.jorek.eu/wiki/doku.php?id=model500_501_555#equations):	     	    	     	  
                    amat(var_u,var_vpar) =   amat(var_u,var_vpar)   &
                          + fact_conservative_u * ( & 
                              - BigR * vpar * F0 * r0_p                          * (v_x * u0_x + v_y * u0_y) * xjac * theta * tstep &
                              - BigR**2 * r0 * (vpar_x * ps0_y - vpar_y * ps0_x) * (v_x * u0_x + v_y * u0_y) * xjac * theta * tstep &
                              - BigR**2 * vpar * (r0_x * ps0_y - r0_y * ps0_x)   * (v_x * u0_x + v_y * u0_y) * xjac * theta * tstep &
                                                   )

                    ! New term coming from -(\partial_t \rho + \nabla \cdot (\rho \mathbf{v})) \mathbf{v} in RHS of momentum equation
                    ! (see wiki: https://www.jorek.eu/wiki/doku.php?id=model500_501_555#equations):	 
                    amat_n(var_u,var_vpar) = - BigR * r0 * F0 * vpar_p * (v_x * u0_x + v_y * u0_y) * xjac * theta * tstep * fact_conservative_u
                  endif ! /with_vpar

                  if (with_neutrals) then
                    amat(var_u,var_rhon) = -(1.d0 - delta_n_convection) * BigR**3 * (r0 * rhon * Sion_T) * (v_x * u0_x + v_y * u0_y) * xjac * theta * tstep 
                  endif

                  !###################################################################################################
                  !#  Current Definition Equation                                                                    #
                  !###################################################################################################

                  amat(var_zj,var_zj) = v * zj / BigR * xjac
                  amat(var_zj,var_psi) = (v_x * psi_x + v_y * psi_y ) / BigR * xjac

                  !###################################################################################################
                  !#  Vorticity Definition Equation                                                                  #
                  !###################################################################################################

                  amat(var_w,var_w) =  v * w * BigR * xjac                                
                  amat(var_w,var_u) = (v_x * u_x + v_y * u_y) * BigR * xjac              

                  !###################################################################################################
                  !#  Density Equation                                                                               #
                  !###################################################################################################


                  amat(var_rho,var_psi) =-((D_par+D_par_sc_num*tau_sc)-D_prof) * BigR * BB2_psi/ BB2**2 * Bgrad_rho_star     * Bgrad_rho     * xjac * theta * tstep &
                                        + ((D_par+D_par_sc_num*tau_sc)-D_prof) * BigR / BB2             * Bgrad_rho_star_psi * Bgrad_rho     * xjac * theta * tstep &
                                        + ((D_par+D_par_sc_num*tau_sc)-D_prof) * BigR / BB2             * Bgrad_rho_star     * Bgrad_rho_psi * xjac * theta * tstep &
                                        + v * Vpar0 * (r0_s * psi_t - r0_t * psi_s)                                           * theta * tstep &
                                        + v * r0 * (vpar0_s * psi_t - vpar0_t * psi_s)                                        * theta * tstep &
             
                                        + tgnum_rho * 0.25d0 / BigR * vpar0**2                                                                          &
                                                  * (r0_x * psi_y - r0_y * psi_x)                                                                     &
                                                  * ( v_x * ps0_y -  v_y * ps0_x                   )                   * xjac * theta * tstep * tstep &
                                        + tgnum_rho * 0.25d0 / BigR * vpar0**2                                                                          &
                                                  * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                                                  &
                                                  * ( v_x * psi_y -  v_y * psi_x                   )                   * xjac * theta * tstep * tstep


                  amat_k(var_rho,var_psi) =-((D_par+D_par_sc_num*tau_sc)-D_prof) * BigR * BB2_psi/ BB2**2 * Bgrad_rho_k_star * Bgrad_rho     * xjac * theta * tstep &
                                          + ((D_par+D_par_sc_num*tau_sc)-D_prof) * BigR / BB2             * Bgrad_rho_k_star * Bgrad_rho_psi * xjac * theta * tstep &

                                          + tgnum_rho * 0.25d0 / BigR * vpar0**2                                                                        &
                                                  * (r0_x * psi_y - r0_y * psi_x)                                                                     &
                                                  * (                            + F0 / BigR * v_p)                    * xjac * theta * tstep * tstep

                  amat(var_rho,var_u) =- v * BigR**2 * ( r0_s * u_t - r0_t * u_s)                                        * theta * tstep &
                                       - v * 2.d0 * BigR * r0 * u_y                                               * xjac * theta * tstep &
                                       
                                       + tgnum_rho * 0.25d0 * BigR**3 * (r0_x * u_y  - r0_y * u_x)                                                     &
                                                                    * ( v_x * u0_y - v_y  * u0_x)                * xjac * theta * tstep * tstep      &
                                       + tgnum_rho * 0.25d0 * BigR**3 * (r0_x * u0_y - r0_y * u0_x)                                                    &
                                                                    * ( v_x * u_y  - v_y  * u_x)                 * xjac * theta * tstep * tstep 

                  amat(var_rho,var_rho) = v * rho * BigR * (1.d0 + zeta)                                            * xjac   &
                          - v * BigR**2 * ( rho_s * u0_t - rho_t * u0_s)                                       * theta * tstep &
                          - v * 2.d0 * BigR * rho * u0_y                                                * xjac * theta * tstep &
                          + ((D_par+D_par_sc_num*tau_sc)-D_prof) * BigR / BB2 * Bgrad_rho_star * Bgrad_rho_rho                * xjac * theta * tstep &
                          + D_prof * BigR  * (v_x*rho_x + v_y*rho_y )                                   * xjac * theta * tstep &
                          + v * Vpar0 * (rho_s * ps0_t - rho_t * ps0_s)                                        * theta * tstep &
                          + v * rho * (vpar0_s * ps0_t - vpar0_t * ps0_s)                                      * theta * tstep &
                          + v * rho * F0 / BigR * vpar0_p                                               * xjac * theta * tstep &

                          - v * 2.d0 * tauIC*2. * (rho_y * Ti0 + rho*Ti0_y) * BigR                         * xjac * theta * tstep &

                          - v * rho * rn0       * BigR * Sion_T                                         * xjac * theta * tstep &
                          + v * rho * 2.d0 * r0 * BigR * Srec_T                                         * xjac * theta * tstep &

                          + D_perp_num * (v_xx + v_x/BigR + v_yy)*(rho_xx + rho_x/BigR + rho_yy)  * BigR * xjac * theta * tstep &

                          + tgnum_rho * 0.25d0 * BigR**3 * (rho_x * u0_y - rho_y * u0_x)                                &
                                                       * ( v_x  * u0_y - v_y   * u0_x) * xjac * theta * tstep * tstep &

                          + tgnum_rho * 0.25d0 / BigR * vpar0**2                                                        &
                                    * (rho_x * ps0_y - rho_y * ps0_x )                             &
                                    * ( v_x * ps0_y -  v_y * ps0_x   ) * xjac * theta * tstep * tstep

                  amat_k(var_rho,var_rho) = + ((D_par+D_par_sc_num*tau_sc)-D_prof) * BigR / BB2 * Bgrad_rho_k_star * Bgrad_rho_rho          * xjac * theta * tstep &
 
                         + tgnum_rho * 0.25d0 / BigR * vpar0**2                                                        &
                                    * (rho_x * ps0_y - rho_y * ps0_x                  )                              &
                                    * (                              + F0 / BigR * v_p) * xjac * theta * tstep * tstep

                  amat_n(var_rho,var_rho) = + ((D_par+D_par_sc_num*tau_sc)-D_prof) * BigR / BB2 * Bgrad_rho_star   * Bgrad_rho_rho_n        * xjac * theta * tstep &
                              + v * F0 / BigR * Vpar0 * rho_p                                             * xjac * theta * tstep &

                          + tgnum_rho * 0.25d0 / BigR * vpar0**2                                                        &
                                    * (                              + F0 / BigR * rho_p)                             &
                                    * ( v_x * ps0_y -  v_y * ps0_x                      ) * xjac * theta * tstep * tstep

                  amat_kn(var_rho,var_rho) = + ((D_par+D_par_sc_num*tau_sc)-D_prof) * BigR / BB2 * Bgrad_rho_k_star * Bgrad_rho_rho_n       * xjac * theta * tstep &
                                 + D_prof * BigR  * ( v_p*rho_p /BigR**2 )                                * xjac * theta * tstep &

                          + tgnum_rho * 0.25d0 / BigR * vpar0**2                                                        &
                                    * ( + F0 / BigR * rho_p)                                                          &
                                    * ( + F0 / BigR * v_p) * xjac * theta * tstep * tstep

                  if ( with_TiTe ) then
                    amat(var_rho,var_Ti) = - v * 2.d0 * tauIC*2. * (Ti_y * r0 + Ti*r0_y) * BigR * xjac * theta * tstep
                    amat(var_rho,var_Te) = - v * BigR * r0 * rn0 * dSion_dT * Te   * xjac * theta * tstep &
                                           + v * BigR * r0 * r0 *  dSrec_dT * Te   * xjac * theta * tstep
                  else ! (with_TiTe)
                    amat(var_rho,var_T)  = - v * 2.d0 * tauIC * (T_y  * r0 + T *r0_y) * BigR * xjac * theta * tstep &
                                           - v * BigR * r0 * rn0 * dSion_dT * T    * xjac * theta * tstep &
                                           + v * BigR * r0 * r0 *  dSrec_dT * T    * xjac * theta * tstep 
                  end if ! (with_TiTe)

                  if ( with_vpar ) then
                    amat(var_rho,var_vpar) = + v * F0 / BigR * Vpar * r0_p                * xjac * theta * tstep &
                                             + v * Vpar * (r0_s * ps0_t - r0_t * ps0_s)          * theta * tstep &
                                             + v * r0 * (vpar_s * ps0_t - vpar_t * ps0_s)        * theta * tstep &
                                         
                                             + tgnum_rho * 0.25d0 / BigR * 2.d0*vpar0*vpar                                            &
                                                  * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                                &
                                                  * ( v_x * ps0_y -  v_y * ps0_x                   ) * xjac * theta * tstep * tstep 
  
                    amat_k(var_rho,var_vpar) = + tgnum_rho * 0.25d0 / BigR * 2.d0*vpar0*vpar                            &
                                     * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                               &
                                     * (                            + F0 / BigR * v_p) * xjac * theta * tstep * tstep 
  
                    amat_n(var_rho,var_vpar) = + v * r0 * F0 / BigR * vpar_p                 * xjac * theta * tstep
                  end if ! (with_vpar)

                  if (with_neutrals) then
                    amat(var_rho,var_rhon) = - BigR * v * r0 * Sion_T * rhon            * xjac * theta * tstep
                  endif

                  !###################################################################################################
                  !#  Parallel Velocity Equation                                                                     #
                  !###################################################################################################
                  
                  if ( with_vpar ) then
  
                    amat(var_vpar,var_psi) = v * r0_corr * vpar0 / BigR * (ps0_x * psi_x + ps0_y * psi_y) * xjac * (1.d0 + zeta) &
  
                              + v * (P0_s * psi_t - P0_t * psi_s)                                       * theta * tstep &
  
                              + 0.5d0 * r0 * vpar0**2 * BB2     * (psi_x * v_y  - psi_y * v_x)   * xjac * theta * tstep &
                              + 0.5d0 * r0 * vpar0**2 * BB2_psi * (ps0_x * v_y  - ps0_y * v_x)   * xjac * theta * tstep &
                              + 0.5d0 * v  * vpar0**2 * BB2     * (psi_x * r0_y - psi_y * r0_x)  * xjac * theta * tstep &
                              + 0.5d0 * v  * vpar0**2 * BB2_psi * (ps0_x * r0_y - ps0_y * r0_x)  * xjac * theta * tstep &
                              - 0.5d0 * v  * vpar0**2 * BB2_psi * F0 / BigR * r0_p               * xjac * theta * tstep &

                             ! New terms coming from -(\partial_t \rho + \nabla \cdot (\rho \mathbf{v})) \mathbf{v} in RHS of momentum equation
                             ! (see wiki: https://www.jorek.eu/wiki/doku.php?id=model500_501_555#equations):
                              + fact_conservative_u * ( & 
                                  + v * r0 * (vpar0_x * psi_y - vpar0_y * psi_x) * vpar0 * BB2     * xjac * theta * tstep &
                                  + v * vpar0 * (r0_x * psi_y - r0_y * psi_x)    * vpar0 * BB2     * xjac * theta * tstep &
                                  - v * (r0_x_hat * u0_y - r0_y_hat * u0_x)      * vpar0 * BB2_psi * xjac * theta * tstep &
                                  + v * F0 / BigR * (r0 * vpar0_p + r0_p * vpar0)* vpar0 * BB2_psi * xjac * theta * tstep &
                                  + v * r0 * (vpar0_x * ps0_y - vpar0_y * ps0_x) * vpar0 * BB2_psi * xjac * theta * tstep &
                                  + v * vpar0 * (r0_x * ps0_y - r0_y * ps0_x)    * vpar0 * BB2_psi * xjac * theta * tstep &
                                                        ) &
  
                              ! Not to be included in conservative form
                              + v*(particle_source(ms,mt) + source_pellet)*vpar0* BB2_psi * BigR * xjac * theta * tstep &
                                 * (1.d0 - fact_conservative_u)  &

                              + (1.d0 - delta_n_convection) * (  &  
                                + v *(r0 * rn0 * Sion_T) * vpar0 * BB2_psi * BigR      * xjac * theta * tstep &
                                - v *(r0 * r0  * Srec_T) * vpar0 * BB2_psi * BigR      * xjac * theta * tstep &
                                ) &

                              + tgnum_vpar * 0.25d0 * r0 * Vpar0**2 * BB2 &
                                        * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac) / BigR  &
                                        * (-(psi_s * v_t     - psi_t * v_s)    /xjac)  * xjac * theta * tstep*tstep &
  
                              + tgnum_vpar * 0.25d0 * r0 * Vpar0**2 * BB2 &
                                        * (-(psi_s * vpar0_t - psi_t * vpar0_s)/xjac) / BigR  &
                                        * (-(ps0_s * v_t     - ps0_t * v_s)    /xjac)  * xjac * theta * tstep*tstep &
  
                              + tgnum_vpar * 0.25d0 * v  * Vpar0**2 * BB2 * (1.d0 - fact_conservative_u) &
                                        * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac) / BigR  &
                                        * (-(psi_s * r0_t    - psi_t * r0_s)   /xjac)  * xjac * theta * tstep*tstep &
  
                              + tgnum_vpar * 0.25d0 * v  * Vpar0**2 * BB2 * (1.d0 - fact_conservative_u) &
                                        * (-(psi_s * vpar0_t - psi_t * vpar0_s)/xjac) / BigR  &
                                        * (-(ps0_s * r0_t    - ps0_t * r0_s)   /xjac)  * xjac * theta * tstep*tstep & 

!=============================== New TG_num terms==================================
                               + tgnum_vpar * 0.25d0 * vpar0 * Vpar0**2 * BB2 * fact_conservative_u &
                                         * (-(ps0_s * r0_t - ps0_t * r0_s)/xjac) / BigR  &
                                         * (-(psi_s * v_t     - psi_t * v_s)    /xjac) * xjac * theta * tstep*tstep &
        
                               + tgnum_vpar * 0.25d0 * vpar0 * Vpar0**2 * BB2 * fact_conservative_u &
                                         * (-(psi_s * r0_t - psi_t * r0_s)/xjac) / BigR  &
                                         * (-(ps0_s * v_t     - ps0_t * v_s)    /xjac) * xjac * theta * tstep*tstep
!===============================End of new TG_num terms============================
  
                    if (normalized_velocity_profile) then
                      amat(var_vpar,var_psi) = amat(var_vpar,var_psi)  - (visco_par + visco_par_sc_num * tau_sc) * (v_x * Vt_x_psi   + v_y * Vt_y_psi) * BigR * xjac * theta * tstep      
                    else
                      amat(var_vpar,var_psi) = amat(var_vpar,var_psi)  - (visco_par + visco_par_sc_num * tau_sc) * 2.d0 * PI * F0 * (v_x * Omega_tor_x_psi + v_y * Omega_tor_y_psi) * BigR * xjac * theta * tstep 
                    endif
  
                    amat_k(var_vpar,var_psi) = - 0.5d0 * r0 * vpar0**2 * BB2_psi * F0 / BigR * v_p    * xjac * theta * tstep 
                    amat(var_vpar,var_u) = 0.d0

                    !---------------------------------------- NEO
                    if ( NEO ) then
                      amat(var_vpar,var_psi) = amat(var_vpar,var_psi) &
                           - v * amu_neo_prof(ms,mt)*BB2/(Btheta2+epsil)*(r0 * (psi_x*u0_x + psi_y*u0_y) + tauIC*2.*(psi_x*Pi0_x + psi_y*Pi0_y) &
                           + aki_neo_prof(ms,mt) * tauIC*2. * r0 * (psi_x*Ti0_x + psi_y*Ti0_y)) * BigR * xjac * theta * tstep                   &
                           - v * amu_neo_prof(ms,mt) * (-Btheta2_psi)*BB2/((Btheta2+epsil)**2) * (r0*(ps0_x*u0_x+ps0_y*u0_y)                         &
                                                        + tauIC*2.*(ps0_x*Pi0_x+ps0_y*Pi0_y)                                                    &                     
                           + aki_neo_prof(ms,mt) * tauIC*2. * r0 * (ps0_x*Ti0_x + ps0_y*Ti0_y)) * BigR * xjac * theta * tstep
  
                      amat(var_vpar,var_u) =  amat(var_vpar,var_u) &
                           -v * amu_neo_prof(ms,mt)*BB2/(Btheta2+epsil) * r0 * (ps0_x*u_x + ps0_y*u_y) * BigR * xjac * theta * tstep 
                    endif
                    !---------------------------------------- NEO
  
                    ! New term coming from -(\partial_t \rho + \nabla \cdot (\rho \mathbf{v})) \mathbf{v} in RHS of momentum equation
                    ! (see wiki: https://www.jorek.eu/wiki/doku.php?id=model500_501_555#equations):
                    amat(var_vpar,var_u) = amat(var_vpar,var_u) - v * (r0_x_hat * u_y - r0_y_hat * u_x) * vpar0 * BB2 * theta * xjac * tstep * fact_conservative_u
  
                    amat(var_vpar,var_rho) = + v * (rho_s * (Ti0+Te0)     * ps0_t - rho_t * (Ti0+Te0)     * ps0_s) * theta * tstep &
                                + v * (rho   * (Ti0_s+Te0_s) * ps0_t - rho   * (Ti0_t+Te0_t) * ps0_s) * theta * tstep &
                                + v * F0 / BigR * rho * (Ti0_p+Te0_p)                          * xjac * theta * tstep &
  
                              + 0.5d0 * rho * vpar0**2 * BB2 * (ps0_s * v_t   - ps0_t * v_s)    * theta * tstep &
                              + 0.5d0 * v   * vpar0**2 * BB2 * (ps0_s * rho_t - ps0_t * rho_s)  * theta * tstep &
  
                       ! New terms coming from -(\partial_t \rho + \nabla \cdot (\rho \mathbf{v})) \mathbf{v} in RHS of momentum equation
                       ! (see wiki: https://www.jorek.eu/wiki/doku.php?id=model500_501_555#equations):
                              + fact_conservative_u * ( & 
                                  + v * rho * vpar0 * F0**2 / BigR * xjac * (1.d0 + zeta)  &
                                  - v * (rho_x_hat * u0_y - rho_y_hat * u0_x)       * vpar0 * BB2 * theta * xjac * tstep &   
                                  + v * F0 / BigR * rho * vpar0_p                   * vpar0 * BB2 * theta * xjac * tstep &
                                  + v * rho * (vpar0_x * ps0_y - vpar0_y * ps0_x)   * vpar0 * BB2 * theta * xjac * tstep &
                                  + v * vpar0 * (rho_x * ps0_y - rho_y * ps0_x)     * vpar0 * BB2 * theta * xjac * tstep &
                                                        ) &

                              + (1.d0 - delta_n_convection) * (  &

                              + v *(rho * rn0       * Sion_T) * vpar0 * BB2 * BigR         * xjac * theta * tstep &
                              - v *(2.d0 * r0 * rho * Srec_T) * vpar0 * BB2 * BigR         * xjac * theta * tstep &
                              ) &

                              + tgnum_vpar * 0.25d0 * rho * Vpar0**2 * BB2 &
                                        * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac + F0 / BigR * vpar0_p) / BigR  &
                                        * (-(ps0_s * v_t     - ps0_t * v_s)    /xjac          )  * xjac * theta * tstep*tstep &
  
                              + tgnum_vpar * 0.25d0 * v * Vpar0**2 * BB2 * (1.d0-fact_conservative_u) &
                                        * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac + F0 / BigR * vpar0_p) / BigR  &
                                        * (-(ps0_s * rho_t   - ps0_t * rho_s)  /xjac           ) * xjac * theta *tstep*tstep &
!=============================== New TG_num terms==================================

                               + tgnum_vpar * 0.25d0 * vpar0 * Vpar0**2 * BB2 * fact_conservative_u &
                                         * (-(ps0_s * rho_t - ps0_t * rho_s)/xjac                          ) / BigR  &
                                         * (-(ps0_s * v_t     - ps0_t * v_s)    /xjac)  * xjac * theta * tstep*tstep

!===============================End of new TG_num terms============================
  
                    amat_k(var_vpar,var_rho) = - 0.5d0 * rho * vpar0**2 * BB2 * F0 / BigR * v_p       * xjac * theta * tstep &
  
                                + tgnum_vpar * 0.25d0 * rho * Vpar0**2 * BB2 &
                                          * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac + F0 / BigR * vpar0_p) / BigR  &
                                          * (                                          + F0 / BigR * v_p)  *  xjac * theta * tstep * tstep &
!=============================== New TG_num terms==================================

                                + tgnum_vpar * 0.25d0 * vpar0 * Vpar0**2 * BB2 * fact_conservative_u &
                                          * (-(ps0_s * rho_t - ps0_t * rho_s)/xjac) / BigR  &
                                          * (+ F0 / BigR * v_p) * xjac * theta * tstep*tstep

!===============================End of new TG_num terms============================
  
                    amat_n(var_vpar,var_rho) = + v * F0 / BigR * rho_p * (Ti0+Te0)                    * xjac * theta * tstep &
                                - 0.5d0 * v   * vpar0**2 * BB2 * F0 / BigR * rho_p       * xjac * theta * tstep &
                         ! New term coming from -(\partial_t \rho + \nabla \cdot (\rho \mathbf{v})) \mathbf{v} in RHS of momentum equation
                         ! (see wiki: https://www.jorek.eu/wiki/doku.php?id=model500_501_555#equations):	 
                                + v * vpar0 * F0 / BigR * rho_p * vpar0 * BB2          * theta * xjac * tstep * fact_conservative_u & 
  
                                + tgnum_vpar * 0.25d0 * v * Vpar0**2 * BB2 * (1.d0-fact_conservative_u) &
                                          * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac + F0 / BigR * vpar0_p) / BigR  &
                                          * (                                          + F0 / BigR * rho_p)* xjac * theta * tstep*tstep &
!=============================== New TG_num terms==================================

                                + tgnum_vpar * 0.25d0 * vpar0 * Vpar0**2 * BB2 * fact_conservative_u &
                                          * (+ F0 / BigR * rho_p) / BigR  &
                                          * (-(ps0_s * v_t     - ps0_t * v_s)    /xjac) * xjac * theta * tstep*tstep

!===============================End of new TG_num terms============================
                    
!=============================== New TG_num terms==================================

                      amat_kn(var_vpar,var_rho) = + tgnum_vpar * 0.25d0 * vpar0 * Vpar0**2 * BB2 * fact_conservative_u &
                                          * (+ F0 / BigR * rho_p) / BigR  &
                                          * (+ F0 / BigR * v_p) * xjac * theta * tstep*tstep

!===============================End of new TG_num terms============================

                    if ( with_TiTe ) then ! (with_TiTe) ********************************************
                      amat(var_vpar,var_Ti)   = + v * (Ti_s * r0   * ps0_t - Ti_t * r0   * ps0_s)        * theta * tstep &
                                                + v * (Ti   * r0_s * ps0_t - Ti   * r0_t * ps0_s)        * theta * tstep &
                                                + v * F0 / BigR * Ti * r0_p                       * xjac * theta * tstep

                      amat_n(var_vpar,var_Ti) = + v * F0 / BigR * Ti_p * R0                       * xjac * theta * tstep

                      amat(var_vpar,var_Te)   = + v * (Te_s * r0   * ps0_t - Te_t * r0   * ps0_s)        * theta * tstep  &
                                                + v * (Te   * r0_s * ps0_t - Te   * r0_t * ps0_s)        * theta * tstep  &
                                                + v * F0 / BigR * Te * r0_p                       * xjac * theta * tstep  &

                                                + (1.d0 - delta_n_convection) * (  &
                                                  + v *(r0 * rn0 * dSion_dT *Te) * vpar0 * BB2 * BigR  * xjac * theta * tstep &
                                                  - v *(r0 * r0  * dSrec_dT *Te) * vpar0 * BB2 * BigR  * xjac * theta * tstep &
                                                  )                                               
    

                      amat_n(var_vpar,var_Te) = + v * F0 / BigR * Te_p * r0                       * xjac * theta * tstep
                    else ! (with_TiTe = .f.), i.e. with single temperature *******************************
                      amat(var_vpar,var_T)    = + v * (T_s  * r0   * ps0_t - T_t  * r0   * ps0_s)        * theta * tstep  &
                                                + v * (T    * r0_s * ps0_t - T    * r0_t * ps0_s)        * theta * tstep  &
                                                + v * F0 / BigR * T  * r0_p                       * xjac * theta * tstep  &

                                                + (1.d0 - delta_n_convection) * (  &
                                                  + v *(r0 * rn0 * dSion_dT * T) * vpar0 * BB2 * BigR  * xjac * theta * tstep &
                                                  - v *(r0 * r0  * dSrec_dT * T) * vpar0 * BB2 * BigR  * xjac * theta * tstep &
                                                  )  

                      amat_n(var_vpar,var_T)  = + v * F0 / BigR * T_p  * r0                       * xjac * theta * tstep
                    end if ! (with_TiTe) ***********************************************************

 
                    amat(var_vpar,var_vpar) = v * Vpar * r0_corr * F0**2 / BigR * xjac * (1.d0 + zeta) &
  
                       ! New terms coming from -(\partial_t \rho + \nabla \cdot (\rho \mathbf{v})) \mathbf{v} in RHS of momentum equation
                       ! (see wiki: https://www.jorek.eu/wiki/doku.php?id=model500_501_555#equations):
                            + fact_conservative_u * ( & 
                                - v * (r0_x_hat * u0_y - r0_y_hat * u0_x) * vpar         * BB2 * xjac * theta * tstep &
                                + 2.0 * v * vpar0 * F0 / BigR * r0_p      * vpar         * BB2 * xjac * theta * tstep &
                                + v * r0 * F0 /BigR * vpar0_p             * vpar         * BB2 * xjac * theta * tstep & 
                                + v * r0 * (vpar_x * ps0_y - vpar_y * ps0_x) * vpar0     * BB2 * xjac * theta * tstep &
                                + v * r0 * (vpar0_x * ps0_y - vpar0_y * ps0_x) * vpar    * BB2 * xjac * theta * tstep &
                                + 2.0 * v * vpar0 * (r0_x * ps0_y - r0_y * ps0_x) * vpar * BB2 * xjac * theta * tstep &
                                                    ) &

                             ! Not to be included in conservative form
                            + v*(particle_source(ms,mt) + source_pellet)*vpar*BB2 * BigR * xjac * theta * tstep &
                               *(1.d0 - fact_conservative_u) &
  
                            + r0 * vpar0 * vpar * BB2 * (ps0_s * v_t - ps0_t * v_s)             * theta * tstep &
                            + v  * vpar0 * vpar * BB2 * (ps0_s * r0_t - ps0_t * r0_s)           * theta * tstep &
                            - v  * vpar0 * vpar * BB2 * F0 / BigR * r0_p                 * xjac * theta * tstep &
  
                            + (1.d0 - delta_n_convection) * (  &

                              + v *(r0 * rn0 * Sion_T) * vpar * BB2 * BigR               * xjac * theta * tstep   &
                              - v *(r0 * r0  * Srec_T) * vpar * BB2 * BigR               * xjac * theta * tstep   &
                              ) &

                            + visco_par_num * (v_xx + v_x/BigR + v_yy)*(vpar_xx + vpar_x/BigR + vpar_yy) * BigR * xjac * theta * tstep&
  
                            + tgnum_vpar * 0.5d0 * r0 * Vpar * Vpar0 * BB2 &
                                      * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac + F0 / BigR * vpar0_p) / BigR                      &
                                      * (-(ps0_s * v_t     - ps0_t * v_s)    /xjac                      )  * xjac * theta * tstep*tstep  &
                            + tgnum_vpar * 0.5d0 * v * Vpar * Vpar0 * BB2 * (1.d0 - fact_conservative_u) &
                                      * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac + F0 / BigR * vpar0_p) / BigR                      &
                                      * (-(ps0_s * r0_t    - ps0_t * r0_s)   /xjac + F0 / BigR * r0_p)  * xjac * theta * tstep*tstep  &
                            + tgnum_vpar * 0.25d0 * r0 * Vpar0**2 * BB2 &
                                      * (-(ps0_s * vpar_t - ps0_t * vpar_s)/xjac                   ) / BigR                           &
                                      * (-(ps0_s * v_t    - ps0_t * v_s)   /xjac                   )  * xjac * theta * tstep*tstep    &
                            + tgnum_vpar * 0.25d0 * v * Vpar0**2 * BB2 * (1.d0 - fact_conservative_u) &
                                      * (-(ps0_s * vpar_t - ps0_t * vpar_s)/xjac                   ) / BigR                           &
                                      * (-(ps0_s * r0_t   - ps0_t * r0_s)  /xjac + F0 / BigR * r0_p)  * xjac * theta * tstep*tstep    & 
!=============================== New TG_num terms==================================

                            + tgnum_vpar * 0.75d0 * Vpar * Vpar0**2 * BB2 * fact_conservative_u &
                                      * (-(ps0_s * r0_t - ps0_t * r0_s)/xjac + F0 / BigR * r0_p) / BigR             &
                                      * (-(ps0_s * v_t     - ps0_t * v_s)    /xjac )  * xjac * theta * tstep*tstep  

!===============================End of new TG_num terms============================
  
                    if (normalized_velocity_profile) then
                      amat(var_vpar,var_vpar) = amat(var_vpar,var_vpar) + (visco_par + visco_par_sc_num * tau_sc) * (v_x * Vpar_x + v_y * Vpar_y) * BigR        * xjac  * theta * tstep 
                    else
                      amat(var_vpar,var_vpar) = amat(var_vpar,var_vpar) + (visco_par + visco_par_sc_num * tau_sc) * F0**2 / BigR**2 * (v_x * (Vpar_x - 2*vpar/BigR) + v_y * Vpar_y) * BigR * xjac  * theta * tstep 
                    endif
  
                    amat_k(var_vpar,var_vpar) = - r0 * vpar0 * vpar * BB2 * F0 / BigR * v_p               * xjac * theta * tstep             &
       
                              + tgnum_vpar * 0.5d0 * r0 * Vpar * Vpar0 * BB2 &
                                        * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac + F0 / BigR * vpar0_p) / BigR                     &
                                        * (                                          + F0 / BigR * v_p)  * xjac * theta * tstep*tstep  &
                              + tgnum_vpar * 0.25d0 * r0 * Vpar0**2 * BB2 &
                                        * (-(ps0_s * vpar_t - ps0_t * vpar_s)/xjac                  ) / BigR                           &
                                        * (                      + F0 / BigR * v_p)  * xjac * theta * tstep*tstep    &
!=============================== New TG_num terms==================================

                              + tgnum_vpar * 0.75d0 * Vpar * Vpar0**2 * BB2 * fact_conservative_u &
                                        * (-(ps0_s * r0_t - ps0_t * r0_s)/xjac + F0 / BigR * r0_p) / BigR             &
                                        * (+ F0 / BigR * v_p)  * xjac * theta * tstep*tstep  

!===============================End of new TG_num terms============================
  
                    amat_n(var_vpar,var_vpar) = &
                        ! New term coming from -(\partial_t \rho + \nabla \cdot (\rho \mathbf{v})) \mathbf{v} in RHS of momentum equation
                        ! (see wiki: https://www.jorek.eu/wiki/doku.php?id=model500_501_555#equations):	    
                              + v * r0 * F0 /BigR * vpar_p * vpar0 * BB2 * xjac * theta * tstep * fact_conservative_u &
  
                              + tgnum_vpar * 0.25d0 * r0 * Vpar0**2 * BB2  &
                                        * (                                        + F0 / BigR * vpar_p) / BigR                        &
                                        * (-(ps0_s * v_t    - ps0_t * v_s)   /xjac                     )  * xjac * theta * tstep*tstep &
                              + tgnum_vpar * 0.25d0 * v * Vpar0**2 * BB2 * (1.d0 - fact_conservative_u) &
                                        * (                                        + F0 / BigR * vpar_p) / BigR                        &
                                        * (-(ps0_s * r0_t   - ps0_t * r0_s)  /xjac + F0 / BigR * r0_p)  * xjac * theta * tstep*tstep 
  
                    amat_kn(var_vpar,var_vpar) = &
                               + tgnum_vpar * 0.25d0 * r0 * Vpar0**2 * BB2 &
                                         * (                                        + F0 / BigR * vpar_p) / BigR                        &
                                         * (                                        + F0 / BigR * v_p)  * xjac * theta * tstep*tstep
  
                    if ( NEO ) then
                      amat(var_vpar,var_rho) = amat(var_vpar,var_rho) - v * amu_neo_prof(ms,mt)*BB2/(Btheta2+epsil) &
                            * (rho * (ps0_x*u0_x + ps0_y*u0_y) + tauIC*2.*(ps0_x * (rho_x*Ti0 + rho*Ti0_x) + ps0_y*(rho_y*Ti0 + rho*Ti0_y)) &
                            + aki_neo_prof(ms,mt) * tauIC*2. * rho*(ps0_x*Ti0_x + ps0_y*Ti0_y) - rho * Vpar0 * Btheta2) * BigR * xjac * tstep * theta
                      
                      if ( with_TiTe ) then ! (with_TiTe) ******************************************
                        amat(var_vpar,var_Ti) = amat(var_vpar,var_Ti) -v*amu_neo_prof(ms,mt)*BB2/(Btheta2+epsil)           &
                                  * (tauIC*2. * (ps0_x * (r0_x*Ti + r0*Ti_x) + ps0_y*(r0_y*Ti + r0*Ti_y)) &
                                  + aki_neo_prof(ms,mt) * tauIC*2. * r0 * (ps0_x*Ti_x + ps0_y*Ti_y)) * BigR * xjac * tstep * theta
                      else ! (with_TiTe = .f.), i.e. with single temperature *****************************
                        amat(var_vpar,var_T)  = amat(var_vpar,var_T)  -v*amu_neo_prof(ms,mt)*BB2/(Btheta2+epsil)           &
                                  * (tauIC * (ps0_x * (r0_x*T  + r0*T_x ) + ps0_y*(r0_y*T  + r0*T_y )) &
                                  + aki_neo_prof(ms,mt) * tauIC * r0 * (ps0_x*T_x  + ps0_y*T_y )) * BigR * xjac * tstep * theta
                      end if ! (with_TiTe) *********************************************************
                    
                      amat(var_vpar,var_vpar) = amat(var_vpar,var_vpar) + v * amu_neo_prof(ms,mt) * BB2 * Btheta2/(Btheta2+epsil) * r0 * vpar * BigR * xjac * tstep * theta 
                    endif

                    if (with_neutrals) then
                      amat(var_vpar,var_rhon) = (1.d0 - delta_n_convection) * v *(r0 * rhon * Sion_T) * vpar0 * BB2 * BigR * xjac * theta * tstep
                    endif

                  end if ! (with_vpar)

                  if ( with_TiTe ) then ! (with_TiTe) **********************************************
  
                    !###################################################################################################
                    !#  Ion Energy Equation                                                                            #
                    !###################################################################################################
                    Bgrad_T_star_psi = ( v_x   * psi_y - v_y   * psi_x  ) / BigR
                    Bgrad_Ti_psi     = ( Ti0_x * psi_y - Ti0_y * psi_x )  / BigR
                    Bgrad_Ti_Ti      = ( Ti_x  * ps0_y - Ti_y  * ps0_x )  / BigR
                    Bgrad_Ti_Ti_n    = ( F0 / BigR * Ti_p) / BigR

                    amat(var_Ti,var_psi) = - (ZKi_par_T-ZKi_prof) * BigR * BB2_psi / BB2**2 * Bgrad_T_star     * Bgrad_Ti     * xjac * theta * tstep &
                                + (ZKi_par_T-ZKi_prof) * BigR / BB2              * Bgrad_T_star_psi * Bgrad_Ti     * xjac * theta * tstep &
                                + (ZKi_par_T-ZKi_prof) * BigR / BB2              * Bgrad_T_star     * Bgrad_Ti_psi * xjac * theta * tstep &
  
                              + v * r0  * Vpar0 * (Ti0_s * psi_t - Ti0_t * psi_s)                                   * theta * tstep &
                              + v * Ti0 * Vpar0 * (r0_s * psi_t - r0_t * psi_s)                                     * theta * tstep &
                              + v * r0  * GAMMA * Ti0 * (vpar0_s * psi_t - vpar0_t * psi_s)                         * theta * tstep &
  
                              !===================== Additional terms from friction terms============
                              - v * ((GAMMA - 1.) / BigR) * vpar0**2 * (psi_x * ps0_x + psi_y * ps0_y)&
                                  * (r0*rn0*Sion_T)                                                    * xjac * theta * tstep &
                              !==============================End of friction terms=================
 
                           + tgnum_Ti* 0.25d0 / BigR * vpar0**2                                                         &
                                     * Ti0 * (r0_x * psi_y - r0_y * psi_x)                                              &
                                     * ( v_x * ps0_y -  v_y * ps0_x                  ) * xjac * theta * tstep * tstep   &
                           + tgnum_Ti* 0.25d0 / BigR * vpar0**2                                                         &
                                     * r0 * (Ti0_x * psi_y - Ti0_y * psi_x)                                             &
                                     * ( v_x * ps0_y -  v_y * ps0_x                  ) * xjac * theta * tstep * tstep   & 
                           + tgnum_Ti* 0.25d0 / BigR * vpar0**2                                                         &
                                     * Ti0 * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                           &
                                     * ( v_x * psi_y -  v_y * psi_x ) * xjac * theta * tstep * tstep                    &
                           + tgnum_Ti* 0.25d0 / BigR * vpar0**2                                                         &
                                     * r0 * (Ti0_x * ps0_y - Ti0_y * ps0_x + F0 / BigR * Ti0_p)                         &
                                     * ( v_x * psi_y -  v_y * psi_x ) * xjac * theta * tstep * tstep
  
                    amat_k(var_Ti,var_psi) = - (ZKi_par_T-ZKi_prof) * BigR * BB2_psi / BB2**2 * Bgrad_T_k_star * Bgrad_Ti     * xjac * theta * tstep &
                                  + (ZKi_par_T-ZKi_prof) * BigR / BB2              * Bgrad_T_k_star * Bgrad_Ti_psi * xjac * theta * tstep &
    
                          + tgnum_Ti* 0.25d0 / BigR * vpar0**2                                                       &
                                    * Ti0 * (r0_x * psi_y - r0_y * psi_x)                                            &
                                    * (                            + F0 / BigR * v_p) * xjac * theta * tstep * tstep &
                          + tgnum_Ti* 0.25d0 / BigR * vpar0**2                                                       &
                                    * r0 * (Ti0_x * psi_y - Ti0_y * psi_x)                                           &
                                    * (                            + F0 / BigR * v_p) * xjac * theta * tstep * tstep
  
  
                    amat(var_Ti,var_u) = - v * r0 * BigR**2  * ( Ti0_x * u_y - Ti0_y * u_x)         * xjac * theta * tstep &
                                - v * Ti0 * BigR**2 * ( r0_x * u_y - r0_y * u_x)           * xjac * theta * tstep &
                                - v * r0 * 2.d0* GAMMA * BigR * Ti0 * u_y                  * xjac * theta * tstep &
  
                               !===================== Additional terms from friction terms============
                                - v * BigR**3 * (GAMMA - 1.) * (u_x * u0_x + u_y * u0_y)  &
                                    * (r0*rn0*Sion_T)                                * xjac * theta * tstep &
                               !==============================End of friction terms===================

                           + tgnum_Ti* 0.25d0 * BigR**2 * Ti0* (r0_x * u_y - r0_y * u_x)                &
                                              * ( v_x * u0_y - v_y * u0_x) * xjac * theta*tstep*tstep  &
                           + tgnum_Ti* 0.25d0 * BigR**2 * r0* (Ti0_x * u_y - Ti0_y * u_x)                &
                                              * ( v_x * u0_y - v_y * u0_x) * xjac * theta*tstep*tstep  &
                           + tgnum_Ti* 0.25d0 * BigR**2 * Ti0* (r0_x * u0_y - r0_y * u0_x)              &
                                              * ( v_x * u_y - v_y * u_x) * xjac * theta*tstep*tstep    &
                           + tgnum_Ti* 0.25d0 * BigR**2 * r0* (Ti0_x * u0_y - Ti0_y * u0_x)              &
                                              * ( v_x * u_y - v_y * u_x) * xjac * theta*tstep*tstep 
  
  
                    amat(var_Ti,var_rho) = v * rho * Ti0_corr   * BigR * xjac * (1.d0 + zeta)     &
                              - v * rho * BigR**2 * ( Ti0_s * u0_t - Ti0_t * u0_s)                        * theta * tstep &
                              - v * Ti0 * BigR**2 * ( rho_s * u0_t - rho_t * u0_s)                        * theta * tstep &
                              - v * rho * 2.d0* GAMMA * BigR * Ti0 * u0_y                          * xjac * theta * tstep &
                              + v * rho * F0 / BigR * Vpar0 * Ti0_p                                * xjac * theta * tstep &
  
                              + v * rho * Vpar0 * (Ti0_s  * ps0_t - Ti0_t * ps0_s)                        * theta * tstep &
                              + v * Ti0 * Vpar0 * (rho_s * ps0_t  - rho_t * ps0_s)                        * theta * tstep & 
  
                              + v * rho * GAMMA * Ti0 * (vpar0_s * ps0_t - vpar0_t * ps0_s)        * theta * tstep &
                              + v * rho * GAMMA * Ti0 * F0 / BigR * vpar0_p                 * xjac * theta * tstep &
  
                              ! Energy exchange term
                              - v * BigR * ddTi_e_drho * rho                                * xjac * theta * tstep &

                             !===================== Additional terms from friction terms============
                              - v * BigR * ((GAMMA - 1.)/2.) * vpar0**2 * BB2 * (rho*rn0*Sion_T) * xjac * theta * tstep &
                              - v * BigR * ((GAMMA - 1.)/2.) * vv2            * (rho*rn0*Sion_T) * xjac * theta * tstep &
                             !==============================End of friction terms=================


                           + tgnum_Ti* 0.25d0 * BigR**2 * Ti0* (rho_x * u0_y - rho_y * u0_x)         &
                                     * ( v_x * u0_y - v_y * u0_x) * xjac * theta*tstep*tstep         &
                           + tgnum_Ti* 0.25d0 * BigR**2 * rho * (Ti0_x * u0_y - Ti0_y * u0_x)        &
                                     * ( v_x * u0_y - v_y * u0_x) * xjac* theta*tstep*tstep          &
                           + tgnum_Ti* 0.25d0 / BigR * vpar0**2                                      &
                                     * Ti0 * (rho_x * ps0_y - rho_y * ps0_x )                        &
                                     * ( v_x * ps0_y -  v_y * ps0_x ) * xjac * theta * tstep * tstep &
                           + tgnum_Ti* 0.25d0 / BigR * vpar0**2                                      &
                                     * rho * (Ti0_x * ps0_y - Ti0_y * ps0_x + F0 / BigR * Ti0_p)     &
                                     * ( v_x * ps0_y -  v_y * ps0_x ) * xjac * theta * tstep * tstep
  
                    amat_n(var_Ti,var_rho) = + v * Ti0  * F0 / BigR * Vpar0 * rho_p   * xjac * theta * tstep    &
  
                           + tgnum_Ti* 0.25d0 / BigR * vpar0**2                                                       &
                                     * Ti0 * (                              + F0 / BigR * rho_p)                      &
                                     * ( v_x * ps0_y -  v_y * ps0_x                  ) * xjac * theta * tstep * tstep
  
                    amat_k(var_Ti,var_rho) = + tgnum_Ti* 0.25d0 / BigR * vpar0**2                                                &
                                     * Ti0 * (rho_x * ps0_y - rho_y * ps0_x                    )                      &
                                     * (                            + F0 / BigR * v_p) * xjac * theta * tstep * tstep &
                                + tgnum_Ti* 0.25d0 / BigR * vpar0**2                                                  &
                                     * rho * (Ti0_x * ps0_y - Ti0_y * ps0_x + F0 / BigR * Ti0_p)                      &
                                     * (                            + F0 / BigR * v_p) * xjac * theta * tstep * tstep
  
                    amat_kn(var_Ti,var_rho) = + tgnum_Ti* 0.25d0 / BigR * vpar0**2                 &
                                     * Ti0 * (+ F0 / BigR * rho_p)                      &
                                     * (      + F0 / BigR * v_p) * xjac * theta * tstep * tstep
  
  
                    amat(var_Ti,var_Ti) = v * r0_corr * Ti  * BigR * xjac * (1.d0 + zeta)     &
                              - v * r0 * BigR**2  * ( Ti_s  * u0_t - Ti_t  * u0_s)             * theta * tstep &
                              - v * Ti  * BigR**2 * ( r0_s * u0_t - r0_t * u0_s)               * theta * tstep &
  
                              - v * r0 * 2.d0* GAMMA * BigR * Ti * u0_y                 * xjac * theta * tstep &
  
                              + v * r0 * Vpar0 * (Ti_s  * ps0_t - Ti_t  * ps0_s)               * theta * tstep &
                              + v * Ti * Vpar0 * (r0_s * ps0_t - r0_t * ps0_s)                 * theta * tstep & 
  
                              + v * r0 * GAMMA * Ti * (vpar0_s * ps0_t - vpar0_t * ps0_s)      * theta * tstep &
                              + v * r0 * GAMMA * Ti * F0 / BigR * vpar0_p               * xjac * theta * tstep &
  
                              ! Energy exchange term
                              - v * BigR * ddTi_e_dTi * Ti                              * xjac * theta * tstep &
  
                              + v * Ti * F0 / BigR * Vpar0 * r0_p                       * xjac * theta * tstep &
  
                              + (ZKi_par_T-ZKi_prof) * BigR / BB2 * Bgrad_T_star * Bgrad_Ti_Ti * xjac * theta * tstep &
                              + ZKi_prof * BigR * ( v_x*Ti_x + v_y*Ti_y )                      * xjac * theta * tstep &
  
                              + dZKi_par_dT * Ti * BigR / BB2 * Bgrad_T_star * Bgrad_Ti       * xjac * theta * tstep &
    
                              + ZK_i_perp_num * (v_xx + v_x/BigR + v_yy)*(Ti_xx + Ti_x/BigR + Ti_yy) * BigR * xjac * theta * tstep &
  
  !!!!                            -v * Te * (gamma-1.d0) * deta_dT_ohm * (zj0 / BigR)**2.d0 * BigR * xjac * theta * tstep &
  
                              + tgnum_Ti* 0.25d0 * BigR**2 * Ti* (r0_x * u0_y - r0_y * u0_x)         &
                                        * ( v_x * u0_y - v_y * u0_x) * xjac * theta * tstep * tstep  &
                              + tgnum_Ti* 0.25d0 * BigR**2 * r0* (Ti_x * u0_y - Ti_y * u0_x)         &
                                        * ( v_x * u0_y - v_y * u0_x) * xjac * theta * tstep * tstep  &
                              + tgnum_Ti* 0.25d0 / BigR * vpar0**2                                                       &
                                        * Ti * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                          &
                                        * ( v_x * ps0_y -  v_y * ps0_x                  ) * xjac * theta * tstep * tstep &
                              + tgnum_Ti* 0.25d0 / BigR * vpar0**2                                                       &
                                        * r0 * (Ti_x * ps0_y - Ti_y * ps0_x             )                                &
                                        * ( v_x * ps0_y -  v_y * ps0_x                  ) * xjac * theta * tstep * tstep 
  
                    amat_k(var_Ti,var_Ti) = + (ZKi_par_T-ZKi_prof) * BigR / BB2 * Bgrad_T_k_star * Bgrad_Ti_Ti * xjac * theta * tstep  &
                                  + dZKi_par_dT * Ti     * BigR / BB2 * Bgrad_T_k_star * Bgrad_Ti    * xjac * theta * tstep  &
  
                                + tgnum_Ti* 0.25d0 / BigR * vpar0**2                                                      &
                                    * Ti * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                               &
                                    * (                                 + F0 / BigR * v_p) * xjac * theta * tstep * tstep &
                                + tgnum_Ti* 0.25d0 / BigR * vpar0**2                                                      &
                                    * r0 * (Ti_x * ps0_y - Ti_y * ps0_x                  )                                &
                                    * (                                + F0 / BigR * v_p) * xjac * theta * tstep * tstep
  
                    amat_n(var_Ti,var_Ti) = + (ZKi_par_T-ZKi_prof) * BigR / BB2 * Bgrad_T_star   * Bgrad_Ti_Ti_n  * xjac * theta * tstep &
  
                                + v * r0 * F0 / BigR * Vpar0 * Ti_p                               * xjac * theta * tstep &
    
                                + tgnum_Ti* 0.25d0 / BigR * vpar0**2                       & 
                                  * r0 * ( + F0 / BigR * Ti_p) * ( v_x * ps0_y -  v_y * ps0_x ) * xjac * theta * tstep * tstep
  
                    amat_kn(var_Ti,var_Ti) = + (ZKi_par_T-ZKi_prof) * BigR / BB2 * Bgrad_T_k_star * Bgrad_Ti_Ti_n * xjac * theta * tstep &
                                   + ZKi_prof * BigR   * (v_p*Ti_p /BigR**2 )                           * xjac * theta * tstep &
  
                                + tgnum_Ti* 0.25d0 / BigR * vpar0**2 &
                                  * r0 * ( + F0 / BigR * Ti_p) * ( + F0 / BigR * v_p)          * xjac * theta * tstep * tstep
  
                    if ( with_vpar ) then
                      amat(var_Ti,var_vpar) = + v * r0 * F0 / BigR * Vpar * Ti0_p                            * xjac * theta * tstep &
                                + v * Ti0 * F0 / BigR * Vpar * r0_p                              * xjac * theta * tstep &
    
                                + v * r0  * Vpar * (Ti0_s * ps0_t - Ti0_t * ps0_s)                      * theta * tstep &
                                + v * Ti0 * Vpar * (r0_s * ps0_t - r0_t * ps0_s)                        * theta * tstep & 
    
                                + v * r0 * GAMMA * Ti0 * (vpar_s * ps0_t - vpar_t * ps0_s)              * theta * tstep &

                               !===================== Additional terms from friction terms============
                                - v * BigR *(GAMMA - 1.) * vpar0 * Vpar * BB2 * (r0*rn0*Sion_T) * xjac * theta * tstep &
                               !==============================End of friction terms=================
                               !============================Behold, the parallel viscous heating terms!=============
                               - (GAMMA - 1.) * v * BigR * visco_par_heating * 2.d0 * (vpar_x*vpar0_x + vpar_y*vpar0_y) * xjac * theta * tstep  &
                               - (GAMMA - 1.) * vpar0 * BigR * visco_par_heating    * (vpar_x*v_x     + vpar_y*v_y)     * xjac * theta * tstep  &
                               - (GAMMA - 1.) * vpar * BigR * visco_par_heating    * (vpar0_x*v_x     + vpar0_y*v_y)    * xjac * theta * tstep  &
                               !==========================End of viscous heating terms==============================


                                + tgnum_Ti* 0.25d0 / BigR * 2.d0 * vpar0*vpar &
                                      * Ti0 * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                               &
                                      * ( v_x * ps0_y -  v_y * ps0_x                        ) * xjac * theta * tstep * tstep &
                                + tgnum_Ti* 0.25d0 / BigR * 2.d0 * vpar0*vpar &
                                      * r0 * (Ti0_x * ps0_y - Ti0_y * ps0_x + F0 / BigR * Ti0_p)                             &
                                      * ( v_x * ps0_y -  v_y * ps0_x                        ) * xjac * theta * tstep * tstep 
    
                      amat_k(var_Ti,var_vpar) =  &
                            + tgnum_Ti* 0.25d0 / BigR * 2.d0 * vpar0*vpar &
                                      * Ti0 * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                          &
                                      * (                            + F0 / BigR * v_p) * xjac * theta * tstep * tstep &
                            + tgnum_Ti* 0.25d0 / BigR * 2.d0 * vpar0*vpar &
                                      * r0 * (Ti0_x * ps0_y - Ti0_y * ps0_x + F0 / BigR * Ti0_p)                          &
                                      * (                            + F0 / BigR * v_p) * xjac * theta * tstep * tstep 
    
                      amat_n(var_Ti,var_vpar) = + v * r0 * GAMMA * Ti0 * F0 / BigR * vpar_p          * xjac * theta * tstep
                    end if ! (with_vpar)
  
                    amat(var_Ti,var_Te) = - v * BigR * ddTi_e_dTe * Te                           * xjac * theta * tstep &

                                          !===================== Additional terms from friction terms============
                                          - v * BigR * ((GAMMA - 1.)/2.) * vpar0**2 * BB2 &
                                              * (r0*rn0*dSion_dT) * Te * xjac * theta * tstep &
                                          - v * BigR * ((GAMMA - 1.)/2.) * vv2 &
                                              * (r0*rn0*dSion_dT) * Te * xjac * theta * tstep 
                                          !==============================End of friction terms=================

                    if (with_neutrals) then
                      !===================== Additional terms from friction terms============
                      amat(var_Ti,var_rhon) = - v * BigR * ((GAMMA - 1.)/2.) * vpar0**2 * BB2 * (r0*rhon*Sion_T) * xjac * theta * tstep &
                                              - v * BigR * ((GAMMA - 1.)/2.) * vv2            * (r0*rhon*Sion_T) * xjac * theta * tstep 
                      !==============================End of friction terms=================

                    endif
  
                    !###################################################################################################
                    !#  Electron Energy Equation                                                                       #
                    !###################################################################################################
                    
                    Bgrad_T_star_psi = ( v_x   * psi_y - v_y   * psi_x  ) / BigR
                    Bgrad_Te_psi     = ( Te0_x * psi_y - Te0_y * psi_x )  / BigR
                    Bgrad_Te_Te      = ( Te_x  * ps0_y - Te_y  * ps0_x )  / BigR
                    Bgrad_Te_Te_n    = ( F0 / BigR * Te_p) / BigR
                    
                    amat(var_Te,var_psi) = - (ZKe_par_T-ZKe_prof) * BigR * BB2_psi / BB2**2 * Bgrad_T_star     * Bgrad_Te     * xjac * theta * tstep &
                                + (ZKe_par_T-ZKe_prof) * BigR / BB2              * Bgrad_T_star_psi * Bgrad_Te     * xjac * theta * tstep &
                                + (ZKe_par_T-ZKe_prof) * BigR / BB2              * Bgrad_T_star     * Bgrad_Te_psi * xjac * theta * tstep &
  
                              + v * r0  * Vpar0 * (Te0_s * psi_t - Te0_t * psi_s)                                   * theta * tstep &
                              + v * Te0 * Vpar0 * (r0_s  * psi_t - r0_t * psi_s)                                    * theta * tstep &
                              + v * r0  * GAMMA * Te0 * (vpar0_s * psi_t - vpar0_t * psi_s)                         * theta * tstep &
  
                           + tgnum_Te * 0.25d0 / BigR * vpar0**2                                                         &
                                     * Te0 * (r0_x * psi_y - r0_y * psi_x)                                              &
                                     * ( v_x * ps0_y -  v_y * ps0_x                  ) * xjac * theta * tstep * tstep   &
                           + tgnum_Te * 0.25d0 / BigR * vpar0**2                                                         &
                                     * r0 * (Te0_x * psi_y - Te0_y * psi_x)                                             &
                                     * ( v_x * ps0_y -  v_y * ps0_x                  ) * xjac * theta * tstep * tstep   & 
                           + tgnum_Te * 0.25d0 / BigR * vpar0**2                                                         &
                                     * Te0 * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                           &
                                     * ( v_x * psi_y -  v_y * psi_x ) * xjac * theta * tstep * tstep                    &
                           + tgnum_Te * 0.25d0 / BigR * vpar0**2                                                         &
                                     * r0 * (Te0_x * ps0_y - Te0_y * ps0_x + F0 / BigR * Te0_p)                         &
                                     * ( v_x * psi_y -  v_y * psi_x ) * xjac * theta * tstep * tstep
  
                    amat_k(var_Te,var_psi) = - (ZKe_par_T-ZKe_prof) * BigR * BB2_psi / BB2**2 * Bgrad_T_k_star * Bgrad_Te     * xjac * theta * tstep &
                                  + (ZKe_par_T-ZKe_prof) * BigR / BB2              * Bgrad_T_k_star * Bgrad_Te_psi * xjac * theta * tstep &
    
                          + tgnum_Te * 0.25d0 / BigR * vpar0**2                                                       &
                                    * Te0 * (r0_x * psi_y - r0_y * psi_x)                                            &
                                    * (                            + F0 / BigR * v_p) * xjac * theta * tstep * tstep &
                          + tgnum_Te * 0.25d0 / BigR * vpar0**2                                                       &
                                    * r0 * (Te0_x * psi_y - Te0_y * psi_x)                                           &
                                    * (                            + F0 / BigR * v_p) * xjac * theta * tstep * tstep
  
  
                    amat(var_Te,var_u) = - v * r0 * BigR**2  * ( Te0_x * u_y - Te0_y * u_x)         * xjac * theta * tstep &
                                - v * Te0 * BigR**2 * ( r0_x * u_y - r0_y * u_x)           * xjac * theta * tstep &
                                - v * r0 * 2.d0* GAMMA * BigR * Te0 * u_y                  * xjac * theta * tstep &
  
                           + tgnum_Te * 0.25d0 * BigR**2 * Te0* (r0_x * u_y - r0_y * u_x)               &
                                              * ( v_x * u0_y - v_y * u0_x) * xjac * theta*tstep*tstep  &
                           + tgnum_Te * 0.25d0 * BigR**2 * r0* (Te0_x * u_y - Te0_y * u_x)              &
                                              * ( v_x * u0_y - v_y * u0_x) * xjac * theta*tstep*tstep  &
                           + tgnum_Te * 0.25d0 * BigR**2 * Te0* (r0_x * u0_y - r0_y * u0_x)             &
                                              * ( v_x * u_y - v_y * u_x) * xjac * theta*tstep*tstep    &
                           + tgnum_Te * 0.25d0 * BigR**2 * r0* (Te0_x * u0_y - Te0_y * u0_x)            &
                                              * ( v_x * u_y - v_y * u_x) * xjac * theta*tstep*tstep 
  
                    amat(var_Te,var_zj) = - v * (gamma-1.d0) * eta_T_ohm * 2.d0 * zj * zj0/(BigR**2.d0) * BigR * xjac * theta * tstep
  
                    amat(var_Te,var_rho) = v * rho * Te0_corr   * BigR * xjac * (1.d0 + zeta)     &
                              - v * rho * BigR**2 * ( Te0_s * u0_t - Te0_t * u0_s)                        * theta * tstep &
                              - v * Te0 * BigR**2 * ( rho_s * u0_t - rho_t * u0_s)                        * theta * tstep &
                              - v * rho * 2.d0* GAMMA * BigR * Te0 * u0_y                          * xjac * theta * tstep &
                              + v * rho * F0 / BigR * Vpar0 * Te0_p                                * xjac * theta * tstep &
  
                              + v * rho * Vpar0 * (Te0_s  * ps0_t - Te0_t * ps0_s)                        * theta * tstep &
                              + v * Te0 * Vpar0 * (rho_s * ps0_t  - rho_t * ps0_s)                        * theta * tstep & 
  
                              + v * rho * GAMMA * Te0 * (vpar0_s * ps0_t - vpar0_t * ps0_s)        * theta * tstep &
                              + v * rho * GAMMA * Te0 * F0 / BigR * vpar0_p                 * xjac * theta * tstep &
                              ! Energy exchange term
                              - v * BigR * ddTe_i_drho * rho                                * xjac * theta * tstep &

                              + v * BigR * rho * rn0 * ksiion * Sion_T                             * xjac * theta * tstep &
                              + v * BigR * rho * rn0_corr * LradDrays_T                            * xjac * theta * tstep &
                              + v * BigR * rho * 2d0 * r0_corr * LradDcont_T                       * xjac * theta * tstep &
                              + v * BigR * rho * frad_bg                                           * xjac * theta * tstep &  
  
                           + tgnum_Te * 0.25d0 * BigR**2 * Te0* (rho_x * u0_y - rho_y * u0_x)         &
                                     * ( v_x * u0_y - v_y * u0_x) * xjac * theta*tstep*tstep         &
                           + tgnum_Te * 0.25d0 * BigR**2 * rho * (Te0_x * u0_y - Te0_y * u0_x)        &
                                     * ( v_x * u0_y - v_y * u0_x) * xjac* theta*tstep*tstep          &
                           + tgnum_Te * 0.25d0 / BigR * vpar0**2                                      &
                                     * Te0 * (rho_x * ps0_y - rho_y * ps0_x )                        &
                                     * ( v_x * ps0_y -  v_y * ps0_x ) * xjac * theta * tstep * tstep &
                           + tgnum_Te * 0.25d0 / BigR * vpar0**2                                      &
                                     * rho * (Te0_x * ps0_y - Te0_y * ps0_x + F0 / BigR * Te0_p)     &
                                     * ( v_x * ps0_y -  v_y * ps0_x ) * xjac * theta * tstep * tstep
  
                    amat_n(var_Te,var_rho) = + v * Te0  * F0 / BigR * Vpar0 * rho_p   * xjac * theta * tstep    &
  
                           + tgnum_Te * 0.25d0 / BigR * vpar0**2                                                       &
                                     * Te0 * (                              + F0 / BigR * rho_p)                      &
                                     * ( v_x * ps0_y -  v_y * ps0_x                  ) * xjac * theta * tstep * tstep
  
                    amat_k(var_Te,var_rho) = + tgnum_Te * 0.25d0 / BigR * vpar0**2                                                &
                                     * Te0 * (rho_x * ps0_y - rho_y * ps0_x                    )                      &
                                     * (                            + F0 / BigR * v_p) * xjac * theta * tstep * tstep &
                                + tgnum_Te * 0.25d0 / BigR * vpar0**2                                                  &
                                     * rho * (Te0_x * ps0_y - Te0_y * ps0_x + F0 / BigR * Te0_p)                      &
                                     * (                            + F0 / BigR * v_p) * xjac * theta * tstep * tstep
  
                    amat_kn(var_Te,var_rho) = + tgnum_Te * 0.25d0 / BigR * vpar0**2                 &
                                     * Te0 * (+ F0 / BigR * rho_p)                      &
                                     * (      + F0 / BigR * v_p) * xjac * theta * tstep * tstep
  
                    amat(var_Te,var_Ti) = - v * BigR * ddTe_i_dTi * Ti                            * xjac * theta * tstep
  
  
                    amat(var_Te,var_Te) = v * r0_corr * Te  * BigR * xjac * (1.d0 + zeta)     &
                              - v * r0 * BigR**2  * ( Te_s * u0_t - Te_t  * u0_s)              * theta * tstep &
                              - v * Te  * BigR**2 * ( r0_s * u0_t - r0_t * u0_s)               * theta * tstep &
  
                              - v * r0 * 2.d0* GAMMA * BigR * Te * u0_y                 * xjac * theta * tstep &
  
                              + v * r0 * Vpar0 * (Te_s * ps0_t - Te_t  * ps0_s)                * theta * tstep &
                              + v * Te * Vpar0 * (r0_s * ps0_t - r0_t * ps0_s)                 * theta * tstep & 
  
                              + v * r0 * GAMMA * Te * (vpar0_s * ps0_t - vpar0_t * ps0_s)      * theta * tstep &
                              + v * r0 * GAMMA * Te * F0 / BigR * vpar0_p               * xjac * theta * tstep &
  
                              + v * Te * F0 / BigR * Vpar0 * r0_p                       * xjac * theta * tstep &
  
                              ! Energy exchange term
                              - v * BigR * ddTe_i_dTe * Te                              * xjac * theta * tstep &
  
                              + (ZKe_par_T-ZKe_prof) * BigR / BB2 * Bgrad_T_star * Bgrad_Te_Te * xjac * theta * tstep &
                              + ZKe_prof * BigR * ( v_x*Te_x + v_y*Te_y )                      * xjac * theta * tstep &
  
                              + dZKe_par_dT * Te * BigR / BB2 * Bgrad_T_star * Bgrad_Te       * xjac * theta * tstep &
    
                              + ZK_e_perp_num * (v_xx + v_x/BigR + v_yy)*(Te_xx + Te_x/BigR + Te_yy) * BigR * xjac * theta * tstep &
  
                              - v * Te * (gamma-1.d0) * deta_dT_ohm * (zj0 / BigR)**2.d0 * BigR * xjac * theta * tstep &

                              + v * BigR * r0 * rn0 * ksiion * dSion_dT * T                   * xjac * theta * tstep &
                              + v * BigR * T * r0_corr * rn0_corr * dLradDrays_dT             * xjac * theta * tstep &
                              + v * BigR * T * r0_corr * r0_corr  * dLradDcont_dT             * xjac * theta * tstep &
                              + v * BigR * T * r0_corr * dfrad_bg_dT                          * xjac * theta * tstep &

  
                              + tgnum_Te * 0.25d0 * BigR**2 * Te* (r0_x * u0_y - r0_y * u0_x)         &
                                        * ( v_x * u0_y - v_y * u0_x) * xjac * theta * tstep * tstep  &
                              + tgnum_Te * 0.25d0 * BigR**2 * r0* (Te_x * u0_y - Te_y * u0_x)         &
                                        * ( v_x * u0_y - v_y * u0_x) * xjac * theta * tstep * tstep  &
                              + tgnum_Te * 0.25d0 / BigR * vpar0**2                                                       &
                                        * Te * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                          &
                                        * ( v_x * ps0_y -  v_y * ps0_x                  ) * xjac * theta * tstep * tstep &
                              + tgnum_Te * 0.25d0 / BigR * vpar0**2                                                       &
                                        * r0 * (Te_x * ps0_y - Te_y * ps0_x             )                                &
                                        * ( v_x * ps0_y -  v_y * ps0_x                  ) * xjac * theta * tstep * tstep 
  
                    amat_k(var_Te,var_Te) = + (ZKe_par_T-ZKe_prof) * BigR / BB2 * Bgrad_T_k_star * Bgrad_Te_Te * xjac * theta * tstep  &
                                  + dZKe_par_dT * Te     * BigR / BB2 * Bgrad_T_k_star * Bgrad_Te    * xjac * theta * tstep  &
  
                                + tgnum_Te * 0.25d0 / BigR * vpar0**2                                                      &
                                    * Te * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                               &
                                    * (                                 + F0 / BigR * v_p) * xjac * theta * tstep * tstep &
                                + tgnum_Te * 0.25d0 / BigR * vpar0**2                                                      &
                                    * r0 * (Te_x * ps0_y - Te_y * ps0_x                  )                                &
                                    * (                                + F0 / BigR * v_p) * xjac * theta * tstep * tstep
  
                    amat_n(var_Te,var_Te) = + (ZKe_par_T-ZKe_prof) * BigR / BB2 * Bgrad_T_star   * Bgrad_Te_Te_n  * xjac * theta * tstep &
  
                                + v * r0 * F0 / BigR * Vpar0 * Te_p                               * xjac * theta * tstep &
    
                                + tgnum_Te * 0.25d0 / BigR * vpar0**2                       & 
                                  * r0 * ( + F0 / BigR * Te_p) * ( v_x * ps0_y -  v_y * ps0_x ) * xjac * theta * tstep * tstep
  
                    amat_kn(var_Te,var_Te) = + (ZKe_par_T-ZKe_prof) * BigR / BB2 * Bgrad_T_k_star * Bgrad_Te_Te_n * xjac * theta * tstep &
                                   + ZKe_prof * BigR   * (v_p*Te_p /BigR**2 )                           * xjac * theta * tstep &
  
                                + tgnum_Te * 0.25d0 / BigR * vpar0**2 &
                                  * r0 * ( + F0 / BigR * Te_p) * ( + F0 / BigR * v_p)           * xjac * theta * tstep * tstep
   
                    if ( with_vpar ) then
                      amat(var_Te,var_vpar) = + v * r0 * F0 / BigR * Vpar * Te0_p                             * xjac * theta * tstep &
                                  + v * Te0 * F0 / BigR * Vpar * r0_p                             * xjac * theta * tstep &
    
                                  + v * r0  * Vpar * (Te0_s * ps0_t - Te0_t * ps0_s)                      * theta * tstep &
                                  + v * Te0 * Vpar * (r0_s * ps0_t  - r0_t * ps0_s)                       * theta * tstep & 
    
                                  + v * r0 * GAMMA * Te0 * (vpar_s * ps0_t - vpar_t * ps0_s)              * theta * tstep &
      
                                  + tgnum_Te * 0.25d0 / BigR * 2.d0 * vpar0*vpar &
                                      * Te0 * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                               &
                                      * ( v_x * ps0_y -  v_y * ps0_x                        ) * xjac * theta * tstep * tstep &
                                  + tgnum_Te * 0.25d0 / BigR * 2.d0 * vpar0*vpar &
                                      * r0 * (Te0_x * ps0_y - Te0_y * ps0_x + F0 / BigR * Te0_p)                             &
                                      * ( v_x * ps0_y -  v_y * ps0_x                        ) * xjac * theta * tstep * tstep 
    
                      amat_k(var_Te,var_vpar) =  &
                            + tgnum_Te * 0.25d0 / BigR * 2.d0 * vpar0*vpar &
                                      * Te0 * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                          &
                                      * (                            + F0 / BigR * v_p) * xjac * theta * tstep * tstep &
                            + tgnum_Te * 0.25d0 / BigR * 2.d0 * vpar0*vpar &
                                      * r0 * (Te0_x * ps0_y - Te0_y * ps0_x + F0 / BigR * Te0_p)                          &
                                      * (                            + F0 / BigR * v_p) * xjac * theta * tstep * tstep 
    
                      amat_n(var_Te,var_vpar) = + v * r0 * GAMMA * Te0 * F0 / BigR * vpar_p          * xjac * theta * tstep
                    end if ! (with_vpar)
                    if (with_neutrals) then
                      amat(var_Te,var_rhon) = + v * BigR * r0 * rhon * ksiion * Sion_T               * xjac * theta * tstep &
                                              + v * BigR * rhon * r0 * LradDrays_T                   * xjac * theta * tstep 
                    endif
                    
                  else ! (with_TiTe = .f.), i.e. with single temperature *********************************
   
                    !###################################################################################################
                    !#  Electron + Ion Energy Equation                                                                 #
                    !###################################################################################################
                    
                    Bgrad_T_star_psi = ( v_x   * psi_y - v_y   * psi_x  ) / BigR
                    Bgrad_T_psi      = ( T0_x  * psi_y - T0_y  * psi_x )  / BigR
                    Bgrad_T_T        = ( T_x   * ps0_y - T_y   * ps0_x )  / BigR
                    Bgrad_T_T_n      = ( F0 / BigR * T_p) / BigR
                    
                    amat(var_T,var_psi) = - (ZK_par_T-ZK_prof) * BigR * BB2_psi / BB2**2 * Bgrad_T_star     * Bgrad_T     * xjac * theta * tstep &
                                          + (ZK_par_T-ZK_prof) * BigR / BB2              * Bgrad_T_star_psi * Bgrad_T     * xjac * theta * tstep &
                                          + (ZK_par_T-ZK_prof) * BigR / BB2              * Bgrad_T_star     * Bgrad_T_psi * xjac * theta * tstep &
  
                                          + v * r0  * Vpar0 * (T0_s * psi_t - T0_t * psi_s)                                      * theta * tstep &
                                          + v * T0 * Vpar0 * (r0_s  * psi_t - r0_t * psi_s)                                      * theta * tstep &
                                          + v * r0  * GAMMA * T0 * (vpar0_s * psi_t - vpar0_t * psi_s)                           * theta * tstep &
                                         !===================== Additional terms from friction terms============
                                          - v * ((GAMMA - 1.) / BigR) * vpar0**2 * (psi_x * ps0_x + psi_y * ps0_y)&
                                              * (r0*rn0*Sion_T)                                                           * xjac * theta * tstep &
                                         !==============================End of friction terms=================
  
                          + tgnum_T * 0.25d0 / BigR * vpar0**2                                                        &
                                    * T0 * (r0_x * psi_y - r0_y * psi_x)                                              &
                                    * ( v_x * ps0_y -  v_y * ps0_x                  ) * xjac * theta * tstep * tstep  &
                          + tgnum_T * 0.25d0 / BigR * vpar0**2                                                        &
                                    * r0 * (T0_x * psi_y - T0_y * psi_x)                                              &
                                    * ( v_x * ps0_y -  v_y * ps0_x                  ) * xjac * theta * tstep * tstep  &
                          + tgnum_T * 0.25d0 / BigR * vpar0**2                                                        &
                                    * T0 * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                           &
                                    * ( v_x * psi_y -  v_y * psi_x ) * xjac * theta * tstep * tstep                   &
                          + tgnum_T * 0.25d0 / BigR * vpar0**2                                                        &
                                    * r0 * (T0_x * ps0_y - T0_y * ps0_x + F0 / BigR * T0_p)                           &
                                    * ( v_x * psi_y -  v_y * psi_x ) * xjac * theta * tstep * tstep
  
                    amat_k(var_T,var_psi) = - (ZK_par_T-ZK_prof) * BigR * BB2_psi / BB2**2 * Bgrad_T_k_star * Bgrad_T     * xjac * theta * tstep &
                                            + (ZK_par_T-ZK_prof) * BigR / BB2              * Bgrad_T_k_star * Bgrad_T_psi * xjac * theta * tstep &
    
                          + tgnum_T * 0.25d0 / BigR * vpar0**2                                                        &
                                    * T0 * (r0_x * psi_y - r0_y * psi_x)                                              &
                                    * (                            + F0 / BigR * v_p)  * xjac * theta * tstep * tstep &
                          + tgnum_T * 0.25d0 / BigR * vpar0**2                                                        &
                                    * r0 * (T0_x * psi_y - T0_y * psi_x)                                              &
                                    * (                            + F0 / BigR * v_p)  * xjac * theta * tstep * tstep
  
  
                    amat(var_T,var_u) = - v * r0 * BigR**2 * ( T0_x * u_y - T0_y * u_x)           * xjac * theta * tstep &
                                        - v * T0 * BigR**2 * ( r0_x * u_y - r0_y * u_x)           * xjac * theta * tstep &
                                        - v * r0 * 2.d0* GAMMA * BigR * T0 * u_y                  * xjac * theta * tstep &
                                       !===================== Additional terms from friction terms============
                                        - v * BigR**3 * (GAMMA - 1.) * (u_x * u0_x + u_y * u0_y)  &
                                            * (r0*rn0*Sion_T)                                     * xjac * theta * tstep &
                                       !==============================End of friction terms===================
  
                          + tgnum_T * 0.25d0 * BigR**2 * T0* (r0_x * u_y - r0_y * u_x)                                &
                                             * ( v_x * u0_y - v_y * u0_x)              * xjac * theta * tstep * tstep &
                          + tgnum_T * 0.25d0 * BigR**2 * r0* (T0_x * u_y - T0_y * u_x)                                &
                                             * ( v_x * u0_y - v_y * u0_x)              * xjac * theta * tstep * tstep &
                          + tgnum_T * 0.25d0 * BigR**2 * T0* (r0_x * u0_y - r0_y * u0_x)                              &
                                             * ( v_x * u_y - v_y * u_x)                * xjac * theta * tstep * tstep &
                          + tgnum_T * 0.25d0 * BigR**2 * r0* (T0_x * u0_y - T0_y * u0_x)                              &
                                             * ( v_x * u_y - v_y * u_x)                * xjac * theta * tstep * tstep 
  
                    amat(var_T,var_zj)  = - v * (gamma-1.d0) * eta_T_ohm * 2.d0 * zj * zj0/(BigR**2.d0) * BigR * xjac * theta * tstep
  
                    amat(var_T,var_rho) =   v * rho * T0_corr   * BigR * xjac * (1.d0 + zeta)     &
                                          - v * rho * BigR**2 * ( T0_s  * u0_t - T0_t  * u0_s)                        * theta * tstep &
                                          - v * T0  * BigR**2 * ( rho_s * u0_t - rho_t * u0_s)                        * theta * tstep &
                                          - v * rho * 2.d0* GAMMA * BigR * T0 * u0_y                           * xjac * theta * tstep &
                                          + v * rho * F0 / BigR * Vpar0 * T0_p                                 * xjac * theta * tstep &
                                          
                                          + v * rho * Vpar0 * (T0_s  * ps0_t  - T0_t  * ps0_s)                        * theta * tstep &
                                          + v * T0  * Vpar0 * (rho_s * ps0_t  - rho_t * ps0_s)                        * theta * tstep & 
                                          
                                          + v * rho * GAMMA * T0 * (vpar0_s * ps0_t - vpar0_t * ps0_s)                * theta * tstep &
                                          + v * rho * GAMMA * T0 * F0 / BigR * vpar0_p                         * xjac * theta * tstep &
  
                                          + v * BigR * rho * rn0 * ksiion * Sion_T                             * xjac * theta * tstep &
                                          + v * BigR * rho * rn0_corr * LradDrays_T                            * xjac * theta * tstep &
                                          + v * BigR * rho * 2d0 * r0_corr * LradDcont_T                       * xjac * theta * tstep &
                                          + v * BigR * rho * frad_bg                                           * xjac * theta * tstep &

                          !===================== Additional terms from friction terms============
                          - v * BigR * ((GAMMA - 1.)/2.) * vpar0**2 * BB2 * (rho*rn0*Sion_T) * xjac * theta * tstep &
                          - v * BigR * ((GAMMA - 1.)/2.) * vv2            * (rho*rn0*Sion_T) * xjac * theta * tstep &
                          !==============================End of friction terms=================

                          + tgnum_T * 0.25d0 * BigR**2 * T0* (rho_x * u0_y - rho_y * u0_x)                            &
                                    * ( v_x * u0_y - v_y * u0_x)                       * xjac * theta * tstep * tstep &
                          + tgnum_T * 0.25d0 * BigR**2 * rho * (T0_x * u0_y - T0_y * u0_x)                            &
                                    * ( v_x * u0_y - v_y * u0_x)                       * xjac * theta * tstep * tstep &
                          + tgnum_T * 0.25d0 / BigR * vpar0**2                                                        &
                                    * T0 * (rho_x * ps0_y - rho_y * ps0_x )                                           &
                                    * ( v_x * ps0_y -  v_y * ps0_x )                   * xjac * theta * tstep * tstep &
                          + tgnum_T * 0.25d0 / BigR * vpar0**2                                                        &
                                    * rho * (T0_x * ps0_y - T0_y * ps0_x + F0 / BigR * T0_p)                          &
                                    * ( v_x * ps0_y -  v_y * ps0_x )                   * xjac * theta * tstep * tstep
  
                    amat_n(var_T,var_rho) = + v * T0  * F0 / BigR * Vpar0 * rho_p      * xjac * theta * tstep         &
  
                          + tgnum_T * 0.25d0 / BigR * vpar0**2                                                        &
                                    * T0 * (                              + F0 / BigR * rho_p)                        &
                                    * ( v_x * ps0_y -  v_y * ps0_x                  )  * xjac * theta * tstep * tstep
  
                    amat_k(var_T,var_rho) =                                                                           &
                          + tgnum_T * 0.25d0 / BigR * vpar0**2                                                        &
                                    * T0 * (rho_x * ps0_y - rho_y * ps0_x                    )                        &
                                    * (                            + F0 / BigR * v_p)  * xjac * theta * tstep * tstep &
                          + tgnum_T * 0.25d0 / BigR * vpar0**2                                                        &
                                    * rho * (T0_x * ps0_y - T0_y * ps0_x + F0 / BigR * T0_p)                          &
                                    * (                            + F0 / BigR * v_p)  * xjac * theta * tstep * tstep
  
                    amat_kn(var_T,var_rho) =                                                                          &
                          + tgnum_T * 0.25d0 / BigR * vpar0**2                                                        &
                                    * T0 * (+ F0 / BigR * rho_p)                                                      &
                                    * (     + F0 / BigR * v_p)                         * xjac * theta * tstep * tstep
  
                    amat(var_T,var_T) = v * r0_corr * T  * BigR * xjac * (1.d0 + zeta)     &
                                      - v * r0 * BigR**2 * ( T_s  * u0_t - T_t  * u0_s)               * theta * tstep &
                                      - v * T  * BigR**2 * ( r0_s * u0_t - r0_t * u0_s)               * theta * tstep &
                                      
                                      - v * r0 * 2.d0* GAMMA * BigR * T * u0_y                 * xjac * theta * tstep &
                                      
                                      + v * r0 * Vpar0 * (T_s  * ps0_t - T_t  * ps0_s)                * theta * tstep &
                                      + v * T  * Vpar0 * (r0_s * ps0_t - r0_t * ps0_s)                * theta * tstep & 
                                      
                                      + v * r0 * GAMMA * T * (vpar0_s * ps0_t - vpar0_t * ps0_s)      * theta * tstep &
                                      + v * r0 * GAMMA * T * F0 / BigR * vpar0_p               * xjac * theta * tstep &
                                      
                                      + v * T * F0 / BigR * Vpar0 * r0_p                       * xjac * theta * tstep &
                                      
                                      + (ZK_par_T-ZK_prof) * BigR / BB2 * Bgrad_T_star * Bgrad_T_T * xjac * theta * tstep &
                                      + ZK_prof * BigR * ( v_x*T_x + v_y*T_y )                      *xjac * theta * tstep &
                                      
                                      + dZK_par_dT * T * BigR / BB2 * Bgrad_T_star * Bgrad_T       * xjac * theta * tstep &
                                      
                                      + ZK_perp_num * (v_xx + v_x/BigR + v_yy)*(T_xx + T_x/BigR + T_yy) * BigR * xjac * theta * tstep &
                                      
                                      - v * T * (gamma-1.d0) * deta_dT_ohm * (zj0 / BigR)**2.d0         * BigR * xjac * theta * tstep &
  
                                      + v * BigR * r0 * rn0           * ksiion * dSion_dT * T  * xjac * theta * tstep &
                                      + v * BigR * T * r0_corr * rn0_corr * dLradDrays_dT      * xjac * theta * tstep &
                                      + v * BigR * T * r0_corr * r0_corr  * dLradDcont_dT      * xjac * theta * tstep &
                                      + v * BigR * T * r0_corr * dfrad_bg_dT                   * xjac * theta * tstep &

                                     !===================== Additional terms from friction terms============
                                      - v * BigR * ((GAMMA - 1.)/2.) * vpar0**2 * BB2 &
                                          * (r0*rn0*dSion_dT) * T * xjac * theta * tstep &
                                      - v * BigR * ((GAMMA - 1.)/2.) * vv2 &
                                          * (r0*rn0*dSion_dT) * T * xjac * theta * tstep &
                                     !==============================End of friction terms=================

                          + tgnum_T * 0.25d0 * BigR**2 * T* (r0_x * u0_y - r0_y * u0_x)                               &
                                    * ( v_x * u0_y - v_y * u0_x)                       * xjac * theta * tstep * tstep &
                          + tgnum_T * 0.25d0 * BigR**2 * r0* (T_x * u0_y - T_y * u0_x)                                &
                                    * ( v_x * u0_y - v_y * u0_x)                       * xjac * theta * tstep * tstep &
                          + tgnum_T * 0.25d0 / BigR * vpar0**2                                                        &
                                    * T * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                            &
                                    * ( v_x * ps0_y -  v_y * ps0_x                  )  * xjac * theta * tstep * tstep &
                          + tgnum_T * 0.25d0 / BigR * vpar0**2                                                        &
                                    * r0 * (T_x * ps0_y - T_y * ps0_x               )                                 &
                                    * ( v_x * ps0_y -  v_y * ps0_x                  )  * xjac * theta * tstep * tstep 
  
                    amat_k(var_T,var_T) = + (ZK_par_T-ZK_prof) * BigR / BB2 * Bgrad_T_k_star * Bgrad_T_T * xjac * theta * tstep  &
                                          + dZK_par_dT * T     * BigR / BB2 * Bgrad_T_k_star * Bgrad_T   * xjac * theta * tstep  &
  
                          + tgnum_T * 0.25d0 / BigR * vpar0**2                                                            &
                                    * T * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                                &
                                    * (                                 + F0 / BigR * v_p) * xjac * theta * tstep * tstep &
                          + tgnum_T * 0.25d0 / BigR * vpar0**2                                                            &
                                    * r0 * (T_x * ps0_y - T_y * ps0_x                  )                                  &
                                    * (                                + F0 / BigR * v_p)  * xjac * theta * tstep * tstep
  
                    amat_n(var_T,var_T) = + (ZK_par_T-ZK_prof) * BigR / BB2 * Bgrad_T_star   * Bgrad_T_T_n  * xjac * theta * tstep &
  
                                          + v * r0 * F0 / BigR * Vpar0 * T_p                                * xjac * theta * tstep &
    
                          + tgnum_T * 0.25d0 / BigR * vpar0**2                                                        &
                                    * r0 * ( + F0 / BigR * T_p) * ( v_x * ps0_y - v_y * ps0_x ) * xjac * theta * tstep * tstep
  
                    amat_kn(var_T,var_T) = + (ZK_par_T-ZK_prof) * BigR / BB2 * Bgrad_T_k_star * Bgrad_T_T_n * xjac * theta * tstep &
                                           + ZK_prof * BigR   * (v_p*T_p /BigR**2 )                         * xjac * theta * tstep &
  
                          + tgnum_T * 0.25d0 / BigR * vpar0**2                                                        &
                                    * r0 * ( + F0 / BigR * T_p) * ( + F0 / BigR * v_p) * xjac * theta * tstep * tstep
   
                    if ( with_vpar ) then
                      amat(var_T,var_vpar) = + v * r0 * F0 / BigR * Vpar * T0_p                             * xjac * theta * tstep &
                                             + v * T0 * F0 / BigR * Vpar * r0_p                             * xjac * theta * tstep &
                                             
                                             + v * r0 * Vpar * (T0_s * ps0_t - T0_t * ps0_s)                       * theta * tstep &
                                             + v * T0 * Vpar * (r0_s * ps0_t - r0_t * ps0_s)                       * theta * tstep & 
                                             
                                             + v * r0 * GAMMA * T0 * (vpar_s * ps0_t - vpar_t * ps0_s)             * theta * tstep &

                            !===================== Additional terms from friction terms============
                            - v * BigR *(GAMMA - 1.) * vpar0 * Vpar * BB2 * (r0*rn0*Sion_T) * xjac * theta * tstep &
                            !==============================End of friction terms=================
                            !============================Behold, the parallel viscous heating terms!=============
                            - (GAMMA - 1.) * v * BigR * visco_par_heating * 2.d0 * (vpar_x*vpar0_x + vpar_y*vpar0_y) * xjac * theta * tstep  &
                            - (GAMMA - 1.) * vpar0 * BigR * visco_par_heating    * (vpar_x*v_x     + vpar_y*v_y)     * xjac * theta * tstep  &
                            - (GAMMA - 1.) * vpar * BigR * visco_par_heating    * (vpar0_x*v_x     + vpar0_y*v_y)    * xjac * theta * tstep  &
                            !==========================End of viscous heating terms==============================

     
                          + tgnum_T * 0.25d0 / BigR * 2.d0 * vpar0*vpar                                               &
                                    * T0 * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                           &
                                    * ( v_x * ps0_y -  v_y * ps0_x                   ) * xjac * theta * tstep * tstep &
                          + tgnum_T * 0.25d0 / BigR * 2.d0 * vpar0*vpar                                               &
                                    * r0 * (T0_x * ps0_y - T0_y * ps0_x + F0 / BigR * T0_p)                           &
                                    * ( v_x * ps0_y -  v_y * ps0_x                   ) * xjac * theta * tstep * tstep 
    
                      amat_k(var_T,var_vpar) =                                                                        &
                          + tgnum_T * 0.25d0 / BigR * 2.d0 * vpar0*vpar                                               &
                                    * T0 * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                           &
                                    * (                            + F0 / BigR * v_p)  * xjac * theta * tstep * tstep &
                          + tgnum_T * 0.25d0 / BigR * 2.d0 * vpar0*vpar                                               &
                                    * r0 * (T0_x * ps0_y - T0_y * ps0_x + F0 / BigR * T0_p)                           &
                                    * (                            + F0 / BigR * v_p)  * xjac * theta * tstep * tstep 
    
                      amat_n(var_T,var_vpar) = + v * r0 * GAMMA * T0 * F0 / BigR * vpar_p          * xjac * theta * tstep

                      if (with_neutrals) then
                        amat(var_T,var_rhon) = + v * BigR * r0 * rhon * ksiion * Sion_T              * xjac * theta * tstep &
                                               + v * BigR * rhon * r0_corr * LradDrays_T             * xjac * theta * tstep &
                              !===================== Additional terms from friction terms============
                              - v * BigR * ((GAMMA - 1.)/2.) * vpar0**2 * BB2 * (r0*rhon*Sion_T)     * xjac * theta * tstep &
                              - v * BigR * ((GAMMA - 1.)/2.) * vv2            * (r0*rhon*Sion_T)     * xjac * theta * tstep 
                              !==============================End of friction terms=================

                      endif

                    end if ! (with_vpar)
                    
                  end if ! (with_TiTe) *************************************************************
 
                  !###################################################################################################
                  !#  Neutral density equation                                                                       #
                  !###################################################################################################

                  if (with_neutrals) then              

                    amat(var_rhon,var_psi) = + delta_n_convection*(                                                                                   &
                                             + v * rn0   * (vpar0_s * psi_t - vpar0_t * psi_s)         * theta * tstep &
                                             + v * Vpar0 * (rn0_s   * psi_t - rn0_t   * psi_s)         * theta * tstep )

                    amat(var_rhon,var_u) = delta_n_convection*(                                                        &
                                  + v * BigR**2 * ( rn0_s * u_t - rn0_t * u_s)                         * theta * tstep &
                                  + v * 2.d0 * BigR * rn0 * u_y                                 * xjac * theta * tstep )

                    amat(var_rhon,var_rho) = + BigR * v * rn0 * Sion_T * rho                    * xjac * theta * tstep &
                                - BigR * v * 2d0 * r0 * rho * Srec_T                            * xjac * theta * tstep 

                 ! We do not include the term coming from div(rhon * v_star_i) because they are prop. to rho_n/rho, because they may cause problems
                 ! in areas where rho is small.   
              
                    if (with_TiTe) then
                      amat(var_rhon,var_Te) = + BigR * v * r0 * rn0 * dSion_dT * Te   * xjac * theta * tstep &
                                              - BigR * v * r0 * r0  * dSrec_dT * Te   * xjac * theta * tstep       
                    else
                      amat(var_rhon,var_T) = + BigR * v * r0 * rn0 * dSion_dT * T     * xjac * theta * tstep &
                                             - BigR * v * r0 * r0  * dSrec_dT * T     * xjac * theta * tstep       
                    endif  ! with_TiTe

                    if (with_vpar) then 

                      amat(var_rhon,var_vpar) = + delta_n_convection * ( v * F0 / BigR * Vpar * rn0_p     * xjac * theta * tstep &
                                                + v * Vpar * (rn0_s * ps0_t - rn0_t * ps0_s)                     * theta * tstep &
                                                + v * rn0 * (vpar_s * ps0_t - vpar_t * ps0_s)                    * theta * tstep )

                      amat_n(var_rhon,var_vpar) = + delta_n_convection * v * rn0 * F0 / BigR * vpar_p     * xjac * theta * tstep  

                    endif ! with_vpar

                    amat(var_rhon,var_rhon) = + v * rhon * BigR * xjac * (1.d0 + zeta)   &

                              + delta_n_convection*(                                                                                     &
                                - v * BigR**2 * ( rhon_s * u0_t - rhon_t * u0_s)                       * theta * tstep &
                                - v * 2.d0 * BigR * rhon * u0_y                                 * xjac * theta * tstep &
                                + v * rhon * (vpar0_s * ps0_t - vpar0_t * ps0_s)                       * theta * tstep &
                                + v * Vpar0 * (rhon_s * ps0_t - rhon_t * ps0_s)                        * theta * tstep &
                                + v * F0 / BigR * rhon * vpar0_p                                * xjac * theta * tstep ) &
                                   
                     + BigR * (Dn0x * rhon_x * v_x + Dn0y * rhon_y * v_y)                       * xjac * theta * tstep &   
                     + BigR * v * r0 * rhon* Sion_T                                             * xjac * theta * tstep &
                     + Dn_perp_num * (v_xx + v_x/BigR + v_yy)*(rhon_xx + rhon_x/BigR + rhon_yy) * BigR * xjac * theta * tstep 

                     amat_n(var_rhon,var_rhon)  = + delta_n_convection*(                                                             &
                                                  + v * F0 / BigR * Vpar0 * rhon_p                                * xjac * theta * tstep )

                     amat_kn(var_rhon,var_rhon) = + BigR * ( + Dn0p * rhon_p * v_p/BigR**2)                   * xjac * theta * tstep    
                  endif ! with_neutrals 
 
                  !###################################################################################################
                  !# end equations                                                                                   #
                  !###################################################################################################

                  if (use_fft) then
                    index_kl = n_var*(n_order+1)*(k-1) + n_var*(l-1) + 1
                  else
                    index_kl = n_tor_local*n_var*(n_order+1)*(k-1) + n_tor_local*n_var*(l-1) + in - n_tor_start +1
                  endif

                  ! --- Fill up the matrix
                  if (use_fft) then
 
                    do kl = 1, n_var
                      do ij = 1, n_var

                        ELM_p(mp,index_kl+kl-1,ij)  =  ELM_p(mp,index_kl+kl-1,ij) + wst * amat(ij,kl)
                        ELM_n(mp,index_kl+kl-1,ij)  =  ELM_n(mp,index_kl+kl-1,ij) + wst * amat_n(ij,kl)
                        ELM_k(mp,index_kl+kl-1,ij)  =  ELM_k(mp,index_kl+kl-1,ij) + wst * amat_k(ij,kl)
                        ELM_kn(mp,index_kl+kl-1,ij) =  ELM_kn(mp,index_kl+kl-1,ij) + wst * amat_kn(ij,kl)
                      
                      enddo
                    enddo

                  else

                    do kl = 1, n_var
                      do ij = 1, n_var

                        ELM(index_ij+(ij-1)*(n_tor_local),index_kl+(kl-1)*(n_tor_local)) = &
                        ELM(index_ij+(ij-1)*(n_tor_local),index_kl+(kl-1)*(n_tor_local))   &
                          + (amat(ij,kl) + amat_k(ij,kl) + amat_n(ij,kl) + amat_kn(ij,kl)) * wst
                      
                      enddo
                    enddo

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
!--- THIS IS ONLY USED FOR DIAGNOSTIC PURPOSES-------------------------
!--- ELM structure is re-used to plot separate terms in vtk ----------- 
if (present(get_terms)) then

  do i_term=1, max_terms

    do j=1, n_vertex_max*n_var*(n_order+1)
    
      in_fft = ELM_p(1:n_plane,i_term, j)
#ifdef USE_FFTW
      call dfftw_execute_dft_r2c(fftw_plan, in_fft, out_fft)
#else
      call my_fft(in_fft, out_fft, n_plane)
#endif
        
      index = n_tor*(j-1) + 1
      ELM(i_term,index) = real(out_fft(1))
    
      do k=2,(n_tor+1)/2
        index = n_tor*(j-1) + 2*(k-1)
        ELM(i_term,index)   =   real(out_fft(k))
        ELM(i_term,index+1) = - imag(out_fft(k))
      enddo
    
    enddo
    
    do j=1, n_vertex_max*n_var*(n_order+1)
    
      in_fft = ELM_k(1:n_plane,i_term,j)
#ifdef USE_FFTW
      call dfftw_execute_dft_r2c(fftw_plan, in_fft, out_fft)
#else
      call my_fft(in_fft, out_fft, n_plane)
#endif
      
      index = n_tor*(j-1) + 1
      ik    = 1
      ELM(i_term,index) = ELM(i_term,index) + imag(out_fft(1)) * float(mode(ik))
    
      do k=2,(n_tor+1)/2
        ik    = max(2*(k-1),1)
        index = n_tor*(j-1) + 2*(k-1)
        ELM(i_term,index)   = ELM(i_term,index)   + imag(out_fft(k)) * float(mode(ik))
        ELM(i_term,index+1) = ELM(i_term,index+1) + real(out_fft(k)) * float(mode(ik))
      enddo
    
    enddo

  enddo ! maxterms
!----------------------------------------------------------------------

else

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

endif !--- get terms

return


CONTAINS

! subroutine that calculates shock-capturing stabilization related terms
subroutine calculate_sc_quantities()
! Define total pressure as (assuming neutrals are at same Temperature as ions)
! P_tot = r0 * (Ti0 + Te0) + rn0 * Ti0         and   r0 * T0 + 0.5d0 * rn0 * T0
if ( with_TiTe ) then 
  Ptot     = Pi0   + Pe0    + rn0 * Ti0
  Ptot_corr= r0_corr * (Ti0_corr + Te0_corr) + rn0_corr * Ti0_corr
  Ptot_p   = Pi0_p + Pe0_p  + rn0 * Ti0_p + rn0_p * Ti0
  Ptot_x   = Pi0_x + Pe0_x  + rn0 * Ti0_x + rn0_x * Ti0
  Ptot_y   = Pi0_y + Pe0_y  + rn0 * Ti0_y + rn0_y * Ti0
else
  Ptot     = P0    + 0.5d0 * rn0 * T0
  Ptot_corr=  r0_corr * T0_corr + 0.5d0 * rn0_corr * T0_corr
  Ptot_p   = P0_p  + 0.5d0 * (rn0 * T0_p + rn0_p * T0)
  Ptot_x   = P0_x  + 0.5d0 * (rn0 * T0_x + rn0_x * T0)
  Ptot_y   = P0_y  + 0.5d0 * (rn0 * T0_y + rn0_y * T0)
endif

d_p = 0.d0

! approximate residual in the density equation: \nabla \cdot (\rho \boldsymbol{v})
R_rho = + BigR**2 * ( r0_x * u0_y - r0_y * u0_x) &
        + 2.d0 * BigR * r0 * u0_y                &
        - F0 / BigR * Vpar0 * r0_p               &
        - Vpar0 * (r0_x * ps0_y - r0_y * ps0_x)  &
        - F0 / BigR * r0 * vpar0_p               &
        - r0 * (vpar0_x * ps0_y - vpar0_y * ps0_x)

! approximate residual in the neutrals density equation: \nabla \cdot (\rho_n  \boldsymbol{v})
R_rhon =  delta_n_convection*( &
       +  BigR**2 * ( rn0_x * u0_y - rn0_y * u0_x)  &
       +  2.d0 * BigR * rn0 * u0_y                  &
       -  rn0 * (vpar0_x * ps0_y - vpar0_y * ps0_x) &
       -  Vpar0 * (rn0_x * ps0_y - rn0_y * ps0_x)   &
       -  F0 / BigR * Vpar0 * rn0_p                 &
       -  F0 / BigR * rn0 * vpar0_p )

! approximate residual in the pressure equations: \boldsymbol{v} \cdot \nabla p + \gamma p \nabla \boldsymbol{v}
if ( with_TiTe ) then ! (with_TiTe)

  R_pi = + r0 * BigR**2 * ( Ti0_x * u0_y - Ti0_y * u0_x)                       &
         + Ti0 * BigR**2 * ( r0_x * u0_y - r0_y * u0_x)                        &
         
         + r0 * Ti0 * 2.d0* GAMMA * BigR * u0_y                                &
         
         - r0 * F0 / BigR * Vpar0 * Ti0_p                                      &
         - Ti0 * F0 / BigR * Vpar0 * r0_p                                      &
         
         - r0 * Vpar0 * (Ti0_x * ps0_y - Ti0_y * ps0_x)                        &
         - Ti0 * Vpar0 * (r0_x * ps0_y - r0_y * ps0_y)                         &
         
         - r0 * Ti0 * GAMMA * (vpar0_x * ps0_y - vpar0_y * ps0_x)              &
         - r0 * Ti0 * GAMMA * F0 / BigR * vpar0_p

  R_pe = + r0 * BigR**2  * (Te0_x * u0_y - Te0_y * u0_x)                       &
         + Te0 * BigR**2 * ( r0_x * u0_y -  r0_y * u0_x)                       &
                                                                        
         + r0 * Te0 * 2.d0* GAMMA * BigR * u0_y                                &
                                                                        
         - r0 * F0 / BigR * Vpar0 * Te0_p                                      &
         - Te0 * F0 / BigR * Vpar0 * r0_p                                      &
                                                                        
         - r0 * Vpar0 * (Te0_x * ps0_y - Te0_y * ps0_x)                        &
         - Te0 * Vpar0 * (r0_x * ps0_y - r0_y * ps0_x)                         &
                                                                        
         - r0 * Te0 * GAMMA * (vpar0_x * ps0_y - vpar0_y * ps0_x)              &
         - r0 * Te0 * GAMMA * F0 / BigR * vpar0_p 
                                                             
  ! 1/BigR removes the factor R from the integrand in (R dR)                 
  d_p = (Ti0 + Te0) * R_rho / BigR +  (R_pi + R_pe) / BigR  + Ti0 * R_rhon / BigR

else

  R_p = + r0 * BigR**2 * (T0_x  * u0_y - T0_y * u0_x)                          &
        + T0 * BigR**2 * ( r0_x * u0_y - r0_y * u0_x)                          &
                                                                   
        + r0 * T0 * 2.d0* GAMMA * BigR * u0_y                                  &
        
        - r0 * F0 / BigR * Vpar0 * T0_p                                        &
        - T0 * F0 / BigR * Vpar0 * r0_p                                        &
        
        - r0 * Vpar0 * (T0_x * ps0_y - T0_y * ps0_x)                           &
        - T0 * Vpar0 * (r0_x * ps0_y - r0_y * ps0_x)                           &
                                                                    
        - r0 * T0 * GAMMA * (vpar0_x * ps0_y - vpar0_y * ps0_x)                &
        - r0 * T0 * GAMMA * F0 / BigR * vpar0_p
                
  ! 1/BigR removes the factor R from the integrand in (R dR)
  d_p = T0 * R_rho / BigR + R_p / BigR + T0 * R_rhon / BigR

endif

! Shock-detector term based on the total pressure gradient
f_p = dsqrt( Ptot_x*Ptot_x + Ptot_y*Ptot_y + Ptot_p*Ptot_p/ (BigR*BigR) ) / Ptot_corr * h_e
! Estimation of the numerical stabilization coefficient
tau_sc = h_e * h_e * abs(d_p) / Ptot_corr * f_p

! Use of source terms to increase the stabilization coefficients
if(add_sources_in_sc)then
  src_rho = (particle_source(ms,mt) + source_pellet) &
          + r0_corr * rn0_corr * Sion_T              &
          - r0_corr * r0_corr  * Srec_T

  src_rhon = - r0_corr * rn0_corr * Sion_T      &
           + r0_corr * r0_corr  * Srec_T        &
           + source_neutral_drift

  if ( with_TiTe ) then ! (with_TiTe)
    src_pi  = heat_source_i(ms,mt)                                      &
            + ((GAMMA - 1.)/2.) * vpar0**2 * BB2 * (r0_corr*rn0*Sion_T) &
            + ((GAMMA - 1.)/2.) * vv2 * ((r0_corr*rn0*Sion_T))

    src_pe  = heat_source_e(ms,mt)                           &
             + (gamma-1.d0) * eta_T_ohm * (zj0 / BigR)**2.d0 &
             - ksiion  * r0_corr * rn0_corr * Sion_T         &
             - r0_corr * rn0_corr * LradDrays_T              &
             - r0_corr * r0_corr  * LradDcont_T              &
             - r0_corr * frad_bg
            
    s_p = (Ti0 + Te_0) * src_rho + (src_pi + src_pe) + Ti0 * src_rhon
  else
    src_p   =  heat_source(ms,mt)                                         &
            +  (gamma-1.d0) * eta_T_ohm * (zj0 / BigR)**2.d0              &
            +  ((GAMMA - 1.)/2.) * vpar0**2 * BB2 * (r0_corr*rn0*Sion_T)  &
            +  ((GAMMA - 1.)/2.) * vv2 * ((r0_corr*rn0*Sion_T))           &
            -  ksiion  * r0_corr * rn0_corr * Sion_T                      &
            -  r0_corr * rn0_corr * LradDrays_T                           &
            -  r0_corr * r0_corr  * LradDcont_T                           &
            -  r0_corr * frad_bg

    s_p = T0 * src_rho + src_p + T0 * src_rhon
  endif

  tau_sc = h_e * h_e * (abs(s_p) + abs(d_p)) / Ptot_corr * f_p

endif

! Updates in the physical diffsivities to locally add numerical stabilization.
visco_T = visco_T + visco_sc_num  * tau_sc
D_prof  = D_prof  + D_perp_sc_num * tau_sc
if ( with_TiTe ) then
  ZKi_prof  = ZKi_prof  + ZK_i_perp_sc_num * tau_sc
  ZKi_par_T = ZKi_par_T + ZK_i_par_sc_num  * tau_sc
  ZKe_prof  = ZKe_prof  + ZK_e_perp_sc_num * tau_sc
  ZKe_par_T = ZKe_par_T + ZK_e_par_sc_num  * tau_sc
else
  ZK_prof  = ZK_prof   + ZK_perp_sc_num * tau_sc
  ZK_par_T = ZK_par_T  + ZK_par_sc_num  * tau_sc
endif
Dn0x = Dn0x + Dn_pol_sc_num * tau_sc
Dn0y = Dn0y + Dn_pol_sc_num * tau_sc
Dn0p = Dn0p + Dn_p_sc_num   * tau_sc

end subroutine calculate_sc_quantities


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
