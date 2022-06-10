module mod_elt_matrix_fft

  implicit none

contains

#include "corr_neg_include.f90"

subroutine element_matrix_fft(element, nodes, xpoint2, xcase2, R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint, Z_xpoint, &
                              ELM, RHS, tid, ELM_p, ELM_n, ELM_k, ELM_kn, RHS_p, RHS_k,                               &
                              eq_g, eq_s, eq_t, eq_p, eq_ss, eq_st, eq_tt, delta_g, delta_s, delta_t,                 &
                              i_tor_min, i_tor_max, aux_nodes, ELM_pnn, get_terms)
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
use corr_neg
use mod_neutral_source
use mod_bootstrap_functions
use mod_atomic_coeff_deuterium, only : atomic_coeff_deuterium
use mod_impurity, only: radiation_function, radiation_function_linear
use mod_sources
use mod_model_settings


implicit none

type (type_element)       :: element 
type (type_node)          :: nodes(n_vertex_max)     ! fluid variables
type (type_node),optional :: aux_nodes(n_vertex_max) ! particle moments

#define DIM0 n_tor*n_vertex_max*n_degrees*n_var

integer, intent(in)            :: tid
integer, intent(in)            :: i_tor_min, i_tor_max
logical, intent(in), optional  :: get_terms   

real*8, dimension (DIM0,DIM0)  :: ELM
real*8, dimension (DIM0)       :: RHS

integer    :: i, j, ms, mt, mp, k, l, index_ij, index_kl, index, index_k, index_m, m, ik, xcase2, i_inj, n_spi_tmp
integer    :: n_tor_start, n_tor_end, n_tor_local
integer    :: in, im, ij1, ij2, ij3, ij4, ij5, ij6, ij7, ij8, kl1, kl2, kl3, kl4, kl5, kl6, kl7, kl8, ij, kl
real*8     :: wst, xjac, xjac_s, xjac_t, xjac_x, xjac_y, BigR, r2, phi, delta_phi
real*8     :: current_source(n_gauss,n_gauss),particle_source(n_gauss,n_gauss),heat_source(n_gauss,n_gauss)
real*8     :: R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint(2), Z_xpoint(2), dj_dpsi, dj_dz, source_pellet, source_volume
real*8     :: Bgrad_rho_star,     Bgrad_rho,     Bgrad_T_star,  Bgrad_T, BB2
real*8     :: Bgrad_rho_star_psi, Bgrad_rho_psi, Bgrad_rho_rho, Bgrad_T_star_psi, Bgrad_T_psi, Bgrad_T_T, BB2_psi
real*8     :: Bgrad_rho_rho_n, Bgrad_T_T_n, Bgrad_rho_k_star, Bgrad_T_k_star, ZKpar_T, dZKpar_dT
real*8     :: D_prof, ZK_prof, psi_norm, theta, zeta, delta_u_x, delta_u_y, delta_ps_x, delta_ps_y
real*8     :: rhs_ij(n_var), rhs_ij_k(n_var)
real*8     :: amat(n_var,n_var), amat_k(n_var,n_var), amat_n(n_var,n_var), amat_kn(n_var,n_var)

real*8     :: v, v_x, v_y, v_s, v_t, v_p, v_ss, v_st, v_tt, v_xx, v_xy, v_yy
real*8     :: ps0, ps0_x, ps0_y, ps0_p, ps0_s, ps0_t, ps0_ss, ps0_tt, ps0_st, ps0_xx, ps0_yy, ps0_xy
real*8     :: zj0, zj0_x, zj0_y, zj0_p, zj0_s, zj0_t
real*8     :: u0, u0_x, u0_y, u0_p, u0_s, u0_t, u0_ss, u0_tt, u0_st, u0_xx, u0_xy, u0_yy
real*8     :: w0, w0_x, w0_y, w0_p, w0_s, w0_t, w0_ss, w0_st, w0_tt, w0_xx, w0_xy, w0_yy
real*8     :: r0, r0_x, r0_y, r0_p, r0_s, r0_t, r0_ss, r0_st, r0_tt, r0_xx, r0_xy, r0_yy, r0_hat, r0_x_hat, r0_y_hat, r0_corr
real*8     :: T0, T0_x, T0_y, T0_p, T0_s, T0_t, T0_ss, T0_st, T0_tt, T0_xx, T0_xy, T0_yy, T0_corr, dT0_corr_dT
real*8     :: psi, psi_x, psi_y, psi_p, psi_s, psi_t, psi_ss, psi_st, psi_tt, psi_xx, psi_xy, psi_yy
real*8     :: zj, zj_x, zj_y, zj_p, zj_s, zj_t, zj_ss, zj_st, zj_tt
real*8     :: u, u_x, u_y, u_p, u_s, u_t, u_ss, u_st, u_tt, u_xx, u_xy, u_yy
real*8     :: w, w_x, w_y, w_p, w_s, w_t, w_ss, w_st, w_tt, w_xx, w_xy, w_yy
real*8     :: rho, rho_x, rho_y, rho_s, rho_t, rho_p, rho_hat, rho_x_hat, rho_y_hat, rho_ss, rho_st, rho_tt, rho_xx, rho_xy, rho_yy
real*8     :: T, T_x, T_y, T_s, T_t, T_p, T_ss, T_st, T_tt, T_xx, T_xy, T_yy
real*8     :: Ti0, Ti0_x, Ti0_y, Te0, Te0_x, Te0_y
real*8     :: zTi, zTi_x, zTi_y, zTe, zTe_x, zTe_y, zn_x, zn_y
real*8     :: rn0, rn0_s, rn0_t, rn0_x, rn0_y, rn0_p, rn0_ss, rn0_st, rn0_tt, rn0_xx, rn0_yy
real*8     :: rn0_hat, rn0_x_hat, rn0_y_hat, rn0_corr
real*8     :: rhon, rhon_s, rhon_t, rhon_p, rhon_x, rhon_y, rhon_ss, rhon_tt, rhon_st, rhon_xx, rhon_yy, rhon_xy
real*8	   :: Jb_0 , Jb, dn0x, dn0y, dn0p
real*8     :: Vpar, Vpar_x, Vpar_y, Vpar_p, Vpar_s, Vpar_t, Vpar_ss, Vpar_st, Vpar_tt, Vpar_xx, Vpar_yy, Vpar_xy
real*8     :: P0, P0_s, P0_t, P0_x, P0_y, P0_p, P0_ss, P0_st, P0_tt, P0_xx, P0_xy, P0_yy
real*8     :: P0_x_rho, P0_xx_rho, P0_y_rho, P0_yy_rho, P0_xy_rho
real*8     :: P0_x_T,   P0_xx_T,   P0_y_T,   P0_yy_T,   P0_xy_T
real*8     :: Vpar0, Vpar0_s, Vpar0_t, Vpar0_p, Vpar0_x, Vpar0_y, Vpar0_ss, Vpar0_st, Vpar0_tt, Vpar0_xx, Vpar0_yy,Vpar0_xy
real*8     :: BigR_x, vv2, eta_T, visco_T, deta_dT, d2eta_d2T, dvisco_dT, d2visco_dT2, visco_num_T, eta_num_T, W_dia, W_dia_rho, W_dia_T
real*8     :: eta_T_ohm, deta_dT_ohm
real*8     :: ZK_par_num, T0_ps0_x, T_ps0_x, T0_psi_x, T0_ps0_y, T_ps0_y, T0_psi_y, v_ps0_x, v_psi_x, v_ps0_y, v_psi_y
real*8     :: TG_num1, TG_num2, TG_num5, TG_num6, TG_num7

real*8     :: Vt0,Omega_tor0_x,Omega_tor0_y,Vt0_x,Vt0_y
real*8     :: V_source(n_gauss,n_gauss), Vt_x_psi, Vt_y_psi, Omega_tor_x_psi, Omega_tor_y_psi
real*8     :: dV_dpsi_source(n_gauss,n_gauss),dV_dz_source(n_gauss,n_gauss)
real*8     :: dV_dpsi2,dV_dz2,dV_dpsi_dz,dV_dpsi3,dV_dpsi_dz2,dV_dpsi2_dz
real*8     :: eq_zne(n_gauss,n_gauss), eq_zTe(n_gauss,n_gauss)
real*8     :: dn_dpsi(n_gauss,n_gauss),dn_dz,dn_dpsi2,dn_dz2,dn_dpsi_dz,dn_dpsi3,dn_dpsi_dz2,dn_dpsi2_dz
real*8     :: dT_dpsi(n_gauss,n_gauss),dT_dz,dT_dpsi2,dT_dz2,dT_dpsi_dz,dT_dpsi3,dT_dpsi_dz2,dT_dpsi2_dz
logical    :: xpoint2, use_fft
real*8     :: Btheta2, epsil, Btheta2_psi
real*8, dimension(n_gauss,n_gauss)    :: amu_neo_prof, aki_neo_prof
! neutral source
real*8     :: source_neutral
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
real*8     :: Te_corr_eV, Te_eV                               ! Temperature used in radiation rate
real*8     :: ne_SI                                           ! Electron density used in radiation rate

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

real*8, dimension(n_tor,n_plane) :: HHZ, HHZ_p, HHZ_pp

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
zk_par_num = 0.d0

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
heat_source     = 0.d0
V_source        = 0.d0
dV_dpsi_source  = 0.d0
dV_dz_source    = 0.d0
eq_zne          = 0.d0
eq_zTe          = 0.d0         

amu_neo_prof   = 0.d0
aki_neo_prof   = 0.d0

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

    if ( NEO ) then 
      if (num_neo_file) then
        call neo_coef( xpoint2, xcase2, y_g(ms,mt), Z_xpoint, eq_g(1,1,ms,mt),psi_axis,psi_bnd, amu_neo_prof(ms,mt), aki_neo_prof(ms,mt))
      else
        amu_neo_prof(ms,mt) = amu_neo_const
        aki_neo_prof(ms,mt) = aki_neo_const
      endif
    endif

  enddo
enddo

eq_zTe = eq_zTe / 2.d0  ! electron temperature  

!--------------------------------------------------- sum over the Gaussian integration points
do i=1,n_vertex_max
  do j=1,n_degrees

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

          Vt0   = V_source(ms,mt)
          if (normalized_velocity_profile) then
            Vt0_x = dV_dpsi_source(ms,mt)*ps0_x
            Vt0_y = dV_dz_source(ms,mt)+dV_dpsi_source(ms,mt)*ps0_y
          else
            Omega_tor0_x = dV_dpsi_source(ms,mt)*ps0_x
            Omega_tor0_y = dV_dz_source(ms,mt)+dV_dpsi_source(ms,mt)*ps0_y
          endif

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
          r0_corr = corr_neg_dens(r0)
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

          rn0    = eq_g(mp,8,ms,mt)
          rn0_x  = (   y_t(ms,mt) * eq_s(mp,8,ms,mt) - y_s(ms,mt) * eq_t(mp,8,ms,mt) ) / xjac    
          rn0_y  = ( - x_t(ms,mt) * eq_s(mp,8,ms,mt) + x_s(ms,mt) * eq_t(mp,8,ms,mt) ) / xjac   
          rn0_p  = eq_p(mp,8,ms,mt)                                                             
          rn0_s  = eq_s(mp,8,ms,mt)                                                             
          rn0_t  = eq_t(mp,8,ms,mt)                                                             
          rn0_ss = eq_ss(mp,8,ms,mt)                                                            
          rn0_st = eq_st(mp,8,ms,mt)                                                            
          rn0_tt = eq_tt(mp,8,ms,mt)  

          rn0_corr = corr_neg_dens(rn0, (/ 0.d-5, 1.d-5 /)) ! Correction for negative rn0 ...
     
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

          T0    = eq_g(mp,6,ms,mt)
          T0_x  = (   y_t(ms,mt) * eq_s(mp,6,ms,mt) - y_s(ms,mt) * eq_t(mp,6,ms,mt) ) / xjac
          T0_y  = ( - x_t(ms,mt) * eq_s(mp,6,ms,mt) + x_s(ms,mt) * eq_t(mp,6,ms,mt) ) / xjac
          T0_p  = eq_p(mp,6,ms,mt)
          T0_s  = eq_s(mp,6,ms,mt)
          T0_t  = eq_t(mp,6,ms,mt)
          T0_ss = eq_ss(mp,6,ms,mt)
          T0_tt = eq_tt(mp,6,ms,mt)
          T0_st = eq_st(mp,6,ms,mt)

          T0_corr     = corr_neg_temp(T0) ! For use in eta(T), visco(T), ...
          dT0_corr_dT = dcorr_neg_temp_dT(T0) ! Improve the correction

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
            dZKpar_dT = ZK_par * (2.5d0)  * T0_corr**(+1.5d0) * T_0**(-2.5d0) * dT0_corr_dT
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

          eta_num_T   = eta_num                         ! hyperresistivity
          visco_num_T = visco_num                       ! hyperviscosity

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

          Dn0x = D_neutral_x      
          Dn0y = D_neutral_y      
          Dn0p = D_neutral_p    
  
          if (use_pellet) then

            call pellet_source2(pellet_amplitude,pellet_R,pellet_Z,pellet_psi,pellet_phi, &
                                pellet_radius, pellet_delta_psi, pellet_sig, pellet_length, pellet_ellipse, pellet_theta, &
                                x_g(ms,mt),y_g(ms,mt), ps0, phi, r0_corr, T0_corr/2.d0, &
                                central_density, pellet_particles, pellet_density, total_pellet_volume, &
                                source_pellet, source_volume)
          endif

          call atomic_coeff_deuterium(0.5d0*T0, Sion_T, dSion_dT, Srec_T, dSrec_dT,        &
                                      LradDcont_T, dLradDcont_dT, LradDrays_T, dLradDrays_dT, r0 )

          ! --- Transform derivatives on Te to derivatives in total T
          dSion_dT      = dSion_dT      / 2.d0
          dSrec_dT      = dSrec_dT      / 2.d0
          dLradDrays_dT = dLradDrays_dT / 2.d0
          dLradDcont_dT = dLradDcont_dT / 2.d0

          !-------------------------------------------
          ! --- Normalisation of the ionization energy cost for Deuterium
          !------------------------------------------- 
      
          ksiion = central_density * 1.d20 * ksi_ion
      
          !--------------------------------------------------------
          ! --- Source of neutrals, e.g. from MGI/SPI
          !--------------------------------------------------------

          source_neutral = 0.d0                   
     
          call total_neutral_source(x_g(ms,mt),y_g(ms,mt),phi,ps0,source_neutral)
     
          source_neutral = max(source_neutral,0.)
      
         !-----------------------------------------------------------------
         ! --- Radiation from background impurity, using ADAS (by default)
         !-----------------------------------------------------------------

          ne_SI = r0_corr * 1.d20 * central_density !electron density (SI)
          Te_corr_eV = T0_corr/(2.d0*EL_CHG*MU_ZERO*central_density*1.d20)  ! Te in eV

          Te_eV = T0/(2.d0*EL_CHG*MU_ZERO*central_density*1.d20)  ! Te in eV, uncorrected

          if (use_imp_adas) then  ! use open adas by default
            frad_bg = 0. 
            dfrad_bg_dT = 0.
            do i_imp =1, n_adas
              r_imp = nimp_bg(i_imp) / (1.d20 * central_density)  ! Background impurity density in JU     
              if (ne_SI > ne_SI_min .and. Te_eV > Te_eV_min .and. r_imp > 0) then
                Lrad_imp = 0.0
                dLrad_imp_dT = 0.0
                call radiation_function_linear(imp_adas(i_imp),imp_cor(i_imp),log10(ne_SI),   & 
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

            Bgrad_T          = ( F0 / BigR * T0_p +  T0_x * ps0_y - T0_y * ps0_x ) / BigR

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
            !#  equation 1   (induction equation)                                                              #
            !###################################################################################################

            rhs_ij(1) = v * eta_T  * (zj0 - current_source(ms,mt) - Jb)/ BigR * xjac * tstep * factor(1,1) & ! Resistive
                      + v * (ps0_s * u0_t - ps0_t * u0_s)                            * tstep * factor(1,2) & ! Poisson bracket
                      - v * F0 / BigR  * u0_p                                 * xjac * tstep * factor(1,2) &
                      + eta_num_T * (v_x * zj0_x + v_y * zj0_y)               * xjac * tstep * factor(1,3) & ! Hyper-resistivity

                      - v * tauIC/(r0_corr*BB2) * F0**2/BigR**2 * (ps0_s * p0_t - ps0_t * p0_s) * tstep * factor(1,4) & ! Diamagnetic
                      + v * tauIC/(r0_corr*BB2) * F0**3/BigR**3 * p0_p        * xjac * tstep * factor(1,4) &

                      + zeta * v * delta_g(mp,1,ms,mt) / BigR                 * xjac* factor(1,5) ! Time variation 


            !###################################################################################################
            !#  equation 2   (perpendicular momentum equation)                                                 #
            !###################################################################################################

            rhs_ij(2) =  - 0.5d0 * vv2 * (v_x * r0_y_hat - v_y * r0_x_hat)   * xjac * tstep * factor(2,1) & ! Convection
                         - r0_hat * BigR**2 * w0 * (v_s * u0_t - v_t * u0_s)        * tstep * factor(2,1) &
                         + v * (ps0_s * zj0_t - ps0_t * zj0_s )                     * tstep * factor(2,2) & ! JXB
                         - visco_T * BigR * (v_x * w0_x + v_y * w0_y)        * xjac * tstep * factor(2,3) & ! Newtonian visco.
                         - v * F0 / BigR * zj0_p                             * xjac * tstep * factor(2,2) &
                         + BigR**2 * (v_s * p0_t - v_t * p0_s)                      * tstep * factor(2,4) & ! grad P

                         - visco_num_T * (v_xx + v_x/Bigr + v_yy)*(w0_xx + w0_x/Bigr + w0_yy) * xjac * tstep * factor(2,5) & ! Hyper-viscosity 

                         - TG_num2 * 0.25d0 * r0_hat * BigR**3 * (w0_x * u0_y - w0_y * u0_x)  &
                                   * ( v_x * u0_y - v_y * u0_x) * xjac * tstep * tstep * factor(2,6) & ! TG stabilization

                         - v * tauIC * BigR**4 * (p0_s * w0_t - p0_t * w0_s)        * tstep * factor(2,7) & ! Diamagnetic

                         - tauIC * BigR**3 * p0_y * (v_x* u0_x + v_y * u0_y) * xjac * tstep * factor(2,7) &

                         - v * tauIC * BigR**4 * (u0_xy * (p0_xx - p0_yy) - p0_xy * (u0_xx - u0_yy) ) * xjac * tstep * factor(2,7) &

                         ! --- Diamagnetic viscosity
                         + dvisco_dT * bigR * W_dia * (v_x*T0_x + v_y*T0_y)    * xjac * tstep  * factor(2,8) &
                         + visco_T   * bigR * W_dia * (v_xx + v_x/bigR + v_yy) * xjac * tstep  * factor(2,8) &

                         + BigR**3 * (particle_source(ms,mt) + source_pellet) * (v_x * u0_x + v_y * u0_y) * xjac* tstep * factor(2,9) & ! External density source
                     
                         + (1.d0 - delta_n_convection) * (   &
                               + BigR**3*(r0_corr*rn0_corr*Sion_T)*(v_x * u0_x + v_y * u0_y)  * xjac * tstep * factor(2,10) & ! Density source from ionization
                               - BigR**3*(r0_corr*r0_corr *Srec_T)*(v_x * u0_x + v_y * u0_y)  * xjac * tstep * factor(2,11) & ! Density sink by recombination
                           )  &

                         - zeta * BigR * r0_hat * (v_x * delta_u_x + v_y * delta_u_y) * xjac* factor(2,12) ! Variation
            
            if (NEO) then
              rhs_ij(2) = rhs_ij(2)  + amu_neo_prof(ms,mt)*BB2/((Btheta2+epsil)**2)   * factor(2,13)        & ! Neoclassical
                        * (ps0_x * v_x + ps0_y * v_y) *                                      &  
                          ( r0 * (ps0_x * u0_x + ps0_y * u0_y)                               &
                            + tauIC * (ps0_x * P0_x + ps0_y * P0_y)                          &
                            + aki_neo_prof(ms,mt) * tauIC * r0 * (ps0_x*T0_x + ps0_y*T0_y)   &
                            - r0 * Vpar0 * Btheta2 )     * xjac * tstep * BigR 
            endif

            !###################################################################################################
            !#  equation 3   (current definition)                                                              #
            !###################################################################################################

            rhs_ij(3) = - ( v_x * ps0_x  + v_y * ps0_y + v*zj0) / BigR * xjac* factor(3,1) 

            !###################################################################################################
            !#  equation 4   (vorticity definition)                                                            #
            !###################################################################################################

            rhs_ij(4) = - ( v_x * u0_x   + v_y * u0_y  + v*w0)  * BigR * xjac * factor(4,1) 

            !###################################################################################################
            !#  equation 5   (density equation)                                                                #
            !###################################################################################################

            rhs_ij(5)  = v * BigR * (particle_source(ms,mt) + source_pellet)                      * xjac * tstep * factor(5,1) & ! External density source
                       + v * BigR**2 * ( r0_s * u0_t - r0_t * u0_s)                                      * tstep * factor(5,2) & ! Convection perp.
                       + v * 2.d0 * BigR * r0 * u0_y                                              * xjac * tstep * factor(5,3) & 
                       - (D_par-D_prof) * BigR / BB2 * Bgrad_rho_star * Bgrad_rho                 * xjac * tstep * factor(5,4) & ! Diffusion paral.
                       - D_prof * BigR  * (v_x*r0_x + v_y*r0_y                                  ) * xjac * tstep * factor(5,5) & ! Diffusion perp.
                       - v * F0 / BigR * Vpar0 * r0_p                                             * xjac * tstep * factor(5,6) & ! Convection paral.
                       - v * Vpar0 * (r0_s * ps0_t - r0_t * ps0_s)                                       * tstep * factor(5,6) &
                       - v * F0 / BigR * r0 * vpar0_p                                             * xjac * tstep * factor(5,3) &
                       - v * r0 * (vpar0_s * ps0_t - vpar0_t * ps0_s)                                    * tstep * factor(5,3) &

                       + v * 2.d0 * tauIC * p0_y * BigR                                           * xjac * tstep * factor(5,7) & ! Diamagnetic

                       + v * r0_corr * rn0_corr * BigR * Sion_T                                   * xjac * tstep * factor(5,8) & ! Source from ionization
                       - v * r0_corr * r0_corr  * BigR * Srec_T                                   * xjac * tstep * factor(5,9) & ! Sink by recombination
                       
                       + zeta * v * delta_g(mp,5,ms,mt) * BigR                                    * xjac         * factor(5,10)   & ! Variation

                       - D_perp_num * (v_xx + v_x/Bigr + v_yy)*(r0_xx + r0_x/Bigr + r0_yy) * BigR * xjac * tstep * factor(5,11) & !Hyper particle diffusitivity

                       - TG_num5 * 0.25d0 * BigR**3 * (r0_x * u0_y - r0_y * u0_x)                                               & 
                                                    * ( v_x * u0_y - v_y * u0_x) * xjac * tstep * tstep          * factor(5,12) & ! TG stabilization
                       - TG_num5 * 0.25d0 / BigR * vpar0**2                                                                     &
                                 * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                              * factor(5,12) &
                                 * ( v_x * ps0_y -  v_y * ps0_x                   ) * xjac * tstep * tstep 


            rhs_ij_k(5) = - (D_par-D_prof) * BigR / BB2 * Bgrad_rho_k_star * Bgrad_rho            * xjac * tstep * factor(5,4) &
                          - D_prof * BigR  * (                  v_p*r0_p /BigR**2 )               * xjac * tstep * factor(5,5) &

                       - TG_num5 * 0.25d0 / BigR * vpar0**2 * factor(5,12) &
                                 * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                                            &
                                 * (                            + F0 / BigR * v_p) * xjac * tstep * tstep 

            !###################################################################################################
            !#  equation 6   (energy equation)                                                                 #
            !###################################################################################################

            rhs_ij(6) =  v * BigR * heat_source(ms,mt)                                    * xjac * tstep * factor(6, 1)& ! External heat source
            
                       + v * r0 * BigR**2 * ( T0_s * u0_t - T0_t * u0_s)                         * tstep * factor(6, 2)& ! vperp*gradP
                       + v * T0 * BigR**2 * ( r0_s * u0_t - r0_t * u0_s)                         * tstep * factor(6, 2)&

                       + v * r0 * T0 * 2.d0* GAMMA * BigR * u0_y                          * xjac * tstep * factor(6, 3)& !gamma*P*divergence(vperp)

                       - v * r0 * F0 / BigR * Vpar0 * T0_p                                * xjac * tstep * factor(6, 4)& !vpar*gradP
                       - v * T0 * F0 / BigR * Vpar0 * r0_p                                * xjac * tstep * factor(6, 4)&

                       - v * r0 * Vpar0 * (T0_s * ps0_t - T0_t * ps0_s)                          * tstep * factor(6, 4)&
                       - v * T0 * Vpar0 * (r0_s * ps0_t - r0_t * ps0_s)                          * tstep * factor(6, 4)&

                       - v * r0 * T0 * GAMMA * (vpar0_s * ps0_t - vpar0_t * ps0_s)               * tstep * factor(6, 5)& !gamma*P*divergence(vpar)
                       - v * r0 * T0 * GAMMA * F0 / BigR * vpar0_p                        * xjac * tstep * factor(6, 5)&

                       - (ZKpar_T-ZK_prof) * BigR / BB2 * Bgrad_T_star * Bgrad_T     * xjac      * tstep * factor(6, 6)& ! Diffusion paral.
                       - ZK_prof * BigR * (v_x*T0_x + v_y*T0_y                     ) * xjac      * tstep * factor(6, 7)& ! Diffusion perp.
 
                       - ZK_perp_num  *  (v_xx + v_x/Bigr + v_yy)*(T0_xx + T0_x/Bigr + T0_yy) * BigR * xjac * tstep * factor(6,8)& ! Hyper perp. diffusitivity

                       - ZK_par_num * (v_ps0_x  * ps0_y - v_ps0_y  * ps0_x)                                            &
                                    * (T0_ps0_x * ps0_y - T0_ps0_y * ps0_x)               * xjac * tstep * factor(6, 9)&  ! Hyper paral. diffusitivity

                       - v * BigR * ksiion * r0_corr * rn0_corr * Sion_T                  * xjac * tstep * factor(6,10)& ! Energy sink by ionization

!===================== Additional terms from friction terms============
                       + v * BigR * ((GAMMA - 1.)/2.) * vpar0**2 * BB2 * (r0_corr*rn0*Sion_T) * xjac * tstep * factor(6,11)& ! Friction
                       + v * BigR * ((GAMMA - 1.)/2.) * vv2 * ((r0_corr*rn0*Sion_T))          * xjac * tstep * factor(6,11)&
!==============================End of friction terms=================

                       + v * (gamma-1.d0) * eta_T_ohm * (zj0 / BigR)**2.d0         * BigR  * xjac * tstep * factor(6,12) & ! Source from Ohmic heating
                       - v * BigR * r0_corr * rn0_corr * LradDrays_T                       * xjac * tstep * factor(6,13) & ! Sink by line radiation
                       - v * BigR * r0_corr * r0_corr  * LradDcont_T                       * xjac * tstep * factor(6,14) & ! Sink by Brem. radiation
                       - v * BigR * r0_corr * frad_bg                                      * xjac * tstep * factor(6,15) & ! Sink by background impurity rad.

                       - TG_num6 * 0.25d0 * BigR**3 * T0 * (r0_x * u0_y - r0_y * u0_x)                                 &
                                          * ( v_x * u0_y - v_y * u0_x) * xjac * tstep * tstep            * factor(6,16)& ! TG stabilization
                       - TG_num6 * 0.25d0 * BigR**3 * r0 * (T0_x * u0_y - T0_y * u0_x)                                 &
                                          * ( v_x * u0_y - v_y * u0_x) * xjac * tstep * tstep            * factor(6,16)&
                       - TG_num6 * 0.25d0 / BigR * vpar0**2                                                            &
                                 * T0 * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                               &
                                 * ( v_x * ps0_y -  v_y * ps0_x             ) * xjac * tstep * tstep     * factor(6,16)&
                       - TG_num6 * 0.25d0 / BigR * vpar0**2                                                            &
                                 * r0 * (T0_x * ps0_y - T0_y * ps0_x + F0 / BigR * T0_p)                               &
                                 * ( v_x * ps0_y -  v_y * ps0_x             ) * xjac * tstep * tstep     * factor(6,16)&

                       + zeta * v * r0 * delta_g(mp,6,ms,mt) * BigR                       * xjac         * factor(6,17)& ! Variation
                       + zeta * v * T0 * delta_g(mp,5,ms,mt) * BigR                       * xjac         * factor(6,17)


            rhs_ij_k(6) = - (ZKpar_T-ZK_prof) * BigR / BB2 * Bgrad_T_k_star * Bgrad_T * xjac * tstep     * factor(6,6)&
                          - ZK_prof * BigR * (                + v_p*T0_p /BigR**2 )   * xjac * tstep     * factor(6,7)&

                         - TG_num6 * 0.25d0 / BigR * vpar0**2 &
                                 * T0 * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                         &
                                 * (                                 + F0 / BigR * v_p) * xjac * tstep * tstep * factor(6,16)  &
                         - TG_num6 * 0.25d0 / BigR * vpar0**2 &
                                 * r0 * (T0_x * ps0_y - T0_y * ps0_x + F0 / BigR * T0_p)                         &
                                 * (                                 + F0 / BigR * v_p) * xjac * tstep * tstep * factor(6,16)

            !###################################################################################################
            !#  equation 7   (parallel velocity equation)                                                      #
            !###################################################################################################

            rhs_ij(7) = - v * F0 / BigR * P0_p                                             * xjac * tstep * factor(7,1) & ! gradP*B
                        - v * (P0_s * ps0_t - P0_t * ps0_s)                                       * tstep * factor(7,1) &

                     - v*(particle_source(ms,mt) + source_pellet) * vpar0 * BB2   * BigR * xjac * tstep * factor(7,2) & ! External particle source

                     - 0.5d0 * r0 * vpar0**2 * BB2 * (ps0_s * v_t - ps0_t * v_s)                * tstep * factor(7,3) & ! Convection paral.
                     - 0.5d0 * v  * vpar0**2 * BB2 * (ps0_s * r0_t - ps0_t * r0_s)              * tstep * factor(7,3) &
                     + 0.5d0 * v  * vpar0**2 * BB2 * F0 / BigR * r0_p                    * xjac * tstep * factor(7,3) &

                     - visco_par_num * (v_xx + v_x/Bigr + v_yy)*(vpar0_xx + vpar0_x/Bigr + vpar0_yy) * BigR * xjac * tstep * factor(7,4) & ! Hyper viscosity

                     + zeta * v * delta_g(mp,7,ms,mt) * R0 * F0**2 / BigR                        * xjac * factor(7,5) & ! Time variation
                     + zeta * v * r0 * vpar0 * (ps0_x * delta_ps_x + ps0_y * delta_ps_y) / BigR  * xjac * factor(7,5) &

             - TG_NUM7 * 0.25d0 * r0 * Vpar0**2 * BB2 &
                       * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac + F0 / BigR * vpar0_p) / BigR   &
                       * (-(ps0_s * v_t     - ps0_t * v_s)    /xjac                      ) * xjac * tstep * tstep * factor(7,6) & ! TG stabilization
             - TG_NUM7 * 0.25d0 * v  * Vpar0**2 * BB2 &
                       * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac + F0 / BigR * vpar0_p) / BigR   &
                       * (-(ps0_s * r0_t    - ps0_t * r0_s)   /xjac + F0 / BigR * r0_p)    * xjac * tstep * tstep * factor(7,6) &
                 
                    + (1.d0 - delta_n_convection) * (     &
                    - v *(r0_corr * rn0_corr * Sion_T) * vpar0 * BB2 * BigR                        * xjac * tstep * factor(7,7) & ! Ionization 
                    + v *(r0_corr * r0_corr  * Srec_T) * vpar0 * BB2 * BigR                        * xjac * tstep * factor(7,8) & ! Recombination
                    )
 
            if (normalized_velocity_profile) then
              rhs_ij(7) = rhs_ij(7) - visco_par * (v_x * (vpar0_x-Vt0_x) + v_y * (vpar0_y-Vt0_y)) * BigR* xjac * tstep * factor(7,9) ! Initial vpar
            else
              rhs_ij(7) = rhs_ij(7) - visco_par * (v_x * (vpar0_x * F0**2 / BigR**2 - 2 * vpar0 * F0**2 / BigR**3  - 2 * PI * F0 * Omega_tor0_x )  & 
                                                 + v_y * (vpar0_y * F0**2 / BigR**2 -2 * PI * F0 * Omega_tor0_y) ) * BigR* xjac * tstep * factor(7,9) 
            endif
 
            if (NEO) then
              rhs_ij(7) =  rhs_ij(7)  + amu_neo_prof(ms,mt)*BB2/(Btheta2+epsil)  * factor(7,10) & ! Neoclassical
                        * v * ( r0 * (ps0_x * u0_x + ps0_y * u0_y)                &
                              + tauIC   * (ps0_x * P0_x + ps0_y * P0_y)                &
                              + aki_neo_prof(ms,mt) * tauIC * r0 * (ps0_x * T0_x + ps0_y * T0_y) - r0 * Vpar0 * Btheta2) * xjac * tstep * BigR
            endif
 
            rhs_ij_k(7) = + 0.5d0 * r0 * vpar0**2 * BB2 * F0 / BigR * v_p                     * xjac * tstep * factor(7,3) &

               - TG_NUM7 * 0.25d0 * r0 * Vpar0**2 * BB2 &
                         * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac + F0 / BigR * vpar0_p) / BigR  &
                         * (                                          + F0 / BigR * v_p)  * xjac * tstep * tstep * factor(7,6) 

            !###################################################################################################
            !#  equation 8   (neutral density equation)                                                        #
            !###################################################################################################
             rhs_ij(8) = BigR * (- Dn0x * rn0_x * v_x - Dn0y * rn0_y * v_y )   * xjac * tstep     * factor(8,1)         & ! Diffusion         
                      
                      + delta_n_convection*(                                                              &  
                      + v * BigR**2 * ( rn0_s * u0_t - rn0_t * u0_s)                              * tstep &
                      + v * 2.d0 * BigR * rn0 * u0_y                                       * xjac * tstep &
                      - v * rn0 * (vpar0_s * ps0_t - vpar0_t * ps0_s)                      * tstep &
                      - v * Vpar0 * (rn0_s * ps0_t - rn0_t * ps0_s)                        * tstep &
                      - v * F0 / BigR * Vpar0 * rn0_p                                      * xjac * tstep &
                      - v * F0 / BigR * rn0 * vpar0_p                                      * xjac * tstep &
                      )                                                                             * factor(8,2)       & ! Convection

                  - BigR * v * r0_corr * rn0_corr * Sion_T                                 * xjac * tstep * factor(8,3) &  ! Ionization sink
                  + BigR * v * r0_corr * r0_corr  * Srec_T                                 * xjac * tstep * factor(8,4) &  ! Recombination source
                  + BigR * v * source_neutral                                              * xjac * tstep * factor(8,5) &  ! External neutral source
                  - Dn_perp_num * (v_xx + v_x/Bigr + v_yy)*(rn0_xx + rn0_x/Bigr + rn0_yy)  * BigR * xjac * tstep * factor(8,6) & ! Hyper diffusitivity

                  + v * delta_g(mp,8,ms,mt) * BigR * xjac * zeta* factor(8,7) ! Variation

             rhs_ij_k(8) = BigR * ( - Dn0p * rn0_p * v_p/BigR**2)   * xjac * tstep   * factor(8,1)     
                       
            !###################################################################################################
            !#  RHS equations end                                                                                  #
            !###################################################################################################

            if (use_fft) then
              index_ij =       n_var*n_degrees*(i-1) +       n_var*(j-1) + 1
            else
              index_ij = n_tor_local*n_var*n_degrees*(i-1) + n_tor_local * n_var * (j-1) + im - n_tor_start +1 
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
              ij8 = index_ij + 7

              RHS_p(mp,ij1) = RHS_p(mp,ij1) + rhs_ij(1) * wst
              RHS_p(mp,ij2) = RHS_p(mp,ij2) + rhs_ij(2) * wst
              RHS_p(mp,ij3) = RHS_p(mp,ij3) + rhs_ij(3) * wst
              RHS_p(mp,ij4) = RHS_p(mp,ij4) + rhs_ij(4) * wst
              RHS_p(mp,ij5) = RHS_p(mp,ij5) + rhs_ij(5) * wst
              RHS_p(mp,ij6) = RHS_p(mp,ij6) + rhs_ij(6) * wst
              RHS_p(mp,ij7) = RHS_p(mp,ij7) + rhs_ij(7) * wst
              RHS_p(mp,ij8) = RHS_p(mp,ij8) + rhs_ij(8) * wst

              RHS_k(mp,ij5) = RHS_k(mp,ij5) + rhs_ij_k(5) * wst
              RHS_k(mp,ij6) = RHS_k(mp,ij6) + rhs_ij_k(6) * wst
              RHS_k(mp,ij7) = RHS_k(mp,ij7) + rhs_ij_k(7) * wst
              RHS_k(mp,ij8) = RHS_k(mp,ij8) + rhs_ij_k(8) * wst

              !--- THIS IS ONLY USED FOR DIAGNOSTIC PURPOSES-------------------------
              !--- ELM structures are re-used to plot separate terms in vtk ---------
              if ( present(get_terms) ) then
                ELM_p(mp,i_term,ij1) = ELM_p(mp,i_term,ij1) + rhs_ij(1) * wst
                ELM_p(mp,i_term,ij2) = ELM_p(mp,i_term,ij2) + rhs_ij(2) * wst
                ELM_p(mp,i_term,ij3) = ELM_p(mp,i_term,ij3) + rhs_ij(3) * wst
                ELM_p(mp,i_term,ij4) = ELM_p(mp,i_term,ij4) + rhs_ij(4) * wst
                ELM_p(mp,i_term,ij5) = ELM_p(mp,i_term,ij5) + rhs_ij(5) * wst
                ELM_p(mp,i_term,ij6) = ELM_p(mp,i_term,ij6) + rhs_ij(6) * wst
                ELM_p(mp,i_term,ij7) = ELM_p(mp,i_term,ij7) + rhs_ij(7) * wst
                ELM_p(mp,i_term,ij8) = ELM_p(mp,i_term,ij8) + rhs_ij(8) * wst

                ELM_k(mp,i_term,ij5) = ELM_k(mp,i_term,ij5) + rhs_ij_k(5) * wst
                ELM_k(mp,i_term,ij6) = ELM_k(mp,i_term,ij6) + rhs_ij_k(6) * wst
                ELM_k(mp,i_term,ij7) = ELM_k(mp,i_term,ij7) + rhs_ij_k(7) * wst
                ELM_k(mp,i_term,ij8) = ELM_k(mp,i_term,ij8) + rhs_ij_k(8) * wst
              endif
              !----------------------------------------------------------------------
            else
              ij1 = index_ij
              ij2 = index_ij + 1*n_tor_local
              ij3 = index_ij + 2*n_tor_local
              ij4 = index_ij + 3*n_tor_local
              ij5 = index_ij + 4*n_tor_local
              ij6 = index_ij + 5*n_tor_local
              ij7 = index_ij + 6*n_tor_local
              ij8 = index_ij + 7*n_tor_local

              RHS(ij1) = RHS(ij1) + (rhs_ij(1) + rhs_ij_k(1)) * wst
              RHS(ij2) = RHS(ij2) + (rhs_ij(2) + rhs_ij_k(2)) * wst
              RHS(ij3) = RHS(ij3) + (rhs_ij(3) + rhs_ij_k(3)) * wst
              RHS(ij4) = RHS(ij4) + (rhs_ij(4) + rhs_ij_k(4)) * wst
              RHS(ij5) = RHS(ij5) + (rhs_ij(5) + rhs_ij_k(5)) * wst
              RHS(ij6) = RHS(ij6) + (rhs_ij(6) + rhs_ij_k(6)) * wst
              RHS(ij7) = RHS(ij7) + (rhs_ij(7) + rhs_ij_k(7)) * wst
              RHS(ij8) = RHS(ij8) + (rhs_ij(8) + rhs_ij_k(8)) * wst

              !--- THIS IS ONLY USED FOR DIAGNOSTIC PURPOSES-------------------------
              !--- ELM structure is re-used to plot separate terms in vtk -----------             
              if ( present(get_terms) ) then
                ELM(i_term,ij1) = ELM(i_term,ij1) + (rhs_ij(1) + rhs_ij_k(1)) * wst
                ELM(i_term,ij2) = ELM(i_term,ij2) + (rhs_ij(2) + rhs_ij_k(2)) * wst
                ELM(i_term,ij3) = ELM(i_term,ij3) + (rhs_ij(3) + rhs_ij_k(3)) * wst
                ELM(i_term,ij4) = ELM(i_term,ij4) + (rhs_ij(4) + rhs_ij_k(4)) * wst
                ELM(i_term,ij5) = ELM(i_term,ij5) + (rhs_ij(5) + rhs_ij_k(5)) * wst
                ELM(i_term,ij6) = ELM(i_term,ij6) + (rhs_ij(6) + rhs_ij_k(6)) * wst
                ELM(i_term,ij7) = ELM(i_term,ij7) + (rhs_ij(7) + rhs_ij_k(7)) * wst
                ELM(i_term,ij8) = ELM(i_term,ij8) + (rhs_ij(8) + rhs_ij_k(8)) * wst
              endif
              !----------------------------------------------------------------------
            endif
            
            enddo ! max_terms_loop for diagnostic purposes

            if ( present(get_terms) ) cycle

            do k=1,n_vertex_max

              do l=1,n_degrees

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

                  u    = psi    ;  zj    = psi    ;  w    = psi    ; rho    = psi    ;  T    = psi    ; vpar    = psi   ; rhon   = psi
                  u_x  = psi_x  ;  zj_x  = psi_x  ;  w_x  = psi_x  ; rho_x  = psi_x  ;  T_x  = psi_x  ; vpar_x  = psi_x ; rhon_x = psi_x
                  u_y  = psi_y  ;  zj_y  = psi_y  ;  w_y  = psi_y  ; rho_y  = psi_y  ;  T_y  = psi_y  ; vpar_y  = psi_y ; rhon_y = psi_y
                  u_p  = psi_p  ;  zj_p  = psi_p  ;  w_p  = psi_p  ; rho_p  = psi_p  ;  T_p  = psi_p  ; vpar_p  = psi_p ; rhon_p = psi_p
                  u_s  = psi_s  ;  zj_s  = psi_s  ;  w_s  = psi_s  ; rho_s  = psi_s  ;  T_s  = psi_s  ; vpar_s  = psi_s ; rhon_s = psi_s
                  u_t  = psi_t  ;  zj_t  = psi_t  ;  w_t  = psi_t  ; rho_t  = psi_t  ;  T_t  = psi_t  ; vpar_t  = psi_t ; rhon_t = psi_t
                  u_ss = psi_ss ;  zj_ss = psi_ss ;  w_ss = psi_ss ; rho_ss = psi_ss ;  T_ss = psi_ss ; vpar_ss = psi_ss; rhon_ss = psi_ss
                  u_tt = psi_tt ;  zj_tt = psi_tt ;  w_tt = psi_tt ; rho_tt = psi_tt ;  T_tt = psi_tt ; vpar_tt = psi_tt; rhon_tt = psi_tt
                  u_st = psi_st ;  zj_st = psi_st ;  w_st = psi_st ; rho_st = psi_st ;  T_st = psi_st ; vpar_st = psi_st; rhon_st = psi_st

                  u_xx = psi_xx ;                    w_xx = psi_xx ; rho_xx = psi_xx ;  T_xx = psi_xx ; vpar_xx = psi_xx ; rhon_xx = psi_xx
                  u_yy = psi_yy ;                    w_yy = psi_yy ; rho_yy = psi_yy ;  T_yy = psi_yy ; vpar_yy = psi_yy ; rhon_yy = psi_yy
                  u_xy = psi_xy ;                    w_xy = psi_xy ; rho_xy = psi_xy ;  T_xy = psi_xy ; vpar_xy = psi_xy ; rhon_xy = psi_xy

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

                  amat(1,1) = v * psi / BigR * xjac * (1.d0 + zeta)                                                   &
                            - v * (psi_s * u0_t - psi_t * u0_s)                                       * theta * tstep &
                            + v * tauIC/(r0_corr*BB2) * F0**2/BigR**2 * (psi_s * p0_t - psi_t * p0_s) * theta * tstep &
                     - v * tauIC/(r0_corr*BB2**2) * BB2_psi * F0**2/BigR**2 * (ps0_x*p0_y - ps0_y*p0_x) * xjac * theta * tstep &
                     + v * tauIC/(r0_corr*BB2**2) * BB2_psi * F0**3/BigR**3 * p0_p           * xjac * theta * tstep

                  amat(1,2) = -  v * (ps0_s * u_t - ps0_t * u_s)                             * theta * tstep

                  amat_n(1,2) = + F0 / BigR * v * u_p * xjac                                 * theta * tstep

                  amat(1,3) = - eta_num_T * (v_x * zj_x + v_y * zj_y)                 * xjac * theta * tstep  &
                            - eta_T * v * zj / BigR                                 * xjac * theta * tstep

                  amat(1,5) = + v * tauIC/(r0_corr*BB2) * F0**2/BigR**2 * T0  * (ps0_s * rho_t - ps0_t * rho_s) * theta * tstep &
                              + v * tauIC/(r0_corr*BB2) * F0**2/BigR**2 * rho * (ps0_s * T0_t  - ps0_t * T0_s)  * theta * tstep &
                              - v * tauIC/(r0_corr*BB2) * F0**3/BigR**3 * rho * T0_p  * xjac          * theta * tstep &

                              - v * tauIC * rho /(r0_corr**2 * BB2) * F0**2/BigR**2 * (ps0_s * p0_t - ps0_t * p0_s) * theta * tstep &
                              + v * tauIC * rho /(r0_corr**2 * BB2) * F0**3/BigR**3 * p0_p * xjac                   * theta * tstep 

                  amat_n(1,5) = - v * tauIC/(r0_corr*BB2) * F0**3/BigR**3 * T0  * rho_p            * xjac * theta * tstep 

                  amat(1,6) = - deta_dT * v * T * (zj0 - current_source(ms,mt) - Jb)/ BigR * xjac         * theta * tstep &
                            + v * tauIC/(r0_corr*BB2) * F0**2/BigR**2 * r0 * (ps0_s * T_t  - ps0_t * T_s) * theta * tstep &
                            + v * tauIC/(r0_corr*BB2) * F0**2/BigR**2 * T  * (ps0_s * r0_t - ps0_t * r0_s)* theta * tstep &
                            - v * tauIC/(r0_corr*BB2) * F0**3/BigR**3 * T  * r0_p * xjac                  * theta * tstep 

                  amat_n(1,6) = - v * tauIC/(r0_corr*BB2) * F0**3/BigR**3 * r0 * T_p               * xjac * theta * tstep 

                  !###################################################################################################
                  !#  equation 2   (perpendicular momentum equation)                                                 #
                  !###################################################################################################

                  amat(2,1) = - v * (psi_s * zj0_t - psi_t * zj0_s)                          * theta * tstep

                  if (NEO) then
                    amat(2,1) = amat(2,1) &
                         - amu_neo_prof(ms,mt)*BB2/((Btheta2+epsil)**2)*(psi_x*v_x+psi_y*v_y)  &
                         * (r0 * (ps0_x*u0_x + ps0_y*u0_y) + tauIC * (ps0_x*P0_x + ps0_y*P0_y) &
                         + aki_neo_prof(ms,mt) * tauIC * r0 * (ps0_x*T0_x + ps0_y*T0_y) -r0 * Vpar0 * Btheta2) * BigR * xjac * theta * tstep &                         
                         - amu_neo_prof(ms,mt) * BB2/((Btheta2+epsil)**2) * (ps0_x*v_x + ps0_y*v_y) * (r0*(psi_x*u0_x + psi_y*u0_y) &
                              + tauIC * (psi_x*P0_x + psi_y*P0_y)                                     &
                         + aki_neo_prof(ms,mt) * tauIC * r0 * (psi_x*T0_x + psi_y*T0_y)) * BigR * xjac * theta * tstep &

                         ! ========= linearization of 1/(Btheta2**i) , i=2 or 1
                         - amu_neo_prof(ms,mt) * BB2 * (-2.d0*Btheta2_psi)/((Btheta2+epsil)**3) * (ps0_x*v_x + ps0_y*v_y)  &
                                                     * (r0 * (ps0_x*u0_x + ps0_y*u0_y) + tauIC * (ps0_x*P0_x + ps0_y*P0_y) &
                         + aki_neo_prof(ms,mt) * tauIC * r0 * (ps0_x*T0_x + ps0_y*T0_y)) * BigR * xjac * theta * tstep     &
                         + amu_neo_prof(ms,mt) * BB2 * (-Btheta2_psi)/((Btheta2+epsil)**2) * r0 * vpar0 * (ps0_x*v_x + ps0_y*v_y) &
                                               * BigR * xjac * tstep * theta 
                  endif

                  amat(2,2) = - BigR**3 * r0_corr * (v_x * u_x + v_y * u_y) * xjac * (1.d0 + zeta)                                 &
                            + r0_hat * BigR**2 * w0 * (v_s * u_t  - v_t  * u_s)                              * theta * tstep &
                            + BigR**2 * (u_x * u0_x + u_y * u0_y) * (v_x * r0_y_hat - v_y * r0_x_hat) * xjac * theta * tstep &
                            
                            + tauIC * BigR**3 * p0_y * (v_x* u_x + v_y * u_y)                         * xjac * theta * tstep &
                            
                            + v * tauIC * BigR**4 * (u_xy * (p0_xx - p0_yy) - p0_xy * (u_xx - u_yy))  * xjac * theta * tstep &
                            
                            - BigR**3 * (particle_source(ms,mt)+source_pellet) * (v_x * u_x + v_y * u_y) * xjac * theta * tstep &
                            
                            + TG_num2 * 0.25d0 * r0_hat * BigR**3 * (w0_x * u_y - w0_y * u_x)       &
                                      * ( v_x * u0_y - v_y * u0_x) * xjac * theta * tstep * tstep   &
                                      
                            + TG_num2 * 0.25d0 * r0_hat * BigR**3 * (w0_x * u0_y - w0_y * u0_x)     &
                                      * ( v_x * u_y - v_y * u_x)   * xjac * theta * tstep * tstep   &
                  
                           + (1.d0 - delta_n_convection) * (  &
                           - BigR**3 * (r0_corr*rn0_corr*Sion_T) * (v_x * u_x + v_y * u_y) * xjac * theta * tstep &
                           + BigR**3 * (r0_corr* r0_corr*Srec_T) * (v_x * u_x + v_y * u_y) * xjac * theta * tstep &
                           )
   
                  if ( NEO ) then
                    amat(2,2) = amat(2,2) &
                              - amu_neo_prof(ms,mt)*BB2/((Btheta2+epsil)**2) * (ps0_x*v_x + ps0_y*v_y) * r0 *(ps0_x*u_x + ps0_y*u_y) &
                              * BigR * xjac * theta * tstep 
                  endif
                  !---------------------------------------- NEO

                  amat(2,3)   = - v * (ps0_s * zj_t  - ps0_t * zj_s)              * theta * tstep

                  amat_n(2,3) = + F0 / BigR * v * zj_p  * xjac                    * theta * tstep

                  amat(2,4) = r0_hat * BigR**2 * w  * ( v_s * u0_t - v_t * u0_s)  * theta * tstep  &
                            + BigR * ( v_x * w_x + v_y * w_y) * visco_T  * xjac   * theta * tstep  &
                            + v * tauIC * BigR**4 * (p0_s * w_t - p0_t * w_s)     * theta * tstep  & 

                            + visco_num_T * (v_xx + v_x/BigR + v_yy)*(w_xx + w_x/BigR + w_yy) * xjac * theta * tstep    &

                            + TG_num2 * 0.25d0 * r0_hat * BigR**3 * (w_x * u0_y - w_y * u0_x)     &
                                    * ( v_x * u0_y - v_y * u0_x) * xjac * theta * tstep * tstep

                  amat(2,5) = + 0.5d0 * vv2 * (v_x * rho_y_hat - v_y * rho_x_hat)   * xjac * theta * tstep &
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
                                      * ( v_x * u0_y - v_y * u0_x) * xjac * theta * tstep * tstep        &

                        + (1.d0 - delta_n_convection) * (  &
                         - BigR**3 * (rho* rn0_corr * Sion_T)        * (v_x * u0_x + v_y * u0_y) * xjac * theta * tstep &
                         + BigR**3 * (rho * 2.d0 * r0_corr * Srec_T) * (v_x * u0_x + v_y * u0_y) * xjac * theta * tstep &
                           )
                           
                  if ( NEO ) then
                    amat(2,5) = amat(2,5) &
                         - amu_neo_prof(ms,mt)*BB2/((Btheta2+epsil)**2) * (ps0_x*v_x + ps0_y*v_y) * (rho*(ps0_x*u0_x+ps0_y*u0_y) &
                         + tauIC * (ps0_x*(rho_x*T0 + rho*T0_x) + ps0_y*(rho_y*T0 + rho*T0_y))                                   &
                         + aki_neo_prof(ms,mt) * tauIC * rho * (ps0_x*T0_x + ps0_y*T0_y) - rho * Vpar0 * Btheta2) * BigR * xjac * tstep * theta
                  endif

                  amat(2,6) = - BigR**2 * (v_s * r0_t * T   - v_t * r0_s * T)           * theta * tstep  &
                              - BigR**2 * (v_s * r0   * T_t - v_t * r0   * T_s)         * theta * tstep  &
                            + dvisco_dT * T * ( v_x * w0_x + v_y * w0_y ) * BigR * xjac * theta * tstep  &

                            + v * tauIC * BigR**4 * r0 * (T_s * w0_t - T_t * w0_s)      * theta * tstep &
                            + v * tauIC * BigR**4 * T  * (r0_s * w0_t - r0_t * w0_s)    * theta * tstep &
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
                            - dvisco_dT     * bigR * W_dia   * (v_x*T_x  + v_y*T_y )    * xjac * theta * tstep  &
                            
                           + (1 - delta_n_convection) * (  &
                           - BigR**3 * (r0_corr * rn0_corr * dSion_dT * T) * (v_x * u0_x + v_y * u0_y) * xjac * theta * tstep  &
                           + BigR**3 * (r0_corr * r0_corr  * dSrec_dT * T) * (v_x * u0_x + v_y * u0_y) * xjac * theta * tstep  &
                           )

                  !---------------------------------------- NEO
                  if ( NEO ) then
                    amat(2,6) = amat(2,6) - amu_neo_prof(ms,mt)*BB2/((Btheta2+epsil)**2)*(ps0_x*v_x + ps0_y*v_y) &
                              * (tauIC*(ps0_x*(r0_x*T+r0*T_x)+ps0_y*(r0_y*T+r0*T_y))                             &
                              + aki_neo_prof(ms,mt) *tauIC * r0 *(ps0_x*T_x + ps0_y*T_y)) * BigR * xjac * tstep * theta
                    
                    amat(2,7) = amu_neo_prof(ms,mt)*BB2/(Btheta2+epsil)*r0*vpar*(ps0_x*v_x+ps0_y*v_y) &
                              * BigR * xjac * tstep * theta 
                  endif
                  !---------------------------------------- NEO

                  amat(2,8) =  - (1.d0 - delta_n_convection) * BigR**3 * (r0 * rhon * Sion_T) * (v_x * u0_x + v_y * u0_y) * xjac * theta * tstep 

                  !###################################################################################################
                  !#  equation 3   (current definition)                                                              #
                  !###################################################################################################

                  amat(3,3) = v * zj / BigR * xjac                                
                  amat(3,1) = (v_x * psi_x + v_y * psi_y ) / BigR * xjac         

                  !###################################################################################################
                  !#  equation 4   (vorticity definition)                                                            #
                  !###################################################################################################

                  amat(4,4) =  v * w * BigR * xjac                                
                  amat(4,2) = (v_x * u_x + v_y * u_y) * BigR * xjac              

                  !###################################################################################################
                  !#  equation 5   (density equation)                                                                #
                  !###################################################################################################


                  amat(5,1) =-(D_par-D_prof) * BigR * BB2_psi/ BB2**2 * Bgrad_rho_star     * Bgrad_rho     * xjac * theta * tstep &
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


                  amat_k(5,1) = - (D_par-D_prof) * BigR * BB2_psi/ BB2**2 * Bgrad_rho_k_star * Bgrad_rho     * xjac * theta * tstep &
                                + (D_par-D_prof) * BigR / BB2             * Bgrad_rho_k_star * Bgrad_rho_psi * xjac * theta * tstep &

                              + TG_num5 * 0.25d0 / BigR * vpar0**2                                                            &
                                      * (r0_x * psi_y - r0_y * psi_x)                                                         &
                                      * (                            + F0 / BigR * v_p) * xjac * theta * tstep * tstep

                  amat(5,2) =- v * BigR**2 * ( r0_s * u_t - r0_t * u_s)                                        * theta * tstep &
                             - v * 2.d0 * BigR * r0 * u_y                                               * xjac * theta * tstep &

                             + TG_num5 * 0.25d0 * BigR**3 * (r0_x * u_y  - r0_y * u_x)                                      &
                                                          * ( v_x * u0_y - v_y  * u0_x) * xjac * theta * tstep * tstep      &
                             + TG_num5 * 0.25d0 * BigR**3 * (r0_x * u0_y - r0_y * u0_x)                                     &
                                                          * ( v_x * u_y  - v_y  * u_x)  * xjac * theta * tstep * tstep 

                  amat(5,5) = v * rho * BigR * (1.d0 + zeta)                                            * xjac   &
                          - v * BigR**2 * ( rho_s * u0_t - rho_t * u0_s)                                       * theta * tstep &
                          - v * 2.d0 * BigR * rho * u0_y                                                * xjac * theta * tstep &
                          + (D_par-D_prof) * BigR / BB2 * Bgrad_rho_star * Bgrad_rho_rho                * xjac * theta * tstep &
                          + D_prof * BigR  * (v_x*rho_x + v_y*rho_y )                                   * xjac * theta * tstep &
                          + v * Vpar0 * (rho_s * ps0_t - rho_t * ps0_s)                                        * theta * tstep &
                          + v * rho * (vpar0_s * ps0_t - vpar0_t * ps0_s)                                      * theta * tstep &
                          + v * rho * F0 / BigR * vpar0_p                                               * xjac * theta * tstep &

                          - v * 2.d0 * tauIC * (rho_y * T0 + rho*T0_y) * BigR                           * xjac * theta * tstep &

                         - v * rho * rn0_corr       * BigR * Sion_T * xjac * theta * tstep &
                         + v * rho * 2.d0 * r0_corr * BigR * Srec_T * xjac * theta * tstep &

                          + D_perp_num * (v_xx + v_x/BigR + v_yy)*(rho_xx + rho_x/BigR + rho_yy)  * BigR * xjac * theta * tstep &

                          + TG_num5 * 0.25d0 * BigR**3 * (rho_x * u0_y - rho_y * u0_x)                                &
                                                       * ( v_x  * u0_y - v_y   * u0_x) * xjac * theta * tstep * tstep &

                          + TG_num5 * 0.25d0 / BigR * vpar0**2                                                        &
                                    * (rho_x * ps0_y - rho_y * ps0_x )                             &
                                    * ( v_x * ps0_y -  v_y * ps0_x   ) * xjac * theta * tstep * tstep

                  amat_k(5,5) = + (D_par-D_prof) * BigR / BB2 * Bgrad_rho_k_star * Bgrad_rho_rho          * xjac * theta * tstep &
 
                         + TG_num5 * 0.25d0 / BigR * vpar0**2                                                        &
                                    * (rho_x * ps0_y - rho_y * ps0_x                  )                              &
                                    * (                              + F0 / BigR * v_p) * xjac * theta * tstep * tstep

                  amat_n(5,5) = + (D_par-D_prof) * BigR / BB2 * Bgrad_rho_star   * Bgrad_rho_rho_n        * xjac * theta * tstep &
                              + v * F0 / BigR * Vpar0 * rho_p                                             * xjac * theta * tstep &

                          + TG_num5 * 0.25d0 / BigR * vpar0**2                                                        &
                                    * (                              + F0 / BigR * rho_p)                             &
                                    * ( v_x * ps0_y -  v_y * ps0_x                      ) * xjac * theta * tstep * tstep

                  amat_kn(5,5) = + (D_par-D_prof) * BigR / BB2 * Bgrad_rho_k_star * Bgrad_rho_rho_n       * xjac * theta * tstep &
                                 + D_prof * BigR  * ( v_p*rho_p /BigR**2 )                                * xjac * theta * tstep &

                          + TG_num5 * 0.25d0 / BigR * vpar0**2                                                        &
                                    * ( + F0 / BigR * rho_p)                                                          &
                                    * ( + F0 / BigR * v_p) * xjac * theta * tstep * tstep

                  amat(5,6) = - v * 2.d0 * tauIC * (T_y * r0 + T*r0_y) * BigR                             * xjac * theta * tstep &
                            - v * BigR * r0_corr * rn0_corr * dSion_dT * T                                * xjac * theta * tstep &
                            + v * BigR * r0_corr * r0_corr *  dSrec_dT * T                                * xjac * theta * tstep  

                  amat(5,7) = + v * F0 / BigR * Vpar * r0_p                                               * xjac * theta * tstep &
                              + v * Vpar * (r0_s * ps0_t - r0_t * ps0_s)                                         * theta * tstep &
                              + v * r0 * (vpar_s * ps0_t - vpar_t * ps0_s)                                       * theta * tstep &

                              + TG_num5 * 0.25d0 / BigR * 2.d0*vpar0*vpar                                   &
                                   * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                       &
                                   * ( v_x * ps0_y -  v_y * ps0_x                   ) * xjac * theta * tstep * tstep 

                  amat_k(5,7) = + TG_num5 * 0.25d0 / BigR * 2.d0*vpar0*vpar                                 &
                                   * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                       &
                                   * (                            + F0 / BigR * v_p) * xjac * theta * tstep * tstep 

                  amat_n(5,7) = + v * r0 * F0 / BigR * vpar_p                                           * xjac * theta * tstep

                  amat(5,8) = - BigR * v * r0_corr * Sion_T * rhon                                      * xjac * theta * tstep

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
                  
                  amat(6,1) = - (ZKpar_T-ZK_prof) * BigR * BB2_psi / BB2**2 * Bgrad_T_star   * Bgrad_T       * xjac * theta * tstep &
                              + (ZKpar_T-ZK_prof) * BigR / BB2              * Bgrad_T_star_psi * Bgrad_T     * xjac * theta * tstep &
                              + (ZKpar_T-ZK_prof) * BigR / BB2              * Bgrad_T_star     * Bgrad_T_psi * xjac * theta * tstep &

                            + v * r0 * Vpar0 * (T0_s * psi_t - T0_t * psi_s)                                     * theta * tstep &
                            + v * T0 * Vpar0 * (r0_s * psi_t - r0_t * psi_s)                                     * theta * tstep &
                            + v * r0 * GAMMA * T0 * (vpar0_s * psi_t - vpar0_t * psi_s)                          * theta * tstep &
!===================== Additional terms from friction terms============
                            - v * ((GAMMA - 1.) / BigR) * vpar0**2 * (psi_x * ps0_x + psi_y * ps0_y)&
                                * (r0_corr*rn0*Sion_T)                                                    * xjac * theta * tstep &
!==============================End of friction terms=================

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

                  amat_k(6,1) = - (ZKpar_T-ZK_prof) * BigR * BB2_psi / BB2**2 * Bgrad_T_k_star * Bgrad_T     * xjac * theta * tstep &
                                + (ZKpar_T-ZK_prof) * BigR / BB2              * Bgrad_T_k_star * Bgrad_T_psi * xjac * theta * tstep &
  
                        + TG_num6 * 0.25d0 / BigR * vpar0**2                                                       &
                                  * T0 * (r0_x * psi_y - r0_y * psi_x)                                             &
                                  * (                            + F0 / BigR * v_p) * xjac * theta * tstep * tstep &
                        + TG_num6 * 0.25d0 / BigR * vpar0**2                                                       &
                                  * r0 * (T0_x * psi_y - T0_y * psi_x)                                             &
                                  * (                            + F0 / BigR * v_p) * xjac * theta * tstep * tstep


                  amat(6,2) = - v * r0 * BigR**2 * ( T0_x * u_y - T0_y * u_x)           * xjac * theta * tstep &
                              - v * T0 * BigR**2 * ( r0_x * u_y - r0_y * u_x)           * xjac * theta * tstep &
                              - v * r0 * 2.d0* GAMMA * BigR * T0 * u_y                  * xjac * theta * tstep &
!===================== Additional terms from friction terms============
                              - v * BigR**3 * (GAMMA - 1.) * (u_x * u0_x + u_y * u0_y)  &
                                  * (r0_corr*rn0*Sion_T)                                * xjac * theta * tstep &
!==============================End of friction terms=================

                         + TG_num6 * 0.25d0 * BigR**2 * T0* (r0_x * u_y - r0_y * u_x)                &
                                            * ( v_x * u0_y - v_y * u0_x) * xjac * theta*tstep*tstep  &
                         + TG_num6 * 0.25d0 * BigR**2 * r0* (T0_x * u_y - T0_y * u_x)                &
                                            * ( v_x * u0_y - v_y * u0_x) * xjac * theta*tstep*tstep  &
                         + TG_num6 * 0.25d0 * BigR**2 * T0* (r0_x * u0_y - r0_y * u0_x)              &
                                            * ( v_x * u_y - v_y * u_x) * xjac * theta*tstep*tstep    &
                         + TG_num6 * 0.25d0 * BigR**2 * r0* (T0_x * u0_y - T0_y * u0_x)              &
                                            * ( v_x * u_y - v_y * u_x) * xjac * theta*tstep*tstep 

                  amat(6,3) = - v * (gamma-1.d0) * eta_T_ohm * 2.d0 * zj * zj0/(BigR**2.d0) * BigR * xjac * theta * tstep

                  amat(6,5) = v * rho * T0   * BigR * xjac * (1.d0 + zeta)     &
                            - v * rho * BigR**2 * ( T0_s * u0_t - T0_t * u0_s)                          * theta * tstep &
                            - v * T0  * BigR**2 * ( rho_s * u0_t - rho_t * u0_s)                        * theta * tstep &
                            - v * rho * 2.d0* GAMMA * BigR * T0 * u0_y                           * xjac * theta * tstep &
                            + v * rho * F0 / BigR * Vpar0 * T0_p                                 * xjac * theta * tstep &

                            + v * rho * Vpar0 * (T0_s  * ps0_t - T0_t * ps0_s)                          * theta * tstep &
                            + v * T0  * Vpar0 * (rho_s * ps0_t - rho_t * ps0_s)                         * theta * tstep & 

                            + v * rho * GAMMA * T0 * (vpar0_s * ps0_t - vpar0_t * ps0_s)        * theta * tstep &
                            + v * rho * GAMMA * T0 * F0 / BigR * vpar0_p                 * xjac * theta * tstep &

                           + v * BigR * rho * rn0_corr * ksiion * Sion_T                      * xjac * theta * tstep &
                           + v * BigR * rho * rn0_corr * LradDrays_T                          * xjac * theta * tstep &
                           + v * BigR * rho * 2d0 * r0_corr * LradDcont_T                * xjac * theta * tstep &
                           + v * BigR * rho * frad_bg                                    * xjac * theta * tstep &
!===================== Additional terms from friction terms============
                            - v * BigR * ((GAMMA - 1.)/2.) * vpar0**2 * BB2 * (rho*rn0*Sion_T) * xjac * theta * tstep &
                            - v * BigR * ((GAMMA - 1.)/2.) * vv2            * (rho*rn0*Sion_T) * xjac * theta * tstep &
!==============================End of friction terms=================

                         + TG_num6 * 0.25d0 * BigR**2 * T0* (rho_x * u0_y - rho_y * u0_x)      &
                                   * ( v_x * u0_y - v_y * u0_x) * xjac * theta*tstep*tstep     &
                         + TG_num6 * 0.25d0 * BigR**2 * rho * (T0_x * u0_y - T0_y * u0_x)      &
                                   * ( v_x * u0_y - v_y * u0_x) * xjac* theta*tstep*tstep      &
                         + TG_num6 * 0.25d0 / BigR * vpar0**2                                  &
                                   * T0 * (rho_x * ps0_y - rho_y * ps0_x )                     &
                                   * ( v_x * ps0_y -  v_y * ps0_x ) * xjac * theta * tstep * tstep&
                         + TG_num6 * 0.25d0 / BigR * vpar0**2 &
                                   * rho * (T0_x * ps0_y - T0_y * ps0_x + F0 / BigR * T0_p)                    &
                                   * ( v_x * ps0_y -  v_y * ps0_x ) * xjac * theta * tstep * tstep

                  amat_n(6,5) = + v * T0  * F0 / BigR * Vpar0 * rho_p                * xjac * theta * tstep    &

                         + TG_num6 * 0.25d0 / BigR * vpar0**2 &
                                   * T0 * (                              + F0 / BigR * rho_p)                  &
                                   * ( v_x * ps0_y -  v_y * ps0_x                  ) * xjac * theta * tstep * tstep

                  amat_k(6,5) = + TG_num6 * 0.25d0 / BigR * vpar0**2 &
                                   * T0 * (rho_x * ps0_y - rho_y * ps0_x                    )                  &
                                   * (                            + F0 / BigR * v_p) * xjac * theta * tstep * tstep&
                              + TG_num6 * 0.25d0 / BigR * vpar0**2 &
                                   * rho * (T0_x * ps0_y - T0_y * ps0_x + F0 / BigR * T0_p)                    &
                                   * (                            + F0 / BigR * v_p) * xjac * theta * tstep * tstep

                  amat_kn(6,5) = + TG_num6 * 0.25d0 / BigR * vpar0**2                &
                                   * T0 * (+ F0 / BigR * rho_p)                      &
                                   * (     + F0 / BigR * v_p) * xjac * theta * tstep * tstep


                  amat(6,6) = v * r0_corr * T   * BigR * xjac * (1.d0 + zeta)     &
                            - v * r0 * BigR**2 * ( T_s  * u0_t - T_t  * u0_s)               * theta * tstep &
                            - v * T  * BigR**2 * ( r0_s * u0_t - r0_t * u0_s)               * theta * tstep &

                            - v * r0 * 2.d0* GAMMA * BigR * T * u0_y                 * xjac * theta * tstep &

                            + v * r0 * Vpar0 * (T_s  * ps0_t - T_t  * ps0_s)                * theta * tstep &
                            + v * T  * Vpar0 * (r0_s * ps0_t - r0_t * ps0_s)                * theta * tstep & 

                            + v * r0 * GAMMA * T * (vpar0_s * ps0_t - vpar0_t * ps0_s)      * theta * tstep &
                            + v * r0 * GAMMA * T * F0 / BigR * vpar0_p               * xjac * theta * tstep &

                            + v * T * F0 / BigR * Vpar0 * r0_p                       * xjac * theta * tstep &

                            + (ZKpar_T-ZK_prof) * BigR / BB2 * Bgrad_T_star * Bgrad_T_T * xjac * theta * tstep &
                            + ZK_prof * BigR * ( v_x*T_x + v_y*T_y )                    * xjac * theta * tstep &

                            + dZKpar_dT * T * BigR / BB2 * Bgrad_T_star * Bgrad_T           * xjac * theta * tstep &
  
                            + ZK_perp_num * (v_xx + v_x/BigR + v_yy)*(T_xx + T_x/BigR + T_yy) * BigR * xjac * theta * tstep &

                            -v * T * (gamma-1.d0) * deta_dT_ohm * (zj0 / BigR)**2.d0 * BigR * xjac * theta * tstep &

                            + v * BigR * r0_corr * rn0_corr * ksiion * dSion_dT * T         * xjac * theta * tstep &

                            + v * BigR * T * r0_corr * rn0_corr * dLradDrays_dT             * xjac * theta * tstep &
                            + v * BigR * T * r0_corr * r0_corr  * dLradDcont_dT             * xjac * theta * tstep &
                            + v * BigR * T * r0_corr * dfrad_bg_dT                          * xjac * theta * tstep &
!===================== Additional terms from friction terms============
                            - v * BigR * ((GAMMA - 1.)/2.) * vpar0**2 * BB2 &
                                * (r0_corr*rn0*dSion_dT) * T * xjac * theta * tstep &
                            - v * BigR * ((GAMMA - 1.)/2.) * vv2 &
                                * (r0_corr*rn0*dSion_dT) * T * xjac * theta * tstep &
!==============================End of friction terms=================

                            + TG_num6 * 0.25d0 * BigR**2 * T* (r0_x * u0_y - r0_y * u0_x)         &
                                      * ( v_x * u0_y - v_y * u0_x) * xjac * theta * tstep * tstep &
                            + TG_num6 * 0.25d0 * BigR**2 * r0* (T_x * u0_y - T_y * u0_x)          &
                                      * ( v_x * u0_y - v_y * u0_x) * xjac * theta * tstep * tstep &
                            + TG_num6 * 0.25d0 / BigR * vpar0**2                                                       &
                                      * T * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                           &
                                      * ( v_x * ps0_y -  v_y * ps0_x                  ) * xjac * theta * tstep * tstep &
                            + TG_num6 * 0.25d0 / BigR * vpar0**2                                                       &
                                      * r0 * (T_x * ps0_y - T_y * ps0_x               )                                &
                                      * ( v_x * ps0_y -  v_y * ps0_x                  ) * xjac * theta * tstep * tstep 

                  amat_k(6,6) = + (ZKpar_T-ZK_prof) * BigR / BB2 * Bgrad_T_k_star * Bgrad_T_T * xjac * theta * tstep  &
                                + dZKpar_dT * T     * BigR / BB2 * Bgrad_T_k_star * Bgrad_T   * xjac * theta * tstep  &

                              + TG_num6 * 0.25d0 / BigR * vpar0**2                                                    &
                                  * T * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                              &
                                  * (                                + F0 / BigR * v_p) * xjac * theta * tstep * tstep&
                              + TG_num6 * 0.25d0 / BigR * vpar0**2                                                    &
                                  * r0 * (T_x * ps0_y - T_y * ps0_x                  )                                &
                                  * (                                + F0 / BigR * v_p) * xjac * theta * tstep * tstep

                  amat_n(6,6) = + (ZKpar_T-ZK_prof) * BigR / BB2 * Bgrad_T_star   * Bgrad_T_T_n  * xjac * theta * tstep &

                              + v * r0 * F0 / BigR * Vpar0 * T_p                               * xjac * theta * tstep &
  
                              + TG_num6 * 0.25d0 / BigR * vpar0**2                       & 
                                * r0 * ( + F0 / BigR * T_p) * ( v_x * ps0_y -  v_y * ps0_x ) * xjac * theta * tstep * tstep

                  amat_kn(6,6) = + (ZKpar_T-ZK_prof) * BigR / BB2 * Bgrad_T_k_star * Bgrad_T_T_n * xjac * theta * tstep &
                                 + ZK_prof * BigR   * (v_p*T_p /BigR**2 )                        * xjac * theta * tstep &

                              + TG_num6 * 0.25d0 / BigR * vpar0**2 &
                                * r0 * ( + F0 / BigR * T_p) * ( + F0 / BigR * v_p)          * xjac * theta * tstep * tstep

                  amat(6,7) = + v * r0 * F0 / BigR * Vpar * T0_p                            * xjac * theta * tstep &
                            + v * T0 * F0 / BigR * Vpar * r0_p                              * xjac * theta * tstep &

                            + v * r0 * Vpar * (T0_s * ps0_t - T0_t * ps0_s)                        * theta * tstep &
                            + v * T0 * Vpar * (r0_s * ps0_t - r0_t * ps0_s)                        * theta * tstep & 

                            + v * r0 * GAMMA * T0 * (vpar_s * ps0_t - vpar_t * ps0_s)       * theta * tstep        &
!===================== Additional terms from friction terms============
                            - v * BigR *(GAMMA - 1.) * vpar0 * Vpar * BB2 * (r0_corr*rn0*Sion_T) * xjac * theta * tstep &
!==============================End of friction terms=================
  
                            + TG_num6 * 0.25d0 / BigR * 2.d0 * vpar0*vpar &
                                  * T0 * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                          &
                                  * ( v_x * ps0_y -  v_y * ps0_x                        ) * xjac * theta * tstep * tstep &
                            + TG_num6 * 0.25d0 / BigR * 2.d0 * vpar0*vpar &
                                  * r0 * (T0_x * ps0_y - T0_y * ps0_x + F0 / BigR * T0_p)                          &
                                  * ( v_x * ps0_y -  v_y * ps0_x                        ) * xjac * theta * tstep * tstep 

                  amat_k(6,7) =  &
                        + TG_num6 * 0.25d0 / BigR * 2.d0 * vpar0*vpar &
                                  * T0 * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                          &
                                  * (                            + F0 / BigR * v_p) * xjac * theta * tstep * tstep &
                        + TG_num6 * 0.25d0 / BigR * 2.d0 * vpar0*vpar &
                                  * r0 * (T0_x * ps0_y - T0_y * ps0_x + F0 / BigR * T0_p)                          &
                                  * (                            + F0 / BigR * v_p) * xjac * theta * tstep * tstep 

                  amat_n(6,7) = + v * r0 * GAMMA * T0 * F0 / BigR * vpar_p          * xjac * theta * tstep

                  amat(6,8) = + v * BigR * r0_corr * rhon * ksiion * Sion_T         * xjac * theta * tstep &
                              + v * BigR * rhon * r0_corr * LradDrays_T             * xjac * theta * tstep & 
!===================== Additional terms from friction terms============
                              - v * BigR * ((GAMMA - 1.)/2.) * vpar0**2 * BB2 * (r0_corr*rhon*Sion_T) * xjac * theta * tstep &
                              - v * BigR * ((GAMMA - 1.)/2.) * vv2            * (r0_corr*rhon*Sion_T) * xjac * theta * tstep 
!==============================End of friction terms=================

                  !###################################################################################################
                  !#  equation 7   (parallel velocity equation)                                                      #
                  !###################################################################################################

                  amat(7,1) = v * r0 * vpar0 / BigR * (ps0_x * psi_x + ps0_y * psi_y) * xjac * (1.d0 + zeta) &

                            + v * (P0_s * psi_t - P0_t * psi_s)                                       * theta * tstep &

                            + 0.5d0 * r0 * vpar0**2 * BB2     * (psi_x * v_y  - psi_y * v_x)   * xjac * theta * tstep &
                            + 0.5d0 * r0 * vpar0**2 * BB2_psi * (ps0_x * v_y  - ps0_y * v_x)   * xjac * theta * tstep &
                            + 0.5d0 * v  * vpar0**2 * BB2     * (psi_x * r0_y - psi_y * r0_x)  * xjac * theta * tstep &
                            + 0.5d0 * v  * vpar0**2 * BB2_psi * (ps0_x * r0_y - ps0_y * r0_x)  * xjac * theta * tstep &
                            - 0.5d0 * v  * vpar0**2 * BB2_psi * F0 / BigR * r0_p               * xjac * theta * tstep &

                            + v*(particle_source(ms,mt) + source_pellet)*vpar0* BB2_psi * BigR * xjac * theta * tstep &
  
                            + (1.d0 - delta_n_convection) * (  &  
                              + v *(r0_corr * rn0_corr * Sion_T) * vpar0 * BB2_psi * BigR    * xjac * theta * tstep  &
                              - v *(r0_corr * r0_corr  * Srec_T) * vpar0 * BB2_psi * BigR    * xjac * theta * tstep  &
                              ) &

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
                    amat(7,1) = amat(7,1)  - visco_par * (v_x * Vt_x_psi   + v_y * Vt_y_psi) * BigR * xjac * theta * tstep      
                  else
                    amat(7,1) = amat(7,1)  - visco_par * 2.d0 * PI * F0 * (v_x * Omega_tor_x_psi + v_y * Omega_tor_y_psi) * BigR * xjac * theta * tstep 
                  endif

                  amat_k(7,1) = - 0.5d0 * r0 * vpar0**2 * BB2_psi * F0 / BigR * v_p    * xjac * theta * tstep 

                  !---------------------------------------- NEO
                  if ( NEO ) then
                    amat(7,1) = amat(7,1) &
                         - v * amu_neo_prof(ms,mt)*BB2/(Btheta2+epsil)*(r0 * (psi_x*u0_x + psi_y*u0_y) + tauIC*(psi_x*P0_x + psi_y*P0_y) &
                         + aki_neo_prof(ms,mt) * tauIC * r0 * (psi_x*T0_x + psi_y*T0_y)) * BigR * xjac * theta * tstep                 &
                         - v * amu_neo_prof(ms,mt) * (-Btheta2_psi)*BB2/(Btheta2**2) * (r0*(ps0_x*u0_x+ps0_y*u0_y)                     &
                                                      + tauIC*(ps0_x*P0_x+ps0_y*P0_y)                                                  &                     
                         + aki_neo_prof(ms,mt) * tauIC * r0 * (ps0_x*T0_x + ps0_y*T0_y)) * BigR * xjac * theta * tstep

                    amat(7,2) = -v * amu_neo_prof(ms,mt)*BB2/(Btheta2+epsil) * r0 * (ps0_x*u_x + ps0_y*u_y) * BigR * xjac * theta * tstep 
                  endif
                  !---------------------------------------- NEO

                  amat(7,5) = + v * (rho_s * T0 * ps0_t - rho_t * T0 * ps0_s)               * theta * tstep &
                            + v * (rho * T0_s * ps0_t - rho * T0_t * ps0_s)                 * theta * tstep &
                            + v * F0 / BigR * rho * T0_p                             * xjac * theta * tstep &

                            + 0.5d0 * rho * vpar0**2 * BB2 * (ps0_s * v_t - ps0_t * v_s)    * theta * tstep &
                            + 0.5d0 * v   * vpar0**2 * BB2 * (ps0_s * rho_t - ps0_t * rho_s)* theta * tstep &

                            + (1.d0 - delta_n_convection) * (  &

                            + v *(rho * rn0_corr       * Sion_T) * vpar0 * BB2 * BigR         * xjac * theta * tstep &
                            - v *(2.d0 * r0_corr * rho * Srec_T) * vpar0 * BB2 * BigR         * xjac * theta * tstep &
                            ) &

                            + TG_NUM7 * 0.25d0 * rho * Vpar0**2 * BB2 &
                                      * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac + F0 / BigR * vpar0_p) / BigR  &
                                      * (-(ps0_s * v_t     - ps0_t * v_s)    /xjac          )  * xjac * theta * tstep*tstep &

                            + TG_NUM7 * 0.25d0 * v * Vpar0**2 * BB2 &
                                      * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac + F0 / BigR * vpar0_p) / BigR  &
                                      * (-(ps0_s * rho_t   - ps0_t * rho_s)  /xjac           ) * xjac * theta * tstep*tstep

                  amat_k(7,5) = - 0.5d0 * rho * vpar0**2 * BB2 * F0 / BigR * v_p       * xjac * theta * tstep &

                              + TG_NUM7 * 0.25d0 * rho * Vpar0**2 * BB2 &
                                        * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac + F0 / BigR * vpar0_p) / BigR  &
                                        * (                                          + F0 / BigR * v_p)  * xjac * theta * tstep*tstep

                  amat_n(7,5) = + v * F0 / BigR * rho_p * T0                           * xjac * theta * tstep &
                              - 0.5d0 * v   * vpar0**2 * BB2 * F0 / BigR * rho_p       * xjac * theta * tstep &

                              + TG_NUM7 * 0.25d0 * v * Vpar0**2 * BB2 &
                                        * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac + F0 / BigR * vpar0_p) / BigR  &
                                        * (                                          + F0 / BigR * rho_p)* xjac * theta * tstep*tstep
                  
                  amat(7,6) = + v * (T_s * r0 * ps0_t - T_t * r0 * ps0_s)                   * theta * tstep &
                            + v * (T * r0_s * ps0_t - T * r0_t * ps0_s)                     * theta * tstep &
                            + v * F0 / BigR * T * r0_p                               * xjac * theta * tstep &
                        
                            + (1.d0 - delta_n_convection) * (  &
                              + v *(r0_corr * rn0_corr * dSion_dT * T) * vpar0 * BB2 * BigR           * xjac * theta * tstep &
                              - v *(r0_corr * r0_corr  * dSrec_dT * T) * vpar0 * BB2 * BigR           * xjac * theta * tstep &
                               )                       * xjac * theta * tstep

                  amat_n(7,6) = + v * F0 / BigR * T_p * R0                           * xjac * theta * tstep

                  amat(7,7) = v * Vpar * r0_corr * F0**2 / BigR * xjac * (1.d0 + zeta) &

                          + v*(particle_source(ms,mt) + source_pellet)*vpar*BB2 * BigR * xjac * theta * tstep &

                          + r0 * vpar0 * vpar * BB2 * (ps0_s * v_t - ps0_t * v_s)             * theta * tstep &
                          + v  * vpar0 * vpar * BB2 * (ps0_s * r0_t - ps0_t * r0_s)           * theta * tstep &
                          - v  * vpar0 * vpar * BB2 * F0 / BigR * r0_p                 * xjac * theta * tstep &

                          + (1.d0 - delta_n_convection) * (  &

                            + v *(r0_corr * rn0_corr * Sion_T) * vpar * BB2 * BigR               * xjac * theta * tstep   &
                            - v *(r0_corr * r0_corr  * Srec_T) * vpar * BB2 * BigR               * xjac * theta * tstep   &
                            ) &

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
                    amat(7,7) = amat(7,7) + visco_par * (v_x * Vpar_x + v_y * Vpar_y) * BigR        * xjac  * theta * tstep 
                  else
                    amat(7,7) = amat(7,7) + visco_par * F0**2 / BigR**2 * (v_x * (Vpar_x - 2*vpar/BigR) + v_y * Vpar_y) * BigR * xjac  * theta * tstep 
                  endif

                  amat_k(7,7) = - r0 * vpar0 * vpar * BB2 * F0 / BigR * v_p               * xjac * theta * tstep             &
     
                            + TG_NUM7 * 0.5d0 * r0 * Vpar * Vpar0 * BB2 &
                                      * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac + F0 / BigR * vpar0_p) / BigR                     &
                                      * (                                          + F0 / BigR * v_p)  * xjac * theta * tstep*tstep  &
                            + TG_NUM7 * 0.25d0 * r0 * Vpar0**2 * BB2 &
                                      * (-(ps0_s * vpar_t - ps0_t * vpar_s)/xjac                  ) / BigR                           &
                                      * (                                        + F0 / BigR * v_p)  * xjac * theta * tstep*tstep

                  amat_n(7,7) = &

                            + TG_NUM7 * 0.25d0 * r0 * Vpar0**2 * BB2 &
                                      * (                                        + F0 / BigR * vpar_p) / BigR                        &
                                      * (-(ps0_s * v_t    - ps0_t * v_s)   /xjac                     )  * xjac * theta * tstep*tstep &
                            + TG_NUM7 * 0.25d0 * v * Vpar0**2 * BB2 &
                                      * (                                        + F0 / BigR * vpar_p) / BigR                        &
                                      * (-(ps0_s * r0_t   - ps0_t * r0_s)  /xjac + F0 / BigR * r0_p)  * xjac * theta * tstep*tstep 

                  amat_kn(7,7) = &
                             + TG_NUM7 * 0.25d0 * r0 * Vpar0**2 * BB2 &
                                       * (                                        + F0 / BigR * vpar_p) / BigR                        &
                                       * (                                        + F0 / BigR * v_p)  * xjac * theta * tstep*tstep

                  if ( NEO ) then
                    amat(7,5) = amat(7,5) - v * amu_neo_prof(ms,mt)*BB2/(Btheta2+epsil) &
                          * (rho * (ps0_x*u0_x + ps0_y*u0_y) + tauIC*(ps0_x * (rho_x*T0 + rho*T0_x) + ps0_y*(rho_y*T0 + rho*T0_y)) &
                          + aki_neo_prof(ms,mt) * tauIC * rho*(ps0_x*T0_x + ps0_y*T0_y) - rho * Vpar0 * Btheta2) * BigR * xjac * tstep * theta
                    
                    amat(7,6) = amat(7,6) -v*amu_neo_prof(ms,mt)*BB2/(Btheta2+epsil)           &
                              * (tauIC * (ps0_x * (r0_x*T + r0*T_x) + ps0_y*(r0_y*T + r0*T_y)) &
                              + aki_neo_prof(ms,mt) * tauIC * r0 * (ps0_x*T_x + ps0_y*T_y)) * BigR * xjac * tstep * theta
                  
                    amat(7,7) = amat(7,7) + v * amu_neo_prof(ms,mt) * BB2 * r0 * vpar * BigR * xjac * tstep * theta 
                  endif

                  amat(7,8) = + (1.d0 - delta_n_convection) * v *(r0_corr * rhon * Sion_T) * vpar0 * BB2 * BigR         * xjac * theta * tstep

                  !###################################################################################################
                  !#  equation 8  (neutral density equation)                                                         #
                  !###################################################################################################
              
                  amat(8,1) = + delta_n_convection*(                                                                                   &
                                + v * rn0   * (vpar0_s * psi_t - vpar0_t * psi_s)                    * theta * tstep &
                                + v * Vpar0 * (rn0_s   * psi_t - rn0_t   * psi_s)                    * theta * tstep )

                  amat(8,2) = delta_n_convection*(                                                                                     &
                                + v * BigR**2 * ( rn0_s * u_t - rn0_t * u_s)                         * theta * tstep &
                                + v * 2.d0 * BigR * rn0 * u_y                                 * xjac * theta * tstep )

                  amat(8,5) = + BigR * v * rn0_corr * Sion_T * rho                            * xjac * theta * tstep &
                              - BigR * v * 2d0 * r0_corr * rho * Srec_T                       * xjac * theta * tstep 

               ! We do not include the term coming from div(rhon * v_star_i) because they are prop. to rho_n/rho, because they may cause problems
               ! in areas where rho is small.                 
         
                  amat(8,6) = + BigR * v * r0_corr * rn0_corr * dSion_dT * T                  * xjac * theta * tstep &
                              - BigR * v * r0_corr * r0_corr  * dSrec_dT * T                  * xjac * theta * tstep       
                 
                  amat(8,7) = + delta_n_convection * ( v * F0 / BigR * Vpar * rn0_p           * xjac * theta * tstep &
                                + v * Vpar * (rn0_s * ps0_t - rn0_t * ps0_s)                         * theta * tstep &
                                + v * rn0 * (vpar_s * ps0_t - vpar_t * ps0_s)                        * theta * tstep )

                  amat_n(8,7) = + delta_n_convection * ( &
                                + v * rn0 * F0 / BigR * vpar_p                                * xjac * theta * tstep ) 

                  amat(8,8) = + v * rhon * BigR * xjac * (1.d0 + zeta)   &

                            + delta_n_convection*(                                                                                     &
                              - v * BigR**2 * ( rhon_s * u0_t - rhon_t * u0_s)                       * theta * tstep &
                              - v * 2.d0 * BigR * rhon * u0_y                                 * xjac * theta * tstep &
                              + v * rhon * (vpar0_s * ps0_t - vpar0_t * ps0_s)                       * theta * tstep &
                              + v * Vpar0 * (rhon_s * ps0_t - rhon_t * ps0_s)                        * theta * tstep &
                              + v * F0 / BigR * rhon * vpar0_p                                * xjac * theta * tstep ) &
                                   
                   + BigR * (Dn0x * rhon_x * v_x + Dn0y * rhon_y * v_y)                       * xjac * theta * tstep &   
                   + BigR * v * r0_corr * rhon* Sion_T                                        * xjac * theta * tstep &
                   + Dn_perp_num * (v_xx + v_x/BigR + v_yy)*(rhon_xx + rhon_x/BigR + rhon_yy) * BigR * xjac * theta * tstep 

                   amat_n(8,8) = + delta_n_convection*(                                                             &
                               + v * F0 / BigR * Vpar0 * rhon_p                                * xjac * theta * tstep )

                   amat_kn(8,8) = + BigR * ( + Dn0p * rhon_p * v_p/BigR**2)                   * xjac * theta * tstep    
                   
                  !###################################################################################################
                  !# end equations                                                                                   #
                  !###################################################################################################

                  ! --- Fill up the matrix
                  if (use_fft) then

                    index_kl = n_var*n_degrees*(k-1) + n_var*(l-1) + 1
 
                    do kl = 1, n_var
                      do ij = 1, n_var

                        ELM_p(mp,index_kl+kl-1,ij)  =  ELM_p(mp,index_kl+kl-1,ij) + wst * amat(ij,kl)
                        ELM_n(mp,index_kl+kl-1,ij)  =  ELM_n(mp,index_kl+kl-1,ij) + wst * amat_n(ij,kl)
                        ELM_k(mp,index_kl+kl-1,ij)  =  ELM_k(mp,index_kl+kl-1,ij) + wst * amat_k(ij,kl)
                        ELM_kn(mp,index_kl+kl-1,ij) =  ELM_kn(mp,index_kl+kl-1,ij) + wst * amat_kn(ij,kl)
                      
                      enddo
                    enddo

                  else

                    index_kl = n_tor_local*n_var*n_degrees*(k-1) + n_tor_local*n_var*(l-1) + in - n_tor_start +1

                    do kl = 1, n_var
                      do ij = 1, n_var

                        ELM(index_ij+(ij-1)*(n_tor_local),index_kl+(kl-1)*(n_tor_local)) = &
                        ELM(index_ij+(ij-1)*(n_tor_local),index_kl+(kl-1)*(n_tor_local))   &
                          + (amat(ij,kl) + amat_k(ij,kl) + amat_n(ij,kl) + amat_kn(ij,kl)) * wst
                      
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

    if (present(get_terms)) cycle

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
        enddo

      enddo

    endif ! apply fft (or not)

  enddo ! j loop n_degrees
enddo ! i loop (n_vertex)


if (.NOT. use_fft) return

!--- THIS IS ONLY USED FOR DIAGNOSTIC PURPOSES-------------------------
!--- ELM structure is re-used to plot separate terms in vtk ----------- 
if (present(get_terms)) then

  do i_term=1, max_terms

    do j=1, n_vertex_max*n_var*n_degrees
    
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
    
    do j=1, n_vertex_max*n_var*n_degrees
    
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

endif ! get_terms

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
