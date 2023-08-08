module mod_elt_matrix_fft

  implicit none

contains

#include "corr_neg_include.f90"

subroutine element_matrix_fft(element, nodes, xpoint2, xcase2, R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint, Z_xpoint, &
                              ELM, RHS, tid, ELM_p, ELM_n, ELM_k, ELM_kn, RHS_p, RHS_k,                               &
                              eq_g, eq_s, eq_t, eq_p, eq_ss, eq_st, eq_tt, delta_g, delta_s, delta_t,                 &
                              i_tor_min, i_tor_max, aux_nodes, ELM_pnn)
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
use diffusivities, only: get_dperp, get_zk_iperp, get_zk_eperp
use equil_info, only : get_psi_n
use corr_neg
use mod_neutral_source
use mod_bootstrap_functions
use mod_sources

implicit none

type (type_element)       :: element 
type (type_node)          :: nodes(n_vertex_max)     ! fluid variables
type (type_node),optional :: aux_nodes(n_vertex_max) ! particle moments


#define DIM0 n_tor*n_vertex_max*n_degrees*n_var

integer, intent(in)            :: tid
integer, intent(in)            :: i_tor_min, i_tor_max

real*8, dimension (DIM0,DIM0)  :: ELM
real*8, dimension (DIM0)       :: RHS

integer    :: i, j, ms, mt, mp, k, l, index_ij, index_kl, index, index_k, index_m, m, ik, xcase2
integer    :: n_tor_start, n_tor_end, n_tor_local, n_tor_loop
integer    :: in, im, ij1, ij2, ij3, ij4, ij5, ij6, ij7, ij8, kl1, kl2, kl3, kl4, kl5, kl6, kl7, kl8, ij, kl
real*8     :: wst, xjac, xjac_s, xjac_t, xjac_x, xjac_y, BigR, r2, phi, delta_phi
real*8     :: current_source(n_gauss,n_gauss),particle_source(n_gauss,n_gauss),heat_source_i(n_gauss,n_gauss),heat_source_e(n_gauss,n_gauss)
real*8     :: R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint(2), Z_xpoint(2), dj_dpsi, dj_dz, source_pellet, source_volume
real*8     :: Bgrad_rho_star,     Bgrad_rho,     Bgrad_T_star,  Bgrad_Ti, Bgrad_Te, BB2
real*8     :: Bgrad_rho_star_psi, Bgrad_rho_psi, Bgrad_rho_rho, Bgrad_T_star_psi, Bgrad_Ti_psi, Bgrad_Ti_Ti, Bgrad_Te_psi, Bgrad_Te_Te, BB2_psi
real*8     :: Bgrad_rho_k_star, Bgrad_T_k_star, Bgrad_Ti_Ti_n, Bgrad_Te_Te_n, Bgrad_rho_rho_n
real*8     :: ZKi_par_T, dZKi_par_dT, ZKe_par_T, dZKe_par_dT
real*8     :: D_prof, ZKi_prof, ZKe_prof, psi_norm, theta, zeta, delta_u_x, delta_u_y, delta_ps_x, delta_ps_y
real*8     :: rhs_ij(n_var), rhs_ij_k(n_var)
real*8     :: amat(n_var,n_var), amat_k(n_var,n_var), amat_n(n_var,n_var), amat_kn(n_var,n_var)

real*8     :: v, v_x, v_y, v_s, v_t, v_p, v_ss, v_st, v_tt, v_xx, v_xy, v_yy
real*8     :: ps0, ps0_x, ps0_y, ps0_p, ps0_s, ps0_t, ps0_ss, ps0_tt, ps0_st, ps0_xx, ps0_yy, ps0_xy
real*8     :: zj0, zj0_x, zj0_y, zj0_p, zj0_s, zj0_t
real*8     :: u0, u0_x, u0_y, u0_p, u0_s, u0_t, u0_ss, u0_tt, u0_st, u0_xx, u0_xy, u0_yy
real*8     :: w0, w0_x, w0_y, w0_p, w0_s, w0_t, w0_ss, w0_st, w0_tt, w0_xx, w0_xy, w0_yy
real*8     :: r0, r0_x, r0_y, r0_p, r0_s, r0_t, r0_ss, r0_st, r0_tt, r0_xx, r0_xy, r0_yy, r0_hat, r0_x_hat, r0_y_hat, r0_corr
real*8     :: Ti0, Ti0_x, Ti0_y, Ti0_p, Ti0_s, Ti0_t, Ti0_ss, Ti0_st, Ti0_tt, Ti0_xx, Ti0_xy, Ti0_yy, Ti0_corr, dTi0_corr_dT
real*8     :: Te0, Te0_x, Te0_y, Te0_p, Te0_s, Te0_t, Te0_ss, Te0_st, Te0_tt, Te0_xx, Te0_xy, Te0_yy, Te0_corr, dTe0_corr_dT
real*8     :: psi, psi_x, psi_y, psi_p, psi_s, psi_t, psi_ss, psi_st, psi_tt, psi_xx, psi_xy, psi_yy
real*8     :: zj, zj_x, zj_y, zj_p, zj_s, zj_t, zj_ss, zj_st, zj_tt
real*8     :: u, u_x, u_y, u_p, u_s, u_t, u_ss, u_st, u_tt, u_xx, u_xy, u_yy
real*8     :: w, w_x, w_y, w_p, w_s, w_t, w_ss, w_st, w_tt, w_xx, w_xy, w_yy
real*8     :: rho, rho_x, rho_y, rho_s, rho_t, rho_p, rho_hat, rho_x_hat, rho_y_hat, rho_ss, rho_st, rho_tt, rho_xx, rho_xy, rho_yy
real*8     :: Ti, Ti_x, Ti_y, Ti_s, Ti_t, Ti_p, Ti_ss, Ti_st, Ti_tt, Ti_xx, Ti_xy, Ti_yy
real*8     :: Te, Te_x, Te_y, Te_s, Te_t, Te_p, Te_ss, Te_st, Te_tt, Te_xx, Te_xy, Te_yy
real*8	   :: zTi, zTi_x, zTi_y, zTe, zTe_x, zTe_y, zn_x, zn_y, Jb_0 , Jb
real*8     :: Vpar, Vpar_x, Vpar_y, Vpar_p, Vpar_s, Vpar_t, Vpar_ss, Vpar_st, Vpar_tt, Vpar_xx, Vpar_yy, Vpar_xy
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
real*8     :: TG_num1, TG_num2, TG_num5, TG_num6, TG_num7, TG_num8

real*8     :: Vt0,Omega_tor0_x,Omega_tor0_y,Vt0_x,Vt0_y
real*8     :: V_source(n_gauss,n_gauss), Vt_x_psi, Vt_y_psi, Omega_tor_x_psi, Omega_tor_y_psi

real*8     :: dV_dpsi_source(n_gauss,n_gauss), dV_dz_source(n_gauss,n_gauss)
real*8     :: dV_dpsi2,dV_dz2,dV_dpsi_dz,dV_dpsi3,dV_dpsi_dz2,dV_dpsi2_dz
real*8     :: eq_zne(n_gauss,n_gauss),  eq_zTi(n_gauss,n_gauss), eq_zTe(n_gauss,n_gauss)
real*8     :: dn_dpsi(n_gauss,n_gauss), dn_dz,  dn_dpsi2,  dn_dz2,  dn_dpsi_dz,  dn_dpsi3,  dn_dpsi_dz2,  dn_dpsi2_dz
real*8     :: dTi_dpsi(n_gauss,n_gauss),dTi_dz, dTi_dpsi2, dTi_dz2, dTi_dpsi_dz, dTi_dpsi3, dTi_dpsi_dz2, dTi_dpsi2_dz
real*8     :: dTe_dpsi(n_gauss,n_gauss),dTe_dz, dTe_dpsi2, dTe_dz2, dTe_dpsi_dz, dTe_dpsi3, dTe_dpsi_dz2, dTe_dpsi2_dz

logical    :: xpoint2, use_fft
real*8     :: Btheta2, epsil, Btheta2_psi
real*8, dimension(n_gauss,n_gauss)    :: amu_neo_prof, aki_neo_prof
real*8     :: t_norm

real*8     :: in_fft(1:n_plane)
complex*16 :: out_fft(1:n_plane)
integer*8  :: plan

integer    :: i_v, i_loc, j_loc

!   -Ion-electron energy transfer
real*8     :: nu_e_bg, lambda_e_bg, dTi_e, dTe_i
real*8     :: dnu_e_bg_dTi, dnu_e_bg_dTe
real*8     :: dnu_e_bg_drho, dnu_e_bg_drhon
real*8     :: ddTi_e_dTi, ddTi_e_dTe, ddTi_e_drho, ddTi_e_drhon
real*8     :: ddTe_i_dTi, ddTe_i_dTe, ddTe_i_drho, ddTe_i_drhon
real*8     :: Te_corr_eV, dTe_corr_eV_dT                      ! Electron temperature in eV
real*8     :: ne_SI                                          ! Electron density in SI unit
real*8     :: dr0_corr_dn

#define DIM1 n_plane
#define DIM2 1:n_vertex_max*n_var*n_degrees

real*8, dimension(DIM1, DIM2, DIM2) :: ELM_p
real*8, dimension(DIM1, DIM2, DIM2) :: ELM_n
real*8, dimension(DIM1, DIM2, DIM2) :: ELM_k
real*8, dimension(DIM1, DIM2, DIM2) :: ELM_kn
real*8, dimension(DIM1, DIM2)       :: RHS_p
real*8, dimension(DIM1, DIM2)       :: RHS_k
real*8, dimension(DIM1, DIM2, DIM2) :: ELM_pnn

real*8, dimension(n_gauss,n_gauss)    :: x_g, x_s, x_t, x_ss, x_st, x_tt
real*8, dimension(n_gauss,n_gauss)    :: y_g, y_s, y_t, y_ss, y_st, y_tt

real*8, dimension(n_plane,n_var,n_gauss,n_gauss) :: eq_g, eq_s, eq_t, eq_p, eq_ss, eq_st, eq_tt
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

epsil=1.d-3

! --- Taylor-Galerkin Stabilisation coefficients
TG_num1    = TGNUM(1); TG_num2    = TGNUM(2); TG_num5    = TGNUM(5); TG_num6    = TGNUM(6); TG_num7    = TGNUM(7); TG_num8    = TGNUM(8);

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
V_source        = 0.d0
dV_dpsi_source  = 0.d0
dV_dz_source    = 0.d0
eq_zne          = 0.d0
eq_zTi          = 0.d0         
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
  
    call sources(xpoint2, xcase2, y_g(ms,mt), Z_xpoint, eq_g(1,1,ms,mt),psi_axis,psi_bnd, &
                 particle_source(ms,mt),heat_source_i(ms,mt),heat_source_e(ms,mt))

    ! Source of parallel velocity
    if ( ( abs(V_0) .ge. 1.e-12 ) .or. ( num_rot ) ) then
      call velocity(xpoint2, xcase2, y_g(ms,mt), Z_xpoint, eq_g(1,1,ms,mt), psi_axis, psi_bnd, V_source(ms,mt), &
                    dV_dpsi_source(ms,mt),dV_dz_source(ms,mt),dV_dpsi2,dV_dz2,dV_dpsi_dz,dV_dpsi3,dV_dpsi_dz2, dV_dpsi2_dz)
    endif

    call density(xpoint2, xcase2, y_g(ms,mt), Z_xpoint, eq_g(1,1,ms,mt),psi_axis,psi_bnd,eq_zne(ms,mt), &
                 dn_dpsi(ms,mt),dn_dz,dn_dpsi2,dn_dz2,dn_dpsi_dz,dn_dpsi3,dn_dpsi_dz2, dn_dpsi2_dz)

    call temperature_i(xpoint2, xcase2, y_g(ms,mt), Z_xpoint, eq_g(1,1,ms,mt),psi_axis,psi_bnd, eq_zTi(ms,mt), &
                     dTi_dpsi(ms,mt),dTi_dz,dTi_dpsi2,dTi_dz2,dTi_dpsi_dz,dTi_dpsi3,dTi_dpsi_dz2, dTi_dpsi2_dz)

    call temperature_e(xpoint2, xcase2, y_g(ms,mt), Z_xpoint, eq_g(1,1,ms,mt),psi_axis,psi_bnd,eq_zTe(ms,mt), &
                     dTe_dpsi(ms,mt),dTe_dz,dTe_dpsi2,dTe_dz2,dTe_dpsi_dz,dTe_dpsi3,dTe_dpsi_dz2, dTe_dpsi2_dz)

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

!--------------------------------------------------- sum over the Gaussian integration points
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

          Vpar0    = eq_g(mp,var_Vpar,ms,mt)
          Vpar0_x  = (   y_t(ms,mt) * eq_s(mp,var_Vpar,ms,mt) - y_s(ms,mt) * eq_t(mp,var_Vpar,ms,mt) ) / xjac
          Vpar0_y  = ( - x_t(ms,mt) * eq_s(mp,var_Vpar,ms,mt) + x_s(ms,mt) * eq_t(mp,var_Vpar,ms,mt) ) / xjac
          Vpar0_p  = eq_p(mp,var_Vpar,ms,mt)
          Vpar0_s  = eq_s(mp,var_Vpar,ms,mt)
          Vpar0_t  = eq_t(mp,var_Vpar,ms,mt)
          Vpar0_ss = eq_ss(mp,var_Vpar,ms,mt)
          Vpar0_st = eq_st(mp,var_Vpar,ms,mt)
          Vpar0_tt = eq_tt(mp,var_Vpar,ms,mt)

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

          P0   = Pi0   + Pe0
          P0_s = Pi0_s + Pe0_s
          P0_t = Pi0_t + Pe0_t
          P0_p = Pi0_p + Pe0_p
          P0_x = Pi0_x + Pe0_x
          P0_y = Pi0_y + Pe0_y

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

          ! --- Temperature dependent resistivity
          if ( eta_T_dependent .and. Te0_corr <= T_max_eta) then
            eta_T     = eta   * (corr_neg_temp1(Te0)/Te_0)**(-1.5d0)
            deta_dT   = - eta   * (1.5d0)  * Te0_corr**(-2.5d0) * Te_0**(1.5d0)
            d2eta_d2T =   eta   * (3.75d0) * Te0_corr**(-3.5d0) * Te_0**(1.5d0)
          else if ( eta_T_dependent .and. Te0_corr > T_max_eta) then
            eta_T     = eta   * (T_max_eta/Te_0)**(-1.5d0)
            deta_dT   = 0.
            d2eta_d2T = 0.     
          else
            eta_T     = eta
            deta_dT   = 0.d0
            d2eta_d2T = 0.d0
          end if

          if ( eta_T_dependent .and.  xpoint2 .and. (Te0 .lt. T_min) ) then
              eta_T     = eta    * (T_min/Te_0)**(-1.5d0)
              deta_dT   = 0.d0
              d2eta_d2T = 0.d0
          end if

          ! --- Eta for ohmic heating
          if ( eta_T_dependent .and. Te0_corr <= T_max_eta_ohm) then
            eta_T_ohm     = eta_ohmic   * (Te0_corr/Te_0)**(-1.5d0)
            deta_dT_ohm   = - eta_ohmic   * (1.5d0)  * Te0_corr**(-2.5d0) * Te_0**(1.5d0)
          else if ( eta_T_dependent .and. Te0_corr > T_max_eta_ohm) then
            eta_T_ohm     = eta_ohmic   * (T_max_eta_ohm/Te_0)**(-1.5d0)
            deta_dT_ohm   = 0.    
          else
            eta_T_ohm     = eta_ohmic
            deta_dT_ohm   = 0.d0
          end if

          ! --- Temperature dependent viscosity
          if ( visco_T_dependent ) then
            visco_T     =   visco * (Te0_corr/Te_0)**(-1.5d0)
            dvisco_dT   = - visco * (1.5d0)  * Te0_corr**(-2.5d0) * Te_0**(1.5d0)
            d2visco_dT2 =   visco * (3.75d0) * Te0_corr**(-3.5d0) * Te_0**(1.5d0)
            if ( xpoint2 .and. (Te0 .lt. T_min) ) then
              visco_T     = visco  * (T_min/Te_0)**(-1.5d0)
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
            ZKi_par_T   = ZK_i_par * (Ti0_corr/Ti_0)**(+2.5d0)              ! temperature dependent parallel conductivity
            dZKi_par_dT = ZK_i_par * (2.5d0)  * Ti0_corr**(+1.5d0) * Ti_0**(-2.5d0) * dTi0_corr_dT
            if (ZKi_par_T .gt. ZK_par_max) then
              ZKi_par_T   = Zk_par_max
              dZKi_par_dT = 0.d0
            endif
            if ( xpoint2 .and. (Ti0 .lt. Ti_min_ZKpar) ) then
              ZKi_par_T   = ZK_i_par * (Ti_min_ZKpar/Ti_0)**(+2.5d0)
              dZKi_par_dT = 0.d0
            endif

            ZKe_par_T   = ZK_e_par * (Te0_corr/Te_0)**(+2.5d0)              ! temperature dependent parallel conductivity
            dZKe_par_dT = ZK_e_par * (2.5d0)  * Te0_corr**(+1.5d0) * Te_0**(-2.5d0) * dTe0_corr_dT
            if (ZKe_par_T .gt. ZK_par_max) then
              ZKe_par_T   = Zk_par_max
              dZKe_par_dT = 0.d0
            endif
            if ( xpoint2 .and. (Te0 .lt. Te_min_ZKpar) ) then
              ZKe_par_T   = ZK_e_par * (Te_min_ZKpar/Te_0)**(+2.5d0)
              dZKe_par_dT = 0.d0
            endif
          else
            ZKi_par_T   = ZK_i_par                                            ! parallel conductivity
            dZKi_par_dT = 0.d0
            ZKe_par_T   = ZK_e_par                                            ! parallel conductivity
            dZKe_par_dT = 0.d0
          endif

          ! --- Temperature dependent hyper-resistivity
          if ( eta_num_T_dependent ) then
            eta_num_T     =   eta_num   * (Te0_corr/Te_0)**(-3.d0)
            deta_num_dT   = - eta_num   * (3.d0)  * Te0_corr**(-4.d0) * Te_0**(3.d0)
            if ( xpoint2 .and. (Te0 .lt. T_min) ) then
              eta_num_T     = eta_num    * (T_min/Te_0)**(-3.d0)
              deta_num_dT   = 0.d0
            endif
          else
            eta_num_T     = eta_num
            deta_num_dT   = 0.d0
          end if
   
          ! --- Temperature dependent hyper-viscosity
          if ( visco_num_T_dependent ) then
            visco_num_T     =   visco_num   * (Te0_corr/Te_0)**(-3.d0)
            dvisco_num_dT   = - visco_num   * (3.d0)  * Te0_corr**(-4.d0) * Te_0**(3.d0)
            if ( xpoint2 .and. (Te0 .lt. T_min) ) then
              visco_num_T     = visco_num    * (T_min/Te_0)**(-3.d0)
              dvisco_num_dT   = 0.d0
            endif
          else
            visco_num_T     = visco_num
            dvisco_num_dT   = 0.d0
          end if


          ! --- Diamagnetic viscosity
          if (Wdia) then
            W_dia = + tauIC*2. /r0_corr    * (pi0_xx + pi0_x/bigR + pi0_yy) &
                    - tauIC*2. /r0_corr**2 * (r0_x*pi0_x + r0_y*pi0_y)
          else
            W_dia = 0.d0
          endif

   !--------------------------------------------------------
   ! --- Ion-electron energy transfer
   !--------------------------------------------------------

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
        
            dTe_i    = nu_e_bg * (Ti0_corr - Te0_corr) * r0_corr
            dTi_e    = -dTe_i
        
            !Calculating the density and temperature derivative for amats
            !We negelect the coulomb log's dericatives due to their smallness
        
            dnu_e_bg_dTi    = -1.5*MASS_ELECTRON*nu_e_bg*dTi0_corr_dT / (MASS_ELECTRON*Ti0_corr + MASS_PROTON*central_mass*Te0_corr)
            dnu_e_bg_dTe    = -1.5*MASS_PROTON*central_mass*nu_e_bg*dTe0_corr_dT &
                              / (MASS_ELECTRON*Ti0_corr + MASS_PROTON*central_mass*Te0_corr)
  
            dnu_e_bg_drho   = nu_e_bg * dr0_corr_dn / r0_corr
  
            ddTe_i_dTi      = dnu_e_bg_dTi * (Ti0_corr - Te0_corr) * r0_corr + nu_e_bg * dTi0_corr_dT * r0_corr
            ddTe_i_dTe      = dnu_e_bg_dTe * (Ti0_corr - Te0_corr) * r0_corr - nu_e_bg * dTe0_corr_dT * r0_corr
            ddTe_i_drho     = dnu_e_bg_drho * (Ti0_corr - Te0_corr) * r0_corr &
                              + nu_e_bg * (Ti0_corr - Te0_corr) * dr0_corr_dn
        
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


          psi_norm = get_psi_n( ps0, y_g(ms,mt))

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
            
            zTi   = eq_zTi(ms,mt)        
            zTi_x = dTi_dpsi(ms,mt) * ps0_x
            zTi_y = dTi_dpsi(ms,mt) * ps0_y
            zTe   = eq_zTe(ms,mt)  
            zTe_x = dTe_dpsi(ms,mt) * ps0_x 
            zTe_y = dTe_dpsi(ms,mt) * ps0_y
            zn_x  = dn_dpsi(ms,mt)  * ps0_x
            zn_y  = dn_dpsi(ms,mt)  * ps0_y

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

          D_prof   = get_dperp (psi_norm)
          ZKi_prof = get_zk_iperp(psi_norm)
          ZKe_prof = get_zk_eperp(psi_norm)

          ! --- Increase diffusivity if very small density/temperature
          if (xpoint2) then
            if (r0 .lt. D_prof_neg_thresh)  then
              D_prof  = D_prof_neg
            endif
            if (Te0 .lt. ZK_prof_neg_thresh) then
              ZKe_prof = ZK_prof_neg
            endif
            if (Ti0 .lt. ZK_prof_neg_thresh) then
              ZKi_prof = ZK_prof_neg
            endif
          endif

          phi       = 2.d0*PI*float(mp-1)/float(n_plane) / float(n_period)
          delta_phi = 2.d0*PI/float(n_plane) / float(n_period)

          source_pellet = 0.d0
          source_volume = 0.d0

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

            BB2              = (F0*F0 + ps0_x * ps0_x + ps0_y * ps0_y )/BigR**2
            Btheta2          = (ps0_x * ps0_x + ps0_y * ps0_y )/BigR**2

            v_ps0_x  = v_xx  * ps0_y - v_xy  * ps0_x + v_x  * ps0_xy - v_y * ps0_xx
            v_ps0_y  = v_xy  * ps0_y - v_yy  * ps0_x + v_x  * ps0_yy - v_y * ps0_xy

            !###################################################################################################
            !#  equation 1   (induction equation)                                                              #
            !###################################################################################################

            rhs_ij(1) = v * eta_T  * (zj0 - current_source(ms,mt) - Jb)/ BigR * xjac * tstep &
                      + v * (ps0_s * u0_t - ps0_t * u0_s)                            * tstep &
                      - v * F0 / BigR  * u0_p                                 * xjac * tstep &
                      + eta_num_T * (v_x * zj0_x + v_y * zj0_y)               * xjac * tstep &

                      - v * tauIC*2./(r0_corr*BB2) * F0**2/BigR**2 * (ps0_s * pe0_t - ps0_t * pe0_s) * tstep &
                      + v * tauIC*2./(r0_corr*BB2) * F0**3/BigR**3 * pe0_p    * xjac * tstep &

                      + zeta * v * delta_g(mp,1,ms,mt) / BigR                 * xjac


            !###################################################################################################
            !#  equation 2   (perpendicular momentum equation)                                                 #
            !###################################################################################################

            rhs_ij(2) =  - 0.5d0 * vv2 * (v_x * r0_y_hat - v_y * r0_x_hat)   * xjac * tstep &
                         - r0_hat * BigR**2 * w0 * (v_s * u0_t - v_t * u0_s)        * tstep &
                         + v * (ps0_s * zj0_t - ps0_t * zj0_s )                     * tstep &

                         - visco_T * BigR * (v_x * w0_x + v_y * w0_y)        * xjac * tstep &

                         - v * F0 / BigR * zj0_p                             * xjac * tstep &
                         + BigR**2 * (v_s * p0_t - v_t * p0_s)                      * tstep &

                         - visco_num_T * (v_xx + v_x/Bigr + v_yy)*(w0_xx + w0_x/Bigr + w0_yy) * xjac * tstep &

                         - TG_num2 * 0.25d0 * r0_hat * BigR**3 * (w0_x * u0_y - w0_y * u0_x) &
                                   * ( v_x * u0_y - v_y * u0_x) * xjac * tstep * tstep &

                         - v * tauIC*2. * BigR**4 * (pi0_s * w0_t - pi0_t * w0_s)       * tstep &

                         - tauIC*2. * BigR**3 * pi0_y * (v_x* u0_x + v_y * u0_y) * xjac * tstep &

                         - v * tauIC*2. * BigR**4 * (u0_xy * (pi0_xx - pi0_yy) - pi0_xy * (u0_xx - u0_yy) ) * xjac * tstep &

                         ! --- Diamagnetic viscosity
                         + dvisco_dT * bigR * W_dia * (v_x*Ti0_x + v_y*Ti0_y)  * xjac * tstep  &
                         + visco_T   * bigR * W_dia * (v_xx + v_x/bigR + v_yy) * xjac * tstep  &

                         + BigR**3 * (particle_source(ms,mt) + source_pellet) * (v_x * u0_x + v_y * u0_y) * xjac* tstep &
                     
                         - zeta * BigR * r0_hat * (v_x * delta_u_x + v_y * delta_u_y) * xjac
            
            !------------------------------------------------------------------------ NEO
            if (NEO) then
              rhs_ij(2) = rhs_ij(2)  + amu_neo_prof(ms,mt)*BB2/((Btheta2+epsil)**2)             &
                        * (ps0_x * v_x + ps0_y * v_y) *                                         &  
                          ( r0 * (ps0_x * u0_x + ps0_y * u0_y)                                  &
                            + tauIC*2. * (ps0_x * Pi0_x + ps0_y * Pi0_y)                        &
                            + aki_neo_prof(ms,mt) * tauIC*2. * r0 * (ps0_x*Ti0_x + ps0_y*Ti0_y) &
                            - r0 * Vpar0 * Btheta2)     * xjac * tstep * BigR 
            endif
            !------------------------------------------------------------------------ NEO

            !###################################################################################################
            !#  equation 3   (current definition)                                                              #
            !###################################################################################################

            rhs_ij(3) = - ( v_x * ps0_x  + v_y * ps0_y + v*zj0) / BigR * xjac

            !###################################################################################################
            !#  equation 4   (vorticity definition)                                                            #
            !###################################################################################################

            rhs_ij(4) = - ( v_x * u0_x   + v_y * u0_y  + v*w0)  * BigR * xjac 

            !###################################################################################################
            !#  equation 5   (density equation)                                                                #
            !###################################################################################################

            rhs_ij(5)  = v * BigR * (particle_source(ms,mt) + source_pellet)                      * xjac * tstep &
                       + v * BigR**2 * ( r0_s * u0_t - r0_t * u0_s)                                      * tstep &
                       + v * 2.d0 * BigR * r0 * u0_y                                              * xjac * tstep &

                       - (D_par-D_prof) * BigR / BB2 * Bgrad_rho_star * Bgrad_rho                 * xjac * tstep &
                       - D_prof * BigR  * (v_x*r0_x + v_y*r0_y                                  ) * xjac * tstep &

                       - v * F0 / BigR * Vpar0 * r0_p                                             * xjac * tstep &
                       - v * Vpar0 * (r0_s * ps0_t - r0_t * ps0_s)                                       * tstep &
                       - v * F0 / BigR * r0 * vpar0_p                                             * xjac * tstep &
                       - v * r0 * (vpar0_s * ps0_t - vpar0_t * ps0_s)                                    * tstep &

                       + v * 2.d0 * tauIC*2. * pi0_y * BigR                                       * xjac * tstep &
                       
                       + zeta * v * delta_g(mp,5,ms,mt) * BigR                                    * xjac         &

                       - D_perp_num * (v_xx + v_x/Bigr + v_yy)*(r0_xx + r0_x/Bigr + r0_yy) * BigR * xjac * tstep &

                       - TG_num5 * 0.25d0 * BigR**3 * (r0_x * u0_y - r0_y * u0_x)                                &
                                                    * ( v_x * u0_y - v_y * u0_x) * xjac * tstep * tstep          &
                       - TG_num5 * 0.25d0 / BigR * vpar0**2                                                      &
                                 * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                              &
                                 * ( v_x * ps0_y -  v_y * ps0_x                   ) * xjac * tstep * tstep


            rhs_ij_k(5) = - (D_par-D_prof) * BigR / BB2 * Bgrad_rho_k_star * Bgrad_rho            * xjac * tstep &
                          - D_prof * BigR  * (                  v_p*r0_p /BigR**2 )               * xjac * tstep &

                       - TG_num5 * 0.25d0 / BigR * vpar0**2 &
                                 * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                              &
                                 * (                            + F0 / BigR * v_p) * xjac * tstep * tstep

            !###################################################################################################
            !#  equation 6   (ion energy equation)                                                             #
            !###################################################################################################

            rhs_ij(6) =  v * BigR * heat_source_i(ms,mt)                                  * xjac * tstep &
            
                       + v * r0 * BigR**2 * ( Ti0_s * u0_t - Ti0_t * u0_s)                       * tstep &
                       + v * Ti0 * BigR**2 * ( r0_s * u0_t - r0_t * u0_s)                        * tstep &

                       + v * r0 * Ti0 * 2.d0* GAMMA * BigR * u0_y                         * xjac * tstep &

                       - v * r0 * F0 / BigR * Vpar0 * Ti0_p                               * xjac * tstep &
                       - v * Ti0 * F0 / BigR * Vpar0 * r0_p                               * xjac * tstep &

                       - v * r0 * Vpar0 * (Ti0_s * ps0_t - Ti0_t * ps0_s)                        * tstep &
                       - v * Ti0 * Vpar0 * (r0_s * ps0_t - r0_t * ps0_s)                         * tstep &

                       - v * r0 * Ti0 * GAMMA * (vpar0_s * ps0_t - vpar0_t * ps0_s)              * tstep &
                       - v * r0 * Ti0 * GAMMA * F0 / BigR * vpar0_p                       * xjac * tstep &

                       - (ZKi_par_T-ZKi_prof) * BigR / BB2 * Bgrad_T_star * Bgrad_Ti      * xjac * tstep &
                       - ZKi_prof * BigR * (v_x*Ti0_x + v_y*Ti0_y                   )     * xjac * tstep &
 
                       - ZK_perp_num  *  (v_xx + v_x/Bigr + v_yy)*(Ti0_xx + Ti0_x/Bigr + Ti0_yy) * BigR * xjac * tstep &

                       - TG_num6 * 0.25d0 * BigR**3 * Ti0 * (r0_x * u0_y - r0_y * u0_x)         &
                                          * ( v_x * u0_y - v_y * u0_x) * xjac * tstep * tstep   &
                       - TG_num6 * 0.25d0 * BigR**3 * r0 * (Ti0_x * u0_y - Ti0_y * u0_x)        &
                                          * ( v_x * u0_y - v_y * u0_x) * xjac * tstep * tstep   &

                       - TG_num6 * 0.25d0 / BigR * vpar0**2 &
                                 * Ti0 * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                        &
                                 * ( v_x * ps0_y -  v_y * ps0_x                        ) * xjac * tstep * tstep  &
                       - TG_num6 * 0.25d0 / BigR * vpar0**2 &
                                 * r0 * (Ti0_x * ps0_y - Ti0_y * ps0_x + F0 / BigR * Ti0_p)                      &
                                 * ( v_x * ps0_y -  v_y * ps0_x                        ) * xjac * tstep * tstep  &

                       + zeta * v * r0_corr  * delta_g(mp,6,ms,mt) * BigR                     * xjac &
                       + zeta * v * Ti0_corr * delta_g(mp,5,ms,mt) * BigR                     * xjac &
                       ! Energy exchange term
                       + v * BigR * dTi_e                                                * xjac * tstep                       


            rhs_ij_k(6) = - (ZKi_par_T-ZKi_prof) * BigR / BB2 * Bgrad_T_k_star * Bgrad_Ti * xjac * tstep &
                          - ZKi_prof * BigR * (                + v_p*Ti0_p /BigR**2 )     * xjac * tstep &

                         - TG_num6 * 0.25d0 / BigR * vpar0**2 &
                                 * Ti0 * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                        &
                                 * (                                  + F0 / BigR * v_p) * xjac * tstep * tstep  &
                         - TG_num6 * 0.25d0 / BigR * vpar0**2 &
                                 * r0 * (Ti0_x * ps0_y - Ti0_y * ps0_x + F0 / BigR * Ti0_p)                      &
                                 * (                                   + F0 / BigR * v_p) * xjac * tstep * tstep

            !###################################################################################################
            !#  equation 7   (parallel velocity equation)                                                      #
            !###################################################################################################

            rhs_ij(7) = - v * F0 / BigR * P0_p                                             * xjac * tstep &
                        - v * (P0_s * ps0_t - P0_t * ps0_s)                                       * tstep &

                     - v*(particle_source(ms,mt) + source_pellet) * vpar0 * BB2   * BigR * xjac * tstep &

                     - 0.5d0 * r0 * vpar0**2 * BB2 * (ps0_s * v_t - ps0_t * v_s)                * tstep &
                     - 0.5d0 * v  * vpar0**2 * BB2 * (ps0_s * r0_t - ps0_t * r0_s)              * tstep &
                     + 0.5d0 * v  * vpar0**2 * BB2 * F0 / BigR * r0_p                    * xjac * tstep &

                     - visco_par_num * (v_xx + v_x/Bigr + v_yy)*(vpar0_xx + vpar0_x/Bigr + vpar0_yy) * BigR * xjac * tstep &

                     + zeta * v * delta_g(mp,7,ms,mt) * r0_corr * F0**2 / BigR                        * xjac &
                     + zeta * v * r0_corr * vpar0 * (ps0_x * delta_ps_x + ps0_y * delta_ps_y) / BigR  * xjac &

             - TG_NUM7 * 0.25d0 * r0 * Vpar0**2 * BB2 &
                       * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac + F0 / BigR * vpar0_p) / BigR  &
                       * (-(ps0_s * v_t     - ps0_t * v_s)    /xjac                      ) * xjac * tstep * tstep &
             - TG_NUM7 * 0.25d0 * v  * Vpar0**2 * BB2 &
                       * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac + F0 / BigR * vpar0_p) / BigR  &
                       * (-(ps0_s * r0_t    - ps0_t * r0_s)   /xjac + F0 / BigR * r0_p)    * xjac * tstep * tstep
                  
            if (normalized_velocity_profile) then
              rhs_ij(7) = rhs_ij(7) - visco_par * (v_x * (vpar0_x-Vt0_x) + v_y * (vpar0_y-Vt0_y)) * BigR* xjac * tstep 
            else
              rhs_ij(7) = rhs_ij(7) - visco_par * (v_x * (vpar0_x * F0**2 / BigR**2 - 2.d0 * vpar0 * F0**2 / BigR**3  - 2.d0 * PI * F0 * Omega_tor0_x ) & 
                                                 + v_y * (vpar0_y * F0**2 / BigR**2 - 2.d0 * PI * F0 * Omega_tor0_y) ) * BigR* xjac * tstep 
            endif
 
            if (NEO) then
              rhs_ij(7) =  rhs_ij(7)  + amu_neo_prof(ms,mt)*BB2/(Btheta2+epsil)  &
                        * v * ( r0 * (ps0_x * u0_x + ps0_y * u0_y)               &
                              + tauIC*2.   * (ps0_x * Pi0_x + ps0_y * Pi0_y)     &
                              + aki_neo_prof(ms,mt) * tauIC*2. * r0 * (ps0_x * Ti0_x + ps0_y * Ti0_y) - r0 * Vpar0 * Btheta2) * xjac * tstep * BigR
            endif
 
            rhs_ij_k(7) = + 0.5d0 * r0 * vpar0**2 * BB2 * F0 / BigR * v_p                 * xjac * tstep &

               - TG_NUM7 * 0.25d0 * r0 * Vpar0**2 * BB2 &
                         * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac + F0 / BigR * vpar0_p) / BigR  &
                         * (                                          + F0 / BigR * v_p)  * xjac * tstep * tstep 

            !###################################################################################################
            !#  equation 8   (electron energy equation)                                                        #
            !###################################################################################################

            rhs_ij(8) =  v * BigR * heat_source_e(ms,mt)                                    * xjac * tstep &
            
                       + v * r0 * BigR**2  * (Te0_s * u0_t - Te0_t * u0_s)                        * tstep &
                       + v * Te0 * BigR**2 * ( r0_s * u0_t -  r0_t * u0_s)                        * tstep &

                       + v * r0 * Te0 * 2.d0* GAMMA * BigR * u0_y                         * xjac * tstep &

                       - v * r0 * F0 / BigR * Vpar0 * Te0_p                               * xjac * tstep &
                       - v * Te0 * F0 / BigR * Vpar0 * r0_p                               * xjac * tstep &

                       - v * r0 * Vpar0 * (Te0_s * ps0_t - Te0_t * ps0_s)                        * tstep &
                       - v * Te0 * Vpar0 * (r0_s * ps0_t - r0_t * ps0_s)                         * tstep &

                       - v * r0 * Te0 * GAMMA * (vpar0_s * ps0_t - vpar0_t * ps0_s)              * tstep &
                       - v * r0 * Te0 * GAMMA * F0 / BigR * vpar0_p                       * xjac * tstep &

                       - (ZKe_par_T-ZKe_prof) * BigR / BB2 * Bgrad_T_star * Bgrad_Te      * xjac * tstep &
                       - ZKe_prof * BigR * (v_x*Te0_x + v_y*Te0_y                   )     * xjac * tstep &
 
                       - ZK_perp_num  *  (v_xx + v_x/Bigr + v_yy)*(Te0_xx + Te0_x/Bigr + Te0_yy) * BigR * xjac * tstep &

                       + v * (gamma-1.d0) * eta_T_ohm * (zj0 / BigR)**2.d0         * BigR  * xjac * tstep  &

                       - TG_num8 * 0.25d0 * BigR**3 * Te0 * (r0_x * u0_y - r0_y * u0_x)         &
                                          * ( v_x * u0_y - v_y * u0_x) * xjac * tstep * tstep   &
                       - TG_num8 * 0.25d0 * BigR**3 * r0 * (Te0_x * u0_y - Te0_y * u0_x)        &
                                          * ( v_x * u0_y - v_y * u0_x) * xjac * tstep * tstep   &

                       - TG_num8 * 0.25d0 / BigR * vpar0**2 &
                                 * Te0 * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                        &
                                 * ( v_x * ps0_y -  v_y * ps0_x                        ) * xjac * tstep * tstep  &
                       - TG_num8 * 0.25d0 / BigR * vpar0**2 &
                                 * r0 * (Te0_x * ps0_y - Te0_y * ps0_x + F0 / BigR * Te0_p)                      &
                                 * ( v_x * ps0_y -  v_y * ps0_x                        ) * xjac * tstep * tstep  &

                       + zeta * v * r0_corr  * delta_g(mp,8,ms,mt) * BigR                     * xjac &
                       + zeta * v * Te0_corr * delta_g(mp,5,ms,mt) * BigR                     * xjac &
                       ! Energy exchange term
                       + v * BigR * dTe_i                                                * xjac * tstep


            rhs_ij_k(8) = - (ZKe_par_T-ZKe_prof) * BigR / BB2 * Bgrad_T_k_star * Bgrad_Te * xjac * tstep &
                          - ZKe_prof * BigR * (                + v_p*Te0_p /BigR**2 )     * xjac * tstep &

                         - TG_num8 * 0.25d0 / BigR * vpar0**2 &
                                 * Te0 * (r0_x * ps0_y - r0_y * ps0_x  + F0 / BigR * r0_p)                        &
                                 * (                                   + F0 / BigR * v_p) * xjac * tstep * tstep  &
                         - TG_num8 * 0.25d0 / BigR * vpar0**2 &
                                 * r0 * (Te0_x * ps0_y - Te0_y * ps0_x + F0 / BigR * Te0_p)                       &
                                 * (                                   + F0 / BigR * v_p) * xjac * tstep * tstep

            !###################################################################################################
            !#  RHS equations end                                                                              #
            !###################################################################################################

            if (use_fft) then
              index_ij = n_var*n_degrees*(i-1) +       n_var*(j-1) + 1
            else
              index_ij = n_tor_local*n_var*n_degrees*(i-1) + n_tor_local*n_var*(j-1) + im - n_tor_start +1 
            endif


            ! --- Fill up the rhs
            if (use_fft) then

              do ij = 1, n_var
                RHS_p(mp,index_ij+ij-1) = RHS_p(mp,index_ij+ij-1) + rhs_ij(ij)   * wst
                RHS_k(mp,index_ij+ij-1) = RHS_k(mp,index_ij+ij-1) + rhs_ij_k(ij) * wst
              enddo

            else

              do ij = 1, n_var
                RHS(index_ij+(ij-1)*(n_tor_local)) = RHS(index_ij+(ij-1)*(n_tor_local)) &
                                                   + (rhs_ij(ij) + rhs_ij_k(ij)) * wst
              enddo

            endif

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

                  u    = psi    ;  zj    = psi    ;  w    = psi    ; rho    = psi    ;  Ti    = psi    ; vpar    = psi   ; Te   = psi
                  u_x  = psi_x  ;  zj_x  = psi_x  ;  w_x  = psi_x  ; rho_x  = psi_x  ;  Ti_x  = psi_x  ; vpar_x  = psi_x ; Te_x = psi_x
                  u_y  = psi_y  ;  zj_y  = psi_y  ;  w_y  = psi_y  ; rho_y  = psi_y  ;  Ti_y  = psi_y  ; vpar_y  = psi_y ; Te_y = psi_y
                  u_p  = psi_p  ;  zj_p  = psi_p  ;  w_p  = psi_p  ; rho_p  = psi_p  ;  Ti_p  = psi_p  ; vpar_p  = psi_p ; Te_p = psi_p
                  u_s  = psi_s  ;  zj_s  = psi_s  ;  w_s  = psi_s  ; rho_s  = psi_s  ;  Ti_s  = psi_s  ; vpar_s  = psi_s ; Te_s = psi_s
                  u_t  = psi_t  ;  zj_t  = psi_t  ;  w_t  = psi_t  ; rho_t  = psi_t  ;  Ti_t  = psi_t  ; vpar_t  = psi_t ; Te_t = psi_t
                  u_ss = psi_ss ;  zj_ss = psi_ss ;  w_ss = psi_ss ; rho_ss = psi_ss ;  Ti_ss = psi_ss ; vpar_ss = psi_ss; Te_ss = psi_ss
                  u_tt = psi_tt ;  zj_tt = psi_tt ;  w_tt = psi_tt ; rho_tt = psi_tt ;  Ti_tt = psi_tt ; vpar_tt = psi_tt; Te_tt = psi_tt
                  u_st = psi_st ;  zj_st = psi_st ;  w_st = psi_st ; rho_st = psi_st ;  Ti_st = psi_st ; vpar_st = psi_st; Te_st = psi_st

                  u_xx = psi_xx ;                    w_xx = psi_xx ; rho_xx = psi_xx ;  Ti_xx = psi_xx ; vpar_xx = psi_xx ; Te_xx = psi_xx
                  u_yy = psi_yy ;                    w_yy = psi_yy ; rho_yy = psi_yy ;  Ti_yy = psi_yy ; vpar_yy = psi_yy ; Te_yy = psi_yy
                  u_xy = psi_xy ;                    w_xy = psi_xy ; rho_xy = psi_xy ;  Ti_xy = psi_xy ; vpar_xy = psi_xy ; Te_xy = psi_xy

                  rho_hat   = BigR**2 * rho
                  rho_x_hat = 2.d0 * BigR * BigR_x  * rho + BigR**2 * rho_x
                  rho_y_hat = BigR**2 * rho_y

                  Btheta2_psi  = 2.d0 * (psi_x * ps0_x + psi_y * ps0_y ) /BigR**2
                  Bgrad_rho_star_psi = ( v_x  * psi_y - v_y  * psi_x ) / BigR
                  Bgrad_rho_psi      = ( r0_x * psi_y - r0_y * psi_x ) / BigR
                  Bgrad_rho_rho      = ( rho_x * ps0_y - rho_y * ps0_x ) / BigR
                  Bgrad_rho_rho_n    = ( F0 / BigR * rho_p ) / BigR
                  BB2_psi            = 2.d0 * (psi_x * ps0_x + psi_y * ps0_y ) /BigR**2

                  Pi0_x_rho  = rho_x*Ti0 + rho*Ti0_x
                  Pi0_xx_rho = rho_xx*Ti0 + 2.0*rho_x*Ti0_x + rho*Ti0_xx
                  Pi0_y_rho  = rho_y*Ti0 + rho*Ti0_y
                  Pi0_yy_rho = rho_yy*Ti0 + 2.0*rho_y*Ti0_y + rho*Ti0_yy
                  Pi0_xy_rho = rho_xy*Ti0 + rho*Ti0_xy + rho_x*Ti0_y + rho_y*Ti0_x
                 
                  Pi0_x_Ti  = r0_x*Ti  + r0*Ti_x
                  Pi0_xx_Ti = r0_xx*Ti + 2.0*r0_x*Ti_x + r0*Ti_xx
                  Pi0_y_Ti  = r0_y*Ti  + r0*Ti_y
                  Pi0_yy_Ti = r0_yy*Ti + 2.0*r0_y*Ti_y + r0*Ti_yy
                  Pi0_xy_Ti = r0_xy*Ti + r0*Ti_xy + r0_x*Ti_y + r0_y*Ti_x
  
                  if (Wdia) then
                    W_dia_rho = - tauIC*2. *     rho/r0_corr**2 * (pi0_xx     + pi0_x    /bigR + pi0_yy    ) &
                                + tauIC*2.          /r0_corr    * (pi0_xx_rho + pi0_x_rho/bigR + pi0_yy_rho) &
                                + tauIC*2. * 2.0*rho/r0_corr**3 * (r0_x *pi0_x     + r0_y *pi0_y    )    &
                                - tauIC*2.          /r0_corr**2 * (rho_x*pi0_x     + rho_y*pi0_y    )    &
                                - tauIC*2.          /r0_corr**2 * (r0_x *pi0_x_rho + r0_y *pi0_y_rho)
                    W_dia_Ti  = + tauIC*2.          /r0_corr    * (pi0_xx_Ti + pi0_x_Ti/bigR + pi0_yy_Ti) &
                                - tauIC*2.          /r0_corr**2 * (r0_x  *pi0_x_Ti + r0_y  *pi0_y_Ti)
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
                  !#  equation 1   (induction equation)                                                              #
                  !###################################################################################################

                  amat(1,1) = v * psi / BigR * xjac * (1.d0 + zeta)                                                                 &
                            - v * (psi_s * u0_t - psi_t * u0_s)                                                     * theta * tstep &

                     + v * tauIC*2./(r0_corr*BB2) * F0**2/BigR**2 * (psi_s * pe0_t - psi_t * pe0_s)                 * theta * tstep &
                     - v * tauIC*2./(r0_corr*BB2**2) * BB2_psi * F0**2/BigR**2 * (ps0_x*pe0_y - ps0_y*pe0_x) * xjac * theta * tstep &
                     + v * tauIC*2./(r0_corr*BB2**2) * BB2_psi * F0**3/BigR**3 * pe0_p                       * xjac * theta * tstep

                  amat(1,2) = -  v * (ps0_s * u_t - ps0_t * u_s)                                                    * theta * tstep
                                                                                                          
                  amat_n(1,2) = + F0 / BigR * v * u_p * xjac                                                        * theta * tstep
                                                                                                          
                  amat(1,3) = - eta_num_T * (v_x * zj_x + v_y * zj_y)                                        * xjac * theta * tstep &
                              - eta_T * v * zj / BigR                                                        * xjac * theta * tstep

                  amat(1,5) = + v * tauIC*2./(r0_corr*BB2) * F0**2/BigR**2 * Te0 * (ps0_s * rho_t - ps0_t * rho_s)  * theta * tstep &
                              + v * tauIC*2./(r0_corr*BB2) * F0**2/BigR**2 * rho * (ps0_s * Te0_t - ps0_t * Te0_s)  * theta * tstep &
                              - v * tauIC*2./(r0_corr*BB2) * F0**3/BigR**3 * rho * Te0_p                     * xjac * theta * tstep &

                              - v * tauIC*2. * rho /(r0_corr**2 * BB2) * F0**2/BigR**2 * (ps0_s * pe0_t - ps0_t * pe0_s) * theta * tstep &
                              + v * tauIC*2. * rho /(r0_corr**2 * BB2) * F0**3/BigR**3 * pe0_p                    * xjac * theta * tstep 

                  amat_n(1,5) = - v * tauIC*2./(r0_corr*BB2) * F0**3/BigR**3 * Te0  * rho_p                  * xjac * theta * tstep 
                                                                                                             
                  amat(1,8) = - deta_dT * v * Te * (zj0 - current_source(ms,mt) - Jb)/ BigR                  * xjac * theta * tstep &
                              - deta_num_dT * Te * (v_x * zj0_x + v_y * zj0_y)                               * xjac * theta * tstep &
                            + v * tauIC*2./(r0_corr*BB2) * F0**2/BigR**2 * r0 * (ps0_s * Te_t  - ps0_t * Te_s)      * theta * tstep &
                            + v * tauIC*2./(r0_corr*BB2) * F0**2/BigR**2 * Te * (ps0_s * r0_t - ps0_t * r0_s)       * theta * tstep &
                            - v * tauIC*2./(r0_corr*BB2) * F0**3/BigR**3 * Te * r0_p                         * xjac * theta * tstep 

                  amat_n(1,8) = - v * tauIC*2./(r0_corr*BB2) * F0**3/BigR**3 * r0 * Te_p                     * xjac * theta * tstep 

                  !###################################################################################################
                  !#  equation 2   (perpendicular momentum equation)                                                 #
                  !###################################################################################################

                  amat(2,1) = - v * (psi_s * zj0_t - psi_t * zj0_s)                          * theta * tstep

                  ! ------------------------------------------------------ NEO
                  if (NEO) then
                    amat(2,1) = amat(2,1) &
                         - amu_neo_prof(ms,mt)*BB2/((Btheta2+epsil)**2)*(psi_x*v_x+psi_y*v_y)  &
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

                  amat(2,2) = - BigR**3 * r0_corr * (v_x * u_x + v_y * u_y) * xjac * (1.d0 + zeta)                                &
                            + r0_hat * BigR**2 * w0 * (v_s * u_t  - v_t  * u_s)                                   * theta * tstep &
                            + BigR**2 * (u_x * u0_x + u_y * u0_y) * (v_x * r0_y_hat - v_y * r0_x_hat)      * xjac * theta * tstep &

                            + tauIC*2. * BigR**3 * pi0_y * (v_x* u_x + v_y * u_y)                          * xjac * theta * tstep &

                            + v * tauIC*2. * BigR**4 * (u_xy * (pi0_xx - pi0_yy) - pi0_xy * (u_xx - u_yy)) * xjac * theta * tstep &

                            - BigR**3 * (particle_source(ms,mt)+source_pellet) * (v_x * u_x + v_y * u_y)   * xjac * theta * tstep &

                            + TG_num2 * 0.25d0 * r0_hat * BigR**3 * (w0_x * u_y - w0_y * u_x)                                     &
                                      * ( v_x * u0_y - v_y * u0_x)                                 * xjac * theta * tstep * tstep &

                            + TG_num2 * 0.25d0 * r0_hat * BigR**3 * (w0_x * u0_y - w0_y * u0_x)                                   &
                                      * ( v_x * u_y - v_y * u_x)                                   * xjac * theta * tstep * tstep

                  if ( NEO ) then
                    amat(2,2) = amat(2,2) &
                              - amu_neo_prof(ms,mt)*BB2/((Btheta2+epsil)**2) * (ps0_x*v_x + ps0_y*v_y) * r0 *(ps0_x*u_x + ps0_y*u_y) &
                                                                                                    * BigR * xjac * theta * tstep 
                  endif
                  !---------------------------------------- NEO

                  amat(2,3)   = - v * (ps0_s * zj_t  - ps0_t * zj_s)              * theta * tstep

                  amat_n(2,3) = + F0 / BigR * v * zj_p  * xjac                    * theta * tstep

                  amat(2,4) = r0_hat * BigR**2 * w  * ( v_s * u0_t - v_t * u0_s)  * theta * tstep &
                            + BigR * ( v_x * w_x + v_y * w_y) * visco_T  * xjac   * theta * tstep &
                            + v * tauIC*2. * BigR**4 * (pi0_s * w_t - pi0_t * w_s)* theta * tstep &

                            + visco_num_T * (v_xx + v_x/BigR + v_yy)*(w_xx + w_x/BigR + w_yy) * xjac * theta * tstep &

                            + TG_num2 * 0.25d0 * r0_hat * BigR**3 * (w_x * u0_y - w_y * u0_x)     &
                                    * ( v_x * u0_y - v_y * u0_x) * xjac * theta * tstep * tstep

                  amat(2,5) = + 0.5d0 * vv2 * (v_x * rho_y_hat - v_y * rho_x_hat)            * xjac * theta * tstep &
                            + rho_hat * BigR**2 * w0 * (v_s * u0_t   - v_t * u0_s)                  * theta * tstep &
                            - BigR**2 * (v_s * rho_t * (Ti0+Te0)     - v_t * rho_s * (Ti0+Te0))     * theta * tstep &
                            - BigR**2 * (v_s * rho   * (Ti0_t+Te0_t) - v_t * rho   * (Ti0_s+Te0_s)) * theta * tstep &

                            + v * tauIC*2. * BigR**4 * Ti0  * (rho_s * w0_t - rho_t * w0_s)         * theta * tstep &
                            + v * tauIC*2. * BigR**4 * rho * (Ti0_s  * w0_t - Ti0_t  * w0_s)        * theta * tstep &
                            + tauIC*2. * BigR**3 * (Ti0_y * rho + Ti0 * rho_y) * (v_x* u0_x + v_y * u0_y) * xjac * theta * tstep      &
                            + v * tauIC*2. * BigR**4 * ( (u0_xy * (rho_xx*Ti0 + 2.d0*rho_x*Ti0_x + rho*Ti0_xx                          &
                                  -  rho_yy*Ti0 - 2.d0*rho_y*Ti0_y - rho*Ti0_yy))                           &
                                  - (rho_xy*Ti0 + rho_x*Ti0_y + rho_y*Ti0_x + rho*Ti0_xy) * (u0_xx - u0_yy)  )   &
                                                      * xjac * theta * tstep                             &
                            ! --- Diamagnetic viscosity
                            - dvisco_dT * bigR * W_dia_rho * (v_x*Ti0_x + v_y*Ti0_y)  * xjac * theta * tstep  &
                            - visco_T   * bigR * W_dia_rho * (v_xx + v_x/bigR + v_yy) * xjac * theta * tstep  &
                                
                            + TG_num2 * 0.25d0 * rho_hat * BigR**3 * (w0_x * u0_y - w0_y * u0_x)         &
                                      * ( v_x * u0_y - v_y * u0_x) * xjac * theta * tstep * tstep        
                           
                  if ( NEO ) then
                    amat(2,5) = amat(2,5) &
                         - amu_neo_prof(ms,mt)*BB2/((Btheta2+epsil)**2) * (ps0_x*v_x + ps0_y*v_y) * (rho*(ps0_x*u0_x+ps0_y*u0_y) &
                         + tauIC*2. * (ps0_x*(rho_x*Ti0 + rho*Ti0_x) + ps0_y*(rho_y*Ti0 + rho*Ti0_y))                                   &
                         + aki_neo_prof(ms,mt) * tauIC*2. * rho * (ps0_x*Ti0_x + ps0_y*Ti0_y) - rho * Vpar0 * Btheta2) * BigR * xjac * tstep * theta
                  endif

                  amat(2,6) = - BigR**2 * (v_s * r0_t * Ti   - v_t * r0_s * Ti)           * theta * tstep  &
                              - BigR**2 * (v_s * r0   * Ti_t - v_t * r0   * Ti_s)         * theta * tstep  &

                            + v * tauIC*2. * BigR**4 * r0 * (Ti_s * w0_t - Ti_t * w0_s)                      * theta * tstep &
                            + v * tauIC*2. * BigR**4 * Ti  * (r0_s * w0_t - r0_t * w0_s)                     * theta * tstep &
                            + tauIC*2. * BigR**3 * (r0_y * Ti + r0 * Ti_y) * (v_x* u0_x + v_y * u0_y) * xjac * theta * tstep &
                            + v * tauIC*2. * BigR**4 * ( (u0_xy * (Ti_xx*r0 + 2.d0*Ti_x*r0_x + Ti*r0_xx                      &
                                                                 - Ti_yy*r0 - 2.d0*Ti_y*r0_y - Ti*r0_yy))                    &
                                                       - (Ti_xy * r0 + Ti_x*r0_y + Ti_y*r0_x + Ti*r0_xy) * (u0_xx - u0_yy))  &
                                                                                                      * xjac * theta * tstep &

                            - dvisco_dT     * bigR * W_dia_Ti * (v_x*Ti0_x + v_y*Ti0_y)  * xjac * theta * tstep  &
                            - visco_T       * bigR * W_dia_Ti * (v_xx + v_x/bigR + v_yy) * xjac * theta * tstep  &
                            - dvisco_dT     * bigR * W_dia   * (v_x*Ti_x  + v_y*Ti_y )   * xjac * theta * tstep
                            
                  if ( NEO ) then
                    amat(2,6) = amat(2,6) - amu_neo_prof(ms,mt)*BB2/((Btheta2+epsil)**2)*(ps0_x*v_x + ps0_y*v_y) &
                              * (tauIC*2.*(ps0_x*(r0_x*Ti+r0*Ti_x) + ps0_y*(r0_y*Ti+r0*Ti_y))                       &
                              + aki_neo_prof(ms,mt) *tauIC*2. * r0 *(ps0_x*Ti_x + ps0_y*Ti_y)) * BigR * xjac * theta * tstep 
                    
                    amat(2,7) = amu_neo_prof(ms,mt)*BB2 * Btheta2 /((Btheta2+epsil)**2)*r0*vpar*(ps0_x*v_x+ps0_y*v_y) &
                              * BigR * xjac * tstep * theta 
                  endif


                  amat(2,8) = - BigR**2 * (v_s * r0_t * Te   - v_t * r0_s * Te)           * theta * tstep  &
                              - BigR**2 * (v_s * r0   * Te_t - v_t * r0   * Te_s)         * theta * tstep  &

                              + dvisco_dT * Te * ( v_x * w0_x + v_y * w0_y ) * BigR * xjac * theta * tstep  &

                              - d2visco_dT2*Te * bigR * W_dia   * (v_x*Ti0_x + v_y*Ti0_y)  * xjac * theta * tstep  &
                              - dvisco_dT*Te   * bigR * W_dia   * (v_xx + v_x/bigR + v_yy) * xjac * theta * tstep  

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

                          - v * 2.d0 * tauIC*2. * (rho_y * Ti0 + rho*Ti0_y) * BigR                      * xjac * theta * tstep &

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

                  amat(5,6) = - v * 2.d0 * tauIC*2. * (Ti_y * r0 + Ti*r0_y) * BigR                        * xjac * theta * tstep

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

                  !###################################################################################################
                  !#  equation 6   (ion energy equation)                                                             #
                  !###################################################################################################
                  Bgrad_T_star_psi = ( v_x   * psi_y - v_y   * psi_x  ) / BigR
                  Bgrad_Ti_psi     = ( Ti0_x * psi_y - Ti0_y * psi_x )  / BigR
                  Bgrad_Ti_Ti      = ( Ti_x  * ps0_y - Ti_y  * ps0_x )  / BigR
                  Bgrad_Ti_Ti_n    = ( F0 / BigR * Ti_p) / BigR
                  
                  amat(6,1) = - (ZKi_par_T-ZKi_prof) * BigR * BB2_psi / BB2**2 * Bgrad_T_star     * Bgrad_Ti     * xjac * theta * tstep &
                              + (ZKi_par_T-ZKi_prof) * BigR / BB2              * Bgrad_T_star_psi * Bgrad_Ti     * xjac * theta * tstep &
                              + (ZKi_par_T-ZKi_prof) * BigR / BB2              * Bgrad_T_star     * Bgrad_Ti_psi * xjac * theta * tstep &

                            + v * r0  * Vpar0 * (Ti0_s * psi_t - Ti0_t * psi_s)                                   * theta * tstep &
                            + v * Ti0 * Vpar0 * (r0_s * psi_t - r0_t * psi_s)                                     * theta * tstep &
                            + v * r0  * GAMMA * Ti0 * (vpar0_s * psi_t - vpar0_t * psi_s)                         * theta * tstep &

                         + TG_num6 * 0.25d0 / BigR * vpar0**2                                                         &
                                   * Ti0 * (r0_x * psi_y - r0_y * psi_x)                                              &
                                   * ( v_x * ps0_y -  v_y * ps0_x                  ) * xjac * theta * tstep * tstep   &
                         + TG_num6 * 0.25d0 / BigR * vpar0**2                                                         &
                                   * r0 * (Ti0_x * psi_y - Ti0_y * psi_x)                                             &
                                   * ( v_x * ps0_y -  v_y * ps0_x                  ) * xjac * theta * tstep * tstep   & 
                         + TG_num6 * 0.25d0 / BigR * vpar0**2                                                         &
                                   * Ti0 * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                           &
                                   * ( v_x * psi_y -  v_y * psi_x ) * xjac * theta * tstep * tstep                    &
                         + TG_num6 * 0.25d0 / BigR * vpar0**2                                                         &
                                   * r0 * (Ti0_x * ps0_y - Ti0_y * ps0_x + F0 / BigR * Ti0_p)                         &
                                   * ( v_x * psi_y -  v_y * psi_x ) * xjac * theta * tstep * tstep

                  amat_k(6,1) = - (ZKi_par_T-ZKi_prof) * BigR * BB2_psi / BB2**2 * Bgrad_T_k_star * Bgrad_Ti     * xjac * theta * tstep &
                                + (ZKi_par_T-ZKi_prof) * BigR / BB2              * Bgrad_T_k_star * Bgrad_Ti_psi * xjac * theta * tstep &
  
                        + TG_num6 * 0.25d0 / BigR * vpar0**2                                                       &
                                  * Ti0 * (r0_x * psi_y - r0_y * psi_x)                                            &
                                  * (                            + F0 / BigR * v_p) * xjac * theta * tstep * tstep &
                        + TG_num6 * 0.25d0 / BigR * vpar0**2                                                       &
                                  * r0 * (Ti0_x * psi_y - Ti0_y * psi_x)                                           &
                                  * (                            + F0 / BigR * v_p) * xjac * theta * tstep * tstep


                  amat(6,2) = - v * r0 * BigR**2  * ( Ti0_x * u_y - Ti0_y * u_x)         * xjac * theta * tstep &
                              - v * Ti0 * BigR**2 * ( r0_x * u_y - r0_y * u_x)           * xjac * theta * tstep &
                              - v * r0 * 2.d0* GAMMA * BigR * Ti0 * u_y                  * xjac * theta * tstep &

                         + TG_num6 * 0.25d0 * BigR**2 * Ti0* (r0_x * u_y - r0_y * u_x)                &
                                            * ( v_x * u0_y - v_y * u0_x) * xjac * theta*tstep*tstep  &
                         + TG_num6 * 0.25d0 * BigR**2 * r0* (Ti0_x * u_y - Ti0_y * u_x)                &
                                            * ( v_x * u0_y - v_y * u0_x) * xjac * theta*tstep*tstep  &
                         + TG_num6 * 0.25d0 * BigR**2 * Ti0* (r0_x * u0_y - r0_y * u0_x)              &
                                            * ( v_x * u_y - v_y * u_x) * xjac * theta*tstep*tstep    &
                         + TG_num6 * 0.25d0 * BigR**2 * r0* (Ti0_x * u0_y - Ti0_y * u0_x)              &
                                            * ( v_x * u_y - v_y * u_x) * xjac * theta*tstep*tstep 


                  amat(6,5) = v * rho * Ti0_corr   * BigR * xjac * (1.d0 + zeta)     &
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

                         + TG_num6 * 0.25d0 * BigR**2 * Ti0* (rho_x * u0_y - rho_y * u0_x)         &
                                   * ( v_x * u0_y - v_y * u0_x) * xjac * theta*tstep*tstep         &
                         + TG_num6 * 0.25d0 * BigR**2 * rho * (Ti0_x * u0_y - Ti0_y * u0_x)        &
                                   * ( v_x * u0_y - v_y * u0_x) * xjac* theta*tstep*tstep          &
                         + TG_num6 * 0.25d0 / BigR * vpar0**2                                      &
                                   * Ti0 * (rho_x * ps0_y - rho_y * ps0_x )                        &
                                   * ( v_x * ps0_y -  v_y * ps0_x ) * xjac * theta * tstep * tstep &
                         + TG_num6 * 0.25d0 / BigR * vpar0**2                                      &
                                   * rho * (Ti0_x * ps0_y - Ti0_y * ps0_x + F0 / BigR * Ti0_p)     &
                                   * ( v_x * ps0_y -  v_y * ps0_x ) * xjac * theta * tstep * tstep

                  amat_n(6,5) = + v * Ti0  * F0 / BigR * Vpar0 * rho_p   * xjac * theta * tstep    &

                         + TG_num6 * 0.25d0 / BigR * vpar0**2                                                       &
                                   * Ti0 * (                              + F0 / BigR * rho_p)                      &
                                   * ( v_x * ps0_y -  v_y * ps0_x                  ) * xjac * theta * tstep * tstep

                  amat_k(6,5) = + TG_num6 * 0.25d0 / BigR * vpar0**2                                                &
                                   * Ti0 * (rho_x * ps0_y - rho_y * ps0_x                    )                      &
                                   * (                            + F0 / BigR * v_p) * xjac * theta * tstep * tstep &
                              + TG_num6 * 0.25d0 / BigR * vpar0**2                                                  &
                                   * rho * (Ti0_x * ps0_y - Ti0_y * ps0_x + F0 / BigR * Ti0_p)                      &
                                   * (                            + F0 / BigR * v_p) * xjac * theta * tstep * tstep

                  amat_kn(6,5) = + TG_num6 * 0.25d0 / BigR * vpar0**2                 &
                                   * Ti0 * (+ F0 / BigR * rho_p)                      &
                                   * (      + F0 / BigR * v_p) * xjac * theta * tstep * tstep


                  amat(6,6) = v * r0_corr * Ti  * BigR * xjac * (1.d0 + zeta)     &
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
  
                            + ZK_perp_num * (v_xx + v_x/BigR + v_yy)*(Ti_xx + Ti_x/BigR + Ti_yy) * BigR * xjac * theta * tstep &

!!!!                            -v * Te * (gamma-1.d0) * deta_dT_ohm * (zj0 / BigR)**2.d0 * BigR * xjac * theta * tstep &

                            + TG_num6 * 0.25d0 * BigR**2 * Ti* (r0_x * u0_y - r0_y * u0_x)         &
                                      * ( v_x * u0_y - v_y * u0_x) * xjac * theta * tstep * tstep  &
                            + TG_num6 * 0.25d0 * BigR**2 * r0* (Ti_x * u0_y - Ti_y * u0_x)         &
                                      * ( v_x * u0_y - v_y * u0_x) * xjac * theta * tstep * tstep  &
                            + TG_num6 * 0.25d0 / BigR * vpar0**2                                                       &
                                      * Ti * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                          &
                                      * ( v_x * ps0_y -  v_y * ps0_x                  ) * xjac * theta * tstep * tstep &
                            + TG_num6 * 0.25d0 / BigR * vpar0**2                                                       &
                                      * r0 * (Ti_x * ps0_y - Ti_y * ps0_x             )                                &
                                      * ( v_x * ps0_y -  v_y * ps0_x                  ) * xjac * theta * tstep * tstep 

                  amat_k(6,6) = + (ZKi_par_T-ZKi_prof) * BigR / BB2 * Bgrad_T_k_star * Bgrad_Ti_Ti * xjac * theta * tstep  &
                                + dZKi_par_dT * Ti     * BigR / BB2 * Bgrad_T_k_star * Bgrad_Ti    * xjac * theta * tstep  &

                              + TG_num6 * 0.25d0 / BigR * vpar0**2                                                      &
                                  * Ti * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                               &
                                  * (                                 + F0 / BigR * v_p) * xjac * theta * tstep * tstep &
                              + TG_num6 * 0.25d0 / BigR * vpar0**2                                                      &
                                  * r0 * (Ti_x * ps0_y - Ti_y * ps0_x                  )                                &
                                  * (                                + F0 / BigR * v_p) * xjac * theta * tstep * tstep

                  amat_n(6,6) = + (ZKi_par_T-ZKi_prof) * BigR / BB2 * Bgrad_T_star   * Bgrad_Ti_Ti_n  * xjac * theta * tstep &

                              + v * r0 * F0 / BigR * Vpar0 * Ti_p                               * xjac * theta * tstep &
  
                              + TG_num6 * 0.25d0 / BigR * vpar0**2                       & 
                                * r0 * ( + F0 / BigR * Ti_p) * ( v_x * ps0_y -  v_y * ps0_x ) * xjac * theta * tstep * tstep

                  amat_kn(6,6) = + (ZKi_par_T-ZKi_prof) * BigR / BB2 * Bgrad_T_k_star * Bgrad_Ti_Ti_n * xjac * theta * tstep &
                                 + ZKi_prof * BigR   * (v_p*Ti_p /BigR**2 )                           * xjac * theta * tstep &

                              + TG_num6 * 0.25d0 / BigR * vpar0**2 &
                                * r0 * ( + F0 / BigR * Ti_p) * ( + F0 / BigR * v_p)          * xjac * theta * tstep * tstep

                  amat(6,7) = + v * r0 * F0 / BigR * Vpar * Ti0_p                            * xjac * theta * tstep &
                            + v * Ti0 * F0 / BigR * Vpar * r0_p                              * xjac * theta * tstep &

                            + v * r0  * Vpar * (Ti0_s * ps0_t - Ti0_t * ps0_s)                      * theta * tstep &
                            + v * Ti0 * Vpar * (r0_s * ps0_t - r0_t * ps0_s)                        * theta * tstep & 

                            + v * r0 * GAMMA * Ti0 * (vpar_s * ps0_t - vpar_t * ps0_s)              * theta * tstep &
  
                            + TG_num6 * 0.25d0 / BigR * 2.d0 * vpar0*vpar &
                                  * Ti0 * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                               &
                                  * ( v_x * ps0_y -  v_y * ps0_x                        ) * xjac * theta * tstep * tstep &
                            + TG_num6 * 0.25d0 / BigR * 2.d0 * vpar0*vpar &
                                  * r0 * (Ti0_x * ps0_y - Ti0_y * ps0_x + F0 / BigR * Ti0_p)                             &
                                  * ( v_x * ps0_y -  v_y * ps0_x                        ) * xjac * theta * tstep * tstep 

                  amat_k(6,7) =  &
                        + TG_num6 * 0.25d0 / BigR * 2.d0 * vpar0*vpar &
                                  * Ti0 * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                          &
                                  * (                            + F0 / BigR * v_p) * xjac * theta * tstep * tstep &
                        + TG_num6 * 0.25d0 / BigR * 2.d0 * vpar0*vpar &
                                  * r0 * (Ti0_x * ps0_y - Ti0_y * ps0_x + F0 / BigR * Ti0_p)                          &
                                  * (                            + F0 / BigR * v_p) * xjac * theta * tstep * tstep 

                  amat_n(6,7) = + v * r0 * GAMMA * Ti0 * F0 / BigR * vpar_p          * xjac * theta * tstep

                  amat(6,8) = - v * BigR * ddTi_e_dTe * Te                           * xjac * theta * tstep

                  !###################################################################################################
                  !#  equation 7   (parallel velocity equation)                                                      #
                  !###################################################################################################

                  amat(7,1) = v * r0_corr * vpar0 / BigR * (ps0_x * psi_x + ps0_y * psi_y) * xjac * (1.d0 + zeta) &

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
                    amat(7,1) = amat(7,1)  - visco_par * (v_x * Vt_x_psi   + v_y * Vt_y_psi) * BigR * xjac * theta * tstep      
                  else
                    amat(7,1) = amat(7,1)  - visco_par * 2.d0 * PI * F0 * (v_x * Omega_tor_x_psi + v_y * Omega_tor_y_psi) * BigR * xjac * theta * tstep 
                  endif

                  amat_k(7,1) = - 0.5d0 * r0 * vpar0**2 * BB2_psi * F0 / BigR * v_p    * xjac * theta * tstep 

                  !---------------------------------------- NEO
                  if ( NEO ) then
                    amat(7,1) = amat(7,1) &
                         - v * amu_neo_prof(ms,mt)*BB2/(Btheta2+epsil)*(r0 * (psi_x*u0_x + psi_y*u0_y) + tauIC*2.*(psi_x*Pi0_x + psi_y*Pi0_y) &
                         + aki_neo_prof(ms,mt) * tauIC*2. * r0 * (psi_x*Ti0_x + psi_y*Ti0_y)) * BigR * xjac * theta * tstep                   &
                         - v * amu_neo_prof(ms,mt) * (-Btheta2_psi)*BB2/((Btheta2+epsil)**2) * (r0*(ps0_x*u0_x+ps0_y*u0_y)                         &
                                                      + tauIC*2.*(ps0_x*Pi0_x+ps0_y*Pi0_y)                                                    &                     
                         + aki_neo_prof(ms,mt) * tauIC*2. * r0 * (ps0_x*Ti0_x + ps0_y*Ti0_y)) * BigR * xjac * theta * tstep

                    amat(7,2) =  amat(7,2) &
                         -v * amu_neo_prof(ms,mt)*BB2/(Btheta2+epsil) * r0 * (ps0_x*u_x + ps0_y*u_y) * BigR * xjac * theta * tstep 
                  endif
                  !---------------------------------------- NEO

                  amat(7,5) = + v * (rho_s * (Ti0+Te0)     * ps0_t - rho_t * (Ti0+Te0)     * ps0_s) * theta * tstep &
                              + v * (rho   * (Ti0_s+Te0_s) * ps0_t - rho   * (Ti0_t+Te0_t) * ps0_s) * theta * tstep &
                              + v * F0 / BigR * rho * (Ti0_p+Te0_p)                          * xjac * theta * tstep &

                            + 0.5d0 * rho * vpar0**2 * BB2 * (ps0_s * v_t   - ps0_t * v_s)    * theta * tstep &
                            + 0.5d0 * v   * vpar0**2 * BB2 * (ps0_s * rho_t - ps0_t * rho_s)  * theta * tstep &

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

                  amat_n(7,5) = + v * F0 / BigR * rho_p * (Ti0+Te0)                    * xjac * theta * tstep &
                              - 0.5d0 * v   * vpar0**2 * BB2 * F0 / BigR * rho_p       * xjac * theta * tstep &

                              + TG_NUM7 * 0.25d0 * v * Vpar0**2 * BB2 &
                                        * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac + F0 / BigR * vpar0_p) / BigR  &
                                        * (                                          + F0 / BigR * rho_p)* xjac * theta * tstep*tstep
                  
                  amat(7,6) = + v * (Ti_s * r0   * ps0_t - Ti_t * r0   * ps0_s)               * theta * tstep &
                              + v * (Ti   * r0_s * ps0_t - Ti   * r0_t * ps0_s)               * theta * tstep &
                              + v * F0 / BigR * Ti * r0_p                              * xjac * theta * tstep     
                        
                  amat_n(7,6) = + v * F0 / BigR * Ti_p * R0                            * xjac * theta * tstep

                  amat(7,7) = v * Vpar * r0_corr * F0**2 / BigR * xjac * (1.d0 + zeta) &

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
                          * (rho * (ps0_x*u0_x + ps0_y*u0_y) + tauIC*2.*(ps0_x * (rho_x*Ti0 + rho*Ti0_x) + ps0_y*(rho_y*Ti0 + rho*Ti0_y)) &
                          + aki_neo_prof(ms,mt) * tauIC*2. * rho*(ps0_x*Ti0_x + ps0_y*Ti0_y) - rho * Vpar0 * Btheta2) * BigR * xjac * tstep * theta
                    
                    amat(7,6) = amat(7,6) -v*amu_neo_prof(ms,mt)*BB2/(Btheta2+epsil)                  &
                              * (tauIC*2. * (ps0_x * (r0_x*Ti + r0*Ti_x) + ps0_y*(r0_y*Ti + r0*Ti_y)) &
                              + aki_neo_prof(ms,mt) * tauIC*2. * r0 * (ps0_x*Ti_x + ps0_y*Ti_y)) * BigR * xjac * tstep * theta
                  
                    amat(7,7) = amat(7,7) + v * amu_neo_prof(ms,mt) * BB2 * Btheta2/(Btheta2+epsil) * r0 * vpar * BigR * xjac * tstep * theta 
                  endif


                  amat(7,8)   = + v * (Te_s * r0   * ps0_t - Te_t * r0   * ps0_s)                   * theta * tstep  &
                                + v * (Te   * r0_s * ps0_t - Te   * r0_t * ps0_s)                   * theta * tstep  &
                                + v * F0 / BigR * Te * r0_p                                  * xjac * theta * tstep 
  
                  amat_n(7,8) = + v * F0 / BigR * Te_p * r0                                  * xjac * theta * tstep
                 

                  !###################################################################################################
                  !#  equation 8   (electron energy equation)                                                        #
                  !###################################################################################################
                  Bgrad_T_star_psi = ( v_x   * psi_y - v_y   * psi_x  ) / BigR
                  Bgrad_Te_psi     = ( Te0_x * psi_y - Te0_y * psi_x )  / BigR
                  Bgrad_Te_Te      = ( Te_x  * ps0_y - Te_y  * ps0_x )  / BigR
                  Bgrad_Te_Te_n    = ( F0 / BigR * Te_p) / BigR
                  
                  amat(8,1) = - (ZKe_par_T-ZKe_prof) * BigR * BB2_psi / BB2**2 * Bgrad_T_star     * Bgrad_Te     * xjac * theta * tstep &
                              + (ZKe_par_T-ZKe_prof) * BigR / BB2              * Bgrad_T_star_psi * Bgrad_Te     * xjac * theta * tstep &
                              + (ZKe_par_T-ZKe_prof) * BigR / BB2              * Bgrad_T_star     * Bgrad_Te_psi * xjac * theta * tstep &

                            + v * r0  * Vpar0 * (Te0_s * psi_t - Te0_t * psi_s)                                   * theta * tstep &
                            + v * Te0 * Vpar0 * (r0_s  * psi_t - r0_t * psi_s)                                    * theta * tstep &
                            + v * r0  * GAMMA * Te0 * (vpar0_s * psi_t - vpar0_t * psi_s)                         * theta * tstep &

                         + TG_num8 * 0.25d0 / BigR * vpar0**2                                                         &
                                   * Te0 * (r0_x * psi_y - r0_y * psi_x)                                              &
                                   * ( v_x * ps0_y -  v_y * ps0_x                  ) * xjac * theta * tstep * tstep   &
                         + TG_num8 * 0.25d0 / BigR * vpar0**2                                                         &
                                   * r0 * (Te0_x * psi_y - Te0_y * psi_x)                                             &
                                   * ( v_x * ps0_y -  v_y * ps0_x                  ) * xjac * theta * tstep * tstep   & 
                         + TG_num8 * 0.25d0 / BigR * vpar0**2                                                         &
                                   * Te0 * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                           &
                                   * ( v_x * psi_y -  v_y * psi_x ) * xjac * theta * tstep * tstep                    &
                         + TG_num8 * 0.25d0 / BigR * vpar0**2                                                         &
                                   * r0 * (Te0_x * ps0_y - Te0_y * ps0_x + F0 / BigR * Te0_p)                         &
                                   * ( v_x * psi_y -  v_y * psi_x ) * xjac * theta * tstep * tstep

                  amat_k(8,1) = - (ZKe_par_T-ZKe_prof) * BigR * BB2_psi / BB2**2 * Bgrad_T_k_star * Bgrad_Te     * xjac * theta * tstep &
                                + (ZKe_par_T-ZKe_prof) * BigR / BB2              * Bgrad_T_k_star * Bgrad_Te_psi * xjac * theta * tstep &
  
                        + TG_num8 * 0.25d0 / BigR * vpar0**2                                                       &
                                  * Te0 * (r0_x * psi_y - r0_y * psi_x)                                            &
                                  * (                            + F0 / BigR * v_p) * xjac * theta * tstep * tstep &
                        + TG_num8 * 0.25d0 / BigR * vpar0**2                                                       &
                                  * r0 * (Te0_x * psi_y - Te0_y * psi_x)                                           &
                                  * (                            + F0 / BigR * v_p) * xjac * theta * tstep * tstep


                  amat(8,2) = - v * r0 * BigR**2  * ( Te0_x * u_y - Te0_y * u_x)         * xjac * theta * tstep &
                              - v * Te0 * BigR**2 * ( r0_x * u_y - r0_y * u_x)           * xjac * theta * tstep &
                              - v * r0 * 2.d0* GAMMA * BigR * Te0 * u_y                  * xjac * theta * tstep &

                         + TG_num8 * 0.25d0 * BigR**2 * Te0* (r0_x * u_y - r0_y * u_x)               &
                                            * ( v_x * u0_y - v_y * u0_x) * xjac * theta*tstep*tstep  &
                         + TG_num8 * 0.25d0 * BigR**2 * r0* (Te0_x * u_y - Te0_y * u_x)              &
                                            * ( v_x * u0_y - v_y * u0_x) * xjac * theta*tstep*tstep  &
                         + TG_num8 * 0.25d0 * BigR**2 * Te0* (r0_x * u0_y - r0_y * u0_x)             &
                                            * ( v_x * u_y - v_y * u_x) * xjac * theta*tstep*tstep    &
                         + TG_num8 * 0.25d0 * BigR**2 * r0* (Te0_x * u0_y - Te0_y * u0_x)            &
                                            * ( v_x * u_y - v_y * u_x) * xjac * theta*tstep*tstep 

                  amat(8,3) = - v * (gamma-1.d0) * eta_T_ohm * 2.d0 * zj * zj0/(BigR**2.d0) * BigR * xjac * theta * tstep

                  amat(8,5) = v * rho * Te0_corr   * BigR * xjac * (1.d0 + zeta)     &
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


                         + TG_num8 * 0.25d0 * BigR**2 * Te0* (rho_x * u0_y - rho_y * u0_x)         &
                                   * ( v_x * u0_y - v_y * u0_x) * xjac * theta*tstep*tstep         &
                         + TG_num8 * 0.25d0 * BigR**2 * rho * (Te0_x * u0_y - Te0_y * u0_x)        &
                                   * ( v_x * u0_y - v_y * u0_x) * xjac* theta*tstep*tstep          &
                         + TG_num8 * 0.25d0 / BigR * vpar0**2                                      &
                                   * Te0 * (rho_x * ps0_y - rho_y * ps0_x )                        &
                                   * ( v_x * ps0_y -  v_y * ps0_x ) * xjac * theta * tstep * tstep &
                         + TG_num8 * 0.25d0 / BigR * vpar0**2                                      &
                                   * rho * (Te0_x * ps0_y - Te0_y * ps0_x + F0 / BigR * Te0_p)     &
                                   * ( v_x * ps0_y -  v_y * ps0_x ) * xjac * theta * tstep * tstep

                  amat_n(8,5) = + v * Te0  * F0 / BigR * Vpar0 * rho_p   * xjac * theta * tstep    &

                         + TG_num8 * 0.25d0 / BigR * vpar0**2                                                       &
                                   * Te0 * (                              + F0 / BigR * rho_p)                      &
                                   * ( v_x * ps0_y -  v_y * ps0_x                  ) * xjac * theta * tstep * tstep

                  amat_k(8,5) = + TG_num8 * 0.25d0 / BigR * vpar0**2                                                &
                                   * Te0 * (rho_x * ps0_y - rho_y * ps0_x                    )                      &
                                   * (                            + F0 / BigR * v_p) * xjac * theta * tstep * tstep &
                              + TG_num8 * 0.25d0 / BigR * vpar0**2                                                  &
                                   * rho * (Te0_x * ps0_y - Te0_y * ps0_x + F0 / BigR * Te0_p)                      &
                                   * (                            + F0 / BigR * v_p) * xjac * theta * tstep * tstep

                  amat_kn(8,5) = + TG_num8 * 0.25d0 / BigR * vpar0**2                 &
                                   * Te0 * (+ F0 / BigR * rho_p)                      &
                                   * (      + F0 / BigR * v_p) * xjac * theta * tstep * tstep

                  amat(8,6) = - v * BigR * ddTe_i_dTi * Ti                            * xjac * theta * tstep


                  amat(8,8) = v * r0_corr * Te  * BigR * xjac * (1.d0 + zeta)     &
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
  
                            + ZK_perp_num * (v_xx + v_x/BigR + v_yy)*(Te_xx + Te_x/BigR + Te_yy) * BigR * xjac * theta * tstep &

                            - v * Te * (gamma-1.d0) * deta_dT_ohm * (zj0 / BigR)**2.d0 * BigR * xjac * theta * tstep &

                            + TG_num8 * 0.25d0 * BigR**2 * Te* (r0_x * u0_y - r0_y * u0_x)         &
                                      * ( v_x * u0_y - v_y * u0_x) * xjac * theta * tstep * tstep  &
                            + TG_num8 * 0.25d0 * BigR**2 * r0* (Te_x * u0_y - Te_y * u0_x)         &
                                      * ( v_x * u0_y - v_y * u0_x) * xjac * theta * tstep * tstep  &
                            + TG_num8 * 0.25d0 / BigR * vpar0**2                                                       &
                                      * Te * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                          &
                                      * ( v_x * ps0_y -  v_y * ps0_x                  ) * xjac * theta * tstep * tstep &
                            + TG_num8 * 0.25d0 / BigR * vpar0**2                                                       &
                                      * r0 * (Te_x * ps0_y - Te_y * ps0_x             )                                &
                                      * ( v_x * ps0_y -  v_y * ps0_x                  ) * xjac * theta * tstep * tstep 

                  amat_k(8,8) = + (ZKe_par_T-ZKe_prof) * BigR / BB2 * Bgrad_T_k_star * Bgrad_Te_Te * xjac * theta * tstep  &
                                + dZKe_par_dT * Te     * BigR / BB2 * Bgrad_T_k_star * Bgrad_Te    * xjac * theta * tstep  &

                              + TG_num8 * 0.25d0 / BigR * vpar0**2                                                      &
                                  * Te * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                               &
                                  * (                                 + F0 / BigR * v_p) * xjac * theta * tstep * tstep &
                              + TG_num8 * 0.25d0 / BigR * vpar0**2                                                      &
                                  * r0 * (Te_x * ps0_y - Te_y * ps0_x                  )                                &
                                  * (                                + F0 / BigR * v_p) * xjac * theta * tstep * tstep

                  amat_n(8,8) = + (ZKe_par_T-ZKe_prof) * BigR / BB2 * Bgrad_T_star   * Bgrad_Te_Te_n  * xjac * theta * tstep &

                              + v * r0 * F0 / BigR * Vpar0 * Te_p                               * xjac * theta * tstep &
  
                              + TG_num8 * 0.25d0 / BigR * vpar0**2                       & 
                                * r0 * ( + F0 / BigR * Te_p) * ( v_x * ps0_y -  v_y * ps0_x ) * xjac * theta * tstep * tstep

                  amat_kn(8,8) = + (ZKe_par_T-ZKe_prof) * BigR / BB2 * Bgrad_T_k_star * Bgrad_Te_Te_n * xjac * theta * tstep &
                                 + ZKe_prof * BigR   * (v_p*Te_p /BigR**2 )                           * xjac * theta * tstep &

                              + TG_num8 * 0.25d0 / BigR * vpar0**2 &
                                * r0 * ( + F0 / BigR * Te_p) * ( + F0 / BigR * v_p)           * xjac * theta * tstep * tstep
 
                  amat(8,7) = + v * r0 * F0 / BigR * Vpar * Te0_p                             * xjac * theta * tstep &
                              + v * Te0 * F0 / BigR * Vpar * r0_p                             * xjac * theta * tstep &

                              + v * r0  * Vpar * (Te0_s * ps0_t - Te0_t * ps0_s)                      * theta * tstep &
                              + v * Te0 * Vpar * (r0_s * ps0_t  - r0_t * ps0_s)                       * theta * tstep & 

                              + v * r0 * GAMMA * Te0 * (vpar_s * ps0_t - vpar_t * ps0_s)              * theta * tstep &
  
                              + TG_num8 * 0.25d0 / BigR * 2.d0 * vpar0*vpar &
                                  * Te0 * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                               &
                                  * ( v_x * ps0_y -  v_y * ps0_x                        ) * xjac * theta * tstep * tstep &
                              + TG_num8 * 0.25d0 / BigR * 2.d0 * vpar0*vpar &
                                  * r0 * (Te0_x * ps0_y - Te0_y * ps0_x + F0 / BigR * Te0_p)                             &
                                  * ( v_x * ps0_y -  v_y * ps0_x                        ) * xjac * theta * tstep * tstep 

                  amat_k(8,7) =  &
                        + TG_num8 * 0.25d0 / BigR * 2.d0 * vpar0*vpar &
                                  * Te0 * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                          &
                                  * (                            + F0 / BigR * v_p) * xjac * theta * tstep * tstep &
                        + TG_num8 * 0.25d0 / BigR * 2.d0 * vpar0*vpar &
                                  * r0 * (Te0_x * ps0_y - Te0_y * ps0_x + F0 / BigR * Te0_p)                          &
                                  * (                            + F0 / BigR * v_p) * xjac * theta * tstep * tstep 

                  amat_n(8,7) = + v * r0 * GAMMA * Te0 * F0 / BigR * vpar_p          * xjac * theta * tstep

 
                  !###################################################################################################
                  !# end equations                                                                                   #
                  !###################################################################################################

                  if (use_fft) then
                    index_kl = n_var*n_degrees*(k-1) + n_var*(l-1) + 1
                  else
                    index_kl = n_tor_local*n_var*n_degrees*(k-1) + n_tor_local*n_var*(l-1) + in - n_tor_start +1
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
        enddo

      enddo

    endif ! apply fft (or not)

  enddo ! j loop n_degrees
enddo ! i loop (n_vertex)

if (.NOT. use_fft) return

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
