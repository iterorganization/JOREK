module mod_elt_matrix_fft
implicit none
contains
subroutine element_matrix_fft(element, nodes, xpoint2, xcase2, R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint, Z_xpoint, ELM, RHS, tid, &
  ELM_p, ELM_n, ELM_k, ELM_kn, RHS_p, RHS_k,  eq_g, eq_s, eq_t, eq_p, eq_ss, eq_st, eq_tt, delta_g, delta_s, delta_t, i_tor_min, i_tor_max, aux_nodes)
!---------------------------------------------------------------
! calculates the matrix contribution of one element
!
! TO BE CHECKED (see element_matrix)
!
!---------------------------------------------------------------
use constants
use mod_parameters
use data_structure
use gauss
use basis_at_gaussian
use phys_module
use tr_module
use diffusivities, only: get_dperp, get_zkperp
use corr_neg
use equil_info, only: get_psi_n
use mod_semianalytical
use mod_equations
use mod_chi
use mod_sources

implicit none
 
type (type_element), intent(in)           :: element
type (type_node),    intent(in)           :: nodes(n_vertex_max)
type (type_node),    intent(in), optional :: aux_nodes(n_vertex_max) ! particle moments

#define DIM0 n_tor*n_vertex_max*(n_order+1)*n_var

real*8, dimension (DIM0,DIM0)       :: ELM
real*8, dimension (DIM0)            :: RHS
integer                , intent(in) :: tid
integer                , intent(in) :: i_tor_min, i_tor_max

integer    :: i, j, ms, mt, mp, k, l, index_ij, index_kl, index, index_k, index_m, m, ik, xcase2
integer    :: in, im, ij1, ij2, ij3, ij4, ij5, ij6, kl1, kl2, kl3, kl4, kl5, kl6, n_tor_start, n_tor_end, n_tor_local
real*8     :: wst, xjac, xjac_x, xjac_y, xjac_s, xjac_t, BigR, r2, phi, eps_cyl
real*8     :: R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint(2), Z_xpoint(2), dj_dpsi, dj_dz
real*8     :: Bgrad_rho_star,     Bgrad_rho,     Bgrad_T_star,  Bgrad_T, BB2
real*8     :: Bgrad_rho_star_psi, Bgrad_rho_psi, Bgrad_rho_rho, Bgrad_T_star_psi, Bgrad_T_psi, Bgrad_T_T, BB2_psi
real*8     :: Bgrad_rho_rho_n, Bgrad_T_T_n, Bgrad_rho_k_star, Bgrad_T_k_star
real*8     :: D_prof, ZK_prof, psi_norm

real*8     :: x_p_x, x_p_y, y_p_x, y_p_y, v_px, v_py, u_px, u_py
real*8     :: v, v_x, v_y, v_s, v_t, v_p, v_ss, v_st, v_tt, v_xx, v_yy, v_xs, v_ys, v_xt, v_yt, v_xy
real*8     :: ps0, ps0_x, ps0_y, ps0_p,ps0_s,ps0_t,  zj0, zj0_x, zj0_y, zj0_p, zj0_s, zj0_t
real*8     :: u0, u0_x, u0_y, u0_p, u0_s, u0_t,  w0, w0_x, w0_y, w0_p, w0_s, w0_t
real*8     :: r0, r0_x, r0_y, r0_p, r0_s, r0_t,  r0_hat, r0_x_hat, r0_y_hat, T0, T0_x, T0_y, T0_p, T0_s, T0_t
real*8     :: psi, psi_x, psi_y, psi_p, psi_s, psi_t, psi_ss, psi_st, psi_tt, psi_xs, psi_ys, psi_xt, psi_yt, psi_xx, psi_yy
real*8     :: zj, zj_x, zj_y, zj_p, zj_s, zj_t
real*8     :: u, u_x, u_y, u_p, u_s, u_t, w, w_x, w_y, w_p, w_s, w_t, w_xx, w_yy
real*8     :: rho, rho_x, rho_y, rho_s, rho_t, rho_p, rho_hat, rho_x_hat, rho_y_hat, T, T_x, T_y, T_s, T_t, T_p
real*8     :: w0_xs, w0_xt, w0_ys, w0_yt, w0_xx, w0_yy, w0_xy, w0_ss, w0_tt, w0_st, P0, P0_s, P0_t, P0_x, P0_y
real*8     :: BigR_x, vv2, eta_T, visco_T, deta_dT, d2eta_d2T, dvisco_dT
real*8     :: theta, zeta, reta
logical    :: xpoint2, use_fft

real*8, dimension(4) :: rhs_ij_1, rhs_ij_2, rhs_ij_3, rhs_ij_4, rhs_ij_5, rhs_ij_6
real*8, dimension(4) :: amat_11, amat_12, amat_13, amat_16, amat_21, amat_22, amat_23, amat_24, amat_25, amat_26, amat_31, amat_33
real*8, dimension(4) :: amat_42, amat_44, amat_51, amat_52, amat_55, amat_61, amat_62, amat_63, amat_65, amat_66

integer*8  :: plan
real*8     :: in_fft(1:n_plane)
complex*16 :: out_fft(1:n_plane)
INTEGER    :: FFTW_FORWARD,  FFTW_BACKWARD, FFTW_ESTIMATE
PARAMETER (FFTW_FORWARD=-1,FFTW_BACKWARD=+1, FFTW_ESTIMATE=64)

#define DIM1 n_plane
#define DIM2 1:n_vertex_max*n_var*(n_order+1)

real*8, dimension(DIM1, DIM2, DIM2) :: ELM_p
real*8, dimension(DIM1, DIM2, DIM2) :: ELM_n
real*8, dimension(DIM1, DIM2, DIM2) :: ELM_k
real*8, dimension(DIM1, DIM2, DIM2) :: ELM_kn
real*8, dimension(DIM1, DIM2)       :: RHS_p
real*8, dimension(DIM1, DIM2)       :: RHS_k

real*8, dimension(n_plane,n_gauss,n_gauss) :: x_g, x_s, x_t, x_p, x_ss, x_st, x_tt, x_sp, x_tp, x_pp
real*8, dimension(n_plane,n_gauss,n_gauss) :: y_g, y_s, y_t, y_p, y_ss, y_st, y_tt, y_sp, y_tp, y_pp
real*8, dimension(n_plane,n_gauss,n_gauss) :: current_source, particle_source, heat_source
real*8, dimension(n_gauss, n_gauss)        :: s_norm

real*8, dimension(n_plane,n_var,n_gauss,n_gauss) :: eq_g, eq_s, eq_t, eq_p, eq_ss, eq_st, eq_tt, eq_pp, eq_sp, eq_tp
real*8, dimension(n_plane,n_var,n_gauss,n_gauss) :: delta_g, delta_s, delta_t, delta_p

real*8, dimension(:,:,:,:,:), pointer :: eq
real*8, dimension(n_var)              :: eq_px, eq_py

real*8, dimension(4), parameter  :: vr = (/1.,0.,1.,0./), vp = (/0.,1.,0.,1./), dur = (/1.,1.,0.,0./), dup = (/0.,0.,1.,1./)

eq => thread_eq(tid)%eq

ELM_p = 0.d0; ELM_n = 0.d0; ELM_k = 0.d0; ELM_kn = 0.d0
RHS_p = 0.d0; RHS_k = 0.d0
ELM   = 0.d0; RHS   = 0.d0

! --- Take time evolution parameters from phys_module
theta = time_evol_theta
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

if (use_fft) then
  ! In case of FFT, don't loop over toroidal harmonics:
  n_tor_start = 1
  n_tor_end   = 1
else
  n_tor_start = i_tor_min
  n_tor_end   = i_tor_max
end if

n_tor_local = n_tor_end - n_tor_start + 1

!---------------------------------------------------- value of (x,y) and derivatives on Gaussian points
#ifdef altcs
psieq = 0.d0; psieq_s = 0.d0; psieq_t = 0.d0
#endif
x_g  = 0.d0; x_s   = 0.d0; x_t   = 0.d0; x_p = 0.d0; x_st  = 0.d0; x_ss  = 0.d0; x_tt  = 0.d0; x_sp = 0.d0; x_tp = 0.d0; x_pp = 0.d0;
y_g  = 0.d0; y_s   = 0.d0; y_t   = 0.d0; y_p = 0.d0; y_st  = 0.d0; y_ss  = 0.d0; y_tt  = 0.d0; y_sp = 0.d0; y_tp = 0.d0; y_pp = 0.d0;
eq_g = 0.d0; eq_s  = 0.d0; eq_t  = 0.d0; eq_st = 0.d0; eq_ss = 0.d0; eq_tt = 0.d0;
eq_p = 0.d0; eq_pp = 0.d0; eq_sp = 0.d0; eq_tp = 0.d0

delta_g = 0.d0; delta_s = 0.d0; delta_t = 0.d0; delta_p = 0.d0

eq = 0.d0

current_source  = 0.d0
particle_source = 0.d0
heat_source     = 0.d0

do i=1,n_vertex_max
  do j=1,n_order+1
    do ms=1, n_gauss
      do mt=1, n_gauss
#ifdef altcs
        psieq(ms,mt)   = psieq(ms,mt)   + nodes(i)%psi_eq(j)*element%size(i,j)*H(i,j,ms,mt)
        psieq_s(ms,mt) = psieq_s(ms,mt) + nodes(i)%psi_eq(j)*element%size(i,j)*H_s(i,j,ms,mt)
        psieq_t(ms,mt) = psieq_t(ms,mt) + nodes(i)%psi_eq(j)*element%size(i,j)*H_t(i,j,ms,mt)
#endif
        s_norm(ms, mt) = s_norm(ms, mt) + nodes(i)%s_eq(j)*element%size(i,j)*H(i,j,ms,mt)

        do mp=1,n_plane
          do in=1,n_coord_tor
            x_g(mp,ms,mt)  = x_g(mp,ms,mt)  + nodes(i)%x(in,j,1)*element%size(i,j)*H(i,j,ms,mt)   *HZ_coord(in,mp)
            x_s(mp,ms,mt)  = x_s(mp,ms,mt)  + nodes(i)%x(in,j,1)*element%size(i,j)*H_s(i,j,ms,mt) *HZ_coord(in,mp)
            x_t(mp,ms,mt)  = x_t(mp,ms,mt)  + nodes(i)%x(in,j,1)*element%size(i,j)*H_t(i,j,ms,mt) *HZ_coord(in,mp)
            x_p(mp,ms,mt)  = x_p(mp,ms,mt)  + nodes(i)%x(in,j,1)*element%size(i,j)*H(i,j,ms,mt)   *HZ_coord_p(in,mp)
            x_ss(mp,ms,mt) = x_ss(mp,ms,mt) + nodes(i)%x(in,j,1)*element%size(i,j)*H_ss(i,j,ms,mt)*HZ_coord(in,mp)
            x_st(mp,ms,mt) = x_st(mp,ms,mt) + nodes(i)%x(in,j,1)*element%size(i,j)*H_st(i,j,ms,mt)*HZ_coord(in,mp)
            x_tt(mp,ms,mt) = x_tt(mp,ms,mt) + nodes(i)%x(in,j,1)*element%size(i,j)*H_tt(i,j,ms,mt)*HZ_coord(in,mp)
            x_sp(mp,ms,mt) = x_sp(mp,ms,mt) + nodes(i)%x(in,j,1)*element%size(i,j)*H_s(i,j,ms,mt) *HZ_coord_p(in,mp)
            x_tp(mp,ms,mt) = x_tp(mp,ms,mt) + nodes(i)%x(in,j,1)*element%size(i,j)*H_t(i,j,ms,mt) *HZ_coord_p(in,mp)
            x_pp(mp,ms,mt) = x_pp(mp,ms,mt) + nodes(i)%x(in,j,1)*element%size(i,j)*H(i,j,ms,mt)   *HZ_coord_pp(in,mp)

            y_g(mp,ms,mt)  = y_g(mp,ms,mt)  + nodes(i)%x(in,j,2)*element%size(i,j)*H(i,j,ms,mt)   *HZ_coord(in,mp)
            y_s(mp,ms,mt)  = y_s(mp,ms,mt)  + nodes(i)%x(in,j,2)*element%size(i,j)*H_s(i,j,ms,mt) *HZ_coord(in,mp)
            y_t(mp,ms,mt)  = y_t(mp,ms,mt)  + nodes(i)%x(in,j,2)*element%size(i,j)*H_t(i,j,ms,mt) *HZ_coord(in,mp)
            y_p(mp,ms,mt)  = y_p(mp,ms,mt)  + nodes(i)%x(in,j,2)*element%size(i,j)*H(i,j,ms,mt)   *HZ_coord_p(in,mp)
            y_ss(mp,ms,mt) = y_ss(mp,ms,mt) + nodes(i)%x(in,j,2)*element%size(i,j)*H_ss(i,j,ms,mt)*HZ_coord(in,mp)
            y_st(mp,ms,mt) = y_st(mp,ms,mt) + nodes(i)%x(in,j,2)*element%size(i,j)*H_st(i,j,ms,mt)*HZ_coord(in,mp)
            y_tt(mp,ms,mt) = y_tt(mp,ms,mt) + nodes(i)%x(in,j,2)*element%size(i,j)*H_tt(i,j,ms,mt)*HZ_coord(in,mp)
            y_sp(mp,ms,mt) = y_sp(mp,ms,mt) + nodes(i)%x(in,j,2)*element%size(i,j)*H_s(i,j,ms,mt) *HZ_coord_p(in,mp)
            y_tp(mp,ms,mt) = y_tp(mp,ms,mt) + nodes(i)%x(in,j,2)*element%size(i,j)*H_t(i,j,ms,mt) *HZ_coord_p(in,mp)
            y_pp(mp,ms,mt) = y_pp(mp,ms,mt) + nodes(i)%x(in,j,2)*element%size(i,j)*H(i,j,ms,mt)   *HZ_coord_pp(in,mp)
          end do

          do k=1,n_var
            do in=1,n_tor
              eq_g(mp,k,ms,mt) = eq_g(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H(i,j,ms,mt)  * HZ(in,mp)
              eq_s(mp,k,ms,mt) = eq_s(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_s(i,j,ms,mt)* HZ(in,mp)
              eq_t(mp,k,ms,mt) = eq_t(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_t(i,j,ms,mt)* HZ(in,mp)
              eq_p(mp,k,ms,mt) = eq_p(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H(i,j,ms,mt)  * HZ_p(in,mp)
              eq_pp(mp,k,ms,mt) = eq_pp(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H(i,j,ms,mt)   * HZ_pp(in,mp)
              eq_sp(mp,k,ms,mt) = eq_sp(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_s(i,j,ms,mt) * HZ_p(in,mp)
              eq_tp(mp,k,ms,mt) = eq_tp(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_t(i,j,ms,mt) * HZ_p(in,mp)
              eq_ss(mp,k,ms,mt) = eq_ss(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_ss(i,j,ms,mt)* HZ(in,mp)
              eq_st(mp,k,ms,mt) = eq_st(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_st(i,j,ms,mt)* HZ(in,mp)
              eq_tt(mp,k,ms,mt) = eq_tt(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_tt(i,j,ms,mt)* HZ(in,mp)

              delta_g(mp,k,ms,mt) = delta_g(mp,k,ms,mt) + nodes(i)%deltas(in,j,k) * element%size(i,j) * H(i,j,ms,mt)   * HZ(in,mp)
              delta_s(mp,k,ms,mt) = delta_s(mp,k,ms,mt) + nodes(i)%deltas(in,j,k) * element%size(i,j) * H_s(i,j,ms,mt) * HZ(in,mp)
              delta_t(mp,k,ms,mt) = delta_t(mp,k,ms,mt) + nodes(i)%deltas(in,j,k) * element%size(i,j) * H_t(i,j,ms,mt) * HZ(in,mp)
              delta_p(mp,k,ms,mt) = delta_p(mp,k,ms,mt) + nodes(i)%deltas(in,j,k) * element%size(i,j) * H(i,j,ms,mt)   * HZ_p(in,mp)
            enddo
          enddo

          if (keep_current_prof) then
            do in=1,n_tor
              current_source(mp,ms,mt) = current_source(mp,ms,mt) + nodes(i)%j_source(in,j)*element%size(i,j)*H(i,j,ms,mt)*HZ(in,mp)
            end do
          end if

#ifdef altcs
          call sources(xpoint2, xcase2, y_g(mp,ms,mt), Z_xpoint, s_norm(ms,mt),0.0,1.0,particle_source(mp,ms,mt),heat_source(mp,ms,mt))
#else
          call sources(xpoint2, xcase2, y_g(mp,ms,mt), Z_xpoint, s_norm(ms,mt), & !eq_g(mp,var_Psi,ms,mt),
                       0.0,1.0,particle_source(mp,ms,mt),heat_source(mp,ms,mt))
#endif
        enddo
      enddo
    enddo
  enddo
enddo

if (eta .ne. 0.d0) then
  reta = eta_ohmic/eta
else
  reta = 0.d0
end if

! changes deltas for variable time steps
delta_g = delta_g*tstep/tstep_prev
delta_s = delta_s*tstep/tstep_prev
delta_t = delta_t*tstep/tstep_prev

!--------------------------------------------------- sum over the Gaussian integration points
do ms=1, n_gauss
  do mt=1, n_gauss
    wst = wgauss(ms)*wgauss(mt)

#ifdef altcs
    psieq_x = ( y_t(1,ms,mt)*psieq_s(ms,mt) - y_s(1,ms,mt)*psieq_t(ms,mt))/xjac
    psieq_y = (-x_t(1,ms,mt)*psieq_s(ms,mt) + x_s(1,ms,mt)*psieq_t(ms,mt))/xjac
#endif

    do mp = 1, n_plane
      phi = 2.d0*pi*float(mp-1)/float(n_plane*n_period)
      xjac    = x_s(mp,ms,mt)*y_t(mp,ms,mt)  - x_t(mp,ms,mt)*y_s(mp,ms,mt)
      xjac_x  = (x_ss(mp,ms,mt)*y_t(mp,ms,mt)**2 - y_ss(mp,ms,mt)*x_t(mp,ms,mt)*y_t(mp,ms,mt) - 2.d0*x_st(mp,ms,mt)*y_s(mp,ms,mt)*y_t(mp,ms,mt) &
              + y_st(mp,ms,mt)*(x_s(mp,ms,mt)*y_t(mp,ms,mt) + x_t(mp,ms,mt)*y_s(mp,ms,mt))                                                      &
              + x_tt(mp,ms,mt)*y_s(mp,ms,mt)**2 - y_tt(mp,ms,mt)*x_s(mp,ms,mt)*y_s(mp,ms,mt)) / xjac
      xjac_y  = (y_tt(mp,ms,mt)*x_s(mp,ms,mt)**2 - x_tt(mp,ms,mt)*y_s(mp,ms,mt)*x_s(mp,ms,mt) - 2.d0*y_st(mp,ms,mt)*x_t(mp,ms,mt)*x_s(mp,ms,mt) &
              + x_st(mp,ms,mt)*(y_t(mp,ms,mt)*x_s(mp,ms,mt) + y_s(mp,ms,mt)*x_t(mp,ms,mt))                                                      &
              + y_ss(mp,ms,mt)*x_t(mp,ms,mt)**2 - x_ss(mp,ms,mt)*y_t(mp,ms,mt)*x_t(mp,ms,mt)) / xjac
      x_p_x = (x_sp(mp,ms,mt)*y_t(mp,ms,mt) - x_tp(mp,ms,mt)*y_s(mp,ms,mt))/xjac
      x_p_y = (x_tp(mp,ms,mt)*x_s(mp,ms,mt) - x_sp(mp,ms,mt)*x_t(mp,ms,mt))/xjac
      y_p_x = (y_sp(mp,ms,mt)*y_t(mp,ms,mt) - y_tp(mp,ms,mt)*y_s(mp,ms,mt))/xjac
      y_p_y = (y_tp(mp,ms,mt)*x_s(mp,ms,mt) - y_sp(mp,ms,mt)*x_t(mp,ms,mt))/xjac
      BigR = x_g(mp,ms,mt)

      ! Values at current time step (u^n)
      eq(1:n_var,0,0,0,1) = eq_g(mp,:,ms,mt)
      eq(1:n_var,1,0,0,1) = (y_t(mp,ms,mt)*eq_s(mp,:,ms,mt) - y_s(mp,ms,mt)*eq_t(mp,:,ms,mt))/xjac
      eq(1:n_var,0,1,0,1) = (-x_t(mp,ms,mt)*eq_s(mp,:,ms,mt) + x_s(mp,ms,mt)*eq_t(mp,:,ms,mt))/xjac
      eq(1:n_var,0,0,1,1) = eq_p(mp,:,ms,mt) - eq(1:n_var,1,0,0,1)*x_p(mp,ms,mt) - eq(1:n_var,0,1,0,1)*y_p(mp,ms,mt)
      eq(1:n_var,2,0,0,1) = (eq_ss(mp,:,ms,mt)*y_t(mp,ms,mt)**2 - 2.d0*eq_st(mp,:,ms,mt)*y_s(mp,ms,mt)*y_t(mp,ms,mt) &
                          + eq_tt(mp,:,ms,mt)*y_s(mp,ms,mt)**2                                                       &
                          + eq_s(mp,:,ms,mt)*(y_st(mp,ms,mt)*y_t(mp,ms,mt) - y_tt(mp,ms,mt)*y_s(mp,ms,mt))           &
                          + eq_t(mp,:,ms,mt)*(y_st(mp,ms,mt)*y_s(mp,ms,mt) - y_ss(mp,ms,mt)*y_t(mp,ms,mt)))/xjac**2  &
                          - xjac_x*(eq_s(mp,:,ms,mt)*y_t(mp,ms,mt) - eq_t(mp,:,ms,mt)*y_s(mp,ms,mt))/xjac**2
      eq(1:n_var,0,2,0,1) = (eq_ss(mp,:,ms,mt)*x_t(mp,ms,mt)**2 - 2.d0*eq_st(mp,:,ms,mt)*x_s(mp,ms,mt)*x_t(mp,ms,mt) &
                          + eq_tt(mp,:,ms,mt)*x_s(mp,ms,mt)**2                                                       &
                          + eq_s(mp,:,ms,mt)*(x_st(mp,ms,mt)*x_t(mp,ms,mt) - x_tt(mp,ms,mt)*x_s(mp,ms,mt))           &
                          + eq_t(mp,:,ms,mt)*(x_st(mp,ms,mt)*x_s(mp,ms,mt) - x_ss(mp,ms,mt)*x_t(mp,ms,mt)))/xjac**2  &
                          - xjac_y*(-eq_s(mp,:,ms,mt)*x_t(mp,ms,mt) + eq_t(mp,:,ms,mt)*x_s(mp,ms,mt))/xjac**2
      eq(1:n_var,1,1,0,1) = (-eq_ss(mp,:,ms,mt)*y_t(mp,ms,mt)*x_t(mp,ms,mt) - eq_tt(mp,:,ms,mt)*x_s(mp,ms,mt)*y_s(mp,ms,mt) &
                          + eq_st(mp,:,ms,mt)*(y_s(mp,ms,mt)*x_t(mp,ms,mt) + y_t(mp,ms,mt)*x_s(mp,ms,mt))                   &
                          - eq_s(mp,:,ms,mt)*(x_st(mp,ms,mt)*y_t(mp,ms,mt) - x_tt(mp,ms,mt)*y_s(mp,ms,mt))                  &
                          - eq_t(mp,:,ms,mt)*(x_st(mp,ms,mt)*y_s(mp,ms,mt) - x_ss(mp,ms,mt)*y_t(mp,ms,mt)))/xjac**2         &
                          - xjac_x*(-eq_s(mp,:,ms,mt)*x_t(mp,ms,mt) + eq_t(mp,:,ms,mt)*x_s(mp,ms,mt))/xjac**2
      eq_px               = (y_t(mp,ms,mt)*eq_sp(mp,:,ms,mt) - y_s(mp,ms,mt)*eq_tp(mp,:,ms,mt))/xjac
      eq_py               = (-x_t(mp,ms,mt)*eq_sp(mp,:,ms,mt) + x_s(mp,ms,mt)*eq_tp(mp,:,ms,mt))/xjac
      ! Second derivatives wrt phi not implemented in FFT, and not necessary for stabilization (2nd derivatives only appear in hyperdissipation terms)
!      eq(1:n_var,0,0,2,1) = eq_pp(mp,:,ms,mt) - x_pp(mp,ms,mt)*eq(1:n_var,1,0,0,1) - 2.d0*(x_p(mp,ms,mt)*eq_px + y_p(mp,ms,mt)*eq_py)            &
!                          - y_pp(mp,ms,mt)*eq(1:n_var,0,1,0,1) + 2.d0*(x_p(mp,ms,mt)*x_p_x*eq(1:n_var,1,0,0,1) + x_p(mp,ms,mt)*y_p_x*eq(1:n_var,0,1,0,1)&
!                          + y_p(mp,ms,mt)*x_p_y*eq(1:n_var,1,0,0,1) + y_p(mp,ms,mt)*y_p_y*eq(1:n_var,0,1,0,1)) + x_p(mp,ms,mt)**2*eq(1:n_var,2,0,0,1)   &
!                          + 2.d0*x_p(mp,ms,mt)*y_p(mp,ms,mt)*eq(1:n_var,1,1,0,1) + y_p(mp,ms,mt)**2*eq(1:n_var,0,2,0,1)
      eq(1:n_var,1,0,1,1) = eq_px - x_p_x*eq(1:n_var,1,0,0,1) - x_p(mp,ms,mt)*eq(1:n_var,2,0,0,1) - y_p_x*eq(1:n_var,0,1,0,1) &
                          - y_p(mp,ms,mt)*eq(1:n_var,1,1,0,1)
      eq(1:n_var,0,1,1,1) = eq_py - x_p_y*eq(1:n_var,1,0,0,1) - x_p(mp,ms,mt)*eq(1:n_var,1,1,0,1) - y_p_y*eq(1:n_var,0,1,0,1) &
                          - y_p(mp,ms,mt)*eq(1:n_var,0,2,0,1)

      ! Increments since previous time step (delta_u^(n-1))
      eq(n_var+1:2*n_var,0,0,0,1) = delta_g(mp,:,ms,mt)
      eq(n_var+1:2*n_var,1,0,0,1) = (y_t(mp,ms,mt)*delta_s(mp,:,ms,mt) - y_s(mp,ms,mt)*delta_t(mp,:,ms,mt))/xjac
      eq(n_var+1:2*n_var,0,1,0,1) = (-x_t(mp,ms,mt)*delta_s(mp,:,ms,mt) + x_s(mp,ms,mt)*delta_t(mp,:,ms,mt))/xjac
      eq(n_var+1:2*n_var,0,0,1,1) = delta_p(mp,:,ms,mt) - eq(n_var+1:2*n_var,1,0,0,1)*x_p(mp,ms,mt) - eq(n_var+1:2*n_var,0,1,0,1)*y_p(mp,ms,mt)

      ! Vacuum scalar magnetic potential (chi) and field (grad chi)
      eq(2*n_var+3,:,:,:,1) = get_chi(x_g(mp,ms,mt),y_g(mp,ms,mt),phi)
      
      do m=2,4
        eq(1:n_var,:,:,:,m) = eq(1:n_var,:,:,:,1)
        eq(n_var+1:2*n_var,:,:,:,m) = eq(n_var+1:2*n_var,:,:,:,1)
        eq(2*n_var+3,:,:,:,m) = eq(2*n_var+3,:,:,:,1)
      end do
      
      eq(2*n_var+4,0,0,0,:) = x_g(mp,ms,mt); eq(2*n_var+4,1,0,0,:) = 1.d0 ! Cylindrical R coordinate

      psi_norm = s_norm(ms,mt)

      ! The Psi in the equations differs by a factor of F0 from the normal JOREK Psi
      eq(1,:,:,:,:) = eq(1,:,:,:,:)/F0
      eq(n_var+1,:,:,:,:) = eq(n_var+1,:,:,:,:)/F0
      eq(3,:,:,:,:) = eq(3,:,:,:,:)/F0
      eq(n_var+3,:,:,:,:) = eq(n_var+3,:,:,:,:)/F0

      eq(2*n_var+5,0,0,0,:) = get_dperp(psi_norm)       ! D_perp
      eq(2*n_var+6,0,0,0,:) = get_zkperp(psi_norm)      ! k_perp
      eq(2*n_var+7,0,0,0,:) = particle_source(mp,ms,mt)   ! S_rho
      eq(2*n_var+8,0,0,0,:) = heat_source(mp,ms,mt)       ! S_e
      eq(2*n_var+9,0,0,0,:) = current_source(mp,ms,mt)/F0 ! S_j
      ! Resistivity
      if (eta_T_dependent) then
        eq(2*n_var+10,0,0,0,:) = eta*(corr_neg_temp(eq(6,0,0,0,1))/T_0)**(-1.5d0)               ! eta
        eq(2*n_var+11,0,0,0,:) = -1.5d0*eta*corr_neg_temp(eq(6,0,0,0,1))**(-2.5d0)*T_0**(1.5d0) ! deta/dT
      else
        eq(2*n_var+10,0,0,0,:) = eta
        eq(2*n_var+11,0,0,0,:) = 0.d0
      end if
      ! Viscosity
      if (visco_T_dependent) then
        eq(2*n_var+12,0,0,0,:) = visco*(corr_neg_temp(eq(6,0,0,0,1))/T_0)**(-1.5d0)               ! visco
        eq(2*n_var+13,0,0,0,:) = -1.5d0*visco*corr_neg_temp(eq(6,0,0,0,1))**(-2.5d0)*T_0**(1.5d0) ! dvisco/dT
      else
        eq(2*n_var+12,0,0,0,:) = visco
        eq(2*n_var+13,0,0,0,:) = 0.d0
      end if
      ! Parallel thermal conductivity
      if (zkpar_T_dependent) then
        eq(2*n_var+14,0,0,0,:) = zk_par*(corr_neg_temp(eq(6,0,0,0,1))/T_0)**(2.5d0)               ! k_par
        eq(2*n_var+15,0,0,0,:) = 2.5d0*zk_par*corr_neg_temp(eq(6,0,0,0,1))**(1.5d0)*T_0**(-2.5d0) ! dk_par_dT
        if (eq(2*n_var+14,0,0,0,1) .gt. zk_par_max) then
          eq(2*n_var+14,0,0,0,:) = zk_par_max
          eq(2*n_var+15,0,0,0,:) = 0.d0
        end if
      else
        eq(2*n_var+14,0,0,0,:) = zk_par
        eq(2*n_var+15,0,0,0,:) = 0.d0
      end if

      ! Auxiliary variables (aux)
#ifdef DEBUG
      eq(2*n_var+16,0,0,0,:) = eval(thread_eq(tid)%aBv2seq); eq(2*n_var+16,1,0,0,:) = eval(thread_eq(tid)%aBv2xseq)
      eq(2*n_var+16,0,1,0,:) = eval(thread_eq(tid)%aBv2yseq); eq(2*n_var+16,0,0,1,:) = eval(thread_eq(tid)%aBv2pseq)
      eq(2*n_var+17,0,0,0,:) = eval(thread_eq(tid)%aB2seq)
#else
#include "aux_unreadable.h"
#endif

      do i=1,n_vertex_max
        do j=1,n_order+1
          do im=n_tor_start,n_tor_end
            if (use_fft) then
              index_ij =       n_var*(n_order+1)*(i-1) +       n_var*(j-1) + 1
            else
              index_ij = n_tor_local*n_var*(n_order+1)*(i-1) + n_tor_local*n_var*(j-1) + im - n_tor_start +1 
            end if
            
            ! Test function (v)
            eq(2*n_var+1,0,0,0,:) =  H(i,j,ms,mt)*element%size(i,j)*vr
            eq(2*n_var+1,1,0,0,:) = (y_t(mp,ms,mt)*h_s(i,j,ms,mt) - y_s(mp,ms,mt)*h_t(i,j,ms,mt))*element%size(i,j)*vr/xjac
            eq(2*n_var+1,0,1,0,:) = (-x_t(mp,ms,mt)*h_s(i,j,ms,mt) + x_s(mp,ms,mt)*h_t(i,j,ms,mt))*element%size(i,j)*vr/xjac
            eq(2*n_var+1,0,0,1,:) = H(i,j,ms,mt)*element%size(i,j)*vp - eq(2*n_var+1,1,0,0,:)*x_p(mp,ms,mt) - eq(2*n_var+1,0,1,0,:)*y_p(mp,ms,mt)
            eq(2*n_var+1,2,0,0,:) = (h_ss(i,j,ms,mt)*y_t(mp,ms,mt)**2 - 2.d0*h_st(i,j,ms,mt)*y_s(mp,ms,mt)*y_t(mp,ms,mt)                       &
	                              + h_tt(i,j,ms,mt)*y_s(mp,ms,mt)**2                                                                           &
	                              + h_s(i,j,ms,mt)*(y_st(mp,ms,mt)*y_t(mp,ms,mt) - y_tt(mp,ms,mt)*y_s(mp,ms,mt))                               &
                                  + h_t(i,j,ms,mt)*(y_st(mp,ms,mt)*y_s(mp,ms,mt) - y_ss(mp,ms,mt)*y_t(mp,ms,mt)))*element%size(i,j)*vr/xjac**2 &
                                  - xjac_x*(h_s(i,j,ms,mt)*y_t(mp,ms,mt) - h_t(i,j,ms,mt)*y_s(mp,ms,mt))*element%size(i,j)*vr/xjac**2
            eq(2*n_var+1,0,2,0,:) = (h_ss(i,j,ms,mt)*x_t(mp,ms,mt)**2 - 2.d0*h_st(i,j,ms,mt)*x_s(mp,ms,mt)*x_t(mp,ms,mt)                    &
                                  + h_tt(i,j,ms,mt)*x_s(mp,ms,mt)**2                                                                        &
                                  + h_s(i,j,ms,mt)*(x_st(mp,ms,mt)*x_t(mp,ms,mt) - x_tt(mp,ms,mt)*x_s(mp,ms,mt))                            &
                                  + h_t(i,j,ms,mt)*(x_st(mp,ms,mt)*x_s(mp,ms,mt) - x_ss(mp,ms,mt)*x_t(mp,ms,mt)))*element%size(i,j)*vr/xjac**2 &
           	                      - xjac_y*(-h_s(i,j,ms,mt)*x_t(mp,ms,mt) + h_t(i,j,ms,mt)*x_s(mp,ms,mt))*element%size(i,j)*vr/xjac**2
            eq(2*n_var+1,1,1,0,:) = (-h_ss(i,j,ms,mt)*y_t(mp,ms,mt)*x_t(mp,ms,mt) - h_tt(i,j,ms,mt)*x_s(mp,ms,mt)*y_s(mp,ms,mt)                &
                                  + h_st(i,j,ms,mt)*(y_s(mp,ms,mt)*x_t(mp,ms,mt) + y_t(mp,ms,mt)*x_s(mp,ms,mt))                                &
                                  - h_s(i,j,ms,mt)*(x_st(mp,ms,mt)*y_t(mp,ms,mt) - x_tt(mp,ms,mt)*y_s(mp,ms,mt))                               &
                                  - h_t(i,j,ms,mt)*(x_st(mp,ms,mt)*y_s(mp,ms,mt) - x_ss(mp,ms,mt)*y_t(mp,ms,mt)))*element%size(i,j)*vr/xjac**2 &
                                  - xjac_x*(-h_s(i,j,ms,mt)*x_t(mp,ms,mt) + h_t(i,j,ms,mt)*x_s(mp,ms,mt))*element%size(i,j)*vr/xjac**2
            v_px                  = (y_t(mp,ms,mt)*h_s(i,j,ms,mt) - y_s(mp,ms,mt)*h_t(i,j,ms,mt))*element%size(i,j)/xjac
            v_py                  = (-x_t(mp,ms,mt)*h_s(i,j,ms,mt) + x_s(mp,ms,mt)*h_t(i,j,ms,mt))*element%size(i,j)/xjac
            ! Second derivatives wrt phi not implemented in FFT,
            ! and not necessary for stabilization (2nd derivatives only appear in hyperdissipation terms)
!            eq(2*n_var+1,0,0,2,:) = -H(i,j,ms,mt)*element%size(i,j)*vr - x_pp(mp,ms,mt)*eq(2*n_var+1,1,0,0,:) - 2.d0*(x_p(mp,ms,mt)*v_px &
!                                  + y_p(mp,ms,mt)*v_py)*vp - y_pp(mp,ms,mt)*eq(2*n_var+1,0,1,0,:) + 2.d0*(x_p(mp,ms,mt)*x_p_x*eq(2*n_var+1,1,0,0,:) &
!                                  + x_p(mp,ms,mt)*y_p_x*eq(2*n_var+1,0,1,0,:) + y_p(mp,ms,mt)*x_p_y*eq(2*n_var+1,1,0,0,:) &
!                                  + y_p(mp,ms,mt)*y_p_y*eq(2*n_var+1,0,1,0,:)) + x_p(mp,ms,mt)**2*eq(2*n_var+1,2,0,0,:) &
!                                  + 2.d0*x_p(mp,ms,mt)*y_p(mp,ms,mt)*eq(2*n_var+1,1,1,0,:) + y_p(mp,ms,mt)**2*eq(2*n_var+1,0,2,0,:)
            eq(2*n_var+1,1,0,1,:) = v_px*vp - x_p_x*eq(2*n_var+1,1,0,0,:) - x_p(mp,ms,mt)*eq(2*n_var+1,2,0,0,:) - y_p_x*eq(2*n_var+1,0,1,0,:) &
                                  - y_p(mp,ms,mt)*eq(2*n_var+1,1,1,0,:)
            eq(2*n_var+1,0,1,1,:) = v_py*vp - x_p_y*eq(2*n_var+1,1,0,0,:) - y_p(mp,ms,mt)*eq(2*n_var+1,0,2,0,:) - y_p_y*eq(2*n_var+1,0,1,0,:) &
                                  - x_p(mp,ms,mt)*eq(2*n_var+1,1,1,0,:)
            
            
#ifdef DEBUG
            rhs_ij_1 = eval(thread_eq(tid)%rhs1seq)*BigR*xjac
            rhs_ij_2 = eval(thread_eq(tid)%rhs2seq)*BigR*xjac
            rhs_ij_3 = eval(thread_eq(tid)%rhs3seq)*BigR*xjac
            rhs_ij_4 = eval(thread_eq(tid)%rhs4seq)*BigR*xjac
            rhs_ij_5 = eval(thread_eq(tid)%rhs5seq)*BigR*xjac
            rhs_ij_6 = eval(thread_eq(tid)%rhs6seq)*BigR*xjac
#else
#include "rhs_unreadable.h"

            rhs_ij_1 = rhs_ij_1*BigR*xjac
            rhs_ij_2 = rhs_ij_2*BigR*xjac
            rhs_ij_3 = rhs_ij_3*BigR*xjac
            rhs_ij_4 = rhs_ij_4*BigR*xjac
            rhs_ij_5 = rhs_ij_5*BigR*xjac
            rhs_ij_6 = rhs_ij_6*BigR*xjac
#endif

            ij1 = index_ij
            ij2 = index_ij + 1*(n_tor_end - n_tor_start + 1)
            ij3 = index_ij + 2*(n_tor_end - n_tor_start + 1)
            ij4 = index_ij + 3*(n_tor_end - n_tor_start + 1)
            ij5 = index_ij + 4*(n_tor_end - n_tor_start + 1)
            ij6 = index_ij + 5*(n_tor_end - n_tor_start + 1)
            
            ! --- Fill up the matrix
            if (use_fft) then
              RHS_p(mp,ij1) = RHS_p(mp,ij1) + rhs_ij_1(1)*wst
              RHS_p(mp,ij2) = RHS_p(mp,ij2) + rhs_ij_2(1)*wst
              RHS_p(mp,ij3) = RHS_p(mp,ij3) + rhs_ij_3(1)*wst
              RHS_p(mp,ij4) = RHS_p(mp,ij4) + rhs_ij_4(1)*wst
              RHS_p(mp,ij5) = RHS_p(mp,ij5) + rhs_ij_5(1)*wst
              RHS_p(mp,ij6) = RHS_p(mp,ij6) + rhs_ij_6(1)*wst

              RHS_k(mp,ij1) = RHS_k(mp,ij1) + rhs_ij_1(2)*wst
              RHS_k(mp,ij2) = RHS_k(mp,ij2) + rhs_ij_2(2)*wst
              RHS_k(mp,ij3) = RHS_k(mp,ij3) + rhs_ij_3(2)*wst
              RHS_k(mp,ij4) = RHS_k(mp,ij4) + rhs_ij_4(2)*wst
              RHS_k(mp,ij5) = RHS_k(mp,ij5) + rhs_ij_5(2)*wst
              RHS_k(mp,ij6) = RHS_k(mp,ij6) + rhs_ij_6(2)*wst
            else
              RHS(ij1) = RHS(ij1) + (rhs_ij_1(1)*HZ(im,mp) + rhs_ij_1(2)*HZ_p(im,mp))*wst
              RHS(ij2) = RHS(ij2) + (rhs_ij_2(1)*HZ(im,mp) + rhs_ij_2(2)*HZ_p(im,mp))*wst
              RHS(ij3) = RHS(ij3) + (rhs_ij_3(1)*HZ(im,mp) + rhs_ij_3(2)*HZ_p(im,mp))*wst
              RHS(ij4) = RHS(ij4) + (rhs_ij_4(1)*HZ(im,mp) + rhs_ij_4(2)*HZ_p(im,mp))*wst
              RHS(ij5) = RHS(ij5) + (rhs_ij_5(1)*HZ(im,mp) + rhs_ij_5(2)*HZ_p(im,mp))*wst
              RHS(ij6) = RHS(ij6) + (rhs_ij_6(1)*HZ(im,mp) + rhs_ij_6(2)*HZ_p(im,mp))*wst
            endif
            
            do k=1,n_vertex_max
              do l=1,n_order+1
                do in = n_tor_start,n_tor_end
                  if (use_fft) then
                    index_kl =       n_var*(n_order+1)*(k-1) +       n_var*(l-1) + 1
                  else
                    index_kl = n_tor_local*n_var*(n_order+1)*(k-1) + n_tor_local*n_var*(l-1) + in - n_tor_start +1
                  endif
                  
                  ! Unknown increments to next time step (delta u^n)
                  eq(2*n_var+2,0,0,0,:) = H(k,l,ms,mt)*element%size(k,l)*dur
                  eq(2*n_var+2,1,0,0,:) = (y_t(mp,ms,mt)*h_s(k,l,ms,mt) - y_s(mp,ms,mt)*h_t(k,l,ms,mt))*element%size(k,l)*dur/xjac
                  eq(2*n_var+2,0,1,0,:) = (-x_t(mp,ms,mt)*h_s(k,l,ms,mt) + x_s(mp,ms,mt)*h_t(k,l,ms,mt))*element%size(k,l)*dur/xjac
                  eq(2*n_var+2,0,0,1,:) = H(k,l,ms,mt)*element%size(k,l)*dup - eq(2*n_var+2,1,0,0,:)*x_p(mp,ms,mt) - eq(2*n_var+2,0,1,0,:)*y_p(mp,ms,mt)
                  eq(2*n_var+2,2,0,0,:) = (h_ss(k,l,ms,mt)*y_t(mp,ms,mt)**2 - 2.d0*h_st(k,l,ms,mt)*y_s(mp,ms,mt)*y_t(mp,ms,mt)                        &
                                        + h_tt(k,l,ms,mt)*y_s(mp,ms,mt)**2                                                                            &
                                        + h_s(k,l,ms,mt)*(y_st(mp,ms,mt)*y_t(mp,ms,mt) - y_tt(mp,ms,mt)*y_s(mp,ms,mt))                                &
                                        + h_t(k,l,ms,mt)*(y_st(mp,ms,mt)*y_s(mp,ms,mt) - y_ss(mp,ms,mt)*y_t(mp,ms,mt)))*element%size(k,l)*dur/xjac**2 &	
                                        - xjac_x*(h_s(k,l,ms,mt)*y_t(mp,ms,mt) - h_t(k,l,ms,mt)*y_s(mp,ms,mt))*element%size(k,l)*dur/xjac**2
                  eq(2*n_var+2,0,2,0,:) = (h_ss(k,l,ms,mt)*x_t(mp,ms,mt)**2 - 2.d0*h_st(k,l,ms,mt)*x_s(mp,ms,mt)*x_t(mp,ms,mt)                        &
                                        + h_tt(k,l,ms,mt)*x_s(mp,ms,mt)**2                                                                            &
                                        + h_s(k,l,ms,mt)*(x_st(mp,ms,mt)*x_t(mp,ms,mt) - x_tt(mp,ms,mt)*x_s(mp,ms,mt))                                &
                                        + h_t(k,l,ms,mt)*(x_st(mp,ms,mt)*x_s(mp,ms,mt) - x_ss(mp,ms,mt)*x_t(mp,ms,mt)))*element%size(k,l)*dur/xjac**2 &
                                        - xjac_y*(-h_s(k,l,ms,mt)*x_t(mp,ms,mt) + h_t(k,l,ms,mt)*x_s(mp,ms,mt))*element%size(k,l)*dur/xjac**2
                  eq(2*n_var+2,1,1,0,:) = (-h_ss(k,l,ms,mt)*y_t(mp,ms,mt)*x_t(mp,ms,mt) - h_tt(k,l,ms,mt)*x_s(mp,ms,mt)*y_s(mp,ms,mt)                 &
     	                                + h_st(k,l,ms,mt)*(y_s(mp,ms,mt)*x_t(mp,ms,mt)  + y_t(mp,ms,mt)*x_s(mp,ms,mt))                                &
                                        - h_s(k,l,ms,mt)*(x_st(mp,ms,mt)*y_t(mp,ms,mt) - x_tt(mp,ms,mt)*y_s(mp,ms,mt))                                &
                                        - h_t(k,l,ms,mt)*(x_st(mp,ms,mt)*y_s(mp,ms,mt) - x_ss(mp,ms,mt)*y_t(mp,ms,mt)))*element%size(k,l)*dur/xjac**2 &
                                        - xjac_x*(-h_s(k,l,ms,mt)*x_t(mp,ms,mt) + h_t(k,l,ms,mt)*x_s(mp,ms,mt))*element%size(k,l)*dur/xjac**2
                  u_px                  = (y_t(mp,ms,mt)*h_s(k,l,ms,mt) - y_s(mp,ms,mt)*h_t(k,l,ms,mt))*element%size(k,l)/xjac
                  u_py                  = (-x_t(mp,ms,mt)*h_s(k,l,ms,mt) + x_s(mp,ms,mt)*h_t(k,l,ms,mt))*element%size(k,l)/xjac
                  ! Second derivatives wrt phi not implemented in FFT, 
                  ! and not necessary for stabilization (2nd derivatives only appear in hyperdissipation terms)
!                  eq(2*n_var+2,0,0,2,:) = -H(k,l,ms,mt)*element%size(k,l)*dur - x_pp(mp,ms,mt)*eq(2*n_var+2,1,0,0,:) - 2.d0*(x_p(mp,ms,mt)*u_px &
!                                        + y_p(mp,ms,mt)*u_py)*dup - y_pp(mp,ms,mt)*eq(2*n_var+2,0,1,0,:) + 2.d0*(x_p(mp,ms,mt)*x_p_x*eq(2*n_var+2,1,0,0,:) &
!                                        + x_p(mp,ms,mt)*y_p_x*eq(2*n_var+2,0,1,0,:) + y_p(mp,ms,mt)*x_p_y*eq(2*n_var+2,1,0,0,:) &
!                                        + y_p(mp,ms,mt)*y_p_y*eq(2*n_var+2,0,1,0,:)) + x_p(mp,ms,mt)**2*eq(2*n_var+2,2,0,0,:)             &
!                                        + 2.d0*x_p(mp,ms,mt)*y_p(mp,ms,mt)*eq(2*n_var+2,1,1,0,:) + y_p(mp,ms,mt)**2*eq(2*n_var+2,0,2,0,:)
                  eq(2*n_var+2,1,0,1,:) = u_px*dup - x_p_x*eq(2*n_var+2,1,0,0,:) - x_p(mp,ms,mt)*eq(2*n_var+2,2,0,0,:) - y_p_x*eq(2*n_var+2,0,1,0,:) &
                                        - y_p(mp,ms,mt)*eq(2*n_var+2,1,1,0,:)
                  eq(2*n_var+2,0,1,1,:) = u_py*dup - x_p_y*eq(2*n_var+2,1,0,0,:) - y_p(mp,ms,mt)*eq(2*n_var+2,0,2,0,:) - y_p_y*eq(2*n_var+2,0,1,0,:) &
                                        - x_p(mp,ms,mt)*eq(2*n_var+2,1,1,0,:)
                  
                  
#ifdef DEBUG
!---------------------------------------------------------------- equation 1
                  amat_11 = eval(thread_eq(tid)%amat11seq)*BigR*xjac/F0
                  amat_12 = eval(thread_eq(tid)%amat12seq)*BigR*xjac
                  amat_13 = eval(thread_eq(tid)%amat13seq)*BigR*xjac/F0
                  amat_16 = eval(thread_eq(tid)%amat16seq)*BigR*xjac

!---------------------------------------------------------------- equation 2
                  amat_21 = eval(thread_eq(tid)%amat21seq)*BigR*xjac/F0
                  amat_22 = eval(thread_eq(tid)%amat22seq)*BigR*xjac
                  amat_23 = eval(thread_eq(tid)%amat23seq)*BigR*xjac/F0
                  amat_24 = eval(thread_eq(tid)%amat24seq)*BigR*xjac
                  amat_25 = eval(thread_eq(tid)%amat25seq)*BigR*xjac
                  amat_26 = eval(thread_eq(tid)%amat26seq)*BigR*xjac

!---------------------------------------------------------------- equation 3
                  amat_31 = eval(thread_eq(tid)%amat31seq)*BigR*xjac/F0
                  amat_33 = eval(thread_eq(tid)%amat33seq)*BigR*xjac/F0

!---------------------------------------------------------------- equation 4
                  amat_42 = eval(thread_eq(tid)%amat42seq)*BigR*xjac
                  amat_44 = eval(thread_eq(tid)%amat44seq)*BigR*xjac
                 
!---------------------------------------------------------------- equation 5
                  amat_51 = eval(thread_eq(tid)%amat51seq)*BigR*xjac/F0
                  amat_52 = eval(thread_eq(tid)%amat52seq)*BigR*xjac
                  amat_55 = eval(thread_eq(tid)%amat55seq)*BigR*xjac
                 
!---------------------------------------------------------------- equation 6
                  amat_61 = eval(thread_eq(tid)%amat61seq)*BigR*xjac/F0
                  amat_62 = eval(thread_eq(tid)%amat62seq)*BigR*xjac
                  amat_63 = eval(thread_eq(tid)%amat63seq)*BigR*xjac/F0
                  amat_65 = eval(thread_eq(tid)%amat65seq)*BigR*xjac
                  amat_66 = eval(thread_eq(tid)%amat66seq)*BigR*xjac
#else
#include "amat_unreadable.h"

                  amat_11 = amat_11*BigR*xjac/F0; amat_12 = amat_12*BigR*xjac; amat_13 = amat_13*BigR*xjac/F0; amat_16 = amat_16*BigR*xjac
                  amat_21 = amat_21*BigR*xjac/F0; amat_22 = amat_22*BigR*xjac; amat_23 = amat_23*BigR*xjac/F0; amat_24 = amat_24*BigR*xjac
                  amat_25 = amat_25*BigR*xjac; amat_26 = amat_26*BigR*xjac
                  amat_31 = amat_31*BigR*xjac/F0; amat_33 = amat_33*BigR*xjac/F0
                  amat_42 = amat_42*BigR*xjac; amat_44 = amat_44*BigR*xjac
                  amat_51 = amat_51*BigR*xjac/F0; amat_52 = amat_52*BigR*xjac; amat_55 = amat_55*BigR*xjac
                  amat_61 = amat_61*BigR*xjac/F0; amat_62 = amat_62*BigR*xjac; amat_63 = amat_63*BigR*xjac/F0; amat_65 = amat_65*BigR*xjac
                  amat_66 = amat_66*BigR*xjac
#endif

                  kl1 = index_kl
                  kl2 = index_kl + 1*(n_tor_end - n_tor_start + 1)
                  kl3 = index_kl + 2*(n_tor_end - n_tor_start + 1)
                  kl4 = index_kl + 3*(n_tor_end - n_tor_start + 1)
                  kl5 = index_kl + 4*(n_tor_end - n_tor_start + 1)
                  kl6 = index_kl + 5*(n_tor_end - n_tor_start + 1)
                  
                  if (use_fft) then
                    ELM_p(mp,ij1,kl1)  =  ELM_p(mp,ij1,kl1) + wst*amat_11(1)
                    ELM_p(mp,ij1,kl2)  =  ELM_p(mp,ij1,kl2) + wst*amat_12(1)
                    ELM_p(mp,ij1,kl3)  =  ELM_p(mp,ij1,kl3) + wst*amat_13(1)
                    ELM_p(mp,ij1,kl6)  =  ELM_p(mp,ij1,kl6) + wst*amat_16(1)
                    
                    ELM_k(mp,ij1,kl1)  =  ELM_k(mp,ij1,kl1) + wst*amat_11(2)
                    ELM_k(mp,ij1,kl2)  =  ELM_k(mp,ij1,kl2) + wst*amat_12(2)
                    ELM_k(mp,ij1,kl3)  =  ELM_k(mp,ij1,kl3) + wst*amat_13(2)
                    ELM_k(mp,ij1,kl6)  =  ELM_k(mp,ij1,kl6) + wst*amat_16(2)
                    
                    ELM_n(mp,ij1,kl1)  =  ELM_n(mp,ij1,kl1) + wst*amat_11(3)
                    ELM_n(mp,ij1,kl2)  =  ELM_n(mp,ij1,kl2) + wst*amat_12(3)
                    ELM_n(mp,ij1,kl3)  =  ELM_n(mp,ij1,kl3) + wst*amat_13(3)
                    ELM_n(mp,ij1,kl6)  =  ELM_n(mp,ij1,kl6) + wst*amat_16(3)
                    
                    ELM_kn(mp,ij1,kl1)  =  ELM_kn(mp,ij1,kl1) + wst*amat_11(4)
                    ELM_kn(mp,ij1,kl2)  =  ELM_kn(mp,ij1,kl2) + wst*amat_12(4)
                    ELM_kn(mp,ij1,kl3)  =  ELM_kn(mp,ij1,kl3) + wst*amat_13(4)
                    ELM_kn(mp,ij1,kl6)  =  ELM_kn(mp,ij1,kl6) + wst*amat_16(4)


                    ELM_p(mp,ij2,kl1)  =  ELM_p(mp,ij2,kl1) + wst*amat_21(1)
                    ELM_p(mp,ij2,kl2)  =  ELM_p(mp,ij2,kl2) + wst*amat_22(1)
                    ELM_p(mp,ij2,kl3)  =  ELM_p(mp,ij2,kl3) + wst*amat_23(1)
                    ELM_p(mp,ij2,kl4)  =  ELM_p(mp,ij2,kl4) + wst*amat_24(1)
                    ELM_p(mp,ij2,kl5)  =  ELM_p(mp,ij2,kl5) + wst*amat_25(1)
                    ELM_p(mp,ij2,kl6)  =  ELM_p(mp,ij2,kl6) + wst*amat_26(1)
                    
                    ELM_k(mp,ij2,kl1)  =  ELM_k(mp,ij2,kl1) + wst*amat_21(2)
                    ELM_k(mp,ij2,kl2)  =  ELM_k(mp,ij2,kl2) + wst*amat_22(2)
                    ELM_k(mp,ij2,kl3)  =  ELM_k(mp,ij2,kl3) + wst*amat_23(2)
                    ELM_k(mp,ij2,kl4)  =  ELM_k(mp,ij2,kl4) + wst*amat_24(2)
                    ELM_k(mp,ij2,kl5)  =  ELM_k(mp,ij2,kl5) + wst*amat_25(2)
                    ELM_k(mp,ij2,kl6)  =  ELM_k(mp,ij2,kl6) + wst*amat_26(2)
                    
                    ELM_n(mp,ij2,kl1)  =  ELM_n(mp,ij2,kl1) + wst*amat_21(3)
                    ELM_n(mp,ij2,kl2)  =  ELM_n(mp,ij2,kl2) + wst*amat_22(3)
                    ELM_n(mp,ij2,kl3)  =  ELM_n(mp,ij2,kl3) + wst*amat_23(3)
                    ELM_n(mp,ij2,kl4)  =  ELM_n(mp,ij2,kl4) + wst*amat_24(3)
                    ELM_n(mp,ij2,kl5)  =  ELM_n(mp,ij2,kl5) + wst*amat_25(3)
                    ELM_n(mp,ij2,kl6)  =  ELM_n(mp,ij2,kl6) + wst*amat_26(3)
                    
                    ELM_kn(mp,ij2,kl1)  =  ELM_kn(mp,ij2,kl1) + wst*amat_21(4)
                    ELM_kn(mp,ij2,kl2)  =  ELM_kn(mp,ij2,kl2) + wst*amat_22(4)
                    ELM_kn(mp,ij2,kl3)  =  ELM_kn(mp,ij2,kl3) + wst*amat_23(4)
                    ELM_kn(mp,ij2,kl4)  =  ELM_kn(mp,ij2,kl4) + wst*amat_24(4)
                    ELM_kn(mp,ij2,kl5)  =  ELM_kn(mp,ij2,kl5) + wst*amat_25(4)
                    ELM_kn(mp,ij2,kl6)  =  ELM_kn(mp,ij2,kl6) + wst*amat_26(4)


                    ELM_p(mp,ij3,kl1)  =  ELM_p(mp,ij3,kl1) + wst*amat_31(1)
                    ELM_p(mp,ij3,kl3)  =  ELM_p(mp,ij3,kl3) + wst*amat_33(1)
                    
                    ELM_k(mp,ij3,kl1)  =  ELM_k(mp,ij3,kl1) + wst*amat_31(2)
                    ELM_k(mp,ij3,kl3)  =  ELM_k(mp,ij3,kl3) + wst*amat_33(2)
                    
                    ELM_n(mp,ij3,kl1)  =  ELM_n(mp,ij3,kl1) + wst*amat_31(3)
                    ELM_n(mp,ij3,kl3)  =  ELM_n(mp,ij3,kl3) + wst*amat_33(3)
                    
                    ELM_kn(mp,ij3,kl1)  =  ELM_kn(mp,ij3,kl1) + wst*amat_31(4)
                    ELM_kn(mp,ij3,kl3)  =  ELM_kn(mp,ij3,kl3) + wst*amat_33(4)


                    ELM_p(mp,ij4,kl2)  =  ELM_p(mp,ij4,kl2) + wst*amat_42(1)
                    ELM_p(mp,ij4,kl4)  =  ELM_p(mp,ij4,kl4) + wst*amat_44(1)
                    
                    ELM_k(mp,ij4,kl2)  =  ELM_k(mp,ij4,kl2) + wst*amat_42(2)
                    ELM_k(mp,ij4,kl4)  =  ELM_k(mp,ij4,kl4) + wst*amat_44(2)
                    
                    ELM_n(mp,ij4,kl2)  =  ELM_n(mp,ij4,kl2) + wst*amat_42(3)
                    ELM_n(mp,ij4,kl4)  =  ELM_n(mp,ij4,kl4) + wst*amat_44(3)
                    
                    ELM_kn(mp,ij4,kl2)  =  ELM_kn(mp,ij4,kl2) + wst*amat_42(4)
                    ELM_kn(mp,ij4,kl4)  =  ELM_kn(mp,ij4,kl4) + wst*amat_44(4)

 
                    ELM_p(mp,ij5,kl1)  =  ELM_p(mp,ij5,kl1)  + wst*amat_51(1)
                    ELM_p(mp,ij5,kl2)  =  ELM_p(mp,ij5,kl2)  + wst*amat_52(1)
                    ELM_p(mp,ij5,kl5)  =  ELM_p(mp,ij5,kl5)  + wst*amat_55(1)
                    
                    ELM_k(mp,ij5,kl1)  =  ELM_k(mp,ij5,kl1)  + wst*amat_51(2)
                    ELM_k(mp,ij5,kl2)  =  ELM_k(mp,ij5,kl2)  + wst*amat_52(2)
                    ELM_k(mp,ij5,kl5)  =  ELM_k(mp,ij5,kl5)  + wst*amat_55(2)
                    
                    ELM_n(mp,ij5,kl1)  =  ELM_n(mp,ij5,kl1)  + wst*amat_51(3)
                    ELM_n(mp,ij5,kl2)  =  ELM_n(mp,ij5,kl2)  + wst*amat_52(3)
                    ELM_n(mp,ij5,kl5)  =  ELM_n(mp,ij5,kl5)  + wst*amat_55(3)
                    
                    ELM_kn(mp,ij5,kl1)  =  ELM_kn(mp,ij5,kl1)  + wst*amat_51(4)
                    ELM_kn(mp,ij5,kl2)  =  ELM_kn(mp,ij5,kl2)  + wst*amat_52(4)
                    ELM_kn(mp,ij5,kl5)  =  ELM_kn(mp,ij5,kl5)  + wst*amat_55(4)


                    ELM_p(mp,ij6,kl1)  =  ELM_p(mp,ij6,kl1)  + wst*amat_61(1)
                    ELM_p(mp,ij6,kl2)  =  ELM_p(mp,ij6,kl2)  + wst*amat_62(1)
                    ELM_p(mp,ij6,kl3)  =  ELM_p(mp,ij6,kl3)  + wst*amat_63(1)
                    ELM_p(mp,ij6,kl5)  =  ELM_p(mp,ij6,kl5)  + wst*amat_65(1)
                    ELM_p(mp,ij6,kl6)  =  ELM_p(mp,ij6,kl6)  + wst*amat_66(1)
                    
                    ELM_k(mp,ij6,kl1)  =  ELM_k(mp,ij6,kl1)  + wst*amat_61(2)
                    ELM_k(mp,ij6,kl2)  =  ELM_k(mp,ij6,kl2)  + wst*amat_62(2)
                    ELM_k(mp,ij6,kl3)  =  ELM_k(mp,ij6,kl3)  + wst*amat_63(2)
                    ELM_k(mp,ij6,kl5)  =  ELM_k(mp,ij6,kl5)  + wst*amat_65(2)
                    ELM_k(mp,ij6,kl6)  =  ELM_k(mp,ij6,kl6)  + wst*amat_66(2)
                    
                    ELM_n(mp,ij6,kl1)  =  ELM_n(mp,ij6,kl1)  + wst*amat_61(3)
                    ELM_n(mp,ij6,kl2)  =  ELM_n(mp,ij6,kl2)  + wst*amat_62(3)
                    ELM_n(mp,ij6,kl3)  =  ELM_n(mp,ij6,kl3)  + wst*amat_63(3)
                    ELM_n(mp,ij6,kl5)  =  ELM_n(mp,ij6,kl5)  + wst*amat_65(3)
                    ELM_n(mp,ij6,kl6)  =  ELM_n(mp,ij6,kl6)  + wst*amat_66(3)
                    
                    ELM_kn(mp,ij6,kl1)  =  ELM_kn(mp,ij6,kl1)  + wst*amat_61(4)
                    ELM_kn(mp,ij6,kl2)  =  ELM_kn(mp,ij6,kl2)  + wst*amat_62(4)
                    ELM_kn(mp,ij6,kl3)  =  ELM_kn(mp,ij6,kl3)  + wst*amat_63(4)
                    ELM_kn(mp,ij6,kl5)  =  ELM_kn(mp,ij6,kl5)  + wst*amat_65(4)
                    ELM_kn(mp,ij6,kl6)  =  ELM_kn(mp,ij6,kl6)  + wst*amat_66(4)
                  else
                    ELM(ij1,kl1) = ELM(ij1,kl1) + (amat_11(1)*HZ(im,mp)*HZ(in,mp)   + amat_11(2)*HZ_p(im,mp)*HZ(in,mp) &
                                                +  amat_11(3)*HZ(im,mp)*HZ_p(in,mp) + amat_11(4)*HZ_p(im,mp)*HZ_p(in,mp))*wst
                    ELM(ij1,kl2) = ELM(ij1,kl2) + (amat_12(1)*HZ(im,mp)*HZ(in,mp)   + amat_12(2)*HZ_p(im,mp)*HZ(in,mp) &
                                                +  amat_12(3)*HZ(im,mp)*HZ_p(in,mp) + amat_12(4)*HZ_p(im,mp)*HZ_p(in,mp))*wst
                    ELM(ij1,kl3) = ELM(ij1,kl3) + (amat_13(1)*HZ(im,mp)*HZ(in,mp)   + amat_13(2)*HZ_p(im,mp)*HZ(in,mp) &
                                                +  amat_13(3)*HZ(im,mp)*HZ_p(in,mp) + amat_13(4)*HZ_p(im,mp)*HZ_p(in,mp))*wst
                    ELM(ij1,kl6) = ELM(ij1,kl6) + (amat_16(1)*HZ(im,mp)*HZ(in,mp)   + amat_16(2)*HZ_p(im,mp)*HZ(in,mp) &
                                                +  amat_16(3)*HZ(im,mp)*HZ_p(in,mp) + amat_16(4)*HZ_p(im,mp)*HZ_p(in,mp))*wst

                    ELM(ij2,kl1) = ELM(ij2,kl1) + (amat_21(1)*HZ(im,mp)*HZ(in,mp)   + amat_21(2)*HZ_p(im,mp)*HZ(in,mp) &
                                                +  amat_21(3)*HZ(im,mp)*HZ_p(in,mp) + amat_21(4)*HZ_p(im,mp)*HZ_p(in,mp))*wst
                    ELM(ij2,kl2) = ELM(ij2,kl2) + (amat_22(1)*HZ(im,mp)*HZ(in,mp)   + amat_22(2)*HZ_p(im,mp)*HZ(in,mp) &
                                                +  amat_22(3)*HZ(im,mp)*HZ_p(in,mp) + amat_22(4)*HZ_p(im,mp)*HZ_p(in,mp))*wst
                    ELM(ij2,kl3) = ELM(ij2,kl3) + (amat_23(1)*HZ(im,mp)*HZ(in,mp)   + amat_23(2)*HZ_p(im,mp)*HZ(in,mp) &
                                                +  amat_23(3)*HZ(im,mp)*HZ_p(in,mp) + amat_23(4)*HZ_p(im,mp)*HZ_p(in,mp))*wst
                    ELM(ij2,kl4) = ELM(ij2,kl4) + (amat_24(1)*HZ(im,mp)*HZ(in,mp)   + amat_24(2)*HZ_p(im,mp)*HZ(in,mp) &
                                                +  amat_24(3)*HZ(im,mp)*HZ_p(in,mp) + amat_24(4)*HZ_p(im,mp)*HZ_p(in,mp))*wst
                    ELM(ij2,kl5) = ELM(ij2,kl5) + (amat_25(1)*HZ(im,mp)*HZ(in,mp)   + amat_25(2)*HZ_p(im,mp)*HZ(in,mp) &
                                                +  amat_25(3)*HZ(im,mp)*HZ_p(in,mp) + amat_25(4)*HZ_p(im,mp)*HZ_p(in,mp))*wst
                    ELM(ij2,kl6) = ELM(ij2,kl6) + (amat_26(1)*HZ(im,mp)*HZ(in,mp)   + amat_26(2)*HZ_p(im,mp)*HZ(in,mp) &
                                                +  amat_26(3)*HZ(im,mp)*HZ_p(in,mp) + amat_26(4)*HZ_p(im,mp)*HZ_p(in,mp))*wst

                    ELM(ij3,kl1) = ELM(ij3,kl1) + (amat_31(1)*HZ(im,mp)*HZ(in,mp)   + amat_31(2)*HZ_p(im,mp)*HZ(in,mp) &
                                                +  amat_31(3)*HZ(im,mp)*HZ_p(in,mp) + amat_31(4)*HZ_p(im,mp)*HZ_p(in,mp))*wst
                    ELM(ij3,kl3) = ELM(ij3,kl3) + (amat_33(1)*HZ(im,mp)*HZ(in,mp)   + amat_33(2)*HZ_p(im,mp)*HZ(in,mp) &
                                                +  amat_33(3)*HZ(im,mp)*HZ_p(in,mp) + amat_33(4)*HZ_p(im,mp)*HZ_p(in,mp))*wst

                    ELM(ij4,kl2) = ELM(ij4,kl2) + (amat_42(1)*HZ(im,mp)*HZ(in,mp)   + amat_42(2)*HZ_p(im,mp)*HZ(in,mp) &
                                                +  amat_42(3)*HZ(im,mp)*HZ_p(in,mp) + amat_42(4)*HZ_p(im,mp)*HZ_p(in,mp))*wst
                    ELM(ij4,kl4) = ELM(ij4,kl4) + (amat_44(1)*HZ(im,mp)*HZ(in,mp)   + amat_44(2)*HZ_p(im,mp)*HZ(in,mp) &
                                                +  amat_44(3)*HZ(im,mp)*HZ_p(in,mp) + amat_44(4)*HZ_p(im,mp)*HZ_p(in,mp))*wst

                    ELM(ij5,kl1) = ELM(ij5,kl1) + (amat_51(1)*HZ(im,mp)*HZ(in,mp)   + amat_51(2)*HZ_p(im,mp)*HZ(in,mp) &
                                                +  amat_51(3)*HZ(im,mp)*HZ_p(in,mp) + amat_51(4)*HZ_p(im,mp)*HZ_p(in,mp))*wst
                    ELM(ij5,kl2) = ELM(ij5,kl2) + (amat_52(1)*HZ(im,mp)*HZ(in,mp)   + amat_52(2)*HZ_p(im,mp)*HZ(in,mp) &
                                                +  amat_52(3)*HZ(im,mp)*HZ_p(in,mp) + amat_52(4)*HZ_p(im,mp)*HZ_p(in,mp))*wst
                    ELM(ij5,kl5) = ELM(ij5,kl5) + (amat_55(1)*HZ(im,mp)*HZ(in,mp)   + amat_55(2)*HZ_p(im,mp)*HZ(in,mp) &
                                                +  amat_55(3)*HZ(im,mp)*HZ_p(in,mp) + amat_55(4)*HZ_p(im,mp)*HZ_p(in,mp))*wst

                    ELM(ij6,kl1) = ELM(ij6,kl1) + (amat_61(1)*HZ(im,mp)*HZ(in,mp)   + amat_61(2)*HZ_p(im,mp)*HZ(in,mp) &
                                                +  amat_61(3)*HZ(im,mp)*HZ_p(in,mp) + amat_61(4)*HZ_p(im,mp)*HZ_p(in,mp))*wst
                    ELM(ij6,kl2) = ELM(ij6,kl2) + (amat_62(1)*HZ(im,mp)*HZ(in,mp)   + amat_62(2)*HZ_p(im,mp)*HZ(in,mp) &
                                                +  amat_62(3)*HZ(im,mp)*HZ_p(in,mp) + amat_62(4)*HZ_p(im,mp)*HZ_p(in,mp))*wst
                    ELM(ij6,kl3) = ELM(ij6,kl3) + (amat_63(1)*HZ(im,mp)*HZ(in,mp)   + amat_63(2)*HZ_p(im,mp)*HZ(in,mp) &
                                                +  amat_63(3)*HZ(im,mp)*HZ_p(in,mp) + amat_63(4)*HZ_p(im,mp)*HZ_p(in,mp))*wst
                    ELM(ij6,kl5) = ELM(ij6,kl5) + (amat_65(1)*HZ(im,mp)*HZ(in,mp)   + amat_65(2)*HZ_p(im,mp)*HZ(in,mp) &
                                                +  amat_65(3)*HZ(im,mp)*HZ_p(in,mp) + amat_65(4)*HZ_p(im,mp)*HZ_p(in,mp))*wst
                    ELM(ij6,kl6) = ELM(ij6,kl6) + (amat_66(1)*HZ(im,mp)*HZ(in,mp)   + amat_66(2)*HZ_p(im,mp)*HZ(in,mp) &
                                                +  amat_66(3)*HZ(im,mp)*HZ_p(in,mp) + amat_66(4)*HZ_p(im,mp)*HZ_p(in,mp))*wst
                  end if
                end do ! in loop (n_tor, or not...)
              end do ! l loop (n_order+1)
            end do ! k loop (n_vertex)
          end do ! im loop (n_tor, or not...)
        end do ! j loop
      end do ! i loop
    end do ! mp loop (n_plane)
  end do ! ms loop
end do ! mt loop


if (use_fft) then
  do i=1,n_vertex_max*n_var*(n_order+1)
    do j=1, n_vertex_max*n_var*(n_order+1)
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

          if ((l .ge. 0) .and. (l .le. n_plane/2)) then
            ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   + real(out_fft(l+1))
            ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   - imag(out_fft(l+1))
            ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) - imag(out_fft(l+1))
            ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) - real(out_fft(l+1))
          else if ((l .lt. 0) .and. (abs(l) .le. n_plane/2)) then
            ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   + real(out_fft(abs(l)+1))
            ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   + imag(out_fft(abs(l)+1))
            ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) + imag(out_fft(abs(l)+1))
            ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) - real(out_fft(abs(l)+1))
          end if

          l = (k-1) - (m-1)

          if ((l .ge. 0) .and. (l .le. n_plane/2)) then
            ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   + real(out_fft(l+1))
            ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   - imag(out_fft(l+1))
            ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) + imag(out_fft(l+1))
            ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) + real(out_fft(l+1))
          else if ((l .lt. 0) .and. (abs(l) .le. n_plane/2)) then
            ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   + real(out_fft(abs(l)+1))
            ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   + imag(out_fft(abs(l)+1))
            ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) - imag(out_fft(abs(l)+1))
            ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) + real(out_fft(abs(l)+1))
          end if
        end do
      end do
    end do
  end do

  do i=1,n_vertex_max*n_var*(n_order+1)
    do j=1, n_vertex_max*n_var*(n_order+1)
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

            if ((l .ge. 0) .and. (l .le. n_plane/2)) then
              ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   + imag(out_fft(l+1)) * float(mode(im))
              ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   + real(out_fft(l+1)) * float(mode(im))
              ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) + real(out_fft(l+1)) * float(mode(im))
              ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) - imag(out_fft(l+1)) * float(mode(im))
            else if ((l .lt. 0) .and. (abs(l) .le. n_plane/2)) then
              ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   - imag(out_fft(abs(l)+1)) * float(mode(im))
              ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   + real(out_fft(abs(l)+1)) * float(mode(im))
              ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) + real(out_fft(abs(l)+1)) * float(mode(im))
              ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) + imag(out_fft(abs(l)+1)) * float(mode(im))
            end if

            l = (k-1) - (m-1)

            if ((l .ge. 0) .and. (l .le. n_plane/2)) then
              ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   - imag(out_fft(l+1)) * float(mode(im))
              ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   - real(out_fft(l+1)) * float(mode(im))
              ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) + real(out_fft(l+1)) * float(mode(im))
              ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) - imag(out_fft(l+1)) * float(mode(im))
            else if ((l .lt. 0) .and. (abs(l) .le. n_plane/2)) then
              ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   + imag(out_fft(abs(l)+1)) * float(mode(im))
              ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   - real(out_fft(abs(l)+1)) * float(mode(im))
              ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) + real(out_fft(abs(l)+1)) * float(mode(im))
              ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) + imag(out_fft(abs(l)+1)) * float(mode(im))
            end if
          end do
        end do
      end if
    end do
  end do

  do i=1,n_vertex_max*n_var*(n_order+1)
    do j=1, n_vertex_max*n_var*(n_order+1)
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

            if ((l .ge. 0) .and. (l .le. n_plane/2)) then
              ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   + imag(out_fft(l+1)) * float(mode(ik))
              ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   + real(out_fft(l+1)) * float(mode(ik))
              ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) + real(out_fft(l+1)) * float(mode(ik))
              ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) - imag(out_fft(l+1)) * float(mode(ik))
            else if ((l .lt. 0) .and. (abs(l) .le. n_plane/2)) then
              ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   - imag(out_fft(abs(l)+1)) * float(mode(ik))
              ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   + real(out_fft(abs(l)+1)) * float(mode(ik))
              ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) + real(out_fft(abs(l)+1)) * float(mode(ik))
              ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) + imag(out_fft(abs(l)+1)) * float(mode(ik))
            end if

            l = (k-1) - (m-1)

            if ((l .ge. 0) .and. (l .le. n_plane/2)) then
              ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   + imag(out_fft(l+1)) * float(mode(ik))
              ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   + real(out_fft(l+1)) * float(mode(ik))
              ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) - real(out_fft(l+1)) * float(mode(ik))
              ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) + imag(out_fft(l+1)) * float(mode(ik))
            else if ((l .lt. 0) .and. (abs(l) .le. n_plane/2)) then
              ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   - imag(out_fft(abs(l)+1)) * float(mode(ik))
              ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   + real(out_fft(abs(l)+1)) * float(mode(ik))
              ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) - real(out_fft(abs(l)+1)) * float(mode(ik))
              ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) - imag(out_fft(abs(l)+1)) * float(mode(ik))
            end if
          end do
        end do
      end if
    end do
  end do

  do i=1,n_vertex_max*n_var*(n_order+1)
    do j=1, n_vertex_max*n_var*(n_order+1)
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

            if ((l .ge. 0) .and. (l .le. n_plane/2)) then
              ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   - real(out_fft(l+1)) * float(mode(im)) * float(mode(ik))
              ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   + imag(out_fft(l+1)) * float(mode(im)) * float(mode(ik))
              ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) + imag(out_fft(l+1)) * float(mode(im)) * float(mode(ik))
              ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) + real(out_fft(l+1)) * float(mode(im)) * float(mode(ik))
            else if ((l .lt. 0) .and. (abs(l) .le. n_plane/2)) then
              ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   - real(out_fft(abs(l)+1)) * float(mode(im)) * float(mode(ik))
              ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   - imag(out_fft(abs(l)+1)) * float(mode(im)) * float(mode(ik))
              ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) - imag(out_fft(abs(l)+1)) * float(mode(im)) * float(mode(ik))
              ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) + real(out_fft(abs(l)+1)) * float(mode(im)) * float(mode(ik))
            end if

            l = (k-1) - (m-1)

            if ((l .ge. 0) .and. (l .le. n_plane/2)) then
              ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   + real(out_fft(l+1)) * float(mode(im)) * float(mode(ik))
              ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   - imag(out_fft(l+1)) * float(mode(im)) * float(mode(ik))
              ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) + imag(out_fft(l+1)) * float(mode(im)) * float(mode(ik))
              ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) + real(out_fft(l+1)) * float(mode(im)) * float(mode(ik))
            else if ((l .lt. 0) .and. (abs(l) .le. n_plane/2)) then
              ELM(index_k,  index_m  ) = ELM(index_k,  index_m)   + real(out_fft(abs(l)+1)) * float(mode(im)) * float(mode(ik))
              ELM(index_k+1,index_m  ) = ELM(index_k+1,index_m)   + imag(out_fft(abs(l)+1)) * float(mode(im)) * float(mode(ik))
              ELM(index_k,  index_m+1) = ELM(index_k,  index_m+1) - imag(out_fft(abs(l)+1)) * float(mode(im)) * float(mode(ik))
              ELM(index_k+1,index_m+1) = ELM(index_k+1,index_m+1) + real(out_fft(abs(l)+1)) * float(mode(im)) * float(mode(ik))
            end if
          end do
        end do
      end if
    end do
  end do

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
    end do
  end do

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
    end do
  end do
else
  return
end if

return
end subroutine element_matrix_fft


subroutine my_fft(in_fft,out_fft,n)

implicit none
      
real*8     :: in_fft(*)
complex*16 :: out_fft(*)
integer    :: n
      
real*8     :: tmp_fft(2*n+2)
integer    :: i
      
tmp_fft(1:n) = in_fft(1:n)
      
call RFT2(tmp_fft,n,1)
      
do i=1,n
  out_fft(i) = cmplx(tmp_fft(2*i-1),tmp_fft(2*i))
enddo
      
return
end subroutine my_fft
end module mod_elt_matrix_fft
