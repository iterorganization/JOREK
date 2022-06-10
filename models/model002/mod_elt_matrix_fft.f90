module mod_elt_matrix_fft

  implicit none

contains

#include "corr_neg_include.f90"

subroutine element_matrix_fft(element, nodes, xpoint2, xcase2, R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint, Z_xpoint, ELM, RHS, tid, &
                              ELM_p, ELM_n, ELM_k, ELM_kn, RHS_p, RHS_k,  eq_g, eq_s, eq_t, eq_p, eq_ss, eq_st, eq_tt, delta_g, delta_s, delta_t, &
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
use corr_neg
use mod_sources

implicit none

type (type_element)       :: element
type (type_node)          :: nodes(n_vertex_max)
type (type_node),optional :: aux_nodes(n_vertex_max)

#define DIM0 n_tor*n_vertex_max*n_degrees*n_var

integer, intent(in)            :: tid
integer, intent(in)            :: i_tor_min, i_tor_max

real*8, dimension (DIM0,DIM0)  :: ELM
real*8, dimension (DIM0)       :: RHS

integer    :: i, j, ms, mt, mp, k, l, index_ij, index_kl, index, index_k, index_m, m, ik, xcase2
integer    :: n_tor_start, n_tor_end, n_tor_local, n_tor_loop
integer    :: in, im, ij, kl
real*8     :: wst, xjac, xjac_s, xjac_t, xjac_x, xjac_y, BigR, r2, phi, delta_phi
real*8     :: particle_source(n_gauss,n_gauss), heat_source(n_gauss,n_gauss)
real*8     :: R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint(2), Z_xpoint(2)
real*8     :: D_prof, theta, zeta, delta_u_x, delta_u_y
real*8     :: rhs_ij(n_var), rhs_ij_k(n_var)
real*8     :: amat(n_var,n_var), amat_k(n_var,n_var), amat_n(n_var,n_var), amat_kn(n_var,n_var)

real*8     :: v, v_x, v_y, v_s, v_t, v_p, v_ss, v_st, v_tt, v_xx, v_xy, v_yy
real*8     :: u0, u0_x, u0_y, u0_p, u0_s, u0_t, u0_ss, u0_tt, u0_st, u0_xx, u0_xy, u0_yy
real*8     :: w0, w0_x, w0_y, w0_p, w0_s, w0_t, w0_ss, w0_st, w0_tt, w0_xx, w0_xy, w0_yy
real*8     :: r0, r0_x, r0_y, r0_p, r0_s, r0_t, r0_ss, r0_st, r0_tt, r0_xx, r0_xy, r0_yy, r0_hat, r0_x_hat, r0_y_hat, r0_corr
real*8     :: psi, psi_x, psi_y, psi_p, psi_s, psi_t, psi_ss, psi_st, psi_tt, psi_xx, psi_xy, psi_yy
real*8     :: u, u_x, u_y, u_p, u_s, u_t, u_ss, u_st, u_tt, u_xx, u_xy, u_yy
real*8     :: w, w_x, w_y, w_p, w_s, w_t, w_ss, w_st, w_tt, w_xx, w_xy, w_yy
real*8     :: rho, rho_x, rho_y, rho_s, rho_t, rho_p, rho_hat, rho_x_hat, rho_y_hat, rho_ss, rho_st, rho_tt, rho_xx, rho_xy, rho_yy
real*8     :: P0,  P0_s,  P0_t,  P0_x,  P0_y,  P0_p
real*8     :: BigR_x, vv2

real*8     :: dr0_corr_dn 
real*8     :: P0_x_rho, P0_xx_rho, P0_y_rho, P0_yy_rho, P0_xy_rho

logical    :: xpoint2, use_fft
real*8     :: in_fft(1:n_plane)
complex*16 :: out_fft(1:n_plane)
integer*8  :: plan

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

particle_source = 0.d0

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

    ! --- The sources are made with R instead of psi
    call sources(xpoint2, xcase2, y_g(ms,mt), Z_xpoint, x_g(ms,mt),R_begin,R_end, &
                 particle_source(ms,mt),heat_source(ms,mt))

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

          P0    = r0    
          P0_x  = r0_x  
          P0_y  = r0_y  
          P0_s  = r0_s  
          P0_t  = r0_t  
          P0_p  = r0_p  

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

          ! --- Increase diffusivity if very small density/temperature
          D_prof = D_perp(1)
          if (r0 .lt. D_prof_neg_thresh)  then
            D_prof  = D_prof_neg
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


            !###################################################################################################
            !#  Perpendicular Momentum Equation                                                                #
            !###################################################################################################

            rhs_ij(var_u) =  - 0.5d0 * vv2 * (v_x * r0_y_hat - v_y * r0_x_hat)                              * xjac * tstep &
                         - r0_hat * BigR**2 * w0 * (v_s * u0_t - v_t * u0_s)                                       * tstep &
                         - visco   * BigR * (v_x * w0_x + v_y * w0_y)                                       * xjac * tstep &
                         + BigR**2 * (v_s * p0_t - v_t * p0_s)                                                     * tstep &
                         - visco_num   * (v_xx + v_x/Bigr + v_yy)*(w0_xx + w0_x/Bigr + w0_yy) * xjac * tstep &
                         + BigR**3 * particle_source(ms,mt) * (v_x * u0_x + v_y * u0_y) * xjac* tstep &
                         - zeta * BigR * r0_hat * (v_x * delta_u_x + v_y * delta_u_y) * xjac
            

            !###################################################################################################
            !#  Vorticity Definition Equation                                                                  #
            !###################################################################################################

            rhs_ij(var_w) = - ( v_x * u0_x   + v_y * u0_y  + v*w0)  * BigR * xjac 

            !###################################################################################################
            !#  Density Equation                                                                               #
            !###################################################################################################

            rhs_ij(var_rho)  = v * BigR * particle_source(ms,mt)                                  * xjac * tstep &
                       + v * BigR**2 * ( r0_s * u0_t - r0_t * u0_s)                                      * tstep &
                       + v * 2.d0 * BigR * r0 * u0_y                                              * xjac * tstep &
                       - D_prof * BigR  * (v_x*r0_x + v_y*r0_y                                  ) * xjac * tstep &
                       + zeta * v * delta_g(mp,var_rho,ms,mt) * BigR                              * xjac         &
                       - D_perp_num * (v_xx + v_x/Bigr + v_yy)*(r0_xx + r0_x/Bigr + r0_yy) * BigR * xjac * tstep 


            rhs_ij_k(var_rho) = - D_prof * BigR  * (                  v_p*r0_p /BigR**2 )         * xjac * tstep 


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

                  u    = psi    ;  w    = psi    ; rho    = psi    
                  u_x  = psi_x  ;  w_x  = psi_x  ; rho_x  = psi_x  
                  u_y  = psi_y  ;  w_y  = psi_y  ; rho_y  = psi_y  
                  u_p  = psi_p  ;  w_p  = psi_p  ; rho_p  = psi_p  
                  u_s  = psi_s  ;  w_s  = psi_s  ; rho_s  = psi_s  
                  u_t  = psi_t  ;  w_t  = psi_t  ; rho_t  = psi_t  
                  u_ss = psi_ss ;  w_ss = psi_ss ; rho_ss = psi_ss 
                  u_tt = psi_tt ;  w_tt = psi_tt ; rho_tt = psi_tt 
                  u_st = psi_st ;  w_st = psi_st ; rho_st = psi_st 

                  u_xx = psi_xx ;  w_xx = psi_xx ; rho_xx = psi_xx
                  u_yy = psi_yy ;  w_yy = psi_yy ; rho_yy = psi_yy
                  u_xy = psi_xy ;  w_xy = psi_xy ; rho_xy = psi_xy

                  rho_hat   = BigR**2 * rho
                  rho_x_hat = 2.d0 * BigR * BigR_x  * rho + BigR**2 * rho_x
                  rho_y_hat = BigR**2 * rho_y

                  P0_x_rho  = rho_x 
                  P0_xx_rho = rho_xx
                  P0_y_rho  = rho_y 
                  P0_yy_rho = rho_yy
                  P0_xy_rho = rho_xy
                 

                  !###################################################################################################
                  !#  Perpendicular Momentum Equation                                                                #
                  !###################################################################################################

                  amat(var_u,var_u) = - BigR**3 * r0_corr * (v_x * u_x + v_y * u_y) * xjac * (1.d0 + zeta)                                &
                                    + r0_hat * BigR**2 * w0 * (v_s * u_t  - v_t  * u_s)                                   * theta * tstep &
                                    + BigR**2 * (u_x * u0_x + u_y * u0_y) * (v_x * r0_y_hat - v_y * r0_x_hat)      * xjac * theta * tstep &
                                    - BigR**3 * particle_source(ms,mt) * (v_x * u_x + v_y * u_y)                   * xjac * theta * tstep 

                  amat(var_u,var_w) = r0_hat * BigR**2 * w  * ( v_s * u0_t - v_t * u0_s)                                   * theta * tstep &
                                    + BigR * ( v_x * w_x + v_y * w_y) * visco    * xjac                                    * theta * tstep &
                                    + visco_num   * (v_xx + v_x/BigR + v_yy)*(w_xx + w_x/BigR + w_yy)               * xjac * theta * tstep 

                  amat(var_u,var_rho) = + 0.5d0 * vv2 * (v_x * rho_y_hat - v_y * rho_x_hat)                         * xjac * theta * tstep &
                                      + rho_hat * BigR**2 * w0 * (v_s * u0_t   - v_t * u0_s)                               * theta * tstep &
                                      - BigR**2 * (v_s * rho_t - v_t * rho_s )                                             * theta * tstep 
                                                                                                                           

                  !###################################################################################################
                  !#  Vorticity Definition Equation                                                                  #
                  !###################################################################################################

                  amat(var_w,var_w) =  v * w * BigR * xjac                                
                  amat(var_w,var_u) = (v_x * u_x + v_y * u_y) * BigR * xjac              

                  !###################################################################################################
                  !#  Density Equation                                                                               #
                  !###################################################################################################


                  amat(var_rho,var_u) =- v * BigR**2 * ( r0_s * u_t - r0_t * u_s)                                        * theta * tstep &
                                       - v * 2.d0 * BigR * r0 * u_y                                               * xjac * theta * tstep 

                  amat(var_rho,var_rho) = v * rho * BigR * (1.d0 + zeta)                                                * xjac   &
                          - v * BigR**2 * ( rho_s * u0_t - rho_t * u0_s)                                        * theta * tstep &
                          - v * 2.d0 * BigR * rho * u0_y                                                 * xjac * theta * tstep &
                          + D_prof * BigR  * (v_x*rho_x + v_y*rho_y )                                    * xjac * theta * tstep &
                          + D_perp_num * (v_xx + v_x/BigR + v_yy)*(rho_xx + rho_x/BigR + rho_yy)  * BigR * xjac * theta * tstep 

                  amat_kn(var_rho,var_rho) = + D_prof * BigR  * ( v_p*rho_p /BigR**2 )                   * xjac * theta * tstep 






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
