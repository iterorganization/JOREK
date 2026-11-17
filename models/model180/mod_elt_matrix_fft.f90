module mod_elt_matrix_fft
implicit none
contains
subroutine element_matrix_fft(element, nodes, xpoint2, xcase2, R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint, Z_xpoint, ELM, RHS, tid, &
  ELM_p, ELM_n, ELM_k, ELM_kn, RHS_p, RHS_k,  eq_g, eq_s, eq_t, eq_p, eq_ss, eq_st, eq_tt, delta_g, delta_s, delta_t, i_tor_min, i_tor_max, aux_nodes, ELM_pnn)
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
use mod_elt_matrix

implicit none
 
type (type_element), intent(in)           :: element
type (type_node),    intent(in)           :: nodes(n_vertex_max)
type (type_node),    intent(in), optional :: aux_nodes(n_vertex_max)

#define DIM0 n_tor*n_vertex_max*(n_order+1)*n_var

!real*8, dimension (DIM0,DIM0)       :: ELM
!real*8, dimension (DIM0)            :: RHS
real*8, dimension(:,:), allocatable :: ELM
real*8, dimension(:),   allocatable :: RHS
integer                , intent(in) :: tid
integer                , intent(in) :: i_tor_min, i_tor_max

integer    :: i, j, ms, mt, mp, k, l, index_ij, index_kl, index, index_k, index_m, m, ik, xcase2
integer    :: in, im, ij1, ij2, ij3, ij4, ij5, ij6, kl1, kl2, kl3, kl4, kl5, kl6
real*8     :: wst, xjac, xjac_x, xjac_y, xjac_s, xjac_t, BigR, r2, phi, eps_cyl
real*8     :: current_source(n_gauss,n_gauss),particle_source(n_gauss,n_gauss),heat_source(n_gauss,n_gauss)
real*8     :: R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint(2), Z_xpoint(2), dj_dpsi, dj_dz
real*8     :: Bgrad_rho_star,     Bgrad_rho,     Bgrad_T_star,  Bgrad_T, BB2
real*8     :: Bgrad_rho_star_psi, Bgrad_rho_psi, Bgrad_rho_rho, Bgrad_T_star_psi, Bgrad_T_psi, Bgrad_T_T, BB2_psi
real*8     :: Bgrad_rho_rho_n, Bgrad_T_T_n, Bgrad_rho_k_star, Bgrad_T_k_star
real*8     :: D_prof, ZK_prof, psi_norm
real*8     :: rhs_ij_1,   rhs_ij_2,   rhs_ij_3,   rhs_ij_4,   rhs_ij_5,   rhs_ij_6
real*8     :: rhs_ij_5_k, rhs_ij_6_k
real*8     :: rhs_stab_1, rhs_stab_2, rhs_stab_3, rhs_stab_4, rhs_stab_5, rhs_stab_6

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
real*8     :: amat_11, amat_12, amat_21, amat_22, amat_23, amat_24, amat_25, amat_26, amat_33, amat_31, amat_44, amat_42
real*8     :: amat_51, amat_52, amat_55, amat_61, amat_62, amat_66, amat_16, amat_13
real*8     :: amat_12_n, amat_23_n, amat_51_k, amat_55_kn, amat_55_k, amat_55_n, amat_61_k, amat_66_kn, amat_66_k, amat_66_n
real*8     :: amat_stab_11, amat_stab_12, amat_stab_13, amat_stab_14 ,amat_stab_21,amat_stab_22, amat_stab_23, amat_stab_24
real*8     :: amat_stab_31, amat_stab_32, amat_stab_33, amat_stab_34 ,amat_stab_41,amat_stab_42, amat_stab_43, amat_stab_44
real*8     :: theta, zeta, delta_u_x, delta_u_y
logical    :: xpoint2

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
real*8, dimension(DIM1, DIM2, DIM2) :: ELM_pnn

real*8, dimension(n_gauss,n_gauss)    :: x_g, x_s, x_t
real*8, dimension(n_gauss,n_gauss)    :: x_ss, x_st, x_tt
real*8, dimension(n_gauss,n_gauss)    :: y_g, y_s, y_t
real*8, dimension(n_gauss,n_gauss)    :: y_ss, y_st, y_tt

real*8, dimension(n_plane,n_eq_var,n_gauss,n_gauss) :: eq_g, eq_s, eq_t
real*8, dimension(n_plane,n_eq_var,n_gauss,n_gauss) :: eq_p
real*8, dimension(n_plane,n_eq_var,n_gauss,n_gauss) :: eq_ss, eq_st, eq_tt
real*8, dimension(n_plane,n_eq_var,n_gauss,n_gauss) :: delta_g, delta_s, delta_t

! FFT version not implemented yet
call element_matrix(element, nodes, xpoint2, xcase2, R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint, Z_xpoint, ELM, RHS, tid, i_tor_min, i_tor_max)

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
