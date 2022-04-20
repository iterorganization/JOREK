module mod_boundary_matrix_open
  implicit none
contains

subroutine boundary_matrix_open(vertex, direction, element, nodes, xpoint2, xcase2, R_axis, Z_axis, psi_axis, &
                                psi_bnd, R_xpoint, Z_xpoint, ELM, RHS, i_tor_min, i_tor_max)
!---------------------------------------------------------------------
! calculates the matrix contribution of the boundaries of one element
! implements the natural boundary conditions
!---------------------------------------------------------------------
use constants
use mod_parameters
use data_structure
use gauss
use basis_at_gaussian
use phys_module
use corr_neg
use mod_interp
use diffusivities, only: get_dperp, get_zkperp

implicit none

type (type_element)   :: element
type (type_node)      :: nodes(n_vertex_max)        ! the two nodes containing the boundary nodes
integer, intent(in)   :: i_tor_min   
integer, intent(in)   :: i_tor_max   

real*8     :: x_g(n_gauss), x_s(n_gauss), x_t(n_gauss), x_ss(n_gauss)
real*8     :: y_g(n_gauss), y_s(n_gauss), y_t(n_gauss), y_ss(n_gauss)

real*8     :: eq_g(n_plane,n_var,n_gauss), eq_s(n_plane,n_var,n_gauss), eq_p(n_plane,n_var,n_gauss)
real*8     :: eq_t(n_plane,n_var,n_gauss), eq_ss(n_plane,n_var,n_gauss)
real*8     :: delta_g(n_plane,n_var,n_gauss), delta_s(n_plane,n_var,n_gauss)

real*8     :: ELM(n_vertex_max*n_var*n_degrees*n_tor,n_vertex_max*n_var*n_degrees*n_tor)
real*8     :: RHS(n_vertex_max*n_var*n_degrees*n_tor)
real*8     :: rhs_ij(n_var), amat(n_var,n_var)

integer    :: vertex(2), direction(2), direction_perp(2), bnd_type1, bnd_type2
integer    :: i, j, j2, j3, ms, mt, mp, k, l, l2, l3, index_ij, index_kl, index, xcase2, is
integer    :: in, im, ij1, ij2, ij3, ij4, ij5, ij6, ij7, ij8, kl1, kl2, kl3, kl4, kl5, kl6, kl7, kl8, i_var, k_var
real*8     :: ws, xjac,  dl, BigR, phi, eps_cyl, Btot
real*8     :: R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint(2), Z_xpoint(2)
real*8     :: rhs_ij_5, rhs_ij_6, rhs_ij_7, rhs_ij_8
real*8     :: theta, zeta, Zbig, BB2, bdotn, factor, psi_ss, vpar_ss
real*8     :: R_inside, Z_inside, R_mid, Z_mid, R_cnt, Z_cnt, normal(2), normal_direction(2)
real*8     :: normal_sign, normal_sign3

real*8     :: v, v_x, v_y, v_s, v_p, v_ss, v_xx, v_yy, v_xs, v_ys
real*8     :: ps0, ps0_s, ps0_t, ps0_x, ps0_y, Vpar0, r0_corr, T0_corr, Ti0_corr, Te0_corr, cs0  
real*8     :: psi, psi_s, psi_t, vpar, T, Ti, Te, cs_T, cs_Ti, cs_Te
real*8     :: T0,   T0_s,  T0_t, T0_p
real*8     :: Ti0, Ti0_s, Ti0_t, Ti0_x, Ti0_y, Ti0_p
real*8     :: Te0, Te0_s, Te0_t, Te0_x, Te0_y, Te0_p
real*8     :: r0, r0_s, r0_t, r0_p, r0_x, r0_y, rho, rho_s, rho_t, rho_x, rho_y
real*8     :: c_1, c_2, c_3, c_angle, neutral_source
real*8     :: element_size_ij, element_size_kl, element_size_perp
real*8     :: grad_t(2), B0_R, B0_Z, factor_cs_bnd_integral
logical    :: xpoint2
integer    :: n_tor_local 
logical    :: apply_natural_bc(0:n_var)

type (type_node)         :: tmp_node


return
end subroutine

end module mod_boundary_matrix_open
