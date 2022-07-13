module mod_elt_matrix
  implicit none
contains

subroutine element_matrix(element,nodes, xpoint2, xcase2, R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint, Z_xpoint, ELM, RHS, tid, i_tor_min, i_tor_max, aux_nodes)
  !---------------------------------------------------------------
  ! calculates the matrix contribution of one element
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

implicit none

type (type_element), intent(in)           :: element
type (type_node),    intent(in)           :: nodes(n_vertex_max)
type (type_node),    intent(in), optional :: aux_nodes(n_vertex_max)

real*8, dimension (:,:), allocatable  :: ELM
real*8, dimension (:)  , allocatable  :: RHS
integer, intent(in) :: tid, i_tor_min, i_tor_max

integer    :: i, j, ms, mt, mp, k, l, index_ij, index_kl, index, xcase2
integer    :: in, im, ij1, ij2, ij3, ij4, ij5, ij6, kl1, kl2, kl3, kl4, kl5, kl6
real*8     :: wst,  xjac, xjac_x, xjac_y, x_p_x, x_p_y, y_p_x, y_p_y, BigR, phi
real*8     :: R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint(2), Z_xpoint(2)
real*8     :: psi_norm, reta, zeta, theta
real*8     :: v_px, v_py, u_px, u_py

real*8, dimension(4) :: rhs_ij_1, rhs_ij_3
real*8, dimension(4) :: amat_11, amat_13, amat_22, amat_33, amat_44, amat_55, amat_66

logical    :: xpoint2
integer    :: n_tor_local

real*8, dimension(n_plane,n_gauss,n_gauss) :: x_g, x_s, x_t, x_p, x_ss, x_st, x_tt, x_sp, x_tp, x_pp
real*8, dimension(n_plane,n_gauss,n_gauss) :: y_g, y_s, y_t, y_p, y_ss, y_st, y_tt, y_sp, y_tp, y_pp
#ifdef altcs
real*8, dimension(n_gauss,n_gauss) :: psieq, psieq_s, psieq_t
real*8                             :: psieq_x, psieq_y
#endif

real*8, dimension(:,:,:,:) , pointer :: eq_g, eq_s, eq_t
real*8, dimension(:,:,:,:) , pointer :: eq_st, eq_ss, eq_tt
real*8, dimension(:,:,:,:) , pointer :: eq_p, eq_pp, eq_sp, eq_tp

real*8, dimension(:,:,:,:,:), pointer :: eq
real*8, dimension(n_var)            :: eq_px, eq_py

real*8, dimension(n_gauss,n_gauss) :: press_gvec
real*8, dimension(n_dim+1,n_plane,n_gauss,n_gauss) :: B_gvec

eq_g    => thread_struct(tid)%eq_g
eq_s    => thread_struct(tid)%eq_s
eq_t    => thread_struct(tid)%eq_t
eq_p    => thread_struct(tid)%eq_p
eq_pp   => thread_struct(tid)%eq_pp
eq_sp   => thread_struct(tid)%eq_sp
eq_tp   => thread_struct(tid)%eq_tp
eq_ss   => thread_struct(tid)%eq_ss
eq_st   => thread_struct(tid)%eq_st
eq_tt   => thread_struct(tid)%eq_tt

eq => thread_eq(tid)%eq

ELM = 0.d0
RHS = 0.d0

! --- Take time evolution parameters from phys_module
theta = time_evol_theta
! change zeta for variable dt
zeta  = time_evol_zeta * 2.0d0 * tstep / (tstep + tstep_prev)

!---------------------------------------------------- value of (x,y) and derivatives on Gaussian points
#ifdef altcs
psieq = 0.d0; psieq_s = 0.d0; psieq_t = 0.d0
#endif
x_g  = 0.d0; x_s   = 0.d0; x_t   = 0.d0; x_p = 0.d0; x_st  = 0.d0; x_ss  = 0.d0; x_tt  = 0.d0; x_sp = 0.d0; x_tp = 0.d0; x_pp = 0.d0;
y_g  = 0.d0; y_s   = 0.d0; y_t   = 0.d0; y_p = 0.d0; y_st  = 0.d0; y_ss  = 0.d0; y_tt  = 0.d0; y_sp = 0.d0; y_tp = 0.d0; y_pp = 0.d0;
eq_g = 0.d0; eq_s  = 0.d0; eq_t  = 0.d0; eq_st = 0.d0; eq_ss = 0.d0; eq_tt = 0.d0;
eq_p = 0.d0; eq_pp = 0.d0; eq_sp = 0.d0; eq_tp = 0.d0

eq = 0.d0
press_gvec = 0.d0; B_gvec = 0.d0

do i=1,n_vertex_max
 do j=1,n_order+1

   do ms=1, n_gauss
     do mt=1, n_gauss

       press_gvec(ms,mt) = press_gvec(ms,mt) + nodes(i)%pressure(j)*element%size(i,j)*H(i,j,ms,mt)

#ifdef altcs
       psieq(ms,mt)   = psieq(ms,mt)   + nodes(i)%psi_eq(j)*element%size(i,j)*H(i,j,ms,mt)
       psieq_s(ms,mt) = psieq_s(ms,mt) + nodes(i)%psi_eq(j)*element%size(i,j)*H_s(i,j,ms,mt)
       psieq_t(ms,mt) = psieq_t(ms,mt) + nodes(i)%psi_eq(j)*element%size(i,j)*H_t(i,j,ms,mt)
#endif

       do mp=1,n_plane
         do in=1,n_coord_tor
           x_g(mp,ms,mt)  = x_g(mp,ms,mt)  + nodes(i)%x(in,j,1) * element%size(i,j) * H(i,j,ms,mt)    * HZ_coord(in,mp)
           x_s(mp,ms,mt)  = x_s(mp,ms,mt)  + nodes(i)%x(in,j,1) * element%size(i,j) * H_s(i,j,ms,mt)  * HZ_coord(in,mp)
           x_t(mp,ms,mt)  = x_t(mp,ms,mt)  + nodes(i)%x(in,j,1) * element%size(i,j) * H_t(i,j,ms,mt)  * HZ_coord(in,mp)
           x_p(mp,ms,mt)  = x_p(mp,ms,mt)  + nodes(i)%x(in,j,1) * element%size(i,j) * H(i,j,ms,mt)    * HZ_coord_p(in,mp)
           x_ss(mp,ms,mt) = x_ss(mp,ms,mt) + nodes(i)%x(in,j,1) * element%size(i,j) * H_ss(i,j,ms,mt) * HZ_coord(in,mp)
           x_st(mp,ms,mt) = x_st(mp,ms,mt) + nodes(i)%x(in,j,1) * element%size(i,j) * H_st(i,j,ms,mt) * HZ_coord(in,mp)
           x_tt(mp,ms,mt) = x_tt(mp,ms,mt) + nodes(i)%x(in,j,1) * element%size(i,j) * H_tt(i,j,ms,mt) * HZ_coord(in,mp)
           x_sp(mp,ms,mt) = x_sp(mp,ms,mt) + nodes(i)%x(in,j,1) * element%size(i,j) * H_s(i,j,ms,mt)  * HZ_coord_p(in,mp)
           x_tp(mp,ms,mt) = x_tp(mp,ms,mt) + nodes(i)%x(in,j,1) * element%size(i,j) * H_t(i,j,ms,mt)  * HZ_coord_p(in,mp)
           x_pp(mp,ms,mt) = x_pp(mp,ms,mt) + nodes(i)%x(in,j,1) * element%size(i,j) * H(i,j,ms,mt)    * HZ_coord_pp(in,mp)

           y_g(mp,ms,mt)  = y_g(mp,ms,mt)  + nodes(i)%x(in,j,2) * element%size(i,j) * H(i,j,ms,mt)    * HZ_coord(in,mp)
           y_s(mp,ms,mt)  = y_s(mp,ms,mt)  + nodes(i)%x(in,j,2) * element%size(i,j) * H_s(i,j,ms,mt)  * HZ_coord(in,mp)
           y_t(mp,ms,mt)  = y_t(mp,ms,mt)  + nodes(i)%x(in,j,2) * element%size(i,j) * H_t(i,j,ms,mt)  * HZ_coord(in,mp)
           y_p(mp,ms,mt)  = y_p(mp,ms,mt)  + nodes(i)%x(in,j,2) * element%size(i,j) * H(i,j,ms,mt)    * HZ_coord_p(in,mp)
           y_ss(mp,ms,mt) = y_ss(mp,ms,mt) + nodes(i)%x(in,j,2) * element%size(i,j) * H_ss(i,j,ms,mt) * HZ_coord(in,mp)
           y_st(mp,ms,mt) = y_st(mp,ms,mt) + nodes(i)%x(in,j,2) * element%size(i,j) * H_st(i,j,ms,mt) * HZ_coord(in,mp)
           y_tt(mp,ms,mt) = y_tt(mp,ms,mt) + nodes(i)%x(in,j,2) * element%size(i,j) * H_tt(i,j,ms,mt) * HZ_coord(in,mp)
           y_sp(mp,ms,mt) = y_sp(mp,ms,mt) + nodes(i)%x(in,j,2) * element%size(i,j) * H_s(i,j,ms,mt)  * HZ_coord_p(in,mp)
           y_tp(mp,ms,mt) = y_tp(mp,ms,mt) + nodes(i)%x(in,j,2) * element%size(i,j) * H_t(i,j,ms,mt)  * HZ_coord_p(in,mp)
           y_pp(mp,ms,mt) = y_pp(mp,ms,mt) + nodes(i)%x(in,j,2) * element%size(i,j) * H(i,j,ms,mt)    * HZ_coord_pp(in,mp)

           B_gvec(1,mp,ms,mt) = B_gvec(1,mp,ms,mt) + nodes(i)%b_field(in,j,1)*element%size(i,j)*H(i,j,ms,mt)*HZ_coord(in,mp)
           B_gvec(2,mp,ms,mt) = B_gvec(2,mp,ms,mt) + nodes(i)%b_field(in,j,2)*element%size(i,j)*H(i,j,ms,mt)*HZ_coord(in,mp)
           B_gvec(3,mp,ms,mt) = B_gvec(3,mp,ms,mt) + nodes(i)%b_field(in,j,3)*element%size(i,j)*H(i,j,ms,mt)*HZ_coord(in,mp)
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
           enddo
         enddo
       enddo
     enddo
   enddo
 enddo
enddo

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

     ! Values
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
     eq(1:n_var,0,0,2,1) = eq_pp(mp,:,ms,mt) - x_pp(mp,ms,mt)*eq(1:n_var,1,0,0,1) - 2.d0*(x_p(mp,ms,mt)*eq_px + y_p(mp,ms,mt)*eq_py)                   &
                         - y_pp(mp,ms,mt)*eq(1:n_var,0,1,0,1) + 2.d0*(x_p(mp,ms,mt)*x_p_x*eq(1:n_var,1,0,0,1) + x_p(mp,ms,mt)*y_p_x*eq(1:n_var,0,1,0,1)&
                         + y_p(mp,ms,mt)*x_p_y*eq(1:n_var,1,0,0,1) + y_p(mp,ms,mt)*y_p_y*eq(1:n_var,0,1,0,1)) + x_p(mp,ms,mt)**2*eq(1:n_var,2,0,0,1)   &
                         + 2.d0*x_p(mp,ms,mt)*y_p(mp,ms,mt)*eq(1:n_var,1,1,0,1) + y_p(mp,ms,mt)**2*eq(1:n_var,0,2,0,1)
     eq(1:n_var,1,0,1,1) = eq_px - x_p_x*eq(1:n_var,1,0,0,1) - x_p(mp,ms,mt)*eq(1:n_var,2,0,0,1) - y_p_x*eq(1:n_var,0,1,0,1) - y_p(mp,ms,mt)*eq(1:n_var,1,1,0,1)
     eq(1:n_var,0,1,1,1) = eq_py - x_p_y*eq(1:n_var,1,0,0,1) - x_p(mp,ms,mt)*eq(1:n_var,1,1,0,1) - y_p_y*eq(1:n_var,0,1,0,1) - y_p(mp,ms,mt)*eq(1:n_var,0,2,0,1)

     eq(n_var+3,:,:,:,1) = get_chi(x_g(mp,ms,mt),y_g(mp,ms,mt),phi)  ! Vacuum scalar magnetic potential (chi) and field (grad chi)
     eq(n_var+4,0,0,0,1) = x_g(mp,ms,mt); eq(n_var+4,1,0,0,1) = 1.d0 ! Cylindrical R coordinate
     
     eq(n_var+5,0,0,0,1) = mu_zero*press_gvec(ms,mt)  ! Pressure, as imported from GVEC
     eq(n_var+6:n_var+8,0,0,0,1) = B_gvec(:,mp,ms,mt) ! Magnetic field, as imported from GVEC
     
     psi_norm = get_psi_n(eq(1,0,0,0,1), y_g(mp,ms,mt))
     
     ! The Psi in the equations differs by a factor of F0 from the normal JOREK Psi
     eq(1,:,:,:,1) = eq(1,:,:,:,1)/F0
     eq(3,:,:,:,1) = eq(3,:,:,:,1)/F0

     ! Auxiliary variables (aux)
#ifdef DEBUG
     eq(n_var+9,0,0,0,:) = eval(thread_eq(tid)%aBv2seq);  eq(n_var+9,1,0,0,:) = eval(thread_eq(tid)%aBv2xseq)
     eq(n_var+9,0,1,0,:) = eval(thread_eq(tid)%aBv2yseq); eq(n_var+9,0,0,1,:) = eval(thread_eq(tid)%aBv2pseq)
#else
#include "aux_unreadable.h"
#endif
     
     n_tor_local = i_tor_max - i_tor_min + 1
     do i=1,n_vertex_max

       do j=1,n_order+1

         do im=i_tor_min,i_tor_max

           index_ij = n_tor_local*n_var*(n_order+1)*(i-1) + n_tor_local*n_var*(j-1) + im - i_tor_min + 1   ! index in the ELM matrix

           ! Test function (v)
           eq(n_var+1,0,0,0,1) =  H(i,j,ms,mt)*element%size(i,j)*HZ(im,mp)
           eq(n_var+1,1,0,0,1) = (y_t(mp,ms,mt)*h_s(i,j,ms,mt) - y_s(mp,ms,mt)*h_t(i,j,ms,mt))*element%size(i,j)*HZ(im,mp)/xjac
           eq(n_var+1,0,1,0,1) = (-x_t(mp,ms,mt)*h_s(i,j,ms,mt) + x_s(mp,ms,mt)*h_t(i,j,ms,mt))*element%size(i,j)*HZ(im,mp)/xjac
           eq(n_var+1,0,0,1,1) = H(i,j,ms,mt)*element%size(i,j)*HZ_p(im,mp) - eq(n_var+1,1,0,0,1)*x_p(mp,ms,mt) - eq(n_var+1,0,1,0,1)*y_p(mp,ms,mt)
           eq(n_var+1,2,0,0,1) = (h_ss(i,j,ms,mt)*y_t(mp,ms,mt)**2 - 2.d0*h_st(i,j,ms,mt)*y_s(mp,ms,mt)*y_t(mp,ms,mt)                             &
	                           + h_tt(i,j,ms,mt)*y_s(mp,ms,mt)**2                                                                                 &
	                           + h_s(i,j,ms,mt)*(y_st(mp,ms,mt)*y_t(mp,ms,mt) - y_tt(mp,ms,mt)*y_s(mp,ms,mt))                                     &
                               + h_t(i,j,ms,mt)*(y_st(mp,ms,mt)*y_s(mp,ms,mt) - y_ss(mp,ms,mt)*y_t(mp,ms,mt)))*element%size(i,j)*HZ(im,mp)/xjac**2&
                               - xjac_x*(h_s(i,j,ms,mt)*y_t(mp,ms,mt) - h_t(i,j,ms,mt)*y_s(mp,ms,mt))*element%size(i,j)*HZ(im,mp)/xjac**2
           eq(n_var+1,0,2,0,1) = (h_ss(i,j,ms,mt)*x_t(mp,ms,mt)**2 - 2.d0*h_st(i,j,ms,mt)*x_s(mp,ms,mt)*x_t(mp,ms,mt)                             &
                               + h_tt(i,j,ms,mt)*x_s(mp,ms,mt)**2                                                                                 &
                               + h_s(i,j,ms,mt)*(x_st(mp,ms,mt)*x_t(mp,ms,mt) - x_tt(mp,ms,mt)*x_s(mp,ms,mt))                                     &
                               + h_t(i,j,ms,mt)*(x_st(mp,ms,mt)*x_s(mp,ms,mt) - x_ss(mp,ms,mt)*x_t(mp,ms,mt)))*element%size(i,j)*HZ(im,mp)/xjac**2&
           	                   - xjac_y*(-h_s(i,j,ms,mt)*x_t(mp,ms,mt) + h_t(i,j,ms,mt)*x_s(mp,ms,mt))*element%size(i,j)*HZ(im,mp)/xjac**2
           eq(n_var+1,1,1,0,1) = (-h_ss(i,j,ms,mt)*y_t(mp,ms,mt)*x_t(mp,ms,mt) - h_tt(i,j,ms,mt)*x_s(mp,ms,mt)*y_s(mp,ms,mt)                      &
       	                       + h_st(i,j,ms,mt)*(y_s(mp,ms,mt)*x_t(mp,ms,mt) + y_t(mp,ms,mt)*x_s(mp,ms,mt))                                      &
                               - h_s(i,j,ms,mt)*(x_st(mp,ms,mt)*y_t(mp,ms,mt) - x_tt(mp,ms,mt)*y_s(mp,ms,mt))                                     &
                               - h_t(i,j,ms,mt)*(x_st(mp,ms,mt)*y_s(mp,ms,mt) - x_ss(mp,ms,mt)*y_t(mp,ms,mt)))*element%size(i,j)*HZ(im,mp)/xjac**2&
                               - xjac_x*(-h_s(i,j,ms,mt)*x_t(mp,ms,mt) + h_t(i,j,ms,mt)*x_s(mp,ms,mt))*element%size(i,j)*HZ(im,mp)/xjac**2
           v_px                = (y_t(mp,ms,mt)*h_s(i,j,ms,mt) - y_s(mp,ms,mt)*h_t(i,j,ms,mt))*element%size(i,j)*HZ_p(im,mp)/xjac
           v_py                = (-x_t(mp,ms,mt)*h_s(i,j,ms,mt) + x_s(mp,ms,mt)*h_t(i,j,ms,mt))*element%size(i,j)*HZ_p(im,mp)/xjac
           eq(n_var+1,0,0,2,1) = H(i,j,ms,mt)*element%size(i,j)*HZ_pp(im,mp) - x_pp(mp,ms,mt)*eq(n_var+1,1,0,0,1) - 2.d0*(x_p(mp,ms,mt)*v_px &
                               + y_p(mp,ms,mt)*v_py) - y_pp(mp,ms,mt)*eq(n_var+1,0,1,0,1) + 2.d0*(x_p(mp,ms,mt)*x_p_x*eq(n_var+1,1,0,0,1) &
                               + x_p(mp,ms,mt)*y_p_x*eq(n_var+1,0,1,0,1) + y_p(mp,ms,mt)*x_p_y*eq(n_var+1,1,0,0,1) + y_p(mp,ms,mt)*y_p_y*eq(n_var+1,0,1,0,1)) &
                               + x_p(mp,ms,mt)**2*eq(n_var+1,2,0,0,1) + 2.d0*x_p(mp,ms,mt)*y_p(mp,ms,mt)*eq(n_var+1,1,1,0,1) &
                               + y_p(mp,ms,mt)**2*eq(n_var+1,0,2,0,1)
           eq(n_var+1,1,0,1,1) = v_px - x_p_x*eq(n_var+1,1,0,0,1) - x_p(mp,ms,mt)*eq(n_var+1,2,0,0,1) - y_p_x*eq(n_var+1,0,1,0,1) - y_p(mp,ms,mt)*eq(n_var+1,1,1,0,1)
           eq(n_var+1,0,1,1,1) = v_py - x_p_y*eq(n_var+1,1,0,0,1) - y_p(mp,ms,mt)*eq(n_var+1,0,2,0,1) - y_p_y*eq(n_var+1,0,1,0,1) - x_p(mp,ms,mt)*eq(n_var+1,1,1,0,1)

#ifdef DEBUG
           rhs_ij_1 = eval(thread_eq(tid)%rhs1seq)*BigR*xjac
           rhs_ij_3 = eval(thread_eq(tid)%rhs3seq)*BigR*xjac
#else
#include "rhs_unreadable.h"

           rhs_ij_1 = rhs_ij_1*BigR*xjac
           rhs_ij_3 = rhs_ij_3*BigR*xjac
#endif

           ij1 = index_ij
           ij2 = index_ij + 1*n_tor_local
           ij3 = index_ij + 2*n_tor_local
           ij4 = index_ij + 3*n_tor_local
           ij5 = index_ij + 4*n_tor_local
           ij6 = index_ij + 5*n_tor_local

           RHS(ij1) = RHS(ij1) + rhs_ij_1(1)*wst
           RHS(ij3) = RHS(ij3) + rhs_ij_3(1)*wst

           do k=1,n_vertex_max

             do l=1,n_order+1

               do in=i_tor_min,i_tor_max

                 ! Unknown increments to next time step (delta u^n)
                 eq(n_var+2,0,0,0,1) = H(k,l,ms,mt)*element%size(k,l)*HZ(in,mp)
                 eq(n_var+2,1,0,0,1) = (y_t(mp,ms,mt)*h_s(k,l,ms,mt) - y_s(mp,ms,mt)*h_t(k,l,ms,mt))*element%size(k,l)*HZ(in,mp)/xjac
                 eq(n_var+2,0,1,0,1) = (-x_t(mp,ms,mt)*h_s(k,l,ms,mt) + x_s(mp,ms,mt)*h_t(k,l,ms,mt))*element%size(k,l)*HZ(in,mp)/xjac
                 eq(n_var+2,0,0,1,1) = H(k,l,ms,mt)*element%size(k,l)*HZ_p(in,mp) - eq(n_var+2,1,0,0,1)*x_p(mp,ms,mt) - eq(n_var+2,0,1,0,1)*y_p(mp,ms,mt)
                 eq(n_var+2,2,0,0,1) = (h_ss(k,l,ms,mt)*y_t(mp,ms,mt)**2 - 2.d0*h_st(k,l,ms,mt)*y_s(mp,ms,mt)*y_t(mp,ms,mt)                             &
                                     + h_tt(k,l,ms,mt)*y_s(mp,ms,mt)**2                                                                                 &
                                     + h_s(k,l,ms,mt)*(y_st(mp,ms,mt)*y_t(mp,ms,mt) - y_tt(mp,ms,mt)*y_s(mp,ms,mt))                                     &
                                     + h_t(k,l,ms,mt)*(y_st(mp,ms,mt)*y_s(mp,ms,mt) - y_ss(mp,ms,mt)*y_t(mp,ms,mt)))*element%size(k,l)*HZ(in,mp)/xjac**2&	
                                     - xjac_x*(h_s(k,l,ms,mt)*y_t(mp,ms,mt) - h_t(k,l,ms,mt)*y_s(mp,ms,mt))*element%size(k,l)*HZ(in,mp)/xjac**2
                 eq(n_var+2,0,2,0,1) = (h_ss(k,l,ms,mt)*x_t(mp,ms,mt)**2 - 2.d0*h_st(k,l,ms,mt)*x_s(mp,ms,mt)*x_t(mp,ms,mt)                             &
                                     + h_tt(k,l,ms,mt)*x_s(mp,ms,mt)**2                                                                                 &
                                     + h_s(k,l,ms,mt)*(x_st(mp,ms,mt)*x_t(mp,ms,mt) - x_tt(mp,ms,mt)*x_s(mp,ms,mt))                                     &
                                     + h_t(k,l,ms,mt)*(x_st(mp,ms,mt)*x_s(mp,ms,mt) - x_ss(mp,ms,mt)*x_t(mp,ms,mt)))*element%size(k,l)*HZ(in,mp)/xjac**2&
                                     - xjac_y*(-h_s(k,l,ms,mt)*x_t(mp,ms,mt) + h_t(k,l,ms,mt)*x_s(mp,ms,mt))*element%size(k,l)*HZ(in,mp)/xjac**2
                 eq(n_var+2,1,1,0,1) = (-h_ss(k,l,ms,mt)*y_t(mp,ms,mt)*x_t(mp,ms,mt) - h_tt(k,l,ms,mt)*x_s(mp,ms,mt)*y_s(mp,ms,mt)                      &
     	                             + h_st(k,l,ms,mt)*(y_s(mp,ms,mt)*x_t(mp,ms,mt)  + y_t(mp,ms,mt)*x_s(mp,ms,mt))                                     &
                                     - h_s(k,l,ms,mt)*(x_st(mp,ms,mt)*y_t(mp,ms,mt) - x_tt(mp,ms,mt)*y_s(mp,ms,mt))                                     &
                                     - h_t(k,l,ms,mt)*(x_st(mp,ms,mt)*y_s(mp,ms,mt) - x_ss(mp,ms,mt)*y_t(mp,ms,mt)))*element%size(k,l)*HZ(in,mp)/xjac**2&
                                     - xjac_x*(-h_s(k,l,ms,mt)*x_t(mp,ms,mt) + h_t(k,l,ms,mt)*x_s(mp,ms,mt))*element%size(k,l)*HZ(in,mp)/xjac**2
                 u_px                = (y_t(mp,ms,mt)*h_s(k,l,ms,mt) - y_s(mp,ms,mt)*h_t(k,l,ms,mt))*element%size(k,l)*HZ_p(in,mp)/xjac
                 u_py                = (-x_t(mp,ms,mt)*h_s(k,l,ms,mt) + x_s(mp,ms,mt)*h_t(k,l,ms,mt))*element%size(k,l)*HZ_p(in,mp)/xjac
                 eq(n_var+2,0,0,2,1) = H(k,l,ms,mt)*element%size(k,l)*HZ_pp(in,mp) - x_pp(mp,ms,mt)*eq(n_var+2,1,0,0,1) - 2.d0*(x_p(mp,ms,mt)*u_px &
                                     + y_p(mp,ms,mt)*u_py) - y_pp(mp,ms,mt)*eq(n_var+2,0,1,0,1) + 2.d0*(x_p(mp,ms,mt)*x_p_x*eq(n_var+2,1,0,0,1) &
                                     + x_p(mp,ms,mt)*y_p_x*eq(n_var+2,0,1,0,1) + y_p(mp,ms,mt)*x_p_y*eq(n_var+2,1,0,0,1) &
                                     + y_p(mp,ms,mt)*y_p_y*eq(n_var+2,0,1,0,1)) + x_p(mp,ms,mt)**2*eq(n_var+2,2,0,0,1)             &
                                     + 2.d0*x_p(mp,ms,mt)*y_p(mp,ms,mt)*eq(n_var+2,1,1,0,1) + y_p(mp,ms,mt)**2*eq(n_var+2,0,2,0,1)
                 eq(n_var+2,1,0,1,1) = u_px - x_p_x*eq(n_var+2,1,0,0,1) - x_p(mp,ms,mt)*eq(n_var+2,2,0,0,1) - y_p_x*eq(n_var+2,0,1,0,1) &
                                     - y_p(mp,ms,mt)*eq(n_var+2,1,1,0,1)
                 eq(n_var+2,0,1,1,1) = u_py - x_p_y*eq(n_var+2,1,0,0,1) - y_p(mp,ms,mt)*eq(n_var+2,0,2,0,1) - y_p_y*eq(n_var+2,0,1,0,1) &
                                     - x_p(mp,ms,mt)*eq(n_var+2,1,1,0,1)
                 
                 index_kl = n_tor_local*n_var*(n_order+1)*(k-1) + n_tor_local*n_var*(l-1) + in - i_tor_min + 1   ! index in the ELM matrix
                 
#ifdef DEBUG
!---------------------------------------------------------------- equation 1
                 amat_11 = eval(thread_eq(tid)%amat11seq)*BigR*xjac/F0
                 amat_13 = eval(thread_eq(tid)%amat13seq)*BigR*xjac/F0

!---------------------------------------------------------------- equation 2
                 amat_22 = eval(thread_eq(tid)%amat22seq)*BigR*xjac

!---------------------------------------------------------------- equation 3
                 amat_33 = eval(thread_eq(tid)%amat33seq)*BigR*xjac/F0

!---------------------------------------------------------------- equation 4
                 amat_44 = eval(thread_eq(tid)%amat44seq)*BigR*xjac
                 
!---------------------------------------------------------------- equation 5
                 amat_55 = eval(thread_eq(tid)%amat55seq)*BigR*xjac
                 
!---------------------------------------------------------------- equation 6
                 amat_66 = eval(thread_eq(tid)%amat66seq)*BigR*xjac
#else
#include "amat_unreadable.h"

                 amat_11 = amat_11*BigR*xjac/F0
                 amat_13 = amat_13*BigR*xjac/F0
                 amat_22 = amat_22*BigR*xjac
                 amat_33 = amat_33*BigR*xjac/F0
                 amat_44 = amat_44*BigR*xjac
                 amat_55 = amat_55*BigR*xjac
                 amat_66 = amat_66*BigR*xjac
#endif

                 kl1 = index_kl
                 kl2 = index_kl + 1*n_tor_local
                 kl3 = index_kl + 2*n_tor_local
                 kl4 = index_kl + 3*n_tor_local
                 kl5 = index_kl + 4*n_tor_local
                 kl6 = index_kl + 5*n_tor_local

                 ELM(ij1,kl1) =  ELM(ij1,kl1) + wst*amat_11(1)
                 ELM(ij1,kl3) =  ELM(ij1,kl3) + wst*amat_13(1)

                 ELM(ij2,kl2) =  ELM(ij2,kl2) + wst*amat_22(1)

                 ELM(ij3,kl3) =  ELM(ij3,kl3) + wst*amat_33(1)

                 ELM(ij4,kl4) =  ELM(ij4,kl4) + wst*amat_44(1)

                 ELM(ij5,kl5) =  ELM(ij5,kl5) + wst*amat_55(1)

                 ELM(ij6,kl6) =  ELM(ij6,kl6) + wst*amat_66(1)
               enddo
             enddo
           enddo

         enddo
       enddo

     enddo
     
   enddo

 enddo
enddo

return

end subroutine element_matrix
end module mod_elt_matrix
