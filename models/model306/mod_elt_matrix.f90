module mod_elt_matrix

  implicit none

contains

subroutine element_matrix(element, nodes, xpoint2, xcase2, R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint, Z_xpoint, ELM,RHS, tid, i_tor_min, i_tor_max)
!---------------------------------------------------------------
! calculates the matrix contribution of one element
! using the equations for model306, which is the same as model303, but solves
! two additional current equations, which sum to the total ECCD current.  
! These currents are injected inside the magnetic islands.
! version to be debugged 28/2/2014
! Author:                                                                    
!   Jane Pratt - J.L.Pratt@differ.nl                                
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
use mod_sources

implicit none

type (type_element)   :: element
type (type_node)      :: nodes(n_vertex_max)

real*8, dimension (:,:), allocatable  :: ELM
real*8, dimension (:)  , allocatable  :: RHS
integer, intent(in)                   :: tid, i_tor_min, i_tor_max

integer    :: i, j, ms, mt, mp, k, l, index_ij, index_kl, index, xcase2, n_tor_local
integer    :: in, im, ij1, ij2, ij3, ij4, ij5, ij6, ij7, kl1, kl2, kl3, kl4, kl5, kl6, kl7
real*8     :: wst, xjac, xjac_s, xjac_t, xjac_x, xjac_y, BigR, r2, phi, delta_phi, eps_cyl
real*8     :: current_source(n_gauss,n_gauss), particle_source(n_gauss,n_gauss), heat_source(n_gauss,n_gauss)
real*8     :: source_volume, source_pellet, source_pellet2
real*8     :: R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint(2), Z_xpoint(2), dj_dpsi, dj_dz
real*8     :: Bgrad_rho_star,     Bgrad_rho,     Bgrad_T_star,  Bgrad_T, BB2
real*8     :: Bgrad_rho_star_psi, Bgrad_rho_psi, Bgrad_rho_rho, Bgrad_T_star_psi, Bgrad_T_psi, Bgrad_T_T, BB2_psi
real*8     :: rhs_ij_1,   rhs_ij_2,   rhs_ij_3,   rhs_ij_4,   rhs_ij_5,   rhs_ij_6, rhs_ij_7
real*8     :: ZK_prof, D_prof, psi_norm, theta, zeta, delta_u_x, delta_u_y, delta_ps_x, delta_ps_y, ZKpar_T, dZKpar_dT

#ifdef altcs
real*8     :: psieq_R, psieq_Z  ! more accurate current
#endif
real*8     :: v, v_x, v_y, v_s, v_t, v_p, v_ss, v_st, v_tt, v_xx, v_yy, v_xy
real*8     :: ps0, ps0_x, ps0_y, ps0_p,ps0_s,ps0_t, ps0_ss, ps0_st, ps0_tt, ps0_xx, ps0_xy, ps0_yy
real*8     :: zj0, zj0_x, zj0_y, zj0_p, zj0_s, zj0_t
real*8     :: u0, u0_x, u0_y, u0_p, u0_s, u0_t, u0_ss, u0_st, u0_tt, u0_xx, u0_xy, u0_yy 
real*8     :: w0, w0_x, w0_y, w0_p, w0_s, w0_t, w0_ss, w0_st, w0_tt, w0_xx, w0_xy, w0_yy
real*8     :: Vpar0, Vpar0_x, Vpar0_y, Vpar0_p, Vpar0_s, Vpar0_t, Vpar0_ss, Vpar0_st, Vpar0_tt, vpar0_xx, vpar0_xy, vpar0_yy
real*8     :: r0, r0_x, r0_y, r0_p, r0_s, r0_t, r0_ss, r0_st, r0_tt, r0_xx, r0_xy, r0_yy, r0_hat, r0_x_hat, r0_y_hat
real*8     :: T0, T0_x, T0_y, T0_p, T0_s, T0_t, T0_ss, T0_st, T0_tt, T0_xx, T0_xy, T0_yy
real*8     :: psi, psi_x, psi_y, psi_p, psi_s, psi_t, psi_ss, psi_st, psi_tt, psi_xx, psi_yy, psi_xy
real*8     :: zj, zj_x, zj_y, zj_p, zj_s, zj_t, zj_ss, zj_st, zj_tt
real*8     :: vpar, vpar_x, vpar_y, vpar_s, vpar_t, vpar_p, vpar_ss, vpar_st, vpar_tt, vpar_xx, vpar_xy, vpar_yy
real*8     :: u, u_x, u_y, u_p, u_s, u_t, u_ss, u_st, u_tt, u_xx, u_xy, u_yy
real*8     :: w, w_x, w_y, w_p, w_s, w_t, w_ss, w_st, w_tt, w_xx, w_xy, w_yy
real*8     :: rho, rho_x, rho_y, rho_s, rho_t, rho_p, rho_ss, rho_st, rho_tt, rho_xx, rho_xy, rho_yy, rho_hat, rho_x_hat, rho_y_hat
real*8     :: T, T_x, T_y, T_s, T_t, T_p, T_ss, T_st, T_tt, T_xx, T_xy, T_yy
real*8     :: P0, P0_x, P0_y, P0_s, P0_t, P0_ss, P0_st, P0_tt, P0_p, P0_xx, P0_xy, P0_yy
real*8     :: BigR_x, vv2, eta_T, visco_T, deta_dT, d2eta_d2T, dvisco_dT, visco_num_T, eta_num_T
real*8     :: amat_11, amat_12, amat_21, amat_22, amat_23, amat_24, amat_25, amat_26, amat_33, amat_31, amat_44, amat_42
real*8     :: amat_51, amat_52, amat_55, amat_56, amat_57, amat_61, amat_62, amat_65, amat_66, amat_67, amat_16, amat_13
real*8     :: amat_71, amat_72, amat_75, amat_76, amat_77, amat_15
real*8     :: ZK_par_num, T0_ps0_x, T_ps0_x, T0_psi_x, T0_ps0_y, T_ps0_y, T0_psi_y, v_ps0_x, v_psi_x, v_ps0_y, v_psi_y
real*8     :: TG_num1, TG_num2, TG_num5, TG_num6, TG_num7
logical    :: xpoint2
!==================MB: velocity profile is kept by a source which compensating diffusion
real*8     :: Vt0,Vt0_x,Vt0_y
real*8     :: V_source(n_gauss,n_gauss)
real*8     :: dV_dpsi_source(n_gauss,n_gauss),dV_dz_source(n_gauss,n_gauss)
real*8     :: dV_dpsi2,dV_dz2,dV_dpsi_dz,dV_dpsi3,dV_dpsi_dz2,dV_dpsi2_dz
!=======================================
real*8     :: eq_zne(n_gauss,n_gauss), eq_zTe(n_gauss,n_gauss)
real*8     :: dn_dpsi,dn_dz,dn_dpsi2,dn_dz2,dn_dpsi_dz,dn_dpsi3,dn_dpsi_dz2,dn_dpsi2_dz
real*8     :: dT_dpsi,dT_dz,dT_dpsi2,dT_dz2,dT_dpsi_dz,dT_dpsi3,dT_dpsi_dz2,dT_dpsi2_dz
real*8     :: w00_xx, w00_yy 
!======================================= NEO
real*8     :: amat_27, Btheta2
real*8     :: epsil, Btheta2_psi
real*8, dimension(n_gauss,n_gauss)    :: amu_neo_prof, aki_neo_prof
!======================================= NEO

!bookmark (alterations and additions for 2 ECCD currents)
integer    :: kl8,ij8,kl9,ij9
real*8     :: zero1,zero2,zero3
real*8     :: rhs_ij_8,rhs_ij_9
real*8     :: Bgrad_jec1_star,  Bgrad_jec1,Bdot_jec1_star
real*8     :: Bgrad_jec1_star_psi, Bgrad_jec1_psi, Bgrad_jec1_jec
real*8     :: Bgrad_jec2_star,  Bgrad_jec2,Bdot_jec_star
real*8     :: Bgrad_jec2_star_psi, Bgrad_jec2_psi, Bgrad_jec2_jec
real*8     :: ec_source(n_gauss,n_gauss),radd
real*8     :: amat_81, amat_88, amat_18
real*8     :: amat_91, amat_99, amat_19
real*8     :: jec10, jec10_x, jec10_y, jec10_p, jec10_s, jec10_t
real*8     :: jec1, jec1_x, jec1_y, jec1_p, jec1_s, jec1_t
real*8     :: jec20, jec20_x, jec20_y, jec20_p, jec20_s, jec20_t
real*8     :: jec2, jec2_x, jec2_y, jec2_p, jec2_s, jec2_t

!real*8     :: BB1, BB1_psi
!real*8                             :: jbs,pnot_x !bootstrap current
!real*8, dimension(n_gauss,n_gauss) :: pnot
!!



real*8, dimension(n_gauss,n_gauss)    :: x_g, x_s, x_t
real*8, dimension(n_gauss,n_gauss)    :: x_ss, x_st, x_tt
real*8, dimension(n_gauss,n_gauss)    :: y_g, y_s, y_t
real*8, dimension(n_gauss,n_gauss)    :: y_ss, y_st, y_tt
#ifdef altcs
real*8, dimension(n_gauss,n_gauss)    :: psieq,psieq_s,psieq_t
#endif

real*8, dimension(:,:,:,:) , pointer :: eq_g, eq_s, eq_t
real*8, dimension(:,:,:,:) , pointer :: eq_p
real*8, dimension(:,:,:,:) , pointer :: eq_ss, eq_st, eq_tt   
real*8, dimension(:,:,:,:) , pointer :: delta_g, delta_s, delta_t

eq_g    => thread_struct(tid)%eq_g   
eq_s    => thread_struct(tid)%eq_s   
eq_t    => thread_struct(tid)%eq_t   
eq_p    => thread_struct(tid)%eq_p   
eq_ss   => thread_struct(tid)%eq_ss  
eq_st   => thread_struct(tid)%eq_st  
eq_tt   => thread_struct(tid)%eq_tt  
delta_g => thread_struct(tid)%delta_g
delta_s => thread_struct(tid)%delta_s
delta_t => thread_struct(tid)%delta_t

ELM = 0.d0
RHS = 0.d0
!======================================= NEO
epsil=1.d-3
!======================================= NEO

zk_par_num = 0.d0
amat_27    = 0.d0

! --- Taylor-Galerkin Stabilisation coefficients
TG_num1    = TGNUM(1); TG_num2    = TGNUM(2); TG_num5    = TGNUM(5); TG_num6    = TGNUM(6); TG_num7    = TGNUM(7);

! --- Take time evolution parameters from phys_module
theta = time_evol_theta
!zeta  = time_evol_zeta
! change zeta for variable dt
zeta  = time_evol_zeta * 2.0d0 * tstep / (tstep + tstep_prev)

!---------------------------------------------------- value of (x,y) and derivatives on Gaussian points
x_g  = 0.d0; x_s  = 0.d0; x_t  = 0.d0; x_st  = 0.d0; x_ss  = 0.d0; x_tt  = 0.d0;
y_g  = 0.d0; y_s  = 0.d0; y_t  = 0.d0; y_st  = 0.d0; y_ss  = 0.d0; y_tt  = 0.d0;
eq_g = 0.d0; eq_s = 0.d0; eq_t = 0.d0; eq_st = 0.d0; eq_ss = 0.d0; eq_tt = 0.d0; eq_p = 0.d0;

#ifdef altcs
psieq   = 0.d0; psieq_s = 0.d0; psieq_t = 0.d0
#endif
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

#ifdef altcs
       psieq(ms,mt)    = psieq(ms,mt)    + nodes(i)%psi_eq(j)   *element%size(i,j) * H(i,j,ms,mt)
       psieq_s(ms,mt)  = psieq_s(ms,mt)  + nodes(i)%psi_eq(j)   *element%size(i,j) * H_s(i,j,ms,mt)
       psieq_t(ms,mt)  = psieq_t(ms,mt)  + nodes(i)%psi_eq(j)   *element%size(i,j) * H_t(i,j,ms,mt)
#endif
       do mp=1,n_plane

         do k=1,n_var

           do in=1,n_tor

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

!!$if ( NEO ) then 
!!$   if (num_neo_file) then
!!$      write (*,*) 'profiles read from ', neo_file
!!$   else
!!$      write (*,*) 'constant neoclass coeffs: ', amu_neo_const, aki_neo_const 
!!$   endif
!!$endif
do ms=1, n_gauss
  do mt=1, n_gauss

  if (keep_current_prof) then
#ifdef altcs
    call current(xpoint2,xcase2,x_g(ms,mt),y_g(ms,mt),Z_xpoint,psieq(ms,mt),psi_axis,psi_bnd,current_source(ms,mt))
#else
    call current(xpoint2, xcase2, x_g(ms,mt),y_g(ms,mt), Z_xpoint,eq_g(1,1,ms,mt),psi_axis,psi_bnd,current_source(ms,mt))
#endif
  end if
! eccurrent_source is defined in this module, below.
  call eccurrent_source(x_g(ms,mt),y_g(ms,mt),ec_source(ms,mt))
  call sources(xpoint2,xcase2,y_g(ms,mt),Z_xpoint,eq_g(1,1,ms,mt),psi_axis,psi_bnd,particle_source(ms,mt),heat_source(ms,mt))
!=========================================MB :velocity profile
       if ( (abs(V_0) .ge. 1.e-12) .or. (num_rot)) then 
        call velocity(xpoint2, xcase2, y_g(ms,mt), Z_xpoint, eq_g(1,1,ms,mt), psi_axis, psi_bnd, V_source(ms,mt), &
                      dV_dpsi_source(ms,mt),dV_dz_source(ms,mt),dV_dpsi2,dV_dz2,dV_dpsi_dz,dV_dpsi3,dV_dpsi_dz2, dV_dpsi2_dz)
       endif
!======================================MB
    call density(xpoint2, xcase2, y_g(ms,mt), Z_xpoint, eq_g(1,1,ms,mt),psi_axis,psi_bnd,eq_zne(ms,mt), &
                 dn_dpsi,dn_dz,dn_dpsi2,dn_dz2,dn_dpsi_dz,dn_dpsi3,dn_dpsi_dz2, dn_dpsi2_dz)

    call temperature(xpoint2, xcase2, y_g(ms,mt), Z_xpoint, eq_g(1,1,ms,mt),psi_axis,psi_bnd,eq_zTe(ms,mt), &
                     dT_dpsi,dT_dz,dT_dpsi2,dT_dz2,dT_dpsi_dz,dT_dpsi3,dT_dpsi_dz2, dT_dpsi2_dz)

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

eq_zTe = eq_zTe / 2.d0	! electron temperature	


!--------------------------------------------------- sum over the Gaussian integration points
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

   !Fprof = Fprofile(ms,mt)
   !Fprof_R = (   y_t(ms,mt) * Fprofile_s(ms,mt)  - y_s(ms,mt) *Fprofile_t(ms,mt) ) / xjac
   !Fprof_Z = ( - x_t(ms,mt) * Fprofile_s(ms,mt)  + x_s(ms,mt) *Fprofile_t(ms,mt) ) / xjac


   eps_cyl = 1.d0          ! for cylinder geometry : epscyl = eps

   do mp = 1, n_plane

! magnetic potential psi
     ps0    = eq_g(mp,1,ms,mt)
     ps0_x  = (   y_t(ms,mt) * eq_s(mp,1,ms,mt) - y_s(ms,mt) * eq_t(mp,1,ms,mt) ) / xjac
     ps0_y  = ( - x_t(ms,mt) * eq_s(mp,1,ms,mt) + x_s(ms,mt) * eq_t(mp,1,ms,mt) ) / xjac
     ps0_p  = eq_p(mp,1,ms,mt)
     ps0_s  = eq_s(mp,1,ms,mt)
     ps0_t  = eq_t(mp,1,ms,mt)
     ps0_ss = eq_ss(mp,1,ms,mt)
     ps0_tt = eq_tt(mp,1,ms,mt)
     ps0_st = eq_st(mp,1,ms,mt)

! velocity stream function
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

! current
     zj0   = eq_g(mp,3,ms,mt)
     zj0_x = (   y_t(ms,mt) * eq_s(mp,3,ms,mt) - y_s(ms,mt) * eq_t(mp,3,ms,mt) ) / xjac
     zj0_y = ( - x_t(ms,mt) * eq_s(mp,3,ms,mt) + x_s(ms,mt) * eq_t(mp,3,ms,mt) ) / xjac
     zj0_p = eq_p(mp,3,ms,mt)
     zj0_s = eq_s(mp,3,ms,mt)
     zj0_t = eq_t(mp,3,ms,mt)
!=======================================MB:parallel velocity 
     Vt0   = V_source(ms,mt)
     Vt0_x = dV_dpsi_source(ms,mt)*ps0_x
     Vt0_y = dV_dz_source(ms,mt)+dV_dpsi_source(ms,mt)*ps0_y

!=======================================MB

! vorticitiy
     w0    = eq_g(mp,4,ms,mt)
     w0_x  = (   y_t(ms,mt) * eq_s(mp,4,ms,mt) - y_s(ms,mt) * eq_t(mp,4,ms,mt) ) / xjac
     w0_y  = ( - x_t(ms,mt) * eq_s(mp,4,ms,mt) + x_s(ms,mt) * eq_t(mp,4,ms,mt) ) / xjac
     w0_p  = eq_p(mp,4,ms,mt)
     w0_s  = eq_s(mp,4,ms,mt)
     w0_t  = eq_t(mp,4,ms,mt)
     w0_ss = eq_ss(mp,4,ms,mt)
     w0_tt = eq_tt(mp,4,ms,mt)
     w0_st = eq_st(mp,4,ms,mt)

! density
     r0    = eq_g(mp,5,ms,mt)
     r0_x  = (   y_t(ms,mt) * eq_s(mp,5,ms,mt) - y_s(ms,mt) * eq_t(mp,5,ms,mt) ) / xjac
     r0_y  = ( - x_t(ms,mt) * eq_s(mp,5,ms,mt) + x_s(ms,mt) * eq_t(mp,5,ms,mt) ) / xjac
     r0_p  = eq_p(mp,5,ms,mt)
     r0_s  = eq_s(mp,5,ms,mt)
     r0_t  = eq_t(mp,5,ms,mt)
     r0_ss = eq_ss(mp,5,ms,mt)
     r0_st = eq_st(mp,5,ms,mt)
     r0_tt = eq_tt(mp,5,ms,mt)

     r0_hat   = BigR**2. * abs(r0)
     r0_x_hat = 2.d0 * BigR * BigR_x  * r0 + BigR**2. * r0_x
     r0_y_hat = BigR**2. * r0_y

! temperature
     T0    = eq_g(mp,6,ms,mt)
     T0_x  = (   y_t(ms,mt) * eq_s(mp,6,ms,mt) - y_s(ms,mt) * eq_t(mp,6,ms,mt) ) / xjac
     T0_y  = ( - x_t(ms,mt) * eq_s(mp,6,ms,mt) + x_s(ms,mt) * eq_t(mp,6,ms,mt) ) / xjac
     T0_p  = eq_p(mp,6,ms,mt)
     T0_s  = eq_s(mp,6,ms,mt)
     T0_t  = eq_t(mp,6,ms,mt)
     T0_ss = eq_ss(mp,6,ms,mt)
     T0_tt = eq_tt(mp,6,ms,mt)
     T0_st = eq_st(mp,6,ms,mt)

! parallel velocity
     Vpar0    = eq_g(mp,7,ms,mt)
     Vpar0_x  = (   y_t(ms,mt) * eq_s(mp,7,ms,mt) - y_s(ms,mt) * eq_t(mp,7,ms,mt) ) / xjac
     Vpar0_y  = ( - x_t(ms,mt) * eq_s(mp,7,ms,mt) + x_s(ms,mt) * eq_t(mp,7,ms,mt) ) / xjac
     Vpar0_p  = eq_p(mp,7,ms,mt)
     Vpar0_s  = eq_s(mp,7,ms,mt)
     Vpar0_t  = eq_t(mp,7,ms,mt)
     Vpar0_ss = eq_ss(mp,7,ms,mt)
     Vpar0_st = eq_st(mp,7,ms,mt)
     Vpar0_tt = eq_tt(mp,7,ms,mt)

!ECCD current 1
     jec10    = eq_g(mp,8,ms,mt)
     jec10_x  = (   y_t(ms,mt) * eq_s(mp,8,ms,mt) - y_s(ms,mt) *eq_t(mp,8,ms,mt)) / xjac
     jec10_y  = ( - x_t(ms,mt) * eq_s(mp,8,ms,mt) + x_s(ms,mt) *eq_t(mp,8,ms,mt)) / xjac
     jec10_p  = eq_p(mp,8,ms,mt)
     jec10_s  = eq_s(mp,8,ms,mt)
     jec10_t  = eq_t(mp,8,ms,mt)

!ECCD current 2
     jec20    = eq_g(mp,9,ms,mt)
     jec20_x  = (   y_t(ms,mt) * eq_s(mp,9,ms,mt) - y_s(ms,mt) *eq_t(mp,9,ms,mt))/ xjac
     jec20_y  = ( - x_t(ms,mt) * eq_s(mp,9,ms,mt) + x_s(ms,mt) *eq_t(mp,9,ms,mt))/ xjac
     jec20_p  = eq_p(mp,9,ms,mt)
     jec20_s  = eq_s(mp,9,ms,mt)
     jec20_t  = eq_t(mp,9,ms,mt)

! pressure
     P0    = abs(r0 * T0)
     P0_x  = r0_x * T0 + r0 * T0_x
     P0_y  = r0_y * T0 + r0 * T0_y
     P0_s  = r0_s * T0 + r0 * T0_s
     P0_t  = r0_t * T0 + r0 * T0_t
     P0_p  = r0_p * T0 + r0 * T0_p
     P0_ss = r0_ss * T0 + 2.d0 * r0_s * T0_s + r0 * T0_ss
     P0_tt = r0_tt * T0 + 2.d0 * r0_t * T0_t + r0 * T0_tt
     P0_st = r0_st * T0 + r0_s * T0_t + r0_t * T0_s + r0 * T0_st

! second derivatives:
     ps0_xx = (ps0_ss * y_t(ms,mt)**2. - 2.d0*ps0_st * y_s(ms,mt)*y_t(ms,mt) + ps0_tt * y_s(ms,mt)**2. &
             + ps0_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                            &
             + ps0_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) ) )    / xjac**2.             &
             - xjac_x * (ps0_s* y_t(ms,mt) - ps0_t * y_s(ms,mt))  / xjac**2.

     ps0_yy = (ps0_ss * x_t(ms,mt)**2. - 2.d0*ps0_st * x_s(ms,mt)*x_t(ms,mt) + ps0_tt * x_s(ms,mt)**2. &
             + ps0_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                              &
             + ps0_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) ) )    / xjac**2.              &
             - xjac_y * (- ps0_s * x_t(ms,mt) + ps0_t * x_s(ms,mt) )  / xjac**2.

     ps0_xy = (- ps0_ss * y_t(ms,mt)*x_t(ms,mt) - ps0_tt * x_s(ms,mt)*y_s(ms,mt)                     &
              + ps0_st * (y_s(ms,mt)*x_t(ms,mt)  + y_t(ms,mt)*x_s(ms,mt)  )                          &
              - ps0_s  * (x_st(ms,mt)*y_t(ms,mt) - x_tt(ms,mt)*y_s(ms,mt) )                          &
              - ps0_t * (x_st(ms,mt)*y_s(ms,mt)  - x_ss(ms,mt)*y_t(ms,mt) ) )  / xjac**2.            &
              - xjac_x * (- ps0_s * x_t(ms,mt) + ps0_t * x_s(ms,mt) )   / xjac**2.

     u0_xx = (u0_ss * y_t(ms,mt)**2. - 2.d0*u0_st * y_s(ms,mt)*y_t(ms,mt) + u0_tt * y_s(ms,mt)**2     &
           + u0_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                                &
           + u0_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) ) )      / xjac**2.              &
           - xjac_x * (u0_s * y_t(ms,mt) - u0_t * y_s(ms,mt)) / xjac**2.

     u0_yy = (u0_ss * x_t(ms,mt)**2. - 2.d0*u0_st * x_s(ms,mt)*x_t(ms,mt) + u0_tt * x_s(ms,mt)**2.    &
           + u0_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                                &
           + u0_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) ) )      / xjac**2.              &
           - xjac_y * (- u0_s * x_t(ms,mt) + u0_t * x_s(ms,mt) ) / xjac**2.

     u0_xy = (- u0_ss * y_t(ms,mt)*x_t(ms,mt) - u0_tt * x_s(ms,mt)*y_s(ms,mt)                        &
              + u0_st * (y_s(ms,mt)*x_t(ms,mt)  + y_t(ms,mt)*x_s(ms,mt)  )                           &
              - u0_s  * (x_st(ms,mt)*y_t(ms,mt) - x_tt(ms,mt)*y_s(ms,mt) )                           &
              - u0_t * (x_st(ms,mt)*y_s(ms,mt)  - x_ss(ms,mt)*y_t(ms,mt) )  )  / xjac**2.            &
              - xjac_x * (- u0_s * x_t(ms,mt) + u0_t * x_s(ms,mt) )   / xjac**2.

!     w0_xx = (w0_ss * y_t(ms,mt)**2 - 2.d0*w0_st * y_s(ms,mt)*y_t(ms,mt) + w0_tt * y_s(ms,mt)**2     &
!            + w0_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                              &
!            + w0_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) ) )     / xjac**2              &
!            - xjac_x * (w0_s* y_t(ms,mt) - w0_t * y_s(ms,mt))  / xjac**2

!     w0_yy = (w0_ss * x_t(ms,mt)**2 - 2.d0*w0_st * x_s(ms,mt)*x_t(ms,mt) + w0_tt * x_s(ms,mt)**2     &
!            + w0_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                              &
!            + w0_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) ) )     / xjac**2              &
!            - xjac_y * (- w0_s * x_t(ms,mt) + w0_t * x_s(ms,mt) )  / xjac**2

     w0_xy = (- w0_ss * y_t(ms,mt)*x_t(ms,mt) - w0_tt * x_s(ms,mt)*y_s(ms,mt)                        &
              + w0_st * (y_s(ms,mt)*x_t(ms,mt)  + y_t(ms,mt)*x_s(ms,mt)  )                           &
              - w0_s  * (x_st(ms,mt)*y_t(ms,mt) - x_tt(ms,mt)*y_s(ms,mt) )                           &
              - w0_t  * (x_st(ms,mt)*y_s(ms,mt) - x_ss(ms,mt)*y_t(ms,mt) ) )  / xjac**2.              &
              - xjac_x * (- w0_s * x_t(ms,mt) + w0_t * x_s(ms,mt) )   / xjac**2.

!----------------- simplified version of 2nd derivatives (for some unknown reason this is more stable!)
     w0_xx = (w0_ss * y_t(ms,mt)**2 - 2.d0*w0_st * y_s(ms,mt)*y_t(ms,mt) + w0_tt * y_s(ms,mt)**2. )     / xjac**2.
     w0_yy = (w0_ss * x_t(ms,mt)**2 - 2.d0*w0_st * x_s(ms,mt)*x_t(ms,mt) + w0_tt * x_s(ms,mt)**2. )     / xjac**2.

     r0_xx = (r0_ss * y_t(ms,mt)**2 - 2.d0*r0_st * y_s(ms,mt)*y_t(ms,mt) + r0_tt * y_s(ms,mt)**2.     &
	    + r0_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                                  &
	    + r0_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) ) )    / xjac**2.                  &
            - xjac_x * (r0_s* y_t(ms,mt) - r0_t * y_s(ms,mt))  / xjac**2.

     r0_yy = (r0_ss * x_t(ms,mt)**2 - 2.d0*r0_st * x_s(ms,mt)*x_t(ms,mt) + r0_tt * x_s(ms,mt)**2.    &
            + r0_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                              &
	    + r0_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) ) )       / xjac**2.               &	
            - xjac_y * (- r0_s * x_t(ms,mt) + r0_t * x_s(ms,mt) )  / xjac**2.

     r0_xy = (- r0_ss * y_t(ms,mt)*x_t(ms,mt) - r0_tt * x_s(ms,mt)*y_s(ms,mt)                        &
     	      + r0_st * (y_s(ms,mt)*x_t(ms,mt)  + y_t(ms,mt)*x_s(ms,mt)  )                           &
              - r0_s  * (x_st(ms,mt)*y_t(ms,mt) - x_tt(ms,mt)*y_s(ms,mt) )                           &
	          - r0_t * (x_st(ms,mt)*y_s(ms,mt)  - x_ss(ms,mt)*y_t(ms,mt) ) )  / xjac**2.             &	
              - xjac_x * (- r0_s * x_t(ms,mt) + r0_t * x_s(ms,mt) )   / xjac**2.

     T0_xx = (T0_ss * y_t(ms,mt)**2 - 2.d0*T0_st * y_s(ms,mt)*y_t(ms,mt) + T0_tt * y_s(ms,mt)**2.    &
	    + T0_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                                  &
	    + T0_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) ) )    / xjac**2.                  &	
            - xjac_x * (T0_s * y_t(ms,mt) - T0_t * y_s(ms,mt))  / xjac**2.

     T0_yy = (T0_ss * x_t(ms,mt)**2 - 2.d0*T0_st * x_s(ms,mt)*x_t(ms,mt) + T0_tt * x_s(ms,mt)**2.    &
            + T0_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                              &
	        + T0_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) ) )       / xjac**2.           &	
            - xjac_y * (- T0_s * x_t(ms,mt) + T0_t * x_s(ms,mt) )  / xjac**2.

     T0_xy = (- T0_ss * y_t(ms,mt)*x_t(ms,mt) - T0_tt * x_s(ms,mt)*y_s(ms,mt)                        &
     	      + T0_st * (y_s(ms,mt)*x_t(ms,mt)  + y_t(ms,mt)*x_s(ms,mt)  )                           &
              - T0_s  * (x_st(ms,mt)*y_t(ms,mt) - x_tt(ms,mt)*y_s(ms,mt) )                           &
	          - T0_t  * (x_st(ms,mt)*y_s(ms,mt) - x_ss(ms,mt)*y_t(ms,mt) )  )  / xjac**2.            &
              - xjac_x * (- T0_s * x_t(ms,mt) + T0_t * x_s(ms,mt) )   / xjac**2.

     vpar0_xx = (vpar0_ss * y_t(ms,mt)**2 - 2.d0*vpar0_st * y_s(ms,mt)*y_t(ms,mt) + vpar0_tt * y_s(ms,mt)**2. &
               + vpar0_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                                 &
               + vpar0_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) ) )   / xjac**2.                  &
               - xjac_x * (vpar0_s * y_t(ms,mt) - vpar0_t * y_s(ms,mt)) / xjac**2.

     vpar0_yy = (vpar0_ss * x_t(ms,mt)**2 - 2.d0*vpar0_st * x_s(ms,mt)*x_t(ms,mt) + vpar0_tt * x_s(ms,mt)**2. &
               + vpar0_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                                 &
               + vpar0_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) ) )      / xjac**2.               &
               - xjac_y * (- vpar0_s * x_t(ms,mt) + vpar0_t * x_s(ms,mt) ) / xjac**2.

     vpar0_xy = (- vpar0_ss * y_t(ms,mt)*x_t(ms,mt) - vpar0_tt * x_s(ms,mt)*y_s(ms,mt)                       &
                 + vpar0_st * (y_s(ms,mt)*x_t(ms,mt)  + y_t(ms,mt)*x_s(ms,mt)  )                             &
                 - vpar0_s  * (x_st(ms,mt)*y_t(ms,mt) - x_tt(ms,mt)*y_s(ms,mt) )                             &
                 - vpar0_t * (x_st(ms,mt)*y_s(ms,mt)  - x_ss(ms,mt)*y_t(ms,mt) )  )  / xjac**2.              &
                 - xjac_x * (- vpar0_s * x_t(ms,mt) + vpar0_t * x_s(ms,mt) )   / xjac**2.

     P0_xx = r0_xx * T0 + 2.d0 * r0_x * T0_x + r0 * T0_xx
     P0_yy = r0_yy * T0 + 2.d0 * r0_y * T0_y + r0 * T0_yy
     P0_xy = r0_xy * T0 + r0_x * T0_y + r0_y * T0_x + r0 * T0_xy

     T0_ps0_x = T0_xx * ps0_y - T0_xy * ps0_x + T0_x * ps0_xy - T0_y * ps0_xx
     T0_ps0_y = T0_xy * ps0_y - T0_yy * ps0_x + T0_x * ps0_yy - T0_y * ps0_xy

     delta_u_x = (   y_t(ms,mt) * delta_s(mp,2,ms,mt) - y_s(ms,mt) * delta_t(mp,2,ms,mt) ) / xjac
     delta_u_y = ( - x_t(ms,mt) * delta_s(mp,2,ms,mt) + x_s(ms,mt) * delta_t(mp,2,ms,mt) ) / xjac

     delta_ps_x = (   y_t(ms,mt) * delta_s(mp,1,ms,mt) - y_s(ms,mt) * delta_t(mp,1,ms,mt) ) / xjac
     delta_ps_y = ( - x_t(ms,mt) * delta_s(mp,1,ms,mt) + x_s(ms,mt) * delta_t(mp,1,ms,mt) ) / xjac

     ! --- Temperature dependent resistivity
     if ( eta_T_dependent .and. T0 <= T_max_eta) then
       eta_T     = eta   * (abs(T0)/T_0)**(-1.5d0)
       deta_dT   = - eta   * (1.5d0)  * abs(T0)**(-2.5d0) * T_0**(1.5d0)
       d2eta_d2T =   eta   * (3.75d0) * abs(T0)**(-3.5d0) * T_0**(1.5d0)
     else if ( eta_T_dependent .and. T0 > T_max_eta) then
       eta_T     = eta   * (T_max_eta/T_0)**(-1.5d0)
       deta_dT   = 0.
       d2eta_d2T = 0.     
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
      
     ! --- Temperature dependent viscosity
     if ( visco_T_dependent ) then
       visco_T   = visco * (abs(T0)/T_0)**(-1.5d0)
       dvisco_dT = - visco * (1.5d0)  * abs(T0)**(-2.5d0) * T_0**(1.5d0)
       if ( xpoint2 .and. (T0 .lt. T_min) ) then
         visco_T   = visco  * (T_min/T_0)**(-1.5d0)
         dvisco_dT = 0.d0
       endif
     else
       visco_T   = visco
       dvisco_dT = 0.d0
     end if
     
     ! --- Temperature dependent parallel heat diffusivity
     if ( ZKpar_T_dependent ) then
       ZKpar_T   = ZK_par * (abs(T0)/T_0)**(+2.5d0)              ! temperature dependent parallel conductivity
       dZKpar_dT = ZK_par * (2.5d0)  * abs(T0)**(+1.5d0) * T_0**(-2.5d0)
       if (ZKpar_T .gt. ZK_par_max) then
         ZKpar_T   = Zk_par_max
         dZKpar_dT = 0.d0
       endif
       if ( xpoint2 .and. (T0 .lt. T_min) ) then
         ZKpar_T   = ZK_par * (T_min/T_0)**(+2.5d0)
         dZKpar_dT = 0.d0
       endif
     else
       ZKpar_T   = ZK_par                                            ! parallel conductivity
       dZKpar_dT = 0.d0
     endif

     eta_num_T   = eta_num                                           ! hyperresistivity
     visco_num_T = visco_num                                         ! hyperviscosity

     psi_norm = get_psi_n(ps0, y_g(ms,mt))

     D_prof  = get_dperp (psi_norm)
     ZK_prof = get_zkperp(psi_norm)

     if (xpoint2) then
       if (r0 .lt. 0.d0)  then
         D_prof  = D_prof_neg  ! JET : 1.d-4; ITER :  4.d-3
       endif
       if (T0 .lt. 0.d0) then
         ZK_prof = ZK_prof_neg  ! JET : 1.d-3; ITER : 2.d-2 
       endif
     endif
     
     phi       = 2.d0*PI*float(mp-1)/float(n_plane) / float(n_period)
     delta_phi = 2.d0*PI/float(n_plane) / float(n_period)

     source_pellet = 0.d0
     source_volume = 0.d0

     if (use_pellet) then

       call pellet_source2(pellet_amplitude,pellet_R,pellet_Z,pellet_psi,pellet_phi, &
                           pellet_radius, pellet_delta_psi, pellet_sig, pellet_length, pellet_ellipse, pellet_theta, &
                           x_g(ms,mt),y_g(ms,mt), ps0, phi, r0, T0/2.d0, &
                           central_density, pellet_particles, pellet_density, total_pellet_volume, &
                           source_pellet, source_volume)
     endif


     !if (use_pellet) then
  ! 
       !call pellet_source2(pellet_amplitude,pellet_R,pellet_Z,pellet_psi,pellet_phi, &
   !                        pellet_radius, pellet_delta_psi, pellet_sig, pellet_length, &
   !                        x_g(ms,mt),y_g(ms,mt), ps0, phi, eq_zne(ms,mt),eq_zTe(ms,mt), &
   !                        central_density, pellet_particles, pellet_density, total_pellet_volume, &
   !                        source_pellet, source_volume)
   !  endif

     n_tor_local = i_tor_max - i_tor_min + 1     
     do i=1,n_vertex_max

       do j=1,n_degrees

         do im=i_tor_min, i_tor_max

           index_ij = n_tor_local*n_var*n_degrees*(i-1) + n_tor_local * n_var * (j-1) + im - i_tor_min + 1  ! index in the ELM matrix

           v   =  H(i,j,ms,mt) * element%size(i,j) * HZ(im,mp)
           v_x = (  y_t(ms,mt) * h_s(i,j,ms,mt) - y_s(ms,mt) * h_t(i,j,ms,mt) ) * element%size(i,j) / xjac * HZ(im,mp)
           v_y = (- x_t(ms,mt) * h_s(i,j,ms,mt) + x_s(ms,mt) * h_t(i,j,ms,mt) ) * element%size(i,j) / xjac * HZ(im,mp)

           v_s = h_s(i,j,ms,mt) * element%size(i,j) * HZ(im,mp)
           v_t = h_t(i,j,ms,mt) * element%size(i,j) * HZ(im,mp)
           v_p = H(i,j,ms,mt)   * element%size(i,j) * HZ_p(im,mp)

           v_ss = h_ss(i,j,ms,mt) * element%size(i,j) * HZ(im,mp)
           v_tt = h_tt(i,j,ms,mt) * element%size(i,j) * HZ(im,mp)
           v_st = h_st(i,j,ms,mt) * element%size(i,j) * HZ(im,mp)
	   
	       v_xx = (v_ss * y_t(ms,mt)**2. - 2.d0*v_st * y_s(ms,mt)*y_t(ms,mt) + v_tt * y_s(ms,mt)**2. &	        
	            + v_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                           &	   
	            + v_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) ) )  / xjac**2.             &		
	            - xjac_x * (v_s * y_t(ms,mt) - v_t * y_s(ms,mt)) / xjac**2.

	       v_yy = (v_ss * x_t(ms,mt)**2 - 2.d0*v_st * x_s(ms,mt)*x_t(ms,mt) + v_tt * x_s(ms,mt)**2. &	        
		        + v_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                          &	   
	            + v_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) ) )     / xjac**2.         &		
		        - xjac_y * (- v_s * x_t(ms,mt) + v_t * x_s(ms,mt) ) / xjac**2.

           v_xy = (- v_ss * y_t(ms,mt)*x_t(ms,mt) - v_tt * x_s(ms,mt)*y_s(ms,mt)                    &
     	        + v_st * (y_s(ms,mt)*x_t(ms,mt)  + y_t(ms,mt)*x_s(ms,mt)  )                         &        
                - v_s  * (x_st(ms,mt)*y_t(ms,mt) - x_tt(ms,mt)*y_s(ms,mt) )                         &	   
	           - v_t  * (x_st(ms,mt)*y_s(ms,mt) - x_ss(ms,mt)*y_t(ms,mt) )  )  / xjac**2.           &		
                - xjac_x * (- v_s * x_t(ms,mt) + v_t * x_s(ms,mt) )   / xjac**2.   	

           Bgrad_rho_star = ( F0 / BigR * v_p  +  v_x  * ps0_y - v_y  * ps0_x ) / BigR    ! F0 due to absence of normalisation
           Bgrad_rho      = ( F0 / BigR * r0_p +  r0_x * ps0_y - r0_y * ps0_x ) / BigR    ! F0 due to absence of normalisation
           Bgrad_T_star   = ( F0 / BigR * v_p  +  v_x  * ps0_y - v_y  * ps0_x ) / BigR    ! F0 due to absence of normalisation
           Bgrad_T        = ( F0 / BigR * T0_p +  T0_x * ps0_y - T0_y * ps0_x ) / BigR    ! F0 due to absence of normalisation

           BB2 = (F0*F0 + ps0_x * ps0_x + ps0_y * ps0_y )/BigR**2.
           Btheta2=(ps0_x * ps0_x + ps0_y * ps0_y )/BigR
!           Btheta2_psi  = 2.d0 * (psi_x * ps0_x + psi_y * ps0_y ) /BigR**2
           
           v_ps0_x  = v_xx  * ps0_y - v_xy  * ps0_x + v_x  * ps0_xy - v_y * ps0_xx
           v_ps0_y  = v_xy  * ps0_y - v_yy  * ps0_x + v_x  * ps0_yy - v_y * ps0_xy

!###################################################################################################
!#  equation 1   (induction equation)                                                              #
!###################################################################################################

           rhs_ij_1 =   v * eta_T  * (zj0 - current_source(ms,mt)-(jec20+jec10))/ BigR  * xjac * tstep &
                      + v * (ps0_s * u0_t - ps0_t * u0_s)                        * tstep &
                      - v * eps_cyl * F0 / BigR  * u0_p                   * xjac * tstep &

                      - v * tauIC/(r0*BB2) * F0**2./BigR**2. * (ps0_s * p0_t - ps0_t * p0_s) * tstep &		      
                      + v * tauIC/(r0*BB2) * F0**3./BigR**3. * eps_cyl * p0_p * xjac * tstep &

                      + eta_num_T * (v_x * zj0_x + v_y * zj0_y)           * xjac * tstep &
                      + zeta * v * delta_g(mp,1,ms,mt) / BigR             * xjac

!###################################################################################################
!#  equation 2   (perpendicular momentum equation)                                                 #
!###################################################################################################

           rhs_ij_2 = - 0.5d0 * vv2 * (v_x * r0_y_hat - v_y * r0_x_hat)   * xjac * tstep &
                      - r0_hat * BigR**2 * w0 * (v_s * u0_t - v_t * u0_s)        * tstep &
                      + v * (ps0_s * zj0_t - ps0_t * zj0_s )                     * tstep &
                      - visco_T * BigR * (v_x * w0_x + v_y * w0_y)        * xjac * tstep &
                      - v * eps_cyl * F0 / BigR * zj0_p                   * xjac * tstep &
                      + BigR**2 * (v_s * p0_t - v_t * p0_s)                      * tstep &

                      - visco_num_T * (v_xx + v_x/Bigr + v_yy)*(w0_xx + w0_x/Bigr + w0_yy) * xjac * tstep &

                      - TG_num2 * 0.25d0 * r0_hat * BigR**3 * (w0_x * u0_y - w0_y * u0_x) &
                                * ( v_x * u0_y - v_y * u0_x) * xjac * tstep * tstep       &
 			       
!****                      +  tauIC * BigR**4 * w0 * (p0_s * v_t - p0_t * v_s)        * tstep &
!****                      -  4.d0 * tauIC * BigR**3 * v * w0 * p0_y           * xjac * tstep &

                      - v * tauIC * BigR**4 * (p0_s * w0_t - p0_t * w0_s)        * tstep &

		      - tauIC * BigR**3 * p0_y * (v_x* u0_x + v_y * u0_y) * xjac * tstep &
! FO to be verified:
!                      -3*v*tauIC* BigR**3 * (p0_x*u0_xy - p0_y*u0_xx) * xjac * tstep &
		      		      
		      - v * tauIC * BigR**4 * (u0_xy * (p0_xx - p0_yy) - p0_xy * (u0_xx - u0_yy) ) * xjac * tstep &
		      
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
!#  equation 3                                                                                     #
!###################################################################################################

           rhs_ij_3 = - ( v_x * ps0_x  + v_y * ps0_y + v*zj0) / BigR * xjac

!###################################################################################################
!#  equation 4                                                                                     #
!###################################################################################################

           rhs_ij_4 = - ( v_x * u0_x   + v_y * u0_y  + v*w0)  * BigR * xjac 

!###################################################################################################
!#  equation 5 (density equation)                                                                  #
!###################################################################################################

           rhs_ij_5 = v * BigR * (particle_source(ms,mt) + source_pellet)                      * xjac * tstep &
                    + v * BigR**2 * ( r0_s * u0_t - r0_t * u0_s)                                      * tstep &
                    + v * 2.d0 * BigR * r0 * u0_y                                              * xjac * tstep &
                    - (D_par-D_prof) * BigR / BB2 * Bgrad_rho_star * Bgrad_rho                 * xjac * tstep &
                    - D_prof * BigR  * (v_x*r0_x + v_y*r0_y + v_p*r0_p * eps_cyl**2 /BigR**2 ) * xjac * tstep &
                    - v * F0 / BigR * Vpar0 * r0_p                                             * xjac * tstep &
                    - v * Vpar0 * (r0_s * ps0_t - r0_t * ps0_s)                                       * tstep &
                    - v * F0 / BigR * r0 * vpar0_p                                             * xjac * tstep &
                    - v * r0 * (vpar0_s * ps0_t - vpar0_t * ps0_s)                                    * tstep &
		    
		    + v * 2.d0 * tauIC * p0_y * BigR                                           * xjac * tstep &
                    
                    - D_perp_num * (v_xx + v_x/Bigr + v_yy)*(r0_xx + r0_x/Bigr + r0_yy) * BigR * xjac * tstep &  
		                                                                  
                    - TG_num5 * 0.25d0 * BigR**3 * (r0_x * u0_y - r0_y * u0_x)                                &
                                                 * ( v_x * u0_y - v_y * u0_x) * xjac * tstep * tstep          &
                    - TG_num5 * 0.25d0 / BigR * vpar0**2 &
                              * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                              &
                              * ( v_x * ps0_y -  v_y * ps0_x + F0 / BigR * v_p) * xjac * tstep * tstep        &
		    
                    + zeta * v * delta_g(mp,5,ms,mt) * BigR                                    * xjac 

!###################################################################################################
!#  equation 6 (energy  equation)                                                                  #
!###################################################################################################

           rhs_ij_6 = v * BigR * heat_source(ms,mt)                                    * xjac * tstep &

                    + v * r0 * BigR**2 * ( T0_s * u0_t - T0_t * u0_s)                         * tstep &
                    + v * T0 * BigR**2 * ( r0_s * u0_t - r0_t * u0_s)                         * tstep &

                    + v * r0 * 2.d0* GAMMA * BigR * T0 * u0_y                          * xjac * tstep &

                    - v * r0 * F0 / BigR * Vpar0 * T0_p                                * xjac * tstep &
                    - v * T0 * F0 / BigR * Vpar0 * r0_p                                * xjac * tstep &

                    - v * r0 * Vpar0 * (T0_s * ps0_t - T0_t * ps0_s)                          * tstep &
                    - v * T0 * Vpar0 * (r0_s * ps0_t - r0_t * ps0_s)                          * tstep &

                    - v * GAMMA * r0 * T0 * (vpar0_s * ps0_t - vpar0_t * ps0_s)               * tstep &
                    - v * GAMMA * r0 * T0 * F0 / BigR * vpar0_p                        * xjac * tstep &

                    - (ZKpar_T-ZK_prof) * BigR / BB2 * Bgrad_T_star * Bgrad_T          * xjac * tstep &
                    - ZK_prof * BigR * (v_x*T0_x + v_y*T0_y + v_p*T0_p /BigR**2 )      * xjac * tstep &

                    - ZK_perp_num  *  (v_xx + v_x/Bigr + v_yy)*(T0_xx + T0_x/Bigr + T0_yy) * BigR * xjac * tstep &
                     
                    - ZK_par_num * (v_ps0_x  * ps0_y - v_ps0_y  * ps0_x) &
                                 * (T0_ps0_x * ps0_y - T0_ps0_y * ps0_x)               * xjac * tstep &

                    - TG_num6 * 0.25d0 * BigR**3 * T0 * (r0_x * u0_y - r0_y * u0_x) &
                                       * ( v_x * u0_y - v_y * u0_x) * xjac * tstep * tstep  &

                    - TG_num6 * 0.25d0 * BigR**3 * r0 * (T0_x * u0_y - T0_y * u0_x) &
                                       * ( v_x * u0_y - v_y * u0_x) * xjac * tstep * tstep  &

                    - TG_num6 * 0.25d0 / BigR * vpar0**2 &
                              * T0 * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                         &
                              * ( v_x * ps0_y -  v_y * ps0_x + F0 / BigR * v_p) * xjac * tstep * tstep        &

                    - TG_num6 * 0.25d0 / BigR * vpar0**2 &
                              * r0 * (T0_x * ps0_y - T0_y * ps0_x + F0 / BigR * T0_p)                         &
                              * ( v_x * ps0_y -  v_y * ps0_x + F0 / BigR * v_p) * xjac * tstep * tstep        &

                    + zeta * v * r0 * delta_g(mp,6,ms,mt) * BigR                       * xjac &
                    + zeta * v * T0 * delta_g(mp,5,ms,mt) * BigR                       * xjac

!###################################################################################################
!#  equation 7 (parallel velocity  equation)                                                       #
!###################################################################################################

           rhs_ij_7 = - v * F0 / BigR * P0_p                                              * xjac * tstep &
                      - v * (P0_s * ps0_t - P0_t * ps0_s)                                        * tstep &
!==============================================MB:(Vt0_x)beg.=>end:(Vt0_y)======parallel velocity source:
                      - visco_par * (v_x * (vpar0_x-Vt0_x) + v_y * (vpar0_y-Vt0_y)) * BigR* xjac * tstep &
!=================================MB: end

		      		      
		      - v * (particle_source(ms,mt) + source_pellet) * vpar0 * BB2 * BigR * xjac * tstep &
		      
		      - 0.5d0 * r0 * vpar0**2 * BB2 * (ps0_s * v_t - ps0_t * v_s)                * tstep &
		      + 0.5d0 * r0 * vpar0**2 * BB2 * F0 / BigR * v_p                     * xjac * tstep &

                      - 0.5d0 * v * vpar0**2 * BB2 * (ps0_s * r0_t - ps0_t * r0_s)                * tstep &
                      + 0.5d0 * v * vpar0**2 * BB2 * F0 / BigR * r0_p                      * xjac * tstep &

                      - visco_par_num * (v_xx + v_x/Bigr + v_yy)*(vpar0_xx + vpar0_x/Bigr + vpar0_yy) * BigR * xjac * tstep &
		      		      
                      + zeta * v * delta_g(mp,7,ms,mt) * R0 * F0**2 / BigR                * xjac         &
		      + zeta * v * r0 * vpar0 * (ps0_x * delta_ps_x + ps0_y * delta_ps_y) / BigR * xjac  &

                      - TG_NUM7 * 0.25d0 * r0 * Vpar0**2 * BB2 &
                             * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac + F0 / BigR * vpar0_p) / BigR  &
                             * (-(ps0_s * v_t     - ps0_t * v_s)    /xjac + F0 / BigR * v_p)  * xjac * tstep * tstep &

                      - TG_NUM7 * 0.25d0 * v  * Vpar0**2 * BB2 &
                             * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac + F0 / BigR * vpar0_p) / BigR  &
                             * (-(ps0_s * r0_t    - ps0_t * r0_s)   /xjac + F0 / BigR * r0_p)  * xjac * tstep * tstep 

! ------------------------------------------------------ NEO
           if (NEO) then
              rhs_ij_7 =  rhs_ij_7  + amu_neo_prof(ms,mt)*BB2/(Btheta2+epsil)*v*(r0*(ps0_x*u0_x+ps0_y*u0_y)+&
                   tauIC*(ps0_x*P0_x+ps0_y*P0_y) &
                   +aki_neo_prof(ms,mt)*tauIC*r0*(ps0_x*T0_x+ps0_y*T0_y)&
                   -r0*Vpar0*Btheta2)*xjac*tstep*BigR
           endif
! ------------------------------------------------------ NEO

!###################################################################################################
!#  equation 8 (first ECCD current equation)
!###################################################################################################

           Bgrad_jec1    = (F0 / BigR*jec10_p + jec10_x*ps0_y -jec10_y*ps0_x)/BigR

           rhs_ij_8 = - v * BigR * ec_source(ms,mt) *xjac* tstep &
                    - v * jec10 * nu_jec1_fast * BigR* xjac * tstep &
                    + v * JJ_par * BigR/sqrt(BB2) * Bgrad_jec1 * xjac * tstep &
                    + v * zeta * delta_g(mp,8,ms,mt) * BigR * xjac

!###################################################################################################
!#  equation 9 (second ECCD current equation)
!###################################################################################################

           Bgrad_jec2    = (F0 / BigR*jec20_p + jec20_x*ps0_y -jec20_y*ps0_x)/BigR

           rhs_ij_9 = v * BigR * ec_source(ms,mt) *xjac* tstep &
                    - v * jec20 * nu_jec2_fast * BigR* xjac * tstep &
                    + v * JJ_par * BigR/sqrt(BB2) * Bgrad_jec2 * xjac * tstep &
                    + v * zeta * delta_g(mp,9,ms,mt) * BigR * xjac


!###################################################################################################
!#  RHS equations end                                                                                  #
!###################################################################################################

           ij1 = index_ij
           ij2 = index_ij + 1*n_tor_local
           ij3 = index_ij + 2*n_tor_local
           ij4 = index_ij + 3*n_tor_local
           ij5 = index_ij + 4*n_tor_local
           ij6 = index_ij + 5*n_tor_local
           ij7 = index_ij + 6*n_tor_local
           ij8 = index_ij + 7*n_tor_local
           ij9 = index_ij + 8*n_tor_local

           RHS(ij1) = RHS(ij1) + rhs_ij_1 * wst
           RHS(ij2) = RHS(ij2) + rhs_ij_2 * wst
           RHS(ij3) = RHS(ij3) + rhs_ij_3 * wst
           RHS(ij4) = RHS(ij4) + rhs_ij_4 * wst
           RHS(ij5) = RHS(ij5) + rhs_ij_5 * wst
           RHS(ij6) = RHS(ij6) + rhs_ij_6 * wst
           RHS(ij7) = RHS(ij7) + rhs_ij_7 * wst
           RHS(ij8) = RHS(ij8) + rhs_ij_8 * wst
           RHS(ij9) = RHS(ij9) + rhs_ij_9 * wst

           do k=1,n_vertex_max

             do l=1,n_degrees

               do in = i_tor_min, i_tor_max

                 psi   = H(k,l,ms,mt) * element%size(k,l) * HZ(in,mp)

                 psi_x = (   y_t(ms,mt) * h_s(k,l,ms,mt) - y_s(ms,mt) * h_t(k,l,ms,mt) ) / xjac    &
                              * element%size(k,l) * HZ(in,mp)

                 psi_y = ( - x_t(ms,mt) * h_s(k,l,ms,mt) + x_s(ms,mt) * h_t(k,l,ms,mt) )  / xjac   &
                              * element%size(k,l) * HZ(in,mp)

                 psi_p  = H(k,l,ms,mt)    * element%size(k,l) * HZ_p(in,mp)
                 psi_s  = h_s(k,l,ms,mt)  * element%size(k,l) * HZ(in,mp)
                 psi_t  = h_t(k,l,ms,mt)  * element%size(k,l) * HZ(in,mp)
                 psi_ss = h_ss(k,l,ms,mt) * element%size(k,l) * HZ(in,mp)
                 psi_tt = h_tt(k,l,ms,mt) * element%size(k,l) * HZ(in,mp)
                 psi_st = h_st(k,l,ms,mt) * element%size(k,l) * HZ(in,mp)

                 psi_xx = (psi_ss * y_t(ms,mt)**2 - 2.d0*psi_st * y_s(ms,mt)*y_t(ms,mt) + psi_tt * y_s(ms,mt)**2  &	        
		        + psi_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                              &	   
	                + psi_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) ) )    / xjac**2               &		
		        - xjac_x * (psi_s * y_t(ms,mt) - psi_t * y_s(ms,mt)) / xjac**2

	         psi_yy = (psi_ss * x_t(ms,mt)**2 - 2.d0*psi_st * x_s(ms,mt)*x_t(ms,mt) + psi_tt * x_s(ms,mt)**2  &	        
		        + psi_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                              &	   
	                + psi_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) ) )    / xjac**2               &		
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


! additions for the two ECCD current equations
jec1   = psi
jec1_x = psi_x
jec1_y = psi_y
jec1_p = psi_p
jec1_s = psi_s
jec1_t = psi_t

jec2   = psi
jec2_x = psi_x
jec2_y = psi_y
jec2_p = psi_p
jec2_s = psi_s
jec2_t = psi_t
!!

                 u_xx = psi_xx ;                                    rho_xx = psi_xx ;  T_xx = psi_xx ; vpar_xx = psi_xx
                 u_yy = psi_yy ;                                    rho_yy = psi_yy ;  T_yy = psi_yy ; vpar_yy = psi_yy
                 u_xy = psi_xy ;                                    rho_xy = psi_xy ;  T_xy = psi_xy ; vpar_xy = psi_xy

                 w_xx = (psi_ss * y_t(ms,mt)**2 - 2.d0*psi_st * y_s(ms,mt)*y_t(ms,mt) + psi_tt * y_s(ms,mt)**2 ) / xjac**2	        
	         w_yy = (psi_ss * x_t(ms,mt)**2 - 2.d0*psi_st * x_s(ms,mt)*x_t(ms,mt) + psi_tt * x_s(ms,mt)**2 ) / xjac**2
		 w_xy = psi_xy	        

                 rho_hat   = BigR**2 * rho
                 rho_x_hat = 2.d0 * BigR * BigR_x  * rho + BigR**2 * rho_x
                 rho_y_hat = BigR**2 * rho_y
                  Btheta2_psi  = 2.d0 * (psi_x * ps0_x + psi_y * ps0_y ) /BigR**2

                 index_kl = n_tor_local*n_var*n_degrees*(k-1) + n_tor_local * n_var * (l-1) + in - i_tor_min + 1  ! index in the ELM matrix

!###################################################################################################
!#  equation 1   (induction equation)                                                              #
!###################################################################################################

                 amat_11 = v * psi / BigR * xjac * (1.d0 + zeta)                                     &
                         - v * (psi_s * u0_t - psi_t * u0_s)                        * theta * tstep  &

                         + v * tauIC/(r0*BB2) * F0**2/BigR**2 * (psi_s * p0_t - psi_t * p0_s) * theta * tstep 
			 
			 ! term with BB2 still missing
		      

                 amat_12 = - v * (ps0_s * u_t - ps0_t * u_s)                        * theta * tstep  &
                           + eps_cyl * F0 / BigR * v * u_p * xjac                   * theta * tstep

                 amat_13 = - eta_num_T * (v_x * zj_x + v_y * zj_y)           * xjac * theta * tstep  &
                           - eta_T * v * zj / BigR                           * xjac * theta * tstep
			   
	         amat_15 = + v * tauIC/(r0*BB2) * F0**2/BigR**2 * T0  * (ps0_s * rho_t - ps0_t * rho_s) * theta * tstep &
		           + v * tauIC/(r0*BB2) * F0**2/BigR**2 * rho * (ps0_s * T0_t  - ps0_t * T0_s)  * theta * tstep &
			   
                           - v * tauIC/(r0*BB2) * F0**3/BigR**3 * eps_cyl * T0  * rho_p * xjac * theta * tstep &
                           - v * tauIC/(r0*BB2) * F0**3/BigR**3 * eps_cyl * rho * T0_p  * xjac * theta * tstep &
			   
                      - v * tauIC * rho /(r0**2 * BB2) * F0**2/BigR**2 * (ps0_s * p0_t - ps0_t * p0_s) * theta * tstep &		      
                      + v * tauIC * rho /(r0**2 * BB2) * F0**3/BigR**3 * eps_cyl * p0_p * xjac         * theta * tstep 


                amat_16 = - deta_dT * v * T * (zj0 - current_source(ms,mt)-(jec20+jec10)) / BigR * xjac * theta * tstep &
		 
		        + v * tauIC/(r0*BB2) * F0**2/BigR**2 * r0 * (ps0_s * T_t  - ps0_t * T_s) * theta * tstep &
		        + v * tauIC/(r0*BB2) * F0**2/BigR**2 * T  * (ps0_s * r0_t - ps0_t * r0_s)* theta * tstep &
			   
                        - v * tauIC/(r0*BB2) * F0**3./BigR**3. * eps_cyl * r0 * T_p  * xjac * theta * tstep &
                        - v * tauIC/(r0*BB2) * F0**3./BigR**3. * eps_cyl * T  * r0_p * xjac * theta * tstep 


!! additions for the two ECCD current equations
                 amat_18 =  v * eta_T * jec1 / BigR     *xjac * theta * tstep
                 amat_19 =  v * eta_T * jec2 / BigR     *xjac * theta * tstep
!!

!---------------------------------------------------------------- equation 1
!                 amat_11 = v * psi / BigR * xjac * (1.d0 + zeta)                                        &
!                         + eta_T * (psi_x * v_x + psi_y * v_y) / BigR * xjac           * theta * tstep  &
!                         - v * (psi_s * u0_t - psi_t * u0_s)                           * theta * tstep  &
!                         + v * deta_dT * (T0_x * psi_x  + T0_y * psi_y ) / BigR * xjac * theta * tstep

!                 amat_12 = -  v * (ps0_s * u_t - ps0_t * u_s)                 * theta * tstep  &
!                           +  eps_cyl * F0 / BigR * v * u_p * xjac            * theta * tstep

!                 amat_13 = - eta_num * (v_s * zj_t + v_t * zj_s)              * theta * tstep

!                 amat_16 = - deta_dT * v * T * (zj0 - current_source(ms,mt)-jec0) / BigR * xjac * theta * tstep

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

                 amat_22 = - BigR * r0_hat * (v_x * u_x + v_y * u_y) * xjac * (1.d0 + zeta)                                 &
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

                 amat_23 = - v * (ps0_s * zj_t  - ps0_t * zj_s)                           * theta * tstep  &
                           + eps_cyl * F0 / BigR * v * zj_p  * xjac                       * theta * tstep

                 amat_24 = r0_hat * BigR**2 * w  * ( v_s * u0_t - v_t * u0_s)  * theta * tstep  &
                         + BigR * ( v_x * w_x + v_y * w_y) * visco_T  * xjac   * theta * tstep  &
                         + v * tauIC * BigR**4 * (p0_s * w_t - p0_t * w_s)              * theta * tstep &

                         + visco_num_T * (v_xx + v_x/BigR + v_yy)*(w_xx + w_x/BigR + w_yy) * xjac * theta * tstep    &

                         + TG_num2 * 0.25d0 * r0_hat * BigR**3 * (w_x * u0_y - w_y * u0_x)     &
                                   * ( v_x * u0_y - v_y * u0_x) * xjac * theta * tstep * tstep


                 amat_25 = + 0.5d0 * vv2 * (v_x * rho_y_hat - v_y * rho_x_hat)   * xjac * theta * tstep &
                           + rho_hat * BigR**2 * w0 * (v_s * u0_t - v_t * u0_s)         * theta * tstep &
                           - BigR**2 * (v_s * rho_t * T0   - v_t * rho_s * T0  )        * theta * tstep &
                           - BigR**2 * (v_s * rho   * T0_t - v_t * rho   * T0_s)        * theta * tstep &

                           + v * tauIC * BigR**4 * T0  * (rho_s * w0_t - rho_t * w0_s)  * theta * tstep &
                           + v * tauIC * BigR**4 * rho * (T0_s  * w0_t - T0_t  * w0_s)  * theta * tstep &

		           + tauIC * BigR**3 * (T0_y * rho + T0 * rho_y) * (v_x* u0_x + v_y * u0_y) * xjac * theta * tstep &

		           + v * tauIC * BigR**4 * ( (u0_xy * (rho_xx*T0 + 2.d0*rho_x*T0_x + rho*T0_xx           &
			                                    -  rho_yy*T0 - 2.d0*rho_y*T0_y - rho*T0_yy))         &			                  
			                           - (rho_xy * T0 + rho_x*T0_y + rho_y*T0_x + rho*T0_xy) * (u0_xx - u0_yy)  )   &
						 * xjac * theta * tstep                                          &

                           + TG_num2 * 0.25d0 * rho_hat * BigR**3 * (w0_x * u0_y - w0_y * u0_x) &
                                    * ( v_x * u0_y - v_y * u0_x) * xjac * theta * tstep * tstep
! FO to be verified:
!                      +3*v*tauIC* BigR**3 * ((rho*T0_x+rho_x*T0)*u0_xy - (rho*T0_y+rho_y*T0)*u0_xx) * xjac * theta * tstep 

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
                           + dvisco_dT * T * ( v_x * w0_x + v_y * w0_y ) * BigR * xjac  * theta * tstep &

                           + v * tauIC * BigR**4 * r0 * (T_s  * w0_t - T_t  * w0_s)  * theta * tstep &
                           + v * tauIC * BigR**4 * T  * (r0_s * w0_t - r0_t * w0_s)  * theta * tstep &

		           + tauIC * BigR**3 * (r0_y * T + r0 * T_y) * (v_x* u0_x + v_y * u0_y) * xjac * theta * tstep &

		           + v * tauIC * BigR**4 * ( (u0_xy * (T_xx*r0 + 2.d0*T_x*r0_x + T*r0_xx         &
			                                     - T_yy*r0 - 2.d0*T_y*r0_y - T*r0_yy))       &			                  
			                           - (T_xy * r0 + T_x*r0_y + T_y*r0_x + T*r0_xy) * (u0_xx - u0_yy)  )         &
						 * xjac * theta * tstep 
! FO to be verified:
!                          +3*v*tauIC* BigR**3 * ((r0*T_x+r0_x*T)*u0_xy - (r0*T_y+r0_y*T)*u0_xx) * xjac * theta * tstep 

!---------------------------------------- NEO
                 if ( NEO ) then
                    amat_26 = amat_26 -amu_neo_prof(ms,mt)*BB2/((Btheta2+epsil)**2)*(ps0_x*v_x+ps0_y*v_y)*&
                         (tauIC*(ps0_x*(r0_x*T+r0*T_x)+ps0_y*(r0_y*T+r0*T_y)) &
                         +aki_neo_prof(ms,mt)*tauIC*r0*(ps0_x*T_x+ps0_y*T_y))*BigR*xjac*tstep*theta

                    amat_27= amu_neo_prof(ms,mt)*BB2/(Btheta2+epsil)*r0*vpar*(ps0_x*v_x+ps0_y*v_y)&
                         *BigR*xjac*tstep*theta 
                 endif

!---------------------------------------- NEO

!###################################################################################################
!#  equation 3                                                                                     #
!###################################################################################################

                 amat_33 = v * zj / BigR * xjac                             
                 amat_31 = (v_x * psi_x + v_y * psi_y ) / BigR * xjac         

!###################################################################################################
!#  equation 4                                                                                     #
!###################################################################################################

                 amat_44 =  v * w * BigR * xjac                                
                 amat_42 = (v_x * u_x + v_y * u_y) * BigR * xjac               

!###################################################################################################
!#  equation 5    continuity equation (density)                                                    #
!###################################################################################################

                 Bgrad_rho_star_psi = ( v_x  * psi_y - v_y  * psi_x ) / BigR
                 Bgrad_rho_psi      = ( r0_x * psi_y - r0_y * psi_x ) / BigR
                 Bgrad_rho_rho      = ( F0 / BigR * rho_p +  rho_x * ps0_y - rho_y * ps0_x ) / BigR
                 BB2_psi            = 2.d0 * (psi_x * ps0_x + psi_y * ps0_y ) /BigR**2

                 amat_51 = - (D_par-D_prof) * BigR * BB2_psi / BB2**2 * Bgrad_rho_star     * Bgrad_rho     * xjac * theta * tstep &
                           + (D_par-D_prof) * BigR / BB2              * Bgrad_rho_star_psi * Bgrad_rho     * xjac * theta * tstep &
                           + (D_par-D_prof) * BigR / BB2              * Bgrad_rho_star     * Bgrad_rho_psi * xjac * theta * tstep &
                           + v * Vpar0 * (r0_s * psi_t - r0_t * psi_s)                                            * theta * tstep &
                           + v * r0 * (vpar0_s * psi_t - vpar0_t * psi_s)                                         * theta * tstep &

                           + TG_num5 * 0.25d0 / BigR * vpar0**2                                                              &
                                     * (r0_x * psi_y - r0_y * psi_x)                                                         &
                                     * ( v_x * ps0_y -  v_y * ps0_x + F0 / BigR * v_p) * xjac * theta * tstep * tstep        &

                           + TG_num5 * 0.25d0 / BigR * vpar0**2                                                              &
                                     * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                                      &
                                     * ( v_x * psi_y -  v_y * psi_x                   ) * xjac * theta * tstep * tstep


                 amat_52 = - v * BigR**2 * ( r0_s * u_t - r0_t * u_s)                                     * theta * tstep &
                           - v * 2.d0 * BigR * r0 * u_y                                            * xjac * theta * tstep &

                           + TG_num5 * 0.25d0 * BigR**3 * (r0_x * u_y  - r0_y * u_x)                                      &
                                                        * ( v_x * u0_y - v_y  * u0_x) * xjac * theta * tstep * tstep      &

                           + TG_num5 * 0.25d0 * BigR**3 * (r0_x * u0_y - r0_y * u0_x)                                     &
                                                        * ( v_x * u_y  - v_y  * u_x)  * xjac * theta * tstep * tstep 
			   
                 amat_55 = v * rho * BigR * xjac * (1.d0 + zeta)  &
                         - v * BigR**2 * ( rho_s * u0_t - rho_t * u0_s)                                       * theta * tstep &
                         - v * 2.d0 * BigR * rho * u0_y                                                * xjac * theta * tstep &
                         + (D_par-D_prof) * BigR / BB2 * Bgrad_rho_star * Bgrad_rho_rho                * xjac * theta * tstep &
                         + D_prof * BigR  * (v_x*rho_x + v_y*rho_y + v_p*rho_p * eps_cyl**2 /BigR**2 ) * xjac * theta * tstep &
                         + v * F0 / BigR * Vpar0 * rho_p                                               * xjac * theta * tstep &
                         + v * Vpar0 * (rho_s * ps0_t - rho_t * ps0_s)                                        * theta * tstep &
                         + v * rho * (vpar0_s * ps0_t - vpar0_t * ps0_s)                                      * theta * tstep &
                         + v * rho * F0 / BigR * vpar0_p                                               * xjac * theta * tstep &

                         - v * 2.d0 * tauIC * (rho_y * T0 + rho*T0_y) * BigR                           * xjac * theta * tstep &
                                                 
                         + D_perp_num     * (v_xx + v_x/BigR + v_yy)*(rho_xx + rho_x/BigR + rho_yy)   * BigR * xjac * theta * tstep &

                         + TG_num5 * 0.25d0 * BigR**3 * (rho_x * u0_y - rho_y * u0_x)                                &
                                                      * ( v_x  * u0_y - v_y   * u0_x) * xjac * theta * tstep * tstep &

                         + TG_num5 * 0.25d0 / BigR * vpar0**2                                                &
                              * (rho_x * ps0_y - rho_y * ps0_x + F0 / BigR * rho_p)                          &
                              * ( v_x * ps0_y -  v_y * ps0_x   + F0 / BigR * v_p) * xjac * theta * tstep * tstep


		 amat_56 = - v * 2.d0 * tauIC * (T_y * r0 + T*r0_y) * BigR                           * xjac * theta * tstep 

                 amat_57 = + v * F0 / BigR * Vpar * r0_p                                             * xjac * theta * tstep &
                           + v * Vpar * (r0_s * ps0_t - r0_t * ps0_s)                                       * theta * tstep &
                           + v * r0 * (vpar_s * ps0_t - vpar_t * ps0_s)                                     * theta * tstep &
                           + v * r0 * F0 / BigR * vpar_p                                             * xjac * theta * tstep &

                           + TG_num5 * 0.25d0 / BigR * 2.d0*vpar0*vpar                                 &
                              * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                       &
                              * ( v_x * ps0_y -  v_y * ps0_x + F0 / BigR * v_p) * xjac * theta * tstep * tstep 

!###################################################################################################
!#  equation 6   energy equation                                                                   #
!###################################################################################################

                 Bgrad_T_star_psi = ( v_x  * psi_y - v_y  * psi_x  ) / BigR
                 Bgrad_T_psi      = ( T0_x * psi_y - T0_y * psi_x )  / BigR
                 Bgrad_T_T        = ( F0 / BigR * T_p +  T_x * ps0_y - T_y * ps0_x ) / BigR          ! F0 due to absence of normalisation

                 T_ps0_x = T_xx * ps0_y - T_xy * ps0_x + T_x * ps0_xy - T_y * ps0_xx
                 T_ps0_y = T_xy * ps0_y - T_yy * ps0_x + T_x * ps0_yy - T_y * ps0_xy
 
                 T0_psi_x = T0_xx * psi_y - T0_xy * psi_x + T0_x * psi_xy - T0_y * psi_xx
                 T0_psi_y = T0_xy * psi_y - T0_yy * psi_x + T0_x * psi_yy - T0_y * psi_xy
                 
                 v_psi_x = v_xx * psi_y - v_xy * psi_x + v_x * psi_xy - v_y * psi_xx
                 v_psi_y = v_xy * psi_y - v_yy * psi_x + v_x * psi_yy - v_y * psi_xy


                 amat_61 = - (ZKpar_T-ZK_prof) * BigR * BB2_psi / BB2**2 * Bgrad_T_star * Bgrad_T     * xjac * theta * tstep &
                           + (ZKpar_T-ZK_prof) * BigR / BB2     * Bgrad_T_star_psi      * Bgrad_T     * xjac * theta * tstep &
                           + (ZKpar_T-ZK_prof) * BigR / BB2     * Bgrad_T_star          * Bgrad_T_psi * xjac * theta * tstep &
                           + v * r0 * Vpar0 * (T0_s * psi_t - T0_t * psi_s)                                 * theta * tstep &
                           + v * T0 * Vpar0 * (r0_s * psi_t - r0_t * psi_s)                                 * theta * tstep &
                           + v * r0 * GAMMA * T0 * (vpar0_s * psi_t - vpar0_t * psi_s)                      * theta * tstep &

                           + ZK_par_num * (v_psi_x  * ps0_y - v_psi_y  * ps0_x + v_ps0_x * psi_y - v_ps0_y * psi_x)          &
                                        * (T0_ps0_x * ps0_y - T0_ps0_y * ps0_x)                       * xjac * theta * tstep &
                           + ZK_par_num * (T0_psi_x * ps0_y - T0_psi_y * ps0_x + T0_ps0_x * psi_y - T0_ps0_y * psi_x)        &
                                        * (v_ps0_x  * ps0_y - v_ps0_y  * ps0_x)                       * xjac * theta * tstep &

                           + TG_num6 * 0.25d0 / BigR * vpar0**2                                                       &
                                     * T0 * (r0_x * psi_y - r0_y * psi_x)                                             &
                                     * ( v_x * ps0_y -  v_y * ps0_x + F0 / BigR * v_p) * xjac * theta * tstep * tstep &
                     
                           + TG_num6 * 0.25d0 / BigR * vpar0**2                                                       &
                                     * r0 * (T0_x * psi_y - T0_y * psi_x)                                             &
                                     * ( v_x * ps0_y -  v_y * ps0_x + F0 / BigR * v_p) * xjac * theta * tstep * tstep &
                     
                           + TG_num6 * 0.25d0 / BigR * vpar0**2                                                       &
                                     * T0 * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                          &
                                     * ( v_x * psi_y -  v_y * psi_x ) * xjac * theta * tstep * tstep                  &
                     
                           + TG_num6 * 0.25d0 / BigR * vpar0**2                                                       &
                                     * r0 * (T0_x * ps0_y - T0_y * ps0_x + F0 / BigR * T0_p)                          &
                                     * ( v_x * psi_y -  v_y * psi_x ) * xjac * theta * tstep * tstep


                 amat_62 = - v * r0 * BigR**2 * ( T0_s * u_t - T0_t * u_s)                              * theta * tstep &
		           - v * T0 * BigR**2 * ( r0_s * u_t - r0_t * u_s)                              * theta * tstep &
                           - v * r0 * 2.d0* GAMMA * BigR * T0 * u_y                              * xjac * theta * tstep &

                           + TG_num6 * 0.25d0 * BigR**2 * T0* (r0_x * u_y - r0_y * u_x)       &
                                     * ( v_x * u0_y - v_y * u0_x) * xjac * theta*tstep*tstep  &

                           + TG_num6 * 0.25d0 * BigR**2 * r0* (T0_x * u_y - T0_y * u_x)       &
                                     * ( v_x * u0_y - v_y * u0_x) * xjac * theta*tstep*tstep  &

                           + TG_num6 * 0.25d0 * BigR**2 * T0* (r0_x * u0_y - r0_y * u0_x)     &
                                     * ( v_x * u_y - v_y * u_x) * xjac * theta*tstep*tstep    &

                           + TG_num6 * 0.25d0 * BigR**2 * r0* (T0_x * u0_y - T0_y * u0_x)     &
                                     * ( v_x * u_y - v_y * u_x) * xjac * theta*tstep*tstep 


                amat_65 =   v * rho * T0   * BigR * xjac * (1.d0 + zeta)    &

		           - v * rho * BigR**2 * ( T0_s  * u0_t - T0_t  * u0_s)                        * theta * tstep &
		           - v * T0  * BigR**2 * ( rho_s * u0_t - rho_t * u0_s)                        * theta * tstep &

                           - v * rho * 2.d0* GAMMA * BigR * T0 * u0_y                           * xjac * theta * tstep &

                           + v * rho * F0 / BigR * Vpar0 * T0_p                                 * xjac * theta * tstep &
                           + v * T0  * F0 / BigR * Vpar0 * rho_p                                * xjac * theta * tstep &

                           + v * rho * Vpar0 * (T0_s  * ps0_t - T0_t  * ps0_s)                         * theta * tstep &
                           + v * T0  * Vpar0 * (rho_s * ps0_t - rho_t * ps0_s)                         * theta * tstep &

                           + v * rho * GAMMA * T0 * (vpar0_s * ps0_t - vpar0_t * ps0_s)                * theta * tstep &
                           + v * rho * GAMMA * T0 * F0 / BigR * vpar0_p                         * xjac * theta * tstep &

                           + TG_num6 * 0.25d0 * BigR**2 * T0* (rho_x * u0_y - rho_y * u0_x)      &
                                     * ( v_x * u0_y - v_y * u0_x) * xjac * theta*tstep*tstep     &

                           + TG_num6 * 0.25d0 * BigR**2 * rho * (T0_x * u0_y - T0_y * u0_x)      &
                                     * ( v_x * u0_y - v_y * u0_x) * xjac* theta*tstep*tstep      &

                           + TG_num6 * 0.25d0 / BigR * vpar0**2 &
                              * T0 * (rho_x * ps0_y - rho_y * ps0_x + F0 / BigR * rho_p)                      &
                              * ( v_x * ps0_y -  v_y * ps0_x + F0 / BigR * v_p) * xjac * theta * tstep * tstep&

                           + TG_num6 * 0.25d0 / BigR * vpar0**2 &
                              * rho * (T0_x * ps0_y - T0_y * ps0_x + F0 / BigR * T0_p)                        &
                              * ( v_x * ps0_y -  v_y * ps0_x + F0 / BigR * v_p) * xjac * theta * tstep * tstep


                 amat_66 =   v * abs(r0) * T   * BigR * xjac * (1.d0 + zeta)    &

                           - v * r0 * BigR**2 * ( T_s  * u0_t - T_t  * u0_s)                       * theta * tstep &
                           - v * T  * BigR**2 * ( r0_s * u0_t - r0_t * u0_s)                       * theta * tstep &

                           - v * r0 * 2.d0* GAMMA * BigR * T * u0_y                         * xjac * theta * tstep &

                           + v * r0 * F0 / BigR * Vpar0 * T_p                               * xjac * theta * tstep &
                           + v * T * F0  / BigR * Vpar0 * r0_p                              * xjac * theta * tstep &

                           + v * r0 * Vpar0 * (T_s  * ps0_t - T_t  * ps0_s)                        * theta * tstep &
                           + v * T  * Vpar0 * (r0_s * ps0_t - r0_t * ps0_s)                        * theta * tstep &

                           + v * r0 * GAMMA * T * (vpar0_s * ps0_t - vpar0_t * ps0_s)              * theta * tstep &
                           + v * r0 * GAMMA * T * F0 / BigR * vpar0_p                       * xjac * theta * tstep &

                           + (ZKpar_T-ZK_prof) * BigR / BB2 * Bgrad_T_star * Bgrad_T_T      * xjac * theta * tstep &
                           + ZK_prof * BigR * (v_x*T_x + v_y*T_y + v_p*T_p /BigR**2 )       * xjac * theta * tstep &

                           + dZKpar_dT * T * BigR / BB2 * Bgrad_T_star * Bgrad_T            * xjac * theta * tstep &

                           + ZK_perp_num * (v_xx + v_x/BigR + v_yy)*(T_xx + T_x/BigR + T_yy) * BigR * xjac * theta * tstep &

                           + TG_num6 * 0.25d0 * BigR**2 * T* (r0_x * u0_y - r0_y * u0_x)         &
                                     * ( v_x * u0_y - v_y * u0_x) * xjac * theta * tstep * tstep &

                           + TG_num6 * 0.25d0 * BigR**2 * r0* (T_x * u0_y - T_y * u0_x)          &
                                     * ( v_x * u0_y - v_y * u0_x) * xjac * theta * tstep * tstep &

                           + TG_num6 * 0.25d0 / BigR * vpar0**2 &
                              * T * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                          &
                              * ( v_x * ps0_y -  v_y * ps0_x + F0 / BigR * v_p) * xjac * theta * tstep * tstep&

                           + TG_num6 * 0.25d0 / BigR * vpar0**2 &
                              * r0 * (T_x * ps0_y - T_y * ps0_x + F0 / BigR * T_p)                            &
                              * ( v_x * ps0_y -  v_y * ps0_x + F0 / BigR * v_p) * xjac * theta * tstep * tstep


                 amat_67 = + v * r0 * F0 / BigR * Vpar * T0_p                               * xjac * theta * tstep &
		           + v * T0 * F0 / BigR * Vpar * r0_p                               * xjac * theta * tstep &

                           + v * r0 * Vpar * (T0_s * ps0_t - T0_t * ps0_s)                         * theta * tstep &
                           + v * T0 * Vpar * (r0_s * ps0_t - r0_t * ps0_s)                         * theta * tstep &

                           + v * r0 * GAMMA * T0 * (vpar_s * ps0_t - vpar_t * ps0_s)               * theta * tstep &
                           + v * r0 * GAMMA * T0 * F0 / BigR * vpar_p                       * xjac * theta * tstep &

                           + TG_num6 * 0.25d0 / BigR * 2.d0 * vpar0*vpar &
                              * T0 * (r0_x * ps0_y - r0_y * ps0_x + F0 / BigR * r0_p)                          &
                              * ( v_x * ps0_y -  v_y * ps0_x + F0 / BigR * v_p) * xjac * theta * tstep * tstep &

                           + TG_num6 * 0.25d0 / BigR * 2.d0 * vpar0*vpar &
                              * r0 * (T0_x * ps0_y - T0_y * ps0_x + F0 / BigR * T0_p)                          &
                              * ( v_x * ps0_y -  v_y * ps0_x + F0 / BigR * v_p) * xjac * theta * tstep * tstep 


!###################################################################################################
!#  equation 7   parallel velocity equation                                                        #
!###################################################################################################

                 amat_71 = v * r0 * vpar0 / BigR * (ps0_x * psi_x + ps0_y * psi_y) * xjac * (1.d0 + zeta) &
		         
			 + v * (P0_s * psi_t - P0_t * psi_s)                                            * theta * tstep &

                         + 0.5d0 * r0 * vpar0**2 * BB2 * (psi_s * v_t - psi_t * v_s)                    * theta * tstep &
                         + 0.5d0 * v  * vpar0**2 * BB2 * (psi_s * r0_t - psi_t * r0_s)                  * theta * tstep &

                         + v * (particle_source(ms,mt) + source_pellet) * vpar0 * BB2_psi * BigR * xjac * theta * tstep &

                         + TG_NUM7 * 0.25d0 * r0 * Vpar0**2 * BB2 &
                                   * (-(psi_s * vpar0_t - psi_t * vpar0_s)/xjac) / BigR  &
                                   * (-(psi_s * v_t     - psi_t * v_s)    /xjac)  * xjac * theta * tstep*tstep &
             
                         + TG_NUM7 * 0.25d0 * v  * Vpar0**2 * BB2 &
                                   * (-(psi_s * vpar0_t - psi_t * vpar0_s)/xjac) / BigR  &
                                   * (-(psi_s * r0_t    - psi_t * r0_s)   /xjac)  * xjac * theta * tstep*tstep 

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
                           + v * F0 / BigR * (rho_p * T0 + rho * T0_p)              * xjac * theta * tstep &

		           + 0.5d0 * rho * vpar0**2 * BB2 * (ps0_s * v_t - ps0_t * v_s)    * theta * tstep &
		           - 0.5d0 * rho * vpar0**2 * BB2 * F0 / BigR * v_p         * xjac * theta * tstep &
                           + 0.5d0 * v   * vpar0**2 * BB2 * (ps0_s * rho_t - ps0_t * rho_s)* theta * tstep &
                           - 0.5d0 * v   * vpar0**2 * BB2 * F0 / BigR * rho_p       * xjac * theta * tstep &

                           + TG_NUM7 * 0.25d0 * rho * Vpar0**2 * BB2 &
                                     * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac + F0 / BigR * vpar0_p) / BigR  &
                                     * (-(ps0_s * v_t     - ps0_t * v_s)    /xjac + F0 / BigR * v_p)  * xjac * theta * tstep*tstep &
               
                           + TG_NUM7 * 0.25d0 * v * Vpar0**2 * BB2 &
                                     * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac + F0 / BigR * vpar0_p) / BigR  &
                                     * (-(ps0_s * rho_t   - ps0_t * rho_s)  /xjac + F0 / BigR * rho_p)* xjac * theta * tstep*tstep

                 amat_76 = + v * (T_s * r0 * ps0_t - T_t * r0 * ps0_s)                     * theta * tstep &
                           + v * (T * r0_s * ps0_t - T * r0_t * ps0_s)                     * theta * tstep &
                           + v * F0 / BigR * (T_p * r0 + T * r0_p)                  * xjac * theta * tstep

                 amat_77 = v * Vpar * abs(r0) * F0**2 / BigR * xjac * (1.d0 + zeta) &
                         + visco_par * (v_x * Vpar_x + v_y * Vpar_y) * BigR        * xjac  * theta * tstep &

                         + v * (particle_source(ms,mt) + source_pellet)*vpar*BB2 * BigR  * xjac * theta * tstep &
		      
                         + r0 * vpar0 * vpar * BB2 * (ps0_s * v_t - ps0_t * v_s)                * theta * tstep &
                         - r0 * vpar0 * vpar * BB2 * F0 / BigR * v_p                     * xjac * theta * tstep &
                         + v  * vpar0 * vpar * BB2 * (ps0_s * r0_t - ps0_t * r0_s)              * theta * tstep &
                         - v  * vpar0 * vpar * BB2 * F0 / BigR * r0_p                    * xjac * theta * tstep &
                     
                         + visco_par_num * (v_xx + v_x/BigR + v_yy)*(vpar_xx + vpar_x/BigR + vpar_yy) * BigR * xjac * theta * tstep &

                         + TG_NUM7 * 0.5d0 * r0 * Vpar * Vpar0 * BB2 &
                                   * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac + F0 / BigR * vpar0_p) / BigR                     &
                                   * (-(ps0_s * v_t     - ps0_t * v_s)    /xjac + F0 / BigR * v_p)  * xjac * theta * tstep*tstep  &
             
                         + TG_NUM7 * 0.5d0 * v * Vpar * Vpar0 * BB2 &
                                   * (-(ps0_s * vpar0_t - ps0_t * vpar0_s)/xjac + F0 / BigR * vpar0_p) / BigR                     &
                                   * (-(ps0_s * r0_t    - ps0_t * r0_s)   /xjac + F0 / BigR * r0_p)  * xjac * theta * tstep*tstep &
             
                         + TG_NUM7 * 0.25d0 * r0 * Vpar0**2 * BB2 &
                                   * (-(ps0_s * vpar_t - ps0_t * vpar_s)/xjac + F0 / BigR * vpar_p) / BigR                        &
                                   * (-(ps0_s * v_t    - ps0_t * v_s)   /xjac + F0 / BigR * v_p)  * xjac * theta * tstep*tstep    &
             
                         + TG_NUM7 * 0.25d0 * v * Vpar0**2 * BB2 &
                                   * (-(ps0_s * vpar_t - ps0_t * vpar_s)/xjac + F0 / BigR * vpar_p) / BigR                        &
                                   * (-(ps0_s * r0_t   - ps0_t * r0_s)  /xjac + F0 / BigR * r0_p)  * xjac * theta * tstep*tstep 

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
!#  equation 8 first ECCD current
!###################################################################################################

       Bgrad_jec1_psi      = ( jec10_x * psi_y - jec10_y * psi_x )  /BigR
       Bgrad_jec1_jec      = ( F0 / BigR * jec1_p +  jec1_x * ps0_y -jec1_y *ps0_x ) / BigR

       amat_81 =  -JJ_par*BigR/sqrt(BB2)       *v*Bgrad_jec1_psi*xjac*theta*tstep &
       +JJ_par*BigR*sqrt(BB2)*BB2_psi /(2.*BB2**2)  *v *Bgrad_jec1*xjac*theta*tstep

       amat_88 =  v * jec1 * BigR * xjac * (1.d0 + zeta) &
                + v * nu_jec1_fast * BigR * jec1 * xjac * theta * tstep &
                - v * JJ_par*BigR/sqrt(BB2)* Bgrad_jec1_jec * xjac * theta *tstep

!###################################################################################################
!#  equation 9 second ECCD current
!###################################################################################################

       Bgrad_jec2_psi      = ( jec20_x * psi_y - jec20_y * psi_x )  /BigR
       Bgrad_jec2_jec      = ( F0 / BigR * jec2_p +  jec2_x * ps0_y -jec2_y *ps0_x )/ BigR

       amat_91 =  -JJ_par*BigR/sqrt(BB2)       *v*Bgrad_jec2_psi*xjac*theta*tstep  &
       +JJ_par*BigR*sqrt(BB2)*BB2_psi /(2.*BB2**2)  *v*Bgrad_jec2*xjac*theta*tstep

       amat_99 =  v * jec2 * BigR * xjac * (1.d0 + zeta) &
                + v * nu_jec2_fast * BigR * jec2 * xjac * theta * tstep &
                - v * JJ_par*BigR/sqrt(BB2)* Bgrad_jec2_jec * xjac * theta *tstep

!###################################################################################################
!#  end LHS equations                                                                              #
!###################################################################################################


                 kl1 = index_kl
                 kl2 = index_kl + 1*n_tor_local
                 kl3 = index_kl + 2*n_tor_local
                 kl4 = index_kl + 3*n_tor_local
                 kl5 = index_kl + 4*n_tor_local
                 kl6 = index_kl + 5*n_tor_local
                 kl7 = index_kl + 6*n_tor_local
                 kl8 = index_kl + 7*n_tor_local
                 kl9 = index_kl + 8*n_tor_local

                 ELM(ij1,kl1) =  ELM(ij1,kl1) + wst * amat_11
                 ELM(ij1,kl2) =  ELM(ij1,kl2) + wst * amat_12
                 ELM(ij1,kl3) =  ELM(ij1,kl3) + wst * amat_13
                 ELM(ij1,kl5) =  ELM(ij1,kl5) + wst * amat_15
                 ELM(ij1,kl6) =  ELM(ij1,kl6) + wst * amat_16
                 ELM(ij1,kl8) =  ELM(ij1,kl8) + wst * amat_18
                 ELM(ij1,kl9) =  ELM(ij1,kl9) + wst * amat_19

                 ELM(ij2,kl1) =  ELM(ij2,kl1) + wst * amat_21
                 ELM(ij2,kl2) =  ELM(ij2,kl2) + wst * amat_22
                 ELM(ij2,kl3) =  ELM(ij2,kl3) + wst * amat_23
                 ELM(ij2,kl4) =  ELM(ij2,kl4) + wst * amat_24
                 ELM(ij2,kl5) =  ELM(ij2,kl5) + wst * amat_25
                 ELM(ij2,kl6) =  ELM(ij2,kl6) + wst * amat_26
!---------------------------------------- NEO
                 ELM(ij2,kl7) =  ELM(ij2,kl7) + wst * amat_27
!---------------------------------------- NEO

                 ELM(ij3,kl1) =  ELM(ij3,kl1) + wst * amat_31
                 ELM(ij3,kl3) =  ELM(ij3,kl3) + wst * amat_33

                 ELM(ij4,kl2) =  ELM(ij4,kl2) + wst * amat_42
                 ELM(ij4,kl4) =  ELM(ij4,kl4) + wst * amat_44

                 ELM(ij5,kl1) =  ELM(ij5,kl1) + wst * amat_51
                 ELM(ij5,kl2) =  ELM(ij5,kl2) + wst * amat_52
                 ELM(ij5,kl5) =  ELM(ij5,kl5) + wst * amat_55
                 ELM(ij5,kl6) =  ELM(ij5,kl6) + wst * amat_56
                 ELM(ij5,kl7) =  ELM(ij5,kl7) + wst * amat_57

                 ELM(ij6,kl1) =  ELM(ij6,kl1) + wst * amat_61
                 ELM(ij6,kl2) =  ELM(ij6,kl2) + wst * amat_62
                 ELM(ij6,kl5) =  ELM(ij6,kl5) + wst * amat_65
                 ELM(ij6,kl6) =  ELM(ij6,kl6) + wst * amat_66
                 ELM(ij6,kl7) =  ELM(ij6,kl7) + wst * amat_67

                 ELM(ij7,kl1) =  ELM(ij7,kl1) + wst * amat_71
                 ELM(ij7,kl2) =  ELM(ij7,kl2) + wst * amat_72
                 ELM(ij7,kl5) =  ELM(ij7,kl5) + wst * amat_75
                 ELM(ij7,kl6) =  ELM(ij7,kl6) + wst * amat_76
                 ELM(ij7,kl7) =  ELM(ij7,kl7) + wst * amat_77

                 ELM(ij8,kl1) =  ELM(ij8,kl1) + wst * amat_81
                 ELM(ij8,kl8) =  ELM(ij8,kl8) + wst * amat_88

                 ELM(ij9,kl1) =  ELM(ij9,kl1) + wst * amat_91
                 ELM(ij9,kl9) =  ELM(ij9,kl9) + wst * amat_99

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

subroutine eccurrent_source(R,Z,zjz)
!-----------------------------------------------------------------------
! routine to calculate the ECCD injected at a given position R,Z
! Currently this injects the current in a radial annulus in the poloidal plane.
! Eventually a more sophisticated current injection is planned.
!-----------------------------------------------------------------------
use phys_module

implicit none

real*8  :: R,Z,zjz
real*8  :: rad,r1,z1

zjz=0.d0

!calculating the radius of the poloidal cross-section
r1=R-R_geo
rad=r1**2.+Z**2.
rad=sqrt(rad)

! positioning a blob
!zjz=exp(-.5*(r1-jec_pos1)**2./jec_width**2.)*&
!    exp(-.5*(Z-jec_pos2)**2./jec_width2**2.) +&
!    exp(-.5*(r1-jec_pos3)**2./jec_width**2.)*&
!    exp(-.5*(Z-jec_pos4)**2./jec_width2**2.)
!
!zjz=zjz*jecamp

! setting the applied ECCD current normally-distributed on a radial
! annulus of the poloidal cross-section.
zjz=jecamp*exp(-.5*(rad-jec_pos1)**2./jec_width**2.)
!zjz=1.d0

! for debugging purposes, a solid annulus of current
!if(rad.gt.jw1.AND.rad.lt.jw2) then
!zjz=jecamp*(rad-jw1)/(jw2-jw1)
!elseif(rad.gt.jw2.AND.rad.lt.jw3) then
!zjz=jecamp*(rad-jw3)/(jw2-jw3)
!endif

!zjz=jecamp

return
end subroutine eccurrent_source

end module mod_elt_matrix
