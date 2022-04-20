module mod_elt_matrix
  implicit none
contains

subroutine element_matrix(element,nodes, xpoint2, xcase2, R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint, Z_xpoint, ELM, RHS, tid, i_tor_min, i_tor_max)
!---------------------------------------------------------------
! calculates the matrix contribution of one element
!---------------------------------------------------------------
use mod_parameters
use data_structure
use gauss
use basis_at_gaussian
use phys_module

implicit none

type (type_element)   :: element
type (type_node)      :: nodes(n_vertex_max)

real*8, dimension (:,:), allocatable  :: ELM
real*8, dimension (:)  , allocatable  :: RHS
integer, intent(in)                   :: tid, i_tor_min, i_tor_max

integer    :: i, j, ms, mt, mp, k, l, index_ij, index_kl, index, xcase2
integer    :: in, im, ij1, ij2, kl1, kl2
real*8     :: wst, xjac, xjac_x, xjac_y
real*8     :: R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint(2), Z_xpoint(2), psi_norm
real*8     :: rhs_ij_1,   rhs_ij_2
real*8     :: theta, zeta, delta_u_x, delta_u_y

real*8     :: v,  v_x,  v_y,  v_s,  v_t,  v_ss,  v_st,  v_tt,  v_xx,  v_yy
real*8     :: u0, u0_x, u0_y, u0_s, u0_t, u0_ss, u0_st, u0_tt
real*8     :: w0, w0_x, w0_y, w0_s, w0_t, w0_ss, w0_st, w0_tt, w0_xx, w0_yy, w0_xy
real*8     :: u,  u_x,  u_y,  u_s,  u_t,  u_ss,  u_st,  u_tt
real*8     :: w,  w_x,  w_y,  w_s,  w_t,  w_ss,  w_st,  w_tt,  w_xx, w_yy
real*8     :: amat_11, amat_12, amat_21, amat_22

real*8     :: R,Z,WW, WW_R, WW_Z, WW_RR, WW_ZZ, WW_RZ
logical    :: xpoint2

real*8, dimension(n_gauss,n_gauss)    :: x_g, x_s, x_t
real*8, dimension(n_gauss,n_gauss)    :: x_ss, x_st, x_tt
real*8, dimension(n_gauss,n_gauss)    :: y_g, y_s, y_t
real*8, dimension(n_gauss,n_gauss)    :: y_ss, y_st, y_tt

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

! --- Take time evolution parameters from phys_module
theta = time_evol_theta
!zeta  = time_evol_zeta
! change zeta for variable dt
zeta  = time_evol_zeta * 2.0d0 * tstep / (tstep + tstep_prev)

!---------------------------------------------------- value of (x,y) and derivatives on Gaussian points
x_g  = 0.d0; x_s  = 0.d0; x_t  = 0.d0; x_ss  = 0.d0; x_st  = 0.d0; x_tt  = 0.d0;
y_g  = 0.d0; y_s  = 0.d0; y_t  = 0.d0; y_ss  = 0.d0; y_st  = 0.d0; y_tt  = 0.d0;
eq_g = 0.d0; eq_s = 0.d0; eq_t = 0.d0; eq_ss = 0.d0; eq_st = 0.d0; eq_tt = 0.d0;

delta_g = 0.d0; delta_s = 0.d0; delta_t = 0.d0

mp = 1
in = 1

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

        do k=1,n_var

          eq_g(mp,k,ms,mt) = eq_g(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H(i,j,ms,mt)  
          eq_s(mp,k,ms,mt) = eq_s(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_s(i,j,ms,mt)
          eq_t(mp,k,ms,mt) = eq_t(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_t(i,j,ms,mt)

          eq_ss(mp,k,ms,mt) = eq_ss(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_ss(i,j,ms,mt)
          eq_st(mp,k,ms,mt) = eq_st(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_st(i,j,ms,mt)
          eq_tt(mp,k,ms,mt) = eq_tt(mp,k,ms,mt) + nodes(i)%values(in,j,k) * element%size(i,j) * H_tt(i,j,ms,mt)

          delta_g(mp,k,ms,mt) = delta_g(mp,k,ms,mt) + nodes(i)%deltas(in,j,k) * element%size(i,j) * H(i,j,ms,mt)  
          delta_s(mp,k,ms,mt) = delta_s(mp,k,ms,mt) + nodes(i)%deltas(in,j,k) * element%size(i,j) * H_s(i,j,ms,mt)
          delta_t(mp,k,ms,mt) = delta_t(mp,k,ms,mt) + nodes(i)%deltas(in,j,k) * element%size(i,j) * H_t(i,j,ms,mt) 

        enddo

      enddo
    enddo 
      
  enddo
enddo

! changes deltas for variable time steps
delta_g = delta_g * tstep / tstep_prev
delta_s = delta_s * tstep / tstep_prev
delta_t = delta_t * tstep / tstep_prev

!--------------------------------------------------- sum over the Gaussian integration points
do ms=1, n_gauss

  do mt=1, n_gauss

    wst = wgauss(ms)*wgauss(mt)

    xjac    = x_s(ms,mt)*y_t(ms,mt)  - x_t(ms,mt)*y_s(ms,mt)
    
    xjac_x  = (x_ss(ms,mt)*y_t(ms,mt)**2 - y_ss(ms,mt)*x_t(ms,mt)*y_t(ms,mt) - 2.d0*x_st(ms,mt)*y_s(ms,mt)*y_t(ms,mt)  &           
	    + y_st(ms,mt)*(x_s(ms,mt)*y_t(ms,mt) + x_t(ms,mt)*y_s(ms,mt))                                              &
	    + x_tt(ms,mt)*y_s(ms,mt)**2 - y_tt(ms,mt)*x_s(ms,mt)*y_s(ms,mt)) / xjac
	   
    xjac_y  = (y_tt(ms,mt)*x_s(ms,mt)**2 - x_tt(ms,mt)*y_s(ms,mt)*x_s(ms,mt) - 2.d0*y_st(ms,mt)*x_t(ms,mt)*x_s(ms,mt)  &           
            + x_st(ms,mt)*(y_t(ms,mt)*x_s(ms,mt) + y_s(ms,mt)*x_t(ms,mt))                                              &
            + y_ss(ms,mt)*x_t(ms,mt)**2 - x_ss(ms,mt)*y_t(ms,mt)*x_t(ms,mt)) / xjac

    u0    = eq_g(mp,1,ms,mt)
    u0_x  = (   y_t(ms,mt) * eq_s(mp,1,ms,mt) - y_s(ms,mt) * eq_t(mp,1,ms,mt) ) / xjac
    u0_y  = ( - x_t(ms,mt) * eq_s(mp,1,ms,mt) + x_s(ms,mt) * eq_t(mp,1,ms,mt) ) / xjac
    u0_s  = eq_s(mp,1,ms,mt)
    u0_t  = eq_t(mp,1,ms,mt)

    w0    = eq_g(mp,2,ms,mt)
    w0_x  = (   y_t(ms,mt) * eq_s(mp,2,ms,mt) - y_s(ms,mt) * eq_t(mp,2,ms,mt) ) / xjac
    w0_y  = ( - x_t(ms,mt) * eq_s(mp,2,ms,mt) + x_s(ms,mt) * eq_t(mp,2,ms,mt) ) / xjac
    w0_s  = eq_s(mp,2,ms,mt)
    w0_t  = eq_t(mp,2,ms,mt)
    w0_ss = eq_ss(mp,2,ms,mt)
    w0_tt = eq_tt(mp,2,ms,mt)
    w0_st = eq_st(mp,2,ms,mt)

    w0_xx = (w0_ss * y_t(ms,mt)**2 - 2.d0 * w0_st * y_s(ms,mt)*y_t(ms,mt) + w0_tt * y_s(ms,mt)**2  &             
           + w0_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                             &        
           + w0_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) ) )    / xjac**2              &             
           - xjac_x * (w0_s* y_t(ms,mt) - w0_t * y_s(ms,mt))  / xjac**2

    w0_yy = (w0_ss * x_t(ms,mt)**2 - 2.d0 * w0_st * x_s(ms,mt)*x_t(ms,mt) + w0_tt * x_s(ms,mt)**2  &             
           + w0_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                             &        
           + w0_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) ) )       / xjac**2           &             
           - xjac_y * (- w0_s * x_t(ms,mt) + w0_t * x_s(ms,mt) )  / xjac**2

    w0_xy = (- w0_ss * y_t(ms,mt)*x_t(ms,mt) - w0_tt * x_s(ms,mt)*y_s(ms,mt)                     &
     	      + w0_st * (y_s(ms,mt)*x_t(ms,mt)  + y_t(ms,mt)*x_s(ms,mt)  )                        &        
              - w0_s  * (x_st(ms,mt)*y_t(ms,mt) - x_tt(ms,mt)*y_s(ms,mt) )                        &	   
	      - w0_t * (x_st(ms,mt)*y_s(ms,mt)  - x_ss(ms,mt)*y_t(ms,mt) )  )  / xjac**2          &		
           - xjac_x * (- w0_s * x_t(ms,mt) + w0_t * x_s(ms,mt) )   / xjac**2
	   
	       
    delta_u_x = (   y_t(ms,mt) * delta_s(mp,1,ms,mt) - y_s(ms,mt) * delta_t(mp,1,ms,mt) ) / xjac
    delta_u_y = ( - x_t(ms,mt) * delta_s(mp,1,ms,mt) + x_s(ms,mt) * delta_t(mp,1,ms,mt) ) / xjac

    do i=1,n_vertex_max

      do j=1,n_degrees

        index_ij = n_var*n_degrees*(i-1) + n_var * (j-1) + 1   ! index in the ELM matrix

        v   =  H(i,j,ms,mt) * element%size(i,j) 
        v_x = (  y_t(ms,mt) * h_s(i,j,ms,mt) - y_s(ms,mt) * h_t(i,j,ms,mt) ) * element%size(i,j) / xjac 
        v_y = (- x_t(ms,mt) * h_s(i,j,ms,mt) + x_s(ms,mt) * h_t(i,j,ms,mt) ) * element%size(i,j) / xjac 

        v_s  = h_s(i,j,ms,mt)  * element%size(i,j) 
        v_t  = h_t(i,j,ms,mt)  * element%size(i,j)
        v_ss = h_ss(i,j,ms,mt) * element%size(i,j) 
        v_tt = h_tt(i,j,ms,mt) * element%size(i,j) 
        v_st = h_st(i,j,ms,mt) * element%size(i,j)

	v_xx = (v_ss * y_t(ms,mt)**2 - 2.d0*v_st * y_s(ms,mt)*y_t(ms,mt) + v_tt * y_s(ms,mt)**2   &	        
	      + v_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                          &	   
	      + v_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) ) )  / xjac**2             &		
	      - xjac_x * (v_s * y_t(ms,mt) - v_t * y_s(ms,mt)) / xjac**2

	v_yy = (v_ss * x_t(ms,mt)**2 - 2.d0*v_st * x_s(ms,mt)*x_t(ms,mt) + v_tt * x_s(ms,mt)**2   &	        
	      + v_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                          &	   
	      + v_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) ) )     / xjac**2          &		
	      - xjac_y * (- v_s * x_t(ms,mt) + v_t * x_s(ms,mt) ) / xjac**2


        rhs_ij_1 = - w0    * (v_s * u0_t - v_t * u0_s)                   * tstep &
                   - visco * (v_x * w0_x + v_y * w0_y)            * xjac * tstep &
                   - visco_num * (v_xx + v_yy)*(w0_xx + w0_yy)    * xjac * tstep & 
		   
		   - 0.25d0 * (w0_x * u0_y - w0_y * u0_x) * ( v_x * u0_y - v_y * u0_x) * xjac * tstep * tstep &
		   
                   - zeta  * (v_x * delta_u_x + v_y * delta_u_y)  * xjac  
           
        rhs_ij_2 = 0.d0 

        ij1 = index_ij
        ij2 = index_ij + 1

        RHS(ij1) = RHS(ij1) + rhs_ij_1 * wst
        RHS(ij2) = RHS(ij2) + rhs_ij_2 * wst

        do k=1,n_vertex_max

          do l=1,n_degrees

            u    = h(k,l,ms,mt)    * element%size(k,l)
            u_s  = h_s(k,l,ms,mt)  * element%size(k,l)
            u_t  = h_t(k,l,ms,mt)  * element%size(k,l) 
            u_ss = h_ss(k,l,ms,mt) * element%size(k,l) 
            u_tt = h_tt(k,l,ms,mt) * element%size(k,l) 
            u_st = h_st(k,l,ms,mt) * element%size(k,l)
		 
            u_x = (   y_t(ms,mt) * h_s(k,l,ms,mt) - y_s(ms,mt) * h_t(k,l,ms,mt) ) / xjac * element%size(k,l) 
            u_y = ( - x_t(ms,mt) * h_s(k,l,ms,mt) + x_s(ms,mt) * h_t(k,l,ms,mt) ) / xjac * element%size(k,l)

            w = u ; w_x = u_x ; w_y = u_y ; w_s = u_s ; w_t = u_t; w_ss = u_ss ; w_st = u_st ; w_tt = u_tt 

            w_xx = (w_ss * y_t(ms,mt)**2 - 2.d0*w_st * y_s(ms,mt)*y_t(ms,mt) + w_tt * y_s(ms,mt)**2  &	        
		  + w_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                         &	   
	          + w_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) ) )    / xjac**2          &		
		  - xjac_x * (w_s * y_t(ms,mt) - w_t * y_s(ms,mt)) / xjac**2

	    w_yy = (w_ss * x_t(ms,mt)**2 - 2.d0*w_st * x_s(ms,mt)*x_t(ms,mt) + w_tt * x_s(ms,mt)**2  &	        
		  + w_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                         &	   
	          + w_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) ) )    / xjac**2          &		
		  - xjac_y * (- w_s * x_t(ms,mt) + w_t * x_s(ms,mt) ) / xjac**2
		
            index_kl = n_var*n_degrees*(k-1) + n_var * (l-1) + 1   ! index in the ELM matrix

!---------------------------------------------------------------- equation 1
 		      
            amat_11 = - (v_x * u_x + v_y * u_y) * xjac * (1.d0 + zeta)                 &  
                      + w0 * (v_s * u_t - v_t * u_s)                   * theta * tstep &
		      + 0.25d0 * (w0_x * u_y - w0_y * u_x)   * ( v_x * u0_y - v_y * u0_x) * xjac * theta * tstep * tstep &
		      + 0.25d0 * (w0_x * u0_y - w0_y * u0_x) * ( v_x * u_y  - v_y * u_x)  * xjac * theta * tstep * tstep 
     
            amat_12 = + w * (v_s * u0_t - v_t * u0_s)                  * theta * tstep &
                      + visco * (v_x * w_x + v_y * w_y)         * xjac * theta * tstep &	      
                      + visco_num * (v_xx + v_yy)*(w_xx + w_yy) * xjac * theta * tstep &
		      + 0.25d0 * (w_x * u0_y - w_y * u0_x) * ( v_x * u0_y - v_y * u0_x) * xjac * theta * tstep * tstep 
		      
!---------------------------------------------------------------- equation 2
            amat_22 =  v * w * xjac                            
            amat_21 = (v_x * u_x + v_y * u_y) * xjac            

            kl1 = index_kl
            kl2 = index_kl + 1

            ELM(ij1,kl1) =  ELM(ij1,kl1) + wst * amat_11
            ELM(ij1,kl2) =  ELM(ij1,kl2) + wst * amat_12

            ELM(ij2,kl1) =  ELM(ij2,kl1) + wst * amat_21
            ELM(ij2,kl2) =  ELM(ij2,kl2) + wst * amat_22

          enddo
        enddo

      enddo
    enddo
  enddo
enddo

return
end subroutine element_matrix
end module mod_elt_matrix
