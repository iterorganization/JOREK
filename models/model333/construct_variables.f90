!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------- Compute the RZ-coordinates and the Jacobians ------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
subroutine ELM_build_RZ_and_Jacobians(element, nodes, ms, mt)
!DEC$ ATTRIBUTES FORCEINLINE :: ELM_build_RZ_and_Jacobians

  ! --- Modules
  use mod_parameters    
  use basis_at_gaussian
  use equation_variables
  use data_structure
  
  implicit none
  
  ! --- Routine Variables
  type (type_element)	      :: element
  type (type_node)	      :: nodes(n_vertex_max)
  integer		      :: ms, mt
  
  ! --- Internal Variables
  integer		      :: i, j
      
  ! --- Empty before integration
  x_g = 0.d0 ; x_s = 0.d0 ; x_t = 0.d0; x_st  = 0.d0; x_ss  = 0.d0; x_tt  = 0.d0
  y_g = 0.d0 ; y_s = 0.d0 ; y_t = 0.d0; y_st  = 0.d0; y_ss  = 0.d0; y_tt  = 0.d0

  ! --- Integrate
  do i=1,n_vertex_max
    do j=1,n_degrees

      x_g  = x_g  + nodes(i)%x(1,j,1) * element%size(i,j) * H   (i,j,ms,mt)
      x_s  = x_s  + nodes(i)%x(1,j,1) * element%size(i,j) * H_s (i,j,ms,mt)
      x_t  = x_t  + nodes(i)%x(1,j,1) * element%size(i,j) * H_t (i,j,ms,mt)

      x_ss = x_ss + nodes(i)%x(1,j,1) * element%size(i,j) * H_ss(i,j,ms,mt)
      x_st = x_st + nodes(i)%x(1,j,1) * element%size(i,j) * H_st(i,j,ms,mt)
      x_tt = x_tt + nodes(i)%x(1,j,1) * element%size(i,j) * H_tt(i,j,ms,mt)

      y_g  = y_g  + nodes(i)%x(1,j,2) * element%size(i,j) * H   (i,j,ms,mt)
      y_s  = y_s  + nodes(i)%x(1,j,2) * element%size(i,j) * H_s (i,j,ms,mt)
      y_t  = y_t  + nodes(i)%x(1,j,2) * element%size(i,j) * H_t (i,j,ms,mt)

      y_ss = y_ss + nodes(i)%x(1,j,2) * element%size(i,j) * H_ss(i,j,ms,mt)
      y_st = y_st + nodes(i)%x(1,j,2) * element%size(i,j) * H_st(i,j,ms,mt)
      y_tt = y_tt + nodes(i)%x(1,j,2) * element%size(i,j) * H_tt(i,j,ms,mt)
    
    enddo
  enddo
  
  R    = x_g
  R_x  = 1.d0
  
  ! --- Jacobians
  xjac    = x_s*y_t - x_t*y_s

  xjac_x  = (x_ss* y_t**2 - y_ss*x_t*y_t - 2.d0*x_st*y_s*y_t   &       
	   + y_st*(x_s*y_t + x_t*y_s)			       &
	   + x_tt* y_s**2 - y_tt*x_s*y_s) / xjac
	 
  xjac_y  = (y_tt* x_s**2 - x_tt*y_s*x_s - 2.d0*y_st*x_t*x_s   &       
	   + x_st*(y_t*x_s + y_s*x_t)			       &
	   + y_ss* x_t**2 - x_ss*y_t*x_t) / xjac

  return

end subroutine ELM_build_RZ_and_Jacobians








!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------- Compute the variables for the equations --------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
subroutine ELM_build_variables(element, nodes, ms, mt, i_plane)
!DEC$ ATTRIBUTES FORCEINLINE :: ELM_build_variables

  ! --- Modules
  use mod_parameters    
  use basis_at_gaussian
  use equation_variables
  use data_structure
  use phys_module
  use corr_neg
  
  implicit none
  
  ! --- Routine variables
  type (type_element)	      :: element
  type (type_node)	      :: nodes(n_vertex_max)
  integer		      :: ms, mt, i_plane
  
  ! --- Internal variables
  integer		      :: i, j, k, i_tor
      
  ! --- Empty before integration
  ps0	= 0.d0; ps0_s	= 0.d0; ps0_t	= 0.d0; ps0_ss	 = 0.d0; ps0_tt   = 0.d0; ps0_st   = 0.d0; ps0_p   = 0.d0; ps0_pp   = 0.d0
  u0	= 0.d0; u0_s	= 0.d0; u0_t	= 0.d0; u0_ss	 = 0.d0; u0_tt    = 0.d0; u0_st    = 0.d0; u0_p    = 0.d0; u0_pp    = 0.d0
  zj0	= 0.d0; zj0_s	= 0.d0; zj0_t	= 0.d0; zj0_ss	 = 0.d0; zj0_tt   = 0.d0; zj0_st   = 0.d0; zj0_p   = 0.d0; zj0_pp   = 0.d0
  w0	= 0.d0; w0_s	= 0.d0; w0_t	= 0.d0; w0_ss	 = 0.d0; w0_tt    = 0.d0; w0_st    = 0.d0; w0_p    = 0.d0; w0_pp    = 0.d0! --- n=0 variables
  r0	= 0.d0; r0_s	= 0.d0; r0_t	= 0.d0; r0_ss	 = 0.d0; r0_tt    = 0.d0; r0_st    = 0.d0; r0_p    = 0.d0; r0_pp    = 0.d0; r00 = 0.d0
  T0	= 0.d0; T0_s	= 0.d0; T0_t	= 0.d0; T0_ss    = 0.d0; T0_tt    = 0.d0; T0_st    = 0.d0; T0_p    = 0.d0; T0_pp    = 0.d0; T00 = 0.d0
  Vpar0 = 0.d0; Vpar0_s = 0.d0; Vpar0_t = 0.d0; Vpar0_ss = 0.d0; Vpar0_tt = 0.d0; Vpar0_st = 0.d0; Vpar0_p = 0.d0; Vpar0_pp = 0.d0
  delta_g = 0.d0 ; delta_s = 0.d0 ; delta_t = 0.d0

  ! --- Integrate
  do i =1,n_vertex_max
    do j=1,n_degrees
      
      ! --- Axisymmetric variables for localised sources
      r00            = r00        + nodes(i)%values(1    ,j,5) * element%size(i,j) * H   (i,j,ms,mt) * HZ   (1    ,i_plane)
      T00            = T00        + nodes(i)%values(1    ,j,6) * element%size(i,j) * H   (i,j,ms,mt) * HZ   (1    ,i_plane)
      
      do i_tor =1,n_tor

	! --- Variable 1
	ps0	     = ps0	  + nodes(i)%values(i_tor,j,1) * element%size(i,j) * H   (i,j,ms,mt) * HZ   (i_tor,i_plane)
	ps0_s	     = ps0_s	  + nodes(i)%values(i_tor,j,1) * element%size(i,j) * H_s (i,j,ms,mt) * HZ   (i_tor,i_plane)
	ps0_t	     = ps0_t	  + nodes(i)%values(i_tor,j,1) * element%size(i,j) * H_t (i,j,ms,mt) * HZ   (i_tor,i_plane)
	ps0_ss	     = ps0_ss	  + nodes(i)%values(i_tor,j,1) * element%size(i,j) * H_ss(i,j,ms,mt) * HZ   (i_tor,i_plane)
	ps0_tt	     = ps0_tt	  + nodes(i)%values(i_tor,j,1) * element%size(i,j) * H_tt(i,j,ms,mt) * HZ   (i_tor,i_plane)
	ps0_st	     = ps0_st	  + nodes(i)%values(i_tor,j,1) * element%size(i,j) * H_st(i,j,ms,mt) * HZ   (i_tor,i_plane)
	ps0_p	     = ps0_p	  + nodes(i)%values(i_tor,j,1) * element%size(i,j) * H   (i,j,ms,mt) * HZ_p (i_tor,i_plane)
	ps0_pp	     = ps0_pp	  + nodes(i)%values(i_tor,j,1) * element%size(i,j) * H   (i,j,ms,mt) * HZ_pp(i_tor,i_plane)

	! --- Variable 2
	u0	     = u0	  + nodes(i)%values(i_tor,j,2) * element%size(i,j) * H   (i,j,ms,mt) * HZ   (i_tor,i_plane)
	u0_s	     = u0_s	  + nodes(i)%values(i_tor,j,2) * element%size(i,j) * H_s (i,j,ms,mt) * HZ   (i_tor,i_plane)
	u0_t	     = u0_t	  + nodes(i)%values(i_tor,j,2) * element%size(i,j) * H_t (i,j,ms,mt) * HZ   (i_tor,i_plane)
	u0_ss	     = u0_ss	  + nodes(i)%values(i_tor,j,2) * element%size(i,j) * H_ss(i,j,ms,mt) * HZ   (i_tor,i_plane)
	u0_tt	     = u0_tt	  + nodes(i)%values(i_tor,j,2) * element%size(i,j) * H_tt(i,j,ms,mt) * HZ   (i_tor,i_plane)
	u0_st	     = u0_st	  + nodes(i)%values(i_tor,j,2) * element%size(i,j) * H_st(i,j,ms,mt) * HZ   (i_tor,i_plane)
	u0_p	     = u0_p	  + nodes(i)%values(i_tor,j,2) * element%size(i,j) * H   (i,j,ms,mt) * HZ_p (i_tor,i_plane)
	u0_pp	     = u0_pp	  + nodes(i)%values(i_tor,j,2) * element%size(i,j) * H   (i,j,ms,mt) * HZ_pp(i_tor,i_plane)

	! --- Variable 3
	zj0	     = zj0	  + nodes(i)%values(i_tor,j,3) * element%size(i,j) * H   (i,j,ms,mt) * HZ   (i_tor,i_plane)
	zj0_s	     = zj0_s	  + nodes(i)%values(i_tor,j,3) * element%size(i,j) * H_s (i,j,ms,mt) * HZ   (i_tor,i_plane)
	zj0_t	     = zj0_t	  + nodes(i)%values(i_tor,j,3) * element%size(i,j) * H_t (i,j,ms,mt) * HZ   (i_tor,i_plane)
	zj0_ss	     = zj0_ss	  + nodes(i)%values(i_tor,j,3) * element%size(i,j) * H_ss(i,j,ms,mt) * HZ   (i_tor,i_plane)
	zj0_tt	     = zj0_tt	  + nodes(i)%values(i_tor,j,3) * element%size(i,j) * H_tt(i,j,ms,mt) * HZ   (i_tor,i_plane)
	zj0_st	     = zj0_st	  + nodes(i)%values(i_tor,j,3) * element%size(i,j) * H_st(i,j,ms,mt) * HZ   (i_tor,i_plane)
	zj0_p	     = zj0_p	  + nodes(i)%values(i_tor,j,3) * element%size(i,j) * H   (i,j,ms,mt) * HZ_p (i_tor,i_plane)
	zj0_pp	     = zj0_pp	  + nodes(i)%values(i_tor,j,3) * element%size(i,j) * H   (i,j,ms,mt) * HZ_pp(i_tor,i_plane)

	! --- Variable 4
	w0	     = w0	  + nodes(i)%values(i_tor,j,4) * element%size(i,j) * H   (i,j,ms,mt) * HZ   (i_tor,i_plane)
	w0_s	     = w0_s	  + nodes(i)%values(i_tor,j,4) * element%size(i,j) * H_s (i,j,ms,mt) * HZ   (i_tor,i_plane)
	w0_t	     = w0_t	  + nodes(i)%values(i_tor,j,4) * element%size(i,j) * H_t (i,j,ms,mt) * HZ   (i_tor,i_plane)
	w0_ss	     = w0_ss	  + nodes(i)%values(i_tor,j,4) * element%size(i,j) * H_ss(i,j,ms,mt) * HZ   (i_tor,i_plane)
	w0_tt	     = w0_tt	  + nodes(i)%values(i_tor,j,4) * element%size(i,j) * H_tt(i,j,ms,mt) * HZ   (i_tor,i_plane)
	w0_st	     = w0_st	  + nodes(i)%values(i_tor,j,4) * element%size(i,j) * H_st(i,j,ms,mt) * HZ   (i_tor,i_plane)
	w0_p	     = w0_p	  + nodes(i)%values(i_tor,j,4) * element%size(i,j) * H   (i,j,ms,mt) * HZ_p (i_tor,i_plane)
	w0_pp	     = w0_pp	  + nodes(i)%values(i_tor,j,4) * element%size(i,j) * H   (i,j,ms,mt) * HZ_pp(i_tor,i_plane)

	! --- Variable 5
	r0	     = r0	  + nodes(i)%values(i_tor,j,5) * element%size(i,j) * H   (i,j,ms,mt) * HZ   (i_tor,i_plane)
	r0_s	     = r0_s	  + nodes(i)%values(i_tor,j,5) * element%size(i,j) * H_s (i,j,ms,mt) * HZ   (i_tor,i_plane)
	r0_t	     = r0_t	  + nodes(i)%values(i_tor,j,5) * element%size(i,j) * H_t (i,j,ms,mt) * HZ   (i_tor,i_plane)
	r0_ss	     = r0_ss	  + nodes(i)%values(i_tor,j,5) * element%size(i,j) * H_ss(i,j,ms,mt) * HZ   (i_tor,i_plane)
	r0_tt	     = r0_tt	  + nodes(i)%values(i_tor,j,5) * element%size(i,j) * H_tt(i,j,ms,mt) * HZ   (i_tor,i_plane)
	r0_st	     = r0_st	  + nodes(i)%values(i_tor,j,5) * element%size(i,j) * H_st(i,j,ms,mt) * HZ   (i_tor,i_plane)
	r0_p	     = r0_p	  + nodes(i)%values(i_tor,j,5) * element%size(i,j) * H   (i,j,ms,mt) * HZ_p (i_tor,i_plane)
	r0_pp	     = r0_pp	  + nodes(i)%values(i_tor,j,5) * element%size(i,j) * H   (i,j,ms,mt) * HZ_pp(i_tor,i_plane)

	! --- Variable 6
	T0	     = T0	  + nodes(i)%values(i_tor,j,6) * element%size(i,j) * H   (i,j,ms,mt) * HZ   (i_tor,i_plane)
	T0_s	     = T0_s	  + nodes(i)%values(i_tor,j,6) * element%size(i,j) * H_s (i,j,ms,mt) * HZ   (i_tor,i_plane)
	T0_t	     = T0_t	  + nodes(i)%values(i_tor,j,6) * element%size(i,j) * H_t (i,j,ms,mt) * HZ   (i_tor,i_plane)
	T0_ss	     = T0_ss	  + nodes(i)%values(i_tor,j,6) * element%size(i,j) * H_ss(i,j,ms,mt) * HZ   (i_tor,i_plane)
	T0_tt	     = T0_tt	  + nodes(i)%values(i_tor,j,6) * element%size(i,j) * H_tt(i,j,ms,mt) * HZ   (i_tor,i_plane)
	T0_st	     = T0_st	  + nodes(i)%values(i_tor,j,6) * element%size(i,j) * H_st(i,j,ms,mt) * HZ   (i_tor,i_plane)
	T0_p	     = T0_p	  + nodes(i)%values(i_tor,j,6) * element%size(i,j) * H   (i,j,ms,mt) * HZ_p (i_tor,i_plane)
	T0_pp	     = T0_pp	  + nodes(i)%values(i_tor,j,6) * element%size(i,j) * H   (i,j,ms,mt) * HZ_pp(i_tor,i_plane)
  
	! --- Variable 7
	Vpar0	     = Vpar0	  + nodes(i)%values(i_tor,j,7) * element%size(i,j) * H   (i,j,ms,mt) * HZ   (i_tor,i_plane)
	Vpar0_s      = Vpar0_s    + nodes(i)%values(i_tor,j,7) * element%size(i,j) * H_s (i,j,ms,mt) * HZ   (i_tor,i_plane)
	Vpar0_t      = Vpar0_t    + nodes(i)%values(i_tor,j,7) * element%size(i,j) * H_t (i,j,ms,mt) * HZ   (i_tor,i_plane)
	Vpar0_ss     = Vpar0_ss   + nodes(i)%values(i_tor,j,7) * element%size(i,j) * H_ss(i,j,ms,mt) * HZ   (i_tor,i_plane)
	Vpar0_tt     = Vpar0_tt   + nodes(i)%values(i_tor,j,7) * element%size(i,j) * H_tt(i,j,ms,mt) * HZ   (i_tor,i_plane)
	Vpar0_st     = Vpar0_st   + nodes(i)%values(i_tor,j,7) * element%size(i,j) * H_st(i,j,ms,mt) * HZ   (i_tor,i_plane)
	Vpar0_p      = Vpar0_p    + nodes(i)%values(i_tor,j,7) * element%size(i,j) * H   (i,j,ms,mt) * HZ_p (i_tor,i_plane)
	Vpar0_pp     = Vpar0_pp	  + nodes(i)%values(i_tor,j,7) * element%size(i,j) * H   (i,j,ms,mt) * HZ_pp(i_tor,i_plane)

	! --- Deltas
	do k=1,n_var
	  delta_g(k) = delta_g(k) + nodes(i)%deltas(i_tor,j,k) * element%size(i,j) * H   (i,j,ms,mt) * HZ   (i_tor,i_plane)
	  delta_s(k) = delta_s(k) + nodes(i)%deltas(i_tor,j,k) * element%size(i,j) * H_s (i,j,ms,mt) * HZ   (i_tor,i_plane)
	  delta_t(k) = delta_t(k) + nodes(i)%deltas(i_tor,j,k) * element%size(i,j) * H_t (i,j,ms,mt) * HZ   (i_tor,i_plane)
	enddo			
  
      enddo
    enddo
  enddo

  ! changes deltas for variable time steps
  delta_g = delta_g * tstep / tstep_prev
  delta_s = delta_s * tstep / tstep_prev
  delta_t = delta_t * tstep / tstep_prev

  ! --- Variable 1
  ps0_x    = get_deriv_x (ps0_s, ps0_t)
  ps0_y    = get_deriv_y (ps0_s, ps0_t)
  ps0_xx   = get_deriv_xx(ps0_s, ps0_t, ps0_ss, ps0_st, ps0_tt)
  ps0_yy   = get_deriv_yy(ps0_s, ps0_t, ps0_ss, ps0_st, ps0_tt)
  ps0_xy   = get_deriv_xy(ps0_s, ps0_t, ps0_ss, ps0_st, ps0_tt)

  ! --- Variable 2
  u0_x	   = get_deriv_x (u0_s, u0_t)
  u0_y	   = get_deriv_y (u0_s, u0_t)
  u0_xx    = get_deriv_xx(u0_s, u0_t, u0_ss, u0_st, u0_tt)
  u0_yy    = get_deriv_yy(u0_s, u0_t, u0_ss, u0_st, u0_tt)
  u0_xy    = get_deriv_xy(u0_s, u0_t, u0_ss, u0_st, u0_tt)
  vv2	   = R**2 *  ( u0_x * u0_x + u0_y *u0_y  )
  
  ! --- Variable 3
  zj0_x	   = get_deriv_x (zj0_s, zj0_t)
  zj0_y	   = get_deriv_y (zj0_s, zj0_t)
  zj0_xx   = get_deriv_xx(zj0_s, zj0_t, zj0_ss, zj0_st, zj0_tt)
  zj0_yy   = get_deriv_yy(zj0_s, zj0_t, zj0_ss, zj0_st, zj0_tt)
  zj0_xy   = get_deriv_xy(zj0_s, zj0_t, zj0_ss, zj0_st, zj0_tt)
  
  ! --- Variable 4
  w0_x	   = get_deriv_x (w0_s, w0_t)
  w0_y	   = get_deriv_y (w0_s, w0_t)
  w0_xx    = get_deriv_xx(w0_s, w0_t, w0_ss, w0_st, w0_tt)
  w0_yy    = get_deriv_yy(w0_s, w0_t, w0_ss, w0_st, w0_tt)
  w0_xy    = get_deriv_xy(w0_s, w0_t, w0_ss, w0_st, w0_tt)
  
  ! --- Variable 5
  r0_corr  = corr_neg_dens(r0)
  r0_corr2 = corr_neg_dens(r0, (/0.5,0.5/) ) ! A second one specially for the diamagnetic terms
  r0_x	   = get_deriv_x (r0_s, r0_t)
  r0_y	   = get_deriv_y (r0_s, r0_t)
  r0_xx    = get_deriv_xx(r0_s, r0_t, r0_ss, r0_st, r0_tt)
  r0_yy    = get_deriv_yy(r0_s, r0_t, r0_ss, r0_st, r0_tt)
  r0_xy    = get_deriv_xy(r0_s, r0_t, r0_ss, r0_st, r0_tt)
  r0_hat   = R**2 * r0
  r0_x_hat = 2.d0 * R * R_x  * r0 + R**2 * r0_x
  r0_y_hat = R**2 * r0_y
  
  ! --- Variable 6
  T0_corr  = corr_neg_temp(T0) ! For use in eta(T), visco(T), ...
  T0_x	   = get_deriv_x (T0_s, T0_t)
  T0_y	   = get_deriv_y (T0_s, T0_t)
  T0_xx    = get_deriv_xx(T0_s, T0_t, T0_ss, T0_st, T0_tt)
  T0_yy    = get_deriv_yy(T0_s, T0_t, T0_ss, T0_st, T0_tt)
  T0_xy    = get_deriv_xy(T0_s, T0_t, T0_ss, T0_st, T0_tt)

  ! --- Variable 7
  Vpar0_x  = get_deriv_x (Vpar0_s, Vpar0_t)
  Vpar0_y  = get_deriv_y (Vpar0_s, Vpar0_t)
  Vpar0_xx = get_deriv_xx(Vpar0_s, Vpar0_t, Vpar0_ss, Vpar0_st, Vpar0_tt)
  Vpar0_yy = get_deriv_yy(Vpar0_s, Vpar0_t, Vpar0_ss, Vpar0_st, Vpar0_tt)
  Vpar0_xy = get_deriv_xy(Vpar0_s, Vpar0_t, Vpar0_ss, Vpar0_st, Vpar0_tt)

  ! --- Deltas
  delta_u_x  = (   y_t * delta_s(2) - y_s * delta_t(2) ) / xjac
  delta_u_y  = ( - x_t * delta_s(2) + x_s * delta_t(2) ) / xjac
  delta_ps_x = (   y_t * delta_s(1) - y_s * delta_t(1) ) / xjac
  delta_ps_y = ( - x_t * delta_s(1) + x_s * delta_t(1) ) / xjac
  
  ! --- Pressure
  P0	   = r0    * T0
  P0_x     = r0_x  * T0 + r0 * T0_x
  P0_y     = r0_y  * T0 + r0 * T0_y
  P0_s     = r0_s  * T0 + r0 * T0_s
  P0_t     = r0_t  * T0 + r0 * T0_t
  P0_p     = r0_p  * T0 + r0 * T0_p
  P0_pp    = r0_pp * T0 + r0 * T0_pp + 2.d0 * r0_p * T0_p
  P0_xx    = r0_xx * T0 + r0 * T0_xx + 2.d0 * r0_x * T0_x
  P0_yy    = r0_yy * T0 + r0 * T0_yy + 2.d0 * r0_y * T0_y
  P0_xy    = r0_xy * T0 + r0 * T0_xy + r0_x * T0_y + r0_y * T0_x
  
  ! --- Magnetic field amplitude (squared)
  BB2	   = (F0*F0 + ps0_x * ps0_x + ps0_y * ps0_y )/R**2
  
  
  return

contains
! --- Function to compute derivative with respect to x
real*8 function get_deriv_x(VAR_s, VAR_t)
  use equation_variables
  real*8, intent(in) :: VAR_s, VAR_t
  get_deriv_x = (   y_t * VAR_s - y_s * VAR_t ) / xjac
end function get_deriv_x

! --- Function to compute derivative with respect to y
real*8 function get_deriv_y(VAR_s, VAR_t)
  use equation_variables
  real*8, intent(in) :: VAR_s, VAR_t
  get_deriv_y = ( - x_t * VAR_s + x_s * VAR_t ) / xjac
end function get_deriv_y

! --- Function to compute derivative with respect to xx
real*8 function get_deriv_xx(VAR_s, VAR_t, VAR_ss, VAR_st, VAR_tt)
  use equation_variables
  real*8, intent(in) :: VAR_s, VAR_t, VAR_ss, VAR_st, VAR_tt
  get_deriv_xx = (VAR_ss * y_t**2 - 2.d0*VAR_st * y_s*y_t + VAR_tt * y_s**2  & 	    
	        + VAR_s  * (y_st*y_t - y_tt*y_s )			     &    
	        + VAR_t  * (y_st*y_s - y_ss*y_t ) )        / xjac**2         & 	
	        - xjac_x * (VAR_s * y_t - VAR_t * y_s)     / xjac**2
end function get_deriv_xx

! --- Function to compute derivative with respect to yy
real*8 function get_deriv_yy(VAR_s, VAR_t, VAR_ss, VAR_st, VAR_tt)
  use equation_variables
  real*8, intent(in) :: VAR_s, VAR_t, VAR_ss, VAR_st, VAR_tt
  get_deriv_yy = (VAR_ss * x_t**2 - 2.d0*VAR_st * x_s*x_t + VAR_tt * x_s**2  & 	    
	        + VAR_s * (x_st*x_t - x_tt*x_s )			     &    
	        + VAR_t * (x_st*x_s - x_ss*x_t ) )         / xjac**2         & 	
	        - xjac_y * (- VAR_s * x_t + VAR_t * x_s )  / xjac**2
end function get_deriv_yy

! --- Function to compute derivative with respect to yy
real*8 function get_deriv_xy(VAR_s, VAR_t, VAR_ss, VAR_st, VAR_tt)
  use equation_variables
  real*8, intent(in) :: VAR_s, VAR_t, VAR_ss, VAR_st, VAR_tt
  get_deriv_xy = (- VAR_ss * y_t*x_t - VAR_tt * x_s*y_s		         &
	          + VAR_st * (y_s*x_t  + y_t*x_s  )  		         &        
	          - VAR_s  * (x_st*y_t - x_tt*y_s )  		         &    
	          - VAR_t  * (x_st*y_s - x_ss*y_t )  )       / xjac**2   & 	
	          - xjac_x * (- VAR_s * x_t + VAR_t * x_s )  / xjac**2
end function get_deriv_xy
end subroutine ELM_build_variables







!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------- Compute the diffusivities and source ---------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
subroutine ELM_build_diffusivities_and_sources(element, nodes, xpoint2, xcase2, R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint, Z_xpoint, i_plane)
!DEC$ ATTRIBUTES FORCEINLINE :: ELM_build_diffusivities_and_sources

  ! --- Modules
  use mod_parameters    
  use basis_at_gaussian
  use phys_module
  use equation_variables
  use data_structure
  use diffusivities, only: get_dperp, get_zkperp
  use pellet_module
  use mod_bootstrap_functions
  use equil_info, only : get_psi_n
  use mod_sources
  
  implicit none
  
  ! --- Routine variables
  type (type_element)	      :: element
  type (type_node)	      :: nodes(n_vertex_max)
  logical		      :: xpoint2
  integer		      :: xcase2
  integer		      :: i_plane
  real*8		      :: R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint(2), Z_xpoint(2)
  
  ! --- Internal variables
  real*8		      :: psi_norm
  real*8		      :: V_source, dV_dpsi2, dV_dz2, dV_dpsi_dz, dV_dpsi3,dV_dpsi_dz2, dV_dpsi2_dz
  real*8		      :: Ti0, Ti0_x, Ti0_y, Te0, Te0_x, Te0_y
  real*8		      :: zTi, zTi_x, zTi_y, zTe, zTe_x, zTe_y, zn_x, zn_y
  real*8		      :: Jb_0
  real*8		      :: rho_norm
  real*8		      :: above_prof
      
  
  ! -------------------------------------
  ! --- Temperature dependent resistivity
  ! -------------------------------------
  if ( eta_T_dependent .and. T0_corr <= T_max_eta) then
    eta_T     =   eta   * (T0_corr / T_0)**(-1.5d0)
    deta_dT   = - eta	* (1.5d0)  * T0_corr**(-2.5d0) * T_0**(1.5d0)
    d2eta_d2T =   eta	* (3.75d0) * T0_corr**(-3.5d0) * T_0**(1.5d0)
  else if ( eta_T_dependent .and. T0_corr > T_max_eta) then
     eta_T     = eta   * (T_max_eta/T_0)**(-1.5d0)
     deta_dT   = 0.
     d2eta_d2T = 0.     
  else
    eta_T     = eta
    deta_dT   = 0.d0
    d2eta_d2T = 0.d0
  end if
  

  ! -------------------------
  ! --- Eta for ohmic heating
  ! -------------------------
  if ( eta_T_dependent .and. T0_corr <= T_max_eta_ohm) then
    eta_T_ohm     = eta_ohmic   * (T0_corr/T_0)**(-1.5d0)
    deta_dT_ohm   = - eta_ohmic   * (1.5d0)  * T0_corr**(-2.5d0) * T_0**(1.5d0)
  else if ( eta_T_dependent .and. T0_corr > T_max_eta_ohm) then
    eta_T_ohm     = eta_ohmic   * (T_max_eta_ohm/T_0)**(-1.5d0)
    deta_dT_ohm   = 0.    
  else
    eta_T_ohm     = eta_ohmic
    deta_dT_ohm   = 0.d0
  end if  
   

  ! -----------------------------------
  ! --- Temperature dependent viscosity
  ! -----------------------------------
  if ( visco_T_dependent ) then       
    visco_T     =   visco * (T0_corr/T_0)**(-1.5d0)
    dvisco_dT   = - visco * (1.5d0)  * T0_corr**(-2.5d0) * T_0**(1.5d0)
    d2visco_dT2 =   visco * (3.75d0) * T0_corr**(-3.5d0) * T_0**(1.5d0)
  else
    visco_T     = visco
    dvisco_dT   = 0.d0
    d2visco_dT2 = 0.d0
  end if
  visco_parr  = visco_par  
  
  
  ! -------------------------------------------------------------
  ! --- D_perp and K_perp profiles (for fixed pedestal gradients)
  ! -------------------------------------------------------------
  ! --- First need psi_norm
  psi_norm = get_psi_n(ps0, y_g)

  ! --- Call Diff functions (same as before, but with additional profile if D_perp(10)=1.d0 or ZK_perp(10) = 1.d0)
  D_prof = get_dperp (ps0, psi_norm, psi_axis, psi_bnd, y_g, Z_xpoint)
  K_prof = get_zkperp(ps0, psi_norm, psi_axis, psi_bnd, y_g, Z_xpoint)
  
  ! --- Increase diffusivity if very small density/temperature
  if (xpoint2) then
    if (r0 .lt. D_prof_neg_thresh)  then
      D_prof = D_prof_neg
    endif
    if (T0 .lt. ZK_prof_neg_thresh) then
      K_prof = ZK_prof_neg
    endif
  endif
  
  ! -----------------------------------------------------
  ! --- Parallel conductivity profiles (Braginskii model)
  ! -----------------------------------------------------
  if ( ZKpar_T_dependent ) then
    K_par    = ZK_par * (T0_corr/T_0)**(+2.5d0)
    dK_par   = ZK_par * (2.5d0)  * T0_corr**(+1.5d0) * T_0**(-2.5d0)
    if (K_par .gt. ZK_par_max) then
      K_par  = Zk_par_max
      dK_par = 0.d0
    endif
  else
    K_par  = ZK_par
    dK_par = 0.d0
  endif
  
 
  ! -------------------------
  ! --- Hyper diffusivitities
  ! -------------------------
  eta_numm	 = eta_num		! hyper-resistivity
  visco_numm	 = visco_num		! hyper-viscosity
  visco_par_numm = visco_par_num	! hyper-viscosity
  D_perp_numm	 = D_perp_num		! hyper-diffusivity
  K_perp_numm	 = ZK_perp_num		! hyper-conductivity

  
  ! ------------------------------------------------
  ! --- Taylor Galerkin (TG2) stabilisation switches
  ! ------------------------------------------------
  TG_num1 = tgnum(1);
  TG_num2 = tgnum(2);
  TG_num5 = tgnum(5);
  TG_num6 = tgnum(6);
  TG_num7 = tgnum(7);
  
  
  ! ---------------------
  ! --- Bootstrap current
  ! ---------------------
  if (bootstrap) then
    ! --- Full Sauter formula
    Ti0   = T0   / 2.d0 ; Te0	= T0   / 2.d0
    Ti0_x = T0_x / 2.d0 ; Te0_x = T0_x / 2.d0
    Ti0_y = T0_y / 2.d0 ; Te0_y = T0_y / 2.d0
    call bootstrap_current(R, y_g,                               &
                           R_axis,   Z_axis,   psi_axis,         &
			   R_xpoint, Z_xpoint, psi_bnd, psi_norm,&
			   ps0, ps0_x, ps0_y,                    &
			   r0,  r0_x,  r0_y,                     &
			   Ti0, Ti0_x, Ti0_y,                    &
			   Te0, Te0_x, Te0_y,                  Jb)
    ! --- Full Sauter formula for initial profiles
    call density(    xpoint2, xcase2, y_g, Z_xpoint, ps0,psi_axis,psi_bnd, &
    		     zn,dn_dpsi,  dn_dz, dn_dpsi2, dn_dz2, dn_dpsi_dz, dn_dpsi3, dn_dpsi_dz2, dn_dpsi2_dz)
    call temperature(xpoint2, xcase2, y_g, Z_xpoint, ps0,psi_axis,psi_bnd, &
    		     zT,dT_dpsi,dT_dz,dT_dpsi2,dT_dz2,dT_dpsi_dz,dT_dpsi3,dT_dpsi_dz2,dT_dpsi2_dz)
    zTi   = zT / 2.d0             
    zTi_x = dT_dpsi * ps0_x / 2.d0
    zTi_y = dT_dpsi * ps0_y / 2.d0
    zTe   = zTi  
    zTe_x = zTi_x
    zTe_y = zTi_y
    zn_x  = dn_dpsi * ps0_x
    zn_y  = dn_dpsi * ps0_y
    call bootstrap_current(R, y_g,                               &
                           R_axis,   Z_axis,   psi_axis,         &
			   R_xpoint, Z_xpoint, psi_bnd, psi_norm,&
			   ps0, ps0_x, ps0_y,                    &
			   zn,  zn_x,  zn_y,                     &
			   zTi, zTi_x, zTi_y,                    &
			   zTe, zTe_x, zTe_y,                  Jb_0)
    ! --- Subtract the initial equilibrium part
    Jb = Jb - Jb_0
  else
    Jb = 0.d0
  endif
  
  
  ! ------------------------------------------------------
  ! --- Diamagnetic terms, avoid problems at the target...
  ! ------------------------------------------------------
  tau_IC = tauIC
  
  ! -------------------------
  ! --- Neoclassical rotation
  ! -------------------------
  epsil   = 1.d-3
  Btheta2 = (ps0_x**2 + ps0_y**2) / R**2
    amu_neo_prof   = 0.d0
    aki_neo_prof   = 0.d0
  if ( NEO ) then 
    if (num_neo_file) then
      call neo_coef(xpoint2, xcase2, y_g, Z_xpoint, ps0, psi_axis, psi_bnd, amu_neo_prof,          &
        aki_neo_prof)
    else
       amu_neo_prof = amu_neo_const
       aki_neo_prof = aki_neo_const
    endif
  endif
  
  ! -------------------------------------------------------------------
  ! --- Heating, current and particle source (the same for all i_plane)
  ! -------------------------------------------------------------------
  if (i_plane .eq. 1) then
    ! --- Current source

    current_source = 0.
    if (keep_current_prof) &
      call current(xpoint2, xcase2, x_g,y_g, Z_xpoint, ps0,psi_axis,psi_bnd,current_source)

    ! --- Density and Temperature source
    call sources(xpoint2, xcase2, y_g, Z_xpoint, ps0,psi_axis,psi_bnd,particle_source,heat_source)
    
    ! --- New source profile: source with exactly the same profile as the initial equilibirum profiles.
    call density(    xpoint2, xcase2, y_g, Z_xpoint, ps0,psi_axis,psi_bnd, &
    		     zn,dn_dpsi,  dn_dz, dn_dpsi2, dn_dz2, dn_dpsi_dz, dn_dpsi3, dn_dpsi_dz2, dn_dpsi2_dz)
    call temperature(xpoint2, xcase2, y_g, Z_xpoint, ps0,psi_axis,psi_bnd, &
    		     zT,dT_dpsi,dT_dz,dT_dpsi2,dT_dz2,dT_dpsi_dz,dT_dpsi3,dT_dpsi_dz2,dT_dpsi2_dz)
    
    ! --- Toroidal momentum source (NBI)
    dV_dpsi_source = 0.d0
    dV_dz_source   = 0.d0
    if ( ( abs(V_0) .ge. 1.e-12 ) .or. ( num_rot ) ) then
      call velocity(xpoint2, xcase2, y_g, z_xpoint, ps0, psi_axis, psi_bnd, V_source,               &
        dV_dpsi_source, dV_dz_source, dV_dpsi2, dV_dz2, dV_dpsi_dz, dV_dpsi3,dV_dpsi_dz2,           &
        dV_dpsi2_dz)
    end if
    if (normalized_velocity_profile) then
      Vt0_x = dV_dpsi_source * ps0_x
      Vt0_y = dV_dz_source + dV_dpsi_source * ps0_y
    else
      Omega_tor0_x = dV_dpsi_source * ps0_x
      Omega_tor0_y = dV_dz_source + dV_dpsi_source * ps0_y
    end if
    
    ! --- Pellet Source
    source_pellet = 0.d0
    source_volume = 0.d0
    if (use_pellet) then
      call pellet_source2(pellet_amplitude,pellet_R,pellet_Z,pellet_psi,pellet_phi, &
    			  pellet_radius, pellet_delta_psi, pellet_sig, pellet_length, pellet_ellipse, pellet_theta, &
    			  x_g,y_g, ps0, phi, zn, zT/2.0, &
    			  central_density, pellet_particles, pellet_density, total_pellet_volume, &
    			  source_pellet, source_volume)
    endif
    
    ! --- Total density source 
    total_rho_source = particle_source + source_pellet
  endif
  
  ! -------------------------------------------------------------------
  ! --- Renormalise all MHD parameters
  ! -------------------------------------------------------------------
  
  ! --- (ie. if input values are given in physical units) Not done for sources yet...
  if (renormalise) then
    if (central_density .gt. 1.d10) then
      rho_norm    = central_density         * central_mass * mass_proton
    else
      rho_norm    = central_density * 1.d20 * central_mass * mass_proton
    endif
    eta_T       = eta_T       * sqrt(rho_norm / mu_zero )
    deta_dT     = deta_dT     * sqrt(rho_norm / mu_zero )
    eta_T_ohm   = eta_T_ohm   * sqrt(rho_norm / mu_zero )
    deta_dT_ohm = deta_dT_ohm * sqrt(rho_norm / mu_zero )
    visco_T     = visco_T     * sqrt(mu_zero  / rho_norm)
    dvisco_dT   = dvisco_dT   * sqrt(mu_zero  / rho_norm)
    visco_parr  = visco_par   * sqrt(mu_zero  / rho_norm)
    D_prof      = D_prof      * sqrt(mu_zero  * rho_norm)
    K_prof      = K_prof      * sqrt(mu_zero  / rho_norm)
    K_par       = K_par       * sqrt(mu_zero  / rho_norm)
    dK_par      = dK_par      * sqrt(mu_zero  / rho_norm)
    tau_IC      = tau_IC      / sqrt(mu_zero  * rho_norm)
  endif


  ! ------------------------------------------------------
  ! --- Diamagnetic viscosity
  ! ------------------------------------------------------
  if (Wdia) then
    W_dia = + tau_IC /r0_corr2    * (p0_xx + p0_x/R + p0_yy) &
            - tau_IC /r0_corr2**2 * (r0_x*p0_x + r0_y*p0_y)
  else
    W_dia = 0.d0
  endif
  

  ! -------------------------------------------------------------------
  ! --- SOL sources to stabilise diamagnetic terms (not applied by default!)
  ! -------------------------------------------------------------------
  
  ! --- Make sure SOL density/temperature stays levelled
  if (.false.) then
    if (psi_norm .gt. 1.0) then
      call density(    xpoint2, xcase2, y_g, Z_xpoint, ps0,psi_axis,psi_bnd, &
    		       zn,dn_dpsi,  dn_dz, dn_dpsi2, dn_dz2, dn_dpsi_dz, dn_dpsi3, dn_dpsi_dz2, dn_dpsi2_dz)
      call temperature(xpoint2, xcase2, y_g, Z_xpoint, ps0,psi_axis,psi_bnd, &
    		       zT,dT_dpsi,dT_dz,dT_dpsi2,dT_dz2,dT_dpsi_dz,dT_dpsi3,dT_dpsi_dz2,dT_dpsi2_dz)
      ! --- If equilibrium density comes below half of initial SOL value, we fill it up
      if (r00 .lt. 0.5*zn) total_rho_source = 0.5 * (0.5*zn-r00) / tstep
      if (T00 .lt. 0.5*zT) heat_source = 0.5 * (0.5*zT-T00) / tstep
    endif
  endif
  
  ! --- Make sure SOL density stays above zero (nonlinearly)
  ! --- This is a bit more complex, but negative density is the worst enemy of diamagnetic terms
  ! --- If we see that density becomes less than 30% of initial SOL value, and if this
  ! --- is happening due to the nonlinear ballooning mode, then we fill up the loss
  if (.false.) then
    ! --- Usually this happens outside the separatrix
    if (psi_norm .gt. 1.0) then
      call density(    xpoint2, xcase2, y_g, Z_xpoint, ps0,psi_axis,psi_bnd, &
    		       zn,dn_dpsi,  dn_dz, dn_dpsi2, dn_dz2, dn_dpsi_dz, dn_dpsi3, dn_dpsi_dz2, dn_dpsi2_dz)
      if ( (r0 .lt. 0.3*zn) .and. (r00 .gt. 0.3*zn) ) then
    	total_rho_source = 0.5 * (0.3*zn-r0) / tstep
      endif
    endif
    ! --- But it can also happen in the pedestal, where you might want to use different values
    if ( (psi_norm .gt. 0.95) .and. (psi_norm .le. 1.0) ) then
      call density(    xpoint2, xcase2, y_g, Z_xpoint, ps0,psi_axis,psi_bnd, &
    		       zn,dn_dpsi,  dn_dz, dn_dpsi2, dn_dz2, dn_dpsi_dz, dn_dpsi3, dn_dpsi_dz2, dn_dpsi2_dz)
      if ( (r0 .lt. 0.3*zn) .and. (r00 .gt. 0.3*zn) ) then
    	total_rho_source = 0.5 * (0.3*zn-r0) / tstep
      endif
    endif
  endif

  ! -------------------------------------------------------------------
  ! --- Pedestal sources to get ELM cycles (not applied by default!)
  ! -------------------------------------------------------------------
  
  ! --- Just a source localised in the pedestal
  if (.false.) then
    ! --- Density
    call density(    xpoint2, xcase2, y_g, Z_xpoint, ps0,psi_axis,psi_bnd, &
                     zn,dn_dpsi,  dn_dz, dn_dpsi2, dn_dz2, dn_dpsi_dz, dn_dpsi3, dn_dpsi_dz2, dn_dpsi2_dz)
    if (r00 .lt. 1.3*zn) then ! should not exceed 30% above equilibrium initial value
      total_rho_source = particlesource * (0.5d0 - 0.5d0 * tanh( (psi_norm - 1.0d0)/0.003) ) * (0.5d0 - 0.5d0 * tanh(-(psi_norm - 0.60d0)/0.1) ) &
                       + 1.d-4 * (0.5d0 - 0.5d0 * tanh(-(psi_norm - 0.998d0)/0.003) )
    endif
    ! --- Temperature
    call temperature(xpoint2, xcase2, y_g, Z_xpoint, ps0,psi_axis,psi_bnd, &
                     zT,dT_dpsi,dT_dz,dT_dpsi2,dT_dz2,dT_dpsi_dz,dT_dpsi3,dT_dpsi_dz2,dT_dpsi2_dz)
    if (T00 .lt. 1.3*zT) then ! should not exceed 30% above equilibrium initial value
      heat_source      = heatsource     * (0.5d0 - 0.5d0 * tanh( (psi_norm - 1.0d0)/0.003) ) * (0.5d0 - 0.5d0 * tanh(-(psi_norm - 0.60d0)/0.1) ) &
                       + 4.d-7 * (0.5d0 - 0.5d0 * tanh(-(psi_norm - 0.998d0)/0.003) )
    endif
  endif
  
  ! --- A source adapted to the actual profile and the profile you want to reach in the end (here based on the initial profiles)
  if (.false.) then
    total_rho_source = 0.d0
    heat_source      = 0.d0
    ! --- We do this only inside the separatrix
    if ( ( (ps0-psi_axis)/(psi_bnd-psi_axis) .lt. 1.0) .and. (y_g .gt. Z_xpoint(1)) ) then
      ! --- This determines how far above the initial equilibrium profiles you want to go
      ! --- eg. 0.4 means 40% above, 0.0 means exactly the initial profile, not above
      above_prof = 1.d0 + 0.4d0 * (0.5d0 - 0.5d0 * tanh(-(psi_norm - 0.6d0)/0.3) )
      call density(    xpoint2, xcase2, y_g, Z_xpoint, ps0,psi_axis,psi_bnd, &
                       zn,dn_dpsi,  dn_dz, dn_dpsi2, dn_dz2, dn_dpsi_dz, dn_dpsi3, dn_dpsi_dz2, dn_dpsi2_dz)
      call temperature(xpoint2, xcase2, y_g, Z_xpoint, ps0,psi_axis,psi_bnd, &
                       zT,dT_dpsi,dT_dz,dT_dpsi2,dT_dz2,dT_dpsi_dz,dT_dpsi3,dT_dpsi_dz2,dT_dpsi2_dz)
      ! --- Note that here, the value of particlesource and heatsource determines how long
      ! --- it will take for the target profile to be recovered. eg. if 0.1, it will take
      ! --- 10 JOREK times. If 0.001 it will take 1000 JOREK times
      if (r00 .lt. above_prof*zn) total_rho_source = particlesource * (above_prof*zn-r00) / tstep
      if (T00 .lt. above_prof*zT) heat_source      = heatsource     * (above_prof*zT-T00) / tstep
    else
      total_rho_source = 0.d0
      heat_source      = 0.d0
    endif
  endif

  return

end subroutine ELM_build_diffusivities_and_sources










!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
!-------------------------------- Compute the basis functions (ie. the Bezier polynomials) ------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
subroutine ELM_build_basis_functions(element, nodes, ms, mt, i_plane, i_vertex, i_order, i_tor, &
				     vv, vv_s,  vv_t,	     vv_p,  vv_x,  vv_y,	       &
        				 vv_ss, vv_tt, vv_st, vv_pp, vv_xx, vv_yy, vv_xy        )
!DEC$ ATTRIBUTES FORCEINLINE :: ELM_build_basis_functions

  ! --- Modules
  use mod_parameters    
  use basis_at_gaussian
  use equation_variables
  use data_structure
  
  implicit none
  
  ! --- Routine variables
  type (type_element)	      :: element
  type (type_node)	      :: nodes(n_vertex_max)
  integer		      :: ms, mt, i_plane, i_vertex, i_order, i_tor
  real*8		      :: vv, vv_s,  vv_t,	  vv_p,  vv_x,  vv_y
  real*8		      ::     vv_ss, vv_tt, vv_st, vv_pp, vv_xx, vv_yy, vv_xy
  
  ! --- Internal variables
  real*8, dimension(n_tor,n_plane) :: HHZ, HHZ_p, HHZ_pp
  
  ! --- Toroidal functions	      
  if (n_tor .gt. 3) then
    HHZ   (i_tor,i_plane) = 1.d0
    HHZ_p (i_tor,i_plane) = 1.d0
    HHZ_pp(i_tor,i_plane) = 1.d0
  else
    HHZ   (i_tor,i_plane) = HZ   (i_tor,i_plane)
    HHZ_p (i_tor,i_plane) = HZ_p (i_tor,i_plane)
    HHZ_pp(i_tor,i_plane) = HZ_pp(i_tor,i_plane)
  endif
  
  ! --- Test functions  	      
  vv	= H   (i_vertex,i_order,ms,mt) * element%size(i_vertex,i_order) * HHZ	(i_tor,i_plane)
  vv_s  = H_s (i_vertex,i_order,ms,mt) * element%size(i_vertex,i_order) * HHZ	(i_tor,i_plane)
  vv_t  = H_t (i_vertex,i_order,ms,mt) * element%size(i_vertex,i_order) * HHZ	(i_tor,i_plane)
  vv_p  = H   (i_vertex,i_order,ms,mt) * element%size(i_vertex,i_order) * HHZ_p (i_tor,i_plane)
  vv_pp = H   (i_vertex,i_order,ms,mt) * element%size(i_vertex,i_order) * HHZ_pp(i_tor,i_plane)

  vv_ss = H_ss(i_vertex,i_order,ms,mt) * element%size(i_vertex,i_order) * HHZ	(i_tor,i_plane)
  vv_tt = H_tt(i_vertex,i_order,ms,mt) * element%size(i_vertex,i_order) * HHZ	(i_tor,i_plane)
  vv_st = H_st(i_vertex,i_order,ms,mt) * element%size(i_vertex,i_order) * HHZ	(i_tor,i_plane)

  vv_x = (  y_t * vv_s - y_s * vv_t ) / xjac
  vv_y = (- x_t * vv_s + x_s * vv_t ) / xjac

  vv_xx = (vv_ss * y_t**2 - 2.d0*vv_st * y_s*y_t + vv_tt * y_s**2   &	      
	 + vv_s * (y_st*y_t - y_tt*y_s )			    &	 
	 + vv_t * (y_st*y_s - y_ss*y_t ) )    / xjac**2 	    &	  
	 - xjac_x * (vv_s * y_t - vv_t * y_s) / xjac**2

  vv_yy = (vv_ss * x_t**2 - 2.d0*vv_st * x_s*x_t + vv_tt * x_s**2   &	      
	 + vv_s * (x_st*x_t - x_tt*x_s )			    &	 
	 + vv_t * (x_st*x_s - x_ss*x_t ) )	/ xjac**2	    &	  
	 - xjac_y * (- vv_s * x_t + vv_t * x_s ) / xjac**2

  vv_xy = (- vv_ss * y_t*x_t - vv_tt * x_s*y_s  		&
	   + vv_st * (y_s*x_t  + y_t*x_s  )			&	 
	   - vv_s  * (x_st*y_t - x_tt*y_s )			&	   
	   - vv_t  * (x_st*y_s - x_ss*y_t )  )       / xjac**2  &		
	   - xjac_x * (- vv_s * x_t + vv_t * x_s )   / xjac**2  	      
  
  return

end subroutine ELM_build_basis_functions


