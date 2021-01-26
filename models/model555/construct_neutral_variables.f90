!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------- Compute the variables for the equations --------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
subroutine ELM_build_neutral_variables(element, nodes, ms, mt, i_plane)
!DEC$ ATTRIBUTES FORCEINLINE :: ELM_build_neutral_variables

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
  w0	= 0.d0; w0_s	= 0.d0; w0_t	= 0.d0; w0_ss	 = 0.d0; w0_tt    = 0.d0; w0_st    = 0.d0; w0_p    = 0.d0; w0_pp    = 0.d0
  r0	= 0.d0; r0_s	= 0.d0; r0_t	= 0.d0; r0_ss	 = 0.d0; r0_tt    = 0.d0; r0_st    = 0.d0; r0_p    = 0.d0; r0_pp    = 0.d0
  T0	= 0.d0; T0_s	= 0.d0; T0_t	= 0.d0; T0_ss    = 0.d0; T0_tt    = 0.d0; T0_st    = 0.d0; T0_p    = 0.d0; T0_pp    = 0.d0
  Vpar0 = 0.d0; Vpar0_s = 0.d0; Vpar0_t = 0.d0; Vpar0_ss = 0.d0; Vpar0_tt = 0.d0; Vpar0_st = 0.d0; Vpar0_p = 0.d0; Vpar0_pp = 0.d0
  rn0	= 0.d0; rn0_s	= 0.d0; rn0_t	= 0.d0; rn0_ss	 = 0.d0; rn0_tt   = 0.d0; rn0_st   = 0.d0; rn0_p   = 0.d0; rn0_pp   = 0.d0
  delta_g = 0.d0 ; delta_s = 0.d0 ; delta_t = 0.d0

  ! --- Integrate
  do i =1,n_vertex_max
    do j=1,n_order+1
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

	! --- Variable 8
	rn0	     = rn0	  + nodes(i)%values(i_tor,j,8) * element%size(i,j) * H   (i,j,ms,mt) * HZ   (i_tor,i_plane)
	rn0_s	     = rn0_s	  + nodes(i)%values(i_tor,j,8) * element%size(i,j) * H_s (i,j,ms,mt) * HZ   (i_tor,i_plane)
	rn0_t	     = rn0_t	  + nodes(i)%values(i_tor,j,8) * element%size(i,j) * H_t (i,j,ms,mt) * HZ   (i_tor,i_plane)
	rn0_ss	     = rn0_ss	  + nodes(i)%values(i_tor,j,8) * element%size(i,j) * H_ss(i,j,ms,mt) * HZ   (i_tor,i_plane)
	rn0_tt	     = rn0_tt	  + nodes(i)%values(i_tor,j,8) * element%size(i,j) * H_tt(i,j,ms,mt) * HZ   (i_tor,i_plane)
	rn0_st	     = rn0_st	  + nodes(i)%values(i_tor,j,8) * element%size(i,j) * H_st(i,j,ms,mt) * HZ   (i_tor,i_plane)
	rn0_p	     = rn0_p	  + nodes(i)%values(i_tor,j,8) * element%size(i,j) * H   (i,j,ms,mt) * HZ_p (i_tor,i_plane)
	rn0_pp	     = rn0_pp	  + nodes(i)%values(i_tor,j,8) * element%size(i,j) * H   (i,j,ms,mt) * HZ_pp(i_tor,i_plane)
  
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
  ps0_x    = (   y_t * ps0_s - y_s * ps0_t ) / xjac
  ps0_y    = ( - x_t * ps0_s + x_s * ps0_t ) / xjac
  ps0_xx   = (ps0_ss * y_t**2 - 2.d0*ps0_st * y_s*y_t + ps0_tt * y_s**2  & 	    
	      + ps0_s * (y_st*y_t - y_tt*y_s )			         &    
	      + ps0_t * (y_st*y_s - y_ss*y_t ) )       / xjac**2         & 	
	    - xjac_x * (ps0_s * y_t - ps0_t * y_s)     / xjac**2
  ps0_yy   = (ps0_ss * x_t**2 - 2.d0*ps0_st * x_s*x_t + ps0_tt * x_s**2  & 	    
	      + ps0_s * (x_st*x_t - x_tt*x_s )			         &    
	      + ps0_t * (x_st*x_s - x_ss*x_t ) )       / xjac**2         & 	
	    - xjac_y * (- ps0_s * x_t + ps0_t * x_s )  / xjac**2
  ps0_xy   = (- ps0_ss * y_t*x_t - ps0_tt * x_s*y_s		         &
	      + ps0_st * (y_s*x_t  + y_t*x_s  )  		         &        
	      - ps0_s  * (x_st*y_t - x_tt*y_s )  		         &    
	      - ps0_t  * (x_st*y_s - x_ss*y_t )  )     / xjac**2         & 	
	    - xjac_x * (- ps0_s * x_t + ps0_t * x_s )  / xjac**2

  ! --- Variable 2
  u0_x     = (   y_t * u0_s - y_s * u0_t ) / xjac
  u0_y     = ( - x_t * u0_s + x_s * u0_t ) / xjac
  u0_xx    = (u0_ss * y_t**2 - 2.d0*u0_st * y_s*y_t + u0_tt * y_s**2  & 	    
	      + u0_s * (y_st*y_t - y_tt*y_s )			      &    
	      + u0_t * (y_st*y_s - y_ss*y_t ) )      / xjac**2        & 	
	    - xjac_x * (u0_s * y_t - u0_t * y_s)     / xjac**2
  u0_yy    = (u0_ss * x_t**2 - 2.d0*u0_st * x_s*x_t + u0_tt * x_s**2  & 	    
	      + u0_s * (x_st*x_t - x_tt*x_s )			      &    
	      + u0_t * (x_st*x_s - x_ss*x_t ) )      / xjac**2        & 	
	    - xjac_y * (- u0_s * x_t + u0_t * x_s )  / xjac**2
  u0_xy    = (- u0_ss * y_t*x_t - u0_tt * x_s*y_s		      &
	      + u0_st * (y_s*x_t  + y_t*x_s  )  		      &        
	      - u0_s  * (x_st*y_t - x_tt*y_s )  		      &    
	      - u0_t  * (x_st*y_s - x_ss*y_t )  )    / xjac**2        & 	
	    - xjac_x * (- u0_s * x_t + u0_t * x_s )  / xjac**2
  vv2	   = R**2 *  ( u0_x * u0_x + u0_y *u0_y  )
  
  ! --- Variable 3
  zj0_x    = (   y_t * zj0_s - y_s * zj0_t ) / xjac
  zj0_y    = ( - x_t * zj0_s + x_s * zj0_t ) / xjac
  zj0_xx   = (zj0_ss * y_t**2 - 2.d0*zj0_st * y_s*y_t + zj0_tt * y_s**2  & 	    
	      + zj0_s * (y_st*y_t - y_tt*y_s )			         &    
	      + zj0_t * (y_st*y_s - y_ss*y_t ) )       / xjac**2         & 	
	    - xjac_x * (zj0_s * y_t - zj0_t * y_s)     / xjac**2
  zj0_yy   = (zj0_ss * x_t**2 - 2.d0*zj0_st * x_s*x_t + zj0_tt * x_s**2  & 	    
	      + zj0_s * (x_st*x_t - x_tt*x_s )			         &    
	      + zj0_t * (x_st*x_s - x_ss*x_t ) )       / xjac**2         & 	
	    - xjac_y * (- zj0_s * x_t + zj0_t * x_s )  / xjac**2
  zj0_xy   = (- zj0_ss * y_t*x_t - zj0_tt * x_s*y_s		         &
	      + zj0_st * (y_s*x_t  + y_t*x_s  )  		         &        
	      - zj0_s  * (x_st*y_t - x_tt*y_s )  		         &    
	      - zj0_t  * (x_st*y_s - x_ss*y_t )  )     / xjac**2         & 	
	    - xjac_x * (- zj0_s * x_t + zj0_t * x_s )  / xjac**2
  
  ! --- Variable 4
  w0_x     = (   y_t * w0_s - y_s * w0_t ) / xjac
  w0_y     = ( - x_t * w0_s + x_s * w0_t ) / xjac
  w0_xx    = (w0_ss * y_t**2 - 2.d0*w0_st * y_s*y_t + w0_tt * y_s**2  & 	    
	      + w0_s * (y_st*y_t - y_tt*y_s )			      &    
	      + w0_t * (y_st*y_s - y_ss*y_t ) )      / xjac**2        & 	
	    - xjac_x * (w0_s * y_t - w0_t * y_s)     / xjac**2
  w0_yy    = (w0_ss * x_t**2 - 2.d0*w0_st * x_s*x_t + w0_tt * x_s**2  & 	    
	      + w0_s * (x_st*x_t - x_tt*x_s )			      &    
	      + w0_t * (x_st*x_s - x_ss*x_t ) )      / xjac**2        & 	
	    - xjac_y * (- w0_s * x_t + w0_t * x_s )  / xjac**2
  w0_xy    = (- w0_ss * y_t*x_t - w0_tt * x_s*y_s		      &
	      + w0_st * (y_s*x_t  + y_t*x_s  )  		      &        
	      - w0_s  * (x_st*y_t - x_tt*y_s )  		      &    
	      - w0_t  * (x_st*y_s - x_ss*y_t )  )    / xjac**2        & 	
	    - xjac_x * (- w0_s * x_t + w0_t * x_s )  / xjac**2
  
  ! --- Variable 5
  r0_x     = (   y_t * r0_s - y_s * r0_t ) / xjac
  r0_y     = ( - x_t * r0_s + x_s * r0_t ) / xjac
  r0_xx    = (r0_ss * y_t**2 - 2.d0*r0_st * y_s*y_t + r0_tt * y_s**2  & 	    
	      + r0_s * (y_st*y_t - y_tt*y_s )			      &    
	      + r0_t * (y_st*y_s - y_ss*y_t ) )      / xjac**2        & 	
	    - xjac_x * (r0_s * y_t - r0_t * y_s)     / xjac**2
  r0_yy    = (r0_ss * x_t**2 - 2.d0*r0_st * x_s*x_t + r0_tt * x_s**2  & 	    
	      + r0_s * (x_st*x_t - x_tt*x_s )			      &    
	      + r0_t * (x_st*x_s - x_ss*x_t ) )      / xjac**2        & 	
	    - xjac_y * (- r0_s * x_t + r0_t * x_s )  / xjac**2
  r0_xy    = (- r0_ss * y_t*x_t - r0_tt * x_s*y_s		      &
	      + r0_st * (y_s*x_t  + y_t*x_s  )  		      &        
	      - r0_s  * (x_st*y_t - x_tt*y_s )  		      &    
	      - r0_t  * (x_st*y_s - x_ss*y_t )  )    / xjac**2        & 	
	    - xjac_x * (- r0_s * x_t + r0_t * x_s )  / xjac**2
  r0_hat   = R**2 * r0
  r0_x_hat = 2.d0 * R * R_x  * r0 + R**2 * r0_x
  r0_y_hat = R**2 * r0_y
  
  ! --- Variable 6
  T_corr    = corr_neg_temp(T0) ! For use in eta(T), visco(T), ...
  T0_x      = (   y_t * T0_s  - y_s * T0_t ) / xjac
  T0_y      = ( - x_t * T0_s  + x_s * T0_t ) / xjac
  T0_xx     = (T0_ss * y_t**2 - 2.d0*T0_st * y_s*y_t + T0_tt * y_s**2	& 	    
	      + T0_s * (y_st*y_t - y_tt*y_s )			        &    
	      + T0_t * (y_st*y_s - y_ss*y_t ) )        / xjac**2        & 	
	     - xjac_x * (T0_s * y_t - T0_t * y_s)      / xjac**2
  T0_yy     = (T0_ss * x_t**2 - 2.d0*T0_st * x_s*x_t + T0_tt * x_s**2   & 	    
	      + T0_s * (x_st*x_t - x_tt*x_s )			        &    
	      + T0_t * (x_st*x_s - x_ss*x_t ) )        / xjac**2        & 	
	     - xjac_y * (- T0_s * x_t + T0_t * x_s )   / xjac**2
  T0_xy     = (- T0_ss * y_t*x_t - T0_tt * x_s*y_s		        &
	       + T0_st * (y_s*x_t  + y_t*x_s  )  		        &        
	       - T0_s  * (x_st*y_t - x_tt*y_s )  		        &    
	       - T0_t  * (x_st*y_s - x_ss*y_t )  )     / xjac**2        & 	
	       - xjac_x * (- T0_s * x_t + T0_t * x_s ) / xjac**2

  
  ! --- Variable 7
  Vpar0_x  = (   y_t * Vpar0_s - y_s * Vpar0_t ) / xjac
  Vpar0_y  = ( - x_t * Vpar0_s + x_s * Vpar0_t ) / xjac
  Vpar0_xx = (Vpar0_ss * y_t**2 - 2.d0*Vpar0_st * y_s*y_t + Vpar0_tt * y_s**2  &	 
	      + Vpar0_s * (y_st*y_t - y_tt*y_s )			       &	
	      + Vpar0_t * (y_st*y_s - y_ss*y_t ) )         / xjac**2           &      
	    - xjac_x * (Vpar0_s * y_t - Vpar0_t * y_s)     / xjac**2
  Vpar0_yy = (Vpar0_ss * x_t**2 - 2.d0*Vpar0_st * x_s*x_t + Vpar0_tt * x_s**2  &	 
	      + Vpar0_s * (x_st*x_t - x_tt*x_s )			       &	
	      + Vpar0_t * (x_st*x_s - x_ss*x_t ) )         / xjac**2           &      
	    - xjac_y * (- Vpar0_s * x_t + Vpar0_t * x_s )  / xjac**2
  Vpar0_xy = (- Vpar0_ss * y_t*x_t - Vpar0_tt * x_s*y_s 		       &
	      + Vpar0_st * (y_s*x_t  + y_t*x_s  )  		               &	    
	      - Vpar0_s  * (x_st*y_t - x_tt*y_s )  		               &	
	      - Vpar0_t  * (x_st*y_s - x_ss*y_t )  )       / xjac**2           &      
	    - xjac_x * (- Vpar0_s * x_t + Vpar0_t * x_s )  / xjac**2
  
  ! --- Variable 8
  rn0_x    = (   y_t * rn0_s - y_s * rn0_t ) / xjac
  rn0_y    = ( - x_t * rn0_s + x_s * rn0_t ) / xjac
  rn0_xx   = (rn0_ss * y_t**2 - 2.d0*rn0_st * y_s*y_t + rn0_tt * y_s**2  & 	    
	      + rn0_s * (y_st*y_t - y_tt*y_s )			         &    
	      + rn0_t * (y_st*y_s - y_ss*y_t ) )       / xjac**2         & 	
	    - xjac_x * (rn0_s * y_t - rn0_t * y_s)     / xjac**2
  rn0_yy   = (rn0_ss * x_t**2 - 2.d0*rn0_st * x_s*x_t + rn0_tt * x_s**2  & 	    
	      + rn0_s * (x_st*x_t - x_tt*x_s )			         &    
	      + rn0_t * (x_st*x_s - x_ss*x_t ) )       / xjac**2         & 	
	    - xjac_y * (- rn0_s * x_t + rn0_t * x_s )  / xjac**2
  rn0_xy   = (- rn0_ss * y_t*x_t - rn0_tt * x_s*y_s		         &
	      + rn0_st * (y_s*x_t  + y_t*x_s  )  		         &        
	      - rn0_s  * (x_st*y_t - x_tt*y_s )  		         &    
	      - rn0_t  * (x_st*y_s - x_ss*y_t )  )     / xjac**2         & 	
	    - xjac_x * (- rn0_s * x_t + rn0_t * x_s )  / xjac**2
  
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
  BB2	    = (F0*F0 + ps0_x * ps0_x + ps0_y * ps0_y )/R**2
  
  return

end subroutine ELM_build_neutral_variables







!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------- Compute the diffusivities and source ---------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
subroutine ELM_build_neutral_diffusivities_and_sources(element, nodes, xpoint2, xcase2, R_axis, Z_axis, psi_axis, psi_bnd, Z_xpoint, i_plane)
!DEC$ ATTRIBUTES FORCEINLINE :: ELM_build_neutral_diffusivities_and_sources

  ! --- Modules
  use mod_parameters    
  use basis_at_gaussian
  use phys_module
  use equation_variables
  use data_structure
  use diffusivities, only: get_dperp, get_zkperp
  use corr_neg
  
  implicit none
  
  ! --- Routine variables
  type (type_element)		:: element
  type (type_node)		:: nodes(n_vertex_max)
  logical			:: xpoint2
  integer			:: xcase2
  integer			:: i_plane
  real*8			:: R_axis, Z_axis, psi_axis, psi_bnd, Z_xpoint(2)
  
  real*8			:: coef_ion_1, coef_ion_2, coef_ion_3, S_ion_puiss
      
 
  ! --------------------------------------------
  ! --- Ionization cost for temperature equation
  ! --------------------------------------------
  coef_ion_3 = 2.738d-24*n_zero                   
  coef_ion_2 = 0.232d0                            
  coef_ion_1 = 1.886d-30*(n_zero)**(1.5d0)        
  S_ion_puiss = 3.9d-1                            
  ksiion = ksi_ion * n_zero
  
  S_ion = coef_ion_1*((coef_ion_3/corr_neg_temp(T0,(/1.d-5,0.3/)))**S_ion_puiss)*1/(coef_ion_2+coef_ion_3/corr_neg_temp(T0,(/1.d-5,0.3/)))*exp(-coef_ion_3/corr_neg_temp(T0,(/1.d-5,0.3/)))			   
  S_ion_T = coef_ion_1*exp(-coef_ion_3/corr_neg_temp(T0,(/1.d-5,0.3/)))*((coef_ion_3/corr_neg_temp(T0,(/1.d-5,0.3/)))**0.39d0)*1/(corr_neg_temp(T0,(/1.d-5,0.3/))*(coef_ion_2*corr_neg_temp(T0,(/1.d-5,0.3/))+coef_ion_3)**2.0d0) &   
  	   * (coef_ion_3*((coef_ion_2+1)*corr_neg_temp(T0,(/1.d-5,0.3/))+coef_ion_3)-0.39d0*(corr_neg_temp(T0,(/1.d-5,0.3/))*(coef_ion_2*corr_neg_temp(T0,(/1.d-5,0.3/))+coef_ion_3)))  		  
   
  
  ! ----------------------------
  ! --- Neutral diffusion coeffs
  ! ----------------------------
  Dn0x = D_neutral_x	   
  Dn0y = D_neutral_y	   
  Dn0p = D_neutral_p	   


  ! -------------------------------------------------------------------
  ! --- Heating, current and particle source (the same for all i_plane)
  ! -------------------------------------------------------------------
  if (i_plane .eq. 1) then
    ! --- Neutral density source
    phi        = 2.d0*PI * float(i_plane-1) / float(n_plane) / float(n_period)
    source_neutral = 0.d0			
    call neutral_source(ns_amplitude, ns_R, ns_Z, ns_phi, ns_radius, ns_sig, ns_length, x_g, y_g, phi, source_neutral)				    
  endif
  
  return

end subroutine ELM_build_neutral_diffusivities_and_sources

