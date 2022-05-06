module mod_elt_matrix_fft
contains

  ! --- Include all the routines directly for runtime efficiency
  INCLUDE "construct_variables.f90"
  INCLUDE "equations.f90"
  INCLUDE "equations_numm.f90"


  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------ Calculates the matrix contribution of one element ---------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------------------------------------------
  subroutine element_matrix_fft(element, nodes, xpoint2, xcase2, R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint, Z_xpoint, &
                                ELM, RHS, tid, ELM_p, ELM_n, ELM_k, ELM_kn, RHS_p, RHS_k,                               &
                                eq_g, eq_s, eq_t, eq_p, eq_ss, eq_st, eq_tt, delta_g_arg, delta_s_arg, delta_t_arg,     &
                                i_tor_min, i_tor_max, aux_nodes, ELM_pnn)

    ! --- Modules
    use equation_variables
    use constants
    use mod_parameters
    use data_structure
    use gauss
    use basis_at_gaussian
    use phys_module
    use tr_module 
    use profiles, only: interpolProf
    use diffusivities, only: get_dperp, get_zkperp
    use corr_neg
    use pellet_module
    use mod_elm_apply_fft
    use mod_bootstrap_functions
    use vacuum

    implicit none
    
    ! --- Structures
    type (type_element)        :: element
    type (type_node)           :: nodes(n_vertex_max)
    type (type_node), optional :: aux_nodes(n_vertex_max)

    type (type_surface_list)   :: flux_list

    ! --- Matrix elements and toroidal functions
    integer, intent(in)        :: tid
    integer, intent(in)        :: i_tor_min, i_tor_max
#define DIM0 n_tor*n_vertex_max*n_degrees*n_var
#define DIM1 n_plane
#define DIM2 1:n_vertex_max*n_var*n_degrees

    real*8, dimension (DIM0,DIM0)	     	:: ELM
    real*8, dimension (DIM0)	     		:: RHS
    real*8, dimension(DIM1,DIM2,DIM2)    	:: ELM_p, ELM_n, ELM_k, ELM_kn
    real*8, dimension(DIM1,DIM2)  	        :: RHS_p, RHS_k
    real*8, dimension(DIM1, DIM2, DIM2)         :: ELM_pnn

! The following buffers are not used by this model:
    real*8, dimension(n_plane,n_var,n_gauss,n_gauss) :: eq_g, eq_s, eq_t
    real*8, dimension(n_plane,n_var,n_gauss,n_gauss) :: eq_p
    real*8, dimension(n_plane,n_var,n_gauss,n_gauss) :: eq_ss, eq_st, eq_tt
    real*8, dimension(n_plane,n_var,n_gauss,n_gauss) :: delta_g_arg, delta_s_arg, delta_t_arg
    
    ! --- Indexes
    integer    :: i_ij, ij_tmp
    integer    :: i_kl, kl_tmp
    integer    :: ms, mt
    integer    :: i_plane
    integer    :: i_order,  i_vertex, i_tor
    integer    :: j_order,  j_vertex, j_tor
    integer    :: n_tor_loop, n_tor_loop2
    integer    :: index_ij, index_kl
    
    ! --- Routine variables (Xpoint and axis)
    logical    :: xpoint2
    integer    :: xcase2
    real*8     :: R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint(2), Z_xpoint(2)
    
    ! --- Integration weight
    real*8     :: wst
    logical    :: use_fft   
    integer    :: n_tor_start, n_tor_end 
        
    ! --- Initialise rhs and lhs terms
    rhs_tmp  = 0.d0; rhs_k_tmp  = 0.d0
    amat_tmp = 0.d0; amat_k_tmp = 0.d0; amat_n_tmp = 0.d0; amat_kn_tmp = 0.d0
    
    ! --- Matrix elements arrays
    ELM_p = 0.d0
    ELM_n = 0.d0
    ELM_k = 0.d0
    ELM_kn = 0.d0
    RHS_p = 0.d0
    RHS_k = 0.d0
    
    ELM = 0.d0; RHS = 0.d0
        
    ! --- Take time evolution parameters from phys_module
    theta = time_evol_theta
    !zeta  = time_evol_zeta
    ! change zeta for variable dt
    zeta  = time_evol_zeta * 2.0d0 * tstep / (tstep + tstep_prev)

    ! for cylinder geometry : epscyl = eps
    eps_cyl = 1.d0

    ! --- Do we need to use the FFT or non-FFT version?
    if ( (i_tor_min == 1) .and. (i_tor_max == n_tor) ) then
      ! In case of global matrix construction:
      use_fft = n_tor > 3 
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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!! Begin integration loop over Gaussian integration points !!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do ms =1,n_gauss
      do mt =1,n_gauss
      
    	wst = wgauss(ms)*wgauss(mt)

    	call ELM_build_RZ_and_Jacobians(element, nodes, ms, mt)

    	do i_plane =1,n_plane

    	  call ELM_build_variables(element, nodes, ms, mt, i_plane)
          
    	  call ELM_build_diffusivities_and_sources(element, nodes, xpoint2, xcase2, R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint, Z_xpoint, i_plane)

    	  ! --- Now the equations, first the RHS
    	  do i_vertex =1,n_vertex_max

    	    do i_order =1,n_degrees

    	      do i_tor =n_tor_start, n_tor_end

    		! --- Index in the ELM matrix	    
    		if (use_fft) then
    		  index_ij =       n_var*n_degrees*(i_vertex-1) +       n_var*(i_order-1) + 1
    		else
    		  index_ij = (n_tor_end - n_tor_start +1)*n_var*n_degrees*(i_vertex-1) + (n_tor_end - n_tor_start +1)*n_var*(i_order-1) + & 
                             i_tor - n_tor_start +1
    		endif
		
	        ! --- Build test functions (which we choose to be the basis functions)
		call ELM_build_basis_functions(element, nodes, ms, mt, i_plane, i_vertex, i_order, i_tor, &
					       v, v_s,  v_t,	    v_p,  v_x,  v_y,			  &
						  v_ss, v_tt, v_st, v_pp, v_xx, v_yy, v_xy		  )
		
		rhs_tmp   = 0.d0
		rhs_k_tmp = 0.d0
		call ELM_main_rhs_1  	(rhs_tmp, rhs_k_tmp)
		call ELM_main_rhs_2  	(rhs_tmp, rhs_k_tmp)
		call ELM_main_rhs_3  	(rhs_tmp, rhs_k_tmp)
		call ELM_main_rhs_4  	(rhs_tmp, rhs_k_tmp)
		call ELM_main_rhs_5  	(rhs_tmp, rhs_k_tmp)
		call ELM_main_rhs_6  	(rhs_tmp, rhs_k_tmp)
		call ELM_main_rhs_7  	(rhs_tmp, rhs_k_tmp)
		call ELM_main_rhs_2_numm(rhs_tmp, rhs_k_tmp)
		call ELM_main_rhs_5_numm(rhs_tmp, rhs_k_tmp)
		call ELM_main_rhs_6_numm(rhs_tmp, rhs_k_tmp)
		call ELM_main_rhs_7_numm(rhs_tmp, rhs_k_tmp)
    		

    		! --- Fill up the matrix
    		if (use_fft) then
    		  do i_ij =1,n_var
		    ij_tmp = index_ij + (i_ij-1)*(n_tor_end - n_tor_start +1)
		    RHS_p(i_plane,ij_tmp) = RHS_p(i_plane,ij_tmp) + rhs_tmp  (i_ij) * wst
		    RHS_k(i_plane,ij_tmp) = RHS_k(i_plane,ij_tmp) + rhs_k_tmp(i_ij) * wst
		  enddo
    		else
    		  do i_ij =1,n_var
		    ij_tmp = index_ij + (i_ij-1)*(n_tor_end - n_tor_start +1)
    		    RHS(ij_tmp) = RHS(ij_tmp) + (rhs_tmp(i_ij) + rhs_k_tmp(i_ij)) * wst
		  enddo
    		endif

    		! --- And the LHS (linearised part)
    		do j_vertex =1,n_vertex_max

    		  do j_order =1,n_degrees

    		    do j_tor =n_tor_start, n_tor_end 

    		      ! --- Index in the ELM matrix
    		      if (use_fft) then
    			index_kl =       n_var*n_degrees*(j_vertex-1) +       n_var*(j_order-1) + 1
    		      else
    			index_kl = (n_tor_end - n_tor_start +1)*n_var*n_degrees*(j_vertex-1) + & 
                                   (n_tor_end - n_tor_start +1)*n_var*(j_order-1) + j_tor - n_tor_start +1
    		      endif

		      ! --- Build basis functions
		      call ELM_build_basis_functions(element, nodes, ms, mt, i_plane, j_vertex, j_order, j_tor,  &
                                                     psi, psi_s,  psi_t,          psi_p,  psi_x,  psi_y,         &
				                          psi_ss, psi_tt, psi_st, psi_pp, psi_xx, psi_yy, psi_xy )
    		      
		      u    = psi; u_x    = psi_x; u_y    = psi_y; u_p    = psi_p; u_s    = psi_s; u_t    = psi_t
		      zj   = psi; zj_x   = psi_x; zj_y   = psi_y; zj_p   = psi_p; zj_s   = psi_s; zj_t   = psi_t
		      w    = psi; w_x    = psi_x; w_y    = psi_y; w_p    = psi_p; w_s    = psi_s; w_t    = psi_t
		      rho  = psi; rho_x  = psi_x; rho_y  = psi_y; rho_p  = psi_p; rho_s  = psi_s; rho_t  = psi_t
		      T    = psi; T_x    = psi_x; T_y    = psi_y; T_p    = psi_p; T_s    = psi_s; T_t    = psi_t
		      Vpar = psi; Vpar_x = psi_x; Vpar_y = psi_y; Vpar_p = psi_p; Vpar_s = psi_s; Vpar_t = psi_t
		      
		      u_ss    = psi_ss; u_tt	= psi_tt; u_st    = psi_st
		      zj_ss   = psi_ss; zj_tt	= psi_tt; zj_st   = psi_st
		      w_ss    = psi_ss; w_tt	= psi_tt; w_st    = psi_st
   		      rho_ss  = psi_ss; rho_tt  = psi_tt; rho_st  = psi_st
   		      T_ss    = psi_ss; T_tt	= psi_tt; T_st    = psi_st
   		      Vpar_ss = psi_ss; Vpar_tt = psi_tt; Vpar_st = psi_st
                      
		      u_xx    = psi_xx; u_yy    = psi_yy; u_xy    = psi_xy; u_pp    = psi_pp
		      zj_xx   = psi_xx; zj_yy   = psi_yy; zj_xy   = psi_xy; zj_pp   = psi_pp
		      w_xx    = psi_xx; w_yy    = psi_yy; w_xy    = psi_xy; w_pp    = psi_pp
		      rho_xx  = psi_xx; rho_yy  = psi_yy; rho_xy  = psi_xy; rho_pp  = psi_pp
		      T_xx    = psi_xx; T_yy    = psi_yy; T_xy    = psi_xy; T_pp    = psi_pp
		      Vpar_xx = psi_xx; Vpar_yy = psi_yy; Vpar_xy = psi_xy; Vpar_pp = psi_pp
                      
		      BB2_psi = 2.d0 * (psi_x * ps0_x + psi_y * ps0_y ) /R**2
   
    		      amat_tmp    = 0.d0
    		      amat_k_tmp  = 0.d0
    		      amat_n_tmp  = 0.d0
    		      amat_kn_tmp = 0.d0
		      call ELM_main_lhs_1     (amat_tmp, amat_k_tmp, amat_n_tmp, amat_kn_tmp)
    		      call ELM_main_lhs_2     (amat_tmp, amat_k_tmp, amat_n_tmp, amat_kn_tmp)
    		      call ELM_main_lhs_3     (amat_tmp, amat_k_tmp, amat_n_tmp, amat_kn_tmp)
    		      call ELM_main_lhs_4     (amat_tmp, amat_k_tmp, amat_n_tmp, amat_kn_tmp)
    		      call ELM_main_lhs_5     (amat_tmp, amat_k_tmp, amat_n_tmp, amat_kn_tmp)
    		      call ELM_main_lhs_6     (amat_tmp, amat_k_tmp, amat_n_tmp, amat_kn_tmp)
    		      call ELM_main_lhs_7     (amat_tmp, amat_k_tmp, amat_n_tmp, amat_kn_tmp)
    		      call ELM_main_lhs_2_numm(amat_tmp, amat_k_tmp, amat_n_tmp, amat_kn_tmp)
    		      call ELM_main_lhs_5_numm(amat_tmp, amat_k_tmp, amat_n_tmp, amat_kn_tmp)
    		      call ELM_main_lhs_6_numm(amat_tmp, amat_k_tmp, amat_n_tmp, amat_kn_tmp)
    		      call ELM_main_lhs_7_numm(amat_tmp, amat_k_tmp, amat_n_tmp, amat_kn_tmp)
		      
    		      ! --- Fill up the matrix
    		      if (use_fft) then
    		  	do i_ij =1,n_var
		  	  ij_tmp = index_ij + (i_ij-1)*(n_tor_end - n_tor_start +1)
    		  	  do i_kl =1,n_var
		  	    kl_tmp = index_kl + (i_kl-1)*(n_tor_end - n_tor_start +1)
    			    ELM_p (i_plane,ij_tmp,kl_tmp) = ELM_p (i_plane,ij_tmp,kl_tmp) + wst * amat_tmp   (i_ij,i_kl)
    			    ELM_k (i_plane,ij_tmp,kl_tmp) = ELM_k (i_plane,ij_tmp,kl_tmp) + wst * amat_k_tmp (i_ij,i_kl)
    			    ELM_n (i_plane,ij_tmp,kl_tmp) = ELM_n (i_plane,ij_tmp,kl_tmp) + wst * amat_n_tmp (i_ij,i_kl)
    			    ELM_kn(i_plane,ij_tmp,kl_tmp) = ELM_kn(i_plane,ij_tmp,kl_tmp) + wst * amat_kn_tmp(i_ij,i_kl)
		  	  enddo
		  	enddo
    		      else
    		  	do i_ij =1,n_var
		  	  ij_tmp = index_ij + (i_ij-1)*(n_tor_end - n_tor_start+1)
    		  	  do i_kl =1,n_var
		  	    kl_tmp = index_kl + (i_kl-1)*(n_tor_end - n_tor_start +1)
    			    ELM(ij_tmp,kl_tmp) = ELM(ij_tmp,kl_tmp) + (amat_tmp(i_ij,i_kl) + amat_k_tmp(i_ij,i_kl) + amat_n_tmp(i_ij,i_kl) + amat_kn_tmp(i_ij,i_kl)) * wst
		  	  enddo
		  	enddo
    		      endif


    		    
    		    enddo ! inner n_tor_loop
    		  enddo   ! inner n_degrees
    		enddo	  ! inner n_vertex_max

    	      enddo	  ! outer n_tor_loop
    	    enddo	  ! outer n_degrees
    	  enddo 	  ! outer n_vertex_max

    	enddo		  ! n_plane

      enddo		  ! n_gauss
    enddo		  ! n_gauss




    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!! Apply FFT !!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (use_fft) then
      call ELM_apply_fft(RHS, RHS_p, RHS_k, ELM, ELM_p, ELM_n, ELM_k, ELM_kn, tid)
    endif
    
    return
  end subroutine element_matrix_fft



end module mod_elt_matrix_fft
