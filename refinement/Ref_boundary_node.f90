!********************************************************************************************
! function to determine the position and vectors of a new boundary node, given the position *
! of the two parent nodes                                                                   *
!********************************************************************************************

subroutine Ref_boundary_node(node_list, element_list, lambda,mu,iref, inew)

 use mod_parameters
 use data_structure
 use mod_basisfunctions

implicit none

 type (type_node_list)    :: node_list
 type (type_element_list) :: element_list

 real*8, dimension(n_vertex_max,n_degrees)	:: H, H_s, H_t, H_st
 real*8, dimension(2)				:: P, dP_ds, dP_dt, d2P_dsdt, u, v, w 
 real*8 					:: lambda, mu 
 real*8						:: Somme_X, Somme_Y, dX_ds, dY_ds, dX_dt, dY_dt, &
                  				   d2X_dsdt, d2Y_dsdt
 real*8                                         :: h_u,h_v,h_w
 real*8,  dimension(n_tor)			:: Psi, dPsi_ds,dPsi_dt, d2Psi_dsdt
 integer, dimension(n_vertex_max)		:: pr	
 integer                                        :: iref, inew					   
 integer					:: i,j,k,l,i_tor,i_var


  do j = 1, n_vertex_max	  
       pr(j) = element_list%element(iref)%vertex(j)
  end do 

 
       call BasisFunctions(lambda,mu,H,H_s,H_t,H_st)

       Somme_X  = 0.
       Somme_Y  = 0.
       dX_ds    = 0.
       dy_ds    = 0.
       dx_dt    = 0.
       dy_dt    = 0.
       d2X_dsdt = 0.
       d2Y_dsdt = 0.
       
          

       
       do k = 1, n_vertex_max	     

            do l = 1, n_degrees
	    
	         Somme_X = Somme_X + node_list%node(pr(k))%x(1,l,1) * H(k,l)	 &
		     * element_list%element(iref)%size(k,l)
	         Somme_Y = Somme_Y + node_list%node(pr(k))%x(1,l,2) * H(k,l)	 &
             	     * element_list%element(iref)%size(k,l)
                
	         dX_ds = dx_ds + node_list%node(pr(k))%x(1,l,1) * H_s(k,l)    &
		     * element_list%element(iref)%size(k,l)
	         dy_ds = dy_ds + node_list%node(pr(k))%x(1,l,2) * H_s(k,l)    &
             	     * element_list%element(iref)%size(k,l)	       
            
	         dx_dt = dx_dt + node_list%node(pr(k))%x(1,l,1) * H_t(k,l)    &
		     * element_list%element(iref)%size(k,l)
	         dy_dt = dy_dt + node_list%node(pr(k))%x(1,l,2) * H_t(k,l)    &
             	     * element_list%element(iref)%size(k,l)
             	       
	         d2X_dsdt = d2X_dsdt + node_list%node(pr(k))%x(1,l,1) * H_st(k,l)    &
		     * element_list%element(iref)%size(k,l)
	         d2Y_dsdt = d2Y_dsdt + node_list%node(pr(k))%x(1,l,2) * H_st(k,l)    &
             	     * element_list%element(iref)%size(k,l)	       
            enddo
       enddo  
       
           
        !**********************
        !   Position          *
        !**********************

       P(1)        = Somme_X   !position of shouted point 
       P(2)        = Somme_Y   !position of shouted point

       
       dP_ds(1)    = dX_ds     !vecteur "U"	     
       dP_ds(2)    = dy_ds     !vecteur "U"   
 
       dP_dt(1)    = dx_dt     !vecteur "V"   
       dP_dt(2)    = dy_dt     !vecteur "V"   
       
       d2P_dsdt(1) = d2X_dsdt  !vecteur "W"
       d2P_dsdt(2) = d2Y_dsdt  !vecteur "W"

    
        !**********************
        !      Vector "U"     *
        !**********************
         
       u(1) = (1./3.)*dX_ds     !vecteur "U"
       u(2) = (1./3.)*dy_ds     !vecteur "U"		       
       h_u  =1.
       
    
       !**********************
       !      Vector "V"     *
       !**********************
       
       v(1) = (1./3.)*dx_dt     !vecteur "V"   
       v(2) = (1./3.)*dy_dt     !vecteur "V" 
       h_v  =1.


        !**********************
        !      Vector "W"     *
        !**********************
       
       w(1) = ((1./3.)**2)*d2X_dsdt  !vecteur "W"
       w(2) = ((1./3.)**2)*d2Y_dsdt  !vecteur "W"
       h_w  = h_u*h_v 
     

       
      ! write(*,*) 'P',inew, P(1),P(2)
      ! write(*,*) 'u',inew, u(1),u(2)
      ! write(*,*) 'v',inew, v(1),v(2) 
      ! write(*,*) 'w',inew, w(1),w(2)
      ! write(*,*) ' '
          
      

      !***********************************************************
      !  Position and and vectors of New nodes  (P, u, v, w)     *
      !***********************************************************

         
       do j = 1, 2 !  directions "s" and "t"				     
            node_list%node(inew)%x(1,1,j) = P(j)	       
            node_list%node(inew)%x(1,2,j) = u(j) 
            node_list%node(inew)%x(1,3,j) = v(j) 
            node_list%node(inew)%x(1,4,j) = w(j)
                
       end do

     !****************************************
     !     update values                     *
     !****************************************
    
    do i_var=1,n_var 
      
       Psi = 0.
       dPsi_ds = 0.
       dPsi_dt = 0.
       d2Psi_dsdt = 0.

       do i_tor = 1, n_tor
            
            do k = 1, n_vertex_max	     

               	do l = 1, n_degrees      
                 

                 Psi(i_tor) = Psi(i_tor) + node_list%node(pr(k))%values(i_tor,l,i_var)* H(k,l) &
                     *element_list%element(iref)%size(k,l)
       
                 dPsi_ds(i_tor) = dPsi_ds(i_tor) + node_list%node(pr(k))%values(i_tor,l,i_var) * H_s(k,l) &
                     *element_list%element(iref)%size(k,l)

                 dPsi_dt(i_tor) = dPsi_dt(i_tor) + node_list%node(pr(k))%values(i_tor,l,i_var) * H_t(k,l) &
                     *element_list%element(iref)%size(k,l)
           
                 d2Psi_dsdt(i_tor) = d2Psi_dsdt(i_tor) + node_list%node(pr(k))%values(i_tor,l,i_var) &
		     * H_st(k,l) *element_list%element(iref)%size(k,l)	     
 
 		
               	end do 
		 
            end do 
	    		     
       
       
       node_list%node(inew)%values(i_tor,1,i_var)	=  Psi(i_tor)
       node_list%node(inew)%values(i_tor,2,i_var) 	= (dPsi_ds(i_tor)) / (3.*h_u)
       node_list%node(inew)%values(i_tor,3,i_var)	= (dPsi_dt(i_tor)) / (3.*h_v)
       node_list%node(inew)%values(i_tor,4,i_var)	= (d2Psi_dsdt(i_tor)) / (9.*h_w)	    
 
      
       end do
  
    enddo !(i_var)  

  return
 
end subroutine Ref_boundary_node
