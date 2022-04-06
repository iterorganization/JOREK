!**************************************************************************
! subroutine to add a node to the node list
! P1        : One of the parent nodes (index)
! P2        : Second parent node index
! i12 == 0  : It is a (new) frontiere node
! i12 == +1 : It is a new contrainte node
! i12 == -1 : It is a new node but not contrainte
! inew     -> index of the new node
!**************************************************************************


subroutine Ref_Add_Node(node_list, element_list,lambda, mu, iref,iside, i12, inew,counter)
 use data_structure
 use mod_parameters
 use mod_basisfunctions
 
 implicit none

 type (type_node_list)    :: node_list
 type (type_element_list) :: element_list 

 real*8, dimension(n_vertex_max,n_degrees)	:: H, H_s, H_t, H_st
 real*8, dimension(2)				:: P, dP_ds, dP_dt, d2P_dsdt, u, v, w 
 real*8						:: Somme_X, Somme_Y, dX_ds, dY_ds, dX_dt, dY_dt, &
                 				   d2X_dsdt, d2Y_dsdt
 real*8, dimension(n_tor)			:: Psi, dPsi_ds,dPsi_dt, d2Psi_dsdt						      
 real*8						:: lambda, mu
 real*8			                 	:: h_u, h_v, h_w
 integer, dimension(n_vertex_max)		:: pr 
 integer 					:: i12, inew, iref,iside,idir,np,P1,P2,counter
 integer                                        :: i,j,k,l,i_tor,i_var

!**********************************************************************************
! Determination of spatial variables at the new node of local coordinates         *
! (s,t) = (lambda,mu) on the element iref                                        *
!**********************************************************************************

  np   = node_list%n_nodes
  inew = np+1
     
    if (i12==0) then ! It is a boundary node

       call Ref_boundary_node(node_list, element_list, lambda, mu, iref, inew)

    
    else             ! It is not a boundary node
    
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
  end if
   

!****************************************************************************
!                  Genealogy of the new nodes
!****************************************************************************    

  P1 = 0
  P2 = 0

  
  if(mu==0.) then
       P1 = element_list%element(iref)%vertex(1) 
       P2 = element_list%element(iref)%vertex(2) 
  end if
   
  if(lambda==1.) then
       P1 = element_list%element(iref)%vertex(2) 
       P2 = element_list%element(iref)%vertex(3) 
  end if	   	

  if(mu==1.) then
       P1 = element_list%element(iref)%vertex(3) 
       P2 = element_list%element(iref)%vertex(4) 
  end if                	    

  if(lambda==0.) then
       P1 = element_list%element(iref)%vertex(4) 
       P2 = element_list%element(iref)%vertex(1) 
  end if     
         	              

  node_list%node(inew)%parents(1)  = P1
  node_list%node(inew)%parents(2)  = P2
  node_list%node(inew)%parent_elem = iref	! Save the element which contains P1 and P2
  node_list%node(inew)%ref_lambda  = lambda     ! Local coordinates of the new node in the element iref
  node_list%node(inew)%ref_mu	   = mu	
  
  if (i12 .eq. +1) then                                         ! Constrained Node
       
         node_list%node(inew)%constrained = .true.
       
       
  elseif (i12 .eq. 0) then                          	        ! boundary nodes
         node_list%node(inew)%boundary     =  2
         node_list%node(inew)%constrained  = .false.
      
       
  elseif (i12 .eq.-1) then                                       !center node
    
         node_list%node(inew)%constrained = .false.
      
      
  endif
   
         node_list%n_nodes = node_list%n_nodes + 1
         counter =counter +1
        
         element_list%element(iref)%contain_node(counter)= inew
  return


end subroutine Ref_Add_Node
