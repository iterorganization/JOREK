!*****************************************************************************
!        Modifications of the stiffness matrix and RHS                       *
!*****************************************************************************
subroutine Chgmt_node(i_element, element,nodes,element_father,nodes_father,ELM,RHS,node_out)



 use mod_parameters
 use data_structure
 use mod_basisfunctions

implicit none

 type (type_element)   :: element
 type (type_element)   :: element_father
 type (type_node)      :: nodes(n_vertex_max)
 type (type_node)      :: nodes_father(n_vertex_max)

 real*8, dimension(n_vertex_max*n_degrees,&
                 n_vertex_max*n_degrees)        :: ELM, ELM_bis, ELM_tr, C_matrix

 real*8, dimension( n_vertex_max*n_degrees)     :: RHS,RHS_bis
 real*8, dimension(4,n_degrees)                 :: H, H_s, H_t, H_st
 real*8, dimension(2,4)                         :: c, dc_ds, dc_dt, d2c_dsdt
 real*8                                         :: lambda, mu
 integer, dimension(n_vertex_max)               :: pr, pos_node_constrained
 integer, dimension(n_vertex_max)               :: node_out
 integer, dimension(2,2)                        :: pos_parent, parent
 real*8                                         :: h_u, h_v,h_w
 integer                                        :: i_element,index_elm
 integer                                        :: pos1, pos2,n_constrained,i_constrained
 integer                                        :: i, j,k,l,n, p,prj
 



  do j = 1, n_vertex_max           
       node_out(j)=0
  enddo 
 
  n_constrained = 0
  i_constrained = 0
  pos_node_constrained = 0
 
  
  do j = 1, n_vertex_max
       
        if ((.not. nodes(j)%constrained)) then
	
          node_out(j)=element%vertex(j)    
        
        else

	 !*****************************************************************************
         ! Processing  "constrained nodes"
         !*****************************************************************************
   
       !if (i_constrained==0) then
         !write(*,*) '**************************************'
         !write(*,*) '*           changement nodes         *'
         !write(*,*) '**************************************'
       !endif
	    i_constrained = 1
	    n_constrained = n_constrained + 1
	    pos_node_constrained(n_constrained) = j                           ! Position of the constrained node in 'i_element'
	     
   
	    
	    parent(n_constrained,1) = nodes(j)%parents(1)
	    parent(n_constrained,2) = nodes(j)%parents(2)
	    do k = 1, n_vertex_max   
	         
		 if(element%vertex(k) == parent(n_constrained,1)) then
		      pos_parent(n_constrained,1) = k
		      pos_parent(n_constrained,2) = j		      
		   elseif (element%vertex(k) == parent(n_constrained,2)) then
	              parent(n_constrained,1) = nodes(j)%parents(2)
	              parent(n_constrained,2) = nodes(j)%parents(1)	
		      pos_parent(n_constrained,1) = k                          !position of the node that belongs to        "i_element"
		      pos_parent(n_constrained,2) = j                          !position of the node that  not  belongs to  "i_element"      
		 endif    
       
	    enddo    
    
       node_out(j)=parent(n_constrained,2)	    
       endif	    
  enddo




  

!*****************************************************************************
!                  Connectivity matrix C_matrix
!*****************************************************************************
 
  if(i_constrained == 1) then

       C_matrix = 0.
	
       do i = 1, n_vertex_max*n_degrees
	
	    C_matrix(i,i) = 1.
	
       enddo


       do j = 1, n_constrained 
  	    
	    prj = pr(pos_node_constrained(j))
        
            !******************************************************************
       	    ! Coefficients of the linear relation between the values at       *
            ! constrained node  and those of both parents                     *
	    !******************************************************************
	    
          
            lambda = nodes(pos_node_constrained(j))%ref_lambda
	    mu     = nodes(pos_node_constrained(j))%ref_mu
	    index_elm  = nodes(pos_node_constrained(j))%parent_elem
 
	    call basisfunctions(lambda, mu, H, H_s, H_t, H_st)
	    
          
	          	   

	    h_u =1
	    h_v =1
	    h_w =h_u*h_v

          
            do k = 1, n_vertex_max           
                 
		 pr(k) = element_father%vertex(k)		 
 
		 do p = 1, 2
		      if(pr(k)==parent(j,p)) then
		       
        	           do l = 1, n_degrees
                                c(p,l) 	 	= (H(k,l)*element_father%size(k,l))
		                dc_ds(p,l) 	= (H_s(k,l)*element_father%size(k,l)) / (3.*h_u)
		                dc_dt(p,l)	= (H_t(k,l)*element_father%size(k,l)) / (3.*h_v)
		                d2c_dsdt(p,l)	= (H_st(k,l)*element_father%size(k,l))/ (9.*h_w)

	                   enddo 
	              endif 
                   
		 enddo
            enddo 



	    
	    
!******************************************************************************************************************************
!         Modifications of the stiffness matrix and RHS                                                                       *
!******************************************************************************************************************************
	    	    	    	   
   
	     

  
	    Pos1 = (pos_parent(j,1)-1)*n_degrees +1                       ! Position of parent node that is in element
            Pos2 = (pos_parent(j,2)-1)*n_degrees +1                       ! Position of parent node that is outside of element
            

          
	    
	         C_matrix(Pos2,:)   = 0.
		 C_matrix(Pos2+1,:) = 0.
	         C_matrix(Pos2+2,:) = 0.
		 C_matrix(Pos2+3,:) = 0.
		 
	
          
	    
	    do k = 1, n_degrees
	         
	    	 C_matrix(Pos2,Pos1+k-1) = c(1,k) 		 
	    	 C_matrix(Pos2,Pos2+k-1) = c(2,k) 		
	    	 
	    	 C_matrix(Pos2+1,Pos1+k-1) = dc_ds(1,k)  	
	    	 C_matrix(Pos2+1,Pos2+k-1) = dc_ds(2,k) 		 
	    	     
	    	 C_matrix(Pos2+2,Pos1+k-1) = dc_dt(1,k) 	
	    	 C_matrix(Pos2+2,Pos2+k-1) = dc_dt(2,k) 	

            	 C_matrix(Pos2+3,Pos1+k-1) = d2c_dsdt(1,k)
	    	 C_matrix(Pos2+3,Pos2+k-1) = d2c_dsdt(2,k) 
                  
            enddo            
     
       enddo

 

        
     

       ELM_bis = 0.
       ELM_tr  = 0.
       RHS_bis = 0.  
            

            !***********************************************************************
            !Right multprlication of the stiffness matrix  by C_matrix             *
            !***********************************************************************

       do i = 1, n_vertex_max*n_degrees
            do k = 1, n_vertex_max*n_degrees
                 do n = 1, n_vertex_max*n_degrees
                      ELM_bis(i,k) = ELM_bis(i,k) + ELM(i,n)*C_matrix(n,k)
                 enddo     
            enddo
       enddo

            !*************************************************************************
            !Left multprlication of the stiffness matrix by the transposed C_matrix  *
            !*************************************************************************


       do i = 1, n_vertex_max*n_degrees
            do k = 1, n_vertex_max*n_degrees
                 do n = 1, n_vertex_max*n_degrees
                      ELM_tr(i,k) = ELM_tr(i,k) + C_matrix(n,i)*ELM_bis(n,k)
                 enddo     
            enddo
       enddo


            !************************************************************************
            !Left multprlication of the load vector by the transposed C_matrix      *
            !************************************************************************


       do i = 1, n_vertex_max*n_degrees
            do k = 1, n_vertex_max*n_degrees
		 RHS_bis(i) = RHS_bis(i) + C_matrix(k,i)*RHS(k) 	 
	    enddo  
       enddo	    

       ELM= ELM_tr
       RHS = RHS_bis
      
      
      
  endif
   

  return


end subroutine Chgmt_node
