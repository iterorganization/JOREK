subroutine update_values(element_list, node_list, rhs_vec)
!-----------------------------------------------------------------------
! subroutine adds the delta_values in RHS to the values in the node_list
!-----------------------------------------------------------------------

use data_structure, only: type_element_list, type_node_list, type_RHS
use phys_module, only: keep_n0_const, treat_axis
use mod_parameters, only: n_tor, n_var, n_degrees, n_vertex_max
use mod_basisfunctions
use mod_axis_treatment

implicit none

! --- Routine parameters
type (type_element_list), intent(in) :: element_list
type (type_node_list),    intent(inout) :: node_list
type(type_RHS),           intent(inout) :: rhs_vec

! --- local variables
real*8, dimension(4,n_degrees)	 :: H, H_s, H_t, H_st
real*8	:: lambda, mu
real*8, dimension(n_tor) :: Psi, dPsi_ds,dPsi_dt, d2Psi_dsdt
real*8, dimension(n_tor) :: Delt,Delt_ds,Delt_dt,Delt_dsdt
real*8	:: h_u, h_v, h_w   
integer, dimension(n_vertex_max)  :: Pr
integer, dimension(2)    :: parent
integer :: index_elm,l,i_tor,ivar
integer :: i, j, k, in, index_node, index, i_tor_min
integer :: id
real*8  :: stored_dofs(1:n_degrees*n_var*n_tor), new_dofs(1:4), old_dofs(1:4)

i_tor_min = 1
if ( keep_n0_const ) i_tor_min = 2 ! Keep equilibrium unchanged during the run

do i = 1, node_list%n_nodes
  if((.not. node_list%node(i)%constrained) ) then

  ! We need to transform the new dof to old ones on the axis.
  ! The respective RHS entries on the axis (which are shared)
  ! are updated during the transformation. At the end of the loop,
  ! we recover RHS entries so that they same can be used for all the axis nodes.
    if(treat_axis .and. node_list%node(i)%axis_node)then
      do k=1,n_var
        do in=1,n_tor
          do j=1,n_degrees
             index_node = node_list%node(i)%index(j)
             index = n_tor*n_var * (index_node - 1) + n_tor*(k-1) + in
             new_dofs(j) = rhs_vec%val(index)
             id = n_degrees*(n_tor)*(k-1) + n_degrees*(in-1) + j
             stored_dofs(id) = rhs_vec%val(index)
          enddo
          call new_to_old_dofs_on_the_axis(node_list, i, new_dofs, old_dofs)
          do j=1,n_degrees
             index_node = node_list%node(i)%index(j)
             index = n_tor*n_var * (index_node - 1) + n_tor*(k-1) + in
             rhs_vec%val(index) = old_dofs(j)
          enddo
        enddo
      enddo
    endif       

    do j=1,n_degrees
  
      index_node = node_list%node(i)%index(j)
  
#ifdef JECCD
  ! the n=0 component of eccd current should never be frozen when keep_n0_const=true
      do k=1,n_var-1
        do in=i_tor_min,n_tor
          index = n_tor*n_var * (index_node - 1) + n_tor*(k-1) + in
          if (index .gt. 0) then
            node_list%node(i)%values(in,j,k) = node_list%node(i)%values(in,j,k)+ rhs_vec%val(index)
            node_list%node(i)%deltas(in,j,k) = rhs_vec%val(index)
          endif  ! index gt 0
        enddo  ! in over n_tor
      enddo !k over nvar
  
  ! for final variable is the eccd current
      do in=1,n_tor
        index = n_tor*n_var * (index_node - 1) + n_tor*(n_var-1) + in
        if (index .gt. 0) then
          node_list%node(i)%values(in,j,n_var) = node_list%node(i)%values(in,j,n_var)+rhs_vec%val(index)
          node_list%node(i)%deltas(in,j,n_var) = rhs_vec%val(index)
        endif  ! index gt 0
      enddo  ! in over n_tor
  
#else
      do k=1,n_var
        do in=i_tor_min,n_tor
          index = n_tor*n_var * (index_node - 1) + n_tor*(k-1) + in
          if (index .gt. 0) then
            node_list%node(i)%values(in,j,k) = node_list%node(i)%values(in,j,k) + rhs_vec%val(index)
            node_list%node(i)%deltas(in,j,k) = rhs_vec%val(index)
          endif  ! index gt 0
        enddo  ! in over n_tor
      enddo !k over nvar
#endif
    enddo
   
    ! recover RHS entries.
    if(treat_axis .and. node_list%node(i)%axis_node)then
      do k=1,n_var
        do in=i_tor_min,n_tor
          do j=1,n_degrees
            index_node = node_list%node(i)%index(j)
            index = n_tor*n_var * (index_node - 1) + n_tor*(k-1) + in
            id = n_degrees*(n_tor)*(k-1) + n_degrees*(in-1) + j
            rhs_vec%val(index) = stored_dofs(id)
          enddo
        enddo
      enddo
    endif
    
  endif
!   write(*,'(i5,20e12.4)') i,node_list%node(i)%values(1,:,2),node_list%node(i)%values(2,:,2)

enddo
  !stop
do i = 1, node_list%n_nodes
  if((node_list%node(i)%constrained) ) then   

    lambda = node_list%node(i)%ref_lambda
	  mu     = node_list%node(i)%ref_mu
	  index_elm = node_list%node(i)%parent_elem
	  parent(1) = node_list%node(i)%parents(1)
	  parent(2) = node_list%node(i)%parents(2)
        
	  call basisfunctions(lambda, mu, H, H_s, H_t, H_st)
	   
    do j = 1, n_vertex_max           
      Pr(j) = element_list%element(index_elm)%vertex(j)    
    enddo 
      
	  h_u =1.
	  h_v =1. 
	  h_w =h_u*h_v 

#ifdef JECCD
! the n=0 component of eccd current should never be frozen when keep_n0_const=true
! this bit of code separates out that final equation.
! update values and deltas for first n_var-1 variables normally
    do ivar=1,n_var-1
      Psi = 0.
      dPsi_ds = 0.
      dPsi_dt = 0.
      d2Psi_dsdt = 0.

      Delt=0.
      Delt_ds=0.
      Delt_dt=0.
      Delt_dsdt=0.

      do i_tor = i_tor_min, n_tor
        do k = 1, n_vertex_max        
          Pr(k) = element_list%element(index_elm)%vertex(k)
          if((Pr(k)==parent(1)).or.(Pr(k)==parent(2))) then
            do l = 1, n_degrees
              !  Values      *
              Psi(i_tor) = Psi(i_tor) +node_list%node(Pr(k))%values(i_tor,l,ivar)* H(k,l) &
                     *element_list%element(index_elm)%size(k,l)
              dPsi_ds(i_tor) = dPsi_ds(i_tor) +node_list%node(Pr(k))%values(i_tor,l,ivar) * H_s(k,l) &
                     *element_list%element(index_elm)%size(k,l)
              dPsi_dt(i_tor) = dPsi_dt(i_tor) +node_list%node(Pr(k))%values(i_tor,l,ivar) * H_t(k,l) &
                     *element_list%element(index_elm)%size(k,l)
              d2Psi_dsdt(i_tor) = d2Psi_dsdt(i_tor) +node_list%node(Pr(k))%values(i_tor,l,ivar) &
                     * H_st(k,l) *element_list%element(index_elm)%size(k,l)     

              !  Deltas      *
              Delt(i_tor) = Delt(i_tor) +node_list%node(Pr(k))%deltas(i_tor,l,ivar)* H(k,l) &
                     *element_list%element(index_elm)%size(k,l)
              Delt_ds(i_tor) = Delt_ds(i_tor) +node_list%node(Pr(k))%deltas(i_tor,l,ivar) * H_s(k,l) &
                     *element_list%element(index_elm)%size(k,l)
              Delt_dt(i_tor) =  Delt_dt(i_tor) +node_list%node(Pr(k))%deltas(i_tor,l,ivar) * H_t(k,l) &
                     *element_list%element(index_elm)%size(k,l)
              Delt_dsdt(i_tor) = Delt_dsdt(i_tor)  +node_list%node(Pr(k))%deltas(i_tor,l,ivar) &
                     * H_st(k,l) *element_list%element(index_elm)%size(k,l)     
            enddo !(l)
          endif
        enddo !(k)                

        ! Values      *
        node_list%node(i)%values(i_tor,1,ivar)   = (Psi(i_tor))
        node_list%node(i)%values(i_tor,2,ivar)   = (dPsi_ds(i_tor)) / (3.*h_u)
        node_list%node(i)%values(i_tor,3,ivar)   = (dPsi_dt(i_tor)) / (3.*h_v)
        node_list%node(i)%values(i_tor,4,ivar)   = (d2Psi_dsdt(i_tor)) / (9.*h_w)        

        ! Deltas      *
        node_list%node(i)%deltas(i_tor,1,ivar)   = (Delt(i_tor) )
        node_list%node(i)%deltas(i_tor,2,ivar)   = (Delt_ds(i_tor)) / (3.*h_u)
        node_list%node(i)%deltas(i_tor,3,ivar)   = (Delt_dt(i_tor))/ (3.*h_v)
        node_list%node(i)%deltas(i_tor,4,ivar)   = (Delt_dsdt(i_tor))/ (9.*h_w)
      enddo!(i_tor)
    enddo !(ivar) 

!final eccd current variable, update values and deltas with i_tor_min (above) always
!equal to 1

    Psi = 0.
    dPsi_ds = 0.
    dPsi_dt = 0.
    d2Psi_dsdt = 0.

    Delt=0.
    Delt_ds=0.
    Delt_dt=0.
    Delt_dsdt=0.

    do i_tor = 1, n_tor
      do k = 1, n_vertex_max
        Pr(k) = element_list%element(index_elm)%vertex(k)
        if((Pr(k)==parent(1)).or.(Pr(k)==parent(2))) then
          do l = 1, n_degrees
            ! Values
            Psi(i_tor) = Psi(i_tor)+node_list%node(Pr(k))%values(i_tor,l,n_var)* H(k,l) &
                     *element_list%element(index_elm)%size(k,l)
            dPsi_ds(i_tor) = dPsi_ds(i_tor)+node_list%node(Pr(k))%values(i_tor,l,n_var) * H_s(k,l) &
                     *element_list%element(index_elm)%size(k,l)
            dPsi_dt(i_tor) = dPsi_dt(i_tor)+node_list%node(Pr(k))%values(i_tor,l,n_var) * H_t(k,l) &
                     *element_list%element(index_elm)%size(k,l)
            d2Psi_dsdt(i_tor) = d2Psi_dsdt(i_tor)+node_list%node(Pr(k))%values(i_tor,l,n_var) &
                     * H_st(k,l) *element_list%element(index_elm)%size(k,l)

            ! Deltas 
            Delt(i_tor) = Delt(i_tor)+node_list%node(Pr(k))%deltas(i_tor,l,n_var)* H(k,l) &
                     *element_list%element(index_elm)%size(k,l)
            Delt_ds(i_tor) = Delt_ds(i_tor)+node_list%node(Pr(k))%deltas(i_tor,l,n_var) * H_s(k,l) &
                     *element_list%element(index_elm)%size(k,l)
            Delt_dt(i_tor) =  Delt_dt(i_tor)+node_list%node(Pr(k))%deltas(i_tor,l,n_var) * H_t(k,l) &
                     *element_list%element(index_elm)%size(k,l)
            Delt_dsdt(i_tor) = Delt_dsdt(i_tor)+node_list%node(Pr(k))%deltas(i_tor,l,n_var) &
                     * H_st(k,l) *element_list%element(index_elm)%size(k,l)
          enddo !(l)
        endif
      enddo !(k)                

      ! Values  2 
      node_list%node(i)%values(i_tor,1,n_var)   = (Psi(i_tor))
      node_list%node(i)%values(i_tor,2,n_var)   = (dPsi_ds(i_tor)) / (3.*h_u)
      node_list%node(i)%values(i_tor,3,n_var)   = (dPsi_dt(i_tor)) / (3.*h_v)
      node_list%node(i)%values(i_tor,4,n_var)   = (d2Psi_dsdt(i_tor)) / (9.*h_w)

      ! Deltas  2   *
      node_list%node(i)%deltas(i_tor,1,n_var)   = (Delt(i_tor) )
      node_list%node(i)%deltas(i_tor,2,n_var)   = (Delt_ds(i_tor)) / (3.*h_u)
      node_list%node(i)%deltas(i_tor,3,n_var)   = (Delt_dt(i_tor))/ (3.*h_v)
      node_list%node(i)%deltas(i_tor,4,n_var)   = (Delt_dsdt(i_tor))/ (9.*h_w)
    enddo!(i_tor)

#else

    !***************************************************
    !     update values and deltas                     *
    !***************************************************
    do ivar=1,n_var

      Psi = 0.
	    dPsi_ds = 0.
	    dPsi_dt = 0.
	    d2Psi_dsdt = 0.

      Delt=0.
      Delt_ds=0.
      Delt_dt=0.
      Delt_dsdt=0.  
  
      do i_tor = i_tor_min, n_tor     
        do k = 1, n_vertex_max
          Pr(k) = element_list%element(index_elm)%vertex(k)    
          if((Pr(k)==parent(1)).or.(Pr(k)==parent(2))) then
            
            do l = 1, n_degrees
              !  Values

              Psi(i_tor) = Psi(i_tor) + node_list%node(Pr(k))%values(i_tor,l,ivar)* H(k,l) &
                     *element_list%element(index_elm)%size(k,l)
       
              dPsi_ds(i_tor) = dPsi_ds(i_tor) + node_list%node(Pr(k))%values(i_tor,l,ivar) * H_s(k,l) &
                     *element_list%element(index_elm)%size(k,l)

              dPsi_dt(i_tor) = dPsi_dt(i_tor) + node_list%node(Pr(k))%values(i_tor,l,ivar) * H_t(k,l) &
                     *element_list%element(index_elm)%size(k,l)
           
              d2Psi_dsdt(i_tor) = d2Psi_dsdt(i_tor) + node_list%node(Pr(k))%values(i_tor,l,ivar) &
		                 * H_st(k,l) *element_list%element(index_elm)%size(k,l)	    
                           
              !  Deltas
 		      
              Delt(i_tor) = Delt(i_tor) + node_list%node(Pr(k))%deltas(i_tor,l,ivar)* H(k,l) &
                     *element_list%element(index_elm)%size(k,l)
       
              Delt_ds(i_tor) = Delt_ds(i_tor) + node_list%node(Pr(k))%deltas(i_tor,l,ivar) * H_s(k,l) &
                     *element_list%element(index_elm)%size(k,l)

              Delt_dt(i_tor) =  Delt_dt(i_tor) + node_list%node(Pr(k))%deltas(i_tor,l,ivar) * H_t(k,l) &
                     *element_list%element(index_elm)%size(k,l)
           
              Delt_dsdt(i_tor) = Delt_dsdt(i_tor)  + node_list%node(Pr(k))%deltas(i_tor,l,ivar) &
		                 * H_st(k,l) *element_list%element(index_elm)%size(k,l)	     

            enddo !(l)
          endif
	      enddo !(k)                
        !  Values

        node_list%node(i)%values(i_tor,1,ivar)	= (Psi(i_tor))
        node_list%node(i)%values(i_tor,2,ivar) 	= (dPsi_ds(i_tor)) / (3.*h_u)
        node_list%node(i)%values(i_tor,3,ivar)	= (dPsi_dt(i_tor)) / (3.*h_v)
        node_list%node(i)%values(i_tor,4,ivar)	= (d2Psi_dsdt(i_tor)) / (9.*h_w)	    
 
       !  Deltas
        node_list%node(i)%deltas(i_tor,1,ivar)	= (Delt(i_tor) )
        node_list%node(i)%deltas(i_tor,2,ivar) 	= (Delt_ds(i_tor)) / (3.*h_u)
        node_list%node(i)%deltas(i_tor,3,ivar)	= (Delt_dt(i_tor))/ (3.*h_v)
        node_list%node(i)%deltas(i_tor,4,ivar)	= (Delt_dsdt(i_tor))/ (9.*h_w)
     
      enddo!(i_tor)
    enddo !(ivar) 
#endif
  endif
enddo !(i)

return
end subroutine update_values
