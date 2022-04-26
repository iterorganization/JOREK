subroutine define_final_grid(node_list, element_list, flux_list, &
                             xcase, n_grids, stpts, nwpts)
!------------------------------------------------------------------------------------------
! subroutine defines the new nodes and elements of the final grid 
!------------------------------------------------------------------------------------------

use constants
use tr_module 
use data_structure
use grid_xpoint_data
use mod_interp
use phys_module, only: write_ps, force_central_node, SDN_threshold, fix_axis_nodes, treat_axis
use equil_info

implicit none

! --- Routine parameters
type (type_surface_list)    , intent(in)    :: flux_list
type (type_node_list)       , intent(inout) :: node_list
type (type_element_list)    , intent(inout) :: element_list
type (type_strategic_points), intent(in)    :: stpts
type (type_new_points)      , intent(in)    :: nwpts
integer,                      intent(in)    :: n_grids(12), xcase


! --- local variables
type (type_node_list),    pointer :: newnode_list
type (type_element_list), pointer :: newelement_list

integer             :: i, j, j2, k, l, my_id, ifail, n_tmp
integer             :: i_elm1, i_vertex1, i_node1, i_node_save
integer             :: i_elm2, i_vertex2, i_node2
integer             :: n_loop, n_loop2, n_start_connect
integer             :: n_psi, n_tht_mid, n_tht_mid2
integer             :: n_flux, n_tht,   n_open,   n_outer,   n_inner    
integer             :: n_private,   n_up_priv,   n_leg,   n_up_leg
integer             :: index
integer             :: n_start_open, n_start_outer, n_start_inner
integer             :: n_start_private, n_start_up_priv
integer             :: n_xpoint_1, n_xpoint_2, n_xpoint_3, n_jump
integer             :: iv, ivp, node_iv, node_ivp, ielm_out
real*8, allocatable :: xp(:),yp(:)
real*8              :: RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss
real*8              :: ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss
real*8              :: PSg1,dPSg1_dr,dPSg1_ds,dPSg1_drs,dPSg1_drr,dPSg1_dss
real*8              :: R1, Z1, s_out, t_out, R_out, Z_out, RZ_jac, dRZ_jac_dR, dRZ_jac_dZ, PSI_R, PSI_Z, PSI_RR, PSI_ZZ, PSI_RZ
real*8              :: R0,Z0, RP,ZP, dR0, dZ0, dRP, dZP, size_0, size_p, denom
character*4         :: label
logical, parameter  :: plot_grid = .true.


write(*,*) '*****************************************'
write(*,*) '* X-point grid : Define the final grid  *'
write(*,*) '*****************************************'

n_flux     = n_grids(1); n_tht     = n_grids(2)
n_open     = n_grids(3); n_outer   = n_grids(4); n_inner = n_grids(5)
n_private  = n_grids(6); n_up_priv = n_grids(7)
n_leg      = n_grids(8); n_up_leg  = n_grids(9)
n_tht_mid  = n_grids(10)
n_tht_mid2 = n_tht-n_tht_mid

!------------------------------------------------------------------------------------------------------------------------!
!************************************************************************************************************************!
!************************************************************************************************************************!
!***************************************** First part: define the new nodes  ********************************************!
!************************************************************************************************************************!
!************************************************************************************************************************!
!------------------------------------------------------------------------------------------------------------------------!


!-------------------------------- Allocate data structures for new nodes and initialize them
allocate(newnode_list)
call tr_register_mem(sizeof(newnode_list),"newnode_list")
newnode_list%n_nodes = 0
newnode_list%n_dof   = 0
index = 0
do i = 1, n_nodes_max
  newnode_list%node(i)%x           = 0.d0
  newnode_list%node(i)%values      = 0.d0
  newnode_list%node(i)%deltas      = 0.d0
  newnode_list%node(i)%index       = 0
  newnode_list%node(i)%boundary    = 0
  newnode_list%node(i)%parents     = 0
  newnode_list%node(i)%parent_elem = 0
  newnode_list%node(i)%ref_lambda  = 0.d0
  newnode_list%node(i)%ref_mu      = 0.d0
  newnode_list%node(i)%constrained = .false.
end do




!-------------------------------------------------------------------------------------------!
!-------------------------------------- The Xpoints ----------------------------------------!
!-------------------------------------------------------------------------------------------!

! THIS ADDS FOUR NODES AT EACH XPOINTS, PLEASE SEE create_x_node FOR MORE DETAILS
if (xcase .eq. LOWER_XPOINT) then
  call create_x_node(node_list, element_list, newnode_list, nwpts, stpts, &
                     LOWER_XPOINT, ES%R_axis, ES%Z_axis, ES%R_xpoint, ES%Z_xpoint, ES%i_elm_xpoint, ES%s_xpoint, ES%t_xpoint)
endif 
if (xcase .eq. UPPER_XPOINT) then
  call create_x_node(node_list, element_list, newnode_list, nwpts, stpts, &
                     UPPER_XPOINT, ES%R_axis, ES%Z_axis, ES%R_xpoint, ES%Z_xpoint, ES%i_elm_xpoint, ES%s_xpoint, ES%t_xpoint)
endif 
if ( (xcase .eq. DOUBLE_NULL) .and. ( (ES%active_xpoint .eq. LOWER_XPOINT) .or. (ES%active_xpoint .eq. SYMMETRIC_XPOINT) ) ) then ! Put lower Xpoint first
  call create_x_node(node_list, element_list, newnode_list, nwpts, stpts, &
                     LOWER_XPOINT, ES%R_axis, ES%Z_axis, ES%R_xpoint, ES%Z_xpoint, ES%i_elm_xpoint, ES%s_xpoint, ES%t_xpoint)
  call create_x_node(node_list, element_list, newnode_list, nwpts, stpts, &
                     UPPER_XPOINT, ES%R_axis, ES%Z_axis, ES%R_xpoint, ES%Z_xpoint, ES%i_elm_xpoint, ES%s_xpoint, ES%t_xpoint)
endif 
if ( (xcase .eq. DOUBLE_NULL) .and. ( ES%active_xpoint .eq. UPPER_XPOINT )                                          ) then ! Put upper Xpoint first
  call create_x_node(node_list, element_list, newnode_list, nwpts, stpts, &
                     UPPER_XPOINT, ES%R_axis, ES%Z_axis, ES%R_xpoint, ES%Z_xpoint, ES%i_elm_xpoint, ES%s_xpoint, ES%t_xpoint)
  call create_x_node(node_list, element_list, newnode_list, nwpts, stpts, &
                     LOWER_XPOINT, ES%R_axis, ES%Z_axis, ES%R_xpoint, ES%Z_xpoint, ES%i_elm_xpoint, ES%s_xpoint, ES%t_xpoint)
endif 
index = newnode_list%n_nodes




!-------------------------------------------------------------------------------------------!
!-------------------------------- The closed field lines -----------------------------------!
!-------------------------------------------------------------------------------------------!

do i=1,n_flux                 
  n_loop = n_tht-1
  if (xcase .eq. DOUBLE_NULL) n_loop = n_loop-1  ! For double-null, n_tht_mid and n_tht_mid+1 are the same lines
  do k=1, n_loop

    j = k
    if ((xcase .eq. DOUBLE_NULL) .and. (j .gt. n_tht_mid)) j = j+1  ! For double-null, n_tht_mid and n_tht_mid+1 are the same lines
    
    index = index + 1
    call create_new_node(node_list, element_list, newnode_list, index, i, j, nwpts)

    if (i .eq. 1) then   !------------------------------------ magnetic axis : special case
      newnode_list%node(index)%x(1,3,:) = 0.d0
      newnode_list%node(index)%x(1,4,:) = 0.d0
    endif

  enddo
enddo
newnode_list%n_nodes = index




!-------------------------------------------------------------------------------------------!
!--------------- The open part (or sandwich part, in case of double-null) ------------------!
!-------------------------------------------------------------------------------------------!

       !- (the one around the main plasma)
       !- Starting from right strike point, all around to left strike point (if lower Xpoint is the main one)
       !- The other way round if upper Xpoint is the main one
       !- We add 4 points for each Xpoint, either on the left or on the right 
       !- (see routine create_x_node for more info)
       
n_start_open = newnode_list%n_nodes + 1
if (xcase .eq. LOWER_XPOINT) n_loop = 2*(n_leg-1)    + n_tht 
if (xcase .eq. UPPER_XPOINT) n_loop = 2*(n_up_leg-1) + n_tht
if ( (xcase .eq. DOUBLE_NULL) .and. ( (ES%active_xpoint .eq. LOWER_XPOINT) .or. (ES%active_xpoint .eq. SYMMETRIC_XPOINT) ) ) n_loop = 2*(n_leg-1)    + n_tht - 1 
if ( (xcase .eq. DOUBLE_NULL) .and. (  ES%active_xpoint .eq. UPPER_XPOINT)                                          ) n_loop = 2*(n_up_leg-1) + n_tht - 1
n_loop2 = n_flux+n_open+1
if (xcase .eq. DOUBLE_NULL) n_loop2 = n_flux+n_open               
n_tmp = n_tht
if (xcase .eq. DOUBLE_NULL) n_tmp = n_tht-1  ! n_tht_mid and n_tht_mid+1 are the same for double null          
if (ES%active_xpoint .ne. SYMMETRIC_XPOINT) then  ! ignore if symmetric double-null
  
  do i=n_flux+1,n_loop2
    
    do l=1,n_loop                 
    
      !-------------------------------- CASE 1 : Lower Xpoint is the main one
      if ( (ES%active_xpoint .eq. LOWER_XPOINT) .or. (ES%active_xpoint .eq. SYMMETRIC_XPOINT) ) then
        
        !---- First part : right lower strike line
        if (l .le. n_leg-1) then 
          j = n_tht + n_leg + l
          index = index + 1
          call create_new_node(node_list, element_list, newnode_list, index, i, j, nwpts)
        endif    
        
        !---- Second part : around the main plasma
        if ( (l .gt. n_leg-1) .and. (l .le. (n_leg-1)+n_tmp) ) then      
          j = l-(n_leg-1)
          if ( ((j .gt. 1) .and. (j .lt. n_tmp)) .or. (i .ne. n_flux+1) ) then ! Don't put the Xpoint twice
            if ( (xcase .eq. DOUBLE_NULL) .and. (j .gt. n_tht_mid) ) j = j+1  ! n_tht_mid and n_tht_mid+1 are the same for double null
            index = index + 1
            call create_new_node(node_list, element_list, newnode_list, index, i, j, nwpts)    
          endif    
        endif    
        
        !---- Third part : left lower strike line
        if (l .gt. (n_leg-1)+n_tmp) then      
          j = l-(n_leg-1)-n_tmp
          j = n_tht + n_leg - j
          index = index + 1
          call create_new_node(node_list, element_list, newnode_list, index, i, j, nwpts)    
        endif
        
      endif


      !-------------------------------- CASE 2 : Upper Xpoint is the main one
      if ( ES%active_xpoint .eq. UPPER_XPOINT ) then
        
        !---- First part : left upper strike line
        if (l .le. n_up_leg-1) then      
          j = n_tht + 2*n_leg + l
          index = index + 1
          call create_new_node(node_list, element_list, newnode_list, index, i, j, nwpts)
        endif    
        
        !---- Second part : around the main plasma
        if ( (l .gt. n_up_leg-1) .and. (l .le. (n_up_leg-1)+n_tmp) ) then    
          j = l-(n_up_leg-1)
          if ( ((j .gt. 1) .and. (j .lt. n_tmp)) .or. (i .ne. n_flux+1) ) then ! Don't put the Xpoint twice
            if ( (xcase .eq. DOUBLE_NULL) .and. (j .gt. n_tht_mid) ) j = j+1  ! n_tht_mid and n_tht_mid+1 are the same for double null
            index = index + 1
            call create_new_node(node_list, element_list, newnode_list, index, i, j, nwpts) 
          endif   
        endif    
        
        !---- Third part : right upper strike line
        if (l .gt. (n_up_leg-1)+n_tmp) then      
          j = l-(n_up_leg-1)-n_tmp
          j = n_tht + 2*n_leg + 2*n_up_leg - j
          index = index + 1
          call create_new_node(node_list, element_list, newnode_list, index, i, j, nwpts)    
        endif    
        
      endif
      
      if ( (xcase .ne. DOUBLE_NULL) .and. (i .eq. n_loop2) ) newnode_list%node(index)%boundary = newnode_list%node(index)%boundary + 2
      if ( (l .eq. 1) .or. (l .eq. n_loop)                 ) newnode_list%node(index)%boundary = newnode_list%node(index)%boundary + 1
      
    enddo
  enddo
endif
newnode_list%n_nodes = index





!-------------------------------------------------------------------------------------------!
!------------------------ The outer part, in case of double-null ---------------------------!
!-------------------------------------------------------------------------------------------!
       
n_start_outer = newnode_list%n_nodes + 1
if (xcase .eq. DOUBLE_NULL) then
  if ( (ES%active_xpoint .eq. LOWER_XPOINT) .or. (ES%active_xpoint .eq. SYMMETRIC_XPOINT) ) n_tmp = n_tht_mid 
  if (  ES%active_xpoint .eq. UPPER_XPOINT                                         ) n_tmp = n_tht_mid2
  n_loop = (n_leg-1) + (n_up_leg-1) + n_tmp
  do i=n_flux+n_open+1,n_flux+n_open+n_outer+1                
    do l=1,n_loop                 
    
      !---- First part : right lower strike line
      if (l .le. n_leg-1) then 
        j = n_tht + n_leg + l
        index = index + 1
        call create_new_node(node_list, element_list, newnode_list, index, i, j, nwpts)
      endif    

      !---- Second part : along the outer side of the separatrix
      if ( (l .gt. n_leg-1) .and. (l .le. (n_leg-1)+n_tmp) ) then        
        j = l-(n_leg-1)
        if ( (ES%active_xpoint .eq. LOWER_XPOINT) .or. (ES%active_xpoint .eq. SYMMETRIC_XPOINT) ) then 
          if ( (j .lt. n_tmp) .or. (i .ne. n_flux+n_open+1) ) then ! Don't put the Xpoint twice
            index = index + 1
            call create_new_node(node_list, element_list, newnode_list, index, i, j, nwpts)
          endif    
        else 
          if ( (j .gt. 1) .or. (i .ne. n_flux+n_open+1) ) then ! Don't put the Xpoint twice
            j = j+n_tht_mid
            index = index + 1
            call create_new_node(node_list, element_list, newnode_list, index, i, j, nwpts)    
          endif  
        endif  
      endif    
      
      !---- Third part : right upper strike line
      if (l .gt. (n_leg-1)+n_tmp) then      
        j = l-(n_leg-1)-n_tmp
        j = n_tht + 2*n_leg + 2*n_up_leg - j
        index = index + 1
        call create_new_node(node_list, element_list, newnode_list, index, i, j, nwpts)    
      endif    
          
      if ( (i .eq. n_flux+n_open+n_outer+1) ) newnode_list%node(index)%boundary = newnode_list%node(index)%boundary + 2
      if ( (l .eq. 1) .or. (l .eq. n_loop) )  newnode_list%node(index)%boundary = newnode_list%node(index)%boundary + 1
        
    enddo
  enddo
  newnode_list%n_nodes = index
endif





!-------------------------------------------------------------------------------------------!
!------------------------ The inner part, in case of double-null ---------------------------!
!-------------------------------------------------------------------------------------------!
       
n_start_inner = newnode_list%n_nodes + 1
if (xcase .eq. DOUBLE_NULL) then
  if ( (ES%active_xpoint .eq. LOWER_XPOINT) .or. (ES%active_xpoint .eq. SYMMETRIC_XPOINT) ) n_tmp = n_tht_mid2
  if (  ES%active_xpoint .eq. UPPER_XPOINT                                         ) n_tmp = n_tht_mid
  n_loop = (n_leg-1) + (n_up_leg-1) + n_tmp
  do k=n_flux+n_open+n_outer+1,n_flux+n_open+n_outer+n_inner+1                
    do l=1,n_loop                 
    
      i = k
      if (i .eq. n_flux+n_open+n_outer+1) i = n_flux+n_open+1  !Doing the second separatrix first
      
      !---- First part : left upper strike line
      if (l .le. n_up_leg-1) then 
        j = n_tht + 2*n_leg + l
        index = index + 1
        call create_new_node(node_list, element_list, newnode_list, index, i, j, nwpts)
      endif    

      !---- Second part : along the outer side of the separatrix
      if ( (l .gt. n_up_leg-1) .and. (l .le. (n_up_leg-1)+n_tmp) ) then          
        j = l-(n_up_leg-1)
        if ( (ES%active_xpoint .eq. LOWER_XPOINT) .or. (ES%active_xpoint .eq. SYMMETRIC_XPOINT) ) then
          if ( (j .gt. 1) .or. (i .ne. n_flux+n_open+1) ) then ! Don't put the Xpoint twice
            j = j+n_tht_mid
            index = index + 1
            call create_new_node(node_list, element_list, newnode_list, index, i, j, nwpts)    
          endif  
        else
          if ( (j .lt. n_tmp) .or. (i .ne. n_flux+n_open+1) ) then ! Don't put the Xpoint twice
            index = index + 1
            call create_new_node(node_list, element_list, newnode_list, index, i, j, nwpts)    
          endif  
        endif  
      endif          
      
      !---- Third part : left lower strike line
      if (l .gt. (n_up_leg-1)+n_tmp) then      
        j = l-(n_up_leg-1)-n_tmp
        j = n_tht + n_leg - j
        index = index + 1
        call create_new_node(node_list, element_list, newnode_list, index, i, j, nwpts)    
      endif    
          
      if ( (i .eq. n_flux+n_open+n_outer+n_inner+1) ) newnode_list%node(index)%boundary = newnode_list%node(index)%boundary + 2
      if ( (l .eq. 1) .or. (l .eq. n_loop) )          newnode_list%node(index)%boundary = newnode_list%node(index)%boundary + 1
        
    enddo
  enddo
  newnode_list%n_nodes = index
endif





!-------------------------------------------------------------------------------------------!
!------------------------------- The private parts -----------------------------------------!
!-------------------------------------------------------------------------------------------!
       
!-------------------------------- The lower private
n_start_private = newnode_list%n_nodes + 1
if(xcase .ne. UPPER_XPOINT) then
  index = newnode_list%n_nodes
  do i=n_flux+n_open+n_outer+n_inner+1+1,n_flux+n_open+n_outer+n_inner+n_private+1              
    do k=1,2  ! Two loops for the two legs
      do l=1, n_leg                         

        if (k .eq. 1) j = n_tht + n_leg + l
        if (k .eq. 2) j = n_tht + n_leg - l + 1

        if ( (k .ne. 2) .or. (l .ne. 1) ) then   ! Do the middle line only once

          index = index + 1
          call create_new_node(node_list, element_list, newnode_list, index, i, j, nwpts)    

          if (i .eq. n_flux+n_open+n_outer+n_inner+n_private+1)    newnode_list%node(index)%boundary = newnode_list%node(index)%boundary + 2
          if ( (l .eq. 1) .or. ((l .eq. n_leg) .and. (k .eq. 2)) ) newnode_list%node(index)%boundary = newnode_list%node(index)%boundary + 1

        endif

      enddo
    enddo
  enddo
endif
newnode_list%n_nodes = index

!-------------------------------- The upper private
n_start_up_priv = newnode_list%n_nodes + 1
if(xcase .ne. LOWER_XPOINT) then
  index = newnode_list%n_nodes
  do i=n_flux+n_open+n_outer+n_inner+n_private+1+1,n_flux+n_open+n_outer+n_inner+n_private+n_up_priv+1
    do k=1,2  ! Two loops for the two legs
      do l=1, n_up_leg                         

        if (k .eq. 1) j = n_tht + 2*n_leg + l
        if (k .eq. 2) j = n_tht + 2*n_leg + 2*n_up_leg - l + 1

        if ( (k .ne. 2) .or. (l .ne. 1) ) then   ! Do the middle line only once

          index = index + 1
          call create_new_node(node_list, element_list, newnode_list, index, i, j, nwpts)    

          if (i .eq. n_flux+n_open+n_outer+n_inner+n_private+n_up_priv+1)  newnode_list%node(index)%boundary = newnode_list%node(index)%boundary + 2
          if ( (l .eq. 1) .or. ((l .eq. n_up_leg) .and. (k .eq. 2)) )      newnode_list%node(index)%boundary = newnode_list%node(index)%boundary + 1

        endif

      enddo
    enddo
  enddo
endif
newnode_list%n_nodes = index




!-------------------------------------------------------------------------------------------!
!--------------------------------- Plot the nodes ------------------------------------------!
!-------------------------------------------------------------------------------------------!
write(*,*) '                 Definition of nodes complete : '
write(*,*) '                     number of nodes = ',newnode_list%n_nodes

!-------------------------------- Plot all the nodes
if (plot_grid .and. write_ps) then
  call nframe(21,11,1,1.0,5.0,-1.8,2.2,' ',1,'R',1,'Z',1)
  call plot_flux_surfaces(node_list,element_list,flux_list,.false.,1,.true.,xcase)
  call tr_allocate(xp,1,index,"xp")
  call tr_allocate(yp,1,index,"yp")
  do i=1,newnode_list%n_nodes
    xp(i) = newnode_list%node(i)%x(1,1,1)
    yp(i) = newnode_list%node(i)%x(1,1,2)
  enddo
  call lplot(1,1,421,xp,yp,-newnode_list%n_nodes,1,'R',1,'Z',1,'nodes',5)
  call tr_deallocate(xp,"xp")
  call tr_deallocate(yp,"yp")
endif

!-------------------------------- Plot the boundary nodes only
if (plot_grid .and. write_ps) then
  call nframe(21,11,1,1.0,5.0,-1.8,2.2,' ',1,'R',1,'Z',1)
  call plot_flux_surfaces(node_list,element_list,flux_list,.false.,1,.true.,xcase)
  index=0
  do i=1,newnode_list%n_nodes
    if(newnode_list%node(i)%boundary .ne. 0) then
      index=index+1
    endif
  enddo
  call tr_allocate(xp,1,index,"xp")
  call tr_allocate(yp,1,index,"yp")
  index=0
  do i=1,newnode_list%n_nodes
    if(newnode_list%node(i)%boundary .ne. 0) then
      index=index+1
      xp(index) = newnode_list%node(i)%x(1,1,1)
      yp(index) = newnode_list%node(i)%x(1,1,2)
    endif
  enddo
  call lplot(1,1,421,xp,yp,-index,1,'R',1,'Z',1,'nodes',5)
  call tr_deallocate(xp,"xp")
  call tr_deallocate(yp,"yp")
endif

!-------------------------------- Plot the divertor nodes only
if (plot_grid .and. write_ps) then
  call nframe(21,11,1,1.0,5.0,-1.8,2.2,' ',1,'R',1,'Z',1)
  call plot_flux_surfaces(node_list,element_list,flux_list,.false.,1,.true.,xcase)
  index=0
  do i=1,newnode_list%n_nodes
    if(newnode_list%node(i)%boundary .eq. 1) then
      index=index+1
    endif
  enddo
  call tr_allocate(xp,1,index,"xp")
  call tr_allocate(yp,1,index,"yp")
  index=0
  do i=1,newnode_list%n_nodes
    if(newnode_list%node(i)%boundary .eq. 1) then
      index=index+1
      xp(index) = newnode_list%node(i)%x(1,1,1)
      yp(index) = newnode_list%node(i)%x(1,1,2)
    endif
  enddo
  call lplot(1,1,421,xp,yp,-index,1,'R',1,'Z',1,'nodes',5)
  call tr_deallocate(xp,"xp")
  call tr_deallocate(yp,"yp")
endif

!-------------------------------- Plot the open flux surface nodes only
if (plot_grid .and. write_ps ) then
  call nframe(21,11,1,1.0,5.0,-1.8,2.2,' ',1,'R',1,'Z',1)
  call plot_flux_surfaces(node_list,element_list,flux_list,.false.,1,.true.,xcase)
  index=0
  do i=1,newnode_list%n_nodes
    if(newnode_list%node(i)%boundary .eq. 2) then
      index=index+1
    endif
  enddo
  call tr_allocate(xp,1,index,"xp")
  call tr_allocate(yp,1,index,"yp")
  index=0
  do i=1,newnode_list%n_nodes
    if(newnode_list%node(i)%boundary .eq. 2) then
      index=index+1
      xp(index) = newnode_list%node(i)%x(1,1,1)
      yp(index) = newnode_list%node(i)%x(1,1,2)
    endif
  enddo
  call lplot(1,1,421,xp,yp,-index,1,'R',1,'Z',1,'nodes',5)
  call tr_deallocate(xp,"xp")
  call tr_deallocate(yp,"yp")
endif

!-------------------------------- Plot the corner nodes only
if (plot_grid .and. write_ps ) then
  call nframe(21,11,1,1.0,5.0,-1.8,2.2,' ',1,'R',1,'Z',1)
  call plot_flux_surfaces(node_list,element_list,flux_list,.false.,1,.true.,xcase)
  index=0
  do i=1,newnode_list%n_nodes
    if(newnode_list%node(i)%boundary .eq. 3) then
      index=index+1
    endif
  enddo
  call tr_allocate(xp,1,index,"xp")
  call tr_allocate(yp,1,index,"yp")
  index=0
  do i=1,newnode_list%n_nodes
    if(newnode_list%node(i)%boundary .eq. 3) then
      index=index+1
      xp(index) = newnode_list%node(i)%x(1,1,1)
      yp(index) = newnode_list%node(i)%x(1,1,2)
    endif
  enddo
  call lplot(1,1,421,xp,yp,-index,1,'R',1,'Z',1,'nodes',5)
  call tr_deallocate(xp,"xp")
  call tr_deallocate(yp,"yp")
endif

!-------------------------------- Verify that only the Xpoint and the magnetic axis appear more than once
!do i=n_tht-1,newnode_list%n_nodes
!  do j=i+1,newnode_list%n_nodes
!    if((newnode_list%node(i)%x(1,1,1) .eq. newnode_list%node(j)%x(1,1,1))      &
!      .and. (newnode_list%node(i)%x(1,1,2) .eq. newnode_list%node(j)%x(1,1,2)) ) then
!      write(*,*)'Found ij',i,j
!      write(*,*)'RZ',newnode_list%node(i)%x(1,1,1),newnode_list%node(i)%x(1,1,2)
!    endif
!  enddo
!enddo
















!------------------------------------------------------------------------------------------------------------------------!
!************************************************************************************************************************!
!************************************************************************************************************************!
!*************************************** Second part: define the new elements  ******************************************!
!************************************************************************************************************************!
!************************************************************************************************************************!
!------------------------------------------------------------------------------------------------------------------------!


!-------------------------------------------------------------------------------------------!
!------------------------------ Define the elements ----------------------------------------!
!-------------------------------------------------------------------------------------------!

       
!-------------------------------- Allocate data structures for new elements and initialize them
index = 0
allocate(newelement_list)
call tr_register_mem(sizeof(newelement_list),"newelement_list")
newelement_list%n_elements = 0
do i = 1, n_elements_max
  newelement_list%element(i)%vertex       = 0
  newelement_list%element(i)%neighbours   = 0
  newelement_list%element(i)%size         = 0.d0
  newelement_list%element(i)%father       = 0
  newelement_list%element(i)%n_sons       = 0
  newelement_list%element(i)%n_gen        = 0
  newelement_list%element(i)%sons         = 0
  newelement_list%element(i)%contain_node = 0
  newelement_list%element(i)%nref         = 0
end do


!-------------------------------- The closed region
n_tmp = 4
if (xcase .eq. DOUBLE_NULL) n_tmp = 8
do i=1,n_flux
  n_loop = n_tht-1
  if (xcase .eq. DOUBLE_NULL) n_loop = n_tht-2
  do l=1, n_loop

    j = l
    index = index + 1
    newelement_list%element(index)%size = 1.d0

    newelement_list%element(index)%vertex(1) = n_tmp + (i-1)*n_loop + j
    newelement_list%element(index)%vertex(2) = n_tmp + (i  )*n_loop + j
    newelement_list%element(index)%vertex(3) = n_tmp + (i  )*n_loop + j + 1
    newelement_list%element(index)%vertex(4) = n_tmp + (i-1)*n_loop + j + 1

    if (j .eq. n_loop) then
      newelement_list%element(index)%vertex(4)  = n_tmp + (i-1)*n_loop + 1
      newelement_list%element(index)%vertex(3)  = n_tmp + (i  )*n_loop + 1
    endif

    ! Connect with open (or sandwich) region (or outer and inner in case of symmetric double-null)
    if (i .eq. n_flux) then
      if ( (ES%active_xpoint .eq. LOWER_XPOINT) .or.(ES%active_xpoint .eq. SYMMETRIC_XPOINT) ) then
        n_xpoint_1 = n_start_open + n_leg-1 -1
      else 
        n_xpoint_1 = n_start_open + n_up_leg-1 -1
      endif
      if (  ES%active_xpoint .eq. SYMMETRIC_XPOINT  ) n_xpoint_1 = n_start_outer + n_leg-1
      if ( (ES%active_xpoint .eq. SYMMETRIC_XPOINT  ) .and. (l .gt. n_tht_mid-1) ) then
        n_xpoint_1 = n_start_inner + n_up_leg-1 -1
        j = j - n_tht_mid+1
      endif
      newelement_list%element(index)%vertex(2) = n_xpoint_1 + j - 1 
      newelement_list%element(index)%vertex(3) = n_xpoint_1 + j
      if (l .eq. 1) then ! Special case for element starting at Xpoint
        newelement_list%element(index)%vertex(2) = 2 ! We need to start at NODE2 of the lower Xpoint (or NODE6=NODE2 of the upper Xpoint)
      endif
      if (l .eq. n_loop) then ! Special case for element arriving at Xpoint
        newelement_list%element(index)%vertex(3) = 3 ! We need to finish at NODE3 of the lower Xpoint (or NODE7=NODE3 of the upper Xpoint)
      endif
      if ( (ES%active_xpoint .eq. SYMMETRIC_XPOINT) .and. (l .eq. n_tht_mid-1) ) then ! Special case for element passing through 2nd Xpoint (for symmetric only)
        newelement_list%element(index)%vertex(3) = 7 ! We need to arrive at NODE7 of the upper Xpoint
      endif
      if ( (ES%active_xpoint .eq. SYMMETRIC_XPOINT) .and. (l .eq. n_tht_mid)   ) then ! Special case for element passing through 2nd Xpoint (for symmetric only)
        newelement_list%element(index)%vertex(2) = 6 ! We need to leave from NODE6 of the upper Xpoint
      endif
    endif  
      
  enddo
enddo
newelement_list%n_elements = index

if ( newelement_list%n_elements > n_elements_max ) then
  write(*,*) 'ERROR in define_final_grid: hard-coded parameter n_elements_max is too small'
  stop
end if

!-------------------------------- The open (or sandwich) region (between the two separatrices)
if (xcase .eq. LOWER_XPOINT) then
  n_xpoint_1 = n_leg-1                   ! First time through Xpoint
  n_xpoint_2 = n_leg-1 + n_tht-1         ! Second time through Xpoint
  n_loop     = 2*(n_leg-1) + n_tht-1
endif 
if ((xcase .eq. DOUBLE_NULL ) .and. ( (ES%active_xpoint .eq. LOWER_XPOINT) .or. (ES%active_xpoint .eq. SYMMETRIC_XPOINT) ) ) then
  n_xpoint_1 = n_leg-1                   ! First time through first Xpoint
  n_xpoint_2 = n_leg-1 + n_tht-2         ! Second time through first Xpoint
  n_xpoint_3 = n_leg-1 + n_tht_mid-1     ! Going through second Xpoint
  n_loop     = 2*(n_leg-1) + n_tht-2
endif 
if ((xcase .eq. DOUBLE_NULL ) .and. ( ES%active_xpoint .eq. UPPER_XPOINT )) then
  n_xpoint_1 = n_up_leg-1                ! First time through first Xpoint
  n_xpoint_2 = n_up_leg-1 + n_tht-2      ! Second time through first Xpoint
  n_xpoint_3 = n_up_leg-1 + n_tht_mid-1  ! Going through second Xpoint
  n_loop     = 2*(n_up_leg-1) + n_tht-2
endif 
if (xcase .eq. UPPER_XPOINT) then
  n_xpoint_1 = n_up_leg-1                ! First time through first Xpoint
  n_xpoint_2 = n_up_leg-1 + n_tht-1      ! Second time through first Xpoint
  n_loop     = 2*(n_up_leg-1) + n_tht-1
endif 
if ( ES%active_xpoint .ne. SYMMETRIC_XPOINT ) then ! ignore if symmetric double-null  
  do i=1,n_open
    do l=1, n_loop

      j = l
      j2 = l
      index = index + 1
      newelement_list%element(index)%size = 1.d0

      ! Note the -2 stands for the two Xpoints at the separatrix
      ! Note that for n_loop elements, there are n_loop+1 nodes
      newelement_list%element(index)%vertex(1) = n_start_open + (i-1)*(n_loop+1)-2 + j - 1
      newelement_list%element(index)%vertex(4) = n_start_open + (i-1)*(n_loop+1)-2 + j
      newelement_list%element(index)%vertex(2) = n_start_open + (i  )*(n_loop+1)-2 + j - 1
      newelement_list%element(index)%vertex(3) = n_start_open + (i  )*(n_loop+1)-2 + j
      
      ! Connect with closed region
      if (i .eq. 1) then
        if (l .gt. n_xpoint_1) j = j-1 ! First time through Xpoint
        if (l .gt. n_xpoint_2) j = j-1 ! Second time through Xpoint
        newelement_list%element(index)%vertex(1) = n_start_open + j - 1
        newelement_list%element(index)%vertex(4) = n_start_open + j
        if ( (l .eq. n_xpoint_1) .or. (l .eq. n_xpoint_2) ) then ! Special case for element arriving at Xpoint
          newelement_list%element(index)%vertex(4) = 4 ! We need to arrive at NODE8=NODE4 of the lower Xpoint (or NODE4 of the upper Xpoint)
        endif
        if ( (l .eq. n_xpoint_1+1) .or. (l .eq. n_xpoint_2+1) ) then ! Special case for element leaving Xpoint
          newelement_list%element(index)%vertex(1) = 1 ! We need to leave at NODE1 of the lower Xpoint (or NODE5=NODE1 of the upper Xpoint)
        endif
      endif

      ! Connect with outer and inner regions
      if ( (i .eq. n_open) .and. (xcase .eq. DOUBLE_NULL) ) then
        if ( (ES%active_xpoint .eq. LOWER_XPOINT) .or. (ES%active_xpoint .eq. SYMMETRIC_XPOINT) ) then
          if (l .le. n_xpoint_3) then
            newelement_list%element(index)%vertex(2) = n_start_outer + j2 - 1
            newelement_list%element(index)%vertex(3) = n_start_outer + j2
          else
            j2 = j2-n_xpoint_3
            newelement_list%element(index)%vertex(2) = n_start_inner + n_up_leg-1 - 1 + j2 - 1
            newelement_list%element(index)%vertex(3) = n_start_inner + n_up_leg-1 - 1 + j2
          endif
        else
          if (l .le. n_xpoint_3) then
            newelement_list%element(index)%vertex(2) = n_start_inner + j2 - 1
            newelement_list%element(index)%vertex(3) = n_start_inner + j2
          else
            j2 = j2-n_xpoint_3
            newelement_list%element(index)%vertex(2) = n_start_outer + n_leg-1 - 1 + j2 - 1
            newelement_list%element(index)%vertex(3) = n_start_outer + n_leg-1 - 1 + j2
          endif
        endif
        if (l .eq. n_xpoint_3) then ! Special case for element arriving at second Xpoint
          newelement_list%element(index)%vertex(3) = 7 ! We need to arrive at NODE7 of the upper Xpoint (or NODE3=NODE7 of the lower Xpoint)
        endif
        if (l .eq. n_xpoint_3+1) then ! Special case for element leaving second Xpoint
          newelement_list%element(index)%vertex(2) = 6 ! We need to leave at NODE6 of the upper Xpoint (or NODE2=NODE6 of the lower Xpoint)
        endif
      endif

    enddo
  enddo
endif
newelement_list%n_elements = index

if ( newelement_list%n_elements > n_elements_max ) then
  write(*,*) 'ERROR in define_final_grid: hard-coded parameter n_elements_max is too small'
  stop
end if

!-------------------------------- The outer region
if (xcase .eq. DOUBLE_NULL) then
  if ( (ES%active_xpoint .eq. LOWER_XPOINT) .or. (ES%active_xpoint .eq. SYMMETRIC_XPOINT) ) then
    n_xpoint_1 = (n_leg-1)+(n_tht_mid-1)
    n_loop     = (n_leg-1)+(n_up_leg-1)+(n_tht_mid-1)
  else
    n_xpoint_1 = (n_leg-1)
    n_loop     = (n_leg-1)+(n_up_leg-1)+(n_tht_mid2-1)
  endif 
  if (ES%active_xpoint .eq. SYMMETRIC_XPOINT) then
    n_xpoint_1 = (n_leg-1)+(n_tht_mid-1)
    n_xpoint_2 = (n_leg-1)
    n_loop     = (n_leg-1)+(n_up_leg-1)+(n_tht_mid-1)
  endif 
  do i=1,n_outer
    do l=1, n_loop

      j = l
      index = index + 1
      newelement_list%element(index)%size = 1.d0

      newelement_list%element(index)%vertex(1) = n_start_outer + (i-1)*(n_loop+1)-1 + j - 1
      newelement_list%element(index)%vertex(4) = n_start_outer + (i-1)*(n_loop+1)-1 + j
      newelement_list%element(index)%vertex(2) = n_start_outer + (i  )*(n_loop+1)-1 + j - 1
      newelement_list%element(index)%vertex(3) = n_start_outer + (i  )*(n_loop+1)-1 + j
      
      ! Connect with open region (or closed and private for symmetric double-null)
      if (i .eq. 1) then
        if (j .gt. n_xpoint_1) j = j-1
        newelement_list%element(index)%vertex(1) = n_start_outer + j - 1
        newelement_list%element(index)%vertex(4) = n_start_outer + j
        if (l .eq. n_xpoint_1)   then ! Special case for element arriving at second Xpoint
          newelement_list%element(index)%vertex(4) = 8 ! We need to arrive at NODE8=NODE4 of the second Xpoint
        endif
        if (l .eq. n_xpoint_1+1) then ! Special case for element leaving second Xpoint
          newelement_list%element(index)%vertex(1) = 5 ! We need to leave at NODE1 of the second Xpoint
        endif
        if ( (l .eq. n_xpoint_2)   .and. (ES%active_xpoint .eq. SYMMETRIC_XPOINT) ) then ! Special case for element arriving at first Xpoint of symmetric double-null
          newelement_list%element(index)%vertex(4) = 4 ! We need to arrive at NODE8=NODE4 of the first Xpoint
        endif
        if ( (l .eq. n_xpoint_2+1) .and. (ES%active_xpoint .eq. SYMMETRIC_XPOINT) ) then ! Special case for element leaving first Xpoint of symmetric double-null
          newelement_list%element(index)%vertex(1) = 1 ! We need to leave at NODE1 of the first Xpoint
        endif
      endif

    enddo
  enddo
endif
newelement_list%n_elements = index

if ( newelement_list%n_elements > n_elements_max ) then
  write(*,*) 'ERROR in define_final_grid: hard-coded parameter n_elements_max is too small'
  stop
end if

!-------------------------------- The inner region
if (xcase .eq. DOUBLE_NULL) then
  if ( (ES%active_xpoint .eq. LOWER_XPOINT) .or. (ES%active_xpoint .eq. SYMMETRIC_XPOINT) ) then
    n_xpoint_1 = (n_up_leg-1)
    n_loop     = (n_up_leg-1)+(n_leg-1)+(n_tht_mid2-1)
  else
    n_xpoint_1 = (n_up_leg-1)+(n_tht_mid-1)
    n_loop     = (n_up_leg-1)+(n_leg-1)+(n_tht_mid-1)
  endif 
  if (ES%active_xpoint .eq. SYMMETRIC_XPOINT) then
    n_xpoint_1 = (n_up_leg-1)
    n_xpoint_2 = (n_up_leg-1)+(n_tht_mid2-1)
    n_loop     = (n_up_leg-1)+(n_leg-1)+(n_tht_mid2-1)
  endif 
  do i=1,n_inner
    do l=1, n_loop

      j = l
      index = index + 1
      newelement_list%element(index)%size = 1.d0

      newelement_list%element(index)%vertex(1) = n_start_inner + (i-1)*(n_loop+1)-1 + j - 1
      newelement_list%element(index)%vertex(4) = n_start_inner + (i-1)*(n_loop+1)-1 + j
      newelement_list%element(index)%vertex(2) = n_start_inner + (i  )*(n_loop+1)-1 + j - 1
      newelement_list%element(index)%vertex(3) = n_start_inner + (i  )*(n_loop+1)-1 + j
      
      ! Connect with open region
      if (i .eq. 1) then
        if (j .gt. n_xpoint_1) j = j-1
        newelement_list%element(index)%vertex(1) = n_start_inner + j - 1
        newelement_list%element(index)%vertex(4) = n_start_inner + j
        if (l .eq. n_xpoint_1) then ! Special case for element arriving at second Xpoint
          newelement_list%element(index)%vertex(4) = 8 ! We need to arrive at NODE4 of the second Xpoint
        endif
        if (l .eq. n_xpoint_1+1) then! Special case for element leaving second Xpoint
          newelement_list%element(index)%vertex(1) = 5 ! We need to leave at NODE5=NODE1 of the second Xpoint
        endif
        if ( (l .eq. n_xpoint_2)   .and. (ES%active_xpoint .eq. SYMMETRIC_XPOINT) ) then ! Special case for element arriving at first Xpoint of symmetric double-null
          newelement_list%element(index)%vertex(4) = 4 ! We need to arrive at NODE4 of the first Xpoint
        endif
        if ( (l .eq. n_xpoint_2+1) .and. (ES%active_xpoint .eq. SYMMETRIC_XPOINT) ) then ! Special case for element leaving first Xpoint of symmetric double-null
          newelement_list%element(index)%vertex(1) = 1 ! We need to leave at NODE5=NODE1 of the first Xpoint
        endif
      endif

    enddo
  enddo
endif
newelement_list%n_elements = index

if ( newelement_list%n_elements > n_elements_max ) then
  write(*,*) 'ERROR in define_final_grid: hard-coded parameter n_elements_max is too small'
  stop
end if

!-------------------------------- The lower private region
if (xcase .ne. UPPER_XPOINT) then
  n_loop          = 2*(n_leg-1)
  n_xpoint_1      = n_leg-1
  n_start_connect = n_start_open
  if (  xcase .eq. LOWER_XPOINT                                                                                       ) n_jump = n_tht-3
  if ( (xcase .eq. DOUBLE_NULL) .and. ( (ES%active_xpoint .eq. LOWER_XPOINT) .or. (ES%active_xpoint .eq. SYMMETRIC_XPOINT) ) ) n_jump = n_tht-4
  if ( (xcase .eq. DOUBLE_NULL) .and. ( ES%active_xpoint .eq. UPPER_XPOINT )                                          ) then
    n_start_connect  =   n_start_outer
    n_jump           = - n_start_outer - (n_leg-1) + n_start_inner + (n_up_leg-1) + (n_tht_mid-1) -1
  endif
  if (ES%active_xpoint .eq. SYMMETRIC_XPOINT) then ! Symmetric double-node
    n_start_connect  =   n_start_outer
    n_jump           = - n_start_outer - (n_leg-1) + n_start_inner + (n_up_leg-1) + (n_tht_mid2-1) -1
  endif
  do i=1,n_private
    do l=1, n_loop

      j = l
      index = index + 1
      newelement_list%element(index)%size = 1.d0

      newelement_list%element(index)%vertex(1) = n_start_private + (i-1)*(n_loop+1) + j - 1
      newelement_list%element(index)%vertex(4) = n_start_private + (i-1)*(n_loop+1) + j
      if (i .gt. 1) then
        newelement_list%element(index)%vertex(2) = n_start_private + (i-2)*(n_loop+1) + j - 1
        newelement_list%element(index)%vertex(3) = n_start_private + (i-2)*(n_loop+1) + j
      endif
      
      ! Connect with open (or outer and inner) region
      if (i .eq. 1) then
        if (l .gt. n_xpoint_1) j = j+n_jump
        newelement_list%element(index)%vertex(2) = n_start_connect + j - 1
        newelement_list%element(index)%vertex(3) = n_start_connect + j
        if (l .eq. n_xpoint_1) then ! Special case for element arriving at lower Xpoint
          newelement_list%element(index)%vertex(3) = 3 ! We need to arrive at NODE7=NODE3 of the Xpoint
          if ( (xcase .eq. DOUBLE_NULL) .and. (ES%active_xpoint .eq. UPPER_XPOINT) ) newelement_list%element(index)%vertex(3) = 7 
        endif
        if (l .eq. n_xpoint_1+1) then! Special case for element leaving lower Xpoint
          newelement_list%element(index)%vertex(2) = 2 ! We need to leave at NODE6=NODE2 of the Xpoint
          if ( (xcase .eq. DOUBLE_NULL) .and. (ES%active_xpoint .eq. UPPER_XPOINT) ) newelement_list%element(index)%vertex(2) = 6 
        endif
      endif

    enddo
  enddo
endif
newelement_list%n_elements = index

if ( newelement_list%n_elements > n_elements_max ) then
  write(*,*) 'ERROR in define_final_grid: hard-coded parameter n_elements_max is too small'
  stop
end if

!-------------------------------- The upper private region
if (xcase .ne. LOWER_XPOINT) then
  n_loop          = 2*(n_up_leg-1)
  n_xpoint_1      = n_up_leg-1
  n_start_connect = n_start_open
  if (  xcase .eq. UPPER_XPOINT                                                                                       ) n_jump = n_tht-3
  if ( (xcase .eq. DOUBLE_NULL) .and. (  ES%active_xpoint .eq. UPPER_XPOINT)                                          ) n_jump = n_tht-4
  if ( (xcase .eq. DOUBLE_NULL) .and. ( (ES%active_xpoint .eq. LOWER_XPOINT) .or. (ES%active_xpoint .eq. SYMMETRIC_XPOINT) ) ) then ! Warning - note that we include the symmetric double-node here
    n_start_connect  =   n_start_inner
    n_jump           = - n_start_inner - (n_up_leg-1) + n_start_outer + (n_leg-1) + (n_tht_mid-1) -1
  endif
  do i=1,n_up_priv
    do l=1, n_loop

      j = l
      index = index + 1
      newelement_list%element(index)%size = 1.d0

      newelement_list%element(index)%vertex(1) = n_start_up_priv + (i-1)*(n_loop+1) + j - 1
      newelement_list%element(index)%vertex(4) = n_start_up_priv + (i-1)*(n_loop+1) + j
      if (i .gt. 1) then
        newelement_list%element(index)%vertex(2) = n_start_up_priv + (i-2)*(n_loop+1) + j - 1
        newelement_list%element(index)%vertex(3) = n_start_up_priv + (i-2)*(n_loop+1) + j
      endif
      
      ! Connect with open (or outer and inner) region
      if (i .eq. 1) then
        if (l .gt. n_xpoint_1) j = j+n_jump
        newelement_list%element(index)%vertex(2) = n_start_connect + j - 1
        newelement_list%element(index)%vertex(3) = n_start_connect + j
        if (l .eq. n_xpoint_1) then ! Special case for element arriving at upper Xpoint
          newelement_list%element(index)%vertex(3) = 3 ! We need to arrive at NODE3 of the Xpoint
          if ( (xcase .eq. DOUBLE_NULL) .and. ( (ES%active_xpoint .eq. LOWER_XPOINT) .or. (ES%active_xpoint .eq. SYMMETRIC_XPOINT) ) ) newelement_list%element(index)%vertex(3) = 7 
        endif
        if (l .eq. n_xpoint_1+1) then! Special case for element leaving upper Xpoint
          newelement_list%element(index)%vertex(2) = 2 ! We need to leave at NODE2 of the Xpoint
          if ( (xcase .eq. DOUBLE_NULL) .and. ( (ES%active_xpoint .eq. LOWER_XPOINT) .or. (ES%active_xpoint .eq. SYMMETRIC_XPOINT) ) ) newelement_list%element(index)%vertex(2) = 6 
        endif
      endif

    enddo
  enddo
endif
newelement_list%n_elements = index

if ( newelement_list%n_elements > n_elements_max ) then
  write(*,*) 'ERROR in define_final_grid: hard-coded parameter n_elements_max is too small'
  stop
end if


!-------------------------------------------------------------------------------------------!
!------------------------------- Plot the new elements -------------------------------------!
!-------------------------------------------------------------------------------------------!
write(*,*) '                 Definition of elements complete : '
write(*,*) '                     number of elements = ',newelement_list%n_elements

!-------------------------------- Plot the elements' nodes
if (plot_grid .and. write_ps ) then
  call nframe(21,11,1,1.0,5.0,-1.8,2.2,' ',1,'R',1,'Z',1)
  call plot_flux_surfaces(node_list,element_list,flux_list,.false.,1,.true.,xcase)
  k = 4*index
  call tr_allocate(xp,1,k,"xp")
  call tr_allocate(yp,1,k,"yp")
  do j=1,newelement_list%n_elements
    do i=1,4
      k = newelement_list%element(j)%vertex(i)
      xp(4*(j-1)+i) = newnode_list%node(k)%x(1,1,1)
      yp(4*(j-1)+i) = newnode_list%node(k)%x(1,1,2)
    enddo
  enddo
  k = 4*index
  call lplot(1,1,421,xp,yp,-k,1,'R',1,'Z',1,'nodes',5)
  call tr_deallocate(xp,"xp")
  call tr_deallocate(yp,"yp")
endif

!----------------------------------- Print a python file that plots a cross with the 4 nodes of each element
if (plot_grid .and. write_ps) then
  n_loop = newelement_list%n_elements
  open(101,file='plot_elements.py')
    write(101,'(A)')        '#!/usr/bin/env python'
    write(101,'(A)')        'import numpy as N'
    write(101,'(A)')        'import pylab'
    write(101,'(A)')        'def main():'
    write(101,'(A,i6,A)')    ' r = N.zeros(',4*n_loop,')'
    write(101,'(A,i6,A)')    ' z = N.zeros(',4*n_loop,')'
    do j=1,n_loop
      do i=1,2
        index = newelement_list%element(j)%vertex(i)
        write(101,'(A,i6,A,f15.4)') ' r[',4*(j-1)+2*i-2,'] = ',newnode_list%node(index)%x(1,1,1)
        write(101,'(A,i6,A,f15.4)') ' z[',4*(j-1)+2*i-2,'] = ',newnode_list%node(index)%x(1,1,2)
        index = newelement_list%element(j)%vertex(i+2)
        write(101,'(A,i6,A,f15.4)') ' r[',4*(j-1)+2*i-1,'] = ',newnode_list%node(index)%x(1,1,1)
        write(101,'(A,i6,A,f15.4)') ' z[',4*(j-1)+2*i-1,'] = ',newnode_list%node(index)%x(1,1,2)
      enddo
    enddo
    write(101,'(A,i6,A)')    ' for i in range (0,',n_loop*2,'):'
    write(101,'(A)')        '  pylab.plot(r[2*i:2*i+2],z[2*i:2*i+2], "r")'
    do j=1,n_loop
      do i=1,4
        index = newelement_list%element(j)%vertex(i)
        write(101,'(A,i6,A,f15.4)') ' r[',4*(j-1)+i-1,'] = ',newnode_list%node(index)%x(1,1,1)
        write(101,'(A,i6,A,f15.4)') ' z[',4*(j-1)+i-1,'] = ',newnode_list%node(index)%x(1,1,2)
      enddo
    enddo
    write(101,'(A,i6,A)')    ' for i in range (0,',n_loop,'):'
    write(101,'(A)')        '  pylab.plot(r[4*i:4*i+4],z[4*i:4*i+4], "b")'
    write(101,'(A)')        ' pylab.axis("equal")'
    write(101,'(A)')        ' pylab.show()'
    write(101,'(A)')        ' '
    write(101,'(A)')        'main()'
  close(101)
endif




!-------------------------------------------------------------------------------------------!
!------------------- Adjust size of elements to get better match ---------------------------!
!-------------------------------------------------------------------------------------------!


do k=1, newelement_list%n_elements   ! fill in the size of the elements
  do iv = 1, 4                    ! over 4 sides of an element

    ivp = mod(iv,4)   + 1         ! vertex with index one higher
    node_iv  = newelement_list%element(k)%vertex(iv)
    node_ivp = newelement_list%element(k)%vertex(ivp) 

    if ((iv .eq. 1) .or. (iv .eq. 3)) then
      R0 = newnode_list%node(node_iv )%X(1,1,1)  ; dR0 = newnode_list%node(node_iv )%X(1,2,1)
      Z0 = newnode_list%node(node_iv )%X(1,1,2)  ; dZ0 = newnode_list%node(node_iv )%X(1,2,2)
      RP = newnode_list%node(node_ivp)%X(1,1,1)  ; dRP = newnode_list%node(node_ivp)%X(1,2,1)
      ZP = newnode_list%node(node_ivp)%X(1,1,2)  ; dZP = newnode_list%node(node_ivp)%X(1,2,2)
    else
      R0 = newnode_list%node(node_iv )%X(1,1,1)  ; dR0 = newnode_list%node(node_iv )%X(1,3,1)
      Z0 = newnode_list%node(node_iv )%X(1,1,2)  ; dZ0 = newnode_list%node(node_iv )%X(1,3,2)
      RP = newnode_list%node(node_ivp)%X(1,1,1)  ; dRP = newnode_list%node(node_ivp)%X(1,3,1)
      ZP = newnode_list%node(node_ivp)%X(1,1,2)  ; dZP = newnode_list%node(node_ivp)%X(1,3,2)
    endif

    size_0 = 1.d0
    size_p = 1.d0
    denom = ( dRP * dZ0 - dR0 * dZP)
    size_0 = sign(sqrt((R0-RP)**2 + (Z0-ZP)**2) /3.d0, dR0 * (RP-R0) + dZ0 * (ZP-Z0) )
    size_P = sign(sqrt((R0-RP)**2 + (Z0-ZP)**2) /3.d0, dRP * (R0-RP) + dZP * (Z0-ZP) )

    if ((R0-RP)**2 + (Z0-ZP)**2 .eq. 0.d0) then
      size_0 = 1.d0
      size_P = 1.d0
    endif

    if ((iv .eq. 1) .or. (iv .eq. 3)) then
      newelement_list%element(k)%size(iv,2)  = size_0
      newelement_list%element(k)%size(ivp,2) = size_p
    else
      newelement_list%element(k)%size(iv,3)  = size_0
      newelement_list%element(k)%size(ivp,3) = size_p
    endif

  enddo

  do iv=1,4
    newelement_list%element(k)%size(iv,1) = 1.d0
    newelement_list%element(k)%size(iv,4) = newelement_list%element(k)%size(iv,2) * newelement_list%element(k)%size(iv,3)
  enddo

  newelement_list%element(k)%father     = 0
  newelement_list%element(k)%n_sons     = 0
  element_list%element(Index)%sons(:)   = 0
enddo




!-------------------------------------------------------------------------------------------!
!--------------------------- Fill in the values into the new grid --------------------------!
!-------------------------------------------------------------------------------------------!


!-------------------------------- Fill in PSI-values at each node
do i=1,newnode_list%n_nodes

  R1 = newnode_list%node(i)%x(1,1,1)
  Z1 = newnode_list%node(i)%x(1,1,2)

  call find_RZ(node_list,element_list,R1,Z1,R_out,Z_out,ielm_out,s_out,t_out,ifail)

  if (ifail .ne. 0) then
    write(*,'(A,2f15.4)')'Warning! did not find node one previous grid!',R1,Z1
    write(*,*)'Unable to extract psi information, the grid might be flawed.'
  endif
  
  call interp_RZ(node_list,element_list,ielm_out,s_out,t_out, &
                 RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss, &
                 ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)

  call interp(node_list,element_list,ielm_out,1,1,s_out,t_out,PSg1,dPSg1_dr,dPSg1_ds,dPSg1_drs,dPSg1_drr,dPSg1_dss)

  RZ_jac     = dRRg1_dr * dZZg1_ds - dRRg1_ds * dZZg1_dr
  dRZ_jac_dR = (dRRg1_drr* dZZg1_ds**2 - dZZg1_drr*dRRg1_ds*dZZg1_ds - 2.d0*dRRg1_drs*dZZg1_dr*dZZg1_ds   &
              + dZZg1_drs*(dRRg1_dr*dZZg1_ds + dRRg1_ds*dZZg1_dr)                                         &
              + dRRg1_dss* dZZg1_dr**2 - dZZg1_dss*dRRg1_dr*dZZg1_dr) / RZ_jac
  dRZ_jac_dZ = (dZZg1_dss* dRRg1_dr**2 - dRRg1_dss*dZZg1_dr*dRRg1_dr - 2.d0*dZZg1_drs*dRRg1_ds*dRRg1_dr   &
              + dRRg1_drs*(dZZg1_ds*dRRg1_dr + dZZg1_dr*dRRg1_ds)                                         &
              + dZZg1_drr* dRRg1_ds**2 - dRRg1_drr*dZZg1_ds*dRRg1_ds) / RZ_jac

  PSI_R  = (   dZZg1_ds * dPSg1_dr - dZZg1_dr * dPSg1_ds ) / RZ_jac
  PSI_Z  = ( - dRRg1_ds * dPSg1_dr + dRRg1_dr * dPSg1_ds ) / RZ_jac
  PSI_RR = (dPSg1_drr * dZZg1_ds**2 - 2.d0*dPSg1_drs * dZZg1_dr*dZZg1_ds + dPSg1_dss * dZZg1_dr**2     &
           + dPSg1_dr * (dZZg1_drs*dZZg1_ds - dZZg1_dss*dZZg1_dr )                                     &
           + dPSg1_ds * (dZZg1_drs*dZZg1_dr - dZZg1_drr*dZZg1_ds ) )  / RZ_jac**2                      &
           - dRZ_jac_dR * (dPSg1_dr * dZZg1_ds - dPSg1_ds * dZZg1_dr) / RZ_jac**2
  PSI_ZZ = (dPSg1_drr * dRRg1_ds**2 - 2.d0*dPSg1_drs * dRRg1_dr*dRRg1_ds + dPSg1_dss * dRRg1_dr**2     &
           + dPSg1_dr * (dRRg1_drs*dRRg1_ds - dRRg1_dss*dRRg1_dr )                                     &
           + dPSg1_ds * (dRRg1_drs*dRRg1_dr - dRRg1_drr*dRRg1_ds ) )     / RZ_jac**2                   &
           - dRZ_jac_dZ * (- dPSg1_dr * dRRg1_ds + dPSg1_ds * dRRg1_dr ) / RZ_jac**2
  PSI_RZ = (- dPSg1_drr * dZZg1_ds*dRRg1_ds - dPSg1_dss * dRRg1_dr*dZZg1_dr                  &
           + dPSg1_drs * (dZZg1_dr*dRRg1_ds  + dZZg1_ds*dRRg1_dr  )                          &
           - dPSg1_dr  * (dRRg1_drs*dZZg1_ds - dRRg1_dss*dZZg1_dr )                          &
           - dPSg1_ds  * (dRRg1_drs*dZZg1_dr - dRRg1_drr*dZZg1_ds )  )     / RZ_jac**2       &
           - dRZ_jac_dR * (- dPSg1_dr * dRRg1_ds + dPSg1_ds * dRRg1_dr )   / RZ_jac**2

  newnode_list%node(i)%values(1,1,1) = PSg1
  newnode_list%node(i)%values(1,2,1) = PSI_R * newnode_list%node(i)%x(1,2,1) + PSI_Z * newnode_list%node(i)%x(1,2,2)
  newnode_list%node(i)%values(1,3,1) = PSI_R * newnode_list%node(i)%x(1,3,1) + PSI_Z * newnode_list%node(i)%x(1,3,2)
  newnode_list%node(i)%values(1,4,1) = PSI_RR * newnode_list%node(i)%x(1,2,1) * newnode_list%node(i)%x(1,3,1) &
                                     + PSI_RZ * newnode_list%node(i)%x(1,2,1) * newnode_list%node(i)%x(1,3,2) &
                                     + PSI_RZ * newnode_list%node(i)%x(1,3,1) * newnode_list%node(i)%x(1,2,2) &
                                     + PSI_ZZ * newnode_list%node(i)%x(1,2,2) * newnode_list%node(i)%x(1,3,2) &
                                     + PSI_R  * newnode_list%node(i)%x(1,4,1)                               &
                                     + PSI_Z  * newnode_list%node(i)%x(1,4,2)

  if (newnode_list%node(i)%boundary .eq. 2) newnode_list%node(i)%values(1,3,1) = 0.d0

enddo

!-------------------------------- Empty Xpoints
newnode_list%node(1)%values(1,2:n_degrees,1) = 0.d0
newnode_list%node(2)%values(1,2:n_degrees,1) = 0.d0
newnode_list%node(3)%values(1,2:n_degrees,1) = 0.d0
newnode_list%node(4)%values(1,2:n_degrees,1) = 0.d0
if (xcase .eq. DOUBLE_NULL) then
  newnode_list%node(5)%values(1,2:n_degrees,1) = 0.d0
  newnode_list%node(6)%values(1,2:n_degrees,1) = 0.d0
  newnode_list%node(7)%values(1,2:n_degrees,1) = 0.d0
  newnode_list%node(8)%values(1,2:n_degrees,1) = 0.d0
endif

!-------------------------------- Empty Axis
if (xcase .ne. DOUBLE_NULL) then
  do j=5,4+n_tht-1
    newnode_list%node(j)%values(1,2:n_degrees,1) = 0.d0
  enddo
else
  do j=9,8+n_tht-2
    newnode_list%node(j)%values(1,2:n_degrees,1) = 0.d0
  enddo
endif

!-------------------------------- Empty old nodes/elements
do i=1,node_list%n_nodes
  node_list%node(i)%x        = 0.d0
  node_list%node(i)%values   = 0.d0
  node_list%node(i)%index    = 0
  node_list%node(i)%boundary = 0
enddo
node_list%n_nodes = 0

do i=1,element_list%n_elements
  element_list%element(i)%vertex     = 0
  element_list%element(i)%size       = 0.d0
  element_list%element(i)%neighbours = 0
enddo

!---------------------------- copy new grid into nodes/elements
element_list%n_elements = newelement_list%n_elements
element_list%element(1:element_list%n_elements) = newelement_list%element(1:element_list%n_elements)

! --- This is the old way
!node_list%n_nodes = newnode_list%n_nodes
!node_list%node(1:node_list%n_nodes) = newnode_list%node(1:node_list%n_nodes)

! --- Now, we define only the nodes that belong to elements! (this gets rid of potential orphan nodes, which the matrix doesn't like, obviously...)
node_list%n_nodes = 4+n_tht-1
node_list%node(1:node_list%n_nodes) = newnode_list%node(1:node_list%n_nodes)
if (xcase .eq. DOUBLE_NULL) then
  node_list%n_nodes = 8+n_tht-2
  node_list%node(1:node_list%n_nodes) = newnode_list%node(1:node_list%n_nodes)
endif
do i_elm1 = 1,element_list%n_elements
  do i_vertex1 = 1,n_vertex_max
    i_node1 = newelement_list%element(i_elm1)%vertex(i_vertex1)
    if ( ((i_node1.gt.8+n_tht-2).and.(xcase .eq. DOUBLE_NULL)) .or. ((i_node1.gt.4+n_tht-1).and.(xcase .ne. DOUBLE_NULL)) ) then
      i_node_save = 0
      do i_elm2 = 1,i_elm1-1
        do i_vertex2 = 1,n_vertex_max
          i_node2 = newelement_list%element(i_elm2)%vertex(i_vertex2)
          if (i_node2 .eq. i_node1) then
            i_node_save = i_node1
            exit
          endif
        enddo
        if (i_node_save .ne. 0) exit
      enddo
      if (i_node_save .eq. 0) then
        node_list%n_nodes = node_list%n_nodes + 1
        node_list%node(node_list%n_nodes) = newnode_list%node(i_node1)
        newnode_list%node(i_node1)%boundary = node_list%n_nodes ! using "boundary" to save new node index, since newnode_list will be scrapped
        element_list%element(i_elm1)%vertex(i_vertex1) = node_list%n_nodes
      else
        element_list%element(i_elm1)%vertex(i_vertex1) = newnode_list%node(i_node_save)%boundary
      endif
    endif
  enddo
enddo

!-------------------------------- Combine multiple nodes at axis and Xpoints
! --- Note: it's very important that we do this after copying the nodes and after eliminating the orphan nodes!
index = 0
do i=1,newnode_list%n_nodes

  node_list%node(i)%axis_node = .false.
  node_list%node(i)%axis_dof  = 0
  if (xcase .ne. DOUBLE_NULL) then
    if ((i .ge. 5) .and. (i .le. 4+n_tht-1)) then
       node_list%node(i)%axis_node = .true.
       node_list%node(i)%axis_dof  = 2
    endif
  else
    if ((i .ge. 9) .and. (i .le. 8+n_tht-2)) then
       node_list%node(i)%axis_node = .true.
       node_list%node(i)%axis_dof  = 2
    endif
  endif

  do k=1,n_degrees

    index = index + 1
    node_list%node(i)%index(k) = index

    ! Remove all but one node at axis
    if (force_central_node) then
      if (xcase .ne. DOUBLE_NULL) then
        if ((i .gt. 5) .and. (i .le. 4+n_tht-1) .and. (k.eq.1)) then
          node_list%node(i)%index(k) = node_list%node(5)%index(1)
          index = index - 1
        endif
      else
        if ((i .gt. 9) .and. (i .le. 8+n_tht-2) .and. (k.eq.1)) then
          node_list%node(i)%index(k) = node_list%node(9)%index(1)
          index = index - 1
        endif
      endif
    endif
    
    ! Share 4 degrees of freedom for all nodes on the grid axis and flag the axis nodes. ONLY FOR C1-elements at the moment!
    if (treat_axis) then
      if (xcase .ne. DOUBLE_NULL) then
        if ((i .gt. 5) .and. (i .le. 4+n_tht-1) .and. (k.le.n_order+1)) then
          node_list%node(i)%index(k) = node_list%node(5)%index(k)
          index = index - 1
        endif
      else
        if ((i .gt. 9) .and. (i .le. 8+n_tht-2) .and. (k.le.n_order+1)) then
          node_list%node(i)%index(k) = node_list%node(9)%index(k)
          index = index - 1
        endif
      endif
    endif
    
    ! Remove all but one node at first Xpoint
    if ((i .eq. 2).and.(k.eq.1)) then
      node_list%node(i)%index(k) = node_list%node(1)%index(k)
      index = index - 1
    endif
    if ((i .eq. 2).and.(k.eq.3)) then
      node_list%node(i)%index(k) = node_list%node(1)%index(k)
      index = index - 1
    endif
    if ((i .eq. 3).and.(k.eq.1)) then
      node_list%node(i)%index(k) = node_list%node(1)%index(k)
      index = index - 1
    endif
    if ((i .eq. 3).and.(k.eq.2)) then
      node_list%node(i)%index(k) = node_list%node(2)%index(k)
      index = index - 1
    endif
    if ((i .eq. 4).and.(k.eq.1)) then
      node_list%node(i)%index(k) = node_list%node(1)%index(k)
      index = index - 1
    endif
    if ((i .eq. 4).and.(k.eq.2)) then
      node_list%node(i)%index(k) = node_list%node(1)%index(k)
      index = index - 1
    endif
    if ((i .eq. 4).and.(k.eq.3)) then
      node_list%node(i)%index(k) = node_list%node(3)%index(k)
      index = index - 1
    endif
  
    ! Remove all but one node at second Xpoint
    if (xcase .eq. DOUBLE_NULL) then
      if ((i .eq. 6).and.(k.eq.1)) then
        node_list%node(i)%index(k) = node_list%node(5)%index(k)
        index = index - 1
      endif
      if ((i .eq. 6).and.(k.eq.3)) then
        node_list%node(i)%index(k) = node_list%node(5)%index(k)
        index = index - 1
      endif
      if ((i .eq. 7).and.(k.eq.1)) then
        node_list%node(i)%index(k) = node_list%node(5)%index(k)
        index = index - 1
      endif
      if ((i .eq. 7).and.(k.eq.2)) then
        node_list%node(i)%index(k) = node_list%node(6)%index(k)
        index = index - 1
      endif
      if ((i .eq. 8).and.(k.eq.1)) then
        node_list%node(i)%index(k) = node_list%node(5)%index(k)
        index = index - 1
      endif
      if ((i .eq. 8).and.(k.eq.2)) then
        node_list%node(i)%index(k) = node_list%node(5)%index(k)
        index = index - 1
      endif
      if ((i .eq. 8).and.(k.eq.3)) then
        node_list%node(i)%index(k) = node_list%node(7)%index(k)
        index = index - 1
      endif
    endif
  
  enddo  
  node_list%node(i)%constrained = .false.
enddo

if (fix_axis_nodes) then
  do k=1, element_list%n_elements
    do iv=1,4
      j = element_list%element(k)%vertex(iv)
      if (node_list%node(j)%axis_node) then
        element_list%element(k)%size(iv,3) = 0.d0
        element_list%element(k)%size(iv,4) = 0.d0
      endif
    enddo
  enddo
endif




!----temporary, needs to be completed, neighbour information
do i=1, element_list%n_elements
  element_list%element(i)%father = 0
  element_list%element(i)%n_sons = 0
  element_list%element(i)%sons(:) = 0
enddo

call tr_unregister_mem(sizeof(newnode_list),"newnode_list")
deallocate(newnode_list)
call tr_unregister_mem(sizeof(newelement_list),"newelement_list")
deallocate(newelement_list)

my_id = 0 !Now we want the output...
call find_axis(my_id,node_list,element_list,ES%psi_axis,ES%R_axis,ES%Z_axis,ES%i_elm_axis,ES%s_axis,ES%t_axis,ifail)
call find_xpoint(my_id,node_list,element_list,ES%psi_xpoint,ES%R_xpoint,ES%Z_xpoint,ES%i_elm_xpoint,ES%s_xpoint,ES%t_xpoint,xcase,ifail)

return
end subroutine 
