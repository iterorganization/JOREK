subroutine finish_grid(node_list, element_list, newnode_list, newelement_list, n_grids, include_axis, include_xpoint, include_psi)
!------------------------------------------------------------------------------------------
! subroutine defines the new nodes and elements of the final grid 
!------------------------------------------------------------------------------------------

use constants
use tr_module 
use data_structure
use grid_xpoint_data
use phys_module, only: xcase, RZ_grid_inside_wall, force_central_node, fix_axis_nodes, R_geo, Z_geo, xpoint, n_pol, treat_axis
use mod_eqdsk_tools
use mod_interp, only: interp_RZ, interp
use mod_element_rtree

implicit none

! --- Routine parameters
type (type_node_list)       , intent(inout) :: node_list
type (type_element_list)    , intent(inout) :: element_list
type (type_node_list)       , intent(inout) :: newnode_list
type (type_element_list)    , intent(inout) :: newelement_list
integer                     , intent(in)    :: n_grids(12)
logical                     , intent(in)    :: include_axis, include_xpoint, include_psi


integer             :: i, j, j2, k, l, my_id, ifail, n_tmp
integer             :: i_elm1, i_vertex1, i_node1, i_node_save
integer             :: i_elm2, i_vertex2, i_node2
integer             :: i_elm_xpoint(2), i_elm_axis
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
integer             :: i_save, j_save, ii, jj, count
integer             :: nR_eqdsk, nZ_eqdsk, ier
real*8, allocatable :: R_eqdsk(:),Z_eqdsk(:),psi_eqdsk(:,:)
real*8              :: psi1, psi2, psi3, psi4, psi_tmp, psi
real*8              :: psi_left, psi_right
real*8              :: deriv_left, deriv_right
real*8              :: R_elm, Z_elm
real*8              :: RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss
real*8              :: ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss
real*8              :: PSg1,dPSg1_dr,dPSg1_ds,dPSg1_drs,dPSg1_drr,dPSg1_dss
real*8              :: psi_xpoint(2), R_xpoint(2), Z_xpoint(2), s_xpoint(2), t_xpoint(2)
real*8              :: psi_axis, R_axis, Z_axis, s_axis, t_axis
real*8              :: R1, Z1, s_out, t_out, R_out, Z_out, RZ_jac, dRZ_jac_dR, dRZ_jac_dZ, PSI_R, PSI_Z, PSI_RR, PSI_ZZ, PSI_RZ
real*8              :: R0,Z0, RP,ZP, dR0, dZ0, dRP, dZP, size_0, size_p, denom
character*4         :: label
logical             :: normal_eqdsk, normal_eqdsk_wall
logical, parameter  :: plot_grid = .false.


write(*,*) '*****************************************'
write(*,*) '* X-point grid : Finalise grid          *'
write(*,*) '*****************************************'



n_flux = n_grids(1)
n_tht  = n_grids(2)
if (n_tht .eq. 0) n_tht = n_pol


!-------------------------------------------------------------------------------------------!
!------------------- Adjust size of elements to get better match ---------------------------!
!-------------------------------------------------------------------------------------------!
write(*,*) '                 Definition of elements size '

! --- renormalise vectors (in case we are extending from grid_polar_bezier)
do i=1,newnode_list%n_nodes
  do k=2,3
    size_0 = sqrt( newnode_list%node(i)%X(1,k,1)**2 + newnode_list%node(i)%X(1,k,2)**2 )
    if (size_0 .lt. 1.d-16) cycle ! axis
    newnode_list%node(i)%X(1,k,1) = newnode_list%node(i)%X(1,k,1) / size_0
    newnode_list%node(i)%X(1,k,2) = newnode_list%node(i)%X(1,k,2) / size_0
  enddo
enddo

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
enddo




!-------------------------------------------------------------------------------------------!
!--------------------------- Fill in the values into the new grid --------------------------!
!-------------------------------------------------------------------------------------------!
write(*,*) '                 Fill in psi-values '
psi    = 0.d0 ; PSI_R  = 0.d0 ; PSI_Z  = 0.d0
PSI_RR = 0.d0 ; PSI_ZZ = 0.d0 ; PSI_RZ = 0.d0
if (include_psi) then
  ier = 0
  if (RZ_grid_inside_wall) then
    call get_eqdsk_style(normal_eqdsk, normal_eqdsk_wall, ier)
    if (ier .ne. 0) then
      write(*,*)'Problem reading eqdsk file',eqdsk_filename
      write(*,*)'Aborting...'
      stop
    endif
    call get_eqdsk_dimensions(normal_eqdsk, nR_eqdsk, nZ_eqdsk, n_wall, ier)
    if (ier .eq. 0) then
      allocate(R_eqdsk(nR_eqdsk), Z_eqdsk(nZ_eqdsk), psi_eqdsk(nR_eqdsk,nZ_eqdsk))
      call get_data_from_eqdsk(normal_eqdsk, normal_eqdsk_wall, nR_eqdsk, nZ_eqdsk, R_eqdsk, Z_eqdsk, psi_eqdsk, n_wall, R_wall(1:n_wall), Z_wall(1:n_wall), ier)
    endif
  endif
  
  do i=1,newnode_list%n_nodes
  
    R1 = newnode_list%node(i)%x(1,1,1)
    Z1 = newnode_list%node(i)%x(1,1,2)
  
    ifail = 99
    if (n_flux .ne. 0) call find_RZ(node_list,element_list,R1,Z1,R_out,Z_out,ielm_out,s_out,t_out,ifail)
  
    if (ifail .ne. 0) then
      if ( (RZ_grid_inside_wall) .and. (ier .eq. 0) ) then
        ! --- psi values from eqdsk
        call interpolate_psi_from_eqdsk_grid(nR_eqdsk, nZ_eqdsk, R_eqdsk, Z_eqdsk, psi_eqdsk, R1, Z1, psi, psi_R, psi_Z)
      else
        write(*,'(A,2f15.4)')'Warning! did not find node one previous grid!',R1,Z1
        write(*,*)'Unable to extract psi information, the grid might be flawed.'
      endif
    else
      call interp_RZ(node_list,element_list,ielm_out,s_out,t_out, &
                     RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss, &
                     ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)
      call interp(node_list,element_list,ielm_out,1,1,s_out,t_out,PSg1,dPSg1_dr,dPSg1_ds,dPSg1_drs,dPSg1_drr,dPSg1_dss)
      RZ_jac  = dRRg1_dr * dZZg1_ds - dRRg1_ds * dZZg1_dr
      dRZ_jac_dR = (dRRg1_drr* dZZg1_ds**2 - dZZg1_drr*dRRg1_ds*dZZg1_ds - 2.d0*dRRg1_drs*dZZg1_dr*dZZg1_ds   &
                  + dZZg1_drs*(dRRg1_dr*dZZg1_ds + dRRg1_ds*dZZg1_dr)                                         &
                  + dRRg1_dss* dZZg1_dr**2 - dZZg1_dss*dRRg1_dr*dZZg1_dr) / RZ_jac
      dRZ_jac_dZ = (dZZg1_dss* dRRg1_dr**2 - dRRg1_dss*dZZg1_dr*dRRg1_dr - 2.d0*dZZg1_drs*dRRg1_ds*dRRg1_dr   &
                  + dRRg1_drs*(dZZg1_ds*dRRg1_dr + dZZg1_dr*dRRg1_ds)                                         &
                  + dZZg1_drr* dRRg1_ds**2 - dRRg1_drr*dZZg1_ds*dRRg1_ds) / RZ_jac
  
      psi    = PSg1
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
    endif
  
    newnode_list%node(i)%values(1,1,1) = psi
    newnode_list%node(i)%values(1,2,1) = PSI_R * newnode_list%node(i)%x(1,2,1) + PSI_Z * newnode_list%node(i)%x(1,2,2)
    newnode_list%node(i)%values(1,3,1) = PSI_R * newnode_list%node(i)%x(1,3,1) + PSI_Z * newnode_list%node(i)%x(1,3,2)
    newnode_list%node(i)%values(1,4,1) = PSI_RR * newnode_list%node(i)%x(1,2,1) * newnode_list%node(i)%x(1,3,1) &
                                       + PSI_RZ * newnode_list%node(i)%x(1,2,1) * newnode_list%node(i)%x(1,3,2) &
                                       + PSI_RZ * newnode_list%node(i)%x(1,3,1) * newnode_list%node(i)%x(1,2,2) &
                                       + PSI_ZZ * newnode_list%node(i)%x(1,2,2) * newnode_list%node(i)%x(1,3,2) &
                                       + PSI_R  * newnode_list%node(i)%x(1,4,1)                               &
                                       + PSI_Z  * newnode_list%node(i)%x(1,4,2)
  
    !if (newnode_list%node(i)%boundary .eq. 2) newnode_list%node(i)%values(1,3,1) = 0.d0 ! this is ok only if bnd 2 is aligned to surface!
  
  enddo
else
  ! --- This is very dirty, but if you don't know the grid, it's difficult to guess...
  do i=1,newnode_list%n_nodes
    if (node_list%node(i)%boundary .eq. 0) then
      newnode_list%node(i)%values(1,1,1) = -1.0 + sqrt( (newnode_list%node(i)%x(1,1,1)-R_geo)**2 + (newnode_list%node(i)%x(1,1,2)-Z_geo)**2 )
    else
      newnode_list%node(i)%values(1,1,1) = 0.0
    endif
    newnode_list%node(i)%values(1,2,1) = 0.d0
    newnode_list%node(i)%values(1,3,1) = 0.d0
    newnode_list%node(i)%values(1,4,1) = 0.d0
  enddo
endif

!-------------------------------------------------------------------------------------------!
!--------------------------- Fill in the values into the new grid --------------------------!
!-------------------------------------------------------------------------------------------!
write(*,*) '                 Copy new grid into old one '

!-------------------------------- Empty Xpoints
if (include_xpoint) then
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
endif

!-------------------------------- Empty Axis
if (include_axis) then
  if (include_xpoint) then
    if (xcase .ne. DOUBLE_NULL) then
      do j=5,4+n_tht-1
        newnode_list%node(j)%values(1,2:n_degrees,1) = 0.d0
      enddo
    else
      do j=9,8+n_tht-2
        newnode_list%node(j)%values(1,2:n_degrees,1) = 0.d0
      enddo
    endif
  else
    do j=2,n_tht
      newnode_list%node(j)%values(1,2:n_degrees,1) = 0.d0
    enddo
  endif
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

! --- This is the old way
!element_list%n_elements = newelement_list%n_elements
!element_list%element(1:element_list%n_elements) = newelement_list%element(1:element_list%n_elements)
!node_list%n_nodes = newnode_list%n_nodes
!node_list%node(1:node_list%n_nodes) = newnode_list%node(1:node_list%n_nodes)

! --- Make sure we copy only elements that have nodes
element_list%n_elements = 0
do i_elm1 = 1,newelement_list%n_elements
  if (newelement_list%element(i_elm1)%vertex(1) .eq. 0) cycle
  if (newelement_list%element(i_elm1)%vertex(2) .eq. 0) cycle
  if (newelement_list%element(i_elm1)%vertex(3) .eq. 0) cycle
  if (newelement_list%element(i_elm1)%vertex(4) .eq. 0) cycle
  element_list%n_elements = element_list%n_elements + 1
  element_list%element(element_list%n_elements) = newelement_list%element(i_elm1)
enddo

! --- Now, we define only the nodes that belong to elements! (this gets rid of potential orphan nodes, which the matrix doesn't like, obviously...)
node_list%n_nodes = 0
if (include_xpoint) then
  node_list%n_nodes = 4+n_tht-1
  node_list%node(1:node_list%n_nodes) = newnode_list%node(1:node_list%n_nodes)
  if (xcase .eq. DOUBLE_NULL) then
    node_list%n_nodes = 8+n_tht-2
    node_list%node(1:node_list%n_nodes) = newnode_list%node(1:node_list%n_nodes)
  endif
else
  node_list%n_nodes = n_tht
  node_list%node(1:node_list%n_nodes) = newnode_list%node(1:node_list%n_nodes)
endif
do i_elm1 = 1,element_list%n_elements
  do i_vertex1 = 1,n_vertex_max
    i_node1 = newelement_list%element(i_elm1)%vertex(i_vertex1)
    if ( (i_node1 .le. 8+n_tht-2) .and. (xcase .eq. DOUBLE_NULL) .and. (include_xpoint) ) cycle
    if ( (i_node1 .le. 4+n_tht-1) .and. (xcase .ne. DOUBLE_NULL) .and. (include_xpoint) ) cycle
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
  enddo
enddo
rtree_initialized = .false.
call populate_element_rtree(node_list, element_list)

!-------------------------------------------------------------------------------------------!
!------------------------------ Define nodes index in the matrix ---------------------------!
!-------------------------------------------------------------------------------------------!
write(*,*) '                 Definition of nodes index '

! --- Note: it's very important that we do this after copying the nodes and after eliminating the orphan nodes!

!-------------------------------- Combine multiple nodes at axis and Xpoints
index = 0
do i=1,node_list%n_nodes

  node_list%node(i)%axis_node = .false.
  node_list%node(i)%axis_dof  = 0  
  if (include_axis) then
    if (include_xpoint) then
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
    else
      if (i .le. n_tht) then
         node_list%node(i)%axis_node = .true.
         node_list%node(i)%axis_dof  = 3
      endif   
    endif
  endif


  do k=1,n_degrees

    index = index + 1
    node_list%node(i)%index(k) = index

    ! Remove all but one node at axis
    if (force_central_node .and. include_axis) then
      if (include_xpoint) then
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
      else
        if ((i .gt. 1) .and. (i .le. n_tht) .and. (k.eq.1)) then
          node_list%node(i)%index(k) = node_list%node(1)%index(1)
          index = index - 1
        endif
      endif
    endif
    
    ! Share 4 degrees of freedom for all nodes on the grid axis and flag the axis nodes.    
    if (treat_axis .and. include_axis) then
      if (include_xpoint) then
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
      else
        if ((i .gt. 1) .and. (i .le. n_tht) .and. (k.le.n_order+1)) then
          node_list%node(i)%index(k) = node_list%node(1)%index(k)
          index = index - 1
        endif
      endif
    endif
    
    ! Remove all but one node at first Xpoint
    if (include_xpoint) then
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
    endif
  
  enddo  
  node_list%node(i)%constrained = .false.
enddo

if (fix_axis_nodes .and. include_axis) then
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



!----temporary, needs to be completed, neighbour and boundary information
call update_neighbours_basic(element_list,node_list)
call update_boundary_types(element_list,node_list, include_xpoint)
call update_boundary_types_final(element_list,node_list)

my_id = 0 !Now we want the output...
call find_axis(my_id,node_list,element_list,psi_axis,R_axis,Z_axis,i_elm_axis,s_axis,t_axis,ifail)
if (xpoint) call find_xpoint(my_id,node_list,element_list,psi_xpoint,R_xpoint,Z_xpoint,i_elm_xpoint,s_xpoint,t_xpoint,xcase,ifail)


!----------------------------------- Print a python file that plots a cross with the 4 nodes of each element
if (plot_grid) then
  n_tmp = newelement_list%n_elements
  open(101,file='plot_finished_grid.py')
    write(101,'(A)')         '#!/usr/bin/env python'
    write(101,'(A)')         'import numpy as N'
    write(101,'(A)')         'import pylab'
    write(101,'(A)')         'def main():'
    write(101,'(A,i6,A)')    ' r = N.zeros(',4*n_tmp,')'
    write(101,'(A,i6,A)')    ' z = N.zeros(',4*n_tmp,')'
    do j=1,n_tmp
      do i=1,2
        index = newelement_list%element(j)%vertex(i)
        write(101,'(A,i6,A,f15.4)') ' r[',4*(j-1)+2*i-2,'] = ',newnode_list%node(index)%x(1,1,1)
        write(101,'(A,i6,A,f15.4)') ' z[',4*(j-1)+2*i-2,'] = ',newnode_list%node(index)%x(1,1,2)
        index = newelement_list%element(j)%vertex(i+2)
        write(101,'(A,i6,A,f15.4)') ' r[',4*(j-1)+2*i-1,'] = ',newnode_list%node(index)%x(1,1,1)
        write(101,'(A,i6,A,f15.4)') ' z[',4*(j-1)+2*i-1,'] = ',newnode_list%node(index)%x(1,1,2)
      enddo
    enddo
    write(101,'(A,i6,A)')    ' for i in range (0,',n_tmp*2,'):'
    write(101,'(A)')         '  if (i%2 == 0):'
    write(101,'(A)')         '   pylab.plot(r[2*i:2*i+2],z[2*i:2*i+2], "r")'
    write(101,'(A)')         '  else:'
    write(101,'(A)')         '   pylab.plot(r[2*i:2*i+2],z[2*i:2*i+2], "g")'
    do j=1,n_tmp
      do i=1,4
        index = newelement_list%element(j)%vertex(i)
        write(101,'(A,i6,A,f15.4)') ' r[',4*(j-1)+i-1,'] = ',newnode_list%node(index)%x(1,1,1)
        write(101,'(A,i6,A,f15.4)') ' z[',4*(j-1)+i-1,'] = ',newnode_list%node(index)%x(1,1,2)
      enddo
    enddo
    write(101,'(A,i6,A)')    ' for i in range (0,',n_tmp,'):'
    write(101,'(A)')         '  pylab.plot(r[4*i:4*i+4],z[4*i:4*i+4], "b")'
    write(101,'(A)')         ' pylab.axis("equal")'
    write(101,'(A)')         ' pylab.show()'
    write(101,'(A)')         ' '
    write(101,'(A)')         'main()'
  close(101)
endif

if (RZ_grid_inside_wall .and. include_psi) deallocate(R_eqdsk, Z_eqdsk, psi_eqdsk)



return
end subroutine finish_grid





