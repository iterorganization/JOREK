!> module with grid conversion routines to create bi-quintic grid from bi-cubic grid
module mod_grid_conversions




contains










!> This routine just sets the n_order>=5 sizes with the basic rules
subroutine set_high_order_sizes(element_list)


use mod_parameters
use data_structure
use mod_node_indices


implicit none

! --- Routine variables
type(type_element_list), intent(inout) :: element_list   !< list of elements with element information

! --- Local variables
integer :: i_elm, i_vertex, k, l, node_index
integer :: node_indices( (n_order+1)/2, (n_order+1)/2 )
real*8  :: size_u, size_v

if (n_order .lt. 5) return

! --- calculate node_indices
call calculate_node_indices(node_indices)

do i_elm = 1,element_list%n_elements
  do i_vertex = 1,n_vertex_max
    size_u = element_list%element(i_elm)%size(i_vertex,2)
    size_v = element_list%element(i_elm)%size(i_vertex,3)
    ! --- C2 continuity
    element_list%element(i_elm)%size(i_vertex,5) = 1.d0
    element_list%element(i_elm)%size(i_vertex,6) = 1.d0
    element_list%element(i_elm)%size(i_vertex,7) = size_v
    element_list%element(i_elm)%size(i_vertex,8) = size_u
    element_list%element(i_elm)%size(i_vertex,9) = 1.d0
    ! --- C3 and above
    ! --- Loop over all indices
    do k = 1,(n_order+1)/2
      do l = 1,(n_order+1)/2
        node_index = node_indices(k,l)
        if (node_index .le. 9) cycle ! we want only derivatives >=3
        if ( (k .eq. 1) .or. (k .eq. 3) ) then
          element_list%element(i_elm)%size(i_vertex,node_index) = size_v
        elseif ( (l .eq. 1) .or. (l .eq. 3) ) then
          element_list%element(i_elm)%size(i_vertex,node_index) = size_u
        else
          element_list%element(i_elm)%size(i_vertex,node_index) = size_u*size_v
        endif
      enddo
    enddo
  enddo
enddo

return

end subroutine set_high_order_sizes











!> This routine just sets the sizes of _t derivatives of n_order>=5 to zero
subroutine set_high_order_sizes_on_axis(node_list,element_list)


use mod_parameters
use data_structure
use mod_node_indices
use phys_module, only: fix_axis_nodes


implicit none

! --- Routine variables
type(type_node_list),    intent(inout) :: node_list      !< list of nodes with grid information
type(type_element_list), intent(inout) :: element_list   !< list of elements with element information

! --- Local variables
integer :: i_elm, i_vertex, i_node, k, l, node_index
integer :: node_indices( (n_order+1)/2, (n_order+1)/2 )

if (.not. fix_axis_nodes) return
if (n_order .lt. 5) return

! --- calculate node_indices
call calculate_node_indices(node_indices)

do i_elm = 1,element_list%n_elements
  do i_vertex = 1,n_vertex_max
    i_node = element_list%element(i_elm)%vertex(i_vertex)
    if (node_list%node(i_node)%axis_node) then
      element_list%element(i_elm)%size(i_vertex,6) = 0.d0
      element_list%element(i_elm)%size(i_vertex,7) = 0.d0
      element_list%element(i_elm)%size(i_vertex,8) = 0.d0
      element_list%element(i_elm)%size(i_vertex,9) = 0.d0
      ! --- Loop over all indices (note node_indices have already been calculated above...)
      do k = 1,(n_order+1)/2
        do l = 2,(n_order+1)/2 ! start t-index from 2 to keep only the pure s-derivatives
          node_index = node_indices(k,l)
          if (node_index .le. 9) cycle ! we want only derivatives >=3
          element_list%element(i_elm)%size(i_vertex,node_index) = 0.d0
        enddo
      enddo
    endif
  enddo
enddo


return

end subroutine set_high_order_sizes_on_axis








!> Routine aligns the second derivative control points (ie. i and j)
!> This is not an exact alignment, it looks at the psi-value at s=0.1
!> and compared to the psi-value on the node itself, and moves the control point 
!> until convergence is obtained towards the right psi-value.
subroutine align_2nd_derivatives(node_list,element_list, newnode_list,newelement_list)


use constants
use tr_module
use mod_parameters
use data_structure
use mod_interp
use phys_module, only: R_geo, Z_geo

implicit none

! --- Routine variables
type(type_node_list),    intent(inout) :: node_list,    newnode_list      !< list of nodes with grid information
type(type_element_list), intent(inout) :: element_list, newelement_list   !< list of elements with element information

! --- Local variables
integer, allocatable :: n_parents(:)         ! for each node, want the number of parent elements
integer, allocatable :: node_parents(:,:)    ! for each node, want to know the 4 parent elements
integer, allocatable :: parent_elm_node(:,:) ! for each node, want to know the corresonding vertex for the 4 parent elements
integer              :: i_node, i_elm, i_vertex, i_vertex1, i_vertex2, i
integer              :: i_node_u, i_node_v
integer              :: index
real*8               :: size_u_min, size_v_min
real*8               :: size_tmp
real*8               :: point1(2), point2(2)
real*8               :: distance1, distance2
real*8               :: direction, offset
real*8               :: R_find, Z_find
real*8               :: psi_node, psi_node1, psi_node2
real*8               :: s_base, t_base, psi_base
real*8               :: psi_prev, psi_next
real*8               :: size_base, size_record
real*8               :: size_min, size_max
logical              :: align_i_not_j
integer              :: i_or_j
integer              :: iter, n_tries, i_min
integer              :: n_alignment_loops, i_align
integer, parameter   :: n_range = 10
real*8               :: psi_range(n_range), size_range(n_range), diff_min
real*8, parameter    :: tolerance = 1.d-5 ! convergence tolerance in absolute psi

! --- First, get an approximation of the 2nd derivative alignment, it's a good starting point
! --- Note, this is actually a very good guess in most cases, particularly if
! --- the grid has fine resolution. The further realignment, below, is probably
! --- not even needed...
call approximate_2nd_derivatives(newnode_list,newelement_list)

allocate( n_parents      (  newnode_list%n_nodes) )
allocate( node_parents   (8,newnode_list%n_nodes) )
allocate( parent_elm_node(8,newnode_list%n_nodes) )
n_parents       = 0
node_parents    = 0
parent_elm_node = 0

! --- Find parent elements
do i_node = 1, newnode_list%n_nodes
  n_parents(i_node) = 0
  do i_elm = 1, newelement_list%n_elements
    do i_vertex = 1, n_vertex_max
      if (newelement_list%element(i_elm)%vertex(i_vertex) .eq. i_node) then
        n_parents(i_node) = n_parents(i_node) + 1
        node_parents   (n_parents(i_node),i_node) = i_elm
        parent_elm_node(n_parents(i_node),i_node) = i_vertex
        exit
      endif
    enddo
  enddo
enddo


! --- We look at the alignment at 10% away from node, or whatever distance you choose here
offset = 0.1

! --- Because we change the size at one end of the element side, and then later
! --- on (while doing another node), we might re-align the size at the order end
! --- of that same element, it might modify the first end alignment a bit, so we
! --- need to do this a few time
n_alignment_loops = 2
n_tries = 2 ! number of tries for each node
do i_align = 1,n_alignment_loops

  ! --- Align nodes one by one
  do i_node=1,newnode_list%n_nodes

    ! --- We only care about one side, the one aligned to a psi-surface, so we need to find which side this is
    ! --- psi on node
    R_find = newnode_list%node(i_node)%x(1,1,1)
    Z_find = newnode_list%node(i_node)%x(1,1,2)
    psi_node = get_psi_on_point(node_list, element_list, R_find, Z_find)
    if (psi_node .eq. 1.d10) cycle
    ! --- just take the first parent element
    i = 1
    i_elm = node_parents(i,i_node)
    ! --- Always set vertex1 to be along u (and i), and vertex2 along v (and j)
    if (parent_elm_node(i,i_node) .eq. 1) then
      i_vertex = 1 ; i_vertex1 = 2 ; i_vertex2 = 4
    elseif (parent_elm_node(i,i_node) .eq. 2) then
      i_vertex = 2 ; i_vertex1 = 1 ; i_vertex2 = 3
    elseif (parent_elm_node(i,i_node) .eq. 3) then
      i_vertex = 3 ; i_vertex1 = 4 ; i_vertex2 = 2
    elseif (parent_elm_node(i,i_node) .eq. 4) then
      i_vertex = 4 ; i_vertex1 = 3 ; i_vertex2 = 1
    endif
    ! --- psi on node1
    i_node_u  = newelement_list%element(i_elm)%vertex(i_vertex1)
    R_find    = newnode_list%node(i_node_u)%x(1,1,1)
    Z_find    = newnode_list%node(i_node_u)%x(1,1,2)
    psi_node1 = get_psi_on_point(node_list, element_list, R_find, Z_find)
    ! --- psi on node2
    i_node_v  = newelement_list%element(i_elm)%vertex(i_vertex2)
    R_find    = newnode_list%node(i_node_v)%x(1,1,1)
    Z_find    = newnode_list%node(i_node_v)%x(1,1,2)
    psi_node2 = get_psi_on_point(node_list, element_list, R_find, Z_find)
    if ( (psi_node1 .eq. 1.d10) .and. (psi_node2 .eq. 1.d10) ) cycle
    ! --- compare values
    if ( abs(psi_node1-psi_node) .lt. abs(psi_node2-psi_node) ) then
      align_i_not_j = .true.
      i_or_j        = 5
    else
      align_i_not_j = .false.
      i_or_j        = 6
    endif
    ! --- Set the (s,t) coords where we will do our alignment
    if (parent_elm_node(i,i_node) .eq. 1) then
      if (align_i_not_j) then
        s_base = offset ; t_base = 0.0
      else
        s_base = 0.0    ; t_base = offset
      endif
    elseif (parent_elm_node(i,i_node) .eq. 2) then
      if (align_i_not_j) then
        s_base = 1.0-offset ; t_base = 0.0
      else
        s_base = 1.0        ; t_base = offset
      endif
    elseif (parent_elm_node(i,i_node) .eq. 3) then
      if (align_i_not_j) then
        s_base = 1.0-offset ; t_base = 1.0
      else
        s_base = 1.0        ; t_base = 1.0-offset
      endif
    elseif (parent_elm_node(i,i_node) .eq. 4) then
      if (align_i_not_j) then
        s_base = offset ; t_base = 1.0
      else
        s_base = 0.0    ; t_base = 1.0-offset
      endif
    endif
 
    ! --- Get current psi_value at base
    call interp_RZ(newnode_list,newelement_list,i_elm,s_base,t_base,R_find,Z_find)
    psi_base = get_psi_on_point(node_list, element_list, R_find, Z_find)
    if (psi_base .eq. 1.d10) cycle
    if ( abs(psi_base-psi_node) .lt. tolerance ) cycle
    ! --- Save element size
    size_base = newelement_list%element(i_elm)%size(i_vertex,i_or_j)
    size_record = size_base
    ! --- Start with a range of 30%
    size_min = 0.7*size_base
    size_max = 1.3*size_base
    do iter = 1,n_tries
      diff_min = 1.d10
      do i=1,n_range
        size_range(i) = size_min + float(i-1)/float(n_range-1) * (size_max - size_min)
        newelement_list%element(i_elm)%size(i_vertex,i_or_j) = size_range(i)
        call interp_RZ(newnode_list,newelement_list,i_elm,s_base,t_base,R_find,Z_find)
        psi_range(i) = get_psi_on_point(node_list, element_list, R_find, Z_find)
        if ( abs(psi_range(i)-psi_node) .lt. diff_min ) then
          diff_min = abs(psi_range(i)-psi_node)
          i_min = i
        endif
      enddo
      ! --- Determine next range
      if (i_min .eq. 1) then
        size_min = 0.7 * size_min
        size_max = size_range(2)
      elseif (i_min .eq. n_range) then
        size_min = size_range(n_range-1)
        size_max = 1.3*size_max
      else
        size_min = size_range(i_min-1)
        size_max = size_range(i_min+1)
      endif
      !write(*,'(A,4i6,3f18.7)')'alignment',i_align,i_node,iter,i_min,psi_node,psi_base,psi_range(i_min)
      if (diff_min .lt. abs(psi_base-psi_node) ) size_base = size_range(i_min)
      if (diff_min .lt. tolerance ) exit
    enddo
    !write(*,'(A,2i6,2f18.7)')'alignment',i_align,i_node,size_record,size_base
    newelement_list%element(i_elm)%size(i_vertex,i_or_j) = size_base

    ! --- Size need to be the same on all parent elements!
    do i = 2,n_parents(i_node)
      i_elm    = node_parents(i,i_node)
      i_vertex = parent_elm_node(i,i_node)
      newelement_list%element(i_elm)%size(i_vertex,i_or_j) = size_base
    enddo

  enddo ! end i_node loop

enddo ! end alignment_loops


deallocate(n_parents      )
deallocate(node_parents   )
deallocate(parent_elm_node)



return

end subroutine align_2nd_derivatives


















!> This routine approximates the second derivative control points (ie. vectors i and j)
!> such that the curvature changes linearly between two vertices of an element
!> Therefore, we assume that the _ss derivative is not known on the nodes, but
!> by assuming that, at s=0.5, the _ss value is 0.5*(node1_ss + node2_ss), this
!> gives us an approximation for the vectors i and j (Please see paper for the formula)
subroutine approximate_2nd_derivatives(node_list,element_list)


use constants
use tr_module
use mod_parameters
use data_structure
use phys_module, only: R_geo, Z_geo

implicit none

! --- Routine variables
type(type_node_list),    intent(inout) :: node_list      !< list of nodes with grid information
type(type_element_list), intent(inout) :: element_list   !< list of elements with element information

! --- Local variables
integer, allocatable :: n_parents(:)         ! for each node, want the number of parent elements
integer, allocatable :: node_parents(:,:)    ! for each node, want to know the 4 parent elements
integer, allocatable :: parent_elm_node(:,:) ! for each node, want to know the corresonding vertex for the 4 parent elements
integer              :: i_node, i_elm, i_vertex, i_vertex1, i_vertex2, i
integer              :: i_node_u, i_node_v
integer              :: index
real*8               :: size_u_min, size_v_min
real*8               :: size_tmp
real*8               :: point1(2), point2(2)
real*8               :: distance1, distance2
real*8               :: direction

allocate( n_parents      (  node_list%n_nodes) )
allocate( node_parents   (8,node_list%n_nodes) )
allocate( parent_elm_node(8,node_list%n_nodes) )
n_parents       = 0
node_parents    = 0
parent_elm_node = 0

! --- Find parent elements
do i_node = 1, node_list%n_nodes
  n_parents(i_node) = 0
  do i_elm = 1, element_list%n_elements
    do i_vertex = 1, n_vertex_max
      if (element_list%element(i_elm)%vertex(i_vertex) .eq. i_node) then
        n_parents(i_node) = n_parents(i_node) + 1
        node_parents   (n_parents(i_node),i_node) = i_elm
        parent_elm_node(n_parents(i_node),i_node) = i_vertex
        exit
      endif
    enddo
  enddo
enddo




! --- After definition of u,v,w, we are ready for i and j.
do i_node = 1, node_list%n_nodes
  
  ! ------------
  ! --- Vector i
  ! ------------
  
  ! --- METHOD: Assume curvature in the middle of edge is 0.5(curv_00+curv_05) 
  ! --- Find the next node
  i = 1
  i_elm = node_parents(i,i_node)
  if (parent_elm_node(i,i_node) .eq. 1) then
    i_vertex1 = 1 ; i_vertex2 = 2
    i_node_u = element_list%element(i_elm)%vertex(2)
  elseif (parent_elm_node(i,i_node) .eq. 2) then
    i_vertex1 = 2 ; i_vertex2 = 1
    i_node_u = element_list%element(i_elm)%vertex(1)
  elseif (parent_elm_node(i,i_node) .eq. 3) then
    i_vertex1 = 3 ; i_vertex2 = 4
    i_node_u = element_list%element(i_elm)%vertex(4)
  elseif (parent_elm_node(i,i_node) .eq. 4) then
    i_vertex1 = 4 ; i_vertex2 = 3
    i_node_u = element_list%element(i_elm)%vertex(3)
  endif
  ! --- Formulation: h_i*\vec{i} = -0.25 * ( h_00^u*\vec{u}_00 + h_50^u*\vec{u}_50 )
  node_list%node(i_node)%X(1,5,1) = -0.25 * ( element_list%element(i_elm)%size(i_vertex1,2)*node_list%node(i_node  )%x(1,2,1) &
                                             +element_list%element(i_elm)%size(i_vertex2,2)*node_list%node(i_node_u)%x(1,2,1) )
  node_list%node(i_node)%X(1,5,2) = -0.25 * ( element_list%element(i_elm)%size(i_vertex1,2)*node_list%node(i_node  )%x(1,2,2) &
                                             +element_list%element(i_elm)%size(i_vertex2,2)*node_list%node(i_node_u)%x(1,2,2) )
  ! --- Take the average from both sides of the node
  if (n_parents(i_node) .gt. 1) then
    do i = 2,n_parents(i_node)
      i_elm = node_parents(i,i_node)
      i_vertex2 = 0
      if ( (i_vertex1 .eq. 1) .or. (i_vertex1 .eq. 4) ) then
        if (parent_elm_node(i,i_node) .eq. 2) then
          i_vertex1 = 2 ; i_vertex2 = 1
          i_node_u = element_list%element(i_elm)%vertex(1)
        elseif (parent_elm_node(i,i_node) .eq. 3) then
          i_vertex1 = 3 ; i_vertex2 = 4
          i_node_u = element_list%element(i_elm)%vertex(4)
        endif
      elseif ( (i_vertex1 .eq. 2) .or. (i_vertex1 .eq. 3) ) then
        if (parent_elm_node(i,i_node) .eq. 1) then
          i_vertex1 = 1 ; i_vertex2 = 2
          i_node_u = element_list%element(i_elm)%vertex(2)
        elseif (parent_elm_node(i,i_node) .eq. 4) then
          i_vertex1 = 4 ; i_vertex2 = 3
          i_node_u = element_list%element(i_elm)%vertex(3)
        endif
      endif
      if (i_vertex2 .eq. 0) then
        cycle
      else
        ! --- Formulation: h_i*\vec{i} = -0.25 * ( h_00^u*\vec{u}_00 + h_50^u*\vec{u}_50 )
        node_list%node(i_node)%X(1,5,1) = 0.5 * node_list%node(i_node)%X(1,5,1) &
                                          -0.5*0.25 * ( element_list%element(i_elm)%size(i_vertex1,2)*node_list%node(i_node  )%x(1,2,1) &
                                                       +element_list%element(i_elm)%size(i_vertex2,2)*node_list%node(i_node_u)%x(1,2,1) )
        node_list%node(i_node)%X(1,5,2) = 0.5 * node_list%node(i_node)%X(1,5,2) &
                                          -0.5*0.25 * ( element_list%element(i_elm)%size(i_vertex1,2)*node_list%node(i_node  )%x(1,2,2) &
                                                       +element_list%element(i_elm)%size(i_vertex2,2)*node_list%node(i_node_u)%x(1,2,2) )
        exit
      endif
    enddo
  endif
  
  ! ------------
  ! --- Vector j
  ! ------------
  
  ! --- METHOD: Assume curvature in the middle of edge is 0.5(curv_00+curv_05) 
  ! --- Find the next node
  i = 1
  i_elm = node_parents(i,i_node)
  if (parent_elm_node(i,i_node) .eq. 1) then
    i_vertex1 = 1 ; i_vertex2 = 4
    i_node_v = element_list%element(i_elm)%vertex(4)
  elseif (parent_elm_node(i,i_node) .eq. 2) then
    i_vertex1 = 2 ; i_vertex2 = 3
    i_node_v = element_list%element(i_elm)%vertex(3)
  elseif (parent_elm_node(i,i_node) .eq. 3) then
    i_vertex1 = 3 ; i_vertex2 = 2
    i_node_v = element_list%element(i_elm)%vertex(2)
  elseif (parent_elm_node(i,i_node) .eq. 4) then
    i_vertex1 = 4 ; i_vertex2 = 1
    i_node_v = element_list%element(i_elm)%vertex(1)
  endif
  ! --- Formulation: h_j*\vec{j} = -0.25 * ( h_00^v*\vec{v}_00 + h_05^v*\vec{v}_05 )
  node_list%node(i_node)%X(1,6,1) = -0.25 * ( element_list%element(i_elm)%size(i_vertex1,3)*node_list%node(i_node  )%x(1,3,1) &
                                             +element_list%element(i_elm)%size(i_vertex2,3)*node_list%node(i_node_v)%x(1,3,1) )
  node_list%node(i_node)%X(1,6,2) = -0.25 * ( element_list%element(i_elm)%size(i_vertex1,3)*node_list%node(i_node  )%x(1,3,2) &
                                             +element_list%element(i_elm)%size(i_vertex2,3)*node_list%node(i_node_v)%x(1,3,2) )
  ! --- Take the average from both sides of the node
  if (n_parents(i_node) .gt. 1) then
    do i = 2,n_parents(i_node)
      i_elm = node_parents(i,i_node)
      i_vertex2 = 0
      if ( (i_vertex1 .eq. 1) .or. (i_vertex1 .eq. 2) ) then
        if (parent_elm_node(i,i_node) .eq. 3) then
          i_vertex1 = 3 ; i_vertex2 = 2
          i_node_v = element_list%element(i_elm)%vertex(2)
        elseif (parent_elm_node(i,i_node) .eq. 4) then
          i_vertex1 = 4 ; i_vertex2 = 1
          i_node_v = element_list%element(i_elm)%vertex(1)
        endif
      elseif ( (i_vertex1 .eq. 3) .or. (i_vertex1 .eq. 4) ) then
        if (parent_elm_node(i,i_node) .eq. 1) then
          i_vertex1 = 1 ; i_vertex2 = 4
          i_node_v = element_list%element(i_elm)%vertex(4)
        elseif (parent_elm_node(i,i_node) .eq. 2) then
          i_vertex1 = 2 ; i_vertex2 = 3
          i_node_v = element_list%element(i_elm)%vertex(3)
        endif
      endif
      if (i_vertex2 .eq. 0) then
        cycle
      else
        ! --- Formulation: h_i*\vec{i} = -0.25 * ( h_00^u*\vec{u}_00 + h_50^u*\vec{u}_50 )
        node_list%node(i_node)%X(1,6,1) = 0.5 * node_list%node(i_node)%X(1,6,1) &
                                          -0.5*0.25 * ( element_list%element(i_elm)%size(i_vertex1,3)*node_list%node(i_node  )%x(1,3,1) &
                                                       +element_list%element(i_elm)%size(i_vertex2,3)*node_list%node(i_node_v)%x(1,3,1) )
        node_list%node(i_node)%X(1,6,2) = 0.5 * node_list%node(i_node)%X(1,6,2) &
                                          -0.5*0.25 * ( element_list%element(i_elm)%size(i_vertex1,3)*node_list%node(i_node  )%x(1,3,2) &
                                                       +element_list%element(i_elm)%size(i_vertex2,3)*node_list%node(i_node_v)%x(1,3,2) )
        exit
      endif
    enddo
  endif
  

enddo

deallocate(n_parents      )
deallocate(node_parents   )
deallocate(parent_elm_node)

return

end subroutine approximate_2nd_derivatives













































subroutine transform_to_bi_quintic(node_list,element_list)


use constants
use tr_module
use mod_parameters
use data_structure
use phys_module, only: R_geo, Z_geo

implicit none

! --- Routine variables
type(type_node_list),    intent(inout) :: node_list      !< list of nodes with grid information
type(type_element_list), intent(inout) :: element_list   !< list of elements with element information

! --- Local variables
integer, allocatable          :: n_parents(:)         ! for each node, want the number of parent elements
integer, allocatable          :: node_parents(:,:)    ! for each node, want to know the 4 parent elements
integer, allocatable          :: parent_elm_node(:,:) ! for each node, want to know the corresonding vertex for the 4 parent elements
integer                       :: i_node, i_elm, i_vertex, i_vertex1, i_vertex2, i
integer                       :: i_node_u, i_node_v
integer                       :: index
real*8                        :: size_u_min, size_v_min
real*8                        :: size_tmp
real*8                        :: point1(2), point2(2)
real*8                        :: distance1, distance2
real*8                        :: direction
real*8                        :: scale_uv, scale_ij, scale_wk

allocate( n_parents      (  node_list%n_nodes) )
allocate( node_parents   (4,node_list%n_nodes) )
allocate( parent_elm_node(4,node_list%n_nodes) )
n_parents       = 0
node_parents    = 0
parent_elm_node = 0

! --- Find parent elements
do i_node = 1, node_list%n_nodes

  n_parents(i_node) = 0
  do i_elm = 1, element_list%n_elements
  
    do i_vertex = 1, n_vertex_max
    
      if (element_list%element(i_elm)%vertex(i_vertex) .eq. i_node) then
        n_parents(i_node) = n_parents(i_node) + 1
        node_parents   (n_parents(i_node),i_node) = i_elm
        parent_elm_node(n_parents(i_node),i_node) = i_vertex
        exit
      endif
    
    enddo
    
    if ( (node_list%node(i_node)%boundary .eq. 2) .and. (n_parents(i_node) .eq. 2) ) exit
    if ( (node_list%node(i_node)%boundary .ne. 2) .and. (n_parents(i_node) .eq. 4) ) exit
  
  enddo
  
enddo


! --- For cubic elements, we scale by 1/3, for quintic elements, by 1/5
scale_uv = 0.2
! scale by 1/10 of the element side?
scale_ij = 0.001 !1.0
! scale by 1/10 of the element side?
scale_wk = 0.0
  

! --- Redefine all vectors at each node
do i_node = 1, node_list%n_nodes

  ! --- Note: all vectors are always of unit size
  ! --- it's the element size that changes


  ! -------------------
  ! --- Vectors u and v
  ! -------------------

  ! --- We must freeze the u/v vector sizes, which means we 
  ! --- must make sure the vector size is based on the smallest element
  
  ! --- Loop over each parent element and get minimal size
  size_u_min = 1.d15
  size_v_min = 1.d15
  do i=1,n_parents(i_node)
    i_elm = node_parents(i,i_node)
    if (parent_elm_node(i,i_node) .eq. 1) then
      i_node_u = element_list%element(i_elm)%vertex(2)
      i_node_v = element_list%element(i_elm)%vertex(4)
    elseif (parent_elm_node(i,i_node) .eq. 2) then
      i_node_u = element_list%element(i_elm)%vertex(1)
      i_node_v = element_list%element(i_elm)%vertex(3)
    elseif (parent_elm_node(i,i_node) .eq. 3) then
      i_node_u = element_list%element(i_elm)%vertex(4)
      i_node_v = element_list%element(i_elm)%vertex(2)
    elseif (parent_elm_node(i,i_node) .eq. 4) then
      i_node_u = element_list%element(i_elm)%vertex(3)
      i_node_v = element_list%element(i_elm)%vertex(1)
    endif
    distance1 = distance_nodes(node_list,i_node,i_node_u)
    distance2 = distance_nodes(node_list,i_node,i_node_v)
    if (distance1 .ne. 0.d0) size_u_min = min(size_u_min, distance1)
    if (distance2 .ne. 0.d0) size_v_min = min(size_v_min, distance2)
  enddo
  if (size_u_min .eq. 1.d15) size_u_min = 0.d0 
  if (size_v_min .eq. 1.d15) size_v_min = 0.d0 
  
  ! --- Set the direction u and v
  ! --- Positive direction of u is always from node-(1->2) === (4->3)
  ! --- Positive direction of v is always from node-(1->4) === (2->3)
  ! --- First, we set a temporary size of the vectors u,v to make sure they
  ! --- are smaller than the element side
  call normalise_vector(node_list,i_node,2)
  call normalise_vector(node_list,i_node,3)
  node_list%node(i_node)%X(1,2,1) = node_list%node(i_node)%X(1,2,1) * size_u_min * scale_uv
  node_list%node(i_node)%X(1,2,2) = node_list%node(i_node)%X(1,2,2) * size_u_min * scale_uv
  node_list%node(i_node)%X(1,3,1) = node_list%node(i_node)%X(1,3,1) * size_v_min * scale_uv
  node_list%node(i_node)%X(1,3,2) = node_list%node(i_node)%X(1,3,2) * size_v_min * scale_uv
  ! --- Then, we check the direction by locating the tip of the vector relative to the nodes
  do i=1,n_parents(i_node)
    i_elm = node_parents(i,i_node)
    if (parent_elm_node(i,i_node) .eq. 1) then
      i_node_u = element_list%element(i_elm)%vertex(2)
      i_node_v = element_list%element(i_elm)%vertex(4)
    elseif (parent_elm_node(i,i_node) .eq. 2) then
      i_node_u = element_list%element(i_elm)%vertex(1)
      i_node_v = element_list%element(i_elm)%vertex(3)
    elseif (parent_elm_node(i,i_node) .eq. 3) then
      i_node_u = element_list%element(i_elm)%vertex(4)
      i_node_v = element_list%element(i_elm)%vertex(2)
    elseif (parent_elm_node(i,i_node) .eq. 4) then
      i_node_u = element_list%element(i_elm)%vertex(3)
      i_node_v = element_list%element(i_elm)%vertex(1)
    endif
    ! --- compare size of side with distance from tip of vector u
    distance1 = distance_nodes(node_list,i_node,i_node_u)
    point1(1) = node_list%node(i_node)%X(1,1,1) + node_list%node(i_node)%X(1,2,1) ! tip of vector u
    point1(2) = node_list%node(i_node)%X(1,1,2) + node_list%node(i_node)%X(1,2,2) ! tip of vector u
    distance2 = distance_node_point(node_list,i_node_u,point1)
    ! --- reverse vector u if needed
    if (distance2 .lt. distance1) then
      if ( (parent_elm_node(i,i_node) .eq. 2) .or. (parent_elm_node(i,i_node) .eq. 3) ) then
        node_list%node(i_node)%X(1,2,1) = - node_list%node(i_node)%X(1,2,1)
        node_list%node(i_node)%X(1,2,2) = - node_list%node(i_node)%X(1,2,2)
      endif
    endif
    ! --- compare size of side with distance from tip of vector v
    distance1 = distance_nodes(node_list,i_node,i_node_v)
    point1(1) = node_list%node(i_node)%X(1,1,1) + node_list%node(i_node)%X(1,3,1) ! tip of vector v
    point1(2) = node_list%node(i_node)%X(1,1,2) + node_list%node(i_node)%X(1,3,2) ! tip of vector v
    distance2 = distance_node_point(node_list,i_node_v,point1)
    ! --- reverse vector v if needed
    if (distance2 .lt. distance1) then
      if ( (parent_elm_node(i,i_node) .eq. 3) .or. (parent_elm_node(i,i_node) .eq. 4) ) then
        node_list%node(i_node)%X(1,3,1) = - node_list%node(i_node)%X(1,3,1)
        node_list%node(i_node)%X(1,3,2) = - node_list%node(i_node)%X(1,3,2)
      endif
    endif
  enddo
  ! --- Finally, we set u,v back to unit vectors
  call normalise_vector(node_list,i_node,2)
  call normalise_vector(node_list,i_node,3)
  
  ! --- SIZE CONDITION: h_u and h_v must be the same on all parent elements (ie. breaks Bezier to become Hermite)
  ! --- Set the element size for nodes u and v
  do i=1,n_parents(i_node)
    i_elm = node_parents(i,i_node)
    i_vertex = parent_elm_node(i,i_node)
    if (i_vertex .eq. 1) then
      element_list%element(i_elm)%size(i_vertex,2) = + size_u_min * scale_uv
      element_list%element(i_elm)%size(i_vertex,3) = + size_v_min * scale_uv
    elseif (i_vertex .eq. 2) then
      element_list%element(i_elm)%size(i_vertex,2) = - size_u_min * scale_uv
      element_list%element(i_elm)%size(i_vertex,3) = + size_v_min * scale_uv
    elseif (i_vertex .eq. 3) then
      element_list%element(i_elm)%size(i_vertex,2) = - size_u_min * scale_uv
      element_list%element(i_elm)%size(i_vertex,3) = - size_v_min * scale_uv
    elseif (i_vertex .eq. 4) then
      element_list%element(i_elm)%size(i_vertex,2) = + size_u_min * scale_uv
      element_list%element(i_elm)%size(i_vertex,3) = - size_v_min * scale_uv
    endif
  enddo
  
  
  ! ------------
  ! --- Vector w
  ! ------------
  
  ! --- There is no rule for vector w, so we just set it as
  ! --- w = u+v, of length scale_wk*(u+v). But the condition comes on the element size
  
  ! --- Vector definition
  node_list%node(i_node)%X(1,4,1) = scale_wk * ( node_list%node(i_node)%X(1,2,1) + node_list%node(i_node)%X(1,3,1) )
  node_list%node(i_node)%X(1,4,2) = scale_wk * ( node_list%node(i_node)%X(1,2,2) + node_list%node(i_node)%X(1,3,2) )
  ! --- SIZE CONDITION: h_w = h_u*h_v
  do i = 1,n_parents(i_node)
    i_elm    = node_parents(i,i_node)
    i_vertex = parent_elm_node(i,i_node)
    element_list%element(i_elm)%size(i_vertex,4) = element_list%element(i_elm)%size(i_vertex,2) * element_list%element(i_elm)%size(i_vertex,3)
  enddo

enddo


! VECTORS i AND j NEED TO DETERMINE THE CURVATURE! THEY NEED TO BE SET DEPENDING ON ELEMENT SHAPE!!!
! VECTORS i AND j NEED TO DETERMINE THE CURVATURE! THEY NEED TO BE SET DEPENDING ON ELEMENT SHAPE!!!
! VECTORS i AND j NEED TO DETERMINE THE CURVATURE! THEY NEED TO BE SET DEPENDING ON ELEMENT SHAPE!!!
! VECTORS i AND j NEED TO DETERMINE THE CURVATURE! THEY NEED TO BE SET DEPENDING ON ELEMENT SHAPE!!!
! VECTORS i AND j NEED TO DETERMINE THE CURVATURE! THEY NEED TO BE SET DEPENDING ON ELEMENT SHAPE!!!
! VECTORS i AND j NEED TO DETERMINE THE CURVATURE! THEY NEED TO BE SET DEPENDING ON ELEMENT SHAPE!!!
! VECTORS i AND j NEED TO DETERMINE THE CURVATURE! THEY NEED TO BE SET DEPENDING ON ELEMENT SHAPE!!!
! VECTORS i AND j NEED TO DETERMINE THE CURVATURE! THEY NEED TO BE SET DEPENDING ON ELEMENT SHAPE!!!
! VECTORS i AND j NEED TO DETERMINE THE CURVATURE! THEY NEED TO BE SET DEPENDING ON ELEMENT SHAPE!!!
! VECTORS i AND j NEED TO DETERMINE THE CURVATURE! THEY NEED TO BE SET DEPENDING ON ELEMENT SHAPE!!!

! --- After definition of u,v,w, we are ready for i and j.
do i_node = 1, node_list%n_nodes
  
  ! ------------
  ! --- Vector i
  ! ------------
  
  ! --- Like w, we define i as i = u+v, of length scale_ij*(u+v)
  ! --- The real condition comes on the element size
  ! --- Which simply that the size should be exactly the same on all elements
  ! --- we set the size to be small, scale_ij*(size_u+size_v) of the first node available,
  ! --- with scale_i = 0.1 ?
  
  ! --- Vector definition
  node_list%node(i_node)%X(1,5,1) = 0.d0 ! node_list%node(i_node)%X(1,2,1) + node_list%node(i_node)%X(1,3,1)
  node_list%node(i_node)%X(1,5,2) = 0.d0 ! node_list%node(i_node)%X(1,2,2) + node_list%node(i_node)%X(1,3,2)
  ! --- SIZE CONDITION: h_i must be the same on all parent nodes
  i_elm = node_parents(1,i_node)
  size_tmp = 0.5 * ( abs(element_list%element(i_elm)%size(1,2)) + abs(element_list%element(i_elm)%size(1,3)) )

  ! --- METHOD: Assume curvature in the middle of edge is 0.5(curv_00+curv_05) 
  ! --- Find the next node
  i = n_parents(1)
  i_elm = node_parents(i,i_node)
  if (parent_elm_node(i,i_node) .eq. 1) then
    i_vertex1 = 1 ; i_vertex2 = 2
    i_node_u = element_list%element(i_elm)%vertex(2)
  elseif (parent_elm_node(i,i_node) .eq. 2) then
    i_vertex1 = 2 ; i_vertex2 = 1
    i_node_u = element_list%element(i_elm)%vertex(1)
  elseif (parent_elm_node(i,i_node) .eq. 3) then
    i_vertex1 = 3 ; i_vertex2 = 4
    i_node_u = element_list%element(i_elm)%vertex(4)
  elseif (parent_elm_node(i,i_node) .eq. 4) then
    i_vertex1 = 4 ; i_vertex2 = 3
    i_node_u = element_list%element(i_elm)%vertex(3)
  endif
  ! --- Formulation: h_i*\vec{i} = -0.25 * ( h_00^u*\vec{u}_00 + h_50^u*\vec{u}_50 )
  point1(1) = -0.25 * ( element_list%element(i_elm)%size(i_vertex1,2)*node_list%node(i_node  )%x(1,2,1) &
                       +element_list%element(i_elm)%size(i_vertex2,2)*node_list%node(i_node_u)%x(1,2,1) )
  point1(2) = -0.25 * ( element_list%element(i_elm)%size(i_vertex1,2)*node_list%node(i_node  )%x(1,2,2) &
                       +element_list%element(i_elm)%size(i_vertex2,2)*node_list%node(i_node_u)%x(1,2,2) )
  point2(1) = 0.d0 ; point2(2) = 0.d0
  distance1 = distance_points(point1,point2)
  node_list%node(i_node)%X(1,5,1) = point1(1)
  node_list%node(i_node)%X(1,5,2) = point1(2)
  call normalise_vector(node_list,i_node,5)
  size_tmp = distance1
  do i = 1,n_parents(i_node)
    i_elm    = node_parents(i,i_node)
    i_vertex = parent_elm_node(i,i_node)
    element_list%element(i_elm)%size(i_vertex,5) = size_tmp * scale_ij
  enddo
  
  
  ! ------------
  ! --- Vector j
  ! ------------
  
  ! --- METHOD: Assume curvature in the middle of edge is 0.5(curv_00+curv_05) 
  ! --- Find the next node
  i = n_parents(1)
  i_elm = node_parents(i,i_node)
  if (parent_elm_node(i,i_node) .eq. 1) then
    i_vertex1 = 1 ; i_vertex2 = 4
    i_node_v = element_list%element(i_elm)%vertex(4)
  elseif (parent_elm_node(i,i_node) .eq. 2) then
    i_vertex1 = 2 ; i_vertex2 = 3
    i_node_v = element_list%element(i_elm)%vertex(3)
  elseif (parent_elm_node(i,i_node) .eq. 3) then
    i_vertex1 = 3 ; i_vertex2 = 2
    i_node_v = element_list%element(i_elm)%vertex(2)
  elseif (parent_elm_node(i,i_node) .eq. 4) then
    i_vertex1 = 4 ; i_vertex2 = 1
    i_node_v = element_list%element(i_elm)%vertex(1)
  endif
  ! --- Formulation: h_j*\vec{j} = -0.25 * ( h_00^v*\vec{v}_00 + h_05^v*\vec{v}_05 )
  point1(1) = -0.25 * ( element_list%element(i_elm)%size(i_vertex1,3)*node_list%node(i_node  )%x(1,3,1) &
                       +element_list%element(i_elm)%size(i_vertex2,3)*node_list%node(i_node_v)%x(1,3,1) )
  point1(2) = -0.25 * ( element_list%element(i_elm)%size(i_vertex1,3)*node_list%node(i_node  )%x(1,3,2) &
                       +element_list%element(i_elm)%size(i_vertex2,3)*node_list%node(i_node_v)%x(1,3,2) )
  point2(1) = 0.d0 ; point2(2) = 0.d0
  distance1 = distance_points(point1,point2)
  node_list%node(i_node)%X(1,6,1) = element_list%element(i_elm)%size(i_vertex1,2) * node_list%node(i_node)%X(1,2,1) !point1(1)
  node_list%node(i_node)%X(1,6,2) = element_list%element(i_elm)%size(i_vertex1,2) * node_list%node(i_node)%X(1,2,2) !point1(2)
  call normalise_vector(node_list,i_node,6)
  ! --- SIZE CONDITION: h_j must be the same on all parent nodes
  size_tmp = distance1
  do i = 1,n_parents(i_node)
    i_elm    = node_parents(i,i_node)
    i_vertex = parent_elm_node(i,i_node)
    element_list%element(i_elm)%size(i_vertex,6) = size_tmp * scale_ij
  enddo
  
  
  ! ------------
  ! --- Vector m
  ! ------------
  
  ! --- m needs to be across from u, so we simply set m = v, of unit size
  ! --- and we set the size the same as that of the v vector
  
  ! --- Vector definition
  node_list%node(i_node)%X(1,7,1) = node_list%node(i_node)%X(1,3,1)
  node_list%node(i_node)%X(1,7,2) = node_list%node(i_node)%X(1,3,2)
  call normalise_vector(node_list,i_node,7)
  ! --- SIZE CONDITION: h_m must opposite on either side of the two parent nodes (like v)
  do i = 1,n_parents(i_node)
    i_elm    = node_parents(i,i_node)
    i_vertex = parent_elm_node(i,i_node)
    element_list%element(i_elm)%size(i_vertex,7) = element_list%element(i_elm)%size(i_vertex,3)
  enddo
  
  
  ! ------------
  ! --- Vector n
  ! ------------
  
  ! --- n needs to be across from v, so we simply set n = u, of unit size
  ! --- and we set the size the same as that of the u vector
  
  ! --- Vector definition
  node_list%node(i_node)%X(1,8,1) = node_list%node(i_node)%X(1,2,1)
  node_list%node(i_node)%X(1,8,2) = node_list%node(i_node)%X(1,2,2)
  call normalise_vector(node_list,i_node,8)
  ! --- SIZE CONDITION: h_n must opposite on either side of the two parent nodes (like u)
  do i = 1,n_parents(i_node)
    i_elm    = node_parents(i,i_node)
    i_vertex = parent_elm_node(i,i_node)
    element_list%element(i_elm)%size(i_vertex,8) = element_list%element(i_elm)%size(i_vertex,2)
  enddo
  
  
  ! ------------
  ! --- Vector k
  ! ------------
  
  ! --- Because it is very similar to vector w, we simply set k = w
  ! --- and we set the size the same as that of the u vector
  
  ! --- Vector definition
  node_list%node(i_node)%X(1,9,1) = node_list%node(i_node)%X(1,4,1)
  node_list%node(i_node)%X(1,9,2) = node_list%node(i_node)%X(1,4,2)
  ! --- SIZE CONDITION: h_k = h_i*h_j       ! OLD: (h_u*h_v)**2 = h_w**2
  do i = 1,n_parents(i_node)
    i_elm    = node_parents(i,i_node)
    i_vertex = parent_elm_node(i,i_node)
    element_list%element(i_elm)%size(i_vertex,9) = element_list%element(i_elm)%size(i_vertex,5) *  element_list%element(i_elm)%size(i_vertex,6)
  enddo

  ! --- IMPORTANT NOTE: We assume that "force_central_node" is used
  ! --- This means we identify axis nodes by checking the %index(1)=1 nodes
  ! --- On these noes, the poloidal vector sizes should be zero, ie. vectors v,w,j,n,k
  ! --- IMPORTANT NOTE: I'm not sure about n and k!
  if (node_list%node(i_node)%index(1) .eq. 1) then
    do i = 1,n_parents(i_node)
      i_elm    = node_parents(i,i_node)
      i_vertex = parent_elm_node(i,i_node)
      element_list%element(i_elm)%size(i_vertex,3) = 0.d0
      element_list%element(i_elm)%size(i_vertex,4) = 0.d0
      element_list%element(i_elm)%size(i_vertex,6) = 0.d0
      !element_list%element(i_elm)%size(i_vertex,8) = 0.d0
      !element_list%element(i_elm)%size(i_vertex,9) = 0.d0
    enddo
  endif


enddo





! YOU NEED TO TAKE CARE OF THE XPOINT AS WELL!!!
! YOU NEED TO TAKE CARE OF THE XPOINT AS WELL!!!
! YOU NEED TO TAKE CARE OF THE XPOINT AS WELL!!!
! YOU NEED TO TAKE CARE OF THE XPOINT AS WELL!!!
! YOU NEED TO TAKE CARE OF THE XPOINT AS WELL!!!
! YOU NEED TO TAKE CARE OF THE XPOINT AS WELL!!!
! YOU NEED TO TAKE CARE OF THE XPOINT AS WELL!!!
! YOU NEED TO TAKE CARE OF THE XPOINT AS WELL!!!
! YOU NEED TO TAKE CARE OF THE XPOINT AS WELL!!!

! --- Redefine node indexes in the matrix
index = 0
do i_node = 1, node_list%n_nodes
  
  ! --- Careful with force_axis_nodes
  if (node_list%node(i_node)%index(1) .eq. 1) then
    if (index .eq. 0) index = 1
    do i=2,n_degrees
      node_list%node(i_node)%index(i) = index + i-1
    enddo
    index = index + n_degrees-1
  else
    do i=1,n_degrees
      node_list%node(i_node)%index(i) = index + i
    enddo
    index = index + n_degrees
  endif
enddo


return

end subroutine transform_to_bi_quintic








pure real*8 function distance_nodes(node_list,i_node1,i_node2)
use data_structure
implicit none
type(type_node_list), intent(in) :: node_list
integer,              intent(in) :: i_node1,i_node2
distance_nodes = sqrt(  (node_list%node(i_node2)%X(1,1,1) - node_list%node(i_node1)%X(1,1,1))**2 &
                      + (node_list%node(i_node2)%X(1,1,2) - node_list%node(i_node1)%X(1,1,2))**2 )
end function distance_nodes

pure real*8 function distance_node_point(node_list,i_node,point)
use data_structure
implicit none
type(type_node_list), intent(in) :: node_list
integer,              intent(in) :: i_node
real*8,               intent(in) :: point(2)
distance_node_point = sqrt(  (point(1) - node_list%node(i_node)%X(1,1,1))**2 &
                           + (point(2) - node_list%node(i_node)%X(1,1,2))**2 )
end function distance_node_point

pure real*8 function distance_points(point1,point2)
use data_structure
implicit none
real*8, intent(in) :: point1(2),point2(2)
distance_points = sqrt(  (point2(1) - point1(1))**2 + (point2(2) - point1(2))**2 )
end function distance_points

pure subroutine normalise_vector(node_list,i_node,i_vector)
use data_structure
implicit none
type(type_node_list), intent(inout) :: node_list
integer,              intent(in)    :: i_node, i_vector
real*8                              :: vector_size
  vector_size = sqrt( node_list%node(i_node)%X(1,i_vector,1)**2.0 + node_list%node(i_node)%X(1,i_vector,2)**2.0 )
  if (vector_size .eq. 0.d0) return ! yes, because with fix_axis, we set v and w to zero!
  node_list%node(i_node)%X(1,i_vector,1) = node_list%node(i_node)%X(1,i_vector,1) / vector_size
  node_list%node(i_node)%X(1,i_vector,2) = node_list%node(i_node)%X(1,i_vector,2) / vector_size
end subroutine normalise_vector

real*8 function get_psi_on_point(node_list, element_list, R_find, Z_find)
use data_structure
use mod_interp
implicit none
type(type_node_list),    intent(in) :: node_list
type(type_element_list), intent(in) :: element_list 
real*8,                  intent(in) :: R_find, Z_find
integer :: ifail, i_elm_out
real*8  :: R_out,Z_out,s_out,t_out
real*8  :: psi, dd1,dd2,dd3,dd4,dd5
  get_psi_on_point = 1.d10
  call find_RZ(node_list,element_list,R_find,Z_find,R_out,Z_out,i_elm_out,s_out,t_out,ifail)
  if (ifail .ne. 0) return
  call interp(node_list,element_list,i_elm_out,1,1,s_out,t_out,psi,dd1,dd2,dd3,dd4,dd5)
  get_psi_on_point = psi
end function get_psi_on_point





end module
