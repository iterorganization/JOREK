subroutine join_grid_patches(node_list1, element_list1, node_list2, element_list2, newnode_list, newelement_list, xcase)

use constants
use tr_module 
use data_structure
use grid_xpoint_data
use phys_module, only:   tokamak_device

implicit none

! --- Routine parameters
type (type_node_list)       , intent(in)    :: node_list1
type (type_element_list)    , intent(in)    :: element_list1
type (type_node_list)       , intent(in)    :: node_list2
type (type_element_list)    , intent(in)    :: element_list2
type (type_node_list)       , intent(inout) :: newnode_list
type (type_element_list)    , intent(inout) :: newelement_list
integer,                      intent(in)    :: xcase ! It's important to have this as a variable, we might be patching grids outside the main part

! --- local variables
integer             :: i, j, k, index, n_loop, n_start
integer             :: n_zip_nodes
integer             :: zip_nodes(n_nodes_max,2)
integer             :: new_index(n_nodes_max)
real*8              :: R1, Z1
real*8              :: R2, Z2
real*8              :: distance
real*8, parameter   :: tolerance = 1.d-5 ! 0.01mm should be largely enough
logical             :: avoid_this_one
logical, parameter  :: plot_grid = .true.
logical, parameter  :: debug = .false.

write(*,*) '*****************************************'
write(*,*) '* X-point grid inside wall :            *'
write(*,*) '*****************************************'
write(*,*) '                 Joining grid patches together'



! --- Allocate data structures for new nodes and initialize them
newnode_list%n_nodes = 0
newnode_list%n_dof   = 0
do i = 1, n_nodes_max
  newnode_list%node(i)%x           = 0.d0
  newnode_list%node(i)%values      = 0.d0
  newnode_list%node(i)%deltas      = 0.d0
  newnode_list%node(i)%index       = 0
  newnode_list%node(i)%boundary    = 0
  newnode_list%node(i)%axis_node   = .false.
  newnode_list%node(i)%axis_dof    = 0
  newnode_list%node(i)%parents     = 0
  newnode_list%node(i)%parent_elem = 0
  newnode_list%node(i)%ref_lambda  = 0.d0
  newnode_list%node(i)%ref_mu      = 0.d0
  newnode_list%node(i)%constrained = .false.
end do

! --- Allocate data structures for new elements and initialize them
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


! --- We find the common nodes to both patches
n_zip_nodes = 0
do i = 1,node_list1%n_nodes
  if (node_list1%node(i)%boundary .eq. 0) cycle
  if ( (xcase .gt. 0          ) .and. (i .le. 4) ) cycle ! xpoint is patched automatically
  if ( (xcase .eq. DOUBLE_NULL) .and. (i .le. 8) ) cycle
  do j=1,node_list2%n_nodes
    if (node_list2%node(j)%boundary .eq. 0) cycle
    if ( (xcase .gt. 0          ) .and. (j .le. 4) ) cycle ! xpoint is patched automatically
    if ( (xcase .eq. DOUBLE_NULL) .and. (j .le. 8) ) cycle
    R1 = node_list1%node(i)%x(1,1,1); R2 = node_list2%node(j)%x(1,1,1)
    Z1 = node_list1%node(i)%x(1,1,2); Z2 = node_list2%node(j)%x(1,1,2)
    distance = sqrt( (R2-R1)**2 + (Z2-Z1)**2 )
    if (distance .lt. tolerance) then
      n_zip_nodes = n_zip_nodes + 1
      zip_nodes(n_zip_nodes,1) = i
      zip_nodes(n_zip_nodes,2) = j
    endif
  enddo
enddo
if (debug) then
  if (xcase .eq. 0) then
    write(*,'(A,i4,A)')' Found ',n_zip_nodes,' zip nodes'
  else
    write(*,'(A,i4,A)')' Found ',n_zip_nodes+1,' zip nodes'
  endif
endif

! --- Copy 1st patch into new grid
do i=1,node_list1%n_nodes
  newnode_list%node(i)%x         = node_list1%node(i)%x          
  newnode_list%node(i)%values    = node_list1%node(i)%values     
  newnode_list%node(i)%deltas    = node_list1%node(i)%deltas     
  newnode_list%node(i)%boundary  = node_list1%node(i)%boundary   
  newnode_list%node(i)%index     = node_list1%node(i)%index
  newnode_list%node(i)%axis_node = node_list1%node(i)%axis_node
  newnode_list%node(i)%axis_dof  = node_list1%node(i)%axis_dof
enddo
newnode_list%n_nodes = node_list1%n_nodes
do i=1,element_list1%n_elements
  newelement_list%element(i)%vertex     = element_list1%element(i)%vertex      
  newelement_list%element(i)%neighbours = element_list1%element(i)%neighbours  
  newelement_list%element(i)%size       = element_list1%element(i)%size        
enddo
newelement_list%n_elements = element_list1%n_elements

! --- Then add the 2nd patch at the end
! --- First the nodes
if (xcase .gt. 0) then
  n_start = 5 ! avoid the xpoints at the beginning
  if (xcase .eq. DOUBLE_NULL) n_start = 9
else
  n_start = 1
endif
do i=1,n_start-1
  new_index(i) = i
enddo
do i=n_start,node_list2%n_nodes
  avoid_this_one = .false.
  do j=1,n_zip_nodes
    if (zip_nodes(j,2) .eq. i) then
      new_index(i) = zip_nodes(j,1)
      avoid_this_one = .true.
      exit
    endif
  enddo
  if (avoid_this_one) cycle
  newnode_list%n_nodes = newnode_list%n_nodes + 1
  new_index(i) = newnode_list%n_nodes
  newnode_list%node(newnode_list%n_nodes)%x         = node_list2%node(i)%x          
  newnode_list%node(newnode_list%n_nodes)%values    = node_list2%node(i)%values     
  newnode_list%node(newnode_list%n_nodes)%deltas    = node_list2%node(i)%deltas     
  newnode_list%node(newnode_list%n_nodes)%boundary  = node_list2%node(i)%boundary   
  newnode_list%node(newnode_list%n_nodes)%index     = node_list2%node(i)%index
  newnode_list%node(newnode_list%n_nodes)%axis_node = node_list2%node(i)%axis_node
  newnode_list%node(newnode_list%n_nodes)%axis_dof  = node_list2%node(i)%axis_dof
enddo
! --- Then the elements
do i=1,element_list2%n_elements
  newelement_list%n_elements = newelement_list%n_elements + 1
  do j=1,n_vertex_max
    k = element_list2%element(i)%vertex(j)
    newelement_list%element(newelement_list%n_elements)%vertex(j) = new_index(k)
  enddo
  newelement_list%element(newelement_list%n_elements)%neighbours  = element_list2%element(i)%neighbours  
  newelement_list%element(newelement_list%n_elements)%size        = element_list2%element(i)%size        
enddo



!----------------------------------- Print a python file that plots a cross with the 4 nodes of each element
if (plot_grid) then
  n_loop = newelement_list%n_elements
  open(101,file='plot_new_patch.py')
    write(101,'(A)')         '#!/usr/bin/env python'
    write(101,'(A)')         'import numpy as N'
    write(101,'(A)')         'import pylab'
    write(101,'(A)')         'def main():'
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
    write(101,'(A)')         '  pylab.plot(r[2*i:2*i+2],z[2*i:2*i+2], "r")'
    do j=1,n_loop
      do i=1,4
        index = newelement_list%element(j)%vertex(i)
        write(101,'(A,i6,A,f15.4)') ' r[',4*(j-1)+i-1,'] = ',newnode_list%node(index)%x(1,1,1)
        write(101,'(A,i6,A,f15.4)') ' z[',4*(j-1)+i-1,'] = ',newnode_list%node(index)%x(1,1,2)
      enddo
    enddo
    write(101,'(A,i6,A)')    ' for i in range (0,',n_loop,'):'
    write(101,'(A)')         '  pylab.plot(r[4*i:4*i+4],z[4*i:4*i+4], "b")'
    write(101,'(A)')         ' pylab.axis("equal")'
    write(101,'(A)')         ' pylab.show()'
    write(101,'(A)')         ' '
    write(101,'(A)')         'main()'
  close(101)
endif


return
end subroutine join_grid_patches



