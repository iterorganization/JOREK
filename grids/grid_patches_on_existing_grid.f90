subroutine grid_patches_on_existing_grid(node_list, element_list)
!-----------------------------------------------------------------------
! subroutine defines a flux surface aligned finite element grid
! inclduing a single x-point
!-----------------------------------------------------------------------

use constants
use tr_module 
use data_structure
use grid_xpoint_data
use mod_export_restart
use mod_eqdsk_tools
use mod_element_rtree

! --- Input parameters
use phys_module, only:     n_wall_blocks, xpoint

implicit none

! --- Routine parameters
type (type_node_list),    intent(inout) :: node_list
type (type_element_list), intent(inout) :: element_list

! --- local variables
type (type_surface_list) :: flux_list, sep_list

type (type_node_list),    pointer :: node_list_tmp,    node_list_tmp2,    node_list_new
type (type_element_list), pointer :: element_list_tmp, element_list_tmp2, element_list_new

integer             :: my_id, i_ext
integer             :: n_grids(12)
integer             :: n_seg_prev
real*8              :: seg_prev(n_seg_max)
integer             :: n_loop, i, j, k, index
logical, parameter  :: plot_grid = .true.

my_id  = 1 ! Just don't want the printout...

write(*,*) ' '
write(*,*) ' '
write(*,*) '*************************************'
write(*,*) '*************************************'
write(*,*) '*************************************'
write(*,*) '* Adding extension patches to grid  *'
write(*,*) '*************************************'
write(*,*) '*************************************'
write(*,*) '*************************************'
write(*,*) ' '
write(*,*) ' '


!-------------------------------------------------------------------------------------------!
!---------------------------- Initialise some internal data --------------------------------!
!-------------------------------------------------------------------------------------------!








!-------------------------------- Redefine neighbours
if (sum(element_list%element(1)%neighbours) .eq. 0) then
  call update_neighbours_basic(element_list,node_list)
endif

! --- Allocate data structures for new nodes and initialize them
allocate(node_list_tmp,node_list_tmp2,node_list_new)
call tr_register_mem(sizeof(node_list_tmp),"node_list_tmp")
call tr_register_mem(sizeof(node_list_tmp2),"node_list_tmp2")
call tr_register_mem(sizeof(node_list_new),"node_list_new")
! --- Allocate data structures for new elements and initialize them
allocate(element_list_tmp,element_list_tmp2,element_list_new)
call tr_register_mem(sizeof(element_list_tmp),"element_list_tmp")
call tr_register_mem(sizeof(element_list_tmp2),"element_list_tmp2")
call tr_register_mem(sizeof(element_list_new),"element_list_new")
! --- Copy structure into work structure
call copy_node_structure(node_list_new, element_list_new, node_list, element_list)

!-------------------------------- Now the wall extension
do i_ext = 1,n_wall_blocks
  call define_extension_patch(node_list_new, element_list_new, node_list_tmp, element_list_tmp, n_seg_prev, seg_prev, i_ext)
  call update_neighbours_basic(element_list_tmp,node_list_tmp)
  call update_boundary_types  (element_list_tmp,node_list_tmp, 0)
  call join_grid_patches(node_list_new,  element_list_new, &
                         node_list_tmp,  element_list_tmp, &
                         node_list_tmp2, element_list_tmp2, .false.)
  call copy_node_structure(node_list_new, element_list_new, node_list_tmp2, element_list_tmp2)
  call update_neighbours_basic(element_list_new,node_list_new)
  call update_boundary_types  (element_list_new,node_list_new, 0)
enddo



! --- Finalise grid (element size, nodes index etc.)
if (xpoint) then
  call finish_grid(node_list, element_list, node_list_new, element_list_new, n_grids)
else
  call finish_grid_no_xpoint_no_axis(node_list, element_list, node_list_new, element_list_new, n_grids)
endif

call export_restart(node_list, element_list, 'jorek_restart')




!----------------------------------- Print a python file that plots a cross with the 4 nodes of each element
if (plot_grid) then
  n_loop = element_list%n_elements
  open(101,file='plot_elements.py')
    write(101,'(A)')         '#!/usr/bin/env python'
    write(101,'(A)')         'import numpy as N'
    write(101,'(A)')         'import pylab'
    write(101,'(A)')         'def main():'
    write(101,'(A,i6,A)')    ' r = N.zeros(',4*n_loop,')'
    write(101,'(A,i6,A)')    ' z = N.zeros(',4*n_loop,')'
    do j=1,n_loop
      do i=1,2
        index = element_list%element(j)%vertex(i)
        write(101,'(A,i6,A,f15.4)') ' r[',4*(j-1)+2*i-2,'] = ',node_list%node(index)%x(1,1,1)
        write(101,'(A,i6,A,f15.4)') ' z[',4*(j-1)+2*i-2,'] = ',node_list%node(index)%x(1,1,2)
        index = element_list%element(j)%vertex(i+2)
        write(101,'(A,i6,A,f15.4)') ' r[',4*(j-1)+2*i-1,'] = ',node_list%node(index)%x(1,1,1)
        write(101,'(A,i6,A,f15.4)') ' z[',4*(j-1)+2*i-1,'] = ',node_list%node(index)%x(1,1,2)
      enddo
    enddo
    write(101,'(A,i6,A)')    ' for i in range (0,',n_loop,'):'
    write(101,'(A)')         '  pylab.plot(r[4*i:4*i+2],z[4*i:4*i+2], "r")'
    write(101,'(A)')         '  pylab.plot(r[4*i+2:4*i+4],z[4*i+2:4*i+4], "g")'
    do j=1,n_loop
      do i=1,4
        index = element_list%element(j)%vertex(i)
        write(101,'(A,i6,A,f15.4)') ' r[',4*(j-1)+i-1,'] = ',node_list%node(index)%x(1,1,1)
        write(101,'(A,i6,A,f15.4)') ' z[',4*(j-1)+i-1,'] = ',node_list%node(index)%x(1,1,2)
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




! --- Deallocate data structures for new nodes and initialize them
deallocate(node_list_tmp,node_list_tmp2,node_list_new)
call tr_unregister_mem(sizeof(node_list_tmp),"node_list_tmp")
call tr_unregister_mem(sizeof(node_list_tmp2),"node_list_tmp2")
call tr_unregister_mem(sizeof(node_list_new),"node_list_new")

! --- Deallocate data structures for new elements and initialize them
deallocate(element_list_tmp,element_list_tmp2,element_list_new)
call tr_register_mem(sizeof(element_list_tmp),"element_list_tmp")
call tr_register_mem(sizeof(element_list_tmp2),"element_list_tmp2")
call tr_register_mem(sizeof(element_list_new),"element_list_new")



return
end subroutine grid_patches_on_existing_grid







