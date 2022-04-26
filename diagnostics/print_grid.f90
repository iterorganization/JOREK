subroutine print_grid(node_list,element_list,boundary_list)
!----------------------------------------------------------------
! plot the grid of finite elements with the correct curved edges
!----------------------------------------------------------------
use data_structure
use mod_parameters, only : n_order

implicit none

type (type_node_list)        :: node_list
type (type_element_list)     :: element_list
type (type_bnd_element_list) :: boundary_list
integer :: i

write(*,*) '**************************************************'
write(*,*) '*           finite element grid                  *'
write(*,*) '**************************************************'

if (n_order .gt. 3) then
  write(*,*)'WARNING:'
  write(*,*)'This routine is not addapted to n_order>5'
  write(*,*)'Please use export_restart instead.'
  return
endif

write(*,*) ' node_list : n_nodes =',node_list%n_nodes
write(*,'(A68)') '   i,    x,      y,      ux,     uy,     vx,     vy,     wx,     wy'
do i=1,node_list%n_nodes
  write(*,'(i5,8f8.3)') i,node_list%node(i)%x(1,1,1:2),node_list%node(i)%x(1,2,1:2), &
                         node_list%node(i)%x(1,3,1:2),node_list%node(i)%x(1,4,1:2)
enddo
write(*,*)
write(*,*) ' element_list : n_elements =',element_list%n_elements
write(*,'(A88)') '    i,  iv1, iv2, iv3, iv4,  hu1,    hv1,    hu2,    hv2,    hu3,    hv3,    hu4,    hv4,'
do i=1,element_list%n_elements
  write(*,'(5i5,8f12.5)') i,element_list%element(i)%vertex(1:4)   , &
                          element_list%element(i)%size(1,2:3), &
                          element_list%element(i)%size(2,2:3), &
                          element_list%element(i)%size(3,2:3), &
                          element_list%element(i)%size(4,2:3)
enddo

write(*,*)
write(*,*) ' boundary_list : n_boundary=',boundary_list%n_bnd_elements
write(*,'(A88)') '    i,  iv1, iv2, idir1, idir1, idir2, idir2, elm, side, h1,  hv1,  h2,  hv2'
do i=1,boundary_list%n_bnd_elements
  write(*,'(9i5,4f12.5)') i,boundary_list%bnd_element(i)%vertex(1:2)   , &
                          boundary_list%bnd_element(i)%direction(1,1:2), &
                          boundary_list%bnd_element(i)%direction(2,1:2), &
                          boundary_list%bnd_element(i)%element, &
                          boundary_list%bnd_element(i)%side, &
                          boundary_list%bnd_element(i)%size(1,1:2), boundary_list%bnd_element(i)%size(2,1:2) 
enddo

write(*,*)
write(*,*) ' node_list : index '
do i=1,node_list%n_nodes
  write(*,'(5i5,4e16.8)') i,node_list%node(i)%index(:),node_list%node(i)%values(1,:,6)
enddo
write(*,*) ' done print_grid'
return
end
