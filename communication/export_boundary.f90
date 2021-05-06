!> Export boundary information for STARWALL
subroutine export_boundary(node_list, boundary_list, bnd_node_list)

  use data_structure 
  use mod_parameters, only: n_coord_tor, n_coord_period

  implicit none

  ! --- Routine parameters.
  type(type_node_list),        intent(in) :: node_list     !< List of grid nodes.
  type(type_bnd_element_list), intent(in) :: boundary_list !< List of boundary elements.
  type(type_bnd_node_list),    intent(in) :: bnd_node_list !< List of boundary nodes.

  ! --- Local variables.
  integer :: i, iv1, iv2, ib1, ib2, idir1, idir2
  integer :: bnd_file_version = 2
  open(22,file='boundary.txt')

  write(22,'(5i7)') boundary_list%n_bnd_elements, bnd_node_list%n_bnd_nodes, n_coord_tor, n_coord_period, bnd_file_version

  do i=1,boundary_list%n_bnd_elements

    iv1   = boundary_list%bnd_element(i)%vertex(1)
    iv2   = boundary_list%bnd_element(i)%vertex(2)
    ib1   = boundary_list%bnd_element(i)%bnd_vertex(1)
    ib2   = boundary_list%bnd_element(i)%bnd_vertex(2)
    idir1 = boundary_list%bnd_element(i)%direction(1,2)
    idir2 = boundary_list%bnd_element(i)%direction(2,2)

    write(22,'(3i7,999e16.8)') i, ib1, ib2,      &
      node_list%node(iv1)%x(:,1,1),             &
      node_list%node(iv1)%x(:,1,2),             &
      node_list%node(iv1)%x(:,idir1,1),         &
      node_list%node(iv1)%x(:,idir1,2),         &
      boundary_list%bnd_element(i)%size(1,1:2), &
      node_list%node(iv2)%x(:,1,1),             &
      node_list%node(iv2)%x(:,1,2),             &
      node_list%node(iv2)%x(:,idir1,1),         &
      node_list%node(iv2)%x(:,idir1,2),         &
      boundary_list%bnd_element(i)%size(2,1:2) 

  end do

  do i=1, bnd_node_list%n_bnd_nodes
    write(22,'(4i7)') i,bnd_node_list%bnd_node(i)%index_starwall,bnd_node_list%bnd_node(i)%n_dof
  enddo

  close(22)
  
  return
end subroutine export_boundary
