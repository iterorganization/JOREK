subroutine log_grid_info(verbose, node_list, element_list, dir_in, filename_appendix_in)

use data_structure

implicit none

! --- Routine parameters
logical,                    intent(in)  :: verbose
type(type_node_list),       intent(in)  :: node_list
type(type_element_list),    intent(in)  :: element_list
character(len=*), optional, intent(in)  :: dir_in
character(len=*), optional, intent(in)  :: filename_appendix_in

! --- Local variables
integer, parameter :: maxbnd = 200
integer :: i, j, boundary_types(0:maxbnd), n_axis
real*8  :: Rmin, Rmax, Zmin, Zmax, x(2)
character(len=1024) :: dir, filename_appendix, filename

write(*,*)
write(*,*) 'Number of grid nodes    / Max: ', node_list%n_nodes, size(node_list%node,1)
write(*,*) 'Number of grid elements / Max: ', element_list%n_elements, size(element_list%element,1)
write(*,*)

boundary_types(:) = 0
Rmin = 1.d99
Rmax = -1.d99
Zmin = 1.d99
Zmax = -1.d99
n_axis = 0

do i = 1, node_list%n_nodes
  
  ! --- Detect how many nodes of each boundary type exist
  if ( (node_list%node(i)%boundary < 0) .or. (node_list%node(i)%boundary > maxbnd) ) then
    write(*,*) 'Skipping node ', i, ' with unusual boundary value ', node_list%node(i)%boundary
  else
    boundary_types(node_list%node(i)%boundary) = boundary_types(node_list%node(i)%boundary) + 1
  end if
  
  ! --- Determine geometrical region covered with nodes
  Rmin = min( Rmin, node_list%node(i)%x(1,1,1) )
  Rmax = max( Rmax, node_list%node(i)%x(1,1,1) )
  Zmin = min( Zmin, node_list%node(i)%x(1,1,2) )
  Zmax = max( Zmax, node_list%node(i)%x(1,1,2) )
  
  ! --- Count nodes on axis
  if ( node_list%node(i)%axis_node ) n_axis = n_axis + 1
end do

! --- Write out base information
write(*,*) 'R range of grid nodes: ', Rmin, Rmax
write(*,*) 'Z range of grid nodes: ', Zmin, Zmax
write(*,*) 'Number of nodes per boundary type:'
do i = 0, maxbnd
  if (boundary_types(i)/=0) write(*,*) 'type', i,': ',boundary_types(i),' nodes'
end do

! --- Write out nodes sorted by boundary types into ascii files (if verbose)
if ( verbose ) then
  dir               = './'
  filename_appendix = '.dat'
  if ( present(dir_in)               ) dir               = dir_in
  if ( present(filename_appendix_in) ) filename_appendix = filename_appendix_in
  write(*,*)
  write(*,*) 'Writing node coordinates sorted by boundary type into the following output files:'
  write(*,*) trim(DIR) // '/nodes_by_boundary_type_*' // trim(filename_appendix)
  do i = 0, maxbnd
    if (boundary_types(i)/=0) then
      write(filename,'(2a,i0.8,a)') trim(DIR), '/nodes_by_boundary_type_', i, trim(filename_appendix)
      open(400+i, file=filename, status='replace', form='formatted', action='write')
    end if
  end do
  do i = 1, node_list%n_nodes
      write(400+node_list%node(i)%boundary,*) node_list%node(i)%x(1,1,1:2)
  end do
  do i = 0, maxbnd
    if (boundary_types(i)/=0) close(400+i)
  end do
end if

! --- Write out all nodes
if ( verbose ) then
  write(*,*)
  write(*,*) 'Writing out all node coordinates to:'
  write(*,*) trim(DIR) // '/nodes_*' // trim(filename_appendix)
  write(filename,'(3a)') trim(DIR), '/nodes', trim(filename_appendix)
  open(400, file=filename, status='replace', form='formatted', action='write')
  do i = 1, node_list%n_nodes
      write(400,*) node_list%node(i)%x(1,1,1:2)
      write(400,*)
      write(400,*)
  end do
  close(400)
end if

! --- Write out the "center" of all elements
if ( verbose ) then
  write(*,*)
  write(*,*) 'Writing out all element center coordinates to:'
  write(*,*) trim(DIR) // '/element_centers' // trim(filename_appendix)
  write(filename,'(3a)') trim(DIR), '/element_centers', trim(filename_appendix)
  open(400, file=filename, status='replace', form='formatted', action='write')
  do i = 1, element_list%n_elements
      x = 0.d0
      do j = 1, n_vertex_max
        x = x + node_list%node(element_list%element(i)%vertex(j))%x(1,1,1:2) / real(n_vertex_max)
      end do
      write(400,*) x
      write(400,*)
      write(400,*)
  end do
  close(400)
  write(*,*)
end if

! --- Write out all nodes
if ( verbose ) then
  write(*,*)
  write(*,*) 'Writing out all node indices to:'
  write(*,*) trim(DIR) // '/nodes_*' // trim(filename_appendix)
  write(filename,'(3a)') trim(DIR), '/indices', trim(filename_appendix)
  open(400, file=filename, status='replace', form='formatted', action='write')
  do i = 1, node_list%n_nodes
      write(400,*) node_list%node(i)%index(:)
      write(400,*)
      write(400,*)
  end do
  close(400)
end if

end subroutine log_grid_info
