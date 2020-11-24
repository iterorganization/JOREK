!< Create a grid from parameters n_R, n_Z, n_radial, n_pol
subroutine initial_grid(node_list,element_list,bnd_node_list,bnd_elm_list)
  use phys_module
  use data_structure
  use mpi_mod
  use mod_boundary, only: boundary_from_grid
  use mod_element_rtree, only: populate_element_rtree
  implicit none
  type(type_node_list), intent(inout) :: node_list
  type(type_element_list), intent(inout) :: element_list
  type(type_bnd_node_list), intent(inout) :: bnd_node_list
  type(type_bnd_element_list), intent(inout) :: bnd_elm_list
  integer :: ierr

  element_list%n_elements      = 0
  bnd_elm_list%n_bnd_elements  = 0
  bnd_node_list%n_bnd_nodes    = 0
  node_list%n_nodes            = 0
  ! --- Define the boundary of the initial grid
  call define_boundary()
  
  if ((n_R > 0) .and. (n_Z > 0) .and. (n_radial > 0)) then
    
    call grid_bezier_square_polar(n_R, n_Z, n_radial, R_begin, R_end, Z_begin, Z_end, R_geo,   &
                                  Z_geo, amin, fbnd, fpsi, mf, .true., node_list, element_list)

  else if ((n_R > 0) .and. (n_Z > 0) ) then
    
    call grid_bezier_square(n_R, n_Z, R_begin, R_end, Z_begin, Z_end, .true., node_list, element_list)
    
  else if ((n_radial > 0) .and. (n_pol > 0) ) then
    
    call grid_polar_bezier(R_geo, Z_geo, amin, 0.d0, 0.d0, fbnd, fpsi, mf, n_radial, n_pol, node_list, element_list)
    
  else
    write(*,*) ' FATAL : no valid combination of grid-sizes specified'
    call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    stop
  end if 
  
  if ( freeboundary .and. freeb_change_indices ) call exchange_indices_for_vacuum(node_list, 0) ! typically call this
  ! routine from my_id 0 only
  
  ! --- Determine boundary information from the grid
  call boundary_from_grid(node_list, element_list, bnd_node_list, bnd_elm_list, .false.)

  call populate_element_rtree(node_list, element_list)

end subroutine initial_grid
