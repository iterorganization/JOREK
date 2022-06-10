module mod_initial_grid



implicit none



contains



!< Create a grid from parameters n_R, n_Z, n_radial, n_pol
subroutine initial_grid(node_list, element_list, bnd_node_list, bnd_elm_list, my_id, n_cpu)
  use phys_module
  use data_structure
  use mpi_mod
  use mod_boundary, only: boundary_from_grid
  use mod_element_rtree, only: populate_element_rtree
  use mod_exchange_indices
  
  implicit none
  
  type(type_node_list),         intent(inout) :: node_list
  type(type_element_list),      intent(inout) :: element_list
  type(type_bnd_node_list),     intent(inout) :: bnd_node_list
  type(type_bnd_element_list),  intent(inout) :: bnd_elm_list
  integer,                      intent(in)    :: my_id
  integer,                      intent(in)    :: n_cpu

  integer :: ierr

  element_list%n_elements      = 0
  bnd_elm_list%n_bnd_elements  = 0
  node_list%n_nodes            = 0
  
  if (my_id == 0) then
    
    ! --- Define the boundary of the initial grid
    call define_boundary()
    
    if ((n_R > 0) .and. (n_Z > 0) .and. RZ_grid_inside_wall) then
      
      call grid_inside_wall(n_R,n_Z,R_begin,R_end,Z_begin,Z_end,fbnd,node_list,element_list)
      
    else if ((n_R > 0) .and. (n_Z > 0) .and. (n_radial > 0)) then
      
      call grid_bezier_square_polar(n_R, n_Z, n_radial, R_begin, R_end, Z_begin, Z_end, R_geo,   &
        Z_geo, amin, fbnd, fpsi, mf, .true., node_list, element_list)
      
    else if ((n_R > 0) .and. (n_Z > 0) ) then
      
      call grid_bezier_square(n_R, n_Z, R_begin, R_end, Z_begin, Z_end, .true., node_list,       &
        element_list)
      
    else if ((n_radial > 0) .and. (n_pol > 0) ) then
      
      call grid_polar_bezier(R_geo, Z_geo, amin, 0.d0, 0.d0, fbnd, fpsi, mf, n_radial, n_pol,    &
        node_list, element_list)

    else
      write(*,*) ' FATAL : no valid combination of grid-sizes specified'
      call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
      stop
    end if 

    ! --- Optional: Add patches to an existing grid imported from restart file
    if ( extend_existing_grid .and. (n_flux .le. 0) ) &
        call grid_patches_on_existing_grid(node_list, element_list)

    if ( freeboundary .and. (n_flux==0) .and. freeb_change_indices ) call exchange_indices(node_list, my_id, n_cpu, .false.)
    
    ! --- Determine boundary information from the grid
    call boundary_from_grid(node_list, element_list, bnd_node_list, bnd_elm_list, .false.)
    call populate_element_rtree(node_list, element_list)

    call tr_debug_write("JMAIN:Def_grid elt_list",element_list%n_elements)
    call tr_debug_write("JMAIN:Def_grid node_list",node_list%n_nodes)
    call tr_debug_write("JMAIN:Def_grid bnd_elt_list",bnd_elm_list%n_bnd_elements)
    
  end if

end subroutine initial_grid



end module mod_initial_grid
