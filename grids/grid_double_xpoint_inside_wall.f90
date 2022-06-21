subroutine grid_double_xpoint_inside_wall(node_list, element_list)
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
use reorder_and_clean_flux_surfaces

! --- Input parameters
use phys_module, only:     n_flux, n_open, n_tht, n_outer, n_inner, n_private, n_leg, n_up_priv, n_up_leg,      &
                           n_leg_out, n_up_leg_out, n_wall_blocks,     &
                           SIG_closed, SIG_theta, SIG_open, SIG_outer, SIG_inner, SIG_private, SIG_up_priv,     &
                           SIG_theta_up, SIG_leg_0, SIG_leg_1, SIG_up_leg_0, SIG_up_leg_1,                      &
                           dPSI_open, dPSI_outer, dPSI_inner, dPSI_private, dPSI_up_priv,                       &
                           xcase, force_horizontal_Xline
use equil_info

implicit none

! --- Routine parameters
type (type_node_list),    intent(inout) :: node_list
type (type_element_list), intent(inout) :: element_list

! --- local variables
type (type_surface_list) :: flux_list, sep_list

type (type_node_list),    pointer :: node_list_inner_leg,    node_list_outer_leg
type (type_element_list), pointer :: element_list_inner_leg, element_list_outer_leg

type (type_node_list),    pointer :: node_list_tmp,    node_list_tmp2,    node_list_new
type (type_element_list), pointer :: element_list_tmp, element_list_tmp2, element_list_new

type (type_strategic_points) , pointer     :: stpts
type (type_new_points)       , pointer     :: nwpts

real*8              :: psi_xpoint(2)
integer             :: n_psi
integer             :: i_elm_find(8), ifail
integer             :: i_ext
real*8              :: psi_bnd, psi_bnd2
real*8              :: sigmas(17)
integer             :: n_grids(12)
integer             :: n_seg_prev
real*8              :: seg_prev(n_seg_max)
integer             :: n_loop, i, j, index
logical, parameter  :: plot_grid = .true.
character*2         :: char_patch
character*256       :: filename

write(*,*) ' '
write(*,*) ' '
write(*,*) '*************************************'
write(*,*) '*************************************'
write(*,*) '*************************************'
write(*,*) '*          X-point grid             *'
write(*,*) '*************************************'
write(*,*) '*************************************'
write(*,*) '*************************************'
write(*,*) ' '
write(*,*) ' '


allocate(stpts)
allocate(nwpts)
call tr_register_mem(sizeof(stpts),"stpts",CAT_GRID)
call tr_register_mem(sizeof(nwpts),"nwpts",CAT_GRID)
!-------------------------------------------------------------------------------------------!
!---------------------------- Initialise some internal data --------------------------------!
!-------------------------------------------------------------------------------------------!









!-------------------------------- Redefine neighbours
if (sum(element_list%element(1)%neighbours) .eq. 0) then
  call update_neighbours_basic(element_list,node_list)
endif

!-------------------------------- Reset some parameters if they are inconsistent with XCASE
if (xcase .eq. LOWER_XPOINT) then
  n_outer   = 0
  n_inner   = 0
  n_up_priv = 0
  n_up_leg  = 0
endif
if (xcase .eq. UPPER_XPOINT) then
  n_outer   = 0
  n_inner   = 0
  n_private = 0
  n_leg     = 0
endif
if ( (xcase .eq. DOUBLE_NULL) .and. (mod(n_tht,2) .ne. 0) )  n_tht = n_tht + 1
if ( (xcase .ne. DOUBLE_NULL) .and. (mod(n_tht,2) .eq. 0) )  n_tht = n_tht + 1

!-------------------------------- Check consistency of grid inputs
if (xcase .eq. DOUBLE_NULL) then
  if ( (n_outer .le. 1) .or. (n_inner .le. 1) .or. (n_private .le. 1) .or. (n_up_priv .le. 1) ) then
    write(*,'(A)')'Unfortunately, there are minimal requirements for some of the grid inputs.'
    write(*,'(A)')'These are as follows:'
    write(*,'(A)')'n_outer   > 1'
    write(*,'(A)')'n_inner   > 1'
    write(*,'(A)')'n_private > 1'
    write(*,'(A)')'n_up_priv > 1'
    write(*,'(A)')'Please changes accordingly. Aborting...'
    stop
  endif
endif

!-------------------------------- Build up some arrays to send as routine parameters (avoid long lists...)
sigmas(1)  = SIG_closed  ; sigmas(2)  = SIG_theta
sigmas(3)  = SIG_open    ; sigmas(4)  = SIG_outer   ; sigmas(5)  = SIG_inner
sigmas(6)  = SIG_private ; sigmas(7)  = SIG_up_priv
sigmas(8)  = SIG_leg_0   ; sigmas(9)  = SIG_leg_1
sigmas(10) = SIG_up_leg_0; sigmas(11) = SIG_up_leg_1
sigmas(12) = dPSI_open   ; sigmas(13) = dPSI_outer  ; sigmas(14) = dPSI_inner
sigmas(15) = dPSI_private; sigmas(16) = dPSI_up_priv
if ( SIG_theta_up .eq. 999.d0 ) then
  sigmas(17) = SIG_theta
else
  sigmas(17) = SIG_theta_up
endif

n_grids(1) = n_flux   ; n_grids(2) = n_tht
n_grids(3) = n_open   ; n_grids(4) = n_outer  ; n_grids(5) = n_inner
n_grids(6) = n_private; n_grids(7) = n_up_priv
n_grids(8) = n_leg    ; n_grids(9) = n_up_leg
n_grids(10)= n_leg_out; n_grids(11)= n_up_leg_out
!n_grids(10)= 0 ! keep for n_tht_outer, which determines the angle of second Xpoint

!-------------------------------------------------------------------------------------------!
!----------------------------- Find MagAxis and Xpoint -------------------------------------!
!-------------------------------------------------------------------------------------------!
psi_xpoint(1) = ES%psi_xpoint(1)
psi_xpoint(2) = ES%psi_xpoint(2)

psi_bnd  = 0.d0
psi_bnd2 = 0.d0
if(xcase .eq. LOWER_XPOINT) psi_bnd = psi_xpoint(1)
if(xcase .eq. UPPER_XPOINT) psi_bnd = psi_xpoint(2)
if(xcase .eq. DOUBLE_NULL ) then
  if(ES%active_xpoint .eq. UPPER_XPOINT) then
    psi_bnd  = psi_xpoint(2)
    psi_bnd2 = psi_xpoint(1)
  else
    psi_bnd  = psi_xpoint(1)
    psi_bnd2 = psi_xpoint(2)  
  endif
  ! If we have a symmetric double-null, force the single separatrix
  if (ES%active_xpoint .eq. SYMMETRIC_XPOINT) then
    psi_xpoint(1)  = (psi_xpoint(1)+psi_xpoint(2))/2.d0
    psi_xpoint(2)  = psi_xpoint(1)
    psi_bnd  = psi_xpoint(1)
    psi_bnd2 = psi_bnd  
    n_grids(3) = 0
  endif
endif



!-------------------------------------------------------------------------------------------!
!--------------- Define the flux values on which grid will be aligned ----------------------!
!-------------------------------------------------------------------------------------------!

!-------------------------------- Define number of psi values and allocate flux_list structure
n_psi           = n_flux   + n_open   + n_outer   + n_inner   + n_private   + n_up_priv + 1   ! this includes the magnetic axis
flux_list%n_psi = n_psi - 1
call tr_allocate(flux_list%psi_values,1,flux_list%n_psi,"flux_list%psi_values",CAT_GRID)

!-------------------------------- Allocate sep_list structure (that's for plotting only)
sep_list%n_psi =3
if(xcase .eq. DOUBLE_NULL) sep_list%n_psi =6
if (allocated(sep_list%psi_values)) call tr_deallocate(sep_list%psi_values,"sep_list%psi_values",CAT_GRID)
call tr_allocate(sep_list%psi_values,1,sep_list%n_psi,"sep_list%psi_values",CAT_GRID)

!-------------------------------- Call the routine
call define_flux_values(node_list, element_list, flux_list, sep_list, xcase, psi_xpoint, n_grids, sigmas)

call plot_flux_surfaces(node_list,element_list,flux_list,.true.,1,.true.,xcase)
call plot_flux_surfaces(node_list,element_list,sep_list,.false.,1,.true.,xcase)

if (allocated(sep_list%flux_surfaces))     deallocate(sep_list%flux_surfaces)



!-------------------------------------------------------------------------------------------!
!-------- Find all strategic points (leg corners, strike points and private middles) -------!
!-------------------------------------------------------------------------------------------!

!-------------------------------- Call the routine
call reorder_flux_surfaces(node_list, element_list, flux_list, .true., ifail)
call clean_surfaces(node_list,element_list,flux_list,n_grids)
call find_strategic_points_advanced(node_list, element_list, flux_list, xcase, force_horizontal_Xline, n_grids, stpts)


!-------------------------------------------------------------------------------------------!
!-------------- Find new grid points by crossing polar and radial coordinates --------------!
!-------------------------------------------------------------------------------------------!

! --- Allocate data structures for new nodes and initialize them
allocate(node_list_inner_leg,node_list_outer_leg,node_list_tmp,node_list_tmp2,node_list_new)
call tr_register_mem(sizeof(node_list_inner_leg),"node_list_inner_leg")
call tr_register_mem(sizeof(node_list_outer_leg),"node_list_outer_leg")
call tr_register_mem(sizeof(node_list_tmp),"node_list_tmp")
call tr_register_mem(sizeof(node_list_tmp2),"node_list_tmp2")
call tr_register_mem(sizeof(node_list_new),"node_list_new")

! --- Allocate data structures for new elements and initialize them
allocate(element_list_inner_leg,element_list_outer_leg,element_list_tmp,element_list_tmp2,element_list_new)
call tr_register_mem(sizeof(element_list_inner_leg),"element_list_inner_leg")
call tr_register_mem(sizeof(element_list_outer_leg),"element_list_outer_leg")
call tr_register_mem(sizeof(element_list_tmp),"element_list_tmp")
call tr_register_mem(sizeof(element_list_tmp2),"element_list_tmp2")
call tr_register_mem(sizeof(element_list_new),"element_list_new")


!-------------------------------- Respline flux surfaces
!call respline_flux_surfaces(node_list,element_list,flux_list)



!-------------------------------- First the Central part
call define_central_grid(node_list, element_list, flux_list, &
                         xcase, n_grids, stpts, sigmas, nwpts, node_list_new,  element_list_new)
call update_neighbours_basic(element_list_new,node_list_new)
call update_boundary_types  (element_list_new,node_list_new, 1)

!-------------------------------- Then the lower legs
if (xcase .ne. UPPER_XPOINT) then
  !-------------------------------- The lower inner leg
  call define_leg_grid(node_list, element_list, node_list_inner_leg, element_list_inner_leg, &
                       flux_list, xcase, n_grids, stpts, sigmas, nwpts, 1)
  call update_neighbours_basic(element_list_inner_leg,node_list_inner_leg)
  call update_boundary_types  (element_list_inner_leg,node_list_inner_leg, 2)

  !-------------------------------- The lower outer leg
  call define_leg_grid(node_list, element_list, node_list_outer_leg, element_list_outer_leg, &
                       flux_list, xcase, n_grids, stpts, sigmas, nwpts, 2)
  call update_neighbours_basic(element_list_outer_leg,node_list_outer_leg)
  call update_boundary_types  (element_list_outer_leg,node_list_outer_leg, 2)

  !-------------------------------- Join the lower legs together
  call join_grid_patches(node_list_inner_leg, element_list_inner_leg, &
                         node_list_outer_leg, element_list_outer_leg, &
                         node_list_tmp,       element_list_tmp, xcase)
  call update_neighbours_basic(element_list_tmp,node_list_tmp)
  call update_boundary_types  (element_list_tmp,node_list_tmp, 1)

  !-------------------------------- Join lower legs with main part
  call join_grid_patches(node_list_new,  element_list_new,  &
                         node_list_tmp,  element_list_tmp,  &
                         node_list_tmp2, element_list_tmp2, xcase)
  call update_neighbours_basic(element_list_tmp2,node_list_tmp2)
  call update_boundary_types  (element_list_tmp2,node_list_tmp2, 1)
  call copy_node_structure(node_list_new, element_list_new, node_list_tmp2, element_list_tmp2)
endif

!-------------------------------- Then the Upper legs
if (xcase .ne. LOWER_XPOINT) then
  !-------------------------------- The Upper inner leg
  call define_leg_grid(node_list, element_list, node_list_inner_leg, element_list_inner_leg, &
                       flux_list, xcase, n_grids, stpts, sigmas, nwpts, 3)
  call update_neighbours_basic(element_list_inner_leg,node_list_inner_leg)
  call update_boundary_types  (element_list_inner_leg,node_list_inner_leg, 2)
  
  !-------------------------------- The Upper outer leg
  call define_leg_grid(node_list, element_list, node_list_outer_leg, element_list_outer_leg, &
                       flux_list, xcase, n_grids, stpts, sigmas, nwpts, 4)
  call update_neighbours_basic(element_list_outer_leg,node_list_outer_leg)
  call update_boundary_types  (element_list_outer_leg,node_list_outer_leg, 2)
  
  !-------------------------------- Join the upper legs together
  call join_grid_patches(node_list_inner_leg, element_list_inner_leg, &
                         node_list_outer_leg, element_list_outer_leg, &
                         node_list_tmp,       element_list_tmp, xcase)
  call update_neighbours_basic(element_list_tmp,node_list_tmp)
  call update_boundary_types  (element_list_tmp,node_list_tmp, 1)
  
  !-------------------------------- Join upper legs with rest of the grid
  call join_grid_patches(node_list_new,  element_list_new, &
                         node_list_tmp,  element_list_tmp, &
                         node_list_tmp2, element_list_tmp2, xcase)
  call copy_node_structure(node_list_new, element_list_new, node_list_tmp2, element_list_tmp2)
  call update_neighbours_basic(element_list_new,node_list_new)
  call update_boundary_types  (element_list_new,node_list_new, 1)
endif
call temporary_element_sizes(node_list_new, element_list_new)
call export_restart(node_list_new, element_list_new, 'grid_no_patch')

!-------------------------------- Now the wall extension
do i_ext = 1,n_wall_blocks
  call define_extension_patch(node_list_new, element_list_new, node_list_tmp, element_list_tmp, n_seg_prev, seg_prev, i_ext)
  call update_neighbours_basic(element_list_tmp,node_list_tmp)
  call update_boundary_types  (element_list_tmp,node_list_tmp, 0)
  ! --- create restart file for vtk plots BEG
  if (i_ext .lt. 10) write(char_patch,'(i1)') i_ext
  if (i_ext .ge. 10) write(char_patch,'(i2)') i_ext
  write(filename,'(A10,A)')'grid_patch',trim(char_patch)
  call temporary_element_sizes(node_list_tmp, element_list_tmp)
  call export_restart(node_list_tmp, element_list_tmp, filename)
  ! --- create restart file for vtk plots END
  call join_grid_patches(node_list_new,  element_list_new, &
                         node_list_tmp,  element_list_tmp, &
                         node_list_tmp2, element_list_tmp2, xcase)
  call copy_node_structure(node_list_new, element_list_new, node_list_tmp2, element_list_tmp2)
  call update_neighbours_basic(element_list_new,node_list_new)
  call update_boundary_types  (element_list_new,node_list_new, 0)
enddo



!-------------------------------- Finalise grid (element size, nodes index etc.)
call finish_grid(node_list, element_list, node_list_new, element_list_new, n_grids, .true., .true., .true.)
call export_restart(node_list, element_list, 'jorek_grid')




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
deallocate(node_list_inner_leg,node_list_outer_leg,node_list_tmp,node_list_tmp2,node_list_new)
call tr_unregister_mem(sizeof(node_list_inner_leg),"node_list_inner_leg")
call tr_unregister_mem(sizeof(node_list_outer_leg),"node_list_outer_leg")
call tr_unregister_mem(sizeof(node_list_tmp),"node_list_tmp")
call tr_unregister_mem(sizeof(node_list_tmp2),"node_list_tmp2")
call tr_unregister_mem(sizeof(node_list_new),"node_list_new")

! --- Deallocate data structures for new elements and initialize them
deallocate(element_list_inner_leg,element_list_outer_leg,element_list_tmp,element_list_tmp2,element_list_new)
call tr_register_mem(sizeof(element_list_inner_leg),"element_list_inner_leg")
call tr_register_mem(sizeof(element_list_outer_leg),"element_list_outer_leg")
call tr_register_mem(sizeof(element_list_tmp),"element_list_tmp")
call tr_register_mem(sizeof(element_list_tmp2),"element_list_tmp2")
call tr_register_mem(sizeof(element_list_new),"element_list_new")



!-------------------------------- Deallocate data
deallocate(stpts)
deallocate(nwpts)
call tr_unregister_mem(sizeof(stpts),"stpts",CAT_GRID)
call tr_unregister_mem(sizeof(nwpts),"nwpts",CAT_GRID)

return
end subroutine grid_double_xpoint_inside_wall










subroutine copy_node_structure(node_list, element_list, newnode_list, newelement_list)

  use data_structure
  
  implicit none
  
  ! --- Routine variables
  type (type_node_list)       , intent(inout) :: node_list
  type (type_element_list)    , intent(inout) :: element_list
  type (type_node_list)       , intent(in)    :: newnode_list
  type (type_element_list)    , intent(in)    :: newelement_list
  
  ! --- copy new grid into nodes/elements
  node_list%n_nodes = newnode_list%n_nodes
  node_list%node(1:node_list%n_nodes) = newnode_list%node(1:node_list%n_nodes)
  
  element_list%n_elements = newelement_list%n_elements
  element_list%element(1:element_list%n_elements) = newelement_list%element(1:element_list%n_elements)
  
  return

end subroutine copy_node_structure


