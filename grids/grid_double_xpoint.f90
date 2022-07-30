subroutine grid_double_xpoint(node_list, element_list)
!-----------------------------------------------------------------------
! subroutine defines a flux surface aligned finite element grid
! inclduing a single x-point
!-----------------------------------------------------------------------

use constants
use tr_module 
use data_structure
use grid_xpoint_data
use mod_neighbours, only: update_neighbours

! --- Input parameters
use phys_module, only:     n_flux, n_open, n_tht, n_outer, n_inner, n_private, n_leg, n_up_priv, n_up_leg,      &
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

type (type_strategic_points) , pointer     :: stpts
type (type_new_points)       , pointer     :: nwpts

real*8              :: psi_xpoint(2)
integer             :: n_psi
integer             :: i_elm_find(8), ifail
real*8              :: psi_bnd, psi_bnd2
real*8              :: sigmas(17)
integer             :: n_grids(12)

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

!-------------------------------- Build up some arrays to send as routine parameters (avoid long lists...)
sigmas(1)  = SIG_closed  ; sigmas(2)  = SIG_theta
sigmas(3)  = SIG_open    ; sigmas(4)  = SIG_outer   ; sigmas(5)  = SIG_inner
sigmas(6)  = SIG_private ; sigmas(7)  = SIG_up_priv
sigmas(8)  = SIG_leg_0   ; sigmas(9)  = SIG_leg_1
sigmas(10) = SIG_up_leg_0; sigmas(11) = SIG_up_leg_1
sigmas(12) = dPSI_open   ; sigmas(13) = dPSI_outer  ; sigmas(14) = dPSI_inner
sigmas(15) = dPSI_private; sigmas(16) = dPSI_up_priv
if ( SIG_theta_up .eq. 999.d0 ) then
  sigmas(17) = SIG_theta ! if SIG_theta_up is not specified it is set to SIG_theta
else
  sigmas(17) = SIG_theta_up
endif

n_grids(1) = n_flux   ; n_grids(2) = n_tht
n_grids(3) = n_open   ; n_grids(4) = n_outer  ; n_grids(5) = n_inner
n_grids(6) = n_private; n_grids(7) = n_up_priv
n_grids(8) = n_leg    ; n_grids(9) = n_up_leg
n_grids(10)= 0 ! keep for n_tht_outer, which determines the angle of second Xpoint

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
    psi_xpoint(1) = (psi_xpoint(1)+psi_xpoint(2))/2.d0
    psi_xpoint(2) = psi_xpoint(1)
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
call find_strategic_points(node_list, element_list, flux_list, xcase, force_horizontal_Xline, n_grids, stpts)


!-------------------------------------------------------------------------------------------!
!-------------- Find new grid points by crossing polar and radial coordinates --------------!
!-------------------------------------------------------------------------------------------!

!-------------------------------- Call the routine
call define_new_grid_points(node_list, element_list, flux_list, xcase, n_grids, stpts, sigmas, nwpts)



!-------------------------------------------------------------------------------------------!
!---------------------- Define the final grid (new nodes and new elements) -----------------!
!-------------------------------------------------------------------------------------------!

!-------------------------------- Call the routine
call define_final_grid(node_list, element_list, flux_list, &
                       xcase, n_grids, stpts, nwpts)


!-------------------------------- Deallocate data
deallocate(stpts)
deallocate(nwpts)
call tr_unregister_mem(sizeof(stpts),"stpts",CAT_GRID)
call tr_unregister_mem(sizeof(nwpts),"nwpts",CAT_GRID)

call update_neighbours(node_list,element_list, force_rtree_initialize=.true.)
return
end subroutine grid_double_xpoint
