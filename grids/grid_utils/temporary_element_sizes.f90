!> Subroutine defines basic element sizes until detailed sizes are set
subroutine temporary_element_sizes(newnode_list, newelement_list)

use constants
use tr_module 
use data_structure
use grid_xpoint_data

implicit none

! --- Routine parameters
type (type_node_list)       , intent(inout) :: newnode_list
type (type_element_list)    , intent(inout) :: newelement_list

! --- local variables
integer :: k, iv, ivp, node_iv, node_ivp
real*8  :: R0, Z0, dR0, dZ0
real*8  :: RP, ZP, dRP, dZP
real*8  :: size_0, size_p, denom


! --- Temporary, just for the vtk plots, will be overwritten in "finish_grid.f90"
do k=1, newelement_list%n_elements   ! fill in the size of the elements
  do iv = 1, 4                    ! over 4 sides of an element

    ivp = mod(iv,4)   + 1         ! vertex with index one higher
    node_iv  = newelement_list%element(k)%vertex(iv)
    node_ivp = newelement_list%element(k)%vertex(ivp) 

    if ((iv .eq. 1) .or. (iv .eq. 3)) then
      R0 = newnode_list%node(node_iv )%X(1,1,1)  ; dR0 = newnode_list%node(node_iv )%X(1,2,1)
      Z0 = newnode_list%node(node_iv )%X(1,1,2)  ; dZ0 = newnode_list%node(node_iv )%X(1,2,2)
      RP = newnode_list%node(node_ivp)%X(1,1,1)  ; dRP = newnode_list%node(node_ivp)%X(1,2,1)
      ZP = newnode_list%node(node_ivp)%X(1,1,2)  ; dZP = newnode_list%node(node_ivp)%X(1,2,2)
    else
      R0 = newnode_list%node(node_iv )%X(1,1,1)  ; dR0 = newnode_list%node(node_iv )%X(1,3,1)
      Z0 = newnode_list%node(node_iv )%X(1,1,2)  ; dZ0 = newnode_list%node(node_iv )%X(1,3,2)
      RP = newnode_list%node(node_ivp)%X(1,1,1)  ; dRP = newnode_list%node(node_ivp)%X(1,3,1)
      ZP = newnode_list%node(node_ivp)%X(1,1,2)  ; dZP = newnode_list%node(node_ivp)%X(1,3,2)
    endif

    size_0 = 1.d0
    size_p = 1.d0
    denom = ( dRP * dZ0 - dR0 * dZP)
    size_0 = sign(sqrt((R0-RP)**2 + (Z0-ZP)**2) /3.d0, dR0 * (RP-R0) + dZ0 * (ZP-Z0) )
    size_P = sign(sqrt((R0-RP)**2 + (Z0-ZP)**2) /3.d0, dRP * (R0-RP) + dZP * (Z0-ZP) )

    if ((R0-RP)**2 + (Z0-ZP)**2 .eq. 0.d0) then
      size_0 = 1.d0
      size_P = 1.d0
    endif

    if ((iv .eq. 1) .or. (iv .eq. 3)) then
      newelement_list%element(k)%size(iv,2)  = size_0
      newelement_list%element(k)%size(ivp,2) = size_p
    else
      newelement_list%element(k)%size(iv,3)  = size_0
      newelement_list%element(k)%size(ivp,3) = size_p
    endif

  enddo

  do iv=1,4
    newelement_list%element(k)%size(iv,1) = 1.d0
    newelement_list%element(k)%size(iv,4) = newelement_list%element(k)%size(iv,2) * newelement_list%element(k)%size(iv,3)
  enddo

enddo


return
end subroutine temporary_element_sizes





