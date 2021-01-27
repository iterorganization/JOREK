subroutine create_x_node(node_list, element_list, newnode_list, nwpts, stpts, &
                         Xup, R_axis, Z_axis, R_xpoint, Z_xpoint, i_elm_xpoint, s_xpoint, t_xpoint)
!------------------------------------------------------------------------------------------
! subroutine defines a new node with the given index
! There are 8 elements at the Xpoint, so we need 8 nodes
! Xup =1 for lower Xpoint and =2 for upper Xpoint
! Here is a picture of an Xpoint:
!
!
!                               vertical
!                                 line
!                            (towards axis)
!                                  \           separatrix
!                        *         \         * 
!                         * NODE-3 \ NODE-2 *
!                          *       \       *
!                           *      \      *
!                            *     \     *
!                             *    \    *
!                     NODE-4   *   \   *   NODE-1
!                               *  \  *
!                                * \ *
!                                 *\*
!                  ----------------\-------------horizontal line (or near horizontal)
!                                 *\*
!                                * \ *
!                               *  \  *
!                     NODE-5   *   \   *   NODE-8
!                             *    \    *
!                            *     \     *
!                           *      \      *
!                          *       \       *
!                         * NODE-6 \ NODE-7 *
!                        *         \         *
!                                  \          separatrix
!                                            
! Note that opposite nodes are exactly the same, only with opposite sign
! Hence, only nodes 1->4 will be added to the list 
!------------------------------------------------------------------------------------------

use tr_module 
use data_structure
use grid_xpoint_data
use mod_interp

implicit none

! --- Routine parameters
type (type_node_list)       , intent(inout) :: node_list
type (type_node_list)       , intent(inout) :: newnode_list
type (type_element_list)    , intent(inout) :: element_list
type (type_strategic_points), intent(in)    :: stpts
type (type_new_points)      , intent(in)    :: nwpts
integer,                      intent(in)    :: Xup, i_elm_xpoint(2)
real*8,                       intent(in)    :: R_axis, Z_axis, R_xpoint(2), Z_xpoint(2), s_xpoint(2), t_xpoint(2)

! --- local variables
integer             :: i, l, m, index
real*8              :: R_cub1d(4), Z_cub1d(4), dR_dt, dZ_dt, RZ_jac, PSI_R, PSI_Z, tmp1, tmp2
real*8              :: EJAC, RX, RY, SX, SY, CRR, CZZ, CRZ, alpha1, alpha2, alpha_max, alpha_min
real*8              :: RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss
real*8              :: ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss
real*8              :: PSg1,dPSg1_dr,dPSg1_ds,dPSg1_drs,dPSg1_drr,dPSg1_dss
real*8              :: angle
real*8,external     :: root

index = newnode_list%n_nodes  
do i=1,4

  index = index + 1
  
  newnode_list%node(index)%x(1,4,:) = 0.d0
  newnode_list%node(index)%boundary = 0
  
  if (Xup .eq. 1) angle = stpts%angle_LowerRight
  if (Xup .eq. 2) angle = stpts%angle_UpperRight
  
  call interp(node_list,element_list,i_elm_xpoint(Xup),1,1,s_xpoint(Xup),t_xpoint(Xup),&
                                     PSg1,dPSg1_dr,dPSg1_ds,dPSg1_drs,dPSg1_drr,dPSg1_dss)
  call interp_RZ(node_list,element_list,i_elm_xpoint(Xup),s_xpoint(Xup),t_xpoint(Xup), &
                                        RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss, &
                                        ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)
  
  EJAC = dRRg1_dr*dZZg1_ds - dRRg1_ds * dZZg1_dr
  RY = - dRRg1_ds / EJAC
  RX =   dZZg1_ds / EJAC
  SY =   dRRg1_dr / EJAC
  SX = - dZZg1_dr / EJAC

  CRR = (dPSg1_drr*RX*RX + 2.*dPSg1_drs*RX*SX + dPSg1_dss*SX*SX)/2.
  CZZ = (dPSg1_drr*RY*RY + 2.*dPSg1_drs*RY*SY + dPSg1_dss*SY*SY)/2.
  CRZ = (dPSg1_drr*RX*RY + dPSg1_drs*(RX*SY + RY*SX) + dPSg1_dss*SX*SY)/2.

  alpha1 = root(CRR,CRZ,CZZ,CRZ*CRZ-4.d0*CRR*CZZ,+1.d0)
  alpha2 = root(CRR,CRZ,CZZ,CRZ*CRZ-4.d0*CRR*CZZ,-1.d0)

  alpha_max = max(alpha1,alpha2)
  alpha_min = min(alpha1,alpha2)

  newnode_list%node(index)%x(1,1,:) = (/ R_xpoint(Xup), Z_xpoint(Xup) /)
  
  if (i .eq. 1) then
    newnode_list%node(index)%x(1,2,:) = (/ cos(angle),sin(angle) /)
    newnode_list%node(index)%x(1,3,:) = (/ alpha_max,1.d0 /) / sqrt(alpha_max**2 + 1.d0)
  endif
  if (i .eq. 2) then
    newnode_list%node(index)%x(1,2,:) = (/ R_xpoint(Xup)-R_axis,Z_xpoint(Xup)-Z_axis/)/sqrt((R_xpoint(Xup)-R_axis)**2+(Z_xpoint(Xup)-Z_axis)**2)
    newnode_list%node(index)%x(1,3,:) = (/ alpha_max,1.d0 /) / sqrt(alpha_max**2 + 1.d0)
  endif
  if (i .eq. 3) then
    newnode_list%node(index)%x(1,2,:) = (/ R_xpoint(Xup)-R_axis,Z_xpoint(Xup)-Z_axis/)/sqrt((R_xpoint(Xup)-R_axis)**2+(Z_xpoint(Xup)-Z_axis)**2)
    newnode_list%node(index)%x(1,3,:) = (/ alpha_min,1.d0 /) / sqrt(alpha_min**2 + 1.d0)
  endif
  if (i .eq. 4) then
    newnode_list%node(index)%x(1,2,:) = (/ cos(angle),sin(angle) /)
    newnode_list%node(index)%x(1,3,:) = (/ alpha_min,1.d0 /) / sqrt(alpha_min**2 + 1.d0)
  endif

enddo
newnode_list%n_nodes = index


return
end subroutine create_x_node
