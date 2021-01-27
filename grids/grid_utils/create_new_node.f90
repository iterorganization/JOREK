subroutine create_new_node(node_list, element_list, newnode_list, index, i2, j2, nwpts)
!------------------------------------------------------------------------------------------
! subroutine defines a new node with the given index, the node is the
! intersection between a flux surface and a polar line, at indices i2 and j2 respectively
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
type (type_new_points)      , intent(in)    :: nwpts
integer,                      intent(in)    :: i2, j2, index

! --- local variables
integer             :: m
real*8              :: R_cub1d(4), Z_cub1d(4), dR_dt, dZ_dt, RZ_jac, PSI_R, PSI_Z, tmp1, tmp2
real*8              :: RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss
real*8              :: ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss
real*8              :: PSg1,dPSg1_dr,dPSg1_ds,dPSg1_drs,dPSg1_drr,dPSg1_dss


  m = nwpts%k_cross(i2,j2)
  
  R_cub1d = (/ nwpts%R_polar(m,1,j2), 3.d0/2.d0 *(nwpts%R_polar(m,2,j2)-nwpts%R_polar(m,1,j2)), &
               nwpts%R_polar(m,4,j2), 3.d0/2.d0 *(nwpts%R_polar(m,4,j2)-nwpts%R_polar(m,3,j2))  /)
  Z_cub1d = (/ nwpts%Z_polar(m,1,j2), 3.d0/2.d0 *(nwpts%Z_polar(m,2,j2)-nwpts%Z_polar(m,1,j2)), &
               nwpts%Z_polar(m,4,j2), 3.d0/2.d0 *(nwpts%Z_polar(m,4,j2)-nwpts%Z_polar(m,3,j2)) /)
  
  call CUB1D(R_cub1d(1), R_cub1d(2), R_cub1d(3), R_cub1d(4),nwpts%t_tht(i2,j2),tmp1, dR_dt)
  call CUB1D(Z_cub1d(1), Z_cub1d(2), Z_cub1d(3), Z_cub1d(4),nwpts%t_tht(i2,j2),tmp2, dZ_dt)

  call interp_RZ(node_list,element_list,nwpts%ielm_flux(i2,j2),nwpts%s_flux(i2,j2),nwpts%t_flux(i2,j2), &
                 RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss, &
                 ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)

  call interp(node_list,element_list,nwpts%ielm_flux(i2,j2),1,1,nwpts%s_flux(i2,j2),nwpts%t_flux(i2,j2),&
                 PSg1,dPSg1_dr,dPSg1_ds,dPSg1_drs,dPSg1_drr,dPSg1_dss)

  RZ_jac  = DRRg1_dr * dZZg1_ds - dRRg1_ds * dZZg1_dr
  
  PSI_R = (   dPSg1_dr * dZZg1_ds - dPSg1_ds * dZZg1_dr ) / RZ_jac
  PSI_Z = ( - dPSg1_dr * dRRg1_ds + dPSg1_ds * dRRg1_dr ) / RZ_jac
  
  newnode_list%node(index)%x(1,1,:) = (/ nwpts%RR_new(i2,j2), nwpts%ZZ_new(i2,j2) /)
  newnode_list%node(index)%x(1,2,:) = (/ dR_dt, dZ_dt /)   / sqrt(dR_dt**2 + dZ_dt**2)
  newnode_list%node(index)%x(1,3,:) = (/ -PSI_Z, +PSI_R /) / sqrt(PSI_R**2 + PSI_Z**2)
  newnode_list%node(index)%x(1,4,:) = 0.d0
  newnode_list%node(index)%boundary = 0

return
end subroutine create_new_node






subroutine create_new_node_polar(newnode_list, index, n_seg, n_nodes, i1, i2, R_polar1, Z_polar1, R_polar2, Z_polar2, R_point, Z_point)
!------------------------------------------------------------------------------------------
! subroutine defines a new node with the given index, the node is the
! intersection between two polar lines, at indices i2 and j2 respectively
!------------------------------------------------------------------------------------------

use tr_module 
use data_structure
use grid_xpoint_data

implicit none

! --- Routine parameters
type (type_node_list)       , intent(inout) :: newnode_list
integer,                      intent(in)    :: index, n_seg, n_nodes, i1, i2
real*8,                       intent(in)    :: R_polar1(n_nodes-1,4,n_seg), Z_polar1(n_nodes-1,4,n_seg)
real*8,                       intent(in)    :: R_polar2(n_seg-1,4,n_nodes), Z_polar2(n_seg-1,4,n_nodes)
real*8,                       intent(in)    :: R_point, Z_point

! --- local variables
integer :: j1, j2
real*8  :: R_cub1d1(4), Z_cub1d1(4)
real*8  :: R_cub1d2(4), Z_cub1d2(4)
real*8  :: ss1, ss2
real*8  :: dR_dt1, dZ_dt1
real*8  :: dR_dt2, dZ_dt2
real*8  :: tmp1, tmp2


  if (i1 .eq. 0) then
    j1 = 1
    ss1 = -1.d0
  else
    j1 = i1
    ss1 = 1.d0
  endif
  if (i2 .eq. 0) then
    j2 = 1
    ss2 = -1.d0
  else
    j2 = i2
    ss2 = 1.d0
  endif
  
  call from_polar_to_cubic(R_polar1(j2,1:4,j1+1),R_cub1d1)
  call from_polar_to_cubic(Z_polar1(j2,1:4,j1+1),Z_cub1d1)
  call from_polar_to_cubic(R_polar2(j1,1:4,j2+1),R_cub1d2)
  call from_polar_to_cubic(Z_polar2(j1,1:4,j2+1),Z_cub1d2)
  
  call CUB1D(R_cub1d1(1), R_cub1d1(2), R_cub1d1(3), R_cub1d1(4),ss1,tmp1, dR_dt1)
  call CUB1D(Z_cub1d1(1), Z_cub1d1(2), Z_cub1d1(3), Z_cub1d1(4),ss1,tmp2, dZ_dt1)
  call CUB1D(R_cub1d2(1), R_cub1d2(2), R_cub1d2(3), R_cub1d2(4),ss2,tmp1, dR_dt2)
  call CUB1D(Z_cub1d2(1), Z_cub1d2(2), Z_cub1d2(3), Z_cub1d2(4),ss2,tmp2, dZ_dt2)

  newnode_list%node(index)%x(1,1,:) = (/ R_point, Z_point /)
  newnode_list%node(index)%x(1,2,:) = (/ dR_dt2, dZ_dt2 /)   / sqrt(dR_dt2**2 + dZ_dt2**2)
  newnode_list%node(index)%x(1,3,:) = (/ dR_dt1, dZ_dt1 /)   / sqrt(dR_dt1**2 + dZ_dt1**2)
  newnode_list%node(index)%x(1,4,:) = 0.d0
  newnode_list%node(index)%boundary = 0

return
end subroutine create_new_node_polar
