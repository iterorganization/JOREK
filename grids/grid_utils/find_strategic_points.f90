subroutine find_strategic_points(node_list, element_list, flux_list, xcase, force_horizontal_Xline, n_grids, stpts)
!----------------------------------------------------------------------------------------
! subroutine finds all the strategic points on the legs (Leg corners, strike points etc.)
!----------------------------------------------------------------------------------------

use constants
use tr_module 
use data_structure
use grid_xpoint_data
use phys_module, only: tokamak_device
use mod_interp, only: interp_RZ
use equil_info

implicit none

! --- Routine parameters
type (type_surface_list),     intent(inout) :: flux_list
type (type_node_list),        intent(inout) :: node_list
type (type_element_list),     intent(inout) :: element_list
type (type_strategic_points), intent(inout) :: stpts
integer,                      intent(in)    :: n_grids(12) 
integer,                      intent(in)    :: xcase
logical,                      intent(in)    :: force_horizontal_Xline


! --- local variables
integer  :: i, k, l, i_elm, i_surf, i_find, i_elm_find(8) , i_max, ifound, ier
integer  :: n_flux,   n_open,   n_outer,   n_inner,   n_private,   n_up_priv  
real*8   :: RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss
real*8   :: ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss
real*8   :: rr1, ss1, s_find(8), t_find(8)
real*8   :: tht_x, tht_x1, tht_x2
real*8   :: angle_LowerCorner 
real*8   :: angle_UpperCorner
real*8   :: R_wall_max

n_flux    = n_grids(1)
n_open    = n_grids(3); n_outer   = n_grids(4); n_inner = n_grids(5)
n_private = n_grids(6); n_up_priv = n_grids(7)

write(*,*) '****************************************'
write(*,*) '* X-point grid : Find strategic points *'
write(*,*) '****************************************'



!--------------------------------------------------------------------------------------------
!-------------------------------- Define all strategic points -------------------------------
!--------------------------------------------------------------------------------------------
if(xcase .eq. UPPER_XPOINT) then
  stpts%RLeftCorn_LowerInnerLeg  = 0.d0;    stpts%ZLeftCorn_LowerInnerLeg  = 0.d0
  stpts%RRightCorn_LowerInnerLeg = 0.d0;    stpts%ZRightCorn_LowerInnerLeg = 0.d0
  stpts%RLeftCorn_LowerOuterLeg  = 0.d0;    stpts%ZLeftCorn_LowerOuterLeg  = 0.d0
  stpts%RRightCorn_LowerOuterLeg = 0.d0;    stpts%ZRightCorn_LowerOuterLeg = 0.d0
  stpts%RStrike_LowerInnerLeg    = 0.d0;    stpts%ZStrike_LowerInnerLeg    = 0.d0
  stpts%RStrike_LowerOuterLeg    = 0.d0;    stpts%ZStrike_LowerOuterLeg    = 0.d0
else
  stpts%RLeftCorn_LowerInnerLeg  = 999.d0;  stpts%ZLeftCorn_LowerInnerLeg  = 1.d10
  stpts%RRightCorn_LowerInnerLeg = 1.d10;   stpts%ZRightCorn_LowerInnerLeg = 999.d0
  stpts%RLeftCorn_LowerOuterLeg  = -1.d10;  stpts%ZLeftCorn_LowerOuterLeg  = 999.d0
  stpts%RRightCorn_LowerOuterLeg = -1.d10;  stpts%ZRightCorn_LowerOuterLeg = 1.d10
  stpts%RStrike_LowerInnerLeg    = 999.d0;  stpts%ZStrike_LowerInnerLeg    = 1.d10
  stpts%RStrike_LowerOuterLeg    = 999.d0;  stpts%ZStrike_LowerOuterLeg    = 1.d10
endif

if(xcase .eq. LOWER_XPOINT) then
  stpts%RLeftCorn_UpperInnerLeg  = 0.d0;    stpts%ZLeftCorn_UpperInnerLeg  = 0.d0
  stpts%RRightCorn_UpperInnerLeg = 0.d0;    stpts%ZRightCorn_UpperInnerLeg = 0.d0
  stpts%RLeftCorn_UpperOuterLeg  = 0.d0;    stpts%ZLeftCorn_UpperOuterLeg  = 0.d0
  stpts%RRightCorn_UpperOuterLeg = 0.d0;    stpts%ZRightCorn_UpperOuterLeg = 0.d0
  stpts%RStrike_UpperInnerLeg    = 0.d0;    stpts%ZStrike_UpperInnerLeg    = 0.d0
  stpts%RStrike_UpperOuterLeg    = 0.d0;    stpts%ZStrike_UpperOuterLeg    = 0.d0
else
  stpts%RLeftCorn_UpperInnerLeg  = 999.d0;  stpts%ZLeftCorn_UpperInnerLeg  = -1.d10
  stpts%RRightCorn_UpperInnerLeg = 1.d10;   stpts%ZRightCorn_UpperInnerLeg = 999.d0
  stpts%RLeftCorn_UpperOuterLeg  = -1.d10;  stpts%ZLeftCorn_UpperOuterLeg  = 999.d0
  stpts%RRightCorn_UpperOuterLeg = -1.d10;  stpts%ZRightCorn_UpperOuterLeg = -1.d10
  stpts%RStrike_UpperInnerLeg    = 999.d0;  stpts%ZStrike_UpperInnerLeg    = -1.d10
  stpts%RStrike_UpperOuterLeg    = 999.d0;  stpts%ZStrike_UpperOuterLeg    = -1.d10
endif

if(xcase .ne. DOUBLE_NULL) then
  stpts%RSecondStrike_InnerLeg   = 0.d0;    stpts%ZSecondStrike_InnerLeg   = 0.d0       
  stpts%RSecondStrike_OuterLeg   = 0.d0;    stpts%ZSecondStrike_OuterLeg   = 0.d0       
else
  if ( ES%active_xpoint .eq. LOWER_XPOINT ) then
    stpts%RSecondStrike_InnerLeg = 999.d0;  stpts%ZSecondStrike_InnerLeg   = 1.d10   
    stpts%RSecondStrike_OuterLeg = 999.d0;  stpts%ZSecondStrike_OuterLeg   = 1.d10   
  else
    stpts%RSecondStrike_InnerLeg = 999.d0;  stpts%ZSecondStrike_InnerLeg   = -1.d10
    stpts%RSecondStrike_OuterLeg = 999.d0;  stpts%ZSecondStrike_OuterLeg   = -1.d10  
  endif
  stpts%ZLimit_UpperMastWall =  1.5d0 
  stpts%ZLimit_LowerMastWall = -1.5d0
endif

stpts%RMiddle_LowerPrivate       = 999.d0;  stpts%ZMiddle_LowerPrivate     = 1.d10
stpts%RMiddle_UpperPrivate       = 999.d0;  stpts%ZMiddle_UpperPrivate     = -1.d1



!--------------------------------------------------------------------------------------------
!------------------- Now find all the points for standard equilibria ------------------------
!--------------------------------------------------------------------------------------------

if (xcase .ne. DOUBLE_NULL) then
  ! ---------------------------------- The last open flux surface (SOL boundary)
  i_surf = n_flux+n_open 
  do k=1,flux_list%flux_surfaces(i_surf)%n_pieces    
    do l=1,3,2
      
      rr1   = flux_list%flux_surfaces(i_surf)%s(l,k)
      ss1   = flux_list%flux_surfaces(i_surf)%t(l,k)
      i_elm = flux_list%flux_surfaces(i_surf)%elm(k)

      call interp_RZ(node_list,element_list,i_elm,rr1,ss1,RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss, &
                                                          ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)

      if (xcase .eq. LOWER_XPOINT) then
        if ((ZZg1 .lt. stpts%ZLeftCorn_LowerInnerLeg) .and. (RRg1 .lt. ES%R_xpoint(1))) then
          stpts%RLeftCorn_LowerInnerLeg = RRg1
          stpts%ZLeftCorn_LowerInnerLeg = ZZg1
        endif
        if ((ZZg1 .lt. stpts%ZRightCorn_LowerOuterLeg) .and. (RRg1 .gt. ES%R_xpoint(1))) then
          stpts%RRightCorn_LowerOuterLeg = RRg1
          stpts%ZRightCorn_LowerOuterLeg = ZZg1
        endif
      else
        if ((ZZg1 .gt. stpts%ZLeftCorn_UpperInnerLeg) .and. (RRg1 .lt. ES%R_xpoint(2))) then
          stpts%RLeftCorn_UpperInnerLeg = RRg1
          stpts%ZLeftCorn_UpperInnerLeg = ZZg1
        endif
        if ((ZZg1 .gt. stpts%ZRightCorn_UpperOuterLeg) .and. (RRg1 .gt. ES%R_xpoint(2))) then
          stpts%RRightCorn_UpperOuterLeg = RRg1
          stpts%ZRightCorn_UpperOuterLeg = ZZg1
        endif
      endif
      
    enddo
  enddo

else ! xcase == DOUBLE_NULL
  ! ---------------------------------- The last open flux surface (SOL boundary) on outer board (LFS)
  i_surf = n_flux+n_open+n_outer  
  do k=1,flux_list%flux_surfaces(i_surf)%n_pieces    
    do l=1,3,2
      
      rr1   = flux_list%flux_surfaces(i_surf)%s(l,k)
      ss1   = flux_list%flux_surfaces(i_surf)%t(l,k)
      i_elm = flux_list%flux_surfaces(i_surf)%elm(k)

      call interp_RZ(node_list,element_list,i_elm,rr1,ss1,RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss, &
                                                          ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)

      if(tokamak_device(1:4) .ne. 'MAST') then
        if ((ZZg1 .lt. stpts%ZRightCorn_LowerOuterLeg) .and. (RRg1 .gt. ES%R_xpoint(1))) then
          stpts%RRightCorn_LowerOuterLeg = RRg1
          stpts%ZRightCorn_LowerOuterLeg = ZZg1
        endif
        if ((ZZg1 .gt. stpts%ZRightCorn_UpperOuterLeg) .and. (RRg1 .gt. ES%R_xpoint(2))) then
          stpts%RRightCorn_UpperOuterLeg = RRg1
          stpts%ZRightCorn_UpperOuterLeg = ZZg1
        endif
      else
        if ((ZZg1 .lt. ES%Z_xpoint(1)) .and. (RRg1 .gt. ES%R_xpoint(1)) .and. (RRg1 .gt. stpts%RRightCorn_LowerOuterLeg)) then
          stpts%RRightCorn_LowerOuterLeg = RRg1
          stpts%ZRightCorn_LowerOuterLeg = ZZg1
        endif
        if ((ZZg1 .gt. ES%Z_xpoint(2)) .and. (RRg1 .gt. ES%R_xpoint(2)) .and. (RRg1 .gt. stpts%RRightCorn_UpperOuterLeg)) then
          stpts%RRightCorn_UpperOuterLeg = RRg1
          stpts%ZRightCorn_UpperOuterLeg = ZZg1
        endif
      endif
      
    enddo
  enddo

  ! ---------------------------------- The last open flux surface (SOL boundary) on inner board (HFS)
  i_surf = n_flux+n_open+n_outer+n_inner 
  do k=1,flux_list%flux_surfaces(i_surf)%n_pieces    
    do l=1,3,2
      
      rr1   = flux_list%flux_surfaces(i_surf)%s(l,k)
      ss1   = flux_list%flux_surfaces(i_surf)%t(l,k)
      i_elm = flux_list%flux_surfaces(i_surf)%elm(k)

      call interp_RZ(node_list,element_list,i_elm,rr1,ss1,RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss, &
                                                          ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)

      if ((ZZg1 .lt. stpts%ZLeftCorn_LowerInnerLeg) .and. (RRg1 .lt. ES%R_xpoint(1))) then
        stpts%RLeftCorn_LowerInnerLeg = RRg1
        stpts%ZLeftCorn_LowerInnerLeg = ZZg1
      endif
      if ((ZZg1 .gt. stpts%ZLeftCorn_UpperInnerLeg) .and. (RRg1 .lt. ES%R_xpoint(2))) then
        stpts%RLeftCorn_UpperInnerLeg = RRg1
        stpts%ZLeftCorn_UpperInnerLeg = ZZg1
      endif
      
    enddo
  enddo
endif

! ---------------------------------- The last open flux surface (Private boundary) under lower X-point 
if (xcase .ne. UPPER_XPOINT) then
  i_surf = n_flux+n_open+n_outer+n_inner+n_private  
  do k=1,flux_list%flux_surfaces(i_surf)%n_pieces
    do l=1,3,2
    
      rr1   = flux_list%flux_surfaces(i_surf)%s(l,k)
      ss1   = flux_list%flux_surfaces(i_surf)%t(l,k)
      i_elm = flux_list%flux_surfaces(i_surf)%elm(k)

      call interp_RZ(node_list,element_list,i_elm,rr1,ss1,RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss, &
                                                          ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)

      if ((RRg1 .lt. stpts%RRightCorn_LowerInnerLeg) .and. (ZZg1 .lt. ES%Z_xpoint(1))) then
        stpts%RRightCorn_LowerInnerLeg = RRg1
        stpts%ZRightCorn_LowerInnerLeg = ZZg1
      endif
      if ((RRg1 .gt. stpts%RLeftCorn_LowerOuterLeg) .and. (ZZg1 .lt. ES%Z_xpoint(1))) then
        stpts%RLeftCorn_LowerOuterLeg = RRg1
        stpts%ZLeftCorn_LowerOuterLeg = ZZg1
      endif
    
    enddo
  enddo
endif

! ---------------------------------- The last open flux surface (Private boundary) above upper X-point 
if (xcase .ne. LOWER_XPOINT) then
  i_surf = n_flux+n_open+n_outer+n_inner+n_private+n_up_priv 
  do k=1,flux_list%flux_surfaces(i_surf)%n_pieces
    do l=1,3,2
    
      rr1   = flux_list%flux_surfaces(i_surf)%s(l,k)
      ss1   = flux_list%flux_surfaces(i_surf)%t(l,k)
      i_elm = flux_list%flux_surfaces(i_surf)%elm(k)

      call interp_RZ(node_list,element_list,i_elm,rr1,ss1,RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss, &
                                                          ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)

      if ((RRg1 .lt. stpts%RRightCorn_UpperInnerLeg) .and. (ZZg1 .gt. ES%Z_xpoint(2))) then
        stpts%RRightCorn_UpperInnerLeg = RRg1
        stpts%ZRightCorn_UpperInnerLeg = ZZg1
      endif
      if ((RRg1 .gt. stpts%RLeftCorn_UpperOuterLeg) .and. (ZZg1 .gt. ES%Z_xpoint(2))) then
        stpts%RLeftCorn_UpperOuterLeg = RRg1
        stpts%ZLeftCorn_UpperOuterLeg = ZZg1
      endif
    
    enddo
  enddo
endif  

! ---------------------------------- Find line from axis to lower X-point and get intersection ZMiddle_LowerPrivate with last private surface
if (xcase .ne. UPPER_XPOINT) then
  tht_x = atan2(ES%Z_xpoint(1)-ES%Z_axis,ES%R_xpoint(1)-ES%R_axis)
  if (tht_x .lt. 0.d0) tht_x = tht_x + 2.d0 * PI
  i_surf = n_flux+n_open+n_outer+n_inner+n_private
  call find_theta_surface(node_list,element_list,flux_list,i_surf,tht_x,ES%R_axis,ES%Z_axis,i_elm_find,s_find,t_find,i_find)

  do i=1,i_find
    call interp_RZ(node_list,element_list,i_elm_find(i),s_find(i),t_find(i),RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss, &
                                                                            ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)
    stpts%RMiddle_LowerPrivate = RRg1
    stpts%ZMiddle_LowerPrivate = ZZg1
    if (stpts%ZMiddle_LowerPrivate .le. ES%Z_xpoint(1)) exit
  enddo
endif

! ---------------------------------- Find line from axis to upper X-point and get intersection ZMiddle_UpperPrivate with last private surface
if (xcase .ne. LOWER_XPOINT) then
  tht_x = atan2(ES%Z_xpoint(2)-ES%Z_axis,ES%R_xpoint(2)-ES%R_axis)
  if (tht_x .lt. 0.d0) tht_x = tht_x + 2.d0 * PI
  i_surf = n_flux+n_open+n_outer+n_inner+n_private+n_up_priv
  call find_theta_surface(node_list,element_list,flux_list,i_surf,tht_x,ES%R_axis,ES%Z_axis,i_elm_find,s_find,t_find,i_find)

  do i=1,i_find
    call interp_RZ(node_list,element_list,i_elm_find(i),s_find(i),t_find(i),RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss, &
                                                                            ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)
    stpts%RMiddle_UpperPrivate = RRg1
    stpts%ZMiddle_UpperPrivate = ZZg1
    if (stpts%ZMiddle_UpperPrivate .ge. ES%Z_xpoint(2)) exit
  enddo
endif

! ---------------------------------- Find lower strike points
if (xcase .ne. UPPER_XPOINT) then
  if ( (xcase .eq. DOUBLE_NULL) .and. ( ES%active_xpoint .eq. UPPER_XPOINT ) ) then
    i_surf = n_flux + n_open
  else
    i_surf = n_flux
  endif
  do k=1,flux_list%flux_surfaces(i_surf)%n_pieces
    do l=1,3,2
    
      rr1   = flux_list%flux_surfaces(i_surf)%s(l,k)
      ss1   = flux_list%flux_surfaces(i_surf)%t(l,k)
      i_elm = flux_list%flux_surfaces(i_surf)%elm(k)

      call interp_RZ(node_list,element_list,i_elm,rr1,ss1,RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss, &
                                                          ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)

      if ((ZZg1 .lt. stpts%ZStrike_LowerInnerLeg) .and. (RRg1 .lt. ES%R_xpoint(1))) then
        stpts%RStrike_LowerInnerLeg = RRg1
        stpts%ZStrike_LowerInnerLeg = ZZg1
      endif
      if ((ZZg1 .lt. stpts%ZStrike_LowerOuterLeg) .and. (RRg1 .gt. ES%R_xpoint(1))) then
        stpts%RStrike_LowerOuterLeg = RRg1
        stpts%ZStrike_LowerOuterLeg = ZZg1
      endif

    enddo
  enddo
endif

! ---------------------------------- Find upper strike points
if (xcase .ne. LOWER_XPOINT) then
  if ( (xcase .eq. DOUBLE_NULL) .and. ( ES%active_xpoint .eq. LOWER_XPOINT ) ) then
    i_surf = n_flux + n_open
  else
    i_surf = n_flux
  endif
  do k=1,flux_list%flux_surfaces(i_surf)%n_pieces
    do l=1,3,2
    
      rr1   = flux_list%flux_surfaces(i_surf)%s(l,k)
      ss1   = flux_list%flux_surfaces(i_surf)%t(l,k)
      i_elm = flux_list%flux_surfaces(i_surf)%elm(k)

      call interp_RZ(node_list,element_list,i_elm,rr1,ss1,RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss, &
                                                          ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)

      if ((ZZg1 .gt. stpts%ZStrike_UpperInnerLeg) .and. (RRg1 .lt. ES%R_xpoint(2))) then
        stpts%RStrike_UpperInnerLeg = RRg1
        stpts%ZStrike_UpperInnerLeg = ZZg1
      endif
      if ((ZZg1 .gt. stpts%ZStrike_UpperOuterLeg) .and. (RRg1 .gt. ES%R_xpoint(2))) then
        stpts%RStrike_UpperOuterLeg = RRg1
        stpts%ZStrike_UpperOuterLeg = ZZg1
      endif

    enddo
  enddo
endif

! ---------------------------------- Find strike points of second separatrix
if (xcase .eq. DOUBLE_NULL) then
  i_surf = n_flux+n_open
  do k=1,flux_list%flux_surfaces(i_surf)%n_pieces
    do l=1,3,2
    
      rr1   = flux_list%flux_surfaces(i_surf)%s(l,k)
      ss1   = flux_list%flux_surfaces(i_surf)%t(l,k)
      i_elm = flux_list%flux_surfaces(i_surf)%elm(k)

      call interp_RZ(node_list,element_list,i_elm,rr1,ss1,RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss, &
                                                          ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)

      if ( ES%active_xpoint .eq. LOWER_XPOINT ) then
        if ((ZZg1 .lt. stpts%ZSecondStrike_InnerLeg) .and. (RRg1 .lt. ES%R_xpoint(1))) then
          stpts%RSecondStrike_InnerLeg = RRg1
          stpts%ZSecondStrike_InnerLeg = ZZg1
        endif
        if ((ZZg1 .lt. stpts%ZSecondStrike_OuterLeg) .and. (RRg1 .gt. ES%R_xpoint(1))) then
          stpts%RSecondStrike_OuterLeg = RRg1
          stpts%ZSecondStrike_OuterLeg = ZZg1
        endif
      else
        if ((ZZg1 .gt. stpts%ZSecondStrike_InnerLeg) .and. (RRg1 .lt. ES%R_xpoint(2))) then
          stpts%RSecondStrike_InnerLeg = RRg1
          stpts%ZSecondStrike_InnerLeg = ZZg1
        endif
        if ((ZZg1 .gt. stpts%ZSecondStrike_OuterLeg) .and. (RRg1 .gt. ES%R_xpoint(2))) then
          stpts%RSecondStrike_OuterLeg = RRg1
          stpts%ZSecondStrike_OuterLeg = ZZg1
        endif
      endif
      
    enddo
  enddo
endif

!----------------------------------- Define the lines separating the central and upper/lower parts of the grid
if (xcase .ne. UPPER_XPOINT) tht_x1 = atan2(ES%Z_xpoint(1)-ES%Z_axis,ES%R_xpoint(1)-ES%R_axis)
if (xcase .eq. UPPER_XPOINT) tht_x1 = atan2(ES%Z_xpoint(2)-ES%Z_axis,ES%R_xpoint(2)-ES%R_axis)
if (xcase .eq. DOUBLE_NULL ) tht_x2 = atan2(ES%Z_xpoint(2)-ES%Z_axis,ES%R_xpoint(2)-ES%R_axis)
if (tht_x1 .lt. 0.d0) tht_x1 = tht_x1 + 2.d0*PI
if (tht_x2 .lt. 0.d0) tht_x2 = tht_x2 + 2.d0*PI

! --- First define the lower angle
if(xcase .ne. UPPER_XPOINT) then
  
  angle_LowerCorner       = atan2(stpts%ZLeftCorn_LowerInnerLeg-ES%Z_xpoint(1),stpts%RLeftCorn_LowerInnerLeg-ES%R_xpoint(1))
  stpts%angle_LowerLeft   = tht_x1 + 1.5d0*PI
  stpts%angle_LowerRight  = tht_x1 + 0.5d0*PI
  ! --- Depending on the equilibrium, it may be better to have a horizontal line... (eg. near double-null at JET)
  if(force_horizontal_Xline) then
    stpts%angle_LowerLeft   = PI
    stpts%angle_LowerRight  = 0.d0
  endif
  if (angle_LowerCorner .lt. 0.d0)          angle_LowerCorner       = angle_LowerCorner       + 2.d0*PI
  if (stpts%angle_LowerLeft   .gt. 2.d0*PI) stpts%angle_LowerLeft   = stpts%angle_LowerLeft   - 2.d0*PI
  if (stpts%angle_LowerRight  .gt. 2.d0*PI) stpts%angle_LowerRight  = stpts%angle_LowerRight  - 2.d0*PI
  
  if (stpts%angle_LowerLeft .gt. angle_LowerCorner) then   ! check if Limit_LowerInnerLeg is above ZLeftCorn_LowerInnerLeg : if not adjust angle
    stpts%angle_LowerLeft  = angle_LowerCorner - 0.1
    stpts%angle_LowerRight = stpts%angle_LowerLeft - PI
    if (stpts%angle_LowerRight  .lt. 0.d0) stpts%angle_LowerRight  = stpts%angle_LowerRight + 2.d0*PI
  endif
  
endif

! --- Define the upper angle
if(xcase .ne. LOWER_XPOINT) then

  angle_UpperCorner       = atan2(stpts%ZLeftCorn_UpperInnerLeg-ES%Z_xpoint(2),stpts%RLeftCorn_UpperInnerLeg-ES%R_xpoint(2))
  stpts%angle_UpperLeft   = tht_x2 + 0.5d0*PI; if(xcase .eq. UPPER_XPOINT) stpts%angle_UpperLeft   = tht_x1 + 0.5d0*PI
  stpts%angle_UpperRight  = tht_x2 + 1.5d0*PI; if(xcase .eq. UPPER_XPOINT) stpts%angle_UpperRight  = tht_x1 + 1.5d0*PI
  ! --- Depending on the equilibrium, it may be better to have a horizontal line... (eg. near double-null at JET)
  if(force_horizontal_Xline) then
    stpts%angle_UpperLeft   = PI
    stpts%angle_UpperRight  = 0.d0
  endif
  if (angle_UpperCorner .lt. 0.d0)          angle_UpperCorner       = angle_UpperCorner + 2.d0*PI
  if (stpts%angle_UpperLeft   .gt. 2.d0*PI) stpts%angle_UpperLeft   = stpts%angle_UpperLeft   - 2.d0*PI
  if (stpts%angle_UpperRight  .gt. 2.d0*PI) stpts%angle_UpperRight  = stpts%angle_UpperRight  - 2.d0*PI
  
  if (stpts%angle_UpperLeft .lt. angle_UpperCorner) then   ! check if Limit_UpperInnerLeg is below ZLeftCorn_UpperInnerLeg : if not adjust angle
    stpts%angle_UpperLeft  = angle_UpperCorner + 0.1
    stpts%angle_UpperRight = stpts%angle_UpperLeft + PI
    if (stpts%angle_UpperRight  .gt. 2.d0*PI) stpts%angle_UpperRight  = stpts%angle_UpperRight  - 2.d0*PI
  endif
  
endif

! --- Then find the lower line
if(xcase .ne. UPPER_XPOINT) then
  
  i_max = n_flux + n_open + n_outer
  call find_theta_surface(node_list,element_list,flux_list,i_max,stpts%angle_LowerRight,ES%R_xpoint(1),ES%Z_xpoint(1),i_elm_find,s_find,t_find,i_find)
  if(i_find .eq. 0) return
  call interp_RZ(node_list,element_list,i_elm_find(1),s_find(1),t_find(1),&
                 RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss,    &
                 ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)

  stpts%RLimit_LowerOuterLeg = RRg1
  stpts%ZLimit_LowerOuterLeg = ZZg1

  i_max = n_flux + n_open + n_outer + n_inner
  call find_theta_surface(node_list,element_list,flux_list,i_max,stpts%angle_LowerLeft,ES%R_xpoint(1),ES%Z_xpoint(1),i_elm_find,s_find,t_find,i_find)
  if(i_find .eq. 0) return
  call interp_RZ(node_list,element_list,i_elm_find(1),s_find(1),t_find(1),&
                 RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss,    &
                 ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)

  stpts%RLimit_LowerInnerLeg = RRg1
  stpts%ZLimit_LowerInnerLeg = ZZg1

endif

! --- And the upper line
if(xcase .ne. LOWER_XPOINT) then
  
  i_max = n_flux + n_open + n_outer
  call find_theta_surface(node_list,element_list,flux_list,i_max,stpts%angle_UpperRight,ES%R_xpoint(2),ES%Z_xpoint(2),i_elm_find,s_find,t_find,i_find)
  if(i_find .eq. 0) return
  call interp_RZ(node_list,element_list,i_elm_find(1),s_find(1),t_find(1),&
                 RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss,    &
                 ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)

  stpts%RLimit_UpperOuterLeg = RRg1
  stpts%ZLimit_UpperOuterLeg = ZZg1

  i_max = n_flux + n_open + n_outer + n_inner
  call find_theta_surface(node_list,element_list,flux_list,i_max,stpts%angle_UpperLeft,ES%R_xpoint(2),ES%Z_xpoint(2),i_elm_find,s_find,t_find,i_find)
  if(i_find .eq. 0) return
  call interp_RZ(node_list,element_list,i_elm_find(1),s_find(1),t_find(1),&
                 RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss,    &
                 ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)

  stpts%RLimit_UpperInnerLeg = RRg1
  stpts%ZLimit_UpperInnerLeg = ZZg1

endif



!--------------------------------------------------------------------------------------------
!------------------- Now find points related to MAST equilibria -----------------------------
!--------------------------------------------------------------------------------------------

if( (tokamak_device(1:4) .eq. 'MAST') .and. (xcase .eq. DOUBLE_NULL) ) then

  ! ---------------------------------- The last open flux surface on outer board that intersects MAST lower Wall
  ifound = 0
  R_wall_max = 1.73d0
  do i=n_flux+n_open, n_flux+n_open+n_outer 
    do k=1,flux_list%flux_surfaces(i)%n_pieces    
      do l=1,3,2
        
        rr1   = flux_list%flux_surfaces(i)%s(l,k)
        ss1   = flux_list%flux_surfaces(i)%t(l,k)
        i_elm = flux_list%flux_surfaces(i)%elm(k)

        call interp_RZ(node_list,element_list,i_elm,rr1,ss1,RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss, &
                                                            ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)

        if ((ZZg1 .lt. stpts%ZLimit_LowerMastWall) .and. (RRg1 .gt. R_wall_max)) then
          stpts%RLimit_LowerMastWall = RRg1
          stpts%ZLimit_LowerMastWall = ZZg1
          ifound = 1
        endif
        
      enddo
    enddo
    stpts%i_surf_wall_low = i 
    if (ifound .eq. 1) exit
  enddo

  ! ---------------------------------- The last open flux surface on outer board that intersects MAST upper Wall
  ifound = 0
  R_wall_max = 1.723d0
  do i=n_flux+n_open, n_flux+n_open+n_outer 
    do k=1,flux_list%flux_surfaces(i)%n_pieces    
      do l=1,3,2
        
        rr1   = flux_list%flux_surfaces(i)%s(l,k)
        ss1   = flux_list%flux_surfaces(i)%t(l,k)
        i_elm = flux_list%flux_surfaces(i)%elm(k)

        call interp_RZ(node_list,element_list,i_elm,rr1,ss1,RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss, &
                                                            ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)

        if ((ZZg1 .gt. stpts%ZLimit_UpperMastWall) .and. (RRg1 .gt. R_wall_max)) then
          stpts%RLimit_UpperMastWall = RRg1
          stpts%ZLimit_UpperMastWall = ZZg1
          ifound = 1
        endif
        
      enddo
    enddo
    stpts%i_surf_wall_up = i 
    if (ifound .eq. 1) exit
  enddo

  ! ---------------------------------- The last lower wall surface of MAST
  call find_theta_surface(node_list,element_list,flux_list,stpts%i_surf_wall_low,stpts%angle_LowerRight,ES%R_xpoint(1),ES%Z_xpoint(1),i_elm_find,s_find,t_find,i_find)
  if(i_find .eq. 0) return
  call interp_RZ(node_list,element_list,i_elm_find(1),s_find(1),t_find(1),&
                 RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss,    &
                 ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)

  stpts%RLimit_LowerMastWallBox = RRg1
  stpts%ZLimit_LowerMastWallBox = ZZg1

  ! ---------------------------------- The last upper wall surface of MAST
  call find_theta_surface(node_list,element_list,flux_list,stpts%i_surf_wall_up,stpts%angle_UpperRight,ES%R_xpoint(2),ES%Z_xpoint(2),i_elm_find,s_find,t_find,i_find)
  if(i_find .eq. 0) return
  call interp_RZ(node_list,element_list,i_elm_find(1),s_find(1),t_find(1),&
                 RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss,    &
                 ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)

  stpts%RLimit_UpperMastWallBox = RRg1
  stpts%ZLimit_UpperMastWallBox = ZZg1

endif



!--------------------------------------------------------------------------------------------
!------------------- And print the output ---------------------------------------------------
!--------------------------------------------------------------------------------------------

write(*,'(A)')                  ' _________________________________________________________'
write(*,'(A)')                  '|                                                         |'
write(*,'(A)')                  '| LEG POINTS (R,Z)                                        |'
write(*,'(A)')                  '|_________________________________________________________|'
write(*,'(A)')                  '|                                                         |'

if (xcase .ne. UPPER_XPOINT) then
  write(*,'(A)')                '| Lower Legs : -------------------------------------------|'
  write(*,'(A,F5.2,A,F5.2,A)')  '|   Upper corner of the Lower Inner Leg  : (',stpts%RLimit_LowerInnerLeg,     ', ', stpts%ZLimit_LowerInnerLeg,     ') |'
  write(*,'(A,F5.2,A,F5.2,A)')  '|   Left corner  of the Lower Inner Leg  : (',stpts%RLeftCorn_LowerInnerLeg,  ', ', stpts%ZLeftCorn_LowerInnerLeg,  ') |'
  write(*,'(A,F5.2,A,F5.2,A)')  '|   Strike point of the Lower Inner Leg  : (',stpts%RStrike_LowerInnerLeg,    ', ', stpts%ZStrike_LowerInnerLeg,    ') |'
  write(*,'(A,F5.2,A,F5.2,A)')  '|   Right corner of the Lower Inner Leg  : (',stpts%RRightCorn_LowerInnerLeg, ', ', stpts%ZRightCorn_LowerInnerLeg, ') |'
  write(*,'(A,F5.2,A,F5.2,A)')  '|   Middle of the Lower Private surface  : (',stpts%RMiddle_LowerPrivate,     ', ', stpts%ZMiddle_LowerPrivate ,    ') |'
  write(*,'(A,F5.2,A,F5.2,A)')  '|   Left corner  of the Lower Outer Leg  : (',stpts%RLeftCorn_LowerOuterLeg,  ', ', stpts%ZLeftCorn_LowerOuterLeg,  ') |'
  write(*,'(A,F5.2,A,F5.2,A)')  '|   Strike point of the Lower Outer Leg  : (',stpts%RStrike_LowerOuterLeg,    ', ', stpts%ZStrike_LowerOuterLeg,    ') |'
  write(*,'(A,F5.2,A,F5.2,A)')  '|   Right corner of the Lower Outer Leg  : (',stpts%RRightCorn_LowerOuterLeg, ', ', stpts%ZRightCorn_LowerOuterLeg, ') |'
  write(*,'(A,F5.2,A,F5.2,A)')  '|   Upper corner of the Lower Outer Leg  : (',stpts%RLimit_LowerOuterLeg,     ', ', stpts%ZLimit_LowerOuterLeg,     ') |'
endif

if (xcase .ne. LOWER_XPOINT) then
  write(*,'(A)')                '| Upper Legs : -------------------------------------------|'
  write(*,'(A,F5.2,A,F5.2,A)')  '|   Lower corner of the Upper Inner Leg  : (',stpts%RLimit_UpperInnerLeg,     ', ', stpts%ZLimit_UpperInnerLeg,     ') |'
  write(*,'(A,F5.2,A,F5.2,A)')  '|   Left corner  of the Upper Inner Leg  : (',stpts%RLeftCorn_UpperInnerLeg,  ', ', stpts%ZLeftCorn_UpperInnerLeg,  ') |'
  write(*,'(A,F5.2,A,F5.2,A)')  '|   Strike point of the Upper Inner Leg  : (',stpts%RStrike_UpperInnerLeg,    ', ', stpts%ZStrike_UpperInnerLeg,    ') |'
  write(*,'(A,F5.2,A,F5.2,A)')  '|   Right corner of the Upper Inner Leg  : (',stpts%RRightCorn_UpperInnerLeg, ', ', stpts%ZRightCorn_UpperInnerLeg, ') |'
  write(*,'(A,F5.2,A,F5.2,A)')  '|   Middle of the Upper Private surface  : (',stpts%RMiddle_UpperPrivate,     ', ', stpts%ZMiddle_UpperPrivate ,    ') |'
  write(*,'(A,F5.2,A,F5.2,A)')  '|   Left corner  of the Upper Outer Leg  : (',stpts%RLeftCorn_UpperOuterLeg,  ', ', stpts%ZLeftCorn_UpperOuterLeg,  ') |'
  write(*,'(A,F5.2,A,F5.2,A)')  '|   Strike point of the Upper Outer Leg  : (',stpts%RStrike_UpperOuterLeg,    ', ', stpts%ZStrike_UpperOuterLeg,    ') |'
  write(*,'(A,F5.2,A,F5.2,A)')  '|   Right corner of the Upper Outer Leg  : (',stpts%RRightCorn_UpperOuterLeg, ', ', stpts%ZRightCorn_UpperOuterLeg, ') |'
  write(*,'(A,F5.2,A,F5.2,A)')  '|   Lower corner of the Upper Outer Leg  : (',stpts%RLimit_UpperOuterLeg,     ', ', stpts%ZLimit_UpperOuterLeg,     ') |'
endif

if (xcase .eq. DOUBLE_NULL) then
  write(*,'(A)')                '| Secondary Strike Points : ------------------------------|'
  write(*,'(A,F5.2,A,F5.2,A)')  '|   Left  Strike point of 2nd separatrix : (',stpts%RSecondStrike_InnerLeg,   ', ', stpts%ZSecondStrike_InnerLeg,   ') |'
  write(*,'(A,F5.2,A,F5.2,A)')  '|   Right Strike point of 2nd separatrix : (',stpts%RSecondStrike_OuterLeg,   ', ', stpts%ZSecondStrike_OuterLeg,   ') |'
endif

if(tokamak_device(1:4) .eq. 'MAST') then
  write(*,'(A)')                '| Last Outer Surface Intersecting MAST Wall : ------------|'
  write(*,'(A,F5.2,A,F5.2,A)')  '|   Intersection with lower MAST wall    : (',stpts%RLimit_LowerMastWall,     ', ', stpts%ZLimit_LowerMastWall,     ') |'
  write(*,'(A,F5.2,A,F5.2,A)')  '|   Intersection with upper MAST wall    : (',stpts%RLimit_UpperMastWall,     ', ', stpts%ZLimit_UpperMastWall,     ') |'
  write(*,'(A,F5.2,A,F5.2,A)')  '|   Intersection with lower Xpoint line  : (',stpts%RLimit_LowerMastWallBox,  ', ', stpts%ZLimit_LowerMastWallBox,  ') |'
  write(*,'(A,F5.2,A,F5.2,A)')  '|   Intersection with upper Xpoint line  : (',stpts%RLimit_UpperMastWallBox,  ', ', stpts%ZLimit_UpperMastWallBox,  ') |'
endif

write(*,'(A)')                  '|_________________________________________________________|'


return
end subroutine find_strategic_points

























subroutine find_strategic_points_advanced(node_list, element_list, flux_list, xcase, force_horizontal_Xline, n_grids, stpts)
!----------------------------------------------------------------------------------------
! subroutine finds all the strategic points on the legs (Leg corners, strike points etc.)
!----------------------------------------------------------------------------------------

use constants
use tr_module 
use data_structure
use grid_xpoint_data
use phys_module, only:   tokamak_device
use mod_interp, only: interp_RZ
use equil_info

implicit none

! --- Routine parameters
type (type_surface_list),     intent(inout) :: flux_list
type (type_node_list),        intent(inout) :: node_list
type (type_element_list),     intent(inout) :: element_list
type (type_strategic_points), intent(inout) :: stpts
integer,                      intent(in)    :: n_grids(12) 
integer,                      intent(in)    :: xcase
logical,                      intent(in)    :: force_horizontal_Xline


! --- local variables
integer  :: i, k, l, i_elm, i_surf, i_find, i_elm_find(8) , i_max, ifound, ier
integer  :: n_flux,   n_open,   n_outer,   n_inner,   n_private,   n_up_priv  
real*8   :: RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss
real*8   :: ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss
real*8   :: rr1, ss1, s_find(8), t_find(8)
real*8   :: tht_x, tht_x1, tht_x2
real*8   :: angle_LowerCorner 
real*8   :: angle_UpperCorner
real*8   :: R_wall_max

integer  :: count_lower, count_upper, count
integer  :: edge_piece(4), edge_side(4)
integer  :: i_part
real*8   :: R_tmp(4), Z_tmp(4)
real*8   :: distance, distance_min

logical, parameter :: debug = .false.

n_flux    = n_grids(1)
n_open    = n_grids(3); n_outer   = n_grids(4); n_inner = n_grids(5)
n_private = n_grids(6); n_up_priv = n_grids(7)
tht_x  = 0.d0
tht_x1 = 0.d0
tht_x2 = 0.d0

write(*,*) '****************************************'
write(*,*) '* X-point grid : Find strategic points *'
write(*,*) '****************************************'



!--------------------------------------------------------------------------------------------
!-------------------------------- Define all strategic points -------------------------------
!--------------------------------------------------------------------------------------------
if(xcase .eq. UPPER_XPOINT) then
  stpts%RLeftCorn_LowerInnerLeg  = 0.d0;    stpts%ZLeftCorn_LowerInnerLeg  = 0.d0
  stpts%RRightCorn_LowerInnerLeg = 0.d0;    stpts%ZRightCorn_LowerInnerLeg = 0.d0
  stpts%RLeftCorn_LowerOuterLeg  = 0.d0;    stpts%ZLeftCorn_LowerOuterLeg  = 0.d0
  stpts%RRightCorn_LowerOuterLeg = 0.d0;    stpts%ZRightCorn_LowerOuterLeg = 0.d0
  stpts%RStrike_LowerInnerLeg    = 0.d0;    stpts%ZStrike_LowerInnerLeg    = 0.d0
  stpts%RStrike_LowerOuterLeg    = 0.d0;    stpts%ZStrike_LowerOuterLeg    = 0.d0
else
  stpts%RLeftCorn_LowerInnerLeg  = 999.d0;  stpts%ZLeftCorn_LowerInnerLeg  = 1.d10
  stpts%RRightCorn_LowerInnerLeg = 1.d10;   stpts%ZRightCorn_LowerInnerLeg = 999.d0
  stpts%RLeftCorn_LowerOuterLeg  = -1.d10;  stpts%ZLeftCorn_LowerOuterLeg  = 999.d0
  stpts%RRightCorn_LowerOuterLeg = -1.d10;  stpts%ZRightCorn_LowerOuterLeg = 1.d10
  stpts%RStrike_LowerInnerLeg    = 999.d0;  stpts%ZStrike_LowerInnerLeg    = 1.d10
  stpts%RStrike_LowerOuterLeg    = 999.d0;  stpts%ZStrike_LowerOuterLeg    = 1.d10
endif

if(xcase .eq. LOWER_XPOINT) then
  stpts%RLeftCorn_UpperInnerLeg  = 0.d0;    stpts%ZLeftCorn_UpperInnerLeg  = 0.d0
  stpts%RRightCorn_UpperInnerLeg = 0.d0;    stpts%ZRightCorn_UpperInnerLeg = 0.d0
  stpts%RLeftCorn_UpperOuterLeg  = 0.d0;    stpts%ZLeftCorn_UpperOuterLeg  = 0.d0
  stpts%RRightCorn_UpperOuterLeg = 0.d0;    stpts%ZRightCorn_UpperOuterLeg = 0.d0
  stpts%RStrike_UpperInnerLeg    = 0.d0;    stpts%ZStrike_UpperInnerLeg    = 0.d0
  stpts%RStrike_UpperOuterLeg    = 0.d0;    stpts%ZStrike_UpperOuterLeg    = 0.d0
else
  stpts%RLeftCorn_UpperInnerLeg  = 999.d0;  stpts%ZLeftCorn_UpperInnerLeg  = -1.d10
  stpts%RRightCorn_UpperInnerLeg = 1.d10;   stpts%ZRightCorn_UpperInnerLeg = 999.d0
  stpts%RLeftCorn_UpperOuterLeg  = -1.d10;  stpts%ZLeftCorn_UpperOuterLeg  = 999.d0
  stpts%RRightCorn_UpperOuterLeg = -1.d10;  stpts%ZRightCorn_UpperOuterLeg = -1.d10
  stpts%RStrike_UpperInnerLeg    = 999.d0;  stpts%ZStrike_UpperInnerLeg    = -1.d10
  stpts%RStrike_UpperOuterLeg    = 999.d0;  stpts%ZStrike_UpperOuterLeg    = -1.d10
endif

if(xcase .ne. DOUBLE_NULL) then
  stpts%RSecondStrike_InnerLeg   = 0.d0;    stpts%ZSecondStrike_InnerLeg   = 0.d0
  stpts%RSecondStrike_OuterLeg   = 0.d0;    stpts%ZSecondStrike_OuterLeg   = 0.d0
else ! xcase == DOUBLE_NULL
  if ( ES%active_xpoint .eq. LOWER_XPOINT ) then
    stpts%RSecondStrike_InnerLeg = 999.d0;  stpts%ZSecondStrike_InnerLeg   = 1.d10
    stpts%RSecondStrike_OuterLeg = 999.d0;  stpts%ZSecondStrike_OuterLeg   = 1.d10
  else
    stpts%RSecondStrike_InnerLeg = 999.d0;  stpts%ZSecondStrike_InnerLeg   = -1.d10
    stpts%RSecondStrike_OuterLeg = 999.d0;  stpts%ZSecondStrike_OuterLeg   = -1.d10
  endif
  stpts%ZLimit_UpperMastWall =  1.5d0 
  stpts%ZLimit_LowerMastWall = -1.5d0
endif

stpts%RMiddle_LowerPrivate       = 999.d0;  stpts%ZMiddle_LowerPrivate     = 1.d10
stpts%RMiddle_UpperPrivate       = 999.d0;  stpts%ZMiddle_UpperPrivate     = -1.d1



!--------------------------------------------------------------------------------------------
!------------------- Now find all the points for standard equilibria ------------------------
!--------------------------------------------------------------------------------------------


! ---------------------------------- Find main strike points for symmetric case
if ( (xcase .eq. DOUBLE_NULL) .and. (ES%active_xpoint .eq. SYMMETRIC_XPOINT) ) then
  if (debug) write(*,*)'looking for strike points of up-down symmetric DN'
  count_lower = 0
  count_upper = 0
  i_surf = n_flux
  do i_part=1,flux_list%flux_surfaces(i_surf)%n_parts
    edge_piece(1) = flux_list%flux_surfaces(i_surf)%parts_index(i_part)
    edge_side(1)  = 1
    edge_piece(2) = flux_list%flux_surfaces(i_surf)%parts_index(i_part+1)-1
    edge_side(2)  = 3
    do l=1,2
      rr1   = flux_list%flux_surfaces(i_surf)%s(edge_side(l),edge_piece(l))
      ss1   = flux_list%flux_surfaces(i_surf)%t(edge_side(l),edge_piece(l))
      i_elm = flux_list%flux_surfaces(i_surf)%elm(edge_piece(l))
      call interp_RZ(node_list,element_list,i_elm,rr1,ss1,RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss, &
                                                          ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)
      if (ZZg1 .lt. ES%Z_axis) then
        count_lower = count_lower + 1
        R_tmp(count_lower) = RRg1
        Z_tmp(count_lower) = ZZg1
      else
        count_upper = count_upper + 1
        R_tmp(2+count_upper) = RRg1
        Z_tmp(2+count_upper) = ZZg1
      endif
    enddo
  enddo
  if (count_lower+count_upper .ne. 4) then
    write(*,*)'Something wrong, we are supposed to find 4 strikes in symmetric cases.',count_lower+count_upper
    write(*,*)'Aborting...'
    stop
  endif
  if (R_tmp(1) .lt. R_tmp(2)) then
    stpts%RStrike_LowerInnerLeg = R_tmp(1)
    stpts%ZStrike_LowerInnerLeg = Z_tmp(1)
    stpts%RStrike_LowerOuterLeg = R_tmp(2)
    stpts%ZStrike_LowerOuterLeg = Z_tmp(2)
  else
    stpts%RStrike_LowerInnerLeg = R_tmp(2)
    stpts%ZStrike_LowerInnerLeg = Z_tmp(2)
    stpts%RStrike_LowerOuterLeg = R_tmp(1)
    stpts%ZStrike_LowerOuterLeg = Z_tmp(1)
  endif
  if (R_tmp(3) .lt. R_tmp(4)) then
    stpts%RStrike_UpperInnerLeg = R_tmp(3)
    stpts%ZStrike_UpperInnerLeg = Z_tmp(3)
    stpts%RStrike_UpperOuterLeg = R_tmp(4)
    stpts%ZStrike_UpperOuterLeg = Z_tmp(4)
  else
    stpts%RStrike_UpperInnerLeg = R_tmp(4)
    stpts%ZStrike_UpperInnerLeg = Z_tmp(4)
    stpts%RStrike_UpperOuterLeg = R_tmp(3)
    stpts%ZStrike_UpperOuterLeg = Z_tmp(3)
  endif
  if (debug) write(*,*)'FOUND THE MAIN STRIKE POINTS:' 
  if (debug) write(*,*)'Lower Inner',stpts%RStrike_LowerInnerLeg,stpts%ZStrike_LowerInnerLeg
  if (debug) write(*,*)'Lower Outer',stpts%RStrike_LowerOuterLeg,stpts%ZStrike_LowerOuterLeg
  if (debug) write(*,*)'Upper Inner',stpts%RStrike_UpperInnerLeg,stpts%ZStrike_UpperInnerLeg
  if (debug) write(*,*)'Upper Outer',stpts%RStrike_UpperOuterLeg,stpts%ZStrike_UpperOuterLeg
else
! ---------------------------------- Find main strike points  
  if (debug) write(*,*)'looking for main strike points'
  count = 0
  i_surf = n_flux
  do i_part=1,flux_list%flux_surfaces(i_surf)%n_parts
    edge_piece(1) = flux_list%flux_surfaces(i_surf)%parts_index(i_part)
    edge_side(1)  = 1
    edge_piece(2) = flux_list%flux_surfaces(i_surf)%parts_index(i_part+1)-1
    edge_side(2)  = 3
    do l=1,2
      rr1   = flux_list%flux_surfaces(i_surf)%s(edge_side(l),edge_piece(l))
      ss1   = flux_list%flux_surfaces(i_surf)%t(edge_side(l),edge_piece(l))
      i_elm = flux_list%flux_surfaces(i_surf)%elm(edge_piece(l))
      call interp_RZ(node_list,element_list,i_elm,rr1,ss1,RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss, &
                                                          ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)
      if ( (xcase .eq. LOWER_XPOINT) .and. (ZZg1 .lt. ES%Z_axis)) then
        count = count + 1
        R_tmp(count) = RRg1
        Z_tmp(count) = ZZg1
      endif
      if ( (xcase .eq. UPPER_XPOINT) .and. (ZZg1 .gt. ES%Z_axis)) then
        count = count + 1
        R_tmp(count) = RRg1
        Z_tmp(count) = ZZg1
      endif
      if (xcase .eq. DOUBLE_NULL ) then
        if ( (ES%active_xpoint .eq. LOWER_XPOINT) .and. (ZZg1 .lt. ES%Z_axis) ) then
          count = count + 1
          R_tmp(count) = RRg1
          Z_tmp(count) = ZZg1
        endif
        if ( (ES%active_xpoint .eq. UPPER_XPOINT) .and. (ZZg1 .gt. ES%Z_axis) ) then
          count = count + 1
          R_tmp(count) = RRg1
          Z_tmp(count) = ZZg1
        endif
      endif
    enddo
    !if ( (xcase .eq. LOWER_XPOINT) .and. ( (count .eq. 2) .or. (count .eq. 4) ) ) then
    !  if (.not. (     ((R_tmp(count  ) .lt. ES%R_xpoint(1)) .and. (R_tmp(count-1) .gt. ES%R_xpoint(1))) &
    !            .or. ((R_tmp(count-1) .lt. ES%R_xpoint(1)) .and. (R_tmp(count  ) .gt. ES%R_xpoint(1))) ) ) then
    !   count = count - 2
    !  endif
    !endif
    !if ( (xcase .eq. UPPER_XPOINT) .and. ( (count .eq. 2) .or. (count .eq. 4) ) ) then
    !  if (.not. (     ((R_tmp(count  ) .lt. ES%R_xpoint(2)) .and. (R_tmp(count-1) .gt. ES%R_xpoint(2))) &
    !            .or. ((R_tmp(count-1) .lt. ES%R_xpoint(2)) .and. (R_tmp(count  ) .gt. ES%R_xpoint(2))) ) ) then
    !   count = count - 2
    !  endif
    !endif
  enddo
  if (count .ne. 2) then
    write(*,*)'Something wrong, we are supposed to find 2 strikes points.',count
    write(*,*)'Aborting...'
    stop
  endif
  if (Z_tmp(1) .lt. ES%Z_axis) then
    if (R_tmp(1) .lt. R_tmp(2)) then
      stpts%RStrike_LowerInnerLeg = R_tmp(1)
      stpts%ZStrike_LowerInnerLeg = Z_tmp(1)
      stpts%RStrike_LowerOuterLeg = R_tmp(2)
      stpts%ZStrike_LowerOuterLeg = Z_tmp(2)
    else
      stpts%RStrike_LowerInnerLeg = R_tmp(2)
      stpts%ZStrike_LowerInnerLeg = Z_tmp(2)
      stpts%RStrike_LowerOuterLeg = R_tmp(1)
      stpts%ZStrike_LowerOuterLeg = Z_tmp(1)
    endif
    if (debug) write(*,*)'FOUND THE MAIN STRIKE POINTS:' 
    if (debug) write(*,*)'Lower Inner',stpts%RStrike_LowerInnerLeg,stpts%ZStrike_LowerInnerLeg
    if (debug) write(*,*)'Lower Outer',stpts%RStrike_LowerOuterLeg,stpts%ZStrike_LowerOuterLeg
  else
    if (R_tmp(1) .lt. R_tmp(2)) then
      stpts%RStrike_UpperInnerLeg = R_tmp(1)
      stpts%ZStrike_UpperInnerLeg = Z_tmp(1)
      stpts%RStrike_UpperOuterLeg = R_tmp(2)
      stpts%ZStrike_UpperOuterLeg = Z_tmp(2)
    else
      stpts%RStrike_UpperInnerLeg = R_tmp(2)
      stpts%ZStrike_UpperInnerLeg = Z_tmp(2)
      stpts%RStrike_UpperOuterLeg = R_tmp(1)
      stpts%ZStrike_UpperOuterLeg = Z_tmp(1)
    endif
    if (debug) write(*,*)'FOUND THE MAIN STRIKE POINTS:' 
    if (debug) write(*,*)'Upper Inner',stpts%RStrike_UpperInnerLeg,stpts%ZStrike_UpperInnerLeg
    if (debug) write(*,*)'Upper Outer',stpts%RStrike_UpperOuterLeg,stpts%ZStrike_UpperOuterLeg
  endif
! ---------------------------------- Second main strike points  
  if (xcase .eq. DOUBLE_NULL) then
    if (debug) write(*,*)'looking for secondary strike points'
    count = 0
    i_surf = n_flux + n_open
    do i_part=1,flux_list%flux_surfaces(i_surf)%n_parts
      edge_piece(1) = flux_list%flux_surfaces(i_surf)%parts_index(i_part)
      edge_side(1)  = 1
      edge_piece(2) = flux_list%flux_surfaces(i_surf)%parts_index(i_part+1)-1
      edge_side(2)  = 3
      do l=1,2
        rr1   = flux_list%flux_surfaces(i_surf)%s(edge_side(l),edge_piece(l))
        ss1   = flux_list%flux_surfaces(i_surf)%t(edge_side(l),edge_piece(l))
        i_elm = flux_list%flux_surfaces(i_surf)%elm(edge_piece(l))
        call interp_RZ(node_list,element_list,i_elm,rr1,ss1,RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss, &
                                                            ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)
        if (xcase .eq. DOUBLE_NULL) then
          if ( (ES%active_xpoint .eq. LOWER_XPOINT) .and. (ZZg1 .gt. ES%Z_axis) ) then
            count = count + 1
            R_tmp(count) = RRg1
            Z_tmp(count) = ZZg1
          endif
          if ( (ES%active_xpoint .eq. UPPER_XPOINT) .and. (ZZg1 .lt. ES%Z_axis) ) then
            count = count + 1
            R_tmp(count) = RRg1
            Z_tmp(count) = ZZg1
          endif
        endif
      enddo
      !if ( (xcase .eq. LOWER_XPOINT) .and. ( (count .eq. 2) .or. (count .eq. 4) ) ) then
      !  if (.not. (     ((R_tmp(count  ) .lt. ES%R_xpoint(1)) .and. (R_tmp(count-1) .gt. ES%R_xpoint(1))) &
      !             .or. ((R_tmp(count-1) .lt. ES%R_xpoint(1)) .and. (R_tmp(count  ) .gt. ES%R_xpoint(1))) ) ) then
      !    count = count - 2
      !  endif
      !endif
      !if ( (xcase .eq. UPPER_XPOINT) .and. ( (count .eq. 2) .or. (count .eq. 4) ) ) then
      !  if (.not. (     ((R_tmp(count  ) .lt. ES%R_xpoint(2)) .and. (R_tmp(count-1) .gt. ES%R_xpoint(2))) &
      !             .or. ((R_tmp(count-1) .lt. ES%R_xpoint(2)) .and. (R_tmp(count  ) .gt. ES%R_xpoint(2))) ) ) then
      !    count = count - 2
      !  endif
      !endif
    enddo
    if (count .ne. 2) then
      write(*,*)'Something wrong, we are supposed to find 2 secondary strikes points.',count
      write(*,*)'Aborting...'
      stop
    endif
    if (Z_tmp(1) .lt. ES%Z_axis) then
      if (R_tmp(1) .lt. R_tmp(2)) then
        stpts%RStrike_LowerInnerLeg = R_tmp(1)
        stpts%ZStrike_LowerInnerLeg = Z_tmp(1)
        stpts%RStrike_LowerOuterLeg = R_tmp(2)
        stpts%ZStrike_LowerOuterLeg = Z_tmp(2)
      else
        stpts%RStrike_LowerInnerLeg = R_tmp(2)
        stpts%ZStrike_LowerInnerLeg = Z_tmp(2)
        stpts%RStrike_LowerOuterLeg = R_tmp(1)
        stpts%ZStrike_LowerOuterLeg = Z_tmp(1)
      endif
      if (debug) write(*,*)'FOUND THE MAIN STRIKE POINTS:' 
      if (debug) write(*,*)'Lower Inner',stpts%RStrike_LowerInnerLeg,stpts%ZStrike_LowerInnerLeg
      if (debug) write(*,*)'Lower Outer',stpts%RStrike_LowerOuterLeg,stpts%ZStrike_LowerOuterLeg
    else
      if (R_tmp(1) .lt. R_tmp(2)) then
        stpts%RStrike_UpperInnerLeg = R_tmp(1)
        stpts%ZStrike_UpperInnerLeg = Z_tmp(1)
        stpts%RStrike_UpperOuterLeg = R_tmp(2)
        stpts%ZStrike_UpperOuterLeg = Z_tmp(2)
      else
        stpts%RStrike_UpperInnerLeg = R_tmp(2)
        stpts%ZStrike_UpperInnerLeg = Z_tmp(2)
        stpts%RStrike_UpperOuterLeg = R_tmp(1)
        stpts%ZStrike_UpperOuterLeg = Z_tmp(1)
      endif
      if (debug) write(*,*)'FOUND THE MAIN STRIKE POINTS:' 
      if (debug) write(*,*)'Upper Inner',stpts%RStrike_UpperInnerLeg,stpts%ZStrike_UpperInnerLeg
      if (debug) write(*,*)'Upper Outer',stpts%RStrike_UpperOuterLeg,stpts%ZStrike_UpperOuterLeg
    endif
  endif
endif

! ---------------------------------- Find strike points of 2nd separatrix
if ( (xcase .eq. DOUBLE_NULL) .and. ( ES%active_xpoint .ne. SYMMETRIC_XPOINT ) ) then
  if (debug) write(*,*)'looking for strike points of secondary separatrix on main target'
  i_surf = n_flux + n_open
  count = 0
  do i_part=1,flux_list%flux_surfaces(i_surf)%n_parts
    edge_piece(1) = flux_list%flux_surfaces(i_surf)%parts_index(i_part)
    edge_side(1)  = 1
    edge_piece(2) = flux_list%flux_surfaces(i_surf)%parts_index(i_part+1)-1
    edge_side(2)  = 3
    do l=1,2
      rr1   = flux_list%flux_surfaces(i_surf)%s(edge_side(l),edge_piece(l))
      ss1   = flux_list%flux_surfaces(i_surf)%t(edge_side(l),edge_piece(l))
      i_elm = flux_list%flux_surfaces(i_surf)%elm(edge_piece(l))
      call interp_RZ(node_list,element_list,i_elm,rr1,ss1,RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss, &
                                                          ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)
      if ( (ES%active_xpoint .eq. LOWER_XPOINT) .and. (ZZg1 .lt. ES%Z_axis) ) then
        count = count + 1
        R_tmp(count) = RRg1
        Z_tmp(count) = ZZg1
      endif
      if ( (ES%active_xpoint .eq. UPPER_XPOINT) .and. (ZZg1 .gt. ES%Z_axis) ) then
        count = count + 1
        R_tmp(count) = RRg1
        Z_tmp(count) = ZZg1
      endif
    enddo
  enddo
  if (count .ne. 2) then
    write(*,*)'Something wrong, we are supposed to find 2 secondary strikes on main target.',count
    write(*,*)'Aborting...'
    stop
  endif
  if (R_tmp(1) .lt. R_tmp(2)) then
    stpts%RSecondStrike_InnerLeg = R_tmp(1)
    stpts%ZSecondStrike_InnerLeg = Z_tmp(1)
    stpts%RSecondStrike_OuterLeg = R_tmp(2)
    stpts%ZSecondStrike_OuterLeg = Z_tmp(2)
  else
    stpts%RSecondStrike_InnerLeg = R_tmp(2)
    stpts%ZSecondStrike_InnerLeg = Z_tmp(2)
    stpts%RSecondStrike_OuterLeg = R_tmp(1)
    stpts%ZSecondStrike_OuterLeg = Z_tmp(1)
  endif
  if (debug) write(*,*)'FOUND THE SECONDARY STRIKE POINTS:' 
  if (debug) write(*,*)'Inner',stpts%RSecondStrike_InnerLeg,stpts%ZSecondStrike_InnerLeg
  if (debug) write(*,*)'Outer',stpts%RSecondStrike_OuterLeg,stpts%ZSecondStrike_OuterLeg
endif


! ---------------------------------- The last open flux surface (outer SOL boundary)
if (xcase .eq. DOUBLE_NULL) then
  if (debug) write(*,*)'looking for outer SOL limits'
  i_surf = n_flux+n_open+n_outer  
  do i=1,flux_list%flux_surfaces(i_surf)%n_parts
    distance_min = 1.d10
    do k=flux_list%flux_surfaces(i_surf)%parts_index(i),flux_list%flux_surfaces(i_surf)%parts_index(i+1)-1
      rr1   = flux_list%flux_surfaces(i_surf)%s(1,k)
      ss1   = flux_list%flux_surfaces(i_surf)%t(1,k)
      i_elm = flux_list%flux_surfaces(i_surf)%elm(k)
      call interp_RZ(node_list,element_list,i_elm,rr1,ss1,RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss, &
                                                          ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)
      if ( (RRg1 .gt. ES%R_axis) .and. (abs(ZZg1-ES%Z_axis) .lt. distance_min) ) then
        distance_min = abs(ZZg1-ES%Z_axis)
        i_part = i
      endif
    enddo
  enddo
  edge_piece(1) = flux_list%flux_surfaces(i_surf)%parts_index(i_part)
  edge_side(1)  = 1
  edge_piece(2) = flux_list%flux_surfaces(i_surf)%parts_index(i_part+1)-1
  edge_side(2)  = 3
  do l=1,2
    rr1   = flux_list%flux_surfaces(i_surf)%s(edge_side(l),edge_piece(l))
    ss1   = flux_list%flux_surfaces(i_surf)%t(edge_side(l),edge_piece(l))
    i_elm = flux_list%flux_surfaces(i_surf)%elm(edge_piece(l))
    call interp_RZ(node_list,element_list,i_elm,rr1,ss1,RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss, &
                                                        ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)
    if (ZZg1 .gt. ES%Z_axis) then
      stpts%RRightCorn_UpperOuterLeg = RRg1
      stpts%ZRightCorn_UpperOuterLeg = ZZg1
    else
      stpts%RRightCorn_LowerOuterLeg = RRg1
      stpts%ZRightCorn_LowerOuterLeg = ZZg1
    endif
  enddo 
  if (debug) write(*,*)'FOUND THE OUTER LIMITS:' 
  if (debug) write(*,*)'Upper',stpts%RRightCorn_UpperOuterLeg,stpts%ZRightCorn_UpperOuterLeg
  if (debug) write(*,*)'Lower',stpts%RRightCorn_LowerOuterLeg,stpts%ZRightCorn_LowerOuterLeg
endif

! ---------------------------------- The last open flux surface (inner SOL boundary)
if (xcase .eq. DOUBLE_NULL) then
  if (debug) write(*,*)'looking for inner SOL limits'
  i_surf = n_flux+n_open+n_outer+n_inner
  do i=1,flux_list%flux_surfaces(i_surf)%n_parts
    distance_min = 1.d10
    do k=flux_list%flux_surfaces(i_surf)%parts_index(i),flux_list%flux_surfaces(i_surf)%parts_index(i+1)-1
      rr1   = flux_list%flux_surfaces(i_surf)%s(1,k)
      ss1   = flux_list%flux_surfaces(i_surf)%t(1,k)
      i_elm = flux_list%flux_surfaces(i_surf)%elm(k)
      call interp_RZ(node_list,element_list,i_elm,rr1,ss1,RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss, &
                                                          ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)
      if ( (RRg1 .lt. ES%R_axis) .and. (abs(ZZg1-ES%Z_axis) .lt. distance_min) ) then
        distance_min = abs(ZZg1-ES%Z_axis)
        i_part = i
      endif
    enddo
  enddo
  edge_piece(1) = flux_list%flux_surfaces(i_surf)%parts_index(i_part)
  edge_side(1)  = 1
  edge_piece(2) = flux_list%flux_surfaces(i_surf)%parts_index(i_part+1)-1
  edge_side(2)  = 3
  do l=1,2
    rr1   = flux_list%flux_surfaces(i_surf)%s(edge_side(l),edge_piece(l))
    ss1   = flux_list%flux_surfaces(i_surf)%t(edge_side(l),edge_piece(l))
    i_elm = flux_list%flux_surfaces(i_surf)%elm(edge_piece(l))
    call interp_RZ(node_list,element_list,i_elm,rr1,ss1,RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss, &
                                                        ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)
    if (ZZg1 .gt. ES%Z_axis) then
      stpts%RLeftCorn_UpperInnerLeg = RRg1
      stpts%ZLeftCorn_UpperInnerLeg = ZZg1
    else
      stpts%RLeftCorn_LowerInnerLeg = RRg1
      stpts%ZLeftCorn_LowerInnerLeg = ZZg1
    endif
  enddo 
  if (debug) write(*,*)'FOUND THE INNER LIMITS:' 
  if (debug) write(*,*)'Upper',stpts%RLeftCorn_UpperInnerLeg,stpts%ZLeftCorn_UpperInnerLeg
  if (debug) write(*,*)'Lower',stpts%RLeftCorn_LowerInnerLeg,stpts%ZLeftCorn_LowerInnerLeg
endif


! ---------------------------------- The last open flux surface (SOL boundary single null)
if (xcase .ne. DOUBLE_NULL) then
  if (debug) write(*,*)'looking for SOL limits'
  i_surf = n_flux + n_open
  count = 0
  do i_part=1,flux_list%flux_surfaces(i_surf)%n_parts
    edge_piece(1) = flux_list%flux_surfaces(i_surf)%parts_index(i_part)
    edge_side(1)  = 1
    edge_piece(2) = flux_list%flux_surfaces(i_surf)%parts_index(i_part+1)-1
    edge_side(2)  = 3
    do l=1,2
      rr1   = flux_list%flux_surfaces(i_surf)%s(edge_side(l),edge_piece(l))
      ss1   = flux_list%flux_surfaces(i_surf)%t(edge_side(l),edge_piece(l))
      i_elm = flux_list%flux_surfaces(i_surf)%elm(edge_piece(l))
      call interp_RZ(node_list,element_list,i_elm,rr1,ss1,RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss, &
                                                          ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)
      if ( (xcase .eq. LOWER_XPOINT) .and. (ZZg1 .lt. ES%Z_axis) ) then
        count = count + 1
        R_tmp(count) = RRg1
        Z_tmp(count) = ZZg1
      endif
      if ( (xcase .eq. UPPER_XPOINT) .and. (ZZg1 .gt. ES%Z_axis) ) then
        count = count + 1
        R_tmp(count) = RRg1
        Z_tmp(count) = ZZg1
      endif
    enddo
    !if ( (xcase .eq. LOWER_XPOINT) .and. ( (count .eq. 2) .or. (count .eq. 4) ) ) then
    !  if (.not. (     ((R_tmp(count  ) .lt. ES%R_xpoint(1)) .and. (R_tmp(count-1) .gt. ES%R_xpoint(1))) &
    !             .or. ((R_tmp(count-1) .lt. ES%R_xpoint(1)) .and. (R_tmp(count  ) .gt. ES%R_xpoint(1))) ) ) then
    !    count = count - 2
    !  endif
    !endif
    !if ( (xcase .eq. UPPER_XPOINT) .and. ( (count .eq. 2) .or. (count .eq. 4) ) ) then
    !  if (.not. (     ((R_tmp(count  ) .lt. ES%R_xpoint(2)) .and. (R_tmp(count-1) .gt. ES%R_xpoint(2))) &
    !             .or. ((R_tmp(count-1) .lt. ES%R_xpoint(2)) .and. (R_tmp(count  ) .gt. ES%R_xpoint(2))) ) ) then
    !    count = count - 2
    !  endif
    !endif
  enddo
  if (count .ne. 2) then
    write(*,*)'Something wrong, we are supposed to find 2 SOL boundaries in single null.',count
    write(*,*)'Aborting...'
    stop
  endif
  if (xcase .eq. LOWER_XPOINT) then
    if (R_tmp(1) .lt. R_tmp(2)) then
      stpts%RLeftCorn_LowerInnerLeg  = R_tmp(1)
      stpts%ZLeftCorn_LowerInnerLeg  = Z_tmp(1)
      stpts%RRightCorn_LowerOuterLeg = R_tmp(2)
      stpts%ZRightCorn_LowerOuterLeg = Z_tmp(2)
    else
      stpts%RLeftCorn_LowerInnerLeg  = R_tmp(2)
      stpts%ZLeftCorn_LowerInnerLeg  = Z_tmp(2)
      stpts%RRightCorn_LowerOuterLeg = R_tmp(1)
      stpts%ZRightCorn_LowerOuterLeg = Z_tmp(1)
    endif
    if (debug) write(*,*)'FOUND THE SOL BOUNDARY:' 
    if (debug) write(*,*)'Inner',stpts%RLeftCorn_LowerInnerLeg,stpts%ZLeftCorn_LowerInnerLeg
    if (debug) write(*,*)'Outer',stpts%RRightCorn_LowerOuterLeg,stpts%ZRightCorn_LowerOuterLeg
  else
    if (R_tmp(1) .lt. R_tmp(2)) then
      stpts%RLeftCorn_UpperInnerLeg  = R_tmp(1)
      stpts%ZLeftCorn_UpperInnerLeg  = Z_tmp(1)
      stpts%RRightCorn_UpperOuterLeg = R_tmp(2)
      stpts%ZRightCorn_UpperOuterLeg = Z_tmp(2)
    else
      stpts%RLeftCorn_UpperInnerLeg  = R_tmp(2)
      stpts%ZLeftCorn_UpperInnerLeg  = Z_tmp(2)
      stpts%RRightCorn_UpperOuterLeg = R_tmp(1)
      stpts%ZRightCorn_UpperOuterLeg = Z_tmp(1)
    endif
    if (debug) write(*,*)'FOUND THE SOL BOUNDARY:' 
    if (debug) write(*,*)'Inner',stpts%RLeftCorn_UpperInnerLeg,stpts%ZLeftCorn_UpperInnerLeg
    if (debug) write(*,*)'Outer',stpts%RRightCorn_UpperOuterLeg,stpts%ZRightCorn_UpperOuterLeg
  endif
endif


! ---------------------------------- The last open flux surface (Private boundary) under lower X-point 
if (xcase .ne. UPPER_XPOINT) then
  if (debug) write(*,*)'looking for lower private limits'
  i_surf = n_flux+n_open+n_outer+n_inner+n_private  
  count = 0
  do i_part=1,flux_list%flux_surfaces(i_surf)%n_parts
    edge_piece(1) = flux_list%flux_surfaces(i_surf)%parts_index(i_part)
    edge_side(1)  = 1
    edge_piece(2) = flux_list%flux_surfaces(i_surf)%parts_index(i_part+1)-1
    edge_side(2)  = 3
    do l=1,2
      rr1   = flux_list%flux_surfaces(i_surf)%s(edge_side(l),edge_piece(l))
      ss1   = flux_list%flux_surfaces(i_surf)%t(edge_side(l),edge_piece(l))
      i_elm = flux_list%flux_surfaces(i_surf)%elm(edge_piece(l))
      call interp_RZ(node_list,element_list,i_elm,rr1,ss1,RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss, &
                                                          ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)
      if (ZZg1 .lt. ES%Z_axis) then
        count = count + 1
        R_tmp(count) = RRg1
        Z_tmp(count) = ZZg1
      endif
    enddo
    !if ( (count.gt. 0) .and. (abs(R_tmp(count)-R_tmp(count-1)) .lt. 1.d-7) .and. (abs(Z_tmp(count)-Z_tmp(count-1)) .lt. 1.d-7) ) then
    !  count = count - 2
    !endif
    !if ( (count .eq. 2) .or. (count .eq. 4) ) then
    !  if (.not. (     ((R_tmp(count  ) .lt. ES%R_xpoint(1)) .and. (R_tmp(count-1) .gt. ES%R_xpoint(1))) &
    !             .or. ((R_tmp(count-1) .lt. ES%R_xpoint(1)) .and. (R_tmp(count  ) .gt. ES%R_xpoint(1))) ) ) then
    !    count = count - 2
    !  endif
    !endif
  enddo
  if (count .ne. 2) then
    write(*,*)'Something wrong, we are supposed to find 2 lower private boundaries.',count
    write(*,*)'Aborting...'
    stop
  endif
  if (R_tmp(1) .lt. R_tmp(2)) then
    stpts%RRightCorn_LowerInnerLeg = R_tmp(1)
    stpts%ZRightCorn_LowerInnerLeg = Z_tmp(1)
    stpts%RLeftCorn_LowerOuterLeg  = R_tmp(2)
    stpts%ZLeftCorn_LowerOuterLeg  = Z_tmp(2)
  else
    stpts%RRightCorn_LowerInnerLeg = R_tmp(2)
    stpts%ZRightCorn_LowerInnerLeg = Z_tmp(2)
    stpts%RLeftCorn_LowerOuterLeg  = R_tmp(1)
    stpts%ZLeftCorn_LowerOuterLeg  = Z_tmp(1)
  endif
  if (debug) write(*,*)'FOUND THE LOWER PRIVATE BOUNDARY:' 
  if (debug) write(*,*)'Inner',stpts%RRightCorn_LowerInnerLeg,stpts%ZRightCorn_LowerInnerLeg
  if (debug) write(*,*)'Outer',stpts%RLeftCorn_LowerOuterLeg,stpts%ZLeftCorn_LowerOuterLeg
endif


! ---------------------------------- The last open flux surface (Private boundary) above upper X-point 
if (xcase .ne. LOWER_XPOINT) then
  if (debug) write(*,*)'looking for upper private limits'
  i_surf = n_flux+n_open+n_outer+n_inner+n_private+n_up_priv 
  count = 0
  do i_part=1,flux_list%flux_surfaces(i_surf)%n_parts
    edge_piece(1) = flux_list%flux_surfaces(i_surf)%parts_index(i_part)
    edge_side(1)  = 1
    edge_piece(2) = flux_list%flux_surfaces(i_surf)%parts_index(i_part+1)-1
    edge_side(2)  = 3
    do l=1,2
      rr1   = flux_list%flux_surfaces(i_surf)%s(edge_side(l),edge_piece(l))
      ss1   = flux_list%flux_surfaces(i_surf)%t(edge_side(l),edge_piece(l))
      i_elm = flux_list%flux_surfaces(i_surf)%elm(edge_piece(l))
      call interp_RZ(node_list,element_list,i_elm,rr1,ss1,RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss, &
                                                          ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)
      if (ZZg1 .gt. ES%Z_axis) then
        count = count + 1
        R_tmp(count) = RRg1
        Z_tmp(count) = ZZg1
      endif
    enddo
    !if ( (count .ge. 2) .and. (abs(R_tmp(count)-R_tmp(count-1)) .lt. 1.d-7) .and. (abs(Z_tmp(count)-Z_tmp(count-1)) .lt. 1.d-7) ) then
    !  count = count - 2
    !endif
    !if ( (count .eq. 4) .or. (count .eq. 2) ) then
    !  if (.not. (     ((R_tmp(count  ) .lt. ES%R_xpoint(2)) .and. (R_tmp(count-1) .gt. ES%R_xpoint(2))) &
    !             .or. ((R_tmp(count-1) .lt. ES%R_xpoint(2)) .and. (R_tmp(count  ) .gt. ES%R_xpoint(2))) ) ) then
    !    count = count - 2
    !  endif
    !endif
  enddo
  if (count .ne. 2) then
    write(*,*)'Something wrong, we are supposed to find 2 upper private boundaries.',count
    write(*,*)'Aborting...'
    stop
  endif
  if (R_tmp(1) .lt. R_tmp(2)) then
    stpts%RRightCorn_UpperInnerLeg = R_tmp(1)
    stpts%ZRightCorn_UpperInnerLeg = Z_tmp(1)
    stpts%RLeftCorn_UpperOuterLeg  = R_tmp(2)
    stpts%ZLeftCorn_UpperOuterLeg  = Z_tmp(2)
  else
    stpts%RRightCorn_UpperInnerLeg = R_tmp(2)
    stpts%ZRightCorn_UpperInnerLeg = Z_tmp(2)
    stpts%RLeftCorn_UpperOuterLeg  = R_tmp(1)
    stpts%ZLeftCorn_UpperOuterLeg  = Z_tmp(1)
  endif
  if (debug) write(*,*)'FOUND THE UPPER PRIVATE BOUNDARY:' 
  if (debug) write(*,*)'Inner',stpts%RLeftCorn_UpperOuterLeg,stpts%ZRightCorn_UpperInnerLeg
  if (debug) write(*,*)'Outer',stpts%RRightCorn_LowerOuterLeg,stpts%ZLeftCorn_UpperOuterLeg
endif


! ---------------------------------- Find line from axis to lower X-point and get intersection ZMiddle_LowerPrivate with last private surface
if (xcase .ne. UPPER_XPOINT) then
  if (debug) write(*,*)'looking for lower axis-xpoint line'
  tht_x = atan2(ES%Z_xpoint(1)-ES%Z_axis,ES%R_xpoint(1)-ES%R_axis)
  if (tht_x .lt. 0.d0) tht_x = tht_x + 2.d0 * PI
  i_surf = n_flux+n_open+n_outer+n_inner+n_private
  call find_theta_surface(node_list,element_list,flux_list,i_surf,tht_x,ES%R_axis,ES%Z_axis,i_elm_find,s_find,t_find,i_find)

  do i=1,i_find
    call interp_RZ(node_list,element_list,i_elm_find(i),s_find(i),t_find(i),RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss, &
                                                                            ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)
    stpts%RMiddle_LowerPrivate = RRg1
    stpts%ZMiddle_LowerPrivate = ZZg1
    if (stpts%ZMiddle_LowerPrivate .le. ES%Z_xpoint(1)) exit
  enddo
  if (debug) write(*,*)'FOUND THE LOWER PRIVATE MIDDLE:'
  if (debug) write(*,*)'middle:',stpts%RMiddle_LowerPrivate,stpts%ZMiddle_LowerPrivate
endif


! ---------------------------------- Find line from axis to upper X-point and get intersection ZMiddle_UpperPrivate with last private surface
if (xcase .ne. LOWER_XPOINT) then
  if (debug) write(*,*)'looking for upper axis-xpoint line'
  tht_x = atan2(ES%Z_xpoint(2)-ES%Z_axis,ES%R_xpoint(2)-ES%R_axis)
  if (tht_x .lt. 0.d0) tht_x = tht_x + 2.d0 * PI
  i_surf = n_flux+n_open+n_outer+n_inner+n_private+n_up_priv
  call find_theta_surface(node_list,element_list,flux_list,i_surf,tht_x,ES%R_axis,ES%Z_axis,i_elm_find,s_find,t_find,i_find)

  do i=1,i_find
    call interp_RZ(node_list,element_list,i_elm_find(i),s_find(i),t_find(i),RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss, &
                                                                            ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)
    stpts%RMiddle_UpperPrivate = RRg1
    stpts%ZMiddle_UpperPrivate = ZZg1
    if (stpts%ZMiddle_UpperPrivate .ge. ES%Z_xpoint(2)) exit
  enddo
  if (debug) write(*,*)'FOUND THE UPPER PRIVATE MIDDLE:'
  if (debug) write(*,*)'middle:',stpts%RMiddle_UpperPrivate,stpts%ZMiddle_UpperPrivate
endif


!----------------------------------- Define the lines separating the central and upper/lower parts of the grid
if (xcase .ne. UPPER_XPOINT) tht_x1 = atan2(ES%Z_xpoint(1)-ES%Z_axis,ES%R_xpoint(1)-ES%R_axis)
if (xcase .eq. UPPER_XPOINT) tht_x1 = atan2(ES%Z_xpoint(2)-ES%Z_axis,ES%R_xpoint(2)-ES%R_axis)
if (xcase .eq. DOUBLE_NULL ) tht_x2 = atan2(ES%Z_xpoint(2)-ES%Z_axis,ES%R_xpoint(2)-ES%R_axis)
if (tht_x1 .lt. 0.d0) tht_x1 = tht_x1 + 2.d0*PI
if (tht_x2 .lt. 0.d0) tht_x2 = tht_x2 + 2.d0*PI

! --- First define the lower angle
if(xcase .ne. UPPER_XPOINT) then
  if (debug) write(*,*)'looking for angle of lower horizontal line'
  angle_LowerCorner       = atan2(stpts%ZLeftCorn_LowerInnerLeg-ES%Z_xpoint(1),stpts%RLeftCorn_LowerInnerLeg-ES%R_xpoint(1))
  stpts%angle_LowerLeft   = tht_x1 + 1.5d0*PI
  stpts%angle_LowerRight  = tht_x1 + 0.5d0*PI
  ! --- Depending on the equilibrium, it may be better to have a horizontal line... (eg. near double-null at JET)
  if(force_horizontal_Xline) then
    stpts%angle_LowerLeft   = PI
    stpts%angle_LowerRight  = 0.d0
  endif
  if (angle_LowerCorner .lt. 0.d0)          angle_LowerCorner       = angle_LowerCorner       + 2.d0*PI
  if (stpts%angle_LowerLeft   .gt. 2.d0*PI) stpts%angle_LowerLeft   = stpts%angle_LowerLeft   - 2.d0*PI
  if (stpts%angle_LowerRight  .gt. 2.d0*PI) stpts%angle_LowerRight  = stpts%angle_LowerRight  - 2.d0*PI
  
  if (stpts%angle_LowerLeft .gt. angle_LowerCorner) then   ! check if Limit_LowerInnerLeg is above ZLeftCorn_LowerInnerLeg : if not adjust angle
    stpts%angle_LowerLeft  = angle_LowerCorner - 0.1
    stpts%angle_LowerRight = stpts%angle_LowerLeft - PI
    if (stpts%angle_LowerRight  .lt. 0.d0) stpts%angle_LowerRight  = stpts%angle_LowerRight + 2.d0*PI
  endif
  if (debug) write(*,*)'FOUND ANGLE OF THE LOWER LINE:'
  if (debug) write(*,*)'angle (left):',stpts%angle_LowerLeft
endif

! --- Define the upper angle
if(xcase .ne. LOWER_XPOINT) then
  if (debug) write(*,*)'looking for angle of upper horizontal line'
  angle_UpperCorner       = atan2(stpts%ZLeftCorn_UpperInnerLeg-ES%Z_xpoint(2),stpts%RLeftCorn_UpperInnerLeg-ES%R_xpoint(2))
  stpts%angle_UpperLeft   = tht_x2 + 0.5d0*PI; if(xcase .eq. UPPER_XPOINT) stpts%angle_UpperLeft   = tht_x1 + 0.5d0*PI
  stpts%angle_UpperRight  = tht_x2 + 1.5d0*PI; if(xcase .eq. UPPER_XPOINT) stpts%angle_UpperRight  = tht_x1 + 1.5d0*PI
  ! --- Depending on the equilibrium, it may be better to have a horizontal line... (eg. near double-null at JET)
  if(force_horizontal_Xline) then
    stpts%angle_UpperLeft   = PI
    stpts%angle_UpperRight  = 0.d0
  endif
  if (angle_UpperCorner .lt. 0.d0)          angle_UpperCorner       = angle_UpperCorner + 2.d0*PI
  if (stpts%angle_UpperLeft   .gt. 2.d0*PI) stpts%angle_UpperLeft   = stpts%angle_UpperLeft   - 2.d0*PI
  if (stpts%angle_UpperRight  .gt. 2.d0*PI) stpts%angle_UpperRight  = stpts%angle_UpperRight  - 2.d0*PI
  
  if (stpts%angle_UpperLeft .lt. angle_UpperCorner) then   ! check if Limit_UpperInnerLeg is below ZLeftCorn_UpperInnerLeg : if not adjust angle
    stpts%angle_UpperLeft  = angle_UpperCorner + 0.1
    stpts%angle_UpperRight = stpts%angle_UpperLeft + PI
    if (stpts%angle_UpperRight  .gt. 2.d0*PI) stpts%angle_UpperRight  = stpts%angle_UpperRight  - 2.d0*PI
  endif
  if (debug) write(*,*)'FOUND ANGLE OF THE UPPER LINE:'
  if (debug) write(*,*)'angle (left):',stpts%angle_UpperLeft
endif

! --- Then find the lower line
if(xcase .ne. UPPER_XPOINT) then
  if (debug) write(*,*)'looking for lower leg limits'
  i_max = n_flux + n_open + n_outer
  call find_theta_surface(node_list,element_list,flux_list,i_max,stpts%angle_LowerRight,ES%R_xpoint(1),ES%Z_xpoint(1),i_elm_find,s_find,t_find,i_find)
  if(i_find .eq. 0) return
  call interp_RZ(node_list,element_list,i_elm_find(1),s_find(1),t_find(1),&
                 RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss,    &
                 ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)

  stpts%RLimit_LowerOuterLeg = RRg1
  stpts%ZLimit_LowerOuterLeg = ZZg1

  i_max = n_flux + n_open + n_outer + n_inner
  call find_theta_surface(node_list,element_list,flux_list,i_max,stpts%angle_LowerLeft,ES%R_xpoint(1),ES%Z_xpoint(1),i_elm_find,s_find,t_find,i_find)
  if(i_find .eq. 0) return
  call interp_RZ(node_list,element_list,i_elm_find(1),s_find(1),t_find(1),&
                 RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss,    &
                 ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)

  stpts%RLimit_LowerInnerLeg = RRg1
  stpts%ZLimit_LowerInnerLeg = ZZg1
  if (debug) write(*,*)'FOUND LOWER LEG LIMITS:'
  if (debug) write(*,*)'inner:',stpts%RLimit_LowerInnerLeg,stpts%ZLimit_LowerInnerLeg
  if (debug) write(*,*)'outer:',stpts%RLimit_LowerOuterLeg,stpts%ZLimit_LowerOuterLeg
endif

! --- And the upper line
if(xcase .ne. LOWER_XPOINT) then
  if (debug) write(*,*)'looking for upper leg limits'
  i_max = n_flux + n_open + n_outer
  call find_theta_surface(node_list,element_list,flux_list,i_max,stpts%angle_UpperRight,ES%R_xpoint(2),ES%Z_xpoint(2),i_elm_find,s_find,t_find,i_find)
  if(i_find .eq. 0) return
  call interp_RZ(node_list,element_list,i_elm_find(1),s_find(1),t_find(1),&
                 RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss,    &
                 ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)

  stpts%RLimit_UpperOuterLeg = RRg1
  stpts%ZLimit_UpperOuterLeg = ZZg1

  i_max = n_flux + n_open + n_outer + n_inner
  call find_theta_surface(node_list,element_list,flux_list,i_max,stpts%angle_UpperLeft,ES%R_xpoint(2),ES%Z_xpoint(2),i_elm_find,s_find,t_find,i_find)
  if(i_find .eq. 0) return
  call interp_RZ(node_list,element_list,i_elm_find(1),s_find(1),t_find(1),&
                 RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss,    &
                 ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)

  stpts%RLimit_UpperInnerLeg = RRg1
  stpts%ZLimit_UpperInnerLeg = ZZg1
  if (debug) write(*,*)'FOUND UPPER LEG LIMITS:'
  if (debug) write(*,*)'inner:',stpts%RLimit_UpperInnerLeg,stpts%ZLimit_UpperInnerLeg
  if (debug) write(*,*)'outer:',stpts%RLimit_UpperOuterLeg,stpts%ZLimit_UpperOuterLeg
endif



!--------------------------------------------------------------------------------------------
!------------------- And print the output ---------------------------------------------------
!--------------------------------------------------------------------------------------------

write(*,'(A)')                  ' _________________________________________________________'
write(*,'(A)')                  '|                                                         |'
write(*,'(A)')                  '| LEG POINTS (R,Z)                                        |'
write(*,'(A)')                  '|_________________________________________________________|'
write(*,'(A)')                  '|                                                         |'

if (xcase .ne. UPPER_XPOINT) then
  write(*,'(A)')                '| Lower Legs : -------------------------------------------|'
  write(*,'(A,F5.2,A,F5.2,A)')  '|   Upper corner of the Lower Inner Leg  : (',stpts%RLimit_LowerInnerLeg,     ', ', stpts%ZLimit_LowerInnerLeg,     ') |'
  write(*,'(A,F5.2,A,F5.2,A)')  '|   Left corner  of the Lower Inner Leg  : (',stpts%RLeftCorn_LowerInnerLeg,  ', ', stpts%ZLeftCorn_LowerInnerLeg,  ') |'
  write(*,'(A,F5.2,A,F5.2,A)')  '|   Strike point of the Lower Inner Leg  : (',stpts%RStrike_LowerInnerLeg,    ', ', stpts%ZStrike_LowerInnerLeg,    ') |'
  write(*,'(A,F5.2,A,F5.2,A)')  '|   Right corner of the Lower Inner Leg  : (',stpts%RRightCorn_LowerInnerLeg, ', ', stpts%ZRightCorn_LowerInnerLeg, ') |'
  write(*,'(A,F5.2,A,F5.2,A)')  '|   Middle of the Lower Private surface  : (',stpts%RMiddle_LowerPrivate,     ', ', stpts%ZMiddle_LowerPrivate ,    ') |'
  write(*,'(A,F5.2,A,F5.2,A)')  '|   Left corner  of the Lower Outer Leg  : (',stpts%RLeftCorn_LowerOuterLeg,  ', ', stpts%ZLeftCorn_LowerOuterLeg,  ') |'
  write(*,'(A,F5.2,A,F5.2,A)')  '|   Strike point of the Lower Outer Leg  : (',stpts%RStrike_LowerOuterLeg,    ', ', stpts%ZStrike_LowerOuterLeg,    ') |'
  write(*,'(A,F5.2,A,F5.2,A)')  '|   Right corner of the Lower Outer Leg  : (',stpts%RRightCorn_LowerOuterLeg, ', ', stpts%ZRightCorn_LowerOuterLeg, ') |'
  write(*,'(A,F5.2,A,F5.2,A)')  '|   Upper corner of the Lower Outer Leg  : (',stpts%RLimit_LowerOuterLeg,     ', ', stpts%ZLimit_LowerOuterLeg,     ') |'
endif

if (xcase .ne. LOWER_XPOINT) then
  write(*,'(A)')                '| Upper Legs : -------------------------------------------|'
  write(*,'(A,F5.2,A,F5.2,A)')  '|   Lower corner of the Upper Inner Leg  : (',stpts%RLimit_UpperInnerLeg,     ', ', stpts%ZLimit_UpperInnerLeg,     ') |'
  write(*,'(A,F5.2,A,F5.2,A)')  '|   Left corner  of the Upper Inner Leg  : (',stpts%RLeftCorn_UpperInnerLeg,  ', ', stpts%ZLeftCorn_UpperInnerLeg,  ') |'
  write(*,'(A,F5.2,A,F5.2,A)')  '|   Strike point of the Upper Inner Leg  : (',stpts%RStrike_UpperInnerLeg,    ', ', stpts%ZStrike_UpperInnerLeg,    ') |'
  write(*,'(A,F5.2,A,F5.2,A)')  '|   Right corner of the Upper Inner Leg  : (',stpts%RRightCorn_UpperInnerLeg, ', ', stpts%ZRightCorn_UpperInnerLeg, ') |'
  write(*,'(A,F5.2,A,F5.2,A)')  '|   Middle of the Upper Private surface  : (',stpts%RMiddle_UpperPrivate,     ', ', stpts%ZMiddle_UpperPrivate ,    ') |'
  write(*,'(A,F5.2,A,F5.2,A)')  '|   Left corner  of the Upper Outer Leg  : (',stpts%RLeftCorn_UpperOuterLeg,  ', ', stpts%ZLeftCorn_UpperOuterLeg,  ') |'
  write(*,'(A,F5.2,A,F5.2,A)')  '|   Strike point of the Upper Outer Leg  : (',stpts%RStrike_UpperOuterLeg,    ', ', stpts%ZStrike_UpperOuterLeg,    ') |'
  write(*,'(A,F5.2,A,F5.2,A)')  '|   Right corner of the Upper Outer Leg  : (',stpts%RRightCorn_UpperOuterLeg, ', ', stpts%ZRightCorn_UpperOuterLeg, ') |'
  write(*,'(A,F5.2,A,F5.2,A)')  '|   Lower corner of the Upper Outer Leg  : (',stpts%RLimit_UpperOuterLeg,     ', ', stpts%ZLimit_UpperOuterLeg,     ') |'
endif

if ( (xcase .eq. DOUBLE_NULL) .and. (ES%active_xpoint .ne. SYMMETRIC_XPOINT) ) then
  write(*,'(A)')                '| Secondary Strike Points : ------------------------------|'
  write(*,'(A,F5.2,A,F5.2,A)')  '|   Left  Strike point of 2nd separatrix : (',stpts%RSecondStrike_InnerLeg,   ', ', stpts%ZSecondStrike_InnerLeg,   ') |'
  write(*,'(A,F5.2,A,F5.2,A)')  '|   Right Strike point of 2nd separatrix : (',stpts%RSecondStrike_OuterLeg,   ', ', stpts%ZSecondStrike_OuterLeg,   ') |'
endif

write(*,'(A)')                  '|_________________________________________________________|'


return
end subroutine find_strategic_points_advanced
