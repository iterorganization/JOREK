!> Subroutine defines the new grid_points on a divertor leg from crossing of polar and radial coordinate lines
subroutine define_leg_grid(node_list, element_list, newnode_list, newelement_list, flux_list, &
                           xcase, n_grids, stpts, sigmas, nwpts, which_leg)

use constants
use tr_module 
use data_structure
use grid_xpoint_data
use phys_module, only:   tokamak_device, SDN_threshold
use py_plots_grids
use mod_interp, only: interp_RZ
use equil_info

implicit none

! --- Routine parameters
type (type_surface_list)    , intent(inout) :: flux_list
type (type_node_list)       , intent(inout) :: node_list
type (type_element_list)    , intent(inout) :: element_list
type (type_node_list)       , intent(inout) :: newnode_list
type (type_element_list)    , intent(inout) :: newelement_list
type (type_strategic_points), intent(inout) :: stpts
type (type_new_points)      , intent(inout) :: nwpts
integer,                      intent(inout) :: n_grids(12)
integer,                      intent(in)    :: xcase
real*8,                       intent(in)    :: sigmas(17)
integer,                      intent(in)    :: which_leg ! 1=LowerInner, 2=OuterInner, 3=LowerUpper, 4=OuterUpper

! --- local variables
real*8, allocatable :: delta(:)
real*8, allocatable :: R_polar(:,:,:),Z_polar(:,:,:)
real*8, allocatable :: R_polar_smooth(:,:,:),Z_polar_smooth(:,:,:)
integer             :: i, j, k, index, i_sep, pieces, i_elm
integer             :: n_surf_tot, i_surf_tmp
integer,allocatable :: i_flux(:)
integer             :: n_flux, n_tht,  n_open,   n_outer,   n_inner
integer             :: n_private,   n_up_priv,   n_leg,    n_up_leg,   n_leg_out,    n_up_leg_out
integer             :: ifail, my_id
integer             :: n_xpoint_1, n_xpoint_2
integer             :: n_loop, n_tmp
real*8              :: R_cub1d(4), Z_cub1d(4)
real*8              :: length
real*8              :: R1,dR1_dr,dR1_ds,dR1_drs,dR1_drr,dR1_dss
real*8              :: Z1,dZ1_dr,dZ1_ds,dZ1_drs,dZ1_drr,dZ1_dss
real*8              :: diff_min, diff
integer             :: i_surf(3)
real*8              :: R_beg(3), Z_beg(3)
real*8              :: R_end(3), Z_end(3)
real*8              :: R_beg_tmp, Z_beg_tmp
real*8              :: R_end_tmp, Z_end_tmp
integer             :: i_elm_find(8), i_find
real*8              :: s_find(8), t_find(8)
real*8              :: SIG_theta, sig_tmp
real*8              :: SIG_leg_0, SIG_leg_1
real*8              :: SIG_up_leg_0, SIG_up_leg_1
real*8              :: SIG_0, SIG_1
real*8              :: distance_leg, leg_cut
real*8              :: tht_x, rr1, ss1
real*8              :: R_tmp(4), Z_tmp(4)
real*8              :: R_min_tmp, Z_min_tmp
integer             :: count, i_part
integer             :: edge_piece(4), edge_side(4)
integer             :: n_seg
real*8, allocatable :: seg(:), R_seg(:), Z_seg(:)
real*8, allocatable :: R_seg_surf(:,:), Z_seg_surf(:,:)
real*8, allocatable :: R_seg_smooth(:,:), Z_seg_smooth(:,:)
logical             :: xpoint_surface
character*256       :: plot_filename
character*1         :: char_tmp
logical, parameter  :: plot_grid = .true.


write(*,*) '*****************************************'
write(*,*) '* X-point grid inside wall :            *'
write(*,*) '*****************************************'
write(*,*) '                 Define leg part of grid',which_leg

SIG_theta    = sigmas(2) 
SIG_leg_0    = sigmas(8) ; SIG_leg_1    = sigmas(9) 
SIG_up_leg_0 = sigmas(10); SIG_up_leg_1 = sigmas(11)

n_flux    = n_grids(1);  n_tht        = n_grids(2)
n_open    = n_grids(3);  n_outer      = n_grids(4); n_inner = n_grids(5)
n_private = n_grids(6);  n_up_priv    = n_grids(7)
n_leg     = n_grids(8);  n_up_leg     = n_grids(9)
n_leg_out = n_grids(10); n_up_leg_out = n_grids(11)

nwpts%k_cross = 0


!-------------------------------- Allocate data structures for new nodes and initialize them
newnode_list%n_nodes = 0
newnode_list%n_dof   = 0
do i = 1, n_nodes_max
  newnode_list%node(i)%x           = 0.d0
  newnode_list%node(i)%values      = 0.d0
  newnode_list%node(i)%deltas      = 0.d0
  newnode_list%node(i)%index       = 0
  newnode_list%node(i)%boundary    = 0
  newnode_list%node(i)%parents     = 0
  newnode_list%node(i)%parent_elem = 0
  newnode_list%node(i)%ref_lambda  = 0.d0
  newnode_list%node(i)%ref_mu      = 0.d0
  newnode_list%node(i)%constrained = .false.
end do

!-------------------------------- Allocate data structures for new elements and initialize them
newelement_list%n_elements = 0
do i = 1, n_elements_max
  newelement_list%element(i)%vertex       = 0
  newelement_list%element(i)%neighbours   = 0
  newelement_list%element(i)%size         = 0.d0
  newelement_list%element(i)%father       = 0
  newelement_list%element(i)%n_sons       = 0
  newelement_list%element(i)%n_gen        = 0
  newelement_list%element(i)%sons         = 0
  newelement_list%element(i)%contain_node = 0
  newelement_list%element(i)%nref         = 0
end do



!------------------------------------------------------------------------------------------------------------------------!
!************************************************************************************************************************!
!************************************************************************************************************************!
!*********************************** First part: find extrapolation points  *********************************************!
!************************************************************************************************************************!
!************************************************************************************************************************!
!------------------------------------------------------------------------------------------------------------------------!



!---------------------------------------!
!------- Extrapolation points ----------!
!---------------------------------------!


! --- Determine end points of surfaces

! --- lower or upper leg?
if (which_leg .le. 2) then
  R_beg(1) = stpts%RMiddle_LowerPrivate; Z_beg(1) = stpts%ZMiddle_LowerPrivate
  R_beg(2) = ES%R_xpoint(1);                Z_beg(2) = ES%Z_xpoint(1)
  i_surf(1) = n_flux + n_open + n_outer + n_inner + n_private
  if ( (xcase .ne. DOUBLE_NULL) .or. ( (xcase .eq. DOUBLE_NULL) .and. ( (ES%active_xpoint .eq. LOWER_XPOINT) .or. (ES%active_xpoint .eq. SYMMETRIC_XPOINT) ) ) ) then
    i_surf(2) = n_flux
  else
    i_surf(2) = n_flux + n_open
  endif
  ! --- inner or outer leg?
  if (which_leg .eq. 1) then
    R_end(1) = stpts%RRightCorn_LowerInnerLeg; Z_end(1) = stpts%ZRightCorn_LowerInnerLeg
    R_end(2) = stpts%RStrike_LowerInnerLeg;    Z_end(2) = stpts%ZStrike_LowerInnerLeg
    R_beg(3) = stpts%RLimit_LowerInnerLeg;     Z_beg(3) = stpts%ZLimit_LowerInnerLeg
    R_end(3) = stpts%RLeftCorn_LowerInnerLeg;  Z_end(3) = stpts%ZLeftCorn_LowerInnerLeg
    i_surf(3) = n_flux + n_open + n_outer + n_inner
    n_seg = n_leg
    if ( (xcase .eq. DOUBLE_NULL) .and. ( ES%active_xpoint .eq. UPPER_XPOINT ) ) then
      n_surf_tot = n_inner + n_private + 1
      allocate(i_flux(n_surf_tot))
      do i=1,n_private
        i_flux(i) = n_flux + n_open + n_outer + n_inner + n_private - i + 1
      enddo
      i_flux(n_private+1) = i_surf(2)
      i_sep = n_private+1
      do i=1,n_inner
        i_flux(n_private+1+i) = n_flux + n_open + n_outer + i
      enddo
    else
      n_surf_tot = n_inner + n_open + n_private + 1
      allocate(i_flux(n_surf_tot))
      do i=1,n_private
        i_flux(i) = n_flux + n_open + n_outer + n_inner + n_private - i + 1
      enddo
      i_flux(n_private+1) = i_surf(2)
      i_sep = n_private+1
      do i=1,n_open
        i_flux(n_private+1+i) = n_flux + i
      enddo
      do i=1,n_inner
        i_flux(n_private+n_open+1+i) = n_flux + n_open + n_outer + i
      enddo
    endif
    n_xpoint_1 = 2 ! please see "create_x_node.f90" if confused
    n_xpoint_2 = 1
  else
    R_end(1) = stpts%RLeftCorn_LowerOuterLeg;  Z_end(1) = stpts%ZLeftCorn_LowerOuterLeg
    R_end(2) = stpts%RStrike_LowerOuterLeg;    Z_end(2) = stpts%ZStrike_LowerOuterLeg
    R_beg(3) = stpts%RLimit_LowerOuterLeg;     Z_beg(3) = stpts%ZLimit_LowerOuterLeg
    R_end(3) = stpts%RRightCorn_LowerOuterLeg; Z_end(3) = stpts%ZRightCorn_LowerOuterLeg
    i_surf(3) = n_flux + n_open + n_outer
    n_seg = n_leg
    if (n_leg_out .gt. 0) n_seg = n_leg_out
    if ( (xcase .eq. DOUBLE_NULL) .and. ( ES%active_xpoint .eq. UPPER_XPOINT ) ) then
      n_surf_tot = n_outer + n_private + 1
      allocate(i_flux(n_surf_tot))
      do i=1,n_private
        i_flux(i) = n_flux + n_open + n_outer + n_inner + n_private - i + 1
      enddo
      i_flux(n_private+1) = i_surf(2)
      i_sep = n_private+1
      do i=1,n_outer
        i_flux(n_private+1+i) = n_flux + n_open + i
      enddo
    else
      n_surf_tot = n_outer + n_open + n_private + 1
      allocate(i_flux(n_surf_tot))
      do i=1,n_private
        i_flux(i) = n_flux + n_open + n_outer + n_inner + n_private - i + 1
      enddo
      i_flux(n_private+1) = i_surf(2)
      i_sep = n_private+1
      do i=1,n_open
        i_flux(n_private+1+i) = n_flux + i
      enddo
      do i=1,n_outer
        i_flux(n_private+n_open+1+i) = n_flux + n_open + i
      enddo
    endif
    n_xpoint_1 = 3 ! please see "create_x_node.f90" if confused
    n_xpoint_2 = 4
  endif
  SIG_0 = SIG_leg_0
  SIG_1 = SIG_leg_1
else
  i_surf(1) = n_flux + n_open + n_outer + n_inner + n_private + n_up_priv
  R_beg(1) = stpts%RMiddle_UpperPrivate; Z_beg(1) = stpts%ZMiddle_UpperPrivate
  R_beg(2) = ES%R_xpoint(2);             Z_beg(2) = ES%Z_xpoint(2)
  if ( (xcase .ne. DOUBLE_NULL) .or. ( (xcase .eq. DOUBLE_NULL) .and. ( (ES%active_xpoint .eq. UPPER_XPOINT) .or. (ES%active_xpoint .eq. SYMMETRIC_XPOINT) ) ) ) then
    i_surf(2) = n_flux
  else
    i_surf(2) = n_flux + n_open
  endif
  if (which_leg .eq. 3) then
    R_end(1) = stpts%RRightCorn_UpperInnerLeg; Z_end(1) = stpts%ZRightCorn_UpperInnerLeg
    R_end(2) = stpts%RStrike_UpperInnerLeg;    Z_end(2) = stpts%ZStrike_UpperInnerLeg
    R_beg(3) = stpts%RLimit_UpperInnerLeg;     Z_beg(3) = stpts%ZLimit_UpperInnerLeg
    R_end(3) = stpts%RLeftCorn_UpperInnerLeg;  Z_end(3) = stpts%ZLeftCorn_UpperInnerLeg
    i_surf(3) = n_flux + n_open + n_outer + n_inner
    n_seg = n_up_leg
    if ( ES%active_xpoint .eq. UPPER_XPOINT  ) then
      n_surf_tot = n_inner + n_open + n_up_priv + 1
      allocate(i_flux(n_surf_tot))
      do i=1,n_up_priv
        i_flux(i) = n_flux + n_open + n_outer + n_inner + n_private + n_up_priv - i + 1
      enddo
      i_flux(n_up_priv+1) = i_surf(2)
      i_sep = n_up_priv+1
      do i=1,n_open
        i_flux(n_up_priv+1+i) = n_flux + i
      enddo
      do i=1,n_inner
        i_flux(n_up_priv+n_open+1+i) = n_flux + n_open + n_outer +i
      enddo
    else
      n_surf_tot = n_inner + n_up_priv + 1
      allocate(i_flux(n_surf_tot))
      do i=1,n_up_priv
        i_flux(i) = n_flux + n_open + n_outer + n_inner + n_private + n_up_priv - i + 1
      enddo
      i_flux(n_up_priv+1) = i_surf(2)
      i_sep = n_up_priv+1
      do i=1,n_inner
        i_flux(n_up_priv+1+i) = n_flux + n_open + n_outer +i
      enddo
    endif
    ! --- please see "create_x_node.f90" if confused
    if (xcase .eq. UPPER_XPOINT) then
      n_xpoint_1 = 3
      n_xpoint_2 = 4
    else
      n_xpoint_1 = 7
      n_xpoint_2 = 8
    endif
  ! --- inner or outer leg?
  else
    R_end(1) = stpts%RLeftCorn_UpperOuterLeg;  Z_end(1) = stpts%ZLeftCorn_UpperOuterLeg
    R_end(2) = stpts%RStrike_UpperOuterLeg;    Z_end(2) = stpts%ZStrike_UpperOuterLeg
    R_beg(3) = stpts%RLimit_UpperOuterLeg;     Z_beg(3) = stpts%ZLimit_UpperOuterLeg
    R_end(3) = stpts%RRightCorn_UpperOuterLeg; Z_end(3) = stpts%ZRightCorn_UpperOuterLeg
    i_surf(3) = n_flux + n_open + n_outer
    n_seg = n_up_leg
    if (n_up_leg_out .gt. 0) n_seg = n_up_leg_out
    if ( ES%active_xpoint .eq. UPPER_XPOINT ) then
      n_surf_tot = n_outer + n_open + n_up_priv + 1
      allocate(i_flux(n_surf_tot))
      do i=1,n_up_priv
        i_flux(i) = n_flux + n_open + n_outer + n_inner + n_private + n_up_priv - i + 1
      enddo
      i_flux(n_up_priv+1) = i_surf(2)
      i_sep = n_up_priv+1
      do i=1,n_open
        i_flux(n_up_priv+1+i) = n_flux + i
      enddo
      do i=1,n_outer
        i_flux(n_up_priv+n_open+1+i) = n_flux + n_open + i
      enddo
    else
      n_surf_tot = n_outer + n_up_priv + 1
      allocate(i_flux(n_surf_tot))
      do i=1,n_up_priv
        i_flux(i) = n_flux + n_open + n_outer + n_inner + n_private + n_up_priv - i + 1
      enddo
      i_flux(n_up_priv+1) = i_surf(2)
      i_sep = n_up_priv+1
      do i=1,n_outer
        i_flux(n_up_priv+1+i) = n_flux + n_open + i
      enddo
    endif
    ! --- please see "create_x_node.f90" if confused
    if (xcase .eq. UPPER_XPOINT) then
      n_xpoint_1 = 2
      n_xpoint_2 = 1
    else
      n_xpoint_1 = 6
      n_xpoint_2 = 5
    endif
  endif
  SIG_0 = SIG_up_leg_0
  SIG_1 = SIG_up_leg_1
endif

! --- Define the segmentation array
allocate(seg(n_seg), R_seg(n_seg), Z_seg(n_seg))
allocate(R_seg_surf(n_surf_tot,n_seg), Z_seg_surf(n_surf_tot,n_seg))
allocate(R_seg_smooth(3,n_seg), Z_seg_smooth(3,n_seg))

! --- We cut the end to ensure we remain inside the domain...
distance_leg = sqrt( (R_end(2)-R_beg(2))**2 + (Z_end(2)-Z_beg(2))**2 )
leg_cut = 1.d0 - 1.d-3 / distance_leg ! we cut by ~1mmm

! --- Segment the private, separatrix and SOL surfaces
do i = 1,n_surf_tot
  seg = 0
  ! --- Change sig near the Xpoint to avoid strong angles
  if (i .lt. i_sep) then
    sig_tmp = SIG_1 + (1.3d0*SIG_1 - SIG_1) * (real(i-1)/real(i_sep-1))**6 ! large exponent to have fast decay...
  else if (i .gt. i_sep) then
    sig_tmp = SIG_1 + (1.3d0*SIG_1 - SIG_1) * (real(n_surf_tot-i)/real(n_surf_tot-i_sep))**6
  else if (i .eq. i_sep) then
    sig_tmp = 1.3d0*SIG_1
  endif
  call meshac2(n_seg,seg,0.d0,1.d0,sig_tmp,SIG_0,0.6d0,1.0d0)
  seg = leg_cut*seg ! we cut the end to ensure we remain inside the domain...
  xpoint_surface = .false.
  if (i .eq. i_sep) xpoint_surface = .true.
  ! --- The end points
  if (i .eq. 1) then
    R_beg_tmp = R_beg(1) ; R_end_tmp = R_end(1)
    Z_beg_tmp = Z_beg(1) ; Z_end_tmp = Z_end(1)
  elseif (i .eq. i_sep) then
    R_beg_tmp = R_beg(2) ; R_end_tmp = R_end(2)
    Z_beg_tmp = Z_beg(2) ; Z_end_tmp = Z_end(2)
  elseif (i .eq. n_surf_tot) then
    R_beg_tmp = R_beg(3) ; R_end_tmp = R_end(3)
    Z_beg_tmp = Z_beg(3) ; Z_end_tmp = Z_end(3)
  else
    ! --- The Beg point
    if (i .lt. i_sep) then
      tht_x = atan2(Z_beg(1)-Z_beg(2),R_beg(1)-R_beg(2))
    else
      tht_x = atan2(Z_beg(3)-Z_beg(2),R_beg(3)-R_beg(2))
    endif
    if (tht_x .lt. 0.d0) tht_x = tht_x + 2.d0 * PI
    call find_theta_surface(node_list,element_list,flux_list,i_flux(i),tht_x,R_beg(2),Z_beg(2),i_elm_find,s_find,t_find,i_find)
    do k=1,i_find
      call interp_RZ(node_list,element_list,i_elm_find(k),s_find(k),t_find(k),R1,dR1_dr,dR1_ds,dR1_drs,dR1_drr,dR1_dss, &
                                                                              Z1,dZ1_dr,dZ1_ds,dZ1_drs,dZ1_drr,dZ1_dss)
      R_beg_tmp = R1
      Z_beg_tmp = Z1
      if ( (Z_beg_tmp .le. ES%Z_xpoint(1)) .and. (which_leg .le. 2) .and. (i .lt. i_sep) ) exit
      if ( (Z_beg_tmp .ge. ES%Z_xpoint(2)) .and. (which_leg .gt. 2) .and. (i .lt. i_sep) ) exit
      if ( (R_beg_tmp .le. ES%R_xpoint(1)) .and. (which_leg .eq. 1) .and. (i .gt. i_sep) ) exit
      if ( (R_beg_tmp .ge. ES%R_xpoint(1)) .and. (which_leg .eq. 2) .and. (i .gt. i_sep) ) exit
      if ( (R_beg_tmp .le. ES%R_xpoint(2)) .and. (which_leg .eq. 3) .and. (i .gt. i_sep) ) exit
      if ( (R_beg_tmp .ge. ES%R_xpoint(2)) .and. (which_leg .eq. 4) .and. (i .gt. i_sep) ) exit
    enddo
    ! --- The End point
    count = 0
    do i_part=1,flux_list%flux_surfaces(i_flux(i))%n_parts
      edge_piece(1) = flux_list%flux_surfaces(i_flux(i))%parts_index(i_part)
      edge_side(1)  = 1
      edge_piece(2) = flux_list%flux_surfaces(i_flux(i))%parts_index(i_part+1)-1
      edge_side(2)  = 3
      do k=1,2
        rr1   = flux_list%flux_surfaces(i_flux(i))%s(edge_side(k),edge_piece(k))
        ss1   = flux_list%flux_surfaces(i_flux(i))%t(edge_side(k),edge_piece(k))
        i_elm = flux_list%flux_surfaces(i_flux(i))%elm(edge_piece(k))
        call interp_RZ(node_list,element_list,i_elm,rr1,ss1,R1,dR1_dr,dR1_ds,dR1_drs,dR1_drr,dR1_dss, &
                                                            Z1,dZ1_dr,dZ1_ds,dZ1_drs,dZ1_drr,dZ1_dss)
        if (which_leg .eq. 1) then
          !if ( (R1 .lt. ES%R_xpoint(1)) .and. (Z1 .lt. ES%Z_axis) ) then
          if (Z1 .lt. ES%Z_axis) then
            count = count + 1
            R_tmp(count) = R1
            Z_tmp(count) = Z1
          endif
        endif
        if (which_leg .eq. 2) then
          !if ( (R1 .gt. ES%R_xpoint(1)) .and. (Z1 .lt. ES%Z_axis) ) then
          if (Z1 .lt. ES%Z_axis) then
            count = count + 1
            R_tmp(count) = R1
            Z_tmp(count) = Z1
          endif
        endif
        if (which_leg .eq. 3) then
          !if ( (R1 .lt. ES%R_xpoint(2)) .and. (Z1 .gt. ES%Z_axis) ) then
          if (Z1 .gt. ES%Z_axis) then
            count = count + 1
            R_tmp(count) = R1
            Z_tmp(count) = Z1
          endif
        endif
        if (which_leg .eq. 4) then
          !if ( (R1 .gt. ES%R_xpoint(2)) .and. (Z1 .gt. ES%Z_axis) ) then
          if (Z1 .gt. ES%Z_axis) then
            count = count + 1
            R_tmp(count) = R1
            Z_tmp(count) = Z1
          endif
        endif
      enddo
      if (count .ge. 2) then
        if ( (abs(R_tmp(count)-R_tmp(count-1)) .lt. 1.d-7) .and. (abs(Z_tmp(count)-Z_tmp(count-1)) .lt. 1.d-7) ) then
          count = count - 2
        endif
      endif
    enddo
    if (count .ne. 1) then
      if ( (which_leg .eq. 1) .or. (which_leg .eq. 3) ) then
        if (R_tmp(1) .gt. R_tmp(2)) then
          R_tmp(1) = R_tmp(2)
          Z_tmp(1) = Z_tmp(2)
        endif
      else
        if (R_tmp(1) .lt. R_tmp(2)) then
          R_tmp(1) = R_tmp(2)
          Z_tmp(1) = Z_tmp(2)
        endif
      endif
      count = count - 1
    endif
    if (count .ne. 1) then
      diff_min = 1.d10
      do j=1,count
        if (which_leg .le. 2) then
          diff = sqrt( (R_tmp(j)-ES%R_xpoint(1))**2 + (Z_tmp(j)-ES%Z_xpoint(1))**2 )
        else
          diff = sqrt( (R_tmp(j)-ES%R_xpoint(2))**2 + (Z_tmp(j)-ES%Z_xpoint(2))**2 )
        endif
        if (diff .lt. diff_min) then
          diff_min = diff
          R_min_tmp = R_tmp(j)
          Z_min_tmp = Z_tmp(j)
        endif
      enddo
      R_tmp(1) = R_min_tmp
      Z_tmp(1) = Z_min_tmp
    endif
    R_end_tmp = R_tmp(1)
    Z_end_tmp = Z_tmp(1)
  endif
  ! --- Segment surface
  call segment_surface_length(node_list,element_list,flux_list%flux_surfaces(i_flux(i)), &
                              R_beg_tmp, Z_beg_tmp, R_end_tmp, Z_end_tmp, n_seg, seg, R_seg, Z_seg, xpoint_surface)
  do j = 1,n_seg
    R_seg_surf(i,j) = R_seg(j)
    Z_seg_surf(i,j) = Z_seg(j)
  enddo
enddo

! --- We know where the X-point is, it's safer to use it directly
R_seg_surf(1,1)          = R_beg(1)
Z_seg_surf(1,1)          = Z_beg(1)
R_seg_surf(i_sep,1)      = R_beg(2)
Z_seg_surf(i_sep,1)      = Z_beg(2)
R_seg_surf(n_surf_tot,1) = R_beg(3)
Z_seg_surf(n_surf_tot,1) = Z_beg(3)

! --- Using polar at every node is safer if target is uneven, but it looks nicer with simpler polar
do i = 1,3
  if (i .eq. 1) i_surf_tmp = 1
  if (i .eq. 2) i_surf_tmp = i_sep
  if (i .eq. 3) i_surf_tmp = n_surf_tot
  do j = 1,n_seg
    R_seg_smooth(i,j) = R_seg_surf(i_surf_tmp,j)
    Z_seg_smooth(i,j) = Z_seg_surf(i_surf_tmp,j)
  enddo
enddo



!----------------------------------- Print a python file that plots a cross with the 4 nodes of each element
if (plot_grid) then
  open(101,file='plot_leg.py')
    write(101,'(A)')                '#!/usr/bin/env python'
    write(101,'(A)')                'import numpy as N'
    write(101,'(A)')                'import pylab'
    write(101,'(A)')                'def main():'
    do k=1,n_surf_tot
      write(101,'(A,i6,A)')         ' r = N.zeros(',n_seg,')'
      write(101,'(A,i6,A)')         ' z = N.zeros(',n_seg,')'
      do i=1,n_seg
        write(101,'(A,i6,A,f15.4)') ' r[',i-1,'] = ',R_seg_surf(k,i)
        write(101,'(A,i6,A,f15.4)') ' z[',i-1,'] = ',Z_seg_surf(k,i)
      enddo
      write(101,'(A)')              ' pylab.plot(r,z, "b-x")'
      write(101,'(A)')              ' pylab.plot(r[0],z[0], "rx")'
    enddo
    write(101,'(A)')                ' pylab.axis("equal")'
    write(101,'(A)')                ' pylab.show()'
    write(101,'(A)')                ' '
    write(101,'(A)')                'main()'
  close(101)
endif



if (plot_grid) then
  plot_filename = 'plot_leg_surf.py'
  call print_py_plot_prepare_plot(plot_filename)
  open(101,file=plot_filename,position='append')
  do k=1,n_surf_tot
    write(101,'(A,i6,A)')         ' r = N.zeros(',n_seg,')'
    write(101,'(A,i6,A)')         ' z = N.zeros(',n_seg,')'
    do i=1,n_seg
      write(101,'(A,i6,A,f15.4)') ' r[',i-1,'] = ',R_seg_surf(k,i)
      write(101,'(A,i6,A,f15.4)') ' z[',i-1,'] = ',Z_seg_surf(k,i)
    enddo
    write(101,'(A)')              ' pylab.plot(r,z, "b-x")'
    write(101,'(A)')              ' pylab.plot(r[0],z[0], "rx")'
  enddo
  close(101)
  call print_py_plot_ordered_flux_surfaces(plot_filename, node_list, element_list, flux_list, 'r', .false.)
  call print_py_plot_wall(plot_filename)
  call print_py_plot_finish_plot(plot_filename)
endif


!------------------------------------------------------------------------------------------------------------------------!
!************************************************************************************************************************!
!************************************************************************************************************************!
!*********************************** Second part: find crossings of lines ***********************************************!
!************************************************************************************************************************!
!************************************************************************************************************************!
!------------------------------------------------------------------------------------------------------------------------!
write(*,*) '                 Find crossings between coordinate lines'





!--------------------------------------------------------------------------!
!------------------- Construct polar coordinate lines ---------------------!
!--------------------------------------------------------------------------!

allocate (delta(n_seg))
allocate (R_polar_smooth(3,4,n_seg),Z_polar_smooth(3,4,n_seg))
allocate (R_polar(n_surf_tot-1,4,n_seg),Z_polar(n_surf_tot-1,4,n_seg))

! --- The safe polar coords
do i=1,n_seg
  call create_polar_lines_simple(n_surf_tot, R_seg_surf(1:n_surf_tot,i), Z_seg_surf(1:n_surf_tot,i), R_polar(1:n_surf_tot-1,1:4,i) , Z_polar(1:n_surf_tot-1,1:4,i) )
enddo

! --- The smooth polar coords
delta = 0.1
delta(1) = 0.d0  
delta(2) = 0.05d0
call create_polar_lines(n_seg, R_seg_smooth(1,1:n_seg), Z_seg_smooth(1,1:n_seg), &
                               R_seg_smooth(2,1:n_seg), Z_seg_smooth(2,1:n_seg), &
                               R_seg_smooth(3,1:n_seg), Z_seg_smooth(3,1:n_seg), &
                               delta, R_polar_smooth, Z_polar_smooth)
deallocate (delta)

do i=1,3
  do j=1,n_seg
    nwpts%R_polar(i,1:4,j) = R_polar_smooth(i,1:4,j)
    nwpts%Z_polar(i,1:4,j) = Z_polar_smooth(i,1:4,j)
  enddo
enddo


!--------------------------------------------------------------------------!
!--------- Find grid_points from crossing of coordinate lines -------------!
!--------------------------------------------------------------------------!


!----------------------------------- The main part (without legs, outer and inner)
do i=1,n_surf_tot
  i_surf_tmp = i_flux(i)
  do j=1, n_seg
    
    ! --- First try with the smooth polar lines
    do k=1,3
      nwpts%R_polar(k,1:4,j) = R_polar_smooth(k,1:4,j)
      nwpts%Z_polar(k,1:4,j) = Z_polar_smooth(k,1:4,j)
    enddo
    do k=1,3     ! 3 line pieces per coordinate line

      call from_polar_to_cubic(R_polar_smooth(k,1:4,j),R_cub1d)
      call from_polar_to_cubic(Z_polar_smooth(k,1:4,j),Z_cub1d)
      call find_crossing(node_list,element_list,flux_list,i_surf_tmp,R_cub1d,Z_cub1d, &
                         nwpts%RR_new(i,j),nwpts%ZZ_new(i,j),nwpts%ielm_flux(i,j),&
                         nwpts%s_flux(i,j),nwpts%t_flux(i,j),nwpts%t_tht(i,j),ifail, .true.)
      ! --- Readjust to make sure we are inside element.
      if ( (nwpts%s_flux(i,j) .lt. 0.d0) .or. (nwpts%s_flux(i,j) .gt. 1.d0) .or. (nwpts%t_flux(i,j) .lt. 0.d0) .or. (nwpts%t_flux(i,j) .gt. 1.d0) ) then
        if (nwpts%s_flux(i,j) .lt. 0.d0) nwpts%s_flux(i,j) = 0.d0
        if (nwpts%s_flux(i,j) .gt. 1.d0) nwpts%s_flux(i,j) = 1.d0
        if (nwpts%t_flux(i,j) .lt. 0.d0) nwpts%t_flux(i,j) = 0.d0
        if (nwpts%t_flux(i,j) .gt. 1.d0) nwpts%t_flux(i,j) = 1.d0
        call interp_RZ(node_list,element_list,nwpts%ielm_flux(i,j),nwpts%s_flux(i,j),nwpts%t_flux(i,j),&
                       R1,dR1_dr,dR1_ds,dR1_drs,dR1_drr,dR1_dss,    &
                       Z1,dZ1_dr,dZ1_ds,dZ1_drs,dZ1_drr,dZ1_dss)
        nwpts%RR_new(i,j) = R1
        nwpts%ZZ_new(i,j) = Z1
      endif
      
      if (ifail .eq. 0) then
        nwpts%k_cross(i,j) = k
        exit
      endif

    enddo

    ! --- If the smooth polar lines don't work, try the safe ones
    if (ifail .ne. 0) then
      do k=1,n_surf_tot-1
        nwpts%R_polar(k,1:4,j) = R_polar(k,1:4,j)
        nwpts%Z_polar(k,1:4,j) = Z_polar(k,1:4,j)
      enddo
      do k=1,n_surf_tot-1  ! 3 line pieces per coordinate line

        call from_polar_to_cubic(R_polar(k,1:4,j),R_cub1d)
        call from_polar_to_cubic(Z_polar(k,1:4,j),Z_cub1d)
        call find_crossing(node_list,element_list,flux_list,i_surf_tmp,R_cub1d,Z_cub1d, &
                           nwpts%RR_new(i,j),nwpts%ZZ_new(i,j),nwpts%ielm_flux(i,j),&
                           nwpts%s_flux(i,j),nwpts%t_flux(i,j),nwpts%t_tht(i,j),ifail,.true.)
        ! --- Readjust to make sure we are inside element.
        if ( (nwpts%s_flux(i,j) .lt. 0.d0) .or. (nwpts%s_flux(i,j) .gt. 1.d0) .or. (nwpts%t_flux(i,j) .lt. 0.d0) .or. (nwpts%t_flux(i,j) .gt. 1.d0) ) then
          if (nwpts%s_flux(i,j) .lt. 0.d0) nwpts%s_flux(i,j) = 0.d0
          if (nwpts%s_flux(i,j) .gt. 1.d0) nwpts%s_flux(i,j) = 1.d0
          if (nwpts%t_flux(i,j) .lt. 0.d0) nwpts%t_flux(i,j) = 0.d0
          if (nwpts%t_flux(i,j) .gt. 1.d0) nwpts%t_flux(i,j) = 1.d0
          call interp_RZ(node_list,element_list,nwpts%ielm_flux(i,j),nwpts%s_flux(i,j),nwpts%t_flux(i,j),&
                         R1,dR1_dr,dR1_ds,dR1_drs,dR1_drr,dR1_dss,    &
                         Z1,dZ1_dr,dZ1_ds,dZ1_drs,dZ1_drr,dZ1_dss)
          nwpts%RR_new(i,j) = R1
          nwpts%ZZ_new(i,j) = Z1
        endif
        
        if (ifail .eq. 0) then
          nwpts%k_cross(i,j) = k
          exit
        endif

      enddo

      if (ifail .ne. 0) then
        write(*,'(A,I6,I6,I6,F20.10)') ' WARNING node not found for leg grid : ',ifail,i,j
      endif
    endif
      
  enddo
enddo


if (plot_grid) then
  open(101,file='plot_leg_nodes.py')
    write(101,'(A)')                '#!/usr/bin/env python'
    write(101,'(A)')                'import numpy as N'
    write(101,'(A)')                'import pylab'
    write(101,'(A)')                'def main():'
    do i=1,n_surf_tot
      do j=1,n_seg
        write(101,'(A,f15.4)') ' r = ',nwpts%RR_new(i,j)
        write(101,'(A,f15.4)') ' z = ',nwpts%ZZ_new(i,j)
        write(101,'(A)')            ' pylab.plot(r,z, "bx")'
      enddo
    enddo
    write(101,'(A)')                ' pylab.axis("equal")'
    write(101,'(A)')                ' pylab.show()'
    write(101,'(A)')                ' '
    write(101,'(A)')                'main()'
  close(101)
endif














!------------------------------------------------------------------------------------------------------------------------!
!************************************************************************************************************************!
!************************************************************************************************************************!
!***************************************** Third part: define the new nodes  ********************************************!
!************************************************************************************************************************!
!************************************************************************************************************************!
!------------------------------------------------------------------------------------------------------------------------!
write(*,*) '                 Defining new nodes'



!-------------------------------------------------------------------------------------------!
!-------------------------------------- The Xpoints ----------------------------------------!
!-------------------------------------------------------------------------------------------!

! THIS ADDS FOUR NODES AT EACH XPOINTS, PLEASE SEE create_x_node FOR MORE DETAILS
if (xcase .eq. LOWER_XPOINT) then
  call create_x_node(node_list, element_list, newnode_list, nwpts, stpts, &
                     LOWER_XPOINT, ES%R_axis, ES%Z_axis, ES%R_xpoint, ES%Z_xpoint, ES%i_elm_xpoint, ES%s_xpoint, ES%t_xpoint)
endif 
if (xcase .eq. UPPER_XPOINT) then
  call create_x_node(node_list, element_list, newnode_list, nwpts, stpts, &
                     UPPER_XPOINT, ES%R_axis, ES%Z_axis, ES%R_xpoint, ES%Z_xpoint, ES%i_elm_xpoint, ES%s_xpoint, ES%t_xpoint)
endif 
if ( (xcase .eq. DOUBLE_NULL) .and. ( (ES%active_xpoint .eq. LOWER_XPOINT) .or. (ES%active_xpoint .eq. SYMMETRIC_XPOINT) ) ) then ! Put lower Xpoint first
  call create_x_node(node_list, element_list, newnode_list, nwpts, stpts, &
                     LOWER_XPOINT, ES%R_axis, ES%Z_axis, ES%R_xpoint, ES%Z_xpoint, ES%i_elm_xpoint, ES%s_xpoint, ES%t_xpoint)
  call create_x_node(node_list, element_list, newnode_list, nwpts, stpts, &
                     UPPER_XPOINT, ES%R_axis, ES%Z_axis, ES%R_xpoint, ES%Z_xpoint, ES%i_elm_xpoint, ES%s_xpoint, ES%t_xpoint)
endif 
if ( (xcase .eq. DOUBLE_NULL) .and. ( ES%active_xpoint .eq. UPPER_XPOINT ) ) then ! Put upper Xpoint first
  call create_x_node(node_list, element_list, newnode_list, nwpts, stpts, &
                     UPPER_XPOINT, ES%R_axis, ES%Z_axis, ES%R_xpoint, ES%Z_xpoint, ES%i_elm_xpoint, ES%s_xpoint, ES%t_xpoint)
  call create_x_node(node_list, element_list, newnode_list, nwpts, stpts, &
                     LOWER_XPOINT, ES%R_axis, ES%Z_axis, ES%R_xpoint, ES%Z_xpoint, ES%i_elm_xpoint, ES%s_xpoint, ES%t_xpoint)
endif 
index = newnode_list%n_nodes



!----------------------------------------------------------------------------!
!-------------------------------- The leg -----------------------------------!
!----------------------------------------------------------------------------!

do i=1,n_surf_tot                 
  do j=1, n_seg

    index = index + 1
    call create_new_node(node_list, element_list, newnode_list, index, i, j, nwpts)

  enddo
enddo
newnode_list%n_nodes = index




!------------------------------------------------------------------------------------------------------------------------!
!************************************************************************************************************************!
!************************************************************************************************************************!
!*************************************** Fourth part: define the new elements  ******************************************!
!************************************************************************************************************************!
!************************************************************************************************************************!
!------------------------------------------------------------------------------------------------------------------------!
write(*,*) '                 Defining new elements'


!-------------------------------- The closed region
n_tmp = 4 ! because we put the Xpoints first
if (xcase .eq. DOUBLE_NULL) n_tmp = 8 ! because we put the Xpoints first
index = 0
do i=1,n_surf_tot-1
  do j=1, n_seg-1

    index = index + 1
    newelement_list%element(index)%size = 1.d0

    if (which_leg .eq. 1) then
      newelement_list%element(index)%vertex(1) = n_tmp + (i-1)*n_seg + j
      newelement_list%element(index)%vertex(2) = n_tmp + (i  )*n_seg + j
      newelement_list%element(index)%vertex(3) = n_tmp + (i  )*n_seg + j + 1
      newelement_list%element(index)%vertex(4) = n_tmp + (i-1)*n_seg + j + 1
      ! Connect with Xpoint properly
      if ( (i .eq. i_sep-1) .and. (j .eq. 1) ) then
        newelement_list%element(index)%vertex(2) = n_xpoint_1
      endif  
      if ( (i .eq. i_sep) .and. (j .eq. 1) ) then
        newelement_list%element(index)%vertex(1) = n_xpoint_2
      endif  
    endif  
    if (which_leg .eq. 2) then
      newelement_list%element(index)%vertex(1) = n_tmp + (i-1)*n_seg + j + 1
      newelement_list%element(index)%vertex(2) = n_tmp + (i  )*n_seg + j + 1
      newelement_list%element(index)%vertex(3) = n_tmp + (i  )*n_seg + j
      newelement_list%element(index)%vertex(4) = n_tmp + (i-1)*n_seg + j
      ! Connect with Xpoint properly
      if ( (i .eq. i_sep-1) .and. (j .eq. 1) ) then
        newelement_list%element(index)%vertex(3) = n_xpoint_1
      endif  
      if ( (i .eq. i_sep) .and. (j .eq. 1) ) then
        newelement_list%element(index)%vertex(4) = n_xpoint_2
      endif  
    endif
    if (which_leg .eq. 3) then
      newelement_list%element(index)%vertex(1) = n_tmp + (i-1)*n_seg + j + 1
      newelement_list%element(index)%vertex(2) = n_tmp + (i  )*n_seg + j + 1
      newelement_list%element(index)%vertex(3) = n_tmp + (i  )*n_seg + j
      newelement_list%element(index)%vertex(4) = n_tmp + (i-1)*n_seg + j
      ! Connect with Xpoint properly
      if ( (i .eq. i_sep-1) .and. (j .eq. 1) ) then
        newelement_list%element(index)%vertex(3) = n_xpoint_1
      endif  
      if ( (i .eq. i_sep) .and. (j .eq. 1) ) then
        newelement_list%element(index)%vertex(4) = n_xpoint_2
      endif  
    endif  
    if (which_leg .eq. 4) then
      newelement_list%element(index)%vertex(1) = n_tmp + (i-1)*n_seg + j
      newelement_list%element(index)%vertex(2) = n_tmp + (i  )*n_seg + j
      newelement_list%element(index)%vertex(3) = n_tmp + (i  )*n_seg + j + 1
      newelement_list%element(index)%vertex(4) = n_tmp + (i-1)*n_seg + j + 1
      ! Connect with Xpoint properly
      if ( (i .eq. i_sep-1) .and. (j .eq. 1) ) then
        newelement_list%element(index)%vertex(2) = n_xpoint_1
      endif  
      if ( (i .eq. i_sep) .and. (j .eq. 1) ) then
        newelement_list%element(index)%vertex(1) = n_xpoint_2
      endif  
    endif
      
  enddo
enddo
newelement_list%n_elements = index


!----------------------------------- Print a python file that plots a cross with the 4 nodes of each element
if (plot_grid) then
  write(char_tmp,'(i1)')which_leg
  plot_filename = 'plot_leg_elements'//char_tmp//'.py'
  n_loop = newelement_list%n_elements
  open(101,file=plot_filename)
    write(101,'(A)')         '#!/usr/bin/env python'
    write(101,'(A)')         'import numpy as N'
    write(101,'(A)')         'import pylab'
    write(101,'(A)')         'def main():'
    write(101,'(A,i6,A)')    ' r = N.zeros(',4*n_loop,')'
    write(101,'(A,i6,A)')    ' z = N.zeros(',4*n_loop,')'
    do j=1,n_loop
      do i=1,2
        index = newelement_list%element(j)%vertex(i)
        write(101,'(A,i6,A,f15.4)') ' r[',4*(j-1)+2*i-2,'] = ',newnode_list%node(index)%x(1,1,1)
        write(101,'(A,i6,A,f15.4)') ' z[',4*(j-1)+2*i-2,'] = ',newnode_list%node(index)%x(1,1,2)
        index = newelement_list%element(j)%vertex(i+2)
        write(101,'(A,i6,A,f15.4)') ' r[',4*(j-1)+2*i-1,'] = ',newnode_list%node(index)%x(1,1,1)
        write(101,'(A,i6,A,f15.4)') ' z[',4*(j-1)+2*i-1,'] = ',newnode_list%node(index)%x(1,1,2)
      enddo
    enddo
    write(101,'(A,i6,A)')    ' for i in range (0,',n_loop*2,'):'
    write(101,'(A)')         '  pylab.plot(r[2*i:2*i+2],z[2*i:2*i+2], "r")'
    do j=1,n_loop
      do i=1,4
        index = newelement_list%element(j)%vertex(i)
        write(101,'(A,i6,A,f15.4)') ' r[',4*(j-1)+i-1,'] = ',newnode_list%node(index)%x(1,1,1)
        write(101,'(A,i6,A,f15.4)') ' z[',4*(j-1)+i-1,'] = ',newnode_list%node(index)%x(1,1,2)
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




deallocate(seg, R_seg, Z_seg, R_seg_surf, Z_seg_surf, R_seg_smooth, Z_seg_smooth)
deallocate(R_polar,Z_polar)
deallocate(i_flux)

return
end subroutine define_leg_grid




























subroutine segment_surface_length(node_list,element_list,surface, R_beg, Z_beg, R_end, Z_end, n_seg, seg, R_seg, Z_seg, xpoint_surface)

  use data_structure
  use phys_module, only: xcase
  use mod_interp, only: interp_RZ
  
  implicit none
  
  ! --- Routine variables
  type (type_node_list),    intent(in)   :: node_list
  type (type_element_list), intent(in)   :: element_list
  type (type_surface),      intent(in)   :: surface
  real*8,                   intent(in)   :: R_beg, Z_beg, R_end, Z_end
  integer,                  intent(in)   :: n_seg
  real*8,                   intent(in)   :: seg(n_seg) ! segment from 0 to 1 (usually from meshac)
  real*8,                   intent(out)  :: R_seg(n_seg), Z_seg(n_seg)
  logical,                  intent(in)   :: xpoint_surface
  
  ! --- Internal variables
  type (type_surface_list) :: surface_list_tmp
  integer :: i_part, i_piece, i_elm, i_find, j_find, i_seg, diff_pieces_min, diff_pieces, i_RZ
  integer :: i_part_beg, i_piece_beg
  integer :: i_part_end, i_piece_end
  integer :: i_elm_find(10)
  integer :: i_start, i_stop, increment
  real*8  :: s_find(10),t_find(10), st_find(10)
  real*8  :: diff_min_beg, diff_beg
  real*8  :: diff_min_end, diff_end
  real*8  :: st_beg, st_end, st_prev
  real*8  :: st_start, st_stop
  real*8  :: st1, st2
  real*8  :: rr, ss, st, s_tmp, t_tmp
  real*8  :: rr1, ss1, drr1, dss1
  real*8  :: rr2, ss2, drr2, dss2
  real*8  :: dR1_dt, dZ1_dt, dR2_dt, dZ2_dt
  real*8  :: R1,dR1_dr,dR1_ds,dR1_drs,dR1_drr,dR1_dss
  real*8  :: Z1,dZ1_dr,dZ1_ds,dZ1_drs,dZ1_drr,dZ1_dss
  real*8  :: R2,dR2_dr,dR2_ds,dR2_drs,dR2_drr,dR2_dss
  real*8  :: Z2,dZ2_dr,dZ2_ds,dZ2_drs,dZ2_drr,dZ2_dss
  real*8  :: R3,dR3_dr,dR3_ds,dR3_drs,dR3_drr,dR3_dss
  real*8  :: Z3,dZ3_dr,dZ3_ds,dZ3_drs,dZ3_drr,dZ3_dss
  real*8  :: surface_length, length, length_sum, length_seg
  real*8, parameter :: tol_find = 5.d-4
  
  ! --- Find the corresponding end points on each surface
  allocate(surface_list_tmp%psi_values(1))
  allocate(surface_list_tmp%flux_surfaces(1))
  surface_list_tmp%n_psi = 1
  surface_list_tmp%flux_surfaces(1)%n_pieces = 1
  i_part_end = 0
  i_part_beg = 0
  ! --- First the end point
  diff_min_end = 1.d10
  do i_part = 1,surface%n_parts
  
    do i_piece = surface%parts_index(i_part), surface%parts_index(i_part+1)-1
    
      ! --- End points
      rr    = surface%s(1,i_piece)
      ss    = surface%t(1,i_piece)
      i_elm = surface%elm(i_piece)
      call interp_RZ(node_list,element_list,i_elm,rr,ss,R1,dR1_dr,dR1_ds,dR1_drs,dR1_drr,dR1_dss, &
                                                        Z1,dZ1_dr,dZ1_ds,dZ1_drs,dZ1_drr,dZ1_dss)
      rr    = surface%s(3,i_piece)
      ss    = surface%t(3,i_piece)
      i_elm = surface%elm(i_piece)
      call interp_RZ(node_list,element_list,i_elm,rr,ss,R2,dR2_dr,dR2_ds,dR2_drs,dR2_drr,dR2_dss, &
                                                        Z2,dZ2_dr,dZ2_ds,dZ2_drs,dZ2_drr,dZ2_dss)
      
      surface_list_tmp%flux_surfaces(1)%s(:,1) = surface%s(:,i_piece)
      surface_list_tmp%flux_surfaces(1)%t(:,1) = surface%t(:,i_piece)
      surface_list_tmp%flux_surfaces(1)%elm(1) = surface%elm(i_piece)
      
      ! --- End point
      do i_RZ=1,2
        i_find = 0
        if (i_RZ .eq. 1) then
          if ( (min(Z1,Z2) .le. Z_end) .and. (Z_end .le. max(Z1,Z2)) ) then
            call find_Z_surface(node_list,element_list,surface_list_tmp,1,Z_end,i_elm_find,s_find,t_find,st_find,i_find)
          endif
        else
          if ( (min(R1,R2) .le. R_end) .and. (R_end .le. max(R1,R2)) ) then
            call find_R_surface(node_list,element_list,surface_list_tmp,1,R_end,i_elm_find,s_find,t_find,st_find,i_find)
          endif
        endif
        do j_find = 1,i_find
          rr    = s_find(j_find)
          ss    = t_find(j_find)
          i_elm = i_elm_find(j_find)
          call interp_RZ(node_list,element_list,i_elm,rr,ss,R3,dR3_dr,dR3_ds,dR3_drs,dR3_drr,dR3_dss, &
                                                            Z3,dZ3_dr,dZ3_ds,dZ3_drs,dZ3_drr,dZ3_dss)
          diff_end = sqrt( (R3-R_end)**2 + (Z3-Z_end)**2 )
          if (diff_end .lt. diff_min_end) then
            diff_min_end = diff_end
            i_part_end  = i_part
            i_piece_end = i_piece
            st_end      = st_find(j_find)
          endif
        enddo
      enddo
    
    enddo
  enddo
  ! --- Then beg point
  diff_min_beg = 1.d10
  diff_pieces_min = 10000000
  do i_part = 1,surface%n_parts
  
    do i_piece = surface%parts_index(i_part), surface%parts_index(i_part+1)-1
    
      ! --- End points
      rr    = surface%s(1,i_piece)
      ss    = surface%t(1,i_piece)
      i_elm = surface%elm(i_piece)
      call interp_RZ(node_list,element_list,i_elm,rr,ss,R1,dR1_dr,dR1_ds,dR1_drs,dR1_drr,dR1_dss, &
                                                        Z1,dZ1_dr,dZ1_ds,dZ1_drs,dZ1_drr,dZ1_dss)
      rr    = surface%s(3,i_piece)
      ss    = surface%t(3,i_piece)
      i_elm = surface%elm(i_piece)
      call interp_RZ(node_list,element_list,i_elm,rr,ss,R2,dR2_dr,dR2_ds,dR2_drs,dR2_drr,dR2_dss, &
                                                        Z2,dZ2_dr,dZ2_ds,dZ2_drs,dZ2_drr,dZ2_dss)
      
      surface_list_tmp%flux_surfaces(1)%s(:,1) = surface%s(:,i_piece)
      surface_list_tmp%flux_surfaces(1)%t(:,1) = surface%t(:,i_piece)
      surface_list_tmp%flux_surfaces(1)%elm(1) = surface%elm(i_piece)
      
      ! --- Beg point
      i_find = 0
      do i_RZ=1,2
        if (i_RZ .eq. 1) then
          if ( (min(Z1,Z2) .le. Z_beg) .and. (Z_beg .le. max(Z1,Z2)) ) then
            call find_Z_surface(node_list,element_list,surface_list_tmp,1,Z_beg,i_elm_find,s_find,t_find,st_find,i_find)
          endif
        else
          if ( (min(R1,R2) .le. R_beg) .and. (R_beg .le. max(R1,R2)) ) then
            call find_R_surface(node_list,element_list,surface_list_tmp,1,R_beg,i_elm_find,s_find,t_find,st_find,i_find)
          endif
        endif
        do j_find = 1,i_find
          rr    = s_find(j_find)
          ss    = t_find(j_find)
          i_elm = i_elm_find(j_find)
          call interp_RZ(node_list,element_list,i_elm,rr,ss,R3,dR3_dr,dR3_ds,dR3_drs,dR3_drr,dR3_dss, &
                                                            Z3,dZ3_dr,dZ3_ds,dZ3_drs,dZ3_drr,dZ3_dss)
          diff_beg = sqrt( (R3-R_beg)**2 + (Z3-Z_beg)**2 )
          if ( (.not. xpoint_surface) .and. (diff_beg .le. diff_min_beg) ) then
            diff_min_beg = diff_beg
            i_part_beg  = i_part
            i_piece_beg = i_piece
            st_beg      = st_find(j_find)
          endif
          if ( (xpoint_surface) .and. (diff_beg .le. tol_find) ) then
            if (i_part .eq. i_part_end) then
              diff_pieces = abs(i_piece-i_piece_end)
              if (diff_pieces .lt. diff_pieces_min) then
                diff_pieces_min = diff_pieces
                diff_min_beg = diff_beg
                i_part_beg  = i_part
                i_piece_beg = i_piece
                st_beg      = st_find(j_find)
              endif
            endif
          endif
        enddo
      enddo
    
    enddo
  enddo
  deallocate(surface_list_tmp%psi_values)
  deallocate(surface_list_tmp%flux_surfaces)
  
  ! --- both end points should be on the same part
  if (i_part_end .ne. i_part_beg) then
    write(*,*)'Problem segmenting surface: end points on different surface parts!',i_part_end,i_part_beg
    stop
  endif
  
  ! --- Get total segment length
  i_start   = i_piece_beg
  i_stop    = i_piece_end
  st_start  = st_beg
  st_stop   = st_end
  if (i_piece_end .gt. i_piece_beg) then
    increment = +1
  else
    increment = -1
  endif
  
  surface_length = 0.d0
  do i_piece = i_start, i_stop, increment
    if (i_piece .eq. i_start) then
      st1 = st_start
      st2 = real(increment)
    elseif (i_piece .eq. i_stop) then
      st1 = -real(increment)
      st2 = st_end
    else
      st1 = -real(increment)
      st2 = real(increment)
    endif
    
    rr1  = surface%s(1,i_piece);   ss1  = surface%t(1,i_piece)
    drr1 = surface%s(2,i_piece);   dss1 = surface%t(2,i_piece)
    rr2  = surface%s(3,i_piece);   ss2  = surface%t(3,i_piece)
    drr2 = surface%s(4,i_piece);   dss2 = surface%t(4,i_piece)
    i_elm = surface%elm(i_piece)
    call interp_RZ(node_list,element_list,i_elm,rr1,ss1,R1,dR1_dr,dR1_ds,dR1_drs,dR1_drr,dR1_dss, &
                                                        Z1,dZ1_dr,dZ1_ds,dZ1_drs,dZ1_drr,dZ1_dss)
    call interp_RZ(node_list,element_list,i_elm,rr2,ss2,R2,dR2_dr,dR2_ds,dR2_drs,dR2_drr,dR2_dss, &
                                                        Z2,dZ2_dr,dZ2_ds,dZ2_drs,dZ2_drr,dZ2_dss)
    dR1_dt = dR1_dr * drr1 + dR1_ds * dss1
    dZ1_dt = dZ1_dr * drr1 + dZ1_ds * dss1
    dR2_dt = dR2_dr * drr2 + dR2_ds * dss2
    dZ2_dt = dZ2_dr * drr2 + dZ2_ds * dss2
  
    call curve_length(R1, dR1_dt, R2, dR2_dt, Z1, dZ1_dt, Z2, dZ2_dt, st1, st2, length)
    surface_length = surface_length + length
    
  enddo
  !write(*,*)'surface length:',surface_length
  
  ! --- Now that we know the surface length, segment it
  do i_seg = 1,n_seg
    
    length_seg = seg(i_seg) * surface_length
    
    length_sum = 0.d0
    do i_piece = i_start, i_stop, increment
      if (i_piece .eq. i_start) then
        st1 = st_start
        st2 = real(increment)
      elseif (i_piece .eq. i_stop) then
        st1 = -real(increment)
        st2 = st_end
      else
        st1 = -real(increment)
        st2 = real(increment)
      endif
      
      rr1  = surface%s(1,i_piece);   ss1  = surface%t(1,i_piece)
      drr1 = surface%s(2,i_piece);   dss1 = surface%t(2,i_piece)
      rr2  = surface%s(3,i_piece);   ss2  = surface%t(3,i_piece)
      drr2 = surface%s(4,i_piece);   dss2 = surface%t(4,i_piece)
      i_elm = surface%elm(i_piece)
      call interp_RZ(node_list,element_list,i_elm,rr1,ss1,R1,dR1_dr,dR1_ds,dR1_drs,dR1_drr,dR1_dss, &
                                                          Z1,dZ1_dr,dZ1_ds,dZ1_drs,dZ1_drr,dZ1_dss)
      call interp_RZ(node_list,element_list,i_elm,rr2,ss2,R2,dR2_dr,dR2_ds,dR2_drs,dR2_drr,dR2_dss, &
                                                          Z2,dZ2_dr,dZ2_ds,dZ2_drs,dZ2_drr,dZ2_dss)
      dR1_dt = dR1_dr * drr1 + dR1_ds * dss1
      dZ1_dt = dZ1_dr * drr1 + dZ1_ds * dss1
      dR2_dt = dR2_dr * drr2 + dR2_ds * dss2
      dZ2_dt = dZ2_dr * drr2 + dZ2_ds * dss2
    
      call curve_length(R1, dR1_dt, R2, dR2_dt, Z1, dZ1_dt, Z2, dZ2_dt, st1, st2, length)
      
      ! --- We found the right piece
      if (length_sum + length .ge. length_seg) then
        length_seg = length_seg - length_sum
        st = length_seg/length
        st = (st2-st1)*st + st1
        call CUB1D(R1, dR1_dt, R2, dR2_dt, st, R3,dR3_dr)
        call CUB1D(Z1, dZ1_dt, Z2, dZ2_dt, st, Z3,dZ3_dr)
        R_seg(i_seg) = R3
        Z_seg(i_seg) = Z3
        st_prev = st
        exit
      else
        length_sum = length_sum + length
      endif
      
    enddo
    
  enddo
  
  
  
  
  return

end subroutine segment_surface_length
















subroutine curve_length(R1, dR1, R2, dR2, Z1, dZ1, Z2, dZ2, s_beg, s_end, length)

  use gauss

  implicit none
  
  ! --- Routine parameters
  real*8, intent(in)  :: R1, dR1, R2, dR2
  real*8, intent(in)  :: Z1, dZ1, Z2, dZ2
  real*8, intent(in)  :: s_beg, s_end
  real*8, intent(out) :: length
  
  ! --- Internal parameters
  integer :: i_gauss
  real*8  :: s_norm, ss, ds, integrand, ws
  real*8  :: R_tmp,  Z_tmp
  real*8  :: dR_ds,  dZ_ds
  real*8  :: dR_dsn, dZ_dsn
  
  
  ! --- Integrate curve on gaussian points
  length = 0.d0
  do i_gauss = 1,n_gauss
    s_norm = xgauss(i_gauss)
    ws = wgauss(i_gauss)
    ss = (s_end-s_beg)*s_norm + s_beg
    ds = (s_end-s_beg)
    call CUB1D(R1, dR1, R2, dR2, ss, R_tmp, dR_ds)
    call CUB1D(Z1, dZ1, Z2, dZ2, ss, Z_tmp, dZ_ds)
    dR_dsn = dR_ds * ds
    dZ_dsn = dZ_ds * ds
    integrand = sqrt( dR_dsn**2 + dZ_dsn**2 )
    length = length + integrand * ws
  enddo

  return

end subroutine curve_length

