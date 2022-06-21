!> Subroutine defines the new grid_points from crossing of polar and radial coordinate lines
subroutine define_central_grid(node_list, element_list, flux_list, &
                               xcase, n_grids, stpts, sigmas, nwpts, newnode_list, newelement_list)

use constants
use tr_module 
use data_structure
use grid_xpoint_data
use phys_module, only:   tokamak_device, n_tht_equidistant, SDN_threshold
use mod_interp, only: interp_RZ
use equil_info

implicit none

! --- Routine parameters
type (type_surface_list)    , intent(inout) :: flux_list
type (type_node_list)       , intent(inout) :: node_list
type (type_element_list)    , intent(inout) :: element_list
type (type_strategic_points), intent(inout) :: stpts
type (type_new_points)      , intent(inout) :: nwpts
integer,                      intent(inout) :: n_grids(12)
integer,                      intent(in)    :: xcase
real*8,                       intent(in)    :: sigmas(17)
type (type_node_list)       , intent(inout) :: newnode_list
type (type_element_list)    , intent(inout) :: newelement_list

! --- local variables
real*8, allocatable :: s_tmp(:), s_tmp2(:), theta_sep(:), theta_beg(:)
real*8, allocatable :: delta1(:),delta2(:)
real*8, allocatable :: R_polar(:,:,:),Z_polar(:,:,:)
real*8, allocatable :: R_beg(:),Z_beg(:)
integer             :: i, j, k, l, m, j2
integer             :: index
integer             :: n_flux_mid
integer             :: i_sep1, i_sep2, i_max, n_start, i_surf
integer             :: n_psi, n_tht_2, n_tht_mid, n_tht_mid2
integer             :: n_flux, n_tht,  n_open,   n_outer,   n_inner
integer             :: n_private,   n_up_priv,   n_leg,   n_up_leg
integer             :: ifail
integer             :: n_xpoint_1, n_xpoint_2, n_xpoint_3
integer             :: n_loop, n_loop2, n_tmp
integer             :: n_start_open, n_start_outer, n_start_inner
integer             :: n_start_private, n_start_up_priv
real*8              :: R_cub1d(4), Z_cub1d(4)
real*8              :: tht_x1, tht_x2
real*8              :: psi_norm, psi_bnd
real*8              :: RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss
real*8              :: ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss
integer             :: i_elm_find(8), i_find
real*8              :: s_find(8), t_find(8), st_find(8)
real*8              :: diff_min
real*8              :: SIG_theta, SIG_theta_up
real*8              :: SIG_leg_0, SIG_leg_1
real*8              :: SIG_up_leg_0, SIG_up_leg_1
real*8              :: SIG_0, SIG_1
real*8              :: bgf_tht
real*8              :: Zbeg, Zend, tht_SOL
real*8              :: scale_out_points
logical, parameter  :: plot_grid = .true.

write(*,*) '*****************************************'
write(*,*) '* X-point grid inside wall :            *'
write(*,*) '*****************************************'
write(*,*) '                 Define central part of grid'

if(xcase .eq. LOWER_XPOINT) psi_bnd = ES%psi_xpoint(1)
if(xcase .eq. UPPER_XPOINT) psi_bnd = ES%psi_xpoint(2)
if(xcase .eq. DOUBLE_NULL ) then
  if (ES%active_xpoint .eq. UPPER_XPOINT) then
    psi_bnd  = ES%psi_xpoint(2)
  else
    psi_bnd  = ES%psi_xpoint(1)
  endif
  ! If we have a symmetric double-null, force the single separatrix
  if (ES%active_xpoint .eq. SYMMETRIC_XPOINT) then
    psi_bnd  = (ES%psi_xpoint(1)+ES%psi_xpoint(2))/2.d0
  endif
endif

SIG_theta    = sigmas(2) ; SIG_theta_up = sigmas(17)
SIG_leg_0    = sigmas(8) ; SIG_leg_1    = sigmas(9) 
SIG_up_leg_0 = sigmas(10); SIG_up_leg_1 = sigmas(11)

n_flux    = n_grids(1); n_tht     = n_grids(2)
n_open    = n_grids(3); n_outer   = n_grids(4); n_inner = n_grids(5)
n_private = n_grids(6); n_up_priv = n_grids(7)
n_leg     = n_grids(8); n_up_leg  = n_grids(9)

nwpts%k_cross = 0

! this is the background of the distribution for the theta grid-points
bgf_tht  = 0.8d0!0.97d0


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



!-------------------------------------------------------------------------------------------!
!------- Extrapolation points for first part of the grid (everything except legs) ----------!
!-------------------------------------------------------------------------------------------!

!-------------------------------- define the polar coordinate
n_psi   = n_flux + n_open + n_outer + n_inner + n_private + n_up_priv + 1   ! this includes the magnetic axis
n_tht_2 = n_tht + 2*n_leg + 2*n_up_leg

! --- Find middle surface
diff_min = 1.d10
n_flux_mid = 0
do i=1,n_flux
  psi_norm = (flux_list%psi_values(i)-ES%psi_axis)/(psi_bnd-ES%psi_axis)
  if (abs(psi_norm-0.5) .lt. diff_min) then
    diff_min = abs(psi_norm-0.5)
    n_flux_mid = i
  endif
enddo

allocate(theta_sep(n_tht))
allocate(theta_beg(n_tht))
allocate(s_tmp (n_tht))
allocate(s_tmp2(n_tht))
  
if (xcase .ne. DOUBLE_NULL) then 
  if (xcase .eq. LOWER_XPOINT) tht_x1 = atan2(ES%Z_xpoint(1)-ES%Z_axis,ES%R_xpoint(1)-ES%R_axis)
  if (xcase .eq. UPPER_XPOINT) tht_x1 = atan2(ES%Z_xpoint(2)-ES%Z_axis,ES%R_xpoint(2)-ES%R_axis)
  n_tht_mid = n_tht/2
  
  ! s_tmp is for the separatrix, where we focus grid-points near the X-points
  ! s_tmp2 is for the grid centre, which we want to be equidistant (ie. not focused at the X-points)
  call meshac2(n_tht,s_tmp, 0.d0,1.d0,SIG_theta,SIG_theta,bgf_tht,1.0d0)
  do j=1,n_tht
    theta_sep(j) = tht_x1 + 2.d0 * PI * s_tmp(j)
    theta_beg(j) = tht_x1 + 2.d0 * PI * real(j-1)/real(n_tht-1)
  enddo
else ! xcase == DOUBLE_NULL
  if ( (ES%active_xpoint .eq. LOWER_XPOINT) .or. (ES%active_xpoint .eq. SYMMETRIC_XPOINT) ) then
    tht_x1 = atan2(ES%Z_xpoint(1)-ES%Z_axis,ES%R_xpoint(1)-ES%R_axis)
    tht_x2 = atan2(ES%Z_xpoint(2)-ES%Z_axis,ES%R_xpoint(2)-ES%R_axis)
  else
    tht_x2 = atan2(ES%Z_xpoint(1)-ES%Z_axis,ES%R_xpoint(1)-ES%R_axis)
    tht_x1 = atan2(ES%Z_xpoint(2)-ES%Z_axis,ES%R_xpoint(2)-ES%R_axis)
  endif
  if (tht_x1 .lt. 0.d0) tht_x1 = tht_x1 + 2.d0 * PI
  if (tht_x2 .lt. 0.d0) tht_x2 = tht_x2 + 2.d0 * PI
  
  ! Spread out points evenly (outer angle between tht_x1 and tht_x2 is usually bigger than inner angle)
  if ( (ES%active_xpoint .eq. LOWER_XPOINT ) .or. (ES%active_xpoint .eq. SYMMETRIC_XPOINT) ) then
    ! n_tht_mid corresponds to the secondary X-point
    n_tht_mid = int(n_tht * (2.d0*PI - (tht_x1 - tht_x2)) / (2.d0*PI))
    ! Make sure n_tht_mid is odd and save it to n_grids for later use
    ! it needs to be odd, because we go through each Xpoint twice:
    ! from 1 to n_tht_mid, and then from n_tht_mid+1 to n_tht
    if(mod(n_tht_mid,2) .eq. 0) n_tht_mid = n_tht_mid + 1
    
    call meshac2(n_tht_mid,s_tmp, 0.d0,1.d0,SIG_theta,SIG_theta_up,bgf_tht,1.0d0)
    call meshac2(n_tht_mid,s_tmp2,0.d0,1.d0,999.0,    999.0,       bgf_tht,1.0d0)
    do j=1,n_tht_mid
      theta_sep(j) = (tht_x1-2.d0*PI) + (tht_x2-(tht_x1-2.d0*PI)) * s_tmp(j)
      theta_beg(j) = (tht_x1-2.d0*PI) + (tht_x2-(tht_x1-2.d0*PI)) * s_tmp2(j)
    enddo
  
    n_tht_mid2 = n_tht-n_tht_mid
    call meshac2(n_tht_mid2,s_tmp, 0.d0,1.d0,SIG_theta_up,SIG_theta,bgf_tht,1.0d0)
    call meshac2(n_tht_mid2,s_tmp2,0.d0,1.d0,999.0,    999.0,       bgf_tht,1.0d0)
    do j=n_tht_mid+1,n_tht
      theta_sep(j) = tht_x2 + (tht_x1-tht_x2) * s_tmp(j-n_tht_mid)
      theta_beg(j) = tht_x2 + (tht_x1-tht_x2) * s_tmp2(j-n_tht_mid)
    enddo
  else
    n_tht_mid = int(n_tht * (tht_x2 - tht_x1) / (2.d0*PI))
    ! Make sure n_tht_mid is odd and save it to n_grids for later use
    if(mod(n_tht_mid,2) .eq. 0) n_tht_mid = n_tht_mid + 1
    
    call meshac2(n_tht_mid,s_tmp, 0.d0,1.d0,SIG_theta_up,SIG_theta,bgf_tht,1.0d0)
    call meshac2(n_tht_mid,s_tmp2,0.d0,1.d0,999.0,       999.0,    bgf_tht,1.0d0)
    do j=1,n_tht_mid
      theta_sep(j) = tht_x1 + (tht_x2-tht_x1) * s_tmp(j-n_tht_mid)
      theta_beg(j) = tht_x1 + (tht_x2-tht_x1) * s_tmp2(j-n_tht_mid)
    enddo

    n_tht_mid2 = n_tht-n_tht_mid
    call meshac2(n_tht_mid2,s_tmp, 0.d0,1.d0,SIG_theta,SIG_theta_up,bgf_tht,1.0d0)
    call meshac2(n_tht_mid2,s_tmp2,0.d0,1.d0,999.0,    999.0,       bgf_tht,1.0d0)
    do j=n_tht_mid+1,n_tht
      theta_sep(j) = (tht_x2-2.d0*PI) + (tht_x1-(tht_x2-2.d0*PI)) * s_tmp(j-n_tht_mid)
      theta_beg(j) = (tht_x2-2.d0*PI) + (tht_x1-(tht_x2-2.d0*PI)) * s_tmp2(j-n_tht_mid)
    enddo
  
  endif      
endif
if (.not. n_tht_equidistant) theta_beg = theta_sep

call tr_deallocate(s_tmp, "s_tmp", CAT_GRID)
call tr_deallocate(s_tmp2,"s_tmp2",CAT_GRID)
do j=1,n_tht
  if (theta_sep(j) .lt. 0.d0)    theta_sep(j) = theta_sep(j) + 2.d0 * PI
  if (theta_sep(j) .gt. 2.d0*PI) theta_sep(j) = theta_sep(j) - 2.d0 * PI
  if (theta_beg(j) .lt. 0.d0)    theta_beg(j) = theta_beg(j) + 2.d0 * PI
  if (theta_beg(j) .gt. 2.d0*PI) theta_beg(j) = theta_beg(j) - 2.d0 * PI
enddo


!------------------------------------- find crossing with middle surface
i_max = n_flux_mid
do j=1,n_tht
  call find_theta_surface(node_list,element_list,flux_list,i_max,theta_beg(j),ES%R_axis,ES%Z_axis,i_elm_find,s_find,t_find,i_find)
  if(i_find .eq. 0) return
  
  do k=1,i_find
    call interp_RZ(node_list,element_list,i_elm_find(k),s_find(k),t_find(k),&
                   RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss,    &
                   ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)
    if(     ( (xcase .eq. LOWER_XPOINT) .and. (ZZg1 .gt. ES%Z_xpoint(1)) ) &
       .or. ( (xcase .eq. UPPER_XPOINT) .and. (ZZg1 .lt. ES%Z_xpoint(2)) ) &
       .or. ( (xcase .eq. DOUBLE_NULL ) .and. (ZZg1 .gt. ES%Z_xpoint(1)) .and. (ZZg1 .lt. ES%Z_xpoint(2)) ) ) then
      nwpts%R_mid(j) = RRg1
      nwpts%Z_mid(j) = ZZg1
      exit
    endif
  enddo
  
enddo


!------------------------------------- find crossing with separatrix
i_sep1  = n_flux
i_sep2  = n_flux + n_open
if (xcase .ne. DOUBLE_NULL) then
  do j=2,n_tht-1
    call find_theta_surface(node_list,element_list,flux_list,i_sep1,theta_sep(j),ES%R_axis,ES%Z_axis,i_elm_find,s_find,t_find,i_find)
    if(i_find .eq. 0) return
    
    do k=1,i_find
      call interp_RZ(node_list,element_list,i_elm_find(k),s_find(k),t_find(k),&
                     RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss,    &
                     ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)
      if(     ( (xcase .eq. LOWER_XPOINT) .and. (ZZg1 .ge. ES%Z_xpoint(1)) ) &
         .or. ( (xcase .eq. UPPER_XPOINT) .and. (ZZg1 .le. ES%Z_xpoint(2)) ) ) then
        nwpts%R_sep(j) = RRg1
        nwpts%Z_sep(j) = ZZg1
        exit
      endif
    enddo
    
  enddo
else ! xcase == DOUBLE_NULL
  do j=2,n_tht-1
    if((j .ne. n_tht_mid) .and. (j .ne. n_tht_mid+1)) then
      if (theta_sep(j) .ge. pi) then
        if ( (ES%active_xpoint .eq. LOWER_XPOINT) .or. (ES%active_xpoint .eq. SYMMETRIC_XPOINT) ) then
          call find_theta_surface(node_list,element_list,flux_list,i_sep1,theta_sep(j),ES%R_axis,ES%Z_axis,i_elm_find,s_find,t_find,i_find)
        else
          call find_theta_surface(node_list,element_list,flux_list,i_sep2,theta_sep(j),ES%R_axis,ES%Z_axis,i_elm_find,s_find,t_find,i_find)
        endif
        if(i_find .eq. 0) return

        do k=1,i_find
          call interp_RZ(node_list,element_list,i_elm_find(k),s_find(k),t_find(k),&
                         RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss,    &
                         ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)
          if(ZZg1 .ge. ES%Z_xpoint(1)) then
            nwpts%R_sep(j) = RRg1
            nwpts%Z_sep(j) = ZZg1
          endif
        enddo
      else
        if ( (ES%active_xpoint .eq. LOWER_XPOINT) .or. (ES%active_xpoint .eq. SYMMETRIC_XPOINT) ) then
          call find_theta_surface(node_list,element_list,flux_list,i_sep2,theta_sep(j),ES%R_axis,ES%Z_axis,i_elm_find,s_find,t_find,i_find)
        else
          call find_theta_surface(node_list,element_list,flux_list,i_sep1,theta_sep(j),ES%R_axis,ES%Z_axis,i_elm_find,s_find,t_find,i_find)
        endif
        if(i_find .eq. 0) return
        
        do k=1,i_find
          call interp_RZ(node_list,element_list,i_elm_find(k),s_find(k),t_find(k),&
                         RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss,    &
                         ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)
          if(ZZg1 .le. ES%Z_xpoint(2)) then
            nwpts%R_sep(j) = RRg1
            nwpts%Z_sep(j) = ZZg1    
          endif
        enddo
      endif
    endif
  enddo
endif

if (xcase .eq. LOWER_XPOINT) then
  nwpts%R_sep(1)             = ES%R_xpoint(1) ! this one is known - safer...
  nwpts%Z_sep(1)             = ES%Z_xpoint(1) ! this one is known - safer...
  nwpts%R_sep(n_tht)         = ES%R_xpoint(1) ! this one is known - safer...
  nwpts%Z_sep(n_tht)         = ES%Z_xpoint(1) ! this one is known - safer...
endif
if (xcase .eq. UPPER_XPOINT) then
  nwpts%R_sep(1)             = ES%R_xpoint(2) ! this one is known - safer...
  nwpts%Z_sep(1)             = ES%Z_xpoint(2) ! this one is known - safer...
  nwpts%R_sep(n_tht)         = ES%R_xpoint(2) ! this one is known - safer...
  nwpts%Z_sep(n_tht)         = ES%Z_xpoint(2) ! this one is known - safer...
endif
if (xcase .eq. DOUBLE_NULL) then
  if ( (ES%active_xpoint .eq. LOWER_XPOINT) .or. (ES%active_xpoint .eq. SYMMETRIC_XPOINT) ) then
    nwpts%R_sep(1)           = ES%R_xpoint(1) ! this one is known - safer...
    nwpts%Z_sep(1)           = ES%Z_xpoint(1) ! this one is known - safer...
    nwpts%R_sep(n_tht)       = ES%R_xpoint(1) ! this one is known - safer...
    nwpts%Z_sep(n_tht)       = ES%Z_xpoint(1) ! this one is known - safer...
    nwpts%R_sep(n_tht_mid)   = ES%R_xpoint(2) ! this one is known - safer...
    nwpts%Z_sep(n_tht_mid)   = ES%Z_xpoint(2) ! this one is known - safer...
    nwpts%R_sep(n_tht_mid+1) = ES%R_xpoint(2) ! this one is known - safer...
    nwpts%Z_sep(n_tht_mid+1) = ES%Z_xpoint(2) ! this one is known - safer...
  else
    nwpts%R_sep(1)           = ES%R_xpoint(2) ! this one is known - safer...
    nwpts%Z_sep(1)           = ES%Z_xpoint(2) ! this one is known - safer...
    nwpts%R_sep(n_tht)       = ES%R_xpoint(2) ! this one is known - safer...
    nwpts%Z_sep(n_tht)       = ES%Z_xpoint(2) ! this one is known - safer...
    nwpts%R_sep(n_tht_mid)   = ES%R_xpoint(1) ! this one is known - safer...
    nwpts%Z_sep(n_tht_mid)   = ES%Z_xpoint(1) ! this one is known - safer...
    nwpts%R_sep(n_tht_mid+1) = ES%R_xpoint(1) ! this one is known - safer...
    nwpts%Z_sep(n_tht_mid+1) = ES%Z_xpoint(1) ! this one is known - safer...
  endif
endif

!------------------------------------ find crossing with last fluxsurface 
do j=1,n_tht

  if (nwpts%Z_sep(j) .le. ES%Z_axis) then
  
    Zbeg = ES%Z_xpoint(1)
    Zend = ES%Z_xpoint(1) + 0.3 * (ES%Z_axis - ES%Z_xpoint(1))
    scale_out_points = (nwpts%Z_sep(j) - Zend)/(Zbeg - Zend)
    if (scale_out_points .lt. 0.d0) scale_out_points = 0.d0
    if (scale_out_points .gt. 1.d0) scale_out_points = 1.d0
    if ( (ES%active_xpoint .eq. LOWER_XPOINT) .or. (ES%active_xpoint .eq. SYMMETRIC_XPOINT) ) then
      if (j .gt. n_tht_mid) then
        tht_SOL = 1.0*PI + (stpts%angle_LowerLeft  - 1.0*PI) * scale_out_points
        i_max = n_flux + n_open + n_outer + n_inner
      else
        if (stpts%angle_LowerRight .gt. PI) then
          tht_SOL = 1.98*PI + (stpts%angle_LowerRight - 1.98*PI) * scale_out_points
        else
          tht_SOL = 0.02*PI + (stpts%angle_LowerRight - 0.02*PI) * scale_out_points
        endif
        i_max = n_flux + n_open + n_outer
      endif
      call find_theta_surface(node_list,element_list,flux_list,i_max,tht_SOL,nwpts%R_sep(j),nwpts%Z_sep(j),i_elm_find,s_find,t_find,i_find)    
    elseif ( (xcase .eq. DOUBLE_NULL) .and. (ES%active_xpoint .eq. UPPER_XPOINT) ) then
      if (j .gt. n_tht_mid) then
        if (stpts%angle_LowerRight .gt. PI) then
          tht_SOL = 1.98*PI + (stpts%angle_LowerRight - 1.98*PI) * scale_out_points
        else
          tht_SOL = 0.02*PI + (stpts%angle_LowerRight - 0.02*PI) * scale_out_points
        endif
        i_max = n_flux + n_open + n_outer
      else
        tht_SOL = 1.0*PI + (stpts%angle_LowerLeft  - 1.0*PI) * scale_out_points
        i_max = n_flux + n_open + n_outer + n_inner
      endif
      call find_theta_surface(node_list,element_list,flux_list,i_max,tht_SOL,nwpts%R_sep(j),nwpts%Z_sep(j),i_elm_find,s_find,t_find,i_find)    
    else    
      i_max = n_flux + n_open
      call find_theta_surface(node_list,element_list,flux_list,i_max,theta_sep(j),ES%R_axis,ES%Z_axis,i_elm_find,s_find,t_find,i_find)    
    endif

    do k=1,i_find
      call interp_RZ(node_list,element_list,i_elm_find(k),s_find(k),t_find(k), &
                     RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss, &
                     ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)

      if(    (xcase .eq. UPPER_XPOINT)   &
        .or. (     ( (ES%active_xpoint .eq. LOWER_XPOINT) .or. (ES%active_xpoint .eq. SYMMETRIC_XPOINT) ) &
             .and. (    ( (RRg1 .gt. ES%R_xpoint(1)) .and. (j.lt.n_tht_mid) ) &
                   .or. ( (RRg1 .lt. ES%R_xpoint(1)) .and. (j.gt.n_tht_mid) ) )  ) &
        .or. (     ( (xcase .eq. DOUBLE_NULL) .and. (ES%active_xpoint .eq. UPPER_XPOINT) ) &
             .and. (    ( (RRg1 .lt. ES%R_xpoint(1)) .and. (j.le.n_tht_mid) ) &
                   .or. ( (RRg1 .gt. ES%R_xpoint(1)) .and. (j.gt.n_tht_mid) ) ) ) ) then

        nwpts%R_max(j)             = RRg1
        nwpts%Z_max(j)             = ZZg1
        nwpts%RR_new(i_max+1,j)    = RRg1
        nwpts%ZZ_new(i_max+1,j)    = ZZg1
        nwpts%ielm_flux(i_max+1,j) = i_elm_find(k)
        nwpts%s_flux(i_max+1,j)    = s_find(k)
        nwpts%t_flux(i_max+1,j)    = t_find(k)
        nwpts%t_tht(i_max+1,j)     = 1.d0
        nwpts%k_cross(i_max+1,j)   = 3
        exit
      endif
    enddo
    
    if (     ( (ES%active_xpoint .eq. LOWER_XPOINT) .or. (ES%active_xpoint .eq. SYMMETRIC_XPOINT) ) &
       .and. ((j .eq. 1) .or. (j .eq. n_tht))                                                   ) then
      nwpts%R_max(1)              = stpts%RLimit_LowerOuterLeg ! this one is known - safer...
      nwpts%Z_max(1)              = stpts%ZLimit_LowerOuterLeg ! this one is known - safer...
      nwpts%RR_new(i_max+1,1)     = stpts%RLimit_LowerOuterLeg ! this one is known - safer...
      nwpts%ZZ_new(i_max+1,1)     = stpts%ZLimit_LowerOuterLeg ! this one is known - safer...
      nwpts%R_max(n_tht)          = stpts%RLimit_LowerInnerLeg ! this one is known - safer...
      nwpts%Z_max(n_tht)          = stpts%ZLimit_LowerInnerLeg ! this one is known - safer...
      nwpts%RR_new(i_max+1,n_tht) = stpts%RLimit_LowerInnerLeg ! this one is known - safer...
      nwpts%ZZ_new(i_max+1,n_tht) = stpts%ZLimit_LowerInnerLeg ! this one is known - safer...
    endif
    if (     ( (xcase .eq. DOUBLE_NULL) .and. (ES%active_xpoint .eq. UPPER_XPOINT) )  &
       .and. ((j .eq. 1) .or. (j .eq. n_tht))                            ) then
      nwpts%R_max(n_tht_mid)            = stpts%RLimit_LowerInnerLeg ! this one is known - safer...
      nwpts%Z_max(n_tht_mid)            = stpts%ZLimit_LowerInnerLeg ! this one is known - safer...
      nwpts%RR_new(i_max+1,n_tht_mid)   = stpts%RLimit_LowerInnerLeg ! this one is known - safer...
      nwpts%ZZ_new(i_max+1,n_tht_mid)   = stpts%ZLimit_LowerInnerLeg ! this one is known - safer...
      nwpts%R_max(n_tht_mid+1)          = stpts%RLimit_LowerOuterLeg ! this one is known - safer...
      nwpts%Z_max(n_tht_mid+1)          = stpts%ZLimit_LowerOuterLeg ! this one is known - safer...
      nwpts%RR_new(i_max+1,n_tht_mid+1) = stpts%RLimit_LowerOuterLeg ! this one is known - safer...
      nwpts%ZZ_new(i_max+1,n_tht_mid+1) = stpts%ZLimit_LowerOuterLeg ! this one is known - safer...
    endif

  else
        
    Zbeg = ES%Z_xpoint(2)
    Zend = ES%Z_xpoint(2) + 0.3 * (ES%Z_axis - ES%Z_xpoint(2))
    scale_out_points = (nwpts%Z_sep(j) - Zend)/(Zbeg - Zend)
    if (scale_out_points .lt. 0.d0) scale_out_points = 0.d0
    if (scale_out_points .gt. 1.d0) scale_out_points = 1.d0
    if (    ( (j .gt. n_tht_mid) .and. (xcase .eq. DOUBLE_NULL ) .and. ( (ES%active_xpoint .eq. LOWER_XPOINT) .or. (ES%active_xpoint .eq. SYMMETRIC_XPOINT) ) ) & 
       .or. ( (j .lt. n_tht_mid) .and. (xcase .eq. DOUBLE_NULL ) .and. (  ES%active_xpoint .eq. UPPER_XPOINT                                         ) ) & 
       .or. ( (j .lt. n_tht_mid) .and. (xcase .eq. UPPER_XPOINT) )                               ) then
      tht_SOL = 1.0*PI + (stpts%angle_UpperLeft  - 1.0*PI) * scale_out_points
      i_max = n_flux + n_open + n_outer + n_inner
      call find_theta_surface(node_list,element_list,flux_list,i_max,tht_SOL,nwpts%R_sep(j),nwpts%Z_sep(j),i_elm_find,s_find,t_find,i_find)    
    endif 
    if (    ( (j .le. n_tht_mid) .and. (xcase .eq. DOUBLE_NULL ) .and. ( (ES%active_xpoint .eq. LOWER_XPOINT) .or. (ES%active_xpoint .eq. SYMMETRIC_XPOINT) ) ) &
       .or. ( (j .gt. n_tht_mid) .and. (xcase .eq. DOUBLE_NULL ) .and. (  ES%active_xpoint .eq. UPPER_XPOINT                                         ) ) & 
       .or. ( (j .gt. n_tht_mid) .and. (xcase .eq. UPPER_XPOINT) )                               ) then
      if (stpts%angle_UpperRight .gt. PI) then
        tht_SOL = 1.98*PI + (stpts%angle_UpperRight - 1.98*PI) * scale_out_points
      else
        tht_SOL = 0.02*PI + (stpts%angle_UpperRight - 0.02*PI) * scale_out_points
      endif
      i_max = n_flux + n_open + n_outer
      call find_theta_surface(node_list,element_list,flux_list,i_max,tht_SOL,nwpts%R_sep(j),nwpts%Z_sep(j),i_elm_find,s_find,t_find,i_find)    
    endif
    if (xcase .eq. LOWER_XPOINT) then
      call find_theta_surface(node_list,element_list,flux_list,i_max,theta_sep(j),ES%R_axis,ES%Z_axis,i_elm_find,s_find,t_find,i_find)
    endif

    do k=1,i_find
      call interp_RZ(node_list,element_list,i_elm_find(k),s_find(k),t_find(k),&
                     RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss,    &
                     ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)

      if(    (xcase .eq. LOWER_XPOINT)     &
        .or. (     ( (xcase .eq. DOUBLE_NULL) .and. ( (ES%active_xpoint .eq. LOWER_XPOINT) .or. (ES%active_xpoint .eq. SYMMETRIC_XPOINT) ) )   &
             .and. (    ( (RRg1 .ge. ES%R_xpoint(2)) .and. (j .le. n_tht_mid) ) &
                   .or. ( (RRg1 .lt. ES%R_xpoint(2)) .and. (j .gt. n_tht_mid) ) )  ) &
        .or. (     ( ES%active_xpoint .eq. UPPER_XPOINT )   &
             .and. (    ( (RRg1 .ge. ES%R_xpoint(2)) .and. (j .gt. n_tht_mid) ) &
                   .or. ( (RRg1 .lt. ES%R_xpoint(2)) .and. (j .lt. n_tht_mid) ) )  ) ) then

        nwpts%R_max(j)             = RRg1
        nwpts%Z_max(j)             = ZZg1
        nwpts%RR_new(i_max+1,j)    = RRg1
        nwpts%ZZ_new(i_max+1,j)    = ZZg1
        nwpts%ielm_flux(i_max+1,j) = i_elm_find(k)
        nwpts%s_flux(i_max+1,j)    = s_find(k)
        nwpts%t_flux(i_max+1,j)    = t_find(k)
        nwpts%t_tht(i_max+1,j)     = 1.d0
        nwpts%k_cross(i_max+1,j)   = 3
      endif
    enddo
    
    if (     ( ES%active_xpoint .eq. UPPER_XPOINT ) & 
       .and. ( (j .eq. 1) .or. (j .eq. n_tht_2) )                         ) then
      nwpts%R_max(n_tht)            = stpts%RLimit_UpperOuterLeg ! this one is known - safer...
      nwpts%Z_max(n_tht)            = stpts%ZLimit_UpperOuterLeg ! this one is known - safer...
      nwpts%RR_new(i_max+1,n_tht)   = stpts%RLimit_UpperOuterLeg ! this one is known - safer...
      nwpts%ZZ_new(i_max+1,n_tht)   = stpts%ZLimit_UpperOuterLeg ! this one is known - safer...
      nwpts%R_max(1)                = stpts%RLimit_UpperInnerLeg ! this one is known - safer...
      nwpts%Z_max(1)                = stpts%ZLimit_UpperInnerLeg ! this one is known - safer...
      nwpts%RR_new(i_max+1,1)       = stpts%RLimit_UpperInnerLeg ! this one is known - safer...
      nwpts%ZZ_new(i_max+1,1)       = stpts%ZLimit_UpperInnerLeg ! this one is known - safer...
    endif
    
    if (     ( (xcase .eq. DOUBLE_NULL) .and. ( (ES%active_xpoint .eq. LOWER_XPOINT) .or. (ES%active_xpoint .eq. SYMMETRIC_XPOINT) ) ) & 
       .and. ( (j .eq. n_tht_mid) .or. (j .eq. n_tht_mid+1) )            ) then
      nwpts%R_max(n_tht_mid)            = stpts%RLimit_UpperOuterLeg ! this one is known - safer...
      nwpts%Z_max(n_tht_mid)            = stpts%ZLimit_UpperOuterLeg ! this one is known - safer...
      nwpts%RR_new(i_max+1,n_tht_mid)   = stpts%RLimit_UpperOuterLeg ! this one is known - safer...
      nwpts%ZZ_new(i_max+1,n_tht_mid)   = stpts%ZLimit_UpperOuterLeg ! this one is known - safer...
      nwpts%R_max(n_tht_mid+1)          = stpts%RLimit_UpperInnerLeg ! this one is known - safer...
      nwpts%Z_max(n_tht_mid+1)          = stpts%ZLimit_UpperInnerLeg ! this one is known - safer...
      nwpts%RR_new(i_max+1,n_tht_mid+1) = stpts%RLimit_UpperInnerLeg ! this one is known - safer...
      nwpts%ZZ_new(i_max+1,n_tht_mid+1) = stpts%ZLimit_UpperInnerLeg ! this one is known - safer...
    endif

  endif

enddo







!------------------------------------------------------------------------------------------------------------------------!
!************************************************************************************************************************!
!************************************************************************************************************************!
!*********************************** Second part: find crossings of lines ***********************************************!
!************************************************************************************************************************!
!************************************************************************************************************************!
!------------------------------------------------------------------------------------------------------------------------!
write(*,*) '                 Find crossings between coordinate lines'





!--------------------------------------------------------------------------!
!--------------- First construct polar coordinate lines -------------------!
!--------------------------------------------------------------------------!

allocate (R_beg (n_tht),Z_beg (n_tht))
allocate (delta1(n_tht),delta2(n_tht))
allocate (R_polar(5,4,n_tht),Z_polar(5,4,n_tht))

! --- These deltas determine how close to the mid-point you place your two middle control points, 10% by default.
! --- But for example, close to the X-point, you want them to be very close to the X-lines, so you go to a smaller values.
! --- On the X-lines themselves, of course you want the polar lines to be exactly on the X-point, and so it needs to be zero.
delta1 = 0.1
delta2 = 0.1
delta2(1) = 0.d0   ; delta2(n_tht  ) = 0.d0
delta2(2) = 0.05d0 ; delta2(n_tht-1) = 0.05d0
delta2(n_tht_mid  ) = 0.d0   ; delta2(n_tht_mid+1) = 0.d0
delta2(n_tht_mid-1) = 0.05d0 ; delta2(n_tht_mid+2) = 0.05d0

do j=1,n_tht
  R_beg (j) = ES%R_axis
  Z_beg (j) = ES%Z_axis
enddo

call create_polar_lines_4(n_tht,       R_beg,       Z_beg, &
                                 nwpts%R_mid, nwpts%Z_mid, &
                                 nwpts%R_sep, nwpts%Z_sep, &
                                 nwpts%R_max, nwpts%Z_max, &
                                 delta1, delta2, R_polar, Z_polar)
do i=1,5
  do j=1,n_tht
    nwpts%R_polar(i,1:4,j) = R_polar(i,1:4,j)
    nwpts%Z_polar(i,1:4,j) = Z_polar(i,1:4,j)
  enddo
enddo

deallocate (R_beg ,Z_beg )
deallocate (delta1,delta2)







!--------------------------------------------------------------------------!
!--------- Find grid_points from crossing of coordinate lines -------------!
!--------------------------------------------------------------------------!


!----------------------------------- The magnetic axis
do j=1, n_tht
  nwpts%RR_new(1,j)    = ES%R_axis
  nwpts%ZZ_new(1,j)    = ES%Z_axis
  nwpts%ielm_flux(1,j) = ES%i_elm_axis
  nwpts%s_flux(1,j)    = ES%s_axis
  nwpts%t_flux(1,j)    = ES%t_axis
  nwpts%t_tht(1,j)     = -1.d0          ! expressed in cubic Hermite (-1<t<+1)
enddo
nwpts%k_cross(1,:) = 1

!----------------------------------- The main part (without legs, outer and inner)
do i=1,n_flux+n_open+n_outer+n_inner
  n_start = 1
  n_loop  = n_tht
  if ( (i .gt. n_flux+n_open) .and. (i .le. n_flux+n_open+n_outer) ) n_start = 1
  if ( (i .gt. n_flux+n_open) .and. (i .le. n_flux+n_open+n_outer) ) n_loop  = n_tht_mid
  if ( (i .gt. n_flux+n_open+n_outer) ) n_start = n_tht_mid+1
  if ( (i .gt. n_flux+n_open+n_outer) ) n_loop  = n_tht
  do j=n_start, n_loop
    
    do k=1,5     ! 5 line pieces per coordinate line

      R_cub1d = (/ R_polar(k,1,j), 3.d0/2.d0 *(R_polar(k,2,j)-R_polar(k,1,j)), &
                   R_polar(k,4,j), 3.d0/2.d0 *(R_polar(k,4,j)-R_polar(k,3,j))  /)
      Z_cub1d = (/ Z_polar(k,1,j), 3.d0/2.d0 *(Z_polar(k,2,j)-Z_polar(k,1,j)), &
                   Z_polar(k,4,j), 3.d0/2.d0 *(Z_polar(k,4,j)-Z_polar(k,3,j)) /)

      call find_crossing(node_list,element_list,flux_list,i,R_cub1d,Z_cub1d, &
                         nwpts%RR_new(i+1,j),nwpts%ZZ_new(i+1,j),nwpts%ielm_flux(i+1,j),&
                         nwpts%s_flux(i+1,j),nwpts%t_flux(i+1,j),nwpts%t_tht(i+1,j),ifail, .true.)
      ! --- Readjust to make sure we are inside element.
      if ( (nwpts%s_flux(i+1,j) .lt. 0.d0) .or. (nwpts%s_flux(i+1,j) .gt. 1.d0) .or. (nwpts%t_flux(i+1,j) .lt. 0.d0) .or. (nwpts%t_flux(i+1,j) .gt. 1.d0) ) then
        if (nwpts%s_flux(i+1,j) .lt. 0.d0) nwpts%s_flux(i+1,j) = 0.d0
        if (nwpts%s_flux(i+1,j) .gt. 1.d0) nwpts%s_flux(i+1,j) = 1.d0
        if (nwpts%t_flux(i+1,j) .lt. 0.d0) nwpts%t_flux(i+1,j) = 0.d0
        if (nwpts%t_flux(i+1,j) .gt. 1.d0) nwpts%t_flux(i+1,j) = 1.d0
        call interp_RZ(node_list,element_list,nwpts%ielm_flux(i+1,j),nwpts%s_flux(i+1,j),nwpts%t_flux(i+1,j),&
                       RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss,    &
                       ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)
        nwpts%RR_new(i+1,j) = RRg1
        nwpts%ZZ_new(i+1,j) = ZZg1
      endif
      
      if (ifail .eq. 0) then
        nwpts%k_cross(i+1,j) = k
        exit
      endif

    enddo

    if (ifail .ne. 0) then
      write(*,'(A,I6,I6,I6,F20.10)') ' WARNING node not found for central grid (without legs) : ',ifail,i,j
    endif
      
  enddo
enddo











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



!-------------------------------------------------------------------------------------------!
!-------------------------------- The closed field lines -----------------------------------!
!-------------------------------------------------------------------------------------------!

n_loop = n_tht-1
if (xcase .eq. DOUBLE_NULL) n_loop = n_loop-1  ! For double-null, n_tht_mid and n_tht_mid+1 are the same lines
do i=1,n_flux                 
  do k=1, n_loop

    j = k
    if ((xcase .eq. DOUBLE_NULL) .and. (j .gt. n_tht_mid)) j = j+1  ! For double-null, n_tht_mid and n_tht_mid+1 are the same lines
    
    index = index + 1
    call create_new_node(node_list, element_list, newnode_list, index, i, j, nwpts)

    if (i .eq. 1) then   !------------------------------------ magnetic axis : special case
      newnode_list%node(index)%x(1,3,:) = 0.d0
      newnode_list%node(index)%x(1,4,:) = 0.d0
    endif

  enddo
enddo
newnode_list%n_nodes = index



!-------------------------------------------------------------------------------------------!
!--------------- The open part (or sandwich part, in case of double-null) ------------------!
!-------------------------------------------------------------------------------------------!

       !- (the one around the main plasma)
       !- Starting from right strike point, all around to left strike point (if lower Xpoint is the main one)
       !- The other way round if upper Xpoint is the main one
       !- We add 4 points for each Xpoint, either on the left or on the right 
       !- (see routine create_x_node for more info)
       
n_start_open = newnode_list%n_nodes + 1
if (xcase .eq. LOWER_XPOINT) n_loop = n_tht 
if (xcase .eq. UPPER_XPOINT) n_loop = n_tht
if ( (xcase .eq. DOUBLE_NULL) .and. ( (ES%active_xpoint .eq. LOWER_XPOINT) .or. (ES%active_xpoint .eq. SYMMETRIC_XPOINT) ) ) n_loop = n_tht - 1 
if ( (xcase .eq. DOUBLE_NULL) .and. (  ES%active_xpoint .eq. UPPER_XPOINT )                                         ) n_loop = n_tht - 1
n_loop2 = n_flux+n_open+1
if (xcase .eq. DOUBLE_NULL) n_loop2 = n_flux+n_open
n_tmp = n_tht
if (xcase .eq. DOUBLE_NULL) n_tmp = n_tht-1  ! n_tht_mid and n_tht_mid+1 are the same for double null          
if (ES%active_xpoint .ne. SYMMETRIC_XPOINT ) then ! ignore if symmetric double-null
  
  do i=n_flux+1,n_loop2
    
    do l=1,n_loop
    
      j = l
      !-------------------------------- CASE 1 : Lower Xpoint is the main one
      if ( (xcase .eq. LOWER_XPOINT) .or. (xcase .eq. UPPER_XPOINT) .or. ( (xcase .eq. DOUBLE_NULL) .and. ( (ES%active_xpoint .eq. LOWER_XPOINT) .or. (ES%active_xpoint .eq. SYMMETRIC_XPOINT) ) ) ) then
        if ( ((j .gt. 1) .and. (j .lt. n_tmp)) .or. (i .ne. n_flux+1) ) then ! Don't put the Xpoint twice
          if ( (xcase .eq. DOUBLE_NULL) .and. (j .gt. n_tht_mid) ) j = j+1  ! n_tht_mid and n_tht_mid+1 are the same for double null
          index = index + 1
          call create_new_node(node_list, element_list, newnode_list, index, i, j, nwpts)    
        endif    
      endif
      
    enddo
  enddo
endif
newnode_list%n_nodes = index





!-------------------------------------------------------------------------------------------!
!------------------------ The outer part, in case of double-null ---------------------------!
!-------------------------------------------------------------------------------------------!
       
n_start_outer = newnode_list%n_nodes + 1
if (xcase .eq. DOUBLE_NULL) then
  if ( (ES%active_xpoint .eq. LOWER_XPOINT) .or. (ES%active_xpoint .eq. SYMMETRIC_XPOINT) ) n_loop = n_tht_mid 
  if (  ES%active_xpoint .eq. UPPER_XPOINT                                         ) n_loop = n_tht_mid2
  do i=n_flux+n_open+1,n_flux+n_open+n_outer+1
    do l=1,n_loop
    
      j = l
      if ( (ES%active_xpoint .eq. LOWER_XPOINT) .or. (ES%active_xpoint .eq. SYMMETRIC_XPOINT) ) then 
        if ( (j .lt. n_tmp) .or. (i .ne. n_flux+n_open+1) ) then ! Don't put the Xpoint twice
          index = index + 1
          call create_new_node(node_list, element_list, newnode_list, index, i, j, nwpts)
        endif    
      else 
        if ( (j .gt. 1) .or. (i .ne. n_flux+n_open+1) ) then ! Don't put the Xpoint twice
          j = j+n_tht_mid
          index = index + 1
          call create_new_node(node_list, element_list, newnode_list, index, i, j, nwpts)    
        endif  
      endif  
      
    enddo
  enddo
  newnode_list%n_nodes = index
endif





!-------------------------------------------------------------------------------------------!
!------------------------ The inner part, in case of double-null ---------------------------!
!-------------------------------------------------------------------------------------------!
       
n_start_inner = newnode_list%n_nodes + 1
if (xcase .eq. DOUBLE_NULL) then
  if ( (ES%active_xpoint .eq. LOWER_XPOINT) .or. (ES%active_xpoint .eq. SYMMETRIC_XPOINT) ) n_loop = n_tht_mid2
  if (  ES%active_xpoint .eq. UPPER_XPOINT                                         ) n_loop = n_tht_mid
  do k=n_flux+n_open+n_outer+1,n_flux+n_open+n_outer+n_inner+1                
    do l=1,n_loop
    
      i = k
      if (i .eq. n_flux+n_open+n_outer+1) i = n_flux+n_open+1  !Doing the second separatrix first
      
      j = l
      if ( (ES%active_xpoint .eq. LOWER_XPOINT) .or. (ES%active_xpoint .eq. SYMMETRIC_XPOINT) ) then
        if ( (j .gt. 1) .or. (i .ne. n_flux+n_open+1) ) then ! Don't put the Xpoint twice
          j = j+n_tht_mid
          index = index + 1
          call create_new_node(node_list, element_list, newnode_list, index, i, j, nwpts)    
        endif  
      else
        if ( (j .lt. n_tmp) .or. (i .ne. n_flux+n_open+1) ) then ! Don't put the Xpoint twice
          index = index + 1
          call create_new_node(node_list, element_list, newnode_list, index, i, j, nwpts)    
        endif  
      endif  
      
    enddo
  enddo
  newnode_list%n_nodes = index
endif















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
n_loop = n_tht-1
if (xcase .eq. DOUBLE_NULL) n_loop = n_tht-2
index = 0
do i=1,n_flux
  do l=1, n_loop

    j = l
    index = index + 1
    newelement_list%element(index)%size = 1.d0

    newelement_list%element(index)%vertex(1) = n_tmp + (i-1)*n_loop + j
    newelement_list%element(index)%vertex(2) = n_tmp + (i  )*n_loop + j
    newelement_list%element(index)%vertex(3) = n_tmp + (i  )*n_loop + j + 1
    newelement_list%element(index)%vertex(4) = n_tmp + (i-1)*n_loop + j + 1

    if (j .eq. n_loop) then
      newelement_list%element(index)%vertex(4)  = n_tmp + (i-1)*n_loop + 1
      newelement_list%element(index)%vertex(3)  = n_tmp + (i  )*n_loop + 1
    endif

    ! Connect with open (or sandwich) region (or outer and inner in case of symmetric double-null)
    if (i .eq. n_flux) then
      if ( (ES%active_xpoint .eq. LOWER_XPOINT) .or. (ES%active_xpoint .eq. SYMMETRIC_XPOINT) ) then
        n_xpoint_1 = n_start_open -1
      else 
        n_xpoint_1 = n_start_open -1
      endif
      if ( ES%active_xpoint .eq. SYMMETRIC_XPOINT) n_xpoint_1 = n_start_outer
      if ((ES%active_xpoint .eq. SYMMETRIC_XPOINT) .and. (l .gt. n_tht_mid-1) ) then
        n_xpoint_1 = n_start_inner -1
        j = j - n_tht_mid+1
      endif
      newelement_list%element(index)%vertex(2) = n_xpoint_1 + j - 1 
      newelement_list%element(index)%vertex(3) = n_xpoint_1 + j
      if (l .eq. 1) then ! Special case for element starting at Xpoint
        newelement_list%element(index)%vertex(2) = 2 ! We need to start at NODE2 of the lower Xpoint (or NODE6=NODE2 of the upper Xpoint)
      endif
      if (l .eq. n_loop) then ! Special case for element arriving at Xpoint
        newelement_list%element(index)%vertex(3) = 3 ! We need to finish at NODE3 of the lower Xpoint (or NODE7=NODE3 of the upper Xpoint)
      endif
      if ( (ES%active_xpoint .eq. SYMMETRIC_XPOINT) .and. (l .eq. n_tht_mid-1) ) then ! Special case for element passing through 2nd Xpoint (for symmetric only)
        newelement_list%element(index)%vertex(3) = 7 ! We need to arrive at NODE7 of the upper Xpoint
      endif
      if ( (ES%active_xpoint .eq. SYMMETRIC_XPOINT) .and. (l .eq. n_tht_mid  ) ) then ! Special case for element passing through 2nd Xpoint (for symmetric only)
        newelement_list%element(index)%vertex(2) = 6 ! We need to leave from NODE6 of the upper Xpoint
      endif
    endif  
      
  enddo
enddo
newelement_list%n_elements = index


!-------------------------------- The open (or sandwich) region (between the two separatrices)
if (xcase .eq. LOWER_XPOINT) then
  n_xpoint_1 = 0               ! First time through Xpoint
  n_xpoint_2 = n_tht-1         ! Second time through Xpoint
  n_loop     = n_tht-1
endif 
if ((xcase .eq. DOUBLE_NULL) .and. ( (ES%active_xpoint .eq. LOWER_XPOINT) .or. (ES%active_xpoint .eq. SYMMETRIC_XPOINT) )) then
  n_xpoint_1 = 0               ! First time through first Xpoint
  n_xpoint_2 = n_tht-2         ! Second time through first Xpoint
  n_xpoint_3 = n_tht_mid-1     ! Going through second Xpoint
  n_loop     = n_tht-2
endif 
if ((xcase .eq. DOUBLE_NULL) .and. ( ES%active_xpoint .eq. UPPER_XPOINT                                          )) then
  n_xpoint_1 = 0            ! First time through first Xpoint
  n_xpoint_2 = n_tht-2      ! Second time through first Xpoint
  n_xpoint_3 = n_tht_mid-1  ! Going through second Xpoint
  n_loop     = n_tht-2
endif 
if (xcase .eq. UPPER_XPOINT) then
  n_xpoint_1 = 0            ! First time through first Xpoint
  n_xpoint_2 = n_tht-1      ! Second time through first Xpoint
  n_loop     = n_tht-1
endif 
if (ES%active_xpoint .ne. SYMMETRIC_XPOINT) then ! ignore if symmetric double-null
  do i=1,n_open
    do l=1, n_loop

      j = l
      j2 = l
      index = index + 1
      newelement_list%element(index)%size = 1.d0

      ! Note the -2 stands for the two Xpoints at the separatrix
      ! Note that for n_loop elements, there are n_loop+1 nodes
      newelement_list%element(index)%vertex(1) = n_start_open + (i-1)*(n_loop+1)-2 + j - 1
      newelement_list%element(index)%vertex(4) = n_start_open + (i-1)*(n_loop+1)-2 + j
      newelement_list%element(index)%vertex(2) = n_start_open + (i  )*(n_loop+1)-2 + j - 1
      newelement_list%element(index)%vertex(3) = n_start_open + (i  )*(n_loop+1)-2 + j
      
      ! Connect with closed region
      if (i .eq. 1) then
        if (l .gt. n_xpoint_1) j = j-1 ! First time through Xpoint
        if (l .gt. n_xpoint_2) j = j-1 ! Second time through Xpoint
        newelement_list%element(index)%vertex(1) = n_start_open + j - 1
        newelement_list%element(index)%vertex(4) = n_start_open + j
        if ( (l .eq. n_xpoint_1) .or. (l .eq. n_xpoint_2) ) then ! Special case for element arriving at Xpoint
          newelement_list%element(index)%vertex(4) = 4 ! We need to arrive at NODE8=NODE4 of the lower Xpoint (or NODE4 of the upper Xpoint)
        endif
        if ( (l .eq. n_xpoint_1+1) .or. (l .eq. n_xpoint_2+1) ) then ! Special case for element leaving Xpoint
          newelement_list%element(index)%vertex(1) = 1 ! We need to leave at NODE1 of the lower Xpoint (or NODE5=NODE1 of the upper Xpoint)
        endif
      endif

      ! Connect with outer and inner regions
      if ( (i .eq. n_open) .and. (xcase .eq. DOUBLE_NULL) ) then
        if ( (ES%active_xpoint .eq. LOWER_XPOINT) .or. (ES%active_xpoint .eq. SYMMETRIC_XPOINT) ) then
          if (l .le. n_xpoint_3) then
            newelement_list%element(index)%vertex(2) = n_start_outer + j2 - 1
            newelement_list%element(index)%vertex(3) = n_start_outer + j2
          else
            j2 = j2-n_xpoint_3
            newelement_list%element(index)%vertex(2) = n_start_inner - 1 + j2 - 1
            newelement_list%element(index)%vertex(3) = n_start_inner - 1 + j2
          endif
        else
          if (l .le. n_xpoint_3) then
            newelement_list%element(index)%vertex(2) = n_start_inner + j2 - 1
            newelement_list%element(index)%vertex(3) = n_start_inner + j2
          else
            j2 = j2-n_xpoint_3
            newelement_list%element(index)%vertex(2) = n_start_outer - 1 + j2 - 1
            newelement_list%element(index)%vertex(3) = n_start_outer - 1 + j2
          endif
        endif
        if (l .eq. n_xpoint_3) then ! Special case for element arriving at second Xpoint
          newelement_list%element(index)%vertex(3) = 7 ! We need to arrive at NODE7 of the upper Xpoint (or NODE3=NODE7 of the lower Xpoint)
        endif
        if (l .eq. n_xpoint_3+1) then ! Special case for element leaving second Xpoint
          newelement_list%element(index)%vertex(2) = 6 ! We need to leave at NODE6 of the upper Xpoint (or NODE2=NODE6 of the lower Xpoint)
        endif
      endif

    enddo
  enddo
endif
newelement_list%n_elements = index


!-------------------------------- The outer region
if (xcase .eq. DOUBLE_NULL) then
  if ( (ES%active_xpoint .eq. LOWER_XPOINT) .or. (ES%active_xpoint .eq. SYMMETRIC_XPOINT) ) then
    n_xpoint_1 = (n_tht_mid-1)
    n_loop     = (n_tht_mid-1)
  else
    n_xpoint_1 = 0
    n_loop     = (n_tht_mid2-1)
  endif 
  if (ES%active_xpoint .eq. SYMMETRIC_XPOINT) then
    n_xpoint_1 = (n_tht_mid-1)
    n_xpoint_2 = 0
    n_loop     = (n_tht_mid-1)
  endif 
  do i=1,n_outer
    do l=1, n_loop

      j = l
      index = index + 1
      newelement_list%element(index)%size = 1.d0

      newelement_list%element(index)%vertex(1) = n_start_outer + (i-1)*(n_loop+1)-1 + j
      newelement_list%element(index)%vertex(4) = n_start_outer + (i-1)*(n_loop+1)-1 + j + 1
      newelement_list%element(index)%vertex(2) = n_start_outer + (i  )*(n_loop+1)-1 + j
      newelement_list%element(index)%vertex(3) = n_start_outer + (i  )*(n_loop+1)-1 + j + 1
      
      ! Connect with open region (or closed and private for symmetric double-null)
      if (i .eq. 1) then
        newelement_list%element(index)%vertex(1) = n_start_outer + j - 1
        newelement_list%element(index)%vertex(4) = n_start_outer + j
        if (l .eq. n_xpoint_1)   then ! Special case for element arriving at second Xpoint
          newelement_list%element(index)%vertex(4) = 8 ! We need to arrive at NODE8=NODE4 of the second Xpoint
        endif
        if ( (l .eq. n_xpoint_2+1) .and. (ES%active_xpoint .eq. SYMMETRIC_XPOINT) ) then ! Special case for element leaving first Xpoint of symmetric double-null
          newelement_list%element(index)%vertex(1) = 1 ! We need to leave at NODE1 of the first Xpoint
        endif
      endif

    enddo
  enddo
endif
newelement_list%n_elements = index


!-------------------------------- The inner region
if (xcase .eq. DOUBLE_NULL) then
  if ( (ES%active_xpoint .eq. LOWER_XPOINT) .or. (ES%active_xpoint .eq. SYMMETRIC_XPOINT) ) then
    n_xpoint_1 = 0
    n_loop     = (n_tht_mid2-1)
  else
    n_xpoint_1 = (n_tht_mid-1)
    n_loop     = (n_tht_mid-1)
  endif 
  if (ES%active_xpoint .eq. SYMMETRIC_XPOINT) then
    n_xpoint_1 = 0
    n_xpoint_2 = (n_tht_mid2-1)
    n_loop     = (n_tht_mid2-1)
  endif 
  do i=1,n_inner
    do l=1, n_loop

      j = l
      index = index + 1
      newelement_list%element(index)%size = 1.d0

      newelement_list%element(index)%vertex(1) = n_start_inner + (i-1)*(n_loop+1)-1 + j - 1
      newelement_list%element(index)%vertex(4) = n_start_inner + (i-1)*(n_loop+1)-1 + j
      newelement_list%element(index)%vertex(2) = n_start_inner + (i  )*(n_loop+1)-1 + j - 1
      newelement_list%element(index)%vertex(3) = n_start_inner + (i  )*(n_loop+1)-1 + j
      
      ! Connect with open region
      if (i .eq. 1) then
        if (j .gt. n_xpoint_1) j = j-1
        newelement_list%element(index)%vertex(1) = n_start_inner + j - 1
        newelement_list%element(index)%vertex(4) = n_start_inner + j
        if (l .eq. n_xpoint_1) then ! Special case for element arriving at second Xpoint
          newelement_list%element(index)%vertex(4) = 8 ! We need to arrive at NODE4 of the second Xpoint
        endif
        if (l .eq. n_xpoint_1+1) then! Special case for element leaving second Xpoint
          newelement_list%element(index)%vertex(1) = 5 ! We need to leave at NODE5=NODE1 of the second Xpoint
        endif
        if ( (l .eq. n_xpoint_2)   .and. (ES%active_xpoint .eq. SYMMETRIC_XPOINT) ) then ! Special case for element arriving at first Xpoint of symmetric double-null
          newelement_list%element(index)%vertex(4) = 4 ! We need to arrive at NODE4 of the first Xpoint
        endif
        if ( (l .eq. n_xpoint_2+1) .and. (ES%active_xpoint .eq. SYMMETRIC_XPOINT) ) then ! Special case for element leaving first Xpoint of symmetric double-null
          newelement_list%element(index)%vertex(1) = 1 ! We need to leave at NODE5=NODE1 of the first Xpoint
        endif
      endif

    enddo
  enddo
endif
newelement_list%n_elements = index


!----------------------------------- Print a python file that plots a cross with the 4 nodes of each element
if (plot_grid) then
  n_loop = newelement_list%n_elements
  open(101,file='plot_central_grid_elements.py')
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



deallocate(theta_sep, theta_beg)
deallocate (R_polar,Z_polar)

return
end subroutine define_central_grid
