!> Subroutine defines the new grid_points from crossing of polar and radial coordinate lines
subroutine find_extrapolation_points_central_part(node_list, element_list, flux_list, xcase,   &
                                                  n_grids, stpts, sigmas, nwpts)

use constants
use tr_module 
use data_structure
use grid_xpoint_data
use phys_module, only:   tokamak_device
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

! --- local variables
real*8, allocatable :: s_tmp(:), theta_sep(:)
real*8, allocatable :: xp(:),yp(:)
integer             :: i, j, k, m, i_elm_find(8), i_find
integer             :: i_sep1, i_sep2, i_max, n_loop, n_start, i_surf
integer             :: n_psi, n_tht_2, n_tht_mid, n_tht_mid2
integer             :: n_flux, n_tht,   n_open,   n_outer,   n_inner    
integer             :: n_private,   n_up_priv,   n_leg,   n_up_leg
integer             :: npl, ifail, my_id
real*8              :: delta, ss, tmp1, tmp2
real*8              :: R_cub1d(4), Z_cub1d(4)
real*8              :: tht_x1, tht_x2
real*8              :: RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss
real*8              :: ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss
real*8              :: PSg1,dPSg1_dr,dPSg1_ds,dPSg1_drs,dPSg1_drr,dPSg1_dss
real*8              :: s_find(8), t_find(8), st_find(8)
real*8              :: R_beg, R_end
real*8              :: Z_beg, Z_end
real*8              :: SIG_theta, SIG_theta_up
real*8              :: SIG_leg_0, SIG_leg_1
real*8              :: SIG_up_leg_0, SIG_up_leg_1
real*8              :: SIG_0, SIG_1
real*8              :: bgf_tht
logical, parameter  :: plot_grid = .true.

write(*,*) '*****************************************'
write(*,*) '* X-point grid : Define new grid points *'
write(*,*) '*****************************************'
write(*,*) '                 Define extrapolation points'

SIG_theta    = sigmas(2) ; SIG_theta_up = sigmas(17)
SIG_leg_0    = sigmas(8) ; SIG_leg_1    = sigmas(9) 
SIG_up_leg_0 = sigmas(10); SIG_up_leg_1 = sigmas(11)

n_flux    = n_grids(1); n_tht     = n_grids(2)
n_open    = n_grids(3); n_outer   = n_grids(4); n_inner = n_grids(5)
n_private = n_grids(6); n_up_priv = n_grids(7)
n_leg     = n_grids(8); n_up_leg  = n_grids(9)

nwpts%k_cross = 0

bgf_tht  = 0.8d0!0.97d0


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

call tr_allocate(theta_sep,1,n_tht_2,"theta_sep",CAT_GRID)
  
if (xcase .ne. DOUBLE_NULL) then 
  if (xcase .eq. LOWER_XPOINT) tht_x1 = atan2(ES%Z_xpoint(1)-ES%Z_axis,ES%R_xpoint(1)-ES%R_axis)
  if (xcase .eq. UPPER_XPOINT) tht_x1 = atan2(ES%Z_xpoint(2)-ES%Z_axis,ES%R_xpoint(2)-ES%R_axis)
  n_tht_mid = n_tht/2
  
  call tr_allocate(s_tmp,1,n_tht,"s_tmp",CAT_GRID)
  s_tmp = 0
  call meshac2(n_tht,s_tmp,0.d0,1.d0,SIG_theta,SIG_theta,bgf_tht,1.0d0)

  do j=1,n_tht
    theta_sep(j) = tht_x1 + 2.d0 * PI * s_tmp(j)
    if (theta_sep(j) .lt. 0.d0)    theta_sep(j) = theta_sep(j) + 2.d0 * PI
    if (theta_sep(j) .gt. 2.d0*PI) theta_sep(j) = theta_sep(j) - 2.d0 * PI
  enddo
  call tr_deallocate(s_tmp,"s_tmp",CAT_GRID)
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
  !write(*,'(A,2f15.4)') ' angles : ',tht_x1,tht_x2
  
  ! Spread out points evenly (outer angle between tht_x1 and tht_x2 is usually bigger than inner angle)
  if ( (ES%active_xpoint .eq. LOWER_XPOINT ) .or. (ES%active_xpoint .eq. SYMMETRIC_XPOINT) ) then
    n_tht_mid = int(n_tht * (2.d0*PI - (tht_x1 - tht_x2)) / (2.d0*PI))
    ! Make sure n_tht_mid is odd and save it to n_grids for later use
    if(mod(n_tht_mid,2) .eq. 0) n_tht_mid = n_tht_mid + 1
    n_grids(10) = n_tht_mid
    
    call tr_allocate(s_tmp,1,n_tht_mid,"s_tmp",CAT_GRID)
    s_tmp = 0
    call meshac2(n_tht_mid,s_tmp,0.d0,1.d0,SIG_theta,SIG_theta,bgf_tht,1.0d0)

    do j=1,n_tht_mid
      theta_sep(j) = (tht_x1-2.d0*PI) + (tht_x2-(tht_x1-2.d0*PI)) * s_tmp(j)
      if (theta_sep(j) .lt. 0.d0)    theta_sep(j) = theta_sep(j) + 2.d0 * PI
      if (theta_sep(j) .gt. 2.d0*PI) theta_sep(j) = theta_sep(j) - 2.d0 * PI
    enddo
    call tr_deallocate(s_tmp,"s_tmp",CAT_GRID)
  
    n_tht_mid2 = n_tht-n_tht_mid
    call tr_allocate(s_tmp,1,n_tht_mid2,"s_tmp",CAT_GRID)
    s_tmp = 0
    call meshac2(n_tht_mid2,s_tmp,0.d0,1.d0,SIG_theta,SIG_theta,bgf_tht,1.0d0)

    do j=n_tht_mid+1,n_tht
      theta_sep(j) = tht_x2 + (tht_x1-tht_x2) * s_tmp(j-n_tht_mid)
      if (theta_sep(j) .lt. 0.d0)    theta_sep(j) = theta_sep(j) + 2.d0 * PI
      if (theta_sep(j) .gt. 2.d0*PI) theta_sep(j) = theta_sep(j) - 2.d0 * PI
    enddo
    call tr_deallocate(s_tmp,"s_tmp",CAT_GRID)
  else
    n_tht_mid = int(n_tht * (tht_x2 - tht_x1) / (2.d0*PI))
    ! Make sure n_tht_mid is odd and save it to n_grids for later use
    if(mod(n_tht_mid,2) .eq. 0) n_tht_mid = n_tht_mid + 1
    n_grids(10) = n_tht_mid
    
    call tr_allocate(s_tmp,1,n_tht_mid,"s_tmp",CAT_GRID)
    s_tmp = 0
    call meshac2(n_tht_mid,s_tmp,0.d0,1.d0,SIG_theta,SIG_theta,bgf_tht,1.0d0)
    do j=1,n_tht_mid
      theta_sep(j) = tht_x1 + (tht_x2-tht_x1) * s_tmp(j)
      if (theta_sep(j) .lt. 0.d0)    theta_sep(j) = theta_sep(j) + 2.d0 * PI
      if (theta_sep(j) .gt. 2.d0*PI) theta_sep(j) = theta_sep(j) - 2.d0 * PI
    enddo
    call tr_deallocate(s_tmp,"s_tmp",CAT_GRID)

    n_tht_mid2 = n_tht-n_tht_mid
    call tr_allocate(s_tmp,1,n_tht_mid2,"s_tmp",CAT_GRID)
    s_tmp = 0
    call meshac2(n_tht_mid2,s_tmp,0.d0,1.d0,SIG_theta,SIG_theta,bgf_tht,1.0d0)
    do j=n_tht_mid+1,n_tht
      theta_sep(j) = (tht_x2-2.d0*PI) + (tht_x1-(tht_x2-2.d0*PI)) * s_tmp(j-n_tht_mid)
      if (theta_sep(j) .lt. 0.d0)    theta_sep(j) = theta_sep(j) + 2.d0 * PI
      if (theta_sep(j) .gt. 2.d0*PI) theta_sep(j) = theta_sep(j) - 2.d0 * PI
    enddo
    call tr_deallocate(s_tmp,"s_tmp",CAT_GRID)
  
  endif      
endif

!------------------------------------- find crossing with separatrix
i_sep1  = n_flux
i_sep2  = n_flux + n_open
if (xcase .ne. DOUBLE_NULL) then
  do j=2,n_tht-1
    call find_theta_surface(node_list,element_list,flux_list,i_sep1,theta_sep(j),ES%R_axis,ES%Z_axis,i_elm_find,s_find,t_find,i_find)
    if(i_find .eq. 0) return
    
    call interp_RZ(node_list,element_list,i_elm_find(1),s_find(1),t_find(1),&
                   RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss,    &
                   ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)
    if( ((ZZg1 .lt. ES%Z_xpoint(1)) .and. (xcase .eq. LOWER_XPOINT)) .or. ((ZZg1 .gt. ES%Z_xpoint(2)) .and. (xcase .eq. UPPER_XPOINT)) ) then
      call interp_RZ(node_list,element_list,i_elm_find(2),s_find(2),t_find(2),&
                     RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss,  &
                     ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)
    endif

    nwpts%R_sep(j) = RRg1
    nwpts%Z_sep(j) = ZZg1
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

        call interp_RZ(node_list,element_list,i_elm_find(1),s_find(1),t_find(1),&
                       RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss,    &
                       ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)
        if(ZZg1 .lt. ES%Z_xpoint(1)) then
          call interp_RZ(node_list,element_list,i_elm_find(2),s_find(2),t_find(2),&
                         RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss,  &
                         ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)
        endif
        
        nwpts%R_sep(j) = RRg1
        nwpts%Z_sep(j) = ZZg1
      else
        if ( (ES%active_xpoint .eq. LOWER_XPOINT) .or. (ES%active_xpoint .eq. SYMMETRIC_XPOINT) ) then
          call find_theta_surface(node_list,element_list,flux_list,i_sep2,theta_sep(j),ES%R_axis,ES%Z_axis,i_elm_find,s_find,t_find,i_find)
        else
          call find_theta_surface(node_list,element_list,flux_list,i_sep1,theta_sep(j),ES%R_axis,ES%Z_axis,i_elm_find,s_find,t_find,i_find)
        endif
        if(i_find .eq. 0) return

        call interp_RZ(node_list,element_list,i_elm_find(1),s_find(1),t_find(1),&
                       RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss,    &
                       ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)
        if(ZZg1 .gt. ES%Z_xpoint(2)) then
          call interp_RZ(node_list,element_list,i_elm_find(2),s_find(2),t_find(2),&
                         RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss,  &
                         ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)
        endif

        nwpts%R_sep(j) = RRg1
        nwpts%Z_sep(j) = ZZg1    
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
if (xcase .eq. DOUBLE_NULL ) then
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
  
    if ( (ES%active_xpoint .eq. LOWER_XPOINT) .or. (ES%active_xpoint .eq. SYMMETRIC_XPOINT) ) then
      if (j .gt. n_tht_mid) then
        nwpts%Z_max(j) = nwpts%Z_sep(j) + (stpts%ZLimit_LowerInnerLeg - ES%Z_xpoint(1)) * ((nwpts%Z_sep(j) - ES%Z_axis)/(ES%Z_xpoint(1) - ES%Z_axis))**2
        i_max = n_flux + n_open + n_outer + n_inner
      else
        nwpts%Z_max(j) = nwpts%Z_sep(j) + (stpts%ZLimit_LowerOuterLeg - ES%Z_xpoint(1)) * ((nwpts%Z_sep(j) - ES%Z_axis)/(ES%Z_xpoint(1) - ES%Z_axis))**2
        i_max = n_flux + n_open + n_outer
      endif      
      call find_Z_surface(node_list,element_list,flux_list,i_max,nwpts%Z_max(j),i_elm_find,s_find,t_find,st_find,i_find)    
    elseif ( (xcase .eq. DOUBLE_NULL) .and. ( ES%active_xpoint .eq. UPPER_XPOINT ) ) then
      if (j .gt. n_tht_mid) then
        nwpts%Z_max(j) = nwpts%Z_sep(j) + (stpts%ZLimit_LowerOuterLeg - ES%Z_xpoint(1)) * ((nwpts%Z_sep(j) - ES%Z_axis)/(ES%Z_xpoint(1) - ES%Z_axis))**2
        i_max = n_flux + n_open + n_outer
      else
        nwpts%Z_max(j) = nwpts%Z_sep(j) + (stpts%ZLimit_LowerInnerLeg - ES%Z_xpoint(1)) * ((nwpts%Z_sep(j) - ES%Z_axis)/(ES%Z_xpoint(1) - ES%Z_axis))**2
        i_max = n_flux + n_open + n_outer + n_inner
      endif      
      call find_Z_surface(node_list,element_list,flux_list,i_max,nwpts%Z_max(j),i_elm_find,s_find,t_find,st_find,i_find)    
    else    
      i_max = n_flux + n_open
      call find_theta_surface(node_list,element_list,flux_list,i_max,theta_sep(j),ES%R_axis,ES%Z_axis,i_elm_find,s_find,t_find,i_find)    
    endif

    call interp_RZ(node_list,element_list,i_elm_find(1),s_find(1),t_find(1), &
                   RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss, &
                   ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)

    if(    (xcase .eq. UPPER_XPOINT)  &
      .or. (     ( (ES%active_xpoint .eq. LOWER_XPOINT) .or. (ES%active_xpoint .eq. SYMMETRIC_XPOINT) ) &
           .and. (    ( (RRg1 .gt. ES%R_xpoint(1)) .and. (j.lt.n_tht_mid) ) &
                 .or. ( (RRg1 .lt. ES%R_xpoint(1)) .and. (j.gt.n_tht_mid) ) ) ) &
      .or. (     ( (xcase .eq. DOUBLE_NULL) .and. ( ES%active_xpoint .eq. UPPER_XPOINT ) ) &
           .and. (    ( (RRg1 .lt. ES%R_xpoint(1)) .and. (j.le.n_tht_mid) ) &
                 .or. ( (RRg1 .gt. ES%R_xpoint(1)) .and. (j.gt.n_tht_mid) ) ) ) ) then

      nwpts%R_max(j)             = RRg1
      nwpts%Z_max(j)             = ZZg1
      nwpts%RR_new(i_max+1,j)    = RRg1
      nwpts%ZZ_new(i_max+1,j)    = ZZg1
      nwpts%ielm_flux(i_max+1,j) = i_elm_find(1)
      nwpts%s_flux(i_max+1,j)    = s_find(1)
      nwpts%t_flux(i_max+1,j)    = t_find(1)
      nwpts%t_tht(i_max+1,j)     = 1.d0
      nwpts%k_cross(i_max+1,j)   = 3
      
    else

      call interp_RZ(node_list,element_list,i_elm_find(2),s_find(2),t_find(2),&
                     RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss,    &
                     ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)
      nwpts%R_max(j)             = RRg1
      nwpts%Z_max(j)             = ZZg1
      nwpts%RR_new(i_max+1,j)    = RRg1
      nwpts%ZZ_new(i_max+1,j)    = ZZg1
      nwpts%ielm_flux(i_max+1,j) = i_elm_find(2)
      nwpts%s_flux(i_max+1,j)    = s_find(2)
      nwpts%t_flux(i_max+1,j)    = t_find(2)
      nwpts%t_tht(i_max+1,j)     = 1.d0
      nwpts%k_cross(i_max+1,j)   = 3

    endif
    
    if (     ( (ES%active_xpoint .eq. LOWER_XPOINT) .or. (ES%active_xpoint .eq. SYMMETRIC_XPOINT) ) &
       .and. ((j .eq. 1) .or. (j .eq. n_tht)) ) then
      nwpts%R_max(1)              = stpts%RLimit_LowerOuterLeg ! this one is known - safer...
      nwpts%Z_max(1)              = stpts%ZLimit_LowerOuterLeg ! this one is known - safer...
      nwpts%RR_new(i_max+1,1)     = stpts%RLimit_LowerOuterLeg ! this one is known - safer...
      nwpts%ZZ_new(i_max+1,1)     = stpts%ZLimit_LowerOuterLeg ! this one is known - safer...
      nwpts%R_max(n_tht)          = stpts%RLimit_LowerInnerLeg ! this one is known - safer...
      nwpts%Z_max(n_tht)          = stpts%ZLimit_LowerInnerLeg ! this one is known - safer...
      nwpts%RR_new(i_max+1,n_tht) = stpts%RLimit_LowerInnerLeg ! this one is known - safer...
      nwpts%ZZ_new(i_max+1,n_tht) = stpts%ZLimit_LowerInnerLeg ! this one is known - safer...
    endif
    if (     ((xcase .eq. DOUBLE_NULL) .and. ( ES%active_xpoint .eq. UPPER_XPOINT ) )  &
       .and. ((j .eq. 1) .or. (j .eq. n_tht)) ) then
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
        
    if (    ( (j .gt. n_tht_mid) .and. (xcase .eq. DOUBLE_NULL ) .and. ( (ES%active_xpoint .eq. LOWER_XPOINT) .or. (ES%active_xpoint .eq. SYMMETRIC_XPOINT) ) ) & 
       .or. ( (j .lt. n_tht_mid) .and. (xcase .eq. DOUBLE_NULL ) .and. (  ES%active_xpoint .eq. UPPER_XPOINT)                                          ) & 
       .or. ( (j .lt. n_tht_mid) .and. (xcase .eq. UPPER_XPOINT) ) ) then
      nwpts%Z_max(j) = nwpts%Z_sep(j) + (stpts%ZLimit_UpperInnerLeg - ES%Z_xpoint(2)) * ((nwpts%Z_sep(j) - ES%Z_axis)/(ES%Z_xpoint(2) - ES%Z_axis))**2
      i_max = n_flux + n_open + n_outer + n_inner
      call find_Z_surface(node_list,element_list,flux_list,i_max,nwpts%Z_max(j),i_elm_find,s_find,t_find,st_find,i_find)
    endif 
    if (    ( (j .le. n_tht_mid) .and. (xcase .eq. DOUBLE_NULL ) .and. ( (ES%active_xpoint .eq. LOWER_XPOINT) .or. (ES%active_xpoint .eq. SYMMETRIC_XPOINT) ) ) &
       .or. ( (j .gt. n_tht_mid) .and. (xcase .eq. DOUBLE_NULL ) .and. (  ES%active_xpoint .eq. UPPER_XPOINT )                                         ) & 
       .or. ( (j .gt. n_tht_mid) .and. (xcase .eq. UPPER_XPOINT) ) ) then
      nwpts%Z_max(j) = nwpts%Z_sep(j) + (stpts%ZLimit_UpperOuterLeg - ES%Z_xpoint(2)) * ((nwpts%Z_sep(j) - ES%Z_axis)/(ES%Z_xpoint(2) - ES%Z_axis))**2
      i_max = n_flux + n_open + n_outer
      call find_Z_surface(node_list,element_list,flux_list,i_max,nwpts%Z_max(j),i_elm_find,s_find,t_find,st_find,i_find)
    endif
    if (xcase .eq. LOWER_XPOINT) then
      call find_theta_surface(node_list,element_list,flux_list,i_max,theta_sep(j),ES%R_axis,ES%Z_axis,i_elm_find,s_find,t_find,i_find)
    endif

    call interp_RZ(node_list,element_list,i_elm_find(1),s_find(1),t_find(1),&
                   RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss,    &
                   ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)


    if(    (xcase .eq. LOWER_XPOINT)  &
      .or. (     ( (xcase .eq. DOUBLE_NULL) .and. ( (ES%active_xpoint .eq. LOWER_XPOINT) .or. (ES%active_xpoint .eq. SYMMETRIC_XPOINT) ) )  &
           .and. (    ( (RRg1 .ge. ES%R_xpoint(2)) .and. (j .le. n_tht_mid) ) &
                 .or. ( (RRg1 .lt. ES%R_xpoint(2)) .and. (j .gt. n_tht_mid) ) ) ) &
      .or. (     (  ES%active_xpoint .eq. UPPER_XPOINT                                                                            )  &
           .and. (    ( (RRg1 .ge. ES%R_xpoint(2)) .and. (j .gt. n_tht_mid) ) &
                 .or. ( (RRg1 .lt. ES%R_xpoint(2)) .and. (j .lt. n_tht_mid) ) ) ) ) then

      nwpts%R_max(j)             = RRg1
      nwpts%Z_max(j)             = ZZg1
      nwpts%RR_new(i_max+1,j)    = RRg1
      nwpts%ZZ_new(i_max+1,j)    = ZZg1
      nwpts%ielm_flux(i_max+1,j) = i_elm_find(1)
      nwpts%s_flux(i_max+1,j)    = s_find(1)
      nwpts%t_flux(i_max+1,j)    = t_find(1)
      nwpts%t_tht(i_max+1,j)     = 1.d0
      nwpts%k_cross(i_max+1,j)   = 3
      
    else

      call interp_RZ(node_list,element_list,i_elm_find(2),s_find(2),t_find(2),&
                     RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss,    &
                     ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)
      nwpts%R_max(j)             = RRg1
      nwpts%Z_max(j)             = ZZg1
      nwpts%RR_new(i_max+1,j)    = RRg1
      nwpts%ZZ_new(i_max+1,j)    = ZZg1
      nwpts%ielm_flux(i_max+1,j) = i_elm_find(2)
      nwpts%s_flux(i_max+1,j)    = s_find(2)
      nwpts%t_flux(i_max+1,j)    = t_find(2)
      nwpts%t_tht(i_max+1,j)     = 1.d0
      nwpts%k_cross(i_max+1,j)   = 3

    endif
    
    if (     (  ES%active_xpoint .eq. UPPER_XPOINT ) & 
       .and. ( (j .eq. 1) .or. (j .eq. n_tht_2) ) ) then
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
       .and. ( (j .eq. n_tht_mid) .or. (j .eq. n_tht_mid+1) ) ) then
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






!----------------------------------- Print a python file that plots the bound points
if (plot_grid) then
  open(100,file='plot_bound_points.py')
    write(100,'(A)')         '#!/usr/bin/env python'
    write(100,'(A)')         'import numpy as N'
    write(100,'(A)')         'import pylab'
    write(100,'(A)')         'def main():'
    write(100,'(A,i6,A)')     ' r = N.zeros(',2*n_tht_2+2*n_leg+2*n_up_leg,')'
    write(100,'(A,i6,A)')     ' z = N.zeros(',2*n_tht_2+2*n_leg+2*n_up_leg,')'
    do j=1,n_tht_2
      write(100,'(A,i6,A,f15.4)') ' r[',j-1,'] = ',nwpts%R_sep(j)
      write(100,'(A,i6,A,f15.4)') ' z[',j-1,'] = ',nwpts%Z_sep(j)
    enddo
    do j=1,n_tht_2
      write(100,'(A,i6,A,f15.4)') ' r[',n_tht_2+j-1,'] = ',nwpts%R_max(j)
      write(100,'(A,i6,A,f15.4)') ' z[',n_tht_2+j-1,'] = ',nwpts%Z_max(j)
    enddo
    do j=1,2*n_leg
      write(100,'(A,i6,A,f15.4)') ' r[',2*n_tht_2+j-1,'] = ',nwpts%R_min(j)
      write(100,'(A,i6,A,f15.4)') ' z[',2*n_tht_2+j-1,'] = ',nwpts%Z_min(j)
    enddo
    do j=1,2*n_up_leg
      write(100,'(A,i6,A,f15.4)') ' r[',2*n_tht_2+2*n_leg+j-1,'] = ',nwpts%R_min(2*n_leg+j)
      write(100,'(A,i6,A,f15.4)') ' z[',2*n_tht_2+2*n_leg+j-1,'] = ',nwpts%Z_min(2*n_leg+j)
    enddo
    write(100,'(A,i6,A)')     ' for i in range (0,',2*n_tht_2+2*n_leg+2*n_up_leg,'):'
    write(100,'(A)')         '  pylab.plot(r[i:i+1],z[i:i+1], "r.")'
    write(100,'(A)')         ' pylab.axis("equal")'
    write(100,'(A)')         ' pylab.show()'
    write(100,'(A)')         ' '
    write(100,'(A)')         'main()'
  close(100)
endif

!----------------------------------- Print a python file that plots the extrapolation points
if (plot_grid) then
  open(101,file='plot_extra_points.py')
    write(101,'(A)')         '#!/usr/bin/env python'
    write(101,'(A)')         'import numpy as N'
    write(101,'(A)')         'import pylab'
    write(101,'(A)')         'def main():'
    write(101,'(A,i6,A)')     ' r = N.zeros(',3*n_tht_2,')'
    write(101,'(A,i6,A)')     ' z = N.zeros(',3*n_tht_2,')'
    do j=1,n_tht
      write(101,'(A,i6,A,f15.4)') ' r[',3*(j-1),'] = ',ES%R_axis
      write(101,'(A,i6,A,f15.4)') ' z[',3*(j-1),'] = ',ES%Z_axis
      write(101,'(A,i6,A,f15.4)') ' r[',3*(j-1)+1,'] = ',nwpts%R_sep(j)
      write(101,'(A,i6,A,f15.4)') ' z[',3*(j-1)+1,'] = ',nwpts%Z_sep(j)
      write(101,'(A,i6,A,f15.4)') ' r[',3*(j-1)+2,'] = ',nwpts%RR_new(n_psi,j)
      write(101,'(A,i6,A,f15.4)') ' z[',3*(j-1)+2,'] = ',nwpts%ZZ_new(n_psi,j)
    enddo
    do j=n_tht+1,n_tht_2
      write(101,'(A,i6,A,f15.4)') ' r[',3*(j-1),'] = ',nwpts%R_min(j)
      write(101,'(A,i6,A,f15.4)') ' z[',3*(j-1),'] = ',nwpts%Z_min(j)
      write(101,'(A,i6,A,f15.4)') ' r[',3*(j-1)+1,'] = ',nwpts%R_sep(j)
      write(101,'(A,i6,A,f15.4)') ' z[',3*(j-1)+1,'] = ',nwpts%Z_sep(j)
      write(101,'(A,i6,A,f15.4)') ' r[',3*(j-1)+2,'] = ',nwpts%R_max(j)
      write(101,'(A,i6,A,f15.4)') ' z[',3*(j-1)+2,'] = ',nwpts%Z_max(j)
    enddo
    write(101,'(A,i6,A)')     ' for i in range (0,',n_tht_2,'):'
    write(101,'(A)')         '  pylab.plot(r[3*i:3*i+3],z[3*i:3*i+3], "r")'
    write(101,'(A)')         '  pylab.plot(r[3*i:3*i+3],z[3*i:3*i+3], "r+")'
    write(101,'(A)')         ' pylab.axis("equal")'
    write(101,'(A)')         ' pylab.show()'
    write(101,'(A)')         ' '
    write(101,'(A)')         'main()'
  close(101)
endif


return
end subroutine find_extrapolation_points_central_part
