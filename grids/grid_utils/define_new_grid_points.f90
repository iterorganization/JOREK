!> Subroutine defines the new grid_points from crossing of polar and radial coordinate lines
subroutine define_new_grid_points(node_list, element_list, flux_list, xcase, n_grids, stpts, sigmas, nwpts)

use constants
use tr_module 
use data_structure
use grid_xpoint_data
use phys_module, only: tokamak_device, write_ps
use mod_interp
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
  if ( (ES%active_xpoint .eq. LOWER_XPOINT) .or. (ES%active_xpoint .eq. SYMMETRIC_XPOINT) ) then
    n_tht_mid = int(n_tht * (2.d0*PI - (tht_x1 - tht_x2)) / (2.d0*PI))
    ! Make sure n_tht_mid is odd and save it to n_grids for later use
    if(mod(n_tht_mid,2) .eq. 0) n_tht_mid = n_tht_mid + 1
    n_grids(10) = n_tht_mid
    
    call tr_allocate(s_tmp,1,n_tht_mid,"s_tmp",CAT_GRID)
    s_tmp = 0
    call meshac2(n_tht_mid,s_tmp,0.d0,1.d0,SIG_theta,SIG_theta_up,bgf_tht,1.0d0)

    do j=1,n_tht_mid
      theta_sep(j) = (tht_x1-2.d0*PI) + (tht_x2-(tht_x1-2.d0*PI)) * s_tmp(j)
      if (theta_sep(j) .lt. 0.d0)    theta_sep(j) = theta_sep(j) + 2.d0 * PI
      if (theta_sep(j) .gt. 2.d0*PI) theta_sep(j) = theta_sep(j) - 2.d0 * PI
    enddo
    call tr_deallocate(s_tmp,"s_tmp",CAT_GRID)
  
    n_tht_mid2 = n_tht-n_tht_mid
    call tr_allocate(s_tmp,1,n_tht_mid2,"s_tmp",CAT_GRID)
    s_tmp = 0
    call meshac2(n_tht_mid2,s_tmp,0.d0,1.d0,SIG_theta_up,SIG_theta,bgf_tht,1.0d0)

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
    call meshac2(n_tht_mid,s_tmp,0.d0,1.d0,SIG_theta_up,SIG_theta,bgf_tht,1.0d0)
    do j=1,n_tht_mid
      theta_sep(j) = tht_x1 + (tht_x2-tht_x1) * s_tmp(j)
      if (theta_sep(j) .lt. 0.d0)    theta_sep(j) = theta_sep(j) + 2.d0 * PI
      if (theta_sep(j) .gt. 2.d0*PI) theta_sep(j) = theta_sep(j) - 2.d0 * PI
    enddo
    call tr_deallocate(s_tmp,"s_tmp",CAT_GRID)

    n_tht_mid2 = n_tht-n_tht_mid
    call tr_allocate(s_tmp,1,n_tht_mid2,"s_tmp",CAT_GRID)
    s_tmp = 0
    call meshac2(n_tht_mid2,s_tmp,0.d0,1.d0,SIG_theta,SIG_theta_up,bgf_tht,1.0d0)
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
    elseif ( (xcase .eq. DOUBLE_NULL) .and. (ES%active_xpoint .eq. UPPER_XPOINT) ) then
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

    ! --- Readjust s,t (sometimes outside element after find_Z_surface and find_R_surface...)
    if (s_find(1) .lt. 0.d0) s_find(1) = 0.d0
    if (s_find(1) .gt. 1.d0) s_find(1) = 1.d0
    if (t_find(1) .lt. 0.d0) t_find(1) = 0.d0
    if (t_find(1) .gt. 1.d0) t_find(1) = 1.d0
      
    call interp_RZ(node_list,element_list,i_elm_find(1),s_find(1),t_find(1), &
                   RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss, &
                   ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)

    if(    (xcase .eq. UPPER_XPOINT)  &
      .or. (     ( (ES%active_xpoint .eq. LOWER_XPOINT) .or. (ES%active_xpoint .eq. SYMMETRIC_XPOINT) ) &
           .and. (    ( (RRg1 .gt. ES%R_xpoint(1)) .and. (j.lt.n_tht_mid) ) &
                 .or. ( (RRg1 .lt. ES%R_xpoint(1)) .and. (j.gt.n_tht_mid) ) ) ) &
      .or. (     ( (xcase .eq. DOUBLE_NULL) .and. (ES%active_xpoint .eq. UPPER_XPOINT) ) &
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

      ! --- Readjust s,t (sometimes outside element after find_Z_surface and find_R_surface...)
      if (s_find(2) .lt. 0.d0) s_find(2) = 0.d0
      if (s_find(2) .gt. 1.d0) s_find(2) = 1.d0
      if (t_find(2) .lt. 0.d0) t_find(2) = 0.d0
      if (t_find(2) .gt. 1.d0) t_find(2) = 1.d0
      
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
    if (     ((xcase .eq. DOUBLE_NULL) .and. ( ES%active_xpoint .eq. UPPER_XPOINT ) )  &
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
        
    if (    ( (j .gt. n_tht_mid) .and. (xcase .eq. DOUBLE_NULL) .and. ( (ES%active_xpoint .eq. LOWER_XPOINT) .or. (ES%active_xpoint .eq. SYMMETRIC_XPOINT) ) ) & 
       .or. ( (j .lt. n_tht_mid) .and. (xcase .eq. DOUBLE_NULL) .and. (  ES%active_xpoint .eq. UPPER_XPOINT                                         ) ) & 
       .or. ( (j .lt. n_tht_mid) .and. (xcase .eq. UPPER_XPOINT)                                                                                    ) ) then
      nwpts%Z_max(j) = nwpts%Z_sep(j) + (stpts%ZLimit_UpperInnerLeg - ES%Z_xpoint(2)) * ((nwpts%Z_sep(j) - ES%Z_axis)/(ES%Z_xpoint(2) - ES%Z_axis))**2
      i_max = n_flux + n_open + n_outer + n_inner
      call find_Z_surface(node_list,element_list,flux_list,i_max,nwpts%Z_max(j),i_elm_find,s_find,t_find,st_find,i_find)
    endif 
    if (    ( (j .le. n_tht_mid) .and. (xcase .eq. DOUBLE_NULL) .and. ( (ES%active_xpoint .eq. LOWER_XPOINT) .or. (ES%active_xpoint .eq. SYMMETRIC_XPOINT) ) ) &
       .or. ( (j .gt. n_tht_mid) .and. (xcase .eq. DOUBLE_NULL) .and. (  ES%active_xpoint .eq. UPPER_XPOINT                                         ) ) & 
       .or. ( (j .gt. n_tht_mid) .and. (xcase .eq. UPPER_XPOINT)                                                                                    ) ) then
      nwpts%Z_max(j) = nwpts%Z_sep(j) + (stpts%ZLimit_UpperOuterLeg - ES%Z_xpoint(2)) * ((nwpts%Z_sep(j) - ES%Z_axis)/(ES%Z_xpoint(2) - ES%Z_axis))**2
      i_max = n_flux + n_open + n_outer
      call find_Z_surface(node_list,element_list,flux_list,i_max,nwpts%Z_max(j),i_elm_find,s_find,t_find,st_find,i_find)
    endif
    if (xcase .eq. LOWER_XPOINT) then
      call find_theta_surface(node_list,element_list,flux_list,i_max,theta_sep(j),ES%R_axis,ES%Z_axis,i_elm_find,s_find,t_find,i_find)
    endif

    ! --- Readjust s,t (sometimes outside element after find_Z_surface and find_R_surface...)
    if (s_find(1) .lt. 0.d0) s_find(1) = 0.d0
    if (s_find(1) .gt. 1.d0) s_find(1) = 1.d0
    if (t_find(1) .lt. 0.d0) t_find(1) = 0.d0
    if (t_find(1) .gt. 1.d0) t_find(1) = 1.d0
      
    call interp_RZ(node_list,element_list,i_elm_find(1),s_find(1),t_find(1),&
                   RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss,    &
                   ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)


    if(    (xcase .eq. LOWER_XPOINT) &
      .or. (     ( (xcase .eq. DOUBLE_NULL) .and. ( (ES%active_xpoint .eq. LOWER_XPOINT) .or. (ES%active_xpoint .eq. SYMMETRIC_XPOINT) ) ) &
           .and. (    ( (RRg1 .ge. ES%R_xpoint(2)) .and. (j .le. n_tht_mid) ) &
                 .or. ( (RRg1 .lt. ES%R_xpoint(2)) .and. (j .gt. n_tht_mid) ) ) ) &
      .or. (     ( ES%active_xpoint .eq. UPPER_XPOINT                                                                             ) &
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

      ! --- Readjust s,t (sometimes outside element after find_Z_surface and find_R_surface...)
      if (s_find(2) .lt. 0.d0) s_find(2) = 0.d0
      if (s_find(2) .gt. 1.d0) s_find(2) = 1.d0
      if (t_find(2) .lt. 0.d0) t_find(2) = 0.d0
      if (t_find(2) .gt. 1.d0) t_find(2) = 1.d0
      
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
    
    if (     ( ES%active_xpoint .eq. UPPER_XPOINT ) & 
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




!--------------------------------------------------------------------------!
!------- Extrapolation points for second part of the grid (Legs) ----------!
!--------------------------------------------------------------------------!

!------------------------------ Intersections with private surfaces
do i=1,4
  
  if (xcase .ne. UPPER_XPOINT) then 
    if (i .eq. 1) then
      n_loop  = n_leg
      n_start = n_tht
      i_surf  = n_flux+n_open+n_outer+n_inner+n_private
      R_beg   = stpts%RRightCorn_LowerInnerLeg
      R_end   = stpts%RMiddle_LowerPrivate
      SIG_0   = SIG_leg_0
      SIG_1   = SIG_leg_1
    endif
    if (i .eq. 2) then
      n_loop  = n_leg
      n_start = n_tht + n_leg
      i_surf  = n_flux+n_open+n_outer+n_inner+n_private
      R_beg   = stpts%RLeftCorn_LowerOuterLeg
      R_end   = stpts%RMiddle_LowerPrivate
      SIG_0   = SIG_leg_0
      SIG_1   = SIG_leg_1
    endif
    if (i .eq. 3) then
      if(xcase .eq. DOUBLE_NULL) then
        n_loop  = n_up_leg
        n_start = n_tht + n_leg + n_leg
        i_surf  = n_flux+n_open+n_outer+n_inner+n_private+n_up_priv
        R_beg   = stpts%RRightCorn_UpperInnerLeg
        R_end   = stpts%RMiddle_UpperPrivate
        SIG_0   = SIG_up_leg_0
        SIG_1   = SIG_up_leg_1
      else
        exit
      endif
    endif
    if (i .eq. 4) then
      if(xcase .eq. DOUBLE_NULL) then
        n_loop  = n_up_leg
        n_start = n_tht + n_leg + n_leg + n_up_leg
        i_surf  = n_flux+n_open+n_outer+n_inner+n_private+n_up_priv
        R_beg   = stpts%RLeftCorn_UpperOuterLeg
        R_end   = stpts%RMiddle_UpperPrivate
        SIG_0   = SIG_up_leg_0
        SIG_1   = SIG_up_leg_1
      else
        exit
      endif
    endif
  else
    if (i .eq. 1) then
      n_loop  = n_up_leg
      n_start = n_tht
      i_surf  = n_flux+n_open+n_up_priv
      R_beg   = stpts%RRightCorn_UpperInnerLeg
      R_end   = stpts%RMiddle_UpperPrivate
      SIG_0   = SIG_up_leg_0
      SIG_1   = SIG_up_leg_1
    endif
    if (i .eq. 2) then
      n_loop  = n_up_leg
      n_start = n_tht + n_up_leg
      i_surf  = n_flux+n_open+n_up_priv
      R_beg   = stpts%RLeftCorn_UpperOuterLeg
      R_end   = stpts%RMiddle_UpperPrivate
      SIG_0   = SIG_up_leg_0
      SIG_1   = SIG_up_leg_1
    endif
    if (i .eq. 3) exit
    if (i .eq. 4) exit
  endif
  
  call tr_allocate(s_tmp,1,n_loop,"s_tmp",CAT_GRID)
  s_tmp = 0
  call meshac2(n_loop,s_tmp,0.d0,1.d0,SIG_0,SIG_1,0.6d0,1.0d0)
  do j=1,n_loop

    nwpts%R_min(n_start + j) = R_beg + (R_end-R_beg) * s_tmp(j)

    call find_R_surface(node_list,element_list,flux_list,i_surf,nwpts%R_min(n_start+j),i_elm_find,s_find,t_find,st_find,i_find)

    do k=1,i_find
 
      ! --- Readjust s,t (sometimes outside element after find_Z_surface and find_R_surface...)
      if (s_find(k) .lt. 0.d0) s_find(k) = 0.d0
      if (s_find(k) .gt. 1.d0) s_find(k) = 1.d0
      if (t_find(k) .lt. 0.d0) t_find(k) = 0.d0
      if (t_find(k) .gt. 1.d0) t_find(k) = 1.d0
      
      call interp_RZ(node_list,element_list,i_elm_find(k),s_find(k),t_find(k),&
                     RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss,    &
                     ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)

      if ( (xcase .eq. LOWER_XPOINT) .and. (ZZg1 .le. ES%Z_xpoint(1)) ) exit
      if ( (xcase .eq. UPPER_XPOINT) .and. (ZZg1 .ge. ES%Z_xpoint(2)) ) exit
      if ( (xcase .eq. DOUBLE_NULL ) .and. (ZZg1 .le. ES%Z_xpoint(1)) .and. (i .le. 2) ) exit
      if ( (xcase .eq. DOUBLE_NULL ) .and. (ZZg1 .ge. ES%Z_xpoint(2)) .and. (i .ge. 3) ) exit

    enddo

    nwpts%Z_min(n_start + j) = ZZg1

  enddo
  call tr_deallocate(s_tmp,"s_tmp",CAT_GRID)
  if ( (xcase .ne. UPPER_XPOINT) .and. (i .eq. 1) ) nwpts%Z_min(n_start+1) = stpts%ZRightCorn_LowerInnerLeg ! this one is known - safer...
  if ( (xcase .ne. UPPER_XPOINT) .and. (i .eq. 2) ) nwpts%Z_min(n_start+1) = stpts%ZLeftCorn_LowerOuterLeg  ! this one is known - safer...
  if ( (xcase .eq. DOUBLE_NULL ) .and. (i .eq. 3) ) nwpts%Z_min(n_start+1) = stpts%ZRightCorn_UpperInnerLeg ! this one is known - safer...
  if ( (xcase .eq. DOUBLE_NULL ) .and. (i .eq. 4) ) nwpts%Z_min(n_start+1) = stpts%ZLeftCorn_UpperOuterLeg  ! this one is known - safer...
  if ( (xcase .eq. UPPER_XPOINT) .and. (i .eq. 1) ) nwpts%Z_min(n_start+1) = stpts%ZRightCorn_UpperInnerLeg ! this one is known - safer...
  if ( (xcase .eq. UPPER_XPOINT) .and. (i .eq. 2) ) nwpts%Z_min(n_start+1) = stpts%ZLeftCorn_UpperOuterLeg  ! this one is known - safer...
enddo

!------------------------------ Intersection with last surface on wall for MAST
if(tokamak_device(1:4) .eq. 'MAST') then
  do i=1,2
    
    if (i .eq. 1) then
      n_loop  = n_leg
      n_start = 0
      i_surf  = stpts%i_surf_wall_low
      Z_beg   = stpts%ZLimit_LowerMastWall
      Z_end   = stpts%ZLimit_LowerMastWallBox
      SIG_0   = SIG_leg_0
      SIG_1   = SIG_leg_1
    endif
    if (i .eq. 2) then
      n_loop  = n_up_leg
      n_start = n_leg
      i_surf  = stpts%i_surf_wall_up
      Z_beg   = stpts%ZLimit_UpperMastWall
      Z_end   = stpts%ZLimit_UpperMastWallBox
      SIG_0   = SIG_up_leg_0
      SIG_1   = SIG_up_leg_1
    endif
    
    call tr_allocate(s_tmp,1,n_loop,"s_tmp",CAT_GRID)
    s_tmp = 0
    call meshac2(n_loop,s_tmp,0.d0,1.d0,SIG_0,SIG_1,0.6d0,1.0d0)
    do j=1,n_loop

      nwpts%Z_wall(n_start + j) = Z_beg + (Z_end-Z_beg) * s_tmp(j)

      call find_Z_surface(node_list,element_list,flux_list,i_surf,nwpts%Z_wall(n_start + j),i_elm_find,s_find,t_find,st_find,i_find)

      do k=1,i_find

        call interp_RZ(node_list,element_list,i_elm_find(k),s_find(k),t_find(k),&
                       RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss,    &
                       ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)

        if ( (i .eq. 1) .and. (RRg1 .gt. ES%R_xpoint(1)) ) exit
        if ( (i .eq. 2) .and. (RRg1 .gt. ES%R_xpoint(2)) ) exit

      enddo

      nwpts%R_wall(n_start + j) = RRg1

    enddo
    call tr_deallocate(s_tmp,"s_tmp",CAT_GRID)
  enddo
endif

!------------------------------ Intersections with separatrices
do i=1,4
  
  if (xcase .ne. UPPER_XPOINT) then 
    if (i .eq. 1) then
      n_loop  = n_leg
      n_start = n_tht
      i_surf  = n_flux+n_open+n_outer+n_inner
      Z_beg   = stpts%ZLeftCorn_LowerInnerLeg
      Z_end   = stpts%ZLimit_LowerInnerLeg
      SIG_0   = SIG_leg_0
      SIG_1   = SIG_leg_1
    endif
    if (i .eq. 2) then
      n_loop  = n_leg
      n_start = n_tht + n_leg
      i_surf  = n_flux+n_open+n_outer
      if(tokamak_device(1:4) .eq. 'MAST') then
        R_beg   = stpts%RRightCorn_LowerOuterLeg
        R_end   = stpts%RLimit_LowerOuterLeg
      else
        Z_beg   = stpts%ZRightCorn_LowerOuterLeg
        Z_end   = stpts%ZLimit_LowerOuterLeg
      endif
      SIG_0   = SIG_leg_0
      SIG_1   = SIG_leg_1
    endif
    if (i .eq. 3) then
      if(xcase .eq. DOUBLE_NULL) then
        n_loop  = n_up_leg
        n_start = n_tht + n_leg + n_leg
        i_surf  = n_flux+n_open+n_outer+n_inner
        Z_beg   = stpts%ZLeftCorn_UpperInnerLeg
        Z_end   = stpts%ZLimit_UpperInnerLeg
        SIG_0   = SIG_up_leg_0
        SIG_1   = SIG_up_leg_1
      else
        exit
      endif
    endif
    if (i .eq. 4) then
      if(xcase .eq. DOUBLE_NULL) then
        n_loop  = n_up_leg
        n_start = n_tht + n_leg + n_leg + n_up_leg
        i_surf  = n_flux+n_open+n_outer
        if(tokamak_device(1:4) .eq. 'MAST') then
          R_beg   = stpts%RRightCorn_UpperOuterLeg
          R_end   = stpts%RLimit_UpperOuterLeg
        else
          Z_beg   = stpts%ZRightCorn_UpperOuterLeg
          Z_end   = stpts%ZLimit_UpperOuterLeg
        endif
        SIG_0   = SIG_up_leg_0
        SIG_1   = SIG_up_leg_1
      else
        exit
      endif
    endif
  else
    if (i .eq. 1) then
      n_loop  = n_up_leg
      n_start = n_tht
      i_surf  = n_flux+n_open
      Z_beg   = stpts%ZLeftCorn_UpperInnerLeg
      Z_end   = stpts%ZLimit_UpperInnerLeg
      SIG_0   = SIG_up_leg_0
      SIG_1   = SIG_up_leg_1
    endif
    if (i .eq. 2) then
      n_loop  = n_up_leg
      n_start = n_tht + n_up_leg
      i_surf  = n_flux+n_open
      Z_beg   = stpts%ZRightCorn_UpperOuterLeg
      Z_end   = stpts%ZLimit_UpperOuterLeg
      SIG_0   = SIG_up_leg_0
      SIG_1   = SIG_up_leg_1
    endif
    if (i .eq. 3) exit
    if (i .eq. 4) exit
  endif
  
  call tr_allocate(s_tmp,1,n_loop,"s_tmp",CAT_GRID)
  s_tmp = 0
  call meshac2(n_loop,s_tmp,0.d0,1.d0,SIG_0,SIG_1,0.6d0,1.0d0)
  do j=1,n_loop

    if ((i .eq. 1) .or. (i .eq. 3) .or. (tokamak_device(1:4) .ne. 'MAST')) then
      nwpts%Z_max(n_start + j) = Z_beg + (Z_end-Z_beg) * s_tmp(j)
      call find_Z_surface(node_list,element_list,flux_list,i_surf,nwpts%Z_max(n_start+j),i_elm_find,s_find,t_find,st_find,i_find)
    else
      nwpts%R_max(n_start + j) = R_beg + (R_end-R_beg) * s_tmp(j)
      call find_R_surface(node_list,element_list,flux_list,i_surf,nwpts%R_max(n_start+j),i_elm_find,s_find,t_find,st_find,i_find)
    endif


    do k=1,i_find

      ! --- Readjust s,t (sometimes outside element after find_Z_surface and find_R_surface...)
      if (s_find(k) .lt. 0.d0) s_find(k) = 0.d0
      if (s_find(k) .gt. 1.d0) s_find(k) = 1.d0
      if (t_find(k) .lt. 0.d0) t_find(k) = 0.d0
      if (t_find(k) .gt. 1.d0) t_find(k) = 1.d0
      
      call interp_RZ(node_list,element_list,i_elm_find(k),s_find(k),t_find(k),&
                     RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss,    &
                     ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)

      if   ( (xcase .eq. UPPER_XPOINT) .and. (i .eq. 1) .and. (RRg1 .le. ES%R_xpoint(2)) ) exit
      if   ( (xcase .eq. UPPER_XPOINT) .and. (i .eq. 2) .and. (RRg1 .ge. ES%R_xpoint(2)) ) exit
      if   ( (xcase .ne. UPPER_XPOINT) .and. (i .eq. 1) .and. (RRg1 .le. ES%R_xpoint(1)) ) exit
      if   ( (xcase .ne. UPPER_XPOINT) .and. (i .eq. 3) .and. (RRg1 .le. ES%R_xpoint(2)) ) exit
      if (tokamak_device(1:4) .ne. 'MAST') then
        if ( (xcase .ne. UPPER_XPOINT) .and. (i .eq. 2) .and. (RRg1 .ge. ES%R_xpoint(1)) ) exit
        if ( (xcase .ne. UPPER_XPOINT) .and. (i .eq. 4) .and. (RRg1 .ge. ES%R_xpoint(2)) ) exit
      else
        if ( (xcase .ne. UPPER_XPOINT) .and. (i .eq. 2) .and. (ZZg1 .lt. ES%Z_xpoint(1)) ) exit
        if ( (xcase .ne. UPPER_XPOINT) .and. (i .eq. 4) .and. (ZZg1 .ge. ES%Z_xpoint(2)) ) exit
      endif

    enddo

    if ((i .eq. 1) .or. (i .eq. 3) .or. (tokamak_device(1:4) .ne. 'MAST')) then
      nwpts%R_max(n_start + j) = RRg1
    else
      nwpts%Z_max(n_start + j) = ZZg1
    endif

  enddo
  call tr_deallocate(s_tmp,"s_tmp",CAT_GRID)
  if ( (xcase .eq. UPPER_XPOINT) .and. (i .eq. 1) ) nwpts%R_max(n_start+n_loop) = nwpts%R_max(n_tht)       ! this one is known - safer...
  if ( (xcase .eq. UPPER_XPOINT) .and. (i .eq. 2) ) nwpts%R_max(n_start+n_loop) = nwpts%R_max(1)           ! this one is known - safer...
  if ( (xcase .eq. DOUBLE_NULL ) .and. (i .eq. 3) ) nwpts%R_max(n_start+n_loop) = nwpts%R_max(n_tht_mid+1) ! this one is known - safer...
  if ( (xcase .eq. DOUBLE_NULL ) .and. (i .eq. 4) ) nwpts%R_max(n_start+n_loop) = nwpts%R_max(n_tht_mid)   ! this one is known - safer...
  if ( (xcase .eq. DOUBLE_NULL ) .and. ( ES%active_xpoint .eq. UPPER_XPOINT ) .and. (i .eq. 1) ) &
                                                    nwpts%R_max(n_start+n_loop) = nwpts%R_max(n_tht_mid)   ! this one is known - safer...
  if ( (xcase .eq. DOUBLE_NULL ) .and. ( ES%active_xpoint .eq. UPPER_XPOINT ) .and. (i .eq. 2) ) &
                                                    nwpts%R_max(n_start+n_loop) = nwpts%R_max(n_tht_mid+1) ! this one is known - safer...
  if ( (xcase .eq. DOUBLE_NULL ) .and. ( ES%active_xpoint .eq. UPPER_XPOINT ) .and. (i .eq. 3) ) &
                                                    nwpts%R_max(n_start+n_loop) = nwpts%R_max(1)           ! this one is known - safer...
  if ( (xcase .eq. DOUBLE_NULL ) .and. ( ES%active_xpoint .eq. UPPER_XPOINT ) .and. (i .eq. 4) ) &
                                                    nwpts%R_max(n_start+n_loop) = nwpts%R_max(n_tht)       ! this one is known - safer...
  if ( (xcase .eq. UPPER_XPOINT) .and. (i .eq. 1) ) nwpts%R_max(n_start+n_loop) = nwpts%R_max(1)           ! this one is known - safer...
  if ( (xcase .eq. UPPER_XPOINT) .and. (i .eq. 2) ) nwpts%R_max(n_start+n_loop) = nwpts%R_max(n_tht)       ! this one is known - safer...
enddo

!------------------------------ Intersections with open surfaces
do i=1,4
  
  if (xcase .ne. UPPER_XPOINT) then 
    if (i .eq. 1) then
      n_loop  = n_leg
      n_start = n_tht
      i_surf  = n_flux
      if((xcase .eq. DOUBLE_NULL) .and. ( ES%active_xpoint .eq. UPPER_XPOINT )) i_surf  = n_flux+n_open
      R_beg   = stpts%RStrike_LowerInnerLeg
      R_end   = ES%R_xpoint(1)
      SIG_0   = SIG_leg_0
      SIG_1   = SIG_leg_1
    endif
    if (i .eq. 2) then
      n_loop  = n_leg
      n_start = n_tht + n_leg
      i_surf  = n_flux
      if((xcase .eq. DOUBLE_NULL) .and. ( ES%active_xpoint .eq. UPPER_XPOINT )) i_surf  = n_flux+n_open
      R_beg   = stpts%RStrike_LowerOuterLeg
      R_end   = ES%R_xpoint(1)
      SIG_0   = SIG_leg_0
      SIG_1   = SIG_leg_1
    endif
    if (i .eq. 3) then
      if(xcase .eq. DOUBLE_NULL) then
        n_loop  = n_up_leg
        n_start = n_tht + n_leg + n_leg
        i_surf  = n_flux+n_open
        if((xcase .eq. DOUBLE_NULL) .and. ( ES%active_xpoint .eq. UPPER_XPOINT )) i_surf  = n_flux
        R_beg   = stpts%RStrike_UpperInnerLeg
        R_end   = ES%R_xpoint(2)
        SIG_0   = SIG_up_leg_0
        SIG_1   = SIG_up_leg_1
      else
        exit
      endif
    endif
    if (i .eq. 4) then
      if(xcase .eq. DOUBLE_NULL) then
        n_loop  = n_up_leg
        n_start = n_tht + n_leg + n_leg + n_up_leg
        i_surf  = n_flux+n_open
        if((xcase .eq. DOUBLE_NULL) .and. ( ES%active_xpoint .eq. UPPER_XPOINT )) i_surf  = n_flux
        R_beg   = stpts%RStrike_UpperOuterLeg
        R_end   = ES%R_xpoint(2)
        SIG_0   = SIG_up_leg_0
        SIG_1   = SIG_up_leg_1
      else
        exit
      endif
    endif
  else
    if (i .eq. 1) then
      n_loop  = n_up_leg
      n_start = n_tht
      i_surf  = n_flux
      R_beg   = stpts%RStrike_UpperInnerLeg
      R_end   = ES%R_xpoint(2)
      SIG_0   = SIG_up_leg_0
      SIG_1   = SIG_up_leg_1
    endif
    if (i .eq. 2) then
      n_loop  = n_up_leg
      n_start = n_tht + n_up_leg
      i_surf  = n_flux
      R_beg   = stpts%RStrike_UpperOuterLeg
      R_end   = ES%R_xpoint(2)
      SIG_0   = SIG_up_leg_0
      SIG_1   = SIG_up_leg_1
    endif
    if (i .eq. 3) exit
    if (i .eq. 4) exit
  endif
  
  call tr_allocate(s_tmp,1,n_loop,"s_tmp",CAT_GRID)
  s_tmp = 0
  call meshac2(n_loop,s_tmp,0.d0,1.d0,SIG_0,SIG_1,0.6d0,1.0d0)
  do j=1,n_loop

    nwpts%R_sep(n_start + j) = R_beg + (R_end-R_beg) * s_tmp(j)

    call find_R_surface(node_list,element_list,flux_list,i_surf,nwpts%R_sep(n_start+j),i_elm_find,s_find,t_find,st_find,i_find)

    do k=1,i_find

      ! --- Readjust s,t (sometimes outside element after find_Z_surface and find_R_surface...)
      if (s_find(k) .lt. 0.d0) s_find(k) = 0.d0
      if (s_find(k) .gt. 1.d0) s_find(k) = 1.d0
      if (t_find(k) .lt. 0.d0) t_find(k) = 0.d0
      if (t_find(k) .gt. 1.d0) t_find(k) = 1.d0
      
      call interp_RZ(node_list,element_list,i_elm_find(k),s_find(k),t_find(k),&
                     RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss,    &
                     ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)
      
      if ( (xcase .eq. LOWER_XPOINT) .and. (ZZg1 .le. ES%Z_xpoint(1)) ) exit
      if ( (xcase .eq. UPPER_XPOINT) .and. (ZZg1 .ge. ES%Z_xpoint(2)) ) exit
      if ( (xcase .eq. DOUBLE_NULL ) .and. (ZZg1 .le. ES%Z_xpoint(1)) .and. (i .le. 2) ) exit
      if ( (xcase .eq. DOUBLE_NULL ) .and. (ZZg1 .ge. ES%Z_xpoint(2)) .and. (i .ge. 3) ) exit

    enddo

    nwpts%Z_sep(n_start + j) = ZZg1

  enddo
  call tr_deallocate(s_tmp,"s_tmp",CAT_GRID)
  if ( (xcase .ne. UPPER_XPOINT) .and. (i .le. 2) ) nwpts%Z_sep(n_start+n_loop) = ES%Z_xpoint(1) ! this one is known - safer...
  if ( (xcase .eq. DOUBLE_NULL ) .and. (i .ge. 3) ) nwpts%Z_sep(n_start+n_loop) = ES%Z_xpoint(2) ! this one is known - safer...
  if   (xcase .eq. UPPER_XPOINT)                    nwpts%Z_sep(n_start+n_loop) = ES%Z_xpoint(2) ! this one is known - safer...
enddo

!if ( write_ps ) then
!  call lincol(2)
!  call lplot6(1,1,nwpts%R_max,nwpts%Z_max,-(n_tht_3),' ')
!  call lincol(0)
!  call lplot6(1,1,nwpts%R_sep,nwpts%Z_sep,-(n_tht_3),' ')
!endif






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

!------------------------------ Central part (without legs)
do j=1,n_tht

  delta = 0.1

  if ( (j .eq. 1) .or. (j .eq. n_tht) )       delta = 0.d0
  if ( (j .eq. 2) .or. (j .eq. n_tht - 1) )   delta = 0.05d0
  if ( ((j .eq. n_tht_mid)   .or. (j .eq. n_tht_mid+1)) .and. (xcase .eq. DOUBLE_NULL) ) delta = 0.d0
  if ( ((j .eq. n_tht_mid-1) .or. (j .eq. n_tht_mid+2)) .and. (xcase .eq. DOUBLE_NULL) ) delta = 0.05d0

  nwpts%R_polar(1,1,j) = ES%R_axis
  nwpts%R_polar(1,4,j) = delta * ES%R_axis + (1.d0 - delta) * nwpts%R_sep(j)
  nwpts%R_polar(1,2,j) = ( 2.d0 * nwpts%R_polar(1,1,j)  +         nwpts%R_polar(1,4,j) ) / 3.d0
  nwpts%R_polar(1,3,j) = (        nwpts%R_polar(1,1,j)  +  2.d0 * nwpts%R_polar(1,4,j) ) / 3.d0

  nwpts%Z_polar(1,1,j) = ES%Z_axis
  nwpts%Z_polar(1,4,j) = delta * ES%Z_axis + (1.d0 - delta) * nwpts%Z_sep(j)
  nwpts%Z_polar(1,2,j) = ( 2.d0 * nwpts%Z_polar(1,1,j)  +         nwpts%Z_polar(1,4,j) ) / 3.d0
  nwpts%Z_polar(1,3,j) = (        nwpts%Z_polar(1,1,j)  +  2.d0 * nwpts%Z_polar(1,4,j) ) / 3.d0

  nwpts%R_polar(3,1,j) = nwpts%R_max(j)
  nwpts%R_polar(3,4,j) = delta * nwpts%R_max(j) + (1.d0 - delta) * nwpts%R_sep(j)
  nwpts%R_polar(3,2,j) = ( 2.d0 * nwpts%R_polar(3,1,j)  +         nwpts%R_polar(3,4,j) ) / 3.d0
  nwpts%R_polar(3,3,j) = (        nwpts%R_polar(3,1,j)  +  2.d0 * nwpts%R_polar(3,4,j) ) / 3.d0

  nwpts%Z_polar(3,1,j) = nwpts%Z_max(j)
  nwpts%Z_polar(3,4,j) = delta * nwpts%Z_max(j) + (1.d0 - delta) * nwpts%Z_sep(j)
  nwpts%Z_polar(3,2,j) = ( 2.d0 * nwpts%Z_polar(3,1,j)  +         nwpts%Z_polar(3,4,j) ) / 3.d0
  nwpts%Z_polar(3,3,j) = (        nwpts%Z_polar(3,1,j)  +  2.d0 * nwpts%Z_polar(3,4,j) ) / 3.d0

  nwpts%R_polar(2,1,j) = nwpts%R_polar(1,4,j)
  nwpts%R_polar(2,4,j) = nwpts%R_polar(3,4,j)
  nwpts%R_polar(2,2,j) = ( nwpts%R_polar(2,1,j) +  2.d0 * nwpts%R_sep(j) ) / 3.d0
  nwpts%R_polar(2,3,j) = ( nwpts%R_polar(2,4,j) +  2.d0 * nwpts%R_sep(j) ) / 3.d0

  nwpts%Z_polar(2,1,j) = nwpts%Z_polar(1,4,j)
  nwpts%Z_polar(2,4,j) = nwpts%Z_polar(3,4,j)
  nwpts%Z_polar(2,2,j) = ( nwpts%Z_polar(2,1,j) +  2.d0 * nwpts%Z_sep(j) ) / 3.d0
  nwpts%Z_polar(2,3,j) = ( nwpts%Z_polar(2,4,j) +  2.d0 * nwpts%Z_sep(j) ) / 3.d0

enddo

!------------------------------ Lower leg
if(xcase .ne. UPPER_XPOINT) then
  do j=1,2*n_leg

    n_start = n_tht
    delta = 0.2
    if ( (j .eq. n_leg)   .or. (j .eq. 2*n_leg) )    delta = 0.d0
    if ( (j .eq. n_leg-1) .or. (j .eq. 2*n_leg-1) )  delta = 0.05d0

    if ( (j .le. n_leg) .or. (tokamak_device(1:4) .ne. 'MAST') ) then
      nwpts%R_polar(1,1,n_start+j) = nwpts%R_min(n_start+j)
      nwpts%R_polar(1,4,n_start+j) = delta * nwpts%R_min(n_start+j) + (1.d0 - delta) * nwpts%R_sep(n_start+j)
      nwpts%R_polar(1,2,n_start+j) = ( 2.d0 * nwpts%R_polar(1,1,n_start+j)  +         nwpts%R_polar(1,4,n_start+j) ) / 3.d0
      nwpts%R_polar(1,3,n_start+j) = (        nwpts%R_polar(1,1,n_start+j)  +  2.d0 * nwpts%R_polar(1,4,n_start+j) ) / 3.d0

      nwpts%Z_polar(1,1,n_start+j) = nwpts%Z_min(n_start+j)
      nwpts%Z_polar(1,4,n_start+j) = delta * nwpts%Z_min(n_start+j) + (1.d0 - delta) * nwpts%Z_sep(n_start+j)
      nwpts%Z_polar(1,2,n_start+j) = ( 2.d0 * nwpts%Z_polar(1,1,n_start+j)  +         nwpts%Z_polar(1,4,n_start+j) ) / 3.d0
      nwpts%Z_polar(1,3,n_start+j) = (        nwpts%Z_polar(1,1,n_start+j)  +  2.d0 * nwpts%Z_polar(1,4,n_start+j) ) / 3.d0

      nwpts%R_polar(3,1,n_start+j) = nwpts%R_max(n_start+j)
      nwpts%R_polar(3,4,n_start+j) = delta * nwpts%R_max(n_start+j) + (1.d0 - delta) * nwpts%R_sep(n_start+j)
      nwpts%R_polar(3,2,n_start+j) = ( 2.d0 * nwpts%R_polar(3,1,n_start+j)  +         nwpts%R_polar(3,4,n_start+j) ) / 3.d0
      nwpts%R_polar(3,3,n_start+j) = (        nwpts%R_polar(3,1,n_start+j)  +  2.d0 * nwpts%R_polar(3,4,n_start+j) ) / 3.d0

      nwpts%Z_polar(3,1,n_start+j) = nwpts%Z_max(n_start+j)
      nwpts%Z_polar(3,4,n_start+j) = delta * nwpts%Z_max(n_start+j) + (1.d0 - delta) * nwpts%Z_sep(n_start+j)
      nwpts%Z_polar(3,2,n_start+j) = ( 2.d0 * nwpts%Z_polar(3,1,n_start+j)  +         nwpts%Z_polar(3,4,n_start+j) ) / 3.d0
      nwpts%Z_polar(3,3,n_start+j) = (        nwpts%Z_polar(3,1,n_start+j)  +  2.d0 * nwpts%Z_polar(3,4,n_start+j) ) / 3.d0

      nwpts%R_polar(2,1,n_start+j) = nwpts%R_polar(1,4,n_start+j)
      nwpts%R_polar(2,4,n_start+j) = nwpts%R_polar(3,4,n_start+j)
      nwpts%R_polar(2,2,n_start+j) = ( nwpts%R_polar(2,1,n_start+j) +  2.d0 * nwpts%R_sep(n_start+j) ) / 3.d0
      nwpts%R_polar(2,3,n_start+j) = ( nwpts%R_polar(2,4,n_start+j) +  2.d0 * nwpts%R_sep(n_start+j) ) / 3.d0

      nwpts%Z_polar(2,1,n_start+j) = nwpts%Z_polar(1,4,n_start+j)
      nwpts%Z_polar(2,4,n_start+j) = nwpts%Z_polar(3,4,n_start+j)
      nwpts%Z_polar(2,2,n_start+j) = ( nwpts%Z_polar(2,1,n_start+j) +  2.d0 * nwpts%Z_sep(n_start+j) ) / 3.d0
      nwpts%Z_polar(2,3,n_start+j) = ( nwpts%Z_polar(2,4,n_start+j) +  2.d0 * nwpts%Z_sep(n_start+j) ) / 3.d0
    else
      nwpts%R_polar(1,1,n_start+j) = nwpts%R_min(n_start+j)
      nwpts%R_polar(1,4,n_start+j) = delta * nwpts%R_min(n_start+j) + (1.d0 - delta) * nwpts%R_sep(n_start+j)
      nwpts%R_polar(1,2,n_start+j) = ( 2.d0 * nwpts%R_polar(1,1,n_start+j)  +        nwpts%R_polar(1,4,n_start+j) ) / 3.d0
      nwpts%R_polar(1,3,n_start+j) = (        nwpts%R_polar(1,1,n_start+j)  +  2.d0 * nwpts%R_polar(1,4,n_start+j) ) / 3.d0

      nwpts%Z_polar(1,1,n_start+j) = nwpts%Z_min(n_start+j)
      nwpts%Z_polar(1,4,n_start+j) = delta * nwpts%Z_min(n_start+j) + (1.d0 - delta) * nwpts%Z_sep(n_start+j)
      nwpts%Z_polar(1,2,n_start+j) = ( 2.d0 * nwpts%Z_polar(1,1,n_start+j)  +        nwpts%Z_polar(1,4,n_start+j) ) / 3.d0
      nwpts%Z_polar(1,3,n_start+j) = (        nwpts%Z_polar(1,1,n_start+j)  +  2.d0 * nwpts%Z_polar(1,4,n_start+j) ) / 3.d0

      nwpts%R_polar(3,1,n_start+j) = nwpts%R_wall(j-n_leg)
      nwpts%R_polar(3,4,n_start+j) = delta * nwpts%R_wall(j-n_leg) + (1.d0 - delta) * nwpts%R_sep(n_start+j)
      nwpts%R_polar(3,2,n_start+j) = ( 2.d0 * nwpts%R_polar(3,1,n_start+j)  +        nwpts%R_polar(3,4,n_start+j) ) / 3.d0
      nwpts%R_polar(3,3,n_start+j) = (        nwpts%R_polar(3,1,n_start+j)  +  2.d0 * nwpts%R_polar(3,4,n_start+j) ) / 3.d0

      nwpts%Z_polar(3,1,n_start+j) = nwpts%Z_wall(j-n_leg)
      nwpts%Z_polar(3,4,n_start+j) = delta * nwpts%Z_wall(j-n_leg) + (1.d0 - delta) * nwpts%Z_sep(n_start+j)
      nwpts%Z_polar(3,2,n_start+j) = ( 2.d0 * nwpts%Z_polar(3,1,n_start+j)  +        nwpts%Z_polar(3,4,n_start+j) ) / 3.d0
      nwpts%Z_polar(3,3,n_start+j) = (        nwpts%Z_polar(3,1,n_start+j)  +  2.d0 * nwpts%Z_polar(3,4,n_start+j) ) / 3.d0

      nwpts%R_polar(2,1,n_start+j) = nwpts%R_polar(1,4,n_start+j)
      nwpts%R_polar(2,4,n_start+j) = nwpts%R_polar(3,4,n_start+j)
      nwpts%R_polar(2,2,n_start+j) = ( nwpts%R_polar(2,1,n_start+j) +  2.d0 * nwpts%R_sep(n_start+j) ) / 3.d0
      nwpts%R_polar(2,3,n_start+j) = ( nwpts%R_polar(2,4,n_start+j) +  2.d0 * nwpts%R_sep(n_start+j) ) / 3.d0

      nwpts%Z_polar(2,1,n_start+j) = nwpts%Z_polar(1,4,n_start+j)
      nwpts%Z_polar(2,4,n_start+j) = nwpts%Z_polar(3,4,n_start+j)
      nwpts%Z_polar(2,2,n_start+j) = ( nwpts%Z_polar(2,1,n_start+j) +  2.d0 * nwpts%Z_sep(n_start+j) ) / 3.d0
      nwpts%Z_polar(2,3,n_start+j) = ( nwpts%Z_polar(2,4,n_start+j) +  2.d0 * nwpts%Z_sep(n_start+j) ) / 3.d0    

      nwpts%R_polar(4,1,n_start+j) = nwpts%R_max(n_start+j)
      nwpts%R_polar(4,4,n_start+j) = nwpts%R_wall(j-n_leg)
      nwpts%R_polar(4,2,n_start+j) = ( 2.d0 * nwpts%R_polar(4,1,n_start+j)  +        nwpts%R_polar(4,4,n_start+j) ) / 3.d0
      nwpts%R_polar(4,3,n_start+j) = (        nwpts%R_polar(4,1,n_start+j)  +  2.d0 * nwpts%R_polar(4,4,n_start+j) ) / 3.d0

      nwpts%Z_polar(4,1,n_start+j) = nwpts%Z_max(n_start+j)
      nwpts%Z_polar(4,4,n_start+j) = nwpts%Z_wall(j-n_leg)
      nwpts%Z_polar(4,2,n_start+j) = ( 2.d0 * nwpts%Z_polar(4,1,n_start+j)  +        nwpts%Z_polar(4,4,n_start+j) ) / 3.d0
      nwpts%Z_polar(4,3,n_start+j) = (        nwpts%Z_polar(4,1,n_start+j)  +  2.d0 * nwpts%Z_polar(4,4,n_start+j) ) / 3.d0
    endif

  enddo
endif

!------------------------------ Upper leg
if(xcase .ne. LOWER_XPOINT) then
  do j=1,2*n_up_leg

    delta = 0.2
    if ( (j .eq. n_up_leg)   .or. (j .eq. 2*n_up_leg) )    delta = 0.d0
    if ( (j .eq. n_up_leg-1) .or. (j .eq. 2*n_up_leg-1) )  delta = 0.05d0
    if(xcase .eq. UPPER_XPOINT) n_start = n_tht
    if(xcase .eq. DOUBLE_NULL ) n_start = n_tht+2*n_leg

    if ( (j .le. n_up_leg) .or. (tokamak_device(1:4) .ne. 'MAST') ) then
      nwpts%R_polar(1,1,n_start+j) = nwpts%R_min(n_start+j)
      nwpts%R_polar(1,4,n_start+j) = delta * nwpts%R_min(n_start+j) + (1.d0 - delta) * nwpts%R_sep(n_start+j)
      nwpts%R_polar(1,2,n_start+j) = ( 2.d0 * nwpts%R_polar(1,1,n_start+j)  +         nwpts%R_polar(1,4,n_start+j) ) / 3.d0
      nwpts%R_polar(1,3,n_start+j) = (        nwpts%R_polar(1,1,n_start+j)  +  2.d0 * nwpts%R_polar(1,4,n_start+j) ) / 3.d0

      nwpts%Z_polar(1,1,n_start+j) = nwpts%Z_min(n_start+j)
      nwpts%Z_polar(1,4,n_start+j) = delta * nwpts%Z_min(n_start+j) + (1.d0 - delta) * nwpts%Z_sep(n_start+j)
      nwpts%Z_polar(1,2,n_start+j) = ( 2.d0 * nwpts%Z_polar(1,1,n_start+j)  +         nwpts%Z_polar(1,4,n_start+j) ) / 3.d0
      nwpts%Z_polar(1,3,n_start+j) = (        nwpts%Z_polar(1,1,n_start+j)  +  2.d0 * nwpts%Z_polar(1,4,n_start+j) ) / 3.d0

      nwpts%R_polar(3,1,n_start+j) = nwpts%R_max(n_start+j)
      nwpts%R_polar(3,4,n_start+j) = delta * nwpts%R_max(n_start+j) + (1.d0 - delta) * nwpts%R_sep(n_start+j)
      nwpts%R_polar(3,2,n_start+j) = ( 2.d0 * nwpts%R_polar(3,1,n_start+j)  +         nwpts%R_polar(3,4,n_start+j) ) / 3.d0
      nwpts%R_polar(3,3,n_start+j) = (        nwpts%R_polar(3,1,n_start+j)  +  2.d0 * nwpts%R_polar(3,4,n_start+j) ) / 3.d0

      nwpts%Z_polar(3,1,n_start+j) = nwpts%Z_max(n_start+j)
      nwpts%Z_polar(3,4,n_start+j) = delta * nwpts%Z_max(n_start+j) + (1.d0 - delta) * nwpts%Z_sep(n_start+j)
      nwpts%Z_polar(3,2,n_start+j) = ( 2.d0 * nwpts%Z_polar(3,1,n_start+j)  +         nwpts%Z_polar(3,4,n_start+j) ) / 3.d0
      nwpts%Z_polar(3,3,n_start+j) = (        nwpts%Z_polar(3,1,n_start+j)  +  2.d0 * nwpts%Z_polar(3,4,n_start+j) ) / 3.d0

      nwpts%R_polar(2,1,n_start+j) = nwpts%R_polar(1,4,n_start+j)
      nwpts%R_polar(2,4,n_start+j) = nwpts%R_polar(3,4,n_start+j)
      nwpts%R_polar(2,2,n_start+j) = ( nwpts%R_polar(2,1,n_start+j) +  2.d0 * nwpts%R_sep(n_start+j) ) / 3.d0
      nwpts%R_polar(2,3,n_start+j) = ( nwpts%R_polar(2,4,n_start+j) +  2.d0 * nwpts%R_sep(n_start+j) ) / 3.d0

      nwpts%Z_polar(2,1,n_start+j) = nwpts%Z_polar(1,4,n_start+j)
      nwpts%Z_polar(2,4,n_start+j) = nwpts%Z_polar(3,4,n_start+j)
      nwpts%Z_polar(2,2,n_start+j) = ( nwpts%Z_polar(2,1,n_start+j) +  2.d0 * nwpts%Z_sep(n_start+j) ) / 3.d0
      nwpts%Z_polar(2,3,n_start+j) = ( nwpts%Z_polar(2,4,n_start+j) +  2.d0 * nwpts%Z_sep(n_start+j) ) / 3.d0
    else
      nwpts%R_polar(1,1,n_start+j) = nwpts%R_min(n_start+j)
      nwpts%R_polar(1,4,n_start+j) = delta * nwpts%R_min(n_start+j) + (1.d0 - delta) * nwpts%R_sep(n_start+j)
      nwpts%R_polar(1,2,n_start+j) = ( 2.d0 * nwpts%R_polar(1,1,n_start+j)  +         nwpts%R_polar(1,4,n_start+j) ) / 3.d0
      nwpts%R_polar(1,3,n_start+j) = (        nwpts%R_polar(1,1,n_start+j)  +  2.d0 * nwpts%R_polar(1,4,n_start+j) ) / 3.d0

      nwpts%Z_polar(1,1,n_start+j) = nwpts%Z_min(n_start+j)
      nwpts%Z_polar(1,4,n_start+j) = delta * nwpts%Z_min(n_start+j) + (1.d0 - delta) * nwpts%Z_sep(n_start+j)
      nwpts%Z_polar(1,2,n_start+j) = ( 2.d0 * nwpts%Z_polar(1,1,n_start+j)  +         nwpts%Z_polar(1,4,n_start+j) ) / 3.d0
      nwpts%Z_polar(1,3,n_start+j) = (        nwpts%Z_polar(1,1,n_start+j)  +  2.d0 * nwpts%Z_polar(1,4,n_start+j) ) / 3.d0

      nwpts%R_polar(3,1,n_start+j) = nwpts%R_wall(j-n_up_leg+n_leg)
      nwpts%R_polar(3,4,n_start+j) = delta * nwpts%R_wall(j-n_up_leg+n_leg) + (1.d0 - delta) * nwpts%R_sep(n_start+j)
      nwpts%R_polar(3,2,n_start+j) = ( 2.d0 * nwpts%R_polar(3,1,n_start+j)  +         nwpts%R_polar(3,4,n_start+j) ) / 3.d0
      nwpts%R_polar(3,3,n_start+j) = (        nwpts%R_polar(3,1,n_start+j)  +  2.d0 * nwpts%R_polar(3,4,n_start+j) ) / 3.d0

      nwpts%Z_polar(3,1,n_start+j) = nwpts%Z_wall(j-n_up_leg+n_leg)
      nwpts%Z_polar(3,4,n_start+j) = delta * nwpts%Z_wall(j-n_up_leg+n_leg) + (1.d0 - delta) * nwpts%Z_sep(n_start+j)
      nwpts%Z_polar(3,2,n_start+j) = ( 2.d0 * nwpts%Z_polar(3,1,n_start+j)  +         nwpts%Z_polar(3,4,n_start+j) ) / 3.d0
      nwpts%Z_polar(3,3,n_start+j) = (        nwpts%Z_polar(3,1,n_start+j)  +  2.d0 * nwpts%Z_polar(3,4,n_start+j) ) / 3.d0

      nwpts%R_polar(2,1,n_start+j) = nwpts%R_polar(1,4,n_start+j)
      nwpts%R_polar(2,4,n_start+j) = nwpts%R_polar(3,4,n_start+j)
      nwpts%R_polar(2,2,n_start+j) = ( nwpts%R_polar(2,1,n_start+j) +  2.d0 * nwpts%R_sep(n_start+j) ) / 3.d0
      nwpts%R_polar(2,3,n_start+j) = ( nwpts%R_polar(2,4,n_start+j) +  2.d0 * nwpts%R_sep(n_start+j) ) / 3.d0

      nwpts%Z_polar(2,1,n_start+j) = nwpts%Z_polar(1,4,n_start+j)
      nwpts%Z_polar(2,4,n_start+j) = nwpts%Z_polar(3,4,n_start+j)
      nwpts%Z_polar(2,2,n_start+j) = ( nwpts%Z_polar(2,1,n_start+j) +  2.d0 * nwpts%Z_sep(n_start+j) ) / 3.d0
      nwpts%Z_polar(2,3,n_start+j) = ( nwpts%Z_polar(2,4,n_start+j) +  2.d0 * nwpts%Z_sep(n_start+j) ) / 3.d0

      nwpts%R_polar(4,1,n_start+j) = nwpts%R_max(n_start+j)
      nwpts%R_polar(4,4,n_start+j) = nwpts%R_wall(j-n_up_leg+n_leg)
      nwpts%R_polar(4,2,n_start+j) = ( 2.d0 * nwpts%R_polar(4,1,n_start+j)  +         nwpts%R_polar(4,4,n_start+j) ) / 3.d0
      nwpts%R_polar(4,3,n_start+j) = (        nwpts%R_polar(4,1,n_start+j)  +  2.d0 * nwpts%R_polar(4,4,n_start+j) ) / 3.d0

      nwpts%Z_polar(4,1,n_start+j) = nwpts%Z_max(n_start+j)
      nwpts%Z_polar(4,4,n_start+j) = nwpts%Z_wall(j-n_up_leg+n_leg)
      nwpts%Z_polar(4,2,n_start+j) = ( 2.d0 * nwpts%Z_polar(4,1,n_start+j)  +         nwpts%Z_polar(4,4,n_start+j) ) / 3.d0
      nwpts%Z_polar(4,3,n_start+j) = (        nwpts%Z_polar(4,1,n_start+j)  +  2.d0 * nwpts%Z_polar(4,4,n_start+j) ) / 3.d0
    endif

  enddo
endif


!------------------------------ Plot coordinate lines
if ( write_ps ) then
  call lincol(3)
  npl = 11
  call tr_allocate(xp,1,npl,"xp",CAT_GRID)
  call tr_allocate(yp,1,npl,"yp",CAT_GRID)
  do j=1,n_tht_2

    n_loop = 3
    if(tokamak_device(1:4) .eq. 'MAST') then
      if ((j .gt. n_tht+n_leg) .and. (j .le. n_tht+2*n_leg)) n_loop = 4
      if ((j .gt. n_tht+2*n_leg+n_up_leg) .and. (j .le. n_tht+2*n_leg+2*n_up_leg)) n_loop = 4
    endif
    do m=1,n_loop

      do k=1,npl
        ss = -1. + 2.*float(k-1)/float(npl-1)

        R_cub1d = (/ nwpts%R_polar(m,1,j), 3.d0/2.d0 *(nwpts%R_polar(m,2,j)-nwpts%R_polar(m,1,j)), &
            nwpts%R_polar(m,4,j), 3.d0/2.d0 *(nwpts%R_polar(m,4,j)-nwpts%R_polar(m,3,j))  /)
        Z_cub1d = (/ nwpts%Z_polar(m,1,j), 3.d0/2.d0 *(nwpts%Z_polar(m,2,j)-nwpts%Z_polar(m,1,j)), &
            nwpts%Z_polar(m,4,j), 3.d0/2.d0 *(nwpts%Z_polar(m,4,j)-nwpts%Z_polar(m,3,j)) /)

        call CUB1D(R_cub1d(1), R_cub1d(2), R_cub1d(3), R_cub1d(4),ss,xp(k), tmp1)
        call CUB1D(Z_cub1d(1), Z_cub1d(2), Z_cub1d(3), Z_cub1d(4),ss,yp(k), tmp2)
      enddo
      call lincol(3)
      write(51,*) ' .1 setlinewidth'
      call lplot6(1,1,xp,yp,-npl,' ')
      write(51,*) ' stroke'

    enddo

    call lincol(0)

  enddo
  call tr_deallocate(xp,"xp",CAT_GRID)
  call tr_deallocate(yp,"yp",CAT_GRID)
endif





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
do i=1,n_flux+n_open
  do j=1, n_tht
    
    do k=1,3     ! 3 line pieces per coordinate line

      R_cub1d = (/ nwpts%R_polar(k,1,j), 3.d0/2.d0 *(nwpts%R_polar(k,2,j)-nwpts%R_polar(k,1,j)), &
                   nwpts%R_polar(k,4,j), 3.d0/2.d0 *(nwpts%R_polar(k,4,j)-nwpts%R_polar(k,3,j))  /)
      Z_cub1d = (/ nwpts%Z_polar(k,1,j), 3.d0/2.d0 *(nwpts%Z_polar(k,2,j)-nwpts%Z_polar(k,1,j)), &
                   nwpts%Z_polar(k,4,j), 3.d0/2.d0 *(nwpts%Z_polar(k,4,j)-nwpts%Z_polar(k,3,j)) /)

      call find_crossing(node_list,element_list,flux_list,i,R_cub1d,Z_cub1d, &
                      nwpts%RR_new(i+1,j),nwpts%ZZ_new(i+1,j),nwpts%ielm_flux(i+1,j),nwpts%s_flux(i+1,j),nwpts%t_flux(i+1,j),nwpts%t_tht(i+1,j),ifail,.false.)
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
      ! We avoid last surface, points are already known...
      if ( (i .eq. n_flux+n_open) .and. (xcase .ne. DOUBLE_NULL) ) then
        nwpts%k_cross(i+1,j) = 3
        nwpts%RR_new(i+1,j)  = nwpts%R_max(j)
        nwpts%ZZ_new(i+1,j)  = nwpts%Z_max(j)
        write(*,*) ' WARNING node not found for last open flux surface -> using RZ_max '
      else
        write(*,'(A,I6,I6,I6,F20.10)') ' WARNING node not found for central grid (without legs) : ',ifail,i,j,theta_sep(j)
      endif
    endif
      
  enddo
enddo

!----------------------------------- The main inner/outer parts (without legs)
if (xcase .eq. DOUBLE_NULL) then
  do i=n_flux+n_open+1,n_flux+n_open+n_outer+n_inner
    if (     ( ( (ES%active_xpoint .eq. LOWER_XPOINT) .or. (ES%active_xpoint .eq. SYMMETRIC_XPOINT) ) .and. (i .lt. n_flux+n_open+n_outer+1) )  &
        .or. ( (  ES%active_xpoint .eq. UPPER_XPOINT                                         ) .and. (i .ge. n_flux+n_open+n_outer+1) )  ) then
      n_start = 1
      n_loop  = n_tht_mid
    else
      n_start = n_tht_mid + 1
      n_loop  = n_tht
    endif
    do j=n_start, n_loop
      
      do k=1,3     ! 3 line pieces per coordinate line

        R_cub1d = (/ nwpts%R_polar(k,1,j), 3.d0/2.d0 *(nwpts%R_polar(k,2,j)-nwpts%R_polar(k,1,j)), &
                     nwpts%R_polar(k,4,j), 3.d0/2.d0 *(nwpts%R_polar(k,4,j)-nwpts%R_polar(k,3,j))  /)
        Z_cub1d = (/ nwpts%Z_polar(k,1,j), 3.d0/2.d0 *(nwpts%Z_polar(k,2,j)-nwpts%Z_polar(k,1,j)), &
                     nwpts%Z_polar(k,4,j), 3.d0/2.d0 *(nwpts%Z_polar(k,4,j)-nwpts%Z_polar(k,3,j)) /)

        call find_crossing(node_list,element_list,flux_list,i,R_cub1d,Z_cub1d, &
                        nwpts%RR_new(i+1,j),nwpts%ZZ_new(i+1,j),nwpts%ielm_flux(i+1,j),nwpts%s_flux(i+1,j),nwpts%t_flux(i+1,j),nwpts%t_tht(i+1,j),ifail,.false.)
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
        ! We avoid last surface, points are already known...
        if ( (i .eq. n_flux+n_open+n_outer) .or. (i .eq. n_flux+n_open+n_outer+n_inner) ) then
          nwpts%k_cross(i+1,j) = 3
          nwpts%RR_new(i+1,j)  = nwpts%R_max(j)
          nwpts%ZZ_new(i+1,j)  = nwpts%Z_max(j)
          write(*,*) ' WARNING node not found for last open flux surface -> using RZ_max '
        else
          write(*,'(A,I6,I6,I6,F20.10)') ' WARNING node not found for central grid (without legs) : ',ifail,i,j,theta_sep(j)
        endif
      endif
        
    enddo
  enddo
endif

!----------------------------------- The legs and private parts
if(xcase .ne. DOUBLE_NULL) then
  do i=n_flux,n_psi-1          ! With the private parts
    do j=n_tht+1, n_tht_2
      do k=1,4      ! 3 line pieces per coordinate line (4 for MAST)

        R_cub1d = (/ nwpts%R_polar(k,1,j), 3.d0/2.d0 *(nwpts%R_polar(k,2,j)-nwpts%R_polar(k,1,j)), &
                     nwpts%R_polar(k,4,j), 3.d0/2.d0 *(nwpts%R_polar(k,4,j)-nwpts%R_polar(k,3,j))  /)
        Z_cub1d = (/ nwpts%Z_polar(k,1,j), 3.d0/2.d0 *(nwpts%Z_polar(k,2,j)-nwpts%Z_polar(k,1,j)), &
                     nwpts%Z_polar(k,4,j), 3.d0/2.d0 *(nwpts%Z_polar(k,4,j)-nwpts%Z_polar(k,3,j)) /)

        call find_crossing(node_list,element_list,flux_list,i,R_cub1d,Z_cub1d, &
                           nwpts%RR_new(i+1,j),nwpts%ZZ_new(i+1,j),nwpts%ielm_flux(i+1,j),nwpts%s_flux(i+1,j),nwpts%t_flux(i+1,j),nwpts%t_tht(i+1,j),ifail,.false.)
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
        ! We avoid last surface, points are already known...
        if (i .eq. n_psi-1) then
          nwpts%k_cross(i+1,j) = 1
          nwpts%RR_new(i+1,j)  = nwpts%R_min(j)
          nwpts%ZZ_new(i+1,j)  = nwpts%Z_min(j)
          write(*,*) ' WARNING node not found for last private surface -> using RZ_max '
        else
          if (xcase .eq. LOWER_XPOINT) write(*,'(A,i6,i6,i6)') ' WARNING node not found for lower part of grid : ',ifail,i,j
          if (xcase .eq. UPPER_XPOINT) write(*,'(A,i6,i6,i6)') ' WARNING node not found for upper part of grid : ',ifail,i,j
        endif
      endif

    enddo
  enddo
else  ! xcase == DOUBLE_NULL
  do i=n_flux,n_flux+n_open        ! The sandwich parts
    if ( ES%active_xpoint .eq. UPPER_XPOINT ) then
      n_start = n_tht+2*n_leg+1
      n_loop  = n_tht_2
    else
      n_start = n_tht+1
      n_loop  = n_tht_2-2*n_up_leg
    endif
    do j=n_start,n_loop
      do k=1,4      ! 3 line pieces per coordinate line (4 for MAST)

        R_cub1d = (/ nwpts%R_polar(k,1,j), 3.d0/2.d0 *(nwpts%R_polar(k,2,j)-nwpts%R_polar(k,1,j)), &
                     nwpts%R_polar(k,4,j), 3.d0/2.d0 *(nwpts%R_polar(k,4,j)-nwpts%R_polar(k,3,j))  /)
        Z_cub1d = (/ nwpts%Z_polar(k,1,j), 3.d0/2.d0 *(nwpts%Z_polar(k,2,j)-nwpts%Z_polar(k,1,j)), &
                     nwpts%Z_polar(k,4,j), 3.d0/2.d0 *(nwpts%Z_polar(k,4,j)-nwpts%Z_polar(k,3,j)) /)

        call find_crossing(node_list,element_list,flux_list,i,R_cub1d,Z_cub1d, &
                           nwpts%RR_new(i+1,j),nwpts%ZZ_new(i+1,j),nwpts%ielm_flux(i+1,j),nwpts%s_flux(i+1,j),nwpts%t_flux(i+1,j),nwpts%t_tht(i+1,j),ifail,.false.)
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

      if (ifail .ne. 0) write(*,'(A,i6,i6,i6)') ' WARNING node not found for the sandwich part of grid : ',ifail,i,j

    enddo
  enddo
  do i=n_flux+n_open,n_psi-n_private-n_up_priv-1          ! Outer/Inner parts
    do j=n_tht+1, n_tht_2
      if (      (i .gt. n_flux+n_open+n_outer) .and. (i .gt. n_flux+n_open)                                             &
          .and. ( (j .gt. n_tht+2*n_leg+n_up_leg) .or. ( (j .gt. n_tht+n_leg)   .and. (j .le. n_tht+2*n_leg) )          ) ) cycle
      if (       (i .le. n_flux+n_open+n_outer) .and. (i .gt. n_flux+n_open)                                            &
          .and. ( (j .le. n_tht+n_leg)            .or. ( (j .gt. n_tht+2*n_leg) .and. (j .le. n_tht+2*n_leg+n_up_leg) ) ) ) cycle
      do k=1,4      ! 3 line pieces per coordinate line (4 for MAST)

        R_cub1d = (/ nwpts%R_polar(k,1,j), 3.d0/2.d0 *(nwpts%R_polar(k,2,j)-nwpts%R_polar(k,1,j)), &
                     nwpts%R_polar(k,4,j), 3.d0/2.d0 *(nwpts%R_polar(k,4,j)-nwpts%R_polar(k,3,j))  /)
        Z_cub1d = (/ nwpts%Z_polar(k,1,j), 3.d0/2.d0 *(nwpts%Z_polar(k,2,j)-nwpts%Z_polar(k,1,j)), &
                     nwpts%Z_polar(k,4,j), 3.d0/2.d0 *(nwpts%Z_polar(k,4,j)-nwpts%Z_polar(k,3,j)) /)

        call find_crossing(node_list,element_list,flux_list,i,R_cub1d,Z_cub1d, &
                           nwpts%RR_new(i+1,j),nwpts%ZZ_new(i+1,j),nwpts%ielm_flux(i+1,j),nwpts%s_flux(i+1,j),nwpts%t_flux(i+1,j),nwpts%t_tht(i+1,j),ifail,.false.)
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
        ! We avoid last surface, points are already known...
        if ( (i .eq. n_flux+n_open+n_outer) .or. (i .eq. n_flux+n_open+n_outer+n_inner) ) then
          nwpts%k_cross(i+1,j) = 3
          nwpts%RR_new(i+1,j)  = nwpts%R_max(j)
          nwpts%ZZ_new(i+1,j)  = nwpts%Z_max(j)
          write(*,*) ' WARNING node not found for last open flux surface -> using RZ_max '
        else
          write(*,'(A,i6,i6,i6)') ' WARNING node not found for the outer/inner part of the legs : ',ifail,i,j
        endif
      endif

    enddo
  enddo
  do i=n_psi-n_private-n_up_priv,n_psi-n_up_priv-1          ! Lower private
    do j=n_tht+1, n_tht_2-2*n_up_leg
      do k=1,4      ! 3 line pieces per coordinate line (4 for MAST)

        R_cub1d = (/ nwpts%R_polar(k,1,j), 3.d0/2.d0 *(nwpts%R_polar(k,2,j)-nwpts%R_polar(k,1,j)), &
                     nwpts%R_polar(k,4,j), 3.d0/2.d0 *(nwpts%R_polar(k,4,j)-nwpts%R_polar(k,3,j))  /)
        Z_cub1d = (/ nwpts%Z_polar(k,1,j), 3.d0/2.d0 *(nwpts%Z_polar(k,2,j)-nwpts%Z_polar(k,1,j)), &
                     nwpts%Z_polar(k,4,j), 3.d0/2.d0 *(nwpts%Z_polar(k,4,j)-nwpts%Z_polar(k,3,j)) /)

        call find_crossing(node_list,element_list,flux_list,i,R_cub1d,Z_cub1d, &
                           nwpts%RR_new(i+1,j),nwpts%ZZ_new(i+1,j),nwpts%ielm_flux(i+1,j),nwpts%s_flux(i+1,j),nwpts%t_flux(i+1,j),nwpts%t_tht(i+1,j),ifail,.false.)
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
        ! We avoid last surface, points are already known...
        if (i .eq. n_psi-n_up_priv-1) then
          nwpts%k_cross(i+1,j) = 1
          nwpts%RR_new(i+1,j)  = nwpts%R_min(j)
          nwpts%ZZ_new(i+1,j)  = nwpts%Z_min(j)
          write(*,*) ' WARNING node not found for last lower private surface -> using RZ_max '
        else
          write(*,'(A,i6,i6,i6)') ' WARNING node not found lower private part of the grid : ',ifail,i,j
        endif
      endif

    enddo
  enddo
  do i=n_psi-n_up_priv,n_psi-1          !Upper private 
    do j=n_tht+2*n_leg+1, n_tht_2
      do k=1,4      ! 3 line pieces per coordinate line (4 for MAST)

        R_cub1d = (/ nwpts%R_polar(k,1,j), 3.d0/2.d0 *(nwpts%R_polar(k,2,j)-nwpts%R_polar(k,1,j)), &
                     nwpts%R_polar(k,4,j), 3.d0/2.d0 *(nwpts%R_polar(k,4,j)-nwpts%R_polar(k,3,j))  /)
        Z_cub1d = (/ nwpts%Z_polar(k,1,j), 3.d0/2.d0 *(nwpts%Z_polar(k,2,j)-nwpts%Z_polar(k,1,j)), &
                     nwpts%Z_polar(k,4,j), 3.d0/2.d0 *(nwpts%Z_polar(k,4,j)-nwpts%Z_polar(k,3,j)) /)

        call find_crossing(node_list,element_list,flux_list,i,R_cub1d,Z_cub1d, &
                           nwpts%RR_new(i+1,j),nwpts%ZZ_new(i+1,j),nwpts%ielm_flux(i+1,j),nwpts%s_flux(i+1,j),nwpts%t_flux(i+1,j),nwpts%t_tht(i+1,j),ifail,.false.)
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
        ! We avoid last surface, points are already known...
        if (i .eq. n_psi-1) then
          nwpts%k_cross(i+1,j) = 1
          nwpts%RR_new(i+1,j)  = nwpts%R_min(j)
          nwpts%ZZ_new(i+1,j)  = nwpts%Z_min(j)
          write(*,*) ' WARNING node not found for last upper private surface -> using RZ_max '
        else
          write(*,'(A,i6,i6,i6)') ' WARNING node not found upper private part of the grid : ',ifail,i,j
        endif
      endif

    enddo
  enddo
endif




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

!----------------------------------- Print a python file that plots the new grid points
if (plot_grid) then
  open(102,file='plot_new_points.py')
    write(102,'(A)')        '#!/usr/bin/env python'
    write(102,'(A)')        'import numpy as N'
    write(102,'(A)')        'import pylab'
    write(102,'(A)')        'def main():'
    write(102,'(A,i6,A)')     ' r = N.zeros(',(n_psi)*n_tht_2,')'
    write(102,'(A,i6,A)')     ' z = N.zeros(',(n_psi)*n_tht_2,')'
    do i=1,n_psi
      do j=1,n_tht_2
      write(102,'(A,i6,A,f15.4)') ' r[',(i-1)*(n_tht_2)+j-1,'] = ',nwpts%RR_new(i,j)
      write(102,'(A,i6,A,f15.4)') ' z[',(i-1)*(n_tht_2)+j-1,'] = ',nwpts%ZZ_new(i,j)
      enddo
    enddo
    write(102,'(A,i6,A,i6,A)')' pylab.plot(r[0:',(n_psi)*n_tht_2,'],z[0:',(n_psi)*n_tht_2,'], "r.")'
    write(102,'(A)')        ' pylab.axis("equal")'
    write(102,'(A)')        ' pylab.show()'
    write(102,'(A)')        ' '
    write(102,'(A)')        'main()'
  close(102)
endif

return
end subroutine define_new_grid_points
