subroutine grid_xpoint(node_list, element_list, n_flux, n_open, n_private, n_leg, n_tht, &
  SIG_open, SIG_closed, SIG_private, SIG_theta, SIG_leg_0, SIG_leg_1, dPSI_open, dPSI_private, xcase)
!-----------------------------------------------------------------------
! subroutine defines a flux surface aligned finite element grid
! inclduing a single x-point
!-----------------------------------------------------------------------

use constants
use tr_module 
use data_structure
use mod_neighbours, only: update_neighbours
use mod_interp
use phys_module, only: force_central_node, write_ps, fix_axis_nodes, treat_axis

implicit none

! --- Routine parameters
type (type_node_list),    intent(inout) :: node_list
type (type_element_list), intent(inout) :: element_list
integer,                  intent(in)    :: n_flux, n_open, n_private, n_leg, n_tht, xcase
real*8,                   intent(in)    :: SIG_open, SIG_closed, SIG_private, SIG_theta, SIG_leg_0, SIG_leg_1
real*8,                   intent(in)    :: dPSI_open, dPSI_private

! --- local variables
type (type_surface_list) :: flux_list

type (type_node_list),    pointer :: newnode_list
type (type_element_list), pointer :: newelement_list

real*8, allocatable :: s_values(:), theta_sep(:), R_sep(:), Z_sep(:), R_max(:), Z_max(:), R_min(:), Z_min(:),s_tmp(:)
real*8              :: psi_axis, R_axis, Z_axis, s_axis, t_axis, R_xpoint(2), Z_xpoint(2), s_xpoint(2), t_xpoint(2), psi_xpoint(2)
real*8              :: s_find(8), t_find(8), st_find(8), tht_x, theta, delta, ss, tmp1, tmp2
real*8              :: RRg1,dRRg1_dr,dRRg1_ds
real*8              :: ZZg1,dZZg1_dr,dZZg1_ds
real*8              :: PSg1,dPSg1_dr,dPSg1_ds,dPSg1_drs,dPSg1_drr,dPSg1_dss
real*8,allocatable  :: R_polar(:,:,:),Z_polar(:,:,:),xout(:),xp(:),yp(:)
real*8              :: R_cub1d(4), Z_cub1d(4), dR_dt, dZ_dt, RZ_jac, PSI_R, PSI_Z
real*8, allocatable :: RR_new(:,:),ZZ_new(:,:),s_flux(:,:),t_flux(:,:),t_tht(:,:)
integer,allocatable :: ielm_flux(:,:), keep(:,:,:), k_cross(:,:)
integer             :: i, j, k, l, m, n_psi, n_flux_2, n_open_2, n_tht_2, n_psi_2, i2, j2
integer             :: n_private_2, i_surf, n_pieces
integer             :: i_elm_axis, i_elm_xpoint(2), i_elm_find(8), i_sep, i_max, i_find, npl, ifail
integer             :: node, index, node_start, index_xpoint, n_xpoint, j_start, j_end
integer             :: iv, ivp, node_iv, node_ivp, i_elm, n_leg_2, n_tht_3
integer             :: my_id, ielm_out
real*8              :: Rmid, Zmid, R0,Z0, RP,ZP, dR0, dZ0, dRP, dZP, size_0, size_p, denom
real*8              :: R1, Z1, s_out, t_out, R_out, Z_out
real*8              :: EJAC, RX, RY, SX, SY, CRR, CZZ, CRZ, alpha1, alpha2, alpha_max, alpha_min
real*8              :: RL1, RL2, RL3, RL4, RL5, RL6, RL7, RL8, RL9, ZL1, ZL2, ZL3, ZL4, ZL5, ZL6, ZL7, ZL8, ZL9
real*8              :: angle_L1, angle_L8, angle_L9, rr1, ss1
logical             :: xpoint
real*8,external     :: root
character*4         :: label
integer             :: i_elm1, i_vertex1, i_node1, i_node_save
integer             :: i_elm2, i_vertex2, i_node2

xpoint = .true.
my_id  = 0

write(*,*) '*************************************'
write(*,*) '*          X-point grid             *'
write(*,*) '*************************************'


call find_axis(my_id,node_list,element_list,psi_axis,R_axis,Z_axis,i_elm_axis,s_axis,t_axis,ifail)

call find_xpoint(my_id,node_list,element_list,psi_xpoint,R_xpoint,Z_xpoint,i_elm_xpoint,s_xpoint,t_xpoint,xcase,ifail)

!-------------------------------- double values to find the mid_points
n_flux_2    = 2 * (n_flux - 1)
n_open_2    = 2 * n_open
n_tht_2     = 2 * (n_tht - 1) + 1
n_private_2 = 2 * n_private
n_leg_2     = 2 * n_leg - 1

write(*,*) ' n_flux,   n_open,   n_tht   : ', n_flux,   n_open,   n_tht
write(*,*) ' n_flux_2, n_open_2, n_tht_2 : ', n_flux_2, n_open_2, n_tht_2

flux_list%n_psi = n_flux_2 + n_open_2 + n_private_2

call tr_allocate(flux_list%psi_values,1,flux_list%n_psi,"flux_list%psi_values",CAT_GRID)
call tr_allocate(s_values,1,n_flux_2+n_open_2+n_private_2,"s_values",CAT_GRID)

call tr_allocate(s_tmp,1,n_flux_2+1,"s_tmp",CAT_GRID)
s_tmp = 0
call meshac2(n_flux_2+1,s_tmp,1.d0,9999.d0,SIG_closed,9999.d0,0.2d0,1.0d0)

do i=1,n_flux_2
  s_values(i) = s_tmp(i+1)
  flux_list%psi_values(i) =  psi_axis + s_values(i)**2 * (psi_xpoint(1) - psi_axis)
enddo

call tr_deallocate(s_tmp,"s_tmp",CAT_GRID); call tr_allocate(s_tmp,1,n_open_2+1,"s_tmp",CAT_GRID)
s_tmp = 0
call meshac2(n_open_2+1,s_tmp,0.d0,9999.d0,SIG_open,9999.d0,0.6d0,1.0d0)

do i=1,n_open_2
  s_values(i+n_flux_2)             =  1.d0 + dPSI_open*s_tmp(i+1)
  flux_list%psi_values(i+n_flux_2) =  psi_axis + s_values(i+n_flux_2)**2 * (psi_xpoint(1) - psi_axis)
enddo

call tr_deallocate(s_tmp,"s_tmp",CAT_GRID); call tr_allocate(s_tmp,1,n_private_2+1,"s_tmp",CAT_GRID)
s_tmp = 0
call meshac2(n_private_2+1,s_tmp,0.d0,9999.d0,SIG_private,9999.d0,0.6d0,1.0d0)

do i=1,n_private_2
  s_values(i+n_flux_2+n_open_2)             =  1.d0 - dPSI_private*s_tmp(i+1)
  flux_list%psi_values(i+n_flux_2+n_open_2) =  psi_axis + s_values(i+n_flux_2+n_open_2)**2 * (psi_xpoint(1) - psi_axis)
enddo


call find_flux_surfaces(my_id,xpoint,xcase,node_list,element_list,flux_list)
!call q_profile(node_list,element_list,flux_list,psi_axis,psi_xpoint,Z_xpoint)

call plot_flux_surfaces(node_list,element_list,flux_list,.true.,1,.true.,xcase)


!-------------------------------------- store some data for the new grid

n_psi   = n_flux   + n_open   + n_private
n_psi_2 = n_flux_2 + n_open_2 + n_private_2 + 1                            ! this includes the magnetic axis

n_tht_3 = n_tht_2 + 2*n_leg_2

call tr_allocate(RR_new,1,n_psi_2,1,n_tht_3,"RR_new",CAT_GRID)
call tr_allocate(ZZ_new,1,n_psi_2,1,n_tht_3,"ZZ_new",CAT_GRID)
call tr_allocate(ielm_flux,1,n_psi_2,1,n_tht_3,"ielm_flux",CAT_GRID)
call tr_allocate(s_flux,1,n_psi_2,1,n_tht_3,"s_flux",CAT_GRID)
call tr_allocate(t_flux,1,n_psi_2,1,n_tht_3,"t_flux",CAT_GRID)
call tr_allocate(t_tht,1,n_psi_2,1,n_tht_3,"t_tht",CAT_GRID)

!------------------------------------- find some points on the legs
RL1 = 999.; ZL1 = 1.d10        ! left  bottom point of inner leg (i.e on last flux_surface)
RL4 = 999.; ZL4 = 1.d10        ! right bottom point of outer leg (i.e on last flux_surface)

i_surf = n_flux_2+n_open_2

do k=1,flux_list%flux_surfaces(i_surf)%n_pieces

  rr1   = flux_list%flux_surfaces(i_surf)%s(1,k)
  ss1   = flux_list%flux_surfaces(i_surf)%t(1,k)
  i_elm = flux_list%flux_surfaces(i_surf)%elm(k)

  call interp_RZ(node_list,element_list,i_elm,rr1,ss1,RRg1,ZZg1)

  if ((ZZg1 .lt. ZL1) .and. (RRg1 .lt. R_xpoint(1))) then
    RL1 = RRg1
    ZL1 = ZZg1
  endif
  if ((ZZg1 .lt. ZL4) .and. (RRg1 .gt. R_xpoint(1))) then
    RL4 = RRg1
    ZL4 = ZZg1
  endif

  rr1   = flux_list%flux_surfaces(i_surf)%s(3,k)
  ss1   = flux_list%flux_surfaces(i_surf)%t(3,k)
  i_elm = flux_list%flux_surfaces(i_surf)%elm(k)

  call interp_RZ(node_list,element_list,i_elm,rr1,ss1,RRg1,ZZg1)

  if ((ZZg1 .lt. ZL1) .and. (RRg1 .lt. R_xpoint(1))) then
    RL1 = RRg1
    ZL1 = ZZg1
  endif
  if ((ZZg1 .lt. ZL4) .and. (RRg1 .gt. R_xpoint(1))) then
    RL4 = RRg1
    ZL4 = ZZg1
  endif

enddo

RL2 =  1.d10; ZL2 = 999.      ! right bottom point of inner leg (i.e on last private flux_surface)
RL3 = -1.d10; ZL3 = 999.      ! left  bottom point of outer leg (i.e on last private flux_surface)

i_surf = n_psi_2 - 1

do k=1,flux_list%flux_surfaces(i_surf)%n_pieces

  rr1   = flux_list%flux_surfaces(i_surf)%s(1,k)
  ss1   = flux_list%flux_surfaces(i_surf)%t(1,k)
  i_elm = flux_list%flux_surfaces(i_surf)%elm(k)

  call interp_RZ(node_list,element_list,i_elm,rr1,ss1,RRg1,ZZg1)

  if ((RRg1 .lt. RL2) .and. (ZZg1 .lt. Z_xpoint(1))) then
    RL2 = RRg1
    ZL2 = ZZg1
  endif
  if ((RRg1 .gt. RL3)  .and. (ZZg1 .lt. Z_xpoint(1))) then
    RL3 = RRg1
    ZL3 = ZZg1
  endif

  rr1   = flux_list%flux_surfaces(i_surf)%s(3,k)
  ss1   = flux_list%flux_surfaces(i_surf)%t(3,k)
  i_elm = flux_list%flux_surfaces(i_surf)%elm(k)

  call interp_RZ(node_list,element_list,i_elm,rr1,ss1,RRg1,ZZg1)

  if ((RRg1 .lt. RL2) .and. (ZZg1 .lt. Z_xpoint(1))) then
    RL2 = RRg1
    ZL2 = ZZg1
  endif
  if ((RRg1 .gt. RL3)  .and. (ZZg1 .lt. Z_xpoint(1))) then
    RL3 = RRg1
    ZL3 = ZZg1
  endif
enddo

tht_x = atan2(Z_xpoint(1)-Z_axis,R_xpoint(1)-R_axis)
if (tht_x .lt. 0.d0) tht_x = tht_x + 2.d0 * PI

call find_theta_surface(node_list,element_list,flux_list,flux_list%n_psi,tht_x,R_axis,Z_axis,i_elm_find,s_find,t_find,i_find)

do i=1,i_find
   call interp_RZ(node_list,element_list,i_elm_find(i),s_find(i),t_find(i),RL5,ZL5)

   if (ZL5 .le. Z_xpoint(1)) exit
enddo

RL6 = 999.; ZL6 = 1.d10  ! left  strike point (i.e. on separatrix)
RL7 = 999.; ZL7 = 1.d10  ! right strike point (i.e on separatrix)

i_surf = n_flux_2

do k=1,flux_list%flux_surfaces(i_surf)%n_pieces

  rr1   = flux_list%flux_surfaces(i_surf)%s(1,k)
  ss1   = flux_list%flux_surfaces(i_surf)%t(1,k)
  i_elm = flux_list%flux_surfaces(i_surf)%elm(k)

  call interp_RZ(node_list,element_list,i_elm,rr1,ss1,RRg1,ZZg1)

  if ((ZZg1 .lt. ZL6) .and. (RRg1 .lt. R_xpoint(1))) then
    RL6 = RRg1
    ZL6 = ZZg1
  endif
  if ((ZZg1 .lt. ZL7)  .and. (RRg1 .gt. R_xpoint(1))) then
    RL7 = RRg1
    ZL7 = ZZg1
  endif

  rr1   = flux_list%flux_surfaces(i_surf)%s(3,k)
  ss1   = flux_list%flux_surfaces(i_surf)%t(3,k)
  i_elm = flux_list%flux_surfaces(i_surf)%elm(k)

  call interp_RZ(node_list,element_list,i_elm,rr1,ss1,RRg1,ZZg1)

  if ((ZZg1 .lt. ZL6) .and. (RRg1 .lt. R_xpoint(1))) then
    RL6 = RRg1
    ZL6 = ZZg1
  endif
  if ((ZZg1 .lt. ZL7)  .and. (RRg1 .gt. R_xpoint(1))) then
    RL7 = RRg1
    ZL7 = ZZg1
  endif

enddo
write(*,*) ' LEG points (1) : ',RL1,ZL1
write(*,*) ' LEG points (2) : ',RL2,ZL2
write(*,*) ' LEG points (3) : ',RL3,ZL3
write(*,*) ' LEG points (4) : ',RL4,ZL4
write(*,*) ' LEG points (5) : ',RL5,ZL5
write(*,*) ' LEG points (6) : ',RL6,ZL6
write(*,*) ' LEG points (7) : ',RL7,ZL7


!-------------------------------- define the polar coordinate
i_sep = n_flux_2
i_max = n_flux_2 + n_open_2

tht_x = atan2(Z_xpoint(1)-Z_axis,R_xpoint(1)-R_axis)

call tr_allocate(theta_sep,1,n_tht_3,"theta_sep",CAT_GRID)
call tr_allocate(R_sep,1,n_tht_3,"R_sep",CAT_GRID)
call tr_allocate(Z_sep,1,n_tht_3,"Z_sep",CAT_GRID)

call tr_deallocate(s_tmp,"s_tmp",CAT_GRID); call tr_allocate(s_tmp,1,n_tht_2,"s_tmp",CAT_GRID)
s_tmp = 0
call meshac2(n_tht_2,s_tmp,0.d0,1.d0,SIG_theta,SIG_theta,0.8d0,1.0d0)

!call meshac3(n_tht_2,s_tmp,0.d0,-tht_x/(2.d0*PI),1.d0,SIG_theta,0.05d0,SIG_theta,0.8d0,1.0d0)

do j=1,n_tht_2

  theta_sep(j) = tht_x + 2.d0 * PI * s_tmp(j)

  if (theta_sep(j) .lt. 0.d0) theta_sep(j) = theta_sep(j) + 2.d0 * PI

enddo

!------------------------------------- find crossing with separatrix
do j=1,n_tht_2

!  write(*,*) ' theta_sep : ',j,theta_sep(j)

  call find_theta_surface(node_list,element_list,flux_list,i_sep,theta_sep(j),R_axis,Z_axis,i_elm_find,s_find,t_find,i_find)

  call interp_RZ(node_list,element_list,i_elm_find(1),s_find(1),t_find(1),&
                 RRg1,ZZg1)

  R_sep(j) = RRg1
  Z_sep(j) = ZZg1
!  write(*,'(A,i5,2e16.8)') ' R_sep, Z_sep : ',j,R_sep(j),Z_sep(j)
enddo

R_sep(1) = R_xpoint(1)
Z_sep(1) = Z_xpoint(1)
R_sep(n_tht_2) = R_xpoint(1)
Z_sep(n_tht_2) = Z_xpoint(1)

angle_L1 = atan2(ZL1-Z_xpoint(1),RL1-R_xpoint(1))
if (angle_L1 .lt. 0.d0) angle_L1 = angle_L1 + 2.d0*PI

angle_L8 = tht_x + 1.5d0*PI
angle_L9 = tht_x + 0.5d0*PI

if (tht_x + 1.5d0*PI .gt. angle_L1) then   ! check if L8 is above L1 : if not adjust angle
  angle_L8 = angle_L1 - 0.1
  angle_L9 = angle_L8 - PI
endif

!------------------------------------ find crossing of last fluxsurface
call tr_allocate(R_max,1,n_tht_3,"R_max",CAT_GRID)
call tr_allocate(Z_max,1,n_tht_3,"Z_max",CAT_GRID)
call tr_allocate(R_min,1,n_tht_3,"R_min",CAT_GRID)
call tr_allocate(Z_min,1,n_tht_3,"Z_min",CAT_GRID)

call find_theta_surface(node_list,element_list,flux_list,i_max,angle_L9,R_xpoint(1),Z_xpoint(1),i_elm_find,s_find,t_find,i_find)
call interp_RZ(node_list,element_list,i_elm_find(1),s_find(1),t_find(1),RL9,ZL9)

write(*,'(A,2e16.8)') ' L9 (outside): ',RL9,ZL9  ! outside wall crossing at x-point level

call find_theta_surface(node_list,element_list,flux_list,i_max,angle_L8,R_xpoint(1),Z_xpoint(1),i_elm_find,s_find,t_find,i_find)
call interp_RZ(node_list,element_list,i_elm_find(1),s_find(1),t_find(1),RL8,ZL8)

write(*,'(A,2e16.8)') ' L8 (inside): ',RL8,ZL8  ! inside wall crossing at x-point level


do j=1,n_tht_2

  if (Z_sep(j) .le. Z_axis) then

!     Z_max(j) = Z_sep(j)

    if (j .gt. n_tht_2/2) then
      Z_max(j) = Z_sep(j) + (ZL8 - Z_xpoint(1)) * ((Z_sep(j) - Z_axis)/(Z_xpoint(1) - Z_axis))**2
    else
      Z_max(j) = Z_sep(j) + (ZL9 - Z_xpoint(1)) * ((Z_sep(j) - Z_axis)/(Z_xpoint(1) - Z_axis))**2
    endif
    
    call find_Z_surface(node_list,element_list,flux_list,i_max,Z_max(j),i_elm_find,s_find,t_find,st_find,i_find)

    call interp_RZ(node_list,element_list,i_elm_find(1),s_find(1),t_find(1),RRg1,ZZg1)


    if ( ((RRg1 .ge. R_xpoint(1)).and.(j .lt. n_tht_2/2)) .or. ((RRg1 .lt. R_xpoint(1)).and.(j .gt. n_tht_2/2)) ) then

      R_max(j)             = RRg1
      RR_new(n_psi_2,j)    = RRg1
      ZZ_new(n_psi_2,j)    = ZZg1
      ielm_flux(n_psi_2,j) = i_elm_find(1)
      s_flux(n_psi_2,j)    = s_find(1)
      t_flux(n_psi_2,j)    = t_find(1)
      t_tht(n_psi_2,j)     = 1.d0
      
    else

      call interp_RZ(node_list,element_list,i_elm_find(2),s_find(2),t_find(2),RRg1,ZZg1)
      R_max(j)             = RRg1
      RR_new(n_psi_2,j)    = RRg1
      ZZ_new(n_psi_2,j)    = ZZg1
      ielm_flux(n_psi_2,j) = i_elm_find(2)
      s_flux(n_psi_2,j)    = s_find(2)
      t_flux(n_psi_2,j)    = t_find(2)
      t_tht(n_psi_2,j)     = 1.d0

    endif

  else

    theta = theta_sep(j) !atan2(Z_sep(j)-Z_axis,R_sep(i) - R_axis)
    if (theta .lt. 0.d0) theta = theta + 2.d0 * PI

    call find_theta_surface(node_list,element_list,flux_list,i_max,theta,R_axis,Z_axis,i_elm_find,s_find,t_find,i_find)

    call interp_RZ(node_list,element_list,i_elm_find(1),s_find(1),t_find(1),RRg1,ZZg1)

    R_max(j) = RRg1
    Z_max(j) = ZZg1

    R_max(j)             = RRg1
    RR_new(n_psi_2,j)    = RRg1
    ZZ_new(n_psi_2,j)    = ZZg1
    ielm_flux(n_psi_2,j) = i_elm_find(1)
    s_flux(n_psi_2,j)    = s_find(1)
    t_flux(n_psi_2,j)    = t_find(1)
    t_tht(n_psi_2,j)     = 1.d0

  endif

enddo

!------------------------------ second part of the grid below the x-point
call tr_deallocate(s_tmp,"s_tmp",CAT_GRID); call tr_allocate(s_tmp,1,n_leg_2,"s_tmp",CAT_GRID)
s_tmp = 0
call meshac2(n_leg_2,s_tmp,0.d0,1.d0,SIG_leg_0,SIG_leg_1,0.6d0,1.0d0)

do j=1,n_leg_2

  R_min(n_tht_2 + j) = RL2 + (RL5-RL2) * s_tmp(j)

  call find_R_surface(node_list,element_list,flux_list,n_psi_2-1,R_min(n_tht_2+j),i_elm_find,s_find,t_find,st_find,i_find)

  do i=1,i_find

    call interp_RZ(node_list,element_list,i_elm_find(i),s_find(i),t_find(i),RRg1,ZZg1)

    if (ZZg1 .le. Z_xpoint(1)) exit

  enddo

  Z_min(n_tht_2 + j) = ZZg1

enddo

do j=1,n_leg_2

  R_min(n_tht_2 + n_leg_2 + j) = RL3 + (RL5-RL3) * s_tmp(j)

  call find_R_surface(node_list,element_list,flux_list,n_psi_2-1,R_min(n_tht_2+n_leg_2+j),i_elm_find,s_find,t_find,st_find,i_find)

  do i=1,i_find

    call interp_RZ(node_list,element_list,i_elm_find(i),s_find(i),t_find(i),RRg1,ZZg1)

    if (ZZg1 .le. Z_xpoint(1)) exit

  enddo

  Z_min(n_tht_2 + n_leg_2 + j) = ZZg1

enddo


do j=1,n_leg_2

  Z_max(n_tht_2 + j) = ZL1 + (Z_max(n_tht_2) -ZL1) * s_tmp(j)

  call find_Z_surface(node_list,element_list,flux_list,i_max,Z_max(n_tht_2+j),i_elm_find,s_find,t_find,st_find,i_find)

  do i=1,i_find

    call interp_RZ(node_list,element_list,i_elm_find(i),s_find(i),t_find(i),RRg1,ZZg1)

    if (RRg1 .le. R_xpoint(1)) exit

  enddo

  R_max(n_tht_2 + j) = RRg1

enddo
R_max(n_tht_2+n_leg_2) = R_max(n_tht_2)   ! this one is known

do j=1,n_leg_2

  Z_max(n_tht_2 + n_leg_2 + j) = ZL4 + (Z_max(1) - ZL4) * s_tmp(j)

  call find_Z_surface(node_list,element_list,flux_list,i_max,Z_max(n_tht_2+n_leg_2+j),i_elm_find,s_find,t_find,st_find,i_find)

  do i=1,i_find

    call interp_RZ(node_list,element_list,i_elm_find(i),s_find(i),t_find(i),RRg1,ZZg1)

    if (RRg1 .ge. R_xpoint(1)) exit

  enddo

  R_max(n_tht_2 + n_leg_2 + j) = RRg1

enddo
R_max(n_tht_2+2*n_leg_2) = R_max(1)   ! this one is known


do j=1,n_leg_2                        ! inside leg

  R_sep(n_tht_2 + j) = RL6 + (R_xpoint(1) - RL6) * s_tmp(j)
  call find_R_surface(node_list,element_list,flux_list,i_sep,R_sep(n_tht_2+j),i_elm_find,s_find,t_find,st_find,i_find)

  do i=1,i_find

    call interp_RZ(node_list,element_list,i_elm_find(i),s_find(i),t_find(i),RRg1,ZZg1)

    if (ZZg1 .le. Z_xpoint(1)) exit

  enddo

  Z_sep(n_tht_2 + j) = ZZg1

enddo
Z_sep(n_tht_2+n_leg_2) = Z_xpoint(1) ! this one is known

do j=1,n_leg_2

  R_sep(n_tht_2 + n_leg_2 + j) = RL7 + (R_xpoint(1) - RL7) * s_tmp(j)

  call find_R_surface(node_list,element_list,flux_list,i_sep,R_sep(n_tht_2+n_leg_2+j),i_elm_find,s_find,t_find,st_find,i_find)

  do i=1,i_find

    call interp_RZ(node_list,element_list,i_elm_find(i),s_find(i),t_find(i),RRg1,ZZg1)

    if (ZZg1 .le. Z_xpoint(1)) exit

  enddo

  Z_sep(n_tht_2 + n_leg_2 + j) = ZZg1

enddo
Z_sep(n_tht_2+2*n_leg_2) = Z_xpoint(1) ! this one is known

if ( write_ps ) then
  call lincol(2)
  call lplot6(1,1,R_max,Z_max,-(n_tht_2+2*n_leg_2),' ')
  call lincol(0)
  call lplot6(1,1,R_sep,Z_sep,-(n_tht_2+2*n_leg_2),' ')
endif


!------------------------------ interpolation points are known, construct polar coordinate lines
n_pieces=3
call tr_allocate(R_polar,1,n_pieces,1,4,1,n_tht_2+2*n_leg_2,"R_polar",CAT_GRID)
call tr_allocate(Z_polar,1,n_pieces,1,4,1,n_tht_2+2*n_leg_2,"Z_polar",CAT_GRID)

do j=1,n_tht_2

  delta = 0.1

  if ((j .eq. 1) .or. (j .eq. n_tht_2))       delta = 0.d0
  if ((j .eq. 2) .or. (j .eq. n_tht_2 - 1))   delta = 0.05d0

  R_polar(1,1,j) = R_axis
  R_polar(1,4,j) = delta * R_axis + (1.d0 - delta) * R_sep(j)
  R_polar(1,2,j) = ( 2.d0 * R_polar(1,1,j)  +         R_polar(1,4,j) ) / 3.d0
  R_polar(1,3,j) = (        R_polar(1,1,j)  +  2.d0 * R_polar(1,4,j) ) / 3.d0

  Z_polar(1,1,j) = Z_axis
  Z_polar(1,4,j) = delta * Z_axis + (1.d0 - delta) * Z_sep(j)
  Z_polar(1,2,j) = ( 2.d0 * Z_polar(1,1,j)  +         Z_polar(1,4,j) ) / 3.d0
  Z_polar(1,3,j) = (        Z_polar(1,1,j)  +  2.d0 * Z_polar(1,4,j) ) / 3.d0

  R_polar(3,1,j) = R_max(j)
  R_polar(3,4,j) = delta * R_max(j) + (1.d0 - delta) * R_sep(j)
  R_polar(3,2,j) = ( 2.d0 * R_polar(3,1,j)  +         R_polar(3,4,j) ) / 3.d0
  R_polar(3,3,j) = (        R_polar(3,1,j)  +  2.d0 * R_polar(3,4,j) ) / 3.d0

  Z_polar(3,1,j) = Z_max(j)
  Z_polar(3,4,j) = delta * Z_max(j) + (1.d0 - delta) * Z_sep(j)
  Z_polar(3,2,j) = ( 2.d0 * Z_polar(3,1,j)  +         Z_polar(3,4,j) ) / 3.d0
  Z_polar(3,3,j) = (        Z_polar(3,1,j)  +  2.d0 * Z_polar(3,4,j) ) / 3.d0

  R_polar(2,1,j) = R_polar(1,4,j)
  R_polar(2,4,j) = R_polar(3,4,j)
  R_polar(2,2,j) = ( R_polar(2,1,j) +  2.d0 * R_sep(j) ) / 3.d0
  R_polar(2,3,j) = ( R_polar(2,4,j) +  2.d0 * R_sep(j) ) / 3.d0

  Z_polar(2,1,j) = Z_polar(1,4,j)
  Z_polar(2,4,j) = Z_polar(3,4,j)
  Z_polar(2,2,j) = ( Z_polar(2,1,j) +  2.d0 * Z_sep(j) ) / 3.d0
  Z_polar(2,3,j) = ( Z_polar(2,4,j) +  2.d0 * Z_sep(j) ) / 3.d0

enddo

do j=1,2*n_leg_2

  delta = 0.2
  if ((j .eq. n_leg_2)   .or. (j .eq. 2*n_leg_2))    delta = 0.d0
  if ((j .eq. n_leg_2-1) .or. (j .eq. 2*n_leg_2-1))  delta = 0.05d0

  R_polar(1,1,n_tht_2+j) = R_min(n_tht_2+j)
  R_polar(1,4,n_tht_2+j) = delta * R_min(n_tht_2+j) + (1.d0 - delta) * R_sep(n_tht_2+j)
  R_polar(1,2,n_tht_2+j) = ( 2.d0 * R_polar(1,1,n_tht_2+j)  +         R_polar(1,4,n_tht_2+j) ) / 3.d0
  R_polar(1,3,n_tht_2+j) = (        R_polar(1,1,n_tht_2+j)  +  2.d0 * R_polar(1,4,n_tht_2+j) ) / 3.d0

  Z_polar(1,1,n_tht_2+j) = Z_min(n_tht_2+j)
  Z_polar(1,4,n_tht_2+j) = delta * Z_min(n_tht_2+j) + (1.d0 - delta) * Z_sep(n_tht_2+j)
  Z_polar(1,2,n_tht_2+j) = ( 2.d0 * Z_polar(1,1,n_tht_2+j)  +         Z_polar(1,4,n_tht_2+j) ) / 3.d0
  Z_polar(1,3,n_tht_2+j) = (        Z_polar(1,1,n_tht_2+j)  +  2.d0 * Z_polar(1,4,n_tht_2+j) ) / 3.d0

  R_polar(3,1,n_tht_2+j) = R_max(n_tht_2+j)
  R_polar(3,4,n_tht_2+j) = delta * R_max(n_tht_2+j) + (1.d0 - delta) * R_sep(n_tht_2+j)
  R_polar(3,2,n_tht_2+j) = ( 2.d0 * R_polar(3,1,n_tht_2+j)  +         R_polar(3,4,n_tht_2+j) ) / 3.d0
  R_polar(3,3,n_tht_2+j) = (        R_polar(3,1,n_tht_2+j)  +  2.d0 * R_polar(3,4,n_tht_2+j) ) / 3.d0

  Z_polar(3,1,n_tht_2+j) = Z_max(n_tht_2+j)
  Z_polar(3,4,n_tht_2+j) = delta * Z_max(n_tht_2+j) + (1.d0 - delta) * Z_sep(n_tht_2+j)
  Z_polar(3,2,n_tht_2+j) = ( 2.d0 * Z_polar(3,1,n_tht_2+j)  +         Z_polar(3,4,n_tht_2+j) ) / 3.d0
  Z_polar(3,3,n_tht_2+j) = (        Z_polar(3,1,n_tht_2+j)  +  2.d0 * Z_polar(3,4,n_tht_2+j) ) / 3.d0

  R_polar(2,1,n_tht_2+j) = R_polar(1,4,n_tht_2+j)
  R_polar(2,4,n_tht_2+j) = R_polar(3,4,n_tht_2+j)
  R_polar(2,2,n_tht_2+j) = ( R_polar(2,1,n_tht_2+j) +  2.d0 * R_sep(n_tht_2+j) ) / 3.d0
  R_polar(2,3,n_tht_2+j) = ( R_polar(2,4,n_tht_2+j) +  2.d0 * R_sep(n_tht_2+j) ) / 3.d0

  Z_polar(2,1,n_tht_2+j) = Z_polar(1,4,n_tht_2+j)
  Z_polar(2,4,n_tht_2+j) = Z_polar(3,4,n_tht_2+j)
  Z_polar(2,2,n_tht_2+j) = ( Z_polar(2,1,n_tht_2+j) +  2.d0 * Z_sep(n_tht_2+j) ) / 3.d0
  Z_polar(2,3,n_tht_2+j) = ( Z_polar(2,4,n_tht_2+j) +  2.d0 * Z_sep(n_tht_2+j) ) / 3.d0

enddo

if ( write_ps ) then
  call lincol(3)

  npl = 11
  call tr_allocate(xout,1,2,"xout",CAT_GRID)
  call tr_allocate(xp,1,npl,"xp",CAT_GRID)
  call tr_allocate(yp,1,npl,"yp",CAT_GRID)
  do j=1,n_tht_2+2*n_leg_2
    
    do m=1,n_pieces

      do k=1,npl
        ss = -1. + 2.*float(k-1)/float(npl-1)

        R_cub1d = (/ R_polar(m,1,j), 3.d0/2.d0 *(R_polar(m,2,j)-R_polar(m,1,j)), &
                   R_polar(m,4,j), 3.d0/2.d0 *(R_polar(m,4,j)-R_polar(m,3,j))  /)
        Z_cub1d = (/ Z_polar(m,1,j), 3.d0/2.d0 *(Z_polar(m,2,j)-Z_polar(m,1,j)), &
                   Z_polar(m,4,j), 3.d0/2.d0 *(Z_polar(m,4,j)-Z_polar(m,3,j)) /)

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
  call tr_deallocate(xout,"xout",CAT_GRID)
  call tr_deallocate(xp,"xp",CAT_GRID)
  call tr_deallocate(yp,"yp",CAT_GRID)
endif
!----------------------------------- find grid_points from crossing of coordinate lines

do j=1, n_tht_2          ! the magnetic axis

  RR_new(1,j)    = R_axis
  ZZ_new(1,j)    = Z_axis
  ielm_flux(1,j) = i_elm_axis
  s_flux(1,j)    = s_axis
  t_flux(1,j)    = t_axis
  t_tht(1,j)     = -1.d0          ! expressed in cubic Hermite (-1<t<+1)

enddo


call tr_allocate(k_cross,1,n_flux_2+n_open_2+n_private_2+1,1,n_tht_2+2*n_leg_2,"k_cross",CAT_GRID)

k_cross = 0

k_cross(1,:) = 1

do i=1,n_flux_2+n_open_2

  do j=1, n_tht_2

    do k=1,n_pieces       ! 3 line pieces per coordinate line

      R_cub1d = (/ R_polar(k,1,j), 3.d0/2.d0 *(R_polar(k,2,j)-R_polar(k,1,j)), &
                   R_polar(k,4,j), 3.d0/2.d0 *(R_polar(k,4,j)-R_polar(k,3,j))  /)
      Z_cub1d = (/ Z_polar(k,1,j), 3.d0/2.d0 *(Z_polar(k,2,j)-Z_polar(k,1,j)), &
                   Z_polar(k,4,j), 3.d0/2.d0 *(Z_polar(k,4,j)-Z_polar(k,3,j)) /)

      call find_crossing(node_list,element_list,flux_list,i,R_cub1d,Z_cub1d, &
                       RR_new(i+1,j),ZZ_new(i+1,j),ielm_flux(i+1,j),s_flux(i+1,j),t_flux(i+1,j),t_tht(i+1,j),ifail,.false.)


      if (ifail .eq. 0) then
        call interp(node_list,element_list,ielm_flux(i+1,j),1,1,s_flux(i+1,j),t_flux(i+1,j),&
                    PSg1,dPSg1_dr,dPSg1_ds,dPSg1_drs,dPSg1_drr,dPSg1_dss)

        k_cross(i+1,j) = k
        exit
      endif

    enddo

    if (ifail .ne. 0) write(*,*) ' WARNING (1) node not found (1) : ',ifail,i,j

  enddo

enddo

do i=n_flux_2,n_flux_2+n_open_2+n_private_2

  do j=n_tht_2+1, n_tht_2 + 2*n_leg_2

    do k=1,n_pieces       ! 3 line pieces per coordinate line

      R_cub1d = (/ R_polar(k,1,j), 3.d0/2.d0 *(R_polar(k,2,j)-R_polar(k,1,j)), &
                   R_polar(k,4,j), 3.d0/2.d0 *(R_polar(k,4,j)-R_polar(k,3,j))  /)
      Z_cub1d = (/ Z_polar(k,1,j), 3.d0/2.d0 *(Z_polar(k,2,j)-Z_polar(k,1,j)), &
                   Z_polar(k,4,j), 3.d0/2.d0 *(Z_polar(k,4,j)-Z_polar(k,3,j)) /)

      call find_crossing(node_list,element_list,flux_list,i,R_cub1d,Z_cub1d, &
                         RR_new(i+1,j),ZZ_new(i+1,j),ielm_flux(i+1,j),s_flux(i+1,j),t_flux(i+1,j),t_tht(i+1,j),ifail,.false.)

      if (ifail .eq. 0) then
        call interp(node_list,element_list,ielm_flux(i+1,j),1,1,s_flux(i+1,j),t_flux(i+1,j),&
                    PSg1,dPSg1_dr,dPSg1_ds,dPSg1_drs,dPSg1_drr,dPSg1_dss)

        k_cross(i+1,j) = k
        exit
      endif

     enddo

     if (ifail .ne. 0) write(*,*) ' WARNING node not found (2) : ',ifail,i,j

  enddo

enddo

!***********************************************************************
!*     define the new nodes and finite elements    (nodes first)       *
!***********************************************************************

! --- Allocate data structures for new nodes and elements and initialize them.
allocate(newnode_list)
call tr_register_mem(sizeof(newnode_list),"newnode_list",CAT_GRID)
allocate(newelement_list)
call tr_register_mem(sizeof(newelement_list),"newelement_list",CAT_GRID)
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

do i=1,n_flux-1                 !------------------------ the closed field lines
  do j=1, n_tht-1

    i2 = 2*(i-1) + 1
    j2 = 2*(j-1) + 1

    k = k_cross(i2,j2)

    R_cub1d = (/ R_polar(k,1,j2), 3.d0/2.d0 *(R_polar(k,2,j2)-R_polar(k,1,j2)), &
                 R_polar(k,4,j2), 3.d0/2.d0 *(R_polar(k,4,j2)-R_polar(k,3,j2))  /)
    Z_cub1d = (/ Z_polar(k,1,j2), 3.d0/2.d0 *(Z_polar(k,2,j2)-Z_polar(k,1,j2)), &
                 Z_polar(k,4,j2), 3.d0/2.d0 *(Z_polar(k,4,j2)-Z_polar(k,3,j2)) /)

    call CUB1D(R_cub1d(1), R_cub1d(2), R_cub1d(3), R_cub1d(4),t_tht(i2,j2),tmp1, dR_dt)
    call CUB1D(Z_cub1d(1), Z_cub1d(2), Z_cub1d(3), Z_cub1d(4),t_tht(i2,j2),tmp2, dZ_dt)

    call interp_RZ(node_list,element_list,ielm_flux(i2,j2),s_flux(i2,j2),t_flux(i2,j2), &
                   RRg1,dRRg1_dr,dRRg1_ds,ZZg1,dZZg1_dr,dZZg1_ds)

    call interp(node_list,element_list,ielm_flux(i2,j2),1,1,s_flux(i2,j2),t_flux(i2,j2),&
                   PSg1,dPSg1_dr,dPSg1_ds,dPSg1_drs,dPSg1_drr,dPSg1_dss)

    RZ_jac  = DRRg1_dr * dZZg1_ds - dRRg1_ds * dZZg1_dr

    PSI_R = (   dPSg1_dr * dZZg1_ds - dPSg1_ds * dZZg1_dr ) / RZ_jac
    PSI_Z = ( - dPSg1_dr * dRRg1_ds + dPSg1_ds * dRRg1_dr ) / RZ_jac


    node  = (n_tht-1)*(i-1) + j
    index = node

    newnode_list%node(index)%x(1,1,:) = (/ RR_new(i2,j2), ZZ_new(i2,j2) /)
    newnode_list%node(index)%x(1,2,:) = (/ dR_dt, dZ_dt /) / sqrt(dR_dt**2 + dZ_dt**2)
    newnode_list%node(index)%boundary = 0

    if (i .eq. 1) then   !------------------------------------ magnetic axis : special case
      newnode_list%node(index)%x(1,3,:) = 0.d0
      newnode_list%node(index)%x(1,4,:) = 0.d0
    else
      newnode_list%node(index)%x(1,3,:) = (/ -PSI_Z, +PSI_R /) / sqrt(PSI_R**2 + PSI_Z**2)
      newnode_list%node(index)%x(1,4,:) = 0.d0
    endif

  enddo
enddo
newnode_list%n_nodes = node

!----------------------------------------- add multiple nodes at the x-point
index_xpoint = newnode_list%n_nodes + 1

n_xpoint=4

call interp(node_list,element_list,i_elm_xpoint(1),1,1,s_xpoint(1),t_xpoint(1),&
                   PSg1,dPSg1_dr,dPSg1_ds,dPSg1_drs,dPSg1_drr,dPSg1_dss)
call interp_RZ(node_list,element_list,i_elm_xpoint(1),s_xpoint(1),t_xpoint(1), &
                   RRg1,dRRg1_dr,dRRg1_ds,ZZg1,dZZg1_dr,dZZg1_ds)

EJAC = dRRg1_dr*dZZg1_ds - dRRg1_ds * dZZg1_dr
RY = - dRRg1_ds / EJAC
RX =   dZZg1_ds / EJAC
SY =   dRRg1_dr / EJAC
SX = - dZZg1_dr / EJAC

CRR = (dPSg1_drr*RX*RX + 2.*dPSg1_drs*RX*SX + dPSg1_dss*SX*SX)/2.
CZZ = (dPSg1_drr*RY*RY + 2.*dPSg1_drs*RY*SY + dPSg1_dss*SY*SY)/2.
CRZ = (dPSg1_drr*RX*RY + dPSg1_drs*(RX*SY + RY*SX) + dPSg1_dss*SX*SY)/2.

write(*,'(A,3e16.8)') ' CRR, CZZ, CRZ : ',CRR,CZZ,CRZ
write(*,*) ' xpoint node_index : ',index_xpoint

alpha1 = root(CRR,CRZ,CZZ,CRZ*CRZ-4.d0*CRR*CZZ,+1.d0)
alpha2 = root(CRR,CRZ,CZZ,CRZ*CRZ-4.d0*CRR*CZZ,-1.d0)

alpha_max = max(alpha1,alpha2)
alpha_min = min(alpha1,alpha2)

write(*,*) ' ALPHA : ',alpha_max,alpha_min

newnode_list%node(index_xpoint)%x(1,1,:) = (/ R_xpoint(1), Z_xpoint(1) /)
newnode_list%node(index_xpoint)%x(1,2,:) = (/ cos(angle_L9),sin(angle_L9) /)
newnode_list%node(index_xpoint)%x(1,3,:) = (/ alpha_max,1.d0 /) / sqrt(alpha_max**2 + 1.d0)
newnode_list%node(index_xpoint)%x(1,4,:) = 0.d0
newnode_list%node(index_xpoint)%boundary = 0

newnode_list%node(index_xpoint+1)%x(1,1,:) = (/ R_xpoint(1), Z_xpoint(1) /)
newnode_list%node(index_xpoint+1)%x(1,2,:) = (/ R_xpoint(1)-R_axis,Z_xpoint(1)-Z_axis/)/sqrt((R_xpoint(1)-R_axis)**2+(Z_xpoint(1)-Z_axis)**2)
newnode_list%node(index_xpoint+1)%x(1,3,:) = (/ alpha_max,1.d0 /) / sqrt(alpha_max**2 + 1.d0)
newnode_list%node(index_xpoint+1)%x(1,4,:) = 0.d0
newnode_list%node(index_xpoint+1)%boundary = 0

newnode_list%node(index_xpoint+2)%x(1,1,:) = (/ R_xpoint(1), Z_xpoint(1) /)
newnode_list%node(index_xpoint+2)%x(1,2,:) = (/ R_xpoint(1)-R_axis,Z_xpoint(1)-Z_axis/)/sqrt((R_xpoint(1)-R_axis)**2+(Z_xpoint(1)-Z_axis)**2)
newnode_list%node(index_xpoint+2)%x(1,3,:) = (/ alpha_min,1.d0 /) / sqrt(alpha_min**2 + 1.d0)
newnode_list%node(index_xpoint+2)%x(1,4,:) = 0.d0
newnode_list%node(index_xpoint+2)%boundary = 0

newnode_list%node(index_xpoint+3)%x(1,1,:) = (/ R_xpoint(1), Z_xpoint(1) /)
newnode_list%node(index_xpoint+3)%x(1,2,:) = (/ -cos(angle_L8),-sin(angle_L8) /)
newnode_list%node(index_xpoint+3)%x(1,3,:) = (/ alpha_min,1.d0 /) / sqrt(alpha_min**2 + 1.d0)
newnode_list%node(index_xpoint+3)%x(1,4,:) = 0.d0
newnode_list%node(index_xpoint+3)%boundary = 0

newnode_list%n_nodes = newnode_list%n_nodes + n_xpoint

index = newnode_list%n_nodes

do i=n_flux,n_flux+n_open           !--------------------------- nodes on the open field lines

  j_start = 1; j_end = n_tht   ! skip first and last point on separatrix (x-points already added)
  if (i.eq. n_flux) then
    j_start = 2
    j_end   = n_tht-1
  endif

  do j=j_start, j_end

    i2 = 2*(i-1) + 1
    j2 = 2*(j-1) + 1

    k = k_cross(i2,j2)

    R_cub1d = (/ R_polar(k,1,j2), 3.d0/2.d0 *(R_polar(k,2,j2)-R_polar(k,1,j2)), &
                 R_polar(k,4,j2), 3.d0/2.d0 *(R_polar(k,4,j2)-R_polar(k,3,j2))  /)
    Z_cub1d = (/ Z_polar(k,1,j2), 3.d0/2.d0 *(Z_polar(k,2,j2)-Z_polar(k,1,j2)), &
                 Z_polar(k,4,j2), 3.d0/2.d0 *(Z_polar(k,4,j2)-Z_polar(k,3,j2)) /)

    if ((j .eq. 1) .or. (j .eq. n_tht)) then
      R_cub1d = (/ R_xpoint(1), 0.5d0*(R_max(j2)-R_xpoint(1)), R_max(j2), 0.5d0*(R_max(j2)-R_xpoint(1))  /)
      Z_cub1d = (/ Z_xpoint(1), 0.5d0*(Z_max(j2)-Z_xpoint(1)), Z_max(j2), 0.5d0*(Z_max(j2)-Z_xpoint(1))  /)
    endif

    call CUB1D(R_cub1d(1), R_cub1d(2), R_cub1d(3), R_cub1d(4),t_tht(i2,j2),tmp1, dR_dt)
    call CUB1D(Z_cub1d(1), Z_cub1d(2), Z_cub1d(3), Z_cub1d(4),t_tht(i2,j2),tmp2, dZ_dt)

    call interp_RZ(node_list,element_list,ielm_flux(i2,j2),s_flux(i2,j2),t_flux(i2,j2), &
                   RRg1,dRRg1_dr,dRRg1_ds,ZZg1,dZZg1_dr,dZZg1_ds)

    call interp(node_list,element_list,ielm_flux(i2,j2),1,1,s_flux(i2,j2),t_flux(i2,j2),&
                   PSg1,dPSg1_dr,dPSg1_ds,dPSg1_drs,dPSg1_drr,dPSg1_dss)

    RZ_jac  = DRRg1_dr * dZZg1_ds - dRRg1_ds * dZZg1_dr

    PSI_R = (   dPSg1_dr * dZZg1_ds - dPSg1_ds * dZZg1_dr ) / RZ_jac
    PSI_Z = ( - dPSg1_dr * dRRg1_ds + dPSg1_ds * dRRg1_dr ) / RZ_jac

    index = index + 1

    newnode_list%node(index)%x(1,1,:) = (/ RR_new(i2,j2), ZZ_new(i2,j2) /)
    newnode_list%node(index)%x(1,2,:) = (/ dR_dt, dZ_dt /)   / sqrt(dR_dt**2 + dZ_dt**2)
    newnode_list%node(index)%x(1,3,:) = (/ -PSI_Z, +PSI_R /) / sqrt(PSI_R**2 + PSI_Z**2)
    newnode_list%node(index)%x(1,4,:) = 0.d0
    newnode_list%node(index)%boundary = 0

    if (i .eq. n_flux+n_open)                   newnode_list%node(index)%boundary = 2

  enddo
enddo

newnode_list%n_nodes = newnode_list%n_nodes + (n_open+1) * n_tht - 2

index = newnode_list%n_nodes

do j=1, n_leg                         !--------------------------- nodes on right leg

  do k=1,n_open+n_private+1           !--------------------------- nodes on right leg


    if (k .le. n_private) then
      i = n_flux + n_open + n_private - k + 1
    else
      i = n_flux + (k-n_private) - 1
    endif

    if (.not. ( (j .eq. n_leg) .and. (i .le. n_flux+n_open))) then  ! exclude horizontal line (already there)

      i2 = 2*(i-1) + 1
      j2 = n_tht_2 + n_leg_2 + 2*(j-1) + 1

      m = k_cross(i2,j2)

      R_cub1d = (/ R_polar(m,1,j2), 3.d0/2.d0 *(R_polar(m,2,j2)-R_polar(m,1,j2)), &
                   R_polar(m,4,j2), 3.d0/2.d0 *(R_polar(m,4,j2)-R_polar(m,3,j2))  /)
      Z_cub1d = (/ Z_polar(m,1,j2), 3.d0/2.d0 *(Z_polar(m,2,j2)-Z_polar(m,1,j2)), &
                   Z_polar(m,4,j2), 3.d0/2.d0 *(Z_polar(m,4,j2)-Z_polar(m,3,j2)) /)

      call CUB1D(R_cub1d(1), R_cub1d(2), R_cub1d(3), R_cub1d(4),t_tht(i2,j2),tmp1, dR_dt)
      call CUB1D(Z_cub1d(1), Z_cub1d(2), Z_cub1d(3), Z_cub1d(4),t_tht(i2,j2),tmp2, dZ_dt)

      call interp_RZ(node_list,element_list,ielm_flux(i2,j2),s_flux(i2,j2),t_flux(i2,j2), &
                     RRg1,dRRg1_dr,dRRg1_ds,ZZg1,dZZg1_dr,dZZg1_ds)

      call interp(node_list,element_list,ielm_flux(i2,j2),1,1,s_flux(i2,j2),t_flux(i2,j2),&
                     PSg1,dPSg1_dr,dPSg1_ds,dPSg1_drs,dPSg1_drr,dPSg1_dss)

      RZ_jac  = DRRg1_dr * dZZg1_ds - dRRg1_ds * dZZg1_dr

      PSI_R = (   dPSg1_dr * dZZg1_ds - dPSg1_ds * dZZg1_dr ) / RZ_jac
      PSI_Z = ( - dPSg1_dr * dRRg1_ds + dPSg1_ds * dRRg1_dr ) / RZ_jac

      index = index + 1

      newnode_list%node(index)%x(1,1,:) = (/ RR_new(i2,j2), ZZ_new(i2,j2) /)
      newnode_list%node(index)%x(1,2,:) = (/ dR_dt, dZ_dt /)   / sqrt(dR_dt**2 + dZ_dt**2)
      newnode_list%node(index)%x(1,3,:) = (/ -PSI_Z, +PSI_R /) / sqrt(PSI_R**2 + PSI_Z**2)
      newnode_list%node(index)%x(1,4,:) = 0.d0
      newnode_list%node(index)%boundary = 0

      if ((k.eq.1) .or. (k .eq. n_open+n_private+1))  newnode_list%node(index)%boundary = newnode_list%node(index)%boundary + 2
      if  (j .eq. 1)                                  newnode_list%node(index)%boundary = newnode_list%node(index)%boundary + 1

   endif

  enddo
enddo

do l=1, n_leg-1                       !--------------------------- nodes on left leg

  j = n_leg - l

  do k=1,n_open+n_private+1           !--------------------------- nodes on left leg

    if (k .le. n_private) then
      i = n_flux + n_open + n_private - k + 1
    else
      i = n_flux + (k-n_private) - 1
    endif

    i2 = 2*(i-1) + 1
    j2 = n_tht_2 + 2*(j-1) + 1

    m = k_cross(i2,j2)

    R_cub1d = (/ R_polar(m,1,j2), 3.d0/2.d0 *(R_polar(m,2,j2)-R_polar(m,1,j2)), &
                 R_polar(m,4,j2), 3.d0/2.d0 *(R_polar(m,4,j2)-R_polar(m,3,j2))  /)
    Z_cub1d = (/ Z_polar(m,1,j2), 3.d0/2.d0 *(Z_polar(m,2,j2)-Z_polar(m,1,j2)), &
                 Z_polar(m,4,j2), 3.d0/2.d0 *(Z_polar(m,4,j2)-Z_polar(m,3,j2)) /)

    call CUB1D(R_cub1d(1), R_cub1d(2), R_cub1d(3), R_cub1d(4),t_tht(i2,j2),tmp1, dR_dt)
    call CUB1D(Z_cub1d(1), Z_cub1d(2), Z_cub1d(3), Z_cub1d(4),t_tht(i2,j2),tmp2, dZ_dt)

    call interp_RZ(node_list,element_list,ielm_flux(i2,j2),s_flux(i2,j2),t_flux(i2,j2), &
                   RRg1,dRRg1_dr,dRRg1_ds,ZZg1,dZZg1_dr,dZZg1_ds)

    call interp(node_list,element_list,ielm_flux(i2,j2),1,1,s_flux(i2,j2),t_flux(i2,j2),&
                   PSg1,dPSg1_dr,dPSg1_ds,dPSg1_drs,dPSg1_drr,dPSg1_dss)

    RZ_jac  = DRRg1_dr * dZZg1_ds - dRRg1_ds * dZZg1_dr

    PSI_R = (   dPSg1_dr * dZZg1_ds - dPSg1_ds * dZZg1_dr ) / RZ_jac
    PSI_Z = ( - dPSg1_dr * dRRg1_ds + dPSg1_ds * dRRg1_dr ) / RZ_jac

    index = index + 1

    newnode_list%node(index)%x(1,1,:) = (/ RR_new(i2,j2), ZZ_new(i2,j2) /)
    newnode_list%node(index)%x(1,2,:) = (/ dR_dt, dZ_dt /)   / sqrt(dR_dt**2 + dZ_dt**2)
    newnode_list%node(index)%x(1,3,:) = (/ -PSI_Z, +PSI_R /) / sqrt(PSI_R**2 + PSI_Z**2)
    newnode_list%node(index)%x(1,4,:) = 0.d0
    newnode_list%node(index)%boundary = 0


    if ((k.eq.1) .or. (k .eq. n_open+n_private+1))  newnode_list%node(index)%boundary = newnode_list%node(index)%boundary + 2
    if  (j .eq. 1)                                  newnode_list%node(index)%boundary = newnode_list%node(index)%boundary + 1

 enddo
enddo

newnode_list%n_nodes = index

write(*,*) ' definition of nodes completed ',newnode_list%n_nodes

if ( write_ps ) call nframe(11,11,1,2.0,3.0,-2.0,-1.0,' ',1,'R',1,'Z',1)
call plot_flux_surfaces(node_list,element_list,flux_list,.false.,1,.true.,xcase)


call tr_allocate(xp,1,index,"xp",CAT_GRID)
call tr_allocate(yp,1,index,"yp",CAT_GRID)

do i=1,newnode_list%n_nodes
  xp(i) = newnode_list%node(i)%x(1,1,1)
  yp(i) = newnode_list%node(i)%x(1,1,2)
!  write(label,'(i4)') i
!  call dlch(int(90.+900.*(xp(i)-2.5)/1.),int(77.+645.*(yp(i)+1.71)/0.71),label,4,1)
enddo

if ( write_ps ) call lplot(1,1,421,xp,yp,-newnode_list%n_nodes,1,'R',1,'Z',1,'nodes',5)
call tr_deallocate(xp,"xp",CAT_GRID)
call tr_deallocate(yp,"yp",CAT_GRID)

!***********************************************************************
!*                   define the new elements                           *
!***********************************************************************

call tr_allocate(keep,1,n_psi*n_tht*2,1,4,1,2,"keep",CAT_GRID)

do i=1,n_flux-1

  do j=1, n_tht-1

    index = (n_tht-1)*(i-1) + j

    newelement_list%element(index)%vertex(1) = (i-1)*(n_tht-1) + j
    newelement_list%element(index)%vertex(2) = (i  )*(n_tht-1) + j
    newelement_list%element(index)%vertex(3) = (i  )*(n_tht-1) + j + 1
    newelement_list%element(index)%vertex(4) = (i-1)*(n_tht-1) + j + 1

    newelement_list%element(index)%size = 1.d0

    keep(index,1,:) = (/ 2*i  , 2*j-1 /)             ! node numbers of the mid-points in the original 'regular' grid
    keep(index,2,:) = (/ 2*i+1, 2*j   /)
    keep(index,3,:) = (/ 2*i  , 2*j+1 /)
    keep(index,4,:) = (/ 2*i-1, 2*j   /)

    if (j .eq. n_tht-1) then
      newelement_list%element(index)%vertex(4)  = (i-1)*(n_tht-1) + 1
      newelement_list%element(index)%vertex(3)  = (i  )*(n_tht-1) + 1
    endif

    if (i .eq. n_flux-1) then
      newelement_list%element(index)%vertex(3) = (i  )*(n_tht-1) + j + n_xpoint
      newelement_list%element(index)%vertex(2) = (i  )*(n_tht-1) + j + n_xpoint - 1

      if (j.eq.1) then
        newelement_list%element(index)%vertex(2)  = index_xpoint + 1
      endif

      if (j.eq.n_tht-1) then
        newelement_list%element(index)%vertex(3)  = index_xpoint + 2
      endif

    endif

  enddo

enddo

newelement_list%n_elements = (n_flux - 1)*(n_tht - 1)

if ( newnode_list%n_nodes > n_nodes_max ) then
  write(*,*) 'ERROR in grid_xpoint: hard-coded parameter n_nodes_max is too small'
  stop
else if ( newelement_list%n_elements > n_elements_max ) then
  write(*,*) 'ERROR in grid_xpoint: hard-coded parameter n_elements_max is too small'
  stop
end if

node_start = (n_flux-1) * (n_tht-1) + n_xpoint

do i=1,n_open

  do j=1, n_tht-1

    index = newelement_list%n_elements + (n_tht-1)*(i-1) + j

    if ( i .eq. 1) then

      newelement_list%element(index)%vertex(1) = node_start + j - 1
      newelement_list%element(index)%vertex(2) = node_start + (n_tht-1) + j - 1
      newelement_list%element(index)%vertex(3) = node_start + (n_tht-1) + j
      newelement_list%element(index)%vertex(4) = node_start + j

      if (j.eq.1)       newelement_list%element(index)%vertex(1) = node_start - n_xpoint + 1
      if (j.eq.n_tht-1) newelement_list%element(index)%vertex(4) = node_start

    else

      newelement_list%element(index)%vertex(1) = node_start - 1 + (i-1)*n_tht + j - 1
      newelement_list%element(index)%vertex(2) = node_start - 1 + (i  )*n_tht + j - 1
      newelement_list%element(index)%vertex(3) = node_start - 1 + (i  )*n_tht + j
      newelement_list%element(index)%vertex(4) = node_start - 1 + (i-1)*n_tht + j

    endif

    keep(index,1,:) = (/ n_flux_2 + 2*i,   2*j-1 /)             ! node numbers of the mid-points in the original 'regular' grid
    keep(index,2,:) = (/ n_flux_2 + 2*i+1, 2*j   /)
    keep(index,3,:) = (/ n_flux_2 + 2*i,   2*j+1 /)
    keep(index,4,:) = (/ n_flux_2 + 2*i-1, 2*j   /)

    newelement_list%element(index)%size = 1.d0

  enddo

enddo

newelement_list%n_elements = newelement_list%n_elements + (n_open) * (n_tht-1)

if ( newnode_list%n_nodes > n_nodes_max ) then
  write(*,*) 'ERROR in grid_xpoint: hard-coded parameter n_nodes_max is too small'
  stop
else if ( newelement_list%n_elements > n_elements_max ) then
  write(*,*) 'ERROR in grid_xpoint: hard-coded parameter n_elements_max is too small'
  stop
end if

node_start   = (n_flux-1)*(n_tht-1) + 4 + (n_open+1) * n_tht - 2

index = newelement_list%n_elements

do i=1, n_open+n_private

  do j=1, 2*n_leg - 2

    index = index + 1

    if (j .lt. n_leg) then

      newelement_list%element(index)%vertex(1) = node_start + (j-1)*(n_open+n_private+1) + i
      newelement_list%element(index)%vertex(2) = node_start + (j-1)*(n_open+n_private+1) + i+1
      newelement_list%element(index)%vertex(3) = node_start + (j  )*(n_open+n_private+1) + i+1
      newelement_list%element(index)%vertex(4) = node_start + (j  )*(n_open+n_private+1) + i

      if (i .lt. n_private) then
        keep(index,1,:) = (/ n_flux_2 + n_open_2 + n_private_2 + 1 - 2*(i-1)-1, n_tht_2 + n_leg_2 + 2*(j-1)+1 /)
        keep(index,2,:) = (/ n_flux_2 + n_open_2 + n_private_2 + 1 - 2*(i  )  , n_tht_2 + n_leg_2 + 2*(j-1)+2 /)
        keep(index,3,:) = (/ n_flux_2 + n_open_2 + n_private_2 + 1 - 2*(i  )+1, n_tht_2 + n_leg_2 + 2*(j  )+1 /)
        keep(index,4,:) = (/ n_flux_2 + n_open_2 + n_private_2 + 1 - 2*(i-1)  , n_tht_2 + n_leg_2 + 2*(j  )   /)
      elseif (i .eq. n_private) then
        keep(index,1,:) = (/ n_flux_2 + n_open_2 + n_private_2 + 1 - 2*(i-1)-1, n_tht_2 + n_leg_2 + 2*(j-1)+1 /)
        keep(index,2,:) = (/ n_flux_2 + 1                                     , n_tht_2 + n_leg_2 + 2*(j-1)+2 /)
        keep(index,3,:) = (/ n_flux_2 + n_open_2 + n_private_2 + 1 - 2*(i  )+1, n_tht_2 + n_leg_2 + 2*(j  )+1 /)
        keep(index,4,:) = (/ n_flux_2 + n_open_2 + n_private_2 + 1 - 2*(i-1)  , n_tht_2 + n_leg_2 + 2*(j  )   /)
      else
        keep(index,1,:) = (/ n_flux_2 + 2*(i-n_private)  , n_tht_2 + n_leg_2 + 2*(j-1)+1 /)
        keep(index,2,:) = (/ n_flux_2 + 2*(i-n_private)+1, n_tht_2 + n_leg_2 + 2*(j-1)+2 /)
        keep(index,3,:) = (/ n_flux_2 + 2*(i-n_private),   n_tht_2 + n_leg_2 + 2*(j  )+1 /)
        keep(index,4,:) = (/ n_flux_2 + 2*(i-n_private)-1, n_tht_2 + n_leg_2 + 2*(j  )   /)
      endif

    elseif (j .ge. n_leg) then

      newelement_list%element(index)%vertex(1) = node_start + (j-1)*(n_open+n_private+1) + i   - n_open - 1
      newelement_list%element(index)%vertex(2) = node_start + (j-1)*(n_open+n_private+1) + i+1 - n_open - 1
      newelement_list%element(index)%vertex(3) = node_start + (j  )*(n_open+n_private+1) + i+1 - n_open - 1
      newelement_list%element(index)%vertex(4) = node_start + (j  )*(n_open+n_private+1) + i   - n_open - 1

      if (i .lt. n_private) then
        keep(index,1,:) = (/ n_flux_2 + n_open_2 + n_private_2 + 1 - 2*(i-1)-1, n_tht_2 + n_leg_2 - 2*(j-n_leg)    /)
        keep(index,2,:) = (/ n_flux_2 + n_open_2 + n_private_2 + 1 - 2*(i  )  , n_tht_2 + n_leg_2 - 2*(j-n_leg) - 1 /)
        keep(index,3,:) = (/ n_flux_2 + n_open_2 + n_private_2 + 1 - 2*(i  )+1, n_tht_2 + n_leg_2 - 2*(j-n_leg) - 2 /)
        keep(index,4,:) = (/ n_flux_2 + n_open_2 + n_private_2 + 1 - 2*(i-1)  , n_tht_2 + n_leg_2 - 2*(j-n_leg) - 1 /)
      elseif (i .eq. n_private) then
        keep(index,1,:) = (/ n_flux_2 + n_open_2 + n_private_2 + 1 - 2*(i-1)-1, n_tht_2 + n_leg_2 - 2*(j-n_leg)    /)
        keep(index,2,:) = (/ n_flux_2 + 1                                     , n_tht_2 + n_leg_2 - 2*(j-n_leg) - 1 /)
        keep(index,3,:) = (/ n_flux_2 + n_open_2 + n_private_2 + 1 - 2*(i  )+1, n_tht_2 + n_leg_2 - 2*(j-n_leg) - 2 /)
        keep(index,4,:) = (/ n_flux_2 + n_open_2 + n_private_2 + 1 - 2*(i-1)  , n_tht_2 + n_leg_2 - 2*(j-n_leg) - 1 /)
      else
        keep(index,1,:) = (/ n_flux_2 + 2*(i-n_private)  , n_tht_2 + n_leg_2 - 2*(j-n_leg)   /)
        keep(index,2,:) = (/ n_flux_2 + 2*(i-n_private)+1, n_tht_2 + n_leg_2 - 2*(j-n_leg)-1 /)
        keep(index,3,:) = (/ n_flux_2 + 2*(i-n_private)  , n_tht_2 + n_leg_2 - 2*(j-n_leg)-2 /)
        keep(index,4,:) = (/ n_flux_2 + 2*(i-n_private)-1, n_tht_2 + n_leg_2 - 2*(j-n_leg)-1 /)
      endif

    endif

    if (j .eq. n_leg - 1) then

      if (i .gt. n_private) then
        newelement_list%element(index)%vertex(3) = node_start - n_open*n_tht + 1 + (i-n_private-1)*n_tht
        newelement_list%element(index)%vertex(4) = node_start - n_open*n_tht + 1 + (i-n_private-2)*n_tht
      endif
      if (i .eq. n_private) then
         newelement_list%element(index)%vertex(3) = index_xpoint + 2
      endif
      if (i .eq. n_private+1) then
         newelement_list%element(index)%vertex(4) = index_xpoint + 3
      endif
    endif

    if (j .eq. n_leg ) then

      if (i .gt. n_private) then
        newelement_list%element(index)%vertex(1) = node_start + (i-n_private-n_open-1)*n_tht
        newelement_list%element(index)%vertex(2) = node_start + (i-n_private-n_open  )*n_tht
      else
        newelement_list%element(index)%vertex(1) = node_start + (j-1)*(n_open+n_private+1) + i
        newelement_list%element(index)%vertex(2) = node_start + (j-1)*(n_open+n_private+1) + i+1
      endif
      if (i.eq. n_private)   then
        newelement_list%element(index)%vertex(2) = index_xpoint + 1
      endif
      if (i.eq. n_private+1) then
        newelement_list%element(index)%vertex(1) = index_xpoint
      endif
    endif

  enddo

enddo

newelement_list%n_elements = newelement_list%n_elements + (n_open+n_private)*(2*n_leg-2)

if ( newnode_list%n_nodes > n_nodes_max ) then
  write(*,*) 'ERROR in grid_xpoint: hard-coded parameter n_nodes_max is too small'
  stop
else if ( newelement_list%n_elements > n_elements_max ) then
  write(*,*) 'ERROR in grid_xpoint: hard-coded parameter n_elements_max is too small'
  stop
end if

!if ( write_ps ) then
!  allocate(xp(index),yp(index))

!  call lincol(0)
!  do i=1,newelement_list%n_elements

!    xp(i) = 0.
!    yp(i) = 0.

!    do k=1,4
!      iv= newelement_list%element(i)%vertex(k)

!      xp(i) = xp(i) +  newnode_list%node(iv)%x(1,1,1) /4.
!      yp(i) = yp(i) +  newnode_list%node(iv)%x(1,1,2) /4.
!    enddo

!    write(label,'(i4)') i

!    call dlch(int(90.+900.*(xp(i)-2.5)/1.),int(77.+645.*(yp(i)+1.71)/0.71),label,4,1)

!  enddo

!  deallocate(xp,yp)
!endif

!***********************************************************************
!*             adjust size of elements to get better match             *
!***********************************************************************


do k=1, newelement_list%n_elements   ! fill in the size of the elements

  do iv = 1, 4                    ! over 4 sides of an element

    ivp = mod(iv,4)   + 1         ! vertex with index one higher

    i2 = keep(k,iv,1)             ! index of the mid-point nodes in the original regular grids
    j2 = keep(k,iv,2)

    node_iv  = newelement_list%element(k)%vertex(iv)
    node_ivp = newelement_list%element(k)%vertex(ivp) 

!    write(*,'(A,8i5)') ' midpoint node : ',k,iv,i2,j2

    Rmid = RR_new(i2,j2)
    Zmid = ZZ_new(i2,j2)

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

!    if (abs(denom) .gt. 1.d+8) then

!      size_0 = + 4.d0 / 3.d0 * (dZp * ( R0 + Rp - 2.d0*Rmid) - dRp * ( Z0 + ZP - 2.d0*Zmid) )/ ( dRP * dZ0 - dR0 * dZP)
!      size_P = - 4.d0 / 3.d0 * (dZ0 * ( R0 + Rp - 2.d0*Rmid) - dR0 * ( Z0 + ZP - 2.d0*Zmid) )/ ( dRP * dZ0 - dR0 * dZP)

!    else

      size_0 = sign(sqrt((R0-RP)**2 + (Z0-ZP)**2) /3.d0, dR0 * (RP-R0) + dZ0 * (ZP-Z0) )
      size_P = sign(sqrt((R0-RP)**2 + (Z0-ZP)**2) /3.d0, dRP * (R0-RP) + dZP * (Z0-ZP) )

      if ((R0-RP)**2 + (Z0-ZP)**2 .eq. 0.d0) then
        size_0 = 1.d0
        size_P = 1.d0
      endif

!    endif

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

  newelement_list%element(k)%father     = 0
  newelement_list%element(k)%n_sons     = 0
  element_list%element(Index)%sons(:)   = 0
enddo

!call plot_flux_surfaces(node_list,element_list,flux_list,.true.,1)

!***********************************************************************
!*             fill in the values into the new grid                    *
!***********************************************************************

index = 0
do i=1,newnode_list%n_nodes

  newnode_list%node(i)%axis_node = .false.
  newnode_list%node(i)%axis_dof  = 0

  if (i .lt. n_tht) then
    newnode_list%node(i)%axis_node = .true.
    newnode_list%node(i)%axis_dof  = 2
  endif

  do k=1,n_order+1

    index = index + 1
    newnode_list%node(i)%index(k) = index

    if ((force_central_node) .and. (i .gt. 1) .and. (i .lt. n_tht) .and. (k.eq.1)) then
      newnode_list%node(i)%index(k) = newnode_list%node(1)%index(1)
      index = index - 1
    endif
    ! Share 4 degrees of freedom for all nodes on the grid axis.    
    if ((treat_axis) .and. (i .gt. 1) .and. (i .lt. n_tht) .and. (k.le.n_order+1)) then
      newnode_list%node(i)%index(k) = newnode_list%node(1)%index(k)
      index = index - 1
    endif    
    if ((i .eq. index_xpoint+1).and.(k.eq.1)) then
      newnode_list%node(i)%index(k) = newnode_list%node(index_xpoint)%index(k)
      index = index - 1
    endif
    if ((i .eq. index_xpoint+1).and.(k.eq.3)) then
      newnode_list%node(i)%index(k) = newnode_list%node(index_xpoint)%index(k)
      index = index - 1
    endif
    if ((i .eq. index_xpoint+2).and.(k.eq.1)) then
      newnode_list%node(i)%index(k) = newnode_list%node(index_xpoint)%index(k)
      index = index - 1
    endif
    if ((i .eq. index_xpoint+2).and.(k.eq.2)) then
      newnode_list%node(i)%index(k) = newnode_list%node(index_xpoint+1)%index(k)
      index = index - 1
    endif
    if ((i .eq. index_xpoint+3).and.(k.eq.1)) then
      newnode_list%node(i)%index(k) = newnode_list%node(index_xpoint)%index(k)
      index = index - 1
    endif
    if ((i .eq. index_xpoint+3).and.(k.eq.2)) then
      newnode_list%node(i)%index(k) = newnode_list%node(index_xpoint)%index(k)
      index = index - 1
    endif
    if ((i .eq. index_xpoint+3).and.(k.eq.3)) then
      newnode_list%node(i)%index(k) = newnode_list%node(index_xpoint+2)%index(k)
      index = index - 1
    endif
  enddo
  
  newnode_list%node(i)%constrained = .false.
enddo

if (fix_axis_nodes) then
  do k=1, newelement_list%n_elements
    do iv=1,4
      j = newelement_list%element(k)%vertex(iv)
      if (newnode_list%node(j)%axis_node) then
        newelement_list%element(k)%size(iv,3) = 0.d0
        newelement_list%element(k)%size(iv,4) = 0.d0
      endif
    enddo
  enddo
endif

do i=1,newnode_list%n_nodes

  R1 = newnode_list%node(i)%x(1,1,1)
  Z1 = newnode_list%node(i)%x(1,1,2)

  call find_RZ(node_list,element_list,R1,Z1,R_out,Z_out,ielm_out,s_out,t_out,ifail)

  call interp_RZ(node_list,element_list,ielm_out,s_out,t_out, &
                 RRg1,dRRg1_dr,dRRg1_ds,ZZg1,dZZg1_dr,dZZg1_ds)

  call interp(node_list,element_list,ielm_out,1,1,s_out,t_out,PSg1,dPSg1_dr,dPSg1_ds,dPSg1_drs,dPSg1_drr,dPSg1_dss)

  RZ_jac  = dRRg1_dr * dZZg1_ds - dRRg1_ds * dZZg1_dr
  PSI_R  = (   dZZg1_ds * dPSg1_dr - dZZg1_dr * dPSg1_ds ) / RZ_jac
  PSI_Z  = ( - dRRg1_ds * dPSg1_dr + dRRg1_dr * dPSg1_ds ) / RZ_jac

  newnode_list%node(i)%values(1,1,1) = PSg1
  newnode_list%node(i)%values(1,2,1) = PSI_R * newnode_list%node(i)%x(1,2,1) + PSI_Z * newnode_list%node(i)%x(1,2,2)
  newnode_list%node(i)%values(1,3,1) = PSI_R * newnode_list%node(i)%x(1,3,1) + PSI_Z * newnode_list%node(i)%x(1,3,2)
  newnode_list%node(i)%values(1,4,1) = PSI_R * newnode_list%node(i)%x(1,4,1) + PSI_Z * newnode_list%node(i)%x(1,4,2)

  if (newnode_list%node(i)%boundary .eq. 2) newnode_list%node(i)%values(1,3,1) = 0.d0

enddo

newnode_list%node(index_xpoint  )%values(1,2:4,1) = 0.d0
newnode_list%node(index_xpoint+1)%values(1,2:4,1) = 0.d0
newnode_list%node(index_xpoint+2)%values(1,2:4,1) = 0.d0
newnode_list%node(index_xpoint+3)%values(1,2:4,1) = 0.d0

do i=1,n_tht - 1
  newnode_list%node(i)%values(1,2:4,1) = 0.d0
enddo

!----------------------------- empty old nodes/elements
do i=1,node_list%n_nodes
  node_list%node(i)%x        = 0.d0
  node_list%node(i)%values   = 0.d0
  node_list%node(i)%index    = 0
  node_list%node(i)%boundary = 0
enddo
node_list%n_nodes = 0

do i=1,element_list%n_elements
  element_list%element(i)%vertex     = 0
  element_list%element(i)%size       = 0.d0
  element_list%element(i)%neighbours = 0
enddo

!---------------------------- copy new grid into nodes/elements
node_list%n_nodes = newnode_list%n_nodes
node_list%node(1:node_list%n_nodes) = newnode_list%node(1:node_list%n_nodes)

element_list%n_elements = newelement_list%n_elements
element_list%element(1:element_list%n_elements) = newelement_list%element(1:element_list%n_elements)


!----temporary, needs to be completed, neighbour information
do i=1, element_list%n_elements
  element_list%element(i)%father = 0
  element_list%element(i)%n_sons = 0
  element_list%element(i)%sons(:) = 0
enddo

call tr_unregister_mem(sizeof(newnode_list),"newnode_list",CAT_GRID)
deallocate(newnode_list)
call tr_unregister_mem(sizeof(newelement_list),"newelement_list",CAT_GRID)
deallocate(newelement_list)

call find_axis(my_id,node_list,element_list,psi_axis,R_axis,Z_axis,i_elm_axis,s_axis,t_axis,ifail)

call find_xpoint(my_id,node_list,element_list,psi_xpoint,R_xpoint,Z_xpoint,i_elm_xpoint,s_xpoint,t_xpoint,xcase,ifail)

call tr_deallocate(s_values,"s_values",CAT_GRID)
call tr_deallocate(theta_sep,"theta_sep",CAT_GRID)
call tr_deallocate(R_sep,"R_sep",CAT_GRID)
call tr_deallocate(Z_sep,"Z_sep",CAT_GRID)
call tr_deallocate(R_max,"R_max",CAT_GRID)
call tr_deallocate(Z_max,"Z_max",CAT_GRID)
call tr_deallocate(R_min,"R_min",CAT_GRID)
call tr_deallocate(Z_min,"Z_min",CAT_GRID)
call tr_deallocate(s_tmp,"s_tmp",CAT_GRID)
call tr_deallocate(R_polar,"R_polar",CAT_GRID)
call tr_deallocate(Z_polar,"Z_polar",CAT_GRID)
call tr_deallocate(RR_new,"RR_new",CAT_GRID)
call tr_deallocate(ZZ_new,"ZZ_new",CAT_GRID)
call tr_deallocate(s_flux,"s_flux",CAT_GRID)
call tr_deallocate(t_flux,"t_flux",CAT_GRID)
call tr_deallocate(t_tht,"t_tht",CAT_GRID)
call tr_deallocate(ielm_flux,"ielm_flux",CAT_GRID)
call tr_deallocate(keep,"keep",CAT_GRID)
call tr_deallocate(k_cross,"k_cross",CAT_GRID)

call update_neighbours(node_list,element_list, force_rtree_initialize=.true.)
return
end subroutine grid_xpoint
