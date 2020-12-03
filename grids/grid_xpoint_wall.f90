subroutine grid_xpoint_wall(node_list, element_list, n_flux, n_open, n_private, n_leg, n_tht, n_ext,&
  SIG_open, SIG_closed, SIG_private, SIG_theta, SIG_leg_0, SIG_leg_1, dPSI_open, dPSI_private)
!-----------------------------------------------------------------------
! subroutine defines a flux surface aligned finite element grid
! including a single x-point.
! Add shape of exetrnal wall, align to wall close to it.
!-----------------------------------------------------------------------

use data_structure
use tr_module 
use gauss
use basis_at_gaussian
use phys_module, only:   n_limiter, R_limiter, Z_limiter, write_ps, fix_axis_nodes, force_central_node
use mod_neighbours, only: update_neighbours
use mod_interp

implicit none

! --- Routine parameters
type (type_node_list),    intent(inout) :: node_list
type (type_element_list), intent(inout) :: element_list
integer,                  intent(in)    :: n_flux, n_open, n_private, n_leg, n_tht
real*8,                   intent(in)    :: SIG_open, SIG_closed, SIG_private, SIG_theta, SIG_leg_0, SIG_leg_1
real*8,                   intent(in)    :: dPSI_open, dPSI_private

! --- local variables
type (type_surface_list) :: flux_list

type (type_node_list), pointer    :: newnode_list
type (type_element_list), pointer :: newelement_list
type (type_element)                   :: element
type (type_node)                      :: nodes(n_vertex_max)

real*8, allocatable :: s_values(:), theta_sep(:), R_sep(:), Z_sep(:), R_max(:), Z_max(:), R_min(:), Z_min(:),s_tmp(:)
real*8, allocatable :: R_wall_max(:), Z_wall_max(:), T_wall_par(:), R_wall_min(:), Z_wall_min(:)
real*8              :: psi_axis, R_axis, Z_axis, s_axis, t_axis
real*8              :: R_xpoint(2), Z_xpoint(2), s_xpoint(2), t_xpoint(2), psi_xpoint(2)
real*8              :: PI, s_find(8), t_find(8), st_find(8), tht_x, theta, delta, ss, tmp1, tmp2, tan_max, tht_bnd, tht_ext
real*8              :: RRg1,dRRg1_dr,dRRg1_ds
real*8              :: ZZg1,dZZg1_dr,dZZg1_ds
real*8              :: PSg1,dPSg1_dr,dPSg1_ds,dPSg1_drs,dPSg1_drr,dPSg1_dss
real*8,allocatable  :: R_polar(:,:,:),Z_polar(:,:,:),xout(:),xp(:),yp(:), R_strike(:,:), Z_strike(:,:)
real*8              :: R_cub1d(4), Z_cub1d(4), dR_dt, dZ_dt, RZ_jac, PSI_R, PSI_Z, Rw, Rw2, Zw, Zw2, Tw, tan12
real*8, allocatable :: RR_new(:,:),ZZ_new(:,:),s_flux(:,:),t_flux(:,:),t_tht(:,:), R_wall(:), Z_wall(:)
integer,allocatable :: ielm_flux(:,:), k_cross(:,:), elm_left(:), elm_right(:)
integer             :: my_id, i, j, k, l, m, n_psi, n_theta, i2, j2, n_total, n_tht_3
integer             :: i_surf, n_pieces, n_wall, i_flux, n_ext, j2rev
integer             :: i_elm_axis, i_elm_xpoint(2), i_elm_find(8), i_sep, i_max, i_find, npl, ifail
integer             :: node, index, node_start, index_xpoint, n_xpoint, j_start, j_end
integer             :: iv, ivp, node_iv, node_ivp, i_elm, ielm_out
integer             :: ms, mt, inode, iter, index_keep, xcase
integer             :: index_ext1, index_ext2, index_ext_leg1, index_ext_leg2, index_leg1, index_leg2
real*8              :: Rmid, Zmid, R0,Z0, RP,ZP, dR0, dZ0, dRP, dZP, size_0, size_p, denom
real*8              :: R1, Z1, s_out, t_out, R_out, Z_out, R_private, Z_private, tht1, tht2, tht_min, tht_max
real*8              :: EJAC, RX, RY, SX, SY, CRR, CZZ, CRZ, alpha1, alpha2, alpha_max, alpha_min
real*8              :: RL5, RL8, RL9, RL10, RL11, ZL5, ZL8, ZL9, ZL10, ZL11, angle_L8, angle_L9, rr1, ss1
real*8,allocatable  :: psi_gaussians(:,:), angle_gaussians(:,:), s_equidistant(:)
real*8,allocatable  :: Aspline(:), Bspline(:), Cspline(:), Dspline(:)
real*8              :: x_g(n_gauss,n_gauss), y_g(n_gauss,n_gauss), psi_g(n_gauss,n_gauss), xmin(2), xmax(2)
real*8              :: abltg(3), t_node, sign_psi
real*8              :: Rtmp, Ztmp, dRtmp, dZtmp, Rtmp2, Ztmp2, dRtmp2, dZtmp2, Rtmp3, Ztmp3, dRtmp3, dZtmp3
real*8              :: t2, t3, t_delta, t_total, xl_axis, theta_axis
logical             :: xpoint, extend
real*8,external     :: root, spwert
character*4         :: label
integer             :: i_elm1, i_vertex1, i_node1, i_node_save, iv1, iv2, iv3, iv4
integer             :: i_elm2, i_vertex2, i_node2
integer             :: n_remove_elements, n_remove_nodes, remove_elements(100), remove_nodes(100), newnode_index(n_nodes_max), skip_index

xpoint = .true.
extend = .true.;   if (n_ext .lt. 1) extend = .false.
xcase  = 1
my_id  = 0

PI = 2.d0 * asin(1.d0)

write(*,*) '*************************************'
write(*,*) '*    X-point(wall) grid             *'
write(*,*) '*************************************'

open(23,file='wall.txt')
read(23,*)
read(23,*) n_wall
write(*,*) ' reading outer wall : ',n_wall
allocate(R_wall(n_wall),Z_wall(n_wall))
read(23,*)
do i=1,n_wall
  read(23,*) R_wall(i),Z_wall(i)
enddo
close(23)

call find_axis(my_id,node_list,element_list,psi_axis,R_axis,Z_axis,i_elm_axis,s_axis,t_axis,ifail)

call find_xpoint(my_id,node_list,element_list,psi_xpoint,R_xpoint,Z_xpoint,i_elm_xpoint,s_xpoint,t_xpoint,xcase,ifail)




! -------------------------------------------------------------------------------------------
! --------------------------------- Define Flux Values --------------------------------------
! -------------------------------------------------------------------------------------------


flux_list%n_psi = n_flux - 1 + n_open + n_private            ! excludes the axis

allocate(flux_list%psi_values(flux_list%n_psi))
allocate(s_values(n_flux+n_open+n_private),psi_gaussians(4,n_flux+n_open+n_private))

allocate(s_tmp(n_flux),Aspline(n_flux),Bspline(n_flux),Cspline(n_flux),Dspline(n_flux),s_equidistant(n_flux))
s_tmp = 0
call meshac2(n_flux,s_tmp,1.d0,9999.d0,SIG_closed,9999.d0,0.2d0,1.0d0)

do i=1,n_flux
  s_equidistant(i) = float(i-1)/float(n_flux-1)
enddo
call spline(n_flux,s_equidistant,s_tmp,0.d0,0.d0,2,Aspline,Bspline,Cspline,Dspline)

do i=1,n_flux-1
  s_values(i) = s_tmp(i+1)
  flux_list%psi_values(i) =  psi_axis + s_values(i)**2 * (psi_xpoint(1) - psi_axis)
enddo

deallocate(s_tmp,Aspline,Bspline,Cspline,Dspline,s_equidistant)
allocate(s_tmp(n_open+1),Aspline(n_open+1),Bspline(n_open+1),Cspline(n_open+1),Dspline(n_open+1),s_equidistant(n_open+1))
s_tmp = 0
call meshac2(n_open+1,s_tmp,0.d0,9999.d0,SIG_open,9999.d0,0.6d0,1.0d0)

do i=1,n_open+1
  s_equidistant(i) = float(i-1)/float(n_open)
enddo
call spline(n_open+1,s_equidistant,s_tmp,0.d0,0.d0,2,Aspline,Bspline,Cspline,Dspline)

do i=1,n_open
  s_values(i+n_flux-1)             =  1.d0 + dPSI_open*s_tmp(i+1)
  flux_list%psi_values(i+n_flux-1) =  psi_axis + s_values(i+n_flux-1)**2 * (psi_xpoint(1) - psi_axis)
enddo

deallocate(s_tmp,Aspline,Bspline,Cspline,Dspline,s_equidistant)
allocate(s_tmp(n_private+1),Aspline(n_private+1),Bspline(n_private+1), &
                            Cspline(n_private+1),Dspline(n_private+1),s_equidistant(n_private+1))

s_tmp = 0
call meshac2(n_private+1,s_tmp,0.d0,9999.d0,SIG_private,9999.d0,0.6d0,1.0d0)

do i=1,n_private+1
  s_equidistant(i) = float(i-1)/float(n_private)
enddo
call spline(n_private+1,s_equidistant,s_tmp,0.d0,0.d0,2,Aspline,Bspline,Cspline,Dspline)

do i=1,n_private
  s_values(i+n_flux-1+n_open)             =  1.d0 - dPSI_private*s_tmp(i+1)
  flux_list%psi_values(i+n_flux-1+n_open) =  psi_axis + s_values(i+n_flux-1+n_open)**2 * (psi_xpoint(1) - psi_axis)
enddo

call find_flux_surfaces(my_id,xpoint,xcase,node_list,element_list,flux_list)

call plot_flux_surfaces(node_list,element_list,flux_list,.true.,1,psi_xpoint,R_xpoint,Z_xpoint,xpoint,xcase)

if ( write_ps ) then 
  call lincol(3)
  call lplot6(21,11,R_wall,Z_wall,-n_wall,'first wall')
  call lincol(0)
endif




! -------------------------------------------------------------------------------------------
! --------------------------------- Find Strategic Points -----------------------------------
! -------------------------------------------------------------------------------------------

!-------------------------------------- store some data for the new grid
n_psi   = n_flux + n_open + n_private          ! includes the axis
n_theta = n_tht + 2*n_leg                      ! includes center and legs

n_total = n_psi + n_ext                        ! add number of non-aligned surfaces

allocate(RR_new(n_total,n_theta),ZZ_new(n_total,n_theta))
allocate(ielm_flux(n_total,n_theta),s_flux(n_total,n_theta))
allocate(t_flux(n_total,n_theta),t_tht(n_total,n_theta))

tht_x = atan2(Z_xpoint(1)-Z_axis,R_xpoint(1)-R_axis)
if (tht_x .lt. 0.d0) tht_x = tht_x + 2.d0 * PI

!-------------------------------- define the polar coordinate
i_sep = n_flux - 1               ! separatrix surface (index in fluxsurface list, excludes axis)
i_max = n_flux + n_open - 1      ! last open surface (non-private) idem

allocate(theta_sep(n_theta),R_sep(n_theta),Z_sep(n_theta))

deallocate(s_tmp); allocate(s_tmp(n_tht))
s_tmp = 0
call meshac2(n_tht,s_tmp,0.d0,1.d0,SIG_theta,SIG_theta,0.8d0,1.0d0)

do j=1,n_tht
  theta_sep(j) = tht_x + 2.d0 * PI * s_tmp(j)
  if (theta_sep(j) .lt. 0.d0)      theta_sep(j) = theta_sep(j) + 2.d0 * PI
  if (theta_sep(j) .gt. 2.d0 * PI) theta_sep(j) = theta_sep(j) - 2.d0 * PI
enddo

!------------------------------------- assume straight line through x-point (but not necesarily horizontal)
!                                      find crossing with first wall
write(*,*) ' tht_x : ',tht_x
tht_max = PI/2.d0+tht_x
write(*,*) ' L9 : ',tht_max
call find_wall_crossing(R_wall,Z_wall,n_wall,R_xpoint(1),Z_xpoint(1),tht_max,RL9,ZL9,Tw)

tht_max = 3.d0*PI/2.d0+tht_x
write(*,*) ' L8 : ',tht_max
call find_wall_crossing(R_wall,Z_wall,n_wall,R_xpoint(1),Z_xpoint(1),tht_max,RL8,ZL8,Tw)

write(*,'(A,2e16.8)') ' L8 (inside wall crossing) : ',RL8,ZL8  ! inside wall crossing at x-point level
write(*,'(A,2e16.8)') ' L9 (outside wall crossing): ',RL9,ZL9  ! outside wall crossing at x-point level




! -------------------------------------------------------------------------------------------
! --------------------------------- Define New Grid Points ----------------------------------
! -------------------------------------------------------------------------------------------

!------------------------------------- find crossing with separatrix
do j=2,n_tht-1

  call find_theta_surface(node_list,element_list,flux_list,i_sep,theta_sep(j),R_axis,Z_axis,i_elm_find,s_find,t_find,i_find)

  call interp_RZ(node_list,element_list,i_elm_find(1),s_find(1),t_find(1),RRg1,ZZg1)

  R_sep(j) = RRg1
  Z_sep(j) = ZZg1

enddo

R_sep(1) = R_xpoint(1)
Z_sep(1) = Z_xpoint(1)
R_sep(n_tht) = R_xpoint(1)
Z_sep(n_tht) = Z_xpoint(1)

if ( write_ps ) then
  call lincol(2)
  call lplot6(21,11,R_sep,Z_sep,-n_tht,' ')
  call lincol(0)
endif

!------------------------------------ find crossing with the first wall
allocate(R_wall_max(n_tht+2*n_leg),Z_wall_max(n_tht+2*n_leg),T_wall_par(n_tht+2*n_leg))
allocate(R_wall_min(n_tht+2*n_leg),Z_wall_min(n_tht+2*n_leg))

do j=1,n_tht

  if (Z_sep(j) .le. Z_axis) then

    if (j .gt. n_tht/2) then
      Z_wall_max(j) = Z_sep(j) + (ZL8 - Z_xpoint(1)) * ((Z_sep(j) - Z_axis)/(Z_xpoint(1) - Z_axis))**2
      tht_max = PI
    else
      Z_wall_max(j) = Z_sep(j) + (ZL9 - Z_xpoint(1)) * ((Z_sep(j) - Z_axis)/(Z_xpoint(1) - Z_axis))**2
      tht_max = 0.d0
    endif

    call find_wall_crossing(R_wall,Z_wall,n_wall,R_sep(j),Z_wall_max(j),tht_max,Rw,Zw,Tw)
  else
    tht_max = theta_sep(j)
    call find_wall_crossing(R_wall,Z_wall,n_wall,R_sep(j),Z_sep(j),tht_max,Rw,Zw,Tw)
  endif

  R_wall_max(j) = Rw
  Z_wall_max(j) = Zw
  T_wall_par(j) = Tw
enddo

angle_L8 = tht_x + 1.5d0*PI
angle_L9 = tht_x + 0.5d0*PI
if (angle_L8 .gt. PI) angle_L8 = angle_L8 - 2.d0*PI
if (angle_L9 .gt. PI) angle_L9 = angle_L9 - 2.d0*PI

R_wall_max(1)     = RL9
Z_wall_max(1)     = ZL9
R_wall_max(n_tht) = RL8
Z_wall_max(n_tht) = ZL8
!T_wall_par(1)     = angle_L9 + PI/2.d0
!T_wall_par(n_tht) = angle_L8 + PI/2.d0

if ( write_ps ) then
  call lincol(6)  
  call lplot6(21,11,R_wall_max,Z_wall_max,-n_tht,' ')
  call lincol(0)
endif

!------------------------------------ find crossing of last fluxsurface
n_tht_3 = n_tht + 2*n_leg

call tr_allocate(R_max,1,n_tht_3,"R_max")
call tr_allocate(Z_max,1,n_tht_3,"Z_max")
call tr_allocate(R_min,1,n_tht_3,"R_min")
call tr_allocate(Z_min,1,n_tht_3,"Z_min")

call find_theta_surface(node_list,element_list,flux_list,i_max,angle_L9,R_xpoint(1),Z_xpoint(1), &
                        i_elm_find,s_find,t_find,i_find)

call interp_RZ(node_list,element_list,i_elm_find(1),s_find(1),t_find(1),RL10,ZL10)

write(*,'(A,2e16.8)') ' L10 (last flux surface level x-point, outside): ',RL10,ZL10 ! outside wall crossing at x-point level

call find_theta_surface(node_list,element_list,flux_list,i_max,angle_L8,R_xpoint(1),Z_xpoint(1), &
                        i_elm_find,s_find,t_find,i_find)

call interp_RZ(node_list,element_list,i_elm_find(1),s_find(1),t_find(1),RL11,ZL11)

write(*,'(A,2e16.8)') ' L11 (last flux surface level x-point, inside): ',RL11,ZL11  ! inside wall crossing at x-point level


!------------------------------ second part of the grid below the x-point

call find_theta_surface(node_list,element_list,flux_list,flux_list%n_psi,tht_x,R_axis,Z_axis,i_elm_find,s_find,t_find,i_find)

do i=1,i_find
   call interp_RZ(node_list,element_list,i_elm_find(i),s_find(i),t_find(i),RL5,ZL5)

   if (ZL5 .le. Z_xpoint(1)) exit
enddo

write(*,'(A,2e16.8)') ' L5 (last private flux surface, below x-point): ',RL5,ZL5

!------------------------------ find wall crossing of open and private field lines (i.e strike points)

allocate(R_strike(n_open+n_private+1,2),Z_strike(n_open+n_private+1,2))

R_strike(:,1) = 1.d10
R_strike(:,2) = 0.d0
Z_strike      = 0.d0

do i=1, n_open + n_private + 1

  i_flux = n_flux -1 + i - 1

  do k=1,n_wall-1

    Rw = R_wall(k+1)
    Zw = Z_wall(k+1)
    theta = atan2(Z_wall(k)-Zw, R_wall(k)-Rw)

    if (R_wall(k) .eq. R_wall(k+1)) then

      call find_R_surface(node_list,element_list,flux_list,i_flux,Rw,i_elm_find,s_find,t_find,st_find,i_find)

    elseif (Z_wall(k) .eq. Z_wall(k+1)) then

      call find_Z_surface(node_list,element_list,flux_list,i_flux,Rw,i_elm_find,s_find,t_find,st_find,i_find)

    else

       call find_theta_surface(node_list,element_list,flux_list,i_flux,theta,Rw,Zw,i_elm_find,s_find,t_find,i_find)

    endif

    do j=1,i_find

      call interp_RZ(node_list,element_list,i_elm_find(j),s_find(j),t_find(j),RRg1,ZZg1)

      if   ( (((R_wall(k)-RRg1)*(R_wall(k+1)-RRg1) .le. 0.d0) .or. (R_wall(k) .eq. R_wall(k+1))) &
      
           .and. ((Z_wall(k)-ZZg1)*(Z_wall(k+1)-ZZg1) .le. 0.d0) )  then

        if ((RRg1 .le. R_xpoint(1)) .and. (ZZg1 .le. Z_axis)) then
            R_strike(i,1) = RRg1
            Z_strike(i,1) = ZZg1
!            write(*,'(A,i3,2f8.4)') ' INNER strike point : ',i,RRg1,ZZg1
        elseif ((RRg1 .gt. R_xpoint(1)) .and. (ZZg1 .le. Z_axis)) then
             R_strike(i,2) = RRg1
             Z_strike(i,2) = ZZg1
!             write(*,'(A,i3,2f8.4)') ' OUTER strike point : ',i,RRg1,ZZg1
         endif

      endif

    enddo
  enddo
enddo

if ( write_ps ) then
  call lincol(2)
  call lplot(1,1,461,R_strike(:,1),Z_strike(:,1),-(n_open+n_private+1),1,'Nodes',5,'X',1,'Y',1)
  call lincol(3)
  call lplot(1,1,461,R_strike(:,2),Z_strike(:,2),-(n_open+n_private+1),1,'Nodes',5,'X',1,'Y',1)
  call lincol(0)
endif

deallocate(s_tmp); allocate(s_tmp(n_leg))
s_tmp = 0
call meshac2(n_leg,s_tmp,0.d0,1.d0,SIG_leg_0,SIG_leg_1,0.6d0,1.0d0)

!----------------------------- inner leg, private side
do j=1,n_leg

  R_min(n_tht + j) = R_strike(n_open+n_private+1,1) + (RL5-R_strike(n_open+n_private+1,1)) * s_tmp(j)

  call find_R_surface(node_list,element_list,flux_list,n_psi-1,R_min(n_tht+j),i_elm_find,s_find,t_find,st_find,i_find)

  do i=1,i_find

    call interp_RZ(node_list,element_list,i_elm_find(i),s_find(i),t_find(i),RRg1,ZZg1)

    if (ZZg1 .le. Z_xpoint(1)) exit

  enddo

  Z_min(n_tht + j) = ZZg1

enddo

!----------------------------- outer leg, private side
do j=1,n_leg

  R_min(n_tht + n_leg + j) = R_strike(n_open+n_private+1,2) + (RL5-R_strike(n_open+n_private+1,2)) * s_tmp(j)

  call find_R_surface(node_list,element_list,flux_list,n_psi-1,R_min(n_tht+n_leg+j),i_elm_find,s_find,t_find,st_find,i_find)

  do i=1,i_find

    call interp_RZ(node_list,element_list,i_elm_find(i),s_find(i),t_find(i),RRg1,ZZg1)

    if (ZZg1 .le. Z_xpoint(1)) exit

  enddo

  Z_min(n_tht + n_leg + j) = ZZg1

enddo



!----------------------------- inner leg
do j=1,n_leg

  Z_max(n_tht + j) = Z_strike(n_open+1,1) + (ZL11 - Z_strike(n_open+1,1)) * s_tmp(j)

  call find_Z_surface(node_list,element_list,flux_list,i_max,Z_max(n_tht+j),i_elm_find,s_find,t_find,st_find,i_find)

  do i=1,i_find

    call interp_RZ(node_list,element_list,i_elm_find(i),s_find(i),t_find(i),RRg1,ZZg1)

    if (RRg1 .le. R_xpoint(1)) exit

  enddo

  R_max(n_tht + j) = RRg1

enddo

!----------------------------- outer leg, SOL
do j=1,n_leg

  Z_max(n_tht + n_leg + j) = Z_strike(n_open+1,2) + (ZL10 - Z_strike(n_open+1,2)) * s_tmp(j)

  call find_Z_surface(node_list,element_list,flux_list,i_max,Z_max(n_tht+n_leg+j),i_elm_find,s_find,t_find,st_find,i_find)

  do i=1,i_find

    call interp_RZ(node_list,element_list,i_elm_find(i),s_find(i),t_find(i),RRg1,ZZg1)

    if (RRg1 .ge. R_xpoint(1)) exit

  enddo

  R_max(n_tht + n_leg + j) = RRg1

enddo

!--------------------------- along the speratrix of on the inner leg
do j=1,n_leg

  R_sep(n_tht + j) =  R_strike(1,1) + (R_xpoint(1) - R_strike(1,1)) * s_tmp(j)

  call find_R_surface(node_list,element_list,flux_list,i_sep,R_sep(n_tht+j),i_elm_find,s_find,t_find,st_find,i_find)

  do i=1,i_find

    call interp_RZ(node_list,element_list,i_elm_find(i),s_find(i),t_find(i),RRg1,ZZg1)

    if (ZZg1 .le. Z_xpoint(1)) exit

  enddo

  Z_sep(n_tht  + j) = ZZg1

enddo
Z_sep(n_tht+n_leg) = Z_xpoint(1) ! this one is known

!--------------------------- along the speratrix of on the outer leg
do j=1,n_leg

  R_sep(n_tht + n_leg + j) = R_strike(1,2) + (R_xpoint(1) - R_strike(1,2)) * s_tmp(j)
  call find_R_surface(node_list,element_list,flux_list,i_sep,R_sep(n_tht+n_leg+j),i_elm_find,s_find,t_find,st_find,i_find)

  do i=1,i_find

    call interp_RZ(node_list,element_list,i_elm_find(i),s_find(i),t_find(i),RRg1,ZZg1)

    if (ZZg1 .le. Z_xpoint(1)) exit

  enddo

  Z_sep(n_tht + n_leg + j) = ZZg1

enddo
Z_sep(n_tht+2*n_leg) = Z_xpoint(1) ! this one is known

s_tmp = 0.d0
call meshac2(n_leg,s_tmp,0.d0,1.d0,0.3d0,0.3d0,0.6d0,1.0d0)

do j=1,n_leg
  R_wall_max(j+n_tht) = R_wall_max(n_tht) + 0.8d0*(R_max(n_tht+1) - R_wall_max(n_tht)) * float(j-1)/float(n_leg-1)
  Z_wall_max(j+n_tht) = Z_wall_max(n_tht) + 0.8d0*(Z_max(n_tht+1) - Z_wall_max(n_tht)) * float(j-1)/float(n_leg-1) 
  tht_max = PI
  call find_wall_crossing(R_wall,Z_wall,n_wall,R_max(n_tht+n_leg-j+1)+0.1,Z_wall_max(j+n_tht),tht_max,Rw,Zw,Tw)
  R_wall_max(j+n_tht) = Rw
  Z_wall_max(j+n_tht) = Zw
  T_wall_par(j+n_tht) = Tw
enddo
do j=1,n_leg
  R_wall_max(j+n_tht+n_leg) = R_wall_max(1) + 0.8d0*(R_max(n_tht+n_leg+1)-R_wall_max(1)) * float(j-1)/float(n_leg-1)
  Z_wall_max(j+n_tht+n_leg) = Z_wall_max(1) + 0.8d0*(Z_max(n_tht+n_leg+1)-Z_wall_max(1)) * float(j-1)/float(n_leg-1)
  tht_max = 0.d0
  call find_wall_crossing(R_wall,Z_wall,n_wall,R_max(j+n_tht+n_leg)-0.1,Z_wall_max(j+n_tht+n_leg),tht_max,Rw,Zw,Tw)
  R_wall_max(j+n_tht+n_leg) = Rw
  Z_wall_max(j+n_tht+n_leg) = Zw
  T_wall_par(j+n_tht+n_leg) = Tw
enddo


call plot_flux_surfaces(node_list,element_list,flux_list,.true.,1,psi_xpoint,R_xpoint,Z_xpoint,xpoint,xcase)

if ( write_ps ) then
  call lincol(3)
  call lplot6(21,11,R_wall,Z_wall,-n_wall,'first wall')
  call lincol(0)

  call lincol(2)
  call lplot6(1,1,R_max,Z_max,-(n_tht+2*n_leg),' ')
  call lincol(5)
  call lplot6(1,1,R_min(n_tht+1),Z_min(n_tht+1),-(2*n_leg),' ')
  call lincol(4)
  call lplot6(1,1,R_sep,Z_sep,-(n_tht+2*n_leg),' ')
  call lincol(0)
endif

!do i=1,n_tht+2*n_leg
!  write(*,'(A,i4,6f10.6)') ' R_max,Z_max : ',i,R_max(i),Z_max(i),R_sep(i),Z_sep(i),R_min(i),Z_min(i)
!enddo

!------------------------------ interpolation points are known, construct polar coordinate lines
n_pieces=3
allocate(R_polar(n_pieces,4,n_tht+2*n_leg),Z_polar(n_pieces,4,n_tht+2*n_leg))

do j=1,n_tht

  delta = 0.08

  if ((j .eq. 1) .or. (j .eq. n_tht))       delta = 0.d0
  if ((j .eq. 2) .or. (j .eq. n_tht - 1))   delta = 0.05d0

  theta_axis = tht_x + 2.*PI*(j-1)/(n_tht-1)

  R_polar(1,1,j) = R_axis
  R_polar(1,4,j) = delta * R_axis + (1.d0 - delta) * R_sep(j)
  R_polar(1,2,j) = ( 2.d0 * R_polar(1,1,j)  +         R_polar(1,4,j) ) / 3.d0
  R_polar(1,3,j) = (        R_polar(1,1,j)  +  2.d0 * R_polar(1,4,j) ) / 3.d0

  Z_polar(1,1,j) = Z_axis
  Z_polar(1,4,j) = delta * Z_axis + (1.d0 - delta) * Z_sep(j)
  Z_polar(1,2,j) = ( 2.d0 * Z_polar(1,1,j)  +         Z_polar(1,4,j) ) / 3.d0
  Z_polar(1,3,j) = (        Z_polar(1,1,j)  +  2.d0 * Z_polar(1,4,j) ) / 3.d0

  xl_axis = sqrt((R_polar(1,2,j) -  R_axis)**2 + (Z_polar(1,2,j) -  Z_axis)**2 )

  R_polar(1,2,j) = R_axis +  0.8 * xl_axis * cos(theta_axis)
  Z_polar(1,2,j) = Z_axis +  0.8 * xl_axis * sin(theta_axis)

  R_polar(3,1,j) = R_wall_max(j)
  R_polar(3,4,j) = delta * R_wall_max(j) + (1.d0 - delta) * R_sep(j)
  R_polar(3,2,j) = ( 2.d0 * R_polar(3,1,j)  +         R_polar(3,4,j) ) / 3.d0
  R_polar(3,3,j) = (        R_polar(3,1,j)  +  2.d0 * R_polar(3,4,j) ) / 3.d0

  Z_polar(3,1,j) = Z_wall_max(j)
  Z_polar(3,4,j) = delta * Z_wall_max(j) + (1.d0 - delta) * Z_sep(j)
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



do j=1,2*n_leg

  delta = 0.2

  if ((j .eq. n_leg)   .or. (j .eq. 2*n_leg))    delta = 0.d0
  if ((j .eq. n_leg-1) .or. (j .eq. 2*n_leg-1))  delta = 0.05d0
  if ((j .eq. n_leg-2) .or. (j .eq. 2*n_leg-2))  delta = 0.1d0

  if (j .le. n_leg) then
    j2 = n_leg -j + 1
  else
    j2 = 2*n_leg - (j-n_leg) + 1
  endif

  R_polar(1,1,n_tht+j) = R_min(n_tht+j)
  R_polar(1,4,n_tht+j) = delta * R_min(n_tht+j) + (1.d0 - delta) * R_sep(n_tht+j)
  R_polar(1,2,n_tht+j) = ( 2.d0 * R_polar(1,1,n_tht+j)  +         R_polar(1,4,n_tht+j) ) / 3.d0
  R_polar(1,3,n_tht+j) = (        R_polar(1,1,n_tht+j)  +  2.d0 * R_polar(1,4,n_tht+j) ) / 3.d0

  Z_polar(1,1,n_tht+j) = Z_min(n_tht+j)
  Z_polar(1,4,n_tht+j) = delta * Z_min(n_tht+j) + (1.d0 - delta) * Z_sep(n_tht+j)
  Z_polar(1,2,n_tht+j) = ( 2.d0 * Z_polar(1,1,n_tht+j)  +         Z_polar(1,4,n_tht+j) ) / 3.d0
  Z_polar(1,3,n_tht+j) = (        Z_polar(1,1,n_tht+j)  +  2.d0 * Z_polar(1,4,n_tht+j) ) / 3.d0

  R_polar(3,1,n_tht+j) = R_wall_max(n_tht+j2)
  R_polar(3,4,n_tht+j) = delta * R_wall_max(n_tht+j2) + (1.d0 - delta) * R_sep(n_tht+j)
  R_polar(3,2,n_tht+j) = ( 2.d0 * R_polar(3,1,n_tht+j)  +         R_polar(3,4,n_tht+j) ) / 3.d0
  R_polar(3,3,n_tht+j) = (        R_polar(3,1,n_tht+j)  +  2.d0 * R_polar(3,4,n_tht+j) ) / 3.d0

  Z_polar(3,1,n_tht+j) = Z_wall_max(n_tht+j2)
  Z_polar(3,4,n_tht+j) = delta * Z_wall_max(n_tht+j2) + (1.d0 - delta) * Z_sep(n_tht+j)
  Z_polar(3,2,n_tht+j) = ( 2.d0 * Z_polar(3,1,n_tht+j)  +         Z_polar(3,4,n_tht+j) ) / 3.d0
  Z_polar(3,3,n_tht+j) = (        Z_polar(3,1,n_tht+j)  +  2.d0 * Z_polar(3,4,n_tht+j) ) / 3.d0

  R_polar(2,1,n_tht+j) = R_polar(1,4,n_tht+j)
  R_polar(2,4,n_tht+j) = R_polar(3,4,n_tht+j)
  R_polar(2,2,n_tht+j) = ( R_polar(2,1,n_tht+j) +  2.d0 * R_sep(n_tht+j) ) / 3.d0
  R_polar(2,3,n_tht+j) = ( R_polar(2,4,n_tht+j) +  2.d0 * R_sep(n_tht+j) ) / 3.d0

  Z_polar(2,1,n_tht+j) = Z_polar(1,4,n_tht+j)
  Z_polar(2,4,n_tht+j) = Z_polar(3,4,n_tht+j)
  Z_polar(2,2,n_tht+j) = ( Z_polar(2,1,n_tht+j) +  2.d0 * Z_sep(n_tht+j) ) / 3.d0
  Z_polar(2,3,n_tht+j) = ( Z_polar(2,4,n_tht+j) +  2.d0 * Z_sep(n_tht+j) ) / 3.d0

#ifdef D3D_WALL
  if (j .eq. 1) then
    R_polar(3,1,n_tht+j) = R_wall_max(n_tht+j2)
    R_polar(3,4,n_tht+j) = R_wall_max(n_tht+j2) 
    R_polar(3,2,n_tht+j) = R_wall_max(n_tht+j2) 
    R_polar(3,3,n_tht+j) = R_wall_max(n_tht+j2) 

    Z_polar(3,4,n_tht+j) = -1.223 

    Z_polar(3,2,n_tht+j) = (2.d0 * Z_polar(3,1,n_tht+j) +        Z_polar(3,4,n_tht+j)) / 3.d0
    Z_polar(3,3,n_tht+j) = (       Z_polar(3,1,n_tht+j) + 2.d0 * Z_polar(3,4,n_tht+j)) / 3.d0

    R_polar(2,4,n_tht+j) = R_wall_max(n_tht+j2) 
    Z_polar(2,4,n_tht+j) = -1.223 

  endif

  if ((j .ge. 2) .and. (j .le. 5)) then

    R_polar(3,4,n_tht+j) =  0.25d0 * R_polar(3,4,n_tht+j) + 0.75d0 * R_polar(3,4,n_tht+j-1)
    
    R_polar(3,2,n_tht+j) = (2.0d0 * R_polar(3,1,n_tht+j) +         R_polar(3,4,n_tht+j)) / 3.d0
    R_polar(3,3,n_tht+j) = (        R_polar(3,1,n_tht+j) + 2.0d0 * R_polar(3,4,n_tht+j)) / 3.d0

    Z_polar(3,4,n_tht+j) = -1.223 

    Z_polar(3,2,n_tht+j) = (2.d0 * Z_polar(3,1,n_tht+j) +        Z_polar(3,4,n_tht+j)) / 3.d0
    Z_polar(3,3,n_tht+j) = (       Z_polar(3,1,n_tht+j) + 2.d0 * Z_polar(3,4,n_tht+j)) / 3.d0

    R_polar(2,4,n_tht+j) = R_polar(3,4,n_tht+j) 
    Z_polar(2,4,n_tht+j) = Z_polar(3,4,n_tht+j) 

  endif
#endif

enddo

#ifdef AUG_WALL

R_polar(3,2,n_tht+1) =  1.294 
Z_polar(3,2,n_tht+1) = -0.955

R_polar(3,3,n_tht+1) =  1.292 
Z_polar(3,3,n_tht+1) = -0.983

R_polar(3,4,n_tht+1) =  1.287 
Z_polar(3,4,n_tht+1) = -1.001

R_polar(2,4,n_tht+1) = R_polar(3,4,n_tht+1) 
Z_polar(2,4,n_tht+1) = Z_polar(3,4,n_tht+1) 

Z_polar(2,3,n_tht+1) = 2.d0 * Z_polar(2,4,n_tht+1) - Z_polar(3,3,n_tht+1)  
R_polar(2,3,n_tht+1) = 2.d0 * R_polar(2,4,n_tht+1) - R_polar(3,3,n_tht+1)  

j_end = n_leg/2

do j= 2, j_end-1

  R_polar(3,:,n_tht+j) =  real(j-1,8)/real(j_end-1,8) * R_polar(3,:,n_tht+j_end) + real(j_end-j,8)/real(j_end-1,8) * R_polar(3,:,n_tht+1)
  Z_polar(3,:,n_tht+j) =  real(j-1,8)/real(j_end-1,8) * Z_polar(3,:,n_tht+j_end) + real(j_end-j,8)/real(j_end-1,8) * Z_polar(3,:,n_tht+1)
  
  R_polar(2,4,n_tht+j) = R_polar(3,4,n_tht+j) 
  Z_polar(2,4,n_tht+j) = Z_polar(3,4,n_tht+j) 

  Z_polar(2,3,n_tht+j) = 2.d0 * Z_polar(2,4,n_tht+j) - Z_polar(3,3,n_tht+j)  
  R_polar(2,3,n_tht+j) = 2.d0 * R_polar(2,4,n_tht+j) - R_polar(3,3,n_tht+j)  

enddo

R_polar(3,:,n_tht+2) = 0.5d0 * R_polar(3,:,n_tht+1) + 0.5d0 * R_polar(3,:,n_tht+3)
Z_polar(3,:,n_tht+2) = 0.5d0 * Z_polar(3,:,n_tht+1) + 0.5d0 * Z_polar(3,:,n_tht+3)
  
R_polar(2,4,n_tht+2) = R_polar(3,4,n_tht+2) 
Z_polar(2,4,n_tht+2) = Z_polar(3,4,n_tht+2) 

Z_polar(2,3,n_tht+2) = 2.d0 * Z_polar(2,4,n_tht+2) - Z_polar(3,3,n_tht+2)  
R_polar(2,3,n_tht+2) = 2.d0 * R_polar(2,4,n_tht+2) - R_polar(3,3,n_tht+2)  
#endif

if ( write_ps ) then
  call lincol(3)

  npl = 11
  allocate(xout(2),xp(npl),yp(npl))
  do j=1,n_tht+2*n_leg

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
      call lincol(m)
      write(51,*) ' .1 setlinewidth'
      call lplot6(1,1,xp,yp,-npl,' ')
      write(51,*) ' stroke'
      
    enddo
    
    call lincol(0)
    
  enddo
  deallocate(xp,yp,xout)
endif

!----------------------------------- find grid_points from crossing of coordinate lines
RR_new = 0.d0
ZZ_new = 0.d0


do j=1, n_tht          ! the magnetic axis

  RR_new(1,j)    = R_axis
  ZZ_new(1,j)    = Z_axis
  ielm_flux(1,j) = i_elm_axis
  s_flux(1,j)    = s_axis
  t_flux(1,j)    = t_axis
  t_tht(1,j)     = -1.d0          ! expressed in cubic Hermite (-1<t<+1)

enddo


allocate(k_cross(n_flux+n_open+n_private+1,n_theta))

k_cross = 0

k_cross(1,:) = 1

do i=1, n_flux + n_open - 1

  do j=1, n_tht

    do k=1,n_pieces       ! 3 line pieces per coordinate line

      R_cub1d = (/ R_polar(k,1,j), 3.d0/2.d0 *(R_polar(k,2,j)-R_polar(k,1,j)), &
                   R_polar(k,4,j), 3.d0/2.d0 *(R_polar(k,4,j)-R_polar(k,3,j))  /)
      Z_cub1d = (/ Z_polar(k,1,j), 3.d0/2.d0 *(Z_polar(k,2,j)-Z_polar(k,1,j)), &
                   Z_polar(k,4,j), 3.d0/2.d0 *(Z_polar(k,4,j)-Z_polar(k,3,j)) /)

      call find_crossing(node_list,element_list,flux_list,i,R_cub1d,Z_cub1d, &
                       RR_new(i+1,j),ZZ_new(i+1,j),ielm_flux(i+1,j),s_flux(i+1,j),t_flux(i+1,j),t_tht(i+1,j),ifail,.false.)

      call CUB1D(R_cub1d(1), R_cub1d(2), R_cub1d(3), R_cub1d(4),t_tht(i+1,j),tmp1, dR_dt)
      call CUB1D(Z_cub1d(1), Z_cub1d(2), Z_cub1d(3), Z_cub1d(4),t_tht(i+1,j),tmp2, dZ_dt)

      if (ifail .eq. 0) then
        call interp(node_list,element_list,ielm_flux(i+1,j),1,1,s_flux(i+1,j),t_flux(i+1,j),&
                    PSg1,dPSg1_dr,dPSg1_ds,dPSg1_drs,dPSg1_drr,dPSg1_dss)

        k_cross(i+1,j) = k
        exit
      endif

    enddo

  enddo

enddo

!DIR$ NOVECTOR
do i=n_flux-1, n_flux - 1  + n_open + n_private
!DIR$ NOVECTOR
  do j=n_tht+1, n_tht + 2*n_leg
!DIR$ NOVECTOR
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

  enddo

enddo

!***********************************************************************
!*     define the new nodes and finite elements    (nodes first)       *
!***********************************************************************

! --------------------------------------------------------------------------------------
! --------------------------------- Define Final Grid ----------------------------------
! --------------------------------------------------------------------------------------

! --- Allocate data structures for new nodes and elements and initialize them.
allocate(newnode_list)
call tr_register_mem(sizeof(newnode_list),"newnode_list")
allocate(newelement_list)
call tr_register_mem(sizeof(newelement_list),"newelement_list")

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

    k = k_cross(i,j)

    R_cub1d = (/ R_polar(k,1,j), 3.d0/2.d0 *(R_polar(k,2,j)-R_polar(k,1,j)), &
                 R_polar(k,4,j), 3.d0/2.d0 *(R_polar(k,4,j)-R_polar(k,3,j))  /)
    Z_cub1d = (/ Z_polar(k,1,j), 3.d0/2.d0 *(Z_polar(k,2,j)-Z_polar(k,1,j)), &
                 Z_polar(k,4,j), 3.d0/2.d0 *(Z_polar(k,4,j)-Z_polar(k,3,j)) /)

    call CUB1D(R_cub1d(1), R_cub1d(2), R_cub1d(3), R_cub1d(4),t_tht(i,j),tmp1, dR_dt)
    call CUB1D(Z_cub1d(1), Z_cub1d(2), Z_cub1d(3), Z_cub1d(4),t_tht(i,j),tmp2, dZ_dt)

    call interp_RZ(node_list,element_list,ielm_flux(i,j),s_flux(i,j),t_flux(i,j), &
                   RRg1,dRRg1_dr,dRRg1_ds,ZZg1,dZZg1_dr,dZZg1_ds)

    call interp(node_list,element_list,ielm_flux(i,j),1,1,s_flux(i,j),t_flux(i,j),&
                   PSg1,dPSg1_dr,dPSg1_ds,dPSg1_drs,dPSg1_drr,dPSg1_dss)

    RZ_jac  = DRRg1_dr * dZZg1_ds - dRRg1_ds * dZZg1_dr

    PSI_R = (   dPSg1_dr * dZZg1_ds - dPSg1_ds * dZZg1_dr ) / RZ_jac
    PSI_Z = ( - dPSg1_dr * dRRg1_ds + dPSg1_ds * dRRg1_dr ) / RZ_jac

    node  = (n_tht-1)*(i-1) + j
    index = node

    newnode_list%node(index)%x(1,:) = (/ RR_new(i,j), ZZ_new(i,j) /)
    newnode_list%node(index)%x(2,:) = (/ dR_dt, dZ_dt /) / sqrt(dR_dt**2 + dZ_dt**2)
    newnode_list%node(index)%boundary = 0

    if (i .eq. 1) then   !------------------------------------ magnetic axis : special case
      newnode_list%node(index)%x(3,:) = 0.d0
      newnode_list%node(index)%x(4,:) = 0.d0
    else
      newnode_list%node(index)%x(3,:) = (/ -PSI_Z, +PSI_R /) / sqrt(PSI_R**2 + PSI_Z**2)
      newnode_list%node(index)%x(4,:) = 0.d0
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

newnode_list%node(index_xpoint)%x(1,:) = (/ R_xpoint(1), Z_xpoint(1) /)
newnode_list%node(index_xpoint)%x(2,:) = (/ cos(angle_L9),sin(angle_L9) /)
newnode_list%node(index_xpoint)%x(3,:) = (/ alpha_max,1.d0 /) / sqrt(alpha_max**2 + 1.d0)
newnode_list%node(index_xpoint)%x(4,:) = 0.d0
newnode_list%node(index_xpoint)%boundary = 0

newnode_list%node(index_xpoint+1)%x(1,:) = (/ R_xpoint(1), Z_xpoint(1) /)
newnode_list%node(index_xpoint+1)%x(2,:) = (/ R_xpoint(1)-R_axis,Z_xpoint(1)-Z_axis/)/sqrt((R_xpoint(1)-R_axis)**2+(Z_xpoint(1)-Z_axis)**2)
newnode_list%node(index_xpoint+1)%x(3,:) = (/ alpha_max,1.d0 /) / sqrt(alpha_max**2 + 1.d0)
newnode_list%node(index_xpoint+1)%x(4,:) = 0.d0
newnode_list%node(index_xpoint+1)%boundary = 0

newnode_list%node(index_xpoint+2)%x(1,:) = (/ R_xpoint(1), Z_xpoint(1) /)
newnode_list%node(index_xpoint+2)%x(2,:) = (/ R_xpoint(1)-R_axis,Z_xpoint(1)-Z_axis/)/sqrt((R_xpoint(1)-R_axis)**2+(Z_xpoint(1)-Z_axis)**2)
newnode_list%node(index_xpoint+2)%x(3,:) = (/ alpha_min,1.d0 /) / sqrt(alpha_min**2 + 1.d0)
newnode_list%node(index_xpoint+2)%x(4,:) = 0.d0
newnode_list%node(index_xpoint+2)%boundary = 0

newnode_list%node(index_xpoint+3)%x(1,:) = (/ R_xpoint(1), Z_xpoint(1) /)
newnode_list%node(index_xpoint+3)%x(2,:) = (/ -cos(angle_L8),-sin(angle_L8) /)
newnode_list%node(index_xpoint+3)%x(3,:) = (/ alpha_min,1.d0 /) / sqrt(alpha_min**2 + 1.d0)
newnode_list%node(index_xpoint+3)%x(4,:) = 0.d0
newnode_list%node(index_xpoint+3)%boundary = 0

newnode_list%n_nodes = newnode_list%n_nodes + n_xpoint

index = newnode_list%n_nodes

do i=n_flux-1,n_flux-1+n_open           !--------------------------- nodes on the open field lines

  j_start = 1; j_end = n_tht   ! skip first and last point on separatrix (x-points already added)
  if (i.eq. n_flux-1) then
    j_start = 2
    j_end   = n_tht-1
  endif

  do j=j_start, j_end

    k = k_cross(i+1,j)

    R_cub1d = (/ R_polar(k,1,j), 3.d0/2.d0 *(R_polar(k,2,j)-R_polar(k,1,j)), &
                 R_polar(k,4,j), 3.d0/2.d0 *(R_polar(k,4,j)-R_polar(k,3,j))  /)
    Z_cub1d = (/ Z_polar(k,1,j), 3.d0/2.d0 *(Z_polar(k,2,j)-Z_polar(k,1,j)), &
                 Z_polar(k,4,j), 3.d0/2.d0 *(Z_polar(k,4,j)-Z_polar(k,3,j)) /)

    call CUB1D(R_cub1d(1), R_cub1d(2), R_cub1d(3), R_cub1d(4),t_tht(i+1,j),tmp1, dR_dt)
    call CUB1D(Z_cub1d(1), Z_cub1d(2), Z_cub1d(3), Z_cub1d(4),t_tht(i+1,j),tmp2, dZ_dt)

    call interp_RZ(node_list,element_list,ielm_flux(i+1,j),s_flux(i+1,j),t_flux(i+1,j), &
                   RRg1,dRRg1_dr,dRRg1_ds,ZZg1,dZZg1_dr,dZZg1_ds)

    call interp(node_list,element_list,ielm_flux(i+1,j),1,1,s_flux(i+1,j),t_flux(i+1,j),&
                   PSg1,dPSg1_dr,dPSg1_ds,dPSg1_drs,dPSg1_drr,dPSg1_dss)

    RZ_jac  = DRRg1_dr * dZZg1_ds - dRRg1_ds * dZZg1_dr

    PSI_R = (   dPSg1_dr * dZZg1_ds - dPSg1_ds * dZZg1_dr ) / RZ_jac
    PSI_Z = ( - dPSg1_dr * dRRg1_ds + dPSg1_ds * dRRg1_dr ) / RZ_jac

    index = index + 1

    newnode_list%node(index)%x(1,:) = (/ RR_new(i+1,j), ZZ_new(i+1,j) /)
    newnode_list%node(index)%x(2,:) = (/ dR_dt, dZ_dt /)   / sqrt(dR_dt**2 + dZ_dt**2)
    newnode_list%node(index)%x(3,:) = (/ -PSI_Z, +PSI_R /) / sqrt(PSI_R**2 + PSI_Z**2)
    newnode_list%node(index)%x(4,:) = 0.d0
    newnode_list%node(index)%boundary = 0

    if (i .eq. n_flux-1+n_open) then
      newnode_list%node(index)%boundary = 2
    endif
  enddo
enddo

index_ext1 = index

newnode_list%n_nodes = newnode_list%n_nodes + (n_open+1) * n_tht - 2

index = newnode_list%n_nodes

index_leg2 = index+1

do j=1, n_leg                         !--------------------------- nodes on outer leg

  j2 = n_tht + n_leg + j

  do k=1,n_open+n_private+1

    if (k .le. n_private) then
      i = n_flux - 1 + n_open + n_private - k + 1
    else
      i = n_flux - 1 + (k-n_private) - 1
    endif

    if (.not. ( (j .eq. n_leg) .and. (i .le. n_flux+n_open-1))) then  ! exclude horizontal line (already there)

      m = k_cross(i+1,j2)

      R_cub1d = (/ R_polar(m,1,j2), 3.d0/2.d0 *(R_polar(m,2,j2)-R_polar(m,1,j2)), &
                   R_polar(m,4,j2), 3.d0/2.d0 *(R_polar(m,4,j2)-R_polar(m,3,j2))  /)
      Z_cub1d = (/ Z_polar(m,1,j2), 3.d0/2.d0 *(Z_polar(m,2,j2)-Z_polar(m,1,j2)), &
                   Z_polar(m,4,j2), 3.d0/2.d0 *(Z_polar(m,4,j2)-Z_polar(m,3,j2)) /)

      call CUB1D(R_cub1d(1), R_cub1d(2), R_cub1d(3), R_cub1d(4),t_tht(i+1,j2),tmp1, dR_dt)
      call CUB1D(Z_cub1d(1), Z_cub1d(2), Z_cub1d(3), Z_cub1d(4),t_tht(i+1,j2),tmp2, dZ_dt)

      call interp_RZ(node_list,element_list,ielm_flux(i+1,j2),s_flux(i+1,j2),t_flux(i+1,j2), &
                     RRg1,dRRg1_dr,dRRg1_ds,ZZg1,dZZg1_dr,dZZg1_ds)

      call interp(node_list,element_list,ielm_flux(i+1,j2),1,1,s_flux(i+1,j2),t_flux(i+1,j2),&
                     PSg1,dPSg1_dr,dPSg1_ds,dPSg1_drs,dPSg1_drr,dPSg1_dss)

      RZ_jac = DRRg1_dr * dZZg1_ds - dRRg1_ds * dZZg1_dr

      PSI_R = (   dPSg1_dr * dZZg1_ds - dPSg1_ds * dZZg1_dr ) / RZ_jac
      PSI_Z = ( - dPSg1_dr * dRRg1_ds + dPSg1_ds * dRRg1_dr ) / RZ_jac

      index = index + 1

      newnode_list%node(index)%x(1,:) = (/ RR_new(i+1,j2), ZZ_new(i+1,j2) /)
      newnode_list%node(index)%x(2,:) = (/ dR_dt, dZ_dt /)   / sqrt(dR_dt**2 + dZ_dt**2)
      newnode_list%node(index)%x(3,:) = (/ -PSI_Z, +PSI_R /) / sqrt(PSI_R**2 + PSI_Z**2)
      newnode_list%node(index)%x(4,:) = 0.d0
      newnode_list%node(index)%boundary = 0

      if ((k.eq.1) .or. (k .eq. n_open+n_private+1)) then
        newnode_list%node(index)%boundary = newnode_list%node(index)%boundary + 2
      endif
      if (j.eq.1) then
        newnode_list%node(index)%boundary = newnode_list%node(index)%boundary + 1
      endif

   endif

  enddo
enddo

index_leg1 = index+1

do l=1, n_leg-1                       !--------------------------- nodes on inner leg

  j  = n_leg - l
  j2 = n_tht + j

  do k=1,n_open+n_private+1

    if (k .le. n_private) then
      i = n_flux -1 + n_open + n_private - k + 1
    else
      i = n_flux -1 + (k-n_private) - 1
    endif

    m = k_cross(i+1,j2)

    R_cub1d = (/ R_polar(m,1,j2), 3.d0/2.d0 *(R_polar(m,2,j2)-R_polar(m,1,j2)), &
                 R_polar(m,4,j2), 3.d0/2.d0 *(R_polar(m,4,j2)-R_polar(m,3,j2))  /)
    Z_cub1d = (/ Z_polar(m,1,j2), 3.d0/2.d0 *(Z_polar(m,2,j2)-Z_polar(m,1,j2)), &
                 Z_polar(m,4,j2), 3.d0/2.d0 *(Z_polar(m,4,j2)-Z_polar(m,3,j2)) /)

    call CUB1D(R_cub1d(1), R_cub1d(2), R_cub1d(3), R_cub1d(4),t_tht(i+1,j2),tmp1, dR_dt)
    call CUB1D(Z_cub1d(1), Z_cub1d(2), Z_cub1d(3), Z_cub1d(4),t_tht(i+1,j2),tmp2, dZ_dt)

    call interp_RZ(node_list,element_list,ielm_flux(i+1,j2),s_flux(i+1,j2),t_flux(i+1,j2), &
                   RRg1,dRRg1_dr,dRRg1_ds,ZZg1,dZZg1_dr,dZZg1_ds)

    call interp(node_list,element_list,ielm_flux(i+1,j2),1,1,s_flux(i+1,j2),t_flux(i+1,j2),&
                   PSg1,dPSg1_dr,dPSg1_ds,dPSg1_drs,dPSg1_drr,dPSg1_dss)

    RZ_jac  = DRRg1_dr * dZZg1_ds - dRRg1_ds * dZZg1_dr

    PSI_R = (   dPSg1_dr * dZZg1_ds - dPSg1_ds * dZZg1_dr ) / RZ_jac
    PSI_Z = ( - dPSg1_dr * dRRg1_ds + dPSg1_ds * dRRg1_dr ) / RZ_jac

    index = index + 1

    newnode_list%node(index)%x(1,:) = (/ RR_new(i+1,j2), ZZ_new(i+1,j2) /)
    newnode_list%node(index)%x(2,:) = (/ dR_dt, dZ_dt /)   / sqrt(dR_dt**2 + dZ_dt**2)
    newnode_list%node(index)%x(3,:) = (/ -PSI_Z, +PSI_R /) / sqrt(PSI_R**2 + PSI_Z**2)
    newnode_list%node(index)%x(4,:) = 0.d0
    newnode_list%node(index)%boundary = 0

    if ((k.eq.1) .or. (k .eq. n_open+n_private+1)) then
      newnode_list%node(index)%boundary = newnode_list%node(index)%boundary + 2
    endif
    if  (j.eq.1) then
      newnode_list%node(index)%boundary = newnode_list%node(index)%boundary + 1
    endif

 enddo
enddo

newnode_list%n_nodes = index
 
if (extend) then

  allocate(elm_left(n_ext),elm_right(n_ext))

  index_ext2 = index + 1

  sign_psi = 1.d0
  if (psi_axis .gt. psi_xpoint(1)) sign_psi = -1.d0

  do i=1,n_ext
    do j=1,n_tht

      k = k_cross(n_flux+n_open,j)

      if (k .eq. 2) then
        
        t2 = 1.d0 - t_tht(n_flux+n_open,j)
        t3 = 2.d0

        t_total = t2 + t3

        t_delta = t_total * real(i,8)/real(n_ext,8)

        t_node = t_tht(n_flux+n_open,j) + t_delta

        if (t_node .gt. 1.d0) then

          k = 3

          t_node = 1.d0 - (t_delta - t2)
        
        endif

      else

        t_node = -1.d0 + float(n_ext-i)/float(n_ext) * (t_tht(n_flux+n_open,j)+1.d0)

      endif  


      R_cub1d = (/ R_polar(k,1,j), 3.d0/2.d0 *(R_polar(k,2,j)-R_polar(k,1,j)), &
                   R_polar(k,4,j), 3.d0/2.d0 *(R_polar(k,4,j)-R_polar(k,3,j))  /)
      Z_cub1d = (/ Z_polar(k,1,j), 3.d0/2.d0 *(Z_polar(k,2,j)-Z_polar(k,1,j)), &
                   Z_polar(k,4,j), 3.d0/2.d0 *(Z_polar(k,4,j)-Z_polar(k,3,j)) /)


      call CUB1D(R_cub1d(1), R_cub1d(2), R_cub1d(3), R_cub1d(4), t_node, Rtmp, dR_dt)
      call CUB1D(Z_cub1d(1), Z_cub1d(2), Z_cub1d(3), Z_cub1d(4), t_node, Ztmp, dZ_dt)

      tht_bnd = atan2(sign_psi*newnode_list%node(index_ext1-n_tht+j)%x(3,2),&
                      sign_psi*newnode_list%node(index_ext1-n_tht+j)%x(3,1))

      if (T_wall_par(j) .gt. PI)  T_wall_par(j) = T_wall_par(j) - 2.d0*PI
      if (T_wall_par(j) .lt.-PI)  T_wall_par(j) = T_wall_par(j) + 2.d0*PI
      if (tht_bnd - T_wall_par(j) .gt. PI)  tht_bnd = tht_bnd - 2.d0*PI
      if (tht_bnd - T_wall_par(j) .lt.-PI)  tht_bnd = tht_bnd + 2.d0*PI
      if (T_wall_par(j) - tht_bnd .gt. PI)  tht_bnd = tht_bnd + 2.d0*PI
      if (T_wall_par(j) - tht_bnd .lt.-PI)  tht_bnd = tht_bnd - 2.d0*PI

      tht_ext = tht_bnd + (T_wall_par(j) - tht_bnd) * float(i)/float(n_ext)

      index = index+1
      newnode_list%node(index)%x(1,:) = (/ Rtmp, Ztmp /)
      newnode_list%node(index)%x(2,:) = (/ dR_dt, dZ_dt /) / sqrt(dR_dt**2 + dZ_dt**2)
      newnode_list%node(index)%x(3,:) = (/ cos(tht_ext), sin(tht_ext) /)
      newnode_list%node(index)%x(4,:) = 0.d0
      newnode_list%node(index)%boundary = 0

      if (i .eq. n_ext) then
        newnode_list%node(index)%boundary = 5
        newnode_list%node(index_ext1-n_tht+j)%boundary = 0
      endif

    enddo

    if (i.eq.1) then
      index_ext_leg1 = index+1
    endif

    do j=2,n_leg                                  ! left leg
      j2    = n_tht + j 
      j2rev = n_tht + n_leg - j + 1
      index = index+1
      dRtmp = R_wall_max(j2) - R_max(j2rev)
      dZtmp = Z_wall_max(j2) - Z_max(j2rev)

      dRtmp = R_wall_max(j2) - newnode_list%node(index_leg1 + (j-1)*(n_open+n_private+1) - 1)%x(1,1)
      dZtmp = Z_wall_max(j2) - newnode_list%node(index_leg1 + (j-1)*(n_open+n_private+1) - 1)%x(1,2)
      Rtmp = newnode_list%node(index_leg1 + (j-1)*(n_open+n_private+1) - 1)%x(1,1) + dRtmp * float(i)/float(n_ext)
      Ztmp = newnode_list%node(index_leg1 + (j-1)*(n_open+n_private+1) - 1)%x(1,2) + dZtmp * float(i)/float(n_ext)

      tht_bnd = atan2(newnode_list%node(index_leg1+(j-1)*(n_open+n_private+1)-1)%x(3,2),&
                      newnode_list%node(index_leg1+(j-1)*(n_open+n_private+1)-1)%x(3,1))

      if (T_wall_par(j2) .gt. PI)            T_wall_par(j2) = T_wall_par(j2) - 2.d0*PI
      if (T_wall_par(j2) - tht_bnd .gt. PI)  T_wall_par(j2) = T_wall_par(j2) - 2.d0*PI
      if (tht_bnd - T_wall_par(j2) .gt. PI)  tht_bnd = tht_bnd - 2.d0*PI

      tht_ext = tht_bnd + (T_wall_par(j2) - tht_bnd) * float(i)/float(n_ext)

      if (j .gt. 2) then
        dRtmp2 = R_wall_max(j2-1) - newnode_list%node(index_leg1 + (j-2)*(n_open+n_private+1) - 1)%x(1,1)
        dZtmp2 = Z_wall_max(j2-1) - newnode_list%node(index_leg1 + (j-2)*(n_open+n_private+1) - 1)%x(1,2)
        Rtmp2 = newnode_list%node(index_leg1 + (j-2)*(n_open+n_private+1) - 1)%x(1,1) + dRtmp2 * float(i)/float(n_ext)
        Ztmp2 = newnode_list%node(index_leg1 + (j-2)*(n_open+n_private+1) - 1)%x(1,2) + dZtmp2 * float(i)/float(n_ext)
      else
        Rtmp2 = Rtmp
        Ztmp2 = Ztmp
      endif

      if (j .lt. n_leg) then
        dRtmp3 = R_wall_max(j2+1) - newnode_list%node(index_leg1 + (j)*(n_open+n_private+1) - 1)%x(1,1)
        dZtmp3 = Z_wall_max(j2+1) - newnode_list%node(index_leg1 + (j)*(n_open+n_private+1) - 1)%x(1,2)
        Rtmp3 = newnode_list%node(index_leg1 + (j)*(n_open+n_private+1) - 1)%x(1,1) + dRtmp3 * float(i)/float(n_ext)
        Ztmp3 = newnode_list%node(index_leg1 + (j)*(n_open+n_private+1) - 1)%x(1,2) + dZtmp3 * float(i)/float(n_ext)
      else
        Rtmp3 = Rtmp
        Ztmp3 = Ztmp
      endif

      tht_ext = atan2(Ztmp2-Ztmp3,Rtmp2-Rtmp3)

      newnode_list%node(index)%x(1,:) = (/ Rtmp, Ztmp /)
      newnode_list%node(index)%x(2,:) = (/ dRtmp, dZtmp /) / sqrt(dRtmp*dRtmp+dZtmp*dZtmp)
      newnode_list%node(index)%x(3,:) = (/ cos(tht_ext),  sin(tht_ext) /)
      newnode_list%node(index)%x(4,:) = 0.d0
      newnode_list%node(index)%boundary = 0

      if ((dRtmp**2 + dZtmp**2) .lt. 1d-8) then
        newnode_list%node(index)%x(2,:) = (/ -sin(tht_ext),  cos(tht_ext) /)
      endif

      if ((i .eq. n_ext) .and. (j .ne. n_leg)) then
        newnode_list%node(index)%boundary = 5
        newnode_list%node(index_leg1 + (j-1)*(n_open+n_private+1) - 1)%boundary = 0
      endif
      if ((i .eq. n_ext) .and. (j .eq. n_leg)) then
        newnode_list%node(index)%boundary = 9
      endif
    enddo

    if (i.eq.1) then
      index_ext_leg2 = index+1
    endif

    do j=2,n_leg                                 ! right leg
      j2    = n_tht + n_leg   + j
      j2rev = n_tht + 2*n_leg - j + 1 
      index = index+1

      dRtmp = R_wall_max(j2) - newnode_list%node(index_leg2 + (n_leg-j+1)*(n_open+n_private+1) - 1)%x(1,1)
      dZtmp = Z_wall_max(j2) - newnode_list%node(index_leg2 + (n_leg-j+1)*(n_open+n_private+1) - 1)%x(1,2)
      Rtmp = newnode_list%node(index_leg2 + (n_leg-j+1)*(n_open+n_private+1) - 1)%x(1,1) + dRtmp * float(i)/float(n_ext)
      Ztmp = newnode_list%node(index_leg2 + (n_leg-j+1)*(n_open+n_private+1) - 1)%x(1,2) + dZtmp * float(i)/float(n_ext)

      tht_bnd = atan2(newnode_list%node(index_leg2 + (n_leg-j+1)*(n_open+n_private+1) - 1)%x(3,2),&
                      newnode_list%node(index_leg2 + (n_leg-j+1)*(n_open+n_private+1) - 1)%x(3,1))

      if (T_wall_par(j2) .gt. PI)            T_wall_par(j2) = T_wall_par(j2) - 2.d0*PI
      if (T_wall_par(j2) - tht_bnd .gt. PI)  T_wall_par(j2) = T_wall_par(j2) - 2.d0*PI
      if (tht_bnd - T_wall_par(j2) .gt. PI)  tht_bnd = tht_bnd - 2.d0*PI

      tht_ext = tht_bnd + (T_wall_par(j2) - tht_bnd) * float(i)/float(n_ext)

      if (j .gt. 2) then
        dRtmp2 = R_wall_max(j2-1) - newnode_list%node(index_leg2 + (n_leg-j+2)*(n_open+n_private+1) - 1)%x(1,1)
        dZtmp2 = Z_wall_max(j2-1) - newnode_list%node(index_leg2 + (n_leg-j+2)*(n_open+n_private+1) - 1)%x(1,2)
        Rtmp2 = newnode_list%node(index_leg2 + (n_leg-j+2)*(n_open+n_private+1) - 1)%x(1,1) + dRtmp2 * float(i)/float(n_ext)
        Ztmp2 = newnode_list%node(index_leg2 + (n_leg-j+2)*(n_open+n_private+1) - 1)%x(1,2) + dZtmp2 * float(i)/float(n_ext)
      else
        Rtmp2 = Rtmp
        Ztmp2 = Ztmp
      endif
      if (j .lt. n_leg) then
        dRtmp3 = R_wall_max(j2+1) - newnode_list%node(index_leg2 + (n_leg-j)*(n_open+n_private+1) - 1)%x(1,1)
        dZtmp3 = Z_wall_max(j2+1) - newnode_list%node(index_leg2 + (n_leg-j)*(n_open+n_private+1) - 1)%x(1,2)
        Rtmp3 = newnode_list%node(index_leg2 + (n_leg-j)*(n_open+n_private+1) - 1)%x(1,1) + dRtmp3 * float(i)/float(n_ext)
        Ztmp3 = newnode_list%node(index_leg2 + (n_leg-j)*(n_open+n_private+1) - 1)%x(1,2) + dZtmp3 * float(i)/float(n_ext)
      else
        Rtmp3 = Rtmp
        Ztmp3 = Ztmp
      endif

      tht_ext = atan2(Ztmp2-Ztmp3,Rtmp2-Rtmp3)

      newnode_list%node(index)%x(1,:) = (/ Rtmp, Ztmp /)
      newnode_list%node(index)%x(2,:) = (/ dRtmp, dZtmp /) / sqrt(dRtmp*dRtmp+dZtmp*dZtmp)
      newnode_list%node(index)%x(3,:) = (/ cos(tht_ext),  sin(tht_ext) /)
      newnode_list%node(index)%x(4,:) = 0.d0
      newnode_list%node(index)%boundary = 0

      if ((dRtmp**2 + dZtmp**2) .lt. 1d-8) then
        newnode_list%node(index)%x(2,:) = (/ -sin(tht_ext),  cos(tht_ext) /)
      endif

      if ((i .eq. n_ext) .and. (j .ne. n_leg)) then
        newnode_list%node(index)%boundary = 5
        newnode_list%node(index_leg2 + (n_leg-j+1)*(n_open+n_private+1) - 1)%boundary = 0
      endif
      if ((i .eq. n_ext) .and. (j .eq. n_leg)) then
        newnode_list%node(index)%boundary = 9
      endif

    enddo
  enddo

  newnode_list%n_nodes = index
endif

write(*,*) ' definition of nodes completed ',newnode_list%n_nodes

!if ( write_ps ) then
  !call nframe(11,11,1,2.5,3.5,-2.0,-1.0,' ',1,'R',1,'Z',1)
  !call plot_flux_surfaces(node_list,element_list,flux_list,.true.,1)
  !
  !allocate(xp(index),yp(index))
  !do i=1,newnode_list%n_nodes
  !  xp(i) = newnode_list%node(i)%x(1,1)
  !  yp(i) = newnode_list%node(i)%x(1,2)
  !enddo
  !call lplot(1,1,421,xp,yp,-newnode_list%n_nodes,1,'R',1,'Z',1,'nodes',5)
  !deallocate(xp,yp)
!endif

!***********************************************************************
!*                   define the new elements                           *
!***********************************************************************

do i=1,n_flux-1

  do j=1, n_tht-1

    index = (n_tht-1)*(i-1) + j

    newelement_list%element(index)%vertex(1) = (i-1)*(n_tht-1) + j
    newelement_list%element(index)%vertex(2) = (i  )*(n_tht-1) + j
    newelement_list%element(index)%vertex(3) = (i  )*(n_tht-1) + j + 1
    newelement_list%element(index)%vertex(4) = (i-1)*(n_tht-1) + j + 1

    do iv = 1,4
      newelement_list%element(index)%size(iv,:) = (/ 1.d0, 0.01d0, 0.01d0, 0.0d0 /)
    enddo

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
  write(*,*) 'ERROR in grid_xpoint_wall: hard-coded parameter n_nodes_max is too small'
  stop
else if ( newelement_list%n_elements > n_elements_max ) then
  write(*,*) 'ERROR in grid_xpoint_wall: hard-coded parameter n_elements_max is too small'
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

    do iv = 1,4
      newelement_list%element(index)%size(iv,:) = (/ 1.d0, 0.01d0, 0.01d0, 0.0d0 /)
    enddo

  enddo

enddo

newelement_list%n_elements = newelement_list%n_elements + (n_open) * (n_tht-1)

if ( newnode_list%n_nodes > n_nodes_max ) then
  write(*,*) 'ERROR in grid_xpoint_wall: hard-coded parameter n_nodes_max is too small'
  stop
else if ( newelement_list%n_elements > n_elements_max ) then
  write(*,*) 'ERROR in grid_xpoint_wall: hard-coded parameter n_elements_max is too small'
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

    elseif (j .ge. n_leg) then

      newelement_list%element(index)%vertex(1) = node_start + (j-1)*(n_open+n_private+1) + i   - n_open - 1
      newelement_list%element(index)%vertex(2) = node_start + (j-1)*(n_open+n_private+1) + i+1 - n_open - 1
      newelement_list%element(index)%vertex(3) = node_start + (j  )*(n_open+n_private+1) + i+1 - n_open - 1
      newelement_list%element(index)%vertex(4) = node_start + (j  )*(n_open+n_private+1) + i   - n_open - 1

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

    do iv = 1,4
      newelement_list%element(index)%size(iv,:) = (/ 1.d0, 0.01d0, 0.01d0, 0.0d0 /)
    enddo

  enddo

enddo

if (extend) then
  do i=1,n_ext
    do j=1, n_tht-1
      index = index + 1
      if (i.eq.1) then
        newelement_list%element(index)%vertex(1) = index_ext1 - n_tht + j
        newelement_list%element(index)%vertex(4) = index_ext1 - n_tht + j + 1
      else
        newelement_list%element(index)%vertex(1) = index_ext2 + (n_tht+2*(n_leg-1))*(i-2) + j - 1
        newelement_list%element(index)%vertex(4) = index_ext2 + (n_tht+2*(n_leg-1))*(i-2) + j
      endif
      newelement_list%element(index)%vertex(2) = index_ext2 + (n_tht+2*(n_leg-1))*(i-1) + j - 1
      newelement_list%element(index)%vertex(3) = index_ext2 + (n_tht+2*(n_leg-1))*(i-1) + j
    enddo
  enddo

  do i=1, n_ext
    do j=1, n_leg-1                      ! left leg
      index = index + 1
      if (i .eq.1) then
        if (j .eq. 1) then
          newelement_list%element(index)%vertex(1) = index_ext1 
          newelement_list%element(index)%vertex(4) = index_leg1 + (j  )*(n_open+n_private+1) - 1
          newelement_list%element(index)%vertex(2) = index_ext2 + n_tht - 1
          newelement_list%element(index)%vertex(3) = index_ext2 + n_tht 
        else
          newelement_list%element(index)%vertex(1) = index_leg1 + (j-1)*(n_open+n_private+1) - 1 
          newelement_list%element(index)%vertex(4) = index_leg1 + (j  )*(n_open+n_private+1) - 1
          newelement_list%element(index)%vertex(2) = index_ext2 + n_tht - 2 + j
          newelement_list%element(index)%vertex(3) = index_ext2 + n_tht - 1 + j
        endif
      else
        if (j .eq. 1) then
          newelement_list%element(index)%vertex(1) = index_ext2 + n_tht-1 + (i-2)*(n_tht+2*(n_leg-1))
          newelement_list%element(index)%vertex(4) = index_ext2 + n_tht-1 + (i-2)*(n_tht+2*(n_leg-1)) + 1
          newelement_list%element(index)%vertex(2) = index_ext2 + n_tht-1 + (i-1)*(n_tht+2*(n_leg-1))
          newelement_list%element(index)%vertex(3) = index_ext2 + n_tht-1 + (i-1)*(n_tht+2*(n_leg-1)) + 1
        else
          newelement_list%element(index)%vertex(1) = index_ext2 + n_tht-1 + (i-2)*(n_tht+2*(n_leg-1)) + j-1
          newelement_list%element(index)%vertex(4) = index_ext2 + n_tht-1 + (i-2)*(n_tht+2*(n_leg-1)) + j
          newelement_list%element(index)%vertex(2) = index_ext2 + n_tht-1 + (i-1)*(n_tht+2*(n_leg-1)) + j-1
          newelement_list%element(index)%vertex(3) = index_ext2 + n_tht-1 + (i-1)*(n_tht+2*(n_leg-1))+ j
        endif

      endif

      if (i .eq. n_ext) then
        newnode_list%node(newelement_list%element(index)%vertex(1))%boundary = 0
        newnode_list%node(newelement_list%element(index)%vertex(4))%boundary = 0
        newnode_list%node(newelement_list%element(index)%vertex(2))%boundary = 5
        newnode_list%node(newelement_list%element(index)%vertex(3))%boundary = 5
        if (j .eq. n_leg-1)then
          newnode_list%node(newelement_list%element(index)%vertex(3))%boundary = 9 !3  ! LEFT LEG
          newnode_list%node(newelement_list%element(index)%vertex(4))%boundary = 4 !3  ! LEFT LEG
        endif
      endif
      if ((j .eq. n_leg-1) .and. (i .ne. n_ext)) then
          newnode_list%node(newelement_list%element(index)%vertex(4))%boundary = 4 !1  ! LEFT LEG
          newnode_list%node(newelement_list%element(index)%vertex(3))%boundary = 4 !1  ! LEFT LEG
      endif

      if (j .eq. n_leg-1) elm_left(i) = index

    enddo
  enddo

  do i=1,n_ext
    do j=1, n_leg-1                    ! right leg
      index = index + 1
      if (i.eq.1) then
        if (j.eq. 1) then
          newelement_list%element(index)%vertex(4) = index_ext1 - n_tht + 1
          newelement_list%element(index)%vertex(1) = index_leg2 + (n_leg-j)*(n_open+n_private+1) - 1
          newelement_list%element(index)%vertex(3) = index_ext2
          newelement_list%element(index)%vertex(2) = index_ext_leg2 + j - 1 
        else
          newelement_list%element(index)%vertex(4) = index_leg2 + (n_leg-j+1)*(n_open+n_private+1) - 1 
          newelement_list%element(index)%vertex(1) = index_leg2 + (n_leg-j  )*(n_open+n_private+1) - 1
          newelement_list%element(index)%vertex(3) = index_ext_leg2 + j - 2
          newelement_list%element(index)%vertex(2) = index_ext_leg2 + j - 1
        endif
      else
        if (j.eq. 1) then
          newelement_list%element(index)%vertex(1) = index_ext_leg2 + (i-2)*(n_tht+2*(n_leg-1))
          newelement_list%element(index)%vertex(2) = index_ext_leg2 + (i-1)*(n_tht+2*(n_leg-1)) 
          newelement_list%element(index)%vertex(3) = index_ext2 + (n_tht+2*(n_leg-1))*(i-1) + j - 1
          newelement_list%element(index)%vertex(4) = index_ext2 + (n_tht+2*(n_leg-1))*(i-2) + j - 1
        else
          newelement_list%element(index)%vertex(4) = index_ext_leg2 + (i-2)*(n_tht+2*(n_leg-1)) + j - 2 
          newelement_list%element(index)%vertex(3) = index_ext_leg2 + (i-1)*(n_tht+2*(n_leg-1)) + j - 2
          newelement_list%element(index)%vertex(1) = index_ext_leg2 + (i-2)*(n_tht+2*(n_leg-1)) + j - 1
          newelement_list%element(index)%vertex(2) = index_ext_leg2 + (i-1)*(n_tht+2*(n_leg-1)) + j - 1
        endif
      endif

      if (i.eq. n_ext) then
        newnode_list%node(newelement_list%element(index)%vertex(4))%boundary = 0
        newnode_list%node(newelement_list%element(index)%vertex(1))%boundary = 0
        newnode_list%node(newelement_list%element(index)%vertex(3))%boundary = 5
        newnode_list%node(newelement_list%element(index)%vertex(2))%boundary = 5
        if (j .eq. n_leg-1)then
          newnode_list%node(newelement_list%element(index)%vertex(1))%boundary = 4 ! 3  ! RIGHT LEG
          newnode_list%node(newelement_list%element(index)%vertex(2))%boundary = 9 ! 3  ! RIGHT LEG
        endif
      endif
 
      if ((j .eq. n_leg-1) .and. (i .ne. n_ext)) then
          newnode_list%node(newelement_list%element(index)%vertex(1))%boundary = 4  ! 1  ! RIGHT LEG
          newnode_list%node(newelement_list%element(index)%vertex(2))%boundary = 4  ! 1  ! RIGHT LEG
      endif

      if (j .eq. n_leg-1) elm_right(i) = index

    enddo
  enddo

  newelement_list%n_elements = newelement_list%n_elements + n_ext * (n_tht-1) + 2*n_ext*(n_leg-1)
  
  if ( newnode_list%n_nodes > n_nodes_max ) then
    write(*,*) 'ERROR in grid_xpoint_wall: hard-coded parameter n_nodes_max is too small'
    stop
  else if ( newelement_list%n_elements > n_elements_max ) then
    write(*,*) 'ERROR in grid_xpoint_wall: hard-coded parameter n_elements_max is too small'
    stop
  end if

endif

!newelement_list%n_elements = newelement_list%n_elements + (n_open+n_private)*(2*n_leg-2)
newelement_list%n_elements = index

if ( newnode_list%n_nodes > n_nodes_max ) then
  write(*,*) 'ERROR in grid_xpoint_wall: hard-coded parameter n_nodes_max is too small'
  stop
else if ( newelement_list%n_elements > n_elements_max ) then
  write(*,*) 'ERROR in grid_xpoint_wall: hard-coded parameter n_elements_max is too small'
  stop
end if

write(*,*) ' definition of elements completed'

!***********************************************************************
!*             adjust size of elements to get better match             *
!***********************************************************************

do k=1, newelement_list%n_elements   ! fill in the size of the elements

  do iv = 1, 4                    ! over 4 sides of an element

    ivp = mod(iv,4)   + 1         ! vertex with index one higher

    node_iv  = newelement_list%element(k)%vertex(iv)
    node_ivp = newelement_list%element(k)%vertex(ivp) 

    if ((iv .eq. 1) .or. (iv .eq. 3)) then
      R0 = newnode_list%node(node_iv )%X(1,1)  ; dR0 = newnode_list%node(node_iv )%X(2,1)
      Z0 = newnode_list%node(node_iv )%X(1,2)  ; dZ0 = newnode_list%node(node_iv )%X(2,2)
      RP = newnode_list%node(node_ivp)%X(1,1)  ; dRP = newnode_list%node(node_ivp)%X(2,1)
      ZP = newnode_list%node(node_ivp)%X(1,2)  ; dZP = newnode_list%node(node_ivp)%X(2,2)
    else
      R0 = newnode_list%node(node_iv )%X(1,1)  ; dR0 = newnode_list%node(node_iv )%X(3,1)
      Z0 = newnode_list%node(node_iv )%X(1,2)  ; dZ0 = newnode_list%node(node_iv )%X(3,2)
      RP = newnode_list%node(node_ivp)%X(1,1)  ; dRP = newnode_list%node(node_ivp)%X(3,1)
      ZP = newnode_list%node(node_ivp)%X(1,2)  ; dZP = newnode_list%node(node_ivp)%X(3,2)
    endif

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

  newelement_list%element(k)%father     = 0
  newelement_list%element(k)%n_sons     = 0
  newelement_list%element(Index)%sons(:)= 0
enddo


call plot_flux_surfaces(node_list,element_list,flux_list,.true.,1,psi_xpoint,R_xpoint,Z_xpoint,xpoint,xcase)

!***********************************************************************
!*             fill in the values into the new grid                    *
!***********************************************************************

index = 0
do i=1,newnode_list%n_nodes

  newnode_list%node(i)%axis_node = .false.
  if ( fix_axis_nodes .and. (i .le. n_tht) ) newnode_list%node(i)%axis_node = .true.

  do k=1,n_order+1

    index = index + 1
    newnode_list%node(i)%index(k) = index

    if ((force_central_node) .and. (i .gt. 1) .and. (i .le. n_tht) .and. (k.eq.1)) then
      newnode_list%node(i)%index(k) = newnode_list%node(1)%index(1)
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

  R1 = newnode_list%node(i)%x(1,1)
  Z1 = newnode_list%node(i)%x(1,2)

  call find_RZ(node_list,element_list,R1,Z1,R_out,Z_out,ielm_out,s_out,t_out,ifail)

  call interp_RZ(node_list,element_list,ielm_out,s_out,t_out, &
                 RRg1,dRRg1_dr,dRRg1_ds,ZZg1,dZZg1_dr,dZZg1_ds)

  call interp(node_list,element_list,ielm_out,1,1,s_out,t_out,PSg1,dPSg1_dr,dPSg1_ds,dPSg1_drs,dPSg1_drr,dPSg1_dss)

  RZ_jac  = dRRg1_dr * dZZg1_ds - dRRg1_ds * dZZg1_dr
  PSI_R  = (   dZZg1_ds * dPSg1_dr - dZZg1_dr * dPSg1_ds ) / RZ_jac
  PSI_Z  = ( - dRRg1_ds * dPSg1_dr + dRRg1_dr * dPSg1_ds ) / RZ_jac

  newnode_list%node(i)%values(1,1,1) = PSg1
  newnode_list%node(i)%values(1,2,1) = PSI_R * newnode_list%node(i)%x(2,1) + PSI_Z * newnode_list%node(i)%x(2,2)
  newnode_list%node(i)%values(1,3,1) = PSI_R * newnode_list%node(i)%x(3,1) + PSI_Z * newnode_list%node(i)%x(3,2)
  newnode_list%node(i)%values(1,4,1) = PSI_R * newnode_list%node(i)%x(4,1) + PSI_Z * newnode_list%node(i)%x(4,2)

  if (newnode_list%node(i)%boundary .eq. 2) newnode_list%node(i)%values(1,3,1) = 0.d0

enddo

newnode_list%node(index_xpoint  )%values(1,2:4,1) = 0.d0
newnode_list%node(index_xpoint+1)%values(1,2:4,1) = 0.d0
newnode_list%node(index_xpoint+2)%values(1,2:4,1) = 0.d0
newnode_list%node(index_xpoint+3)%values(1,2:4,1) = 0.d0

do j=1,n_tht - 1
  newnode_list%node(i)%values(1,2:4,1) = 0.d0
enddo

n_remove_elements = 0
n_remove_nodes    = 0

do i=1, newelement_list%n_elements
  do j=1, n_vertex_max

    iv = newelement_list%element(i)%vertex(j)

    if (newnode_list%node(iv)%boundary .eq. 9) then  ! remove the small edge "triangles"
      write(*,*) 'removing element : ',i      
      remove_elements(n_remove_elements+1) = i
      n_remove_elements = n_remove_elements + 1

      remove_nodes(n_remove_nodes+1)       = iv
      n_remove_nodes = n_remove_nodes + 1

      iv1 = newelement_list%element(i)%vertex(j)
      iv2 = newelement_list%element(i)%vertex(mod(j,4)+1)
      iv3 = newelement_list%element(i)%vertex(mod(j+1,4)+1)
      iv4 = newelement_list%element(i)%vertex(mod(j+2,4)+1)

      newnode_list%node(iv2)%boundary = 99  ! temporary value to be reset to 9 below
      newnode_list%node(iv3)%boundary = 99  ! temporary value to be reset to 9 below
      newnode_list%node(iv4)%boundary = 99  ! temporary value to be reset to 9 below

    endif
  enddo
enddo
!----------------------------- grid optimisation

!call align_grid(node_List,element_list,new_nodelist,new_element_list,psi_gaussians)


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

!---------------------------- copy new grid into nodes/elements, optional remove some nodes/elements
node_list%n_nodes = 0
skip_index        = 0
do i = 1, newnode_list%n_nodes
  if (.not. any(remove_nodes(1:n_remove_nodes) == i)) then
    newnode_index(i)  = node_list%n_nodes+1 
    node_list%node(node_list%n_nodes+1) = newnode_list%node(i) 
    node_list%node(node_list%n_nodes+1)%index = newnode_list%node(i)%index - skip_index
    node_list%n_nodes = node_list%n_nodes + 1
  else 
    skip_index = skip_index + newnode_list%node(i)%index(4) - newnode_list%node(i)%index(1) + 1
    write(*,*) ' skip_index : ',skip_index
  endif
enddo

element_list%n_elements = 0
do i = 1, newelement_list%n_elements
  if (.not. any(remove_elements(1:n_remove_elements) == i)) then
    element_list%element(element_list%n_elements+1) = newelement_list%element(i)    
    element_list%n_elements                         = element_list%n_elements + 1
  endif
enddo

do i=1, node_list%n_nodes
  if (node_list%node(i)%boundary == 99) node_list%node(i)%boundary = 9  
enddo

do i = 1, element_list%n_elements
  do iv=1, n_vertex_max
    inode = element_list%element(i)%vertex(iv)
    element_list%element(i)%vertex(iv) = newnode_index(inode)  
  enddo
enddo

if ( newnode_list%n_nodes > n_nodes_max ) then
  write(*,*) 'ERROR in grid_xpoint_wall: hard-coded parameter n_nodes_max is too small'
  stop
else if ( element_list%n_elements > n_elements_max ) then
  write(*,*) 'ERROR in grid_xpoint_wall: hard-coded parameter n_elements_max is too small'
  stop
end if

!----temporary, needs to be completed, neighbour information
do i=1, element_list%n_elements
  element_list%element(i)%father = 0
  element_list%element(i)%n_sons = 0
  element_list%element(i)%sons(:) = 0
enddo

deallocate(newnode_list, newelement_list)
deallocate(s_values,theta_sep,R_sep,Z_sep,R_max,Z_max,R_min,Z_min,s_tmp)
deallocate(R_polar,Z_polar)
deallocate(RR_new,ZZ_new,s_flux,t_flux,t_tht)
deallocate(ielm_flux,k_cross)

call update_neighbours(node_list,element_list, force_rtree_initialize=.true.)

write(*,*) ' completed grid_xpoint_wall'

return
end subroutine grid_xpoint_wall
