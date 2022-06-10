subroutine flux_surface_add_line(node_list,element_list,surface_list,i_elm,j,r_psi,s_psi,dpsi_dr,dpsi_ds)
use data_structure
use mod_interp, only: interp

implicit none

type (type_node_list)    :: node_list
type (type_element_list) :: element_list
type (type_surface_list) :: surface_list

integer :: i_elm, j,node1, node2, node3, node4, iter, n_iter, ifail
real*8  :: psi_value, r_psi(*), s_psi(*), dpsi_dr(*), dpsi_ds(*)
real*8  :: rr1, rr2, ss1, ss2,sgn, drs, xl1 ,xl2, drr1, drr2, dss1, dss2, t
real*8  :: ri, si, dri, dsi, delta_ri, delta_si,psi_test, psi_test2, dl1, dl2
real*8  :: psi_r, psi_s, dummy1, dummy2, dummy3, dummy4, dummy5, dummy6, error_psi
real*8  :: drr1_save, dss1_save, drr2_save, dss2_save

rr1  = r_psi(1)
rr2  = r_psi(2)
ss1  = s_psi(1)
ss2  = s_psi(2)

sgn = 1.d0

if ((rr1 .eq. 0.d0) .and. (dpsi_ds(1) .lt. 0.d0)) then
  sgn = -1.d0
elseif ((rr1 .eq. +1.d0).and. (dpsi_ds(1) .gt. 0.d0)) then
  sgn = -1.d0
endif
if ((ss1 .eq. 0.d0) .and. (dpsi_dr(1) .gt. 0.d0)) then
  sgn = -1.d0
elseif ((ss1 .eq. +1.d0).and. (dpsi_dr(1) .lt. 0.d0)) then
  sgn = -1.d0
endif

drs = sgn * sqrt((rr2-rr1)**2 + (ss2-ss1)**2)

xl1 = 0.5d0 * drs / sqrt(dpsi_dr(1)**2 + dpsi_ds(1)**2)  ! temporary fix; waiting for ideas
xl2 = 0.5d0 * drs / sqrt(dpsi_dr(2)**2 + dpsi_ds(2)**2)

drr1 =   xl1 * dpsi_ds(1)
drr2 =   xl2 * dpsi_ds(2)
dss1 = - xl1 * dpsi_dr(1)
dss2 = - xl2 * dpsi_dr(2)

drr1_save = drr1
dss1_save = dss1
drr2_save = drr2
dss2_save = dss2

!------------------------------------------- attempt a correction step to correct drr1,drr2,dss1,dss2
!
ifail  = 0
n_iter = 11

do iter=1, n_iter

  t = 0.d0                 ! make the curve fit exactly at t=0.

  call CUB1D(rr1, drr1, rr2, drr2, t, ri, dri)
  call CUB1D(ss1, dss1, ss2, dss2, t, si, dsi)

  call interp(node_list,element_list,i_elm,1,1,ri,si,psi_test,psi_r,psi_s,dummy1,dummy2,dummy3)

  delta_ri = 0.5d0 * (surface_list%psi_values(j) - psi_test) / psi_r
  delta_si = 0.5d0 * (surface_list%psi_values(j) - psi_test) / psi_s

  call solveM2(0.25d0*drr1,-0.25d0*drr2,0.25d0*dss1,-0.25d0*dss2,delta_ri,delta_si,dl1,dl2)

  dl1 = max(min(dl1,0.2d0),-0.2d0)
  dl2 = max(min(dl2,0.2d0),-0.2d0)

  drr1 = (1.d0+ dl1) * drr1
  dss1 = (1.d0+ dl1) * dss1
  drr2 = (1.d0+ dl2) * drr2
  dss2 = (1.d0+ dl2) * dss2

  call CUB1D(rr1, drr1, rr2, drr2, t, ri, dri)
  call CUB1D(ss1, dss1, ss2, dss2, t, si, dsi)

  call interp(node_list,element_list,i_elm,1,1,ri,si,psi_test2,psi_r,psi_s,dummy1,dummy2,dummy3)

!  write(*,'(A,i3,12e16.8)') ' CHECK IMPROVEMENT : ',iter,dl1,dl2,surface_list%psi_values(j), &
!                         surface_list%psi_values(j)-psi_test,surface_list%psi_values(j)-psi_test2

  error_psi = abs(surface_list%psi_values(j)-psi_test2)/(abs(surface_list%psi_values(j))+abs(psi_test2)+1.d-8)

  if (error_psi .lt. 1.d-9) exit

  if (abs(surface_list%psi_values(j)-psi_test) .lt. abs(surface_list%psi_values(j)-psi_test2)) then
    if ( error_psi .gt.1.d-6) then
!      write(*,'(A,i3,3e16.8)') ' WARNING : IMPROVEMENT GONE WRONG ',iter,surface_list%psi_values(j),psi_test,psi_test2
      ifail = 99
    endif
  endif

!  endif

  if (iter .eq. n_iter) then
!      write(*,'(A,i3,3e16.8)') ' WARNING : IMPROVEMENT FAILED TO CONVERGE ',iter,surface_list%psi_values(j),psi_test,psi_test2
    ifail = 98
  endif

enddo

if (ifail .ne. 0) then
  r_psi(1) = rr1
  r_psi(2) = rr2
  s_psi(1) = ss1
  s_psi(2) = ss2
  drr1 = drr1_save
  dss1 = dss1_save
  drr2 = drr2_save
  dss2 = dss2_save
endif

! --- Derivatives are normally around 0.5 along the flux line
! --- Refinement sometimes goes wrong and give stupid values above 1.0, even up to 3.0 or 4.0
! --- This causes the flux surface line to wrap on itself, if this happens, use the initial guess
! --- Use a threshold of 0.7
if (abs(drr1) .gt. 0.7) drr1 = drr1_save
if (abs(drr2) .gt. 0.7) drr2 = drr2_save
if (abs(dss1) .gt. 0.7) dss1 = dss1_save
if (abs(dss2) .gt. 0.7) dss2 = dss2_save

surface_list%flux_surfaces(j)%n_pieces                                      = surface_list%flux_surfaces(j)%n_pieces + 1
surface_list%flux_surfaces(j)%elm(surface_list%flux_surfaces(j)%n_pieces)   = i_elm
surface_list%flux_surfaces(j)%s(1:4,surface_list%flux_surfaces(j)%n_pieces) = (/ rr1, drr1, rr2, drr2 /)
surface_list%flux_surfaces(j)%t(1:4,surface_list%flux_surfaces(j)%n_pieces) = (/ ss1, dss1, ss2, dss2 /)

!write(*,'(A,i5,4f10.4)') ' adding line : ',i_elm,rr1,ss1,rr2,ss2

return
end
