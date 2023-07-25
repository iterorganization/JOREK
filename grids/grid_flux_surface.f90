subroutine grid_flux_surface(xpoint,xcase,node_list,element_list,surface_list,n_flux,n_tht,xr1,sig1,xr2,sig2,refinement)
!------------------------------------------------------------------------
! subroutine calculates a new flux surface grid (adapted from HELENA20)
!------------------------------------------------------------------------
use constants
use tr_module
use data_structure
use mod_neighbours, only: update_neighbours
use mod_interp
use phys_module, only: force_central_node, fix_axis_nodes, treat_axis
use equil_info
use mod_grid_conversions

implicit none

! --- Routine parameters
logical,                  intent(in)    :: xpoint, refinement
type (type_node_list),    intent(inout) :: node_list
type (type_element_list), intent(inout) :: element_list
type (type_surface_list), intent(inout) :: surface_list
integer,                  intent(in)    :: xcase
integer,                  intent(in)    :: n_flux
integer,                  intent(in)    :: n_tht
real*8,                   intent(in)    :: xr1, xr2
real*8,                   intent(in)    :: sig1, sig2

! --- local variables
integer            :: nrnew, npnew, i, j, k, i_elm
real*8,allocatable :: RRnew(:,:),ZZnew(:,:),PSInew(:,:)
real*8             :: abltg(3), xtmp
real*8,allocatable :: s_values(:),radius(:),psi_values(:),tht_start(:),tht_end(:)
real*8,allocatable :: sp1(:),sp2(:),sp3(:),sp4(:)
real*8             :: dpsi_ds, dpsi_dss, tht_min, tht_max, rr1, rr2, ss1, ss2
real*8             :: RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss
real*8             :: ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss
real*8             :: RRg2,dRRg2_dr,dRRg2_ds,dRRg2_drs,dRRg2_drr,dRRg2_dss
real*8             :: ZZg2,dZZg2_dr,dZZg2_ds,dZZg2_drs,dZZg2_drr,dZZg2_dss
real*8             :: psg1, dpsg1_dr, dpsg1_ds, dpsg1_drs, dpsg1_drr, dpsg1_dss
real*8             :: tht1, tht2, tht_tmp
real*8             :: dRRg1_dt, dZZg1_dt, dRRg2_dt, dZZg2_dt, rz0, rz1, rz2, drz1, drz2
real*8             :: a0, a1, a2, a3
real*8             :: theta, drr1, dss1, drr2, dss2, t, t2, t3, ri, si, dri, dsi, check
real*8             :: rad2, th_z, th_r, th_rr, th_zz, th_rz, dth_ds, dth_dr, dth_drs, dth_drr, dth_dss
real*8             :: rzjac, ps_z, ps_r, ejac, ptjac, rt, st, dptjac_dr, dptjac_ds, rpt, spt, rtt, stt
real*8             :: dr_dr, dr_dz, ds_dr, ds_dz, dps_drr, dps_dzz, crr_axis, czz_axis, cx, cy
real*8             :: dr_dpt, dz_dpt, dr_dtt, dz_dtt, r_ax, s_ax, tn, tn2, cn
real*8             :: delta_rp, delta_zp, delta_rm, delta_zm, dir_2, dir_3, B_axis, q_axis
real*8             :: psi_bnd
real*8, external   :: spwert
integer            :: ifail, inode, node, index, index0, n_node_start, n_element_start, iv, ivp, ivm
integer            :: my_id, n_index_start, node_iv, node_ivp, node_ivm
integer            :: i_sons, n_max

real*8             :: psi_axis_local, R_axis_local, Z_axis_local, s_axis_local, t_axis_local
real*8             :: psi_xpoint_local(2), R_xpoint_local(2), Z_xpoint_local(2), s_xpoint_local(2), t_xpoint_local(2)
integer            :: i_elm_axis_local, i_elm_xpoint_local(2)

write(*,*) '**************************************'
write(*,*) '*         flux surface grid          *'
write(*,*) '**************************************'

my_id = 0
 
call find_axis(my_id,node_list,element_list,psi_axis_local,R_axis_local,Z_axis_local,i_elm_axis_local,s_axis_local,t_axis_local,ifail)  ! left to print some info
 
if (xpoint) call find_xpoint(my_id,node_list,element_list,psi_xpoint_local,R_xpoint_local,Z_xpoint_local,i_elm_xpoint_local,s_xpoint_local,t_xpoint_local,xcase,ifail) !left to print some info

surface_list%n_psi = n_flux - 1
nrnew              = n_flux
npnew              = n_tht

n_max = n_degrees
if (n_order .gt. 5) n_max = (5+1)**2 / 4 ! we don't care about derivatives >= 3...

call tr_allocate(surface_list%psi_values,1,surface_list%n_psi,"surface_list%psi_values",CAT_GRID)
call tr_allocate(psi_values,1,surface_list%n_psi+1,"psi_values",CAT_GRID)
call tr_allocate(s_values,1,nrnew,"s_values",CAT_GRID)
call tr_allocate(radius,1,nrnew,"radius",CAT_GRID)
call tr_allocate(tht_start,1,n_pieces_max,"tht_start",CAT_GRID)
call tr_allocate(tht_end,1,n_pieces_max,"tht_end",CAT_GRID)
call tr_allocate(RRnew,1,n_max,1,nrnew*npnew,"RRnew",CAT_GRID)
call tr_allocate(ZZnew,1,n_max,1,nrnew*npnew,"ZZnew",CAT_GRID)
call tr_allocate(PSInew,1,n_max,1,nrnew*npnew,"PSInew",CAT_GRID)

s_values = 0.d0
call meshac2(surface_list%n_psi+1,s_values,xr1,xr2,sig1,sig2,0.6d0,1.0d0)

psi_values(1) = ES%psi_axis

psi_bnd = ES%psi_bnd + 1.d-8*(ES%psi_axis-ES%psi_bnd)
if (xpoint) then
  psi_bnd = ES%psi_xpoint(1)
  if (ES%active_xpoint .eq. UPPER_XPOINT) then
    psi_bnd = ES%psi_xpoint(2)
  endif
endif

do i=1,surface_list%n_psi
  radius(i+1)                = float(i)/float(surface_list%n_psi)
  surface_list%psi_values(i) = ES%psi_axis + s_values(i+1)**2 *  (psi_bnd - ES%psi_axis)
  psi_values(i+1)            = surface_list%psi_values(i)
!  write(*,'(A,i5,3f14.6)') ' psi values : ',i,radius(i+1),s_values(i+1),surface_list%psi_values(i)
enddo
radius(1)     = 0.d0
psi_values(1) = ES%psi_axis

RRnew(1:n_max,1:nrnew*npnew)  = 0.d0
ZZnew(1:n_max,1:nrnew*npnew)  = 0.d0
PSInew(1:n_max,1:nrnew*npnew) = 0.d0

call find_flux_surfaces(my_id,xpoint,xcase,node_list,element_list,surface_list)
call plot_flux_surfaces(node_list,element_list,surface_list,.true.,1,xpoint,xcase)

!call q_profile(node_list,element_list,surface_list,ES%psi_axis,ES%psi_xpoint,ES%Z_xpoint)

call tr_allocate(sp1,1,surface_list%n_psi+1,"sp1",CAT_GRID)
call tr_allocate(sp2,1,surface_list%n_psi+1,"sp2",CAT_GRID)
call tr_allocate(sp3,1,surface_list%n_psi+1,"sp3",CAT_GRID)
call tr_allocate(sp4,1,surface_list%n_psi+1,"sp4",CAT_GRID)

call spline(surface_list%n_psi+1,radius,s_values,0.d0,0.d0,0,sp1,sp2,sp3,sp4)

do i=1,surface_list%n_psi

    xtmp     = spwert(surface_list%n_psi+1,radius(i+1),sp1,sp2,sp3,sp4,radius,abltg)
    dpsi_ds  = abltg(1) * (psi_bnd - ES%psi_axis) * 2.d0 * s_values(i+1)
    dpsi_dss = abltg(2) * (psi_bnd - ES%psi_axis) * 2.d0 * s_values(i+1)

    tht_min =  1d20
    tht_max = -1d20

    do k=1, surface_list%flux_surfaces(i)%n_pieces

      rr1  = surface_list%flux_surfaces(i)%s(1,k)
      rr2  = surface_list%flux_surfaces(i)%s(3,k)

      ss1  = surface_list%flux_surfaces(i)%t(1,k)
      ss2  = surface_list%flux_surfaces(i)%t(3,k)

      i_elm = surface_list%flux_surfaces(i)%elm(k)

      call interp_RZ(node_list,element_list,i_elm,rr1,ss1,RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss, &
                                                          ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)
      call interp_RZ(node_list,element_list,i_elm,rr2,ss2,RRg2,dRRg2_dr,dRRg2_ds,dRRg2_drs,dRRg2_drr,dRRg2_dss, &
                                                          ZZg2,dZZg2_dr,dZZg2_ds,dZZg2_drs,dZZg2_drr,dZZg2_dss)

      tht1 = atan2(ZZg1-ES%Z_axis,RRg1-ES%R_axis)
      tht2 = atan2(ZZg2-ES%Z_axis,RRg2-ES%R_axis)

      if (tht1 .lt. 0.d0) tht1 = tht1 + 2.d0*PI
      if (tht2 .lt. 0.d0) tht2 = tht2 + 2.d0*PI

      tht_start(k) = min(tht1,tht2)
      tht_end(k)   = max(tht1,tht2)

      if ((tht_end(k) - tht_start(k)) .gt. 3.d0*PI/4.d0) then
         tht_tmp      = tht_end(k)
         tht_end(k)   = tht_start(k)
         tht_start(k) = tht_tmp - 2.d0*PI
      endif

      tht_min = min(tht_min,tht_start(k))
      tht_max = max(tht_max,tht_end(k))

    enddo

    do j=1, npnew

      theta = 2.d0 * PI * float(j-1)/float(npnew)

      if (theta .gt. tht_max) theta = theta - 2.d0*PI
      if (theta .lt. tht_min) theta = theta + 2.d0*PI

      do k=1, surface_list%flux_surfaces(i)%n_pieces

        if ( (theta .ge. tht_start(k)) .and. (theta .le. tht_end(k)) ) then

          rr1  = surface_list%flux_surfaces(i)%s(1,k);   ss1  = surface_list%flux_surfaces(i)%t(1,k)
          drr1 = surface_list%flux_surfaces(i)%s(2,k);   dss1 = surface_list%flux_surfaces(i)%t(2,k)
          rr2  = surface_list%flux_surfaces(i)%s(3,k);   ss2  = surface_list%flux_surfaces(i)%t(3,k)
          drr2 = surface_list%flux_surfaces(i)%s(4,k);   dss2 = surface_list%flux_surfaces(i)%t(4,k)

          i_elm = surface_list%flux_surfaces(i)%elm(k)

          call interp_RZ(node_list,element_list,i_elm,rr1,ss1,RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss, &
                                                              ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)
          call interp_RZ(node_list,element_list,i_elm,rr2,ss2,RRg2,dRRg2_dr,dRRg2_ds,dRRg2_drs,dRRg2_drr,dRRg2_dss, &
                                                              ZZg2,dZZg2_dr,dZZg2_ds,dZZg2_drs,dZZg2_drr,dZZg2_dss)
          dRRg1_dt = dRRg1_dr * drr1 + dRRg1_ds * dss1
          dZZg1_dt = dZZg1_dr * drr1 + dZZg1_ds * dss1
          dRRg2_dt = dRRg2_dr * drr2 + dRRg2_ds * dss2
          dZZg2_dt = dZZg2_dr * drr2 + dZZg2_ds * dss2

          RZ1  = RRg1     * tan(theta) - ZZg1
          RZ2  = RRg2     * tan(theta) - ZZg2
          dRZ1 = dRRg1_dt * tan(theta) - dZZg1_dt
          dRZ2 = dRRg2_dt * tan(theta) - dZZg2_dt

          RZ0  = ES%R_axis  * tan(theta) - ES%Z_axis

          a3 = (   RZ1 + dRZ1 -   RZ2 + dRZ2 )/4.d0
          a2 = (       - dRZ1         + dRZ2 )/4.d0
          a1 = (-3.d0*RZ1 - dRZ1 + 3.d0*RZ2 - dRZ2 )/4.d0
          a0 = ( 2.d0*RZ1 + dRZ1 + 2.d0*RZ2 - dRZ2 )/4.d0 - RZ0

          call SOLVP3(a0,a1,a2,a3,t,t2,t3,ifail)

          if (abs(t) .le. 1.d0 + 1.d-6) then

            call CUB1D(rr1, drr1, rr2, drr2, t, ri, dri)
            call CUB1D(ss1, dss1, ss2, dss2, t, si, dsi)

            call interp_RZ(node_list,element_list,i_elm,ri,si,RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss, &
                                                              ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)
            call interp_RZ(node_list,element_list,i_elm,ri,si,RRg2,dRRg2_dr,dRRg2_ds,dRRg2_drs,dRRg2_drr,dRRg2_dss, &
                                                              ZZg2,dZZg2_dr,dZZg2_ds,dZZg2_drs,dZZg2_drr,dZZg2_dss)
            call interp(node_list,element_list,i_elm,1,1,ri,si,PSg1,dPSg1_dr,dPSg1_ds,dPSg1_drs,dPSg1_drr,dPSg1_dss)

            check = atan2(ZZg1- ES%Z_axis,RRg1-ES%R_axis)
            if (check .lt. 0.d0) check = check + 2.d0*PI

            RRnew(1,npnew*(i) + j)  = RRg1
            ZZnew(1,npnew*(i) + j)  = ZZg1
            PSInew(1,npnew*(i) + j) = surface_list%psi_values(i)

            RAD2    =   (RRg1 - ES%R_axis)**2 + (ZZg1 - ES%Z_axis)**2
            TH_Z    =   (RRg1 - ES%R_axis) / RAD2
            TH_R    = - (ZZg1 - ES%Z_axis) / RAD2

            TH_RR   = 2.d0*(ZZg1 - ES%Z_axis) * (RRg1 - ES%R_axis) / RAD2**2
            TH_ZZ   = - TH_RR
            TH_RZ   = ( (ZZg1 - ES%Z_axis)**2 - (RRg1 - ES%R_axis)**2 ) / RAD2**2

            dTH_ds  = TH_R * dRRg1_ds + TH_Z * dZZg1_ds
            dTH_dr  = TH_R * dRRg1_dr + TH_Z * dZZg1_dr

            dTH_drr = TH_RR * dRRg1_dr * dRRg1_dr + TH_R * dRRg1_drr + 2.d0* TH_RZ * dRRg1_dr * dZZg1_dr &
                    + TH_ZZ * dZZg1_dr * dZZg1_dr + TH_Z * dZZg1_drr

            dTH_drs = TH_RR * dRRg1_dr * dRRg1_ds + TH_RZ * dRRg1_dr * dZZg1_ds + TH_R * dRRg1_drs &
                    + TH_RZ * dZZg1_dr * dRRg1_ds + TH_ZZ * dZZg1_dr * dZZg1_ds + TH_Z * dZZg1_drs

            dTH_dss = TH_RR * dRRg1_ds * dRRg1_ds + TH_R * dRRg1_dss + 2.d0* TH_RZ * dRRg1_ds * dZZg1_ds &
                    + TH_ZZ * dZZg1_ds * dZZg1_ds + TH_Z * dZZg1_dss

            RZjac   = dRRg1_dr * dZZg1_ds - dRRg1_ds * dZZg1_dr

            PS_R    = (  dZZg1_ds * dPSg1_dr - dZZg1_dr * dPSg1_ds) / RZjac
            PS_Z    = (- dRRg1_ds * dPSg1_dr + dRRg1_dr * dPSg1_ds) / RZjac

            Ejac    =  (PS_R * TH_Z - PS_Z * TH_R)
            PTjac   = (dPSg1_dr * dTH_ds - dPSg1_ds * dTH_dr)

            dPTjac_dr = dPSg1_drr * dTH_ds + dPSg1_dr * dTH_drs - dPSg1_drs * dTH_dr - dPSg1_ds * dTH_drr
            dPTjac_ds = dPSg1_drs * dTH_ds + dPSg1_dr * dTH_dss - dPSg1_dss * dTH_dr - dPSg1_ds * dTH_drs

            RT      = - dPSg1_ds / PTjac
            ST      =   dPSg1_dr / PTjac

            RPT = (-dPTjac_dr * dTH_ds / PTjac**2 + dTH_drs/PTjac) * RT + (- dPTjac_ds * dTH_ds / PTjac**2 + dTH_dss/PTjac)*ST
            SPT = ( dPTjac_dr * dTH_dr / PTjac**2 - dTH_drr/PTjac) * RT + (  dPTjac_ds * dTH_dr / PTjac**2 - dTH_drs/PTjac)*ST

            RTT = ( dPTjac_dr * dPSg1_ds / PTjac**2 - dPSg1_drs/PTjac) * RT + (  dPTjac_ds * dPSg1_ds / PTjac**2 - dPSg1_dss/PTjac)*ST
            STT = (-dPTjac_dr * dPSg1_dr / PTjac**2 + dPSg1_drr/PTjac) * RT + (- dPTjac_ds * dPSg1_dr / PTjac**2 + dPSg1_drs/PTjac)*ST

            ! 2nd order psi-tht derivative
            DR_dpt = - dRRg1_drr * dTH_ds * dPSg1_ds / PTjac**2 + dRRg1_drs * (dTH_ds * dPSg1_dr + dTH_dr * dPSg1_ds)/PTjac**2 &
                     + dRRg1_dr  * RPT    + dRRg1_ds * SPT - dRRg1_dss * dTH_dr * dPSg1_dr / PTjac**2

            DZ_dpt = - dZZg1_drr * dTH_ds * dPSg1_ds / PTjac**2 + dZZg1_drs * (dTH_ds * dPSg1_dr + dTH_dr * dPSg1_ds)/PTjac**2 &
                     + dZZg1_dr  * RPT    + dZZg1_ds * SPT - dZZg1_dss * dTH_dr * dPSg1_dr / PTjac**2

            ! 2nd order tht-tht derivative
            DR_dtt = + dRRg1_drr * dPSg1_ds**2 / PTjac**2 - dRRg1_drs * 2.d0 * dPSg1_ds * dPSg1_dr /PTjac**2 &
                     - dRRg1_dr  * RTT    - dRRg1_ds * STT + dRRg1_dss * dPSg1_dr**2 / PTjac**2

            DZ_dtt = + dZZg1_drr * dPSg1_ds**2 / PTjac**2 - dZZg1_drs * 2.d0 * dPSg1_ds * dPSg1_dr /PTjac**2 &
                     - dZZg1_dr  * RTT    - dZZg1_ds * STT + dZZg1_dss * dPSg1_dr**2 / PTjac**2

            RRnew(1,npnew*(i) + j)  = RRg1
            ZZnew(1,npnew*(i) + j)  = ZZg1

            RRnew(2,npnew*(i) + j) = (  TH_Z / Ejac) /  (2.d0*(nrnew-1)) * dpsi_ds
            ZZnew(2,npnew*(i) + j) = (- TH_R / Ejac) /  (2.d0*(nrnew-1)) * dpsi_ds

            RRnew(3,npnew*(i) + j) = + (- PS_Z / Ejac) /  (npnew/PI)
            ZZnew(3,npnew*(i) + j) = + (  PS_R / Ejac) /  (npnew/PI)

            RRnew(4,npnew*(i) + j) =  dR_dpt /(2.d0*(nrnew-1)*(npnew)/PI)  * dpsi_ds
            ZZnew(4,npnew*(i) + j) =  dZ_dpt /(2.d0*(nrnew-1)*(npnew)/PI)  * dpsi_ds

            PSInew(1,npnew*(i) + j) = PSg1
            PSInew(2,npnew*(i) + j) = dpsi_ds / (2.d0*(nrnew-1))
            PSInew(3,npnew*(i) + j) = 0.d0
            PSInew(4,npnew*(i) + j) = 0.d0
            
            if (n_order .ge. 5) then
              RRnew (5,npnew*(i) + j) = 0.d0
              ZZnew (5,npnew*(i) + j) = 0.d0
              RRnew (6,npnew*(i) + j) = dR_dtt / (npnew/PI)**2
              ZZnew (6,npnew*(i) + j) = dZ_dtt / (npnew/PI)**2
              RRnew (7:n_max,npnew*(i) + j) = 0.d0
              ZZnew (7:n_max,npnew*(i) + j) = 0.d0
              PSInew(5,npnew*(i) + j) = dpsi_dss / (2.d0*(nrnew-1))**2
              PSInew(6:n_max,npnew*(i) + j) = 0.d0
            endif

          else

            write(*,*) ' T TOO BIG : ',T,T2,T3

          endif

        endif

      enddo

    enddo

enddo

!----------------------------------- magnetic axis

r_ax = ES%s_axis
s_ax = ES%t_axis

call interp_RZ(node_list,element_list,ES%i_elm_axis,r_ax,s_ax, RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss, &
                                                            ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)
call interp(node_list,element_list,ES%i_elm_axis,1,1,r_ax,s_ax,PSg1,dPSg1_dr,dPSg1_ds,dPSg1_drs,dPSg1_drr,dPSg1_dss)

ejac  = dRRg1_dr * dZZg1_ds - dRRg1_ds * DZZg1_dr
dr_dZ = - dRRg1_ds / ejac
dr_dR = + dZZg1_ds / ejac
ds_dZ = + dRRg1_dr / ejac
ds_dR = - dZZg1_dr / ejac

dPS_dRR = dPSg1_drr * dr_dR * dr_dR + 2.d0*dPSg1_drs * dr_dR * ds_dR + dPSg1_dss * ds_dR * ds_dR
dPS_dZZ = dPSg1_drr * dr_dZ * dr_dZ + 2.d0*dPSg1_drs * dr_dZ * ds_dZ + dPSg1_dss * ds_dZ * ds_dZ

CRR_axis = dPS_dRR / 2.d0 / (ES%psi_bnd-ES%psi_axis)
CZZ_axis = dPS_dZZ / 2.d0 / (ES%psi_bnd-ES%psi_axis)

B_axis = 1.d0 / ES%R_axis
q_axis = B_axis   /(2.d0*SQRT(CRR_axis*CZZ_axis)) / abs(ES%psi_axis)

write(*,'(A,4f14.8)') ' magnetic axis, q : ',ES%R_axis,ES%Z_axis,ES%psi_axis,q_axis

CX = CRR_axis
CY = CZZ_axis

do j=1,npnew

  inode = j
  RRnew(1,inode) = RRg1
  ZZnew(1,inode) = ZZg1

  theta = float(j-1)/float(npnew) * 2.d0*PI

  TN  = TAN(theta)
  TN2 = TN**2
  CN  = COS(theta)

  if (theta .eq. PI/2.d0) then
    RRnew(2,inode) = 0.d0
    ZZnew(2,inode) = +1.d0/(sqrt(abs(CY))*2.d0*float(nrnew-1))
    RRnew(4,inode) = -1.d0/(sqrt(abs(CY))*2.d0*float(nrnew-1)*float(npnew)/PI)
    ZZnew(4,inode) = 0.d0
  ELSEIF (theta .eq. (3.d0*PI/2.d0)) THEN
    RRnew(2,inode) = 0.d0
    ZZnew(2,inode) = +1.d0/(sqrt(abs(CY))*2.d0*float(nrnew-1))
    RRnew(4,inode) = -1.d0/(sqrt(abs(CY))*2.d0*float(nrnew-1)*float(npnew)/PI)
    ZZnew(4,inode) = 0.d0
  ELSE
    RRnew(2,inode) = + sign(1.d0,CN)/(sqrt(abs(CX+CY*TN2))*2.d0 *float(nrnew-1))
    ZZnew(2,inode) = + abs(TN)/(sqrt(abs(CX+CY*TN2))*2.d0*float(nrnew-1))
    RRnew(4,inode) = - abs(CX+CY*TN2)**(-1.5d0) * CY * abs(TN) / (CN**2 * 2.d0*float(nrnew-1)*float(npnew)/PI)
    ZZnew(4,inode) = + CX * abs(CX + CY*TN2)**(-1.5d0) / (CN*abs(CN) * 2.d0*float(nrnew-1)*float(npnew-1)/PI)
  ENDIF

  IF (theta .gt. PI) THEN
    ZZnew(2,inode) = - ZZnew(2,inode)
    RRnew(4,inode) = - RRnew(4,inode)
  ENDIF
  RRnew(3,inode) = 0.d0
  ZZnew(3,inode) = 0.d0
  PSInew(1,inode) = ES%psi_axis

  RRnew(2,inode) = RRnew(2,inode) * sp2(1)
  RRnew(4,inode) = RRnew(4,inode) * sp2(1)
  ZZnew(2,inode) = ZZnew(2,inode) * sp2(1)
  ZZnew(4,inode) = ZZnew(4,inode) * sp2(1)
  PSInew(2,inode) = PSInew(2,inode) * sp2(1)
  PSInew(4,inode) = PSInew(4,inode) * sp2(1)

  if (n_order .ge. 5) then
    RRnew (5:n_max,inode) = 0.d0
    ZZnew (5:n_max,inode) = 0.d0
    PSInew(5:n_max,inode) = 0.d0
  endif

enddo


!----------------------------- empty old nodes/elements

do i=1,n_nodes_max
  node_list%node(i)%x        = 0.d0
  node_list%node(i)%values   = 0.d0
  node_list%node(i)%index    = 0
  node_list%node(i)%boundary = 0
enddo
node_list%n_nodes = 0

do i=1,n_elements_max
  element_list%element(i)%vertex     = 0
  element_list%element(i)%size       = 0.d0
  element_list%element(i)%neighbours = 0
enddo

!----------------------------- convert from cubic Hermite grid to Bezier grid

element_list%n_elements = (nrnew-1)*npnew
node_list%n_nodes       = nrnew*npnew

if ( node_list%n_nodes > n_nodes_max ) then
  write(*,*) 'ERROR in grid_flux_surface: hard-coded parameter n_nodes_max is too small'
  stop
else if ( element_list%n_elements > n_elements_max ) then
  write(*,*) 'ERROR in grid_flux_surface: hard-coded parameter n_elements_max is too small'
  stop
end if

n_node_start    = 0
n_element_start = 0
n_index_start   = 0

do i=1,nrnew-1

  do j=1,npnew-1
    node  = npnew*(i-1) + j
    index = node
    element_list%element(index)%vertex(1) = (i-1)*npnew + j
    element_list%element(index)%vertex(4) = (i-1)*npnew + j + 1
    element_list%element(index)%vertex(3) = (i  )*npnew + j + 1
    element_list%element(index)%vertex(2) = (i  )*npnew + j
  enddo

  index = npnew*(i-1) + npnew

  element_list%element(index)%vertex(1)  = (i  )*npnew
  element_list%element(index)%vertex(4)  = (i  )*npnew - npnew + 1
  element_list%element(index)%vertex(3)  = (i  )*npnew + 1
  element_list%element(index)%vertex(2)  = (i  )*npnew + npnew

enddo

do i=1,nrnew

  do j=1,npnew

    index0 =                npnew*(i-1) + j
    index  = n_node_start + npnew*(i-1) + j

    node_list%node(index)%X(1,:,1)      = 0.d0
    node_list%node(index)%X(1,:,2)      = 0.d0
    node_list%node(index)%values(1,:,1) = 0.d0

    node_list%node(index)%X(1,1,1)      = RRnew (1,index0)
    node_list%node(index)%X(1,1,2)      = ZZnew (1,index0)
    node_list%node(index)%values(1,1,1) = PSInew(1,index0)

    node_list%node(index)%X(1,2,1)      = RRnew (2,index0) * 2.d0/float(n_order)
    node_list%node(index)%X(1,2,2)      = ZZnew (2,index0) * 2.d0/float(n_order)
    node_list%node(index)%values(1,2,1) = PSInew(2,index0) * 2.d0/float(n_order)

    node_list%node(index)%X(1,3,1)      = RRnew (3,index0) * 2.d0/float(n_order)
    node_list%node(index)%X(1,3,2)      = ZZnew (3,index0) * 2.d0/float(n_order)
    node_list%node(index)%values(1,3,1) = PSInew(3,index0) * 2.d0/float(n_order)

    node_list%node(index)%X(1,4,1)      = RRnew (4,index0) * 4.d0/float(n_order)**2
    node_list%node(index)%X(1,4,2)      = ZZnew (4,index0) * 4.d0/float(n_order)**2
    node_list%node(index)%values(1,4,1) = PSInew(4,index0) * 4.d0/float(n_order)**2

    if (n_order .ge. 5) then
      node_list%node(index)%X(1,5,1)      = RRnew (5,index0)  * 4.d0/float(n_order)**2
      node_list%node(index)%X(1,6,1)      = RRnew (6,index0)  * 4.d0/float(n_order)**2
      node_list%node(index)%X(1,7,1)      = RRnew (7,index0)  * 8.d0/float(n_order)**3
      node_list%node(index)%X(1,8,1)      = RRnew (8,index0)  * 8.d0/float(n_order)**3
      node_list%node(index)%X(1,9,1)      = RRnew (9,index0)  *16.d0/float(n_order)**4
      
      node_list%node(index)%X(1,5,2)      = ZZnew (5,index0)  * 4.d0/float(n_order)**2
      node_list%node(index)%X(1,6,2)      = ZZnew (6,index0)  * 4.d0/float(n_order)**2
      node_list%node(index)%X(1,7,2)      = ZZnew (7,index0)  * 8.d0/float(n_order)**3
      node_list%node(index)%X(1,8,2)      = ZZnew (8,index0)  * 8.d0/float(n_order)**3
      node_list%node(index)%X(1,9,2)      = ZZnew (9,index0)  *16.d0/float(n_order)**4
      
      node_list%node(index)%values(1,5,1) = PSInew(5,index0)  * 4.d0/float(n_order)**2
      node_list%node(index)%values(1,6,1) = PSInew(6,index0)  * 4.d0/float(n_order)**2
      node_list%node(index)%values(1,7,1) = PSInew(7,index0)  * 8.d0/float(n_order)**3
      node_list%node(index)%values(1,8,1) = PSInew(8,index0)  * 8.d0/float(n_order)**3
      node_list%node(index)%values(1,9,1) = PSInew(9,index0)  *16.d0/float(n_order)**4
    endif

    if (i .eq. nrnew) node_list%node(index)%boundary = 2

    node_list%node(index)%axis_node = .false.
    node_list%node(index)%axis_dof  = 0    

    if (i .eq. 1) node_list%node(index)%axis_node = .true.

    if (.not. refinement) then       ! keep original formulation if not using refinement
   
      ! Share 4 degrees of freedom for all nodes on the grid axis and flag the axis nodes.
      if(treat_axis)then

         if(i.eq.1)then
           node_list%node(index)%index(1) = 1
           node_list%node(index)%index(2) = 2
           node_list%node(index)%index(3) = 3
           node_list%node(index)%index(4) = 4
           n_index_start = 4
           node_list%node(index)%axis_node = .true.
           node_list%node(index)%axis_dof  = 3

         else
           do k=1,n_degrees
             node_list%node(index)%index(k) = n_index_start + k
           enddo
           n_index_start = n_index_start + n_degrees
           node_list%node(index)%axis_node = .false.
         endif

      elseif (force_central_node .and. (i.eq.1)) then

        node_list%node(index)%index(1) = 1

        if (j.eq.1) n_index_start = n_index_start + 1

        do k=2,n_degrees
          node_list%node(index)%index(k) = n_index_start + k-1
        enddo
        n_index_start = n_index_start +n_degrees-1

      else
        do k=1,n_degrees
          node_list%node(index)%index(k) = n_index_start + k
        enddo
        n_index_start = n_index_start + n_degrees
      endif
   
    else      ! in case of refinement
  
      do k=1,n_degrees
        node_list%node(index)%index(k) = n_index_start + n_degrees*(index0-1)+k
      enddo
  
          !Neighbours of the element (for refinement procedure)

      if(i==1) then
        element_list%element(Index)%neighbours(4) = 0    
      else
        element_list%element(Index)%neighbours(4) = Index - npnew 
      end if 
    
      if(j==npnew) then
        element_list%element(Index)%neighbours(3) = Index - npnew + 1  
      else
        element_list%element(Index)%neighbours(3) = Index + 1       
      end if
    
      if(i==nrnew-1) then
        element_list%element(Index)%neighbours(2) = 0   
      else   
        element_list%element(Index)%neighbours(2) = Index + npnew 
      end if 
        
      if(j==1) then
        element_list%element(Index)%neighbours(1) = Index + npnew -1  
      else   
        element_list%element(Index)%neighbours(1) = Index -1      
      end if  
      
    endif
  
! Initialization of the genealogy  (for refinement procedure)
  
    element_list%element(Index)%father = 0
    element_list%element(Index)%n_sons = 0
    element_list%element(Index)%sons(:) = 0

  enddo
enddo


do k=n_element_start+1 , element_list%n_elements   ! fill in the size of the elements

 do iv = 1, 4                    ! over 4 corners of an element

   ivp = mod(iv,4)   + 1         ! vertex with index one higher
   ivm = mod(iv+2,4) + 1         ! vertex with index one below

   node_iv  = element_list%element(k)%vertex(iv)
   node_ivp = element_list%element(k)%vertex(ivp)
   node_ivm = element_list%element(k)%vertex(ivm)

   if ((iv .eq. 1) .or. (iv .eq.3)) then

     delta_Rp = node_list%node(node_ivp)%X(1,1,1) - node_list%node(node_iv)%X(1,1,1)
     delta_Zp = node_list%node(node_ivp)%X(1,1,2) - node_list%node(node_iv)%X(1,1,2)
     dir_2    = delta_Rp * node_list%node(node_iv)%X(1,2,1) + delta_Zp * node_list%node(node_iv)%X(1,2,2)

     delta_Rm = node_list%node(node_ivm)%X(1,1,1) - node_list%node(node_iv)%X(1,1,1)
     delta_Zm = node_list%node(node_ivm)%X(1,1,2) - node_list%node(node_iv)%X(1,1,2)
     dir_3    = delta_Rm * node_list%node(node_iv)%X(1,3,1) + delta_Zm * node_list%node(node_iv)%X(1,3,2)

   else

     delta_Rp = node_list%node(node_ivp)%X(1,1,1) - node_list%node(node_iv)%X(1,1,1)
     delta_Zp = node_list%node(node_ivp)%X(1,1,2) - node_list%node(node_iv)%X(1,1,2)
     dir_3    = delta_Rp * node_list%node(node_iv)%X(1,3,1) + delta_Zp * node_list%node(node_iv)%X(1,3,2)

     delta_Rm = node_list%node(node_ivm)%X(1,1,1) - node_list%node(node_iv)%X(1,1,1)
     delta_Zm = node_list%node(node_ivm)%X(1,1,2) - node_list%node(node_iv)%X(1,1,2)
     dir_2    = delta_Rm * node_list%node(node_iv)%X(1,2,1) + delta_Zm * node_list%node(node_iv)%X(1,2,2)

   endif

   if (dir_2 .ne. 0.d0) then
     dir_2 = dir_2 / abs(dir_2)
   else
     dir_2 = 1.d0
   endif
   if (dir_3 .ne. 0.d0) then
     dir_3 = dir_3 / abs(dir_3)
   else
     dir_3 = -1.d0
     if (iv.eq.1) dir_3 = 1.d0              ! admittedly not very elegant
   endif

   element_list%element(k)%size(iv,1) = 1.d0
   element_list%element(k)%size(iv,2) = dir_2
   element_list%element(k)%size(iv,3) = dir_3
   element_list%element(k)%size(iv,4) = element_list%element(k)%size(iv,2) * element_list%element(k)%size(iv,3)
   if (fix_axis_nodes) then
      j = element_list%element(k)%vertex(iv)
      if (node_list%node(j)%axis_node) then
        element_list%element(k)%size(iv,3) = 0.d0
        element_list%element(k)%size(iv,4) = 0.d0
      endif
   endif

!   if ((RR(2,node_iv)**2 + ZZ(2,node_iv)**2) .eq. 0.) element_list%element(k)%size(iv,2) = dir_2
!   if ((RR(3,node_iv)**2 + ZZ(3,node_iv)**2) .eq. 0.) element_list%element(k)%size(iv,3) = dir_3
!   if ((RR(4,node_iv)**2 + ZZ(4,node_iv)**2) .eq. 0.) element_list%element(k)%size(iv,4) = dir_2 * dir_3

!    write(*,'(2i5,12e16.8)') k,iv,element_list%element(k)%size(iv,1:4)

 enddo

enddo

if (n_order .ge. 5) call set_high_order_sizes(element_list)
if ( (n_order .ge. 5) .and. fix_axis_nodes) call set_high_order_sizes_on_axis(node_list,element_list)

call update_neighbours(node_list,element_list, force_rtree_initialize=.true.)
return
end subroutine grid_flux_surface
