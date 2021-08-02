!> Export plasma boundary and profiles for later use in HELENA
subroutine export_helena(node_list,element_list,bnd_elm_list)

use constants
use tr_module 
use mod_parameters
use data_structure
use phys_module
use mod_interp
use equil_info

implicit none

! --- Routine parameters
type (type_node_list),        intent(in)    :: node_list
type (type_element_list),     intent(in)    :: element_list
type (type_bnd_element_list), intent(in)    :: bnd_elm_list

! --- Local variables
type (type_surface_list)     :: surface_list
real*8, allocatable      :: rplot(:), zplot(:)
real*8 :: psi_axis, R_axis, Z_axis, s_axis, t_axis, psi_xpoint(2), R_xpoint(2), Z_xpoint(2), s_xpoint(2), t_xpoint(2)
real*8 :: psi_lim,R_lim,Z_lim
real*8 :: psi_bnd, Rmin, Rmax, rr1, drr1, rr2, drr2, ss1, dss1, ss2, dss2 
real*8 :: t, ri, dri, si, dsi, rplot_tmp, zplot_tmp, s_value
real*8 :: aminor, Rgeo, Bgeo, current,beta_p,beta_t,beta_n, dp_int, zjz_int, sum_dl, q, dl, dp_dpsi
real*8 :: PSgi,dPSgi_dr,dPSgi_ds,dPSgi_drs,dPSgi_drr,dPSgi_dss
real*8 :: R0gi,dR0gi_dr,dR0gi_ds,dR0gi_drs,dR0gi_drr,dR0gi_dss
real*8 :: T0gi,dT0gi_dr,dT0gi_ds,dT0gi_drs,dT0gi_drr,dT0gi_dss
real*8 :: Ti0gi,dTi0gi_dr,dTi0gi_ds,dTi0gi_drs,dTi0gi_drr,dTi0gi_dss
real*8 :: Te0gi,dTe0gi_dr,dTe0gi_ds,dTe0gi_drs,dTe0gi_drr,dTe0gi_dss
real*8 :: ZJgi,dZJgi_dr,dZJgi_ds,dZJgi_drs,dZJgi_drr,dZJgi_dss
real*8 :: RRgi,dRRgi_dr,dRRgi_ds,dRRgi_drs,dRRgi_drr,dRRgi_dss
real*8 :: ZZgi,dZZgi_dr,dZZgi_ds,dZZgi_drs,dZZgi_drr,dZZgi_dss
real*8 :: dRRgi_dt, dZZgi_dt, RZjac, PSI_R, PSI_Z, grad_psi, B_tot2, P0gi, dP0gi_dr,dP0gi_ds, P0_R, P0_Z 
real*8 :: density, density_in, density_out, pressure, pressure_in, pressure_out, heat_src_in, heat_src_out, part_src_in, part_src_out
integer :: i_elm_axis, i_elm_xpoint(2), nplot, i, j, n_bnd, i_elm, k, ip, ig
integer :: node1, node2, node3, node4, ifail, my_id

!--------------------------------------- gaussian points between (-1.,1.)
real*8 :: xgs(4), wgs(4)
data xgs /-0.861136311594053, -0.339981043584856, 0.339981043584856,  0.861136311594053 /
data wgs / 0.347854845137454,  0.652145154862546, 0.652145154862546,  0.347854845137454 /

write(*,*) '***************************************'
write(*,*) '* export_helena                       *'
write(*,*) '***************************************'
my_id = 0

call find_axis(my_id,node_list,element_list,psi_axis,R_axis,Z_axis,i_elm_axis,s_axis,t_axis,ifail)

psi_bnd = 0.d0
Z_xpoint(1) = -99.d0
Z_xpoint(2) = +99.d0
   
if (xpoint) then
  call find_xpoint(my_id,node_list,element_list,psi_xpoint,R_xpoint,Z_xpoint,i_elm_xpoint,s_xpoint,t_xpoint,xcase,ifail)
  if (ifail .ne. 1) then      
    psi_bnd  = psi_xpoint(1)
    if( ES%active_xpoint .eq. UPPER_XPOINT ) then
      psi_bnd = psi_xpoint(2)
    endif
    if(xcase .eq. LOWER_XPOINT) Z_xpoint(2) = +99.d0
    if(xcase .eq. UPPER_XPOINT) Z_xpoint(1) = -99.d0
  else
    Z_xpoint(1) = -99.d0
    Z_xpoint(2) = +99.d0
  endif
endif

call find_limiter(my_id,node_list,element_list,bnd_elm_list,psi_lim,R_lim,Z_lim)
if ( (Z_lim .gt. Z_xpoint(1)) .and. (Z_lim .lt. Z_xpoint(2)) ) then
  psi_bnd = min(psi_lim,psi_bnd)
  write(*,'(A,3f8.3)') ' LIMITER PLASMA ',psi_lim,R_lim,Z_lim
endif
    
surface_list%n_psi =3
if (allocated(surface_list%psi_values)) call tr_deallocate(surface_list%psi_values,"surface_list%psi_values",CAT_GRID)
call tr_allocate(surface_list%psi_values,1,3,"surface_list%psi_values",CAT_GRID)
surface_list%psi_values = 0 ! XL : uninitialised value.

if (xpoint) then
  write(*,*) ' x-point plasma'
  if( ES%active_xpoint .eq. UPPER_XPOINT ) then
    surface_list%psi_values(1) =  psi_axis + 0.95  * (psi_xpoint(2) - psi_axis)
    surface_list%psi_values(2) =  psi_axis + 0.99  * (psi_xpoint(2) - psi_axis)
    surface_list%psi_values(3) =  psi_axis + 0.995 * (psi_xpoint(2) - psi_axis)
  else
    surface_list%psi_values(1) =  psi_axis + 0.95  * (psi_xpoint(1) - psi_axis)
    surface_list%psi_values(2) =  psi_axis + 0.99  * (psi_xpoint(1) - psi_axis)
    surface_list%psi_values(3) =  psi_axis + 0.995 * (psi_xpoint(1) - psi_axis)
  endif
  if(xcase .eq. LOWER_XPOINT) Z_xpoint(2) = +99.d0
  if(xcase .eq. UPPER_XPOINT) Z_xpoint(1) = -99.d0
else
  write(*,*) ' NOT an x-point plasma'
  surface_list%psi_values(1) =  psi_axis + 0.95  * (psi_bnd - psi_axis)
  surface_list%psi_values(2) =  psi_axis + 0.99  * (psi_bnd - psi_axis)
  surface_list%psi_values(3) =  psi_axis + 0.999 * (psi_bnd - psi_axis)
endif

call find_flux_surfaces(my_id,xpoint,xcase,node_list,element_list,surface_list)

nplot = 3

open(11,file='equilibrium.txt')

j=3

psi_bnd = surface_list%psi_values(j)

call tr_allocate(rplot,1,surface_list%flux_surfaces(j)%n_pieces * nplot,"rplot",CAT_GRID)
call tr_allocate(zplot,1,surface_list%flux_surfaces(j)%n_pieces * nplot,"zplot",CAT_GRID)

Rmin =  1.d20
Rmax = -1.d20

n_bnd = 0

do k=1,surface_list%flux_surfaces(j)%n_pieces

  i_elm = surface_list%flux_surfaces(j)%elm(k)

  node1 = element_list%element(i_elm)%vertex(1)
  node2 = element_list%element(i_elm)%vertex(2)
  node3 = element_list%element(i_elm)%vertex(3)
  node4 = element_list%element(i_elm)%vertex(4)

  rr1  = surface_list%flux_surfaces(j)%s(1,k)
  drr1 = surface_list%flux_surfaces(j)%s(2,k)
  rr2  = surface_list%flux_surfaces(j)%s(3,k)
  drr2 = surface_list%flux_surfaces(j)%s(4,k)

  ss1  = surface_list%flux_surfaces(j)%t(1,k)
  dss1 = surface_list%flux_surfaces(j)%t(2,k)
  ss2  = surface_list%flux_surfaces(j)%t(3,k)
  dss2 = surface_list%flux_surfaces(j)%t(4,k)

  do ip = 1, nplot

    t = -0.75 + 1.5*float(ip-1)/float(nplot-1)

    call CUB1D(rr1, drr1, rr2, drr2, t, ri, dri)
    call CUB1D(ss1, dss1, ss2, dss2, t, si, dsi)

    call interp_RZ(node_list,element_list,i_elm,ri,si,rplot_tmp,zplot_tmp)

    if ( (zplot_tmp .ge. Z_xpoint(1)) .and. ((zplot_tmp .le. Z_xpoint(2)) ) ) then
      Rmin = min(rplot_tmp,Rmin)
      Rmax = max(rplot_tmp,Rmax)

      n_bnd = n_bnd + 1

      rplot(n_bnd) = rplot_tmp
      zplot(n_bnd) = zplot_tmp

    endif

  enddo
enddo

write(11,'(8e16.8)') R_axis,Z_axis,F0
write(11,'(8e16.8)') surface_list%psi_values(j),psi_axis!,psi_xpoint(1)
write(11,*)          n_bnd

do i=1,n_bnd
  write(11,*) rplot(i),zplot(i)
enddo
call tr_deallocate(rplot,"rplot",CAT_GRID)
call tr_deallocate(zplot,"zplot",CAT_GRID)

aminor = (Rmax - Rmin) /2.d0
Rgeo   = (Rmax + Rmin) /2.d0
Bgeo   = F0 / Rgeo

write(*,'(A,f8.5,A)') ' amin : ',aminor,' m'
write(*,'(A,f8.5,A)') ' Rgeo : ',Rgeo,' m'
write(*,'(A,f8.5,A)') ' Bgeo : ',Bgeo,' T'

call Integrals(node_list, element_list, R_axis, Z_axis, psi_axis, R_xpoint, Z_xpoint, psi_xpoint, psi_bnd, aminor, &
  Bgeo, current, beta_p, beta_t, beta_n, density, density_in, density_out, pressure, pressure_in,  &
  pressure_out, heat_src_in, heat_src_out, part_src_in, part_src_out)
  
write(11,*) aminor, Rgeo, Bgeo
write(11,*) current,beta_p,beta_t,beta_n

n_flux = 201

surface_list%n_psi =n_flux
if (allocated(surface_list%psi_values)) call tr_deallocate(surface_list%psi_values,"surface_list%psi_values",CAT_GRID)
call tr_allocate(surface_list%psi_values,1,surface_list%n_psi,"surface_list%psi_values",CAT_GRID)

do i=1,n_flux
  s_value = float(i)/float(n_flux)
  surface_list%psi_values(i) =  psi_axis + s_value**2 * (psi_bnd - psi_axis)
enddo

call find_flux_surfaces(my_id,xpoint,xcase,node_list,element_list,surface_list)

write(11,*) n_flux-1

do i=2, surface_list%n_psi

  dp_int   = 0.d0
  zjz_int  = 0.d0
  sum_dl   = 0.d0
  q        = 0.d0

  do k=1, surface_list%flux_surfaces(i)%n_pieces

    do ig = 1, 4

      t = xgs(ig)

      rr1  = surface_list%flux_surfaces(i)%s(1,k)
      drr1 = surface_list%flux_surfaces(i)%s(2,k)
      rr2  = surface_list%flux_surfaces(i)%s(3,k)
      drr2 = surface_list%flux_surfaces(i)%s(4,k)

      ss1  = surface_list%flux_surfaces(i)%t(1,k)
      dss1 = surface_list%flux_surfaces(i)%t(2,k)
      ss2  = surface_list%flux_surfaces(i)%t(3,k)
      dss2 = surface_list%flux_surfaces(i)%t(4,k)

      call CUB1D(rr1, drr1, rr2, drr2, t, ri, dri)
      call CUB1D(ss1, dss1, ss2, dss2, t, si, dsi)

      i_elm = surface_list%flux_surfaces(i)%elm(k)

      call interp(node_list,element_list,i_elm,1,1,ri,si,PSgi,dPSgi_dr,dPSgi_ds,dPSgi_drs,dPSgi_drr,dPSgi_dss)
      call interp(node_list,element_list,i_elm,5,1,ri,si,R0gi,dR0gi_dr,dR0gi_ds,dR0gi_drs,dR0gi_drr,dR0gi_dss)
      if (with_TiTe) then
        call interp(node_list,element_list,i_elm,var_Ti,1,ri,si,Ti0gi,dTi0gi_dr,dTi0gi_ds,dTi0gi_drs,dTi0gi_drr,dTi0gi_dss)
        call interp(node_list,element_list,i_elm,var_Te,1,ri,si,Te0gi,dTe0gi_dr,dTe0gi_ds,dTe0gi_drs,dTe0gi_drr,dTe0gi_dss)
        T0gi     = Ti0gi + Te0gi
        dT0gi_dr = dTi0gi_dr + dTe0gi_dr
        dT0gi_ds = dTi0gi_ds + dTe0gi_ds
      else
        call interp(node_list,element_list,i_elm,var_T,1,ri,si,T0gi,dT0gi_dr,dT0gi_ds,dT0gi_drs,dT0gi_drr,dT0gi_dss)
      endif      
      call interp(node_list,element_list,i_elm,3,1,ri,si,ZJgi,dZJgi_dr,dZJgi_ds,dZJgi_drs,dZJgi_drr,dZJgi_dss)

      call interp_RZ(node_list,element_list,i_elm,ri,si,RRgi,dRRgi_dr,dRRgi_ds,dRRgi_drs,dRRgi_drr,dRRgi_dss, &
                                                        ZZgi,dZZgi_dr,dZZgi_ds,dZZgi_drs,dZZgi_drr,dZZgi_dss)
      
      ! --- Make sure that for flux surfaces at Psi_N < 1, the surface integral is carried out only
      !     over the flux surface segments of the plasma region.
      !     I.e., ignore flux surface segments in the private flux region below the x-point.
      if (xcase .ne. UPPER_XPOINT) then 
        if ( xpoint .and. ((PSgi-psi_axis)/(psi_xpoint(1)-psi_axis) < 1.d0) .and. (ZZgi < z_xpoint(1))) cycle
      endif
      if (xcase .ne. LOWER_XPOINT) then
        if ( xpoint .and. ((PSgi-psi_axis)/(psi_xpoint(2)-psi_axis) < 1.d0) .and. (ZZgi > z_xpoint(2))) cycle
      endif

      dRRgi_dt = dRRgi_dr * dri + dRRgi_ds * dsi
      dZZgi_dt = dZZgi_dr * dri + dZZgi_ds * dsi

      dl = sqrt(dRRgi_dt**2 + dZZgi_dt**2)

      RZjac  = DRRgi_dr * dZZgi_ds - dRRgi_ds * dZZgi_dr

      PSI_R = (   dPSgi_dr * dZZgi_ds - dPSgi_ds * dZZgi_dr ) / RZjac
      PSI_Z = ( - dPSgi_dr * dRRgi_ds + dPSgi_ds * dRRgi_dr ) / RZjac

      grad_psi = sqrt(PSI_R * PSI_R + PSI_Z * PSI_Z)

      B_tot2 =  (F0 / RRgi)**2 + (grad_psi/ RRgi)**2

      P0gi     = R0gi * T0gi
      dP0gi_dr = dR0gi_dr * T0gi + R0gi * dT0gi_dr
      dP0gi_ds = dR0gi_ds * T0gi + R0gi * dT0gi_ds

      P0_R = (   dP0gi_dr * dZZgi_ds - dP0gi_ds * dZZgi_dr ) / RZjac
      P0_Z = ( - dP0gi_dr * dRRgi_ds + dP0gi_ds * dRRgi_dr ) / RZjac

      dp_dpsi = (P0_R * PSI_R + P0_Z * PSI_Z) / (PSI_R**2 + PSI_Z**2)

      sum_dl   = sum_dl  +  wgs(ig) * dl
      dp_int   = dp_int  +  wgs(ig) * dP_dpsi * dl
      zjz_int  = zjz_int +  wgs(ig) * ZJgi /RRgi * dl       ! toroidal current density (not JOREKs J_3)
      q        = q       +  wgs(ig) / (RRgi * grad_psi) * dl

    enddo

  enddo

  if ( sum_dl == 0.d0 ) then
    sum_dl = 1.d99
    write(*,*) 'WARNING: Something went wrong in export_helena. sum_dl==0.'
  end if
  
  write(11,'(8e16.8)') surface_list%psi_values(i),dp_int/sum_dl,zjz_int/sum_dl,F0 * q / (2.d0 * PI)

enddo

close(11)

write(*,*) ' end export_helena'

return
end

