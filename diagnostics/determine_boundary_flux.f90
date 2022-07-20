!> Calculates n.B on the plasma boundary

subroutine determine_boundary_flux(node_list,element_list)

use constants
use tr_module 
use data_structure
use phys_module
use equil_info, only : get_psi_n
use mod_interp
use mod_chi
use mod_parameters


implicit none

! --- Gaussian points between (-1.,1.) for Gauss-integration
real*8, parameter :: xgs(4) = (/-0.861136311594053, -0.339981043584856, 0.339981043584856,  0.861136311594053 /)
real*8, parameter :: wgs(4) = (/ 0.347854845137454,  0.652145154862546, 0.652145154862546,  0.347854845137454 /)

! --- Input parameters.
type (type_node_list),    intent(in)    :: node_list
type (type_element_list), intent(in)    :: element_list

! --- Local variables
integer :: n_int
integer :: i_elm, j, k, n1, n2, n3
real*8  :: p, t,rr1, rr2, drr1, drr2, ss1, ss2, dss1, dss2, ri, si, dri, dsi, dA
real*8  :: RRgi, dRRgi_dr, dRRgi_ds, dRRgi_dp, ZZgi, dZZgi_dr, dZZgi_ds, dZZgi_dp, dRRgi_dt, dZZgi_dt
real*8  :: PSgi, dPSgii_dr, dPSgii_ds, dPSgi_dr, dPSgi_ds, dPSgi_dp, PSI_R, PSI_Z, Psi_p, RZJAC, grad_psi, psi_n
real*8  :: Fgi,dFgi_dr,dFgi_ds,dFgi_drs,dFgi_drr,dFgi_dss
real*8  :: surface_area, sum_dA, sum_dA_abs, B_tot2, delta_phi
real*8  :: dRRgi_drs,dRRgi_drr,dRRgi_dss, dRRgi_drp,dRRgi_dsp,dRRgi_dpp
real*8  :: dZZgi_drs,dZZgi_drr,dZZgi_dss, dZZgi_drp,dZZgi_dsp,dZZgi_dpp, dPSgi_drs,dPSgi_drr,dPSgi_dss
integer :: m, ig1, ig2, i_plane, i_tor
real*8  :: s_phi, c_phi, ndotB, cross_deriv(3), n_perp(3), B_boundary(3)
real*8  :: chi(0:n_order-1,0:n_order-1,0:n_order-1)
real*8  :: ndotB_max=0.0


write(*,*) "*********************************"
write(*,*) "*    Determine Boundary Flux    *"
write(*,*) "*********************************"

delta_phi = 2 * PI / float(n_plane) / float(n_period)
sum_dA = 0.d0
sum_dA_abs = 0.d0
surface_area = 0.d0
do i_elm=(n_flux-2)*n_tht+1, (n_flux-1)*n_tht
  do i_plane=1,n_plane
    do ig1 = 1, 4
      si = 0.5 * (xgs(ig1) + 1.0)
      ri = 1.0

      call interp_RZP(node_list,element_list,i_elm,ri,si,(i_plane-1)*delta_phi,   &
                      RRgi,dRRgi_dr,dRRgi_ds,dRRgi_dp,dRRgi_drs,dRRgi_drr,dRRgi_dss,dRRgi_drp, dRRgi_dsp, dRRgi_dpp, &
                      ZZgi,dZZgi_dr,dZZgi_ds,dZZgi_dp,dZZgi_drs,dZZgi_drr,dZZgi_dss,dZZgi_drp, dZZgi_dsp, dZZgi_dpp)
      do ig2 = 1, 4
        p = (i_plane - 1 + 0.5 * (xgs(ig2) + 1.0)) * delta_phi

        call interp(node_list,element_list,i_elm,1,1,ri,si,PSgi,dPSgi_dr,dPSgi_ds,dPSgi_drs,dPSgi_drr,dPSgi_dss)
        dPSgi_dp = 0.d0
        do i_tor=1,(n_tor-1)/2
          call interp(node_list,element_list,i_elm,1,2*i_tor,ri,si,PSgi,dPSgii_dr,dPSgii_ds,dPSgi_drs,dPSgi_drr,dPSgi_dss)
          dPSgi_dr = dPSgi_dr + dPSgii_dr*cos(mode(2*i_tor)*p)
          dPSgi_ds = dPSgi_ds + dPSgii_ds*cos(mode(2*i_tor)*p)
          dPSgi_dp = dPSgi_dp - mode(2*i_tor)*PSgi*sin(mode(2*i_tor)*p)
          call interp(node_list,element_list,i_elm,1,2*i_tor+1,ri,si,PSgi,dPSgii_dr,dPSgii_ds,dPSgi_drs,dPSgi_drr,dPSgi_dss)
          dPSgi_dr = dPSgi_dr + dPSgii_dr*sin(mode(2*i_tor+1)*p)
          dPSgi_ds = dPSgi_ds + dPSgii_ds*sin(mode(2*i_tor+1)*p)
          dPSgi_dp = dPSgi_dp + mode(2*i_tor+1)*PSgi*cos(mode(2*i_tor+1)*p)
        end do

        call interp_RZP(node_list,element_list,i_elm,ri,si,p,   &
                        RRgi,dRRgi_dr,dRRgi_ds,dRRgi_dp,dRRgi_drs,dRRgi_drr,dRRgi_dss,dRRgi_drp, dRRgi_dsp, dRRgi_dpp, &
                        ZZgi,dZZgi_dr,dZZgi_ds,dZZgi_dp,dZZgi_drs,dZZgi_drr,dZZgi_dss,dZZgi_drp, dZZgi_dsp, dZZgi_dpp)

        chi = get_chi(RRgi, ZZgi, p)
        

        ! Radial component discarded because grid is already flux surface aligned
        dRRgi_dt = dRRgi_ds  ! + dRRgi_dr * dri
        dZZgi_dt = dZZgi_ds  ! + dZZgi_dr * dri
        
        ! Calculate normal to boundary
        n_perp = (/-dZZgi_dt, dRRgi_dt, (dRRgi_dp*dZZgi_dt-dRRgi_dt*dZZgi_dp)/RRgi /)
        n_perp = n_perp / sqrt(sum(n_perp*n_perp))

        ! Calculate surface area contribution from covariant components
        c_phi = cos(p)
        s_phi = sin(p)
        cross_deriv = (/(dRRgi_dt*s_phi*dZZgi_dp-(dRRgi_dp*s_phi+RRgi*c_phi)*dZZgi_dt),   &
                      (-dRRgi_dt*dZZgi_dp*c_phi+(dRRgi_dp*c_phi-RRgi*s_phi)*dZZgi_dt),  &
                      (dRRgi_dt * RRgi)  /)
        dA = sqrt(cross_deriv(1)*cross_deriv(1) + cross_deriv(2)*cross_deriv(2) + cross_deriv(3)*cross_deriv(3)) 

        ! Calculate n.B
        RZjac  = DRRgi_dr * dZZgi_ds - dRRgi_ds * dZZgi_dr
        PSI_R = (   dPSgi_dr * dZZgi_ds - dPSgi_ds * dZZgi_dr ) / RZjac
        PSI_Z = ( - dPSgi_dr * dRRgi_ds + dPSgi_ds * dRRgi_dr ) / RZjac
        Psi_p = dPSgi_dp - Psi_R*dRRgi_dp - Psi_z*dZZgi_dp
        B_boundary = (/ chi(1,0,0) + (Psi_z*chi(0,0,1) - Psi_p*chi(0,1,0))/(F0*RRgi), &
                        chi(0,1,0) - (Psi_R*chi(0,0,1) - Psi_p*chi(1,0,0))/(F0*RRgi), &
                        chi(0,0,1)/RRgi + (Psi_R*chi(0,1,0) - Psi_z*chi(1,0,0))/F0 /)
        ndotB = sum(n_perp*B_boundary)      
        ndotB_max = max(abs(ndotB), ndotB_max)

        ! Factors of 0.5 to convert the integration interval from [-1,1] to [0,1] (poloidal direction) 
        !   and 0.5*delta_phi to convert the integration interval from [-1,1] to [phi_i,phi_i+delta_phi] (toroidal direction)
        surface_area = surface_area +  wgs(ig1) * 0.5 * wgs(ig2) * 0.5 * delta_phi * dA
        sum_dA = sum_dA +  wgs(ig1) * 0.5 * wgs(ig2) * 0.5 * delta_phi * dA * ndotB 
        sum_dA_abs = sum_dA_abs +  wgs(ig1) * 0.5 * wgs(ig2) * 0.5 * delta_phi * dA * abs(ndotB) 
      end do
    end do
  end do
end do

write(*,*) "Max n.B: ", ndotB_max
write(*,*) "Surface area:        ", n_period * surface_area, "m^2"
write(*,*) "Integrated abs(n.B): ", n_period * sum_dA_abs, "Tm^2"
write(*,*) "Total Boundary Flux: ", n_period * sum_dA, "Tm^2"

end subroutine determine_boundary_flux
