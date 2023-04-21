subroutine initial_conditions(my_id,node_list,element_list,bnd_node_list, bnd_elm_list, xpoint2, xcase2)

use constants
use data_structure
use phys_module
use mod_poiss
use equil_info
use mod_sources

implicit none

type (type_node_list)    :: node_list
type (type_element_list) :: element_list
type (type_surface_list) :: surface_list
type (type_bnd_node_list)    :: bnd_node_list
type (type_bnd_element_list) :: bnd_elm_list

integer    :: my_id, i, in, mm, i_elm, ifail, xcase2
real*8     :: amplitude, psi, psi_n, theta
real*8     :: zn, dn_dpsi, dn_dpsi2, dn_dz, dn_dz2, dn_dpsi_dz, dn_dpsi3, dn_dpsi2_dz, dn_dpsi_dz2
real*8     :: zT, dT_dpsi, dT_dpsi2, dT_dz, dT_dz2, dT_dpsi_dz, dT_dpsi3, dT_dpsi2_dz, dT_dpsi_dz2
real*8     :: R, Z, BigR, BigR_s
real*8     :: p_s, p_t, R_s, R_t, Z_s, Z_t, xjac, direction
integer    :: nj
real*8     :: rr,ww, drr_dR, drr_dZ, drr_dR2, drr_dZ2, drr_dRdZ
real*8     :: zjz, dj_dR, dj_dZ, dj_dR_dZ, dj_dR_DR, dj_dZ_dZ
real*8     :: save_zjz, save_dj_dR, save_dj_dZ, save_dj_dR_dZ, save_dj_dR_DR, save_dj_dZ_dZ
logical    :: xpoint2

if (my_id .eq. 0) then
  write(*,*) '***************************************'
  write(*,*) '*      initial conditions  (002)      *'
  write(*,*) '***************************************'
endif

! --- This old projection method should be replaced by a direct solve on the elements
! --- It is much cleaner and more generic
if ( (my_id .eq. 0) .and. (n_order .le. 3) ) then

  do i=1,node_list%n_nodes

    R   = node_list%node(i)%x(1,1,1)
    Z   = node_list%node(i)%x(1,1,2)

    ! --- Density background: profile is made using R instead of psi
    call density    (xpoint2, xcase2, Z, ES%Z_xpoint, R,R_begin,R_end,zn,dn_dpsi,dn_dz,dn_dpsi2,dn_dz2,             &
                                                                      dn_dpsi_dz,dn_dpsi3,dn_dpsi_dz2, dn_dpsi2_dz)
    node_list%node(i)%values(1,1,var_rho) = zn
    node_list%node(i)%values(1,2,var_rho) = dn_dpsi    * node_list%node(i)%values(1,2,1) + dn_dz * node_list%node(i)%x(1,2,2)
    node_list%node(i)%values(1,3,var_rho) = dn_dpsi    * node_list%node(i)%values(1,3,1) + dn_dz * node_list%node(i)%x(1,3,2)
    node_list%node(i)%values(1,4,var_rho) = dn_dpsi    * node_list%node(i)%values(1,4,1) + dn_dz * node_list%node(i)%x(1,4,2) &
                                    + dn_dpsi2   * node_list%node(i)%values(1,2,1) * node_list%node(i)%values(1,3,1)  &
                                    + dn_dz2     * node_list%node(i)%x(1,2,2)        * node_list%node(i)%x(1,3,2)         &
                                    + dn_dpsi_dz * node_list%node(i)%values(1,3,1) * node_list%node(i)%x(1,2,2)         &
                                    + dn_dpsi_dz * node_list%node(i)%values(1,2,1) * node_list%node(i)%x(1,3,2)      

    ! --- Temperature background: profile is made using R instead of psi
    call temperature(xpoint2, xcase2, Z, ES%Z_xpoint, R,R_begin,R_end,zT,dT_dpsi,dT_dz,dT_dpsi2,dT_dz2,             &
                                                                      dT_dpsi_dz,dT_dpsi3,dT_dpsi_dz2, dT_dpsi2_dz)
    node_list%node(i)%values(1,1,var_T) = zT
    node_list%node(i)%values(1,2,var_T) = dT_dpsi    * node_list%node(i)%values(1,2,1) + dT_dz * node_list%node(i)%x(1,2,2)
    node_list%node(i)%values(1,3,var_T) = dT_dpsi    * node_list%node(i)%values(1,3,1) + dT_dz * node_list%node(i)%x(1,3,2)
    node_list%node(i)%values(1,4,var_T) = dT_dpsi    * node_list%node(i)%values(1,4,1) + dT_dz * node_list%node(i)%x(1,4,2) &
                                    + dT_dpsi2   * node_list%node(i)%values(1,2,1) * node_list%node(i)%values(1,3,1)  &
                                    + dT_dz2     * node_list%node(i)%x(1,2,2)        * node_list%node(i)%x(1,3,2)         &
                                    + dT_dpsi_dz * node_list%node(i)%values(1,3,1) * node_list%node(i)%x(1,2,2)         &
                                    + dT_dpsi_dz * node_list%node(i)%values(1,2,1) * node_list%node(i)%x(1,3,2)      

    ! --- Use current ropes to define density blobs
    do nj=1,n_jropes
      rr = sqrt((R-R_jropes(nj))**2 + (Z-Z_jropes(nj))**2)
      drr_dR   = (R-R_jropes(nj)) / rr
      drr_dZ   = (Z-Z_jropes(nj)) / rr
      drr_dR2  = 1./rr - (R-R_jropes(nj)) / rr**2 * drr_dR
      drr_dZ2  = 1./rr - (Z-Z_jropes(nj)) / rr**2 * drr_dZ
      drr_dRdZ = - (R-R_jropes(nj)) / rr**2 * drr_dZ
      ww = w_jropes(nj)
      save_zjz        = 0.d0
      save_dj_dR      = 0.d0
      save_dj_dZ      = 0.d0
      save_dj_dR_dR   = 0.d0
      save_dj_dZ_dZ   = 0.d0
      save_dj_dR_dZ   = 0.d0
      if (rr .lt. 1.00*ww) then
        save_zjz        =                     (1.0 - (rr/ww)**2 )**2
        save_dj_dR      = -4. * rr   /ww**2 * (1.0 - (rr/ww)**2 ) * drr_dR
        save_dj_dZ      = -4. * rr   /ww**2 * (1.0 - (rr/ww)**2 ) * drr_dZ
        save_dj_dR_dR   = -4.        /ww**2 * (1.0 - (rr/ww)**2 ) * drr_dR**2 & 
                          -4. * rr   /ww**2 * (1.0 - (rr/ww)**2 ) * drr_dR2   & 
                          +8. * rr**2/ww**4                       * drr_dR**2 
        save_dj_dZ_dZ   = -4.        /ww**2 * (1.0 - (rr/ww)**2 ) * drr_dZ**2 & 
                          -4. * rr   /ww**2 * (1.0 - (rr/ww)**2 ) * drr_dZ2   & 
                          +8. * rr**2/ww**4                       * drr_dZ**2 
        save_dj_dR_dZ   = -4.        /ww**2 * (1.0 - (rr/ww)**2 ) * drr_dR*drr_dZ & 
                          -4. * rr   /ww**2 * (1.0 - (rr/ww)**2 ) * drr_dRdZ      & 
                          +8. * rr**2/ww**4                       * drr_dR*drr_dZ 
      endif
     
      ! --- Density
      zjz        = save_zjz      * rho_jropes(nj)
      dj_dR      = save_dj_dR    * rho_jropes(nj)
      dj_dZ      = save_dj_dZ    * rho_jropes(nj)
      dj_dR_dR   = save_dj_dR_dR * rho_jropes(nj)
      dj_dZ_dZ   = save_dj_dZ_dZ * rho_jropes(nj)
      dj_dR_dZ   = save_dj_dR_dZ * rho_jropes(nj)

      node_list%node(i)%values(1,1,var_rho) = node_list%node(i)%values(1,1,var_rho) + zjz
     
      node_list%node(i)%values(1,2,var_rho) = node_list%node(i)%values(1,2,var_rho) &
                                            + dj_dR   * node_list%node(i)%x(1,2,1)  &
                                            + dj_dZ   * node_list%node(i)%x(1,2,2)
     
      node_list%node(i)%values(1,3,var_rho) = node_list%node(i)%values(1,3,var_rho) &
                                            + dj_dR   * node_list%node(i)%x(1,3,1)  &
                                            + dj_dZ   * node_list%node(i)%x(1,3,2)
     
      node_list%node(i)%values(1,4,var_rho) = node_list%node(i)%values(1,4,var_rho) &
                                            + dj_dR    * node_list%node(i)%x(1,4,1) &
                                            + dj_dZ    * node_list%node(i)%x(1,4,2) &
                                            + dj_dR_dR * node_list%node(i)%x(1,2,1) * node_list%node(i)%x(1,3,1)    &
                                            + dj_dZ_dZ * node_list%node(i)%x(1,2,2) * node_list%node(i)%x(1,3,2)    &
                                            + dj_dR_dZ * ( node_list%node(i)%x(1,2,1) * node_list%node(i)%x(1,3,2)  &
                                                         + node_list%node(i)%x(1,3,1) * node_list%node(i)%x(1,2,2) )

     
      ! --- Temperature
      zjz        = save_zjz      * T_jropes(nj)
      dj_dR      = save_dj_dR    * T_jropes(nj)
      dj_dZ      = save_dj_dZ    * T_jropes(nj)
      dj_dR_dR   = save_dj_dR_dR * T_jropes(nj)
      dj_dZ_dZ   = save_dj_dZ_dZ * T_jropes(nj)
      dj_dR_dZ   = save_dj_dR_dZ * T_jropes(nj)

      node_list%node(i)%values(1,1,var_T) = node_list%node(i)%values(1,1,var_T) + zjz
     
      node_list%node(i)%values(1,2,var_T) = node_list%node(i)%values(1,2,var_T) &
                                            + dj_dR   * node_list%node(i)%x(1,2,1)  &
                                            + dj_dZ   * node_list%node(i)%x(1,2,2)
     
      node_list%node(i)%values(1,3,var_T) = node_list%node(i)%values(1,3,var_T) &
                                            + dj_dR   * node_list%node(i)%x(1,3,1)  &
                                            + dj_dZ   * node_list%node(i)%x(1,3,2)
     
      node_list%node(i)%values(1,4,var_T) = node_list%node(i)%values(1,4,var_T) &
                                            + dj_dR    * node_list%node(i)%x(1,4,1) &
                                            + dj_dZ    * node_list%node(i)%x(1,4,2) &
                                            + dj_dR_dR * node_list%node(i)%x(1,2,1) * node_list%node(i)%x(1,3,1)    &
                                            + dj_dZ_dZ * node_list%node(i)%x(1,2,2) * node_list%node(i)%x(1,3,2)    &
                                            + dj_dR_dZ * ( node_list%node(i)%x(1,2,1) * node_list%node(i)%x(1,3,2)  &
                                                         + node_list%node(i)%x(1,3,1) * node_list%node(i)%x(1,2,2) )

    enddo

    node_list%node(i)%values(1,:,var_u) = 0.d0 
    node_list%node(i)%values(1,:,var_w) = 0.d0 

    node_list%node(i)%deltas = 0.d0

  enddo

endif ! n_order<=3

! --- Variable projection is better at higher order...
! --- (by the way, we could use this for n_order=3 and remove all the above as well, 
! --- and remove all derivatives from profiles functions, which are not really needed, 
! --- except dn_dpsi and dT_dpsi for current profile...)
if (n_order .ge. 5) then
  call Poisson(my_id,0,node_list,element_list,bnd_node_list,bnd_elm_list, &
               var_psi,var_rho,1, ES%psi_axis,ES%psi_bnd,xpoint2,xcase2,ES%Z_xpoint,freeboundary_equil,refinement,1)
  call Poisson(my_id,0,node_list,element_list,bnd_node_list,bnd_elm_list, &
               var_psi,var_T,1, ES%psi_axis,ES%psi_bnd,xpoint2,xcase2,ES%Z_xpoint,freeboundary_equil,refinement,1)
endif


    
!---------------------------- initialise perturbations
amplitude = 1.d-12
mm = 2

do in=2,n_tor

  if (my_id .eq. 0) then

    do i=1,node_list%n_nodes

      node_list%node(i)%values(in,:,:) = 0.d0

      node_list%node(i)%values(in,1,var_w) = amplitude
      node_list%node(i)%values(in,2,var_w) = 0.d0
      node_list%node(i)%values(in,3,var_w) = 0.d0
      node_list%node(i)%values(in,4,var_w) = 0.d0
      
      node_list%node(i)%deltas = 0.d0

    enddo

  endif

  call Poisson(my_id,1,node_list,element_list,bnd_node_list,bnd_elm_list, &
               var_w,var_u,1, ES%psi_axis,ES%psi_bnd,xpoint2, xcase2,ES%Z_xpoint,freeboundary_equil,refinement,1)
enddo

return
end
