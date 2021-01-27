subroutine initial_conditions(my_id,node_list,element_list,bnd_node_list, bnd_elm_list, xpoint2, xcase2)
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
use constants
use data_structure
use phys_module
use mod_poiss
use equil_info
use mod_interp, only: interp

implicit none

type (type_node_list)    :: node_list
type (type_element_list) :: element_list
type (type_surface_list) :: surface_list
type (type_bnd_node_list)    :: bnd_node_list
type (type_bnd_element_list) :: bnd_elm_list

integer    :: my_id, i, in, mm, i_elm_axis, i_elm_xpoint(2), i_elm, ifail, xcase2
real*8     :: amplitude, psi, psi_n, theta
real*8     :: zn, dn_dpsi, dn_dpsi2, dn_dz, dn_dz2, dn_dpsi_dz, dn_dpsi3, dn_dpsi2_dz, dn_dpsi_dz2
real*8     :: zTi, dTi_dpsi, dTi_dpsi2, dTi_dz, dTi_dz2, dTi_dpsi_dz, dTi_dpsi3, dTi_dpsi2_dz, dTi_dpsi_dz2
real*8     :: zTe, dTe_dpsi, dTe_dpsi2, dTe_dz, dTe_dz2, dTe_dpsi_dz, dTe_dpsi3, dTe_dpsi2_dz, dTe_dpsi_dz2
real*8     :: zFFprime,dFFprime_dpsi,dFFprime_dz, dFFprime_dpsi_dz, dFFprime_dz2, dFFprime_dpsi2
real*8     :: R, Z, BigR, Ti0, Ti0_s, Te0, Te0_s, BigR_s
real*8     :: zjz, dj_dpsi, dj_dR, dj_dZ, dj_dR_dZ, dj_dR_DR, dj_dZ_dZ, dj_dpsi2, dj_dR_dpsi, dj_dZ_dpsi
real*8     :: zp, dp_dpsi, dp_dpsi2, dp_dz, dp_dz2, P_ss, P_st, P_tt, R_out,Z_out,s_out,t_out
real*8     :: ps0_s, ps0_t, p_s, p_t, zj0_s, zj0_t, ps0_x, ps0_y, R_s, R_t, Z_s, Z_t, xjac, direction, Btot
real*8     :: zV, dV_dpsi, dV_dpsi2, dV_dz, dV_dz2, dV_dpsi_dz, dV_dpsi3, dV_dpsi2_dz, dV_dpsi_dz2
real*8     :: Omega, dOmega_dpsi, dOmega_dz, dOmega_dpsi2, dOmega_dz2, dOmega_dpsi_dz
logical    :: xpoint2
real*8     :: theta_bal, theta_pol, gauss_amp, perturbation, R_left, q95
integer    :: mod_n

if (my_id .eq. 0) then
  write(*,*) '***************************************'
  write(*,*) '*      initial conditions (400)       *'
  write(*,*) '***************************************'
endif

if (my_id .eq. 0) then

  do i=1,node_list%n_nodes

    psi = node_list%node(i)%values(1,1,1)
    R   = node_list%node(i)%x(1,1,1)
    Z   = node_list%node(i)%x(1,1,2)

    call density(      xpoint2, xcase2, Z, ES%Z_xpoint, psi,ES%psi_axis,ES%psi_bnd, &
                       zn,dn_dpsi,dn_dz,dn_dpsi2,dn_dz2,dn_dpsi_dz,dn_dpsi3,dn_dpsi_dz2, dn_dpsi2_dz)

    call temperature_i(xpoint2, xcase2, Z, ES%Z_xpoint, psi,ES%psi_axis,ES%psi_bnd, &
    		       zTi,dTi_dpsi,dTi_dz,dTi_dpsi2,dTi_dz2,dTi_dpsi_dz,dTi_dpsi3,dTi_dpsi_dz2,dTi_dpsi2_dz)
    
    call temperature_e(xpoint2, xcase2, Z, ES%Z_xpoint, psi,ES%psi_axis,ES%psi_bnd, &
    		       zTe,dTe_dpsi,dTe_dz,dTe_dpsi2,dTe_dz2,dTe_dpsi_dz,dTe_dpsi3,dTe_dpsi_dz2,dTe_dpsi2_dz)

    call FFprime(      xpoint2, xcase2, Z, ES%Z_xpoint, psi,ES%psi_axis,ES%psi_bnd, &
                       zFFprime,dFFprime_dpsi,dFFprime_dz,dFFprime_dpsi2,dFFprime_dz2, dFFprime_dpsi_dz, .true.)
    
!============================MB
    if ( (abs(V_0) .ge. 1.d-19) .or. (num_rot) ) then
       if (normalized_velocity_profile) then
          call velocity(xpoint2, xcase2, Z, ES%Z_xpoint, psi,ES%psi_axis,ES%psi_bnd,zV,dV_dpsi,dV_dz,dV_dpsi2,dV_dz2, &
               dV_dpsi_dz,dV_dpsi3,dV_dpsi_dz2, dV_dpsi2_dz)
       else
          call velocity(xpoint2, xcase2, Z, ES%Z_xpoint, psi,ES%psi_axis,ES%psi_bnd,Omega,dOmega_dpsi,dOmega_dz,dOmega_dpsi2,dOmega_dz2, &
               dOmega_dpsi_dz,dV_dpsi3,dV_dpsi_dz2, dV_dpsi2_dz)
       endif
    endif
!============================MB

    zp       = zn * (zTi + zTe)
    dp_dpsi  = zn * (dTi_dpsi + dTe_dpsi) + dn_dpsi * (zTi + zTe)
    dp_dpsi2 = zn * (dTi_dpsi2 + dTe_dpsi2) + 2.d0 * dn_dpsi * (dTi_dpsi + dTe_dpsi) + dn_dpsi2 * (zTi + zTe)
    dp_dz    = zn * (dTi_dz + dTe_dz) + dn_dz * (zTi + zTe)
    dp_dz2   = zn * (dTi_dz2 + dTe_dz2) + 2.d0 * dn_dz * (dTi_dz + dTe_dz) + dn_dz2 * (zTi + zTe)				       

    node_list%node(i)%values(1,1,5) = zn
    node_list%node(i)%values(1,2,5) = dn_dpsi  * node_list%node(i)%values(1,2,1) + dn_dz * node_list%node(i)%x(1,2,2)
    node_list%node(i)%values(1,3,5) = dn_dpsi  * node_list%node(i)%values(1,3,1) + dn_dz * node_list%node(i)%x(1,3,2)
    node_list%node(i)%values(1,4,5) = dn_dpsi    * node_list%node(i)%values(1,4,1) + dn_dz * node_list%node(i)%x(1,4,2) &
                                    + dn_dpsi2   * node_list%node(i)%values(1,2,1) * node_list%node(i)%values(1,3,1)  &
                                    + dn_dz2     * node_list%node(i)%x(1,2,2)        * node_list%node(i)%x(1,3,2)         &
                                    + dn_dpsi_dz * node_list%node(i)%values(1,3,1) * node_list%node(i)%x(1,2,2)         &
                                    + dn_dpsi_dz * node_list%node(i)%values(1,2,1) * node_list%node(i)%x(1,3,2)      

    node_list%node(i)%values(1,1,6) = zTi
    node_list%node(i)%values(1,2,6) = dTi_dpsi    * node_list%node(i)%values(1,2,1) + dTi_dz * node_list%node(i)%x(1,2,2)
    node_list%node(i)%values(1,3,6) = dTi_dpsi    * node_list%node(i)%values(1,3,1) + dTi_dz * node_list%node(i)%x(1,3,2)
    node_list%node(i)%values(1,4,6) = dTi_dpsi    * node_list%node(i)%values(1,4,1) + dTi_dz * node_list%node(i)%x(1,4,2) &
                                    + dTi_dpsi2   * node_list%node(i)%values(1,2,1) * node_list%node(i)%values(1,3,1)  &
                                    + dTi_dz2     * node_list%node(i)%x(1,2,2)        * node_list%node(i)%x(1,3,2)         &
                                    + dTi_dpsi_dz * node_list%node(i)%values(1,3,1) * node_list%node(i)%x(1,2,2)         &
                                    + dTi_dpsi_dz * node_list%node(i)%values(1,2,1) * node_list%node(i)%x(1,3,2)      

    node_list%node(i)%values(1,1,8) = zTe
    node_list%node(i)%values(1,2,8) = dTe_dpsi    * node_list%node(i)%values(1,2,1) + dTe_dz * node_list%node(i)%x(1,2,2)
    node_list%node(i)%values(1,3,8) = dTe_dpsi    * node_list%node(i)%values(1,3,1) + dTe_dz * node_list%node(i)%x(1,3,2)
    node_list%node(i)%values(1,4,8) = dTe_dpsi    * node_list%node(i)%values(1,4,1) + dTe_dz * node_list%node(i)%x(1,4,2) &
                                    + dTe_dpsi2   * node_list%node(i)%values(1,2,1) * node_list%node(i)%values(1,3,1)  &
                                    + dTe_dz2     * node_list%node(i)%x(1,2,2)        * node_list%node(i)%x(1,3,2)         &
                                    + dTe_dpsi_dz * node_list%node(i)%values(1,3,1) * node_list%node(i)%x(1,2,2)         &
                                    + dTe_dpsi_dz * node_list%node(i)%values(1,2,1) * node_list%node(i)%x(1,3,2)      

    node_list%node(i)%values(1,1,2) = - tauIC * zp 
    node_list%node(i)%values(1,2,2) = - tauIC * (dp_dpsi  * node_list%node(i)%values(1,2,1) + dp_dz * node_list%node(i)%x(1,2,2))
    node_list%node(i)%values(1,3,2) = - tauIC * (dp_dpsi  * node_list%node(i)%values(1,3,1) + dp_dz * node_list%node(i)%x(1,3,2))
    node_list%node(i)%values(1,4,2) = - tauIC * (dp_dpsi  * node_list%node(i)%values(1,4,1) + dp_dz * node_list%node(i)%x(1,4,2) &
                                    + dp_dpsi2 * node_list%node(i)%values(1,2,1) * node_list%node(i)%values(1,3,1)  &
                                    + dp_dz2   * node_list%node(i)%x(1,2,2)        * node_list%node(i)%x(1,3,2) )

    node_list%node(i)%values(1,1,7) = 0.d0        ! parallel velocity
    node_list%node(i)%values(1,2,7) = 0.d0
    node_list%node(i)%values(1,3,7) = 0.d0
    node_list%node(i)%values(1,4,7) = 0.d0

    if ( (abs(V_0) .ge. 1.d-19) .or. (num_rot) ) then
      if (normalized_velocity_profile) then
        node_list%node(i)%values(1,1,7) = zV
        node_list%node(i)%values(1,2,7) = dV_dpsi  * node_list%node(i)%values(1,2,1) + dV_dz * node_list%node(i)%x(1,2,2)
        node_list%node(i)%values(1,3,7) = dV_dpsi  * node_list%node(i)%values(1,3,1) + dV_dz * node_list%node(i)%x(1,3,2)
        node_list%node(i)%values(1,4,7) = dV_dpsi    * node_list%node(i)%values(1,4,1) + dV_dz * node_list%node(i)%x(1,4,2) &
                                        + dV_dpsi2   * node_list%node(i)%values(1,2,1) * node_list%node(i)%values(1,3,1)  &
                                        + dV_dz2     * node_list%node(i)%x(1,2,2)        * node_list%node(i)%x(1,3,2)         &
                                        + dV_dpsi_dz * node_list%node(i)%values(1,3,1) * node_list%node(i)%x(1,2,2)         &
                                        + dV_dpsi_dz * node_list%node(i)%values(1,2,1) * node_list%node(i)%x(1,3,2) 
      else
        node_list%node(i)%values(1,1,7) = R**2 * Omega
        node_list%node(i)%values(1,2,7) = 2.d0 * R * node_list%node(i)%x(1,2,1) * Omega               &
                                          + R**2 * dOmega_dpsi  * node_list%node(i)%values(1,2,1)   &
                                          + R**2 * dOmega_dz    * node_list%node(i)%x(1,2,2)
        node_list%node(i)%values(1,3,7) = 2.d0 * R * node_list%node(i)%x(1,3,1) * Omega               &
                                          + R**2 * dOmega_dpsi  * node_list%node(i)%values(1,3,1)   &
                                          + R**2 * dOmega_dz * node_list%node(i)%x(1,3,2)
   
        node_list%node(i)%values(1,4,7) =   2.d0 *     node_list%node(i)%x(1,2,1)**2  * Omega &
                                          + 2.d0 * R * node_list%node(i)%x(1,4,1)     * Omega &
                                          + 2.d0 * R * node_list%node(i)%x(1,2,1)     * dOmega_dpsi  * node_list%node(i)%values(1,3,1) &
                                          + 2.d0 * R * node_list%node(i)%x(1,2,1)     * dOmega_dz * node_list%node(i)%x(1,3,2) 
        node_list%node(i)%values(1,4,7) = node_list%node(i)%values(1,4,7) &
                                          + 2.d0 * R * node_list%node(i)%x(1,3,1) * dOmega_dpsi  * node_list%node(i)%values(1,2,1)        &
                                          + R**2 * dOmega_dpsi2   * node_list%node(i)%values(1,3,1) * node_list%node(i)%values(1,2,1)   &
                                          + R**2 * dOmega_dpsi_dz * node_list%node(i)%x(1,3,2)        * node_list%node(i)%values(1,2,1)   &
                                          + R**2 * dOmega_dpsi    * node_list%node(i)%values(1,4,1)
        node_list%node(i)%values(1,4,7) = node_list%node(i)%values(1,4,7) &
                                          + 2.d0 * R * node_list%node(i)%x(1,3,1) * dOmega_dz    * node_list%node(i)%x(1,2,2)      &
                                          + R**2 * dOmega_dpsi_dz * node_list%node(i)%values(1,3,1) * node_list%node(i)%x(1,3,2) &
                                          + R**2 * dOmega_dz2     * node_list%node(i)%x(1,3,2)        * node_list%node(i)%x(1,3,2) &
                                          + R**2 * dOmega_dz * node_list%node(i)%x(1,4,2)
   
        node_list%node(i)%values(1,1,7) = 2.d0 * PI / F0 * node_list%node(i)%values(1,1,7)
        node_list%node(i)%values(1,2,7) = 2.d0 * PI / F0 * node_list%node(i)%values(1,2,7)
        node_list%node(i)%values(1,3,7) = 2.d0 * PI / F0 * node_list%node(i)%values(1,3,7)
        node_list%node(i)%values(1,4,7) = 2.d0 * PI / F0 * node_list%node(i)%values(1,4,7)
      endif
    end if
    
    node_list%node(i)%deltas = 0.d0
    
  enddo

endif

if (tauIC .ne. 0.d0) then
  call Poisson(my_id,2,node_list,element_list,bnd_node_list,bnd_elm_list, &
               2,4,1, ES%psi_axis,ES%psi_bnd,xpoint2, xcase2,ES%Z_xpoint,freeboundary_equil,refinement,1)      ! inverse Poisson
endif


!----------------------------------------- flux boundary perturbation (to be completed, see Marina)
!if (my_id .eq. 0) then
!  do i=1,node_list%n_nodes
!    psi = node_list%node(i)%values(1,1,1)
!    R   = node_list%node(i)%x(1,1,1)
!    Z   = node_list%node(i)%x(1,1,2)
!    theta = atan2(Z-ES%Z_axis, R-ES%R_axis)
!    psi_bnd = 0.d0
!    if (xpoint2 .and. (xcase2 .ne. 2)) psi_bnd = psi_xpoint(1)
!    if (xpoint2 .and. (xcase2 .eq. 2)) psi_bnd = psi_xpoint(2)
!    if (xpoint2 .and. (xcase2 .eq. 3) .and. (psi_xpoint(2) .lt. psi_xpoint(1))) psi_bnd = psi_xpoint(2)
!    psi_n = (psi - psi_axis)/(psi_bnd - psi_axis)     
!    if (node_list%node(i)%boundary .ne. 0) then
!      node_list%node(i)%values(2,1,1) =  0.01 * sin(2.d0*theta)
!      node_list%node(i)%values(3,1,1) = -0.01 * cos(2.d0*theta) 
!      node_list%node(i)%values(2,3,1) =  0.01 * cos(2.d0*theta) * 2.d0*PI/float(n_tht)
!      node_list%node(i)%values(3,3,1) = +0.01 * sin(2.d0*theta) * 2.d0*PI/float(n_tht)
!    endif
!  enddo
!call poisson(my_id,-2,node_list,element_list,1,3,2,xpoint2, xcase2) 
!call poisson(my_id,-2,node_list,element_list,1,3,3,xpoint2, xcase2) 
!endif    

!---------------------------- initialise perturbations
amplitude = 1.d-12
mm = 2

do in=2,n_tor

  if (my_id .eq. 0) then

    do i=1,node_list%n_nodes

      node_list%node(i)%values(in,:,:) = 0.d0

      psi = node_list%node(i)%values(1,1,1)
      R   = node_list%node(i)%x(1,1,1)
      Z   = node_list%node(i)%x(1,1,2)
      psi_n = (psi - ES%psi_axis)/(ES%psi_bnd - ES%psi_axis)
      if (xpoint2) then
        if ((psi_n .lt. 1.d0) .and. (Z .lt. ES%Z_xpoint(1)) .and. (xcase2 .ne. 2)) then
          psi_n = 2.d0 - psi_n
        endif
        if ((psi_n .lt. 1.d0) .and. (Z .gt. ES%Z_xpoint(2)) .and. (xcase2 .ne. 1)) then
          psi_n = 2.d0 - psi_n
        endif
      endif

      !node_list%node(i)%values(in,1,4) = amplitude * psi_n * (1.d0 -psi_n)
      !node_list%node(i)%values(in,2,4) = amplitude * (1. - 2.d0 * psi_n)/(ES%psi_bnd - psi_axis) * node_list%node(i)%values(1,2,1)
      !node_list%node(i)%values(in,3,4) = amplitude * (1. - 2.d0 * psi_n)/(ES%psi_bnd - psi_axis) * node_list%node(i)%values(1,3,1)
      !node_list%node(i)%values(in,4,4) = amplitude * (1. - 2.d0 * psi_n)/(ES%psi_bnd - psi_axis) * node_list%node(i)%values(1,4,1)

      ! --- Trying a ballooning mode perturbation
      q95 = 3.3 ! just estimate for now...
      R_left = 1.95
      if (mod(in,2) .eq. 0) then
        mod_n = in/2 * n_period
      else
        mod_n = (in-1)/2 * n_period
      endif
      gauss_amp = amplitude * exp( -(psi_n-0.90)**2/(2*0.07**2) ) 
      theta_pol = atan2(Z-ES%Z_axis,R-ES%R_axis)
      if (theta_pol .lt. 0.d0) theta_pol = theta_pol + 2*PI
      if (theta_pol .gt. 2*PI) theta_pol = theta_pol - 2*PI
      if (theta_pol .lt. PI) then 
        theta_bal = q95*mod_n*theta_pol       *(ES%R_axis/R**2)
      else
        theta_bal = q95*mod_n*(2*PI-theta_pol)*(ES%R_axis/R**2)
      endif
      perturbation = gauss_amp * cos(theta_bal) * abs(R-R_left)
      node_list%node(i)%values(in,1,4) = perturbation

      node_list%node(i)%deltas = 0.d0

    enddo

  endif

  call Poisson(my_id,1,node_list,element_list,bnd_node_list,bnd_elm_list, &
               4,2,in, ES%psi_axis,ES%psi_bnd,xpoint2, xcase2,ES%Z_xpoint,freeboundary_equil,refinement,1)
enddo

!----------------------------------- fill in parallel velocity at boundary (on open field lines)
do i=1,node_list%n_nodes

  if ((node_list%node(i)%boundary .eq. 1) .or. (node_list%node(i)%boundary .eq. 3)) then

    ps0_s     = node_list%node(i)%values(1,2,1)
    ps0_t     = node_list%node(i)%values(1,3,1)
    R_s       = node_list%node(i)%x(1,2,1)
    R_t       = node_list%node(i)%x(1,3,1)
    Z_s       = node_list%node(i)%x(1,2,2)
    Z_t       = node_list%node(i)%x(1,3,2)

    xjac  =  R_s*Z_t - R_t*Z_s
    ps0_x = (	Z_t * ps0_s - Z_s * ps0_t ) / xjac
    ps0_y = ( - R_t * ps0_s + R_s * ps0_t ) / xjac

    if (tokamak_device(1:4) .eq. 'MAST') then
      if ( (node_list%node(i)%x(1,1,1) .gt. (ES%R_xpoint(1)+ES%R_xpoint(2))/2.d0) ) then
        direction = 1.d0
      else
        direction = -1.d0
      endif
    else
      direction = + ps0_x / abs(ps0_x)		 ! temporary solution for lower x-point only
    endif
    if (xcase2 .eq. 2) direction = -direction
    if ( (xcase2 .eq. 3) .and. (node_list%node(i)%x(1,1,2) .gt. (ES%Z_xpoint(1)+ES%Z_xpoint(2))/2.d0) ) direction = -direction

    BigR = node_list%node(i)%x(1,1,1)
    Btot = sqrt(F0**2 + ps0_x**2 + ps0_y**2) / BigR
    BigR_s = node_list%node(i)%x(1,2,1)
    
    do in=1,n_tor

      Ti0   = node_list%node(i)%values(in,1,6)
      Te0   = node_list%node(i)%values(in,1,8)
      node_list%node(i)%values(in,1,7) = direction / Btot * sqrt(GAMMA * (Ti0 + Te0))

      Ti0_s   = node_list%node(i)%values(in,2,6)
      Te0_s   = node_list%node(i)%values(in,2,8)
      node_list%node(i)%values(in,2,7) = BigR_s / (BigR*Btot) * sqrt(GAMMA * (Ti0 + Te0)) + 0.5d0 / Btot * sqrt(GAMMA / (Ti0 + Te0)) * (Ti0_s + Te0_s)
      node_list%node(i)%values(in,2,7) = direction *  node_list%node(i)%values(in,2,7)

      if(xcase2 .eq. 1) then
        write(*,'(A,8e14.6)') ' Boundary condition (eq): ',BigR,ES%psi_xpoint(1),node_list%node(i)%values(1,1,1),ps0_x,ps0_y, &
			    node_list%node(i)%values(in,1,n_var),BigR/F0 * sqrt(GAMMA*(Ti0 + Te0))
      endif
      if( (xcase2 .eq. 2) .or. ( (xcase2 .eq. 3) .and. (ES%psi_xpoint(2) .lt. ES%psi_xpoint(1)) ) ) then
        write(*,'(A,8e14.6)') ' Boundary condition (eq): ',BigR,ES%psi_xpoint(2),node_list%node(i)%values(1,1,1),ps0_x,ps0_y, &
			    node_list%node(i)%values(in,1,n_var),BigR/F0 * sqrt(GAMMA*(Ti0 + Te0))
      endif

    enddo
  endif
enddo

!call add_pellet(node_list,element_list,5.,0.08,0.03,R_geo-0.78,Z_geo)
!call add_pellet(node_list,element_list,25.,0.08,0.03,R_geo+0.85,Z_geo)

return
end
