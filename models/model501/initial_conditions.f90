subroutine initial_conditions(my_id,node_list,element_list,bnd_node_list, bnd_elm_list, xpoint2, xcase2)
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
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
real*8     :: zrn, drn_dpsi, drn_dpsi2, drn_dz, drn_dz2, drn_dpsi_dz, drn_dpsi3, drn_dpsi2_dz, drn_dpsi_dz2
real*8     :: zT, dT_dpsi, dT_dpsi2, dT_dz, dT_dz2, dT_dpsi_dz, dT_dpsi3, dT_dpsi2_dz, dT_dpsi_dz2
real*8     :: zFFprime,dFFprime_dpsi,dFFprime_dz, dFFprime_dpsi_dz, dFFprime_dz2, dFFprime_dpsi2
real*8     :: R, Z, BigR, T0, BigR_s, T0_s
real*8     :: zjz, dj_dpsi, dj_dR, dj_dZ, dj_dR_dZ, dj_dR_DR, dj_dZ_dZ, dj_dpsi2, dj_dR_dpsi, dj_dZ_dpsi
real*8     :: zp, dp_dpsi, dp_dpsi2, dp_dz, dp_dz2, dp_dpsi_dz, P_ss, P_st, P_tt, R_out,Z_out,s_out,t_out
real*8     :: ps0_s, ps0_t, p_s, p_t, zj0_s, zj0_t,R_s, R_t, ps0_x, ps0_y, Z_s, Z_t, xjac, direction, Btot
logical    :: xpoint2
real*8     :: zV, dV_dpsi, dV_dpsi2, dV_dz, dV_dz2, dV_dpsi_dz, dV_dpsi3, dV_dpsi2_dz, dV_dpsi_dz2
real*8     :: Omega, dOmega_dpsi, dOmega_dz, dOmega_dpsi2, dOmega_dz2, dOmega_dpsi_dz

if (my_id .eq. 0) then
  write(*,*) '***************************************'
  write(*,*) '*      initial conditions  (501)      *'
  write(*,*) '***************************************'
endif

! --- This old projection method should be replaced by a direct solve on the elements
! --- It is much cleaner and more generic
if ( (my_id .eq. 0) .and. (n_order .le. 3) ) then

  do i=1,node_list%n_nodes

    psi = node_list%node(i)%values(1,1,var_psi)
    R   = node_list%node(i)%x(1,1,1)
    Z   = node_list%node(i)%x(1,1,2)
   

    call density(  xpoint2, xcase2, Z, ES%Z_xpoint, psi,ES%psi_axis,ES%psi_bnd,zn,dn_dpsi,dn_dz,dn_dpsi2,dn_dz2,      &
                                                               dn_dpsi_dz,dn_dpsi3,dn_dpsi_dz2, dn_dpsi2_dz)

    call temperature(xpoint2, xcase2, Z, ES%Z_xpoint, psi,ES%psi_axis,ES%psi_bnd,zT,dT_dpsi,dT_dz,dT_dpsi2,dT_dz2,    &
                                                               dT_dpsi_dz,dT_dpsi3,dT_dpsi_dz2, dT_dpsi2_dz)

    call neutral_density(xpoint2, xcase2, Z, ES%Z_xpoint, psi,ES%psi_axis,ES%psi_bnd,zrn,drn_dpsi,drn_dz,drn_dpsi2,drn_dz2,             &
                                                               drn_dpsi_dz,drn_dpsi3,drn_dpsi_dz2, drn_dpsi2_dz)

    call FFprime(   xpoint2, xcase2, Z, ES%Z_xpoint, psi,ES%psi_axis,ES%psi_bnd,zFFprime,dFFprime_dpsi,dFFprime_dz, &
                                                               dFFprime_dpsi2,dFFprime_dz2, dFFprime_dpsi_dz, .true.)

    if ( (abs(V_0) .ge. 1.d-19) .or. (num_rot) ) then
       if (normalized_velocity_profile) then
          call velocity(xpoint2, xcase2, Z, ES%Z_xpoint, psi,ES%psi_axis,ES%psi_bnd,zV,dV_dpsi,dV_dz,dV_dpsi2,dV_dz2, &
               dV_dpsi_dz,dV_dpsi3,dV_dpsi_dz2, dV_dpsi2_dz)
       else
          call velocity(xpoint2, xcase2, Z, ES%Z_xpoint, psi,ES%psi_axis,ES%psi_bnd,Omega,dOmega_dpsi,dOmega_dz,dOmega_dpsi2,dOmega_dz2, &
               dOmega_dpsi_dz,dV_dpsi3,dV_dpsi_dz2, dV_dpsi2_dz)
       endif
    endif

    zp         = zn * zT
    dp_dpsi    = zn * dT_dpsi + dn_dpsi * zT
    dp_dpsi2   = zn * dT_dpsi2 + 2.d0 * dn_dpsi * dT_dpsi + dn_dpsi2 * zT
    dp_dz      = zn * dT_dz + dn_dz * zT
    dp_dz2     = zn * dT_dz2 + 2.d0 * dn_dz * dT_dz + dn_dz2 * zT 							       
    dp_dpsi_dz = zn * dT_dpsi_dz + dn_dz * dT_dpsi + dn_dpsi * dT_dz + dn_dpsi_dz * zT

    node_list%node(i)%values(1,1,var_rho) = zn
    node_list%node(i)%values(1,2,var_rho) = dn_dpsi    * node_list%node(i)%values(1,2,var_psi) + dn_dz * node_list%node(i)%x(1,2,2)
    node_list%node(i)%values(1,3,var_rho) = dn_dpsi    * node_list%node(i)%values(1,3,var_psi) + dn_dz * node_list%node(i)%x(1,3,2)
    node_list%node(i)%values(1,4,var_rho) = dn_dpsi    * node_list%node(i)%values(1,4,var_psi) + dn_dz * node_list%node(i)%x(1,4,2) &
                                    + dn_dpsi2   * node_list%node(i)%values(1,2,var_psi) * node_list%node(i)%values(1,3,var_psi)  &
                                    + dn_dz2     * node_list%node(i)%x(1,2,2)        * node_list%node(i)%x(1,3,2)         &
                                    + dn_dpsi_dz * node_list%node(i)%values(1,3,var_psi) * node_list%node(i)%x(1,2,2)         &
                                    + dn_dpsi_dz * node_list%node(i)%values(1,2,var_psi) * node_list%node(i)%x(1,3,2)      

    node_list%node(i)%values(1,1,var_rhon) = zrn
    node_list%node(i)%values(1,2,var_rhon) = drn_dpsi    * node_list%node(i)%values(1,2,var_psi) + drn_dz * node_list%node(i)%x(1,2,2)
    node_list%node(i)%values(1,3,var_rhon) = drn_dpsi    * node_list%node(i)%values(1,3,var_psi) + drn_dz * node_list%node(i)%x(1,3,2)
    node_list%node(i)%values(1,4,var_rhon) = drn_dpsi    * node_list%node(i)%values(1,4,var_psi) + drn_dz * node_list%node(i)%x(1,4,2) &
                                    + drn_dpsi2   * node_list%node(i)%values(1,2,var_psi) * node_list%node(i)%values(1,3,var_psi)  &
                                    + drn_dz2     * node_list%node(i)%x(1,2,2)        * node_list%node(i)%x(1,3,2)         &
                                    + drn_dpsi_dz * node_list%node(i)%values(1,3,var_psi) * node_list%node(i)%x(1,2,2)         &
                                    + drn_dpsi_dz * node_list%node(i)%values(1,2,var_psi) * node_list%node(i)%x(1,3,2)

    node_list%node(i)%values(1,1,var_T) = zT
    node_list%node(i)%values(1,2,var_T) = dT_dpsi    * node_list%node(i)%values(1,2,var_psi) + dT_dz * node_list%node(i)%x(1,2,2)
    node_list%node(i)%values(1,3,var_T) = dT_dpsi    * node_list%node(i)%values(1,3,var_psi) + dT_dz * node_list%node(i)%x(1,3,2)
    node_list%node(i)%values(1,4,var_T) = dT_dpsi    * node_list%node(i)%values(1,4,var_psi) + dT_dz * node_list%node(i)%x(1,4,2) &
                                    + dT_dpsi2   * node_list%node(i)%values(1,2,var_psi) * node_list%node(i)%values(1,3,var_psi)  &
                                    + dT_dz2     * node_list%node(i)%x(1,2,2)        * node_list%node(i)%x(1,3,2)         &
                                    + dT_dpsi_dz * node_list%node(i)%values(1,3,var_psi) * node_list%node(i)%x(1,2,2)         &
                                    + dT_dpsi_dz * node_list%node(i)%values(1,2,var_psi) * node_list%node(i)%x(1,3,2)      

    node_list%node(i)%values(1,1,var_u) = - tauIC * zp 
    node_list%node(i)%values(1,2,var_u) = - tauIC * (dp_dpsi  * node_list%node(i)%values(1,2,var_psi) + dp_dz * node_list%node(i)%x(1,2,2))
    node_list%node(i)%values(1,3,var_u) = - tauIC * (dp_dpsi  * node_list%node(i)%values(1,3,var_psi) + dp_dz * node_list%node(i)%x(1,3,2))
    node_list%node(i)%values(1,4,var_u) = - tauIC * (dP_dpsi  * node_list%node(i)%values(1,4,var_psi) + dP_dz * node_list%node(i)%x(1,4,2) &
                                    + dP_dpsi2   * node_list%node(i)%values(1,2,var_psi) * node_list%node(i)%values(1,3,var_psi)  &
                                    + dP_dz2     * node_list%node(i)%x(1,2,2)              * node_list%node(i)%x(1,3,2)         &
                                    + dP_dpsi_dz * node_list%node(i)%values(1,3,var_psi) * node_list%node(i)%x(1,2,2)         &
                                    + dP_dpsi_dz * node_list%node(i)%values(1,2,var_psi) * node_list%node(i)%x(1,3,2) )
 
    node_list%node(i)%values(1,:,var_w) = 0.d0        ! vorticity (will be filled just below with inverse Poisson)

    node_list%node(i)%values(1,:,var_Vpar) = 0.d0        ! parallel velocity
    if ( (abs(V_0) .ge. 1.d-19) .or. (num_rot) ) then
      ! set Vpar,0 (JOREK normalized, ie without unit)= parallel velocity given as input profile 
      if (normalized_velocity_profile) then
        node_list%node(i)%values(1,1,var_Vpar) = zV
        node_list%node(i)%values(1,2,var_Vpar) = dV_dpsi    * node_list%node(i)%values(1,2,1) + dV_dz * node_list%node(i)%x(1,2,2)
        node_list%node(i)%values(1,3,var_Vpar) = dV_dpsi    * node_list%node(i)%values(1,3,1) + dV_dz * node_list%node(i)%x(1,3,2)
        node_list%node(i)%values(1,4,var_Vpar) = dV_dpsi    * node_list%node(i)%values(1,4,1) + dV_dz * node_list%node(i)%x(1,4,2) &
                                               + dV_dpsi2   * node_list%node(i)%values(1,2,1) * node_list%node(i)%values(1,3,1)    &
                                               + dV_dz2     * node_list%node(i)%x(1,2,2)      * node_list%node(i)%x(1,3,2)         &
                                               + dV_dpsi_dz * node_list%node(i)%values(1,3,1) * node_list%node(i)%x(1,2,2)         &
                                               + dV_dpsi_dz * node_list%node(i)%values(1,2,1) * node_list%node(i)%x(1,3,2) 
      ! set Vpar,0 = 2piR^2/F0 * Omega_tor,0
      ! where Omega_tor,0 is the angular toroidal (~parallel) rotation given as input profile.
      else
        node_list%node(i)%values(1,1,var_Vpar) = R**2 * Omega
        node_list%node(i)%values(1,2,var_Vpar) = 2.d0 * R * Omega       * node_list%node(i)%x(1,2,1)        &
                                                 + R**2   * dOmega_dpsi * node_list%node(i)%values(1,2,1)   &
                                                 + R**2   * dOmega_dz   * node_list%node(i)%x(1,2,2)
        node_list%node(i)%values(1,3,var_Vpar) = 2.d0 * R * Omega       * node_list%node(i)%x(1,3,1)        &
                                                 + R**2   * dOmega_dpsi * node_list%node(i)%values(1,3,1)   &
                                                 + R**2   * dOmega_dz   * node_list%node(i)%x(1,3,2)
      
        node_list%node(i)%values(1,4,var_Vpar) =   2.d0 *     node_list%node(i)%x(1,2,1)**2  * Omega &
                                                 + 2.d0 * R * node_list%node(i)%x(1,4,1)     * Omega &
                                                 + 2.d0 * R * node_list%node(i)%x(1,2,1)     * dOmega_dpsi  * node_list%node(i)%values(1,3,1) &
                                                 + 2.d0 * R * node_list%node(i)%x(1,2,1)     * dOmega_dz    * node_list%node(i)%x(1,3,2) 
        node_list%node(i)%values(1,4,var_Vpar) = node_list%node(i)%values(1,4,var_Vpar) &
                                                 + 2.d0 * R * dOmega_dpsi    * node_list%node(i)%x(1,3,1)      * node_list%node(i)%values(1,2,1) &
                                                 + R**2     * dOmega_dpsi2   * node_list%node(i)%values(1,3,1) * node_list%node(i)%values(1,2,1) &
                                                 + R**2     * dOmega_dpsi_dz * node_list%node(i)%x(1,3,2)      * node_list%node(i)%values(1,2,1) &
                                                 + R**2     * dOmega_dpsi    * node_list%node(i)%values(1,4,1)
        node_list%node(i)%values(1,4,var_Vpar) = node_list%node(i)%values(1,4,var_Vpar) &
                                                 + 2.d0 * R * dOmega_dz      * node_list%node(i)%x(1,3,1)      * node_list%node(i)%x(1,2,2) &
                                                 + R**2     * dOmega_dpsi_dz * node_list%node(i)%values(1,3,1) * node_list%node(i)%x(1,3,2) &
                                                 + R**2     * dOmega_dz2     * node_list%node(i)%x(1,3,2)      * node_list%node(i)%x(1,3,2) &
                                                 + R**2     * dOmega_dz      * node_list%node(i)%x(1,4,2)
      
        node_list%node(i)%values(1,1,var_Vpar) = 2.d0 * PI / F0 * node_list%node(i)%values(1,1,var_Vpar)
        node_list%node(i)%values(1,2,var_Vpar) = 2.d0 * PI / F0 * node_list%node(i)%values(1,2,var_Vpar)
        node_list%node(i)%values(1,3,var_Vpar) = 2.d0 * PI / F0 * node_list%node(i)%values(1,3,var_Vpar)
        node_list%node(i)%values(1,4,var_Vpar) = 2.d0 * PI / F0 * node_list%node(i)%values(1,4,var_Vpar)
      endif
    endif

    node_list%node(i)%deltas = 0.d0

  enddo

endif ! if n_order<=3

! --- Variable projection is better at higher order...
! --- (by the way, we could use this for n_order=3 and remove all the above as well, 
! --- and remove all derivatives from profiles functions, which are not really needed, 
! --- except dn_dpsi and dT_dpsi for current profile...)
if (n_order .ge. 5) then
  call Poisson(my_id,0,node_list,element_list,bnd_node_list,bnd_elm_list, &
               var_psi,var_rho,1, ES%psi_axis,ES%psi_bnd,xpoint2,xcase2,ES%Z_xpoint,freeboundary_equil,refinement,1)
  call Poisson(my_id,0,node_list,element_list,bnd_node_list,bnd_elm_list, &
               var_psi,var_T,1, ES%psi_axis,ES%psi_bnd,xpoint2,xcase2,ES%Z_xpoint,freeboundary_equil,refinement,1)
  call Poisson(my_id,0,node_list,element_list,bnd_node_list,bnd_elm_list, &
               var_psi,var_Vpar,1, ES%psi_axis,ES%psi_bnd,xpoint2,xcase2,ES%Z_xpoint,freeboundary_equil,refinement,1)
  call Poisson(my_id,0,node_list,element_list,bnd_node_list,bnd_elm_list, &
               var_psi,var_rhon,1, ES%psi_axis,ES%psi_bnd,xpoint2,xcase2,ES%Z_xpoint,freeboundary_equil,refinement,1)
endif



if (tauIC .ne. 0.d0) then
  call Poisson(my_id,2,node_list,element_list,bnd_node_list,bnd_elm_list, &
               var_u,var_w,1, ES%psi_axis,ES%psi_bnd,xpoint2, xcase2,ES%Z_xpoint,freeboundary_equil,refinement,1)      ! inverse Poisson
endif
    
!---------------------------- initialise perturbations
amplitude = 1.d-12
mm = 2

do in=2,n_tor

  if (my_id .eq. 0) then

    do i=1,node_list%n_nodes

      node_list%node(i)%values(in,:,:) = 0.d0

      psi = node_list%node(i)%values(1,1,var_psi)
      Z   = node_list%node(i)%x(1,1,2)
      psi_n = (psi - ES%psi_axis)/(ES%psi_bnd - ES%psi_axis)

      node_list%node(i)%values(in,1,var_w) = amplitude * psi_n * (1.d0 -psi_n)
      node_list%node(i)%values(in,2,var_w) = amplitude * (1. - 2.d0 * psi_n)/(ES%psi_bnd - ES%psi_axis) * node_list%node(i)%values(1,2,var_psi)
      node_list%node(i)%values(in,3,var_w) = amplitude * (1. - 2.d0 * psi_n)/(ES%psi_bnd - ES%psi_axis) * node_list%node(i)%values(1,3,var_psi)
      node_list%node(i)%values(in,4,var_w) = amplitude * (1. - 2.d0 * psi_n)/(ES%psi_bnd - ES%psi_axis) * node_list%node(i)%values(1,4,var_psi)
      
      if (xpoint2 .and. ((psi_n .gt. 1.d0) .or. ((Z .lt. ES%Z_xpoint(1)) .and. (xcase2 .ne. UPPER_XPOINT)) ) ) then
        node_list%node(i)%values(in,1:4,var_w) = 0.d0
      endif
      if (xpoint2 .and. ((psi_n .gt. 1.d0) .or. ((Z .gt. ES%Z_xpoint(2)) .and. (xcase2 .ne. LOWER_XPOINT)) ) ) then
        node_list%node(i)%values(in,1:4,var_w) = 0.d0
      endif

      node_list%node(i)%deltas = 0.d0

    enddo

  endif

  call Poisson(my_id,1,node_list,element_list,bnd_node_list,bnd_elm_list, &
               var_w,var_u,1, ES%psi_axis,ES%psi_bnd,xpoint2, xcase2,ES%Z_xpoint,freeboundary_equil,refinement,1)
enddo

return

! The following seems don't have any meaning since it is after the return, should we delete this
!----------------------------------- fill in parallel velocity at boundary (on open field lines)
if (.not. no_mach1_bc) then
  do i=1,node_list%n_nodes

#ifdef altcs
    node_list%node(i)%psi_eq(:) = node_list%node(i)%values(1,:,var_psi)
#endif

    if ((node_list%node(i)%boundary .eq. 1) .or. (node_list%node(i)%boundary .eq. 3)) then
 
      ps0_s     = node_list%node(i)%values(1,2,var_psi)
      ps0_t     = node_list%node(i)%values(1,3,var_psi)
      R_s       = node_list%node(i)%x(1,2,1)
      R_t       = node_list%node(i)%x(1,3,1)
      Z_s       = node_list%node(i)%x(1,2,2)
      Z_t       = node_list%node(i)%x(1,3,2)
 
      xjac  =  R_s*Z_t - R_t*Z_s
      ps0_x = (   Z_t * ps0_s - Z_s * ps0_t ) / xjac
      ps0_y = ( - R_t * ps0_s + R_s * ps0_t ) / xjac
 
      direction = + ps0_x / abs(ps0_x)		 ! temporary solution for lower x-point only
      if (xcase2 .eq. UPPER_XPOINT) direction = -direction
      if ( (xcase2 .eq. DOUBLE_NULL) .and. (node_list%node(i)%x(1,1,2) .gt. (ES%Z_xpoint(1)+ES%Z_xpoint(2))/2.d0) ) direction = -direction
      if ( (grid_to_wall) .and. (n_wall_blocks .ne. 0) ) direction = 0.d0 ! everything to zero for grid with patches
 
      BigR = node_list%node(i)%x(1,1,1)
      Btot = sqrt(F0**2 + ps0_x**2 + ps0_y**2) / BigR
      BigR_s = node_list%node(i)%x(1,2,1)
 
      T0   = node_list%node(i)%values(1,1,var_T)
      node_list%node(i)%values(1,1,var_Vpar) = direction / Btot * sqrt(GAMMA * T0)
 
      T0_s   = node_list%node(i)%values(1,2,var_T)
      node_list%node(i)%values(1,2,var_Vpar) = BigR_s / (BigR*Btot) * sqrt(GAMMA * T0) + 0.5d0 / Btot * sqrt(GAMMA / T0) * T0_s
      node_list%node(i)%values(1,2,var_Vpar) = direction *  node_list%node(i)%values(1,2,var_Vpar)
 
      if(xcase2 .eq. LOWER_XPOINT) then
        write(*,'(A,8e14.6)') ' Boundary condition (eq): ',BigR,ES%psi_xpoint(1),node_list%node(i)%values(1,1,var_psi),ps0_x,ps0_y, &
        		    node_list%node(i)%values(1,1,var_Vpar),BigR/F0 * sqrt(GAMMA*T0)
      endif
      if( ES%active_xpoint .eq. UPPER_XPOINT ) then
        write(*,'(A,8e14.6)') ' Boundary condition (eq): ',BigR,ES%psi_xpoint(2),node_list%node(i)%values(1,1,var_psi),ps0_x,ps0_y, &
        		    node_list%node(i)%values(1,1,var_Vpar),BigR/F0 * sqrt(GAMMA*T0)
      endif
 
    endif
  enddo
endif

return
end
