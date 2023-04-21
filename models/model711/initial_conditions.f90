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
use mod_F_profile

implicit none

type (type_node_list)    :: node_list
type (type_element_list) :: element_list
type (type_surface_list) :: surface_list
type (type_bnd_node_list)    :: bnd_node_list
type (type_bnd_element_list) :: bnd_elm_list

integer    :: my_id, i, in, mm, i_elm, ifail, xcase2
integer    :: index0, index, n_node_start, n_index_start, j, k, ivar
real*8     :: amplitude, psi, psi_n, theta
real*8     :: zn, dn_dpsi, dn_dpsi2, dn_dz, dn_dz2, dn_dpsi_dz, dn_dpsi3, dn_dpsi2_dz, dn_dpsi_dz2
real*8     :: zTi, dTi_dpsi, dTi_dpsi2, dTi_dz, dTi_dz2, dTi_dpsi_dz, dTi_dpsi3, dTi_dpsi2_dz, dTi_dpsi_dz2
real*8     :: zTe, dTe_dpsi, dTe_dpsi2, dTe_dz, dTe_dz2, dTe_dpsi_dz, dTe_dpsi3, dTe_dpsi2_dz, dTe_dpsi_dz2
real*8     :: R, Z, BigR, T0, BigR_s, T0_s
real*8     :: zjz, dj_dpsi, dj_dR, dj_dZ, dj_dR_dZ, dj_dR_DR, dj_dZ_dZ, dj_dpsi2, dj_dR_dpsi, dj_dZ_dpsi
real*8     :: P_ss, P_st, P_tt, R_out,Z_out,s_out,t_out 
real*8     :: ps0_s, ps0_t, p_s, p_t, zj0_s, zj0_t,R_s, R_t, ps0_x, ps0_y, Z_s, Z_t, xjac, direction, Btot
real*8     :: Omega, dOmega_dpsi, dOmega_dpsi2, zeta, Lam, dLam_dpsi, dLam_dpsi2
real*8     :: zn0, zT0, dn0_dpsi, dT0_dpsi, dn_dR, dn_dR2, dT_dR, dT_dR2, R2sh, rf, rf0
real*8     :: x21, x31, x41, psi2, psi3, psi4
logical    :: xpoint2
real*8     :: F_prof,          dF_dpsi,  dF_dz,  dF_dpsi2,  dF_dz2,  dF_dpsi_dz
real*8     :: FFprime_profile, dFF_dpsi, dFF_dz, dFF_dpsi2, dFF_dz2, dFF_dpsi_dz



if (my_id .eq. 0) then
  write(*,*) '***************************************'
  write(*,*) '*      initial conditions  (711)      *'
  write(*,*) '***************************************'
endif

! --- This old projection method should be replaced by a direct solve on the elements
! --- It is much cleaner and more generic
if ( (my_id .eq. 0) .and. (n_order .le. 3) ) then

  do i=1,node_list%n_nodes

    psi = node_list%node(i)%values(1,1,1)
    R   = node_list%node(i)%x(1,1,1)
    Z   = node_list%node(i)%x(1,1,2)

    call density(    xpoint2, xcase2, Z, ES%Z_xpoint, psi,ES%psi_axis,ES%psi_bnd,zn,dn_dpsi,dn_dz,dn_dpsi2,dn_dz2,             &
                                                               dn_dpsi_dz,dn_dpsi3,dn_dpsi_dz2, dn_dpsi2_dz)

    call temperature_i(xpoint2, xcase2, Z, ES%Z_xpoint, psi,ES%psi_axis,ES%psi_bnd,zTi,dTi_dpsi,dTi_dz,dTi_dpsi2,dTi_dz2, &
                                                               dTi_dpsi_dz,dTi_dpsi3,dTi_dpsi_dz2, dTi_dpsi2_dz)

    call temperature_e(xpoint2, xcase2, Z, ES%Z_xpoint, psi,ES%psi_axis,ES%psi_bnd,zTe,dTe_dpsi,dTe_dz,dTe_dpsi2,dTe_dz2, &
                                                               dTe_dpsi_dz,dTe_dpsi3,dTe_dpsi_dz2, dTe_dpsi2_dz)


    node_list%node(i)%values(1,:,var_uR) = 0.d0    
    node_list%node(i)%values(1,:,var_uZ) = 0.d0

    node_list%node(i)%values(1,:,var_AR) = 0.d0    
    node_list%node(i)%values(1,:,var_AZ) = 0.d0

    node_list%node(i)%values(1,1,var_rho) = zn
    node_list%node(i)%values(1,2,var_rho) = dn_dpsi    * node_list%node(i)%values(1,2,var_A3) + dn_dz * node_list%node(i)%x(1,2,2)
    node_list%node(i)%values(1,3,var_rho) = dn_dpsi    * node_list%node(i)%values(1,3,var_A3) + dn_dz * node_list%node(i)%x(1,3,2)
    node_list%node(i)%values(1,4,var_rho) = dn_dpsi    * node_list%node(i)%values(1,4,var_A3) + dn_dz * node_list%node(i)%x(1,4,2) &
                                          + dn_dpsi2   * node_list%node(i)%values(1,2,var_A3) * node_list%node(i)%values(1,3,var_A3)  &
                                          + dn_dz2     * node_list%node(i)%x(1,2,2)             * node_list%node(i)%x(1,3,2)         &
                                          + dn_dpsi_dz * node_list%node(i)%values(1,3,var_A3) * node_list%node(i)%x(1,2,2)         &
                                          + dn_dpsi_dz * node_list%node(i)%values(1,2,var_A3) * node_list%node(i)%x(1,3,2)      


    node_list%node(i)%values(1,1,var_Ti) = zTi
    node_list%node(i)%values(1,2,var_Ti) = dTi_dpsi    * node_list%node(i)%values(1,2,var_A3) + dTi_dz * node_list%node(i)%x(1,2,2)
    node_list%node(i)%values(1,3,var_Ti) = dTi_dpsi    * node_list%node(i)%values(1,3,var_A3) + dTi_dz * node_list%node(i)%x(1,3,2)
    node_list%node(i)%values(1,4,var_Ti) = dTi_dpsi    * node_list%node(i)%values(1,4,var_A3) + dTi_dz * node_list%node(i)%x(1,4,2) &
                                         + dTi_dpsi2   * node_list%node(i)%values(1,2,var_A3) * node_list%node(i)%values(1,3,var_A3)  &
                                         + dTi_dz2     * node_list%node(i)%x(1,2,2)             * node_list%node(i)%x(1,3,2)         &
                                         + dTi_dpsi_dz * node_list%node(i)%values(1,3,var_A3) * node_list%node(i)%x(1,2,2)         &
                                         + dTi_dpsi_dz * node_list%node(i)%values(1,2,var_A3) * node_list%node(i)%x(1,3,2)      
    
    node_list%node(i)%values(1,1,var_Te) = zTe
    node_list%node(i)%values(1,2,var_Te) = dTe_dpsi    * node_list%node(i)%values(1,2,var_A3) + dTe_dz * node_list%node(i)%x(1,2,2)
    node_list%node(i)%values(1,3,var_Te) = dTe_dpsi    * node_list%node(i)%values(1,3,var_A3) + dTe_dz * node_list%node(i)%x(1,3,2)
    node_list%node(i)%values(1,4,var_Te) = dTe_dpsi    * node_list%node(i)%values(1,4,var_A3) + dTe_dz * node_list%node(i)%x(1,4,2) &
                                         + dTe_dpsi2   * node_list%node(i)%values(1,2,var_A3) * node_list%node(i)%values(1,3,var_A3)  &
                                         + dTe_dz2     * node_list%node(i)%x(1,2,2)             * node_list%node(i)%x(1,3,2)         &
                                         + dTe_dpsi_dz * node_list%node(i)%values(1,3,var_A3) * node_list%node(i)%x(1,2,2)         &
                                         + dTe_dpsi_dz * node_list%node(i)%values(1,2,var_A3) * node_list%node(i)%x(1,3,2)      
    
    node_list%node(i)%values(1,:,var_up) = 0.d0

    node_list%node(i)%deltas = 0.d0

    node_list%node(i)%psi_eq(:) = node_list%node(i)%values(1,:,1)

    ! Fprof_eq was aleady initialised in equilibrium.f90. 
    ! We fill in the values here as well, but anyway, we solve Fprof = Fprof below to ensure that the node values are clean
    ! This makes it 100% certain that all derivatives of Fprofile (when taken from the node values), will be accurate
    ! to the level of our finite elements.
    call F_profile(xpoint2, xcase2, Z, ES%Z_xpoint, psi, ES%psi_axis, ES%psi_bnd, &
                   F_prof,          dF_dpsi,  dF_dz,  dF_dpsi2,  dF_dz2,  dF_dpsi_dz , &
                   FFprime_profile, dFF_dpsi, dFF_dz, dFF_dpsi2, dFF_dz2, dFF_dpsi_dz)
    node_list%node(i)%Fprof_eq(1) =   F_prof
    node_list%node(i)%Fprof_eq(2) =   dF_dpsi  * node_list%node(i)%values(1,2,var_A3) + dF_dz * node_list%node(i)%x(1,2,2)
    node_list%node(i)%Fprof_eq(3) =   dF_dpsi  * node_list%node(i)%values(1,3,var_A3) + dF_dz * node_list%node(i)%x(1,3,2)
    node_list%node(i)%Fprof_eq(4) = dF_dpsi    * node_list%node(i)%values(1,4,var_A3) + dF_dz * node_list%node(i)%x(1,4,2) &
                                  + dF_dpsi2   * node_list%node(i)%values(1,2,var_A3) * node_list%node(i)%values(1,3,var_A3)  &
                                  + dF_dz2     * node_list%node(i)%x(1,2,2)             * node_list%node(i)%x(1,3,2)         &
                                  + dF_dpsi_dz * node_list%node(i)%values(1,3,var_A3) * node_list%node(i)%x(1,2,2)         &
                                  + dF_dpsi_dz * node_list%node(i)%values(1,2,var_A3) * node_list%node(i)%x(1,3,2)      

  enddo

endif ! if n_order<=3

! --- Variable projection is better at higher order...
! --- (by the way, we could use this for n_order=3 and remove all the above as well, 
! --- and remove all derivatives from profiles functions, which are not really needed, 
! --- except dn_dpsi and dT_dpsi for current profile...)
if (n_order .ge. 5) then
  call Poisson(my_id,0,node_list,element_list,bnd_node_list,bnd_elm_list, &
               var_A3,var_rho,1, ES%psi_axis,ES%psi_bnd,xpoint2,xcase2,ES%Z_xpoint,freeboundary_equil,refinement,1)
  call Poisson(my_id,0,node_list,element_list,bnd_node_list,bnd_elm_list, &
               var_A3,var_Ti,1, ES%psi_axis,ES%psi_bnd,xpoint2,xcase2,ES%Z_xpoint,freeboundary_equil,refinement,1)
  call Poisson(my_id,0,node_list,element_list,bnd_node_list,bnd_elm_list, &
               var_A3,var_Te,1, ES%psi_axis,ES%psi_bnd,xpoint2,xcase2,ES%Z_xpoint,freeboundary_equil,refinement,1)
endif



! --- This is the special Poisson for Fprofile (it will not overwrite var_A3)
call Poisson(my_id,0,node_list,element_list,bnd_node_list,bnd_elm_list, &
             var_A3,710,1, ES%psi_axis,ES%psi_bnd,xpoint2, xcase2,ES%Z_xpoint,freeboundary_equil,refinement,1)      ! inverse Poisson


!---------------------------- initialise perturbations
amplitude = 1.d-12
mm = 2

do in=2,n_tor

  if (my_id .eq. 0) then

    do i=1,node_list%n_nodes

      node_list%node(i)%values(in,:,:) = 0.d0

      psi = node_list%node(i)%values(1,1,1)
      Z   = node_list%node(i)%x(1,1,2)
      psi_n = (psi - ES%psi_axis)/(ES%psi_bnd - ES%psi_axis)

     
      ! Initialise perturbation for A3 nonzero n harmonics
      node_list%node(i)%values(in,:,:)= 0.d0

      node_list%node(i)%values(in,1,var_A3) = amplitude * psi_n * (1.d0 -psi_n)
      node_list%node(i)%values(in,2,var_A3) = amplitude * (1. - 2.d0 * psi_n)/(ES%psi_bnd - ES%psi_axis) * node_list%node(i)%values(1,2,var_A3)
      node_list%node(i)%values(in,3,var_A3) = amplitude * (1. - 2.d0 * psi_n)/(ES%psi_bnd - ES%psi_axis) * node_list%node(i)%values(1,3,var_A3)
      node_list%node(i)%values(in,4,var_A3) = amplitude * (1. - 2.d0 * psi_n)/(ES%psi_bnd - ES%psi_axis) * node_list%node(i)%values(1,4,var_A3)

      if (xpoint2 .and. ((psi_n .gt. 1.d0) .or. ((Z .lt. ES%Z_xpoint(1)) .and. (xcase2 .ne. UPPER_XPOINT)) ) ) then
        node_list%node(i)%values(in,1:4,4) = 0.d0
      endif
      if (xpoint2 .and. ((psi_n .gt. 1.d0) .or. ((Z .gt. ES%Z_xpoint(2)) .and. (xcase2 .ne. LOWER_XPOINT)) ) ) then
        node_list%node(i)%values(in,1:4,4) = 0.d0
      endif

      node_list%node(i)%deltas = 0.d0

    enddo

  endif

enddo




return
end
