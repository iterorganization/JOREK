subroutine initial_conditions(my_id,node_list,element_list,bnd_node_list, bnd_elm_list, xpoint2, xcase2)
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
use data_structure
use phys_module
use mod_poiss
use equil_info
use mod_interp, only: interp

implicit none

type (type_node_list)        :: node_list
type (type_element_list)     :: element_list
type (type_surface_list)     :: surface_list
type (type_bnd_node_list)    :: bnd_node_list
type (type_bnd_element_list) :: bnd_elm_list

integer    :: my_id, i, in, mm, i_elm_axis, i_elm_xpoint(2), ifail, i_elm, xcase2
real*8     :: amplitude, psi, psi_n
real*8     :: zn, dn_dpsi, dn_dpsi2, dn_dz, dn_dz2, dn_dpsi_dz, dn_dpsi3, dn_dpsi2_dz, dn_dpsi_dz2
real*8     :: zT, dT_dpsi, dT_dpsi2, dT_dz, dT_dz2, dT_dpsi_dz, dT_dpsi3, dT_dpsi2_dz, dT_dpsi_dz2
real*8     :: zFFprime,dFFprime_dpsi,dFFprime_dz, dFFprime_dpsi_dz, dFFprime_dz2, dFFprime_dpsi2
real*8     :: R, Z, BigR
real*8     :: R_out, Z_out, s_out, t_out, R_lim, Z_lim, s_lim, t_lim, psi_lim
real*8     :: p_s, p_t, p_ss, p_st, p_tt
logical    :: xpoint2

if (my_id .eq. 0) then
  write(*,*) '***************************************'
  write(*,*) '*      initial conditions  (199)      *'
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

    call temperature(xpoint2, xcase2, Z, ES%Z_xpoint, psi,ES%psi_axis,ES%psi_bnd,zT,dT_dpsi,dT_dz,dT_dpsi2,dT_dz2,             &
                                                               dT_dpsi_dz,dT_dpsi3,dT_dpsi_dz2, dT_dpsi2_dz)

    call FFprime(    xpoint2, xcase2, Z, ES%Z_xpoint, psi,ES%psi_axis,ES%psi_bnd,zFFprime,dFFprime_dpsi,dFFprime_dz, &
                                                               dFFprime_dpsi2,dFFprime_dz2, dFFprime_dpsi_dz, .true.)

    node_list%node(i)%values(1,1,5) = zn
    node_list%node(i)%values(1,2,5) = dn_dpsi  * node_list%node(i)%values(1,2,1) + dn_dz * node_list%node(i)%x(1,2,2)
    node_list%node(i)%values(1,3,5) = dn_dpsi  * node_list%node(i)%values(1,3,1) + dn_dz * node_list%node(i)%x(1,3,2)
    node_list%node(i)%values(1,4,5) = dn_dpsi  * node_list%node(i)%values(1,4,1) + dn_dz * node_list%node(i)%x(1,4,2) &
                                    + dn_dpsi2 * node_list%node(i)%values(1,2,1) * node_list%node(i)%values(1,3,1)  &
                                    + dn_dz2   * node_list%node(i)%x(1,2,2)        * node_list%node(i)%x(1,3,2) &
                                    + dn_dpsi_dz * ( &
                                      + node_list%node(i)%values(1,3,1) * node_list%node(i)%x(1,2,2) &
                                      + node_list%node(i)%values(1,2,1) * node_list%node(i)%x(1,3,2) &
                                      )

    node_list%node(i)%values(1,1,6) = zT
    node_list%node(i)%values(1,2,6) = dT_dpsi  * node_list%node(i)%values(1,2,1) + dT_dz * node_list%node(i)%x(1,2,2)
    node_list%node(i)%values(1,3,6) = dT_dpsi  * node_list%node(i)%values(1,3,1) + dT_dz * node_list%node(i)%x(1,3,2)
    node_list%node(i)%values(1,4,6) = dT_dpsi  * node_list%node(i)%values(1,4,1) + dT_dz * node_list%node(i)%x(1,4,2) &
                                    + dT_dpsi2 * node_list%node(i)%values(1,2,1) * node_list%node(i)%values(1,3,1)  &
                                    + dT_dz2   * node_list%node(i)%x(1,2,2)        * node_list%node(i)%x(1,3,2) &
                                    + dT_dpsi_dz * ( &
                                      + node_list%node(i)%values(1,3,1) * node_list%node(i)%x(1,2,2) &
                                      + node_list%node(i)%values(1,2,1) * node_list%node(i)%x(1,3,2) &
                                      )

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

amplitude = 1.d-10
mm = 2

do in=2,n_tor

  if (my_id .eq. 0) then

    do i=1,node_list%n_nodes

      node_list%node(i)%values(in,:,:) = 0.d0

      psi = node_list%node(i)%values(1,1,1)
      Z   = node_list%node(i)%x(1,1,2)

      psi_n = (psi - ES%psi_axis)/(ES%psi_bnd - ES%psi_axis)

      node_list%node(i)%values(in,1,4) = amplitude * psi_n * (1.d0 -psi_n)
      node_list%node(i)%values(in,2,4) = amplitude * (1. - 2.d0 * psi_n)/(ES%psi_bnd - ES%psi_axis) * node_list%node(i)%values(1,2,1)
      node_list%node(i)%values(in,3,4) = amplitude * (1. - 2.d0 * psi_n)/(ES%psi_bnd - ES%psi_axis) * node_list%node(i)%values(1,3,1)
      node_list%node(i)%values(in,4,4) = amplitude * (1. - 2.d0 * psi_n)/(ES%psi_bnd - ES%psi_axis) * node_list%node(i)%values(1,4,1)

      if (xpoint2 .and. ((psi_n .gt. 1.d0) .or. ((Z .lt. ES%Z_xpoint(1)) .and. (xcase2 .ne. UPPER_XPOINT)) ) ) then
        node_list%node(i)%values(in,1:4,4) = 0.d0
      endif
      if (xpoint2 .and. ((psi_n .gt. 1.d0) .or. ((Z .gt. ES%Z_xpoint(2)) .and. (xcase2 .ne. LOWER_XPOINT)) ) ) then
        node_list%node(i)%values(in,1:4,4) = 0.d0
      endif

      node_list%node(i)%deltas = 0.d0

    enddo

  endif

  call Poisson(my_id,1,node_list,element_list,bnd_node_list,bnd_elm_list, &
               4,2,in, ES%psi_axis,ES%psi_bnd,xpoint2, xcase2,ES%Z_xpoint,freeboundary_equil,refinement,1)

enddo


!call add_pellet(node_list,element_list,50.,0.06,0.02,ES%R_axis-0.96,ES%Z_axis)

return
end
