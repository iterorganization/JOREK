subroutine initial_conditions(my_id,node_list,element_list,bnd_node_list, bnd_elm_list, xpoint2, xcase2)
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
use data_structure
use phys_module
use mod_poiss
use mod_interp, only: interp
use constants, only: mu_zero
implicit none

type (type_node_list)        :: node_list
type (type_element_list)     :: element_list
type (type_surface_list)     :: surface_list
type (type_bnd_node_list)    :: bnd_node_list
type (type_bnd_element_list) :: bnd_elm_list

integer    :: my_id, i, i2, j, in, mm, i_elm_axis, i_elm_xpoint(2), ifail, i_elm, xcase2
real*8     :: amplitude, psi, psi_axis
real*8     :: zn, dn_dpsi, dn_dpsi2, dn_dz, dn_dz2, dn_dpsi_dz, dn_dpsi3, dn_dpsi2_dz, dn_dpsi_dz2
real*8     :: zT, dT_dpsi, dT_dpsi2, dT_dz, dT_dz2, dT_dpsi_dz, dT_dpsi3, dT_dpsi2_dz, dT_dpsi_dz2
real*8     :: R_axis, Z_axis, s_axis, t_axis, R, Z, BigR
real*8     :: R_out, Z_out, s_out, t_out, R_lim, Z_lim, s_lim, t_lim, psi_lim
real*8     :: psi_n, psi_bnd
real*8     :: p_s, p_t, p_ss, p_st, p_tt
logical    :: xpoint2
real*8     :: s_norm, s_factor
real*8     :: Z_xpoint(2) =  (/ 0.0, 0.0 /)

if (my_id .eq. 0) then
  write(*,*) '***************************************'
  write(*,*) '*      initial conditions  (083)      *'
  write(*,*) '***************************************'
endif

if (my_id .eq. 0) then
  s_factor = 1.0 / (n_flux - 1)

  do i=1,node_list%n_nodes

    psi = node_list%node(i)%values(1,1,1)
    R   = node_list%node(i)%x(1,1,1)
    Z   = node_list%node(i)%x(1,1,2)
    s_norm = node_list%node(i)%s_eq(1)

    ! Initialise density
    call density(    xpoint2, xcase2, Z, Z_xpoint, s_norm, 0.0, 1.0, zn, dn_dpsi, dn_dz, dn_dpsi2, dn_dz2, dn_dpsi_dz, dn_dpsi3, dn_dpsi_dz2, dn_dpsi2_dz)
    node_list%node(i)%values(1,1,5) = zn
    node_list%node(i)%values(1,2,5) = dn_dpsi * s_factor * 1.0 / 3.0
    node_list%node(i)%values(1,3,5) = 0.d0
    node_list%node(i)%values(1,4,5) = 0.d0
    
    ! Initialise temperature based on density and pressure
    node_list%node(i)%values(1,1,6) = mu_zero*node_list%node(i)%pressure(1) / zn
    node_list%node(i)%values(1,2,6) = mu_zero * (node_list%node(i)%pressure(2) - node_list%node(i)%pressure(1) / zn * node_list%node(i)%values(1,2,5)) / zn   ! dT/ds = (dp/ds - dn/ds*T) / n  ... scale factor is already taken into account
    node_list%node(i)%values(1,3,6) = 0.d0
    node_list%node(i)%values(1,4,6) = 0.d0

  enddo

endif

!---------------------------- initialise perturbations

amplitude = 1.d-10
mm = 2

do in=2,n_tor

  if (my_id .eq. 0) then

    do i=1,node_list%n_nodes

      node_list%node(i)%values(in,:,:) = 0.d0

      node_list%node(i)%deltas = 0.d0

    enddo

  endif

enddo

#ifdef altcs
do i=1,node_list%n_nodes
  node_list%node(i)%psi_eq(:) = node_list%node(i)%values(1,:,1)
end do
#endif

!call add_pellet(node_list,element_list,50.,0.06,0.02,R_axis-0.96,Z_axis)

return
end
