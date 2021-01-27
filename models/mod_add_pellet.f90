module mod_add_pellet
contains
subroutine add_pellet(node_list,element_list,pellet_density,pellet_size,sig_pellet,pellet_R,pellet_Z)
use constants
use tr_module
use data_structure

implicit none

type (type_node_list)    :: node_list
type (type_element_list) :: element_list
type (type_surface_list) :: surface_list

real*8 :: pellet_density,pellet_size,pellet_R,pellet_Z,sig_pellet
real*8 :: radius, dradius_dR, dradius_dZ, dradius_DRR, dradius_DRZ, dradius_DZZ
real*8 :: dR_ds, dR_dt, dR_dsdt, dZ_ds, dZ_dt, dZ_dsdt, dradius_ds, dradius_dt, dradius_dsdt
real*8 :: atn, datn, d2atn, delta_density, delta_density_ds, delta_density_dt, delta_density_dsdt
real*8 :: density_old, density_ds_old, density_dt_old, density_dsdt_old, T_old, T_ds_old, T_dt_old, T_dsdt_old   
integer :: i, in, nfft

real*8, allocatable :: delta_n_phi(:), T_phi(:), T_phi_s(:), T_phi_t(:), T_phi_st(:)
real*8              :: amplitude_n(n_tor), amplitude_T(n_tor), amplitude_T_s(n_tor), amplitude_T_t(n_tor), amplitude_T_st(n_tor), sum_A, phi

write(*,*) '*********************************************'
write(*,*) '* adding (non-axisymmetric) pellet density  *'
write(*,*) '*********************************************'

nfft=1024
call tr_allocate(delta_n_phi,1,nfft+2,"delta_n_phi",CAT_GRID)
call tr_allocate(T_phi,1,nfft+2,"T_phi",CAT_GRID)
call tr_allocate(T_phi_s,1,nfft+2,"T_phi_s",CAT_GRID)
call tr_allocate(T_phi_t,1,nfft+2,"T_phi_t",CAT_GRID)
call tr_allocate(T_phi_st,1,nfft+2,"T_phi_st",CAT_GRID)

do i=1,node_list%n_nodes

  radius = sqrt((node_list%node(i)%x(1,1,1)-pellet_R)**2 + (node_list%node(i)%x(1,1,2)-pellet_Z)**2 )
  
  dradius_dR = (node_list%node(i)%x(1,1,1)-pellet_R)/radius
  dradius_dZ = (node_list%node(i)%x(1,1,2)-pellet_Z)/radius
  
  dradius_DRR = -1d0 / radius**3d0 * (node_list%node(i)%x(1,1,1)-pellet_R)**2 + 1.d0/radius
  dradius_DRZ = -1d0 / radius**3d0 * (node_list%node(i)%x(1,1,2)-pellet_Z)**2 + 1.d0/radius
  dradius_DZZ = -1d0 / radius**3d0 * (node_list%node(i)%x(1,1,1)-pellet_R)*(node_list%node(i)%x(1,1,2)-pellet_Z) 
  
  dR_ds   = node_list%node(i)%x(1,2,1)
  dR_dt   = node_list%node(i)%x(1,3,1)
  dR_dsdt = node_list%node(i)%x(1,4,1)
  dZ_ds   = node_list%node(i)%x(1,2,2)
  dZ_dt   = node_list%node(i)%x(1,3,2)
  dZ_dsdt = node_list%node(i)%x(1,4,2)
  
  dradius_ds = dradius_dR * dR_ds + dradius_dZ * dZ_ds
  dradius_dt = dradius_dR * dR_dt + dradius_dZ * dZ_dt

  dradius_dsdt = dradius_dR  * dR_dsdt + dradius_dZ * dZ_dsdt           &
               + dradius_dRR * dR_ds*dR_dt + dradius_dZZ * dZ_ds*dZ_dt  &
	       + dradius_dRZ * (dR_ds*dZ_dt + dR_dt*dZ_ds)
  
  atn   = (0.5d0 - 0.5d0*tanh((radius - pellet_size)/sig_pellet))
  
  datn  = - 1.d0/cosh((radius - pellet_size)/sig_pellet)**2 / (2.d0 * sig_pellet) 
  d2atn =   1.d0/cosh((radius - pellet_size)/sig_pellet)**2 / (sig_pellet**2) * tanh((radius - pellet_size)/sig_pellet)

  do in=1,nfft
    phi = float(in-1)/float(nfft) * 2.d0 * PI
    
    delta_n_phi(in) = cos(phi/2.d0)**20 
    
    T_phi(in) = node_list%node(i)%values(1,1,6) * node_list%node(i)%values(1,1,5) &
              / ( node_list%node(i)%values(1,1,5) + pellet_density * atn * delta_n_phi(in))
	      
    T_phi_s(in) = (  node_list%node(i)%values(1,2,6) * node_list%node(i)%values(1,1,5)   &
                   + node_list%node(i)%values(1,1,6) * node_list%node(i)%values(1,2,5) ) &
                / ( node_list%node(i)%values(1,1,5) + pellet_density * atn * delta_n_phi(in)) &
	        - ( node_list%node(i)%values(1,1,6) * node_list%node(i)%values(1,1,5)) &
                 / ( node_list%node(i)%values(1,1,5) + pellet_density * atn * delta_n_phi(in))**2 &
	         * (node_list%node(i)%values(1,2,5)  + pellet_density * datn * dradius_ds * delta_n_phi(in))
		 
    T_phi_t(in) = (  node_list%node(i)%values(1,3,6) * node_list%node(i)%values(1,1,5)   &
                   + node_list%node(i)%values(1,1,6) * node_list%node(i)%values(1,3,5) ) &
                / ( node_list%node(i)%values(1,1,5) + pellet_density * atn * delta_n_phi(in)) &
	        - (node_list%node(i)%values(1,1,6) * node_list%node(i)%values(1,1,5)) &
                 / (node_list%node(i)%values(1,1,5) + pellet_density * atn * delta_n_phi(in))**2 &
	         * (node_list%node(i)%values(1,3,5) + pellet_density * datn * dradius_dt * delta_n_phi(in))

    T_phi_st(in) = ( node_list%node(i)%values(1,4,6) * node_list%node(i)%values(1,1,5)         &
                   + node_list%node(i)%values(1,1,6) * node_list%node(i)%values(1,4,5)         &
		   + node_list%node(i)%values(1,2,6) * node_list%node(i)%values(1,3,5)         &
		   + node_list%node(i)%values(1,3,6) * node_list%node(i)%values(1,2,5) )       &
		 / ( node_list%node(i)%values(1,1,5) + pellet_density * atn * delta_n_phi(in)) &
		 
                 - ( node_list%node(i)%values(1,3,6) * node_list%node(i)%values(1,1,5)            &
                   + node_list%node(i)%values(1,1,6) * node_list%node(i)%values(1,3,5) )          &
                 / ( node_list%node(i)%values(1,1,5) + pellet_density * atn * delta_n_phi(in))**2 &
	         * (node_list%node(i)%values(1,2,5)  + pellet_density * datn * dradius_ds * delta_n_phi(in)) &

                 - ( node_list%node(i)%values(1,2,6) * node_list%node(i)%values(1,1,5)            &
                   + node_list%node(i)%values(1,1,6) * node_list%node(i)%values(1,2,5) )          &
                 / ( node_list%node(i)%values(1,1,5) + pellet_density * atn * delta_n_phi(in))**2 &
	         * (node_list%node(i)%values(1,3,5)  + pellet_density * datn * dradius_dt * delta_n_phi(in)) &

                 + 2.d0* ( node_list%node(i)%values(1,1,6) * node_list%node(i)%values(1,1,5) )    &
                 / ( node_list%node(i)%values(1,1,5) + pellet_density * atn * delta_n_phi(in))**3 &
	         * (node_list%node(i)%values(1,2,5)  + pellet_density * datn * dradius_ds * delta_n_phi(in)) &
		 * (node_list%node(i)%values(1,3,5)  + pellet_density * datn * dradius_dt * delta_n_phi(in)) &

                 - ( node_list%node(i)%values(1,1,6) * node_list%node(i)%values(1,1,5) )          &
                 / ( node_list%node(i)%values(1,1,5) + pellet_density * atn * delta_n_phi(in))**2 &
	         * (node_list%node(i)%values(1,4,5)  + pellet_density * datn * dradius_dsdt * delta_n_phi(in)) 
				
  enddo
  
  call rft2(delta_n_phi,nfft,1)
  call rft2(T_phi,nfft,1)
  call rft2(T_phi_s,nfft,1)
  call rft2(T_phi_t,nfft,1)
  call rft2(T_phi_st,nfft,1)

  amplitude_n(1)    = delta_n_phi(1)/real(nfft)
  amplitude_T(1)    = T_phi(1)/real(nfft)
  amplitude_T_s(1)  = T_phi_s(1)/real(nfft)
  amplitude_T_t(1)  = T_phi_t(1)/real(nfft)
  amplitude_T_st(1) = T_phi_st(1)/real(nfft)
  do in=2, n_tor
    amplitude_n(in)    = 2.d0*delta_n_phi(in+1)/real(nfft)
    amplitude_T(in)    = 2.d0*T_phi(in+1)/real(nfft)
    amplitude_T_s(in)  = 2.d0*T_phi_s(in+1)/real(nfft)
    amplitude_T_t(in)  = 2.d0*T_phi_t(in+1)/real(nfft)
    amplitude_T_st(in) = 2.d0*T_phi_st(in+1)/real(nfft)
  enddo

  do in=1,n_tor

    delta_density      = amplitude_n(in) * pellet_density * atn
    delta_density_ds   = amplitude_n(in) * pellet_density * datn * dradius_ds   
    delta_density_dt   = amplitude_n(in) * pellet_density * datn * dradius_dt   
    delta_density_dsdt = amplitude_n(in) * pellet_density * datn * dradius_dsdt + d2atn * dradius_ds * dradius_dt  

    T_old      = node_list%node(i)%values(1,1,6)
    T_ds_old   = node_list%node(i)%values(1,2,6)
    T_dt_old   = node_list%node(i)%values(1,3,6)
    T_dsdt_old = node_list%node(i)%values(1,4,6)

    node_list%node(i)%values(in,1,5) = node_list%node(i)%values(in,1,5) + delta_density 
    node_list%node(i)%values(in,2,5) = node_list%node(i)%values(in,2,5) + delta_density_ds 
    node_list%node(i)%values(in,3,5) = node_list%node(i)%values(in,3,5) + delta_density_dt 
    node_list%node(i)%values(in,4,5) = node_list%node(i)%values(in,4,5) + delta_density_dsdt  
				  
    node_list%node(i)%values(in,1,6) = amplitude_T(in)
    node_list%node(i)%values(in,2,6) = amplitude_T_s(in)
    node_list%node(i)%values(in,3,6) = amplitude_T_t(in)
    node_list%node(i)%values(in,4,6) = amplitude_T_st(in)
                                    
  enddo
enddo

call tr_deallocate(delta_n_phi,"delta_n_phi",CAT_GRID)
call tr_deallocate(T_phi,"T_phi",CAT_GRID)
call tr_deallocate(T_phi_s,"T_phi_s",CAT_GRID)
call tr_deallocate(T_phi_t,"T_phi_t",CAT_GRID)
call tr_deallocate(T_phi_st,"T_phi_st",CAT_GRID)

return
end
end module mod_add_pellet
