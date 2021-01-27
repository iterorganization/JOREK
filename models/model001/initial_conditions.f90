subroutine initial_conditions(my_id,node_list,element_list,bnd_node_list, bnd_elm_list, xpoint2, xcase2)
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
use data_structure
use phys_module
use mod_poiss
implicit none

type (type_node_list)        :: node_list
type (type_element_list)     :: element_list
type (type_surface_list)     :: surface_list
type (type_bnd_node_list)    :: bnd_node_list
type (type_bnd_element_list) :: bnd_elm_list

integer    :: my_id, i , xcase2
real*8     :: amplitude1, amplitude2, R, Z, R1, Z1, R2, Z2, radius1, sigma1, radius2, sigma2, psi_bnd, psi_axis, Z_xpoint(2)
real*8     :: W, W_R, W_Z, W_RR, W_RZ, W_ZZ, R_s, R_t, R_st, Z_s, Z_t, Z_st
logical    :: xpoint2

if (my_id .eq. 0) then
  write(*,*) '***************************************'
  write(*,*) '*      initial conditions  (001)      *'
  write(*,*) '***************************************'
endif

!---------------------------- initialise perturbations

amplitude1 = 1.d0
amplitude2 = 1.d0
sigma1     = 0.64d0
sigma2     = 0.64d0
R1 = 1.5d0
Z1 = 0.d0
R2 = -1.5d0
Z2 = 0.d0

if (my_id .eq. 0) then

  do i=1,node_list%n_nodes

    node_list%node(i)%values(1,:,:) = 0.d0

    R   = node_list%node(i)%x(1,1,1)
    Z   = node_list%node(i)%x(1,1,2)

    radius1 = (R-R1)*(R-R1) + (Z-Z1)*(Z-Z1)
    radius2 = (R-R2)*(R-R2) + (Z-Z2)*(Z-Z2)
    
    W   =   amplitude1 * exp(-radius1/sigma1) +  amplitude2 * exp(-radius2/sigma2)
    W_R = - amplitude1 * 2.*(R-R1)/sigma1 * exp(-radius1/sigma1) &
          - amplitude2 * 2.*(R-R2)/sigma2 * exp(-radius2/sigma2)
    W_Z = - amplitude1 * 2.*(Z-Z1)/sigma1 * exp(-radius1/sigma1)  &
          - amplitude2 * 2.*(Z-Z2)/sigma2 * exp(-radius2/sigma2)
    W_RZ = + amplitude1 * 2.*(Z-Z1)/sigma1 * 2.*(R-R1)/sigma1 * exp(-radius1/sigma1)  &
           + amplitude2 * 2.*(Z-Z2)/sigma2 * 2.*(R-R2)/sigma2 * exp(-radius2/sigma2) 
	   
    W_RR =  - amplitude1 * 2./sigma1 * exp(-radius1/sigma1) + amplitude1 * 2.*(R-R1)/sigma1 * 2.*(R-R1)/sigma1 * exp(-radius1/sigma1) &
            - amplitude2 * 2./sigma2 * exp(-radius2/sigma2) + amplitude2 * 2.*(R-R2)/sigma2 * 2.*(R-R2)/sigma2 * exp(-radius2/sigma2)
	    
    W_ZZ =  - amplitude1 * 2./sigma1 * exp(-radius1/sigma1) + amplitude1 * 2.*(Z-Z1)/sigma1 * 2.*(Z-Z1)/sigma1 * exp(-radius1/sigma1) &
            - amplitude2 * 2./sigma2 * exp(-radius2/sigma2) + amplitude2 * 2.*(Z-Z2)/sigma2 * 2.*(Z-Z2)/sigma2 * exp(-radius2/sigma2)
 
    R_s  = node_list%node(i)%x(1,2,1)
    R_t  = node_list%node(i)%x(1,3,1)
    R_st = node_list%node(i)%x(1,4,1)
    Z_s  = node_list%node(i)%x(1,2,2)
    Z_t  = node_list%node(i)%x(1,3,2)
    Z_st = node_list%node(i)%x(1,4,2)
        
    node_list%node(i)%values(1,1,2) = W
    node_list%node(i)%values(1,2,2) = W_R * R_s + W_z * Z_s                                      
    node_list%node(i)%values(1,3,2) = W_R * R_t + W_z * Z_t 
    node_list%node(i)%values(1,4,2) = W_RR * R_s*R_t + W_R * R_st + W_RZ * (R_s*Z_t+R_t*Z_s) + W_ZZ*Z_s*Z_t + W_Z*Z_st

    node_list%node(i)%deltas = 0.d0
    
  enddo

  call Poisson(my_id,3,node_list,element_list,bnd_node_list,bnd_elm_list, &
               2,1,1, psi_axis,psi_bnd,xpoint2, xcase2,Z_xpoint,freeboundary_equil,refinement,1)

endif

return
end subroutine initial_conditions
