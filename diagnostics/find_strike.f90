!> Locate the position of the strike points.
!!
!! Note: The routine assumes that psi_bnd in the equilibrium datastructure ES is updated already.
subroutine find_strike(node_list, bnd_elm_list)
  
  use data_structure
  use gauss
  use equil_info
  use basis_at_gaussian
  
  implicit none
  
  ! --- Routine parameters
  type(type_node_list),        intent(in)    :: node_list    !< List of grid nodes
  type(type_bnd_element_list), intent(in)    :: bnd_elm_list !< List of grid boundary elements
  
  ! --- Local variables
  integer :: ibndelm, iv1, iv2
  real*8  :: psi1, psi2, coeff, R1, R2, Z1, Z2
  
  ! ### The determination of the strike points should be implemented properly. This is just a rough
  ! ### estimate now. And ES%s_strike is not even calculated this way...
  
  do ibndelm = 1, bnd_elm_list%n_bnd_elements
    
    iv1  = bnd_elm_list%bnd_element(ibndelm)%vertex(1)
    iv2  = bnd_elm_list%bnd_element(ibndelm)%vertex(2)
    
    psi1 = node_list%node(iv1)%values(1,1,1)
    psi2 = node_list%node(iv2)%values(1,1,1)
    
    R1   = node_list%node(iv1)%x(1,1,1)
    Z1   = node_list%node(iv1)%x(1,1,2)
    R2   = node_list%node(iv2)%x(1,1,1)
    Z2   = node_list%node(iv2)%x(1,1,2)
    
    if ( (min(psi1,psi2) < ES%psi_bnd) .and. (max(psi1,psi2) > ES%psi_bnd) ) then
      ES%num_strike = ES%num_strike + 1
      coeff = (ES%psi_bnd - psi1)/(psi2 - psi1)
      ES%R_strike(ES%num_strike) = (1.d0-coeff) * R1 + coeff * R2
      ES%Z_strike(ES%num_strike) = (1.d0-coeff) * Z1 + coeff * Z2
      ES%i_bndelm_strike(ES%num_strike) = ibndelm
    end if
  end do
  
end subroutine find_strike
