module mod_boundary_matrix_open
  implicit none
contains

subroutine boundary_matrix_open(vertex, direction, element, nodes, xpoint2, xcase2, R_axis, Z_axis, psi_axis, &
                                psi_bnd, R_xpoint, Z_xpoint, ELM, RHS, i_tor_min, i_tor_max)

!---------------------------------------------------------------------
! calculates the matrix contribution of the boundaries of one element
! implements the natural boundary conditions
!---------------------------------------------------------------------
use mod_parameters
use data_structure
use gauss
use basis_at_gaussian
use phys_module

implicit none

type (type_element)   :: element
type (type_node)      :: nodes(n_vertex_max)        ! the two nodes containing the boundary nodes
integer, intent(in)   :: i_tor_min   
integer, intent(in)   :: i_tor_max   

real*8     :: ELM(n_vertex_max*n_var*n_degrees*n_tor,n_vertex_max*n_var*n_degrees*n_tor)
real*8     :: RHS(n_vertex_max*n_var*n_degrees*n_tor)
real*8     :: psi_axis, R_axis, Z_axis, psi_bnd, R_xpoint(2), Z_xpoint(2)
integer    :: vertex(2), direction(2), xcase2
logical    :: xpoint2

return
end subroutine boundary_matrix_open
end module mod_boundary_matrix_open
