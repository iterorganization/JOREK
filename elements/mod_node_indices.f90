!> This module takes care of the relation between the 2D-indexing of our FEM vectorial basis,
!> (which is directly related to Bezier control points and (s,t)-derivatives)
!> and the matrix indexing, defined in the structure
!> node_list%node(i_node)%index(i_index)
!> With bi-cubic elements, ordering of the matrix index at each node is trivial, 1-4 nodes.
!> At each node, the 4 control points, with (s,t)-indices (1,1), (1,2), (2,1), (2,2)
!> correspond to matrix indices
!>
!>              (s,t)-indices                             matrix-index                                  (s,t)-derivatives
!>                                                                                             
!>      (1,2) --------------- (2,2)                 (3) -------------------- (4)                dpsi/dt ------------------ d2psi/dsdt
!>        |                     |                    |                        |                    |                        |
!>        |                     |                    |                        |                    |                        |
!>        |                     |        ==>         |                        |        ==>         |                        |
!>        |                     |                    |                        |                    |                        |
!>      (1,1) --------------- (2,1)                 (1) -------------------- (2)                  psi -------------------- dpsi/ds
!>
!> With n_order=5 and higher, this ordering quickly becomes very complex
!> This module contains a set of routines to find the matrix-index given the (s,t)-indices,
!> and vice versa, the (s,t)-indices given the matrix index.
module mod_node_indices

use mod_parameters

contains

!> calculate the node matrix indices as a function of (s,t)-coord indices
subroutine calculate_node_indices(node_indices)

  implicit none

  integer, intent(inout) :: node_indices((n_order+1)/2,(n_order+1)/2)
  integer                :: k, l, i, j, counter

  node_indices = 0
  counter = 1
  ! --- We do it square by square
  do k = 1,(n_order+1)/2
    ! --- For each square, we do i row and j column
    do l = 1,k
        ! --- First the row
        i = k
        j = l
        node_indices(i,j) = counter
        counter = counter + 1
        ! --- Then the column
        i = l
        j = k
        if (i .eq. j) cycle ! don't record corners twice
        node_indices(i,j) = counter
        counter = counter + 1
    enddo
  enddo

  return
end subroutine calculate_node_indices


!> given the node matrix index, find the (s,t)-coord indices of a node
subroutine get_node_coords_from_index(node_indices, index, i, j)

  implicit none

  integer, intent(in)    :: node_indices((n_order+1)/2,(n_order+1)/2), index
  integer, intent(inout) :: i, j
  integer                :: k, l
  logical                :: found

  i = 0 ; j = 0
  found = .false.
  do k = 1,(n_order+1)/2
    do l = 1,(n_order+1)/2
      if (node_indices(k,l) .eq. index) then
        i = k
        j = l
        found = .true.
        exit
      endif
    enddo
    if (found) exit
  enddo

  return
end subroutine get_node_coords_from_index





end module mod_node_indices
