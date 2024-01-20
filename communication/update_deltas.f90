subroutine update_deltas(node_list, deltas)
!---------------------------------------------------------------------
! subroutine to create a local list of delta values
!---------------------------------------------------------------------
use tr_module
use data_structure, only: type_node_list, type_RHS
use mod_parameters, only: n_tor, n_var, n_degrees
!use global_distributed_matrix

implicit none
integer               :: i, j, k, in, index, index_node
type(type_node_list)  :: node_list
type(type_RHS)        :: deltas

if (.not. associated(deltas%val)) then
  call tr_allocatep(deltas%val,int1,node_list%n_dof,"deltas",CAT_FEM)
  deltas%n = node_list%n_dof
  deltas%val(1:deltas%n) = 0.d0
endif

do i = 1, node_list%n_nodes
 if ((.not. node_list%node(i)%constrained)) then
  do j=1,n_degrees

    index_node = node_list%node(i)%index(j)

    do k=1,n_var

      do in=1,n_tor

        index = n_tor*n_var * (index_node - 1) + n_tor*(k-1) + in

        deltas%val(index) = node_list%node(i)%deltas(in,j,k)

      enddo

    enddo

  enddo
 endif
enddo

return
end
