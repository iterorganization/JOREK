!> Add condition for the axis directly in the matrix.
!> This is aimed at stabilising numerical noise on the axis.
module mod_fix_axis_nodes
contains

subroutine fix_nodes_on_axis(node_list, element_list, local_elms, n_local_elms, index_min, index_max, a_mat)

  use mod_assembly, only : boundary_conditions_add_one_entry, boundary_conditions_add_RHS
  use data_structure
  use mod_locate_irn_jcn
  use mod_integer_types

  implicit none

  ! Subroutine parameters
  integer,                            intent(in)    :: local_elms(*)         !< List of local elements
  integer,                            intent(in)    :: n_local_elms          !< Number of local elements
  integer,                            intent(in)    :: index_min, index_max  !< Min/max index of local elements
  type (type_node_list),              intent(in)    :: node_list             !< List of nodes
  type (type_element_list),           intent(in)    :: element_list          !< List of all elements
  type(type_SP_MATRIX)                              :: a_mat
  ! Internal parameters
  real*8                :: zbig
  integer               :: i, in, iv, inode, k, ielm, ilarge2, n_tor_local
  integer               :: index_node, index_node2
  integer(kind=int_all) :: index_large_i
  integer(kind=int_all) :: ijA_position,ijA_position2

  n_tor_local = (a_mat%i_tor_max - a_mat%i_tor_min + 1)

  zbig = 1.d12
  do i=1, n_local_elms

    ielm = local_elms(i)

    do iv=1, n_vertex_max

      inode = element_list%element(ielm)%vertex(iv)

      if (node_list%node(inode)%axis_node) then

        do in=a_mat%i_tor_min, a_mat%i_tor_max
          do k=1, n_var

            ! --- For t-derivative
            index_node = node_list%node(inode)%index(3)
            if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then
              call locate_irn_jcn(index_node,index_node,index_min,index_max,ijA_position,a_mat)
              index_large_i = n_tor_local * n_var * (index_node - 1)
              ilarge2 = ijA_position - 1 + ((k-1)*n_tor_local + in-a_mat%i_tor_min) * n_var*n_tor_local &
                + (k-1)*n_tor_local + in - a_mat%i_tor_min + 1
              a_mat%irn(ilarge2) =  n_tor_local * n_var * (index_node-1) + (k-1)*n_tor_local + in - a_mat%i_tor_min + 1
              a_mat%jcn(ilarge2) =  n_tor_local * n_var * (index_node-1) + (k-1)*n_tor_local + in - a_mat%i_tor_min + 1
              a_mat%val(ilarge2)   = zbig
            end if

            ! --- For cross st-derivative
            index_node = node_list%node(inode)%index(4)
            if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then
              call locate_irn_jcn(index_node,index_node,index_min,index_max,ijA_position,a_mat)
              index_large_i = n_tor_local * n_var * (index_node - 1)
              ilarge2 = ijA_position - 1 + ((k-1)*n_tor_local + in-a_mat%i_tor_min) * n_var*n_tor_local &
                + (k-1)*n_tor_local + in - a_mat%i_tor_min + 1
              a_mat%irn(ilarge2) =  n_tor_local * n_var * (index_node-1) + (k-1)*n_tor_local + in - a_mat%i_tor_min + 1
              a_mat%jcn(ilarge2) =  n_tor_local * n_var * (index_node-1) + (k-1)*n_tor_local + in - a_mat%i_tor_min + 1
              a_mat%val(ilarge2)   = zbig
            end if

          enddo
        enddo

      endif

    enddo  ! n_vertex
  enddo ! n_elements

  return

end subroutine fix_nodes_on_axis

end module mod_fix_axis_nodes
