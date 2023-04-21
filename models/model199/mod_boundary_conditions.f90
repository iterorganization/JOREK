!*******************************************************************************
!* Subroutine: boundary_condition                                              *
!*******************************************************************************
!*                                                                             *
!* Add boundary condition on the matrix.                                       *
!*                                                                             *
!* Parameters:                                                                 *
!*   my_id        - Identifier of the node in MPI_COMM_WORLD                   *
!*   node_list    - List of nodes                                              *
!*   element_list - List of all elements                                       *
!*   local_elms   - List of local elements                                     *
!*   n_local_elms - Number of local elements                                   *
!*   index_min    - Minimal index of local elements                            *
!*   index_max    - Maximal index of local elements                            *
!*   xpoint2      -                                                            *
!*   xcase2       -                                                            *
!*   psi_axis     -                                                            *
!*   psi_bnd      -                                                            *
!*   Z_xpoint     -                                                            *
!*                                                                             *
!*******************************************************************************
module mod_boundary_conditions
implicit none
contains
  subroutine boundary_conditions( my_id, node_list, element_list, bnd_node_list, local_elms,          &
                                  n_local_elms, index_min, index_max, rhs_loc, xpoint2, xcase2,       & 
                                  R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint, Z_xpoint, psi_xpoint,  &  
                                  a_mat)

    use data_structure
    use phys_module, only: F0, GAMMA, keep_n0_const
    use vacuum, only: is_freebound
    use mpi_mod
    use mod_integer_types
    use mod_assembly, only : boundary_conditions_add_one_entry, boundary_conditions_add_RHS
    use mod_node_indices

    implicit none

    ! --- Routine parameters
    integer,                            intent(in)    :: my_id
    type (type_node_list),              intent(in)    :: node_list
    type (type_element_list),           intent(in)    :: element_list
    type (type_bnd_node_list),          intent(in)    :: bnd_node_list
    integer,                            intent(in)    :: local_elms(*)
    integer,                            intent(in)    :: n_local_elms
    integer,                            intent(in)    :: index_min
    integer,                            intent(in)    :: index_max
    real*8,                             intent(inout) :: rhs_loc(*)
    logical,                            intent(in)    :: xpoint2
    integer,                            intent(in)    :: xcase2
    real*8,                             intent(in)    :: R_axis
    real*8,                             intent(in)    :: Z_axis
    real*8,                             intent(in)    :: psi_axis
    real*8,                             intent(in)    :: psi_bnd
    real*8,                             intent(in)    :: R_xpoint(2)
    real*8,                             intent(in)    :: Z_xpoint(2)
    real*8,                             intent(in)    :: psi_xpoint(2)
    type(type_SP_MATRIX)                              :: a_mat

    ! Internal parameters
    real*8                :: zbig, zbig_backup
    integer               :: i, in, iv, inode, k
    integer               :: ielm
    integer               :: index_node
    integer(kind=int_all) :: ijA_position
    integer               :: ilarge2
    integer               :: ierr, n_tor_local
    integer               :: node_indices( (n_order+1)/2, (n_order+1)/2 ), index_tmp, kk, ll, iv_dir

    ! --- calculate node_indices
    call calculate_node_indices(node_indices)

    n_tor_local = a_mat%i_tor_max - a_mat%i_tor_min + 1
    zbig = 1.d12
    zbig_backup = zbig

    do i=1, n_local_elms

       ielm = local_elms(i)

       do iv=1, n_vertex_max

          inode = element_list%element(ielm)%vertex(iv)

          if (node_list%node(inode)%boundary .ne. 0) then

             do in=a_mat%i_tor_min, a_mat%i_tor_max 
               if (keep_n0_const  .and.  in .eq. 1 ) then
                 zbig = 1.d15
               else
                 zbig = zbig_backup
               endif

                do k=1, n_var

                   !------------------------------------ the open field lines (in case of x-point grid)
                   if ((node_list%node(inode)%boundary .eq. 1) .or. (node_list%node(inode)%boundary .eq. 3)) then

                      if ((k .eq. 1) .or. (k .eq. 2) .or. (k .eq. 3) .or. &
                           (k .eq. 4) .or. (k .eq. 5) .or. (k .eq. 6) ) then
 
                       if ( (.not. is_freebound(in,k)) ) then ! apply fixed boundary conditions where necessary

                         ! --- Fix derivatives in one direction
                         iv_dir = 2
                         do kk = 1,(n_order+1)/2
                           if ( (iv_dir .eq. 3) .and. (kk .gt. 1) ) cycle ! do only t-derivatives and node value
                           do ll = 1,(n_order+1)/2
                             if ( (iv_dir .eq. 2) .and. (ll .gt. 1) ) cycle ! do only s-derivatives and node value
                             index_tmp = node_indices(kk,ll)
                             index_node = node_list%node(inode)%index(index_tmp)
                             call boundary_conditions_add_one_entry(                 &
                                    index_node, k, in, index_node, k, in,            &
                                    zbig, index_min, index_max, a_mat)
                           enddo
                         enddo

                       endif
                     endif
                   endif

                   !------------------------------------ wall aligned with fluxsurface (in case of x-point grid)
                   if ((node_list%node(inode)%boundary .eq. 2) .or. (node_list%node(inode)%boundary .eq. 3)) then

                      if ( (.not. is_freebound(in,k)) ) then ! apply fixed boundary conditions where necessary

                         ! --- Fix derivatives in one direction
                         iv_dir = 3
                         do kk = 1,(n_order+1)/2
                           if ( (iv_dir .eq. 3) .and. (kk .gt. 1) ) cycle ! do only t-derivatives and node value
                           do ll = 1,(n_order+1)/2
                             if ( (iv_dir .eq. 2) .and. (ll .gt. 1) ) cycle ! do only s-derivatives and node value
                             index_tmp = node_indices(kk,ll)
                             index_node = node_list%node(inode)%index(index_tmp)
                             call boundary_conditions_add_one_entry(                 &
                                    index_node, k, in, index_node, k, in,            &
                                    zbig, index_min, index_max, a_mat)
                           enddo
                         enddo

                      endif

                   endif

                enddo

             enddo
          endif
       enddo
    enddo

    return
  end subroutine boundary_conditions
end module mod_boundary_conditions
