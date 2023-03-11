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
contains
  subroutine boundary_conditions( my_id, node_list, element_list, bnd_node_list, local_elms,          &
                                  n_local_elms, index_min, index_max, rhs_loc, xpoint2, xcase2,       & 
                                  R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint, Z_xpoint, psi_xpoint, a_mat)

    use mod_assembly, only : boundary_conditions_add_one_entry, boundary_conditions_add_RHS
    use data_structure
    use phys_module, only: F0, GAMMA, freeboundary, keep_n0_const
    use mpi_mod
    use mod_locate_irn_jcn
    use mod_integer_types

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
    real*8                :: zbig, zbig_backup,  T0, Vpar0, bigR, dT0_ds, dVpar0_ds, dBigR_ds
    real*8                :: R_s, R_t, Z_s, Z_t, ps0_s, ps0_t, ps0_x, ps0_y, direction, xjac
    real*8                :: Btot
    real*8                :: grad_psi, u0_s, u0_t, u0_x, u0_y
    integer               :: i, in, iv, inode, k
    integer               :: ielm
    integer               :: index_node, index_node2
    integer(kind=int_all) :: ijA_position,ijA_position2
    integer               :: ilarge2, kv, kT, ku, ilarge_vv, ilarge_vT, ilarge_vus
    integer               :: ilarge_vsvs, ilarge_vsTs, ilarge_vsT
    integer               :: first_tor, last_tor, ierr, n_tor_local

    n_tor_local = a_mat%i_tor_max - a_mat%i_tor_min + 1 
    zbig = 1.d10
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

                         if ( (k .eq. 1) .or. (k .eq. 2) ) then

                            index_node = node_list%node(inode)%index(1)
                            if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then

                               call locate_irn_jcn(index_node,index_node,index_min,index_max,ijA_position,a_mat)

                               ilarge2 = ijA_position - 1 + ((k-1)*n_tor_local + in-a_mat%i_tor_min) * n_var*n_tor_local & 
                                 + (k-1)*n_tor_local +  in - a_mat%i_tor_min +1

                               a_mat%irn(ilarge2) = n_tor_local * n_var * (index_node-1) + (k-1)*n_tor_local + in - a_mat%i_tor_min +1
                               a_mat%jcn(ilarge2) = n_tor_local * n_var * (index_node-1) + (k-1)*n_tor_local + in - a_mat%i_tor_min +1
                               a_mat%val(ilarge2) = zbig

                            endif

                            index_node = node_list%node(inode)%index(2)
                            if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then

                               call locate_irn_jcn(index_node,index_node,index_min,index_max,ijA_position,a_mat)

                               ilarge2 = ijA_position - 1 + ((k-1)*n_tor_local + in-a_mat%i_tor_min) * n_var*n_tor_local & 
                                 + (k-1)*n_tor_local + in - a_mat%i_tor_min +1

                               a_mat%irn(ilarge2) = n_tor_local * n_var * (index_node-1) + (k-1)*n_tor_local + in - a_mat%i_tor_min +1
                               a_mat%jcn(ilarge2) = n_tor_local * n_var * (index_node-1) + (k-1)*n_tor_local + in - a_mat%i_tor_min +1
                               a_mat%val(ilarge2) = zbig

                            endif
                         endif

                      end if

                      !------------------------------------ wall aligned with fluxsurface (in case of x-point grid)

                      if ((node_list%node(inode)%boundary .eq. 2) .or. (node_list%node(inode)%boundary .eq. 3)) then

                         if ( (k .eq. 1) .or. (k .eq. 2) ) then

                            index_node = node_list%node(inode)%index(1)
                            if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then

                               call locate_irn_jcn(index_node,index_node,index_min,index_max,ijA_position,a_mat)

                               ilarge2 = ijA_position - 1 + ((k-1)*n_tor_local + in-a_mat%i_tor_min) * n_var*n_tor_local & 
                                 + (k-1)*n_tor_local + in - a_mat%i_tor_min +1

                               a_mat%irn(ilarge2) = n_tor_local * n_var * (index_node-1) + (k-1)*n_tor_local + in - a_mat%i_tor_min +1
                               a_mat%jcn(ilarge2) = n_tor_local * n_var * (index_node-1) + (k-1)*n_tor_local + in - a_mat%i_tor_min +1
                               a_mat%val(ilarge2) = zbig

                            endif
                            index_node = node_list%node(inode)%index(3)

                            if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then

                               call locate_irn_jcn(index_node,index_node,index_min,index_max,ijA_position,a_mat)

                               ilarge2 = ijA_position - 1 + ((k-1)*n_tor_local + in-a_mat%i_tor_min) * n_var*n_tor_local & 
                                 + (k-1)*n_tor_local + in - a_mat%i_tor_min +1

                               a_mat%irn(ilarge2) = n_tor_local * n_var * (index_node-1) + (k-1)*n_tor_local + in - a_mat%i_tor_min +1
                               a_mat%jcn(ilarge2) = n_tor_local * n_var * (index_node-1) + (k-1)*n_tor_local + in - a_mat%i_tor_min +1
                               a_mat%val(ilarge2) = zbig

                            endif
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
