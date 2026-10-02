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
!*   gmres        - boolean indicating if we are using GMRES method            *
!*   solve_only   - Indicate if we want to perform only solve                  *
!*                                                                             *
!*******************************************************************************
module mod_boundary_conditions
implicit none
contains
  subroutine boundary_conditions( my_id, node_list, element_list, bnd_node_list, local_elms,& 
                                n_local_elms, index_min, index_max, rhs_loc, xpoint2,     &
                                xcase2, R_axis, Z_axis, psi_axis, psi_bnd,                &
                                R_xpoint, Z_xpoint, psi_xpoint, a_mat)

    use mod_assembly, only : boundary_conditions_add_one_entry, boundary_conditions_add_RHS

    use phys_module, only: F0, bc_natural_open, vpar_sbc_enable,&
                           particle_flux_sbc_enable, heat_flux_sbc_enable, &
                           loop_voltage, tstep, central_density, central_mass
    use mod_boundary_ndotB, only: get_vpar_target_for_column, get_vpar_target_dT_for_column
    use mod_model_settings, only: var_Psi, var_Phi, var_zj, var_w, var_rho, var_T, &
                                  var_Vpar, var_Ti, var_Te, n_var
    use vacuum, only: is_freebound
    use constants, only: MU_ZERO, ATOMIC_MASS_UNIT
    use mpi_mod
    use mod_locate_irn_jcn
    use mod_integer_types
    use data_structure

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
    logical,                            intent(in)    :: xpoint2
    integer,                            intent(in)    :: xcase2
    real*8,                             intent(in)    :: R_axis
    real*8,                             intent(in)    :: Z_axis
    real*8,                             intent(in)    :: psi_axis
    real*8,                             intent(in)    :: psi_bnd
    real*8,                             intent(in)    :: R_xpoint(2)
    real*8,                             intent(in)    :: Z_xpoint(2)
    real*8,                             intent(in)    :: psi_xpoint(2)
    real*8,                             intent(inout) :: rhs_loc(*)
    type(type_SP_MATRIX)                              :: a_mat

    ! Internal parameters
    real*8                :: zbig
    integer               :: i, in, iv, inode, k
    integer               :: ielm
    integer               :: index_node
    
    ! v_par SBC variables
    real*8                :: vpar_target, vpar_current, delta_vpar
    real*8                :: ndotB_norm, T_local
    real*8                :: cs, alpha_rad, alpha0_rad, factor_sbc

    zbig = 1.d12
       do i=1, n_local_elms

          ielm = local_elms(i)

          do iv=1, n_vertex_max

             inode = element_list%element(ielm)%vertex(iv)

             if (node_list%node(inode)%boundary .ne. 0) then

                do in=a_mat%i_tor_min, a_mat%i_tor_max

                   do k=1, n_var
                      ! Apply loop voltage to drive psi evolution (n=0 mode at boundary)
                      ! Physics: loop_voltage drives Ohmic current via Faraday's law:
                      !   d(psi)/dt = -V_loop => psi(t) = psi(0) - V_loop * t
                      !  Loop voltage: n=0 mode, psi only - matches model600 pattern
                      if (loop_voltage .ne. 0.d0 .and. k .eq. var_Psi .and. in .eq. 1) then
                          if ( (.not. is_freebound(in, var_psi)) ) then
                              index_node = node_list%node(inode)%index(1)
                              call boundary_conditions_add_RHS(       &
                                        index_node, var_psi, in, index_min, index_max,         &
                                        RHS_loc, zbig*loop_voltage*sqrt(MU_ZERO*central_density*central_mass*ATOMIC_MASS_UNIT*1.d20)*tstep, &
                                        a_mat%i_tor_min, a_mat%i_tor_max)
                          endif
                      endif

                     ! Skip Dirichlet BC for density when particle flux SBC is enabled
                     ! (weak-form natural BC handled in mod_boundary_matrix_open.f90)
                     if (particle_flux_sbc_enable .and. k .eq. var_rho) cycle
                     ! Skip Dirichlet BC for temperature when heat flux SBC is enabled
                     if (heat_flux_sbc_enable .and. k .eq. var_T) cycle

                      !------------------------------------ boundary nodes (types 1, 2, 3)
                      ! Type 1: Open field lines (tokamak divertor targets)
                      ! Type 2: Wall/limiter (stellarator with divertor geometry)
                      ! Type 3: Corner nodes (both target and wall)
                      if ((node_list%node(inode)%boundary .eq. 1) .or. &
                          (node_list%node(inode)%boundary .eq. 2) .or. &
                          (node_list%node(inode)%boundary .eq. 3)) then

                          if ( (.not. is_freebound(in,k)) ) then ! apply fixed boundary conditions where necessary

                            ! Value DOF: constrain for ALL boundary types (1,2,3) — v_par target lives here
                            index_node = node_list%node(inode)%index(1)        
                            call boundary_conditions_add_one_entry(                 &
                                   index_node, k, in, index_node, k, in,            &
                                   zbig, index_min, index_max, a_mat)
                            
                            ! v_par Dirichlet BC: add RHS term for non-zero target
                            ! For n=0 mode, use toroidally-averaged ndotB
                            ! For n>0 modes, use Fourier coefficients of ndotB
                            if (k .eq. var_Vpar) then
                              ! Compute angle-dependent target if SBC enabled
                              if (vpar_sbc_enable) then
                                ! calculate delta vpar and apply, just like in model 600
                                delta_vpar = get_vpar_target_for_column(inode, in) &
                                              - node_list%node(inode)%values(in,1,var_Vpar)

                                call boundary_conditions_add_RHS(                     &
                                        index_node, k, in, index_min, index_max,       &
                                        rhs_loc, zbig * delta_vpar,                    &
                                        a_mat%i_tor_min, a_mat%i_tor_max)

                                 ! add dT derivative entry if using local T
                                 if (var_T .gt. 0) then
                                    call boundary_conditions_add_one_entry(                                &
                                          index_node, var_Vpar, in, index_node, var_T, in,                &
                                          -zbig * get_vpar_target_dT_for_column(inode, in),               &
                                          index_min, index_max, a_mat)
                                 endif
                              endif
                            endif

                            ! index(2) derivative for types 1 and 3 only, not type 2
                            if ((node_list%node(inode)%boundary .eq. 1) .or. (node_list%node(inode)%boundary .eq. 3)) then
                              index_node = node_list%node(inode)%index(2)
                              call boundary_conditions_add_one_entry(                 &
                                    index_node, k, in, index_node, k, in,            &
                                    zbig, index_min, index_max, a_mat)
                            endif
                            
                          endif
                      endif

                      !------------------------------------ wall aligned with fluxsurface (in case of x-point grid)
                      if ((node_list%node(inode)%boundary .eq. 2) .or. (node_list%node(inode)%boundary .eq. 3)) then

                         if ( (.not. is_freebound(in,k)) ) then ! apply fixed boundary conditions where necessary
                         ! free-boundary not yet implemented for 183

                            ! --- constrain second tangential derivative DOF at type 2/3 nodes (X-point corner geometry)
                            index_node = node_list%node(inode)%index(3)

                            call boundary_conditions_add_one_entry(                 &
                                   index_node, k, in, index_node, k, in,            &
                                   zbig, index_min, index_max, a_mat)

                         endif
                      endif
                   enddo  ! k=1,n_var (variables loop)
                enddo  ! in (toroidal modes)
             endif
          enddo
       enddo

    return
  end subroutine boundary_conditions
end module mod_boundary_conditions
