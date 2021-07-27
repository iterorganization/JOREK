module mod_boundary_conditions
contains
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
  !*   index_max    - Maximal index of local elements (                           
  !*   xpoint2      -                                                            *
  !*   xcase2       -                                                            *
  !*   psi_axis     -                                                            *
  !*   psi_bnd      -                                                            *
  !*   Z_xpoint     -                                                            *
  !*   gmres        - boolean indicating if we are using GMRES method            *
  !*   solve_only   - Indicate if we want to perform only solve                  *
  !*                                                                             *
  !*******************************************************************************
  subroutine boundary_conditions( my_id, node_list, element_list, bnd_node_list, local_elms,           &
                                  n_local_elms, index_min, index_max, rhs_loc, xpoint2, xcase2,        & 
                                  R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint, Z_xpoint, psi_xpoint,   &
                                  gmres, solve_only, ijA_index, ijA_size, irn_jcn, irn, jcn, &
                                  A_mat, i_tor_min, i_tor_max )

    use constants
    use data_structure
    use phys_module, only: F0, GAMMA, n_pol, n_tht, freeboundary, RMP_on, psi_RMP_cos, dpsi_RMP_cos_dR, dpsi_RMP_cos_dZ, &
         psi_RMP_sin, dpsi_RMP_sin_dR, dpsi_RMP_sin_dZ, t_now, RMP_growth_rate, RMP_ramp_up_time,  &
         RMP_start_time, tstep, RMP_har_cos, RMP_har_sin, grid_to_wall, n_wall_blocks, keep_n0_const
    USE tr_module
    use mpi_mod
    use mod_locate_irn_jcn
    use mod_integer_types

    implicit none

!!!! WARNING: gmres already defined in phys_module!!! Hence we use phys_module, ONLY...


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
    logical,                            intent(in)    :: gmres
    logical,                            intent(in)    :: solve_only
    integer,                            intent(in)    :: i_tor_min, i_tor_max 
    integer(kind=int_all), allocatable, intent(in)    :: ijA_index(:,:), ijA_size(:), irn_jcn(:,:) 
    integer(kind=int_all), allocatable, intent(inout) :: irn(:), jcn(:) 
    real*8,                allocatable, intent(inout) :: A_mat(:) 

    ! Internal parameters
    real*8                :: zbig, zbig_backup, T0, Vpar0, bigR, dT0_ds, dVpar0_ds, dBigR_ds, psi_1, R_1, Z_1
    real*8                :: R_s, R_t, Z_s, Z_t, ps0_s, ps0_t, ps0_x, ps0_y, direction, xjac
    real*8                :: Btot, alpha, dT0_dt, dVpar0_dt, dBigR_dt, R_inside, Z_inside
    real*8                :: grad_psi, u0_s, u0_t, u0_x, u0_y
    integer               :: i, in, iv, inode, k, ii
    integer               :: ielm
    integer               :: index_node, index_node2
    integer(kind=int_all) :: ijA_position,ijA_position2
    integer               :: ilarge2, kv, kT, ku, ilarge_vv, ilarge_vT, ilarge_vus
    integer               :: ilarge_vsvs, ilarge_vsTs, ilarge_vsT, ilarge_vut, ilarge_vtvt, ilarge_vtTt, ilarge_vtT
    integer               :: ierr, n_tor_local

    !=============== RMP ==============
    real*8, allocatable   :: psi_RMP_cos1(:),dpsi_RMP_cos_dR1(:),dpsi_RMP_cos_dZ1(:)
    real*8, allocatable   :: psi_RMP_sin1(:),dpsi_RMP_sin_dR1(:),dpsi_RMP_sin_dZ1(:)
    real*8                :: Rnode, dRnode_ds, Znode, dZnode_ds, dRnode_dt, dZnode_dt, establish_RMP
    real*8                :: delta_psi_rmp, delta_psi_rmp_dR, delta_psi_rmp_dZ, delta_psi_rmp_ds, delta_psi_rmp_dt, psi_test, sigmo_fonc
    integer               :: ilarge_vp, ilarge_vp2
    integer               :: kp, j, err, itest

    n_tor_local = i_tor_max - i_tor_min + 1
    if (RMP_on .and. (n_tor .ge. 3)) then

       call tr_allocate(psi_RMP_cos1,1, bnd_node_list%n_bnd_nodes,"psi_RMP_cos1",CAT_UNKNOWN)
       call tr_allocate(dpsi_RMP_cos_dR1,1,bnd_node_list%n_bnd_nodes,"dpsi_RMP_cos_dR1",CAT_UNKNOWN)
       call tr_allocate(dpsi_RMP_cos_dZ1,1,bnd_node_list%n_bnd_nodes,"dpsi_RMP_cos_dZ1",CAT_UNKNOWN)
       call tr_allocate(psi_RMP_sin1,1,bnd_node_list%n_bnd_nodes,"psi_RMP_sin1",CAT_UNKNOWN)
       call tr_allocate(dpsi_RMP_sin_dR1,1,bnd_node_list%n_bnd_nodes,"dpsi_RMP_sin_dR1",CAT_UNKNOWN)
       call tr_allocate(dpsi_RMP_sin_dZ1,1,bnd_node_list%n_bnd_nodes,"dpsi_RMP_sin_dZ1",CAT_UNKNOWN)

    psi_test =  node_list%node(bnd_node_list%bnd_node(1)%index_jorek)%values(RMP_har_cos,1,1)
    ! if necessary, replace by:
    ! psi_test =  node_list%node(bnd_node_list%bnd_node(1)%index_jorek)%values(min(RMP_har_cos, n_tor),1,1)
    write (*,*) 'psi_bnd at previous time step', psi_test
    
    if (abs(psi_test) .le. abs(psi_RMP_cos(1))) then
      sigmo_fonc = ( 1.d0 + exp(-RMP_growth_rate*( t_now - RMP_start_time - RMP_ramp_up_time/2.d0 )))**(-1) &
          - ( 1.d0 + exp(-RMP_growth_rate*( 0.d0 - RMP_ramp_up_time/2.d0 )))**(-1) 
      establish_RMP = (RMP_growth_rate*sigmo_fonc*(1-sigmo_fonc)+1.e-6)*tstep 
    else
      establish_RMP = 0.d0
    endif
    ! Other possibility (simpler) : if ( (t_now - RMP_start_time) .ge. 2.2*RMP_ramp_up_time/2.d0 ) then establish_RMP =0.0
  
    do j=1, bnd_node_list%n_bnd_nodes  
      psi_RMP_cos1(j)     = psi_RMP_cos(j)     * establish_RMP
      dpsi_RMP_cos_dR1(j) = dpsi_RMP_cos_dR(j) * establish_RMP
      dpsi_RMP_cos_dZ1(j) = dpsi_RMP_cos_dZ(j) * establish_RMP
      psi_RMP_sin1(j)     = psi_RMP_sin(j)     * establish_RMP
      dpsi_RMP_sin_dR1(j) = dpsi_RMP_sin_dR(j) * establish_RMP
      dpsi_RMP_sin_dZ1(j) = dpsi_RMP_sin_dZ(j) * establish_RMP
    end do

       if (my_id == 0) then

          write (*,*) 'psi_RMP_cos1(1) and derivatives after multiplication in boundary conditions'
          write (*,*) psi_RMP_cos1(1), dpsi_RMP_cos_dR1(1), dpsi_RMP_cos_dZ1(1)
          write (*,*) 'establish_RMP', establish_RMP

       endif

    endif
    !=============== RMP ==============

    zbig = 1.d12
    zbig_backup = zbig
       do i=1, n_local_elms

          ielm = local_elms(i)

          do iv=1, n_vertex_max

             inode = element_list%element(ielm)%vertex(iv)


             ! A crude way of imposing partial regularity at the grid axis
             ! May be improved upon in the future ( see forthcoming paper on full
             ! MHD JOREK )
             !---------------------------------------------------------------------------------------------
             do in=i_tor_min, i_tor_max 
               if (keep_n0_const  .and.  in .eq. 1 ) then
                 zbig = 1.d15
               else
                 zbig = zbig_backup
               endif
                do k=1, n_var
                   ! Restrain the coefficients of the 2nd, 3rd, 4th basis functions
                   ! on axis from changing
                   if ( ( inode <= n_tht .or. ( n_tht < 1 .and. inode <= n_pol ) ).and. 1==1 ) then
                      do ii=2,4
                         index_node = node_list%node(inode)%index(ii)

                            if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then

                               call locate_irn_jcn(index_node,index_node,index_min,index_max,ijA_position,ijA_index, ijA_size, irn_jcn)

                               ilarge2 = ijA_position - 1 + ((k-1)*n_tor_local + in-i_tor_min) *n_var*n_tor_local &
                                                          +  (k-1)*n_tor_local + in - i_tor_min +1

                               irn(ilarge2) =  n_tor_local * n_var * (index_node-1) +(k-1)*n_tor_local + in - i_tor_min +1
                               jcn(ilarge2) =  n_tor_local * n_var * (index_node-1) +(k-1)*n_tor_local + in - i_tor_min +1
                               A_mat(ilarge2)    = zbig

                            end if
                      enddo ! do loop over indices 2-4

                   endif

                enddo
             enddo
             !---------------------------------------------------------------------------------------------


             if (node_list%node(inode)%boundary .ne. 0) then

                do in=i_tor_min, i_tor_max 
               if (keep_n0_const  .and.  in .eq. 1 ) then
                 zbig = 1.d15
               else
                 zbig = zbig_backup
               endif

                   do k=1, n_var

                      !========================================================================
                      ! conditions for direction 1 (s), i.e. boundary types 1, 3, 4, 9
                      ! apply fixed bc for variables k=1,2,3,4
                      ! apply v_par = cs for k=7
                      !========================================================================
                      if     ((node_list%node(inode)%boundary .eq.  1) &
                         .or. (node_list%node(inode)%boundary .eq. 11) &
                         .or. (node_list%node(inode)%boundary .eq.  9) &
                         .or. (node_list%node(inode)%boundary .eq. 19) &
                         .or. (node_list%node(inode)%boundary .eq.  3) &
                         .or. (node_list%node(inode)%boundary .eq.  4)) then

                         !====================================== beginning RMPs at boundary ======================================================
                         !================================== type 1 - boundary: only depends on 's'
                         ! ======================================================================================================================

                         if (RMP_on ) then

                            if ((k.eq.1) .and. ((in.eq.2) .or. (in.eq.3)) .and. (.not. freeboundary)) then
                               ! in .eq. 2 corresponds to cos(n_perturbation)
                               ! in .eq. 3 corresponds to sin(n_perturbation)

                               kp=1    ! variable psi
                               kv=1    ! equation for psi

                               index_node = node_list%node(inode)%index(1)  ! index in RHS (or matrix A not compressed)

                               Rnode     = node_list%node(inode)%x(1,1,1) 
                               dRnode_ds = node_list%node(inode)%x(1,2,1) 
                               Znode     = node_list%node(inode)%x(1,1,2) 
                               dZnode_ds = node_list%node(inode)%x(1,2,2) 

                               if (in.eq.2) then
                                  delta_psi_rmp = psi_RMP_cos1(node_list%node(inode)%boundary_index)
                                  delta_psi_rmp_dR = dpsi_RMP_cos_dR1(node_list%node(inode)%boundary_index)
                                  delta_psi_rmp_dZ = dpsi_RMP_cos_dZ1(node_list%node(inode)%boundary_index)

                               else 
                                  delta_psi_rmp = psi_RMP_sin1(node_list%node(inode)%boundary_index)
                                  delta_psi_rmp_dR = dpsi_RMP_sin_dR1(node_list%node(inode)%boundary_index)
                                  delta_psi_rmp_dZ = dpsi_RMP_sin_dZ1(node_list%node(inode)%boundary_index)

                               endif

                               delta_psi_rmp_ds = delta_psi_rmp_dR * dRnode_ds + delta_psi_rmp_dZ * dZnode_ds

                               if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then

                                  call locate_irn_jcn(index_node,index_node,index_min,index_max,ijA_position,ijA_index, ijA_size, irn_jcn)

                                  !-------- index dans A_mat
                                  ilarge_vp  = ijA_position  - 1 + ((kv-1)*n_tor_local + in-i_tor_min) * n_var*n_tor_local &
                                                                 +  (kp-1)*n_tor_local + in - i_tor_min +1

                                  Rhs_loc(n_tor_local*n_var * (index_node-1) + (kv-1)*n_tor_local + in - i_tor_min +1) = ZBIG * delta_psi_rmp

                                  irn(ilarge_vp) =  n_tor_local * n_var * (index_node-1) + (kv-1)*n_tor_local + in - i_tor_min +1
                                  jcn(ilarge_vp) =  n_tor_local * n_var * (index_node-1) + (kp-1)*n_tor_local + in - i_tor_min +1
                                  A_mat(ilarge_vp)   = ZBIG
                               endif

                               index_node2 = node_list%node(inode)%index(2)

                               if ((index_node2 .ge. index_min) .and. (index_node2 .le. index_max)) then                         
                                  call locate_irn_jcn(index_node2,index_node2,index_min,index_max,ijA_position2,ijA_index, ijA_size, irn_jcn)

                                  ilarge_vp2  = ijA_position2  - 1 + ((kv-1)*n_tor_local + in-i_tor_min) * n_var*n_tor_local &
                                                                   +  (kp-1)*n_tor_local + in - i_tor_min +1

                                  Rhs_loc(n_tor_local*n_var * (index_node2-1) + (kv-1)*n_tor_local + in - i_tor_min +1) = ZBIG * delta_psi_rmp_ds

                                  irn(ilarge_vp2) =  n_tor_local * n_var * (index_node2-1) + (kv-1)*n_tor_local + in - i_tor_min +1
                                  jcn(ilarge_vp2) =  n_tor_local * n_var * (index_node2-1) + (kp-1)*n_tor_local + in - i_tor_min +1
                                  A_mat(ilarge_vp2)   = ZBIG
                               endif
                            endif


                         endif !(end RMP)
                         !======================================= end RMPs ==================================


                         if (                  &
                              ((k .eq. 1) .and. (.not. RMP_on) .and. ( in .ge. 2 ))      &
                              .or. ((k .eq. 1) .and. (in .eq. 1))  & 
                              .or. ((k .eq. 1) .and. (in .ge. 4))  & 
                              .or. (k .eq. 2)  &
                              .or. (k .eq. 3)  &
                              .or. (k .eq. 4)  &
                                !.or. (k .eq. 5) &
                                !.or. (k .eq. 6) &
                                !.or. (k .eq. 7) &
                              ) then


                            index_node = node_list%node(inode)%index(1)
                               if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then

                                  call locate_irn_jcn(index_node,index_node,index_min,index_max,ijA_position,ijA_index, ijA_size, irn_jcn)

                                  ilarge2 = ijA_position - 1 + ((k-1)*n_tor_local + in-i_tor_min) * n_var*n_tor_local & 
                                                             +  (k-1)*n_tor_local + in - i_tor_min +1

                                  irn(ilarge2) =  n_tor_local * n_var * (index_node-1) + (k-1)*n_tor_local + in - i_tor_min +1
                                  jcn(ilarge2) =  n_tor_local * n_var * (index_node-1) + (k-1)*n_tor_local + in - i_tor_min +1
                                  A_mat(ilarge2)   = zbig

                               endif

                            index_node = node_list%node(inode)%index(2)
                               if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then

                                  call locate_irn_jcn(index_node,index_node,index_min,index_max,ijA_position,ijA_index, ijA_size, irn_jcn)

                                  ilarge2 = ijA_position - 1 + ((k-1)*n_tor_local + in-i_tor_min) * n_var*n_tor_local & 
                                                             +  (k-1)*n_tor_local + in - i_tor_min +1

                                  irn(ilarge2) =  n_tor_local * n_var * (index_node-1) + (k-1)*n_tor_local + in - i_tor_min +1
                                  jcn(ilarge2) =  n_tor_local * n_var * (index_node-1) + (k-1)*n_tor_local + in - i_tor_min +1
                                  A_mat(ilarge2)    = zbig

                               endif
                         endif


                         if (k .eq. 7) then

                            index_node  = node_list%node(inode)%index(1)             ! position of value
                            index_node2 = node_list%node(inode)%index(2)             ! position of first deriative

                            T0        = abs(node_list%node(inode)%values(1,1,6))
                            Vpar0     = node_list%node(inode)%values(1,1,7)
                            BigR      = node_list%node(inode)%x(1,1,1)
                            dT0_ds    = node_list%node(inode)%values(1,2,6)
                            dVpar0_ds = node_list%node(inode)%values(1,2,7)
                            dBigR_ds  = node_list%node(inode)%x(1,2,1)

                            ps0_s     = node_list%node(inode)%values(1,2,1)
                            ps0_t     = node_list%node(inode)%values(1,3,1)

                            U0_s      = node_list%node(inode)%values(1,2,2)
                            U0_t      = node_list%node(inode)%values(1,3,2)

                            R_s       = node_list%node(inode)%x(1,2,1)
                            R_t       = node_list%node(inode)%x(1,3,1)
                            Z_s       = node_list%node(inode)%x(1,2,2)
                            Z_t       = node_list%node(inode)%x(1,3,2)

                            xjac  =  R_s*Z_t - R_t*Z_s
                            ps0_x = (   Z_t * ps0_s - Z_s * ps0_t ) / xjac
                            ps0_y = ( - R_t * ps0_s + R_s * ps0_t ) / xjac

                            u0_x = (   Z_t * u0_s - Z_s * u0_t ) / xjac
                            u0_y = ( - R_t * u0_s + R_s * u0_t ) / xjac

                            alpha = ((Z_axis) - Z_xpoint(1))/(R_axis - R_xpoint(1))

                            R_inside = alpha*(node_list%node(inode)%x(1,1,2)-Z_xpoint(1)) &
                                 + node_list%node(inode)%x(1,1,1) + alpha**2 * R_xpoint(1)
                            R_inside = R_inside / (1.d0 + alpha**2)
                            Z_inside = alpha * (R_inside - R_xpoint(1)) + Z_xpoint(1)

                            R_inside = min(max(R_inside,R_xpoint(1)),R_axis)
                            Z_inside = min(max(Z_inside,Z_xpoint(1)),Z_axis)

                            direction = ps0_s * (  (node_list%node(inode)%x(1,1,1)-R_inside)*Z_s &
                                 - (node_list%node(inode)%x(1,1,2)-Z_inside)*R_s)
                            direction = direction / abs(direction)

                            if (xcase2 .eq. UPPER_XPOINT) then
                               direction = -direction
                            else if ((xcase2 .eq. DOUBLE_NULL).and.(node_list%node(inode)%x(1,1,2).gt.Z_axis +0.1).and.(node_list%node(inode)%x(1,1,1).gt.R_xpoint(2))) then
                              direction = -1.
                            else if ((xcase2 .eq. DOUBLE_NULL) .and. (node_list%node(inode)%x(1,1,2).gt.Z_axis +0.1).and.(node_list%node(inode)%x(1,1,1).lt.R_xpoint(2)))then
                              direction = +1.
                            end if

                            ! --- Special field direction for grid patches
                            ! --- (hopefully temporary until Vpar is treated in boundary_matrix_open)
                            if ( (grid_to_wall) .and. (n_wall_blocks .gt. 0) ) then
                              direction = 1
                              if (node_list%node(inode)%boundary .eq. 11) direction = -direction
                              if (node_list%node(inode)%boundary .eq. 19) direction = -direction
                            endif
                            
                            grad_psi = sqrt(ps0_x**2 + ps0_y**2)

                            Btot = sqrt(F0**2 + ps0_x**2 + ps0_y**2) / BigR

                            if (in .eq. 1) then

                               !                write(*,'(A,3e14.6,A,e14.6)') ' Boundary : ',Vpar0, -BigR**2 * u0_s/ps0_s, direction*sqrt(GAMMA*T0)/Btot,&
                               !                                              ' error : ',Vpar0 - BigR**2 * u0_s/ps0_s - direction*sqrt(GAMMA*T0)/Btot

                            endif

                            ku = 2
                            kv = 7
                            kT = 6

                               if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then

                                  call locate_irn_jcn(index_node,index_node, index_min,index_max,ijA_position,ijA_index, ijA_size, irn_jcn)
                                  call locate_irn_jcn(index_node,index_node2,index_min,index_max,ijA_position2,ijA_index, ijA_size, irn_jcn)

                                  ilarge_vv  = ijA_position  - 1 + ((kv-1)*n_tor_local + in-i_tor_min) * n_var*n_tor_local &
                                                                 +  (kv-1)*n_tor_local + in - i_tor_min +1
                                  ilarge_vT  = ijA_position  - 1 + ((kv-1)*n_tor_local + in-i_tor_min) * n_var*n_tor_local &
                                                                 +  (kT-1)*n_tor_local + in - i_tor_min +1
                                  ilarge_vus = ijA_position2 - 1 + ((kv-1)*n_tor_local + in-i_tor_min) * n_var*n_tor_local &
                                                                 +  (ku-1)*n_tor_local + in - i_tor_min +1

                                  irn(ilarge_vv) =  n_tor_local * n_var * (index_node-1) + (kv-1)*n_tor_local + in - i_tor_min +1
                                  jcn(ilarge_vv) =  n_tor_local * n_var * (index_node-1) + (kv-1)*n_tor_local + in - i_tor_min +1
                                  A_mat(ilarge_vv)   =  zbig

                                  irn(ilarge_vT) =  n_tor_local * n_var * (index_node-1) + (kv-1)*n_tor_local + in - i_tor_min +1
                                  jcn(ilarge_vT) =  n_tor_local * n_var * (index_node-1) + (kT-1)*n_tor_local + in - i_tor_min +1
                                  A_mat(ilarge_vT)   = - zbig / Btot * 0.5d0 * GAMMA / sqrt(GAMMA*T0) * direction

                                  irn(ilarge_vus) =  n_tor_local * n_var * (index_node -1) + (kv-1)*n_tor_local + in - i_tor_min +1
                                  jcn(ilarge_vus) =  n_tor_local * n_var * (index_node2-1) + (ku-1)*n_tor_local + in - i_tor_min +1
                                  A_mat(ilarge_vus)   = - zbig * BigR**2 / ps0_s


                                  if (in .eq. 1) then
                                     RHS_loc(n_tor_local*n_var * (index_node-1) + (kv-1)*n_tor_local + in - i_tor_min +1) = &
                                          Zbig * ( - Vpar0 + BigR**2 * U0_s /ps0_s + direction*sqrt(GAMMA*T0) / Btot)
                                  else
                                     RHS_loc(n_tor_local*n_var * (index_node-1) + (kv-1)*n_tor_local + in - i_tor_min +1) = 0.d0
                                  endif

                               endif

                            index_node  = node_list%node(inode)%index(1)
                            index_node2 = node_list%node(inode)%index(2)
                            kv = 7
                            kT = 6

                               if ((index_node2 .ge. index_min) .and. (index_node2 .le. index_max)) then

                                  call locate_irn_jcn(index_node2,index_node,index_min,index_max,ijA_position,ijA_index, ijA_size, irn_jcn)
                                  call locate_irn_jcn(index_node2,index_node2,index_min,index_max,ijA_position2,ijA_index, ijA_size, irn_jcn)

                                  ilarge_vsvs = ijA_position2 - 1 + ((kv-1)*n_tor_local + in-i_tor_min) * n_var*n_tor_local &
                                                                  +  (kv-1)*n_tor_local + in - i_tor_min +1
                                  ilarge_vsTs = ijA_position2 - 1 + ((kv-1)*n_tor_local + in-i_tor_min) * n_var*n_tor_local &
                                                                  +  (kT-1)*n_tor_local + in - i_tor_min +1
                                  ilarge_vsT  = ijA_position  - 1 + ((kv-1)*n_tor_local + in-i_tor_min) * n_var*n_tor_local &
                                                                  +  (kT-1)*n_tor_local + in - i_tor_min +1

                                  irn(ilarge_vsvs) =  n_tor_local * n_var * (index_node2-1) + (kv-1)*n_tor_local + in - i_tor_min +1
                                  jcn(ilarge_vsvs) =  n_tor_local * n_var * (index_node2-1) + (kv-1)*n_tor_local + in - i_tor_min +1
                                  A_mat(ilarge_vsvs)   = zbig

                                  irn(ilarge_vsTs) =  n_tor_local * n_var * (index_node2-1) + (kv-1)*n_tor_local + in - i_tor_min +1
                                  jcn(ilarge_vsTs) =  n_tor_local * n_var * (index_node2-1) + (kT-1)*n_tor_local + in - i_tor_min +1
                                  A_mat(ilarge_vsTs)   = - zbig / Btot * 0.5d0 * GAMMA / sqrt(GAMMA*T0) * direction

                                  irn(ilarge_vsT) =  n_tor_local * n_var * (index_node2-1) + (kv-1)*n_tor_local + in - i_tor_min +1
                                  jcn(ilarge_vsT) =  n_tor_local * n_var * (index_node -1) + (kT-1)*n_tor_local + in - i_tor_min +1
                                  A_mat(ilarge_vsT)   = + zbig / Btot * 0.25d0 * GAMMA**2 / (GAMMA*T0)**(3/2) * dT0_ds * direction

                                  if (in .eq. 1) then
                                     Rhs_loc(n_tor_local * n_var * (index_node2-1) + (kv-1)*n_tor_local + in - i_tor_min +1) = &
                                          Zbig*(-dVpar0_ds +  0.5d0 / Btot * GAMMA / sqrt(GAMMA*T0) * dT0_ds * direction)
                                  else
                                     Rhs_loc(n_tor_local * n_var * (index_node2-1) + (kv-1)*n_tor_local + in - i_tor_min +1) = 0.d0
                                  endif

                               endif
                         end if

                      end if

                      !========================================================================
                      ! conditions for direction 2 (s), i.e. boundary types 5, 9
                      ! apply fixed bc for variables k=1,2,3,4
                      ! apply v_par = cs for k=7
                      !========================================================================
                      if    ((node_list%node(inode)%boundary .eq.     5) &
                           .or. (node_list%node(inode)%boundary .eq. 15) &
                           .or. (node_list%node(inode)%boundary .eq.  9) &
                           .or. (node_list%node(inode)%boundary .eq. 19)) then


                         if (                  &
                              (k .eq. 1)  &
                              .or. (k .eq. 2)  &
                              .or. (k .eq. 3)  &
                              .or. (k .eq. 4)  &
                                !.or. (k .eq. 5) &
                                !.or. (k .eq. 6) &
                                !.or. (k .eq. 7) &
                              ) then


                            index_node = node_list%node(inode)%index(1)
                               if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then

                                  call locate_irn_jcn(index_node,index_node,index_min,index_max,ijA_position,ijA_index, ijA_size, irn_jcn)

                                  ilarge2 = ijA_position - 1 + ((k-1)*n_tor_local + in-i_tor_min) * n_var*n_tor_local &
                                                          + (k-1)*n_tor_local + in - i_tor_min +1

                                  irn(ilarge2) =  n_tor_local * n_var * (index_node-1) + (k-1)*n_tor_local + in - i_tor_min +1
                                  jcn(ilarge2) =  n_tor_local * n_var * (index_node-1) + (k-1)*n_tor_local + in - i_tor_min +1
                                  A_mat(ilarge2)   = zbig

                               endif

                            index_node = node_list%node(inode)%index(3)

                               if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then

                                  call locate_irn_jcn(index_node,index_node,index_min,index_max,ijA_position,ijA_index, ijA_size, irn_jcn)

                                  ilarge2 = ijA_position - 1 + ((k-1)*n_tor_local + in-i_tor_min) * n_var*n_tor_local &
                                                             +  (k-1)*n_tor_local + in - i_tor_min +1

                                  irn(ilarge2) =  n_tor_local * n_var * (index_node-1) + (k-1)*n_tor_local + in - i_tor_min +1
                                  jcn(ilarge2) =  n_tor_local * n_var * (index_node-1) + (k-1)*n_tor_local + in - i_tor_min +1
                                  A_mat(ilarge2)    = zbig

                               endif
                         endif


                         if (k .eq. 7) then

                            index_node  = node_list%node(inode)%index(1)             ! position of value
                            index_node2 = node_list%node(inode)%index(3)             ! position of first deriative

                            T0        = abs(node_list%node(inode)%values(1,1,6));  T0 = max(T0,0.001)
                            Vpar0     = node_list%node(inode)%values(1,1,7)
                            BigR      = node_list%node(inode)%x(1,1,1)
                            dT0_dt    = node_list%node(inode)%values(1,3,6)
                            dVpar0_dt = node_list%node(inode)%values(1,3,7)
                            dBigR_dt  = node_list%node(inode)%x(1,3,1)

                            ps0_s     = node_list%node(inode)%values(1,2,1)
                            ps0_t     = node_list%node(inode)%values(1,3,1)

                            U0_s      = node_list%node(inode)%values(1,2,2)
                            U0_t      = node_list%node(inode)%values(1,3,2)

                            R_s       = node_list%node(inode)%x(1,2,1)
                            R_t       = node_list%node(inode)%x(1,3,1)
                            Z_s       = node_list%node(inode)%x(1,2,2)
                            Z_t       = node_list%node(inode)%x(1,3,2)

                            xjac  =  R_s*Z_t - R_t*Z_s
                            ps0_x = (   Z_t * ps0_s - Z_s * ps0_t ) / xjac
                            ps0_y = ( - R_t * ps0_s + R_s * ps0_t ) / xjac

                            u0_x = (   Z_t * u0_s - Z_s * u0_t ) / xjac
                            u0_y = ( - R_t * u0_s + R_s * u0_t ) / xjac

                            alpha = ((Z_axis) - Z_xpoint(1))/(R_axis - R_xpoint(1))

                            R_inside = alpha*(node_list%node(inode)%x(1,1,2)-Z_xpoint(1)) &
                                 + node_list%node(inode)%x(1,1,1) + alpha**2 * R_xpoint(1)
                            R_inside = R_inside / (1.d0 + alpha**2)
                            Z_inside = alpha * (R_inside - R_xpoint(1)) + Z_xpoint(1)

                            R_inside = min(max(R_inside,R_xpoint(1)),R_axis)
                            Z_inside = min(max(Z_inside,Z_xpoint(1)),Z_axis)

                            direction = ps0_t * (  (node_list%node(inode)%x(1,1,1)-R_inside)*Z_t &
                                 - (node_list%node(inode)%x(1,1,2)-Z_inside)*R_t)
                            direction = direction / abs(direction)

                            ! --- Special field direction for grid patches
                            ! --- (hopefully temporary until Vpar is treated in boundary_matrix_open)
                            if ( (grid_to_wall) .and. (n_wall_blocks .gt. 0) ) then
                              direction = 1
                              if (node_list%node(inode)%boundary .eq. 15) direction = -direction
                              if (node_list%node(inode)%boundary .eq. 19) direction = -direction
                            endif
                            
                            grad_psi = sqrt(ps0_x**2 + ps0_y**2)

                            Btot = sqrt(F0**2 + ps0_x**2 + ps0_y**2) / BigR

                            ku = 2
                            kv = 7
                            kT = 6
                               if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then

                                  call locate_irn_jcn(index_node,index_node, index_min,index_max,ijA_position,ijA_index, ijA_size, irn_jcn)
                                  call locate_irn_jcn(index_node,index_node2,index_min,index_max,ijA_position2,ijA_index, ijA_size, irn_jcn)

                                  ilarge_vv  = ijA_position  - 1 + ((kv-1)*n_tor_local + in-i_tor_min) * n_var*n_tor_local &
                                                                 +  (kv-1)*n_tor_local + in - i_tor_min +1
                                  ilarge_vT  = ijA_position  - 1 + ((kv-1)*n_tor_local + in-i_tor_min) * n_var*n_tor_local &
                                                                 +  (kT-1)*n_tor_local + in - i_tor_min +1
                                  ilarge_vut = ijA_position2 - 1 + ((kv-1)*n_tor_local + in-i_tor_min) * n_var*n_tor_local &
                                                                 +  (ku-1)*n_tor_local + in - i_tor_min +1

                                  irn(ilarge_vv) =  n_tor_local * n_var * (index_node-1) + (kv-1)*n_tor_local + in - i_tor_min +1
                                  jcn(ilarge_vv) =  n_tor_local * n_var * (index_node-1) + (kv-1)*n_tor_local + in - i_tor_min +1
                                  A_mat(ilarge_vv)   =  zbig

                                  irn(ilarge_vT) =  n_tor_local * n_var * (index_node-1) + (kv-1)*n_tor_local + in - i_tor_min +1
                                  jcn(ilarge_vT) =  n_tor_local * n_var * (index_node-1) + (kT-1)*n_tor_local + in - i_tor_min +1
                                  A_mat(ilarge_vT)   = - zbig / Btot * 0.5d0 * GAMMA / sqrt(GAMMA*T0) * direction

                                  irn(ilarge_vut) =  n_tor_local * n_var * (index_node -1) + (kv-1)*n_tor_local + in - i_tor_min +1
                                  jcn(ilarge_vut) =  n_tor_local * n_var * (index_node2-1) + (ku-1)*n_tor_local + in - i_tor_min +1
                                  A_mat(ilarge_vut)   = - zbig * BigR**2 / ps0_t

                                  if (in .eq. 1) then
                                     RHS_loc(n_tor_local*n_var * (index_node-1) + (kv-1)*n_tor_local + in - i_tor_min +1) = &
                                          Zbig * ( - Vpar0 + BigR**2 * U0_t /ps0_t + direction*sqrt(GAMMA*T0) / Btot)

                                  else
                                     RHS_loc(n_tor_local*n_var * (index_node-1) + (kv-1)*n_tor_local + in - i_tor_min +1) = 0.d0
                                  endif
                                  !         write(*,'(A,i6,3e16.8)') ' bc5:',inode, Vpar0,direction*sqrt(GAMMA*T0) / Btot, Vpar0-direction*sqrt(GAMMA*T0) / Btot

                               endif

                            index_node  = node_list%node(inode)%index(1)
                            index_node2 = node_list%node(inode)%index(3)
                            kv = 7
                            kT = 6
                               if ((index_node2 .ge. index_min) .and. (index_node2 .le. index_max)) then

                                  call locate_irn_jcn(index_node2,index_node,index_min,index_max,ijA_position,ijA_index, ijA_size, irn_jcn)
                                  call locate_irn_jcn(index_node2,index_node2,index_min,index_max,ijA_position2,ijA_index, ijA_size, irn_jcn)

                                  ilarge_vtvt = ijA_position2 - 1 + ((kv-1)*n_tor_local + in-i_tor_min) * n_var*n_tor_local &
                                                                  +  (kv-1)*n_tor_local + in - i_tor_min +1
                                  ilarge_vtTt = ijA_position2 - 1 + ((kv-1)*n_tor_local + in-i_tor_min) * n_var*n_tor_local &
                                                                  +  (kT-1)*n_tor_local + in - i_tor_min +1
                                  ilarge_vtT  = ijA_position  - 1 + ((kv-1)*n_tor_local + in-i_tor_min) * n_var*n_tor_local &
                                                                  +  (kT-1)*n_tor_local + in - i_tor_min +1

                                  irn(ilarge_vtvt) = n_tor_local * n_var * (index_node2-1) + (kv-1)*n_tor_local + in - i_tor_min +1
                                  jcn(ilarge_vtvt) = n_tor_local * n_var * (index_node2-1) + (kv-1)*n_tor_local + in - i_tor_min +1
                                  A_mat(ilarge_vtvt)   = zbig

                                  irn(ilarge_vtTt) =  n_tor_local * n_var * (index_node2-1) + (kv-1)*n_tor_local + in - i_tor_min +1
                                  jcn(ilarge_vtTt) =  n_tor_local * n_var * (index_node2-1) + (kT-1)*n_tor_local + in - i_tor_min +1
                                  A_mat(ilarge_vtTt)   = - zbig / Btot * 0.5d0 * GAMMA / sqrt(GAMMA*T0) * direction

                                  irn(ilarge_vtT) =  n_tor_local * n_var * (index_node2-1) + (kv-1)*n_tor_local + in - i_tor_min +1
                                  jcn(ilarge_vtT) =  n_tor_local * n_var * (index_node -1) + (kT-1)*n_tor_local + in - i_tor_min +1
                                  A_mat(ilarge_vtT)   = + zbig / Btot * 0.25d0 * GAMMA**2 / (GAMMA*T0)**(3/2) * dT0_dt * direction

                                  if (in .eq. 1) then
                                     RHS_loc(n_tor_local * n_var * (index_node2-1) + (kv-1)*n_tor_local + in - i_tor_min +1) = &
                                          Zbig*(-dVpar0_dt +  0.5d0 / Btot * GAMMA / sqrt(GAMMA*T0) * dT0_dt * direction)
                                  else
                                     RHS_loc(n_tor_local * n_var * (index_node2-1) + (kv-1)*n_tor_local + in - i_tor_min +1) = 0.d0
                                  endif

                               endif
                         end if

                      end if



                      !------------------------------------ wall aligned with fluxsurface (in case of x-point grid)
                      if    ((node_list%node(inode)%boundary .eq.  2) &
                        .or. (node_list%node(inode)%boundary .eq. 12) &
                        .or. (node_list%node(inode)%boundary .eq.  3)) then

                         !====================================== begining RMPs at boundary ======================================================
                         !================================== type 2 - boundary: only depends on 't'
                         ! ======================================================================================================================

                         if (RMP_on ) then

                            if ((k.eq.1) .and. ((in.eq.2) .or. (in.eq.3)) .and. (.not. freeboundary)) then
                               ! in .eq. 2 corresponds to cos(n_perturbation)
                               ! in .eq. 3 corresponds to sin(n_perturbation)



                               kp=1    ! variable psi
                               kv=1    ! equation for psi

                               index_node = node_list%node(inode)%index(1)  ! index in RHS (or matrix A not compressed)

                               Rnode     = node_list%node(inode)%x(1,1,1) 
                               dRnode_dt = node_list%node(inode)%x(1,3,1) 
                               Znode     = node_list%node(inode)%x(1,1,2) 
                               dZnode_dt = node_list%node(inode)%x(1,3,2) 

                               if (in.eq.2) then
                                  delta_psi_rmp = psi_RMP_cos1(node_list%node(inode)%boundary_index)
                                  delta_psi_rmp_dR = dpsi_RMP_cos_dR1(node_list%node(inode)%boundary_index)
                                  delta_psi_rmp_dZ = dpsi_RMP_cos_dZ1(node_list%node(inode)%boundary_index)

                                  if (node_list%node(inode)%boundary_index == 1 ) then
                                     write (*,*) 'type2_bnd: my_id, psi_RMP_cos1, Rnode, Znode'
                                     write (*,*) my_id, delta_psi_rmp, Rnode, Znode
                                     write (*,*) 'delta_psi_rmp_dR, delta_psi_rmp_dZ'      
                                     write (*,*) delta_psi_rmp_dR, delta_psi_rmp_dZ
                                  endif
                               else 
                                  delta_psi_rmp = psi_RMP_sin1(node_list%node(inode)%boundary_index)
                                  delta_psi_rmp_dR = dpsi_RMP_sin_dR1(node_list%node(inode)%boundary_index)
                                  delta_psi_rmp_dZ = dpsi_RMP_sin_dZ1(node_list%node(inode)%boundary_index)

                               endif

                               delta_psi_rmp_dt = delta_psi_rmp_dR * dRnode_dt + delta_psi_rmp_dZ * dZnode_dt
                               if (in.eq.2) then

                                  if (node_list%node(inode)%boundary_index == 1 ) then
                                     write (*,*) 'delta_psi_rmp_dt', delta_psi_rmp_dt
                                     write (*,*) 'delta_psi_rmp_ds', delta_psi_rmp_ds
                                  endif
                               endif

                               if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then

                                  call locate_irn_jcn(index_node,index_node,index_min,index_max,ijA_position,ijA_index, ijA_size, irn_jcn)

                                  !-------- index dans A_mat
                                  ilarge_vp  = ijA_position  - 1 + ((kv-1)*n_tor_local + in-i_tor_min) * n_var*n_tor_local &
                                                                 +  (kp-1)*n_tor_local + in - i_tor_min +1

                                  Rhs_loc(n_tor_local*n_var * (index_node-1) + (kv-1)*n_tor_local + in - i_tor_min +1) = ZBIG * delta_psi_rmp

                                  irn(ilarge_vp) =  n_tor_local * n_var * (index_node-1) + (kv-1)*n_tor_local + in - i_tor_min +1
                                  jcn(ilarge_vp) =  n_tor_local * n_var * (index_node-1) + (kp-1)*n_tor_local + in - i_tor_min +1
                                  A_mat(ilarge_vp)   = ZBIG

                               endif

                               index_node2 = node_list%node(inode)%index(3)

                               if ((index_node2 .ge. index_min) .and. (index_node2 .le. index_max)) then                         
                                  call locate_irn_jcn(index_node2,index_node2,index_min,index_max,ijA_position2,ijA_index, ijA_size, irn_jcn)

                                  ilarge_vp2  = ijA_position2  - 1 + ((kv-1)*n_tor_local + in-i_tor_min) * n_var*n_tor_local &
                                                                   +  (kp-1)*n_tor_local + in - i_tor_min +1

                                  Rhs_loc(n_tor_local*n_var * (index_node2-1) + (kv-1)*n_tor_local + in - i_tor_min +1) = ZBIG * delta_psi_rmp_dt

                                  irn(ilarge_vp2) =  n_tor_local * n_var * (index_node2-1) + (kv-1)*n_tor_local + in - i_tor_min +1
                                  jcn(ilarge_vp2) =  n_tor_local * n_var * (index_node2-1) + (kp-1)*n_tor_local + in - i_tor_min +1
                                  A_mat(ilarge_vp2)   = ZBIG

                               endif
                            endif
                         endif

                         !======================================= end RMPs ==================================


                         if (                                                      &
                              ((freeboundary) .and. (k .eq. 1) .and. (in .eq. 1))  &               ! exclude condition on psi (freeboundary) except n=0
                              .or. (( .not. freeboundary) .and. (k .eq. 1) .and. (.not. RMP_on) .and. ( in .ge. 2 ))   &
                              .or. (( .not. freeboundary) .and. (k .eq. 1) .and. (in .eq. 1))  &
                              .or. (( .not. freeboundary) .and. (k .eq. 1) .and. (in .ge. 4))  & 
                              .or. (k .eq. 2)    &
                              .or. (k .eq. 3)    &
                              .or. (k .eq. 4)    &
                              .or. (k .eq. 5)    &
                              .or. (k .eq. 6)    &
                              .or. (k .eq. 7)    &
                              ) then

                            index_node = node_list%node(inode)%index(1)
                               if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then

                                  call locate_irn_jcn(index_node,index_node,index_min,index_max,ijA_position,ijA_index, ijA_size, irn_jcn)

                                  ilarge2 = ijA_position - 1 + ((k-1)*n_tor_local + in-i_tor_min) * n_var*n_tor_local &
                                                             +  (k-1)*n_tor_local + in - i_tor_min +1

                                  irn(ilarge2) =  n_tor_local * n_var * (index_node-1) + (k-1)*n_tor_local + in - i_tor_min +1
                                  jcn(ilarge2) =  n_tor_local * n_var * (index_node-1) + (k-1)*n_tor_local + in - i_tor_min +1
                                  A_mat(ilarge2)   = zbig

                               endif
                            index_node = node_list%node(inode)%index(3)

                               if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then

                                  call locate_irn_jcn(index_node,index_node,index_min,index_max,ijA_position,ijA_index, ijA_size, irn_jcn)

                                  ilarge2 = ijA_position - 1 + ((k-1)*n_tor_local + in-i_tor_min) * n_var*n_tor_local &
                                                             +  (k-1)*n_tor_local + in - i_tor_min +1

                                  irn(ilarge2) =  n_tor_local * n_var * (index_node-1) + (k-1)*n_tor_local + in - i_tor_min +1
                                  jcn(ilarge2) =  n_tor_local * n_var * (index_node-1) + (k-1)*n_tor_local + in - i_tor_min +1
                                  A_mat(ilarge2)    = zbig

                               endif
                         endif

                      endif

                      !------------------------------------ Special corners (only for grid with patches)
                      if    ((node_list%node(inode)%boundary .eq. 21) &
                        .or. (node_list%node(inode)%boundary .eq. 20)) then

!====================================== begining RMPs at boundary ======================================================
!================================== type 20-21 - boundary: corners apply on 's' and 't'
! ======================================================================================================================
                       
                        if (RMP_on ) then

                           if ((k.eq.1) .and. ((in.eq.2) .or. (in.eq.3)) .and. (.not. freeboundary)) then
                              ! in .eq. 2 corresponds to cos(n_perturbation)
                              ! in .eq. 3 corresponds to sin(n_perturbation)

                              kp=1    ! variable psi
                              kv=1    ! equation for psi

                              index_node = node_list%node(inode)%index(1)  ! index in RHS (or matrix A not compressed)

                              Rnode     = node_list%node(inode)%x(1,1,1) 
                              dRnode_ds = node_list%node(inode)%x(1,2,1) 
                              dRnode_dt = node_list%node(inode)%x(1,3,1) 
                              Znode     = node_list%node(inode)%x(1,1,2) 
                              dZnode_ds = node_list%node(inode)%x(1,2,2) 
                              dZnode_dt = node_list%node(inode)%x(1,3,2) 


                              if (in.eq.2) then
                                 delta_psi_rmp = psi_RMP_cos1(node_list%node(inode)%boundary_index)
                                 delta_psi_rmp_dR = dpsi_RMP_cos_dR1(node_list%node(inode)%boundary_index)
                                 delta_psi_rmp_dZ = dpsi_RMP_cos_dZ1(node_list%node(inode)%boundary_index)
                                 if (node_list%node(inode)%boundary_index == 1 ) then
                                    write (*,*) 'type2_bnd: my_id, psi_RMP_cos1, Rnode, Znode'
                                    write (*,*) my_id, delta_psi_rmp, Rnode, Znode
                                    write (*,*) 'delta_psi_rmp_dR, delta_psi_rmp_dZ'      
                                    write (*,*) delta_psi_rmp_dR, delta_psi_rmp_dZ
                                 endif
                              else 
                                 delta_psi_rmp = psi_RMP_sin1(node_list%node(inode)%boundary_index)
                                 delta_psi_rmp_dR = dpsi_RMP_sin_dR1(node_list%node(inode)%boundary_index)
                                 delta_psi_rmp_dZ = dpsi_RMP_sin_dZ1(node_list%node(inode)%boundary_index)
                              endif

                              delta_psi_rmp_ds = delta_psi_rmp_dR * dRnode_ds + delta_psi_rmp_dZ * dZnode_ds
                              delta_psi_rmp_dt = delta_psi_rmp_dR * dRnode_dt + delta_psi_rmp_dZ * dZnode_dt
                              if (in.eq.2) then
                                 if (node_list%node(inode)%boundary_index == 1 ) then
                                    write (*,*) 'delta_psi_rmp_dt', delta_psi_rmp_dt
                                    write (*,*) 'delta_psi_rmp_ds', delta_psi_rmp_ds
                                 endif
                              endif

                              if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then
                                 call locate_irn_jcn(index_node,index_node,index_min,index_max,ijA_position,ijA_index, ijA_size, irn_jcn)
                                 ilarge_vp  = ijA_position  - 1 + ((kv-1)*n_tor_local +in-i_tor_min) * n_var*n_tor_local & 
                                                                +  (kp-1)*n_tor_local + in - i_tor_min +1
                                 Rhs_loc(n_tor_local*n_var * (index_node-1) + (kv-1)*n_tor_local + in - i_tor_min +1) = ZBIG * delta_psi_rmp
                                 irn(ilarge_vp) =  n_tor_local * n_var * (index_node-1) + (kv-1)*n_tor_local + in - i_tor_min +1
                                 jcn(ilarge_vp) =  n_tor_local * n_var * (index_node-1) + (kp-1)*n_tor_local + in - i_tor_min +1
                                 A_mat(ilarge_vp)   = ZBIG
                              endif

                              index_node2 = node_list%node(inode)%index(2)

                              if ((index_node2 .ge. index_min) .and. (index_node2 .le. index_max)) then                         
                                 call locate_irn_jcn(index_node2,index_node2,index_min,index_max,ijA_position2,ijA_index, ijA_size, irn_jcn)
                                 ilarge_vp2  = ijA_position2  - 1 + ((kv-1)*n_tor_local +in-i_tor_min) * n_var*n_tor_local & 
                                                                  +  (kp-1)*n_tor_local + in - i_tor_min +1
                                 Rhs_loc(n_tor_local*n_var * (index_node2-1) + (kv-1)*n_tor_local + in - i_tor_min +1) = ZBIG * delta_psi_rmp_ds
                                 irn(ilarge_vp2) =  n_tor_local * n_var * (index_node2-1) + (kv-1)*n_tor_local + in - i_tor_min +1
                                 jcn(ilarge_vp2) =  n_tor_local * n_var * (index_node2-1) + (kp-1)*n_tor_local + in - i_tor_min +1
                                 A_mat(ilarge_vp2)   = ZBIG
                              endif
                              
                              index_node2 = node_list%node(inode)%index(3)

                              if ((index_node2 .ge. index_min) .and. (index_node2 .le. index_max)) then                         
                                 call locate_irn_jcn(index_node2,index_node2,index_min,index_max,ijA_position2,ijA_index, ijA_size, irn_jcn)
                                 ilarge_vp2  = ijA_position2  - 1 + ((kv-1)*n_tor_local +in-i_tor_min) * n_var*n_tor_local & 
                                                                  +  (kp-1)*n_tor_local + in - i_tor_min +1
                                 Rhs_loc(n_tor_local*n_var * (index_node2-1) + (kv-1)*n_tor_local + in - i_tor_min +1) = ZBIG * delta_psi_rmp_dt
                                 irn(ilarge_vp2) =  n_tor_local * n_var * (index_node2-1) + (kv-1)*n_tor_local + in - i_tor_min +1
                                 jcn(ilarge_vp2) =  n_tor_local * n_var * (index_node2-1) + (kp-1)*n_tor_local + in - i_tor_min +1
                                 A_mat(ilarge_vp2)   = ZBIG
                              endif
                              
                           endif
                        endif

                        !======================================= end RMPs ==================================

                        ! decides when the boundary conditions should be applied (for freeboundary and RMP cases)
                        if (                                                      &
                             ((freeboundary) .and. (k .eq. 1) .and. (in .eq. 1))  &               ! exclude condition on psi (freeboundary) except n=0
                             .or. (( .not. freeboundary) .and. (k .eq. 1) .and. (.not. RMP_on) .and. ( in .ge. 2 ))   &
                             .or. (( .not. freeboundary) .and. (k .eq. 1) .and. (in .eq. 1))  &
                             .or. (( .not. freeboundary) .and. (k .eq. 1) .and. (in .ge. 4))  & 
                             .or. (k .eq. 2)    &
                             .or. (k .eq. 3)    &
                             .or. (k .eq. 4)    &
                             .or. (k .eq. 5)    &
                             .or. (k .eq. 6)    &
                             .or. (k .eq. 7)    &
                             ) then

                          index_node = node_list%node(inode)%index(1)
                          if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then
                            call locate_irn_jcn(index_node,index_node,index_min,index_max,ijA_position,ijA_index, ijA_size, irn_jcn)
                            ilarge2 = ijA_position - 1 + ((k-1)*n_tor_local +in-i_tor_min) * n_var*n_tor_local & 
                                                       +  (k-1)*n_tor_local + in - i_tor_min +1
                            irn(ilarge2) =  n_tor_local * n_var * (index_node-1) + (k-1)*n_tor_local + in - i_tor_min +1
                            jcn(ilarge2) =  n_tor_local * n_var * (index_node-1) + (k-1)*n_tor_local + in - i_tor_min +1
                            A_mat(ilarge2)   = zbig
                          endif
                          
                          index_node = node_list%node(inode)%index(2)
                          if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then
                            call locate_irn_jcn(index_node,index_node,index_min,index_max,ijA_position,ijA_index, ijA_size, irn_jcn)
                            ilarge2 = ijA_position - 1 + ((k-1)*n_tor_local +in-i_tor_min) * n_var*n_tor_local &
                                                       +  (k-1)*n_tor_local + in - i_tor_min +1
                            irn(ilarge2) =  n_tor_local * n_var * (index_node-1) + (k-1)*n_tor_local + in - i_tor_min +1
                            jcn(ilarge2) =  n_tor_local * n_var * (index_node-1) + (k-1)*n_tor_local + in - i_tor_min +1
                            A_mat(ilarge2)    = zbig
                          endif

                          index_node = node_list%node(inode)%index(3)
                          if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then
                            call locate_irn_jcn(index_node,index_node,index_min,index_max,ijA_position,ijA_index, ijA_size, irn_jcn)
                            ilarge2 = ijA_position - 1 + ((k-1)*n_tor_local +in-i_tor_min) * n_var*n_tor_local &
                                                       +  (k-1)*n_tor_local + in - i_tor_min +1
                            irn(ilarge2) =  n_tor_local * n_var * (index_node-1) + (k-1)*n_tor_local + in - i_tor_min +1
                            jcn(ilarge2) =  n_tor_local * n_var * (index_node-1) + (k-1)*n_tor_local + in - i_tor_min +1
                            A_mat(ilarge2)    = zbig
                          endif

                        endif

                      endif

                   enddo

                enddo
             endif
          enddo
       enddo

    if (RMP_on) then
       if (allocated(psi_RMP_cos1))         call tr_deallocate(psi_RMP_cos1,"psi_RMP_cos1",CAT_UNKNOWN)
       if (allocated(dpsi_RMP_cos_dR1))     call tr_deallocate(dpsi_RMP_cos_dR1,"dpsi_RMP_cos_dR1",CAT_UNKNOWN)
       if (allocated(dpsi_RMP_cos_dZ1))     call tr_deallocate(dpsi_RMP_cos_dZ1,"dpsi_RMP_cos_dZ1",CAT_UNKNOWN)
       if (allocated(psi_RMP_sin1))         call tr_deallocate(psi_RMP_sin1,"psi_RMP_sin1",CAT_UNKNOWN)
       if (allocated(dpsi_RMP_sin_dR1))     call tr_deallocate(dpsi_RMP_sin_dR1,"dpsi_RMP_sin_dR1",CAT_UNKNOWN)
       if (allocated(dpsi_RMP_sin_dZ1))     call tr_deallocate(dpsi_RMP_sin_dZ1,"dpsi_RMP_sin_dZ1",CAT_UNKNOWN)
    endif

    return
  end subroutine boundary_conditions
end module mod_boundary_conditions
