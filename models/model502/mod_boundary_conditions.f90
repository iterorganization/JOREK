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
!*   index_max    - Maximal index of local elements (                          *
!*   xpoint2      -                                                            *
!*   xcase2       -                                                            *
!*   psi_axis     -                                                            *
!*   psi_bnd      -                                                            *
!*   Z_xpoint     -                                                            *
!*   gmres        - boolean indicating if we are using GMRES method            *
!*   solve_only   - Indicate if we want to perform only solve                  *
!*                                                                             *
!*******************************************************************************
  subroutine boundary_conditions( my_id, node_list, element_list, bnd_node_list, local_elms,& 
                                  n_local_elms, index_min, index_max, rhs_loc, xpoint2,     &
                                  xcase2, R_axis, Z_axis, psi_axis, psi_bnd,                &
                                  R_xpoint, Z_xpoint, psi_xpoint, gmres, solve_only,        & 
                                  ijA_index, ijA_size, irn_jcn, irn, jcn, A_mat, i_tor_min, i_tor_max )

    use mod_assembly, only : boundary_conditions_add_one_entry, boundary_conditions_add_RHS
    use data_structure
    use vacuum, ONLY: is_freebound
    use phys_module, only: F0, GAMMA, freeboundary, RMP_on, psi_RMP_cos, dpsi_RMP_cos_dR, dpsi_RMP_cos_dZ, &
       psi_RMP_sin, dpsi_RMP_sin_dR, dpsi_RMP_sin_dZ, t_now, RMP_growth_rate, RMP_ramp_up_time,            &
       RMP_start_time, tstep, RMP_har_cos, RMP_har_sin, T_min,                                             &
       mach_one_bnd_integral, Vpar_smoothing, vpar_smoothing_coef,                                         &
       Number_RMP_harmonics, RMP_har_cos_spectrum,RMP_har_sin_spectrum, grid_to_wall, n_wall_blocks, keep_n0_const
    USE tr_module
    use mpi_mod
    use mod_locate_irn_jcn
    use mod_basisfunctions
    use mod_interp
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
  real*8,                             intent(inout) :: rhs_loc(*)
  integer,                            intent(in)    :: i_tor_min, i_tor_max 
  integer(kind=int_all), allocatable, intent(in)    :: ijA_index(:,:)
  integer(kind=int_all), allocatable, intent(in)    :: ijA_size(:)
  integer(kind=int_all), allocatable, intent(in)    :: irn_jcn(:,:) 
  integer(kind=int_all), allocatable, intent(inout) :: irn(:)
  integer(kind=int_all), allocatable, intent(inout) :: jcn(:) 
  real*8,                allocatable, intent(inout) :: A_mat(:) 

  ! Internal parameters
  real*8                :: zbig, zbig_backup,  T0, Vpar0, bigR, dT0_ds, dVpar0_ds, dBigR_ds, psi_1, R_1, Z_1
  real*8                :: R_s, R_t, Z, Z_s, Z_t, R_tt, Z_tt, ps0, ps0_s, ps0_t, ps0_tt, ps0_x, ps0_y, direction, xjac
  real*8                :: Btot, alpha, dT0_dt, dVpar0_dt, dBigR_dt, R_inside, Z_inside
  real*8                :: grad_psi, u0_s, u0_t, u0_x, u0_y, element_size_0, element_size_2, element_size_perp
  real*8                :: H1(2,2), H1_s(2,2), H1_ss(2,2)
  integer               :: i, in, iv, iv2, inode, inode2, k
  integer               :: ielm
  integer               :: index_node, index_node2
  integer(kind=int_all) :: ijA_position,ijA_position2
  integer               :: ilarge2, kv, kT, ku, kn, ilarge_vv, ilarge_vT, ilarge_vus, ilarge_vn
  integer               :: ilarge_vsvs, ilarge_vsTs, ilarge_vsT, ilarge_vut, ilarge_vtvt, ilarge_vtTt, ilarge_vtT
  integer               :: ierr
  logical               :: apply_psi_BC, apply_current_BC

  real*8, allocatable   :: psi_RMP_cos1(:),dpsi_RMP_cos_dR1(:),dpsi_RMP_cos_dZ1(:)
  real*8, allocatable   :: psi_RMP_sin1(:),dpsi_RMP_sin_dR1(:),dpsi_RMP_sin_dZ1(:)
  real*8                :: Rnode, dRnode_ds, Znode, dZnode_ds, dRnode_dt, dZnode_dt, establish_RMP
  real*8                :: delta_psi_rmp, delta_psi_rmp_dR, delta_psi_rmp_dZ, delta_psi_rmp_ds, delta_psi_rmp_dt, psi_test, sigmo_fonc
  real*8                :: R_mid, Z_mid, R_center, Z_center, direction2, normal(2), normal_direction(2), grad_s(2), grad_t(2)
  real*8                :: factor, factor_t, c_1, c_2, c_3, bn, cs0, cs0_T, cs0_TT, dl, dl_dt, bn_t, bn_t_abs, hfact_t
  integer               :: ilarge_vp, ilarge_vp2
  integer               :: kp, j, err, itest, i_mid, i_bnd, inode_m, inode_p
  integer               :: n_rmp_harm, N_rmp_har_block_size

  real*8                :: R_out, Z_out, s_elm, t_elm, QR,QR_s,QR_t,QR_st,QR_ss,QR_tt,QZ,QZ_s,QZ_t,QZ_st,QZ_ss,QZ_tt
  real*8                :: QPs0,QPs0_s,QPs0_t,QPs0_st,QPs0_ss,QPs0_tt
  integer               :: ifail, i_elm

  RMPspectrum: if (RMP_on .and. (n_tor .ge. 3)) then !*****
  
! for the moment it's done in a way that all RMP harmonics follow each other,i.e. n=2,n=3,n=4... 
! if you want for example n=2 and n=4 RMP you should consider n=2,3,4, but put zeros at the boundary in the input file for n=3 RMP
! example: ntor=13 and nperiod=1(so taking into account, toroidal numbers n=0,1,2....6) and  n=2 and n=3 are toroidal numbers of RMPs, 
! so Number_RMP_harmonics=2, RMP_har_cos_spectrum(1)=4,RMP_har_sin_spectrum(1)=5,RMP_har_cos_spectrum(2)=6,RMP_har_sin_spectrum(2)=7.  
  
    call tr_allocate(psi_RMP_cos1,1, bnd_node_list%n_bnd_nodes*Number_RMP_harmonics,"psi_RMP_cos1",CAT_UNKNOWN)
    call tr_allocate(dpsi_RMP_cos_dR1,1,bnd_node_list%n_bnd_nodes*Number_RMP_harmonics,"dpsi_RMP_cos_dR1",CAT_UNKNOWN)
    call tr_allocate(dpsi_RMP_cos_dZ1,1,bnd_node_list%n_bnd_nodes*Number_RMP_harmonics,"dpsi_RMP_cos_dZ1",CAT_UNKNOWN)
    call tr_allocate(psi_RMP_sin1,1,bnd_node_list%n_bnd_nodes*Number_RMP_harmonics,"psi_RMP_sin1",CAT_UNKNOWN)
    call tr_allocate(dpsi_RMP_sin_dR1,1,bnd_node_list%n_bnd_nodes*Number_RMP_harmonics,"dpsi_RMP_sin_dR1",CAT_UNKNOWN)
    call tr_allocate(dpsi_RMP_sin_dZ1,1,bnd_node_list%n_bnd_nodes*Number_RMP_harmonics,"dpsi_RMP_sin_dZ1",CAT_UNKNOWN)
    N_rmp_har_block_size=bnd_node_list%n_bnd_nodes
    
    psi_test =  node_list%node(bnd_node_list%bnd_node(1)%index_jorek)%values(RMP_har_cos_spectrum(1),1,1)
    ! if necessary, replace by:
    ! psi_test =  node_list%node(bnd_node_list%bnd_node(1)%index_jorek)%values(min(RMP_har_cos_spectrum(1), n_tor),1,1)
    write (*,*) 'psi_bnd at previous time step', psi_test
    
    if (abs(psi_test) .le. abs(psi_RMP_cos(1))) then
      sigmo_fonc = ( 1.d0 + exp(-RMP_growth_rate*( t_now - RMP_start_time - RMP_ramp_up_time/2.d0 )))**(-1) &
          - ( 1.d0 + exp(-RMP_growth_rate*( 0.d0 - RMP_ramp_up_time/2.d0 )))**(-1) 
      establish_RMP = (RMP_growth_rate*sigmo_fonc*(1-sigmo_fonc)+1.e-6)*tstep 
    else
      establish_RMP = 0.d0
    endif
    ! Other possibility (simpler) : if ( (t_now - RMP_start_time) .ge. 2.2*RMP_ramp_up_time/2.d0 ) then establish_RMP =0.0
  
    do j=1, bnd_node_list%n_bnd_nodes*Number_RMP_harmonics  
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

  end if RMPspectrum

  zbig = 1.d12
  zbig_backup = zbig

     do i=1, n_local_elms !===============================do elements

        ielm = local_elms(i)

        i_bnd = 0
        do iv=1, n_vertex_max !==========================do vertex
          inode = element_list%element(ielm)%vertex(iv)
          if (node_list%node(inode)%boundary .ne. 0) i_bnd = i_bnd + 1
        enddo
        if (i_bnd .lt. 2) cycle           

        i_mid = 0; R_mid = 0.d0; Z_mid = 0.d0; R_center = 0.d0; Z_center = 0.d0

        do iv=1, n_vertex_max !==========================do vertex
          inode = element_list%element(ielm)%vertex(iv)
          if (node_list%node(inode)%boundary .ne. 0) then
            i_mid = i_mid + 1
            R_mid = R_mid + node_list%node(inode)%x(1,1,1)     ! mid point on boundary (approx.)
            Z_mid = Z_mid + node_list%node(inode)%x(1,1,2)     ! mid point on boundary (approx.)
          endif
          R_center = R_center + node_list%node(inode)%x(1,1,1)       ! center point within element (approx.)
          Z_center = Z_center + node_list%node(inode)%x(1,1,2)       ! center point within element (approx.)
        enddo
        R_mid    = R_mid    / real(i_mid,8)
        Z_mid    = Z_mid    / real(i_mid,8)
        R_center = R_center / real(n_vertex_max,8)
        Z_center = Z_center / real(n_vertex_max,8)
         
        normal_direction = (/R_mid - R_center, Z_mid - Z_center /) / norm2((/R_mid - R_center, Z_mid - Z_center /))

        do iv=1, n_vertex_max !==========================do vertex

           inode = element_list%element(ielm)%vertex(iv)

           if (node_list%node(inode)%boundary .ne. 0) then !==================if boundary nodes

             R_mid = node_list%node(inode)%x(1,1,1)
             Z_mid = node_list%node(inode)%x(1,1,2)

             inode_p = element_list%element(ielm)%vertex(mod(iv  ,4) + 1)
             inode_m = element_list%element(ielm)%vertex(mod(iv+2,4) + 1)

             if (node_list%node(inode_p)%boundary .eq. 0) then 
               R_center = node_list%node(inode_p)%x(1,1,1)
               Z_center = node_list%node(inode_p)%x(1,1,2)
               iv2 = mod(iv+2,4) + 1                                  ! the index of the other boundary vertex
               inode2 = inode_m
             elseif (node_list%node(inode_m)%boundary .eq. 0) then 
               R_center = node_list%node(inode_m)%x(1,1,1)
               Z_center = node_list%node(inode_m)%x(1,1,2)
               iv2 = mod(iv  ,4) + 1                                  ! the index of the other boundary vertex
               inode2 = inode_p
             endif

             normal_direction = (/R_mid - R_center, Z_mid - Z_center /) / norm2((/R_mid - R_center, Z_mid - Z_center /))

              do in=i_tor_min, i_tor_max  !========================do n_tor
                if (keep_n0_const  .and.  in .eq. 1 ) then
                  zbig = 1.d15
                else
                  zbig = zbig_backup
                endif

!              do n_rmp_harm=1, Number_RMP_harmonics !===========do RMP harmonics

                 do k=1, n_var ! ================================do variables
                 
                                                                                      !-----(General for all bnd types)
                   !------------ Decide when Psi or Current need BCs --------------------------------------------------                      
                   !----Psi
                   apply_psi_BC = .false.
                   if (k == 1) then                        
                     if ( (RMP_on) .and. (in .lt. RMP_har_cos_spectrum(1))                    )   apply_psi_BC = .true.
                     if ( (RMP_on) .and. (in .gt. RMP_har_sin_spectrum(Number_RMP_harmonics)) )   apply_psi_BC = .true.
                     if ( (.not. RMP_on) .and. (in .ge. 2)              )                         apply_psi_BC = .true.
                     if (              in .eq. 1                        )                         apply_psi_BC = .true.
                     if (           is_freebound(in,k)                  )                         apply_psi_BC = .false.                     
                   endif
                      
                   !----Current
                   apply_current_BC = .false.
                   if (k == 3) then
                     if ( .not. is_freebound(in,k) )   apply_current_BC = .true.
                   endif
                   !---------------------------------------------------------------------------------------------------

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
!=========================================================RMP spectrum========================================Marina
                          do n_rmp_harm=1, Number_RMP_harmonics !===========do RMP harmonics

                          if ((k.eq.1) .and. ((in.eq.RMP_har_cos_spectrum(n_rmp_harm)) .or. (in.eq.RMP_har_sin_spectrum(n_rmp_harm))) &
                              .and. (.not. freeboundary)) then
                             ! in .eq. RMP_har_cos corresponds to cos(n_perturbation)
                             ! in .eq. RMP_har_sin corresponds to sin(n_perturbation)

                             kp=1    ! variable psi
                             kv=1    ! equation for psi

                             index_node = node_list%node(inode)%index(1)  ! index in RHS (or matrix A not compressed)
                                                
                             Rnode     = node_list%node(inode)%x(1,1,1) 
                             dRnode_ds = node_list%node(inode)%x(1,2,1) 
                             Znode     = node_list%node(inode)%x(1,1,2) 
                             dZnode_ds = node_list%node(inode)%x(1,2,2) 
                          
                             if (in.eq.RMP_har_cos_spectrum(n_rmp_harm)) then
                                delta_psi_rmp = psi_RMP_cos1(node_list%node(inode)%boundary_index +N_rmp_har_block_size*(n_rmp_harm-1))
                                delta_psi_rmp_dR = dpsi_RMP_cos_dR1(node_list%node(inode)%boundary_index+N_rmp_har_block_size*(n_rmp_harm-1))
                                delta_psi_rmp_dZ = dpsi_RMP_cos_dZ1(node_list%node(inode)%boundary_index+N_rmp_har_block_size*(n_rmp_harm-1))

                             else 
                                delta_psi_rmp = psi_RMP_sin1(node_list%node(inode)%boundary_index+N_rmp_har_block_size*(n_rmp_harm-1))
                                delta_psi_rmp_dR = dpsi_RMP_sin_dR1(node_list%node(inode)%boundary_index+N_rmp_har_block_size*(n_rmp_harm-1))
                                delta_psi_rmp_dZ = dpsi_RMP_sin_dZ1(node_list%node(inode)%boundary_index+N_rmp_har_block_size*(n_rmp_harm-1))

                             endif

                             delta_psi_rmp_ds = delta_psi_rmp_dR * dRnode_ds + delta_psi_rmp_dZ * dZnode_ds
                             call boundary_conditions_add_one_entry(   &
                                  index_node, kv, in,                  &
                                  index_node, kp, in,                  &
                                  zbig, solve_only, gmres,             &
                                  index_min, index_max,                & 
                                  ijA_index, ijA_size, irn_jcn,        & 
                                  irn, jcn, A_mat, i_tor_min, i_tor_max)

                                call boundary_conditions_add_RHS(      &
                                     index_node, kv, in,               &
                                     index_min, index_max,             &
                                     RHS_loc, ZBIG * delta_psi_rmp,    &
                                     i_tor_min, i_tor_max)
                             
                             index_node2 = node_list%node(inode)%index(2)

                             call boundary_conditions_add_one_entry(   &
                                  index_node2, kv, in,                 &
                                  index_node2, kp, in,                 &
                                  zbig, solve_only, gmres,             &
                                  index_min, index_max,                & 
                                  ijA_index, ijA_size, irn_jcn,        & 
                                  irn, jcn, A_mat, i_tor_min, i_tor_max)
                                call boundary_conditions_add_RHS(      &
                                     index_node2, kv, in,              &
                                     index_min, index_max,             &
                                     RHS_loc, ZBIG * delta_psi_rmp_ds, &
                                     i_tor_min, i_tor_max)
                          endif
                       enddo  !(end RMP harmonics)   
                       endif !(end RMP)
!======================================= end RMPs ==================================

                      
                       if (        apply_psi_BC      &
                              .or. apply_current_BC  &
                              .or. (k .eq. 2)        &
                              .or. (k .eq. 4)        &
                             !.or. (k .eq. 5)        &
                             !.or. (k .eq. 6)        &
                             !.or. (k .eq. 7)        &
                             !.or. (k .eq. 8)        &
                              ) then


                            index_node = node_list%node(inode)%index(1)

                            call boundary_conditions_add_one_entry(   &
                                 index_node, k, in,                   &
                                 index_node, k, in,                   &
                                 zbig, solve_only, gmres,             &
                                 index_min, index_max,                & 
                                 ijA_index, ijA_size, irn_jcn,        & 
                                 irn, jcn, A_mat, i_tor_min, i_tor_max)

                            index_node = node_list%node(inode)%index(2)

                            call boundary_conditions_add_one_entry(   &
                                 index_node, k, in,                   &
                                 index_node, k, in,                   &
                                 zbig, solve_only, gmres,             &
                                 index_min, index_max,                & 
                                 ijA_index, ijA_size, irn_jcn,        & 
                                 irn, jcn, A_mat, i_tor_min, i_tor_max)
                         endif

                         if ( (k .eq. 7) .and. (.not. mach_one_bnd_integral) ) then

                            index_node  = node_list%node(inode)%index(1)             ! position of value
                            index_node2 = node_list%node(inode)%index(2)             ! position of first deriative

                            T0        = max(node_list%node(inode)%values(1,1,6), T_min)
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

                            if (xcase2 .eq. 2) then
                               direction = -direction
                            else if ((xcase2 .eq. 3).and.(node_list%node(inode)%x(1,1,2).gt.Z_axis +0.1).and.(node_list%node(inode)%x(1,1,1).gt.R_xpoint(2))) then
                              direction = -1.
                            else if ((xcase2 .eq. 3) .and. (node_list%node(inode)%x(1,1,2).gt.Z_axis +0.1).and.(node_list%node(inode)%x(1,1,1).lt.R_xpoint(2)))then
                              direction = +1.
                            end if
                            
                            ! --- Special field direction for grid patches
                            ! --- (hopefully temporary until Vpar is treated in boundary_matrix_open)
                            if ( (grid_to_wall) .and. (n_wall_blocks .gt. 0) ) then
                              direction = 1
                              if (node_list%node(inode)%boundary .eq. 11) direction = -direction
                              if (node_list%node(inode)%boundary .eq. 19) direction = -direction
                            endif

                            grad_t = (/ -Z_s,   R_s /) / xjac
                          
                            normal     = dot_product(grad_t,normal_direction) * grad_t      ! outward pointing normal
                            normal     = normal / norm2(normal)
                            direction  = sign(1.d0,dot_product((/ps0_y,-ps0_x/),normal))
                            
                            grad_psi = sqrt(ps0_x**2 + ps0_y**2)

                            Btot = sqrt(F0**2 + ps0_x**2 + ps0_y**2) / BigR
                            
!                            if (in .eq. 1) then                              
!                              write(*,'(i3,A,3e14.6,A,e14.6)') node_list%node(inode)%boundary, &
!                                                              ' Boundary (s): ',Vpar0, -BigR**2 * u0_s/ps0_s, direction*sqrt(GAMMA*T0)/Btot, &
!                                                              ' error : ',Vpar0 - BigR**2 * u0_s/ps0_s - direction*sqrt(GAMMA*T0)/Btot
!                            endif

                            ku = 2
                            kv = 7
                            kT = 6

                            call boundary_conditions_add_one_entry(   &
                                 index_node, kv, in,                  &
                                 index_node, kv, in,                  &
                                 zbig, solve_only, gmres,             &
                                 index_min, index_max,                & 
                                 ijA_index, ijA_size, irn_jcn,        & 
                                 irn, jcn, A_mat, i_tor_min, i_tor_max)

                            call boundary_conditions_add_one_entry(   &
                                 index_node, kv, in,                  &
                                 index_node, kT, in,                  &
                                 - zbig / Btot * 0.5d0 * GAMMA        &
                                 / sqrt(GAMMA*T0) * direction         &
                                 , solve_only, gmres,                 &
                                 index_min, index_max,                & 
                                 ijA_index, ijA_size, irn_jcn,        & 
                                 irn, jcn, A_mat, i_tor_min, i_tor_max)

                            call boundary_conditions_add_one_entry(   &
                                 index_node,  kv, in,                 &
                                 index_node2, ku, in,                 &
                                 - zbig * BigR**2 / ps0_s,            &
                                 solve_only, gmres,                   &
                                 index_min, index_max,                & 
                                 ijA_index, ijA_size, irn_jcn,        & 
                                 irn, jcn, A_mat, i_tor_min, i_tor_max)

                            if (in .eq. 1) then
                               call boundary_conditions_add_RHS(              &
                                    index_node, kv, in,                       &
                                    index_min, index_max,                     &
                                    RHS_loc,                                  &
                                    Zbig * (-Vpar0 + BigR**2 * U0_s/ps0_s +   & 
                                           direction * sqrt(GAMMA*T0) / Btot),& 
                                    i_tor_min, i_tor_max)

                            else
                               call boundary_conditions_add_RHS(              &
                                    index_node, kv, in,                       &
                                    index_min, index_max,                     &
                                    RHS_loc, 0.d0,                            &
                                    i_tor_min, i_tor_max)

                            endif
                            index_node  = node_list%node(inode)%index(1)
                            index_node2 = node_list%node(inode)%index(2)
                            kv = 7
                            kT = 6

                            call boundary_conditions_add_one_entry(   &
                                 index_node2, kv, in,                 &
                                 index_node2, kv, in,                 &
                                 zbig, solve_only, gmres,             &
                                 index_min, index_max,                & 
                                 ijA_index, ijA_size, irn_jcn,        & 
                                 irn, jcn, A_mat, i_tor_min, i_tor_max)

                            call boundary_conditions_add_one_entry(   &
                                 index_node2, kv, in,                 &
                                 index_node2, kT, in,                 &
                                 - zbig / Btot * 0.5d0 * GAMMA        &
                                 / sqrt(GAMMA*T0) * direction,        &
                                 solve_only, gmres,                   &
                                 index_min, index_max,                & 
                                 ijA_index, ijA_size, irn_jcn,        & 
                                 irn, jcn, A_mat, i_tor_min, i_tor_max)

                            call boundary_conditions_add_one_entry(   &
                                 index_node2,  kv, in,                &
                                 index_node,   kT, in,                &
                                 + zbig / Btot * 0.25d0 * GAMMA**2    &
                                 / (GAMMA*T0)**(3/2) * dT0_ds *       &
                                 direction,                           &
                                 solve_only, gmres,                   &
                                 index_min, index_max,                & 
                                 ijA_index, ijA_size, irn_jcn,        & 
                                 irn, jcn, A_mat, i_tor_min, i_tor_max)

                            if (in .eq. 1) then
                               call boundary_conditions_add_RHS(                  &
                                    index_node2, kv, in,                          &
                                    index_min, index_max,                         &
                                    RHS_loc,                                      &
                                    Zbig*(-dVpar0_ds +  0.5d0 / Btot *            &
                                    GAMMA / sqrt(GAMMA*T0) * dT0_ds * direction), &
                                    i_tor_min, i_tor_max)
                            else
                               call boundary_conditions_add_RHS(                  &
                                    index_node2, kv, in,                          &
                                    index_min, index_max,                         &
                                    RHS_loc, 0.d0,                                &
                                    i_tor_min, i_tor_max)

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

!====================================== begining RMPs at boundary ======================================================
!================================== type 2 - boundary: only depends on 't'
! ======================================================================================================================
                       
                       if (RMP_on ) then
                          do n_rmp_harm=1, Number_RMP_harmonics !===========do RMP harmonics

                          if ((k.eq.1) .and. ((in.eq.RMP_har_cos_spectrum(n_rmp_harm)) .or. (in.eq.RMP_har_sin_spectrum(n_rmp_harm))) .and. (.not. freeboundary)) then
                             ! in .eq. RMP_har_cos  corresponds to cos(n_perturbation)
                             ! in .eq. RMP_har_sin   corresponds to sin(n_perturbation)
                                                      
                             kp=1    ! variable psi
                             kv=1    ! equation for psi
                          
                             index_node = node_list%node(inode)%index(1)  ! index in RHS (or matrix A not compressed)
                                                
                             Rnode     = node_list%node(inode)%x(1,1,1) 
                             dRnode_dt = node_list%node(inode)%x(1,3,1) 
                             Znode     = node_list%node(inode)%x(1,1,2) 
                             dZnode_dt = node_list%node(inode)%x(1,3,2) 
                          
                             if (in.eq.RMP_har_cos_spectrum(n_rmp_harm)) then
                                delta_psi_rmp = psi_RMP_cos1(node_list%node(inode)%boundary_index+N_rmp_har_block_size*(n_rmp_harm-1))
                                delta_psi_rmp_dR = dpsi_RMP_cos_dR1(node_list%node(inode)%boundary_index+N_rmp_har_block_size*(n_rmp_harm-1))
                                delta_psi_rmp_dZ = dpsi_RMP_cos_dZ1(node_list%node(inode)%boundary_index+N_rmp_har_block_size*(n_rmp_harm-1))

                                if (node_list%node(inode)%boundary_index == 1 ) then
                                   write (*,*) 'type2_bnd: my_id, psi_RMP_cos1, Rnode, Znode, in'
                                   write (*,*) my_id, delta_psi_rmp, Rnode, Znode,in
                                   write (*,*) 'delta_psi_rmp_dR, delta_psi_rmp_dZ'      
                                   write (*,*) delta_psi_rmp_dR, delta_psi_rmp_dZ
                                endif
                             else 
                                delta_psi_rmp = psi_RMP_sin1(node_list%node(inode)%boundary_index+N_rmp_har_block_size*(n_rmp_harm-1))
                                delta_psi_rmp_dR = dpsi_RMP_sin_dR1(node_list%node(inode)%boundary_index+N_rmp_har_block_size*(n_rmp_harm-1))
                                delta_psi_rmp_dZ = dpsi_RMP_sin_dZ1(node_list%node(inode)%boundary_index+N_rmp_har_block_size*(n_rmp_harm-1))
                             endif

                             delta_psi_rmp_dt = delta_psi_rmp_dR * dRnode_dt + delta_psi_rmp_dZ * dZnode_dt
                             if (in.eq.RMP_har_cos_spectrum(n_rmp_harm)) then

                                if (node_list%node(inode)%boundary_index == 1 ) then
                                   write (*,*) 'delta_psi_rmp_dt', delta_psi_rmp_dt
                                   write (*,*) 'delta_psi_rmp_ds', delta_psi_rmp_ds
                                endif
                             endif

                             call boundary_conditions_add_one_entry( &
                                  index_node, kv, in,                &
                                  index_node, kp, in,                &
                                  zbig, solve_only, gmres,           &
                                  index_min, index_max,              & 
                                  ijA_index, ijA_size, irn_jcn,      & 
                                  irn, jcn, A_mat, i_tor_min, i_tor_max)

                                call boundary_conditions_add_RHS(  &
                                     index_node, kv, in,           &
                                     index_min, index_max,         &
                                     RHS_loc, ZBIG * delta_psi_rmp,&
                                     i_tor_min, i_tor_max)
                             
                             index_node2 = node_list%node(inode)%index(3)

                             call boundary_conditions_add_one_entry( &
                                  index_node2, kv, in,               &
                                  index_node2, kp, in,               &
                                  zbig, solve_only, gmres,           &
                                  index_min, index_max,              & 
                                  ijA_index, ijA_size, irn_jcn,      & 
                                  irn, jcn, A_mat, i_tor_min, i_tor_max)
                                call boundary_conditions_add_RHS(       &
                                     index_node2, kv, in,               &
                                     index_min, index_max,              &
                                     RHS_loc, ZBIG * delta_psi_rmp_dt,  &
                                     i_tor_min, i_tor_max)
                          endif
                        enddo        !(end RMP harmonics)
                        endif        !(end RMPs on)  ==================================
!======================================= end RMPs ==================================


                         if (      apply_psi_BC      &
                              .or. apply_current_BC  &
                              .or. (k .eq. 2)        &
                              .or. (k .eq. 4)        &
                             !.or. (k .eq. 5)        &
                             !.or. (k .eq. 6)        &
                             !.or. (k .eq. 7)        &
                             !.or. (k .eq. 8)        &
                              ) then


                            index_node = node_list%node(inode)%index(1)

                            call boundary_conditions_add_one_entry(   &
                                 index_node,  k, in,                  &
                                 index_node,  k, in,                  &
                                 zbig, solve_only, gmres,             &
                                 index_min, index_max,                & 
                                 ijA_index, ijA_size, irn_jcn,        & 
                                 irn, jcn, A_mat, i_tor_min, i_tor_max)


                            index_node = node_list%node(inode)%index(3)
                            call boundary_conditions_add_one_entry(   &
                                 index_node,  k, in,                  &
                                 index_node,  k, in,                  &
                                 zbig, solve_only, gmres,             &
                                 index_min, index_max,                & 
                                 ijA_index, ijA_size, irn_jcn,        & 
                                 irn, jcn, A_mat, i_tor_min, i_tor_max)

                         endif


                         if ( (k .eq. 7) .and. (.not. mach_one_bnd_integral) ) then

                            call basisfunctions1(0.d0, H1, H1_s, H1_ss)

                            if ( (iv.eq.1) .or. (iv .eq.2)) then
                              element_size_0 =   element_list%element(ielm)%size(iv ,3) * H1_s(1,2)  ! this is a boundary along t
                            else
                              element_size_0 = - element_list%element(ielm)%size(iv ,3) * H1_s(1,2)  ! this is a boundary along t
                            endif

                            if (iv*iv2 .eq. 6) then
                              element_size_perp = - element_list%element(ielm)%size(iv,2) * H1_s(1,2)
                            else
                              element_size_perp = + element_list%element(ielm)%size(iv,2) * H1_s(1,2)
                            endif

                            index_node  = node_list%node(inode)%index(1)             ! position of value
                            index_node2 = node_list%node(inode)%index(3)             ! position of first deriative

                            T0        = max(node_list%node(inode)%values(1,1,6), T_min)
                            Vpar0     = node_list%node(inode)%values(1,1,7)
                            BigR      = node_list%node(inode)%x(1,1,1)

                            dT0_dt    = node_list%node(inode)%values(1,3,6) * element_size_0 
                            dVpar0_dt = node_list%node(inode)%values(1,3,7) * element_size_0 
                            dBigR_dt  = node_list%node(inode)%x(1,3,1)        * element_size_0

                            ps0       = node_list%node(inode)%values(1,1,1)
                            ps0_s     = node_list%node(inode)%values(1,2,1) * element_size_perp
                            ps0_t     = node_list%node(inode)%values(1,3,1) * element_size_0

                            U0_s      = node_list%node(inode)%values(1,2,2) * element_size_perp
                            U0_t      = node_list%node(inode)%values(1,3,2) * element_size_0   

                            R_s       = node_list%node(inode)%x(1,2,1) * element_size_perp
                            R_t       = node_list%node(inode)%x(1,3,1) * element_size_0    
                            Z_s       = node_list%node(inode)%x(1,2,2) * element_size_perp
                            Z_t       = node_list%node(inode)%x(1,3,2) * element_size_0    
                            Z         = node_list%node(inode)%x(1,1,2)
                            
                            ps0_tt = element_list%element(ielm)%size(iv ,1) * node_list%node(inode)%values(1,1,1)  * H1_ss(1,1) &
                                   + element_list%element(ielm)%size(iv ,3) * node_list%node(inode)%values(1,3,1)  * H1_ss(1,2) &
                                   + element_list%element(ielm)%size(iv2,1) * node_list%node(inode2)%values(1,1,1) * H1_ss(2,1) &
                                   + element_list%element(ielm)%size(iv2,3) * node_list%node(inode2)%values(1,3,1) * H1_ss(2,2)

                            R_tt = + element_list%element(ielm)%size(iv ,1) * node_list%node(inode )%x(1,1,1) * H1_ss(1,1)  &
                                   + element_list%element(ielm)%size(iv ,3) * node_list%node(inode )%x(1,3,1) * H1_ss(1,2)  &
                                   + element_list%element(ielm)%size(iv2,1) * node_list%node(inode2)%x(1,1,1) * H1_ss(2,1)  &
                                   + element_list%element(ielm)%size(iv2,3) * node_list%node(inode2)%x(1,3,1) * H1_ss(2,2)  

                            Z_tt = + element_list%element(ielm)%size(iv ,1) * node_list%node(inode )%x(1,1,2) * H1_ss(1,1)  &
                                   + element_list%element(ielm)%size(iv ,3) * node_list%node(inode )%x(1,3,2) * H1_ss(1,2)  &
                                   + element_list%element(ielm)%size(iv2,1) * node_list%node(inode2)%x(1,1,2) * H1_ss(2,1)  &
                                   + element_list%element(ielm)%size(iv2,3) * node_list%node(inode2)%x(1,3,2) * H1_ss(2,2)  
                            
!                            call find_RZ(node_list, element_list, BigR, Z, R_out, Z_out, i_elm, s_elm, t_elm, ifail)
!                            call interp_RZ_2(node_list,element_list,i_elm,s_elm,t_elm,QR,QR_s,QR_t,QR_st,QR_ss,QR_tt,QZ,QZ_s,QZ_t,QZ_st,QZ_ss,QZ_tt)
!                            call interp(node_list,element_list,i_elm,1,1,s_elm, t_elm, QPs0,QPs0_s,QPs0_t,QPs0_st,QPs0_ss,QPs0_tt)
!                            if (ielm .eq. i_elm) then
!                              write(112,'(3i8,20e14.6)') ielm,i_elm,iv,s_elm,t_elm,Qps0,ps0,Qps0_s,ps0_s, Qps0_t,ps0_t, Qps0_tt, ps0_tt
!                              write(113,'(3i8,20e14.6)') ielm,i_elm,iv,s_elm,t_elm,QR,bigR,QR_s,R_s, QR_t,R_t, QR_tt, R_tt
!                              write(114,'(3i8,20e14.6)') ielm,i_elm,iv,s_elm,t_elm,QZ,Z,QZ_s,Z_s, QZ_t,Z_t, QZ_tt, Z_tt
!                            endif

                            xjac  =  R_s*Z_t - R_t*Z_s
                            ps0_x = (   Z_t * ps0_s - Z_s * ps0_t ) / xjac
                            ps0_y = ( - R_t * ps0_s + R_s * ps0_t ) / xjac

                            u0_x = (   Z_t * u0_s - Z_s * u0_t ) / xjac
                            u0_y = ( - R_t * u0_s + R_s * u0_t ) / xjac

                            grad_s = (/  Z_t,  - R_t /) / xjac

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
                          
                            normal    = dot_product(grad_s,normal_direction) * grad_s      ! outward pointing normal
                            normal    = normal / norm2(normal)
                            direction = sign(1.d0,dot_product((/ps0_y,-ps0_x/),normal))
                            
                            grad_psi = sqrt(ps0_x**2 + ps0_y**2)
                            Btot     = sqrt(F0**2 + ps0_x**2 + ps0_y**2) / BigR
                            dl       = sqrt(R_t**2 + Z_t**2)

                            dl_dt    = (R_t*R_tt + Z_t*Z_tt) / dl

                            cs0    =   sqrt(gamma*T0)
                            cs0_T  =   0.5d0  * gamma    / cs0
                            cs0_TT = - 0.25d0 * gamma**2 / cs0**2 

                            bn     = dot_product( (/ps0_y,-ps0_x/), normal ) /  (BigR*Btot)  ! B\B7n/Btot
!                            bn     = ps0_t/(BigR*Btot*dl)
                            bn_t   = 1.d0 / (Btot*dl*BigR) * (ps0_tt - ps0_t * dl_dt /dl )

                            c_1 = vpar_smoothing_coef(1); c_2 = vpar_smoothing_coef(2); c_3 = vpar_smoothing_coef(3)

                            if (vpar_smoothing) then
                               if (c_2 .gt. 0d0) then
                                 factor    = 0.25d0 * ( 1.d0 + tanh( (abs(bn) - c_1) / c_2 ) )**2 - c_3
                                 factor_t  = 0.5d0  * ( 1.d0 + tanh( (abs(bn) - c_1) / c_2 ) )           & 
                                           * (bn_t * bn/abs(bn) /c_2) /(cosh( (abs(bn) - c_1) / c_2 ) )**2
                                else
                                  factor    = tanh(bn/c_1)
                                  factor_t  = bn_t /c_1 / cosh(bn/c_1)**2
                                  direction = 1.d0                            
                                endif                       
                            else
                              factor   = 1.d0
                              factor_t = 0.d0
                            endif

                            Hfact_t   = factor * R_t / BigR  + factor_t

 !                           if ((in .eq. 1) .and. (node_list%node(inode)%boundary .eq. 5)) then                              
 !                            write(111,'(6i8,A,3e14.6,A,20e14.6)') node_list%node(inode)%boundary, ielm, inode, iv, inode2, iv2, &
 !                                                              ' Boundary (t): ',Vpar0, -BigR**2 * u0_t/ps0_t, factor*direction*sqrt(GAMMA*T0)/Btot,&
 !                                                              ' error : ',Vpar0 - BigR**2 * u0_t/ps0_t - factor*direction*sqrt(GAMMA*T0)/Btot                                                               
 !                           endif

                            ku = 2
                            kv = 7
                            kT = 6
                            call boundary_conditions_add_one_entry(             &
                                 index_node, kv, in, index_node, kv, in,        &
                                 zbig,                                          &
                                 solve_only, gmres, index_min, index_max,       & 
                                 ijA_index, ijA_size, irn_jcn, irn, jcn, A_mat, i_tor_min, i_tor_max)

                            call boundary_conditions_add_one_entry(             &
                                 index_node, kv, in, index_node, kT, in,        &
                                 - zbig * factor / Btot * cs0_T * direction,    &
                                 solve_only, gmres, index_min, index_max,       & 
                                 ijA_index, ijA_size, irn_jcn, irn, jcn, A_mat, i_tor_min, i_tor_max)

                            call boundary_conditions_add_one_entry(             &
                                 index_node,  kv, in, index_node2, ku, in,      &
                                 - zbig * factor * BigR**2 * element_size_0 / ps0_t,   &
                                 solve_only, gmres, index_min, index_max,       & 
                                 ijA_index, ijA_size, irn_jcn, irn, jcn, A_mat, i_tor_min, i_tor_max)

                            if (in .eq. 1) then
                               call boundary_conditions_add_RHS(                      &
                                    index_node, kv, in,index_min, index_max, RHS_loc, &
                                    Zbig * ( - Vpar0 + factor*(BigR**2 * U0_t/ps0_t + cs0 * direction) / Btot), &
                                    i_tor_min, i_tor_max)
                            else
                               call boundary_conditions_add_RHS(                       &
                                    index_node, kv, in, index_min, index_max, RHS_loc, &
                                    0.d0,                                              &
                                    i_tor_min, i_tor_max)

                            endif
  
                            index_node  = node_list%node(inode)%index(1)
                            index_node2 = node_list%node(inode)%index(3)
                            kv = 7
                            kT = 6

                            call boundary_conditions_add_one_entry(               &
                                   index_node2, kv, in, index_node2, kv, in,      &
                                   zbig * element_size_0,                         &
                                   solve_only, gmres, index_min, index_max,       & 
                                   ijA_index, ijA_size, irn_jcn, irn, jcn, A_mat, i_tor_min, i_tor_max)

                            call boundary_conditions_add_one_entry(               &
                                   index_node2, kv, in, index_node2, kT, in,      &
                                   - zbig * element_size_0 * factor / Btot * cs0_T * direction,    &
                                   solve_only, gmres, index_min, index_max,       & 
                                   ijA_index, ijA_size, irn_jcn, irn, jcn, A_mat, i_tor_min, i_tor_max)

                            call boundary_conditions_add_one_entry(               &
                                   index_node2, kv, in, index_node,  kT, in,      &
                                   - zbig * factor  / Btot * cs0_TT * dT0_dt * direction  &
                                   - zbig * Hfact_t / Btot * cs0_T           * direction, & 
                                   solve_only, gmres, index_min, index_max,       & 
                                   ijA_index, ijA_size, irn_jcn, irn, jcn, A_mat, i_tor_min, i_tor_max)

                            if (in .eq. 1) then
                              call boundary_conditions_add_RHS(                    &
                                     index_node2, kv, in, index_min, index_max, RHS_loc,              &
                                     Zbig*(-dVpar0_dt + factor  / Btot * cs0_T * dT0_dt * direction   &
                                                      + Hfact_t / Btot * cs0            * direction), &
                                     i_tor_min, i_tor_max)
                            else
                               call boundary_conditions_add_RHS(                   &
                                     index_node2, kv, in, index_min, index_max, RHS_loc, &
                                     0.d0,                                         &
                                     i_tor_min, i_tor_max) 
                            endif

                       endif

                    endif



                    !------------------------------------ wall aligned with fluxsurface (in case of x-point grid)
                    if    ((node_list%node(inode)%boundary .eq.  2) &
                      .or. (node_list%node(inode)%boundary .eq. 12) &
                      .or. (node_list%node(inode)%boundary .eq.  3)) then

!====================================== begining RMPs at boundary ======================================================
!================================== type 2 - boundary: only depends on 't'
! ======================================================================================================================

                       if (RMP_on ) then
                          do n_rmp_harm=1, Number_RMP_harmonics !===========do RMP harmonics

                          if ((k.eq.1) .and. ((in.eq.RMP_har_cos_spectrum(n_rmp_harm)) .or. (in.eq.RMP_har_sin_spectrum(n_rmp_harm))) .and. (.not. freeboundary)) then
                             ! in .eq. RMP_har_cos  corresponds to cos(n_perturbation)
                             ! in .eq. RMP_har_sin   corresponds to sin(n_perturbation)



                             kp=1    ! variable psi
                             kv=1    ! equation for psi

                             index_node = node_list%node(inode)%index(1)  ! index in RHS (or matrix A not compressed)
                                                
                             Rnode     = node_list%node(inode)%x(1,1,1) 
                             dRnode_dt = node_list%node(inode)%x(1,3,1) 
                             Znode     = node_list%node(inode)%x(1,1,2) 
                             dZnode_dt = node_list%node(inode)%x(1,3,2) 
                          
                             if (in.eq.RMP_har_cos_spectrum(n_rmp_harm)) then
                                delta_psi_rmp = psi_RMP_cos1(node_list%node(inode)%boundary_index+N_rmp_har_block_size*(n_rmp_harm-1))
                                delta_psi_rmp_dR = dpsi_RMP_cos_dR1(node_list%node(inode)%boundary_index+N_rmp_har_block_size*(n_rmp_harm-1))
                                delta_psi_rmp_dZ = dpsi_RMP_cos_dZ1(node_list%node(inode)%boundary_index+N_rmp_har_block_size*(n_rmp_harm-1))

                                if (node_list%node(inode)%boundary_index == 1 ) then
                                   write (*,*) 'type2_bnd: my_id, psi_RMP_cos1, Rnode, Znode, in'
                                   write (*,*) my_id, delta_psi_rmp, Rnode, Znode,in
                                   write (*,*) 'delta_psi_rmp_dR, delta_psi_rmp_dZ'
                                   write (*,*) delta_psi_rmp_dR, delta_psi_rmp_dZ
                                endif
                             else 
                                delta_psi_rmp = psi_RMP_sin1(node_list%node(inode)%boundary_index+N_rmp_har_block_size*(n_rmp_harm-1))
                                delta_psi_rmp_dR = dpsi_RMP_sin_dR1(node_list%node(inode)%boundary_index+N_rmp_har_block_size*(n_rmp_harm-1))
                                delta_psi_rmp_dZ = dpsi_RMP_sin_dZ1(node_list%node(inode)%boundary_index+N_rmp_har_block_size*(n_rmp_harm-1))

                             endif

                             delta_psi_rmp_ds = delta_psi_rmp_dR * dRnode_ds + delta_psi_rmp_dZ * dZnode_ds
                             delta_psi_rmp_dt = delta_psi_rmp_dR * dRnode_dt + delta_psi_rmp_dZ * dZnode_dt
                             if (in.eq.RMP_har_cos_spectrum(n_rmp_harm)) then

                                if (node_list%node(inode)%boundary_index == 1 ) then
                                   write (*,*) 'delta_psi_rmp_dt', delta_psi_rmp_dt
                                   write (*,*) 'delta_psi_rmp_ds', delta_psi_rmp_ds
                                endif
                             endif

                             call boundary_conditions_add_one_entry( &
                                  index_node, kv, in,                &
                                  index_node, kp, in,                &
                                  zbig, solve_only, gmres,           &
                                  index_min, index_max,              & 
                                  ijA_index, ijA_size, irn_jcn,      & 
                                  irn, jcn, A_mat, i_tor_min, i_tor_max)

                             call boundary_conditions_add_RHS(  &
                                     index_node, kv, in,           &
                                     index_min, index_max,         &
                                     RHS_loc, ZBIG * delta_psi_rmp,&
                                     i_tor_min, i_tor_max)

                             index_node2 = node_list%node(inode)%index(3)
                             ! --- special case for grid with patches
                             if (node_list%node(inode)%boundary .eq. 12) then
                               index_node2 = node_list%node(inode)%index(2)
                               delta_psi_rmp_dt = delta_psi_rmp_ds
                             endif

                             call boundary_conditions_add_one_entry( &
                                  index_node2, kv, in,               &
                                  index_node2, kp, in,               &
                                  zbig, solve_only, gmres,           &
                                  index_min, index_max,              & 
                                  ijA_index, ijA_size, irn_jcn,      & 
                                  irn, jcn, A_mat, i_tor_min, i_tor_max)
                             
                             call boundary_conditions_add_RHS(          &
                                     index_node2, kv, in,               &
                                     index_min, index_max,              &
                                     RHS_loc, ZBIG * delta_psi_rmp_dt,  &
                                     i_tor_min, i_tor_max)
                          endif

                        enddo        !(end RMP harmonics)
                        endif        !(end RMPs on)  ==================================
!======================================= end RMPs ==================================

                       ! decides when the boundary conditions should be applied (for freeboundary and RMP cases)
                        if (       apply_psi_BC                   &
                              .or. apply_current_BC               &
                              .or. (( k /= 1 ) .and. ( k /= 3 ))  ) then

                          index_node = node_list%node(inode)%index(1)

                          call boundary_conditions_add_one_entry( &
                               index_node,  k,  in,               &
                               index_node,  k,  in,               &
                               zbig, solve_only, gmres,           &
                               index_min, index_max,              & 
                               ijA_index, ijA_size, irn_jcn,      & 
                               irn, jcn, A_mat, i_tor_min, i_tor_max)

                          index_node = node_list%node(inode)%index(3)

                          call boundary_conditions_add_one_entry( &
                               index_node,  k,  in,               &
                               index_node,  k,  in,               &
                               zbig, solve_only, gmres,           &
                               index_min, index_max,              & 
                               ijA_index, ijA_size, irn_jcn,      & 
                               irn, jcn, A_mat, i_tor_min, i_tor_max)

                       endif

                    endif

                    !------------------------------------ Special corners (only for grid with patches)
                    if    ((node_list%node(inode)%boundary .eq. 21) &
                      .or. (node_list%node(inode)%boundary .eq. 20)) then

!====================================== begining RMPs at boundary ======================================================
!================================== type 20-21 - boundary: corners apply on 's' and 't'
! ======================================================================================================================
                       
                       if (RMP_on ) then
                          do n_rmp_harm=1, Number_RMP_harmonics !===========do RMP harmonics

                          if ((k.eq.1) .and. ((in.eq.RMP_har_cos_spectrum(n_rmp_harm)) .or. (in.eq.RMP_har_sin_spectrum(n_rmp_harm))) .and. (.not. freeboundary)) then
                             ! in .eq. RMP_har_cos  corresponds to cos(n_perturbation)
                             ! in .eq. RMP_har_sin   corresponds to sin(n_perturbation)

                             kp=1    ! variable psi
                             kv=1    ! equation for psi

                             index_node = node_list%node(inode)%index(1)  ! index in RHS (or matrix A not compressed)
                                                
                             Rnode     = node_list%node(inode)%x(1,1,1) 
                             dRnode_dt = node_list%node(inode)%x(1,3,1) 
                             Znode     = node_list%node(inode)%x(1,1,2) 
                             dZnode_dt = node_list%node(inode)%x(1,3,2) 
                          
                             if (in.eq.RMP_har_cos_spectrum(n_rmp_harm)) then
                                delta_psi_rmp = psi_RMP_cos1(node_list%node(inode)%boundary_index+N_rmp_har_block_size*(n_rmp_harm-1))
                                delta_psi_rmp_dR = dpsi_RMP_cos_dR1(node_list%node(inode)%boundary_index+N_rmp_har_block_size*(n_rmp_harm-1))
                                delta_psi_rmp_dZ = dpsi_RMP_cos_dZ1(node_list%node(inode)%boundary_index+N_rmp_har_block_size*(n_rmp_harm-1))

                                if (node_list%node(inode)%boundary_index == 1 ) then
                                   write (*,*) 'type2_bnd: my_id, psi_RMP_cos1, Rnode, Znode, in'
                                   write (*,*) my_id, delta_psi_rmp, Rnode, Znode,in
                                   write (*,*) 'delta_psi_rmp_dR, delta_psi_rmp_dZ'
                                   write (*,*) delta_psi_rmp_dR, delta_psi_rmp_dZ
                                endif
                             else 
                                delta_psi_rmp = psi_RMP_sin1(node_list%node(inode)%boundary_index+N_rmp_har_block_size*(n_rmp_harm-1))
                                delta_psi_rmp_dR = dpsi_RMP_sin_dR1(node_list%node(inode)%boundary_index+N_rmp_har_block_size*(n_rmp_harm-1))
                                delta_psi_rmp_dZ = dpsi_RMP_sin_dZ1(node_list%node(inode)%boundary_index+N_rmp_har_block_size*(n_rmp_harm-1))

                             endif

                             delta_psi_rmp_ds = delta_psi_rmp_dR * dRnode_ds + delta_psi_rmp_dZ * dZnode_ds
                             delta_psi_rmp_dt = delta_psi_rmp_dR * dRnode_dt + delta_psi_rmp_dZ * dZnode_dt
                             if (in.eq.RMP_har_cos_spectrum(n_rmp_harm)) then

                                if (node_list%node(inode)%boundary_index == 1 ) then
                                   write (*,*) 'delta_psi_rmp_dt', delta_psi_rmp_dt
                                   write (*,*) 'delta_psi_rmp_ds', delta_psi_rmp_ds
                                endif
                             endif

                             call boundary_conditions_add_one_entry(   &
                                  index_node, kv, in,                  &
                                  index_node, kp, in,                  &
                                  zbig, solve_only, gmres,             &
                                  index_min, index_max,                & 
                                  ijA_index, ijA_size, irn_jcn,        & 
                                  irn, jcn, A_mat, i_tor_min, i_tor_max)

                                call boundary_conditions_add_RHS(      &
                                     index_node, kv, in,               &
                                     index_min, index_max,             &
                                     RHS_loc, ZBIG * delta_psi_rmp,    &
                                     i_tor_min, i_tor_max)

                             index_node2 = node_list%node(inode)%index(2)

                             call boundary_conditions_add_one_entry( &
                                  index_node2, kv, in,               &
                                  index_node2, kp, in,               &
                                  zbig, solve_only, gmres,           &
                                  index_min, index_max,              & 
                                  ijA_index, ijA_size, irn_jcn,      & 
                                  irn, jcn, A_mat, i_tor_min, i_tor_max)
                                call boundary_conditions_add_RHS(       &
                                     index_node2, kv, in,               &
                                     index_min, index_max,              &
                                     RHS_loc, ZBIG * delta_psi_rmp_ds,  &
                                     i_tor_min, i_tor_max)

                             index_node2 = node_list%node(inode)%index(3)

                             call boundary_conditions_add_one_entry(   &
                                  index_node2, kv, in,                 &
                                  index_node2, kp, in,                 &
                                  zbig, solve_only, gmres,             &
                                  index_min, index_max,                & 
                                  ijA_index, ijA_size, irn_jcn,        & 
                                  irn, jcn, A_mat, i_tor_min, i_tor_max)
                                call boundary_conditions_add_RHS(      &
                                     index_node2, kv, in,              &
                                     index_min, index_max,             &
                                     RHS_loc, ZBIG * delta_psi_rmp_dt, &
                                     i_tor_min, i_tor_max)

                          endif
                        enddo        !(end RMP harmonics)
                        endif        !(end RMPs on)  ==================================
                        !======================================= end RMPs ==================================

                       ! decides when the boundary conditions should be applied (for freeboundary and RMP cases)
                        if (       apply_psi_BC                   &
                              .or. apply_current_BC               &
                              .or. (( k /= 1 ) .and. ( k /= 3 ))  ) then

                          index_node = node_list%node(inode)%index(1)

                          call boundary_conditions_add_one_entry(   &
                               index_node,  k,  in,                 &
                               index_node,  k,  in,                 &
                               zbig, solve_only, gmres,             &
                               index_min, index_max,                & 
                               ijA_index, ijA_size, irn_jcn,        & 
                               irn, jcn, A_mat, i_tor_min, i_tor_max)

                          index_node = node_list%node(inode)%index(2)

                          call boundary_conditions_add_one_entry( &
                               index_node,  k,  in,               &
                               index_node,  k,  in,               &
                               zbig, solve_only, gmres,           &
                               index_min, index_max,              & 
                               ijA_index, ijA_size, irn_jcn,      & 
                               irn, jcn, A_mat, i_tor_min, i_tor_max)

                          index_node = node_list%node(inode)%index(3)

                          call boundary_conditions_add_one_entry(   &
                               index_node,  k,  in,                 &
                               index_node,  k,  in,                 &
                               zbig, solve_only, gmres,             &
                               index_min, index_max,                & 
                               ijA_index, ijA_size, irn_jcn,        & 
                               irn, jcn, A_mat, i_tor_min, i_tor_max)

                       endif

                    endif

                 enddo  ! ================================do variables

!              enddo !===========do RMP harmonics , by defalt=1

              enddo !========================do n_tor

           endif !==================if boundary nodes

     enddo  !==========================do vertex

 
    enddo !===============================do elements

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
