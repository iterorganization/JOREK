module mod_boundary_conditions

  !*******************************************************************************
  !************ Define global variables for all internal routines ****************
  !*******************************************************************************
  ! --- ZBIG parameter to make equations "more important" than element_matrix equations
  real*8  		        :: zbig = 1.d10
  ! --- R,Z variables
  real*8			:: R, R_s, R_t, R_inside
  real*8			:: Z, Z_s, Z_t, Z_inside
  real*8			:: alpha
  real*8			:: xjac, zbig_backup
  ! --- Variable numbers
  integer, parameter		:: k_psi  = 1
  integer, parameter		:: k_u    = 2
  integer, parameter		:: k_Vpar = 7
  integer, parameter		:: k_Ti   = 6
  integer, parameter		:: k_rho  = 5
  integer, parameter		:: k_Te   = 8
  ! --- Variables
  real*8			:: ps0,   ps0_s,   ps0_t, ps0_x, ps0_y, grad_psi, Btot
  real*8			:: u0,    u0_s,    u0_t,  u0_x,  u0_y
  real*8			:: Vpar0, Vpar0_s, Vpar0_t
  real*8			:: rho0,  rho0_s,  rho0_t
  real*8			:: Ti0,   Ti0_s,   Ti0_t
  real*8			:: Pi0,   Pi0_s,   Pi0_t
  real*8			:: tau_IC, rho_norm
  real*8			:: lhs_tmp, rhs_tmp
  ! --- Direction of Vpar on target
  real*8			:: direction


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
  subroutine boundary_conditions(my_id, node_list, element_list, bnd_node_list,           &
                                 local_elms, n_local_elms, index_min, index_max, rhs_loc, &
                                 xpoint2, xcase2,                                         &
                                 R_axis, Z_axis, psi_axis,                                &
                                 psi_bnd, R_xpoint, Z_xpoint, psi_xpoint, a_mat)

    use constants
    use data_structure
    use vacuum, ONLY: is_freebound
    use phys_module, only: F0, GAMMA, freeboundary, tokamak_device, U_sheath,               &
                           RMP_on, psi_RMP_cos, dpsi_RMP_cos_dR, dpsi_RMP_cos_dZ,           &
                           psi_RMP_sin, dpsi_RMP_sin_dR, dpsi_RMP_sin_dZ,                   &
                           t_now, RMP_growth_rate, RMP_ramp_up_time,                        &
                           RMP_start_time, tstep, RMP_har_cos, RMP_har_sin,                 &
                           grid_to_wall, n_wall_blocks, keep_n0_const
    USE tr_module
    use mpi_mod
    use mod_integer_types

    implicit none
    !include 'mpif.h'

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

    ! --- Internal parameters
    real*8  :: mach1, dmach1, d2mach1_dTi, d2mach1_dTe, mach_u, dmach_u, dmach_rho
    integer :: i, i_tor, iv, inode, k_var, side
    integer :: ielm
    integer :: ierr
    logical :: apply_dirichlet, apply_on_psi, apply_on_current, on_private, on_inner, on_inner_or_private

    ! --- RMP parameters
    real*8, allocatable :: psi_RMP_cos1(:),dpsi_RMP_cos_dR1(:),dpsi_RMP_cos_dZ1(:)
    real*8, allocatable :: psi_RMP_sin1(:),dpsi_RMP_sin_dR1(:),dpsi_RMP_sin_dZ1(:)
    real*8              :: establish_RMP
    real*8              :: delta_psi_rmp, delta_psi_rmp_dR, delta_psi_rmp_dZ, delta_psi_rmp_ds, delta_psi_rmp_dt, psi_test, sigmo_fonc
    integer             :: ilarge_vp, ilarge_vp2
    integer             :: j, err, itest

    zbig_backup = zbig
    ! -------------------------
    ! --- Retrieve RMP profiles
    if ( RMP_on .and. (n_tor .ge. 3) ) then
      call tr_allocate(psi_RMP_cos1,    1, bnd_node_list%n_bnd_nodes,"psi_RMP_cos1",    CAT_UNKNOWN)
      call tr_allocate(dpsi_RMP_cos_dR1,1, bnd_node_list%n_bnd_nodes,"dpsi_RMP_cos_dR1",CAT_UNKNOWN)
      call tr_allocate(dpsi_RMP_cos_dZ1,1, bnd_node_list%n_bnd_nodes,"dpsi_RMP_cos_dZ1",CAT_UNKNOWN)
      call tr_allocate(psi_RMP_sin1,    1, bnd_node_list%n_bnd_nodes,"psi_RMP_sin1",    CAT_UNKNOWN)
      call tr_allocate(dpsi_RMP_sin_dR1,1, bnd_node_list%n_bnd_nodes,"dpsi_RMP_sin_dR1",CAT_UNKNOWN)
      call tr_allocate(dpsi_RMP_sin_dZ1,1, bnd_node_list%n_bnd_nodes,"dpsi_RMP_sin_dZ1",CAT_UNKNOWN)

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
        psi_RMP_cos1(j)     =  psi_RMP_cos(j)	 * establish_RMP
        dpsi_RMP_cos_dR1(j) = dpsi_RMP_cos_dR(j) * establish_RMP
        dpsi_RMP_cos_dZ1(j) = dpsi_RMP_cos_dZ(j) * establish_RMP
        psi_RMP_sin1(j)     =  psi_RMP_sin(j)	 * establish_RMP
        dpsi_RMP_sin_dR1(j) = dpsi_RMP_sin_dR(j) * establish_RMP
        dpsi_RMP_sin_dZ1(j) = dpsi_RMP_sin_dZ(j) * establish_RMP
      end do
     
      !if (my_id == 0) then
      !  write (*,*) 'psi_RMP_cos1(1) and derivatives after multiplication in boundary conditions'
      !  write (*,*) psi_RMP_cos1(1), dpsi_RMP_cos_dR1(1), dpsi_RMP_cos_dZ1(1)
      !  write (*,*) 'establish_RMP', establish_RMP
      !endif
    
    endif
    ! --- Retrieve RMP profiles (END)
    ! -------------------------------
    
      ! --- Loop on each element
      do i=1, n_local_elms
        ielm = local_elms(i)

        ! --- Take each node of element
        do iv=1, n_vertex_max
          inode = element_list%element(ielm)%vertex(iv)

          ! --- We only care about boundary elements
          if (node_list%node(inode)%boundary .ne. 0) then

            call construct_variables(node_list%node(inode), R_axis, Z_axis, R_xpoint, Z_xpoint, psi_bnd)
            
            do i_tor=a_mat%i_tor_min, a_mat%i_tor_max
              if (keep_n0_const  .and.  i_tor .eq. 1 ) then
                 zbig = 1.d15
               else
                 zbig = zbig_backup
               endif
              
              do k_var=1, n_var

                ! --------------------------------------------------------------------------------------------------------------
                ! ------------------------------------ the targets (in case of x-point grid) -----------------------------------
                ! --------------------------------------------------------------------------------------------------------------
                if    ((node_list%node(inode)%boundary .eq.  1) &
                  .or. (node_list%node(inode)%boundary .eq. 11) &
                  .or. (node_list%node(inode)%boundary .eq.  9) &
                  .or. (node_list%node(inode)%boundary .eq. 19) &
                  .or. (node_list%node(inode)%boundary .eq.  3) &
                  .or. (node_list%node(inode)%boundary .eq.  4)) then
                      
                  ! --- Which side is this? 2 => d/ds, 3 => d/dt
                  side = 2

                  ! --- Field direction
                  if ( (grid_to_wall) .and. (n_wall_blocks .gt. 0) ) then
                    direction = 1
                    if (node_list%node(inode)%boundary .eq. 11) direction = -direction
                    if (node_list%node(inode)%boundary .eq. 19) direction = -direction
                  endif

                  ! ---------------------------------------------
                  ! --- Apply RMP on target (only depends on 's')

                  if (      RMP_on                                                      &
                      .and. (k_var .eq. 1)                                              &
                      .and. ((i_tor.eq.RMP_har_cos) .or. (i_tor.eq.RMP_har_sin))        &
                      .and. (.not. freeboundary)                                        ) then
                                         
                      call apply_RMP_BCs(rhs_loc, node_list%node(inode), side, i_tor,           &
                                         psi_RMP_cos1, dpsi_RMP_cos_dR1, dpsi_RMP_cos_dZ1,      &
                                         psi_RMP_sin1, dpsi_RMP_sin_dR1, dpsi_RMP_sin_dZ1,      &
                                         index_min, index_max, a_mat)

                  endif
                  
                  ! -----------------------------------------------
                  ! --- Dirichlet BCs (or Neumann if commented out)
                  apply_dirichlet = .false.
                  ! --- Determine if we need to apply condition on psi (we don't want to overwrite RMPs)
                  apply_on_psi = .false.
                  if ((k_var .eq. 1) .and. (.not. is_freebound(i_tor,k_var))) then
                    if                        (i_tor .eq. 1)             apply_on_psi = .true.
                    if ( (.not. RMP_on) .and. (i_tor .ge. 2 )          ) apply_on_psi = .true.
                    if ( (RMP_on)       .and. (i_tor .lt. RMP_har_cos) ) apply_on_psi = .true.
                    if ( (RMP_on)       .and. (i_tor .gt. RMP_har_sin) ) apply_on_psi = .true.
                  endif
                  
                  apply_on_current = .false.
                  if ((k_var .eq. 3) .and. (.not. is_freebound(i_tor,k_var))) apply_on_current = .true.
                  
                  ! --- Apply conditions to which variables?
                  if (                                                  &
                           apply_on_psi                                 &
                      .or. apply_on_current                             &
                      .or.( (k_var .eq. 2) .and. (.not. U_sheath) )     &
                      .or.  (k_var .eq. 4)                              &
                      !.or.  (k_var .eq. 7)                              &
                      ) apply_dirichlet = .true.


                  ! --- Apply Dirichlet if required
                  if (apply_dirichlet) then
                    call apply_Dirichlet_BCs(node_list%node(inode), side, k_var,i_tor, index_min, index_max, a_mat)
                  endif

                  ! --------------
                  ! --- Mach-1 BCs
                  if (k_var .eq. k_Vpar) then
                    call apply_Mach1_BCs(rhs_loc, node_list%node(inode), side, i_tor, index_min, index_max, a_mat)
                    if (U_sheath) call apply_U_sheath (rhs_loc, node_list%node(inode), side, i_tor, index_min,index_max, a_mat)
                  endif

                endif

                
                ! --------------------------------------------------------------------------------------------------------------
                ! ------------------------- the non-targets open field-lines (for grid_xpoint_wall) ----------------------------
                ! --------------------------------------------------------------------------------------------------------------
                if    ((node_list%node(inode)%boundary .eq.  5) &
                  .or. (node_list%node(inode)%boundary .eq. 15) &
                  .or. (node_list%node(inode)%boundary .eq.  9) &
                  .or. (node_list%node(inode)%boundary .eq. 19)) then

                  ! --- Which side is this? 2 => d/ds, 3 => d/dt
                  side = 3

                  ! --- Field direction
                  if ( (grid_to_wall) .and. (n_wall_blocks .gt. 0) ) then
                    direction = 1
                    if (node_list%node(inode)%boundary .eq. 15) direction = -direction
                    if (node_list%node(inode)%boundary .eq. 19) direction = -direction
                  endif                 

                  ! ---------------------------------------------
                  ! --- Apply RMP on target (only depends on 's')
                  if (      RMP_on                                                      &
                      .and. (k_var .eq. 1)                                              &
                      .and. ((i_tor.eq.RMP_har_cos) .or. (i_tor.eq.RMP_har_sin))        &
                      .and. (.not. freeboundary)                                        ) then
                                         
                      call apply_RMP_BCs(rhs_loc, node_list%node(inode), side, i_tor,           &
                                         psi_RMP_cos1, dpsi_RMP_cos_dR1, dpsi_RMP_cos_dZ1,      &
                                         psi_RMP_sin1, dpsi_RMP_sin_dR1, dpsi_RMP_sin_dZ1,      &
                                         index_min, index_max, a_mat)
                  endif
                  
                  ! -----------------------------------------------
                  ! --- Dirichlet BCs (or Neumann if commented out)
                  apply_dirichlet = .false.
                  ! --- Determine if we need to apply condition on psi (we don't want to overwrite RMPs)
                  apply_on_psi = .false.
                  if ((k_var .eq. 1) .and. (.not. is_freebound(i_tor,k_var))) then
                    if                        (i_tor .eq. 1)             apply_on_psi = .true.
                    if ( (.not. RMP_on) .and. (i_tor .ge. 2 )          ) apply_on_psi = .true.
                    if ( (RMP_on)       .and. (i_tor .lt. RMP_har_cos) ) apply_on_psi = .true.
                    if ( (RMP_on)       .and. (i_tor .gt. RMP_har_sin) ) apply_on_psi = .true.
                  endif
                  
                  apply_on_current = .false.
                  if ((k_var .eq. 3) .and. (.not. is_freebound(i_tor,k_var))) apply_on_current = .true.
                  
                  ! --- Apply conditions to which variables?
                  if (                                                  &
                           apply_on_psi                                 &
                      .or. apply_on_current                             &
                      .or.( (k_var .eq. 2) .and. (.not. U_sheath) )     &
                      .or.  (k_var .eq. 4)                              &
                      ) apply_dirichlet = .true.


                  ! --- Apply Dirichlet if required
                  if (apply_dirichlet) then
                    call apply_Dirichlet_BCs(node_list%node(inode), side, k_var,i_tor, index_min, index_max, a_mat)
                  endif

                  ! --------------
                  ! --- Mach-1 BCs
                  if (k_var .eq. k_Vpar) then
                    call apply_Mach1_BCs(rhs_loc, node_list%node(inode), side, i_tor, index_min, index_max, a_mat)
                    if (U_sheath) call apply_U_sheath (rhs_loc, node_list%node(inode), side, i_tor, index_min,index_max, a_mat)
                  endif

                endif

                
                ! ----------------------------------------------------------------------------------------------------
                ! ------------------------------------ the flux-surface boundaries -----------------------------------
                ! ----------------------------------------------------------------------------------------------------
                if    ((node_list%node(inode)%boundary .eq.  2) &
                  .or. (node_list%node(inode)%boundary .eq. 12) &
                  .or. (node_list%node(inode)%boundary .eq.  3)) then

                  ! --- Which side is this? 2 => d/ds, 3 => d/dt
                  side = 3
                  if (node_list%node(inode)%boundary .eq. 12) side = 2
                  
                  ! ---------------------------------------------
                  ! --- Apply RMP on target (only depends on 't')

                  if (      RMP_on                                                      &
                      .and. (k_var .eq. 1)                                              &
                      .and. ((i_tor.eq.RMP_har_cos) .or. (i_tor.eq.RMP_har_sin))        &
                      .and. (.not. freeboundary)                                        ) then
                                         
                      call apply_RMP_BCs(rhs_loc, node_list%node(inode), side, i_tor,           &
                                         psi_RMP_cos1, dpsi_RMP_cos_dR1, dpsi_RMP_cos_dZ1,      &
                                         psi_RMP_sin1, dpsi_RMP_sin_dR1, dpsi_RMP_sin_dZ1,      &
                                         index_min, index_max, a_mat)
                  endif
                  
                  ! -----------------------------------------------
                  ! --- Dirichlet BCs (or Neumann if commented out)
                  apply_dirichlet = .false.
                  
                  ! --- Determine if we are on the private or the inner boundary
                  on_private            = .false.
                  on_inner              = .false.
                  on_inner_or_private   = .false.
                  if ((xcase2 .ne. DOUBLE_NULL) .and. (ps0 .lt. psi_bnd)) then
                    on_inner_or_private = .true.
                    on_private          = .true.
                  endif
                  if  (xcase2 .eq. DOUBLE_NULL) then
                    if ( (Z .lt. (Z_xpoint(1)+Z_xpoint(2))/2.d0) .and. (ps0 .lt. psi_xpoint(1)) )                    on_private = .true.
                    if ( (Z .gt. (Z_xpoint(1)+Z_xpoint(2))/2.d0) .and. (ps0 .lt. psi_xpoint(2)) )                    on_private = .true.
                    if ( (R .lt. (R_xpoint(1)+R_xpoint(2))/2.d0) .and. (ps0 .gt. max(psi_xpoint(2),psi_xpoint(2))) ) on_inner   = .true.
                  endif
                  if (on_private .or. on_inner) on_inner_or_private = .true.
                  
                  ! --- Determine if we need to apply condition on psi (we don't want to overwrite RMPs)
                  apply_on_psi = .false.
                  if ((k_var .eq. 1).and.(.not. is_freebound(i_tor,k_var))) then
                    if                        (i_tor .eq. 1)             apply_on_psi = .true.
                    if ( (.not. RMP_on) .and. (i_tor .ge. 2)           ) apply_on_psi = .true.
                    if ( (RMP_on)         .and. (i_tor .lt. RMP_har_cos) ) apply_on_psi = .true.
                    if ( (RMP_on)         .and. (i_tor .gt. RMP_har_sin) ) apply_on_psi = .true.
                  endif
                  
                  apply_on_current = .false.
                  if ((k_var .eq. 3) .and. (.not. is_freebound(i_tor,k_var))) apply_on_current = .true.
                  
                  ! Apply conditions to which variables and where?
                  if (                                                                  &
                            (apply_on_psi)                                              &
                      .or.  (apply_on_current)                                          &
                      .or.( (k_var .eq. 2) .and. (.not. U_sheath) )                     &
                      .or.  (k_var .eq. 4)                                              &
                      .or.  (k_var .eq. 5)                                              &
                      .or.  (k_var .eq. 6)                                              &
                      !.or.( (k_var .eq. 5) .and. (on_private) )                                & 
                      !.or.( (k_var .eq. 6) .and. (on_private) )                                & 
                      .or.  (k_var .eq. 7)                                              &
                      ) apply_dirichlet = .true.

                  if (apply_dirichlet) then
                    call apply_Dirichlet_BCs(node_list%node(inode), side, k_var,i_tor, index_min, index_max, a_mat)
                    if (U_sheath) call apply_U_sheath (rhs_loc, node_list%node(inode), side, i_tor, index_min,index_max, a_mat)
                  endif

                endif

                ! ----------------------------------------------------------------------------------------------------
                ! ------------------------------------ the special corners -------------------------------------------
                ! ----------------------------------------------------------------------------------------------------
                if (node_list%node(inode)%boundary .eq. 21) then
                  if ( (k_var .eq. k_Ti) .or. (k_var .eq. k_Vpar) ) then
                    side = 2
                    call apply_Dirichlet_BCs(node_list%node(inode), side, k_var,i_tor, index_min, index_max, a_mat)
                    side = 3
                    call apply_Dirichlet_BCs(node_list%node(inode), side, k_var,i_tor, index_min, index_max, a_mat)
                  endif
                endif
                if (node_list%node(inode)%boundary .eq. 20) then
                  side = 2
                  call apply_Dirichlet_BCs(node_list%node(inode), side, k_var,i_tor, index_min, index_max, a_mat)
                  side = 3
                  call apply_Dirichlet_BCs(node_list%node(inode), side, k_var,i_tor, index_min, index_max, a_mat)
                endif
              enddo

            enddo
          endif
        enddo
      enddo
    return
  end subroutine boundary_conditions
  
  
  
  
  
  
  
  
  
  
  !******************************************************************************
  !******************************************************************************
  !************ Routine to construct variables once for all routines ************
  !******************************************************************************
  !******************************************************************************
  subroutine construct_variables(node, R_axis, Z_axis, R_xpoint, Z_xpoint, psi_bnd)

    use constants    
    use data_structure
    use phys_module, only: F0, FF_0, xpoint, xcase, tokamak_device,    &
                           central_mass, mass_proton, central_density, &
                           mu_zero, tauIC, renormalise, grid_to_wall, n_wall_blocks
    use corr_neg
    
    implicit none
    
    ! --- Routine variables
    type (type_node),   intent(in)    :: node
    real*8,             intent(in)    :: R_axis
    real*8,             intent(in)    :: Z_axis
    real*8,             intent(in)    :: R_xpoint(2)
    real*8,             intent(in)    :: Z_xpoint(2)
    real*8,             intent(in)    :: psi_bnd
  
    ! --- Define (R,Z) coords and Jacobian
    R         = node%x(1,1,1)
    R_s       = node%x(1,2,1)
    R_t       = node%x(1,3,1)
    Z         = node%x(1,1,2)
    Z_s       = node%x(1,2,2)
    Z_t       = node%x(1,3,2)
    xjac      = R_s*Z_t - R_t*Z_s
    
    ! --- Define psi variables
    ps0       = node%values(1,1,1)
    ps0_s     = node%values(1,2,1)
    ps0_t     = node%values(1,3,1)
    ps0_x     = (  Z_t*ps0_s - Z_s*ps0_t) / xjac
    ps0_y     = (- R_t*ps0_s + R_s*ps0_t) / xjac
    grad_psi  = sqrt(        ps0_x**2 + ps0_y**2)
    Btot      = sqrt(F0**2 + ps0_x**2 + ps0_y**2) / R
    
    ! --- Define U variables
    U0        = node%values(1,1,2)
    U0_s      = node%values(1,2,2)
    U0_t      = node%values(1,3,2)
    u0_x      = (  Z_t*u0_s - Z_s*u0_t) / xjac
    u0_y      = (- R_t*u0_s + R_s*u0_t) / xjac

    ! --- Define Ti variables
    rho0      = node%values(1,1,5)
    rho0_s    = node%values(1,2,5)
    rho0_t    = node%values(1,3,5)
    
    ! --- Define Ti variables
    Ti0       = corr_neg_temp(node%values(1,1,6))
    Ti0_s     = node%values(1,2,6)
    Ti0_t     = node%values(1,3,6)
    
    ! --- Define Ti variables
    Pi0       = rho0 * Ti0
    Pi0_s     = rho0_s * Ti0 + rho0 * Ti0_s
    Pi0_t     = rho0_t * Ti0 + rho0 * Ti0_t
    
    ! --- Define Vpar variables
    Vpar0     = node%values(1,1,k_Vpar)
    Vpar0_s   = node%values(1,2,k_Vpar)
    Vpar0_t   = node%values(1,3,k_Vpar)

    ! --- Define direction of Vpar on target. Careful, using ps0_x/abs(ps0_x) can be treacherous.
    if ( (.not. grid_to_wall) .or. (n_wall_blocks .eq. 0) ) then
      if (tokamak_device(1:4) .eq. 'MAST') then
        if ( (R .gt. (R_xpoint(1)+R_xpoint(2))/2.d0) ) then
          direction = 1.d0
        else
          direction = -1.d0
        endif
      else
        !direction = + ps0_x / abs(ps0_x)
        alpha = (Z_axis - Z_xpoint(1))/(R_axis - R_xpoint(1))
        R_inside = alpha*(Z-Z_xpoint(1)) + R + alpha**2 * R_xpoint(1)
        R_inside = R_inside / (1.d0 + alpha**2)
        Z_inside = alpha * (R_inside - R_xpoint(1)) + Z_xpoint(1)
        R_inside = min(max(R_inside,R_xpoint(1)),R_axis)
        Z_inside = min(max(Z_inside,Z_xpoint(1)),Z_axis)
        direction = ps0_s * ( (R-R_inside)*Z_s - (Z-Z_inside)*R_s )
        direction = direction / abs(direction)
      endif
      if (xcase .eq. UPPER_XPOINT) then
        direction = -direction
      else if ((xcase .eq. DOUBLE_NULL).and.(Z .gt. Z_axis +0.1) .and. ( R .gt.R_xpoint(2))) then
        direction = -1.
      else if ((xcase .eq. DOUBLE_NULL) .and. (Z .gt. Z_axis +0.1) .and. (R .lt. R_xpoint(2))) then
        direction = +1.
      end if
    endif

    ! --- Diamagnetic term
    tau_IC = tauIC
    
    ! --- Renormalise MHD parameters?
    if (central_density .gt. 1.d10) then
      rho_norm = central_density	 * central_mass * mass_proton
    else
      rho_norm = central_density * 1.d20 * central_mass * mass_proton
    endif
    if (renormalise) tau_IC = tau_IC / sqrt(mu_zero * rho_norm)

    return
  end subroutine construct_variables
  
  
  




  
 
  
  !******************************************************************************
  !******************************************************************************
  !********* Routine to apply RMP perturbation on boundary conditions ***********
  !******************************************************************************
  !******************************************************************************
  subroutine apply_RMP_BCs(RHS_loc, node, side, i_tor,                          &
                           psi_RMP_cos1, dpsi_RMP_cos_dR1, dpsi_RMP_cos_dZ1,    &
                           psi_RMP_sin1, dpsi_RMP_sin_dR1, dpsi_RMP_sin_dZ1,    &
                           index_min, index_max, a_mat)

  
    use mod_parameters
    use data_structure
    use phys_module, only: RMP_har_cos, RMP_har_sin
    use mod_assembly, only : boundary_conditions_add_one_entry, boundary_conditions_add_RHS
    use mod_integer_types
    
    implicit none
    
    ! --- Routine variables
    real*8,                             intent(inout) :: RHS_loc(*)
    type (type_node),                   intent(in)    :: node
    integer,                            intent(in)    :: side ! == 2 for d/ds, == 3 for d/dt
    integer,                            intent(in)    :: i_tor
    real*8,                             intent(in)    :: psi_RMP_cos1(*), dpsi_RMP_cos_dR1(*), dpsi_RMP_cos_dZ1(*)
    real*8,                             intent(in)    :: psi_RMP_sin1(*), dpsi_RMP_sin_dR1(*), dpsi_RMP_sin_dZ1(*)
    integer,                            intent(in)    :: index_min, index_max
    type(type_SP_MATRIX)                              :: a_mat

    ! --- Internal variables
    integer                                           :: index_node,   index_node2
    real*8                                            :: delta_psi_rmp, delta_psi_rmp_dR, delta_psi_rmp_dZ, delta_psi_rmp_dl
    
    ! --- Get psi perturbation and its derivatives
    if (i_tor.eq.RMP_har_cos) then
      delta_psi_rmp    =  psi_RMP_cos1   (node%boundary_index)
      delta_psi_rmp_dR = dpsi_RMP_cos_dR1(node%boundary_index)
      delta_psi_rmp_dZ = dpsi_RMP_cos_dZ1(node%boundary_index)
    else 
      delta_psi_rmp    =  psi_RMP_sin1   (node%boundary_index)
      delta_psi_rmp_dR = dpsi_RMP_sin_dR1(node%boundary_index)
      delta_psi_rmp_dZ = dpsi_RMP_sin_dZ1(node%boundary_index)
    endif
    if (side .eq. 2) then
      delta_psi_rmp_dl = delta_psi_rmp_dR*R_s + delta_psi_rmp_dZ*Z_s
    else
      delta_psi_rmp_dl = delta_psi_rmp_dR*R_t + delta_psi_rmp_dZ*Z_t
    endif
    
    ! --- Get nodes index
    index_node  = node%index(1)
    index_node2 = node%index(side)
    
    ! --- Condition on nodes
    lhs_tmp = ZBIG
    rhs_tmp = ZBIG * delta_psi_rmp
    
    call boundary_conditions_add_one_entry( &
         index_node, k_psi, i_tor,          &
         index_node, k_psi, i_tor,          &
         lhs_tmp, index_min, index_max, a_mat)
      call boundary_conditions_add_RHS(     &
           index_node, k_psi, i_tor,        &
           index_min, index_max,            &
           RHS_loc, rhs_tmp,                &
           a_mat%i_tor_min, a_mat%i_tor_max)
    
    ! --- Condition between nodes (d/ds or d/dt)
    lhs_tmp = ZBIG
    rhs_tmp = ZBIG * delta_psi_rmp_dl
    
    call boundary_conditions_add_one_entry( &
         index_node2, k_psi, i_tor,         &
         index_node2, k_psi, i_tor,         &
         lhs_tmp, index_min, index_max, a_mat)
      call boundary_conditions_add_RHS(     &
           index_node2, k_psi, i_tor,       &
           index_min, index_max,            &
           RHS_loc, rhs_tmp,                &
           a_mat%i_tor_min, a_mat%i_tor_max)
  
    return
  end subroutine apply_RMP_BCs
  
  
  
  
  
  !******************************************************************************
  !******************************************************************************
  !***************** Routine to apply Dirichlet boundary conditions *************
  !******************************************************************************
  !******************************************************************************
  subroutine apply_Dirichlet_BCs(node, side, k_var,i_tor, index_min,index_max, a_mat)
  
    use mod_parameters
    use data_structure
    use phys_module, only: RMP_har_cos, RMP_har_sin
    use mod_assembly, only : boundary_conditions_add_one_entry, boundary_conditions_add_RHS
    use mod_integer_types
    
    implicit none
    
    ! --- Routine variables
    type (type_node),                   intent(in)    :: node
    integer,                            intent(in)    :: side ! == 2 for d/ds, == 3 for d/dt
    integer,                            intent(in)    :: k_var
    integer,                            intent(in)    :: i_tor
    integer,                            intent(in)    :: index_min, index_max
    type(type_SP_MATRIX)                              :: a_mat
    ! --- Internal variables
    integer                                           :: index_node,   index_node2
    
    ! --- Get nodes index
    index_node  = node%index(1)
    index_node2 = node%index(side)
    
    lhs_tmp = ZBIG
    
    ! --- Condition on nodes
    call boundary_conditions_add_one_entry( &
         index_node, k_var, i_tor,          &
         index_node, k_var, i_tor,          &
         lhs_tmp, index_min, index_max, a_mat)

    ! --- Condition between nodes (d/ds or d/dt)
    call boundary_conditions_add_one_entry( &
         index_node2, k_var, i_tor,         &
         index_node2, k_var, i_tor,         &
         lhs_tmp, index_min, index_max, a_mat)
    
    return
  end subroutine apply_Dirichlet_BCs
  
  
  
  
  
  !******************************************************************************
  !******************************************************************************
  !****************** Routine to apply Mach-1 boundary conditions ***************
  !******************************************************************************
  !******************************************************************************
  subroutine apply_Mach1_BCs(rhs_loc, node, side, i_tor, index_min,index_max, a_mat)
  
    use mod_parameters
    use data_structure
    use phys_module, only: GAMMA, tauIC, central_density, mu_zero
    use mod_assembly, only : boundary_conditions_add_one_entry, boundary_conditions_add_RHS
    use mod_integer_types
    
    implicit none
    
    ! --- Routine variables
    real*8,                             intent(inout) :: rhs_loc(*)
    type (type_node),                   intent(in)    :: node
    integer,                            intent(in)    :: side ! == 2 for d/ds, == 3 for d/dt
    integer,                            intent(in)    :: i_tor
    integer,                            intent(in)    :: index_min, index_max
    type(type_SP_MATRIX)                              :: a_mat
    
    ! --- Internal variables
    integer                                           :: index_node,   index_node2
    real*8                                            :: mach1, dmach1_dVpar, dmach1_dpsis, dmach1_dus, dmach1_drho, dmach1_drhos, dmach1_dTi, dmach1_dTis
    real*8                                            :: mach1_ds, mach1_ds_Vpars, mach1_ds_Ti, mach1_ds_Tis
    
    ! --- Define node indices
    index_node  = node%index(1)             ! position of value
    index_node2 = node%index(side)          ! position of first deriative

    ! --- Define equations first, depends on which side we're on (d/ds or d/dt)
    if (side .eq. 2) then
      mach1          =   zbig*Vpar0                                                           & ! Main part
                       - zbig * R**2 * u0_s / ps0_s                                           & ! Vperp part
                       - zbig * R**2 * tau_IC / rho0    * Pi0_s / ps0_s                       & ! Vdia part
                       - zbig * direction / Btot                             * sqrt(GAMMA*Ti0)  ! Sound speed
      ! --- Linearise mach1
      dmach1_dVpar   =   zbig                                                          
      dmach1_dpsis   = + zbig * R**2 * u0_s / ps0_s**2                                        &
                       + zbig * R**2 * tau_IC / rho0    * Pi0_s / ps0_s**2
      dmach1_dus     = - zbig * R**2        / ps0_s
      dmach1_drho    = + zbig * R**2 * tau_IC / rho0**2 * Pi0_s / ps0_s                       &
                       - zbig * R**2 * tau_IC / rho0    * Ti0_s / ps0_s                      
      dmach1_drhos   = - zbig * R**2 * tau_IC / rho0    * Ti0   / ps0_s 
      dmach1_dTi     = - zbig * R**2 * tau_IC / rho0    * rho0_s/ ps0_s                       &
                       - zbig * direction / Btot * 0.5d0  * GAMMA            / sqrt(GAMMA*Ti0)       
      dmach1_dTis    = - zbig * R**2 * tau_IC / rho0    * rho0  / ps0_s                  
      ! --- d(mach1)/ds
      mach1_ds       =   zbig*Vpar0_s                                                         &
                       - zbig * direction / Btot * 0.5d0  * GAMMA    * Ti0_s / sqrt(GAMMA*Ti0)       
      ! --- Linearise d(mach1)/ds
      mach1_ds_Vpars =   zbig
      mach1_ds_Ti    = + zbig * direction / Btot * 0.25d0 * GAMMA**2 * Ti0_s / sqrt(GAMMA*Ti0)       
      mach1_ds_Tis   = - zbig * direction / Btot * 0.5d0  * GAMMA            / sqrt(GAMMA*Ti0)       
    else
      mach1          =   zbig*Vpar0                                                           & ! Main part
                       - zbig * R**2 * u0_t / ps0_t                                           & ! Vperp part
                       - zbig * R**2 * tau_IC / rho0    * Pi0_t / ps0_t                       & ! Vdia part
                       - zbig * direction / Btot                             * sqrt(GAMMA*Ti0)  ! Sound speed
      ! --- Linearise mach1
      dmach1_dVpar   =   zbig                                                          
      dmach1_dpsis   = + zbig * R**2 * u0_t / ps0_t**2                                        &
                       + zbig * R**2 * tau_IC / rho0    * Pi0_t / ps0_t**2
      dmach1_dus     = - zbig * R**2        / ps0_t
      dmach1_drho    = + zbig * R**2 * tau_IC / rho0**2 * Pi0_t / ps0_t                       &
                       - zbig * R**2 * tau_IC / rho0    * Ti0_t / ps0_t                      
      dmach1_drhos   = - zbig * R**2 * tau_IC / rho0    * Ti0   / ps0_t 
      dmach1_dTi     = - zbig * R**2 * tau_IC / rho0    * rho0_t/ ps0_t                       &
                       - zbig * direction / Btot * 0.5d0  * GAMMA            / sqrt(GAMMA*Ti0)       
      dmach1_dTis    = - zbig * R**2 * tau_IC / rho0    * rho0  / ps0_t                  
      ! --- d(mach1)/ds
      mach1_ds       =   zbig*Vpar0_t                                                         &
                       - zbig * direction / Btot * 0.5d0  * GAMMA    * Ti0_t / sqrt(GAMMA*Ti0)       
      ! --- Linearise d(mach1)/ds
      mach1_ds_Vpars =   zbig
      mach1_ds_Ti    = + zbig * direction / Btot * 0.25d0 * GAMMA**2 * Ti0_t / sqrt(GAMMA*Ti0)       
      mach1_ds_Tis   = - zbig * direction / Btot * 0.5d0  * GAMMA            / sqrt(GAMMA*Ti0)       
    endif
    
    
    ! --- Condition on nodes
    lhs_tmp = dmach1_dVpar
    call boundary_conditions_add_one_entry( &
         index_node, k_Vpar, i_tor,         &
         index_node, k_Vpar, i_tor,         &
         lhs_tmp, index_min, index_max, a_mat)

    lhs_tmp = dmach1_dpsis
    call boundary_conditions_add_one_entry( &
         index_node,  k_Vpar, i_tor,        &
         index_node2, k_psi,  i_tor,        &
         lhs_tmp, index_min, index_max, a_mat)

    lhs_tmp = dmach1_dus
    call boundary_conditions_add_one_entry( &
         index_node,  k_Vpar, i_tor,        &
         index_node2, k_u,    i_tor,        &
         lhs_tmp, index_min, index_max, a_mat)

    lhs_tmp = dmach1_drho
    call boundary_conditions_add_one_entry( &
         index_node,  k_Vpar, i_tor,        &
         index_node,  k_rho,  i_tor,        &
         lhs_tmp, index_min, index_max, a_mat)

    lhs_tmp = dmach1_drhos
    call boundary_conditions_add_one_entry( &
         index_node,  k_Vpar, i_tor,        &
         index_node2, k_rho,  i_tor,        &
         lhs_tmp, index_min, index_max, a_mat)

    lhs_tmp = dmach1_dTi
    call boundary_conditions_add_one_entry( &
         index_node,  k_Vpar, i_tor,        &
         index_node,  k_Ti,   i_tor,        &
         lhs_tmp, index_min, index_max, a_mat)

    lhs_tmp = dmach1_dTis
    call boundary_conditions_add_one_entry( &
         index_node,  k_Vpar, i_tor,        &
         index_node2, k_Ti,   i_tor,        &
         lhs_tmp, index_min, index_max, a_mat)

      if (i_tor .eq. 1) then
        rhs_tmp = - mach1
      else
        rhs_tmp = 0.d0
      endif
      call boundary_conditions_add_RHS(       &
           index_node, k_Vpar, i_tor,         &
           index_min, index_max,              &
           RHS_loc, rhs_tmp,                  &
           a_mat%i_tor_min, a_mat%i_tor_max)

    ! --- Condition between nodes (d/ds or d/dt)
    lhs_tmp = mach1_ds_Vpars
    call boundary_conditions_add_one_entry( &
         index_node2, k_Vpar, i_tor,        &
         index_node2, k_Vpar, i_tor,        &
         lhs_tmp, index_min, index_max, a_mat)

    lhs_tmp = mach1_ds_Ti
    call boundary_conditions_add_one_entry( &
         index_node2, k_Vpar, i_tor,        &
         index_node,  k_Ti,   i_tor,        &
         lhs_tmp, index_min, index_max, a_mat)

    lhs_tmp = mach1_ds_Tis
    call boundary_conditions_add_one_entry( &
         index_node2, k_Vpar, i_tor,        &
         index_node2, k_Ti,   i_tor,        &
         lhs_tmp, index_min, index_max, a_mat)

    if (i_tor .eq. 1) then
      rhs_tmp = - mach1_ds
    else
      rhs_tmp = 0.d0
    endif
    
    call boundary_conditions_add_RHS(       &
           index_node2, k_Vpar, i_tor,        &
           index_min, index_max,              &
           RHS_loc, rhs_tmp,                  &
           a_mat%i_tor_min, a_mat%i_tor_max)
  
    return
  end subroutine apply_Mach1_BCs
  
  
  
  
  
  
  
  
  
  
  !******************************************************************************
  !******************************************************************************
  !****************** Routine to apply STANGEBY boundary conditions on U ********
  !******************************************************************************
  !******************************************************************************
  subroutine apply_U_sheath(rhs_loc, node, side, i_tor, index_min,index_max, a_mat)
  
    use mod_parameters
    use data_structure
    use phys_module, only: GAMMA, tauIC, central_density, mu_zero, F0, FF_0
    use mod_assembly, only : boundary_conditions_add_one_entry, boundary_conditions_add_RHS
    use mod_integer_types
    
    implicit none
    
    ! --- Routine variables
    real*8,                             intent(inout) :: rhs_loc(*)
    type (type_node),                   intent(in)    :: node
    integer,                            intent(in)    :: side ! == 2 for d/ds, == 3 for d/dt
    integer,                            intent(in)    :: i_tor
    integer,                            intent(in)    :: index_min, index_max
    type(type_SP_MATRIX)                              :: a_mat
    
    ! --- Internal variables
    integer                                           :: index_node,   index_node2
    real*8                                            :: constmp
    real*8                                            :: sheath_u, dsheath_u, dsheath_T
    real*8                                            :: sheath_ds, dsheath_ds_us, dsheath_ds_Ts
    real*8                                            :: sign_tmp
    
    ! --- Define node indices
    index_node  = node%index(1)             ! position of value
    index_node2 = node%index(side)          ! position of first deriative

    ! --- Define equations before MURGE and non-MURGE fork
    
    ! --- The constant factor
    !constmp   = - tauIC/abs(tauIC) * 3.1d-8/sqrt(mu_zero*rho_norm) ! constant from the formula U = -Te/e * ln(mi/(2pi.me))**0.5  (Guido)
    !constmp   = - 3.1d-8/sqrt(mu_zero*rho_norm)/F0 ! constant from the formula U = -Te/e * ln(mi/(2pi.me))**0.5  (Guido/Stangeby at wall itself, not into plasma)
    constmp   = - 0.7d-9/sqrt(mu_zero*rho_norm)/F0 ! constant from the formula U = -0.07 Te/e (Stangeby at sheath edge into plasma)
    
    ! --- The sign depends on direction of Field and direction of current
    sign_tmp  = FF_0/abs(FF_0) * tauIC/abs(tauIC)
    constmp   = sign_tmp * constmp
    
    sheath_u  = zbig * ( U0 - constmp * Ti0 )
    dsheath_u = zbig
    dsheath_T =-zbig * constmp
    if (side .eq. 2) then
      sheath_ds     = zbig * ( U0_s - constmp * Ti0_s )
      dsheath_ds_us = zbig
      dsheath_ds_Ts =-zbig * constmp
    else
      sheath_ds     = zbig * ( U0_t - constmp * Ti0_t )
      dsheath_ds_us = zbig
      dsheath_ds_Ts =-zbig * constmp
    endif

    ! --- Condition on nodes
    lhs_tmp = dsheath_u
    call boundary_conditions_add_one_entry( &
         index_node, k_u, i_tor,            &
         index_node, k_u, i_tor,            &
         lhs_tmp, index_min, index_max, a_mat)

    lhs_tmp = dsheath_T
    call boundary_conditions_add_one_entry( &
         index_node, k_u,  i_tor,           &
         index_node, k_Ti, i_tor,           &
         lhs_tmp, index_min, index_max, a_mat)

      if (i_tor .eq. 1) then
        rhs_tmp = - sheath_u
      else
        rhs_tmp = 0.d0
      endif
      call boundary_conditions_add_RHS(       &
           index_node, k_u, i_tor,            &
           index_min, index_max,              &
           RHS_loc, rhs_tmp,                  &
           a_mat%i_tor_min, a_mat%i_tor_max)
      
    ! --- Condition between nodes
    lhs_tmp = dsheath_ds_us
    call boundary_conditions_add_one_entry( &
         index_node2, k_u, i_tor,           &
         index_node2, k_u, i_tor,           &
         lhs_tmp, index_min, index_max, a_mat)

    lhs_tmp = dsheath_ds_Ts
    call boundary_conditions_add_one_entry( &
         index_node2, k_u, i_tor,           &
         index_node2, k_Ti, i_tor,          &
         lhs_tmp, index_min, index_max, a_mat)

      if (i_tor .eq. 1) then
        rhs_tmp = - sheath_ds
      else
        rhs_tmp = 0.d0
      endif
      call boundary_conditions_add_RHS(       &
           index_node2, k_u, i_tor,           &
           index_min, index_max,              &
           RHS_loc, rhs_tmp,                  &
           a_mat%i_tor_min, a_mat%i_tor_max)

    return
  end subroutine apply_U_sheath
  
  
  
  
  
  
end module mod_boundary_conditions
