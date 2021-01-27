module mod_boundary_conditions

  !*******************************************************************************
  !************ Define global variables for all internal routines ****************
  !*******************************************************************************
  ! --- ZBIG parameter to make equations "more important" than element_matrix equations
  real*8  		        :: zbig = 1.d10, zbig_backup
  ! --- R,Z variables
  real*8			:: R, R_s, R_t
  real*8			:: Z, Z_s, Z_t
  real*8			:: xjac
  ! --- Variable numbers
  integer, parameter		:: k_psi  = 1
  integer, parameter		:: k_u    = 2
  integer, parameter		:: k_Vpar = 7
  integer, parameter		:: k_Ti   = 6
  integer, parameter		:: k_Te   = 8
  ! --- Variables
  real*8			:: ps0,   ps0_s,   ps0_t, ps0_x, ps0_y, grad_psi, Btot
  real*8			::        u0_s,    u0_t,  u0_x,  u0_y
  real*8			:: Vpar0, Vpar0_s, Vpar0_t
  real*8			:: Ti0,   Ti0_s,   Ti0_t
  real*8			:: Te0,   Te0_s,   Te0_t
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
  subroutine boundary_conditions( my_id, node_list, element_list, bnd_node_list,           &
                                  local_elms, n_local_elms, index_min, index_max, rhs_loc, &
                                  xpoint2, xcase2,                                         &
                                  R_axis, Z_axis, psi_axis,                                &
                                  psi_bnd, R_xpoint, Z_xpoint, psi_xpoint,                 &
                                  gmres, solve_only, ijA_index, ijA_size, irn_jcn,         & 
                                  irn, jcn, A_mat, i_tor_min, i_tor_max )

    use data_structure
    use vacuum, ONLY: is_freebound
    use phys_module, only: F0, GAMMA, freeboundary, tokamak_device,                  &
                           RMP_on, psi_RMP_cos, dpsi_RMP_cos_dR, dpsi_RMP_cos_dZ,    &
                           psi_RMP_sin, dpsi_RMP_sin_dR, dpsi_RMP_sin_dZ,            &
                           t_now, RMP_start_time, tstep, RMP_har_cos, RMP_har_sin,   &
                           RMP_growth_rate, RMP_ramp_up_time, T_min, grid_to_wall,   &
                           n_wall_blocks, keep_n0_const
    use mpi_mod
    use mod_locate_irn_jcn

    implicit none
    !include 'mpif.h'

    ! --- Routine parameters
    integer,                   intent(in)    :: my_id
    type (type_node_list),     intent(in)    :: node_list
    type (type_element_list),  intent(in)    :: element_list
    type (type_bnd_node_list), intent(in)    :: bnd_node_list
    integer,                   intent(in)    :: local_elms(*)
    integer,                   intent(in)    :: n_local_elms
    integer,                   intent(in)    :: index_min
    integer,                   intent(in)    :: index_max
    real*8,                    intent(inout) :: rhs_loc(*)
    logical,                   intent(in)    :: xpoint2
    integer,                   intent(in)    :: xcase2
    real*8,                    intent(in)    :: R_axis
    real*8,                    intent(in)    :: Z_axis
    real*8,                    intent(in)    :: psi_axis
    real*8,                    intent(in)    :: psi_bnd
    real*8,                    intent(in)    :: R_xpoint(2)
    real*8,                    intent(in)    :: Z_xpoint(2)
    real*8,                    intent(in)    :: psi_xpoint(2)
    logical,                   intent(in)    :: gmres
    logical,                   intent(in)    :: solve_only
    integer,                   intent(in)    :: i_tor_min, i_tor_max 
    integer, allocatable,      intent(in)    :: ijA_index(:,:), ijA_size(:), irn_jcn(:,:) 
    integer, allocatable,      intent(inout) :: irn(:), jcn(:) 
    real*8, allocatable,       intent(inout) :: A_mat(:) 

    ! --- Internal parameters
    real*8  :: mach1, dmach1, d2mach1_dTi, d2mach1_dTe, mach_u, dmach_u, dmach_rho
    integer :: i, i_tor, iv, inode, k_var, side
    integer :: index_node, index_node2, ielm, index_tmp
    integer :: ijA_position,ijA_position2, ilarge2, ilarge_vv, ilarge_vTi, ilarge_vTe, ilarge_vus
    integer :: ilarge_vsvs, ilarge_vsTi, ilarge_vsTe, ilarge_vsTis, ilarge_vsTes
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
    if (RMP_on .and. (n_tor .ge. 3)) then
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
        establish_RMP = 0.5d-3 
      else
        establish_RMP = 0.0
      endif
    
      do j=1, bnd_node_list%n_bnd_nodes  
        psi_RMP_cos1(j)     =  psi_RMP_cos(j)    * establish_RMP
        dpsi_RMP_cos_dR1(j) = dpsi_RMP_cos_dR(j) * establish_RMP
        dpsi_RMP_cos_dZ1(j) = dpsi_RMP_cos_dZ(j) * establish_RMP
        psi_RMP_sin1(j)     =  psi_RMP_sin(j)    * establish_RMP
        dpsi_RMP_sin_dR1(j) = dpsi_RMP_sin_dR(j) * establish_RMP
        dpsi_RMP_sin_dZ1(j) = dpsi_RMP_sin_dZ(j) * establish_RMP
      end do
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
	    
	    do i_tor=i_tor_min, i_tor_max
              if (keep_n0_const  .and.  i_tor .eq. 1 ) then
                zbig = 1.d15
              else
                zbig = zbig_backup
              endif

              do k_var=1, n_var

                ! --------------------------------------------------------------------------------------------------------------
                ! ------------------------------------ the targets (in case of x-point grid) -----------------------------------
                ! --------------------------------------------------------------------------------------------------------------
                if     ((node_list%node(inode)%boundary .eq.  1) &
                   .or. (node_list%node(inode)%boundary .eq. 11) &
                   .or. (node_list%node(inode)%boundary .eq.  9) &
                   .or. (node_list%node(inode)%boundary .eq. 19) &
                   .or. (node_list%node(inode)%boundary .eq.  3) &
                   .or. (node_list%node(inode)%boundary .eq.  4)) then
                      
                  ! --- Which side is this? 2 => d/ds, 3 => d/dt
                  side = 2

                  ! --- Special field direction for grid patches
                  ! --- (hopefully temporary until Vpar is treated in boundary_matrix_open)
                  if ( (grid_to_wall) .and. (n_wall_blocks .gt. 0) ) then
                    direction = 1
                    if (node_list%node(inode)%boundary .eq. 11) direction = -direction
                    if (node_list%node(inode)%boundary .eq. 19) direction = -direction
                  endif
                            
                  ! ---------------------------------------------
                  ! --- Apply RMP on target (only depends on 's')
                  if (      RMP_on							&
		      .and. (k_var .eq. 1)						&
		      .and. ((i_tor.eq.RMP_har_cos) .or. (i_tor.eq.RMP_har_sin))	&
		      .and. (.not. freeboundary)					) then
                        		 
		      call apply_RMP_BCs(rhs_loc, node_list%node(inode), side, i_tor,	        &
		                         psi_RMP_cos1, dpsi_RMP_cos_dR1, dpsi_RMP_cos_dZ1,	&
		                         psi_RMP_sin1, dpsi_RMP_sin_dR1, dpsi_RMP_sin_dZ1,	&
                                         index_min,index_max,i_tor_min,i_tor_max,               & 
                                         ijA_index, ijA_size, irn_jcn, irn, jcn, A_mat)

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
                  if (                          &
                           apply_on_psi         & 
                      .or. apply_on_current                             &
                      .or. (k_var .eq. 2)       &
                      .or. (k_var .eq. 4)       &
                      ) apply_dirichlet = .true.

		  ! --- Apply Dirichlet if required
		  if (apply_dirichlet) then
		    call apply_Dirichlet_BCs(node_list%node(inode), side, k_var,i_tor, index_min,index_max, gmres, solve_only, & 
                                             i_tor_min, i_tor_max, ijA_index, ijA_size, irn_jcn, irn, jcn, A_mat)
                  endif

                  ! --------------
                  ! --- Mach-1 BCs
                  if (k_var .eq. 7) then
                    call apply_Mach1_BCs(rhs_loc, node_list%node(inode), side, i_tor, index_min,index_max, gmres, solve_only, & 
                                         i_tor_min, i_tor_max, ijA_index, ijA_size, irn_jcn, irn, jcn, A_mat)
                  endif

                endif

                
                ! --------------------------------------------------------------------------------------------------------------
                ! ------------------------- the non-targets open field-lines (for grid_xpoint_wall) ----------------------------
                ! --------------------------------------------------------------------------------------------------------------
                if    ((node_list%node(inode)%boundary .eq.     5) &
                     .or. (node_list%node(inode)%boundary .eq. 15) &
                     .or. (node_list%node(inode)%boundary .eq.  9) &
                     .or. (node_list%node(inode)%boundary .eq. 19)) then

                  ! --- Which side is this? 2 => d/ds, 3 => d/dt
                  side = 3
                  
                  ! --- Special field direction for grid patches
                  ! --- (hopefully temporary until Vpar is treated in boundary_matrix_open)
                  if ( (grid_to_wall) .and. (n_wall_blocks .gt. 0) ) then
                    direction = 1
                    if (node_list%node(inode)%boundary .eq. 15) direction = -direction
                    if (node_list%node(inode)%boundary .eq. 19) direction = -direction
                  endif
                  
		  ! ---------------------------------------------
		  ! --- Apply RMP on target (only depends on 's')
                  if (      RMP_on							&
		      .and. (k_var .eq. 1)						&
		      .and. ((i_tor.eq.RMP_har_cos) .or. (i_tor.eq.RMP_har_sin))	&
		      .and. (.not. freeboundary)					) then
                        		 
		      call apply_RMP_BCs(rhs_loc, node_list%node(inode), side, i_tor,	        &
		                         psi_RMP_cos1, dpsi_RMP_cos_dR1, dpsi_RMP_cos_dZ1,	&
		                         psi_RMP_sin1, dpsi_RMP_sin_dR1, dpsi_RMP_sin_dZ1,	&
                                         index_min,index_max,i_tor_min,i_tor_max,               & 
                                         ijA_index, ijA_size, irn_jcn, irn, jcn, A_mat)


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
                  if (                          &
                           apply_on_psi         & 
                      .or. apply_on_current                             &
                      .or. (k_var .eq. 2)       &
                      .or. (k_var .eq. 4)       &
                      ) apply_dirichlet = .true.

		  ! --- Apply Dirichlet if required
		  if (apply_dirichlet) then
		    call apply_Dirichlet_BCs(node_list%node(inode), side, k_var,i_tor, index_min,index_max, gmres, solve_only, & 
                                             i_tor_min, i_tor_max, ijA_index, ijA_size, irn_jcn, irn, jcn, A_mat)
                  endif

                  ! --------------
                  ! --- Mach-1 BCs
                  if (k_var .eq. 7) then
                    call apply_Mach1_BCs(rhs_loc, node_list%node(inode), side, i_tor, index_min,index_max, gmres, solve_only, & 
                                         i_tor_min, i_tor_max, ijA_index, ijA_size, irn_jcn, irn, jcn, A_mat)
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
                  
                  ! ---------------------------------------------
                  ! --- Apply RMP on target (only depends on 't')
                  if (      RMP_on                                                      &
                      .and. (k_var .eq. 1)                                              &
                      .and. ((i_tor.eq.RMP_har_cos) .or. (i_tor.eq.RMP_har_sin))        &
                      .and. (.not. freeboundary)                                        ) then
                                         
                      call apply_RMP_BCs(rhs_loc, node_list%node(inode), side, i_tor,           &
                                         psi_RMP_cos1, dpsi_RMP_cos_dR1, dpsi_RMP_cos_dZ1,      &
                                         psi_RMP_sin1, dpsi_RMP_sin_dR1, dpsi_RMP_sin_dZ1,      &
                                         index_min,index_max,i_tor_min,i_tor_max,               & 
                                         ijA_index, ijA_size, irn_jcn, irn, jcn, A_mat)

                  endif
                  
                  ! -----------------------------------------------
                  ! --- Dirichlet BCs (or Neumann if commented out)
                  apply_dirichlet = .false.
                  
                  ! --- Determine if we are on the private or the inner boundary
                  on_private            = .false.
                  on_inner              = .false.
                  on_inner_or_private   = .false.
                  if ((xcase2 .ne. 3) .and. (ps0 .lt. psi_bnd)) then
                    on_inner_or_private = .true.
                    on_private          = .true.
                  endif
                  if  (xcase2 .eq. 3) then
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
                      if ( (RMP_on)       .and. (i_tor .lt. RMP_har_cos) ) apply_on_psi = .true.
                      if ( (RMP_on)       .and. (i_tor .gt. RMP_har_sin) ) apply_on_psi = .true.
                    endif
                  
                  apply_on_current = .false.
                  if ((k_var .eq. 3) .and. (.not. is_freebound(i_tor,k_var))) apply_on_current = .true.
                  
                  ! Apply conditions to which variables and where?
                  if (                                                                  &
                            (apply_on_psi)                                              &
                      .or.  (apply_on_current)                                          &
                      .or.  (k_var .eq. 2)                      &
                      .or.  (k_var .eq. 4)                                              &
                      .or.  (k_var .eq. 5)                                              &
                      .or.  (k_var .eq. 6)                                              &
                      !.or.( (k_var .eq. 5) .and. (on_private) )                                & 
                      !.or.( (k_var .eq. 6) .and. (on_private) )                                & 
                      .or.  (k_var .eq. 7)                                              &
                      .or.  (k_var .eq. 8)                                              &
                      ) apply_dirichlet = .true.

                  if (apply_dirichlet) then
                    call apply_Dirichlet_BCs(node_list%node(inode), side, k_var,i_tor, index_min,index_max, gmres, solve_only, &
                                             i_tor_min, i_tor_max, ijA_index, ijA_size, irn_jcn, irn, jcn, A_mat)
                  endif

                endif

                ! ----------------------------------------------------------------------------------------------------
                ! ------------------------------------ the flux-surface boundaries -----------------------------------
                ! ----------------------------------------------------------------------------------------------------
                if    ((node_list%node(inode)%boundary .eq. 20) &
                  .or. (node_list%node(inode)%boundary .eq. 21)) then

                  ! --- We do both sides! 2 => d/ds, 3 => d/dt
                  side = 2
                  
                  ! ---------------------------------------------
                  ! --- Apply RMP on target (only depends on 't')
                  if (      RMP_on 							&
		      .and. (k_var .eq. 1)						&
		      .and. ((i_tor.eq.RMP_har_cos) .or. (i_tor.eq.RMP_har_sin))	&
		      .and. (.not. freeboundary) 					) then
                        		 
		      call apply_RMP_BCs(rhs_loc, node_list%node(inode), side, i_tor,	        &
		                         psi_RMP_cos1, dpsi_RMP_cos_dR1, dpsi_RMP_cos_dZ1,	&
		                         psi_RMP_sin1, dpsi_RMP_sin_dR1, dpsi_RMP_sin_dZ1,	&
                                         index_min,index_max,i_tor_min,i_tor_max,               & 
                                         ijA_index, ijA_size, irn_jcn, irn, jcn, A_mat)
                  
                  endif
                  
                  ! -----------------------------------------------
                  ! --- Dirichlet BCs (or Neumann if commented out)
                  apply_dirichlet = .false.
                  
                  ! --- Determine if we are on the private or the inner boundary
                  on_private            = .false.
                  on_inner              = .false.
                  on_inner_or_private   = .false.
                  if ((xcase2 .ne. 3) .and. (ps0 .lt. psi_bnd)) then
                    on_inner_or_private = .true.
                    on_private          = .true.
                  endif
                  if  (xcase2 .eq. 3) then
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
                      if ( (RMP_on)       .and. (i_tor .lt. RMP_har_cos) ) apply_on_psi = .true.
                      if ( (RMP_on)       .and. (i_tor .gt. RMP_har_sin) ) apply_on_psi = .true.
                    endif
                  
                  apply_on_current = .false.
                  if ((k_var .eq. 3) .and. (.not. is_freebound(i_tor,k_var))) apply_on_current = .true.
                  
                  ! Apply conditions to which variables and where?
                  if (                                                                  &
                            (apply_on_psi)                                              &
                      .or.  (apply_on_current)                                          &
                      .or.  (k_var .eq. 2)                                              &
                      .or.  (k_var .eq. 4)                                              &
                      .or.  (k_var .eq. 5)                                              &
                      .or.  (k_var .eq. 6)                                              &
                      !.or.( (k_var .eq. 5) .and. (on_private) )                                & 
                      !.or.( (k_var .eq. 6) .and. (on_private) )                                & 
                      .or.  (k_var .eq. 7)                                              &
                      .or.  (k_var .eq. 8)                                              &
                      ) apply_dirichlet = .true.

                  if (apply_dirichlet) then
                    call apply_Dirichlet_BCs(node_list%node(inode), side, k_var,i_tor, index_min,index_max, gmres, solve_only, & 
                                             i_tor_min, i_tor_max, ijA_index, ijA_size, irn_jcn, irn, jcn, A_mat)
                  endif

                  ! --- We do both sides! 2 => d/ds, 3 => d/dt
                  side = 3
                  
                  ! ---------------------------------------------
                  ! --- Apply RMP on target (only depends on 't')
                  if (      RMP_on                                                      &
                      .and. (k_var .eq. 1)                                              &
                      .and. ((i_tor.eq.RMP_har_cos) .or. (i_tor.eq.RMP_har_sin))        &
                      .and. (.not. freeboundary)                                        ) then
                                         
                      call apply_RMP_BCs(rhs_loc, node_list%node(inode), side, i_tor,           &
                                         psi_RMP_cos1, dpsi_RMP_cos_dR1, dpsi_RMP_cos_dZ1,      &
                                         psi_RMP_sin1, dpsi_RMP_sin_dR1, dpsi_RMP_sin_dZ1,      &
                                         index_min,index_max,i_tor_min,i_tor_max,               & 
                                         ijA_index, ijA_size, irn_jcn, irn, jcn, A_mat)
                  endif
                  
                  ! -----------------------------------------------
                  ! --- Dirichlet BCs (or Neumann if commented out)
                  apply_dirichlet = .false.
                  
                  ! --- Determine if we are on the private or the inner boundary
                  on_private            = .false.
                  on_inner              = .false.
                  on_inner_or_private   = .false.
                  if ((xcase2 .ne. 3) .and. (ps0 .lt. psi_bnd)) then
                    on_inner_or_private = .true.
                    on_private          = .true.
                  endif
                  if  (xcase2 .eq. 3) then
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
                      if ( (RMP_on)       .and. (i_tor .lt. RMP_har_cos) ) apply_on_psi = .true.
                      if ( (RMP_on)       .and. (i_tor .gt. RMP_har_sin) ) apply_on_psi = .true.
                    endif
                  
                  apply_on_current = .false.
                  if ((k_var .eq. 3) .and. (.not. is_freebound(i_tor,k_var))) apply_on_current = .true.
                  
                  ! Apply conditions to which variables and where?
                  if (                                                                  &
                            (apply_on_psi)                                              &
                      .or.  (apply_on_current)                                          &
                      .or.  (k_var .eq. 2)                                              &
                      .or.  (k_var .eq. 4)                                              &
                      .or.  (k_var .eq. 5)                                              &
                      .or.  (k_var .eq. 6)                                              &
                      !.or.( (k_var .eq. 5) .and. (on_private) )                                & 
                      !.or.( (k_var .eq. 6) .and. (on_private) )                                & 
                      .or.  (k_var .eq. 7)                                              &
                      .or.  (k_var .eq. 8)                                              &
                      ) apply_dirichlet = .true.

		  if (apply_dirichlet) then
                    call apply_Dirichlet_BCs(node_list%node(inode), side, k_var,i_tor, index_min,index_max, gmres, solve_only, &
                                             i_tor_min, i_tor_max, ijA_index, ijA_size, irn_jcn, irn, jcn, A_mat)
                  endif

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
  
    use data_structure
    use phys_module, only: F0, xpoint, xcase, tokamak_device, T_min
    
    implicit none
    
    ! --- Routine variables
    type (type_node),	intent(in)    :: node
    real*8,		intent(in)    :: R_axis
    real*8,		intent(in)    :: Z_axis
    real*8,		intent(in)    :: R_xpoint(2)
    real*8,		intent(in)    :: Z_xpoint(2)
    real*8,		intent(in)    :: psi_bnd
    
    real*8 :: alpha, R_inside, Z_inside
  
    ! --- Define (R,Z) coords and Jacobian
    R	      = node%x(1,1,1)
    R_s       = node%x(1,2,1)
    R_t       = node%x(1,3,1)
    Z	      = node%x(1,1,2)
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
    U0_s      = node%values(1,2,2)
    U0_t      = node%values(1,3,2)
    u0_x      = (  Z_t*u0_s - Z_s*u0_t) / xjac
    u0_y      = (- R_t*u0_s + R_s*u0_t) / xjac

    ! --- Define Ti variables
    Ti0       = max(node%values(1,1,6),T_min)
    Ti0_s     = node%values(1,2,6)
    Ti0_t     = node%values(1,3,6)
    
    ! --- Define Te variables
    Te0       = max(node%values(1,1,8),T_min)
    Te0_s     = node%values(1,2,8)
    Te0_t     = node%values(1,3,8)
    
    ! --- Define Vpar variables
    Vpar0     = node%values(1,1,k_Vpar)
    Vpar0_s   = node%values(1,2,k_Vpar)
    Vpar0_t   = node%values(1,3,k_Vpar)

    ! --- Define direction of Vpar on target. Careful, using ps0_x/abs(ps0_x) can be treacherous.
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
    if (xcase .eq. 2) then
      direction = -direction
    else if ((xcase .eq. 3).and.(Z .gt. Z_axis +0.1) .and. ( R .gt.R_xpoint(2))) then
      direction = -1.
    else if ((xcase .eq. 3) .and. (Z .gt. Z_axis +0.1) .and. (R .lt. R_xpoint(2))) then
      direction = +1.
    end if

    return
  end subroutine construct_variables
  
  
  
  
  
  
  !******************************************************************************
  !******************************************************************************
  !********* Routine to apply RMP perturbation on boundary conditions ***********
  !******************************************************************************
  !******************************************************************************
  subroutine apply_RMP_BCs(rhs_loc, node, side, i_tor, 				&
		           psi_RMP_cos1, dpsi_RMP_cos_dR1, dpsi_RMP_cos_dZ1,	&
		           psi_RMP_sin1, dpsi_RMP_sin_dR1, dpsi_RMP_sin_dZ1,	&
                           index_min,index_max,i_tor_min, i_tor_max,            & 
                           ijA_index, ijA_size, irn_jcn, irn, jcn, A_mat)
  
    use mod_parameters
    use data_structure
    use phys_module, only: RMP_har_cos, RMP_har_sin
    use mod_locate_irn_jcn
    
    implicit none
    
    ! --- Routine variables
    real*8,	          intent(inout) :: rhs_loc(*)
    type (type_node),     intent(in)    :: node
    integer,	          intent(in)    :: side ! == 2 for d/ds, == 3 for d/dt
    integer,	          intent(in)    :: i_tor
    real*8,	          intent(in)    :: psi_RMP_cos1(*), dpsi_RMP_cos_dR1(*), dpsi_RMP_cos_dZ1(*)
    real*8,	          intent(in)    :: psi_RMP_sin1(*), dpsi_RMP_sin_dR1(*), dpsi_RMP_sin_dZ1(*)
    integer,	          intent(in)    :: index_min, index_max
    integer,              intent(in)    :: i_tor_min, i_tor_max 
    integer, allocatable, intent(in)    :: ijA_index(:,:), ijA_size(:), irn_jcn(:,:) 
    integer, allocatable, intent(inout) :: irn(:), jcn(:) 
    real*8,  allocatable, intent(inout) :: A_mat(:) 
    
    ! --- Internal variables
    integer				:: index_node,   index_node2, index_tmp
    integer				:: ijA_position, ijA_position2
    integer				:: ilarge_vp,    ilarge_vp2
    real*8				:: delta_psi_rmp, delta_psi_rmp_dR, delta_psi_rmp_dZ, delta_psi_rmp_dl
    integer				:: n_tor_local 
    
    n_tor_local = i_tor_max - i_tor_min +1
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
    if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then
      call locate_irn_jcn(index_node,index_node,index_min,index_max,ijA_position,ijA_index, ijA_size, irn_jcn)
      ilarge_vp  = ijA_position  - 1 + ((k_psi-1)*n_tor_local + i_tor-i_tor_min) * n_var*n_tor_local + & 
        (k_psi-1)*n_tor_local + i_tor -i_tor_min + 1
      
      irn(ilarge_vp) = n_tor_local * n_var * (index_node-1) + (k_psi-1)*n_tor_local + i_tor -i_tor_min + 1
      jcn(ilarge_vp) = n_tor_local * n_var * (index_node-1) + (k_psi-1)*n_tor_local + i_tor -i_tor_min + 1
      A_mat(ilarge_vp)   = ZBIG
      
      index_tmp = n_tor_local*n_var * (index_node-1) + (k_psi-1)*n_tor_local + i_tor -i_tor_min + 1
      Rhs_loc(index_tmp) = ZBIG * delta_psi_rmp
    endif
    
    ! --- Condition between nodes (d/ds or d/dt)
    if ((index_node2 .ge. index_min) .and. (index_node2 .le. index_max)) then			      
      call locate_irn_jcn(index_node2,index_node2,index_min,index_max,ijA_position2,ijA_index, ijA_size, irn_jcn)
      
      ilarge_vp2  = ijA_position2 - 1 + ((k_psi-1)*n_tor_local + i_tor-i_tor_min) * n_var*n_tor_local + & 
        (k_psi-1)*n_tor_local + i_tor -i_tor_min + 1

      irn(ilarge_vp2) = n_tor_local * n_var * (index_node2-1) + (k_psi-1)*n_tor_local + i_tor -i_tor_min + 1
      jcn(ilarge_vp2) = n_tor_local * n_var * (index_node2-1) + (k_psi-1)*n_tor_local + i_tor -i_tor_min + 1
      A_mat(ilarge_vp2)   = ZBIG
      
      index_tmp = n_tor_local*n_var * (index_node2-1) + (k_psi-1)*n_tor_local + i_tor -i_tor_min + 1
      Rhs_loc(index_tmp) = ZBIG * delta_psi_rmp_dl
    endif
  
    return
  end subroutine apply_RMP_BCs
  
  
  
  
  
  !******************************************************************************
  !******************************************************************************
  !***************** Routine to apply Dirichlet boundary conditions *************
  !******************************************************************************
  !******************************************************************************
  subroutine apply_Dirichlet_BCs(node, side, k_var,i_tor, index_min,index_max, gmres, solve_only, &
                                 i_tor_min, i_tor_max, ijA_index, ijA_size, irn_jcn, irn, jcn, A_mat)
  
    use mod_parameters
    use data_structure
    use phys_module, only: RMP_har_cos, RMP_har_sin
    use mod_locate_irn_jcn
    
    implicit none
    
    ! --- Routine variables
    type (type_node),     intent(in)    :: node
    integer,	          intent(in)    :: side ! == 2 for d/ds, == 3 for d/dt
    integer,	          intent(in)    :: k_var
    integer,	          intent(in)    :: i_tor
    integer,	          intent(in)    :: index_min, index_max
    logical,	          intent(in)    :: gmres, solve_only
    integer,              intent(in)    :: i_tor_min, i_tor_max 
    integer, allocatable, intent(in)    :: ijA_index(:,:), ijA_size(:), irn_jcn(:,:) 
    integer, allocatable, intent(inout) :: irn(:), jcn(:) 
    real*8,  allocatable, intent(inout) :: A_mat(:) 
    
    ! --- Internal variables
    integer				:: index_node,   index_node2
    integer				:: ijA_position
    integer				:: ilarge
    logical				:: is_local
    integer				:: n_tor_local 
    
    n_tor_local = i_tor_max - i_tor_min +1
    ! --- Get nodes index
    index_node  = node%index(1)
    index_node2 = node%index(side)
    
    ! --- Condition on nodes
    if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then
      call locate_irn_jcn(index_node,index_node,index_min,index_max,ijA_position,ijA_index, ijA_size, irn_jcn)
      ilarge = ijA_position - 1 + ((k_var-1)*n_tor_local + i_tor-i_tor_min) * n_var*n_tor_local + & 
        (k_var-1)*n_tor_local + i_tor -i_tor_min + 1

      irn(ilarge) = n_tor_local * n_var * (index_node-1) + (k_var-1)*n_tor_local + i_tor -i_tor_min + 1
      jcn(ilarge) = n_tor_local * n_var * (index_node-1) + (k_var-1)*n_tor_local + i_tor -i_tor_min + 1
      A_mat(ilarge)   = zbig
    endif
    ! --- Condition between nodes (d/ds)
    if ((index_node2 .ge. index_min) .and. (index_node2 .le. index_max)) then
      call locate_irn_jcn(index_node2,index_node2,index_min,index_max,ijA_position,ijA_index, ijA_size, irn_jcn)
      ilarge = ijA_position - 1 + ((k_var-1)*n_tor_local + i_tor-i_tor_min) * n_var*n_tor_local + & 
        (k_var-1)*n_tor_local + i_tor -i_tor_min + 1

      irn(ilarge) = n_tor_local * n_var * (index_node2-1) + (k_var-1)*n_tor_local + i_tor -i_tor_min + 1
      jcn(ilarge) = n_tor_local * n_var * (index_node2-1) + (k_var-1)*n_tor_local + i_tor -i_tor_min + 1
      A_mat(ilarge)   = zbig
    endif
  
    return
  end subroutine apply_Dirichlet_BCs
  
  
  
  
  
  !******************************************************************************
  !******************************************************************************
  !****************** Routine to apply Mach-1 boundary conditions ***************
  !******************************************************************************
  !******************************************************************************
  subroutine apply_Mach1_BCs(rhs_loc, node, side, i_tor, index_min,index_max, gmres, solve_only, &
                             i_tor_min, i_tor_max, ijA_index, ijA_size, irn_jcn, irn, jcn, A_mat)
  
  
    use mod_parameters
    use data_structure
    use phys_module, only: GAMMA
    use mod_locate_irn_jcn
    
    implicit none
    
    ! --- Routine variables
    real*8,	          intent(inout) :: rhs_loc(*)
    type (type_node),     intent(in)    :: node
    integer,	          intent(in)    :: side ! == 2 for d/ds, == 3 for d/dt
    integer,	          intent(in)    :: i_tor
    integer,	          intent(in)    :: index_min, index_max
    logical,	          intent(in)    :: gmres, solve_only
    integer,              intent(in)    :: i_tor_min, i_tor_max 
    integer, allocatable, intent(in)    :: ijA_index(:,:), ijA_size(:), irn_jcn(:,:) 
    integer, allocatable, intent(inout) :: irn(:), jcn(:) 
    real*8,  allocatable, intent(inout) :: A_mat(:) 
    
    ! --- Internal variables
    integer				:: index_node,   index_node2
    integer				:: index_rhs,    index_rhs2
    integer				:: ijA_position, ijA_position2
    integer				:: ilarge
    logical				:: is_local
    real*8				:: mach1,  dmach1,  d2mach1_dTi, d2mach1_dTe
    real*8				:: mach_u, dmach_u, dmach_rho
    integer				:: ilarge_vv,   ilarge_vTi,   ilarge_vTe,   ilarge_vus
    integer				:: ilarge_vsvs, ilarge_vsTis, ilarge_vsTes, ilarge_vsTi, ilarge_vsTe
    integer				:: n_tor_local 
    
    n_tor_local = i_tor_max - i_tor_min +1
    ! --- Define node indices
    index_node  = node%index(1) 	    ! position of value
    index_node2 = node%index(side) 	    ! position of first deriative
    index_rhs   = n_tor_local*n_var*(index_node-1 ) + (k_Vpar-1)*n_tor_local + i_tor -i_tor_min + 1
    index_rhs2  = n_tor_local*n_var*(index_node2-1) + (k_Vpar-1)*n_tor_local + i_tor -i_tor_min + 1

    ! --- Define equations before MURGE and non-MURGE fork
    mach1       = - zbig / Btot * direction                     * sqrt(GAMMA*(Ti0 + Te0))	     
    dmach1      = - zbig / Btot * direction * 0.5d0  * GAMMA    / sqrt(GAMMA*(Ti0 + Te0))	     
    if (side .eq. 2) then
      d2mach1_dTi = + zbig / Btot * direction * 0.25d0 * GAMMA**2 / (GAMMA*(Ti0 + Te0))**(3/2) * Ti0_s
      d2mach1_dTe = + zbig / Btot * direction * 0.25d0 * GAMMA**2 / (GAMMA*(Ti0 + Te0))**(3/2) * Te0_s
      mach_u      = - zbig * U0_s * R**2 / ps0_s
      dmach_u     = - zbig        * R**2 / ps0_s
      !mach_u      = - zbig * (U0_s + tauIC*   Pi0_s     /rho0   ) * R**2 / ps0_s
      !dmach_u     = - zbig                                        * R**2 / ps0_s
      !dmach_rho   = - zbig * (     - tauIC*   Pi0_s     /rho0**2) * R**2 / ps0_s &
      !              - zbig * (     + tauIC*(Ti0_s+Ti0)  /rho0   ) * R**2 / ps0_s
    else
      d2mach1_dTi = + zbig / Btot * direction * 0.25d0 * GAMMA**2 / (GAMMA*(Ti0 + Te0))**(3/2) * Ti0_t
      d2mach1_dTe = + zbig / Btot * direction * 0.25d0 * GAMMA**2 / (GAMMA*(Ti0 + Te0))**(3/2) * Te0_t
      mach_u      = - zbig * U0_t * R**2 / ps0_t
      dmach_u     = - zbig        * R**2 / ps0_t
      !mach_u      = - zbig * (U0_t + tauIC*   Pi0_t     /rho0   ) * R**2 / ps0_t
      !dmach_u     = - zbig                                        * R**2 / ps0_t
      !dmach_rho   = - zbig * (     - tauIC*   Pi0_t     /rho0**2) * R**2 / ps0_t &
      !              - zbig * (     + tauIC*(Ti0_t+Ti0)  /rho0   ) * R**2 / ps0_t
    endif
    
    ! --- Condition on nodes
    if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then
      call locate_irn_jcn(index_node,index_node, index_min,index_max,ijA_position,ijA_index, ijA_size, irn_jcn)
      call locate_irn_jcn(index_node,index_node2,index_min,index_max,ijA_position2,ijA_index, ijA_size, irn_jcn)

      ilarge_vv  = ijA_position  - 1 + ((k_Vpar-1)*n_tor_local + i_tor-i_tor_min) * n_var*n_tor_local + & 
        (k_Vpar-1)*n_tor_local + i_tor -i_tor_min + 1
      ilarge_vTi = ijA_position  - 1 + ((k_Vpar-1)*n_tor_local + i_tor-i_tor_min) * n_var*n_tor_local + & 
        (k_Ti  -1)*n_tor_local + i_tor -i_tor_min + 1
      ilarge_vTe = ijA_position  - 1 + ((k_Vpar-1)*n_tor_local + i_tor-i_tor_min) * n_var*n_tor_local + & 
        (k_Te  -1)*n_tor_local + i_tor -i_tor_min + 1
      !ilarge_vus = ijA_position2 - 1 + ((k_Vpar-1)*n_tor + i_tor-1) * n_var*n_tor + (k_u   -1)*n_tor + i_tor

      irn(ilarge_vv)  = n_tor_local * n_var * (index_node-1) + (k_Vpar-1)*n_tor_local + i_tor -i_tor_min + 1
      jcn(ilarge_vv)  = n_tor_local * n_var * (index_node-1) + (k_Vpar-1)*n_tor_local + i_tor -i_tor_min + 1
      A_mat(ilarge_vv)    = zbig

      irn(ilarge_vTi) = n_tor_local * n_var * (index_node-1) + (k_Vpar-1)*n_tor_local + i_tor -i_tor_min + 1
      jcn(ilarge_vTi) = n_tor_local * n_var * (index_node-1) + (k_Ti  -1)*n_tor_local + i_tor -i_tor_min + 1
      A_mat(ilarge_vTi)   = dmach1

      irn(ilarge_vTe) = n_tor_local * n_var * (index_node-1) + (k_Vpar-1)*n_tor_local + i_tor -i_tor_min + 1
      jcn(ilarge_vTe) = n_tor_local * n_var * (index_node-1) + (k_Te  -1)*n_tor_local + i_tor -i_tor_min + 1
      A_mat(ilarge_vTe)   = dmach1

      !irn(ilarge_vus) = n_tor * n_var * (index_node -1) + (k_Vpar-1)*n_tor + i_tor
      !jcn(ilarge_vus) = n_tor * n_var * (index_node2-1) + (k_u   -1)*n_tor + i_tor
      !A_mat(ilarge_vus)   = dmach_u

      if (i_tor .eq. 1) then
        !RHS_loc(index_rhs) = - Zbig*Vpar0 - mach_u - mach1
        RHS_loc(index_rhs) = - Zbig*Vpar0 - mach1
      else
        RHS_loc(index_rhs) = 0.d0
      endif
    endif
    
    ! --- Condition between nodes (d/ds or d/dt)
    if ((index_node2 .ge. index_min) .and. (index_node2 .le. index_max)) then
      call locate_irn_jcn(index_node2,index_node, index_min,index_max,ijA_position,ijA_index, ijA_size, irn_jcn)
      call locate_irn_jcn(index_node2,index_node2,index_min,index_max,ijA_position2,ijA_index, ijA_size, irn_jcn)

      ilarge_vsvs  = ijA_position2 - 1 + ((k_Vpar-1)*n_tor_local + i_tor-i_tor_min) * n_var*n_tor_local + & 
        (k_Vpar-1)*n_tor_local + i_tor -i_tor_min + 1
      ilarge_vsTis = ijA_position2 - 1 + ((k_Vpar-1)*n_tor_local + i_tor-i_tor_min) * n_var*n_tor_local + & 
        (k_Ti  -1)*n_tor_local + i_tor -i_tor_min + 1
      ilarge_vsTes = ijA_position2 - 1 + ((k_Vpar-1)*n_tor_local + i_tor-i_tor_min) * n_var*n_tor_local + & 
        (k_Te  -1)*n_tor_local + i_tor -i_tor_min + 1
      ilarge_vsTi  = ijA_position  - 1 + ((k_Vpar-1)*n_tor_local + i_tor-i_tor_min) * n_var*n_tor_local + & 
        (k_Ti  -1)*n_tor_local + i_tor -i_tor_min + 1
      ilarge_vsTe  = ijA_position  - 1 + ((k_Vpar-1)*n_tor_local + i_tor-i_tor_min) * n_var*n_tor_local + & 
        (k_Te  -1)*n_tor_local + i_tor -i_tor_min + 1

      irn(ilarge_vsvs)  = n_tor_local * n_var * (index_node2-1) + (k_Vpar-1)*n_tor_local + i_tor -i_tor_min + 1
      jcn(ilarge_vsvs)  = n_tor_local * n_var * (index_node2-1) + (k_Vpar-1)*n_tor_local + i_tor -i_tor_min + 1
      A_mat(ilarge_vsvs)    = zbig

      irn(ilarge_vsTis) = n_tor_local * n_var * (index_node2-1) + (k_Vpar-1)*n_tor_local + i_tor -i_tor_min + 1
      jcn(ilarge_vsTis) = n_tor_local * n_var * (index_node2-1) + (k_Ti  -1)*n_tor_local + i_tor -i_tor_min + 1
      A_mat(ilarge_vsTis)   = dmach1

      irn(ilarge_vsTes) = n_tor_local * n_var * (index_node2-1) + (k_Vpar-1)*n_tor_local + i_tor -i_tor_min + 1
      jcn(ilarge_vsTes) = n_tor_local * n_var * (index_node2-1) + (k_Te  -1)*n_tor_local + i_tor -i_tor_min + 1
      A_mat(ilarge_vsTes)   = dmach1

      irn(ilarge_vsTi)  = n_tor_local * n_var * (index_node2-1) + (k_Vpar-1)*n_tor_local + i_tor -i_tor_min + 1
      jcn(ilarge_vsTi)  = n_tor_local * n_var * (index_node -1) + (k_Ti  -1)*n_tor_local + i_tor -i_tor_min + 1
      A_mat(ilarge_vsTi)    = d2mach1_dTi

      irn(ilarge_vsTe)  = n_tor_local * n_var * (index_node2-1) + (k_Vpar-1)*n_tor_local + i_tor -i_tor_min + 1
      jcn(ilarge_vsTe)  = n_tor_local * n_var * (index_node -1) + (k_Te  -1)*n_tor_local + i_tor -i_tor_min + 1
      A_mat(ilarge_vsTe)    = d2mach1_dTe

      if (i_tor .eq. 1) then
        if (side .eq. 2) then
          RHS_loc(index_rhs2) = - Zbig*Vpar0_s - dmach1 * (Ti0_s + Te0_s)
        else
          RHS_loc(index_rhs2) = - Zbig*Vpar0_t - dmach1 * (Ti0_t + Te0_t)
        endif
      else
        Rhs_loc(index_rhs2) = 0.d0
      endif
    endif
  
    return
  end subroutine apply_Mach1_BCs
  
  
end module mod_boundary_conditions
