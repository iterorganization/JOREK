program jorek2_connection_flux_aligned
  !-----------------------------------------------------------------------
  ! This routine is a modified version of jorek2_poincare. Field lines are traced around the simulated domain in order
  ! to determine the distance travelled before reaching the simulation boundary, or an arbitrary psi_n flux surface
  ! within the plasma.
  !
  ! Currently the diagnostic calculates the connection length and strike points on the set end boundary.
  ! 
  ! STEPS:
  ! 1.  Read input parameters from connect.nml and distribute start points among MPI tasks.
  ! 2.  Loop over start points
  !   3. Trace field lines around torus for a pre-set number of turns using a pre-set number of steps.
  !      4. Loop over max turns and steps per turn
  !        5. Perform step.
  !        6. Check if element boundary is crossed.
  !        7. If end boundary is crossed - signal to break and start new field line
  !      8. Record turn based data for poincares. 
  ! 9.  Write out poincare 
  ! 10. Write connection length
  ! 11. Write strike points
  !-----------------------------------------------------------------------
  use data_structure
  use phys_module
  use basis_at_gaussian
  use mod_element_rtree, only: populate_element_rtree
  use elements_nodes_neighbours
  use constants
  use mod_import_restart
  use mod_boundary
  use mod_neighbours
  use mod_interp
  use equil_info
  
  implicit none

  
  type (type_bnd_element_list)           :: bnd_elm_list    
  type (type_bnd_node_list)              :: bnd_node_list 
  
  ! --- Poincare data
  real*8,allocatable  :: rp(:), zp(:), R_all(:), Z_all(:), psin_all(:), C_all(:)                        ! arrays for position variables, and connection length for all lines
  real*8,allocatable  :: R_strike(:),  Z_strike(:), P_strike(:)                                         ! position of strike points
  real*8,allocatable  :: C_strike(:),  B_strike(:)                                                      ! connection length, boundary type at strike points
  real*8,allocatable  :: T0_strike(:), T_strike(:)                                                      ! temperature at start and end of fieldline
  real*8,allocatable  :: ZN0_strike(:), ZN_strike(:)                                                    ! density at start and end of fieldline
  real*8,allocatable  :: PS0_strike(:)                                                                  ! flux at starting point
  real*8,allocatable  :: in_domain_strike(:)                                                            ! Determine if field line is inside the domain
  real*8,allocatable  :: R_turn(:,:), Z_turn(:,:), C_turn(:,:), C_turn_tmp(:,:)                         ! position, and connection length of field line after each turn
  real*8,allocatable  :: T_turn(:,:), PSI_turn(:,:), ZN_turn(:,:)                                       ! physical parameters of field line after each turn
  
  ! --- Extra data
  integer   :: ntheta, npsin
  integer   :: my_id, ikeep, n_cpu, ierr, nsend, nrecv, ikeep0, i_line0
  integer   :: thetas_per_cpu, local_theta_start, local_theta_end
  integer   :: nnos, i_var, i_strike, i_strike0
  integer   :: i, j, iside_i, iside_j, ip, i_line, n_lines, i_tor, i_harm, i_var_psi, i_var_n, i_var_T, i_dir, k, m
  integer   :: i_elm, i_elm_start, ifail, i_phi, n_phi, i_turn, n_turns, i_elm_out, i_elm_prev, i_steps, n_turn_max(2)
  integer   :: v_s0_t0, v_s1_t0, v_s1_t1, v_s0_t1                                                        ! vertex indices for local element vertices
  integer   :: e_splus, e_sminus, e_tplus, e_tminus                                                      ! edge indices for local element neighbours
  real*8    :: delta_theta, delta_psin
  real*8    :: R_start, Z_start, phi_start
  real*8    :: R_line, Z_line, s_line, t_line, p_line
  real*8    :: s_mid, t_mid, p_mid
  real*8    :: s_ini, t_ini
  real*8    :: s_out, t_out
  real*8    :: R, R_s, R_t, R_st, R_ss, R_tt, R_in, R_keep, Rmid, Rmid_s, Rmid_t
  real*8    :: Z, Z_s, Z_t, Z_st, Z_ss, Z_tt, Z_in, Z_keep, Zmid, Zmid_s, Zmid_t
  real*8    :: P, P_s, P_t, P_st, P_ss, P_tt
  real*8    :: tol, psi_s, psi_t
  real*8    :: psin_range_min, psin_range_max
  real*8    :: Rmin, Rmax
  real*8    :: Zmin, Zmax
  real*8    :: psin_strike_bnd
  real*8    :: delta_phi, delta_phi_local, delta_phi_step, total_phi
  real*8    :: delta_s, small_delta_s
  real*8    :: delta_t, small_delta_t
  real*8    :: small_delta, dl2, total_length, length_max
  real*8    :: zl1, zl2, partial(2)
  real*8    :: psi_bnd, psi_bnd2
  integer   :: bnd_tmp, bnd_tmp_opp
  real*8    :: s_tmp, t_tmp
  real*8,allocatable  :: RZkeep(:,:),RhoThetakeep(:,:)
  integer   :: status(MPI_STATUS_SIZE)
  integer   :: count_lines, count_lines_tot
  real*8    :: C_average,   C_average_tot
  real*8    :: small_r, theta_pol
  logical, parameter  :: Rtheta_plot = .true.
  integer   :: n_points_max
  real*8    :: R_tmp,Z_tmp
  real*8    :: psi_tmp,psin_tmp, P0_s,P0_t,P0_st,P0_ss,P0_tt
  integer   :: ielm_tmp

  ! Data for find starting point routine
  real*8  :: theta_start, psin_start, previous_r_lower, previous_r_upper                                    

  ! -------------------------------------------------------------------------------------------------
  ! --- First part: Initialisation
  ! -------------------------------------------------------------------------------------------------
  namelist /connect_params/ n_turns, n_phi, ntheta, npsin, phi_start, tol, psin_range_min, psin_range_max, psin_strike_bnd

#define DEBUG

  ! --- MPI initilisation
  call MPI_INIT(IERR)
  call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)             ! id of each MPI proc
  call MPI_COMM_SIZE(MPI_COMM_WORLD, n_cpu, ierr)             ! number of MPI procs
  
  if (my_id .eq. 0) then
    write(*,*) '***************************************'
    write(*,*) '* JOREK2_connection_flux_aligned      *'
    write(*,*) '***************************************'
    write(*,*) 'WARNING: This program is still in development'
  endif
  
  ! --- Initilise data
  call initialise_basis                                     ! define the basis functions at the Gaussian points
  call initialise_parameters(my_id, "__NO_FILENAME__")
  do i_tor=1, n_tor
    mode(i_tor) = + int(i_tor / 2) * n_period
  enddo
  
  call broadcast_phys(my_id)                                ! physics parameters
  call import_restart(node_list, element_list, 'jorek_restart', rst_format, ierr)
  call boundary_from_grid(node_list, element_list, bnd_node_list, bnd_elm_list, .false.)

  ! --- Broadcast accross MPIs
!  call broadcast_elements(my_id, element_list)              ! elements
!  call broadcast_nodes(my_id, node_list)                    ! nodes
  call populate_element_rtree(node_list, element_list)       ! populate tree for all mpi tasks

  ! --- Define element neighbours
  allocate(element_neighbours(4,element_list%n_elements))
  element_neighbours = 0
  do i=1,element_list%n_elements
    do j=i+1,element_list%n_elements
      if (neighbours(node_list,element_list%element(i),element_list%element(j),iside_i,iside_j)) then
        element_neighbours(iside_i,i) = j
        element_neighbours(iside_j,j) = i
      endif
    enddo
  enddo
  
  ! --- Initialise constants
  i_var_psi = 1                             ! magnetic flux index
  i_var_n   = 1                             ! number density index
  i_var_T   = 6                             ! temperature index
  
  v_s0_t0   = 1                             ! the vertex and edge indices follow and anti-clockwise convention
  v_s1_t0   = 2                             ! around the element: https://www.jorek.eu/wiki/doku.php?id=grids
  v_s1_t1   = 3
  v_s0_t1   = 4
  e_splus   = 2
  e_sminus  = 4
  e_tplus   = 3
  e_tminus  = 1

  ! --- Preset parameters and read inputs from file if provided
  n_turns     = 500                         ! number of toroidal turns to follow a fieldline
  n_phi       = 1000                        ! number of steps per toroidal turn
  tol         = 1.d-6                       ! tolerance when stepping from element to element
  phi_start   = PI/4.                       ! starting toroidal angle
  ntheta      = 500                         ! number of poloidal starting points on each flux surface
  npsin       = 25                          ! number of radial starting points on each flux surface
  psin_range_min = 0.d0                     ! Minimum psin value
  psin_range_max = 0.995                    ! Maximym psin value - should be smaller than psin_strike_bnd
  psin_strike_bnd = 1.0                     ! Boundary location outside of which field lines are no longer traced and strike points are recorded
  open(42, file='connect.nml', action='read', status='old', iostat=ierr)
  if ( ierr == 0 ) then
    if (my_id .eq. 0 ) write(*,*) 'Reading parameters from connect.nml namelist.'
    read(42,connect_params)
    close(42)
  end if
  n_lines = 2 * ntheta * npsin              ! number of starting points to allocate
  if (psin_range_max .gt. 0.999) then
    write(*, *) 'Maximum psi normalised cannot be greater than 0.999. The input value is: ', psin_range_max
    stop
  end if

  ! --- Allocate data
  allocate(R_strike(n_lines),Z_strike(n_lines),P_strike(n_lines),C_strike(n_lines),B_strike(n_lines), in_domain_strike(n_lines))
  allocate(T0_strike(n_lines),T_strike(n_lines),ZN0_strike(n_lines),ZN_strike(n_lines),PS0_strike(n_lines))
  allocate(R_all(n_lines),Z_all(n_lines),psin_all(n_lines),C_all(n_lines))
  allocate(R_turn(n_turns+1,2),Z_turn(n_turns+1,2),C_turn(n_turns+1,2),C_turn_tmp(n_turns+1,2))
  allocate(T_turn(n_turns+1,2),PSI_turn(n_turns+1,2),ZN_turn(n_turns+1,2))
  n_points_max = 10000000
  allocate(RZkeep(2,n_points_max),RhoThetakeep(2,n_points_max))

  ! --- Initialise allocated data
  R_all     = 0.d0; Z_all     = 0.d0; psin_all=0.d0; C_all = 0.d0;  
  R_strike  = 0.d0; Z_strike  = 0.d0; P_strike  = 0.d0;  C_strike   = 0.d0; in_domain_strike = 0.d0
  T0_strike = 0.d0; T_strike  = 0.d0; ZN0_strike = 0.d0; ZN_strike  = 0.d0; PS0_strike = 0.d0
  R_turn    = 0.d0; Z_turn    = 0.d0; C_turn  = 0.d0;  C_turn_tmp = 0.d0
  
  ! --- Get domain limits
  Rmin = 1.d20; Rmax = -1.d20; Zmin = 1.d20; Zmax=-1.d20
  do i=1,node_list%n_nodes
    Rmin = min(Rmin,node_list%node(i)%x(1,1,1))
    Rmax = max(Rmax,node_list%node(i)%x(1,1,1))
    Zmin = min(Zmin,node_list%node(i)%x(1,1,2))
    Zmax = max(Zmax,node_list%node(i)%x(1,1,2))
  enddo
  
  ! --- find x-point(s)
  if (xpoint) then
    psi_bnd  = ES%psi_xpoint(1)
    psi_bnd2 = ES%psi_xpoint(2)
    if( ES%active_xpoint .eq. UPPER_XPOINT ) then
      psi_bnd  = ES%psi_xpoint(2)
      psi_bnd2 = ES%psi_xpoint(1)
    endif
    if (xcase .eq. LOWER_XPOINT) psi_bnd2 = psi_bnd
  else
    call find_limiter(my_id, node_list, element_list, bnd_elm_list, psi_bnd, R_start, Z_start)
    psi_bnd2 = psi_bnd
  endif
    
  ! --- The elements our local MPI is looking at
  ! --- Field lines are separated poloidally to ensure that the computational time
  ! --- should be evenly distributed among MPI tasks
  if (mod(ntheta, n_cpu) .ne. 0) then
    write(*,*) "ERROR: n_theta must be exactly divisible by n_cpu! Otherwise poloidal distribution of field lines is non-uniform.", ntheta, n_cpu
    stop 
  endif 
  thetas_per_cpu = ntheta / n_cpu
  local_theta_start = my_id * thetas_per_cpu + 1
  local_theta_end = min(ntheta, (my_id + 1) * thetas_per_cpu)
  if (my_id .eq. 0) write(*, *) thetas_per_cpu
  write(*, *) my_id, local_theta_start, local_theta_end
  
  ! --- Some info print outs
  delta_theta = 2.d0 * PI / ntheta
  delta_psin = (psin_range_max - psin_range_min) / npsin
  write(*,*) ' Local MPI poloidal angles (mpi_id, theta_start, theta_end) : ', my_id, local_theta_start * delta_theta, local_theta_end * delta_theta
  
  ! --- Initialise counters
  i_line      = 0
  i_strike    = 0
  ikeep       = 0
  count_lines = 0
  C_average   = 0.d0
  
  ! -------------------------------------------------------------------------------------------------
  ! --- Second part: Loop over starting points and store connection lengths
  ! -------------------------------------------------------------------------------------------------
  ! --- Loop over poloidal locations
  do i = local_theta_start, local_theta_end
    theta_start = i * delta_theta

    ! --- Print progress
    write(*,'(A,4i6)')' Progress: MPI_id, element, start, end :', my_id, i, local_theta_start, local_theta_end
    
    ! Initialise minor radius bounds for locating starting point
    previous_r_lower = 0.0
    previous_r_upper = 0.05

    ! Loop over radial locations
    do j=1, npsin
      ! Get starting R, Z, phi - phi is set arbitrarily
      psin_start = psin_range_min + j * delta_psin
      call find_starting_element(i_elm_start, psin_start, ES%psi_axis, psi_bnd, theta_start, ES%R_axis, ES%Z_axis, ES%Z_xpoint, s_ini, t_ini, previous_r_lower, previous_r_upper, R_start, Z_start)

      ! ----------------------
      ! --- We have a new line
      ! ----------------------
      i_line = i_line + 1
    
      ! --- Initialise temporary data (for the line)
      R_turn     = 0.d0
      Z_turn     = 0.d0
      C_turn_tmp = 0.d0

      ! --- Do both directions of field line
      do i_dir = -1,1,2
         
        ! --- The toroidal step (with direction)
        delta_phi = 2.d0 * PI * float(i_dir) / float(n_period*n_phi)

        ! --- Initialise data before start
        i_elm = i_elm_start
        s_line = s_ini
        t_line = t_ini
        total_length = 0.d0
        total_phi    = 0.d0

        ! --- Record the first point
        R_all(i_line) = R_start
        Z_all(i_line) = Z_start
        call var_value(i_elm, i_var_psi, s_line, t_line, phi_start, psin_all(i_line))
        psin_all(i_line) = (psin_all(i_line) - ES%psi_axis) / (psi_bnd - ES%psi_axis)  

        R_turn(1,(i_dir+1)/2+1) = R_start
        Z_turn(1,(i_dir+1)/2+1) = Z_start
        C_turn(1,(i_dir+1)/2+1) = 0.d0
    
        ! --- Get variables at start
        call var_value(i_elm, i_var_T, s_line,t_line,phi_start, T_turn  (1,(i_dir+1)/2+1) )
        PSI_turn(1,(i_dir+1)/2+1) = psin_all(i_line)
        call var_value(i_elm, i_var_n, s_line,t_line,phi_start, ZN_turn (1,(i_dir+1)/2+1) )
        
        ! --- Set and store starting location
        R_line = R_start
        Z_line = Z_start
        p_line = phi_start

        ! --- We assume this line will give a strike on the target
        i_strike = i_strike + 1
        ZN0_strike(i_strike) = ZN_turn (1,(i_dir+1)/2+1)
        T0_strike (i_strike) = T_turn  (1,(i_dir+1)/2+1)
        PS0_strike(i_strike) = PSI_turn(1,(i_dir+1)/2+1)
    
        ! --- Loop over toroidal turns
        do i_turn = 1, n_turns
          ! --- Record the maximum number of turns
          n_turn_max((i_dir+1)/2+1) = i_turn
    
          ! --- Loop over toroidal steps
          do i_phi=1,n_phi
            ! --- Loop within one element
            delta_phi_local = 0.d0
            i_steps = 0
            do while ((abs(delta_phi_local) .lt. abs(delta_phi)) .and. (i_steps .lt. 100) )
              ! --- Count element steps (lt 100)
              i_steps = i_steps + 1

              ! --- Find the new position after element step
              delta_phi_step = delta_phi - delta_phi_local
              call step(i_elm,s_line,t_line,p_line,delta_phi_step, delta_s,delta_t,R,Z,R_s,R_t,Z_s,Z_t)
          
              ! --- Take another step from middle point (ie. 1.5*step)
              s_mid = s_line + 0.5d0 * delta_s
              t_mid = t_line + 0.5d0 * delta_t
              p_mid = p_line + 0.5d0 * delta_phi_step
              call step(i_elm,s_mid,t_mid,p_mid,delta_phi_step,delta_s,delta_t,Rmid,Zmid,Rmid_s,Rmid_t,Zmid_s,Zmid_t)
    
              ! --- Step to element boundary, not beyond (we do smaller steps to approach boundary)
              small_delta_s = 1.d0
              if (s_line + delta_s .gt. 1.d0) small_delta_s = (1.d0 - s_line)/delta_s
              if (s_line + delta_s .lt. 0.d0) small_delta_s = abs(s_line/delta_s)
    
              ! --- Step to element boundary, not beyond (we do smaller steps to approach boundary)
              small_delta_t = 1.d0
              if (t_line + delta_t .gt. 1.d0) small_delta_t = (1.d0 - t_line)/delta_t
              if (t_line + delta_t .lt. 0.d0) small_delta_t = abs(t_line/delta_t)
    
              ! --- Do we require a smaller step (< 1.0) 
              small_delta = min(small_delta_s, small_delta_t)
    
              if (small_delta .lt. 1.d0)  then       ! this step is crossing the boundary
                s_mid = s_line + 0.5d0 * small_delta * delta_s
                t_mid = t_line + 0.5d0 * small_delta * delta_t
                p_mid = p_line + 0.5d0 * small_delta * delta_phi_step
          
                call step(i_elm,s_mid,t_mid,p_mid,delta_phi_step,delta_s,delta_t,Rmid,Zmid,Rmid_s,Rmid_t,Zmid_s,Zmid_t)
    
                if (small_delta_s .lt. small_delta_t) then
    
                  if (s_line + delta_s .gt. 1.d0) then     ! crossing boundary 2 or 4 at s=1
                    s_line = 1.0
                    call find_new_element(.true., i_elm, i_elm_prev, &
                                          R_strike(i_strike), Z_strike(i_strike), P_strike(i_strike), C_strike(i_strike), &
                                          T_strike(i_strike), ZN_strike(i_strike), B_strike(i_strike), &
                                          Rmid, Rmid_s, Rmid_t, Zmid_s, Zmid_t, &
                                          total_length, &
                                          s_line, t_line, p_line, small_delta, delta_s, delta_t, delta_phi_step, &
                                          e_splus, e_sminus, v_s1_t0, v_s1_t1, &
                                          dl2)
                  elseif (s_line + delta_s .lt. 0.d0) then ! crossing boundary 2 or 4 at s=0
                    s_line = 0.d0
                    call find_new_element(.true., i_elm, i_elm_prev, &
                                          R_strike(i_strike), Z_strike(i_strike), P_strike(i_strike), C_strike(i_strike), &
                                          T_strike(i_strike), ZN_strike(i_strike), B_strike(i_strike), &
                                          Rmid, Rmid_s, Rmid_t, Zmid_s, Zmid_t, &
                                          total_length, &
                                          s_line, t_line, p_line, small_delta, delta_s, delta_t, delta_phi_step, &
                                          e_sminus, e_splus, v_s0_t1, v_s0_t0, &
                                          dl2)
                  endif
                else
                  if (t_line + delta_t .gt. 1.d0) then  ! crossing boundary 1 or 3 at t=1
                    t_line = 1.0
                    call find_new_element(.false., i_elm, i_elm_prev, &
                                          R_strike(i_strike), Z_strike(i_strike), P_strike(i_strike), C_strike(i_strike), &
                                          T_strike(i_strike), ZN_strike(i_strike), B_strike(i_strike), &
                                          Rmid, Rmid_s, Rmid_t, Zmid_s, Zmid_t, &
                                          total_length, &
                                          s_line, t_line, p_line, small_delta, delta_s, delta_t, delta_phi_step, &
                                          e_tplus, e_tminus, v_s1_t1, v_s0_t1, &
                                          dl2)
                  elseif (t_line + delta_t .lt. 0.d0) then  ! crossing boundary 1 or 3 at t=0
                    t_line = 0.d0
                    call find_new_element(.false., i_elm, i_elm_prev, &
                                          R_strike(i_strike), Z_strike(i_strike), P_strike(i_strike), C_strike(i_strike), &
                                          T_strike(i_strike), ZN_strike(i_strike), B_strike(i_strike), &
                                          Rmid, Rmid_s, Rmid_t, Zmid_s, Zmid_t, &
                                          total_length, &
                                          s_line, t_line, p_line, small_delta, delta_s, delta_t, delta_phi_step, &
                                          e_tminus, e_tplus, v_s0_t0, v_s0_t1, &
                                          dl2)
                  endif
                endif
              else  ! this step remains within the element
                s_line = s_line + delta_s
                t_line = t_line + delta_t
                p_line = p_line + delta_phi_step
          
                dl2 = (Rmid_s**2 + Zmid_s**2)*delta_s**2 + (Rmid_t**2 + Zmid_t**2)*delta_t**2 + Rmid**2 * delta_phi_step**2
#ifdef DEBUG
                if (dl2 .ne. dl2) then
                  write(*, *) 'NaN detected in dl2!', Rmid_s, Zmid_s, delta_s, Rmid_t, Zmid_t, delta_t, Rmid, delta_phi_step
                  stop
                end if
#endif

                small_delta = 1.d0
    
                call interp_RZ(node_list,element_list,i_elm,s_line,t_line,R_in,R_s,R_t,R_st,R_ss,R_tt,Z_in,Z_s,Z_t,Z_st,Z_ss,Z_tt)
    
              endif

              ! --- Reset the toroial step
              delta_phi_local = delta_phi_local + small_delta * delta_phi_step
    
              ! --- Record total lengths
              total_length = total_length + sqrt(abs(dl2))
              total_phi    = total_phi    + small_delta * delta_phi_step

              ! --- Exit if we stepped out of domain
              if (i_elm .eq. 0) exit
              
              call interp(node_list,element_list,i_elm,1,1,s_line,t_line,psi_tmp,P0_s,P0_t,P0_st,P0_ss,P0_tt)
              psin_tmp = (psi_tmp - ES%psi_axis) / (psi_bnd - ES%psi_axis)
              if (psin_tmp .gt. psin_strike_bnd) exit
            enddo  ! end of loop over steps within one element
    
            ! --- Exit if we stepped out of domain
            if ((i_elm .eq. 0) .or. (psin_tmp .gt. psin_strike_bnd)) exit
          enddo    ! end of a 2Pi turn (or before if end of open field line)
    
          ! --- Exit if we stepped out of domain
          if ((i_elm .eq. 0)) exit
    
          ! --- Record Poincare point
          call interp_RZ(node_list,element_list,i_elm,s_line,t_line,R_in,R_s,R_t,R_st,R_ss,R_tt,Z_in,Z_s,Z_t,Z_st,Z_ss,Z_tt)
          R_turn(i_turn+1,(i_dir+1)/2+1) = R_in
          Z_turn(i_turn+1,(i_dir+1)/2+1) = Z_in
          C_turn_tmp(i_turn+1,(i_dir+1)/2+1) = total_length
          call var_value(i_elm, i_var_T, s_line,t_line,p_line, T_turn  (i_turn+1,(i_dir+1)/2+1) )
          call var_value(i_elm, i_var_psi, s_line,t_line,p_line, PSI_turn(i_turn+1,(i_dir+1)/2+1) )
          call var_value(i_elm, i_var_n, s_line,t_line,p_line, ZN_turn (i_turn+1,(i_dir+1)/2+1) )
          
          ! Finish line if it has struck the psin based boundary
          if (psin_tmp .gt. psin_strike_bnd) then
            R_strike(i_strike) = R_in
            Z_strike(i_strike) = Z_in
            P_strike(i_strike) = p_line
            C_strike(i_strike) = total_length
            call var_value(i_elm,i_var_T,s_line,t_line,p_line,T_strike(i_strike))
            call var_value(i_elm,i_var_n,s_line,t_line,p_line,ZN_strike(i_strike))
            exit
          endif
        enddo  ! end of loop over toroidal turns
    
        ! --- Field line still in domain, after n_turn turns
        if ((i_elm .ne. 0) .and. (psin_tmp .lt. psin_strike_bnd)) then
          R_strike(i_strike) = R_in
          Z_strike(i_strike) = Z_in
          P_strike(i_strike) = p_line
          C_strike(i_strike) = total_length     
          in_domain_strike(i_strike) = 1
          call var_value(i_elm,i_var_T,s_line,t_line,p_line,T_strike(i_strike))
          call var_value(i_elm,i_var_n,s_line,t_line,p_line,ZN_strike(i_strike))
        endif
    
        ! --- Record connection length
        if (i_dir .eq. -1) then  
          C_all(i_line)   = total_length
          partial(1)      = total_length
        else
          C_all(i_line)   = C_all(i_line)+total_length!min(C_all(i_line),total_length)
          partial(2)      = total_length
        endif
    
      enddo  ! end of two directions

      ! ------------------------------------------------------
      ! --- Record turn based data for poincare plots --------
      ! ------------------------------------------------------
      ! --- Reverse connection length (not from plasma to position, but from target to position)
      do i_dir=1,2
        do i_turn = 1, n_turn_max(i_dir)
          C_turn(i_turn,i_dir) = partial(i_dir) - c_turn_tmp(i_turn,i_dir)
        enddo
      enddo

      ! --- Keep only field lines starting inside the plasma region
      if ( ( ( (xcase .eq. 1) .and. (Z_turn(1,1) .gt. ES%Z_xpoint(1)) ) &
          .or.( (xcase .eq. 2) .and. (Z_turn(1,1) .lt. ES%Z_xpoint(2)) ) &
          .or.( (xcase .eq. 3) .and. (Z_turn(1,1) .gt. ES%Z_xpoint(1)) .and. (Z_turn(1,1) .lt. ES%Z_xpoint(2)) ) ) ) then

        ! --- Compute averaged connection length
        count_lines = count_lines + 1
        C_average = C_average + C_all(i_line)
  
        ! --- Do both directions
        do i_dir=1,2
          ! --- Do all turns
          do i_turn=1,n_turn_max(i_dir)+1
            if (R_turn(i_turn,i_dir) .gt. 0.d0) then

              ! --- Connection lengths in both directions
              ikeep = ikeep + 1
              small_r   = sqrt( (R_turn(i_turn,i_dir)-ES%R_axis)**2 + (Z_turn(i_turn,i_dir)-ES%Z_axis)**2 )
              theta_pol = atan2(Z_turn(i_turn,i_dir)-ES%Z_axis,R_turn(i_turn,i_dir)-ES%R_axis)
              if (theta_pol .lt. 0.d0) theta_pol = theta_pol + 2.d0*PI
              
              if (ikeep .le. n_points_max) then
                call find_RZ(node_list,element_list,R_turn(i_turn,i_dir),Z_turn(i_turn,i_dir),R_tmp,Z_tmp,ielm_tmp,s_tmp,t_tmp,ifail)
                call interp(node_list,element_list,ielm_tmp,1,1,s_tmp,t_tmp,psi_tmp,P0_s,P0_t,P0_st,P0_ss,P0_tt)
                RhoThetakeep(1,ikeep)      = (psi_tmp-ES%psi_axis)/(psi_bnd-ES%psi_axis)!small_r
                RhoThetakeep(2,ikeep)      = theta_pol / (2.d0*PI)
                RZkeep(1,ikeep)            = R_turn(i_turn,i_dir)
                RZkeep(2,ikeep)            = Z_turn(i_turn,i_dir)
              else
                write(*,*) 'Warning! Exceeded maximum number of points',n_points_max
              endif
            endif
          enddo  
        enddo ! end recording in both directions
      endif
    enddo ! end over loop over radial coordinate
  enddo ! end of loop over poloidal coordinate
  write(*,*)'Finished loop for MPI task: ', my_id
  
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  
  ! --- Print averaged connection length
  call MPI_Reduce(count_lines,count_lines_tot,1,MPI_INTEGER,         MPI_SUM,0,MPI_COMM_WORLD,ierr)
  call MPI_Reduce(C_average,  C_average_tot,  1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  if (my_id .eq. 0) then
    if (count_lines_tot .gt. 0) then
      C_average_tot = C_average_tot / real(count_lines_tot)
      write(*,*)'Averaged connection length:',C_average_tot
    endif
  endif

  ! -------------------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------------------
  ! --- Third part: write poincare plot to VTK file
  ! -------------------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------------------
  write(*,*)'Writing poincare files'
  
  ! --- The local number of points (for mpi_0)
  ikeep0  = ikeep
  call MPI_Reduce(ikeep,nnos,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  if (my_id .eq. 0) write(*,*) ' number of points : ',nnos
  
  ! --- Open file and write headers
  open(21,file='poinc_R-Z.dat',status='replace')
  write(21,*) '#  R                 Z'
  open(22,file='poinc_rho-theta.dat',status='replace')
  write(22,*) '# rho=sqrt(psi_n)'
  write(22,*) '# psi_n=(psi - ES%psi_axis)/(psi_bnd - ES%psi_axis)'
  write(22,*) '#'
  write(22,*) '#  rho               theta'
  
  ! --- Write points for local MPI (id=0)
  if (my_id .eq. 0) then
    write(21,'(2e18.8)') (RZkeep(1,i),RZkeep(2,i), i=1,ikeep0 )
    write(22,'(2e18.8)') (RhoThetakeep(1,i),RhoThetakeep(2,i), i=1,ikeep0 )
    write(21,*)
    write(21,*)
    write(22,*)
    write(22,*)
  endif
  
  ! --- Write points for all other MPIs
  if (my_id .eq. 0) then
    ! --- If this is mpi_0, we receive data from the other MPIs and print it
    do j=1,n_cpu-1
      call mpi_recv(ikeep,1, MPI_INTEGER, j, j, MPI_COMM_WORLD, status, ierr)
      if (ikeep .gt. 0) then
        nrecv = 2*ikeep
        call mpi_recv(RZkeep,nrecv, MPI_DOUBLE_PRECISION, j, j, MPI_COMM_WORLD, status, ierr)
        call mpi_recv(RhoThetakeep,nrecv, MPI_DOUBLE_PRECISION, j, j, MPI_COMM_WORLD, status, ierr)
        write(21,'(2e18.8)') (RZkeep(1,i),RZkeep(2,i), i=1,ikeep0 )
        write(22,'(2e18.8)') (RhoThetakeep(1,i),RhoThetakeep(2,i), i=1,ikeep0 )
        write(21,*)
        write(21,*)
        write(22,*)
        write(22,*)
      endif
    enddo
  else
    ! --- If this is not mpi_0, we send data to the main MPI 0
    call mpi_send(ikeep, 1, MPI_INTEGER, 0, my_id, MPI_COMM_WORLD, ierr)
    if (ikeep .gt. 0) then
      nsend = 2*ikeep
      call mpi_send(RZkeep, nsend,MPI_DOUBLE_PRECISION, 0, my_id, MPI_COMM_WORLD, ierr)
      call mpi_send(RhoThetakeep, nsend,MPI_DOUBLE_PRECISION, 0, my_id, MPI_COMM_WORLD, ierr)
    endif
  endif

  ! -------------------------------------------------------------------------------------------------
  ! --- Fourth part: write connection lengths to file
  ! -------------------------------------------------------------------------------------------------
  if (my_id .eq. 0) write(*,*) 'Writing lines to connections.txt'
  open(23, file='connections.txt', action='write', status='replace', form='formatted', recl=400)
  
  ! --- The local number of lines (for mpi_0)
  i_line0 = i_line
  call MPI_Reduce(i_line,nnos,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  if (my_id .eq. 0) write(*,*) ' number of lines : ',nnos

  ! --- Write points for local MPI (id=0)
  if (my_id .eq. 0) then
    write(23,*) "# R    Z    Psi_N    Connection_length"
    write(23,'(4e16.8)') ( (/ R_all(i), Z_all(i), psin_all(i), C_all(i) /),i=1,i_line0)
  endif
  
  ! --- Write points for all other MPIs
  if (my_id .eq. 0) then
    ! --- If this is mpi_0, we receive data from the other MPIs and print it
    do j=1,n_cpu-1
      call mpi_recv(i_line,1, MPI_INTEGER, j, j, MPI_COMM_WORLD, status, ierr)
      if (i_line .gt. 0) then
        nrecv = i_line
        call mpi_recv(R_all,nrecv, MPI_DOUBLE_PRECISION, j, j, MPI_COMM_WORLD, status, ierr)
        call mpi_recv(Z_all,nrecv, MPI_DOUBLE_PRECISION, j, j, MPI_COMM_WORLD, status, ierr)
        call mpi_recv(psin_all,nrecv, MPI_DOUBLE_PRECISION, j, j, MPI_COMM_WORLD, status, ierr)
        call mpi_recv(C_all,nrecv, MPI_DOUBLE_PRECISION, j, j, MPI_COMM_WORLD, status, ierr)
        write(23,'(4e16.8)') ( (/R_all(i), Z_all(i), psin_all(i), C_all(i) /),i=1,i_line)
      endif
      write(*, *) 'Received MPI task: ', j
    enddo
  else
    ! --- If this is not mpi_0, we send data to the main MPI 0
    call mpi_send(i_line, 1, MPI_INTEGER, 0, my_id, MPI_COMM_WORLD, ierr)
    if (i_line .gt. 0) then
      nsend = i_line
      call mpi_send(R_all, nsend, MPI_DOUBLE_PRECISION, 0, my_id, MPI_COMM_WORLD, ierr)
      call mpi_send(Z_all, nsend, MPI_DOUBLE_PRECISION, 0, my_id, MPI_COMM_WORLD, ierr)
      call mpi_send(psin_all, nsend, MPI_DOUBLE_PRECISION, 0, my_id, MPI_COMM_WORLD, ierr)
      call mpi_send(C_all, nsend, MPI_DOUBLE_PRECISION, 0, my_id, MPI_COMM_WORLD, ierr)
    endif
  endif
  close(23)

  deallocate(RZkeep,RhoThetakeep)

  ! -------------------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------------------
  ! --- Fifth part: write strike points to file
  ! -------------------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------------------
  ! --- Open text file
  write(*,*)'Writing strikes.txt...'
  open(23,file='strikes.txt',status='replace')
  
  ! --- The local number of strikes (for mpi_0)
  i_strike0 = i_strike
  call MPI_Reduce(i_strike,nnos,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  if (my_id .eq. 0) write(*,*) ' number of points : ',nnos
  
  ! --- Write points for local MPI (id=0)
  if (my_id .eq. 0) then
    write(23,*) "# R    Z    phi    Psi_N    T    Connection_length    in_domain"
    write(23,'(7f22.8)') ( (/ R_strike(i), Z_strike(i), P_strike(i), Ps0_strike(i), T0_strike(i), C_strike(i), in_domain_strike(i) /),i=1,i_strike0)
  endif
  
  ! --- Write points for all other MPIs
  if (my_id .eq. 0) then
    ! --- If this is mpi_0, we receive data from the other MPIs and print it
    do j=1,n_cpu-1
      call mpi_recv(i_strike,1, MPI_INTEGER, j, j, MPI_COMM_WORLD, status, ierr)
      if (i_strike .gt. 0) then
        nrecv = i_strike
        call mpi_recv(R_strike,nrecv, MPI_DOUBLE_PRECISION, j, j, MPI_COMM_WORLD, status, ierr)
        call mpi_recv(Z_strike,nrecv, MPI_DOUBLE_PRECISION, j, j, MPI_COMM_WORLD, status, ierr)
        call mpi_recv(P_strike,nrecv, MPI_DOUBLE_PRECISION, j, j, MPI_COMM_WORLD, status, ierr)
        call mpi_recv(Ps0_strike,nrecv, MPI_DOUBLE_PRECISION, j, j, MPI_COMM_WORLD, status, ierr)
        call mpi_recv(T0_strike,nrecv, MPI_DOUBLE_PRECISION, j, j, MPI_COMM_WORLD, status, ierr)
        call mpi_recv(C_strike,nrecv, MPI_DOUBLE_PRECISION, j, j, MPI_COMM_WORLD, status, ierr)
        call mpi_recv(in_domain_strike,nrecv, MPI_DOUBLE_PRECISION, j, j, MPI_COMM_WORLD, status, ierr)
        write(23,'(7f22.8)') ( (/ R_strike(i), Z_strike(i), P_strike(i), Ps0_strike(i), T0_strike(i), C_strike(i), in_domain_strike(i) /),i=1,i_strike)
      endif
    enddo
  else
    ! --- If this is not mpi_0, we send data to the main MPI 0
    call mpi_send(i_strike, 1, MPI_INTEGER, 0, my_id, MPI_COMM_WORLD, ierr)
    if (i_strike .gt. 0) then
      nsend = i_strike
      call mpi_send(R_strike, nsend, MPI_DOUBLE_PRECISION, 0, my_id, MPI_COMM_WORLD, ierr)
      call mpi_send(Z_strike, nsend, MPI_DOUBLE_PRECISION, 0, my_id, MPI_COMM_WORLD, ierr)
      call mpi_send(P_strike, nsend, MPI_DOUBLE_PRECISION, 0, my_id, MPI_COMM_WORLD, ierr)
      call mpi_send(Ps0_strike, nsend, MPI_DOUBLE_PRECISION, 0, my_id, MPI_COMM_WORLD, ierr)
      call mpi_send(T0_strike, nsend, MPI_DOUBLE_PRECISION, 0, my_id, MPI_COMM_WORLD, ierr)
      call mpi_send(C_strike, nsend, MPI_DOUBLE_PRECISION, 0, my_id, MPI_COMM_WORLD, ierr)
      call mpi_send(in_domain_strike, nsend, MPI_DOUBLE_PRECISION, 0, my_id, MPI_COMM_WORLD, ierr)
    endif
  endif
  
  ! Clean and exit
  close(23)
  call MPI_FINALIZE(IERR)                                ! clean up MPI

end program jorek2_connection_flux_aligned

! --------------------------------
! ----- BEGIN SUBROUTINES --------
! --------------------------------
subroutine find_starting_element(i_elm, &
                                 psin_start, psi_axis, psi_bnd, theta_start, &
                                 R_axis, Z_axis, Z_xpoint, &
                                 s_ini, t_ini, &
                                 previous_r_lower, previous_r_upper, &
                                 R_start, Z_start)
  use mod_parameters
  use elements_nodes_neighbours
  use phys_module
  use mod_interp 

  implicit none
  ! Determine the starting point of a new field line on a particular flux surface
  !
  ! 1. Hunting algorithm for upper and lower bounds to psi normalised location, along radial coordinate in terms of psi normalised
  ! 2. Bisection algorithm for location of radial coordinate in terms of psi normalised

  ! Input parameters
  integer, intent(inout) :: i_elm                                                                      ! Output element index
  real*8, intent(in)     :: psin_start, psi_axis, psi_bnd, theta_start                                 ! Required starting psi normalised and poloidal coordinate
  real*8, intent(in)     :: R_axis, Z_axis                                                             ! Plasma axis
  real*8, intent(in)     :: Z_xpoint(2)                                                                ! xpoint Z location
  real*8, intent(inout)  :: s_ini, t_ini                                                               ! Starting s, t values for line
  real*8, intent(inout)  :: previous_r_lower, previous_r_upper                                         ! Previous minor radius bounds of field line from last call
  real*8, intent(inout)  :: R_start, Z_start                                                           ! Output R-Z location of field line

  ! local variables
  real*8                 :: r_lower, r_upper, r_tmp                                                    ! Upper and lower bounds for search
  real*8                 :: BigR_tmp, BigR_upper, BigR_lower, Z_upper, Z_lower, Z_tmp                  ! R-Z Coordinates of upper and lower bound
  real*8                 :: psin_diff, delta_psi_bnd                                                   ! delta psi between upper and lower boundaries and normalising constant for psi
  real*8                 :: psin_upper, psin_lower, psi_tmp, P0_s,P0_t,P0_st,P0_ss,P0_tt               ! psin_upper, psin_lower and dummy variables for interpolation
  logical                :: psin_lower_gt_psin_start, psin_upper_lt_psin_start                         ! truth statements for upper and lower bounds
  logical                :: z_beyond_xpoint                                                            ! logical to test if current bounds are in xpoint region
  integer                :: ifail=0                                                                    ! flag for algorithm failure
  real*8                 :: iter_tol=1.d-8                                                             ! Tolerance for bisection
  integer                :: num_iter=0                                                                 ! Number of iterations
  integer                :: sub_iter=0                                                                 ! Number of sub-iterations in algorithm
  integer                :: max_iterations=100                                                         ! maximum iterations in hunting and bisection algorithms
  
  delta_psi_bnd = psi_bnd - psi_axis
  if (previous_r_lower .ge. previous_r_upper) then
    write(*, *) 'Inconsistent starting values for finding new element: ', previous_r_lower, previous_r_upper
    stop
  end if

  ! Get initial bounds
  r_lower = previous_r_lower
  BigR_tmp = R_axis + r_lower * cos(theta_start)
  Z_tmp = Z_axis + r_lower * sin(theta_start)
  call find_RZ(node_list, element_list, BigR_tmp, Z_tmp, BigR_lower, Z_lower, i_elm, s_ini, t_ini, ifail)
  if (ifail .ne. 0) then
    write(*, *) 'Lower find_RZ failed'
    stop
  endif
  call interp(node_list,element_list,i_elm,1,1,s_ini,t_ini,psin_lower,P0_s,P0_t,P0_st,P0_ss,P0_tt)
  psin_lower = (psin_lower - psi_axis) / delta_psi_bnd
  
  r_upper = previous_r_upper
  BigR_tmp = R_axis + r_upper * cos(theta_start)
  Z_tmp = Z_axis + r_upper * sin(theta_start)
  call find_RZ(node_list, element_list, BigR_tmp, Z_tmp, BigR_upper, Z_upper, i_elm, s_ini, t_ini, ifail)
  if (ifail .ne. 0) then
    write(*, *) 'Upper find_RZ failed'
    stop
  endif
  call interp(node_list,element_list,i_elm,1,1,s_ini,t_ini,psin_upper,P0_s,P0_t,P0_st,P0_ss,P0_tt)
  psin_upper = (psin_upper - psi_axis) / delta_psi_bnd
  
  ! Peform hunting algorithm
  ! Check boundaries are consistent
  if ((psin_lower .gt. psin_upper) .or. (r_lower .gt. r_upper)) then
    write(*,*)  "ERROR: Inconsistent input boundaries for hunting algorithm: ", psin_lower, psin_upper, r_lower, r_upper 
    stop
  endif

  num_iter = 0
  psin_lower_gt_psin_start = (psin_lower .gt. psin_start)
  psin_upper_lt_psin_start = (psin_upper .lt. psin_start)
  do while ((psin_lower_gt_psin_start .or. psin_upper_lt_psin_start)) 
    num_iter = num_iter + 1
    if (num_iter .gt. max_iterations) then
      write(*, *) "ERROR: Maximum iterations exceeded in Hunting algorithm"
      stop
    endif
    
    ! Find new boundaries
    if (psin_lower_gt_psin_start) then
      r_tmp = r_lower
      r_lower = max(0.0, r_lower - 2 * abs(r_upper - r_lower))
      psin_upper = psin_lower
      r_upper = r_tmp
      BigR_tmp = R_axis + r_lower * cos(theta_start)
      Z_tmp = Z_axis + r_lower * sin(theta_start)

      ! Interpolate lower bound psi
      call find_RZ(node_list, element_list, BigR_tmp, Z_tmp, BigR_lower, Z_lower, i_elm, s_ini, t_ini, ifail)
      if (ifail .ne. 0) then
        write(*, *) 'Lower find_RZ failed'
        stop
      endif
      call interp(node_list,element_list,i_elm,1,1,s_ini,t_ini,psin_lower,P0_s,P0_t,P0_st,P0_ss,P0_tt)
      psin_lower = (psin_lower - psi_axis) / delta_psi_bnd
    else if (psin_upper_lt_psin_start) then
      r_tmp = r_upper
      r_upper = r_upper + 2 * abs(r_upper - r_lower)
      psin_lower = psin_upper
      r_lower = r_tmp
      BigR_tmp = R_axis + r_upper * cos(theta_start)
      Z_tmp = Z_axis + r_upper * sin(theta_start)

      ! interpolate upper bound psi
      call find_RZ(node_list, element_list, BigR_tmp, Z_tmp, BigR_upper, Z_upper, i_elm, s_ini, t_ini, ifail)

      ! Correct r_upper if new point is outside of the domain, or in the x-point region, where bisections will fail
      z_beyond_xpoint = (xpoint .and. (((xcase .eq. 1) .and. (Z_upper .lt. Z_xpoint(1)) ) &
          .or.( (xcase .eq. 2) .and. (Z_upper .gt. Z_xpoint(2)) ) &
          .or.( (xcase .eq. 3) .and. ((Z_upper .lt. Z_xpoint(1)) .or. (Z_upper .gt. Z_xpoint(2))))))
      sub_iter = 0
      do while ((ifail .ne. 0) .or. z_beyond_xpoint)
        if (sub_iter .gt. 50) then 
          write(*,*) 'ERROR: Hunting algorithm failed to find point inside the domain'
          stop
        endif
        sub_iter = sub_iter + 1
        r_upper = 0.5 * (r_lower + r_upper)
        BigR_tmp = R_axis + r_upper * cos(theta_start)
        Z_tmp = Z_axis + r_upper * sin(theta_start)
        call find_RZ(node_list, element_list, BigR_tmp, Z_tmp, BigR_upper, Z_upper, i_elm, s_ini, t_ini, ifail)
        z_beyond_xpoint = (xpoint .and. (((xcase .eq. 1) .and. (Z_upper .lt. Z_xpoint(1)) ) &
          .or.( (xcase .eq. 2) .and. (Z_upper .gt. Z_xpoint(2)) ) &
          .or.( (xcase .eq. 3) .and. ((Z_upper .lt. Z_xpoint(1)) .or. (Z_upper .gt. Z_xpoint(2))))))
      enddo
      call interp(node_list,element_list,i_elm,1,1,s_ini,t_ini,psin_upper,P0_s,P0_t,P0_st,P0_ss,P0_tt)
      psin_upper = (psin_upper - psi_axis) / delta_psi_bnd
    endif

    psin_lower_gt_psin_start = (psin_lower .gt. psin_start)
    psin_upper_lt_psin_start = (psin_upper .lt. psin_start)
  enddo
  previous_r_lower = r_lower
  previous_r_upper = r_upper

  ! Perform bisection
  num_iter = 0
  psin_diff = psin_upper - psin_lower
  do while ((psin_diff .gt. iter_tol) .and. num_iter .lt. max_iterations) 
    num_iter = num_iter + 1
    if (num_iter .gt. max_iterations) then
      write(*, *) "ERROR: Maximum iterations exceeded in Bisection algorithm"
      stop
    endif

    r_tmp = 0.5 * (r_upper + r_lower)
    BigR_tmp = R_axis + r_tmp * cos(theta_start)
    Z_tmp = Z_axis + r_tmp * sin(theta_start)
    call find_RZ(node_list, element_list, BigR_tmp, Z_tmp, R_start, Z_start, i_elm, s_ini, t_ini, ifail)
    call interp(node_list,element_list,i_elm,1,1,s_ini,t_ini,psi_tmp,P0_s,P0_t,P0_st,P0_ss,P0_tt)
    psi_tmp = (psi_tmp - psi_axis) / delta_psi_bnd
    
    ! Set new boundaries
    if (psin_start .lt. psi_tmp) then
      r_upper = r_tmp; psin_upper = psi_tmp;
    else
      r_lower = r_tmp; psin_lower = psi_tmp;
    endif

    psin_diff = psin_upper - psin_lower;
  enddo
end

subroutine find_new_element(s_bnd, i_elm, i_elm_prev, &
                            R_strike, Z_strike, P_strike, C_strike, T_strike, ZN_strike, B_strike, &
                            R_mid, R_mid_s, R_mid_t, Z_mid_s, Z_mid_t, &
                            total_length, &
                            s_line, t_line, p_line, small_delta, delta_s, delta_t, delta_phi_step, &
                            element_to_neighbour_idx, neighbour_to_element_idx, vertex_1, vertex_2, &
                            dl2)
  use mod_parameters
  use elements_nodes_neighbours
  use phys_module
  use mod_interp

  implicit none
  ! Determines the new element after performing a step along a field line that has crossed a boundary
  ! 
  ! 1. Find boundary point in original element 
  ! 2. Get neighbour element
  ! 3. Confirm boundary point is the same for both elements
  ! 4. If new element is outside domain, record strike point 
  
  ! Input parameters
  logical, intent(in)     :: s_bnd                                                                    ! Determine whether crossing is from the s coordinate  
  integer, intent(inout)  :: i_elm, i_elm_prev                                                        ! Element index of current element and previous element
  real*8, intent(inout)   :: total_length                                                             ! Total connection length
  real*8, intent(inout)   :: R_strike, Z_strike, P_strike                                             ! Strike locations
  real*8, intent(inout)   :: C_strike, T_strike, ZN_strike, B_strike                                  ! Strike connection length, and physical state
  real*8, intent(in)      :: R_mid, R_mid_s, R_mid_t, Z_mid_s, Z_mid_t                                ! R and Z derivatives for current half step
  integer, intent(in)     :: element_to_neighbour_idx, neighbour_to_element_idx                       ! Indices of current element to neighbour and expected neighbour to current element without orientation change
  integer, intent(in)     :: vertex_1, vertex_2                                                       ! Vertex indices for boundary nodes
  real*8, intent(inout)   :: s_line, t_line, p_line, small_delta, delta_s, delta_t, delta_phi_step    ! Current values of s, t, and phi and deltas from current step
  real*8, intent(out)     :: dl2                                                                      ! Change in distance from step

  ! Local variables
  real*8                  :: R_in, R_out, R_s, R_t, R_st, R_ss, R_tt                                  ! Variables for interpolation of R
  real*8                  :: Z_in, Z_out, Z_s, Z_t, Z_st, Z_ss, Z_tt                                  ! Variables for interpolation of Z
  integer                 :: i_elm_tmp                                                                ! Temporary element to check neighbour element orientation is consistent
  integer                 :: inode1, inode2                                                           ! Nodes corresponding to boundary edge
  integer                 :: i_var_n=5, i_var_T=6;                                                    ! Indices of number density and temperature

  if (s_bnd) then
    t_line = t_line + small_delta * delta_t
  else
    s_line = s_line + small_delta * delta_s
  endif
  p_line = p_line + small_delta * delta_phi_step

  dl2 = (R_mid_s**2 + Z_mid_s**2)*delta_s**2 + (R_mid_t**2 + Z_mid_t**2)*delta_t**2 + R_mid**2 * delta_phi_step**2
  dl2 = dl2 * small_delta**2

  call interp_RZ(node_list,element_list,i_elm,s_line,t_line,R_in,R_s,R_t,R_st,R_ss,R_tt,Z_in,Z_s,Z_t,Z_st,Z_ss,Z_tt)

  i_elm_prev = i_elm
  i_elm = element_neighbours(element_to_neighbour_idx,i_elm_prev)

  if (i_elm .ne. 0) then
    i_elm_tmp  = element_neighbours(neighbour_to_element_idx,i_elm)
    if (i_elm_prev .ne. i_elm_tmp) write(*,*) ' WARNING : CHANGE OF ORIENTATION'

    ! Switch boundary coordinate from 1.0 to 0.0 or vice versa
    if (s_bnd) then
      s_line = abs(1.d0 - s_line)
    else
      t_line = abs(1.d0 - t_line)
    endif

    ! Check if R and Z from both elements are close
    call interp_RZ(node_list,element_list,i_elm,s_line,t_line,R_out,R_s,R_t,R_st,R_ss,R_tt,Z_out,Z_s,Z_t,Z_st,Z_ss,Z_tt)
    if ( (abs(R_in - R_out) .gt. 1.d-8) .or. (abs(Z_in - Z_out) .gt. 1.d-8)) &
      write(*,'(A,2i6,4f12.4)') ' error in element change ',i_elm_prev,i_elm,R_in,R_out,Z_in,Z_out

  else ! crossing an outer boundary
    R_strike = R_in
    Z_strike = Z_in
    P_strike = p_line
    C_strike = total_length + sqrt(abs(dl2))

    call var_value(i_elm_prev,i_var_T,s_line,t_line,p_line,T_strike)
    call var_value(i_elm_prev,i_var_n,s_line,t_line,p_line,ZN_strike)
    
    inode1 = element_list%element(i_elm_prev)%vertex(vertex_1)
    inode2 = element_list%element(i_elm_prev)%vertex(vertex_2)

    if ((node_list%node(inode1)%boundary .ne. 0) .and. (node_list%node(inode2)%boundary .ne. 0)) then
      B_strike = min(node_list%node(inode1)%boundary,node_list%node(inode2)%boundary)
    else
      write(*,*) 'error : leaving domain but not at correct boundary!',inode1,inode2,R_in,Z_in
    endif
  endif
end

! Take step along the traced field line
subroutine step(i_elm,s_in,t_in,p_in,delta_p,delta_s,delta_t,R,Z,R_s,R_t,Z_s,Z_t)
  use mod_parameters
  use elements_nodes_neighbours
  use phys_module
  use mod_interp
  
  implicit none
  
  integer :: i_var_psi, i_elm, i_tor, i_harm
  
  real*8 :: s_in, t_in, p_in, delta_p, delta_s, delta_t
  real*8 :: R_out, Z_out, Rs_out, Rt_out, Zs_out, Zt_out
  real*8 :: R,R_s,R_t,R_st,R_ss,R_tt,Z,Z_s,Z_t,Z_st,Z_ss,Z_tt
  real*8 :: Pcos,Pcos_s,Pcos_t,Pcos_st,Pcos_ss,Pcos_tt, Psin,Psin_s,Psin_t,Psin_st,Psin_ss,Psin_tt
  real*8 :: P0,P0_s,P0_t,P0_st,P0_ss,P0_tt, psi_s, psi_t, Zjac
  
  i_var_psi = 1
  
  call interp_RZ(node_list,element_list,i_elm,s_in,t_in,R,R_s,R_t,Z,Z_s,Z_t)
  
  Zjac = (R_s * Z_t - R_t * Z_s)
#ifdef DEBUG
  if (Zjac .eq. 0.0) then
    write(*, *) 'Jacobian is zero!', i_elm, s_in, t_in, R, R_s, R_t, Z, Z_s, Z_t
    stop
  end if
#endif

  call interp(node_list,element_list,i_elm,i_var_psi,1,s_in,t_in,P0,P0_s,P0_t,P0_st,P0_ss,P0_tt)
  
  psi_s = P0_s 
  psi_t = P0_t 
  
  do i_tor = 1, (n_tor-1)/2
    i_harm = 2*i_tor
  
    call interp(node_list,element_list,i_elm,i_var_psi,i_harm,s_in,t_in,Pcos,Pcos_s,Pcos_t,Pcos_st,Pcos_ss,Pcos_tt)
    psi_s = psi_s + Pcos_s * cos(mode(i_harm)*p_in)
    psi_t = psi_t + Pcos_t * cos(mode(i_harm)*p_in)
  
    call interp(node_list,element_list,i_elm,i_var_psi,i_harm+1,s_in,t_in,Psin,Psin_s,Psin_t,Psin_st,Psin_ss,Psin_tt)
    psi_s = psi_s + Psin_s * sin(mode(i_harm+1)*p_in)
    psi_t = psi_t + Psin_t * sin(mode(i_harm+1)*p_in)
  enddo
  
  delta_s =   psi_t * R / (Zjac * F0) * delta_p
  delta_t = - psi_s * R / (Zjac * F0) * delta_p
  
  return
end

subroutine var_value(i_elm,i_var,s_in,t_in,p_in,value_out)
  use mod_parameters
  use elements_nodes_neighbours
  use phys_module
  use mod_interp
  
  implicit none
  
  integer :: i_var, i_elm, i_tor, i_harm
  
  real*8 :: s_in, t_in, p_in
  real*8 :: Pcos,Pcos_s,Pcos_t,Pcos_st,Pcos_ss,Pcos_tt, Psin,Psin_s,Psin_t,Psin_st,Psin_ss,Psin_tt
  real*8 :: P0,P0_s,P0_t,P0_st,P0_ss,P0_tt
  real*8 :: value_out
  
  call interp(node_list,element_list,i_elm,i_var,1,s_in,t_in,P0,P0_s,P0_t,P0_st,P0_ss,P0_tt)
  value_out = P0
  
  do i_tor = 1, (n_tor-1)/2
    i_harm = 2*i_tor
  
    call interp(node_list,element_list,i_elm,i_var,i_harm,s_in,t_in,Pcos,Pcos_s,Pcos_t,Pcos_st,Pcos_ss,Pcos_tt)
    value_out = value_out + Pcos * cos(mode(i_harm)*p_in)
  
    call interp(node_list,element_list,i_elm,i_var,i_harm+1,s_in,t_in,Psin,Psin_s,Psin_t,Psin_st,Psin_ss,Psin_tt)
    value_out = value_out + Psin * sin(mode(i_harm+1)*p_in)
  enddo

  return
end
