program jorek2_connection_fmhd
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  use data_structure
  use phys_module
  use basis_at_gaussian
  use elements_nodes_neighbours
  use constants
  use mod_import_restart
  use mod_neighbours
  use mod_interp
  use equil_info
  use mod_element_rtree
  use mpi
  
  implicit none
  
  ! --- Poincare data
  real*8,allocatable    :: rp(:), zp(:), R_all(:), Z_all(:)
  real*4,allocatable    :: C_all_4(:)
  real*8,allocatable    :: R_turn(:,:), Z_turn(:,:), C_turn(:,:)
  real*8,allocatable    :: T_turn(:,:), PSI_turn(:,:), ZN_turn(:,:)
  
  ! --- Extra data
  integer               :: my_id, ikeep, n_cpu, ierr, nsend, nrecv, ikeep0, inode1, inode2, i_line0
  integer               :: elm_start, elm_end, elm_delta, local_elm_start, local_elm_end
  integer               :: i_elm_half, cpu_fraction
  integer               :: nnos, n_scalars, ivtk, i_var
  integer               :: i, j, iside_i, iside_j, i_line, n_lines, i_tor, i_harm, i_var_psi, i_dir, k, m
  integer               :: i_elm, ifail, i_phi, n_phi, i_turn, n_turns, i_elm_out, n_turn_max(2)
  real*8                :: R_start, Z_start, P_start
  real*8                :: R_line, Z_line, s_line, t_line, p_line
  real*8                :: xyz(3), xyz_prev(3)
  real*8                :: s_out, t_out
  real*8                :: R, R_s, R_t, R_st, R_ss, R_tt, R_in, R_out
  real*8                :: Z, Z_s, Z_t, Z_st, Z_ss, Z_tt, Z_in, Z_out
  real*8                :: P, P_s, P_t, P_st, P_ss, P_tt
  real*8                :: psi, psi_s, psi_t, psi_st, psi_ss, psi_tt, psi_R, psi_Z
  real*8                :: Rmin, Rmax
  real*8                :: Zmin, Zmax
  real*8                :: delta_phi
  real*8                :: zl1, zl2, partial(2)
  real*8                :: value_out
  real*4,allocatable    :: RZkeep(:,:),scalars(:,:)
  real*4                :: ZERO
  integer               :: status(MPI_STATUS_SIZE)
  character             :: buffer*80, lf*1, str1*12, str2*12
  character*12, allocatable :: scalar_names(:)
  integer               :: count_lines, count_lines_tot
  real*8                :: C_average, C_average_tot, C_min, C_all, total_length
  real*8                :: small_r, theta_pol
  logical               :: Rtheta_plot
  logical               :: psitheta_plot
  logical               :: ignore_SOL
  integer               :: n_points_max
  real*8                :: P0_s,P0_t,P0_st,P0_ss,P0_tt
  integer               :: n_points_start
  real*8                :: Rstart_min, Rstart_max
  real*8                :: Zstart_min, Zstart_max
  real*8                :: phi_start, progress
  real*8                :: xjac, dpsi_dR, dpsi_dZ, phi, R_half, Z_half
  real*8                :: Fprof,Fprof_s,Fprof_t,Fprof_st,Fprof_ss,Fprof_tt
  real*8                :: P0
  real*8                :: AR0,AR0_s,AR0_t,AR0_st,AR0_ss,AR0_tt,AR0_R,AR0_Z,AR0_p
  real*8                :: AZ0,AZ0_s,AZ0_t,AZ0_st,AZ0_ss,AZ0_tt,AZ0_R,AZ0_Z,AZ0_p
  real*8                :: A30,A30_s,A30_t,A30_st,A30_ss,A30_tt,A30_R,A30_Z,A30_p
  real*8                :: BR, BZ, Bp
  real*8                :: HHZ(n_tor), HHZ_p(n_tor)

  namelist /poincare_params/ &
      n_points_start,        &! Number of points along a line to start field-lines
      Rstart_min,            &! R-start of line along which field lines will be started, default: R_axis
      Zstart_min,            &! Z-start of line along which field lines will be started, default: Z_axis
      Rstart_max,            &! R-end   of line along which field lines will be started, default: R_max of grid
      Zstart_max,            &! Z-end   of line along which field lines will be started, default: Z_axis
      n_turns,               &! Number of toroidal turns per each field line (not, going in two directions), default: 500
      n_phi,                 &! Toroidal discretisation, number of steps per toroidal turn, default: 1000
      phi_start,             &! Toroidal angle were Poincare plot is made, default: 0.0
      n_points_max,          &! Maximum number of points in Poincare plot, default: 10000000
      ignore_SOL,            &! Ignore lines that started outside the separatrix, default: false
      Rtheta_plot,           &! Plot with (x-axis,y-axis) = (amin,theta), default, false
      psitheta_plot           ! Plot with (x-axis,y-axis) = (psi_n,theta), default: false
 
  
  if (my_id .eq. 0) write(*,*) '***************************************'
  if (my_id .eq. 0) write(*,*) '* JOREK2_poincare                     *'
  if (my_id .eq. 0) write(*,*) '***************************************'
  
  
  ! -------------------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------------------
  ! --- First part: Initialisation
  ! -------------------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------------------
  
  ! --- MPI initilisation
  call MPI_INIT(IERR)
  !required=MPI_THREAD_MULTIPLE
  !call MPI_Init_thread(required,provided,StatInfo)
  call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)      ! id of each MPI proc
  call MPI_COMM_SIZE(MPI_COMM_WORLD, n_cpu, ierr)      ! number of MPI procs
  
  ! --- Initilise data
  call initialise_parameters(my_id, "__NO_FILENAME__")
  do i_tor=1, n_tor
    mode(i_tor) = + int(i_tor / 2) * n_period
  enddo
  call broadcast_elements(my_id, element_list)                ! elements
  call broadcast_nodes(my_id, node_list)                      ! nodes
  call populate_element_rtree(node_list, element_list)
  call broadcast_phys(my_id)                                  ! physics parameters
  call broadcast_equil_state(my_id)                           ! equil_state
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  call import_restart(node_list, element_list, 'jorek_restart', rst_format, ierr, .true.)
  call initialise_basis                                       ! define the basis functions at the Gaussian points
  
  ! --- Broadcast accross MPIs
  call broadcast_elements(my_id, element_list)                ! elements
  call broadcast_nodes(my_id, node_list)                      ! nodes
  call broadcast_phys(my_id)                                  ! physics parameters
  
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
  
  ! --- Initialise variables
  n_scalars = 4
  ZERO      = 0.
  i_var_psi = 1                                 ! the index of the magnetic flux variable
  allocate(scalar_names(n_scalars))
  scalar_names  = (/ 'length_tot  ','length_min  ','psi_start   ','T_start     '/)
  
  ! --- Get domain limits
  Rmin = 1.d20; Rmax = -1.d20; Zmin = 1.d20; Zmax=-1.d20
  do i=1,node_list%n_nodes
    Rmin = min(Rmin,node_list%node(i)%x(1,1,1))
    Rmax = max(Rmax,node_list%node(i)%x(1,1,1))
    Zmin = min(Zmin,node_list%node(i)%x(1,1,2))
    Zmax = max(Zmax,node_list%node(i)%x(1,1,2))
  enddo
  
  ! --- Initialise toroidal variables
  n_points_start         = 100       ! Number of points along a line to start field-lines
  Rstart_min             = ES%R_axis ! R-start of line along which field lines will be started, default: R_axis
  Zstart_min             = ES%Z_axis ! Z-start of line along which field lines will be started, default: Z_axis
  Rstart_max             = Rmax      ! R-end   of line along which field lines will be started, default: R_max of grid
  Zstart_max             = ES%Z_axis ! Z-end   of line along which field lines will be started, default: Z_axis
  n_turns                = 500       ! Number of toroidal turns per each field line (not, going in two directions), default: 500
  n_phi                  = 1000      ! Toroidal discretisation, number of steps per toroidal turn, default: 1000
  phi_start              = 0.0       ! Toroidal angle were Poincare plot is made, default: 0.0
  n_points_max           = 10000000  ! Maximum number of points in Poincare plot, default: 10000000
  ignore_SOL             = .false.   ! Ignore lines that started outside the separatrix, default: false
  Rtheta_plot            = .false.   ! Plot with (x-axis,y-axis) = (amin,theta), default, false
  psitheta_plot          = .false.   ! Plot with (x-axis,y-axis) = (psi_n,theta), default: false
    
  ! --- Read parameters from namelist file 'vtk.nml' if it exists
  open(42, file='poincare.nml', action='read', status='old', iostat=ierr)
  if ( ierr == 0 ) then
    if (my_id .eq. 0) write(*,*) 'Reading parameters from poincare.nml namelist.'
    read(42,poincare_params)
    close(42)
  end if

  if (my_id .eq. 0) then
    write(*,*)
    write(*,*) '-----------'
    write(*,*) 'Parameters:'
    write(*,*) '-----------'
    write(*,*) 'n_points_start  :',  n_points_start     
    write(*,*) 'Rstart_min      :',  Rstart_min         
    write(*,*) 'Zstart_min      :',  Zstart_min         
    write(*,*) 'Rstart_max      :',  Rstart_max         
    write(*,*) 'Zstart_max      :',  Zstart_max         
    write(*,*) 'n_turns         :',  n_turns            
    write(*,*) 'n_phi           :',  n_phi              
    write(*,*) 'phi_start       :',  phi_start          
    write(*,*) 'n_points_max    :',  n_points_max       
    write(*,*) 'ignore_SOL      :',  ignore_SOL         
    write(*,*) 'Rtheta_plot     :',  Rtheta_plot        
    write(*,*) 'psitheta_plot   :',  psitheta_plot     
  endif 
    
  ! --- Number of lines
  n_lines = n_points_start

  ! --- Allocate data
  allocate(R_all(n_lines),Z_all(n_lines))
  allocate(R_turn(n_turns+1,2),Z_turn(n_turns+1,2),C_turn(n_turns+1,2))
  allocate(T_turn(n_turns+1,2),PSI_turn(n_turns+1,2),ZN_turn(n_turns+1,2))
  allocate(RZkeep(2,n_points_max),scalars(n_points_max,n_scalars))

  ! --- Initialise allocated data
  R_all     = 0.d0; Z_all     = 0.d0; 
  R_turn    = 0.d0; Z_turn    = 0.d0; C_turn    = 0.d0
  
  ! --- The elements we are looking at (note element indexing starts at axis going gradually outwards)
  elm_start = 1
  elm_end   = n_points_start
  
  ! --- The elements our local MPI is looking at
  elm_delta = (elm_end - elm_start) / n_cpu + 1
  local_elm_start = elm_start + my_id*elm_delta
  local_elm_end   = min(elm_end,elm_start+(my_id+1)*elm_delta - 1)
  
  ! --- Initialise counters
  i_line      = 0
  ikeep       = 0
  count_lines = 0
  C_average   = 0.d0
  
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  if (my_id .eq. 0) write(*,*)'Starting main loop...'

  ! -------------------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------------------
  ! --- Second part: Loop over starting points and get poincare points
  ! -------------------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------------------
  
  ! --- Start going over all local elements
  do i = local_elm_start, local_elm_end
    
    ! --- Print progress
    if (local_elm_end .ne. local_elm_start) then
      progress = real(i-local_elm_start)/real(local_elm_end-local_elm_start)*100.d0
    else
      progress = 100.0
    endif
    if (mod(i,10).eq.0) &
      write(*,'(A,5i6,A)')' Progress: MPI_id, element, start, end :', my_id, i, local_elm_start, local_elm_end,int(progress),'%'
    
    i_line = i_line + 1
    
    ! --- Initialise temporary data (for the line)
    R_turn     = 0.d0
    Z_turn     = 0.d0
    
    ! --- Do both directions of field line
    do i_dir = -1,1,2
    
      ! --- The toroidal step (with direction)
      delta_phi = 2.d0 * PI * float(i_dir) / float(n_period*n_phi)
      
    
      R_start = Rstart_min + float(i-1)/float(n_points_start-1) * (Rstart_max - Rstart_min)
      Z_start = Zstart_min + float(i-1)/float(n_points_start-1) * (Zstart_max - Zstart_min)
      P_start = phi_start
      R_turn(1,(i_dir+1)/2+1) = R_start
      Z_turn(1,(i_dir+1)/2+1) = Z_start
      C_turn(1,(i_dir+1)/2+1) = 0.d0
      R_line = R_start
      Z_line = Z_start
      phi = phi_start
      total_length = 0.d0

      call find_RZ(node_list,element_list,R_start,Z_start,R_out,Z_out,i_elm,s_out,t_out,ifail)
      if ( ifail .ne. 0 ) cycle

      call var_value(i_elm, var_T,   s_out,t_out,P_start, T_turn  (1,(i_dir+1)/2+1) )
      call var_value(i_elm, var_A3,  s_out,t_out,P_start, PSI_turn(1,(i_dir+1)/2+1) )
      call var_value(i_elm, var_rho, s_out,t_out,P_start, ZN_turn (1,(i_dir+1)/2+1) )

      xyz(1) = R_start * cos(phi)
      xyz(2) = Z_start
      xyz(3) = R_start * sin(phi)
      xyz_prev = xyz
      
      ! --- Loop over toroidal turns
      do i_turn = 1, n_turns
        
        ! --- Record the maximum number of turns
        n_turn_max((i_dir+1)/2+1) = i_turn
      
        ! --- Perform the field line tracing.
        do j = 1, n_phi
        
          ! --- toroidal functions
          HHZ  (1) = 1.d0
          HHZ_p(1) = 0.d0
          do i_tor=1,(n_tor-1)/2
            HHZ  (2*i_tor)   = +cos(mode(2*i_tor)   * phi )
            HHZ  (2*i_tor+1) = +sin(mode(2*i_tor+1) * phi )
            HHZ_p(2*i_tor)   = -sin(mode(2*i_tor)   * phi ) * float(mode(2*i_tor))
            HHZ_p(2*i_tor+1) = +cos(mode(2*i_tor+1) * phi ) * float(mode(2*i_tor+1))
          enddo
          
          ! -----------------
          ! --- half step ---
          ! -----------------
          call find_RZ(node_list,element_list,R_line,Z_line,R_out,Z_out,i_elm,s_out,t_out,ifail)
          if ( ifail .eq. 0 ) then
            ! --- RZ-coords
            call interp_RZ(node_list,element_list,i_elm,s_out,t_out,R,R_s,R_t,Z,Z_s,Z_t)
            xjac = (R_s * Z_t - R_t * Z_s)
#ifdef fullmhd
            ! --- B-variables
            AR0_p  = 0.d0 ; AR0_R  = 0.d0 ; AR0_Z  = 0.d0
            AZ0_p  = 0.d0 ; AZ0_R  = 0.d0 ; AZ0_Z  = 0.d0
            A30_p  = 0.d0 ; A30_R  = 0.d0 ; A30_Z  = 0.d0
            do i_tor = 1, n_tor
              call interp(node_list,element_list,i_elm,var_AR, i_tor,s_out,t_out,AR0,AR0_s,AR0_t,AR0_st,AR0_ss,AR0_tt)
              call interp(node_list,element_list,i_elm,var_AZ, i_tor,s_out,t_out,AZ0,AZ0_s,AZ0_t,AZ0_st,AZ0_ss,AZ0_tt)
              call interp(node_list,element_list,i_elm,var_A3, i_tor,s_out,t_out,A30,A30_s,A30_t,A30_st,A30_ss,A30_tt)
              AR0_p  = AR0_p  + AR0 * HHZ_p(i_tor)
              AZ0_p  = AZ0_p  + AZ0 * HHZ_p(i_tor)
              A30_p  = A30_p  + A30 * HHZ_p(i_tor)
              if ((xjac .gt. 1.d-6)) then  ! avoid the axis
                AR0_R  = AR0_R  + (   Z_t * AR0_s - Z_s * AR0_t ) / xjac * HHZ(i_tor)
                AR0_Z  = AR0_Z  + ( - R_t * AR0_s + R_s * AR0_t ) / xjac * HHZ(i_tor)
                AZ0_R  = AZ0_R  + (   Z_t * AZ0_s - Z_s * AZ0_t ) / xjac * HHZ(i_tor)
                AZ0_Z  = AZ0_Z  + ( - R_t * AZ0_s + R_s * AZ0_t ) / xjac * HHZ(i_tor)
                A30_R  = A30_R  + (   Z_t * A30_s - Z_s * A30_t ) / xjac * HHZ(i_tor)
                A30_Z  = A30_Z  + ( - R_t * A30_s + R_s * A30_t ) / xjac * HHZ(i_tor)
              endif
            enddo

            ! --- Magnetic field
            call interp(node_list,element_list,i_elm,710,1,s_out,t_out,Fprof,Fprof_s,Fprof_t,Fprof_st,Fprof_ss,Fprof_tt)
            BR = ( A30_Z - AZ0_p )/ R
            BZ = ( AR0_p - A30_R )/ R
            Bp = ( AZ0_R - AR0_Z ) + Fprof / R
! reduced-MHD
#else
            psi_R  = 0.d0 ; psi_Z  = 0.d0
            do i_tor = 1, n_tor
              call interp(node_list,element_list,i_elm,var_psi, i_tor,s_out,t_out,psi,psi_s,psi_t,psi_st,psi_ss,psi_tt)
              if ((xjac .gt. 1.d-6)) then  ! avoid the axis
                psi_R  = psi_R  + (   Z_t * psi_s - Z_s * psi_t ) / xjac * HHZ(i_tor)
                psi_Z  = psi_Z  + ( - R_t * psi_s + R_s * psi_t ) / xjac * HHZ(i_tor)
              endif
            enddo
            BR = + psi_Z / R
            BZ = - psi_R / R
            Bp = F0 / R
#endif
            
            R_half = R_line + R*delta_phi/2. / Bp * BR
            Z_half = Z_line + R*delta_phi/2. / Bp * BZ
          else
            cycle
          endif
          
          ! -----------------
          ! --- full step ---
          ! -----------------
          call find_RZ(node_list,element_list,R_half,Z_half,R_out,Z_out,i_elm,s_out,t_out,ifail)
          if ( ifail .eq. 0 ) then
            ! --- RZ-coords
            call interp_RZ(node_list,element_list,i_elm,s_out,t_out,R,R_s,R_t,Z,Z_s,Z_t)
            xjac = (R_s * Z_t - R_t * Z_s)
#ifdef fullmhd
            ! --- B-variables
            AR0_p  = 0.d0 ; AR0_R  = 0.d0 ; AR0_Z  = 0.d0
            AZ0_p  = 0.d0 ; AZ0_R  = 0.d0 ; AZ0_Z  = 0.d0
            A30_p  = 0.d0 ; A30_R  = 0.d0 ; A30_Z  = 0.d0
            do i_tor = 1, n_tor
              call interp(node_list,element_list,i_elm,var_AR, i_tor,s_out,t_out,AR0,AR0_s,AR0_t,AR0_st,AR0_ss,AR0_tt)
              call interp(node_list,element_list,i_elm,var_AZ, i_tor,s_out,t_out,AZ0,AZ0_s,AZ0_t,AZ0_st,AZ0_ss,AZ0_tt)
              call interp(node_list,element_list,i_elm,var_A3, i_tor,s_out,t_out,A30,A30_s,A30_t,A30_st,A30_ss,A30_tt)
              AR0_p  = AR0_p  + AR0 * HHZ_p(i_tor)
              AZ0_p  = AZ0_p  + AZ0 * HHZ_p(i_tor)
              A30_p  = A30_p  + A30 * HHZ_p(i_tor)
              if ((xjac .gt. 1.d-6)) then  ! avoid the axis
                AR0_R  = AR0_R  + (   Z_t * AR0_s - Z_s * AR0_t ) / xjac * HHZ(i_tor)
                AR0_Z  = AR0_Z  + ( - R_t * AR0_s + R_s * AR0_t ) / xjac * HHZ(i_tor)
                AZ0_R  = AZ0_R  + (   Z_t * AZ0_s - Z_s * AZ0_t ) / xjac * HHZ(i_tor)
                AZ0_Z  = AZ0_Z  + ( - R_t * AZ0_s + R_s * AZ0_t ) / xjac * HHZ(i_tor)
                A30_R  = A30_R  + (   Z_t * A30_s - Z_s * A30_t ) / xjac * HHZ(i_tor)
                A30_Z  = A30_Z  + ( - R_t * A30_s + R_s * A30_t ) / xjac * HHZ(i_tor)
              endif
            enddo

            ! --- Magnetic field
            call interp(node_list,element_list,i_elm,710,1,s_out,t_out,Fprof,Fprof_s,Fprof_t,Fprof_st,Fprof_ss,Fprof_tt)
            BR = ( A30_Z - AZ0_p )/ R
            BZ = ( AR0_p - A30_R )/ R
            Bp = ( AZ0_R - AR0_Z ) + Fprof / R
! reduced-MHD
#else
            psi_R  = 0.d0 ; psi_Z  = 0.d0
            do i_tor = 1, n_tor
              call interp(node_list,element_list,i_elm,var_psi, i_tor,s_out,t_out,psi,psi_s,psi_t,psi_st,psi_ss,psi_tt)
              if ((xjac .gt. 1.d-6)) then  ! avoid the axis
                psi_R  = psi_R  + (   Z_t * psi_s - Z_s * psi_t ) / xjac * HHZ(i_tor)
                psi_Z  = psi_Z  + ( - R_t * psi_s + R_s * psi_t ) / xjac * HHZ(i_tor)
              endif
            enddo
            BR = + psi_Z / R
            BZ = - psi_R / R
            Bp = F0 / R
#endif
            
            R_line = R_line + R*delta_phi / Bp * BR
            Z_line = Z_line + R*delta_phi / Bp * BZ
          else
            cycle
          endif
          
          ! --- record connection length
          xyz(1) = R_line * cos(phi)
          xyz(2) = Z_line
          xyz(3) = R_line * sin(phi)
          total_length = total_length + sqrt((xyz(1)-xyz_prev(1))**2 + (xyz(2)-xyz_prev(2))**2 + (xyz(3)-xyz_prev(3))**2)
          xyz_prev = xyz
          
          ! --- Step in toroidal direction
          phi = phi + delta_phi
          
        end do ! n_phi
      
        ! -----------------
        ! --- record    ---
        ! -----------------
        call find_RZ(node_list,element_list,R_line,Z_line,R_out,Z_out,i_elm,s_out,t_out,ifail)
        if ( ifail /= 0 ) cycle
        R_turn(i_turn+1,(i_dir+1)/2+1) = R_line
        Z_turn(i_turn+1,(i_dir+1)/2+1) = Z_line
        C_turn(i_turn+1,(i_dir+1)/2+1) = total_length
        call var_value(i_elm, var_T,   s_out,t_out,P_start, T_turn  (i_turn+1,(i_dir+1)/2+1) )
        call var_value(i_elm, var_A3,  s_out,t_out,P_start, PSI_turn(i_turn+1,(i_dir+1)/2+1) )
        call var_value(i_elm, var_rho, s_out,t_out,P_start, ZN_turn (i_turn+1,(i_dir+1)/2+1) )
          
      end do ! n_turns
      
    end do ! i_dir
    
    ! --- Record lines that started inside separatrix only?
    if (     ( (FF_0 .gt. 0.d0) .and. (PSI_turn(1,1) .lt. ES%psi_bnd) ) &
        .or. ( (FF_0 .lt. 0.d0) .and. (PSI_turn(1,1) .gt. ES%psi_bnd) ) &
        .or. (.not. ignore_SOL)                                     )then
  
      ! --- Compute averaged connection length
      C_all = C_turn(n_turns+1,1) + C_turn(n_turns+1,2)
      count_lines = count_lines + 1
      C_average = C_average + C_all
      
      ! --- Do both directions
      do i_dir=1,2
  
        ! --- Do all turns
        do i_turn=1,n_turn_max(i_dir)+1
  
          ! --- Safety check
          if (R_turn(i_turn,i_dir) .gt. 0.d0) then
  
            ! --- Connection lengths in both directions
            C_min = C_all - C_turn(i_turn,i_dir)
  
            ikeep = ikeep + 1

            ! --- Polar coordinates
            small_r   = sqrt( (R_turn(i_turn,i_dir)-ES%R_axis)**2 + (Z_turn(i_turn,i_dir)-ES%Z_axis)**2 )
            theta_pol = atan2(Z_turn(i_turn,i_dir)-ES%Z_axis,R_turn(i_turn,i_dir)-ES%R_axis)
            if (theta_pol .lt. 0.d0) theta_pol = theta_pol + 2.d0*PI

            ! --- Save point
            if (ikeep .le. n_points_max) then
              if (Rtheta_plot) then
                RZkeep(1,ikeep)          = small_r
                RZkeep(2,ikeep)          = theta_pol / (2.d0*PI)
              elseif (psitheta_plot) then
                call find_RZ(node_list,element_list,R_turn(i_turn,i_dir),Z_turn(i_turn,i_dir),R_out,Z_out,i_elm,s_out,t_out,ifail)
                call interp(node_list,element_list,i_elm,1,1,s_out,t_out,psi,P0_s,P0_t,P0_st,P0_ss,P0_tt)
                RZkeep(1,ikeep)          = (psi-ES%psi_axis)/(ES%psi_bnd-ES%psi_axis)
                RZkeep(2,ikeep)          = theta_pol / (2.d0*PI)
              else
                RZkeep(1,ikeep)          = R_turn(i_turn,i_dir)
                RZkeep(2,ikeep)          = Z_turn(i_turn,i_dir)
              endif
              ! --- Variables on point
              scalars(ikeep,1:n_scalars) = (/ C_all, C_min, PSI_turn(1,i_dir), T_turn(1,i_dir) /)
            else
              write(*,*) 'Warning! Exceeded maximum number of points',n_points_max
            endif
  
          endif
        enddo
      
      enddo ! end recording in both directions
  
    endif
  end do

  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  
  if (my_id .eq. 0) write(*,*)'Finished loop!!!'
  
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
  
  if (my_id .eq. 0) write(*,*)'Writing connection.vtk'
  
  ! --- An arbitrary unit number for the VTK output file
  ivtk = 22
  
  ! --- The local number of points (for mpi_0)
  ikeep0  = ikeep
  call MPI_Reduce(ikeep,nnos,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  if (my_id .eq. 0) write(*,*) ' number of points : ',nnos
  
  ! --- Line feed character
  lf = char(10)
  
  ! --- Open file and write headers
  if (my_id .eq. 0) then
#ifdef IBM_MACHINE
    open(unit=ivtk,file='connection.vtk',form='unformatted',access='stream',status='replace')
#else
    open(unit=ivtk,file='connection.vtk',form='unformatted',access='stream',convert='BIG_ENDIAN',status='replace')
#endif
    buffer = '# vtk DataFile Version 3.0'//lf                                             ; write(ivtk) trim(buffer)
    buffer = 'vtk output'//lf                                                             ; write(ivtk) trim(buffer)
    buffer = 'BINARY'//lf                                                                 ; write(ivtk) trim(buffer)
    buffer = 'DATASET UNSTRUCTURED_GRID'//lf//lf                                          ; write(ivtk) trim(buffer)
    write(str1(1:12),'(i12)') nnos
    buffer = 'POINTS '//str1//'  float'//lf                                               ; write(ivtk) trim(buffer)
  endif
  
  ! --- Write points for local MPI (id=0)
  if (my_id .eq. 0) then
    write(ivtk) ( (/RZkeep(1,i), RZkeep(2,i), ZERO /),i=1,ikeep0)
  endif
  
  ! --- Write points for all other MPIs
  if (my_id .eq. 0) then
    ! --- If this is mpi_0, we receive data from the other MPIs and print it
    do j=1,n_cpu-1
      call mpi_recv(ikeep,1, MPI_INTEGER, j, j, MPI_COMM_WORLD, status, ierr)
      if (ikeep .gt. 0) then
        nrecv = 2*ikeep
        call mpi_recv(RZkeep,nrecv, MPI_DOUBLE_PRECISION, j, j, MPI_COMM_WORLD, status, ierr)
        write(ivtk) ( (/RZkeep(1,i), RZkeep(2,i), ZERO /),i=1,ikeep)
      endif
    enddo
  else
    ! --- If this is not mpi_0, we send data to the main MPI 0
    call mpi_send(ikeep, 1, MPI_INTEGER, 0, my_id, MPI_COMM_WORLD, ierr)
    if (ikeep .gt. 0) then
      nsend = 2*ikeep
      call mpi_send(RZkeep, nsend,MPI_DOUBLE_PRECISION, 0, my_id, MPI_COMM_WORLD, ierr)
    endif
  endif
  
  ! --- Now for the variables on points (headers first)
  if (my_id .eq. 0) then
    write(str1(1:12),'(i12)') nnos
    buffer = lf//lf//'POINT_DATA '//str1//lf                                              ; write(ivtk) trim(buffer)
  endif
  
  ! --- Loop over all variables
  do i_var =1, n_scalars

    ! --- More headers
    if (my_id .eq. 0) then
      buffer = 'SCALARS '//scalar_names(i_var)//' float'//lf                              ; write(ivtk) trim(buffer)
      buffer = 'LOOKUP_TABLE default'//lf                                                 ; write(ivtk) trim(buffer)
    endif
  
    ! --- Write data for local MPI (id=0)
    if (my_id .eq. 0) then
      write(ivtk) (scalars(i,i_var),i=1,ikeep0)
    endif
  
    ! --- Write data for all other MPIs
    if (my_id .eq. 0) then
      ! --- If this is mpi_0, we receive data from the other MPIs and print it
      do j=1,n_cpu-1
        call mpi_recv(ikeep,1, MPI_INTEGER, j, j, MPI_COMM_WORLD, status, ierr)
        if (ikeep .gt. 0) then
          nrecv = ikeep
          call mpi_recv(scalars(1:ikeep,i_var),nrecv, MPI_DOUBLE_PRECISION, j, j, MPI_COMM_WORLD, status, ierr)
          write(ivtk) (scalars(i,i_var),i=1,ikeep)
        endif
      enddo
    else
      ! --- If this is not mpi_0, we send data to the main MPI 0
      call mpi_send(ikeep, 1, MPI_INTEGER, 0, my_id, MPI_COMM_WORLD, ierr)
      if (ikeep .gt. 0) then
        nsend = ikeep
        call mpi_send(scalars(1:ikeep,i_var), nsend, MPI_DOUBLE_PRECISION, 0, my_id, MPI_COMM_WORLD, ierr)
      endif
    endif
  
  enddo
  
  close(ivtk)
  
  deallocate(RZkeep,scalars,scalar_names)
  
  
  call MPI_FINALIZE(IERR)                                ! clean up MPI

end




































! --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
! --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
! --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
! --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
! --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
! --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------








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
