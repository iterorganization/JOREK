!> Program to put JOREK data into IMAS IDSs
program jorek2_IDS

#ifdef USE_IMAS
  use ids_schemas 
  use ids_routines
  use equil_info, only: print_equil_state

  use mod_jorek2IMAS 
  use constants
  use nodes_elements
  use data_structure
  use mod_import_restart
  use phys_module
  use mod_boundary,        only: boundary_from_grid
  use vacuum
  use vacuum_response,     only: get_vacuum_response, broadcast_starwall_response, init_wall_currents
  use vacuum_equilibrium,  only: import_external_fields
  use mod_impurity, only: init_imp_adas 
  use basis_at_gaussian, only: initialise_basis
  
  implicit none
  
  character(len=1024):: user, database, passive_coil_geo_file, active_coil_geo_file, URI
  character(len=64) :: file_name, name_proj, dd_version_maj, backend, str_shot, str_run
  integer :: shot_number, run_number, i_begin, i_end, i_step, i_jump_steps, i_fmt
  integer :: ierr, idx, stat_mhd, stat_core, stat_rad, stat_eq, n_grid, stat, stat_wall, n_phi_PFC_wall=32
  integer :: stat_pass, stat_act, stat_sum, stat_dis, stat_vac, stat_spi, stat_trans
  integer :: summary_disruption_type, summary_disruption_vertical_displacement
  logical :: first_step, file_exists, rad_only_projections_h5, overwrite_entry
  logical :: export_JOREK_variables, export_radiation, export_1d_profiles, export_equilibrium
  logical :: export_wall, export_pf_passive, export_pf_active, export_summary, export_disruption
  logical :: export_field_extension, new_entry, export_spi, export_transport
  real*8  :: rho0, fact_time, time_SI, wall_thickness
  real*8, allocatable :: res0D(:), res_bnd(:,:,:)
  real*8, allocatable :: avg(:,:)    ! average like expr_avg_list + q_prof + rho_tor
  real*8, allocatable :: q_prof(:), rho_tor(:)


  integer   :: my_id, my_id_n, my_id_master, ierr2
  integer   :: i_rank(n_tor), n_mpi, n_mpi_n, n_mpi_master, m_mpi, n_masters, n_mpi_trans, my_id_trans
  integer   :: MPI_COMM_N, MPI_GROUP_MASTER, MPI_GROUP_WORLD, MPI_COMM_MASTER, MPI_COMM_TRANS
  integer   :: required,provided,StatInfo, resultlength
  integer*4 :: rank, comm_size 
  character(len=MPI_MAX_PROCESSOR_NAME) :: name
  character(len=1000) :: simulation_description

  type(ids_plasma_profiles), target   :: plasma_profiles_ids1, plasma_profiles_vac_extension
  type(ids_plasma_profiles)           :: plasma_profiles_ids
  type(ids_equilibrium)   :: equilibrium_ids
  type(ids_summary)       :: summary_ids
  type(ids_radiation)     :: radiation_ids
  type(ids_wall), target  :: wall_ids
  type(ids_pf_passive)    :: pf_passive
  type(ids_pf_active)     :: pf_active
  type(ids_disruption)    :: disruption_ids
  type(ids_spi)           :: spi_ids
  type(ids_plasma_transport), target :: transport_ids
  type(t_rect_grid_params):: rect_grid_params
  type(t_PFC_triang_grid) :: PFC_wall_grid
  type(t_expr_list) :: expr_avg_list

  namelist /imas_params/ URI, shot_number, run_number, user, database, i_begin, i_end, &
                         export_JOREK_variables, export_radiation, export_1d_profiles, n_grid, &
                         export_equilibrium, rad_only_projections_h5, export_wall,   &
                         export_pf_passive, export_pf_active, passive_coil_geo_file, &
                         active_coil_geo_file, wall_thickness, export_disruption,    &
                         export_summary, overwrite_entry, i_jump_steps,              &
                         simulation_description, export_field_extension,             &
                         rect_grid_params, backend, export_spi, export_transport,    &
                         summary_disruption_type, summary_disruption_vertical_displacement

  ! --- Necessary initialization ------------------
  ! --- MPI initialization (for wall current reconstruction)
#ifdef FUNNELED
    required = MPI_THREAD_FUNNELED
#else
    required = MPI_THREAD_MULTIPLE
#endif

  call MPI_Init_thread(required, provided, StatInfo)

  call init_threads()  ! on some systems init_threads needs to come after mpi_init_thread

  ! --- Determine number of MPI procs
  call MPI_COMM_SIZE(MPI_COMM_WORLD, comm_size, ierr)
  n_mpi = comm_size
  
  if (n_mpi > 1) then
    write(*,*) '  jorek2_IDS needs to be adapated for several MPI processes'
    write(*,*) '  please run with 1 MPI process for the moment'
    write(*,*) n_mpi
    stop
  endif

  ! --- Determine ID of each MPI proc
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  my_id = rank

  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  CALL MPI_GET_PROCESSOR_NAME (name,resultlength,ierr)
  call MPI_Barrier(MPI_COMM_WORLD,ierr)

  ! --- Initialize mode and mode_type arrays
  call det_modes()
  
  ! --- Initialize the basis functions
  call initialise_basis()
  
  ! --- Preset namelist input parameters
  call preset_parameters()

  ! --- Initialize and broadcast input parameters
  call initialise_and_broadcast_parameters(my_id, "__NO_FILENAME__", .false.)

  ! --- Ensure that aux_node_list is associated
  if (.not. associated(aux_node_list)) allocate(aux_node_list)
  
  ! --- Initialize the vacuum part.
  call vacuum_init(my_id, freeboundary_equil, freeboundary, resistive_wall)

  call MPI_Barrier(MPI_COMM_WORLD,ierr)

  ! --- Initialize ADAS
#if (defined WITH_Neutrals) || (defined WITH_Impurities)
    ! --- Read ADAS data and generate coronal equilibrium if needed
    call init_imp_adas(0)
#else
    if (use_imp_adas .and. (nimp_bg(1) > 0.d0)) then
      call init_imp_adas(0)
    endif
#endif
  ! ------------------ end initialization ------------------------
  
  ! --- Preset parameters for this program
  URI         = ''                        !< Optional complete IMAS URI; overrides the fields below
  backend     = 'hdf5'                    !< Name of the backend to store the data (mdsplus,hdf5...)
  database    = 'test'                    !< Name of the database to export the results
  shot_number = 111112;   run_number=1;   
  i_begin     = 0                         !< Starting restart file index
  i_end       = 99999                     !< Ending restart file index
  i_jump_steps= 1                         !< Jump this many steps to read the next restart file
  export_JOREK_variables = .true.
  export_field_extension = .false.        !< Export magnetic field also into the vacuum region? (freeboundary only)
  export_radiation     = .false.
  export_1d_profiles   = .false. 
  export_equilibrium   = .false.
  export_wall          = .false.
  export_pf_passive    = .false.
  export_pf_active     = .false.
  export_summary       = .false.
  export_disruption    = .false.
  export_spi           = .false.
  export_transport     = .false.      
  summary_disruption_type = -999       !< Optional summary/disruption/type/index
  summary_disruption_vertical_displacement = -999  !< Optional summary/disruption/vertical_displacement/value
  rad_only_projections_h5 = .false.    !< use only *.h5 projection files for radiation IDS (single jorek_restart.h5 still needed)
  overwrite_entry      = .false.       !< If true, it overwrites the shot/run even if it already exists in the database.
                                       !  Otherwise it appends the IDSs to the existing entry
  n_grid               = 100           !< Number of points used for 1D and 2D profiles  
  wall_thickness       = 0.06          !< Thickness used for the STARWALL thin wall (default value is for ITER)
  passive_coil_geo_file= 'None'
  active_coil_geo_file = 'None'
  simulation_description = 'JOREK simulation'

  call getenv('USER',user)
  
  ! --- Read parameters from namelist file 'imas.nml' if it exists
  open(42, file='imas.nml', action='read', status='old', iostat=ierr)
  if ( ierr == 0 ) then
    write(*,*) 'Reading parameters from imas.nml namelist.'
    read(42, imas_params)
    close(42)
  end if

  ! --- Use a URI supplied in the namelist in preference to its component fields.
  if (len_trim(URI) == 0) then
    write(str_run,  '(I0)') run_number
    write(str_shot, '(I0)') shot_number
    write(dd_version_maj, '(I0)') al_dd_major_version

    URI = "imas:" // trim(backend) // "?user="     // trim(user)     // ";pulse="    // TRIM(str_shot)      // &
          ";run=" // TRIM(str_run) // ";database=" // trim(database) // ";version=" //  TRIM(dd_version_maj)
  endif

  write(*,*) ' Exporting to URI = '//trim(URI)
  new_entry = .true.
  if (overwrite_entry) then
    call imas_open(URI, FORCE_CREATE_PULSE, idx, stat)
  else
    ! --- Try to open shot number if it exists
    call imas_open(URI, OPEN_PULSE, idx, stat)
    new_entry = .false.

    if (stat /= 0) then  ! --- Create a new shot if it doesn't exist
      write(*,*) '  Shot/run number did not exist, creating new one...'
      call imas_open(URI, FORCE_CREATE_PULSE, idx, stat)
      new_entry = .true.
    endif
  endif

  if (export_radiation .and. use_marker)  then
    allocate( aux_node_list )
  end if

  ! --- Time normalization
  rho0       = central_density * 1.d20 * central_mass * ATOMIC_MASS_UNIT
  fact_time  = sqrt( mu_zero * rho0 )

  first_step = .true.

  call initialise_postproc_settings(first_step, n_grid)
  
  ! --- Loop over
  do i_step = i_begin, i_end, i_jump_steps
 
    ! --- Cycle when required files don't exist 
    if (rad_only_projections_h5 .and. export_radiation) then
      write(name_proj,'(a,i5.5,a)') 'projections', i_step, '.h5'  ! This formatting should be improved
      inquire (file=trim(name_proj), exist=file_exists)
      if (.not. file_exists) then
        write(*,*) 'Warning:',name_proj,'is not found'
        cycle
      endif
      file_name = 'jorek_restart' 
      fact_time = 1.d0  ! Projection times are already in SI units...
    else
      i_fmt = restart_file_exists(i_step)  ! -1 if it does not exist, otherwise index for restart file digit format
      if (i_fmt < 0) cycle
      write(file_name,rst_file_ind_fmt(i_fmt)) 'jorek', i_step
      write(name_proj,rst_file_ind_fmt(i_fmt)) 'projections', i_step
    endif

    ! --- Import restart file
    write(*,*)
    write(*,'(a,i9.9,a)') '#################### STEP ', i_step, ' ####################'
    write(*,*)
    call import_restart(node_list, element_list, file_name, rst_format, ierr)
    call print_equil_state(.true.)
    call boundary_from_grid(node_list, element_list, bnd_node_list, bnd_elm_list, output_bnd_elements)
    if (ierr /=0 ) then
       write(*,*) '  Could not read the JOREK restart file'
       stop
    endif
    time_SI = t_start * fact_time

    ! --- Check whether jorek2_IDS has been compiled with a sufficient number of harmonics
    if (n_tor < n_tor_restart) then
      write(*,*) '  You must compile jorek2_IDS at least with n_tor = ', n_tor_restart
      write(*,*) '  Otherwise you are cutting information for the IDSs'
      stop
    endif

    ! --- Read STARWALL response to export wall currents for wall_IDS
    if (first_step .and. freeboundary .and. &
        (export_wall .or. export_pf_passive .or. export_pf_active .or. export_field_extension)) then
      call get_vacuum_response(my_id, node_list, bnd_elm_list, bnd_node_list, freeboundary_equil,    &
           resistive_wall)
      call import_external_fields('coil_field.dat', my_id)
      if ( .not. wall_curr_initialized ) call init_wall_currents(my_id, resistive_wall)
   endif

    ! --- Compute integrals 
    if (export_equilibrium .or. export_summary .or. export_disruption)  then
      call get_0D_data(first_step, res0D)
    endif

    ! --- Compute flux average profiles
    if (export_equilibrium .or. export_1d_profiles .or. export_radiation)  then
      call get_average_data(first_step, n_grid, expr_avg_list, avg)
    end if

    ! --- Get boundary data and create grid for the JOREK boundary
    if (export_wall .or. export_transport) then
      call get_boundary_data(n_phi=n_phi_PFC_wall, result=res_bnd)
      if (first_step) call triangulate_thin_wall_from_bnd(bnd_points=res_bnd, PFC_wall=PFC_wall_grid)
    end if

    ! --- Fill and export a plasma profiles IDS with the JOREK variables
    if (export_JOREK_variables)  call fill_profiles_w_JOREK_var(first_step, time_SI, plasma_profiles_ids1)  

    ! --- Fill and export a plasma_profiles IDS
    if (export_1d_profiles)  call fill_plasma_profiles_IDS(first_step, time_SI, expr_avg_list, avg, plasma_profiles_ids, n_grid)  


    if (export_equilibrium) call fill_equilibrium_IDS(first_step, time_SI, n_grid, res0D, expr_avg_list, avg, equilibrium_ids, rect_grid_params)  
    if (export_summary)     call fill_summary_IDS(first_step, time_SI, res0D, summary_ids, simulation_description, &
                                                   summary_disruption_type, summary_disruption_vertical_displacement)
    if (export_disruption)  call fill_disruption_IDS(first_step, time_SI, res0D, disruption_ids) 
    if (export_transport)   call fill_transport_IDS(first_step, time_SI, transport_ids, res_bnd, PFC_wall_grid)


    ! --- Extend the fields to vacuum is possible
    if (export_field_extension) then
      if (freeboundary) then
        if (export_equilibrium) then
          call fill_fields_vacuum_extension(first_step, time_SI, plasma_profiles_vac_extension, rect_grid_params, equilibrium_ids) 
        else
          call fill_fields_vacuum_extension(first_step, time_SI, plasma_profiles_vac_extension, rect_grid_params) 
        endif
      else
        write(*,*) 'ERROR: You need the freeboundary extension to calculate magnetic fields in the vacuuum'
        stop
      endif
    endif

    ! --- Fill and export a wall IDS
    if (export_wall)  call fill_wall_IDS(first_step, time_SI, wall_thickness, wall_ids, res_bnd, PFC_wall_grid)  

    ! --- Fill and export a pf_passive IDS
    if (export_pf_passive)  call fill_pf_passive_IDS(first_step, time_SI, pf_passive, passive_coil_geo_file)  

    ! --- Fill and export a pf_active IDS
    if (export_pf_active)   call fill_pf_active_IDS(first_step, time_SI, pf_active, active_coil_geo_file)  

    ! --- Fill and export a radiation IDS
    if (export_radiation) then
      if (use_marker) then
        call import_hdf5_restart_aux(aux_node_list, trim(name_proj) // '.h5', rst_format, ierr)
        if (ierr /= 0) then
          write(*,*) ' Could not open projections file where radiation is stored'
          stop
        endif
      endif
      call fill_radiation_IDS(first_step, time_SI, expr_avg_list, avg, n_grid, radiation_ids)  
    end if
    ! --- Fill and export an SPI IDS
    if (export_spi)   call fill_spi_IDS(first_step, time_SI, spi_ids) 

    stat_mhd = 1;   stat_core = 1;   stat_rad = 1;   stat_eq  = 1;   stat_wall = 1;
    stat_pass= 1;   stat_act  = 1;   stat_sum = 1;   stat_dis = 1;   stat_vac  = 1;
    stat_spi = 1;   stat_trans= 1;

    ! --- Put IDSs into database
    if (first_step .and. new_entry) then  
      if (export_1d_profiles)      call ids_put(idx,'plasma_profiles',plasma_profiles_ids,stat_core)
      if (export_JOREK_variables)  call ids_put(idx,'plasma_profiles/1',plasma_profiles_ids1,stat_mhd)
      if (export_field_extension)  call ids_put(idx,'plasma_profiles/2',plasma_profiles_vac_extension,stat_vac)
      if (export_equilibrium)      call ids_put(idx,'equilibrium',equilibrium_ids,stat_eq)
      if (export_radiation)        call ids_put(idx,'radiation',radiation_ids,stat_rad)
      if (export_wall)             call ids_put(idx,'wall',wall_ids,stat_wall)
      if (export_pf_passive)       call ids_put(idx,'pf_passive',pf_passive,stat_pass)
      if (export_pf_active)        call ids_put(idx,'pf_active',pf_active,stat_act)
      if (export_summary)          call ids_put(idx,'summary',summary_ids,stat_sum)
      if (export_disruption)       call ids_put(idx,'disruption',disruption_ids,stat_dis)
      if (export_spi)              call ids_put(idx,'spi',spi_ids,stat_spi)
      if (export_transport)        call ids_put(idx,'plasma_transport',transport_ids,stat_trans)
    else
      if (export_1d_profiles)      call ids_put_slice(idx,'plasma_profiles',plasma_profiles_ids,stat_core)
      if (export_JOREK_variables)  call ids_put_slice(idx,'plasma_profiles/1',plasma_profiles_ids1,stat_mhd)
      if (export_field_extension)  call ids_put_slice(idx,'plasma_profiles/2',plasma_profiles_vac_extension,stat_vac)
      if (export_equilibrium)      call ids_put_slice(idx,'equilibrium',equilibrium_ids,stat_eq)
      if (export_radiation)        call ids_put_slice(idx,'radiation',radiation_ids,stat_rad)
      if (export_wall)             call ids_put_slice(idx,'wall',wall_ids,stat_wall)
      if (export_pf_passive)       call ids_put_slice(idx,'pf_passive',pf_passive,stat_pass)
      if (export_pf_active)        call ids_put_slice(idx,'pf_active',pf_active,stat_act)
      if (export_summary)          call ids_put_slice(idx,'summary',summary_ids,stat_sum)
      if (export_disruption)       call ids_put_slice(idx,'disruption',disruption_ids,stat_dis)
      if (export_spi)              call ids_put_slice(idx,'spi',spi_ids,stat_spi)
      if (export_transport)        call ids_put_slice(idx,'plasma_transport',transport_ids,stat_trans)
    endif

    if (export_JOREK_variables.and. (stat_mhd==0 ))  write(*,*) '    JOREK variables exported to plasma profiles IDS'
    if (export_field_extension.and.(stat_vac==0 ))   write(*,*) '    Vacuum extension exported'
    if (export_1d_profiles   .and. (stat_core==0))   write(*,*) '    1D profiles exported to plasma profiles IDS'
    if (export_equilibrium   .and. (stat_eq==0  ))   write(*,*) '    Equlibrium IDS exported'
    if (export_radiation     .and. (stat_rad==0 ))   write(*,*) '    Radiation IDS exported'
    if (export_wall          .and. (stat_wall==0 ))  write(*,*) '    Wall IDS exported'
    if (export_pf_passive    .and. (stat_pass==0 ))  write(*,*) '    Pf passive IDS exported'
    if (export_pf_active     .and. (stat_act==0 ))   write(*,*) '    Pf active IDS exported'
    if (export_summary       .and. (stat_sum==0 ))   write(*,*) '    Summary IDS exported'
    if (export_disruption    .and. (stat_dis==0 ))   write(*,*) '    Disruption IDS exported'
    if (export_SPI           .and. (stat_spi==0 ))   write(*,*) '    SPI IDS exported'
    if (export_transport     .and. (stat_trans==0 )) write(*,*) '    Plasma transport IDS exported'

    if (export_JOREK_variables.and. (stat_mhd/=0 ))  write(*,*) '    Problem saving JOREK variables to plasma profiles IDS'
    if (export_field_extension.and.(stat_vac/=0 ))   write(*,*) '    Problem saving vacuum extension'
    if (export_1d_profiles   .and. (stat_core/=0))   write(*,*) '    Problem saving 1D profiles to plasma profiles IDS'
    if (export_equilibrium   .and. (stat_eq/=0  ))   write(*,*) '    Problem saving Equlibrium IDS'
    if (export_radiation     .and. (stat_rad/=0 ))   write(*,*) '    Problem saving Radiation IDS'
    if (export_wall          .and. (stat_wall/=0))   write(*,*) '    Problem saving wall IDS'
    if (export_pf_passive    .and. (stat_pass/=0))   write(*,*) '    Problem saving PF passive IDS'
    if (export_pf_active     .and. (stat_act/=0 ))   write(*,*) '    Problem saving PF active IDS'
    if (export_summary       .and. (stat_sum/=0 ))   write(*,*) '    Problem saving summary IDS'
    if (export_disruption    .and. (stat_dis/=0 ))   write(*,*) '    Problem saving disruption IDS'
    if (export_spi           .and. (stat_spi/=0 ))   write(*,*) '    Problem saving SPI IDS'
    if (export_transport     .and. (stat_trans/=0 )) write(*,*) '    Problem saving Plasma transport IDS'

    first_step = .false.

  enddo

  call imas_close(idx) 
 
#else

  write(*,*) 'Error: jorek2_IDS must be compiled with IMAS (USE_IMAS=1)'
  write(*,*) 'You will also have to load the IMAS module, in case you have not done it'
  write(*,*) '     module load IMAS                                                   '

#endif

end program jorek2_IDS
