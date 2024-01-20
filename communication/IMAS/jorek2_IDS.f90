!> Program to put JOREK data into IMAS IDSs
program jorek2_IDS

#ifdef USE_IMAS
  use ids_schemas 
  use ids_routines, only: imas_open_env, &
     imas_create_env, imas_close, ids_get, ids_put

  use mod_jorek2IMAS 
  use nodes_elements
  use data_structure
  use mod_import_restart
  use phys_module, only : rst_format
  use basis_at_gaussian, only: initialise_basis
  
  implicit none
  
  character(len=200):: user, database
  character(len=64) :: file_name
  integer :: shot_number, run_number, i_begin, i_end, i_step
  integer :: ierr, idx, stat
  logical :: first_step
  logical :: export_MHD, export_radiation

  namelist /imas_params/ shot_number, run_number, user, database, i_begin, i_end, &
                         export_mhd, export_radiation

  ! --- Necessary initialization ------------------
  ! --- Initialize mode and mode_type arrays
  call det_modes()
  
  ! --- Initialize the basis functions
  call initialise_basis()
  
  ! --- Preset namelist input parameters
  call preset_parameters()

  ! --- Read input file
  call initialise_parameters(0, "__NO_FILENAME__")
  ! -----------------------------------------------
  
  ! --- Preset parameters for this program
  database    = 'test'                    ! Name of the database to export the results
  shot_number = 111112;   run_number=1;   
  i_begin     = 0                         ! Starting restart file index
  i_end       = 99999                     ! Ending restart file index
  export_MHD       = .true.
  export_radiation = .false.

  call getenv('USER',user)
  
  ! --- Read parameters from namelist file 'imas.nml' if it exists
  open(42, file='imas.nml', action='read', status='old', iostat=ierr)
  if ( ierr == 0 ) then
    write(*,*) 'Reading parameters from imas.nml namelist.'
    read(42, imas_params)
    close(42)
  end if

  ! --- Try to open shot and number if it exists
  call imas_open_env( 'ids', shot_number,run_number,idx,user,database,'3',stat)! 3 is the database version  
  
  if (stat /= 0) then  ! --- Create a new shot if it doesn't exist
    write(*,*) '  Shot/run number did not exist, creating new one...'
    call imas_create_env('ids',shot_number,run_number, 0,0,idx,user,database,'3') 
  endif

  if (export_radiation)  allocate( aux_node_list )

  first_step = .true.

  ! --- Loop over
  do i_step = i_begin, i_end
  
    ! --- Check whether the restart file exists
    if (.not. restart_file_exists(i_step)) cycle

    ! --- Import restart file
    write(*,*)
    write(*,'(a,i5.5,a)') '#################### STEP ', i_step, ' ####################'
    write(*,*)
    write(file_name,'(a,i5.5)') 'jorek', i_step
    call import_restart(node_list, element_list, file_name, rst_format, ierr)
    if (ierr /=0 ) then
       write(*,*) '  Could not read the restart file'
       stop
    endif

    ! --- Fill and export an MHD IDS
    if (export_mhd)  call fill_mhd_IDS(first_step, idx)  

    ! --- Fill and export a radiation IDS
    if (export_radiation) then
      write(file_name,'(a,i5.5,a)') 'projections', i_step, '.h5'
      call import_hdf5_restart_aux(aux_node_list, file_name, rst_format, ierr)
      if (ierr /= 0) then
        write(*,*) ' Could not open projections file were radiation is stored'
        stop
      endif
      call fill_radiation_IDS(first_step, idx)  
    endif

    first_step = .false.

  enddo

  call imas_close(idx) 
 
#else

  write(*,*) 'Error: jorek2_IDS must be compiled with IMAS (USE_IMAS=1)'
  write(*,*) 'You will also have to load the IMAS module, in case you have not done it'
  write(*,*) '     module load IMAS                                                   '

#endif

end program jorek2_IDS
