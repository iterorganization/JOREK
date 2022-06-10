!> Program for post-processing jorek data.
!!
!! Documentation at http://jorek.eu/wiki
program jorek2_postproc
  
  use nodes_elements, only: node_list, element_list
  use phys_module
  use mod_new_diag
  use parse_commands, only: read_command, type_command
  use exec_commands,  only: exec_command, specific_help
  use settings,       only: set_setting
  use basis_at_gaussian, only: initialise_basis
  
  implicit none
  
  type(type_command) :: command
  integer            :: ierr
  
  ! --- Initialize mode and mode_type arrays
  call det_modes()
  
  ! --- Initialize the basis functions
  call initialise_basis()
  
  ! --- Initialize the new_diag package.
  call init_new_diag(.false.)
  
  ! --- Preset namelist input parameters
  call preset_parameters()
  
  ! --- Preset some parameters
  call set_setting('units',           '0',     ierr, 'Calculate quantities in which units (1=JOREK, 0=SI)')
  call set_setting('loop_units',      '0',     ierr, 'Use which units for time-loops (1=JOREK, 0=SI)'     )
  call set_setting('linepoints',      '200',   ierr, 'Number of points along a line e.g. for pol_line'    )
  call set_setting('tor_points',      '200',   ierr, 'Number of toroidal points e.g. for tor_line'        )
  call set_setting('surfaces',        '100',   ierr, 'number for flux surfaces e.g. for qprofile'         )
  call set_setting('verbose',         'false', ierr, 'print detailed information?'                        )
  call set_setting('debug',           'false', ierr, 'debugging output for postproc?'                     )
  call set_setting('nsmallsteps',     '3',     ierr, 'numerical parameter for field line tracing'         )
  call set_setting('nmaxsteps',       '2500',  ierr, 'numerical parameter for field line tracing'         )
  call set_setting('deltaphi',        '0.3',   ierr, 'numerical parameter for field line tracing'         )
  call set_setting('rad_range_min',   '0.001', ierr, 'numerical parameter for field line tracing'         )
  call set_setting('rad_range_max',   '0.999', ierr, 'numerical parameter for field line tracing'         )
  call set_setting('nTht',            '32',    ierr, 'numerical parameter for field line tracing'         )
  call set_setting('nsub_bnd',        '4',     ierr, 'points per boundary element for boundary_quantities')
  call set_setting('nsub_vtk',        '4',     ierr, 'number of element subdivisions for vtk plots       ')
  call set_setting('vtk_phi_value',   '0',     ierr, 'phi value at which the 2D vtk is computed          ')
  call set_setting('only_itor',      '-1',     ierr, 'select a single toroidal harmonic for calculations ')
  
  ! --- Print getting started information
  call specific_help('getting_started')
  
  do ! (main loop: Read, parse, and execute one command after the other)
    
    ! --- Read and parse a command line
    call read_command(command, ierr)
    if ( ierr /= 0 ) then
      write(*,*) 'ERROR in read_command'
      cycle
    else if ( command%n_args < 0 ) then
      cycle
    end if
    
    ! --- Exit upon user request
    if ( (command%args(0) == 'exit') .or. (command%args(0) == 'quit') ) stop
    
    ! --- Execute the command
    call exec_command(command, .true., ierr)
    if ( ierr /= 0 ) then
      write(*,*) 'ERROR in exec_command'
      cycle
    end if
    
  end do ! (main loop)
  
end program jorek2_postproc
