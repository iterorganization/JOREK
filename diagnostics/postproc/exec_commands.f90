!> Module for execution of user commands (used by jorek2_postproc)
module exec_commands
  
  use constants
  use mod_parameters
  use phys_module
  use data_structure
  use equil_info
  use nodes_elements
  use mod_boundary
  use basis_at_gaussian
  use mod_new_diag
  use domains
  use parse_commands
  use settings
  use convert_character
  use postproc_help
  use mod_log_params
  use mod_import_restart
  use mod_interp
  use mod_poloidal_currents 
  
  
  
  
  implicit none
  
  
  
  
  character(len=11), parameter, private :: DIR = './postproc/' !< Output goes into this directory!
  
  integer, parameter :: NORMAL_MODE = 1 !< Normal mode
  integer, parameter :: LOOP_S_MODE = 2 !< Mode started by 'for step' and ended by 'done' commands
  integer, parameter :: LOOP_T_MODE = 3 !< Mode started by 'for time' and ended by 'done' commands
  integer :: exec_mode = NORMAL_MODE    !< Current operation mode (NORMAL_MODE or LOOP_X_MODE)
  integer :: loop_min_step              !< Smallest timestep index of for loop
  integer :: loop_max_step              !< Largest timestep index of for loop
  integer :: loop_inc_step              !< Timestep index step width of for loop
  real*8  :: loop_min_time              !< Smallest time of for loop
  real*8  :: loop_max_time              !< Largest time of for loop
  real*8  :: loop_inc_time              !< Timestep time interval width of for loop  
  
  integer, parameter :: MAX_QUEUE_LENGTH  = 9999         !< Maximum length of command queue
  integer            :: n_queued_commands = 0            !< Number of commands in the queue
  type(type_command) :: command_queue(MAX_QUEUE_LENGTH)  !< Queued commands
  
  logical,             private, save :: input_loaded  = .false. !< Has an input file been loaded?
  logical,             private, save :: step_imported = .false. !< Has a restart file been imported?
  logical,             private, save :: dir_created   = .false. !< Postproc directory created?
  logical,             private, save :: verbose
  logical,             private, save :: debug
  type(t_equil_state), private, save :: eq !< Equilibrium state; updated when time step is loaded
  type(t_expr_list),   private, save :: expr_list
  real*8, allocatable, private, save :: result(:,:,:,:), res2d(:,:,:), res1d(:,:), res0d(:), sum(:)
  complex*16, allocatable, private, save :: cp(:,:,:,:)
  real*8,              private, save :: time_now !< Time of current restart file in selected units
  
  ! --- used by average_h5 command:
  real*8, allocatable, private, save :: values(:,:,:,:)
  real*8,                       save :: weight, total_weight
  
  
  
  private
  public exec_command, general_help, specific_help, clean_up
  
  
  
  
  
  save
  
  
  
  
  
  contains
  
  
  
  
  
  !> Execute a command
  subroutine exec_command(command, first_step, ierr)
    
    ! --- Routine parameters
    type(type_command), intent(in)  :: command     !< Command to be executed
    logical,            intent(in)  :: first_step  !< First time step of a for loop?
    integer,            intent(out) :: ierr        !< Error flag
    
    ierr = 0
    
    if ( .not. dir_created ) then
      call system('mkdir -p '//DIR)
      dir_created = .true.
    end if
    
    if ( get_setting('debug',ierr) == 'true' ) then
      write(*,'(a)') 'Exec_command was called with:'
      call print_command(command)
      write(*,*)
    end if
    
    ! --- Free shared module arrays
    call clean_up()
    
    ! --- In normal mode, some commands are directly executed
    if ( exec_mode == NORMAL_MODE ) then
      
      if ( ( .not. input_loaded ) .and. ( trim(command%args(0)) == 'for' ) ) then
        write(*,*) 'ERROR: No namelist input file loaded.'
        call specific_help('namelist')
        ierr = 1
        return
      end if
      
      verbose = get_log_setting('verbose', ierr)
      debug   = get_log_setting('debug', ierr)
      
      if ( verbose .or. debug ) then
        write(*,*)
        write(*,*) '- Executing "'//trim(command%args(0))//'"'
      end if
      
      select case ( trim(command%args(0)) )
        case ( 'average' )
          call average(command, first_step, ierr)
        case ( 'average_h5' )
          call average_h5(command, first_step, ierr)
        case ( 'int2d' )
          call int2d(command, first_step, ierr)
        case ( 'int3d' )
          call int3d(command, first_step, ierr)
        case ( 'equil_params' )
          call equil_params(command, first_step, ierr)
        case ( 'energy_spectrum' )
          call energy_spectrum(command, first_step, ierr)
        case ( 'expressions' )
          call expressions(command, ierr)
        case ( 'expressions_int' )
          call expressions_int(command, ierr)
        case ( 'fluxsurfaces' )
          call fluxsurfaces(command, ierr)
        case ( 'for' )
          call loop_start(command, ierr)
        case ( 'four2d' )
          call four2d(command, ierr)
        case ( 'gourdon' )
          call gourdon(command, first_step, ierr)
        case ( 'grid' )
          call grid(command, ierr)
        case ( 'help' )
          call help(command, ierr)
        case ( 'I_halo_TPF' )
          call I_halo_TPF(command, first_step, ierr)
        case ( 'jnorm_bnd_curr' )
          call jnorm_bnd_RZ(command, ierr)
        case ( 'jorek-units' )
          call select_jorek_units(command, ierr)
        case ( 'pol_line' )
          call pol_line(command, first_step, ierr)
        case ( 'int_along_pol_line' )
          call int_along_pol_line(command, first_step, ierr)
        case ( 'tor_line' )
          call tor_line(command, first_step, ierr)
        case ( 'rectangle' )
          call rectangle(command, first_step, ierr)
        case ( 'mark_coords' )
          call mark_coords(command, ierr)
        case ( 'midplane' )
          call midplane(command, first_step, ierr)
        case ( 'namelist' )
          call load_namelist(command, ierr)
        case ( 'params' )
          call log_parameters(0, .false.)
        case ( 'point' )
          call point(command, first_step, ierr)
        case ( 'qprofile' )
          call qprofile(command, first_step, ierr)
         case ( 'q_at_psin' )
          call q_at_given_psin(command, first_step, ierr)
        case ( 'separatrix' )
          call separatrix(command, ierr)
        case ( 'set' )
          call set(command, ierr)
        case ( 'si-units' )
          call select_si_units(command, ierr)
        case ( 'spi-state' )
          call spi_state(command, first_step, ierr)
        case ( 'timesteps' )
          call timesteps 
        case default
          write(*,*) 'Command "', trim(command%args(0)), '" does not exist'
          call general_help() 
      end select
      
    ! --- In loop mode, commands are queued and afterwards executed for each timestep separately
    else if (( exec_mode == LOOP_S_MODE ) .or. ( exec_mode == LOOP_T_MODE )) then
      
      select case ( trim(command%args(0)) )
        case ( 'expressions', 'expressions_int', 'mark_coords', 'int2d', 'int3d','midplane', 'average', 'point',      &
          'pol_line', 'int_along_pol_line', 'tor_line', 'equil_params', 'qprofile',        &
          'q_at_psin', 'fluxsurfaces', 'separatrix', 'set', 'four2d', 'gourdon', 'jorek-units',         & 
          'jnorm_bnd_curr', 'si-units', 'grid', 'rectangle', 'energy_spectrum', 'average_h5', &
          'I_halo_TPF', 'spi-state')
          call add_to_command_queue(command, ierr)
        case ( 'help' )
          call help(command, ierr)
        case ( 'done' )
          call loop_end(command, ierr)
        case default
          write(*,*) 'Command "', trim(command%args(0)), '" unknown or invalid inside a loop'
          call general_help() 
      end select
        
    end if
    
  end subroutine exec_command
  
  
  
  
  
  !> Finalize a command once a loop has finished (not needed for most commands)
  subroutine finalize_command(command, first_step, ierr)
    
    ! --- Routine parameters
    type(type_command), intent(in)  :: command     !< Command to be executed
    logical,            intent(in)  :: first_step  !< First time step of a for loop?
    integer,            intent(out) :: ierr        !< Error flag
    
    ierr = 0
    
    verbose = get_log_setting('verbose', ierr)
    debug   = get_log_setting('debug', ierr)
    
    if ( debug ) then
      write(*,'(a)') 'Finalize_command was called with:'
      call print_command(command)
      write(*,*)
    end if
    
    select case ( trim(command%args(0)) )
      case ( 'average_h5' )
        call average_h5_finalize(command, first_step, ierr)
      case default
        if (debug) write(*,*) 'Command does not need finalize.'
    end select
    
  end subroutine finalize_command
  
  
  
  
  
  !> Loads a time step from a restart file if the restart file exists
  subroutine load_step(istep, ierr)
    
    ! --- Routine parameters
    integer,            intent(in)  :: istep !< Load this time step
    integer,            intent(out) :: ierr !< Error flag
    
    character(len=64) :: file_name
    logical           :: file_exists
    integer           :: ifail
    
    ierr = 0
    
    write(file_name,'(a,i5.5)') 'jorek', istep
    if ( rst_hdf5 .ne. 0 ) then
      inquire (file=trim(file_name)//'.h5', exist=file_exists)
    else
      inquire (file=trim(file_name)//'.rst', exist=file_exists)
    end if
    if ( .not. file_exists ) then
      ierr = 42
      return
    end if
    
    write(*,*)
    write(*,'(a,i5.5,a)') '#################### TIME STEP ', istep, ' ####################'
    write(*,*)
    
    ! --- Load the restart file
    call import_restart(node_list, element_list, file_name, rst_format, ierr, .true.)
    if ( ierr /= 0 ) return
    call boundary_from_grid(node_list, element_list, bnd_node_list, bnd_elm_list, .false.)
    
    ! --- Locate magnetic axis and X-point.
    call update_equil_state(node_list, element_list, bnd_elm_list, xpoint, xcase, eq)
    
    t_now         = t_start
    index_now     = index_start
    step_imported = .true.
    
    ! (not elegant, admittedly... but guarantees consistent time normalization:)
    call eval_expr(eq, get_int_setting('units', ierr), exprs('t',1),                               &
      pol_pos(node_list,element_list,eq,R=eq%R_axis,Z=eq%Z_axis), tor_pos(phi=0.d0), result, ierr)
    time_now = result(1,1,1,1)
    
  end subroutine load_step
  
  
  
  
  
  !> Is called at the start of a for loop.
  !!
  !! It switches the module behaviour (exec_mode) to LOOP_MODE, i.e., all commands inside the for
  !! loop are not executed immediately but collected in a command queue. All commands of the
  !! command queue are executed for each time step after the 'done' command of the for loop
  !! has been entered.
  subroutine loop_start(command, ierr)
    
    ! --- Routine parameters
    type(type_command), intent(in)     :: command !< Command to be executed
    integer,            intent(out)    :: ierr    !< Error flag
    
    ierr = 0
    
    ! --- Some checks.
    call check_args(command%n_args,ierr,3,5,7);  if ( ierr /= 0 ) return

    ! --- loop through steps
    if ( trim(command%args(1)) == 'step' ) then

      loop_min_step = to_int(command%args(2),ierr)
      if ( ierr /= 0 ) then
        write(*,*) 'ERROR: Wrong syntax for "for" statement.'
        call specific_help('for')
        return
      end if

      if ( command%n_args == 3 ) then ! for step <XXX> do
        if ( trim(command%args(3)) /= 'do' ) then
          write(*,*) 'ERROR: Wrong syntax for "for" statement.'
          call specific_help('for')
          return
        end if
        loop_max_step = loop_min_step
        loop_inc_step = 1
      else if ( command%n_args == 5 ) then! for step <XXX> to <YYY> do
        if ( trim(command%args(3)) /= 'to' ) then
          write(*,*) 'ERROR: Wrong syntax for "for" statement.'
          call specific_help('for')
          return
        end if
        loop_max_step = to_int(command%args(4),ierr)
        if ( ierr /= 0 ) then
          write(*,*) 'ERROR: Wrong syntax for "for" statement.'
          call specific_help('for')
          return
        end if
        if ( trim(command%args(5)) /= 'do' ) then
          write(*,*) 'ERROR: Wrong syntax for "for" statement.'
          call specific_help('for')
          return
        end if
        loop_inc_step = 1
      else if ( command%n_args == 7 ) then! for step <XXX> to <YYY> by <ZZZ> do
        if ( trim(command%args(3)) /= 'to' ) then
          write(*,*) 'ERROR: Wrong syntax for "for" statement.'
          call specific_help('for')
          return
        end if
        loop_max_step = to_int(command%args(4),ierr)
        if ( ierr /= 0 ) then
          write(*,*) 'ERROR: Wrong syntax for "for" statement.'
          call specific_help('for')
          return
        end if
        if ( trim(command%args(5)) /= 'by' ) then
          write(*,*) 'ERROR: Wrong syntax for "for" statement.'
          call specific_help('for')
          return
        end if
        loop_inc_step = to_int(command%args(6),ierr)
        if ( ierr /= 0 ) then
          write(*,*) 'ERROR: Wrong syntax for "for" statement.'
          call specific_help('for')
          return
        end if
        if ( trim(command%args(7)) /= 'do' ) then
          write(*,*) 'ERROR: Wrong syntax for "for" statement.'
          call specific_help('for')
          return
        end if
      end if
    
      exec_mode = LOOP_S_MODE

    else if ( trim(command%args(1)) == 'time' ) then

      loop_min_time = to_float(command%args(2),ierr)
      if ( ierr /= 0 ) then
        write(*,*) 'ERROR: Wrong syntax for "for" statement.'
        call specific_help('for')
        return
      end if

      if ( command%n_args == 3 ) then ! for time <XXX> do
        if ( trim(command%args(3)) /= 'do' ) then
          write(*,*) 'ERROR: Wrong syntax for "for" statement.'
          call specific_help('for')
          return
        end if
        loop_max_time = loop_min_time
        loop_inc_time = 1.d0
      else if ( command%n_args == 5 ) then
        write(*,*) 'ERROR: Wrong syntax for "for" statement.'
        write(*,*) '       "time" mode needs one or three arguments.'
        call specific_help('for')
        return
      else if ( command%n_args == 7 ) then! for step <XXX> to <YYY> by <ZZZ> do
        if ( trim(command%args(3)) /= 'to' ) then
          write(*,*) 'ERROR: Wrong syntax for "for" statement.'
          call specific_help('for')
          return
        end if
        loop_max_time = to_float(command%args(4),ierr)
        if ( ierr /= 0 ) then
          write(*,*) 'ERROR: Wrong syntax for "for" statement.'
          call specific_help('for')
          return
        end if
        if ( trim(command%args(5)) /= 'by' ) then
          write(*,*) 'ERROR: Wrong syntax for "for" statement.'
          call specific_help('for')
          return
        end if
        loop_inc_time = to_float(command%args(6),ierr)
        if ( ierr /= 0 ) then
          write(*,*) 'ERROR: Wrong syntax for "for" statement.'
          call specific_help('for')
          return
        end if
        if ( trim(command%args(7)) /= 'do' ) then
          write(*,*) 'ERROR: Wrong syntax for "for" statement.'
          call specific_help('for')
          return
        end if
      end if

      ! --- Check if values are okay
      if ( loop_max_time .lt. loop_min_time ) then
        write(*,*) 'ERROR: Minimum value has to be larger than maximum value.'
        call specific_help('for')
        return
      end if
      if ( loop_inc_time .lt. 0.d0 ) then
        write(*,*) 'ERROR: Increment needs to be a positive value.'
        call specific_help('for')
        return
      end if

      exec_mode = LOOP_T_MODE

    else
      call report_error('for', 'Wrong syntax.', command)
      return
    end if
    
  end subroutine loop_start
  
  
  
  
  
  !> Is called upon the 'done' command of a for loop and executes all commands enclosed in the loop.
  subroutine loop_end(command, ierr)
    
    include 'omp_lib.h'
    
    ! --- Routine parameters
    type(type_command), intent(in)     :: command !< Command to be executed
    integer,            intent(out)    :: ierr    !< Error flag
    
    ! --- Local variables
    integer              :: jcmd, istep, load_error, n_avail, itime, n_time, iavail, f_avail
    logical              :: first_step, file_exists   ! Is true for the first timestep loaded in the for-loop
    real*8               :: time_loop
    character(len=256)   :: filename
    integer, allocatable :: selected_steps(:)
    real*8 , allocatable :: selected_times(:)
    integer, allocatable :: available_steps(:)
    
    ierr = 0
    
    if ( n_queued_commands == 0 ) then
      write(*,*) 'WARNING: No commands in for-loop.'
      exec_mode = NORMAL_MODE
      return
    end if
    
    first_step = .true.
    if ( exec_mode == LOOP_S_MODE ) then
      exec_mode = NORMAL_MODE
      do istep = loop_min_step, loop_max_step, loop_inc_step
        write(*,*) istep
        call load_step(istep, load_error)
        if ( load_error /= 0 ) cycle
      
        do jcmd = 1, n_queued_commands
        
          call exec_command(command_queue(jcmd), first_step, ierr)  
          if ( ierr /= 0 ) then
            write(*,*) 'ERROR executing the following command (ignoring it):'
            call print_command(command_queue(jcmd))
            call specific_help(command_queue(jcmd)%args(0))
            ierr = 0
          end if
        
        end do
      
        first_step = .false.
      end do
    else
      exec_mode = NORMAL_MODE
      write(*,*) '--------------------------------------------------'
      write(*,*) '   Selects restart files for choosen time points'
      write(*,*) '--------------------------------------------------'

      ! --- Gets xtime from jorek_restart
      call import_restart(node_list, element_list, "jorek_restart", rst_format, ierr, .true.)
      allocate(available_steps(index_start+1))

      ! --- Gets list of available restart files
      n_avail=0
      do istep = 0, index_start
        if ( rst_hdf5 .ne. 0 ) then
          write (filename,'(a, i5.5, a)') 'jorek', istep, '.h5'
          inquire (file=filename, exist=file_exists)
        else
          write (filename,'(a, i5.5, a)') 'jorek', istep, '.rst'
          inquire (file=filename, exist=file_exists)
        end if

        if ( file_exists ) then
          n_avail=n_avail+1					
          available_steps(n_avail) = istep
        end if
      end do

      ! --- Find for each selected time a corresponding step number
      if (loop_max_time .gt. xtime(index_start)) then
        loop_max_time=xtime(index_start)
      end if
      if (loop_min_time .gt. xtime(index_start)) then
        loop_min_time=xtime(index_start)
      end if
      n_time = int(((loop_max_time-loop_min_time)/loop_inc_time))+1
      allocate(selected_steps(n_time))
      f_avail = 1 !< Thats the index of the first available file that will be compared
      do itime = 1,n_time
        time_loop = loop_min_time+loop_inc_time*(itime-1)        
        ! --- Go through list of available (and relevant) restart files
        !     and compare their xtime with time_loop
        do iavail = f_avail, n_avail
          if ( iavail .gt. f_avail ) then
            if ( abs(xtime(available_steps(iavail)) - time_loop) .gt.  &
                 abs(xtime( selected_steps(itime)) - time_loop)      ) then
              exit
            end if
          end if
          selected_steps(itime)=available_steps(iavail)
          !f_avail=iavail+1
        end do
        write(*,'(a,f13.6)') 'For time point: ', time_loop
        write(*,'(a, i5.5, a,f13.6)') 'Selected Step: ', selected_steps(itime),&
        ' at t=',xtime(selected_steps(itime))
        write(*,*) '****************************'
      end do

      loop_min_step=selected_steps(1)
      loop_max_step=selected_steps(size(selected_steps))

      write(*,*) '--------------------------------------------------'

      do istep = 1, size(selected_steps)
        call load_step( selected_steps(istep), load_error)
        if ( load_error /= 0 ) cycle
      
        do jcmd = 1, n_queued_commands
        
          call exec_command(command_queue(jcmd), first_step, ierr)  
          if ( ierr /= 0 ) then
            write(*,*) 'ERROR executing the following command (ignoring it):'
            call print_command(command_queue(jcmd))
            call specific_help(command_queue(jcmd)%args(0))
            ierr = 0
          end if
        
        end do
      
        first_step = .false.
      end do
    end if 


    do jcmd = 1, n_queued_commands
      call finalize_command(command_queue(jcmd), first_step, ierr)
    end do
    
    if ( first_step ) then
      write(*,'(a,i5.5,a,i5.5,a)') 'WARNING: There were no restart files for steps ',              &
        loop_min_step, ' to ', loop_max_step, '.'
    end if
    
    n_queued_commands = 0
    
  end subroutine loop_end
  
  
  
  
  
  !> Add a command to the command queue (for loop)
  subroutine add_to_command_queue(command, ierr)
    
    ! --- Routine parameters
    type(type_command), intent(in)     :: command !< Command to be executed
    integer,            intent(out)    :: ierr    !< Error flag
    
    ierr = 0
    
    if ( n_queued_commands + 1 > MAX_QUEUE_LENGTH) then
      write(*,*) 'ERROR: Too many commands in the command queue of the for loop.'
      ierr = 1
      return
    end if
    
    n_queued_commands = n_queued_commands + 1
    command_queue(n_queued_commands) = command
    
  end subroutine add_to_command_queue
  
  
  
  
  
  !> Retrieve a command from the command queue (currently not used!)
  subroutine get_from_command_queue(command, ierr)
    
    ! --- Routine parameters
    type(type_command), intent(out)    :: command !< Command to be executed
    integer,            intent(out)    :: ierr    !< Error flag
    
    ierr = 0
    
    if ( n_queued_commands < 1 ) then
      write(*,*) 'ERROR: Cannot get a command from an empty queue.'
      ierr = 1
      return
    end if
    
    n_queued_commands = n_queued_commands - 1
    command = command_queue(1)
    command_queue(1:n_queued_commands) = command_queue(2:n_queued_commands+1)
    
  end subroutine get_from_command_queue
  
  
  
  
  
  !> Clean up global module arrays.
  subroutine clean_up()
    if ( allocated(result) ) deallocate(result)
    if ( allocated(res2d)  ) deallocate(res2d)
    if ( allocated(res1d)  ) deallocate(res1d)
    if ( allocated(res0d)  ) deallocate(res0d)
    if ( allocated(cp)     ) deallocate(cp)
  end subroutine clean_up
  
  
  
  
  
  !> Implements the 'set variable value' command
  subroutine set(command, ierr)
  
    ! --- Routine parameters
    type(type_command), intent(in)     :: command !< Command to be executed
    integer,            intent(out)    :: ierr    !< Error flag
      
    ! --- Local variables
    character(len=128)   :: name
    character(len=1024)  :: value
    
    ierr = 0
    
    ! --- Transformation of input data
    !### check param count
    
    if ( command%n_args == 0 ) then
      call print_settings()
    else
      name  = command%args(1)
      value = command%args(2)
      call set_setting(name, value, ierr)
    end if
    
  end subroutine set
  
  
  
  
  
  !> List all existing restart files
  subroutine timesteps()
    
    character(len=256)  :: filename
    logical             :: file_exists
    integer             :: i
    
    write(*,'(a)') 'Available restart files:'
    do i = 0, 99999
      write (filename,'(a, i5.5, a)') 'jorek', i, '.rst'
      inquire (file=filename, exist=file_exists)
      if (file_exists) write(*,'(i6)',advance='no') i
      write (filename,'(a, i5.5, a)') 'jorek', i, '.h5'
      inquire (file=filename, exist=file_exists)
      if (file_exists) write(*,'(i6)',advance='no') i
    end do
    write(*,*)
    
  end subroutine timesteps
  
  
  
  
  
  !> Load a specific namelist input file
  subroutine load_namelist(command, ierr)
    
    use phys_module     
    
    ! --- Routine parameters
    type(type_command), intent(in)  :: command     !< Command to be executed
    integer,            intent(out) :: ierr        !< Error flag
    
    ! --- Local variables
    character(len=1024) ::  filename
    logical             ::  file_exists
    
    ierr = 0
    
    ! --- Some checks.
    call check_args(command%n_args,ierr,1);  if ( ierr /= 0 ) return
    
    filename = trim(command%args(1))
    inquire (file=filename, exist=file_exists)
      
    ! --- Read the input namelist file
    if (file_exists) then
      call initialise_parameters(0, filename)
      input_loaded = .true.
    else
      ierr = 1
      write(*,*) 'ERROR: input file "', trim(filename), '" does not exist.'
      call specific_help('namelist')
    end if

  end subroutine load_namelist
  
  
  
  
  
  !> Check if a restart file has already been imported
  subroutine check_step_imported(ierr)
    integer, intent(out) :: ierr !< Error flag
    
    ierr = 0
    if ( .not. step_imported ) then
      ierr = 1
      write(*,*) 'ERROR: No restart file has been imported yet. Use the "for" loop:'
      call specific_help('for')
    end if
    
  end subroutine check_step_imported
  
  
  
  
  
  !> Check if expressions have been selected.
  subroutine check_exprs_selected(ierr)
    integer, intent(out) :: ierr !< Error flag
    
    ierr = 0
    if ( expr_list%n_expr <= 0 ) then
      ierr = 1
      write(*,*) 'ERROR: No physical expressions selected yet. Use command "expressions":'
      call specific_help('expressions')
    end if
    
  end subroutine check_exprs_selected
  
  
  
  
  
  !> Check number of arguments correct.
  subroutine check_args(nargs,ierr,ok1,ok2,ok3,ok4)
    integer, intent(in)              :: nargs               !< Number of arguments.
    integer, intent(inout)           :: ierr                !< Error code.
    integer, intent(in),    optional :: ok1, ok2, ok3, ok4
    
    ierr = 1
    if ( present(ok1) .and. (nargs==ok1) ) ierr = 0
    if ( present(ok2) .and. (nargs==ok2) ) ierr = 0
    if ( present(ok3) .and. (nargs==ok3) ) ierr = 0
    if ( present(ok4) .and. (nargs==ok4) ) ierr = 0
    if ( ierr /= 0 ) write(*,*) 'ERROR: Wrong number of parameters for command.'
  end subroutine check_args
  
  
  
  
  
  !> Print usage information
  subroutine help(command, ierr)
  
    ! --- Routine parameters
    type(type_command), intent(in)  :: command !< Command to be executed
    integer,            intent(out) :: ierr    !< Error flag
    
    ! --- Local variables
    integer :: i
    
    ierr = 0
    
    if ( command%n_args == 0 ) then
      call general_help()
    else
      do i = 1, command%n_args
        call specific_help(command%args(i))
      end do
    end if
    
  end subroutine help
  
  
  
  
  
  !> Report an error message.
  subroutine report_error(routine, message, command)
    
    ! --- Routine parameters.
    character(len=*),             intent(in) :: routine
    character(len=*),             intent(in) :: message
    type(type_command), optional, intent(in) :: command
    
    write(*,*)
    write(*,*) 'ERROR in '//trim(routine)//': '//trim(message)
    write(*,*)
    
    if ( present(command) ) call specific_help(trim(command%args(0)))
    
  end subroutine report_error
  
  
  
  
  
  !> List or Select Available Expressions.
  subroutine expressions(command, ierr)
    
    ! --- Routine parameters
    type(type_command), intent(in)  :: command     !< Command to be executed
    integer,            intent(out) :: ierr        !< Error flag
    
    ierr = 0
    
    if ( command%n_args == 0 ) then
      
      call print_exprs(exprs_all)
      
    else
      
      expr_list = exprs(command%args(1:command%n_args), command%n_args)
      call print_exprs(expr_list,.true.)
      
    end if
    
  end subroutine expressions
  
  


  !> List or Select Available Expressions.
  subroutine expressions_int(command, ierr)
    
    ! --- Routine parameters
    type(type_command), intent(in)  :: command     !< Command to be executed
    integer,            intent(out) :: ierr        !< Error flag
    
    ierr = 0
    
    if ( command%n_args == 0 ) then
      
      call print_exprs(exprs_all_int)
      
    else
      
    expr_list = exprs_int(command%args(1:command%n_args), command%n_args)
    call print_exprs(expr_list,.true.)
       
    end if
    
  end subroutine expressions_int
 
  
  
  !> Mark some expressions as coordinates.
  subroutine mark_coords(command, ierr)
    
    ! --- Routine parameters
    type(type_command), intent(in)  :: command     !< Command to be executed
    integer,            intent(out) :: ierr        !< Error flag
    
    ! --- Local variables
    integer :: n_coord
    
    ierr = 0
    
    call check_args(command%n_args,ierr,1);  if ( ierr /= 0 ) return
    n_coord = to_int(command%args(1), ierr); if ( ierr /= 0 ) return
    expr_list%n_coord = n_coord
    
  end subroutine mark_coords
  
  
  
  
  
  !> Auxilliary routine for file name construction.
  character(len=64) function step_range_string(min_step, max_step)

    integer, intent(in) :: min_step, max_step
    character(len=2)     :: prefix

    if ( exec_mode .eq. LOOP_S_MODE )  then 
      prefix='_s'
    else
      prefix='_t'
    end if

    if ( min_step /= max_step ) then
      write(step_range_string,'(a,i5.5,a,i5.5)') prefix, min_step, '..', max_step
    else
      write(step_range_string,'(a,i5.5)') prefix, min_step
    end if

  end function step_range_string
  
  
  
  
  
  !> Read several .h5 files and create an "average" one (average of absolute values).
  subroutine average_h5(command, first_step, ierr)
    
    use mod_export_restart
    
    ! --- Routine parameters
    type(type_command), intent(in)  :: command     !< Command to be executed
    logical,            intent(in)  :: first_step  !< First time step of a for loop?
    integer,            intent(out) :: ierr        !< Error flag
    
    ! --- Local variables
    integer :: i
    integer, save :: prev_index
    
    ierr = 0
    
    ! --- Some checks
    call check_args(command%n_args,ierr,0);  if ( ierr /= 0 ) return
    call check_step_imported(ierr);          if ( ierr /= 0 ) return
    
    weight = 0.d0
    if (first_step) then
      allocate(values(n_tor,n_order+1,n_var,node_list%n_nodes))
      values = 0.d0
    else
      weight = xtime(index_now)-xtime(prev_index)
    end if
    total_weight = total_weight + weight
    prev_index = index_now
    
    do i = 1, node_list%n_nodes
      values(:,1,:,i) = values(:,1,:,i) + abs(node_list%node(i)%values(:,1,:)) * weight
      ! Since the purpose of this postproc command is to visualize the localization of a
      ! particular mode activity, we take the time average over the absolute value. In the
      ! Bezier representation, this we have to throw away the degrees of freedomn 2,3,4 in
      ! doing so, since their respective basis functions are not only positive but change
      ! sign. Consequently, the absolute value of them cannot be represented in the same
      ! basis.
      ! As a result, the created h5 file contains the time average of the absolute values
      ! with a limited resolution since only the first dof is kept on each node.
    end do
    
    ! see also average_h5_finalize below, which writes out the result
    
  end subroutine average_h5
  
  
  
  
  
  !> Read several .h5 files and create an "average" one (average of absolute values).
  subroutine average_h5_finalize(command, first_step, ierr)
    
    use mod_export_restart
    
    ! --- Routine parameters
    type(type_command), intent(in)  :: command     !< Command to be executed
    logical,            intent(in)  :: first_step  !< First time step of a for loop?
    integer,            intent(out) :: ierr        !< Error flag
    
    integer :: i
    
    ierr = 0
    if ( .not. allocated(values) ) then
      write(*,*) 'average_h5_finalize called, but values not allocated!'
      ierr = 99
      return
    end if
    
    ! copy back for writing out
    do i = 1, node_list%n_nodes
      node_list%node(i)%values(:,:,:) = values(:,:,:,i) / total_weight
    end do
    call export_restart(node_list, element_list, 'jorek99999')
    deallocate(values)
    
  end subroutine average_h5_finalize
  
  
  
  
  
  !> Evaluate expressions at a single point.
  subroutine energy_spectrum(command, first_step, ierr)
    
    ! --- Routine parameters
    type(type_command), intent(in)  :: command     !< Command to be executed
    logical,            intent(in)  :: first_step  !< First time step of a for loop?
    integer,            intent(out) :: ierr        !< Error flag
    
    ! --- Local variables
    real*8              :: tmin, tmax, s(n_tor,2), t((n_tor+1)/2,2), left(n_tor,2), right(n_tor,2)
    real*8              :: tleft, tright
    integer             :: i_file, i, j
    character(len=1024) :: filename
    
    ierr = 0
    
    ! --- Some checks
    call check_args(command%n_args,ierr,2);  if ( ierr /= 0 ) return
    call check_step_imported(ierr);          if ( ierr /= 0 ) return
    
    ! --- Preparation
    tmin  = to_float(command%args(1), ierr); if ( ierr /= 0 ) return
    tmax  = to_float(command%args(2), ierr); if ( ierr /= 0 ) return
    
    write(filename,'(9a)') DIR, 'energyspectrum_tmin', trim(real2str(tmin,'(f12.4)')), '_tmax', &
      trim(real2str(tmax,'(f12.4)')), trim(step_range_string(index_now,index_now)), '.dat'
    
    if ( (tmin<xtime(1)) .or. (xtime(index_now)<tmax) .or. (tmax<=tmin) ) then
      write(*,*) 'ERROR in energy_spectrum: Specified time window is invalid.'
      return
    end if
    
    s(:,:) = 0.d0
    
    ! --- sum up for integration (loop over intervals between)
    do i = 2, index_now
      
      ! --- Left value for summation
      if ( (xtime(i-1) < tmin) .and. (tmin <= xtime(i)) ) then
        ! (we are at the beginning of the [tmin,tmax] interval)
        left (:,:) = ( energies(:,:,i-1) * (xtime(i)-tmin) + energies(:,:,i) * (tmin-xtime(i-1)) ) / (xtime(i)-xtime(i-1))
        tleft      = tmin
      else if ( tmin <= xtime(i-1) ) then
        ! (we are somewhere in the middle of the [tmin,tmax] interval)
        left (:,:) = energies(:,:,i-1)
        tleft      = xtime(i-1)
      else
        cycle ! present time point can be skipped
      end if
      
      ! --- Right value for summation
      if ( (xtime(i-1) < tmax) .and. (tmax <= xtime(i)) ) then
        ! (we are at the end of the [tmin,tmax] interval)
        right(:,:) = ( energies(:,:,i-1) * (xtime(i)-tmax) + energies(:,:,i) * (tmax-xtime(i-1)) ) / (xtime(i)-xtime(i-1))
        tright     = tmax
      else if ( xtime(i) < tmax ) then
        ! (we are somewhere in the middle of the [tmin,tmax] interval)
        right(:,:) = energies(:,:,i)
        tright     = xtime(i)
      else
        cycle ! present time point can be skipped
      end if
      
      s(:,:) = s(:,:) + 0.5d0 * ( left(:,:) + right(:,:) ) * ( tright - tleft )
      
    end do
    
    s(:,:) = s(:,:) / (tmax-tmin) ! normalize integral to interval length
    
    ! --- combine cosine and sine components
    t(:,:) = 0.d0
    t(1,:) = s(1,:)
    do i = 1, (n_tor-1)/2
      t(i+1,:) = s(2*i,:) + s(2*i+1,:)
    end do
    
    ! --- write to ascii file
    i_file=133
    open(i_file, file=trim(filename), form='formatted', status='replace', iostat=ierr)
    write(i_file,*) '# energy spectrum'
    write(i_file,*) '# toroidal mode number | magnetic energy spectrum | kinetic energy spectrum'
    do i = 0, (n_tor-1)/2
      write(i_file,'(i7,2es25.15)') i, t(i+1,:)
    end do
    close(i_file)
    
  end subroutine energy_spectrum
  
  
  
  
  
  !> Evaluate expressions at a single point.
  subroutine point(command, first_step, ierr)
    
    ! --- Routine parameters
    type(type_command), intent(in)  :: command     !< Command to be executed
    logical,            intent(in)  :: first_step  !< First time step of a for loop?
    integer,            intent(out) :: ierr        !< Error flag
    
    ! --- Local variables
    real*8  :: R, Z, phi
    integer :: units
    character(len=1024) :: filename
    
    ierr = 0
    
    ! --- Some checks
    call check_args(command%n_args,ierr,2,3);  if ( ierr /= 0 ) return
    call check_step_imported(ierr);          if ( ierr /= 0 ) return
    call check_exprs_selected(ierr);         if ( ierr /= 0 ) return
    
    ! --- Preparation
    R     = to_float(command%args(1), ierr); if ( ierr /= 0 ) return
    Z     = to_float(command%args(2), ierr); if ( ierr /= 0 ) return
    units = get_int_setting('units', ierr)
    
    if (command%n_args == 3) then ! local values
      
      phi = to_float(command%args(3), ierr); if ( ierr /= 0 ) return
      
      write(filename,'(9a)') DIR, 'exprs_at_R', trim(real2str(R)), '_Z', trim(real2str(Z)), '_p',  &
        trim(real2str(phi)), trim(step_range_string(loop_min_step,loop_max_step)), '.dat'
      
      call eval_expr(eq, units, expr_list, pol_pos(node_list,element_list,eq,R=R,Z=Z),             &
        tor_pos(phi=phi), result, ierr)
      
      call reduce_result_to_0d(ierr, result, res0d, 1, 1, 1)
      
    else ! toroidally averaged values
      
      write(filename,'(9a)') DIR, 'exprs_at_R', trim(real2str(R)), '_Z', trim(real2str(Z)),        &
        '_toroidally-averaged', trim(step_range_string(loop_min_step,loop_max_step)), '.dat'
      
      call eval_expr(eq, units, expr_list, pol_pos(node_list,element_list,eq,R=R,Z=Z),             &
        tor_pos(nphi=4*n_plane), result, ierr)
      
      call apply_four_filter(result, simple_filter(n=0), expr_list%n_coord, ierr)
      call reduce_result_to_0d(ierr, result, res0d, 1, 1, 1)
      
    end if
    
    call write_ascii_0d(ierr, eq, expr_list, res0d, FORM_TABLE, header=first_step,                 &
      filename=trim(filename), append=(.not.first_step), blanks=.false.)
    
  end subroutine point
  
  
  
  
  
  !> Evaluate expressions on the midplane.
  subroutine midplane(command, first_step, ierr)
    
    ! --- Routine parameters
    type(type_command), intent(in)  :: command     !< Command to be executed
    logical,            intent(in)  :: first_step  !< First time step of a for loop?
    integer,            intent(out) :: ierr        !< Error flag
    
    ! --- Local variables
    integer :: units, npts, side
    character(len=1024) :: filename, comment, s
    
    ierr = 0
    
    ! --- Some checks
    call check_args(command%n_args,ierr,0,1);  if ( ierr /= 0 ) return
    call check_step_imported(ierr);            if ( ierr /= 0 ) return
    call check_exprs_selected(ierr);           if ( ierr /= 0 ) return
    
    units = get_int_setting('units', ierr)
    npts  = get_int_setting('linepoints', ierr)
    
    if ( command%n_args == 0 ) then
      s = 'midplane'
      side = BOTH_SIDES
    else if ( command%args(1) == 'outer' ) then
      s = 'outer-midplane'
      side = LOWFIELD_SIDE
    else if ( command%args(1) == 'inner' ) then
      s = 'inner-midplane'
      side = HIGHFIELD_SIDE
    else
      write(*,*) 'WARNING: Illegal parameter for command "midplane".'
      ierr = 1
      return
    end if
    
    write(filename,'(4a)') DIR, 'exprs_'//trim(s)//                                                &
      trim(step_range_string(loop_min_step,loop_max_step)), '.dat'
    
    write(comment,'(a,i6.6)') 'time step #', index_now
    
    call midplane_profile(node_list, element_list, eq, units, expr_list, res1d, side, npts,        &
      ierr, filename=trim(filename), append=(.not.first_step), comment=trim(comment) )
    
  end subroutine midplane
  
  
  
  
  
  !> Expressions along a line in the poloidal plane.
  subroutine pol_line(command, first_step, ierr)
    
    ! --- Routine parameters
    type(type_command), intent(in)  :: command     !< Command to be executed
    logical,            intent(in)  :: first_step  !< First time step of a for loop?
    integer,            intent(out) :: ierr        !< Error flag
    
    ! --- Local variables
    real*8  :: Rstart, Zstart, Rend, Zend, phi
    integer :: units, npts
    character(len=1024) :: filename, comment
    
    ierr = 0
    
    ! --- Some checks
    call check_args(command%n_args,ierr,5);  if ( ierr /= 0 ) return
    call check_step_imported(ierr);          if ( ierr /= 0 ) return
    call check_exprs_selected(ierr);         if ( ierr /= 0 ) return
    
    ! --- Preparation
    Rstart = to_float(command%args(1), ierr); if ( ierr /= 0 ) return
    Zstart = to_float(command%args(2), ierr); if ( ierr /= 0 ) return
    Rend   = to_float(command%args(3), ierr); if ( ierr /= 0 ) return
    Zend   = to_float(command%args(4), ierr); if ( ierr /= 0 ) return
    phi    = to_float(command%args(5), ierr); if ( ierr /= 0 ) return
    units  = get_int_setting('units', ierr);      if ( ierr /= 0 ) return
    npts   = get_int_setting('linepoints', ierr); if ( ierr /= 0 ) return
    
    write(filename,'(15a)') DIR, 'exprs_along_line_R', trim(real2str(Rstart)), '..',               &
      trim(real2str(Rend)), '_Z', trim(real2str(Zstart)), '..', trim(real2str(Zend)), '_p',        &
      trim(real2str(phi)), trim(step_range_string(loop_min_step,loop_max_step)), '.dat'
    
    write(comment,'(a,i6.6)') 'time step #', index_now
    
    call pol_lineout(node_list, element_list, eq, units, expr_list, res1d, phi, Rstart, Zstart,    &
      Rend, Zend, npts, ierr, filename, append=(.not.first_step), comment=trim(comment) )
    
  end subroutine pol_line
  
   !> Integrate expressions along a line in the poloidal plane.
  subroutine int_along_pol_line(command, first_step, ierr)

    ! --- Routine parameters
    type(type_command), intent(in)  :: command     !< Command to be executed
    logical,            intent(in)  :: first_step  !< First time step of a for loop?
    integer,            intent(out) :: ierr        !< Error flag

    ! --- Local variables
    real*8  :: Rstart, Zstart, Rend, Zend, phi
    integer :: units, npts
    character(len=1024) :: filename, comment

    ierr = 0

    ! --- Some checks
    call check_args(command%n_args,ierr,5);  if ( ierr /= 0 ) return
    call check_step_imported(ierr);          if ( ierr /= 0 ) return
    call check_exprs_selected(ierr);         if ( ierr /= 0 ) return

    ! --- Preparation
    Rstart = to_float(command%args(1), ierr); if ( ierr /= 0 ) return
    Zstart = to_float(command%args(2), ierr); if ( ierr /= 0 ) return
    Rend   = to_float(command%args(3), ierr); if ( ierr /= 0 ) return
    Zend   = to_float(command%args(4), ierr); if ( ierr /= 0 ) return
    phi    = to_float(command%args(5), ierr); if ( ierr /= 0 ) return
    units  = get_int_setting('units', ierr);      if ( ierr /= 0 ) return
    npts   = get_int_setting('linepoints', ierr); if ( ierr /= 0 ) return

    write(filename,'(15a)') DIR, 'integrate_exprs_along_line_R', trim(real2str(Rstart)), '..',             &
      trim(real2str(Rend)), '_Z', trim(real2str(Zstart)), '..', trim(real2str(Zend)), '_p',                &
      trim(real2str(phi)), trim(step_range_string(loop_min_step,loop_max_step)), '.dat'

    write(comment,'(a,i6.6)') 'time step #', index_now

    call int_along_pol_lineout(node_list, element_list, eq, units, expr_list, sum, phi, Rstart, Zstart,    &
      Rend, Zend, npts, ierr, filename, append=(.not.first_step), comment=trim(comment) )

  end subroutine int_along_pol_line
 
  
  
  
  !> Expressions along a toroidal line.
  subroutine tor_line(command, first_step, ierr)
    
    ! --- Routine parameters
    type(type_command), intent(in)  :: command     !< Command to be executed
    logical,            intent(in)  :: first_step  !< First time step of a for loop?
    integer,            intent(out) :: ierr        !< Error flag
    
    ! --- Local variables
    real*8  :: R, Z, phi_start, phi_end
    integer :: units, npts
    character(len=1024) :: filename, comment
    
    ierr = 0
    
    ! --- Some checks
    call check_args(command%n_args,ierr,4);  if ( ierr /= 0 ) return
    call check_step_imported(ierr);          if ( ierr /= 0 ) return
    call check_exprs_selected(ierr);         if ( ierr /= 0 ) return
    
    ! --- Preparation
    R         = to_float(command%args(1), ierr); if ( ierr /= 0 ) return
    Z         = to_float(command%args(2), ierr); if ( ierr /= 0 ) return
    phi_start = to_float(command%args(3), ierr); if ( ierr /= 0 ) return
    phi_end   = to_float(command%args(4), ierr); if ( ierr /= 0 ) return
    units     = get_int_setting('units', ierr);      if ( ierr /= 0 ) return
    npts      = get_int_setting('linepoints', ierr); if ( ierr /= 0 ) return
    
    write(filename,'(15a)') DIR, 'exprs_along_line_R', trim(real2str(R)), '_Z', trim(real2str(Z)), &
      '_p', trim(real2str(phi_start)), '..', trim(real2str(phi_end)),                              &
      trim(step_range_string(loop_min_step,loop_max_step)), '.dat'
    
    write(comment,'(a,i6.6)') 'time step #', index_now
    
    call tor_lineout(node_list, element_list, eq, units, expr_list, res1d, phi_start, phi_end, R,  &
      Z, npts, ierr, filename, append=(.not.first_step), comment=trim(comment) )
    
  end subroutine tor_line
  
  
  
  
  
  !> Expressions in a rectangular area.
  subroutine rectangle(command, first_step, ierr)
    
    use mod_position, only: pol_pos, tor_pos
    
    ! --- Routine parameters
    type(type_command), intent(in)  :: command     !< Command to be executed
    logical,            intent(in)  :: first_step  !< First time step of a for loop?
    integer,            intent(out) :: ierr        !< Error flag
    
    ! --- Local variables
    real*8  :: Rmin, Rmax, Zmin, Zmax, phi
    integer :: nR, nZ, units
    character(len=1024) :: filename, comment
    
    ierr = 0
    
    ! --- Some checks
    call check_args(command%n_args,ierr,7);  if ( ierr /= 0 ) return
    call check_step_imported(ierr);          if ( ierr /= 0 ) return
    call check_exprs_selected(ierr);         if ( ierr /= 0 ) return
    
    ! --- Preparation
    Rmin      = to_float(command%args(1), ierr); if ( ierr /= 0 ) return
    Rmax      = to_float(command%args(2), ierr); if ( ierr /= 0 ) return
    nR        = to_int  (command%args(3), ierr); if ( ierr /= 0 ) return
    Zmin      = to_float(command%args(4), ierr); if ( ierr /= 0 ) return
    Zmax      = to_float(command%args(5), ierr); if ( ierr /= 0 ) return
    nZ        = to_int  (command%args(6), ierr); if ( ierr /= 0 ) return
    phi       = to_float(command%args(7), ierr); if ( ierr /= 0 ) return
    units     = get_int_setting('units', ierr);      if ( ierr /= 0 ) return
    
    write(filename,'(15a)') DIR, 'exprs_Rmin', trim(real2str(Rmin)), '_Rmax', trim(real2str(Rmax)),&
      '_Zmin', trim(real2str(Zmin)), '_Zmax', trim(real2str(Zmax)), '_phi', trim(real2str(phi)),   &
      trim(step_range_string(index_now,index_now)), '.h5'
      
    comment = 'Output produced by jorek2_postproc command "rectangle"'
    
    call eval_expr(eq, units, expr_list,                                                           &
       pol_pos(node_list,element_list,eq,Rmin=Rmin,Rmax=Rmax,nR=nR,Zmin=Zmin,Zmax=Zmax,nZ=nZ),     &
       tor_pos(phi=phi), result, ierr)
    
    call reduce_result_to_2d(ierr, result, res2d, i1=1)
    call write_hdf5_2d(ierr, expr_list, res2d, trim(filename))
    
    if ( allocated(result) ) deallocate(result)
    if ( allocated(res2d ) ) deallocate(res2d )
    
  end subroutine rectangle
  
  
  
  
  
  !> Toroidally and poloidally averaged expressions.
  subroutine average(command, first_step, ierr)
    
    ! --- Routine parameters
    type(type_command), intent(in)  :: command     !< Command to be executed
    logical,            intent(in)  :: first_step  !< First time step of a for loop?
    integer,            intent(out) :: ierr        !< Error flag
    
    ! --- Local variables
    integer :: units, npts, nsmall
    character(len=1024) :: filename, comment
    type(t_pol_pos_list), save :: pol_pos_list
    type(t_tor_pos_list), save :: tor_pos_list
    
    ierr = 0
    
    ! --- Some checks
    call check_args(command%n_args,ierr,0);  if ( ierr /= 0 ) return
    call check_step_imported(ierr);          if ( ierr /= 0 ) return
    call check_exprs_selected(ierr);         if ( ierr /= 0 ) return
    
    units = get_int_setting('units', ierr)
    npts  = get_int_setting('surfaces', ierr)
    nsmall= get_int_setting('nsmallsteps', ierr)
    
    write(filename,'(4a)') DIR, 'exprs_averaged',                                                  &
      trim(step_range_string(loop_min_step,loop_max_step)), '.dat'
    
    ! ### is nTht and nphi really chosen well???
    pol_pos_list = pol_pos(node_list, element_list, eq, nPsiN=npts, nTht=max(150,6*n_plane),                &
      nsmallsteps=nsmall)
    tor_pos_list = tor_pos(nphi=max(n_plane,2))
    
    call eval_expr(eq, units, expr_list, pol_pos_list, tor_pos_list, result, ierr)
    call apply_four_filter(result, simple_filter(m=0,n=0), expr_list%n_coord, ierr)
    call reduce_result_to_1d(ierr, result, res1d, i1=1, i2=1)
    
    write(comment,'(a,i6.6)') 'time step #', index_now
    
    call write_ascii_1d(ierr, eq, expr_list, res1d, FORM_TABLE, header=.true.,                     &
      filename=trim(filename), append=(.not.first_step), blanks=.true., comment=trim(comment))
    
  end subroutine average
  
  
  
  
  
  !> Output equilibrium information.
  subroutine equil_params(command, first_step, ierr)
    
    ! --- Routine parameters
    type(type_command), intent(in)  :: command     !< Command to be executed
    logical,            intent(in)  :: first_step  !< First time step of a for loop?
    integer,            intent(out) :: ierr        !< Error flag
    
    ! --- Local variables
    integer :: units, i_file
    character(len=1024) :: filename, status, access
    
    ierr = 0
    
    ! --- Some checks
    call check_args(command%n_args,ierr,0);  if ( ierr /= 0 ) return
    call check_step_imported(ierr);          if ( ierr /= 0 ) return
    
    units = get_int_setting('units', ierr)
    
    write(filename,'(4a)') DIR, 'equil_params',                                                    &
      trim(step_range_string(loop_min_step,loop_max_step)), '.dat'
    
    status = 'replace'
    access = 'sequential'
    if ( .not. first_step ) then
      status = 'old'
      access = 'append'
    end if
    i_file=133
    open(i_file, file=trim(filename), form='formatted', status=trim(status), access=trim(access),  &
        iostat=ierr)
    
    if ( first_step ) then
      write(i_file,'(a)') '# time                R_axis              Z_axis              '//       &
        'Psi_axis            R_xpoint(1)         Z_xpoint(1)         Psi_xpoint(1)       '//       &
        'R_xpoint(2)         Z_xpoint(2)         Psi_xpoint(2)       R_lim               '//       &
        'Z_lim               Psi_lim             Psi_bnd'
    end if
    
    write(i_file,'(es20.13,33f20.16)') time_now, eq%R_axis, eq%Z_axis, eq%Psi_axis, eq%R_xpoint(1),&
      eq%Z_xpoint(1), eq%Psi_xpoint(1), eq%R_xpoint(2), eq%Z_xpoint(2), eq%Psi_xpoint(2), eq%R_lim,&
      eq%Z_lim, eq%Psi_lim, eq%Psi_bnd
    
    close(i_file)
    
  end subroutine equil_params




  
  !> Write out SPI related information
  subroutine spi_state(command, first_step, ierr)
    
    ! --- Routine parameters
    type(type_command), intent(in)  :: command     !< Command to be executed
    logical,            intent(in)  :: first_step  !< First time step of a for loop?
    integer,            intent(out) :: ierr        !< Error flag
    
    ! --- Local variables
    integer :: units, i_file, i
    character(len=1024) :: filename, status, access
    real*8 :: R_av, Z_av, phi_av, shard_atoms_left, atoms_left, abl_tot, xx, yy, zz
    
    ierr = 0
    
    ! --- Some checks
    call check_args(command%n_args,ierr,0);  if ( ierr /= 0 ) return
    call check_step_imported(ierr);          if ( ierr /= 0 ) return
    
    units = get_int_setting('units', ierr)
    
    write(filename,'(4a)') DIR, 'spi', trim(step_range_string(loop_min_step,loop_max_step)), '.dat'
    
    status = 'replace'
    access = 'sequential'
    if ( .not. first_step ) then
      status = 'old'
      access = 'append'
    end if
    i_file=133
    open(i_file, file=trim(filename), form='formatted', status=trim(status), access=trim(access),  &
        iostat=ierr)
    
    xx         = 0.d0
    yy         = 0.d0
    zz         = 0.d0
    atoms_left = 0.d0
    abl_tot    = 0.d0
    
    do i = 1, n_SPI
      shard_atoms_left = 4./3.*PI*pellets(i)%spi_radius**3 * pellet_density * 1.d20
      
      atoms_left = atoms_left + shard_atoms_left
      
      xx = xx + pellets(i)%spi_R * cos(pellets(i)%spi_phi) * shard_atoms_left
      yy = yy - pellets(i)%spi_R * sin(pellets(i)%spi_phi) * shard_atoms_left
      zz = zz + pellets(i)%spi_Z                           * shard_atoms_left
      
      abl_tot = abl_tot + pellets(i)%spi_abl
    end do
    
    if ( atoms_left /= 0.d0 ) then
      R_av   = sqrt( xx**2 + yy**2 ) / atoms_left
      Z_av   = Z_av                  / atoms_left
      phi_av = -atan2( yy, xx )      / atoms_left
      if ( phi_av < 0.d0 ) phi_av = phi_av + 2.d0*PI
    else
      R_av   = 0.d0
      Z_av   = 0.d0
      phi_av = 0.d0
    end if
    
    if ( first_step ) then
      write(i_file,'(a)') '# time                R_average           Z_average           '//       &
        'phi_average         atoms_left          ablation_rate'
    end if
    
    write(i_file,'(es20.13,3f20.16,3es20.12)') time_now, R_av, Z_av, phi_av, atoms_left, abl_tot
    
    close(i_file)
    
  end subroutine spi_state




  
  !> Output integrated poloidal current that is normal to the boudary and
  !! toroidal peaking factor (TPF)
  subroutine I_halo_TPF(command, first_step, ierr)
    
    ! --- Routine parameters
    type(type_command), intent(in)  :: command     !< Command to be executed
    logical,            intent(in)  :: first_step  !< First time step of a for loop?
    integer,            intent(out) :: ierr        !< Error flag
    
    ! --- Local variables
    integer :: i_file
    character(len=1024) :: filename, status, access
    real*8 :: I_halo, TPF  
    
    ierr = 0
    
    ! --- Some checks
    call check_args(command%n_args,ierr,0);  if ( ierr /= 0 ) return
    call check_step_imported(ierr);          if ( ierr /= 0 ) return
    
    write(filename,'(4a)') DIR, 'I_halo_TPF',  &
       trim(step_range_string(loop_min_step,loop_max_step)), '.dat'
    
    status = 'replace'
    access = 'sequential'
    if ( .not. first_step ) then
      status = 'old'
      access = 'append'
    end if
    i_file=133

    open(i_file, file=trim(filename), form='formatted', status=trim(status), access=trim(access),  &
        iostat=ierr)
    
    if ( first_step ) then
      write(i_file,'(a)') '#               time            I_halo [MA]             TPF'
    end if
    
    call integrated_normal_bnd_curr(node_list, bnd_node_list, bnd_elm_list, I_halo, TPF)
 
    write(i_file,'(3es20.9)') time_now, I_halo, TPF 
    
    close(i_file)

  end subroutine I_halo_TPF
  
  
  
  
  
  
  !> Output 2d integrals.
  subroutine int2d(command, first_step, ierr)
    
    ! --- Routine parameters
    type(type_command), intent(in)  :: command     !< Command to be executed
    logical,            intent(in)  :: first_step  !< First time step of a for loop?
    integer,            intent(out) :: ierr        !< Error flag
    
    ! --- Local variables
    integer :: units, i_file
    character(len=1024) :: filename, status, access
    real*8 :: aminor, Bgeo, current, beta_p, beta_t, beta_n, density, density_in, density_out,     &
      pressure, pressure_in, pressure_out, heat_src_in, heat_src_out, part_src_in, part_src_out
    real*8 :: fact_mu_zero, fact_ne
    
    ierr = 0
    
    ! --- Some checks
    call check_args(command%n_args,ierr,0);  if ( ierr /= 0 ) return
    call check_step_imported(ierr);          if ( ierr /= 0 ) return
    
    units = get_int_setting('units', ierr)
    
    write(filename,'(4a)') DIR, 'int2d',                                                           &
      trim(step_range_string(loop_min_step,loop_max_step)), '.dat'
    
    status = 'replace'
    access = 'sequential'
    if ( .not. first_step ) then
      status = 'old'
      access = 'append'
    end if
    i_file=133
    open(i_file, file=trim(filename), form='formatted', status=trim(status), access=trim(access),  &
        iostat=ierr)
    
    if ( first_step ) then
      write(i_file,'(a)') '#               time            pressure         pressure_in        ' //&
        'pressure_out             density          density_in         density_out              ' //&
        'beta_n              beta_t              beta_p            current'
    end if
    
    if ( units == SI_UNITS ) then
      fact_mu_zero = MU_zero
      fact_ne      = central_density * 1.d20
    else
      fact_mu_zero = 1.d0
      fact_ne      = 1.d0
    end if
    
    call integrals(node_list, element_list, eq%R_axis, eq%Z_axis, eq%psi_axis, eq%R_xpoint,        &
      eq%Z_xpoint, eq%psi_xpoint, eq%psi_lim, aminor, Bgeo, current, beta_p, beta_t, beta_n,       &
      density, density_in, density_out, pressure, pressure_in, pressure_out, heat_src_in,          &
      heat_src_out, part_src_in, part_src_out)
    
    write(i_file,'(33es20.13)') time_now, pressure/fact_mu_zero, pressure_in/fact_mu_zero,         &
      pressure_out/fact_mu_zero, density*fact_ne, density_in*fact_ne, density_out*fact_ne, beta_n, &
      beta_t, beta_p, current
    
    close(i_file)
    
  end subroutine int2d
 



 
  
  !> Output 3d integrals.
  subroutine int3D(command, first_step, ierr)
    
    use mod_integrals3D
    use mpi_mod

    ! --- Routine parameters
    type(type_command), intent(in)  :: command     !< Command to be executed
    logical,            intent(in)  :: first_step  !< First time step of a for loop?
    integer,            intent(out) :: ierr        !< Error flag
    
    ! --- Local variables
    integer :: i_file, i, units, my_id
    integer :: required, provided, StatInfo
    character(len=1024) :: filename, status, access
    real*8, allocatable :: res(:)
    character(len=23)   :: s

    ierr = 0
    my_id=0

    ! --- Initialize MPI
#ifdef FUNNELED
    required = MPI_THREAD_FUNNELED
#else
    required = MPI_THREAD_MULTIPLE
#endif

    if (first_step)  call MPI_Init_thread(required, provided, StatInfo)
   
    ! --- Some checks
    call check_args(command%n_args,ierr,0,1);  if ( ierr /= 0 ) return
    call check_step_imported(ierr);            if ( ierr /= 0 ) return
    call check_exprs_selected(ierr);           if ( ierr /= 0 ) return
    units = get_int_setting('units', ierr)

    allocate(res(expr_list%n_expr+1))
    res = 0.d0   
 
    write(filename,'(4a)') DIR, 'integrals3D',  &
       trim(step_range_string(loop_min_step,loop_max_step)), '.dat'
    
    status = 'replace'
    access = 'sequential'
    if ( .not. first_step ) then
      status = 'old'
      access = 'append'
    end if
    i_file=133

    open(i_file, file=trim(filename), form='formatted', status=trim(status), access=trim(access),  &
        iostat=ierr)
    
    if ( first_step ) then
      write(i_file,'(a)',advance='no') '# time                   '
      do i = 1, expr_list%n_expr
        s = trim(expr_list%expr(i)%name)
        write(i_file,'(a)',advance='no') s
      end do
      write(i_file,'(a)')
    end if
    close(i_file)
 
   call int3d_new(my_id, node_list, element_list, bnd_node_list, bnd_elm_list, expr_list, res, units)        

   call write_ascii_0d(ierr, eq, expr_list, res, FORM_TABLE, header=.false.,                   &
     filename=filename, append=.true., blanks=.false.)
   
  end subroutine int3D
  
  
  
 


  !> Output current density normal to the jorek boundary as a function of Rbnd
  !! and Zbnd
  subroutine jnorm_bnd_RZ(command, ierr)
    
    ! --- Routine parameters
    type(type_command), intent(in)  :: command     !< Command to be executed
    integer,            intent(out) :: ierr        !< Error flag
    
    ! --- Local variables
    logical   :: bool_si_units 
    integer   :: i_plane, units
    
    ierr = 0
    
    ! --- Some checks
    call check_args(command%n_args,ierr,1);  if ( ierr /= 0 ) return
    call check_step_imported(ierr);          if ( ierr /= 0 ) return
    units = get_int_setting('units', ierr)

    i_plane  = to_int(command%args(1), ierr); if ( ierr /= 0 ) return

    if ((i_plane <= 0) .or. (i_plane > n_plane) ) then
      write(*,*) 'Incorrect i_plane, note that    0 < i_plane <= n_plane'
      return
    endif

    if ( units == SI_UNITS ) then
      bool_si_units = .true. 
    else
      bool_si_units = .false.
    end if
 
    call normal_bnd_curr(node_list, element_list, bnd_node_list, &
                            bnd_elm_list, i_plane, bool_si_units) 

  end subroutine jnorm_bnd_RZ
  
  


 
  !> Output the q-profile as a function of Psi_N
  recursive subroutine qprofile(command, first_step, ierr)
  
    ! --- Routine parameters
    type(type_command), intent(in)  :: command     !< Command to be executed
    logical,            intent(in)  :: first_step  !< First time step of a for loop?
    integer,            intent(out) :: ierr        !< Error flag

    
    ! --- Local variables
    integer                  :: k, k2, npts
    real*8, allocatable      :: q(:), rad(:)
    type (type_surface_list) :: surface_list
    character(len=1024)      :: filename, comment
    type(t_expr_list)        :: tmp_expr_list
    
    ierr = 0
    
    ! --- Some checks
    call check_args(command%n_args,ierr,0);  if ( ierr /= 0 ) return
    call check_step_imported(ierr);          if ( ierr /= 0 ) return
    
    npts  = get_int_setting('surfaces', ierr)
    
    write(filename,'(4a)') DIR, 'qprofile', trim(step_range_string(loop_min_step,loop_max_step)),  &
      '.dat'
    
    ! --- Find flux surfaces and determine q-profile
    surface_list%n_psi = npts
    allocate( surface_list%psi_values(npts), q(npts), rad(npts) )
    do k = 1, npts
      surface_list%psi_values(k) = eq%psi_axis + (eq%psi_bnd - eq%psi_axis) * real(k-1)/real(npts-1)
    end do
    call find_flux_surfaces(0,xpoint, xcase, node_list, element_list, surface_list)
    call determine_q_profile(node_list, element_list, surface_list, eq%psi_axis, eq%psi_xpoint,    &
      eq%Z_xpoint, q, rad)
    
    ! --- Write out q-profile versus Psi_n
    tmp_expr_list%n_expr = 0
    tmp_expr_list%expr(1)%name = 'Psi_n'
    tmp_expr_list%expr(2)%name = 'q'
    write(comment,'(a,i6.6)') 'time step #', index_now
    allocate(res1d(npts-2,2))
    do k2 = 1, npts-2
      k = k2 + 1 ! to avoid first and last point of q-profile which often is bad
      res1d(k2,:) = (/ get_psi_n(eq, surface_list%psi_values(k)), q(k) /)
    end do
    
    call write_ascii_1d(ierr, eq, tmp_expr_list, res1d, FORM_TABLE, header=.true.,                 &
      filename=trim(filename), append=(.not.first_step), blanks=.true., comment=trim(comment))
    
    ! --- Clean up.
    if ( allocated(surface_list%psi_values)    ) deallocate(surface_list%psi_values)
    if ( allocated(surface_list%flux_surfaces) ) deallocate(surface_list%flux_surfaces)
    if ( allocated(q)                          ) deallocate(q)
    if ( allocated(rad)                        ) deallocate(rad)
    
  end subroutine qprofile




  
  !> Output q over time at a certain psi_n
  subroutine q_at_given_psin(command, first_step, ierr)
    
    ! --- Routine parameters
    type(type_command), intent(in)  :: command     !< Command to be executed
    logical,            intent(in)  :: first_step  !< First time step of a for loop?
    integer,            intent(out) :: ierr        !< Error flag
    
    ! --- Local variables
    integer :: i_file
    character(len=1024) :: filename, status, access
    real*8 :: t_norm, psin, q_psin(2), rad(2) 
    
    type (type_surface_list) :: surface_list
 
    ierr = 0
    
    ! --- Some checks
    call check_args(command%n_args,ierr,1);  if ( ierr /= 0 ) return
    call check_step_imported(ierr);          if ( ierr /= 0 ) return
    
    psin  = to_float(command%args(1), ierr); if ( ierr /= 0 ) return
   
    write(filename,'(5a)') DIR, 'q_at_psin_', trim(real2str(psin,'(f12.4)')), &
       trim(step_range_string(loop_min_step,loop_max_step)), '.dat'
    
    status = 'replace'
    access = 'sequential'
    if ( .not. first_step ) then
      status = 'old'
      access = 'append'
    end if
    i_file=133

    open(i_file, file=trim(filename), form='formatted', status=trim(status), access=trim(access),  &
        iostat=ierr)
    
    if ( first_step ) then
      write(i_file,'(a)') '#               time            q'
    end if
    
    ! --- Find flux surfaces and determine q-profile
    surface_list%n_psi = 2 
    allocate( surface_list%psi_values(2) )
    surface_list%psi_values(1) = eq%psi_axis + (eq%psi_bnd - eq%psi_axis) * 0.2d0
    surface_list%psi_values(2) = eq%psi_axis + (eq%psi_bnd - eq%psi_axis) * psin

    call find_flux_surfaces(0,xpoint, xcase, node_list, element_list, surface_list)
    call determine_q_profile(node_list, element_list, surface_list, eq%psi_axis, eq%psi_xpoint,    &
      eq%Z_xpoint, q_psin, rad)
    
    write(i_file,'(2es20.13)') time_now, q_psin(2) 
    
    close(i_file)

    ! --- Clean up.
    if ( allocated(surface_list%psi_values)    ) deallocate(surface_list%psi_values)
    if ( allocated(surface_list%flux_surfaces) ) deallocate(surface_list%flux_surfaces)
    
  end subroutine q_at_given_psin
  
  
  

  
  !> Output the flux surfaces.
  recursive subroutine fluxsurfaces(command, ierr)
  
    ! --- Routine parameters
    type(type_command), intent(in)  :: command     !< Command to be executed
    integer,            intent(out) :: ierr        !< Error flag
    
    
    ! --- Local variables
    integer                  :: i, j, i_elm, npts, ip, nplot, i_file
    type (type_surface_list) :: surface_list
    character(len=1024)      :: filename, comment
    type(t_expr_list)        :: tmp_expr_list
    real*8                   :: psi_min, psi_max, psi_min2, psi_max2, ss1, dss1, ss2, dss2, tt1,   &
      dtt1, tt2, dtt2, u, si, dsi, ti, dti, R, R_s, R_t, R_st, R_ss, R_tt, Z, Z_s, Z_t, Z_st, Z_ss,&
      Z_tt
    
    ierr = 0
    
    ! --- Some checks
    call check_args(command%n_args,ierr,0);  if ( ierr /= 0 ) return
    call check_step_imported(ierr);          if ( ierr /= 0 ) return
    
    npts  = get_int_setting('surfaces', ierr)
    
    write(filename,'(4a)') DIR, 'fluxsurfaces', trim(step_range_string(index_start,index_start)),  &
      '.dat'
    
    ! --- Find minimum and maximum psi values
    psi_min=+1.d99
    psi_max=-1.d99
    do i_elm = 1, element_list%n_elements
      call psi_minmax(node_list, element_list, i_elm, psi_min2, psi_max2)
      psi_min = min(psi_min,psi_min2)
      psi_max = max(psi_max,psi_max2)
    end do
    psi_min = psi_min + 0.001*(psi_max-psi_min)
    psi_max = psi_max - 0.001*(psi_max-psi_min)
    
    ! --- Find flux surfaces
    surface_list%n_psi = npts
    allocate( surface_list%psi_values(npts) )
    do i = 1, npts
      surface_list%psi_values(i) = psi_min + (psi_max-psi_min) * real(i-1)/real(npts-1)
    end do
    call find_flux_surfaces(0,xpoint, xcase, node_list, element_list, surface_list)
    
    ! --- Write out flux surfaces
    nplot  = 5
    i_file = 111
    call open_ascii_file(ierr, i_file, filename, .false.)
    do i = 1, npts
      
      ! --- Loop over all segments of this flux surface
      do j=1,surface_list%flux_surfaces(i)%n_pieces
        
        ! --- Bezier element, in which the current flux surface segment is located
        i_elm = surface_list%flux_surfaces(i)%elm(j)
        ss1  = surface_list%flux_surfaces(i)%s(1,j)
        dss1 = surface_list%flux_surfaces(i)%s(2,j)
        ss2  = surface_list%flux_surfaces(i)%s(3,j)
        dss2 = surface_list%flux_surfaces(i)%s(4,j)
        
        tt1  = surface_list%flux_surfaces(i)%t(1,j)
        dtt1 = surface_list%flux_surfaces(i)%t(2,j)
        tt2  = surface_list%flux_surfaces(i)%t(3,j)
        dtt2 = surface_list%flux_surfaces(i)%t(4,j)
        
        ! --- Loop over nplot points in a flux surface segment
        do ip = 1, nplot
          u = -1. + 2.*float(ip-1)/float(nplot-1)
          
          ! --- Determine s and t values of the current point inside element i_elm
          call CUB1D(ss1, dss1, ss2, dss2, u, si, dsi)
          call CUB1D(tt1, dtt1, tt2, dtt2, u, ti, dti)
          
          ! --- Determine (R,Z)-coordinates of the current point on the current flux surface
          call interp_RZ(node_list, element_list, i_elm, si, ti, R, R_s, R_t, R_st, R_ss, R_tt, &
            Z, Z_s, Z_t, Z_st, Z_ss, Z_tt)
            
          ! --- Write out the (R,Z)-coordinates
          write(i_file,'(2ES16.7)') R, Z
        end do
        
        write(i_file,*)
        write(i_file,*)
      
      end do
      
    end do
    
    close(i_file)
    
    ! --- Clean up.
    if ( allocated(surface_list%psi_values)    ) deallocate(surface_list%psi_values)
    if ( allocated(surface_list%flux_surfaces) ) deallocate(surface_list%flux_surfaces)
    
  end subroutine fluxsurfaces
  
  
  
  
  
  !> Output the separatrix.
  recursive subroutine separatrix(command, ierr)
  
    ! --- Routine parameters
    type(type_command), intent(in)  :: command     !< Command to be executed
    integer,            intent(out) :: ierr        !< Error flag
    
    ! --- Local variables
    integer                  :: i, j, i_elm, npts, ip, nplot, i_file
    type (type_surface_list) :: surface_list
    character(len=1024)      :: filename, comment
    type(t_expr_list)        :: tmp_expr_list
    real*8                   :: psi_min, psi_max, psi_min2, psi_max2, ss1, dss1, ss2, dss2, tt1,   &
      dtt1, tt2, dtt2, u, si, dsi, ti, dti, R, R_s, R_t, R_st, R_ss, R_tt, Z, Z_s, Z_t, Z_st, Z_ss,&
      Z_tt
    
    ierr = 0
    
    ! --- Some checks
    call check_args(command%n_args,ierr,0);  if ( ierr /= 0 ) return
    call check_step_imported(ierr);          if ( ierr /= 0 ) return
    
    write(filename,'(4a)') DIR, 'separatrix', trim(step_range_string(index_start,index_start)),    &
      '.dat'
    
    ! --- Find flux surfaces
    npts = 1
    surface_list%n_psi = 1
    allocate( surface_list%psi_values(1) )
    surface_list%psi_values(1) = eq%psi_bnd
    call find_flux_surfaces(0,xpoint, xcase, node_list, element_list, surface_list)
    
    ! --- Write out flux surfaces
    nplot  = 5
    i_file = 111
    call open_ascii_file(ierr, i_file, filename, .false.)
    i = 1
    ! --- Loop over all segments of this flux surface
    do j=1,surface_list%flux_surfaces(i)%n_pieces
      
      ! --- Bezier element, in which the current flux surface segment is located
      i_elm = surface_list%flux_surfaces(i)%elm(j)
      ss1  = surface_list%flux_surfaces(i)%s(1,j)
      dss1 = surface_list%flux_surfaces(i)%s(2,j)
      ss2  = surface_list%flux_surfaces(i)%s(3,j)
      dss2 = surface_list%flux_surfaces(i)%s(4,j)
      
      tt1  = surface_list%flux_surfaces(i)%t(1,j)
      dtt1 = surface_list%flux_surfaces(i)%t(2,j)
      tt2  = surface_list%flux_surfaces(i)%t(3,j)
      dtt2 = surface_list%flux_surfaces(i)%t(4,j)
      
      ! --- Loop over nplot points in a flux surface segment
      do ip = 1, nplot
        u = -1. + 2.*float(ip-1)/float(nplot-1)
        
        ! --- Determine s and t values of the current point inside element i_elm
        call CUB1D(ss1, dss1, ss2, dss2, u, si, dsi)
        call CUB1D(tt1, dtt1, tt2, dtt2, u, ti, dti)
        
        ! --- Determine (R,Z)-coordinates of the current point on the current flux surface
        call interp_RZ(node_list, element_list, i_elm, si, ti, R, R_s, R_t, R_st, R_ss, R_tt,      &
          Z, Z_s, Z_t, Z_st, Z_ss, Z_tt)
          
        ! --- Write out the (R,Z)-coordinates
        write(i_file,'(2ES16.7)') R, Z
      end do
      
      write(i_file,*)
      write(i_file,*)
    
    end do
    
    close(i_file)
    
    ! --- Clean up.
    if ( allocated(surface_list%psi_values)    ) deallocate(surface_list%psi_values)
    if ( allocated(surface_list%flux_surfaces) ) deallocate(surface_list%flux_surfaces)
    
  end subroutine separatrix
  
  
  
  
  
  !> Perform a 2D Fourier analysis (in straight field line coordinates).
  subroutine four2d(command, ierr)
    
    ! --- Routine parameters
    type(type_command), intent(in)  :: command     !< Command to be executed
    integer,            intent(out) :: ierr        !< Error flag
    
    ! --- Local variables
    integer :: units, npts, nsmall, nmaxstep, n_thetastar
    real*8  :: radial_range(2), delta_phi
    character(len=1024) :: filename_start
    type(t_pol_pos_list), save :: pol_pos_list
    type(t_tor_pos_list), save :: tor_pos_list
    
    ierr = 0
    
    ! --- Some checks
    call check_args(command%n_args,ierr,0);  if ( ierr /= 0 ) return
    call check_step_imported(ierr);          if ( ierr /= 0 ) return
    call check_exprs_selected(ierr);         if ( ierr /= 0 ) return
    
    units  = get_int_setting('units', ierr)
    npts   = get_int_setting('surfaces', ierr)
    nsmall = get_int_setting('nsmallsteps', ierr)
    nmaxstep = get_int_setting('nmaxsteps', ierr)
    delta_phi = get_float_setting('deltaphi', ierr)
    radial_range(1) = get_float_setting('rad_range_min', ierr)
    radial_range(2) = get_float_setting('rad_range_max', ierr)
    n_thetastar = get_int_setting('nTht', ierr)
    

    write(filename_start,'(3a)') DIR, 'exprs_four2d', trim(step_range_string(index_now,index_now))
    write(*,*) 'Input parameters set:'
    write(*,*) 'units        =', units
    write(*,*) 'surfaces     =', npts
    write(*,*) 'nsmallsteps  =', nsmall
    write(*,*) 'nmaxsteps    =', nmaxstep
    write(*,*) 'deltaphi     =', delta_phi
    write(*,*) 'rad_range    =', radial_range
    write(*,*) 'n_thetastar  =', n_thetastar
    
    call fourier_analysis(node_list, element_list, eq, units, expr_list, cp, npts, ierr,           &
      filename_start, OUTP_ABS_VALUE, nsmallsteps=nsmall, nmaxsteps=nmaxstep, deltaphi=delta_phi,  &
      rad_range=radial_range, nTht=n_thetastar)
    
  end subroutine four2d
  
  
  
  
  
  !> Output magnetic field for the Gourdon code.
  subroutine gourdon(command, first_step, ierr)
    
    ! --- Routine parameters
    type(type_command), intent(in)  :: command     !< Command to be executed
    logical,            intent(in)  :: first_step  !< First time step of a for loop?
    integer,            intent(out) :: ierr        !< Error flag
    
    ! --- Local variables
    integer :: units, npts, n_R, n_Z, n_phi, i, j, k
    character(len=1024) :: filename
    real*8 :: R_min, R_max, Z_min, Z_max, R_max2, Z_max2, phi_max2, fact_phi, fact_btor, fact_bpol,&
      tmp
    real*8, allocatable :: field(:,:,:,:)
    type(t_pol_pos_list), save :: pol_pos_list
    type(t_tor_pos_list), save :: tor_pos_list
    type(t_expr_list),    save :: tmp_expr_list
    
    ierr = 0
    
    ! --- Some checks
    call check_args(command%n_args,ierr,10); if ( ierr /= 0 ) return
    call check_step_imported(ierr);          if ( ierr /= 0 ) return
    
    ! --- Preparation
    R_min = to_float(command%args(1), ierr); if ( ierr /= 0 ) return
    R_max = to_float(command%args(2), ierr); if ( ierr /= 0 ) return
    n_R   = to_int  (command%args(3), ierr); if ( ierr /= 0 ) return
    Z_min = to_float(command%args(4), ierr); if ( ierr /= 0 ) return
    Z_max = to_float(command%args(5), ierr); if ( ierr /= 0 ) return
    n_Z   = to_int  (command%args(6), ierr); if ( ierr /= 0 ) return
    n_phi = to_int  (command%args(7), ierr); if ( ierr /= 0 ) return
    !   --- The following three parameters are correction factors which allow to
    !     transform the field between different coordinate systems.
    fact_phi  = to_float(command%args(8),  ierr); if ( ierr /= 0 ) return
    fact_btor = to_float(command%args(9),  ierr); if ( ierr /= 0 ) return
    fact_bpol = to_float(command%args(10), ierr); if ( ierr /= 0 ) return
    
    write(filename,'(3a)') DIR, 'gourdon', trim(step_range_string(index_now,index_now))
    
    ! --- Take into account that the last points are not included in Gourdon format!
    R_max2   = R_max   - (R_max-R_min) / real(n_R)
    Z_max2   = Z_max   - (Z_max-Z_min) / real(n_Z)
    phi_max2 = 2.d0*pi - (2.d0*pi)     / real(n_phi)    * fact_phi
    
    ! --- Make sure that positions outside the JOREK domain get bfield=0.
    tmp = expr_outside_value
    expr_outside_value = 0.d0
    
    ! --- Calculate field components.
    if ( first_step ) then ! (Positions remain unchanged for all time steps, compute only once)
      call create_pol_pos(pol_pos_list, ierr, node_list, element_list, eq, Rmin=R_min, Rmax=R_max2,&
        nR=n_R, Zmin=Z_min, Zmax=Z_max2, nZ=n_Z)
      tor_pos_list  = tor_pos(phistart=0.d0, phiend=phi_max2, nphi=n_phi)
      tmp_expr_list = exprs((/'B_tor', 'B_R  ', 'B_Z  '/), 3)
    end if
    call eval_expr(eq, JOREK_UNITS, tmp_expr_list, pol_pos_list, tor_pos_list, result, ierr)
    if ( fact_btor /= 1.d0 ) result(:,:,:,1  ) = result(:,:,:,1  ) * fact_btor
    if ( fact_bpol /= 1.d0 ) result(:,:,:,2:3) = result(:,:,:,2:3) * fact_bpol
    allocate(field(3,0:n_R-1,0:n_Z-1,0:n_phi-1))
    !    --- Array transform for file output...
    do k = 0, n_Z-1
      do j = 0, n_R-1
        do i = 1, 3
          field(i,j,k,0:n_phi-1) = result(1:n_phi,j+1,k+1,i)
        end do
      end do
    end do
    
    ! --- Write out the data.
    open(122, file=filename, status='replace', action='write', form='unformatted', iostat=ierr)
    if (ierr /= 0) then
      write(*,*) 'ERROR in routine gourdon: Creating file "', trim(filename), '" failed.'
      return
    end if
    write(122) real(n_phi), real(n_R), real(n_Z), real(3), real(1), real(n_phi)
    write(122) (R_max+R_min)/2.d0, (Z_max+Z_min)/2.d0, (R_max-R_min)/2.d0, (Z_max-Z_min)/2.d0
    do k = 0, n_phi - 1
      write(122) field(:,:,:,k)
    end do
    close(122)
    
    ! --- Clean up.
    expr_outside_value = tmp ! restore previous value
    if ( allocated(field ) ) deallocate(field)
    
  end subroutine gourdon
  
  
  
  
  
  !> Select JOREK normalized units.
  subroutine select_jorek_units(command, ierr)
    
    ! --- Routine parameters
    type(type_command), intent(in)  :: command     !< Command to be executed
    integer,            intent(out) :: ierr        !< Error flag
    
    ierr = 0
    
    ! --- Some checks
    call check_args(command%n_args,ierr,0); if ( ierr /= 0 ) return
    
    call set_setting('units', '0', ierr)
    
  end subroutine select_jorek_units
  
  
  
  
  
  !> Select SI units.
  subroutine select_si_units(command, ierr)
    
    ! --- Routine parameters
    type(type_command), intent(in)  :: command     !< Command to be executed
    integer,            intent(out) :: ierr        !< Error flag
    
    ierr = 0
    
    ! --- Some checks
    call check_args(command%n_args,ierr,0); if ( ierr /= 0 ) return
    
    call set_setting('units', '1', ierr)
    
  end subroutine select_si_units
  
  
  
  
  
  !> Output the computational grid.
  subroutine grid(command, ierr)
    
    ! --- Routine parameters
    type(type_command), intent(in)  :: command     !< Command to be executed
    integer,            intent(out) :: ierr        !< Error flag
    
    ierr = 0
    
    ! --- Some checks
    call check_args(command%n_args,ierr,0); if ( ierr /= 0 ) return
    
    call plot_grid(node_list, element_list, bnd_elm_list, bnd_node_list, .true., .false.,          &
      trim(step_range_string(index_now,index_now)))
    
    call system('mv '//'grid_'//trim(step_range_string(index_now,index_now))//'.dat '//DIR)
    
  end subroutine grid
  
  
  
  
  
end module exec_commands
