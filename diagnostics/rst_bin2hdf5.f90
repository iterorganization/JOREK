!> Program to convert a JOREK2 BINARY restart file into HDF5 restart file
program RST_convert_bin2hdf5

  use data_structure
  use phys_module
  ! Argument parsing
  use cla
  use mod_import_restart
  use mod_export_restart
  use mod_boundary, only: boundary_from_grid
  use vacuum
  use vacuum_response
  use vacuum_equilibrium
  use nodes_elements
  use mpi_mod

  implicit none

  integer :: ierr, i, provided, StatInfo, my_id
  character(len=80) :: filein, fileout
  logical :: verbose, file_exists

#ifndef USE_HDF5
#error " Should be compiled with -DUSE_HDF5"
#endif

  call MPI_Init_thread(MPI_THREAD_FUNNELED, provided, StatInfo)
  call init_threads()
  call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)
  
  call initialise_parameters(0, '__NO_FILENAME__')
  call vacuum_init(0, freeboundary_equil, freeboundary, resistive_wall)
  call update_time_evol_params()

  ! Parse command line arguments
  call cla_init
  call pcla_register('filename', 'name of the restart file to convert',  cla_char, 'jorek_restart.rst')    
  call cla_register('-v','--verbose','enable verbose output', cla_flag,'n')
  call cla_validate("rst_bin2hdf5")
  call cla_get('filename',filein)
  verbose = cla_key_present('--verbose')

  ! Create output filename
  fileout = filein(1:index(filein,'.rst',.true.)) // 'h5' ! .true. searches backwards

  ! --- Initialize mode and mode_type arrays
  call det_modes()
  rst_format = 0

  ! --- Check for presence of the restart file
  inquire(file=filein, exist=file_exists)
  if (.not. file_exists) then
    write(*,*) "File " // trim(filein) // " not found", ""
    call cla_help('rst_bin2hdf5')
    call exit(1)
  endif

  ! --- Read the restart binary file
  if (verbose) write (6,*) " =============> rst_bin2hdf5 for filename = ",filein
  call import_binary_restart(node_list, element_list, filein, rst_format, ierr)

  index_now = index_start
  t_now     = t_start
  visco     = visco_rst
  visco_par = visco_par_rst
  eta       = eta_rst

  call boundary_from_grid(node_list, element_list, bnd_node_list, bnd_elm_list, .false.)
  
  if ( freeboundary ) then
    call get_vacuum_response(0, node_list, bnd_elm_list, bnd_node_list, freeboundary_equil,  &
      resistive_wall)
    call update_response(my_id, tstep, resistive_wall)
    call import_external_fields('coil_field.dat', 0)
    if ( .not. wall_curr_initialized ) call init_wall_currents(0, resistive_wall)
  end if

  ! -- Write the HDF5 restart file
  if (verbose) write (6,*) " =============> rst_bin2hdf5, write HDF5 file = ",fileout
  call export_hdf5_restart(node_list, element_list, fileout)

end program RST_convert_bin2hdf5
