program process_hdf5_jorek

  USE HDF5 ! This module contains all necessary modules
  use hdf5_io_module

  implicit none
  include 'mpif.h'

  ! --- Executable arguments variables
  integer :: n_arg
  character*250 :: str_input, str_input2

  ! --- HDF5 variables including filename
  character*250  :: filename, filename_data, filename_path, hdf5_dataname
  integer        :: filename_min, filename_max, filename_step, filename_index
  logical        :: filename_list, hdf5_compact
  character*250  :: filename_data_multiple
  integer(HID_T) :: file_id, data_id
  integer        :: error
  
  ! --- Picture variables
  character*250 :: filename_picture
  character*50  :: char_min, char_max
  character*1   :: rgb(3)
  integer       :: colormap
  logical       :: image_output
  integer       :: itmp, icnt
  
  ! --- MPI variables
  integer :: n_cpu, my_id, ierr, nrecv, nsend
  integer :: status(MPI_STATUS_SIZE)
  
  ! --- Data variables
  integer :: n_tor, n_period, n_elements, n_nodes, n_var, n_data, index_now, jorek_model
  real*8  :: phi_angle
  real*8  :: normalise_min, normalise_max, normalise_data
  real*8  :: central_density, F0, delta_t, time_now
  integer :: i_var(100) ! use a large buffer
  integer :: i_times, n_times
  integer, allocatable :: elm_vertex(:,:)
  real*8,  allocatable :: R_axis_t(:), Z_axis_t(:), psi_axis_t(:), psi_bnd_t(:), Z_xpoint_t(:,:)
  real*8,  allocatable :: modes(:), elm_size(:,:,:), RZ_nodes(:,:,:), values(:,:,:,:), deltas(:,:,:,:), time(:)
  real*8,  allocatable :: HZ(:), HZ_p(:), HZ_pp(:)
  real*8,  allocatable :: my_var(:,:,:,:), my_var3D(:,:,:), my_var3D_time(:,:,:,:,:), my_var_time(:,:,:,:,:)
  real*8,  allocatable :: time_record(:)
  real*8,  allocatable :: Rgrid(:,:), Zgrid(:,:)
  real*8,  allocatable :: radius_grid(:,:), theta_grid(:,:), normalised_radius_grid(:,:)
  real*8,  allocatable :: Xgrid3D(:), Ygrid3D(:), Zgrid3D(:)
  character*50 :: variable_names(100) ! use a large buffer
  
  ! --- Grid variables
  integer :: nx_3D, ny_3D, nz_3D
  integer :: nx_pix, ny_pix
  real*8  :: resolution, limit_buffer
  real*8  :: Rmin, Rmax, Zmin, Zmax
  real*8  :: Rpix_min, Rpix_max
  real*8  :: Xpix_min, Xpix_max
  real*8  :: Ypix_min, Ypix_max
  real*8  :: Zpix_min, Zpix_max
  logical :: save_pixels, use_pixel_file, use_3D_grid, use_r_theta
  real*8  :: radius_max, radius_find, radius, theta
  real*8, parameter    :: PI=3.141592653589793
  integer :: pix_start, pix_end, pix_delta
  integer :: local_pix_start, local_pix_end
  integer, allocatable :: i_elm_save(:,:)
  real*8,  allocatable :: s_save(:,:), t_save(:,:)
  integer, allocatable :: Xindex3D(:), Yindex3D(:),Zindex3D(:), index3D(:,:,:)

  ! --- VTK variables
  integer                 :: ivtk, vtk_n_cells, vtk_points_per_cell
  integer(kind=C_INT8_T)  :: uint8
  integer(kind=C_INT32_T) :: int32
  !integer :: uint8 ! changed because not using Mitsuba here and compile errors...
  !integer :: int32
  real*4                  :: float32
  real*4,allocatable      :: vtk_xyz (:,:), vtk_scalars(:,:), vtk_scalars3D(:,:,:,:)
  integer,allocatable     :: vtk_cells(:,:)
  integer                 :: etype
  character               :: buffer*80, lf*1, str1*10, str2*10, str3*3
  
  ! --- Photon Emissivity Coeff (PEC) variables
  integer               :: PEC_size, PEC_index_Ne, PEC_index_Te, PEC_index, k_pec
  character*50          :: line
  real*8,allocatable    :: PEC_dens(:), PEC_temp(:), PEC(:)
  logical               :: D_alpha
  real*8                :: D_alpha_emission

  ! --- Other variables
  integer :: i, j, k, l, i_elm, i_tor, i_data, i_check, j_check, k_check, n_tor_min, n_tor_max, k_tor
  real*8  :: xjac
  real*8  :: R_tmp, R_out, RR_s, RR_t, s_out
  real*8  :: Z_tmp, Z_out, ZZ_s, ZZ_t, t_out
  real*8  :: x_tmp3D, y_tmp3D, z_tmp3D
  real*8  :: psi, psi_s, psi_t, psi_ss, psi_tt, psi_st, psi_p, delta_psi, psi_R, psi_Z
  real*8  :: phi, phi_s, phi_t, phi_ss, phi_tt, phi_st, phi_p, delta_phi, phi_R, phi_Z
  real*8  :: pp,  pp_s,  pp_t,  pp_ss,  pp_tt,  pp_st,  pp_p,  delta_pp
  real*8  :: rho, rho_n, Te, Ti
  real*8  :: BR, BZ, Bp, ER, EZ, Ep, Epar
  real*8  :: R_axis, Z_axis, psi_axis, psi_bnd, Z_xpoint(2), psi_max
  integer :: ifail
  real*8  :: progress
  logical :: txt_data, bin_data, hdf5_data, fourier
  real*8  :: rho_0, mu_0, t_norm, eV2Joules, B_tot
  real*8  :: psi_norm, custom_sig, custom_mu
  real*8, parameter :: D_alpha_norm   = 1.d+17
  real*8, parameter :: D_alpha_thresh = 0.0005
  
  !***********************************************************************
  !***********************************************************************
  !********** SECTION: Initialisation & user-input ***********************
  !***********************************************************************
  !***********************************************************************

  ! --- MPI initilisation
  call MPI_INIT(IERR)
  call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)      ! id of each MPI proc
  call MPI_COMM_SIZE(MPI_COMM_WORLD, n_cpu, ierr)      ! number of MPI procs
  

  if (my_id .eq. 0) then
    write(*,*) '******************************************************************************'
    write(*,*) '* fortran_process_jorek.f90                                                  *'
    write(*,*) '* This code requires a JOREK hdf5 restart file                               *'
    write(*,*) '* use "-h" to get options                                                    *'
    write(*,*) '******************************************************************************'
  endif
  
  ! --- Initialise default arguments
  filename       = 'jorek_restart.h5'    ! name of JOREK hdf5 file
  filename_path  = ''                    ! name path to input files
  filename_min   = -1                    ! if using a list, this is the first JOREK file index 
  filename_max   = -1                    ! if using a list, this is the last  JOREK file index
  filename_step  = 1                     ! if using a list, this is the step between JOREK file indices
  filename_list  = .false.               ! if using a list, just set an internal flag...
  filename_data  = 'data_jorek.h5'       ! name of file where data is saved
  save_pixels    = .false.               ! Save grid pixel locations into file to make subsequents runs (much) faster
  use_pixel_file = .false.               ! Read the grid pixel locations from file to make subsequents runs (much) faster (need to run with -save_pixels at least once)
  image_output   = .false.               ! Create an image of the variable you have extracted
  colormap       = 1                     ! Colormap for image (1=heat, 2=rainbow)
  phi_angle      = 0.d0                  ! Toroidal angle where you want your poloidal slice [0,2pi)
  resolution     = 0.01                  ! Poloidal resolution (in meters)
  limit_buffer   = 0.03                  ! The grid will by default be a little larger than the exact JOREK domain
  use_3D_grid    = .false.               ! Save data on a 3D-grid boxed around plasma domain
  use_r_theta    = .false.               ! Save data on (r,theta) coords instead of (R,Z) coords
  i_var(1)       = 5                     ! Variable to extract (default is density rho, ie. JOREK variable 5)
  n_data         = 0                     ! Number of physics quantities to be saved
  i_tor          = -1                    ! Toroidal mode number (-1 for all spectrum)
  hdf5_data      = .true.                ! Data format is HDF5 (by default)
  hdf5_compact   = .false.               ! If extracting multiple time sclices, write a single file (not one file/timestep)
  txt_data       = .false.               ! Data format is txt (default is hdf5)
  bin_data       = .false.               ! Data format is binary ASCII (default is hdf5)
  fourier        = .false.               ! Output data as separated n=0 and n!=0 modes
  D_alpha        = .false.               ! User requests D_alpha
  
  ! --- Initialise variable names
  do i=1,100
    variable_names(i)  = ''
  enddo
  
  ! --- Names of JOREK variables will be fully define once jorek_model is known
  variable_names(1)  = 'psi'
  
  ! --- We start all non-physics-model variables at 20 (to leave space for various JOREK models)
  variable_names(21) = 'BR'
  variable_names(22) = 'BZ'
  variable_names(23) = 'Bp'
  variable_names(24) = 'ER'
  variable_names(25) = 'EZ'
  variable_names(26) = 'Ep'
  variable_names(27) = 'Epar'
  variable_names(28) = 'dpsi_dt'
  variable_names(29) = 'psi_norm'
  variable_names(30) = 'D_alpha' ! requires my_pec.dat !!!
  variable_names(31) = 'custom'


  !***********************************************************************
  !***********************************************************************
  !********** SECTION: Find which JOREK model is used ********************
  !***********************************************************************
  !***********************************************************************

  ! --- Before going through all arguments, we need to know what the JOREK model is
  ! --- Loop through arguments to find the name of the jorek file(s)
  n_arg=command_argument_count()
  do i = 1,n_arg
    call get_command_argument(i,str_input)
    ! --- All other options
    if (trim(str_input) .eq. '-jorek_file') then
      call get_command_argument(i+1,str_input2)
      filename = trim(str_input2)
      if (my_id .eq. 0) write(*,*) 'Using data hdf5 file: ',trim(filename)
    endif
    if (trim(str_input) .eq. '-jorek_file_min') then
      call get_command_argument(i+1,str_input2)
      read(str_input2,'(i5)')filename_min
      if (my_id .eq. 0) write(*,'(A,i0.5,A)') ' Using the first data hdf5 file: jorek',filename_min,'.h5'
    endif
    if (trim(str_input) .eq. '-jorek_file_max') then
      call get_command_argument(i+1,str_input2)
      read(str_input2,'(i5)')filename_max
      if (my_id .eq. 0) write(*,'(A,i0.5,A)') ' Using the last  data hdf5 file: jorek',filename_max,'.h5'
    endif
    if (trim(str_input) .eq. '-jorek_file_step') then
      call get_command_argument(i+1,str_input2)
      read(str_input2,'(i5)')filename_step
      if (my_id .eq. 0) write(*,'(A,i0.5,A)') ' Using increment step for files list: ',filename_step
    endif
    if (trim(str_input) .eq. '-jorek_path') then
      call get_command_argument(i+1,str_input2)
      filename_path = trim(str_input2)
      if (my_id .eq. 0) write(*,*) 'Using specific path for data files: ',trim(filename_path)
    endif
  enddo

  ! --- If using a list of files, make sure min and max are both here.
  if ( (filename_min .ne. -1) .or. (filename_max .ne. -1) ) then
      if ( (my_id .eq. 0) .and. (filename_min .eq. -1) ) then
        write(*,*)'Warning: it seems you specified the last jorek file index, but not the first!'
        write(*,*)'aborting...'
        stop
      endif
      if ( (my_id .eq. 0) .and. (filename_max .eq. -1) ) then
        write(*,*)'Warning: it seems you specified the first jorek file index, but not the last!'
        write(*,*)'aborting...'
        stop
      endif
  endif

  ! --- If using a list of files, set the internal flag
  if ( (filename_min .ne. -1) .and. (filename_max .ne. -1) ) then
    filename_list = .true.
    save_pixels   = .true.
  endif

  ! --- If not using a list of files, we make a list with one entry only
  if (.not. filename_list) then
    filename_min = 1
    filename_max = 1
  endif

  ! --- Initialize FORTRAN interface.
  CALL h5open_f(error)

  ! --- Use the first jorek file
  filename_index = filename_min

  ! --- If doing a list of files, set the filename
  if (filename_list) then
    write(char_min,'(i0.5)')filename_index
    write(filename,'(A,A5,A,A3)')trim(filename_path),'jorek',trim(char_min),'.h5'
  else
    ! --- Include filename path
    if (trim(filename_path) .ne. '') then
      write(filename,'(A,A)')trim(filename_path),trim(filename)
    endif 
  endif

  ! --- Open the JOREK data file.
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  if (my_id .eq. 0) write(*,*)'Opening JOREK data file ',trim(filename)
  call MPI_Barrier(MPI_COMM_WORLD,ierr)

  ! --- Open jorek h5 file.
  call h5fopen_f (trim(filename), H5F_ACC_RDONLY_F, file_id, error)
  if (error .ne. 0) then
    if (my_id .eq. 0) write(*,*)'Failed to open file ',trim(filename)
    if (my_id .eq. 0) write(*,*)'aborting...'
    stop
  endif
  
  ! --- Define names of JOREK variables depending on model
  ! --- IMPORTANT NOTE: the JOREK variable names are not saved in the HDF5 file
  ! ---                 However, it it were, then this extractor program could 
  ! ---                 in principle use exactly the same names, which would be more coherent...
  call HDF5_integer_reading(file_id,jorek_model,'jorek_model')
  call define_jorek_variable_names(jorek_model,variable_names)
  
  ! --- Close the JOREK data file.
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  CALL h5fclose_f(file_id, error)


  !***********************************************************************
  !***********************************************************************
  !********** SECTION: user-inputs ***************************************
  !***********************************************************************
  !***********************************************************************

  ! --- Get all arguments
  do i = 1,n_arg
    call get_command_argument(i,str_input)
    ! --- Asking for help
    if (trim(str_input) .eq. '-h') then
      if (my_id .eq. 0) then
        if (i .lt. n_arg) then
          call get_command_argument(i+1,str_input)
          if ( (trim(str_input) .eq. 'variables') .or. (trim(str_input) .eq. 'variable') ) then
            write(*,*) ''
            write(*,*) 'Available variables are:'
            write(*,*) '  "psi":      poloidal flux (magnetic potential)'
            write(*,*) '  "Phi":      electric potential'
            write(*,*) '  "j":        toroidal current'
            write(*,*) '  "w":        toroidal vorticity'
            write(*,*) '  "rho":      density'
            write(*,*) '  "T":        temperature'
            write(*,*) '  "Ti":       ion temperature (only for model-400)'
            write(*,*) '  "Te":       electron temperature (only for model-400)'
            write(*,*) '  "Vpar":     parallel velocity (only for model >= 300)'
            write(*,*) '  "B":        magnetic field components (BR,BZ,Bp)'
            write(*,*) '  "E":        electric field components (ER,EZ,Ep)'
            write(*,*) '  "BR":       horizontal magnetic field'
            write(*,*) '  "BZ":       vertical magnetic field'
            write(*,*) '  "Bp":       toroidal magnetic field'
            write(*,*) '  "ER":       horizontal electric field'
            write(*,*) '  "EZ":       vertical electric field'
            write(*,*) '  "Ep":       toroidal electric field'
            write(*,*) '  "Epar":     parallel electric field'
            write(*,*) '  "dpsi_dt":  time derivative of psi'
            write(*,*) '  "psi_norm": normalised psi'
            write(*,*) '  "D_alpha":  D_alpha emissions (only for model-500s)'
            write(*,*) '  "custom":   user-defined variable (hard-coded)'
            write(*,*) '  (NOTE: default variable is "rho")'
            write(*,*) '  (NOTE: "psi-norm" jumps to 100.0 outside JOREK domain)'
          else
            write(*,*) ''
            write(*,*) 'Unknown -h flag ',trim(str_input)
            write(*,*) 'Try "-h variables"'
          endif
        else
          write(*,*) ''
          write(*,*) 'Here are the available options:'
          write(*,*) '  -h                     : print this...'
          write(*,*) '  -h variables           : print help about available variables'
          write(*,*) '  -jorek_file <filename> : use the JOREK hdf5 file named <filename> to generate'
          write(*,*) '                           the data (default is "jorek_restart.h5").'
          write(*,*) '                           Instead of using a single file, you can also loop'
          write(*,*) '                           through a list of files using the following options'
          write(*,*) '  -jorek_file_min <filename1> : Integer XXXXX of the first JOREK hdf5 file'
          write(*,*) '                                named "jorekXXXXX.h5" in your list of files.'
          write(*,*) '  -jorek_file_max <filename2> : Integer YYYYY of the last JOREK hdf5 file'
          write(*,*) '                                named "jorekYYYYY.h5" in your list of files.'
          write(*,*) '  -jorek_file_step <interval> : Integer ZZ to increment along the file list'
          write(*,*) '                                like a do loop "do i=XXXXX,YYYYY,Z"'
          write(*,*) '                                default is 1.'
          write(*,*) '                         NOTE : if processing thousands of files with a large'
          write(*,*) '                                3D grid, the pixel_file.txt may take a while'
          write(*,*) '                                to read for every file. Using a list will'
          write(*,*) '                                read the pixel_file.txt only once for all.'
          write(*,*) '  -jorek_path <filepath> : Optional! You can use a path to point to the data'
          write(*,*) '                           files, mostly useful for file lists. For a single'
          write(*,*) '                           you can simply include the path inside'
          write(*,*) '                           the option "-jorek_file <filename>"'
          write(*,*) '                           WARNING: make sure your path and filenames'
          write(*,*) '                           are coherent, eg. path should end with "/"'
          write(*,*) '  -data_file <filename>  : name of text file where data will be saved'
          write(*,*) '                           (default is "data_jorek.txt")'
          write(*,*) '  -txt_data              : Data format is txt (default is hdf5)'
          write(*,*) '  -bin_data              : Data format is binary ASCII (default is hdf5)'
          write(*,*) '                           (only valid for 3D-grids, uses Mitsuba2 format)'
          write(*,*) '  -hdf5_compact          : If extracting multiple time-steps, write only'
          write(*,*) '                           one hdf5 file for all timesteps (default is false)'
          write(*,*) '  -no_data               : Do not extract data to file'
          write(*,*) '  -save_pixels           : save the location of pixels in the JOREK domain to'
          write(*,*) '                           the file "saved_pixels.dat" (default is "false")'
          write(*,*) '  -use_pixel_file        : use the data file named "saved_pixels.dat" to get'
          write(*,*) '                           the location of pixels in the JOREK domain - much'
          write(*,*) '                           much faster - (default is "false") to use this'
          write(*,*) '                           option, you first need to run once with -save_pixels'
          write(*,*) '  -save_image            : save the extracted data into image named'
          write(*,*) '                           "jorek_image.ppm" (default is "false")'
          write(*,*) '  -colorbar <id>         : use id=1 for heat colorbar and id=2 for rainbow'
          write(*,*) '                           colorbar (default is 1)'
          write(*,*) '  -phi <angle>           : select the toroidal angle [0,2pi) at which you want'
          write(*,*) '                           your poloidal slice (default is phi=0)'
          write(*,*) '  -fourier               : instead of a poloidal slice at a given phi angle,'
          write(*,*) '                           including all perturbations, this will output the'
          write(*,*) '                           n=0 and each toroidal mode (sin/cos) perturbation'
          write(*,*) '                           separately, for each variable'
          write(*,*) '  -resolution <res>      : select the resolution of the rectangular grid in'
          write(*,*) '                           meters (default is 0.01m)'
          write(*,*) '  -limit_buffer <buff>   : select the buffer [m] of the rectangular grid with'
          write(*,*) '                           respect to the exact JOREK domainmeters (the default'
          write(*,*) '                           is 0.03m outside the JOREK domain). <buff> can also'
          write(*,*) '                           be negative (ie. using only points that are inside'
          write(*,*) '                           the JOREK domain min/max limits)'
          write(*,*) '  -3D_grid               : Instead of a 2D poloidal grid at a selected phi-angle'
          write(*,*) '                           use a 3D grid (x,y,z) boxed around simulation domain'
          write(*,*) '                           with individual dimensions (nx,ny,nz) each determined'
          write(*,*) '                           by the "resolution" parameter'
          write(*,*) '  -r_theta               : Instead of a 2D (R,Z)-grid, use (r,theta) coords'
          write(*,*) '  -variable <var>        : select the variable you want to extract'
          write(*,*) '                           (note, you can extract many variables at once).'
          write(*,*) '                           To view available variables, run with:'
          write(*,*) '                           "-h variables"'
          write(*,*) '  -mode <n>              : select toroidal mode you want to extract from file'
          write(*,*) '                           n>=0 will select a given mode number'
          write(*,*) '                           n=-1 will give you the whole spectrum'
          write(*,*) '                           default is n=-1'
          write(*,*) '***************************'
          write(*,*) '***************************'
          write(*,*) '***************************'
          write(*,'(A)') 'Here is an example of how to run the code to extract data at 4 different phi planes.'
          write(*,'(A)') 'Note that the first command (to define the pixel points) only needs to be run once'
          write(*,'(A)') 'After this, data can be extracted for any phi-angle, and for any jorek#####.h5 file from the same simulation'
          write(*,'(A)') '>> mpirun -np 16 ./fortran_process -jorek_file jorek01234.h5 -save_pixels -resolution 0.005'
          write(*,'(A)') '>> mpirun -np 16 ./fortran_process -jorek_file jorek01234.h5 -use_pixel_file -phi 0.00000 -variable rho -variable BR -variable BZ -variable Bp -data_file processed_data1.h5'
          write(*,'(A)') '>> mpirun -np 16 ./fortran_process -jorek_file jorek01234.h5 -use_pixel_file -phi 0.31415 -variable rho -variable BR -variable BZ -variable Bp -data_file processed_data2.h5'
          write(*,'(A)') '>> mpirun -np 16 ./fortran_process -jorek_file jorek01234.h5 -use_pixel_file -phi 0.62830 -variable rho -variable BR -variable BZ -variable Bp -data_file processed_data3.h5'
          write(*,'(A)') '>> mpirun -np 16 ./fortran_process -jorek_file jorek01234.h5 -use_pixel_file -phi 0.94245 -variable rho -variable BR -variable BZ -variable Bp -data_file processed_data4.h5'
        endif
      endif
      call MPI_FINALIZE(IERR)
      stop
    endif
    if (trim(str_input) .eq. '-data_file') then
      call get_command_argument(i+1,str_input2)
      filename_data = trim(str_input2)
      if (my_id .eq. 0) write(*,*) 'Writing data to file: ',trim(filename_data)
    endif
    if (trim(str_input) .eq. '-txt_data') then
      hdf5_data = .false.
      txt_data  = .true.
      if (my_id .eq. 0) write(*,*) 'Data format will be .txt'
    endif
    if (trim(str_input) .eq. '-bin_data') then
      hdf5_data = .false.
      bin_data  = .true.
      if (my_id .eq. 0) write(*,*) 'Data format will be .ascii.bin'
    endif
    if (trim(str_input) .eq. '-hdf5_compact') then
      hdf5_compact = .true.
      if (my_id .eq. 0) write(*,*) 'Using compact HDF5 format'
    endif
    if (trim(str_input) .eq. '-no_data') then
      hdf5_data = .false.
      bin_data  = .false.
      txt_data  = .false.
      if (my_id .eq. 0) write(*,*) 'Data will not be saved'
    endif
    if (trim(str_input) .eq. '-save_pixels') then
      save_pixels = .true.
      if (my_id .eq. 0) write(*,*) 'Saving grid pixel locations into "saved_pixels.dat"'
    endif
    if (trim(str_input) .eq. '-use_pixel_file') then
      use_pixel_file = .true.
      if (my_id .eq. 0) write(*,*) 'Reading grid pixel locations from "saved_pixels.dat"'
    endif
    if (trim(str_input) .eq. '-save_image') then
      image_output = .true.
    endif
    if (trim(str_input) .eq. '-colorbar') then
      call get_command_argument(i+1,str_input2)
      read(str_input2,'(i8)')colormap
      if (my_id .eq. 0) then
        if (colormap .eq. 1) write(*,*) 'Using heat colorbar '
        if (colormap .eq. 2) write(*,*) 'Using rainbow colorbar '
      endif
    endif
    if (trim(str_input) .eq. '-phi') then
      call get_command_argument(i+1,str_input2)
      read(str_input2,'(f14.6)')phi_angle
      if (my_id .eq. 0) write(*,*) 'Getting data at toroidal angle : ',phi_angle
    endif
    if (trim(str_input) .eq. '-fourier') then
      fourier = .true.
    endif
    if (trim(str_input) .eq. '-resolution') then
      call get_command_argument(i+1,str_input2)
      read(str_input2,'(f14.6)')resolution
      if (my_id .eq. 0) write(*,*) 'Getting data with resolution : ',resolution
    endif
    if (trim(str_input) .eq. '-limit_buffer') then
      call get_command_argument(i+1,str_input2)
      read(str_input2,'(f14.6)')limit_buffer
      if (my_id .eq. 0) write(*,*) 'Getting data with buffer limit : ',limit_buffer
    endif
    if (trim(str_input) .eq. '-3D_grid') then
      use_3D_grid = .true.
    endif
    if (trim(str_input) .eq. '-r_theta') then
      use_r_theta = .true.
    endif
    if (trim(str_input) .eq. '-variable') then
      n_data = n_data + 1
      call get_command_argument(i+1,str_input2)
      ! --- JOREK variables
      do j=1,20 ! at most 20 JOREK physics model variables
        if (trim(str_input2) .eq. variable_names(j) ) i_var(n_data) = j
      enddo
      ! --- non-physics-model variables
      if (trim(str_input2) .eq. 'BR'   ) i_var(n_data) = 21
      if (trim(str_input2) .eq. 'BZ'   ) i_var(n_data) = 22
      if (trim(str_input2) .eq. 'Bp'   ) i_var(n_data) = 23
      if (trim(str_input2) .eq. 'B'    ) then
        i_var(n_data) = 21
        n_data = n_data + 1
        i_var(n_data) = 22
        n_data = n_data + 1
        i_var(n_data) = 23
      endif
      if (trim(str_input2) .eq. 'ER'   ) i_var(n_data) = 24
      if (trim(str_input2) .eq. 'EZ'   ) i_var(n_data) = 25
      if (trim(str_input2) .eq. 'Ep'   ) i_var(n_data) = 26
      if (trim(str_input2) .eq. 'Epar' ) i_var(n_data) = 27
      if (trim(str_input2) .eq. 'E'    ) then
        i_var(n_data) = 24
        n_data = n_data + 1
        i_var(n_data) = 25
        n_data = n_data + 1
        i_var(n_data) = 26
      endif
      if (trim(str_input2) .eq. 'dpsi_dt' ) i_var(n_data) = 28
      if (trim(str_input2) .eq. 'psi_norm') i_var(n_data) = 29
      if (trim(str_input2) .eq. 'D_alpha' ) i_var(n_data) = 30
      if (trim(str_input2) .eq. 'D_alpha' ) D_alpha       = .true.
      if (trim(str_input2) .eq. 'custom'  ) i_var(n_data) = 31
      if (my_id .eq. 0) write(*,*) 'Getting variable : ',trim(str_input2)
    endif
    if (trim(str_input) .eq. '-mode') then
      call get_command_argument(i+1,str_input2)
      read(str_input2,'(i8)')i_tor
      if (my_id .eq. 0) then
        if (i_tor .eq. -1) write(*,*) 'Extracting complete Fourier spectrum '
        if (i_tor .ge. 0 ) write(*,*) 'Extracting Fourier mode : ',i_tor
      endif
    endif
  enddo
  n_data = max(1,n_data)
  if (txt_data .and. (trim(filename_data) .eq. 'data_jorek.h5') ) filename_data = 'data_jorek.txt'
  if (bin_data .and. (trim(filename_data) .eq. 'data_jorek.h5') ) filename_data = 'data_jorek.ascii.bin'
  if (bin_data .and. (.not. use_3D_grid) ) then
    if (my_id .eq. 0) write(*,*)'Warning, binary data format is reserved for 3D-grids...'
    if (my_id .eq. 0) write(*,*)'Defaulting back to hdf5 format...'
    txt_data  = .false.
    bin_data  = .false.
    hdf5_data = .true.
    filename_data = 'data_jorek.h5'
  endif
  
  ! --- Choose if you want to output an image or a pixel file or use an existing pixel file
  if (save_pixels .and. use_pixel_file) then
    if (my_id .eq. 0) write(*,*)'Warning, you cannot both read and write a pixel file...'
    if (my_id .eq. 0) write(*,*)'Defaulting to reading a pixel file instead of overwriting...'
    save_pixels = .false.
  endif

  ! --- Fourier doesn't make sense on a 3D_grid
  if (use_3D_grid .and. fourier) then
    if (my_id .eq. 0) write(*,*)'Warning, you cannot use Fourier modes with a 3D-grid...'
    if (my_id .eq. 0) write(*,*)'Defaulting to Fourier=False...'
    fourier = .false.
  endif

  ! --- (r,theta)-coords doesn't make sense on a 3D_grid
  if (use_3D_grid .and. use_r_theta) then
    if (my_id .eq. 0) write(*,*)'Warning, you cannot use (r,theta)-coords with a 3D-grid...'
    if (my_id .eq. 0) write(*,*)'Defaulting to use_r_theta=False...'
    use_r_theta = .false.
  endif


  !***********************************************************************
  !***********************************************************************
  !********** SECTION: Read JOREK data from HDF5 *************************
  !***********************************************************************
  !***********************************************************************

  ! --- Loop over each JOREK data file in the list
  do filename_index = filename_min,filename_max,filename_step

    ! --- If doing a list of files, set the filename
    if (filename_list) then
      write(char_min,'(i0.5)')filename_index
      write(filename,'(A,A5,A,A3)')trim(filename_path),'jorek',trim(char_min),'.h5'
    else
      ! --- Include filename path
      if (trim(filename_path) .ne. '') then
        write(filename,'(A,A)')trim(filename_path),trim(filename)
      endif 
    endif

    ! --- Open the JOREK data file.
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    if (my_id .eq. 0) write(*,*)''
    if (my_id .eq. 0) write(*,*)'-------------------------------------------'
    if (my_id .eq. 0) write(*,*)'Opening JOREK data file ',trim(filename)
    call MPI_Barrier(MPI_COMM_WORLD,ierr)

    ! --- Open jorek h5 file.
    call h5fopen_f (trim(filename), H5F_ACC_RDONLY_F, file_id, error)
    if (error .ne. 0) then
      if (my_id .eq. 0) write(*,*)'Failed to open file ',trim(filename)
      if (my_id .eq. 0) write(*,*)'aborting...'
      stop
    endif
    
    ! --- Get toroidal resolution
    call HDF5_integer_reading(file_id, n_tor,    'n_tor')
    call HDF5_integer_reading(file_id, n_period, 'n_period')
    ! --- Convert to JOREK definition of modes
    if (i_tor .ge. 0) then
      i_tor = 2 * i_tor / n_period
    endif
    if (my_id .eq. 0) write(*,*)''
    if (my_id .eq. 0) write(*,*)'Toroidal resolution'
    if (my_id .eq. 0) write(*,*)'number of modes = ',(n_tor-1)/2
    if (my_id .eq. 0) write(*,*)'n_period =',n_period
    if (my_id .eq. 0) write(*,*)'modes included:'
    do i = 0,(n_tor-1)/2
      if (my_id .eq. 0) write(*,'(A,i3)',advance='no')' ',i*n_period
    enddo
    if (my_id .eq. 0) write(*,*)' '
    
    ! --- Define mode numbers
    if (allocated(modes)) deallocate(modes)
    allocate(modes(n_tor))
    do i = 1,n_tor
      modes(i) = int((i)/2) * n_period
    enddo
    
    ! --- Get poloidal grid data
    call HDF5_integer_reading(file_id,n_elements,'n_elements')
    call HDF5_integer_reading(file_id,n_nodes,   'n_nodes')
    call HDF5_integer_reading(file_id,n_var,     'n_var')
    if (allocated(elm_vertex)) deallocate(elm_vertex)
    if (allocated(elm_size)) deallocate(elm_size)
    if (allocated(RZ_nodes)) deallocate(RZ_nodes)
    if (allocated(values)) deallocate(values)
    if (allocated(deltas)) deallocate(deltas)
    allocate(elm_vertex(n_elements,4))
    allocate(elm_size(n_elements,4,4))
    allocate(RZ_nodes(n_nodes,4,2))
    allocate(values(n_nodes,n_tor,4,n_var))
    allocate(deltas(n_nodes,n_tor,4,n_var))
    call HDF5_array2D_reading_int(file_id,elm_vertex, 'vertex')
    call HDF5_array3D_reading    (file_id,elm_size,   'size')
    call HDF5_array3D_reading(file_id,RZ_nodes,       'x')
    call HDF5_array4D_reading(file_id,values,         'values')
    call HDF5_array4D_reading(file_id,deltas,         'deltas')
    
    ! --- Initialise toroidal basis functions (ie. Fourier modes)
    if (allocated(HZ)) deallocate(HZ)
    if (allocated(HZ_p)) deallocate(HZ_p)
    if (allocated(HZ_pp)) deallocate(HZ_pp)
    allocate(HZ(n_tor),HZ_p(n_tor),HZ_pp(n_tor))
    call initialise_toroidal_basis(phi_angle, n_period, n_tor, modes, HZ, HZ_p, HZ_pp)
    if (fourier) then
      n_tor_min = 1
      n_tor_max = n_tor
    else
      n_tor_min = 1
      n_tor_max = 1
    endif
    
    ! --- Get axis, bnd, density and field
    call HDF5_integer_reading(file_id,index_now,'index_now')
    if (index_now .gt. 0) then
      if (allocated(R_axis_t)) deallocate(R_axis_t)
      if (allocated(Z_axis_t)) deallocate(Z_axis_t)
      if (allocated(psi_axis_t)) deallocate(psi_axis_t)
      if (allocated(psi_bnd_t )) deallocate(psi_bnd_t)
      if (allocated(Z_xpoint_t)) deallocate(Z_xpoint_t)
      allocate(R_axis_t(index_now), Z_axis_t(index_now), psi_axis_t(index_now), psi_bnd_t(index_now), Z_xpoint_t(index_now,2))
      call HDF5_array1D_reading(file_id,R_axis_t, 'R_axis_t')
      call HDF5_array1D_reading(file_id,Z_axis_t, 'Z_axis_t')
      call HDF5_array1D_reading(file_id,psi_axis_t, 'psi_axis_t')
      call HDF5_array1D_reading(file_id,psi_bnd_t, 'psi_bnd_t')
      call HDF5_array2D_reading(file_id,Z_xpoint_t, 'Z_xpoint_t')
      R_axis   = R_axis_t(index_now)
      Z_axis   = Z_axis_t(index_now)
      psi_axis = psi_axis_t(index_now)
      psi_bnd  = psi_bnd_t(index_now)
      Z_xpoint(1:2) = Z_xpoint_t(index_now,1:2)
    else
      R_axis   = 10.0
      Z_axis   = 0.0
      psi_axis = 0.d0
      psi_bnd  = 1.d0
      Z_xpoint(1:2) = (/-99.0,+99.0/)
    endif
    if (psi_axis .eq. psi_bnd) then
      psi_axis = 0.d0
      psi_bnd  = 1.d0
    endif
    call HDF5_real_reading(file_id,central_density,'central_density')
    if (central_density .lt. 1.d15) central_density = central_density * 1.d20
    call HDF5_real_reading(file_id,F0,'F0')
    call HDF5_real_reading(file_id,delta_t,'tstep')
    if (my_id .eq. 0) write(*,*)'Finished reading data'
    rho_0     = 3.32d-27 * central_density
    mu_0      = 1.257e-6
    t_norm    = sqrt(rho_0*mu_0)
    eV2Joules = 1.602176487d-19
    
    ! --- Get JOREK time and normalise it
    if (index_now .gt. 0) then
      if (allocated(time)) deallocate(time)
      allocate(time(index_now))
      call HDF5_array1D_reading(file_id,time,           'xtime')
      time_now = time(index_now) * t_norm
    else
      time_now = 0.d0
    endif
  
    ! --- Close the JOREK data file.
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    CALL h5fclose_f(file_id, error)

    ! --- If user requests D_alpha, we need jorek_model500 and my_pec.dat file
    if (D_alpha) then
      ! --- D_alpha makes no sense for harmonics
      if (fourier) then
        if (my_id .eq. 0) then
          write(*,*) 'Variable "D_alpha" only for full toroidal variables.'
          write(*,*) '(ie. do not use "-fourier" option)'
          write(*,*) 'Aborting...'
        endif
        call MPI_FINALIZE(ierr)
        stop
      endif
      ! --- Make sure this is jorek_model500
      if ( (jorek_model .ne. 500) .and. (jorek_model .ne. 545) ) then
        if (my_id .eq. 0) then
          write(*,*) 'You can only request variable "D_alpha" with jorek-model 500.'
          write(*,*) 'Your model:',jorek_model
          write(*,*) 'Aborting...'
        endif
        call MPI_FINALIZE(ierr)
        stop
      endif
      ! --- If user requests D_alpha, we need the my_pec.dat file (only once, not for every file...)
      if (filename_index .eq. filename_min) then
        open(123, file="./my_pec.dat", action='read', iostat=ierr)
        if ( ierr .ne. 0 ) then
          if (my_id .eq. 0) then
            write(*,*) 'Failed to open PEC data file "./my_pec.dat" with error ',ierr
            write(*,*) 'If you request variable "D_alpha", you need the file "my_pec.dat".'
            write(*,*) 'Aborting...'
          endif
          call MPI_FINALIZE(ierr)
          stop
        ! --- Read data file with PEC(Ne,Te) grid (Photon Emissivity Coefficients)  !
        else
          ! --- File should start with a comment line
          read(123,'(A)') line
          ! --- Second line should be "np"
          read(123,'(A)') line
          ! --- Third line should be the value of "np"
          read(123,'(A)') line
          read(line,*) PEC_size
          ! --- Allocate vectors
          allocate(PEC_dens(PEC_size),PEC_temp(PEC_size),PEC(PEC_size*PEC_size))
          ! --- Then comes the density profile
          read(123,'(A)') line
          do i=1, PEC_size
            read(123,'(A)') line
            read(line,*) PEC_dens(i)
          enddo
          ! --- Then comes the temperature profile
          read(123,'(A)') line
          do i=1, PEC_size
            read(123,'(A)') line
            read(line,*) PEC_temp(i)
          enddo
          ! --- Then comes the PEC profiles
          read(123,'(A)') line
          do i=1, PEC_size*PEC_size
            read(123,'(A)') line
            read(line,*) PEC(i)
          enddo
          close(123)
        endif
      endif
    endif
    
    !***********************************************************************
    !***********************************************************************
    !********** SECTION: Setup grid parameters *****************************
    !***********************************************************************
    !***********************************************************************
  
    ! --- Get/Set grid dimensions
    if (filename_index .eq. filename_min) then
      if (use_pixel_file) then
        ! --- Open image file with PPM P3 format
        open(unit=21,file='saved_pixels.dat', ACTION = 'read')
        if (my_id .eq. 0) write(*,*) 'Reading pixel file : saved_pixels.dat'
        ! --- header with pixel size and domain
        if (use_3D_grid) then
          read(21,'(3i8)')  nx_3D, ny_3D, nz_3D
          read(21,'(6e21.11)')   Xpix_min, Xpix_max, Ypix_min, Ypix_max, Zpix_min, Zpix_max
        else
          read(21,'(2i8)')  nx_pix, ny_pix
          read(21,'(4e21.11)')   Rpix_min, Rpix_max, Zpix_min, Zpix_max
        endif
      else
        ! --- Determine boundaries of grid
        if (my_id .eq. 0) write(*,*)'Finding minmax of domain'
        Rpix_min = +1.d10
        Rpix_max = -1.d10
        Zpix_min = +1.d10
        Zpix_max = -1.d10
        do i_elm = 1,n_elements
          call RZ_minmax(n_elements, n_nodes, elm_vertex, elm_size, RZ_nodes, i_elm, Rmin, Rmax, Zmin, Zmax)
          if (Rmin .lt. Rpix_min) Rpix_min = Rmin
          if (Zmin .lt. Zpix_min) Zpix_min = Zmin
          if (Rmax .gt. Rpix_max) Rpix_max = Rmax
          if (Zmax .gt. Zpix_max) Zpix_max = Zmax
        enddo
        ! --- Take few cm extra?
        Rpix_min = Rpix_min - limit_buffer
        Zpix_min = Zpix_min - limit_buffer
        Rpix_max = Rpix_max + limit_buffer
        Zpix_max = Zpix_max + limit_buffer
        ! --- Grid size
        nx_pix = int( (Rpix_max-Rpix_min)/resolution + 1 )
        ny_pix = int( (Zpix_max-Zpix_min)/resolution + 1 )
        ! --- Add Y-dimention if doing 3D
        if (use_3D_grid) then
          Xpix_min =-Rpix_max
          Xpix_max = Rpix_max
          Ypix_min = Xpix_min
          Ypix_max = Xpix_max
          nx_3D = 2*nx_pix
          ny_3D = 2*nx_pix
          nz_3D = ny_pix
        endif
      endif
    endif
    ! --- Allocate data
    if (filename_index .eq. filename_min) then
      if (use_3D_grid) then
        if (my_id .eq. 0) write(*,*)'Doing 3D grid domain'
        if (my_id .eq. 0) write(*,*)'   X-min-max:',Xpix_min,Xpix_max
        if (my_id .eq. 0) write(*,*)'   Y-min-max:',Ypix_min,Ypix_max
        if (my_id .eq. 0) write(*,*)'   Z-min-max:',Zpix_min,Zpix_max
        if (my_id .eq. 0) write(*,*)'   (nX,nY,nZ):',nx_3D,ny_3D,nz_3D
        ! --- We use ny_pix for the parallelisation, with one vector to save indexing
        nx_pix = 1
        ny_pix = nx_3D*ny_3D*nz_3D
        allocate(Xindex3D(ny_pix),Yindex3D(ny_pix),Zindex3D(ny_pix),index3D(nx_3D,ny_3D,nz_3D))
        ny_pix = 0
        do i = 1,nx_3D
          do j = 1,ny_3D
            do k = 1,nz_3D
              ny_pix = ny_pix + 1
              index3D(i,j,k)   = ny_pix
              Xindex3D(ny_pix) = i
              Yindex3D(ny_pix) = j
              Zindex3D(ny_pix) = k
            enddo
          enddo
        enddo
        allocate(Xgrid3D(ny_pix),Ygrid3D(ny_pix),Zgrid3D(ny_pix))
        allocate(my_var(1,ny_pix,1,n_data))
      else
        if (my_id .eq. 0) write(*,*)'Doing rectangular grid domain'
        if (my_id .eq. 0) write(*,*)'   R-min-max:',Rpix_min,Rpix_max
        if (my_id .eq. 0) write(*,*)'   Z-min-max:',Zpix_min,Zpix_max
        if (my_id .eq. 0) write(*,*)'   nR:',nx_pix
        if (my_id .eq. 0) write(*,*)'   nZ:',ny_pix
        allocate(Rgrid(nx_pix,ny_pix),Zgrid(nx_pix,ny_pix))
        if (use_r_theta) allocate(radius_grid(nx_pix,ny_pix),theta_grid(nx_pix,ny_pix),normalised_radius_grid(nx_pix,ny_pix))
        if (fourier) then
          allocate(my_var(nx_pix,ny_pix,n_tor,n_data))
        else
          allocate(my_var(nx_pix,ny_pix,1,n_data))
        endif
      endif
      if ( (save_pixels) .or. (use_pixel_file) ) then
        allocate(i_elm_save(nx_pix,ny_pix), s_save(nx_pix,ny_pix), t_save(nx_pix,ny_pix))
      endif
    endif
    
    ! --- Get the position data from file
    if ( (use_pixel_file) .and. (filename_index .eq. filename_min) ) then
      if (use_3D_grid) then
        ! --- Read pixel file
        if (my_id .eq. 0) write(*,*)'Now reading data from saved_pixels.dat...'
        ny_pix = 0
        do i = 1,nx_3D
          if (my_id .eq. 0) then
            progress = 1.d2 * float(i-1) / float(nx_3D)
            progress = max(0.d0,progress)
            progress = min(1.d2,progress)
            write(*,"(' Processing  ... ',I0,'%',A,$)") int(progress),CHAR(13)
          endif
          do j = 1,ny_3D
            do k = 1,nz_3D
              ny_pix = ny_pix + 1
              read(21,'(4i8,5e21.11)') i_check, j_check, k_check, &
                                       i_elm_save(1,ny_pix), s_save(1,ny_pix), t_save(1,ny_pix), &
                                       Xgrid3D(ny_pix), Ygrid3D(ny_pix), Zgrid3D(ny_pix)
              !too slow to do if-statement here in 3D
              !if ( (i_check .ne. i) .or. (j_check .ne. j) .or. (k_check .ne. k) ) then
              !  if (my_id .eq. 0) write(*,*)'WARNING! Cannot make sense out of file saved_pixels.dat'
              !  if (my_id .eq. 0) write(*,*)'         You should rerun with save_pixel option'
              !  if (my_id .eq. 0) write(*,*)'         Aborting'
              !  call MPI_FINALIZE(IERR)
              !  stop
              !endif
            enddo
          enddo
        enddo
      else
        ! --- Read pixel file
        do i = 1,nx_pix
          do j = 1,ny_pix
            if (use_r_theta) then
              read(21,'(3i8,7e21.11)') i_check, j_check, i_elm_save(i,j), s_save(i,j), t_save(i,j), Rgrid(i,j), Zgrid(i,j), radius_grid(i,j), theta_grid(i,j), normalised_radius_grid(i,j)
            else
              read(21,'(3i8,4e21.11)') i_check, j_check, i_elm_save(i,j), s_save(i,j), t_save(i,j), Rgrid(i,j), Zgrid(i,j)
            endif
            if ( (i_check .ne. i) .or. (j_check .ne. j) ) then
              if (my_id .eq. 0) write(*,*)'WARNING! Cannot make sense out of file saved_pixels.dat'
              if (my_id .eq. 0) write(*,*)'         You should rerun with save_pixel option'
              if (my_id .eq. 0) write(*,*)'         Aborting'
              call MPI_FINALIZE(IERR)
              stop
            endif
          enddo
        enddo
      endif
      close(21)
    endif
      
    !***********************************************************************
    !***********************************************************************
    !********** SECTION: MPI parallelise on ny_pix parameter ***************
    !***********************************************************************
    !***********************************************************************
  
    ! --- The elements we are looking at (note element indexing starts at axis going gradually outwards)
    ! --- We parallelise on ny_pix, both in 2D and 3D (in 3D, ny_pix is the total number of points)
    pix_start = 1
    pix_end   = ny_pix
    
    pix_delta = (pix_end - pix_start) / n_cpu
    local_pix_start = pix_start + my_id*pix_delta + 1
    local_pix_end   = min(pix_end,pix_start+(my_id+1)*pix_delta)
    if (my_id .eq. 0      ) local_pix_start = local_pix_start - 1
    if (my_id .eq. n_cpu-1) local_pix_end   = ny_pix
    
    ! --- Some info print outs
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    if (use_3D_grid) then
      if (my_id .eq. 0) write(*,*) ' Total number of points in 3D grid : ',ny_pix
    else
      if (my_id .eq. 0) write(*,*) ' Total number of horizontal lines : ',ny_pix
    endif
    write(*,'(A,3i8)') ' Local MPI elements (mpi_id, line_start, line_end) : ', my_id, local_pix_start, local_pix_end
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    
    !***********************************************************************
    !***********************************************************************
    !********** SECTION: Interpolate JOREK data on grid ********************
    !***********************************************************************
    !***********************************************************************
  
    ! --- Now get the data on grid
    if (my_id .eq. 0) write(*,*)'Extracting data...'
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    icnt = 0
    do j = local_pix_end,local_pix_start,-1
      ! --- Progress
      if (my_id .eq. 0) then
        progress = 1.d2 - 1.d2 * float(j-local_pix_start) / float(pix_delta)
        progress = max(0.d0,progress)
        progress = min(1.d2,progress)
        write(*,"(' Processing  ... ',I0,'%',A,$)") int(progress),CHAR(13)
      endif
        
      do i = 1,nx_pix
        
        ! --- Get data at RZ-location
        if (use_3D_grid) then
          if ( (use_pixel_file) .or. (filename_index .gt. filename_min) ) then
            i_elm = i_elm_save(1,j)
            s_out = s_save(1,j)
            t_out = t_save(1,j)
            call interp_RZ(elm_vertex, elm_size, RZ_nodes, i_elm, n_elements, n_nodes, s_out, t_out, &
                           R_tmp, RR_s, RR_t, Z_tmp, ZZ_s, ZZ_t)
            ! --- Transform (x,y,z) -> (R,Z,phi)
            R_out = R_tmp
            Z_out = Z_tmp
            x_tmp3D = Xgrid3D(j)
            y_tmp3D = Ygrid3D(j)
            z_tmp3D = Zgrid3D(j)
            phi_angle = atan2(y_tmp3D,x_tmp3D)
            call initialise_toroidal_basis(phi_angle, n_period, n_tor, modes, HZ, HZ_p, HZ_pp)
            ifail = 0
            if (i_elm .eq. 0) ifail = 99
          else
            x_tmp3D = Xpix_min + float(Xindex3D(j)-1)/float(nx_3D-1) * (Xpix_max-Xpix_min)
            y_tmp3D = Ypix_min + float(Yindex3D(j)-1)/float(ny_3D-1) * (Ypix_max-Ypix_min)
            z_tmp3D = Zpix_min + float(Zindex3D(j)-1)/float(nz_3D-1) * (Zpix_max-Zpix_min)
            Xgrid3D(j) = x_tmp3D
            Ygrid3D(j) = y_tmp3D
            Zgrid3D(j) = z_tmp3D
            ! --- Transform (x,y,z) -> (R,Z,phi)
            R_tmp = sqrt( x_tmp3D**2 + y_tmp3D**2 )
            Z_tmp = z_tmp3D
            phi_angle = atan2(y_tmp3D,x_tmp3D)
            call initialise_toroidal_basis(phi_angle, n_period, n_tor, modes, HZ, HZ_p, HZ_pp)
            call find_RZ(n_elements, n_nodes, elm_vertex, elm_size, RZ_nodes, &
                         R_tmp, Z_tmp, R_out, Z_out, i_elm, s_out, t_out, ifail)
            call interp_RZ(elm_vertex, elm_size, RZ_nodes, i_elm, n_elements, n_nodes, s_out, t_out, &
                           R_out, RR_s, RR_t, Z_out, ZZ_s, ZZ_t)
          endif
        else ! not use_3D_grid
          if ( (use_pixel_file) .or. (filename_index .gt. filename_min) ) then
            i_elm = i_elm_save(i,j)
            s_out = s_save(i,j)
            t_out = t_save(i,j)
            call interp_RZ(elm_vertex, elm_size, RZ_nodes, i_elm, n_elements, n_nodes, s_out, t_out, &
                           R_tmp, RR_s, RR_t, Z_tmp, ZZ_s, ZZ_t)
            R_out = R_tmp
            Z_out = Z_tmp
            ifail = 0
            if (i_elm .eq. 0) ifail = 99
          else
            ! --- With (r,theta)-coords, it's a bit more complicated
            if (use_r_theta) then
              ! --- Estimation of maximum radius based on domain size
              radius_max = 2.0 * max(Rpix_max-Rpix_min,Zpix_max-Zpix_min)
              ! --- Define theta
              theta = 2.0 * PI * float(j-1)/float(ny_pix)
              ! --- Only do the search once for each theta
              if (i .eq. 1) then
                call find_radius_with_angle(n_elements, n_nodes, elm_vertex, elm_size, RZ_nodes, R_axis, Z_axis, theta, radius_max, radius_find)
              endif
              ! --- Define radius
              radius = radius_find * float(i-1)/float(nx_pix-1)
              ! --- Get RZ-coords
              R_tmp = R_axis + radius * cos(theta)
              Z_tmp = Z_axis + radius * sin(theta)
              radius_grid(i,j)            = radius
              theta_grid(i,j)             = theta
              normalised_radius_grid(i,j) = radius / radius_find
            else
              R_tmp = Rpix_min + float(i-1)/float(nx_pix-1) * (Rpix_max-Rpix_min)
              Z_tmp = Zpix_min + float(j-1)/float(ny_pix-1) * (Zpix_max-Zpix_min)
            endif
            Rgrid(i,j) = R_tmp
            Zgrid(i,j) = Z_tmp
            call find_RZ(n_elements, n_nodes, elm_vertex, elm_size, RZ_nodes, &
                         R_tmp, Z_tmp, R_out, Z_out, i_elm, s_out, t_out, ifail)
            call interp_RZ(elm_vertex, elm_size, RZ_nodes, i_elm, n_elements, n_nodes, s_out, t_out, &
                           R_out, RR_s, RR_t, Z_out, ZZ_s, ZZ_t)
          endif
        endif
        ! --- Save Pixel file?
        if (save_pixels .or. (filename_list .and. (.not. use_pixel_file) .and. (filename_index .eq. filename_min)) ) then
          if (ifail .eq. 0) then
            i_elm_save(i,j) = i_elm
            s_save(i,j)     = s_out
            t_save(i,j)     = t_out
          else
            i_elm_save(i,j) = 0
            s_save(i,j)     = 0.d0
            t_save(i,j)     = 0.d0
          endif
        endif
        ! --- Now get data at location
        if (ifail .eq. 0) then
          xjac = RR_s*ZZ_t - RR_t*ZZ_s
          do k_tor = n_tor_min, n_tor_max
            if (fourier) then
              call interp_fourier(elm_vertex, elm_size, values, deltas, i_elm, n_elements, n_nodes, &
                                  1, n_var, k_tor, n_tor, s_out, t_out, &
                                  psi, psi_s, psi_t, psi_ss, psi_tt, psi_st, psi_p, delta_psi)
              call interp_fourier(elm_vertex, elm_size, values, deltas, i_elm, n_elements, n_nodes, &
                                  2, n_var, k_tor, n_tor, s_out, t_out, &
                                  phi, phi_s, phi_t, phi_ss, phi_tt, phi_st, phi_p, delta_phi)
            else
              call interp(elm_vertex, elm_size, values, deltas, i_elm, n_elements, n_nodes, &
                          1, n_var, i_tor, n_tor, HZ, HZ_p, s_out, t_out, &
                          psi, psi_s, psi_t, psi_ss, psi_tt, psi_st, psi_p, delta_psi)
              call interp(elm_vertex, elm_size, values, deltas, i_elm, n_elements, n_nodes, &
                          2, n_var, i_tor, n_tor, HZ, HZ_p, s_out, t_out, &
                          phi, phi_s, phi_t, phi_ss, phi_tt, phi_st, phi_p, delta_phi)
            endif
            psi_R = ( + ZZ_t*psi_s - ZZ_s*psi_t ) / xjac
            psi_Z = ( - RR_t*psi_s + RR_s*psi_t ) / xjac
            B_tot = sqrt( F0**2 + psi_R**2 + psi_Z**2) / R_out
            phi_R = ( + ZZ_t*phi_s - ZZ_s*phi_t ) / xjac
            phi_Z = ( - RR_t*phi_s + RR_s*phi_t ) / xjac
            do k = 1,n_data
              if (i_var(k) .le. n_var) then
                if (fourier) then
                  call interp_fourier(elm_vertex, elm_size, values, deltas, i_elm, n_elements, n_nodes, i_var(k), &
                                      n_var, k_tor, n_tor, s_out, t_out, &
                                      pp, pp_s, pp_t, pp_ss, pp_tt, pp_st, pp_p, delta_pp)
                else
                  call interp(elm_vertex, elm_size, values, deltas, i_elm, n_elements, n_nodes, i_var(k), &
                              n_var, i_tor, n_tor, HZ, HZ_p, s_out, t_out, &
                              pp, pp_s, pp_t, pp_ss, pp_tt, pp_st, pp_p, delta_pp)
                endif
                ! --- Normal variable
                my_var(i,j,k_tor,k) = pp
              endif
              if (trim(variable_names(i_var(k))) .eq. 'Phi'     ) my_var(i,j,k_tor,k) = - F0 * my_var(i,j,k_tor,k) / t_norm
              if (trim(variable_names(i_var(k))) .eq. 'j'       ) my_var(i,j,k_tor,k) = - my_var(i,j,k_tor,k) / R_out / mu_0
              if (trim(variable_names(i_var(k))) .eq. 'rho'     ) my_var(i,j,k_tor,k) = max(0.d0,central_density * my_var(i,j,k_tor,k))
              if (trim(variable_names(i_var(k))) .eq. 'T'       ) my_var(i,j,k_tor,k) = max(0.d0,my_var(i,j,k_tor,k) / (central_density * mu_0) / eV2Joules)
              if (trim(variable_names(i_var(k))) .eq. 'Ti'      ) my_var(i,j,k_tor,k) = max(0.d0,my_var(i,j,k_tor,k) / (central_density * mu_0) / eV2Joules)
              if (trim(variable_names(i_var(k))) .eq. 'Te'      ) my_var(i,j,k_tor,k) = max(0.d0,my_var(i,j,k_tor,k) / (central_density * mu_0) / eV2Joules)
              if (trim(variable_names(i_var(k))) .eq. 'rho_n'   ) my_var(i,j,k_tor,k) = max(0.d0,central_density * my_var(i,j,k_tor,k))
              if (trim(variable_names(i_var(k))) .eq. 'Vpar'    ) my_var(i,j,k_tor,k) = my_var(i,j,k_tor,k) * B_tot / t_norm
              if (trim(variable_names(i_var(k))) .eq. 'BR'      ) my_var(i,j,k_tor,k) = psi_Z / R_out
              if (trim(variable_names(i_var(k))) .eq. 'BZ'      ) my_var(i,j,k_tor,k) = -psi_R / R_out
              if (trim(variable_names(i_var(k))) .eq. 'Bp'      ) my_var(i,j,k_tor,k) = F0 / R_out
              if (trim(variable_names(i_var(k))) .eq. 'ER'      ) my_var(i,j,k_tor,k) = - F0 * phi_R / t_norm
              if (trim(variable_names(i_var(k))) .eq. 'EZ'      ) my_var(i,j,k_tor,k) = - F0 * phi_Z / t_norm
              if (trim(variable_names(i_var(k))) .eq. 'Ep'      ) my_var(i,j,k_tor,k) = (- F0 * phi_p / R_out - delta_psi / delta_t / R_out ) / t_norm
              if (trim(variable_names(i_var(k))) .eq. 'dpsi_dt' ) my_var(i,j,k_tor,k) = - delta_psi / delta_t / R_out / t_norm
              if (trim(variable_names(i_var(k))) .eq. 'Epar'    ) then
                BR   = psi_Z / R_out
                BZ   = -psi_R / R_out
                Bp   = F0 / R_out
                ER   = - F0 * phi_R / t_norm
                EZ   = - F0 * phi_Z / t_norm
                Ep   = (- F0 * phi_p / R_out - delta_psi / delta_t / R_out ) / t_norm
                Epar = (BR*ER + BZ*EZ + Bp*Ep) / sqrt(BR**2 + BZ**2 + Bp**2)
                my_var(i,j,k_tor,k) = Epar
              endif
              psi_norm   = (psi - psi_axis) / (psi_bnd - psi_axis)
                if ((psi_norm .lt. 1.d0) .and. (Z_out .lt. Z_xpoint(1)) ) psi_norm = 2.d0 - psi_norm
                if ((psi_norm .lt. 1.d0) .and. (Z_out .gt. Z_xpoint(2)) ) psi_norm = 2.d0 - psi_norm
                if ( psi_norm .lt. 0.d0) psi_norm = 0.d0
              if (trim(variable_names(i_var(k))) .eq. 'psi_norm') my_var(i,j,k_tor,k) = psi_norm
              if (trim(variable_names(i_var(k))) .eq. 'D_alpha'  ) then
                ! --- Need density, neutrals and temperature
                call interp(elm_vertex, elm_size, values, deltas, i_elm, n_elements, n_nodes, &
                            5, n_var, i_tor, n_tor, HZ, HZ_p, s_out, t_out, &
                            rho, pp_s, pp_t, pp_ss, pp_tt, pp_st, pp_p, delta_pp)
                call interp(elm_vertex, elm_size, values, deltas, i_elm, n_elements, n_nodes, &
                            9, n_var, i_tor, n_tor, HZ, HZ_p, s_out, t_out, &
                            rho_n, pp_s, pp_t, pp_ss, pp_tt, pp_st, pp_p, delta_pp)
                call interp(elm_vertex, elm_size, values, deltas, i_elm, n_elements, n_nodes, &
                            6, n_var, i_tor, n_tor, HZ, HZ_p, s_out, t_out, &
                            Ti, pp_s, pp_t, pp_ss, pp_tt, pp_st, pp_p, delta_pp)
                call interp(elm_vertex, elm_size, values, deltas, i_elm, n_elements, n_nodes, &
                            8, n_var, i_tor, n_tor, HZ, HZ_p, s_out, t_out, &
                            Te, pp_s, pp_t, pp_ss, pp_tt, pp_st, pp_p, delta_pp)
                rho   = max(0.d0, rho * central_density)
                rho_n = max(0.d0, rho_n * central_density)
                Te    = max(0.d0, Te/2.0 / (central_density * mu_0 * eV2Joules) )
                ! --- Calculate PEC(Ne,Te)
                do k_pec=2,PEC_size
                  if (rho .lt. PEC_dens(k_pec)) then
                    if (abs(PEC_dens(k_pec)-rho) .lt. abs(PEC_dens(k_pec-1)-rho)) then
                      PEC_index_Ne = k_pec
                    else
                      PEC_index_Ne = k_pec-1        
                    endif
                    exit
                  endif
                  !if ( (k_pec .eq. PEC_size) .and. (my_id .eq. 0) ) &
                  !  write(*,'(A,3e)') 'Warning! no PEC found for density :',rho,PEC_dens(1),PEC_dens(k_pec)
                enddo
                do k_pec=2,PEC_size
                  if (Te .lt. PEC_temp(k_pec)) then
                    if (abs(PEC_temp(k_pec)-Te) .lt. abs(PEC_temp(k_pec-1)-Te)) then
                      PEC_index_Te = k_pec
                    else
                      PEC_index_Te = k_pec-1        
                    endif
                    exit
                  endif
                  !if ( (k_pec .eq. PEC_size) .and. (my_id .eq. 0) ) &
                  !  write(*,'(A,3e)') 'Warning! no PEC found for temperature:',Te,PEC_temp(1),PEC_temp(k_pec)
                enddo
                PEC_index = (PEC_index_Ne-1)*PEC_size + PEC_index_Te
                ! --- Integrate Emissivity
                my_var(i,j,k_tor,k) = rho_n*rho*PEC(PEC_index)
                !my_var(i,j,k_tor,k) = rho * Te**2 * rho_n * (0.5 - 0.5*tanh(-(psi_norm-0.97)/0.001) )
                !my_var(i,j,k_tor,k) = rho**2 * Te**2 * rho_n * (0.5 - 0.5*tanh(-(psi_norm-0.90)/0.03) ) 
                !my_var(i,j,k_tor,k) = rho * Ti * max(0.0005,rho_n/1.d20) * (0.5 - 0.5 * tanh(-(psi_norm-0.95)/0.05) ) * (1.0+abs(Z_out)) * R_out
                !my_var(i,j,k_tor,k) = min(D_alpha_thresh,my_var(i,j,k_tor,k) / D_alpha_norm)
              endif
              if (trim(variable_names(i_var(k))) .eq. 'custom'  ) then
                custom_mu  = 0.97
                custom_sig = 0.05
                my_var(i,j,k_tor,k) = exp(-((psi_norm-custom_mu)/custom_sig)**2.0/2.0) / custom_sig / 2.0/3.1415
              endif
            enddo ! n_data
          enddo ! n_tor
        else
          ! --- If we are outside the simulation domain
          ! --- All variables are set to zero, except psi_norm
          ! --- which is used as a flag, by being set to psi_norm=100.0
          my_var(i,j,:,:) = 0.d0
          do k = 1,n_data
            if (trim(variable_names(i_var(k))) .eq. 'psi_norm') then
              my_var(i,j,n_tor_min:n_tor_max,k) = 1.d2
            endif
          enddo
        endif
        
      enddo
    enddo
    
    !***********************************************************************
    !***********************************************************************
    !********** SECTION: MPI gather data on 0th process ********************
    !***********************************************************************
    !***********************************************************************
  
    ! --- Wait here before sending/collecting data
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    
    ! --- Gather MPI data to process 0
    if (my_id .eq. 0) then
      ! --- If this is mpi_0, we receive data from the other MPIs and print it
      do j=1,n_cpu-1
        call mpi_recv(local_pix_start,1, MPI_INTEGER, j, j, MPI_COMM_WORLD, status, ierr)
        call mpi_recv(local_pix_end,  1, MPI_INTEGER, j, j, MPI_COMM_WORLD, status, ierr)
        nrecv = (local_pix_end-local_pix_start+1)*nx_pix
        if (nrecv .gt. 0) then
          if (save_pixels) then
            call mpi_recv(i_elm_save(:,local_pix_start:local_pix_end),nrecv, MPI_INTEGER,          j, j, MPI_COMM_WORLD, status, ierr)
            call mpi_recv(s_save    (:,local_pix_start:local_pix_end),nrecv, MPI_DOUBLE_PRECISION, j, j, MPI_COMM_WORLD, status, ierr)
            call mpi_recv(t_save    (:,local_pix_start:local_pix_end),nrecv, MPI_DOUBLE_PRECISION, j, j, MPI_COMM_WORLD, status, ierr)
          endif
          if (use_3D_grid) then
            call mpi_recv(Xgrid3D(local_pix_start:local_pix_end),nrecv, MPI_DOUBLE_PRECISION, j, j, MPI_COMM_WORLD, status, ierr)
            call mpi_recv(Ygrid3D(local_pix_start:local_pix_end),nrecv, MPI_DOUBLE_PRECISION, j, j, MPI_COMM_WORLD, status, ierr)
            call mpi_recv(Zgrid3D(local_pix_start:local_pix_end),nrecv, MPI_DOUBLE_PRECISION, j, j, MPI_COMM_WORLD, status, ierr)
          else
            call mpi_recv(Rgrid(:,local_pix_start:local_pix_end),nrecv, MPI_DOUBLE_PRECISION, j, j, MPI_COMM_WORLD, status, ierr)
            call mpi_recv(Zgrid(:,local_pix_start:local_pix_end),nrecv, MPI_DOUBLE_PRECISION, j, j, MPI_COMM_WORLD, status, ierr)
            if (use_r_theta) then
              call mpi_recv(radius_grid           (:,local_pix_start:local_pix_end),nrecv, MPI_DOUBLE_PRECISION, j, j, MPI_COMM_WORLD, status, ierr)
              call mpi_recv(theta_grid            (:,local_pix_start:local_pix_end),nrecv, MPI_DOUBLE_PRECISION, j, j, MPI_COMM_WORLD, status, ierr)
              call mpi_recv(normalised_radius_grid(:,local_pix_start:local_pix_end),nrecv, MPI_DOUBLE_PRECISION, j, j, MPI_COMM_WORLD, status, ierr)
            endif
          endif
          do i=1,n_data
            do k_tor = n_tor_min,n_tor_max
              call mpi_recv(my_var(:,local_pix_start:local_pix_end,k_tor,i),nrecv, MPI_DOUBLE_PRECISION, j, j, MPI_COMM_WORLD, status, ierr)
            enddo
          enddo
        endif
      enddo
    else
      ! --- If this is not mpi_0, we send data to the main MPI 0
      call mpi_send(local_pix_start, 1, MPI_INTEGER, 0, my_id, MPI_COMM_WORLD, ierr)
      call mpi_send(local_pix_end,   1, MPI_INTEGER, 0, my_id, MPI_COMM_WORLD, ierr)
      nsend = (local_pix_end-local_pix_start+1)*nx_pix
      if (nsend .gt. 0) then
        if (save_pixels) then
          call mpi_send(i_elm_save(:,local_pix_start:local_pix_end), nsend, MPI_INTEGER,          0, my_id, MPI_COMM_WORLD, ierr)
          call mpi_send(s_save    (:,local_pix_start:local_pix_end), nsend, MPI_DOUBLE_PRECISION, 0, my_id, MPI_COMM_WORLD, ierr)
          call mpi_send(t_save    (:,local_pix_start:local_pix_end), nsend, MPI_DOUBLE_PRECISION, 0, my_id, MPI_COMM_WORLD, ierr)
        endif
        if (use_3D_grid) then
          call mpi_send(Xgrid3D(local_pix_start:local_pix_end), nsend, MPI_DOUBLE_PRECISION, 0, my_id, MPI_COMM_WORLD, ierr)
          call mpi_send(Ygrid3D(local_pix_start:local_pix_end), nsend, MPI_DOUBLE_PRECISION, 0, my_id, MPI_COMM_WORLD, ierr)
          call mpi_send(Zgrid3D(local_pix_start:local_pix_end), nsend, MPI_DOUBLE_PRECISION, 0, my_id, MPI_COMM_WORLD, ierr)
        else
          call mpi_send(Rgrid(:,local_pix_start:local_pix_end), nsend, MPI_DOUBLE_PRECISION, 0, my_id, MPI_COMM_WORLD, ierr)
          call mpi_send(Zgrid(:,local_pix_start:local_pix_end), nsend, MPI_DOUBLE_PRECISION, 0, my_id, MPI_COMM_WORLD, ierr)
          if (use_r_theta) then
            call mpi_send(radius_grid           (:,local_pix_start:local_pix_end), nsend, MPI_DOUBLE_PRECISION, 0, my_id, MPI_COMM_WORLD, ierr)
            call mpi_send(theta_grid            (:,local_pix_start:local_pix_end), nsend, MPI_DOUBLE_PRECISION, 0, my_id, MPI_COMM_WORLD, ierr)
            call mpi_send(normalised_radius_grid(:,local_pix_start:local_pix_end), nsend, MPI_DOUBLE_PRECISION, 0, my_id, MPI_COMM_WORLD, ierr)
          endif
        endif
        do i=1,n_data
          do k_tor = n_tor_min,n_tor_max
            call mpi_send(my_var(:,local_pix_start:local_pix_end,k_tor,i), nsend, MPI_DOUBLE_PRECISION, 0, my_id, MPI_COMM_WORLD, ierr)
          enddo
        enddo
      endif
    endif
    
    !***********************************************************************
    !***********************************************************************
    !********** SECTION: Write Pixel File **********************************
    !***********************************************************************
    !***********************************************************************
  
    ! --- Write pixel-file
    if (my_id .eq. 0) then
      
      ! --- Write pixel file
      if (save_pixels .and. (filename_index .eq. filename_min)) then
        ! --- Open image file with PPM P3 format
        open(unit=2,file='saved_pixels.dat', ACTION = 'write')
        write(*,*) 'Now writing pixel file : saved_pixels.dat'
        if (use_3D_grid) then
          ! --- header with pixel size and domain
          write(2,'(3i8)')  nx_3D, ny_3D, nz_3D
          write(2,'(6e21.11)')   Xpix_min, Xpix_max, Ypix_min, Ypix_max, Zpix_min, Zpix_max
          ! --- Write to image file
          ny_pix = 0
          do i = 1,nx_3D
            do j = 1,ny_3D
              do k = 1,nz_3D
                ny_pix = ny_pix + 1
                write(2,'(4i8,5e21.11)') i, j, k, &
                                         i_elm_save(1,ny_pix), s_save(1,ny_pix), t_save(1,ny_pix), &
                                         Xgrid3D(ny_pix), Ygrid3D(ny_pix), Zgrid3D(ny_pix)
              enddo
            enddo
          enddo
        else
          ! --- header with pixel size and domain
          write(2,'(2i8)')  nx_pix, ny_pix
          write(2,'(4e21.11)')   Rpix_min, Rpix_max, Zpix_min, Zpix_max
          ! --- Write to image file
          if (use_r_theta) then
            do i = 1,nx_pix
              do j = 1,ny_pix
                write(2,'(3i8,7e21.11)') i, j, i_elm_save(i,j), s_save(i,j), t_save(i,j), Rgrid(i,j), Zgrid(i,j), radius_grid(i,j), theta_grid(i,j), normalised_radius_grid(i,j)
              enddo
            enddo
          else
            do i = 1,nx_pix
              do j = 1,ny_pix
                write(2,'(3i8,4e21.11)') i, j, i_elm_save(i,j), s_save(i,j), t_save(i,j), Rgrid(i,j), Zgrid(i,j)
              enddo
            enddo
          endif
        endif
        close(2)
      endif
  
    endif ! my_id=0
      
    !***********************************************************************
    !***********************************************************************
    !********** SECTION: Write HDF5 Data File ******************************
    !***********************************************************************
    !***********************************************************************
  
    ! --- Write data-file
    if ( (my_id .eq. 0) .and. (hdf5_data) ) then
      
      ! --- Write a compact HDF5 file with all the timesteps
      if (hdf5_compact .and. (filename_min .ne. filename_max) ) then
        ! --- Actions at the first step
        if (filename_index .eq. filename_min) then
          ! --- Allocate data
          n_times = (filename_max-filename_min)/filename_step + 1
          allocate(time_record(n_times))
          if (use_3D_grid) then
            allocate(my_var3D_time(nx_3D,ny_3D,nz_3D,n_data,n_times))
          else
            allocate(my_var_time(nx_pix,ny_pix,n_tor_max-n_tor_min+1,n_data,n_times))
          endif
        endif
        ! --- Record time-dependent data
        i_times = (filename_index - filename_min) / filename_step + 1
        time_record(i_times) = time_now
        if (use_3D_grid) then
          do k = 1,n_data
            do i=1,ny_pix
              my_var3D_time(Xindex3D(i),Yindex3D(j),Zindex3D(i),k,i_times) = my_var(1,i,1,k)
            enddo
          enddo
        else
          if (fourier) then
            do k = 1,n_data
              do i=1,nx_pix
                do j=1,ny_pix
                  do k_tor = n_tor_min, n_tor_max
                    my_var_time(i,j,k_tor,k,i_times) = my_var(i,j,k_tor,k)
                  enddo
                enddo
              enddo
            enddo
          else
            do k = 1,n_data
              do i=1,nx_pix
                do j=1,ny_pix
                  my_var_time(i,j,1,k,i_times) = my_var(i,j,1,k)
                enddo
              enddo
            enddo
          endif
        endif

        ! --- Actions at the last step
        if (filename_index .eq. filename_max) then
          ! --- Open jorek h5 file.
          call HDF5_create(trim(filename_data),data_id,error)
          if (error .ne. 0) then
            if (my_id .eq. 0) write(*,*)'Failed to open file ',trim(filename_data)
            if (my_id .eq. 0) write(*,*)'aborting...'
            stop
          endif
          write(*,*) 'Now writing data file : ', trim(filename_data)
          ! --- Save grid
          if (use_3D_grid) then
            call HDF5_integer_saving(data_id,nx_3D,'nX')
            call HDF5_integer_saving(data_id,ny_3D,'nY')
            call HDF5_integer_saving(data_id,nz_3D,'nZ')
            call HDF5_real_saving   (data_id,Xpix_min,'Xmin')
            call HDF5_real_saving   (data_id,Xpix_max,'Xmax')
            call HDF5_real_saving   (data_id,Ypix_min,'Ymin')
            call HDF5_real_saving   (data_id,Ypix_max,'Ymax')
            call HDF5_real_saving   (data_id,Zpix_min,'Zmin')
            call HDF5_real_saving   (data_id,Zpix_max,'Zmax')
            allocate(my_var3D(nx_3D,ny_3D,nz_3D))
            do i=1,ny_pix
              my_var3D(Xindex3D(i),Yindex3D(j),Zindex3D(i)) = Xgrid3D(i)
            enddo
            call HDF5_array3D_saving(data_id,my_var3D,nx_3D,ny_3D,nz_3D,'Xgrid(nX,nY,nZ)')
            do i=1,ny_pix
              my_var3D(Xindex3D(i),Yindex3D(j),Zindex3D(i)) = Ygrid3D(i)
            enddo
            call HDF5_array3D_saving(data_id,my_var3D,nx_3D,ny_3D,nz_3D,'Ygrid(nX,nY,nZ)')
            do i=1,ny_pix
              my_var3D(Xindex3D(i),Yindex3D(j),Zindex3D(i)) = Zgrid3D(i)
            enddo
            call hdf5_array3d_saving(data_id,my_var3d,nx_3d,ny_3d,nz_3d,'zgrid(nx,ny,nz)')
            deallocate(my_var3d)
          else
            call HDF5_integer_saving(data_id,nx_pix,'nR')
            call HDF5_integer_saving(data_id,ny_pix,'nZ')
            call HDF5_real_saving   (data_id,Rpix_min,'Rmin')
            call HDF5_real_saving   (data_id,Rpix_max,'Rmax')
            call HDF5_real_saving   (data_id,Zpix_min,'Zmin')
            call HDF5_real_saving   (data_id,Zpix_max,'Zmax')
            call HDF5_array2D_saving(data_id,Rgrid,nx_pix,ny_pix,'Rgrid(nR,nZ)')
            call HDF5_array2D_saving(data_id,Zgrid,nx_pix,ny_pix,'Zgrid(nR,nZ)')
            if (use_r_theta) then
              call HDF5_array2D_saving(data_id,radius_grid,           nx_pix,ny_pix,'radius_grid(nR,nZ)')
              call HDF5_array2D_saving(data_id,theta_grid,            nx_pix,ny_pix,'theta_grid(nR,nZ)')
              call HDF5_array2D_saving(data_id,normalised_radius_grid,nx_pix,ny_pix,'normalised_radius_grid(nR,nZ)')
            endif
          endif

          ! --- save recorded time-dependent data
          call HDF5_integer_saving(data_id,n_times,'n_times')
          call HDF5_array1D_saving(data_id,time_record(:), n_times,'time(n_times)')
          if (use_3D_grid) then
            do k = 1,n_data
              write(hdf5_dataname,'(A,A)')trim(variable_names(i_var(k))),'(nx,ny,nz,n_times)'
              call hdf5_array4d_saving(data_id,my_var3D_time(:,:,:,k,:),nx_3d,ny_3d,nz_3d,n_times,hdf5_dataname)
            enddo
          else
            if (fourier) then
              do k = 1,n_data
                write(hdf5_dataname,'(A,A)')trim(variable_names(i_var(k))),'(nR,nZ,n_tor,n_times)'
                call HDF5_array4D_saving(data_id,my_var_time(:,:,:,k,:), nx_pix,ny_pix,n_tor_max-n_tor_min+1,n_times,hdf5_dataname)
              enddo
            else
              do k = 1,n_data
                write(hdf5_dataname,'(A,A)')trim(variable_names(i_var(k))),'(nR,nZ,n_times)'
                call HDF5_array3D_saving(data_id,my_var_time(:,:,1,k,:), nx_pix,ny_pix,n_times,hdf5_dataname)
              enddo
            endif
          endif

          ! --- Deallocate data
          deallocate(time_record)
          if (use_3D_grid) then
            deallocate(my_var3D_time)
          else
            deallocate(my_var_time)
          endif
          ! --- Close HDF5 file
          call HDF5_close(data_id)
        endif

      ! --- Write one HDF5 file per timesteps
      else

        ! --- Determine filename in case of list
        if (filename_list) then
          write(char_min,'(i0.5)')filename_index
          write(filename_data,'(A10,A,A3)')'data_jorek',trim(char_min),'.h5'
        endif
        
        ! --- Open jorek h5 file.
        call HDF5_create(trim(filename_data),data_id,error)
        if (error .ne. 0) then
          if (my_id .eq. 0) write(*,*)'Failed to open file ',trim(filename_data)
          if (my_id .eq. 0) write(*,*)'aborting...'
          stop
        endif
        write(*,*) 'Now writing data file : ', trim(filename_data)
        call HDF5_real_saving   (data_id,time_now,'time_now')
        if (use_3D_grid) then
          call HDF5_integer_saving(data_id,nx_3D,'nX')
          call HDF5_integer_saving(data_id,ny_3D,'nY')
          call HDF5_integer_saving(data_id,nz_3D,'nZ')
          call HDF5_real_saving   (data_id,Xpix_min,'Xmin')
          call HDF5_real_saving   (data_id,Xpix_max,'Xmax')
          call HDF5_real_saving   (data_id,Ypix_min,'Ymin')
          call HDF5_real_saving   (data_id,Ypix_max,'Ymax')
          call HDF5_real_saving   (data_id,Zpix_min,'Zmin')
          call HDF5_real_saving   (data_id,Zpix_max,'Zmax')
          allocate(my_var3D(nx_3D,ny_3D,nz_3D))
          do i=1,ny_pix
            my_var3D(Xindex3D(i),Yindex3D(j),Zindex3D(i)) = Xgrid3D(i)
          enddo
          call HDF5_array3D_saving(data_id,my_var3D,nx_3D,ny_3D,nz_3D,'Xgrid(nX,nY,nZ)')
          do i=1,ny_pix
            my_var3D(Xindex3D(i),Yindex3D(j),Zindex3D(i)) = Ygrid3D(i)
          enddo
          call HDF5_array3D_saving(data_id,my_var3D,nx_3D,ny_3D,nz_3D,'Ygrid(nX,nY,nZ)')
          do i=1,ny_pix
            my_var3D(Xindex3D(i),Yindex3D(j),Zindex3D(i)) = Zgrid3D(i)
          enddo
          call HDF5_array3D_saving(data_id,my_var3D,nx_3D,ny_3D,nz_3D,'Zgrid(nX,nY,nZ)')
          do k = 1,n_data
            do i=1,ny_pix
              my_var3D(Xindex3D(i),Yindex3D(j),Zindex3D(i)) = my_var(1,i,1,k)
            enddo
            call HDF5_array3D_saving(data_id,my_var3D,nx_3D,ny_3D,nz_3D,trim(variable_names(i_var(k))))
          enddo
          deallocate(my_var3D)
        else
          call HDF5_integer_saving(data_id,nx_pix,'nR')
          call HDF5_integer_saving(data_id,ny_pix,'nZ')
          call HDF5_real_saving   (data_id,Rpix_min,'Rmin')
          call HDF5_real_saving   (data_id,Rpix_max,'Rmax')
          call HDF5_real_saving   (data_id,Zpix_min,'Zmin')
          call HDF5_real_saving   (data_id,Zpix_max,'Zmax')
          call HDF5_array2D_saving(data_id,Rgrid,nx_pix,ny_pix,'Rgrid(nR,nZ)')
          call HDF5_array2D_saving(data_id,Zgrid,nx_pix,ny_pix,'Zgrid(nR,nZ)')
          if (use_r_theta) then
            call HDF5_array2D_saving(data_id,radius_grid,           nx_pix,ny_pix,'radius_grid(nR,nZ)')
            call HDF5_array2D_saving(data_id,theta_grid,            nx_pix,ny_pix,'theta_grid(nR,nZ)')
            call HDF5_array2D_saving(data_id,normalised_radius_grid,nx_pix,ny_pix,'normalised_radius_grid(nR,nZ)')
          endif
          if (fourier) then
            do k = 1,n_data
              call HDF5_array3D_saving(data_id,my_var(:,:,:,k), nx_pix,ny_pix,n_tor_max-n_tor_min,trim(variable_names(i_var(k))))
            enddo
          else
            do k = 1,n_data
              call HDF5_array2D_saving(data_id,my_var(:,:,1,k), nx_pix,ny_pix,trim(variable_names(i_var(k))))
            enddo
          endif
        endif
        call HDF5_close(data_id)

      endif ! hdf5_compact

    endif ! my_id=0
      
    !***********************************************************************
    !***********************************************************************
    !********** SECTION: Write TEXT Data File ******************************
    !***********************************************************************
    !***********************************************************************
  
    ! --- Write data-file
    if ( (my_id .eq. 0) .and. (txt_data) ) then
      
      ! --- Determine filename in case of list
      if (filename_list) then
        write(char_min,'(i0.5)')filename_index
        write(filename_data,'(A10,A,A4)')'data_jorek',trim(char_min),'.txt'
      endif

      open(unit=2,file=trim(filename_data), ACTION = 'write')
      write(*,*) 'Now writing data file : ', trim(filename_data)
      if (use_3D_grid) then
        write(2,'(A21)')    'time_now'
        write(2,'(e21.11)')  time_now
        write(2,'(A8)',advance='no')  'nX'
        write(2,'(A8)',advance='no')  'nY'
        write(2,'(A8)')  'nZ'
        write(2,'(3i8)')  nx_3D, ny_3D, nz_3D
        write(2,'(A21)',advance='no')  'Xmin'
        write(2,'(A21)',advance='no')  'Xmax'
        write(2,'(A21)',advance='no')  'Ymin'
        write(2,'(A21)',advance='no')  'Ymax'
        write(2,'(A21)',advance='no')  'Zmin'
        write(2,'(A21)')               'Zmax'
        write(2,'(6e21.11)')   Xpix_min, Xpix_max, Ypix_min, Ypix_max, Zpix_min, Zpix_max
        write(2,'(A21)',advance='no')  'X'
        write(2,'(A21)',advance='no')  'Y'
        write(2,'(A21)',advance='no')  'Z'
        do k = 1,n_data
          write(2,'(A21)',advance='no')  trim(variable_names(i_var(k)))
        enddo
        write(2,'(A)')  ''
        ! --- Write to text file
        ny_pix = 0
        do i = 1,nx_3D
          progress = 1.d2 * float(i-1) / float(nx_3D)
          progress = max(0.d0,progress)
          progress = min(1.d2,progress)
          write(*,"(' Processing  ... ',I0,'%',A,$)") int(progress),CHAR(13)
          do j = 1,ny_3D
            do k = 1,nz_3D
              write(2,'(3e21.11)',advance='no') Xgrid3D(index3D(i,j,k)),Ygrid3D(index3D(i,j,k)),Zgrid3D(index3D(i,j,k))
              ny_pix = ny_pix + 1
              do l = 1,n_data
                  write(2,'(e21.11)',advance='no') my_var(1,ny_pix,1,l)
              enddo
              write(2,'(A)') ''
            enddo
          enddo
        enddo
      else ! not 3D_grid
        write(2,'(A21)')    'time_now'
        write(2,'(e21.11)')  time_now
        write(2,'(A8)',advance='no')  'nR'
        write(2,'(A8)')  'nZ'
        write(2,'(2i8)')  nx_pix, ny_pix
        write(2,'(A21)',advance='no')  'Rmin'
        write(2,'(A21)',advance='no')  'Rmax'
        write(2,'(A21)',advance='no')  'Zmin'
        write(2,'(A21)')               'Zmax'
        write(2,'(4e21.11)')   Rpix_min, Rpix_max, Zpix_min, Zpix_max
        write(2,'(A21)',advance='no')  'R'
        write(2,'(A21)',advance='no')  'Z'
        if (use_r_theta) then
          write(2,'(A21)',advance='no')  'radius'
          write(2,'(A21)',advance='no')  'theta'
          write(2,'(A21)',advance='no')  'normalised_radius'
        endif
        do k = 1,n_data
          do k_tor = n_tor_min, n_tor_max
            if (fourier) then
              if (k_tor .eq. 1) then
                write(char_min,'(A,A)')trim(variable_names(i_var(k))),' (n=0)'
              else
                write(char_min,'(A,A,i2)')trim(variable_names(i_var(k))),' (n=',n_period*((k_tor-mod(k_tor,2)))/2
                if (mod(k_tor,2) .eq. 0) write(char_min,'(A,A)')trim(char_min),' cos)'
                if (mod(k_tor,2) .eq. 1) write(char_min,'(A,A)')trim(char_min),' sin)'
              endif
              write(2,'(A21)',advance='no')  trim(char_min)
            else
              write(2,'(A21)',advance='no')  trim(variable_names(i_var(k)))
            endif
          enddo
        enddo
        write(2,'(A)')  ''
        ! --- Write to text file
        do i = 1,nx_pix
          progress = 1.d2 * float(i-1) / float(ny_pix)
          progress = max(0.d0,progress)
          progress = min(1.d2,progress)
          write(*,"(' Processing  ... ',I0,'%',A,$)") int(progress),CHAR(13)
          do j = 1,ny_pix
            write(2,'(2e21.11)',advance='no') Rgrid(i,j),Zgrid(i,j)
            if (use_r_theta) write(2,'(3e21.11)',advance='no') radius_grid(i,j),theta_grid(i,j),normalised_radius_grid(i,j)
            do k = 1,n_data
              do k_tor = n_tor_min, n_tor_max
                write(2,'(e21.11)',advance='no') my_var(i,j,k_tor,k)
              enddo
            enddo
            write(2,'(A)') ''
          enddo
        enddo
      endif
      write(2,'(A)') ' '
      close(2)
    endif ! my_id=0
      
    !***********************************************************************
    !***********************************************************************
    !********** SECTION: Write ASCII-BINARY File (Mitsuba2 format) *********
    !***********************************************************************
    !***********************************************************************
  
    ! --- Write data-file
    if ( (my_id .eq. 0) .and. (bin_data) ) then
      ! --- Copy data in usable format
      allocate(vtk_scalars3D(n_data,nx_3D,ny_3D,nz_3D))
      do i=1,ny_pix
        do i_data=1,n_data
          vtk_scalars3D(i_data,Xindex3D(i),Yindex3D(i),Zindex3D(i)) = my_var(1,i,1,i_data)
        enddo
      enddo
      ! --- Write one file per data
      do i_data = 1,n_data
        ! --- Determine filename in case of list
        if (filename_list) then
          write(char_min,'(i0.5)')filename_index
          write(filename_data_multiple,'(A10,A,A1,A,A10)')'jorek_data',trim(char_min),'_',trim(variable_names(i_var(i_data))),'.ascii.bin'
        else
          write(filename_data_multiple,'(A11,A,A10)')'jorek_data_',trim(variable_names(i_var(i_data))),'.ascii.bin'
        endif
        ! --- Open file and write
        open(unit=2,file=trim(filename_data_multiple), ACTION = 'write', form="unformatted", status='replace', access='stream')
        write(*,*) 'Now writing data file : ', trim(filename_data_multiple)
        str3 = 'VOL'    ; write(2) str3  ! for Mitsuba2 (guess this means "VOLUME")
        uint8 = 3       ; write(2) uint8 ! for Mitsuba2 (File format version, currently 3)
        int32 = 1       ; write(2) int32 ! for Mitsuba2: [=1->float32] [=2->float16(NotSupported)] [=3->uint8(0..255)] [=4->DenseQuantized]
        int32 = nx_3D   ; write(2) int32
        int32 = ny_3D   ; write(2) int32
        int32 = nz_3D   ; write(2) int32
        int32 = 1       ; write(2) int32
        float32 = Xpix_min ; write(2) float32
        float32 = Ypix_min ; write(2) float32
        float32 = Zpix_min ; write(2) float32
        float32 = Xpix_max ; write(2) float32
        float32 = Ypix_max ; write(2) float32
        float32 = Zpix_max ; write(2) float32
        do k = 1,nz_3D
          progress = 1.d2 * float(k-1) / float(nz_3D)
          progress = max(0.d0,progress)
          progress = min(1.d2,progress)
          write(*,"(' Processing  ... ',I0,'%',A,$)") int(progress),CHAR(13)
          do j = 1,ny_3D
            do i = 1,nx_3D
              write(2) vtk_scalars3D(i_data,i,j,k)
            enddo
          enddo
        enddo
        close(2)
      enddo
      deallocate(vtk_scalars3D)
    endif ! my_id=0
  
    !***********************************************************************
    !***********************************************************************
    !********** SECTION: Write Image File **********************************
    !***********************************************************************
    !***********************************************************************
  
    ! --- Write image-files
    if (my_id .eq. 0) then
      
      ! --- Write image file
      if (image_output) then
  
        ! --- In 3D, we write a .vtk file
        if (use_3D_grid) then
          write(*,*) 'Now writing VTK file : jorek_3D_grid.vtk'
  
          ! --- Determine filename in case of list
          if (filename_list) then
            write(char_min,'(i0.5)')filename_index
            write(filename_data,'(A5,A,A12)')'jorek',trim(char_min),'_3D_data.vtk'
          else
            write(filename_data,'(A17)')'jorek_3D_data.vtk'
          endif

          ! --- Open file
          ivtk = 101 ! just the file write number
#ifdef IBM_MACHINE
          open(unit=ivtk,file=filename_data,form='unformatted',access='stream',status='replace')
#else
          open(unit=ivtk,file=filename_data,form='unformatted',access='stream',convert='BIG_ENDIAN',status='replace')
#endif
          ! --- Headers
          lf = char(10) ! line feed character
          buffer = '# vtk DataFile Version 3.0'//lf                                             ; write(ivtk) trim(buffer)
          buffer = 'vtk output'//lf                                                             ; write(ivtk) trim(buffer)
          buffer = 'BINARY'//lf                                                                 ; write(ivtk) trim(buffer)
          buffer = 'DATASET UNSTRUCTURED_GRID'//lf//lf                                          ; write(ivtk) trim(buffer)
          
          ! --- POINTS SECTION
          write(*,*) '                     : writing points...'
          allocate(vtk_xyz(3,ny_pix), vtk_scalars(n_data,ny_pix))
          do i=1,ny_pix
            vtk_xyz(1,i) = Xgrid3D(i)
            vtk_xyz(2,i) = Ygrid3D(i)
            vtk_xyz(3,i) = Zgrid3D(i)
            do j=1,n_data
              vtk_scalars(j,i) = my_var(1,i,1,j)
            enddo
          enddo
          write(str1(1:10),'(i10)') ny_pix
          buffer = 'POINTS '//str1//'  float'//lf                                               ; write(ivtk) trim(buffer)
          write(ivtk) ((vtk_xyz(i,j),i=1,3),j=1,ny_pix)
          
          ! --- CELLS SECTION
          write(*,*) '                     : writing cells...'
          vtk_points_per_cell = 8
          vtk_n_cells = (nx_3D-1)*(ny_3D-1)*(nz_3D-1)
          allocate(vtk_cells(vtk_points_per_cell,vtk_n_cells))
          vtk_n_cells = 0
          do l=1,ny_pix
            if (Xindex3D(l) .ge. nx_3D-0) cycle
            if (Yindex3D(l) .ge. ny_3D-0) cycle
            if (Zindex3D(l) .ge. nz_3D-0) cycle
            vtk_n_cells = vtk_n_cells + 1
            i = Xindex3D(l) ; j = Yindex3D(l) ; k = Zindex3D(l)
            vtk_cells(1,vtk_n_cells) = index3D(i  ,j  ,k  ) - 1
            vtk_cells(2,vtk_n_cells) = index3D(i+1,j  ,k  ) - 1
            vtk_cells(3,vtk_n_cells) = index3D(i+1,j+1,k  ) - 1
            vtk_cells(4,vtk_n_cells) = index3D(i  ,j+1,k  ) - 1
            vtk_cells(5,vtk_n_cells) = index3D(i  ,j  ,k+1) - 1
            vtk_cells(6,vtk_n_cells) = index3D(i+1,j  ,k+1) - 1
            vtk_cells(7,vtk_n_cells) = index3D(i+1,j+1,k+1) - 1
            vtk_cells(8,vtk_n_cells) = index3D(i  ,j+1,k+1) - 1
          enddo
          write(str1(1:10),'(i10)') vtk_n_cells                          ! number of elements (cells)
          write(str2(1:10),'(i10)') vtk_n_cells*(1+vtk_points_per_cell)  ! size of element list (+1 because first entry is # of entry in line)
          buffer = lf//lf//'CELLS '//str1//' '//str2//lf                                        ; write(ivtk) trim(buffer)
          write(ivtk) (vtk_points_per_cell,(vtk_cells(i,j),i=1,vtk_points_per_cell),j=1,vtk_n_cells)
          
          ! --- CELL_TYPES SECTION
          write(*,*) '                     : writing cell types...'
          etype = 12  ! for vtk_quad
          write(str1(1:10),'(i10)') vtk_n_cells                          ! number of elements (cells)
          buffer = lf//lf//'CELL_TYPES'//str1//lf                                               ; write(ivtk) trim(buffer)
          write(ivtk) (etype,i=1,vtk_n_cells)
          
          ! --- POINT_DATA SECTION
          write(*,*) '                     : writing data on points...'
          write(str1(1:10),'(i10)') ny_pix
          buffer = lf//lf//'POINT_DATA '//str1//lf                                              ; write(ivtk) trim(buffer)
          do i =1, n_data
            buffer = 'SCALARS '//variable_names(i_var(i))//' float'//lf                         ; write(ivtk) trim(buffer)
            buffer = 'LOOKUP_TABLE default'//lf                                                 ; write(ivtk) trim(buffer)
            write(ivtk) (vtk_scalars(i,j),j=1,ny_pix)
          enddo
          
          close(ivtk)
          deallocate(vtk_xyz,vtk_scalars,vtk_cells)
          write(*,*) '                     : finished'
  
        ! --- In 2D, we write a simple pixel .ppm image
        else
          ! --- One file per variable
          do l=1,n_data 
            ! --- Get min/max of data (for normalisation and filename)
            normalise_min = +1.d10
            normalise_max = -1.d10
            do j = 1,ny_pix
              do i = 1,nx_pix
                normalise_min = min(my_var(i,j,1,l),normalise_min)
                normalise_max = max(my_var(i,j,1,l),normalise_max)
              enddo
            enddo
            do j = 1,ny_pix
              do i = 1,nx_pix
                my_var(i,j,1,l) = (my_var(i,j,1,l)-normalise_min) / (normalise_max-normalise_min)
              enddo
            enddo
            ! --- Open image file with PPM P3 format
            ! --- Determine filename in case of list
            if (filename_list) then
              write(filename_data,'(i0.5)')filename_index
              write(char_min,'(sp,e10.3)')normalise_min
              write(char_max,'(sp,e10.3)')normalise_max
              write(filename_picture,'(A5,A,A1,A,A5,A,A5,A,A4)')'jorek',trim(filename_data),'_',trim(variable_names(i_var(l))),'_min_',trim(char_min),'_max_',trim(char_max),'.ppm'
            else
              write(char_min,'(sp,e10.3)')normalise_min
              write(char_max,'(sp,e10.3)')normalise_max
              write(filename_picture,'(A6,A,A5,A,A5,A,A4)')'jorek_',trim(variable_names(i_var(l))),'_min_',trim(char_min),'_max_',trim(char_max),'.ppm'
            endif
            open(unit=2,file=trim(filename_picture),status='unknown')
            write(*,*) 'Now writing PPM (P3) file : ', trim(filename_picture)
            ! --- header
            write(2,'(A)') 'P3'
            write(2,'(2(1x,i4),'' 255 '')')  nx_pix, ny_pix
            ! --- Write to image file
            do j = ny_pix,1,-1
              progress = 1.d2 - 1.d2 * float(j-1) / float(ny_pix)
              progress = max(0.d0,progress)
              progress = min(1.d2,progress)
              write(*,"(' Processing  ... ',I0,'%',A,$)") int(progress),CHAR(13)
              do i = 1,nx_pix
                call get_color(colormap, my_var(i,j,1,l), rgb)
                do k = 1, 3
                  itmp = ichar(rgb(k))
                  icnt = icnt + 4
                  if (icnt .LT. 60) then
                    write(2,fmt='(1x,i3,$)') itmp ! "$" is not standard.
                  else
                    write(2,fmt='(1x,i3)') itmp
                    icnt = 0
                  endif
                enddo
              enddo
            enddo
            write(2,'(A)') ' '
            close(2)
          enddo
        endif ! 3D_grid
      endif ! image_output
     
    endif ! my_id=0
    
    ! --- Wait here before going to next data file
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    
  enddo ! Loop over list of JOREK data files

  !***********************************************************************
  !***********************************************************************
  !********** SECTION: Allocation and MPI cleanup ************************
  !***********************************************************************
  !***********************************************************************
  
  ! --- Cleanup
  if (my_id .eq. 0) write(*,*) 'Finished. Cleaning up...'
  deallocate(modes, elm_vertex, elm_size, RZ_nodes, values, deltas, HZ, HZ_p, HZ_pp)
  deallocate(my_var)
  if (use_3D_grid) then
    deallocate(Xindex3D,Yindex3D,Zindex3D,index3D)
    deallocate(Xgrid3D,Ygrid3D,Zgrid3D)
  else
    deallocate(Rgrid,Zgrid)
    if (use_r_theta) deallocate(radius_grid,theta_grid,normalised_radius_grid) 
  endif
  if ( save_pixels .or. use_pixel_file .or. filename_list ) deallocate(i_elm_save, s_save, t_save)
  
  ! --- Close FORTRAN interface.
  CALL h5close_f(error)

  ! --- Close MPI
  call MPI_FINALIZE(IERR)                                ! clean up MPI
  
END PROGRAM process_hdf5_jorek



















!***********************************************************************
!***********************************************************************
!********** SECTION-ROUTINES: JOREK polynomial routines ****************
!***********************************************************************
!***********************************************************************
  

! --- Given density, return red/green/blue coefs for colormap
subroutine get_color(colormap, density, rgb)

  implicit none
  
  ! --- Routine parameters
  integer, intent(in)      :: colormap
  real*8,  intent(in)      :: density
  character*1, intent(out) :: rgb(3)
  
  ! --- Internal parameters
  real*8  ::  red,  gre,  blu
  integer :: ired, igre, iblu
  

  ! --- Heat colorbar from white to yellow to red to black
  if (colormap .eq. 1) then
    if (density .lt. 1.d0/3.d0) then
      red = 1.d0     
      gre = 1.d0
      blu = 1.d0 - density * 3.d0
    elseif (density .lt. 2.d0/3.d0) then
      red = 1.d0
      gre = 1.d0 - (density-1.d0/3.d0) * 3.d0
      blu = 0.d0
    else
      red = 1.d0 - (density-2.d0/3.d0) * 3.d0
      gre = 0.d0
      blu = 0.d0
    endif
  else
  ! --- Rainbow colorbar from blue to green to white to yellow to red to black
    if (density .lt. 1.d0/5.d0) then
      red = 0.d0   
      gre = density * 5.d0
      blu = 1.d0 - density * 5.d0
    elseif (density .lt. 2.d0/5.d0) then
      red = (density-1.d0/5.d0) * 5.d0
      gre = 1.d0  
      blu = (density-1.d0/5.d0) * 5.d0
    elseif (density .lt. 3.d0/5.d0) then
      red = 1.d0
      gre = 1.d0  
      blu = 1.d0 - (density-2.d0/5.d0) * 5.d0
    elseif (density .lt. 4.d0/5.d0) then
      red = 1.d0
      gre = 1.d0 - (density-3.d0/5.d0) * 5.d0
      blu = 0.d0
    else
      red = 1.d0 - (density-4.d0/5.d0) * 5.d0
      gre = 0.d0
      blu = 0.d0
    endif
  endif
  
  ! --- Convert to 255 bit map  
  ired = int(red * 255.0D+00)
  igre = int(gre * 255.0D+00)
  iblu = int(blu * 255.0D+00)
  if (ired .GT. 255) ired = 255
  if (igre .GT. 255) igre = 255
  if (iblu .GT. 255) iblu = 255
  if (ired .LT.   0) ired =   0
  if (igre .LT.   0) igre =   0
  if (iblu .LT.   0) iblu =   0
  rgb(1) = char(ired)
  rgb(2) = char(igre)
  rgb(3) = char(iblu)
      

  return
end subroutine get_color









! --- Define JOREK variable names depending on JOREK model
subroutine define_jorek_variable_names(jorek_model,variable_names)

  implicit none

  ! --- Routine variables
  integer,      intent(in)    :: jorek_model
  character*50, intent(inout) :: variable_names(100)

  ! --- Internal variables
  integer :: i

  do i=1,20
    variable_names(1)  = ''
  enddo

  if (jorek_model .eq. 002) then
    variable_names(1)  = 'Phi'
    variable_names(2)  = 'w'
    variable_names(3)  = 'rho'
  endif

  if (jorek_model .eq. 003) then
    variable_names(1)  = 'Phi'
    variable_names(2)  = 'w'
    variable_names(3)  = 'rho'
    variable_names(3)  = 'T'
  endif

  if (jorek_model .ge. 199) then
    variable_names(1)  = 'psi'
    variable_names(2)  = 'Phi'
    variable_names(3)  = 'j'
    variable_names(4)  = 'w'
    variable_names(5)  = 'rho'
    variable_names(6)  = 'T'
  endif

  if ( (jorek_model .ge. 300) .and. (jorek_model .lt. 400) ) then
    variable_names(1)  = 'psi'
    variable_names(2)  = 'Phi'
    variable_names(3)  = 'j'
    variable_names(4)  = 'w'
    variable_names(5)  = 'rho'
    variable_names(6)  = 'T'
    variable_names(7)  = 'Vpar'
  endif

  if ( (jorek_model .ge. 400) .and. (jorek_model .lt. 500) ) then
    variable_names(1)  = 'psi'
    variable_names(2)  = 'Phi'
    variable_names(3)  = 'j'
    variable_names(4)  = 'w'
    variable_names(5)  = 'rho'
    variable_names(6)  = 'Ti'
    variable_names(7)  = 'Vpar'
    variable_names(8)  = 'Te'
  endif

  if (jorek_model .eq. 500) then
    variable_names(1)  = 'psi'
    variable_names(2)  = 'Phi'
    variable_names(3)  = 'j'
    variable_names(4)  = 'w'
    variable_names(5)  = 'rho'
    variable_names(6)  = 'T'
    variable_names(7)  = 'Vpar'
    variable_names(8)  = 'rho_n'
  endif

  if (jorek_model .eq. 545) then
    variable_names(1)  = 'psi'
    variable_names(2)  = 'Phi'
    variable_names(3)  = 'j'
    variable_names(4)  = 'w'
    variable_names(5)  = 'rho'
    variable_names(6)  = 'Ti'
    variable_names(7)  = 'Vpar'
    variable_names(8)  = 'Te'
    variable_names(9)  = 'rho_n'
  endif

  return

end subroutine define_jorek_variable_names








! --- Initialise toroidal basis functions (ie. Fourier modes)
subroutine initialise_toroidal_basis(phi_angle, n_period, n_tor, modes, HZ, HZ_p, HZ_pp)

  implicit none
  
  ! --- Routine variables
  integer, intent(in)    :: n_period, n_tor
  real*8,  intent(in)    :: phi_angle, modes(n_tor)
  real*8,  intent(inout) :: HZ(n_tor), HZ_p(n_tor), HZ_pp(n_tor)
  
  ! --- Internal variables
  integer :: i
  
  ! --- Initialise
  HZ	= 0.d0
  HZ_p  = 0.d0
  HZ_pp = 0.d0
  
  HZ   (1) = 1.0
  HZ_p (1) = 0.0
  HZ_pp(1) = 0.0

  do i = 1, (n_tor-1)/2
    HZ   (2*i)   =                           cos(modes(2*i  )*phi_angle)
    HZ_p (2*i)   = - real(modes(2*i  ))    * sin(modes(2*i  )*phi_angle)
    HZ_pp(2*i)   = - real(modes(2*i  ))**2 * cos(modes(2*i  )*phi_angle)
    HZ   (2*i+1) =                           sin(modes(2*i+1)*phi_angle)
    HZ_p (2*i+1) = + real(modes(2*i+1))    * cos(modes(2*i+1)*phi_angle)
    HZ_pp(2*i+1) = - real(modes(2*i+1))**2 * sin(modes(2*i+1)*phi_angle)
  enddo
  
  return

end subroutine










! --- Subroutine calculates the interpolation within one element (i_elm) for a given position (s,t) in the local coordinates
subroutine interp(elm_vertex, elm_size, values, deltas, i_elm, n_elm, n_nodes, i_var, n_var, i_tor, n_tor, HZ, HZ_p, s, t, P, P_s, P_t, P_ss, P_tt, P_st, P_p, delta_P)
  
  implicit none
  
  ! --- Routine parameters
  integer,		    intent(in)  :: elm_vertex(n_elm,4)
  real*8,		    intent(in)  :: elm_size(n_elm,4,4)
  real*8,		    intent(in)  :: values(n_nodes,n_tor,4,n_var)
  real*8,		    intent(in)  :: deltas(n_nodes,n_tor,4,n_var)
  integer,		    intent(in)  :: i_elm
  integer,		    intent(in)  :: n_elm
  integer,		    intent(in)  :: n_nodes
  integer,		    intent(in)  :: i_var
  integer,		    intent(in)  :: n_var
  integer,		    intent(in)  :: i_tor
  integer,		    intent(in)  :: n_tor
  real*8,		    intent(in)  :: HZ(n_tor), HZ_p(n_tor)
  real*8,		    intent(in)  :: s
  real*8,		    intent(in)  :: t
  real*8,		    intent(out) :: P, P_s, P_t, P_ss, P_tt, P_st, P_p, delta_P
  
  ! --- Local variables
  real*8 :: G(4,4), G_s(4,4), G_t(4,4), G_st(4,4), G_ss(4,4), G_tt(4,4)
  integer :: kv, iv, kf, k_tor, n_tor_min, n_tor_max
    
  call basisfunctions(s,t,G(1:4,1:4),G_s(1:4,1:4),G_t(1:4,1:4),G_st(1:4,1:4),G_ss(1:4,1:4),G_tt(1:4,1:4))
  
  P = 0.d0; P_s = 0.d0; P_t = 0.d0; P_st = 0.d0; P_ss = 0.d0; P_tt = 0.d0; P_p = 0.d0; delta_P = 0.d0
  
  do kv = 1,4 ! 4 vertices
  
    iv = elm_vertex(i_elm,kv) ! the node number
  
    do kf = 1,4 ! 4 basis functions
      if (i_tor .ne. -1) then
        n_tor_min = i_tor
        n_tor_max = i_tor+1
      else
        n_tor_min = 1
        n_tor_max = n_tor
      endif
      do k_tor = n_tor_min,n_tor_max
        delta_P = delta_P + deltas(iv,k_tor,kf,i_var) * elm_size(i_elm,kv,kf) * G   (kv,kf) * HZ  (k_tor)
        P       = P       + values(iv,k_tor,kf,i_var) * elm_size(i_elm,kv,kf) * G   (kv,kf) * HZ  (k_tor)
        P_s     = P_s     + values(iv,k_tor,kf,i_var) * elm_size(i_elm,kv,kf) * G_s (kv,kf) * HZ  (k_tor)
        P_t     = P_t     + values(iv,k_tor,kf,i_var) * elm_size(i_elm,kv,kf) * G_t (kv,kf) * HZ  (k_tor)
        P_st    = P_st    + values(iv,k_tor,kf,i_var) * elm_size(i_elm,kv,kf) * G_st(kv,kf) * HZ  (k_tor)
        P_ss    = P_ss    + values(iv,k_tor,kf,i_var) * elm_size(i_elm,kv,kf) * G_ss(kv,kf) * HZ  (k_tor)
        P_tt    = P_tt    + values(iv,k_tor,kf,i_var) * elm_size(i_elm,kv,kf) * G_tt(kv,kf) * HZ  (k_tor)
        P_p     = P_p     + values(iv,k_tor,kf,i_var) * elm_size(i_elm,kv,kf) * G_tt(kv,kf) * HZ_p(k_tor)
      enddo
    enddo

  enddo

  return
end subroutine interp










! --- Subroutine calculates the interpolation within one element (i_elm) for a given position (s,t) in the local coordinates
subroutine interp_fourier(elm_vertex, elm_size, values, deltas, i_elm, n_elm, n_nodes, i_var, n_var, i_tor, n_tor, s, t, P, P_s, P_t, P_ss, P_tt, P_st, P_p, delta_P)
  
  implicit none
  
  ! --- Routine parameters
  integer,		    intent(in)  :: elm_vertex(n_elm,4)
  real*8,		    intent(in)  :: elm_size(n_elm,4,4)
  real*8,		    intent(in)  :: values(n_nodes,n_tor,4,n_var)
  real*8,		    intent(in)  :: deltas(n_nodes,n_tor,4,n_var)
  integer,		    intent(in)  :: i_elm
  integer,		    intent(in)  :: n_elm
  integer,		    intent(in)  :: n_nodes
  integer,		    intent(in)  :: i_var
  integer,		    intent(in)  :: n_var
  integer,		    intent(in)  :: i_tor
  integer,		    intent(in)  :: n_tor
  real*8,		    intent(in)  :: s
  real*8,		    intent(in)  :: t
  real*8,		    intent(out) :: P, P_s, P_t, P_ss, P_tt, P_st, P_p, delta_P
  
  ! --- Local variables
  real*8 :: G(4,4), G_s(4,4), G_t(4,4), G_st(4,4), G_ss(4,4), G_tt(4,4)
  integer :: kv, iv, kf, k_tor, n_tor_min, n_tor_max
    
  call basisfunctions(s,t,G(1:4,1:4),G_s(1:4,1:4),G_t(1:4,1:4),G_st(1:4,1:4),G_ss(1:4,1:4),G_tt(1:4,1:4))
  
  P = 0.d0; P_s = 0.d0; P_t = 0.d0; P_st = 0.d0; P_ss = 0.d0; P_tt = 0.d0; P_p = 0.d0; delta_P = 0.d0
  
  do kv = 1,4 ! 4 vertices
  
    iv = elm_vertex(i_elm,kv) ! the node number
  
    do kf = 1,4 ! 4 basis functions
      k_tor = i_tor
      delta_P = delta_P + deltas(iv,k_tor,kf,i_var) * elm_size(i_elm,kv,kf) * G   (kv,kf)
      P       = P       + values(iv,k_tor,kf,i_var) * elm_size(i_elm,kv,kf) * G   (kv,kf)
      P_s     = P_s     + values(iv,k_tor,kf,i_var) * elm_size(i_elm,kv,kf) * G_s (kv,kf)
      P_t     = P_t     + values(iv,k_tor,kf,i_var) * elm_size(i_elm,kv,kf) * G_t (kv,kf)
      P_st    = P_st    + values(iv,k_tor,kf,i_var) * elm_size(i_elm,kv,kf) * G_st(kv,kf)
      P_ss    = P_ss    + values(iv,k_tor,kf,i_var) * elm_size(i_elm,kv,kf) * G_ss(kv,kf)
      P_tt    = P_tt    + values(iv,k_tor,kf,i_var) * elm_size(i_elm,kv,kf) * G_tt(kv,kf)
      P_p     = P_p     + values(iv,k_tor,kf,i_var) * elm_size(i_elm,kv,kf) * G_tt(kv,kf)
    enddo

  enddo

  return
end subroutine interp_fourier










! --- Subroutine calculates the interpolation within one element (i_elm) for a given position (s,t) in the local coordinates
subroutine interp_RZ(elm_vertex, elm_size, RZ_nodes, i_elm, n_elm, n_nodes, s, t, R, R_s, R_t, Z, Z_s, Z_t)

  implicit none
  
  ! --- Routine parameters
  integer,		    intent(in)  :: elm_vertex(n_elm,4)
  real*8,		    intent(in)  :: elm_size(n_elm,4,4)
  real*8,		    intent(in)  :: RZ_nodes(n_nodes,4,2)
  integer,		    intent(in)  :: i_elm
  integer,		    intent(in)  :: n_elm
  integer,		    intent(in)  :: n_nodes
  real*8,		    intent(in)  :: s,t
  real*8,		    intent(out) :: R, R_s, R_t, Z, Z_s, Z_t
  
  ! --- Local variables
  real*8  :: G(4,4), G_s(4,4), G_t(4,4), G_st(4,4), G_ss(4,4), G_tt(4,4)
  real*8  :: xx1, xx2, ss
  integer :: kv, iv, kf
  
  call basisfunctions(s,t,G,G_s,G_t,G_st,G_ss,G_tt)
  
  R = 0.d0; R_s = 0.d0; R_t = 0.d0
  Z = 0.d0; Z_s = 0.d0; Z_t = 0.d0
  
  do kv = 1,4 ! 4 vertices
  
    iv = elm_vertex(i_elm,kv) ! the node number
  
    do kf = 1,4 ! 4 basis functions
      
      xx1 = RZ_nodes(iv,kf,1)
      xx2 = RZ_nodes(iv,kf,2)
      ss  = elm_size(i_elm,kv,kf)
      
      R    = R    + xx1 * ss * G(kv,kf)
      R_s  = R_s  + xx1 * ss * G_s(kv,kf)
      R_t  = R_t  + xx1 * ss * G_t(kv,kf)
  
      Z    = Z    + xx2 * ss * G(kv,kf)
      Z_s  = Z_s  + xx2 * ss * G_s(kv,kf)
      Z_t  = Z_t  + xx2 * ss * G_t(kv,kf)
  
    enddo
  
  enddo
  
  return
end subroutine interp_RZ
  

  







! --- Find RZ position in the JOREK domain (ie. determine emelent and s,t coords)
subroutine find_RZ(n_elm, n_nodes, elm_vertex, elm_size, RZ_nodes, R_find, Z_find, R_out, Z_out, ielm_out, s_out, t_out, ifail)

  implicit none
  
  ! --- Constants
  integer, parameter :: niter	  = 20  		 !< Maximum number of Newton iterations
  real*8,  parameter :: tolf	  = 1.d-6		 !< Tolerance for spatial distance
  real*8,  parameter :: tolx	  = 1.d-15		 !< Tolerance for iteration step width
  real*8,  parameter :: delta	  = 0.05d0		 !< Maximum number of Newton iterations
  
  ! --- Routine parameters
  integer, intent(in)	 :: n_elm
  integer, intent(in)	 :: n_nodes
  integer, intent(in)	 :: elm_vertex(n_elm,4)
  real*8,  intent(in)	 :: elm_size(n_elm,4,4)
  real*8,  intent(in)	 :: RZ_nodes(n_nodes,4,2)
  real*8,  intent(in)	 :: R_find, Z_find
  real*8,  intent(out)   :: R_out, Z_out, s_out, t_out
  integer, intent(out)   :: ielm_out, ifail
  
  integer :: i, j, k, iv, istart
  real*8  :: Rmin, Rmax, Zmin, Zmax, temp, dis, RR, RR_s, RR_t, ZZ, ZZ_s, ZZ_t, x(2), fvec(2), p(2)
  
  ielm_out = 0
  ifail    = 99
  
  L_EL: do k = 1, n_elm
    
    call RZ_minmax(n_elm, n_nodes, elm_vertex, elm_size, RZ_nodes, k, Rmin, Rmax, Zmin, Zmax)
    
    if ( (R_find > Rmin - delta) .and. (R_find < Rmax + delta) .and. (Z_find > Zmin - delta) .and.   &
       (Z_find < Zmax + delta) ) then ! (If the element could be relevant, proceed:)
      
      L_ST: do istart = 1, 5
  	
  	! Try up to five different starting positions inside the element:
  	if (istart == 1) then
  	  x(:) = (/ 0.50d0, 0.50d0 /)
  	else if (istart == 2) then
  	  x(:) = (/ 0.23d0, 0.23d0 /)
  	else if (istart == 3) then
  	  x(:) = (/ 0.77d0, 0.77d0 /)
  	else if (istart == 4) then
  	  x(:) = (/ 0.77d0, 0.23d0 /)
  	else if (istart == 5) then
  	  x(:) = (/ 0.23d0, 0.77d0 /)
  	end if
  	
  	do i = 1, niter
  	  
          call interp_RZ(elm_vertex, elm_size, RZ_nodes, k, n_elm, n_nodes, x(1), x(2), RR, RR_s, RR_t, ZZ, ZZ_s, ZZ_t)
  	  
  	  fvec(:) = (/ RR - R_find, ZZ - Z_find /)
  	  
  	  if (sqrt(sum(fvec**2)) <= tolf) then
  	    ielm_out = k
  	    exit L_EL
  	  endif
  	  
  	  dis  = ZZ_t * RR_s - RR_t * ZZ_s
  	  if (dis == 0.d0) exit L_ST
  	  
  	  p(:) = (/ RR_t * fvec(2) - ZZ_t * fvec(1), ZZ_s * fvec(1) - RR_s * fvec(2) /) / dis
  	  
  	  p(:) = max( min(p(:),     +0.25d0), -0.25d0 ) ! (limit iteration step size)
  	  x(:) = max( min(x(:)+p(:),+1.00d0), -0.00d0 ) ! (restict s and t to valid range)
  	  
  	  if (sqrt(sum(p**2)) <= tolx) then
  	    ielm_out  = k
  	    exit L_EL
  	  end if
  	  
  	end do
      end do L_ST
      
    end if
    
  end do L_EL
  
  if ( ielm_out /= 0 ) then
    s_out     = x(1)
    t_out     = x(2)
    R_out     = RR
    Z_out     = ZZ
    ifail     = 0
  end if
  
  return
end subroutine find_RZ
  
  
  
  
! --- Find minmax values of R,Z on a given element
subroutine RZ_minmax(n_elm, n_nodes, elm_vertex, elm_size, RZ_nodes, i_elm, Rmin, Rmax, Zmin, Zmax)

  implicit none
  
  ! --- Routine parameters
  integer, intent(in)	 :: n_elm
  integer, intent(in)	 :: n_nodes
  integer, intent(in)	 :: elm_vertex(n_elm,4)
  real*8,  intent(in)	 :: elm_size(n_elm,4,4)
  real*8,  intent(in)	 :: RZ_nodes(n_nodes,4,2)
  integer, intent(in)	 :: i_elm
  real*8,  intent(inout) :: Rmin, Rmax, Zmin, Zmax
  
  ! --- Internal parameters
  real*8  :: psimin, psimax, psma, psmi, psmima, psim, psimr, psip, psipr
  real*8  :: aa, bb, cc, det, r, dummy
  real*8,external :: root
  integer :: iv, n, im, n1, n2
  real*8  :: s,t,P,P_s,P_t,P_st,P_ss,P_tt, ss1, ss2, ss3
  integer :: k
  
  do k=1,2
  
    psimin = 1d10
    psimax =-1d10
  
    do iv= 1, 4
  
      im = mod(iv,4) + 1
      n1 = elm_vertex(i_elm,iv)
      n2 = elm_vertex(i_elm,im)
      
      ss1 = elm_size(i_elm,iv,1)
  
      if ((iv .eq. 1) .or. (iv .eq. 3)) THEN
  
  	ss2 = elm_size(i_elm,iv,2)
  
  	PSIM  =  RZ_nodes(n1,1,k) * ss1      ! PSI(1,n1)
  	PSIMR =  RZ_nodes(n1,2,k) * ss2 * 3.d0/2.d0 ! PSI(3,n1)
  	PSIP  =  RZ_nodes(n2,1,k) * ss1      ! PSI(1,n2)
  	PSIPR =  RZ_nodes(n2,2,k) * ss2 * 3.d0/2.d0 ! PSI(3,n2)
  
      elseif ((iv .eq. 2) .or. (iv .eq. 4)) then
  
  	ss3 = elm_size(i_elm,iv,3)
  	
  	PSIM  =   RZ_nodes(n1,1,k) * ss1	     ! PSI(1,n1)
  	PSIMR =   RZ_nodes(n1,3,k) * ss3 * 3.d0/2.d0 ! PSI(2,n1)
  	PSIP  =   RZ_nodes(n2,1,k) * ss1	     ! PSI(1,n2)
  	PSIPR =   RZ_nodes(n2,3,k) * ss3 * 3.d0/2.d0 ! PSI(2,n2)
  
      endif
  
      if ((PSIM .eq. PSIP) .and. (PSIMR .eq. 0.d0) .and. (PSIPR .eq. 0.d0)) then
  
  	psimin = min(psimin,PSIM)
  	psimax = max(psimax,PSIM)
  
      else
  
  	PSMA = MAX(PSIM,PSIP)
  	PSMI = MIN(PSIM,PSIP)
  	AA =  3.d0 * (PSIM + PSIMR - PSIP + PSIPR ) / 4.d0
  	BB =  ( - PSIMR + PSIPR ) / 2.d0
  	CC =  ( - 3.d0*PSIM - PSIMR + 3.d0*PSIP - PSIPR) / 4.d0
  	DET = BB**2 - 4.d0*AA*CC
  	IF (DET .GE. 0.d0) THEN
  	  R = ROOT(AA,BB,CC,DET,1.d0)
  	  IF (ABS(R) .GT. 1.d0) THEN
  	    R = ROOT(AA,BB,CC,DET,-1.d0)
  	  ENDIF
  	  IF (ABS(R) .LE. 1.d0) THEN
  	    CALL CUB1D(PSIM,PSIMR,PSIP,PSIPR,R,PSMIMA,DUMMY)
  	    psma = max(psma,psmima)
  	    psmi = min(psmi,psmima)
  	  ENDIF
  	ENDIF
  	psimin = min(psimin,psmi)
  	psimax = max(psimax,psma)
  
      endif
  
    enddo
  
    if (k.eq.1) then
      Rmin = psimin
      Rmax = psimax
    else
      Zmin = psimin
      Zmax = psimax
    endif

  ENDDO

  RETURN
END subroutine RZ_minmax
  









subroutine find_radius_with_angle(n_elm, n_nodes, elm_vertex, elm_size, RZ_nodes, R_axis, Z_axis, theta, radius_max, radius_find)

  implicit none
  
  ! --- Routine parameters
  integer, intent(in)    :: n_elm
  integer, intent(in)    :: n_nodes
  integer, intent(in)    :: elm_vertex(n_elm,4)
  real*8,  intent(in)    :: elm_size(n_elm,4,4)
  real*8,  intent(in)    :: RZ_nodes(n_nodes,4,2)
  real*8,  intent(in)    :: R_axis, Z_axis
  real*8,  intent(in)    :: theta, radius_max
  real*8,  intent(inout) :: radius_find 
  
  ! --- Internal parameters
  integer :: iter, n_iter
  real*8  :: distance
  real*8  :: R_prev,  Z_prev
  real*8  :: R_left,  Z_left
  real*8  :: R_right, Z_right
  real*8  :: R_find,  Z_find
  real*8  :: R_out,   Z_out
  real*8  :: s_out,   t_out
  integer :: ielm_out, ifail
  real*8  :: accuracy

  ! --- Brute force search of domain edge
  R_left  = R_axis  
  Z_left  = Z_axis  
  R_right = R_axis + radius_max * cos(theta)
  Z_right = Z_axis + radius_max * sin(theta)
  R_find  = 0.5 * (R_left+R_right)
  Z_find  = 0.5 * (Z_left+Z_right)
  R_prev = R_find
  Z_prev = Z_find

  radius_find = 0.d0
  accuracy = 1.d-6
  n_iter   = 50
  do iter=1,n_iter
    call find_RZ(n_elm, n_nodes, elm_vertex, elm_size, RZ_nodes, R_find, Z_find, R_out, Z_out, ielm_out, s_out, t_out, ifail)
    if (ifail .eq. 0) then
      R_left  = R_out
      Z_left  = Z_out
    else
      R_right = R_find
      Z_right = Z_find
    endif
    R_find  = 0.5 * (R_left+R_right)
    Z_find  = 0.5 * (Z_left+Z_right)
    distance = sqrt( (R_prev-R_find)**2 + (Z_prev-Z_find)**2 )
    if (distance .lt. accuracy) then
      !radius_find = sqrt( (R_find-R_axis)**2 + (Z_find-Z_axis)**2 )
      radius_find = sqrt( (R_left-R_axis)**2 + (Z_left-Z_axis)**2 )
      exit
    endif
    R_prev = R_find
    Z_prev = Z_find
    if ( (iter .eq. n_iter) .and. (radius_find .eq. 0.d0) ) then
      write(*,*)'WARNING: Failed to find radius of domain given theta'
    endif
  enddo

  return 
end subroutine






subroutine cub1D(X1,X1S,X2,X2S,S,X,XS)
  !-----------------------------------------------------------------------
  ! CUBIC HERMITE INTERPOLATION IN ONE DIMENSION
  !-----------------------------------------------------------------------
  implicit none
  real*8   :: X1,X1S,X2,X2S,S,X,XS
  real*8   :: H0M,H0P,H1M,H1P,H0MS,H0PS,H1MS,H1PS

  H0M  =  (S-1.)**2 *(S+2.) * 0.25
  H0MS =  (S-1.)*(S+2.)/2. + (S-1.)**2 * 0.25
  H0P  = -(S+1.)**2 *(S-2.) * 0.25
  H0PS = -(S+1.)*(S-2.)/2. - (S+1.)**2 * 0.25
  H1M  =  (S-1.)**2 *(S+1.) * 0.25
  H1MS =  (S-1.)*(S+1.)/2. + (S-1.)**2 * 0.25
  H1P  =  (S+1.)**2 *(S-1.) * 0.25
  H1PS =  (S+1.)*(S-1.)/2. + (S+1.)**2 * 0.25 

  X  = X1*H0M  + X1S*H1M +  X2*H0P  + X2S*H1P
  XS = X1*H0MS + X1S*H1MS + X2*H0PS + X2S*H1PS

  return
end subroutine cub1D




function root(A,B,C,D,SGN)
  !---------------------------------------------------------------------
  ! THIS FUNCTION GIVES BETTER ROOTS OF QUADRATICS BY AVOIDING
  ! CANCELLATION OF SMALLER ROOT
  !---------------------------------------------------------------------
  implicit none
  real*8 :: root,a, b, c, d, sgn

  if(((B .EQ. 0.D0) .and. (D .EQ. 0.D0)) .or. (A .eq. 0.D0)) then
   root = 1.d20 ! ill defined
   return
  endif
  
  if (B*SGN .GE. 0.d0) then
    root = -2.d0*C/(B+SGN*SQRT(D))
  else
    ROOT = (-B + SGN*SQRT(D)) / (2.d0 * A)
  endif
  return
end function root





! --- Subroutine which defines the basis functions derived from a mixed Bezier/Cubic finite element representation.
subroutine basisfunctions(s, t, H, H_s, H_t, H_st, H_ss, H_tt)

  implicit none
  
  ! --- Routine parameters
  real*8, intent(in)  :: s	    !< s-coordinate in the element
  real*8, intent(in)  :: t	    !< t-coordinate in the element
  real*8, intent(out) :: H(4,4)     !< Basis functions
  real*8, intent(out) :: H_s(4,4)   !< Basis functions derived with respect to s
  real*8, intent(out) :: H_t(4,4)   !< Basis functions derived with respect to t
  real*8, intent(out) :: H_st(4,4)  !< Basis functions derived with respect to s and t
  real*8, intent(out) :: H_ss(4,4)  !< Basis functions derived two times with respect to s
  real*8, intent(out) :: H_tt(4,4)  !< Basis functions derived two times with respect to t
  
  !---------------------------------------------------------- vertex (1)
  H(1,1)   =(-1.d0 + s)**2*(1.d0 + 2.d0*s)*(-1.d0 + t)**2*(1.d0 + 2.d0*t)
  H_s(1,1) =6.d0*(-1.d0 + s)*s*(-1.d0 + t)**2*(1.d0 + 2.d0*t)
  H_t(1,1) =6.d0*(-1.d0 + s)**2*(1.d0 + 2.d0*s)*(-1.d0 + t)*t
  H_st(1,1)=36.d0*(-1.d0 + s)*s*(-1.d0 + t)*t
  H_ss(1,1)=6.d0*(-1.d0 + 2.d0*s)*(-1.d0 + t)**2*(1.d0 + 2.d0*t)
  H_tt(1,1)=6.d0*(-1.d0 + s)**2*(1.d0 + 2.d0*s)*(-1.d0 + 2.d0*t)
  
  H(1,2)   =3.d0*(-1.d0 + s)**2*s*(-1.d0 + t)**2*(1.d0 + 2.d0*t)
  H_s(1,2) =3.d0*(-1.d0 + s)*(-1.d0 + 3.d0*s)*(-1.d0 + t)**2*(1.d0 + 2.d0*t)
  H_t(1,2) =18.d0*(-1.d0 + s)**2*s*(-1.d0 + t)*t
  H_st(1,2)=18.d0*(1.d0 - 4.d0*s + 3.d0*s**2)*(-1.d0 + t)*t
  H_ss(1,2)=6.d0*(-2.d0 + 3.d0*s)*(-1.d0 + t)**2*(1.d0 + 2.d0*t)
  H_tt(1,2)=18.d0*(-1.d0 + s)**2*s*(-1.d0 + 2.d0*t)

  H(1,3)   =3.d0*(-1.d0 + s)**2*(1.d0 + 2.d0*s)*(-1.d0 + t)**2*t
  H_s(1,3) =18.d0*(-1.d0 + s)*s*(-1.d0 + t)**2*t
  H_t(1,3) =3.d0*(-1.d0 + s)**2*(1.d0 + 2.d0*s)*(-1.d0 + t)*(-1.d0 + 3.d0*t)
  H_st(1,3)=18.d0*(-1.d0 + s)*s*(1.d0 - 4.d0*t + 3.d0*t**2)
  H_ss(1,3)=18.d0*(-1.d0 + 2.d0*s)*(-1.d0 + t)**2*t
  H_tt(1,3)=6.d0*(-1.d0 + s)**2*(1.d0 + 2.d0*s)*(-2.d0 + 3.d0*t)
  
  H(1,4)   =9.d0*(-1.d0 + s)**2*s*(-1.d0 + t)**2*t
  H_s(1,4) =9.d0*(-1.d0 + s)*(-1.d0 + 3.d0*s)*(-1.d0 + t)**2*t
  H_t(1,4) =9.d0*(-1.d0 + s)**2*s*(-1.d0 + t)*(-1.d0 + 3.d0*t)
  H_st(1,4)=9.d0*(1.d0 - 4.d0*s + 3.d0*s**2)*(1.d0 - 4.d0*t + 3.d0*t**2)
  H_ss(1,4)=18.d0*(-2.d0 + 3.d0*s)*(-1.d0 + t)**2*t
  H_tt(1,4)=18.d0*(-1.d0 + s)**2*s*(-2.d0 + 3.d0*t)
  
  !---------------------------------------------------------- vertex (2)
  H(2,1)   =-(s**2*(-3.d0 + 2.d0*s)*(-1.d0 + t)**2*(1.d0 + 2.d0*t))
  H_s(2,1) =-6.d0*(-1.d0 + s)*s*(-1.d0 + t)**2*(1.d0 + 2.d0*t)
  H_t(2,1) =-6.d0*s**2*(-3.d0 + 2.d0*s)*(-1.d0 + t)*t
  H_st(2,1)=-36.d0*(-1.d0 + s)*s*(-1.d0 + t)*t
  H_ss(2,1)=-6.d0*(-1.d0 + 2.d0*s)*(-1.d0 + t)**2*(1.d0 + 2.d0*t)
  H_tt(2,1)=-6.d0*s**2*(-3.d0 + 2.d0*s)*(-1.d0 + 2.d0*t)
  
  H(2,2)   =-3.d0*(-1.d0 + s)*s**2*(-1.d0 + t)**2*(1.d0 + 2.d0*t)
  H_s(2,2) =-3.d0*s*(-2.d0 + 3.d0*s)*(-1.d0 + t)**2*(1.d0 + 2.d0*t)
  H_t(2,2) =-18.d0*(-1.d0 + s)*s**2*(-1.d0 + t)*t
  H_st(2,2)=-18.d0*s*(-2.d0 + 3.d0*s)*(-1.d0 + t)*t
  H_ss(2,2)=-6.d0*(-1.d0 + 3.d0*s)*(-1.d0 + t)**2*(1 + 2.d0*t)
  H_tt(2,2)=-18.d0*(-1.d0 + s)*s**2*(-1.d0 + 2.d0*t)
  
  H(2,3)   =-3.d0*s**2*(-3.d0 + 2.d0*s)*(-1.d0 + t)**2*t
  H_s(2,3) =-18.d0*(-1.d0 + s)*s*(-1.d0 + t)**2*t
  H_t(2,3) =3.d0*s**2*(-3.d0 + 2.d0*s)*(1.d0 - 3.d0*t)*(-1.d0 + t)
  H_st(2,3)=-18.d0*(-1.d0 + s)*s*(1.d0 - 4.d0*t + 3.d0*t**2)
  H_ss(2,3)=6.d0*(3.d0 - 6.d0*s)*(-1.d0 + t)**2*t
  H_tt(2,3)=-6.d0*s**2*(-3.d0 + 2.d0*s)*(-2.d0 + 3.d0*t)
  
  H(2,4)   =-9.d0*(-1.d0 + s)*s**2*(-1.d0 + t)**2*t
  H_s(2,4) =-9.d0*s*(-2.d0 + 3.d0*s)*(-1.d0 + t)**2*t
  H_t(2,4) =9.d0*(-1.d0 + s)*s**2*(1.d0 - 3.d0*t)*(-1.d0 + t)
  H_st(2,4)=-9.d0*s*(-2.d0 + 3.d0*s)*(1.d0 - 4.d0*t + 3.d0*t**2)
  H_ss(2,4)=18.d0*(1.d0 - 3.d0*s)*(-1.d0 + t)**2*t
  H_tt(2,4)=-18.d0*(-1.d0 + s)*s**2*(-2.d0 + 3.d0*t)
  
  !---------------------------------------------------------- vertex (3)
  H(3,1)   =s**2*(-3.d0 + 2.d0*s)*t**2*(-3.d0 + 2.d0*t)
  H_s(3,1) =6.d0*(-1.d0 + s)*s*t**2*(-3.d0 + 2.d0*t)
  H_t(3,1) =6.d0*s**2*(-3.d0 + 2.d0*s)*(-1.d0 + t)*t
  H_st(3,1)=36.d0*(-1.d0 + s)*s*(-1.d0 + t)*t
  H_ss(3,1)=6.d0*(-1.d0 + 2*s)*t**2*(-3.d0 + 2.d0*t)
  H_tt(3,1)=6.d0*s**2*(-3.d0 + 2.d0*s)*(-1.d0 + 2.d0*t)
  
  H(3,2)   =3.d0*(-1.d0 + s)*s**2*t**2*(-3.d0 + 2.d0*t)
  H_s(3,2) =3.d0*s*(-2 + 3.d0*s)*t**2*(-3.d0 + 2.d0*t)
  H_t(3,2) =18.d0*(-1.d0 + s)*s**2*(-1.d0 + t)*t
  H_st(3,2)=18.d0*s*(-2.d0 + 3.d0*s)*(-1.d0 + t)*t
  H_ss(3,2)=6.d0*(-1.d0 + 3.d0*s)*t**2*(-3.d0 + 2.d0*t)
  H_tt(3,2)=18.d0*(-1.d0 + s)*s**2*(-1.d0 + 2.d0*t)
  
  H(3,3)   =3.d0*s**2*(-3.d0 + 2.d0*s)*(-1.d0 + t)*t**2
  H_s(3,3) =18.d0*(-1.d0 + s)*s*(-1.d0 + t)*t**2
  H_t(3,3) =3.d0*s**2*(-3.d0 + 2.d0*s)*t*(-2.d0 + 3.d0*t)
  H_st(3,3)=18.d0*(-1.d0 + s)*s*t*(-2.d0 + 3.d0*t)
  H_ss(3,3)=18.d0*(-1.d0 + 2.d0*s)*(-1.d0 + t)*t**2
  H_tt(3,3)=6.d0*s**2*(-3.d0 + 2.d0*s)*(-1.d0 + 3.d0*t)
  
  H(3,4)   =9.d0*(-1.d0 + s)*s**2*(-1.d0 + t)*t**2
  H_s(3,4) =9.d0*s*(-2.d0 + 3.d0*s)*(-1.d0 + t)*t**2
  H_t(3,4) =9.d0*(-1.d0 + s)*s**2*t*(-2.d0 + 3.d0*t)
  H_st(3,4)=9.d0*s*(-2.d0 + 3.d0*s)*t*(-2.d0 + 3.d0*t)
  H_ss(3,4)=18.d0*(-1.d0 + 3.d0*s)*(-1.d0 + t)*t**2
  H_tt(3,4)=18.d0*(-1.d0 + s)*s**2*(-1.d0 + 3.d0*t)
  
  !---------------------------------------------------------- vertex (4)
  H(4,1)   =-((-1.d0 + s)**2*(1.d0 + 2.d0*s)*t**2*(-3.d0 + 2.d0*t))
  H_s(4,1) =-6.d0*(-1.d0 + s)*s*t**2*(-3.d0 + 2.d0*t)
  H_t(4,1) =-6.d0*(-1.d0 + s)**2*(1.d0 + 2.d0*s)*(-1.d0 + t)*t
  H_st(4,1)=-36.d0*(-1.d0 + s)*s*(-1.d0 + t)*t
  H_ss(4,1)=-6.d0*(-1.d0 + 2.d0*s)*t**2*(-3.d0 + 2.d0*t)
  H_tt(4,1)=-6.d0*(-1.d0 + s)**2*(1.d0 + 2.d0*s)*(-1.d0 + 2.d0*t)
  
  H(4,2)   =-3.d0*(-1.d0 + s)**2*s*t**2*(-3.d0 + 2.d0*t)
  H_s(4,2) =3.d0*(1.d0 - 3*s)*(-1.d0 + s)*t**2*(-3.d0 + 2.d0*t)
  H_t(4,2) =-18.d0*(-1.d0 + s)**2*s*(-1.d0 + t)*t
  H_st(4,2)=-18.d0*(1.d0 - 4.d0*s + 3.d0*s**2)*(-1.d0 + t)*t
  H_ss(4,2)=-6.d0*(-2.d0 + 3.d0*s)*t**2*(-3.d0 + 2.d0*t)
  H_tt(4,2)=6.d0*(-1.d0 + s)**2*s*(3.d0 - 6.d0*t)
  
  H(4,3)   =-3.d0*(-1.d0 + s)**2*(1.d0 + 2.d0*s)*(-1.d0 + t)*t**2
  H_s(4,3) =-18.d0*(-1.d0 + s)*s*(-1.d0 + t)*t**2
  H_t(4,3) =-3.d0*(-1.d0 + s)**2*(1.d0 + 2.d0*s)*t*(-2.d0 + 3.d0*t)
  H_st(4,3)=-18.d0*(-1.d0 + s)*s*t*(-2.d0 + 3.d0*t)
  H_ss(4,3)=-18.d0*(-1.d0 + 2.d0*s)*(-1.d0 + t)*t**2
  H_tt(4,3)=-6.d0*(-1.d0 + s)**2*(1.d0 + 2.d0*s)*(-1 + 3.d0*t)

  H(4,4)   =-9.d0*(-1.d0 + s)**2*s*(-1.d0 + t)*t**2
  H_s(4,4) =9.d0*(1.d0 - 3.d0*s)*(-1.d0 + s)*(-1.d0 + t)*t**2
  H_t(4,4) =-9.d0*(-1.d0 + s)**2*s*t*(-2.d0 + 3.d0*t)
  H_st(4,4)=-9.d0*(1.d0 - 4.d0*s + 3.d0*s**2)*t*(-2.d0 + 3.d0*t)
  H_ss(4,4)=-18.d0*(-2.d0 + 3.d0*s)*(-1.d0 + t)*t**2
  H_tt(4,4)=18.d0*(-1.d0 + s)**2*s*(1.d0 - 3.d0*t)
  
  return
end subroutine basisfunctions
  
