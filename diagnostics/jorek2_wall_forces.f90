!< This program reads jorek_restart files and a "starwall_response.dat" file
!! and calculates the total wall forces on the STARWALL wall. It involves the 
!! computation of expensive volume integrals in the plasma, so run it
!! in paralllel with several MPI processes for large cases. 
!! To select the given JOREK restart files create a file named wall_forces.nml
!! with 3 integers in the first line (istart, iend, delta_step)
!!   istart     = 1st  restart file index
!!   iend       = Last restart file index
!!   delta_step = index jump between the restart files
!! The forces are exported to the total_wall_forces.dat file
!! Note that to calculate the require plasma fields n_plane typically should be
!! more than 60 (even for 2D)
!! See more documentation in
!!    https://www.jorek.eu/wiki/doku.php?id=jorek2_wall_forces
program jorek2_wall_forces 

  use mod_vacuum_fields,   only: total_wall_forces 
  use data_structure
  use phys_module
  use mod_parameters
  use nodes_elements
  use mod_boundary,            only: boundary_from_grid
  use vacuum
  use vacuum_response,     only: get_vacuum_response, broadcast_starwall_response, init_wall_currents
  use vacuum_equilibrium,  only: import_external_fields
  use mod_import_restart
  use mod_element_rtree, only: populate_element_rtree
  use basis_at_gaussian, only: initialise_basis
  use tr_module

#ifdef USE_HDF5
  use hdf5
  use hdf5_io_module
  use matio_module, only: timestamp
#endif
  use mpi_mod

  use, intrinsic :: iso_c_binding
  use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
                                            stdout=>output_unit, &
                                            stderr=>error_unit
  implicit none
  character(len=60), parameter ::      FMT = "(1I5.5,'..', 1I5.5,'_' , 1f5.3 )" 
  integer   :: my_id, my_id_n, my_id_master, ierr, ierr2
  integer   :: i_rank(n_tor), n_cpu, n_cpu_n, n_cpu_master, m_cpu, n_masters, n_cpu_trans, my_id_trans
  integer   :: MPI_COMM_N, MPI_GROUP_MASTER, MPI_GROUP_WORLD, MPI_COMM_MASTER, MPI_COMM_TRANS
  integer   :: required,provided,StatInfo
  integer   :: istep, delta_step, istart, iend
  integer*4 :: rank, comm_size 
  logical   :: first_step

  real*8    :: Fx, Fy, Fz, scale_fact=1.d5

  character*17                          :: file_in
  character*20                          :: fact
  character(len=MPI_MAX_PROCESSOR_NAME) :: name
  integer :: resultlength
  namelist /wall_forces/ istart, iend, delta_step, scale_fact 
 
  !***********************************************************************
  !*                  intialisation (copied from jorek2_main)            *
  !***********************************************************************
#ifdef FUNNELED
  required = MPI_THREAD_FUNNELED
#else
  required = MPI_THREAD_MULTIPLE
#endif

  call MPI_Init_thread(required, provided, StatInfo)

  call init_threads()  ! on some systems init_threads needs to come after mpi_init_thread
  
  ! --- Determine number of MPI procs
  call MPI_COMM_SIZE(MPI_COMM_WORLD, comm_size, ierr)
  n_cpu = comm_size
  
  ! --- Determine ID of each MPI proc
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  my_id = rank
  
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  CALL MPI_GET_PROCESSOR_NAME (name,resultlength,ierr)
  write(*,'(A,I5,2A)') '  #MPI id, ProcessorName ', rank, ': ', name
  call MPI_Barrier(MPI_COMM_WORLD,ierr)

  ! --- Initialize mode and mode_type arrays
  call det_modes()
 
  ! --- Preset input parameters to reasonable defaults, then read the input file.
  call initialise_and_broadcast_parameters(my_id, "__NO_FILENAME__")
  
  ! --- Initialize the vacuum part.
  call vacuum_init(my_id, freeboundary_equil, freeboundary, resistive_wall)
  
  call MPI_Barrier(MPI_COMM_WORLD,ierr)

  ! --- Define the basis functions at the Gaussian points
  call initialise_basis()

  first_step = .true.
  if (my_id==0) then
     open(25, file='wall_forces.nml', action='read', status='old', iostat=ierr)
     if (ierr==0) then
        read(25,wall_forces)
     end if
     if (scale_fact > 1.d3 .or. scale_fact < 1.d-3) scale_fact = 1.01d0
     write(*,*) 'scale factor:', scale_fact
     write(fact,FMT) istart, iend, scale_fact
     open(87,file=trim('total_wall_forces_')//trim(fact)//'.dat',action='write')
     write(87,*) '#Step  time(norm)   time(ms)      Fx(N)        Fy(N)          Fz(N)'
  end if

  call MPI_BCAST(     istart,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(       iend,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST( delta_step,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST( scale_fact,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)

  ! --- Loop over restart files
  do istep = istart, iend, delta_step 

    write(file_in,'(A5,i5.5)') 'jorek', istep

    if ( my_id == 0 ) then
      call import_restart(node_list, element_list, file_in, rst_format, ierr2)
    endif
    call MPI_BCAST(ierr2,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    if ( ierr2 /= 0 ) cycle

    call broadcast_phys(my_id)  
    call broadcast_elements(my_id, element_list)                ! elements
    call broadcast_nodes(my_id, node_list)                      ! nodes
  
    if (.not. freeboundary) then
      write(*,*) ' **** Fatal: jorek2_wall_forces needs freeboundary simulations ****'
      stop
    endif

    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    ! --- Fill the vacuum response matrices for freeboundary computations
    if (first_step) then
    
       ! --- Determine boundary information from the grid
       if ( my_id == 0 ) call boundary_from_grid(node_list, element_list, bnd_node_list, bnd_elm_list, output_bnd_elements)
       call broadcast_boundary(my_id, bnd_elm_list, bnd_node_list)
   
       call broadcast_vacuum(my_id, resistive_wall)

       call get_vacuum_response(my_id, node_list, bnd_elm_list, bnd_node_list, freeboundary_equil,    &
        resistive_wall)
       call import_external_fields('coil_field.dat', my_id)
       if ( .not. wall_curr_initialized ) call init_wall_currents(my_id, resistive_wall)
       first_step = .false.
    endif

    call MPI_Barrier(MPI_COMM_WORLD,ierr)
  
    !***********************************************************************
    !*              end intialisation                                      *
    !***********************************************************************
  
    if (.not. starwall_equil_coils) then
      write(*,*) ' **** Fatal: jorek2_wall_forces needs freeboundary n=0 and '
      write(*,*) '             starwall PF coils                             '
      stop
    endif 
  
    ! --- FORCES ---
    call total_wall_forces(my_id, node_list, element_list, scale_fact, Fx, Fy, Fz)
 
    if (my_id==0) then
      write(87,'(I5.5,5ES14.6)') istep, t_start, t_start*sqrt_mu0_rho0*1.d3, Fx, Fy, Fz 
    endif
  
    if (my_id==0) then
      write(*,*) ' Fx = ', Fx
      write(*,*) ' Fy = ', Fy
      write(*,*) ' Fz = ', Fz
    endif

  enddo

  close(87)

  call MPI_FINALIZE(IERR)                                ! clean up MPI

end program jorek2_wall_forces 

