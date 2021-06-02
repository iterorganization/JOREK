!< This program reads jorek_restart files and a "starwall_response.dat" file
!! and calculates the magnetic field at arbitrary xyz points. It involves the 
!! computation of expensive volume integrals in the plasma, so run it
!! in paralllel with several MPI processes for large cases. 
!! As input you must provide the points in cartesian coordinates in the file
!! xyz.nml in the format
!! 
!!   n (number_of_points)
!!   x_1       y_1       z_1  
!!   x_2       y_2       z_2  
!!    :         :         :   
!!    :         :         :   
!!   x_n       y_n       z_n  
!!
!! DO NOT specify the points exactly at the STARWALL's coil/wall
!! triangles. Otherwise singularities at those points may occur!
!!
!! To select the given JOREK restart files create a file named steps.nml
!! with 3 integers in the first line (istart, iend, delta_step)
!!   istart     = 1st  restart file index
!!   iend       = Last restart file index
!!   delta_step = index jump between the restart files
!!
!! The fields are exported in the file "fields_xyz.dat" and are also separated in 
!! coil/wall/plasma contributions. If you are running with freeboundary=.f. then
!! you only get the fields produced by the plasma currents
!!
!! Note that to calculate the require plasma fields n_plane typically should be
!! more than 60 (even for 2D)
!!
!! See more documentation here
!!   https://www.jorek.eu/wiki/doku.php?id=jorek2_fields_xyz
program jorek2_fields_xyz 

  use mod_vacuum_fields
  use mod_plasma_response,  only: plasma_fields_at_xyz
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
 
  integer   :: my_id, my_id_n, my_id_master, ierr, ierr2
  integer   :: i_rank(n_tor), n_cpu, n_cpu_n, n_cpu_master, m_cpu, n_masters, n_cpu_trans, my_id_trans
  integer   :: MPI_COMM_N, MPI_GROUP_MASTER, MPI_GROUP_WORLD, MPI_COMM_MASTER, MPI_COMM_TRANS
  integer   :: required,provided,StatInfo
  integer   :: istep, delta_step, istart, iend, np, i
  integer*4 :: rank, comm_size 
  logical   :: first_step

  real*8               :: bx,      by,      bz
  real*8, allocatable  :: bx_c(:), by_c(:), bz_c(:)
  real*8, allocatable  :: bx_w(:), by_w(:), bz_w(:)
  real*8, allocatable  :: bx_p(:), by_p(:), bz_p(:)
  real*8, allocatable  ::    x(:),    y(:),    z(:)

  character*17      :: file_in
  character(len=MPI_MAX_PROCESSOR_NAME) :: name
  integer :: resultlength
 
 
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

  open(26,file='xyz.nml',action='read',iostat=ierr)
  if (ierr/=0) then
    write(*,*) 'Could not read xyz.nml'
    stop
  endif
  read(26,*) np 
  allocate(x(np), y(np), z(np))
  do i=1, np
    read(26,*) x(i), y(i), z(i)
  enddo
  close(26)

  open(25,file='steps.nml',action='read',iostat=ierr)
  if (ierr/=0) then
    write(*,*) 'Could not read steps.nml'
    stop
  endif
  read(25,*) istart, iend, delta_step   
  close(25)

  if (my_id==0) then
    open(87,file='fields_xyz.dat',action='write')
    write(87,'(a)', advance='no') '#Step  '
    write(87,'(a)', advance='no') 'time(norm)    '
    write(87,'(a)', advance='no') 'time(ms)      '
    write(87,'(a)', advance='no') 'x(m)          '
    write(87,'(a)', advance='no') 'y(m)          '
    write(87,'(a)', advance='no') 'z(m)          '
    write(87,'(a)', advance='no') 'Bx(T)         '
    write(87,'(a)', advance='no') 'By(T)         '
    write(87,'(a)', advance='no') 'Bz(T)         '
    write(87,'(a)', advance='no') 'Bx_p(T)       '
    write(87,'(a)', advance='no') 'By_p(T)       '
    write(87,'(a)', advance='no') 'Bz_p(T)       '
    write(87,'(a)', advance='no') 'Bx_w(T)       '
    write(87,'(a)', advance='no') 'By_w(T)       '
    write(87,'(a)', advance='no') 'Bz_w(T)       '
    write(87,'(a)', advance='no') 'Bx_c(T)       '
    write(87,'(a)', advance='no') 'By_c(T)       '
    write(87,'(a)')               'Bz_c(T)       '
  endif

  allocate(bx_c(np), by_c(np), bz_c(np))
  allocate(bx_w(np), by_w(np), bz_w(np))
  allocate(bx_p(np), by_p(np), bz_p(np))

  ! --- Loop over restart files
  do istep = istart, iend, delta_step 

    write(file_in,'(A5,i5.5)') 'jorek', istep

    if ( my_id == 0 ) then
      call import_restart(node_list, element_list, file_in, rst_format, ierr)
    endif

    call MPI_BCAST(ierr,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr2)
    if ( ierr /= 0 ) cycle

    call broadcast_phys(my_id)  
    call broadcast_elements(my_id, element_list)                ! elements
    call broadcast_nodes(my_id, node_list)                      ! nodes
  
    if (.not. freeboundary) then
      write(*,*) ' **** Fatal: jorek2_fields_xyz needs freeboundary simulations ****'
      stop
    endif
  
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    
    ! --- Determine boundary information from the grid
    if ( my_id == 0 ) call boundary_from_grid(node_list, element_list, bnd_node_list, bnd_elm_list, output_bnd_elements)
    call broadcast_boundary(my_id, bnd_elm_list, bnd_node_list)
   
    if (freeboundary)  call broadcast_vacuum(my_id, resistive_wall)

    ! --- Fill the vacuum response matrices for freeboundary computations
    if (first_step .and. freeboundary) then
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
  
    if (freeboundary) then
      call wall_fields_at_xyz(my_id, x, y, z, bx_w, by_w, bz_w)
      if (sr%ncoil > 0) then
        if (.not. starwall_equil_coils) then
          if (my_id==0) then
            write(*,*) '******************************************************'
            write(*,*) '*** WARNING: You are running without equilibrium coils'
            write(*,*) '***    n=0 component of PF coils fields is missing!   '
            write(*,*) '******************************************************'
          endif
        endif
        call coil_fields_at_xyz(my_id, x, y, z, bx_c, by_c, bz_c)
      else
        if (my_id==0) then
          write(*,*) '******************************************************'
          write(*,*) '*** WARNING: You are running without coils.           '
          write(*,*) '***    coil fields are set to 0 !                     '
          write(*,*) '******************************************************'
        endif
        bx_c = 0.d0;    by_c = 0.d0;    bz_c = 0.d0;
      endif
    else
      if (my_id==0) then
        write(*,*) '******************************************************'
        write(*,*) '*** WARNING: You are running with freeboundary=.f.    '
        write(*,*) '***    coils and wall fields are set to 0 !           '
        write(*,*) '******************************************************'
      endif
      bx_c = 0.d0;    by_c = 0.d0;    bz_c = 0.d0;
      bx_w = 0.d0;    by_w = 0.d0;    bz_w = 0.d0;
    endif 

    call plasma_fields_at_xyz(my_id, node_list,element_list, x, y, z, bx_p, by_p, bz_p)
 
    if (my_id==0) then
      do i=1, np
        bx = bx_p(i) + bx_w(i) + bx_c(i)
        by = by_p(i) + by_w(i) + by_c(i)
        bz = bz_p(i) + bz_w(i) + bz_c(i)
        write(87,'(I5.5,17ES14.6)') istep, t_start, t_start*sqrt_mu0_rho0*1.d3,x(i),y(i),z(i), bx, by, bz, &
                                    bx_p(i), by_p(i), bz_p(i), bx_w(i), by_w(i), bz_w(i), bx_c(i), by_c(i), bz_c(i)      
      enddo
    endif
    
    write(87,*) ' '
    write(87,*) ' '

  enddo

  close(87)

  call MPI_FINALIZE(IERR)                                ! clean up MPI

end program jorek2_fields_xyz 

