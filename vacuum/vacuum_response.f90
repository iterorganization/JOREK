!> Implements the interplay of the plasma with a conducting wall.
!!
!! The plasma-wall interaction is characterized by vacuum response matrices which are calculated by
!! the STARWALL code and imported into JOREK by the routine read_starwall_response(). The matrices
!! allow to express the magnetic field parallel to the interface (boundary of the JOREK domain) in
!! terms of the poloidal flux at the interface and the wall currents. The vacuum
!! response enters into the boundary integral of the current equation which vanishes for fixed
!! boundary conditions. The boundary integral is implemented in the routine
!! vacuum_boundary_integral().
!!
!! @note The variable s in a boundary element may correspond to s or t of the respective
!! 2D element depending on element orientation.
module vacuum_response
  
  use vacuum
  
  implicit none
  
  
  
  contains
  
  
  
  !> Is the filehandle associated with a formatted file?
  logical function is_formatted(filehandle)
    implicit none
    ! --- Routine parameters
    integer,          intent(in) :: filehandle
    ! --- Local variables
    character(len=64) :: format_type
    
    inquire(filehandle,form=format_type)
    is_formatted = ( trim(format_type) == 'FORMATTED' )
  end function is_formatted
  
  
  
  !> Get the vacuum response for an ideal or resistive wall.
  subroutine get_vacuum_response(my_id, node_list, bnd_elm_list, bnd_node_list,                    &
    freeboundary_equil, resistive_wall)
  
    use mpi_mod
    use constants
    use mod_parameters, only: n_tor
    use data_structure, only: type_node_list, type_bnd_element_list, type_bnd_node_list
    use phys_module,    only: central_mass, central_density

    implicit none    
    
    integer,                     intent(in) :: my_id              !< MPI proc ID
    type(type_node_list),        intent(in) :: node_list          !< List of boundary nodes
    type(type_bnd_element_list), intent(in) :: bnd_elm_list       !< List of boundary elements
    type(type_bnd_node_list),    intent(in) :: bnd_node_list      !< List of boundary nodes
    logical,                     intent(in) :: freeboundary_equil !< Use free boundary equilibrium?
    logical,                     intent(in) :: resistive_wall     !< Resistive or ideal wall?

    integer :: i,j, ierr, dim
    logical :: exists

    ! --- Determine total number of boundary degrees of freedom per harmonic (skipping duplicates).
    do i=1, bnd_node_list%n_bnd_nodes
      exists = .false.
      do j=1, i-1
        if (bnd_node_list%bnd_node(i)%index_jorek .eq. bnd_node_list%bnd_node(j)%index_jorek) then
          exists  = .true.
          exit
        end if
      end do
      if (.not. exists) then
        n_dof_bnd = n_dof_bnd + bnd_node_list%bnd_node(i)%n_dof ! Number of boundary degrees of freedom per harmonic
      end if
    end do

    if (my_id == 0) write(*,'(A,i5)') '   total number of degrees of freedom on the boundary : ',n_dof_bnd
    
    ! --- Write out the boundary information for STARWALL.
    if (my_id == 0) call export_boundary(node_list, bnd_elm_list, bnd_node_list)
    
    if ( vacuum_debug  ) call log_starwall_response(my_id, sr) 
    call read_starwall_response(my_id, sr,'starwall-response.dat',bnd_elm_list%n_bnd_elements)
    call broadcast_starwall_response(my_id, sr)
    if( sr%ncoil .gt. 0 ) call distribute_coil_names

    ! --- Set the "wall resistivity" to be used inside JOREK (actually it is the normalized thin wall resistivity)
    if ( (sr%file_version == 1)  .and. (my_id == 0) ) then
      write(*,*) 'Remark: STARWALL response file_version==1 means that wall_resistivity is specified in the JOREK namelist file.'
      write(*,*) '        Thus, the input parameter wall_resistivity_fact is ignored (WARNING).'
    else
      if (my_id ==0) then
          write(*,*) 'Remark: STARWALL response file_version>=2 means that eta_thin_w is specified in the STARWALL input.'
          write(*,*) '        The JOREK variable wall_resistivity is automatically calculated from it.'
          write(*,*) '        Thus, the input parameter wall_resistivity is ignored (WARNING).'
      end if
      wall_resistivity = wall_resistivity_fact * sr%eta_thin_w * &
        sqrt( central_density * 1.d20 * central_mass * mass_proton / mu_zero )
    end if
   
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    call log_starwall_response(my_id, sr)
              
  end subroutine get_vacuum_response



  !> Read an integer parameter from the STARWALL response file.
  integer function read_intparam(filehandle, parameter_name)

    implicit none

    ! --- Routine parameters
    integer,          intent(in) :: filehandle
    character(len=*), intent(in) :: parameter_name

    ! --- Local variables
    character(len=12) :: marker
    character(len=24) :: name
    integer           :: ierr

    if ( is_formatted(filehandle) ) then
      read(filehandle,'(A12,A24,I12)',iostat=ierr) marker, name, read_intparam
    else
      read(filehandle, iostat=ierr) marker, name, read_intparam
    end if

    if ( (ierr /= 0) .or. (trim(adjustl(marker)) /= '#@intparam') .or. (trim(adjustl(name)) /= trim(parameter_name)) ) then
      write(*,*) 'ERROR: Could not read parameter "', trim(parameter_name) ,'"from STARWALL response.'
      stop
    end if

    if ( vacuum_debug ) write(*,'(3x,"Read: ",A24,"=",I12)',iostat=ierr) name,read_intparam

  end function read_intparam



  !> Read an array from the STARWALL respone file
  subroutine read_array(filehandle, array_name, dim, int1d, int2d, float1d, float2d)

    implicit none

    ! --- Routine parameters
    integer, intent(in) :: filehandle
    character(len=*), intent(in) :: array_name
    integer, intent(in) :: dim(2)
    integer, allocatable, optional, intent(inout)   :: int1d(:)
    integer, allocatable, optional, intent(inout)   :: int2d(:,:)
    real*8,  allocatable, optional, intent(inout)   :: float1d(:)
    real*8,  allocatable, optional, intent(inout)   :: float2d(:,:)

    ! --- Local variables
    character(len=12) :: marker
    integer           :: nd, d(2), ierr
    character(len=24) :: name, datatype, requested_type
    logical           :: error

    if ( present(int1d) .or. present(int2d) ) then
      requested_type = 'int'
    else
      requested_type = 'float'
    end if

    if ( is_formatted(filehandle) ) then
      read(filehandle,'(A12,A24,I12,A24,2I12)',iostat=ierr) marker, name, nd, datatype, d
    else
      read(filehandle,iostat=ierr) marker, name, nd, datatype, d
    end if
    marker   = adjustl(marker)
    name     = adjustl(name)
    datatype = adjustl(datatype)

    error = ( ierr /= 0 ) .or. ( trim(marker) /= '#@array' ) .or. ( trim(name) /= trim(array_name) )                     &
      .or. ( dim(1) /= d(1) ) .or. ( dim(2) /= d(2) ) .or. ( trim(datatype) /= trim(requested_type) )

    if ( error ) then
      write(*,*) 'ERROR: Could not read array ', trim(array_name), ' from STARWALL response.'
      stop
    end if

    if ( present(int1d) ) then

      if ( allocated(int1d) ) deallocate( int1d )
      allocate( int1d(dim(1)) )
      if ( is_formatted(filehandle) ) then
        read(filehandle,*) int1d(:)
      else
        read(filehandle) int1d(:)
      end if

    else if ( present(int2d) ) then

      if ( allocated(int2d) ) deallocate( int2d )
      allocate( int2d(dim(1),dim(2)) )
      if ( is_formatted(filehandle) ) then
        read(filehandle,*) int2d(:,:)
      else
        read(filehandle) int2d(:,:)
      end if

    else if ( present(float1d) ) then

      if ( allocated(float1d) ) deallocate( float1d )
      allocate( float1d(dim(1)) )
      if ( is_formatted(filehandle) ) then
        read(filehandle,*) float1d(:)
      else
        read(filehandle) float1d(:)
      end if

    else if ( present(float2d) ) then
      
      if ( allocated(float2d) ) deallocate( float2d )
      allocate( float2d(dim(1),dim(2)) )
      if ( is_formatted(filehandle) ) then
        read(filehandle,'(4ES24.16)') float2d(:,:)
      else
        read(filehandle) float2d(:,:)
      end if

    end if

    if ( vacuum_debug ) write(*,'(3x,"Read: ",A24,"> type ",A," size ",2I7)')name, trim(datatype), d(1:nd)

  end subroutine read_array

  

  !> Read an integer parameter from a STARWALL response file (MPI I/O).
  integer function read_intparam_parallel(filehandle, parameter_name,disp)

    use mpi_mod
    implicit none

    ! --- Routine parameters
    integer,                       intent(in)    :: filehandle
    character(len=*),              intent(in)    :: parameter_name
    integer(kind=MPI_OFFSET_KIND), intent(inout) :: disp !< Present location in the file

    ! --- Local variables
    integer, dimension (MPI_STATUS_SIZE) :: status
    character(len=12) :: marker
    character(len=24) :: name
    integer           :: err

    call MPI_FILE_READ(filehandle, marker,int(sizeof(marker),4),MPI_CHARACTER,status,err)
    call MPI_FILE_READ(filehandle, name,  int(sizeof(name),4)  ,MPI_CHARACTER,status,err)
    call MPI_FILE_READ(filehandle, read_intparam_parallel, 1,MPI_INTEGER,status,err)
    disp = disp + sizeof(marker) + sizeof(name) + sizeof(10)! any integer 

    if ( (err /= 0) .or. (trim(adjustl(marker)) /= '#@intparam') .or. (trim(adjustl(name)) /= trim(parameter_name)) ) then
      write(*,*) 'ERROR: Could not read parameter "', trim(parameter_name) ,'"from STARWALL response.'
      stop
    end if
    if ( vacuum_debug ) write(*,'(3x,"Read: ",A24,"=",I12)',iostat=err) name, read_intparam_parallel

  end function read_intparam_parallel  
 

 
  
  !> Read an array from the STARWALL respone file (MPI I/O).
  subroutine read_array_not_distr(filehandle, array_name, dim, disp, int1d, int2d, float1d, float2d, char1d)

    use mpi_mod
    implicit none

    ! --- Routine parameters
    integer,                        intent(in)      :: filehandle
    character(len=*),               intent(in)      :: array_name
    integer,                        intent(in)      :: dim(2)
    integer(kind=MPI_OFFSET_KIND),  intent(inout)   :: disp !< Present location in the file
    integer, allocatable, optional, intent(inout)   :: int1d(:)
    integer, allocatable, optional, intent(inout)   :: int2d(:,:)
    real*8,  allocatable, optional, intent(inout)   :: float1d(:)
    real*8,  allocatable, optional, intent(inout)   :: float2d(:,:)
    character(len=*), allocatable, optional, intent(inout)   :: char1d(:)

    ! --- Local variables
    integer, dimension (MPI_STATUS_SIZE) :: status
    character(len=12) :: marker
    integer           :: nd, d(2), ierr, err, str_len
    character(len=24) :: name, datatype, requested_type
    logical           :: error
    character, allocatable :: tmp_char1d(:)
    
    if ( present(int1d) .or. present(int2d) ) then
      requested_type = 'int'
    elseif (present(char1d))then
      requested_type = 'char'
    else
      requested_type = 'float'
    end if

    call MPI_FILE_READ(filehandle,  marker,   int(sizeof(marker),4),    MPI_CHARACTER,status,ierr)
    call MPI_FILE_READ(filehandle,  name,     int(sizeof(name),4),      MPI_CHARACTER,status,ierr)
    call MPI_FILE_READ(filehandle,  nd,       1,                        MPI_INTEGER  ,status,ierr)
    call MPI_FILE_READ(filehandle,  datatype, int(sizeof(datatype),4),  MPI_CHARACTER,status,ierr)
    call MPI_FILE_READ(filehandle,  d,        2,                        MPI_INTEGER  ,status,ierr)
    disp = disp + sizeof(marker) + sizeof(name) + sizeof(nd) + sizeof(datatype) + sizeof(d) 

    marker   = adjustl(marker)
    name     = adjustl(name)
    datatype = adjustl(datatype)

    error = ( ierr /= 0 ) .or. ( trim(marker) /= '#@array' ) .or. ( trim(name) /= trim(array_name) ) &
      .or. ( dim(1) /= d(1) ) .or. ( dim(2) /= d(2) ) .or. ( trim(datatype) /= trim(requested_type) )

    if ( error ) then
      write(*,*) 'ERROR: Could not read array ', trim(array_name), ' from STARWALL response.'
      stop
    end if
 
    if ( present(int1d) ) then

      if ( allocated(int1d) ) deallocate( int1d ); allocate( int1d(dim(1)) )
      call MPI_FILE_READ(filehandle, int1d, dim(1), MPI_INTEGER  ,status,ierr)
      disp = disp + sizeof(int1d)

    else if ( present(int2d) ) then

      if ( allocated(int2d) ) deallocate( int2d ); allocate( int2d(dim(1),dim(2)) )
      call MPI_FILE_READ(filehandle, int2d, dim(1)*dim(2), MPI_INTEGER  ,status,ierr)
      disp = disp + sizeof(int2d)

    else if ( present(float1d) ) then

      if ( allocated(float1d) ) deallocate( float1d ); allocate( float1d(dim(1)) )
      call MPI_FILE_READ(filehandle, float1d, dim(1), MPI_DOUBLE_PRECISION  ,status,ierr)
      disp = disp + sizeof(float1d)

    else if ( present(float2d) ) then

      if ( allocated(float2d) ) deallocate( float2d ); allocate( float2d(dim(1),dim(2)) )
      call MPI_FILE_READ(filehandle, float2d, dim(1)*dim(2), MPI_DOUBLE_PRECISION  ,status,ierr)
      disp = disp + sizeof(float2d)

    else if (present(char1d)) then
      str_len = dim(1) * sizeof(char1d(1))
      if ( allocated(tmp_char1d) ) deallocate( tmp_char1d ); allocate( tmp_char1d(str_len))
      if ( allocated(char1d) ) deallocate( char1d ); allocate( char1d(dim(1)))
      call MPI_FILE_READ(filehandle, tmp_char1d, str_len, MPI_CHARACTER  ,status,ierr)
      char1d = transfer(tmp_char1d, char1d)
      deallocate(tmp_char1d)
      disp = disp + str_len
    end if

    if ( vacuum_debug ) write(*,'(3x,"Read: ",A24,"> type ",A," size ",2I7)')name, trim(datatype), d(1:nd)

  end subroutine read_array_not_distr
  
  
  
  !> Read an array from the STARWALL respone file
  subroutine read_array_par(filehandle, array_name, dim, disp, my_id, float2d, row_wise)

    use mpi_mod
    implicit none

    ! --- Routine parameters
    integer,                        intent(in)      :: filehandle
    character(len=*),               intent(in)      :: array_name
    integer,                        intent(in)      :: dim(2)
    integer(kind=MPI_OFFSET_KIND),  intent(inout)   :: disp !< Present location in the file
    integer,                        intent(in)      :: my_id
    type(t_distrib_mat),            intent(inout)   :: float2d
    logical,                        intent(in)      :: row_wise ! if  .true. -  rowwise reading; .false. - columnwise    

    ! --- Local variables
    integer, dimension (MPI_STATUS_SIZE) :: status
    character(len=12) :: marker
    integer           :: nd,d(2),loc_sizes(2),loc_starts(2),ierr,err,ntasks
    character(len=24) :: name, datatype
    logical           :: error
    integer           :: my_subarray
    integer           :: num_read_elements, n_step_read, step_read
    integer(KIND=8)   :: local_num_elements, num_read_elements_const, i, i_ind, j_ind
    
    if (my_id==0) then
      call MPI_FILE_READ(filehandle, marker,   int(sizeof(marker), 4),   MPI_CHARACTER, status, ierr)
      call MPI_FILE_READ(filehandle, name,     int(sizeof(name), 4),     MPI_CHARACTER, status, ierr)
      call MPI_FILE_READ(filehandle, nd,       1,                        MPI_INTEGER,   status, ierr)
      call MPI_FILE_READ(filehandle, datatype, int(sizeof(datatype), 4), MPI_CHARACTER, status, ierr)
      call MPI_FILE_READ(filehandle, d,        2,                        MPI_INTEGER,   status, ierr)

      disp = disp + sizeof(marker) + sizeof(name) + sizeof(nd) + sizeof(datatype) + sizeof(d)
      marker   = adjustl(marker)
      name     = adjustl(name)
      datatype = adjustl(datatype)

      error = ( ierr /= 0 ) .or. ( trim(marker) /= '#@array' ) .or. (trim(name) /= trim(array_name) ) &
        .or. ( dim(1) /= d(1) ) .or. ( dim(2) /= d(2) ) .or. (trim(datatype) /= 'float' )

      if ( error ) then
        write(*,*) 'ERROR: Could not read array ', trim(array_name), ' from STARWALL response.'
        stop
      end if
    end if 
          
    call MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, ierr)  
    call MPI_BCAST(disp, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    if(row_wise) then ! If we read matrix rowwise
      float2d%step       = dim(1) / ntasks
      loc_starts(:)      = (/ my_id*float2d%step, 0 /)
      loc_sizes(:)       = (/ float2d%step, dim(2) /)
      if (my_id==ntasks-1) loc_sizes(1) = float2d%step + dim(1)-ntasks*float2d%step
      float2d%row_wise   = .true.
    else             ! If we read matrix columnwise
      float2d%step       = dim(2) / ntasks
      loc_starts(:)      = (/ 0, my_id*float2d%step /)
      loc_sizes(:)       = (/ dim(1), float2d%step /)
      if (my_id==ntasks-1) loc_sizes(2) = float2d%step + dim(2)-ntasks*float2d%step
      float2d%row_wise   = .false.
    endif

    local_num_elements = int(loc_sizes(1),8) * int(loc_sizes(2),8)

    call  alloc_distr(my_id, float2d, dim , row_wise)

    call MPI_TYPE_CREATE_SUBARRAY(2,dim,loc_sizes,loc_starts,MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,my_subarray,ierr)
    call MPI_Type_commit(my_subarray,ierr)
    call MPI_File_set_view(filehandle, disp, MPI_DOUBLE_PRECISION, my_subarray, "native", MPI_INFO_NULL, ierr)

    ! -----------------------------------
    ! The following is a work around for an Intel MPI limitation not allowing to read more than 2 GB
    ! of data per MPI task in each call corresponding to 268435455 double precision values.
    num_read_elements_const = 250000000

    if (num_read_elements_const>local_num_elements) then
      num_read_elements = local_num_elements
      n_step_read       = 1
    else
      n_step_read       = local_num_elements/num_read_elements_const + 1
      num_read_elements = num_read_elements_const
    end if

    call MPI_AllREDUCE(MPI_IN_PLACE,n_step_read,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)

    i_ind=1
    j_ind=1

    do i = 1, n_step_read
      call MPI_File_read_all(filehandle, float2d%loc_mat(i_ind, j_ind),num_read_elements, &
        MPI_DOUBLE_PRECISION, status, ierr)

      ! Calculation of starting indices of each slice in local array     
      j_ind = int(num_read_elements,8)*i/loc_sizes(1) + 1
      i_ind = int(num_read_elements,8)*i-int((j_ind-1),8)*int(loc_sizes(1),8) + 1

      ! If number of reading step is not eaqual for each MPI task.         
      if (num_read_elements_const*int((i+1),8)>local_num_elements) then
        num_read_elements = local_num_elements-i*num_read_elements_const
        if (num_read_elements <= 0) then
          num_read_elements = 0
          j_ind = 1
          i_ind = 1
        endif
      endif
    end do
    ! End of workaround for Intel MPI limitation
    ! -----------------------------------
      
    !The following line can replace the workaround in the future
    !call MPI_File_read_all(filehandle, float2d%loc_mat, loc_sizes(1)*loc_sizes(2), MPI_DOUBLE_PRECISION,status, ierr)

    disp = disp + sizeof(float2d%loc_mat(1,1))*dim(1)*dim(2)

    ! Return to ordinary file view
    call MPI_File_set_view(filehandle, disp, MPI_BYTE, MPI_BYTE,"native",MPI_INFO_NULL, ierr)

    if ( vacuum_debug  .and. (my_id==0) ) write(*,'(3x,"Read: ",A24,"> type ",A," size",2I7)')name,&
      trim(datatype), d(1:nd)

  end subroutine read_array_par
  
  !===========================================================================================================
  
  !> Read the STARWALL response matrices from a single file.
  !!
  !! file_version 1: Original
  !! file_version 2: Includes eta_thin_w
  !! file_version 3: Includes additional coil information
  subroutine read_starwall_response(my_id, sr, filename, n_bnd)

    use constants
    use mod_parameters, only: n_tor, n_period
    use mpi_mod

    implicit none

    ! --- Routine parameters
    integer,                   intent(in)    :: my_id
    type(t_starwall_response), intent(inout) :: sr
    character(len=*),          intent(in)    :: filename
    integer,                   intent(in)    :: n_bnd !< Number of boundary elements for consistency check

    ! --- Local variables
    integer            :: filehandle
    character(len=512) :: comment, tmp_int2str
    integer            :: i, j, i_starw, n, is_sin, err, i_tmp, ier
    real*8             :: r_tmp
    real*8, allocatable:: tmp(:)
    integer, dimension (MPI_STATUS_SIZE) :: status
    integer(kind=MPI_OFFSET_KIND)        :: disp !< Present location in the file
    real*8             :: test_sum
    integer            :: loc_sizes(2),ntasks

    disp = 0
    call MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, err)

    ! --- Open file
    if (my_id ==0 ) write(*,*) 'Trying to open response file...'

    call MPI_FILE_OPEN(MPI_COMM_WORLD, 'starwall-response.dat', MPI_MODE_RDONLY, MPI_INFO_NULL, filehandle, err)
    if ( (err == 0) .and. (my_id == 0) ) then
      call MPI_FILE_READ(filehandle, i_tmp, 1, MPI_INTEGER,          status, err)
      call MPI_FILE_READ(filehandle, r_tmp, 1, MPI_DOUBLE_PRECISION, status, err)
      disp = disp + sizeof(42) + sizeof(42.d0)
      if ( (i_tmp/=42) .or. (r_tmp/=42.d0) ) err=-42         
    end if
    call MPI_BCAST(err, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ier )
    if (err/=0) call MPI_FILE_CLOSE(filehandle, err)

    if ( (err /= 0) .and. (my_id == 0) ) then
      write(*,*) '  ... failed.'
      write(*,*) 'ERROR: STARWALL response file (',trim(filename),') could not be opened.'
      stop
    end if
    if (my_id ==0) write(*,*) '  ... succeeded.'

    ! --- Read data from STARWALL response file
    if (my_id == 0) then

      call MPI_FILE_READ(filehandle, comment,int(sizeof(comment),4),MPI_CHARACTER,status,err)
      disp = disp + sizeof(comment)

      sr%file_version = read_intparam_parallel(filehandle, 'file_version', disp)
      if ( sr%file_version > 5 ) then
        write(*,*) 'ERROR: STARWALL response file version ', sr%file_version, ' is not supported.'
        stop
      end if

      sr%n_bnd  = read_intparam_parallel(filehandle, 'n_bnd' , disp)
      if ( n_bnd /= sr%n_bnd ) then
        write(*,*) 'ERROR: The number of boundary elements in the STARWALL response file is different from your grid.'
        stop
      end if
      
      sr%nd_bez = read_intparam_parallel(filehandle, 'nd_bez', disp)

      if ( sr%file_version >= 5 ) then
        sr%nv       = read_intparam_parallel(filehandle, 'nv',        disp)
        sr%n_points = read_intparam_parallel(filehandle, 'n_points',  disp)
      endif

      sr%ncoil  = read_intparam_parallel(filehandle, 'ncoil' , disp)
      sr%npot_w = read_intparam_parallel(filehandle, 'npot_w', disp)
      sr%n_w    = read_intparam_parallel(filehandle, 'n_w'   , disp)
      sr%ntri_w = read_intparam_parallel(filehandle, 'ntri_w', disp)

      if ( sr%file_version >= 5 ) then
        sr%iwall    = read_intparam_parallel(filehandle, 'iwall',    disp)
        sr%nwu      = read_intparam_parallel(filehandle, 'nwu',      disp)
        sr%nwv      = read_intparam_parallel(filehandle, 'nwv' ,     disp)
        sr%mn_w     = read_intparam_parallel(filehandle, 'mn_w',     disp)
        sr%max_mn_w = read_intparam_parallel(filehandle, 'MAX_MN_W', disp)

        call read_array_not_distr(filehandle, 'm_w',         (/sr%max_mn_w,0/),  disp,  int1d=sr%m_w)
        call read_array_not_distr(filehandle, 'n_w_fourier', (/sr%max_mn_w,0/),  disp,  int1d=sr%n_w_fourier)
        call read_array_not_distr(filehandle, 'rc_w',        (/sr%max_mn_w,0/),  disp,  float1d=sr%rc_w)
        call read_array_not_distr(filehandle, 'rs_w',        (/sr%max_mn_w,0/),  disp,  float1d=sr%rs_w)
        call read_array_not_distr(filehandle, 'zc_w',        (/sr%max_mn_w,0/),  disp,  float1d=sr%zc_w)
        call read_array_not_distr(filehandle, 'zs_w',        (/sr%max_mn_w,0/),  disp,  float1d=sr%zs_w)
      endif

      sr%n_tor  = read_intparam_parallel(filehandle, 'n_tor' , disp)
      sr%n_tor0 = sr%n_tor

      call read_array_not_distr(filehandle, 'i_tor', (/sr%n_tor,0/), disp, int1d=sr%i_tor)

      if ( sr%file_version >= 3 ) then
        
        sr%ntri_c                  = read_intparam_parallel(filehandle, 'ntri_c',          disp)
        sr%n_pol_coils             = read_intparam_parallel(filehandle, 'n_pol_coils',     disp)
        sr%n_rmp_coils             = read_intparam_parallel(filehandle, 'n_rmp_coils',     disp)
        sr%n_voltage_coils         = read_intparam_parallel(filehandle, 'n_voltage_coils', disp)
        sr%n_diag_coils            = read_intparam_parallel(filehandle, 'n_diag_coils',    disp)
        
        if (sr%n_voltage_coils > 0 ) then
          write(*,*) 'ERROR: voltage_coils not yet implemented.'
          stop
        end if
      
        if ( sr%ncoil /= sr%n_pol_coils + sr%n_rmp_coils + sr%n_voltage_coils + sr%n_diag_coils) then
          write(*,*) 'ERROR: STARWALL response is inconsistent: ncoil does not match sum.'
          stop
        end if
 
        sr%ind_start_pol_coils     = read_intparam_parallel(filehandle, 'ind_start_pol_coils',     disp)
        sr%ind_start_rmp_coils     = read_intparam_parallel(filehandle, 'ind_start_rmp_coils',     disp)
        sr%ind_start_voltage_coils = read_intparam_parallel(filehandle, 'ind_start_voltage_coils', disp)
        sr%ind_start_diag_coils    = read_intparam_parallel(filehandle, 'ind_start_diag_coils',    disp)

        if ( sr%ncoil > 0 ) then
          call read_array_not_distr(filehandle, 'jtri_c',        (/sr%ncoil,0/),  disp,  int1d=sr%jtri_c)
          call read_array_not_distr(filehandle, 'x_coil',        (/sr%ntri_c,3/), disp,  float2d=sr%x_coil)
          call read_array_not_distr(filehandle, 'y_coil',        (/sr%ntri_c,3/), disp,  float2d=sr%y_coil)
          call read_array_not_distr(filehandle, 'z_coil',        (/sr%ntri_c,3/), disp,  float2d=sr%z_coil)
          call read_array_not_distr(filehandle, 'phi_coil',      (/sr%ntri_c,3/), disp,  float2d=sr%phi_coil)
          call read_array_not_distr(filehandle, 'eta_thin_coil', (/sr%ntri_c,0/), disp,  float1d=sr%eta_thin_coil)
          call read_array_not_distr(filehandle, 'coil_resist',   (/sr%ncoil,0/),  disp,  float1d=sr%coil_resist)
          if (sr%file_version > 3)  then
            call read_array_not_distr(filehandle, 'coil_name',     (/sr%ncoil,0/),  disp,  char1d=sr%coil_name)
           endif
        end if
        
        if ( (sr%ncoil .gt. 0) .and. sr%file_version .le. 3) then
         write(*,*) "Coil names not yet supported by starwall (ver <=3). Generic names used."
         if (allocated(sr%coil_name)) deallocate(sr%coil_name)
         allocate(sr%coil_name(sr%ncoil))

         do i=1, sr%n_diag_coils
           write(tmp_int2str,'(I3.3)'), i
           sr%coil_name(sr%ind_start_diag_coils + i - 1) = trim("Dia_"//trim(tmp_int2str))
         end do

         do i=1, sr%n_pol_coils
           write(tmp_int2str,'(I3.3)'), i
           sr%coil_name(sr%ind_start_pol_coils + i - 1) = trim("PF_"//trim(tmp_int2str))
         end do
         
         do i=1, sr%n_rmp_coils
           write(tmp_int2str,'(I3.3)'), i
           sr%coil_name(sr%ind_start_rmp_coils + i - 1) = trim("RMP_"//trim(tmp_int2str))
         end do
       endif
     end if

      ! --- eta_thin_w is only part of the STARWALL response file since file_version 2
      if ( sr%file_version >= 2 ) then
        allocate(tmp(1))
        call read_array_not_distr(filehandle, 'eta_thin_w', (/1,0/), disp, float1d=tmp)
        sr%eta_thin_w = tmp(1)
        deallocate(tmp)
      else
        sr%eta_thin_w = 0.
      end if

      call read_array_not_distr(filehandle, 'yy', (/sr%n_w,0/), disp, float1d=sr%d_yy)
   
    end if

    if(my_id==0) call memory_prediction(sr%n_w, sr%nd_bez, ntasks)

    call MPI_BCAST(sr%n_w,    1, MPI_INTEGER, 0, MPI_COMM_WORLD, err)
    call MPI_BCAST(sr%nd_bez, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, err)

    ! Last parameter defines how matrix should be distributed: .false. - rowwise; .true. - columnwise
    call read_array_par    (filehandle, 'ye',       (/sr%n_w,sr%nd_bez/),    disp, my_id,  sr%a_ye,     .false.)    
    call read_array_par    (filehandle, 'ey',       (/sr%nd_bez,sr%n_w/),    disp, my_id,  sr%a_ey,     .true. )
    call read_array_par    (filehandle, 'ee',       (/sr%nd_bez,sr%nd_bez/), disp, my_id,  sr%a_ee,     .true. )
    call read_array_par    (filehandle, 's_ww',     (/sr%n_w,sr%n_w/),       disp, my_id,  sr%s_ww,     .true. )
    call read_array_par    (filehandle, 's_ww_inv', (/sr%n_w,sr%n_w/),       disp, my_id,  sr%s_ww_inv, .true. )

    if(my_id == 0) then
      call read_array_not_distr(filehandle, 'xyzpot_w', (/sr%npot_w,3/), disp, float2d=sr%xyzpot_w)
      call read_array_not_distr(filehandle, 'jpot_w',   (/sr%ntri_w,3/), disp, int2d=sr%jpot_w)
      if ( sr%file_version >= 5 ) then
        call read_array_not_distr(filehandle, 'phi0_w',   (/sr%ntri_w,3/), disp, float2d=sr%phi0_w)
      endif
    end if

    call MPI_FILE_CLOSE(filehandle, err)

    if ( vacuum_debug .and. (my_id==0) ) write(*,*) 'Finished reading vacuum response.'

    ! --- Import normalization
    sr%a_ee%loc_mat(:,:) = sr%a_ee%loc_mat(:,:) * 2.d0*PI
    sr%a_ye%loc_mat(:,:) = sr%a_ye%loc_mat(:,:) * 2.d0*PI
    if ( vacuum_debug .and. (my_id==0) ) write(*,*) 'Applied import normalization.'

    ! --- Compute ideal-wall and no-wall response matrices.
    
    call  alloc_distr(my_id, sr%a_nw,(/sr%nd_bez, sr%nd_bez/), .true.)

    sr%a_nw%loc_mat(:,:) = sr%a_ee%loc_mat(:,:)
    call matrix_multiplication(my_id,sr%a_ey,mat2=sr%a_ye, res_mat=sr%a_id )
    sr%a_id%loc_mat(:,:) = sr%a_ee%loc_mat(:,:) - sr%a_id%loc_mat(:,:)

    ! --- Transform STARWALL harmonics to account for periodicity
    if ( my_id==0 ) then
      j = 0
      do i = 1, sr%n_tor
        i_starw = sr%i_tor(i)
        n       = i_starw / 2
        is_sin  = i_starw - 2 * n
        i_starw = 2 * n/n_period + is_sin
        if ( (mod(n, n_period) /= 0) .or. (i_starw < 1) .or. (i_starw > n_tor) ) then
          write(*,*) 'WARNING: STARWALL harmonic has no JOREK equivalent!'
          write(*,*) 'i_starw    =', sr%i_tor(i)
          write(*,*) 'n_period   =', n_period
          write(*,*) 'n_tor      =', n_tor
        else
          j = j + 1
          sr%i_tor(j) = i_starw
        end if
      end do
     
      sr%n_tor = j
      if ( vacuum_debug) write(*,*) 'End of routine read_starwall_response.'
    end if

  end subroutine read_starwall_response
  
  
  
  !> Routine for multiplying two (distributed) matrices.
  subroutine matrix_multiplication(my_id, mat1, mat2, mat2_not_distr, res_mat, res_mat_not_distr)
    
    use mpi_mod
  
    implicit none
  
    ! --- Routine parameters
    integer,                        intent(in)      :: my_id
    type(t_distrib_mat),            intent(in)      :: mat1
    type(t_distrib_mat), optional,  intent(in)      :: mat2
    real*8, allocatable, optional,  intent(inout)   :: mat2_not_distr(:,:)
    type(t_distrib_mat), optional,  intent(inout)   :: res_mat
    real*8, allocatable, optional,  intent(inout)   :: res_mat_not_distr(:,:)
  
    ! --- Local variables
    integer             :: ntasks, ierr, i, k, j, z, length, glob_index_i, glob_index_j 
    real*8, allocatable :: tmp(:,:)
    real*8              :: sum_element
    
    call MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, ierr)
  
    ! --- Multiplication mat1(distributed)*mat2(distributed) =res_mat(distributed)
    if ( present(res_mat) .and. present(mat2) ) then
  
      call  alloc_distr(my_id, res_mat, (/mat1%dim(1), mat2%dim(2)/), .true.)       
      res_mat%loc_mat=0.0
    
      if ( res_mat%row_wise .and. mat1%row_wise .and. (.not. mat2%row_wise) ) then
        
       if(mat1%dim(2) /= mat2%dim(1) ) then
          write(*,*) "ERROR: Dimensions of  multiplied matrices (1) #col =", mat1%dim(2), &
                     " #row =", mat2%dim(1) , "STOP (1)"
          stop
       endif   
 
        ! loop over all tasks
        do i = 1, ntasks
    
          ! Broadcasting size of the distributed column
          length = size(mat2%loc_mat,2)      
          call MPI_BCAST(length, 1, MPI_INTEGER, i-1, MPI_COMM_WORLD, ierr )
    
          ! temporal matrix for receiving part of the distributed matrix and
          ! multiply it
          if ( allocated(tmp) ) deallocate(tmp)
          allocate( tmp(size(mat2%loc_mat,1),length) )
    
          ! part of second matrix broadcasted to all tasks
          if(my_id == i-1)  tmp    = mat2%loc_mat  
          if(my_id == i-1 ) length = size(tmp)
          call MPI_BCAST(length, 1, MPI_INTEGER, i-1, MPI_COMM_WORLD, ierr)
          call MPI_BCAST( tmp, length, MPI_DOUBLE_PRECISION, i-1, MPI_COMM_WORLD, ierr)
       
          ! Matrix -matrix multiplications
          glob_index_j=(i-1)*mat2%step
          do z = 1, size(mat1%loc_mat,1)
            do k = 1, size(tmp,2)
    
              sum_element=0.0 
!$omp parallel private (j) shared (sum_element,mat1,z,k,tmp) 
!$omp do reduction (+:sum_element)
              do j = 1, size(mat1%loc_mat,2)
                sum_element = sum_element + mat1%loc_mat(z,j) * tmp(j,k)         
              end do
!$omp end do 
!$omp end parallel
              
              res_mat%loc_mat(z,k+glob_index_j) = sum_element  
        
            end do
          end do
    
        end do     
      
      else
        write(*,*) 'ERROR in matrix_multiplication: Combination of rowwise/columnwise not implemented (1).'
        stop
      end if 
    
    ! --- Multiplication mat1(distributed)*mat2(distributed) =res_mat(not_distributed)
    else if ( present(mat2) .and. present(res_mat_not_distr) ) then
      
      if (allocated(res_mat_not_distr)) deallocate(res_mat_not_distr)
      allocate( res_mat_not_distr(mat1%dim(1), mat2%dim(2)) )  

      if ( mat1%row_wise .and. (.not. mat2%row_wise) ) then
       
        if(mat1%dim(2) /= mat2%dim(1) ) then
          write(*,*) "ERROR: Dimensions of  multiplied matrices (2) #col =", mat1%dim(2), &
                     " #row =", mat2%dim(1) , "STOP (2)"
          stop
        endif
    
        res_mat_not_distr=0.0   
    
        ! loop over all tasks
        do i = 1, ntasks
    
          ! Broadcasting size of the distributed column
          length = size(mat2%loc_mat,2)
        
          call MPI_BCAST(length, 1, MPI_INTEGER, i-1, MPI_COMM_WORLD, ierr )
    
          ! temporal matrix for receiving part of the distributed matrix and
          ! multiply it
          if ( allocated(tmp) ) deallocate(tmp)
          allocate( tmp(size(mat2%loc_mat,1),length) )
    
          ! part of second matrix broadcasted to all tasks
          if(my_id == i-1)  tmp    = mat2%loc_mat
          if(my_id == i-1 ) length = size(tmp)
          call MPI_BCAST(length, 1, MPI_INTEGER, i-1, MPI_COMM_WORLD, ierr)
          call MPI_BCAST(tmp, length, MPI_DOUBLE_PRECISION, i-1, MPI_COMM_WORLD,ierr)
    
          ! Matrix -matrix multiplications
    
          glob_index_j = (i-1)*mat2%step
          glob_index_i = my_id*mat2%step
          do z = 1, size(mat1%loc_mat,1)
            do k = 1, size(tmp,2)
    
              sum_element=0.0

!$omp parallel private (j) shared (sum_element,mat1,z,k,tmp) 
!$omp do reduction (+:sum_element)
              do j = 1, size(mat1%loc_mat,2)
                sum_element = sum_element + mat1%loc_mat(z,j) * tmp(j,k)
              end do
!$omp end do 
!$omp end parallel

              res_mat_not_distr(z + glob_index_i, k + glob_index_j) = sum_element
    
            end do
          end do
    
        end do
      else
        write(*,*) 'ERROR in matrix_multiplication: Combination of rowwise/columnwise not implemented (2).'
        stop
      end if
     
    ! --- Multiplication mat1(distributed)*mat2(not_distributed) =res_mat(not_distributed)
    else if(present(mat2_not_distr) .AND. present(res_mat_not_distr)) then

       if (allocated(res_mat_not_distr)) deallocate(res_mat_not_distr)
       allocate( res_mat_not_distr(mat1%dim(1), size(mat2_not_distr,2)) )

       if(mat1%dim(2) /= size(mat2_not_distr,1) ) then
          write(*,*) "ERROR: Dimensions of  multiplied matrices (3) #col =", mat1%dim(2), &
                     " #row =", size(mat2_not_distr,1) , "STOP (3)"
          stop
       endif
       
      !Check distribution. Result matrix not distributed, first matrix row wise, second matrix not distributed
      if (mat1%row_wise) then
    
        res_mat_not_distr=0.0
        
        ! Matrix -matrix multiplications     
        glob_index_i = my_id*mat1%step
    
        do z = 1, size(mat1%loc_mat,1)
          do k = 1, size(mat2_not_distr,2)
    
            sum_element=0.0

!$omp parallel private (j) shared (sum_element,mat1,z,k,tmp,mat2_not_distr) 
!$omp do reduction (+:sum_element) 
            do j = 1, size(mat1%loc_mat,2)
              sum_element = sum_element + mat1%loc_mat(z,j) *mat2_not_distr(j,k)
            end do
!$omp end do 
!$omp end parallel

            res_mat_not_distr(z + glob_index_i, k) = sum_element
    
          end do
        end do
        
      else
        write(*,*) 'ERROR in matrix_multiplication: Combination of rowwise/columnwise not implemented (3).'
        stop
      end if
    
    else
 
      write(*,*) 'ERROR in matrix_multiplication: Unsupported set of matrices provided.'
      stop

    end if
  
  end subroutine matrix_multiplication
  
  
  
  !> Broadcast the STARWALL response matrices to the other MPI procs.
  subroutine broadcast_starwall_response(my_id, sr)
    
    use mpi_mod
    
    implicit none
    
    
    ! --- Routine parameters
    integer,                   intent(in)    :: my_id
    type(t_starwall_response), intent(inout) :: sr
    
    ! --- Local parameters
    integer :: ierr
    real*8    :: test_sum, loc_sum    

    if ( vacuum_debug .and. (my_id == 0) ) write(*,*) my_id, 'Entering broadcast_starwall_response.'
    
    ! --- Broadcast parameters.
    call MPI_bcast(sr%file_version,            1, MPI_INTEGER,  0, MPI_COMM_WORLD, ierr)
    call MPI_bcast(sr%n_bnd,                   1, MPI_INTEGER,  0, MPI_COMM_WORLD, ierr)
    call MPI_bcast(sr%nd_bez,                  1, MPI_INTEGER,  0, MPI_COMM_WORLD, ierr)
    call MPI_bcast(sr%ncoil,                   1, MPI_INTEGER,  0, MPI_COMM_WORLD, ierr)
    call MPI_bcast(sr%npot_w,                  1, MPI_INTEGER,  0, MPI_COMM_WORLD, ierr)
    call MPI_bcast(sr%n_w,                     1, MPI_INTEGER,  0, MPI_COMM_WORLD, ierr)
    call MPI_bcast(sr%ntri_w,                  1, MPI_INTEGER,  0, MPI_COMM_WORLD, ierr)
    call MPI_bcast(sr%n_tor,                   1, MPI_INTEGER,  0, MPI_COMM_WORLD, ierr)
    call MPI_bcast(sr%n_tor0,                  1, MPI_INTEGER,  0, MPI_COMM_WORLD, ierr)

    call MPI_bcast(sr%nv,                      1, MPI_INTEGER,  0, MPI_COMM_WORLD, ierr)
    call MPI_bcast(sr%n_points,                1, MPI_INTEGER,  0, MPI_COMM_WORLD, ierr)
    call MPI_bcast(sr%iwall,                   1, MPI_INTEGER,  0, MPI_COMM_WORLD, ierr)
    call MPI_bcast(sr%nwu,                     1, MPI_INTEGER,  0, MPI_COMM_WORLD, ierr)
    call MPI_bcast(sr%nwv,                     1, MPI_INTEGER,  0, MPI_COMM_WORLD, ierr)
    call MPI_bcast(sr%mn_w,                    1, MPI_INTEGER,  0, MPI_COMM_WORLD, ierr)
    call MPI_bcast(sr%max_mn_w,                1, MPI_INTEGER,  0, MPI_COMM_WORLD, ierr)

    call MPI_bcast(sr%ntri_c,                  1, MPI_INTEGER,  0, MPI_COMM_WORLD, ierr)
    call MPI_bcast(sr%n_pol_coils,             1, MPI_INTEGER,  0, MPI_COMM_WORLD, ierr)
    call MPI_bcast(sr%n_rmp_coils,             1, MPI_INTEGER,  0, MPI_COMM_WORLD, ierr)
    call MPI_bcast(sr%n_voltage_coils,         1, MPI_INTEGER,  0, MPI_COMM_WORLD, ierr)
    call MPI_bcast(sr%n_diag_coils,            1, MPI_INTEGER,  0, MPI_COMM_WORLD, ierr)
    call MPI_bcast(sr%ind_start_pol_coils,     1, MPI_INTEGER,  0, MPI_COMM_WORLD, ierr)
    call MPI_bcast(sr%ind_start_rmp_coils,     1, MPI_INTEGER,  0, MPI_COMM_WORLD, ierr)
    call MPI_bcast(sr%ind_start_voltage_coils, 1, MPI_INTEGER,  0, MPI_COMM_WORLD, ierr)
    call MPI_bcast(sr%ind_start_diag_coils,    1, MPI_INTEGER,  0, MPI_COMM_WORLD, ierr)
    call MPI_bcast(sr%eta_thin_w,              1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    
    n_dof_starwall = sr%nd_bez
    n_wall_curr    = sr%n_w
    
    ! --- Allocate matrices.
    if ( my_id /= 0 ) then
      if (allocated(sr%i_tor)   ) deallocate(sr%i_tor);    allocate(sr%i_tor(sr%n_tor))
      if (allocated(sr%d_yy)    ) deallocate(sr%d_yy);     allocate(sr%d_yy(sr%n_w))
      if (allocated(sr%xyzpot_w)) deallocate(sr%xyzpot_w); allocate(sr%xyzpot_w(sr%npot_w,3))
      if (allocated(sr%jpot_w)  ) deallocate(sr%jpot_w);   allocate(sr%jpot_w(sr%ntri_w,3))

      if ( sr%file_version>=5 ) then 
        if (allocated(sr%m_w)        ) deallocate(sr%m_w);          allocate(        sr%m_w(sr%max_mn_w))
        if (allocated(sr%n_w_fourier)) deallocate(sr%n_w_fourier);  allocate(sr%n_w_fourier(sr%max_mn_w))

        if (allocated(sr%rc_w)       ) deallocate(sr%rc_w);         allocate(sr%rc_w(sr%max_mn_w))
        if (allocated(sr%rs_w)       ) deallocate(sr%rs_w);         allocate(sr%rs_w(sr%max_mn_w))
        if (allocated(sr%zc_w)       ) deallocate(sr%zc_w);         allocate(sr%zc_w(sr%max_mn_w))
        if (allocated(sr%zs_w)       ) deallocate(sr%zs_w);         allocate(sr%zs_w(sr%max_mn_w))
        if (allocated(sr%phi0_w)     ) deallocate(sr%phi0_w);       allocate(sr%phi0_w(sr%ntri_w,3))
      endif

        
      if ( sr%ncoil > 0 ) then
        if (allocated(sr%jtri_c)       ) deallocate(sr%jtri_c);        allocate(sr%jtri_c(sr%ncoil))
        if (allocated(sr%x_coil)       ) deallocate(sr%x_coil);        allocate(sr%x_coil(sr%ntri_c,3))
        if (allocated(sr%y_coil)       ) deallocate(sr%y_coil);        allocate(sr%y_coil(sr%ntri_c,3))
        if (allocated(sr%z_coil)       ) deallocate(sr%z_coil);        allocate(sr%z_coil(sr%ntri_c,3))
        if (allocated(sr%phi_coil)     ) deallocate(sr%phi_coil);      allocate(sr%phi_coil(sr%ntri_c,3))
        if (allocated(sr%eta_thin_coil)) deallocate(sr%eta_thin_coil); allocate(sr%eta_thin_coil(sr%ntri_c))
        if (allocated(sr%coil_resist)  ) deallocate(sr%coil_resist);   allocate(sr%coil_resist(sr%ncoil))
        if (allocated(sr%coil_name)  )   deallocate(sr%coil_name);     allocate(sr%coil_name(sr%ncoil))
      end if
    end if
    
    ! --- Broadcast matrices.
    call MPI_bcast(sr%i_tor,    sr%n_tor,                MPI_INTEGER,          0, MPI_COMM_WORLD, ierr)
    call MPI_bcast(sr%d_yy,     sr%n_w,                  MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_bcast(sr%xyzpot_w, sr%npot_w*3,             MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_bcast(sr%jpot_w,   sr%ntri_w*3,             MPI_INTEGER,          0, MPI_COMM_WORLD, ierr)

    if ( sr%file_version>=5 ) then 
      call MPI_bcast(sr%m_w,            sr%max_mn_w,     MPI_INTEGER,          0, MPI_COMM_WORLD, ierr)
      call MPI_bcast(sr%n_w_fourier,    sr%max_mn_w,     MPI_INTEGER,          0, MPI_COMM_WORLD, ierr)
      call MPI_bcast(sr%rc_w,           sr%max_mn_w,     MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call MPI_bcast(sr%rs_w,           sr%max_mn_w,     MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call MPI_bcast(sr%zc_w,           sr%max_mn_w,     MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call MPI_bcast(sr%zs_w,           sr%max_mn_w,     MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call MPI_bcast(sr%phi0_w,         sr%ntri_w*3,     MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    endif

    if ( sr%ncoil > 0 ) then
      call MPI_bcast(sr%jtri_c,   sr%ncoil,              MPI_INTEGER,          0, MPI_COMM_WORLD, ierr)
      call MPI_bcast(sr%x_coil,   sr%ntri_c*3,           MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call MPI_bcast(sr%y_coil,   sr%ntri_c*3,           MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call MPI_bcast(sr%z_coil,   sr%ntri_c*3,           MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call MPI_bcast(sr%phi_coil,      sr%ntri_c*3,      MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call MPI_bcast(sr%eta_thin_coil, sr%ntri_c,        MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call MPI_bcast(sr%coil_resist,   sr%ncoil,         MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call MPI_bcast(sr%coil_name,sr%ncoil*COIL_NAME_LEN,MPI_CHARACTER       , 0, MPI_COMM_WORLD, ierr)
    end if
    
    if ( vacuum_debug ) then

      loc_sum = sum(sr%a_ye%loc_mat) + sum(sr%a_ey%loc_mat) + sum(sr%a_ee%loc_mat) + &
                sum(sr%a_nw%loc_mat) + sum(sr%s_ww%loc_mat) + sum(sr%s_ww_inv%loc_mat)   

      if (my_id==0) then

        loc_sum =   loc_sum + sum(sr%i_tor) + sum(sr%d_yy)+ sum(sr%xyzpot_w)     &
                  + sum(sr%jpot_w) + sr%n_bnd + sr%nd_bez + sr%ncoil + sr%npot_w &
                  + sr%n_w + sr%ntri_w + sr%n_tor + sr%eta_thin_w                &
                  + sr%file_version + sr%ntri_c + sr%n_pol_coils                 &
                  + sr%n_rmp_coils + sr%n_voltage_coils + sr%n_diag_coils        &
                  + sr%ind_start_pol_coils + sr%ind_start_rmp_coils              &
                  + sr%ind_start_voltage_coils + sr%ind_start_diag_coils
      end if
  
      test_sum=0.0
      call MPI_Reduce(loc_sum, test_sum, 1, MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
 
      if(my_id == 0)  write(*,'("Checksum",I4,ES24.16)') my_id, test_sum
      
      if(my_id == 0)  write(*,*) my_id, 'Exiting broadcast_starwall_response.'
   
    end if  
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  end subroutine broadcast_starwall_response
  
  

  !> Write out information about the STARWALL response matrices.
  subroutine log_starwall_response(my_id, sr)

    use mod_parameters, only: n_period
    use mpi_mod

    implicit none

    integer,                   intent(in) :: my_id
    type(t_starwall_response), intent(in) :: sr

    real*8  :: test_sum
    integer :: ierr

    32 format(3x,77('-'))
    33 format(3x,a,i8)
    34 format(3x,'sum(',a,')=',es24.16)
    35 format(3x,'sum(',a,')=',i24)
    36 format(3x,'sum(',a,')= ---not allocated---')
    37 format(3x,a,es25.15)
  
    if (my_id == 0) then
         write(*,*)
       write(*,32)
       write(*,33) 'STARWALL RESPONSE INFORMATION:'
       write(*,32)
       write(*,33) 'file_version            =', sr%file_version
       write(*,33) 'n_bnd                   =', sr%n_bnd
       write(*,33) 'nd_bez                  =', sr%nd_bez
       write(*,33) 'nv                      =', sr%nv    
       write(*,33) 'n_points                =', sr%n_points
       write(*,33) 'ncoil                   =', sr%ncoil
       write(*,33) 'npot_w                  =', sr%npot_w
       write(*,33) 'n_w                     =', sr%n_w
       write(*,33) 'ntri_w                  =', sr%ntri_w
       write(*,33) 'iwall                   =', sr%iwall  
       write(*,33) 'nwu                     =', sr%nwu    
       write(*,33) 'nwv                     =', sr%nwv    
       write(*,33) 'mn_w                    =', sr%mn_w   
       write(*,33) 'max_mn_w                =', sr%max_mn_w   
       write(*,33) 'n_tor                   =', sr%n_tor
       write(*,33) 'n_tor0                  =', sr%n_tor0
       write(*,33) 'n_pol_coils             =', sr%n_pol_coils
       write(*,33) 'n_rmp_coils             =', sr%n_rmp_coils
       write(*,33) 'n_voltage_coils         =', sr%n_voltage_coils
       write(*,33) 'n_diag_coils            =', sr%n_diag_coils
       write(*,33) 'ind_start_pol_coils     =', sr%ind_start_pol_coils
       write(*,33) 'ind_start_rmp_coils     =', sr%ind_start_rmp_coils
       write(*,33) 'ind_start_voltage_coils =', sr%ind_start_voltage_coils
       write(*,33) 'ind_start_diag_coils    =', sr%ind_start_diag_coils
      if ( sr%file_version >= 2) write(*,37) 'eta_thin_w        =', sr%eta_thin_w
      if (allocated(sr%i_tor)) write(*,33) 'i_tor               ='//trim(modes_to_str(sr%i_tor,sr%n_tor,n_period))
    end if 

    if ( vacuum_debug ) then
      
      if (my_id == 0) then
        write(*,32)
        if (allocated(sr%i_tor        )) then; write(*,35) 'i_tor        ',sum(sr%i_tor         ); else; write(*,36) 'i_tor        '; end if
        if (allocated(sr%d_yy         )) then; write(*,34) 'd_yy         ',sum(sr%d_yy          ); else; write(*,36) 'd_yy         '; end if
      end if 

      test_sum=0.0
      if (allocated(sr%a_ye%loc_mat)) then
        call MPI_Reduce(sum(sr%a_ye%loc_mat), test_sum, 1, MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        if (my_id == 0) write(*,34) 'a_ye         ', test_sum 
      else
        if (my_id == 0)  write(*,36) 'a_ye         '
      end if            
 
      test_sum=0.0
      if (allocated(sr%a_ey%loc_mat)) then 
        call MPI_Reduce(sum(sr%a_ey%loc_mat), test_sum, 1, MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        if (my_id == 0) write(*,34) 'a_ey         ', test_sum 
      else
        if (my_id == 0)  write(*,36) 'a_ey         ' 
      end if

      test_sum=0.0
      if (allocated(sr%a_ee%loc_mat)) then
        call MPI_Reduce(sum(sr%a_ee%loc_mat), test_sum, 1, MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        if (my_id == 0) write(*,34) 'a_ee         ', test_sum
      else
        if (my_id == 0) write(*,36) 'a_ee         '
      end if

      test_sum=0.0
      if (allocated(sr%a_id%loc_mat)) then
        call MPI_Reduce(sum(sr%a_id%loc_mat), test_sum, 1, MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        if (my_id == 0) write(*,34) 'a_id         ', test_sum
      else
        if (my_id == 0) write(*,36) 'a_id         '
      end if

      test_sum=0.0
      if (allocated(sr%a_nw%loc_mat)) then
        call MPI_Reduce(sum(sr%a_nw%loc_mat), test_sum, 1, MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        if (my_id == 0) write(*,34) 'a_nw         ', test_sum
      else
        if (my_id == 0) write(*,36) 'a_nw         '
      end if
      
      test_sum=0.0
      if (allocated(sr%s_ww%loc_mat)) then
        call MPI_Reduce(sum(sr%s_ww%loc_mat), test_sum, 1, MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        if (my_id == 0) write(*,34) 's_ww         ', test_sum
      else
        if (my_id == 0) write(*,36) 's_nw         '
      end if

      test_sum=0.0
      if (allocated(sr%s_ww_inv%loc_mat)) then
        call MPI_Reduce(sum(sr%s_ww_inv%loc_mat), test_sum, 1, MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        if (my_id == 0) write(*,34) 's_ww_inv     ', test_sum
      else
        if (my_id == 0) write(*,36) 's_nw_inv         '
      end if
        
      if (my_id == 0) then
        if (allocated(sr%xyzpot_w     )) then; write(*,34) 'xyzpot_w     ',sum(sr%xyzpot_w      ); else; write(*,36) 'xyzpot_w     '; end if
        if (allocated(sr%jpot_w       )) then; write(*,35) 'jpot_w       ',sum(sr%jpot_w        ); else; write(*,36) 'jpot_w       '; end if
        if (allocated(sr%jtri_c       )) then; write(*,35) 'jtri_c       ',sum(sr%jtri_c        ); else; write(*,36) 'jtri_c       '; end if
        if (allocated(sr%x_coil       )) then; write(*,34) 'x_coil       ',sum(sr%x_coil        ); else; write(*,36) 'x_coil       '; end if
        if (allocated(sr%y_coil       )) then; write(*,34) 'y_coil       ',sum(sr%y_coil        ); else; write(*,36) 'y_coil       '; end if
        if (allocated(sr%z_coil       )) then; write(*,34) 'z_coil       ',sum(sr%z_coil        ); else; write(*,36) 'z_coil       '; end if
        if (allocated(sr%phi_coil     )) then; write(*,34) 'phi_coil     ',sum(sr%phi_coil      ); else; write(*,36) 'phi_coil     '; end if
        if (allocated(sr%eta_thin_coil)) then; write(*,34) 'eta_thin_coil',sum(sr%eta_thin_coil ); else; write(*,36) 'eta_thin_coil'; end if
        if (allocated(sr%coil_resist  )) then; write(*,34) 'coil_resist  ',sum(sr%coil_resist   ); else; write(*,36) 'coil_resist  '; end if
      end if
    end if

    if (my_id==0) then
      write(*,32)
      write(*,*)
    end if
   
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    
  end subroutine log_starwall_response
  
  
  
  !> Write out resistive-wall data as a VTK-file.
  !!
  !! Scalar quantities (at wall triangle nodes):
  !! * pot_w: Wall current potentials
  !! * dpot_w: Change of wall current potentials in previous time-step
  !!
  !! Vector quantities (at triangles):
  !! * jsurf_w: Surface currents on the wall
  subroutine write_wall_vtk(index, resistive_wall, my_id)
    
    use phys_module, only: nout
    use constants
    use mpi_mod    

    implicit none
    
    ! --- Routine parameters
    integer, intent(in) :: index !< Time step index
    logical, intent(in) :: resistive_wall
    integer, intent(in) :: my_id
    
    ! --- Local variables
    real*8              :: phi1, phi2, phi3, r1(3), r2(3), r3(3), r21(3), r32(3), r13(3), r21_cross_r32(3)
    real*8              :: nm(3), n13(3), n32(3), n21(3), j_lin(3), j13, j32, j21
    real*8              :: ephi12(3),ephi13(3),ephi23(3), jsides(3)
    real*8              :: pol13(3), pol32(3), pol21(3)
    real*8              :: mid12(3), mid23(3), mid13(3), angle12, angle13, angle23
    real*8              :: Iw_net_tor 
    integer             :: filehandle = 60, i, maxcurr_pos
    logical             :: Iphi_max, Ipol_max, jphi_lin, jpol_lin
    character(len=18)   :: filename
    real*8, allocatable :: tripot_w(:), dtripot_w(:)
    integer :: ierr

    if ( mod(index,nout) /= 0 ) return

    if(my_id == 0) then
      if ( mod(index,nout) /= 0 ) return
      
      ! --- Preset values
      Iphi_max = .false.
      jphi_lin = .false.
      Ipol_max = .false.
      jpol_lin = .false.    
      
      ! --- VTK file header
      write(filename,'(A,I5.5,A)') 'wallcurr.',index,'.vtk'
      open(filehandle, file=filename, status='replace', action='write')
      140 format(a)
      141 format(a,i8,a)
      142 format(3es16.8)
      143 format(a,2i8)
      144 format(4i8)
      write(filehandle,140) '# vtk DataFile Version 2.0'
      write(filehandle,140) 'testdata'
      write(filehandle,140) 'ASCII'
      write(filehandle,140) 'DATASET POLYDATA'
      
      ! --- Triangle node positions
      write(filehandle,141) 'POINTS', sr%npot_w, ' float'
      do i = 1, sr%npot_w
        write(filehandle,142) sr%xyzpot_w(i,:)
      end do
      
      ! --- Node indices corresponding to triangles
      write(filehandle,143) 'POLYGONS', sr%ntri_w, sr%ntri_w * 4
      do i = 1, sr%ntri_w
        write(filehandle,144) 3, sr%jpot_w(i,:) - 1
      end do
      
      ! --- Wall current potentials
      write(filehandle,141) 'POINT_DATA', sr%npot_w
      write(filehandle,140) 'SCALARS pot_w float'
      write(filehandle,140) 'LOOKUP_TABLE default'
    end if
    
    call reconstruct_triangle_potentials(tripot_w, wall_curr, my_id, Iw_net_tor)

    if(my_id == 0) then
      do i = 1, sr%npot_w
        write(filehandle,142) tripot_w(i)
      end do
      
      ! --- Change of wall current potentials in previous time-step
      write(filehandle,140) 'SCALARS dpot_w float'
      write(filehandle,140) 'LOOKUP_TABLE default'
    end if

    call reconstruct_triangle_potentials(dtripot_w, dwall_curr, my_id)

    if(my_id == 0) then
  
      do i = 1, sr%npot_w
        write(filehandle,142) dtripot_w(i)
      end do

      deallocate(dtripot_w)
      
      ! --- Cell data variables
      write(filehandle,141) 'CELL_DATA', sr%ntri_w
     
      ! --- Maximum toroidal current flwoing on a triangle (kA)
      if (Iphi_max) then
          
        write(filehandle,140) 'SCALARS Iphi_max(kA) float'
        write(filehandle,140) 'LOOKUP_TABLE default'
         
        do i = 1, sr%ntri_w
         
          ! --- Wall potentials at triangle nodes
          phi1   = tripot_w(sr%jpot_w(i,1))
          phi2   = tripot_w(sr%jpot_w(i,2))
          phi3   = tripot_w(sr%jpot_w(i,3))
          
          ! --- Positions of triangle nodes
          r1(:)  = sr%xyzpot_w(sr%jpot_w(i,1),:)
          r2(:)  = sr%xyzpot_w(sr%jpot_w(i,2),:)
          r3(:)  = sr%xyzpot_w(sr%jpot_w(i,3),:)
          
          ! --- Middle points on triangle sides
          mid12(:) =   (r1(:) + r2(:)) * 0.5d0
          mid13(:) =   (r1(:) + r3(:)) * 0.5d0
          mid23(:) =   (r2(:) + r3(:)) * 0.5d0      
               
          ! --- Toroidal angles on middle points, JOREK coordinate system  y --> 1, x --> 3
          angle12 = atan2(-mid12(1),mid12(3))
          angle13 = atan2(-mid13(1),mid13(3))
          angle23 = atan2(-mid23(1),mid23(3))
          
          ! --- Toroidal basis vectors on middle points, in clock-wise direction looking from above the torus
          ephi12(:) = (/- cos(angle12), 0.d0, - sin(angle12) /)
          ephi13(:) = (/- cos(angle13), 0.d0, - sin(angle13) /)
          ephi23(:) = (/- cos(angle23), 0.d0, - sin(angle23) /)
         
          ! --- Quantities needed to calculate the current density vector 
          r21(:) = r1(:)-r2(:)
          r32(:) = r2(:)-r3(:)
          r13(:) = r3(:)-r1(:)
         
          r21_cross_r32(:) = (/ r21(2)*r32(3) - r21(3)*r32(2), r21(3)*r32(1) - r21(1)*r32(3),          &
            r21(1)*r32(2) - r21(2)*r32(1) /)
            
          !--- current density vector in kA/m, Merkel 2015
          j_lin = ( phi1*(r3-r2)+phi2*(r1-r3)+phi3*(r2-r1) ) / sqrt(sum(r21_cross_r32**2))
          j_lin = j_lin / mu_zero * 1.d-3  
            
          !--- Vector normal to triangle surface
          nm(:)  = r21_cross_r32(:) / sqrt(sum(r21_cross_r32**2))
          
          !--- Normal vectors to triangle sides
          n13(:) = (/ r13(2)*nm(3) - r13(3)*nm(2), r13(3)*nm(1) - r13(1)*nm(3),          &
            r13(1)*nm(2) - r13(2)*nm(1) /)
          n21(:) = (/ r21(2)*nm(3) - r21(3)*nm(2), r21(3)*nm(1) - r21(1)*nm(3),          &
            r21(1)*nm(2) - r21(2)*nm(1) /)
          n32(:) = (/ r32(2)*nm(3) - r32(3)*nm(2), r32(3)*nm(1) - r32(1)*nm(3),          &
            r32(1)*nm(2) - r32(2)*nm(1) /)
          
          !--- Normalized vectors normal to triangle sides
          n13(:) = n13(:)/ sqrt(sum(n13**2))
          n21(:) = n21(:)/ sqrt(sum(n21**2))
          n32(:) = n32(:)/ sqrt(sum(n32**2))
          
          !--- Toroidal current flowing per side of the triangle 
          j13    = abs(dot_product(n13,ephi13)) * sqrt(sum(r13**2)) * dot_product(ephi13,j_lin)
          j21    = abs(dot_product(n21,ephi12)) * sqrt(sum(r21**2)) * dot_product(ephi12,j_lin)
          j32    = abs(dot_product(n32,ephi23)) * sqrt(sum(r32**2)) * dot_product(ephi23,j_lin)
          
          jsides(:) = (/ j13, j21, j32 /)
              
          !--- Find the triangle side with maximum absolute current
          maxcurr_pos   = maxloc(abs(jsides),1)
          
          !--- Maximum toroidal current flowing on a triangle side, in kA
          write(filehandle,142) jsides(maxcurr_pos) 
        
        end do
      end if
      
      ! --- Maximum poloidal current flwoing on a triangle (kA)
      if (Ipol_max) then
     
        write(filehandle,140) 'SCALARS Ipol_max(kA) float'
        write(filehandle,140) 'LOOKUP_TABLE default'
         
        do i = 1, sr%ntri_w
         
          ! --- Wall potentials at triangle nodes
          phi1   = tripot_w(sr%jpot_w(i,1))
          phi2   = tripot_w(sr%jpot_w(i,2))
          phi3   = tripot_w(sr%jpot_w(i,3))
          
          ! --- Positions of triangle nodes
          r1(:)  = sr%xyzpot_w(sr%jpot_w(i,1),:)
          r2(:)  = sr%xyzpot_w(sr%jpot_w(i,2),:)
          r3(:)  = sr%xyzpot_w(sr%jpot_w(i,3),:)
          
          ! --- Middle points on triangle sides
          mid12(:) =   (r1(:) + r2(:)) * 0.5d0
          mid13(:) =   (r1(:) + r3(:)) * 0.5d0
          mid23(:) =   (r2(:) + r3(:)) * 0.5d0      
               
          ! --- Toroidal angles on middle points, JOREK coordinate system  y --> 1, x --> 3
          angle12 = atan2(-mid12(1),mid12(3))
          angle13 = atan2(-mid13(1),mid13(3))
          angle23 = atan2(-mid23(1),mid23(3))
          
          ! --- Toroidal basis vectors on middle points, in clock-wise direction looking from above the torus
          ephi12(:) = (/- cos(angle12), 0.d0, - sin(angle12) /)
          ephi13(:) = (/- cos(angle13), 0.d0, - sin(angle13) /)
          ephi23(:) = (/- cos(angle23), 0.d0, - sin(angle23) /)
         
          ! --- Quantities needed to calculate the current density vector 
          r21(:) = r1(:)-r2(:)
          r32(:) = r2(:)-r3(:)
          r13(:) = r3(:)-r1(:)
         
          r21_cross_r32(:) = (/ r21(2)*r32(3) - r21(3)*r32(2), r21(3)*r32(1) - r21(1)*r32(3),          &
            r21(1)*r32(2) - r21(2)*r32(1) /)
            
          !--- current density vector in kA/m, Merkel 2015
          j_lin = ( phi1*(r3-r2)+phi2*(r1-r3)+phi3*(r2-r1) ) / sqrt(sum(r21_cross_r32**2))
          j_lin = j_lin / mu_zero * 1.d-3  
            
          !--- Vector normal to triangle surface
          nm(:)  = r21_cross_r32(:) / sqrt(sum(r21_cross_r32**2))
          
          !--- Poloidal vectors by triangle side (e_phi x n), (at middle points)
          pol13(:) = (/ ephi13(2)*nm(3) - ephi13(3)*nm(2), ephi13(3)*nm(1) - ephi13(1)*nm(3),          &
            ephi13(1)*nm(2) - ephi13(2)*nm(1) /)
          pol21(:) = (/ ephi12(2)*nm(3) - ephi12(3)*nm(2), ephi12(3)*nm(1) - ephi12(1)*nm(3),          &
            ephi12(1)*nm(2) - ephi12(2)*nm(1) /)
          pol32(:) = (/ ephi23(2)*nm(3) - ephi23(3)*nm(2), ephi23(3)*nm(1) - ephi23(1)*nm(3),          &
            ephi23(1)*nm(2) - ephi23(2)*nm(1) /)
          
          !--- Normalized poloidal vectors by triangle side
          pol13(:) = pol13(:)/ sqrt(sum(pol13**2))
          pol21(:) = pol21(:)/ sqrt(sum(pol21**2))
          pol32(:) = pol32(:)/ sqrt(sum(pol32**2))
          
          !--- Normal vectors to triangle sides
          n13(:) = (/ r13(2)*nm(3) - r13(3)*nm(2), r13(3)*nm(1) - r13(1)*nm(3),          &
            r13(1)*nm(2) - r13(2)*nm(1) /)
          n21(:) = (/ r21(2)*nm(3) - r21(3)*nm(2), r21(3)*nm(1) - r21(1)*nm(3),          &
            r21(1)*nm(2) - r21(2)*nm(1) /)
          n32(:) = (/ r32(2)*nm(3) - r32(3)*nm(2), r32(3)*nm(1) - r32(1)*nm(3),          &
            r32(1)*nm(2) - r32(2)*nm(1) /)
          
          !--- Normalized vectors normal to triangle sides
          n13(:) = n13(:)/ sqrt(sum(n13**2))
          n21(:) = n21(:)/ sqrt(sum(n21**2))
          n32(:) = n32(:)/ sqrt(sum(n32**2))
          
          !--- Poloidal current flowing per side of the triangle 
          j13    = abs(dot_product(n13,pol13)) * sqrt(sum(r13**2)) * dot_product(pol13,j_lin)
          j21    = abs(dot_product(n21,pol21)) * sqrt(sum(r21**2)) * dot_product(pol21,j_lin)
          j32    = abs(dot_product(n32,pol32)) * sqrt(sum(r32**2)) * dot_product(pol32,j_lin)
          
          jsides(:) = (/ j13, j21, j32 /)
              
          !--- Find the triangle side with maximum absolute current
          maxcurr_pos   = maxloc(abs(jsides),1)
          
          !--- Maximum poloidal current flowing on a triangle side, in kA
          write(filehandle,142) jsides(maxcurr_pos) 
        
        end do
      end if
      
      ! --- Linear toroidal current density per triangle (kA/m)
      if (jphi_lin) then
     
        write(filehandle,140) 'SCALARS jphi_lin(kA/m) float'
        write(filehandle,140) 'LOOKUP_TABLE default'
         
        do i = 1, sr%ntri_w
         
          ! --- Wall potentials at triangle nodes
          phi1   = tripot_w(sr%jpot_w(i,1))
          phi2   = tripot_w(sr%jpot_w(i,2))
          phi3   = tripot_w(sr%jpot_w(i,3))
          
          ! --- Positions of triangle nodes
          r1(:)  = sr%xyzpot_w(sr%jpot_w(i,1),:)
          r2(:)  = sr%xyzpot_w(sr%jpot_w(i,2),:)
          r3(:)  = sr%xyzpot_w(sr%jpot_w(i,3),:)
          
          ! --- Middle points on triangle sides
          mid12(:) =   (r1(:) + r2(:)) * 0.5d0
          mid13(:) =   (r1(:) + r3(:)) * 0.5d0
          mid23(:) =   (r2(:) + r3(:)) * 0.5d0      
               
          ! --- Toroidal angles on middle points, JOREK coordinate system  y --> 1, x --> 3
          angle12 = atan2(-mid12(1),mid12(3))
          angle13 = atan2(-mid13(1),mid13(3))
          angle23 = atan2(-mid23(1),mid23(3))
          
          ! --- Toroidal basis vectors on middle points, in clock-wise direction looking from above the torus
          ephi12(:) = (/- cos(angle12), 0.d0, - sin(angle12) /)
          ephi13(:) = (/- cos(angle13), 0.d0, - sin(angle13) /)
          ephi23(:) = (/- cos(angle23), 0.d0, - sin(angle23) /)
         
          ! --- Quantities needed to calculate the current density vector 
          r21(:) = r1(:)-r2(:)
          r32(:) = r2(:)-r3(:)
          r13(:) = r3(:)-r1(:)
         
          r21_cross_r32(:) = (/ r21(2)*r32(3) - r21(3)*r32(2), r21(3)*r32(1) - r21(1)*r32(3),          &
            r21(1)*r32(2) - r21(2)*r32(1) /)
            
          !--- current density vector in kA/m, Merkel 2015
          j_lin = ( phi1*(r3-r2)+phi2*(r1-r3)+phi3*(r2-r1) ) / sqrt(sum(r21_cross_r32**2))
          j_lin = j_lin / mu_zero * 1.d-3  
          
          !--- Toroidal linear density current in each triangle side
          j13    =  dot_product(ephi13,j_lin)
          j21    =  dot_product(ephi12,j_lin)
          j32    =  dot_product(ephi23,j_lin)
          
          jsides(:) = (/ j13, j21, j32 /)
              
          !--- Find the triangle side with maximum current
          maxcurr_pos   = maxloc(abs(jsides),1)
          
          !--- Maximum toroidal linear density current flowing on a triangle side, in kA/m
          write(filehandle,142) jsides(maxcurr_pos) 
        
        end do
      end if
      
      ! --- Poloidal linear density current flowing on a triangle (kA/m)
      if (jpol_lin) then
     
        write(filehandle,140) 'SCALARS jpol_lin(kA/m) float'
        write(filehandle,140) 'LOOKUP_TABLE default'
         
        do i = 1, sr%ntri_w
         
          ! --- Wall potentials at triangle nodes
          phi1   = tripot_w(sr%jpot_w(i,1))
          phi2   = tripot_w(sr%jpot_w(i,2))
          phi3   = tripot_w(sr%jpot_w(i,3))
          
          ! --- Positions of triangle nodes
          r1(:)  = sr%xyzpot_w(sr%jpot_w(i,1),:)
          r2(:)  = sr%xyzpot_w(sr%jpot_w(i,2),:)
          r3(:)  = sr%xyzpot_w(sr%jpot_w(i,3),:)
          
          ! --- Middle points on triangle sides
          mid12(:) =   (r1(:) + r2(:)) * 0.5d0
          mid13(:) =   (r1(:) + r3(:)) * 0.5d0
          mid23(:) =   (r2(:) + r3(:)) * 0.5d0      
               
          ! --- Toroidal angles on middle points, JOREK coordinate system  y --> 1, x --> 3
          angle12 = atan2(-mid12(1),mid12(3))
          angle13 = atan2(-mid13(1),mid13(3))
          angle23 = atan2(-mid23(1),mid23(3))
          
          ! --- Toroidal basis vectors on middle points, in clock-wise direction looking from above the torus
          ephi12(:) = (/- cos(angle12), 0.d0, - sin(angle12) /)
          ephi13(:) = (/- cos(angle13), 0.d0, - sin(angle13) /)
          ephi23(:) = (/- cos(angle23), 0.d0, - sin(angle23) /)
         
          ! --- Quantities needed to calculate the current density vector 
          r21(:) = r1(:)-r2(:)
          r32(:) = r2(:)-r3(:)
          r13(:) = r3(:)-r1(:)
         
          r21_cross_r32(:) = (/ r21(2)*r32(3) - r21(3)*r32(2), r21(3)*r32(1) - r21(1)*r32(3),          &
            r21(1)*r32(2) - r21(2)*r32(1) /)
            
          !--- current density vector in kA/m, Merkel 2015
          j_lin = ( phi1*(r3-r2)+phi2*(r1-r3)+phi3*(r2-r1) ) / sqrt(sum(r21_cross_r32**2))
          j_lin = j_lin / mu_zero * 1.d-3  
            
          !--- Vector normal to triangle surface
          nm(:)  = r21_cross_r32(:) / sqrt(sum(r21_cross_r32**2))
          
          !--- Poloidal vectors by triangle side (e_phi x n), (at middle points)
          pol13(:) = (/ ephi13(2)*nm(3) - ephi13(3)*nm(2), ephi13(3)*nm(1) - ephi13(1)*nm(3),          &
            ephi13(1)*nm(2) - ephi13(2)*nm(1) /)
          pol21(:) = (/ ephi12(2)*nm(3) - ephi12(3)*nm(2), ephi12(3)*nm(1) - ephi12(1)*nm(3),          &
            ephi12(1)*nm(2) - ephi12(2)*nm(1) /)
          pol32(:) = (/ ephi23(2)*nm(3) - ephi23(3)*nm(2), ephi23(3)*nm(1) - ephi23(1)*nm(3),          &
            ephi23(1)*nm(2) - ephi23(2)*nm(1) /)
          
          !--- Normalized poloidal vectors by triangle side
          pol13(:) = pol13(:)/ sqrt(sum(pol13**2))
          pol21(:) = pol21(:)/ sqrt(sum(pol21**2))
          pol32(:) = pol32(:)/ sqrt(sum(pol32**2))
          
          !--- Poloidal current flowing per side of the triangle 
          j13    =  dot_product(pol13,j_lin)
          j21    =  dot_product(pol21,j_lin)
          j32    =  dot_product(pol32,j_lin)
          
          jsides(:) = (/ j13, j21, j32 /)
              
          !--- Find the triangle side with maximum absolute current
          maxcurr_pos   = maxloc(abs(jsides),1)
          
          !--- Maximum poloidal linear density current flowing on a triangle, in kA/m
          write(filehandle,142) jsides(maxcurr_pos) 
        
        end do
      end if
      
      ! --- Current density: Contribution of single valued potentials  
      write(filehandle,140) 'VECTORS jw_single_val(MA/m) float'

      do i = 1, sr%ntri_w
        ! --- Wall potential at triangle nodes
        phi1   = tripot_w(sr%jpot_w(i,1))
        phi2   = tripot_w(sr%jpot_w(i,2))
        phi3   = tripot_w(sr%jpot_w(i,3))
        ! --- Position of triangle nodes
        r1(:)  = sr%xyzpot_w(sr%jpot_w(i,1),:)
        r2(:)  = sr%xyzpot_w(sr%jpot_w(i,2),:)
        r3(:)  = sr%xyzpot_w(sr%jpot_w(i,3),:)
        r21(:) = r1(:)-r2(:)
        r32(:) = r2(:)-r3(:)
        r21_cross_r32(:) = (/ r21(2)*r32(3) - r21(3)*r32(2), r21(3)*r32(1) - r21(1)*r32(3),          &
          r21(1)*r32(2) - r21(2)*r32(1) /)
          
        ! Exports the linear wall density current in MA/m  
        write(filehandle,142) ( phi1*(r3-r2)+phi2*(r1-r3)+phi3*(r2-r1) ) / sqrt(sum(r21_cross_r32**2)) &
                              / mu_zero * 1.d-6
      end do

      if (sr%file_version >= 5) then

        ! --- Current density: Contribution of net wall potential 
        write(filehandle,140) 'VECTORS jw_net(MA/m) float'
  
        do i = 1, sr%ntri_w
          ! --- Wall potential at triangle nodes
          phi1   = Iw_net_tor * sr%phi0_w(i,1) 
          phi2   = Iw_net_tor * sr%phi0_w(i,2) 
          phi3   = Iw_net_tor * sr%phi0_w(i,3)
          ! --- Position of triangle nodes
          r1(:)  = sr%xyzpot_w(sr%jpot_w(i,1),:)
          r2(:)  = sr%xyzpot_w(sr%jpot_w(i,2),:)
          r3(:)  = sr%xyzpot_w(sr%jpot_w(i,3),:)
          r21(:) = r1(:)-r2(:)
          r32(:) = r2(:)-r3(:)
          r21_cross_r32(:) = (/ r21(2)*r32(3) - r21(3)*r32(2), r21(3)*r32(1) - r21(1)*r32(3),          &
            r21(1)*r32(2) - r21(2)*r32(1) /)
            
          ! Exports the linear wall density current in MA/m  
          write(filehandle,142) ( phi1*(r3-r2)+phi2*(r1-r3)+phi3*(r2-r1) ) / sqrt(sum(r21_cross_r32**2)) &
                                / mu_zero * 1.d-6
        end do

        ! --- Total wall current density
        write(filehandle,140) 'VECTORS jw_tot(MA/m) float'
  
        do i = 1, sr%ntri_w
          ! --- Wall potential at triangle nodes
          phi1   = Iw_net_tor * sr%phi0_w(i,1) + tripot_w(sr%jpot_w(i,1)) 
          phi2   = Iw_net_tor * sr%phi0_w(i,2) + tripot_w(sr%jpot_w(i,2))
          phi3   = Iw_net_tor * sr%phi0_w(i,3) + tripot_w(sr%jpot_w(i,3))
          ! --- Position of triangle nodes
          r1(:)  = sr%xyzpot_w(sr%jpot_w(i,1),:)
          r2(:)  = sr%xyzpot_w(sr%jpot_w(i,2),:)
          r3(:)  = sr%xyzpot_w(sr%jpot_w(i,3),:)
          r21(:) = r1(:)-r2(:)
          r32(:) = r2(:)-r3(:)
          r21_cross_r32(:) = (/ r21(2)*r32(3) - r21(3)*r32(2), r21(3)*r32(1) - r21(1)*r32(3),          &
            r21(1)*r32(2) - r21(2)*r32(1) /)
            
          ! Exports the linear wall density current in MA/m  
          write(filehandle,142) ( phi1*(r3-r2)+phi2*(r1-r3)+phi3*(r2-r1) ) / sqrt(sum(r21_cross_r32**2)) &
                                / mu_zero * 1.d-6
        end do
      endif    ! file version 5

      
      ! --- Close file, clean up
      close(filehandle)
  
    end if  ! my_id==0

    if ( allocated(tripot_w) ) deallocate( tripot_w )

  end subroutine write_wall_vtk
  
  

  !> Implements the boundary integral in the current equation which vanishes for fixed boundary
  !! conditions.
  !!
  !! The magnetic field parallel to the interface (boundary of JOREK computational domain) is
  !! expressed by the STARWALL vacuum response in terms of the poloidal magnetic field at the
  !! interface (ideal and resistive wall) and the wall currents (resistive wall).
  subroutine vacuum_boundary_integral(my_id, bnd_node_list, node_list,bnd_elm_list,               &
    freeboundary_equil, resistive_wall, index_min, index_max, rhs_loc, A_mat, tstep, index_now,   &
    irn, jcn, n_matrix_block_size, ijA_index, ijA_size, irn_jcn, i_tor_min, i_tor_max)

    use data_structure, only: type_node_list, type_bnd_node_list,type_bnd_element_list, type_bnd_element
    use mod_parameters,     only: n_plane, n_var, n_tor
    use gauss,          only: n_gauss, xgauss, wgauss
    use global_distributed_matrix, only: irn_glob, jcn_glob, a_glob, ndof_glob,det_row_col,       &
      det_sparse_pos, det_sparse_pos_block
    use basis_at_gaussian, only: H1, H1_s, HZ
    use phys_module, only: t_now, t_start
    use mpi_mod
    use mod_integer_types

    implicit none

    ! --- Routine parameters
    integer,                            intent(in)    :: my_id                !< MPI process ID
    type(type_node_list),               intent(in)    :: node_list            !< List of grid nodes
    type(type_bnd_node_list),           intent(in)    :: bnd_node_list        !< List of boundary grid nodes
    type(type_bnd_element_list),        intent(in)    :: bnd_elm_list         !< List of boundary elements
    logical,                            intent(in)    :: freeboundary_equil   !< Use free boundary equilibrium?
    logical,                            intent(in)    :: resistive_wall       !< Resistive or ideal wall?
    integer,                            intent(in)    :: index_min, index_max !< Responsibility of MPI proc
    real*8,                             intent(inout) :: rhs_loc(:)           !< Part of RHS of MPI proc 
    real*8,                allocatable, intent(inout) :: A_mat(:)             !< Distributed global or harmonic matrix 
    real*8,                             intent(in)    :: tstep                !< delta t, timestep
    integer,                            intent(in)    :: index_now            !< Current timestep index
    integer,                            intent(in)    :: i_tor_min, i_tor_max !< Toroidal harmonics 
    integer(kind=int_all), allocatable, intent(inout) :: irn(:), jcn(:)  
    integer(kind=int_all),              intent(in)    :: n_matrix_block_size
    integer(kind=int_all), allocatable, intent(in)    :: ijA_index(:,:)
    integer(kind=int_all), allocatable, intent(in)    :: ijA_size(:)
    integer(kind=int_all), allocatable, intent(in)    :: irn_jcn(:,:)

    ! --- Local variables
    real*8, allocatable :: psibnd_vec(:)    ! Vector of the values of Psi at the boundary
    real*8, allocatable :: dpsibnd_vec(:)   ! Vector of the values of deltaPsi at the boundary
    real*8, allocatable :: psibnd_coils(:)  ! Vector of the values of Psi_coil at the boundary
    real*8   :: amat_contrib, rhs_contrib   ! Vacuum response contribution to lhs and rhs
    real*8   :: testfunc_l                  ! j^*_l in documentation
    real*8   :: basfunc_i                   ! b_i in documentation
    real*8   :: dA                          ! factor from definition of dA
    real*8   :: x_s(n_gauss), y_s(n_gauss)  ! values of dR/ds and dZ/ds at Gaussian points
    real*8   :: common_prefactor
    integer  :: m_bndelem                   ! Boundary element index
    type(type_bnd_element) :: bndelem_m     ! Boundary element corresponding to index m_bndelem
    integer  :: ms                          ! Gauss point index
    integer  :: m_plane                     ! Toroidal plane index
    integer  :: sparsepos_jp, sparsepos_pp  ! Position of lhs contribution in the sparse matrix
    integer  :: blockpos_jp, blockpos_pp    ! Position of respective block in sparse matrix
    !   --- Test function related quantities
    integer  :: l_vertex, l_dof, l_dir, l_node, l_node_bnd, l_index, l_tor, l_row_j, l_row_psi
    real*8   :: l_size
    !   --- Quantities related to the boundary dof at which response is calculated
    integer  :: i_vertex, i_dof, i_dir, i_node, i_node_bnd, i_index, i_starwall, i_tor, i_resp, i_resp_0
    real*8   :: i_size
    !   --- Quantities related to the boundary dof contributing to the response
    integer  :: j_dof, j_dir, j_node, j_node_bnd, j_index, j_starwall, j_tor, j_col_psi, j_resp

    integer  :: i_resp_old

    integer  :: i_start_pf, i_end_pf     ! Indices for PF coils
   
#ifdef __GFORTRAN__
    real*8 :: wgauss_copy(4)
#endif
    real*8 :: t_elaps_start, t_elaps_end !### timing ###
    logical, save  :: PF_perturbation = .true.
    integer  :: ierr,i
    real*8, allocatable :: rhs_contrib_arr(:)

    if ( sr%n_tor == 0 ) then
      write(*,*) 'Skipping vacuum_boundary_integral since sr%n_tor==0.'
      return
    end if
    
    if ( vacuum_debug ) t_elaps_start = MPI_WTIME()  ! for timing

    if ( vacuum_debug ) write(*,*) my_id, 'Before:', sum(abs(rhs_loc)),sum(abs(A_mat))

    ! --- Determine vectors of the psi and deltapsi boundary values.
    call det_psibnd_vec(bnd_node_list, node_list, psibnd_vec, dpsibnd_vec, psibnd_coils)

    ! --- Calculate current source term for PF coil currents imposition
    if ( (sr%ncoil /= 0) ) then    !avoid it when using COIL_FIELD
      call coil_current_source(my_id)
    else
      if (.not. allocated (Y_coils0)) allocate(Y_coils0(n_wall_curr))
      Y_coils0(:) = 0.d0
    endif

    ! --- Update the derived response matrices
    call update_response(my_id, tstep, freeboundary_equil, resistive_wall)

    ! --- Perform the time-stepping for the wall currents.
    if ( resistive_wall .and.  (sr%n_tor>0) ) call evolve_wall_currents(my_id, psibnd_vec, dpsibnd_vec)

    ! --- Update old_dpsibnd_vec  (used for updating the wall currents in the next iteration)
    old_dpsibnd_vec(:) = dpsibnd_vec(:)

    if ( vacuum_debug ) then
      write(*,*) my_id, 'psibnd_vec:  ', sum(abs(psibnd_vec)), sum(psibnd_vec)
      write(*,*) my_id, 'dpsibnd_vec: ', sum(abs(dpsibnd_vec)), sum(dpsibnd_vec)
    end if

    call boundary_check(my_id)
    allocate(rhs_contrib_arr(n_dof_starwall))
    rhs_contrib_arr=0.0

    do i = response_m_f%ind_start,response_m_f%ind_end
      rhs_contrib_arr(i) =  sum(response_m_f%loc_mat(i-my_id*response_m_f%step, :) *  wall_curr(:) )  &
        +sum(response_m_g%loc_mat(i-my_id*response_m_g%step, :) * dwall_curr(:) )  &
        -sum(response_m_v%loc_mat(i-my_id*response_m_v%step, :) *   Y_coils0(:) )
    enddo
    call MPI_ALLReduce(MPI_IN_PLACE, rhs_contrib_arr, size(rhs_contrib_arr),MPI_DOUBLE_PRECISION, &
                       MPI_SUM, MPI_COMM_WORLD, ierr)

#ifdef __GFORTRAN__
    wgauss_copy(1:4) = wgauss(1:4)
#endif
    ! --- Sum over boundary elements
    !$omp parallel do                                                                                         &
    !$omp default(none)                                                                                       &    
    !$omp shared(my_id,A_mat, rhs_loc, bnd_elm_list, bnd_node_list, node_list, index_min, index_max,          &
    !$omp   response_m_e, response_m_f, response_m_g, response_m_h, response_m_j, H1, HZ, sr,                 &
    !$omp   i_tor_min, i_tor_max, irn, jcn, ijA_index, ijA_size, irn_jcn,                                     & 
    !$omp   bext_tan, I_coils, wall_curr, dwall_curr, psibnd_vec, dpsibnd_vec,psibnd_coils,                   &
#ifdef __GFORTRAN__
    !$omp   wgauss_copy,                                                                                      &
#endif
    !$omp   starwall_equil_coils, response_m_v, Y_coils0,rhs_contrib_arr, n_matrix_block_size,resistive_wall) &
    !$omp private(m_bndelem, bndelem_m, x_s, y_s, l_vertex, l_dof, l_node,l_dir, l_node_bnd,                  &
    !$omp   l_index, l_size, l_tor, l_row_j, l_row_psi, ms, dA, m_plane,common_prefactor,                     &
    !$omp   testfunc_l, i_vertex, i_dof, i_node, i_dir, i_node_bnd, i_index, i_size, i_starwall,              &
    !$omp   i_tor, i_resp, i_resp_old, i_resp_0, basfunc_i, j_node_bnd, j_dof, j_node, j_dir,                 &
    !$omp   j_index, j_starwall, j_tor, j_resp, j_col_psi, sparsepos_jp, sparsepos_pp,                        &
    !$omp   amat_contrib, rhs_contrib, blockpos_jp, blockpos_pp, ierr     )                                   &
    !$omp schedule(dynamic,1) collapse(4)
    L_MB: do m_bndelem = 1, bnd_elm_list%n_bnd_elements

      ! --- Select a test function (the weak form equation must hold for every test function)
      L_LS: do l_tor = i_tor_min, i_tor_max ! (loop over toroidal harmonics)
        L_LV: do l_vertex = 1, 2 ! (loop over nodes in element m_bndelem)
          L_LD: do l_dof = 1, 2 ! (loop over node dofs)

            bndelem_m = bnd_elm_list%bnd_element(m_bndelem)

            ! --- Determine the values of R,s and Z,s at the Gaussian points.
            call det_coord_bnd(bndelem_m, node_list, R_S=x_s, Z_S=y_s)

            l_node      = bndelem_m%vertex(l_vertex)
            l_dir       = bndelem_m%direction(l_vertex,l_dof)
            l_node_bnd  = bndelem_m%bnd_vertex(l_vertex)
            l_index     = node_list%node(l_node)%index(l_dir)
            l_size      = bndelem_m%size(l_vertex,l_dof)

            if ( (l_index < index_min) .or. (l_index > index_max) ) cycle ! This MPI proc responsible?

            ! --- Determine the row in the main matrix.
            l_row_psi = det_row_col(l_index, ivar_psi, l_tor, i_tor_min, i_tor_max)
            l_row_j   = det_row_col(l_index, ivar_j,   l_tor, i_tor_min, i_tor_max)

            ! --- Sum over boundary dofs at which response is calculated
            L_IV: do i_vertex = 1, 2 ! (loop over nodes in element m_bndelem)

              i_node      = bndelem_m%vertex(i_vertex)
              i_node_bnd  = bndelem_m%bnd_vertex(i_vertex)

              L_ID: do i_dof = 1, 2 ! (loop over node dofs)

                i_dir       = bndelem_m%direction(i_vertex,i_dof)
                i_index     = node_list%node(i_node)%index(i_dir)
                i_size      = bndelem_m%size(i_vertex,i_dof)

                L_IS: do i_starwall = 1, sr%n_tor ! (loop over STARWALL harmonics)

                  i_tor    = sr%i_tor(i_starwall)

                  if ( (i_tor < i_tor_min) .or. (i_tor > i_tor_max) ) cycle    

                  i_resp_old   = response_index(i_node_bnd,i_starwall,i_dof)

                  i_resp   = (bnd_node_list%bnd_node(i_node_bnd)%index_starwall(1) - 1)*sr%n_tor0 &
                           + bnd_node_list%bnd_node(i_node_bnd)%n_dof*(i_starwall-1)              &
                           + bnd_node_list%bnd_node(i_node_bnd)%index_starwall(i_dof)             &
                           - bnd_node_list%bnd_node(i_node_bnd)%index_starwall(1) + 1

!                  if (i_resp_old .ne. i_resp) write(*,'(A,8i5)') 'PANIC! : ',i_node, i_starwall, &
!                  i_dof,bnd_node_list%bnd_node(i_node_bnd)%index_starwall,i_resp_old, i_resp

                  i_resp_0 = response_index_eq(i_node_bnd,i_dof)

                  ! --- Loop over Gaussian points -- integration in s-direction
                  common_prefactor = 0.
                  L_MS: do ms = 1, n_gauss

                    ! --- Integration factor from the definition of dA:
                    !     int dA = sum_{m_bndelem} int ds int dphi  sqrt{(R,s)^2 + (Z,s)^2}
                    dA = sqrt(x_s(ms)**2 + y_s(ms)**2)

                    ! --- Loop over toroidal planes -- integration in phi-direction
                    L_MP: do m_plane = 1, n_plane

                      ! --- Evaluate test function at current position
                      testfunc_l = H1(l_vertex,l_dof,ms) * l_size * HZ(l_tor,m_plane)

                      ! --- Determine basis function
                      basfunc_i = H1(i_vertex,i_dof,ms) * i_size * HZ(i_tor,m_plane)

#ifdef __GFORTRAN__
                      common_prefactor = common_prefactor + wgauss_copy(ms) * dA * testfunc_l * basfunc_i
#else
                      common_prefactor = common_prefactor + wgauss(ms) * dA * testfunc_l * basfunc_i
#endif

                    end do L_MP

                  end do L_MS

                  ! --- Sum over boundary dofs contributing to the response
                  L_JB: do j_node_bnd = 1, bnd_node_list%n_bnd_nodes ! (loop over boundary nodes)

                    j_node      = bnd_node_list%bnd_node(j_node_bnd)%index_jorek

                    L_JD: do j_dof = 1, 2 ! (loop over node dofs)

                      j_dir       = bnd_node_list%bnd_node(j_node_bnd)%direction(j_dof)
                      j_index     = node_list%node(j_node)%index(j_dir)

                      L_JS: do j_starwall = 1, sr%n_tor ! (loop over STARWALL harmonics)

                        j_tor  = sr%i_tor(j_starwall)
                  
                        if ( (j_tor < i_tor_min) .or. (j_tor > i_tor_max) ) cycle    

                        j_resp   = (bnd_node_list%bnd_node(j_node_bnd)%index_starwall(1) - 1)*sr%n_tor0 &
                                 + bnd_node_list%bnd_node(j_node_bnd)%n_dof*(j_starwall-1)              &
                                 + bnd_node_list%bnd_node(j_node_bnd)%index_starwall(j_dof)             &
                                 - bnd_node_list%bnd_node(j_node_bnd)%index_starwall(1) + 1



                        ! --- Option to switch off mode coupling due to a 3D  wall
                        if ( vacuum_decouple_modes .and. (j_tor /= i_tor) ) cycle

                        ! --- Determine the column in the main matrix
                        j_col_psi = det_row_col(j_index, ivar_psi, j_tor, i_tor_min, i_tor_max)

                        ! --- Determine the position in the sparse matrix data structure
                        !     which corresponds to the matrix entry at  l_row_j, j_col_psi.
                        sparsepos_jp = det_sparse_pos(l_row_j,   j_col_psi, index_min, n_matrix_block_size, ijA_index, ijA_size, irn_jcn)
                        sparsepos_pp = det_sparse_pos(l_row_psi, j_col_psi, index_min, n_matrix_block_size, ijA_index, ijA_size, irn_jcn)

                        ! --- Vacuum response contribution to the lhs of the current equation
                        amat_contrib = - common_prefactor * response_m_e(i_resp, j_resp)
                        !$omp atomic
                        A_mat(sparsepos_jp) = A_mat(sparsepos_jp) + amat_contrib

                      end do L_JS
                    end do L_JD
                  end do L_JB

                  ! --- Contribution of vacuum response to the rhs of the current equation

                  rhs_contrib = sum( response_m_h(i_resp,:) * psibnd_vec(:))                     &
                              + sum( response_m_j(i_resp,:) * dpsibnd_vec(:))

                  if ( (l_tor == 1) .and. (sr%i_tor(1) == 1) .and. (.not. starwall_equil_coils)) &
                    rhs_contrib = rhs_contrib - sum( bext_tan(i_resp_0, :) * I_coils(:) )        &
                                - sum( response_m_h(i_resp,:) *   psibnd_coils(:) )

                  if ( resistive_wall )rhs_contrib=rhs_contrib+rhs_contrib_arr(i_resp)

                  rhs_contrib = rhs_contrib * common_prefactor
                  !$omp atomic
                  rhs_loc(l_row_j) = rhs_loc(l_row_j) + rhs_contrib

                end do L_IS
              end do L_ID
            end do L_IV

          end do L_LD
        end do L_LV
      end do L_LS

    end do L_MB
    !$omp end parallel do

    ! --- timing
    if ( vacuum_debug ) then
      t_elaps_end = MPI_WTIME()
      write(*,'(I4,A,F10.7,A)') my_id, '  Elapsed time vacuum_boundary_integral', t_elaps_end - t_elaps_start, ' s'
    end if
    
    if ( vacuum_debug ) write(*,*) my_id, 'After:', sum(abs(rhs_loc)),sum(abs(A_mat))

    if ( allocated(psibnd_vec ) ) deallocate( psibnd_vec  )
    if ( allocated(dpsibnd_vec) ) deallocate( dpsibnd_vec )

  end subroutine vacuum_boundary_integral
  
  

  !> Determine the values of \f$R\f$, \f$Z\f$, \f$\partial R/\partial s\f$, and
  !! \f$ \partial Z/\partial s \f$ on the Gaussian points of a given boundary element.
  subroutine det_coord_bnd(bndelem, node_list, R, Z, R_s, Z_s)
    
    use gauss,             only: n_gauss, xgauss, wgauss
    use data_structure,    only: type_node, type_bnd_element, type_node_list
    use basis_at_gaussian, only: H1, H1_s, HZ
    
    implicit none
    
    ! --- Routine parameters
    type(type_bnd_element), intent(in)  :: bndelem       !< Boundary element to be considered.
    type(type_node_list),   intent(in)  :: node_list     !< List of grid nodes
    real*8, optional,       intent(out) :: R(n_gauss)    !< Values of R on Gaussian points
    real*8, optional,       intent(out) :: Z(n_gauss)    !< Values of Z on Gaussian points
    real*8, optional,       intent(out) :: R_s(n_gauss)  !< Values of R,s on Gaussian points
    real*8, optional,       intent(out) :: Z_s(n_gauss)  !< Values of Z,s on Gaussian points
    
    ! --- Local variables
    integer         :: k_vertex, k_dof, k_node, k_dir
    real*8          :: k_size
    type(type_node) :: node_k
    
    if ( present(R  ) ) R   = 0.d0
    if ( present(Z  ) ) Z   = 0.d0
    if ( present(R_s) ) R_s = 0.d0
    if ( present(Z_s) ) Z_s = 0.d0
    
    do k_vertex = 1, 2
      do k_dof = 1, 2
        k_node      = bndelem%vertex(k_vertex)
        k_dir       = bndelem%direction(k_vertex,k_dof)
        k_size      = bndelem%size(k_vertex,k_dof)
        node_k      = node_list%node(k_node)
        if ( present(R  ) ) R  (:)  = R  (:)  + node_k%x(1,k_dir,1) * k_size * H1  (k_vertex,k_dof,:)
        if ( present(Z  ) ) Z  (:)  = Z  (:)  + node_k%x(1,k_dir,2) * k_size * H1  (k_vertex,k_dof,:)
        if ( present(R_s) ) R_s(:)  = R_s(:)  + node_k%x(1,k_dir,1) * k_size * H1_s(k_vertex,k_dof,:)
        if ( present(Z_s) ) Z_s(:)  = Z_s(:)  + node_k%x(1,k_dir,2) * k_size * H1_s(k_vertex,k_dof,:)
      end do
    end do
    
  end subroutine det_coord_bnd
  
  
  
  !> Determine vectors of the psi and deltapsi values at the boundary.
  subroutine det_psibnd_vec(bnd_node_list, node_list, psibnd_vec, dpsibnd_vec, psibnd_coils)

    use data_structure, only: type_node_list, type_bnd_node_list

    implicit none

    ! --- Routine parameters
    type(type_node_list),     intent(in)  :: node_list      !< List of grid nodes
    type(type_bnd_node_list), intent(in)  :: bnd_node_list  !< List of boundary grid nodes
    real*8, allocatable,      intent(out) :: psibnd_vec(:)  !< Vector of Psi boundary values
    real*8, allocatable,      intent(out) :: dpsibnd_vec(:) !< Vector of deltaPsi boundary values
    real*8, allocatable, optional, intent(out) :: psibnd_coils(:)!< Vector of Psi_coil boundary values

    ! --- Local variables
    integer :: jnode, jnode_glob, j_starwall, jtor, jbas, jdir, j_resp, j_resp_0, j_resp_old

    if ( allocated(psibnd_vec) ) deallocate(psibnd_vec)
    allocate( psibnd_vec(n_dof_starwall) )
    psibnd_vec(:) = 0.d0

    if ( allocated(dpsibnd_vec) ) deallocate(dpsibnd_vec)
    allocate( dpsibnd_vec(n_dof_starwall) )
    dpsibnd_vec(:) = 0.d0

    if ( present(psibnd_coils) ) then
      if ( allocated(psibnd_coils) ) deallocate(psibnd_coils)
      allocate( psibnd_coils(n_dof_starwall) )
      psibnd_coils(:) = 0.d0
    end if

    ! --- Determine vector of (delta)psi boundary values.
    do jnode = 1, bnd_node_list%n_bnd_nodes       ! loop over nodes

      jnode_glob = bnd_node_list%bnd_node(jnode)%index_jorek

      do j_starwall = 1, sr%n_tor     ! loop over STARWALL harmonics

        jtor = sr%i_tor(j_starwall)   ! (mode corresponding to STARWALL harmonic)

        do jbas = 1, 2                ! loop over basis functions

          jdir         = bnd_node_list%bnd_node(jnode)%direction(jbas)
          j_resp_old   = response_index(jnode,j_starwall,jbas)
          j_resp_0     = response_index_eq(jnode,jbas)

          j_resp = (bnd_node_list%bnd_node(jnode)%index_starwall(1) - 1)*sr%n_tor0 &
                 +  bnd_node_list%bnd_node(jnode)%n_dof*(j_starwall-1) &
                 + (bnd_node_list%bnd_node(jnode)%index_starwall(jbas)-bnd_node_list%bnd_node(jnode)%index_starwall(1)) + 1

!          if (j_resp_old .ne. j_resp) write(*,'(A4i5)') 'PANIC jresp: ',j_resp_old,j_resp, &
!           bnd_node_list%bnd_node(jnode)%index_starwall(1), bnd_node_list%bnd_node(jnode)%index_starwall(jbas)

          psibnd_vec ( j_resp ) = node_list%node(jnode_glob)%values(jtor, jdir, ivar_psi)
          dpsibnd_vec( j_resp ) = node_list%node(jnode_glob)%deltas(jtor, jdir, ivar_psi)

          if ( (present(psibnd_coils)) .and. (allocated(I_coils)) .and. (jtor==1) .and. (.not. starwall_equil_coils) ) then
            j_resp_0 = 2*(jnode-1) + jbas
            psibnd_coils( j_resp ) = sum( bext_psi(j_resp_0,:) * I_coils(:) )
          end if

        end do
      end do
    end do
    
  end subroutine det_psibnd_vec
  
  
  
  !> Initialize the currents in the resistive wall and the external coil currents.
  subroutine init_wall_currents(my_id, resistive_wall)
    
    use nodes_elements, only: bnd_node_list, node_list
    use phys_module, only: index_now, index_start, nstep
    use mpi_mod
    
    implicit none
    
    ! --- Routine parameters
    integer, intent(in) :: my_id
    logical, intent(in) :: resistive_wall
    
    ! --- Local variables
    real*8, allocatable :: psibnd_vec(:)    ! Vector of Psi values at the boundary
    real*8, allocatable :: dpsibnd_vec(:)   ! Vector of deltaPsi values at the boundary
    real*8, allocatable :: wall_and_coil_curr(:)
    real*8, allocatable :: tmp_coil_curr(:)
    integer             :: k, k2, global_index, ierr
    
    call det_psibnd_vec(bnd_node_list, node_list, psibnd_vec, dpsibnd_vec)
    
    if ( resistive_wall ) then
      
      if ( (allocated(wall_curr)) .and. (.not. wall_curr_initialized) ) &
        deallocate(wall_curr)
      if ( (allocated(dwall_curr)) .and. (.not. wall_curr_initialized) ) &
        deallocate(dwall_curr)
      
      if ( .not. allocated(wall_curr) ) then
        allocate( wall_curr(n_wall_curr) )
        do k = 1, n_wall_curr
          wall_curr(k) = 0.d0! - sum( sr%a_ye(k,:) * psibnd_vec(:) )
        end do
      end if
      
      if ( .not. allocated(dwall_curr) ) then
        allocate( dwall_curr(n_wall_curr) )
        do k = 1, n_wall_curr
          dwall_curr(:) = 0.d0!- sum( sr%a_ye(k,:) * dpsibnd_vec(:) )
        end do
      end if         
      
    end if
    
    ! --- Also initialize the old_dpsibnd_vec
    if ( allocated(old_dpsibnd_vec) ) deallocate(old_dpsibnd_vec)
    allocate( old_dpsibnd_vec(n_dof_starwall) )
    old_dpsibnd_vec(:) = 0.d0
    
    call write_wall_vtk(0, resistive_wall, my_id)
    deallocate( psibnd_vec, dpsibnd_vec )

    ! --- Initialize net toroidal wall current (for plot_live_data)
    if ( resistive_wall .and. (index_now>0) .and. (index_start+nstep>0)) then
      if ( .not. allocated(net_tor_wall_curr) ) then
        allocate( net_tor_wall_curr(index_start+nstep) )
        net_tor_wall_curr(:) = 0.d0
      end if
      k2 = sr%ncoil + 1 
      if ( (k2 >= sr%s_ww%ind_start) .and. (k2 <= sr%s_ww%ind_end) ) then
        global_index = my_id*sr%s_ww%step
        net_tor_wall_curr(index_now) = sum(sr%s_ww%loc_mat(k2 - global_index,:) * wall_curr(:))
      endif
      call MPI_ALLReduce(MPI_IN_PLACE, net_tor_wall_curr(index_now),1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    end if
    
    ! --- Initialize coil currents (for plot_live_data)
    if ( resistive_wall .and. (sr%ncoil > 0) .and. (index_now>0) .and. (index_start+nstep>0)) then
      ! Calculate coil currents for all coil types
      if  (.not. allocated(tmp_coil_curr) ) allocate(tmp_coil_curr(sr%ncoil))
      tmp_coil_curr = 0.d0
      do k = 1, sr%ncoil
        if ( (k >= sr%s_ww%ind_start) .and. (k <= sr%s_ww%ind_end) ) then
          global_index = my_id*sr%s_ww%step
          tmp_coil_curr(k) = sum(sr%s_ww%loc_mat(k - global_index,:) * wall_curr(:))
        endif
      end do
      call MPI_ALLReduce(MPI_IN_PLACE, tmp_coil_curr,sr%ncoil,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

      ! Distribute the coil currents to the different coil types
      if ( sr%n_diag_coils > 0 ) then
        if ( (.not. allocated(diag_coil_curr)) .and. (index_start+nstep >0) ) then
          allocate( diag_coil_curr(index_start+nstep,sr%n_diag_coils) )
          diag_coil_curr(:,:) = 0.d0
        end if
        diag_coil_curr(index_now,:) =  tmp_coil_curr(sr%ind_start_diag_coils:sr%ind_start_diag_coils + sr%n_diag_coils -1)
      endif
      
      if ( sr%n_rmp_coils > 0 ) then
        if ( (.not. allocated(rmp_coil_curr)) .and. (index_start+nstep >0) ) then
          allocate( rmp_coil_curr(index_start+nstep,sr%n_rmp_coils) )
          rmp_coil_curr(:,:) = 0.d0
        end if
        rmp_coil_curr(index_now,:) =  tmp_coil_curr(sr%ind_start_rmp_coils : sr%ind_start_rmp_coils + sr%n_rmp_coils -1) 
      end if
      
      if ( sr%n_pol_coils > 0 ) then
        if ( (.not. allocated(pf_coil_curr)) .and. (index_start+nstep >0) ) then
          allocate( pf_coil_curr(index_start+nstep,sr%n_pol_coils) )
          pf_coil_curr(:,:) = 0.d0
        end if
        pf_coil_curr(index_now,:) =  tmp_coil_curr(sr%ind_start_pol_coils:sr%ind_start_pol_coils + sr%n_pol_coils -1) 
      end if

      ! Voltage coils not yet implemented
      ! if ( sr%n_voltage_coils > 0 ) then
      !  if ( (.not. allocated(voltage_coil_curr)) .and. (index_start+nstep >0) ) then
      !    allocate( voltage_coil_curr(index_start+nstep,sr%n_voltage_coils) )
      !    voltage_coil_curr(:,:) = 0.d0
      !  end if
      !  voltage_coil_curr(index_now,:) =  tmp_coil_curr(sr%ind_start_voltage_coils:sr%ind_start_voltage_coils + sr%n_voltage_coils -1) 
      !  endif
      !end if
    endif
          
 
    if ( vacuum_debug .and. (my_id == 0) ) write(*,*) 'Wall currents initialized.'
    wall_curr_initialized = .true.
    
  end subroutine init_wall_currents
  
  
  
  !> Imposing PF coil currents as a current source term
  subroutine coil_current_source(my_id)
    
    use constants
    use phys_module, only: t_now, index_now, tstep, nstep, index_start
    use equil_info, only: ES
    use profiles
    use mpi_mod    

    implicit none
   
    ! --- Routine parameters
    integer, intent(in)       :: my_id  
    integer                   :: i, ierr
    real*8, allocatable       :: potentials_real_0(:)
    real*8, save, allocatable :: delta_Icoils_0(:)
    real*8, save              :: t_last, Z_p = 1.d10
    logical, save             :: initialized=.false.
    real*8                    :: z_ref_inter, Z_n
    if( my_id == 0 ) write(*,*) ' Imposing PF coil currents with a current source term '

    ! --- Calculate the difference between the input file coil currents and the ones coming from the restart file
    ! This is used below to avoid violent jumps in I_coils that can happen when restarting from a freeboudary_equil with feedback
    if(.not. initialized) then
      if (allocated(delta_Icoils_0)) deallocate(delta_Icoils_0)
      allocate(delta_Icoils_0(n_coils))
      delta_Icoils_0 = 0.d0
      do i=1, n_coils
        if ((vert_FB_amp(i) /= 0.d0) .or. (vert_FB_amp_ts(i)/=0.d0) ) then
          delta_Icoils_0(i) = I_coils(i) - interpolProf(coil_curr_time_trace(i)%time, coil_curr_time_trace(i)%curr, coil_curr_time_trace(i)%len, t_now)
        end if
      end do
      initialized = .true.
    end if
    

    ! --- Calculate the specified coil currents at present time 
    do i=1, n_coils
      I_coils(i) = interpolProf(coil_curr_time_trace(i)%time, coil_curr_time_trace(i)%curr, coil_curr_time_trace(i)%len, t_now) 
    end do

    ! -- During timestepping apply vertical feedback by the PF coils which were activated
    ! Add the feedback current to the difference from the restart file
    if (.not. allocated(vert_FB_response) .and. index_start+nstep>0) then
      allocate(vert_FB_response(index_start+nstep,4))
      vert_FB_response = 0.d0
    endif
    z_ref_inter = interpolProf(Z_axis_ref_ts%time, Z_axis_ref_ts%position ,  Z_axis_ref_ts%len, t_now)
    vert_FB_response(index_now,4) = z_ref_inter
    if ( t_now>start_VFB_ts .and. sum( abs(vert_FB_amp_ts(1:n_pf_coils)) )>1.e-6 .and. sr%i_tor(1) ==  1 ) then
      Z_n = ES%Z_axis
      if ( Z_p > 1.d9 ) Z_p = Z_n
      dZ_axis_integral = dZ_axis_integral + ( Z_n - z_ref_inter )*tstep
      if ( (t_now-t_last)>vert_FB_tact ) then 
        if (my_id==0) write(*,*) 'Vertical feedback active'
        t_last = t_now
        vert_FB_response(index_now,1:3) =  &
            (/  1.d3  * vert_FB_gain(1)* (Z_n - z_ref_inter), &
                1.d5  * vert_FB_gain(2)* (Z_n - Z_p )/ tstep, &
                1.d-3 * vert_FB_gain(3)* dZ_axis_integral /)
        delta_Icoils_0(1:n_pf_coils) = delta_Icoils_0(1:n_pf_coils) +  sum(vert_FB_response(index_now,1:3))*vert_FB_amp_ts(1:n_pf_coils)
      endif
      Z_p = Z_n
    endif

    

    ! -- Apply coil current limits
    do i=1, n_coils
      if ( abs( I_coils(i) + delta_Icoils_0(i)) .gt. I_coils_max(i) ) &
          delta_Icoils_0(i) = ( I_coils(i) + delta_Icoils_0(i) ) / abs( I_coils(i) + delta_Icoils_0(i))&
                                  *I_coils_max(i) - ( I_coils(i) + delta_Icoils_0(i) )
      I_coils(i) = I_coils(i)  + delta_Icoils_0(i)
    end do


    if (.not. allocated (Y_coils0)) allocate(Y_coils0(n_wall_curr))
    if (.not. allocated (potentials_real_0)) allocate(potentials_real_0(n_wall_curr))
    ! --- Transform real currents into starwall currents
    Y_coils0          = 0.d0
    potentials_real_0 = 0.d0
    potentials_real_0(1:n_coils) = I_coils(1:n_coils) * mu_zero
    do i = 1, n_wall_curr
      if ( (i>=sr%s_ww_inv%ind_start) .and. (i<=sr%s_ww_inv%ind_end) ) then
        Y_coils0(i) = sum(sr%s_ww_inv%loc_mat(i-my_id*sr%s_ww_inv%step,:)*potentials_real_0(:))
      end if
    end do
    call MPI_ALLReduce(MPI_IN_PLACE, Y_coils0, size(Y_coils0),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

  end subroutine coil_current_source
  
  
  
  !> Perform the time-evolution of the wall currents (resistive wall).
  subroutine evolve_wall_currents(my_id, psibnd_vec, dpsibnd_vec)
    
    use phys_module, only: index_now, index_start, nstep, resistive_wall
    use mpi_mod
    
    implicit none
    
    ! --- Routine parameters
    integer, intent(in) :: my_id
    real*8,  intent(in) :: psibnd_vec (n_dof_starwall) !< Vector of Psi boundary values
    real*8,  intent(in) :: dpsibnd_vec(n_dof_starwall) !< Vector of deltaPsi boundary values
    
    ! --- Local variables
    integer              :: k, k2, global_index
    integer              :: ierr
    real *8, allocatable :: dwall_curr_2(:)
    real *8, allocatable :: tmp_coil_curr(:)

    if  (.not. allocated(dwall_curr_2)) allocate (dwall_curr_2(n_wall_curr))
    if ( vacuum_debug ) write(*,*) 'wall_curr(before)', sum(abs(wall_curr)), sum(wall_curr)    
    
    if ( (.not. allocated(old_dpsibnd_vec)) .or. (size(old_dpsibnd_vec,1)/= n_dof_starwall) ) then
      if ( allocated(old_dpsibnd_vec) ) deallocate(old_dpsibnd_vec)
      allocate( old_dpsibnd_vec(n_dof_starwall) )
    end if    

     do k = 1, n_wall_curr

       ! --- contribution from not distributed matrices
       dwall_curr_2(k) = response_d_b(k) * (wall_curr(k) - Y_coils0(k)) &
         + response_d_c(k) * dwall_curr(k)
       
       ! --- contribution from distributed matrices
       dwall_curr(k) = sum(response_m_a%loc_mat(k,:) * dpsibnd_vec(response_m_a%ind_start:response_m_a%ind_end)) &
         + sum(response_m_d%loc_mat(k,:) * old_dpsibnd_vec(response_m_d%ind_start:response_m_d%ind_end) )

     end do

    call MPI_ALLReduce(MPI_IN_PLACE, dwall_curr,size(dwall_curr),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

    dwall_curr(:)= dwall_curr(:)+dwall_curr_2(:)

    if (index_now <= 1 ) dwall_curr = 0.d0
    wall_curr(:) = wall_curr(:) + dwall_curr(:)

    call write_wall_vtk(index_now, resistive_wall, my_id)

    if ( vacuum_debug .and. resistive_wall ) then
      call log_wall_curr(my_id)
      !call log_coil_curr()
    end if

    ! --- Extract net toroidal wall current such that it can be written to the macroscopic_vars.dat file (e.g., for ./util/plot_live_data.sh)
    if ( (.not. allocated(net_tor_wall_curr)) .and. (index_start+nstep >0) ) then
      allocate( net_tor_wall_curr(index_start+nstep) )
      net_tor_wall_curr(:) = 0.d0
    end if
    if (index_now>0) then
      k2 = sr%ncoil + 1
      if ( (k2 >= sr%s_ww%ind_start) .and. (k2 <= sr%s_ww%ind_end) ) then
        global_index = my_id*sr%s_ww%step
        net_tor_wall_curr(index_now) = sum(sr%s_ww%loc_mat(k2 - global_index,:) * wall_curr(:))
      endif
    endif
    call MPI_ALLReduce(MPI_IN_PLACE, net_tor_wall_curr(index_now),1 ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

    ! --- Extract coil currents such that they can be written to the macroscopic_vars.dat file (e.g., for ./util/plot_live_data.sh)
    if (sr%ncoil > 0) then
      ! Calculate coil currents for all coil types
      if  (.not. allocated(tmp_coil_curr) ) allocate(tmp_coil_curr(sr%ncoil))
      tmp_coil_curr = 0.d0
      do k = 1, sr%ncoil
        if ( (k >= sr%s_ww%ind_start) .and. (k <= sr%s_ww%ind_end) ) then
          global_index = my_id*sr%s_ww%step
          tmp_coil_curr(k) = sum(sr%s_ww%loc_mat(k - global_index,:) * wall_curr(:))
        endif
      end do
      call MPI_ALLReduce(MPI_IN_PLACE, tmp_coil_curr,sr%ncoil,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

      ! Distribute the coil currents to the different coil types
      if ( sr%n_diag_coils > 0 ) then
        if ( (.not. allocated(diag_coil_curr)) .and. (index_start+nstep >0) ) then
          allocate( diag_coil_curr(index_start+nstep,sr%n_diag_coils) )
          diag_coil_curr(:,:) = 0.d0
        end if
        if (index_now>0) then
          diag_coil_curr(index_now,:) =  tmp_coil_curr(sr%ind_start_diag_coils:sr%ind_start_diag_coils + sr%n_diag_coils -1) 
        endif
      end if
      
      if ( sr%n_rmp_coils > 0 ) then
        if ( (.not. allocated(rmp_coil_curr)) .and. (index_start+nstep >0) ) then
          allocate( rmp_coil_curr(index_start+nstep,sr%n_rmp_coils) )
          rmp_coil_curr(:,:) = 0.d0
        end if
        if (index_now>0) then
          rmp_coil_curr(index_now,:) =  tmp_coil_curr(sr%ind_start_rmp_coils : sr%ind_start_rmp_coils + sr%n_rmp_coils -1) 
        endif
      end if
      
      if ( sr%n_pol_coils > 0 ) then
        if ( (.not. allocated(pf_coil_curr)) .and. (index_start+nstep >0) ) then
          allocate( pf_coil_curr(index_start+nstep,sr%n_pol_coils) )
          pf_coil_curr(:,:) = 0.d0
        end if
        if (index_now>0) then
          pf_coil_curr(index_now,:) =  tmp_coil_curr(sr%ind_start_pol_coils:sr%ind_start_pol_coils + sr%n_pol_coils -1) 
        endif
      end if

      ! Voltage coils not yet implemented
      !if ( sr%n_voltage_coils > 0 ) then
      !  if ( (.not. allocated(voltage_coil_curr)) .and. (index_start+nstep >0) ) then
      !    allocate( voltage_coil_curr(index_start+nstep,sr%n_voltage_coils) )
      !    voltage_coil_curr(:,:) = 0.d0
      !  end if
      !  if (index_now>0) then
      !    voltage_coil_curr(index_now,:) =  tmp_coil_curr(sr%ind_start_voltage_coils:sr%ind_start_voltage_coils + sr%n_voltage_coils -1) 
      !  endif
      !end if
    endif


    if ( vacuum_debug .and. (my_id ==0) ) write(*,*) 'wall_curr(after)', sum(abs(wall_curr)), sum(wall_curr)
    
  end subroutine evolve_wall_currents
  
  
  
  !> Reconstruct the potential values at the wall triangle nodes.
  !!
  !! To reconstruct the physical wall potentials from the ones we work 
  !! with we need to muliply them by the similarity transform matrix.
  !! The physical wall potentials include different types of potentials,
  !! and the way they are ordered in the array is
  !!
  !!   (I_coil_1, I_coil_2, ..., I_coil_ncoil, Iw_net_tor, Potw_1, Potw_2, ..., Potw_npotw-1)
  !!  
  !! where I_coil are the coil currents, Iw_net is the net wall current
  !! and Potw are the single valued wall potentials.  
  subroutine reconstruct_triangle_potentials(tripot_w, wall_curr, my_id, Iw_net_tor)
    
    use mpi_mod

    implicit none
    
    ! --- Routine parameters
    real*8, allocatable, intent(inout) :: tripot_w(:)
    real*8, allocatable, intent(in)    :: wall_curr(:)
    integer,             intent(in)    :: my_id
    real*8,  optional,   intent(inout) :: Iw_net_tor 

    ! --- Local variables
    integer              :: i, j, ierr, global_index, ntasks
    integer              :: count=1
    real*8, allocatable  :: pot_tmp(:)

    if ( allocated(tripot_w) ) deallocate(tripot_w); allocate( tripot_w(sr%npot_w) )    
    if ( allocated(pot_tmp)  ) deallocate(pot_tmp);  allocate(  pot_tmp(sr%npot_w) )    
    tripot_w = 0.d0
    pot_tmp  = 0.d0

    ! --- multiply by the similarity transform matrix to get the physical wall potentials
    if ( allocated(wall_curr) ) then
      global_index = my_id*sr%s_ww%step

      do i = 1, sr%npot_w
        j = i + sr%ncoil
        if ( (j >= sr%s_ww%ind_start) .and. (j <= sr%s_ww%ind_end) ) then
          pot_tmp(i) = sum(sr%s_ww%loc_mat(j - global_index,:) * wall_curr(:))
        end if  
      end do
     
    end if

    call MPI_AllREDUCE(MPI_IN_PLACE,pot_tmp,size(pot_tmp),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

    if (present(Iw_net_tor)) then
      Iw_net_tor = pot_tmp(1)
    endif

    ! --- Correct indexing (1st potential is the net curret potential)
    ! --- Shift indexing so 1st wall node corresponds to 1st single valued potential (and so on)
    do i = 1, sr%npot_w-1
      tripot_w(i) = pot_tmp(i+1)
    enddo

    tripot_w(sr%npot_w) = 0.d0 ! STARWALL's BC for the single valued potentials (last potential is 0)
    !!! Net current potential contribution to be added soon!

    deallocate(pot_tmp) 

  end subroutine reconstruct_triangle_potentials
  
  
  
  
  !> Write wall current potentials to logfile.
  subroutine log_wall_curr(my_id)
    
    implicit none
    
    integer,   intent(in) :: my_id

    ! --- Local variables
    real*8, allocatable :: tripot_w(:)
    
    call reconstruct_triangle_potentials(tripot_w, wall_curr, my_id)
    if (my_id==0) write(*,'(" Wall current potentials (min, max): ",ES16.8,"...",ES16.8)') minval(tripot_w), maxval(tripot_w)
    deallocate( tripot_w )
    
  end subroutine log_wall_curr
  
  
  
  !> Update vacuum response
  !! This is necessary
  !! - right after the start or restart of the code
  !! - when wall resistivity, tstep, or some other parameters have changed.
  subroutine update_response(my_id, tstep, freeboundary_equil, resistive_wall)

    use phys_module, only: time_evol_theta, time_evol_zeta
    use mpi_mod

    implicit none

    ! --- Routine parameters
    integer,                     intent(in) :: my_id              !< MPI task ID
    real*8,                      intent(in) :: tstep              !< delta t,timestep
    logical,                     intent(in) :: freeboundary_equil !< Use free boundary equilibrium?
    logical,                     intent(in) :: resistive_wall     !< Resistive or ideal wall?

    ! --- Local variables
    integer :: i, j, k, j2, k2
    real*8, allocatable :: tmp_d_s(:)
    real*8  :: theta, zeta
    logical :: update_required
    integer :: ierr,loc_size,ntasks
    real*8  :: test_sum,test_sum2

! --- Local variables to store the previous values of some parameters.
    real*8,  save :: old_thick   = 0.0
    real*8,  save :: old_res     = 0.0
    real*8,  save :: old_tstep   = 0.0
    real*8,  save :: old_theta   = 0.0
    real*8,  save :: old_zeta    = 0.0
    logical, save :: old_reswall = .false.

    if ( sr%n_tor == 0 ) then
      write(*,*) 'Remark: Routine update_response is not doing anything since sr%n_tor==0.'
      return
    end if

    theta = time_evol_theta
    zeta  = time_evol_zeta

    call MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, ierr)

! --- Update response matrices only, if parameter values changed or matrices not allocated
    update_required = ( old_res     /= wall_resistivity       ) &
                 .or. ( old_tstep   /= tstep                  ) &
                 .or. ( old_theta   /= theta                  ) &
                 .or. ( old_zeta    /= zeta                   ) &
                 .or. ( old_reswall .neqv. resistive_wall     ) &
                 .or. ( .not. allocated(response_m_a%loc_mat) ) &
                 .or. ( .not. allocated(response_d_b)         ) &
                 .or. ( .not. allocated(response_d_c)         ) &
                 .or. ( .not. allocated(response_m_d%loc_mat) ) &
                 .or. ( .not. allocated(response_m_e)         ) &
                 .or. ( .not. allocated(response_m_f%loc_mat) ) &
                 .or. ( .not. allocated(response_m_g%loc_mat) ) &
                 .or. ( .not. allocated(response_m_h)         ) &
                 .or. ( .not. allocated(response_m_j)         ) &
                 .or. ( .not. allocated(response_m_k)         ) &
                 .or. ( .not. allocated(response_m_l)         ) &
                 .or. ( .not. allocated(response_m_v%loc_mat) ) &
                 .or. ( .not. allocated(response_m_eq)        )

    if ( update_required ) then
      ! --- Remember parameter values.
      old_res     = wall_resistivity
      old_tstep   = tstep
      old_theta   = theta
      old_zeta    = zeta
      old_reswall = resistive_wall

      ! --- Allocate matrices if required
      if ( .not. allocated(response_m_eq)) allocate( response_m_eq(sr%nd_bez/sr%n_tor, sr%nd_bez/sr%n_tor) )
      if ( .not. allocated(response_d_b) ) allocate( response_d_b(n_wall_curr) )
      if ( .not. allocated(response_d_c) ) allocate( response_d_c(n_wall_curr) )
      if ( .not. allocated(response_m_h) ) allocate( response_m_h(n_dof_starwall, n_dof_starwall) )
      if ( .not. allocated(response_m_k) ) allocate( response_m_k(n_wall_curr, sr%ncoil) )   

      if( .not. allocated (response_m_a%loc_mat) ) &
            call  alloc_distr(my_id, response_m_a, (/n_wall_curr, n_dof_starwall/), .false.)
      if( .not. allocated (response_m_d%loc_mat) ) &
            call  alloc_distr(my_id, response_m_d, (/n_wall_curr, n_dof_starwall/), .false.)
      if( .not. allocated (response_m_f%loc_mat) ) &
            call  alloc_distr(my_id, response_m_f, (/n_dof_starwall, n_wall_curr/), .true.)
      if( .not. allocated (response_m_g%loc_mat) ) & 
            call  alloc_distr(my_id, response_m_g, (/n_dof_starwall, n_wall_curr/), .true.)
      if( .not. allocated (response_m_v%loc_mat) ) &
            call  alloc_distr(my_id, response_m_v, (/n_dof_starwall, n_wall_curr/), .true.)

      ! --- Derived response matrix for equilibrium (extract n=0 part from STARWALL EE matrix)
      response_m_eq = 0.d0

      do j = 1, sr%nd_bez/sr%n_tor0, 2
        j2 = (j-1)*sr%n_tor0 + 1
        do k = 1, sr%nd_bez/sr%n_tor0, 2
          k2 = (k-1)*sr%n_tor0 + 1

           if ( sr%a_ee%row_wise ) then

             if ( (j2>=sr%a_ee%ind_start) .and. (j2<=sr%a_ee%ind_end) ) then
               response_m_eq(j,k:k+1) = sr%a_ee%loc_mat(j2-sr%a_ee%step*my_id,k2:k2+1)
             end if

             if( ((j2+1)>=sr%a_ee%ind_start) .and. ((j2+1)<=sr%a_ee%ind_end) ) then
               response_m_eq(j+1,k:k+1) = sr%a_ee%loc_mat((j2+1)-sr%a_ee%step*my_id,k2:k2+1)
             end if

           end if
        
        end do
      end do

      call MPI_AllREDUCE(MPI_IN_PLACE,response_m_eq,size(response_m_eq),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)

      ! --- Derived response matrices for time-evolution
      if ( resistive_wall) then

        allocate( tmp_d_s(n_wall_curr) )

        tmp_d_s(:) = 1.d0 + zeta + tstep * theta * wall_resistivity * sr%d_yy(:)

        do j = 1, size(response_m_a%loc_mat,2) 
          response_m_a%loc_mat(:,j) = -(1.d0+zeta) * sr%a_ye%loc_mat(:,j) / tmp_d_s(:)
        end do

        response_d_b(:) = - tstep * wall_resistivity * sr%d_yy(:) / tmp_d_s(:)
        response_d_c(:) = zeta / tmp_d_s(:)

        do j = 1, size(response_m_d%loc_mat,2)
          response_m_d%loc_mat(:,j) = zeta * sr%a_ye%loc_mat(:,j) / tmp_d_s(:)
        end do
        
        ! --- m_e = a_ee - matmul(a_ey, m_a)
        call matrix_multiplication(my_id,sr%a_ey,mat2=response_m_a, res_mat_not_distr=response_m_e)
     
        do j = 1, size(sr%a_ee%loc_mat,1)
          response_m_e(j+my_id*sr%a_ee%step,:) = response_m_e(j+my_id*sr%a_ee%step,:) + sr%a_ee%loc_mat(j,:)
        end do
        call MPI_AllREDUCE(MPI_IN_PLACE,response_m_e,size(response_m_e),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)

        do k = 1, n_wall_curr
          response_m_f%loc_mat(:,k) = sr%a_ey%loc_mat(:,k) * ( 1.d0 + response_d_b(k) )
        end do

        do k = 1, n_wall_curr
          response_m_g%loc_mat(:,k) = sr%a_ey%loc_mat(:,k) * response_d_c(k)
        end do

        response_m_h =0.0
        response_m_h(sr%a_ee%ind_start:sr%a_ee%ind_end,:) = sr%a_ee%loc_mat(:,:)
        call MPI_AllREDUCE(MPI_IN_PLACE,response_m_h,size(response_m_h),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)

        ! --- m_j =  matmul(a_ey, m_d)
        call matrix_multiplication(my_id,sr%a_ey,mat2=response_m_d,res_mat_not_distr=response_m_j)
        call MPI_AllREDUCE(MPI_IN_PLACE,response_m_j,size(response_m_j),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)

        response_m_k=0.0        
        do k = 1, n_wall_curr
          do j = 1, sr%ncoil
            if(sr%s_ww%ind_start<=j .AND. sr%s_ww%ind_end>=j) then   
              response_m_k(k,j) = -tstep * sr%d_yy(k) * sr%s_ww%loc_mat(j,k)
            end if
          end do
        end do
        call MPI_AllREDUCE(MPI_IN_PLACE,response_m_k,size(response_m_k),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)

        call matrix_multiplication(my_id,sr%a_ey,mat2_not_distr=response_m_k,res_mat_not_distr=response_m_l)
        call MPI_AllREDUCE(MPI_IN_PLACE,response_m_l,size(response_m_l),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)

        do k = 1, n_wall_curr
          response_m_v%loc_mat(:,k) = sr%a_ey%loc_mat(:,k) * ( response_d_b(k) )
        end do
        
        deallocate( tmp_d_s )

      else ! (Ideal wall)

        response_m_a%loc_mat(:,:) = 0.d0
        response_d_b(:)   = 0.d0
        response_d_c(:)   = 0.d0
        response_m_d%loc_mat(:,:) = 0.d0

        response_m_e=0.0
        response_m_e(sr%a_id%ind_start:sr%a_id%ind_end,:) =sr%a_id%loc_mat(:,:)
        call MPI_AllREDUCE(MPI_IN_PLACE,response_m_e,size(response_m_e),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)

        response_m_f%loc_mat(:,:) = 0.d0
        response_m_g%loc_mat(:,:) = 0.d0

        response_m_h =0.0
        response_m_h(sr%a_id%ind_start:sr%a_id%ind_end,:) = sr%a_id%loc_mat(:,:)
        call MPI_AllREDUCE(MPI_IN_PLACE,response_m_h,size(response_m_h),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)

        response_m_j(:,:) = 0.d0
        response_m_k(:,:) = 0.d0 !####
        response_m_l(:,:) = 0.d0 !####
        response_m_v%loc_mat(:,:) = 0.d0 !####

      end if

      if ( vacuum_debug ) then
        if(my_id == 0)   write(*,*) 'DEBUG: Checksums'

        test_sum=0.0
        test_sum2=0.0
        call MPI_AllREDUCE(sum(abs(response_m_a%loc_mat)),test_sum2,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_AllREDUCE(sum(response_m_a%loc_mat),test_sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
        if(my_id == 0)   write(*,*) 'm_a:', test_sum2, test_sum
        if(my_id == 0)   write(*,*) 'd_b:', sum(abs(response_d_b)), sum(response_d_b)
        if(my_id == 0)   write(*,*) '1+d_b:',sum(abs(1.d0+response_d_b)), sum(1.d0+response_d_b)
        if(my_id == 0)   write(*,*) 'd_c:', sum(abs(response_d_c)), sum(response_d_c)

        test_sum=0.0
        test_sum2=0.0
        call MPI_AllREDUCE(sum(abs(response_m_d%loc_mat)),test_sum2,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_AllREDUCE(sum(response_m_d%loc_mat),test_sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
        if(my_id == 0)   write(*,*) 'm_d:', test_sum2, test_sum
        if(my_id == 0)   write(*,*) 'm_e:', sum(abs(response_m_e)), sum(response_m_e)

        test_sum=0.0
        test_sum2=0.0
        call MPI_AllREDUCE(sum(abs(response_m_f%loc_mat)),test_sum2,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_AllREDUCE(sum(response_m_f%loc_mat),test_sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
        if(my_id == 0)   write(*,*) 'm_f:',test_sum2, test_sum

        test_sum=0.0
        test_sum2=0.0
        call MPI_AllREDUCE(sum(abs(response_m_g%loc_mat)),test_sum2,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_AllREDUCE(sum(response_m_g%loc_mat),test_sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
        if(my_id == 0)   write(*,*) 'm_g:',test_sum2, test_sum
        if(my_id == 0)   write(*,*) 'm_h:', sum(abs(response_m_h)), sum(response_m_h)
        if(my_id == 0)   write(*,*) 'm_j:', sum(abs(response_m_j)), sum(response_m_j)
        if(my_id == 0)   write(*,*) 'm_k:', sum(abs(response_m_k)), sum(response_m_k)
        if(my_id == 0)   write(*,*) 'm_l:', sum(abs(response_m_l)), sum(response_m_l)
    
        test_sum=0.0
        test_sum2=0.0
        call MPI_AllREDUCE(sum(abs(response_m_v%loc_mat)),test_sum2,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_AllREDUCE(sum(response_m_v%loc_mat),test_sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
        if(my_id == 0)   write(*,*) 'm_v:',test_sum2, test_sum
        if(my_id == 0)   write(*,*) 'm_eq:', sum(abs(response_m_eq)), sum(response_m_eq)
        if(my_id == 0)   write(*,*) 'END: Checksums'
      end if

    end if  

   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  end subroutine update_response
  
  
  !> Read a STARWALL response matrix from a file.
  subroutine read_response_matrix( matrix, dim_expected, filename, err )
    
    ! --- Routine parameters
    real*8, allocatable, intent(inout) :: matrix(:,:)     !< Matrix to be read
    integer,             intent(in)    :: dim_expected(2) !< Matrix dimension expected
    character(len=*),    intent(in)    :: filename        !< Filename to read from
    integer,             intent(out)   :: err             !< Error code
    
    ! --- Local variables
    integer :: dim(2), i, j, i2, j2
    
    err = 0
    
    open(42, FILE=trim(filename), status='old', action='read', iostat=err)
    if ( err /= 0 ) return

    read(42,*) dim

    if ( (dim(1) /= dim_expected(1)) .or. (dim(2) /= dim_expected(2)) ) then
      write(*,*) 'FATAL ERROR: Matrix dimension not as expected. Different resolutions?'
      write(*,'(1x,A,2I7)') 'dim_expected=',dim_expected
      stop
    end if
    
    if ( allocated(matrix) ) deallocate(matrix)
    allocate( matrix(dim(1),dim(2)) )
    matrix = 0.d0
    
    do i = 1, dim(1)
      do j = 1, dim(2)
        
        read (42,*) i2, j2, matrix(i,j)
        
        if ( ( i2 /= i ) .or. ( j2 /= j ) ) then
          write(*,*) 'FATAL ERROR: Matrix indices not as expected. Different resolutions?'
          stop
        end if
      
      end do
    end do
      
    close(42)

  end subroutine read_response_matrix
  
  
  
  !> Read a diagonal STARWALL response matrix from a file.
  subroutine read_response_diagonal( diagonal, dim_expected, filename, err )
    
    ! --- Routine parameters
    real*8, allocatable, intent(inout) :: diagonal(:)     !< Matrix to be read
    integer,             intent(in)    :: dim_expected    !< Matrix dimension expected
    character(len=*),    intent(in)    :: filename        !< Filename to read from
    integer,             intent(out)   :: err             !< Error code
    
    ! --- Local variables
    integer :: dim, i, i2
    
    open(42, FILE=trim(filename), status='old', action='read', iostat=err)
    if ( err /= 0 ) return

    read(42,*) dim

    if ( dim /= dim_expected ) then
      write(*,*) 'FATAL ERROR: Matrix dimension not as expected. Different resolutions?'
      stop
    end if
    
    if ( allocated(diagonal) ) deallocate(diagonal)
    allocate( diagonal(dim) )
    diagonal = 0.d0
    
    do i = 1, dim
        
      read (42,*) i2, diagonal(i)
      
      if ( i2 /= i ) then
        write(*,*) 'FATAL ERROR: Matrix indices not as expected. Different resolutions?'
        stop
      end if
        
    end do
      
    close(42)
    
  end subroutine read_response_diagonal
  
  
  
  !> Determine the index in the response matrix for a certain boundary degree of freedom.
  integer recursive function response_index(inode, i_starwall, ibas)
    
    ! --- Routine parameters
    integer, intent(in)    :: inode      !< Boundary index of the node
    integer, intent(in)    :: i_starwall !< STARWALL harmonic
    integer, intent(in)    :: ibas       !< Basis function (1 or 2)
    
    if ( (i_starwall < 0) .or. (i_starwall > sr%n_tor) ) then
      write(*,*) 'response_index: illegal value i_starwall=', i_starwall
      stop
    end if
    response_index = 2*sr%n_tor0*(inode-1) + 2*(i_starwall-1) + ibas
    
    if ( response_index < 1 ) then
      write(*,*) 'FATAL: RESPONSE_INDEX < 1 DETECTED'
      stop
    end if
    
  end function response_index
  
  
  
  !> Determine the index in the response matrix for a certain boundary degree of freedom.
  integer recursive function response_index_eq(inode, ibas)
    
    ! --- Routine parameters
    integer, intent(in)    :: inode      !< Boundary index of the node
    integer, intent(in)    :: ibas       !< Basis function (1 or 2)
    
    response_index_eq = 2*(inode-1) + ibas
    
    if ( response_index_eq < 1 ) then
      write(*,*) 'FATAL: RESPONSE_INDEX_EQ < 1 DETECTED'
      stop
    end if
    
  end function response_index_eq
  
  
  
  !> Return a description for a given toroidal mode index (mode number, cos/sin).
  character(len=12) function mode_to_str(i_tor, n_period)
    
    use phys_module, only: mode, mode_type
    
    implicit none
    
    ! --- Routine parameters
    integer, intent(in) :: i_tor     !< Toroidal mode index
    integer, intent(in) :: n_period  !< Periodicity
    
    ! --- Local variables
    integer           :: n                ! Toroidal mode number
    character(len=3)  :: typ              ! sin or cos
    character(len=30) :: i_tor_str, n_str ! Character string representations for i_tor and n
    

    write(i_tor_str,'(I10)') i_tor   ! toroidal mode index
    n   = mode(i_tor)                ! toroidal mode number
    write(n_str,'(I10)') n           !  -"-
    typ = mode_type(i_tor)           ! sin or cos
    
    mode_to_str = trim(adjustl(i_tor_str))//' (n='//trim(adjustl(n_str))//' '//trim(typ)//')'
    
  end function mode_to_str
  
  
  
  !> Return a description of several toroidal modes.
  character(len=1400) function modes_to_str(i_tors, n_tor, n_period)
  
    implicit none
    
    ! --- Routine parameters
    integer, intent(in) :: n_tor         !< Dimension of i_tors
    integer, intent(in) :: i_tors(n_tor) !< Toroidal mode numbers
    integer, intent(in) :: n_period      !< Periodicity
    
    ! --- Local variables
    integer :: i
    
    modes_to_str = ''
    do i = 1, n_tor
      modes_to_str = trim(modes_to_str)//' '//mode_to_str(i_tors(i), n_period)
      if ( i == n_tor ) exit
      modes_to_str = trim(modes_to_str)//', '
    end do
    
    modes_to_str = adjustl(modes_to_str)
    
  end function modes_to_str


  subroutine memory_prediction(n_w, nd_bez, ntasks)
  
    use mpi_mod
    implicit none

    ! --- Routine parameters
    integer, intent(in) :: n_w           !< Number of wall 
    integer, intent(in) :: nd_bez        !< Number of bezier elements
    integer, intent(in) :: ntasks        !< Number of MPI tasks

    ! --- Local variables
    real*8    :: tot_mem, gbyte_conv, real_n_w, real_nd_bez

    real_n_w = real(n_w,8); real_nd_bez = real(nd_bez,8) 

    100 format(80('-'))
    101 format(5x, a,'(',i6,',',i6,') = ',E11.3E2,' GB total / ',E11.3E2, ' GB per MPI task')
    102 format(10x,'TOTAL MEMORY CONSUMPTION PREDICTION = ', E12.3E2, ' GB total /', E11.3E2, 'GB per MPI task') 
    
    !1 DP converted to GB 
    gbyte_conv = 8.d0/1024.d0/1024.d0/1024.d0
 
    write(*,100)  
    write(*,*) 'Predicted memory consumption in the JOREK-STARWALL part:'
    write(*,*)
    write(*,101)'s_ww        ', n_w,      n_w,      real_n_w    * real_n_w    * gbyte_conv, real_n_w    * real_n_w    * gbyte_conv / ntasks
    write(*,101)'s_ww_inv    ', n_w,      n_w,      real_n_w    * real_n_w    * gbyte_conv, real_n_w    * real_n_w    * gbyte_conv / ntasks
    write(*,101)'a_ye        ', n_w,      nd_bez,   real_n_w    * real_nd_bez * gbyte_conv, real_n_w    * real_nd_bez * gbyte_conv / ntasks
    write(*,101)'a_ey        ', nd_bez,   n_w,      real_nd_bez * real_n_w    * gbyte_conv, real_nd_bez * real_n_w    * gbyte_conv / ntasks
    write(*,101)'a_ee        ', nd_bez,   nd_bez,   real_nd_bez * real_nd_bez * gbyte_conv, real_nd_bez * real_nd_bez * gbyte_conv / ntasks
    write(*,101)'a_id        ', nd_bez,   nd_bez,   real_nd_bez * real_nd_bez * gbyte_conv, real_nd_bez * real_nd_bez * gbyte_conv / ntasks
    write(*,101)'a_nw        ', nd_bez,   nd_bez,   real_nd_bez * real_nd_bez * gbyte_conv, real_nd_bez * real_nd_bez * gbyte_conv / ntasks
    write(*,101)'response_m_a', n_w,      nd_bez,   real_n_w    * real_nd_bez * gbyte_conv, real_n_w    * real_nd_bez * gbyte_conv / ntasks
    write(*,101)'response_m_d', n_w,      nd_bez,   real_n_w    * real_nd_bez * gbyte_conv, real_n_w    * real_nd_bez * gbyte_conv / ntasks
    write(*,101)'response_m_f', nd_bez,   n_w,      real_nd_bez * real_n_w    * gbyte_conv, real_nd_bez * real_n_w    * gbyte_conv / ntasks
    write(*,101)'response_m_g', nd_bez,   n_w,      real_nd_bez * real_n_w    * gbyte_conv, real_nd_bez * real_n_w    * gbyte_conv / ntasks
    write(*,101)'response_m_f', nd_bez,   n_w,      real_nd_bez * real_n_w    * gbyte_conv, real_nd_bez * real_n_w    * gbyte_conv / ntasks
    write(*,101)'response_m_h', nd_bez,   nd_bez,   real_nd_bez * real_nd_bez * gbyte_conv, real_nd_bez * real_nd_bez * gbyte_conv / ntasks
    write(*,101)'response_m_j', nd_bez,   nd_bez,   real_nd_bez * real_nd_bez * gbyte_conv, real_nd_bez * real_nd_bez * gbyte_conv / ntasks
    write(*,101)'response_m_e', nd_bez,   nd_bez,   real_nd_bez * real_nd_bez * gbyte_conv, real_nd_bez * real_nd_bez * gbyte_conv / ntasks

    tot_mem =  (2.d0*real_n_w**2 + 7.d0*real_n_w*real_nd_bez + 6.d0*real_nd_bez**2) * gbyte_conv

    write(*,*)
    write(*,102) tot_mem, tot_mem / ntasks
    write(*,100)

  end subroutine memory_prediction 

  subroutine distribute_coil_names
    use vacuum

    if ((sr%n_diag_coils .gt. 0) .and. (.not. allocated(diag_coil_name)) )  then
      allocate( diag_coil_name(sr%n_diag_coils) )
      diag_coil_name(:) = sr%coil_name(sr%ind_start_diag_coils:sr%ind_start_diag_coils +  sr%n_diag_coils -1)
    end if

    if (( sr%n_rmp_coils .gt. 0) .and. (.not. allocated( rmp_coil_name)) ) then
      allocate( rmp_coil_name(sr%n_rmp_coils) )
      rmp_coil_name(:)  = sr%coil_name(  sr%ind_start_rmp_coils:sr%ind_start_rmp_coils +  sr%n_rmp_coils  -1)
    end if
    
    if (( sr%n_pol_coils .gt. 0) .and. (.not. allocated(  pf_coil_name)) )  then
      allocate( pf_coil_name(sr%n_pol_coils) )
      pf_coil_name(:)   = sr%coil_name(  sr%ind_start_pol_coils:sr%ind_start_pol_coils +  sr%n_pol_coils  -1)
    end if
    
  end subroutine distribute_coil_names
  

end module vacuum_response

  
