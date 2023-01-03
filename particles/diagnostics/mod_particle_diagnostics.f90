!> This module contains some routines for calculating diagnostics on particles.
module mod_particle_diagnostics
use mod_io_actions
use mod_particle_sim
use mod_fields
use mod_particle_types
use hdf5
implicit none
private
public write_particle_diagnostics, calculate_particle_diagnostics

!> Cannot use HDF5 types here because these are invalid before h5open_f is called
!> (I think, did not take the chance)
integer, parameter :: REAL4 = 1, INT4 = 2, REAL8 = 3
integer, parameter :: n_vars = 15
character(len=7)  :: var_names(n_vars) = ["e      ", "k      ", "mu     ", &
  "psi_n  ", "psi_bar", "p_phi  ", "weight ", "lost   ", "q      ", "region ", &
  "theta  ", "phi    ", "R      ", "Z      ","i_elm  "]
integer, parameter :: var_types(n_vars) = [REAL8, REAL8, REAL4, REAL4, REAL4, REAL8, REAL4, INT4, INT4, INT4, REAL4, REAL4, REAL4, REAL4,INT4]
integer, parameter :: n_real8_var      = 3 !count(var_types .eq. REAL8)
integer, parameter :: n_real4_var      = 8 !count(var_types .eq. REAL4)
integer, parameter :: n_int4_var       = 4 !count(var_types .eq. INT4)
! HDF5 does not support booleans, use INT4

integer(HSIZE_T), parameter :: chunk_size(2) = [50000_HSIZE_T,1_HSIZE_T] !< particle, time


!> Action to calculate diagnostics and write this to an HDF5 file
!> in extensible (in the time-dimension) datasets
type, extends(io_action) :: write_particle_diagnostics
  integer(HID_T) :: file_id !< file identifier
  logical :: append = .false. !< Append to existing file
  logical, private :: init_done = .false. !< True after we have checked for an existing file
  integer, allocatable :: only(:) !< Do not output all variables (list of index into var_names)
  !< if not allocated skip it
contains
  procedure :: do => do_write_particle_diagnostics
end type write_particle_diagnostics
interface write_particle_diagnostics
  module procedure new_write_particle_diagnostics
end interface write_particle_diagnostics


contains
!> Constructor
function new_write_particle_diagnostics(filename, append, only) result(new)
  type(write_particle_diagnostics) :: new
  character(len=*), intent(in)     :: filename
  logical, intent(in), optional    :: append
  integer, intent(in), optional    :: only(:) !< if present do not write every variable
  new%filename = filename
  new%name = "WriteConstantsOfMotion"
  new%log = .true.
  if (present(append)) new%append = append
  if (present(only)) allocate(new%only(1:size(only)), source=only)
end function new_write_particle_diagnostics

!> Action to calculate all of these values and write them to an HDF5 file
subroutine do_write_particle_diagnostics(this, sim, ev)
  use hdf5_io_module
  use mpi
  use mod_event
  class(write_particle_diagnostics), intent(inout) :: this
  type(particle_sim), intent(inout)                :: sim
  type(event), intent(inout), optional             :: ev

  integer(HSIZE_T), parameter :: n_time = 1_HSIZE_T

  integer(HSIZE_T)  :: data_dims(2), time_dims(1), data_maxdims(2), time_maxdims(1)
  integer           :: i, j, n, my_id, n_cpu, ierr, dot_index, var_index
  integer(HID_T)    :: dspace, dset, mem_space
  integer(HID_T)    :: tspace, t_mem_space, tset, group_id, plist
  logical           :: link_exists, file_exists
  character(len=80) :: dataset_name, timeset_name, group_name
  character(len=123):: filename !< len=120 from this%filename + 3
  integer, dimension(:), allocatable           :: particles_per_proc
  real*8, dimension(:,:), allocatable, target  :: real8_var, real8_var_all !< Data storage order: particle index, var
  real*4,  dimension(:,:), allocatable, target :: real4_var, real4_var_all !< Data storage order: particle index, var
  integer, dimension(:,:), allocatable, target :: int4_var, int4_var_all !< Data storage order: particle index, var
  real*4 :: tmpreal4, times(1)
  integer :: tmpint4
  real*8 :: tmpreal8

  ! For psi_Axis, psi_sep
  real*8 :: psi_val, R, Z, s, t, psi_axis, psi_xpoint
  integer :: i_elm, ifail, xcase

  call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)      ! id of each MPI proc
  call MPI_COMM_SIZE(MPI_COMM_WORLD, n_cpu, ierr)      ! number of MPI procs
  if (my_id .eq. 0) call h5open_f(ierr)

  ! Move an old file if this is the initialisation, we are id 0 and we are not appending
  if (.not. this%append .and. my_id .eq. 0 .and. .not. this%init_done) then
    this%init_done = .true. ! only on ID 0, watch out
    ! Move old file if present
    inquire(file=trim(this%filename), exist=file_exists)
    dot_index = scan(this%filename,'.',back=.true.)
    if (file_exists) then
      do i=1,99
        write(filename,'(A,i0.2,A)') this%filename(1:dot_index-1), i, this%filename(dot_index:)
        inquire(file=trim(filename), exist=file_exists)
        if (.not. file_exists) then
          call system("mv '"//trim(this%filename)//"' '"//trim(filename)//"'")
          exit ! loop
        end if
      end do
    end if
  end if
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
      
  if (my_id .eq. 0) then
    call HDF5_open_or_create(this%filename, file_id=this%file_id, ierr=ierr)
    if (ierr .ne. 0) then
      write(*,*) "ERROR: cannot open HDF5 file"
      call MPI_Abort(MPI_COMM_WORLD, -1, ierr)
    end if

    ! Create group to write particle groups in if it does not exist yet
    call h5lexists_f(this%file_id, "/groups", link_exists, ierr)
    if (.not. link_exists) then
      call h5gcreate_f(this%file_id, "/groups", group_id, ierr)
      call h5gclose_f(group_id, ierr)
    end if
    ! assume that if it exists it's a group
  end if


  ! Safety checks
  if (.not. allocated(sim%groups)) return

  ! Preparation
  allocate(particles_per_proc(0:n_cpu-1))
  psi_axis = -1.d0
  psi_xpoint = 0.d0

  ! For each of the groups
  do i=lbound(sim%groups,1),ubound(sim%groups,1)
    ! Find the number of particles on each node
    n = size(sim%groups(i)%particles,1)
    call MPI_Gather(n,1,MPI_INTEGER,&
        particles_per_proc,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    allocate(real8_var(n,n_real8_var))
    allocate(real4_var(n,n_real4_var))
    allocate(int4_var( n,n_int4_var))

    ! Create group to write particles in if it does not exist yet
    if (my_id .eq. 0) then
      write(group_name,'(A,i0.3)') 'groups/', i
      call h5lexists_f(this%file_id, group_name, link_exists, ierr)
      if (.not. link_exists) then
        call h5gcreate_f(this%file_id, group_name, group_id, ierr)
        call h5gclose_f(group_id, ierr)
      end if
      ! assume that if it exists it's a group. This breaks appending to old files, not so bad


      ! Check the timeset existence and properties
      write(timeset_name,'(A,i0.3,A)') 'groups/', i, '/t'
      call h5lexists_f(this%file_id, trim(timeset_name), link_exists, ierr)

      if (link_exists) then
        call h5dopen_f(this%file_id, trim(timeset_name), tset, ierr)
        if (ierr .ne. 0) then
          write(*,*) "Error opening timeset", i
          call MPI_Abort(MPI_COMM_WORLD, -1, ierr)
        else
          call h5dget_space_f(tset, tspace, ierr)
        end if
      else
        call create_constants_time_dataset(this%file_id, trim(timeset_name), &
            tset, tspace)
      end if

      ! Get the current time dimensions
      call h5sget_simple_extent_dims_f(tspace, time_dims, time_maxdims, ierr)
      ! Check that the current time is > the last stored time
      if (time_dims(1) .ge. 1) then
        call HDF5_array1D_reading_r4(this%file_id,times,trim(timeset_name),[time_dims(1)-1])
        if (times(1) .gt. sim%time) then
          write(*,*) "Current time smaller than last stored diagnostic time, aborting."
          call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
        end if
      end if

      ! Extend the dataset by 1 in the time-dimension
      ! After extending, tspace is invalid. Close here already
      call h5sclose_f(tspace, ierr)
      time_dims(1) = time_dims(1) + 1_HSIZE_T
      call h5dset_extent_f(tset, time_dims, ierr)
    end if


    ! Generate data
    call calculate_particle_diagnostics(sim%fields, sim%time, sim%groups(i)%particles, sim%groups(i)%mass, &
        real8_var(:,:), real4_var(:,:), int4_var(:,:), psi_axis, psi_xpoint, dt=sim%groups(i)%dt)

    ! Gather to root
    if (my_id .eq. 0) then
      allocate(real8_var_all(sum(particles_per_proc,1),n_real8_var))
      allocate(real4_var_all(sum(particles_per_proc,1),n_real4_var))
      allocate(int4_var_all( sum(particles_per_proc,1),n_int4_var))
    else
      ! dummy allocations to make the compiler happy
      allocate(real8_var_all(1,n_real8_var))
      allocate(real4_var_all(1,n_real4_var))
      allocate(int4_var_all( 1,n_int4_var))
    end if
    do j=1,n_real8_var
      call MPI_Gatherv(real8_var(:,j), size(real8_var,1), MPI_REAL8, &
        real8_var_all(:,j), particles_per_proc, [(sum(particles_per_proc(0:i-1),1), i=0,n_cpu-1)], &
        MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    end do
    do j=1,n_real4_var
      call MPI_Gatherv(real4_var(:,j), size(real4_var,1), MPI_REAL4, &
        real4_var_all(:,j), particles_per_proc, [(sum(particles_per_proc(0:i-1),1), i=0,n_cpu-1)], &
        MPI_REAL4, 0, MPI_COMM_WORLD, ierr)
    end do
    do j=1,n_int4_var
      call MPI_Gatherv(int4_var(:,j), size(int4_var,1), MPI_INTEGER, &
        int4_var_all(:,j), particles_per_proc, [(sum(particles_per_proc(0:i-1),1), i=0,n_cpu-1)], &
        MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    end do

    if (my_id .eq. 0) then
      ! Loop over variables, write
      do j=1,n_vars
        if (allocated(this%only)) then
          if (.not. any(this%only .eq. j)) cycle ! skip this variable
        end if

        write(dataset_name,'(A,i0.3,A,A)') 'groups/', i, "/", trim(var_names(j))
        call h5lexists_f(this%file_id, trim(dataset_name), link_exists, ierr)
        if (link_exists) then
          call h5dopen_f(this%file_id, trim(dataset_name), dset, ierr)
          if (ierr .ne. 0) then
            write(*,*) "Error opening dataset", i
            call MPI_Abort(MPI_COMM_WORLD, -1, ierr)
          else
            call h5dget_space_f(dset, dspace, ierr)
          end if
        else
          ! Create a dataspace with unlimited dimensions, 
          call h5screate_simple_f(2, [int(sum(particles_per_proc,1),kind=HSIZE_T),0_HSIZE_T], dspace, ierr, &
              maxdims=[H5S_UNLIMITED_F,H5S_UNLIMITED_F])
          ! Create a dataset property list, enable chunking
          call h5pcreate_f(H5P_DATASET_CREATE_F, plist, ierr)
          call h5pset_chunk_f(plist, 2, chunk_size, ierr)
          select case (var_types(j))
            case (REAL8)
              call h5dcreate_f(this%file_id, dataset_name, H5T_IEEE_F64LE, dspace, dset, ierr, plist)
            case (REAL4)
              call h5dcreate_f(this%file_id, dataset_name, H5T_IEEE_F32LE, dspace, dset, ierr, plist)
            case (INT4)
              call h5dcreate_f(this%file_id, dataset_name, H5T_STD_I32LE, dspace, dset, ierr, plist)
            case DEFAULT
              write(*,*) "Unknown variable type for diagnostics"
              call MPI_ABORT(MPI_COMM_WORLD, -1, ierr)
          end select
          call h5pclose_f(plist, ierr)
        end if

        ! Now we must add one to the time dimension for this variable
        call h5sget_simple_extent_dims_f(dspace, data_dims, data_maxdims, ierr)
        data_dims(1) = max(data_dims(1), int(sum(particles_per_proc,1),HSIZE_T)) ! only grow this set
        data_dims(2) = data_dims(2) + 1_HSIZE_T
        if (time_dims(1) .ne. data_dims(2)) then
          write(*,*) "ERROR: data and time series are not of equal length"
          call MPI_Abort(MPI_COMM_WORLD, -1, ierr)
        end if
        call h5dset_extent_f(dset, data_dims, ierr) ! dspace becomes invalid now

        ! Create memory space
        call h5screate_simple_f(2, [int(sum(particles_per_proc,1),HSIZE_T), n_time], mem_space, ierr)
        ! We must create a hyperslab space here because simple does not play well
        ! with chunked HDF5 (but does not tell you this!).
        ! It is also important to get the space belonging to this dataset, instead of
        ! creating a new one.
        call h5dget_space_f(dset, dspace, ierr)
        call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, &
            start=[0_HSIZE_T, data_dims(2)-1_HSIZE_T], &
            count=[1_HSIZE_T,1_HSIZE_T], &
            block=[int(sum(particles_per_proc,1),HSIZE_T),n_time], &
            hdferr=ierr)

        var_index = count(var_types(1:j) .eq. var_types(j))
        select case (var_types(j))
          case (REAL8)
            call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, real8_var_all(:,var_index:var_index), [1_HSIZE_T,1_HSIZE_T], &
                 ierr, file_space_id = dspace, mem_space_id = mem_space)
          case (REAL4)
            call h5dwrite_f(dset, H5T_NATIVE_REAL, real4_var_all(:,var_index:var_index), [1_HSIZE_T,1_HSIZE_T], &
                 ierr, file_space_id = dspace, mem_space_id = mem_space)
          case (INT4)
            call h5dwrite_f(dset, H5T_NATIVE_INTEGER, int4_var_all(:,var_index:var_index), [1_HSIZE_T,1_HSIZE_T], &
                 ierr, file_space_id = dspace, mem_space_id = mem_space)
          case DEFAULT
            write(*,*) "Unknown variable type for diagnostics"
            call MPI_ABORT(MPI_COMM_WORLD, -1, ierr)
        end select

        call h5sclose_f(mem_space, ierr)
        call h5sclose_f(dspace, ierr)
        call h5dclose_f(dset, ierr)
      end do

      ! Add the current time to the timeset
      call h5dget_space_f(tset, tspace, ierr)
      call h5sselect_hyperslab_f(tspace, H5S_SELECT_SET_F, &
          start=[time_dims(1)-1_HSIZE_T], count=[n_time], hdferr=ierr)
      call h5screate_f(H5S_SCALAR_F, t_mem_space, ierr)
      call h5dwrite_f(tset, H5T_NATIVE_REAL, real(sim%time,4), [n_time], ierr, file_space_id=tspace, mem_space_id=t_mem_space)


      call h5sclose_f(t_mem_space, ierr)
      call h5sclose_f(tspace, ierr)
      call h5dclose_f(tset, ierr)
    end if
    deallocate(real8_var_all,real4_var_all,int4_var_all)
    deallocate(real8_var,real4_var,int4_var)
  end do


  ! Write psi_axis and psi_sep into the file
  if (my_id .eq. 0) then
    do i=1,2
      if (i .eq. 1) then
        timeset_name = 'psi_axis'
        psi_val = psi_axis
      else
        timeset_name = 'psi_sep'
        psi_val = psi_xpoint
      end if
      call h5lexists_f(this%file_id, trim(timeset_name), link_exists, ierr)
      if (link_exists) then
        call h5dopen_f(this%file_id, trim(timeset_name), tset, ierr)
        if (ierr .ne. 0) then
          write(*,*) "Error opening ", trim(timeset_name), i
          call MPI_Abort(MPI_COMM_WORLD, -1, ierr)
        else
          call h5dget_space_f(tset, tspace, ierr)
        end if
      else
        call create_constants_time_dataset(this%file_id, trim(timeset_name), &
            tset, tspace)
      end if
      ! Time_dims is set from earlier
      call h5dset_extent_f(tset, time_dims, ierr)
      ! close the tspace and make a new one after setting extent above
      call h5sclose_f(tspace, ierr)
      ! Add the current value to the timeset
      call h5dget_space_f(tset, tspace, ierr)
      call h5sselect_hyperslab_f(tspace, H5S_SELECT_SET_F, &
          start=[time_dims(1)-1_HSIZE_T], count=[n_time], hdferr=ierr)
      call h5screate_f(H5S_SCALAR_F, t_mem_space, ierr)
      call h5dwrite_f(tset, H5T_NATIVE_REAL, real(psi_val,4), [n_time], ierr, file_space_id=tspace, mem_space_id=t_mem_space)
      call h5dclose_f(tset, ierr)
      call h5sclose_f(t_mem_space, ierr)
      call h5sclose_f(tspace, ierr)
    end do


    call h5fclose_f(this%file_id, ierr)
    call h5close_f(ierr)
  end if
end subroutine do_write_particle_diagnostics

!> Create a new dataset for time data (1-d extensible, chunked (required for extensibility))
subroutine create_constants_time_dataset(file_id, dataset_name, dset, dspace)
  integer(HID_T), intent(in)   :: file_id
  character(len=*), intent(in) :: dataset_name
  integer(HID_T), intent(out)  :: dset, dspace
  integer :: ierr
  integer(HID_T) :: crp_list
  integer(HSIZE_T), parameter :: chunk_size(1) = [1000]

  ! Create a dataspace with unlimited dimensions, 
  call h5screate_simple_f(1, [0_HSIZE_T], dspace, ierr, &
      maxdims=[H5S_UNLIMITED_F])
  ! Create a dataset property list, enable chunking
  call h5pcreate_f(H5P_DATASET_CREATE_F, crp_list, ierr)
  call h5pset_chunk_f(crp_list, 1, chunk_size, ierr)
  call h5dcreate_f(file_id, dataset_name, H5T_NATIVE_REAL, dspace, dset, ierr, crp_list)
  call h5pclose_f(crp_list, ierr)
end subroutine create_constants_time_dataset


!> Calculate H, mu, P_phi, Psi_bar, Psi, q, weight for a list of particles.
!> mask is .f. if a particle is lost. Those values in out are set to 0.d0
subroutine calculate_particle_diagnostics(fields, time, particles, mass, real8_stats, real4_stats, int_stats, psi_axis, psi_limit, mask, dt)
  use mod_particle_types
  use phys_module, only: F0, xpoint, xcase
  use constants
  use mod_boris
  use mod_kinetic_relativistic, only: relativistic_kinetic_to_particle, relativistic_kinetic_to_gc
  use mod_gc_relativistic, only: relativistic_gc_to_particle, relativistic_gc_to_gc
  use mod_fields_linear
  use domains
  use equil_info
  class(fields_base), intent(in)                               :: fields
  real*8, intent(in)                                           :: time
  class(particle_base), intent(in), dimension(:)               :: particles
  real*8, intent(in)                                           :: mass
  real*8, dimension(:,:), intent(out)                          :: real8_stats !< List of values (H, mu, P_phi, Rho, weight, theta, phi)
  real*4, dimension(:,:), intent(out)                          :: real4_stats !< List of values (H, mu, P_phi, Rho, weight, theta, phi)
  integer, dimension(:,:), intent(out)                         :: int_stats !< List of values (q, lost, region)
  real*8, intent(out)                                          :: psi_axis, psi_limit
  logical, dimension(:), intent(out), optional                 :: mask !< Mask containing .f. if particle is lost
  real*8, intent(in), optional                                 :: dt !< for leapfrog types perform a half step to get the correct velocity
  real*8                 :: psi, U, E(3), B(3), v_par
  type(particle_gc)      :: particle
  type(particle_kinetic) :: particle_centered ! a particle where position and velocity are known at same time

  integer :: i, ifail
  integer :: domain, i_elm_axis, i_elm_xpoint(2)
  real*8  :: R_axis, Z_axis, s_axis, t_axis
  real*8, dimension(2) :: psi_xpoint, R_xpoint, Z_xpoint, s_xpoint, t_xpoint

  real*8 :: real_stats_tmp(n_real8_var+n_real4_var) !< Temporary storage for real8 and real4 stats
  integer :: i_real8, i_real4, i_tmp, j

  ! Preparation (force my_id to 1 to suppress message)
  call find_axis(1,fields%node_list,fields%element_list,psi_axis,R_axis,Z_axis,i_elm_axis,s_axis,t_axis,ifail)

  if (xpoint) then
    call find_xpoint(1,fields%node_list,fields%element_list,psi_xpoint,R_xpoint,Z_xpoint,i_elm_xpoint,s_xpoint,t_xpoint,xcase,ifail)
    psi_limit  = psi_xpoint(1)
    if(ES%active_xpoint .eq. UPPER_XPOINT) then
      psi_limit = psi_xpoint(2)
    endif
  else
    write(*,*) "WARNING: particle diagnostics with limiter not fully supported yet! Setting psi_limit to 0"
    psi_limit = 0.d0
  endif

  ! Call which_domain once to setup saved values
  domain = which_domain(fields%node_list, fields%element_list, &
      0.d0, 0.d0, &
      0.d0, xpoint, xcase, R_xpoint, Z_xpoint, psi_xpoint, psi_limit, &
      R_axis, Z_axis, psi_axis)


  if (present(mask)) mask = .true.
  real8_stats = 0.d0
  real4_stats = 0.d0
  int_stats  = 0
#ifndef __NVCOMPILER  
#ifdef __GFORTRAN__
  !$omp parallel do default(shared) &
#else
  !$omp parallel do default(none) &
  !$omp shared(particles, fields, int_stats, real4_stats, real8_stats, mask, time, f0, mass, dt, &
  !$omp        xpoint, xcase, R_xpoint, Z_xpoint, psi_xpoint, psi_limit, R_axis, Z_axis, psi_axis) &
#endif
  !$omp private(E, B, psi, U, particle, v_par, domain, particle_centered, real_stats_tmp, i_real8, i_real4, i_tmp, j)
#endif
  do i=1,size(particles,1)
    int_stats(i,4) = particles(i)%i_elm
    if (particles(i)%i_elm .lt. 1) then
      if (present(mask)) mask(i) = .false.
      int_stats(i,1) = 1 ! lost
      select type (particle_in => particles(i))
      type is (particle_kinetic_leapfrog)
        int_stats(i,2) = particle_in%q
      type is (particle_kinetic)
        int_stats(i,2) = particle_in%q
      type is (particle_gc)
        int_stats(i,2) = particle_in%q
      type is (particle_kinetic_relativistic)
        int_stats(i,2) = particle_in%q
     type is (particle_gc_relativistic)
        int_stats(i,2) = particle_in%q
      end select
    else
      call fields%calc_EBpsiU(time, particles(i)%i_elm, particles(i)%st, particles(i)%x(3), E, B, psi, U)
      ! Calculate psi and B in the current particle location (either GC or kinetic)

      select type (particle_in => particles(i))
      type is (particle_kinetic_leapfrog)
        ! Let the conversion function calculate the conserved quantities and a new position
        if (present(dt)) then
          particle_centered = kinetic_leapfrog_to_kinetic(particle_in, E, B, mass, dt)
        else
          particle_centered = kinetic_leapfrog_to_kinetic(particle_in, E, B, mass, 0.d0) ! does no real work
        end if
        particle = kinetic_to_gc(fields%node_list, fields%element_list, particle_centered, B, mass)

        if (particle%i_elm .le. 0) cycle

        ! P_phi (generalized canonical toroidal momentum) at the particle position
        real_stats_tmp(6) = particle_centered%q * EL_CHG * psi + mass * ATOMIC_MASS_UNIT * particle_centered%x(1) * particle_centered%v(3)
      type is (particle_gc)
        v_par    = sign(sqrt(2*(particle%E-particle%mu*norm2(B))*EL_CHG/(mass*ATOMIC_MASS_UNIT)),particle%mu)
        particle = particle_in
        real_stats_tmp(6) = real(particle%q,8) * EL_CHG * psi + mass * ATOMIC_MASS_UNIT * particle%x(1) * v_par * B(3)/norm2(B)
      type is (particle_fieldline)
        particle = particle_in
        real_stats_tmp(6) = 0.d0 ! Since there is no momentum defined for this we just use 0
      type is (particle_kinetic_relativistic)
  ! compute the canonical toroidal momentum P_phi
         real_stats_tmp(6) = real(particle_in%q,8)*EL_CHG*psi - ATOMIC_MASS_UNIT*particle_in%x(1)* &
                                 (particle_in%p(1)*sin(particle_in%x(3))+particle_in%p(2)*cos(particle_in%x(3)))
	! transform the particle into a gc to get E and mu
        call relativistic_kinetic_to_particle(fields%node_list,fields%element_list, particle_in, particle, mass, B)
      type is (particle_gc_relativistic)
  ! compute the canonical toroidal momentum P_phi
        real_stats_tmp(6) = EL_CHG*particle_in%q*psi + ATOMIC_MASS_UNIT*particle_in%x(1)* particle_in%p(1)*B(3)/norm2(B)
  ! transform the particle into a gc to get E and mu
        call relativistic_gc_to_particle(fields%node_list,fields%element_list, particle_in,particle,mass,B) 
      class default
        write(*,*) "ERROR: calculate_particle_diagnostics not implemented for this particle type"
        cycle ! skip this iteration
      end select


      ! Calculate output variables (see var_names on top of the module for ordering)

      ! Total energy (including electric potential at this charge state) at kinetic or GC position
      real_stats_tmp(1) = particle%E + particle%q * F0 * U ! E in [eV] + q [e] * U [V]
      ! Kinetic energy
      real_stats_tmp(2) = particle%E    
      ! Magnetic moment
      real_stats_tmp(3) = particle%mu
      ! Psi_n (normalized psi, at kinetic position)
      real_stats_tmp(4) = (psi-psi_axis)/(psi_limit-psi_axis)
      ! P_phi (generalized canonical toroidal momentum) (see above
      ! skipped, see above
      ! real_stats(i,6) = ...
      ! Psi_bar = P_phi/Ze
      if (particle%q .ne. 0) real_stats_tmp(5) = real_stats_tmp(6) / (particle%q * EL_CHG)
      ! weight
      real_stats_tmp(7) = real(particle%weight, 8)
      ! theta
      real_stats_tmp(8) = atan2(particles(i)%x(2)-Z_axis, particles(i)%x(1)-R_axis)
      ! phi
      real_stats_tmp(9) = particles(i)%x(3)
      ! R
      real_stats_tmp(10) = particles(i)%x(1)
      ! Z
      real_stats_tmp(11) = particles(i)%x(2)           

      ! 1. lost (boolean)
      int_stats(i,1) = 0
      ! 2. q (charge)
      int_stats(i,2) = particle%q
      ! 3. region (enum)
      domain = which_domain(fields%node_list, fields%element_list, &
          particle%x(1), particle%x(2), &
          psi, xpoint, xcase, R_xpoint, Z_xpoint, psi_xpoint, psi_limit, &
          R_axis, Z_axis, psi_axis)
      int_stats(i,3) = domain

      ! Distribute output to real4 or real8 as needed
      i_real8 = 0
      i_real4 = 0
      i_tmp = 0
      do j=1,n_vars
        if (var_types(j) .eq. REAL8) then
          i_tmp = i_tmp + 1
          i_real8 = i_real8 + 1
          real8_stats(i,i_real8) = real_stats_tmp(i_tmp)
        else if (var_types(j) .eq. REAL4) then
          i_tmp = i_tmp + 1
          i_real4 = i_real4 + 1
          real4_stats(i,i_real4) = real(real_stats_tmp(i_tmp), 4)
        end if
      end do
    endif
  enddo
#ifndef __NVCOMPILER  
  !$omp end parallel do
#endif
end subroutine calculate_particle_diagnostics


!> Calculate particles present in specific regions on all particles
!> Performs MPI communications to sum values, returns the value
!> corresponding to all particles on node 0, and the value for each node on this node
!> Regions are: DOMAIN_PLASMA, DOMAIN_SOL, DOMAIN_OUTER_SOL,
!> DOMAIN_UPPER_PRIVATE, DOMAIN_LOWER_PRIVATE
function particles_in_regions(node_list, element_list, particles)
  use constants
  use data_structure
  use phys_module, only: DOMAIN_PLASMA, DOMAIN_SOL, DOMAIN_OUTER_SOL, DOMAIN_UPPER_PRIVATE,        &
      DOMAIN_LOWER_PRIVATE, xpoint, xcase
  use mod_particle_types
  use domains
  use mpi
  use mod_interp, only: interp
  use equil_info
  implicit none

  type(type_node_list), intent(in)     :: node_list
  type(type_element_list), intent(in)  :: element_list
  class(particle_base), intent(in), dimension(:) :: particles

  integer, dimension(DOMAIN_PLASMA:DOMAIN_LOWER_PRIVATE) :: particles_in_regions, tmp
  integer :: i, ifail, my_id
  integer :: domain, i_elm_axis, i_elm_xpoint(2)
  real*8  :: psi, psi_s, psi_t, psi_st, psi_ss, psi_tt
  real*8  :: psi_axis, R_axis, Z_axis, s_axis, t_axis, psi_limit
  real*8, dimension(2) :: psi_xpoint, R_xpoint, Z_xpoint, s_xpoint, t_xpoint


  !! Preparation (force my_id to 1 to suppress message)
  call find_axis(1,node_list,element_list,psi_axis,R_axis,Z_axis,i_elm_axis,s_axis,t_axis,ifail)

  if (xpoint) then
    call find_xpoint(1,node_list,element_list,psi_xpoint,R_xpoint,Z_xpoint,i_elm_xpoint,s_xpoint,t_xpoint,xcase,ifail)
    psi_limit  = psi_xpoint(1)
    if(ES%active_xpoint .eq. UPPER_XPOINT) then
      psi_limit = psi_xpoint(2)
    endif
  else
    psi_limit = 0.d0
  endif

  ! Call which_domain once to setup saved values
  domain = which_domain(node_list, element_list, &
      0.d0, 0.d0, &
      0.d0, xpoint, xcase, R_xpoint, Z_xpoint, psi_xpoint, psi_limit, &
      R_axis, Z_axis, psi_axis)

  tmp = 0
#ifdef __GFORTRAN__
  !$omp parallel do default(shared) & 
#else
  !$omp parallel do default(none)   &
  !$omp shared(node_list, element_list, particles, xpoint, xcase, R_xpoint, Z_xpoint, psi_xpoint, psi_limit, &
  !$omp        R_axis, Z_axis, psi_axis) &
#endif
  !$omp private(domain, psi, psi_s, psi_t, psi_st, psi_ss, psi_tt) &
  !$omp reduction(+:tmp)
  do i=1,size(particles,1)
    if (particles(i)%i_elm .le. 0 .or. particles(i)%i_elm .gt. element_list%n_elements) cycle
    call interp(node_list, element_list, particles(i)%i_elm, 1, 1, &                               ! force i_harm to 1
                particles(i)%st(1), particles(i)%st(2), psi, psi_s, psi_t, psi_st, psi_ss, psi_tt)

    domain = which_domain(node_list, element_list, &
        particles(i)%x(1), particles(i)%x(2), &
        psi, xpoint, xcase, R_xpoint, Z_xpoint, psi_xpoint, psi_limit, &
        R_axis, Z_axis, psi_axis)

    tmp(domain) = tmp(domain) + 1
  end do
  !$omp end parallel do

  ! Save values on nodes
  particles_in_regions = tmp
  ! Mpi communication to get the total answer on node 0
  call MPI_Reduce(tmp, particles_in_regions, DOMAIN_LOWER_PRIVATE-DOMAIN_PLASMA, &
    MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ifail)
  call MPI_Comm_Rank(MPI_COMM_WORLD, my_id, ifail)
end function particles_in_regions
end module mod_particle_diagnostics
