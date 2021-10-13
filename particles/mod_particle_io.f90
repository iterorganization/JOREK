!> Particle input-output module, containing hdf5 data_type and writing routines
!> TODO: add metadata and/or use H5MD format (http://nongnu.org/h5md/h5md.html)
module mod_particle_io
use iso_c_binding
use hdf5_io_module
use hdf5
use mpi
use mod_interp, only: interp_0
use mod_particle_types
use mod_particle_sim
implicit none
private
public write_simulation_hdf5, read_simulation_hdf5, get_simulation_hdf5_time

integer(HSIZE_T), parameter :: particle_type_name_length = 40 !< length of the string used to identify a specific type of particle
character(len=*), parameter :: particle_type_name_field_name = 'particle_type' !< Name of the field containing the particle_type_name
contains

!> Export all particles using HDF5 IO
!> This communicates the particles internally before writing since this gave
!> better performance in our tests. A parallel version is in the git history.
!> Note that reading is still performed in parallel as this is not considered
!> time-critical and simpler to write.
subroutine write_simulation_hdf5(sim, filename)
type(particle_sim)   , intent(in) :: sim
character*(*)        , intent(in) :: filename

integer :: my_id, n_cpu, ierr
integer :: n_here, n_total
integer(HSIZE_T) :: i_here
integer, allocatable, dimension(:) :: particles_per_proc

! For HDF5 writing
integer(HID_T)                :: file, create_file_space, write_file_space, dset, plist ! handles
integer(HID_T)                :: group_id
integer(HID_T)                :: data_type
integer(HID_T)                :: time_space_id, time_set_id
character(len=12)             :: group_name
character(len=particle_type_name_length) :: particle_type_name
integer                       :: i, j, hdferr
type(c_ptr) :: p_ptr
real*8, dimension(:,:), allocatable :: x, v, x_all, v_all, st, st_all
real*8, dimension(:), allocatable   :: Vpar, E, mu, v1, Vpar_all, E_all, mu_all, v1_all
real*4, dimension(:), allocatable   :: weight, weight_all, t_birth, t_birth_all
integer, dimension(:), allocatable  :: i_elm, i_elm_all, i_life, i_life_all
integer, dimension(:), allocatable  :: q, q_all, lost, lost_all

! Preparation
call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)      ! id of each MPI proc
call MPI_COMM_SIZE(MPI_COMM_WORLD, n_cpu, ierr)      ! number of MPI procs
call h5open_f(hdferr)
allocate(particles_per_proc(0:n_cpu-1))
particles_per_proc = 0


if (my_id .eq. 0) then
  ! Create file
  call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file, hdferr)
  if (hdferr .gt. 0) then
    write(*,*) "file open failed:", hdferr
    call MPI_Abort(MPI_COMM_WORLD, -1, ierr)
  endif

  ! Create group to write particle groups in
  call h5gcreate_f(file, "/groups", group_id, hdferr)
  call h5gclose_f(group_id, hdferr)

  ! Write the time
  call HDF5_real_saving(file,sim%time,'/time')
end if


if (allocated(sim%groups)) then
  do i=1,size(sim%groups,1)
    if (.not. allocated(sim%groups(i)%particles)) then
      write(*,*) "WARNING: group ", i, " not allocated, exiting"
      call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
    end if
    ! Find the number of particles on each node
    n_here = size(sim%groups(i)%particles,1)
    call MPI_Gather(n_here,1,MPI_INTEGER,particles_per_proc,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    n_total = sum(particles_per_proc,1) ! set it on other processes as well (to
    ! 0) so they can allocate the useless arrays

    if (my_id .eq. 0) then
      ! Create group to write in
      write(group_name,"(A,i0.3,A)") "/groups/", i, "/"
      call h5gcreate_f(file, group_name, group_id, hdferr)
      call h5gclose_f(group_id, hdferr)
    end if

    ! particle_base properties
    ! x
    allocate(x(3,n_here), x_all(3,n_total))
    do j=1,n_here
      x(:,j) = sim%groups(i)%particles(j)%x
    end do
!    call MPI_Gatherv(x(:,:), 3*n_here, MPI_REAL8, &
!      x_all(:,:), particles_per_proc*3, [(sum(particles_per_proc(1:i),1)*3, i=0,n_cpu-1)], &
!      MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

    call MPI_Gatherv(x(1,:), n_here, MPI_REAL8, &
                     x_all(1,:), particles_per_proc, [(sum(particles_per_proc(1:i),1), i=0,n_cpu-1)], &
                     MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    call MPI_Gatherv(x(2,:), n_here, MPI_REAL8, &
                     x_all(2,:), particles_per_proc, [(sum(particles_per_proc(1:i),1), i=0,n_cpu-1)], &
                     MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    call MPI_Gatherv(x(3,:), n_here, MPI_REAL8, &
                     x_all(3,:), particles_per_proc, [(sum(particles_per_proc(1:i),1), i=0,n_cpu-1)], &
                     MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

    ! st
    allocate(st(2,n_here), st_all(2,n_total))
    do j=1,n_here
      st(:,j) = sim%groups(i)%particles(j)%st
    end do
    call MPI_Gatherv(st(:,:), 2*n_here, MPI_REAL8, &
      st_all(:,:), particles_per_proc*2, [(sum(particles_per_proc(1:i),1)*2, i=0,n_cpu-1)], &
      MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

    ! weight
    allocate(weight(n_here), weight_all(n_total))
    do j=1,n_here
      weight(j) = sim%groups(i)%particles(j)%weight
    end do
    call MPI_Gatherv(weight(:), n_here, MPI_REAL4, &
      weight_all(:), particles_per_proc, [(sum(particles_per_proc(1:i),1), i=0,n_cpu-1)], &
      MPI_REAL4, 0, MPI_COMM_WORLD, ierr)

    ! i_elm
    allocate(i_elm(n_here), i_elm_all(n_total))
    do j=1,n_here
      i_elm(j) = sim%groups(i)%particles(j)%i_elm
    end do
    call MPI_Gatherv(i_elm(:), n_here, MPI_INTEGER, &
      i_elm_all(:), particles_per_proc, [(sum(particles_per_proc(1:i),1), i=0,n_cpu-1)], &
      MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! i_life
    allocate(i_life(n_here), i_life_all(n_total))
    do j=1,n_here
      i_life(j) = sim%groups(i)%particles(j)%i_life
    end do
    call MPI_Gatherv(i_life(:), n_here, MPI_INTEGER, &
      i_life_all(:), particles_per_proc, [(sum(particles_per_proc(1:i),1), i=0,n_cpu-1)], &
      MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! t_birth
    allocate(t_birth(n_here), t_birth_all(n_total))
    do j=1,n_here
      t_birth(j) = sim%groups(i)%particles(j)%t_birth
    end do
    call MPI_Gatherv(t_birth(:), n_here, MPI_REAL4, &
      t_birth_all(:), particles_per_proc, [(sum(particles_per_proc(1:i),1), i=0,n_cpu-1)], &
      MPI_REAL4, 0, MPI_COMM_WORLD, ierr)

    ! Write out stuff depending on particle type
    select type (p => sim%groups(i)%particles)

    type is (particle_kinetic)
      particle_type_name = 'particle_kinetic'

      ! v
      allocate(v(3,n_here), v_all(3,n_total))
      do j=1,n_here
        v(:,j) = p(j)%v
      end do
!      call MPI_Gatherv(v(:,:), 3*n_here, MPI_REAL8, &
!        v_all(:,:), particles_per_proc*3, [(sum(particles_per_proc(1:i),1)*3, i=0,n_cpu-1)], &
!        MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

      call MPI_Gatherv(v(1,:), n_here, MPI_REAL8, &
                       v_all(1,:), particles_per_proc, [(sum(particles_per_proc(1:i),1), i=0,n_cpu-1)], &
                       MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
      call MPI_Gatherv(v(2,:), n_here, MPI_REAL8, &
                       v_all(2,:), particles_per_proc, [(sum(particles_per_proc(1:i),1), i=0,n_cpu-1)], &
                       MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
      call MPI_Gatherv(v(3,:), n_here, MPI_REAL8, &
                       v_all(3,:), particles_per_proc, [(sum(particles_per_proc(1:i),1), i=0,n_cpu-1)], &
                       MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
      ! q
      allocate(q(n_here), q_all(n_total))
      do j=1,n_here
        q(j) = p(j)%q
      end do
      call MPI_Gatherv(q(:), n_here, MPI_INTEGER, &
        q_all(:), particles_per_proc, [(sum(particles_per_proc(1:i),1), i=0,n_cpu-1)], &
        MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

      if (my_id .eq. 0) then
        call HDF5_array2D_saving(file,v_all,3,n_total,group_name//"v")
        call HDF5_array1D_saving_int(file,q_all,n_total,group_name//"q")
      end if
      deallocate(v,q,v_all,q_all)

    type is (particle_kinetic_leapfrog)
      particle_type_name = 'particle_kinetic_leapfrog'

      ! v
      allocate(v(3,n_here), v_all(3,n_total))
      do j=1,n_here
        v(:,j) = p(j)%v
      end do
!      call MPI_Gatherv(v(:,:), 3*n_here, MPI_REAL8, &
!        v_all(:,:), particles_per_proc*3, [(sum(particles_per_proc(1:i),1)*3, i=0,n_cpu-1)], &
!        MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

      call MPI_Gatherv(v(1,:), n_here, MPI_REAL8, &
                       v_all(1,:), particles_per_proc, [(sum(particles_per_proc(1:i),1), i=0,n_cpu-1)], &
                       MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
      call MPI_Gatherv(v(2,:), n_here, MPI_REAL8, &
                       v_all(2,:), particles_per_proc, [(sum(particles_per_proc(1:i),1), i=0,n_cpu-1)], &
                       MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
      call MPI_Gatherv(v(3,:), n_here, MPI_REAL8, &
                       v_all(3,:), particles_per_proc, [(sum(particles_per_proc(1:i),1), i=0,n_cpu-1)], &
                       MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

      ! q
      allocate(q(n_here), q_all(n_total))
      do j=1,n_here
        q(j) = p(j)%q
      end do
      call MPI_Gatherv(q(:), n_here, MPI_INTEGER, &
        q_all(:), particles_per_proc, [(sum(particles_per_proc(1:i),1), i=0,n_cpu-1)], &
        MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

      if (my_id .eq. 0) then
        call HDF5_array2D_saving(file,v_all,3,n_total,group_name//"v")
        call HDF5_array1D_saving_int(file,q_all,n_total,group_name//"q")
      end if
      deallocate(v,q,v_all,q_all)

    type is (particle_gc)
      particle_type_name = 'particle_gc'

      ! E
      allocate(E(n_here), E_all(n_total))
      do j=1,n_here
        E(j) = p(j)%E
      end do
      call MPI_Gatherv(E(:), n_here, MPI_REAL8, &
        E_all(:), particles_per_proc, [(sum(particles_per_proc(1:i),1), i=0,n_cpu-1)], &
        MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

      ! mu
      allocate(mu(n_here), mu_all(n_total))
      do j=1,n_here
        mu(j) = p(j)%mu
      end do
      call MPI_Gatherv(mu(:), n_here, MPI_REAL8, &
        mu_all(:), particles_per_proc, [(sum(particles_per_proc(1:i),1), i=0,n_cpu-1)], &
        MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

      ! q
      allocate(q(n_here), q_all(n_total))
      do j=1,n_here
        q(j) = p(j)%q
      end do
      call MPI_Gatherv(q(:), n_here, MPI_INTEGER, &
        q_all(:), particles_per_proc, [(sum(particles_per_proc(1:i),1), i=0,n_cpu-1)], &
        MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

      if (my_id .eq. 0) then
        call HDF5_array1D_saving(file,E_all,n_total,group_name//"E")
        call HDF5_array1D_saving(file,mu_all,n_total,group_name//"mu")
        call HDF5_array1D_saving_int(file,q_all,n_total,group_name//"q")
      end if
      deallocate(E,mu,q,E_all,mu_all,q_all)

    type is (particle_gc_vpar)
      particle_type_name = 'particle_gc_vpar'

      ! Vpar
      allocate(Vpar(n_here), Vpar_all(n_total))
      do j=1,n_here
        Vpar(j) = p(j)%Vpar
      end do
      call MPI_Gatherv(Vpar(:), n_here, MPI_REAL8, &
        Vpar_all(:), particles_per_proc, [(sum(particles_per_proc(1:i),1), i=0,n_cpu-1)], &
        MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

      ! mu
      allocate(mu(n_here), mu_all(n_total))
      do j=1,n_here
        mu(j) = p(j)%mu
      end do
      call MPI_Gatherv(mu(:), n_here, MPI_REAL8, &
        mu_all(:), particles_per_proc, [(sum(particles_per_proc(1:i),1), i=0,n_cpu-1)], &
        MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

      ! q
      allocate(q(n_here), q_all(n_total))
      do j=1,n_here
        q(j) = p(j)%q
      end do
      call MPI_Gatherv(q(:), n_here, MPI_INTEGER, &
        q_all(:), particles_per_proc, [(sum(particles_per_proc(1:i),1), i=0,n_cpu-1)], &
        MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

      if (my_id .eq. 0) then
        call HDF5_array1D_saving(file,Vpar_all,n_total,group_name//"Vpar")
        call HDF5_array1D_saving(file,mu_all,n_total,group_name//"mu")
        call HDF5_array1D_saving_int(file,q_all,n_total,group_name//"q")
      end if
      deallocate(Vpar,mu,q,Vpar_all,mu_all,q_all)

    type is (particle_fieldline)
      particle_type_name = 'particle_fieldline'
      ! v
      allocate(v1(n_here), v1_all(n_total))
      do j=1,n_here
        v1(j) = p(j)%v
      end do
      call MPI_Gatherv(v1(:), n_here, MPI_REAL8, &
        v1_all(:), particles_per_proc, [(sum(particles_per_proc(1:i),1), i=0,n_cpu-1)], &
        MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
      if (my_id .eq. 0) then
        call HDF5_array1D_saving(file,v1_all,n_total,group_name//"v")
      end if
      deallocate(v1, v1_all)

    type is (particle_kinetic_relativistic)
      particle_type_name = 'particle_kinetic_relativistic'
 
      ! momenta
      allocate(v(3,n_here), v_all(3,n_total))
      do j=1,n_here
        v(:,j) = p(j)%p
      end do
      call MPI_Gatherv(v(:,:), 3*n_here, MPI_REAL8, &
        v_all(:,:), particles_per_proc*3, [(sum(particles_per_proc(1:i),1)*3, i=0,n_cpu-1)], &
        MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

      ! q
      allocate(q(n_here), q_all(n_total))
      do j=1,n_here
        q(j) = p(j)%q
      end do
      call MPI_Gatherv(q(:), n_here, MPI_INTEGER, &
        q_all(:), particles_per_proc, [(sum(particles_per_proc(1:i),1), i=0,n_cpu-1)], &
        MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

      if (my_id .eq. 0) then
        call HDF5_array2D_saving(file,v_all,3,n_total,group_name//"v") 
        call HDF5_array1D_saving_int(file,q_all,n_total,group_name//"q") 
      end if
      deallocate(v,q,v_all,q_all)

    type is (particle_gc_relativistic)
      particle_type_name = 'particle_gc_relativistic'

      ! momenta (parallel momentum and magnetic moment)
      allocate(v(2,n_here), v_all(2,n_total))
      do j=1,n_here
        v(:,j) = p(j)%p
      end do
      call MPI_Gatherv(v(:,:), 2*n_here, MPI_REAL8, &
        v_all(:,:), particles_per_proc*2, [(sum(particles_per_proc(1:i),1)*2, i=0,n_cpu-1)], &
        MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

      ! q
      allocate(q(n_here), q_all(n_total))
      do j=1,n_here
        q(j) = p(j)%q
      end do
      call MPI_Gatherv(q(:), n_here, MPI_INTEGER, &
        q_all(:), particles_per_proc, [(sum(particles_per_proc(1:i),1), i=0,n_cpu-1)], &
        MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

      if (my_id .eq. 0) then
        call HDF5_array2D_saving(file,v_all,2,n_total,group_name//"v") 
        call HDF5_array1D_saving_int(file,q_all,n_total,group_name//"q") 
      end if
      deallocate(v,q,v_all,q_all)

    class default
      write(*,*) "error: missing type name declaration for write"
      call exit(1)
    end select


    ! It is important to do the gathering first, because that is the collective part
    if (my_id .eq. 0) then
      call HDF5_array2D_saving(file,x_all,3,n_total,group_name//"x")
      call HDF5_array2D_saving(file,st_all,2,n_total,group_name//"st")
      call HDF5_array1D_saving_r4(file,weight_all,n_total,group_name//"weight")
      call HDF5_array1D_saving_int(file,i_elm_all,n_total,group_name//"i_elm")
      call HDF5_array1D_saving_int(file,i_life_all,n_total,group_name//"i_life")
      call HDF5_array1D_saving_r4(file,t_birth_all,n_total,group_name//"t_birth")

      call HDF5_char_saving(file,particle_type_name,group_name//"type")
      call HDF5_integer_saving(file,sim%groups(i)%Z,group_name//"Z")
      call HDF5_real_saving(file,sim%groups(i)%mass,group_name//"mass")

      call HDF5_char_saving(file,sim%groups(i)%ad%suffix,group_name//"adas_suffix")
    end if
    deallocate(x,x_all,st,st_all,weight,weight_all,t_birth,t_birth_all,i_elm,i_elm_all,i_life,i_life_all)
  end do
end if

! Close everything
if (my_id .eq. 0) then
  call h5fclose_f(file, hdferr)
  call h5close_f(hdferr)

  write(*,*) "Writing particle output file to ", filename, " succeeded"
end if
end subroutine write_simulation_hdf5




!> Import all particles.
!> Reads the number of particles from a file, determines a
!> particle distribution over all processors and read this many
!> particles per processor.
subroutine read_simulation_hdf5(sim, filename)
use mod_openadas, only: read_adf11
use mod_coronal
type(particle_sim) , intent(inout) :: sim
character*(*)      , intent(in)  :: filename

integer                            :: my_id, n_cpu, ierr, rank, n_particles
integer, allocatable, dimension(:) :: particles_per_proc

! For HDF5 reading
integer(HID_T)    :: file, file_space, mem_space, dset, plist ! handles
integer(HID_T)    :: data_type
integer(HID_T)    :: group_id
integer(HID_T)    :: time_set_id
integer(HSIZE_T)  :: i_here
integer           :: n_here
integer           :: storage_type, max_corder
character(len=12) :: group_name
character(len=particle_type_name_length) :: particle_type_name
integer           :: i, j, n, hdferr, n_alive
integer, allocatable :: n_alive_all(:)
logical           :: exists

type(c_ptr) :: p_ptr
integer*8, dimension(1:2) :: tmp, maxdims
real*8, dimension(:,:), allocatable :: real8_2D
integer*4, dimension(:), allocatable :: int4_1D
real*4, dimension(:), allocatable :: real4_1D
real*8, dimension(:), allocatable :: real8_1D

! Preparation
call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)      ! id of each MPI proc
call MPI_COMM_SIZE(MPI_COMM_WORLD, n_cpu, ierr)      ! number of MPI procs
call h5open_f(hdferr)
allocate(particles_per_proc(0:n_cpu-1))

! Create file property list for parallel access
call h5pcreate_f(H5P_FILE_ACCESS_F, plist, hdferr)
call h5pset_fapl_mpio_f(plist, MPI_COMM_WORLD, MPI_INFO_NULL, hdferr)

! Open the file
call h5fopen_f(filename, H5F_ACC_RDONLY_F, file, hdferr, access_prp=plist)
call h5pclose_f(plist, hdferr)

! read the time
call HDF5_real_reading(file,sim%time,'/time')

! Get the number of groups
call h5gopen_f(file, '/groups/', group_id, hdferr)
call h5gget_info_f(group_id, storage_type, n, max_corder, hdferr)
call h5gclose_f(group_id, hdferr)

! Reallocate groups if necessary
if (allocated(sim%groups)) deallocate(sim%groups)
allocate(sim%groups(n))
do i=1,n
  ! Open the dataset for x
  write(group_name,'(A,i0.3,A)') '/groups/', i, '/'
  call h5dopen_f(file, group_name//"x", dset, hdferr)

  ! Open the file dataspace
  call h5dget_space_f(dset, file_space, hdferr)

  ! Get the number of particles
  call h5sget_simple_extent_ndims_f(file_space, rank, hdferr)
  call h5sget_simple_extent_dims_f(file_space, tmp, maxdims, hdferr)
  n_particles = int(tmp(2),4)
  call h5sclose_f(file_space, hdferr)
  call h5dclose_f(dset, hdferr)

  ! Divide particles over processors
  particles_per_proc(0)         = n_particles/n_cpu + modulo(n_particles, n_cpu)
  particles_per_proc(1:n_cpu-1) = n_particles/n_cpu
  i_here = sum(particles_per_proc(0:my_id-1))
  n_here = particles_per_proc(my_id)

  ! Get the particle type from the attribute
  call HDF5_char_reading(file,particle_type_name,group_name//"type")

  ierr = 0
  select case (trim(particle_type_name))
  case ('particle_kinetic')
    allocate(particle_kinetic::sim%groups(i)%particles(n_here), stat=ierr)
  case ('particle_kinetic_leapfrog')
    allocate(particle_kinetic_leapfrog::sim%groups(i)%particles(n_here), stat=ierr)
  case ('particle_gc')
    allocate(particle_gc::sim%groups(i)%particles(n_here), stat=ierr)
  case ('particle_gc_vpar')
    allocate(particle_gc_vpar::sim%groups(i)%particles(n_here), stat=ierr)
  case ('particle_fieldline')
    allocate(particle_fieldline::sim%groups(i)%particles(n_here), stat=ierr)
  case ('particle_kinetic_relativistic')
    allocate(particle_kinetic_relativistic::sim%groups(i)%particles(n_here), stat=ierr)
  case ('particle_gc_relativistic')
    allocate(particle_gc_relativistic::sim%groups(i)%particles(n_here), stat=ierr)
  case default
    write(*,*) "error: missing type name declaration ", trim(particle_type_name), " for read"
    call exit(1)
  end select
  if (ierr .gt. 0) write(*,"(i3,a,i12,a)") my_id, &
      "unable to allocate particles(", particles_per_proc(my_id), ")"

  call HDF5_integer_reading(file,sim%groups(i)%Z,group_name//"Z")
  call HDF5_real_reading(file,sim%groups(i)%mass,group_name//"mass")
  call HDF5_char_reading(file,sim%groups(i)%ad%suffix,group_name//"adas_suffix")
  if (len_trim(sim%groups(i)%ad%suffix) .gt. 0) then
    sim%groups(i)%ad = read_adf11(my_id,sim%groups(i)%ad%suffix)
    sim%groups(i)%cor = coronal(sim%groups(i)%ad)
  end if

  ! Read base particle attributes
  ! x
  allocate(real8_2D(3,n_here))
  call HDF5_array2D_reading(file, real8_2D, group_name//"x",start=[0_HSIZE_T,i_here])
  do j=1,n_here
    sim%groups(i)%particles(j)%x = real8_2D(:,j)
  end do
  deallocate(real8_2D)

  ! st
  allocate(real8_2D(2,n_here))
  call HDF5_array2D_reading(file, real8_2D, group_name//"st",start=[0_HSIZE_T,i_here])
  do j=1,n_here
    sim%groups(i)%particles(j)%st = real8_2D(:,j)
  end do
  deallocate(real8_2D)

  ! weight
  allocate(real4_1D(n_here))
  call HDF5_array1D_reading_r4(file, real4_1D, group_name//"weight",start=[i_here])
  do j=1,n_here
    sim%groups(i)%particles(j)%weight = real4_1D(j)
  end do
  deallocate(real4_1D)
  ! i_elm
  allocate(int4_1D(n_here))
  call HDF5_array1D_reading_int(file, int4_1D, group_name//"i_elm", start=[i_here])
  do j=1,n_here
    sim%groups(i)%particles(j)%i_elm = int4_1D(j)
  end do
  deallocate(int4_1D)

  ! The following two are relatively new, so might not always be present
  ! i_life
  call h5lexists_f(file, group_name//"i_life", exists, ierr)
  if (exists) then
    allocate(int4_1D(n_here))
    int4_1D = 0 ! preset to 0 in case not present (i.e. old restart files)
    ! will still give a nasty error message, but should work,
    ! though who knows what hdf5 does with our array if reading fails...
    call HDF5_array1D_reading_int(file, int4_1D, group_name//"i_life", start=[i_here])
    do j=1,n_here
      sim%groups(i)%particles(j)%i_life = int4_1D(j)
    end do
    deallocate(int4_1D)
  end if

  ! t_birth
  call h5lexists_f(file, group_name//"t_birth", exists, ierr)
  if (exists) then
    allocate(real4_1D(n_here))
    real4_1D = 0.0 ! preset to 0 in case not present (i.e. old restart files)
    call HDF5_array1D_reading_r4(file, real4_1D, group_name//"t_birth",start=[i_here])
    do j=1,n_here
      sim%groups(i)%particles(j)%t_birth = real4_1D(j)
    end do
    deallocate(real4_1D)
  end if

  select type (p => sim%groups(i)%particles)
  type is (particle_kinetic)
    ! v
    allocate(real8_2D(3,n_here))
    call HDF5_array2D_reading(file, real8_2D, group_name//"v",start=[0_HSIZE_T,i_here])
    do j=1,n_here
      p(j)%v = real8_2D(:,j)
    end do
    deallocate(real8_2D)

    ! q
    allocate(int4_1D(n_here))
    call HDF5_array1D_reading_int(file, int4_1D, group_name//"q", start=[i_here])
    do j=1,n_here
      p(j)%q = int4_1D(j)
    end do
    deallocate(int4_1D)

  type is (particle_kinetic_leapfrog)
    ! v
    allocate(real8_2D(3,n_here))
    call HDF5_array2D_reading(file, real8_2D, group_name//"v",start=[0_HSIZE_T,i_here])
    do j=1,n_here
      p(j)%v = real8_2D(:,j)
    end do
    deallocate(real8_2D)

    ! q
    allocate(int4_1D(n_here))
    call HDF5_array1D_reading_int(file, int4_1D, group_name//"q", start=[i_here])
    do j=1,n_here
      p(j)%q = int4_1D(j)
    end do
    deallocate(int4_1D)

  type is (particle_gc)
    ! E
    allocate(real8_1D(n_here))
    call HDF5_array1D_reading(file, real8_1D, group_name//"E",start=[i_here])
    do j=1,n_here
      p(j)%E = real8_1D(j)
    end do
    deallocate(real8_1D)
    ! mu
    allocate(real8_1D(n_here))
    call HDF5_array1D_reading(file, real8_1D, group_name//"mu",start=[i_here])
    do j=1,n_here
      p(j)%mu = real8_1D(j)
    end do
    deallocate(real8_1D)
    ! q
    allocate(int4_1D(n_here))
    call HDF5_array1D_reading_int(file, int4_1D, group_name//"q", start=[i_here])
    do j=1,n_here
      p(j)%q = int4_1D(j)
    end do
    deallocate(int4_1D)
  
  type is (particle_gc_vpar)
    ! 
    allocate(real8_1D(n_here))
    call HDF5_array1D_reading(file, real8_1D, group_name//"Vpar",start=[i_here])
    do j=1,n_here
      p(j)%Vpar = real8_1D(j)
    end do
    deallocate(real8_1D)
    ! mu
    allocate(real8_1D(n_here))
    call HDF5_array1D_reading(file, real8_1D, group_name//"mu",start=[i_here])
    do j=1,n_here
      p(j)%mu = real8_1D(j)
    end do
    deallocate(real8_1D)
    ! q
    allocate(int4_1D(n_here))
    call HDF5_array1D_reading_int(file, int4_1D, group_name//"q", start=[i_here])
    do j=1,n_here
      p(j)%q = int4_1D(j)
    end do
    deallocate(int4_1D)
  
  type is (particle_fieldline)
    ! v
    allocate(real8_1D(n_here))
    call HDF5_array1D_reading(file, real8_1D, group_name//"v",start=[i_here])
    do j=1,n_here
      p(j)%v = real8_1D(j)
    end do
    deallocate(real8_1D)

  type is (particle_kinetic_relativistic)

    ! momenta [AMU*m/s]
    allocate(real8_2D(3,n_here))
    call HDF5_array2D_reading(file, real8_2D, group_name//"v",start=[0_HSIZE_T,i_here])
    do j=1,n_here
      p(j)%p = real8_2D(:,j)
    end do
    deallocate(real8_2D)

    ! electric charge q
    allocate(int4_1D(n_here))
    call HDF5_array1D_reading_int(file, int4_1D, group_name//"q", start=[i_here])
    do j=1,n_here
      p(j)%q = int4_1D(j)
    end do
    deallocate(int4_1D)    

  type is (particle_gc_relativistic)

    ! momenta (parallel momentum and magnetic moment)
    allocate(real8_2D(2,n_here))
    call HDF5_array2D_reading(file, real8_2D, group_name//"v",start=[0_HSIZE_T,i_here])
    do j=1,n_here
      p(j)%p = real8_2D(:,j)
    end do
    deallocate(real8_2D)

    ! electric charge q
    allocate(int4_1D(n_here))
    call HDF5_array1D_reading_int(file, int4_1D, group_name//"q", start=[i_here])
    do j=1,n_here
      p(j)%q = int4_1D(j)
    end do
    deallocate(int4_1D)    

  end select

  ! Check if the balance between processors is okay, by comparing the number of
  ! alive particles
  n_alive = count(sim%groups(i)%particles(:)%i_elm .gt. 0)
  allocate(n_alive_all(sim%n_cpu))
  call MPI_Gather(n_alive, 1, MPI_INTEGER, n_alive_all, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

  ! Check if the imbalance is not too great
  if (sim%my_id .eq. 0) then
    if (maxval(n_alive_all) .gt. minval(n_alive_all) * 1.5) then
      write(*,*) "WARNING: ", (maxval(n_alive_all)*100)/minval(n_alive_all), '% imbalance between CPUs, counts:'
      write(*,*) n_alive_all
    end if
  end if
  deallocate(n_alive_all)

end do

! Close everything else
call h5fclose_f(file, hdferr)
call h5close_f(hdferr)

end subroutine read_simulation_hdf5


!> Get '/time' from a file. Does not alter the units in any way
!> code works for jorek and particle restart files, and returns values in
!> different units for both
function get_simulation_hdf5_time(filename) result(time)
  character*(*)      , intent(in)  :: filename
  real*8 :: time
  integer(HID_T) :: file
  integer :: hdferr
  call h5open_f(hdferr)
  call h5fopen_f(filename, H5F_ACC_RDONLY_F, file, hdferr)
  call HDF5_real_reading(file,time,'/time')
  call h5fclose_f(file,hdferr)
  call h5close_f(hdferr)
end function get_simulation_hdf5_time
end module mod_particle_io
