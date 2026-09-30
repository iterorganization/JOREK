!> Program for hot tail RE generation during a disruption
!> This code follows Guiding Centers (GC) in an already computed MHD simulation from JOREK
!> More details about the implementation and verification can be found in: L. Puel et al 2026 NF 10.1088/1741-4326/ae2d70
!> The present parameters are used for the study of the baseline H-mode 15 MA ITER case presented in: L. Puel et al 2026 NF 10.1088/1741-4326/aea175
!> It represents a good start to study other cases.
!> See the related documentation in the wiki for details

module testccoll_helpers
  use particle_tracer
  use mod_particle_io
  use mod_particle_diagnostics
  use mod_gc_relativistic
  use hdf5_io_module
  use constants, only: EL_CHG, SPEED_OF_LIGHT, ATOMIC_MASS_UNIT, PI
  use mod_coordinate_transforms, only: vector_cylindrical_to_cartesian
  implicit none

contains

!> Returns the total momentum normalised to m*c and the cosine of the pitch angle
!> of a relativistic guiding centre, from its parallel momentum and magnetic moment.
pure subroutine gc_momentum_pitch(p_par, mu, Bnorm, mass, p_over_mc, xi)
  implicit none
  real*8, intent(in)  :: p_par     !< parallel momentum [AMU m/s]
  real*8, intent(in)  :: mu        !< magnetic moment [(AMU m^2)/(T s^2)]
  real*8, intent(in)  :: Bnorm     !< magnetic field magnitude [T]
  real*8, intent(in)  :: mass      !< marker mass [AMU]
  real*8, intent(out) :: p_over_mc !< total momentum normalised to m*c [1]
  real*8, intent(out) :: xi        !< cosine of the pitch angle [1]

  real*8 :: pnorm

  pnorm     = sqrt(max(2.d0*mu*Bnorm*mass + p_par*p_par, 0.d0))
  p_over_mc = pnorm / (mass * SPEED_OF_LIGHT)
  if (pnorm .gt. 0.d0) then
    xi = p_par / pnorm
  else
    xi = 0.d0
  end if
end subroutine gc_momentum_pitch

!> Draws a marker index uniformly in [1, nprt].
integer function draw_marker_index(nprt) result(idx)
  implicit none
  integer, intent(in) :: nprt !< number of markers to draw from

  real*8 :: rnd

  call random_number(rnd)
  idx = min(int(rnd * real(nprt, 8)) + 1, nprt)
end function draw_marker_index

!> Creates (or truncates) an HDF5 output file on rank 0.
subroutine create_hdf5(fnout)
  implicit none
  character(len=*), intent(in) :: fnout !< name of the file to create

  integer(HID_T) :: file
  integer :: hdferr, ierr, rank

  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  if (rank .eq. 0) then
    call h5open_f(hdferr)
    call h5fcreate_f(trim(fnout), H5F_ACC_TRUNC_F, file, hdferr)
    if (hdferr .ne. 0) then
      write(*,*) "ERROR: could not create HDF5 file ", trim(fnout), " (hdferr =", hdferr, ")"
      call h5close_f(hdferr)
      call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    end if
    call h5fclose_f(file, hdferr)
    call h5close_f(hdferr)
  end if
end subroutine create_hdf5

!> Appends one column (one output time) to a 2D dataset [n_marker, n_time] of an
!> already open HDF5 file. The dataset is stored in single precision, which is
!> ample for the n, T and E marker diagnostics; the conversion is done explicitly
!> so that the in-memory and in-file types always match.
subroutine append_2Dreal_column_to_hdf5(file_id, datasetname, array)
  use hdf5
  implicit none
  integer(HID_T),   intent(in) :: file_id     !< identifier of an open HDF5 file
  character(len=*), intent(in) :: datasetname !< name of the dataset
  real*8,           intent(in) :: array(:)    !< column to append

  integer(HID_T)   :: dset_id, dspace_id, memspace_id, plist_id
  integer(HSIZE_T) :: dims(2), maxdims(2), start(2), count(2), chunk_dims(2)
  integer          :: error
  logical          :: exists
  real*4, allocatable :: buffer(:)

  call h5lexists_f(file_id, datasetname, exists, error)
  if (error .ne. 0) then
    write(*,*) "ERROR: h5lexists_f failed for ", trim(datasetname)
    return
  end if

  if (.not. exists) then
    !> Create an extendable dataset: first dimension = markers, second = time
    dims       = [size(array, KIND=HSIZE_T), 0_HSIZE_T]
    maxdims    = [size(array, KIND=HSIZE_T), H5S_UNLIMITED_F]
    chunk_dims = [size(array, KIND=HSIZE_T), 1_HSIZE_T]

    call h5screate_simple_f(2, dims, dspace_id, error, maxdims)
    call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
    call h5pset_chunk_f(plist_id, 2, chunk_dims, error)
    call h5dcreate_f(file_id, datasetname, H5T_NATIVE_REAL, dspace_id, &
                     dset_id, error, plist_id)
    call h5pclose_f(plist_id, error)
    call h5sclose_f(dspace_id, error)
  else
    call h5dopen_f(file_id, datasetname, dset_id, error)
  end if
  if (error .ne. 0) then
    write(*,*) "ERROR: could not open or create dataset ", trim(datasetname)
    return
  end if

  !> Extend the dataset by one column
  call h5dget_space_f(dset_id, dspace_id, error)
  call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, error)
  call h5sclose_f(dspace_id, error)
  dims(2) = dims(2) + 1_HSIZE_T
  call h5dset_extent_f(dset_id, dims, error)

  !> Select the newly created column
  start = [0_HSIZE_T, dims(2) - 1_HSIZE_T]
  count = [size(array, KIND=HSIZE_T), 1_HSIZE_T]
  call h5dget_space_f(dset_id, dspace_id, error)
  call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, start, count, error)
  call h5screate_simple_f(1, [size(array, KIND=HSIZE_T)], memspace_id, error)

  allocate(buffer(size(array)))
  buffer = real(array, 4)
  call h5dwrite_f(dset_id, H5T_NATIVE_REAL, buffer, [size(array, KIND=HSIZE_T)], error, &
                  mem_space_id=memspace_id, file_space_id=dspace_id)
  deallocate(buffer)

  call h5sclose_f(memspace_id, error)
  call h5sclose_f(dspace_id, error)
  call h5dclose_f(dset_id, error)
end subroutine append_2Dreal_column_to_hdf5

!> Appends a 1D array at the end of an extendable 1D dataset of an already open file.
subroutine append_1Dreal_array_to_hdf5(file_id, datasetname, array)
  use hdf5
  implicit none
  integer(HID_T),   intent(in) :: file_id     !< identifier of an open HDF5 file
  character(len=*), intent(in) :: datasetname !< name of the dataset
  real*8,           intent(in) :: array(:)    !< data to append

  integer(HID_T)   :: dset_id, dspace_id, memspace_id, plist_id
  integer(HSIZE_T) :: dims(1), maxdims(1), offset(1), count(1), chunk_dims(1)
  integer          :: error
  logical          :: exists

  call h5lexists_f(file_id, datasetname, exists, error)
  if (error .ne. 0) then
    write(*,*) "ERROR: h5lexists_f failed for ", trim(datasetname)
    return
  end if

  if (.not. exists) then
    !> Create a dataset with an unlimited first dimension
    dims(1)       = 0
    maxdims(1)    = H5S_UNLIMITED_F
    chunk_dims(1) = max(1, min(1024, size(array)))

    call h5screate_simple_f(1, dims, dspace_id, error, maxdims)
    call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
    call h5pset_chunk_f(plist_id, 1, chunk_dims, error)
    call h5dcreate_f(file_id, datasetname, H5T_NATIVE_DOUBLE, &
                     dspace_id, dset_id, error, plist_id)
    call h5pclose_f(plist_id, error)
    call h5sclose_f(dspace_id, error)
  else
    call h5dopen_f(file_id, datasetname, dset_id, error)
  end if
  if (error .ne. 0) then
    write(*,*) "ERROR: could not open or create dataset ", trim(datasetname)
    return
  end if

  !> Nothing to append: the dataset exists and is left unchanged
  if (size(array) .eq. 0) then
    call h5dclose_f(dset_id, error)
    return
  end if

  !> Extend the dataset by size(array) entries
  call h5dget_space_f(dset_id, dspace_id, error)
  call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, error)
  call h5sclose_f(dspace_id, error)

  offset(1) = dims(1)
  count(1)  = size(array, kind=HSIZE_T)
  dims(1)   = offset(1) + count(1)
  call h5dset_extent_f(dset_id, dims, error)

  !> Select the new entries and write
  call h5dget_space_f(dset_id, dspace_id, error)
  call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, count, error)
  call h5screate_simple_f(1, count, memspace_id, error)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, array, count, error, &
                  memspace_id, dspace_id)

  call h5sclose_f(memspace_id, error)
  call h5sclose_f(dspace_id, error)
  call h5dclose_f(dset_id, error)
end subroutine append_1Dreal_array_to_hdf5

!> Appends a scalar at the end of an extendable 1D dataset of an already open file.
subroutine append_scalar_to_hdf5(file_id, datasetname, scalar)
  use hdf5
  implicit none
  integer(HID_T),   intent(in) :: file_id     !< identifier of an open HDF5 file
  character(len=*), intent(in) :: datasetname !< name of the dataset
  real*8,           intent(in) :: scalar      !< value to append

  real*8 :: array(1)

  array(1) = scalar
  call append_1Dreal_array_to_hdf5(file_id, datasetname, array)
end subroutine append_scalar_to_hdf5

!> Appends values at the end of an allocatable array, allocating it if needed.
subroutine append_array(array, values)
  implicit none
  real*8, allocatable, intent(inout) :: array(:)  !< array to grow
  real*8,              intent(in)    :: values(:) !< values to append

  real*8, allocatable :: tmp(:)
  integer :: n, m

  if (.not. allocated(array)) then
    allocate(array(size(values)))
    array = values
  else
    n = size(array)
    m = size(values)
    allocate(tmp(n + m))
    tmp(1:n)     = array
    tmp(n+1:n+m) = values
    call move_alloc(tmp, array)
  end if
end subroutine append_array

!> Appends the entries of source selected by mask at the end of array, and resets
!> the selected entries of source to zero.
subroutine append_filtered(array, source, mask)
  implicit none
  real*8, allocatable, intent(inout) :: array(:)  !< array to grow
  real*8,              intent(inout) :: source(:) !< source values, selected entries are reset
  logical,             intent(in)    :: mask(:)   !< selection mask

  real*8, allocatable :: to_append(:)
  integer :: count_true, i, j

  count_true = count(mask)
  if (count_true .eq. 0) return

  allocate(to_append(count_true))
  j = 0
  do i = 1, size(source)
    if (mask(i)) then
      j = j + 1
      to_append(j) = source(i)
      source(i)    = 0.d0
    end if
  end do

  call append_array(array, to_append)
end subroutine append_filtered

!> Splits the markers evenly over the MPI ranks and allocates the per-rank
!> collision-field diagnostic buffers.
subroutine MPI_tracker_init(nprt, rank, n_mpi, ierr, start_iprt, end_iprt, nprt_MPI, &
                            n_write_MPI, T_write_MPI, E_write_MPI)
  implicit none
  !> Input
  integer, intent(in) :: nprt !< total number of markers
  !> In/Out
  integer, intent(inout) :: rank, n_mpi, ierr
  integer, intent(inout) :: start_iprt, end_iprt, nprt_MPI
  real*8, allocatable, intent(inout) :: n_write_MPI(:), T_write_MPI(:), E_write_MPI(:)

  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, n_mpi, ierr)

  if (mod(nprt, n_mpi) .ne. 0) then
    if (rank .eq. 0) then
      write(*,*) "ERROR: the marker number (", nprt, ") is not divisible by the &
                 &number of MPI processes (", n_mpi, ")"
    end if
    !> MPI_Abort instead of STOP: a STOP on a single rank leaves the other ranks
    !> waiting in the next collective call
    call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  end if

  start_iprt = 1 + rank * (nprt / n_mpi)
  end_iprt   = min((rank + 1) * (nprt / n_mpi), nprt)
  nprt_MPI   = end_iprt - start_iprt + 1

  allocate(T_write_MPI(nprt_MPI), E_write_MPI(nprt_MPI), n_write_MPI(nprt_MPI))
  T_write_MPI = 0.d0
  E_write_MPI = 0.d0
  n_write_MPI = 0.d0

  write(*,*) "> MPI process / start_part, end_part :", rank, start_iprt, end_iprt
end subroutine MPI_tracker_init

!> Gathers the per-marker collision fields on rank 0 and appends one time column
!> to the diagnostic file. The file and the HDF5 library are opened once per call.
subroutine write_collparam_hdf5(fnout_collparam, n_write, T_write, E_write, &
                                n_write_MPI, T_write_MPI, E_write_MPI, nprt_MPI, time)
  implicit none
  character(len=*), intent(in)    :: fnout_collparam !< output file name
  real*8,           intent(inout) :: n_write(:), T_write(:), E_write(:)         !< global buffers (rank 0)
  real*8,           intent(inout) :: n_write_MPI(:), T_write_MPI(:), E_write_MPI(:) !< per-rank buffers, reset on exit
  integer,          intent(in)    :: nprt_MPI !< number of markers on this rank
  real*8,           intent(in)    :: time     !< current time [s]

  integer, allocatable :: particles_per_proc(:), displs(:)
  integer(HID_T) :: file_id
  integer :: ierr, rank, n_cpus, i, error

  call MPI_Comm_size(MPI_COMM_WORLD, n_cpus, ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

  !> recvcounts and displs are only meaningful on the root, but they must be
  !> defined everywhere: they are read by the MPI implementation on every rank
  allocate(particles_per_proc(0:n_cpus-1), displs(0:n_cpus-1))
  particles_per_proc = 0
  displs             = 0

  call MPI_Gather(nprt_MPI, 1, MPI_INTEGER, &
                  particles_per_proc, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  if (rank .eq. 0) then
    do i = 1, n_cpus - 1
      displs(i) = displs(i-1) + particles_per_proc(i-1)
    end do
  end if

  call MPI_Gatherv(n_write_MPI, size(n_write_MPI, 1), MPI_REAL8, &
                   n_write, particles_per_proc, displs, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
  call MPI_Gatherv(T_write_MPI, size(T_write_MPI, 1), MPI_REAL8, &
                   T_write, particles_per_proc, displs, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
  call MPI_Gatherv(E_write_MPI, size(E_write_MPI, 1), MPI_REAL8, &
                   E_write, particles_per_proc, displs, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

  if (rank .eq. 0) then
    call h5open_f(error)
    call h5fopen_f(trim(fnout_collparam), H5F_ACC_RDWR_F, file_id, error)
    if (error .ne. 0) then
      write(*,*) "ERROR: unable to open HDF5 file ", trim(fnout_collparam)
    else
      call append_2Dreal_column_to_hdf5(file_id, "n", n_write)
      call append_2Dreal_column_to_hdf5(file_id, "T", T_write)
      call append_2Dreal_column_to_hdf5(file_id, "E", E_write)
      call append_scalar_to_hdf5(file_id, "time", time)
      call h5fclose_f(file_id, error)
    end if
    call h5close_f(error)
  end if

  deallocate(particles_per_proc, displs)

  n_write_MPI = 0.d0
  T_write_MPI = 0.d0
  E_write_MPI = 0.d0
end subroutine write_collparam_hdf5

!> Gathers the markers lost since the previous call on rank 0, appends them to the
!> loss file and empties the per-rank loss buffers.
subroutine write_loss_to_hdf5(fnout_losses, t_losses, p_losses, xi_losses, &
                              R_losses, Z_losses, Phi_losses, Weight_losses, &
                              iprt_losses)
  use hdf5
  implicit none
  character(len=*), intent(in) :: fnout_losses !< output file name

  real*8, allocatable, intent(inout) :: t_losses(:), p_losses(:), xi_losses(:)
  real*8, allocatable, intent(inout) :: R_losses(:), Z_losses(:)
  real*8, allocatable, intent(inout) :: Phi_losses(:), Weight_losses(:)
  real*8, allocatable, intent(inout) :: iprt_losses(:)

  ! MPI
  integer :: ierr, rank, nprocs
  integer :: sendcount, total_size, root_size, i
  integer, allocatable :: recvcounts(:), displs(:)

  ! Global arrays (meaningful on rank 0 only)
  real*8, allocatable :: t_all(:), p_all(:), xi_all(:)
  real*8, allocatable :: R_all(:), Z_all(:)
  real*8, allocatable :: Phi_all(:), Weight_all(:)
  real*8, allocatable :: iprt_all(:)

  ! HDF5
  integer(HID_T) :: file_id
  integer :: error

  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)

  sendcount = size(t_losses)

  allocate(recvcounts(nprocs), displs(nprocs))

  !> The counts are one integer per rank, so an Allgather is cheap and keeps the
  !> displacements consistent on every rank
  call MPI_Allgather(sendcount, 1, MPI_INTEGER, &
                     recvcounts, 1, MPI_INTEGER, &
                     MPI_COMM_WORLD, ierr)

  displs(1) = 0
  do i = 2, nprocs
    displs(i) = displs(i-1) + recvcounts(i-1)
  end do
  total_size = sum(recvcounts)

  !> The gathered size is reported because it is the quantity that drives the
  !> root allocation below, and therefore the memory footprint of this routine.
  if (rank .eq. 0) then
    write(*,*) 'INFO: writing ', total_size, ' loss records (largest rank buffer:', &
               maxval(recvcounts), ')'
    flush(6)
  end if

  !> Only rank 0 writes, so only rank 0 needs the full arrays. The other ranks
  !> allocate a single element so that the receive buffer address stays valid.
  if (rank .eq. 0) then
    root_size = max(total_size, 1)
  else
    root_size = 1
  end if
  allocate(t_all(root_size), p_all(root_size), xi_all(root_size))
  allocate(R_all(root_size), Z_all(root_size), Phi_all(root_size))
  allocate(Weight_all(root_size), iprt_all(root_size))

  call MPI_Gatherv(t_losses, sendcount, MPI_REAL8, &
                   t_all, recvcounts, displs, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
  call MPI_Gatherv(p_losses, sendcount, MPI_REAL8, &
                   p_all, recvcounts, displs, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
  call MPI_Gatherv(xi_losses, sendcount, MPI_REAL8, &
                   xi_all, recvcounts, displs, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
  call MPI_Gatherv(R_losses, sendcount, MPI_REAL8, &
                   R_all, recvcounts, displs, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
  call MPI_Gatherv(Z_losses, sendcount, MPI_REAL8, &
                   Z_all, recvcounts, displs, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
  call MPI_Gatherv(Phi_losses, sendcount, MPI_REAL8, &
                   Phi_all, recvcounts, displs, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
  call MPI_Gatherv(Weight_losses, sendcount, MPI_REAL8, &
                   Weight_all, recvcounts, displs, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
  call MPI_Gatherv(iprt_losses, sendcount, MPI_REAL8, &
                   iprt_all, recvcounts, displs, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

  !> The file is written even when no loss occurred, so that the datasets always
  !> exist and the post-processing does not have to special-case an empty run
  if (rank .eq. 0) then
    call h5open_f(error)
    if (error .ne. 0) write(*,*) "ERROR: h5open_f failed"

    !> Open the file once and pass the identifier to every append
    call h5fopen_f(trim(fnout_losses), H5F_ACC_RDWR_F, file_id, error)
    if (error .ne. 0) then
      write(*,*) "ERROR: unable to open HDF5 file:"
      write(*,*) trim(fnout_losses)
      write(*,*) "HDF5 error =", error
    else
      call append_1Dreal_array_to_hdf5(file_id, "t",      t_all(1:total_size))
      call append_1Dreal_array_to_hdf5(file_id, "p",      p_all(1:total_size))
      call append_1Dreal_array_to_hdf5(file_id, "xi",     xi_all(1:total_size))
      call append_1Dreal_array_to_hdf5(file_id, "R",      R_all(1:total_size))
      call append_1Dreal_array_to_hdf5(file_id, "Z",      Z_all(1:total_size))
      call append_1Dreal_array_to_hdf5(file_id, "Phi",    Phi_all(1:total_size))
      call append_1Dreal_array_to_hdf5(file_id, "Weight", Weight_all(1:total_size))
      call append_1Dreal_array_to_hdf5(file_id, "iprt",   iprt_all(1:total_size))
      call h5fclose_f(file_id, error)
    end if

    call h5close_f(error)
  end if

  deallocate(t_all, p_all, xi_all)
  deallocate(R_all, Z_all, Phi_all)
  deallocate(Weight_all, iprt_all)
  deallocate(recvcounts, displs)

  deallocate(t_losses, p_losses, xi_losses)
  deallocate(R_losses, Z_losses, Phi_losses)
  deallocate(Weight_losses, iprt_losses)

  allocate(t_losses(0), p_losses(0), xi_losses(0))
  allocate(R_losses(0), Z_losses(0), Phi_losses(0))
  allocate(Weight_losses(0), iprt_losses(0))
end subroutine write_loss_to_hdf5

!> Correlated Russian roulette. A marker that has just fallen below p_cutoff is a
!> candidate: it survives with probability l/k with its weight multiplied by k/l,
!> otherwise it is killed and replaced by a copy of a marker above p_far, the two
!> copies sharing the donor weight.
!>
!> The survival draws come from a barrel of k slots containing exactly l ones, so
!> that exactly l out of every k candidates survive. The barrel click therefore
!> advances on candidates only.
subroutine correlated_russian_roulette(particles, fields, nprt, time, mass, p_cutoff, p_far, &
                                       p_before_reweight, prt_age, age_split, k, l, &
                                       use_psi_res, use_p_far_sweep, psi_res, delta_psi_res, delta_p_far)
  implicit none
  class(particle_base), dimension(:), intent(inout) :: particles
  class(fields_base),                 intent(in)    :: fields
  real*8,                             intent(inout) :: p_before_reweight(nprt)
  integer*4,                          intent(inout) :: prt_age(nprt)
  integer,                            intent(in)    :: nprt, k, l, age_split
  real*8,                             intent(in)    :: p_cutoff, p_far, time, mass
  logical,                            intent(in)    :: use_psi_res, use_p_far_sweep
  real*8,                             intent(in)    :: psi_res, delta_psi_res, delta_p_far

  !> Number of draws, in units of nprt, after which the search is declared failed
  integer, parameter :: MAX_DRAWS_PER_MARKER = 20

  integer, allocatable :: barrel(:)
  integer*4 :: i_click, iprt, new_iprt, i_loop_counter, ierr
  real*8    :: new_weight, pnorm, xi, p_far_tmp
  real*8    :: E(3), B(3), psi_old, psi_new, psi_res_tmp, U
  logical   :: is_new_part

  allocate(barrel(k))
  call generate_barrel(barrel, k, l)
  i_click = 1

  select type (prt => particles)
  type is (particle_gc_relativistic)
    do iprt = 1, nprt

      !> Defensive guard: calling the field interpolation with i_elm < 1 indexes
      !> element_list out of bounds. All invalid markers are recycled at the end
      !> of each time step, so this branch should never be taken.
      if (prt(iprt)%i_elm .lt. 1 .or. &
          prt(iprt)%i_elm .gt. fields%element_list%n_elements) then
        p_before_reweight(iprt) = 0.d0
        prt_age(iprt)           = prt_age(iprt) + 1
        cycle
      end if

      psi_res_tmp = psi_res
      p_far_tmp   = p_far

      !> Only B and psi are needed here, so the E = eta*J_par normalisation is
      !> skipped.
      call fields%calc_EBpsiU(time, prt(iprt)%i_elm, prt(iprt)%st, prt(iprt)%x(3), &
                              E, B, psi_old, U, E_par_norm=.false.)
      call gc_momentum_pitch(prt(iprt)%p(1), prt(iprt)%p(2), norm2(B), mass, pnorm, xi)

      if ((pnorm .le. p_cutoff) .and. (p_before_reweight(iprt) .gt. p_cutoff)) then
        if (barrel(i_click) .eq. 1) then
          !> Survive: the weight is increased by the inverse survival probability
          prt(iprt)%weight = prt(iprt)%weight * real(k, 8) / real(l, 8)
        else
          !> Killed: look for a marker above p_far to duplicate
          is_new_part    = .false.
          i_loop_counter = 0
          do while (.not. is_new_part)
            i_loop_counter = i_loop_counter + 1
            if (i_loop_counter .gt. MAX_DRAWS_PER_MARKER*nprt) then
              !> The largest momentum actually present in the population is the
              !> quantity that decides whether the condition pnorm > p_far_tmp can
              !> be met at all. Once the field goes three dimensional the markers
              !> cool, the population maximum falls below p_far, and no draw can
              !> ever succeed. Reporting it turns this abort into a diagnosis.
              write(*,*) ">>> ERROR : the reweighting loop reached ", MAX_DRAWS_PER_MARKER, &
                         "*nprt draws without finding a duplication candidate."
              write(*,*) "    p_far required        :", p_far_tmp
              write(*,*) "    largest p/(m c) present:", maxval(p_before_reweight)
              write(*,*) "    psi_old, psi_res_tmp  :", psi_old, psi_res_tmp
              write(*,*) "    The population has cooled below p_far, so the condition is"
              write(*,*) "    unsatisfiable. Enable use_p_far_sweep, or lower p_far and"
              write(*,*) "    p_cutoff_recycle to match the momenta actually present."
              call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
            end if

            new_iprt = draw_marker_index(nprt)
            if (new_iprt .eq. iprt) cycle
            if (prt(new_iprt)%i_elm .lt. 1) cycle
            if (prt(new_iprt)%i_elm .gt. fields%element_list%n_elements) cycle

            call fields%calc_EBpsiU(time, prt(new_iprt)%i_elm, prt(new_iprt)%st, &
                                    prt(new_iprt)%x(3), E, B, psi_new, U, E_par_norm=.false.)
            call gc_momentum_pitch(prt(new_iprt)%p(1), prt(new_iprt)%p(2), norm2(B), mass, pnorm, xi)

            !> Duplication conditions
            is_new_part = (pnorm .gt. p_far_tmp)            .and. &
                          (prt(new_iprt)%weight .gt. 0.d0)  .and. &
                          (prt_age(new_iprt) .gt. age_split)      !< avoids clustering
            if (use_psi_res) then
              !> Keeps the duplication local in psi, which limits the spatial
              !> displacement introduced by the reweighting
              is_new_part = is_new_part .and. (abs(psi_old - psi_new) .lt. psi_res_tmp)
            end if

            !> Relax the search window every nprt unsuccessful draws
            if (mod(i_loop_counter, nprt) .eq. 0) then
              if (use_psi_res)     psi_res_tmp = psi_res_tmp + delta_psi_res
              if (use_p_far_sweep) p_far_tmp   = p_far_tmp   - delta_p_far
            end if
          end do

          prt(iprt)%x       = prt(new_iprt)%x
          prt(iprt)%p       = prt(new_iprt)%p
          prt(iprt)%st      = prt(new_iprt)%st
          prt(iprt)%i_elm   = prt(new_iprt)%i_elm
          prt(iprt)%t_birth = prt(new_iprt)%t_birth
          prt(iprt)%i_life  = prt(new_iprt)%i_life + 1

          !> The donor weight is shared between the two copies
          new_weight           = 0.5d0 * prt(new_iprt)%weight
          prt(iprt)%weight     = new_weight
          prt(new_iprt)%weight = new_weight

          prt_age(iprt)     = -1
          prt_age(new_iprt) = 0
        end if

        !> The barrel click advances on candidates only, so that exactly l out of
        !> every k candidates survive
        i_click = i_click + 1
        if (i_click .gt. k) then
          call generate_barrel(barrel, k, l)
          i_click = 1
        end if
      end if

      prt_age(iprt)           = prt_age(iprt) + 1
      p_before_reweight(iprt) = pnorm
    end do
  end select

  deallocate(barrel)
end subroutine correlated_russian_roulette

!> Fills a barrel of k slots with exactly l ones placed at random positions.
!> The Fisher-Yates shuffle is used because the naive "swap slot i with a random
!> slot" loop does not produce a uniformly distributed permutation.
subroutine generate_barrel(barrel, k, l)
  implicit none
  integer, intent(in)    :: k       !< barrel size
  integer, intent(in)    :: l       !< number of surviving slots
  integer, intent(inout) :: barrel(k)

  integer :: i, j, tmp
  real*8  :: rnd

  barrel      = 0
  barrel(1:l) = 1

  do i = k, 2, -1
    call random_number(rnd)
    j = min(int(rnd * real(i, 8)) + 1, i)
    tmp       = barrel(i)
    barrel(i) = barrel(j)
    barrel(j) = tmp
  end do
end subroutine generate_barrel

!> Replaces marker iprt by a copy of a randomly drawn marker above p_cutoff, the
!> two copies sharing the donor weight. This keeps the marker number constant
!> after a loss without changing the total weight carried by the population.
subroutine recycle_lost(particles, fields, nprt, iprt, time, mass, p_cutoff,  &
                        p_before_reweight, min_weight, &
                        n_write_MPI, T_write_MPI, E_write_MPI, ifail)
  use, intrinsic :: ieee_arithmetic
  implicit none
  class(particle_base), dimension(:), intent(inout) :: particles
  class(fields_base),                 intent(in)    :: fields
  integer,                            intent(in)    :: nprt  !< number of local markers
  integer,                            intent(in)    :: iprt  !< marker to be replaced
  real*8,                             intent(in)    :: time, mass
  real*8,                             intent(in)    :: p_cutoff, min_weight
  real*8,                             intent(inout) :: p_before_reweight(nprt)
  real*8,                             intent(inout) :: n_write_MPI(:), T_write_MPI(:), E_write_MPI(:)
  integer,                            intent(out)   :: ifail !< 0 on success, 1 if no donor was found

  !> Draws after which the momentum condition on the donor is dropped
  integer, parameter :: RELAX_DRAWS_PER_MARKER = 3
  !> Draws after which the search is declared failed
  integer, parameter :: MAX_DRAWS_PER_MARKER = 20

  integer :: i_draw, new_iprt
  real*8  :: pnorm, xi, new_weight
  real*8  :: E(3), B(3), psi, U
  logical :: accepted

  ifail    = 0
  accepted = .false.
  new_iprt = 0
  i_draw   = 0
  pnorm    = 0.d0

  select type (prt => particles)
  type is (particle_gc_relativistic)
    do while (.not. accepted)
      !> The counter is incremented on every draw. Incrementing it only for
      !> candidates that pass the validity test makes the loop non-terminating
      !> when no marker passes it.
      i_draw = i_draw + 1
      if (i_draw .gt. MAX_DRAWS_PER_MARKER*nprt) then
        ifail = 1
        return
      end if

      new_iprt = draw_marker_index(nprt)

      if (new_iprt .eq. iprt)                     cycle
      !> Both bounds are tested: interp_PRZ indexes element_list%element(i_elm)
      !> with no bound check of its own, so a stale index is an out-of-range
      !> access rather than a detectable error.
      if (prt(new_iprt)%i_elm .lt. 1)             cycle
      if (prt(new_iprt)%i_elm .gt. fields%element_list%n_elements) cycle
      if (prt(new_iprt)%weight .le. min_weight)   cycle
      if (ieee_is_nan(prt(new_iprt)%p(1)))        cycle

      !> Only B is needed here, so the E = eta*J_par normalisation is skipped
      call fields%calc_EBpsiU(time, prt(new_iprt)%i_elm, prt(new_iprt)%st, &
                              prt(new_iprt)%x(3), E, B, psi, U, E_par_norm=.false.)
      call gc_momentum_pitch(prt(new_iprt)%p(1), prt(new_iprt)%p(2), norm2(B), mass, pnorm, xi)

      if (i_draw .lt. RELAX_DRAWS_PER_MARKER*nprt) then
        accepted = (pnorm .gt. p_cutoff)
      else
        !> Relax the momentum condition rather than stall
        accepted = .true.
      end if
    end do

    prt(iprt)%x       = prt(new_iprt)%x
    prt(iprt)%p       = prt(new_iprt)%p
    prt(iprt)%st      = prt(new_iprt)%st
    prt(iprt)%i_elm   = prt(new_iprt)%i_elm
    prt(iprt)%t_birth = prt(new_iprt)%t_birth
    prt(iprt)%i_life  = prt(new_iprt)%i_life + 1

    n_write_MPI(iprt) = n_write_MPI(new_iprt)
    T_write_MPI(iprt) = T_write_MPI(new_iprt)
    E_write_MPI(iprt) = E_write_MPI(new_iprt)

    !> The donor weight is shared between the two copies
    new_weight           = 0.5d0 * prt(new_iprt)%weight
    prt(iprt)%weight     = new_weight
    prt(new_iprt)%weight = new_weight

    p_before_reweight(iprt) = pnorm
  end select
end subroutine recycle_lost

end module testccoll_helpers

!############################################################################################
!#############################         SIMULATION LOOP        ###############################
!############################################################################################
program hot_tail_track

  use particle_tracer
  use mod_event
  use mod_particle_io
  use mod_particle_diagnostics
  use mod_fields_linear
  use mod_fields_hermite_birkhoff
  use mod_gc_relativistic
  use mod_ccoll_relativistic
  use mod_radreactforce
  use mod_kinetic_relativistic
  use mod_project_particles
  use mod_impurity, only: init_imp_adas
  use testccoll_helpers
  use, intrinsic :: ieee_arithmetic
  use constants
  use mod_bessel, only : bessel_k2exp
  use mpi

  implicit none

!---------------------------------------!
! Set up the simulation variables
!> dt_marker is the marker time step in seconds. It is deliberately not called
!> tstep, which in JOREK denotes the MHD time step in JOREK units.
real*8      :: dt_marker, deltat, deltat_write, deltat_reweight, duration, mass
real*8      :: w_jump, total_weight, range_weight, min_weight, sigma_psi

real*8      :: timesteps(1)
real*8      :: time_startsim, target_time, next_time_write, next_time_reweight
integer*4   :: nstep, nstep_write, nstep_reweight
!> n_mpi instead of size: the latter shadows the intrinsic SIZE
integer     :: rank, n_mpi, ierr
integer     :: start_iprt, end_iprt, n_jump, n_invalid, n_recycle_fail
integer     :: n_elm_max, n_stale, n_stale_total
integer     :: n_loss_buffer, n_loss_buffer_max, max_loss_per_rank
integer     :: k, l, age_split, iseed, n_seed, i_seed
integer, allocatable  :: seed_array(:)
real*8,dimension(2)   :: Rbound, Zbound, Phibound, Ekinbound, Pitchbound
logical     :: restart_particles, average_coll_fields, file_exists
type(event) :: fieldreader

! Variables for PDF initialization
type(pcg32_rng)                         :: rng_pcg32
integer                                 :: n_int_pdf_param, n_real_pdf_param
integer                                 :: n_int_weight_param,n_real_weight_param
integer                                 :: n_variables
integer                                 :: n_int_pdf_to_part_coord_param, n_real_pdf_to_part_coord_param
integer                                 :: n_int_gdf_param,n_real_gdf_param
real*8                                  :: pdf_upper_bound, gdf_upper_bound, charge, sup_pdf_safety_factor
real*8                                  :: psi_axis, R0_axis, Z0_axis, s_axis, t_axis
integer*4                               :: i_elm, ifail
real*8,dimension(2)                     :: Rbox,Zbox
real*8,dimension(2)                     :: Pbound
real*8,dimension(6,2)                   :: phase_space_bounds
integer,dimension(:),allocatable        :: int_pdf_param, int_weight_param, int_gdf_param
integer,dimension(:),allocatable        :: int_pdf_to_part_coord_param
real*8,dimension(:),allocatable         :: real_pdf_to_part_coord_param
real*8,dimension(:),allocatable         :: real_pdf_param, real_weight_param,real_gdf_param
procedure(real_f),pointer               :: pdf_to_use         => NULL()
procedure(real_f),pointer               :: weight_to_use      => NULL()
procedure(real_f),pointer               :: gdf_to_use         => NULL()
procedure(real_arr_inout_s),pointer     :: gdf_sampler_to_use => NULL()
procedure(part_inout_s),pointer         :: pdf_to_part_coord  => NULL()

character(len=200) :: fnout_collparam, fnout_losses !< File where the collision / losses outputs are written

real*8, allocatable :: n_write(:), T_write(:), E_write(:)
real*8, allocatable :: n_write_MPI(:), T_write_MPI(:), E_write_MPI(:)
real*8, allocatable :: t_losses(:), p_losses(:), xi_losses(:), R_losses(:), Z_losses(:), Phi_losses(:), Weight_losses(:), iprt_losses(:)
real*8, allocatable :: p_before_reweight(:)
integer*4, allocatable :: prt_age(:)
real*8, allocatable :: t_losses_tmp(:), p_losses_tmp(:), xi_losses_tmp(:), R_losses_tmp(:), Z_losses_tmp(:), Phi_losses_tmp(:), Weight_losses_tmp(:), iprt_losses_tmp(:)
!> lost          : the marker left the domain and must be written to the loss file
!> needs_recycle : the marker must be replaced by a clone (physical loss, numerical
!>                 jump or invalid state). A numerical jump is NOT a physical loss.
logical, allocatable :: lost(:), needs_recycle(:)

integer*4 :: iprt, nprt, nprt_MPI, istep, chargenum, n_lost

real*8  :: E(3), B(3), psi, U, p_after_RK4, delta_jump, jump_coll, E_par_max, t_push
real*8  :: p_before, xi_before, xi_after_RK4, R_before, Z_before, Phi_before
real*8  :: pnorm
real*8, dimension(:), allocatable    :: ni
real*8  :: the, ne, ti
real*8  :: p_cutoff, p_far, p_cutoff_recycle
logical :: use_psi_res, use_p_far_sweep
real*8  :: psi_res, delta_psi_res, delta_p_far

type(write_particle_diagnostics)                  :: diag

! For timing
real*8 :: t0, t1, t_wall0, t_wall1
! Collision data
type(ccoll_data) :: dat
!---------------------------------------!

!<<<< INPUT : Simulation options >>>>!
dt_marker     = 10.e-9           !< Marker time step [s]
duration      = 3.0e-3           !< How long markers are traced [s]
nprt          = 30000            !< Number of markers
time_startsim = 2.0e-3           !< Time to start simulation (not accounted if restart particles)
nstep_write = 100                !< Number of writting time step for diagnostic (no counting the t=0 one)
nstep_reweight = 100             !< Number of reweighting steps using the russian roulette method
average_coll_fields = .TRUE.     !< Average the fields observed by a particle between two writting times (usefull when E is very noisy!) <!> Must be improved: at t = t - dt/2 | no reweighting consideration
restart_particles = .FALSE.      !< Use a 'part_restart.h5' file to initialize markers
iseed         = 0                !< 0: seed the RNG from the OS entropy (results differ from run to run)
                                 !< >0: reproducible per-rank seed derived from iseed. Note that gfortran
                                 !<     derives the per-thread states from the master state, so results
                                 !<     still depend on the number of OpenMP threads.

! Initialization and reweighting parameters
Rbound          = [6.3692-1.8110, 6.3692+1.8110]      !< R coordinate initialization window
Zbound          = [-3.,4.]           !< Z coodinate initialization window
Phibound        = [TWOPI, 2.*TWOPI]   !< Toroidal angle initialization window
Ekinbound       = [10., 410.0d+3]     !< kinetic energy initialization window [eV]
Pitchbound      = [1.d-1, PI-1.d-1]   !< Pitch angle initialization window [rad] (sampled uniformly in its cosine)
total_weight    = 8.4d22      !< Total weight of markers
range_weight    = 1.d15           !< Range of weights
sigma_psi       = 1000.0       !< Standard deviation of Pi function for spatial initialization (JOREK units of Psi)
jump_coll       = 3.0          !< Tolerance on momentum variation between two timesteps for the collision push (in order to avoid unphysical jumps). Jumpi if: p(t)>p_th and p(t+dt) > jump_coll*p_th
E_par_max       = 1.d4         !< Upper bound assumed on |E_par| [V/m], used to build the RK4 jump detection threshold
p_cutoff = 0.2       !< cutoff momentum for the correlated russian roulette (set a negative value if no reweighting is wanted)
p_far = 0.6          !< far momentum over which duplication of markers is allowed for the correlated russian roulette
k = 10               !< Size of the barrel for the russian roulette
l = 1                !< (k,l) define the save probability : P_sav = l/k
age_split = 1        !< marker minimum age for being duplicated
use_psi_res     = .FALSE.  !< True: activate the condition of duplication close to psi_res from the death of the killed marker (may help to avoid clustering)
use_p_far_sweep = .FALSE.   !< True: allows the p_far to be reduced by delta_p_far if no candidate is found after nprt draws (activate if finding duplication candidate is too hard/impossible)
psi_res         = 0.5   !< if use_psi_res=.TRUE., define the psi distance for the duplication condition
delta_psi_res   = 0.5   !< if not duplication candidate found after nprt draws, psi_res = psi_res + delta_psi_res
delta_p_far     = 0.1   !< if use_p_far_sweep=.TRUE., if not duplication candidate found after nprt draws, p_far = p_far - delta_p_far
p_cutoff_recycle = 0.6  !< cutoff momentum for the recycling of losses
max_loss_per_rank = 200000 !< Maximum number of buffered loss records per rank before an early
                           !< flush is forced. write_loss_to_hdf5 gathers every rank's buffer on
                           !< rank 0, so the root allocation is n_mpi * this * 8 arrays * 8 bytes,
                           !< i.e. 0.41 GB at 32 ranks with the default. Without a cap the buffers
                           !< grow with the loss rate times the number of steps between two
                           !< scheduled writes, which reaches several GB during a thermal quench.


! Naming of output files
write(fnout_collparam, "(A, F0.0, A, F0.0, A, I0, A, F0.0, A)") &
        "Jcoll_K", Ekinbound(1)*1e-3, "_", Ekinbound(2)*1e-3, &
        "_np", nprt, "_dt", dt_marker*1.e+9, "ns.h5"

write(fnout_losses, "(A, F0.0, A, F0.0, A, I0, A, F0.0, A)") &
        "Jlosses_K", Ekinbound(1)*1e-3, "_", Ekinbound(2)*1e-3, &
        "_np", nprt, "_dt", dt_marker*1.e+9, "ns.h5"

!> Electron mass in atomic mass units
mass = MASS_ELECTRON / ATOMIC_MASS_UNIT
chargenum = -1
!<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>!

! Initialization Collision Operator and Ions/Imp
call sim%initialize(num_groups=1)
call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
call MPI_Comm_size(MPI_COMM_WORLD, n_mpi, ierr)

!> Random number generator seeding
if (iseed .eq. 0) then
  call random_seed()
else
  call random_seed(size = n_seed)
  allocate(seed_array(n_seed))
  seed_array = [ ( iseed + 37*(i_seed-1) + 1013*rank, i_seed = 1, n_seed ) ]
  call random_seed(put = seed_array)
  deallocate(seed_array)
end if

call init_imp_adas(0)

!> The 'ccolldata' file must already be present
inquire(file='ccolldata', exist=file_exists)
if (.not. file_exists) then
  if (rank .eq. 0) then
    write(*,*) "ERROR: the tabulated collision data file 'ccolldata' was not found"
    write(*,*) "       in the run directory. It has to be generated once, by running"
    write(*,*) "       the program ccoll_generate_L0L1 (particles/util/ccoll_generate_L0L1.f90)."
  end if
  call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
end if
call ccoll_init('ccolldata', dat)
allocate(ni(ubound(dat%Ii, dim=1)))

! Define simulation and writting steps
nstep  = nint(duration/dt_marker)
deltat = duration / nstep
deltat_write = duration / nstep_write
deltat_reweight = duration / (nstep_reweight + 1.0)
next_time_write = 0.
next_time_reweight = deltat_reweight
timesteps = [deltat]
n_jump = 0
n_invalid = 0
n_recycle_fail = 0
n_stale_total = 0
w_jump = 0.d0
!> Global weight floor, set once from the total weight. Markers below it are
!> replaced, and they are not accepted as recycling donors.
min_weight = total_weight / range_weight

! Allocate and initialize collisions parameters
allocate(T_write(nprt), E_write(nprt), n_write(nprt))
T_write = 0.d0
E_write = 0.d0
n_write = 0.d0

! MPI Initialization
call MPI_tracker_init(nprt, rank, n_mpi, ierr, start_iprt, end_iprt, nprt_MPI, &
                      n_write_MPI, T_write_MPI, E_write_MPI)
call MPI_Barrier(MPI_COMM_WORLD, ierr)

if (rank .eq. 0) then
   write(*,*) "  "
   write(*,*) "####################### Inputs #######################"
   write(*,*) "dt_marker = ", dt_marker
   write(*,*) "duration =", duration
   write(*,*) "nprt =", nprt
   write(*,*) "time_startsim =", time_startsim
   write(*,*) "nstep_write =", nstep_write
   write(*,*) "nstep_reweight =", nstep_reweight
   write(*,*) "average_coll_fields =", average_coll_fields
   write(*,*) "restart_particles =", restart_particles
   write(*,*) "iseed =", iseed
   write(*,*) "Rbound =", Rbound
   write(*,*) "Zbound =", Zbound
   write(*,*) "Phibound =", Phibound
   write(*,*) "Ekinbound =", Ekinbound
   write(*,*) "Pitchbound =", Pitchbound
   write(*,*) "total_weight =", total_weight
   write(*,*) "range_weight =", range_weight
   write(*,*) "sigma_psi =", sigma_psi
   write(*,*) "jump_coll =", jump_coll
   write(*,*) "E_par_max =", E_par_max
   write(*,*) "p_cutoff =", p_cutoff
   write(*,*) "p_far =", p_far
   write(*,*) "k =", k
   write(*,*) "l =", l
   write(*,*) "age_split =", age_split
   write(*,*) "use_psi_res =", use_psi_res
   write(*,*) "use_p_far_sweep =", use_p_far_sweep
   write(*,*) "psi_res =", psi_res
   write(*,*) "delta_psi_res =", delta_psi_res
   write(*,*) "delta_p_far =", delta_p_far
   write(*,*) "p_cutoff_recycle =", p_cutoff_recycle
   write(*,*) "max_loss_per_rank =", max_loss_per_rank
   write(*,*) "#################################################"
   write(*,*) "  "
end if
call MPI_Barrier(MPI_COMM_WORLD, ierr)

! Set up the diagnostics output
diag = write_particle_diagnostics(filename='diag.h5', only=[5,7,8,12,13,14,15,16,17], append=.false.)
call create_hdf5(fnout_collparam)
call create_hdf5(fnout_losses)

!#### Initialization of markers: from a PDF or a restart file ####
if (restart_particles) then
   !> Restore the markers from 'part_restart.h5'. read_simulation_hdf5 reallocates
   !> sim%groups(1)%particles to this rank's share of the marker count stored in
   !> the file and restores x, st, p, weight, i_elm, i_life, t_birth and q, as
   !> well as sim%time and the group mass. 
   if (rank .eq. 0) then
      write(*,*) '*******************************'
      write(*,*) 'Reading particle restart file'
   end if
   call read_simulation_hdf5(sim, 'part_restart.h5')
   if (rank .eq. 0) then
      write(*,*) 'Time = ', sim%time
      write(*,*) '*******************************'
   end if
   time_startsim = sim%time

   !> Begin the simulation at the restored time slice by loading the field there
   fieldreader = event(read_jorek_fields_interp_linear(i=last_file_before_time(sim%time), stop_at_end=.true.))
   events = [fieldreader]
   call with(sim, events)

else
   ! Initialize markers using PDF Initializer
   sim%time = time_startsim
   sim%groups(1)%mass = mass

   ! Begin simulation at new time slice by initializing the field at that point
   fieldreader = event(read_jorek_fields_interp_linear(i=last_file_before_time(sim%time), stop_at_end=.true.))
   !fieldreader = event(read_jorek_fields_interp_hermite_birkhoff(i=last_file_before_time(sim%time)))
   events = [fieldreader]
   call with(sim,events)

   ! Use relativistic GCs
   allocate(particle_gc_relativistic::sim%groups(1)%particles(nprt_MPI))
   ! Initialize each sub-population of markers for each MPI process
   n_variables     = 6
   charge          = chargenum * 1.d0

   ! Magnetic Axis coordinates
   call find_axis(rank, sim%fields%node_list, sim%fields%element_list, psi_axis, R0_axis, Z0_axis, i_elm, s_axis, t_axis, ifail)

   ! Extra real parameter of pdf and weights
   n_real_pdf_param        = 1
   n_real_weight_param     = 1
   n_real_pdf_to_part_coord_param = 1
   allocate(real_pdf_param(n_real_pdf_param))
   allocate(real_weight_param(n_real_weight_param))
   allocate(real_pdf_to_part_coord_param(n_real_pdf_to_part_coord_param))
   real_pdf_param(1) = sim%groups(1)%mass
   real_weight_param(1) = sim%groups(1)%mass
   real_pdf_to_part_coord_param(1) = sim%groups(1)%mass

   ! Extra integer parameter of pdf and weights
   n_int_pdf_param         = 0
   n_int_weight_param      = 0
   n_int_pdf_to_part_coord_param  = 0
   allocate(int_pdf_param(n_int_pdf_param))
   allocate(int_weight_param(n_int_weight_param))
   allocate(int_pdf_to_part_coord_param(n_int_pdf_to_part_coord_param))

   ! Creation of phase space bounds
   sup_pdf_safety_factor   = 1d0
   call domain_bounding_box(sim%fields%node_list,sim%fields%element_list,Rbox(1),Rbox(2),Zbox(1),Zbox(2))
   if(Rbox(1).ge.Rbound(1)) Rbound(1) = Rbox(1)
   if((Rbox(2).lt.Rbound(2)).and.((Rbox(2)-Rbound(1)).gt.0.d0)) Rbound(2) = Rbox(2)
   if((Zbox(1).gt.0.d0).and.(Zbound(1).ge.Zbox(1))) Zbound(1) = Zbox(1)
   if((Zbox(1).lt.0.d0).and.(Zbound(1).lt.Zbox(1))) Zbound(1) = Zbox(1)
   if((Zbox(2).gt.0.d0).and.(Zbound(2).ge.Zbox(2)).and.((Zbox(2)-Zbound(1)).gt.0.d0)) Zbound(2) = Zbox(2)
   if((Zbox(2).lt.0.d0).and.(Zbound(2).lt.Zbox(2)).and.((Zbox(2)-Zbound(1)).gt.0.d0)) Zbound(2) = Zbox(2)
   Pbound = sim%groups(1)%mass*SPEED_OF_LIGHT*sqrt(((EL_CHG*Ekinbound/(ATOMIC_MASS_UNIT*sim%groups(1)%mass*SPEED_OF_LIGHT**2))+1.d0)**2-1.d0)
   phase_space_bounds(:,1) = [Rbound(1),Zbound(1),Phibound(1),Pbound(1),Pitchbound(1),charge]
   phase_space_bounds(:,2) = [Rbound(2),Zbound(2),Phibound(2),Pbound(2),Pitchbound(2),charge]

   ! PDF
   pdf_to_use        => pdf_uniform
   ! SUP PDF
   pdf_upper_bound   = sup_pdf_uniform(n_variables,&
   phase_space_bounds(1:n_variables,1),phase_space_bounds(1:n_variables,2),&
   n_real_pdf_param,real_pdf_param,n_int_pdf_param,int_pdf_param)
   ! WEIGHT
   weight_to_use     => weight_MJ
   ! Coordinate computation
   pdf_to_part_coord => spherical_p_cartesian_q_to_relativistic_gc
   ! GDF
   n_real_gdf_param = 0; n_int_gdf_param = 0;
   gdf_to_use         => gdf_uniform_phase
   gdf_sampler_to_use => gdf_uniform_sampler
   gdf_upper_bound    = sup_gdf_uniform_phase(n_variables,&
   phase_space_bounds(1:n_variables,1),phase_space_bounds(1:n_variables,2), &
   n_real_pdf_param,real_pdf_param,n_int_pdf_param,int_pdf_param)

   ! Initialization
   if (rank .eq. 0) then
      write(*,*) "... initialising guiding center in phase space"
   end if
   call initialise_particles_in_phase_space(n_variables,sim%groups(1)%particles,sim%fields,rng_pcg32,&
   pdf_to_use,weight_to_use,gdf_to_use,gdf_sampler_to_use,pdf_upper_bound,gdf_upper_bound,&
   pdf_to_part_coord,sim%groups(1)%mass,sim%time,phase_space_bounds,n_real_pdf_param,real_pdf_param,&
   n_int_pdf_param,int_pdf_param,n_real_weight_param,real_weight_param,n_int_weight_param,&
   int_weight_param,n_real_gdf_param,real_gdf_param,n_int_gdf_param,int_gdf_param,&
   n_real_pdf_to_part_coord_param,real_pdf_to_part_coord_param,&
   n_int_pdf_to_part_coord_param,int_pdf_to_part_coord_param)

   do iprt=1, nprt_MPI
      if (ieee_is_nan(sim%groups(1)%particles(iprt)%weight) .OR. &
          (.NOT. ieee_is_finite(sim%groups(1)%particles(iprt)%weight)) .OR. &
          (abs(sim%groups(1)%particles(iprt)%weight) .le. 0.d0)) then
         sim%groups(1)%particles(iprt)%weight = 0.d0
         sim%groups(1)%particles(iprt)%i_elm = 0
      end if
   end do
   call adjust_particle_weights(sim%groups(1)%particles, total_weight)
end if

! Allocate the reweighting and loss bookkeeping
allocate(p_before_reweight(nprt_MPI), prt_age(nprt_MPI))
allocate(t_losses(0), p_losses(0), xi_losses(0), R_losses(0), Z_losses(0), Phi_losses(0), Weight_losses(0), iprt_losses(0))
allocate(t_losses_tmp(nprt_MPI))
allocate(p_losses_tmp(nprt_MPI))
allocate(xi_losses_tmp(nprt_MPI))
allocate(R_losses_tmp(nprt_MPI))
allocate(Z_losses_tmp(nprt_MPI))
allocate(Phi_losses_tmp(nprt_MPI))
allocate(Weight_losses_tmp(nprt_MPI))
allocate(iprt_losses_tmp(nprt_MPI))
allocate(lost(nprt_MPI), needs_recycle(nprt_MPI))
t_losses_tmp = 0d0
p_losses_tmp = 0d0
xi_losses_tmp = 0d0
R_losses_tmp = 0d0
Z_losses_tmp = 0d0
Phi_losses_tmp = 0d0
Weight_losses_tmp = 0d0
iprt_losses_tmp = 0d0
lost          = .false.
needs_recycle = .false.

! First census for the reweighting method
select type (prt => sim%groups(1)%particles)
type is (particle_gc_relativistic)
   do iprt=1,nprt_MPI
      if (prt(iprt)%i_elm .lt. 1) then
         p_before_reweight(iprt) = 0.d0
      else
         call sim%fields%calc_EBpsiU(sim%time, prt(iprt)%i_elm, prt(iprt)%st, prt(iprt)%x(3), &
                                     E, B, psi, U, E_par_norm=.FALSE.)
         call gc_momentum_pitch(prt(iprt)%p(1), prt(iprt)%p(2), norm2(B), sim%groups(1)%mass, &
                                p_before_reweight(iprt), xi_before)
      end if
      prt_age(iprt) = 10
   end do
end select
call MPI_Barrier(MPI_COMM_WORLD, ierr)

! Replace the markers that did not survive the initialisation, then compute the
! first collision fields
select type (prt => sim%groups(1)%particles)
type is (particle_gc_relativistic)
   !> On a restart the marker weights come from the restart file and are not to
   !> be second-guessed, so only an invalid element index calls for a replacement.
   if (restart_particles) then
      needs_recycle = (prt(1:nprt_MPI)%i_elm .lt. 1)
   else
      needs_recycle = (prt(1:nprt_MPI)%i_elm .lt. 1) .OR. (prt(1:nprt_MPI)%weight .lt. min_weight)
   end if
   do iprt=1,nprt_MPI
      if (needs_recycle(iprt)) then
         call recycle_lost(sim%groups(1)%particles, sim%fields, nprt_MPI, iprt, sim%time, &
                           sim%groups(1)%mass, p_cutoff_recycle, p_before_reweight, &
                           min_weight, n_write_MPI, T_write_MPI, E_write_MPI, ifail)
         if (ifail .ne. 0) then
            write(*,*) "ERROR: no valid marker left to recycle marker", iprt, "on rank", rank
            write(*,*) "       The whole local population is invalid; check the initialisation window."
            call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
         end if
         needs_recycle(iprt) = .false.
      end if
      call sim%fields%calc_EBpsiU(sim%time, prt(iprt)%i_elm, prt(iprt)%st, prt(iprt)%x(3), &
                                  E, B, psi, U, E_par_norm=.TRUE.)
      call sim%fields%calc_NjTj(sim%time, prt(iprt)%i_elm, prt(iprt)%st, prt(iprt)%x(3), &
                                dat%m_i_over_m_imp, ne, the, ni, ti)
      n_write_MPI(iprt) = ne
      T_write_MPI(iprt) = the*K_BOLTZ/EL_CHG
      E_write_MPI(iprt) = dot_product(E,B)/norm2(B)
   end do
end select
needs_recycle = .false.

! Set dpsi/dt=0
!call sim%fields%set_flag_dpsidt(.true.)
call check_and_fix_timesteps(timesteps, events)
target_time = next_event_at(sim, events)

!> Threshold used to detect a non-physical momentum jump during the RK4 push.
!> The parallel electric field cannot change p/(m c) by more than
!> e * E_par_max * dt / (m_e c) over one step, so a larger variation signals a
!> numerical artefact. Numerically this is 586.7 * E_par_max * dt.
delta_jump = EL_CHG * E_par_max * dt_marker / (mass * ATOMIC_MASS_UNIT * SPEED_OF_LIGHT)

! Define Start CPU and wall clock times
call cpu_time(t0)
t_wall0 = MPI_Wtime()
! #######################################################################################
! ############################### START SIMULATION LOOP #################################
! #######################################################################################
if (rank .eq. 0) then
   write(*,*) ">>>>>>>>>>>>>> START SIMULATION LOOP <<<<<<<<<<<<<<"
end if
do istep=1,nstep
   !> Diag output
   if (sim%time - time_startsim .ge. next_time_write) then
      call with(sim, diag)    ! tracking diagnostic
      call write_simulation_hdf5(sim, 'part_restart.h5')   ! restart file for particles
      if (n_recycle_fail .gt. 0) then
         write(*,*) 'WARNING: failed recycling attempts so far =', n_recycle_fail, &
                    ', invalid markers at push entry =', n_invalid, ' / rank =', rank
         flush(6)
      end if
      call write_collparam_hdf5(fnout_collparam, n_write, T_write, E_write, n_write_MPI, T_write_MPI, E_write_MPI, nprt_MPI, sim%time)   ! field diagnostic
      call write_loss_to_hdf5(fnout_losses, t_losses, p_losses, xi_losses, R_losses, Z_losses, Phi_losses, Weight_losses, iprt_losses)  ! losses diagnostic
   end if

   select type (prt => sim%groups(1)%particles)
   type is (particle_gc_relativistic)

      !$omp parallel do default(shared) &
      !$omp firstprivate(ni) &
      !$omp private(pnorm, xi_before) &
      !$omp private(the, ti, ne, E, B, psi, U, p_after_RK4, xi_after_RK4, t_push) &
      !$omp private(p_before, R_before, Z_before, Phi_before) &
      !$omp private(iprt) &
      !$omp reduction(+:n_jump, w_jump, n_invalid)
      do iprt = 1, nprt_MPI

         if (prt(iprt)%i_elm .lt. 1 .or. &
             prt(iprt)%i_elm .gt. sim%fields%element_list%n_elements) then
            n_invalid = n_invalid + 1
            needs_recycle(iprt) = .TRUE.
            cycle
         end if

         ! Keep data information before the push
         !> Only B is needed here, so the E = eta*J_par normalisation is skipped
         call sim%fields%calc_EBpsiU(sim%time, prt(iprt)%i_elm, prt(iprt)%st, prt(iprt)%x(3), &
                                     E, B, psi, U, E_par_norm=.FALSE.)
         call gc_momentum_pitch(prt(iprt)%p(1), prt(iprt)%p(2), norm2(B), sim%groups(1)%mass, &
                                p_before, xi_before)
         R_before = prt(iprt)%x(1)
         Z_before = prt(iprt)%x(2)
         Phi_before = prt(iprt)%x(3)

         !>>> RK4 push <<<
         t_push = sim%time
         call runge_kutta_fixed_dt_gc_push_jorek_radreact(sim%fields,t_push,timesteps(1), &
                                                          sim%groups(1)%mass, prt(iprt), E_par_norm=.TRUE.)

         if ((prt(iprt)%i_elm .lt. 1) .OR. (ieee_is_nan(prt(iprt)%p(1))) &
                                      .OR. (ieee_is_nan(prt(iprt)%p(2)))) then
            !> Physical loss: register it and flag the marker for recycling
            t_losses_tmp(iprt) = sim%time
            p_losses_tmp(iprt) = p_before
            xi_losses_tmp(iprt) = xi_before
            R_losses_tmp(iprt) = R_before
            Z_losses_tmp(iprt) = Z_before
            Phi_losses_tmp(iprt) = Phi_before
            Weight_losses_tmp(iprt)  = prt(iprt)%weight
            iprt_losses_tmp(iprt)  = real(iprt + start_iprt - 1)
            lost(iprt)  = .TRUE.
            needs_recycle(iprt) = .TRUE.
         else
            ! Jump Check for the RK4 push
            !> E is needed below only when the collision fields are averaged
            call sim%fields%calc_EBpsiU(sim%time, prt(iprt)%i_elm, prt(iprt)%st, prt(iprt)%x(3), &
                                        E, B, psi, U, E_par_norm=average_coll_fields)
            call gc_momentum_pitch(prt(iprt)%p(1), prt(iprt)%p(2), norm2(B), sim%groups(1)%mass, &
                                   p_after_RK4, xi_after_RK4)
            if ((ABS(p_after_RK4-p_before) .gt. delta_jump) .OR. (prt(iprt)%weight .lt. 0.d0)) then
               !> Numerical artefact, not a physical loss: the marker is replaced
               !> at the end of the step and is NOT written to the loss file.
               !> Registering it there would pollute the loss diagnostic.
               n_jump = n_jump + 1
               w_jump = w_jump + prt(iprt)%weight
               needs_recycle(iprt) = .TRUE.
            else
               !>>> Collision Operator <<<
               call ccoll_gc_relativistic_push(dat, prt(iprt), sim%fields, sim%groups(1)%mass, sim%time, timesteps(1), jump_coll=jump_coll)
               !> The collision push sets i_elm = 0 when it has to discard a marker
               if (prt(iprt)%i_elm .lt. 1) then
                  t_losses_tmp(iprt) = sim%time
                  p_losses_tmp(iprt) = p_after_RK4
                  xi_losses_tmp(iprt) = xi_after_RK4
                  R_losses_tmp(iprt) = prt(iprt)%x(1)
                  Z_losses_tmp(iprt) = prt(iprt)%x(2)
                  Phi_losses_tmp(iprt) = prt(iprt)%x(3)
                  Weight_losses_tmp(iprt) = prt(iprt)%weight
                  iprt_losses_tmp(iprt) = real(iprt + start_iprt - 1)
                  lost(iprt) = .TRUE.
                  needs_recycle(iprt) = .TRUE.
               end if
            end if

            if (average_coll_fields .AND. (.NOT. needs_recycle(iprt))) then
               call sim%fields%calc_NjTj(sim%time, prt(iprt)%i_elm, prt(iprt)%st, prt(iprt)%x(3), dat%m_i_over_m_imp, ne, the, ni, ti)
               n_write_MPI(iprt) = n_write_MPI(iprt) + (ne / (deltat_write/timesteps(1)))
               T_write_MPI(iprt) = T_write_MPI(iprt) + (the*K_BOLTZ/EL_CHG / (deltat_write/timesteps(1)))
               E_write_MPI(iprt) = E_write_MPI(iprt) + dot_product(E,B)/norm2(B) / (deltat_write/timesteps(1))
            end if
         end if

         ! Write if time >= next_time_write
         !> Fields for
         if ((.NOT. average_coll_fields) .AND. (.NOT. needs_recycle(iprt)) &
             .AND. (sim%time - time_startsim + timesteps(1) .ge. next_time_write)) then
            call sim%fields%calc_EBpsiU(sim%time, prt(iprt)%i_elm, prt(iprt)%st, prt(iprt)%x(3), E, B, psi, U, E_par_norm=.TRUE.)
            call sim%fields%calc_NjTj(sim%time, prt(iprt)%i_elm, prt(iprt)%st, prt(iprt)%x(3), dat%m_i_over_m_imp, ne, the, ni, ti)
            n_write_MPI(iprt) = ne
            T_write_MPI(iprt) = the*K_BOLTZ/EL_CHG
            E_write_MPI(iprt) = dot_product(E,B)/norm2(B)  ! Parallel component
         end if

      end do
      !$omp end parallel do
   end select

   !> Recycle the markers flagged during the push.
   if (any(needs_recycle)) then
      do iprt = 1, nprt_MPI
         if (.NOT. needs_recycle(iprt)) cycle
         call recycle_lost(sim%groups(1)%particles, sim%fields, nprt_MPI, iprt, sim%time, &
                           sim%groups(1)%mass, p_cutoff_recycle, p_before_reweight, &
                           min_weight, n_write_MPI, T_write_MPI, E_write_MPI, ifail)
         if (ifail .ne. 0) then
            !> No donor available. 
            n_recycle_fail = n_recycle_fail + 1
         else
            needs_recycle(iprt) = .false.
         end if
      end do
   end if

   ! Gather losses
   n_lost = count(lost)
   if (n_lost .gt. 0) then   ! if a loss occured during the last timestep
      call append_filtered(t_losses, t_losses_tmp, lost)
      call append_filtered(p_losses, p_losses_tmp, lost)
      call append_filtered(xi_losses, xi_losses_tmp, lost)
      call append_filtered(R_losses, R_losses_tmp, lost)
      call append_filtered(Z_losses, Z_losses_tmp, lost)
      call append_filtered(Phi_losses, Phi_losses_tmp, lost)
      call append_filtered(Weight_losses, Weight_losses_tmp, lost)
      call append_filtered(iprt_losses, iprt_losses_tmp, lost)
      lost  = .FALSE.
   end if

   !> Bound the loss buffers, which grow between two scheduled writes and are all
   !> gathered on rank 0. 
   if (mod(istep, 20) .eq. 0) then
      n_loss_buffer = size(t_losses)
      call MPI_Allreduce(n_loss_buffer, n_loss_buffer_max, 1, MPI_INTEGER, MPI_MAX, &
                         MPI_COMM_WORLD, ierr)
      if (n_loss_buffer_max .gt. max_loss_per_rank) then
         if (rank .eq. 0) then
            write(*,*) 'INFO: loss buffer reached ', n_loss_buffer_max, &
                       ' records on the largest rank, flushing early to bound memory'
            flush(6)
         end if
         call write_loss_to_hdf5(fnout_losses, t_losses, p_losses, xi_losses, R_losses, &
                                 Z_losses, Phi_losses, Weight_losses, iprt_losses)
      end if
   end if

   ! Set the next writting time
   if ((sim%time - time_startsim) .ge. next_time_write) then
      if (istep .lt. nstep) then
         next_time_write = next_time_write + deltat_write
      end if
      if (rank .eq. 0) then
         write(*,*) "Progression :", NINT(100*(sim%time - time_startsim)/duration), " %"
      end if
   end if

   ! Update the sim%time and read fields if necessary
   sim%time = sim%time + timesteps(1)
   if (sim%time >= target_time) then
      call with(sim,events, at=target_time)
      target_time = next_event_at(sim, events)
      call MPI_Barrier(MPI_COMM_WORLD, ierr)

      !> A field reload replaces the node and element lists. If the grid changed
      !> between the two restart files, for instance through refinement, the
      !> element indices carried by the markers refer to the previous grid and
      !> are no longer valid. interp_PRZ performs no bound check, so such an
      !> index is an out-of-range memory access rather than a reported error.
      !> Stale markers are invalidated here and replaced by the recycling pass.
      n_elm_max = sim%fields%element_list%n_elements
      n_stale   = 0
      do iprt = 1, nprt_MPI
         if (sim%groups(1)%particles(iprt)%i_elm .gt. n_elm_max) then
            sim%groups(1)%particles(iprt)%i_elm = 0
            needs_recycle(iprt) = .TRUE.
            n_stale = n_stale + 1
         end if
      end do
      n_stale_total = n_stale_total + n_stale
      if (n_stale .gt. 0) then
         write(*,*) 'WARNING: the grid changed at the field reload. n_elements =', n_elm_max, &
                    ', stale marker indices invalidated =', n_stale, ' / rank =', rank
         flush(6)
      end if
   end if

   ! Russian roulette
   if (((sim%time - time_startsim) .ge. next_time_reweight) .AND. (istep .ne. nstep)) then
      call correlated_russian_roulette(sim%groups(1)%particles, sim%fields, nprt_MPI, sim%time, sim%groups(1)%mass, p_cutoff, p_far, &
                                       p_before_reweight, prt_age, age_split, k, l, &
                                       use_psi_res, use_p_far_sweep, psi_res, delta_psi_res, delta_p_far)
      next_time_reweight = next_time_reweight + deltat_reweight
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
   end if

   ! Last step output
   if (istep .eq. nstep) then
      !> Diag output
      call with(sim, diag)
      call write_simulation_hdf5(sim, 'part_restart.h5')
      call write_collparam_hdf5(fnout_collparam, n_write, T_write, E_write, n_write_MPI, T_write_MPI, E_write_MPI, nprt_MPI, sim%time)
      call write_loss_to_hdf5(fnout_losses, t_losses, p_losses, xi_losses, R_losses, Z_losses, Phi_losses, Weight_losses, iprt_losses)  ! losses diagnostic
   end if

end do
write(*,*) "Progression : 100 % // MPI process :", rank
! #####################################################################################
! ############################### END SIMULATION LOOP #################################
! #####################################################################################
! Finalisation MPI
call MPI_Barrier(MPI_COMM_WORLD, ierr)
if (rank .eq. 0) then
   write(*,*) ">>>>>>>>>>>>>> END SIMULATION LOOP <<<<<<<<<<<<<<"
end if
call cpu_time(t1)
t_wall1 = MPI_Wtime()
!> cpu_time returns the CPU time summed over the OpenMP threads, so the wall clock
!> time from MPI_Wtime is the relevant figure for the run duration
write(*,*) 'Wall time: ', t_wall1-t_wall0, ' s / rank = ', rank
write(*,*) 'CPU time (all threads): ', t1-t0, ' s / rank = ', rank
write(*,*) 'Nb jumps: ', n_jump, ' / rank = ', rank
write(*,*) 'Weight jumps: ', w_jump, ' / rank = ', rank
write(*,*) 'Nb invalid markers found at push entry: ', n_invalid, ' / rank = ', rank
write(*,*) 'Nb failed recycling attempts: ', n_recycle_fail, ' / rank = ', rank
write(*,*) 'Nb stale marker indices after field reloads: ', n_stale_total, ' / rank = ', rank

! Finalize the simulation
deallocate (n_write, T_write, E_write, n_write_MPI, T_write_MPI, E_write_MPI)
deallocate (t_losses_tmp, p_losses_tmp, xi_losses_tmp, R_losses_tmp)
deallocate (Z_losses_tmp, Phi_losses_tmp, Weight_losses_tmp, iprt_losses_tmp)
deallocate (lost, needs_recycle, p_before_reweight, prt_age, ni)
call ccoll_deallocate(dat)

call sim%finalize

contains

function weight_MJ(nx,x,st,time,i_elm,fields,x_min,x_max,&
  n_real_param,real_param,n_int_param,int_param) result(weight)
    use constants,  only: PI, SPEED_OF_LIGHT, K_BOLTZ, ATOMIC_MASS_UNIT
    use mod_fields, only: fields_base
    implicit none
    !> Inputs:
    integer,intent(in)                          :: nx,i_elm,n_real_param
    integer,intent(in)                          :: n_int_param
    integer,dimension(:),allocatable,intent(in) :: int_param
    real*8,intent(in)                           :: time
    real*8,dimension(nx),intent(in)             :: x,x_min,x_max
    real*8,dimension(2),intent(in)              :: st
    class(fields_base),intent(in)               :: fields
    real*8,dimension(:),allocatable,intent(in)  :: real_param
    !> Local variables:
    real*8                          :: gamma        ! relativistic Lorentz factor of the marker
    real*8                          :: n_e, T_e, ti ! electron density [m^-3], electron/ion temperatures [K]
    real*8                          :: psi, U       ! poloidal flux and electric potential at the marker
    real*8,dimension(3)             :: B_field, E_field
    real*8                          :: theta        ! normalised temperature kT_e/(m c^2)
    real*8                          :: log_norm     ! log of the theta-dependent normalisation theta*K_2(1/theta)
    real*8                          :: chi2         ! (psi-psi_axis)^2/(2 sigma_psi^2); Pi(psi) = exp(-chi2)
    real*8                          :: log_arg      ! full exponent of the (unnormalised) weight, including 1/Pi(psi)
    real*8,dimension(:),allocatable :: ni           ! ion density per species [m^-3]
    !> Output:
    real*8 :: weight

    ! Asymptotic-series coefficients of K_2(z) for z = 1/theta -> +infinity (DLMF 10.40.2):
    ! K_2(1/theta) ~ sqrt(pi*theta/2) exp(-1/theta) [1 + c1*theta + c2*theta^2 + O(theta^3)]
    ! The leading term alone underestimates K_2 by ~10% near 30 keV; two terms reach machine precision.
    real*8,parameter :: c1 = 15.0d0/8.0d0
    real*8,parameter :: c2 = 105.0d0/128.0d0
    ! Largest exponent kept before the marker is treated as lying outside the sampling
    ! envelope Pi(psi). Set well below log(HUGE) ~ 709 to leave headroom for the n_e*p^2
    ! prefactor; it is only reached ~35 sigma from the axis, i.e. never for sampled markers.
    real*8,parameter :: MAX_LOG_ARG = 600.0d0

    if (i_elm .gt. 0) then
      allocate(ni(ubound(dat%Ii, dim=1)))
      !> Local field and profile evaluation at the marker position
      call fields%calc_NjTj(time, i_elm, st, x(3), dat%m_i_over_m_imp, n_e, T_e, ni, ti)
      !> Only psi is needed here, so the E = eta*J_par normalisation is skipped.
      !> This function is evaluated many times per accepted marker by the rejection
      !> sampler, so the saving is significant at initialisation.
      call fields%calc_EBpsiU(time, i_elm, st, x(3), E_field, B_field, psi, U, E_par_norm=.FALSE.)

      !> Lorentz factor
      gamma = sqrt(1.0d0 + (x(4)/(real_param(1)*SPEED_OF_LIGHT))**2)

      !> Normalised electron temperature theta = k T_e/(m c^2), with T_e in [K].
      theta = T_e*K_BOLTZ / (real_param(1)*ATOMIC_MASS_UNIT*SPEED_OF_LIGHT**2)

      !> Log of the temperature-dependent normalisation theta*K_2(1/theta)
      log_norm = log(theta) + 0.5d0*log(PI*theta/2.0d0) &
               + log(1.0d0 + c1*theta + c2*theta*theta)

      !> Squared poloidal-flux distance to the axis
      chi2 = (psi - psi_axis)**2 / (2.0d0*sigma_psi**2)

      !> Maxwell-Juttner weight (unnormalised), reweighted by 1/Pi(psi)
      log_arg = (1.0d0 - gamma)/theta - log_norm + chi2
      if (log_arg >= MAX_LOG_ARG) then
         ! Marker effectively outside the sampling envelope: its sampling probability is ~0,
         ! so it is a numerical artifact and gets zero weight. Use huge(1.0d0) instead if a
         ! finite cap better matches your reweighting scheme.
         weight = 0.0d0
      else
         weight = n_e * x(4)*x(4) * exp(log_arg)
      end if
      deallocate(ni)
    else
      weight = 0.d0
    end if
end function weight_MJ

function pdf_uniform(nx,x,st,time,i_elm,fields,x_min,x_max,&
  n_real_param,real_param,n_int_param,int_param)
    use constants,  only: PI
    use mod_fields, only: fields_base
    implicit none
    !> Inputs:
    integer,intent(in)                          :: nx,i_elm,n_real_param
    integer,intent(in)                          :: n_int_param
    integer,dimension(:),allocatable,intent(in) :: int_param
    real*8,intent(in)                           :: time
    real*8,dimension(nx),intent(in)             :: x,x_min,x_max
    real*8,dimension(2),intent(in)              :: st
    class(fields_base),intent(in)               :: fields
    real*8,dimension(:),allocatable,intent(in)  :: real_param
    !> Variables
    real*8                                      :: n_e, T_e, ti
    real*8,dimension(:),allocatable             :: ni
    real*8                                      :: psi,U
    real*8,dimension(3)                         :: B_field,E_field
    real*8                                      :: Pi_function, weight_test
    !> Outputs:
    real*8 :: pdf_uniform

    if (i_elm .gt. 0) then
      !> Evalutate pdf. Only psi is needed, so the E = eta*J_par normalisation is skipped.
      call fields%calc_EBpsiU(time,i_elm,st,x(3),E_field,B_field,psi,U,E_par_norm=.FALSE.)
      Pi_function = exp(-(psi-psi_axis)**2 / (2*sigma_psi**2))

      allocate(ni(ubound(dat%Ii, dim=1)))
      call fields%calc_NjTj(time, i_elm, st, x(3), dat%m_i_over_m_imp, n_e, T_e, ni, ti)
      deallocate(ni)

      weight_test =  weight_MJ(nx,x,st,time,i_elm,fields,x_min,x_max,&
      n_real_param,real_param,n_int_param,int_param)
      if (weight_test .lt. (n_e/range_weight)) then
         pdf_uniform = 0.0d0
      else
         pdf_uniform = 3.d0/((x_max(1)**2-x_min(1)**2)*(x_max(2)-x_min(2))* &
            (x_max(3)-x_min(3))*(x_max(4)**3-x_min(4)**3)* &
            (cos(x_min(5))-cos(x_max(5)))*PI) * Pi_function
      end if
   else
      pdf_uniform = 0.0d0
   end if
end function pdf_uniform

function sup_pdf_uniform(nx,x_min,x_max,n_real_param,real_param,&
n_int_param,int_param) result(sup_pdf)
   use constants,  only: PI
   use mod_fields, only: fields_base
   implicit none
   !> Inputs:
   integer,intent(in)                          :: nx,n_real_param
   integer,intent(in)                          :: n_int_param
   integer,dimension(:),allocatable,intent(in) :: int_param
   real*8,dimension(nx),intent(in)             :: x_min,x_max
   real*8,dimension(:),allocatable,intent(in)  :: real_param
   !> Outputs:
   real*8 :: sup_pdf
   !> Evalutate the upper extremum of the pdf
   sup_pdf = 3.d0/((x_max(1)**2-x_min(1)**2)*(x_max(2)-x_min(2))* &
   (x_max(3)-x_min(3))*(x_max(4)**3-x_min(4)**3)* &
   (cos(x_min(5))-cos(x_max(5)))*PI)
end function sup_pdf_uniform

function sup_gdf_uniform_phase(nx,x_min,x_max,n_real_param,real_param,&
n_int_param,int_param) result(sup_gdf)
   use constants,  only: PI
   use mod_fields, only: fields_base
   implicit none
   !> Inputs:
   integer,intent(in)                          :: nx,n_real_param
   integer,intent(in)                          :: n_int_param
   integer,dimension(:),allocatable,intent(in) :: int_param
   real*8,dimension(nx),intent(in)             :: x_min,x_max
   real*8,dimension(:),allocatable,intent(in)  :: real_param
   !> Outputs:
   real*8 :: sup_gdf
   !> Evalutate the upper extremum of the pdf
   sup_gdf = 3.d0/((x_max(1)**2-x_min(1)**2)*(x_max(2)-x_min(2))* &
   (x_max(3)-x_min(3))*(x_max(4)**3-x_min(4)**3)* &
   (cos(x_min(5))-cos(x_max(5)))*PI)
end function sup_gdf_uniform_phase

subroutine gdf_uniform_sampler(nx,x,st,time,i_elm,fields,&
x_min,x_max,n_real_param,real_param,n_int_param,int_param)
   use mod_fields, only: fields_base
   implicit none
   !> Inputs:
   integer,intent(in)                          :: nx,n_real_param
   integer,intent(in)                          :: n_int_param
   integer,dimension(:),allocatable,intent(in) :: int_param
   real*8,intent(in)                           :: time
   real*8,dimension(nx),intent(in)             :: x_min,x_max
   class(fields_base),intent(in)               :: fields
   real*8,dimension(:),allocatable,intent(in)  :: real_param
   !> Inputs-Outputs:
   integer,intent(inout)                       :: i_elm
   real*8,dimension(2),intent(inout)           :: st
   real*8,dimension(nx),intent(inout)          :: x
   !> Variables:
   integer :: ifail
   real*8  :: R_found, Z_found, s_found, t_found !< separate outputs, see below
   !> Compute new particle position in phase space
   x(1) = sqrt(x_min(1)**2 + (x_max(1)**2 - x_min(1)**2)*x(1))
   x(2:3) = x_min(2:3) + (x_max(2:3)-x_min(2:3))*x(2:3)
   x(4) = x_min(4) + (x_max(4)-x_min(4))*x(4)
   x(5) = acos(cos(x_min(5))+(cos(x_max(5))-cos(x_min(5)))*x(5))
   x(6) = x_min(6) + (x_max(6)-x_min(6))*x(6)
   !> find RZ coordinates. x(1) and x(2) must not serve as both the input and the
   !> output of find_RZ: F2018 15.5.2.13 forbids it, and find_RZ_general writes
   !> its outputs inside the loop over candidate elements. Separate locals are
   !> used, primed with the sampled position because find_RZ_general leaves them
   !> unassigned on failure, so a rejected draw keeps the point it was drawn at.
   R_found = x(1)
   Z_found = x(2)
   s_found = st(1)
   t_found = st(2)
   call find_RZ(fields%node_list,fields%element_list,x(1),x(2),&
   R_found,Z_found,i_elm,s_found,t_found,ifail)
   x(1)  = R_found
   x(2)  = Z_found
   st(1) = s_found
   st(2) = t_found
end subroutine gdf_uniform_sampler

function gdf_uniform_phase(nx,x,st,time,i_elm,fields,x_min,x_max,&
n_real_param,real_param,n_int_param,int_param) result(gdf)
   use constants,  only: PI
   use mod_fields, only: fields_base
   implicit none
   !> Inputs:
   integer,intent(in)                          :: nx,i_elm,n_real_param
   integer,intent(in)                          :: n_int_param
   integer,dimension(:),allocatable,intent(in) :: int_param
   real*8,intent(in)                           :: time
   real*8,dimension(nx),intent(in)             :: x,x_min,x_max
   real*8,dimension(2),intent(in)              :: st
   class(fields_base),intent(in)               :: fields
   real*8,dimension(:),allocatable,intent(in)  :: real_param
   !> Outputs:
   real*8 :: gdf
   !> Evalutate pdf
   gdf = 3.d0/((x_max(1)**2-x_min(1)**2)*(x_max(2)-x_min(2))* &
   (x_max(3)-x_min(3))*(x_max(4)**3-x_min(4)**3)* &
   (cos(x_min(5))-cos(x_max(5)))*PI)
end function gdf_uniform_phase

subroutine spherical_p_cartesian_q_to_relativistic_gc(p_inout,&
n_x,x,time,fields,n_real_param,real_param,n_int_param,int_param)
   use mod_particle_types,        only: particle_base
   use mod_particle_types,        only: particle_gc_relativistic
   use mod_fields,                only: fields_base
   implicit none
   !> Inputs-Outputs:
   class(particle_base),intent(inout) :: p_inout
   !> Inputs:
   class(fields_base),intent(in)               :: fields
   integer,intent(in)                          :: n_x,n_real_param,n_int_param
   integer,dimension(:),allocatable,intent(in) :: int_param
   real*8,intent(in)                           :: time
   real*8,dimension(n_x),intent(in)            :: x
   real*8,dimension(:),allocatable,intent(in)  :: real_param
   !> variables
   real*8              :: psi,U, normB
   real*8,dimension(3) :: B_field,E_field

   select type (p=>p_inout)
   type is (particle_gc_relativistic)
   !> Only B is needed here, so the E = eta*J_par normalisation is skipped
   call fields%calc_EBpsiU(time,p%i_elm,p%st,p%x(3),&
   E_field,B_field,psi,U,E_par_norm=.FALSE.)
   normB = norm2(B_field)
   p%p = x(4)*[cos(x(5)),&
   (x(4)*((sin(x(5)))**2))/(2d0*real_param(1)*normB)]
   p%q = int(x(6),kind=1)
   end select
end subroutine spherical_p_cartesian_q_to_relativistic_gc

end program hot_tail_track