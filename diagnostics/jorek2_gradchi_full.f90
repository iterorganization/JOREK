!> Evaluates the full grad(chi), the vacuum magnetic field, at arbitrary Cartesian points.
!! It reads a restart file, locates each point in the mesh and calls get_chi, which adds
!! the finite-element correction of chi (or the external field) to the Dommaschk part,
!! depending on the build (USE_DOMM, USE_EXT_FIELD).
!!
!! Points outside the mesh get i_elm = -1. Their value is the Dommaschk part alone in a
!! USE_DOMM build, where this is exact, and NaN otherwise.
!!
!! Usage:  jorek2_gradchi_full [--xyz xyz.dat] [--out fields_xyz.dat]
!!                             [--restart jorek_restart] < namelist
!! File formats are described in mod_gradchi_xyz. Three extra columns give the element
!! index and the local coordinates s, t of each point.
program jorek2_gradchi_full

#ifndef USE_DOMM
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
#endif
  use data_structure
  use phys_module
  use mod_parameters,       only: n_order
  use mod_chi,              only: init_chi_basis, get_chi, get_chi_domm
  use nodes_elements
  use mod_boundary,         only: boundary_from_grid
  use mod_import_restart
  use mod_element_rtree,    only: populate_element_rtree
  use basis_at_gaussian,    only: initialise_basis
  use mod_gradchi_xyz
  use mpi_mod

  implicit none

  integer             :: my_id, ierr, ierr_bcast, required, provided, np, i
  integer             :: i_elm, ifail, checked_elms, n_outside
  real*8              :: R, phi, R_out, Z_out, s, t, gx, gy, gz
  real*8              :: chi(0:n_order-1,0:n_order-1,0:n_order-1)
  real*8, allocatable :: xp(:), yp(:), zp(:)
  logical             :: found
  character(len=1024) :: xyz_path, out_path, restart_name

#ifdef FUNNELED
  required = MPI_THREAD_FUNNELED
#else
  required = MPI_THREAD_MULTIPLE
#endif
  call MPI_Init_thread(required, provided, ierr)
  call init_threads()
  call MPI_Comm_rank(MPI_COMM_WORLD, my_id, ierr)

  xyz_path     = 'xyz.dat'
  out_path     = 'fields_xyz.dat'
  restart_name = 'jorek_restart'
  call parse_gradchi_args(my_id, 'jorek2_gradchi_full', xyz_path, out_path, restart_name)

  ! --- Read the namelist and the restart file
  call det_modes()
  call initialise_and_broadcast_parameters(my_id, "__NO_FILENAME__", .false.)
  call initialise_basis()

  if (my_id == 0) then
    ! import_restart crashes on a missing file, so check first
    if (rst_hdf5 == 0) then
      inquire(file=trim(restart_name)//'.rst', exist=found)
    else
      inquire(file=trim(restart_name)//'.h5', exist=found)
    end if
    if (.not. found) call gradchi_abort('Restart file not found: '//trim(restart_name))
    call import_restart(node_list, element_list, trim(restart_name), rst_format, ierr)
  end if
  call MPI_Bcast(ierr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr_bcast)
  if (ierr /= 0) call gradchi_abort('Could not read restart file '//trim(restart_name))

  call init_chi_basis()
  call broadcast_phys(my_id)
  call broadcast_elements(my_id, element_list)
  call broadcast_nodes(my_id, node_list)
  if (my_id == 0) call boundary_from_grid(node_list, element_list, bnd_node_list, bnd_elm_list, .false.)
  call broadcast_boundary(my_id, bnd_elm_list, bnd_node_list)
  call populate_element_rtree(node_list, element_list)

  if (my_id == 0) then
    write(*,*) '--- jorek2_gradchi_full ---'
    write(*,*) 'F0      = ', F0
    write(*,*) 'restart = ', trim(restart_name)
    write(*,*) 'input   = ', trim(xyz_path)
    write(*,*) 'output  = ', trim(out_path)

    call read_xyz_points(xyz_path, np, xp, yp, zp)
    write(*,*) 'Number of points: ', np

    open(87, file=trim(out_path), action='write', status='replace', iostat=ierr)
    if (ierr /= 0) call gradchi_abort('Could not open output file '//trim(out_path))
    write(87,'(a)') gradchi_header//'  i_elm  s  t'

    n_outside = 0
    do i = 1, np
      R   = sqrt(xp(i)**2 + yp(i)**2)
      phi = atan2(-yp(i), xp(i))

      i_elm = 0
      call find_RZP(node_list, element_list, R, zp(i), phi, R_out, Z_out, i_elm, s, t, ifail, checked_elms)

      if (ifail == 0 .and. i_elm > 0) then
        chi = get_chi(R, zp(i), phi, node_list, element_list, i_elm, s, t)
        call gradchi_cartesian(chi, xp(i), yp(i), gx, gy, gz)
      else
        n_outside = n_outside + 1
        i_elm = -1; s = 0.d0; t = 0.d0
#ifdef USE_DOMM
        chi = get_chi_domm(R, zp(i), phi)
        call gradchi_cartesian(chi, xp(i), yp(i), gx, gy, gz)
#else
        gx = ieee_value(gx, ieee_quiet_nan); gy = gx; gz = gx
#endif
      end if

      write(87,'(6ES24.15,I10,2ES24.15)') xp(i), yp(i), zp(i), gx, gy, gz, i_elm, s, t
    end do
    close(87)

    write(*,*) 'Output written to ', trim(out_path)
    if (n_outside > 0) write(*,*) 'WARNING: ', n_outside, ' points are outside the mesh (i_elm = -1)'
  end if

  call MPI_Finalize(ierr)

end program jorek2_gradchi_full
