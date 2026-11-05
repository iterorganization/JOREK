!> Evaluates grad(chi), the vacuum magnetic field, from the Dommaschk potentials alone
!! at arbitrary Cartesian points. No restart file is needed, only the namelist (stdin)
!! with F0 and domm_file, so the tool is cheap enough to be called from other codes.
!! CARIDDI uses it this way to get grad(chi) on its control surface for the
!! stellarator free-boundary coupling.
!!
!! For stellarator models it requires a build with USE_DOMM=1 and domm_file set,
!! otherwise the Dommaschk sum is not part of chi. Use jorek2_gradchi_full to include
!! the finite-element correction of chi. For tokamak models it returns F0*grad(phi).
!!
!! Usage:  jorek2_gradchi_domm [--xyz xyz.dat] [--out fields_xyz.dat] < namelist
!! File formats are described in mod_gradchi_xyz.
program jorek2_gradchi_domm

  use phys_module,     only: F0, R_domm
  use mod_parameters,  only: n_order
  use mod_chi,         only: init_chi_basis, get_chi_domm
  use mod_gradchi_xyz
  use mpi_mod

  implicit none

  integer             :: my_id, ierr, required, provided, np, i
  real*8              :: R, phi, gx, gy, gz
  real*8              :: chi(0:n_order-1,0:n_order-1,0:n_order-1)
  real*8, allocatable :: xp(:), yp(:), zp(:)
  character(len=1024) :: xyz_path, out_path

#ifdef FUNNELED
  required = MPI_THREAD_FUNNELED
#else
  required = MPI_THREAD_MULTIPLE
#endif
  call MPI_Init_thread(required, provided, ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, my_id, ierr)

  xyz_path = 'xyz.dat'
  out_path = 'fields_xyz.dat'
  call parse_gradchi_args(my_id, 'jorek2_gradchi_domm', xyz_path, out_path)

  ! --- Read the namelist and the Dommaschk coefficients
  call det_modes()
  call initialise_and_broadcast_parameters(my_id, "__NO_FILENAME__", .false.)
  call init_chi_basis()

#if STELLARATOR_MODEL && !defined(USE_DOMM)
  call gradchi_abort('jorek2_gradchi_domm needs a build with USE_DOMM=1, use jorek2_gradchi_full instead')
#endif

  if (my_id == 0) then
    write(*,*) '--- jorek2_gradchi_domm ---'
    write(*,*) 'F0     = ', F0
    write(*,*) 'R_domm = ', R_domm
    write(*,*) 'input  = ', trim(xyz_path)
    write(*,*) 'output = ', trim(out_path)

    call read_xyz_points(xyz_path, np, xp, yp, zp)
    write(*,*) 'Number of points: ', np

    open(87, file=trim(out_path), action='write', status='replace', iostat=ierr)
    if (ierr /= 0) call gradchi_abort('Could not open output file '//trim(out_path))
    write(87,'(a)') gradchi_header
    do i = 1, np
      R   = sqrt(xp(i)**2 + yp(i)**2)
      phi = atan2(-yp(i), xp(i))
      chi = get_chi_domm(R, zp(i), phi)
      call gradchi_cartesian(chi, xp(i), yp(i), gx, gy, gz)
      write(87,gradchi_fmt) xp(i), yp(i), zp(i), gx, gy, gz
    end do
    close(87)

    write(*,*) 'Output written to ', trim(out_path)
  end if

  call MPI_Finalize(ierr)

end program jorek2_gradchi_domm
