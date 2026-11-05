!> Shared input/output routines for the grad(chi) diagnostics jorek2_gradchi_domm and
!! jorek2_gradchi_full. Both evaluate the gradient of the vacuum scalar potential chi
!! (B_vac = grad(chi)) at arbitrary Cartesian points.
!!
!! Command line (all optional):
!!   --xyz <file>      input points            (default: xyz.dat)
!!   --out <file>      output file             (default: fields_xyz.dat)
!!   --restart <name>  restart file, no .h5    (default: jorek_restart, jorek2_gradchi_full only)
!!
!! Input file: an optional first line with the number of points, then one point per line
!!   x y z
!! in metres. Empty lines and lines starting with '#' are skipped.
!!
!! Output file: one header line starting with '#', then one line per input point, in input order
!!   x y z gradchi_x gradchi_y gradchi_z [extra columns]
!! The gradient is in Tesla. Coordinates follow the JOREK convention x = R cos(phi),
!! y = -R sin(phi), z = Z.
module mod_gradchi_xyz

  use, intrinsic :: iso_fortran_env, only: output_unit
  use mpi_mod

  implicit none

  private
  public :: parse_gradchi_args, read_xyz_points, gradchi_cartesian, gradchi_abort
  public :: gradchi_header, gradchi_fmt

  !> Header and number format of the output file (15 significant digits)
  character(len=*), parameter :: gradchi_header = '#  x(m)  y(m)  z(m)  gradchi_x(T)  gradchi_y(T)  gradchi_z(T)'
  character(len=*), parameter :: gradchi_fmt    = '(6ES24.15)'

contains

  !> Read --xyz, --out and (if allowed) --restart from the command line
  subroutine parse_gradchi_args(my_id, prog, xyz_path, out_path, restart_name)
    integer,                    intent(in)    :: my_id
    character(len=*),           intent(in)    :: prog
    character(len=*),           intent(inout) :: xyz_path, out_path
    character(len=*), optional, intent(inout) :: restart_name

    character(len=1024) :: arg
    integer             :: narg, iarg, ios

    narg = command_argument_count()
    iarg = 1
    do while (iarg <= narg)
      call get_command_argument(iarg, arg)
      select case (trim(arg))
      case ('--xyz', '--out', '--restart')
        if (iarg >= narg) call usage(my_id, prog, 'Missing value after '//trim(arg), present(restart_name))
        if (trim(arg) == '--restart' .and. .not. present(restart_name)) &
          call usage(my_id, prog, 'Unknown argument: --restart', .false.)
        select case (trim(arg))
        case ('--xyz');     call get_command_argument(iarg+1, xyz_path)
        case ('--out');     call get_command_argument(iarg+1, out_path)
        case ('--restart'); call get_command_argument(iarg+1, restart_name)
        end select
        iarg = iarg + 2
      case ('--help', '-h')
        call usage(my_id, prog, '', present(restart_name))
      case default
        call usage(my_id, prog, 'Unknown argument: '//trim(arg), present(restart_name))
      end select
    end do

    ! Empty the output file now, so that a run that stops early never leaves the
    ! results of a previous run behind
    if (my_id == 0) then
      open(87, file=trim(out_path), action='write', status='replace', iostat=ios)
      if (ios /= 0) call gradchi_abort('Could not open output file '//trim(out_path))
      close(87)
    end if
  end subroutine parse_gradchi_args

  subroutine usage(my_id, prog, message, with_restart)
    integer,          intent(in) :: my_id
    character(len=*), intent(in) :: prog, message
    logical,          intent(in) :: with_restart

    if (my_id == 0) then
      if (len_trim(message) > 0) write(*,*) 'ERROR: ', trim(message)
      if (with_restart) then
        write(*,*) 'Usage: ', prog, ' [--xyz file] [--out file] [--restart name] < namelist'
      else
        write(*,*) 'Usage: ', prog, ' [--xyz file] [--out file] < namelist'
      end if
    end if
    call gradchi_abort('')
  end subroutine usage

  !> Read the evaluation points (call on one MPI task only)
  subroutine read_xyz_points(path, np, xp, yp, zp)
    character(len=*),    intent(in)  :: path
    integer,             intent(out) :: np
    real*8, allocatable, intent(out) :: xp(:), yp(:), zp(:)

    integer, parameter  :: unit = 26
    character(len=2048) :: line
    real*8              :: x, y, z
    integer             :: ios, n_header, ip, pass
    logical             :: first

    open(unit, file=trim(path), action='read', status='old', iostat=ios)
    if (ios /= 0) call gradchi_abort('Could not open input file '//trim(path))

    ! Pass 1 counts the points, pass 2 stores them
    do pass = 1, 2
      rewind(unit)
      first = .true.
      ip    = 0
      do
        read(unit, '(A)', iostat=ios) line
        if (ios /= 0) exit
        line = adjustl(line)
        if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
        read(line, *, iostat=ios) x, y, z
        if (ios /= 0) then
          ! A single integer on the first data line is the optional point count
          read(line, *, iostat=ios) n_header
          if (first .and. ios == 0) then
            first = .false.
            cycle
          end if
          call gradchi_abort('Could not parse point in '//trim(path)//': '//trim(line))
        end if
        first = .false.
        ip = ip + 1
        if (pass == 2) then
          xp(ip) = x; yp(ip) = y; zp(ip) = z
        end if
      end do
      if (pass == 1) then
        np = ip
        if (np == 0) call gradchi_abort('No points found in '//trim(path))
        allocate(xp(np), yp(np), zp(np))
      end if
    end do
    close(unit)
  end subroutine read_xyz_points

  !> Convert the chi derivatives returned by mod_chi (d/dR, d/dZ, d/dphi) into a
  !! Cartesian gradient at the point (x, y)
  subroutine gradchi_cartesian(chi, x, y, gx, gy, gz)
    real*8, intent(in)  :: chi(0:,0:,0:)
    real*8, intent(in)  :: x, y
    real*8, intent(out) :: gx, gy, gz

    real*8 :: R, g_R, g_Z, g_phi

    R     = sqrt(x**2 + y**2)
    g_R   = chi(1,0,0)
    g_Z   = chi(0,1,0)
    g_phi = chi(0,0,1) / R
    ! With x = R cos(phi), y = -R sin(phi): e_R = (x, y)/R, e_phi = (y, -x)/R
    gx = g_R * x/R + g_phi * y/R
    gy = g_R * y/R - g_phi * x/R
    gz = g_Z
  end subroutine gradchi_cartesian

  !> Print an error and stop all MPI tasks with a non-zero exit code
  subroutine gradchi_abort(message)
    character(len=*), intent(in) :: message
    integer :: ierr

    if (len_trim(message) > 0) write(*,*) 'ERROR: ', trim(message)
    flush(output_unit)
    call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  end subroutine gradchi_abort

end module mod_gradchi_xyz
