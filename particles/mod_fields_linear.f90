!> Module for linearly interpolating in time between values and deltas
!> in JOREK restart files. Contains an action to read the fields.
module mod_fields_linear
use data_structure
use mod_particle_sim
use mod_event
use mod_fields
use mod_interp
implicit none
private
public jorek_fields_interp_linear, read_jorek_fields_interp_linear, last_file_before_time

!> Action to read in the fields into sim%fields
type, extends(action) :: read_jorek_fields_interp_linear
  character(len=80) :: basename = 'jorek' !< Comes before the file number or extension
  integer :: i = 0 !< Number of the restart file to read. Set to -1 to not include. Corresponds to the index of NOW
  integer :: rst_format = 0 !< Format of restart file if .rst type
  logical :: stop_at_end = .true. !< Whether to stop the simulation at the end of the file list
  contains
    procedure :: do => do_read
end type read_jorek_fields_interp_linear
 
interface read_jorek_fields_interp_linear
  module procedure new_read_jorek_fields_interp_linear
end interface read_jorek_fields_interp_linear

!> store enough data for linear interpolation
!> in time.
!> This does not really work for changing grids in time!
!> Use at your own peril in that case. (e.g. i_elm and s,t for a specific spatial position might depend on time)
!>
!> The reason behind not using deltas is that we do not have to alter much code
!> and can import two restarts which are not consecutive and still interpolate.
type, extends(fields_base) :: jorek_fields_interp_linear
  real*8 :: time_now  = 0.d0 !< Time of current restart file (SI units)
  real*8 :: time_prev = 0.d0!< Time of previous restart file (SI units)
  contains
    procedure :: interp_PRZ => do_interp_PRZ_1
    procedure :: interp_PRZ_2 => do_interp_PRZ_2
end type jorek_fields_interp_linear
contains

!> Interpolate a variable at a specific position (with phi), with first derivatives only
pure subroutine do_interp_PRZ_1(this, time, i_elm, i_v, n_v, s, t, phi, P, P_s, P_t, P_phi, P_time, R, R_s, R_t, Z, Z_s, Z_t)
  use mod_interp
  use constants, only: mu_zero, mass_proton
  use phys_module, only: tstep, central_mass, central_density
  use mod_linear, only: linear_interp_differentials
  use mod_linear, only: linear_interp_differentials_dt
  class(jorek_fields_interp_linear),  intent(in)  :: this
  real*8,                   intent(in)  :: time !< Time at which to calculate this variable
  integer,                  intent(in)  :: i_elm
  integer,                  intent(in)  :: n_v, i_v(n_v)
  real*8,                   intent(in)  :: s, t, phi
  real*8,                   intent(out) :: P(n_v), P_s(n_v), P_t(n_v), P_phi(n_v), P_time(n_v)
  real*8,                   intent(out) :: R, R_s, R_t, Z, Z_s, Z_t

  real*8                 :: df, dt
  real*8, dimension(n_v) :: Pd, Pd_s, Pd_t, Pd_phi
  real*8                 :: t_jorek
  
  ! JOREK time step in seconds
  t_jorek = tstep*sqrt(mu_zero * mass_proton * central_mass * central_density * 1.d20)

  P_time = 0.d0
  
  !> interpolate values
  call interp_PRZ(this%node_list,this%element_list,i_elm,i_v,n_v,s,t,phi,P, P_s, P_t, P_phi, R,R_s,R_t,Z,Z_s,Z_t)

  !> interpolate differentials
  if(t_jorek .gt. 0.d0) then
    call interp_PRZ(this%node_list,this%element_list,i_elm,i_v,n_v,s,t,phi, &
      Pd,Pd_s,Pd_t,Pd_phi,R,R_s,R_t,Z,Z_s,Z_t,deltas=.true.)
    if(abs(this%time_now-this%time_prev) .gt. 1d-10 .and. .not. this%static) then
      !> compute time fraction df
      dt = 1.d0/(this%time_now - this%time_prev)
      df = (this%time_now - time)*dt
      !> apply linear interpolation
      P     = linear_interp_differentials(n_v,P,Pd,df)
      P_s   = linear_interp_differentials(n_v,P_s,Pd_s,df)
      P_t   = linear_interp_differentials(n_v,P_t,Pd_t,df)
      P_phi = linear_interp_differentials(n_v,P_phi,Pd_phi,df)
    else
      dt = 1.d0/t_jorek
    endif
    !> compute time derivative
    P_time = linear_interp_differentials_dt(n_v,Pd,dt) 
  endif

end subroutine do_interp_PRZ_1

!> This procedure interpolates a variable, its first and second order derivatives in space
!> and first derivatives in time. Only first order spatial derivatives in s and t
!> are also derived in time.
!> inputs:
!>   this:    (jorek_fields_interp_linear) interpolation class
!>   time:    (real8) interpolation time
!>   i_elm:   (integer) mesh element index
!>   n_v:     (integer) number of jorek fields to be interpolated
!>   i_v:     (integer)(n_v) indices of the jorek fields to be interpolated
!>   s,t,phi: (real8) interpolation position in local (s,t) coordinates and
!>            toroidal angle (phi)
!> outputs:
!>   P:                (real8)(n_v) interpolated values
!>   P_s,P_t,P_time:   (real8)(n_v) interpolated first order derivatives
!>   P_ss,P_st,P_tt:   (real8)(n_v) interpolated second order derivatives
!>   P_sphi,P_tphi:    (real8)(n_v) interpolated second order derivatives
!>   P_stime,P_ttime:  (real8)(n_v) interpolated spatial-time cross derivatives
!>   R,Z:              (real8) global coordinates
!>   R_s,R_t,Z_s,Z_t:: (real8) global coordinates first order derivatives
!>   R_ss,R_st,R_tt:   (real8) R second order derivatives
!>   Z_ss,Z_st,Z_tt:   (real8) Z second order derivatives
pure subroutine do_interp_PRZ_2(this,time,i_elm,i_v,n_v,s,t,phi, &
  P,P_s,P_t,P_phi,P_time,P_ss,P_st,P_tt,P_sphi,P_tphi,P_stime,   &
  P_ttime,R,R_s,R_t,R_ss,R_st,R_tt,Z,Z_s,Z_t,Z_ss,Z_st,Z_tt)
  !> load module and functions
  use mod_interp
  use constants, only: mu_zero,mass_proton
  use phys_module, only: tstep,central_mass,central_density
  use mod_linear, only: linear_interp_differentials
  use mod_linear, only: linear_interp_differentials_dt
  implicit none
  !> declare input variables
  class(jorek_fields_interp_linear), intent(in) :: this
  real(kind=8), intent(in)                      :: s, t, phi, time
  integer, intent(in)                           :: i_elm, n_v
  integer, dimension(n_v), intent(in)           :: i_v
  !> declare output variables
  real(kind=8), intent(out) :: R, R_s, R_t, R_ss, R_st, R_tt
  real(kind=8), intent(out) :: Z, Z_s, Z_t, Z_ss, Z_st, Z_tt
  real(kind=8), dimension(n_v), intent(out) :: P, P_s, P_t, P_phi, P_time
  real(kind=8), dimension(n_v), intent(out) :: P_ss, P_st, P_tt, P_sphi
  real(kind=8), dimension(n_v), intent(out) :: P_tphi, P_stime, P_ttime
  !> declare internal variables
  real(kind=8) :: t_jorek
  real(kind=8) :: time_fraction, inverse_time_interval
  !> additional values
  real(kind=8), dimension(n_v) :: P_phiphi
  !> additional differentials
  real(kind=8), dimension(n_v) :: dP, dP_s, dP_t, dP_phi, dP_st, dP_ss, dP_tt
  real(kind=8), dimension(n_v) :: dP_sphi, dP_tphi, dP_phiphi
  
  ! JOREK time step in seconds
  t_jorek = tstep*sqrt(mu_zero*central_density*mass_proton*central_mass*1.d20)

  P_time = 0.d0
  P_stime = 0.d0
  P_ttime = 0.d0
  
  !> interpolate variables
  call interp_PRZ(this%node_list,this%element_list,i_elm,i_v,n_v,  &
    s,t,phi,P,P_s,P_t,P_phi,P_st,P_ss,P_tt,P_sphi,P_tphi,P_phiphi, &
    R,R_s,R_t,R_st,R_ss,R_tt,Z,Z_s,Z_t,Z_st,Z_ss,Z_tt,&
    deltas=.false.)

  !> interpolate differentials    
  if(t_jorek .gt. 0.d0) then
    call interp_PRZ(this%node_list,this%element_list,i_elm,i_v,n_v, &
      s,t,phi,dP,dP_s,dP_t,dP_phi,dP_st,dP_ss,dP_tt,dP_sphi,        &
      dP_tphi,dP_phiphi,R,R_s,R_t,R_st,R_ss,R_tt,                   &
      Z,Z_s,Z_t,Z_st,Z_ss,Z_tt,deltas=.true.)
    if((this%time_now-this%time_prev) .gt. 1.d-10 .and. .not. this%static) then
       inverse_time_interval = 1.d0/(this%time_now-this%time_prev)
       !> compute time fraction
       time_fraction = inverse_time_interval*(this%time_now - time)
       !> apply linear interpolation
       P      = linear_interp_differentials(n_v,P,dP,time_fraction)
       P_s    = linear_interp_differentials(n_v,P_s,dP_s,time_fraction)
       P_t    = linear_interp_differentials(n_v,P_t,dP_t,time_fraction)
       P_phi  = linear_interp_differentials(n_v,P_phi,dP_phi,time_fraction)
       P_st   = linear_interp_differentials(n_v,P_st,dP_st,time_fraction)
       P_ss   = linear_interp_differentials(n_v,P_ss,dP_ss,time_fraction)
       P_tt   = linear_interp_differentials(n_v,P_tt,dP_tt,time_fraction)
       P_sphi = linear_interp_differentials(n_v,P_sphi,dP_sphi,time_fraction)
       P_tphi = linear_interp_differentials(n_v,P_tphi,dP_tphi,time_fraction)
    else
       inverse_time_interval = 1.d0/t_jorek
    endif
    !> compute time derivatives
    P_time  = linear_interp_differentials_dt(n_v,dP,inverse_time_interval)
    P_stime = linear_interp_differentials_dt(n_v,dP_s,inverse_time_interval)
    P_ttime = linear_interp_differentials_dt(n_v,dP_t,inverse_time_interval)
 endif
 
end subroutine do_interp_PRZ_2

!> Constructor to allow for optional and default variables
function new_read_jorek_fields_interp_linear(basename, i, rst_format, stop_at_end) result(new)
  character(len=*), intent(in), optional :: basename
  integer, intent(in), optional :: i
  integer, intent(in), optional :: rst_format
  logical, intent(in), optional :: stop_at_end
  type(read_jorek_fields_interp_linear) :: new
  if (present(basename)) new%basename = basename
  if (present(i)) new%i = i
  if (present(rst_format)) new%rst_format = rst_format
  if (present(stop_at_end)) new%stop_at_end = stop_at_end
  new%name = "ReadJorekFieldsInterpLinear"
  new%log = .true.
end function new_read_jorek_fields_interp_linear


!> Find the number of the latest restart file < time (SI units)
!> Perform a bisection method of all the jorek$num.h5 files in the directory
!> Perhaps better to use xtime if this is always present, combined with a filter
!> for all step numbers that are in the current directory
function last_file_before_time(time) result(file_number)
  use phys_module
  use mpi
  real*8, intent(in) :: time
  integer :: file_number
  integer :: i, ierr, my_id, u
  character(len=5) :: my_id_s, num_s
  integer, dimension(:), allocatable :: filenums_tmp, filenums
  integer :: n, i_lower, i_guess, i_upper, io
  real*8 :: t_norm, t_lower, t_guess, t_upper

  call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)
  t_norm = sqrt(mu_zero * mass_proton * central_mass * central_density * 1.d20) ! 1 jorek time unit in seconds
  if (my_id .eq. 0) then
    write(*,*) "Looking for jorek restart file just before time ", time

    ! Get list of filenumbers
    write(my_id_s,"(i0.5)") my_id
    call execute_command_line("ls jorek[0-9]*.h5 | grep -o '[0-9]\{5\}' > .jorek_filenums."//my_id_s)
    open(newunit=u,file=".jorek_filenums."//my_id_s)
    allocate(filenums_tmp(100000)) ! assumes 5-digit numbers
    n=0
    do i=1,100000
      read(u,*,iostat=io) filenums_tmp(i)
      if (io/=0) exit
      n = n + 1
    end do
    allocate(filenums(n))
    filenums(:) = filenums_tmp(1:n)
    deallocate(filenums_tmp)
    close(u, status='delete')

    if (n .le. 0) then
      write(*,*) "No files found!"
      file_number = 0
      return
    end if

    ! Calculate upper and lower bounds
    i_lower = 1 ! index into filenumber array
    write(num_s,'(i0.5)') filenums(i_lower)
    t_lower = get_jorek_hdf5_time('jorek'//num_s//'.h5')*t_norm
    i_upper = n
    write(num_s,'(i0.5)') filenums(i_upper)
    t_upper = get_jorek_hdf5_time('jorek'//num_s//'.h5')*t_norm
    i_guess = nint((time-t_lower)/(t_upper-t_lower)*real(i_upper - i_lower)) + i_lower

    do i=1,20
      if (i_guess .le. 1 .or. i_guess .gt. n) then
        if (my_id .eq. 0) write(*,*) "ERROR: requested time out of range"
        exit
      end if
      if (i_guess .eq. i_lower .or. i_guess .eq. i_upper) exit

      write(num_s,'(i0.5)') filenums(i_guess)
      t_guess = get_jorek_hdf5_time('jorek'//num_s//'.h5')*t_norm
      if (my_id .eq. 0) write(*,"(i5,A,g14.7,A,i5,A,g14.7,A,i5,A,g14.7,A)") i_lower, " (", t_lower, &
        ")    ", i_guess, " (", t_guess, &
        ")    ", i_upper, " (", t_upper, ")    "
      ! Based on the value of t_guess, replace either the lower or upper bound
      if (t_guess .le. time) then
        t_lower = t_guess
        i_lower = i_guess
      else
        t_upper = t_guess
        i_upper = i_guess
      end if
      i_guess = i_lower + (i_upper-i_lower)/2
    end do
    file_number = filenums(i_lower)
    write(*,"(A,i5,A,g14.7,A,g14.7,A)") 'Selected ', i_lower, " (", t_lower, ') as last file before (', time, ')'
  end if
  call MPI_Bcast(file_number, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
end function last_file_before_time

!> Get '/t_now' from a file. Does not alter the units in any way
function get_jorek_hdf5_time(filename) result(time)
  use hdf5
  use hdf5_io_module
  character*(*)      , intent(in)  :: filename
  real*8 :: time
  integer(HID_T) :: file
  integer :: hdferr
  call h5open_f(hdferr)
  call h5fopen_f(filename, H5F_ACC_RDONLY_F, file, hdferr)
  call HDF5_real_reading(file,time,'/t_now')
  call h5fclose_f(file,hdferr)
  call h5close_f(hdferr)
end function get_jorek_hdf5_time



!> Read jorek fields from the next restart file
subroutine do_read(this, sim, ev)
  use mod_import_restart
  use phys_module
  use mpi
  use mod_neighbours
  use nodes_elements
  class(read_jorek_fields_interp_linear), intent(inout) :: this
  type(particle_sim), intent(inout) :: sim
  type(event), intent(inout), optional :: ev
  character(len=80) :: restart_file
  integer :: i, ierr, my_id
  logical :: file_exists, next_file_found

  real*8 :: t_norm
  t_norm = sqrt(mu_zero * mass_proton * central_mass * central_density * 1.d20) ! 1 jorek time unit in seconds

  call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)

  ! Check that the right fields are allocated in sim and allocate if needed
  if (allocated(sim%fields)) then
    select type (f => sim%fields)
    type is (jorek_fields_interp_linear) ! do nothing
    class default
      write(*,*) "WARNING: wrong type of fields in particle%sim, reallocating"
      deallocate(sim%fields)
      allocate(jorek_fields_interp_linear::sim%fields)
    end select
  else
    allocate(jorek_fields_interp_linear::sim%fields)
  end if

  if (.not. associated(sim%fields%node_list))    sim%fields%node_list    => node_list
  if (.not. associated(sim%fields%element_list)) sim%fields%element_list => element_list

  
  ! Continue for jorek_fields_interp_linear
  select type (f => sim%fields)
  type is (jorek_fields_interp_linear)

    if (my_id .eq. 0) then
      ! Read only one file
      if (this%i .eq. -1 .or. f%static) then
        if (this%i .eq. -1) then
          write(restart_file,'(A,A)') trim(this%basename), '_restart.h5'
        else
          write(restart_file,'(A,i5.5,A)') trim(this%basename), this%i, '.h5'
        end if
        inquire(file=trim(restart_file), exist=file_exists)
        if (file_exists) then
          call import_hdf5_restart(f%node_list,f%element_list,restart_file,this%rst_format,ierr)
          f%static = .true.
        else
          if (my_id .eq. 0) write(*,*) "ERROR: file ", trim(restart_file), " does not exist"
          call exit(1)
        end if
      else ! Linearly interpolating case
        ! If nothing has been loaded (i.e. fields%time_prev = 0.d0) load the initial file
        if (abs(f%time_prev) .lt. 1.d-50) then
          write(restart_file,'(A,i5.5,A)') trim(this%basename), this%i, '.h5'
          inquire(file=trim(restart_file), exist=file_exists)
          if (file_exists) then
            call import_hdf5_restart(f%node_list,f%element_list,trim(restart_file),this%rst_format,ierr)
            if (ierr .ne. 0) then
              if (my_id .eq. 0) write(*,*) "ERROR: cannot open restart file"
              call exit(1)
            else
              f%time_now = t_start*t_norm ! set by import_hdf5_restart
              ! Set sim%time to this time also, to start at the right point
              if (sim%time .gt. 1d-16) then ! check if this is the right file if we have already set a time
                if (sim%time .le. f%time_now) then
                  if (my_id .eq. 0) write(*,*) "ERROR: restart file read that is too far in the future"
                end if
              else ! otherwise set the time to the time of this file
                sim%time = f%time_now
              end if
              if (my_id .eq. 0) write(*,"(A,f9.8,A)") "Read initial restart file, set t=", sim%time, " [s]"
            end if
          else
            if (my_id .eq. 0) write(*,*) "ERROR: cannot read initial file ", trim(restart_file)
            call exit(1)
          end if
        end if
        
        ! Find the following file (next timestep number)
        next_file_found=.false.
        do i=this%i+1,this%i+20 ! check 20 files ahead
          write(restart_file,'(A,i5.5,A)') trim(this%basename), i, '.h5'
          inquire(file=trim(restart_file), exist=file_exists)
          if (file_exists) then
            next_file_found=.true.
            call merge_restart(f%node_list, f%element_list, trim(restart_file), this%rst_format,my_id, ierr)
            if (ierr .ne. 0) then
              if (my_id .eq. 0) write(*,*) "ERROR: cannot open restart file"
              call exit(1)
            else
              f%time_prev = f%time_now
              f%time_now = t_start*t_norm ! Set by merge_restart
              if (my_id .eq. 0 .and. i-this%i .gt. 1) write(*,"(i2,A)") i-this%i, " JOREK steps between restarts"
              this%i=i
              ! set the time to run this event at next
              if (my_id .eq. 0) write(*,"(A,f9.8,A)") " Read next restart file, values until t=", f%time_now, " [s]"
              if (present(ev)) then
                ev%start = f%time_now
                if (my_id .eq. 0) write(*,*) "Set time for next restart file read to ", ev%start
              end if
              exit ! the file-finding loop
            end if
          endif
        enddo
        if (.not. next_file_found .and. this%stop_at_end) then
          if (my_id .eq. 0) write(*,*) "WARNING: cannot find any next restart files. Stopping."
          call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        if (.not. next_file_found .and. .not. this%stop_at_end) then
          if (my_id .eq. 0) write(*,*) "WARNING: cannot find any next restart files. Continuing with &
          the last values without time-dependence"
          f%static = .true.
        end if
      end if

    end if

    ! Communicate the new fields to all processes
    call broadcast_elements(my_id, f%element_list)
    call broadcast_nodes(my_id, f%node_list)
    call update_neighbours(f%node_list, f%element_list) ! needs to be done on every process to have an RTree everywhere
    call MPI_Bcast(f%time_prev, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(f%time_now, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(sim%time, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(t_start, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr) ! might be needed in some cases (static with number)
    call MPI_Bcast(f%static, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    if (present(ev)) then
      call MPI_Bcast(ev%start, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    end if
  class default
    if (my_id .eq. 0) write(*,*) "ERROR, do_read called with wrong sim%fields"
  end select
end subroutine do_read



!> Import a binary restart file and merges it with the values currently known
!> This can then be used to interpolate linearly between any two restart files
subroutine merge_restart(node_list,element_list, restart_file, format_rst,my_id, ierr)
  use data_structure
  use phys_module
  use mod_import_restart
  implicit none

  ! --- Routine parameters
  type(type_node_list),    intent(inout) :: node_list
  type(type_element_list), intent(inout) :: element_list
  character(len=*),        intent(in)    :: restart_file !< Filename of new restart file to import
  integer,                 intent(out)   :: ierr
  integer,                 intent(in)    :: format_rst !< Restart file format
  integer,                 intent(in)    :: my_id

  ! --- Internal variables
  real*8, allocatable, dimension(:,:,:,:) :: values
  integer :: inode
  real*8 :: tstart_old

  ! Save the old values to calculate the new deltas
  allocate(values(n_tor,n_order+1,n_var,node_list%n_nodes))
  !$omp parallel do default(shared) private(inode)
  do inode=1,node_list%n_nodes
    values(:,:,:,inode) = node_list%node(inode)%values(:,:,:)
  enddo
  !$omp end parallel do
  tstart_old = t_start

  ! Import new values
  call import_hdf5_restart(node_list,element_list, restart_file, format_rst, ierr)

  ! Calculate deltas as values_new - values_old
  !$omp parallel do default(shared) private(inode)
  do inode=1,node_list%n_nodes
    node_list%node(inode)%deltas = node_list%node(inode)%values - values(:,:,:,inode)
  enddo
  !$omp end parallel do

  ! Set timestep to time between restart files
  tstep = t_start - tstart_old

  deallocate(values)
end subroutine merge_restart
end module mod_fields_linear
