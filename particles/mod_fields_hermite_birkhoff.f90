!> Module for Hermite-Birkhoff spline interpolation in time
!> in JOREK restart files. Contains an action to read the fields.
!>
!> This module calculates a finite-difference approximation to the time-derivative
!> and stores this in the deltas.
!>
!> We create a ring buffer https://en.wikipedia.org/wiki/Circular_buffer to
!> store the node_lists we will interpolate from.
!> This contains the values at 3 or 4 different timesteps, with derivatives
!> calculated for all (i.e. 1 or 2) intermediate timesteps.
!> We implement this as a starting index into the array for the oldest stored
!> timestep, plus a length. We also store the times.
!>
!> WARNING: There is no guarantee that fields%node_list contains the values at the
!> current timestep. It will contain the most recently read node_list.
!> The only use for this is to get the structure of the grid,
!> which is assumed not to change. For all other uses please call the interp_PRZ
!> function.
module mod_fields_hermite_birkhoff
use data_structure
use mod_particle_sim
use mod_event
use mod_fields
use mod_interp
use mod_fields_linear
use mod_hermite_birkhoff
implicit none
private
public jorek_fields_interp_hermite_birkhoff, read_jorek_fields_interp_hermite_birkhoff
public node_list_deltas_to_values

!> Action to read in the fields into sim%fields
type, extends(action) :: read_jorek_fields_interp_hermite_birkhoff
  character(len=80) :: basename = 'jorek' !< Comes before the file number or extension
  integer :: i = 0 !< Number of the restart file to read. Set to -1 to not include. Corresponds to the index of NOW
  integer :: rst_format = 0 !< Format of restart file if .rst type
  logical :: stop_at_end = .true. !< Whether to stop the simulation at the end of the file list
  contains
    procedure :: do => do_read
end type
interface read_jorek_fields_interp_hermite_birkhoff
  module procedure new_read_jorek_fields_interp_hermite_birkhoff
end interface read_jorek_fields_interp_hermite_birkhoff

!> Store enough data for interpolation in time.
!> This does not work for changing grids in time!
!> Use at your own peril in that case. (e.g. i_elm and s,t for a specific spatial position might depend on time)
!>
!> We use a kind of ring buffer to store the node lists
integer, parameter :: NL = 4 !< number of node_lists
type, extends(fields_base) :: jorek_fields_interp_hermite_birkhoff
  type(type_node_list), allocatable, dimension(:)    :: node_lists    !< Ring buffer of node lists
  type(type_element_list), allocatable, dimension(:) :: element_lists !< Ring buffer of element lists
  real*8, dimension(NL) :: t !< Time at each restart file (SI units)

  integer :: start = 1 !< Index in the ring buffer of the oldest file
  integer :: len = 0
  contains
    procedure :: interp_PRZ => do_interp_PRZ_1
    procedure :: interp_PRZ_2 => do_interp_PRZ_2
    procedure, nopass :: mask !< Map number to 1..NL
    procedure :: ind !< Return (array of) indices into 1..NL (start+i etc)
end type jorek_fields_interp_hermite_birkhoff
 
contains
!> The correct indexing operation for our ring buffer
pure function mask(i)
  integer, intent(in) :: i
  integer :: mask
  mask = mod(i-1,NL)+1
end function mask
!> Convert an index in the ring buffer into an array index
elemental function ind(f, j)
  class(jorek_fields_interp_hermite_birkhoff), intent(in) :: f
  integer, intent(in) :: j !< fortran-style index (1-based)
  integer :: ind !< Array index for node_lists etc corresponding to element j
  ind = mask(f%start+j-1)
end function ind


!> Interpolate a variable at a specific position (with phi), with first derivatives only
pure subroutine do_interp_PRZ_1(this, time, i_elm, i_v, n_v, s, t, phi, P, P_s, P_t, P_phi, P_time, R, R_s, R_t, Z, Z_s, Z_t)
  class(jorek_fields_interp_hermite_birkhoff),  intent(in)  :: this
  real*8,                   intent(in)  :: time !< Time at which to calculate this variable
  integer,                  intent(in)  :: i_elm
  integer,                  intent(in)  :: n_v, i_v(n_v)
  real*8,                   intent(in)  :: s, t, phi
  real*8,                   intent(out) :: P(n_v), P_s(n_v), P_t(n_v), P_phi(n_v), P_time(n_v)
  real*8,                   intent(out) :: R, R_s, R_t, Z, Z_s, Z_t

  real*8, dimension(n_v,4,4) :: V !< n_var, (P, P_s, P_t, P_phi), (interp, interp, delta, delta)
  integer :: i1, i2, j

  !> initialise deltas to zero
  V(:,:,3:4) = 0.d0

  ! Determine between which sets to interpolate by selecting
  ! t_i <= t_now and taking the last true and first false value.
  do j=1,this%len
    if (this%t(this%ind(j)) .gt. time) exit ! the loop
  end do
  i1 = this%ind(j-1)
  i2 = this%ind(j)
#ifdef DEBUG
  if (j .eq. this%len) i2 = 99999 ! intentional segfault for debugging
#endif

  call interp_PRZ(this%node_lists(i1),  this%element_lists(i1),  i_elm,i_v,n_v,s,t,phi, &
    V(:,1,1), V(:,2,1), V(:,3,1), V(:,4,1),   R,R_s,R_t,Z,Z_s,Z_t)
  call interp_PRZ(this%node_lists(i2),  this%element_lists(i2),  i_elm,i_v,n_v,s,t,phi, &
       V(:,1,2), V(:,2,2), V(:,3,2), V(:,4,2),   R,R_s,R_t,Z,Z_s,Z_t)
  call interp_PRZ(this%node_lists(i1),  this%element_lists(i1),  i_elm,i_v,n_v,s,t,phi, &
    V(:,1,3), V(:,2,3), V(:,3,3), V(:,4,3),   R,R_s,R_t,Z,Z_s,Z_t,deltas=.true.)
  call interp_PRZ(this%node_lists(i2),  this%element_lists(i2),  i_elm,i_v,n_v,s,t,phi, &
    V(:,1,4), V(:,2,4), V(:,3,4), V(:,4,4),   R,R_s,R_t,Z,Z_s,Z_t,deltas=.true.)

  call HB_interp(this%t(i1), this%t(i2), n_v, V(:,1,1), V(:,1,2), V(:,1,3), V(:,1,4), time, P)
  call HB_interp(this%t(i1), this%t(i2), n_v, V(:,2,1), V(:,2,2), V(:,2,3), V(:,2,4), time, P_s)
  call HB_interp(this%t(i1), this%t(i2), n_v, V(:,3,1), V(:,3,2), V(:,3,3), V(:,3,4), time, P_t)
  call HB_interp(this%t(i1), this%t(i2), n_v, V(:,4,1), V(:,4,2), V(:,4,3), V(:,4,4), time, P_phi)
  call HB_interp_dt(this%t(i1), this%t(i2), n_v, V(:,1,1), V(:,1,2), V(:,1,3), V(:,1,4), time, P_time)
end subroutine do_interp_PRZ_1

!> This procedure interpolates a variable, its first and second order derivatives in space
!> and first order derivatives in time. Only first order spatial derivatives in
!> s and t are also derived in time.
!> inputs:
!>   this:    (jorek_fields_interp_hermite_birkhoff) interpolation class
!>   time:    (real8) current time
!>   i_elm:   (integer) element index
!>   i_v:     (integer)(n_v) array of jorek field indices
!>   n_v:     (integer) number of jorek fields to interpolate
!>   s,t,phi: (real8) interpolation position in local (s,t) coordinates and
!>            toroidal angle (phi)
!> outputs:
!>  P:               (real8)(n_v) interpolated quantities
!>  P_s,P_t,P_time:  (real8)(n_v) interpolated quantity first order derivatives
!>  P_ss,P_st,P_tt:  (real8)(n_v) interpolated quantity second order derivatives
!>  P_sphi,P_tphi:   (real8)(n_v) interpolated quantitiy second order derivatives
!>  P_stime,P_ttime: (real8)(n_v) interpolated quantity cross time-space derivatives
!>  R,Z:             (real8) radial and vertical position
!>  R_s,R_t,Z_s,Z_t: (real8) radial and vertical position first order derivatives
!>  R_ss,R_st,R_tt:  (real8) R second order derivatives
!>  Z_ss,Z_st,Z_tt:  (real8) Z second order derivatives
pure subroutine do_interp_PRZ_2(this,time,i_elm,i_v,n_v,s,t,phi,       &
  P,P_s,P_t,P_phi,P_time,P_ss,P_st,P_tt,P_sphi,P_tphi,P_stime,P_ttime, &
  R,R_s,R_t,R_ss,R_st,R_tt,Z,Z_s,Z_t,Z_ss,Z_st,Z_tt)
  implicit none
  !> declare input variables
  class(jorek_fields_interp_hermite_birkhoff), intent(in) :: this
  real(kind=8), intent(in)                                :: time, s, t, phi
  integer, intent(in)                                     :: i_elm,n_v
  integer, dimension(n_v), intent(in)                     :: i_v
  !> declare output variables
  real(kind=8), intent(out) :: R, R_s, R_t, R_ss, R_st, R_tt
  real(kind=8), intent(out) :: Z, Z_s, Z_t, Z_ss, Z_st, Z_tt
  real(kind=8), dimension(n_v), intent(out) :: P, P_s, P_t, P_phi, P_time
  real(kind=8), dimension(n_v), intent(out) :: P_ss, P_st, P_tt, P_sphi, P_tphi
  real(kind=8), dimension(n_v), intent(out) :: P_stime, P_ttime
  !> declare internal variables
  !> values array: first index: jorek field to interpolate
  !>               second index: (P,P_s,P_t,P_phi,P_st,P_ss,P_tt,P_sphi,P_tphi,P_phiphi)
  !>               third index:  (value_restart_1,value_restart_2,deltas_1,deltas_2)
  real(kind=8), dimension(n_v,10,4) :: values
  integer                           :: i1,i2,j

  !> initialise deltas to zero
  values(:,:,3:4) = 0.d0

  !> Determine between which restart files the interpolation occurs
  !> selection: t_i <= t_now and taking the last true and first false values.
  do j=1,this%len
     if(this%t(this%ind(j)).gt.time) exit !< found the first restart after time
  enddo
  i1 = this%ind(j-1) !< earlier restart
  i2 = this%ind(j)   !< later restart
#ifdef DEBUG
    if(j.eq.this%len) i2=9999999 !< intentional segfault for debugging
#endif

    !> field spatial interpolation
    !>   first restart values
    call interp_PRZ(this%node_lists(i1),this%element_lists(i1),i_elm,i_v,n_v, &
      s,t,phi,values(:,1,1),values(:,2,1),values(:,3,1),values(:,4,1),        &
      values(:,5,1),values(:,6,1),values(:,7,1),values(:,8,1),                &
      values(:,9,1),values(:,10,1),R,R_s,R_t,R_st,R_ss,R_tt,                  &
      Z,Z_s,Z_t,Z_st,Z_ss,Z_tt) 
    !>   second restart values                                     
    call interp_PRZ(this%node_lists(i2),this%element_lists(i2),i_elm,i_v,n_v, &
      s,t,phi,values(:,1,2),values(:,2,2),values(:,3,2),values(:,4,2),        &
      values(:,5,2),values(:,6,2),values(:,7,2),values(:,8,2),                &
      values(:,9,2),values(:,10,2),R,R_s,R_t,R_st,R_ss,R_tt,                  &
      Z,Z_s,Z_t,Z_st,Z_ss,Z_tt)  
    !> first restart derivatives                                    
    call interp_PRZ(this%node_lists(i1),this%element_lists(i1),i_elm,i_v,n_v, &
      s,t,phi,values(:,1,3),values(:,2,3),values(:,3,3),values(:,4,3),        &
      values(:,5,3),values(:,6,3),values(:,7,3),values(:,8,3),                &
      values(:,9,3),values(:,10,3),R,R_s,R_t,R_st,R_ss,R_tt,                  &
      Z,Z_s,Z_t,Z_st,Z_ss,Z_tt,deltas=.true.)
    !> second restart derivatives
    call interp_PRZ(this%node_lists(i2),this%element_lists(i2),i_elm,i_v,n_v, &
      s,t,phi,values(:,1,4),values(:,2,4),values(:,3,4),values(:,4,4),        &
      values(:,5,4),values(:,6,4),values(:,7,4),values(:,8,4),                &
      values(:,9,4),values(:,10,4),R,R_s,R_t,R_st,R_ss,R_tt,                  &
      Z,Z_s,Z_t,Z_st,Z_ss,Z_tt,deltas=.true.) 
 
    !> compute field time interpolations
    call HB_interp(this%t(i1),this%t(i2),n_v,values(:,1,1),values(:,1,2), &
      values(:,1,3),values(:,1,4),time,P)
    call HB_interp(this%t(i1),this%t(i2),n_v,values(:,2,1),values(:,2,2), &
      values(:,2,3),values(:,2,4),time,P_s)
    call HB_interp(this%t(i1),this%t(i2),n_v,values(:,3,1),values(:,3,2), &
      values(:,3,3),values(:,3,4),time,P_t)
    call HB_interp(this%t(i1),this%t(i2),n_v,values(:,4,1),values(:,4,2), &
      values(:,4,3),values(:,4,4),time,P_phi)
    call HB_interp(this%t(i1),this%t(i2),n_v,values(:,5,1),values(:,5,2), &
      values(:,5,3),values(:,5,4),time,P_st)
    call HB_interp(this%t(i1),this%t(i2),n_v,values(:,6,1),values(:,6,2), &
      values(:,6,3),values(:,6,4),time,P_ss)
    call HB_interp(this%t(i1),this%t(i2),n_v,values(:,7,1),values(:,7,2), &
      values(:,7,3),values(:,7,4),time,P_tt)
    call HB_interp(this%t(i1),this%t(i2),n_v,values(:,8,1),values(:,8,2), &
      values(:,8,3),values(:,8,4),time,P_sphi)
    call HB_interp(this%t(i1),this%t(i2),n_v,values(:,9,1),values(:,9,2), &
      values(:,9,3),values(:,9,4),time,P_tphi)
    
    !> compute field time derivatives
    call HB_interp_dt(this%t(i1),this%t(i2),n_v,values(:,1,1),values(:,1,2), &
      values(:,1,3),values(:,1,4),time,P_time)
    call HB_interp_dt(this%t(i1),this%t(i2),n_v,values(:,2,1),values(:,2,2), &
      values(:,2,3),values(:,2,4),time,P_stime)
    call HB_interp_dt(this%t(i1),this%t(i2),n_v,values(:,3,1),values(:,3,2), &
      values(:,3,3),values(:,4,3),time,P_ttime)
    
end subroutine do_interp_PRZ_2

!> Constructor to allow for optional variables
function new_read_jorek_fields_interp_hermite_birkhoff(basename, i, rst_format, stop_at_end) result(new)
  character(len=*), intent(in), optional :: basename
  integer, intent(in), optional :: i
  integer, intent(in), optional :: rst_format
  logical, intent(in), optional :: stop_at_end
  type(read_jorek_fields_interp_hermite_birkhoff) :: new
  if (present(basename)) new%basename = basename
  if (present(i)) new%i = i
  if (present(rst_format)) new%rst_format = rst_format
  if (present(stop_at_end)) new%stop_at_end = stop_at_end
  new%name = "ReadJorekFieldsInterpHermiteBirkhoff"
  new%log = .true.
  !< check if a single restart is required and exit the program if this is the case
  if(i.eq.(-1)) then 
    write(*,*) "ERROR: Hermite-Birkhoff interpolation requires several restart files."
    write(*,*) "For a single restart file, use the linear interpolation."
    stop
  endif
end function new_read_jorek_fields_interp_hermite_birkhoff

!> Read jorek fields from the next restart file.
!>
!> This routine needs to handle setting up the first two files and reading the
!> next file in the sequence.
!>
!> For the initial read we read the first file, containing info about
!> timestep i (values) and i-1 (deltas). This is transformed into two node
!> lists in the ring buffer.

!> Any subsequent read is simple. We read the next file i(jend)+N (try N=2,1,3+)
!> containing i+N (values) and i+N-1 (deltas) (N>=1).
!> If N == 1 we have overlap and only need to store the values in the third slot
!> in the ring buffer. If N > 1 we can store two node lists and the buffer is
!> full.
!> Now we can calculate the derivatives of the middle one or two fields. These
!> we store in the deltas (changing the meaning from delta to derivative!).
!> The oldest element in the ring buffer can now be reclaimed and we can
!> interpolate from the next oldest element onwards.
!> The next read action is set for the time when we run out of derivatives, i.e.
!> t(jend-1).
subroutine do_read(this, sim, ev)
  use phys_module, only: central_mass, central_density, t_start, tstep
  use constants, only: mu_zero, mass_proton
  use mpi
  use mod_neighbours
  class(read_jorek_fields_interp_hermite_birkhoff), intent(inout) :: this
  type(particle_sim), intent(inout) :: sim
  type(event), intent(inout), optional :: ev
  character(len=80) :: restart_file
  integer :: i, j, k, di, ierr, my_id, i1, is(3), dt(2)
  logical :: file_exists, next_file_found

  real*8 :: t_norm, invdet
  t_norm = sqrt(mu_zero * mass_proton * central_mass * central_density * 1.d20) ! 1 jorek time unit in seconds

  call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)

  ! Check that the right fields are allocated in sim and allocate if needed
  if (allocated(sim%fields)) then
    select type (f => sim%fields)
    type is (jorek_fields_interp_hermite_birkhoff) ! do nothing
    class default
      write(*,*) "WARNING: wrong type of fields in particle%sim, reallocating"
      deallocate(sim%fields)
      allocate(jorek_fields_interp_hermite_birkhoff::sim%fields)
    end select
  else
    allocate(jorek_fields_interp_hermite_birkhoff::sim%fields)
  end if
  if (.not. associated(sim%fields%node_list))    allocate(sim%fields%node_list)
  if (.not. associated(sim%fields%element_list)) allocate(sim%fields%element_list)
  
  ! Continue for jorek_fields_interp_hermite_birkhoff
  select type (f => sim%fields)
  type is (jorek_fields_interp_hermite_birkhoff)
    if (.not. allocated(f%node_lists))    allocate(f%node_lists(NL))
    if (.not. allocated(f%element_lists)) allocate(f%element_lists(NL))

    ! If nothing has been loaded load the initial file
    if (f%len .eq. 0) then
      this%i = this%i-1 ! trick to reuse normal reading code
      call read_next_file(this, f, i, prefer_plus_2=.false.)
      this%i = i
      call update_neighbours(f%node_list, f%element_list)
      call append_to_fields(f, f%node_list, f%element_list, t_start*t_norm, &
        tstep*t_norm, from_deltas=.true.)
      ! Set sim%time to this time also, to start at the right point
      if (sim%time .gt. 1d-16) then ! check if this is the right file if we have already set a time
        if (sim%time .lt. f%t(f%ind(f%len))) then
          if (my_id .eq. 0) then
            write(*,*) "ERROR: restart file read that is too far in the future"
            call exit(1)
          end if
        end if
      else ! otherwise set the time to the time of this file
        sim%time = f%t(f%ind(f%len))
      end if

      ! Now read the next file
      call read_next_file(this, f, i, prefer_plus_2=.true.)
      call update_neighbours(f%node_list, f%element_list)
      call append_to_fields(f, f%node_list, f%element_list, t_start*t_norm, &
        tstep*t_norm, from_deltas=i - this%i .ge. 2)
      ! note that t_start is set in import_hdf5_restart instead of t_now
      this%i = i ! update index of latest file read
      ! Now we have 3 or 4 time points in our list

      ! Finally we need to calculate the derivatives of the middle points (1 or
      ! 2 points)
      do j=2,f%len-1
        call interp_derivatives(f, j)
      end do

      ! Now we can remove the first element
      f%start = f%ind(2)
      f%len = f%len-1
      ! And we finish with 2 or 3 timesteps in the buffer

      if (my_id .eq. 0) write(*,"(A,f9.8,A)") "Read initial restart file, set t=", sim%time, " [s]"
      if (present(ev)) then
        ev%start = f%t(f%ind(f%len-1)) ! read again at next-to-last file
        if (my_id .eq. 0) write(*,*) "Set time for next restart file read to ", ev%start
      end if


    else

      ! Incremental read. We are now roughly at f%t(f%ind(f%len-1)) (i.e. time
      ! of next-to-last step)

      ! Drop old steps. All except the last 2
      f%start = f%ind(f%len-1)
      f%len = 2


      ! Do an incremental read
      call read_next_file(this, f, i, prefer_plus_2=.true.)
      call update_neighbours(f%node_list, f%element_list)
      if (my_id .eq. 0) write(*,"(i2,A)") i-this%i, " JOREK steps between restarts"
      call append_to_fields(f, f%node_list, f%element_list, t_start*t_norm, &
        tstep*t_norm, from_deltas=i - this%i .ge. 2)
        ! note that t_start is set in import_hdf5_restart instead of t_now
      this%i=i ! set index of last-read file

      !> interpolate fields with midpoint rule if more than one restart is used
      do j=2,f%len-1
        call interp_derivatives(f, j)
      end do

      ! set the time to run this event at next
      if (my_id .eq. 0) write(*,"(A,f9.8,A)") " Read next restart file, values until t=", f%t(f%ind(f%len-1)), " [s]"
      if (present(ev)) then
        ev%start = f%t(f%ind(f%len-1)) ! read again at next-to-last file
        if (my_id .eq. 0) write(*,*) "Set time for next restart file read to ", ev%start
      end if
    end if

  class default
    if (my_id .eq. 0) write(*,*) "ERROR, do_read called with wrong sim%fields"
  end select
end subroutine do_read

!> Starting from i+2 (if prefer_plus_2=.true.), i+1, i+N search for a next file and read it into
!> f%node_list, f%element_list.
!> Performs MPI communication to get values from root process to every other process.
!> Also broadcasted are tstep and t_start
subroutine read_next_file(this, f, i_found, prefer_plus_2)
  use mod_import_restart
  use phys_module
  use mpi
  class(read_jorek_fields_interp_hermite_birkhoff), intent(inout) :: this
  class(jorek_fields_interp_hermite_birkhoff), intent(inout) :: f
  integer, intent(out) :: i_found
  logical, optional, intent(in) :: prefer_plus_2 !< Set to true if we need to
  !< check for the presence of this%i+2 first

  character(len=80) :: restart_file
  integer :: i, j, k, di, ierr, my_id
  logical :: file_exists, next_file_found, flip_i12 = .false.

  call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)
  if (present(prefer_plus_2)) then
    if (prefer_plus_2) flip_i12 = .true.
  end if

  if (my_id .eq. 0) then
    ! Find the following file (next timestep number)
    next_file_found=.false.
    do di=1,20
      if (di .le. 2 .and. flip_i12) then ! flip first and second element to search to save effort reading
        i = this%i + (3-di)
      else
        i = this%i + di
      end if
      write(restart_file,'(A,i5.5,A)') trim(this%basename), i, '.h5'
      inquire(file=trim(restart_file), exist=file_exists)
      if (file_exists) then
        next_file_found=.true.

        ! Read the next restart file
        call import_hdf5_restart(f%node_list,f%element_list,trim(restart_file),this%rst_format,ierr)
        call broadcast_elements(my_id, f%element_list)
        call broadcast_nodes(my_id, f%node_list)
        call broadcast_phys(my_id) ! we only really use tstep and t_start but this is simpler to write
        if (ierr .ne. 0) then
          if (my_id .eq. 0) write(*,*) "ERROR: cannot open restart file"
          call exit(1)
        else
          exit ! the file-finding loop
        end if
      endif
    enddo
    if (.not. next_file_found) then
      if (my_id .eq. 0) write(*,*) "WARNING: cannot find any next restart files. Stopping."
      call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
      i_found = 0
    else
      i_found = i
    end if
    call MPI_Bcast(i_found, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  else
    call broadcast_elements(my_id, f%element_list)
    call broadcast_nodes(my_id, f%node_list)
    call broadcast_phys(my_id) ! we only really use tstep and t_start but this is simpler to write
    call MPI_Bcast(i_found, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  end if
end subroutine read_next_file

!> Interpolate derivatives at a time j by nonuniform finite differences
!> j must be between 2 and f%len-1
pure subroutine interp_derivatives(f, j)
  class(jorek_fields_interp_hermite_birkhoff), intent(inout) :: f
  integer, intent(in) :: j !< index of time to interpolate at. Must be between
  !< 2 and f%len-1
  integer :: k, is(3)
  real*8 :: invdet, dt(2)

  is = f%ind([j-1,j,j+1])
  dt = [f%t(is(2))-f%t(is(1)), f%t(is(3))-f%t(is(2))]
  invdet = 1.d0/(dt(1)*(dt(1)+dt(2))*dt(2))
  do k=1,f%node_lists(is(2))%n_nodes
    f%node_lists(is(2))%node(k)%deltas = (&
        -dt(2)**2*f%node_lists(is(1))%node(k)%values &
      + (dt(2)**2-dt(1)**2)*f%node_lists(is(2))%node(k)%values &
      + dt(1)**2*f%node_lists(is(3))%node(k)%values) &
      * invdet
  end do
end subroutine interp_derivatives

!> Append a node and element list to the fields array
pure subroutine append_to_fields(f, node_list, element_list, t, dt, from_deltas)
  class(jorek_fields_interp_hermite_birkhoff), intent(inout) :: f
  type(type_node_list), intent(in) :: node_list
  type(type_element_list), intent(in) :: element_list
  real*8, intent(in) :: t, dt !< Time and timestep in SI units
  logical, intent(in) :: from_deltas
  integer :: i1, i2

  if (.not. from_deltas) then
    i1 = f%ind(f%len+1)
    f%node_lists(i1)    = node_list
    f%element_lists(i1) = element_list
    f%t(i1)             = t
    f%len = f%len+1
  else
    i1 = f%ind(f%len+1)
    i2 = f%ind(f%len+2)
    f%node_lists(i1)         = node_list
    f%element_lists(i1)      = element_list
    call node_list_deltas_to_values(f%node_lists(i1))
    f%node_lists(i2)    = node_list
    f%element_lists(i2) = element_list
    f%t(i1) = t - dt
    f%t(i2) = t
    f%len = f%len+2
  end if
end subroutine append_to_fields

!> Subtract deltas from values to get the values at the previous timestep.
!> Deltas = f(i) - f(i-1) => f(i-1) = f(i) - deltas
pure subroutine node_list_deltas_to_values(node_list)
  type(type_node_list), intent(inout) :: node_list
  integer :: i
  do i=1,node_list%n_nodes
    node_list%node(i)%values = node_list%node(i)%values - node_list%node(i)%deltas
  end do
end subroutine node_list_deltas_to_values
end module mod_fields_hermite_birkhoff
