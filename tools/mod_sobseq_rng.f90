!> Wrapper in the mod_rng type for [[mod_sobseq]]
module mod_sobseq_rng
  use mod_sobseq
  use mod_rng
  implicit none
  private
  public type_rng, sobseq_rng

  type, extends(type_rng) :: sobseq_rng
    private
    integer :: n_streams
    integer :: i_stream
    type(sobol_state), dimension(:), allocatable :: state
    contains
      procedure, public :: initialize => initialize_sobseq_rng
      procedure, public :: next => next_sobseq_rng
      procedure, public :: jump_ahead => jump_ahead_sobseq_rng
  end type
  interface sobseq_rng
    module procedure empty_sobseq_rng
  end interface sobseq_rng

  ! Direction numbers taken from Joe and Kuo
  integer, parameter, dimension(1:8)   :: s = (/1,2,3,3,4,4,5,5/)
  integer, parameter, dimension(1:8)   :: a = (/0,1,1,2,1,4,2,4/)
  integer, parameter, dimension(5,1:8) :: m = reshape((/1,0,0,0,0, &
                                                1,3,0,0,0, &
                                                1,3,1,0,0, &
                                                1,1,1,0,0, &
                                                1,1,3,3,0, &
                                                1,3,5,13,0,&
                                                1,1,5,5,17,&
                                                1,1,5,5,5/), (/5,8/))
  public :: ilog2_b_ceil
contains
  type(sobseq_rng) function empty_sobseq_rng()
  end function empty_sobseq_rng

  subroutine initialize_sobseq_rng(rng, n_dims, seed, n_streams, i_stream, ierr)
    use mpi
    implicit none
    class(sobseq_rng), intent(inout) :: rng
    integer, intent(in)  :: n_dims !< Dimension of the generated output vector (not used)
    integer, intent(in)  :: seed !< Seed for the RNG if required
    integer, intent(in)  :: n_streams !< Number of output streams needed (not used)
    integer, intent(in)  :: i_stream !< Index of this output stream (sequence number)
    integer, intent(out), optional :: ierr !< Error code. If present returns, otherwise calls MPI_ABORT

    integer :: i, ifail, mpi_err
    real*8 :: DUMMY_REAL

    ifail = 0
    if (n_dims .le. 0) then
      write(*,*) "Number of dimensions must be greater than 0"
      ifail = 1
    else if (n_dims .gt. 8) then
      write(*,*) "Number of dimensions greater than 8 not implemented yet"
      ifail = 6
    else if (n_streams .le. 0) then
      write(*,*) "Number of streams must be greater than zero"
      ifail = 3
    else if (i_stream .le. 0) then
      write(*,*) "Index of this stream must be greater than zero"
      ifail = 4
    else if (i_stream .gt. n_streams) then
      write(*,*) "Index of this stream must be less than number of streams"
      ifail = 5
    else if (n_streams .lt. 2**ilog2_b_ceil(n_streams)) then
      write(*,*) "Number of streams must be a power of 2"
      ifail = 6
    endif
    if (present(ierr)) ierr = ifail
    if (ifail .gt. 0) then
      if (present(ierr)) then
        return
      else
        call MPI_ABORT(MPI_COMM_WORLD, ierr, mpi_err)
      end if
    end if

    if (allocated(rng%state)) deallocate(rng%state)
    allocate(rng%state(n_dims))
    rng%n_streams = n_streams
    rng%i_stream = i_stream

    ! Seed with default values
    do i=1,n_dims
      call rng%state(i)%initialize(s(i), a(i), m(:,i), stride=ilog2_b_ceil(n_streams))
      DUMMY_REAL = rng%state(i)%skip_ahead(i_stream)
    end do
  end subroutine initialize_sobseq_rng

  !> Generate the next value in the sobol series for n_streams.
  !> This routine contains some extra logic to handle n_streams not being a
  !> power of 2, in which case n_streams values need to be selected out of 2^ilog2_b_ceil(n_streams)
  subroutine next_sobseq_rng(rng, out)
    implicit none
    class(sobseq_rng), intent(inout)  :: rng
    real*8, dimension(:), intent(out) :: out

    integer :: i, N

    if (rng%n_streams .eq. 1) then
      do i=1,size(rng%state,1)
        out(i) = rng%state(i)%next()
      end do
    else
      N = ilog2_b_ceil(rng%n_streams)
      if (rng%n_streams .lt. 2**N) then
        ! We have to select n_streams values from 2^N values
        write(*,*) "ERROR: selecting ", rng%n_streams, " values from ", 2**N, " values not implemented"
        !call exit(137)
      else
        do i=1,size(rng%state,1)
          out(i) = rng%state(i)%next_strided()
        end do
      end if
    end if
  end subroutine next_sobseq_rng

  !> Jump ahead many steps fast
  subroutine jump_ahead_sobseq_rng(rng, delta)
    implicit none
    class(sobseq_rng), intent(inout) :: rng
    integer, intent(in)              :: delta
    real*8 :: out
    integer :: i
    do i=1,size(rng%state,1)
      out = rng%state(i)%jump_ahead(delta)
    end do
  end subroutine jump_ahead_sobseq_rng

  !> Integer logarithm base 2 rounded up (i.e. the smallest n in 2^n >= val)
  function ilog2_b_ceil(val) result(res)
    integer, intent(in) :: val
    integer             :: res
    integer             :: tmp

    res = -1
    ! Negative values not allowed
    if (val < 1) return

    tmp = val
    do while (tmp > 0)
      res = res + 1
      tmp = rshift(tmp, 1)
    end do
    ! round up if not an exact fit
    if (val .gt. 2**res) res = res + 1
  end function ilog2_b_ceil
end module mod_sobseq_rng
