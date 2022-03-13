!> Wrapper module for the pcg32 rng in the mod_rng type
module mod_pcg32_rng
use mod_pcg32
use mod_rng
implicit none
private
public pcg32_rng

type, extends(type_rng) :: pcg32_rng
#ifndef __NVCOMPILER 
  private
#endif
    type(pcg_state_setseq_64) :: state = pcg_state_setseq_64( &
      state = 234617468, inc = 123674237) !< mash-they-keyboard random
  contains
    procedure :: initialize => initialize_pcg32_rng
    procedure :: next => next_pcg32_rng
    procedure :: jump_ahead => jump_ahead_pcg32_rng
end type

contains
  subroutine initialize_pcg32_rng(rng, n_dims, seed, n_streams, i_stream, ierr)
    use mpi
    implicit none
    class(pcg32_rng), intent(inout) :: rng
    integer, intent(in)  :: n_dims !< Dimension of the generated output vector (not used)
    integer, intent(in)  :: seed !< Seed for the RNG if required
    integer, intent(in)  :: n_streams !< Number of output streams needed (not used)
    integer, intent(in)  :: i_stream !< Index of this output stream (sequence number)
    integer, intent(out), optional :: ierr !< Error code. Exit on error if not present
    integer :: ifail, mpi_err

    ifail = 0
    if (n_dims .le. 0) ifail = 1
    if (seed .eq. 0) ifail = 2

    if (present(ierr)) ierr = ifail

    if (ifail .gt. 0) then
      if (present(ierr)) then
        return
      else
        call MPI_ABORT(MPI_COMM_WORLD, ierr, mpi_err)
      end if
    end if

    call pcg32_srandom_r(rng%state, int(seed, 8), int(i_stream, 8))
  end subroutine initialize_pcg32_rng

  subroutine next_pcg32_rng(rng, out)
    implicit none
    class(pcg32_rng), intent(inout) :: rng
    real*8, dimension(:), intent(out) :: out
    call pcg32_random_doubles_r(rng%state, out)
  end subroutine next_pcg32_rng

  subroutine jump_ahead_pcg32_rng(rng, delta)
    class(pcg32_rng), intent(inout) :: rng
    integer, intent(in) :: delta
    call pcg32_jumpahead(rng%state, delta)
  end subroutine jump_ahead_pcg32_rng
end module mod_pcg32_rng
