!> the module mod_rng_mpi_tools contains variables
!> and procedures needed for using rng in parallel environments
module mod_rng_mpi_tools
implicit none
private
public :: init_parallel_rngs
!> Variables and datatypes -------------------------------------
!> Interfaces --------------------------------------------------
contains
!> Procedures --------------------------------------------------
!> initialise random number generators in parallel environment
!> inputs:
!>   rngs:     (type_rng) random number generators to be initialised
!>   n_rngs:   (integer) number of random number generators
!>   rng_size: (integer) size of the random number vectors
!>   n_tasks:  (integer) number of MPI tasks
!>   my_id:    (integer) id of the current task
!>   ifail:    (integer) failing code: 0-success 
!> outputs:
!>   rngs: (type_rng) initialised random number generators
!>   ifail:    (integer) failing code: 0-success
subroutine init_parallel_rngs(rngs,n_rngs,rng_size,n_tasks,my_id,ifail)
  use mpi
  use mod_rng, only: type_rng
  implicit none
  !> inputs:
  integer,intent(in) :: n_rngs,rng_size,n_tasks,my_id
  !> inputs-outputs:
  class(type_rng),dimension(:),allocatable,intent(inout) :: rngs
  integer,intent(inout) :: ifail
  !> variables
  integer :: ii,n_streams,seq
  integer,dimension(n_tasks) :: seeds
  real*8,dimension(n_tasks)  :: rands
  !> initialisations
  !> obtain different seeds from the master task and broadcast them
  if(my_id.eq.0) then
    call random_number(rands); seeds = floor(huge(0.d0)*rands)
  endif
  call MPI_Bcast(seeds,n_tasks,MPI_INTEGER,0,MPI_COMM_WORLD,ifail)
  n_streams = n_tasks*n_rngs
  !> initialise the rng
  do ii=1,n_rngs
    seq = my_id*n_rngs + ii
    call rngs(ii)%initialize(rng_size,seeds(my_id+1),n_streams,seq,ifail)
    if(ifail.ne.0) then
      write(*,'(/A)') "Error initialisation rngs: abort!"
      call MPI_abort(MPI_COMM_WORLD,-1,ifail)
    endif
  enddo
end subroutine init_parallel_rngs
!>--------------------------------------------------------------
end module mod_rng_mpi_tools
