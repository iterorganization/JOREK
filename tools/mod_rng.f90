!> Module containing abstract type for multi-dimensional (Q)RNGs
!> subclass type_rng to implement other generators
module mod_rng
  implicit none
  private
  public type_rng
  public setup_shared_rngs

  type, abstract :: type_rng
    contains
      procedure (initialize), deferred :: initialize
      procedure (next),       deferred :: next
      procedure (jump_ahead), deferred :: jump_ahead
  end type

  interface
    subroutine initialize(rng, n_dims, seed, n_streams, i_stream, ierr)
      import :: type_rng
      implicit none
      class(type_rng), intent(inout) :: rng
      integer, intent(in)  :: n_dims !< Dimension of the generated output vector
      integer, intent(in)  :: seed !< Seed for the RNG if required
      integer, intent(in)  :: n_streams !< Number of output streams needed
      integer, intent(in)  :: i_stream !< Index of this output stream (1<=i_stream<=n_streams)
      integer, intent(out), optional :: ierr !< Error code. If present, return on error, otherwise call mpi_abort
    end subroutine initialize

    subroutine next(rng, out)
      import :: type_rng
      implicit none
      class(type_rng), intent(inout) :: rng
      real*8, dimension(:), intent(out) :: out
    end subroutine next

    subroutine jump_ahead(rng, delta)
      import :: type_rng
      implicit none
      class(type_rng), intent(inout) :: rng
      integer, intent(in) :: delta
    end subroutine jump_ahead
  end interface

contains


!> Setup many coordinated RNGs on MPI ranks and openmp threads in array rngs
subroutine setup_shared_rngs(n_dim, seed, rng_type, rngs)
  use mpi_mod
  !$ use omp_lib
  integer, intent(in)                                     :: n_dim !< number of dimensions of the RNG
  integer, intent(in)                                     :: seed !< seed value for the RNG, typically from [[mod_random_seed]].  Only used on mpi ID 0, discarded on others
  class(type_rng), intent(in)                             :: rng_type !< what type of rng do we want?
  class(type_rng), dimension(:), allocatable, intent(out) :: rngs !< output array of RNGs, one per openmp thread

  integer :: i, n_stream, n_streams_total
  integer, allocatable :: n_streams(:)
  integer :: my_id, n_cpu, ierr

  call MPI_Comm_Rank(MPI_COMM_WORLD, my_id, ierr)
  call MPI_Comm_Size(MPI_COMM_WORLD, n_cpu, ierr)

  call MPI_Bcast(seed, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

  n_stream = 1
  !$ n_stream = omp_get_max_threads()

  ! Get the total number of streams necessary
  ! this supports having a nonuniform number of openmp streams per mpi process. why not.
  allocate(n_streams(n_cpu))
  call MPI_AllGather(n_stream, 1, MPI_INTEGER, n_streams, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
  n_streams_total = sum(n_streams,1)

  allocate(rngs(n_stream),source=rng_type)
  do i=1,n_stream
    call rngs(i)%initialize(n_dim, seed, n_streams_total, sum(n_streams(1:my_id),1)+i)
  end do 
end subroutine setup_shared_rngs
end module mod_rng
