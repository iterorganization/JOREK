! The module mod_mpi_tools contains general pourposes mpi
! variables and procedures
module mod_mpi_tools
implicit none

private
public :: init_mpi_threads
public :: finalize_mpi_threads
public :: get_mpi_wtime

contains

! Procedures ------------------------------------------------
! the init_mpi_threads procedure initialises the mpi 
! comm-world with threads and returns the task id 
! and the total number of tasks
! outputs:
!   my_id:          (integer) rank identifier
!   n_tasks:        (integer) total number of tasks
!   ierr:           (integer) =0 : success
!   start_time_out: (real8)(optional) start time
subroutine init_mpi_threads(my_id,n_tasks,ierr,start_time_out)
  use mpi
  implicit none
  ! outputs:
  integer,intent(out)         :: my_id,n_tasks,ierr
  real*8,intent(out),optional :: start_time_out
  ! variables
  logical :: initialised
  integer :: required,provided
  real*8  :: start_time
  integer :: result_length
  character(len=MPI_MAX_PROCESSOR_NAME) :: proc_name

  ! set type of MPI thread from macro
#ifdef FUNNELED
    required = MPI_THREAD_FUNNELED
#elif SERIALIZED
    required = MPI_THREAD_SERIALIZED
#else
    required = MPI_THREAD_MULTIPLE
#endif

  ! check if mpi_init has already been called
  call MPI_Initialized(initialised,ierr)
  if(.not.initialised) call MPI_Init_thread(required,provided,ierr)
  if(ierr.ne.0) write(*,*) "Error: ",ierr," in MPI_Init_thread"
  if((my_id.eq.0).and.(provided.ne.required)) write(*,*) &
  "WARNING: provided(",provided,"different than required (",&
  required,")"
  call MPI_Comm_rank(MPI_COMM_WORLD,my_id,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,n_tasks,ierr)
  ! Synchronise all tasks
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  ! Extract wall time: accurate up to network latency 
  !(fine for time measured in seconds)
  start_time = MPI_Wtime()
  if(present(start_time_out)) start_time_out = start_time
  call MPI_Get_processor_name(proc_name,result_length,ierr)
  write(*,'(A,I5,A,2A)') "MPI ID: ",my_id," Processor name: ",proc_name
  
end subroutine init_mpi_threads

!> Finilize MPI threads and return time if required
!> inputs
!>   ierr: (integer) if !=0 an error occurred
!> outptus:
!>   ierr: (integer) if !=0 an error occurred
!>   final_time_out: (real8)(optional) process time before finalizing
subroutine finalize_mpi_threads(ierr,final_time_out)
  use mpi
  implicit none
  integer,intent(inout)       :: ierr
  real*8,intent(out),optional :: final_time_out
  logical                     :: initialised
  call MPI_Initialized(initialised,ierr)
  if(initialised) then
    if(present(final_time_out)) final_time_out = MPI_Wtime()
    call MPI_Finalize(ierr)  
  endif
end subroutine finalize_mpi_threads

!> return the mpi elapsed time on the calling task
function get_mpi_wtime()
  use mpi
  implicit none
  real*8 :: get_mpi_wtime
  get_mpi_wtime = MPI_Wtime()
end function get_mpi_wtime

!------------------------------------------------------------
end module mod_mpi_tools
