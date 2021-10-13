!> Module containing a datatype for simulation parameters
module mod_particle_sim
use mod_particle_types
use mod_fields
use mod_openadas
use mod_coronal
use basis_at_gaussian
implicit none
private
public particle_group, particle_sim

!> A group of particles, implemented as an allocatable array.
!> It must contain particles of the same species (charge number).
type :: particle_group
  integer :: Z !< Atomic number of al particles in the group (-1 for electrons, 0 for fieldline-following)
  real*8  :: mass !< Mass of all the particles in the group
  type(ADF11_all) :: ad !< OPEN-ADAS datafiles for this species
  type(coronal) :: cor !< (coronal) equilibrium pre-calculation
  class(particle_base), dimension(:), allocatable :: particles
  real*8 :: dt !< timestep (if fixed for all particles in this group)
end type particle_group

!> Particle simulation type, containing all variables pertaining to a simulation.
type :: particle_sim
  real*8                                          :: time = 0.d0 !< time of the simulation. Only accurate when in events with sync or at
  !< the start of the simulation
  class(fields_base), allocatable                 :: fields
  logical                                         :: stop_now = .false.
  real*8                                          :: t_norm !< JOREK normalisation factor
  type(particle_group), dimension(:), allocatable :: groups
  !< MPI settings
  integer :: my_id = 0
  integer :: n_cpu = 1 ! if not initialized, act as if there is no mpi
  real*8  :: wtime_start !< Clock time at the start of the program
contains
  procedure :: finalize
  procedure :: initialize
end type particle_sim

contains
!> Actions to perform when setting up a simulation
subroutine initialize(sim, num_groups, skip_jorek2help)
  use mpi
  use mod_parameters, only: n_tor, n_period
  use phys_module, only: mode, central_mass, central_density
  use basis_at_gaussian, only: initialise_basis
  use constants, only: MU_ZERO, MASS_PROTON
  use data_structure, only: init_threads, nbthreads
  !$ use omp_lib
  class(particle_sim), intent(inout) :: sim !< why is this class() and not type()?
  integer, intent(in) :: num_groups
  logical, optional :: skip_jorek2help
  integer :: required, provided, ierr, i_tor
  character(len=MPI_MAX_PROCESSOR_NAME) :: name
  integer :: resultlength, nthreads
  logical :: my_skip_help

#ifdef FUNNELED
  required = MPI_THREAD_FUNNELED
#else
  required = MPI_THREAD_MULTIPLE
#endif

  call MPI_Init_thread(required, provided, ierr)
  if (ierr .ne. 0) write(*,*) "Error ", ierr, " in MPI_Init_thread"
  call MPI_COMM_RANK(MPI_COMM_WORLD, sim%my_id, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, sim%n_cpu, ierr)
  ! Synchronize clocks
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  sim%wtime_start = MPI_Wtime() ! accurate up to the network latency (fine for times measured in seconds)

  if (provided .ne. required .and. sim%my_id .eq. 0) write(*,*) "WARNING: provided(", provided, ") != required(", required, ")"
  allocate(sim%groups(num_groups))
  call MPI_GET_PROCESSOR_NAME(name,resultlength,ierr)
  write(*,'(A,I5,2A)') '#MPI id, ProcessorName ', sim%my_id, ': ', name
  
  call init_threads()

  if (present(skip_jorek2help)) then
    if (sim%my_id .eq. 0 .and. .not. skip_jorek2help) call jorek2help(sim%n_cpu, nbthreads)
  else
    if (sim%my_id .eq. 0) call jorek2help(sim%n_cpu, nbthreads)
  end if

  ! Initialise mode numbers
  call det_modes()

  ! Initialise parameters
  call initialise_and_broadcast_parameters(sim%my_id, "__NO_FILENAME__")

  ! Broadcast physics parameters
  call broadcast_phys(sim%my_id)

  ! Set up normalisation factors
  sim%t_norm = sqrt(MU_ZERO * central_mass * MASS_PROTON * central_density * 1.d20)

  ! Initialise the gaussian points at basis functions
  call initialise_basis
end subroutine

!> Actions to perform when stopping the simulation.
subroutine finalize(sim)
  use mod_startup_teardown, only: jorek_finalize => finalize
  class(particle_sim), intent(in) :: sim !< why is this class() and not type()?
  integer :: ierr
  if (sim%stop_now) then
    write(*,"(A,g14.6,A)") "INFO: Stop requested at ", sim%time, " , exiting"
  else
    write(*,"(A,g14.6,A)") "INFO: End of events at ", sim%time, " , exiting"
  end if
  call MPI_Finalize(ierr)
end subroutine
end module mod_particle_sim
