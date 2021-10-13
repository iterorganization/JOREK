module mod_startup_teardown
use data_structure, only: nbthreads
use mpi
implicit none
contains

!> Initialize solvers, parameters, MPI, threads etc.
subroutine initialise(my_id, n_cpu, skip_help)
  use tr_module, only: tr_meminit
  use mod_clock, only: clck_init
  use data_structure, only: init_threads
  use basis_at_gaussian
  use phys_module, only: gmres
  use mod_openadas, only : read_adf11
  use mod_impurity, only : init_imp_adas

#if ((defined WITH_Neutrals) && (!defined WITH_Impurities))
  use mod_neutral_source
#endif
#ifdef WITH_Impurities
  use mod_injection_source
#endif

#include "r3_info.h"
! Necessary for dependency reasons... should clean that up a bit and create a module
  integer, intent(out) :: my_id, n_cpu
  logical, optional, intent(in) :: skip_help
  integer :: ierr
  integer :: required, provided
  character(len=MPI_MAX_PROCESSOR_NAME) :: name
  integer :: resultlength

  interface
    subroutine set_trap_sigterm() bind(C)
    end subroutine set_trap_sigterm
  end interface

  ! --- Initialise MPI / threaded MPI
#ifdef FUNNELED
  required = MPI_THREAD_FUNNELED
#else
  required = MPI_THREAD_MULTIPLE
#endif
  call MPI_Init_thread(required, provided, ierr)
  if (ierr .ne. 0) then
      write(*,*) 'Error initializing MPI', ierr
      stop
  end if

  call init_threads()
  
  ! --- Determine number of MPI procs, ID of this proc
  call MPI_COMM_SIZE(MPI_COMM_WORLD, n_cpu, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)
  
  ! --- Process command line arguments
  if (present(skip_help)) then
    if ( my_id == 0 .and. .not. skip_help) call jorek2help(n_cpu, nbthreads)
  else
    if ( my_id == 0) call jorek2help(n_cpu, nbthreads)
  end if
  
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  call MPI_Get_processor_name(name,resultlength,ierr)
  write(*,'(A,I5,2A)') '  #MPI id, ProcessorName ', my_id, ': ', name
  call MPI_Barrier(MPI_COMM_WORLD,ierr)

  ! --- Initialise memory tracing
  call tr_meminit(my_id, n_cpu)

  ! --- Initialise timing
  call clck_init()
  call r3_info_init()
  
  ! --- Initialize mode and mode_type arrays
  call det_modes()
  
  ! --- Remove file STOP_NOW if it exists
  if ( my_id == 0 ) then
    open(42, file='STOP_NOW', iostat=ierr)
    if ( ierr == 0 ) close(42, status='delete')
  end if

  ! --- Set a signal handler for SIGTERM
  call set_trap_sigterm()

  ! --- Did MPI load correctly?
  if (required .ne. provided) then
    write(*,*) 'FATAL : MPI_THREAD_MULTIPLE (provided < required)', my_id, required, provided
    call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    stop
  end if

  ! --- MURGE with ntor=1 doesn't work up to now because i_tor is not allocated correctly
  if (n_tor == 1) then
    gmres     = .false.
  end if

#if (defined WITH_Neutrals) || (defined WITH_Impurities)
  ! --- Read ADAS data and generate coronal equilibrium if needed
  call init_imp_adas(my_id)
#endif

  
  ! --- Define the basis functions at the Gaussian points
  call initialise_basis()

  ! --- Initialise ppplib plotting library
  if (my_id == 0) call begplt('jorek2.ps')
end subroutine initialise



!> Verify that we are not doing stupid things. Run this after loading parameters
!> from the input file.
subroutine sanity_checks(my_id, n_cpu)
  use mumps_module
  use pastix_module
  use wsmp_module
  use mod_parameters, only: n_tor, n_plane
  use phys_module

  integer :: ierr
  integer, intent(in) :: my_id, n_cpu

  if ( (.not. use_mumps) .and. (.not. use_pastix) .and. (.not. use_wsmp) ) then
    write(*,*) ' FATAL : specify a valid solver'
    call MPI_Abort(MPI_COMM_WORLD, 2, ierr)
    stop
  else if ( n_plane < 2*(n_tor-1) ) then
    write(*,*) ' FATAL: n_plane >= 2 * (n_tor-1) required to avoid aliasing.'
    call MPI_Abort(MPI_COMM_WORLD, 3, ierr)
    stop
#ifndef USE_FFTW
  else if ( ( n_tor >= n_tor_fft_thresh ) .and. ( iand(n_plane,n_plane-1) /= 0 ) ) then
    write(*,*) ' FATAL: If n_tor >= n_tor_fft_thresh, n_plane must be a power of 2.'
    write(*,*) ' Hint: USE_FFTW removes this constraint.'
    call MPI_Abort(MPI_COMM_WORLD, 4, ierr)
    stop
#endif
  else if ( gmres .and. (nstep > 0) .and. (mod(n_cpu,(n_tor-1)/2+1) /= 0) ) then
    write(*,'(A,i4,A,i4,A)') ' FATAL : need a multiple of ',(n_tor-1)/2+1,' cpus for ',            &
      (n_tor-1)/2+1,' harmonics'
    call MPI_Abort(MPI_COMM_WORLD, 5, ierr)
    stop
  else if ( use_mumps ) then
#ifndef USE_MUMPS
    write(*,*) 'FATAL : use_mumps=.true. requires USE_MUMPS=1 in Makefile.inc'
    call MPI_Abort(MPI_COMM_WORLD, 6, ierr)
    stop
#endif
  else if ( use_pastix ) then
#ifndef USE_PASTIX
     write(*,*) 'FATAL : use_pastix=.true. requires USE_PASTIX=1 in Makefile.inc'
     call MPI_Abort(MPI_COMM_WORLD, 7, ierr)
     stop
#endif
  else if ( use_wsmp ) then
#ifndef USE_WSMP
    write(*,*) 'FATAL : use_wsmp=.true. requires USE_WSMP=1 in Makefile.inc'
    call MPI_Abort(MPI_COMM_WORLD, 9, ierr)
    stop
#endif
#ifdef USE_BLOCK
    write(*,*) 'FATAL : USE_BLOCK=1 in Makefile.inc is currently not possible with use_wsmp'
    call MPI_Abort(MPI_COMM_WORLD, 10, ierr)
    stop
#endif
      if ( .not. restart ) then
      write(*,*) 'FATAL : use_wsmp is currently not supported for the equilibrium'
      call MPI_Abort(MPI_COMM_WORLD, 11, ierr)
      stop
    end if
    if ( use_pastix ) then
      write(*,*) 'FATAL : you should only select one of use_wsmp or use_pastix'
      call MPI_Abort(MPI_COMM_WORLD, 12, ierr)
      stop
    end if
  end if
  if ( iand(n_plane,n_plane-1) /= 0 ) then
    write(*,*) 'WARNING: n_plane is not a power of two. This might be inefficient.'
    write(*,*) '  When using FFTW, it is possible to run like this, but it might not be fast.'
  end if
  if ( (nbthreads > 24) .and. (my_id == 0) ) then
    write(*,*) 'WARNING: You are using more than 24 OpenMP threads which might be inefficient.'
    write(*,*) '  Consider testing, whether you get better performance by increasing the number'
    write(*,*) '  of MPI tasks and reducing the number of OpenMP threads in the jobscript.'
  end if
#ifndef USE_BLOCK
  write(*,*) 'WARNING: You are not using USE_BLOCK=1 which might be inefficient.'
  write(*,*) '  Consider setting USE_BLOCK=1 in your Makefile.inc'
#endif
#ifndef USE_FFTW
  write(*,*) 'WARNING: You are not using USE_FFTW=1 which might be inefficient.'
  write(*,*) '  Consider setting USE_FFTW=1 in your Makefile.inc'
#endif
end subroutine sanity_checks


subroutine finalize(my_id)
  use phys_module, only: xtime, energies, energies2, energies3, energies4, fftw_plan
  use tr_module, only: tr_deallocate, CAT_UNKNOWN

  integer, intent(in) :: my_id
  integer :: ierr

  if (my_id == 0) then
    if (allocated(energies))  call tr_deallocate(energies,"energies",CAT_UNKNOWN)
    if (allocated(xtime))     call tr_deallocate(xtime,"xtime",CAT_UNKNOWN)

#ifdef JECCD
    if (allocated(energies2)) call tr_deallocate(energies2,"energies2",CAT_UNKNOWN)
    if (allocated(energies3)) call tr_deallocate(energies3,"energies3",CAT_UNKNOWN)
#ifdef JEC2DIAG
    if (allocated(energies4)) call tr_deallocate(energies4,"energies4",CAT_UNKNOWN)
#endif
#endif
  endif
  
#ifdef USE_FFTW
  call dfftw_destroy_plan(fftw_plan)
#endif

  call r3_info_summary()                                 ! timing
  call MPI_Finalize(ierr)                                ! clean up MPI

end subroutine finalize
end module mod_startup_teardown
