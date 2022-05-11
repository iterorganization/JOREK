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
subroutine sanity_checks(my_id, n_cpu, mpi_required, mpi_provided)
  use mumps_module
  use pastix_module
  use wsmp_module
  use mod_parameters, only: n_tor, n_plane
  use phys_module

  integer :: ierr, i
  integer, intent(in) :: my_id, n_cpu
  integer :: nsolvers=0
  logical :: solvers(4), solvers_eq(3)
  integer :: mpi_required, mpi_provided

  ! WARNING for axis treatment
  if(treat_axis .and. (fix_axis_nodes .or. force_central_node))then
    write(*,*) 'WARNING :'
    write(*,*) 'If using treat_axis = .true. then'
    write(*,*) 'fix_axis_nodes and force_central_node both MUST be .false.'
    write(*,*) 'Setting fix_axis_nodes and force_central_node to .false.'    
    force_central_node  = .false.
    fix_axis_nodes      = .false.
  endif
  if(treat_axis .and. (n_order .gt. 3))then
    write(*,*) 'WARNING :'
    write(*,*) 'treat_axis = .true. is not possible'
    write(*,*) 'at the moment with n_order>3, please use fix_axis_nodes instead.'
    call MPI_FINALIZE(IERR) 
    stop
  endif

  ! WARNING for freeboundary with n_order>3
  if(freeboundary .and. (n_order .gt. 3))then
    write(*,*) 'WARNING :'
    write(*,*) 'freeboundary = .true. is not possible'
    write(*,*) 'at the moment with n_order>3, aborting.'
    call MPI_FINALIZE(IERR) 
    stop
  endif

#if (!defined(USE_PASTIX))&&(!defined(USE_PASTIX6))
  if (use_pastix.or.use_pastix_eq) then
    write(*,*) ' FATAL : use_pastix requires defined USE_PASTIX or USE_PASTIX6'
    call MPI_Abort(MPI_COMM_WORLD, 3, ierr)
    stop
  endif
#endif

#ifndef USE_MUMPS
  if (use_mumps.or.use_mumps_eq) then
    write(*,*) ' FATAL : use_mumps requires defined USE_MUMPS'
    call MPI_Abort(MPI_COMM_WORLD, 3, ierr)
    stop
  endif
#endif

#ifndef USE_WSMP
  if (use_wsmp) then
    write(*,*) ' FATAL : use_wsmp requires defined USE_WSMP'
    call MPI_Abort(MPI_COMM_WORLD, 3, ierr)
    stop
  endif
#endif

#ifndef USE_STRUMPACK
  if (use_strumpack.or.use_strumpack_eq) then
    write(*,*) ' FATAL : use_strumpack requires defined USE_STRUMPACK'
    call MPI_Abort(MPI_COMM_WORLD, 3, ierr)
    stop
  endif
#endif
  
  ! --- Check solver consistency
  solvers = (/use_mumps,use_pastix,use_wsmp,use_strumpack/)
  nsolvers = 0
  do i=1,size(solvers)
    if (solvers(i)) nsolvers = nsolvers + 1
  enddo
  if (nsolvers.eq.0) then
    write(*,*) ' FATAL : specify a valid solver'
    call MPI_Abort(MPI_COMM_WORLD, 3, ierr)
    stop
  elseif (nsolvers.gt.1) then
    write(*,*) ' FATAL : specify only one solver'
    call MPI_Abort(MPI_COMM_WORLD, 3, ierr)
    stop
  endif
  
  solvers_eq = (/use_mumps_eq,use_pastix_eq,use_strumpack_eq/)
  nsolvers = 0
  do i=1,size(solvers_eq)
    if (solvers_eq(i)) nsolvers = nsolvers + 1
  enddo
  if (nsolvers.eq.0) then
    write(*,*) ' FATAL : specify a valid equilibrium solver'
    call MPI_Abort(MPI_COMM_WORLD, 3, ierr)
    stop
  elseif (nsolvers.gt.1) then
    write(*,*) ' FATAL : specify only one equilibrium solver'
    call MPI_Abort(MPI_COMM_WORLD, 3, ierr)
    stop
  endif

  ! --- Some checks not to waste any cpu time
  if ( (n_tor < 1) .or. (mod(n_tor,2) == 0) ) then
    write(*,*) 'FATAL : Hard-coded parameter n_tor has an illegal value', n_tor
    call MPI_Abort(MPI_COMM_WORLD, 23, ierr)
    stop
  else if ( (n_coord_tor > 1) .and. (rst_hdf5_version .eq. 1) ) then
    write(*,*) 'FATAL : Hard-coded parameter n_coord_tor > 1 is only possible with rst_hdf5_version > 1'
    call MPI_Abort(MPI_COMM_WORLD, 23, ierr)
  else if ( n_period<1 ) then
    write(*,*) 'FATAL : Hard-coded parameter n_period has an illegal value', n_period
    call MPI_Abort(MPI_COMM_WORLD, 24, ierr)
    stop
  else if ( n_elements_max<1 ) then
    write(*,*) 'FATAL : Hard-coded parameter n_elements_max has an illegal value', n_elements_max
    call MPI_Abort(MPI_COMM_WORLD, 25, ierr)
    stop
  else if ( n_nodes_max<1 ) then
    write(*,*) 'FATAL : Hard-coded parameter n_nodes_max has an illegal value', n_nodes_max
    call MPI_Abort(MPI_COMM_WORLD, 25, ierr)
    stop
  else if ( n_boundary_max<1 ) then
    write(*,*) 'FATAL : Hard-coded parameter n_boundary_max has an illegal value', n_boundary_max
    call MPI_Abort(MPI_COMM_WORLD, 25, ierr)
    stop
  else if ( n_pieces_max<1 ) then
    write(*,*) 'FATAL : Hard-coded parameter n_pieces_max has an illegal value', n_pieces_max
    call MPI_Abort(MPI_COMM_WORLD, 25, ierr)
    stop
  else if ( n_vertex_max/=4 ) then
    write(*,*) 'WARNING : hard-coded parameter n_vertex_max /= 4', n_vertex_max
    call MPI_Abort(MPI_COMM_WORLD, 25, ierr)
    stop
  else if (mpi_required .ne. mpi_provided) then
    write(*,*) 'FATAL : MPI_THREAD_MULTIPLE (provided < required)', my_id, mpi_required, mpi_provided
    call MPI_Abort(MPI_COMM_WORLD, 2, ierr)
    stop
  else if ( mod(n_tor,2) == 0 ) then
    write(*,*) ' FATAL: n_tor must be an uneven number.'
    call MPI_Abort(MPI_COMM_WORLD, 4, ierr)
    stop
  else if ( n_plane < 2*(n_tor-1) ) then
    write(*,*) ' FATAL: n_plane >= 2 * (n_tor-1) required to avoid aliasing.'
    call MPI_Abort(MPI_COMM_WORLD, 4, ierr)
    stop
#ifndef USE_FFTW
  else if ( ( n_tor >= n_tor_fft_thresh ) .and. ( iand(n_plane,n_plane-1) /= 0 ) ) then
    write(*,*) ' FATAL: If n_tor >= n_tor_fft_thresh, n_plane must be a power of 2.'
    write(*,*) ' Hint: USE_FFTW removes this constraint.'
    call MPI_Abort(MPI_COMM_WORLD, 5, ierr)
    stop
#endif
  else if ( n_tor_fft_thresh < 2 ) then
    write(*,*) ' FATAL: n_tor_fft_thresh < 2 presently not allowed. Will cause problems for n_tor=1.'
    call MPI_Abort(MPI_COMM_WORLD, 5, ierr)
    stop
  else if ( use_wsmp ) then
#ifdef USE_BLOCK
    write(*,*) 'FATAL : USE_BLOCK=1 in Makefile.inc is currently not possible with use_wsmp'
    call MPI_Abort(MPI_COMM_WORLD, 11, ierr)
    stop
#endif
    if ( .not. restart ) then
      write(*,*) 'FATAL : use_wsmp is currently not supported for the equilibrium'
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
  if ( ( tauIC .ne. 0.d0 ) .and. ( jorek_model == 401 ) ) then
    write(*,*) 'WARNING: tauIC in model401 has been modified to match model303. '
    write(*,*) '         tauIC should be = m_{ion} / ( e * F0 * sqrt_mu0_rho0 * (1. + T_i/T_e) )'
  endif
  if (abs(visco_par-visco_par_heating)/(visco_par+visco_par_heating+1.d-12) > 1.d-6) then
    write(*,*) 'WARNING: The viscosity visco_par and the viscosity used for viscous heating '
    write(*,*) '  visco_par_heating are not the same. No problem if you know what you are doing,  ' 
    write(*,*) '  but with this setup you are not conserving energy.   '
  endif

  if (abs(eta-eta_ohmic)/(eta+eta_ohmic+1.d-12) > 1.d-6) then
    write(*,*) 'WARNING: The resistivity eta and the resistivity used for Ohmic heating '
    write(*,*) '  eta_ohm are not the same. No problem if you know what you are doing,  ' 
    write(*,*) '  but with this setup you are not conserving energy.   '
  endif
  if (abs(T_max_eta-T_max_eta_ohm)/(T_max_eta+T_max_eta_ohm) > 1.d-6) then
    write(*,*) 'WARNING: T_max_eta and T_max_eta_ohm are not the same, which breaks  &
        energy conservation. No problem if you know what you are doing (a good reason to &
	do this could be to avoid spurious Ohmic heating in the plasma core).'
  end if
  if ((T_min_neg .lt. 0.d0) .or. (rho_min_neg .lt. 0.d0)) then
	write(*,*) 'WARNING: You did not specify T_min_neg and/or rho_min_neg for the correction of negative temperatures and densities.  & 
	   The lower values of the equilibrium profiles (T_1 and/or rho_1) will be used instead.'
	write(*,*) 'For instance, try in your input file: rho_min_neg = 1.d-3 and T_min_neg = 4.02d-4 !=2.01d-5*central_density*Tmin_ev (with central_density = 1 and Tmin_eV= 20 eV)'    
  endif

#ifndef USE_BLOCK
  write(*,*) 'WARNING: You are not using USE_BLOCK=1 which might be inefficient.'
  write(*,*) '  Consider setting USE_BLOCK=1 in your Makefile.inc'
#endif
#ifndef USE_FFTW
  write(*,*) 'WARNING: You are not using USE_FFTW=1 which might be inefficient.'
  write(*,*) '  Consider setting USE_FFTW=1 in your Makefile.inc'
#endif
  if (use_pastix .and. use_BLR_compression) then
    write(*,*) 'WARNING: PaStiX versions before 6.x do not support BLR compression.'
    write(*,*) '  No compression will be used in this run.'
  endif
  
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
