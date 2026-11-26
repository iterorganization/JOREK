module mod_startup_teardown
use data_structure, only: nbthreads
use mpi
implicit none
contains

!> Initialize solvers, parameters, MPI, threads etc.
subroutine initialise(my_id, n_mpi, skip_help)
  use tr_module, only: tr_meminit
  use mod_clock, only: clck_init
  use data_structure, only: init_threads
  use basis_at_gaussian
  use mod_openadas, only : read_adf11
  use mod_impurity, only : init_imp_adas
  use mod_plasma_functions, only: initialise_reference_parameters

#if ((defined WITH_Neutrals) && (!defined WITH_Impurities))
  use mod_neutral_source
#endif
#ifdef WITH_Impurities
  use mod_injection_source
#endif

#include "r3_info.h"
! Necessary for dependency reasons... should clean that up a bit and create a module
  integer, intent(out) :: my_id, n_mpi
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
  call MPI_COMM_SIZE(MPI_COMM_WORLD, n_mpi, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)
  
  ! --- Process command line arguments
  if (present(skip_help)) then
    if ( my_id == 0 .and. .not. skip_help) call jorek2help(n_mpi, nbthreads)
  else
    if ( my_id == 0) call jorek2help(n_mpi, nbthreads)
  end if
  
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  call MPI_Get_processor_name(name,resultlength,ierr)
  write(*,'(A,I5,2A)') '  #MPI id, ProcessorName ', my_id, ': ', name
  call MPI_Barrier(MPI_COMM_WORLD,ierr)

  ! --- Initialise memory tracing
  call tr_meminit(my_id, n_mpi)

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

  ! --- Initialize derived reference parameters
  call initialise_reference_parameters()
  
  ! --- Define the basis functions at the Gaussian points
  call initialise_basis()

  ! --- Initialise ppplib plotting library
  if (my_id == 0) call begplt('jorek2.ps')
end subroutine initialise


!> Ensure that the prescribed density is compatible with a model that does not evolve density.
subroutine check_flat_density_profile(my_id)
  use phys_module, only: LOWER_XPOINT

  integer, intent(in) :: my_id
  integer, parameter :: n_profile_points = 201
  real*8, parameter :: tolerance = 1.d-5

  integer :: i, ierr
  real*8 :: density_profile, derivatives(8), psi_n
  real*8 :: max_deviation, psi_n_max_deviation

  max_deviation = 0.d0
  psi_n_max_deviation = 0.d0

  do i = 1, n_profile_points
    psi_n = 1.5d0 * dble(i-1) / dble(n_profile_points-1)
    call density(.false., LOWER_XPOINT, 0.d0, (/0.d0, 0.d0/), psi_n, 0.d0, 1.d0, density_profile, &
                 derivatives(1), derivatives(2), derivatives(3), derivatives(4), &
                 derivatives(5), derivatives(6), derivatives(7), derivatives(8))

    if (abs(density_profile - 1.d0) .gt. max_deviation) then
      max_deviation = abs(density_profile - 1.d0)
      psi_n_max_deviation = psi_n
    endif
  enddo

  if (max_deviation .le. tolerance) return

  if (my_id .eq. 0) then
    write(*,*) 'ERROR:'
    write(*,*) '  with_rho = .false. requires density() to return a flat profile of 1.'
    write(*,*) '  Maximum deviation from 1: ', max_deviation
    write(*,*) '  At normalized flux coordinate: ', psi_n_max_deviation
  endif
  call MPI_FINALIZE(ierr)
  stop
end subroutine check_flat_density_profile



!> Verify that we are not doing stupid things. Run this after loading parameters
!> from the input file.
subroutine sanity_checks(my_id, n_mpi, mpi_required, mpi_provided)
  use mod_parameters, only: n_tor, n_plane, n_period, n_coord_period
  use phys_module
  use gauss

  integer :: ierr, i
  integer, intent(in) :: my_id, n_mpi
  integer :: nsolvers=0
  logical :: solvers(4), solvers_eq(3)
  integer :: mpi_required, mpi_provided

  ! Check for a non-flat density profile when density is not evolved
  if (.not. with_rho) call check_flat_density_profile(my_id)

  ! WARNING for axis treatment
  if(treat_axis .and. (fix_axis_nodes .or. force_central_node))then
    if (my_id .eq. 0) then
      write(*,*) 'WARNING:'
      write(*,*) '  If using treat_axis = .true. then'
      write(*,*) '  fix_axis_nodes and force_central_node both MUST be .false.'
      write(*,*) '  Setting fix_axis_nodes and force_central_node to .false.'
    endif
    force_central_node  = .false.
    fix_axis_nodes      = .false.
  endif
  if(treat_axis .and. (n_order .gt. 3))then
    if (my_id .eq. 0) then
      write(*,*) 'WARNING:'
      write(*,*) '  treat_axis = .true. is not possible'
      write(*,*) '  at the moment with n_order>3, please use fix_axis_nodes instead.'
    endif
    call MPI_FINALIZE(IERR) 
    stop
  endif

  ! WARNING for freeboundary with n_order>3
  if(freeboundary .and. (n_order .gt. 3))then
    if (my_id .eq. 0) then
      write(*,*) 'WARNING:'
      write(*,*) '  freeboundary = .true. is not possible'
      write(*,*) '  at the moment with n_order>3, aborting.'
    endif
    call MPI_FINALIZE(IERR) 
    stop
  endif

  ! WARNING for invalid values of n_order (needs to be odd)
  if( mod(n_order+1,2) .ne. 0 )then
    if (my_id .eq. 0) then
      write(*,*) 'WARNING:'
      write(*,*) '  n_order needs to be an odd integer'
    endif
    call MPI_FINALIZE(IERR) 
    stop
  endif

  ! WARNING for n_order>=7 if auto-generated mod_basisfunctions.f90 and gauss.f90 have not been created
  if( (n_order .ge. 7) .and. (n_gauss .le. 8) )then
    if (my_id .eq. 0) then
      write(*,*) 'WARNING:'
      write(*,*) '  if you are using n_order>=7 (ie. G3-Bezier or higher)'
      write(*,*) '  you need to generate the routines mod_basisfunctions.f90 and gauss.f90'
      write(*,*) '  using the code ./util/generate_codes_for_norder_gt_7.py'
      write(*,*) '  please see'
      write(*,*) '  https://www.jorek.eu/wiki/doku.php?id=gn_grid_tutorial#gn-continuous_grid_tutorial'
      write(*,*) '  for instructions.'
    endif
    call MPI_FINALIZE(IERR) 
    stop
  endif

#if (!defined(USE_PASTIX))&&(!defined(USE_PASTIX6))
  if (use_pastix .or. use_pastix_eq .or. (use_particles .and. use_pastix_prj)) then
     write(*,*) ' FATAL : use_pastix/use_pastix_eq/use_pastix_prj requires defined USE_PASTIX or USE_PASTIX6'
     call MPI_Abort(MPI_COMM_WORLD, 3, ierr)
     stop
  end if
#endif

#ifndef USE_MUMPS
  if (use_mumps.or.use_mumps_eq .or. (use_particles .and. use_mumps_prj)) then
     write(*,*) ' FATAL : use_mumps/use_mumps_eq/use_mumps_prj requires defined USE_MUMPS'
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
  if ( use_strumpack .or. use_strumpack_eq .or. (use_particles .and. use_strumpack_prj)) then
     write(*,*) ' FATAL : use_strumpack/use_strumpack_eq/use_strumpack_prj requires defined USE_STRUMPACK'
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
  else if ( (n_coord_tor < 1) .or. (mod(n_coord_tor,2) == 0) ) then
    write(*,*) 'FATAL : Hard-coded parameter n_coord_tor has an illegal value', n_coord_tor
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
  else if ( mod(n_period*n_plane, n_coord_period) /= 0 ) then
    ! The toroidal integrals are sums over n_period*n_plane equally spaced planes. Unless n_coord_period divides
    ! that number, aliasing couples toroidal harmonics of different field-period mode families.
    write(*,*) ' FATAL: n_coord_period must divide n_period*n_plane, otherwise plane aliasing couples'
    write(*,*) '        the field-period mode families.', n_coord_period, n_period, n_plane
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
  if ( (iand(n_plane,n_plane-1) /= 0)  .and. (my_id .eq. 0)) then
    write(*,*) 'WARNING: n_plane is not a power of two. This might be inefficient.'
    write(*,*) '  When using FFTW, it is possible to run like this, but it might not be fast.'
  end if
  if ( (nbthreads > 24) .and. (my_id == 0)  .and. (my_id .eq. 0)) then
    write(*,*) 'WARNING: You are using more than 24 OpenMP threads which might be inefficient.'
    write(*,*) '  Consider testing, whether you get better performance by increasing the number'
    write(*,*) '  of MPI tasks and reducing the number of OpenMP threads in the jobscript.'
  end if
  if ( ( tauIC .ne. 0.d0 ) .and. ( jorek_model == 401 )  .and. (my_id .eq. 0)) then
    write(*,*) 'WARNING: tauIC in model401 has been modified to match model303. '
    write(*,*) '         tauIC should be = m_{ion} / ( e * F0 * sqrt_mu0_rho0 * (1. + T_i/T_e) )'
  endif
  if ((abs(visco_par-visco_par_heating)/(visco_par+visco_par_heating+1.d-12) > 1.d-6) .and. (my_id .eq. 0)) then
    write(*,*) 'WARNING: The viscosity visco_par and the viscosity used for viscous heating'
    write(*,*) '  visco_par_heating are not the same. No problem if you know what you`re doing,' 
    write(*,*) '  but with this setup you are not conserving energy.   '
  endif
  if ((abs(visco-visco_heating)/(visco+visco_heating+1.d-12) > 1.d-6) .and. (my_id .eq. 0)) then
    write(*,*) 'WARNING: The viscosity visco and the viscosity used for viscous heating '
    write(*,*) '  visco_heating are not the same. No problem if you know what you`re doing,' 
    write(*,*) '  but with this setup you are not conserving energy.   '
  endif

  if ((abs(eta-eta_ohmic)/(eta+eta_ohmic+1.d-12) > 1.d-6) .and. (my_id .eq. 0)) then
    write(*,*) 'WARNING: The resistivity eta and the resistivity used for Ohmic heating '
    write(*,*) '  eta_ohm are not the same. No problem if you know what you are doing,  ' 
    write(*,*) '  but with this setup you are not conserving energy.   '
  endif
  if ((abs(T_max_eta-T_max_eta_ohm)/(T_max_eta+T_max_eta_ohm) > 1.d-6) .and. (my_id .eq. 0)) then
    write(*,*) 'WARNING: T_max_eta and T_max_eta_ohm are not the same, which breaks energy'
    write(*,*) '  conservation. No problem if you know what you are doing (a good reason  '
    write(*,*) '  to do this could be to avoid spurious Ohmic heating in the plasma core).'
  end if
  if (((T_min_neg .lt. 0.d0) .or. (rho_min_neg .lt. 0.d0)) .and. (my_id .eq. 0)) then
    write(*,*) 'WARNING: You did not specify T_min_neg and/or rho_min_neg for the correction' 
    write(*,*) '  of negative temperatures and densities. The lower values of the equilibrium'
    write(*,*) '  profiles (T_1 and/or rho_1) will be used instead.'
    write(*,*) '  For instance, try in your input file: rho_min_neg = 1.d-3' 
    write(*,*) '  and T_min_neg = 4.02d-4 !=2.01d-5*central_density*Tmin_ev' 
    write(*,*) '  (with central_density = 1 and Tmin_eV= 20 eV)' 
  endif
#ifdef WITH_Impurities
  if ((D_prof_imp_neg_thresh .gt. -1.d3) .and. (my_id .eq. 0)) then
    write(*,*) 'WARNING: You are using a value for D_prof_imp_neg_thresh that is likely to activate the correction for negative impurity density.' 
    write(*,*) '  No problem if you know what you are doing, but this could lead to convergence issues'
    write(*,*) '  in particular at the beginning of impurity injection when nimp is oscillating around zero.'
  endif
#endif

#ifndef USE_BLOCK
  if (my_id .eq. 0) then
    write(*,*) 'WARNING: You are not using USE_BLOCK=1 which might be inefficient.'
    write(*,*) '  Consider setting USE_BLOCK=1 in your Makefile.inc'
  endif
#endif
#ifndef USE_FFTW
  if (my_id .eq. 0) then
    write(*,*) 'WARNING: You are not using USE_FFTW=1 which might be inefficient.'
    write(*,*) '  Consider setting USE_FFTW=1 in your Makefile.inc'
  endif
#endif
  if (use_pastix .and. use_BLR_compression .and. (my_id .eq. 0) ) then
    write(*,*) 'WARNING: PaStiX versions before 6.x do not support BLR compression.'
    write(*,*) '  No compression will be used in this run.'
  endif
  
end subroutine sanity_checks


subroutine finalize(my_id)
  use phys_module, only: xtime, energies, energies2, energies3, energies4, fftw_plan
  use tr_module, only: tr_deallocate, CAT_UNKNOWN
  use nodes_elements
  use data_structure, only: dealloc_node_list

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

if (allocated(node_list%node)) call dealloc_node_list(node_list)
if (allocated(aux_node_list%node)) call dealloc_node_list(aux_node_list)

  call r3_info_summary()                                 ! timing
  call MPI_Finalize(ierr)                                ! clean up MPI

end subroutine finalize
end module mod_startup_teardown
