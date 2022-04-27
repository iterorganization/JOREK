!> JOREK 2.0 -- Solves the (reduced) MHD equations in 3D toroidal geometry.
!!
!! By using this code, you accept the code license available at:
!! https://www.jorek.eu/JOREK-license.pdf
!!
!! - solvers implemented:
!!   - MUMPS
!!   - PastiX
!!   - STRUMPACK
!!   - GMRES (+MUMPS, PastiX or STRUMPACK preconditioner)
!!
!! - required libraries :
!!   - MPI
!!   - MUMPS
!!   - PastiX
!!   - SCOTCH (metis)
!!   - FFTW
!!   - SCALAPACK (BLACS)
!!   - LAPACK, BLAS
!!   - PPPLIB
!!
!! @author Guido Huysmans (Euratom / CEA Association) and the whole JOREK Team
!! @date 18-7-2008
program JOREK2

  use constants
  use mumps_module
  use pastix_module
  use wsmp_module
  use data_structure
  use phys_module
  use mod_parameters
  use mod_log_params
  use global_distributed_matrix
  use nodes_elements
  use pellet_module
  use equil_info
  use mod_boundary,            only: boundary_from_grid
  use vacuum
  use vacuum_response,     only: get_vacuum_response, update_response, init_wall_currents, I_coils
  use vacuum_equilibrium,  only: import_external_fields
  use live_data
  use mod_bootstrap_functions
  use construct_matrix_mod, only : construct_matrix
  use mod_global_matrix_structure
  use mod_import_restart
  use mod_export_restart
  use mod_element_rtree, only: populate_element_rtree
  use mod_interp
  use basis_at_gaussian, only: initialise_basis
  use mod_expression, only: exprs_all_int, init_expr
  use mod_integrals3D
  use mod_openadas, only : read_adf11
  use mod_atomic_coeff_deuterium, only: ad_deuterium 
#ifdef USE_STRUMPACK
  use strumpack_module
#endif
#ifdef USE_PASTIX6
  use mod_pastix, only: pastix_finalize
#endif
  use preconditioner_module
  use mod_distribute_preconditioner
  use direct_construction_mod
  use centralization_mod
  use mod_exchange_indices
  use mod_gmres, only: gmres_driver
  use mod_startup_teardown
  use mod_initial_grid
  use mod_flux_grid

! these write additional live data (global data) used when an ECCD current is applied)
#ifdef JECCD
  use live_data2,          only: init_live_data2, write_live_data2, finalize_live_data2
  use live_data3,          only: init_live_data3, write_live_data3, finalize_live_data3
#ifdef JEC2DIAG
  use live_data4,          only: init_live_data4, write_live_data4, finalize_live_data4
#endif
#endif

  use solve_mat_n
  use tr_module
  use mod_clock
#ifdef USE_HDF5
  use hdf5
  use hdf5_io_module
#endif
  use mpi_mod
  use mod_impurity, only: init_imp_adas
#ifdef USE_BICGSTAB
  use mod_bicgstab, only: bicgstab_driver, bicgstab_finalize
#endif  

  use, intrinsic :: iso_c_binding
  use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
                                            stdout=>output_unit, &
                                            stderr=>error_unit
  
  implicit none

#ifdef USE_FFTW
  include 'fftw3.f03'
#endif
  
#include "r3_info.h"
  
  interface
    subroutine equilibrium(my_id,node_list,element_list,bnd_node_list,bnd_elm_list,xpoint2,xcase2, nice_q)
      use data_structure
      integer(kind=4),             intent(in)    :: my_id
      integer(kind=4),             intent(in)    :: xcase2
      type (type_node_list),       intent(inout) :: node_list
      type (type_element_list),    intent(inout) :: element_list
      type (type_bnd_node_list)   ,intent(inout) :: bnd_node_list    
      type (type_bnd_element_list),intent(inout) :: bnd_elm_list    
      logical(kind=4),             intent(in)    :: xpoint2
      logical(kind=4),             intent(in)    :: nice_q
    end subroutine equilibrium

    subroutine set_trap_sigterm() bind(C)
    end subroutine set_trap_sigterm
    logical function sigterm_called() bind(C)
    end function sigterm_called
  end interface
  
  type (type_surface_list) :: surface_list
  real*8                   :: W_mag(n_tor), W_kin(n_tor), growth_mag, growth_kin, growth_mag0, growth_kin0
#ifdef JECCD
  real*8                   :: A_tem(n_tor), A_den(n_tor), A_jen(n_tor), A_jec(n_tor),A_jec1(n_tor), A_jec2(n_tor)
#endif
  real*8                   :: t_matrix, t_send, t_solve
  type(clcktype)           :: t_itstart, t0, t1
  real*8                   :: mindelta, maxdelta
  integer                  :: my_id, my_id_n, my_id_master
  integer                  :: istep,jstep,ierr,i,itor,inode, i_elm_axis, i_elm_xpoint(2)
  integer                  :: n_local_ELMs
  integer                  :: i_rank(n_tor), n_cpu, n_cpu_n, n_cpu_master, m_cpu, n_masters, n_cpu_trans, my_id_trans
  integer                  :: iter_gmres
  integer                  :: MPI_COMM_N, MPI_GROUP_MASTER, MPI_GROUP_WORLD, MPI_COMM_MASTER, MPI_COMM_TRANS
  character*8              :: label, itlabel
  character*14             :: fileout
  integer                  :: mpi_required,mpi_provided,StatInfo
  integer, allocatable     :: local_elms(:), index_min(:), index_max(:)
  real*8                   :: zjz, E_min, E_max
  logical                  :: solve_only, to_quit, freeb_equil2
  integer*4                :: rank, comm_size 
  real*8                   :: zn,  dn_dpsi,  dn_dz,  dn_dpsi2,  dn_dz2,  dn_dpsi_dz,  dn_dpsi3,  dn_dpsi_dz2,  dn_dpsi2_dz
  real*8                   :: zT,  dT_dpsi,  dT_dz,  dT_dpsi2,  dT_dz2,  dT_dpsi_dz,  dT_dpsi3,  dT_dpsi_dz2,  dT_dpsi2_dz
  real*8                   :: zTi, dTi_dpsi, dTi_dz, dTi_dpsi2, dTi_dz2, dTi_dpsi_dz, dTi_dpsi3, dTi_dpsi_dz2, dTi_dpsi2_dz
  real*8                   :: zTe, dTe_dpsi, dTe_dz, dTe_dpsi2, dTe_dz2, dTe_dpsi_dz, dTe_dpsi3, dTe_dpsi_dz2, dTe_dpsi2_dz
  real*8                   :: zFFprime, dFFprime_dpsi, dFFprime_dz, dFFprime_dpsi_dz,dFFprime_dpsi2,dFFprime_dz2
  real*8                   :: Rp, Zp, R_out,Z_out,s_out,t_out,P_s,P_t,P_st,P_ss,P_tt, psi
  real*8                   :: Rp_start, Rp_end, density_tot,density_in,density_out,pressure_tot,pressure_in,pressure_out,Bgeo
  real*8,allocatable       :: xp(:), yp1(:), yp2(:), yp3(:)
  real*8,allocatable       :: res(:) 
  integer                  :: nplot, iplot, i_elm, ifail, ivar, iter_big, n_aa, iter_prev, n_since_update, n_spi_begin
  logical                  :: is_local, file_exists
  integer                  :: i_elem, inode1, i_order, index_node1
  type (type_element)      :: element
  integer                  :: index_size, id_elements
  integer                  :: list_to_be_refined(n_ref_list), n_to_be_refined    
  REAL*8                   :: max_time, min_time, tsecond
  integer, allocatable     :: tab_n_local_elems(:)
  real*8                   :: t_this, sum_deltas
! =================== plot NEO coeffs ==================
  real*8                   :: amu_neo_node, aki_neo_node
  real*8,allocatable       :: mu_neo(:), ki_neo(:)
! ======================================================
#ifdef USE_FFTW
  real*8     :: in_fft(1:n_plane)
  complex*16 :: out_fft(1:n_plane)
#endif

  real*8  :: DUMMY_REAL(1:1)
  integer :: DUMMY_INT (1:1)
  character(len=MPI_MAX_PROCESSOR_NAME) :: name
  integer :: resultlength

  integer :: holder
  integer :: getpid

  logical :: input_treat_axis
  
  call init_expr()
  allocate(res(exprs_all_int%n_expr+1))
  res = 0.d0   
  
  !***********************************************************************
  !*                  intialisation                                      *
  !***********************************************************************
  ! --- Initialize OpenMP threads before MPI_init
  !call init_threads()
  
  ! --- Initialise MPI / threaded MPI
#ifdef FUNNELED
  mpi_required = MPI_THREAD_FUNNELED
#else
  mpi_required = MPI_THREAD_MULTIPLE
#endif
#ifdef STAN_FLAG
mpi_required = 0
#endif
  call MPI_Init_thread(mpi_required, mpi_provided, StatInfo)

  call init_threads()  ! on some systems init_threads needs to come after mpi_init_thread
  
  ! --- Determine number of MPI procs
  call MPI_COMM_SIZE(MPI_COMM_WORLD, comm_size, ierr)
  n_cpu = comm_size
  
  ! --- Determine ID of each MPI proc
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  my_id = rank
  
  ! --- Process command line arguments
  if ( my_id == 0 ) call jorek2help(n_cpu, nbthreads)
  
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  CALL MPI_GET_PROCESSOR_NAME (name,resultlength,ierr)
  write(*,'(A,I5,2A)') '  #MPI id, ProcessorName ', rank, ': ', name
  call MPI_Barrier(MPI_COMM_WORLD,ierr)

  ! --- Initialise memory tracing
  call tr_meminit(my_id, n_cpu)

  ! --- Initialise timing
  call clck_init()
  call r3_info_init ()
  
  ! --- Initialize mode and mode_type arrays
  call det_modes()
  
  ! --- Remove file STOP_NOW if it exists
  if ( my_id == 0 ) then
    open(42, file='STOP_NOW', iostat=ierr)
    if ( ierr == 0 ) close(42, status='delete')
  end if

  ! --- Set a signal handler for SIGTERM
  call set_trap_sigterm()

  ! --- Preset some solver variables
  pastix_initialised = .false.
  pastix_analysed    = .false.
#ifdef USE_STRUMPACK  
  spss_initialized = .false.
  spss_analyzed    = .false.
#endif  
  
  ! --- Preset input parameters to reasonable defaults, then read the input file.
  call initialise_and_broadcast_parameters(my_id, "__NO_FILENAME__")
  
  ! --- Initialize the vacuum part.
  call vacuum_init(my_id, freeboundary_equil, freeboundary, resistive_wall)

  if (nstep .gt. 0)   call check_preconditioner_consistency
  
  ! --- Initialize live data file which will be filled during the code run
  if ( my_id == 0 ) call init_live_data()
#ifdef JECCD
  if ( my_id == 0 ) call init_live_data2()
  if ( my_id == 0 ) call init_live_data3()
#ifdef JEC2DIAG
  if ( my_id == 0 ) call init_live_data4()
#endif
#endif
  
  ! --- Initialise ppplib plotting library
  if (my_id == 0 .and. write_ps)  call begplt('jorek2.ps')
  
  ! --- Define the basis functions at the Gaussian points
  call initialise_basis()
  
#if (defined WITH_Neutrals) || (defined WITH_Impurities)
  ! --- Read ADAS data and generate coronal equilibrium if needed
  call init_imp_adas(my_id)
#else
  if (use_imp_adas .and. (nimp_bg(1) > 0.d0)) then
    call init_imp_adas(my_id)
  endif
#endif

  ! --- Write out all parameters defined in parameters and the namelist input file.
  call log_parameters(my_id)
 
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  
  call sanity_checks(my_id, n_cpu, mpi_required, mpi_provided)

  call tr_print_memsize("InitStep")

  !***********************************************************************
  !*                  read restart file                                  *
  !***********************************************************************
  
  input_treat_axis = treat_axis   ! store the value from the input file

  if ( restart .and. (my_id == 0) ) then
    
    call import_restart(node_list, element_list, 'jorek_restart', rst_format, ierr)
    if ( ierr /= 0 ) stop

    ! for variable time step Gears method
    if ( index_now <= 1 ) then
      tstep_prev = tstep
    else
      tstep_prev = xtime(index_start) - xtime(index_start-1)
    end if

    ! check consistency for axis treatment with restart file
    if (input_treat_axis .neqv. treat_axis) then
      write(*,*) 'WARNING: Axis treatments set via input file is not the same as that in the restart file.'
      write(*,*) 'You are trying to restart the simulation with treat_axis = ', input_treat_axis
      write(*,*) 'Earlier treat_axis was set to = ', treat_axis
      write(*,*) 'STOP' 
      stop      
    endif

    ! --- Write live data for previous time-steps
    if ( .not. bench_without_plot ) then
      do index_now = 1, index_start
        call write_live_data(index_now)
        call write_live_data_vacuum(index_now)
#ifdef JECCD
        call write_live_data2(index_now)
        call write_live_data3(index_now)
#ifdef JEC2DIAG
        call write_live_data4(index_now)
#endif
#endif
      end do
    end if
    
    ! --- Optional: Redo flux aligned grid (DOES NOT WORK CURRENTLY)
    if (regrid) then
      if (xpoint)  then
        if ( (xcase .ge. UPPER_XPOINT) .or. (RZ_grid_inside_wall) ) then
          if (grid_to_wall) then
            call grid_double_xpoint_inside_wall(node_list, element_list)
          else
            call grid_double_xpoint(node_list, element_list)
          endif
        else
	  call grid_xpoint(node_list,element_list,n_flux,n_open,n_private,n_leg,n_tht,  &
                           SIG_open,SIG_closed,SIG_private,SIG_theta,SIG_leg_0,SIG_leg_1,dPSI_open,dPSI_private, xcase)
        endif
      else
        call grid_flux_surface(xpoint,xcase, node_list, element_list, surface_list, n_flux, n_tht, xr1,  &
                               sig1, xr2, sig2, refinement)
      end if
      
    end if
    
    if ( freeboundary .and. freeb_change_indices) call exchange_indices(node_list, my_id, n_cpu, .false.)
    
 end if !   if ( restart .and. (my_id == 0) ) then

  ! This is necessary for the parallel vacuum version during the code restart 
  if(restart) then
    call broadcast_phys(my_id)  
    if(freeboundary) call broadcast_vacuum(my_id, resistive_wall)
  end if
  call populate_element_rtree(node_list, element_list)
  
  !***********************************************************************
  !*                  define grid / equilibrium                          *
  !***********************************************************************
  
  if_not_restart: if (.not. restart) then
    call tr_resetfile()
    
    call initial_grid(node_list, element_list, bnd_node_list, bnd_elm_list, my_id, n_cpu)
    
    ! --- Synchronizing MPI processes avoid deadlock issues on some machine
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    
    ! --- Send boundary elements and nodes to other MPI procs
    call broadcast_boundary(my_id,bnd_elm_list,bnd_node_list)
    
    ! --- Fill the vacuum response matrices for freeboundary computations
    if ( freeboundary_equil .and. (n_flux .eq. 0)) then
      call get_vacuum_response(my_id, node_list, bnd_elm_list, bnd_node_list, freeboundary_equil,  &
        resistive_wall)
      call update_response(my_id,tstep, freeboundary_equil, resistive_wall)
      call import_external_fields('coil_field.dat', my_id)
      call set_coil_curr_time_trace()
      if ( (.not. restart) .or. (.not. wall_curr_initialized) ) call init_wall_currents(my_id, resistive_wall)
    else
      freeb_equil2        = freeboundary_equil
      freeboundary_equil  = .false.
    end if
    
    ! --- Plot the grid  
    if ( (my_id == 0) .and. (.not. bench_without_plot) ) then
      call plot_grid(node_list,element_list,bnd_elm_list,bnd_node_list,.true.,.false.,'initial')
    end if
    
    ! --- Check sanity of grid
    if (.not. RZ_grid_inside_wall) call check_grid(my_id, node_list, element_list)

#ifdef USE_MUMPS
    ! --- Initialize MUMPS solver (used for equilibrium)
    call MPI_COMM_GROUP(MPI_COMM_WORLD,MPI_GROUP_WORLD,ierr)
    call MPI_GROUP_INCL(MPI_GROUP_WORLD,1,[0],MPI_GROUP_MUMPS_EQUIL,ierr)
    call MPI_COMM_CREATE(MPI_COMM_WORLD,MPI_GROUP_MUMPS_EQUIL,MPI_COMM_MUMPS_EQUIL,ierr)
    if (my_id == 0) call initialise_mumps(MPI_COMM_MUMPS_EQUIL)
#endif

    ! --- Compute the plasma equilibrium
    if (equil) then
      call equilibrium(my_id,node_list,element_list,bnd_node_list,bnd_elm_list,xpoint,xcase, .true.) 
      if (export_for_nemec) then
        if(my_id ==0 ) call export_nemec(node_list, element_list, xpoint, xcase)
      endif
      if (my_id == 0) call update_equil_state(my_id,node_list, element_list, bnd_elm_list, xpoint, xcase)
    end if ! if (equil) then

    ! --- Determine a flux surface aligned grid and re-calculate the equilibrium on it
    if (n_flux > 1) then

      call flux_grid(node_list, element_list, bnd_node_list, bnd_elm_list, my_id, n_cpu)
      
      if ( freeb_equil2) then
        freeboundary_equil = .true.
        call get_vacuum_response(my_id, node_list, bnd_elm_list, bnd_node_list, freeboundary_equil,  &
          resistive_wall)
        call update_response(my_id,tstep, freeboundary_equil, resistive_wall)
        call import_external_fields('coil_field.dat', my_id)
        call set_coil_curr_time_trace()
        if ( (.not. restart) .or. (.not. wall_curr_initialized) ) call init_wall_currents(my_id, resistive_wall)
      end if
      
      ! --- Compute the plasma equilibrium
      call equilibrium(my_id, node_list, element_list, bnd_node_list, bnd_elm_list, xpoint,xcase, .false.)

    end if ! if (n_flux > 1) then
 
    if (my_id == 0) then
          
      ! --- Update the status of the equilibrium
      call update_equil_state(my_id,node_list, element_list, bnd_elm_list, xpoint, xcase)
      
      ! --- Set initial conditions for time-evolution
      call initial_conditions(my_id,node_list,element_list,bnd_node_list, bnd_elm_list, xpoint,xcase)

!      call remove_centre(node_list,element_list,n_tht,67*(n_tht-1))

      ! --- Determine initial energies
      call energy(W_mag,W_kin)
      write(*,'(A,12e16.8)') ' initial energies : ', W_mag, W_kin

    end if ! (my_id == 0)
    
#ifdef USE_MUMPS
    ! --- Clean up this instance of mumps (used for equilibrium)
    mumps_par%JOB = -2
    if (my_id == 0) call DMUMPS(mumps_par)
#endif
    ! -- For PaStiX solver before version 6.x
    if (allocated(pastix_perm_vars))  call tr_deallocate(pastix_perm_vars,"pastix_perm_vars",CAT_UNKNOWN)
    if (allocated(pastix_iperm_vars)) call tr_deallocate(pastix_iperm_vars,"pastix_iperm_vars",CAT_UNKNOWN)

  end if if_not_restart
  
  ! --- Print some grid information
  if ( my_id == 0 ) call log_grid_info(.false., node_list, element_list)
  
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  
  ! --- Determine boundary information from the grid
  if ( my_id == 0 ) call boundary_from_grid(node_list, element_list, bnd_node_list, bnd_elm_list, output_bnd_elements)
  call broadcast_boundary(my_id, bnd_elm_list, bnd_node_list)
  
  ! --- Fill the vacuum response matrices for freeboundary computations
  if ( freeboundary ) then
    call get_vacuum_response(my_id, node_list, bnd_elm_list, bnd_node_list, freeboundary_equil,    &
      resistive_wall)
    call update_response(my_id,tstep, freeboundary_equil, resistive_wall)
    call import_external_fields('coil_field.dat', my_id)
    call set_coil_curr_time_trace()
    call read_Z_axis_profile() 
    if ( (.not. restart) .or. (.not. wall_curr_initialized) ) call init_wall_currents(my_id, resistive_wall)
  end if
  
  call tr_print_memsize("AfterEquilibrium")

  if (RMP_on) then
     if (my_id == 0) then
        call read_RMP_profiles(bnd_node_list)
     endif

  endif
  
  ! --- Broadcast grid information and input parameters to other MPI procs
  call broadcast_elements(my_id, element_list)                ! elements
  if (RMP_on) then
     call broadcast_RMP_profiles(my_id, bnd_node_list)        ! psi_RMP profiles
  endif

  call broadcast_nodes(my_id, node_list)                      ! nodes

  ! Let every mpi proc calculate this
  call populate_element_rtree(node_list, element_list)

  call broadcast_phys(my_id)                                  ! physics parameters

  ! --- Broadcast equil_state: This is needed because find_axis depends on the axis
  ! --- from the previous time-step, which is only read by my_id=0 from the restart file
  call broadcast_equil_state(my_id)                           ! equil_state

  if ( freeboundary ) call broadcast_vacuum(my_id, resistive_wall)
  n_AA = 0  
  do inode = 1, node_list%n_nodes  
    n_AA = max(n_AA,node_list%node(inode)%index(4))  
  end do
  mumps_par%n = n_AA

  ! --- Load deuterium ADAS data if required
  if (deuterium_adas) ad_deuterium =  read_adf11(my_id,'96_h') 

   ! --- Initialize FFTW
#ifdef USE_FFTW
  call dfftw_plan_dft_r2c_1d(fftw_plan,n_plane,in_fft,out_fft,FFTW_PATIENT)
#endif
 
  call tr_debug_write("JMAIN:End_init elt_list",element_list%n_elements)
  call tr_debug_write("JMAIN:End_init bnd_elt_list",bnd_elm_list%n_bnd_elements)
  call tr_debug_write("JMAIN:End_init node_list",node_list%n_nodes)
  call tr_debug_write("JMAIN:End_init nAA",n_AA)

  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  
  call update_equil_state(my_id,node_list, element_list, bnd_elm_list, xpoint, xcase)
  if ( my_id == 0 ) then
    call print_equil_state(.true.)
    call save_special_points('special_equilibrium_points.dat', .false., ierr)
  end if

  !***********************************************************************
  !*                 end of initilisation/equilibrium                    *
  !***********************************************************************
  
  t_now     = t_start      ! t_now: current time in the simulation
  
  if (nstep > 0) then
    
    !*******************************************************
    !*      create groups /communicators		   *
    !* MPI_COMM_N      : group for each harmonic	   *
    !* MPI_COMM_TRANS  : Transversal communicator	   *
    !*   (ie : all first proc of MPI_COMM_N, all second,   *
    !*         all third...)				   *
    !* MPI_COMM_MASTER : group of masters of each harmonic *
    !*  		 (i.e id=0 from each MPI_COMM_N)   *
    !*******************************************************
    if (gmres) then
    
       call create_communicators(my_id_n, n_cpu_n, MPI_COMM_N, my_id_master, n_masters, &
                                 MPI_COMM_MASTER, MPI_COMM_TRANS)
       m_cpu = n_cpu_n
       write(*,*) "my_id, my_id_n", my_id, my_id_n

       call distribute_modes
       
    else
       my_id_n = my_id
       MPI_COMM_N = MPI_COMM_WORLD
       m_cpu = n_cpu
    endif


    !***********************************************************************
    !*  	  distribute nodes and elements over cpu's		   *
    !***********************************************************************
    index_size  = n_cpu
    id_elements = my_id

    call tr_allocate(local_elms,1,element_list%n_elements,"local_elms",CAT_FEM)
    call tr_allocate(index_min,1,index_size,"index_min",CAT_FEM)
    call tr_allocate(index_max,1,index_size,"index_max",CAT_FEM)
    call tr_allocate(local_index_start,1,n_cpu,"local_index_start",CAT_FEM)
    call tr_allocate(local_index_end,1,n_cpu,"local_index_end",CAT_FEM)

    !
    ! Construct index_min, index_max and local_elems
    !
    call distribute_nodes_elements(id_elements,m_cpu,index_size,node_list,element_list,.false.,local_elms, & 
         n_local_elms,ndof_glob,index_min,index_max, restart, freeboundary)    

    node_list%n_dof = ndof_glob
    local_index_start = index_min
    local_index_end   = index_max
    ! Build ijA_index, ijA_size and irn_jcn

    call global_matrix_structure(my_id,my_id_n,node_List,element_list,bnd_elm_list, freeboundary,&
         local_elms,n_local_elms,index_min(id_elements+1),index_max(id_elements+1),              & 
         ijA_index, ijA_size, irn_jcn, irn_glob, jcn_glob, 1, n_tor,                 &
         n_glob, nz_glob, ndof_glob, n_matrix_block_size)

    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    if ( freeboundary .and. ( sr%n_tor /= 0 ) ) then 
      call global_matrix_structure_vacuum(node_list, bnd_node_list, index_min(my_id+1), index_max(my_id+1),& 
           1, n_tor, irn_glob, jcn_glob, n_matrix_block_size, ijA_index, ijA_size, irn_jcn) 
    endif

    if ((gmres).and.(my_id_n.eq.0)) call map_row_index(ndof_glob)
    if (use_mumps) then
       if (.not. gmres) then
         call initialise_mumps(MPI_COMM_WORLD)    ! start MUMPS sparse matrix solver all cpus
       else
         call initialise_mumps(MPI_COMM_N)        ! start MUMPS sparse matrix solver on local groups
       endif
    endif

  endif ! (nstep >0)
  
  ! --- Export a restart file before the first timestep
  if ( (my_id == 0) .and. (.not. restart) ) then
    if ( freeboundary .and. freeb_change_indices ) call exchange_indices(node_list, my_id, n_cpu, .true.)
    fileout = 'jorek00000'
    call export_restart(node_list, element_list, fileout)
    if ( freeboundary .and. freeb_change_indices ) call exchange_indices(node_list, my_id, n_cpu, .false.)
  end if
  
  if ( ( my_id == 0 ) .and. ( (node_list%n_nodes > n_nodes_max+1000)                               &
    .or. (element_list%n_elements > n_elements_max+1000) ) ) then
    write(*,*) 'WARNING: n_nodes_max and/or n_elements_max is too large. This wastes memory.'
    write(*,*) '  n_nodes_max,    n_nodes    =', n_nodes_max,    node_list%n_nodes
    write(*,*) '  n_elements_max, n_elements =', n_elements_max, element_list%n_elements
    write(*,*) '  Note: for the equilibrium calculation higher values might be needed depending'
    write(*,*) '  on the resolution of your initial grid. In that case, you can run with reduced'
    write(*,*) '  values after restarting.'
  end if
  
  !***********************************************************************
  !***********************************************************************
  !*                          time stepping                              *
  !***********************************************************************
  !***********************************************************************
  
  if (nstep > 0) call update_deltas(my_id, node_list) ! create list of delta values in local_matrix module

  iter_gmres     = iter_precon
  iter_big       = gmres_max_iter
  iter_prev      = 0
  n_since_update = 0

  call tr_print_memsize("BeforeTimeStepping")
  call r3_info_print (-2, -2, 'INITIALIZATION')    ! timing
  
  if (.not. associated(aux_node_list)) allocate(aux_node_list) ! information of particle moments is stored in aux_list

  index_now = index_start  ! index_now: Index of current timestep

  jstep_loop: do jstep = 1, 10 ! Go through the different values of the tstep_n and nstep_n arrays
  istep_loop: do istep = 1, nstep_n(jstep)
    call clck_time_barrier(t_itstart)
    t0 = t_itstart

    flush stdout
    call tr_debug_write("JMAIN:Index_now",index_now)

    index_now = index_now + 1
    
    tstep = tstep_n(jstep)
    ! start from t=0 
    if ( index_now <= 1 ) tstep_prev = tstep
    
    if ( freeboundary ) call update_response(my_id,tstep, freeboundary_equil, resistive_wall)

    ! ---- For now running the jorek2_main should not include aux inputs
    aux_node_list%n_nodes = 0
    do i = 1, size(aux_node_list%node, 1)
      aux_node_list%node(i)%values = 0.d0
      aux_node_list%node(i)%deltas = 0.d0
    enddo

    if ( my_id == 0 ) then
      write(*,*) '******************************************************'
      write(*,'(A17,3i7,2f14.5,A)') ' *   time step : ',jstep,istep,index_now,tstep,tstep_prev,'  *'
      write(*,*) '******************************************************'
    end if

    ! --- Initialise the buffers needed by OpenMP threads. The values of n_tor, 
    ! --- n_plane, n_var have to remain the same until the end of the program.
    call new_thread_buffers()
    
    call update_equil_state(my_id,node_list, element_list, bnd_elm_list, xpoint, xcase)
    if ( my_id == 0 ) call print_equil_state(.false.)

    ! --- Prepare minor radius and q-,ft-,B-splines for bootstrap current
    minRad = 0.0
    if (bootstrap) then
      call bootstrap_find_minRad(node_list, element_list, ES%R_axis, ES%Z_axis, ES%psi_axis, ES%psi_bnd)
      call bootstrap_get_q_and_ft_splines(node_list, element_list, ES%psi_axis, ES%psi_xpoint, ES%R_xpoint, ES%Z_xpoint)
    endif
    
    call tr_debug_write("JMAIN:Find_axis_R",ES%R_axis)
    call tr_debug_write("JMAIN:Find_axis_Z",ES%Z_axis)
    call clck_time_barrier(t1)
    call clck_ldiff(t0,t1,tsecond)
!    if (my_id .eq. 0) then
!       write(*,FMT_TIMING)  my_id, '# Elapsed time init_time_step :',tsecond
!    end if

    ! Build the matrix 
    call clck_time_barrier(t0)
    if (gmres) then
      ! Matrix analysis and factorization in the preconditioner is re-done...
      ! ... in the first step of a simulation (also when restarting)
      ! ... when tstep changes
      ! ... when the previous time steps took too many iterations
      solve_only = (istep > 1) .and. ((iter_gmres+iter_prev <= 2*iter_precon) .and. (n_since_update < max_steps_noUpdate))
      if (solve_only) then 
        n_since_update = n_since_update + 1
      else
        n_since_update = 0
      endif
      !if ( my_id == 0 ) write(*,*) 'solve_only: ', solve_only
    endif
    
    if (use_pellet) then	    ! calculating the pellet_volume (total_pellet_volume)
      pellet_volume = PI * pellet_radius**2 * 2.d0 * PI * pellet_R * (pellet_phi/PI)
      call int3d_new(my_id, node_list, element_list, bnd_node_list, bnd_elm_list, exprs_all_int, res, 1)
    endif
    call tr_debug_write("JMAIN:Debconstruct_n_elms",n_local_elms)

    ! --- The following is for parallel debugging only

    !holder = 0;
    !write(*,*) "my_id", my_id, "PID", getpid(), "Host", name

    !do while (holder == 0)
    !  call sleep(5)
    !end do

    ! --- End of parallel debugging section 

    !--------- Constructing Global Matrix
    call construct_matrix(my_id, MPI_COMM_N, my_id_n, MPI_COMM_MASTER, my_id_master, local_elms,   &
         n_local_ELms, index_min(my_id+1), index_max(my_id+1), xpoint, xcase, ES%R_axis, ES%Z_axis,&
         ES%psi_axis, ES%psi_bnd, ES%R_xpoint, ES%Z_xpoint, ES%psi_xpoint, 1, n_tor,   &
         n_glob, nz_glob, ndof_glob, n_matrix_block_size, A_glob, rhs_glob, irn_glob, jcn_glob, ijA_index, ijA_size,    &
         irn_jcn, harmonic_matrix=.false.)


    call clck_time_barrier(t1)
    if (my_id .eq. 0) then
      call clck_ldiff(t0,t1,tsecond)
      write(*,FMT_TIMING) my_id, '# Elapsed time in construct global matrix :',tsecond
    endif     

    if (.not. gmres) then

      if (use_mumps) then
#ifdef USE_MUMPS
        call solve_mumps_all(my_id)
#endif
      elseif (use_strumpack) then
#ifdef USE_STRUMPACK
        call solve_strumpack_all(n_cpu,my_id,index_min(my_id+1),index_max(my_id+1))
#endif
      elseif (use_pastix) then
#if defined(USE_PASTIX) || defined(USE_PASTIX6)     
         call solve_pastix_all(n_cpu,my_id,index_min(my_id+1),index_max(my_id+1))
#endif
      endif

    else

      if (.not. solve_only) then

#ifndef DIRECT_CONSTRUCTION
        call clck_time(t0)
         ! --- Extract harmonic matrix from global matrix via MPI communication
        call distribute_harmonics(my_id,my_id_n,n_cpu)
        if(my_id_n.eq.0) call distribute_vector(rhs_glob,mumps_par%rhs,MPI_COMM_MASTER)
        call MPI_Barrier(MPI_COMM_WORLD,ierr)
        call clck_time_barrier(t1)
        call clck_ldiff(t0,t1,tsecond)
        if (my_id .eq. 0) then
          write(*,FMT_TIMING) my_id, '# Elapsed time distribute :',tsecond
        end if
#else 

         call clck_time_barrier(t0) 
         ! --- Direct construction of harmonic matrix
         call direct_construction_harmonic(my_id, my_id_n, m_cpu, n_cpu, MPI_COMM_N, MPI_COMM_MASTER, my_id_master, & 
              node_list, element_list, bnd_elm_list, bnd_node_list, xpoint, xcase, restart, freeboundary, .true.)
        call MPI_Barrier(MPI_COMM_WORLD,ierr)
        call clck_time_barrier(t1) 

        if (my_id .eq. 0) then
          call clck_ldiff(t0,t1,tsecond)
          write(*,FMT_TIMING) my_id, '# Elapsed time in construct harmonic matrix :',tsecond
        endif     

        call clck_time_barrier(t0) 
        ! --- Centralize the harmonic matrix on the master task of the MPI group (if needed)
        call centralization_harmonic(my_id, my_id_n, n_cpu_n, MPI_COMM_N)
        call MPI_Barrier(MPI_COMM_WORLD,ierr)
  
        call clck_time_barrier(t1) 

        if (my_id .eq. 0) then
          call clck_ldiff(t0,t1,tsecond)
          write(*,FMT_TIMING) my_id, '# Elapsed time in centralizing the matrix:',tsecond
        endif     

#endif

      else

        if(my_id_n.eq.0) call distribute_vector(rhs_glob,mumps_par%rhs,MPI_COMM_MASTER)
      endif

       ! --- Free the buffers needed by OpenMP threads (ELM-RHS etc.)
       call del_thread_buffers()

       call clck_time(t0)
      if (use_strumpack) then 
#ifdef USE_STRUMPACK
        call solve_matrix_n_spk(my_id,MPI_COMM_N,MPI_COMM_MASTER,solve_only)
#endif
      else
#if defined(USE_PASTIX) || defined(USE_MUMPS)
        call solve_matrix_n(my_id,MPI_COMM_N,MPI_COMM_MASTER,solve_only) ! factorise preconditioning matrices
#endif
#if defined(USE_PASTIX6)
        call solve_matrix_n_ptx(my_id,MPI_COMM_N,MPI_COMM_MASTER,solve_only)
#endif
      endif


      call clck_time_barrier(t1)
      call clck_ldiff(t0,t1,tsecond)
      if (my_id .eq. 0) then
        write(*,FMT_TIMING) my_id, '# Elapsed time first solve :',tsecond
      end if
    endif

    call clck_time(t0)
    if (gmres) then
      iter_prev = iter_gmres
      iter_gmres = gmres_max_iter

#ifdef USE_BICGSTAB
      call bicgstab_driver(irn_glob, jcn_glob, a_glob, deltas, rhs_glob, iter_gmres, gmres_tol, MPI_COMM_WORLD, MPI_COMM_N, MPI_COMM_MASTER)
#else
      call gmres_driver(my_id,my_id_n,MPI_COMM_N,MPI_COMM_MASTER,iter_gmres)
#endif

    endif
    call clck_time_barrier(t1)
    call clck_ldiff(t0,t1,tsecond)
    if (my_id .eq. 0) then
      write(*,FMT_TIMING)  my_id, '# Elapsed time gmres/solve :',tsecond
    end if

    call clck_time(t0)
    if ( (gmres .and. (iter_gmres .lt. iter_big)) .or. (.not.gmres) ) then

      if (use_pellet) then
        pellet_volume = total_pellet_volume
        call update_pellet(my_id,node_list,element_list)

        if (my_id == 0) then
          xtime_pellet_R(index_now)         = pellet_R
          xtime_pellet_Z(index_now)         = pellet_Z
          xtime_pellet_psi(index_now)       = pellet_psi
          xtime_pellet_particles(index_now) = pellet_particles
          xtime_phys_ablation(index_now)    = phys_ablation
        endif

      endif

#if (defined WITH_Neutrals) || (defined WITH_Impurities)
       if (using_spi) then
         n_spi_begin = 1
         do i = 1, n_inj !< Do one update for each injection location
           if (t_now >= t_ns(i)) call update_spi(my_id,node_list,element_list,i,n_spi_begin)
           n_spi_begin = n_spi_begin + n_spi(i)
         end do
       end if
#endif


      call update_values(my_id,element_list,node_list,deltas)         ! add solution to node values
      call update_deltas(my_id,node_list)

      t_now = t_now + tstep

      ! save previous time step
      tstep_prev = tstep


    else
      if ( my_id == 0 ) then
        write(*,*)
        write(*,'(a,i6.6,a)') '>>>>> NO CONVERGENCE AFTER ', iter_gmres, ' ITERATIONS. ABORTING <<<<<'
        write(*,*)
      end if
      index_now = index_now - 1 ! Undo the time step
      exit jstep_loop
    end if
    call clck_time_barrier(t1)
    call clck_ldiff(t0,t1,tsecond)
    if (my_id .eq. 0) then
      write(*,FMT_TIMING)  my_id, '#  Elapsed time Final Update:',tsecond
    end if

    !-------------------------------------------------------- adapt time step (in progress...)
    mindelta = minval(deltas); maxdelta = maxval(deltas);

    if (gmres .and. adaptive_time) then        ! experimental
       if (iter_gmres .ge. iter_big) then
          tstep = tstep /2.d0
          write(*,*) my_id,' REDUCTION TIMESTEP : ',tstep
       elseif (max(abs(mindelta),abs(maxdelta)) .gt. 0.05) then
          !	 tstep = tstep /2.d0
          !	 iter_gmres = 99999
          !	 write(*,*) my_id,' REDUCTION TIMESTEP : ',tstep
       elseif (max(abs(mindelta),abs(maxdelta)) .lt. 0.001) then
          !	 tstep = tstep * 2.d0
          !	 iter_gmres = 99999
          !	 write(*,*) my_id,' INCREASE TIMESTEP : ',tstep
       endif
    endif

    !--------------------------------------------------------- energies
    if ( (my_id == 0) .and. (.not. bench_without_plot) ) then

       call energy(W_mag,W_kin)

       R_axis_t(index_now)       = ES%R_axis
       Z_axis_t(index_now)       = ES%Z_axis
       psi_axis_t(index_now)     = ES%psi_axis
       R_xpoint_t(index_now,:)   = ES%R_xpoint(:)
       Z_xpoint_t(index_now,:)   = ES%Z_xpoint(:)
       psi_xpoint_t(index_now,:) = ES%psi_xpoint(:)
       R_bnd_t(index_now)        = ES%R_bnd
       Z_bnd_t(index_now)        = ES%Z_bnd
       psi_bnd_t(index_now)      = ES%psi_bnd

       xtime(index_now) = t_now
       energies(1:n_tor,1,index_now) = W_mag(1:n_tor)
       energies(1:n_tor,2,index_now) = W_kin(1:n_tor)

#ifdef JECCD
       call temp(node_list,element_list,A_tem,A_den,A_jen,A_jec,A_jec1,A_jec2)
       write(*,'(A,12e16.8)') ' current energies2 : ',A_tem,A_den
       write(*,'(A,12e16.8)') ' current energies3 : ',A_jen,A_jec
#ifdef JEC2DIAG
       write(*,'(A,12e16.8)') ' current energies4 : ',A_jec1,A_jec2
#endif

       energies2(1:n_tor,1,index_now) = A_tem(1:n_tor)
       energies2(1:n_tor,2,index_now) = A_den(1:n_tor)

       energies3(1:n_tor,1,index_now) = A_jen(1:n_tor)
       energies3(1:n_tor,2,index_now) = A_jec(1:n_tor)

#ifdef JEC2DIAG
       energies4(1:n_tor,1,index_now) = A_jec1(1:n_tor)
       energies4(1:n_tor,2,index_now) = A_jec2(1:n_tor)
#endif

       write(*,*) ' exiting current energies '
#endif

       ! --- Output some information about the current timestep
       130 format(1x,a,i5.5,a,es10.3,a)
       131 format(1x,a,2(2(es10.2,' ...',es10.2,',')))
       132 format(1x,'-------------------------------------------------------------------')
       133 format(1x,a,2(es10.2,' at ',i10,','))
       write(*,*)
       write(*,132)
       write(*,130) 'After step ', istep, ' (t_now=', t_now, '):'
       write(*,132)
       write(*,133) 'min,max deltas  =', mindelta, minloc(deltas), maxdelta, maxloc(deltas)
       write(*,131) 'W_mag,_kin      =', W_mag(1), W_mag(n_tor), W_kin(1), W_kin(n_tor)
       Growth_mag  = 0.d0; Growth_kin  = 0.d0; Growth_mag0 = 0.d0; Growth_kin0 = 0.d0
       if (index_now > index_start+1) then
         Growth_mag  = 0.5d0*log(abs(energies(n_tor,1,index_now)/energies(n_tor,1,index_now-1)))/ tstep
         if (energies(n_tor,2,index_now-1) .gt. 0.d0) then
           Growth_kin = 0.5d0*log(abs(energies(n_tor,2,index_now)/energies(n_tor,2,index_now-1)))/ tstep
         else
           Growth_kin = 0.d0
         endif
         Growth_mag0 = 0.5d0*log(abs(energies(1,1,index_now)/energies(1,1,index_now-1)))/ tstep
         if (energies(1,2,index_now-1) .gt. 0.d0) then
           Growth_kin0 = 0.5d0*log(abs(energies(1,2,index_now)/energies(1,2,index_now-1)))/ tstep
         else
           Growth_kin0 = 0.d0
         endif
         write(*,131) 'Growth_mag,_kin =', Growth_mag0, Growth_mag, Growth_kin0, Growth_kin
       endif
       write(*,132)
       write(*,*)
    endif   !--- my_id=0

    call int3d_new(my_id, node_list, element_list, bnd_node_list, bnd_elm_list, exprs_all_int, res, 1)

    if (my_id .eq. 0 ) then
      ! --- Output energies and growth_rates to text files during the code run
      call write_live_data(index_now)
      call write_live_data_vacuum(index_now)

#ifdef JECCD
      call write_live_data2(index_now)
      call write_live_data3(index_now)
#ifdef JEC2DIAG
      call write_live_data4(index_now)
#endif
#endif
    endif

    call clck_time_barrier(t1)
    call clck_ldiff(t0,t1,tsecond)
    if (my_id .eq. 0) then
       write(*,FMT_TIMING)  my_id, '#  Elapsed time Diagnostics :',tsecond
    end if
    !---------------------------------------------------------timing
    if ( istep == 1 ) then
       call r3_info_print (-3, -2, 'ITERATION	 1')
    else
       call r3_info_print (istep, -2, 'ITERATION')
    endif
    write(itlabel,'(I8)') istep
    call tr_print_memsize("AfterIter"//itlabel)
    
    ! --- Write a restart file every nout timesteps
    if ( (my_id == 0) .and. (mod(index_now,nout) == 0) ) then
      if ( freeboundary .and. freeb_change_indices ) call exchange_indices(node_list, my_id, n_cpu, .true.)
      write(fileout,'(A5,i5.5)') 'jorek',index_now
      call export_restart(node_list, element_list, fileout)
      if ( freeboundary .and. freeb_change_indices ) call exchange_indices(node_list, my_id, n_cpu, .false.)
    endif
    
    ! --- Exit the code if a file "STOP_NOW" exists in the run directory.
    inquire(file='STOP_NOW', exist=file_exists)
    if ( file_exists ) then
      if ( my_id == 0 ) then
        write(*,*)
        write(*,*) '>>>>> FOUND FILE STOP_NOW: EXITING THE CODE <<<<<'
        write(*,*)
      end if
      exit jstep_loop
    end if

    ! --- Redo LU decomposition if a file "REDO_LU" exists in the run directory.
    inquire(file='REDO_LU', exist=file_exists)
    if ( file_exists ) then
      if ( my_id == 0 ) then
        write(*,*)
        write(*,*) '>>>>> FOUND FILE REDO_LU: NEXT STEP DO AN LU DECOMPOSITION <<<<<'
        write(*,*)
        open(42, file='REDO_LU', iostat=ierr)
        if ( ierr == 0 ) close(42, status='delete')
      end if
      iter_prev  = iter_precon + 1
      iter_gmres = iter_precon + 1
    end if


    ! --- Exit the code if SIGTERM has been called on any node
    call MPI_ALLReduce(sigterm_called(), to_quit, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
    if (to_quit) then ! only present on id 0
      if ( my_id == 0 ) then
        write(*,*)
        write(*,*) ">>>>> SIGTERM RECEIVED: EXITING THE CODE <<<<<"
        write(*,*)
      end if
      exit jstep_loop
    end if
    
    ! --- Exit the code if NaNs are detected.
    if ( allocated(deltas) ) then
      sum_deltas = sum(deltas)
      if ( sum_deltas /= sum_deltas ) then
        write(*,*)
        write(*,*) '>>>>> NaNs DETECTED: EXITING THE CODE <<<<<'
        write(*,*)
        exit jstep_loop
      end if
    end if
    
    call clck_time_barrier(t1)
    call clck_ldiff(t_itstart,t1,tsecond)
    if (my_id .eq. 0) then
       write(*,FMT_TIMING)  my_id, '# Elapsed time ITERATION :',tsecond
    end if

  enddo istep_loop
  enddo jstep_loop
  
  !***********************************************************************
  !*                         cleanup  (solvers)                          *
  !***********************************************************************

  if (nstep .gt.0) then
#ifdef USE_MUMPS
    if (use_mumps) then
      mumps_par%JOB = -2                            ! clean up this instance of mumps
      call DMUMPS(mumps_par)
    endif
#endif

#ifdef USE_STRUMPACK
    if (use_strumpack) then
      call strumpack_finalize(MPI_COMM_WORLD)
    endif
#endif

#ifdef USE_BICGSTAB
    if (gmres) call bicgstab_finalize()
#endif

#ifdef USE_PASTIX6
    if (use_pastix) then
      call pastix_finalize()
    endif
#endif

#if defined(USE_PASTIX)
    if (use_pastix) then
      ! -- For PaStiX solver before version 6.x
      pastix_iparm(2)     = 7                       ! Clean-up
      pastix_iparm(3)     = 7

      if (.not. gmres) then
        call pastix_fortran(pastix_data,MPI_COMM_WORLD,mumps_par%n,DUMMY_INT,DUMMY_INT,DUMMY_REAL, &
          pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,1,pastix_iparm,pastix_dparm)
      elseif ( (.not. pastix_smp_only) .or. (pastix_smp_only .and. (my_id_n .eq.0))  ) then
        call pastix_fortran(pastix_data,MPI_COMM_N,mumps_par%n, &
          DUMMY_INT,DUMMY_INT,DUMMY_REAL,                       &
          pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,1,pastix_iparm,pastix_dparm)
      endif
    endif
#endif

#ifdef USE_WSMP
    if (use_wsmp) then
      call PWGSMP__deallocate()
    endif
#endif
    
  endif
  
  ! --- Close open files
  if ( my_id == 0 ) call finalize_live_data()
#ifdef JECCD
  if ( my_id == 0 ) call finalize_live_data2()
  if ( my_id == 0 ) call finalize_live_data3()
#ifdef JEC2DIAG
  if ( my_id == 0 ) call finalize_live_data4()
#endif
#endif

  !***********************************************************************
  !*                          plots etc.                                 *
  !***********************************************************************

  if (my_id .eq. 0)  then
    if ( freeboundary .and. freeb_change_indices ) call exchange_indices(node_list, my_id, n_cpu, .true.)
    fileout = 'jorek_restart'
    call export_restart(node_list, element_list, fileout)
    if ( freeboundary .and. freeb_change_indices ) call exchange_indices(node_list, my_id, n_cpu, .false.)
    if ( write_ps ) then
      if (.not. bench_without_plot) then
        do ivar=1,n_var
          call plot_solution(node_list,element_list,ivar,-1,1,variable_names(ivar))
        enddo

        do i=1,n_tor,2
          write(label,'(A4,i3,A1)') '(n =',((i-1)/2)*n_period,')'

          do ivar=1,n_var
            if ((ivar .ne. 3) .and. (ivar .ne. 4)) then
              call plot_solution(node_list,element_list,ivar,i,1,variable_names(ivar)//label)
            endif
          enddo

        enddo
      endif

      if (index_now .gt. 1) then

        E_min =  1.d20
        E_max = -1.d20
        E_max = max(E_max,maxval(energies(1,2,1:index_now)))
        E_min = min(E_min,minval(energies(1,2,1:index_now)))
        do i=2,n_tor
          E_max = max(E_max,maxval(energies(i,1,1:index_now)))
          E_min = min(E_min,minval(energies(i,1,1:index_now)))
          E_max = max(E_max,maxval(energies(i,2,1:index_now)))
          E_min = min(E_min,minval(energies(i,2,1:index_now)))
        enddo

        call nframe(1,1,2,xtime(1),xtime(index_now),E_min,E_max,'energies',7,'time',4,' ',1)

        do i=1,n_tor
          if (mod(i,2) .eq. 0) then
            call lincol(mod(i/2,10))
          else
            call lincol(mod((i-1)/2,10))
          endif
          call lplot(1,1,2,xtime(1:index_now),energies(i,1,1:index_now),-index_now,1,'Magnetic Energie',16,'time',4,'Emag',4)
          call lincol(4)
          if (n_tor .eq. 3) call lincol(2)
          call lplot(1,1,2,xtime(1:index_now),energies(i,2,1:index_now),-index_now,1,'Kinetic Energie',15,'time',4,'Ekin',4)
        enddo
        call lincol(3)
        call lplot(1,1,2,xtime(1:index_now),energies(1,2,1:index_now),-index_now,1,'Kinetic Energie',15,'time',4,'Ekin',4)
        call lincol(0)
      endif

!---------------------------------------------- plot equilibrium current profile (to be removed)

      nplot = 501
      call tr_allocate(xp,1,nplot,"xp",CAT_GRID)
      call tr_allocate(yp1,1,nplot,"yp1",CAT_GRID)
      call tr_allocate(yp2,1,nplot,"yp2",CAT_GRID)
      call tr_allocate(yp3,1,nplot,"yp3",CAT_GRID)
      ! ---- plot neoclassical coefficients -----
      if (NEO) then
        call tr_allocate(mu_neo,1,nplot,"mu_neo",CAT_GRID)
        call tr_allocate(ki_neo,1,nplot,"ki_neo",CAT_GRID)
      endif
      iplot = 0

      Rp_start = ES%R_axis - amin*2.d0
      Rp_end   = ES%R_axis + amin*2.d0

      Zp = ES%Z_axis

      do i=1,nplot

        Rp =  Rp_start + float(i-1)/float(nplot-1) * (Rp_end - Rp_start)

        call find_RZ(node_list,element_list,Rp,Zp,R_out,Z_out,i_elm,s_out,t_out,ifail)

        if (ifail .eq. 0) then

    	  call interp(node_list,element_list,i_elm,1,1,s_out,t_out,psi,P_s,P_t,P_st,P_ss,P_tt)

    	  call density(    xpoint,xcase, Zp, ES%Z_xpoint, psi,ES%psi_axis,ES%psi_bnd,	       &
    	     zn,dn_dpsi,dn_dz,dn_dpsi2,dn_dz2,dn_dpsi_dz,dn_dpsi3,dn_dpsi_dz2,dn_dpsi2_dz)
    	  if (with_TiTe) then	     
    	    call temperature_i(xpoint,xcase, Zp, ES%Z_xpoint, psi,ES%psi_axis,ES%psi_bnd, &
    	      zTi,dTi_dpsi,dTi_dz,dTi_dpsi2,dTi_dz2,dTi_dpsi_dz,dTi_dpsi3,dTi_dpsi_dz2,dTi_dpsi2_dz)			   
    	    call temperature_e(xpoint,xcase, Zp, ES%Z_xpoint, psi,ES%psi_axis,ES%psi_bnd, &
    	      zTe,dTe_dpsi,dTe_dz,dTe_dpsi2,dTe_dz2,dTe_dpsi_dz,dTe_dpsi3,dTe_dpsi_dz2,dTe_dpsi2_dz)	     
            zT = zTi + zTe
    	    dT_dpsi = dTi_dpsi + dTe_dpsi	    
    	  else
            call temperature(xpoint,xcase, Zp, ES%Z_xpoint, psi,ES%psi_axis,ES%psi_bnd, &
    	      zT,dT_dpsi,dT_dz,dT_dpsi2,dT_dz2,dT_dpsi_dz,dT_dpsi3,dT_dpsi_dz2,dT_dpsi2_dz)
    	  endif
    	  call FFprime(    xpoint,xcase, Zp, ES%Z_xpoint, psi,ES%psi_axis,ES%psi_bnd,	       &
    	       zFFprime,dFFprime_dpsi,dFFprime_dz,dFFprime_dpsi2,dFFprime_dz2,dFFprime_dpsi_dz, .true.)

          if (NEO) then
            if (num_neo_file) then
              call neo_coef (xpoint, xcase, Zp, ES%Z_xpoint, psi, ES%psi_axis,ES%psi_bnd, &
                amu_neo_node, aki_neo_node)
            endif
          endif

          zjz	= (zFFprime - Rp*Rp * (zn * dT_dpsi + dn_dpsi * zT)) / Rp

          iplot = iplot + 1

          xp(iplot)  = Rp
          yp1(iplot) = zFFprime / Rp
          yp2(iplot) = zjz
          yp3(iplot) = - Rp*Rp * (zn * dT_dpsi + dn_dpsi * zT) / Rp

          if (NEO) then
            if ( num_neo_file) then
              mu_neo(iplot) = amu_neo_node
              ki_neo(iplot) = aki_neo_node
              write(*,'(A,8e16.8)') ' profiles : ',xp(iplot),psi,ES%psi_axis,ES%psi_xpoint,mu_neo(iplot),ki_neo(iplot)
            endif
          endif

        endif

      enddo

      call lplot6(1,1,xp,yp2,iplot,' ')
      call lincol(1)
      call lplot6(1,1,xp,yp1,-iplot,' ')
      call lincol(2)
      call lplot6(1,1,xp,yp3,-iplot,' ')
      call lincol(0)
      if (NEO) then
        if ( num_neo_file) then
          call lplot6(1,1,xp,mu_neo,iplot,' ')
          call lincol(1)
          call lplot6(1,1,xp,ki_neo,iplot,' ')
          call lincol(0)
        end if
      endif
      call finplt 					 ! close plot file
    endif !  write_ps
!  cll export_POV(node_list,element_list,3,1)	       ! export to POVray native bezier patch format
#ifdef fullmhd
    write(*,*) ' '
    write(*,*) 'Warning: Export to helena is not adapted for full MHD'
    write(*,*) ' '
#else
    call export_helena(node_list,element_list,bnd_elm_list)
#endif
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

  call r3_info_summary ()                                ! timing
  call MPI_FINALIZE(IERR)                                ! clean up MPI

end program JOREK2


! --- The following comment block defines the start page of the Doxygen code documentation ---
!
!> \mainpage
!!
!! \section JOREK About JOREK
!!
!! JOREK solves the (reduced) MHD equations in 3D toroidal geometry using a discretization with
!! Bezier finite elements in the poloidal plane and a Fourier expansion in toroidal direction.
!!
!! \section DOCU About this documentation
!!
!! This documentation covers the code structure, i.e.,
!! - <a href="dirs.html">directories</a>,
!! - <a href="files.html">source files</a>,
!! - <a href="annotated.html">Fortran modules</a>,
!! - <a href="globals_func.html">subroutines and functions</a>, and
!! - routine parameters.
!!
!! <img src="JOREK_DOC.jpeg" style="border:1px solid black;float:right"/>
!! Background information on the following topics is covered in a seperate documentation
!! available in the docu/ folder of the JOREK repository as LaTeX source or
!! <a href="JOREK_DOC.pdf">online as PDF</a>:
!! - Implemented equations
!! - Numerical methods
!! - Using the version control system Subversion
!! - Compiling and running the code
!! - Bezier elements
!! - Free boundary extension
!! - ...
!! 
!! \section Repository Browse Subversion Repository online
!! 
!! The <a href="http://subversion.apache.org/">Subversion</a> repository of JOREK can
!! be <a href="https://gforge.inria.fr/scm/viewvc.php?view=rev&root=aster">browsed online</a>
!! after logging in if you already have an account for it. On the Server,
!! <a href="http://www.viewvc.org/">ViewVC</a> is installed for this purpose.
!!
!! \section Doxygen Code documentation with Doxygen
!!
!! This documentation is generated directly from the source code using
!! <a href="http://www.doxygen.org">Doxygen</a>.
!! For this to work properly, you should follow some
!! simple rules when writing comments in the code.
!!
!! - An example for the proper documentation of a subroutine:                               \n
!!                                                                                          \n
!! <code>
!! !\> Brief documentation for the test routine                                             \n
!! !!                                                                                       \n
!! !! More details on the functionality of the                                              \n
!! !! test routines can be added like this.                                                 \n
!! !!                                                                                       \n
!! subroutine test(i,r)                                                                     \n
!!                                                                                          \n
!!   ! --- Routine parameters                                                               \n
!!   integer, intent(in)  :: i !< Information on parameter i                                \n
!!   real*8,  intent(out) :: r !< Information regarding parameter r                         \n
!!                                                                                          \n
!!   ! --- Local variables                                                                  \n
!!   integer :: k   ! Information on k (Local variables are not documented by Doxygen)      \n
!!                                                                                          \n
!!   ...                                                                                    \n
!!                                                                                          \n
!! end subroutine test                                                                      \n
!! </code>
!!
!! - The documentation of a module works in the same way:                                   \n
!!                                                                                          \n
!! <code>
!! !\> Brief documentation for the test module                                              \n
!! !!                                                                                       \n
!! !! More details on the functionality of the                                              \n
!! !! test module can be added like this.                                                   \n
!! !!                                                                                       \n
!! module some_test_module                                                                  \n
!!                                                                                          \n
!!  implicit none                                                                           \n
!!                                                                                          \n
!!  !\> \@name Rectangular Grid                                                             \n
!!  integer :: n_R               !< Number of grid points in R-direction                    \n
!!  integer :: n_Z               !< Number of grid points in Z-direction                    \n
!!                                                                                          \n
!!  !\> \@name Polar Grid                                                                   \n
!!  !! Parameters defining a non flux-aligned polar grid in the poloidal plane.             \n
!!  integer :: n_radial          !< Number of radial grid points                            \n
!!  integer :: n_pol             !< Number of poloidal grid points                          \n
!!                                                                                          \n
!!  contains                                                                                \n
!!                                                                                          \n
!!  !\> Description for routine                                                             \n
!!  subroutine test()                                                                       \n
!!  ...                                                                                     \n
!! </code>
!! Note, that the \@name command allows to create groups of variables. An example for this
!! is the module ::phys_module.
!!
!! - HTML markup is possible, e.g., \<b\>some text\</b\> will display in bold face:
!!   <b>some text</b>
!!
!! - For further information, refer to the
!!   <a href="http://www.stack.nl/~dimitri/doxygen/manual.html">Doxygen manual</a>.
!!
