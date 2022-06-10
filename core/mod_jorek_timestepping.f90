module mod_jorek_timestepping
use mod_event
use mod_particle_sim
use iso_c_binding ! for fftw03.f03
use mod_parameters, only: n_plane
use data_structure, only: type_bnd_element_list, type_bnd_node_list !< store these in jorek_timestep_action

! Solvers
use mumps_module
use pastix_module
use wsmp_module

#ifdef USE_STRUMPACK
use strumpack_module
#endif
use preconditioner_module
use mod_distribute_preconditioner
use direct_construction_mod
use centralization_mod

use equil_info

implicit none

private
public jorek_timestep_action, new_jorek_timestep_action

#ifdef USE_FFTW
  include 'fftw3.f03'
#endif

type, extends(action) :: jorek_timestep_action
  integer :: istep !< index in timestep size array from namelist (not jorek timestep number!)

  logical :: setup_done = .false. !< have we set up the solvers etc?
#ifdef USE_FFTW
  real*8     :: in_fft(1:n_plane)
  complex*16 :: out_fft(1:n_plane)
#endif

  type(type_bnd_element_list) :: bnd_elm_list !< List of boundary elements
  type(type_bnd_node_list)    :: bnd_node_list !< List of boundary nodes.

  type (t_equil_state) :: eq !< Information about the equilibrium

  ! GMRES communicators (TODO: move to separate type?)
  ! if no gmres is used mpi_comm_n and my_id_n are equal to the global ones (mpi_comm_world and my_id)
  integer :: my_id_n, n_cpu_n !< id, count of procs in local comm (for one harmonic)
  integer :: my_id_trans, n_cpu_trans !< id, count of procs in transverse comm (all first of MPI_COMM_N, all second, all third)
  integer :: my_id_master
  integer :: MPI_COMM_N !< group for each harmonic (see [[gmres_setup]])
  integer :: MPI_COMM_TRANS !< transversal groups (i.e. every 1st of MPI_COMM_N, every 2nd of MPI_COMM_N etc) (see [[gmres_setup]])
  integer :: MPI_COMM_MASTER !< Every first of MPI_COMM_N (see [[gmres_setup]])
  integer :: MPI_GROUP_WORLD
  integer :: MPI_GROUP_MASTER !< subset of MPI_COMM_WORLD corresponding to MPI_COMM_MASTER
  integer, allocatable :: i_tor(:) !< toroidal harmonic solved by this process
  integer, allocatable :: local_elms(:), index_min(:), index_max(:) !< division of work across processes
  integer :: n_local_elms
  integer :: n_AA !< number of nonzeros

  ! when to recalculate the preconditioner
  logical :: prec_needed = .true.
  integer :: iter_gmres, iter_prev

  ! Coupling data for in construct_matrix
  type(type_node_list), pointer :: auxiliary_node_list => null()

  ! Optionally update the start time of anohter event
  type(event), pointer :: extra_event => null()
contains
  procedure :: do => do_jorek_timestep
end type
interface jorek_timestep_action
  module procedure new_jorek_timestep_action
end interface

contains


function new_jorek_timestep_action(auxiliary_node_list) result(new)
  type(jorek_timestep_action) :: new
  type(type_node_list), intent(in), target,  optional :: auxiliary_node_list
!  if (present(auxiliary_node_list)) new%auxiliary_node_list => auxiliary_node_list
  new%istep = 1
  new%name = "JOREK timestep"
  new%log = .true.
end function new_jorek_timestep_action


!> Initialise parameters for the solvers used in JOREK and perform
!> other preliminary work to be done once before solving
subroutine setup_solvers(this, sim)
  use phys_module
  use nodes_elements
  use mod_clock,            only: clck_init
  use live_data,            only: init_live_data
  use mod_live_data_core,   only: write_live_data_all
  use mpi_mod
  use tr_module
  use global_distributed_matrix
  use mod_boundary,         only: boundary_from_grid
  use mod_global_matrix_structure
  use gmres_setup,          only: gmres_setup_jorek
  use vacuum
  use vacuum_response,      only: get_vacuum_response, update_response, init_wall_currents, I_coils
  use vacuum_equilibrium,   only: import_external_fields
  use mod_startup_teardown, only: sanity_checks
  use mod_log_params,       only: log_parameters

  implicit none

  class(jorek_timestep_action), intent(inout) :: this
  type(particle_sim),           intent(inout) :: sim

  real*8 :: psi_xpoint(2), R_xpoint, Z_xpoint, s_xpoint, t_xpoint
  integer :: i_elm_xpoint, ifail
  real*8 :: psi_lim, R_lim, Z_lim

  integer :: index_size, id_elements !< number of elements locally
  integer :: inode, ierr, i, block_size, n_masters

  write(*,*) 'setting up the solvers'
  call tr_meminit(sim%my_id, sim%n_cpu)
  call tr_resetfile()

  ! Initialize clock
  call clck_init()
  call r3_info_init()

  ! Initialise mode and mode_type arrays
  call det_modes()

  ! Initialise the data writing 
  if (sim%my_id .eq. 0) then
    call init_live_data()

    if (restart) then
      do i = 1, index_start
       call write_live_data_all(i)
!      call write_live_data_vacuum(index_now, diag_coil_curr)
      end do
    endif
  endif

  ! --- Preset some solver variables
  pastix_initialised = .false.
  pastix_analysed    = .false.

  ! MURGE with ntor=1 doesn't work
  if (n_tor .eq. 1) then
    gmres     = .false.
  end if

  ! --- Initialize the vacuum part.
  call vacuum_init(sim%my_id, freeboundary_equil, freeboundary, resistive_wall)

  ! --- Set time-stepping scheme
  call update_time_evol_params()

  ! --- Write out all parameters defined in parameters and the namelist input file.
  call log_parameters(sim%my_id)

  ! Warn on doing stupid stuff
  call sanity_checks(sim%my_id, sim%n_cpu, 7, 7) ! #### the 7, 7 is just a dummy that needs to be removed later on; the sanity_checks should anyway not be part of setup_solvers in the end (to be addressed in a separate pull request) @TODO
  if (nstep .gt. 0)   call check_preconditioner_consistency

  ! Initialise the boundary element and node list
  if (sim%my_id .eq. 0) then
    call boundary_from_grid(sim%fields%node_list, sim%fields%element_list, bnd_node_list, bnd_elm_list, output_bnd_elements)
  endif

  call broadcast_boundary(sim%my_id, bnd_elm_list, bnd_node_list)

  ! --- Fill the vacuum response matrices for freeboundary computations
  if ( freeboundary ) then
  
    call get_vacuum_response(sim%my_id, sim%fields%node_list, bnd_elm_list, bnd_node_list, freeboundary_equil, resistive_wall)

    call update_response(sim%my_id,get_tstep_n(1), freeboundary_equil, resistive_wall)
    
    call import_external_fields('coil_field.dat', sim%my_id)
    
    call set_coil_curr_time_trace()
    
    call MPI_BCAST(wall_curr_initialized, 1 , MPI_LOGICAL,          0, MPI_COMM_WORLD, ierr)

    if ( (.not. restart) .or. (.not. wall_curr_initialized) ) then
      call init_wall_currents(sim%my_id, resistive_wall)
    endif
  
  end if


  if (RMP_on) then
     if (sim%my_id == 0) then
        call read_RMP_profiles(bnd_node_list)
     endif
     call broadcast_RMP_profiles(sim%my_id, bnd_node_list)        ! psi_RMP profiles
  endif

  ! nodes, elements, bnd_nodes and phys have already been broadcast
  if ( freeboundary ) call broadcast_vacuum(sim%my_id, resistive_wall)

  this%n_AA = 0
  do inode = 1, sim%fields%node_list%n_nodes  
    this%n_AA = max(this%n_AA,sim%fields%node_list%node(inode)%index(4))  
  end do
  mumps_par%n = this%n_AA

  call update_equil_state(sim%my_id, sim%fields%node_list, sim%fields%element_list, bnd_elm_list, xpoint, xcase)
  this%eq = ES

  if ( sim%my_id == 0 ) then
    call print_equil_state(.true.)
    call save_special_points('special_equilibrium_points.dat', .false., ierr)
  end if

   ! --- Initialize FFTW
#ifdef USE_FFTW
  call dfftw_plan_dft_r2c_1d(fftw_plan,n_plane,this%in_fft,this%out_fft,FFTW_PATIENT)
#endif

  n_masters = (n_tor+1)/2
  if (gmres) then
    ! setup per-harmonic and transverse communicators
    !call gmres_setup_jorek(sim%my_id, sim%n_cpu, this%i_tor, this%my_id_n, this%n_cpu_n, &
    !                       this%my_id_trans, this%n_cpu_trans, this%my_id_master,        &
    !                       this%MPI_COMM_N, this%MPI_COMM_TRANS, this%MPI_COMM_MASTER,   &
    !                       this%MPI_GROUP_MASTER, this%MPI_GROUP_WORLD, my_family_id)
    
    call create_communicators(this%my_id_n, this%n_cpu_n, this%MPI_COMM_N, this%my_id_master, n_masters, &
                             this%MPI_COMM_MASTER, this%MPI_COMM_TRANS)

    call distribute_modes                       
                      
  else
     this%my_id_n    = sim%my_id
     this%n_cpu_n    = sim%n_cpu
     this%MPI_COMM_N = MPI_COMM_WORLD
  endif

  !***********************************************************************
  !*              distribute nodes and elements over cpu's               *
  !***********************************************************************
  index_size  = sim%n_cpu
  id_elements = sim%my_id

  call tr_allocate(this%local_elms,1,sim%fields%element_list%n_elements,"local_elms",CAT_FEM)
  call tr_allocate(this%index_min,1,index_size,"index_min",CAT_FEM)
  call tr_allocate(this%index_max,1,index_size,"index_max",CAT_FEM)
  call tr_allocate(local_index_start,1,sim%n_cpu,"local_index_start",CAT_FEM)
  call tr_allocate(local_index_end,1,sim%n_cpu,"local_index_end",CAT_FEM)
 !
  ! Construct index_min, index_max and local_elems
  !
  call distribute_nodes_elements(id_elements,this%n_cpu_n,index_size,sim%fields%node_list,sim%fields%element_list, .false., &
                                 this%local_elms, this%n_local_elms, ndof_glob, this%index_min, this%index_max, restart, freeboundary)
                                          
  sim%fields%node_list%n_dof = ndof_glob
  local_index_start = this%index_min
  local_index_end   = this%index_max
 
  call update_deltas(sim%my_id, sim%fields%node_list) ! create list of delta values in local_matrix module
  
   ! Build ijA_index, ijA_size and irn_jcn
 
  block_size = n_tor*n_var

  call global_matrix_structure(sim%my_id,this%my_id_n, sim%fields%node_List, sim%fields%element_list, bnd_elm_list, freeboundary,&
                               this%local_elms,this%n_local_elms,this%index_min(id_elements+1),this%index_max(id_elements+1),      &
                               ijA_index, ijA_size, irn_jcn, irn_glob, jcn_glob, 1, n_tor, n_glob, nz_glob, ndof_glob, n_matrix_block_size)                      

  if ( freeboundary .and. ( sr%n_tor /= 0 ) ) then 
    call global_matrix_structure_vacuum(sim%fields%node_list, bnd_node_list, this%index_min(sim%my_id+1), this%index_max(sim%my_id+1), &
                                        1, n_tor, irn_glob, jcn_glob, n_matrix_block_size, ijA_index, ijA_size, irn_jcn) 
  endif
  
  if ((gmres) .and. (this%my_id_n .eq. 0)) call map_row_index(ndof_glob)

  call MPI_Barrier(MPI_COMM_WORLD,ierr)

  if (use_mumps) then
    if (.not. gmres) then
      call initialise_mumps(MPI_COMM_WORLD) ! start MUMPS sparse matrix solver all cpus
    else
      call initialise_mumps(this%MPI_COMM_N) ! start MUMPS sparse matrix solver on local groups
    endif
  endif

  this%iter_gmres  = iter_precon
  this%iter_prev   = 0

  this%setup_done = .true.

  if (.not. associated(aux_node_list)) allocate(aux_node_list) ! information of particle moments is stored in aux_list

end subroutine setup_solvers


!> Destroy memory used by these solvers
!> note that this is not yet explicitly called somewhere!
subroutine cleanup_solvers(this, sim)

  use phys_module, only: gmres, use_mumps, use_pastix
  use mpi_mod,     only: MPI_COMM_WORLD
  class(jorek_timestep_action), intent(inout) :: this
  type(particle_sim), intent(inout)           :: sim

  integer :: DUMMY_INT
  real*8  :: DUMMY_REAL

  if (use_mumps) then

#ifdef USE_MUMPS
    mumps_par%JOB = -2                            ! clean up this instance of mumps
    call DMUMPS(mumps_par)
#endif

  elseif (use_pastix) then

    pastix_iparm(2)     = 7                       ! Clean-up
    pastix_iparm(3)     = 7

    if (.not. gmres) then

      call pastix_fortran(pastix_data,MPI_COMM_WORLD,mumps_par%n,DUMMY_INT,DUMMY_INT,DUMMY_REAL, &
                          pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,1,pastix_iparm,pastix_dparm)

    elseif ( (.not. pastix_smp_only) .or. (pastix_smp_only .and. (this%my_id_n .eq.0))  ) then

      call pastix_fortran(pastix_data,this%MPI_COMM_N,mumps_par%n,&
                          DUMMY_INT,DUMMY_INT,DUMMY_REAL, &
                          pastix_perm_vars,pastix_iperm_vars,mumps_par%rhs,1,pastix_iparm,pastix_dparm)
    endif

  endif

end subroutine cleanup_solvers


!> Perform a single jorek timestep, with timestep size from current time - last time
!> Notes:
!> - STOP_NOW file handling does not work
!> - set_trap_sigterm() is not setup yet
subroutine do_jorek_timestep(this, sim, ev)
  use phys_module
  use nodes_elements
  use mod_clock
  use global_distributed_matrix
  use data_structure,          only: new_thread_buffers, del_thread_buffers
  use mod_bootstrap_functions, only: bootstrap_find_minRad, bootstrap_get_q_and_ft_splines
  use live_data
  use mod_live_data_core,      only: write_live_data_all
  use tr_module,               only: tr_print_memsize, tr_resetfile
  use mod_export_restart
  use construct_matrix_mod
  use solve_mat_n
  use pellet_module
  use vacuum
  use vacuum_response,         only: update_response
  use mod_fields_linear
  use mod_gmres,               only: gmres_driver
#ifdef USE_BICGSTAB
  use mod_bicgstab, only: bicgstab_driver, bicgstab_finalize
#endif
  use mod_expression,          only: exprs_all_int, init_expr
  use mod_integrals3D

#if (defined WITH_Neutrals) && (!defined WITH_Impurities)
  use mod_neutral_source
#endif
#ifdef WITH_Impurities
  use mod_injection_source
#endif

  class(jorek_timestep_action), intent(inout) :: this
  type(particle_sim), intent(inout)           :: sim
  type(event), intent(inout), optional        :: ev

  real*8         :: dt, dt_jorek
  type(clcktype) :: t0, t1, t_itstart
  real*8         :: tsecond
  logical        :: solve_only

  real*8         :: W_mag(n_tor), W_kin(n_tor), growth_mag, growth_kin, growth_mag0, growth_kin0
  real*8         :: density_tot,density_in,density_out,pressure_tot,pressure_in,pressure_out,Bgeo
  real*8         :: kin_par_tot, kin_par_in, kin_par_out, mom_par_tot, mom_par_in, mom_par_out
  real*8, allocatable :: res(:)

  real*8         :: mindelta, maxdelta, sum_deltas
  character*8    :: label, itlabel
  character*14   :: fileout
  integer        :: i, n_spi_begin

  call init_expr()
  allocate(res(exprs_all_int%n_expr+1))
  res = 0.d0  

  ! Get the timestep size
  dt_jorek = get_tstep_n(this%istep)
  if (dt_jorek .eq. 0.d0) then
    write(*,*) "Jorek timestep is 0, assuming end of simulation."
    sim%stop_now = .true.
    return
  end if
  tstep = dt_jorek !< Update the jorek timestep for use in mod_elt_matrix
  !< Update the jorek previous timestep for use in mod_elt_matrix. 
  !< If ommited, certain models (e.g. 710+) will divide by zero. Not fully tested.
  if ( this%istep -1 > 0) then
    tstep_prev = get_tstep_n(this%istep-1) 
  else
    tstep_prev = tstep
    write(*,*) "INFO: tstep_prev set to tstep at first iteration"
  endif
  dt = dt_jorek * sim%t_norm

  if (.not. this%setup_done) then
    t_now = sim%time / sim%t_norm  
    index_now = index_start
    if (sim%my_id .eq. 0) write(*,"(A,f16.8,A,g12.6,A)") "INFO: JOREK timestep: ", dt_jorek, " = ", dt, " s"
    call setup_solvers(this, sim)
  end if

  index_now = index_now + 1 ! we started at 0

  ! Set up the next start time to run this
  if (present(ev)) then
    ev%start = ev%start + dt
    if (sim%my_id .eq. 0) write(*,*) "INFO: scheduling next JOREK event for ", ev%start
  end if

  if (associated(this%extra_event)) then
    this%extra_event%start = ev%start
  endif

  call clck_time_barrier(t_itstart)
  t0 = t_itstart

  if ( freeboundary ) call update_response(sim%my_id,dt_jorek, freeboundary_equil, resistive_wall)
  
  ! --- Initialise the buffers needed by OpenMP threads. The values of n_tor, 
  ! --- n_plane, n_var have to remain the same until the end of the program.
  call new_thread_buffers()

  call update_equil_state(sim%my_id, sim%fields%node_list, sim%fields%element_list, bnd_elm_list, xpoint, xcase )
  this%eq = ES

  if ( sim%my_id == 0 ) call print_equil_state(.false.)
  
  ! --- Prepare minor radius and q-,ft-,B-splines for bootstrap current
  minRad=0.d0
  if (bootstrap) then
    call bootstrap_find_minRad(sim%fields%node_list, sim%fields%element_list, this%eq%R_axis, this%eq%Z_axis, this%eq%psi_axis, this%eq%psi_bnd)

    call bootstrap_get_q_and_ft_splines(sim%fields%node_list, sim%fields%element_list, this%eq%psi_axis, this%eq%psi_xpoint, this%eq%R_xpoint, this%eq%Z_xpoint)
  endif
  
  call clck_time_barrier(t1)
  call clck_ldiff(t0,t1,tsecond)

  ! Build the matrix 
  call clck_time_barrier(t0)
  if (gmres) then
    ! Matrix analysis and factorization in the preconditioner is re-done...
    ! ... in the first step of a simulation (also when restarting)
    ! ... when tstep changes
    ! ... when the previous time steps took too many iterations
    solve_only = (this%istep .gt. 1) .and. (this%iter_gmres+this%iter_prev <= 2*iter_precon)
!    solve_only = (.not. this%prec_needed) .and. (this%iter_gmres+this%iter_prev <= 2*iter_precon)
  endif

  if (use_pellet) then            ! calculating the pellet_volume (total_pellet_volume)
    pellet_volume = PI * pellet_radius**2 * 2.d0 * PI * pellet_R * (pellet_phi/PI)
    call Integrals_3D(sim%my_id, sim%fields%node_list, sim%fields%element_list, density_tot,density_in,density_out,pressure_tot,pressure_in,pressure_out, &
                                                                                kin_par_tot, kin_par_in, kin_par_out, mom_par_tot, mom_par_in, mom_par_out)
    
  endif


  call construct_matrix(sim%my_id, this%MPI_COMM_N, this%my_id_n, this%MPI_COMM_MASTER, this%my_id_master,               &
                        this%local_elms, this%n_local_elms, this%index_min(sim%my_id+1),                                 &
                        this%index_max(sim%my_id+1), xpoint, xcase, this%eq%R_axis, this%eq%Z_axis, this%eq%psi_axis,    &
                        this%eq%psi_bnd, this%eq%R_xpoint, this%eq%Z_xpoint, this%eq%psi_xpoint, 1, n_tor,               &
                        n_glob, nz_glob, ndof_glob, n_matrix_block_size, A_glob, rhs_glob, irn_glob, jcn_glob, ijA_index, ijA_size, &
                        irn_jcn,  harmonic_matrix=.false.)                    
  
  ! --- Free the buffers needed by OpenMP threads (ELM-RHS etc.)
  call del_thread_buffers()

  call clck_time_barrier(t1)
  if (sim%my_id .eq. 0) then
     call clck_ldiff(t0,t1,tsecond)
    write(*,FMT_TIMING) sim%my_id, '# Elapsed time construct_matrix :',tsecond
  endif     

  if (.not. gmres) then
    if (use_mumps) then
#ifdef USE_MUMPS    
      call solve_mumps_all(sim%my_id)
#endif      
    else
#if defined(USE_PASTIX) || defined(USE_PASTIX6)    
      call solve_pastix_all(sim%n_cpu,sim%my_id,this%index_min(sim%my_id+1),this%index_max(sim%my_id+1))
#endif      
    endif
  else
    call clck_time(t0)
    if (.not. solve_only) then
      call distribute_harmonics(sim%my_id,this%my_id_n,sim%n_cpu)
    endif
    if(this%my_id_n.eq.0) call distribute_vector(rhs_glob,mumps_par%rhs,this%MPI_COMM_MASTER)
    
    call clck_time_barrier(t1)
    call clck_ldiff(t0,t1,tsecond)
    if (sim%my_id .eq. 0) write(*,FMT_TIMING) sim%my_id, '# Elapsed time distribute :',tsecond

    call clck_time(t0)
    call solve_matrix_n(sim%my_id,this%MPI_COMM_N,this%MPI_COMM_MASTER,solve_only)    ! factorise preconditioning matrices
      
    call clck_time_barrier(t1)
    call clck_ldiff(t0,t1,tsecond)
    if (sim%my_id .eq. 0) write(*,FMT_TIMING) sim%my_id, '# Elapsed time first solve :',tsecond
  endif

  call clck_time(t0)
  if (gmres) then
    this%iter_prev = this%iter_gmres
    this%iter_gmres = gmres_max_iter
 
#ifdef USE_BICGSTAB
    call bicgstab_driver(irn_glob, jcn_glob, a_glob, deltas, rhs_glob, &
                     this%iter_gmres, gmres_tol, MPI_COMM_WORLD, this%MPI_COMM_N, this%MPI_COMM_MASTER)
#else 
    call gmres_driver(sim%my_id,this%my_id_n,this%MPI_COMM_N,this%MPI_COMM_MASTER,this%iter_gmres)
#endif

  endif
  call clck_time_barrier(t1)
  call clck_ldiff(t0,t1,tsecond)
  if (sim%my_id .eq. 0) write(*,FMT_TIMING) sim%my_id, '# Elapsed time gmres/solve :',tsecond

  call clck_time(t0)
  if ( (gmres .and. (this%iter_gmres .lt. gmres_max_iter)) .or. (.not. gmres) ) then

    ! TODO add if use_pellet

#if (defined WITH_Neutrals) || (defined WITH_Impurities)
    if (using_spi) then
      n_spi_begin = 1
      do i = 1, n_inj !< Do one update for each injection location
        if (t_now >= t_ns(i)) call update_spi(sim%my_id,sim%fields%node_list,sim%fields%element_list,i,n_spi_begin)
        n_spi_begin = n_spi_begin + n_spi(i)
      end do
    end if
#endif    

    call update_values(sim%my_id, sim%fields%element_list, sim%fields%node_list, deltas)         ! add solution to node values
    call update_deltas(sim%my_id, sim%fields%node_list)
    t_now = t_now + dt_jorek
  else
    if ( sim%my_id == 0 ) then
      write(*,*)
      write(*,'(a,i6.6,a)') '>>>>> NO CONVERGENCE AFTER ', this%iter_gmres, ' ITERATIONS. ABORTING <<<<<'
      write(*,*)
    end if
    sim%stop_now = .true.
    return
  end if
  call clck_time_barrier(t1)
  call clck_ldiff(t0,t1,tsecond)
  if (sim%my_id .eq. 0) write(*,FMT_TIMING) sim%my_id, '#  Elapsed time Final Update:',tsecond

  !--------------------------------------------------------- energies
  if (sim%my_id == 0) then
    ! This is a change from jorek2_main, where these quantities are calculated using the old xpoint and axis data
    call update_equil_state(sim%my_id,sim%fields%node_list, sim%fields%element_list, bnd_elm_list, xpoint, xcase)
    this%eq = ES

    call energy(W_mag, W_kin)
    
!    call integrals(sim%fields%node_list, sim%fields%element_list,                                                         &
!        this%eq%R_axis, this%eq%Z_axis, this%eq%psi_axis, this%eq%R_xpoint, this%eq%Z_xpoint,       &
!        this%eq%psi_xpoint, this%eq%psi_bnd, amin, Bgeo, current_t(index_now), beta_p_t(index_now), &
!        beta_t_t(index_now), beta_n_t(index_now), density_tot, density_in_t(index_now),             &
!        density_out_t(index_now), pressure_tot, pressure_in_t(index_now),                           &
!        pressure_out_t(index_now), heat_src_in_t(index_now), heat_src_out_t(index_now),             &
!        part_src_in_t(index_now), part_src_out_t(index_now))

    R_axis_t(index_now)   = this%eq%R_axis
    Z_axis_t(index_now)   = this%eq%Z_axis
    psi_axis_t(index_now) = this%eq%psi_axis

    xtime(index_now)              = t_now
    energies(1:n_tor,1,index_now) = W_mag(1:n_tor)
    energies(1:n_tor,2,index_now) = W_kin(1:n_tor)

    mindelta = minval(deltas); maxdelta = maxval(deltas);
    
    ! --- Output some information about the current timestep
    130 format(1x,a,i5.5,a,es10.3,a)
    131 format(1x,a,2(2(es10.2,' ...',es10.2,',')))
    132 format(1x,'-------------------------------------------------------------------')
    133 format(1x,a,2(es10.2,' at ',i10,','))
    write(*,*)
    write(*,132)
    write(*,130) 'After step ', index_now, ' (t_now=', t_now, '):'
    write(*,132)
    write(*,133) 'min,max deltas  =', mindelta, minloc(deltas), maxdelta, maxloc(deltas)
    write(*,131) 'W_mag,_kin      =', W_mag(1), W_mag(n_tor), W_kin(1), W_kin(n_tor)
    
    Growth_mag  = 0.d0; Growth_kin  = 0.d0; Growth_mag0 = 0.d0; Growth_kin0 = 0.d0

    if (index_now > index_start+1) then
      if (n_tor .gt. 1) then
        Growth_mag  = 0.5d0*log(abs(energies(n_tor,1,index_now)/energies(n_tor,1,index_now-1)))/ tstep
        Growth_kin  = 0.5d0*log(abs(energies(n_tor,2,index_now)/energies(n_tor,2,index_now-1)))/ tstep
      else
        Growth_mag  = 0.d0
        Growth_kin  = 0.d0
      endif
      if (linear_run) then
        Growth_mag0 = 0.d0
        Growth_kin0 = 0.d0
      else
        Growth_mag0 = 0.5d0*log(abs(energies(1,1,index_now)/energies(1,1,index_now-1)))/ tstep
        Growth_kin0 = 0.5d0*log(abs(energies(1,2,index_now)/energies(1,2,index_now-1)))/ tstep
      endif
      write(*,131) 'Growth_mag,_kin =', Growth_mag0, Growth_mag, Growth_kin0, Growth_kin
    endif
    write(*,132)
    write(*,*)
  
  endif ! myid = 0

  call int3d_new(sim%my_id, sim%fields%node_list, sim%fields%element_list, bnd_node_list, bnd_elm_list, exprs_all_int, res, 1)

  if (sim%my_id .eq. 0 ) then
    ! --- Output energies and growth_rates to text files during the code run
    call write_live_data(index_now)
    call write_live_data_vacuum(index_now)
  endif

  call clck_time_barrier(t1)
  call clck_ldiff(t0,t1,tsecond)

  if (sim%my_id .eq. 0) write(*,FMT_TIMING) sim%my_id, '#  Elapsed time Diagnostics :',tsecond
  
  !---------------------------------------------------------timing
  if ( this%istep == 1) then
    call r3_info_print (-3, -2, 'ITERATION    1')
  else
    call r3_info_print (this%istep, -2, 'ITERATION')
  endif
  write(itlabel,'(I8)') this%istep
  call tr_print_memsize("AfterIter"//itlabel)
  
  ! --- Write a restart file every nout timesteps
  if ( (sim%my_id == 0) .and. (mod(index_now,nout) == 0) ) then
    write(fileout,'(A5,i5.5)') 'jorek',index_now
    call export_restart(sim%fields%node_list, sim%fields%element_list, fileout)
  endif
  
  ! --- Exit the code if NaNs are detected.
  if ( allocated(deltas) ) then
    sum_deltas = sum(deltas)
    if ( sum_deltas /= sum_deltas ) then
      write(*,*)
      write(*,*) '>>>>> NaNs DETECTED: EXITING THE CODE <<<<<'
      write(*,*)
      sim%stop_now = .true.
    end if
  end if
  
  call clck_time_barrier(t1)
  call clck_ldiff(t_itstart,t1,tsecond)
  if (sim%my_id .eq. 0) write(*,FMT_TIMING) sim%my_id, '# Elapsed time ITERATION :',tsecond

  this%istep = this%istep + 1

  ! --- Exit the code on end of timestepping spec
  if (get_tstep_n(this%istep) .eq. 0.d0) then
    write(*,*) "Last timestep executed, stopping"
    sim%stop_now = .true.
  end if

  ! Write a restart file on code exit
  if (sim%stop_now .and. sim%my_id .eq. 0) then
    call export_restart(sim%fields%node_list, sim%fields%element_list, 'jorek_restart')
  end if

  select type (fields => sim%fields)
  type is (jorek_fields_interp_linear)
    fields%time_prev = sim%time
    fields%time_now  = sim%time + dt
  end select
end subroutine do_jorek_timestep


!> Calculate the timestep size (in jorek units) at step i by looping through the nstep_n array
function get_tstep_n(i) result(dt)
  use phys_module
  integer, intent(in) :: i
  integer :: j, n_steps_accum
  real*8  :: dt

  n_steps_accum = 0
  dt = 0.d0
  do j=1,size(tstep_n,1)
    n_steps_accum = n_steps_accum + nstep_n(j)
    if (n_steps_accum .ge. i) then
      dt = tstep_n(j)
      if (dt .eq. 0.d0) then
        write(*,*) "ERROR: No time stepping detected. Exiting."
        stop
      endif
      return
    else
      if (j .eq. size(tstep_n,1)) then
        write(*,*) "WARNING: no next timestep detected, returning 0"
      end if
    end if
  end do
end function get_tstep_n
end module mod_jorek_timestepping
