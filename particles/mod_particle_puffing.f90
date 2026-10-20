!> Module for puffing gas into the plasma
!> This module will create new particles at the locations where gas valves will be.

module mod_particle_puffing
  use mod_io_actions, only: io_action
  use mod_sampling
  use mod_particle_types
  use constants, only: TWOPI, K_BOLTZ, ATOMIC_MASS_UNIT
  use mod_rng, only: type_rng, setup_shared_rngs
  use mod_boundary, only: wall_normal_vector
  use mod_atomic_elements 
  use mod_particle_sim
  use mod_event
  use mod_find_rz_nearby, only: find_rz_nearby
  use mod_particle_allocation, only: calc_n_particles_per_mpi
  use phys_module, only: type_valve, valves, part_group_configs, type_puff_ctrl, n_puff_segment_max 
  use mod_particle_group_id, only: matching_part_config_indices, matching_sim_groups_indices
  use mod_particle_create, only: type_part_create_scheme
  use mpi
  
  implicit none

  real*8  :: supers_ratio_puff_default = 1.d-4                 !< if none of the puff_ctrl(i)%supers_..._puff options are set, supers_to_puff will be calculated
                                                               !< as supers_ratio_puff_default * part_group_config(i)%n_particles
                                                               !< In this case this default value overrides the value from preset_parameters.f90

  private
  public  :: particle_puffing

  ! Extend type
  type, extends(io_action) :: particle_puffing
    
    !> Set by user
    type(type_valve)     :: puff_valve                         !< determines the location of the puffing
    integer              :: valve_num                          !< the index of the valve used for this puffing event
    type(type_puff_ctrl) :: puff_ctrl                          !< the controller for this puffing action
    integer              :: target_group                       !< target particle group

    !> Internal
    type(type_part_create_scheme) :: create_scheme             !< super particle create scheme 
    class(type_rng), dimension(:), allocatable :: rng          !< one RNG per openmp thread
    real*8               :: vector_normal(3)                   !< vector_normal of puff_valve
    integer              :: val_i_elm                          !< element index of puff_valve
    real*8               :: val_R, val_Z                       !< R, Z location of puff_valve
    real*8               :: val_s, val_t                       !< s, t location of puff_valve
    real*8, dimension(:), allocatable :: times                 !< timepoints of puffrate (in s)
    real*8, dimension(:), allocatable :: rates                 !< timepoints of puffrate (in s)

    !> Variables required for piecewise linear time dependent puffing control 
    integer              :: current_puff_seg = 0   
    integer              :: last_puff_seg = n_puff_segment_max !< the segment number after the last defined puffing keyframe
    
  contains
    procedure :: do => do_particle_puffing
  end type particle_puffing

  interface particle_puffing
    module procedure new_particle_puffing
  end interface particle_puffing
contains


!> validity checks for user input settings in puff_ctrl
!> and set puffing settings based on these settings
subroutine initialize_settings_from_puff_ctrl(sim, group_num, valve_num, new)
  use mod_particle_create, only: part_create_scheme
  use profiles, only: readProf
  use tr_module
  use mpi_mod
  
  implicit none

  type(particle_sim),     intent(in)     :: sim
  integer,                intent(in)     :: group_num
  integer,                intent(in)     :: valve_num
  type(particle_puffing), intent(inout)  :: new

  integer                                :: i, config_num, file_len, ierr
  character(len=1000)                    :: identifier

  logical :: from_namelist, from_file

  !> variables for piecewise linear
  integer :: times_counter = 0 
  integer :: rates_counter = 0

  config_num = matching_part_config_indices(group_num)

  ! determine the create scheme
  write(identifier,"(A,I2,A,I2,A)") "part_group_config(",config_num,")%puff_ctrl(",valve_num,")"
  new%create_scheme = part_create_scheme(new%puff_ctrl%supers_num_puff, new%puff_ctrl%supers_weight_puff, &
    new%puff_ctrl%supers_ratio_puff, sim%groups(group_num)%n_particles, &
    default=supers_ratio_puff_default, my_id=sim%my_id, identifier=identifier)
  
  !> specific to piecewise linear puff_ctrl ----------------------------

  !> validity checks
  ! check whether namelist values are set
  rates_counter = count(new%puff_ctrl%rates /= -1)
  times_counter = count(new%puff_ctrl%times /= -1)
  from_namelist = rates_counter > 0 .or. times_counter > 0

  ! check whether puff rate file is set
  from_file = trim(new%puff_ctrl%from_file) /= "none"

  !error if both are set simultaneously
  if (from_namelist .and. from_file) then
    if (sim%my_id == 0) write(*,"(A,I2,A,I2,A)") "ERROR: both namelist input (%rates and %times) and %from_file specified for part_group_configs(",config_num,")%puff_ctrl(",valve_num,"). Please use only 1 way to specify input puff rate over time"
    stop
  endif

  !using namelist
  if (from_namelist) then
    if (rates_counter == 0) then
      if (sim%my_id == 0) write(*,"(A,I2,A,I2,A)") "ERROR: you set %times but not %rates for part_group_configs(",config_num,")%puff_ctrl(",valve_num,")."
      stop
    endif

    if (times_counter == 0) then
      if (sim%my_id == 0) write(*,"(A,I2,A,I2,A)") "ERROR: you set %rates but not %times for part_group_configs(",config_num,")%puff_ctrl(",valve_num,")."
      stop
    endif

    if (times_counter /= rates_counter) then
      if (sim%my_id == 0) then
        write(*,"(A,I2,A,I2,A)") "ERROR: Mismatch in part_group_configs(",config_num,")%puff_ctrl(",valve_num,"),"
        write(*,*)               "  between number of entries of %times and %rates. "
      endif
      stop
    endif

    do i=1, times_counter
      if (new%puff_ctrl%rates(i) < 0) then
        if (sim%my_id == 0) write(*,"(A,I2,A,I2,A,I2,A)") "ERROR: part_group_configs(",config_num,")%puff_ctrl(",valve_num,")%rates(",i,") < 0."
        stop
      endif
      if (new%puff_ctrl%times(i) < 0) then
        if (sim%my_id == 0) write(*,"(A,I2,A,I2,A,I2,A)") "ERROR: part_group_configs(",config_num,")%puff_ctrl(",valve_num,")%times(",i,") < 0."
        stop
      endif

      !> check that the set times are increasing
      if (i > 1) then
        if (new%puff_ctrl%times(i-1) >= new%puff_ctrl%times(i)) then
          if (sim%my_id == 0) write(*,"(A,I2,A,I2,A)") "ERROR: inputs in %times for part_group_config(",config_num,")%puff_ctrl(",valve_num,") must be strictly increasing"
          stop
        endif 
      endif
    enddo

    allocate(new%times(times_counter))
    new%times=new%puff_ctrl%times(1:times_counter)
    allocate(new%rates(rates_counter))
    new%rates=new%puff_ctrl%rates(1:rates_counter)
  end if

  !using puff rate from file
  if (from_file) then
    ! read and check on main mpi thread
    if (sim%my_id==0) then
      !loading file
      call readProf(new%times,new%rates,file_len,new%puff_ctrl%from_file)

      !sanity checks
      do i=1, file_len
        if (new%rates(i) < 0) then
          write(*,"(A,I2,A,I2,3A,I5,A)") "ERROR: negative rate in part_group_configs(",config_num,")%puff_ctrl(",valve_num,')%from_file="',trim(new%puff_ctrl%from_file),'" (position ',i,")."
          call MPI_Abort(MPI_COMM_WORLD, 51, IERR)
        endif
        if (new%times(i) < 0) then
          write(*,"(A,I2,A,I2,3A,I5,A)") "ERROR: negative time in part_group_configs(",config_num,")%puff_ctrl(",valve_num,')%from_file="',trim(new%puff_ctrl%from_file),'" (position ',i,")."
          call MPI_Abort(MPI_COMM_WORLD, 52, IERR)
        endif

        !> check that the set times are increaseing
        if (i > 1) then
          if (new%times(i-1) >= new%times(i)) then
            write(*,"(A,I2,A,I2,3A,I5,A)") "ERROR: puff rate time inputs in part_group_configs(",config_num,")%puff_ctrl(",valve_num,')%from_file="',trim(new%puff_ctrl%from_file),'" must be strictly increasing (this is violated at position ',i,")."
            call MPI_Abort(MPI_COMM_WORLD, 53, IERR)
          endif 
        endif
      enddo
    end if

    !broadcast to all other ranks
    call MPI_BCAST(file_len,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    if ( sim%my_id /= 0 ) then
      call tr_allocate(new%times,1,file_len,"new%times",CAT_UNKNOWN)
      call tr_allocate(new%rates,1,file_len,"new%rates",CAT_UNKNOWN)
    end if
    call MPI_BCAST(new%times,file_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(new%rates,file_len,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    
  end if

  !> determine the current puffing segment and the last puffing segment 
  do i=1, size(new%times)
    !> find current puffing segment
    if (sim%time > new%times(i)) new%current_puff_seg = i
  enddo
  new%last_puff_seg=size(new%times)
  ! -------------------------------------------------
end subroutine


!> 1. checks whether a specific valve has been initialized properly, 
!>    i.e. whether it has a valid type, and if it is within the domain
!> 2. if initialized properly, determine its location in the domain
!> and get the wall vector normal 
subroutine initialize_puff_valve(sim, valve_num, new)
  implicit none
  type(particle_sim),     intent(in)             :: sim
  integer,                intent(in)             :: valve_num
  type(particle_puffing), intent(inout)          :: new
  type(type_valve)                               :: valve
  integer                                        :: ifail

  valve = valves(valve_num)

  select case (trim(valve%type))
    !> Check whether the user forgot to set the valve type
    case ('none')
      if(sim%my_id == 0) write(*,"(A,I2,A,I2,A)") "ERROR: Valve, ", valve_num, " has been selected for puffing but no valve()%type has been set"
      stop

    !> For valid valve types, check that it is within the domain ------------
    !> (add a case here when implementing new valve types)
    case ('circ')
      call find_RZ(sim%fields%node_list, sim%fields%element_list, valve%R_valve_loc, valve%Z_valve_loc, new%val_R, new%val_Z, &
      new%val_i_elm, new%val_s, new%val_t ,ifail)

    case ('poly')
      call find_RZ(sim%fields%node_list, sim%fields%element_list, sum(valve%poly_R(1:2))/2.d0, sum(valve%poly_Z(1:2))/2.d0, new%val_R, new%val_Z, &
      new%val_i_elm, new%val_s, new%val_t ,ifail)

    !> Types not listed as a case above are unsupported ---------------------
    case default
      write(*,*) "ERROR: the valve type '", valve%type, "' is currently not supported."
      stop
  end select

  if (ifail .ne. 0) then
    if (sim%my_id .eq. 0) then
      write(*,"(A,I2,A)") "WARNING: The location set for valve ", valve_num, " could"
      write(*,*)          "  not be found, maybe it was placed outside of the grid?"
    endif
    stop
  else
    new%vector_normal = wall_normal_vector(sim%fields%node_list, sim%fields%element_list, new%val_i_elm, new%val_s, new%val_t)
  end if

end subroutine 

function new_particle_puffing(sim, target_group, valve_num, rng, seed) result(new)

  use mod_pcg32_rng,   only: pcg32_rng
  use mod_random_seed, only: random_seed
  
  type(particle_puffing)    :: new

  type(particle_sim), intent(in)  :: sim           
  integer, intent(in)             :: target_group
  integer, intent(in)             :: valve_num           !< the valve number to use for this puffing
    
  class(type_rng), intent(in), optional :: rng  !< random number generator to use (deafult PCG32)
  integer, intent(in), optional         :: seed !< Seed for the RNG (default random_seed() on my_id + bcast)
  integer                               :: my_seed, i, config_num

  config_num = matching_part_config_indices(target_group)
  
  !> check the validity of the given puff valve -----
  call initialize_puff_valve(sim, valve_num, new)

  !> set puff action parameters based on input -----
  new%target_group     = target_group
  new%puff_valve       = valves(valve_num)
  new%valve_num        = valve_num
  new%puff_ctrl        = part_group_configs(config_num)%puff_ctrl(valve_num)

  !> check the validity and initialize settings from puff_ctrl
  call initialize_settings_from_puff_ctrl(sim, target_group, valve_num, new)

  !> allocate random seed for sampling -----
  if (present(seed)) then
    my_seed = seed
  else
    my_seed = random_seed()
  end if
  if (present(rng)) then
    call setup_shared_rngs(n_dim=3, seed=my_seed, rng_type=rng, rngs=new%rng)
  else
    ! default to pcg32_rng for reflection
    call setup_shared_rngs(n_dim=3, seed=my_seed, rng_type=pcg32_rng(), rngs=new%rng)
  end if

end function new_particle_puffing


!< linearly interpolates the puff rate between two determined points (t0, y0), (t1, y1)
!< y = y0 + [(y1 - y0)/(t1 - t0)] * (t - t0)
subroutine calc_puff_rate_linear(this, time, puff_rate_0, puff_rate_1, puff_rate)

  implicit none
  class(particle_puffing), intent(inout) :: this
  real*8, intent(in)                     :: time          !< current simulation time (t)
  real*8, intent(inout)                  :: puff_rate_0   !< the left set value of the linear function (y0)
  real*8, intent(inout)                  :: puff_rate_1   !< the right set value of the linear function (y1)
  real*8                                 :: puff_time_0   !< t0
  real*8                                 :: puff_time_1   !< t1
  real*8, intent(inout)                  :: puff_rate     !< the puff rate at time t (y)

  !> check if the simulation has advanced to the next puffing segment
  ! if some puffing segments are very small (e.g. to create a step function), you might need to skip more than 1 puffing segment, 
  ! so keep on checking until t<t_next_seg or the last segment is reached
  ! do ... enddo just does the loop infinitely until exit is called
  do
    if (this%current_puff_seg /= this%last_puff_seg) then
      if (time > this%times(this%current_puff_seg + 1)) then
        this%current_puff_seg = this%current_puff_seg + 1 
      else
        exit
      endif
    else
      exit
    endif
  enddo

  !> set the bounding values of the current puffing segment
  if (this%current_puff_seg == 0) then 
    !> for t < times(1), we keep puff_rate constant at rates(1)
    puff_rate_0 = this%rates(1)
    puff_rate_1 = this%rates(1) 

    !> the times no longer matter in this case but (times_1 - times_0) has to be non zero 
    puff_time_0 = this%times(1) 
    puff_time_1 = this%times(1) + 1 
  else if (this%current_puff_seg == this%last_puff_seg) then
    !> for t > times(last_puff_seg), we keep puff_rate constant at rates(last_puff_seg)
    puff_rate_0 = this%rates(this%last_puff_seg)
    puff_rate_1 = this%rates(this%last_puff_seg) 

    !> the times no longer matter in this case but (times_1 - times_0) has to be non zero 
    puff_time_0 = this%times(this%last_puff_seg) 
    puff_time_1 = this%times(this%last_puff_seg) + 1 
  else
    !> in general, times and rates are the defined values bounding the segment
    puff_rate_0 = this%rates(this%current_puff_seg)
    puff_rate_1 = this%rates(this%current_puff_seg + 1)
    puff_time_0 = this%times(this%current_puff_seg)
    puff_time_1 = this%times(this%current_puff_seg + 1)
  endif

  puff_rate = puff_rate_0 + (puff_rate_1 - puff_rate_0) / (puff_time_1 - puff_time_0) * (time - puff_time_0)

end subroutine calc_puff_rate_linear

!> Actually puff gas
subroutine do_particle_puffing(this,sim, ev)
  use mpi_mod
  use phys_module, only: tstep, central_mass, central_density
  use constants, only: ATOMIC_MASS_UNIT, MU_ZERO
  ! use mod_particle_create, only: free_particle_indices ! TODO
  ! !$ use omp_lib

  class(particle_puffing) , intent(inout) :: this
  type(particle_sim), intent(inout)       :: sim
  type(event), intent(inout), optional    :: ev 

  integer :: ierr,i_scalar, n_free, j, k, n_group, i_elm, i_elm_new, ifail, i_p, to_puff, supers_to_puff, supers_to_puff_local, i_rng
  logical, allocatable, dimension(:) :: is_free ! TODO rm
  integer, allocatable, dimension(:) :: i_free
  real*8  :: tstep_fluid_si, c, R, Z, phi, s, t
  real*8  :: R_new, Z_new, s_new, t_new, r_valve, theta
  real*8  :: u(5)
  real*8  :: puff_rate, puff_rate_local !< possibly time dependent fueling rate, fueling rate/n_mpi

  integer ::    puffed_this_step_local, all_puffed_this_step
  real*8  ::    puff_weight_local, all_puff_weight

  !> Variables for specific puffing patterns ----
  !> for piecewise linear time dependent puffing
  real*8  :: puff_rate_0, puff_rate_1 

  if (sim%my_id == 0) write(*,'(A,A,A,I1,A)') "------ For Group: ", sim%groups(this%target_group)%id, ", Valve: ", this%valve_num, " ------"


  !> Set R, Z, s, t, and i_elm to the location of the center of the valve -----
  R = this%val_R
  Z = this%val_Z
  s = this%val_s
  t = this%val_t
  i_elm = this%val_i_elm

  !> get SI time ----------------
  tstep_fluid_si = tstep*sqrt((MU_ZERO * CENTRAL_MASS * ATOMIC_MASS_UNIT * CENTRAL_DENSITY * 1.d20))

  ! TODO rm
  !============== Finding free particles !< make into a function?
  allocate(is_free(size(sim%groups(this%target_group)%particles,1))) 
  !$omp parallel do default(none) shared(sim, this, n_free, i_free, is_free) &
  !$omp private(j) schedule(dynamic, 100)
  do j=1,size(sim%groups(this%target_group)%particles,1) !sim%groups(this%target_group)%particles
    is_free(j) = sim%groups(this%target_group)%particles(j)%i_elm .le. 0  !< array T/F is particle is free
  end do
  !$omp end parallel do
  !$omp barrier
  n_free = count(is_free)
  allocate(i_free(n_free))
  k = 1
  do j=1,size(is_free,1)
    if (is_free(j)) then
      i_free(k) = j !< i_free(k) has index of free particle in  sim%groups(this%target_group)%particles(j)
      k = k+1
      !if (sim%my_id .eq. 0) write(*,*) "Adding to the list number: ", j
    end if
  end do
  ! ==================
  
  ! The current set up may only work for puffing Hydrogen
  ! Assuming the incoming gas at T=300C and a diatomic gas
  ! c is the sound speed of the puffed gas, which is given as c = sqrt(gamma*kB*T/m) for an ideal gas
  ! A factor of 2 is used for the mass as we are puffing diatomic molecules
  c = sqrt((7.d0/5.d0)*(300.d0+273.d0)*K_BOLTZ/(2.d0*sim%groups(this%target_group)%mass*ATOMIC_MASS_UNIT))

  !> calculate time dependent puffing rate -----------------------
  !> piecewise linear approach 
  call calc_puff_rate_linear(this, sim%time, puff_rate_0, puff_rate_1, puff_rate)
  puff_rate_local = puff_rate / sim%n_mpi

  !> calculate how many superparticles to initiate ---------------
  supers_to_puff = this%create_scheme%supers_to_create(sim%my_id, tstep_fluid_si * puff_rate)
  supers_to_puff_local = calc_n_particles_per_mpi(supers_to_puff, sim%n_mpi, sim%my_id)

  ! TODO: use function
  ! ! determine indices of free particles
  ! call free_particle_indices(sim%groups(this%target_group)%particles, n_free, i_free, n_needed=supers_to_puff_local)  

  to_puff = supers_to_puff_local
  if (to_puff .ge. n_free) then
    to_puff = n_free
    write(*,"(A, I3, A)") "WARNING: On mpi proc ", sim%my_id, ", not enough free particles available to puff the requested amount."
    write(*,"(A, I12, A, I12)") "  Requested: ", supers_to_puff_local, " | Free (actually puffed): ", n_free
  end if

  !> write out puffing details -------------------------------------
  if (sim%my_id .eq. 0) then
    write(*,"(2X,A12)") "Set-up:     "
    write(*,"(4X,A18, ' = ', I12)")        "puff segment      ", this%current_puff_seg
    write(*,"(4X,A18, ' = ', G13.6)")      "puff_rate_0       " , puff_rate_0
    write(*,"(4X,A18, ' = ', G13.6)")      "puff_rate_1       " , puff_rate_1
    write(*,"(4X,A18, ' = ', G13.6,A)")    "puff_rate         "   , puff_rate, " atoms/s"
    if (trim(this%create_scheme%scheme) == "num") write(*,"(4X,A18, ' = ', I12)") "supers_num_puff   ", this%puff_ctrl%supers_num_puff
    if (trim(this%create_scheme%scheme) == "weight") write(*,"(4X,A18, ' = ', G13.6)") "supers_weight_puff", this%puff_ctrl%supers_weight_puff
    if (trim(this%create_scheme%scheme) == "ratio") write(*,"(4X,A18, ' = ', G13.6)") "supers_ratio_puff ", this%puff_ctrl%supers_ratio_puff
    write(*,*) ""
  endif

  if(puff_rate <= 0) then
    if (sim%my_id .eq. 0) write(*,"(2X,A)") "puff_rate=0 so nothing to be done here. Returning"
    return
  end if
      
  !> Puffing loop ----------------------------------------
  puffed_this_step_local = 0
  puff_weight_local      = 0.d0
  select type (pa => sim%groups(this%target_group)%particles)
  type is (particle_kinetic_leapfrog)
  !> To do: Parallellize loop 
  !> This loop initializes the to be puffed particles and places then in the computational domain.
  !> It counts the amount of marker particles and weight that was initialized.
  !> This loop can be OMP parallel. there was a small bug in the OMP implementation.
  !> The almost done loop is written below.
  !> reduction of puffed_this_step_local,puff_weight_local for diagnostics
  !>
  ! #ifdef __GFORTRAN__
  !  !$omp parallel do default(shared) & ! workaround for Error: �__vtab_mod_pcg32_rng_Pcg32_rng� not specified in enclosing �parallel�
  ! #else
  ! !$omp parallel do default(shared) &
  ! !$omp schedule(dynamic,10) &
  ! !$omp shared(sim, pa, this,i_free,c,                      &
  ! !$omp to_puff,supers_to_puff_local, tstep_fluid_si,puff_rate_local )                        &
  ! !$omp private(i_p, i_rng, j,k,u , R,Z,s,t,R_new,Z_new,s_new,t_new,     &
  ! !$omp  i_elm,i_elm_new,r_valve, theta,                                         &
  ! !$omp ifail)                                                                    &
  ! !$omp reduction(+:puffed_this_step_local,puff_weight_local)
    do j = 1, to_puff
      i_p = i_free(j)

      !> keep doing the following until a valid position is found
      do 
        !  !$ i_rng = omp_get_thread_num()+1
        call this%rng(1)%next(u) !rng(1)

        !> additional valve type specific behaviour should be added in this if block
        if (this%puff_valve%type == "circ") then
          r_valve = this%puff_valve%r_valve*sample_piecewise_linear(2, [0.d0, 1.d0], [1.d0, 0.d0], u(1))
          theta = TWOPI * u(2)
          R_new = this%puff_valve%R_valve_loc + r_valve * cos(theta)
          Z_new = this%puff_valve%Z_valve_loc + r_valve * sin(theta)
        else if (this%puff_valve%type == "poly") then
          s = u(1)
          t = u(2)
          R_new = this%puff_valve%poly_R(1)*(1.d0-s)*(1.d0-t) + this%puff_valve%poly_R(2)*s*(1.d0-t) + this%puff_valve%poly_R(3)*(1.d0-s)*t + this%puff_valve%poly_R(4)*s*t
          Z_new = this%puff_valve%poly_Z(1)*(1.d0-s)*(1.d0-t) + this%puff_valve%poly_Z(2)*s*(1.d0-t) + this%puff_valve%poly_Z(3)*(1.d0-s)*t + this%puff_valve%poly_Z(4)*s*t
        endif
    
        call find_RZ_nearby(sim%fields%node_list, sim%fields%element_list, R, Z, s, t, i_elm, &
        R_new, Z_new, s_new, t_new, i_elm_new, ifail)
        if (ifail .ge. 0) exit
      end do
      R     = R_new
      Z     = Z_new
      s     = s_new
      t     = t_new
      i_elm = i_elm_new
      if (this%puff_valve%phi .lt. 0.d0) then
        pa(i_p)%x(3) = TWOPI*u(3)
      else 
        pa(i_p)%x(3) = this%puff_valve%phi 
      end if
      pa(i_p)%x(1:2)  = [R, Z]
      pa(i_p)%st(1:2) = [s, t]
      pa(i_p)%i_elm   = i_elm
      pa(i_p)%weight  = real(1.d0/to_puff) * tstep_fluid_si * puff_rate_local
   
      pa(i_p)%v       = c * sample_cosine(u(4:5), this%vector_normal)   
      pa(i_p)%q       = 0_1
      if (sim%groups(this%target_group)%particles(i_p)%weight  .le. 1.d-2) then ! if the weight is too low. 
        sim%groups(this%target_group)%particles(i_p)%i_elm = 0
        cycle       
      end if
    puffed_this_step_local = puffed_this_step_local+1
    puff_weight_local      = puff_weight_local + pa(i_p)%weight 
    end do
  ! !$omp end parallel do  
  class default
    write(*,*) 'Particle type not implemented for gas fueling.'
    stop
  end select

! puffed_this_step_local
call MPI_REDUCE(puffed_this_step_local,all_puffed_this_step,1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)   
call MPI_REDUCE(puff_weight_local,all_puff_weight,1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)  
if (sim%my_id .eq. 0) then
  write(*,"(2X,A12)") "Puffed:     "
  write(*,"(4X,A18, ' = ', I12)")    "Superparticles    ", all_puffed_this_step
  write(*,"(4X,A18, ' = ', E13.6)")  "Total Weight      ", all_puff_weight
endif
  
  
end subroutine do_particle_puffing


end module