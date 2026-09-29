!> Module for particle to particle and fluid to particle interactions at the wall.
!> Documentation on the wiki can be found at 
!> https://jorek.eu/wiki/doku.php?id=particles:wall_actions
!> 
!> The main object is the wall_action, in which the wall_action%type is a string
!> determining the type of interaction. Examples of wall_actions are a plasma fluid 
!> species sputtering one particle group (e.g. D plasma sputtering W impurities),
!> or particles from a particular group (e.g. N) reflecting against the wall
!> Every such interaction needs it's own object, and internally the right routines
!> are then called when wall_action%do(sim) is called.
!> 
!> The currently implemented interaction types are: 
!> "self sputter" (e.g. W -> W), "fluid sputter" (e.g. fluid D+ -> W), "reflection" 
!> (e.g. kinetic D -> D), "wall recomb" (e.g. fluid D+ -> D) and "pump" (e.g. kinetic
!> D -> D with lowered weight).
!> "other sputter" (e.g. kinetic N -> W) is not implemented as of yet although some 
!> preparation is already available in the code.  
!> 
!> Eckstein coefficients are used to determine the yield of the interaction and the 
!> resulting energy of the resulting particles. These yields are automatically loaded
!> from the simulation folder based on the original and target species symbols.
!>
!> Actions have global diagnostics (e.g. total particles of group i reflected off the 
!> wall) printed out in the logfile, and sputter diagnostics also have local vtk
!> diagnostics using mod_edge_elements
!>
!> Wall_actions are gathered into wall_act_groups. There are 3 types of groups,
!> - "fluid sputter to one" is a group of "fluid sputter" type wall_actions with
!>   the same target species (e.g. all sputter into "W01"). The yields are precalculated
!>   and there is one create scheme for the whole group, in which the number of particles
!>   are distributed evenly according to their yields. The create scheme for this group is
!>   set by setting it for exactly one of the wall_act_configs in the group. Because there
!>   is a group per target species, there can be multiple of these groups.
!> - "part2self" contains "pump", "self sputter" and "reflection" type wall_actions (in that order),
!>   it needs no create scheme and it should be run after the evolution loop so that particles hitting
!>   the wall are reflected back before being saved or overwritten. This group exists only once.
!>   These actions can be specified to only happen in a certain polygon of the domain. If a particle 
!>   could undergo multiple actions, only the first one is applied (first "pump", then "self sputter", 
!>   then "reflection", and if more than one action of the same type is set (e.g. overlapping pump 
!>   domains), the first applying action in the order as specified by the user is applied)
!> - "other" contains the rest, which currently only is "wall recomb". They all have their own
!>   create scheme. This group exists only once
!> 
!> Limitations:
!> - The incoming angle of the particle/fluid is not taken into account (it is hardcoded 
!>   to 0). Implementing this correctly would require some estimation of surface roughness.
!> - The sampling from and integrating the fluid integrals is done using mod_edge_elements 
!>   rather than using a bezier FE description like the fluid (same is true for wall projections)
!> - A simplified model is used for the energy of sampled particles from the fluid in 
!>   fluid2part actions
!> - Local sputtering/wall recombination (using the polygons) is currently not possible.
module mod_particle_wall_interaction
  use mod_io_actions, only: io_action
  use mod_sampling
  use mod_particle_types
  use mod_eckstein_y_ye
  use constants
  use mod_rng, only: type_rng, setup_shared_rngs
  use mod_boundary, only: wall_normal_vector
  use mod_interp
  use mod_atomic_elements !< chemical elements
  use mod_particle_sim
  use mod_event
  use mod_particle_allocation, only: calc_n_particles_per_mpi
  use equil_info, only:find_xpoint
  use mod_particle_create, only: part_create_scheme, type_part_create_scheme, create_scheme_is_set
  use mod_edge_elements, only: edge_elements, type_cdf_data
  use mpi
  !$ use omp_lib 
  
  implicit none
   
  private
  public :: wall_act_group, wall_actions_from_config, gcd_wall_acts

  ! action containing the wall interaction information for one origin species to one target species
  type, extends(io_action) :: wall_action
    integer           :: origin_group  !< index specifying which group is undergoing this wall interaction. Either particle group number (sim%groups(target_group)) or fluid group number (if it is a fluid to particle interaction type)
    integer           :: target_group  !< which particle group (sim%groups(target_group)) this wall interaction affects
    character(len=20) :: type = "none" !< type of the wall interaction, namely "self sputter" (e.g. W -> W), "fluid sputter" (e.g. fluid D+ -> W), "other sputter" (e.g. kinetic N -> W), "reflection" (e.g. kinetic D -> D), "wall recomb" (e.g. kinetic D+ -> D) or "pump" (e.g. kinetic D -> D with weight reduction)

    ! polygon in which the reaction is done (currently only implemented for part2self actions)
    logical                           :: only_in_polygon = .false. !< whether to execute this wall_action only on the specified polygon (.true.) or on the full domain (.false.)
    real*8, dimension(:), allocatable :: poly_R                    !< R coordinates of the polygon. Make sure to define the polygon in order (the polygon is defined as linesegments drawn from point 1 to point 2 to point 3, etc.)
    real*8, dimension(:), allocatable :: poly_Z                    !< Z coordinates of the polygon.

    ! internal variables to determine which kind of backend function needs to be called, depending on whether the origin group is fluid or not and the target group is the origin group or not
    logical :: part2self = .false., part2other = .false., fluid2part = .false.

    type(eckstein_sputter_yield)          :: yield  !< eckstein coefficients for the wall interaction yield
    type(eckstein_sputtered_energy_coeff) :: energy !< eckstein coefficients for determining energy of the resulting particle
    type(thompson_dist)                   :: E_dist = thompson_dist(E_b = 8.7d0, n=2) !< produces energies in eV (value for W default)
    logical :: use_thompson = .false. !< Use a thompson distribution for the energy of sputtered particles
    logical :: use_Yn_func  = .false. !< Use Ecksteins interpolating functions instead of interpolating manually
    
    class(type_rng), dimension(:), allocatable :: rng !< one RNG per openmp thread

    integer             :: each_nstep_part = -9999991 !< run this particular action at every i_inner_loop = each_nstep_part
    
    !> when the origin group is a fluid species
    integer             :: fluid_Z         = -999     !< Z of this fluid species (e.g. -2 for D)
    type(edge_elements) :: fluid_yield_integral       !< the yield (of the specified interaction type) integrated over f(v) for this fluid species
    type(type_cdf_data) :: res                        !< data on the cumulative distribution function calculated at the integral which is needed when sampling
    real*8              :: domain_integral            !< [# particles] total weight of all created particles in this wall_action for this timestep
    logical             :: yield_calculated=.false.   !< whether the fluid yield integral has already been separately calculated for this timestep (.true.) or not (.false.)

    type(type_part_create_scheme) :: create_scheme    !< super particles create scheme
    real*8  :: weight_factor = 1.d0 !< additional weight factor of the yield, combines density_fraction for fluids and wall_act_config(:)%weight_factor (e.g. useful to split a single plasma fluid into D and T neutrals upon wall recombination, or set finite wall absorption)

    ! diagnostics
    logical             :: do_wall_projection=.true. !< whether to do wall projections for this interaction
    type(edge_elements) :: wall_projection           !< diagnostic to keep track of particle- and heat fluxes, sputtering yields, etc. resolved on the wall (1D for 2D simulations, 2D for 3D simulations)
    integer             :: i_step_diag = 0           !< how many steps have been taken between the previous diagnostic output and now
    integer             :: n_step_diag               !< after how many timesteps the wall projection should be saved (as vtk file in the simulation folder)
    real*8              :: last_diag_time = -9.d99   !< Last time of output of diagnostics
    integer             :: n_project_extra           !< Number of extra projection diagnostics for this interaction on top of n_project_general
    integer             :: n_project_tot             !< Total number of wall projections for this diagnostic (i.e. shorthand for n_project_general + n_project_extra)
    integer             :: n_project_part = -1       !< Number of particle projections (as those need to be MPI reduced)
    real*8              :: delta_t                   !< [s] tstep in SI
 
    logical             :: constructed =.false.      !< whether the constructor has been called (this is used as assert in the do action) 
  contains
    procedure :: calc_fluid_yield
    procedure :: do => do_wall_action
    procedure :: load_eckstein_data
    procedure :: update_delta_t
  end type wall_action

  ! wall action groups host similar wall actions that should run together, for instance because they have to communicate with eachother
  ! as is the case for "fluid sputter to one" (fluid sputter type wall actions making one particle species) that need to communicate the 
  ! sputter yields to distribute the superparticles per fluid species
  type :: wall_act_group
    character(len=80)                            :: group_type="none"    !< "fluid sputter to one" for fluids sputtering into one particle species (so type "fluid sputter", and same target_group_id), "part2self" for particle interactions with themself (currently "reflection", "self sputter"), or "other" for other interactions (currently "wall recomb")
    logical                                      :: contains_part2part   !< whether the group contains some particle -> particle (.true.) interaction(s) or only fluid -> particle (.false.) interactions
    character(len=3)                             :: target_id="non"      !< id of the particle group being made by this wall_act_group (if it is "sputter many to one")
    type(type_part_create_scheme)                :: create_scheme        !< shared super particles create scheme for the whole group (for type "fluid sputter to one")
    type(wall_action), allocatable, dimension(:) :: wall_actions         !< array with wall actions for this group

    logical                                      :: constructed =.false. !< whether the constructor has been called (this is used as assert in the do action) 
  contains
    procedure :: do => do_wall_act_group
  end type wall_act_group
  
  !> indices of different diagnostics in the global diagnostics array which is used for the output file
  !> number of super particles is intentionally stored in a real, to easily handle all diagnostics simultaneously (in omp reductions and in MPI_reduce)
  integer, parameter :: n_global_diagnostics=10, i_wall_part_in=1, i_wall_flux_in=2, i_wall_heat_in=3, i_wall_part_out=4, i_wall_flux_out=5, i_wall_heat_out=6, i_wall_flux_refl=7, i_wall_heat_refl=8, i_removed=9, i_super_killed=10
  
  integer, parameter :: n_project_general=4                    !< number of general projections (on top of the number of interaction type specific interactions)

  real*8, parameter  :: supers_ratio_wall_default = 5.d-4      !< if none of the wall_act_configs(i)%supers_..._wall options are set, supers_to_create will be calculated
                                                               !< as supers_ratio_wall_default * part_group_config(this%target_group)%n_particles
                                                               !< In this case this default value overrides the value from preset_parameters.f90
  
  integer            :: gcd_wall_acts = -999999                !< greatest common divisor of the %each_nstep_part of all wall_actions
contains

!> Constructor for the particle_sputter type, setting the io_action parameters and sputtering parameters.
subroutine construct_wall_action(this, sim, origin_group, config, edge_element_template, origin_is_fluid, fluid_Z, fluid_density_fraction, each_nstep_part, filename, basename, decimal_digits, fractional_digits, rng, input_identifier)
  use mod_pcg32_rng, only: pcg32_rng
  use mod_random_seed, only: random_seed
  use phys_module, only: nout_projection, n_fluid_groups_max, n_part_groups, n_part_groups_max, fluid_configs, type_wall_act_config
  use mod_particle_group_id, only: matching_sim_groups_indices
  use mod_particle_sim, only: group_num_from_id

  implicit none
  type(wall_action),             intent(inout) :: this         !< the new wall_action object. Inout because it may need some settings already
  type(particle_sim),            intent(inout) :: sim
  integer,                       intent(in) :: origin_group    !< config index specifying which group is undergoing this wall interaction. Either particle group number or fluid group number (if it is a fluid to particle interaction type)
  type(type_wall_act_config),    intent(in) :: config          !< wall_act_config to make the wall_action from
  type(edge_elements),           intent(in) :: edge_element_template !< a prepared set of edge elements
  logical,                       intent(in) :: origin_is_fluid       !< whether the origin group is a fluid group (.true.) or a particle group (.false.)
  integer,                       intent(in), optional :: fluid_Z                !< Z of this fluid species (e.g. -2 for D)
  real*8,                        intent(in), optional :: fluid_density_fraction !< fraction of the plasma density of this specific fluid. The density fractions of all used fluid configs should add up to 1
  integer,                       intent(in), optional :: each_nstep_part        !< run this particular particle-particle action at every i_inner_loop = each_nstep_part 
  character(len=*),              intent(in), optional :: filename               !< where to save the diagnostics
  character(len=*),              intent(in), optional :: basename
  integer,                       intent(in), optional :: decimal_digits
  integer,                       intent(in), optional :: fractional_digits
  class(type_rng),               intent(in), optional :: rng                    !< random-number generator to use (default PCG32)
  character(len=*),              intent(in), optional :: input_identifier       !< extra message on stops, to determine which construct_wall_action call had wrong input
 
  character(len=100) :: name, origin_name
  integer :: my_seed, i, j, target_group_loc, supers_num_loc, n_poly_R, n_poly_Z
  real*8  :: supers_weight_loc, supers_ratio_loc, n_particles
  character(len=14), dimension(:), allocatable :: extra_proj_scalar_names !< additional scalar names on top of the normal ones
  character(len=1000) :: msg !< error message
  character(len=1000) :: identifier
  logical :: creates_particles

  ! setting the identifier (it is optional for the sake of using construct_wall_action directly from inside a program rather than through the namelist)
  if(present(input_identifier)) then
    identifier = input_identifier
  else
    identifier = ""
  end if

  ! --- determining the interaction type
  this%type = trim(config%type)

  select case(trim(this%type))
  case("self sputter")
    this%part2self = .true.
  case("fluid sputter")
    this%fluid2part = .true.
  case("other sputter")
    this%part2other = .true.
    write(msg,"(A)") 'type "other sputter" is still to be implemented'
    call wrong_input(msg, sim%my_id, identifier)
  case("reflection")
    this%part2self = .true.
  case("wall recomb")
    this%fluid2part = .true.
  case("pump")
    this%part2self=.true.
  case default
    call wrong_interaction_type(this%type, identifier)
  end select
  
  ! --- general checks on input
  
  ! check whether the origin group (particle or fluid group) is compatible with the interaction
  if (this%fluid2part) then ! type suggests that origin is fluid group
    if(.not. origin_is_fluid) then
      write(msg,"(3A)") "interaction type '",trim(this%type),"' cannot be used for particle species"
      call wrong_input(msg, sim%my_id, identifier)
    end if
  else ! type suggests that origin is particle group
    if(origin_is_fluid) then
      write(msg,"(3A)") "interaction type '",trim(this%type),"' cannot be used for fluid species"
      call wrong_input(msg, sim%my_id, identifier)
    end if
  end if

  if (this%fluid2part) then
    creates_particles = .true.
    if(present(fluid_Z)) then
      ! check whether Z is sensical
      if(fluid_Z == -999) then
        write(msg,"(A,I2,A)") "it seems like you forgot to set the fluid_configs(",origin_group,")%Z"
        call wrong_input(msg, sim%my_id, identifier)
      end if
      if(fluid_Z < lbound(element_symbols,1) .or. fluid_Z > ubound(element_symbols,1)) then
        write(msg,"(A,3I5)") "fluid Z not in bound (Z/min/max)",fluid_Z,lbound(element_symbols,1),ubound(element_symbols,1)
        call wrong_input(msg, sim%my_id, identifier)
      end if
      this%fluid_Z = fluid_Z
    else ! fluid Z must be present
      write(msg,"(A)") "fluid Z must be specified for fluid type interaction"
      call wrong_input(msg, sim%my_id, identifier)
    end if

    if(origin_group < 1 .or. origin_group > n_fluid_groups_max) then
      write(msg,"(A,2I4)") "fluid origin group is not valid (origin_group/max)",origin_group,n_fluid_groups_max
      call wrong_input(msg, sim%my_id, identifier)
    end if
    this%origin_group = origin_group

    if(present(fluid_density_fraction)) then
      !check whether fluid density is sensical
      if(fluid_density_fraction < -1.d-12) then
        write(msg,"(A,es12.2)") "It seems like you forgot to set the fluid density_fraction, because it has a negative value still: ",fluid_density_fraction
        call wrong_input(msg, sim%my_id, identifier)
      end if
      if(fluid_density_fraction > 1.d0 + 1.d-12) then
        write(msg,"(A,es12.2,A)") "Fluid density_fraction=",fluid_density_fraction,", but cannot be > 1, please set the relative density between 0 and 1"
        call wrong_input(msg, sim%my_id, identifier)
      end if
    end if

  else ! origin_group is a particle group
    ! check whether it is in bounds for the matching array
    if(origin_group < 1 .or. origin_group > n_part_groups_max) then
      write(msg,"(A,2I4)") "particle config origin group is not valid (origin_group/max)",origin_group,n_part_groups_max
      call wrong_input(msg, sim%my_id, identifier)
    end if
    
    this%origin_group = matching_sim_groups_indices(origin_group)

    ! check if it is in bounds of sim%groups(this%origin)
    if(this%origin_group < 1 .or. this%origin_group > n_part_groups) then
      write(msg,"(A,2I4)") "sim%groups(this%origin_group) is not valid (this%origin_group/max)",this%origin_group,n_part_groups
      call wrong_input(msg, sim%my_id, identifier)
    end if

    if(this%part2self) then
      creates_particles = .false.
    end if
    if (this%part2other) then
      creates_particles = .true.
    end if
  end if

  !checking and setting target group
  !self interactions default to have the same target as origin if the wall_act_config%target_group_id is the unchanged namelist input value "non"
  if(this%part2self .and. config%target_group_id == "non") then
    this%target_group = this%origin_group
  else
    !the id should be specified
    if(config%target_group_id == "non") then
      write(msg,"(A)") "%target_group_id was not set and has to be set (for this interaction type)"
      call wrong_input(msg, sim%my_id, identifier)
    end if
    !get the corresponding group_num
    this%target_group = group_num_from_id(sim,config%target_group_id)
    !check that a match was found
    if(this%target_group == -1) then
      write(msg,"(3A)") "%target_group_id ",config%target_group_id," is not a valid id that is in use in sim%groups(:). Did you spell it correctly and include this particle group in part_groups_in_use? Error"
      call wrong_input(msg, sim%my_id, identifier)
    else if(this%target_group < 1 .or. this%target_group > n_part_groups_max) then     
      !sanity check on this%target_group, we should never end up here, so there's a bug somewhere if you get this print
      write(msg,"(A,2I3,3A)") "%target_group is not valid (target_group/max): ",target_group_loc,n_part_groups_max," This happened for id=",config%target_group_id," Something strange happened"
      call wrong_input(msg, sim%my_id, identifier)
    end if   
  end if

  !checking whether the user set origin group and target group differently while it is a self interaction
  if (this%part2self) then
    call check_self_type(this, sim%my_id, identifier)
  end if

  if(creates_particles) then ! this interaction requires resulting species particles to be created
    ! setting the creation scheme
    n_particles = sim%groups(this%target_group)%n_particles
    this%create_scheme = part_create_scheme(config%supers_num_wall,config%supers_weight_wall,config%supers_ratio_wall,n_particles,supers_ratio_wall_default,sim%my_id,identifier)    
  else ! creates no particles
    ! check that no create scheme was set
    if(create_scheme_is_set(config%supers_num_wall,config%supers_weight_wall,config%supers_ratio_wall)) then
      write(msg,"(3A)") "you cannot set a create scheme (supers_..._wall) for wall_action type '",trim(this%type),"', as it does not create particles. Problem found"
      call wrong_input(msg, sim%my_id, identifier)
    end if
  end if


  ! checking and setting the weight_factor
  if(config%weight_factor > 1.d0 + 1.d-12) then
    if(sim%my_id == 0) write(*,"(2A)") "WARNING: having a weight_factor > 1 is usually undesireable. Make sure you really want this. weight_factor > 1",identifier
  end if
  if(config%weight_factor < - 1.d-12) then
    write(msg,"(A,es12.2)") "you cannot have negative weight_factor: ",config%weight_factor
    call wrong_input(msg, sim%my_id, identifier)
  end if
  if (trim(this%type)=="pump" .and. abs(config%weight_factor - 1.d0) < 1.d-12) then !if pumping but not pumping
    if(sim%my_id == 0) write(*,"(2A)") "WARNING: you set %type='pump', but left %weight_factor=1.d0 (meaning in effect this will be a strange kind of reflection)",identifier 
  end if
  this%weight_factor = config%weight_factor
  if(this%fluid2part) then !in the backend this%weight_factor also has to take into account fluid_density_fraction
    this%weight_factor = this%weight_factor * fluid_density_fraction
  end if

  ! checking and setting polygon settings
  this%only_in_polygon = config%only_in_polygon
  n_poly_R = count(config%poly_R > -1.d98)
  n_poly_Z = count(config%poly_Z > -1.d98)
  if(config%only_in_polygon) then
    if(this%fluid2part) then
      write(msg,"(A)") "%only_in_polygon is not implemented for fluid2part wall_actions"
      call wrong_input(msg, sim%my_id, identifier)
    end if
    if(n_poly_R /= n_poly_Z) then
      write(msg,"(A,I5,A,I5,A)") "you must specify an equal amount of %poly_R (",n_poly_R," specified) and %poly_Z values (",n_poly_Z," specified)"
      call wrong_input(msg, sim%my_id, identifier)
    end if
    if(n_poly_R == 0) then
      write(msg,"(A)") "you set %only_in_polygon=.t. but you did not specify %poly_R and %poly_Z"
      call wrong_input(msg, sim%my_id, identifier)
    end if
    if(n_poly_R < 3) then
      write(msg,"(A,I5,A)") "you need at least 3 %poly_R and %poly_Z points to define a polygon, but you specified only ",n_poly_R," point(s)"
      call wrong_input(msg, sim%my_id, identifier)
    end if
    
    allocate(this%poly_R(n_poly_R), this%poly_Z(n_poly_R))
    
    do i=1,n_poly_R
      if(config%poly_R(i) < -1.d98 .or. config%poly_Z(i) < -1.d98) then
        write(msg,"(A)") "you need to specify %poly_R and %poly_Z points in order without leaving gaps (i.e. you cannot specify poly_R(2) if you did not specify poly_R(1))"
        call wrong_input(msg, sim%my_id, identifier)
      end if
      this%poly_R(i) = config%poly_R(i)
      this%poly_Z(i) = config%poly_Z(i)
    end do
    
  else
    if(max(n_poly_R,n_poly_Z) > 0) then
      write(msg,"(A)") "mixed messages in input. You set %poly_R and/or %poly_Z value(s) while not setting %only_in_polygon=.true."
      call wrong_input(msg, sim%my_id, identifier)
    end if
  end if

  ! --- diagnostics
  if(this%fluid2part) extra_proj_scalar_names = ["n_e           ","T_e           ","cos_alpha     ","Psi_n         ","fluid_flux    ","fluid_heatflux","fluid_yield   "]
  
  ! if there are no extra projections, set the allocatable to 0
  if(.not. allocated(extra_proj_scalar_names)) allocate(extra_proj_scalar_names(0))
  if(this%n_project_part < 0) this%n_project_part = n_project_general

  call this%load_eckstein_data(sim)

  ! initialising the edge_element objects from the template
  if (.not. allocated(edge_element_template%patch(1)%xyz)) then
    write(msg,"(A)") 'Edge element template needs to be prepared, exiting'
    call wrong_input(msg, sim%my_id, identifier)
  end if
  
  this%wall_projection = edge_element_template
  this%fluid_yield_integral = edge_element_template
  ! Clean up the passed edge elements
  do i=1,size(edge_element_template%patch,1)
    if (allocated(edge_element_template%patch(i)%scalars)) then
      deallocate(this%wall_projection%patch(i)%scalars, &
                 this%fluid_yield_integral%patch(i)%scalars)
    end if
    if (allocated(edge_element_template%patch(i)%scalar_names)) then
      deallocate(this%wall_projection%patch(i)%scalar_names, &
                 this%fluid_yield_integral%patch(i)%scalar_names)
    end if
  end do

  ! settings for the diagnostics
  if(this%fluid2part) then
    write(origin_name,"(I2.2)") this%origin_group
    write(name, "(7A)") spaces2underscore(this%type),"_", trim(origin_name), "_to_", sim%groups(this%target_group)%id
  else
    origin_name = sim%groups(this%origin_group)%id
    !part2self:
    write(name, "(5A)") spaces2underscore(this%type),"_", trim(origin_name)
  end if

  if(this%only_in_polygon) then
    write(name,"(3A)") trim(name),"_at_",trim(config%nametag)
  end if
  
  this%n_step_diag = nout_projection
  this%basename = trim(name)//"_"
  if (present(filename)) this%filename = filename
  if (present(basename)) this%basename = basename
  if (present(decimal_digits)) this%decimal_digits = decimal_digits
  if (present(fractional_digits)) this%fractional_digits = fractional_digits
  this%extension = '.vtk'
  this%name = trim(name)
  this%log = .true.

  ! Set up scalars and scalar names for the diagnostic projections
  this%n_project_extra = size(extra_proj_scalar_names, dim=1)
  this%n_project_tot = n_project_general + this%n_project_extra
  do i=1,size(this%wall_projection%patch,1)
    allocate(this%wall_projection%patch(i)%scalars(size(this%wall_projection%patch(i)%st,2), this%n_project_tot))
    this%wall_projection%patch(i)%scalars(:,:) = 0.d0 ! initialising scalars
    
    allocate(this%wall_projection%patch(i)%scalar_names(this%n_project_tot))
    
    ! defining the scalar names
    associate (sn => this%wall_projection%patch(i)%scalar_names)
      sn(1) = "part_flux"
      sn(2) = "part_heatflux"
      sn(3) = "part_promptflux"
      sn(4) = "part_yield"
      do j=1,this%n_project_extra
        sn(n_project_general+j) = trim(extra_proj_scalar_names(j))
      enddo
    end associate
  end do

  ! Allocate scalars for the fluid_yield_integral
  do i=1,size(this%fluid_yield_integral%patch,1)
    allocate(this%fluid_yield_integral%patch(i)%scalars( &
      size(this%fluid_yield_integral%patch(i)%st,2), 1))
    
    this%fluid_yield_integral%patch(i)%scalars = -1
  end do

  ! --- allocate random seed for sampling
  my_seed = random_seed()
  if (present(rng)) then
    call setup_shared_rngs(n_dim=3, seed=my_seed, rng_type=rng, rngs=this%rng)
  else
    ! default to pcg32_rng
    call setup_shared_rngs(n_dim=3, seed=my_seed, rng_type=pcg32_rng(), rngs=this%rng)
  end if

  !sanity checks and seting when to run the action
  if(present(each_nstep_part)) then
    if(this%fluid2part) then
      write(msg,"(A)") "each_nstep_part is only supported for particle-particle wall_actions but you set it for a fluid-particle wall_action"
      call wrong_input(msg, sim%my_id, identifier)
    endif
    if(each_nstep_part /= -9999991) then
      if(each_nstep_part <= 0) then
        write(msg,"(A)") "each_nstep_part <= 0 which is not allowed"
        call wrong_input(msg, sim%my_id, identifier)
      endif

      this%each_nstep_part = each_nstep_part ! setting it
    endif
  endif

  !constructor finished
  this%constructed = .true.
end subroutine construct_wall_action


!> Set up the wall_act_groups array from the configs of the namelist.
!> There is currently 1 "part2self" group created (or none), 1 "other" group (or none), 
!> and as many "fluid sputter to one" groups as there are unique targets for "fluid sputter"
!> actions. (so between 0 and n_part_groups_max) 
!> The create scheme of the "fluid sputter to one" groups can be set by setting exactly one of the 
!> child wall actions of that group
function wall_actions_from_config(sim, edge_element_template) result(wall_act_groups)
  use phys_module, only: part_group_configs, n_part_groups, n_part_groups_max, type_wall_act_config
  use phys_module, only: fluid_configs, n_fluid_groups_max
  use mod_particle_group_id, only: matching_part_config_indices
  use mod_particle_sim, only: group_num_from_id
  use mod_math_operators, only: gcd

  implicit none

  type(particle_sim),  intent(inout) :: sim
  type(edge_elements), intent(in)    :: edge_element_template !< a prepared set of edge elements
  
  type(wall_act_group), dimension(:), allocatable :: wall_act_groups
  type(type_wall_act_config) :: config
  type(type_part_create_scheme) :: create_scheme

  character(len=1000) :: identifier
  character(len=1000) :: msg !< error message
  integer :: i, j, k, Z, config_num_i, n_fluids, n_groups, idx_group, idx_act, i_target_group, each_nstep_part
  integer :: i_wall_acts, n_wall_acts !< total number of wall_action objects to make
  integer :: i_other, n_other !< number of wall_action objects in the "other" group
  integer :: i_part2self, n_part2self !< number of wall_action objects in the "part2self" group
  integer :: i_pump, n_pump, i_self_sputter, n_self_sputter, i_reflection, n_reflection !< number of wall_action objects of type "reflection"
  real*8  :: density_fraction, n_particles
  real*8  :: density_fraction_sum !< to check wether sum of density fractions is 1
  character(len=3), dimension(n_part_groups_max) :: sputter_target_ids="non" !< the id's of target groups that are being sputtered into
  integer, dimension(n_part_groups_max) :: i_sputter_group, n_sputter_group=0 !< how many wall actions in the sputter group of the id from sputter_target_ids
  logical, dimension(n_part_groups_max) :: group_scheme_set=.false. !< whether the group scheme was already set or not
  integer, dimension(n_part_groups_max,3) :: group_scheme_namelist_idx=-1 !< index of where the group create scheme was set in the namelist (part_group_configs (1) or fluid_configs (2), config_num, wall_act_num)
  integer :: n_sputter_groups=0 !< how many "fluid sputter to one" groups to make
  logical :: new_target !< whether the sputter target was already encountered before and written into sputter_target_ids (.false.), or not (.true.)
  integer :: idx_part2self,idx_other
  integer :: idx_offset !< how many non sputter groups there are

  if(sim%my_id == 0) write(*,*) "determining wall_actions from the namelist configs"

  ! --- determining number of wall_actions and groups necessary
  n_wall_acts    = 0
  n_other        = 0
  n_part2self    = 0
  n_pump         = 0
  n_self_sputter = 0
  n_reflection   = 0
  
  !from particles
  do i=1,n_part_groups
    do j=1,n_part_groups_max
      config = part_group_configs(i)%wall_act_configs(j)
      select case(trim(config%type))
      case("none") 
        cycle
      case("pump")
        n_part2self    = n_part2self    + 1
        n_pump         = n_pump         + 1
      case("self sputter")
        n_part2self    = n_part2self    + 1
        n_self_sputter = n_self_sputter + 1        
      case("reflection")
        n_part2self    = n_part2self    + 1
        n_reflection   = n_reflection   + 1
      case default
        n_other = n_other + 1
      end select

      n_wall_acts = n_wall_acts + 1
    end do
  end do
  
  !from the fluid
  do i=1,n_fluid_groups_max
    do j=1,n_part_groups_max
      config = fluid_configs(i)%wall_act_configs(j)
      select case(trim(config%type))
      case("none") 
        cycle
      case("fluid sputter")
        new_target=.true.
        ! check whether a group already exists for this target
        do k=1,n_part_groups_max
          if(sputter_target_ids(k) == config%target_group_id) then
            n_sputter_group(k) = n_sputter_group(k) + 1
            new_target=.false.
            exit
          end if
        end do

        !if no group exists yet, make it
        if(new_target) then
          !save this target as the first unused id in sputter_target_ids
          do k=1,n_part_groups_max
            if(sputter_target_ids(k) == "non") then
              sputter_target_ids(k) = config%target_group_id
              n_sputter_group(k) = 1
              exit
            end if
          end do

          n_sputter_groups = n_sputter_groups + 1
        end if

        !check and set the part_create_scheme for the group 
        if(config%supers_num_wall /= -1 .or. config%supers_weight_wall /= -1 .or. config%supers_ratio_wall /= -1) then
          do k=1,n_part_groups_max
            if(sputter_target_ids(k) == config%target_group_id) then
              if(group_scheme_set(k)) then
                select case(group_scheme_namelist_idx(k,1))
                case(1)
                  write(identifier,"(A)") "part_group_configs("
                case(2)
                  write(identifier,"(A)") "fluid_configs("
                case default
                  write(identifier,"(A)") "[ERROR in traceback]("
                end select
                write(msg,"(5A,I2,A,I2,A,I2,A,I2,4A)") "you can only specify one create scheme for all fluid_configs sputtering into a particle group, ", &
                                  "but for sputtering into particle group with id=",config%target_group_id," both a ",trim(identifier),group_scheme_namelist_idx(k,2), &
                                  ")%wall_act_configs(",group_scheme_namelist_idx(k,3),")%supers...wall option and a fluid_configs(",i,")%wall_act_configs(",j,")%super...wall ", &
                                  "option was set. This is not allowed, please provide only one create option (%supers...wall) for all fluids sputtering into ",config%target_group_id,"."
                call wrong_input(msg,sim%my_id,"")
              end if
              group_scheme_namelist_idx(k,:) = [2,i,j]
              group_scheme_set(k) = .true.
              exit
            end if
          end do
          ! if(trim(wall_act_groups(idx_group)%create_scheme%scheme == "non")
        end if
      case default
        n_other = n_other + 1
      end select
    
      n_wall_acts = n_wall_acts + 1
    end do
  end do

  if(n_wall_acts > n_part_groups_max**2 + n_part_groups_max*n_fluid_groups_max) then
    write(msg,"(A)") "size of wall_actions is bigger than should be possible?"
    call wrong_input(msg,sim%my_id,"")
  end if

  ! --- setting up the wall_act_groups
  ! determining the number of wall_act_groups to set up
  n_groups = 0
  idx_part2self = -1
  idx_other = -1
  if(n_part2self > 0) then
    n_groups = n_groups + 1
    idx_part2self = n_groups
  end if
  if(n_other > 0) then 
    n_groups = n_groups + 1
    idx_other = n_groups
  end if
  idx_offset = n_groups
  n_groups = n_groups + n_sputter_groups

  ! allocating and setting wall_act_groups settings
  allocate(wall_act_groups(n_groups))
  if(size(wall_act_groups,1) == 0) return ! avoid allocation issues below
  
  wall_act_groups(:)%contains_part2part = .false.
  if (n_part2self > 0) then
    allocate(wall_act_groups(idx_part2self)%wall_actions(n_part2self))
    wall_act_groups(idx_part2self)%group_type = "part2self"
    wall_act_groups(idx_part2self)%contains_part2part = .true.
  end if
  if(n_other > 0) then
    allocate(wall_act_groups(idx_other)%wall_actions(n_other))
    wall_act_groups(idx_other)%group_type = "other"
  end if
  
  do k=1,n_sputter_groups ! for the sputter groups
    i=k+idx_offset
    !generic group settings
    wall_act_groups(i)%group_type = "fluid sputter to one"
    wall_act_groups(i)%target_id = sputter_target_ids(k)

    !setting the creation scheme of the group
    i_target_group = group_num_from_id(sim,wall_act_groups(i)%target_id)
    if(i_target_group > 0) then !check that id was found before asking for the n_particles of this group
      n_particles = sim%groups(i_target_group)%n_particles
    else
      n_particles = 1 !just to avoid crashing here, in the generation of the wall_actions the wrong id will give a helpful error message, but we don't want to duplicate that check
    end if
    if(group_scheme_set(k)) then !if a scheme is set by one of the action in the input, use that for the group
      config_num_i = group_scheme_namelist_idx(k,2) ! config index
      idx_act = group_scheme_namelist_idx(k,3) ! action index
      select case(group_scheme_namelist_idx(k,1))
      case(1)
        config = part_group_configs(config_num_i)%wall_act_configs(idx_act)
        write(identifier,"(A,I2,A,I2,A)") "part_group_configs(",config_num_i,")%wall_act_configs(",idx_act,")"
      case(2)
        write(identifier,"(A,I2,A,I2,A)") "fluid_configs(",config_num_i,")%wall_act_configs(",idx_act,")"
        config = fluid_configs(config_num_i)%wall_act_configs(idx_act)
      case default
        write(msg,"(A,I2,A,I2,A,I2,A)") " bug in setting the creation scheme of wall_action group ",i,": group_scheme_namelist_idx(",k,",1)=",group_scheme_namelist_idx(k,1),", but it should be either 1 or 2..."
        call wrong_input(msg,sim%my_id,"")
      end select
      wall_act_groups(i)%create_scheme = part_create_scheme(config%supers_num_wall,config%supers_weight_wall,config%supers_ratio_wall,n_particles,my_id=sim%my_id,identifier=identifier)
    else ! no scheme was set in the input for this group, so use the default scheme for this group
      write(identifier,"(2A)") 'for wall_action group: "fluid sputter to one" -> ',sputter_target_ids(k)
      wall_act_groups(i)%create_scheme = part_create_scheme(-1,-1.d0,-1.d0,n_particles,default=supers_ratio_wall_default,my_id=sim%my_id,identifier=identifier)
    end if
    
    !allocation of wall_actions themselves
    allocate(wall_act_groups(i)%wall_actions(n_sputter_group(k)))
  end do

  ! --- filling out the wall actions in the groups from the configs
  i_wall_acts        = 0
  i_other            = 0
  i_part2self        = 0
  i_sputter_group(:) = 0
  i_pump             = 0
  i_self_sputter     = 0
  i_reflection       = 0
  
  !from particles
  do i=1,n_part_groups !loop over particle groups
    config_num_i = matching_part_config_indices(i)
    
    each_nstep_part = part_group_configs(config_num_i)%wall_act_each_nstep_part
    if(each_nstep_part /= -9999991) then
      if(each_nstep_part <= 0) then
        write(*,"(A,I2,A)") "part_group_configs(",config_num_i,")%wall_act_each_nstep_part <= 0 which is not allowed. Aborting."
        stop
      endif
      call sim%update_lcm_gcd(each_nstep_part)
      if (gcd_wall_acts == -9999991) then
        gcd_wall_acts = each_nstep_part
      else
        gcd_wall_acts = gcd(gcd_wall_acts,each_nstep_part)
      endif
    else
      if(sim%my_id == 0) write(*,"(A,I2,A)") "part_group_configs(",config_num_i,")%wall_act_each_nstep_part was not set, so particle-particle wall_actions originating from this species will be done once every fluid step"
    endif
        
    do j=1,n_part_groups_max !loop over wall_action configs
      config = part_group_configs(config_num_i)%wall_act_configs(j)
      
      select case(trim(config%type))
      case("none") 
        cycle
      case("pump")
        idx_group      = idx_part2self
        i_part2self    = i_part2self    + 1
        i_pump         = i_pump         + 1
        idx_act        = i_pump        
      case("self sputter")
        idx_group      = idx_part2self
        i_part2self    = i_part2self    + 1
        i_self_sputter = i_self_sputter + 1
        idx_act        = n_pump + i_self_sputter    
      case("reflection")
        idx_group      = idx_part2self
        i_part2self    = i_part2self    + 1
        i_reflection   = i_reflection   + 1
        idx_act        = n_pump + n_self_sputter + i_reflection
      case default
        idx_group = idx_other
        i_other = i_other + 1
        idx_act = i_other
      end select

      ! being here means it is a wall_action that should be used
      i_wall_acts = i_wall_acts + 1
      write(identifier,"(A,I2,A,I2,A,I3,A)") "for input namelist: particle_group_configs(",config_num_i,")%wall_act_configs(",j,"). (This corresponds to wall_action: ",i_wall_acts,")"
      call construct_wall_action(wall_act_groups(idx_group)%wall_actions(idx_act),sim,i,config,edge_element_template,.false.,each_nstep_part=each_nstep_part,input_identifier=identifier)      
    end do
  end do

  !from the fluid
  n_fluids = 0
  density_fraction_sum = 0
  do i=1,n_fluid_groups_max !loop over fluid groups
    Z = fluid_configs(i)%Z
    density_fraction=fluid_configs(i)%density_fraction
    if(Z /= -999) then
      if(density_fraction < -1.d98) then 
        write(msg,"(A,I2,A)") "it seems like you forgot to set fluid_configs(",i,")%density_fraction (it is still < 0, just like its preset)."
        call wrong_input(msg,sim%my_id,"")
      end if
      density_fraction_sum = density_fraction_sum + density_fraction
      n_fluids = n_fluids + 1
      !don't cycle the loop, because it should be checked whether the user set wall actions for a fluid config without having specified Z. This is checked in construct_wall_action
    end if
    do j=1,n_part_groups_max !loop over wall_action configs
      config = fluid_configs(i)%wall_act_configs(j)
      
      select case(trim(config%type))
      case("none") 
        cycle
      case("fluid sputter")
        do k=1,n_part_groups_max
          if(sputter_target_ids(k) == config%target_group_id) then
            idx_group = idx_offset + k
            i_sputter_group(k) = i_sputter_group(k) + 1
            idx_act = i_sputter_group(k)
            exit
          end if
        end do
        !overriding the create scheme to supers_num
        config%supers_num_wall    = 1
        config%supers_ratio_wall  = -1.d0
        config%supers_weight_wall = -1.d0
      case default
        idx_group = idx_other
        i_other = i_other + 1
        idx_act = i_other
      end select

      ! being here means it is a wall_action that should be used
      i_wall_acts = i_wall_acts + 1
      write(identifier,"(A,I2,A,I2,A,I3,A)") "for input namelist: fluid_configs(",i,")%wall_act_configs(",j,"). (This corresponds to wall_action: ",i_wall_acts,")"
      call construct_wall_action(wall_act_groups(idx_group)%wall_actions(idx_act),sim,i,config,edge_element_template,.true.,fluid_Z=Z,fluid_density_fraction=density_fraction,input_identifier=identifier)
    end do
  end do

  !check that the density_fractions add up to 1 if there is at least one fluid_config defined
  if(n_fluids > 0 .and. abs(density_fraction_sum - 1.d0) > 1.d-12) then
    write(msg,"(A,f24.14,A)") "the sum of fluid_configs(:)%density_fraction has to be 1 (of the fluid configs in use, i.e. with specified Z), but it is ",density_fraction_sum,". Please check your fluid_configs(:)%density_fraction."
    call wrong_input(msg,sim%my_id,"")
  end if

  !sanity checks on the pre calculated number of wall actions and the actual number. If this creates an error that must mean there's a bug inside this function somewhere
  if(i_wall_acts /= n_wall_acts) then
    write(msg,"(A,I3,A,I3,A)") "the total amount of wall_actions (",n_wall_acts,") is not equal to the number of initialised wall_actions (",i_wall_acts,"), so something went wrong."
    call wrong_input(msg,sim%my_id,"")
  end if

  if(n_wall_acts /= i_other + i_part2self + sum(i_sputter_group)) then
    write(msg,"(A,I3,A,I3,A)") "the total amount of wall_actions (",n_wall_acts,") is not equal to sum of initialised wall_actions from all groups (",i_other + i_part2self + sum(i_sputter_group),"), so something went wrong."
    call wrong_input(msg,sim%my_id,"")
  end if

  if(n_part2self /= i_pump + i_self_sputter + i_reflection) then
    write(msg,"(A,I3,A,I3,A)") "the total amount of part2self wall_actions (",n_part2self,") is not equal to sum of initialised part2self wall_actions (",i_pump + i_self_sputter + i_reflection,"), so something went wrong."
    call wrong_input(msg,sim%my_id,"")
  end if

  wall_act_groups(:)%constructed =.true.

end function wall_actions_from_config


!> Load eckstein sputtering yields and energy coefficients for this interaction
subroutine load_eckstein_data(this, sim)
  implicit none

  class(wall_action), intent(inout) :: this
  type(particle_sim), intent(in)    :: sim

  integer :: Z_origin, Z_target

  ! determining Z of origin and target (needed to read the correct data file)
  if (this%fluid2part) then
    Z_origin = this%fluid_Z
  else !< origin is particle group
    Z_origin = sim%groups(this%origin_group)%Z
  end if

  Z_target = sim%groups(this%target_group)%Z

  ! setting the yield object
  this%yield%Z_ion    = Z_origin
  this%yield%Z_target = Z_target
  this%yield%use_Yn_func = this%use_Yn_func

  ! reading the yield data
  call this%yield%read()

  if (.not. this%use_thompson) then ! use eckstein coefficients
    ! setting the energy object
    this%energy%Z_ion    = Z_origin
    this%energy%Z_target = Z_target
    this%energy%use_Yn_func = this%use_Yn_func

    !reading the energy data
    call this%energy%read()
  end if
end subroutine load_eckstein_data

!> does everything for the group (wall_actions & any other things like setting the create schemes)
subroutine do_wall_act_group(this, sim, post_evolution)
  implicit none
  
  class(wall_act_group), intent(inout)    :: this
  type(particle_sim),    intent(inout)    :: sim
  logical,               intent(in)       :: post_evolution !< whether this call is after the evolve particle groups call (.true.) or not (.false.)

  integer :: i, n_supers_tot, n_supers_child
  real*8  :: total_yield, yield_fraction

  !> if we're after the evolution, we should run part2part wall_actions
  if(post_evolution) then
    if(this%contains_part2part) then
      call all_acts_in_group(this,sim)
    end if
  else ! we should run the creation wall actions
    select case(trim(this%group_type))
    case("part2self") !< is all already done in post_evolution
      return
    case("other")
      call all_acts_in_group(this,sim)
    case("fluid sputter to one")
      if(sim%my_id == 0) write(*,"(3A)") '===== wall_action group: "fluid sputter to one" -> ',this%target_id,' ===== '

      ! calculating the partial yields and the total_yield
      total_yield = 0.d0
      do i=1,size(this%wall_actions,1)
        call this%wall_actions(i)%calc_fluid_yield(sim)
        total_yield = total_yield + this%wall_actions(i)%domain_integral*this%wall_actions(i)%weight_factor
      end do

      ! determining the number of supers for the group
      n_supers_tot = this%create_scheme%supers_to_create(sim%my_id,total_yield)
  
      if(sim%my_id == 0) then
        write(*,"(A50,' = ',es16.6)") "total sputter yield for this wall_act_group       ",total_yield
        if (trim(this%create_scheme%scheme) == "num") write(*,"(A50,' = ',I12)") "supers_num                                        ", this%create_scheme%supers_num
        if (trim(this%create_scheme%scheme) == "weight") write(*,"(A50,' = ',es16.6)") "supers_weight                                     ", this%create_scheme%supers_weight
        if (trim(this%create_scheme%scheme) == "ratio") write(*,"(A50,' = ',es16.6)") "supers_ratio                                      ", this%create_scheme%supers_ratio
        write(*,"(A50,' = ',I12)") "superparticles to create for this wall_act_group  ",n_supers_tot
      end if

      ! setting the number of supers for the children actions by equal distribution (at least 1)
      do i=1,size(this%wall_actions,1)
        yield_fraction = this%wall_actions(i)%domain_integral*this%wall_actions(i)%weight_factor / total_yield
        n_supers_child = max(1,nint(n_supers_tot * yield_fraction))
        this%wall_actions(i)%create_scheme%supers_num = n_supers_child
      end do

      !run children actions
      call all_acts_in_group(this,sim)

      if(sim%my_id == 0) write(*,"(A)") "===== end wall_action group ===== "
    case default
      call wrong_input("wrong ",sim%my_id,"do_wall_act_group")
    end select
  end if
  
end subroutine do_wall_act_group

!> runs all wall_actions in this group
subroutine all_acts_in_group(this, sim)
  implicit none
  class(wall_act_group), intent(inout)    :: this
  type(particle_sim),    intent(inout)    :: sim

  integer :: i

  do i=1,size(this%wall_actions,1)
    call this%wall_actions(i)%do(sim)
  end do
  
end subroutine all_acts_in_group

!> Perform the wall interaction according to the setting in the wall_action object
!> How and when to run depends on the details of the interaction
subroutine do_wall_action(this, sim, ev)
  use mod_atomic_elements, only: element_symbols
  use mod_parameters, only: n_plane, n_period
  use phys_module, only: use_manual_random_seed
  
  class(wall_action), intent(inout)    :: this
  type(particle_sim), intent(inout)    :: sim
  type(event), intent(inout), optional :: ev   !< this is here so that it is compatible with the event structure, but we can remove it if we want to move away from it

  integer :: i

  ! --- setup  

  !check whether this action should be run right now
  if(.not. this%fluid2part) then ! always run fluid2part actions when called
    !> skip any part2part action if it is not it's time to be run
    if(.not. (mod(sim%istep_inner_loop,this%each_nstep_part)==0 .or. sim%istep_inner_loop==sim%nstep_inner_loop)) return
  endif

  if(sim%my_id == 0) write(*,"(A)") "--- wall_action: "//trim(this%name)//" --- "

  ! check whether the constructor was used (so that all other sanity checks can be done once in the constructor)
  if (.not. this%constructed) then
    write(*,*)'=======================ERROR!!=================================='
    write(*,*)"particle wall_action object of type '"//trim(this%type)//"'"
    write(*,"(A,I2,A,I2, A)") "with origin ",this%origin_group," and target ", this%target_group, " has not finished it's construction"
    write(*,*)'please use the constructor'
    stop
  end if

  ! updating the timestep in SI (as tstep can change)
  call this%update_delta_t()

  ! --- underlying function calls
  if(this%fluid2part) call fluid2part_action(this, sim)

  if(this%part2self) call part2self_action(this, sim)
  
  ! in the future add part2other

  ! --- area for writing the projected diagnostic
  if (sim%istep_inner_loop >= sim%nstep_inner_loop .or. sim%istep_inner_loop == -1) then !if at the last step of inner loop, or outside inner loop
    this%i_step_diag = this%i_step_diag + 1
    if (this%i_step_diag .ge. this%n_step_diag) then
      call write_wall_project_vtk(this, sim)
    end if
  end if
end subroutine do_wall_action


!> Routine to do the fluid to particle wall interactions.
!> Models fluid sputtering ("fluid sputter") and wall 
!> recombination ("wall recomb")
!>
!> The fluid-particle interaction is done by calculating 
!> the yield by integrating over the velocity distribution.
!> Then particles are sampled using these local yields to represent
!> the incoming flux (with the weight already adjusted for particles 
!> resulting from the interaction). 
!> The particles then undergo a particle-particle interaction to 
!> determine their new energy, and the new particles are stored in
!> free slots in the group's particle array
subroutine fluid2part_action(this, sim)
  use mpi_mod
  use mod_atomic_elements, only: element_symbols
  use mod_interp, only: interp_RZ
  use mod_particle_create, only: free_particle_indices
  use mod_particle_types, only: initialize_particle_to_zero
  use mod_edge_elements, only: sample_edge_elements
  
  class(wall_action), intent(inout) :: this
  type(particle_sim), intent(inout) :: sim
  
  type(particle_kinetic_leapfrog) :: particle
  integer :: j, i_p, n_supers, n_supers_loc, i_rng
  integer :: q, Z
  real*8 :: E !< [eV] particle energy  (eV because of eckstein coeffs).
  real*8 :: n_e, T_e, T_i, Te_eV, Ti_eV
  real*8,  allocatable :: xyz_sampled(:,:), st_sampled(:,:), rng_sample(:,:) !< (3,n_supers_loc), (2,n_supers_loc), (6,n_supers_loc)
  integer, allocatable :: i_elm_sampled(:) !< (n_supers_loc)
  logical :: do_main !> whether to do the main calculation (we can't just return early because that breaks the MPI_REDUCE at the end of the subroutine)

  !> For check free particles
  integer, allocatable, dimension(:) :: i_free

  ! diagnostics
  real*8, dimension(n_global_diagnostics) :: diagnostics         !< diagnostics for the global wall loads
  real*8, dimension(n_global_diagnostics) :: diagnostics_all_mpi !< MPI reduced version of diagnostics
  real*8  :: mol_binding_E=3.526d-19, ion_binding_E=2.18d-18  !< ! (J) default values are only true for hydrogen, should be the sum of ionisation energies from 0 to q for impurities.

  do_main=.true.

  diagnostics = 0.d0

  ! calculate the yield if not done so already
  if(.not. this%yield_calculated) then
    call this%calc_fluid_yield(sim)
  end if

  ! check if we need to do anything
  if (this%domain_integral .le. 1d-12) then
    if(sim%my_id == 0) write(*,"(3A,I2,3A)") "fluid2part wall_action ",trim(this%type)," with origin ",this%origin_group," and target_id=",sim%groups(this%target_group)%id," has 0 yield. returning"
    do_main = .false. ! Move along, nothing to do
  else ! normal behaviour
    ! determine how many particles to initialise on this MPI proces
    n_supers = this%create_scheme%supers_to_create(sim%my_id,this%domain_integral*this%weight_factor)
    n_supers_loc = calc_n_particles_per_mpi(n_supers, sim%n_mpi, sim%my_id)

    ! determine indices of free particles
    call free_particle_indices(sim%groups(this%target_group)%particles, i_free, n_needed=n_supers_loc)  
    n_supers_loc = size(i_free,1)

    ! check if we need to do anything on this mpi process
    if(n_supers_loc == 0) then
      write(*,"(A,I2,3A,I2,3A)") "on sim%my_id=",sim%my_id," fluid2part wall_action ",trim(this%type)," with origin ",this%origin_group," and target_id=",sim%groups(this%target_group)%id," will not create any particles. returning"
      do_main=.false. ! Move along, nothing to do
    end if
  end if
  
  if(do_main) then  
    allocate(rng_sample(6,size(i_free)))
    allocate(xyz_sampled(3,size(i_free)))
    allocate(st_sampled(2,size(i_free)))
    allocate(i_elm_sampled(size(i_free)))

    ! We need to properly use all RNGS here to avoid missing numbers
    ! needs default(shared) for gfortran
#ifdef __GFORTRAN__
    !$omp parallel default(shared) &
#else
    !$omp parallel default(none) &
    !$omp shared(this, rng_sample, n_supers_loc) &
#endif 
    !$omp private(i_rng, j)
    i_rng = 1
    !$ i_rng = omp_get_thread_num()+1
    !$omp do schedule(static,1)
    do j=1,n_supers_loc
      call this%rng(i_rng)%next(rng_sample(:,j))
    end do
    !$omp end do
    !$omp end parallel

    q = this%fluid_Z
    if (q .le. 0) q = 1 ! deuterium, tritium special case
    q = min(q, 4) ! limit to 4 for divertor conditions
    Z = this%fluid_Z

    call sample_edge_elements(this%fluid_yield_integral, this%res, 1, n_supers_loc, rng_sample(1:3,:), xyz_sampled, st_sampled, i_elm_sampled)

    if (sim%my_id .eq. 0) then
      write(*,"(A,i8,A,A,A,A,A,i2,3A,i2,A,es16.6,A,es16.6)") "fluid2wall will create ", n_supers," ", element_symbols(sim%groups(this%target_group)%Z),&
        " from ", element_symbols(Z), " in group ", this%target_group, &
      " (ID=",sim%groups(this%target_group)%id,", Z=", sim%groups(this%target_group)%Z, ") with total weight ", this%domain_integral*this%weight_factor, "  particles flux #/s : ", this%domain_integral*this%weight_factor/this%delta_t
    end if

    select type (pa => sim%groups(this%target_group)%particles)
    type is (particle_kinetic_leapfrog)
#ifdef __GFORTRAN__
    !$omp parallel default(shared) &
#else
    !$omp parallel default(none) &
    !$omp shared(this, sim, rng_sample, xyz_sampled, st_sampled, i_elm_sampled, i_free, &
    !$omp q, Z, n_supers_loc, n_supers) &
#endif
    !$omp private(i_rng, j, E, T_e, T_i, Te_eV, Ti_eV, n_e, particle, i_p) &
    !$omp reduction(+:diagnostics)
    i_rng = 1
    !$ i_rng = omp_get_thread_num()+1
    !$omp do schedule(static,1)
    do j=1,n_supers_loc
      !> make a new particle which at the end of the do loop will be written into a free particle in the array
      call initialize_particle_to_zero(particle)

      particle%q = int(q,1)
      particle%i_elm = i_elm_sampled(j)
      if (i_elm_sampled(j) .le. 0) cycle
      particle%st = st_sampled(:,j)
      call interp_RZ(sim%fields%node_list, sim%fields%element_list, i_elm_sampled(j), &
        st_sampled(1,j), st_sampled(2,j), &
        particle%x(1), &
        particle%x(2))
      particle%x(3) = xyz_sampled(3,j) ! phi coordinate from sampling

      !> weight of fluid particle is equally distributed as a fraction of the incoming flux. Such that the sum of all incoming fluid particles,
      !> is the total amount of incoming particles over the edge domain area * delta_t
      !> multiplication with this%weight_fraction will be done in single_self_interaction
      !> note that in the off case there are not enough free particles this is wrong because the sum of n_supers_loc < n_supers then (but stuff starts breaking then anyway)
      !> this could be fixed using an MPI reduce, but that costs communication time for something which is normally never necessary
      particle%weight = this%domain_integral/n_supers

      ! Calculate temperature at this position to determine particle energy
#ifdef WITH_TiTe
      call sim%fields%calc_NeTeTi(sim%time, i_elm_sampled(j), st_sampled(:,j), xyz_sampled(3,j), n_e=n_e, T_e=T_e, T_i=T_i)
      Ti_eV = T_i * K_BOLTZ / EL_CHG
#else
      call sim%fields%calc_NeTe(sim%time, i_elm_sampled(j), st_sampled(:,j), xyz_sampled(3,j), n_e=n_e, T_e=T_e)
#endif
      Te_eV = T_e * K_BOLTZ / EL_CHG

      select case(trim(this%type))
      case("wall recomb")
        ! determine E
#ifdef WITH_TiTe
        call sample_fluid_particle_energy(Te_eV, rng_sample(4:6,j), Z, E, Ti_eV=Ti_eV)
#else
        call sample_fluid_particle_energy(Te_eV, rng_sample(4:6,j), Z, E)
#endif

        ! determine outcoming particle
        call single_self_interaction(this, sim, particle, this%rng(i_rng), diagnostics, E, "reflection")
      case("fluid sputter")
        ! The yield at a specific position is given by
        ! \[
        !   \int_v Y(E) f(v) dv
        ! \]
        ! where $f(v)$ is a maxwellian and $E$ includes the sheath potential and the Bohm outflow
        ! condition additionally.
        !
        ! To now calculate the energy of the sputtered particle we multiply the sputtered energy
        ! coefficient with E of a particle sampled from f(v).
        ! Taking the sputtered energy coefficient * Y as a weight factor and sampling from f(v) will do the trick.
        ! we need to normalize with the sputter yield at that position, which we have calculated before.
        ! Basically this is a weighted average of Y_E(E) * E, weighted with Y(E) f(E)
        
        ! If sampling from the incoming energy distribution function, the
        ! sputtered energy coefficient needs to be reweighed with the sputtering
        ! yield at this energy (since the tail contributes more)
        ! This is commented below since we have simplified the model for now to
        ! work at a fixed energy of 3 q T_e + 2 T_i, so we don't need to do this
        ! anymore. The extension to realistic IEDFs should be done later, so 
        ! I've kept some of the code around.
#ifdef WITH_TiTe
        E = 2 * Ti_eV
#else
        E = 2 * Te_eV !< from the bohm criterion, E = E_sheath_entrance + E_sheath_acceleration = 2 T_i + 3 q T_e, but for now T_i = T_e
        ! so E_sheath_entrance  = 2 T_i = 2 T, and E_sheath_acceleration will be added later
#endif
        call single_self_interaction(this, sim, particle, this%rng(i_rng), diagnostics, E, "self sputter", .true.)

        ! Non-implemented alternative to the above method:
        ! We could sample directly from Y(E) Y_E(E) f(E), but I don't know how to do this generally.
        ! That would have the advantage of better distribution of statistics (more uniform weights).

        !sputtering_yield = this%yield%interp(E, theta)
        ! Workaround if sputtered energy coeff threshold is lower than sputtering
        ! threshold: use sputtered energy coeff just above threshold instead
        ! (note: all this doesn't take into account theta properly)

        !av_yield = fluid_sputtering_yield(this%yield, T_eV, Z, theta)
        ! we could probably avoid the calculation of fluid_sputtering_yield by
        ! using the discretisation we just sampled from (if theta is constant)
        !if (av_yield .le. 1d-18) av_yield = 1d-6 ! does not matter since then sputtering_yield must be 0, just to avoid a NaN below
        
        ! now we weigh the particles with the prevalence of this energy in sputtered particles, i.e. sputtering_yield
        ! over the integral of sputtering_yield, which we calculate (again)
        !particle%weight = &
        !particle%weight * &
          !sputtering_yield / av_yield
      case default
        call wrong_interaction_type(this%type)
      end select

      ! write the created particle to a free slot in the array
      i_p = i_free(j)
      pa(i_p) = particle ! assignment(=+ operator is defined for particle_base as copy, so this works as you would intuitively think
    end do
    !$omp end do
    !$omp end parallel
    class default
      write(*,*) 'Target particle type not implemented for fluid2part actions (origin/target id)', this%origin_group, sim%groups(this%target_group)%id
      stop
    end select
  end if !do_main

  call reduce_global_diag(diagnostics)

  ! write global diagnostics, taylored for fluid to particle
  if (sim%my_id .eq. 0) then
    write(*,'(A,1f14.0)' ) "superparticles created        = ", diagnostics(i_wall_part_out) 
    write(*,'(A,2es16.6)') "particle flux (in/out) [#/s]  = ", this%domain_integral*this%weight_factor/this%delta_t,diagnostics(i_wall_flux_out)/this%delta_t 
    if(trim(this%type) == "wall recomb") then
      write(*,'(A,2es16.6)') "heatflux (in/out) [W]            = ", diagnostics(i_wall_heat_in)/this%delta_t,diagnostics(i_wall_heat_out)/this%delta_t 
      write(*,'(A50,1es16.8)') "atom wall-assisted recombination power [W] = ",      diagnostics(i_wall_flux_in)                                  * ion_binding_E / this%delta_t
      write(*,'(A50,1es16.8)') "molecule wall-assisted recombination power [W] = ", (diagnostics(i_wall_flux_in) - diagnostics(i_wall_flux_refl)) * mol_binding_E / this%delta_t
      write(*,'(A50,1es16.8)') "Power to (fast) reflected atoms [W] = ",             diagnostics(i_wall_heat_refl)                                                / this%delta_t
      !< atom wall-assisted recombination power = recycled flux * 13.6 eV. All ions are neutralized on the wall. This increaes the heat load on the wall ~stangeby2000 p.653
      !< molecule wall-assisted recombination power = thermal desorption flux*2.2 eV. When neutrals on the wall form neutrals, the wall heat load is increased by 2.2 eV per molecule. ~stangeby2000 p.653
      !< power to (fast) reflected atoms = energy retained by reflected neutrals. This energy is not deposited on the wall, thus decreases the plasma heat load.
      !  From ITER PFPO-1 test in 2D : energy_reflected_all > enery_wall_recombi_all >> energy_mol_recombi_all
    else ! fluid sputter
      write(*,'(A,2es16.6)') "heatflux (out) [W]            = ", diagnostics(i_wall_heat_out)/this%delta_t 
    end if
  endif

  this%yield_calculated = .false.

  if(allocated(rng_sample)) deallocate(rng_sample, xyz_sampled, st_sampled, i_elm_sampled)
  if(allocated(i_free)) deallocate(i_free)
end subroutine fluid2part_action


!> calls single_self_interaction() for all particles in the specified this%target_group
!> also prints the global diagnostics to the output file
subroutine part2self_action(this, sim)
  use phys_module, only: use_manual_random_seed
  use mod_polygon, only: inside_polygon
  use mod_interp,  only: interp_RZ
  
  class(wall_action), intent(inout) :: this
  type(particle_sim), intent(inout) :: sim
  
  real*8, dimension(n_global_diagnostics) :: diagnostics         !< diagnostics for the global wall loads
  real*8, dimension(n_global_diagnostics) :: diagnostics_all_mpi !< MPI reduced version of diagnostics

  integer :: j, i_rng
  
  diagnostics = 0.d0

  select type (pa => sim%groups(this%target_group)%particles)
  type is (particle_kinetic_leapfrog)
  if(use_manual_random_seed) then
    !$ call omp_set_schedule(omp_sched_static,10)
  else
    !$ call omp_set_schedule(omp_sched_dynamic,10)
  end if
  
#ifdef __GFORTRAN__
  !$omp parallel default(shared) & ! workaround for Error: ‘__vtab_mod_pcg32_rng_Pcg32_rng’ not specified in enclosing ‘parallel’
#else
  !$omp parallel default(none) &
  !$omp shared(this, sim)      & 
#endif
  !$omp private(i_rng)         &
  !$omp reduction(+:diagnostics)

  i_rng = 1
  !$ i_rng = omp_get_thread_num()+1
  !$omp do schedule(runtime)
  do j = 1,size(sim%groups(this%target_group)%particles,1)
    ! Skip if this particle is not lost in a specific location (i_elm .eq. 0 means lost 'somewhere')
    if (pa(j)%i_elm .ge. 0) cycle

    if(this%only_in_polygon) then
      if (.not. inside_polygon(size(this%poly_R),this%poly_R,this%poly_Z,pa(j)%x(1),pa(j)%x(2))) then !if particle not inside polygon, don't do the action
        cycle
      end if
    end if

    !> Place particle back into domain
    pa(j)%i_elm = -pa(j)%i_elm !reset i_elm to positive value (i_elm, s, t from where particle was lost at the boundary are known from find_RZ_nearby in mod_particle_evolution)
    call interp_RZ(sim%fields%node_list,sim%fields%element_list,pa(j)%i_elm,pa(j)%st(1),pa(j)%st(2),pa(j)%x(1),pa(j)%x(2)) !get corresponding R,Z at boundary
    
    !> do single particle wall interaction
    call single_self_interaction(this, sim, pa(j), this%rng(i_rng), diagnostics)
  end do
  !$omp end do
  !$omp end parallel
  class default
    write(*,*) "part2self_action not implemented for this kinetic type, group id=",sim%groups(this%target_group)%id
    call exit(13)
  end select
  
  call write_global_diag(this, sim, diagnostics)
end subroutine part2self_action


!> The interaction of a single particle with the wall, only affecting that super particle (self sputter or reflect)
subroutine single_self_interaction(this, sim, particle, rng, diagnostics, E_in, type_in, weight_preadjusted)
  use phys_module, only: part_kill_ratio

  implicit none

  class(wall_action),                      intent(inout) :: this
  type(particle_sim),                      intent(in)    :: sim
  type(particle_kinetic_leapfrog),         intent(inout) :: particle           !< particle to undergo interaction
  class(type_rng),                         intent(inout) :: rng                !< RNG object of the current openmp thread
  real*8, dimension(n_global_diagnostics), intent(inout) :: diagnostics        !< diagnostics for the global wall loads
  real*8,            optional,             intent(in)    :: E_in               !< [eV] energy of the incoming particle (if not specified will be determined from particle%v)
  character(len=*),  optional,             intent(in)    :: type_in            !< type of single interaction (either "self sputter" or "reflection") (if not specified will be set to this%type) 
  logical,           optional,             intent(in)    :: weight_preadjusted !< whether the weight was already adjusted beforehand to take the yield into account (true) or not (false, default)

  real*8 :: n_e, T_e, T_i, theta
  real*8 :: E !<[eV] particle energy. E is in [eV] in this subroutine, because of eckstein coeffs.
  real*8 :: vector_normal(3)
  logical :: fast_reflection !< whether the reflection is a fast reflection or a thermal desorption (not that release is instant, but the energy of the reflected particle is different)
  real*8 :: yield, energy_coeff, Te_eV, Ti_eV, fast_reflect_chance, v_new
  real*8 :: u(2), p_kill(1)
  character(len=20) :: local_type !< which single particle interaction to do, used to call self interaction from within fluid2part_action (=type_in if present, else =this%type)
  logical :: skip_yield !< if weight_preadjusted = true, then the yield calculation should be skipped
  
  ! determine the type, this can be different from this%type if single_self_interaction is called from within fluid2part
  if(present(type_in)) then
    local_type = type_in
  else
    local_type = this%type
  end if

  if (trim(local_type) /= "pump") then ! for the pump we want to see the difference between before and after pumping
    ! pre-update weight of simulated particle according to weight factor for correct incoming diagnostics
    particle%weight = this%weight_factor * particle%weight 
  end if

  ! set the incoming particle energy
  if(present(E_in)) then
    E = E_in
  else
    ! calculate the energy associated with the velocity of the particle (in eV)
    E = 0.5d0*sim%groups(this%target_group)%mass*ATOMIC_MASS_UNIT*dot_product(particle%v, particle%v)/EL_CHG !< must be in eV
  end if

  ! determine whether to calculate the yield
  if(present(weight_preadjusted)) then
    skip_yield = weight_preadjusted
  else
    skip_yield = .false.
  end if

  ! use normal vector and velocity of particle to determine incoming angle
  ! cos(theta) = (n . v)/ (||n||.||v||)
  vector_normal = wall_normal_vector(sim%fields%node_list, sim%fields%element_list, particle%i_elm, particle%st(1), particle%st(2))

  ! Hard-code theta to 0 to fix issues with sputtering module at strange angles
  ! the angle calculation should be revisited. Before using theta != 0 the
  ! surface roughness should be estimated, as this gives a distribution of
  ! impact angles as well
  theta = 0.d0 
  !> old theta 
  !theta = acos(dot_product(-vector_normal,particle%v)/norm2(particle%v))*180.d0/PI !< acos gives results in radians
  ! ! theta must be in degrees as the theta_star is also in degrees
  ! if (abs(theta) .gt. 91) then
  !   ! This is like an assert, it cannot really happen... but it does
  !   !!$omp critical
  !   !write(*,*) 'incoming angle warning', theta, vector_normal, particle%v
  !   !!$omp end critical
  ! end if

  if (trim(local_type) /= "pump") then
    ! Update the particle energy from the potential drop in the sheath
#ifdef WITH_TiTe
    call sim%fields%calc_NeTeTi(sim%time, particle%i_elm, particle%st, particle%x(3), n_e=n_e, T_i=T_i, T_e=T_e)
    Ti_eV = T_i * K_BOLTZ / EL_CHG
#else
    call sim%fields%calc_NeTeTi(sim%time, particle%i_elm, particle%st, particle%x(3), n_e=n_e, T_e=T_e)
#endif
    Te_eV = T_e * K_BOLTZ / EL_CHG
  
    E = E + simple_potential_drop(int(particle%q,4),Te_eV)
  end if

  ! store this particle's contribution to incoming particle, heatflux and flux onto the wall
  diagnostics(i_wall_part_in) = diagnostics(i_wall_part_in) + 1
  diagnostics(i_wall_flux_in) = diagnostics(i_wall_flux_in) + particle%weight
  diagnostics(i_wall_heat_in) = diagnostics(i_wall_heat_in) + particle%weight * E *EL_CHG

  !> determining interaction yield and new energy depending on interaction type
  select case (trim(local_type))
  case ("pump")
    yield = this%weight_factor
    
    !> storing this particle's contribution on a 2D edge element patch grid as diagnostic
    call particle_projection_diagnostic(this, sim, particle, E, (1-yield))
  case ("reflection")
    !> a particle can either bounce of the wall (fast_reflection=.true.) or be thermally released
    !> whether a particle reflects directly is determined through eckstein coefficients set for this goal
    fast_reflect_chance = this%yield%interp(E,theta)
    
    call rng%next(u)
    if (u(1) .le. fast_reflect_chance) then
      fast_reflection = .true.
    else
      fast_reflection = .false. 
    end if
    
    !> assume wall saturation (pumping implementation is done separately)
    yield = 1.d0

    !> storing this particle's contribution on a 2D edge element patch grid as diagnostic
    call particle_projection_diagnostic(this, sim, particle, E, yield)

    !> determine new energy
    if (fast_reflection) then
      ! still some energy and momentum can be lost at the reflection against the wall, this is modelled using another set of eckstein coefficients
      energy_coeff = this%energy%interp(E,theta)
      E = energy_coeff * E

      ! since we have wall_flux_in, and wall_flux_in = wall_flux_refl + wall_flux_therm, we also know wall_flux_thermal. Similarly we know wall_heat_thermal
      diagnostics(i_wall_flux_refl)   = diagnostics(i_wall_flux_refl) + particle%weight
      diagnostics(i_wall_heat_refl)   = diagnostics(i_wall_heat_refl) + particle%weight * E * EL_CHG
    else ! thermal release
      E = (800.d0 + 273.d0) *K_BOLTZ/EL_CHG! must be in eV (800 degrees celsius)
    endif  
  case ("self sputter")
    ! use eckstein sputtering coefficients to determine both the sputter yield and resulting energy
    
    if(skip_yield) then
      yield = 1.d0
    else
      yield = this%yield%interp(E,theta)
    end if

    !> exponential self sputtering for yield > 1
    if (yield .gt. 1.d0 + 1.d-12) then
      !$omp critical
      write(*,"(A,f5.0,A,f8.3)") "> 1 self-sputtering detected, E=", E, "yield=", yield
      !$omp end critical
    end if

    !> storing this particle's contribution on a 2D edge element patch grid as diagnostic
    call particle_projection_diagnostic(this, sim, particle, E, yield)

    !> determining the energy of the particle post sputtering
    if (this%use_thompson) then
      call rng%next(u)
      ! Option below to remove the highest 2% of the distribution by clipping u (hacky)
      ! u = min(u, 0.98d0)
      E = sample_dist(this%E_dist, u(1))
    else
      !> avoiding numerical issues with E being too small to calculate energy_coeff
      if (E < this%energy%E_threshold + 1d0) then
        !$omp critical
        write(*,*) "WARNING: E too small for yields",E,this%energy%E_threshold,"setting E to just above threshold, please expand coefficients range"
        !$omp end critical
        E = this%energy%E_threshold + 1d0
      end if

      energy_coeff = this%energy%interp(E,theta)
      E = energy_coeff * E
    end if

  case default
    write(*,*) "ERROR: unknown single_self_interaction type",local_type
  end select
  
  ! update weight of simulated particle after the wall interaction
  if (particle%weight .le. sim%groups(this%target_group)%average_weight * part_kill_ratio .and. yield .le. 1.d0) then
    call rng%next(p_kill)
    if (p_kill(1) .le. (1-yield)) then
      particle%i_elm  = 0 !takes the particle out of active use

      diagnostics(i_super_killed) = diagnostics(i_super_killed) + 1
      diagnostics(i_removed)      = diagnostics(i_removed)      + particle%weight
      
      return ! we explicitly don't want this particle to be updated or counted in the outgoing diagnostics anymore
    
    endif ! else do nothing to the weight of this particle
  else 
    diagnostics(i_removed) = diagnostics(i_removed) + particle%weight*(1-yield)

    particle%weight = yield * particle%weight
  endif

  ! use E from previous section to calculate velocity in one 
  v_new = sqrt(2.d0* E *EL_CHG/(sim%groups(this%target_group)%mass * ATOMIC_MASS_UNIT))
  
  ! give particle a new direction:
  ! Calculate vector normal and select a random vector with a cosine distribution in angle between the normal and itself
  call rng%next(u)
  particle%v =  v_new * sample_cosine(u(1:2),vector_normal) 
  ! [[not sure what this comment is about]] Since it is a neutral the half-step for boris method does not matter at all

  ! wall interactions typically neutralise the particles if they used to have charge
  ! this is assumed to be true for pump surfaces as well, as particles coming from the pump duct into the vessel will likely be neutral
  particle%q = 0_1

  ! after the wall interaction, the particle is now considered a new particle, so update i_life and t_birth
  particle%i_life = particle%i_life + 1
  particle%t_birth = sim%time
  ! For particle-particle sputtering we might want them to have the same identifiers
  ! if so comment the line above

  !> nan check (in fortran, for x=nan, x == x will return false)
  if (any(particle%x .ne. particle%x) .or. E .ne. E .or. particle%weight .ne. particle%weight) then
    !$omp critical
    write(*,*) 'ERROR: removing particle with nans in function single_self_interaction() (x,E,w,i_elm):', particle%x, E, particle%weight, particle%i_elm
    !$omp end critical
    particle%i_elm = 0 ! skip this one since sputtering went wrong
  end if

  ! store this particle's contribution to outgoing particle, heatflux and flux onto the wall
  diagnostics(i_wall_part_out) = diagnostics(i_wall_part_out) + 1
  diagnostics(i_wall_flux_out) = diagnostics(i_wall_flux_out) + particle%weight
  diagnostics(i_wall_heat_out) = diagnostics(i_wall_heat_out) + particle%weight * E *EL_CHG
  
end subroutine single_self_interaction


!> calculates the fluid yield and stores the results in the wall action
subroutine calc_fluid_yield(this,sim)
  use mod_edge_elements, only: integrate_edge_elements

  implicit none
  
  class(wall_action),                      intent(inout) :: this
  type(particle_sim),                      intent(in)    :: sim

  ! updating the timestep in SI (as tstep can change)
  call this%update_delta_t()

  !> this subroutine will calculate the incident ion flux over every fluid species on edge domain
  !> And the resulting yield of created particles (in atoms/m^2 during delta_t)
  call project_sputter_vars_on_edge(this, sim)

  ! determine integral over the domain
  call integrate_edge_elements(this%fluid_yield_integral, 1, this%domain_integral, this%res)

  this%yield_calculated = .true.
end subroutine calc_fluid_yield


!> updates this%delta_t, as that can change when tstep changes during the simulation
subroutine update_delta_t(this)
  use phys_module, only: tstep, central_mass, central_density
  use constants, only: MU_ZERO, ATOMIC_MASS_UNIT

  implicit none
  
  class(wall_action), intent(inout)    :: this
  
  this%delta_t = (tstep*sqrt((MU_ZERO * central_mass * ATOMIC_MASS_UNIT * central_density * 1.d20)))

end subroutine update_delta_t


!> The potential drop from a debye sheath. Could support two-temperature model later
pure function debye_potential_drop(q, T_eV) result(U_drop)
  use constants, only: TWOPI, ATOMIC_MASS_UNIT, MASS_ELECTRON
  use phys_module, only: central_mass
  integer, intent(in) :: q
  real*8, intent(in) :: T_eV !< Local temperature in eV
  real*8 :: T_i, T_e, U_drop
  ! Equal temperatures
  T_i = T_eV
  T_e = T_eV

  !> Potential drop in eV
  U_drop = 0.5d0 * log((TWOPI * MASS_ELECTRON/(central_mass * ATOMIC_MASS_UNIT))*(1.d0+T_i/T_e))
end function debye_potential_drop


!> Calculate the energy gain of a potential drop from a sheath in the simplest model possible
pure function simple_potential_drop(q, T_eV) result(ion_energy)
  integer, intent(in) :: q
  real*8, intent(in) :: T_eV !< Local temperature in eV
  real*8 :: ion_energy !< The energy of the outgoing ion in eV
  
  ion_energy = 3.d0*real(q,8)*T_eV !< for sputtering from fluid perspective. Add the original energy E to this
end function simple_potential_drop


!> Integrate the sputtering yield over the distribution of incoming velocities.
pure function fluid_sputtering_yield(coeff, T_eV, Z, theta) result(yield)
  use gauss
  class(eckstein_coeff_set), intent(in) :: coeff
  real*8,                    intent(in) :: T_eV  !< Plasma temperature in eV
  integer,                   intent(in) :: Z     !< Atomic number of the incoming particles
  real*8,                    intent(in) :: theta !< angle of impact (usually assumed 0) in degrees
  real*8                                :: yield !< The sputter yield in atoms/ion

  real*8 :: U_drop
  integer :: q

  if (Z .le. 0) then
    q = 1
  else
    q = min(Z, 4) ! cap to 4 for divertor conditions
  end if 

  ! We use a simplified model for now! 2 T_i + 3 q T_e
  U_drop = simple_potential_drop(q, T_eV) ! assume particle has full charge
  yield = coeff%interp(2*T_eV + U_drop, theta)

  ! Alternative but unused version not assuming the simplified model:
  
  ! Since we use inverse transform sampling on u to calculate the energy we can
  ! just integrate over u from 0 to 1 to cover the whole distribution.
  ! Do this with n subelements, using gaussian quadrature in each element
  ! the subintervals go from 1/2 to n-1/2 to avoid using 0 and 1, since 1 should
  ! lead to infinity for sampling from a gaussian. Skipping the first part is reasonable
  ! since the sputtering yield will be very low there. For the high energies a maxwellian
  ! is perhaps not even the best approximation so that is probably not so bad either.
  
  ! The returned yield is averaged over the maxwellian at T_eV + the potential drop
  
  ! real*8 :: E
  ! integer :: i, j, k
  ! integer, parameter :: n_interval = 4 !< number of intervals to calculate. (using 1 is already pretty good)
  ! real*8, parameter :: idu = 1.d0/real(n_interval,8) !< interval size
  ! real*8 :: u(3) !< the integration point

  ! yield = 0.d0
  ! if (T_eV .le. 1d-1) return
  ! do i=0,n_interval*n_gauss-1
  !   u(1) = (real(i/n_gauss,8) + xgauss(mod(i,n_gauss)+1))*idu
  !   do j=0,n_interval*n_gauss-1
  !     u(2) = (real(j/n_gauss,8) + xgauss(mod(j,n_gauss)+1))*idu
  !     do k=0,1 ! we only use the sign of this one
  !       u(3) = real(k,8)

  !       ! add to this energy the plasma sheath potential
  !       call sample_fluid_particle_energy(T_eV, u, Z, E)
  !       U_drop = simple_potential_drop(q, T_eV) ! assume particle has full charge
  !       yield = yield + coeff%interp(E + U_drop, theta)
  !     end do
  !   end do
  ! end do
  ! yield = yield / (2*n_interval**2)
end function fluid_sputtering_yield


!> Calculate the flux to and some diagnostics for fluid flux in
!> a period delta_t.
!> 
!> Fluid_yield_integral contains a single scalar, the incoming 
!> fluid flux. Assume all particles are moving at the same
!> velocity, so multiplying with the relative density is enough 
!> to get the flux of a specific species.
!>
!> Assume that the impact angle of all particles is 0
subroutine project_sputter_vars_on_edge(this, sim)
  use mod_atomic_elements, only: atomic_weights
  use phys_module, only: central_mass, xpoint, xcase, min_sheath_angle, gamma, n_vertex_max, bcs, D_par, central_density
  use data_structure, only: type_node
  use diffusivities, only: get_dperp
  use equil_info, only: get_psi_n
  use constants, only: MU_ZERO, ATOMIC_MASS_UNIT
  
  type(wall_action),  intent(inout) :: this
  type(particle_sim), intent(in)    :: sim

  type(type_node) :: node
  
  integer :: q, i, i_patch
  real*8 :: vector_normal(3), cos_alpha, mass_ion, c_s, Gamma_d, Gamma_convective
  real*8 :: T_i, T_e, n_e, n_e_corr, yield, vpar, grad_par_n_e
  real*8, dimension(3) :: E, B, B_hat, grad_n_e
  real*8 :: psi, U, psi_norm
  real*8 :: D_par_si, D_perp_si, D_norm !< normalisation for diffusivity: D_si = D_jor * D_norm
  real*8 :: c_angle !< min_sheath_angle but then in radians, same as in mod_boundary_matrix_open

  real*8 :: psi_axis, R_axis, Z_axis, s_axis, t_axis, psi_xpoint(2), psi_limit, R_xpoint(2), Z_xpoint(2), s_xpoint(2), t_xpoint(2)
  integer :: i_elm_axis, ifail, i_elm_xpoint(2), ierr
  integer :: n_dirichlet_BC_nodes, n_sheath_BC_nodes, n_grad0_BC_nodes, inode, inodes(n_vertex_max)

  c_angle = min_sheath_angle * PI/180.d0

  D_norm = 1.d0/(sqrt(MU_ZERO*central_mass*central_density*1.d20*ATOMIC_MASS_UNIT))

  ! projection diagnostic
  ! Preparation (force my_id to 1 to suppress message)
  ! Note that this does not do proper time interpolation! We should probably
  ! have a proper function on the simulation to obtain those parameters
  ! for a rough estimate it will work however
  !t_xpoint = 0.d0
  !s_xpoint= 0.d0
  call find_axis(1,sim%fields%node_list,sim%fields%element_list,psi_axis,R_axis,Z_axis,i_elm_axis,s_axis,t_axis,ifail)

  if (xpoint) then
    call find_xpoint(1,sim%fields%node_list,sim%fields%element_list,psi_xpoint,R_xpoint,Z_xpoint,i_elm_xpoint,s_xpoint,t_xpoint,xcase,ifail)
    psi_limit  = psi_xpoint(1)
    if((xcase .eq. 2) .or. ((xcase .eq. 3) .and. (psi_xpoint(2) .lt. psi_xpoint(1)))) then
      psi_limit = psi_xpoint(2)
    end if
  else
    if (sim%my_id .eq. 0) then
      write(*,*) "WARNING: limiter config for sputtering unsupported, use at your own risk"
    end if
    psi_limit = 0.d0 ! not really supported
  end if

  ! resetting fluid yield integral scalars
  do i=1,size(this%fluid_yield_integral%patch,1)
    this%fluid_yield_integral%patch(i)%scalars = -1
  end do

  do i_patch = 1, size(this%fluid_yield_integral%patch,1) !< different parts of edge domain
#ifdef __GFORTRAN__
    !$omp parallel do default(shared) &
#else
    !$omp parallel do default(none) &
    !$omp shared(this, sim, gamma, &
    !$omp i_patch, central_mass, psi_axis, psi_limit, c_angle, bcs, D_par, D_norm) &
#endif
    !$omp private(i, n_e, n_e_corr, T_e, vpar, grad_n_e, E, B, psi, U, psi_norm, vector_normal, B_hat, cos_alpha, grad_par_n_e, q, T_i, mass_ion, c_s, Gamma_d, &
    !$omp         Gamma_convective, yield, inodes, inode, node, n_dirichlet_BC_nodes, n_sheath_BC_nodes, n_grad0_BC_nodes, &
    !$omp         D_par_si, D_perp_si, ierr) schedule(static)
    do i = 1, size(this%fluid_yield_integral%patch(i_patch)%xyz, 2) !< over all nodes
      call sim%fields%calc_NeTeTi(sim%time, this%fluid_yield_integral%patch(i_patch)%i_elm_jorek_edge(i), this%fluid_yield_integral%patch(i_patch)%st(:,i), &
        real(this%fluid_yield_integral%patch(i_patch)%xyz(3,i), 8), n_e_corr, T_e, T_i, n_e_raw=n_e, grad_raw_n_e=grad_n_e)

      call sim%fields%calc_vpar(sim%time, this%fluid_yield_integral%patch(i_patch)%i_elm_jorek_edge(i), this%fluid_yield_integral%patch(i_patch)%st(:,i), &
        real(this%fluid_yield_integral%patch(i_patch)%xyz(3,i), 8), vpar)
      
      call sim%fields%calc_EBpsiU(sim%time, this%fluid_yield_integral%patch(i_patch)%i_elm_jorek_edge(i), &
           this%fluid_yield_integral%patch(i_patch)%st(:,i), &
           real(this%fluid_yield_integral%patch(i_patch)%xyz(3,i), 8), &
           E, B, psi, U)
      
      !> normal vector calculation (wall_normal_vector gives inward pointing normal, so vector_normal is outward pointing)
      vector_normal = - wall_normal_vector(sim%fields%node_list, sim%fields%element_list, &
          this%fluid_yield_integral%patch(i_patch)%i_elm_jorek_edge(i), &
          this%fluid_yield_integral%patch(i_patch)%st(1,i), &
          this%fluid_yield_integral%patch(i_patch)%st(2,i))
      
      !alpha = acos( dot_product(vector_normal,NORM2(B,dim=1))) !< acos is in radians
      ! the flux is given by the velocity along B dot n
      B_hat = B/norm2(B)
      cos_alpha = dot_product(vector_normal,B_hat)

      !convective flux to the wall
      Gamma_convective = n_e * vpar * norm2(B) * cos_alpha
        
      q = 1 ! for calculation of sound speed
      mass_ion = central_mass* ATOMIC_MASS_UNIT !< now we use only the deuterium soundspeed
      ! c_s = sqrt((k_boltz/mass_ion)*(T_e + gamma * T_i)) ! m/s !< gamma *(Te+Ti) in model303 and 307
      c_s = sqrt((k_boltz/mass_ion)*(gamma * (T_i+T_e))) !< IF model =303 / 307
      !<TODO: test c_s is vpar0, as this should account for all models
            
      ! logic to find which boundary condition we should be applying on the particle side to match the fluid side
      inodes = sim%fields%element_list%element(this%fluid_yield_integral%patch(i_patch)%i_elm_jorek_edge(i))%vertex
      n_dirichlet_BC_nodes=0; n_sheath_BC_nodes=0; n_grad0_BC_nodes=0;
      do inode=1,n_vertex_max
        node = sim%fields%node_list%node(inodes(inode))
        if(node%boundary > 0 .and. node%boundary < 20) then
          if(bcs(node%boundary)%dirichlet%rho) then
            n_dirichlet_BC_nodes = n_dirichlet_BC_nodes + 1
          else if(bcs(node%boundary)%natural%rho) then
            n_sheath_BC_nodes    = n_sheath_BC_nodes    + 1
          else
            n_grad0_BC_nodes     = n_grad0_BC_nodes     + 1
          end if
        end if
      enddo

      Gamma_d = 0

      ! calculating the particle flux to the wall
      if (n_sheath_BC_nodes > 0) then ! apply sheath BC -> loss is only convective + c_angle term
        Gamma_d = Gamma_convective + n_e * c_s * c_angle
      else if (n_dirichlet_BC_nodes > 0) then ! apply dirichlet BC -> no c_angle term, but diffusive and convective terms are present
        psi_norm = get_psi_n(psi,z=real(this%fluid_yield_integral%patch(i_patch)%xyz(2,i),8))
        D_par_si = D_par * D_norm
        D_perp_si = get_dperp(psi_norm) * D_norm
        grad_par_n_e=dot_product(b_hat, grad_n_e)
        
        Gamma_d = Gamma_convective - ((D_par_si - D_perp_si)*grad_par_n_e*cos_alpha + D_perp_si*dot_product(grad_n_e,vector_normal))
        
      else if (n_grad0_BC_nodes > 0) then !no BC applied -> grad rho . n = 0, so only loss is convective
        Gamma_d = Gamma_convective
      else
        !$omp critical
        write(*,*) "ERROR: in mod_particle_wall_interaction, it should be impossible to end up in this place in the if, else if tree"
        !$omp end critical
        call MPI_Abort(MPI_COMM_WORLD, 32, ierr)
      end if
      
      ! Assume an impact angle of 0!
      ! need the abs here because we cheat using negative numbers to indicate D, T
      ! cap ionisation level to 4
      q = min(abs(this%fluid_Z), 4)
      select case(trim(this%type))
      case("wall recomb")
        yield = 1.d0 !<assuming complete wall saturation
      case("fluid sputter")
        yield = fluid_sputtering_yield(this%yield, T_e * K_BOLTZ/EL_CHG, q, 0.d0)
      case default
        call wrong_interaction_type(trim(this%type))
      end select

      this%fluid_yield_integral%patch(i_patch)%scalars(i,1) = Gamma_d * this%delta_t * yield !< particles / m^2 in this timestep

      if (this%do_wall_projection) then
        !associate (sc => this%wall_projection%patch(i_patch)%scalars) ! associate is nice to make more readable but cannot be used in OMP before version 4.5 (so not in OneAPI's OMP)
        ! These are all also multiplied by delta_t so we can make an average
        ! over the diagnostics period. Disregard the time in the units.
        this%wall_projection%patch(i_patch)%scalars(i, n_project_general+1) = &
        this%wall_projection%patch(i_patch)%scalars(i, n_project_general+1) + n_e * this%delta_t ! n_e [m^-3]

        this%wall_projection%patch(i_patch)%scalars(i, n_project_general+2) = &
        this%wall_projection%patch(i_patch)%scalars(i, n_project_general+2) + T_e * K_BOLTZ / EL_CHG * this%delta_t ! T_e [eV]

        this%wall_projection%patch(i_patch)%scalars(i, n_project_general+3) = &
        this%wall_projection%patch(i_patch)%scalars(i, n_project_general+3) + cos_alpha * this%delta_t ! dimensionless, cosine of angle between wall normal and fieldline B

        this%wall_projection%patch(i_patch)%scalars(i, n_project_general+4) = &
        this%wall_projection%patch(i_patch)%scalars(i, n_project_general+4) + (psi - psi_axis)/(psi_limit - psi_axis) * this%delta_t ! normalized psi, dimensionless
        
        ! incoming fluid projections, should be similar to incoming particle projections, so this is like a sanity check
        ! number of particles incoming
        this%wall_projection%patch(i_patch)%scalars(i, n_project_general+5) = &
        this%wall_projection%patch(i_patch)%scalars(i, n_project_general+5) + Gamma_d * this%delta_t !< incident particle flux (particles/m^2)

        ! incident energy integrated over delta_t
        ! where we assume the ion energy to be 2 k T_i + 3 q k T_e as in the ! sputtering calculation above
        ! J/m^2
        this%wall_projection%patch(i_patch)%scalars(i, n_project_general+6) = &
        this%wall_projection%patch(i_patch)%scalars(i, n_project_general+6) + &
          Gamma_d * this%delta_t * (2.d0 * k_boltz * T_i + 3.d0 * k_boltz * q * T_e)

        ! sputtering yield in this time interval at this location [particles/m^2]
        this%wall_projection%patch(i_patch)%scalars(i, n_project_general+7) = &
        this%wall_projection%patch(i_patch)%scalars(i, n_project_general+7) + &
          Gamma_d * this%delta_t * yield
        !end associate
      end if
    end do
    !$omp end parallel do
  end do
end subroutine project_sputter_vars_on_edge


!> Sample the energy of a particle with charge Z_ion in the plasma (before the sheath)
!> from the local temperature
!>
!> For the ions treated as a fluid we make the assumption that they travel at the
!> background plasma sound speed.
!>
!> The energy is determined by the criterion $<v_par> > c_s$ with $c_s$ the background
!> plasma sound speed. That leads to the factor sqrt(m_ion/central_mass)
!>
!> The calculation here proceeds as follows:
!> 1. Calculate total energy from chi_squared(3) distribution, equal to chi_squared(2) + chi_squared(1)
!> 2. Calculate ratio of perpendicular and total energies muB = perp/(perp + par) (see https://en.wikipedia.org/wiki/Chi-squared_distribution#Relation_to_other_distributions)
!> 3. Calculate new parallel energy from the square of E +- cs, with + or - 50/50
!> 4. Add all energies together
!> 5. Correct for atomic weight, assuming all velocities are central_mass velocities
pure subroutine sample_fluid_particle_energy(Te_eV, u, Z_ion, E, E_threshold, Ti_eV)
  use phys_module, only: central_mass, gamma
  use mod_sampling, only: sample_chi_squared_3
  use mod_atomic_elements, only: atomic_weights

  real*8, intent(in)             :: Te_eV !< Electron temperature in eV
  real*8, intent(in)             :: u(3) !< random numbers for sampling
  integer, intent(in)            :: Z_ion
  real*8, intent(out)            :: E !< Energy in eV
  real*8, intent(in), optional   :: E_threshold !< Theshold energy in eV, not to sample particles below this energy
  real*8, intent(in), optional  :: Ti_eV !< ion temperature (eV), optional

  real*8                         :: beta, v, Ti_eV_local, fact

  ! Use Ti if provided, otherwise fall back to Te
  if (present(Ti_eV)) then
    Ti_eV_local = Ti_eV
  else
    Ti_eV_local = Te_eV
  end if

  ! Sample an energy at the local temperature
  E = Ti_eV_local*0.5d0*sample_chi_squared_3(u(1)) ! in eV

  ! Solve now for u = 1-sqrt(1-x) (CDF of beta(1,1/2) distribution)
  beta = 2.d0*u(2)-u(2)**2
  ! this is also the ratio between perpendicular and total energies
  ! the parallel energy is then given by
  ! E*(1-beta)
  ! and we take the square root of that to get a parallel velocity
  ! the direction of this is either + or - with 50/50 probability.
  ! Add the soundspeed (positive) to this and calculate the new energy
  ! v = sqrt(2E/m) (+ or - with 50/50 prob)
  v = sign(sqrt(2.d0*E*EL_CHG*(1.d0-beta)/(central_mass*ATOMIC_MASS_UNIT)), u(3)-0.5d0) ! m/s
  ! the sound speed is sqrt(k (1+gamma) T/m) = sqrt(T_eV*EL_CHG/m)
  ! gamma=1 is assumed, valid for cold dense plasma
  v = v + sqrt(gamma*(Te_eV+Ti_eV_local)*EL_CHG/(central_mass*ATOMIC_MASS_UNIT)) ! m/s
  E = E * beta + 0.5d0 * central_mass*ATOMIC_MASS_UNIT * v**2 / EL_CHG

  E = E*sqrt(atomic_weights(Z_ion)/central_mass) ! correct for atomic weight
end subroutine sample_fluid_particle_energy


!> find in which edge element patch index the particle is
function elm_in_patch(i_elm, edge_element_obj) result(i_patch)
  implicit none
  integer, intent(in) :: i_elm
  type(edge_elements), intent(in) :: edge_element_obj
  integer :: i_patch

  logical :: found

  found = .false.
  i_patch = -1

  do i_patch = 1,size(edge_element_obj%patch,1)
    ! if i_elm in the i_elm list of this edge domain exit the loop
    ! Note that this has issues at sharp corners, where particles may be
    ! lost in a different patch but at the same element number!
    if (any(i_elm .eq. edge_element_obj%patch(i_patch)%i_elm_jorek_edge(:))) then
      found = .true.
      exit
    endif
  end do
  ! i_patch should now be the first patch with correct element number, unless it wasn't found
  
  if (.not. found) then
    i_patch = -1 ! impossible number
    return
  end if

end function elm_in_patch


!> adds this particle's contribution to the wall_projection diagnostic tool
subroutine particle_projection_diagnostic(this, sim, particle, E, sputtering_yield)
  use phys_module, only: n_period, n_plane
  use mod_edge_elements, only: find_edge_element

  implicit none

  class(wall_action),              intent(inout) :: this
  class(particle_sim),             intent(in)    :: sim
  type(particle_kinetic_leapfrog), intent(in)    :: particle !< particle to undergo interaction
  real*8,                          intent(in)    :: E !< old energy of particle in eV
  real*8,                          intent(in)    :: sputtering_yield
  
  integer :: k, i_patch
 
  !> Prompt loss calculation
  integer :: is_prompt_loss
  real*8 :: Efield(3), B(3), pot, psi
  
  integer :: i_edge_elm, i_edge_nodes(4)
  real*8 :: area(4), dphi
  !> for mpi_reduce of particle contributions
  integer :: toroidal_offset !< Number of elements in the toroidal direction
  
  if(.not. this%do_wall_projection) return

  !> find in which patch the particle is lost
  i_patch = elm_in_patch(particle%i_elm, this%fluid_yield_integral)
  if (i_patch < 0) then
    !$omp critical
    write(*,"(A,I8,5es15.5)") "ERROR: in particle_self_reflection elm_in_patch, particle lost to somewhere unknown i_elm,s,t,R,Z,phi",particle%i_elm,particle%st,particle%x
    !$omp end critical
    return
  end if

  !> Write several diagnostics for the particle-particle sputtering
    ! the projection of a variable into the edge elements is simply a weighted addition to four points around an element
    ! Calculate the weight factors first and then store the relevant diagnostics
    ! find the corner point of the edge element we'll add the diagnostics to
  i_edge_elm = find_edge_element(this%wall_projection%patch(i_patch), particle%i_elm, particle%st(1), particle%st(2), particle%x(3))
  if (i_edge_elm .le. 0) then
    !$omp critical
    write(*,*) "ERROR: cannot find edge element for particle lost in this patch", particle%i_elm, particle%x(1), particle%x(2), i_edge_elm
    ! call flush(6)
    !$omp end critical
    return
  end if
  ! the weighting is done by inverse area
  ! 3-------|-----------4
  ! |   2   |   k=1     |
  ! |       |           |
  ! --------X------------
  ! |   4   |     3     |
  ! 1-------|-----------2
  ! in real space. i.e. calculate for each of the four quadrants above the surface area of the element
  ! and give them a fraction opposite area / total each.
  !
  ! The integrals are simple, since the elements are linear. It is given by
  ! \[
  !   \int_{l_0}^{l_1} \int_{\phi_0}^{\phi_1} R dl dphi
  ! \]
  ! The phi-integral drops out since it does not depend on l (they are orthogonal)
  ! and the other integral can be simplified since dl is along a straight line.
  ! this has as answer: 
  ! \[
  !   \left(r_0 l + \frac{1}{2} l^2 \frac{dr}{dl}\right) * (\phi_1 - \phi_0)
  ! \]
  ! with dr/dl = delta r / delta l (i.e. bounded between 0 and 1), 1 for purely outwards.
  !this%wall_projection%patch(i_patch)%scalars(index_node,5) = E * particle%weight
  toroidal_offset = this%wall_projection%patch(i_patch)%nsub_toroidal*n_plane
  if (toroidal_offset .eq. 1) toroidal_offset = 0 ! special case for fully axisymmetric
  i_edge_nodes = [i_edge_elm, i_edge_elm+1, &
      i_edge_elm + toroidal_offset,  &
      i_edge_elm + toroidal_offset + 1]

  ! area = r_0 l + (r_1-r_0) l / 2 = (r_1 + r_0) l / 2
  ! The indices k are as above shown, i.e. of the area opposite the node
  ! this is related to the edge nodes as
  ! 1 <-> 4 and 2 <-> 3, so 5-i
  do k=1,4
    if (i_edge_nodes(5-k) .gt. size(this%wall_projection%patch(i_patch)%xyz(1,:))) then
      write(*,*) "ERROR indexing problem in mod_wall_actioning",k,i_edge_elm,toroidal_offset, i_edge_nodes(5-k), size(this%wall_projection%patch(i_patch)%xyz(1,:))
      write(*,*) "ERROR temporary fix: set i_edge_nodes(5-k) = 1"
      i_edge_nodes(5-k) = 1
    end if

    area(k) = (this%wall_projection%patch(i_patch)%xyz(1,i_edge_nodes(5-k)) + particle%x(1)) &
       * norm2(this%wall_projection%patch(i_patch)%xyz(1:2,i_edge_nodes(5-k))-particle%x(1:2), dim=1) * 0.5d0
  end do
  ! multiply with delta-phi part
  ! we assume below that the particle is in this element (as it came from find_edge_element)
  dphi = TWOPI / (n_period * n_plane)
  area(1:2) = area(1:2) * modulo(dphi - particle%x(3), dphi) ! distance from X to top row
  area(3:4) = area(3:4) * modulo(particle%x(3) - dphi, dphi) ! distance from X to bottom row

  ! Multiply by this below (I might be guilty of some premature optimization here)
  is_prompt_loss = 0
  call sim%fields%calc_EBpsiU(sim%time, particle%i_elm, particle%st, particle%x(3), Efield, B, psi, pot)
  ! If the age of this particle is less than an a gyroperiod at the local magnetic field strength
  ! this particle is considered a prompt loss and will be written down below
  if ((sim%time - particle%t_birth) .lt. TWOPI * sim%groups(this%target_group)%mass*ATOMIC_MASS_UNIT/(EL_CHG * norm2(B))) is_prompt_loss = 1

  ! we need to loop here since omp atomic cannot set an array at once
  !associate (sc => this%wall_projection%patch(i_patch)%scalars) ! associate is nice to make more readable but cannot be used in OMP before version 4.5 (so not in OneAPI's OMP)
  do k=1,4
    ! particle flux
    ! (weight/n_period since we are only looking at the flux of one 1/n_period wedge)
    !$omp atomic
    this%wall_projection%patch(i_patch)%scalars(i_edge_nodes(k),1) = &
    this%wall_projection%patch(i_patch)%scalars(i_edge_nodes(k),1) + (particle%weight/n_period) * area(k)/sum(area)**2
    
    ! particle heat flux on edge elements (including sheath potential)
    ! (weight/n_period since we are only looking at the flux of one 1/n_period wedge)
    !$omp atomic
    this%wall_projection%patch(i_patch)%scalars(i_edge_nodes(k),2) = &
    this%wall_projection%patch(i_patch)%scalars(i_edge_nodes(k),2) + (particle%weight/n_period) * E * EL_CHG * area(k)/sum(area)**2
    
    ! particle flux from prompt redeposition (i.e. from particles younger than 2 pi / omega_c)
    ! (weight/n_period since we are only looking at the flux of one 1/n_period wedge)
    !$omp atomic
    this%wall_projection%patch(i_patch)%scalars(i_edge_nodes(k),3) = &
    this%wall_projection%patch(i_patch)%scalars(i_edge_nodes(k),3) + (particle%weight/n_period) * is_prompt_loss * area(k)/sum(area)**2
    
    ! sputtering yield
    ! (weight/n_period since we are only looking at the flux of one 1/n_period wedge)
    !$omp atomic
    this%wall_projection%patch(i_patch)%scalars(i_edge_nodes(k),4) = &
    this%wall_projection%patch(i_patch)%scalars(i_edge_nodes(k),4) + (particle%weight/n_period) * sputtering_yield * area(k)/sum(area)**2
    
  end do
  !end associate
end subroutine particle_projection_diagnostic


subroutine write_wall_project_vtk(this, sim)
  use mpi_mod

  implicit none

  type(wall_action),  intent(inout) :: this
  type(particle_sim), intent(in)    :: sim
  
  integer :: nnos, i, ierr
  real*4, allocatable :: scalars(:,:) !< for mpi_reduce of particle contributions
  character(len=120)  :: filename

  ! if not initialised, setting initial value of this%last_diag_time
  if (this%last_diag_time < 0) then
    this%last_diag_time = sim%time - this%delta_t
  end if

  ! determine filename
  if (len_trim(this%filename) .eq. 0) then
    filename = this%get_filename(sim%time)
  else
    filename = this%filename
  end if
  
  ! time normalising and MPI reducing the quantities
  do i = 1,size(this%wall_projection%patch,1)
    ! Turn all quantities from fluences into fluxes by dividing by the time since the last diagnostics output
    ! some of these (like T_e and n_e) were actually not fluences, but
    ! multiply those by this%delta_t anyway so this normalisation works and we get
    ! a decent time average
    this%wall_projection%patch(i)%scalars(:,:) = &
        this%wall_projection%patch(i)%scalars(:,:) / real(sim%time - this%last_diag_time,4)

    nnos = size(this%wall_projection%patch(i)%scalars,1)
    if (sim%my_id .eq. 0) then
      allocate(scalars(nnos,this%n_project_part))
    else
      allocate(scalars(0,0))
    end if
    
    ! Calculate the sum across mpi procs
    ! this needs to be done for all particle-quantities only (i.e. 1:this%n_project_part)
    call MPI_Reduce(this%wall_projection%patch(i)%scalars(:,1:this%n_project_part), &
        scalars, &
        nnos*this%n_project_part, MPI_REAL4, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    
    if (sim%my_id .eq. 0) then
      this%wall_projection%patch(i)%scalars(:,1:this%n_project_part) = scalars
    end if
    
    deallocate(scalars)
  end do 

  ! writing the projection
  if (sim%my_id .eq. 0) write(*,*) 'Writing wall projection diagnostics to ', trim(filename)
  call this%wall_projection%write_vtk_projection(filename)

  ! Reset diagnostic
  do i = 1,size(this%wall_projection%patch,1)
    this%wall_projection%patch(i)%scalars(:,:) = 0.d0
  end do
  this%i_step_diag = 0
  this%last_diag_time = sim%time
end subroutine write_wall_project_vtk


!> centralised routine to write out what input was wrong, and then stop the program
subroutine wrong_input(message, my_id, identifier)
  implicit none

  character(len=*), intent(in) :: message
  integer,          intent(in) :: my_id
  character(len=*), intent(in) :: identifier

  if(my_id > -1) then
    write(*,"(A,I2,5A)") "ERROR: (MPI ID=",my_id,") ",trim(message), " ", trim(identifier), " Exiting."
  end if

  stop
end subroutine wrong_input


!> centralised routine to throw the error of unsupported type and exit (please change this when you add a new type)
subroutine wrong_interaction_type(type, identifier)
  implicit none
  character(len=*), intent(in) :: type !< type which is not supported (will be trimmed in this subroutine)
  character(len=*), intent(in), optional :: identifier

  write(*,"(A)") "Wall interaction type "//trim(type)//" not supported (mod_wall_interaction.f90)"
  write(*,"(A)") 'Available types: "self sputter", "fluid sputter", "reflection" or "wall recomb" '
  if(present(identifier)) write(*,"(2A)") "Error detected ",trim(identifier)
  stop
end subroutine wrong_interaction_type


!> exit with message if origin_group is not target_group
subroutine check_self_type(this, my_id, identifier)
  implicit none
  type(wall_action), intent(in) :: this
  integer,           intent(in) :: my_id
  character(len=*),  intent(in) :: identifier
  
  character(len=1000) :: msg

  if(this%origin_group /= this%target_group) then
    write(msg,"(A)") "type "//trim(this%type)//" is a self interaction type so origin_group should be target_group. If you want you can leave the target_group_id blank for self interaction types, but you cannot specify another species. Please check your input"
    call wrong_input(msg, my_id, identifier)
  end if
end subroutine check_self_type


!> MPI reduces and writes the normal global diagnostics
!> returns the MPI reduced diagnostics back into diagnostics
subroutine write_global_diag(this,sim,diagnostics)
  implicit none

  type(wall_action),   intent(in) :: this
  type(particle_sim),  intent(in) :: sim
  real*8, dimension(n_global_diagnostics), intent(inout) :: diagnostics !< diagnostics for the global wall loads
  
  integer :: ierr

  call reduce_global_diag(diagnostics)

  ! write standard diagnostics to logfile
  if (sim%my_id .eq. 0) then
    write(*,'(A,2f14.0)' ) "superparticles going (in/out) = ", diagnostics(i_wall_part_in),             diagnostics(i_wall_part_out) 
    write(*,'(A,2es16.6)') "particle flux (in/out) [#/s]  = ", diagnostics(i_wall_flux_in)/this%delta_t,diagnostics(i_wall_flux_out)/this%delta_t 
    write(*,'(A,2es16.6)') "heatflux (in/out) [W]         = ", diagnostics(i_wall_heat_in)/this%delta_t,diagnostics(i_wall_heat_out)/this%delta_t 
    if(trim(this%type) == "pump") then
      write(*,'(A,1f14.0)' ) "pumped superparticles         = ", diagnostics(i_super_killed)
      write(*,'(A,1es16.6)') "pumped particles this step    = ", diagnostics(i_removed)
      write(*,'(A,1es16.6)') "pumped particles flux [#/s]   = ", diagnostics(i_removed)/this%delta_t
    endif
  endif

end subroutine write_global_diag


!> MPI reduces the global diagnostics and returns the reduced version
subroutine reduce_global_diag(diagnostics)
  use mpi_mod

  implicit none

  real*8, dimension(n_global_diagnostics), intent(inout) :: diagnostics !< diagnostics for the global wall loads
  
  real*8, dimension(n_global_diagnostics)                :: diagnostics_all_mpi !< MPI reduced diagnostics for the global wall loads
  integer :: ierr

  ! MPI reduce can be done at once for all diagnostics  
  call MPI_REDUCE(diagnostics, diagnostics_all_mpi, n_global_diagnostics, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  
  ! write MPI reduced diagnostics back to diagnostics
  diagnostics = diagnostics_all_mpi
end subroutine reduce_global_diag


!> first trims and then replaces spaces by underscores (_) in the string
function spaces2underscore(string_in) result(string_out)
  implicit none

  character(len=*), intent(in)  :: string_in !< string in which to replace spaces for _
  character(len=:), allocatable :: string_out

  character(len=:), allocatable :: string
  integer :: i

  string = trim(string_in)

  do i=1,len(string)
    if(string(i:i) == " ") string(i:i) = "_"
  end do

  string_out = string
end function

end module mod_particle_wall_interaction
