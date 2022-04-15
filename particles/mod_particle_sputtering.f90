!> Module for particle-particle and fluid-particle sputtering calculations.
!> This module will read where the simulated particles hit the wall.
!> and it calculates the sputtering from that event on that point

!> in delta_t time some amount of simulated particles hit the wall.
!> calculate on every point the amount of incident particles
!> calculate on the points of incident particles the amount of sputtered particles

!>
!> For the fluid-particle sputtering there are basically 2 approaches to follow
!> the first one is to sample particles from the incoming flux of every species
!> and then calculate the particle-particle sputtering yields. The second one is
!> to calculate the sputtering yield by integrating over the velocity distribution
!> and then sampling first in sputtering yield the spatial location of the particle
!> and then by sampling from the velocity distribution * sputtering yield in some way.
!>
!> The sputtering yield at a specific position is given by
!> \[
!>   \int_v Y(E) f(v) dv
!> \]
!> where $f(v)$ is a maxwellian and $E$ includes the sheath potential and the Bohm outflow
!> condition additionally.
!>
!> To now calculate the energy of the sputtered particle we multiply the sputtered energy
!> coefficient with E of a particle sampled from f(v).
!> Taking the sputtered energy coefficient * Y as a weight factor and sampling from f(v) will do the trick.
!> we need to normalize with the sputter yield at that position, which we have calculated before.
!> Basically this is a weighted average of Y_E(E) * E, weighted with Y(E) f(E).
!> Alternatively we could sample directly from Y(E) Y_E(E) f(E), but I don't know how to do this generally.
!> That would have the advantage of better distribution of statistics (more uniform weights).


!> Sputtering:
!> search for the number of lost particles on all mpi procs
!> Calculate sputter yield and amount of sputtered particles
!>   per particle species
!>   incoming energy and angle leads to sputter yield (and sputtered energy coefficient)
!>   distribute the to-sputter particles proportionally among processes
!>
!>
!> The module comes with a few diagnostics for sputtering
!> particle-particle sputtering:
!>   per incoming group:
!>     * particle flux on edge elements
!>     * particle energy flux on edge elements
!>     * sputtering yield
!> fluid-particle sputtering:
!>   per background species:
!>     * particle flux on edge elements
!>     * energy flux on edge elements
!>     * sputtering yield
!> global:
!>   * T_e [eV]
!> total:
!>   * overall sputtering yield (= sum yields above)
!> leading to 3*n_groups + 3*n_fluids + 1
module mod_particle_sputtering
  use mod_edge_elements
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
  !$ use omp_lib 
  
  implicit none
   
  private
  public :: particle_sputter
  public :: fluid_sputtering_yield !< only public for testing
  public :: sample_fluid_particle_energy

  !> The particle_sputter action encompasses sputtering coefficients and diagnostics
  !> as the main entry point to the sputtering tools. It needs to be initialised with
  !> an edge_domain and a target_group.
  type, extends(io_action) :: particle_sputter
    integer :: target_group !< which group to deposit sputtered particles in
  
    type(eckstein_sputter_yield), allocatable          :: yield(:) !< eckstein coefficients for particle-particle and fluid-particle sputtering
    type(eckstein_sputtered_energy_coeff), allocatable :: energy(:)
    type(thompson_dist) :: E_dist = thompson_dist(E_b = 8.7d0, n=2) !< produces energies in eV (value for W default)
    logical :: use_thompson = .false. !< Use a thompson distribution for the energy of sputtered particles
    logical :: use_Yn_func = .true. !< Use Ecksteins interpolating functions instead of interpolating manually
    
    class(type_rng), dimension(:), allocatable :: rng !< one RNG per openmp thread
    
    type(edge_elements) :: diagnostics
    type(edge_elements) :: fluid_sputter_yield !< the yield integrated over f(v) for every fluid background species
        
    real*8, dimension(:), allocatable  :: background_relative_density
    integer, dimension(:), allocatable :: background_species_Z
    
    integer :: i = 0, n_save = 100 !< used for diagnostics to see when the deposition diagnostic has to be evaluated, written
    integer :: n_sputter = -1 !< number of simulation particles to sputter from fluid-particle sputtering across all processes, across
    !< all fluid species (proportioned by sputtering yield, to get similar weights)

    real*8 :: last_time = 0.d0 !< When did we sputter last (used for integrating fluid flux and for calculating particle fluxes)
    real*8 :: last_diag_time = 0.d0 !< Last time of output of diagnostics

    real*8 :: sputtered_particle_weight_threshold = 1d7 !< Minimum weight of a macroparticle to be sputtered from particle-particle sputtering
  contains
    procedure :: do => do_particle_sputter
    procedure :: load_eckstein_data
  end type particle_sputter
  
  interface particle_sputter
    module procedure new_particle_sputter
  end interface particle_sputter


  ! Diagnostics settings
  integer, parameter :: n_particle_diag = 4 !< number of scalars to print for particle-particle sputtering diagnostics per particle group
  integer, parameter :: n_fluid_diag = 4 !< number of scalars to print for fluid-particle sputtering diagnostics per fluid group
  integer, parameter :: n_general_diag = 5 !< number of scalars to print as general diagnostics
contains

!> Constructor for the particle_sputter type, setting the io_action parameters and sputtering parameters.
function new_particle_sputter(edge_element_template, target_group, n_sputter, background_relative_density, &
        background_species_Z, filename, basename, decimal_digits, fractional_digits, rng, seed) result(new)
  use mod_pcg32_rng, only: pcg32_rng
  use mod_random_seed, only: random_seed
  use mpi_mod
  type(particle_sputter)          :: new
  type(edge_elements), intent(in) :: edge_element_template !< a prepared set of edge elements
  
  integer, intent(in) :: target_group !< which group to add particles to
  integer, intent(in) :: n_sputter !< number of particles to sputter across all processes (incl mpi)
  real*8, dimension(:), intent(in), optional  :: background_relative_density !< relative density of components of assumed divertor plasma
  integer, dimension(:), intent(in), optional :: background_species_Z !< Z of components of assumed divertor plasma
  
  character(len=*), intent(in), optional :: filename
  character(len=*), intent(in), optional :: basename
  integer, intent(in), optional          :: decimal_digits
  integer, intent(in), optional          :: fractional_digits
  class(type_rng), intent(in), optional  :: rng !< random-number generator to use (default PCG32)
  integer, intent(in), optional          :: seed !< Seed for the RNG (default random_seed() on my_id 0 + bcast)
  integer :: u, my_seed, n_stream, n_streams, ierr, i


  new%basename = "dep_"
  if (present(filename)) new%filename = filename
  if (present(basename)) new%basename = basename
  if (present(decimal_digits)) new%decimal_digits = decimal_digits
  if (present(fractional_digits)) new%fractional_digits = fractional_digits
  new%extension = '.vtk'
  new%name = "Particle sputter"
  new%log = .true.

  new%n_sputter = n_sputter
  
  !> make check if background_relative_density and background_species_Z are allocated
  !> otherwise make D plasma
  !> Allocating species for the sputtering from fluid
  new%target_group = target_group
  if (present(background_relative_density) .and. present(background_species_Z)) then
    allocate(new%background_relative_density(size(background_relative_density,1)))
    new%background_relative_density = background_relative_density
    allocate(new%background_species_Z(size(background_species_Z,1)))
    new%background_species_Z = background_species_Z
  end if
  
  if(.not. allocated(new%background_relative_density)) then
    allocate(new%background_relative_density(1),new%background_species_Z(1))
    new%background_species_Z = [-2]
    new%background_relative_density = [1]
    write(*,*) '===========================WARNING=================================='
    write(*,*) 'No plasma (fluid) species initialized, assume a deuterium plasma'
    write(*,*) '===================================================================='
  end if

  if (.not. allocated(edge_element_template%patch(1)%xyz)) then
    write(*,*) 'Edge element template needs to be prepared, exiting'
    stop
  end if
  
  new%diagnostics = edge_element_template
  new%fluid_sputter_yield = edge_element_template
  ! Clean up the passed edge elements
  do i=1,size(edge_element_template%patch,1)
    if (allocated(edge_element_template%patch(i)%scalars)) then
      deallocate(new%diagnostics%patch(i)%scalars, &
                 new%fluid_sputter_yield%patch(i)%scalars)
    end if
    if (allocated(edge_element_template%patch(i)%scalar_names)) then
      deallocate(new%diagnostics%patch(i)%scalar_names, &
                 new%fluid_sputter_yield%patch(i)%scalar_names)
    end if
  end do


  !> allocate random seed for sampling
  if (present(seed)) then
    my_seed = seed
  else
    my_seed = random_seed()
  end if
  if (present(rng)) then
    call setup_shared_rngs(n_dim=3, seed=my_seed, rng_type=rng, rngs=new%rng)
  else
    ! default to pcg32_rng for sputtering
    call setup_shared_rngs(n_dim=3, seed=my_seed, rng_type=pcg32_rng(), rngs=new%rng)
  end if
end function new_particle_sputter



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



!> Load eckstein sputtering yields and energy coefficients for all particle groups
!> and for all background species specified if not already loaded.
subroutine load_eckstein_data(this, sim)
  class(particle_sputter), intent(inout) :: this
  type(particle_sim), intent(in)         :: sim
  integer :: i, n_g

  if (.not. allocated(this%yield) .and. .not. allocated(this%energy)) then
    allocate(this%yield(size(sim%groups,1)+size(this%background_species_Z,1))) 
    allocate(this%energy(size(sim%groups,1)+size(this%background_species_Z,1)))
  
    !< reads eckstein sputter coefficients of all groups (such as tungsten)
    n_g = size(sim%groups,1)
    do i=1,n_g
      this%yield(i)%Z_ion    = sim%groups(i)%Z
      this%yield(i)%Z_target = sim%groups(this%target_group)%Z
      call this%yield(i)%read
      this%yield(i)%use_Yn_func = this%use_Yn_func
     
      if (.not. this%use_thompson) then
        this%energy(i)%Z_ion    = sim%groups(i)%Z
        this%energy(i)%Z_target = sim%groups(this%target_group)%Z
        call this%energy(i)%read
        this%energy(i)%use_Yn_func = this%use_Yn_func
      end if
    end do

    !< reads eckstein sputter coefficients for all background species
    do i=1,size(this%background_species_Z,1)
      this%yield(i+n_g)%Z_ion    = this%background_species_Z(i)
      this%yield(i+n_g)%Z_target = sim%groups(this%target_group)%Z
      call this%yield(i+n_g)%read
      this%yield(i+n_g)%use_Yn_func = this%use_Yn_func
      
      if (.not. this%use_thompson) then
        this%energy(i+n_g)%Z_ion    = this%background_species_Z(i)
        this%energy(i+n_g)%Z_target = sim%groups(this%target_group)%Z
        call this%energy(i+n_g)%read
        this%energy(i+n_g)%use_Yn_func = this%use_Yn_func
      end if
    end do
  end if 
end subroutine load_eckstein_data


!> Perform the sputtering on both fluid and particles hitting the wall.
!> Run this in an event every microsecond or less.
!>
!> Runs in 2 parts. First particle-particle sputtering and then particle-fluid
!> sputtering.
subroutine do_particle_sputter(this, sim, ev)
  use mpi_mod
  use mod_atomic_elements, only: element_symbols
  use mod_parameters, only: n_plane, n_period
  use mod_interp, only: interp_RZ
  class(particle_sputter), intent(inout) :: this
  type(particle_sim), intent(inout)      :: sim
  type(event), intent(inout), optional   :: ev

  integer :: n_fluid_groups, n_particle_groups

  integer :: i, j, k, i_patch, i_scalar, n_samples, ierr
  integer :: q, Z
  real*8 :: E, sputtering_yield, sputtered_energy_coeff, theta, T_eV, integral
  real*8, allocatable :: integral_i(:)
  real*8 :: velocity(3), delta_t, vector_normal(3)
  real*8, allocatable :: xyz_sampled(:,:), st_sampled(:,:), rng_sample(:,:) !< (3,n_samples) , (2,n_samples), (2,n_samples)

  !> Prompt loss calculation
  integer :: is_prompt_loss
  real*8 :: Efield(3), B(3), pot, psi

  character(len=120)  :: filename
  !> For check free particles
  integer, allocatable, dimension(:) :: i_free, i_elm_sampled
  logical, allocatable, dimension(:) :: is_free
  integer :: n_free
  !> For RNG
  real*8 :: u(2)
  integer :: i_rng, n_j, n_j_total
  real*8 :: n_e, T_e

  integer :: n_scalars, n, i_edge_elm, i_edge_nodes(4),i_p, nnos
  character(len=30) :: source_name, target_name, s_tmp
  real*8 :: area(4), av_yield, dphi
  !> For fluid sampling
  integer, allocatable, dimension(:) :: n_samples_fluid !< per background species
  !> for mpi_reduce of particle contributions
  real*4, allocatable :: scalars(:,:) !< size(st,1) by n_particle*3
  integer :: toroidal_offset !< Number of elements in the toroidal direction

  
  if (this%last_time .eq. 0.d0) then
    this%last_time = sim%time
    this%last_diag_time = sim%time
    if (abs(sim%time) .le. 1d-10) then
      this%last_time = 1d-10 ! so we actually do most of the sputtering should we start at 0
    end if
    return ! do nothing on the first run typically
  end if
  delta_t = sim%time - this%last_time
  this%last_time = sim%time

  n_fluid_groups = size(this%background_species_Z,1)
  n_particle_groups = size(sim%groups,1)

  ! load eckstein data if not already loaded
  call this%load_eckstein_data(sim)
  
  !> check if  edge domain is allocate to access the list of i_elm that are included in the sputtering
  !> If it is not allocated, don't sputter 
  if (.not. allocated(this%fluid_sputter_yield%patch)) then
    write(*,*)'=======================ERROR!!=================================='
    write(*,*)'fluid_sputter_yield not initialised. Did you use the constructor?'
    write(*,*)'exiting'
    return
  end if

  if (this%n_sputter .le. 0 .and. sim%my_id .eq. 0) then
    write(*,*) 'Warning: No sputtering quota set, calculating yields only'
  end if

  ! Set up scalars and scalar names for each of the patches if necessary
  n_scalars = n_particle_diag*n_particle_groups + n_fluid_diag*n_fluid_groups + n_general_diag
  do i=1,size(this%diagnostics%patch,1)
    if (.not. allocated(this%diagnostics%patch(i)%scalars)) then
      allocate(this%diagnostics%patch(i)%scalars(size(this%diagnostics%patch(i)%st,2), n_scalars))
      this%diagnostics%patch(i)%scalars(:,:) = 0.d0
    end if
    if (.not. allocated(this%diagnostics%patch(i)%scalar_names)) then
      allocate(this%diagnostics%patch(i)%scalar_names(n_scalars))
      associate (sn => this%diagnostics%patch(i)%scalar_names)
      n_j_total = count(sim%groups(:)%Z .eq. sim%groups(this%target_group)%Z) ! total number of groups of the target species
      if (n_j_total .gt. 1) then
        n_j = count(sim%groups(1:this%target_group)%Z .eq. sim%groups(this%target_group)%Z) ! number of groups of this species so far
        write(s_tmp,'(i1)') n_j
        target_name = trim(element_symbols(sim%groups(this%target_group)%Z))//trim(s_tmp)
      else
        target_name = trim(element_symbols(sim%groups(this%target_group)%Z))
      end if

      do j=1,n_particle_groups
        ! convention is source-target-type (without dashes), with numbers where necessary
        n_j_total = count(sim%groups(:)%Z .eq. sim%groups(j)%Z) ! total number of groups of this species
        if (n_j_total .gt. 1) then
          n_j = count(sim%groups(1:j)%Z .eq. sim%groups(j)%Z) ! number of groups of this species so far
          write(s_tmp,'(i1)') n_j
          source_name = trim(element_symbols(sim%groups(j)%Z))//trim(s_tmp)
        else
          source_name = trim(element_symbols(sim%groups(j)%Z))
        end if
        sn(n_particle_diag*j-3) = trim(source_name)//"flux"
        sn(n_particle_diag*j-2) = trim(source_name)//"heatflux"
        sn(n_particle_diag*j-1) = trim(source_name)//"promptflux"
        sn(n_particle_diag*j-0) = trim(source_name)//trim(target_name)//"yield"
      end do
      n = n_particle_diag*n_particle_groups !< offset
      do j=1,n_fluid_groups
        source_name = trim(element_symbols(this%background_species_Z(j)))//'f'
        sn(n+n_fluid_diag*j-3) = trim(source_name)//"density"
        sn(n+n_fluid_diag*j-2) = trim(source_name)//"flux"
        sn(n+n_fluid_diag*j-1) = trim(source_name)//"heatflux"
        sn(n+n_fluid_diag*j-0) = trim(source_name)//trim(target_name)//"yield"
      end do
      n = n_particle_diag*n_particle_groups + n_fluid_diag*n_fluid_groups
      sn(n+1) = "n_e" ! m^-3
      sn(n+2) = "T_e" ! eV
      sn(n+3) = "cos_alpha" ! cosine of angle between normal and B-field
      sn(n+4) = "Psi_n" ! normalized psi
      sn(n+5) = trim(target_name)//"_yield"
      end associate
    end if
  end do

  ! Allocate scalars for the fluid particle flux, one per background species
  do i=1,size(this%fluid_sputter_yield%patch,1)
    if (.not. allocated(this%fluid_sputter_yield%patch(i)%scalars)) then
      allocate(this%fluid_sputter_yield%patch(i)%scalars( &
        size(this%fluid_sputter_yield%patch(i)%st,2), n_fluid_groups))
    end if
    ! always reset this one to zero since its only used for the current iteration
    this%fluid_sputter_yield%patch(i)%scalars = 0
  end do
    
  
  

  !=============================================PARTICLE PART============================================================
  do i = 1,n_particle_groups ! source particles, i.e. those hitting the wall
    ! For each particle we need the location, the charge and the energy.
    ! instead of selecting type here we will loop first and use functions to get the charge and energy of the particle.
    ! the location requirement is fullfilled by particle_base already

    ! gfortran wants and does not want to have the types in the shared section at the same time.... default(shared) it is
    ! be very very careful however!
#ifdef __GFORTRAN__
    !$omp parallel default(shared) & 
#else
    !$omp parallel default(none) &
    !$omp shared(this, sim, i)   &
#endif
    !$omp private(q, velocity, theta, E, &
    !$omp        sputtering_yield, sputtered_energy_coeff, i_rng, u, i_patch,j, i_edge_nodes, vector_normal, T_eV, &
    !$omp        k, area, i_edge_elm, toroidal_offset, dphi, is_prompt_loss, Efield, B, psi, pot, T_e, n_e)

    i_rng = 1
    !$ i_rng = omp_get_thread_num()+1
    !$omp do schedule(dynamic, 10)
    do j = 1,size(sim%groups(i)%particles,1)
      ! Skip if this particle is not lost in a specific location (i_elm .eq. 0 means lost 'somewhere')
      if (sim%groups(i)%particles(j)%i_elm .ge. 0) cycle !< .not. .lt.!< if this is not a lost particle go to next particle
        
      ! Find out if this particle is lost in any of the edge domains
      do i_patch = 1,size(this%fluid_sputter_yield%patch,1)
        ! if i_elm in the i_elm list of this edge domain exit the loop
        if (any(-sim%groups(i)%particles(j)%i_elm .eq. this%fluid_sputter_yield%patch(i_patch)%i_elm_jorek_edge(:))) exit
        ! Note that this has issues at sharp corners, where particles may be
        ! lost in a different patch but at the same element number!
      end do
      if (i_patch .gt. size(this%fluid_sputter_yield%patch,1)) cycle ! particle not lost in the right area, skip it
        
        
      !> Place particle back into domain
      sim%groups(i)%particles(j)%i_elm = -sim%groups(i)%particles(j)%i_elm

      select type (pa => sim%groups(i)%particles(j))
      type is (particle_kinetic_leapfrog)
        velocity = pa%v
        q = pa%q
      class default
        write(*,*) "particle-particle sputtering pre-calc not implemented for this type, group=", i
        call exit(13)
      end select

      ! Calculate sputter yield, energy and write diagnostics
      associate (pa => sim%groups(i)%particles(j))
        ! use normal vector and velocity of particle to determine incoming angle
        ! cos(theta) = (n . v)/ (||n||.||v||)
        vector_normal = wall_normal_vector(sim%fields%node_list, sim%fields%element_list, pa%i_elm, pa%st(1), pa%st(2))
        theta = acos(dot_product(-vector_normal,velocity)/norm2(velocity))*180.d0/PI !< acos gives results in radians
        ! theta must be in degrees as the theta_star is also in degrees
        if (abs(theta) .gt. 91) then
          ! This is like an assert, it cannot really happen... but it does
          !!$omp critical
          !write(*,*) 'incoming angle warning', theta, vector_normal, velocity
          !!$omp end critical
        end if

        call sim%fields%calc_NeTe(sim%time, pa%i_elm, pa%st, pa%x(3), n_e, T_e)
        T_eV = T_e*K_BOLTZ/EL_CHG
        ! Energy associated with the velocity of the particle
        E = 0.5d0*sim%groups(i)%mass*ATOMIC_MASS_UNIT*dot_product(velocity, velocity)/EL_CHG !< must be in eV
          
        ! Update the particle energy from the potential drop in the sheath
        E = E + simple_potential_drop(q,T_eV)
          
          
          
        !> -------------Sputter yield---------------------------------------------------------------------     
        ! Hard-code theta to 0 to fix issues with sputtering module at strange angles
        ! the angle calculation should be revisited. Before using theta != 0 the
        ! surface roughness should be estimated, as this gives a distribution of
        ! impact angles as well
        theta = 0.d0 
        sputtering_yield = this%yield(i)%interp(E,theta)

        if (sputtering_yield .gt. 1) then
          !!$omp critical
          !write(*,"(A,f5.0,A,f8.3)") "> 1 self-sputtering detected, E=", E, "yield=", sputtering_yield
          !!$omp end critical
        end if

        !> Write several diagnostics for the particle-particle sputtering
        ! the projection of a variable into the edge elements is simply a weighted addition to four points around an element
        ! Calculate the weight factors first and then store the relevant diagnostics
        ! find the corner point of the edge element we'll add the diagnostics to
        i_edge_elm = find_edge_element(this%diagnostics%patch(i_patch), pa%i_elm, pa%st(1), pa%st(2), pa%x(3))
        if (i_edge_elm .le. 0) then
          !!$omp critical
          !write(*,*) "ERROR: cannot find edge element for particle lost in this patch", pa%i_elm, pa%x(1), pa%x(2), i_edge_elm
          !call flush(6)
          !!$omp end critical
          cycle
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
        !this%diagnostics%patch(i_patch)%scalars(index_node,5) = E * pa%weight
        toroidal_offset = this%diagnostics%patch(i_patch)%nsub_toroidal*n_plane
        if (toroidal_offset .eq. 1) toroidal_offset = 0 ! special case for fully axisymmetric
        i_edge_nodes = [i_edge_elm, i_edge_elm+1, &
            i_edge_elm + toroidal_offset,  &
            i_edge_elm + toroidal_offset + 1]

        ! area = r_0 l + (r_1-r_0) l / 2 = (r_1 + r_0) l / 2
        ! The indices k are as above shown, i.e. of the area opposite the node
        ! this is related to the edge nodes as
        ! 1 <-> 4 and 2 <-> 3, so 5-i
        do k=1,4
          area(k) = (this%diagnostics%patch(i_patch)%xyz(1,i_edge_nodes(5-k)) + pa%x(1)) &
             * norm2(this%diagnostics%patch(i_patch)%xyz(1:2,i_edge_nodes(5-k))-pa%x(1:2), dim=1) * 0.5d0
        end do
        ! multiply with delta-phi part
        ! we assume below that the particle is in this element (as it came from find_edge_element)
        dphi = TWOPI / (n_period * n_plane)
        area(1:2) = area(1:2) * modulo(dphi - pa%x(3), dphi) ! distance from X to top row
        area(3:4) = area(3:4) * modulo(pa%x(3) - dphi, dphi) ! distance from X to bottom row

        ! Multiply by this below (I might be guilty of some premature optimization here)
        is_prompt_loss = 0
        call sim%fields%calc_EBpsiU(sim%time, pa%i_elm, pa%st, pa%x(3), Efield, B, psi, pot)
        ! If the age of this particle is less than an a gyroperiod at the local magnetic field strength
        ! this particle is considered a prompt loss and will be written down below
        if ((sim%time - pa%t_birth) .lt. TWOPI * sim%groups(i)%mass*ATOMIC_MASS_UNIT/(EL_CHG * norm2(B))) is_prompt_loss = 1

        ! we need to loop here since omp atomic cannot set an array at once
        if (this%n_save .ge. 1) then
          do k=1,4
            ! these are n_particle_diag = 3 variables
            ! store on all 4 simultaneously, multiplied with weights
            ! * particle flux on edge elements
            !$omp atomic
            this%diagnostics%patch(i_patch)%scalars(i_edge_nodes(k),n_particle_diag*i-3) = &
            this%diagnostics%patch(i_patch)%scalars(i_edge_nodes(k),n_particle_diag*i-3) + pa%weight * area(k)/sum(area)**2
            ! * particle energy flux on edge elements (including sheath potential)
            !$omp atomic
            this%diagnostics%patch(i_patch)%scalars(i_edge_nodes(k),n_particle_diag*i-2) = &
            this%diagnostics%patch(i_patch)%scalars(i_edge_nodes(k),n_particle_diag*i-2) + pa%weight * E * area(k)/sum(area)**2
            ! * particle flux from prompt redeposition (i.e. from particles younger than 2 pi / omega_c)
            !$omp atomic
            this%diagnostics%patch(i_patch)%scalars(i_edge_nodes(k),n_particle_diag*i-1) = &
            this%diagnostics%patch(i_patch)%scalars(i_edge_nodes(k),n_particle_diag*i-1) + pa%weight * is_prompt_loss * area(k)/sum(area)**2
            ! * sputtering yield
            !$omp atomic
            this%diagnostics%patch(i_patch)%scalars(i_edge_nodes(k),n_particle_diag*i-0) = &
            this%diagnostics%patch(i_patch)%scalars(i_edge_nodes(k),n_particle_diag*i-0) + pa%weight * sputtering_yield * area(k)/sum(area)**2
          end do
        end if
        

        ! update weight of simulated particle after sputtering
        pa%weight = real(sputtering_yield,4) * pa%weight 
        !> If the weight gets too low, remove the particle
        if (pa%weight .le. this%sputtered_particle_weight_threshold) then ! this should really be an adaptive method based on something more rigorous
          pa%i_elm = 0
          cycle
        end if

        if (this%use_thompson) then
          call this%rng(i_rng)%next(u)
          E = sample_dist(this%E_dist, u(1))
        else
          sputtered_energy_coeff = this%energy(i)%interp(E,theta)
          E = sputtered_energy_coeff * E
        end if
      end associate
        
        
      !> -------------Set particle velocity/energy---------------------------------------------------------------------
      select type (pa => sim%groups(i)%particles(j))
      type is (particle_kinetic_leapfrog)
        ! give sputtered particle a new direction  
        ! use E from previous section to calculate velocity

        ! Calculate vector normal and select a random vector with a cosine distribution in angle between the normal and itself
        call this%rng(i_rng)%next(u)
        pa%v = sqrt(2.d0* E *EL_CHG/(sim%groups(i)%mass * ATOMIC_MASS_UNIT)) &
                * sample_cosine(u(1:2),vector_normal) 
        ! Since it is a neutral the half-step for boris method does not matter at all
        pa%q = 0_1
        pa%i_life = pa%i_life + 1
        pa%t_birth = sim%time
        ! For particle-particle sputtering we might want them to have the same identifiers
        ! if so comment the line above
        !<------------------------------------------------------------------------------------------------------------------
      class default
        write(*,*) "particle-particle sputtering post-calc not implemented for this type, group=", i
        call exit(13)
      end select

      associate (pa => sim%groups(i)%particles(j))
      if (any(pa%x .ne. pa%x) .or. E .ne. E .or. pa%weight .ne. pa%weight) then
        pa%i_elm = 0 ! skip this one since sputtering went wrong
        ! TODO debug logging? openmp threading though
      end if
      end associate
    end do
    !$omp end do
    !$omp end parallel
  end do




  !=============================================FLUID PART================================================================
  this%background_relative_density = this%background_relative_density/sum(this%background_relative_density, dim =1)
  
  !> this subroutine will calculate the incident ion flux over every fluid species on edge domain
  !> And the sputtering yield per group (in atoms/m^2 during delta_t)
  if (this%n_save .ge. 1) then
    call project_sputter_vars_on_edge(sim, this%background_relative_density, this%background_species_Z, &
        this%yield(n_particle_groups+1:n_particle_groups+n_fluid_groups), this%fluid_sputter_yield, delta_t, this%diagnostics)
  else
    call project_sputter_vars_on_edge(sim, this%background_relative_density, this%background_species_Z, &
        this%yield(n_particle_groups+1:n_particle_groups+n_fluid_groups), this%fluid_sputter_yield, delta_t)
  end if


  !===========================Check for free simulation particles=========================================
  allocate(is_free(size(sim%groups(this%target_group)%particles,1))) 

  ! might be replaced with omp workshare, or just the array expression.
  ! there is an issue with derived type arrays in gfortran though, and this works
  !$omp parallel do default(none) shared(sim, this, n_free, i_free, is_free, i) &
  !$omp private(j) schedule(dynamic, 100)
  do j=1,size(sim%groups(this%target_group)%particles,1)
    is_free(j) = sim%groups(this%target_group)%particles(j)%i_elm .le. 0 
  end do
  !$omp end parallel do
  !$omp barrier
  n_free = count(is_free)
  allocate(i_free(n_free))
  k = 1
  do j=1,size(is_free,1)
    if (is_free(j)) then
      i_free(k) = j
      k = k+1
    end if
  end do
  
  !=======================================================================================================
  ! share the work over all cpus. the last one takes the extra work
  if (this%n_sputter .ge. sim%n_cpu) then
    n_samples = min(this%n_sputter/sim%n_cpu,n_free)
    if (sim%my_id .eq. sim%n_cpu) then
      n_samples = min(mod(this%n_sputter,sim%n_cpu) + this%n_sputter/sim%n_cpu, n_free)
    end if
  else
    if (sim%my_id .lt. this%n_sputter) then
      n_samples = 1
    else
      n_samples = 0
    end if
  end if
 
  if (n_samples .eq. n_free .and. this%n_sputter/sim%n_cpu .gt. n_free) then
    write(*,"(i3,A,i3,A,i8,A,i8)") sim%my_id, 'Warning: could not sputter requested ', sim%n_cpu, 'x ', this%n_sputter/sim%n_cpu, ", bounded to ", n_samples
  end if

  ! Distribute the number of tries by integrated fluid sputtering yield (so the
  ! weights are similar)
  ! Every group gets the same number of tries, since we assume the sputtering yields are relatively similar
  ! (we have a lower density of the impurities but they sputter more efficiently)
  ! the weights are normalised for the discreteness of the number of particles to sputter
  ! we normalize the weights with the integral of the sputter yield later

  allocate(n_samples_fluid(n_fluid_groups),integral_i(n_fluid_groups))
  allocate(rng_sample(3,0),xyz_sampled(3,0),st_sampled(2,0),i_elm_sampled(0))
  do i=1,n_fluid_groups
    ! Sampled particles are already distributed according to sputtering yield, we only need to calculate the energy
    call sample_edge_elements(this%fluid_sputter_yield, i, 0, rng_sample(1:2,0), &
        integral_i(i), xyz_sampled, st_sampled, i_elm_sampled)
    ! This puts the integral into n_samples_fluid
  end do
  n_samples_fluid = floor(n_samples*integral_i/sum(integral_i,1)) ! Normalize to sum 1, then multiply by number of samples to make
  n_samples_fluid(1) = n_samples_fluid(1) + (n_samples - sum(n_samples_fluid)) !  give the extra particles to the first group
  deallocate(rng_sample,xyz_sampled,st_sampled,i_elm_sampled,integral_i)

  k = 0 !< offset for index of i_free particles
  ! Loop over particle sources
  do i = 1, n_fluid_groups
    allocate(rng_sample(3,n_samples_fluid(i)))
    allocate(xyz_sampled(3,n_samples_fluid(i)))
    allocate(st_sampled(2,n_samples_fluid(i)))
    allocate(i_elm_sampled(n_samples_fluid(i)))

    ! We need to properly use all RNGS here to avoid missing numbers
    ! needs default(shared) for gfortran
    !$omp parallel default(none) &
    !$omp shared(this, rng_sample, i, n_samples_fluid) private(i_rng, j)
    i_rng = 1
    !$ i_rng = omp_get_thread_num()+1
    !$omp do schedule(static,1)
    do j=1,n_samples_fluid(i)
      call this%rng(i_rng)%next(rng_sample(:,j))
    end do
    !$omp end do
    !$omp end parallel

    q = this%background_species_Z(i)
    if (q .le. 0) q = 1 ! deuterium, tritium special case
    q = min(q, 4) ! limit to 4 for divertor conditions
    Z = this%background_species_Z(i)

    ! Sampled particles are already distributed according to sputtering yield, we only need to calculate the energy
    call sample_edge_elements(this%fluid_sputter_yield, i, n_samples_fluid(i), rng_sample(1:2,:), &
        integral, xyz_sampled, st_sampled, i_elm_sampled)

    if (integral .le. 1d-12) then
      deallocate(rng_sample, xyz_sampled, st_sampled, i_elm_sampled)
      cycle ! Move along, nothing to sputter here
    end if

    if (sim%my_id .eq. 0 .and. this%n_sputter .gt. 0) then
      write(*,"(A,i3,A,i8,A,A,A,A,A,i1,A,i2,A,g12.4)") "Sputtered ", sim%n_cpu, "x", n_samples_fluid(i), &
        " ", element_symbols(sim%groups(this%target_group)%Z),&
        " from ", element_symbols(Z), " in group ", this%target_group, &
      " (Z=", sim%groups(this%target_group)%Z, ") with total weight ", integral
    end if

#ifdef __GFORTRAN__
    !$omp parallel default(shared) &
#else
    !$omp parallel default(none) &
    !$omp shared(this, sim, i, k, rng_sample, xyz_sampled, st_sampled, i_elm_sampled, n_samples_fluid, i_free, &
    !$omp        integral, delta_t, q, Z, n_particle_groups) &
#endif
    !$omp private(i_rng, j, theta, E, sputtering_yield, av_yield, sputtered_energy_coeff, u, i_p, vector_normal, &
    !$omp         T_e, T_eV, n_e)
    i_rng = 1
    !$ i_rng = omp_get_thread_num()+1
    !$omp do schedule(static,1)
    do j=1,n_samples_fluid(i)
      i_p = i_free(k+j)
      sim%groups(this%target_group)%particles(i_p)%i_elm = i_elm_sampled(j)
      if (i_elm_sampled(j) .le. 0) cycle
      sim%groups(this%target_group)%particles(i_p)%st = st_sampled(:,j)
      call interp_RZ(sim%fields%node_list, sim%fields%element_list, i_elm_sampled(j), &
        st_sampled(1,j), st_sampled(2,j), &
        sim%groups(this%target_group)%particles(i_p)%x(1), &
        sim%groups(this%target_group)%particles(i_p)%x(2))
      sim%groups(this%target_group)%particles(i_p)%x(3) = xyz_sampled(3,j) ! phi coordinate from sampling
  
      !> weight of fluid particle is equally distributed as a fraction of the incoming flux. Such that the sum of all incoming fluid particles ,
      !> is the total amount of incoming particles over the edge domain area * delta_t
      sim%groups(this%target_group)%particles(i_p)%weight = real(1.d0/(n_samples_fluid(i)*sim%n_cpu) * integral,4)
  
      ! Assume the particle is fully ionized! For D, He, Ar this is reasonable
      theta = 0 !< surface roughness leads to an average incoming angle of zero

      ! Calculate temperature at this position
      call sim%fields%calc_NeTe(sim%time, i_elm_sampled(j), st_sampled(:,j), xyz_sampled(3,j), n_e, T_e)
      T_eV = T_e * K_BOLTZ / EL_CHG

      ! sample from energy distribution of the plasma on the edge of the plasma sheath
      ! do not sample energies lower than E_threshold, since they will not sputter anyways
      !call sample_fluid_particle_energy(T_eV, rng_sample(1:3,j), Z, E)
      E = 2 * T_eV
      ! add to this energy the plasma sheath potential
      E = E + simple_potential_drop(q, T_eV)
      
      ! If sampling from the incoming energy distribution function, the
      ! sputtered energy coefficient needs to be reweighed with the sputtering
      ! yield at this energy (since the tail contributes more)
      ! This is uncommented below since we have simplified the model for now to
      ! work at a fixed energy of 3 q T_e + 2 T_i, so we don't need to do this
      ! anymore. The extension to realistic IEDFs should be done later, so 
      ! I've kept some of the code around.

      ! Workaround if sampling leads to positions where the temperature is just
      ! too low we warn here, and set the energy higher
      ! This can happen when the cubic and linear elements mismatch in their
      ! estimation of T
      !if (E .le. this%yield(i + n_particle)%E_threshold) then
        !write(*,*) "Wrong location selected:", E, this%yield(i + n_particle)%E_threshold
      !end if
      !E = max(E, this%yield(i + n_particle)%E_threshold + 1d0) ! add little bit to prevent zeros

      !sputtering_yield = this%yield(i + n_particle)%interp(E, theta)
      ! Workaround if sputtered energy coeff threshold is lower than sputtering
      ! threshold: use sputtered energy coeff just above threshold instead
      ! (note: all this doesn't take into account theta properly)
      if (this%use_thompson) then
        call this%rng(i_rng)%next(u)
        ! Remove the highest 2% of the distribution by clipping u (hacky)
        u = min(u, 0.98d0)
        E = sample_dist(this%E_dist, u(1))
      else
        E = max(E, this%energy(i + n_particle_groups)%E_threshold + 1d0) ! add little bit to prevent zeros
        sputtered_energy_coeff = this%energy(i + n_particle_groups)%interp(E, theta)
        E = sputtered_energy_coeff * E
      end if
      !av_yield = fluid_sputtering_yield(this%yield(i + n_particle), T_eV, Z, theta)
      ! we could probably avoid the calculation of fluid_sputtering_yield by
      ! using the discretisation we just sampled from (if theta is constant)
      !if (av_yield .le. 1d-18) av_yield = 1d-6 ! does not matter since then sputtering_yield must be 0, just to avoid a NaN below
      
      ! now we weigh the particles with the prevalence of this energy in sputtered particles, i.e. sputtering_yield
      ! over the integral of sputtering_yield, which we calculate (again)
      !sim%groups(this%target_group)%particles(i_p)%weight = &
      !sim%groups(this%target_group)%particles(i_p)%weight * &
        !sputtering_yield / av_yield

      ! In the case that sputtering from fluid is really low, we discard it.
      ! we've put the threshold way lower than for self-bombardment to make sure the statistics of the sputtering yield,
      ! still is valid from the fluid perspective
      ! 10 mEv is very low already!
      if (sim%groups(this%target_group)%particles(i_p)%weight .le. 1.d4 .or. E .le. 1d-2) then ! if energy negative the sqrt below will cause trouble. If zero the particle will not enter the domain
        !write(*,*) 'weight or E limit ', E, sputtered_energy_coeff, sputtering_yield, av_yield
        sim%groups(this%target_group)%particles(i_p)%i_elm = 0
        cycle       
      end if
      

      ! Store the result
      select type (pa => sim%groups(this%target_group)%particles)
      type is (particle_kinetic_leapfrog)
        vector_normal = wall_normal_vector(sim%fields%node_list, sim%fields%element_list, pa(i_p)%i_elm, pa(i_p)%st(1), pa(i_p)%st(2))
        call this%rng(i_rng)%next(u)
        pa(i_p)%v = sqrt(2.d0* E *EL_CHG/(sim%groups(1)%mass * ATOMIC_MASS_UNIT)) &
                * sample_cosine(u(1:2), vector_normal)
        ! Since it is a neutral the half-step for boris method does not matter at all
        ! particles come back as neutrals
        pa(i_p)%q = 0_1
        pa(i_p)%i_life = pa(i_p)%i_life + 1 ! This is now really a new particle
        pa(i_p)%t_birth = sim%time
      class default
        write(*,*) 'Particle type not implemented as sputtered particle'
        stop
      end select

      ! NaN checks
      if (     (sim%groups(this%target_group)%particles(i_p)%x(1) .ne. sim%groups(this%target_group)%particles(i_p)%x(1)) &
          .or. (sim%groups(this%target_group)%particles(i_p)%x(2) .ne. sim%groups(this%target_group)%particles(i_p)%x(2)) &
          .or. (sim%groups(this%target_group)%particles(i_p)%x(3) .ne. sim%groups(this%target_group)%particles(i_p)%x(3)) &
          .or. (E .ne. E)                                                                                                 &
          .or. (sim%groups(this%target_group)%particles(i_p)%weight .ne. sim%groups(this%target_group)%particles(i_p)%weight) ) then
       
        sim%groups(this%target_group)%particles(i_p)%i_elm = 0 ! skip this one since sputtering went wrong
       
      end if
    end do
    !$omp end do
    !$omp end parallel

    k = k + n_samples_fluid(i)
    deallocate(rng_sample, xyz_sampled, st_sampled, i_elm_sampled)
  end do
  !<------------------------------------------------------------------------------------

  
  deallocate(i_free, is_free)
  
  !> Part for writing diagnostics
  if (len_trim(this%filename) .eq. 0) then
    filename = this%get_filename(sim%time)
  else
    filename = this%filename
  end if
  
  this%i = this%i+1
  if (this%i .ge. this%n_save) then
    do i = 1,size(this%diagnostics%patch,1)
      ! Turn all quantities from fluences into fluxes by dividing by the time since the last diagnostics output
      ! some of these (like T_e and n_e) were actually not fluences, but
      ! multiply those by delta_t anyway so this normalisation works and we get
      ! a decent time average
      this%diagnostics%patch(i)%scalars(:,:) = &
          this%diagnostics%patch(i)%scalars(:,:) / real(sim%time - this%last_diag_time,4)

      nnos = size(this%diagnostics%patch(i)%scalars,1)
      if (sim%my_id .eq. 0) then
        allocate(scalars(nnos,n_particle_groups*n_particle_diag))
      else
        allocate(scalars(0,0))
      end if
      ! And if necessary calculate the sum across mpi procs
      ! this needs to be done for all particle-quantities only (i.e. 1:n_particle_groups*n_particle_diag)
      call MPI_Reduce(this%diagnostics%patch(i)%scalars(:,1:n_particle_groups*n_particle_diag), &
          scalars, &
          nnos*n_particle_groups*n_particle_diag, MPI_REAL4, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      if (sim%my_id .eq. 0) then
        this%diagnostics%patch(i)%scalars(:,1:n_particle_groups*n_particle_diag) = scalars

        ! sum the sputter yields together to calculate the total
        ! which is at the last variable in diagnostics
        this%diagnostics%patch(i)%scalars(:,n_particle_groups*n_particle_diag + n_fluid_groups*n_fluid_diag + n_general_diag) = &
          ! Particle set
          sum(this%diagnostics%patch(i)%scalars(:,n_particle_diag*1:n_particle_diag*n_particle_groups:n_particle_diag), dim=2) + &
          ! Fluid set
          sum(this%diagnostics%patch(i)%scalars(:, &
            n_particle_diag*n_particle_groups+n_fluid_diag*1:&
            n_particle_diag*n_particle_groups+n_fluid_diag*n_fluid_groups:n_fluid_diag), dim=2)
      end if
      deallocate(scalars)
    end do 


    if (sim%my_id .eq. 0) write(*,*) 'Writing particle sputtering diagnostics to ', trim(filename)
    call this%diagnostics%write(filename)
    
    ! Reset everything to 0
    do i = 1,size(this%diagnostics%patch,1)
      this%diagnostics%patch(i)%scalars(:,:) = 0.d0
    end do
    
    this%i = 0
    this%last_diag_time = sim%time
  end if
end subroutine do_particle_sputter

!> Integrate the sputtering yield over the distribution of incoming velocities.
!> Since we use inverse transform sampling on u to calculate the energy we can
!> just integrate over u from 0 to 1 to cover the whole distribution.
!> Do this with n subelements, using gaussian quadrature in each element
!> the subintervals go from 1/2 to n-1/2 to avoid using 0 and 1, since 1 should
!> lead to infinity for sampling from a gaussian. Skipping the first part is reasonable
!> since the sputtering yield will be very low there. For the high energies a maxwellian
!> is perhaps not even the best approximation so that is probably not so bad either.
!>
!> The returned yield is averaged over the maxwellian at T_eV + the potential drop
pure function fluid_sputtering_yield(coeff, T_eV, Z, theta) result(yield)
  use gauss
  class(eckstein_coeff_set), intent(in) :: coeff
  real*8, intent(in)  :: T_eV !< Plasma temperature in eV
  integer, intent(in) :: Z !< Atomic number of the incoming particles
  real*8, intent(in)  :: theta !< angle of impact (usually assumed 0) in degrees
  real*8 :: yield !< The sputter yield in atoms/ion

  real*8 :: E, U_drop
  integer :: i, j, k
  integer, parameter :: n_interval = 4 !< number of intervals to calculate. (using 1 is already pretty good)
  real*8, parameter :: idu = 1.d0/real(n_interval,8) !< interval size
  real*8 :: u(3) !< the integration point
  integer :: q

  if (Z .le. 0) then
    q = 1
  else
    q = min(Z, 4) ! cap to 4 for divertor conditions
  end if 

  ! We use a simplified model for now! 2 T_i + 3 q T_e
  U_drop = simple_potential_drop(q, T_eV) ! assume particle has full charge
  yield = coeff%interp(2*T_eV + U_drop, theta)

  !!--------------------
  return
  !!--------------------

  yield = 0.d0
  if (T_eV .le. 1d-1) return
  do i=0,n_interval*n_gauss-1
    u(1) = (real(i/n_gauss,8) + xgauss(mod(i,n_gauss)+1))*idu
    do j=0,n_interval*n_gauss-1
      u(2) = (real(j/n_gauss,8) + xgauss(mod(j,n_gauss)+1))*idu
      do k=0,1 ! we only use the sign of this one
        u(3) = real(k,8)

        ! add to this energy the plasma sheath potential
        call sample_fluid_particle_energy(T_eV, u, Z, E)
        U_drop = simple_potential_drop(q, T_eV) ! assume particle has full charge
        yield = yield + coeff%interp(E + U_drop, theta)
      end do
    end do
  end do
  yield = yield / (2*n_interval**2)
end function fluid_sputtering_yield


!> Calculate the flux to and some diagnostics for fluid flux in a period delta_t
!>
!> fluid_sputter_yield contains a single scalar, the incoming fluid flux. Assume all particles are moving at the same
!> velocity, so multiplying with the relative density is enough to get the flux of a specific species
!>
!> Assume that the impact angle of all particles is 0
subroutine project_sputter_vars_on_edge(sim, n_relative, background_species, coeff, fluid_sputter_yield, delta_t, diagnostics)
  use mod_edge_elements, only: edge_elements
  use mod_atomic_elements, only: atomic_weights
  use phys_module, only: central_mass, xpoint, xcase
  use mod_parameters, only: n_plane, n_period

  type(particle_sim), intent(in)                     :: sim
  real*8, dimension(:), intent(inout)                :: n_relative !< array of the fraction of density of different species in the plasma
  integer, dimension(size(n_relative,1)), intent(in) :: background_species !< array of the Z of different species in the plasma
  type(eckstein_sputter_yield), intent(in)           :: coeff(size(n_relative,1)) !< eckstein sputtering yield coefficients
  type(edge_elements), intent(inout)                 :: fluid_sputter_yield !< the sputtering yield from species n
  real*8, intent(in)                                 :: delta_t
  type(edge_elements), intent(inout), optional       :: diagnostics !< diagnostics, saved if present

  integer :: q, i, j, i_patch, n_offset, Z
  real*8 :: vector_normal(3), cos_alpha, mass_ion, c_s, Gamma_d, n_species
  real*8 :: T_i, T_e, n_e, yield
  real*8, dimension(3) :: E, B, B_hat
  real*8 :: m, psi, U

  real*8, parameter :: gamma = 5.d0 / 3.d0 !< Heat capacity ratio, for adiabatic
  real*8 :: psi_axis, R_axis, Z_axis, s_axis, t_axis, psi_xpoint(2), psi_limit, R_xpoint, Z_xpoint, s_xpoint, t_xpoint
  integer :: i_elm_axis, ifail, i_elm_xpoint

  if (present(diagnostics)) then
    ! Preparation (force my_id to 1 to suppress message)
    ! Note that this does not do proper time interpolation! We should probably
    ! have a proper function on the simulation to obtain those parameters
    ! for a rough estimate it will work however
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
  end if

  ! normalize relative densities
  n_relative = n_relative/sum(n_relative, dim =1)
  
  do i_patch = 1, size(fluid_sputter_yield%patch,1) !< different parts of edge domain
#ifdef __GFORTRAN__
    !$omp parallel do default(shared) &
#else
    !$omp parallel do default(none) &
    !$omp shared(fluid_sputter_yield, sim, n_relative, background_species, coeff, diagnostics, delta_t, &
    !$omp        i_patch, central_mass, psi_axis, psi_limit) &
#endif
    !$omp private(i, n_e, T_e, E, B, psi, U, vector_normal, B_hat, cos_alpha, q, T_i, mass_ion, c_s, j, m, n_species, Gamma_d, &
    !$omp         n_offset, yield, Z) schedule(static)
    do i = 1, size(fluid_sputter_yield%patch(i_patch)%xyz, 2) !< over all nodes
      call sim%fields%calc_NeTe(sim%time, fluid_sputter_yield%patch(i_patch)%i_elm_jorek_edge(i), fluid_sputter_yield%patch(i_patch)%st(:,i), &
        real(fluid_sputter_yield%patch(i_patch)%xyz(3,i), 8), n_e, T_e)
      
      call sim%fields%calc_EBpsiU(sim%time, fluid_sputter_yield%patch(i_patch)%i_elm_jorek_edge(i), &
           fluid_sputter_yield%patch(i_patch)%st(:,i), &
           real(fluid_sputter_yield%patch(i_patch)%xyz(3,i), 8), &
           E, B, psi, U)
      
      !> normal vector calculation
      vector_normal = wall_normal_vector(sim%fields%node_list, sim%fields%element_list, &
          fluid_sputter_yield%patch(i_patch)%i_elm_jorek_edge(i), &
          fluid_sputter_yield%patch(i_patch)%st(1,i), &
          fluid_sputter_yield%patch(i_patch)%st(2,i))
      
      !alpha = acos( dot_product(vector_normal,NORM2(B,dim=1))) !< acos is in radians
      ! the flux is given by the velocity along B dot n
      B_hat = B/norm2(B)
      cos_alpha = abs(dot_product(vector_normal,B_hat))
        
      q = 1 ! for calculation of sound speed
      T_i = T_e !< not made for model 400 [K]
      mass_ion = central_mass* ATOMIC_MASS_UNIT !< now we use only the deuterium soundspeed
      c_s = sqrt((k_boltz/mass_ion)*(T_e + gamma * T_i)) ! m/s
      
      do j = 1, size(n_relative,1) !< over all fluid groups
        Z = background_species(j)
        m = atomic_weights(Z) * ATOMIC_MASS_UNIT
        n_species = n_e * n_relative(j)
        
        Gamma_d = n_species * c_s * cos_alpha

        ! Assume an impact angle of 0!
        ! need the abs here because we cheat using negative numbers to indicate D, T
        ! cap ionisation level to 4
        q = min(abs(background_species(j)), 4)
        yield = fluid_sputtering_yield(coeff(j), T_e * K_BOLTZ/EL_CHG, q, 0.d0)
        fluid_sputter_yield%patch(i_patch)%scalars(i,j) = Gamma_d * delta_t * yield !< particles / m^2 in this timestep

        n_offset = size(sim%groups,1)*n_particle_diag
        if (present(diagnostics)) then
          ! Particle density
          diagnostics%patch(i_patch)%scalars(i,n_offset+n_fluid_diag*j-3) = &
          diagnostics%patch(i_patch)%scalars(i,n_offset+n_fluid_diag*j-3) + n_species * delta_t !< particle density (particles/m^3*s) (i.e. time-integrated density)

          ! number of particles incoming
          diagnostics%patch(i_patch)%scalars(i,n_offset+n_fluid_diag*j-2) = &
          diagnostics%patch(i_patch)%scalars(i,n_offset+n_fluid_diag*j-2) + Gamma_d * delta_t !< incident particle flux (particles/m^2)

          ! incident energy integrated over delta_t
          ! where we assume the ion energy to be 2 k T_i + 3 q k T_e as in the ! sputtering calculation above
          ! J/m^2
          diagnostics%patch(i_patch)%scalars(i,n_offset+n_fluid_diag*j-1) = &
          diagnostics%patch(i_patch)%scalars(i,n_offset+n_fluid_diag*j-1) + &
            Gamma_d * delta_t * (2.d0 * k_boltz * T_i + 3.d0 * k_boltz * q * T_e)

          ! sputtering yield in this time interval at this location [particles/m^2]
          diagnostics%patch(i_patch)%scalars(i,n_offset+n_fluid_diag*j-0) = &
          diagnostics%patch(i_patch)%scalars(i,n_offset+n_fluid_diag*j-0) + &
            Gamma_d * delta_t * yield
        end if
      end do

      ! Save electron temperature
      if (present(diagnostics)) then
        n_offset = size(sim%groups,1)*n_particle_diag + size(n_relative,1)*n_fluid_diag
        ! These are all also multiplied by delta_t so we can make an average
        ! over the diagnostics period. Disregard the time in the units.
        diagnostics%patch(i_patch)%scalars(i,n_offset+1) = &
        diagnostics%patch(i_patch)%scalars(i,n_offset+1) + n_e * delta_t ! n_e [m^-3]

        diagnostics%patch(i_patch)%scalars(i,n_offset+2) = &
        diagnostics%patch(i_patch)%scalars(i,n_offset+2) + T_e * K_BOLTZ / EL_CHG * delta_t ! T_e [eV]

        diagnostics%patch(i_patch)%scalars(i,n_offset+3) = &
        diagnostics%patch(i_patch)%scalars(i,n_offset+3) + cos_alpha * delta_t ! dimensionless, angle between normal and fieldline

        diagnostics%patch(i_patch)%scalars(i,n_offset+4) = &
        diagnostics%patch(i_patch)%scalars(i,n_offset+4) + (psi - psi_axis)/(psi_limit - psi_axis) * delta_t ! normalized psi, dimensionless
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
pure subroutine sample_fluid_particle_energy(T_eV, u, Z_ion, E, E_threshold)
  use phys_module, only: central_mass
  use mod_sampling, only: sample_chi_squared_3
  use mod_atomic_elements, only: atomic_weights

  real*8, intent(in)             :: T_eV !< Temperature in eV
  real*8, intent(in)             :: u(3) !< random numbers for sampling
  integer, intent(in)            :: Z_ion
  real*8, intent(out)            :: E !< Energy in eV
  real*8, intent(in), optional   :: E_threshold !< Theshold energy in eV, not to sample particles below this energy

  real*8                         :: beta, v

  ! Sample an energy at the local temperature
  E = T_eV*0.5d0*sample_chi_squared_3(u(1)) ! in eV
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
  v = v + sqrt(T_eV*EL_CHG/(central_mass*ATOMIC_MASS_UNIT)) ! m/s
  E = E * beta + 0.5d0 * central_mass*ATOMIC_MASS_UNIT * v**2 / EL_CHG

  E = E*sqrt(atomic_weights(Z_ion)/central_mass) ! correct for atomic weight
end subroutine sample_fluid_particle_energy
end module mod_particle_sputtering
