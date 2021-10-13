!> Module for puffing gas into the plasma
!> This module will create new particles at the locations where gas valves will be.

module mod_particle_puffing
  use mod_edge_elements
  use mod_io_actions, only: io_action
  use mod_sampling
  use mod_particle_types
  use constants, only: TWOPI, K_BOLTZ, ATOMIC_MASS_UNIT
  use mod_rng, only: type_rng, setup_shared_rngs
  use mod_boundary, only: wall_normal_vector
  use mod_atomic_elements !mod_elements !< Chemical elements
  use mod_particle_sim
  use mod_event
  use mod_find_rz_nearby, only: find_rz_nearby

  
  

  implicit none

  private
  public  :: particle_puffing

  ! Extend type
  type, extends(io_action) :: particle_puffing
   
    class(type_rng), dimension(:), allocatable :: rng  !< one RNG per openmp thread
   
    ! number of simulation particles/s to puff across all processes
    integer :: n_puff = -1 
    ! Average fueling rate: 9.7d22; max fueling rate 18d22
    real*8  :: fueling_rate = -1.d0
    real*8  :: R = -1.d0, Z = -1.d0, phi = -1.d0
    real*8  :: valve_r = -1.d0  !< radius of gas valve
    real*8  :: last_time = 0.d0 !< When did we puff last 
    real*8 :: last_diag_time = 0.d0 !< Last time of output of diagnostics

  contains
    procedure :: do => do_particle_puffing
  end type particle_puffing

  interface particle_puffing
    module procedure new_particle_puffing
  end interface particle_puffing
contains

function new_particle_puffing(n_puff, fueling_rate, valve_r, R, Z, phi, rng, seed) result(new)
  use mod_pcg32_rng,   only: pcg32_rng
  use mod_random_seed, only: random_seed
  
  type(particle_puffing)    :: new

  integer, intent(in)           :: n_puff
  real*8, intent(in)            :: fueling_rate
  real*8, intent(in)            :: valve_r
  real*8, intent(in)            :: R, Z
  real*8, intent(in), optional  :: phi ! If no phi is given axisymmetric puffing will be excecuted.

  class(type_rng), intent(in), optional :: rng !< random number generator to use (deafult PCG32)
  integer, intent(in), optional         :: seed !< Seed for the RNG (default random_seed() on my_id + bcast)
  integer                               :: my_seed
  
  new%n_puff       = n_puff
  new%fueling_rate = fueling_rate
  new%R            = R
  new%Z            = Z
  new%valve_r      = valve_r
  if (present(phi))  new%phi = phi
  


  !> allocate random seed for sampling
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

  ! Actually puff gass
subroutine do_particle_puffing(this,sim, ev)

  class(particle_puffing) , intent(inout) :: this
  type(particle_sim), intent(inout)       :: sim
  type(event), intent(inout), optional    :: ev !<STIJN> is this nececary?

  integer :: i_scalar, n_free, j, k, n_group, i_elm, i_elm_new, ifail, i_p, to_puff
  logical, allocatable, dimension(:) :: is_free
  integer, allocatable, dimension(:) :: i_free
  real*8  :: delta_t, c, R, Z, phi, s, t
  real*8  :: R_new, Z_new, s_new, t_new, r_valve, theta
  real*8  :: vector_normal(3), u(5)


  if (sim%my_id .eq. 0) write(*,*) "Started puffing!"
  
  if (this%last_time .eq. 0.d0) then
    this%last_time = sim%time
    this%last_diag_time = sim%time
    if (abs(sim%time) .le. 1d-10) then
      this%last_time = 1d-10 ! so we actually do most of the sputtering should we start at 0
    end if
    return
  end if
  delta_t = sim%time - this%last_time
  this%last_time = sim%time
  
  if (this%n_puff .le. -1.d-6) then ! 0.d0
    if (sim%my_id .eq. 0) write(*,*) 'No puf quota set, exiting. --- n_puff == 0 this will now stop the program'
    stop
  end if

  allocate(is_free(size(sim%groups(1)%particles,1))) 

  ! might be replaced with omp workshare, or just the array expression.
  ! there is an issue with derived type arrays in gfortran though, and this works
  !!$omp parallel do default(none) shared(sim, this, n_free, i_free, is_free) &
  !!$omp private(j)
  do j=1,size(sim%groups(1)%particles,1)
    is_free(j) = sim%groups(1)%particles(j)%i_elm .le. 0 
  end do
  !!$omp end parallel do
  !!$omp barrier
  n_free = count(is_free)
  allocate(i_free(n_free))
  k = 1
  do j=1,size(is_free,1)
    if (is_free(j)) then
      i_free(k) = j
      k = k+1
      !if (sim%my_id .eq. 0) write(*,*) "Adding to the list number: ", j
    end if
  end do

  n_group = 1   ! Puffing Hydrogen (or actually the element at groups(1)) only
  ! Assuming the incoming gas at T=300K and a diatomic gas
  c = sqrt((7.d0/5.d0)*300.d0*K_BOLTZ/(2.d0*sim%groups(n_group)%mass*ATOMIC_MASS_UNIT))

  call find_RZ(sim%fields%node_list, sim%fields%element_list, this%R, this%Z, R, Z, &
               i_elm, s, t ,ifail)
  if (ifail .ne. 0) then
    if (sim%my_id .eq. 0) write(*,*) "Warning: The valve location for puffing could not be found, maybe it was placed outside of the grid?"
    stop
  end if

  vector_normal = wall_normal_vector(sim%fields%node_list, sim%fields%element_list, &
          i_elm, s, t)

  if (this%n_puff .ge. n_free) then
    write(*,*) "Warning could not puff the requested amount."
    to_puff = n_free
  else
    to_puff = this%n_puff
  end if

  select type (pa => sim%groups(1)%particles)
  type is (particle_kinetic_leapfrog)
    do j = 1, to_puff
      i_p = i_free(j)
      do 
        call this%rng(1)%next(u)
        r_valve = this%valve_r*sample_piecewise_linear(2, [0.d0, 1.d0], [1.d0, 0.d0], u(1))
        theta = TWOPI * u(2)
        R_new = R + r_valve * cos(theta)
        Z_new = Z + r_valve * sin(theta)
        call find_RZ_nearby(sim%fields%node_list, sim%fields%element_list, R, Z, s, t, i_elm, &
        R_new, Z_new, s_new, t_new, i_elm_new, ifail)
        !call find_RZ(sim%fields%node_list, sim%fields%element_list, R_new, Z_new, R, Z, &
        !       i_elm_new, s_new, t_new, ifail)
        if (ifail .ge. 0) exit
      end do
      R     = R_new
      Z     = Z_new
      s     = s_new
      t     = t_new
      i_elm = i_elm_new
      if (this%phi .lt. 0.d0) then
        pa(i_p)%x(3) = TWOPI*u(3)
      else 
        pa(i_p)%x(3) = this%phi 
      end if
      pa(i_p)%x(1:2)  = [R, Z]
      pa(i_p)%st(1:2) = [s, t]
      pa(i_p)%i_elm   = i_elm
      pa(i_p)%weight  = real(1.d0/this%n_puff) * delta_t * this%fueling_rate
      pa(i_p)%v       = c * sample_cosine(u(4:5), vector_normal)   ! <STIJN> Maybe this needs to be isotropic, or 1+cos like?
      pa(i_p)%q       = 0_1
      if (sim%groups(1)%particles(i_p)%weight  .le. 1.d-2) then ! if the weight is too low. 
        sim%groups(1)%particles(i_p)%i_elm = 0
        cycle       
      end if
    end do
  class default
    write(*,*) 'Particle type not impletement for gas fueling.'
    stop
  end select

end subroutine do_particle_puffing
end module