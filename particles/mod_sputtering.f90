!> Very simple sputtering source on a point in the poloidal plane, with energies
!> sampled from the thompson distribution and a fixed source rate.
!>
!> Use [[mod_particle_sputtering]] instead of this one in almost all cases.
module mod_sputtering
  use mod_event
  use mod_sampling
  use mod_particle_types
  implicit none
  private
  public :: simple_sputter


  !> A line source, toroidally symmetric
  !> This is a hacky way to make a source. Better would be to make a new category
  !> which is linked to a specific particle group. That would work better also when
  !> using sources/sinks in different ways. (Also we should generalize the source/sink concept)
  !> For now it will work though, in all groups. Just assume that the sputtered atom
  !> belongs to the group we are looking at
  type, extends(action) :: simple_sputter
    real*8 :: x(2) !< R, Z
    real*8 :: st(2) !< s, t
    integer :: i_elm

    real*8 :: normal(3) !< Surface normal vector

    real*8 :: fill_fraction = 0.25d0 !< How many of the remaining particles to initialize. This is related to the prompt redeposition
    !< fraction and the sputtering timescale.
    real*8 :: source = 1d15 !< How many physical atoms to sputter per second
    real*8 :: sim_source = 1d8 !< How many simulation atoms to sputter per second at most (target)

    type(thompson_dist) :: E_dist = thompson_dist(E_b = 8.7d0, n=2) !< produces energies in eV (value for W default)
    real*8 :: last_time !< When did we sputter last (i.e. how many atoms do we need to make)
  contains
    procedure :: do => do_sputter
  end type simple_sputter

contains

!> We perform the sputtering by searching first for the number of lost particles on all mpi procs
!> We then distribute the to-sputter particles proportionally among processes.
!> Every process creates max(min(nint(n_free*fill_fraction),n_free),0) particles and initializes
!> these with R, Z, s, t, phi = u * TWOPI, i_elm and the sampled velocity.
!> We ignore any variation due to roundoff here, hich should not be important.
subroutine do_sputter(this, sim, ev)
  use mod_particle_sim
  use mpi
  use mod_pcg32_rng
  use mod_random_seed
  use constants, only: TWOPI, EL_CHG, ATOMIC_MASS_UNIT
  !$ use omp_lib
  class(simple_sputter), intent(inout) :: this
  type(particle_sim), intent(inout) :: sim
  type(event), intent(inout), optional :: ev
  integer :: i, j, k, n_free, i_err, n_free_total, n_sputter, seed, ierr, i_rng, n_stream, n_target
  integer, allocatable, dimension(:) :: i_free
  logical, allocatable, dimension(:) :: is_free
  real*8 :: u(4), free_frac

  type(pcg32_rng) :: rng

  do i=1,size(sim%groups,1)
    if (sim%my_id .eq. 0) seed = random_seed()
    call MPI_Bcast(seed, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    allocate(is_free(size(sim%groups(i)%particles,1)))
#ifdef __GFORTRAN__
    !$omp parallel default(shared) & 
#else
    !$omp parallel default(none) &
    !$omp shared(sim, this, n_free, i_free, is_free, n_free_total, i, seed, free_frac, n_sputter, n_target) &
#endif
  !$omp private(i_rng, n_stream, j, k, rng, u, ierr)
    i_rng = 1
    n_stream = 1
    !$ n_stream = omp_get_max_threads()
    !$ i_rng = omp_get_thread_num()+1
    call rng%initialize(4, seed, n_stream, i_rng)

    ! We need a loop here due to a gfortran bug with arrays of derived types
    !$omp do
    do j=1,size(sim%groups(i)%particles,1)
      is_free(j) = sim%groups(i)%particles(j)%i_elm .le. 0
    end do
    !$omp end do
    !$omp barrier
    !$omp single
    n_free = count(is_free)
    allocate(i_free(n_free))
    k = 1
    do j=1,size(sim%groups(i)%particles,1)
      if (is_free(j)) then
        i_free(k) = j
        k = k+1
      end if
    end do
    call MPI_AllReduce(n_free, n_free_total, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    n_free = size(i_free,1) ! because MPI_AllReduce might zero the value
    free_frac = real(n_free,8)/real(size(sim%groups(i)%particles,1),8)
    n_target = ceiling(min(real(n_free,4)*this%fill_fraction, this%sim_source*(sim%time - this%last_time)))
    n_sputter = max(min(n_target,n_free),0)
    if (sim%my_id .eq. 0 .and. (free_frac .le. 0.2)) write(*,*) 'group', i, 'free %', free_frac*100.d0, n_sputter
    !$omp end single

    ! Number of particles to sample is n_free*fill_fraction
    !$omp do
    do j=1,n_sputter
      k = i_free(j)
      call rng%next(u)
      sim%groups(i)%particles(k)%i_elm = this%i_elm
      sim%groups(i)%particles(k)%x = [this%x(1), this%x(2), TWOPI*u(1)]
      sim%groups(i)%particles(k)%st = this%st
      ! The particle weight is given by the number of particles sputtered from
      ! this source in this time divided by the number of
      ! simulation particles sputtered (all cpus). The number of simulation particles is given by
      ! n_free_all*this%fill_fraction
      sim%groups(i)%particles(k)%weight = real(this%source * (sim%time - this%last_time),4) / real(n_sputter,4)

      ! Implementation only for kinetic leapfrog particles currently
      select type (p => sim%groups(i)%particles(k))
      type is (particle_kinetic_leapfrog)
        ! sample_dist produces eV, sample_cosine produces the orientation vector
        ! scale to m/s
        p%v = sqrt(2.d0*sample_dist(this%E_dist, u(2))*EL_CHG/(sim%groups(i)%mass * ATOMIC_MASS_UNIT)) &
            * sample_cosine(u(3:4), this%normal)
        ! Since it is a neutral the half-step for boris method does not matter at all
        p%q = 0_1
      end select
    end do
    !$omp end do
    !$omp end parallel
    deallocate(i_free, is_free)
  end do

  ! Make a marker for the next time to sputter
  this%last_time = sim%time
end subroutine do_sputter
end module mod_sputtering
