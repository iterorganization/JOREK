!> Module containing routines to project a function of particles onto the JOREK
!> finite elements.
!>
!> The general routine allows the user to specify a transformation object
!> to calculate the quantity to be projected.
!> This quantity is multiplied by the particle weight behind the scenes.
!> Included transformation types
!>   - 1 (proj_one)
!>
!> These functions must subclass proj_transform in this module.
!> Other modules defining transformations:
!>   - tools/mod_radiation.f90
!>
!>
!> Projected results can be written to vtk or hdf5 by specifying the appropriate
!> options at creation time of project_particles_base.
!>
!> Projections are written to a node_list. The number of projections is thus limited
!> by the number of JOREK variables. If you need more, use multiple instances.
module mod_project_particles

use mod_io_actions
use data_structure
use mod_particle_sim
use mod_particle_types
use mod_fields
use mod_vtk
use constants, only : el_chg, atomic_mass_unit

implicit none
include 'dmumps_struc.h'        ! MUMPS include files defining its datastructure
private
public projection
public new_projection !< The constructor function is also public since it provides better error handling than the real constructor
public proj_f
public proj_f_interface, proj_one, proj_q, proj_vR, proj_vZ, proj_vPhi, proj_Ekin, proj_Ekin_keV, proj_jR, proj_jZ, proj_jPhi
public write_particle_distribution_to_vtk, write_particle_distribution_to_h5 !< public for testing reasons, please don't use directly
public prepare_mumps_par, prepare_mumps_par_n0, sample_rhs !< public for testing reasons
public DMUMPS_STRUC

interface
  function proj_f_interface(sim, group, particle)
    import particle_sim, particle_base
    type(particle_sim), intent(in) :: sim
    integer, intent(in) :: group
    class(particle_base), intent(in) :: particle !< Input particle
    real*8 :: proj_f_interface !< Value to be projected
  end function proj_f_interface
end interface

!> A wrapper type to contain a projection function and a group to apply it to
type :: proj_f
  procedure(proj_f_interface), nopass, pointer :: f
  integer :: group !< group number to apply this function to
end type proj_f
interface proj_f
  module procedure new_proj_f !< The constructor for this type
end interface proj_f

!> Action to project all particle distributions and save them to vtk.
!> You must use the (new_)projection() constructor to set this up.
!> At construction time the matrix is solved.
type, extends(io_action) :: projection
  type(type_node_list),    pointer :: node_list !< node lists to save particle projections in
  type(type_element_list), pointer :: element_list

  real*8 :: filter                      !< Smoothing factor used for this projection (Laplacian, poloidal plane)
  real*8 :: filter_hyper                !< hyper-smoothing factor used for this projection (double Laplacian, poloidal plane)
  real*8 :: filter_parallel             !< smoothing factor used for this projection (parallel direction)

  real*8 :: filter_n0                   !< Smoothing factor used for this projection for n=0 
  real*8 :: filter_hyper_n0             !< hyper-smoothing factor used for this projection for n=0
  real*8 :: filter_parallel_n0          !< smoothing factor used for this projection (parallel direction) for n=0

  real*8 :: scaling_integral_weights    !< multiplication factor to subtract integral_weights from rhs(n=0)

  logical,public :: do_zonal = .false.  !< solve zonal flow system for n=0 (instead of usual projection)

  !> Output storage (optional)
  type(vtk_grid), allocatable, private :: vtk_grid !< if allocated output to vtk
  logical, public :: to_h5 = .false.    !< Output to hdf5 file

  !> Right-hand side
  type(proj_f), dimension(:),   allocatable :: f   !< List of projection transformations to use (n_proj)
  real*8, dimension(:,:,:,:,:), allocatable :: rhs !< dim (n_order+1,n_vertex_max,n_elements,n_tor,n_proj2)
  !< right-hand side for accumulation during sampling
  !< assumed to be filled by the user. Will be MPI_Reduced (+) before projecting
  !< n_proj + n_proj2 should be less than n_var (extra input will be ignored)
  !< After projection this will be zeroed but not deallocated
  real*8 :: rhs_gather_time = 0.d0 !< Time that the rhs has been integrated over (used for normalisation)
  !< note that this does not really work very well for multiple groups with
  !< different timesteps

  logical :: calc_integrals = .true.    !< Calculate and print integrals of all projected quantities over the entire volume
  !< (for n=0 only)
  !> Note that the integral is just a weighted sum over the node values.
  !> Precalculate the weights during the projection matrix assembly step and store
  !> them here.
  real*8, dimension(:), allocatable     :: integral_weights !< Weights per basis function towards full integral
  !< indexing is as in rhs, and the weight is the volume of the basis function

  real*8 :: area, volume

  real*8, dimension(:,:,:,:,:), allocatable :: rhs_f !< dim (n_order+1,n_vertex_max,n_elements,n_tor,n_proj) storage
  !< location for proj_f output

  integer :: mpi_comm_world    ! mpi communicator of the whole world
  integer :: mpi_comm_n        ! mpi communicator of each toroidal harmonic
  integer :: mpi_comm_master   ! mpi communicator of the group of masters (of each harmonic)
  integer :: my_id             ! mpi id within the world 
  integer :: my_id_n           ! mpi id with comm_n
  integer :: mpi_group_world
  integer :: mpi_group_master
  integer :: n_cpu, m_cpu      ! n_cpu : the total number of cores, m_cpu the number of cores per harmonic
  integer :: n_tor_local       ! 1 or 2 : (1) or (cos,sin)
  integer :: i_tor_local       ! the starting index in the array of toroidal hamonics (as in HZ)
  integer :: n_dof             ! the number of unknowns for (n=0)

  !> Internal variables
  type (DMUMPS_STRUC) :: mumps_par !< matrix is factored by mumps and stored here
contains
  procedure :: do => project
  procedure :: close_mumps => close_mumps
end type projection
interface projection
  module procedure new_projection
end interface projection

contains

!*****************************************************************
!* Example projection functions for generating a right-hand-side *
!*****************************************************************

!> Project the particle density by using transformation function 1
pure function proj_one(sim, group, particle)
  type(particle_sim), intent(in) :: sim
  integer, intent(in) :: group
  class(particle_base), intent(in) :: particle
  real*8 :: proj_one
  proj_one = 1.d0
end function proj_one

pure function proj_vR(sim, group, particle)
  use mod_particle_types, only: particle_kinetic_leapfrog
  type(particle_sim), intent(in) :: sim
  integer, intent(in) :: group
  class(particle_base), intent(in) :: particle
  real*8 :: proj_vR
  select type (p => particle)
  type is (particle_kinetic_leapfrog)
    proj_vR = p%v(1)
  class default
    proj_vR = 0.d0
  end select
end function proj_vR

pure function proj_vZ(sim, group, particle)
  use mod_particle_types, only: particle_kinetic_leapfrog
  type(particle_sim), intent(in) :: sim
  integer, intent(in) :: group
  class(particle_base), intent(in) :: particle
  real*8 :: proj_vZ
  select type (p => particle)
  type is (particle_kinetic_leapfrog)
    proj_vZ = p%v(2)
  class default
    proj_vZ = 0.d0
  end select
end function proj_vZ

pure function proj_vPhi(sim, group, particle)
  use mod_particle_types, only: particle_kinetic_leapfrog
  type(particle_sim), intent(in) :: sim
  integer, intent(in) :: group
  class(particle_base), intent(in) :: particle
  real*8 :: proj_vPhi
  select type (p => particle)
  type is (particle_kinetic_leapfrog)
    proj_vPhi = p%v(3)
  class default
    proj_vPhi = 0.d0
  end select
end function proj_vPhi

!< Energy in joules
pure function proj_Ekin(sim, group, particle) 
  use mod_particle_types, only: particle_kinetic_leapfrog
  type(particle_sim), intent(in) :: sim
  integer, intent(in) :: group
  class(particle_base), intent(in) :: particle
  real*8 :: proj_Ekin
  select type (p => particle)
  type is (particle_kinetic_leapfrog)
    proj_Ekin = sim%groups(group)%mass * atomic_mass_unit * dot_product(p%v,p%v)/2.d0
  class default
    proj_Ekin = 0.d0
  end select
end function proj_Ekin

pure function proj_Ekin_keV(sim, group, particle) 
  use mod_particle_types, only: particle_kinetic_leapfrog
  type(particle_sim), intent(in) :: sim
  integer, intent(in) :: group
  class(particle_base), intent(in) :: particle
  real*8 :: proj_Ekin_keV
  select type (p => particle)
  type is (particle_kinetic_leapfrog)
    proj_Ekin_keV = sim%groups(group)%mass * atomic_mass_unit * dot_product(p%v,p%v)/2.d0 / el_chg /1.d3
  class default
    proj_Ekin_keV = 0.d0
  end select
end function proj_Ekin_keV

!> Project the particle charge by using transformation function q
!> (only valid for particles of type kinetic(_leapfrog) or gc
!>
!> TODO: normalize projection with density. This function calculates
!> the integral of q instead of the mean value.
pure function proj_q(sim, group, particle)
  type(particle_sim), intent(in) :: sim
  integer, intent(in) :: group
  class(particle_base), intent(in) :: particle
  real*8 :: proj_q
  select type (p => particle)
  type is (particle_kinetic)
    proj_q = real(p%q, 8) * el_chg
  type is (particle_kinetic_leapfrog)
    proj_q = real(p%q, 8) * el_chg
  type is (particle_gc)
    proj_q = real(p%q, 8) * el_chg
  class default
    proj_q = 0.d0
  end select
end function proj_q

pure function proj_jR(sim, group, particle)
  type(particle_sim), intent(in) :: sim
  integer, intent(in) :: group
  class(particle_base), intent(in) :: particle
  real*8 :: proj_jR
  select type (p => particle)
  type is (particle_kinetic)
    proj_jR = real(p%q, 8) * el_chg * p%v(1)
  type is (particle_kinetic_leapfrog)
    proj_jR = real(p%q, 8) * el_chg * p%v(1)
  class default
    proj_jR = 0.d0
  end select
end function proj_jR

pure function proj_jZ(sim, group, particle)
  type(particle_sim), intent(in) :: sim
  integer, intent(in) :: group
  class(particle_base), intent(in) :: particle
  real*8 :: proj_jZ
  select type (p => particle)
  type is (particle_kinetic)
    proj_jZ = real(p%q, 8) * el_chg * p%v(2)
  type is (particle_kinetic_leapfrog)
    proj_jZ = real(p%q, 8) * el_chg * p%v(2)
  class default
    proj_jZ = 0.d0
  end select
end function proj_jZ

pure function proj_jPhi(sim, group, particle)
  type(particle_sim), intent(in) :: sim
  integer, intent(in) :: group
  class(particle_base), intent(in) :: particle
  real*8 :: proj_jPhi
  select type (p => particle)
  type is (particle_kinetic)
    proj_jPhi = real(p%q, 8) * el_chg * p%v(3)
  type is (particle_kinetic_leapfrog)
    proj_jPhi = real(p%q, 8) * el_chg * p%v(3)
  class default
    proj_jPhi = 0.d0
  end select
end function proj_jPhi

!> Constructor for projection function
function new_proj_f(f, group)
  type(proj_f) :: new_proj_f
  procedure(proj_f_interface), pointer, intent(in) :: f
  integer, intent(in) :: group
  new_proj_f%f => f
  new_proj_f%group = group
end function new_proj_f
  

!> Constructor for project_particles
!> Be sure to use keyword arguments when initializing, to avoid confusion
function new_projection(node_list, element_list,                                                    &
                        filter,    filter_hyper,    filter_parallel,                                &
                        filter_n0, filter_hyper_n0, filter_parallel_n0,                             &
                        f, do_zonal, to_h5, to_vtk,                                                 &
                        nsub, filename, basename, decimal_digits, fractional_digits, calc_integrals &
                        ) result(new)
  use mpi_mod
  !use mod_parameters, only n_node_max
  type(projection) :: new
  type(type_node_list), intent(in)       :: node_list
  type(type_element_list), intent(in), target    :: element_list
  real*8, intent(in), optional           :: filter,    filter_hyper,    filter_parallel    !< normal, hyper and parallel smoothing
  real*8, intent(in), optional           :: filter_n0, filter_hyper_n0, filter_parallel_n0 !< for n=0
  type(proj_f), intent(in), dimension(:), optional :: f !< Type with function to map over particles before projection
  logical, intent(in), optional          :: do_zonal    !< solve zonal flow  system for n=0 instead of projection (false if omitted)
  logical, intent(in), optional          :: to_h5 !< Write HDF5 output after projecting (false if omitted)
  logical, intent(in), optional          :: to_vtk !< Write vtk output after projecting (false if omitted)
  integer, intent(in), optional          :: nsub !< number of subdivisions of the finite elements
  character(len=*), intent(in), optional :: filename
  character(len=*), intent(in), optional :: basename
  integer, intent(in), optional          :: decimal_digits
  integer, intent(in), optional          :: fractional_digits
  logical, intent(in), optional          :: calc_integrals !< After projecting, calculate and print the integral of each projected quantity

  integer              :: ierr, my_nsub, inode, n_masters, i
  integer, allocatable :: i_tor(:)
  integer              :: i_rank(n_tor), my_id_tmp

  call MPI_Comm_dup(MPI_COMM_WORLD, new%mpi_comm_world, ierr)

  call MPI_COMM_RANK(new%mpi_comm_world, new%my_id, ierr)
  call MPI_COMM_SIZE(new%mpi_comm_world, new%n_cpu, ierr)

  n_masters = (n_tor+1)/2
  if (MOD(new%n_cpu, n_masters) == 0) then
    new%m_cpu = new%n_cpu / (n_masters)
  else
    new%m_cpu = (new%n_cpu - MOD(new%n_cpu, n_masters))/n_masters +1
  end if

  if (allocated(i_tor)) call tr_deallocate(i_tor,"i_tor",CAT_UNKNOWN)
  call tr_allocate(i_tor,1,new%n_cpu,"i_tor",CAT_UNKNOWN)

  do i=1,new%n_cpu
    i_tor(i) = ((i-1) - MOD(i-1, new%m_cpu))/ new%m_cpu  + 1
  enddo

  call MPI_COMM_SPLIT(new%mpi_comm_world,i_tor(new%my_id+1),new%my_id,new%mpi_comm_n,ierr)
  
  do i=1,n_masters
    i_rank(i) = (i-1) * new%m_cpu
  enddo

  call MPI_COMM_GROUP(new%mpi_comm_world,new%mpi_group_world,ierr)
  call MPI_GROUP_INCL(new%mpi_group_world,n_masters,i_rank,new%mpi_group_master,ierr)
  call MPI_COMM_CREATE(new%mpi_comm_world,new%mpi_group_master,new%mpi_comm_master,ierr)
  call MPI_COMM_RANK(new%mpi_comm_n, new%my_id_n, ierr) ! id of this cpu in local comm
  
  if (i_tor(new%my_id+1) .eq. 1) then
    new%n_tor_local = 1
    new%i_tor_local = 1
  else
    new%n_tor_local = 2
    new%i_tor_local = 2*i_tor(new%my_id+1) - 2       ! i_tor_local is the (starting) index in HZ
  endif

  allocate(new%node_list,    source=node_list)
  new%node_list = node_list

  do inode = 1, n_nodes_max
    new%node_list%node(inode)%values = 0.d0
    new%node_list%node(inode)%deltas = 0.d0
  end do

  allocate(new%element_list, source=element_list)
  new%element_list = element_list

  new%filter             = 0.d0
  new%filter_hyper       = 0.d0
  new%filter_parallel    = 0.d0
  new%filter_n0          = 0.d0
  new%filter_hyper_n0    = 0.d0
  new%filter_parallel_n0 = 0.d0

  if (present(filter))             new%filter             = filter
  if (present(filter_hyper))       new%filter_hyper       = filter_hyper
  if (present(filter_parallel))    new%filter_parallel    = filter_parallel
  if (present(filter_n0))          new%filter_n0          = filter_n0
  if (present(filter_hyper_n0))    new%filter_hyper_n0    = filter_hyper_n0
  if (present(filter_parallel_n0)) new%filter_parallel_n0 = filter_parallel_n0

  new%do_zonal = .false.
  if (present(do_zonal)) new%do_zonal = do_zonal

  if (present(to_vtk)) then
    if (to_vtk) then
      my_nsub = 4
      if (present(nsub)) my_nsub = nsub
      ! Precalculate the node positions in the vtk file and the connectivity
      if (new%my_id .eq. 0) then
        new%vtk_grid = vtk_grid(new%node_list, new%element_list, my_nsub)
      endif
      ! We don't set the extension here since this is dynamically set in the do
      ! action (to support both vtk and h5 output)
    end if
  end if
  if (present(to_h5)) new%to_h5 = to_h5

  new%basename = "proj"
  if (present(filename)) new%filename = filename
  if (present(basename)) new%basename = basename
  if (present(decimal_digits)) new%decimal_digits = decimal_digits
  if (present(fractional_digits)) new%fractional_digits = fractional_digits
  if (present(calc_integrals)) new%calc_integrals = calc_integrals
  new%name = "Project"
  new%log = .true.

  if (new%n_tor_local .eq. 1) then

    call prepare_mumps_par_n0(node_list, element_list, new%n_tor_local, new%i_tor_local,            & 
                              new%mpi_comm_world, new%mpi_comm_n, new%mpi_comm_master,              &
                              new%mumps_par, new%area, new%volume,                                  &
                              new%filter_n0, new%filter_hyper_n0, new%filter_parallel_n0,           &
                              integral_weights=new%integral_weights, do_zonal=new%do_zonal)  

    new%n_dof = new%mumps_par%n / 2
  
  else

    call prepare_mumps_par(node_list, element_list, new%n_tor_local, new%i_tor_local,            &
                           new%mpi_comm_world, new%mpi_comm_n, new%mpi_comm_master,              &
                           new%mumps_par, new%filter, new%filter_hyper, new%filter_parallel)

    new%n_dof = new%mumps_par%n / new%n_tor_local

  end if


  if (present(f)) then
    new%f = f
  endif

end function new_projection

subroutine close_mumps(this)
  class(projection), intent(inout) :: this
  this%mumps_par%JOB = -2
  call DMUMPS(this%mumps_par)
end subroutine close_mumps


subroutine project(this, sim, ev)
  use mod_event
  class(projection), intent(inout)     :: this
  type(particle_sim), intent(inout)    :: sim
  type(event), intent(inout), optional :: ev

  ! Sample projection functions for each particle
  call sample_rhs(this, sim)

  ! Project all right-hand sides
  call project_only(this, sim)

  ! Save output if requested
  if (this%to_h5) call save_to_h5(this, sim)
  if (allocated(this%vtk_grid)) then
    call save_to_vtk(this, sim)
  endif

  ! Clean up storage
  if (allocated(this%rhs)) this%rhs = 0.d0
  if (allocated(this%rhs_f)) deallocate(this%rhs_f)
end subroutine project


!> Gather all of the rhs-es into a single matrix and feed it to mumps, and then
!> broadcast the result
subroutine project_only(this, sim)
  use mpi_mod
  use mod_event
  use, intrinsic :: ieee_exceptions
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
  !$ use omp_lib
  class(projection), intent(inout) :: this
  type(particle_sim), intent(inout)    :: sim
  integer :: my_id, ierr
  integer :: in, i_elm, i, j, k, i_var, i_start, i_tor
  integer :: index_large_i, inode, index
  real*8  :: t0, t1, ostart, oend, mmm(3), mmm2(3)
  real*8,  allocatable :: my_rhs(:,:), y_tmp(:)
  integer, allocatable :: recv_counts(:), recv_disp(:)
  integer :: n_rhs, n_rhs_f, i_rhs, n_tor_local, i_tor_local, n_loc_n
  integer :: in_local, in_global, index_n, id_master_in_world, offset
  logical :: halt(size(IEEE_USUAL,1)), found_nan

  call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)
  call cpu_time(t0)
  !$ ostart = omp_get_wtime()

  ! Safety checks
  if (.not. allocated(sim%groups)) return
  
  if (.not. allocated(this%rhs)) then
    n_rhs = 0
  else
    n_rhs = size(this%rhs,5)
  end if

  if (.not. allocated(this%rhs_f)) then
    n_rhs_f = 0
  else 
    n_rhs_f = size(this%rhs_f,5)
  end if

  if (n_rhs .gt. n_var .and. this%my_id_n .eq. 0) then
    write(*,*) 'WARNING: too many rhs-es supplied, skipping after ', n_var
    n_rhs = n_var
    ! n_rhs_f should already be 0 in this case, as done by sample_rhs
  end if

  n_tor_local = this%n_tor_local
  i_tor_local = this%i_tor_local

  this%mumps_par%nrhs = (n_rhs + n_rhs_f)
  this%mumps_par%lrhs = this%mumps_par%n

  if (2*this%n_dof .ne. this%mumps_par%n) then
    write(*,*) 'FATAL : 2*this%n_dof .ne. this%mumps_par%n'
    write(*,*) n_tor_local*this%n_dof,  this%mumps_par%n
  endif

  if (this%my_id_n .eq. 0) then
    ! For some reason gfortran throws an error with the allocated statement here. Instead just reallocate on every call
    !if (allocated(this%mumps_par%rhs) .and. size(this%mumps_par%rhs)/this%mumps_par%n .ne. n_rhs + n_rhs_f) deallocate(this%mumps_par%rhs)
    !if (.not. allocated(this%mumps_par%rhs)) allocate(this%mumps_par%rhs(this%mumps_par%n*(n_rhs+n_rhs_f)))
    allocate(this%mumps_par%rhs(this%mumps_par%n*this%mumps_par%nrhs))
  else
    allocate(this%mumps_par%rhs(0)) ! dummy allocation for MPI
  end if

  this%mumps_par%rhs = 0.d0

  allocate(my_rhs(this%mumps_par%n * this%mumps_par%nrhs, (n_tor+1)/2))
  
  my_rhs = 0.d0

  do i_rhs=1,n_rhs

    i_start =  this%mumps_par%n * (i_rhs-1)
        
    do i_elm=1,this%element_list%n_elements

      do i=1,n_vertex_max

        inode = this%element_list%element(i_elm)%vertex(i)

        do j=1,n_order+1

          index_large_i = 2 * (this%node_list%node(inode)%index(j)-1) + 1 + i_start ! base index in the main matrix + rhs index

          my_rhs(index_large_i,1)   = my_rhs(index_large_i,1)   + this%rhs(j, i, i_elm, 1, i_rhs) !/ this%rhs_gather_time
          my_rhs(index_large_i+1,1) = 0.d0

          do in=2, n_tor, 2

            index_n = in / 2 + 1

            my_rhs(index_large_i,  index_n) = my_rhs(index_large_i,  index_n) + this%rhs(j, i, i_elm, in,   i_rhs) !/ this%rhs_gather_time
            my_rhs(index_large_i+1,index_n) = my_rhs(index_large_i+1,index_n) + this%rhs(j, i, i_elm, in+1, i_rhs) !/ this%rhs_gather_time

          enddo
        enddo

      enddo
    enddo
  enddo
  
  this%rhs_gather_time = 0.d0

  ! Reset the RHS for the next step
  if (allocated(this%rhs)) this%rhs = 0.d0

  do i_rhs=1,n_rhs_f
    ! Fill projection function part

    i_start =  this%mumps_par%n * (n_rhs + i_rhs - 1)  
    
    do i_elm=1,this%element_list%n_elements
        
      do i=1,n_vertex_max
      
        inode = this%element_list%element(i_elm)%vertex(i)
          
        do j=1,n_order+1

          index_large_i = 2*(this%node_list%node(inode)%index(j)-1) + 1 + i_start ! base index in the main matrix + rhs index

          my_rhs(index_large_i,1)   = my_rhs(index_large_i,1)   + this%rhs_f(j, i, i_elm, 1, i_rhs)
          my_rhs(index_large_i+1,1) = 0.

          do in=2, n_tor, 2

            index_n = in / 2 + 1
        
            my_rhs(index_large_i,  index_n) = my_rhs(index_large_i,  index_n) + this%rhs_f(j, i, i_elm, in,   i_rhs)
            my_rhs(index_large_i+1,index_n) = my_rhs(index_large_i+1,index_n) + this%rhs_f(j, i, i_elm, in+1, i_rhs)
            
          enddo
    
        enddo
      enddo

    enddo
  enddo

  ! Gather the RHS's to the root process
  ! cannot do it directly into mumps_par%rhs because this is not allocated in every process
  ! results are gathered on the my_id_n=0 nodes

  call MPI_Reduce(my_rhs(:,1),this%mumps_par%rhs,this%mumps_par%n*this%mumps_par%nrhs, MPI_REAL8, MPI_SUM, 0, this%mpi_comm_world, ierr)

  do in=2, n_tor, 2
    id_master_in_world = in/2 * this%m_cpu
    index_n = in/2 + 1
    call MPI_Reduce(my_rhs(:,index_n),this%mumps_par%rhs,this%mumps_par%n*this%mumps_par%nrhs, MPI_REAL8, MPI_SUM, id_master_in_world, this%mpi_comm_world, ierr)
  enddo

  if ((this%my_id .eq. 0) .and. (this%do_zonal)) then   
    write(*,'(A,3e14.6)') 'check project :',maxval(this%mumps_par%rhs(:)), maxval(this%scaling_integral_weights*this%integral_weights), &
                                            maxval(this%mumps_par%rhs(:) - this%scaling_integral_weights * this%integral_weights)

    this%mumps_par%rhs(1:this%mumps_par%n) = this%mumps_par%rhs(1:this%mumps_par%n) - this%scaling_integral_weights * this%integral_weights(:)
  endif

  ! Compute the solution of Ax=B (B = RHSes)
  this%mumps_par%JOB = 3
  this%mumps_par%icntl(21) = 0 ! solution is available only on host
  this%mumps_par%icntl(4)  = 0 !1 ! print only errors == 1

  ! Disable floating point exceptions in MUMPS
  ! some of the MKL routines make these exceptions on some vectorized
  ! calculations but then don't use the result for a speed increase.
  ! To allow running our code with -fpe0 we need to temporarily disable the
  ! checks, otherwise it'll crash here.
  call ieee_get_halting_mode(IEEE_USUAL, halt)
  call ieee_set_halting_mode(IEEE_USUAL, [.false., .false., .false.])
  call DMUMPS(this%mumps_par)
  call ieee_set_halting_mode(IEEE_USUAL, halt)

  ! collect the solution of all the toroidal harmonics (ntor+1)/2 to my_id=0
  

  if (this%my_id_n .eq. 0) then

    n_loc_n = this%mumps_par%n * this%mumps_par%nrhs
  
    allocate(y_tmp(n_loc_n*(n_tor+1)))            ! allocate only on my_id=0???

    allocate(recv_counts(this%n_cpu/this%m_cpu))
    allocate(recv_disp(this%n_cpu/this%m_cpu))
    
    y_tmp = 0.d0

    recv_counts(1) = n_loc_n
    do i=2,(n_tor+1)/2
      recv_counts(i) = n_loc_n
    enddo

    recv_disp(1) = 0
    do i=2,(n_tor+1)/2
      recv_disp(i) = recv_disp(i-1) + recv_counts(i-1)
    enddo

    call mpi_gatherv(this%mumps_par%rhs, this%mumps_par%n * this%mumps_par%nrhs, MPI_DOUBLE_PRECISION, &
                     y_tmp, recv_counts, recv_disp, MPI_DOUBLE_PRECISION, 0, this%mpi_comm_master,ierr)

  endif

  call MPI_BARRIER(this%mpi_comm_world, ierr)

  ! Write the solution to the node_list
  if (this%my_id .eq. 0) then

    do i_var=1,min(n_rhs+n_rhs_f,n_var)
  
      found_nan = .false.
      
      do i=1,this%node_list%n_nodes

        do k=1,n_order+1
      
          index = this%node_list%node(i)%index(k)

          if (this%do_zonal) then
             this%node_list%node(i)%values(1,k,i_var) = y_tmp(2*(index-1) + 1 + 2*this%n_dof*(i_var-1)) + y_tmp(2*(index-1) + 2 + 2*this%n_dof*(i_var-1))
             this%node_list%node(i)%values(1,k,2)     = y_tmp(2*(index-1) + 2 + 2*this%n_dof*(i_var-1))
          else
             this%node_list%node(i)%values(1,k,i_var) = y_tmp(2*(index-1) + 1 + 2*this%n_dof*(i_var-1))
          endif

          offset = 2*this%n_dof * this%mumps_par%nrhs
          
          do i_tor=2,n_tor,2

            this%node_list%node(i)%values(i_tor,  k,i_var) = y_tmp(2*(index-1) + 1 + offset + 2*this%n_dof*(i_var-1) + (i_tor-2)*this%n_dof * this%mumps_par%nrhs)
            this%node_list%node(i)%values(i_tor+1,k,i_var) = y_tmp(2*(index-1) + 2 + offset + 2*this%n_dof*(i_var-1) + (i_tor-2)*this%n_dof * this%mumps_par%nrhs)
          
          end do
        
        enddo    ! order
        
        ! Check for NaNs in the projection
        if (any(ieee_is_nan(this%node_list%node(i)%values(:,:,i_var)))) then
          found_nan = .true.
        end if
      
      enddo      ! nodes

      if (found_nan) then
        write(*,*) "Found NaNs in projection number ", i_var
      end if

    enddo
  endif

  call MPI_BARRIER(this%mpi_comm_world, ierr)

  call broadcast_nodes(this%my_id, this%node_list)

#ifdef DEBUG
    call cpu_time(t1)
    !$ oend = omp_get_wtime()
    !$ mmm = mpi_minmeanmax(t1-t0)
    !$ mmm2 = mpi_minmeanmax(oend-ostart)
    if (this%my_id .eq. 0) then
      write(*,"(A,3g12.5)") "projection cpu time", mmm
      !$ write(*,"(A,3g12.5)") "projection wall time", mmm2
    end if
#endif

end subroutine project_only


!> Add samples to the right-hand side, stored in `this` by calling this%f(i)%f
!> for every particle and saving the contribution multiplied by each of the basis
!> functions (poloidal and toroidal) in `this%rhs_f`
subroutine sample_rhs(this, sim)
  use mpi_mod
  use mod_event
  use mod_interp, only: mode_moivre, interp_RZ
  use constants, only: PI
  use mod_basisfunctions
  !$ use omp_lib
  class(projection), intent(inout)  :: this
  type(particle_sim), intent(inout) :: sim
  integer :: n_sample !< number of groups to sample
  real*8 :: HP(n_tor)
  real*8 :: v, R_g, R_s, R_t, Z_g, Z_s, Z_t, xjac, x(3), HH(4,4), HH_s(4,4), HH_t(4,4)
  integer :: i_group, m, i, j, im, im_index, i_f
  integer :: index_ij, i_tor
  ! For openmp reduce
  real*8, allocatable :: my_rhs(:,:,:,:,:)


  ! Safety checks
  if (.not. allocated(sim%groups)) return
  
  if (.not. allocated(this%f)) allocate(this%f(0))

  if (allocated(this%rhs)) then
    n_sample = min(size(this%f),n_var-min(size(this%rhs,5),n_var))  ! because we have only n_var storage for now
  else
    n_sample = min(size(this%f),n_var)  ! because we have only n_var storage for now
  end if

  if (n_sample .lt. size(this%f) .and. sim%my_id .eq. 0) then
    write(*,*) 'WARNING: ignoring proj_f after ', n_sample, ' due to lack of output space'
  end if

  if (.not. allocated(this%rhs_f)) then
    if (.not. associated(this%element_list)) then
      write(*,*) "ERROR: Call the constructor for projection or allocate element_list"
      return
    end if

    allocate(this%rhs_f(n_vertex_max,n_order+1,this%element_list%n_elements,n_tor,n_sample))
  
  end if
  
  this%rhs_f = 0.d0


  ! We now write multiple omp regions since this allows us to do reductions per
  ! variable instead of all-at-once, helping with the stack size requirements
  ! there is however some thread-creation overhead

  allocate(my_rhs(n_vertex_max,n_order+1,this%element_list%n_elements,n_tor,1))

  do i_f=1,n_sample

    i_group = this%f(i_f)%group

    if (.not. allocated(sim%groups(i_group)%particles)) cycle
  
    my_rhs = 0.d0

    !$omp parallel do default(none)                                            &
    !$omp shared(this, sim, n_sample, i_group, i_f)                            &
    !$omp private(x, xjac, HH, HH_s, HH_t, R_g, R_s, R_t, Z_g, Z_s, Z_t,       &
    !$omp         m, i, j, i_tor, index_ij, v, HP),                            &
    !$omp reduction(+:my_rhs) schedule(dynamic,10)
    ! Note that the stack per thread can become larger than the 2MB default in
    ! this case. You could try setting the OMP_STACKSIZE environment variable
    ! somewhat larger than the size of my_rhs
    do m=1,size(sim%groups(i_group)%particles,1)

      if (sim%groups(i_group)%particles(m)%i_elm .le. 0) cycle

      x(1:2) = sim%groups(i_group)%particles(m)%st
      x(3)   = sim%groups(i_group)%particles(m)%x(3)

      call interp_RZ(this%node_list,this%element_list,sim%groups(i_group)%particles(m)%i_elm,x(1),x(2),R_g,R_s,R_t,Z_g,Z_s,Z_t)
    
      xjac = R_s*Z_t - R_t*Z_s
      
      call basisfunctions(x(1), x(2), HH, HH_s, HH_t)
      call mode_moivre(x(3), HP)
            
      do i=1,n_vertex_max
        do j=1,n_order+1
       
          v = HH(i,j) * this%element_list%element(sim%groups(i_group)%particles(m)%i_elm)%size(i,j)

          v = v * this%f(i_f)%f(sim, i_group, sim%groups(i_group)%particles(m)) * sim%groups(i_group)%particles(m)%weight
       
          do i_tor = 1, n_tor

            my_rhs(j,i,sim%groups(i_group)%particles(m)%i_elm,i_tor,1) = &
            my_rhs(j,i,sim%groups(i_group)%particles(m)%i_elm,i_tor,1) + HP(i_tor) * v
          
          enddo

        enddo
      enddo
    enddo

    !$omp end parallel do
    this%rhs_f(:,:,:,:,i_f) = my_rhs(:,:,:,:,1)
  enddo
 
  deallocate(my_rhs)

end subroutine sample_rhs


!> Save an already-projected set to a vtk file with current parameters
subroutine save_to_vtk(this, sim)
  use mod_event
  !$ use omp_lib
  class(projection), intent(inout) :: this
  type(particle_sim), intent(inout)    :: sim
  integer :: i, ierr, n_proj
  real*8 :: t0, t1, ostart, oend
  character(len=120) :: filename

  this%extension = '.vtk'
  if (.not. allocated(this%vtk_grid)) then
    write(*,*) "ERROR: Trying to write unprepared vtk file"
    return
  end if

  if (len_trim(this%filename) .eq. 0) then
    filename = this%get_filename(sim%time)
  else
    filename = this%filename
  end if

  call cpu_time(t0)
  !$ ostart = omp_get_wtime()
  if (sim%my_id .eq. 0) then
    ! write only on the host

    n_proj = size(this%f)
    if (allocated(this%rhs)) n_proj = n_proj + size(this%rhs,5)

    call write_particle_distribution_to_vtk(this%node_list, this%element_list, &
      trim(filename), this%vtk_grid%nsub, min(n_proj,n_var), this%vtk_grid%xyz, this%vtk_grid%ien)
      
    write(*,*) "Written projection to ", trim(filename)
  end if
  call cpu_time(t1)
  !$ oend = omp_get_wtime()
#ifdef DEBUG
  if (sim%my_id .eq. 0 ) then
    write(*,"(A,2g12.5)") "writing cpu time", t1-t0
    !$ write(*,"(A,2g12.5)") "writing wall time", oend-ostart
  end if
#endif
end subroutine save_to_vtk




!> Action for projecting all particles and writing output to a hdf5 file
subroutine save_to_h5(this, sim)
  use mpi_mod
  use mod_event
  !$ use omp_lib
  class(projection), intent(inout)  :: this
  type(particle_sim), intent(inout)    :: sim
  integer :: my_id, ierr, n_proj
  character(len=120) :: filename
  real*8 :: t0, t1, ostart, oend

  this%extension = '.h5'

  if (len_trim(this%filename) .eq. 0) then
    filename = this%get_filename(sim%time)
  else
    filename = this%filename
  end if

  call cpu_time(t0)
  !$ ostart = omp_get_wtime()
  call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)
  if (my_id .eq. 0) then
    ! write only on the host
    n_proj = size(this%f)
    if (allocated(this%rhs)) n_proj = n_proj + size(this%rhs,5)
    call write_particle_distribution_to_h5(this%node_list, this%element_list, &
      trim(filename), min(n_proj,n_var), sim%time)

    write(*,*) "Written projection to ", trim(filename)
  end if
  call cpu_time(t1)
  !$ oend = omp_get_wtime()
#ifdef DEBUG
  if (my_id .eq. 0 ) then
    write(*,"(A,2g12.5)") "writing cpu time", t1-t0
    !$ write(*,"(A,2g12.5)") "writing wall time", oend-ostart
  end if
#endif
end subroutine save_to_h5



!> Project particles by weight onto the elements
!> Saves output in node_list%values(1)
!> The projection is done by solving
!> $$\int [p(x) u(x) + \lambda p'(x) u'(x)] dV =
!> \int \sum \delta(x_i-x) w_i u(x) dV $$
!> where p(x) is in the bernstein representation, $u(x)$ are the test functions
!> composed of two basis functions and $w_i$ is the particle weight.
!> A smoothing factor lambda is included
!> x is a vector (R,Z,phi) and dV is r dr dphi
!> divide by 1 or 2pi on both sides (LHS gets 2pi for n=0 mode, 1pi for other modes)
!>
!> See also [project_particles] and [mod_elt_matrix] for reference of the integration method
subroutine prepare_mumps_par(node_list, element_list, n_tor_local, i_tor_local,           &
                             this_mpi_comm_world, this_mpi_comm_n, this_mpi_comm_master,  &
                             mumps_par, filter, filter_hyper, filter_parallel,            &
                             skip_factorisation)
use phys_module, only : F0, TWOPI, mode
use data_structure
use basis_at_gaussian
use mod_basisfunctions
use mpi_mod
use, intrinsic :: ieee_exceptions
implicit none

type (type_node_list), intent(in)    :: node_list !< A copy of the node list which will be used to save variables
type (type_element_list), intent(in) :: element_list
integer, intent(in)                  :: n_tor_local, i_tor_local
integer, intent(in)                  :: this_mpi_comm_world, this_mpi_comm_n, this_mpi_comm_master
type (DMUMPS_STRUC), intent(inout)   :: mumps_par
real*8, intent(in)                   :: filter
real*8, intent(in)                   :: filter_hyper
real*8, intent(in)                   :: filter_parallel
logical, intent(in), optional        :: skip_factorisation

type (type_element)      :: element
type (type_node)         :: nodes(n_vertex_max)

real*8, allocatable      :: ELM(:,:)
real*8     :: wgauss2(n_gauss)
real*8, dimension(n_gauss,n_gauss) :: x_g,   x_s,   x_t,   x_ss,   x_tt,   x_st
real*8, dimension(n_gauss,n_gauss) :: y_g,   y_s,   y_t,   y_ss,   y_tt,   y_st
real*8, dimension(n_gauss,n_gauss) :: psi_g, psi_s, psi_t, psi_ss, psi_tt, psi_st
real*8     :: v, v_s, v_t, v_ss, v_st, v_tt, v_x, v_y, v_xx, v_yy, v_p
real*8     :: p, p_s, p_t, p_ss, p_st, p_tt, p_x, p_y, p_xx, p_yy, p_p
real*8     :: wst, area, volume, xjac, xjac_x, xjac_y, psi_x, psi_y
real*8     :: Bgrad_p, Bgrad_v_star, BB2
integer    :: i, j, k, l, m, in, im, ilarge, index_large_i, index_large_k, inode, knode
integer    :: nz_AA, n_AA, nz_bnd, i_elm, index_ij, index_kl, im_index, in_index, index1
integer    :: ms, mt, mp, my_id, my_id_n, my_id_master, ierr, MPI_COMM_MUMPS
logical    :: apply_dirichlet_condition
logical    :: halt(size(IEEE_USUAL,1))

! We need a separate communicator to be able to run multiple MUMPSes
call MPI_Comm_dup(this_mpi_comm_n, MPI_COMM_MUMPS, ierr)

! Initialise MUMPS
mumps_par%COMM = MPI_COMM_MUMPS
mumps_par%JOB  = -1
mumps_par%SYM  = 0
mumps_par%PAR  = 1

call DMUMPS(mumps_par)
call MPI_COMM_RANK(this_mpi_comm_world,  my_id, ierr)
call MPI_COMM_RANK(mumps_par%COMM,       my_id_n, ierr)

nz_AA = 4 * element_list%n_elements * (n_vertex_max * (n_order+1))**2
n_AA  = 2 * maxval(node_list%node(1:node_list%n_nodes)%index(4))

apply_dirichlet_condition = .false.

nz_bnd = 0
if (apply_dirichlet_condition) then
  do i=1,node_list%n_nodes
    if (node_list%node(i)%boundary .eq. 1) nz_bnd = nz_bnd + 4
    if (node_list%node(i)%boundary .eq. 2) nz_bnd = nz_bnd + 4
    if (node_list%node(i)%boundary .eq. 3) nz_bnd = nz_bnd + 8
    if (node_list%node(i)%boundary .eq. 4) nz_bnd = nz_bnd + 4
    if (node_list%node(i)%boundary .eq. 5) nz_bnd = nz_bnd + 4
    if (node_list%node(i)%boundary .eq. 9) nz_bnd = nz_bnd + 8
  enddo
endif

! Only perform the construction of the matrix on the host
if (my_id_n .eq. 0) then

  allocate(ELM(2*n_vertex_max*(n_order+1),2*n_vertex_max*(n_order+1)))

! Allocate space for elements

  allocate(mumps_par%A(nz_AA+nz_bnd),mumps_par%irn(nz_AA+nz_bnd),mumps_par%jcn(nz_AA+nz_bnd))

  mumps_par%irn = 0
  mumps_par%jcn = 0
  mumps_par%A   = 0.d0

! Copy wgauss into wgauss2 to get around gfortran not recognizing it as a shared
! thing https://groups.google.com/forum/#!topic/comp.lang.fortran/VKhoAm8m9KE
  wgauss2 = wgauss

  write(*,*) '**************************************************'
  write(*,'(A,i2,A)') ' * constructing particle projection matrix (n=',mode(i_tor_local),') *'
  write(*,*) '**************************************************'
  write(*,*) ' n_AA = ',n_AA, n_tor_local
  write(*,'(2i3,A,3e12.4)') my_id,my_id_n,'  filters       ',filter, filter_hyper, filter_parallel

  if (apply_dirichlet_condition) write(*,*) 'applying Dirichlet conditions'

!$omp parallel do default(none) &
!$omp shared(element_list, node_list, n_tor_local, i_tor_local,                     &
!$omp        H, H_s, H_t, H_ss, H_st, H_tt, Hz, Hz_p, mumps_par, wgauss2,           &
!$omp        filter, filter_hyper, filter_parallel, F0, my_id_master)               &
!$omp private(ELM, i_elm, element, nodes, i, j, k, l, ms, mt, in, im, mp,           &
!$omp         x_g, x_s, x_t, x_ss, x_st, x_tt,                                      &
!$omp         y_g, y_s, y_t, y_ss, y_st, y_tt,                                      &
!$omp         psi_g, psi_s, psi_t, psi_ss, psi_tt, psi_st,                          &
!$omp         v, v_s, v_t, v_ss, v_st, v_tt, v_x, v_y, v_xx, v_yy, v_p,             &
!$omp         p, p_s, p_t, p_ss, p_st, p_tt, p_x, p_y, p_xx, p_yy, p_p,             &
!$omp         wst, xjac, xjac_x, xjac_y, psi_x, psi_y, BB2, Bgrad_p, Bgrad_v_star,  &
!$omp         index_ij, index_kl, ilarge, in_index, im_index,                       &
!$omp         inode, index_large_i, knode, index_large_k)                           &
!$omp schedule(static)
do i_elm=1,element_list%n_elements
  
  ELM = 0.d0

  element = element_list%element(i_elm)
  do m=1,n_vertex_max
    nodes(m) = node_list%node(element%vertex(m))
  enddo

  ! Set up gauss points in this element
  x_g = 0.d0;   x_s = 0.d0;   x_t = 0.d0;   x_ss = 0.d0;   x_st = 0.d0;   x_tt = 0.d0
  y_g = 0.d0;   y_s = 0.d0;   y_t = 0.d0;   y_ss = 0.d0;   y_st = 0.d0;   y_tt = 0.d0
  psi_g = 0.d0; psi_s = 0.d0; psi_t = 0.d0; psi_ss = 0.d0; psi_st = 0.d0; psi_tt = 0.d0

  do i=1,n_vertex_max
    do j=1,n_order+1
      do ms=1, n_gauss
        do mt=1, n_gauss
          x_g(ms,mt)  = x_g(ms,mt)  + nodes(i)%x(1,j,1) * element%size(i,j) * H(i,j,ms,mt)
          x_s(ms,mt)  = x_s(ms,mt)  + nodes(i)%x(1,j,1) * element%size(i,j) * H_s(i,j,ms,mt)
          x_t(ms,mt)  = x_t(ms,mt)  + nodes(i)%x(1,j,1) * element%size(i,j) * H_t(i,j,ms,mt)

          x_ss(ms,mt) = x_ss(ms,mt) + nodes(i)%x(1,j,1) * element%size(i,j) * H_ss(i,j,ms,mt)
          x_st(ms,mt) = x_st(ms,mt) + nodes(i)%x(1,j,1) * element%size(i,j) * H_st(i,j,ms,mt)
          x_tt(ms,mt) = x_tt(ms,mt) + nodes(i)%x(1,j,1) * element%size(i,j) * H_tt(i,j,ms,mt)

          y_g(ms,mt)  = y_g(ms,mt)  + nodes(i)%x(1,j,2) * element%size(i,j) * H(i,j,ms,mt)
          y_s(ms,mt)  = y_s(ms,mt)  + nodes(i)%x(1,j,2) * element%size(i,j) * H_s(i,j,ms,mt)
          y_t(ms,mt)  = y_t(ms,mt)  + nodes(i)%x(1,j,2) * element%size(i,j) * H_t(i,j,ms,mt)

          y_ss(ms,mt) = y_ss(ms,mt) + nodes(i)%x(1,j,2) * element%size(i,j) * H_ss(i,j,ms,mt)
          y_st(ms,mt) = y_st(ms,mt) + nodes(i)%x(1,j,2) * element%size(i,j) * H_st(i,j,ms,mt)
          y_tt(ms,mt) = y_tt(ms,mt) + nodes(i)%x(1,j,2) * element%size(i,j) * H_tt(i,j,ms,mt)

          psi_g(ms,mt)  = psi_g(ms,mt)  + nodes(i)%values(1,j,1) * element%size(i,j) * H(i,j,ms,mt)
          psi_s(ms,mt)  = psi_s(ms,mt)  + nodes(i)%values(1,j,1) * element%size(i,j) * H_s(i,j,ms,mt)
          psi_t(ms,mt)  = psi_t(ms,mt)  + nodes(i)%values(1,j,1) * element%size(i,j) * H_t(i,j,ms,mt)

          psi_ss(ms,mt) = psi_ss(ms,mt) + nodes(i)%values(1,j,1) * element%size(i,j) * H_ss(i,j,ms,mt)
          psi_st(ms,mt) = psi_st(ms,mt) + nodes(i)%values(1,j,1) * element%size(i,j) * H_st(i,j,ms,mt)
          psi_tt(ms,mt) = psi_tt(ms,mt) + nodes(i)%values(1,j,1) * element%size(i,j) * H_tt(i,j,ms,mt)
        enddo
      enddo
    enddo
  enddo

  do ms=1, n_gauss
    do mt=1, n_gauss

      wst = wgauss2(ms)*wgauss2(mt)
      xjac =  x_s(ms,mt)*y_t(ms,mt) - x_t(ms,mt)*y_s(ms,mt)

      xjac_x  = (x_ss(ms,mt)*y_t(ms,mt)**2 - y_ss(ms,mt)*x_t(ms,mt)*y_t(ms,mt) - 2.d0*x_st(ms,mt)*y_s(ms,mt)*y_t(ms,mt)   &
              + y_st(ms,mt)*(x_s(ms,mt)*y_t(ms,mt) + x_t(ms,mt)*y_s(ms,mt))                                               &
              + x_tt(ms,mt)*y_s(ms,mt)**2 - y_tt(ms,mt)*x_s(ms,mt)*y_s(ms,mt)) / xjac

      xjac_y  = (y_tt(ms,mt)*x_s(ms,mt)**2 - x_tt(ms,mt)*y_s(ms,mt)*x_s(ms,mt) - 2.d0*y_st(ms,mt)*x_t(ms,mt)*x_s(ms,mt)   &
              + x_st(ms,mt)*(y_t(ms,mt)*x_s(ms,mt) + y_s(ms,mt)*x_t(ms,mt))                                               &
              + y_ss(ms,mt)*x_t(ms,mt)**2 - x_ss(ms,mt)*y_t(ms,mt)*x_t(ms,mt)) / xjac

      psi_x = (  y_t(ms,mt) * psi_s(ms,mt) - y_s(ms,mt) * psi_t(ms,mt)) / xjac
      psi_y = (- x_t(ms,mt) * psi_s(ms,mt) + x_s(ms,mt) * psi_t(ms,mt)) / xjac

      BB2 = 1.d0
      if (filter_parallel .gt. 0.d0) BB2 = (F0*F0 + psi_x * psi_x + psi_y * psi_y )/x_g(ms,mt)**2

      do mp = 1, n_plane

        do i=1,n_vertex_max
          do j=1,n_order+1

            do im = 1, n_tor_local

              im_index = i_tor_local + im - 1   ! i_tor_local is the starting index in HZ

              index_ij = 2*(n_order+1)*(i-1) + 2 * (j-1) + im   ! index in the ELM matrix

              v    = H(i,j,ms,mt)    * element%size(i,j) * HZ(im_index,mp)
              v_s  = H_s(i,j,ms,mt)  * element%size(i,j) * HZ(im_index,mp)
              v_t  = H_t(i,j,ms,mt)  * element%size(i,j) * HZ(im_index,mp)
              v_p  = H(i,j,ms,mt)    * element%size(i,j) * HZ_p(im_index,mp)

              v_ss = H_ss(i,j,ms,mt) * element%size(i,j) * HZ(im_index,mp)
              v_tt = H_tt(i,j,ms,mt) * element%size(i,j) * HZ(im_index,mp)
              v_st = H_st(i,j,ms,mt) * element%size(i,j) * HZ(im_index,mp)

              v_x = (  y_t(ms,mt) * v_s - y_s(ms,mt) * v_t) / xjac
              v_y = (- x_t(ms,mt) * v_s + x_s(ms,mt) * v_t) / xjac

              v_xx = (v_ss * y_t(ms,mt)**2 - 2.d0*v_st * y_s(ms,mt)*y_t(ms,mt) + v_tt * y_s(ms,mt)**2  &
                   + v_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                          &
                   + v_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) ) )  / xjac**2             &
                   - xjac_x * (v_s * y_t(ms,mt) - v_t * y_s(ms,mt)) / xjac**2

              v_yy = (v_ss * x_t(ms,mt)**2 - 2.d0*v_st * x_s(ms,mt)*x_t(ms,mt) + v_tt * x_s(ms,mt)**2  &
                   + v_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                          &
                   + v_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) ) )     / xjac**2          &
                   - xjac_y * (- v_s * x_t(ms,mt) + v_t * x_s(ms,mt) ) / xjac**2


              Bgrad_v_star = 0.d0
              if (filter_parallel .gt. 0.d0) Bgrad_v_star = ( F0 / x_g(ms,mt) * v_p  +  v_x  * psi_y - v_y  * psi_x ) / x_g(ms,mt)


              do k=1,n_vertex_max
                do l=1,n_order+1

                  do in = 1, n_tor_local

                    in_index = i_tor_local + in - 1

                    index_kl = 2*(n_order+1)*(k-1) + 2 * (l-1) + in   ! index in the ELM matrix

                    p   = h(k,l,ms,mt)     * element%size(k,l) * HZ(in_index,mp)
                    p_s = h_s(k,l,ms,mt)   * element%size(k,l) * HZ(in_index,mp)
                    p_t = h_t(k,l,ms,mt)   * element%size(k,l) * HZ(in_index,mp)
                    p_p = h(k,l,ms,mt)     * element%size(k,l) * HZ_p(in_index,mp)

                    p_ss = h_ss(k,l,ms,mt) * element%size(k,l) * HZ(in_index,mp)
                    p_tt = h_tt(k,l,ms,mt) * element%size(k,l) * HZ(in_index,mp)
                    p_st = h_st(k,l,ms,mt) * element%size(k,l) * HZ(in_index,mp)

                    p_x = (  y_t(ms,mt) * p_s - y_s(ms,mt) * p_t) / xjac
                    p_y = (- x_t(ms,mt) * p_s + x_s(ms,mt) * p_t) / xjac

                    p_xx = (p_ss * y_t(ms,mt)**2 - 2.d0*p_st * y_s(ms,mt)*y_t(ms,mt) + p_tt * y_s(ms,mt)**2  &
                         + p_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                          &
                         + p_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) ) )  / xjac**2             &
                         - xjac_x * (p_s * y_t(ms,mt) - p_t * y_s(ms,mt)) / xjac**2

                    p_yy = (p_ss * x_t(ms,mt)**2 - 2.d0*p_st * x_s(ms,mt)*x_t(ms,mt) + p_tt * x_s(ms,mt)**2  &
                         + p_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                          &
                         + p_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) ) )     / xjac**2          &
                         - xjac_y * (- p_s * x_t(ms,mt) + p_t * x_s(ms,mt) ) / xjac**2

                    Bgrad_p = 0.d0
                    if (filter_parallel .gt. 0.d0) Bgrad_p = ( F0 / x_g(ms,mt) * p_p +  p_x * psi_y - p_y * psi_x ) / x_g(ms,mt)

                    ELM(index_ij,index_kl) = ELM(index_ij,index_kl) &
                    
                                           + p * v * xjac * x_g(ms,mt) * wst &

                                           + filter          * (p_x * v_x + p_y * v_y) * xjac * x_g(ms,mt) * wst &

                                           + filter_hyper    * (v_xx + v_x/x_g(ms,mt) + v_yy)*(p_xx + p_x/x_g(ms,mt) + p_yy) * xjac * x_g(ms,mt) * wst &

                                           + filter_parallel * Bgrad_v_star * Bgrad_p / BB2 * xjac * x_g(ms,mt) * wst
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo

  ! Save contribution of this element in MUMPS format
  do i=1,n_vertex_max

    inode = element_list%element(i_elm)%vertex(i)
  
    do j=1,n_order+1

      do im =1, n_tor_local
    
        index_ij = 2*(n_order+1)*(i-1) + 2 * (j-1) + im   ! index in the ELM matrix

        index_large_i = 2*(node_list%node(inode)%index(j)-1) + im   ! base index in the main matrix

        do k=1,n_vertex_max
      
          knode = element_list%element(i_elm)%vertex(k)
        
          do l=1,n_order+1

            do in =1, n_tor_local
        
              index_kl = 2*(n_order+1)*(k-1) + 2 * (l-1) + in   ! index in the ELM matrix

              index_large_k = 2*(node_list%node(knode)%index(l)-1) + in   ! base index in the main matrix

             ! Explicitly calculate the index

              ilarge = in + (l-1) * 2 + (k-1)*2*(n_order+1) &
                      
                     + (im-1) * 2    * n_vertex_max*(n_order+1)       &
                     
                     + (j-1)  * 4 * n_vertex_max*(n_order+1)       &
                     
                     + (i-1)  * 4 * n_vertex_max*(n_order+1)**2    &
                     
                     + (i_elm-1)*(4 * (n_vertex_max*(n_order+1))**2 )

!$omp critical
              mumps_par%irn(ilarge) = index_large_i
              mumps_par%jcn(ilarge) = index_large_k
              mumps_par%A(ilarge)   = ELM(index_ij,index_kl) * TWOPI / real(n_plane,8)
!$omp end critical

            enddo
          enddo
        enddo
      enddo
    enddo
  enddo
enddo
!$omp end parallel do
ilarge = nz_AA

if (apply_dirichlet_condition) then

  do i=1,node_list%n_nodes
    
    if ((node_list%node(i)%boundary .eq. 2) .or. (node_list%node(i)%boundary .eq. 3) .or. &
        (node_list%node(i)%boundary .eq. 5) .or. (node_list%node(i)%boundary .eq. 9)) then

      do j=1,3,2             ! order
        do k=1,2             ! variables

          index1 = node_list%node(i)%index(j)

          ilarge = ilarge + 1

          mumps_par%irn(ilarge) = 2*(index1-1) + k
          mumps_par%jcn(ilarge) = 2*(index1-1) + k
          mumps_par%A(ilarge)   = 1.d12
        enddo
      enddo

    elseif ((node_list%node(i)%boundary .eq. 1) .or. (node_list%node(i)%boundary .eq. 3) .or. &
            (node_list%node(i)%boundary .eq. 4) .or. (node_list%node(i)%boundary .eq. 9)) then

      do j=1,2               ! order
        do k=1,2             ! variables

          index1 = node_list%node(i)%index(j)

          ilarge = ilarge + 1

          mumps_par%irn(ilarge) = 2*(index1-1) + k
          mumps_par%jcn(ilarge) = 2*(index1-1) + k
          mumps_par%A(ilarge)   = 1.d12
        enddo
      enddo
  
    endif
  enddo

endif

nz_AA = ilarge

end if

! Perform the analysis and factorisation with all nodes
mumps_par%JOB       = 4
mumps_par%n         = n_AA
mumps_par%nz        = nz_AA
mumps_par%icntl(2)  = 6 ! print diagnostics, statistics and warnings to stderr
mumps_par%icntl(4)  = 1 ! print errors(1), debug(2), much(3)
mumps_par%icntl(5)  = 0 ! assembled form
mumps_par%icntl(18) = 0 ! centralized input matrix (i.e. only on cpu 0)
mumps_par%icntl(7)  = 7 ! compute symmetric permutation (PORD or SCOTCH autoselect)
mumps_par%icntl(8)  = 8 ! scaling
mumps_par%icntl(14) = 80 ! memory relaxation parameter

if (present(skip_factorisation) .and. skip_factorisation) then
else
  call ieee_get_halting_mode(IEEE_USUAL, halt)
  call ieee_set_halting_mode(IEEE_USUAL, [.false., .false., .false.])
  call DMUMPS(mumps_par)
  call ieee_set_halting_mode(IEEE_USUAL, halt)
endif

end subroutine prepare_mumps_par


subroutine prepare_mumps_par_n0(node_list, element_list, n_tor_local, i_tor_local,           &
                                this_mpi_comm_world, this_mpi_comm_n, this_mpi_comm_master,  &
                                mumps_par, area, volume,                                     &
                                filter, filter_hyper, filter_parallel,                       &
                                integral_weights, skip_factorisation, do_zonal)
use phys_module, only : F0, TWOPI
use data_structure
use basis_at_gaussian
use mod_basisfunctions
use mpi_mod
use, intrinsic :: ieee_exceptions
implicit none

type (type_node_list), intent(in)    :: node_list !< A copy of the node list which will be used to save variables
type (type_element_list), intent(in) :: element_list
integer, intent(in)                  :: n_tor_local, i_tor_local
integer, intent(in)                  :: this_mpi_comm_world, this_mpi_comm_n, this_mpi_comm_master
type (DMUMPS_STRUC), intent(inout)   :: mumps_par
real*8, intent(in)                   :: filter
real*8, intent(in)                   :: filter_hyper
real*8, intent(in)                   :: filter_parallel
logical, intent(in), optional        :: do_zonal
logical, intent(in), optional        :: skip_factorisation
real*8, intent(out), dimension(:), allocatable :: integral_weights !< these will be filled with the volume of each basis function
real*8, intent(out)                  :: area, volume
type (type_element)      :: element
type (type_node)         :: nodes(n_vertex_max)

real*8, allocatable      :: ELM(:,:)
real*8     :: wgauss2(n_gauss)
real*8, dimension(n_gauss,n_gauss) :: x_g,   x_s,   x_t,   x_ss,   x_tt,   x_st
real*8, dimension(n_gauss,n_gauss) :: y_g,   y_s,   y_t,   y_ss,   y_tt,   y_st
real*8, dimension(n_gauss,n_gauss) :: psi_g, psi_s, psi_t, psi_ss, psi_tt, psi_st
real*8     :: v, v_s, v_t, v_ss, v_st, v_tt, v_x, v_y, v_xx, v_yy, v_p
real*8     :: p, p_s, p_t, p_ss, p_st, p_tt, p_x, p_y, p_xx, p_yy, p_p
real*8     :: wst, xjac, xjac_x, xjac_y, psi_x, psi_y
real*8     :: Bgrad_p, Bgrad_v_star, BB2
real*8     :: filter_n0, filter_hyper_n0, filter_parallel_n0, zonal_factor
integer    :: i, j, k, l, m, in, im, ilarge, index_large_i, index_large_k, inode, knode
integer    :: nz_AA, n_AA, nz_bnd, i_elm, index_ij, index_kl, im_index, in_index, index1, index2, index_rhs
integer    :: ms, mt, mp, my_id, my_id_n, my_id_master, ierr, MPI_COMM_MUMPS
logical    :: halt(size(IEEE_USUAL,1))
logical    :: apply_dirichlet_condition, apply_zonal
real*8, dimension(n_vertex_max,n_order+1) :: basisfunction_volume

! We need a separate communicator to be able to run multiple MUMPSes
call MPI_Comm_dup(this_mpi_comm_n, MPI_COMM_MUMPS, ierr)

! Initialise MUMPS
mumps_par%COMM = MPI_COMM_MUMPS
mumps_par%JOB  = -1
mumps_par%SYM  = 0
mumps_par%PAR  = 1

call DMUMPS(mumps_par)
call MPI_COMM_RANK(this_mpi_comm_world,  my_id, ierr)
call MPI_COMM_RANK(mumps_par%COMM,       my_id_n, ierr)

apply_zonal = .false.
if (present(do_zonal)) apply_zonal = do_zonal

apply_dirichlet_condition = .false.
if (apply_zonal)  apply_dirichlet_condition = .true.


nz_AA = 4 * element_list%n_elements * (n_vertex_max * (n_order+1))**2
n_AA =  2 * maxval(node_list%node(1:node_list%n_nodes)%index(4))

nz_bnd = 0
if (apply_dirichlet_condition) then
  do i=1,node_list%n_nodes
    if (node_list%node(i)%boundary .eq. 1) nz_bnd = nz_bnd + 4
    if (node_list%node(i)%boundary .eq. 2) nz_bnd = nz_bnd + 4
    if (node_list%node(i)%boundary .eq. 3) nz_bnd = nz_bnd + 8
    if (node_list%node(i)%boundary .eq. 4) nz_bnd = nz_bnd + 4
    if (node_list%node(i)%boundary .eq. 5) nz_bnd = nz_bnd + 4
    if (node_list%node(i)%boundary .eq. 9) nz_bnd = nz_bnd + 8
  enddo
endif

! Only perform the construction of the matrix on the host
if (my_id_n .eq. 0) then

  allocate(ELM(2*n_vertex_max*(n_order+1),2*n_vertex_max*(n_order+1)))
  allocate(mumps_par%A(nz_AA+nz_bnd),mumps_par%irn(nz_AA+nz_bnd),mumps_par%jcn(nz_AA+nz_bnd))
  
  mumps_par%irn = 0
  mumps_par%jcn = 0
  mumps_par%A   = 0.d0

  allocate(integral_weights(n_AA))
 
  integral_weights = 0.d0
  
! Copy wgauss into wgauss2 to get around gfortran not recognizing it as a shared
! thing https://groups.google.com/forum/#!topic/comp.lang.fortran/VKhoAm8m9KE
  wgauss2 = wgauss

  area   = 0.
  volume = 0.

  filter_n0             = filter 
  filter_hyper_n0       = filter_hyper 
  filter_parallel_n0    = filter_parallel 
  zonal_factor          = 1.d0

  write(*,*) '*************************************************'
  write(*,*) '* constructing particle projection matrix (n=0) *'
  write(*,*) '*************************************************'
  write(*,*)  ' n_AA = ',n_AA
  write(*,'(2I3,A,3e12.4)') my_id, my_id_n,'  filters (n=0) : ',filter_n0, filter_hyper_n0, filter_parallel_n0
  
  if (apply_zonal)               write(*,*) 'using n=0 zonal flow equations'
  if (apply_dirichlet_condition) write(*,*) 'applying Dirichlet conditions'

!$omp parallel do default(none) &
!$omp shared(element_list, node_list, n_tor_local, i_tor_local,                       &
!$omp        apply_dirichlet_condition, zonal_factor, apply_zonal,                    &
!$omp        H, H_s, H_t, H_ss, H_st, H_tt, Hz, Hz_p, mumps_par, wgauss2,             &
!$omp        filter_n0, filter_hyper_n0, filter_parallel_n0, integral_weights, F0, my_id_master) &
!$omp private(ELM, i_elm, element, nodes, i, j, k, l, ms, mt, in, im, mp,           &
!$omp         x_g, x_s, x_t, x_ss, x_st, x_tt,                                      &
!$omp         y_g, y_s, y_t, y_ss, y_st, y_tt,                                      &
!$omp         psi_g, psi_s, psi_t, psi_ss, psi_tt, psi_st,                          &
!$omp         v, v_s, v_t, v_ss, v_st, v_tt, v_x, v_y, v_xx, v_yy, v_p,             &
!$omp         p, p_s, p_t, p_ss, p_st, p_tt, p_x, p_y, p_xx, p_yy, p_p,             &
!$omp         wst, xjac, xjac_x, xjac_y, psi_x, psi_y, BB2, Bgrad_p, Bgrad_v_star,  &
!$omp         index_ij, index_kl, ilarge, in_index, im_index, basisfunction_volume, &
!$omp         inode, index_large_i, knode, index_large_k, index_rhs)                &
!$omp reduction(+:area,volume) schedule(static)
do i_elm=1,element_list%n_elements
  
  ELM = 0.d0

  element = element_list%element(i_elm)
  do m=1,n_vertex_max
    nodes(m) = node_list%node(element%vertex(m))
  enddo

  ! Set up gauss points in this element
  x_g = 0.d0;   x_s = 0.d0;   x_t = 0.d0;   x_ss = 0.d0;   x_st = 0.d0;   x_tt = 0.d0
  y_g = 0.d0;   y_s = 0.d0;   y_t = 0.d0;   y_ss = 0.d0;   y_st = 0.d0;   y_tt = 0.d0
  psi_g = 0.d0; psi_s = 0.d0; psi_t = 0.d0; psi_ss = 0.d0; psi_st = 0.d0; psi_tt = 0.d0

  do i=1,n_vertex_max
    do j=1,n_order+1
      do ms=1, n_gauss
        do mt=1, n_gauss
          x_g(ms,mt)  = x_g(ms,mt)  + nodes(i)%x(1,j,1) * element%size(i,j) * H(i,j,ms,mt)
          x_s(ms,mt)  = x_s(ms,mt)  + nodes(i)%x(1,j,1) * element%size(i,j) * H_s(i,j,ms,mt)
          x_t(ms,mt)  = x_t(ms,mt)  + nodes(i)%x(1,j,1) * element%size(i,j) * H_t(i,j,ms,mt)

          x_ss(ms,mt) = x_ss(ms,mt) + nodes(i)%x(1,j,1) * element%size(i,j) * H_ss(i,j,ms,mt)
          x_st(ms,mt) = x_st(ms,mt) + nodes(i)%x(1,j,1) * element%size(i,j) * H_st(i,j,ms,mt)
          x_tt(ms,mt) = x_tt(ms,mt) + nodes(i)%x(1,j,1) * element%size(i,j) * H_tt(i,j,ms,mt)

          y_g(ms,mt)  = y_g(ms,mt)  + nodes(i)%x(1,j,2) * element%size(i,j) * H(i,j,ms,mt)
          y_s(ms,mt)  = y_s(ms,mt)  + nodes(i)%x(1,j,2) * element%size(i,j) * H_s(i,j,ms,mt)
          y_t(ms,mt)  = y_t(ms,mt)  + nodes(i)%x(1,j,2) * element%size(i,j) * H_t(i,j,ms,mt)

          y_ss(ms,mt) = y_ss(ms,mt) + nodes(i)%x(1,j,2) * element%size(i,j) * H_ss(i,j,ms,mt)
          y_st(ms,mt) = y_st(ms,mt) + nodes(i)%x(1,j,2) * element%size(i,j) * H_st(i,j,ms,mt)
          y_tt(ms,mt) = y_tt(ms,mt) + nodes(i)%x(1,j,2) * element%size(i,j) * H_tt(i,j,ms,mt)

          psi_g(ms,mt)  = psi_g(ms,mt)  + nodes(i)%values(1,j,1) * element%size(i,j) * H(i,j,ms,mt)
          psi_s(ms,mt)  = psi_s(ms,mt)  + nodes(i)%values(1,j,1) * element%size(i,j) * H_s(i,j,ms,mt)
          psi_t(ms,mt)  = psi_t(ms,mt)  + nodes(i)%values(1,j,1) * element%size(i,j) * H_t(i,j,ms,mt)

          psi_ss(ms,mt) = psi_ss(ms,mt) + nodes(i)%values(1,j,1) * element%size(i,j) * H_ss(i,j,ms,mt)
          psi_st(ms,mt) = psi_st(ms,mt) + nodes(i)%values(1,j,1) * element%size(i,j) * H_st(i,j,ms,mt)
          psi_tt(ms,mt) = psi_tt(ms,mt) + nodes(i)%values(1,j,1) * element%size(i,j) * H_tt(i,j,ms,mt)

        enddo
      enddo
    enddo
  enddo

  basisfunction_volume = 0.d0

  do ms=1, n_gauss
    do mt=1, n_gauss

      wst = wgauss2(ms)*wgauss2(mt)
      xjac =  x_s(ms,mt)*y_t(ms,mt) - x_t(ms,mt)*y_s(ms,mt)

      xjac_x  = (x_ss(ms,mt)*y_t(ms,mt)**2 - y_ss(ms,mt)*x_t(ms,mt)*y_t(ms,mt) - 2.d0*x_st(ms,mt)*y_s(ms,mt)*y_t(ms,mt)   &
              + y_st(ms,mt)*(x_s(ms,mt)*y_t(ms,mt) + x_t(ms,mt)*y_s(ms,mt))                                               &
              + x_tt(ms,mt)*y_s(ms,mt)**2 - y_tt(ms,mt)*x_s(ms,mt)*y_s(ms,mt)) / xjac

      xjac_y  = (y_tt(ms,mt)*x_s(ms,mt)**2 - x_tt(ms,mt)*y_s(ms,mt)*x_s(ms,mt) - 2.d0*y_st(ms,mt)*x_t(ms,mt)*x_s(ms,mt)   &
              + x_st(ms,mt)*(y_t(ms,mt)*x_s(ms,mt) + y_s(ms,mt)*x_t(ms,mt))                                               &
              + y_ss(ms,mt)*x_t(ms,mt)**2 - x_ss(ms,mt)*y_t(ms,mt)*x_t(ms,mt)) / xjac

      psi_x = (  y_t(ms,mt) * psi_s(ms,mt) - y_s(ms,mt) * psi_t(ms,mt)) / xjac
      psi_y = (- x_t(ms,mt) * psi_s(ms,mt) + x_s(ms,mt) * psi_t(ms,mt)) / xjac

      BB2 = 1.d0
      if (filter_parallel_n0 .gt. 0.d0) BB2 = (F0*F0 + psi_x * psi_x + psi_y * psi_y )/x_g(ms,mt)

      area   = area   + xjac * wst
      volume = volume + TWOPI * x_g(ms,mt) * xjac * wst

      do i=1,n_vertex_max
        do j=1,n_order+1

          index_ij = 2*(n_order+1)*(i-1) + 2*(j-1) + 1   ! index in the ELM matrix

          v    = H(i,j,ms,mt)    * element%size(i,j) 
          v_s  = H_s(i,j,ms,mt)  * element%size(i,j)
          v_t  = H_t(i,j,ms,mt)  * element%size(i,j)
          v_p  = 0.d0

          v_ss = H_ss(i,j,ms,mt) * element%size(i,j)
          v_tt = H_tt(i,j,ms,mt) * element%size(i,j)
          v_st = H_st(i,j,ms,mt) * element%size(i,j)

          v_x = (  y_t(ms,mt) * v_s - y_s(ms,mt) * v_t) / xjac
          v_y = (- x_t(ms,mt) * v_s + x_s(ms,mt) * v_t) / xjac

          v_xx = (v_ss * y_t(ms,mt)**2 - 2.d0*v_st * y_s(ms,mt)*y_t(ms,mt) + v_tt * y_s(ms,mt)**2  &
               + v_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                          &
               + v_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) ) )  / xjac**2             &
               - xjac_x * (v_s * y_t(ms,mt) - v_t * y_s(ms,mt)) / xjac**2

          v_yy = (v_ss * x_t(ms,mt)**2 - 2.d0*v_st * x_s(ms,mt)*x_t(ms,mt) + v_tt * x_s(ms,mt)**2  &
               + v_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                          &
               + v_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) ) )     / xjac**2          &
               - xjac_y * (- v_s * x_t(ms,mt) + v_t * x_s(ms,mt) ) / xjac**2

          Bgrad_v_star = 0.d0
          if (filter_parallel_n0 .gt. 0.d0) Bgrad_v_star = ( F0 / x_g(ms,mt) * v_p  +  v_x  * psi_y - v_y  * psi_x ) / x_g(ms,mt)

          basisfunction_volume(i,j) = basisfunction_volume(i,j) +  v * TWOPI * x_g(ms,mt) * xjac * wst

          do k=1,n_vertex_max
            do l=1,n_order+1

              index_kl = 2*(n_order+1)*(k-1) + 2*(l-1) + 1   ! index in the ELM matrix

              p   = h(k,l,ms,mt)     * element%size(k,l) 
              p_s = h_s(k,l,ms,mt)   * element%size(k,l)
              p_t = h_t(k,l,ms,mt)   * element%size(k,l)
              p_p = 0.d0

              p_ss = h_ss(k,l,ms,mt) * element%size(k,l)
              p_tt = h_tt(k,l,ms,mt) * element%size(k,l)
              p_st = h_st(k,l,ms,mt) * element%size(k,l)

              p_x = (  y_t(ms,mt) * p_s - y_s(ms,mt) * p_t) / xjac
              p_y = (- x_t(ms,mt) * p_s + x_s(ms,mt) * p_t) / xjac

              p_xx = (p_ss * y_t(ms,mt)**2 - 2.d0*p_st * y_s(ms,mt)*y_t(ms,mt) + p_tt * y_s(ms,mt)**2  &
                   + p_s * (y_st(ms,mt)*y_t(ms,mt) - y_tt(ms,mt)*y_s(ms,mt) )                          &
                   + p_t * (y_st(ms,mt)*y_s(ms,mt) - y_ss(ms,mt)*y_t(ms,mt) ) )  / xjac**2             &
                   - xjac_x * (p_s * y_t(ms,mt) - p_t * y_s(ms,mt)) / xjac**2

              p_yy = (p_ss * x_t(ms,mt)**2 - 2.d0*p_st * x_s(ms,mt)*x_t(ms,mt) + p_tt * x_s(ms,mt)**2  &
                   + p_s * (x_st(ms,mt)*x_t(ms,mt) - x_tt(ms,mt)*x_s(ms,mt) )                          &
                   + p_t * (x_st(ms,mt)*x_s(ms,mt) - x_ss(ms,mt)*x_t(ms,mt) ) )     / xjac**2          &
                   - xjac_y * (- p_s * x_t(ms,mt) + p_t * x_s(ms,mt) ) / xjac**2

              Bgrad_p = 0.d0
              if (filter_parallel_n0 .gt. 0.d0) Bgrad_p = ( F0 / x_g(ms,mt) * p_p +  p_x * psi_y - p_y * psi_x ) / x_g(ms,mt)


              if (apply_zonal) then

                ELM(index_ij,index_kl)     = ELM(index_ij,index_kl)     + p * v        * xjac * x_g(ms,mt) * wst &

                                           + filter_n0       * (p_x * v_x + p_y * v_y) * xjac * x_g(ms,mt) * wst &

                                           + filter_hyper_n0 * (v_xx + v_x/x_g(ms,mt) + v_yy)*(p_xx + p_x/x_g(ms,mt) + p_yy) * xjac * x_g(ms,mt) * wst

                ELM(index_ij,index_kl+1)   = ELM(index_ij,index_kl+1) &

                                           + filter_n0  * (p_x * v_x + p_y * v_y) * xjac * x_g(ms,mt) * wst 


                ELM(index_ij+1,index_kl)   = ELM(index_ij+1,index_kl)   + p * v   * xjac * x_g(ms,mt) * wst * zonal_factor

                ELM(index_ij+1,index_kl+1) = ELM(index_ij+1,index_kl+1) &

                                           - filter_parallel_n0 * Bgrad_v_star * Bgrad_p / BB2 * xjac * x_g(ms,mt) * wst

              else

                ELM(index_ij,index_kl)     = ELM(index_ij,index_kl)     + p * v           * xjac * x_g(ms,mt) * wst &

                                           + filter_n0          * (p_x * v_x + p_y * v_y) * xjac * x_g(ms,mt) * wst &

                                           + filter_hyper_n0    * (v_xx + v_x/x_g(ms,mt) + v_yy)*(p_xx + p_x/x_g(ms,mt) + p_yy) * xjac * x_g(ms,mt) * wst &

                                           + filter_parallel_n0 * Bgrad_v_star * Bgrad_p / BB2 * xjac * x_g(ms,mt) * wst

                ELM(index_ij+1,index_kl+1) = ELM(index_ij+1,index_kl+1)  + p * v      * xjac * x_g(ms,mt) * wst       ! (dummy equation)
                
              endif
              
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo

  ! Save contribution of this element in MUMPS format
  do i=1,n_vertex_max

    inode = element_list%element(i_elm)%vertex(i)
  
    do j=1,n_order+1

      index_rhs = 2*(node_list%node(inode)%index(j)-1) + 1   ! base index in the main matrix

      integral_weights(index_rhs) = integral_weights(index_rhs) + basisfunction_volume(i,j)

      do im =1, 2
    
        index_ij = 2*(n_order+1)*(i-1) + 2 * (j-1) + im   ! index in the ELM matrix

        index_large_i = 2*(node_list%node(inode)%index(j)-1) + im   ! base index in the main matrix

        do k=1,n_vertex_max
      
          knode = element_list%element(i_elm)%vertex(k)
        
          do l=1,n_order+1

            do in =1, 2
        
              index_kl = 2*(n_order+1)*(k-1) + 2 * (l-1) + in   ! index in the ELM matrix

              index_large_k = 2*(node_list%node(knode)%index(l)-1) + in   ! base index in the main matrix

             ! Explicitly calculate the index

              ilarge = in + 2*(l-1) + 2*(k-1)*(n_order+1) &
                      
                     + 2*(im-1)* n_vertex_max*(n_order+1)       &
                     
                     + 4*(j-1) * n_vertex_max*(n_order+1)       &
                     
                     + 4*(i-1) * n_vertex_max*(n_order+1)**2    &
                     
                     + 4*(i_elm-1)*((n_vertex_max*(n_order+1))**2 )

              mumps_par%irn(ilarge) = index_large_i
              mumps_par%jcn(ilarge) = index_large_k
              mumps_par%A(ilarge)   = ELM(index_ij,index_kl) * TWOPI
!              mumps_par%A(ilarge)   = ELM(index_ij,index_kl)

            enddo
          enddo
        enddo
      enddo
    enddo
  enddo
enddo
!$omp end parallel do

! boundary conditions 

ilarge = nz_AA

if (apply_dirichlet_condition) then

  do i=1,node_list%n_nodes
    
    if ((node_list%node(i)%boundary .eq. 2) .or. (node_list%node(i)%boundary .eq. 3) .or. &
        (node_list%node(i)%boundary .eq. 5) .or. (node_list%node(i)%boundary .eq. 9)) then

      do j=1,3,2         ! order
        do k=1,2         ! variables

          index1 = node_list%node(i)%index(j)

          ilarge = ilarge + 1

          mumps_par%irn(ilarge) = 2*(index1-1) + k
          mumps_par%jcn(ilarge) = 2*(index1-1) + k
          mumps_par%A(ilarge)   = 1.d12
        enddo
      enddo

    elseif ((node_list%node(i)%boundary .eq. 1) .or. (node_list%node(i)%boundary .eq. 3) .or. &
            (node_list%node(i)%boundary .eq. 4) .or. (node_list%node(i)%boundary .eq. 9)) then

      do j=1,2           ! order
        do k=1,2         ! variables

          index1 = node_list%node(i)%index(j)

          ilarge = ilarge + 1

          mumps_par%irn(ilarge) = 2*(index1-1) + k
          mumps_par%jcn(ilarge) = 2*(index1-1) + k
          mumps_par%A(ilarge)   = 1.d12
        enddo
      enddo
  
    endif
  enddo

endif

nz_AA = ilarge

write(*,'(A,2e16.8)') 'area volume : ',area, volume
end if

! Perform the analysis and factorisation with all nodes
mumps_par%JOB       = 4
mumps_par%n         = n_AA
mumps_par%nz        = nz_AA
mumps_par%icntl(2)  = 6 ! print diagnostics, statistics and warnings to stderr
mumps_par%icntl(4)  = 1 ! print errors(1), debug(2), much(3)
mumps_par%icntl(5)  = 0 ! assembled form
mumps_par%icntl(18) = 0 ! centralized input matrix (i.e. only on cpu 0)
mumps_par%icntl(7)  = 7 ! compute symmetric permutation (PORD or SCOTCH autoselect)
mumps_par%icntl(8)  = 8 ! scaling
mumps_par%icntl(14) = 80 ! memory relaxation parameter

if (present(skip_factorisation) .and. skip_factorisation) then
else
  call ieee_get_halting_mode(IEEE_USUAL, halt)
  call ieee_set_halting_mode(IEEE_USUAL, [.false., .false., .false.])
  call DMUMPS(mumps_par)
  call ieee_set_halting_mode(IEEE_USUAL, halt)
endif

end subroutine prepare_mumps_par_n0


!> Save a particle distribution in a sort of minimized JOREK restart format.
!> This contains only: n_tor, n_period, n_fields->n_var, n_vertex_max, vertex, x, size, n_elements, values
!> neighbours, RCS_version
!> model=-1 to signify that the variables have not the expected meaning
subroutine write_particle_distribution_to_h5(node_list,element_list,filename,n_fields,time)
use data_structure
use hdf5
use hdf5_io_module
use tr_module
!$ use omp_lib
implicit none

!> Input parameters
type(type_node_list), intent(in)    :: node_list
type(type_element_list), intent(in) :: element_list
character*(*), intent(in)           :: filename
integer, intent(in) :: n_fields !< number of different particle groups to output
real*8, intent(in) :: time

 
#include "version.h"
! --- Local variables
integer :: i
character*50             :: version_control

integer(HID_T)     :: file_id
integer            :: ierr

! type_node, node_list%n_nodes
real(RKIND), allocatable :: t_x(:,:,:,:)                   ! n_coord_tor, n_order+1, n_dim
real(RKIND), allocatable :: t_values(:,:,:,:)              !       n_tor, n_order+1, n_fields

! element, element_list%n_elements
integer,     allocatable :: t_vertex(:,:)                ! n_vertex_max
integer,     allocatable :: t_neighbours(:,:)            ! n_vertex_max
real(RKIND), allocatable :: t_size(:,:,:)                ! n_vertex_max,n_order+1

! type_node, node_list%n_nodes
call tr_allocate(t_x,     1,node_list%n_nodes,1,n_coord_tor,1,n_order+1,1,n_dim,   "node_list%x",     CAT_UNKNOWN)
call tr_allocate(t_values,1,node_list%n_nodes,1,n_tor,      1,n_order+1,1,n_fields,"node_list%values",CAT_UNKNOWN)

! element_list%n_elements
call tr_allocate(t_vertex,    1,element_list%n_elements,1,n_vertex_max,"vertex",CAT_UNKNOWN)
call tr_allocate(t_neighbours,1,element_list%n_elements,1,n_vertex_max,"neighbours",CAT_UNKNOWN)
call tr_allocate(t_size,      1,element_list%n_elements,1,n_vertex_max,1,n_order+1,"size",CAT_UNKNOWN)

do i=1,node_list%n_nodes
   t_x(i,:,:,:)      = node_list%node(i)%x
   t_values(i,:,:,:) = node_list%node(i)%values(:,:,1:n_fields)
end do

do i=1,element_list%n_elements
  t_vertex(i,:)       = element_list%element(i)%vertex
  t_neighbours(i,:)   = element_list%element(i)%neighbours
  t_size(i,:,:)       = element_list%element(i)%size
end do

! -> Create and open HDF5 file
call HDF5_create(trim(filename),file_id,ierr)
if (ierr.ne.0) then
  write(*,*) ' ==> error for opening of HDF5 file',filename
end if
  
! -> Save version of revision control system
write(version_control,'(A)') trim(adjustl(RCS_VERSION))
version_control = trim(adjustl(version_control))
call HDF5_char_saving(file_id,version_control,"RCS_version"//char(0))

! -> Save parameters
call HDF5_integer_saving(file_id,-1,'jorek_model'//char(0)) ! Indicate that this is not a normal JOREK model
call HDF5_integer_saving(file_id,n_fields,'n_var'//char(0))
call HDF5_integer_saving(file_id,n_dim,'n_dim'//char(0))
call HDF5_integer_saving(file_id,n_order,'n_order'//char(0))
call HDF5_integer_saving(file_id,n_tor,'n_tor'//char(0))
call HDF5_integer_saving(file_id,n_coord_tor,'n_coord_tor'//char(0))
call HDF5_integer_saving(file_id,n_period,'n_period'//char(0))
call HDF5_integer_saving(file_id,n_vertex_max,'n_vertex_max'//char(0))
call HDF5_integer_saving(file_id,n_nodes_max,'n_nodes_max'//char(0))
call HDF5_integer_saving(file_id,n_elements_max,'n_elements_max'//char(0))

! -> 
call HDF5_integer_saving(file_id,node_list%n_nodes,'n_nodes'//char(0))
call HDF5_integer_saving(file_id,element_list%n_elements,'n_elements'//char(0))
call HDF5_integer_saving(file_id,node_list%n_dof,'n_dof'//char(0))

call HDF5_array4D_saving(file_id,t_x, &
     node_list%n_nodes,n_coord_tor,n_order+1,n_dim,'x'//char(0))
call HDF5_array4D_saving(file_id,t_values, &
     node_list%n_nodes,n_tor,n_order+1,n_fields,'values'//char(0))

call HDF5_array2D_saving_int(file_id,t_vertex, &
     element_list%n_elements,n_vertex_max,'vertex'//char(0))
call HDF5_array2D_saving_int(file_id,t_neighbours, &
     element_list%n_elements,n_vertex_max,'neighbours'//char(0))
call HDF5_array3D_saving(file_id,t_size, &
     element_list%n_elements,n_vertex_max,n_order+1,'size'//char(0))
call HDF5_real_saving(file_id,time,'t_now'//char(0))

! -> close file
call HDF5_close(file_id)
end subroutine write_particle_distribution_to_h5


!> Helper subroutine to write a particle distribution, saved in variables 1:n,
!> to `filename`. `nsub` is the number of subdivisions to make per element.
!> Should be called only by process 1
subroutine write_particle_distribution_to_vtk(node_list,element_list,filename,nsub,n_fields,xyz,ien)
use data_structure
use phys_module, only: mode
use mod_vtk
use mod_interp
use mod_basisfunctions
!$ use omp_lib
implicit none

!> Input parameters
type(type_node_list), intent(in)    :: node_list
type(type_element_list), intent(in) :: element_list
character*(*), intent(in)           :: filename
integer, intent(in) :: nsub !< Number of subdivisions of each element
integer, intent(in) :: n_fields !< number of different particle groups to output
real*4, intent(in)  :: xyz(:,:)
integer, intent(in) :: ien(:,:)

integer :: nnos, i, j, k, l, m, inode, ivar
real*4, allocatable :: scalars(:,:), vectors(:,:,:)
integer :: n_scalars, n_vectors = 0
character*12, allocatable :: vector_names(:), scalar_names(:)
real*8 :: s, t
real*8 :: P, P_s, P_t, P_st, P_ss, P_tt

integer, parameter :: etype = 9 ! for vtk_quad

n_scalars = n_tor * n_fields
nnos = nsub*nsub*element_list%n_elements

allocate(scalars(nnos,n_scalars),vectors(nnos,3,n_vectors))
allocate(scalar_names(n_scalars),vector_names(n_vectors))

do i=1,n_fields
  write(scalar_names(n_tor*(i-1)+1),'(A,i0.2)') "rho_", i
  do j=1,(n_tor-1)/2
    write(scalar_names(n_tor*(i-1)+2*j),"(A4,i0.2,A4,i0.2)") "rho_", i, "_cos", mode(2*j)
    write(scalar_names(n_tor*(i-1)+2*j+1),"(A4,i0.2,A4,i0.2)") "rho_", i, "_sin", mode(2*j+1)
  end do
end do

scalars = 0.e0
vectors = 0.e0

! Create points for each element
!$omp parallel do default(none) &
!$omp shared(element_list,nsub,node_list,n_fields,scalars) &
!$omp private(i,j,k,l,m,inode,ivar,s,t,P, P_s, P_t, P_st, P_ss, P_tt) schedule(static)
do i=1,element_list%n_elements
  do j=1,n_fields
    do k=1,n_tor
      ivar = (j-1)*n_tor + k
      do l=1,nsub
        s = real(l-1,8)/real(nsub-1,8)
        do m=1,nsub
          t = real(m-1,8)/real(nsub-1,8)
          inode = (i-1)*nsub*nsub+(l-1)*nsub+m
          call interp(node_list, element_list, i, j, k, s, t, P, P_s, P_t, P_st, P_ss, P_tt)
          scalars(inode,ivar) = real(P,4)
        end do
      end do
    end do
  end do
end do ! n_elements
!$omp end parallel do

! ------------- Write to VTK
call write_vtk(filename,xyz, ien, etype, scalar_names, scalars, vector_names, vectors)

end subroutine write_particle_distribution_to_vtk

end module mod_project_particles
