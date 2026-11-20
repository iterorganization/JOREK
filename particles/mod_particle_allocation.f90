!> This module contains functions to allocate sim%groups(:)%particles arrays based on the input from
!> part_group_configs, with the particles in each group distributed evenly across the mpi processes.

!> TODO:
!> In the future this module should use the different initialization schemes defined
!> in mod_particle_group_init.f90 (or similar) once that system has been reworked.
module mod_particle_allocation
use mpi
use mod_particle_sim
use mod_particle_types

implicit none

public allocate_particles_for_sim, allocate_particles_for_group, calc_n_particles_per_mpi, calc_n_particles_per_mpi_array

contains



  !> allocates all particle arrays for a simulation dependent on the configured particle group properties
  subroutine allocate_particles_for_sim(sim)
    use phys_module, only: part_group_configs, n_part_groups
    implicit none
    class(particle_sim), intent(inout)       :: sim
    integer                                  :: i, ierr
    integer,    dimension(:),    allocatable :: n_particles_per_mpi
    character(len=:),            allocatable :: particle_type_str

    do i=1, n_part_groups
      if (.not. allocated(sim%groups(i)%particles)) then

        ! calculate load balancing ------
        n_particles_per_mpi = calc_n_particles_per_mpi_array(int(sim%groups(i)%n_particles), sim%n_mpi)

        ! putting particle_type in the right format for allocate_particle_arrray function
        ! this is so that the function can be used both here and in io
        if (allocated(particle_type_str)) deallocate(particle_type_str)
        particle_type_str = trim(part_group_configs(i)%type)

        ! allocating particle arrays ---------
        call allocate_particles_for_group(sim, i, particle_type_str, n_particles_per_mpi(sim%my_id+1))

        ! initialize them to empty ---------- (temporary, to be replaced by a more general initialization function)
        call initialize_particle_list_to_zero(n_particles_per_mpi(sim%my_id+1),sim%groups(i)%particles,ierr)
      endif !if allocated
    enddo

    if (allocated(n_particles_per_mpi)) deallocate(n_particles_per_mpi)
    if (allocated(particle_type_str)) deallocate(particle_type_str)
    
  end subroutine allocate_particles_for_sim

  !> allocate an empty array of a specified particle type for a specific particle group for each mpi process
  subroutine allocate_particles_for_group(sim, group_num, particle_type_str, n_particles_per_mpi, mpi_comm_loc)
    use mod_particle_types

    implicit none

    type(particle_sim),            intent(inout) :: sim
    integer,                       intent(in)    :: group_num
    character(len=:), allocatable, intent(in)    :: particle_type_str
    integer,                       intent(in)    :: n_particles_per_mpi
    integer, optional,             intent(in)    :: mpi_comm_loc
    integer                                      :: errorcode, ierr, i

    select case (trim(particle_type_str))
    case ("particle_kinetic")
      allocate(particle_kinetic::sim%groups(group_num)%particles(n_particles_per_mpi))
    case ("particle_kinetic_leapfrog")
      allocate(particle_kinetic_leapfrog::sim%groups(group_num)%particles(n_particles_per_mpi))
    case ("particle_gc")
      allocate(particle_gc::sim%groups(group_num)%particles(n_particles_per_mpi))
    case ("particle_gc_vpar")
      allocate(particle_gc_vpar::sim%groups(group_num)%particles(n_particles_per_mpi))
    case ("particle_gc_Qin")
      allocate(particle_gc_Qin::sim%groups(group_num)%particles(n_particles_per_mpi))
    case ("particle_fieldline")
      allocate(particle_fieldline::sim%groups(group_num)%particles(n_particles_per_mpi))
    case ("particle_kinetic_relativistic")
      allocate(particle_kinetic_relativistic::sim%groups(group_num)%particles(n_particles_per_mpi))
    case ("particle_gc_relativistic")
      allocate(particle_gc_relativistic::sim%groups(group_num)%particles(n_particles_per_mpi))
    case default
      write(*,*) "Error: missing type name declaration ",trim(particle_type_str)," for reading: ABORT!"
      if (present(mpi_comm_loc)) then
        call MPI_Abort(mpi_comm_loc,errorcode,ierr)
      else
        stop
      endif
    end select

  end subroutine allocate_particles_for_group

  !> function ran locally on each mpi process which calcuates the number of particles 
  !> to distribute to it from a total of n_particles
  function calc_n_particles_per_mpi(n_particles, n_mpi, my_id) result(n_particles_per_mpi)

    implicit none
    integer, intent(in)            :: n_particles
    integer, intent(in)            :: n_mpi
    integer, intent(in)            :: my_id
    integer                        :: n_particles_per_mpi

    n_particles_per_mpi = n_particles/n_mpi          ! integer division truncates towards 0

    !> the remainder R of the division is then divided evenly amongst the first R processes
    if (my_id .lt. mod(n_particles,n_mpi)) then 
      n_particles_per_mpi = n_particles_per_mpi + 1 
    end if 

  end function calc_n_particles_per_mpi

  !> calculate how to distribute the particles in a group across the mpi processes when importing
  !> returns an array containing the number of particles for each mpi processor 
  function calc_n_particles_per_mpi_array(n_particles_tot, n_mpi) result(n_particles_per_mpi_array)
    implicit none

    integer,  intent(in)                  :: n_particles_tot             !< total number of particles for a particle group
    integer,  intent(in)                  :: n_mpi                       !< number of mpi processes being used
    integer,  dimension(:), allocatable   :: n_particles_per_mpi_array   !< output
    integer                               :: i

    allocate(n_particles_per_mpi_array(n_mpi))

    do i=1, n_mpi
      n_particles_per_mpi_array(i) = calc_n_particles_per_mpi(n_particles_tot, n_mpi, i-1) 
    enddo
    
  end function calc_n_particles_per_mpi_array



end module mod_particle_allocation