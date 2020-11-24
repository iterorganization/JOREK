!> Extension to JOREK of the penning test case.
!> Use here the parameters defined in [[mod_penning_case]] and project them onto a grid
!> to test the accuracy of the pusher in JOREK fields.
module mod_penning_case_jorek
  use mod_penning_case, only: charge, mass, omega_e, omega_b, epsilon
  use data_structure, only: type_node_list, type_element_list
  use mod_parameters, only: n_order
  implicit none
  private
  public jorek_penning_fields
contains

subroutine jorek_penning_fields(node_list, element_list)
  use projection_helpers, only: prepare_mumps_par, calc_rhs_f
  use phys_module, only: F0, central_mass, central_density
  use constants, only: mass_proton, mu_zero, el_chg, atomic_mass_unit
  use mpi
  implicit none

  include 'dmumps_struc.h'        ! MUMPS include files defining its datastructure

  type(type_node_list), intent(inout) :: node_list
  type(type_element_list), intent(inout) :: element_list

  type(DMUMPS_STRUC) :: p
  integer :: i, k, index
  real*8 :: t_norm, qom, B0, Phi0

  !> Note that we have to cheat a little bit here:
  !> We need F0 nonzero for the electric potential, but we need F0 = 0 to not have
  !> a toroidal magnetic field. Set F0 to a very small number, and the potential
  !> to a very high number to get an asymptotically correct answer (at the expense
  !> of accuracy in the field representation of U due to floating-point accuracy)
  real*8, parameter :: F0_CHEAT_FACTOR = 1d-30

  F0 = F0_CHEAT_FACTOR ! Set a very small F0 to have nearly no toroidal magnetic field

  ! Initialize fields in JOREK units
  t_norm  = sqrt(mu_zero * mass_proton * central_mass * central_density * 1.d20) ! 1 jorek time unit in seconds
  qom     = real(charge) * el_chg / (mass * atomic_mass_unit)
  B0      = omega_b/qom ! In T
  Phi0    = epsilon*omega_e**2/qom/2.d0*t_norm ! In JOREK units: E_SI*t_norm

  call prepare_mumps_par(node_list, element_list, 0, 0, MPI_COMM_WORLD, MPI_COMM_WORLD, MPI_COMM_WORLD, &
                         p, filter=0.d0, filter_hyper=0.d0, filter_parallel=0.d0, skip_factorisation=.false.)

  allocate(p%rhs(p%n))

  ! Project manually
  p%JOB = 3
  p%icntl(21) = 0 ! solution is available only on host
  p%icntl(4)  = 1 ! print only errors
  ! Setup RHS by integrating manually
  p%rhs = 0.d0

  ! For U
  ! the shape of the function is fixed, and we can calculate the scale factor after
  ! projecting the function.
  call calc_rhs_f(node_list,element_list,penning_U,p%rhs)
  call DMUMPS(p)
  do i=1,node_list%n_nodes
    do k=1,n_order+1
      index = node_list%node(i)%index(k)
      ! Scale JOREK fields correctly and apply F0_CHEAT_FACTOR
      ! p%rhs contains R^2 - 2 Z^2
      node_list%node(i)%values(1,k,2) = Phi0 * p%rhs(index) / F0_CHEAT_FACTOR
    enddo
  enddo

  ! For psi
  call calc_rhs_f(node_list,element_list,penning_psi,p%rhs)
  call DMUMPS(p)
  do i=1,node_list%n_nodes
    do k=1,n_order+1
      index = node_list%node(i)%index(k)
      ! p^rhs contains -R^2/2
      node_list%node(i)%values(1,k,1) = B0*p%rhs(index)
    enddo
  enddo

  p%JOB=-2
  call DMUMPS(p)

  node_list%node(i)%deltas = 0.d0
end subroutine jorek_penning_fields

!> E = -Grad F0*U, U = Phi0(R^2 - 2 Z^2) (see model001/initial_conditions.f90 for reference)
function penning_U(R, Z)
  real*8, intent(in) :: R, Z
  real*8 :: penning_U
  penning_U = R**2 - 2.d0*Z**2
end function penning_U

!> Projection function for psi to yield the right magnetic field.
!> B has only a z component, so F = 0, Psi(R) = -B0/2 R^2
function penning_psi(R, Z)
  real*8, intent(in) :: R, Z
  real*8 :: penning_psi
  penning_psi = -0.5d0 * R**2
end function penning_psi
end module mod_penning_case_jorek
