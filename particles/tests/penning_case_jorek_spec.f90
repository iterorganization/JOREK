!> Verify that [[mod_penning_case_jorek]] produces the right fields (those of the
!> penning trap) on a few standard grids, with fields%interp.
module penning_case_jorek_spec
use data_structure
use projection_helpers
use fruit
use mod_penning_case
use mod_fields_linear
use mod_interp, only: interp_RZ
use mod_penning_case_jorek
implicit none
contains

!> Actions to perform before any of these tests
subroutine setup_penning_case_jorek_spec
  use basis_at_gaussian, only: initialise_basis
  use phys_module, only: central_mass, central_density
  call initialise_basis !< Calculate the basis functions at the gaussian points
  ! Set these two to have a usable t_norm later
  central_mass = 2.d0 !< D
  central_density = 1.d0 !< 10^20 at core
end subroutine setup_penning_case_jorek_spec

subroutine test_square_10_10
  type(jorek_fields_interp_linear) :: f
  allocate(f%node_list, f%element_list)
  call default_square_grid(f%node_list, f%element_list, 10)
  call jorek_penning_fields(f%node_list, f%element_list)
  call verify_solution(f)
  write(*,*) ' completed test_square_10_10'
end subroutine test_square_10_10

subroutine test_polar_30_32
  type(jorek_fields_interp_linear) :: f
  allocate(f%node_list, f%element_list)
  call default_polar_grid(f%node_list, f%element_list, 32)
  call jorek_penning_fields(f%node_list, f%element_list)
  call verify_solution(f)
  write(*,*) ' completed test_polar_30_32'
end subroutine test_polar_30_32

! flux aligned grid does not work since there is no axis in the domain (and there aren't really flux surfaces in any useful way)




!> Given a fields type, verify the solution at time 0 against the reference case
subroutine verify_solution(f)
  use phys_module, only: F0
  type(jorek_fields_interp_linear), intent(in) :: f
  type(case_penning_cylindrical), parameter :: ref = case_penning_cylindrical()

  real*8 :: E(3), B(3), psi, U
  real*8 :: E_ref(3), B_ref(3)
  real*8 :: R, Z
  integer :: i_elm

  do i_elm=1,f%element_list%n_elements
    call f%calc_EBpsiU(0.d0, i_elm, [0.5d0, 0.5d0], 0.d0, E, B, psi, U)
    call interp_RZ(f%node_list,f%element_list,i_elm,0.5d0,0.5d0,R,Z)

    E_ref = ref%E([R,Z,0.d0], 0.d0)
    call assert_equals(E_ref(1), E(1), 1d-10, 'E_r')
    call assert_equals(E_ref(2), E(2), 3d-11, 'E_z')
    call assert_equals(E_ref(3), E(3), 0.d0, 'E_phi')
    B_ref = ref%B([R,Z,0.d0], 0.d0)
    call assert_equals(B_ref(1), B(1), 2d-11, 'B_r')
    call assert_equals(B_ref(2), B(2), 2d-11, 'B_z')
    call assert_equals(B_ref(3), B(3), 2*F0, 'B_phi') ! includes the cheat factor
  end do
end subroutine verify_solution
end module penning_case_jorek_spec
