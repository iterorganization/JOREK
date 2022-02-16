!> mod_kinetic_relativistic_test contains variables and
!> procedures for testing the mod_kinetic_relativistic
module mod_kinetic_relativistic_test
use fruit
use mod_particle_types
use constants, only: EL_CHG,ATOMIC_MASS_UNIT,SPEED_OF_LIGHT,PI,TWOPI
implicit none

private
public :: run_fruit_kinetic_relativistic

!> Variables ------------------------------------------
integer*1,parameter :: q=-1                  !< electron charge in [A]
real*8,parameter    :: tol_basis=7.5d-15      !< tolerance for success orbital basis
real*8,parameter    :: mass=5.48579909065d-4 !< mass in AMU
real*8,parameter    :: E_kin=1.d1            !< kinetic energy MeV
real*8,parameter    :: I_E=1.25              !< E-field intensity [V/m]
real*8,parameter    :: I_B=3.d0              !< B-field intensity [B]
!> polar and azimuthal angle intervals
real*8,dimension(2),parameter :: thetachi_lowbnd=(/0.d0,0.d0/)
real*8,dimension(2),parameter :: thetachi_uppbnd=(/PI,TWOPI/)
type(particle_kinetic_relativistic) :: particle
!> electric and magnetic field in cartesian coordinates
real*8,dimension(3) :: E_field_cart,B_field_cart
!> orbital basis in cartesian coordinates
real*8,dimension(3) :: T_cart_sol,N_cart_sol,B_cart_sol

!>-----------------------------------------------------

contains

!> Fruit basket ---------------------------------------
!> run_fruit_kinetic_relativistic executes the feature
!> set-up, the tests and the feature tear-down
subroutine run_fruit_kinetic_relativistic()
  implicit none

  write(*,'(/A)') "  ... setting-up: kinetic relativistic tests"
  call setup
  write(*,'(/A)') "  ... running: kinetic relativistic tests"
  call test_relativistic_kinetic_orbital_basis
  write(*,'(/A)') "  ... tearing-up: kinetic relativistic tests"
  call teardown
end subroutine run_fruit_kinetic_relativistic

!> Set-up and tear-down -------------------------------
!> initialize kinetic relativistic test features
subroutine setup()
  use mod_math_operators, only: cross_product
  use mod_gnu_rng,        only: gnu_rng_interval
  implicit none
  !> variables
  real*8 :: gam                   !< relativistic factor
  real*8,dimension(2) :: thetachi !< polar and azimuthal angles

  !> initialise particle
  gam = sqrt(((E_kin*1.d6*EL_CHG)/(mass*ATOMIC_MASS_UNIT*&
  SPEED_OF_LIGHT*SPEED_OF_LIGHT))**2.d0 - 1.d0)
  call gnu_rng_interval(2,thetachi_lowbnd,thetachi_uppbnd,thetachi)
  particle%p = mass*gam*gam*(/sin(thetachi(1))*cos(thetachi(2)),&
  sin(thetachi(1))*sin(thetachi(2)),cos(thetachi(1))/)
  particle%q=q

  !> initialise electring and magnetic field
  call gnu_rng_interval(2,thetachi_lowbnd,thetachi_uppbnd,thetachi)
  E_field_cart = I_E*(/sin(thetachi(1))*cos(thetachi(2)),&
  sin(thetachi(1))*sin(thetachi(2)),cos(thetachi(1))/) !< cartesian E-field
  call gnu_rng_interval(2,thetachi_lowbnd,thetachi_uppbnd,thetachi)
  B_field_cart = I_B*(/sin(thetachi(1))*cos(thetachi(2)),&
  sin(thetachi(1))*sin(thetachi(2)),cos(thetachi(1))/) !< cartesian B-field

  !> compute cartesian orbital basis
  gam = sqrt(1.d0+(dot_product(particle%p,particle%p)/&
  (mass*SPEED_OF_LIGHT))**2.d0)
  T_cart_sol = particle%p/norm2(particle%p)
  N_cart_sol = cross_product(particle%p/(mass*gam),B_field_cart) + &
  E_field_cart - T_cart_sol*(dot_product(T_cart_sol,E_field_cart))
  N_cart_sol = N_cart_sol/norm2(N_cart_sol)
  B_cart_sol = cross_product(particle%p/(mass*gam),&
  cross_product(particle%p/(mass*gam),B_field_cart) + E_field_cart)
  B_cart_sol = B_cart_sol/norm2(B_cart_sol)
end subroutine setup

!> clean-up kinetic realativistic test features
subroutine teardown()
  implicit none
  !> clean up particle
  particle%p = 0.d0; particle%q = 0;
  !> clean-up fields
  E_field_cart = 0.d0; B_field_cart = 0.d0;
  !> clean-up orbital basis
  T_cart_sol = 0.d0; N_cart_sol = 0.d0; B_cart_sol = 0.d0;
end subroutine teardown

!> Tests ----------------------------------------------
!> Test cartesian orbital basis for relativistic 
!> full orbit particles
subroutine test_relativistic_kinetic_orbital_basis()
  use mod_kinetic_relativistic, only: compute_relativistic_kinetic_orbital_basis_cartesian
  implicit none
  !> variables
  real*8,dimension(3) :: T_cart,N_cart,B_cart

  !> compute the orbital basis
  call compute_relativistic_kinetic_orbital_basis_cartesian(&
  particle,mass,E_field_cart,B_field_cart,T_cart,N_cart,B_cart)
  !> test orthonormality
  call test_orthonormality_basis(T_cart,N_cart,B_cart,tol_basis)
  !> test solution
  call assert_equals(T_cart,T_cart_sol,3,tol_basis,& 
  "Error relativistic kinetic orbital basis: T direction mismatch!")
  call assert_equals(N_cart,N_cart_sol,3,tol_basis,& 
  "Error relativistic kinetic orbital basis: N direction mismatch!")
  call assert_equals(B_cart,B_cart_sol,3,tol_basis,& 
  "Error relativistic kinetic orbital basis: B direction mismatch!")
end subroutine test_relativistic_kinetic_orbital_basis

!> Tools ----------------------------------------------
!> check the orthonormality of a basis
subroutine test_orthonormality_basis(v1,v2,v3,tol)
  implicit none
  real*8 :: tol
  real*8,dimension(3),intent(in) :: v1,v2,v3
  call assert_equals(dot_product(v1,v1),1.d0,tol,&
  "Error basis orthonormality (double): v1 is not normalized!")
  call assert_equals(dot_product(v2,v2),1.d0,tol,&
  "Error basis orthonormality (double): v2 is not normalized!")
  call assert_equals(dot_product(v3,v3),1.d0,tol,&
  "Error basis orthonormality (double): v3 is not normalized!")
  call assert_equals(dot_product(v1,v2),0.d0,tol,&
  "Error basis orthonormality (double): v1 and v2 are not orthogonal!")
  call assert_equals(dot_product(v1,v3),0.d0,tol,&
  "Error basis orthonormality (double): v1 and v3 are not orthogonal!")
  call assert_equals(dot_product(v2,v3),0.d0,tol,&
  "Error basis orthonormality (double): v2 and v3 are not orthogonal!")
end subroutine test_orthonormality_basis

!>-----------------------------------------------------
end module mod_kinetic_relativistic_test

