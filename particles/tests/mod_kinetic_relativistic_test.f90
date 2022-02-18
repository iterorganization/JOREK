!> mod_kinetic_relativistic_test contains variables and
!> procedures for testing the mod_kinetic_relativistic
module mod_kinetic_relativistic_test
use fruit
use mod_fields_analytical, only: fields_analytical
use mod_particle_types,    only: particle_kinetic_relativistic
use constants,             only: EL_CHG,ATOMIC_MASS_UNIT,SPEED_OF_LIGHT,PI,TWOPI
implicit none

private
public :: run_fruit_kinetic_relativistic

!> Variables ------------------------------------------
integer*1,parameter :: q=-1                  !< electron charge in [A]
integer,parameter   :: i_elm_zero=0
real*8,parameter    :: tol_real8=5.d-16
real*8,parameter    :: tol_basis=7.5d-15     !< tolerance for success orbital basis
real*8,parameter    :: mass=5.48579909065d-4 !< mass in AMU
real*8,parameter    :: time_sol=0.d0
real*8,dimension(2),parameter :: st_zero=0.d0
!> polar and azimuthal angle intervals
type(fields_analytical)             :: fields_sol
type(particle_kinetic_relativistic) :: particle
integer :: sign_theta_sol
real*8 :: phi_sol
!> orbital basis in cartesian coordinates
real*8,dimension(2) :: p_gc_sol !< guiding center momentum
real*8,dimension(3) :: pThetaChi_sol !< momentum, pitch and gyro angles
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
  call test_momentum_relat_kinetic_to_relat_gc
  call test_kinetic_relat_momentum_spherical_cart
  write(*,'(/A)') "  ... tearing-up: kinetic relativistic tests"
  call teardown
end subroutine run_fruit_kinetic_relativistic

!> Set-up and tear-down -------------------------------
!> initialize kinetic relativistic test features
subroutine setup()
  use mod_math_operators,             only: cross_product
  use mod_coordinate_transforms,      only: vector_cylindrical_to_cartesian
  use mod_gnu_rng,                    only: gnu_rng_interval
  use mod_sampling,                   only: sample_uniform_sphere_corona_rthetaphi
  use mod_particle_common_test_tools, only: EThetaChi_RE_lowbnd
  use mod_particle_common_test_tools, only: EThetaChi_RE_uppbnd
  use mod_particle_common_test_tools, only: RZ0_lowbnd,RZ0_uppbnd
  use mod_particle_common_test_tools, only: BE0_lowbnd,BE0_uppbnd
  use mod_pusher_tools,               only: get_orthonormals
  implicit none
  !> variables
  real*8 :: psi,U,gam,B_norm,rand
  integer,dimension(0) :: int_param
  real*8,dimension(2)  :: p_interval,costheta_interval
  real*8,dimension(3)  :: rand3,B_field,E_field
  real*8,dimension(3)  :: e1_mag_cart,e2_mag_cart,e3_mag_cart
  real*8,dimension(4)  :: real_param

  !> initialise electring and magnetic field
  call random_number(rand)
  sign_theta_sol = -1; if(rand.gt.5.d-1) sign_theta_sol=1
  call random_number(rand); phi_sol=TWOPI*rand
  call gnu_rng_interval(2,RZ0_lowbnd,RZ0_uppbnd,real_param(1:2))
  call gnu_rng_interval(2,BE0_lowbnd,BE0_uppbnd,real_param(3:4))
  call fields_sol%init_fields(0,4,int_param,real_param)
  call fields_sol%calc_EBPsiU(time_sol,i_elm_zero,st_zero,&
  phi_sol,E_field,B_field,psi,U); B_norm = norm2(B_field);
  B_field = vector_cylindrical_to_cartesian(phi_sol,B_field)
  E_field = vector_cylindrical_to_cartesian(phi_sol,E_field)
  e1_mag_cart = B_field/B_norm; 
  call get_orthonormals(e1_mag_cart,e2_mag_cart,e3_mag_cart)

  !> compute momentum and cosinus of the pitch angle interval
  p_interval = (/EThetaChi_RE_lowbnd(1),EThetaChi_RE_uppbnd(1)/)/&
  (mass*SPEED_OF_LIGHT*SPEED_OF_LIGHT)
  p_interval = p_interval*p_interval - (/1.d0,1.d0/)
  p_interval = p_interval*mass*SPEED_OF_LIGHT
  costheta_interval = (/cos(EThetaChi_RE_lowbnd(2)),cos(EThetaChi_RE_uppbnd(2))/)
  !> sampling the particle momentum
  call random_number(rand3)
  pThetaChi_sol = sample_uniform_sphere_corona_rthetaphi(p_interval*p_interval*p_interval,&
  costheta_interval,(/EThetaChi_RE_lowbnd(3),EThetaChi_RE_uppbnd(3)/),rand3)
  pThetaChi_sol = pThetaChi_sol*real(sign_theta_sol,kind=8)
  particle%p = pthetaChi_sol(1)*(e1_mag_cart*cos(pThetaChi_sol(2)) +&
  sin(pThetaChi_sol(2))*(e2_mag_cart*cos(pThetaChi_sol(3))+e3_mag_cart*sin(pThetaChi_sol(3))))  
  particle%q=q

  !> compute gc momentum
  p_gc_sol(1) = dot_product(particle%p,e1_mag_cart)
  p_gc_sol(2) = dot_product(particle%p-p_gc_sol(1)*e1_mag_cart,&
  particle%p-p_gc_sol(1)*e1_mag_cart)/(2.d0*mass*B_norm)

  !> compute cartesian orbital basis
  gam = sqrt(1.d0+(dot_product(particle%p,particle%p)/&
  (mass*SPEED_OF_LIGHT))**2.d0)
  T_cart_sol = particle%p/norm2(particle%p)
  N_cart_sol = cross_product(particle%p/(mass*gam),B_field) + &
  E_field - T_cart_sol*(dot_product(T_cart_sol,E_field))
  N_cart_sol = N_cart_sol/norm2(N_cart_sol)
  B_cart_sol = cross_product(particle%p/(mass*gam),&
  cross_product(particle%p/(mass*gam),B_field) + E_field)
  B_cart_sol = B_cart_sol/norm2(B_cart_sol)
end subroutine setup

!> clean-up kinetic realativistic test features
subroutine teardown()
  implicit none
  !> clean up particle
  p_gc_sol = 0.d0; particle%p = 0.d0; particle%q = 0;
  !> clean-up fields
  call fields_sol%deallocate_fields()
  !> clean-up orbital basis
  T_cart_sol = 0.d0; N_cart_sol = 0.d0; B_cart_sol = 0.d0;
end subroutine teardown

!> Tests ----------------------------------------------
!> test transformation of momentum from relativistic kinetic
!> to relativistic gc particles
subroutine test_momentum_relat_kinetic_to_relat_gc()
  use mod_kinetic_relativistic, only: momentum_relativistic_kinetic_to_relativistic_gc
  implicit none
  !> variables:
  real*8 :: psi,U,B_norm
  real*8,dimension(2) :: p_gc_test
  real*8,dimension(3) :: B,E
  !> compute gc momentum
  call fields_sol%calc_EBPsiU(time_sol,i_elm_zero,st_zero,&
  phi_sol,E,B,psi,U); B_norm = norm2(B); B=B/B_norm
  p_gc_test = momentum_relativistic_kinetic_to_relativistic_gc(&
  mass,particle%p,phi_sol,B_norm,B)
  !> check
  call assert_equals(p_gc_test,p_gc_sol,2,tol_real8,&
  "Error momentum relat. kinetic to relat. gc: gc momenta mismatch!")
end subroutine test_momentum_relat_kinetic_to_relat_gc

!> test the momentum coordinate transform from spherical to cartesian
subroutine test_kinetic_relat_momentum_spherical_cart()
  use mod_kinetic_relativistic, only: kinetic_relativistic_momentum_spherical_to_cart
  implicit none
  !> variables:
  real*8 :: psi,U
  real*8,dimension(3) :: p_kin_test,B,E
  !> transform kinetic particle momentum from spherical to cartesian coord.
  call fields_sol%calc_EBPsiU(time_sol,i_elm_zero,st_zero,&
  phi_sol,E,B,psi,U)
  p_kin_test = kinetic_relativistic_momentum_spherical_to_cart(&
  sign_theta_sol,phi_sol,pThetaChi_sol,B/norm2(B))
  !> checks
  call assert_equals(p_kin_test,particle%p,3,tol_real8,&
  "Error kinetic relat. momentum spherical to cartesian: momenta mismatch!")
end subroutine test_kinetic_relat_momentum_spherical_cart

subroutine test_dummy()
  use mod_kinetic_relativistic
  implicit none
  write(*,'(/A)') "particle kinetic dummy test"
end subroutine test_dummy

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

