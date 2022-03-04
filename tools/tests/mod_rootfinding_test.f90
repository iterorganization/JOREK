!> the mod_rootfinding_test module contains variables and procedures
!> used for testing the method used for finding roots of nonlinear equations
module mod_rootfinding_test
use fruit
use constants,       only: PI,TWOPI
use mod_rootfinding, only: fun_dfun_2d,fun_dfun_3d
implicit none
private
public :: run_fruit_rootfinding

!> Variables and datatypes -------------------------------------------------
type, extends(fun_dfun_2d) :: fun_cossinpar_2d
  real*8,dimension(2) :: x_min,x_max
  real*8              :: ksin,kcos
  real*8,dimension(2) :: a,b,c
  contains
  procedure,pass :: f_df => cossinpar_2d
  procedure,pass :: reset_fun_cossipar_2d
end type fun_cossinpar_2d

type, extends(fun_dfun_3d) :: fun_cossinpar
  real*8,dimension(3) :: x_min,x_max
  real*8,dimension(3) :: ksin,kcos,a,b,c
  contains
  procedure,pass :: f_df => cossinpar
  procedure,pass :: reset_fun_cossipar
end type fun_cossinpar

integer,parameter                :: maxit=100000
real*8,parameter                 :: tol_real8=5.d-10
real*8,parameter                 :: delta=1.3d1
real*8,dimension(3),parameter    :: abc_quadratic_lowbnd=(/-3.d1,-2.d1,-2.45d1/)
real*8,dimension(3),parameter    :: abc_quadratic_uppbnd=(/5.d1,3.d1,3.43d1/)
real*8,dimension(3),parameter    :: acoeffdeltabeta_lowbnd=(/-3.d1,1.d-2,2.45d-2/)
real*8,dimension(3),parameter    :: acoeffdeltabeta_uppbnd=(/5.d1,1.d1,3.43d1/)
real*8,dimension(3),parameter    :: kcos_sol=(/PI/7.d0,TWOPI,5.d0*PI/)
real*8,dimension(3),parameter    :: ksin_sol=(/PI/6.d0,2.d0*PI/3.d0,PI/)
real*8,dimension(3),parameter    :: a_parabolic=(/-5.d0,9.d0,-7.d0/)
real*8,dimension(3),parameter    :: b_parabolic=(/0.d0,-5.d0,1.2d1/)
real*8,dimension(3),parameter    :: c_parabolic=(/4.d0,-1.d1,1.d1/)
real*8,dimension(3),parameter    :: x_min_sol=(/-4.d1,3.d0,-1.2d1/)
real*8,dimension(3),parameter    :: x_max_sol=(/2.d1,5.d1,9.d0/)
type(fun_cossinpar)              :: cossinpar_sol
type(fun_cossinpar_2d)           :: cossinpar_2d_sol
integer                          :: n_no_int_coords_sol
integer,dimension(:),allocatable :: no_int_coords_sol
real*8,dimension(2)              :: y0_cossinpar_2d
real*8,dimension(3)              :: x0_rand,y0_cossinpar
real*8,dimension(3)              :: real_quadratic_coeff_real_roots
real*8,dimension(3)              :: real_cubic_coeff_real_roots
real*8,dimension(3)              :: real_cubic_coeff_cmplx_roots

!> Interfaces --------------------------------------------------------------
!> define the cossinpar constructor by overloading
interface fun_cossinpar_2d
  module procedure init_fun_cossipar_2d
end interface fun_cossinpar_2d
interface fun_cossinpar
  module procedure init_fun_cossipar
end interface fun_cossinpar

contains
!> Fruit basket ------------------------------------------------------------
!> run all set-up, tear-down and test procedures
subroutine run_fruit_rootfinding()
  implicit none
  write(*,*) "  ... setting-up: root finding tests"
  call setup
  call setup_polynomials
  write(*,*) "  ... running: root finding tests"
  call root_quadratic_real_coeff_real_roots
  call root_cubic_real_coeff_real_roots
  call root_cubic_real_coeff_cmplx_roots
  call test_newtons_method_2d_zeros
  call test_newtons_method_2d_y0
  call test_newtons_method_3d_zeros
  call test_newtons_method_3d_y0
  write(*,*) "  ... tearing-down: root finding tests"
  call teardown
end subroutine run_fruit_rootfinding

!> Set-up and tear-down ----------------------------------------------------
!> initialise test features
subroutine setup()
  use mod_gnu_rng, only: gnu_rng_interval
  implicit none
  !> variables:
  real*8,dimension(3)   :: x_loc
  real*8,dimension(3,3) :: J_loc
  !> initialise no integer coordinates
  n_no_int_coords_sol=0; allocate(no_int_coords_sol(n_no_int_coords_sol));
  !> initialise the cossinpar object
  cossinpar_2d_sol = fun_cossinpar_2d(x_min_sol(1:2),x_max_sol(1:2),&
  kcos_sol(1),ksin_sol(1),a_parabolic(1:2),b_parabolic(1:2),c_parabolic(1:2))
  cossinpar_sol = fun_cossinpar(x_min_sol,x_max_sol,kcos_sol,ksin_sol,&
  a_parabolic,b_parabolic,c_parabolic)
  !> set-up a random initial value
  call gnu_rng_interval(3,x_min_sol,x_max_sol,x0_rand)
  !> compute y0
  call gnu_rng_interval(3,x_min_sol,x_max_sol,x_loc)
  call cossinpar_sol%f_df(y0_cossinpar,J_loc,x_loc,&
  n_no_int_coords_sol,no_int_coords_sol)
  call cossinpar_2d_sol%f_df(y0_cossinpar_2d,J_loc,x_loc(1:2),&
  n_no_int_coords_sol,no_int_coords_sol)
end subroutine setup

!> set-up features for finding the root of polynomials
subroutine setup_polynomials()
  use mod_gnu_rng, only: gnu_rng_interval
  implicit none
  real*8 :: Q,R,sig,rand,uppbnd
  real*8,dimension(2) :: rands
  real*8,dimension(3) :: adeltabeta
  !> compute the quadratic coefficients:
  call gnu_rng_interval(3,abc_quadratic_lowbnd,abc_quadratic_uppbnd,&
  real_quadratic_coeff_real_roots)
  !> be sure to have real roots
  if(real_quadratic_coeff_real_roots(1)*real_quadratic_coeff_real_roots(3).gt.0.d0) then
    if(real_quadratic_coeff_real_roots(2)*real_quadratic_coeff_real_roots(2).lt.&
    4.d0*real_quadratic_coeff_real_roots(1)*real_quadratic_coeff_real_roots(3)) then
      call gnu_rng_interval(3,acoeffdeltabeta_lowbnd,acoeffdeltabeta_uppbnd,adeltabeta)
      real_quadratic_coeff_real_roots(2) = sign(1.d0,real_quadratic_coeff_real_roots(2))*(4.d0*&
      real_quadratic_coeff_real_roots(1)*real_quadratic_coeff_real_roots(3) + adeltabeta(2))       
    endif
  endif
  !> compute the coefficients for purely real solutions
  call gnu_rng_interval(3,acoeffdeltabeta_lowbnd,&
  acoeffdeltabeta_uppbnd,adeltabeta)
  real_cubic_coeff_real_roots(1)  = adeltabeta(1)
  !> compute first and 0th order coefficients for real solutions
  call random_number(rands)
  sig = -1.d0; if(rands(1).gt.5.d-1) sig=1.d0;
  uppbnd = (adeltabeta(1)*adeltabeta(1)-9.d0*adeltabeta(2))/3.d0
  real_cubic_coeff_real_roots(2) = uppbnd - adeltabeta(3)*rands(2)
  Q = (adeltabeta(1)*adeltabeta(1)-3.d0*real_cubic_coeff_real_roots(2))/9.d0
  real_cubic_coeff_real_roots(3) = (9.d0*adeltabeta(1)*real_cubic_coeff_real_roots(2)-&
  2.d0*adeltabeta(1)*adeltabeta(1)*adeltabeta(1)+&
  sig*5.4d1*sqrt((Q-adeltabeta(2))**3.d0))/2.7d1
  R = (2.d0*real_cubic_coeff_real_roots(1)**3.d0-&
  9.d0*real_cubic_coeff_real_roots(1)*real_cubic_coeff_real_roots(2)+&
  2.7d1*real_cubic_coeff_real_roots(3))/5.4d1
  call assert_true(((Q-adeltabeta(2)).gt.0.d0).and.((R**2).lt.(Q**3.d0)),&
  "Error set-up cube root finding real: Q and R conditions for real roots not respected!")
  !> compute the coefficients for complex solutions
  call gnu_rng_interval(3,acoeffdeltabeta_lowbnd,&
  acoeffdeltabeta_uppbnd,adeltabeta)
  real_cubic_coeff_cmplx_roots(1) = adeltabeta(1)
  !> compute first and 0th order coefficients for complex solutions
  call random_number(rands)
  sig = -1.d0; if(rands(1).gt.5.d-1) sig=1.d0; 
  uppbnd = (adeltabeta(1)*adeltabeta(1)+9.d0*adeltabeta(2))/3.d0
  real_cubic_coeff_cmplx_roots(2) = uppbnd - adeltabeta(3)*rands(2)
  Q = (adeltabeta(1)*adeltabeta(1)-3.d0*real_cubic_coeff_cmplx_roots(2))/9.d0
  real_cubic_coeff_cmplx_roots(3) = (9.d0*adeltabeta(1)*real_cubic_coeff_cmplx_roots(2)-&
  2.d0*adeltabeta(1)*adeltabeta(1)*adeltabeta(1)+&
  sig*5.4d1*sqrt((Q+adeltabeta(2))**3.d0))/2.7d1
  R = (2.d0*real_cubic_coeff_cmplx_roots(1)**3.d0-&
  9.d0*real_cubic_coeff_cmplx_roots(1)*real_cubic_coeff_cmplx_roots(2)+&
  2.7d1*real_cubic_coeff_cmplx_roots(3))/5.4d1
  call assert_true(((Q+adeltabeta(2)).gt.0.d0).and.((R**2).ge.(Q**3.d0)),&
  "Error set-up cube root finding real: Q and R conditions for complex roots not respected!")
end subroutine setup_polynomials

!> tearing-down test features
subroutine teardown()
  implicit none
  x0_rand=0.d0; y0_cossinpar=0.d0; n_no_int_coords_sol=0;
  real_cubic_coeff_real_roots=0.d0; real_cubic_coeff_cmplx_roots=0.d0;
  if(allocated(no_int_coords_sol)) deallocate(no_int_coords_sol)
  call cossinpar_sol%reset_fun_cossipar
end subroutine teardown

!> Tests -------------------------------------------------------------------
!> test find at least one zero
subroutine test_newtons_method_2d_zeros()
  use mod_rootfinding, only: newtons_method
  implicit none
  !> variables
  integer               :: ierr
  real*8,dimension(2)   :: x_test,y_test
  real*8,dimension(2,2) :: J_test
  !> call the newton method:
  call newtons_method(cossinpar_2d_sol,(/0.d0,0.d0/),x0_rand(1:2),x_test,&
  n_no_int_coords_sol,no_int_coords_sol,ierr,tol_real8,maxit)
  !> checks
  call cossinpar_2d_sol%f_df(y_test,J_test,x_test,&
  n_no_int_coords_sol,no_int_coords_sol)
  call assert_equals(y_test,(/0.d0,0.d0/),2,tol_real8,&
  "Error root finding newton method 2d object: root values are not zeros!")
  call assert_equals(ierr,0,&
  "Error root finding newton method 2d object: maximum number of iterations reach finding zeros!")
end subroutine test_newtons_method_2d_zeros

!> test find at least one zero
subroutine test_newtons_method_2d_y0()
  use mod_rootfinding, only: newtons_method
  implicit none
  !> variables
  integer               :: ierr
  real*8,dimension(2)   :: x_test,y_test
  real*8,dimension(2,2) :: J_test
  !> call the newton method:
  call newtons_method(cossinpar_2d_sol,y0_cossinpar_2d,x0_rand(1:2),x_test,&
  n_no_int_coords_sol,no_int_coords_sol,ierr,tol_real8,maxit)
  !> checks
  call cossinpar_2d_sol%f_df(y_test,J_test,x_test,&
  n_no_int_coords_sol,no_int_coords_sol)
  call assert_equals(y_test,y0_cossinpar_2d,2,tol_real8,&
  "Error root finding newton method 2d object: root values are not zeros!")
  call assert_equals(ierr,0,&
  "Error root finding newton method 2d object: maximum number of iterations reach finding zeros!")
end subroutine test_newtons_method_2d_y0

!> test find at least one zero
subroutine test_newtons_method_3d_zeros()
  use mod_rootfinding, only: newtons_method
  implicit none
  !> variables
  integer               :: ierr
  real*8,dimension(3)   :: x_test,y_test
  real*8,dimension(3,3) :: J_test
  !> call the newton method:
  call newtons_method(cossinpar_sol,(/0.d0,0.d0,0.d0/),x0_rand,x_test,&
  n_no_int_coords_sol,no_int_coords_sol,ierr,tol_real8,maxit)
  !> checks
  call cossinpar_sol%f_df(y_test,J_test,x_test,&
  n_no_int_coords_sol,no_int_coords_sol)
  call assert_equals(y_test,(/0.d0,0.d0,0.d0/),3,tol_real8,&
  "Error root finding newton method 3d object: root values are not zeros!")
  call assert_equals(ierr,0,&
  "Error root finding newton method 3d object: maximum number of iterations reach finding zeros!")
end subroutine test_newtons_method_3d_zeros

!> test find at least one zero
subroutine test_newtons_method_3d_y0()
  use mod_rootfinding, only: newtons_method
  implicit none
  !> variables
  integer               :: ierr
  real*8,dimension(3)   :: x_test,y_test
  real*8,dimension(3,3) :: J_test
  !> call the newton method:
  call newtons_method(cossinpar_sol,y0_cossinpar,x0_rand,x_test,&
  n_no_int_coords_sol,no_int_coords_sol,ierr,tol_real8,maxit)
  !> checks
  call cossinpar_sol%f_df(y_test,J_test,x_test,&
  n_no_int_coords_sol,no_int_coords_sol)
  call assert_equals(y_test,y0_cossinpar,3,tol_real8,&
  "Error root finding newton method 3d object: root values are not zeros!")
  call assert_equals(ierr,0,&
  "Error root finding newton method 3d object: maximum number of iterations reach finding zeros!")
end subroutine test_newtons_method_3d_y0

!> test the finding of roots for quadratic polynomials (real roots only real coeff. only)
subroutine root_quadratic_real_coeff_real_roots()
  use mod_rootfinding, only: root_quadratic
  implicit none
  real*8,dimension(2) :: y,x,zeros
  !> initialisation
  zeros = 0.d0
  !> compute the roots
  x = root_quadratic(real_quadratic_coeff_real_roots)
  y = real_quadratic_coeff_real_roots(1)*x**2.d0 + &
  real_quadratic_coeff_real_roots(2)*x+real_quadratic_coeff_real_roots(3)
  call assert_equals(y,zeros,2,tol_real8,&
  "Error quadratic polynomials root real coeff real solutions: x values are not roots!")
end subroutine root_quadratic_real_coeff_real_roots

!> test the finding of roots for cubic polynomials (real roots only real coeff. only)
subroutine root_cubic_real_coeff_real_roots()
  use mod_rootfinding, only: root_cubic
  implicit none
  complex*16,dimension(3) :: y,x,zeros
  !> initialisation
  zeros = cmplx(0.d0,0.d0)
  !> compute the roots
  x = root_cubic(real_cubic_coeff_real_roots)
  y = x**3.d0 + real_cubic_coeff_real_roots(1)*x**2.d0 + &
  real_cubic_coeff_real_roots(2)*x+real_cubic_coeff_real_roots(3)
  call assert_equals(aimag(x),(/0.d0,0.d0,0.d0/),3,&
  "Error cubic polynomials root real coeff real solutions: x roots not real!")
  call assert_equals(y,zeros,3,tol_real8,&
  "Error cubic polynomials root real coeff real solutions: x values are not roots!")
end subroutine root_cubic_real_coeff_real_roots

!> test the finding of roots for cubic polynomials (complex roots only real coeff. only)
subroutine root_cubic_real_coeff_cmplx_roots()
  use mod_rootfinding, only: root_cubic
  implicit none
  complex*16,dimension(3) :: y,x,zeros
  !> initialisation
  zeros = cmplx(0.d0,0.d0)
  !> compute the roots
  x = root_cubic(real_cubic_coeff_cmplx_roots)
  y = x**3.d0 + real_cubic_coeff_cmplx_roots(1)*x**2.d0 + &
  real_cubic_coeff_cmplx_roots(2)*x+real_cubic_coeff_cmplx_roots(3)
  call assert_true(.not.all(aimag(x).ne.0.d0),&
  "Error cubic polynomials root real coeff complex solutions: x roots not imaginary!")
  call assert_equals(y,zeros,3,tol_real8,&
  "Error cubic polynomials root real coeff complex solutions: x values are not roots!")
end subroutine root_cubic_real_coeff_cmplx_roots

!> Tools -------------------------------------------------------------------
!> initialise function cossinparabola to use for testing
!> inputs:
!>   this:  (fun_cossinpar) function object to be initialised
!>   x_min: (real8)(3) left boundary of the space interval x,y,z
!>   x_max: (real8)(3) right boundary of the space interval x,y,z
!>   ksin   (real8)(3) cosinus mode numbers
!>   kcos   (real8)(3) sinus mode numbers
!>   a:     (real8)(3) parabola x^2 parameters
!>   b:     (real8)(3) parabola x b parameters
!>   c:     (real8)(3) parabola c parameters
!> outputs:
!>   this:  (fun_cossinpar) initialised function object
function init_fun_cossipar(x_min,x_max,ksin,kcos,a,b,c) result(this)
  implicit none
  !> inputs-outputs:
  type(fun_cossinpar) :: this
  !> inputs:
  real*8,dimension(3),intent(in) :: kcos,ksin,a,b,c
  real*8,dimension(3),intent(in) :: x_min,x_max
  !> initialise variables
  this%x_min=x_min; this%x_max=x_max; this%ksin=ksin; 
  this%kcos=kcos; this%a=a; this%b=b; this%c=c;
end function init_fun_cossipar

!> reset function cossinparabola parameters
!> inputs:
!>   this:  (fun_cossinpar) function object to reset
!> outputs:
!>   this:  (fun_cossinpar) reset function object
subroutine reset_fun_cossipar(this)
  implicit none
  !> inputs-outputs:
  class(fun_cossinpar),intent(inout) :: this
  !> reset variables
  this%x_min=0.d0; this%x_max=0.d0; this%ksin=0.d0; 
  this%kcos=0.d0; this%a=0.d0; this%b=0.d0; this%c=0.d0;
end subroutine reset_fun_cossipar

!> function cossinparabola to be used for testing
!> inputs:
!>   this: (cossinpar) cossinpar object
!>   x:    (real8)(3) abscissa
!> outputs:
!>   y:    (real8)(3) ordinate
!>   J:    (real8)(3,3) jacobian
pure subroutine cossinpar(this,f,J,x,n_int_coords,int_coords)
  implicit none
  !> inputs:
  class(fun_cossinpar),intent(inout) :: this
  integer,intent(in)                 :: n_int_coords
  integer,dimension(n_int_coords),intent(in) :: int_coords
  real*8,dimension(3),intent(in)  :: x
  !> outputs:
  real*8,dimension(3),intent(out)   :: f
  real*8,dimension(3,3),intent(out) :: J
  !> variables
  real*8,dimension(3) :: x_norm_val
  !> initialise
  x_norm_val = (x - this%x_min)/(this%x_max - this%x_min)
  !> compute function
  f = (/cos(this%kcos(1)*x_norm_val(1))*sin(this%ksin(1)*&
      x_norm_val(2))*(this%a(1)*x_norm_val(3)*x_norm_val(3)+&
      this%b(1)*x_norm_val(3)+this%c(1)),sin(this%ksin(2)*&
      x_norm_val(1))*(this%a(2)*x_norm_val(2)*x_norm_val(2)+&
      this%b(2)*x_norm_val(2)+this%c(2))*cos(this%kcos(2)*x_norm_val(3)),&
      (this%a(3)*x_norm_val(1)*x_norm_val(1)+this%b(3)*x_norm_val(1)+&
      this%c(3))*sin(this%ksin(3)*x_norm_val(2))*cos(this%kcos(3)*x_norm_val(3))/)
  !> comput jacobian
  J(:,1) = (/-this%kcos(1)*sin(this%kcos(1)*x_norm_val(1))*&
           sin(this%ksin(1)*x_norm_val(2))*(this%a(1)*x_norm_val(3)*&
           x_norm_val(3)+this%b(1)*x_norm_val(3)+this%c(1)),&
           this%ksin(2)*cos(this%ksin(2)*x_norm_val(1))*(this%a(2)*&
           x_norm_val(2)*x_norm_val(2)+this%b(2)*x_norm_val(2)+&
           this%c(2))*cos(this%kcos(2)*x_norm_val(3)),(2.d0*this%a(3)*&
           x_norm_val(1)+this%b(1))*sin(this%ksin(3)*x_norm_val(2))*&
           cos(this%kcos(3)*x_norm_val(3))/)/(this%x_max(1) - this%x_min(1))
  J(:,2) = (/cos(this%kcos(1)*x_norm_val(1))*this%ksin(1)*cos(this%ksin(1)*&
           x_norm_val(2))*(this%a(1)*x_norm_val(3)*x_norm_val(3)+this%b(1)&
           *x_norm_val(3)+this%c(1)),sin(this%ksin(2)*x_norm_val(1))*&
           (2.d0*this%a(2)*x_norm_val(2)+this%b(2))*cos(this%kcos(2)*x_norm_val(3)),&
           (this%a(3)*x_norm_val(1)*x_norm_val(1)+this%b(3)*x_norm_val(1)+&
           this%c(3))*this%ksin(3)*cos(this%ksin(3)*x_norm_val(2))*&
           cos(this%kcos(3)*x_norm_val(3))/)/(this%x_max(2) - this%x_min(2))
  J(:,3) = (/cos(this%kcos(1)*x_norm_val(1))*sin(this%ksin(1)*x_norm_val(2))*&
           (2.d0*this%a(1)*x_norm_val(3)+this%b(1)),-sin(this%ksin(2)*&
           x_norm_val(1))*(this%a(2)*x_norm_val(2)*x_norm_val(2)+this%b(2)*&
           x_norm_val(2)+this%c(2))*this%kcos(2)*sin(this%kcos(2)*x_norm_val(3)),&
           -(this%a(3)*x_norm_val(1)*x_norm_val(1)+this%b(3)*x_norm_val(1)+&
           this%c(3))*sin(this%ksin(3)*x_norm_val(2))*this%kcos(3)*sin(this%kcos(3)*&
           x_norm_val(3))/)/(this%x_max(3) - this%x_min(3))
end subroutine cossinpar

!> initialise function cossinparabola 2d to use for testing
!> inputs:
!>   this:  (fun_cossinpar_2d) function object to be initialised
!>   x_min: (real8)(2) left boundary of the space interval x,y
!>   x_max: (real8)(2) right boundary of the space interval x,y
!>   ksin   (real8) cosinus mode numbers
!>   kcos   (real8) sinus mode numbers
!>   a:     (real8)(2) parabola x^2 parameters
!>   b:     (real8)(2) parabola x b parameters
!>   c:     (real8)(2) parabola c parameters
!> outputs:
!>   this:  (fun_cossinpar_2d) initialised function object
function init_fun_cossipar_2d(x_min,x_max,ksin,kcos,a,b,c) result(this)
  implicit none
  !> inputs-outputs:
  type(fun_cossinpar_2d) :: this
  !> inputs:
  real*8,intent(in)              :: kcos,ksin
  real*8,dimension(2),intent(in) :: x_min,x_max,a,b,c
  !> initialise variables
  this%x_min=x_min; this%x_max=x_max; this%ksin=ksin; 
  this%kcos=kcos; this%a=a; this%b=b; this%c=c;
end function init_fun_cossipar_2d

!> reset function cossinparabola parameters
!> inputs:
!>   this:  (fun_cossinpar_2d) function object to reset
!> outputs:
!>   this:  (fun_cossinpar_2d) reset function object
subroutine reset_fun_cossipar_2d(this)
  implicit none
  !> inputs-outputs:
  class(fun_cossinpar_2d),intent(inout) :: this
  !> reset variables
  this%x_min=0.d0; this%x_max=0.d0; this%ksin=0.d0; 
  this%kcos=0.d0; this%a=0.d0; this%b=0.d0; this%c=0.d0;
end subroutine reset_fun_cossipar_2d

!> function cossinparabola to be used for testing (2d)
!> inputs:
!>   this: (cossinpar_2d) cossinpar object
!>   x:    (real8)(2) abscissa
!> outputs:
!>   y:    (real8)(2) ordinate
!>   J:    (real8)(2,2) jacobian
pure subroutine cossinpar_2d(this,f,J,x,n_int_coords,int_coords)
  implicit none
  !> inputs:
  class(fun_cossinpar_2d),intent(inout) :: this
  integer,intent(in)                    :: n_int_coords
  integer,dimension(n_int_coords),intent(in) :: int_coords
  real*8,dimension(2),intent(in)  :: x
  !> outputs:
  real*8,dimension(2),intent(out)   :: f
  real*8,dimension(2,2),intent(out) :: J
  !> variables
  real*8,dimension(2) :: x_norm_val
  !> initialise
  x_norm_val = (x - this%x_min)/(this%x_max - this%x_min)
  !> compute function
  f = (/cos(this%kcos*x_norm_val(1))*(this%a(1)*x_norm_val(2)*&
      x_norm_val(2)+this%b(1)*x_norm_val(2)+this%c(1)),&
      (this%a(2)*x_norm_val(1)*x_norm_val(1)+this%b(2)*&
      x_norm_val(1)+this%c(2))*sin(this%ksin*x_norm_val(2))/)
  !> comput jacobian
  J(:,1) = (/-1.d0*this%kcos*sin(this%kcos*x_norm_val(1))*(this%a(1)*x_norm_val(2)*&
           x_norm_val(2)+this%b(1)*x_norm_val(2)+this%c(1)),&
           (2.d0*this%a(2)*x_norm_val(1)+this%b(2))*sin(this%ksin*x_norm_val(2))/)/&
           (this%x_max(1) - this%x_min(1))
  J(:,2) = (/cos(this%kcos*x_norm_val(1))*(2.d0*this%a(1)*x_norm_val(2)+this%b(1)),&
           (this%a(2)*x_norm_val(1)*x_norm_val(1)+this%b(2)*x_norm_val(1)+this%c(2))*&
           this%ksin*cos(this%ksin*x_norm_val(2))/)/(this%x_max(2) - this%x_min(2))
end subroutine cossinpar_2d

!>--------------------------------------------------------------------------
end module mod_rootfinding_test

