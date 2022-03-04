!> module mod_math_operators_test contains variables and
!> and procedures for testing the procedures in 
!> mod_math_operators
module mod_math_operators_test
use fruit
implicit none

private
public :: run_fruit_math_operators

!> Variables --------------------------------------------
integer,parameter :: zero_r4=real(0.d0,kind=4)
integer,parameter :: n_vectors=4
real*4,parameter  :: tol_r4=real(1.0d-1,kind=4)
real*8,parameter  :: tol_r8=2.5d-10
real*8,parameter  :: tol_c16=7.50d-9
!> vectors for testing the vector product
real*8,dimension(2) :: a_interval_r8=(/-3.5d1,5.4d1/)
real*8,dimension(2) :: b_interval_r8=(/5.4d1,1.15d2/)
real*4,dimension(2,n_vectors) :: vec_a_2_r4
real*8,dimension(2,n_vectors) :: vec_a_2_r8
real*4,dimension(3,n_vectors) :: vec_a_3_r4,vec_b_3_r4
real*8,dimension(3,n_vectors) :: vec_a_3_r8,vec_b_3_r8
real*4,dimension(2,2,n_vectors) :: matrix_2x2_r4
real*8,dimension(2,2,n_vectors) :: matrix_2x2_r8
real*4,dimension(3,3,n_vectors) :: matrix_3x3_r4
real*8,dimension(3,3,n_vectors) :: matrix_3x3_r8

contains

!> Test basket ------------------------------------------
!> test basket for setting-up, executing and tearing-down
!> all test feature
subroutine run_fruit_math_operators()
  implicit none
  write(*,'(/A)') "  ... setting-up: math operators tests"
  call setup
  write(*,'(/A)') "  ... running: math operators tests"
  call test_cross_product
  call test_solve_2x2_linear_problem
  call test_solve_3x3_linear_problem
  call test_compute_eigenvalues_real_2x2_r8
  call test_compute_eigenvalues_real_3x3_r8
  write(*,'(/A)') "  ... tearing-down: math operators tests"
  call teardown
end subroutine run_fruit_math_operators

!> Set-up and tear-down ---------------------------------
!> Set-up of the math operator test features
subroutine setup()
  use mod_gnu_rng, only: gnu_rng_interval
  implicit none
  !> generate vectors for cross product tests
  call gnu_rng_interval(2,n_vectors,a_interval_r8,vec_a_2_r8)
  vec_a_2_r4 = real(vec_a_2_r8,kind=4)
  call gnu_rng_interval(3,n_vectors,a_interval_r8,vec_a_3_r8)
  call gnu_rng_interval(3,n_vectors,b_interval_r8,vec_b_3_r8)
  vec_a_3_r4 = real(vec_a_3_r8,kind=4)
  vec_b_3_r4 = real(vec_b_3_r8,kind=4)
  call gnu_rng_interval(2,2,n_vectors,b_interval_r8,matrix_2x2_r8)
  matrix_2x2_r4 = real(matrix_2x2_r8,kind=4)
  call gnu_rng_interval(3,3,n_vectors,b_interval_r8,matrix_3x3_r8)
  matrix_3x3_r4 = real(matrix_3x3_r8,kind=4)
end subroutine setup

!> tear-down the math operator test features
subroutine teardown()
  implicit none
  !> set all static vectors to zero
  vec_a_2_r8 = 0.d0;  vec_a_2_r4 = zero_r4;
  vec_a_3_r4 = zero_r4; vec_b_3_r4 = zero_r4;
  vec_a_3_r8 = 0.d0; vec_b_3_r8 = 0.d0;
  matrix_2x2_r8 = 0.d0; matrix_2x2_r4 = zero_r4;
  matrix_3x3_r8 = 0.d0; matrix_3x3_r4 = zero_r4;
end subroutine teardown

!> Tests ------------------------------------------------
!> Test vector product for both single and double precision
subroutine test_cross_product()
  use mod_math_operators, only: cross_product
  implicit none
  integer :: ii
  real*4,dimension(3) :: vec_c_r4
  real*8,dimension(3) :: vec_c_r8

  !> test single precision cross product
  do ii=1,n_vectors
    vec_c_r4 = cross_product(vec_a_3_r4(:,ii),vec_b_3_r4(:,ii))
    call assert_true((abs(dot_product(vec_a_3_r4(:,ii),vec_c_r4)).lt.tol_r4),&
    "Error math operators cross product (float): vectors a and c not orthogonal")
    call assert_true((abs(dot_product(vec_b_3_r4(:,ii),vec_c_r4)).lt.tol_r4),&
    "Error math operators cross product (float): vectors b and c not orthogonal")
  enddo
  do ii=1,n_vectors
    vec_c_r8 = cross_product(vec_a_3_r8(:,ii),vec_b_3_r8(:,ii))
    call assert_equals(dot_product(vec_a_3_r8(:,ii),vec_c_r8),0.d0,tol_r8,&
    "Error math operators cross product (double): vectors a and c not orthogonal")
    call assert_equals(dot_product(vec_b_3_r8(:,ii),vec_c_r8),0.d0,tol_r8,&
    "Error math operators cross product (double): vectors b and c not orthogonal")
  enddo
end subroutine test_cross_product

!> test 2x2 linear problem solver for both single and double precision
subroutine test_solve_2x2_linear_problem()
  use mod_math_operators, only: solve_2x2_linear_problem
  implicit none
  integer :: ii
  real*4,dimension(2) :: x_r4
  real*8,dimension(2) :: x_r8
  !> test single precision 2x2 linear solver
  do ii=1,n_vectors
    call solve_2x2_linear_problem(matrix_2x2_r4(:,:,ii),vec_a_2_r4(:,ii),x_r4)
    call assert_equals(matmul(matrix_2x2_r4(:,:,ii),x_r4),&
    vec_a_2_r4(:,ii),2,tol_r4,&
    "Error math operators solve 2x2 linear problems (float): rhs mismatch!") 
  enddo
  do ii=1,n_vectors
    call solve_2x2_linear_problem(matrix_2x2_r8(:,:,ii),vec_a_2_r8(:,ii),x_r8)
    call assert_equals(matmul(matrix_2x2_r8(:,:,ii),x_r8),&
    vec_a_2_r8(:,ii),2,tol_r8,&
    "Error math operators solve 2x2 linear problems (double): rhs mismatch!") 
  enddo
end subroutine test_solve_2x2_linear_problem

!> test 3x3 linear problem solver for both single and double precision
subroutine test_solve_3x3_linear_problem()
  use mod_math_operators, only: solve_3x3_linear_problem
  implicit none
  integer :: ii
  real*4,dimension(3) :: x_r4
  real*8,dimension(3) :: x_r8
  !> test single precision 2x2 linear solver
  do ii=1,n_vectors
    call solve_3x3_linear_problem(matrix_3x3_r4(:,:,ii),vec_a_3_r4(:,ii),x_r4)
    call assert_equals(matmul(matrix_3x3_r4(:,:,ii),x_r4),&
    vec_a_3_r4(:,ii),2,tol_r4,&
    "Error math operators solve 3x3 linear problems (float): rhs mismatch!") 
  enddo
  do ii=1,n_vectors
    call solve_3x3_linear_problem(matrix_3x3_r8(:,:,ii),vec_a_3_r8(:,ii),x_r8)
    call assert_equals(matmul(matrix_3x3_r8(:,:,ii),x_r8),&
    vec_a_3_r8(:,ii),3,tol_r8,&
    "Error math operators solve 3x3 linear problems (double): rhs mismatch!") 
  enddo
end subroutine test_solve_3x3_linear_problem

!> test method for computing the eigenvalues of a 2x2 double matrix
!> with real values
subroutine test_compute_eigenvalues_real_2x2_r8()
  use mod_math_operators, only: compute_eigenvalues
  implicit none
  !> variables:
  integer :: ii,jj
  real*8,dimension(2)           :: eigv
  real*8,dimension(2,2)         :: A
  real*8,dimension(2,n_vectors) :: det_test,zeros
  !> initialisations
  zeros = 0.d0; det_test = zeros;
  do jj=1,n_vectors
    !> compute eigenvalues
    A = matrix_2x2_r8(:,:,jj)
    call compute_eigenvalues(A,eigv)
    !> compute determinant of A-eigv(ii)*I 
    do ii=1,size(eigv)
      det_test(ii,jj) = (A(1,1)-eigv(ii))*(A(2,2)-eigv(ii))-A(2,1)*A(1,2) 
    enddo
    !> checks
    call assert_equals(det_test,zeros,2,n_vectors,tol_c16,&
    "Error math operators compute eigenvalues real 2x2 r8: determinants not zero!")
  enddo
end subroutine test_compute_eigenvalues_real_2x2_r8

!> test method for computing the eigenvalues of a 3x3 double matrix
!> with real values
subroutine test_compute_eigenvalues_real_3x3_r8()
  use mod_math_operators, only: compute_eigenvalues
  implicit none
  !> variables:
  integer :: ii,jj
  real*8,dimension(3,3)             :: zeros_r8
  complex*16,dimension(3)           :: eigv
  complex*16,dimension(3,3)         :: A
  complex*16,dimension(3,n_vectors) :: det_test,zeros
  !> initialisations
  zeros = cmplx(0.d0,0.d0); zeros_r8 = 0.d0; det_test = zeros;
  do jj=1,n_vectors
    !> compute eigenvalues
    call compute_eigenvalues(matrix_3x3_r8(:,:,jj),eigv)
    !> compute determinant of A-eigv(ii)*I
    A = cmplx(matrix_3x3_r8(:,:,jj),zeros_r8)
    do ii=1,size(eigv)
      det_test(ii,jj) = (A(1,1)-eigv(ii))*((A(2,2)-eigv(ii))*(A(3,3)-eigv(ii))-A(3,2)*A(2,3)) + & 
                     A(1,2)*(A(3,1)*A(2,3)-A(2,1)*(A(3,3)-eigv(ii))) + &
                     A(1,3)*(A(2,1)*A(3,2)-A(3,1)*(A(2,2)-eigv(ii)))
    enddo
    !> checks
    call assert_equals(det_test,zeros,3,n_vectors,tol_c16,&
    "Error math operators compute eigenvalues real 3x3 r8: determinants not zero!")
  enddo
end subroutine test_compute_eigenvalues_real_3x3_r8

!>-------------------------------------------------------

end module mod_math_operators_test
