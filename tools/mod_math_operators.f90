!> This module contains multi-purposes mathematical operators
module mod_math_operators
  implicit none

  private
  public :: cross_product
  public :: solve_2x2_linear_problem
  public :: solve_3x3_linear_problem
  public :: compute_eigenvalues
interface cross_product
  module procedure cross_product_r4
  module procedure cross_product_r8
end interface cross_product

interface solve_2x2_linear_problem
  module procedure solve_2x2_linear_problem_r4
  module procedure solve_2x2_linear_problem_r8
end interface solve_2x2_linear_problem

interface solve_3x3_linear_problem
  module procedure solve_3x3_linear_problem_r4
  module procedure solve_3x3_linear_problem_r8
end interface solve_3x3_linear_problem

interface compute_eigenvalues
module procedure compute_eigenvalues_real_3x3_r8
end interface compute_eigenvalues

contains

!> The vector cross product single precision
pure function cross_product_r4(a, b)
  real*4, dimension(3) :: cross_product_r4
  real*4, dimension(3), intent(in) :: a, b

  cross_product_r4(1) = a(2) * b(3) - a(3) * b(2)
  cross_product_r4(2) = a(3) * b(1) - a(1) * b(3)
  cross_product_r4(3) = a(1) * b(2) - a(2) * b(1)

end function cross_product_r4
  
!> The vector cross product double precision
pure function cross_product_r8(a, b)
  real*8, dimension(3) :: cross_product_r8
  real*8, dimension(3), intent(in) :: a, b

  cross_product_r8(1) = a(2) * b(3) - a(3) * b(2)
  cross_product_r8(2) = a(3) * b(1) - a(1) * b(3)
  cross_product_r8(3) = a(1) * b(2) - a(2) * b(1)

end function cross_product_r8

!> solve a 2x2 linear problem analytically (single precision)
subroutine solve_2x2_linear_problem_r4(A,b,x)
  implicit none
  real*4,dimension(2),intent(in)   :: b
  real*4,dimension(2,2),intent(in) :: A
  real*4,dimension(2),intent(out)  :: x

  x = (/A(2,2)*b(1)-A(1,2)*b(2),&
      A(1,1)*b(2)-A(2,1)*b(1)/)/&
      (A(1,1)*A(2,2)-A(2,1)*A(1,2))

end subroutine solve_2x2_linear_problem_r4

!> solve a 2x2 linear problem analytically (double precision)
subroutine solve_2x2_linear_problem_r8(A,b,x)
  implicit none
  real*8,dimension(2),intent(in)   :: b
  real*8,dimension(2,2),intent(in) :: A
  real*8,dimension(2),intent(out)  :: x

  x = (/A(2,2)*b(1)-A(1,2)*b(2),&
      A(1,1)*b(2)-A(2,1)*b(1)/)/&
      (A(1,1)*A(2,2)-A(2,1)*A(1,2))

end subroutine solve_2x2_linear_problem_r8

!> solve a 3x3 linear problem analytically (single precision)
subroutine solve_3x3_linear_problem_r4(A,b,x)
  implicit none
  real*4,dimension(3),intent(in)   :: b
  real*4,dimension(3,3),intent(in) :: A
  real*4,dimension(3),intent(out)  :: x
  real*4 :: det
  real*4,dimension(3,3) :: invA

  det = A(1,1)*(A(2,2)*A(3,3)-A(3,2)*A(2,3)) + &
        A(1,2)*(A(3,1)*A(2,3)-A(2,1)*A(3,3)) + &
        A(1,3)*(A(2,1)*A(3,2)-A(3,1)*A(2,2))
  invA(:,1) = (/A(2,2)*A(3,3)-A(3,2)*A(2,3),&
              A(3,1)*A(2,3)-A(2,1)*A(3,3),&
              A(2,1)*A(3,2)-A(3,1)*A(2,2)/)
  invA(:,2) = (/A(3,2)*A(1,3)-A(1,2)*A(3,3),&
              A(1,1)*A(3,3)-A(3,1)*A(1,3),&
              A(3,1)*A(1,2)-A(1,1)*A(3,2)/)
  invA(:,3) = (/A(1,2)*A(2,3)-A(2,2)*A(1,3),&
              A(2,1)*A(1,3)-A(1,1)*A(2,3),&
              A(1,1)*A(2,2)-A(2,1)*A(1,2)/)
  x = (/invA(1,1)*b(1)+invA(1,2)*b(2)+invA(1,3)*b(3),&
      invA(2,1)*b(1)+invA(2,2)*b(2)+invA(2,3)*b(3),&
      invA(3,1)*b(1)+invA(3,2)*b(2)+invA(3,3)*b(3)/)/det
  
end subroutine solve_3x3_linear_problem_r4
  
!> solve a 3x3 linear problem analytically (double precision)
subroutine solve_3x3_linear_problem_r8(A,b,x)
  implicit none
  real*8,dimension(3),intent(in)   :: b
  real*8,dimension(3,3),intent(in) :: A
  real*8,dimension(3),intent(out)  :: x
  real*8 :: det
  real*8,dimension(3,3) :: invA

  det = A(1,1)*(A(2,2)*A(3,3)-A(3,2)*A(2,3)) + &
        A(1,2)*(A(3,1)*A(2,3)-A(2,1)*A(3,3)) + &
        A(1,3)*(A(2,1)*A(3,2)-A(3,1)*A(2,2))
  invA(:,1) = (/A(2,2)*A(3,3)-A(3,2)*A(2,3),&
              A(3,1)*A(2,3)-A(2,1)*A(3,3),&
              A(2,1)*A(3,2)-A(3,1)*A(2,2)/)
  invA(:,2) = (/A(3,2)*A(1,3)-A(1,2)*A(3,3),&
              A(1,1)*A(3,3)-A(3,1)*A(1,3),&
              A(3,1)*A(1,2)-A(1,1)*A(3,2)/)
  invA(:,3) = (/A(1,2)*A(2,3)-A(2,2)*A(1,3),&
              A(2,1)*A(1,3)-A(1,1)*A(2,3),&
              A(1,1)*A(2,2)-A(2,1)*A(1,2)/)
  x = (/invA(1,1)*b(1)+invA(1,2)*b(2)+invA(1,3)*b(3),&
      invA(2,1)*b(1)+invA(2,2)*b(2)+invA(2,3)*b(3),&
      invA(3,1)*b(1)+invA(3,2)*b(2)+invA(3,3)*b(3)/)/det
  
end subroutine solve_3x3_linear_problem_r8

!> compute eigenvalues for 3x3 square matrices with real coeff.
subroutine compute_eigenvalues_real_3x3_r8(A,eigv)
  use mod_rootfinding, only: root_cubic
  implicit none
  !> inputs:
  real*8,dimension(3,3),intent(in)    :: A
  !> outputs:
  complex*16,dimension(3),intent(out) :: eigv
  !> variables:
  real*8,dimension(3) :: coeff
  !> compute eigenvalue problem coefficients
  coeff(1) = -1.d0*(A(1,1) + A(2,2) + A(3,3))
  coeff(2) = A(1,1)*A(2,2) + A(1,1)*A(3,3) + A(2,2)*A(3,3) - &
             (A(2,1)*A(1,2) + A(3,2)*A(2,3) + A(3,1)*A(1,3))
  coeff(3) = A(1,1)*A(3,2)*A(2,3) + A(2,1)*A(1,2)*A(3,3) + A(3,1)*A(2,2)*A(1,3) - &  
             (A(1,1)*A(2,2)*A(3,3) + A(2,1)*A(3,2)*A(1,3) + A(3,1)*A(1,2)*A(2,3))
  !> compute eigenvalues
  eigv = root_cubic(coeff)
end subroutine compute_eigenvalues_real_3x3_r8

end module mod_math_operators
