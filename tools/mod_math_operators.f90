!> This module contains multi-purposes mathematical operators
module mod_math_operators

  implicit none
  private
  public cross_product
  
contains
  
!> The vector cross product
#ifdef CUDA_KERNELS
  attributes(host,device) &
#endif
  pure function cross_product(a, b)
  real*8, dimension(3) :: cross_product
  real*8, dimension(3), intent(in) :: a, b

  cross_product(1) = a(2) * b(3) - a(3) * b(2)
  cross_product(2) = a(3) * b(1) - a(1) * b(3)
  cross_product(3) = a(1) * b(2) - a(2) * b(1)

end function cross_product
  
end module mod_math_operators
