!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
!
!>
!> @brief
!> A module for evaluating one-dimensional definite integrals with 
!> the adaptive Simpson rule
!<
module simpson_mod
  use func_mod

  implicit none
 
  integer, private, parameter :: simpson_default_max_depth = 10
  
  public :: simpson_adaptive, simpson_test

  private
  
contains
    
  !> Helper routine for "simpson_adaptive"
  recursive function simpson_helper(f, a, b, eps, S, fa, fb, fc, bottom) result (val)
    implicit none
    real*8, intent(in) :: a,b,eps,s,fa,fb,fc
    integer, intent(in) :: bottom
    real*8 :: val
    real*8 :: c, h, d, e,fd,fe, Sleft, Sright, S2
    procedure(func_real8_1d), pointer :: f

    c = (a + b)/2
    h = b - a
    d = (a + c)/2
    e = (c + b)/2
    fd = f(d) 
    fe = f(e);  
    Sleft = (h/12)*(fa + 4*fd + fc)
    Sright = (h/12)*(fc + 4*fe + fb)
    S2 = Sleft + Sright
        
    if (bottom .lt. 0 .or. abs(S2 - S) .le. eps*abs(S)) then
       val = S2 + (S2 - S)/15;      
    else
       val = simpson_helper(f, a, c, eps, Sleft,  fa, fc, fd, bottom-1) &
            &+ simpson_helper(f, c, b, eps, Sright, fc, fb, fe, bottom-1)  
    end if
  end function simpson_helper
  
  !> @brief Adaptive Simpsons rule for integral
  !<
  !< This function uses recursive splitting of the interval
  !< until either maximum number of intervals is reached or
  !< given relative error tolerance is obtained
  !<
  !< f: pointer to the integrand function (of one variable)
  !< a: lower limit for the interval
  !< b: upper limit for the interval
  !< eps: absolute tolerance
  !< maxDepth: maximum number of splits for the interval (if none given
  !<           then value 10 is used)
  !<
  function simpson_adaptive(f, a, b, eps, maxDepth) result (val)
    implicit none
    procedure(func_real8_1d), pointer :: f
    real*8, intent(in) :: a,b,eps
    integer, optional, intent(in) :: maxDepth
    
    real*8:: val
    integer :: depth
    
    ! Helper variables
    real*8 :: c,h,fa,fb,fc,S
    
    if(.not. present(maxDepth)) then 
       depth = simpson_default_max_depth
    else
       depth = maxDepth
    end if

    c = (a + b)/2
    h = b - a                                      
    fa = f(a) 
    fb = f(b) 
    fc = f(c)                                      
    S = (h/6)*(fa + 4*fc + fb);
    val = simpson_helper(f, a, b, eps, S, fa, fb, fc, depth);        
  end function simpson_adaptive

  !> test routine to integrate a couple of functions
  subroutine simpson_test()
    implicit none
    procedure(func_real8_1d), pointer :: f1_ptr => null(),f2_ptr => null(),f3_ptr => null()
    integer :: i
    integer :: maxdepth
    real*8 :: eps
    real*8 :: pi
    
    f1_ptr => f1
    f2_ptr => f2
    f3_ptr => f3

    pi=4.D0*atan(1.D0)
    maxdepth=100

    write(*,*) "Testing adaptive Simpsons rule."
    write(*,*) "Computing I = int_0^pi f(x) dx for"
    write(*,*) "f1(x)=x**2 with result I1=pi**3/3"
    write(*,*) "f2(x)=cos(x)**3*sin(12*x) with result I2=1096/6435" 
    write(*,*) "f3(x)=exp(x) with result I3=e^pi-1"
    write(*,*) "exact values are I1=",pi**3/3.D0," I2=",1096.D0/6435.D0," I3=",exp(pi)-1.D0

    loop_different_precisions : do i=1,15
       eps=10.D0**(-i)
       write(*,*) "eps=",eps, " I1=",simpson_adaptive(f1_ptr, 0.D0, pi, eps, maxdepth),&
            !&" I2=",simpson_adaptive(f2_ptr, 0.D0, pi, eps, maxdepth),&
            &" I3=",simpson_adaptive(f3_ptr, 0.D0, pi, eps, maxdepth)
    end do loop_different_precisions
    
  contains
    
    function f1(x) result (val)
      real*8, intent(in) :: x
      real*8 :: val
      val = x**2
    end function f1
    
    function f2(x) result (val)
      real*8, intent(in) :: x
      real*8 :: val
      val = cos(x)**3*sin(12*x)
    end function f2

    function f3(x) result (val)
      real*8, intent(in) :: x
      real*8 :: val
      val = exp(x)
    end function f3
    

  end subroutine simpson_test

end module simpson_mod
