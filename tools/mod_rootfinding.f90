!> Module to find roots of functions.
!> We provide drivers for both Newton and Halley's method.
!> These come in two versions, bare versions (_f) taking just a function
!> and object versions (_o) taking a description object containing
!> the required functions and some data.
!>
!> Depending on the number of derivatives present you base your problem
!> class on fun, dfun or ddfun. This ensures that a solver will always
!> be able to call all derivatives, i.e. the Halley's method solver will
!> only operate on objects of class(ddfun)
!>
!> Note that all of these are 1D. An extension to multiple dimensions is
!> straightforward but not yet necessary.
module mod_rootfinding
  implicit none
  private
  public :: newtons_method
  public :: halleys_method
  public :: fun, dfun, ddfun
  public :: root

  !> Base type describing a function with optional private parameters
  !> Provide a guess of the inverse too, to serve as starting point
  type, abstract :: fun
  contains
    procedure(f), pass, deferred :: f
    procedure(inverse_f), pass, deferred :: inverse_f
  end type

  !> Type describing a function with parameters + its derivative
  type, abstract, extends(fun) :: dfun
  contains
    procedure(df), pass, deferred :: df
  end type

  !> Type describing a function with parameters + 2 derivatives
  type, abstract, extends(dfun) :: ddfun
  contains
    procedure(ddf), pass, deferred :: ddf
  end type


  interface
    pure function inverse_f(this, f) result(x)
      import fun
      class(fun), intent(in) :: this
      real*8, intent(in) :: f
      real*8 :: x
    end function inverse_f
    pure function f(this, x)
      import fun
      class(fun), intent(in) :: this
      real*8, intent(in) :: x
      real*8 :: f
    end function f
    pure function df(this, x)
      import dfun
      class(dfun), intent(in) :: this
      real*8, intent(in) :: x
      real*8 :: df
    end function df
    pure function ddf(this, x)
      import ddfun
      class(ddfun), intent(in) :: this
      real*8, intent(in) :: x
      real*8 :: ddf
    end function ddf
  end interface

  interface newtons_method
    module procedure newtons_method_f
    module procedure newtons_method_o
  end interface
  interface halleys_method
    module procedure halleys_method_f
    module procedure halleys_method_o
  end interface
contains

  !> Use newton's method to solve f(x) == y0, starting at x0
  pure subroutine newtons_method_f(f, df, y0, x0, x, ierr)
    real*8, intent(in)   :: y0    !< Intersection to find
    real*8, intent(in)   :: x0    !< Initial value
    real*8, intent(out)  :: x     !< Result value
    integer, intent(out) :: ierr  !< Status code. If == 0 we found a result
    interface
      pure function f(x)
        real*8, intent(in) :: x
        real*8 :: f
      end function f
      pure function df(x)
        real*8, intent(in) :: x
        real*8 :: df
      end function df
    end interface

    integer, parameter :: n_iter = 10
    real*8, parameter  :: tolerance = 1d-10

    integer :: i
    real*8  :: y
    
    ierr = 0
    x = x0
    do i=1,n_iter
      y = f(x) - y0
      x = x - y/df(x)

      if (abs(y) .le. tolerance) return
    end do
    ierr = 1 ! We did not find a root
  end subroutine newtons_method_f


  !> Use newton's method to solve f(x) == y0, starting at x0.
  !> Work from a class(dfun) object containing the functions to be called
  pure subroutine newtons_method_o(fun, y0, x, ierr)
    class(dfun), intent(in) :: fun
    real*8, intent(in)   :: y0    !< Intersection to find
    real*8, intent(out)  :: x     !< Result value
    integer, intent(out) :: ierr  !< Status code. If == 0 we found a result

    integer, parameter :: n_iter = 10
    real*8, parameter  :: tolerance = 1d-10

    integer :: i
    real*8  :: y
    
    ierr = 0
    x = fun%inverse_f(y0)
    do i=1,n_iter
      y = fun%f(x) - y0
      x = x - y/fun%df(x)

      if (abs(y) .le. tolerance) return
    end do
    ierr = 1 ! We did not find a root
  end subroutine newtons_method_o

  !> Use Halley's method to solve f(x) == y0, starting at x0
  pure subroutine halleys_method_f(f, df, ddf, y0, x0, x, ierr)
    real*8, intent(in)   :: y0    !< Intersection to find
    real*8, intent(in)   :: x0    !< Initial value
    real*8, intent(out)  :: x     !< Result value
    integer, intent(out) :: ierr  !< Status code. If == 0 we found a result
    interface
      pure function f(x)
        real*8, intent(in) :: x
        real*8 :: f
      end function f
      pure function df(x)
        real*8, intent(in) :: x
        real*8 :: df
      end function df
      pure function ddf(x)
        real*8, intent(in) :: x
        real*8 :: ddf
      end function ddf
    end interface

    integer, parameter :: n_iter = 10
    real*8, parameter  :: tolerance = 1d-10

    integer :: i
    real*8  :: y, dy, ddy
    
    ierr = 0
    x = x0
    do i=1,n_iter
      y   = f(x) - y0
      if (abs(y) .le. tolerance) return
      dy  = df(x)
      ddy = ddf(x)
      x = x - 2.d0*y*dy/(2.d0*dy*dy - y*ddy)
    end do
    ierr = 1 ! We did not find a root
  end subroutine halleys_method_f

  !> Use Halley's method to solve f(x) == y0, starting at x0
  !> Pass a ddfun object to encapsulate parameters
  pure subroutine halleys_method_o(fun, y0, x, ierr)
    class(ddfun), intent(in) :: fun !< Function description
    real*8, intent(in)   :: y0    !< Intersection to find
    real*8, intent(out)  :: x     !< Result value
    integer, intent(out) :: ierr  !< Status code. If == 0 we found a result

    integer, parameter :: n_iter = 10
    real*8, parameter  :: tolerance = 1d-10

    integer :: i
    real*8  :: y, dy, ddy
    
    ierr = 0
    x = fun%inverse_f(y0)
    do i=1,n_iter
      y   = fun%f(x) - y0
      if (abs(y) .le. tolerance) return
      dy  = fun%df(x)
      ddy = fun%ddf(x)
      x = x - 2.d0*y*dy/(2.d0*dy*dy - y*ddy)
    end do
    ierr = 1 ! We did not find a root
  end subroutine halleys_method_o


  !> Repeated here from solvers/root.f90 to be pure
  pure function root(A,B,C,D,SGN)
  !---------------------------------------------------------------------
  ! THIS FUNCTION GIVES BETTER ROOTS OF QUADRATICS BY AVOIDING
  ! CANCELLATION OF SMALLER ROOT
  ! Solve A x^2 + B x + C = 0
  ! D = B^2 - 4 A C
  !---------------------------------------------------------------------
  implicit none
  real*8, intent(in) :: a, b, c, d, sgn
  real*8 :: root

  if (((B .EQ. 0.D0) .and. (D .EQ. 0.D0)) .or. (A .eq. 0.D0)) then
   root = 1.d20 ! ill defined
   return
  endif
   
  if (B*SGN .GE. 0.d0) then
    root = -2.d0*C/(B+SGN*SQRT(D))
  else
    ROOT = (-B + SGN*SQRT(D)) / (2.d0 * A)
  endif
  return
  end function root
end module mod_rootfinding
