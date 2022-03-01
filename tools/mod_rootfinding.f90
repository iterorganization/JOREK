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
  public :: fun, dfun, ddfun, fun_dfun_3d
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

  !> type descrbing a 3d function and its jacobian 
  type,abstract :: fun_dfun_3d
    contains
    procedure(f_df_3d),pass,deferred :: f_df
  end type fun_dfun_3d

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
    pure subroutine f_df_3d(this,f,J,x)
      import fun_dfun_3d
      class(fun_dfun_3d),intent(inout)  :: this
      real*8,dimension(3),intent(out)   :: f
      real*8,dimension(3,3),intent(out) :: J
      real*8,dimension(3),intent(in)    :: x
    end subroutine f_df_3d
  end interface

  interface newtons_method
    module procedure newtons_method_f
    module procedure newtons_method_o
    module procedure newtons_method_3d_o
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

  !> implement the newton method for finding the roots of a system
  !> of three equations in three variables. A subroutine
  !> complying with the abstract class int_f_df_3d has to be provided
  !> inputs:
  !>   f_df:     (fun_dfun_3d) class computing function and derivative values
  !>   x_0:      (real8)(3) newton method first guess
  !>   y0:       (real8)(3) value of the root to find
  !>   tol_in:   (real8) tolerance
  !>   maxit_in: (integer)  maximum number of iteration
  !> outputs:
  !>   f_df: (f_df_3d) procedure computing function and derivative values
  !>   x:    (real8)(3) system root
  !>   ierr: (integer) 0 for success 1 otherwise
  pure subroutine newtons_method_3d_o(f_df_3d,y0,x0,x,ierr,tol_in,maxit_in)
    implicit none
    !> inputs-outputs:
    class(fun_dfun_3d),intent(inout) :: f_df_3d
    !> inputs:
    integer,intent(in)               :: maxit_in
    real*8,intent(in)                :: tol_in
    real*8,dimension(3),intent(in)   :: y0,x0
    !> outputs:
    integer,intent(out)              :: ierr
    real*8,dimension(3),intent(out)  :: x
    !> variables:
    integer               :: ii
    real*8                :: det
    real*8,dimension(3)   :: y
    real*8,dimension(3,3) :: J,invJ
    !> initialisation
    x=x0; ierr=0;
    do ii=1,maxit_in
     call f_df_3d%f_df(y,J,x)
     if(maxval(abs(y-y0)).lt.tol_in) return
     det = J(1,1)*(J(2,2)*J(3,3)-J(3,2)*J(2,3))+J(1,2)*(J(2,3)*J(3,1)-J(2,1)*J(3,3))+&
     J(1,3)*(J(2,1)*J(3,2)-J(2,2)*J(3,1))
     invJ(:,1) = (/J(2,2)*J(3,3)-J(3,2)*J(2,3),J(3,2)*J(1,3)-J(1,2)*J(3,3),J(1,2)*J(2,3)-J(2,2)*J(1,3)/)
     invJ(:,2) = (/J(2,3)*J(3,2)-J(2,1)*J(3,3),J(1,1)*J(3,3)-J(1,3)*J(3,1),J(1,3)*J(2,1)-J(2,3)*J(1,1)/)
     invJ(:,3) = (/J(2,1)*J(3,2)-J(2,2)*J(3,1),J(3,1)*J(1,2)-J(1,1)*J(3,2),J(1,1)*J(2,2)-J(1,2)*J(2,1)/)
     x = x - matmul(invJ,y-y0)/det
    enddo
  end subroutine newtons_method_3d_o

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
