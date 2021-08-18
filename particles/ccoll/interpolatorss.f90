!##################################################

! This module contains "local" interpolation routines based on equally
! spaced grid points. The user must give the values of the function on
! a certain number of grid points per dimension (1-4) according to the
! precision order (linear=2 points, quadratic=3 points, cubic=4
! points, more details below). He must also give the "normalised" axis
! value for which to interpolate (so that the interpolation interval
! is [0,1]). The "normalised" axis value can be calculated with
! function NORM_IT.

! NOTICE: These routines do not deal with end points and beyond. The
! user must therefore take care himself of preparing function values
! beyond end points (using extrapolation or other methods,
! assumptions).

!##################################################

MODULE interpolatorss
  IMPLICIT NONE

  ! LINEAR = 2 points : lower and upper grid points
  INTEGER, PARAMETER :: NB_LINEAR = 2

  ! QUADRATIC = 3 points : previous, next and second next grid points
  ! Second order Taylor expansion: f(x) = f(x_0) + x Df + 1/2 x^2 D^2f
  ! with forward finite difference estimation of Df and D^2f
  INTEGER, PARAMETER :: NB_QUAD = 3

  ! CUBIC = 4 points : one before previous, previous, next and second
  ! next grid points. According to CATMULL-ROM formula (see Wikipedia)
  INTEGER, PARAMETER :: NB_CUBIC = 4

CONTAINS

  !################################################## 

  !Simple tool to calculate the normalised value on the axis given the
  ! grid points and the value x
  REAL(KIND=8) FUNCTION NORM_IT(x,grid)
    REAL(KIND=8), INTENT(in) :: x
    REAL(KIND=8), INTENT(in) :: grid(2)
    norm_it = (x-grid(1))/(grid(2)-grid(1))
  END FUNCTION NORM_IT

  !##################################################

  REAL(KIND=8) FUNCTION linear_interp(normed_pos,points)
    REAL(KIND=8), INTENT(in) :: points(NB_LINEAR)
    REAL(KIND=8), INTENT(in) :: normed_pos
    linear_interp = points(1)*(1.D0-normed_pos) + points(2)*normed_pos
  END FUNCTION linear_interp

  COMPLEX FUNCTION linear_interp_cpx(normed_pos,points)
    COMPLEX, INTENT(in) :: points(NB_LINEAR)
    REAL(KIND=8), INTENT(in) :: normed_pos
    linear_interp_cpx = points(1)*(1.D0-normed_pos) + points(2)*normed_pos
  END FUNCTION linear_interp_cpx

  REAL(KIND=8) FUNCTION cubic_interp(normed_pos,points)
    REAL(KIND=8), INTENT(in) :: points(NB_CUBIC)
    REAL(KIND=8), INTENT(in) :: normed_pos
    ! !EXPLICIT VERSION
    ! real(kind=8) :: a(NB_CUBIC),np2
    ! np2=normed_pos*normed_pos
    !! NAIVE
    ! a(1) = points(2)
    ! a(2) = points(3) - points(1)
    ! a(3) =-points(4) + points(3) -2.*points(2) + 2.*points(1) 
    ! a(4) = points(4) - points(3) + points(2) - points(1)
    !!CATMULL-ROM
    ! a(1) = points(2)
    ! a(2) = 0.5*(points(3) - points(1))
    ! a(3) = points(1) - 2.5*points(2) + 2.*points(3) - 0.5*points(4)
    ! a(4) = -0.5*points(1) + 1.D05*points(2) - 1.D05*points(3) + 0.5*points(4)
    !cubic_interp = a(1) + a(2)*normed_pos + a(3)*np2 + a(4)*np2*normed_pos
    !! SQUASHED CATMULL VERSION
    cubic_interp = points(2) + 0.5*normed_pos*(&
         points(3) - points(1) + normed_pos*(&
         2.*points(1) - 5.*points(2) + 4.*points(3) - points(4) + normed_pos*(&
         3.*(points(2)-points(3)) + points(4) - points(1))))
  END FUNCTION cubic_interp

  REAL(KIND=8) FUNCTION quad_interp(normed_pos,points)
    REAL(KIND=8), INTENT(IN):: points(NB_QUAD)
    REAL(KIND=8), INTENT(IN):: normed_pos
    quad_interp = points(1) + normed_pos*(&
         points(2)-points(1)+0.5*(normed_pos-1.D0)*(points(1)+points(3)-2*points(2)))
  END FUNCTION quad_interp

  !##################################################

  REAL(KIND=8) FUNCTION bilinear_interp(normed_x,normed_y,points)
    REAL(KIND=8), INTENT(IN) :: normed_x,normed_y
    REAL(KIND=8), INTENT(in) :: points(NB_LINEAR,NB_LINEAR)
    bilinear_interp = bi_interp(normed_x,normed_y,points,NB_LINEAR,linear_interp)
  END FUNCTION bilinear_interp

  REAL(KIND=8) FUNCTION trilinear_interp(normed_x,normed_y,normed_z,points)
    REAL(KIND=8), INTENT(IN) :: normed_x,normed_y,normed_z
    REAL(KIND=8), INTENT(in) :: points(NB_LINEAR,NB_LINEAR,NB_LINEAR)
    trilinear_interp = tri_interp(normed_x,normed_y,normed_z,points,NB_LINEAR,linear_interp)
  END FUNCTION trilinear_interp

  COMPLEX FUNCTION trilinear_interp_cpx(normed_x,normed_y,normed_z,points) !JF
    REAL(KIND=8), INTENT(IN) :: normed_x,normed_y,normed_z
    COMPLEX, INTENT(in) :: points(NB_LINEAR,NB_LINEAR,NB_LINEAR) !JF
    trilinear_interp_cpx = tri_interp_cpx(normed_x,normed_y,normed_z,points,NB_LINEAR,linear_interp_cpx)
  END FUNCTION trilinear_interp_cpx

  REAL(KIND=8) FUNCTION quadrilinear_interp(normed_x,normed_y,normed_z,normed_w,points)
    REAL(KIND=8), INTENT(IN) :: normed_x,normed_y,normed_z,normed_w
    REAL(KIND=8), INTENT(in) :: points(NB_LINEAR,NB_LINEAR,NB_LINEAR,NB_LINEAR)
    quadrilinear_interp = quadri_interp(normed_x,normed_y,normed_z,normed_w,points,NB_LINEAR,linear_interp)
  END FUNCTION quadrilinear_interp

  REAL(KIND=8) FUNCTION bicubic_interp(normed_x,normed_y,points)
    REAL(KIND=8), INTENT(IN) :: normed_x,normed_y
    REAL(KIND=8), INTENT(in) :: points(NB_CUBIC,NB_CUBIC)
    bicubic_interp = bi_interp(normed_x,normed_y,points,NB_CUBIC,cubic_interp)
  END FUNCTION bicubic_interp

  REAL(KIND=8) FUNCTION tricubic_interp(normed_x,normed_y,normed_z,points)
    REAL(KIND=8), INTENT(IN) :: normed_x,normed_y,normed_z
    REAL(KIND=8), INTENT(in) :: points(NB_CUBIC,NB_CUBIC,NB_CUBIC)
    tricubic_interp = tri_interp(normed_x,normed_y,normed_z,points,NB_CUBIC,cubic_interp)
  END FUNCTION tricubic_interp

  REAL(KIND=8) FUNCTION quadricubic_interp(normed_x,normed_y,normed_z,normed_w,points)
    REAL(KIND=8), INTENT(IN) :: normed_x,normed_y,normed_z,normed_w
    REAL(KIND=8), INTENT(in) :: points(NB_CUBIC,NB_CUBIC,NB_CUBIC,NB_CUBIC)
    quadricubic_interp = quadri_interp(normed_x,normed_y,normed_z,normed_w,points,NB_CUBIC,cubic_interp)
  END FUNCTION quadricubic_interp

  REAL(KIND=8) FUNCTION biquad_interp(normed_x,normed_y,points)
    REAL(KIND=8), INTENT(IN) :: normed_x,normed_y
    REAL(KIND=8), INTENT(in) :: points(NB_QUAD,NB_QUAD)
    biquad_interp = bi_interp(normed_x,normed_y,points,NB_QUAD,quad_interp)
  END FUNCTION biquad_interp

  REAL(KIND=8) FUNCTION triquad_interp(normed_x,normed_y,normed_z,points)
    REAL(KIND=8), INTENT(IN) :: normed_x,normed_y,normed_z
    REAL(KIND=8), INTENT(in) :: points(NB_QUAD,NB_QUAD,NB_QUAD)
    triquad_interp = tri_interp(normed_x,normed_y,normed_z,points,NB_QUAD,quad_interp)
  END FUNCTION triquad_interp

  REAL(KIND=8) FUNCTION quadriquad_interp(normed_x,normed_y,normed_z,normed_w,points)
    REAL(KIND=8), INTENT(IN) :: normed_x,normed_y,normed_z,normed_w
    REAL(KIND=8), INTENT(in) :: points(NB_QUAD,NB_QUAD,NB_QUAD,NB_QUAD)
    quadriquad_interp = quadri_interp(normed_x,normed_y,normed_z,normed_w,points,NB_QUAD,quad_interp)
  END FUNCTION quadriquad_interp

  !##################################################

  ! Higher dimension interpolation routines. The procedure is
  ! recursively general.

  REAL(KIND=8) FUNCTION bi_interp(normed_x,normed_y,points,NB_POINTS,which_interp)
    INTEGER, INTENT(in) :: NB_POINTS
    REAL(KIND=8), INTENT(IN) :: normed_x,normed_y
    REAL(KIND=8), INTENT(in) :: points(NB_POINTS,NB_POINTS)
    REAL(KIND=8), EXTERNAL :: which_interp
    REAL(KIND=8) :: compress(NB_POINTS)
    INTEGER :: i
    DO i=1,NB_POINTS
       compress(i) = which_interp(normed_x,points(:,i))
    END DO
    bi_interp = which_interp(normed_y,compress)
  END FUNCTION bi_interp

  COMPLEX FUNCTION bi_interp_cpx(normed_x,normed_y,points,NB_POINTS,which_interp) !JF
    INTEGER, INTENT(in) :: NB_POINTS
    REAL(KIND=8), INTENT(IN) :: normed_x,normed_y
    COMPLEX, INTENT(in) :: points(NB_POINTS,NB_POINTS)
    COMPLEX, EXTERNAL :: which_interp
    COMPLEX :: compress(NB_POINTS)
    INTEGER :: i
    DO i=1,NB_POINTS
       compress(i) = which_interp(normed_x,points(:,i))
    END DO
    bi_interp_cpx = which_interp(normed_y,compress)
  END FUNCTION bi_interp_cpx

  REAL(KIND=8) FUNCTION tri_interp(normed_x,normed_y,normed_z,points,NB_POINTS,which_interp)
    INTEGER, INTENT(IN)::NB_POINTS
    REAL(KIND=8), INTENT(IN) :: normed_x,normed_y,normed_z
    REAL(KIND=8), INTENT(IN) :: points(NB_POINTS,NB_POINTS,NB_POINTS)
    REAL(KIND=8), EXTERNAL :: which_interp
    REAL(KIND=8) :: compress(NB_POINTS)
    INTEGER::i
    DO i=1,NB_POINTS
       compress(i) = bi_interp(normed_x,normed_y,points(:,:,i),NB_POINTS,which_interp)
    END DO
    tri_interp = which_interp(normed_z,compress)
  END FUNCTION tri_interp

  COMPLEX FUNCTION tri_interp_cpx(normed_x,normed_y,normed_z,points,NB_POINTS,which_interp) !JF
    INTEGER, INTENT(IN)::NB_POINTS
    REAL(KIND=8), INTENT(IN) :: normed_x,normed_y,normed_z
    COMPLEX, INTENT(IN) :: points(NB_POINTS,NB_POINTS,NB_POINTS)
    COMPLEX, EXTERNAL :: which_interp
    COMPLEX :: compress(NB_POINTS)
    INTEGER::i
    DO i=1,NB_POINTS
       compress(i) = bi_interp_cpx(normed_x,normed_y,points(:,:,i),NB_POINTS,which_interp)
    END DO
    tri_interp_cpx = which_interp(normed_z,compress)
  END FUNCTION tri_interp_cpx

  REAL(KIND=8) FUNCTION quadri_interp(normed_x,normed_y,normed_z,normed_w,points,NB_POINTS,which_interp)
    INTEGER, INTENT(IN)::NB_POINTS
    REAL(KIND=8), INTENT(IN) :: normed_x,normed_y,normed_z,normed_w
    REAL(KIND=8), INTENT(IN) :: points(NB_POINTS,NB_POINTS,NB_POINTS,NB_POINTS)
    REAL(KIND=8), EXTERNAL :: which_interp
    REAL(KIND=8) :: compress(NB_POINTS)
    INTEGER::i
    DO i=1,NB_POINTS
       compress(i) = tri_interp(normed_x,normed_y,normed_z,points(:,:,:,i),NB_POINTS,which_interp)
    END DO
    quadri_interp = which_interp(normed_w,compress)
  END FUNCTION quadri_interp

  !##################################################

  ! ARRAY INTERPOLATION whenever the arguments has an extra array
  ! dimension (like components of vector field). It assumes the extra
  ! components are stored in the first index

  FUNCTION linear_array_interp(normed_pos,points,array_length)
    INTEGER,INTENT(IN) :: array_length
    REAL(KIND=8), INTENT(in) :: points(ARRAY_LENGTH,NB_LINEAR)
    REAL(KIND=8), INTENT(in) :: normed_pos
    REAL(KIND=8) :: linear_array_interp(ARRAY_LENGTH)
    linear_array_interp = points(:,1)*(1.D0-normed_pos) + points(:,2)*normed_pos
  END FUNCTION linear_array_interp

  FUNCTION cubic_array_interp(normed_pos,points,array_length)
    INTEGER,INTENT(IN) :: array_length
    REAL(KIND=8), INTENT(in) :: points(array_length,NB_CUBIC)
    REAL(KIND=8), INTENT(in) :: normed_pos
    REAL(KIND=8) :: cubic_array_interp(array_length)
    ! !EXPLICIT VERSION
    ! real(kind=8) :: a(NB_CUBIC),np2
    ! np2=normed_pos*normed_pos
    !! NAIVE
    ! a(1) = points(2)
    ! a(2) = points(3) - points(1)
    ! a(3) =-points(4) + points(3) -2.*points(2) + 2.*points(1) 
    ! a(4) = points(4) - points(3) + points(2) - points(1)
    !!CATMULL-ROM
    ! a(1) = points(2)
    ! a(2) = 0.5*(points(3) - points(1))
    ! a(3) = points(1) - 2.5*points(2) + 2.*points(3) - 0.5*points(4)
    ! a(4) = -0.5*points(1) + 1.D05*points(2) - 1.D05*points(3) + 0.5*points(4)
    !cubic_interp = a(1) + a(2)*normed_pos + a(3)*np2 + a(4)*np2*normed_pos
    !! SQUASHED CATMULL VERSION
    cubic_array_interp = points(:,2) + 0.5*normed_pos*(&
         points(:,3) - points(:,1) + normed_pos*(&
         2.*points(:,1) - 5.*points(:,2) + 4.*points(:,3) - points(:,4) + normed_pos*(&
         3.*(points(:,2)-points(:,3)) + points(:,4) - points(:,1))))
  END FUNCTION cubic_array_interp

  FUNCTION quad_array_interp(normed_pos,points,array_length)
    INTEGER,INTENT(IN) :: array_length
    REAL(KIND=8), INTENT(in) :: points(ARRAY_LENGTH,NB_QUAD)
    REAL(KIND=8), INTENT(in) :: normed_pos
    REAL(KIND=8) :: quad_array_interp(ARRAY_LENGTH)
    quad_array_interp = points(:,1) + normed_pos*(&
         points(:,2)-points(:,1)+0.5*(normed_pos-1.D0)*(points(:,1)+points(:,3)-2*points(:,2)))
  END FUNCTION quad_array_interp

  FUNCTION bilinear_array_interp(normed_x,normed_y,points,array_length)
    INTEGER,INTENT(IN) :: array_length
    REAL(KIND=8), INTENT(IN) :: normed_x,normed_y
    REAL(KIND=8), INTENT(in) :: points(ARRAY_LENGTH,NB_LINEAR,NB_LINEAR)
    INTEGER :: i
    REAL(KIND=8) :: compress(array_length,NB_LINEAR)
    REAL(KIND=8) :: bilinear_array_interp(array_length)
    DO i=1,NB_LINEAR
       compress(:,i) = linear_array_interp(normed_x,points(:,:,i),array_length)
    END DO
    bilinear_array_interp = linear_array_interp(normed_y,compress,array_length)
  END FUNCTION bilinear_array_interp

  FUNCTION trilinear_array_interp(normed_x,normed_y,normed_z,points,array_length)
    INTEGER,INTENT(IN) :: array_length
    REAL(KIND=8), INTENT(IN) :: normed_x,normed_y,normed_z
    REAL(KIND=8), INTENT(IN) :: points(array_length,NB_LINEAR,NB_LINEAR,NB_LINEAR)
    INTEGER::i
    REAL(KIND=8) :: compress(array_length,NB_LINEAR)
    REAL(KIND=8) :: trilinear_array_interp(array_length)
    DO i=1,NB_LINEAR
       compress(:,i) = bilinear_array_interp(normed_x,normed_y,points(:,:,:,i),array_length)
    END DO
    trilinear_array_interp = linear_array_interp(normed_z,compress,array_length)
  END FUNCTION trilinear_array_interp

  FUNCTION biquad_array_interp(normed_x,normed_y,points,array_length)
    INTEGER,INTENT(IN) :: array_length
    REAL(KIND=8), INTENT(IN) :: normed_x,normed_y
    REAL(KIND=8), INTENT(in) :: points(ARRAY_LENGTH,NB_QUAD,NB_QUAD)
    INTEGER :: i
    REAL(KIND=8) :: compress(array_length,NB_QUAD)
    REAL(KIND=8) :: biquad_array_interp(array_length)
    DO i=1,NB_QUAD
       compress(:,i) = quad_array_interp(normed_x,points(:,:,i),array_length)
    END DO
    biquad_array_interp = quad_array_interp(normed_y,compress,array_length)
  END FUNCTION biquad_array_interp

  FUNCTION triquad_array_interp(normed_x,normed_y,normed_z,points,array_length)
    INTEGER,INTENT(IN) :: array_length
    REAL(KIND=8), INTENT(IN) :: normed_x,normed_y,normed_z
    REAL(KIND=8), INTENT(IN) :: points(array_length,NB_QUAD,NB_QUAD,NB_QUAD)
    INTEGER::i
    REAL(KIND=8) :: compress(array_length,NB_QUAD)
    REAL(KIND=8) :: triquad_array_interp(array_length)
    DO i=1,NB_QUAD
       compress(:,i) = biquad_array_interp(normed_x,normed_y,points(:,:,:,i),array_length)
    END DO
    triquad_array_interp = quad_array_interp(normed_z,compress,array_length)
  END FUNCTION triquad_array_interp

  FUNCTION bicubic_array_interp(normed_x,normed_y,points,array_length)
    INTEGER,INTENT(IN) :: array_length
    REAL(KIND=8), INTENT(IN) :: normed_x,normed_y
    REAL(KIND=8), INTENT(in) :: points(ARRAY_LENGTH,NB_CUBIC,NB_CUBIC)
    INTEGER :: i
    REAL(KIND=8) :: compress(array_length,NB_CUBIC)
    REAL(KIND=8) :: bicubic_array_interp(array_length)
    DO i=1,NB_CUBIC
       compress(:,i) = cubic_array_interp(normed_x,points(:,:,i),array_length)
    END DO
    bicubic_array_interp = cubic_array_interp(normed_y,compress,array_length)
  END FUNCTION bicubic_array_interp

  FUNCTION tricubic_array_interp(normed_x,normed_y,normed_z,points,array_length)
    INTEGER,INTENT(IN) :: array_length
    REAL(KIND=8), INTENT(IN) :: normed_x,normed_y,normed_z
    REAL(KIND=8), INTENT(IN) :: points(array_length,NB_CUBIC,NB_CUBIC,NB_CUBIC)
    INTEGER::i
    REAL(KIND=8) :: compress(array_length,NB_CUBIC)
    REAL(KIND=8) :: tricubic_array_interp(array_length)
    DO i=1,NB_CUBIC
       compress(:,i) = bicubic_array_interp(normed_x,normed_y,points(:,:,:,i),array_length)
    END DO
    tricubic_array_interp = cubic_array_interp(normed_z,compress,array_length)
  END FUNCTION tricubic_array_interp

  !==================================================
  ! THIS ROUTINE HAS NOT BEEN FULLY TESTED (WORKS IN PRINCIPAL)
  REAL(KIND=8) FUNCTION cubic_irregular_interp(x,xgrid,points)
    REAL(KIND=8), INTENT(IN) :: x
    REAL(KIND=8), INTENT(IN) :: xgrid(NB_CUBIC)
    REAL(KIND=8), INTENT(IN) :: points(NB_CUBIC)

    REAL(KIND=8) :: t
    !REAL(KIND=8) :: h(NB_CUBIC)   ! necessary array for CATMULL-ROM and optimized version
    REAL(KIND=8) :: m(NB_LINEAR)

    t = (x-xgrid(2))/(xgrid(3)-xgrid(2))

    ! TANGEANTS
    ! FINITE DIFFERENCE AVERAGE
    ! m(1) = 0.5*(&
    !      (points(3)-points(2))/(xgrid(3)-xgrid(2))&
    !      +(points(2)-points(1))/(xgrid(2)-xgrid(1)))
    ! m(2) = 0.5*(&
    !      (points(4)-points(3))/(xgrid(4)-xgrid(3))&
    !      +(points(3)-points(2))/(xgrid(3)-xgrid(2)))
    ! 
    ! FINITE DIFFERENCE LARGER
    m(1) = (points(3)-points(1))/(xgrid(3)-xgrid(1))
    m(2) = (points(4)-points(2))/(xgrid(4)-xgrid(2))

    ! CUBIC SPLINE
    cubic_irregular_interp=(1.D0-t)*points(2) + t*points(3) &
         +t*(1.D0-t)*(  (m(1)*(xgrid(3)-xgrid(2))-points(3)+points(2))*(1.D0-t) &
         + (-m(2)*(xgrid(3)-xgrid(2))+points(3)-points(2))*t)

    ! CATUMULL-ROM (CUBIC HERMITE SPLINE)
    ! h(1) = (1.D0+2.*t)*(1.D0-t)**2
    ! h(2) = t*(1.D0-t)**2
    ! h(3) = t**2*(3.-2.*t)
    ! h(4) = t**2*(t-1.D0)
    ! cubic_irregular_interp = h(1)*points(2) &
    !      + h(2)*m(1)*(xgrid(3)-xgrid(2)) &
    !      + h(3)*points(3) &
    !      + h(4)*m(2)*(xgrid(3)-xgrid(2))   

    ! OPTIMIZED VERSION
    ! h(1) = t-1.D0
    ! h(2) = h(1)*h(1)
    ! h(3) = t*t
    ! h(4) = h(3)*h(1)
    ! h(1) = (1.D0+2.*t)*h(2)
    ! h(2) = t*h(2)
    ! h(3) = h(3)*(3.-2.*t)
    ! t = h(2)/(xgrid(4)-xgrid(2))
    ! h(4) = h(4)/(xgrid(3)-xgrid(1))
    ! h(2) = h(1)-h(4)
    ! h(3) = h(3)+t
    ! h(1) = -t
    !cubic_irregular_interp = DOT_PRODUCT(h,points)       
  END FUNCTION cubic_irregular_interp

  real(kind=8) function neville_polynomial_interp(x,xgrid,points,npoints)
    INTEGER,INTENT(IN)::npoints
    REAL(KIND=8),INTENT(IN) :: x
    REAL(KIND=8),INTENT(IN) :: XGRID(NPOINTS),POINTS(NPOINTS)

    INTEGER :: i,j
    REAL(KIND=8) :: polynome(NPOINTS)

    polynome=points

    DO J=1,NPOINTS
       DO I=1,NPOINTS-j
          polynome(i) = ((xgrid(i)-x)*polynome(i+1) + (x-xgrid(i+j))*polynome(i))/(xgrid(i)-xgrid(i+j))
       END DO
    END DO
    NEVILLE_POLYNOMIAL_INTERP=POLYNOME(1)
  END FUNCTION neville_polynomial_interp

  function neville_array_interp(x,xgrid,points,npoints,array_length)
    INTEGER,INTENT(IN)::npoints
    INTEGER,INTENT(IN)::array_length
    REAL(KIND=8),INTENT(IN) :: x
    REAL(KIND=8),INTENT(IN) :: XGRID(NPOINTS)
    REAL(KIND=8),INTENT(IN) :: POINTS(array_length,NPOINTS)

    REAL(KIND=8) neville_array_interp(array_length)
    
    INTEGER :: i,j
    REAL(KIND=8) :: polynome(array_length,NPOINTS)

    polynome=points

    DO J=1,NPOINTS
       DO I=1,NPOINTS-j
          polynome(:,i) = ((xgrid(i)-x)*polynome(:,i+1) + (x-xgrid(i+j))*polynome(:,i))/(xgrid(i)-xgrid(i+j))
       END DO
    END DO
    neville_array_interp=POLYNOME(:,1)
  END FUNCTION neville_array_interp

END MODULE interpolatorss
