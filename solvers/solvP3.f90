SUBROUTINE SOLVP3(C0,C1,C2,C3,X1,X2,X3,IFAIL)
!-----------------------------------------------------------------------
! SOLVES A CUBIC EQUATION WITH A SOLUTION WITH -1.< X < 1
! CN : THE COEFFICIENT OF X**N, X : THE REAL SOLUTION WITH -1.< X < 1.
!-----------------------------------------------------------------------
implicit none
real*8    :: C0, C1, C2, C3, X1, X2, X3, dum, tol, angle
real*8    :: a0, a1, a2, aa, bb, cc, det, pi, p, q, u, v
integer :: ifail
real*8, external :: root

x1    = 99.d0
x2    = 999.d0
x3    = 9999.d0
tol   = 0.d0
ifail = 0

IF ((ABS(C1)+ABS(C2)+ABS(C3)) .EQ. 0.d0) THEN
  IFAIL = 1
  RETURN
ENDIF

!------------------------------------- 2nd order poly for small c3
IF (ABS(C3)/(ABS(C1)+ABS(C2)+ABS(C3)) .LT. 1.d-9) THEN
  AA = C2
  BB = C1
  CC = C0
  DET = BB**2 - 4.d0*AA*CC
  IF (DET.GE.0.d0) THEN
    X1 = ROOT(AA,BB,CC,DET,1.d0)
    IF (ABS(X1).GT. 1.d0 + TOL) THEN
      X1 = ROOT(AA,BB,CC,DET,-1.d0)
    ENDIF
  ELSE
    IFAIL = 1
  ENDIF

ELSE
!------------------------------------- 3rd order poly solution
  PI = 2.d0*ASIN(1.d0)
  A0 = C0 / C3
  A1 = C1 / C3
  A2 = C2 / C3
  P = - (A2**2)/3.d0 + A1
  Q = 2.d0/27.d0*(A2**3) - A2 * A1/3.d0 + A0
  DET = (P/3.d0)**3 + (Q/2.d0)**2
  IF (DET .GE. 0.d0) THEN
    U  = SIGN(1.d0,-Q/2.d0+SQRT(DET))*ABS(-Q/2.d0 + SQRT(DET))**(1.d0/3.d0)
    V  = SIGN(1.d0,-Q/2.d0-SQRT(DET))*ABS(-Q/2.d0 - SQRT(DET))**(1.d0/3.d0)
    X1 =  U + V - A2/3.d0
    IF (ABS(X1) .GE. (1.d0+TOL)) IFAIL = 1
  ELSE
    P = -P
    ANGLE = SIGN(1.d0,P)*ACOS((Q/2.d0)/SQRT(ABS(P)/3.d0)**3)
    X1 = -2.d0*SQRT(ABS(P)/3.d0)*COS(ANGLE/3.d0) - A2/3.d0
    X2 = -2.d0*SQRT(ABS(P)/3.d0)*COS(2.d0*PI/3.d0 - ANGLE/3.d0) - A2/3.d0
    X3 = -2.d0*SQRT(ABS(P)/3.d0)*COS(2.d0*PI/3.d0 + ANGLE/3.d0) - A2/3.d0
  ENDIF
  IF (ABS(X1) .GT. ABS(X2)) THEN
    DUM = X1
    X1 = X2
    X2 = DUM
  ENDIF
  IF (ABS(X2) .GT. ABS(X3)) THEN
    DUM = X2
    X2 = X3
    X3 = DUM
  ENDIF
  IF (ABS(X1) .GT. ABS(X2)) THEN
    DUM = X1
    X1 = X2
    X2 = DUM
  ENDIF
ENDIF
IF (ABS(X1) .GT. (1.d0 + TOL)) IFAIL=1

RETURN
END
