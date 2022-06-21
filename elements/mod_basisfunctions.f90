!> Module providing the different implementations of the basis functions derived from a
!> mixed Bezier/Cubic finite element representation.
!> \see ::basisfunctions and ::basisfunctions1
module mod_basisfunctions
implicit none

private
public :: basisfunctions1, basisfunctions
public :: basisfunctions_T !< Transposed version, for faster interp_PRZ
public :: basisfunctions3

!> One-dimensional basisfunctions with derivatives of order n
interface basisfunctions1
  module procedure basisfunctions_1D_0, basisfunctions_1D_1, basisfunctions_1D_2
end interface basisfunctions1

! Two-dimensional basisfunctions with derivatives of order n
! basisfunctions_2D_1 has only s and t derivatives
! basisfunctions_2D_1p includes the st cross-derivative too
! basisfunctions_2D_2 includes the ss and tt derivatives additionally
interface basisfunctions
  module procedure basisfunctions_2D_0, basisfunctions_2D_1, basisfunctions_2D_1p, basisfunctions_2D_2
end interface basisfunctions

!> Two dimensional basisfunction with derivatives of order n
!> and transposed matrix for better vectorisation
!> basisfunctions_2D_1_T: first order derivatives in s and t
!> basisfunctions_2D_2_T: first and second order derivatives in s and t
interface basisfunctions_T
  module procedure basisfunctions_2D_1_T, basisfunctions_2D_2_T 
end interface basisfunctions_T

contains

!> Subroutine which defines the basis functions in one dimension with no derivatives
pure subroutine basisfunctions_1D_0(s,H)
real*8, intent(in)  :: s          !< s-coordinate in the element (in [0,1])
real*8, intent(out) :: H(2,2)     !< Basis functions
!---------------------------------------------------------- vertex (1)
H(1,1)   =(-1.d0 + s)**2*(1.d0 + 2.d0*s)
H(1,2)   =3.d0*(-1.d0 + s)**2*s
!---------------------------------------------------------- vertex (2)
H(2,1)   =-s**2*(-3.d0 + 2.d0*s)
H(2,2)   =-3.d0*(-1.d0 + s)*s**2
end subroutine basisfunctions_1D_0

!> Subroutine which defines the basis functions in one dimension with first derivatives
pure subroutine basisfunctions_1D_1(s,H,H_s)
real*8, intent(in)  :: s          !< s-coordinate in the element (in [0,1])
real*8, intent(out) :: H(2,2)     !< Basis functions
real*8, intent(out) :: H_s(2,2)   !< Basis functions derived with respect to s
call basisfunctions_1D_0(s,H)
!---------------------------------------------------------- vertex (1)
H_s(1,1) =6.d0*(-1.d0 + s)*s
H_s(1,2) =3.d0*(-1.d0 + s)*(-1.d0 + 3.d0*s)
!---------------------------------------------------------- vertex (2)
H_s(2,1) =-6.d0*(-1.d0 + s)*s
H_s(2,2) =-3.d0*s*(-2.d0 + 3.d0*s)
end subroutine basisfunctions_1D_1

!> Subroutine which defines the basis functions in one dimension with first and second derivatives
pure subroutine basisfunctions_1D_2(s,H,H_s,H_ss)
real*8, intent(in)  :: s          !< s-coordinate in the element (in [0,1])
real*8, intent(out) :: H(2,2)     !< Basis functions
real*8, intent(out) :: H_s(2,2)   !< Basis functions derived with respect to s
real*8, intent(out) :: H_ss(2,2)  !< Basis functions derived two times with respect to s
call basisfunctions_1D_1(s,H,H_s)
!---------------------------------------------------------- vertex (1)
H_ss(1,1)=6.d0*(-1.d0 + 2.d0*s)
H_ss(1,2)=6.d0*(-2.d0 + 3.d0*s)
!---------------------------------------------------------- vertex (2)
H_ss(2,1)=-6.d0*(-1.d0 + 2.d0*s)
H_ss(2,2)=-6.d0*(-1.d0 + 3.d0*s)
end subroutine basisfunctions_1D_2




!> Basisfunctions in 2D, value only.
pure subroutine basisfunctions_2D_0(s, t, H)
implicit none
real*8, intent(in)  :: s          !< s-coordinate in the element [0,1]
real*8, intent(in)  :: t          !< t-coordinate in the element [0,1]
real*8, intent(out) :: H(4,4)     !< Basis functions

!---------------------------------------------------------- vertex (1)
H(1,1)   =(-1.d0 + s)**2*(1.d0 + 2.d0*s)*(-1.d0 + t)**2*(1.d0 + 2.d0*t)
H(1,2)   =3.d0*(-1.d0 + s)**2*s*(-1.d0 + t)**2*(1.d0 + 2.d0*t)
H(1,3)   =3.d0*(-1.d0 + s)**2*(1.d0 + 2.d0*s)*(-1.d0 + t)**2*t
H(1,4)   =9.d0*(-1.d0 + s)**2*s*(-1.d0 + t)**2*t
!---------------------------------------------------------- vertex (2)
H(2,1)   =-(s**2*(-3.d0 + 2.d0*s)*(-1.d0 + t)**2*(1.d0 + 2.d0*t))
H(2,2)   =-3.d0*(-1.d0 + s)*s**2*(-1.d0 + t)**2*(1.d0 + 2.d0*t)
H(2,3)   =-3.d0*s**2*(-3.d0 + 2.d0*s)*(-1.d0 + t)**2*t
H(2,4)   =-9.d0*(-1.d0 + s)*s**2*(-1.d0 + t)**2*t
!---------------------------------------------------------- vertex (3)
H(3,1)   =s**2*(-3.d0 + 2.d0*s)*t**2*(-3.d0 + 2.d0*t)
H(3,2)   =3.d0*(-1.d0 + s)*s**2*t**2*(-3.d0 + 2.d0*t)
H(3,3)   =3.d0*s**2*(-3.d0 + 2.d0*s)*(-1.d0 + t)*t**2
H(3,4)   =9.d0*(-1.d0 + s)*s**2*(-1.d0 + t)*t**2
!---------------------------------------------------------- vertex (4)
H(4,1)   =-((-1.d0 + s)**2*(1.d0 + 2.d0*s)*t**2*(-3.d0 + 2.d0*t))
H(4,2)   =-3.d0*(-1.d0 + s)**2*s*t**2*(-3.d0 + 2.d0*t)
H(4,3)   =-3.d0*(-1.d0 + s)**2*(1.d0 + 2.d0*s)*(-1.d0 + t)*t**2
H(4,4)   =-9.d0*(-1.d0 + s)**2*s*(-1.d0 + t)*t**2
end subroutine basisfunctions_2D_0



pure subroutine basisfunctions_2D_1(s, t, H, H_s, H_t)
implicit none
real*8, intent(in)  :: s          !< s-coordinate in the element [0,1]
real*8, intent(in)  :: t          !< t-coordinate in the element [0,1]
real*8, intent(out) :: H(4,4)     !< Basis functions
real*8, intent(out) :: H_s(4,4)   !< Basis functions derived with respect to s
real*8, intent(out) :: H_t(4,4)   !< Basis functions derived with respect to t

call basisfunctions_2D_0(s, t, H)

!---------------------------------------------------------- vertex (1)
H_s(1,1) =6.d0*(-1.d0 + s)*s*(-1.d0 + t)**2*(1.d0 + 2.d0*t)
H_t(1,1) =6.d0*(-1.d0 + s)**2*(1.d0 + 2.d0*s)*(-1.d0 + t)*t

H_s(1,2) =3.d0*(-1.d0 + s)*(-1.d0 + 3.d0*s)*(-1.d0 + t)**2*(1.d0 + 2.d0*t)
H_t(1,2) =18.d0*(-1.d0 + s)**2*s*(-1.d0 + t)*t

H_s(1,3) =18.d0*(-1.d0 + s)*s*(-1.d0 + t)**2*t
H_t(1,3) =3.d0*(-1.d0 + s)**2*(1.d0 + 2.d0*s)*(-1.d0 + t)*(-1.d0 + 3.d0*t)

H_s(1,4) =9.d0*(-1.d0 + s)*(-1.d0 + 3.d0*s)*(-1.d0 + t)**2*t
H_t(1,4) =9.d0*(-1.d0 + s)**2*s*(-1.d0 + t)*(-1.d0 + 3.d0*t)

!---------------------------------------------------------- vertex (2)
H_s(2,1) =-6.d0*(-1.d0 + s)*s*(-1.d0 + t)**2*(1.d0 + 2.d0*t)
H_t(2,1) =-6.d0*s**2*(-3.d0 + 2.d0*s)*(-1.d0 + t)*t

H_s(2,2) =-3.d0*s*(-2.d0 + 3.d0*s)*(-1.d0 + t)**2*(1.d0 + 2.d0*t)
H_t(2,2) =-18.d0*(-1.d0 + s)*s**2*(-1.d0 + t)*t

H_s(2,3) =-18.d0*(-1.d0 + s)*s*(-1.d0 + t)**2*t
H_t(2,3) =3.d0*s**2*(-3.d0 + 2.d0*s)*(1.d0 - 3.d0*t)*(-1.d0 + t)

H_s(2,4) =-9.d0*s*(-2.d0 + 3.d0*s)*(-1.d0 + t)**2*t
H_t(2,4) =9.d0*(-1.d0 + s)*s**2*(1.d0 - 3.d0*t)*(-1.d0 + t)

!---------------------------------------------------------- vertex (3)
H_s(3,1) =6.d0*(-1.d0 + s)*s*t**2*(-3.d0 + 2.d0*t)
H_t(3,1) =6.d0*s**2*(-3.d0 + 2.d0*s)*(-1.d0 + t)*t

H_s(3,2) =3.d0*s*(-2 + 3.d0*s)*t**2*(-3.d0 + 2.d0*t)
H_t(3,2) =18.d0*(-1.d0 + s)*s**2*(-1.d0 + t)*t

H_s(3,3) =18.d0*(-1.d0 + s)*s*(-1.d0 + t)*t**2
H_t(3,3) =3.d0*s**2*(-3.d0 + 2.d0*s)*t*(-2.d0 + 3.d0*t)

H_s(3,4) =9.d0*s*(-2.d0 + 3.d0*s)*(-1.d0 + t)*t**2
H_t(3,4) =9.d0*(-1.d0 + s)*s**2*t*(-2.d0 + 3.d0*t)

!---------------------------------------------------------- vertex (4)
H_s(4,1) =-6.d0*(-1.d0 + s)*s*t**2*(-3.d0 + 2.d0*t)
H_t(4,1) =-6.d0*(-1.d0 + s)**2*(1.d0 + 2.d0*s)*(-1.d0 + t)*t

H_s(4,2) =3.d0*(1.d0 - 3*s)*(-1.d0 + s)*t**2*(-3.d0 + 2.d0*t)
H_t(4,2) =-18.d0*(-1.d0 + s)**2*s*(-1.d0 + t)*t

H_s(4,3) =-18.d0*(-1.d0 + s)*s*(-1.d0 + t)*t**2
H_t(4,3) =-3.d0*(-1.d0 + s)**2*(1.d0 + 2.d0*s)*t*(-2.d0 + 3.d0*t)

H_s(4,4) =9.d0*(1.d0 - 3.d0*s)*(-1.d0 + s)*(-1.d0 + t)*t**2
H_t(4,4) =-9.d0*(-1.d0 + s)**2*s*t*(-2.d0 + 3.d0*t)
end subroutine basisfunctions_2D_1



!> Transposed from the normal usage for easier vector operations (i.e.  basisfunctions_2D_1.T)
!> - index 1 : counts the variables (p,u,v,w)
!> - index 2 : counts the vertex
!> - the functions are defined on the interval [0,1][0,1]
#ifdef CUDA_KERNELS
  attributes(host,device) &
#endif
pure subroutine basisfunctions_2D_1_T(s, t, H, H_s, H_t)
implicit none
real*8, intent(in)  :: s          !< s-coordinate in the element [0,1]
real*8, intent(in)  :: t          !< t-coordinate in the element [0,1]
real*8, intent(out) :: H(4,4)     !< Basis functions
real*8, intent(out) :: H_s(4,4)   !< Basis functions derived with respect to s
real*8, intent(out) :: H_t(4,4)   !< Basis functions derived with respect to t

!---------------------------------------------------------- vertex (1)
H(1,1)   =(-1.d0 + s)**2*(1.d0 + 2.d0*s)*(-1.d0 + t)**2*(1.d0 + 2.d0*t)
H_s(1,1) =6.d0*(-1.d0 + s)*s*(-1.d0 + t)**2*(1.d0 + 2.d0*t)
H_t(1,1) =6.d0*(-1.d0 + s)**2*(1.d0 + 2.d0*s)*(-1.d0 + t)*t

H(2,1)   =3.d0*(-1.d0 + s)**2*s*(-1.d0 + t)**2*(1.d0 + 2.d0*t)
H_s(2,1) =3.d0*(-1.d0 + s)*(-1.d0 + 3.d0*s)*(-1.d0 + t)**2*(1.d0 + 2.d0*t)
H_t(2,1) =18.d0*(-1.d0 + s)**2*s*(-1.d0 + t)*t

H(3,1)   =3.d0*(-1.d0 + s)**2*(1.d0 + 2.d0*s)*(-1.d0 + t)**2*t
H_s(3,1) =18.d0*(-1.d0 + s)*s*(-1.d0 + t)**2*t
H_t(3,1) =3.d0*(-1.d0 + s)**2*(1.d0 + 2.d0*s)*(-1.d0 + t)*(-1.d0 + 3.d0*t)

H(4,1)   =9.d0*(-1.d0 + s)**2*s*(-1.d0 + t)**2*t
H_s(4,1) =9.d0*(-1.d0 + s)*(-1.d0 + 3.d0*s)*(-1.d0 + t)**2*t
H_t(4,1) =9.d0*(-1.d0 + s)**2*s*(-1.d0 + t)*(-1.d0 + 3.d0*t)

!---------------------------------------------------------- vertex (2)
H(1,2)   =-(s**2*(-3.d0 + 2.d0*s)*(-1.d0 + t)**2*(1.d0 + 2.d0*t))
H_s(1,2) =-6.d0*(-1.d0 + s)*s*(-1.d0 + t)**2*(1.d0 + 2.d0*t)
H_t(1,2) =-6.d0*s**2*(-3.d0 + 2.d0*s)*(-1.d0 + t)*t

H(2,2)   =-3.d0*(-1.d0 + s)*s**2*(-1.d0 + t)**2*(1.d0 + 2.d0*t)
H_s(2,2) =-3.d0*s*(-2.d0 + 3.d0*s)*(-1.d0 + t)**2*(1.d0 + 2.d0*t)
H_t(2,2) =-18.d0*(-1.d0 + s)*s**2*(-1.d0 + t)*t

H(3,2)   =-3.d0*s**2*(-3.d0 + 2.d0*s)*(-1.d0 + t)**2*t
H_s(3,2) =-18.d0*(-1.d0 + s)*s*(-1.d0 + t)**2*t
H_t(3,2) =3.d0*s**2*(-3.d0 + 2.d0*s)*(1.d0 - 3.d0*t)*(-1.d0 + t)

H(4,2)   =-9.d0*(-1.d0 + s)*s**2*(-1.d0 + t)**2*t
H_s(4,2) =-9.d0*s*(-2.d0 + 3.d0*s)*(-1.d0 + t)**2*t
H_t(4,2) =9.d0*(-1.d0 + s)*s**2*(1.d0 - 3.d0*t)*(-1.d0 + t)

!---------------------------------------------------------- vertex (3)
H(1,3)   =s**2*(-3.d0 + 2.d0*s)*t**2*(-3.d0 + 2.d0*t)
H_s(1,3) =6.d0*(-1.d0 + s)*s*t**2*(-3.d0 + 2.d0*t)
H_t(1,3) =6.d0*s**2*(-3.d0 + 2.d0*s)*(-1.d0 + t)*t

H(2,3)   =3.d0*(-1.d0 + s)*s**2*t**2*(-3.d0 + 2.d0*t)
H_s(2,3) =3.d0*s*(-2 + 3.d0*s)*t**2*(-3.d0 + 2.d0*t)
H_t(2,3) =18.d0*(-1.d0 + s)*s**2*(-1.d0 + t)*t

H(3,3)   =3.d0*s**2*(-3.d0 + 2.d0*s)*(-1.d0 + t)*t**2
H_s(3,3) =18.d0*(-1.d0 + s)*s*(-1.d0 + t)*t**2
H_t(3,3) =3.d0*s**2*(-3.d0 + 2.d0*s)*t*(-2.d0 + 3.d0*t)

H(4,3)   =9.d0*(-1.d0 + s)*s**2*(-1.d0 + t)*t**2
H_s(4,3) =9.d0*s*(-2.d0 + 3.d0*s)*(-1.d0 + t)*t**2
H_t(4,3) =9.d0*(-1.d0 + s)*s**2*t*(-2.d0 + 3.d0*t)

!---------------------------------------------------------- vertex (4)
H(1,4)   =-((-1.d0 + s)**2*(1.d0 + 2.d0*s)*t**2*(-3.d0 + 2.d0*t))
H_s(1,4) =-6.d0*(-1.d0 + s)*s*t**2*(-3.d0 + 2.d0*t)
H_t(1,4) =-6.d0*(-1.d0 + s)**2*(1.d0 + 2.d0*s)*(-1.d0 + t)*t

H(2,4)   =-3.d0*(-1.d0 + s)**2*s*t**2*(-3.d0 + 2.d0*t)
H_s(2,4) =3.d0*(1.d0 - 3*s)*(-1.d0 + s)*t**2*(-3.d0 + 2.d0*t)
H_t(2,4) =-18.d0*(-1.d0 + s)**2*s*(-1.d0 + t)*t

H(3,4)   =-3.d0*(-1.d0 + s)**2*(1.d0 + 2.d0*s)*(-1.d0 + t)*t**2
H_s(3,4) =-18.d0*(-1.d0 + s)*s*(-1.d0 + t)*t**2
H_t(3,4) =-3.d0*(-1.d0 + s)**2*(1.d0 + 2.d0*s)*t*(-2.d0 + 3.d0*t)

H(4,4)   =-9.d0*(-1.d0 + s)**2*s*(-1.d0 + t)*t**2
H_s(4,4) =9.d0*(1.d0 - 3.d0*s)*(-1.d0 + s)*(-1.d0 + t)*t**2
H_t(4,4) =-9.d0*(-1.d0 + s)**2*s*t*(-2.d0 + 3.d0*t)
end subroutine basisfunctions_2D_1_T


!> Basisfunctions in 2D with first and cross-derivative
pure subroutine basisfunctions_2D_1p(s,t,H,H_s,H_t,H_st)
implicit none
real*8, intent(in)  :: s          !< s-coordinate in the element [0,1]
real*8, intent(in)  :: t          !< t-coordinate in the element [0,1]
real*8, intent(out) :: H(4,4)     !< Basis functions
real*8, intent(out) :: H_s(4,4)   !< Basis functions derived with respect to s
real*8, intent(out) :: H_t(4,4)   !< Basis functions derived with respect to t
real*8, intent(out) :: H_st(4,4)  !< Basis functions derived with respect to s and t

! --- Local variables
real*8 :: scale1, dscale1_ds, dscale1_dt, dscale1_dsdt
real*8 :: scale2, dscale2_ds, dscale2_dt, dscale2_dsdt
real*8 :: scale3, dscale3_ds, dscale3_dt, dscale3_dsdt
real*8 :: scale4, dscale4_ds, dscale4_dt, dscale4_dsdt

real*8 :: h11, h12, h13, h14, h11_s, h12_s, h13_s, h14_s, h11_t, h12_t, h13_t, h14_t, h11_st, h12_st, h13_st, h14_st
real*8 :: h21, h22, h23, h24, h21_s, h22_s, h23_s, h24_s, h21_t, h22_t, h23_t, h24_t, h21_st, h22_st, h23_st, h24_st
real*8 :: h31, h32, h33, h34, h31_s, h32_s, h33_s, h34_s, h31_t, h32_t, h33_t, h34_t, h31_st, h32_st, h33_st, h34_st
real*8 :: h41, h42, h43, h44, h41_s, h42_s, h43_s, h44_s, h41_t, h42_t, h43_t, h44_t, h41_st, h42_st, h43_st, h44_st

scale1       =      (s-1.d0)**2 * (t-1.d0)**2
dscale1_ds   = 2.d0*(s-1.d0)    * (t-1.d0)**2
dscale1_dt   = 2.d0*(t-1.d0)    * (s-1.d0)**2
dscale1_dsdt = 4.d0*(s-1.d0)    * (t-1.d0)

h11 = (1.d0+2.d0*s)*(1.d0+2.d0*t); h11_s = 2.d0*(1.d0+2.d0*t); h11_t = 2.d0*(1.d0+2.d0*s); h11_st = 4.d0
h12 = 3.d0*s * (1.d0+2.d0*t)     ; h12_s = 3.d0*(1.d0+2.d0*t); h12_t = 6.d0*s            ; h12_st = 6.d0
h13 = 3.d0*t * (1.d0+2.d0*s)     ; h13_s = 6.d0*t            ; h13_t = 3.d0*(1.d0+2.d0*s); h13_st = 6.d0
h14 = 9.d0*s*t                   ; h14_s = 9.d0*t            ; h14_t = 9.d0*s            ; h14_st = 9.d0

H(1,1) = scale1 * h11
H(1,2) = scale1 * h12
H(1,3) = scale1 * h13
H(1,4) = scale1 * h14

H_s(1,1) = dscale1_ds * h11  + scale1 * h11_s
H_s(1,2) = dscale1_ds * h12  + scale1 * h12_s
H_s(1,3) = dscale1_ds * h13  + scale1 * h13_s
H_s(1,4) = dscale1_ds * h14  + scale1 * h14_s

H_t(1,1) = dscale1_dt * h11  + scale1 * h11_t
H_t(1,2) = dscale1_dt * h12  + scale1 * h12_t
H_t(1,3) = dscale1_dt * h13  + scale1 * h13_t
H_t(1,4) = dscale1_dt * h14  + scale1 * h14_t

H_st(1,1) = dscale1_dsdt * h11  + dscale1_dt * h11_s  +  dscale1_ds * h11_t  + scale1 * h11_st
H_st(1,2) = dscale1_dsdt * h12  + dscale1_dt * h12_s  +  dscale1_ds * h12_t  + scale1 * h12_st
H_st(1,3) = dscale1_dsdt * h13  + dscale1_dt * h13_s  +  dscale1_ds * h13_t  + scale1 * h13_st
H_st(1,4) = dscale1_dsdt * h14  + dscale1_dt * h14_s  +  dscale1_ds * h14_t  + scale1 * h14_st

!----------------------------------------------------------------------------
scale2       =       s**2 *      (t-1.d0)**2
dscale2_ds   =  2.d0*s    *      (t-1.d0)**2
dscale2_dt   =       s**2 * 2.d0*(t-1.d0)
dscale2_dsdt =  4.d0*s    *      (t-1.d0)

h21 = (3.d0-2.d0*s)*(1.d0+2.d0*t); h21_s = -2.d0*(1.d0+2.d0*t) ; h21_t =  2.d0*(3.d0-2.d0*s) ; h21_st = -4.d0
h22 = 3.d0*(1.d0-s)*(1.d0+2.d0*t); h22_s = -3.d0*(1.d0+2.d0*t) ; h22_t =  6.d0*(1.d0-s)      ; h22_st = -6.d0
h23 = 3.d0*(3.d0-2.d0*s)*t       ; h23_s = -6.d0*t             ; h23_t =  3.d0*(3.d0-2.d0*s) ; h23_st = -6.d0
h24 = 9.d0*(1.d0-s)* t           ; h24_s = -9.d0*t             ; h24_t =  9.d0*(1.d0-s)      ; h24_st = -9.d0

H(2,1) = scale2 * h21
H(2,2) = scale2 * h22
H(2,3) = scale2 * h23
H(2,4) = scale2 * h24

H_s(2,1) = dscale2_ds * h21 + scale2 * h21_s
H_s(2,2) = dscale2_ds * h22 + scale2 * h22_s
H_s(2,3) = dscale2_ds * h23 + scale2 * h23_s
H_s(2,4) = dscale2_ds * h24 + scale2 * h24_s

H_t(2,1) = dscale2_dt * h21 + scale2 * h21_t
H_t(2,2) = dscale2_dt * h22 + scale2 * h22_t
H_t(2,3) = dscale2_dt * h23 + scale2 * h23_t
H_t(2,4) = dscale2_dt * h24 + scale2 * h24_t

H_st(2,1) = dscale2_dsdt * h21 + dscale2_dt * h21_s +  dscale2_ds * h21_t + scale2 * h21_st
H_st(2,2) = dscale2_dsdt * h22 + dscale2_dt * h22_s +  dscale2_ds * h22_t + scale2 * h22_st
H_st(2,3) = dscale2_dsdt * h23 + dscale2_dt * h23_s +  dscale2_ds * h23_t + scale2 * h23_st
H_st(2,4) = dscale2_dsdt * h24 + dscale2_dt * h24_s +  dscale2_ds * h24_t + scale2 * h24_st


!----------------------------------------------------------------------------
scale3       =      s**2 *      t**2
dscale3_ds   = 2.d0*s    *      t**2
dscale3_dt   =      s**2 * 2.d0*t
dscale3_dsdt = 4.d0*s    *      t

h31 = (3.d0-2.d0*s)*(3.d0-2.d0*t) ; h31_s = -2.d0*(3.d0-2.d0*t) ;  h31_t = -2.d0*(3.d0-2.d0*s) ;  h31_st = +4.d0
h32 = 3.d0*(1.d0-s)*(3.d0-2.d0*t) ; h32_s = -3.d0*(3.d0-2.d0*t) ;  h32_t = -6.d0*(1.d0-s)      ;  h32_st = +6.d0
h33 = 3.d0*(3.d0-2.*s)*(1.d0-t)   ; h33_s = -6.d0*(1.d0-t)      ;  h33_t = -3.d0*(3.d0-2.d0*s) ;  h33_st = +6.d0
h34 = 9.d0*(1.d0-s)*(1.d0-t)      ; h34_s = -9.d0*(1.d0-t)      ;  h34_t = -9.d0*(1.d0-s)      ;  h34_st = +9.d0

H(3,1) = scale3 * h31
H(3,2) = scale3 * h32
H(3,3) = scale3 * h33
H(3,4) = scale3 * h34

H_s(3,1) = dscale3_ds * h31 + scale3 * h31_s
H_s(3,2) = dscale3_ds * h32 + scale3 * h32_s
H_s(3,3) = dscale3_ds * h33 + scale3 * h33_s
H_s(3,4) = dscale3_ds * h34 + scale3 * h34_s

H_t(3,1) = dscale3_dt * h31 + scale3 * h31_t
H_t(3,2) = dscale3_dt * h32 + scale3 * h32_t
H_t(3,3) = dscale3_dt * h33 + scale3 * h33_t
H_t(3,4) = dscale3_dt * h34 + scale3 * h34_t

H_st(3,1) = dscale3_dsdt * h31 + dscale3_dt * h31_s + dscale3_ds * h31_t + scale3 * h31_st
H_st(3,2) = dscale3_dsdt * h32 + dscale3_dt * h32_s + dscale3_ds * h32_t + scale3 * h32_st
H_st(3,3) = dscale3_dsdt * h33 + dscale3_dt * h33_s + dscale3_ds * h33_t + scale3 * h33_st
H_st(3,4) = dscale3_dsdt * h34 + dscale3_dt * h34_s + dscale3_ds * h34_t + scale3 * h34_st

!----------------------------------------------------------------------------
scale4       =       (s-1.d0)**2 *      t**2
dscale4_ds   =  2.d0*(s-1.d0)    *      t**2
dscale4_dt   =       (s-1.d0)**2 * 2.d0*t
dscale4_dsdt =  4.d0*(s-1.d0)    *      t

h41    =  (1.d0+2.d0*s)*(3.d0-2.d0*t);  h41_s = 2.d0*(3.d0-2.d0*t);  h41_t = -2.d0*(1.d0+2.d0*s) ;  h41_st =-4.d0
h42    = 3.d0*s*(3.d0-2.d0*t)        ;  h42_s = 3.d0*(3.d0-2.d0*t);  h42_t = -6.d0*s             ;  h42_st =-6.d0
h43    = 3.d0*(1.d0+2.d0*s)*(1.-t)   ;  h43_s = 6.d0*(1.d0-t)     ;  h43_t = -3.d0*(1.d0+2.d0*s) ;  h43_st =-6.d0
h44    = 9.d0*s*(1.d0-t)             ;  h44_s = 9.d0*(1.d0-t)     ;  h44_t = -9.d0*s             ;  h44_st =-9.d0

H(4,1) = scale4 * h41
H(4,2) = scale4 * h42
H(4,3) = scale4 * h43
H(4,4) = scale4 * h44

H_s(4,1) = dscale4_ds * h41 + scale4 *  h41_s
H_s(4,2) = dscale4_ds * h42 + scale4 *  h42_s
H_S(4,3) = dscale4_ds * h43 + scale4 *  h43_s
H_s(4,4) = dscale4_ds * h44 + scale4 *  h44_s

H_t(4,1) = dscale4_dt * h41  + scale4 * h41_t
H_t(4,2) = dscale4_dt * h42  + scale4 * h42_t
H_t(4,3) = dscale4_dt * h43  + scale4 * h43_t
H_t(4,4) = dscale4_dt * h44  + scale4 * h44_t

H_st(4,1) = dscale4_dsdt * h41 + dscale4_dt * h41_s +  dscale4_ds * h41_t + scale4 * h41_st
H_st(4,2) = dscale4_dsdt * h42 + dscale4_dt * h42_s +  dscale4_ds * h42_t + scale4 * h42_st
H_st(4,3) = dscale4_dsdt * h43 + dscale4_dt * h43_s +  dscale4_ds * h43_t + scale4 * h43_st
H_st(4,4) = dscale4_dsdt * h44 + dscale4_dt * h44_s +  dscale4_ds * h44_t + scale4 * h44_st
end subroutine basisfunctions_2D_1p


!> Basisfunctions with second derivatives in 2D
pure subroutine basisfunctions_2D_2(s, t, H, H_s, H_t, H_st, H_ss, H_tt)
implicit none
real*8, intent(in)  :: s          !< s-coordinate in the element [0,1]
real*8, intent(in)  :: t          !< t-coordinate in the element [0,1]
real*8, intent(out) :: H(4,4)     !< Basis functions
real*8, intent(out) :: H_s(4,4)   !< Basis functions derived with respect to s
real*8, intent(out) :: H_t(4,4)   !< Basis functions derived with respect to t
real*8, intent(out) :: H_st(4,4)  !< Basis functions derived with respect to s and t
real*8, intent(out) :: H_ss(4,4)  !< Basis functions derived two times with respect to s
real*8, intent(out) :: H_tt(4,4)  !< Basis functions derived two times with respect to t

call basisfunctions_2D_1(s, t, H, H_s, H_t)

!---------------------------------------------------------- vertex (1)
H_st(1,1)=36.d0*(-1.d0 + s)*s*(-1.d0 + t)*t
H_ss(1,1)=6.d0*(-1.d0 + 2.d0*s)*(-1.d0 + t)**2*(1.d0 + 2.d0*t)
H_tt(1,1)=6.d0*(-1.d0 + s)**2*(1.d0 + 2.d0*s)*(-1.d0 + 2.d0*t)

H_st(1,2)=18.d0*(1.d0 - 4.d0*s + 3.d0*s**2)*(-1.d0 + t)*t
H_ss(1,2)=6.d0*(-2.d0 + 3.d0*s)*(-1.d0 + t)**2*(1.d0 + 2.d0*t)
H_tt(1,2)=18.d0*(-1.d0 + s)**2*s*(-1.d0 + 2.d0*t)

H_st(1,3)=18.d0*(-1.d0 + s)*s*(1.d0 - 4.d0*t + 3.d0*t**2)
H_ss(1,3)=18.d0*(-1.d0 + 2.d0*s)*(-1.d0 + t)**2*t
H_tt(1,3)=6.d0*(-1.d0 + s)**2*(1.d0 + 2.d0*s)*(-2.d0 + 3.d0*t)

H_st(1,4)=9.d0*(1.d0 - 4.d0*s + 3.d0*s**2)*(1.d0 - 4.d0*t + 3.d0*t**2)
H_ss(1,4)=18.d0*(-2.d0 + 3.d0*s)*(-1.d0 + t)**2*t
H_tt(1,4)=18.d0*(-1.d0 + s)**2*s*(-2.d0 + 3.d0*t)

!---------------------------------------------------------- vertex (2)
H_st(2,1)=-36.d0*(-1.d0 + s)*s*(-1.d0 + t)*t
H_ss(2,1)=-6.d0*(-1.d0 + 2.d0*s)*(-1.d0 + t)**2*(1.d0 + 2.d0*t)
H_tt(2,1)=-6.d0*s**2*(-3.d0 + 2.d0*s)*(-1.d0 + 2.d0*t)

H_st(2,2)=-18.d0*s*(-2.d0 + 3.d0*s)*(-1.d0 + t)*t
H_ss(2,2)=-6.d0*(-1.d0 + 3.d0*s)*(-1.d0 + t)**2*(1 + 2.d0*t)
H_tt(2,2)=-18.d0*(-1.d0 + s)*s**2*(-1.d0 + 2.d0*t)

H_st(2,3)=-18.d0*(-1.d0 + s)*s*(1.d0 - 4.d0*t + 3.d0*t**2)
H_ss(2,3)=6.d0*(3.d0 - 6.d0*s)*(-1.d0 + t)**2*t
H_tt(2,3)=-6.d0*s**2*(-3.d0 + 2.d0*s)*(-2.d0 + 3.d0*t)

H_st(2,4)=-9.d0*s*(-2.d0 + 3.d0*s)*(1.d0 - 4.d0*t + 3.d0*t**2)
H_ss(2,4)=18.d0*(1.d0 - 3.d0*s)*(-1.d0 + t)**2*t
H_tt(2,4)=-18.d0*(-1.d0 + s)*s**2*(-2.d0 + 3.d0*t)

!---------------------------------------------------------- vertex (3)
H_st(3,1)=36.d0*(-1.d0 + s)*s*(-1.d0 + t)*t
H_ss(3,1)=6.d0*(-1.d0 + 2*s)*t**2*(-3.d0 + 2.d0*t)
H_tt(3,1)=6.d0*s**2*(-3.d0 + 2.d0*s)*(-1.d0 + 2.d0*t)

H_st(3,2)=18.d0*s*(-2.d0 + 3.d0*s)*(-1.d0 + t)*t
H_ss(3,2)=6.d0*(-1.d0 + 3.d0*s)*t**2*(-3.d0 + 2.d0*t)
H_tt(3,2)=18.d0*(-1.d0 + s)*s**2*(-1.d0 + 2.d0*t)

H_st(3,3)=18.d0*(-1.d0 + s)*s*t*(-2.d0 + 3.d0*t)
H_ss(3,3)=18.d0*(-1.d0 + 2.d0*s)*(-1.d0 + t)*t**2
H_tt(3,3)=6.d0*s**2*(-3.d0 + 2.d0*s)*(-1.d0 + 3.d0*t)

H_st(3,4)=9.d0*s*(-2.d0 + 3.d0*s)*t*(-2.d0 + 3.d0*t)
H_ss(3,4)=18.d0*(-1.d0 + 3.d0*s)*(-1.d0 + t)*t**2
H_tt(3,4)=18.d0*(-1.d0 + s)*s**2*(-1.d0 + 3.d0*t)

!---------------------------------------------------------- vertex (4)
H_st(4,1)=-36.d0*(-1.d0 + s)*s*(-1.d0 + t)*t
H_ss(4,1)=-6.d0*(-1.d0 + 2.d0*s)*t**2*(-3.d0 + 2.d0*t)
H_tt(4,1)=-6.d0*(-1.d0 + s)**2*(1.d0 + 2.d0*s)*(-1.d0 + 2.d0*t)

H_st(4,2)=-18.d0*(1.d0 - 4.d0*s + 3.d0*s**2)*(-1.d0 + t)*t
H_ss(4,2)=-6.d0*(-2.d0 + 3.d0*s)*t**2*(-3.d0 + 2.d0*t)
H_tt(4,2)=6.d0*(-1.d0 + s)**2*s*(3.d0 - 6.d0*t)

H_st(4,3)=-18.d0*(-1.d0 + s)*s*t*(-2.d0 + 3.d0*t)
H_ss(4,3)=-18.d0*(-1.d0 + 2.d0*s)*(-1.d0 + t)*t**2
H_tt(4,3)=-6.d0*(-1.d0 + s)**2*(1.d0 + 2.d0*s)*(-1 + 3.d0*t)

H_st(4,4)=-9.d0*(1.d0 - 4.d0*s + 3.d0*s**2)*t*(-2.d0 + 3.d0*t)
H_ss(4,4)=-18.d0*(-2.d0 + 3.d0*s)*(-1.d0 + t)*t**2
H_tt(4,4)=18.d0*(-1.d0 + s)**2*s*(1.d0 - 3.d0*t)
end subroutine basisfunctions_2D_2

!> Basisfunctions with second derivatives in 2D and transposed matrix for 
!> better vectorisation
pure subroutine basisfunctions_2D_2_T(s, t, H, H_s, H_t, H_st, H_ss, H_tt)
implicit none
real*8, intent(in)  :: s          !< s-coordinate in the element [0,1]
real*8, intent(in)  :: t          !< t-coordinate in the element [0,1]
real*8, intent(out) :: H(4,4)     !< Basis functions
real*8, intent(out) :: H_s(4,4)   !< Basis functions derived with respect to s
real*8, intent(out) :: H_t(4,4)   !< Basis functions derived with respect to t
real*8, intent(out) :: H_st(4,4)  !< Basis functions derived with respect to s and t
real*8, intent(out) :: H_ss(4,4)  !< Basis functions derived two times with respect to s
real*8, intent(out) :: H_tt(4,4)  !< Basis functions derived two times with respect to t

call basisfunctions_2D_1_T(s, t, H, H_s, H_t)

!---------------------------------------------------------- vertex (1)
H_st(1,1)=36.d0*(-1.d0 + s)*s*(-1.d0 + t)*t
H_ss(1,1)=6.d0*(-1.d0 + 2.d0*s)*(-1.d0 + t)**2*(1.d0 + 2.d0*t)
H_tt(1,1)=6.d0*(-1.d0 + s)**2*(1.d0 + 2.d0*s)*(-1.d0 + 2.d0*t)

H_st(2,1)=18.d0*(1.d0 - 4.d0*s + 3.d0*s**2)*(-1.d0 + t)*t
H_ss(2,1)=6.d0*(-2.d0 + 3.d0*s)*(-1.d0 + t)**2*(1.d0 + 2.d0*t)
H_tt(2,1)=18.d0*(-1.d0 + s)**2*s*(-1.d0 + 2.d0*t)

H_st(3,1)=18.d0*(-1.d0 + s)*s*(1.d0 - 4.d0*t + 3.d0*t**2)
H_ss(3,1)=18.d0*(-1.d0 + 2.d0*s)*(-1.d0 + t)**2*t
H_tt(3,1)=6.d0*(-1.d0 + s)**2*(1.d0 + 2.d0*s)*(-2.d0 + 3.d0*t)

H_st(4,1)=9.d0*(1.d0 - 4.d0*s + 3.d0*s**2)*(1.d0 - 4.d0*t + 3.d0*t**2)
H_ss(4,1)=18.d0*(-2.d0 + 3.d0*s)*(-1.d0 + t)**2*t
H_tt(4,1)=18.d0*(-1.d0 + s)**2*s*(-2.d0 + 3.d0*t)

!---------------------------------------------------------- vertex (2)
H_st(1,2)=-36.d0*(-1.d0 + s)*s*(-1.d0 + t)*t
H_ss(1,2)=-6.d0*(-1.d0 + 2.d0*s)*(-1.d0 + t)**2*(1.d0 + 2.d0*t)
H_tt(1,2)=-6.d0*s**2*(-3.d0 + 2.d0*s)*(-1.d0 + 2.d0*t)

H_st(2,2)=-18.d0*s*(-2.d0 + 3.d0*s)*(-1.d0 + t)*t
H_ss(2,2)=-6.d0*(-1.d0 + 3.d0*s)*(-1.d0 + t)**2*(1 + 2.d0*t)
H_tt(2,2)=-18.d0*(-1.d0 + s)*s**2*(-1.d0 + 2.d0*t)

H_st(3,2)=-18.d0*(-1.d0 + s)*s*(1.d0 - 4.d0*t + 3.d0*t**2)
H_ss(3,2)=6.d0*(3.d0 - 6.d0*s)*(-1.d0 + t)**2*t
H_tt(3,2)=-6.d0*s**2*(-3.d0 + 2.d0*s)*(-2.d0 + 3.d0*t)

H_st(4,2)=-9.d0*s*(-2.d0 + 3.d0*s)*(1.d0 - 4.d0*t + 3.d0*t**2)
H_ss(4,2)=18.d0*(1.d0 - 3.d0*s)*(-1.d0 + t)**2*t
H_tt(4,2)=-18.d0*(-1.d0 + s)*s**2*(-2.d0 + 3.d0*t)

!---------------------------------------------------------- vertex (3)
H_st(1,3)=36.d0*(-1.d0 + s)*s*(-1.d0 + t)*t
H_ss(1,3)=6.d0*(-1.d0 + 2*s)*t**2*(-3.d0 + 2.d0*t)
H_tt(1,3)=6.d0*s**2*(-3.d0 + 2.d0*s)*(-1.d0 + 2.d0*t)

H_st(2,3)=18.d0*s*(-2.d0 + 3.d0*s)*(-1.d0 + t)*t
H_ss(2,3)=6.d0*(-1.d0 + 3.d0*s)*t**2*(-3.d0 + 2.d0*t)
H_tt(2,3)=18.d0*(-1.d0 + s)*s**2*(-1.d0 + 2.d0*t)

H_st(3,3)=18.d0*(-1.d0 + s)*s*t*(-2.d0 + 3.d0*t)
H_ss(3,3)=18.d0*(-1.d0 + 2.d0*s)*(-1.d0 + t)*t**2
H_tt(3,3)=6.d0*s**2*(-3.d0 + 2.d0*s)*(-1.d0 + 3.d0*t)

H_st(4,3)=9.d0*s*(-2.d0 + 3.d0*s)*t*(-2.d0 + 3.d0*t)
H_ss(4,3)=18.d0*(-1.d0 + 3.d0*s)*(-1.d0 + t)*t**2
H_tt(4,3)=18.d0*(-1.d0 + s)*s**2*(-1.d0 + 3.d0*t)

!---------------------------------------------------------- vertex (4)
H_st(1,4)=-36.d0*(-1.d0 + s)*s*(-1.d0 + t)*t
H_ss(1,4)=-6.d0*(-1.d0 + 2.d0*s)*t**2*(-3.d0 + 2.d0*t)
H_tt(1,4)=-6.d0*(-1.d0 + s)**2*(1.d0 + 2.d0*s)*(-1.d0 + 2.d0*t)

H_st(2,4)=-18.d0*(1.d0 - 4.d0*s + 3.d0*s**2)*(-1.d0 + t)*t
H_ss(2,4)=-6.d0*(-2.d0 + 3.d0*s)*t**2*(-3.d0 + 2.d0*t)
H_tt(2,4)=6.d0*(-1.d0 + s)**2*s*(3.d0 - 6.d0*t)

H_st(3,4)=-18.d0*(-1.d0 + s)*s*t*(-2.d0 + 3.d0*t)
H_ss(3,4)=-18.d0*(-1.d0 + 2.d0*s)*(-1.d0 + t)*t**2
H_tt(3,4)=-6.d0*(-1.d0 + s)**2*(1.d0 + 2.d0*s)*(-1 + 3.d0*t)

H_st(4,4)=-9.d0*(1.d0 - 4.d0*s + 3.d0*s**2)*t*(-2.d0 + 3.d0*t)
H_ss(4,4)=-18.d0*(-2.d0 + 3.d0*s)*(-1.d0 + t)*t**2
H_tt(4,4)=18.d0*(-1.d0 + s)**2*s*(1.d0 - 3.d0*t)
end subroutine basisfunctions_2D_2_T

pure subroutine basisfunctions3(s, t, H, H_s, H_t)
implicit none

! --- Routine parameters
real*8, intent(in)  :: s          !< s-coordinate in the element
real*8, intent(in)  :: t          !< t-coordinate in the element
real*8, intent(out) :: H(4,4)     !< Basis functions
real*8, intent(out) :: H_s(4,4)   !< Basis functions derived with respect to s
real*8, intent(out) :: H_t(4,4)   !< Basis functions derived with respect to t

!---------------------------------------------------------- vertex (1)
H(1,1)   =(-1.d0 + s)**2*(1.d0 + 2.d0*s)*(-1.d0 + t)**2*(1.d0 + 2.d0*t)
H_s(1,1) =6.d0*(-1.d0 + s)*s*(-1.d0 + t)**2*(1.d0 + 2.d0*t)
H_t(1,1) =6.d0*(-1.d0 + s)**2*(1.d0 + 2.d0*s)*(-1.d0 + t)*t

H(1,2)   =3.d0*(-1.d0 + s)**2*s*(-1.d0 + t)**2*(1.d0 + 2.d0*t)
H_s(1,2) =3.d0*(-1.d0 + s)*(-1.d0 + 3.d0*s)*(-1.d0 + t)**2*(1.d0 + 2.d0*t)
H_t(1,2) =18.d0*(-1.d0 + s)**2*s*(-1.d0 + t)*t

H(1,3)   =3.d0*(-1.d0 + s)**2*(1.d0 + 2.d0*s)*(-1.d0 + t)**2*t
H_s(1,3) =18.d0*(-1.d0 + s)*s*(-1.d0 + t)**2*t
H_t(1,3) =3.d0*(-1.d0 + s)**2*(1.d0 + 2.d0*s)*(-1.d0 + t)*(-1.d0 + 3.d0*t)

H(1,4)   =9.d0*(-1.d0 + s)**2*s*(-1.d0 + t)**2*t
H_s(1,4) =9.d0*(-1.d0 + s)*(-1.d0 + 3.d0*s)*(-1.d0 + t)**2*t
H_t(1,4) =9.d0*(-1.d0 + s)**2*s*(-1.d0 + t)*(-1.d0 + 3.d0*t)

!---------------------------------------------------------- vertex (2)
H(2,1)   =-(s**2*(-3.d0 + 2.d0*s)*(-1.d0 + t)**2*(1.d0 + 2.d0*t))
H_s(2,1) =-6.d0*(-1.d0 + s)*s*(-1.d0 + t)**2*(1.d0 + 2.d0*t)
H_t(2,1) =-6.d0*s**2*(-3.d0 + 2.d0*s)*(-1.d0 + t)*t

H(2,2)   =-3.d0*(-1.d0 + s)*s**2*(-1.d0 + t)**2*(1.d0 + 2.d0*t)
H_s(2,2) =-3.d0*s*(-2.d0 + 3.d0*s)*(-1.d0 + t)**2*(1.d0 + 2.d0*t)
H_t(2,2) =-18.d0*(-1.d0 + s)*s**2*(-1.d0 + t)*t

H(2,3)   =-3.d0*s**2*(-3.d0 + 2.d0*s)*(-1.d0 + t)**2*t
H_s(2,3) =-18.d0*(-1.d0 + s)*s*(-1.d0 + t)**2*t
H_t(2,3) =3.d0*s**2*(-3.d0 + 2.d0*s)*(1.d0 - 3.d0*t)*(-1.d0 + t)

H(2,4)   =-9.d0*(-1.d0 + s)*s**2*(-1.d0 + t)**2*t
H_s(2,4) =-9.d0*s*(-2.d0 + 3.d0*s)*(-1.d0 + t)**2*t
H_t(2,4) =9.d0*(-1.d0 + s)*s**2*(1.d0 - 3.d0*t)*(-1.d0 + t)

!---------------------------------------------------------- vertex (3)
H(3,1)   =s**2*(-3.d0 + 2.d0*s)*t**2*(-3.d0 + 2.d0*t)
H_s(3,1) =6.d0*(-1.d0 + s)*s*t**2*(-3.d0 + 2.d0*t)
H_t(3,1) =6.d0*s**2*(-3.d0 + 2.d0*s)*(-1.d0 + t)*t

H(3,2)   =3.d0*(-1.d0 + s)*s**2*t**2*(-3.d0 + 2.d0*t)
H_s(3,2) =3.d0*s*(-2 + 3.d0*s)*t**2*(-3.d0 + 2.d0*t)
H_t(3,2) =18.d0*(-1.d0 + s)*s**2*(-1.d0 + t)*t

H(3,3)   =3.d0*s**2*(-3.d0 + 2.d0*s)*(-1.d0 + t)*t**2
H_s(3,3) =18.d0*(-1.d0 + s)*s*(-1.d0 + t)*t**2
H_t(3,3) =3.d0*s**2*(-3.d0 + 2.d0*s)*t*(-2.d0 + 3.d0*t)

H(3,4)   =9.d0*(-1.d0 + s)*s**2*(-1.d0 + t)*t**2
H_s(3,4) =9.d0*s*(-2.d0 + 3.d0*s)*(-1.d0 + t)*t**2
H_t(3,4) =9.d0*(-1.d0 + s)*s**2*t*(-2.d0 + 3.d0*t)

!---------------------------------------------------------- vertex (4)
H(4,1)   =-((-1.d0 + s)**2*(1.d0 + 2.d0*s)*t**2*(-3.d0 + 2.d0*t))
H_s(4,1) =-6.d0*(-1.d0 + s)*s*t**2*(-3.d0 + 2.d0*t)
H_t(4,1) =-6.d0*(-1.d0 + s)**2*(1.d0 + 2.d0*s)*(-1.d0 + t)*t

H(4,2)   =-3.d0*(-1.d0 + s)**2*s*t**2*(-3.d0 + 2.d0*t)
H_s(4,2) =3.d0*(1.d0 - 3*s)*(-1.d0 + s)*t**2*(-3.d0 + 2.d0*t)
H_t(4,2) =-18.d0*(-1.d0 + s)**2*s*(-1.d0 + t)*t

H(4,3)   =-3.d0*(-1.d0 + s)**2*(1.d0 + 2.d0*s)*(-1.d0 + t)*t**2
H_s(4,3) =-18.d0*(-1.d0 + s)*s*(-1.d0 + t)*t**2
H_t(4,3) =-3.d0*(-1.d0 + s)**2*(1.d0 + 2.d0*s)*t*(-2.d0 + 3.d0*t)

H(4,4)   =-9.d0*(-1.d0 + s)**2*s*(-1.d0 + t)*t**2
H_s(4,4) =9.d0*(1.d0 - 3.d0*s)*(-1.d0 + s)*(-1.d0 + t)*t**2
H_t(4,4) =-9.d0*(-1.d0 + s)**2*s*t*(-2.d0 + 3.d0*t)
end subroutine basisfunctions3

end module mod_basisfunctions
