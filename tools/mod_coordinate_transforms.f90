!> Module to calculate coordinate transforms between XYZ and RZPhi coordinate systems.
module mod_coordinate_transforms
  implicit none
  private
  public cartesian_to_cylindrical
  public cylindrical_to_cartesian
  public vector_cartesian_to_cylindrical
  public vector_cylindrical_to_cartesian
  public vector_rotation
  public cartesian_velocity_to_cylindrical
  public transform_derivatives_st_to_RZ, transform_derivatives_RZ_to_st
  public transform_first_derivatives_st_to_RZ
  public transform_second_derivatives_st_to_RZ
  !> interfaces
  interface cylindrical_to_cartesian
    module procedure cylindrical_to_cartesian_real8
    module procedure cylindrical_to_cartesian_real4
  end interface cylindrical_to_cartesian

  !> overload transform derivatives from local to global coordinates
  interface transform_derivatives_st_to_RZ
     module procedure transform_first_derivatives_st_to_RZ, &
       transform_second_derivatives_st_to_RZ
  end interface transform_derivatives_st_to_RZ

  !> overload transfrom derivatives from global to local coorsinates
  interface transform_derivatives_RZ_to_st
     module procedure transform_first_derivatives_RZ_to_st, &
       transform_second_derivatives_RZ_to_st
  end interface transform_derivatives_RZ_to_st

contains
  !> convert a position in xyz coordinates to RZPhi coordinates
  pure function cartesian_to_cylindrical(xyz) result(cyl)
    real*8, intent(in)           :: xyz(3) !< The position in xyz coordinates
    real*8                       :: cyl(3) !< The position in RZPhi coordinates

    cyl(1) = sqrt(xyz(1)*xyz(1) + xyz(2)*xyz(2))
    cyl(2) = xyz(3)
    cyl(3) = atan2(-xyz(2), xyz(1))
  end function cartesian_to_cylindrical

 !> converts a position in RZPhi coordinates to xyz coordinates
  pure function cylindrical_to_cartesian_real8(cyl) result(xyz)
    real*8, intent(in)           :: cyl(3) !< The vector components in RZPhi coordinates
    real*8                       :: xyz(3) !< The vector components in xyz coordinates

    xyz(1) = cyl(1)*cos(-cyl(3))
    xyz(2) = cyl(1)*sin(-cyl(3))
    xyz(3) = cyl(2)
  end function cylindrical_to_cartesian_real8

  !> converts a position in RZPhi coordinates to xyz coordinates
  pure function cylindrical_to_cartesian_real4(cyl) result(xyz)
    real*4, intent(in)           :: cyl(3) !< The vector components in RZPhi coordinates
    real*4                       :: xyz(3) !< The vector components in xyz coordinates

    xyz(1) = cyl(1)*cos(-cyl(3))
    xyz(2) = cyl(1)*sin(-cyl(3))
    xyz(3) = cyl(2)
  end function cylindrical_to_cartesian_real4

  pure function vector_cartesian_to_cylindrical(phi,a) result(b)
    real*8, intent(in)               :: phi !< The local toroidal angle
    real*8, dimension(3), intent(in) :: a   !< The vector components in (ex,ey,ez) basis
    real*8, dimension(3)             :: b   !< The vector components in (eR,eZ,ephi) basis
    real*8, dimension(2)             :: sincosphi

    sincosphi = (/sin(phi),cos(phi)/)
   
    b(1) = a(1)*sincosphi(2) - a(2)*sincosphi(1) 
    b(2) = a(3)
    b(3) = -1.d0*(a(1)*sincosphi(1) + a(2)*sincosphi(2))
  end function vector_cartesian_to_cylindrical  

  !> convert a vector in (eR,eZ,ephi) basis into (ex,ey,ez)  basis
  pure function vector_cylindrical_to_cartesian(phi,a) result(b)
    real*8, intent(in)               :: phi !< The local toroidal angle
    real*8, dimension(3), intent(in) :: a   !< The vector components in (eR,eZ,ephi) basis
    real*8, dimension(3)             :: b   !< The vector components in (ex,ey,ez) basis
    real*8, dimension(2)             :: sincosphi

    sincosphi = (/sin(phi),cos(phi)/)

    b(1) = a(1)*sincosphi(2) - a(3)*sincosphi(1) 
    b(2) = -1.d0*(a(1)*sincosphi(1) + a(3)*sincosphi(2))
    b(3) = a(2)
  end function vector_cylindrical_to_cartesian

  !> multiply a vector with the rotation matrix for angle phi around the z-axis in RZPhi coordinates
  pure function vector_rotation(in, phi) result(out)
    real*8, intent(in) :: in(3) !< Input vector in RZPhi coordinates
    real*8, intent(in) :: phi
    real*8             :: out(3) !< Output vector in RZPhi coordinates

    out(1) = cos(phi) * in(1) - sin(phi) * in(3)
    out(2) = in(2)
    out(3) = sin(phi) * in(1) + cos(phi) * in(3)
  end function vector_rotation

 !> Calculate the corresponding cylindrical expression for a cartesian velocity vector
  pure function cartesian_velocity_to_cylindrical(in, phi) result(out)
    real*8, intent(in) :: in(3) !< Input velocity in xyz coordinates
    real*8, intent(in) :: phi !< toroidal angle of particle (jorek coordinate 3)
    real*8             :: out(3) !< Output velocity in RZPhi coordinates

    out(1) = cos(phi) * in(1) + sin(phi) * in(2)
    out(2) = in(3)
    out(3) = -sin(phi) * in(1) + cos(phi) * in(2)
  end function cartesian_velocity_to_cylindrical

!--------------------------------------------------------------------------
!> This procedure expresses first order derivatives from the local (s,t)
!> to the global (R,Z) reference system.
!> inputs:
!>   n_v: (integer) number of physical quantities
!>   P_s: (real8)(n_v) physical quantities first derivative in s
!>   P_t: (real8)(n_v) physical quantities first derivatives in t
!>   R_s: (real8) major radius first derivative in s
!>   R_z: (real8) major radius first derivative in t
!>   Z_s: (real8) vertical position first derivative in s
!>   Z_t: (real8) vertical position first derivative in t
!> outputs:
!>   P_R: (real8)(n_v) physical quantities first derivative in R
!>   P_Z: (real8)(n_v) physcial quantities first derivative in Z
pure subroutine transform_first_derivatives_st_to_RZ(P_R,P_Z,n_v, &
  P_s,P_t,R_s,R_t,Z_s,Z_t)
  implicit none
  !> input variables
  integer, intent(in)                      :: n_v
  real(kind=8), dimension(n_v), intent(in) :: P_s, P_t
  real(kind=8), intent(in)                 :: R_s, R_t, Z_s, Z_t
  !> output variables
  real(kind=8), dimension(n_v), intent(out) :: P_R, P_Z
  !> internal variables
  real(kind=8) :: inverse_jacobian

  !> compute the inverse of the jacobian
  inverse_jacobian = 1.d0/(R_s*Z_t-R_t*Z_s)

  !> compute derivatives wrt. R and Z
  P_R = (P_s*Z_t-P_t*Z_s)*inverse_jacobian
  P_Z = (P_t*R_s-P_s*R_t)*inverse_jacobian
  
end subroutine transform_first_derivatives_st_to_RZ

!---------------------------------------------------------------------------

!> This procedure expresses second order derivatives from the local (s,t)
!> to the global (R,Z) reference system. Functions are assumed to be Hessian.
!> inputs:
!>   n_v:  (integer) number of physical quantities
!>   P_ss: (real8)(n_v) physical quantities second derivatives in s
!>   P_st: (real8)(n_v) physical quantities cross derivatives in s,t
!>   P_tt: (real8)(n_v) physical quantities second derivative in t
!>   P_R:  (real8)(n_v) physical quantities first derivative in R
!>   P_Z:  (real8)(n_v) physical quantities first derivative in Z
!>   R_s:  (real8) major radius first derivative in s
!>   R_t:  (real8) major radius first derivative in t
!>   R_ss: (real8) major radius second derivative in s
!>   R_st: (real8) major radius cross derivative in s,t
!>   R_tt: (real8) major radius second derivative in t
!>   Z_s:  (real8) vertical position first derivative in s
!>   Z_t:  (real8) vertical position first derivative in t
!>   Z_ss: (real8) vertical position second derivative in s
!>   Z_st: (real8) vertical position cross derivative in s,t
!>   Z_tt: (real8) vertical position second derivative in t
!> outputs:
!>   P_RR: (real8)(n_v) physical quantities second derivatives in R
!>   P_RZ: (real8)(n_v) physical quantities cross derivatives in R,Z
!>   P_ZZ: (real8)(n_v) physical quantities second derivatives in Z
pure subroutine transform_second_derivatives_st_to_RZ(P_RR,P_RZ,P_ZZ,n_v, &
  P_ss,P_st,P_tt,P_R,P_Z,R_s,R_t,R_ss,R_st,R_tt,Z_s,Z_t,Z_ss,Z_st,Z_tt)
  implicit none
  !> input variables
  integer, intent(in) :: n_v
  real(kind=8), dimension(n_v), intent(in) :: P_ss, P_st, P_tt, P_R, P_Z
  real(kind=8), intent(in) :: R_s, R_t, R_ss, R_st, R_tt
  real(kind=8), intent(in) :: Z_s, Z_t, Z_ss, Z_st, Z_tt
  !> output variables
  real(kind=8), dimension(n_v), intent(out) :: P_RR, P_RZ, P_ZZ
  !> internal variables
  integer                      :: i
  real(kind=8), dimension(10)  :: transformation_matrix !< 10:jacobian
  real(kind=8), dimension(n_v) :: RHS_RR, RHS_RZ, RHS_ZZ

  !< compute matrix elements
  transformation_matrix(1:9) = [R_s*R_s,2.d0*R_s*Z_s,Z_s*Z_s, &
     R_t*R_s,R_s*Z_t+Z_s*R_t,Z_s*Z_t,R_t*R_t,                 &
     2.d0*R_t*Z_t,Z_t*Z_t]

  !< compute the inverse of the matrix jacobian
  transformation_matrix(10) = 1.d0/(transformation_matrix(1)*( &
    transformation_matrix(5)*transformation_matrix(9)-         &
    transformation_matrix(8)*transformation_matrix(6))+        &
    transformation_matrix(2)*(transformation_matrix(7)*        &
    transformation_matrix(6)-transformation_matrix(4)*         &
    transformation_matrix(9))+transformation_matrix(3)*(       &
    transformation_matrix(4)*transformation_matrix(8)-         &
    transformation_matrix(5)*transformation_matrix(7)))
    
  !> compute RHS
  RHS_RR = P_ss - R_ss*P_R - Z_ss*P_Z
  RHS_RZ = P_st - R_st*P_R - Z_st*P_Z
  RHS_ZZ = P_tt - R_tt*P_R - Z_tt*P_Z
  
  !> compute second order derivatives
  P_RR = ((transformation_matrix(5)*transformation_matrix(9)-   &
    transformation_matrix(8)*transformation_matrix(6))*RHS_RR + &
    (transformation_matrix(8)*transformation_matrix(3)-         &
    transformation_matrix(2)*transformation_matrix(9))*RHS_RZ + &
    (transformation_matrix(2)*transformation_matrix(6)-         &
    transformation_matrix(5)*transformation_matrix(3))*         &
    RHS_ZZ)*transformation_matrix(10)
  P_RZ = ((transformation_matrix(7)*transformation_matrix(6)-   &
    transformation_matrix(4)*transformation_matrix(9))*RHS_RR + &
    (transformation_matrix(1)*transformation_matrix(9)-         &
    transformation_matrix(7)*transformation_matrix(3))*RHS_RZ + &
    (transformation_matrix(4)*transformation_matrix(3)-         &
    transformation_matrix(1)*transformation_matrix(6))*         &
    RHS_ZZ)*transformation_matrix(10)
  P_ZZ = ((transformation_matrix(4)*transformation_matrix(8)-   &
    transformation_matrix(5)*transformation_matrix(7))*RHS_RR+  &
    (transformation_matrix(2)*transformation_matrix(7)-         &
    transformation_matrix(1)*transformation_matrix(8))*RHS_RZ+  &
    (transformation_matrix(1)*transformation_matrix(5)-         &
    transformation_matrix(2)*transformation_matrix(4))*         &
    RHS_ZZ)*transformation_matrix(10) 
  
end subroutine transform_second_derivatives_st_to_RZ

!---------------------------------------------------------------------------------------

!> This procedure transforms back the first order derivatives from global RZ coordinates
!> to local s,t coordinates. This method is used mainly for tests.
!> inputs:
!>   n_v:  (integer) number of physical quantities
!>   P_R,P_Z: (real8)(n_v) physical quantities first derivatives in R and Z
!>   R_s,R_t: (real8) major radius first derivatives in s and t
!>   Z_s,Z_t: (real) vertical position first derivatives in s and t
!> outputs:
!>   P_s,P_t: (real8)(n_v) physical quantities first derivatives in s and t
pure subroutine transform_first_derivatives_RZ_to_st(P_s,P_t,n_v,P_R,P_Z, &
  R_s,R_t,Z_s,Z_t)
  implicit none
  !> input variables
  integer, intent(in)                      :: n_v
  real(kind=8), dimension(n_v), intent(in) :: P_R, P_Z
  real(kind=8), intent(in)                 :: R_s, R_t, Z_s, Z_t
  !> output varibales
  real(kind=8), dimension(n_v), intent(out) :: P_s, P_t

  P_s = P_R*R_s + P_Z*Z_s
  P_t = P_R*R_t + P_Z*Z_t
  
end subroutine transform_first_derivatives_RZ_to_st
!---------------------------------------------------------------------------------------

!> This procedure transforms back the second order derivatives from global RZ coordinates
!> to local s,t coordinates. This method is used mainly for tests.
!> inputs:
!>   n_v:            (integer) number of physical quantities
!>   P_R,P_Z:        (real8)(n_v) physical quantities first derivatives in R and Z
!>   P_RR,P_RZ,P_ZZ: (real8)(n_v) physical quantities second and cross derivatives in R and Z
!>   R_s,R_t:        (real8) major radius first derivatives in s and t
!>   Z_s,Z_t:        (real8) vertical position first derivatives in s and t
!>   R_ss,R_st,R_tt: (real8) major radius second and cross derivatives in s and t
!>   Z_ss,Z_st,Z_tt: (real8) vertical position second and cross derivatives in s and t
!> outputs:
!>   P_ss,P_st,P_tt: (real)(n_v) physical quantities second and cross derivatives in s and t
pure subroutine transform_second_derivatives_RZ_to_st(P_ss,P_st,P_tt,n_v,P_R,P_Z, &
  P_RR,P_RZ,P_ZZ,R_s,R_t,R_ss,R_st,R_tt,Z_s,Z_t,Z_ss,Z_st,Z_tt)
  implicit none
  !> input variables
  integer, intent(in)                      :: n_v
  real(kind=8), dimension(n_v), intent(in) :: P_RR, P_RZ, P_ZZ, P_R, P_Z
  real(kind=8), intent(in)                 :: R_s, R_t, R_ss, R_st, R_tt
  real(kind=8), intent(in)                 :: Z_s, Z_t, Z_ss, Z_st, Z_tt
  !> output varibales
  real(kind=8), dimension(n_v), intent(out) :: P_ss, P_st, P_tt

  P_ss = P_RR*R_s*R_s + 2.d0*P_RZ*R_s*Z_s + P_ZZ*Z_s*Z_s + P_R*R_ss + P_Z*Z_ss 
  P_st = P_RR*R_t*R_s + P_RZ*(Z_t*R_s+R_t*Z_s) + P_ZZ*Z_t*Z_s + P_R*R_st + P_Z*Z_st 
  P_tt = P_RR*R_t*R_t + 2.d0*P_RZ*R_t*Z_t + P_ZZ*Z_t*Z_t + P_R*R_tt + P_Z*Z_tt 
  
end subroutine transform_second_derivatives_RZ_to_st

!---------------------------------------------------------------------------------------

end module mod_coordinate_transforms
