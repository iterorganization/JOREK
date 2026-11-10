!> This modules contains tools needed by the different
!> particle pushers. Generally, they are mathematical 
!> stand-alone routines which can be used by multiple
!> independent algorithm
module mod_pusher_tools
implicit none
private
!> public procedures
public get_orthonormals
public cayley_transform, approximated_cayley_transform
public gc_position_to_particle, particle_position_to_gc
public cylindrical_displacement, larmor_displacement
public GC_GEOMETRY_INVALID
integer, parameter :: GC_GEOMETRY_INVALID = -1001
contains

!---------------------------------------------------------------------------

!> Get two vectors orthogonal to a given vector.
!> JOREK uses clockwise phi: (eR,eZ,ephi) is right-handed.
!> Retain this perpendicular-basis convention for existing gyrophases.
pure subroutine get_orthonormals(b, e1, e2)
  use mod_math_operators, only: cross_product
  implicit none

  real*8, dimension(3), intent(in)  :: b !< Does not need to be normalized
  real*8, dimension(3), intent(out) :: e1, e2

  ! Take r as a reference vector (this will fail therefore if B is purely in r direction)
  if (b(2) .eq. 0.d0 .and. b(3) .eq. 0.d0) then
    e1 = [0.d0, 1.d0, 0.d0]
    e2 = [0.d0, 0.d0, 1.d0]
  else
    e1 = cross_product(b, [1.d0, 0.d0, 0.d0])
    ! Normalize
    e1 = e1/norm2(e1)
    ! Obtain a second reference vector
    e2 = cross_product(b, e1)
    ! Normalize
    e2 = e2/norm2(e2)
  end if
end subroutine get_orthonormals

!> Embed a physical displacement (metres, local eR/eZ/ephi components).
!> JOREK Cartesian coordinates are (R*cos(phi),-R*sin(phi),Z), so phi is
!> clockwise. The parent angle remains continuous, including across 2*pi.
pure subroutine cylindrical_displacement(x, displacement, x_out, delta_phi, ifail)
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  real*8, intent(in) :: x(3), displacement(3)
  real*8, intent(out) :: x_out(3), delta_phi
  integer, intent(out) :: ifail
  real*8 :: radial, radius

  x_out = x
  delta_phi = 0.d0
  ifail = GC_GEOMETRY_INVALID
  if (.not. all(ieee_is_finite(x))) return
  if (.not. all(ieee_is_finite(displacement))) return
  if (x(1) <= 0.d0) return
  radial = x(1) + displacement(1)
  radius = hypot(radial, displacement(3))
  if (.not. ieee_is_finite(radius)) return
  if (radius <= 0.d0) return ! The cylindrical basis is undefined on the axis.
  delta_phi = atan2(displacement(3), radial)
  x_out = [radius, x(2) + displacement(2), x(3) + delta_phi]
  if (.not. all(ieee_is_finite(x_out))) then
    x_out = x
    delta_phi = 0.d0
    return
  endif
  ifail = 0
end subroutine cylindrical_displacement

!> rho = b cross p / (q e B). Momentum is in amu*m/s; B and p are
!> expressed in the SAME local basis. For JOREK's CLOCKWISE phi, (eR,eZ,ephi)
!> is right-handed: the ordinary array cross product has the physical sign.
!> A neutral has no Larmor displacement. Invalid charged-field data fail explicitly.
pure subroutine larmor_displacement(momentum, charge, B_hat, B_norm, rho, ifail)
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use constants, only: ATOMIC_MASS_UNIT, EL_CHG
  use mod_math_operators, only: cross_product
  real*8, intent(in) :: momentum(3), B_hat(3), B_norm
  integer(kind=1), intent(in) :: charge
  real*8, intent(out) :: rho(3)
  integer, intent(out) :: ifail
  rho = 0.d0
  ifail = GC_GEOMETRY_INVALID
  if (.not. all(ieee_is_finite(momentum))) return
  if (charge == 0) then
    ifail = 0
    return
  endif
  if (.not. all(ieee_is_finite(B_hat)) .or. .not. ieee_is_finite(B_norm)) return
  if (B_norm <= 0.d0 .or. norm2(B_hat) <= 0.d0) return
  rho = (ATOMIC_MASS_UNIT/(real(charge,8)*EL_CHG*B_norm))*cross_product(B_hat,momentum)
  if (.not. all(ieee_is_finite(rho))) then
    rho = 0.d0
    return
  endif
  ifail = 0
end subroutine larmor_displacement

!---------------------------------------------------------------------------

!> This function computes the right-handed Cayley transform of a vector vec multiplied 
!> by a scalar alpha. The Cayley transform is defined as:
!> cayley(alpha*B) = (I-alpha*B)^(-1) * (I+alpha*B)
!> where B is the vector product skew symmetric matrix of the vector vec
!> and I is the identity matrix.
pure function cayley_transform(alpha,vec)
  ! defining input variables
  real(kind=8),intent(in) :: alpha !< multiplicative constant
  real(kind=8),dimension(3),intent(in) :: vec !< vector to be transformed
  ! defining output variables
  real(kind=8),dimension(3,3) :: cayley_transform !< Cayley transform of vec
  ! defining internal variables
  real(kind=8),dimension(3,3) :: A,B !< (I-alpha*B)^(-1) and (I+alpha*B)

! computing (I+alpha*B)
  B(1:3,1) = (/1.d0,-alpha*vec(3),alpha*vec(2)/)
  B(1:3,2) = (/alpha*vec(3),1.d0,-alpha*vec(1)/)
  B(1:3,3) = (/-alpha*vec(2),alpha*vec(1),1.d0/)

! computing (I-alpha*B)^(-1)
  A(1:3,1) = (/(1.0 + alpha*alpha*vec(1)*vec(1)),(alpha*(alpha*vec(1)*vec(2) - vec(3))),&
            (alpha*(alpha*vec(3)*vec(1) + vec(2)))/)
  A(1:3,2) = (/(alpha*(alpha*vec(2)*vec(1) + vec(3))),(1.0 + alpha*alpha*vec(2)*vec(2)),&
            (alpha*(alpha*vec(3)*vec(2) - vec(1)))/)
  A(1:3,3) = (/(alpha*(alpha*vec(3)*vec(1) - vec(2))),(alpha*(alpha*vec(2)*vec(3) + vec(1))),&
            (1.0 + alpha*alpha*vec(3)*vec(3))/)

! computing the Cayley transform
  cayley_transform = (matmul(A,B))/(1.d0 + alpha*alpha*(vec(1)*vec(1)+vec(2)*vec(2)+vec(3)*vec(3)))

end function cayley_transform

!---------------------------------------------------------------------------

!> This function computes the approximated Cayley transform a vector vec 
!> multiplied by a constant alpha as reported in R. Zhang et al., Phys. Plasmas 22 (2015) 044501.
pure function approximated_cayley_transform(alpha,vec)
  ! defining input variables
  real(kind=8),intent(in) :: alpha !< multiplicative constant
  real(kind=8),dimension(3),intent(in) :: vec !< vector to be transformed
  ! defining output variables
  real(kind=8),dimension(3,3) :: approximated_cayley_transform !< approximated Cayley transform of V
  ! internal variables
  real(kind=8) :: coefficient

  ! compute the approximated transfrom coefficient
  coefficient = (2.0*alpha)/(1.0 + alpha*alpha*(dot_product(vec,vec)))

  ! compute the approximated Cayley transform
  approximated_cayley_transform(1:3,1) = &
    (/1.0-coefficient*alpha*(vec(3)*vec(3) + vec(2)*vec(2)),&
    coefficient*(vec(3) + alpha*vec(2)*vec(1)),&
    coefficient*(alpha*vec(3)*vec(1) - vec(2))/)
  approximated_cayley_transform(1:3,2) = &
    (/coefficient*(alpha*vec(1)*vec(2) - vec(3)),&
    1.0-coefficient*alpha*(vec(3)*vec(3) + vec(1)*vec(1)),&
    coefficient*(alpha*vec(3)*vec(2) + vec(1))/)
  approximated_cayley_transform(1:3,3) = &
    (/coefficient*(alpha*vec(1)*vec(3) + vec(2)),&
    coefficient*(alpha*vec(2)*vec(3) - vec(1)),1.0-&
    coefficient*alpha*(vec(2)*vec(2) + vec(1)*vec(1))/)
end function approximated_cayley_transform

!---------------------------------------------------------------------------

!> This subroutine computes the guiding center position
!> from the particle position and momentum
subroutine particle_position_to_gc(node_list,element_list,&
  x_in,st_in,i_elm_in,p_in,q_in,B_hat,B_norm,x_gc_out,st_gc_out,i_elm_out)
  use data_structure
  use constants, only: ATOMIC_MASS_UNIT,EL_CHG
  use mod_find_rz_nearby
  ! declare input variables
  type(type_node_list), intent(in)       :: node_list
  type(type_element_list), intent(in)    :: element_list
  integer(kind=1), intent(in)            :: q_in !< particle charge
  integer(kind=4), intent(in)            :: i_elm_in !< particle element
  real(kind=8), dimension(3), intent(in) :: x_in, p_in !< particle position and momentum
  real(kind=8), dimension(2), intent(in) :: st_in !< particle local coordinates
  real(kind=8), dimension(3), intent(in) :: B_hat !< Magnetic field direction B/B_norm
  real(kind=8), intent(in)               :: B_norm !< magnetic field intensity in [T]
  ! declare output variables
  integer(kind=4), intent(out)            :: i_elm_out !< gc element
  real(kind=8), dimension(2), intent(out) :: st_gc_out !< local gc position s,t
  real(kind=8), dimension(3), intent(out) :: x_gc_out  !< global position in R,Z,phi
  ! declare internal variables
  integer :: ifail
  real*8 :: rho(3), delta_phi

  x_gc_out = x_in
  st_gc_out = st_in
  call larmor_displacement(p_in,q_in,B_hat,B_norm,rho,ifail)
  if (ifail == 0) call cylindrical_displacement(x_in,-rho,x_gc_out,delta_phi,ifail)
  if (ifail /= 0) then
    i_elm_out = 0
    return
  endif

  ! find the local coordinates
  call find_RZ_nearby(node_list,element_list,x_in(1),x_in(2), &
    st_in(1),st_in(2),i_elm_in,x_gc_out(1),x_gc_out(2),       &
    st_gc_out(1),st_gc_out(2),i_elm_out,ifail)
end subroutine particle_position_to_gc

!---------------------------------------------------------------------------
!> This subroutine computes the particle position from the guding center position
subroutine gc_position_to_particle(node_list,element_list, &
  x_gc_in,st_gc_in,i_elm_in,p_gc_in,q_gc_in,B_hat,B_norm,x_out,st_out,i_elm_out)
  use data_structure
  use constants, only: ATOMIC_MASS_UNIT,EL_CHG
  use mod_find_rz_nearby
  ! declare input variables
  type(type_node_list), intent(in)       :: node_list
  type(type_element_list), intent(in)    :: element_list
  integer(kind=1), intent(in)            :: q_gc_in !< gc charge
  integer(kind=4), intent(in)            :: i_elm_in !< gc element
  real(kind=8), dimension(3), intent(in) :: x_gc_in,p_gc_in !< gc position and momentum
  real(kind=8), dimension(2), intent(in) :: st_gc_in !< gc local coordinates
  real(kind=8), dimension(3), intent(in) :: B_hat !< Magnetic field direction B/B_norm
  real(kind=8), intent(in)               :: B_norm !< magnetic field intensity in [T]
  ! declare output variables
  integer(kind=4), intent(out)            :: i_elm_out !< particle element
  real(kind=8), dimension(2), intent(out) :: st_out !< local particle position s,t
  real(kind=8), dimension(3), intent(out) :: x_out  !< global oarticle position in R, Z, phi
  ! declare internal variables
  integer :: ifail
  real*8 :: rho(3), delta_phi

  x_out = x_gc_in
  st_out = st_gc_in
  call larmor_displacement(p_gc_in,q_gc_in,B_hat,B_norm,rho,ifail)
  if (ifail == 0) call cylindrical_displacement(x_gc_in,rho,x_out,delta_phi,ifail)
  if (ifail /= 0) then
    i_elm_out = 0
    return
  endif

  ! find the local coordinates
  call find_RZ_nearby(node_list,element_list,x_gc_in(1),x_gc_in(2), &
    st_gc_in(1),st_gc_in(2),i_elm_in,x_out(1),x_out(2),             &
    st_out(1),st_out(2),i_elm_out,ifail)
end subroutine gc_position_to_particle

!---------------------------------------------------------------------------

end module mod_pusher_tools
