!> Find the element and position at this psi, with poloidal angle theta.
!> This can be useful for experimenting with low-discrepancy sequences.
!> Root finding on the line from R_axis, Z_axis, parametrised by u.
!> R = R_axis + u cos(theta); Z = Z_axis + u sin(theta).
!> The objective function is (psi - psi(u)) where psi(u) = psi(st(RZ(u)))
!> and we calculate st with find_RZ_nearby and RZ as above.
!> Use Newton's method to find the root. This code could be upgraded to use a
!> higher-order method, as we also have the second derivatives.
!> Be careful when using this routine with psi_n > 1, as multiple intersections
!> can exist.
subroutine find_theta_psi(node_list,element_list,psi_minmax,theta,psi,phi,R_axis,Z_axis,i_elm,s,t,R,Z)
use constants
use data_structure
use mod_interp
use mod_find_rz_nearby
use, intrinsic :: ieee_arithmetic, only: ieee_is_nan

implicit none
! Allow using .cross. to mean a right-handed cross product
interface operator(.cross.)
  pure function right_handed_cross_product(a, b) result(cross)
    real*8, dimension(3), intent(in) :: a, b
    real*8, dimension(3) :: cross
  end function right_handed_cross_product
end interface

! --- Routine parameters
type (type_node_list),    intent(in)    :: node_list
type (type_element_list), intent(in)    :: element_list
real*8, dimension(element_list%n_elements,2), intent(in) :: psi_minmax !< list of minima and maxima of psi in these elements
real*8,                   intent(in)    :: theta
real*8,                   intent(in)    :: psi, phi
real*8,                   intent(in)    :: R_axis
real*8,                   intent(in)    :: Z_axis
integer,                  intent(out)   :: i_elm
real*8,                   intent(out)   :: s, t, R, Z

logical, parameter :: verbose = .true.

! --- Internal variables
real*8  :: u, du, R_try, Z_try, s_out, t_out, err
real*8, dimension(1) :: P, P_s, P_t, P_phi
real*8  :: R_s, R_t, Z_s, Z_t
real*8  :: inv_st_jac, psi_R, psi_Z, psi_u
integer :: i_elm_out, ifail, backtrack_step, newton_iter_number
integer, parameter :: num_backtrack_steps = 10 ! to prevent going over the border of the domain
real*8,  parameter :: backtrack_factor = 0.99d0
integer, parameter :: newton_iter_max = 12
integer, parameter :: i_var(1) = [1], n_ivar=1
real*8,  parameter :: tol = 1d-8

integer :: i, j, j1, j2, i_root
logical :: out_of_domain, correct_quadrant
logical, dimension(element_list%n_elements) :: psi_right
real*8, dimension(n_dim) :: x1, x2
real*8 :: u1, u2, lambda, root, r1, r2, r3
real*8 :: s1, t1, s2, t2, psi1, psi2, psi_u1, psi_u2
integer :: ielm1, ielm2
real*8, dimension(1) :: A, B, C, D
real*8, dimension(2) :: xi
real*8, dimension(3) :: intersection
integer  ::  vertex_start, vertex_step

! 0. Preparation
out_of_domain = .false.

! 1. Find the right elements (all elements where the theta and psi lines enter)
psi_right = (psi_minmax(:,1) .lt. psi) .and. (psi_minmax(:,2) .gt. psi)
u = 0.d0
i_elm_out = 0
do i=1,element_list%n_elements
  if (psi_right(i)) then
    ! Calculate crossings of element vertex and line from R_axis, Z_axis
    u1 = 0.d0
    j1 = 0
    j2 = 0
    if(node_list%node(element_list%element(i)%vertex(1))%axis_node) then
       vertex_start=2
       vertex_step=2
    else
       vertex_start=1
       vertex_step=1
    endif
    do j=vertex_start,n_vertex_max,vertex_step
      x1 = node_list%node(element_list%element(i)%vertex(j))%x(1,1,:)
      x2 = node_list%node(element_list%element(i)%vertex(mod(j,n_vertex_max)+1))%x(1,1,:)
      ! Intersection point in homogeneous coordinates
      ! (https://en.wikipedia.org/wiki/Line%E2%80%93line_intersection#Using_homogeneous_coordinates)
      ! See also http://robotics.stanford.edu/~birch/projective/node4.html
      ! equations of the lines:
      ! vertex-line: a1 x + b1 y + c1 = 0 written as (a1, b1, c1)
      !   where c1 = x1y2-y2x1, b1 = x2-y1, a1 = y1-y2 (this is the cross-product of (x1,y1,1) and (x2,y2,1))
      ! theta-line: (a2, b2, c2)
      !   where c1 = x0 sin(theta) - y0 cos(theta), b1 = cos(theta), a1 = -sin(theta)
      ! The intersection point is then given by (vertex-line) x (theta-line) (cross product, normal orientation)
      ! intersection-point: (ap, bp, cp)
      !   where cp = a1b2-a2b1, bp = a2c1-a1c2, ap = b1c2-b2c1
      ! convert this to a real point by taking x = ap/cp, y = bp/cp (perhaps only if cp is larger than some value)
      ! and check if u and lambda (normalized distance from node i) are within bounds
      intersection = ([x1(1),x1(2),1.d0] .cross. [x2(1),x2(2),1.d0]) &
          .cross.  [-sin(theta),cos(theta),R_axis*sin(theta)-Z_axis*cos(theta)] ! intersection point in homogeneous coordinates
      if (abs(intersection(3)) .lt. 1d-8) cycle ! we will probably not stay in the domain, skip
      xi = intersection(1:2)/intersection(3) ! xy of intercept
      u  = norm2(xi-[R_axis,Z_axis])
      lambda = norm2(xi-x1)/norm2(x2-x1)
      if( sign(1.d0,xi(1)-R_axis)*sign(1.d0,cos(theta)) .gt. 0.d0 .and. &
           sign(1.d0,xi(2)-Z_axis)*sign(1.d0,sin(theta)) .gt. 0.d0 ) then
         correct_quadrant = .true.
      else
         correct_quadrant = .false.
      endif
      if (lambda .ge. 0.d0 .and. lambda .le. 1.d0 .and. u .gt. 0.d0 .and. correct_quadrant) then ! if it is a crossing inside the element
        if (u1 .lt. 1d-30) then ! and we have not yet found the first one
          u1 = u ! store values of first crossing
          j1 = j
        else
          u2 = u
          j2 = j
          exit
        end if
      end if
    end do

    ! if this element has two proper crossings
    if (j1 .gt. 0 .and. j2 .gt. 0) then
      ! Calculate coordinates of first crossing
      ielm1=i ! help the search a bit by starting in this element
      call find_RZ(node_list,element_list,R_axis+u1*cos(theta),Z_axis+u1*sin(theta),R,Z,ielm1,s1,t1,ifail)
      ! Calculate coordinates of second crossing
      ielm2=i
      call find_RZ(node_list,element_list,R_axis+u2*cos(theta),Z_axis+u2*sin(theta),R,Z,ielm2,s2,t2,ifail)
      ! Calculate psi and derivative at first crossing
      if (ielm1 .eq. 0 .or. ielm2 .eq. 0) then
        if (verbose) write(*,*) "Cannot find intersection point in element, ielm1=", ielm1, ", ielm2=", ielm2, ", istart=",i, &
          R_axis+u1*cos(theta), Z_axis+u1*sin(theta), &
          R_axis+u2*cos(theta), Z_axis+u2*sin(theta)
        return
      end if
      call interp_PRZ(node_list,element_list,ielm1,i_var,n_ivar,s1,t1,phi, &
          P,P_s,P_t,P_phi,R,R_s,R_t,Z,Z_s,Z_t)
      inv_st_jac = 1.d0/(R_s * Z_t - R_t * Z_s)
      psi_R    = (  P_s(1) * Z_t - P_t(1) * Z_s ) * inv_st_jac
      psi_Z    = (- P_s(1) * R_t + P_t(1) * R_s ) * inv_st_jac
      psi1     = P(1)
      psi_u1   = (psi_R*cos(theta) + psi_Z*sin(theta))
      ! Calculate psi and derivative at second crossing
      call interp_PRZ(node_list,element_list,ielm2,i_var,n_ivar,s2,t2,phi, &
          P,P_s,P_t,P_phi,R,R_s,R_t,Z,Z_s,Z_t)
      inv_st_jac = 1.d0/(R_s * Z_t - R_t * Z_s)
      psi_R    = (  P_s(1) * Z_t - P_t(1) * Z_s ) * inv_st_jac
      psi_Z    = (- P_s(1) * R_t + P_t(1) * R_s ) * inv_st_jac
      psi2     = P(1)
      psi_u2   = (psi_R*cos(theta) + psi_Z*sin(theta))
      ! Get coefficients of the 3rd-order polynomial
      ! v = (u - u1)/(u2-u1) ! between 0 (at u1) and 1 (at u2)
      ! psi = A + B v + C v^2 + D v^3
      A = psi1
      B = psi_u1
      C = 3.d0*(psi2-psi1) - 2.d0*psi_u1 - psi_u2
      D = 2.d0*(psi1-psi2) + psi_u1 + psi_u2
      ! 0 = (A - psi) + ....
      ! Now solve the equation for psi
      call SOLVP3(A-psi,B,C,D,r1,r2,r3,ifail)
      ! Pick the root closest to 0.5
      i_root = minloc(abs([r1,r2,r3]-0.5d0),1)
      root = r1 ! if (i_root .eq. 1) is implied
      if (i_root .eq. 2) root = r2
      if (i_root .eq. 3) root = r3
      u = root*(u2-u1)+u1
      ! check if it is reasonable close to the element
      if (root .gt. -0.5d0 .and. root .lt. 1.5d0) then
        ! Exit the element check loop
        ! set i_elm_out to speed up find_RZ later
        i_elm_out = i
        exit
      end if
    end if
    ! workaround for the one element at the boundary
  end if
end do
if (u .lt. 0.d0) then
  if (verbose) write(*,"(A,g13.6,A,g13.6,A)") "WARNING: no suitable elements found, skipping for psi=", psi, " theta=", theta/PI, "pi", "u=", u
  i_elm = 0
  return
end if


! prepare looping
! Keep trying with smaller u until we find an element
! this helps a bit for particles very close to the edge
do backtrack_step = 0, num_backtrack_steps
  R_try = R_axis + u*cos(theta)
  Z_try = Z_axis + u*sin(theta)
  call find_RZ(node_list,element_list,R_try,Z_try,R,Z,i_elm_out,s_out,t_out,ifail)
  if (ifail .eq. 0) then
    exit
  else
    u = u*backtrack_factor
  end if
end do
if (backtrack_step .gt. num_backtrack_steps .or. ifail .ne. 0 .or. i_elm_out .eq. 0) then
  if (verbose) write(*,*) "Cannot find initial position after ", backtrack_step
  i_elm = 0
  return
end if


! 2. Iterate to find the right value of psi(u)
do newton_iter_number = 1, newton_iter_max
  ! Calculate psi and the derivative in u direction
  call interp_PRZ(node_list,element_list,i_elm_out,i_var,n_ivar,s_out,t_out,phi, &
      P,P_s,P_t,P_phi,R,R_s,R_t,Z,Z_s,Z_t)
  inv_st_jac = 1.d0/(R_s * Z_t - R_t * Z_s)
  psi_R    = (  P_s(1) * Z_t - P_t(1) * Z_s ) * inv_st_jac
  psi_Z    = (- P_s(1) * R_t + P_t(1) * R_s ) * inv_st_jac
  psi_u    = (psi_R*cos(theta) + psi_Z*sin(theta))

  ! Swap variables
  i_elm = i_elm_out
  s = s_out
  t = t_out
  ! Calculate the error, exit if it is low enough
  err = (P(1) - psi)
  if (abs(err) .lt. tol) return

  ! Formulate a new guess (u) based on psi and psi'
  du = -err/psi_u
  ! Make sure it is not negative
  do while (-du .gt. u)
    du = du*backtrack_factor
  end do

  !write(*,"(12g14.6)") psi, P(1)-psi, u, du, i_elm_out, s, t, out_of_domain
  !if (out_of_domain .and. du .gt. 0.d0) du = -0.5d0*du ! u must go down in this case
  ! Calculate R and Z for this guess
  R_try = R_axis + (u+du)*cos(theta)
  Z_try = Z_axis + (u+du)*sin(theta)
  ! Find st based on the old value and the size of this step
  call find_RZ_nearby(node_list, element_list, R, Z, s, t, i_elm, R_try, Z_try, s_out, t_out, i_elm_out, ifail)
  u = u+du
  if (ifail .lt. 0 .or. i_elm_out .eq. 0) then
    out_of_domain = .true.
    i_elm_out = i_elm ! keep going in the fields of the old element
  else
    out_of_domain = .false.
  end if
  if (ieee_is_nan(P(1)) .or. abs(u) .gt. 1d5 .or. abs(P(1)) .gt. 1d5) then
    if (verbose) write(*,"(A,g13.6,A,g13.6,A,3g13.6,A,g13.6,A)") "Position not found for psi=", psi, ", theta=",theta, " last guess=", u, R_try, Z_try, " (psi=",P(1),")"
    i_elm = 0
    return
  end if
end do

! 3. Warn if failed
if (newton_iter_number .gt. newton_iter_max) then
  ! Indicate this by setting i_elm
  if (verbose) write(*,*) "Too many iterations, skipping with an error of ", err
  i_elm = 0
end if
end subroutine find_theta_psi

!> The cross product in a right-handed coordinate system (e.g. XYZ or RPhiZ)
pure function right_handed_cross_product(a, b)
  real*8, dimension(3) :: right_handed_cross_product
  real*8, dimension(3), intent(in) :: a, b

  right_handed_cross_product(1) = a(2) * b(3) - a(3) * b(2)
  right_handed_cross_product(2) = a(3) * b(1) - a(1) * b(3)
  right_handed_cross_product(3) = a(1) * b(2) - a(2) * b(1)
end function right_handed_cross_product
