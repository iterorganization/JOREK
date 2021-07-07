module mod_find_rz_nearby
  implicit none
  private
  public :: find_rz_nearby
contains
!> Optimized subroutine to find st coordinates corresponding to x_new=[R_new, Z_new] using
!! The previous values x_old=[R_old, Z_old], st_old(2), i_elm_old and checking adjacent elements first
!! It first checks the current element, using newton iteration to find st_try corresponding to x_new
!! Once this is found it calls check_element_boundary to see if we have crossed an element boundary
!! If this is the case, the search is restarted within the new element.
!! This procedure continues until the required tolerance is reached, or element_try_max is reached.
!!
!! Example:
!! (o is the starting point, T is the target, x is the guess)
!! |^^^^^^|^^^^^^|  |^^^^^^|^^^^^^|  |^^^^^^|^^^^^^|
!! |1     |2     |  |1     |2     |  |1     |2     |
!! |      |  T   |  |      |  T   |  |      | xT   |
!! |  o   |      |  |  o   |x     |  |  o   |      |
!! |______|______|  |______|______|  |______|______|
!! We start out knowing the old position (R,Z) x_old, the new position (R,Z) x_new,
!! the old position and element number st_old, i_elm_old.
!!
!! We are now in element 1, and element_try_index=1. Using 2D newton iteration
!! with backtracking we find the (s,t) position st_try corresponding to x_new, starting
!! from x_old, in the basis functions of the current element. (Which are not
!! guaranteed to be nice outside of the element)
!! 
!! When we have found the position st_try, we check whether it is outside of the
!! current element. (middle figure) If so, we switch to this element (element 2), perform a linear transform
!! of the coordinates and start the process again.
!! If it is inside of the element and we have found a good st_try the routine is done.
!! 
!! There are some cases where this method does not work, notably near element
!! boundaries and the magnetic axis. Here, overshoot of the newton's method near
!! the element boundary causes errors. If too many elements are crossed, i.e.
!! element_try_max is exceeded, or if too many iterations are needed, i.e.
!! newton_iter_max is exceeded, we stop the routine and call find_RZ, which is
!! extremely slow but works for all positions in the domain.
!! In that case, ifail=2:5 is returned.
!! If ifail=-1 the particle is lost
subroutine find_RZ_nearby(node_list, element_list, R_old, Z_old, s_old, t_old, i_elm_old, &
        R_new, Z_new, s_new, t_new, i_elm_new, ifail)
use data_structure
use mod_neighbours
implicit none
!> Input parameters
type (type_node_list),    intent(in)    :: node_list
type (type_element_list), intent(in)    :: element_list
real*8,                   intent(in)    :: R_old, Z_old !< The old R,Z location
real*8,                   intent(in)    :: R_new, Z_new !< The new R,Z location
real*8,                   intent(in)    :: s_old, t_old !< The old st location (used to compute a guess)
integer,                  intent(in)    :: i_elm_old
real*8,                   intent(out)   :: s_new, t_new !< The found new coordinates
integer,                  intent(out)   :: i_elm_new
integer,                  intent(out)   :: ifail !< if ifail = -1 the position could not be found in the grid.
!< ifail > 0 indicates various other cases

!> Accuracy defaults (tolerances are squared!, units of element size)
real*8,  parameter :: element_tolerance   = 1.d-12 !< Tolerance for finding a position inside an element
integer, parameter :: newton_iter_max     = 8 !< Number of iterations to try

!> Internal variables
integer :: newton_iter_number, i_elm_tmp
real*8 :: inv_st_jac_det, R_s, R_t, Z_s, Z_t
real*8 :: st_step(2), x_step(2), x_tmp(2), st_new(2), x_new(2) ! x_step = (R,Z) of trial position
real*8 :: err2, err2_old, dist(2), fact

! Check if element is valid
if (i_elm_old .lt. 1 .or. i_elm_old .gt. element_list%n_elements) then
  call find_RZ(node_list,element_list,R_new,Z_new,x_step(1),x_step(2),i_elm_new,s_new,t_new,ifail)
  return
end if

! Setup initial values
x_step = [R_old,Z_old] ! start at the current position
i_elm_new = i_elm_old ! start in the current element
st_new = [s_old,t_old] ! start at the old position
x_new = [R_new,Z_new]
! Find the jacobian at the current s and t position
call try_interp(node_list,element_list,i_elm_new,st_new,x_step,R_s,R_t,Z_s,Z_t,inv_st_jac_det)
err2 = dot_product(x_step-x_new,x_step-x_new)
ifail=0


! Newton iteration to find s and t in or out of this element
do newton_iter_number = 1, newton_iter_max
  ! Perform newton iteration by calculating the inverse of the jacobian matrix explicitly
  err2_old = err2

  ! Calculate the trial newton step
  st_step(1) = ( Z_t * (x_new(1)-x_step(1)) - R_t * (x_new(2)-x_step(2))) * inv_st_jac_det
  st_step(2) = (-Z_s * (x_new(1)-x_step(1)) + R_s * (x_new(2)-x_step(2))) * inv_st_jac_det

  ! Limit this step if it goes outside of the element
  dist = merge(1-st_new,st_new,st_step .gt. 0) ! dist = 1-s if step > 0, s if step < 0 (distance to 0 or 1)
  fact = maxval(abs(st_step)/max(dist,1d-30)) ! if fact>=1 we are on the boundary
  ! (it is the overshoot: i.e. how many times we overshoot the boundary with one st_step)

  if (fact .ge. 1.d0-1d-12) then
    st_new = st_new + st_step/fact
#ifdef DEBUG
    call try_interp(node_list,element_list,i_elm_new,st_new,x_tmp,R_s,R_t,Z_s,Z_t,inv_st_jac_det)
#endif
    i_elm_tmp = i_elm_new
    call coord_in_neighbour(node_list,element_list,i_elm_tmp,i_elm_new,st_new)
    if (i_elm_new .lt. 0) then
      call find_RZ(node_list,element_list,x_new(1),x_new(2),x_step(1),x_step(2),i_elm_new,s_new,t_new,ifail)
      if (ifail .ne. 0) i_elm_new = 0
    end if
    if (i_elm_new .eq. 0) then ! No element on that side, particle is lost
      i_elm_new = - i_elm_tmp ! Save position of particle
      ! Calculate new R and Z in x_new
      call try_interp(node_list,element_list,i_elm_tmp,st_new,x_new,R_s,R_t,Z_s,Z_t,inv_st_jac_det)
      ! Set new element-local coordinates for the point on the axis
      s_new = st_new(1)
      t_new = st_new(2)
      ifail = -1
      return
    end if

    call try_interp(node_list,element_list,i_elm_new,st_new,x_step,R_s,R_t,Z_s,Z_t,inv_st_jac_det)
#ifdef DEBUG 
    if (norm2(x_step-x_tmp) .gt. 1d-8) then
      !write(*,*) "ERROR on element edge crossing", x_step, x_tmp, norm2(x_step-x_tmp), &
      !i_elm_new, i_elm_old, i_elm_tmp, "deleting particle"
      !call exit(1)
      i_elm_new = 0
      return
    end if
#endif
  else
    st_new = st_new + st_step
    call try_interp(node_list,element_list,i_elm_new,st_new,x_step,R_s,R_t,Z_s,Z_t,inv_st_jac_det)
  end if
  err2 = dot_product(x_step-x_new,x_step-x_new)
  s_new = st_new(1)
  t_new = st_new(2)

  if (err2 < element_tolerance) exit
enddo


if (isnan(err2)) then
  !write(*,*) "WARNING: NaN encountered after newton iteration, using find_RZ"
  call find_RZ(node_list,element_list,x_new(1),x_new(2),x_step(1),x_step(2),i_elm_new,s_new,t_new,ifail)
  if (ifail .eq. 0) ifail=2
  return
endif
if (newton_iter_number .gt. newton_iter_max) then
  !write(*,"(A,i4,A,i5,A,2g14.6,A,3g14.6)") "WARNING: iteration for st did not converge after", newton_iter_max, " tries in element ", i_elm_new, &
  !" using find_RZ", x_new, "err2(old)/convergence: ", err2, err2_old, err2_old/err2
    !write(*,"(A,2g16.8)") "Find_RZ at ", x_new
  call find_RZ(node_list,element_list,x_new(1),x_new(2),x_step(1),x_step(2),i_elm_new,s_new,t_new,ifail)
  if (ifail .eq. 0) ifail=3
  return
endif
end subroutine find_RZ_nearby


!> Auxiliary subroutine for find_RZ_nearby
pure subroutine try_interp(node_list,element_list,i_elm,st,x,R_s,R_t,Z_s,Z_t,inv_st_jac_det)
use data_structure
use mod_interp
implicit none
!> Input parameters
type (type_node_list),    intent(in)    :: node_list
type (type_element_list), intent(in)    :: element_list
real*8,                   intent(in)    :: st(2)
integer,                  intent(in)    :: i_elm
real*8,                   intent(out)   :: x(2), R_s, R_t, Z_s, Z_t, inv_st_jac_det
real*8 :: jac

call interp_RZ(node_list,element_list,i_elm,st(1),st(2),x(1),R_s,R_t,x(2),Z_s,Z_t)
! Guard against the determinant being close to zero
jac = R_s * Z_t - R_t * Z_s
if (abs(jac) .lt. 1d-8) then
  inv_st_jac_det = sign(1d8, jac) 
else
  inv_st_jac_det = 1.d0/(jac)
end if
end subroutine try_interp
end module mod_find_rz_nearby
