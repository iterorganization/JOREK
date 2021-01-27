module mod_neighbours
  use data_structure
  implicit none
  private
  public :: neighbours, update_neighbours
  public :: coord_in_neighbour
contains
logical function neighbours(node_list,elm1,elm2,inb1,inb2)
!-------------------------------------------------------
! function to check if the two elements elm1 and elm2
! are neighbours. Two elements are neighbours if they
! share a side, i.e. if two nodes are the same.
! inb1 -> the index of the shared neighbour of elm1
! inb2 -> the index of the shared neighbour of elm2
!   note : does not work for unequal sized neighbours
!-------------------------------------------------------
implicit none

type (type_node_list), intent(in) :: node_list
type (type_element), intent(in)   :: elm1, elm2
integer, intent(out)              :: inb1, inb2

integer :: n1(2), n2(2) ! nodes of side i of elm1 or of side j of elm2
integer :: i, j
neighbours = .false.

! First test by node number
do i=1,4 ! Loop over sides of element 1
  n1 = elm1%vertex(mod([i-1,i],4)+1)
  do j=1,4 ! Loop over sides of element 2
    ! node numbers are related to sides as node1=mod(side-1,4)+1, node2=mod(side,4)+1
    n2 = elm2%vertex(mod([j-1,j],4)+1)
    ! If match cross or straight (i.e. 1->2/2->1 or 1->1/2->2)
    ! Four different cases
    if      (n1(1) .eq. n2(2)) then
      neighbours = node_same_pos(node_list, n1(2), n2(1))
    else if (n1(2) .eq. n2(1)) then
      neighbours = node_same_pos(node_list, n1(1), n2(2))
    else if (n1(1) .eq. n2(1)) then
      neighbours = node_same_pos(node_list, n1(2), n2(2))
    else if (n1(2) .eq. n2(2)) then
      neighbours = node_same_pos(node_list, n1(1), n2(1))
    end if
    if (neighbours) then
      inb1 = i
      inb2 = j
      return
    end if
  enddo
enddo
end function neighbours

pure function node_same_pos(node_list, i1, i2) result(same)
  type (type_node_list), intent(in) :: node_list
  integer, intent(in) :: i1, i2
  logical :: same
  real*8, parameter :: tol=1.d-8
  if (i1 .eq. i2) then
    same = .true.
  else
    if (norm2(node_list%node(i1)%x(1,1,1:2)-node_list%node(i2)%x(1,1,1:2)) .lt. tol) then
      same = .true.
    else
      same = .false.
    end if
  end if
end function node_same_pos


!> Find if two elements (elm1, elm2 which is on side1 of elm1) have the same
!> orientation. Element node numbering must always be consecutive when going
!> co or counter-clockwise around the element.
subroutine neighbours_side_co_counter(element_list,elm1,elm2,side1,side2,neighbours,co)
  type (type_element_list), intent(in) :: element_list
  integer, intent(in)               :: elm1, elm2
  integer, intent(in)               :: side1
  integer, intent(out)              :: side2
  logical, intent(out)              :: neighbours
  logical, intent(out)              :: co !< True if elements are both clockwise or
!< counter-clockwise (i.e. they have the same direction)
  integer :: i, n1(2), n2(2)

  co = .false. ! does not mean anything if not neighbours
  neighbours = .false.

  ! Find the side in elm2 pointing to elm1
  ! If elm2 has no neighbour == elm1 use the last one that is 0 (i.e. the one on
  ! the axis itself.)
  side2 = 0
  do i=1,n_vertex_max
    if (element_list%element(elm2)%neighbours(i) .eq. elm1) then
      side2 = i
      exit ! loop
    end if
  end do
  if (side2 .gt. 0) then
    neighbours = .true.

    ! Determine node numbers of the sides
    ! node numbers are related to sides as node1=mod(side-1,4)+1, node2=mod(side,4)+1
    n1 = element_list%element(elm1)%vertex(mod([side1-1,side1],4)+1)
    n2 = element_list%element(elm2)%vertex(mod([side2-1,side2],4)+1)
    ! Find if match cross or straight (i.e. 1->2/2->1 or 1->1/2->2)
    ! We do not need to check the node position since there is only one option
    ! and one of the node numbers must match.
    if      (n1(1) .eq. n2(2)) then
      co = .true.
    else if (n1(2) .eq. n2(1)) then
      co = .true.
    else if (n1(1) .eq. n2(1)) then
      co = .false.
    else if (n1(2) .eq. n2(2)) then
      co = .false.
    end if
  end if
end subroutine neighbours_side_co_counter


!> Convert s,t coordinates on the boundary of element i_from into s,t coordinates
!> on the boundary of i_to.
!> Element side and node numbering:
!> 4______3 1
!> |\  3 /|
!> | \  / |
!> |  \/  | t
!> |4 /\ 2|
!> | /  \ |
!> |/  1 \| 0
!> 1^^^^^^2 
!>
!> 0   s  1
!>
!> On the boundary between elements the following is guaranteed:
!> * One of the local coordinates (s,t) is either 0 or 1 (1: t=0, 2: s=1, 3: t=1, 4: s=0)
!> * The other coordinates x_i, x_j (elements i and j) are related: |dx_i/dx_j| = 1
subroutine coord_in_neighbour(node_list,element_list,i_from,i_to,st)
  type (type_node_list), intent(in)    :: node_list
  type (type_element_list), intent(in) :: element_list
  integer, intent(in)                  :: i_from
  integer, intent(out)                 :: i_to !<  >0 if a neighbour exists on that side, -1 if search
  real*8, intent(inout)                :: st(2)

  integer :: q_from, q_to !< Quadrants
  logical :: nb, co !< .true. if they have the same orientation (co or counter clockwise)
  real*8  :: x

  if (st(1) .gt. st(2)) then ! 1 or 2
    if (1.d0 - st(1) .gt. st(2)) then
      q_from = 1
    else
      q_from = 2
    end if
  else
    if (1.d0 - st(1) .le. st(2)) then
      q_from = 3
    else
      q_from = 4
    end if
  end if

  i_to = element_list%element(i_from)%neighbours(q_from)
  if (i_to .le. 0) return

  ! Check once more that they are neighbours and determine the orientation
  call neighbours_side_co_counter(element_list,i_from,i_to,q_from,q_to,nb,co)
  if (.not. nb .or. q_to .eq. 0) then
    write(*,"(A,i5,A,i5)") "ERROR IN element_list%element(", i_from, ")%neighbours to ", i_to
    i_to = 0
    return
  end if

  ! x is the coordinate along the boundary from vertex i_side to i_side+1
  select case (q_from)
  case (1)
    x = st(1)
  case (2)
    x = st(2)
  case (3)
    x = 1.d0 - st(1)
  case (4)
    x = 1.d0 - st(2)
  end select

  if (co) x = 1.d0-x ! if the vectors along the boundary are antiparallel

  select case (q_to)
  case (1)
    st = [x,0.d0]
  case (2)
    st = [1.d0,x]
  case (3)
    st = [1.d0-x,1.d0]
  case (4)
    st = [0.d0,1.d0-x]
  end select
end subroutine coord_in_neighbour



subroutine update_neighbours(node_list, element_list, force_rtree_initialize)
use mod_element_rtree, only: rtree_initialized, populate_element_rtree, nearby_elements
implicit none

type (type_node_list), intent(in)       :: node_list
type (type_element_list), intent(inout) :: element_list
logical, intent(in), optional           :: force_rtree_initialize !< default false

type (type_element)      :: elm_i, elm_j
integer                  :: inb_i, inb_j, i, j, k
integer                  :: i_elm, j_elm, i_node1, i_node2
real*8                   :: s_i, t_i, R_i, Rs_i, Rt_i, Rst_i, Rss_i, Rtt_i, Z_i, Zs_i, Zt_i, Zst_i,Zss_i,Ztt_i
real*8                   :: s_j, t_j, R_j, Rs_j, Rt_j, Rst_j, Rss_j, Rtt_j, Z_j, Zs_j, Zt_j, Zst_j,Zss_j,Ztt_j
integer, dimension(:), allocatable :: i_nearby

! Be careful here. If the grid changes the information will be incorrect and you
! need to manually call populate_element_rtree
if (present(force_rtree_initialize)) then
  if (force_rtree_initialize) call populate_element_rtree(node_list, element_list)
else if (.not. rtree_initialized) then
  call populate_element_rtree(node_list, element_list)
end if

!$omp parallel default(none) private(i,k,j,i_node1,i_node2,inb_i,inb_j,i_nearby) shared(element_list,node_list)

!$omp do
do i=1, element_list%n_elements
  element_list%element(i)%neighbours(:) = 0
  call nearby_elements(node_list, element_list, i, i_nearby)
  do k=1,size(i_nearby,1)
    j = i_nearby(k)
    if (i .eq. j) cycle

    if  (neighbours(node_list, element_list%element(i), element_list%element(j),inb_i,inb_j)) then
      element_list%element(i)%neighbours(inb_i) = j
    endif
  enddo
enddo
!$omp end do

! Force search for element crossings over the axis
!$omp do
do i=1, element_list%n_elements
  do j=1, n_vertex_max
    i_node1 = element_list%element(i)%vertex(j)
    i_node2 = element_list%element(i)%vertex(mod(j,4)+1)

    ! If the two nodes for this side are in the same position (degenerate element)
    if (norm2(node_list%node(i_node1)%x(1,1,1:2) - node_list%node(i_node2)%x(1,1,1:2)) .lt. 1d-8) then
      ! Force search for the opposite side
      element_list%element(i)%neighbours(j) = -1
      do k=1,n_vertex_max
        ! If any of the other sides was not found force search there as well
        if(element_list%element(i)%neighbours(k).eq.0) then
           element_list%element(i)%neighbours(k) = -1
        endif
      enddo
    endif
  enddo
enddo
!$omp end do
!$omp end parallel
end subroutine update_neighbours
end module mod_neighbours
