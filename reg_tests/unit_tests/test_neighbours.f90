!> Module testing accuracy of calculating sines by applying De Moivre's formula
module test_neighbours
use fruit
use mod_neighbours
use data_structure
use mod_element_rtree, only: rtree_initialized, populate_element_rtree, nearby_elements
implicit none
contains

!> Test the detection of nearby elements with an RTree
subroutine test_nearby_elements
  type(type_node_list), allocatable :: node_list
  type(type_element_list), allocatable :: element_list
  integer, parameter :: nn = 8
  integer, dimension(:), allocatable :: i_nearby
  allocate(node_list, element_list)

  ! Construct a few specific set of elements
  !
  ! 5-------6-------7-------8 x(2) = 1
  ! |       |       |       |
  ! |   1   |   2   |   3   |
  ! |       |       |       |
  ! 1-------2-------3-------4 x(2) = 0
  ! 0       1       2       3
  !x(1)
  call grid_bezier_square(4,2,0.d0,3d0,0.d0,1.d0,.false.,node_list,element_list)

  call populate_element_rtree(node_list, element_list)

  ! Returns in insertion order I think, so everything is fine.
  call nearby_elements(node_list, element_list, 1, i_nearby)
  call assert_equals(2, size(i_nearby), "1 returns 2 elems")
  call assert_equals(1, i_nearby(1), "1 will contain self")
  call assert_equals(2, i_nearby(2), "1 will contain neighbour")

  call nearby_elements(node_list, element_list, 2, i_nearby)
  call assert_equals(3, size(i_nearby), "2 returns 3 elems")
  call assert_equals(1, i_nearby(1), "2 will contain neighbour left") 
  call assert_equals(2, i_nearby(2), "2 will contain self")
  call assert_equals(3, i_nearby(3), "2 will contain neighbour right")

  call nearby_elements(node_list, element_list, 3, i_nearby)
  call assert_equals(2, size(i_nearby), "3 returns 2 elems")
  call assert_equals(2, i_nearby(1), "3 will contain neighbour")
  call assert_equals(3, i_nearby(2), "3 will contain self")
end subroutine test_nearby_elements



!> Test a simple case for the update_neighbours routine
!> TODO:
!> * Test degenerate elements (like axis)
!> * Test multi-node points (like axis, xpoint)
!> * Test elements with different orientation (co/counter)
!>
!> When manually creating grids make sure to set sizes correctly, they are used
!> to determine the element bounding boxes and this in nearby_elements.
subroutine test_update_neighbours
  type(type_node_list), allocatable :: node_list
  type(type_element_list), allocatable :: element_list
  integer, parameter :: nn = 6
  integer :: i, nb(4), inb_i, inb_j
  allocate(node_list, element_list)

  ! Construct a few specific set of elements
  !
  ! 4-------5-------6 x(2) = 1
  ! |       |       |
  ! |   1   |   2   |
  ! |       |       |
  ! 1-------2-------3 x(2) = 0
  ! 0       1       2
  !x(1)
  call grid_bezier_square(3,2,0.d0,2d0,0.d0,1.d0,.false.,node_list,element_list)

  call assert_true(neighbours(node_list, element_list%element(1), element_list%element(2),inb_i,inb_j), "1 must nb 2")
  call assert_equals(2, inb_i, "i on side 2")
  call assert_equals(4, inb_j, "j on side 4")
  call assert_true(neighbours(node_list, element_list%element(2), element_list%element(1),inb_i,inb_j), "2 must nb 1")
  call assert_equals(4, inb_i, "i on side 4")
  call assert_equals(2, inb_j, "j on side 2")

  rtree_initialized = .false. ! Force remake rtree
  call update_neighbours(node_list,element_list)
  nb = [0,2,0,0]
  do i=1,4
    call assert_equals(nb(i), element_list%element(1)%neighbours(i), "element 1 is next to 2")
  end do
  nb = [0,0,0,1]
  do i=1,4
    call assert_equals(nb(i), element_list%element(2)%neighbours(i), "element 2 is next to 1")
  end do
end subroutine test_update_neighbours



subroutine test_coord_in_neighbour_square
  type(type_node_list), allocatable :: node_list
  type(type_element_list), allocatable :: element_list
  real*8 :: st(2)
  integer, parameter :: nn = 6
  real*8, parameter :: tol = 1d-12
  integer :: ito
  allocate(node_list, element_list)

  ! Construct a specific set of elements
  ! 4-------5-------6 x(2) = 1
  ! |       |       |
  ! |   1   |   2   |
  ! |       |       |
  ! 1-------2-------3 x(2) = 0
  ! 0       1       2
  !x(1)
  call grid_bezier_square(3,2,0.d0,2d0,0.d0,1.d0,.false.,node_list,element_list)
  ! All sizes are equal since the elements are square => we can shuffle vertices
  element_list%element(1)%vertex = [1,2,5,4]

  ! 4-------5-------6 (node number)
  ! |4     3|4     3| (vertex number)
  ! |   1   |   2   |
  ! |1     2x1     2|
  ! 1-------2-------3
  element_list%element(2)%vertex = [2,3,6,5]
  rtree_initialized = .false. ! Force remake rtree
  call update_neighbours(node_list,element_list)
  ! Now test 1 position on the boundary and see if the transformation is correct
  st = [1.d0, 0.3d0]
  call coord_in_neighbour(node_list,element_list,1,ito,st)
  call assert_equals(2, ito, "right el 1")
  call assert_equals(0.d0, st(1), tol, "s zero co 1")
  call assert_equals(0.3d0, st(2), tol, "t set co 1")
  call coord_in_neighbour(node_list,element_list,2,ito,st)
  call assert_equals(1.d0, st(1), tol, "s orig co 1")
  call assert_equals(0.3d0, st(2), tol, "t orig co 1")

  ! Rotate the right element and retest
  ! 4-------5-------6 (node number)
  ! |4     3|3     2| (vertex number)
  ! |   1   |   2   |
  ! |1     2x4     1|
  ! 1-------2-------3
  element_list%element(2)%vertex = [3,6,5,2]
  rtree_initialized = .false. ! Force remake rtree
  call update_neighbours(node_list,element_list)
  st = [1.d0, 0.3d0]
  call coord_in_neighbour(node_list,element_list,1,ito,st)
  call assert_equals(2, ito, "right el 2")
  call assert_equals(0.3d0, st(1), tol, "s set co 2")
  call assert_equals(1.d0, st(2), tol, "t one co 2")
  call coord_in_neighbour(node_list,element_list,2,ito,st)
  call assert_equals(1.d0, st(1), tol, "s orig co 2")
  call assert_equals(0.3d0, st(2), tol, "t orig co 2")


  ! Rotate the left element and retest
  ! 4-------5-------6 (node number)
  ! |3     2|3     2| (vertex number)
  ! |   1   |   2   |
  ! |4     1x4     1|
  ! 1-------2-------3
  element_list%element(1)%vertex = [2,5,4,1]
  rtree_initialized = .false. ! Force remake rtree
  call update_neighbours(node_list,element_list)
  st = [0.3d0, 0.d0]
  call coord_in_neighbour(node_list,element_list,1,ito,st)
  call assert_equals(2, ito, "right el 3")
  call assert_equals(0.3d0, st(1), tol, "s same co 3")
  call assert_equals(1.d0, st(2), tol, "t one co 3")
  call coord_in_neighbour(node_list,element_list,2,ito,st)
  call assert_equals(0.3d0, st(1), tol, "s orig co 3")
  call assert_equals(0.0d0, st(2), tol, "t orig co 3")

  ! Reverse orientation of right element and retest
  ! 4-------5-------6 (node number)
  ! |3     2|3     4| (vertex number)
  ! |   1   |   2   |
  ! |4     1x2     1|
  ! 1-------2-------3
  element_list%element(2)%vertex = [3,2,5,6]
  rtree_initialized = .false. ! Force remake rtree
  call update_neighbours(node_list,element_list)
  st = [0.3d0, 0.d0]
  call coord_in_neighbour(node_list,element_list,1,ito,st)
  call assert_equals(2, ito, "right el 4")
  call assert_equals(1d0, st(1), tol, "s 1 counter 4")
  call assert_equals(0.3d0, st(2), tol, "t same counter 4")
  call coord_in_neighbour(node_list,element_list,2,ito,st)
  call assert_equals(0.3d0, st(1), tol, "s orig counter 4")
  call assert_equals(0.0d0, st(2), tol, "t orig counter 4")
end subroutine test_coord_in_neighbour_square
end module test_neighbours
