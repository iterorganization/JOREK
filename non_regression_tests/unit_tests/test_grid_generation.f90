!> Verify the generation of a few simple grids
module test_grid_generation
use fruit
use data_structure
implicit none
contains

!> Test some reasonable properties of our simple grid.
!> note that we should not use knowledge of the grid construction algorithm
!> and we should not use the node ordering or element ordering in any way.
subroutine test_square_grid
!  type(type_node_list) :: node_list
!  type(type_element_list) :: element_list
  use nodes_elements
  integer, parameter :: n = 10
  real*8, parameter :: R_geo = 1.d0, Z_geo = 0.d0, amin = 0.5d0
  character(len=40) :: s
  integer :: i, j, c
  real*8 :: x(4), z(4), area, xt(2), quad_area

  node_list%n_nodes = 0
  element_list%n_elements = 0
  call grid_bezier_square(n, n, R_geo-amin,R_geo+amin, Z_geo-amin, Z_geo+amin, .true., node_list, element_list)

  ! Verify all nodes are in box R_geo-amin, R_geo+amin, Z_geo-amin,Z_geo+amin
  call assert_true(minval(node_list%node(1:node_list%n_nodes)%x(1,1)) .ge. R_geo-amin, 'All nodes r >= R-a')
  call assert_true(maxval(node_list%node(1:node_list%n_nodes)%x(1,1)) .le. R_geo+amin, 'All nodes r <= R+a')

  call assert_true(maxval(node_list%node(1:node_list%n_nodes)%x(1,2)) .ge. Z_geo-amin, 'All nodes z >= R-a')
  call assert_true(maxval(node_list%node(1:node_list%n_nodes)%x(1,2)) .le. Z_geo+amin, 'All nodes z <= R+a')

  ! Verify there are no duplicate node positions
  do i=1,node_list%n_nodes
    do j=1,node_list%n_nodes
      if (i .ne. j) then
        write(s,'(A,i3,A,i3,A)') '(', i, ',', j, ')'
        call assert_false(norm2(node_list%node(i)%x(1,:) - node_list%node(j)%x(1,:)) .le. 1d-20, 'Node in same position ' // s)
      end if
    end do
  end do

  ! Test if every element has valid nodes
  do i=1,element_list%n_elements
    do j=1,n_vertex_max
      write(s,'(A,i3,A,i3,A)') 'element ', i, ', vertex ', j, ')'
      call assert_true(element_list%element(i)%vertex(j) .le. node_list%n_nodes, 'element has vertex' // s)
      call assert_true(element_list%element(i)%vertex(j) .gt. 0, 'element has vertex 2 ' // s)
    end do
  end do

  ! Test every node is used by 1-4 elements
  do i=1,node_list%n_nodes
    c = 0
    do j=1,n_vertex_max
      c = c + count(element_list%element(1:element_list%n_elements)%vertex(j) .eq. i)
      write(s,'(A,i3,A,i3,A)') 'node ', i, ', vertex ', j, ')'
    end do
    call assert_true(c .gt. 0, 'no orphan nodes ' // s)
    call assert_true(c .le. n_vertex_max, 'no overused nodes ' // s)
  end do

  ! Test there are no nodes inside the quadrilateral defined by the element vertices
  ! https://stackoverflow.com/questions/5922027/how-to-determine-if-a-point-is-within-a-quadrilateral
  do i=1,element_list%n_elements
    x = node_list%node(element_list%element(i)%vertex(:))%x(1,1)
    z = node_list%node(element_list%element(i)%vertex(:))%x(1,2)

    ! Calculate area of quad
    xt = [sum(x)/4, sum(z)/4] ! test point in quad
    quad_area = triangle_area([x(1),z(1)], [x(2),z(2)], xt) + &
                triangle_area([x(2),z(2)], [x(3),z(3)], xt) + &
                triangle_area([x(3),z(3)], [x(4),z(4)], xt) + &
                triangle_area([x(4),z(4)], [x(1),z(1)], xt)

    do j=1,node_list%n_nodes
      if (.not. any(element_list%element(i)%vertex .eq. j)) then
        xt = node_list%node(j)%x(1,:)
        ! Calculate sum of areas of triangles between point and consecutive pairs of vertices
        area = triangle_area([x(1),z(1)], [x(2),z(2)], xt) + &
               triangle_area([x(2),z(2)], [x(3),z(3)], xt) + &
               triangle_area([x(3),z(3)], [x(4),z(4)], xt) + &
               triangle_area([x(4),z(4)], [x(1),z(1)], xt)
        write(s,'(A,i3,A,i3)') 'element', i, ', node ', j
        call assert_true(area .gt. quad_area, 'node outside of quad of other element ' // s)
      end if
    end do
  end do
end subroutine test_square_grid

!> Use Heron's formula to calculate the area
function triangle_area(x1, x2, x3) result(area)
  real*8, intent(in), dimension(2) :: x1, x2, x3 ! corners
  real*8 :: s, a, b, c
  real*8 :: area

  a = norm2(x1-x2)
  b = norm2(x2-x3)
  c = norm2(x3-x1)
  s = (a + b + c)/2
  area = sqrt(max(s*(s-a)*(s-b)*(s-c), 0.d0))
end function triangle_area
end module test_grid_generation
