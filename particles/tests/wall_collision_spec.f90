!> This module contains tests for the wall collision algorithm.
module wall_collision_spec
use mod_wall_collision
use constants, only: PI, TWOPI
use fruit

implicit none
contains

!> Test the octree construct that contains the wall triangles
subroutine test_octree
  
  type(octree_triangle)      :: triangles(4)
  type(octree_node)          :: tree 
  type(octree_node), pointer :: node
  integer :: max_depth, i, err
  real(kind = 8), dimension(2,3) :: boundary, emptybox
  real(kind = 8) :: x(3)


  boundary = reshape((/ -2, 2, -2, 2, -2, 2 /), shape(boundary))

  ! Large triangle that spawns half of the xy space
  triangles(1)%triangle_id = 1
  triangles(1)%v0 = (/ -1.5, -1.5, 0.5 /)
  triangles(1)%v1 = (/  1.5, -1.5, 0.5 /)
  triangles(1)%v2 = (/ -1.5,  1.5, 0.5 /)

  ! Small triangle completely within a small box
  triangles(2)%triangle_id = 2
  triangles(2)%v0 = (/ -0.8, 0.3, 1.5 /)
  triangles(2)%v1 = (/ -0.2, 0.3, 1.5 /)
  triangles(2)%v2 = (/ -0.5, 0.8, 1.5 /)

  ! Small triangle with one vertex at the edge of a box
  triangles(3)%triangle_id = 3
  triangles(3)%v0 = (/ -0.8, -0.5, 2.0 /)
  triangles(3)%v1 = (/ -0.2, -0.5, 2.0 /)
  triangles(3)%v2 = (/ -0.5,  0.0, 2.0 /)

  ! Small triangle whose side goes through box corner
  triangles(4)%triangle_id = 4
  triangles(4)%v0 = (/ -0.5, 1.5, 2.0 /)
  triangles(4)%v1 = (/  0.5, 0.5, 2.0 /)
  triangles(4)%v2 = (/  0.5, 1.5, 2.0 /)

  ! First initialization
  max_depth = 2
  call octree_init(triangles, max_depth, boundary, tree, err)
  call assert_equals(0, err, "Initialization successful.")
 
  ! Test that all four triangles can be found in this box
  x = (/ -1.5, 1.5, 1.5 /)
  call octree_find(tree, x, node, err)
  call assert_equals(0, err, "Query successful.")
  err = size(node%contained)
  call assert_equals(4, err, "Queried node contains all triangles.")

  ! No triangles in this box
  x = (/ -1.0, -1.0, -1.0 /)
  call octree_find(tree, x, node, err)
  call assert_equals(0, err, "Query successful.")
  err =size(node%contained)
  call assert_equals(0, err, "Queried empty node.")
  emptybox = node%boundary ! Store for later use

  ! Error when trying to query outside octree
  x = (/ 3.0, 0.0, 0.0 /)
  call octree_find(tree, x, node, err)
  call assert_equals(1, err, "Query fail as point outside octree.")

  ! Free resources and second initialization
  call octree_free(tree)

  max_depth = 3
  call octree_init(triangles, max_depth, boundary, tree, err)
  call assert_equals(0, err, "A new initialization successful.")

  ! A box that fully contains a triangle and corners two others
  x = (/ -0.5, 0.5, 1.5 /)
  call octree_find(tree, x, node, err)
  call assert_equals(0, err, "Query successful.")

  err = 3
  do i = 1,size(node%contained)
     if( any( node%contained(i)%triangle_id .eq. (/2,3,4/) ) ) err = err - 1
  end do
  call assert_equals(0, err, "Queried node fully contains one triangle, vertex of second triangle, and corners side of a third.")

  ! A box completely inside a triangle
  x = (/ -0.5, -0.5, 0.5 /)
  call octree_find(tree, x, node, err)
  call assert_equals(0, err, "Query successful.")
  err = 1
  if( size(node%contained) .eq. 1 .and. node%contained(1)%triangle_id .eq. 1 ) err = 0
  call assert_equals(0, err, "Queried a node completely within a triangle.")

  ! A box which is outside the triangle but inside its bounding box
  x = (/ 1.5, 1.5, 0.5 /)
  call octree_find(tree, x, node, err)
  call assert_equals(0, err, "Query successful.")
  if( size(node%contained) .eq.1 .and.node%contained(1)%triangle_id .eq. 1 ) err =0
  call assert_equals(0, err, "Queried a node outside triangle but within the triangle's bounding box.")

  ! This box was already empty in previous iteration, so it should not have been divided
  x = (/ -1.0, -1.0, -1.0 /)
  call octree_find(tree, x, node, err)
  call assert_equals(0, err, "Query successful.")
  call assert_true( all(node%boundary .eq. emptybox), "No children spawned for already empty node.")

  ! Free resources and third initialization which fails as a triangle is outside the box
  call octree_free(tree)

  boundary  = reshape((/ -2, 2, -2, 2, -1, 1 /), shape(boundary))
  max_depth = 1
  call octree_init(triangles, max_depth, boundary, tree, err)
  call assert_equals(1, err, "Initialization fail as a triangle is outside octree.")

  ! Free resources
  call octree_free(tree)

end subroutine test_octree

!> Test collision between line segment and a triangle
subroutine test_triangle_intersection
  type(octree_triangle) :: triangles(4)
  real*8 :: pxyz(3), qxyz(3), rnd(2), r, t, normal(3), iangle, eangle
  logical :: collided, fail
  integer i, j, n

  ! The line segment begins at origo
  pxyz = (/0.d0, 0.d0, 0.d0/)

  ! Now surround origo with triangles so that there is a guaranteed
  ! hit (miss) in positive (neggative) z direction.
  triangles(1)%triangle_id = 1
  triangles(1)%v0 = (/  1.0, -1.0, 0.0 /)
  triangles(1)%v1 = (/  1.0,  1.0, 0.0 /)
  triangles(1)%v2 = (/  0.0,  0.0, 1.0 /)

  triangles(2)%triangle_id = 2
  triangles(2)%v0 = (/  1.0,  1.0, 0.0 /)
  triangles(2)%v1 = (/ -1.0,  1.0, 0.0 /)
  triangles(2)%v2 = (/  0.0,  0.0, 1.0 /)

  triangles(3)%triangle_id = 3
  triangles(3)%v0 = (/ -1.0, -1.0, 0.0 /)
  triangles(3)%v1 = (/  1.0, -1.0, 0.0 /)
  triangles(3)%v2 = (/  0.0,  0.0, 1.0 /)

  triangles(4)%triangle_id = 4
  triangles(4)%v0 = (/ -1.0,  1.0, 0.0 /)
  triangles(4)%v1 = (/ -1.0, -1.0, 0.0 /)
  triangles(4)%v2 = (/  0.0,  0.0, 1.0 /)

!  normals = transpose( reshape(/ -1.0, 0.0, -1.0, 0.0, -1.0, -1.0, 0.0, 1.0, -1.0, 1.0, 0.0, -1.0 /), (/ size(normals, 2), size(normals, 1) /)))

  ! Generate random directions for the line segment and verify that we have
  ! a collision when its z component is positive and miss when it is negative.
  ! Also checks that the angle of incidence is correct.
  fail = .false.
  n = 100
  r = 1.0
  do j=1,n
     call random_number(rnd)
     rnd(1) = TWOPI * rnd(1)
     rnd(2) = PI * rnd(2)
     qxyz = (/ r * sin(rnd(2)) * cos(rnd(1)), r * sin(rnd(2)) * sin(rnd(1)), r * cos(rnd(2)) /)

     collided = .false.
     do i=1,4
        call mod_wall_collision_intersect(pxyz, qxyz, triangles(i)%v0, &
             triangles(i)%v1, triangles(i)%v2, t, iangle)
        if( t .ge. 0.d0 ) then
           collided = .true.
           normal = (/ -( triangles(i)%v0(1) + triangles(i)%v1(1) ) / 2.0, -( triangles(i)%v0(2) + triangles(i)%v1(2) ) / 2.0, -1.0  /)
           normal = normal / sqrt(2.0)
           eangle = PI - acos( dot_product(qxyz, normal) / r )
        end if
     end do
     if( collided .and. ( cos(rnd(2)) .lt. 0.d0 ) ) fail = .true.
     if( collided .and. .not. fail .and. ( abs(iangle - eangle) .gt. 1.e-6 ) ) fail = .true.
  end do
  call assert_false(fail, "Collision check")

end subroutine test_triangle_intersection


end module wall_collision_spec
