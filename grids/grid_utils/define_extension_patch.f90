!> Subroutine defines the new grid_points from crossing of polar and radial coordinate lines
subroutine define_extension_patch(node_list, element_list, newnode_list, newelement_list, n_seg_prev, seg_prev, i_ext)

use constants
use tr_module 
use data_structure
use grid_xpoint_data
use phys_module, only: tokamak_device, n_ext_block, n_wall_blocks, n_wall_block_points_max, n_ext_equidistant, &
                       n_block_points_left, R_block_points_left, Z_block_points_left, corner_block, &
                       n_block_points_right, R_block_points_right, Z_block_points_right, xcase, &
                       n_limiter, R_limiter, Z_limiter
use py_plots_grids

implicit none

! --- Routine parameters
type (type_node_list)       , intent(inout) :: node_list
type (type_element_list)    , intent(inout) :: element_list
type (type_node_list)       , intent(inout) :: newnode_list
type (type_element_list)    , intent(inout) :: newelement_list
integer,                      intent(inout) :: n_seg_prev
real*8 ,                      intent(inout) :: seg_prev(n_seg_max)
integer,                      intent(in)    :: i_ext

! --- local variables
real*8, allocatable :: delta(:)
real*8, allocatable :: R_polar_radial(:,:,:),Z_polar_radial(:,:,:)
real*8, allocatable :: R_polar_sides (:,:,:),Z_polar_sides (:,:,:)
real*8              :: R_polar_left (n_wall_block_points_max,4),Z_polar_left (n_wall_block_points_max,4)
real*8              :: R_polar_right(n_wall_block_points_max,4),Z_polar_right(n_wall_block_points_max,4)
real*8              :: R_polar_bnd  (n_nodes_max/4,          4),Z_polar_bnd  (n_nodes_max/4,          4)
real*8              :: R_polar_wall (n_wall_max,             4),Z_polar_wall (n_wall_max,             4)
integer             :: i, j, k, l, index, i_sep, pieces, i_node, i_node2, count, i_wall
integer             :: n_tmp, n_start, i_start, i_end
integer             :: i_refine, n_refine, i_save
real*8              :: s_refine_min, s_refine_max
real*8              :: s_refine_min_next, s_refine_max_next
real*8              :: R_cub1d(4), Z_cub1d(4)
real*8              :: length
real*8              :: R1,R2,R3,dR1_dr,dR1_ds,dR1_drs,dR1_drr,dR1_dss, dR3_dr
real*8              :: Z1,Z2,Z3,dZ1_dr,dZ1_ds,dZ1_drs,dZ1_drr,dZ1_dss, dZ3_dr
integer             :: n_bnd, index_bnd(n_nodes_max/4), index_bnd_tmp
real*8              :: polar_length, previous_length, sig_tmp, bgf_tmp
real*8              :: alpha1, alpha2, alpha
integer             :: bnd_type_start
integer             :: i_bnd_beg, i_bnd_end
integer             :: i_bnd_beg_prev, i_bnd_end_prev, index_bnd_prev(n_nodes_max/4)
integer             :: i_side_beg, i_side_end, index_side(n_nodes_max/4)
integer             :: n_lim, index_lim(n_nodes_max/4), ier
integer             :: i_lim_beg, i_lim_end
real*8              :: R_lim(n_wall_max), Z_lim(n_wall_max)
integer             :: n_corner_lim, index_corner_lim(n_nodes_max/4)
integer             :: i_corner_lim_beg, i_corner_lim_end
real*8              :: R_corner_lim(n_wall_max), Z_corner_lim(n_wall_max)
integer             :: i_lim_next, i_lim_prev
integer             :: i_node_next, i_node_prev
integer             :: i_elm_beg
integer             :: direction, wall_direction, wall_direction_corner
logical             :: change_direction
real*8              :: st
real*8              :: bgd_radial, sig_radial1
real*8              :: length_seg, length_tmp, length_sum, length_save
real*8              :: length_left, length_right, length_find
real*8              :: length_bottom, length_top, length_prev
real*8              :: diff_min_beg, diff_min_end, diff, diff_with_grid, diff_min
integer             :: i_elm_find(8), i_find
real*8              :: s_find(8), t_find(8)
integer             :: n_nodes, n_nodes_prev, n_nodes_side
integer             :: n_seg, i_seg
real*8, allocatable :: seg(:),      R_seg(:,:),    Z_seg(:,:), seg_tmp(:), seg_new(:)
real*8, allocatable :: seg_bnd(:),  R_seg_bnd(:),  Z_seg_bnd(:)
real*8, allocatable :: seg_wall(:), R_seg_wall(:), Z_seg_wall(:)
real*8, allocatable ::              R_seg_prev(:), Z_seg_prev(:)
real*8, allocatable :: R_dev_bnd(:),     Z_dev_bnd(:)     ! deviation of bnd  from straight line between end points
real*8, allocatable :: R_dev_wall(:),    Z_dev_wall(:)    ! deviation of wall from straight line between end points
real*8, allocatable :: R_deviation(:,:), Z_deviation(:,:) ! deviation average from straight line between end points
character*256       :: plot_filename
character*1         :: char_tmp
character*2         :: char_tmp2
logical, parameter  :: plot_grid = .true.
real*8,  parameter  :: tolerance = 1.d-14
real*8,  parameter  :: side_tolerance = 0.5d-2 ! 0.5cm?
real*8,  parameter  :: wall_node_proximity_tolerance = 0.5d-2 ! 0.5cm?
real*8,  parameter  :: far_from_wall_tolerance = 2.d-2 ! 2cm?
logical             :: far_from_wall, close_to_grid
logical             :: attached, attached_side
logical             :: found_elm, found_smaller
integer             :: element_direction, i_elm, i_elm_save


write(*,*) '*****************************************'
write(*,*) '* X-point grid inside wall :            *'
write(*,*) '*****************************************'
write(*,*) '                 Define extension patch',i_ext


! --- Avoid xpoint nodes
n_start = 4
if (xcase .eq. DOUBLE_NULL) n_start = 8

!-------------------------------- Allocate data structures for new nodes and initialize them
newnode_list%n_nodes = 0
newnode_list%n_dof   = 0
do i = 1, n_nodes_max
  newnode_list%node(i)%x           = 0.d0
  newnode_list%node(i)%values      = 0.d0
  newnode_list%node(i)%deltas      = 0.d0
  newnode_list%node(i)%index       = 0
  newnode_list%node(i)%boundary    = 0
  newnode_list%node(i)%parents     = 0
  newnode_list%node(i)%parent_elem = 0
  newnode_list%node(i)%ref_lambda  = 0.d0
  newnode_list%node(i)%ref_mu      = 0.d0
  newnode_list%node(i)%constrained = .false.
end do

!-------------------------------- Allocate data structures for new elements and initialize them
newelement_list%n_elements = 0
do i = 1, n_elements_max
  newelement_list%element(i)%vertex       = 0
  newelement_list%element(i)%neighbours   = 0
  newelement_list%element(i)%size         = 0.d0
  newelement_list%element(i)%father       = 0
  newelement_list%element(i)%n_sons       = 0
  newelement_list%element(i)%n_gen        = 0
  newelement_list%element(i)%sons         = 0
  newelement_list%element(i)%contain_node = 0
  newelement_list%element(i)%nref         = 0
end do



!------------------------------------------------------------------------------------------------------------------------!
!************************************************************************************************************************!
!************************************************************************************************************************!
!*********************************** First part: find extrapolation points  *********************************************!
!************************************************************************************************************************!
!************************************************************************************************************************!
!------------------------------------------------------------------------------------------------------------------------!



! --- This determines the background of the radial distribution of points (in meshac2 routine)
! --- we use 999. because we want the concentration only on the side of the grid that the patch attaches to
bgd_radial  = 0.6d0
sig_radial1 = 999.!0.3



!---------------------------------------!
!------- Alignment to existing grid ----!
!---------------------------------------!

! --- First, find out which bnd nodes are our starting/ending points
diff_min_beg = 1.d10
diff_min_end = 1.d10
do i_node = 1,node_list%n_nodes
  if (node_list%node(i_node)%boundary .eq. 0) cycle
  diff = sqrt( (node_list%node(i_node)%x(1,1,1)-R_block_points_left(i_ext,1))**2 &
              +(node_list%node(i_node)%x(1,1,2)-Z_block_points_left(i_ext,1))**2 )
  if (diff .lt. diff_min_beg) then
    diff_min_beg = diff
    i_bnd_beg = i_node
  endif
  diff = sqrt( (node_list%node(i_node)%x(1,1,1)-R_block_points_right(i_ext,1))**2 &
              +(node_list%node(i_node)%x(1,1,2)-Z_block_points_right(i_ext,1))**2 )
  if (diff .lt. diff_min_end) then
    diff_min_end = diff
    i_bnd_end = i_node
  endif
enddo
R_block_points_left (i_ext,1) = node_list%node(i_bnd_beg)%x(1,1,1)
Z_block_points_left (i_ext,1) = node_list%node(i_bnd_beg)%x(1,1,2)
R_block_points_right(i_ext,1) = node_list%node(i_bnd_end)%x(1,1,1)
Z_block_points_right(i_ext,1) = node_list%node(i_bnd_end)%x(1,1,2)

! --- Now step along boundary between these two nodes
call find_next_bnd_node(node_list,element_list,i_bnd_beg,-1,i_node_prev)
call find_next_bnd_node(node_list,element_list,i_bnd_beg,+1,i_node_next)

R1 = node_list%node(i_bnd_end)%x(1,1,1)
Z1 = node_list%node(i_bnd_end)%x(1,1,2)
R2 = node_list%node(i_node_prev)%x(1,1,1)
Z2 = node_list%node(i_node_prev)%x(1,1,2)
R3 = node_list%node(i_node_next)%x(1,1,1)
Z3 = node_list%node(i_node_next)%x(1,1,2)

if ( sqrt( (R1-R3)**2 + (Z1-Z3)**2 ) .lt. sqrt( (R1-R2)**2 + (Z1-Z2)**2 ) ) then
  direction = +1
else
  direction = -1
endif

count = 1
index_bnd(1) = i_bnd_beg
i_node = i_bnd_beg
change_direction = .false.
do i=1,node_list%n_nodes
  call find_next_bnd_node(node_list,element_list,i_node,direction,i_node_next)
  count = count + 1
  index_bnd(count) = i_node_next
  if (i_node_next .eq. i_bnd_end) exit
  if ( (node_list%node(i_node_next)%boundary .eq. 3) .and. (count .ge. 2) ) then
    change_direction = .true.
    exit
  endif
  i_node = i_node_next
enddo
n_nodes = count
if (change_direction) then
  direction = -direction
  count = 1
  index_bnd(1) = i_bnd_beg
  i_node = i_bnd_beg
  do i=1,node_list%n_nodes
    call find_next_bnd_node(node_list,element_list,i_node,direction,i_node_next)
    count = count + 1
    index_bnd(count) = i_node_next
    if (i_node_next .eq. i_bnd_end) exit
    if ( (node_list%node(i_node_next)%boundary .eq. 3) .and. (count .ge. 2) ) then
      write(*,*) 'Extended bnd nodes should not have a corner in the middle. Aborting...'
      stop
    endif
    i_node = i_node_next
  enddo
  n_nodes = count
endif

! --- Now, determine which direction our new nodes will have to be
! --- ie. elm1 with nodes (1,2) should be neighbour of elm2 with nodes (4,3) respectively)
found_elm = .false.
do i_elm=1,element_list%n_elements
  do i_node = 1,n_vertex_max
    if (element_list%element(i_elm)%vertex(i_node) .eq. index_bnd(1)) then
      do i_node2 = 1,n_vertex_max
        if (i_node2 .eq. i_node) cycle
        if (element_list%element(i_elm)%vertex(i_node2) .eq. index_bnd(2)) then
          found_elm = .true.
          i_elm_save = i_elm
          i_node_prev = i_node
          i_node_next = i_node2
          exit
        endif
      enddo
      if (found_elm) exit
    endif
  enddo
  if (found_elm) exit
enddo
if (.not. found_elm) then
  write(*,*) 'Could not find corresponding element at boundary. Aborting...'
  stop
else
  if     ( (i_node_prev .eq. 1) .and. (i_node_next .eq. 2) ) then
    element_direction = 1
  elseif ( (i_node_prev .eq. 2) .and. (i_node_next .eq. 1) ) then
    element_direction = 2
  elseif ( (i_node_prev .eq. 2) .and. (i_node_next .eq. 3) ) then
    element_direction = 3
  elseif ( (i_node_prev .eq. 3) .and. (i_node_next .eq. 2) ) then
    element_direction = 4
  elseif ( (i_node_prev .eq. 3) .and. (i_node_next .eq. 4) ) then
    element_direction = 5
  elseif ( (i_node_prev .eq. 4) .and. (i_node_next .eq. 3) ) then
    element_direction = 6
  elseif ( (i_node_prev .eq. 4) .and. (i_node_next .eq. 1) ) then
    element_direction = 7
  elseif ( (i_node_prev .eq. 1) .and. (i_node_next .eq. 4) ) then
    element_direction = 8
  else
    write(*,*) 'Something very wrong here... Aborting...'
    stop
  endif
endif







!---------------------------------------!
!------- Alignment to the wall ---------!
!---------------------------------------!


! --- Second, find out which wall points are our starting/ending points
if (n_wall .eq. 0) then
  n_wall = n_limiter
  R_wall(1:n_wall) = R_limiter(1:n_wall)
  Z_wall(1:n_wall) = Z_limiter(1:n_wall)
  if (n_wall .eq. 0) then
    write(*,*)'Error getting wall data?'
    return
  endif
endif
diff_min_beg = 1.d10
diff_min_end = 1.d10
do i_wall = 1,n_wall
  diff = sqrt( (R_wall(i_wall)-R_block_points_left(i_ext,n_block_points_left(i_ext)))**2 &
              +(Z_wall(i_wall)-Z_block_points_left(i_ext,n_block_points_left(i_ext)))**2 )
  if (diff .lt. diff_min_beg) then
    diff_min_beg = diff
    i_lim_beg = i_wall
  endif
  diff = sqrt( (R_wall(i_wall)-R_block_points_right(i_ext,n_block_points_right(i_ext)))**2 &
              +(Z_wall(i_wall)-Z_block_points_right(i_ext,n_block_points_right(i_ext)))**2 )
  if (diff .lt. diff_min_end) then
    diff_min_end = diff
    i_lim_end = i_wall
  endif
enddo
if (diff_min_beg .lt. wall_node_proximity_tolerance) then
  R_block_points_left (i_ext,n_block_points_left (i_ext)) = R_wall(i_lim_beg)
  Z_block_points_left (i_ext,n_block_points_left (i_ext)) = Z_wall(i_lim_beg)
endif
if (diff_min_end .lt. wall_node_proximity_tolerance) then
  R_block_points_right(i_ext,n_block_points_right(i_ext)) = R_wall(i_lim_end)
  Z_block_points_right(i_ext,n_block_points_right(i_ext)) = Z_wall(i_lim_end)
endif
far_from_wall = .false.
if (     (diff_min_beg .gt. far_from_wall_tolerance) &
    .or. (diff_min_end .gt. far_from_wall_tolerance) ) far_from_wall = .true.

! --- Make sure we are going in right direction
! --- We always take the shortest route! if you need a very large extension that spans almost all
! --- the wall around the whole plasma, then you need to split the extension into several extensions.
! --- Sorry but this is really the most robust way to do it...
count = 1
i_wall = i_lim_beg
direction = +1
do i=1,n_wall
  i_lim_next = i_wall + direction
  if (i_lim_next .gt. n_wall) i_lim_next = 1
  if (i_lim_next .lt. 1     ) i_lim_next = n_wall
  count = count + 1
  if (i_lim_next .eq. i_lim_end) exit
  i_wall = i_lim_next
enddo
n_tmp = count

count = 1
i_wall = i_lim_beg
direction = -1
do i=1,n_wall
  i_lim_next = i_wall + direction
  if (i_lim_next .gt. n_wall) i_lim_next = 1
  if (i_lim_next .lt. 1     ) i_lim_next = n_wall
  count = count + 1
  if (i_lim_next .eq. i_lim_end) exit
  i_wall = i_lim_next
enddo
if (count .lt. n_tmp) then
  direction = -1
else
  direction = +1
endif
wall_direction = direction

count = 1
index_lim(1) = i_lim_beg
i_wall = i_lim_beg
do i=1,n_wall
  i_lim_next = i_wall + direction
  if (i_lim_next .gt. n_wall) i_lim_next = 1
  if (i_lim_next .lt. 1     ) i_lim_next = n_wall
  count = count + 1
  index_lim(count) = i_lim_next
  if (i_lim_next .eq. i_lim_end) exit
  i_wall = i_lim_next
enddo
n_lim = count

! --- We don't allow going all the way around the wall, if this happens, it means our patch is between two
! --- wall nodes, ie. both ends of the patch are closer to a single wall node than any other nodes
if (n_lim .ge. n_wall-1) n_lim = 0

! --- If all points are far from wall, we just ignore the wall...
if (far_from_wall) n_lim = 1

! --- In some cases, we might want to attach the top of the patch to the grid.
! --- This is rare, but useful for cases with holes like the ITER dome.
close_to_grid = .false.
if (far_from_wall) then
  diff_min_beg = 1.d10
  diff_min_end = 1.d10
  do i_node = 1,node_list%n_nodes
    if (node_list%node(i_node)%boundary .eq. 0) cycle
    diff = sqrt( (node_list%node(i_node)%x(1,1,1)-R_block_points_left(i_ext,n_block_points_left(i_ext)))**2 &
                +(node_list%node(i_node)%x(1,1,2)-Z_block_points_left(i_ext,n_block_points_left(i_ext)))**2 )
    if (diff .lt. diff_min_beg) then
      diff_min_beg = diff
      i_lim_beg = i_node
    endif
    diff = sqrt( (node_list%node(i_node)%x(1,1,1)-R_block_points_right(i_ext,n_block_points_right(i_ext)))**2 &
                +(node_list%node(i_node)%x(1,1,2)-Z_block_points_right(i_ext,n_block_points_right(i_ext)))**2 )
    if (diff .lt. diff_min_end) then
      diff_min_end = diff
      i_lim_end = i_node
    endif
  enddo
  if (diff_min_beg .lt. wall_node_proximity_tolerance) then
    R_block_points_left (i_ext,n_block_points_left (i_ext)) = node_list%node(i_lim_beg)%x(1,1,1)
    Z_block_points_left (i_ext,n_block_points_left (i_ext)) = node_list%node(i_lim_beg)%x(1,1,2)
  endif
  if (diff_min_end .lt. wall_node_proximity_tolerance) then
    R_block_points_right(i_ext,n_block_points_right(i_ext)) = node_list%node(i_lim_end)%x(1,1,1)
    Z_block_points_right(i_ext,n_block_points_right(i_ext)) = node_list%node(i_lim_end)%x(1,1,2)
  endif
  close_to_grid = .true.
  if (     (diff_min_beg .gt. wall_node_proximity_tolerance) &
      .or. (diff_min_end .gt. wall_node_proximity_tolerance) ) close_to_grid = .false.
  
  if (close_to_grid) then
    ! --- Now step along boundary between these two nodes
    call find_next_bnd_node(node_list,element_list,i_lim_beg,-1,i_node_prev)
    call find_next_bnd_node(node_list,element_list,i_lim_beg,+1,i_node_next)
    
    R1 = node_list%node(i_lim_end)%x(1,1,1)
    Z1 = node_list%node(i_lim_end)%x(1,1,2)
    R2 = node_list%node(i_node_prev)%x(1,1,1)
    Z2 = node_list%node(i_node_prev)%x(1,1,2)
    R3 = node_list%node(i_node_next)%x(1,1,1)
    Z3 = node_list%node(i_node_next)%x(1,1,2)
    
    if ( sqrt( (R1-R3)**2 + (Z1-Z3)**2 ) .lt. sqrt( (R1-R2)**2 + (Z1-Z2)**2 ) ) then
      direction = +1
    else
      direction = -1
    endif
    
    count = 1
    index_lim(1) = i_lim_beg
    i_node = i_lim_beg
    change_direction = .false.
    do i=1,node_list%n_nodes
      call find_next_bnd_node(node_list,element_list,i_node,direction,i_node_next)
      count = count + 1
      index_lim(count) = i_node_next
      if (i_node_next .eq. i_lim_end) exit
      if ( (node_list%node(i_node_next)%boundary .eq. 3) .and. (count .ge. 2) ) then
        change_direction = .true.
        exit
      endif
      i_node = i_node_next
    enddo
    n_lim = count
    if (change_direction) then
      direction = -direction
      count = 1
      index_lim(1) = i_lim_beg
      i_node = i_lim_beg
      do i=1,node_list%n_nodes
        call find_next_bnd_node(node_list,element_list,i_node,direction,i_node_next)
        count = count + 1
        index_lim(count) = i_node_next
        if (i_node_next .eq. i_lim_end) exit
        if ( (node_list%node(i_node_next)%boundary .eq. 3) .and. (count .ge. 2) ) then
          write(*,*) 'Extended bnd top nodes should not have a corner in the middle. Aborting...'
          stop
        endif
        i_node = i_node_next
      enddo
      n_lim = count
    endif
  else
    n_lim = 1
  endif
  
  if ( (n_lim .gt. 1) .and. (n_lim .ne. n_nodes) ) then
    write(*,*) 'You are trying to join the bottom and top of the patch to two different'
    write(*,*) 'parts of the grid, but the number of nodes (from previous grid) at either'
    write(*,*) 'ends do not match. Please redefine your patch. Aborting...'
    stop
  endif

endif



!------------------------------------------------------------------------------------------------!
!------- Alignment to existing grid on the right side (to connect both left and right sides) ----!
!------------------------------------------------------------------------------------------------!

! --- Are we joining this block with the grid on the right side? With which nodes?
attached_side = .false.
if (.not. close_to_grid) then
  diff_min_beg = 1.d10
  diff_min_end = 1.d10
  do i_node = 1,node_list%n_nodes
    if (node_list%node(i_node)%boundary .eq. 0) cycle
    diff = sqrt( (node_list%node(i_node)%x(1,1,1)-R_block_points_right(i_ext,1))**2 &
                +(node_list%node(i_node)%x(1,1,2)-Z_block_points_right(i_ext,1))**2 )
    if (diff .lt. diff_min_beg) then
      diff_min_beg = diff
      i_side_beg = i_node
    endif
    diff = sqrt( (node_list%node(i_node)%x(1,1,1)-R_block_points_right(i_ext,n_block_points_right(i_ext)))**2 &
                +(node_list%node(i_node)%x(1,1,2)-Z_block_points_right(i_ext,n_block_points_right(i_ext)))**2 )
    if (diff .lt. diff_min_end) then
      diff_min_end = diff
      i_side_end = i_node
    endif
  enddo
  if (      (diff_min_beg .lt. side_tolerance) &
      .and. (diff_min_end .lt. side_tolerance) &
      .and. (n_block_points_right(i_ext) .eq. 2) ) attached_side = .true.
  ! --- If yes, step along boundary between these two nodes
  if (attached_side) then
    ! --- Precise beg/end points
    R_block_points_right(i_ext,1) = node_list%node(i_side_beg)%x(1,1,1)
    Z_block_points_right(i_ext,1) = node_list%node(i_side_beg)%x(1,1,2)
    R_block_points_right(i_ext,n_block_points_right(i_ext)) = node_list%node(i_side_end)%x(1,1,1)
    Z_block_points_right(i_ext,n_block_points_right(i_ext)) = node_list%node(i_side_end)%x(1,1,2)
    ! --- Now step along boundary between these two nodes
    call find_next_bnd_node(node_list,element_list,i_side_beg,-1,i_node_prev)
    call find_next_bnd_node(node_list,element_list,i_side_beg,+1,i_node_next)
    R1 = node_list%node(i_side_end)%x(1,1,1)
    Z1 = node_list%node(i_side_end)%x(1,1,2)
    R2 = node_list%node(i_node_prev)%x(1,1,1)
    Z2 = node_list%node(i_node_prev)%x(1,1,2)
    R3 = node_list%node(i_node_next)%x(1,1,1)
    Z3 = node_list%node(i_node_next)%x(1,1,2)
    if ( sqrt( (R1-R3)**2 + (Z1-Z3)**2 ) .lt. sqrt( (R1-R2)**2 + (Z1-Z2)**2 ) ) then
      direction = +1
    else
      direction = -1
    endif
    ! --- In the main direction
    count = 1
    index_side(1) = i_side_beg
    i_node = i_side_beg
    change_direction = .false.
    do i=1,node_list%n_nodes
      call find_next_bnd_node(node_list,element_list,i_node,direction,i_node_next)
      count = count + 1
      index_side(count) = i_node_next
      if (i_node_next .eq. i_side_end) exit
      if ( (node_list%node(i_node_next)%boundary .eq. 3) .and. (count .ge. 2) ) then
        change_direction = .true.
        exit
      endif
      i_node = i_node_next
    enddo
    n_nodes_side = count
    ! --- Try the other direction if we failed
    if (change_direction) then
      direction = -direction
      count = 1
      index_side(1) = i_side_beg
      i_node = i_side_beg
      do i=1,node_list%n_nodes
        call find_next_bnd_node(node_list,element_list,i_node,direction,i_node_next)
        count = count + 1
        index_side(count) = i_node_next
        if (i_node_next .eq. i_side_end) exit
        if ( (node_list%node(i_node_next)%boundary .eq. 3) .and. (count .ge. 2) ) then
          write(*,*) 'Extended bnd nodes on right side should not'
          write(*,*) 'have a corner in the middle. Aborting...'
          stop
        endif
        i_node = i_node_next
      enddo
      n_nodes_side = count
    endif
  endif
endif






!-------------------------------------------------------------------------!
!------- Alignment to the wall on "left" side (for corners only) ---------!
!-------------------------------------------------------------------------!


! --- Then, find out which wall points are our starting/ending points, if we want to align on the left side as well
if (corner_block(i_ext) .eq. 1) then
  diff_min_beg = 1.d10
  diff_min_end = 1.d10
  do i_wall = 1,n_wall
    diff = sqrt( (R_wall(i_wall)-R_block_points_left(i_ext,1))**2 &
                +(Z_wall(i_wall)-Z_block_points_left(i_ext,1))**2 )
    if (diff .lt. diff_min_beg) then
      diff_min_beg = diff
      i_corner_lim_beg = i_wall
    endif
    diff = sqrt( (R_wall(i_wall)-R_block_points_left(i_ext,n_block_points_left(i_ext)))**2 &
                +(Z_wall(i_wall)-Z_block_points_left(i_ext,n_block_points_left(i_ext)))**2 )
    if (diff .lt. diff_min_end) then
      diff_min_end = diff
      i_corner_lim_end = i_wall
    endif
  enddo
  if (diff_min_beg .lt. wall_node_proximity_tolerance) then
    R_block_points_left (i_ext,1) = R_wall(i_corner_lim_beg)
    Z_block_points_left (i_ext,1) = Z_wall(i_corner_lim_beg)
  endif
  if (diff_min_end .lt. wall_node_proximity_tolerance) then
    R_block_points_left (i_ext,n_block_points_left (i_ext)) = R_wall(i_corner_lim_end)
    Z_block_points_left (i_ext,n_block_points_left (i_ext)) = Z_wall(i_corner_lim_end)
  endif
  far_from_wall = .false.
  if (     (diff_min_beg .gt. far_from_wall_tolerance) &
      .or. (diff_min_end .gt. far_from_wall_tolerance) ) far_from_wall = .true.
  
  ! --- Make sure we are going in right direction
  ! --- We always take the shortest route! if you need a very large extension that spans almost all
  ! --- the wall around the whole plasma, then you need to split the extension into several extensions.
  ! --- Sorry but this is really the most robust way to do it...
  count = 1
  i_wall = i_corner_lim_beg
  direction = +1
  do i=1,n_wall
    i_lim_next = i_wall + direction
    if (i_lim_next .gt. n_wall) i_lim_next = 1
    if (i_lim_next .lt. 1     ) i_lim_next = n_wall
    count = count + 1
    if (i_lim_next .eq. i_corner_lim_end) exit
    i_wall = i_lim_next
  enddo
  n_tmp = count
  
  count = 1
  i_wall = i_corner_lim_beg
  direction = -1
  do i=1,n_wall
    i_lim_next = i_wall + direction
    if (i_lim_next .gt. n_wall) i_lim_next = 1
    if (i_lim_next .lt. 1     ) i_lim_next = n_wall
    count = count + 1
    if (i_lim_next .eq. i_corner_lim_end) exit
    i_wall = i_lim_next
  enddo
  if (count .lt. n_tmp) then
    direction = -1
  else
    direction = +1
  endif
  wall_direction_corner = direction
  
  count = 1
  index_corner_lim(1) = i_corner_lim_beg
  i_wall = i_corner_lim_beg
  do i=1,n_wall
    i_lim_next = i_wall + direction
    if (i_lim_next .gt. n_wall) i_lim_next = 1
    if (i_lim_next .lt. 1     ) i_lim_next = n_wall
    count = count + 1
    index_corner_lim(count) = i_lim_next
    if (i_lim_next .eq. i_corner_lim_end) exit
    i_wall = i_lim_next
  enddo
  n_corner_lim = count
  
  ! --- We don't allow going all the way around the wall, if this happens, it means our patch is between two
  ! --- wall nodes, ie. both ends of the patch are closer to a single wall node than any other nodes
  if (n_corner_lim .ge. n_wall-1) n_corner_lim = 0
  
  ! --- If all points are far from wall, we just ignore the wall...
  if (far_from_wall) n_corner_lim = 1
endif





!------------------------------------------------------!
!------- Two consecutive patches are attached ---------!
!------------------------------------------------------!


! --- Are we joining this block with the previous one? (note: corner block are not considered)
attached = .false.
if ( (i_ext .gt. 1) .and. (corner_block(i_ext) .ne. 1) ) then
  diff = sqrt( (R_block_points_left(i_ext,1)-R_block_points_right(i_ext-1,1))**2 + (Z_block_points_left(i_ext,1)-Z_block_points_right(i_ext-1,1))**2 )
  ! --- Check difference with any existing grid
  diff_min_beg = 1.d10
  diff_min_end = 1.d10
  if (.not. close_to_grid) then
    do i_node = 1,node_list%n_nodes
      if (node_list%node(i_node)%boundary .eq. 0) cycle
      diff_with_grid = sqrt( (node_list%node(i_node)%x(1,1,1)-R_block_points_left(i_ext,1))**2 &
                            +(node_list%node(i_node)%x(1,1,2)-Z_block_points_left(i_ext,1))**2 )
      if (diff_with_grid .lt. diff_min_beg) then
        diff_min_beg = diff_with_grid
        i_bnd_beg_prev = i_node
      endif
      n_tmp = n_block_points_left(i_ext)
      diff_with_grid = sqrt( (node_list%node(i_node)%x(1,1,1)-R_block_points_left(i_ext,n_tmp))**2 &
                            +(node_list%node(i_node)%x(1,1,2)-Z_block_points_left(i_ext,n_tmp))**2 )
      if (diff_with_grid .lt. diff_min_end) then
        diff_min_end = diff_with_grid
        i_bnd_end_prev = i_node
      endif
    enddo
  endif
  diff_with_grid = diff_min_beg + diff_min_end
  ! --- If we're attached, save connection nodes
  if ( (diff .lt. tolerance) .or. (diff_with_grid .lt. side_tolerance) ) then
    attached = .true.
    ! --- Determine which nodes are attached with previous block (this was already done for diff_with_grid case)
    ! --- First, find out which bnd nodes are our starting/ending points for previous block
    if (diff_with_grid .gt. side_tolerance) then
      attached = .true.
      ! --- Determine which nodes are attached with previous block
      ! --- First, find out which bnd nodes are our starting/ending points for previous block
      diff_min_beg = 1.d10
      diff_min_end = 1.d10
      do i_node = 1,node_list%n_nodes
        if (node_list%node(i_node)%boundary .eq. 0) cycle
        diff = sqrt( (node_list%node(i_node)%x(1,1,1)-R_block_points_right(i_ext-1,1))**2 &
                    +(node_list%node(i_node)%x(1,1,2)-Z_block_points_right(i_ext-1,1))**2 )
        if (diff .lt. diff_min_beg) then
          diff_min_beg = diff
          i_bnd_beg_prev = i_node
        endif
        n_tmp = n_block_points_right(i_ext-1)
        diff = sqrt( (node_list%node(i_node)%x(1,1,1)-R_block_points_right(i_ext-1,n_tmp))**2 &
                    +(node_list%node(i_node)%x(1,1,2)-Z_block_points_right(i_ext-1,n_tmp))**2 )
        if (diff .lt. diff_min_end) then
          diff_min_end = diff
          i_bnd_end_prev = i_node
        endif
      enddo
    endif
    
    ! --- Now step along boundary between these two nodes
    call find_next_bnd_node(node_list,element_list,i_bnd_beg_prev,-1,i_node_prev)
    call find_next_bnd_node(node_list,element_list,i_bnd_beg_prev,+1,i_node_next)
    R1 = node_list%node(i_bnd_end_prev)%x(1,1,1)
    Z1 = node_list%node(i_bnd_end_prev)%x(1,1,2)
    R2 = node_list%node(i_node_prev)%x(1,1,1)
    Z2 = node_list%node(i_node_prev)%x(1,1,2)
    R3 = node_list%node(i_node_next)%x(1,1,1)
    Z3 = node_list%node(i_node_next)%x(1,1,2)
    if ( sqrt( (R1-R3)**2 + (Z1-Z3)**2 ) .lt. sqrt( (R1-R2)**2 + (Z1-Z2)**2 ) ) then
      direction = +1
      bnd_type_start = node_list%node(i_node_next)%boundary
    else
      direction = -1
      bnd_type_start = node_list%node(i_node_prev)%boundary
    endif
    
    count = 1
    index_bnd_prev(1) = i_bnd_beg_prev
    i_node = i_bnd_beg_prev
    change_direction = .false.
    do i=1,node_list%n_nodes
      call find_next_bnd_node(node_list,element_list,i_node,direction,i_node_next)
      count = count + 1
      index_bnd_prev(count) = i_node_next
      if (i_node_next .eq. i_bnd_end_prev) exit
      if ( (node_list%node(i_node)%boundary .ne. bnd_type_start) .and. (count .gt. 2) ) then
        change_direction = .true.
        exit
      endif
      i_node = i_node_next
    enddo
    n_nodes_prev = count
    if (change_direction) then
      direction = -direction
      count = 1
      index_bnd_prev(1) = i_bnd_beg_prev
      i_node = i_bnd_beg_prev
      do i=1,node_list%n_nodes
        call find_next_bnd_node(node_list,element_list,i_node,direction,i_node_next)
        count = count + 1
        index_bnd_prev(count) = i_node_next
        if (i_node_next .eq. i_bnd_end_prev) exit
        if ( (node_list%node(i_node_next)%boundary .eq. 3) .and. (count .ge. 2) ) then
          write(*,*) 'Extended bnd nodes from previous block'
          write(*,*) 'should not have a corner in the middle. Aborting...'
          stop
        endif
        i_node = i_node_next
      enddo
      n_nodes_prev = count
    endif
    if (diff_with_grid .lt. side_tolerance) n_seg_prev = n_nodes_prev
    if (n_nodes_prev .ne. n_seg_prev) then
      write(*,*)'Problem: number of common nodes between blocks does not match'
      stop
    endif
    
    ! ---  Which points were on previous block side?
    if (diff_with_grid .gt. side_tolerance) then
      n_tmp = n_block_points_right(i_ext-1)
      call create_polar_lines_simple(n_tmp, R_block_points_right(i_ext-1,1:n_tmp), Z_block_points_right(i_ext-1,1:n_tmp), R_polar_right(1:n_tmp-1,1:4), Z_polar_right(1:n_tmp-1,1:4))
      length_prev = 0.d0
      do i=1,n_block_points_right(i_ext-1)-1
        call from_polar_to_cubic(R_polar_right(i,1:4),R_cub1d)
        call from_polar_to_cubic(Z_polar_right(i,1:4),Z_cub1d)
        call curve_length(R_cub1d(1), R_cub1d(2), R_cub1d(3), R_cub1d(4), &
                          Z_cub1d(1), Z_cub1d(2), Z_cub1d(3), Z_cub1d(4), -1.d0, 1.d0, length)
        length_prev = length_prev + length
      enddo
      allocate(R_seg_prev(n_seg_prev), Z_seg_prev(n_seg_prev))
      do i=1,n_seg_prev
        i_node = index_bnd_prev(i)
        R_seg_prev(i) = node_list%node(i_node)%x(1,1,1)
        Z_seg_prev(i) = node_list%node(i_node)%x(1,1,2)
      enddo
    else
      allocate(R_seg_prev(n_seg_prev), Z_seg_prev(n_seg_prev))
      do i=1,n_seg_prev
        i_node = index_bnd_prev(i)
        R_seg_prev(i) = node_list%node(i_node)%x(1,1,1)
        Z_seg_prev(i) = node_list%node(i_node)%x(1,1,2)
      enddo
      length_prev = 0.d0
      seg_prev(1) = 0.d0
      do i=2,n_seg_prev
        length = sqrt( (R_seg_prev(i)-R_seg_prev(i-1))**2 + (Z_seg_prev(i)-Z_seg_prev(i-1))**2 )
        length_prev = length_prev + length
        seg_prev(i) = length_prev
      enddo
      seg_prev(1:n_seg_prev) = seg_prev(1:n_seg_prev) / length_prev
    endif
    
  endif
endif







!------------------------------------------------------!
!------- Number of points along left/right sides ------!
!------------------------------------------------------!


! --- Get length of our radial segments on each side
n_tmp = n_block_points_left (i_ext)
call create_polar_lines_simple(n_tmp, R_block_points_left (i_ext,1:n_tmp), Z_block_points_left (i_ext,1:n_tmp), R_polar_left(1:n_tmp-1,1:4) , Z_polar_left(1:n_tmp-1,1:4) )
n_tmp = n_block_points_right(i_ext)
call create_polar_lines_simple(n_tmp, R_block_points_right(i_ext,1:n_tmp), Z_block_points_right(i_ext,1:n_tmp), R_polar_right(1:n_tmp-1,1:4), Z_polar_right(1:n_tmp-1,1:4))
length_right = 0.d0
do i=1,n_block_points_right(i_ext)-1
  call from_polar_to_cubic(R_polar_right(i,1:4),R_cub1d)
  call from_polar_to_cubic(Z_polar_right(i,1:4),Z_cub1d)
  call curve_length(R_cub1d(1), R_cub1d(2), R_cub1d(3), R_cub1d(4), &
                    Z_cub1d(1), Z_cub1d(2), Z_cub1d(3), Z_cub1d(4), -1.d0, 1.d0, length)
  length_right = length_right + length
enddo
length_left = 0.d0
do i=1,n_block_points_left(i_ext)-1
  call from_polar_to_cubic(R_polar_left(i,1:4),R_cub1d)
  call from_polar_to_cubic(Z_polar_left(i,1:4),Z_cub1d)
  call curve_length(R_cub1d(1), R_cub1d(2), R_cub1d(3), R_cub1d(4), &
                    Z_cub1d(1), Z_cub1d(2), Z_cub1d(3), Z_cub1d(4), -1.d0, 1.d0, length)
  length_left = length_left + length
enddo

! --- How many points do we have in radial direction?
if (.not. attached) then
  n_seg = n_ext_block(i_ext)
else
  if (diff_with_grid .gt. side_tolerance) then
    if (n_block_points_left(i_ext) .gt. n_block_points_right(i_ext-1)) then
      if (n_ext_block(i_ext) .le. n_seg_prev+2) then 
        length_tmp  = (1.d0 - seg_prev(n_seg_prev-1)) * length_prev
        length_find = length_left - length_prev
        n_seg = n_seg_prev + int(length_find/length_tmp) + 1
      else
        n_seg = n_ext_block(i_ext)
      endif
    elseif (n_block_points_left(i_ext) .eq. n_block_points_right(i_ext-1)) then
      n_seg = n_seg_prev
    else
      ! ---  Which of these previous block points are we using?
      count = 0
      do i_seg = 1,n_seg_prev
        length_find = seg_prev(i_seg) * length_prev
        if (length_find .gt. length_left) exit
        count = count + 1
      enddo
      n_seg = count
    endif
  else
    n_seg = n_seg_prev
  endif
endif
! --- Check that we have the same number on the left and right if we are attaching both ends
if (attached_side) then
  if (n_nodes_side .ne. n_seg) then
    write(*,*)'You are trying to join both the left and right sides'
    write(*,*)'but the number of nodes on left and right sides are different'
    write(*,*)'left/right node numbers:',n_seg,n_nodes_side
    write(*,*)'Aborting...'
    stop
  endif
endif






!------------------------------------------------------!
!------- Allocate all points --------------------------!
!------------------------------------------------------!


! --- Now we know how many points our grid will have: n_nodes along the grid edge, and n_seg in the other direction, allocate data
allocate(seg(n_seg), R_seg(n_seg,n_nodes), Z_seg(n_seg,n_nodes), seg_tmp(n_seg), seg_new(n_seg))
allocate(seg_bnd (n_nodes),  R_seg_bnd (n_nodes),  Z_seg_bnd (n_nodes))
allocate(seg_wall(n_nodes),  R_seg_wall(n_nodes),  Z_seg_wall(n_nodes))
allocate(R_dev_bnd(n_nodes),   Z_dev_bnd(n_nodes))
allocate(R_dev_wall(n_nodes),  Z_dev_wall(n_nodes))
allocate(R_deviation(n_seg,n_nodes), Z_deviation(n_seg,n_nodes))
allocate(R_polar_radial(n_nodes-1,4,n_seg  ),Z_polar_radial(n_nodes-1,4,n_seg  ))
allocate(R_polar_sides (n_seg-1  ,4,n_nodes),Z_polar_sides (n_seg-1  ,4,n_nodes))




!------------------------------------------------------!
!------- Segmentation in both directions --------------!
!------------------------------------------------------!
seg = 0.d0

! --- Segmentation in radial direction is the input (note this will get changed later on for smooth transition)
if (attached) then
  if ( (n_block_points_left(i_ext) .gt. n_block_points_right(i_ext-1)) .and. (diff_with_grid .gt. side_tolerance) ) then
    do i_seg = 1,n_seg_prev
      seg(i_seg) = (seg_prev(i_seg) * length_prev) / length_left
    enddo
    do i_seg = n_seg_prev+1,n_seg
      seg(i_seg) = seg(n_seg_prev) + real(i_seg-n_seg_prev)/real(n_seg-n_seg_prev) * (1.d0 - seg(n_seg_prev))
    enddo
  else
    seg(1:n_seg) = seg_prev(1:n_seg) / seg_prev(n_seg)
  endif
else
  call meshac2(n_seg,seg,1.d0,9999.d0,sig_radial1,9999.d0,bgd_radial,1.0d0)
endif

! --- Segmentation in other direction will be given by bnd nodes
do i=1,n_nodes
  i_node = index_bnd(i)
  R_seg_bnd(i) = node_list%node(i_node)%x(1,1,1)
  Z_seg_bnd(i) = node_list%node(i_node)%x(1,1,2)
enddo
call create_polar_lines_simple(n_nodes, R_seg_bnd(1:n_nodes), Z_seg_bnd(1:n_nodes), R_polar_bnd(1:n_nodes-1,1:4) , Z_polar_bnd(1:n_nodes-1,1:4) )
length_bottom = 0.d0
do i=1,n_nodes-1
  call from_polar_to_cubic(R_polar_bnd(i,1:4),R_cub1d)
  call from_polar_to_cubic(Z_polar_bnd(i,1:4),Z_cub1d)
  call curve_length(R_cub1d(1), R_cub1d(2), R_cub1d(3), R_cub1d(4), &
                    Z_cub1d(1), Z_cub1d(2), Z_cub1d(3), Z_cub1d(4), -1.d0, 1.d0, length)
  length_bottom = length_bottom + length
enddo
length_tmp = 0.d0
seg_bnd(1) = 0.d0
do i=1,n_nodes-1
  call from_polar_to_cubic(R_polar_bnd(i,1:4),R_cub1d)
  call from_polar_to_cubic(Z_polar_bnd(i,1:4),Z_cub1d)
  call curve_length(R_cub1d(1), R_cub1d(2), R_cub1d(3), R_cub1d(4), &
                    Z_cub1d(1), Z_cub1d(2), Z_cub1d(3), Z_cub1d(4), -1.d0, 1.d0, length)
  length_tmp = length_tmp + length
  seg_bnd(i+1) = length_tmp / length_bottom
enddo
seg_bnd(n_nodes) = 1.d0
seg_wall(1:n_nodes) = seg_bnd(1:n_nodes)
! --- Segmentation might be different at the top if joining to another part of the grid
if (close_to_grid) then
  do i=1,n_lim
    i_node = index_lim(i)
    R_seg_wall(i) = node_list%node(i_node)%x(1,1,1)
    Z_seg_wall(i) = node_list%node(i_node)%x(1,1,2)
  enddo
  call create_polar_lines_simple(n_lim, R_seg_wall(1:n_lim), Z_seg_wall(1:n_lim), R_polar_wall(1:n_lim-1,1:4) , Z_polar_wall(1:n_lim-1,1:4) )
  length_top = 0.d0
  do i=1,n_lim-1
    call from_polar_to_cubic(R_polar_wall(i,1:4),R_cub1d)
    call from_polar_to_cubic(Z_polar_wall(i,1:4),Z_cub1d)
    call curve_length(R_cub1d(1), R_cub1d(2), R_cub1d(3), R_cub1d(4), &
                      Z_cub1d(1), Z_cub1d(2), Z_cub1d(3), Z_cub1d(4), -1.d0, 1.d0, length)
    length_top = length_top + length
  enddo
  length_tmp = 0.d0
  seg_wall(1) = 0.d0
  do i=1,n_lim-1
    call from_polar_to_cubic(R_polar_wall(i,1:4),R_cub1d)
    call from_polar_to_cubic(Z_polar_wall(i,1:4),Z_cub1d)
    call curve_length(R_cub1d(1), R_cub1d(2), R_cub1d(3), R_cub1d(4), &
                      Z_cub1d(1), Z_cub1d(2), Z_cub1d(3), Z_cub1d(4), -1.d0, 1.d0, length)
    length_tmp = length_tmp + length
    seg_wall(i+1) = length_tmp / length_top
  enddo
  seg_wall(n_lim) = 1.d0
endif






!------------------------------------------------------!
!------- Points for polar line along the wall ---------!
!------------------------------------------------------!


! --- Find intermediate points along grid bnd aligned to the wall
if (.not. close_to_grid) then
  if (n_lim .gt. 1) then
    do i=1,n_lim
      R_lim(i) = R_wall(index_lim(i))
      Z_lim(i) = Z_wall(index_lim(i))
    enddo
    ! --- We add a point at the beg/end if we are far from any wall point (end)
    i_lim_prev = index_lim(n_lim) - wall_direction
    if (i_lim_prev .gt. n_wall) i_lim_prev = 1
    if (i_lim_prev .lt. 1     ) i_lim_prev = n_wall
    diff         = sqrt( (R_wall(index_lim(n_lim))-R_block_points_right(i_ext,n_block_points_right(i_ext)))**2 &
                        +(Z_wall(index_lim(n_lim))-Z_block_points_right(i_ext,n_block_points_right(i_ext)))**2 )
    if (diff .gt. wall_node_proximity_tolerance) then
      diff_min_beg = sqrt( (R_wall(index_lim(n_lim))-R_wall(i_lim_prev))**2 &
                          +(Z_wall(index_lim(n_lim))-Z_wall(i_lim_prev))**2 )
      diff         = sqrt( (R_wall(i_lim_prev)  -R_block_points_right(i_ext,n_block_points_right(i_ext)))**2 &
                          +(Z_wall(i_lim_prev)  -Z_block_points_right(i_ext,n_block_points_right(i_ext)))**2 )
      if (diff .gt. diff_min_beg) n_lim = n_lim + 1
      R_lim(n_lim) = R_block_points_right(i_ext,n_block_points_right(i_ext))
      Z_lim(n_lim) = Z_block_points_right(i_ext,n_block_points_right(i_ext))
    endif
    ! --- We add a point at the beg/end if we are far from any wall point (beg)
    i_lim_prev = index_lim(1) - wall_direction
    if (i_lim_prev .gt. n_wall) i_lim_prev = 1
    if (i_lim_prev .lt. 1     ) i_lim_prev = n_wall
    diff         = sqrt( (R_wall(index_lim(1))-R_block_points_left(i_ext,n_block_points_left(i_ext)))**2 &
                        +(Z_wall(index_lim(1))-Z_block_points_left(i_ext,n_block_points_left(i_ext)))**2 )
    if (diff .gt. wall_node_proximity_tolerance) then
      diff_min_beg = sqrt( (R_wall(index_lim(1))-R_wall(i_lim_prev))**2 &
                          +(Z_wall(index_lim(1))-Z_wall(i_lim_prev))**2 )
      diff         = sqrt( (R_wall(i_lim_prev)  -R_block_points_left(i_ext,n_block_points_left(i_ext)))**2 &
                          +(Z_wall(i_lim_prev)  -Z_block_points_left(i_ext,n_block_points_left(i_ext)))**2 )
      if (diff .lt. diff_min_beg) then
        n_lim = n_lim + 1
        do i = n_lim,2,-1
          R_lim(i) = R_lim(i-1)
          Z_lim(i) = Z_lim(i-1)
        enddo
      endif
      R_lim(1) = R_block_points_left(i_ext,n_block_points_left(i_ext))
      Z_lim(1) = Z_block_points_left(i_ext,n_block_points_left(i_ext))
    endif
  else
    n_lim = 2
    R_lim(1) = R_block_points_left(i_ext,n_block_points_left(i_ext))
    Z_lim(1) = Z_block_points_left(i_ext,n_block_points_left(i_ext))
    R_lim(2) = R_block_points_right(i_ext,n_block_points_right(i_ext))
    Z_lim(2) = Z_block_points_right(i_ext,n_block_points_right(i_ext))
  endif
  ! --- Make sure we properly connect to the previous patch
  if ( (attached) .and. (diff_with_grid .gt. side_tolerance) ) then
    if (n_block_points_left(i_ext) .le. n_block_points_right(i_ext-1)) then
      R_lim(1) = R_seg_prev(n_seg)
      Z_lim(1) = Z_seg_prev(n_seg)
    endif
  endif
  call create_polar_lines_simple(n_lim, R_lim(1:n_lim), Z_lim(1:n_lim), R_polar_wall(1:n_lim-1,1:4) , Z_polar_wall(1:n_lim-1,1:4) )
  length_top = 0.d0
  do i=1,n_lim-1
    call from_polar_to_cubic(R_polar_wall(i,1:4),R_cub1d)
    call from_polar_to_cubic(Z_polar_wall(i,1:4),Z_cub1d)
    call curve_length(R_cub1d(1), R_cub1d(2), R_cub1d(3), R_cub1d(4), &
                      Z_cub1d(1), Z_cub1d(2), Z_cub1d(3), Z_cub1d(4), -1.d0, 1.d0, length)
    length_top = length_top + length
  enddo
endif




!------------------------------------------------------!
!------- Segment points along the wall ----------------!
!------------------------------------------------------!

if (.not. close_to_grid) then
  R_seg_wall(1) = R_polar_wall(1,1)
  Z_seg_wall(1) = Z_polar_wall(1,1)
  length_seg = 0.d0
  do i_seg = 2,n_nodes-1
    length_find = seg_wall(i_seg)
    length_sum  = 0.d0
    do i=1,n_lim-1
      call from_polar_to_cubic(R_polar_wall(i,1:4),R_cub1d)
      call from_polar_to_cubic(Z_polar_wall(i,1:4),Z_cub1d)
      call curve_length(R_cub1d(1), R_cub1d(2), R_cub1d(3), R_cub1d(4), &
                        Z_cub1d(1), Z_cub1d(2), Z_cub1d(3), Z_cub1d(4), -1.d0, 1.d0, length)
      if ((length_sum + length)/length_top .ge. length_find) then
        length_tmp = length_find*length_top - length_sum
        st = 2.d0 * length_tmp/length - 1.d0
        call CUB1D(R_cub1d(1), R_cub1d(2), R_cub1d(3), R_cub1d(4), st, R3, dR3_dr)
        call CUB1D(Z_cub1d(1), Z_cub1d(2), Z_cub1d(3), Z_cub1d(4), st, Z3, dZ3_dr)
        R_seg_wall(i_seg) = R3
        Z_seg_wall(i_seg) = Z3
        exit
      endif
      length_sum = length_sum + length
    enddo
  enddo
  R_seg_wall(n_nodes) = R_polar_wall(n_lim-1,4)
  Z_seg_wall(n_nodes) = Z_polar_wall(n_lim-1,4)
endif




!------------------------------------------------------!
!------- Segment points along the sides ---------------!
!------------------------------------------------------!


! --- Find intermediate points along sides
R_seg(1,1)       = R_polar_left (1,1)
Z_seg(1,1)       = Z_polar_left (1,1)
R_seg(1,n_nodes) = R_polar_right(1,1)
Z_seg(1,n_nodes) = Z_polar_right(1,1)
length_seg = 0.d0
do i_seg = 2,n_seg-1
  length_find = seg(i_seg)
  ! --- Left side
  length_sum  = 0.d0
  do i=1,n_block_points_left(i_ext)-1
    call from_polar_to_cubic(R_polar_left(i,1:4),R_cub1d)
    call from_polar_to_cubic(Z_polar_left(i,1:4),Z_cub1d)
    call curve_length(R_cub1d(1), R_cub1d(2), R_cub1d(3), R_cub1d(4), &
                      Z_cub1d(1), Z_cub1d(2), Z_cub1d(3), Z_cub1d(4), -1.d0, 1.d0, length)
    if ((length_sum + length)/length_left .ge. length_find) then
      length_tmp = length_find*length_left - length_sum
      st = 2.d0 * length_tmp/length - 1.d0
      call CUB1D(R_cub1d(1), R_cub1d(2), R_cub1d(3), R_cub1d(4), st, R3, dR3_dr)
      call CUB1D(Z_cub1d(1), Z_cub1d(2), Z_cub1d(3), Z_cub1d(4), st, Z3, dZ3_dr)
      R_seg(i_seg,1) = R3
      Z_seg(i_seg,1) = Z3
      exit
    endif
    length_sum = length_sum + length
  enddo
  ! --- Right side
  length_sum  = 0.d0
  do i=1,n_block_points_right(i_ext)-1
    call from_polar_to_cubic(R_polar_right(i,1:4),R_cub1d)
    call from_polar_to_cubic(Z_polar_right(i,1:4),Z_cub1d)
    call curve_length(R_cub1d(1), R_cub1d(2), R_cub1d(3), R_cub1d(4), &
                      Z_cub1d(1), Z_cub1d(2), Z_cub1d(3), Z_cub1d(4), -1.d0, 1.d0, length)
    if ((length_sum + length)/length_right .ge. length_find) then
      length_tmp = length_find*length_right - length_sum
      st = 2.d0 * length_tmp/length - 1.d0
      call CUB1D(R_cub1d(1), R_cub1d(2), R_cub1d(3), R_cub1d(4), st, R3, dR3_dr)
      call CUB1D(Z_cub1d(1), Z_cub1d(2), Z_cub1d(3), Z_cub1d(4), st, Z3, dZ3_dr)
      R_seg(i_seg,n_nodes) = R3
      Z_seg(i_seg,n_nodes) = Z3
      exit
    endif
    length_sum = length_sum + length
  enddo
enddo
R_seg(n_seg,1)       = R_polar_left (n_block_points_left (i_ext)-1,4)
Z_seg(n_seg,1)       = Z_polar_left (n_block_points_left (i_ext)-1,4)
R_seg(n_seg,n_nodes) = R_polar_right(n_block_points_right(i_ext)-1,4)
Z_seg(n_seg,n_nodes) = Z_polar_right(n_block_points_right(i_ext)-1,4)

! --- Align beginning of radial segments to the existing grid
do i=1,n_nodes
  R_seg(1,i) = R_seg_bnd(i)
  Z_seg(1,i) = Z_seg_bnd(i)
enddo

! --- Align right side to existing grid if attached
if (attached_side) then
  do i=1,n_nodes_side
    i_node = index_side(i)
    R_seg(i,n_nodes) = node_list%node(i_node)%x(1,1,1)
    Z_seg(i,n_nodes) = node_list%node(i_node)%x(1,1,2)
  enddo
  call create_polar_lines_simple(n_nodes_side, R_seg(1:n_nodes_side,n_nodes), Z_seg(1:n_nodes_side,n_nodes), R_polar_right(1:n_nodes_side-1,1:4) , Z_polar_right(1:n_nodes_side-1,1:4) )
  length_right = 0.d0
  do i=1,n_nodes_side-1
    call from_polar_to_cubic(R_polar_right(i,1:4),R_cub1d)
    call from_polar_to_cubic(Z_polar_right(i,1:4),Z_cub1d)
    call curve_length(R_cub1d(1), R_cub1d(2), R_cub1d(3), R_cub1d(4), &
                      Z_cub1d(1), Z_cub1d(2), Z_cub1d(3), Z_cub1d(4), -1.d0, 1.d0, length)
    length_right = length_right + length
  enddo
endif

! --- Align to previous block
if (attached) then
  if (n_seg .gt. n_seg_prev) then
    n_tmp = n_seg_prev
  else
    n_tmp = n_seg
  endif
  do i_seg=1,n_tmp
    R_seg(i_seg,1) = R_seg_prev(i_seg)
    Z_seg(i_seg,1) = Z_seg_prev(i_seg)
  enddo
endif




!---------------------------------------------------------------------------------!
!------- Points for polar line along the wall on the left side of corner ---------!
!---------------------------------------------------------------------------------!

! --- Find intermediate points along grid bnd aligned to the wall
if (corner_block(i_ext) .eq. 1) then
  if (n_corner_lim .gt. 1) then
    do i=1,n_corner_lim
      R_corner_lim(i) = R_wall(index_corner_lim(i))
      Z_corner_lim(i) = Z_wall(index_corner_lim(i))
    enddo
    ! --- We add a point at the beg/end if we are far from any wall point (end)
    i_lim_prev = index_corner_lim(n_corner_lim) - wall_direction_corner
    if (i_lim_prev .gt. n_wall) i_lim_prev = 1
    if (i_lim_prev .lt. 1     ) i_lim_prev = n_wall
    diff         = sqrt( (R_wall(index_corner_lim(n_corner_lim))-R_block_points_left(i_ext,n_block_points_left(i_ext)))**2 &
                        +(Z_wall(index_corner_lim(n_corner_lim))-Z_block_points_left(i_ext,n_block_points_left(i_ext)))**2 )
    if (diff .gt. wall_node_proximity_tolerance) then
      diff_min_beg = sqrt( (R_wall(index_corner_lim(n_corner_lim))-R_wall(i_lim_prev))**2 &
                          +(Z_wall(index_corner_lim(n_corner_lim))-Z_wall(i_lim_prev))**2 )
      diff         = sqrt( (R_wall(i_lim_prev)  -R_block_points_left(i_ext,n_block_points_left(i_ext)))**2 &
                          +(Z_wall(i_lim_prev)  -Z_block_points_left(i_ext,n_block_points_left(i_ext)))**2 )
      if (diff .gt. diff_min_beg) n_corner_lim = n_corner_lim + 1
      R_corner_lim(n_corner_lim) = R_block_points_left(i_ext,n_block_points_left(i_ext))
      Z_corner_lim(n_corner_lim) = Z_block_points_left(i_ext,n_block_points_left(i_ext))
    endif
    ! --- We add a point at the beg/end if we are far from any wall point (beg)
    i_lim_prev = index_corner_lim(1) - wall_direction_corner
    if (i_lim_prev .gt. n_wall) i_lim_prev = 1
    if (i_lim_prev .lt. 1     ) i_lim_prev = n_wall
    diff         = sqrt( (R_wall(index_corner_lim(1))-R_block_points_left(i_ext,1))**2 &
                        +(Z_wall(index_corner_lim(1))-Z_block_points_left(i_ext,1))**2 )
    if (diff .gt. wall_node_proximity_tolerance) then
      diff_min_beg = sqrt( (R_wall(index_corner_lim(1))-R_wall(i_lim_prev))**2 &
                          +(Z_wall(index_corner_lim(1))-Z_wall(i_lim_prev))**2 )
      diff         = sqrt( (R_wall(i_lim_prev)  -R_block_points_left(i_ext,1))**2 &
                          +(Z_wall(i_lim_prev)  -Z_block_points_left(i_ext,1))**2 )
      if (diff .lt. diff_min_beg) then
        n_corner_lim = n_corner_lim + 1
        do i = n_corner_lim,2,-1
          R_corner_lim(i) = R_corner_lim(i-1)
          Z_corner_lim(i) = Z_corner_lim(i-1)
        enddo
      endif
      R_corner_lim(1) = R_block_points_left(i_ext,1)
      Z_corner_lim(1) = Z_block_points_left(i_ext,1)
    endif
  else
    n_corner_lim = 2
    R_corner_lim(1) = R_block_points_left(i_ext,1)
    Z_corner_lim(1) = Z_block_points_left(i_ext,1)
    R_corner_lim(2) = R_block_points_left(i_ext,n_block_points_left(i_ext))
    Z_corner_lim(2) = Z_block_points_left(i_ext,n_block_points_left(i_ext))
  endif
  call create_polar_lines_simple(n_corner_lim, R_corner_lim(1:n_corner_lim), Z_corner_lim(1:n_corner_lim), R_polar_wall(1:n_corner_lim-1,1:4) , Z_polar_wall(1:n_corner_lim-1,1:4) )
  length_left = 0.d0
  do i=1,n_corner_lim-1
    call from_polar_to_cubic(R_polar_wall(i,1:4),R_cub1d)
    call from_polar_to_cubic(Z_polar_wall(i,1:4),Z_cub1d)
    call curve_length(R_cub1d(1), R_cub1d(2), R_cub1d(3), R_cub1d(4), &
                      Z_cub1d(1), Z_cub1d(2), Z_cub1d(3), Z_cub1d(4), -1.d0, 1.d0, length)
    length_left = length_left + length
  enddo
endif




!------------------------------------------------------!
!------- Segment points along the wall for corners ----!
!------------------------------------------------------!

! --- Find intermediate points along sides
if (corner_block(i_ext) .eq. 1) then
  R_seg(1,1)       = R_polar_wall(1,1)
  Z_seg(1,1)       = Z_polar_wall(1,1)
  length_seg = 0.d0
  do i_seg = 2,n_seg-1
    length_find = seg(i_seg)
    length_sum  = 0.d0
    do i=1,n_corner_lim-1
      call from_polar_to_cubic(R_polar_wall(i,1:4),R_cub1d)
      call from_polar_to_cubic(Z_polar_wall(i,1:4),Z_cub1d)
      call curve_length(R_cub1d(1), R_cub1d(2), R_cub1d(3), R_cub1d(4), &
                        Z_cub1d(1), Z_cub1d(2), Z_cub1d(3), Z_cub1d(4), -1.d0, 1.d0, length)
      if ((length_sum + length)/length_left .ge. length_find) then
        length_tmp = length_find*length_left - length_sum
        st = 2.d0 * length_tmp/length - 1.d0
        call CUB1D(R_cub1d(1), R_cub1d(2), R_cub1d(3), R_cub1d(4), st, R3, dR3_dr)
        call CUB1D(Z_cub1d(1), Z_cub1d(2), Z_cub1d(3), Z_cub1d(4), st, Z3, dZ3_dr)
        R_seg(i_seg,1) = R3
        Z_seg(i_seg,1) = Z3
        exit
      endif
      length_sum = length_sum + length
    enddo
  enddo
  R_seg(n_seg,1) = R_polar_wall(n_corner_lim-1,4)
  Z_seg(n_seg,1) = Z_polar_wall(n_corner_lim-1,4)
endif






!----------------------------------- Print a python file that plots a cross with the 4 nodes of each element
if (plot_grid) then
  open(101,file='plot_ext_sides.py')
    write(101,'(A)')                '#!/usr/bin/env python'
    write(101,'(A)')                'import numpy as N'
    write(101,'(A)')                'import pylab'
    write(101,'(A)')                'def main():'
    
    do l=1,2
      if (l.eq.1) k = 1
      if (l.eq.2) k = n_nodes
      write(101,'(A,i6,A)')         ' r = N.zeros(',n_seg,')'
      write(101,'(A,i6,A)')         ' z = N.zeros(',n_seg,')'
      do i=1,n_seg
        write(101,'(A,i6,A,f15.4)') ' r[',i-1,'] = ',R_seg(i,k)
        write(101,'(A,i6,A,f15.4)') ' z[',i-1,'] = ',Z_seg(i,k)
      enddo
      write(101,'(A)')              ' pylab.plot(r,z, "b-x")'
      write(101,'(A)')              ' pylab.plot(r[0],z[0], "rx")'
      write(101,'(A,i3)')           ' n_points = ',n_seg
      write(101,'(A)')              ' pylab.plot(r[n_points-1],z[n_points-1], "rx")'
    enddo
    
    write(101,'(A,i6,A)')           ' r = N.zeros(',n_nodes,')'
    write(101,'(A,i6,A)')           ' z = N.zeros(',n_nodes,')'
    do i=1,n_nodes
      write(101,'(A,i6,A,f15.4)')   ' r[',i-1,'] = ',R_seg_bnd(i)
      write(101,'(A,i6,A,f15.4)')   ' z[',i-1,'] = ',Z_seg_bnd(i)
    enddo
    write(101,'(A)')                ' pylab.plot(r,z, "b-x")'
    write(101,'(A)')                ' pylab.plot(r[0],z[0], "rx")'
    write(101,'(A,i3)')             ' n_points = ',n_nodes
    write(101,'(A)')                ' pylab.plot(r[n_points-1],z[n_points-1], "rx")'
    
    write(101,'(A,i6,A)')           ' r = N.zeros(',n_nodes,')'
    write(101,'(A,i6,A)')           ' z = N.zeros(',n_nodes,')'
    do i=1,n_nodes
      write(101,'(A,i6,A,f15.4)')   ' r[',i-1,'] = ',R_seg_wall(i)
      write(101,'(A,i6,A,f15.4)')   ' z[',i-1,'] = ',Z_seg_wall(i)
    enddo
    write(101,'(A)')                ' pylab.plot(r,z, "b-x")'
    write(101,'(A)')                ' pylab.plot(r[0],z[0], "rx")'
    write(101,'(A,i3)')             ' n_points = ',n_nodes
    write(101,'(A)')                ' pylab.plot(r[n_points-1],z[n_points-1], "rx")'
    
    write(101,'(A)')                ' pylab.axis("equal")'
    write(101,'(A)')                ' pylab.show()'
    write(101,'(A)')                ' '
    write(101,'(A)')                'main()'
  close(101)
endif




!------------------------------------------------------------------------------------------------------------------------!
!************************************************************************************************************************!
!************************************************************************************************************************!
!*********************************** Second part: find crossings of lines ***********************************************!
!************************************************************************************************************************!
!************************************************************************************************************************!
!------------------------------------------------------------------------------------------------------------------------!
write(*,*) '                 Find crossings between coordinate lines'



! --- We need to determine the deviation between our sides and a straight line
do i=1,n_nodes
  ! --- Deviation of bnd nodes
  R1 = R_seg_bnd(1) + seg_bnd(i) * (R_seg_bnd(n_nodes)-R_seg_bnd(1))
  Z1 = Z_seg_bnd(1) + seg_bnd(i) * (Z_seg_bnd(n_nodes)-Z_seg_bnd(1))
  R2 = R_seg_bnd(i)
  Z2 = Z_seg_bnd(i)
  R_dev_bnd(i) = R2-R1
  Z_dev_bnd(i) = Z2-Z1
  ! --- Deviation of wall points
  R1 = R_seg_wall(1) + seg_wall(i) * (R_seg_wall(n_nodes)-R_seg_wall(1))
  Z1 = Z_seg_wall(1) + seg_wall(i) * (Z_seg_wall(n_nodes)-Z_seg_wall(1))
  R2 = R_seg_wall(i)
  Z2 = Z_seg_wall(i)
  R_dev_wall(i) = R2-R1
  Z_dev_wall(i) = Z2-Z1
  ! --- Averaged deviation
  do i_seg=1,n_seg
    R_deviation(i_seg,i) = (1.d0-seg(i_seg)) * R_dev_bnd(i) + seg(i_seg) * R_dev_wall(i)
    Z_deviation(i_seg,i) = (1.d0-seg(i_seg)) * Z_dev_bnd(i) + seg(i_seg) * Z_dev_wall(i)
  enddo
enddo

! --- Define the grid points
do i=2,n_nodes-1
  R_seg(1,i) = R_seg_bnd(i)
  Z_seg(1,i) = Z_seg_bnd(i)
  R_seg(n_seg,i) = R_seg_wall(i)
  Z_seg(n_seg,i) = Z_seg_wall(i)
  do i_seg=2,n_seg-1
    ! --- The straight line
    R_seg(i_seg,i) = R_seg(i_seg,1) + seg_bnd(i) * (R_seg(i_seg,n_nodes)-R_seg(i_seg,1))
    Z_seg(i_seg,i) = Z_seg(i_seg,1) + seg_bnd(i) * (Z_seg(i_seg,n_nodes)-Z_seg(i_seg,1))
    ! --- Plus the deviation
    R_seg(i_seg,i) = R_seg(i_seg,i) + R_deviation(i_seg,i)
    Z_seg(i_seg,i) = Z_seg(i_seg,i) + Z_deviation(i_seg,i)
  enddo
enddo


!----------------------------------- Print a python file that plots a cross with the 4 nodes of each element
if (plot_grid) then
  open(101,file='plot_ext_points.py')
    write(101,'(A)')                '#!/usr/bin/env python'
    write(101,'(A)')                'import numpy as N'
    write(101,'(A)')                'import pylab'
    write(101,'(A)')                'def main():'
    
    do i=1,n_nodes
      write(101,'(A,i6,A)')         ' r = N.zeros(',n_seg,')'
      write(101,'(A,i6,A)')         ' z = N.zeros(',n_seg,')'
      do i_seg=1,n_seg
        write(101,'(A,i6,A,f15.4)') ' r[',i_seg-1,'] = ',R_seg(i_seg,i)
        write(101,'(A,i6,A,f15.4)') ' z[',i_seg-1,'] = ',Z_seg(i_seg,i)
      enddo
      write(101,'(A)')              ' pylab.plot(r,z, "b-x")'
    enddo
    
    do i_seg=1,n_seg
      write(101,'(A,i6,A)')         ' r = N.zeros(',n_nodes,')'
      write(101,'(A,i6,A)')         ' z = N.zeros(',n_nodes,')'
      do i=1,n_nodes
        write(101,'(A,i6,A,f15.4)') ' r[',i-1,'] = ',R_seg(i_seg,i)
        write(101,'(A,i6,A,f15.4)') ' z[',i-1,'] = ',Z_seg(i_seg,i)
      enddo
      write(101,'(A)')              ' pylab.plot(r,z, "b-x")'
    enddo
    
    write(101,'(A)')                ' pylab.axis("equal")'
    write(101,'(A)')                ' pylab.show()'
    write(101,'(A)')                ' '
    write(101,'(A)')                'main()'
  close(101)
endif



! --- Create radial polar lines
do i_seg = 1,n_seg
  call create_polar_lines_simple(n_nodes, R_seg(i_seg,1:n_nodes), Z_seg(i_seg,1:n_nodes), R_polar_radial(1:n_nodes-1,1:4,i_seg) , Z_polar_radial(1:n_nodes-1,1:4,i_seg) )
enddo
! --- Create sides polar lines
do i = 1,n_nodes
  call create_polar_lines_simple(n_seg, R_seg(1:n_seg,i), Z_seg(1:n_seg,i), R_polar_sides(1:n_seg-1,1:4,i) , Z_polar_sides(1:n_seg-1,1:4,i) )
enddo



! --- Now we want to get a smooth transition of radial segmentation between the grid and the new extension patch
i_start = 1
if (attached) i_start = 2
i_end = n_nodes
if (attached_side) i_end = n_nodes-1
do i = i_start,i_end
  ! --- First get the element we have at the bnd
  if (i .le. 2) then
    i_elm = i_elm_save
  else
    found_elm = .false.
    i_elm_save = i_elm
    do i_elm=1,element_list%n_elements
      if (i_elm .eq. i_elm_save) cycle
      do i_node = 1,n_vertex_max
        if (element_list%element(i_elm)%vertex(i_node) .eq. index_bnd(i-1)) then
          found_elm = .true.
          i_elm_save = i_elm
          exit
        endif
      enddo
      if (found_elm) exit
    enddo
    if (.not. found_elm) then
      write(*,*) 'Could not find next element on boundary. Aborting...'
      stop
    endif
    i_elm = i_elm_save
  endif
  ! --- First get all the lengths between the grid boundary and the previous surface (at each bnd_node)
  if (i .eq. 1) then
    if (element_direction .eq. 1) index_bnd_tmp = 4
    if (element_direction .eq. 2) index_bnd_tmp = 3
    if (element_direction .eq. 3) index_bnd_tmp = 1
    if (element_direction .eq. 4) index_bnd_tmp = 4
    if (element_direction .eq. 5) index_bnd_tmp = 2
    if (element_direction .eq. 6) index_bnd_tmp = 1
    if (element_direction .eq. 7) index_bnd_tmp = 3
    if (element_direction .eq. 8) index_bnd_tmp = 2
  else
    if (element_direction .eq. 1) index_bnd_tmp = 3
    if (element_direction .eq. 2) index_bnd_tmp = 4
    if (element_direction .eq. 3) index_bnd_tmp = 4
    if (element_direction .eq. 4) index_bnd_tmp = 1
    if (element_direction .eq. 5) index_bnd_tmp = 1
    if (element_direction .eq. 6) index_bnd_tmp = 2
    if (element_direction .eq. 7) index_bnd_tmp = 2
    if (element_direction .eq. 8) index_bnd_tmp = 3
  endif
  i_node = element_list%element(i_elm)%vertex(index_bnd_tmp)
  previous_length = sqrt(  (node_list%node(i_node)%x(1,1,1)-node_list%node(index_bnd(i))%x(1,1,1))**2 &
                         + (node_list%node(i_node)%x(1,1,2)-node_list%node(index_bnd(i))%x(1,1,2))**2 )
  ! --- Correct length with respective angle
  alpha1 = atan2(node_list%node(i_node)%x(1,1,2)-node_list%node(index_bnd(i))%x(1,1,2),&
                 node_list%node(i_node)%x(1,1,1)-node_list%node(index_bnd(i))%x(1,1,1))
  if (alpha1 .lt. 0.d0) alpha1 = alpha1 + 2.d0 * PI
  if (i .eq. 1) then
    alpha2 = atan2(node_list%node(index_bnd(i+1))%x(1,1,2)-node_list%node(index_bnd(i))%x(1,1,2),&
                   node_list%node(index_bnd(i+1))%x(1,1,1)-node_list%node(index_bnd(i))%x(1,1,1))
  else
    alpha2 = atan2(node_list%node(index_bnd(i-1))%x(1,1,2)-node_list%node(index_bnd(i))%x(1,1,2),&
                   node_list%node(index_bnd(i-1))%x(1,1,1)-node_list%node(index_bnd(i))%x(1,1,1))
  endif
  if (alpha2 .lt. 0.d0) alpha2 = alpha2 + 2.d0 * PI
  alpha = alpha2 - alpha1
  if (alpha .lt. 0.d0   ) alpha = alpha + 2.d0 * PI
  if (alpha .gt. 2.d0*PI) alpha = alpha - 2.d0 * PI
  previous_length = previous_length * abs(sin(alpha))
  ! --- Then get size of all polar lines
  polar_length = 0.d0
  do j=1,n_seg-1
    call from_polar_to_cubic(R_polar_sides(j,1:4,i),R_cub1d)
    call from_polar_to_cubic(Z_polar_sides(j,1:4,i),Z_cub1d)
    call curve_length(R_cub1d(1), R_cub1d(2), R_cub1d(3), R_cub1d(4), &
                      Z_cub1d(1), Z_cub1d(2), Z_cub1d(3), Z_cub1d(4), -1.d0, 1.d0, length)
    polar_length = polar_length + length
    if (j .eq. 1) length_save = length
  enddo
  ! --- Get angle to correct length
  if (i .eq. 1) then
    i_node = index_bnd(i+1)
  else
    i_node = index_bnd(i-1)
  endif
  alpha1 = atan2(node_list%node(i_node)%x(1,1,2)-node_list%node(index_bnd(i))%x(1,1,2),&
                 node_list%node(i_node)%x(1,1,1)-node_list%node(index_bnd(i))%x(1,1,1))
  if (alpha1 .lt. 0.d0) alpha1 = alpha1 + 2.d0 * PI
  alpha2 = atan2(Z_seg(2,i)-node_list%node(index_bnd(i))%x(1,1,2),&
                 R_seg(2,i)-node_list%node(index_bnd(i))%x(1,1,1))
  if (alpha2 .lt. 0.d0) alpha2 = alpha2 + 2.d0 * PI
  alpha = alpha2 - alpha1
  if (alpha .lt. 0.d0   ) alpha = alpha + 2.d0 * PI
  if (alpha .gt. 2.d0*PI) alpha = alpha - 2.d0 * PI
  ! --- In case we do not find good parameters, save an equidistant segmentation
  call meshac2(n_seg,seg_new,1.d0,9999.d0,9999.0,9999.d0,1.d0,1.0d0)
  if (.not. n_ext_equidistant(i_ext)) then
    ! --- Then find the right segmentation to have a smooth transition
    ! --- It is important to do this with high refinment, otherwise the segmentation
    ! --- will be jumpy between each iteration of the main loop [i = i_start,n_nodes]
    n_refine = 4
    s_refine_min_next = 0.d0
    s_refine_max_next = 1.d0
    do i_refine = 1,n_refine
      s_refine_min = s_refine_min_next
      s_refine_max = s_refine_max_next
      n_tmp    = 20 ! we try a few of them, and take the closest one
      diff_min = 1.d10
      sig_tmp  = 0.15d0 ! we take a transition 1/3 of the total length
      if (attached) then
        if (n_block_points_left(i_ext) .gt. n_block_points_right(i_ext-1)) then
          sig_tmp = (sig_tmp * length_prev) / polar_length ! special case
        endif
      endif
      do j = 1,n_tmp
        bgf_tmp = s_refine_min + (s_refine_max - s_refine_min) * real(j)/real(n_tmp+1)
        call meshac2(n_seg,seg_tmp,0.d0,9999.d0,sig_tmp,9999.d0,bgf_tmp,1.0d0)
        length_tmp = polar_length * seg_tmp(2) * abs(sin(alpha))
        diff = abs(length_tmp-previous_length)
        if (diff .lt. diff_min) then
          diff_min = diff
          call meshac2(n_seg,seg_new,0.d0,9999.d0,sig_tmp,9999.d0,bgf_tmp,1.0d0)
          i_save = j
        endif
        ! --- Try the other way around as well
        call meshac2(n_seg,seg_tmp,1.d0,9999.d0,sig_tmp,9999.d0,bgf_tmp,1.0d0)
        length_tmp = polar_length * seg_tmp(2) * abs(sin(alpha))
        diff = abs(length_tmp-previous_length)
        if (diff .lt. diff_min) then
          diff_min = diff
          call meshac2(n_seg,seg_new,1.d0,9999.d0,sig_tmp,9999.d0,bgf_tmp,1.0d0)
          i_save = j
        endif
      enddo
      ! --- New interval for next refinment
      j = max(1,i_save-1)
      s_refine_min_next = s_refine_min + (s_refine_max - s_refine_min) * real(j)/real(n_tmp+1)
      j = min(n_tmp,i_save+1)
      s_refine_max_next = s_refine_min + (s_refine_max - s_refine_min) * real(j)/real(n_tmp+1)
    enddo
  endif
  ! --- In case you know you want a specific extension beyond the previous patch
  if (attached) then
    if ( (n_block_points_left(i_ext) .gt. n_block_points_right(i_ext-1)) .and. (n_ext_block(i_ext) .gt. n_seg_prev+2) ) then
      do j=1,n_seg
        seg(j) = seg(j) + real(i-1)/real(n_nodes-1) * (seg_new(j) - seg(j))
      enddo
    else
      seg(1:n_seg) = seg_new(1:n_seg)
    endif
  else
    seg(1:n_seg) = seg_new(1:n_seg)
  endif
  ! --- Now, we need to resegment these polar lines
  do i_seg = 2,n_seg-1
    length_find = seg(i_seg) * polar_length
    length_sum  = 0.d0
    do j=1,n_seg-1
      call from_polar_to_cubic(R_polar_sides(j,1:4,i),R_cub1d)
      call from_polar_to_cubic(Z_polar_sides(j,1:4,i),Z_cub1d)
      call curve_length(R_cub1d(1), R_cub1d(2), R_cub1d(3), R_cub1d(4), &
                        Z_cub1d(1), Z_cub1d(2), Z_cub1d(3), Z_cub1d(4), -1.d0, 1.d0, length)
      if (length_sum + length .ge. length_find) then
        length_tmp = length_find - length_sum
        st = 2.d0 * length_tmp/length - 1.d0
        call CUB1D(R_cub1d(1), R_cub1d(2), R_cub1d(3), R_cub1d(4), st, R3, dR3_dr)
        call CUB1D(Z_cub1d(1), Z_cub1d(2), Z_cub1d(3), Z_cub1d(4), st, Z3, dZ3_dr)
        R_seg(i_seg,i) = R3
        Z_seg(i_seg,i) = Z3
        exit
      endif
      length_sum = length_sum + length
    enddo
  enddo
  ! --- And finally, we respline the new segments
  call create_polar_lines_simple(n_seg, R_seg(1:n_seg,i), Z_seg(1:n_seg,i), R_polar_sides(1:n_seg-1,1:4,i) , Z_polar_sides(1:n_seg-1,1:4,i) )
enddo
! --- And we need to do the same for the other direction as well
do i_seg = 1,n_seg
  call create_polar_lines_simple(n_nodes, R_seg(i_seg,1:n_nodes), Z_seg(i_seg,1:n_nodes), R_polar_radial(1:n_nodes-1,1:4,i_seg) , Z_polar_radial(1:n_nodes-1,1:4,i_seg) )
enddo









!------------------------------------------------------------------------------------------------------------------------!
!************************************************************************************************************************!
!************************************************************************************************************************!
!***************************************** Third part: define the new nodes  ********************************************!
!************************************************************************************************************************!
!************************************************************************************************************************!
!------------------------------------------------------------------------------------------------------------------------!
write(*,*) '                 Defining new nodes'


index = n_start
do i_seg=1,n_seg
  do i=1,n_nodes
    index = index + 1
    ! --- If extending on the s-side
    if ( (element_direction .eq. 3) .or. (element_direction .eq. 4) .or. (element_direction .eq. 7) .or. (element_direction .eq. 8) )  then
      call create_new_node_polar(newnode_list, index, n_seg, n_nodes, i_seg-1, i-1, R_polar_radial, Z_polar_radial, R_polar_sides, Z_polar_sides, R_seg(i_seg,i), Z_seg(i_seg,i))
    ! --- If extending on the t-side
    else
      call create_new_node_polar(newnode_list, index, n_nodes, n_seg, i-1, i_seg-1, R_polar_sides, Z_polar_sides, R_polar_radial, Z_polar_radial, R_seg(i_seg,i), Z_seg(i_seg,i))
    endif
  enddo
enddo
newnode_list%n_nodes = index


!----------------------------------- Print a python file that plots a cross with the 4 nodes of each element
if (plot_grid) then
  open(101,file='plot_ext_nodes.py')
    write(101,'(A)')                '#!/usr/bin/env python'
    write(101,'(A)')                'import numpy as N'
    write(101,'(A)')                'import pylab'
    write(101,'(A)')                'def main():'
    
    write(101,'(A,i6,A)')           ' r = N.zeros(',newnode_list%n_nodes,')'
    write(101,'(A,i6,A)')           ' z = N.zeros(',newnode_list%n_nodes,')'
    do i=1,newnode_list%n_nodes
      write(101,'(A,i6,A,f15.4)') ' r[',i-1,'] = ',newnode_list%node(i)%x(1,1,1)
      write(101,'(A,i6,A,f15.4)') ' z[',i-1,'] = ',newnode_list%node(i)%x(1,1,2)
    enddo
    write(101,'(A)')              ' pylab.plot(r,z, "bx")'

    write(101,'(A)')                ' pylab.axis("equal")'
    write(101,'(A)')                ' pylab.show()'
    write(101,'(A)')                ' '
    write(101,'(A)')                'main()'
  close(101)
endif




!------------------------------------------------------------------------------------------------------------------------!
!************************************************************************************************************************!
!************************************************************************************************************************!
!*************************************** Fourth part: define the new elements  ******************************************!
!************************************************************************************************************************!
!************************************************************************************************************************!
!------------------------------------------------------------------------------------------------------------------------!
write(*,*) '                 Defining new elements'


!-------------------------------- The closed region
index = 0
do i=1,n_seg-1
  do j=1,n_nodes-1

    index = index + 1
    newelement_list%element(index)%size = 1.d0

    if (element_direction .eq. 1) then
      newelement_list%element(index)%vertex(4) = n_start + (i-1)*n_nodes + j
      newelement_list%element(index)%vertex(1) = n_start + (i  )*n_nodes + j
      newelement_list%element(index)%vertex(2) = n_start + (i  )*n_nodes + j + 1
      newelement_list%element(index)%vertex(3) = n_start + (i-1)*n_nodes + j + 1
    elseif (element_direction .eq. 2) then
      newelement_list%element(index)%vertex(3) = n_start + (i-1)*n_nodes + j
      newelement_list%element(index)%vertex(2) = n_start + (i  )*n_nodes + j
      newelement_list%element(index)%vertex(1) = n_start + (i  )*n_nodes + j + 1
      newelement_list%element(index)%vertex(4) = n_start + (i-1)*n_nodes + j + 1
    elseif (element_direction .eq. 3) then
      newelement_list%element(index)%vertex(1) = n_start + (i-1)*n_nodes + j
      newelement_list%element(index)%vertex(2) = n_start + (i  )*n_nodes + j
      newelement_list%element(index)%vertex(3) = n_start + (i  )*n_nodes + j + 1
      newelement_list%element(index)%vertex(4) = n_start + (i-1)*n_nodes + j + 1
    elseif (element_direction .eq. 4) then
      newelement_list%element(index)%vertex(4) = n_start + (i-1)*n_nodes + j
      newelement_list%element(index)%vertex(3) = n_start + (i  )*n_nodes + j
      newelement_list%element(index)%vertex(2) = n_start + (i  )*n_nodes + j + 1
      newelement_list%element(index)%vertex(1) = n_start + (i-1)*n_nodes + j + 1
    elseif (element_direction .eq. 5) then
      newelement_list%element(index)%vertex(2) = n_start + (i-1)*n_nodes + j
      newelement_list%element(index)%vertex(3) = n_start + (i  )*n_nodes + j
      newelement_list%element(index)%vertex(4) = n_start + (i  )*n_nodes + j + 1
      newelement_list%element(index)%vertex(1) = n_start + (i-1)*n_nodes + j + 1
    elseif (element_direction .eq. 6) then
      newelement_list%element(index)%vertex(1) = n_start + (i-1)*n_nodes + j
      newelement_list%element(index)%vertex(4) = n_start + (i  )*n_nodes + j
      newelement_list%element(index)%vertex(3) = n_start + (i  )*n_nodes + j + 1
      newelement_list%element(index)%vertex(2) = n_start + (i-1)*n_nodes + j + 1
    elseif (element_direction .eq. 7) then
      newelement_list%element(index)%vertex(3) = n_start + (i-1)*n_nodes + j
      newelement_list%element(index)%vertex(4) = n_start + (i  )*n_nodes + j
      newelement_list%element(index)%vertex(1) = n_start + (i  )*n_nodes + j + 1
      newelement_list%element(index)%vertex(2) = n_start + (i-1)*n_nodes + j + 1
    elseif (element_direction .eq. 8) then
      newelement_list%element(index)%vertex(2) = n_start + (i-1)*n_nodes + j
      newelement_list%element(index)%vertex(1) = n_start + (i  )*n_nodes + j
      newelement_list%element(index)%vertex(4) = n_start + (i  )*n_nodes + j + 1
      newelement_list%element(index)%vertex(3) = n_start + (i-1)*n_nodes + j + 1
    endif
      
  enddo
enddo
newelement_list%n_elements = index



!----------------------------------- Print a python file that plots a cross with the 4 nodes of each element
if (plot_grid) then
  n_tmp = newelement_list%n_elements
  if (i_ext .ge. 10) then
    write(char_tmp2,'(i2)')i_ext
    plot_filename = 'plot_extension_'//char_tmp2//'.py'
  else
    write(char_tmp,'(i1)')i_ext
    plot_filename = 'plot_extension_'//char_tmp//'.py'
  endif
  open(101,file=plot_filename)
    write(101,'(A)')         '#!/usr/bin/env python'
    write(101,'(A)')         'import numpy as N'
    write(101,'(A)')         'import pylab'
    write(101,'(A)')         'def main():'
    write(101,'(A,i6,A)')    ' r = N.zeros(',4*n_tmp,')'
    write(101,'(A,i6,A)')    ' z = N.zeros(',4*n_tmp,')'
    do j=1,n_tmp
      do i=1,2
        index = newelement_list%element(j)%vertex(i)
        write(101,'(A,i6,A,f15.4)') ' r[',4*(j-1)+2*i-2,'] = ',newnode_list%node(index)%x(1,1,1)
        write(101,'(A,i6,A,f15.4)') ' z[',4*(j-1)+2*i-2,'] = ',newnode_list%node(index)%x(1,1,2)
        index = newelement_list%element(j)%vertex(i+2)
        write(101,'(A,i6,A,f15.4)') ' r[',4*(j-1)+2*i-1,'] = ',newnode_list%node(index)%x(1,1,1)
        write(101,'(A,i6,A,f15.4)') ' z[',4*(j-1)+2*i-1,'] = ',newnode_list%node(index)%x(1,1,2)
      enddo
    enddo
    write(101,'(A,i6,A)')    ' for i in range (0,',n_tmp,'):'
    write(101,'(A)')         '  pylab.plot(r[4*i:4*i+2],z[4*i:4*i+2], "r")'
    write(101,'(A)')         '  pylab.plot(r[4*i+2:4*i+4],z[4*i+2:4*i+4], "g")'
    do j=1,n_tmp
      do i=1,4
        index = newelement_list%element(j)%vertex(i)
        write(101,'(A,i6,A,f15.4)') ' r[',4*(j-1)+i-1,'] = ',newnode_list%node(index)%x(1,1,1)
        write(101,'(A,i6,A,f15.4)') ' z[',4*(j-1)+i-1,'] = ',newnode_list%node(index)%x(1,1,2)
      enddo
    enddo
    write(101,'(A,i6,A)')    ' for i in range (0,',n_tmp,'):'
    write(101,'(A)')         '  pylab.plot(r[4*i:4*i+4],z[4*i:4*i+4], "b")'
    write(101,'(A)')         ' pylab.axis("equal")'
    write(101,'(A)')         ' pylab.show()'
    write(101,'(A)')         ' '
    write(101,'(A)')         'main()'
  close(101)
endif


! --- Make sure we send the segmentation to the next block!
do i_seg=1,n_seg
  seg_prev(i_seg) = seg(i_seg)
enddo
n_seg_prev = n_seg

deallocate(seg,       R_seg,       Z_seg, seg_tmp, seg_new)
deallocate(seg_bnd,   R_seg_bnd,   Z_seg_bnd)
deallocate(seg_wall,  R_seg_wall,  Z_seg_wall)
deallocate(R_dev_bnd,   Z_dev_bnd)
deallocate(R_dev_wall,  Z_dev_wall)
deallocate(R_deviation, Z_deviation)
deallocate(R_polar_radial,Z_polar_radial)
deallocate(R_polar_sides ,Z_polar_sides )
if (attached) deallocate(R_seg_prev, Z_seg_prev)

return
end subroutine define_extension_patch













subroutine find_next_bnd_node(node_list,element_list,i_node_in,direction,i_node_next)

  use data_structure
  
  implicit none
  
  ! --- Routine variables
  type (type_node_list),    intent(in)   :: node_list
  type (type_element_list), intent(in)   :: element_list
  integer,                  intent(in)   :: i_node_in
  integer,                  intent(in)   :: direction
  integer,                  intent(out)  :: i_node_next
  
  ! --- Internal variables
  integer :: i_elm, i_vertex, i_node, i_vertex_next

  ! --- Find next bnd node along boundary
  i_node_next = 0
  do i_elm = 1,element_list%n_elements
    do i_vertex=1,4
      i_node = element_list%element(i_elm)%vertex(i_vertex)
      if (node_list%node(i_node)%boundary .eq. 0) cycle
      if (i_node .eq. i_node_in) then
        if (direction .eq. +1) then
          i_vertex_next = mod(i_vertex,4) + 1
        else
          i_vertex_next = mod(i_vertex+2,4) + 1
        endif
        i_node = element_list%element(i_elm)%vertex(i_vertex_next)
        if (node_list%node(i_node)%boundary .gt. 0) then
          i_node_next = i_node
          return
        endif
      endif
    enddo
  enddo

  return
end subroutine find_next_bnd_node























