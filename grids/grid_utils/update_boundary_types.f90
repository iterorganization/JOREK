subroutine update_boundary_types(element_list,node_list, across_xpoint)
  
  ! -------------------------------------------------------------------
  ! -------------------------------------------------------------------
  ! --------- IMPORTANT NOTE! THIS ROUTINE ASSUMES THAT NODES ---------
  ! --------- 1-4 ARE XPOINT NODES (AND ALSO 5-8 IF XCASE=3)  ---------
  ! --------- IF YOU ARE BUILDING A PATCH THAT DOES NOT HAVE  ---------
  ! --------- ANY XPOINT IN IT, THEN MAKE SURE THE ELEMENT    ---------
  ! --------- VERTEX INDEXING STARTS AT N_XPOINT+1            ---------
  ! -------------------------------------------------------------------
  ! -------------------------------------------------------------------
  
  use constants
  use mod_parameters
  use data_structure
  use phys_module, only: xcase
  
  implicit none
  
  ! --- Routine variables
  type(type_node_list),    intent(inout) :: node_list      !< list of grid nodes
  type(type_element_list), intent(inout) :: element_list   !< list of finite elements
  integer,                 intent(in)    :: across_xpoint  !< 0=no xpoint at all, 1=we have 4 xpoint nodes, 2=only two
                                                           ! eg, you go from 1 to 4, or you go from 1 to 2 (see create_x_node.f90)
  
  ! --- Internal variables
  integer :: i_elm,       i_vertex,  i_node,  i_node_prev
  integer :: i_elm2,      i_vertex2, i_node2
  integer :: elm_sum
  integer :: i_elm_first, i_vertex_first
  integer :: i_elm_now,   i_vertex_now, i_vertex_next, i_vertex_save
  integer :: i_elm_prev
  integer :: iter, n_xpoints
  logical :: found_first, found_next
  logical, parameter :: debug  = .false.
  logical, parameter :: debug2 = .false.
  integer, parameter :: n_times = 3 ! assume maximum of 2 holes...
  integer            :: i_times
  
  ! --- Some printouts?
  if (debug) then
    write(*,*)'----------------------------------------------------------'
    write(*,*)'---------------- Updating boundary types. ----------------'
    write(*,*)'----------------------------------------------------------'
  endif
  
  if (across_xpoint .gt. 0) then
    if (xcase .eq. DOUBLE_NULL) then
      n_xpoints = 8
    else
      n_xpoints = 4
    endif
  else
    n_xpoints = 0
  endif
  
  ! --- Initialise first!
  do i_node=1,node_list%n_nodes
    node_list%node(i_node)%boundary = 0
  enddo
  
  ! --- There migt be holes, we need to check several times to make sure we don't miss any
  do i_times =1,n_times
    ! --- Boundary check
    ! --- Find first element with a boundary
    !write(*,*)'Looking for first boundary element...'
    i_elm_first = 0
    found_first = .false.
    do i_elm=1,element_list%n_elements
      do i_vertex=1,4
        i_node = element_list%element(i_elm)%vertex(i_vertex)
        if ( (i_node .le. 4) .and. (xcase .ne. DOUBLE_NULL) .and. (across_xpoint .gt. 0) ) cycle
        if ( (i_node .le. 8) .and. (xcase .eq. DOUBLE_NULL) .and. (across_xpoint .gt. 0) ) cycle
        call adjacent_elements(element_list,node_list,i_elm,i_vertex,3,elm_sum)
        ! --- We want a corner node to start with (note this also make it safer if we have axis nodes on our grid)
        ! --- We want an inverted corner when looking for additional holes
        if (     ( (elm_sum .eq. 0) .and. (i_times .eq. 1) ) &
            .or. ( (elm_sum .eq. 2) .and. (i_times .gt. 1) ) ) then
          ! --- Use one of the adjacent nodes to start from
          i_vertex_next = mod(i_vertex+2,4) + 1
          call adjacent_elements(element_list,node_list,i_elm,i_vertex_next,3,elm_sum)
          if (elm_sum .eq. 1) then
            i_vertex_first = i_vertex_next
            i_elm_first    = i_elm
            found_first = .true.
            exit
          else
            ! --- Try the other direction (this might be the xpoint when we have an individual leg alone)
            i_vertex_next = mod(i_vertex,4) + 1
            call adjacent_elements(element_list,node_list,i_elm,i_vertex_next,3,elm_sum)
            if (elm_sum .eq. 1) then
              i_vertex_first = i_vertex_next
              i_elm_first    = i_elm
              found_first = .true.
              exit
            else
              if (debug) write(*,*)'Something strange when looking for the first element'
              if (debug) write(*,*)'Found a type 3 but its neighbours are not bnd...'
            endif  
          endif  
        endif
      enddo
      if (found_first) exit
    enddo
    if (.not. found_first) then
      if (i_times .eq. 1) then
        write(*,*)'Could not find first boundary element. Aborting...'
        return
      else
        cycle
      endif
    endif
    if (debug) write(*,'(A,i6,2f10.3)')'Found first bnd elm    :',i_node,node_list%node(i_node)%x(1,1,1:2)
    
    ! --- Make sure our boundary is coherent (should not last longer than the number of elements)
    found_first  = .false.
    i_elm_now    = i_elm_first
    i_vertex_now = i_vertex_first
    i_elm_prev   = 0
    if (debug) write(*,*)'Going around boundary...'
    do iter=1,element_list%n_elements
      !write(*,*)'iter...',iter
      
      ! --- What type is this boundary?
      call adjacent_elements(element_list,node_list,i_elm_now,i_vertex_now,2,elm_sum)
      if ( (elm_sum .eq. 0) .and. (i_elm_prev .eq. 0) ) then
        write(*,*)'Impossible, we should always start on a non-corner node. Aborting...'
        exit
      endif
      i_node = element_list%element(i_elm_now)%vertex(i_vertex_now)
      if (debug2) write(*,'(A,i6,2f10.3)')'Starting on new element:',i_elm_now
      if (debug2) write(*,'(A,i6,2f10.3)')'On new elm, starting at:',i_node,node_list%node(i_node)%x(1,1,1:2)
      if (elm_sum .eq. 1) node_list%node(i_node)%boundary = 1
      if (elm_sum .eq. 2) node_list%node(i_node)%boundary = 3
      
      ! --- Find the next boundary node on this element
      i_vertex_next = mod(i_vertex_now,4) + 1 ! it should always be the same direction (ie. +1, not -1).
      call adjacent_elements(element_list,node_list,i_elm_now,i_vertex_next,3,elm_sum)
      if (debug2) write(*,'(A,i6,2f10.3,i2)')'Trying next bnd node   :',element_list%element(i_elm_now)%vertex(i_vertex_next),&
                                              node_list%node(element_list%element(i_elm_now)%vertex(i_vertex_next))%x(1,1,1:2),elm_sum
      if (elm_sum .eq. 3) then
        ! --- Something wrong: we went back inside grid. change direction...'
        i_vertex_next = mod(i_vertex_now+2,4) + 1
        call adjacent_elements(element_list,node_list,i_elm_now,i_vertex_next,3,elm_sum)
        if (debug2) write(*,'(A,i6,2f10.3,i2)')'Wrong direction, other :',element_list%element(i_elm_now)%vertex(i_vertex_next),&
                                                node_list%node(element_list%element(i_elm_now)%vertex(i_vertex_next))%x(1,1,1:2),elm_sum
        if (elm_sum .eq. 3) then
          write(*,*)'Something wrong: we went back inside grid. Aborting...'
          exit
        endif
      endif
      i_node_prev = i_node
      i_node = element_list%element(i_elm_now)%vertex(i_vertex_next)
      if (elm_sum .eq. 0) node_list%node(i_node)%boundary = 3
      if (elm_sum .eq. 1) node_list%node(i_node)%boundary = 1
      if (elm_sum .eq. 2) node_list%node(i_node)%boundary = 3
      ! --- On a corner, we need to repeat (careful, xpoints also have only one elm)
      if ( (elm_sum .eq. 0) .and. (i_node .gt. n_xpoints) ) then
        i_vertex_save = i_vertex_next
        i_vertex_next = mod(i_vertex_next,4) + 1 ! it should always be the same direction (ie. +1, not -1)
        call adjacent_elements(element_list,node_list,i_elm_now,i_vertex_next,3,elm_sum)
        if (debug2) write(*,'(A,i6,2f10.3,i2)')'On corner, next node is:',element_list%element(i_elm_now)%vertex(i_vertex_next),&
                                                node_list%node(element_list%element(i_elm_now)%vertex(i_vertex_next))%x(1,1,1:2),elm_sum
        if ( (elm_sum .eq. 3) .or. (element_list%element(i_elm_now)%vertex(i_vertex_next) .eq. i_node_prev) ) then
          ! --- Something wrong: we went back inside grid. change direction...'
          i_vertex_next = mod(i_vertex_save+2,4) + 1
          call adjacent_elements(element_list,node_list,i_elm_now,i_vertex_next,3,elm_sum)
          if (debug2) write(*,'(A,i6,2f10.3,i2)')'Wrong direction, other :',element_list%element(i_elm_now)%vertex(i_vertex_next),&
                                                  node_list%node(element_list%element(i_elm_now)%vertex(i_vertex_next))%x(1,1,1:2),elm_sum
          if (elm_sum .eq. 3) then
            write(*,*)'Something wrong at corner: we went back inside grid. Aborting...'
            exit
          endif
        endif
        i_node = element_list%element(i_elm_now)%vertex(i_vertex_next)
        if (elm_sum .eq. 1) node_list%node(i_node)%boundary = 1
        if (elm_sum .eq. 2) node_list%node(i_node)%boundary = 3
        if (elm_sum .eq. 0) then
          write(*,*)'Impossible, you cannot have two corner nodes on the same element. Aborting...'
          exit
        endif
      endif
      if (debug2) write(*,'(A,i6,2f10.3)')'Found next bnd node    :',i_node,node_list%node(i_node)%x(1,1,1:2)
      i_vertex_now = i_vertex_next
    
      ! --- Find the next element (it should have at least 2 boundary nodes!)
      found_next = .false.
      i_node = element_list%element(i_elm_now)%vertex(i_vertex_now)
      if (i_node .le. n_xpoints) then
        if (across_xpoint .eq. 1) then
          if (mod(i_node,4) .eq. 0) then
            i_node = i_node - 3
          elseif (mod(i_node,4) .eq. 1) then
            i_node = i_node + 3
          else
            if (debug) write(*,*)'Problem: found wrong Xpoint node on our path'
          endif
        elseif (across_xpoint .eq. 2) then
          if (mod(i_node,2) .eq. 0) then
            i_node = i_node - 1
          else
            i_node = i_node + 1
          endif
        endif
      endif
      do i_elm2=1,element_list%n_elements
        if (i_elm2 .eq. i_elm_now) cycle
        if (i_elm2 .eq. i_elm_prev) cycle
        do i_vertex2=1,4
          i_node2 = element_list%element(i_elm2)%vertex(i_vertex2)
          if (i_node2 .eq. i_node) then
            ! --- If this is one, check that the next node is also a boundary
            i_vertex_next = mod(i_vertex2,4) + 1 ! it should always be the same direction (ie. +1, not -1)
            call adjacent_elements(element_list,node_list,i_elm2,i_vertex_next,3,elm_sum)
            if (elm_sum .le. 2) then
              found_next   = .true.
              i_elm_prev   = i_elm_now
              i_elm_now    = i_elm2
              i_vertex_now = i_vertex2
            else
              ! --- Something wrong: change direction...'
              i_vertex_next = mod(i_vertex2+2,4) + 1
              call adjacent_elements(element_list,node_list,i_elm2,i_vertex_next,3,elm_sum)
              if (elm_sum .le. 2) then
                found_next   = .true.
                i_elm_prev   = i_elm_now
                i_elm_now    = i_elm2
                i_vertex_now = i_vertex2
              endif
            endif
            exit
          endif
        enddo
        if (found_next) exit
      enddo
      if (.not. found_next) then
        write(*,*)'Could not find next element. Aborting...',i_node
        exit
      else
        if (debug2) write(*,'(A,i6,2f10.3)')'Found next bnd elm/node:',i_node2,node_list%node(i_node2)%x(1,1,1:2)
      endif
      
      ! --- Check if we looped all around
      if (i_elm_now .eq. i_elm_first) then
        if (debug) write(*,*)'Finished boundary'
        if (i_vertex_now .ne. i_vertex_first) write(*,*)'But finished on wrong node!!!'
        exit
      endif
      
    enddo
  enddo ! n_times
  
  ! --- Now, we need to check non-corner boundaries to see if they are on the s-side or the t-side
  do i_elm=1,element_list%n_elements
    do i_vertex=1,4
      i_node = element_list%element(i_elm)%vertex(i_vertex)
      if (node_list%node(i_node)%boundary .eq. 1) then
        i_vertex2 = mod(i_vertex,4) + 1
        i_node2 = element_list%element(i_elm)%vertex(i_vertex2)
        if (node_list%node(i_node2)%boundary .gt. 0) then
          if (i_vertex+i_vertex2 .eq. 5) node_list%node(i_node)%boundary = 2
          cycle
        endif
        i_vertex2 = mod(i_vertex+2,4) + 1
        i_node2 = element_list%element(i_elm)%vertex(i_vertex2)
        if (node_list%node(i_node2)%boundary .gt. 0) then
          if (i_vertex+i_vertex2 .eq. 5) node_list%node(i_node)%boundary = 2
          cycle
        endif
      endif
    enddo
  enddo
  
  if (debug) then
    open(101,file='plot_bnd_points.py')
      write(101,'(A)')                '#!/usr/bin/env python'
      write(101,'(A)')                'import numpy as N'
      write(101,'(A)')                'import pylab'
      write(101,'(A)')                'def main():'
      do i_node=1,node_list%n_nodes
        write(101,'(A,f15.4)')        ' r = ',node_list%node(i_node)%x(1,1,1)
        write(101,'(A,f15.4)')        ' z = ',node_list%node(i_node)%x(1,1,2)
        if (node_list%node(i_node)%boundary .eq. 1) then
          write(101,'(A)')            ' pylab.plot(r,z, "rx")'
        elseif (node_list%node(i_node)%boundary .eq. 3) then
          write(101,'(A)')            ' pylab.plot(r,z, "gx")'
        elseif (node_list%node(i_node)%boundary .eq. 2) then
          write(101,'(A)')            ' pylab.plot(r,z, "yx")'
        else
          write(101,'(A)')            ' pylab.plot(r,z, "bx")'
        endif
      enddo
      write(101,'(A)')                ' pylab.axis("equal")'
      write(101,'(A)')                ' pylab.show()'
      write(101,'(A)')                ' '
      write(101,'(A)')                'main()'
    close(101)
  endif
  
  return
  
end subroutine update_boundary_types



subroutine adjacent_elements(element_list,node_list,i_elm,i_vertex,side_max, elm_sum)

  use mod_parameters
  use data_structure
  
  implicit none
  
  ! --- Routine variables
  type(type_node_list),    intent(in)  :: node_list        !< list of grid nodes
  type(type_element_list), intent(in)  :: element_list   !< list of finite elements
  integer,                 intent(in)  :: i_elm,i_vertex,side_max
  integer,                 intent(out) :: elm_sum
  
  ! --- Internal variables
  integer :: i_elm2, i_vertex2, i_node, i_node2
  
  elm_sum = 0
  i_node = element_list%element(i_elm)%vertex(i_vertex)
  do i_elm2=1,element_list%n_elements
    if (i_elm2 .eq. i_elm) cycle
    do i_vertex2=1,4
      i_node2 = element_list%element(i_elm2)%vertex(i_vertex2)
      if (i_node2 .eq. i_node) then
        elm_sum = elm_sum + 1
        exit
      endif
    enddo
    if (elm_sum .eq. side_max) exit
  enddo
  
  return
end subroutine adjacent_elements






subroutine update_boundary_types_final(element_list,node_list)
  
  
  use constants
  use mod_parameters
  use data_structure
  use phys_module, only: xcase, use_simple_bnd_types
  use mod_boundary
  
  implicit none
  
  ! --- Routine variables
  type(type_node_list),    intent(inout) :: node_list      !< list of grid nodes
  type(type_element_list), intent(inout) :: element_list   !< list of finite elements
  
  type (type_bnd_node_list)    :: bnd_node_list
  type (type_bnd_element_list) :: bnd_elm_list
  
  ! --- Internal variables
  integer :: i_elm,       i_vertex,  i_node,  i_node_prev, i
  integer :: i_elm2,      i_vertex2, i_node2, i_node_inside, i_node_side, i_node_side2
  integer :: elm_sum
  integer :: i_elm_first, i_vertex_first
  integer :: i_elm_now,   i_vertex_now, i_vertex_next, i_vertex_save
  integer :: i_elm_prev
  integer :: iter, n_xpoints, count
  integer :: i_elm_bnd(4), i_vertex_bnd(4)
  logical :: found_first, found_next, reverse
  logical, parameter :: debug  = .false.
  real*8  :: xjac
  real*8  :: psi_s, psi_t, psi_x, psi_y
  real*8  :: R, R_s, R_t
  real*8  :: Z, Z_s, Z_t
  real*8  :: BR, BZ
  real*8  :: norm_R,  norm_Z
  real*8  :: tang_R,  tang_Z
  real*8  :: tang_R2, tang_Z2
  real*8  :: alpha_Bp, alpha_norm, alpha_tang, alpha_tang2, alpha_tmp, alpha_between1, alpha_between2
  logical :: surface_is_tangent
  real*8, parameter  :: tol_tangent = 3.d-2 !1.5d-2 !1.d-2
  
  ! --- Some printouts?
  if (debug) then
    write(*,*)'----------------------------------------------------------'
    write(*,*)'---------------- Updating boundary types. ----------------'
    write(*,*)'----------------------------------------------------------'
  endif
  
  ! --- Loop over nodes
  do i_node=1,node_list%n_nodes
    if (node_list%node(i_node)%boundary .eq. 0) cycle
    ! --- Count how many elements we have on this node
    count = 0
    do i_elm=1,element_list%n_elements
      do i_vertex=1,4
        i_node2 = element_list%element(i_elm)%vertex(i_vertex)
        if (i_node2 .eq. i_node) then
          count = count + 1
          i_elm_bnd(count)  = i_elm
          i_vertex_bnd(count) = i_vertex
        endif
      enddo
    enddo

    ! --- psi and RZ variables
    R         = node_list%node(i_node)%x(1,1,1)
    R_s       = node_list%node(i_node)%x(1,2,1)
    R_t       = node_list%node(i_node)%x(1,3,1)
    Z         = node_list%node(i_node)%x(1,1,2)
    Z_s       = node_list%node(i_node)%x(1,2,2)
    Z_t       = node_list%node(i_node)%x(1,3,2)
    xjac      =  R_s*Z_t - R_t*Z_s
    psi_s     = node_list%node(i_node)%values(1,2,1)
    psi_t     = node_list%node(i_node)%values(1,3,1)
    psi_x     = (   Z_t * psi_s - Z_s * psi_t ) / xjac
    psi_y     = ( - R_t * psi_s + R_s * psi_t ) / xjac
    
    ! --- Poloidal field
    BR =  psi_y / R
    BZ = -psi_x / R
    alpha_Bp = atan2(BZ,BR)
    if (alpha_Bp .lt. 0.d0) alpha_Bp = alpha_Bp + 2.d0*PI
           

    ! --- If we have an outer corner
    if (count .eq. 1) then
      
      ! --- The insider
      i_node_inside = mod(i_vertex_bnd(1),4) + 1
      i_node_inside = element_list%element(i_elm_bnd(1))%vertex(i_node_inside)
      i_node_side   = mod(i_vertex_bnd(1)+2,4) + 1
      i_node_side   = element_list%element(i_elm_bnd(1))%vertex(i_node_side)
      
      ! --- The tangent vector
      tang_R = node_list%node(i_node_inside)%x(1,1,1) - node_list%node(i_node)%x(1,1,1)
      tang_Z = node_list%node(i_node_inside)%x(1,1,2) - node_list%node(i_node)%x(1,1,2)
      alpha_tang = atan2(tang_Z,tang_R)
      if (alpha_tang .lt. 0.d0) alpha_tang = alpha_tang + 2.d0*PI

      ! --- The other tangent vector
      tang_R = node_list%node(i_node_side)%x(1,1,1) - node_list%node(i_node)%x(1,1,1)
      tang_Z = node_list%node(i_node_side)%x(1,1,2) - node_list%node(i_node)%x(1,1,2)
      alpha_tang2 = atan2(tang_Z,tang_R)
      if (alpha_tang2 .lt. 0.d0) alpha_tang2 = alpha_tang2 + 2.d0*PI
      
      ! --- Swap angles to have good direction
      if (alpha_tang2 .lt. alpha_tang) then
        alpha_norm  = alpha_tang2
        alpha_tang2 = alpha_tang
        alpha_tang  = alpha_norm
      endif
      if (alpha_tang2-alpha_tang .gt. PI) then
        alpha_norm  = alpha_tang2
        alpha_tang2 = alpha_tang
        alpha_tang  = alpha_norm
        alpha_tang2 = alpha_tang2 + 2.d0*PI
      endif
      
      ! --- Now, are we inside or outside
      if (alpha_Bp .lt. alpha_tang) alpha_Bp = alpha_Bp + 2.d0*PI
      if ( (alpha_Bp .gt. alpha_tang) .and. (alpha_Bp .lt. alpha_tang2) ) then
        node_list%node(i_node)%boundary = 19
      else
        node_list%node(i_node)%boundary = 9
      endif
      
    ! --- If we have an inner corner (there we need to let is free, or fix it???)
    elseif (count .eq. 3) then
      
      node_list%node(i_node)%boundary = 21
    
    ! --- Standard boundary
    else
      ! --- The insider and tangent nodes
      i_node_inside = mod(i_vertex_bnd(1),4) + 1
      i_node_inside = element_list%element(i_elm_bnd(1))%vertex(i_node_inside)
      i_node_side   = mod(i_vertex_bnd(1)+2,4) + 1
      i_node_side   = element_list%element(i_elm_bnd(1))%vertex(i_node_side)
      if (node_list%node(i_node_inside)%boundary .ne. 0) then
        i_node_inside = mod(i_vertex_bnd(1)+2,4) + 1
        i_node_inside = element_list%element(i_elm_bnd(1))%vertex(i_node_inside)
        i_node_side   = mod(i_vertex_bnd(1),4) + 1
        i_node_side   = element_list%element(i_elm_bnd(1))%vertex(i_node_side)
      endif
      ! --- The other tangent nodes
      i_node_side2 = mod(i_vertex_bnd(2)+2,4) + 1
      i_node_side2 = element_list%element(i_elm_bnd(2))%vertex(i_node_side2)
      if (node_list%node(i_node_side2)%boundary .eq. 0) then
        i_node_side2 = mod(i_vertex_bnd(2),4) + 1
        i_node_side2 = element_list%element(i_elm_bnd(2))%vertex(i_node_side2)
      endif
      
      ! --- The normal vector
      norm_R = node_list%node(i_node_inside)%x(1,1,1) - node_list%node(i_node)%x(1,1,1)
      norm_Z = node_list%node(i_node_inside)%x(1,1,2) - node_list%node(i_node)%x(1,1,2)

      ! --- The tangent vector
      tang_R = node_list%node(i_node_side)%x(1,1,1) - node_list%node(i_node)%x(1,1,1)
      tang_Z = node_list%node(i_node_side)%x(1,1,2) - node_list%node(i_node)%x(1,1,2)

      ! --- The other tangent vector
      tang_R2 = node_list%node(i_node_side2)%x(1,1,1) - node_list%node(i_node)%x(1,1,1)
      tang_Z2 = node_list%node(i_node_side2)%x(1,1,2) - node_list%node(i_node)%x(1,1,2)

      ! --- The angles
      alpha_norm = atan2(norm_Z,norm_R)
      if (alpha_norm .lt. 0.d0) alpha_norm = alpha_norm + 2.d0*PI
      alpha_tang = atan2(tang_Z,tang_R)
      if (alpha_tang .lt. 0.d0) alpha_tang = alpha_tang + 2.d0*PI
      alpha_tang2 = atan2(tang_Z2,tang_R2)
      if (alpha_tang2 .lt. 0.d0) alpha_tang2 = alpha_tang2 + 2.d0*PI
      if (alpha_tang2 .lt. alpha_tang) then
        alpha_tmp   = alpha_tang
        alpha_tang  = alpha_tang2
        alpha_tang2 = alpha_tmp
      endif
      
      if (alpha_norm .lt. alpha_tang) alpha_norm = alpha_norm + 2.d0*PI
      if (alpha_Bp   .lt. alpha_tang) alpha_Bp   = alpha_Bp   + 2.d0*PI
      
      ! Determine if surface is tangent to boundary
      surface_is_tangent = .false.
      if (abs(alpha_Bp-alpha_tang      ) .lt. tol_tangent*2*PI) surface_is_tangent = .true.
      if (abs(alpha_Bp-alpha_tang -2*PI) .lt. tol_tangent*2*PI) surface_is_tangent = .true.
      if (abs(alpha_Bp-alpha_tang2     ) .lt. tol_tangent*2*PI) surface_is_tangent = .true.
      if (abs(alpha_Bp-alpha_tang2-2*PI) .lt. tol_tangent*2*PI) surface_is_tangent = .true.
      alpha_between1 = alpha_tang2 - PI
      alpha_between2 = alpha_tang
      if (alpha_between2 .lt. alpha_between1) then
        alpha_tmp      = alpha_between1
        alpha_between1 = alpha_between2
        alpha_between2 = alpha_tmp
      endif
      if ( (alpha_between1 .le. alpha_Bp     ) .and. (alpha_Bp      .le. alpha_between2) ) surface_is_tangent = .true.
      if ( (alpha_between1 .le. alpha_Bp-1*PI) .and. (alpha_Bp-1*PI .le. alpha_between2) ) surface_is_tangent = .true.
      if ( (alpha_between1 .le. alpha_Bp-2*PI) .and. (alpha_Bp-2*PI .le. alpha_between2) ) surface_is_tangent = .true.
      
      ! --- If the angle is very small, we want to fix conditions (like along a flux surface)
      if (surface_is_tangent) then
        if (node_list%node(i_node)%boundary .eq. 1) node_list%node(i_node)%boundary = 12
        ! note the bnd 2 is already this
      ! --- If the angle is large, we need Bohm
      else
        reverse = .false.
        if ((alpha_norm .lt. alpha_tang2) .and. (alpha_Bp .lt. alpha_tang2)) reverse = .true.
        if ((alpha_norm .gt. alpha_tang2) .and. (alpha_Bp .gt. alpha_tang2)) reverse = .true.
        
        if (reverse) then
          if (node_list%node(i_node)%boundary .eq. 1) then
            node_list%node(i_node)%boundary = 11
          else
            node_list%node(i_node)%boundary = 15
          endif
        else
          if (node_list%node(i_node)%boundary .ne. 1) then
            node_list%node(i_node)%boundary = 5
          endif
        endif
        
      endif
      
    endif
  enddo
  
  
  ! --- Loop over nodes again to check corners that need to be fixed
  ! --- ie. cannot have Vpar=+Cs on one side of corner, and Vpar=-Cs on the other...
  do i_node=1,node_list%n_nodes
    if ( (node_list%node(i_node)%boundary .ne. 9) .and. (node_list%node(i_node)%boundary .ne. 19) ) cycle
    
    ! --- Count how many elements we have on this node
    count = 0
    do i_elm=1,element_list%n_elements
      do i_vertex=1,4
        i_node2 = element_list%element(i_elm)%vertex(i_vertex)
        if (i_node2 .eq. i_node) then
          count = count + 1
          i_elm_bnd(count)  = i_elm
          i_vertex_bnd(count) = i_vertex
        endif
      enddo
    enddo
    
    
    ! --- The side nodes
    i_node_inside = mod(i_vertex_bnd(1),4) + 1
    i_node_inside = element_list%element(i_elm_bnd(1))%vertex(i_node_inside)
    i_node_side   = mod(i_vertex_bnd(1)+2,4) + 1
    i_node_side   = element_list%element(i_elm_bnd(1))%vertex(i_node_side)
    
    ! --- Check consistency
    if (     (node_list%node(i_node_inside)%boundary .eq.  9) &
        .or. (node_list%node(i_node_inside)%boundary .eq. 19) &
        .or. (node_list%node(i_node_inside)%boundary .eq. 21) &
        ) then
      node_list%node(i_node)%boundary = 20
    endif
    if (     (node_list%node(i_node_side)%boundary .eq.  9) &
        .or. (node_list%node(i_node_side)%boundary .eq. 19) &
        .or. (node_list%node(i_node_side)%boundary .eq. 21) &
        ) then
      node_list%node(i_node)%boundary = 20
    endif
    if (     (node_list%node(i_node_inside)%boundary .eq.  1) &
        .or. (node_list%node(i_node_inside)%boundary .eq.  5) &
        ) then
      if (     (node_list%node(i_node_side)%boundary .eq. 11) &
          .or. (node_list%node(i_node_side)%boundary .eq. 15) &
          ) then
        node_list%node(i_node)%boundary = 20
      endif
    endif
    if (     (node_list%node(i_node_inside)%boundary .eq. 11) &
        .or. (node_list%node(i_node_inside)%boundary .eq. 15) &
        ) then
      if (     (node_list%node(i_node_side)%boundary .eq.  1) &
          .or. (node_list%node(i_node_side)%boundary .eq.  5) &
          ) then
        node_list%node(i_node)%boundary = 20
      endif
    endif
  enddo
  
  !!! RECAP !!! DEFINED BY GUIDO
  ! 1: TARGET,  side 2
  ! 2: TANGENT, side 3
  ! 3: CORNER,  between type-1 and type-2
  ! 4: TARGET,  side 2 (Same as type-1)
  ! 5: TARGET,  side 3
  ! 9: CORNER,  between type-4 and type-5
  
  !!! RECAP !!! DEFINED BY STAN TO INCLUDE FIELD DIRECTION, COMPATIBLE WITH OLD GRIDS
  ! 1: TARGET,  side 2                                   (inward  field)
  !11: TARGET,  side 2                                   (outward field)
  ! 5: TARGET,  side 3                                   (inward  field)
  !15: TARGET,  side 3                                   (outward field)
  ! 9: CORNER, free                                      (inward  field)
  !19: CORNER, free                                      (outward field)
  ! 2: TANGENT, side 3                                   (tangent field)
  !12: TANGENT, side 2                                   (tangent field)
  !20: CORNER, fixed                                     (tangent field)
  !21: CORNER, inverted (3 elements)                     (special case)
  ! NOT USED IN NEW GRID
  ! 3: Corner between type-1 and type-2                            (only for old grid)
  ! 4: Same as type-1                                              (only for old grid)
  
  !!! MAPPING !!! TO GO FROM STAN'S DEFINITION TO GUIDO'S
  !  1 ->  1
  ! 11 ->  1
  !  5 ->  5
  ! 15 ->  5
  !  9 ->  9
  ! 19 ->  9
  !  2 ->  2
  ! 12 ->  not defined!
  ! 20 ->  9 (because not defined!)
  ! 21 ->  not defined!
  !  3 ->  3


  ! --- Convert to Guido's definition
  if (use_simple_bnd_types) then
    do i_node=1,node_list%n_nodes
      if (node_list%node(i_node)%boundary .eq. 0 ) cycle
      if (node_list%node(i_node)%boundary .eq. 11) node_list%node(i_node)%boundary = 1
      if (node_list%node(i_node)%boundary .eq. 15) node_list%node(i_node)%boundary = 5
      if (node_list%node(i_node)%boundary .eq. 19) node_list%node(i_node)%boundary = 9
      if (node_list%node(i_node)%boundary .eq. 12) node_list%node(i_node)%boundary = 4
      if (node_list%node(i_node)%boundary .eq. 20) node_list%node(i_node)%boundary = 9
      if (node_list%node(i_node)%boundary .eq. 21) node_list%node(i_node)%boundary = 9
    enddo
  endif
    


  
  !!! OSOLETE !!! THE INITIAL DEFINITION BY STAN (INCOMPATIBLE WITH GUIDO)
  ! 1: TARGET,  side 2                                   (inward  field)
  ! 2: TANGENT, side 3                                   (tangent field)
  ! 3: CORNER,  free                                     (inward  field)
  ! 4: TANGENT, side 2                                   (tangent field)
  ! 5: TARGET,  side 2                                   (outward field)
  ! 6: TARGET,  side 3                                   (inward  field)
  ! 7: TARGET,  side 3                                   (outward field)
  ! 8: CORNER,  free                                     (outward field)
  ! 9: CORNER, inverted (3 elements)                     (special case)
  !10: CORNER, fixed                                     (tangent field)
  
  !!! MAPPING !!! TO GO FROM STAN'S OLD DEFINITION TO STAN'S NEW DEFINITION
  !  1 ->  1
  !  2 ->  2
  !  3 ->  9
  !  4 -> 12
  !  5 -> 11
  !  6 ->  5
  !  7 -> 15
  !  8 -> 19
  !  9 -> 21
  ! 10 -> 20
  
  return
  
end subroutine update_boundary_types_final










! --- This is the old version, which defines the field direction directly in the nodes types
! --- Keeping it temporatily just in case...
subroutine update_boundary_types_final_old(element_list,node_list)
  
  
  use constants
  use mod_parameters
  use data_structure
  use phys_module, only: xcase
  use mod_boundary
  
  implicit none
  
  ! --- Routine variables
  type(type_node_list),    intent(inout) :: node_list      !< list of grid nodes
  type(type_element_list), intent(inout) :: element_list   !< list of finite elements
  
  type (type_bnd_node_list)    :: bnd_node_list
  type (type_bnd_element_list) :: bnd_elm_list
  
  ! --- Internal variables
  integer :: i_elm,       i_vertex,  i_node,  i_node_prev, i
  integer :: i_elm2,      i_vertex2, i_node2, i_node_inside, i_node_side, i_node_side2
  integer :: elm_sum
  integer :: i_elm_first, i_vertex_first
  integer :: i_elm_now,   i_vertex_now, i_vertex_next, i_vertex_save
  integer :: i_elm_prev
  integer :: iter, n_xpoints, count
  integer :: i_elm_bnd(4), i_vertex_bnd(4)
  logical :: found_first, found_next, reverse
  logical, parameter :: debug  = .false.
  real*8  :: xjac
  real*8  :: psi_s, psi_t, psi_x, psi_y
  real*8  :: R, R_s, R_t
  real*8  :: Z, Z_s, Z_t
  real*8  :: BR, BZ
  real*8  :: norm_R,  norm_Z
  real*8  :: tang_R,  tang_Z
  real*8  :: tang_R2, tang_Z2
  real*8  :: alpha_Bp, alpha_norm, alpha_tang, alpha_tang2, alpha_tmp, alpha_between1, alpha_between2
  logical :: surface_is_tangent
  real*8, parameter  :: tol_tangent = 3.d-2 !1.5d-2 !1.d-2
  
  ! --- Some printouts?
  if (debug) then
    write(*,*)'----------------------------------------------------------'
    write(*,*)'---------------- Updating boundary types. ----------------'
    write(*,*)'----------------------------------------------------------'
  endif
  
  ! --- Loop over nodes
  do i_node=1,node_list%n_nodes
    if (node_list%node(i_node)%boundary .eq. 0) cycle
    ! --- Count how many elements we have on this node
    count = 0
    do i_elm=1,element_list%n_elements
      do i_vertex=1,4
        i_node2 = element_list%element(i_elm)%vertex(i_vertex)
        if (i_node2 .eq. i_node) then
          count = count + 1
          i_elm_bnd(count)  = i_elm
          i_vertex_bnd(count) = i_vertex
        endif
      enddo
    enddo

    ! --- psi and RZ variables
    R         = node_list%node(i_node)%x(1,1,1)
    R_s       = node_list%node(i_node)%x(1,2,1)
    R_t       = node_list%node(i_node)%x(1,3,1)
    Z         = node_list%node(i_node)%x(1,1,2)
    Z_s       = node_list%node(i_node)%x(1,2,2)
    Z_t       = node_list%node(i_node)%x(1,3,2)
    xjac      =  R_s*Z_t - R_t*Z_s
    psi_s     = node_list%node(i_node)%values(1,2,1)
    psi_t     = node_list%node(i_node)%values(1,3,1)
    psi_x     = (   Z_t * psi_s - Z_s * psi_t ) / xjac
    psi_y     = ( - R_t * psi_s + R_s * psi_t ) / xjac
    
    ! --- Poloidal field
    BR =  psi_y / R
    BZ = -psi_x / R
    alpha_Bp = atan2(BZ,BR)
    if (alpha_Bp .lt. 0.d0) alpha_Bp = alpha_Bp + 2.d0*PI
           

    ! --- If we have an outer corner
    if (count .eq. 1) then
      
      ! --- The insider
      i_node_inside = mod(i_vertex_bnd(1),4) + 1
      i_node_inside = element_list%element(i_elm_bnd(1))%vertex(i_node_inside)
      i_node_side   = mod(i_vertex_bnd(1)+2,4) + 1
      i_node_side   = element_list%element(i_elm_bnd(1))%vertex(i_node_side)
      
      ! --- The tangent vector
      tang_R = node_list%node(i_node_inside)%x(1,1,1) - node_list%node(i_node)%x(1,1,1)
      tang_Z = node_list%node(i_node_inside)%x(1,1,2) - node_list%node(i_node)%x(1,1,2)
      alpha_tang = atan2(tang_Z,tang_R)
      if (alpha_tang .lt. 0.d0) alpha_tang = alpha_tang + 2.d0*PI

      ! --- The other tangent vector
      tang_R = node_list%node(i_node_side)%x(1,1,1) - node_list%node(i_node)%x(1,1,1)
      tang_Z = node_list%node(i_node_side)%x(1,1,2) - node_list%node(i_node)%x(1,1,2)
      alpha_tang2 = atan2(tang_Z,tang_R)
      if (alpha_tang2 .lt. 0.d0) alpha_tang2 = alpha_tang2 + 2.d0*PI
      
      ! --- Swap angles to have good direction
      if (alpha_tang2 .lt. alpha_tang) then
        alpha_norm  = alpha_tang2
        alpha_tang2 = alpha_tang
        alpha_tang  = alpha_norm
      endif
      if (alpha_tang2-alpha_tang .gt. PI) then
        alpha_norm  = alpha_tang2
        alpha_tang2 = alpha_tang
        alpha_tang  = alpha_norm
        alpha_tang2 = alpha_tang2 + 2.d0*PI
      endif
      
      ! --- Now, are we inside or outside
      if (alpha_Bp .lt. alpha_tang) alpha_Bp = alpha_Bp + 2.d0*PI
      if ( (alpha_Bp .gt. alpha_tang) .and. (alpha_Bp .lt. alpha_tang2) ) then
        node_list%node(i_node)%boundary = 8
      else
        node_list%node(i_node)%boundary = 3
      endif
      
    ! --- If we have an inner corner (there we need to let is free, or fix it???)
    elseif (count .eq. 3) then
      
      node_list%node(i_node)%boundary = 9
    
    ! --- Standard boundary
    else
      ! --- The insider and tangent nodes
      i_node_inside = mod(i_vertex_bnd(1),4) + 1
      i_node_inside = element_list%element(i_elm_bnd(1))%vertex(i_node_inside)
      i_node_side   = mod(i_vertex_bnd(1)+2,4) + 1
      i_node_side   = element_list%element(i_elm_bnd(1))%vertex(i_node_side)
      if (node_list%node(i_node_inside)%boundary .ne. 0) then
        i_node_inside = mod(i_vertex_bnd(1)+2,4) + 1
        i_node_inside = element_list%element(i_elm_bnd(1))%vertex(i_node_inside)
        i_node_side   = mod(i_vertex_bnd(1),4) + 1
        i_node_side   = element_list%element(i_elm_bnd(1))%vertex(i_node_side)
      endif
      ! --- The other tangent nodes
      i_node_side2 = mod(i_vertex_bnd(2)+2,4) + 1
      i_node_side2 = element_list%element(i_elm_bnd(2))%vertex(i_node_side2)
      if (node_list%node(i_node_side2)%boundary .eq. 0) then
        i_node_side2 = mod(i_vertex_bnd(2),4) + 1
        i_node_side2 = element_list%element(i_elm_bnd(2))%vertex(i_node_side2)
      endif
      
      ! --- The normal vector
      norm_R = node_list%node(i_node_inside)%x(1,1,1) - node_list%node(i_node)%x(1,1,1)
      norm_Z = node_list%node(i_node_inside)%x(1,1,2) - node_list%node(i_node)%x(1,1,2)

      ! --- The tangent vector
      tang_R = node_list%node(i_node_side)%x(1,1,1) - node_list%node(i_node)%x(1,1,1)
      tang_Z = node_list%node(i_node_side)%x(1,1,2) - node_list%node(i_node)%x(1,1,2)

      ! --- The other tangent vector
      tang_R2 = node_list%node(i_node_side2)%x(1,1,1) - node_list%node(i_node)%x(1,1,1)
      tang_Z2 = node_list%node(i_node_side2)%x(1,1,2) - node_list%node(i_node)%x(1,1,2)

      ! --- The angles
      alpha_norm = atan2(norm_Z,norm_R)
      if (alpha_norm .lt. 0.d0) alpha_norm = alpha_norm + 2.d0*PI
      alpha_tang = atan2(tang_Z,tang_R)
      if (alpha_tang .lt. 0.d0) alpha_tang = alpha_tang + 2.d0*PI
      alpha_tang2 = atan2(tang_Z2,tang_R2)
      if (alpha_tang2 .lt. 0.d0) alpha_tang2 = alpha_tang2 + 2.d0*PI
      if (alpha_tang2 .lt. alpha_tang) then
        alpha_tmp   = alpha_tang
        alpha_tang  = alpha_tang2
        alpha_tang2 = alpha_tmp
      endif
      
      if (alpha_norm .lt. alpha_tang) alpha_norm = alpha_norm + 2.d0*PI
      if (alpha_Bp   .lt. alpha_tang) alpha_Bp   = alpha_Bp   + 2.d0*PI
      
      ! Determine if surface is tangent to boundary
      surface_is_tangent = .false.
      if (abs(alpha_Bp-alpha_tang      ) .lt. tol_tangent*2*PI) surface_is_tangent = .true.
      if (abs(alpha_Bp-alpha_tang -2*PI) .lt. tol_tangent*2*PI) surface_is_tangent = .true.
      if (abs(alpha_Bp-alpha_tang2     ) .lt. tol_tangent*2*PI) surface_is_tangent = .true.
      if (abs(alpha_Bp-alpha_tang2-2*PI) .lt. tol_tangent*2*PI) surface_is_tangent = .true.
      alpha_between1 = alpha_tang2 - PI
      alpha_between2 = alpha_tang
      if (alpha_between2 .lt. alpha_between1) then
        alpha_tmp      = alpha_between1
        alpha_between1 = alpha_between2
        alpha_between2 = alpha_tmp
      endif
      if ( (alpha_between1 .le. alpha_Bp     ) .and. (alpha_Bp      .le. alpha_between2) ) surface_is_tangent = .true.
      if ( (alpha_between1 .le. alpha_Bp-1*PI) .and. (alpha_Bp-1*PI .le. alpha_between2) ) surface_is_tangent = .true.
      if ( (alpha_between1 .le. alpha_Bp-2*PI) .and. (alpha_Bp-2*PI .le. alpha_between2) ) surface_is_tangent = .true.
      
      ! --- If the angle is very small, we want to fix conditions (like along a flux surface)
      if (surface_is_tangent) then
        if (node_list%node(i_node)%boundary .eq. 1) node_list%node(i_node)%boundary = 4
        ! note the bnd 2 is already this
      ! --- If the angle is large, we need Bohm
      else
        reverse = .false.
        if ((alpha_norm .lt. alpha_tang2) .and. (alpha_Bp .lt. alpha_tang2)) reverse = .true.
        if ((alpha_norm .gt. alpha_tang2) .and. (alpha_Bp .gt. alpha_tang2)) reverse = .true.
        
        if (reverse) then
          if (node_list%node(i_node)%boundary .eq. 1) then
            node_list%node(i_node)%boundary = 5
          else
            node_list%node(i_node)%boundary = 7
          endif
        else
          if (node_list%node(i_node)%boundary .ne. 1) then
            node_list%node(i_node)%boundary = 6
          endif
        endif
        
      endif
      
    endif
  enddo
  
  
  ! --- Loop over nodes again to check corners
  do i_node=1,node_list%n_nodes
    if ( (node_list%node(i_node)%boundary .ne. 3) .and. (node_list%node(i_node)%boundary .ne. 8) ) cycle
    
    ! --- Count how many elements we have on this node
    count = 0
    do i_elm=1,element_list%n_elements
      do i_vertex=1,4
        i_node2 = element_list%element(i_elm)%vertex(i_vertex)
        if (i_node2 .eq. i_node) then
          count = count + 1
          i_elm_bnd(count)  = i_elm
          i_vertex_bnd(count) = i_vertex
        endif
      enddo
    enddo
    
    
    ! --- The side nodes
    i_node_inside = mod(i_vertex_bnd(1),4) + 1
    i_node_inside = element_list%element(i_elm_bnd(1))%vertex(i_node_inside)
    i_node_side   = mod(i_vertex_bnd(1)+2,4) + 1
    i_node_side   = element_list%element(i_elm_bnd(1))%vertex(i_node_side)
    
    ! --- Check consistency
    if (     (node_list%node(i_node_inside)%boundary .eq. 3) &
        .or. (node_list%node(i_node_inside)%boundary .eq. 8) &
        .or. (node_list%node(i_node_inside)%boundary .eq. 9) &
        ) then
      node_list%node(i_node)%boundary = 10
    endif
    if (     (node_list%node(i_node_side)%boundary .eq. 3) &
        .or. (node_list%node(i_node_side)%boundary .eq. 8) &
        .or. (node_list%node(i_node_side)%boundary .eq. 9) &
        ) then
      node_list%node(i_node)%boundary = 10
    endif
    if (     (node_list%node(i_node_inside)%boundary .eq. 1) &
        .or. (node_list%node(i_node_inside)%boundary .eq. 6) &
        ) then
      if (     (node_list%node(i_node_side)%boundary .eq. 5) &
          .or. (node_list%node(i_node_side)%boundary .eq. 7) &
          ) then
        node_list%node(i_node)%boundary = 10
      endif
    endif
    if (     (node_list%node(i_node_inside)%boundary .eq. 5) &
        .or. (node_list%node(i_node_inside)%boundary .eq. 7) &
        ) then
      if (     (node_list%node(i_node_side)%boundary .eq. 1) &
          .or. (node_list%node(i_node_side)%boundary .eq. 6) &
          ) then
        node_list%node(i_node)%boundary = 10
      endif
    endif
  enddo
  
  !!! RECAP!!!
  ! 1: Non-corner target, with tangent on vector 2, inward field
  ! 2: Non-corner flux-aligned boundary, with tangent on vector 3
  ! 3: Corner target, inward field
  ! 4: Non-corner flux-aligned boundary, with tangent on vector 2
  ! 5: Non-corner target, with tangent on vector 2, outward field
  ! 6: Non-corner target, with tangent on vector 3, inward field
  ! 7: Non-corner target, with tangent on vector 3, outward field
  ! 8: Corner target, outward field
  ! 9: Inverted corner target (3 elements)
  !10: Corner target, tangent field
  
     
  
  return
  
end subroutine update_boundary_types_final_old


