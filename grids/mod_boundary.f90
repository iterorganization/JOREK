!> Allows to generate the boundary-element and boundary-node data structures for a given grid.
module mod_boundary

  implicit none
  private
  public boundary_from_grid, log_bnd_info
  public wall_normal_vector

  
  
  
  
  
  contains
  
  
  
  
  
  !> Routine extracts the boundary information (boundary element and node lists)
  !! from the information stored in the grid (element and node lists).
  !! 
  !! Note: Grid nodes with boundary=3 (located at the edges of the boundary in
  !!       the divertor region) are intentionally added twice to the bnd_node_list.
  subroutine boundary_from_grid(node_list,element_list,bnd_node_list,bnd_elm_list,infos)

    use data_structure

    implicit none

    type (type_node_list),        intent(inout) :: node_list
    type (type_element_list),     intent(in)    :: element_list
    type (type_bnd_node_list),    intent(inout) :: bnd_node_list
    type (type_bnd_element_list), intent(inout) :: bnd_elm_list
    logical,                      intent(in)    :: infos ! If .true., writes bnd nodes and bnd elements in files 'boundary_nodes.dat' and 'boundaru_elements.dat' 

    integer :: i_elem         ! Element index in element_list.
    integer :: iside          ! Side of the element (1...4).
    integer :: iv1, iv2       ! Numbers of the two nodes on the element side (1...4).
    integer :: inode1, inode2, inode  ! Indices of the nodes in the node_list.
    integer :: b1, b2         ! Boundary properties of the two nodes.
    integer :: i, j, index_bnd, i_dof
    
    ! --- Empty the boundary node and element lists.
    bnd_node_list%n_bnd_nodes   = 0
    bnd_elm_list%n_bnd_elements = 0

    ! --- Go through all elements of the grid and find out which are located
    !     at the boundary. Build up the lists of boundary nodes and elements.

    do i_elem = 1, element_list%n_elements ! for each element

      do iside = 1, 4 ! for each side of the element (1...4)

        ! --- Indices of the nodes of the current element, which are located on this side.
        iv1 = iside
        iv2 = mod( iside, 4 ) + 1

        ! --- Indices of these nodes in the node_list.
        inode1 = element_list%element(i_elem)%vertex(iv1)
        inode2 = element_list%element(i_elem)%vertex(iv2)

        ! --- Boundary information of the two nodes.
        b1 = node_list%node(inode1)%boundary
        b2 = node_list%node(inode2)%boundary

        ! --- If the current side of this element is located at the boundary,
        !     add a new boundary element.
        if ( (b1 > 0) .AND. (b2 > 0) ) then
          call add_bnd_elem( i_elem, iv1, iv2, inode1, inode2, iside, b1, b2, &
                             element_list, bnd_node_list, bnd_elm_list )
        end if

      end do

    end do

    ! --- Sort the boundary elements.
    call sort_bnd_elements( bnd_elm_list )


!=================== definition of boundary_index
    do i=1, bnd_node_list%n_bnd_nodes
       inode = bnd_node_list%bnd_node(i)%index_jorek
       node_list%node(inode)%boundary_index = i          ! the index in the bnd_node_list (not response_matrix)
    end do

!=================== define index in the response matrix (including corners whose nodes appear twice)

    index_bnd = 0

    do i=1, bnd_node_list%n_bnd_nodes

      bnd_node_list%bnd_node(i)%index_starwall = (/index_bnd+1,index_bnd+2/)
      bnd_node_list%bnd_node(i)%n_dof = 2

      index_bnd = index_bnd + bnd_node_list%bnd_node(i)%n_dof

    !  write(*,*)  'bnd index : ',i, bnd_node_list%bnd_node(i)%index_starwall

    enddo

    ! --- Output short (.false.) or verbose (.true.) boundary information.
    call log_bnd_info(infos, node_list, bnd_node_list, bnd_elm_list)

  end subroutine boundary_from_grid
  
  
  
  
  
  
  !> Adds a boundary element to the boundary element list.
  subroutine add_bnd_elem( i_elem, iv1, iv2, inode1, inode2, iside, b1, b2, &
    element_list, bnd_node_list, bnd_elm_list)

    use data_structure

    implicit none

    integer, intent(in)    :: i_elem         ! Index of the element in the element_list.
    integer, intent(in)    :: iv1, iv2       ! Index of the nodes of the element i_elem.
    integer, intent(in)    :: inode1, inode2 ! Indices of nodes in node_list.
    integer, intent(in)    :: iside          ! Which side of the element is located at the boundary?
    integer, intent(in)    :: b1, b2         ! Boundary information of the two nodes.
    type (type_element_list),     intent(in)    :: element_list
    type (type_bnd_node_list),    intent(inout) :: bnd_node_list
    type (type_bnd_element_list), intent(inout) :: bnd_elm_list

    integer :: ib1, ib2       ! Indices of the nodes in the bnd_node_list.
    integer :: idir

    ! --- Add the two boundary nodes to the list.
    call add_bnd_node( i_elem, iv1, inode1, iside, b1, ib1, bnd_node_list )
    call add_bnd_node( i_elem, iv2, inode2, iside, b2, ib2, bnd_node_list )

    bnd_elm_list%n_bnd_elements = bnd_elm_list%n_bnd_elements + 1
    
    if ( bnd_elm_list%n_bnd_elements > n_boundary_max ) then
      write(*,*) 'ERROR in mod_boundary:boundary_from_grid: hard-coded parameter n_boundary_max is too small'
      stop
    end if

    ! --- Store vertex indices belonging to the boundary element.
    bnd_elm_list%bnd_element(bnd_elm_list%n_bnd_elements)%vertex = (/ inode1, inode2 /)

    ! --- Store side of the element that is at the boundary.
    bnd_elm_list%bnd_element(bnd_elm_list%n_bnd_elements)%side   = iside

    ! --- Store, which element the boundary element belongs to.
    bnd_elm_list%bnd_element(bnd_elm_list%n_bnd_elements)%element = i_elem

    ! --- Store direction (directly connected to %side)
    bnd_elm_list%bnd_element(bnd_elm_list%n_bnd_elements)%direction(:,1) = 1
    if ( mod(iside,2) == 1 ) then
      idir = 2
    else
      idir = 3
    end if
    bnd_elm_list%bnd_element(bnd_elm_list%n_bnd_elements)%direction(:,2) = idir

    ! --- Store element size.
    bnd_elm_list%bnd_element(bnd_elm_list%n_bnd_elements)%size(1,1) = element_list%element(i_elem)%size(iv1,1)
    bnd_elm_list%bnd_element(bnd_elm_list%n_bnd_elements)%size(2,1) = element_list%element(i_elem)%size(iv2,1)

    bnd_elm_list%bnd_element(bnd_elm_list%n_bnd_elements)%size(1,2) = element_list%element(i_elem)%size(iv1,idir)
    bnd_elm_list%bnd_element(bnd_elm_list%n_bnd_elements)%size(2,2) = element_list%element(i_elem)%size(iv2,idir)

    ! --- Store boundary node indices.
    bnd_elm_list%bnd_element(bnd_elm_list%n_bnd_elements)%bnd_vertex(1) = ib1
    bnd_elm_list%bnd_element(bnd_elm_list%n_bnd_elements)%bnd_vertex(2) = ib2

  end subroutine add_bnd_elem
  
  
  
  
  
  
  !> Adds a node to the bnd_node_list avoiding duplicates (except for boundary=3)
  !! and returns the index of the node in the bnd_node_list as bnd_vertex.
  subroutine add_bnd_node(i_elem, iv, inode, iside, boundary, bnd_vertex, bnd_node_list)

    use data_structure

    implicit none

    integer, intent(in)    :: i_elem         ! Index of the element in the element_list.
    integer, intent(in)    :: iv             ! Index of the node of the element i_elem.
    integer, intent(in)    :: inode          ! Index of the node in the node_list.
    integer, intent(in)    :: iside          ! Which side of the element is located at the boundary?
    integer, intent(in)    :: boundary       ! Boundary property of the node.
    integer, intent(inout) :: bnd_vertex     ! Index of the node in the bnd_node_list.
    type (type_bnd_node_list), intent(inout) :: bnd_node_list

    integer :: i, idir

    ! --- Make sure the node is not in the boundary node list yet (except for boundary=3).
    if ( boundary /= 3 ) then
      do i = 1, bnd_node_list%n_bnd_nodes
        if ( bnd_node_list%bnd_node(i)%index_jorek == inode ) then
          bnd_vertex = i ! Node is already in the list, return its index.
          return
        end if
      end do
    end if

    ! --- Add the node.
    bnd_node_list%n_bnd_nodes = bnd_node_list%n_bnd_nodes + 1
    bnd_vertex = bnd_node_list%n_bnd_nodes

    bnd_node_list%bnd_node(bnd_vertex)%index_jorek    = inode
    bnd_node_list%bnd_node(bnd_vertex)%index_starwall = bnd_vertex

    ! --- Store direction (directly connected to iside)
    bnd_node_list%bnd_node(bnd_vertex)%direction(1) = 1
    if ( mod(iside,2) == 1 ) then
      idir = 2
    else
      idir = 3
    end if
    bnd_node_list%bnd_node(bnd_vertex)%direction(2) = idir

  end subroutine add_bnd_node
  
  
  
  
  
  
  !> Outputs information about the boundary elements and nodes.
  subroutine log_bnd_info(verbose, node_list, bnd_node_list, bnd_elm_list, dir_in, filename_appendix_in)

    use data_structure

    implicit none

    logical, intent(in) :: verbose ! Output verbose information?
    type (type_node_list),        intent(in) :: node_list
    type (type_bnd_node_list),    intent(in) :: bnd_node_list
    type (type_bnd_element_list), intent(in) :: bnd_elm_list
    character(len=*), optional,   intent(in) :: dir_in
    character(len=*), optional,   intent(in) :: filename_appendix_in

    integer             :: i
    character(len=20)   :: s
    character(len=1024) :: dir, filename_appendix
    
    120 format(3x,77('-'))
    121 format(3X,A,I10,A)
    141 format(5X,A,I10,A)
    161 format(7X,A,I10,A)
    182 format(9X,A,15I10)
    183 format(9X,A,15ES10.2)

    write(*,*)
    write(*,120)
    write(*,121) 'BOUNDARY INFORMATION:'
    write(*,120)
    write(*,141) 'n_bnd_elements=', bnd_elm_list%n_bnd_elements
    write(*,141) 'n_bnd_nodes   =', bnd_node_list%n_bnd_nodes
    write(*,120)

    if ( verbose ) then
      
      dir               = './'
      filename_appendix = '.dat'
      if ( present(dir_in)               ) dir               = dir_in
      if ( present(filename_appendix_in) ) filename_appendix = filename_appendix_in
      
!      write(*,*)
!      write(*,120)
!      write(*,141) 'BOUNDARY ELEMENTS:'
!      write(*,120)
!      do i = 1, bnd_elm_list%n_bnd_elements
!        write(s,*) i
!        write(*,161) '#'//trim(adjustl(s))//':'
!        write(*,182) 'vertex        =', bnd_elm_list%bnd_element(i)%vertex
!        write(*,182) 'bnd_vertex    =', bnd_elm_list%bnd_element(i)%bnd_vertex
!        write(*,182) 'direction     =', bnd_elm_list%bnd_element(i)%direction
!        write(*,182) 'element       =', bnd_elm_list%bnd_element(i)%element
!        write(*,182) 'side          =', bnd_elm_list%bnd_element(i)%side
!        write(*,183) 'size          =', bnd_elm_list%bnd_element(i)%size
!      end do
!      write(*,120)
!      write(*,*)
!      
!      write(*,120)
!      write(*,141) 'BOUNDARY NODES:'
!      write(*,120)
!      do i = 1, bnd_node_list%n_bnd_nodes
!        write(s,*) i
!        write(*,161) '#'//trim(adjustl(s))//':'
!        write(*,182) 'index_jorek   =', bnd_node_list%bnd_node(i)%index_jorek
!        write(*,182) 'index_starwall=', bnd_node_list%bnd_node(i)%index_starwall
!        write(*,182) 'direction     =', bnd_node_list%bnd_node(i)%direction
!      end do
!      write(*,120)
!      write(*,*)
      
      write(*,*)
      write(*,*) 'Writing boundary element node coordinates to the following file:'
      write(*,*) trim(DIR) // '/boundary_element_nodes' // trim(filename_appendix)
      open(42, file=trim(DIR) // '/boundary_element_nodes' // trim(filename_appendix), status='replace', action='write')
      do i = 1, bnd_elm_list%n_bnd_elements
        write(42,*) node_list%node( bnd_elm_list%bnd_element(i)%vertex(1) )%x(1,:)
        write(42,*) node_list%node( bnd_elm_list%bnd_element(i)%vertex(2) )%x(1,:)
        write(42,*)
        write(42,*)
      end do
      close(42)
      
      write(*,*)
      write(*,*) 'Writing boundary node coordinates to the following file:'
      write(*,*) trim(DIR) // '/boundary_nodes' // trim(filename_appendix)
      open(42, file=trim(DIR) // '/boundary_nodes' // trim(filename_appendix), status='replace', action='write')
      do i = 1, bnd_node_list%n_bnd_nodes
        write(42,*) node_list%node( bnd_node_list%bnd_node(i)%index_jorek )%x(1,:)
      end do
      close(42)
      
      write(*,*)
      write(*,*) 'Writing boundary element details to the following file:'
      write(*,*) trim(DIR) // '/boundary_element_details' // trim(filename_appendix)
      open(42, file=trim(DIR) // '/boundary_element_details' // trim(filename_appendix), status='replace', action='write')
      write(42,'(a)') '#   bnd_element        vertex1        vertex2    bnd_vertex1    bnd_vertex2   direction1,1   ' // &
        'direction1,2   direction2,1   direction2,2        element           side      boundary1      boundary2     ' // &
        '        R1             Z1             R2             Z2'
      do i = 1, bnd_elm_list%n_bnd_elements
        777 format(13i15,4f15.8)
        write(42,777) i, bnd_elm_list%bnd_element(i)%vertex(:), bnd_elm_list%bnd_element(i)%bnd_vertex(:),   &
          bnd_elm_list%bnd_element(i)%direction(:,:), bnd_elm_list%bnd_element(i)%element,                   &
          bnd_elm_list%bnd_element(i)%side, node_list%node(bnd_elm_list%bnd_element(i)%vertex(1))%boundary,  &
          node_list%node(bnd_elm_list%bnd_element(i)%vertex(2))%boundary,                                    &
          node_list%node(bnd_elm_list%bnd_element(i)%vertex(1))%x(1,:),                                      &
          node_list%node(bnd_elm_list%bnd_element(i)%vertex(1))%x(2,:)
      end do
      close(42)
      
    end if

    write(*,*)

  end subroutine log_bnd_info
  
  
  
  
  
  
  !> Sorts the boundary elements.
  subroutine sort_bnd_elements( bnd_elm_list )

    use data_structure

    implicit none

    integer :: current_vertex, iter, iter_max
    integer :: ibnd_elem
    type(type_bnd_element)      :: bnd_elem
    type(type_bnd_element_list) :: sorted_bnd_element_list
    type (type_bnd_element_list), intent(inout) :: bnd_elm_list
    logical ::found_neighbour

    sorted_bnd_element_list%n_bnd_elements = 0

    current_vertex = bnd_elm_list%bnd_element(1)%vertex(1)

    iter = 0
    iter_max = 2*bnd_elm_list%n_bnd_elements
    found_neighbour = .true.

    L_A: do while ( (bnd_elm_list%n_bnd_elements > 0) .and. found_neighbour)

      iter = iter + 1
      if (iter .gt. iter_max) exit

      found_neighbour = .false.

      do ibnd_elem = 1, bnd_elm_list%n_bnd_elements

        bnd_elem = bnd_elm_list%bnd_element(ibnd_elem)

        if ( bnd_elem%vertex(2) == current_vertex ) call reverse_elem( bnd_elem )

        if ( bnd_elem%vertex(1) == current_vertex ) then
          current_vertex = bnd_elem%vertex(2)
          call add_elem( bnd_elem, sorted_bnd_element_list )
          call remove_elem( ibnd_elem, bnd_elm_list )
          found_neighbour = .true.
          exit
        end if

      end do

    end do L_A

    bnd_elm_list = sorted_bnd_element_list

    if (.not. found_neighbour) then
      write(*,*) 'BOUNDARY WARNING : NO NEIGHBOUR BND ELEMENT!',current_vertex
    endif

  end subroutine sort_bnd_elements
  
  
  
  
  
  
  !> Add the given boundary element to the end of the sorted list.
  subroutine add_elem( bnd_elem, sorted_bnd_element_list )

    use data_structure

    implicit none

    type(type_bnd_element),      intent(in)    :: bnd_elem
    type(type_bnd_element_list), intent(inout) :: sorted_bnd_element_list

    integer :: ibnd_elem

    ibnd_elem = sorted_bnd_element_list%n_bnd_elements + 1
    sorted_bnd_element_list%n_bnd_elements = ibnd_elem
    sorted_bnd_element_list%bnd_element(ibnd_elem) = bnd_elem

  end subroutine add_elem
  
  
  
  
  
  
  !> Reverse the given boundary element, i.e., exchange the nodes.
  subroutine reverse_elem( bnd_elem )

    use data_structure

    implicit none

    type(type_bnd_element), intent(inout) :: bnd_elem

    type(type_bnd_element) :: reversed_elem

    !write(*,*) '############################################################################'
    write(*,*) 'REVERSE_ELEM'
    !write(*,*) '############################################################################'

    reversed_elem%vertex(1) = bnd_elem%vertex(2)
    reversed_elem%vertex(2) = bnd_elem%vertex(1)

    reversed_elem%bnd_vertex(1) = bnd_elem%bnd_vertex(2)
    reversed_elem%bnd_vertex(2) = bnd_elem%bnd_vertex(1)

    reversed_elem%direction(1,:) = bnd_elem%direction(2,:)
    reversed_elem%direction(2,:) = bnd_elem%direction(1,:)

    reversed_elem%size(1,:) = bnd_elem%size(2,:)
    reversed_elem%size(2,:) = bnd_elem%size(1,:)
    
    bnd_elem = reversed_elem

  end subroutine reverse_elem
  
  
  
  
  
  
  !> Remove the ibnd_elem-th boundary element from bnd_elm_list.
  subroutine remove_elem( ibnd_elem, bnd_elm_list )

    use data_structure

    implicit none

    integer, intent(in) :: ibnd_elem
    type (type_bnd_element_list), intent(inout) :: bnd_elm_list

    integer :: nbnd_elem

    nbnd_elem = bnd_elm_list%n_bnd_elements

    bnd_elm_list%bnd_element(ibnd_elem:nbnd_elem-1) = bnd_elm_list%bnd_element(ibnd_elem+1:nbnd_elem)

    bnd_elm_list%n_bnd_elements = nbnd_elem - 1

  end subroutine remove_elem


  !> Get a normal vector from the wall on your element. It is your own responsibility to ensure
  !> that it is actually on the wall ;) We just find the closest edge of the element and calculate the gradient
  !> towards the inside of the element.
  pure function wall_normal_vector(node_list, element_list, i_elm, s, t) result(n)
    use data_structure
    use mod_interp, only: interp_RZ
    type(type_node_list), intent(in)    :: node_list
    type(type_element_list), intent(in) :: element_list
    integer, intent(in)                 :: i_elm
    real*8, intent(in)                  :: s, t
    real*8                              :: n(3)
    real*8 :: R, R_s, R_t, Z, Z_s, Z_t
    logical :: allowed(4)
    integer :: i_side

    ! Calculate the two vectors that span the surface element
    ! One of them is always in the positive phi direction
    ! the poloidal vector is constructed from the local coordinate space
    ! Taking the cross product of those vectors gives the normal vector direction
    ! after normalization we now have the normal vector of the surface where the particle left the grid
        
    ! interp RZ for (Rs,Zs) and (Rt,Zt)
    call interp_RZ(node_list, element_list, i_elm, s, t, R, R_s, R_t, Z, Z_s, Z_t)    
        
    !> if statements for checking the side, to construct a normal vector on the element
    !> Explanation for the first two if statements
    !> if s>t, is the area below line from bottom left to right top.
    !> and then, if (1-s)>t, is the other diagonal line. So this specifies the lower triangle and we appoint to according side to it
    !> ________1
    !> |\  3 /|
    !> | \  / |
    !> |4 \/ 2|
    !> |  /\  | t       
    !> | /  \ |                     
    !> |/  1 \|0                      
    !> 0^^^^^^1 
    !>    s  
    ! First select the side to use.
    ! eligible sides are those with neighbours equal to 0.
    ! if there are multiple the one which is closest to s,t is chosen
    ! the distance is given by t, 1-s, 1-t, s respectively
    allowed = element_list%element(i_elm)%neighbours(:) .eq. 0
    if (.not. any(allowed)) allowed = .true. ! shortcut
    i_side = minloc([t, 1.d0-s, 1.d0-t, s], mask=allowed, dim=1)
    ! The normal can be found by norm(grad_R x) (where x = either s or t depending on the side)
    ! A proper coordinate transfrom has to be caried out i.e. grad_R x = grad_s x J_s^-1, where J_s^-1 = J_R is the Jacobian matrix
    ! Note that the normal vector does not have to be multiplied with the inverse of the determinant of the Jacobian matrix, 
    ! because it is normalized directly after.
    select case (i_side)
    case (1)
      n = [-Z_s, R_s, 0.d0]
    case (2)
      n = [-Z_t, R_t, 0.d0]
    case (3)
      n = [Z_s, -R_s, 0.d0]
    case (4)
      n = [Z_t, -R_t, 0.d0]
    end select
    n = n/norm2(n,1)
  end function wall_normal_vector
end module mod_boundary
