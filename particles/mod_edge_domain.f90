!> mod edge domain, module for making a list of the edge domains that will be used to sputter on
!>
!> The overall process is as follows:
!> loop over all jorek elements (node list and element list)
!> retrieve user defined sides to work with + RZ cutoff
!> check per side if neighbours = 0, keep if so (then we are at the edge)
!> these elements will form a new list of elements on the boundary
!> use function to check if geometrical average is above the cutoff line
!>
!> Now we have a list of all boundary elements we want to use for the edge domains (unsorted in this step)
!> The next step is to sort and connect these
!>
!> Subroutine find adjacent elements
!> Do a sorting loop to find adjacent element via coupling nodes 
!>
!> The main output type is type_edge_domain
!> with size M
!> this contains M edge domains with lists of i_elm i_side
!> from this xyz and ien +scalars
module mod_edge_domain
use data_structure

implicit none
private
public :: type_edge_domain
public :: find_edge_domains
public :: type_edge_element !< internal variable, public for testing/debug purposes

! short description of the edge_domain list and its purpose:
! This mod will produce a set of lists. Each list contains a configurable subset of the jorek elements and their sides that sit on the edge of the grid in the poloidal plane
! The user can specify which parts of the edge in the poloidal plane will be selected to sputter on.
! The elements are split up in lists (subsets). This can be used e.g. to specify which part of the wall is what material.

! the user can specify:
!  -  which sides of elements may be used (sides = (1,2,3,4)
!  -  a RZ-cutoff line, everything above (or below) the line won't be added to the edge domain (see function below)
!  -  if (sharp) corners of the grid should be disconnected (after each corner there will be a new list)

!> Describes a section of the edge of the JOREK domain.
!> i_elm and i_side should be of equal sizes.
!> Note that this should describe a *connected* region, i.e. when traversing
!> i_elm in order every JOREK element needs to have one node in common with its
!> predecessor in the list.
type :: type_edge_domain
  integer, allocatable, dimension(:) :: i_elm  !< JOREK element number (associating with edge element)
  integer, allocatable, dimension(:) :: i_side !< which side of the JOREK element to take into account
  integer, allocatable, dimension(:) :: direction !< ccw or cw elements | 1 = ccw | -1 = cw
end type type_edge_domain

!> type edge element is has 7 parameters on each element describing the desired edge domain.
!> This is an internal type used to store some information necessary for the reconstruction
type type_edge_element
  integer :: i_elm  !< element number in JOREK
  integer :: i_side !< which of the element sides is on the boundary
  integer :: node_1 !< first node on the boundary side
  integer :: node_2 !< second node on the boundary side
  integer :: i_list !< number of edge domain
  integer :: i !< number of element in domain in order
  integer :: direction !< is the element counter-clockwise or clockwise | 1 = ccw | -1 = cw
end type type_edge_element

contains 

!> function that will produce a user specified line where only the boundary elements above or below will be taking into the eddge_domain
!> User selects 2 points in RZ space, the 2 points form a line and only everything below (or right) to the line will be taken as a valid edge_domain
!> this way only the divertor region can be taken into account without having to use all boundary elements
!> ```
!>Z |          /\\\\\\\\\
!>  |         * (R2, Z2)\
!>  |        /\\\\\\\\\\\
!>  |       /\\\\\\\\\\\\
!>  |      /\\\\\\\\\\\\\
!>  |     * (R1,Z1)\\\\\\
!>  |    /\\\\\\\\\\\\\\\
!>  |___/_\_\_\_\_\_\_\_\  R
!> ```
!> flip_ud to choose above or under the line (flip_ud = false, left /  above the line will be cut off)
!> this option will flip the logical value of cut_off
function check_cutoff_edge_domain(x, x1, x2, flip_ud) result(cut_off) !< checks per element if it should be cut off
  real*8, intent(in), dimension(2)   :: x, x1, x2 !< position and 2 reference points
  logical, intent(in)                :: flip_ud !< side selection of the line x1, x2. false -> left/above is discarded
  logical                            :: cut_off !< true if this point should not be included in the edge domain
  real*8                             :: slope, z_offset, cutoff_line, R1, Z1, R2,Z2
 
  R1 = x1(1) ! just to make it more readable
  Z1 = x1(2) ! x1 = [R1, Z1]
  R2 = x2(1) ! x2 = [R1, Z1]
  Z2 = x2(2)
  
  !> x(1,1,1) = R and x(1,1,2) = Z
  !> constructing RZ cutoff line
  !> with form y = ax + b = slope * R + z_ofset
  if (abs(R1 - R2) .le. 1d-10) then
    ! write no cutoff
    ! otherwise division by 0
    ! When RZ points left empty by user, R1 = R2 = 0, then no cut_off
    cut_off = .false.
  else
     slope = (Z2-Z1)/(R2-R1)
     z_offset = (Z1 - slope*R1)
     cutoff_line = (slope)*x(1) + z_offset
     cut_off = x(2) .gt. cutoff_line !< compares the Z coordinates
  endif
  cut_off = cut_off .neqv. flip_ud !< same a .ne. this switches only the value of cut_off if flip_ud is true
  !< if true element will be cut off
end function check_cutoff_edge_domain


!> this subroutine searches for the first element where i_list = 0 as 0 is the minimum is this list
!> then it updates the edge_element%i_list, element_i, edge_element%i to be set up for the next iteration in the sorting loop
!> This first unindexed element is used as the starting point to find adjacent elements from
subroutine find_new_edge_domain_start(i_elm, i_list, current_element, i, i_side, direction, all_found, element_list, node_list)
  use mod_sampling, only: cross_product
  
  integer, intent(in), dimension(:)    :: i_elm, i_side
  integer, intent(inout), dimension(:) :: i_list, i, direction   
  integer, intent(inout)               :: current_element
  logical, intent(out)                 :: all_found
  type(type_node_list), intent(in)     :: node_list
  type(type_element_list), intent(in)  :: element_list

  integer                              :: new_index, node_1, node_2, node_3
  real*8, dimension(3)                 :: c, a_12, b_13
  type(type_node)                      :: node
  type(type_element)                   :: element
  
  new_index = minloc(i_list, dim = 1) !new index is the location in the array of the node for next loop
  if (new_index .lt. lbound(i_list,1)) then
    write(*,*) 'empty list in find_new_edge_domain start, exiting'
    all_found=.true.
    return
  end if
  
  !Error handler: if min(i_list) /= 0, all elements are listed. Go out of the do while loop.
  if (i_list(new_index) .ne. 0) then 
    all_found=.true. !< output that says all domains have been found and labeled.
    return
  else 
    all_found=.false.
    ! i_list(new_index)         = i_list(current_element)+1
    i_list(new_index)         = maxval(i_list, dim = 1)+1
    current_element           = new_index !update element_i
    i(current_element)        = 1      !new indexation for new list
    ! direction(current_element)= 1
    
    !check for absolute direction
    !this will determine the direction for the complete list as it is internally consistent

    node_1 = element_list%element(i_elm(current_element))%vertex(i_side(current_element))
    node_2 = element_list%element(i_elm(current_element))%vertex(modulo(i_side(current_element),4)+1)
    node_3 = element_list%element(i_elm(current_element))%vertex(modulo(modulo(i_side(current_element),4)+1,4)+1)

    !> This part checks for absolute direction ( ccw vs cw)
    !> It takes the first element of the new list, and calculates its absolute direction.
    !> It takes the 2 nodes on the edge side of the element and takes a 3rd node that is next in order
    !> vector a_12 is the vector from node_1 to node_2
    !> vector b_13 is the vector from node_1 to node_3
    !> by taking the cross product of those 2 vectors we know the absolute direction
    !> All the other elements in the list will have the right direction because the relative direction within a list is consistent
    !   4-------3                        _
    !   |     ^ |  a x b  has direction |.| (going out of the screen, thus it is ccw, direction =1)
    !   |    /  |                        -
    !   |  b/   |
    !   |  /    |
    !   | / a   |
    !   1------>2
    !     1
    
    a_12(1) = node_list%node(node_2)%x(1,1,1) - node_list%node(node_1)%x(1,1,1)
    a_12(2) = node_list%node(node_2)%x(1,1,2) - node_list%node(node_1)%x(1,1,2)
    a_12(3) = 0 !< element lies in the RZ plane
    b_13(1) = node_list%node(node_3)%x(1,1,1) - node_list%node(node_1)%x(1,1,1)
    b_13(2) = node_list%node(node_3)%x(1,1,2) - node_list%node(node_1)%x(1,1,2)
    b_13(3) = 0
    
    c =  cross_product(a_12, b_13)

    direction(current_element) = int( c(3)/abs(c(3))) !< results in an integer value of 1 or -1. | 1 = ccw | -1 = cw
  end if
end subroutine find_new_edge_domain_start


!> Entry point for routines in this module. Search for the edges of the domain in node_list, element_list
!> and return these in edge_domain(:).
!> sides, x1, x2, flip_ud and discont_corner are optional parameters to influence the routine.
subroutine find_edge_domains(node_list, element_list, edge_domain, sides, x1, x2, flip_ud, discont_corner)
use mod_neighbours

type(type_node_list), intent(in)                                :: node_list
type(type_element_list), intent(inout)                          :: element_list
type(type_edge_domain), dimension(:) , allocatable, intent(out) :: edge_domain

logical, optional, intent(in)              :: sides(4) !< whether to include a boundary at side i = 1..4. Use these to select divertor plates for instance.
real*8, optional, dimension(2), intent(in) :: x1, x2 !< Points defining a line to exclude stuff above/left of
logical, optional, intent(in)              :: flip_ud !< flip_ud false, everything above cutoff line is disregarded, flip_up true -> everything below
!< has no effect if x1 and x2 are not present
logical, optional, intent(in)              :: discont_corner !< If present and .true. create different edge_domains at hard corners (i.e. where side changes)

integer                :: i, j , i_elm, element_i, element_j, i_side, i_list, element_n, lower_list, higher_list, unique_lists, ni, li, list_number
integer                :: i_node1, i_node2, j_node1, j_node2, j_elm, i_dir
logical                :: allowed
real*8, dimension(2)   :: x
logical                :: cut_off, node_found, all_found, use_cutoff_function, my_flip_ud, my_discont_corner
integer                :: n_edge_max, node_1, node_2, node_3, direction, last_vertex
logical, dimension(:), allocatable  :: lists_used_boolean
integer, dimension(:), allocatable  :: lists_used
real*8, dimension(3) :: c, a_12, b_13

type(type_node)    :: node
type(type_element) :: element
type(type_edge_element), dimension(:), allocatable :: edge_element
type(type_edge_element), dimension(:), allocatable :: resize_edge_element !< temporary variable to resize edge_element

my_discont_corner = .false.
if(present(discont_corner)) my_discont_corner = discont_corner


! Precalculate the maximum possible amount of edge elements
n_edge_max = 0
do i = 1, element_list%n_elements 
  n_edge_max = n_edge_max + count(element_list%element(i)%neighbours == 0 , dim=1)
end do

! Check if neighbours have been calculated properly before this routine was called
if (n_edge_max .eq. 4*element_list%n_elements) then
  write(*,*) 'WARNING: all neighbours == 0'
  write(*,*) 'update_neighbours will be executed automatically'
  
  call update_neighbours(node_list,element_list)
  
  ! There is a new neighbour list, redo maximum possible amount of edge elements
  n_edge_max = 0
  do i = 1, element_list%n_elements 
    n_edge_max = n_edge_max + count(element_list%element(i)%neighbours == 0 , dim=1)
  end do
end if

allocate(edge_element(n_edge_max))
edge_element(:)%i_list = 0

!> do loop over n_elements
!> in every element loop over sides and check if neighbours = 0
!> if neighbours == 0
!> save i_elm, i_side, node_1, node_2
!> save this element as this is an edge element
use_cutoff_function = .false.
if (present(x1) .and. present(x2)) use_cutoff_function = .true.

element_i = 1 ! counter that tracks the number of edge elements for the location in the array
! by usage of element_i, elements with 2 boundaries will be selected twice and will form separate boundary elements
jorek_element : do i_elm=1,element_list%n_elements
  check_sides : do i_side = 1,4 ! defined which sides by user (input arg sides) array
    if (element_list%element(i_elm)%neighbours(i_side) .eq. 0) then ! if we are on the edge of the plasma
      node_1 = element_list%element(i_elm)%vertex(i_side)
      node_2 = element_list%element(i_elm)%vertex(modulo(i_side,4)+1)

      if (use_cutoff_function) then
        ! use cutoff function to determine if i_elm location is below the cutoff line
        x      = (node_list%node(node_1)%x(1,1,:) + node_list%node(node_2)%x(1,1,:))*0.5 !< average R and Z position of the 2 nodes on a side    
        my_flip_ud = .false.
        if (present(flip_ud)) my_flip_ud = flip_ud
        cut_off = check_cutoff_edge_domain(x, x1, x2, my_flip_ud) 
        if (cut_off) cycle ! element is not included in edge_element, check next side
      endif
      
      if (present(sides)) then
        if (.not. sides(i_side)) cycle !< skip this element if sides are given and this is not that one
      end if   

      ! Save the found element
      edge_element(element_i)%i_elm  = i_elm
      edge_element(element_i)%i_side = i_side 
      edge_element(element_i)%node_1 = node_1 ! element%vertex for the correct jorek node
      edge_element(element_i)%node_2 = node_2 ! See mod_neighbours for element structure

      element_i = element_i + 1
    endif

  end do check_sides
end do jorek_element
element_n = element_i-1



! reallocate edge_element array to size 1:element_n
! because initially the array is made too large to make sure all the elements fit in there.
allocate(resize_edge_element(element_n), source=edge_element(1:element_n))
call move_alloc(resize_edge_element,edge_element)


!>--------------------------------------------------------------------------------------Sorting part
!> Sorting loop to find adjacent elements and sort them per edge domain
!> ( per side or user specified region)
!> Depending on the cutoff criteria and/or other user specification
!> multiple disconnected domains can be produced (the amount of lists is variable)


!> JOREK elements can be ccw and cw, the direction of the elements has to be tracked to sort them correctly for the edge element arrays
!> Counter clockwise elements (direction = 1)   | ccw element and cw element, cw element has (direction = -1)
!>     i       j                              i       j
!> |  /\  ||  /\  |                       |  /\  ||  /\  |
!> | /  \ || /  \ |                       | /  \ || /  \ |
!> |/  1 \||/  1 \|                       |/  1 \||/  1 \|
!> 1^^^^^^21^^^^^^2                       1^^^^^^22^^^^^^1

!> This loop goes over all the boundary element
!> per target boundary element it will go through all other boundary element to see if it is adjacent
!> this is done by checking if node_2 of the current element is equal to node_1 of another element,
!> where those 2 elements have the same sides on the boundary (if discont_corner .eq. .true. )
!> Then they have the same number in i_list
!> edge_element%i of next adjacent element is increased by 1
!> element_i is target element, element_j is a possible neighbour
all_found                       = .false. !< for check if all edge elements have been listed
last_vertex                     = -1 !< remembers the last vertex to make sure to the loop cannot bounce between 2 elements

call find_new_edge_domain_start(edge_element%i_elm, edge_element%i_list, element_i, edge_element%i,edge_element%i_side, edge_element%direction, all_found, element_list, node_list ) 


sorting_loop : do while (.not. all_found)

  node_found = .false. ! if an adjacent node is found (without an i_list) in the loop, change node_found to true
  check_for_adjacent_node : do element_j = 1, element_n
    i_dir   = edge_element(element_i)%direction
    i_elm   = edge_element(element_i)%i_elm
    i_node1 = edge_element(element_i)%node_1
    i_node2 = edge_element(element_i)%node_2
    j_elm   = edge_element(element_j)%i_elm
    j_node1 = edge_element(element_j)%node_1
    j_node2 = edge_element(element_j)%node_2
    allowed = (i_elm .ne. j_elm .or. .not. my_discont_corner) .and. (element_i .ne. element_j) ! if this transition is allow due to discont corner rules and not self


    ! if the element we are testing is adjacent to the current element (element_i) and has the same orientation (ccw)
    ! (if i_elm .ne. j_elm because we cannot have sharp corners if the element numbers are different. otherwise if .not.
    ! my_discont_corner we don't care)
    if (j_node1 .eq. i_node2 .and. allowed .and. j_node1 .ne. last_vertex .and. i_dir .eq. 1) then ! ccw to ccw
      node_found = .true.
      edge_element(element_j)%direction = edge_element(element_i)%direction
      last_vertex = edge_element(element_i)%node_2
      exit check_for_adjacent_node
       
    else if (j_node2 .eq. i_node2 .and. allowed .and. j_node2 .ne. last_vertex) then ! ccw to cw
      ! why do we not check i_dir?
      node_found = .true.
      edge_element(element_j)%direction = (-1)*edge_element(element_i)%direction
      last_vertex = edge_element(element_i)%node_2
      exit check_for_adjacent_node
      
    else if (j_node1 .eq. i_node1 .and. allowed .and. j_node1 .ne. last_vertex) then ! cw to ccw
      node_found = .true.
      edge_element(element_j)%direction = (-1)*edge_element(element_i)%direction
      last_vertex = edge_element(element_i)%node_1
      exit check_for_adjacent_node
      
    else if (j_node2 .eq. i_node1 .and. allowed .and. j_node2 .ne. last_vertex .and. i_dir .ne. 1) then ! cw to cw
      ! why do we not check i_dir?
      node_found = .true.
      edge_element(element_j)%direction = edge_element(element_i)%direction
      last_vertex = edge_element(element_i)%node_1
      exit check_for_adjacent_node  
    end if
    
  end do check_for_adjacent_node
          
     
  ! Since we have not found a connection to another element for this node,
  ! this element must be the last in our sub domain, and we move on to the next sub domain if it exists
  if (.not. node_found) then
    call find_new_edge_domain_start(edge_element%i_elm, edge_element%i_list, element_i, edge_element%i,edge_element%i_side,edge_element%direction, all_found, element_list, node_list)
    
    cycle !goes to beginning of while loop, so the next part will be skipped
  end if

  ! now connect the new found element to our previous one
  if (edge_element(element_j)%i_list .eq. 0) then !< i_list = 0, element has not been located as a neighbour of another element ('fresh' element)
    edge_element(element_j)%i_list = edge_element(element_i)%i_list ! so place it in the same list
    edge_element(element_j)%i      = edge_element(element_i)%i +1 ! next to its neighbour, in order
    element_i                      = element_j ! update number for the next loop

  else if (edge_element(element_j)%i_list .ne. edge_element(element_i)%i_list) then
    ! another list has been encountered with elements connected to this one
    ! combine the two lists to 1 list with the lowest i_list
    lower_list  = min(edge_element(element_j)%i_list,edge_element(element_i)%i_list)
    higher_list = max(edge_element(element_j)%i_list,edge_element(element_i)%i_list)
    
    where (edge_element%i_list .eq. edge_element(element_j)%i_list) edge_element%i = edge_element%i + edge_element(element_i)%i      
    where (edge_element%i_list .eq. higher_list)                    edge_element%i_list = lower_list      
    
    ! after connecting find a new place to start the search
    call find_new_edge_domain_start(edge_element%i_elm, edge_element%i_list, element_i, edge_element%i,edge_element%i_side,edge_element%direction, all_found, element_list, node_list) 
    
  else if (edge_element(element_j)%i_list .eq. edge_element(element_i)%i_list .and. edge_element(element_j)%i_list .ne. 0) then
    ! the same list as the current one has been encountered
    ! this means we have come full circle. we could renumber the elements or we can not.
    ! let's not.

    ! in the case of disconnected domains we might need to still find edge elements. Otherwise this function just sets all_found to .true.
    call find_new_edge_domain_start(edge_element%i_elm, edge_element%i_list, element_i, edge_element%i,edge_element%i_side, edge_element%direction, all_found, element_list, node_list)
  endif
end do  sorting_loop



! We now have a list of boundary elements with an index attached to it describing the order of elements in a list.
! We need to split this into many lists, one per connected domain (i.e. between 1 and 4 times the number of disconnected domains
! (which is 1 in jorek))

! the loop will figure out which i_list are used for elements
! it will then count how many lists are used. 

allocate(lists_used_boolean(maxval(edge_element%i_list, dim = 1))) ! maximum number of different lists
do i = 1, maxval(edge_element%i_list, dim = 1)
  lists_used_boolean(i) = any(edge_element%i_list .eq. i, dim = 1)
  ! if the list number occurs in i_lists it will get value = true in lists_used
  ! the position of the index is the number of the list that is used
end do
unique_lists = count(lists_used_boolean, dim=1) ! Amount of different lists in use

allocate(edge_domain(unique_lists)) ! edge_domain has the size of the amount of different lists
allocate(lists_used(unique_lists))
j = 1 ! dummy counter to produce right list number
do i = 1, unique_lists
  if (lists_used_boolean(i)) then 
    lists_used(j) = i
    j = j+1
  end if
end do


!> This loop will construct the output edge_domain
do li = 1,unique_lists !li  = list counter  
  list_number = lists_used(li) !< gives the number of i_list used in domain
  
  allocate(edge_domain(li)%i_elm( count(edge_element%i_list .eq. list_number) ))
  allocate(edge_domain(li)%i_side( count(edge_element%i_list .eq. list_number) ))
  allocate(edge_domain(li)%direction( count(edge_element%i_list .eq. list_number) ))
  
  do i = 1,element_n  
    if (edge_element(i)%i_list .eq. list_number) then
      ! go over all elements and check if this is the next number to place in edge_domain
      ! then place the information of that element in the corresponding location of the output array
        
      edge_domain(list_number)%i_elm(edge_element(i)%i)        = edge_element(i)%i_elm
      edge_domain(list_number)%i_side(edge_element(i)%i)       = edge_element(i)%i_side
      edge_domain(list_number)%direction(edge_element(i)%i)    = edge_element(i)%direction       
    end if  
  end do 
end do
!< Now there are unique_lists amount of edge_domains, each with an ordered list of the elements, sides and direction
end subroutine find_edge_domains
end module mod_edge_domain
