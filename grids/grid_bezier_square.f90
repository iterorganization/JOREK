!> Create a square grid based on the Bezier/modified cubic Hermite representation
subroutine grid_bezier_square(nR,nZ,R_begin,R_end,Z_begin,Z_end,boundary,node_list,element_list)

use mod_parameters
use data_structure
use mod_neighbours, only: update_neighbours
use phys_module, only: n_radial, n_flux, XR_r, XR_z, SIG_r, SIG_z, bgf_r, bgf_z, rect_grid_vac_psi

implicit none

! --- Routine parameters
integer,                 intent(in)    :: nR             !< Number of horizontal nodes (square grid)
integer,                 intent(in)    :: nZ             !< Number of vertical nodes (square grid)
real*8,                  intent(in)    :: R_begin        !< R-min (square grid)
real*8,                  intent(in)    :: R_end          !< R-max (square grid)
real*8,                  intent(in)    :: Z_begin        !< Z-min (square grid)
real*8,                  intent(in)    :: Z_end          !< Z-max (square grid)
logical,                 intent(in)    :: boundary       !< Fill boundary information?
type(type_node_list),    intent(inout) :: node_list      !< list of grid nodes
type(type_element_list), intent(inout) :: element_list   !< list of finite elements

! --- Local variables
integer                  :: m, n, i, j, k, iv, inode_0, inode_p, inode, n_elements, i_element, ip, iuv
integer                  :: n_element_start, n_node_start, n_index_start
real*8                   :: xx_0(n_dim),xx_p(n_dim),uv_0(n_dim),uv_p(n_dim)
real*8, external         :: dlength
real*8, external         :: ddot
real*8, allocatable      :: s_tmp(:), t_tmp(:)
logical                  :: grid_accumulation

write(*,*) '*************************************'
write(*,*) '*       grid_bezier_square          *'
write(*,*) '*************************************'
write(*,*) ' creating a square grid :  n_R, n_Z = ',nR,nZ

n_element_start = element_list%n_elements
n_node_start    = node_list%n_nodes

n_index_start = 0
do i=1,n_node_start
  n_index_start = max(n_index_start,maxval(node_list%node(i)%index(:)))
enddo

write(*,*) ' existing no. of elements : ',n_element_start
write(*,*) ' existing number of nodes : ',n_node_start
write(*,*) ' index_start              : ',n_index_start

grid_accumulation = .false.
if ( (n_radial .eq. 0) .and. (n_flux .eq. 0) ) then
  if (      (XR_r(1) .lt. 1.1) .or. (XR_r(2) .lt. 1.1) &
       .or. (XR_z(1) .lt. 1.1) .or. (XR_z(2) .lt. 1.1)  ) then
    grid_accumulation = .true.
  endif
endif

if (grid_accumulation) then
  call tr_allocate(s_tmp,1,nR+1,"s_tmp",CAT_GRID)
  call tr_allocate(t_tmp,1,nZ+1,"t_tmp",CAT_GRID)
  s_tmp = 0
  t_tmp = 0
  call meshac2(nR+1,s_tmp,XR_r(1),XR_r(2),SIG_r(1),SIG_r(2),bgf_r,1.0d0)
  call meshac2(nZ+1,t_tmp,XR_z(1),XR_z(2),SIG_z(1),SIG_z(2),bgf_z,1.0d0)
endif

inode = 0
do j=1,nZ
 do i=1,nR
    inode = inode + 1
    node_list%node(inode)%boundary = 0
 enddo
enddo

inode = 0
do j=1,nZ
  do i=1,nR

    inode = inode + 1

    ! --- For backward compatibility
    if (grid_accumulation) then
      node_list%node(inode)%x(1,1,1) = R_begin + (R_end - R_begin) * s_tmp(i+1)   ! the position of the node
      node_list%node(inode)%x(1,1,2) = Z_begin + (Z_end - Z_begin) * t_tmp(j+1)
    else
      node_list%node(inode)%x(1,1,1) = R_begin + (R_end - R_begin) * float(i-1)/float(nR-1)   ! the position of the node
      node_list%node(inode)%x(1,1,2) = Z_begin + (Z_end - Z_begin) * float(j-1)/float(nZ-1)
    endif

    node_list%node(inode)%x(1,2,1) = 1.0d0                ! the unit vector u
    node_list%node(inode)%x(1,2,2) = sqrt(1.d0 - node_list%node(inode)%x(1,2,1)**2)

    node_list%node(inode)%x(1,3,1) = 0.d0                ! the unit vector v
    node_list%node(inode)%x(1,3,2) = sqrt(1.d0 - node_list%node(inode)%x(1,3,1)**2)

    node_list%node(inode)%x(1,4,1) = 0.                ! the vector w
    node_list%node(inode)%x(1,4,2) = 0.

    if (boundary) then
      if ((i .eq. 1) .or. (i .eq. nR)) node_list%node(inode)%boundary = node_list%node(inode)%boundary + 2
      if ((j .eq. 1) .or. (j .eq. nZ)) node_list%node(inode)%boundary = node_list%node(inode)%boundary + 1
      ! --- Psi vacuum boundary conditions
      if ((i .eq. 1) .or. (i .eq. nR)) node_list%node(inode)%values(1,1,1) = rect_grid_vac_psi * (node_list%node(inode)%x(1,1,1))**2
      if ((j .eq. 1) .or. (j .eq. nZ)) node_list%node(inode)%values(1,1,1) = rect_grid_vac_psi * (node_list%node(inode)%x(1,1,1))**2
    endif
    ! --- Psi vacuum initial conditions
    node_list%node(inode)%values(1,1,1) = rect_grid_vac_psi * (node_list%node(inode)%x(1,1,1))**2

    do k=1, n_degrees
      node_list%node(inode)%index(k) = n_degrees*(inode-1) + k
    enddo

  enddo
enddo

if (grid_accumulation) then
  call tr_deallocate(s_tmp,"s_tmp",CAT_GRID)
  call tr_deallocate(t_tmp,"t_tmp",CAT_GRID)
endif

node_list%n_nodes = nR*nZ
do i=1,node_list%n_nodes
  node_list%node(i)%axis_node = .false.
enddo

n_elements = (nR-1)*(nZ-1)

if ( node_list%n_nodes > n_nodes_max ) then
  write(*,*) 'ERROR in grid_bezier_square: hard-coded parameter n_nodes_max is too small'
  stop
else if ( n_elements > n_elements_max ) then
  write(*,*) 'ERROR in grid_bezier_square: hard-coded parameter n_elements_max is too small'
  stop
end if

i_element  = 0
do n=1,nZ-1
  do m=1,nR-1                        ! define connectivity of finite element mesh
    i_element = i_element + 1
    element_list%element(i_element)%vertex(1) = (n-1)*nR + m
    element_list%element(i_element)%vertex(2) = (n-1)*nR + m + 1
    element_list%element(i_element)%vertex(3) = (n  )*nR + m + 1
    element_list%element(i_element)%vertex(4) = (n  )*nR + m
  enddo
enddo


if (i_element .ne. n_elements) write(*,*) ' something is wrong ! (1)'

element_list%n_elements = n_elements

do k=1, element_list%n_elements   ! fill in the size of the elements


 do iv=1,4                       ! loop over the vertices

   iuv = mod(iv+1,2)+1           ! the direction vector corresponding to this edge (i)

   inode_0 = element_list%element(k)%vertex(iv)
   xx_0    = node_list%node(inode_0)%x(1,1,:)
   uv_0    = node_list%node(inode_0)%x(1,iuv+1,:)

   ip      = mod(iv,4)+1
   inode_p = element_list%element(k)%vertex(ip)
   xx_p    = node_list%node(inode_p)%x(1,1,:)
   uv_p    = node_list%node(inode_p)%x(1,iuv+1,:)

   element_list%element(k)%size(iv,1)     = 1.
   element_list%element(k)%size(iv,iuv+1) = sign(dlength(xx_p,xx_0),ddot(n_dim,xx_p - xx_0,1,uv_0,1)) /3.d0
   element_list%element(k)%size(ip,iuv+1) = sign(dlength(xx_p,xx_0),ddot(n_dim,xx_0 - xx_p,1,uv_p,1)) /3.d0

!    write(*,*) element_list%element(k)%size(iv,iuv+1),element_list%element(k)%size(ip,iuv+1)

 enddo

 do iv=1,4
   element_list%element(k)%size(iv,4) = element_list%element(k)%size(iv,2) * element_list%element(k)%size(iv,3)
 enddo

enddo

call update_neighbours(node_list,element_list, force_rtree_initialize=.true.)
return
end subroutine grid_bezier_square
