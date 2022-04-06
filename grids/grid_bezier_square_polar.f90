!> Create a square grid surrounded by a transition ring and a polar grid
subroutine grid_bezier_square_polar(nR, nZ, n_radial, R_begin, R_end, Z_begin, Z_end, R_geo, Z_geo,&
  amin, fbnd, fpsi, mf, boundary, node_list, element_list)

use constants
use mod_parameters
use data_structure
use mod_neighbours, only: update_neighbours

implicit none

! --- Routine parameters
integer,                 intent(in)    :: nR             !< Number of horizontal nodes (square grid)
integer,                 intent(in)    :: nZ             !< Number of vertical nodes (square grid)
integer,                 intent(in)    :: n_radial       !< Number of radial nodes (polar grid)
real*8,                  intent(in)    :: R_begin        !< R-min (square grid)
real*8,                  intent(in)    :: R_end          !< R-max (square grid)
real*8,                  intent(in)    :: Z_begin        !< Z-min (square grid)
real*8,                  intent(in)    :: Z_end          !< Z-max (square grid)
real*8,                  intent(in)    :: R_geo          !< R-position of grid center
real*8,                  intent(in)    :: Z_geo          !< Z-position of grid center
real*8,                  intent(in)    :: amin           !< minor radius
real*8,                  intent(in)    :: fbnd(*)        !< Fourier series describing the radius as
                                                         !!   function of the poloidal angle
real*8,                  intent(in)    :: fpsi(*)        !< Fourier series of flux at the boundary
integer,                 intent(in)    :: mf             !< Number of Fourier modes in fbnd and fpsi
logical,                 intent(in)    :: boundary       !< (CURRENTLY UNUSED)
type(type_node_list),    intent(inout) :: node_list      !< list of grid nodes
type(type_element_list), intent(inout) :: element_list   !< list of finite elements

! --- Local variables
real*8  :: radius, x_tmp(n_coord_tor, n_degrees,n_dim), index_tmp(n_degrees)
real*8  :: Rgeo, Zgeo, u_length, angle_start
integer :: n_pol, i, k, n_node_start, n_element_start, n_element_polar, n_polar, n_glue, n_sqr
integer :: i_glue, i_polar, n_index_start
real*8, external :: dlength

write(*,*) '**********************************'
write(*,*) '* creating square/polar grid     *'
write(*,*) '**********************************'
write(*,*) '  n_R, n_Z, n_radial : ',nR,nZ,n_radial

!-------------- start with a square central grid

call grid_bezier_square(nR,nZ,R_begin,R_end,Z_begin,Z_end,.false.,node_list,element_list)

write(*,'(A,i6)') ' number of nodes    : ',node_list%n_nodes
write(*,'(A,i6)') ' number of elements : ',element_list%n_elements

!-------------- adapt square grid such that all edge vectors point outwards

do i=1, nR
  node_list%node(i)%x(1,3,2) = - node_list%node(i)%x(1,3,2)
enddo
do i=1,nR-1
  element_list%element(i)%size(1,3) = -  element_list%element(i)%size(1,3)
  element_list%element(i)%size(2,3) = -  element_list%element(i)%size(2,3)
enddo
do i=1, nZ
  node_list%node((i-1)*nR+1)%x(1,2,1) = - node_list%node((i-1)*nR+1)%x(1,2,1)
enddo
do i=1,nZ-1
  element_list%element((i-1)*(nR-1)+1)%size(1,2) = - element_list%element((i-1)*(nR-1)+1)%size(1,2)
  element_list%element((i-1)*(nR-1)+1)%size(4,2) = - element_list%element((i-1)*(nR-1)+1)%size(4,2)
enddo


!-------------- add additional nodes (adventurous)

!--------------------------------------- south
n_node_start = node_list%n_nodes

do i=1, nR
  node_list%node(n_node_start+i) = node_list%node(i)

  x_tmp = node_list%node(n_node_start+i)%x
  node_list%node(n_node_start+i)%x(1,2,:) = x_tmp(1,3,:)
  node_list%node(n_node_start+i)%x(1,3,:) = x_tmp(1,2,:)

  index_tmp = node_list%node(n_node_start+i)%index
  node_list%node(n_node_start+i)%index(2) = index_tmp(3)
  node_list%node(n_node_start+i)%index(3) = index_tmp(2)

enddo
node_list%n_nodes = node_list%n_nodes + nR
n_node_start      = node_list%n_nodes

!--------------------------------------- east
do i=1, nZ
  node_list%node(n_node_start+i) = node_list%node(i*nR)
enddo
node_list%n_nodes = node_list%n_nodes + nZ
n_node_start      = node_list%n_nodes

!--------------------------------------- north
do i=1, nR
  node_list%node(n_node_start+i) = node_list%node(nR*nZ - i + 1)

  x_tmp = node_list%node(n_node_start+i)%x
  node_list%node(n_node_start+i)%x(1,2,:) = x_tmp(1,3,:)
  node_list%node(n_node_start+i)%x(1,3,:) = x_tmp(1,2,:)

  index_tmp = node_list%node(n_node_start+i)%index
  node_list%node(n_node_start+i)%index(2) = index_tmp(3)
  node_list%node(n_node_start+i)%index(3) = index_tmp(2)

enddo
node_list%n_nodes = node_list%n_nodes + nR
n_node_start      = node_list%n_nodes

!--------------------------------------- west
do i=1, nZ
  node_list%node(n_node_start+i) = node_list%node(nR*NZ - i*nR + 1)
enddo
node_list%n_nodes = node_list%n_nodes + nZ
n_node_start      = node_list%n_nodes

n_index_start = 0
do i=1,n_node_start
  n_index_start = max(n_index_start,maxval(node_list%node(i)%index(:)))
enddo

!-------------------------------------- add index for new degrees of freedom at corner nodes

node_list%node(nR*nZ+1)%index(2)          = n_index_start + 1
node_list%node(nR*nZ+nR)%index(2)         = n_index_start + 2
node_list%node(nR*nZ+nR+1)%index(2)       = n_index_start + 2
node_list%node(nR*nZ+nR+nZ)%index(2)      = n_index_start + 3
node_list%node(nR*nZ+nR+nZ+1)%index(2)    = n_index_start + 3
node_list%node(nR*nZ+2*nR+nZ)%index(2)    = n_index_start + 4
node_list%node(nR*nZ+2*nR+nZ+1)%index(2)  = n_index_start + 4
node_list%node(nR*nZ+2*nR+2*nZ)%index(2)  = n_index_start + 1

!do i=nR*nZ+1,node_list%n_nodes
!  write(*,*) i, node_list%node(i)%index
!enddo


element_list%n_elements = element_list%n_elements + 2*(nR-1)+2*(nZ-1)

if ( node_list%n_nodes > n_nodes_max ) then
  write(*,*) 'ERROR in grid_bezier_square_polar: hard-coded parameter n_nodes_max is too small'
  stop
else if ( element_list%n_elements > n_elements_max ) then
  write(*,*) 'ERROR in grid_bezier_square_polar: hard-coded parameter n_elements_max is too small'
  stop
end if

!-------------- add the polar grid

n_pol    = 2*(nR-1) + 2*(nZ-1)

Rgeo = (R_begin + R_end)/2.d0
Zgeo = (Z_begin + Z_end)/2.d0

if ( ( R_geo - Rgeo > 1.d-10 ) .or. ( Z_geo - Zgeo > 1.d-10 ) ) then
  write(*,*) 'ERROR! Input inconsistent. For grid_bezier_square_polar, the following is required:'
  write(*,*) '  R_geo = (R_begin+R_end)/2 = ', Rgeo
  write(*,*) '  Z_geo = (Z_begin+Z_end)/2 = ', Zgeo
  stop
end if

radius = ( 1.d0 + 2.d0/(n_radial-1.d0) ) * sqrt( (R_end-R_begin)**2 + (Z_end-Z_begin)**2 ) / 2.d0
write(*,*) radius, amin

write(*,'(A,i6)') ' number of nodes    : ',node_list%n_nodes
write(*,'(A,i6)') ' number of elements : ',element_list%n_elements

angle_start = -0.75 * PI

call grid_polar_bezier(Rgeo,Zgeo,amin,radius,angle_start,fbnd,fpsi,mf,n_radial,n_pol,node_list,element_list)

!-------------- adapt the corner nodes
!node_list%node(nR*nZ+1)%x(1,1,2)  = 1.1*node_list%node(nR*nZ+1)%x(1,1,2)
node_list%node(nR*nZ+1)%x(1,2,:)  = node_list%node(nR*nZ + 2*nR + 2*Nz + 1)%x(1,2,:)
node_list%node(nR*nZ+1)%x(1,3,:)  = node_list%node(1)%x(1,2,:)

!node_list%node(nR*nZ+nR)%x(1,1,2) = 1.1 * node_list%node(nR*nZ+nR)%x(1,1,2)
node_list%node(nR*nZ+nR)%x(1,2,:) = node_list%node(nR*nZ + 2*nR + 2*Nz + nR)%x(1,2,:)
node_list%node(nR*nZ+nR)%x(1,3,:) = node_list%node(nR)%x(1,2,:)



!-------------- glue the two grids together (south)


n_polar = nR*nZ + 2*nR + 2*Nz         ! starting node of the polar grid
n_glue  = nR*nZ                       ! starting node of the additional nodes

n_sqr   = (nR-1) * (nZ-1)             ! starting element of the additional elements

n_element_start = (nR-1)*(nZ-1)
n_element_polar = (nR-1)*(nZ-1)  + 2*(nR-1) + 2*(nZ-1)

do i=1,nR-1

  element_list%element(n_element_start+i)%vertex(2) = n_polar + i
  element_list%element(n_element_start+i)%vertex(3) = n_polar + i + 1
  element_list%element(n_element_start+i)%vertex(4) = n_glue  + i + 1
  element_list%element(n_element_start+i)%vertex(1) = n_glue  + i

  u_length = sqrt(node_list%node(n_glue+i)%x(1,2,1)**2 + node_list%node(n_glue+i)%x(1,2,2)**2)
  element_list%element(n_element_start+i)%size(1,2) = dlength(node_list%node(n_polar+i)%x(1,1,:), &
                                                              node_list%node(n_glue+i)%x(1,1,:))/3.d0 / u_length
  element_list%element(n_element_start+i)%size(1,3) = element_list%element(i)%size(1,2)

  u_length = sqrt(node_list%node(n_polar+i)%x(1,2,1)**2 + node_list%node(n_polar+i)%x(1,2,2)**2)
  element_list%element(n_element_start+i)%size(2,2) = - dlength(node_list%node(n_polar+i)%x(1,1,:), &
                                                                node_list%node(n_glue+i)%x(1,1,:))/3.d0 / u_length
  element_list%element(n_element_start+i)%size(2,3) = element_list%element(n_element_polar+i)%size(1,3)

  u_length = sqrt(node_list%node(n_polar+i+1)%x(1,2,1)**2 + node_list%node(n_polar+i+1)%x(1,2,2)**2)
  element_list%element(n_element_start+i)%size(3,2) = - dlength(node_list%node(n_polar+i+1)%x(1,1,:), &
                                                                node_list%node(n_glue+i+1)%x(1,1,:))/3.d0 / u_length
  element_list%element(n_element_start+i)%size(3,3) = element_list%element(n_element_polar+i)%size(4,3)

  u_length = sqrt(node_list%node(n_glue+i+1)%x(1,2,1)**2 + node_list%node(n_glue+i+1)%x(1,2,2)**2)
  element_list%element(n_element_start+i)%size(4,2) = dlength(node_list%node(n_polar+i+1)%x(1,1,:), &
                                                              node_list%node(n_glue+i+1)%x(1,1,:))/3.d0 / u_length
  element_list%element(n_element_start+i)%size(4,3) = element_list%element(i)%size(3,2)

enddo


!-------------- adapt the corner nodes (east)

n_element_start = (nR-1)*(nZ-1) + (nR-1)

!node_list%node(nR*nZ+nR+1)%x(1,1,1)  = 1.1*node_list%node(nR*nZ+nR+1)%x(1,1,1)
node_list%node(nR*nZ+nR+1)%x(1,2,:)  = node_list%node(nR*nZ+nR)%x(1,2,:)
node_list%node(nR*nZ+nR+1)%x(1,3,:)  = node_list%node(nR)%x(1,3,:)

!node_list%node(nR*nZ+nR+nZ)%x(1,1,1) = 1.1*node_list%node(nR*nZ+nR+nZ)%x(1,1,1)
node_list%node(nR*nZ+nR+nZ)%x(1,2,:) = node_list%node(nR*nZ + 2*nR + 2*Nz + nR + nZ - 1)%x(1,2,:)
node_list%node(nR*nZ+nR+nZ)%x(1,3,:) = node_list%node(nR*nZ)%x(1,3,:)

do i=1,nZ-1

  element_list%element(n_element_start+i)%vertex(2) = n_polar + (nR-1) + i
  element_list%element(n_element_start+i)%vertex(3) = n_polar + (nR-1) + i + 1
  element_list%element(n_element_start+i)%vertex(4) = n_glue  + nR     + i + 1
  element_list%element(n_element_start+i)%vertex(1) = n_glue  + nR     + i

  u_length = sqrt(node_list%node(n_glue+nR+i)%x(1,2,1)**2 + node_list%node(n_glue+nR+i)%x(1,2,2)**2)
  element_list%element(n_element_start+i)%size(1,2) = dlength(node_list%node(n_polar+nR+i-1)%x(1,1,:), &
                                                              node_list%node(n_glue+nR+i)%x(1,1,:))/3.d0 / u_length
  element_list%element(n_element_start+i)%size(1,3) = element_list%element(i)%size(1,2)

  u_length = sqrt(node_list%node(n_polar+nR+i-1)%x(1,2,1)**2 + node_list%node(n_polar+nR+i-1)%x(1,2,2)**2)
  element_list%element(n_element_start+i)%size(2,2) = - dlength(node_list%node(n_polar+nR+i-1)%x(1,1,:), &
                                                                node_list%node(n_glue+nR+i)%x(1,1,:))/3.d0 / u_length
  element_list%element(n_element_start+i)%size(2,3) = element_list%element(n_element_polar+nR+i-1)%size(1,3)

  u_length = sqrt(node_list%node(n_polar+nR+i)%x(1,2,1)**2 + node_list%node(n_polar+nR+i)%x(1,2,2)**2)
  element_list%element(n_element_start+i)%size(3,2) = - dlength(node_list%node(n_polar+nR+i)%x(1,1,:), &
                                                                node_list%node(n_glue+i+nR+1)%x(1,1,:))/3.d0 / u_length
  element_list%element(n_element_start+i)%size(3,3) = element_list%element(n_element_polar+nR+i-1)%size(4,3)

  u_length = sqrt(node_list%node(n_glue+nR+i+1)%x(1,2,1)**2 + node_list%node(n_glue+nR+i+1)%x(1,2,2)**2)
  element_list%element(n_element_start+i)%size(4,2) = dlength(node_list%node(n_polar+nR+i)%x(1,1,:), &
                                                              node_list%node(n_glue+nR+i+1)%x(1,1,:))/3.d0 / u_length
  element_list%element(n_element_start+i)%size(4,3) = element_list%element(i)%size(3,2)

enddo

!-------------- adapt the corner nodes (north)

!node_list%node(nR*nZ+nR+nZ+1)%x(1,1,2)  = 1.1*node_list%node(nR*nZ+nR+nZ+1)%x(1,1,2)
node_list%node(nR*nZ+nR+nZ+1)%x(1,2,:)  = node_list%node(nR*nZ+nR+nZ)%x(1,2,:)
node_list%node(nR*nZ+nR+nZ+1)%x(1,3,:)  = node_list%node(nR*nZ)%x(1,2,:)

!node_list%node(nR*nZ+2*nR+nZ)%x(1,1,2) = 1.1*node_list%node(nR*nZ+2*nR+nZ)%x(1,1,2)
node_list%node(nR*nZ+2*nR+nZ)%x(1,2,:) = node_list%node(nR*nZ + 2*nR + 2*Nz + 2*nR + nZ - 2)%x(1,2,:)
node_list%node(nR*nZ+2*nR+nZ)%x(1,3,:) = node_list%node(nR*nZ - nR + 1)%x(1,2,:)

n_element_start = (nR-1)*(nZ-1) + (nR-1) + (nZ-1)

do i=1,nR-1

  element_list%element(n_element_start+i)%vertex(2) = n_polar + (nR-1) + (nZ-1) + i
  element_list%element(n_element_start+i)%vertex(3) = n_polar + (nR-1) + (nZ-1) + i + 1
  element_list%element(n_element_start+i)%vertex(4) = n_glue  + nR + nZ + i + 1
  element_list%element(n_element_start+i)%vertex(1) = n_glue  + nR + nZ + i

  u_length = sqrt(node_list%node(n_glue+nR+nZ+i)%x(1,2,1)**2 + node_list%node(n_glue+nR+nZ+i)%x(1,2,2)**2)
  element_list%element(n_element_start+i)%size(1,2) = dlength(node_list%node(n_polar+nR+nZ+i-2)%x(1,1,:), &
                                                              node_list%node(n_glue+nR+nZ+i)%x(1,1,:))/3.d0 / u_length
  element_list%element(n_element_start+i)%size(1,3) = element_list%element((nR-1)*(nZ-1)-i+1)%size(3,2)

  u_length = sqrt(node_list%node(n_polar+nR+nZ+i-2)%x(1,2,1)**2 + node_list%node(n_polar+nR+nZ+i-2)%x(1,2,2)**2)
  element_list%element(n_element_start+i)%size(2,2) = - dlength(node_list%node(n_polar+nR+nZ+i-2)%x(1,1,:), &
                                                                node_list%node(n_glue+nR+nZ+i)%x(1,1,:))/3.d0 / u_length
  element_list%element(n_element_start+i)%size(2,3) = element_list%element(n_element_polar+nR+nZ+i-2)%size(1,3)

  u_length = sqrt(node_list%node(n_polar+nR+nZ+i-1)%x(1,2,1)**2 + node_list%node(n_polar+nR+nZ+i-1)%x(1,2,2)**2)
  element_list%element(n_element_start+i)%size(3,2) = - dlength(node_list%node(n_polar+nR+nZ+i-1)%x(1,1,:), &
                                                                node_list%node(n_glue+i+nR+nZ+1)%x(1,1,:))/3.d0 / u_length
  element_list%element(n_element_start+i)%size(3,3) = element_list%element(n_element_polar+nR+nZ+i-2)%size(4,3)

  u_length = sqrt(node_list%node(n_glue+nR+nZ+i+1)%x(1,2,1)**2 + node_list%node(n_glue+nR+nZ+i+1)%x(1,2,2)**2)
  element_list%element(n_element_start+i)%size(4,2) = dlength(node_list%node(n_polar+nR+nZ+i-1)%x(1,1,:), &
                                                              node_list%node(n_glue+nR+nZ+i+1)%x(1,1,:))/3.d0 / u_length
  element_list%element(n_element_start+i)%size(4,3) = - element_list%element((nR-1)*(nZ-1)-i)%size(3,2)

  if (i .eq. (nR-1)) element_list%element(n_element_start+i)%size(4,3) = - element_list%element(n_element_start+i)%size(4,3)

enddo

!-------------- adapt the corner nodes (west)

!node_list%node(nR*nZ+2*nR+nZ+1)%x(1,1,1)  = 1.1*node_list%node(nR*nZ+2*nR+nZ+1)%x(1,1,1)
node_list%node(nR*nZ+2*nR+nZ+1)%x(1,2,:)  = node_list%node(nR*nZ+2*nR+nZ)%x(1,2,:)
node_list%node(nR*nZ+2*nR+nZ+1)%x(1,3,:)  = node_list%node(nR*nZ-nR+1)%x(1,3,:)

!node_list%node(nR*nZ+2*nR+2*nZ)%x(1,1,1) = 1.1*node_list%node(nR*nZ+2*nR+2*nZ)%x(1,1,1)
node_list%node(nR*nZ+2*nR+2*nZ)%x(1,2,:) = node_list%node(nR*nZ+1)%x(1,2,:)
node_list%node(nR*nZ+2*nR+2*nZ)%x(1,3,:) = node_list%node(1)%x(1,3,:)

n_element_start = (nR-1)*(nZ-1) + 2*(nR-1) + (nZ-1)

do i=1,nZ-1

  element_list%element(n_element_start+i)%vertex(2) = n_polar + 2*(nR-1) + (nZ-1) + i
  element_list%element(n_element_start+i)%vertex(3) = n_polar + 2*(nR-1) + (nZ-1) + i + 1
  element_list%element(n_element_start+i)%vertex(4) = n_glue  + 2*nR     + nZ     + i + 1
  element_list%element(n_element_start+i)%vertex(1) = n_glue  + 2*nR     + nZ     + i

  if (i .eq. (nZ-1)) element_list%element(n_element_start+i)%vertex(3) = n_polar + 1

  i_glue  = n_glue + nR + nZ + nR + i
  i_polar = n_polar+nR+nZ+nR+i-3

  u_length = sqrt(node_list%node(i_glue)%x(1,2,1)**2 + node_list%node(i_glue)%x(1,2,2)**2)
  element_list%element(n_element_start+i)%size(1,2) = dlength(node_list%node(i_polar)%x(1,1,:), &
                                                              node_list%node(i_glue)%x(1,1,:))/3.d0 / u_length
  element_list%element(n_element_start+i)%size(1,3) = element_list%element((nR-1)*(nZ-1)-i*(nR-1)+1)%size(4,3)

  u_length = sqrt(node_list%node(i_polar)%x(1,2,1)**2 + node_list%node(i_polar)%x(1,2,2)**2)
  element_list%element(n_element_start+i)%size(2,2) = - dlength(node_list%node(i_polar)%x(1,1,:), &
                                                                node_list%node(i_glue)%x(1,1,:))/3.d0 / u_length
  element_list%element(n_element_start+i)%size(2,3) = element_list%element(n_element_polar+nR+nZ+i-2)%size(1,3)

  u_length = sqrt(node_list%node(i_polar+1)%x(1,2,1)**2 + node_list%node(i_polar+1)%x(1,2,2)**2)
  element_list%element(n_element_start+i)%size(3,2) = - dlength(node_list%node(i_polar+1)%x(1,1,:), &
                                                                node_list%node(i_glue+1)%x(1,1,:))/3.d0 / u_length
  element_list%element(n_element_start+i)%size(3,3) = element_list%element(n_element_polar+nR+nZ+i-2)%size(4,3)

  u_length = sqrt(node_list%node(i_glue+1)%x(1,2,1)**2 + node_list%node(i_glue+1)%x(1,2,2)**2)
  element_list%element(n_element_start+i)%size(4,2) = dlength(node_list%node(i_polar+1)%x(1,1,:), &
                                                              node_list%node(i_glue+1)%x(1,1,:))/3.d0 / u_length
  element_list%element(n_element_start+i)%size(4,3) = element_list%element((nR-1)*(nZ-1)-i*(nZ-1)+1)%size(1,3)

  if (i .eq. (nZ-1)) then

    u_length = sqrt(node_list%node(n_polar+1)%x(1,2,1)**2 + node_list%node(n_polar+1)%x(1,2,2)**2)
    element_list%element(n_element_start+i)%size(3,2) = - dlength(node_list%node(n_polar+1)%x(1,1,:), &
                                                                  node_list%node(i_glue+1)%x(1,1,:))/3.d0 / u_length

    u_length = sqrt(node_list%node(i_glue+1)%x(1,2,1)**2 + node_list%node(i_glue+1)%x(1,2,2)**2)
    element_list%element(n_element_start+i)%size(4,2) = dlength(node_list%node(n_polar+1)%x(1,1,:), &
                                                                node_list%node(i_glue+1)%x(1,1,:))/3.d0 / u_length

  endif

enddo

do i=(nR-1)*(nZ-1)+1,(nR-1)*(nZ-1) + 2*(nR-1)+2*(nZ-1)

  do k=1,n_degrees
    element_list%element(i)%size(k,1) = 1.
    element_list%element(i)%size(k,4) = element_list%element(i)%size(k,2) * element_list%element(i)%size(k,3)
  enddo

enddo

!------------------------------- initialise square center
do i=1,n_polar
  node_list%node(i)%values = 0.d0
  node_list%node(i)%values(1,1,1) = node_list%node(n_polar+1)%values(1,1,1)
enddo

do i=1,node_list%n_nodes
  node_list%node(i)%axis_node = .false.
enddo

call update_neighbours(node_list,element_list, force_rtree_initialize=.true.)
return
end subroutine grid_bezier_square_polar
