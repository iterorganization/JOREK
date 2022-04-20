!> Plot the grid of finite elements with the correct curved edges
subroutine plot_grid(node_list, element_list, boundary_list, bnd_node_list, frame, bezier, &
  gridname)

use tr_module 
use mod_parameters
use data_structure
use phys_module, only: write_ps

implicit none

! --- Routine parameters
type (type_node_list),        intent(in)   :: node_list       !< List of grid nodes
type (type_element_list),     intent(in)   :: element_list    !< List of grid elements
type (type_bnd_element_list), intent(in)   :: boundary_list   !< List of boundary elements
type (type_bnd_node_list),    intent(in)   :: bnd_node_list   !< List of boundary nodes
logical,                      intent(in)   :: frame
logical,                      intent(in)   :: bezier
character(len=*),             intent(in)   :: gridname        !< Used for the name of the ascii output file

! --- Local variables
real*8,allocatable :: xp(:,:)
real*8             :: xs(2,n_dim), xx_0(n_dim), uu_0(n_dim), vv_0(n_dim), ww_0(n_dim)
real*8             :: uv_0(n_dim), uv_p(n_dim), xx_p(n_dim) ,xb(n_degrees,n_dim)
real*8             :: xmax, xmin, ymax, ymin, s, huv_0, huv_p, x_length
integer            :: i, j, inode_0, inode_p, iplot, k, ip, np, iuv, idir_0, idir_p
character*3        :: label

write(*,*) '************************************'
write(*,*) '*           plot_grid              *'
write(*,*) '************************************'
write(*,*) ' number of elements          : ',element_list%n_elements
write(*,*) ' number of boundary elements : ',boundary_list%n_bnd_elements
write(*,*) ' number of nodes             : ',node_list%n_nodes
write(*,*) ' bezier                      : ',bezier

call tr_allocate(xp,1,node_list%n_nodes,1,n_dim,"xp",CAT_GRID)

do i=1,node_list%n_nodes
 xp(i,1:n_dim) = node_list%node(i)%x(1,1,1:n_dim)
enddo
xmax = 1.1 * maxval(xp(:,1))
xmin = minval(xp(:,1))
ymax = 1.1 * maxval(xp(:,2))
ymin = minval(xp(:,2))

xmin  = 1.1 * xmin - 0.1 * xmax
ymin  = 1.1 * ymin - 0.1 * ymax

if ( write_ps ) then
  if (frame) call nframe(21,11,1,xmin,xmax,ymin,ymax,'Bezier grid',11,'X',1,'Y',1)
  call lplot(21,11,461,xp(:,1),xp(:,2),-node_list%n_nodes,1,'Nodes',5,'X',1,'Y',1)

!---------------------------------------------- plot unit vectors

  call lincol(2)
  do i=1,node_list%n_nodes

    xs(1,1:n_dim) = node_list%node(i)%x(1,1,1:n_dim)

    if ((node_list%node(i)%boundary .eq. 1) .or. (node_list%node(i)%boundary .eq. 3)) then

      !    write(*,'(A,i5,6e16.8)') ' boundary point (1,3) : ',i,node_list%node(i)%x(1,1,:)

      x_length = sqrt( node_list%node(i)%x(1,2,1)**2 +  node_list%node(i)%x(1,2,2)**2)

      if (x_length .gt. 0.d0) then
        xs(2,1:n_dim) = xs(1,1:n_dim) + 0.1 * node_list%node(i)%x(1,2,1:n_dim) / x_length
        call lplot6(1,1,xs(:,1),xs(:,2),-2,' ')
      endif

    endif

  enddo


  call lincol(3)
  do i=1,node_list%n_nodes

    if ((node_list%node(i)%boundary .eq. 2) .or. (node_list%node(i)%boundary .eq. 3)) then

      !    write(*,'(A,i5,6e16.8)') ' boundary point (2,3) : ',i,node_list%node(i)%x(1,1,:)

      x_length = sqrt( node_list%node(i)%x(1,3,1)**2 +  node_list%node(i)%x(1,3,2)**2)

      if (x_length .gt. 0.d0) then
        xs(1,1:n_dim) = node_list%node(i)%x(1,1,1:n_dim)
        xs(2,1:n_dim) = xs(1,1:n_dim) + 0.1 * node_list%node(i)%x(1,3,1:n_dim) / x_length
        call lplot6(1,1,xs(:,1),xs(:,2),-2,' ')
      endif

    endif

  enddo
  call lincol(0)
endif
!return

!------------------------------ plot the curved boundaries
call tr_deallocate(xp,"xp",CAT_GRID)


np = 11

call tr_allocate(xp,1,np,1,n_dim,"xp",CAT_GRID)

iplot = 0

open(21, file='grid_'//trim(adjustl(gridname))//'.dat') ! Open the ascii file

do k=1, element_list%n_elements

 do i=1,4                         ! over the 4 edges

   iuv = mod(i+1,2)+1             ! the direction vector corresponding to this edge (i)

   inode_0 = element_list%element(k)%vertex(i)
   xx_0    = node_list%node(inode_0)%x(1,1,:)
   uv_0    = node_list%node(inode_0)%x(1,iuv+1,:)
   huv_0   = element_list%element(k)%size(i,iuv+1)

   ip      = mod(i,4)+1
   inode_p = element_list%element(k)%vertex(ip)
   xx_p    = node_list%node(inode_p)%x(1,1,:)
   uv_p    = node_list%node(inode_p)%x(1,iuv+1,:)
   huv_p   = element_list%element(k)%size(ip,iuv+1)

!    write(*,*) xx_0
!    write(*,*) xx_0 + uv_0*huv_0
!    write(*,*) xx_p + uv_p*huv_p
!    write(*,*) xx_p

   do j=1,np
     s = (float(j-1)/float(np-1))
     xb      = 0.d0
     xb(1,:) = xx_0
     xb(2,:) = xx_0+uv_0*huv_0
     xb(3,:) = xx_p+uv_p*huv_p
     xb(4,:) = xx_p

     call bezier_1d(n_dim, s, xb, xp(j,:))
     write(21,'(2ES16.8)') xp(j,:) ! Output the grid information to the ascii file
   enddo
   write(21,*)
   write(21,*)
   
   if ( write_ps ) call lplot6(1,1,xp(1:np,1),xp(1:np,2),-np,' ')
 enddo
enddo

close(21) ! Close the ascii file
if ( write_ps ) then
  call lincol(1)
  call lplot(11,11,461,xp(:,1),xp(:,2),-iplot,1,' ',1,' ',1,' ',1)
  
  call lincol(4)
  
  do k=1, boundary_list%n_bnd_elements

    idir_0 = boundary_list%bnd_element(k)%direction(1,2)

    inode_0 = boundary_list%bnd_element(k)%vertex(1)
    xx_0    = node_list%node(inode_0)%x(1,1,:)
    uv_0    = node_list%node(inode_0)%x(1,idir_0,:)
    huv_0   = boundary_list%bnd_element(k)%size(1,2)

    idir_p = boundary_list%bnd_element(k)%direction(2,2)

    inode_p = boundary_list%bnd_element(k)%vertex(2)
    xx_p    = node_list%node(inode_p)%x(1,1,:)
    uv_p    = node_list%node(inode_p)%x(1,idir_p,:)
    huv_p   = boundary_list%bnd_element(k)%size(2,2)

!    write(*,*) xx_0
!    write(*,*) xx_0 + uv_0*huv_0
!    write(*,*) xx_p + uv_p*huv_p
!    write(*,*) xx_p

    do j=1,np
      s = (float(j-1)/float(np-1))
      xb      = 0.d0
      xb(1,:) = xx_0
      xb(2,:) = xx_0+uv_0*huv_0
      xb(3,:) = xx_p+uv_p*huv_p
      xb(4,:) = xx_p

      call bezier_1d(n_dim, s, xb, xp(j,:))
    enddo

    call lplot6(1,1,xp(1:np,1),xp(1:np,2),-np,' ')
   
  enddo

  call lincol(0)

  call tr_deallocate(xp,"xp",CAT_GRID)

  if (.not. bezier) return

  !------------------------------ plot the Bezier points
  call tr_allocate(xp,1,12*element_list%n_elements,1,n_dim,"xp",CAT_GRID)

  iplot = 0

  do k=1, element_list%n_elements

    do i=1,4

      inode_0 = element_list%element(k)%vertex(i)
      xx_0    = node_list%node(inode_0)%x(1,1,:)
      uu_0    = node_list%node(inode_0)%x(1,2,:)
      vv_0    = node_list%node(inode_0)%x(1,3,:)
      ww_0    = node_list%node(inode_0)%x(1,4,:)

      xp(iplot+1,:) = node_list%node(inode_0)%x(1,1,:) + uu_0(:) * element_list%element(k)%size(i,2)
      xp(iplot+2,:) = node_list%node(inode_0)%x(1,1,:) + vv_0(:) * element_list%element(k)%size(i,3)
      xp(iplot+3,:) = node_list%node(inode_0)%x(1,1,:)  &
          + uu_0(:) * element_list%element(k)%size(i,2) &
          + vv_0(:) * element_list%element(k)%size(i,3) &
          + ww_0(:) * element_list%element(k)%size(i,4)

      iplot = iplot + 3
    enddo
  enddo

  call lincol(1)
  !call lplot(11,11,461,xp(:,1),xp(:,2),-iplot,1,' ',1,' ',1,' ',1)
  call pplot(11,11,xp(:,1),xp(:,2),iplot,1)
endif
call tr_deallocate(xp,"xp",CAT_GRID)


return
end
