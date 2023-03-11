!> Create a square grid based on the Bezier/modified cubic Hermite representation
subroutine grid_inside_wall(n_R,n_Z,R_begin,R_end,Z_begin,Z_end,boundary,node_list,element_list)

  use mod_parameters
  use data_structure
  use mod_export_restart
  use phys_module, only: xshift, n_limiter, R_limiter, Z_limiter, tokamak_device, manipulate_psi_map
  use grid_xpoint_data, only: n_wall, R_wall, Z_wall
  use mod_eqdsk_tools
  use mod_grid_conversions

  implicit none
  
  ! --- Routine parameters
  integer,                 intent(in)    :: n_R             !< Number of horizontal nodes (square grid)
  integer,                 intent(in)    :: n_Z             !< Number of vertical nodes (square grid)
  real*8,                  intent(in)    :: R_begin        !< R-min (square grid)
  real*8,                  intent(in)    :: R_end          !< R-max (square grid)
  real*8,                  intent(in)    :: Z_begin        !< Z-min (square grid)
  real*8,                  intent(in)    :: Z_end          !< Z-max (square grid)
  logical,                 intent(in)    :: boundary       !< Fill boundary information?
  type(type_node_list),    intent(inout) :: node_list      !< list of grid nodes
  type(type_element_list), intent(inout) :: element_list   !< list of finite elements
  
  ! --- Local variables
  integer                  :: m, n, i, j, k, jp1, ii, jj, iv, inode_0, inode_p, inode, n_elements, i_element, ip, iuv, n_loop, index
  integer                  :: n_element_start, n_node_start, n_index_start
  real*8                   :: xx_0(n_dim),xx_p(n_dim),uv_0(n_dim),uv_p(n_dim)
  real*8, external         :: dlength, ddot
  
  ! --- eqdsk variables
  integer          :: nR_eqdsk, nZ_eqdsk, ier
  real,allocatable :: R_eqdsk(:),Z_eqdsk(:),psi_eqdsk(:,:)
  logical          :: normal_eqdsk, normal_eqdsk_wall
  
  ! --- Grid variables
  integer :: i_elm, i_node
  integer :: nR, nZ, n_elm
  integer :: n_tmp, n_tmp2
  real*8  :: psi, psi_R, psi_Z
  real*8  :: R_tmp(4), Z_tmp(4)
  real*8  :: width
  real*8  :: R_elm, Z_elm, diff
  real*8  :: Z_min, zmin, Z_max, zmax
  integer :: ii_n(4), jj_n(4)
  integer, allocatable :: nR_grid(:,:), node_index(:,:,:), elm_node_index(:,:)
  real*8,  allocatable :: R_grid(:,:), Z_grid(:,:), Zlines(:)
  logical, parameter  :: plot_grid = .true.
  real*8  :: amp, Rm, Zm, dRm, dZm, dPsi, RR, ZZ
  
  
  write(*,*) '*************************************'
  write(*,*) '*       grid_inside_wall            *'
  write(*,*) '*************************************'
  write(*,*) ' creating a rectangular grid inside wall:  n_R, n_Z = ',n_R,n_Z
  

  ! -------------------------------------
  ! --- Get equilibrium from eqdsk file
  ! -------------------------------------
  
  call get_eqdsk_style(normal_eqdsk, normal_eqdsk_wall, ier)
  if (ier .ne. 0) then
    write(*,*)'Problem reading eqdsk file',eqdsk_filename
    write(*,*)'Aborting...'
    stop
  endif
  call get_eqdsk_dimensions(normal_eqdsk, nR_eqdsk, nZ_eqdsk, n_wall, ier)
  if (ier .eq. 0) then
    allocate(R_eqdsk(nR_eqdsk), Z_eqdsk(nZ_eqdsk), psi_eqdsk(nR_eqdsk,nZ_eqdsk))
    call get_data_from_eqdsk(normal_eqdsk, normal_eqdsk_wall, nR_eqdsk, nZ_eqdsk, R_eqdsk, Z_eqdsk, psi_eqdsk, n_wall, R_wall(1:n_wall), Z_wall(1:n_wall), ier)
  endif
  
  ! --- Manipulate the Psi from eqdsk
  do i = 1, nR_eqdsk
    RR = R_eqdsk(i)
    do j = 1, nZ_eqdsk
      ZZ = Z_eqdsk(j)
      dPsi = 0.d0
      do k = 1, 5
        amp = manipulate_psi_map(k,1)
        Rm  = manipulate_psi_map(k,2)
        Zm  = manipulate_psi_map(k,3)
        dRm = manipulate_psi_map(k,4)
        dZm = manipulate_psi_map(k,5)
        dPsi = dPsi + amp * exp(-(RR-Rm)**2/dRm**2-(ZZ-Zm)**2/dZm**2)
      end do
      psi_eqdsk(i,j) = psi_eqdsk(i,j) + dPsi
    end do
  end do
  
  ! --- We don't want to use the eqdsk wall, we want the user to be able to modify the wall
  ! --- eg. take pieces out if they are physically irrelevant...
  if (n_limiter .lt. 4) then
    write(*,*)'Problem: you need to provide the wall contour inside your'
    write(*,*)'         input file, using (n_limiter,R_limiter,Z_limiter).'
    write(*,*)'         If you want, you can use the limiter points provided'
    write(*,*)'         for you in the file limiter_file.txt and limiter_file.vtk'
    write(*,*)'         (this is based on the eqdsk wall)'
    write(*,*)'         Aborting...'
    ! --- Just trying to be nice...
    open(101, file='limiter_file.txt')
      write(101,'(A,i6)') 'n_limiter = ',n_wall
      do i=1,n_wall
        write(101,'(A,i6,A,e20.8,A,i6,A,e20.8)') 'R_limiter(',i,') = ',R_wall(i),', Z_limiter(',i,') = ',Z_wall(i)
      enddo
    close(101)
    open(101, file='limiter_file.vtk')
      write(101,'(A)')       '# vtk DataFile Version 1.0'
      write(101,'(A)')       'Line representation of vtk'
      write(101,'(A)')       'ASCII'
      write(101,'(A)')       ''
      write(101,'(A)')       'DATASET POLYDATA'
      write(101,'(A,i6,A)')  'POINTS ',n_wall,' float'
      do i=1,n_wall
        write(101,'(f15.7,A,f15.7,A,f15.7)') R_wall(i),'  ',Z_wall(i),'  ',0.d0
      enddo
      write(101,'(A)')       ''
      n_tmp  = (n_wall+1) / 10
      n_tmp2 = 12*n_tmp
      if (mod(n_wall+1,10) .ne. 0) n_tmp  = n_tmp  + 1
      if (mod(n_wall+1,10) .ne. 0) n_tmp2 = n_tmp2 + mod(n_wall+1,10) + 1
      write(101,'(A,2i6)')   'LINES ',n_tmp,n_tmp2
      do i=1,n_tmp-1
        write(101,'(A,11i6)') '11 ',10*(i-1),10*(i-1)+1,10*(i-1)+2,10*(i-1)+3,10*(i-1)+4,10*(i-1)+5,10*(i-1)+6,10*(i-1)+7,10*(i-1)+8,10*(i-1)+9,10*(i-1)+10
      enddo
      if (mod(n_wall+1,10) .ne. 0) then
        write(101,'(1i6)', advance='no') mod(n_wall+1,10)
        do i=1,mod(n_wall+1,10)-1
          write(101,'(1i6)', advance='no') (n_tmp-1)*10+i-1
        enddo
        write(101,'(1i6)') 0
      else
        i = n_tmp
        write(101,'(A,11i6)') '10 ',10*(i-1),10*(i-1)+1,10*(i-1)+2,10*(i-1)+3,10*(i-1)+4,10*(i-1)+5,10*(i-1)+6,10*(i-1)+7,10*(i-1)+8,0
      endif
    close(101)
    stop
  else
    n_wall = n_limiter
    R_wall(1:n_limiter) = R_limiter(1:n_limiter)
    Z_wall(1:n_limiter) = Z_limiter(1:n_limiter)
  endif
  
  
  
  
  ! --------------------------------------------------------------------------
  ! --- Build a grid inside the wall (not the case below has jumps for MAST-U)
  ! --------------------------------------------------------------------------
  
  
  ! --- Input parameters
  nR = n_R
  nZ = n_Z - mod(n_Z,2) ! This needs to be even!!!
  ! --- nR is the resolution we want in the core, but there can be jumps
  ! --- We start at midplane and go down (or up)
  allocate(nR_grid(nZ+1,2),node_index(4*nR*nZ,4,2),elm_node_index(4*nR,nZ+1))
  allocate(R_grid(4*nR,nZ+1),Z_grid(4*nR,nZ+1)) ! need some margin for nR due to jumps
  allocate(Zlines(nZ/2+1))
  
  if (tokamak_device(1:6) .eq. 'MAST-U') then
    nR_grid(1:nZ+1,2) = 0
    call create_grid_inside_wall_MASTU(nR, nZ, nR_grid(1:nZ+1,1), node_index, Zlines, R_grid, Z_grid, n_elm)
  else if (tokamak_device(1:4) .eq. 'STEP') then
    call create_grid_inside_wall_STEP(nR, nZ, nR_grid, node_index, R_grid, Z_grid, n_elm)
  else
    call create_grid_inside_wall_usual(nR, nZ, nR_grid, node_index, Zlines, R_grid, Z_grid, n_elm)
  endif
  
  
  
  ! --- Convert into jorek grid
  node_list%n_nodes = 0
  do i = 1,nZ+1
    do j = 1,nR_grid(i,1)+nR_grid(i,2)
      node_list%n_nodes = node_list%n_nodes + 1
      elm_node_index(j,i) = node_list%n_nodes
      ! --- Position
      node_list%node(node_list%n_nodes)%x(1,1,1) = R_grid(j,i)
      node_list%node(node_list%n_nodes)%x(1,1,2) = Z_grid(j,i)
    enddo
  enddo
      
  element_list%n_elements = n_elm
  node_list%node(1:n_nodes_max)%values(1,1,1) = 0.d0
  do i_elm = 1,n_elm
    
    ! --- Nodes
    do j = 1,4
      ii_n(j) = node_index(i_elm,j,1)
      jj_n(j) = node_index(i_elm,j,2)
      element_list%element(i_elm)%vertex(j) = elm_node_index(ii_n(j),jj_n(j))
      R_tmp(j) =  R_grid(ii_n(j),jj_n(j))
      Z_tmp(j) =  Z_grid(ii_n(j),jj_n(j))
    enddo
    
    ! --- vectors
    do j = 1,4
      i_node = element_list%element(i_elm)%vertex(j)
      if (node_list%node(i_node)%values(1,1,1) .ne. 0.d0) cycle
      
      ! --- vector u
      if (j .eq. 1) jp1 = 2
      if (j .eq. 2) jp1 = 1
      if (j .eq. 3) jp1 = 4
      if (j .eq. 4) jp1 = 3
      width = sqrt((R_tmp(jp1) - R_tmp(j))**2 + (Z_tmp(jp1) - Z_tmp(j))**2)
      node_list%node(i_node)%x(1,2,1) = (R_tmp(jp1) - R_tmp(j)) / width
      node_list%node(i_node)%x(1,2,2) = (Z_tmp(jp1) - Z_tmp(j)) / width
      ! --- vector v
      if (j .eq. 1) jp1 = 4
      if (j .eq. 2) jp1 = 3
      if (j .eq. 3) jp1 = 2
      if (j .eq. 4) jp1 = 1
      width = sqrt((R_tmp(jp1) - R_tmp(j))**2 + (Z_tmp(jp1) - Z_tmp(j))**2)
      node_list%node(i_node)%x(1,3,1) = (R_tmp(jp1) - R_tmp(j)) / width
      node_list%node(i_node)%x(1,3,2) = (Z_tmp(jp1) - Z_tmp(j)) / width
      ! --- vector w
      node_list%node(i_node)%x(1,4,1) = 0.d0
      node_list%node(i_node)%x(1,4,2) = 0.d0
      
      ! --- boundary
      node_list%node(i_node)%boundary = 0
      
      ! --- matrix index
      do k=1, n_degrees
        node_list%node(i_node)%index(k) = n_degrees*(i_node-1) + k
      enddo
      
      ! --- psi values from eqdsk
      call interpolate_psi_from_eqdsk_grid(nR_eqdsk, nZ_eqdsk, R_eqdsk, Z_eqdsk, psi_eqdsk, R_tmp(j), Z_tmp(j), psi, psi_R, psi_Z)
      
      ! --- values in grid
      node_list%node(i_node)%values(1,1,1) = psi
      node_list%node(i_node)%values(1,2,1) = psi_R * node_list%node(i_node)%x(1,2,1) + psi_Z * node_list%node(i_node)%x(1,2,2)
      node_list%node(i_node)%values(1,3,1) = psi_R * node_list%node(i_node)%x(1,3,1) + psi_Z * node_list%node(i_node)%x(1,3,2)
      node_list%node(i_node)%values(1,4,1) = psi_R * node_list%node(i_node)%x(1,4,1) + psi_Z * node_list%node(i_node)%x(1,4,2)
    enddo
                       
    
  enddo
  


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
      element_list%element(k)%size(iv,iuv+1) = sign(dlength(xx_p,xx_0),ddot(n_dim,xx_p - xx_0,1,uv_0,1)) / float(n_order) 
      element_list%element(k)%size(ip,iuv+1) = sign(dlength(xx_p,xx_0),ddot(n_dim,xx_0 - xx_p,1,uv_p,1)) / float(n_order) 

    enddo

    do iv=1,4
      element_list%element(k)%size(iv,4) = element_list%element(k)%size(iv,2) * element_list%element(k)%size(iv,3)
    enddo
  enddo
  
  if (n_order .ge. 5) then
    call set_high_order_sizes(element_list)
    ! this grid is better without 2nd derivatives...
    call approximate_2nd_derivatives(node_list,element_list)
    do i=1,node_list%n_nodes
      node_list%node(i)%x(1,7:n_degrees,:) = 0.d0
    enddo
  endif
  
  ! --- Update neighbours and boundary
  call update_neighbours_basic(element_list,node_list)
  call update_boundary_types(element_list,node_list, 0)
  
  ! --- MAST-U specific (due to stupid equilibria with nonsense p and ffp profiles...)
  ! --- Shift plasma away from solenoid
  if ( (tokamak_device(1:6) .eq. 'MAST-U') .and. (abs(xshift) .gt. 0.d0) ) then
    write(*,*)'WARNING!!! shifting plasma, use xshift=0.0 if you do not want to!' 
    do inode=1, node_list%n_nodes
      if ( (node_list%node(inode)%boundary .gt. 0) .and. (node_list%node(inode)%x(1,1,1) .lt. 0.37) ) then
      !if ( (node_list%node(inode)%boundary .gt. 0) .and. (abs(node_list%node(inode)%x(1,1,2)) .le. 1.2) ) then
        Zmin  = -1.1!-0.7
        z_min = -0.3!-0.4
        Zmax  = +1.1!+0.7
        z_max = +0.3!+0.4
        R_elm = node_list%node(inode)%x(1,1,1)
        Z_elm = node_list%node(inode)%x(1,1,2)
        diff  = xshift
        if ( (Z_elm .le. z_min) .and. (Z_elm .ge. Zmin) ) then
          diff = (Z_elm-Zmin)/(z_min-Zmin)*diff
        elseif ( (Z_elm .ge. z_max) .and. (Z_elm .le. Zmax) ) then
          diff = (Z_elm-Zmax)/(z_max-Zmax)*diff
        endif
        if ( (Z_elm .lt. Zmin) .or. (Z_elm .gt. Zmax) ) then
          diff = 0.d0
        endif
        node_list%node(inode)%values(1,1,1) = node_list%node(inode)%values(1,1,1) + diff * node_list%node(inode)%values(1,1,1)
      endif
    enddo
  endif

  ! --- Free data
  deallocate(nR_grid, node_index, elm_node_index)
  deallocate(R_grid, Z_grid)
  deallocate(Zlines)
  
  do inode=1, node_list%n_nodes
    node_list%node(inode)%axis_node = .false.
  enddo
  
  ! --- This is just for debug, it could be removed (or not?)
  call export_restart(node_list, element_list, 'grid_initial')
  
  ! --- Also for debug: Print a python file that plots a cross with the 4 nodes of each element
  if (plot_grid) then
    n_loop = element_list%n_elements
    open(101,file='plot_grid_inside_wall.py')
      write(101,'(A)')         '#!/usr/bin/env python'
      write(101,'(A)')         'import numpy as N'
      write(101,'(A)')         'import pylab'
      write(101,'(A)')         'def main():'
      write(101,'(A,i6,A)')    ' r = N.zeros(',4*n_loop,')'
      write(101,'(A,i6,A)')    ' z = N.zeros(',4*n_loop,')'
      do j=1,n_loop
        do i=1,2
          index = element_list%element(j)%vertex(i)
          write(101,'(A,i6,A,f15.4)') ' r[',4*(j-1)+2*i-2,'] = ',node_list%node(index)%x(1,1,1)
          write(101,'(A,i6,A,f15.4)') ' z[',4*(j-1)+2*i-2,'] = ',node_list%node(index)%x(1,1,2)
          index = element_list%element(j)%vertex(i+2)
          write(101,'(A,i6,A,f15.4)') ' r[',4*(j-1)+2*i-1,'] = ',node_list%node(index)%x(1,1,1)
          write(101,'(A,i6,A,f15.4)') ' z[',4*(j-1)+2*i-1,'] = ',node_list%node(index)%x(1,1,2)
        enddo
      enddo
      write(101,'(A,i6,A)')    ' for i in range (0,',n_loop,'):'
      write(101,'(A)')         '  pylab.plot(r[4*i:4*i+2],z[4*i:4*i+2], "r")'
      write(101,'(A)')         '  pylab.plot(r[4*i+2:4*i+4],z[4*i+2:4*i+4], "g")'
      do j=1,n_loop
        do i=1,4
          index = element_list%element(j)%vertex(i)
          write(101,'(A,i6,A,f15.4)') ' r[',4*(j-1)+i-1,'] = ',node_list%node(index)%x(1,1,1)
          write(101,'(A,i6,A,f15.4)') ' z[',4*(j-1)+i-1,'] = ',node_list%node(index)%x(1,1,2)
        enddo
      enddo
      write(101,'(A,i6,A)')    ' for i in range (0,',n_loop,'):'
      write(101,'(A)')         '  pylab.plot(r[4*i:4*i+4],z[4*i:4*i+4], "b")'
      write(101,'(A)')         ' pylab.axis("equal")'
      write(101,'(A)')         ' pylab.show()'
      write(101,'(A)')         ' '
      write(101,'(A)')         'main()'
    close(101)
  endif

  return

end subroutine grid_inside_wall
  






! --- Find the polygon segment a point belongs to
integer function inbetweenSegment(N, rpol, zpol, R, Z)
  
  implicit none
  
  ! --- Input variables
  integer, intent(in) :: N
  real*8,  intent(in) :: rpol(N), zpol(N), R, Z
  
  ! --- Local variables
  integer :: i
  real*8  :: rmin, rmax, zmin, zmax
  
  inbetweenSegment = 0
  do i = 1,N-1
    rmin = min(rpol(i),rpol(i+1))
    rmax = max(rpol(i),rpol(i+1))
    zmin = min(zpol(i),zpol(i+1))
    zmax = max(zpol(i),zpol(i+1))
    if ( (R .ge. rmin) .and. (R .le. rmax) .and. (Z .ge. zmin) .and. (Z .le. zmax) ) then
      inbetweenSegment = i
      exit
    endif
  enddo
  
  return
end function inbetweenSegment




! --- Find intersection between horizontal segment (Rmin,Rmax)@height=Z and polygon
subroutine RintersectPolygon(N, rpol, zpol, Rmin, Rmax, Z, accuracy, new_Rmin, new_Rmax, new_Rmin2, new_Rmax2)

  implicit none
  
  ! --- Input variables
  integer, intent(in)  :: N
  real*8,  intent(in)  :: rpol(N), zpol(N), Rmin, Rmax, Z, accuracy
  real*8,  intent(out) :: new_Rmin, new_Rmax, new_Rmin2, new_Rmax2
  
  ! --- Local variables
  integer, parameter :: n_iter = 1000
  integer :: inside, i1, i2, i, insidePolygon
  real*8  :: R_iter(n_iter), r_tmp, z_tmp
  real*8  :: Rleft, Rmid, Rright
  real*8  :: diff
  
  new_Rmin  = 0.d0
  new_Rmax  = 0.d0
  new_Rmin2 = 0.d0
  new_Rmax2 = 0.d0
  
  ! --- Find entry and exit points of horizontal line (in/out of polygon)
  inside = 0
  do i = 1,n_iter
    R_iter(i) = Rmin + real(i-1)/real(n_iter-1) * (Rmax-Rmin)
  enddo
  z_tmp  = Z
  i1 = -1
  i2 = -1
  do i=1,n_iter
    r_tmp = R_iter(i)
    if (insidePolygon(N, rpol, zpol, r_tmp, z_tmp) .eq. 1) then
      if (inside .eq. 0) then
        i1 = i
        inside = 1
      endif
    endif
    if (insidePolygon(N, rpol, zpol, r_tmp, z_tmp) .eq. 0) then
      if (inside .eq. 1) then
        i2 = i
        exit
      endif
    endif
  enddo
  if (i1 .eq. -1) then
    write(*,*)'Warning: did not find entry'
    return
  endif
  if (i2 .eq. -1) then
    write(*,*)'Warning: did not find exit'
    return
  endif
  
  ! --- Converge to solution at entry point
  Rleft  = R_iter(i1-1)
  Rright = R_iter(i1)
  Rmid   = 0.5*(Rleft+Rright)
  do i = 1,n_iter
    diff = abs(Rmid - 0.5*(Rleft+Rright))
    Rmid = 0.5*(Rleft+Rright)
    if ( (diff .lt. accuracy) .and. (i .gt. 1) ) then
      new_Rmin = Rmid
      exit
    endif
    if (insidePolygon(N, rpol, zpol, Rmid, z_tmp) .eq. 1) then
      Rright = Rmid
    else
      Rleft  = Rmid
    endif
    if (i  .eq. n_iter) then
      write(*,*)'WARNING: RintersectPolygon did not converge in entry search!'
    endif
  enddo
  
  ! --- Converge to solution at exit point
  Rleft  = R_iter(i2-1)
  Rright = R_iter(i2)
  Rmid   = 0.5*(Rleft+Rright)
  do i = 1,n_iter
    diff = abs(Rmid - 0.5*(Rleft+Rright))
    Rmid = 0.5*(Rleft+Rright)
    if ( (diff .lt. accuracy) .and. (i .gt. 1) ) then
      new_Rmax = Rmid
      exit
    endif
    if (insidePolygon(N, rpol, zpol, Rmid, z_tmp) .eq. 1) then
      Rleft = Rmid
    else
      Rright  = Rmid
    endif
    if (i  == n_iter-1) then
      write(*,*)'WARNING: RintersectPolygon did not converge in exit search!'
    endif
  enddo
  
  ! --- We do it again if we have two intersection pairs
  new_Rmin2 = -1.e10
  new_Rmax2 = -1.e10
  
  ! --- Find entry and exit points of horizontal line (in/out of polygon)
  inside = 0
  do i = 1,n_iter
    R_iter(i) = Rmax + real(i-1)/real(n_iter-1) * (Rmin-Rmax)
  enddo
  z_tmp  = Z
  i1 = -1
  i2 = -1
  do i=1,n_iter
    r_tmp = R_iter(i)
    if (insidePolygon(N, rpol, zpol, r_tmp, z_tmp) .eq. 1) then
      if (inside .eq. 0) then
        i1 = i
        inside = 1
      endif
    endif
    if (insidePolygon(N, rpol, zpol, r_tmp, z_tmp) .eq. 0) then
      if (inside .eq. 1) then
        i2 = i
        exit
      endif
    endif
  enddo
  if (i1 .eq. -1) then
    write(*,*)'Warning: did not find entry'
    i1 = 1
  endif
  if (i2 .eq. -1) then
    write(*,*)'Warning: did not find exit'
    i2 = 1
  endif
  i = i1
  i1 = i2
  i2 = i
  if (i1 .eq. 1) i1 = 2
  if (i2 .eq. 1) i2 = 2
  
  ! --- Converge to solution at entry point
  Rleft  = R_iter(i1)
  Rright = R_iter(i1-1)
  Rmid   = 0.5*(Rleft+Rright)
  do i = 1,n_iter
    diff = abs(Rmid - 0.5*(Rleft+Rright))
    Rmid = 0.5*(Rleft+Rright)
    if ( (diff .lt. accuracy) .and. (i .gt. 1) ) then
      new_Rmin2 = Rmid
      exit
    endif
    if (insidePolygon(N, rpol, zpol, Rmid, z_tmp) .eq. 1) then
      Rright = Rmid
    else
      Rleft  = Rmid
    endif
    if (i  .eq. n_iter) then
      write(*,*)'WARNING: RintersectPolygon did not converge in entry search!'
    endif
  enddo
  
  ! --- Converge to solution at exit point
  Rleft  = R_iter(i2)
  Rright = R_iter(i2-1)
  Rmid   = 0.5*(Rleft+Rright)
  do i = 1,n_iter
    diff = abs(Rmid - 0.5*(Rleft+Rright))
    Rmid = 0.5*(Rleft+Rright)
    if ( (diff .lt. accuracy) .and. (i .gt. 1) ) then
      new_Rmax2 = Rmid
      exit
    endif
    if (insidePolygon(N, rpol, zpol, Rmid, z_tmp) .eq. 1) then
      Rleft = Rmid
    else
      Rright  = Rmid
    endif
    if (i  == n_iter-1) then
      write(*,*)'WARNING: RintersectPolygon did not converge in exit search!'
    endif
  enddo
  
  if ( (abs(new_Rmin-new_Rmin2) .lt. 1.d-2) .and. (abs(new_Rmax-new_Rmax2) .lt. 1.d-2)) then
    new_Rmin2 = -1.e10
    new_Rmax2 = -1.e10
  endif
  
  return
end subroutine RintersectPolygon





! --- Check if point is inside polygon
integer function insidePolygon(N, rpol, zpol, r, z)

  implicit none
  
  ! --- Input variables
  integer, intent(in)  :: N
  real*8,  intent(in)  :: rpol(N), zpol(N), r, z
  
  ! --- Local variables
  integer :: counter, i
  real*8  :: r1, z1, r2, z2, xinters
  
  insidePolygon = 0
  
  counter = 0
  r1 = rpol(N)
  z1 = zpol(N)
  do i = 1,N
    r2 = rpol(i)
    z2 = zpol(i)
    if ( (z .ge. min(z1,z2)) .and. (z .le. max(z1,z2)) .and. (r .lt. max(r1,r2)) ) then
      if (z1 .ne. z2) then
        xinters = (z-z1)*(r2-r1)/(z2-z1) + r1
        if ((r .lt. xinters) .and. (z .ne. z2)) then
          counter = counter + 1
        endif
      endif
    endif
    r1 = r2
    z1 = z2
  enddo

  ! --- Outside
  if (mod(counter,2) .eq. 0) then
    insidePolygon = 0
  ! --- Inside
  else
    insidePolygon = 1
  endif
  
  return

end function insidePolygon





















subroutine create_grid_inside_wall_usual(nR, nZ, nR_grid, node_index, Zlines, R_grid, Z_grid, n_elm)


  use grid_xpoint_data, only: n_wall, R_wall, Z_wall
  use phys_module, only : RZ_grid_jump_thres
  
  implicit none
  
  ! --- Input variables
  integer, intent(inout)  :: nR, nZ, n_elm
  integer, intent(inout)  :: nR_grid(nZ+1,2),node_index(4*nR*nZ,4,2)
  real*8,  intent(inout)  :: R_grid(4*nR,nZ+1),Z_grid(4*nR,nZ+1), Zlines(nZ/2+1)
  
  ! --- Local variables
  integer :: i, j
  integer :: nR_save, n_off, i_seg, nR_tmp, i_start_tmp, i_save, j_save, i_elm, i_node, elm_count, i_start
  integer :: iRp1, iZp1
  integer :: inside_boot, n_extra, inbetweenSegment
  real*8  :: Rmin, Rmax, r_min1, r_max1, r_min2, r_max2
  real*8  :: Zmin, Zmax, z_min, z_max, Z_first_jump, jump_threshold
  real*8  :: diff, width
  real*8  :: width1, width_prev1
  real*8  :: width2, width_prev2
  real*8  :: accuracy
  real*8  :: distance, distance_min
  real*8  :: sig_Z
  
  
  write(*,*)'Building grid inside wall'
  
  R_grid   = 0.d0
  Z_grid   = 0.d0
  nR_grid  = 0
  accuracy = 1.d-5
  
  ! --- We cut the domain with horizontal lines
  Rmin = minval(R_wall(1:n_wall))-1.e-3 ; Rmax = maxval(R_wall(1:n_wall))+1.e-3 ! want them slightly outside wall
  Zmin = minval(Z_wall(1:n_wall))+1.e-3 ; Zmax = maxval(Z_wall(1:n_wall))-1.e-3 ! want them slightly inside wall
  write(*,*)'Rminmax:',Rmin, Rmax
  write(*,*)'Zminmax:',Zmin, Zmax
  
  ! --- First lower part
  nR_save = nR
  sig_Z = 0.3
  call meshac2(nZ/2+1,Zlines,1.d0,9999.d0,sig_Z,9999.d0,0.3,1.0d0)
  Zlines = Zlines * Zmin
  width_prev1    = 0.0
  width_prev2    = 0.0
  inside_boot    = 0
  jump_threshold = RZ_grid_jump_thres
  elm_count      = 0
  do i = 1,nZ/2+1
    !write(*,*)'Doing lower horizontal line ',nZ/2+1-i
    n_off = nZ/2-1
    call RintersectPolygon(n_wall, R_wall(1:n_wall), Z_wall(1:n_wall), Rmin, Rmax, Zlines(i), accuracy, r_min1, r_max1, r_min2, r_max2)
    width1 = r_max1-r_min1
    width2 = r_max2-r_min2
    if ( (width2 .eq. 0.d0) .and. (width_prev2 .ne. 0.d0) ) then
      if ( abs(r_min1-R_grid(nR_grid(n_off+i-1,1)+1,n_off+i-1)) .le. abs(r_min1-R_grid(1,n_off+i-1)) ) then
        width1 = 0.d0
        nR_grid(n_off+i,1) = 0
        r_min2 = r_min1
        r_max2 = r_max1
        width2 = r_max2-r_min2
      else
        width2 = 0.d0
        width_prev2 = 0.d0
        nR_grid(n_off+i,2) = 0
      endif
    endif
    
    if (width1 .gt. 0.d0) then
      i_start = 0
      ! Jump of resolution
      nR_grid(n_off+i,1) = nR_grid(n_off+i-1,1)
      if (i .eq. 1) nR_grid(n_off+i,1) = nR
      if ( (width1 .lt. jump_threshold*width_prev1) .and. (i .gt. 1) .and. (i .lt. nZ/2+1) ) then
        diff = width1/width_prev1
        nR_grid(n_off+i,1) = int(real(nR_grid(n_off+i-1,1)) * diff) + 1
        distance_min = 1.d10
        do j = 1,nR_grid(n_off+i-1,1)
          distance = abs(R_grid(j,n_off+i-1)-r_min1)
          if (distance .lt. distance_min) then
            distance_min = distance
            i_start = j-1
          endif
        enddo
        i_start = i_start-1
        if (i_start .lt. 0) i_start = 0
      endif
      do j = 1,nR_grid(n_off+i,1)
        R_grid(j,n_off+i) = r_min1 + real(j-1)/real(nR_grid(n_off+i,1)-1) * (r_max1 - r_min1)
        Z_grid(j,n_off+i) = Zlines(i)
      enddo
      if (i .gt. 1) then
        nR_tmp = min(nR_grid(n_off+i,1),nR_grid(n_off+i-1,1))
        do j = 1,nR_tmp-1
          elm_count = elm_count + 1
          ! R-index
          node_index(elm_count,1,1) = j
          node_index(elm_count,2,1) = j+1
          node_index(elm_count,3,1) = i_start+j+1
          node_index(elm_count,4,1) = i_start+j
          ! Z-index
          node_index(elm_count,1,2) = n_off+i
          node_index(elm_count,2,2) = n_off+i
          node_index(elm_count,3,2) = n_off+i-1
          node_index(elm_count,4,2) = n_off+i-1
        enddo
      endif
      width_prev1 = width1
    endif
    
    if (width2 .gt. 0.d0) then
      i_start = 0
      if (width_prev2 .eq. 0.d0) then
        ! Find start
        distance_min = 1.d10
        do j = 1,nR_grid(n_off+i-1,1)
          distance = abs(R_grid(j,n_off+i-1)-r_min2)
          if (distance .lt. distance_min) then
            distance_min = distance
            i_start = j-1
          endif
        enddo
        i_start = i_start - 2
        if (i_start .lt. nR_grid(n_off+i,1)-1) i_start = i_start + 1
        nR_grid(n_off+i,2) = nR_grid(n_off+i-1,1) - i_start
      else
        ! Jump of resolution
        nR_grid(n_off+i,2) = nR_grid(n_off+i-1,2)
        i_start = nR_grid(n_off+i-1,1)
        if ( (width2 .lt. jump_threshold*width_prev2) .and. (i .gt. 1) .and. (i .lt. nZ/2+1) ) then
          diff = width2/width_prev2
          nR_grid(n_off+i,2) = int(real(nR_grid(n_off+i-1,2)) * diff)
          distance_min = 1.d10
          do j = 1,nR_grid(n_off+i-1,1)
            distance = abs(R_grid(j,n_off+i-1)-r_min2)
            if (distance .lt. distance_min) then
              distance_min = distance
              i_start = j-2
            endif
          enddo
          do j = nR_grid(n_off+i-1,1)+1,nR_grid(n_off+i-1,1)+nR_grid(n_off+i-1,2)
            distance = abs(R_grid(j,n_off+i-1)-r_min2)
            if (distance .lt. distance_min) then
              distance_min = distance
              i_start = j-2
            endif
          enddo
        endif
        nR_tmp = nR_grid(n_off+i-1,1)+nR_grid(n_off+i-1,2)
        if (abs(R_grid(nR_tmp,n_off+i-1)-r_max2) .lt. 1.d-2) then
          nR_grid(n_off+i,2) = nR_tmp - i_start
        endif
      endif
      do j = 1,nR_grid(n_off+i,2)
        R_grid(nR_grid(n_off+i,1)+j,n_off+i) = r_min2 + real(j-1)/real(nR_grid(n_off+i,2)-1) * (r_max2 - r_min2)
        Z_grid(nR_grid(n_off+i,1)+j,n_off+i) = Zlines(i)
      enddo
      if (i .gt. 1) then
        nR_tmp = min(nR_grid(n_off+i,2),nR_grid(n_off+i-1,1)+nR_grid(n_off+i-1,2)-i_start+1)
        if (width_prev2 .eq. 0.d0) nR_tmp = nR_grid(n_off+i,2)
        do j = 1,nR_tmp-1
          elm_count = elm_count + 1
          ! R-index
          node_index(elm_count,1,1) = nR_grid(n_off+i,1)+j
          node_index(elm_count,2,1) = nR_grid(n_off+i,1)+j+1
          node_index(elm_count,3,1) = i_start+j+1
          node_index(elm_count,4,1) = i_start+j
          ! Z-index
          node_index(elm_count,1,2) = n_off+i
          node_index(elm_count,2,2) = n_off+i
          node_index(elm_count,3,2) = n_off+i-1
          node_index(elm_count,4,2) = n_off+i-1
        enddo
      endif
      width_prev2 = width2
    endif
    
    ! --- This may happen if there is a sharp corner at the bottom (eg. ITER)
    if ( (width1 .eq. 0.d0) .and. (width2 .eq. 0.d0) ) then
      if (width_prev1 .gt. 0.d0) then
        nR_grid(n_off+i,1) = nR_grid(n_off+i-1,1)
        do j = 1,nR_grid(n_off+i,1)
          R_grid(j,n_off+i) = R_grid(j,n_off+i-1)
          Z_grid(j,n_off+i) = Zlines(i)
        enddo
        nR_tmp = min(nR_grid(n_off+i,1),nR_grid(n_off+i-1,1))
        do j = 1,nR_tmp-1
          elm_count = elm_count + 1
          ! R-index
          node_index(elm_count,1,1) = j
          node_index(elm_count,2,1) = j+1
          node_index(elm_count,3,1) = i_start+j+1
          node_index(elm_count,4,1) = i_start+j
          ! Z-index
          node_index(elm_count,1,2) = n_off+i
          node_index(elm_count,2,2) = n_off+i
          node_index(elm_count,3,2) = n_off+i-1
          node_index(elm_count,4,2) = n_off+i-1
        enddo
      endif
      if (width_prev2 .gt. 0.d0) then
        nR_grid(n_off+i,2) = nR_grid(n_off+i-1,2)
        i_start = nR_grid(n_off+i-1,1)
        nR_tmp = nR_grid(n_off+i-1,1)+nR_grid(n_off+i-1,2)
        if (abs(R_grid(nR_tmp,n_off+i-1)-r_max2) .lt. 1.d-2) then
          nR_grid(n_off+i,2) = nR_tmp - i_start
        endif
        do j = 1,nR_grid(n_off+i,2)
          R_grid(nR_grid(n_off+i,1)+j,n_off+i) = R_grid(nR_grid(n_off+i,1)+j,n_off+i-1)
          Z_grid(nR_grid(n_off+i,1)+j,n_off+i) = Zlines(i)
        enddo
        nR_tmp = min(nR_grid(n_off+i,2),nR_grid(n_off+i-1,1)+nR_grid(n_off+i-1,2)-i_start+1)
        if (width_prev2 .eq. 0.d0) nR_tmp = nR_grid(n_off+i,2)
        do j = 1,nR_tmp-1
          elm_count = elm_count + 1
          ! R-index
          node_index(elm_count,1,1) = nR_grid(n_off+i,1)+j
          node_index(elm_count,2,1) = nR_grid(n_off+i,1)+j+1
          node_index(elm_count,3,1) = i_start+j+1
          node_index(elm_count,4,1) = i_start+j
          ! Z-index
          node_index(elm_count,1,2) = n_off+i
          node_index(elm_count,2,2) = n_off+i
          node_index(elm_count,3,2) = n_off+i-1
          node_index(elm_count,4,2) = n_off+i-1
        enddo
      endif
    endif
    
  enddo
  write(*,*)'Finished lower part of grid'
  
  ! --- Then upper part
  nR = nR_save
  sig_Z = 0.3
  call meshac2(nZ/2,Zlines,1.d0,9999.d0,sig_Z,9999.d0,0.3,1.0d0)
  Zlines = Zlines * Zmax
  width_prev1    = 0.0
  width_prev2    = 0.0
  inside_boot    = 0
  do i = 1,nZ/2
    !write(*,*)'Doing upper horizontal line ',nZ/2+1-i
    n_off = nZ/2+1
    call RintersectPolygon(n_wall, R_wall(1:n_wall), Z_wall(1:n_wall), Rmin, Rmax, Zlines(i), accuracy, r_min1, r_max1, r_min2, r_max2)
    width1 = r_max1-r_min1
    width2 = r_max2-r_min2
    if ( (width2 .eq. 0.d0) .and. (width_prev2 .ne. 0.d0) ) then
      width1 = 0.d0
      nR_grid(n_off-i,1) = 0
      r_min2 = r_min1
      r_max2 = r_max1
      width2 = r_max2-r_min2
    endif
    
    if (width1 .gt. 0.d0) then
      i_start = 0
      ! Jump of resolution
      nR_grid(n_off-i,1) = nR_grid(n_off-i+1,1)
      if (i .eq. 1) nR_grid(n_off-i,1) = nR
      if ( (width1 .lt. jump_threshold*width_prev1) .and. (i .gt. 1) .and. (i .lt. nZ/2) ) then
        diff = width1/width_prev1
        nR_grid(n_off-i,1) = int(real(nR_grid(n_off-i+1,1)) * diff) + 1
        distance_min = 1.d10
        do j = 1,nR_grid(n_off-i+1,1)
          distance = abs(R_grid(j,n_off-i+1)-r_min1)
          if (distance .lt. distance_min) then
            distance_min = distance
            i_start = j-1
          endif
        enddo
      endif
      do j = 1,nR_grid(n_off-i,1)
        R_grid(j,n_off-i) = r_min1 + real(j-1)/real(nR_grid(n_off-i,1)-1) * (r_max1 - r_min1)
        Z_grid(j,n_off-i) = Zlines(i)
      enddo
      if (i .gt. 1) then
        nR_tmp = min(nR_grid(n_off-i,1),nR_grid(n_off-i+1,1))
        do j = 1,nR_tmp-1
          elm_count = elm_count + 1
          ! R-index
          node_index(elm_count,1,1) = i_start+j
          node_index(elm_count,2,1) = i_start+j+1
          node_index(elm_count,3,1) = j+1
          node_index(elm_count,4,1) = j
          ! Z-index
          node_index(elm_count,1,2) = n_off-i+1
          node_index(elm_count,2,2) = n_off-i+1
          node_index(elm_count,3,2) = n_off-i
          node_index(elm_count,4,2) = n_off-i
        enddo
      endif
      width_prev1 = width1
    endif
    
    if (width2 .gt. 0.d0) then
      i_start = 0
      if (width_prev2 .eq. 0.d0) then
        ! Find start
        distance_min = 1.d10
        do j = 1,nR_grid(n_off-i+1,1)
          distance = abs(R_grid(j,n_off-i+1)-r_min2)
          if (distance .lt. distance_min) then
            distance_min = distance
            i_start = j-1
          endif
        enddo
        if (i_start .lt. nR_grid(n_off-i,1)) i_start = i_start + 1
        nR_grid(n_off-i,2) = nR_grid(n_off-i+1,1) - i_start
      else
        ! Jump of resolution
        nR_grid(n_off-i,2) = nR_grid(n_off-i+1,2)
        i_start = nR_grid(n_off-i+1,1)
        if ( (width2 .lt. jump_threshold*width_prev2) .and. (i .gt. 1) .and. (i .lt. nZ/2) ) then
          diff = width2/width_prev2
          nR_grid(n_off-i,2) = int(real(nR_grid(n_off-i+1,2)) * diff)
          distance_min = 1.d10
          do j = 1,nR_grid(n_off-i+1,1)
            distance = abs(R_grid(j,n_off-i+1)-r_min2)
            if (distance .lt. distance_min) then
              distance_min = distance
              i_start = j-2
            endif
          enddo
          do j = nR_grid(n_off-i+1,1)+1,nR_grid(n_off-i+1,1)+nR_grid(n_off-i+1,2)
            distance = abs(R_grid(j,n_off-i+1)-r_min2)
            if (distance .lt. distance_min) then
              distance_min = distance
              i_start = j-2
            endif
          enddo
        endif
        nR_tmp = nR_grid(n_off-i+1,1)+nR_grid(n_off-i+1,2)
        if (abs(R_grid(nR_tmp,n_off-i+1)-r_max2) .lt. 1.d-2) then
          nR_grid(n_off-i,2) = nR_tmp - i_start
        endif
      endif
      do j = 1,nR_grid(n_off-i,2)
        R_grid(nR_grid(n_off-i,1)+j,n_off-i) = r_min2 + real(j-1)/real(nR_grid(n_off-i,2)-1) * (r_max2 - r_min2)
        Z_grid(nR_grid(n_off-i,1)+j,n_off-i) = Zlines(i)
      enddo
      if (i .gt. 1) then
        nR_tmp = min(nR_grid(n_off-i,2),nR_grid(n_off-i+1,1)+nR_grid(n_off-i+1,2)-i_start+1)
        if (width_prev2 .eq. 0.d0) nR_tmp = nR_grid(n_off-i,2)
        do j = 1,nR_tmp-1
          elm_count = elm_count + 1
          ! R-index
          node_index(elm_count,1,1) = nR_grid(n_off-i,1)+j
          node_index(elm_count,2,1) = nR_grid(n_off-i,1)+j+1
          node_index(elm_count,3,1) = i_start+j+1
          node_index(elm_count,4,1) = i_start+j
          ! Z-index
          node_index(elm_count,1,2) = n_off-i
          node_index(elm_count,2,2) = n_off-i
          node_index(elm_count,3,2) = n_off-i+1
          node_index(elm_count,4,2) = n_off-i+1
        enddo
      endif
      width_prev2 = width2
    endif
    
  enddo
  write(*,*)'Finished upper part of grid'
  
  
  n_elm = elm_count


  return

end subroutine create_grid_inside_wall_usual










subroutine create_grid_inside_wall_MASTU(nR, nZ, nR_grid, node_index, Zlines, R_grid, Z_grid, n_elm)


  use grid_xpoint_data, only: n_wall, R_wall, Z_wall
  
  implicit none
  
  ! --- Input variables
  integer, intent(inout)  :: nR, nZ, n_elm
  integer, intent(inout)  :: nR_grid(nZ+1),node_index(4*nR*nZ,4,2)
  real*8,  intent(inout)  :: R_grid(4*nR,nZ+1),Z_grid(4*nR,nZ+1), Zlines(nZ/2+1)
  
  ! --- Local variables
  integer :: i, j
  integer :: nR_save, n_off, i_seg, nR_tmp, i_start_tmp, i_save, j_save, i_elm, i_node, elm_count, i_start
  integer :: iRp1, iZp1
  integer :: inside_boot, n_extra, inbetweenSegment
  real*8  :: Rmin, Rmax, r_min, r_max, r_min2, r_max2
  real*8  :: Zmin, Zmax, z_min, z_max, Z_first_jump, jump_threshold, z_tmp
  real*8  :: diff
  real*8  :: width, width_prev
  real*8  :: accuracy
  real*8  :: distance, distance_min
  real*8  :: sig_Z
  
  
  write(*,*)'Building grid inside wall'
  
  ! --- Input parameters
  Z_first_jump = 1.0!meters
  R_grid   = 0.d0
  Z_grid   = 0.d0
  nR_grid  = 0
  accuracy = 1.d-5
  
  ! --- We cut the domain with horizontal lines
  Rmin = minval(R_wall(1:n_wall))-1.e-3 ; Rmax = maxval(R_wall(1:n_wall))+1.e-3 ! want them slightly outside wall
  Zmin = minval(Z_wall(1:n_wall))+1.e-3 ; Zmax = maxval(Z_wall(1:n_wall))-1.e-3 ! want them slightly inside wall
  write(*,*)'Rminmax:',Rmin, Rmax
  write(*,*)'Zminmax:',Zmin, Zmax
  
  ! --- First lower part
  nR_save = nR
  do i=1,nZ/2+1
    Zlines(i) = real(i-1)/real(nZ/2) * Zmin
  enddo
  width_prev  = 0.0
  inside_boot = 0
  n_extra     = 0
  do i = 1,nZ/2+1
    !write(*,*)'Doing lower horizontal line ',nZ/2+1-i
    n_off = nZ/2
    call RintersectPolygon(n_wall, R_wall(1:n_wall), Z_wall(1:n_wall), Rmin, Rmax, Zlines(i), accuracy, r_min, r_max, r_min2, r_max2)
    width = r_max-r_min
    ! First jump
    if ( (width .lt. 0.85*width_prev) .and. (Zlines(i) .gt. -Z_first_jump) .and. (i .gt. 1) ) then
      diff = width/width_prev
      nR = int(real(nR) * diff) + 1
    endif
    ! Second jump
    if ( (width .gt. 1.2*width_prev) .and. (i .gt. 1) ) then
      diff = width/width_prev
      n_extra = int(real(nR) * diff) - nR
      inside_boot = 1
    endif
    nR_grid(n_off+i) = nR + n_extra
    if (inside_boot .eq. 1) then
      do j = 1,nR
        R_grid(j,n_off+i) = r_min + (R_grid(j,n_off+i-1) - R_grid(1,n_off+i-1))
        Z_grid(j,n_off+i) = Zlines(i)
      enddo
      do j = 1,n_extra
        R_grid(nR+j,n_off+i) = R_grid(nR,n_off+i) + real(j)/real(n_extra) * (r_max - R_grid(nR,n_off+i))
        Z_grid(nR+j,n_off+i) = Zlines(i)
      enddo
      i_seg = inbetweenSegment(n_wall, R_wall(1:n_wall), Z_wall(1:n_wall), r_max, Zlines(i))
      r_max = min(R_wall(i_seg),R_wall(i_seg+1))
      z_tmp = max(Z_wall(i_seg),Z_wall(i_seg+1))
      nR_grid(n_off+i-1) = nR + n_extra
      do j = 1,n_extra
        R_grid(nR+j,n_off+i-1) = R_grid(nR,n_off+i-1) + real(j)/real(n_extra) * (r_max - R_grid(nR,n_off+i-1))
        Z_grid(nR+j,n_off+i-1) = Z_grid(nR,n_off+i-1) + real(j)/real(n_extra) * (z_tmp - Z_grid(nR,n_off+i-1))
      enddo
      inside_boot = 2
    else
      do j = 1,nR_grid(n_off+i)
        R_grid(j,n_off+i) = r_min + real(j-1)/real(nR_grid(n_off+i)-1) * (r_max - r_min)
        Z_grid(j,n_off+i) = Zlines(i)
      enddo
    endif
    width_prev = width
  enddo
  write(*,*)'Finished lower part of grid'
  
  ! --- Then upper part
  nR = nR_save
  do i=1,nZ/2+1
    Zlines(i) = real(i-1)/real(nZ/2) * Zmax
  enddo
  width_prev  = 0.0
  inside_boot = 0
  n_extra     = 0
  do i = 1,nZ/2
    n_off = nZ/2+1
    !write(*,*)'Doint upper horizontal line ',nZ/2-i+1
    call RintersectPolygon(n_wall, R_wall(1:n_wall), Z_wall(1:n_wall), Rmin, Rmax, Zlines(i+1), accuracy, r_min, r_max, r_min2, r_max2)
    width = r_max-r_min
    ! First jump
    if ( (width .lt. 0.85*width_prev) .and. (Zlines(i+1) .lt. Z_first_jump) .and. (i .gt. 1) ) then
      diff = width/width_prev
      nR = int(real(nR) * diff) + 1
    endif
    ! Second jump
    if ( (width .gt. 1.2*width_prev) .and. (i .gt. 1) ) then
      diff = width/width_prev
      n_extra = int(real(nR) * diff) - nR
      inside_boot = 1
    endif
    nR_grid(n_off-i) = nR + n_extra
    if (inside_boot .eq. 1) then
      do j = 1,nR
        R_grid(j,n_off-i) = r_min + (R_grid(j,n_off-i+1)-R_grid(1,n_off-i+1))
        Z_grid(j,n_off-i) = Zlines(i+1)
      enddo
      do j = 1,n_extra
        R_grid(nR+j,n_off-i) = R_grid(nR,n_off-i) + real(j)/real(n_extra) * (r_max - R_grid(nR,n_off-i))
        Z_grid(nR+j,n_off-i) = Zlines(i+1)
      enddo
      i_seg = inbetweenSegment(n_wall, R_wall(1:n_wall), Z_wall(1:n_wall), r_max, Zlines(i+1))
      r_max = min(R_wall(i_seg),R_wall(i_seg+1))
      z_tmp = min(Z_wall(i_seg),Z_wall(i_seg+1))
      nR_grid(n_off-i+1) = nR + n_extra
      do j = 1,n_extra
        R_grid(nR+j,n_off-i+1) = R_grid(nR,n_off-i+1) + real(j)/real(n_extra) * (r_max - R_grid(nR,n_off-i+1))
        Z_grid(nR+j,n_off-i+1) = Z_grid(nR,n_off-i+1) + real(j)/real(n_extra) * (z_tmp - Z_grid(nR,n_off-i+1))
      enddo
      inside_boot = 2
    else
      do j = 1,nR_grid(n_off-i)
        R_grid(j,n_off-i) = r_min + real(j-1)/real(nR_grid(n_off-i)-1) * (r_max - r_min)
        Z_grid(j,n_off-i) = Zlines(i+1)
      enddo
    endif
    width_prev = width
  enddo
  write(*,*)'Finished upper part of grid'
  
  ! --- Print to file?
  if (.true.) then
    open(21,file='RZ_grid_inside_wall.txt')
    do i = 1,nZ
      !write(*,*)'Doint horizontal line ',nZ-i
      nR_tmp = min(nR_grid(i),nR_grid(i+1))
      do j = 1,nR_tmp-1
        write(21,'(e16.9,A,e16.9,A,e16.9,A,e16.9,A,e16.9,A,e16.9,A,e16.9,A,e16.9)') &
          R_grid(j,i),     ' ',Z_grid(j,i),    ' ', &
          R_grid(j,i+1),   ' ',Z_grid(j,i+1),  ' ', &
          R_grid(j+1,i+1), ' ',Z_grid(j+1,i+1),' ', &
          R_grid(j+1,i),   ' ',Z_grid(j+1,i)
      enddo
    enddo
    close(21)
  endif
    
  ! --- nodes indexing
  elm_count = 0
  do i = 1,nZ
    nR_tmp = min(nR_grid(i),nR_grid(i+1))
    do j = 1,nR_tmp-1
      elm_count = elm_count + 1
      ! R-index
      node_index(elm_count,1,1) = j
      node_index(elm_count,2,1) = j+1
      node_index(elm_count,3,1) = j+1
      node_index(elm_count,4,1) = j
      ! Z-index
      node_index(elm_count,1,2) = i
      node_index(elm_count,2,2) = i
      node_index(elm_count,3,2) = i+1
      node_index(elm_count,4,2) = i+1
    enddo
  enddo

  n_elm = elm_count

  return

end subroutine create_grid_inside_wall_MASTU

subroutine create_grid_inside_wall_STEP(nR, nZ, nR_grid, node_index, R_grid, Z_grid, n_elm)
  ! --- Create an r/z grid for STEP
  ! --- STEP has extended inner and outer divertors
  ! --- which are not handled adequately by create_grid_inside_wall_usual
  ! --- Assume that the outer divertors are larger than the inner divertors
  !
  !           ___
  !   ___    /  /
  !   \  \  /  /
  !    \  \/  /
  !    |      |
  !    |      |   ! A rather poor ASCII art cartoon of STEP
  !    |      |
  !    /  /\  \
  !   /__/  \  \
  !          \__\
  !

  use grid_xpoint_data, only: n_wall, R_wall, Z_wall

  implicit none

  ! --- Input variables
  integer, intent(inout)  :: nR, nZ, n_elm
  integer, intent(inout)  :: nR_grid(nZ+1,2),node_index(4*nR*nZ,4,2)
  real*8,  intent(inout)  :: R_grid(4*nR,nZ+1),Z_grid(4*nR,nZ+1)

  ! --- Local variables
  real*8  :: Rmin, Rmax, r_min1, r_max1, r_min2, r_max2, width1, width2
  real*8  :: lower_inner_z, lower_outer_z, upper_inner_z, upper_outer_z
  real*8  :: lower_split_z, lower_split_inner_r, lower_split_r, lower_split_outer_r
  real*8  :: upper_split_z, upper_split_inner_r, upper_split_r, upper_split_outer_r
  integer :: nR_lower_inner, nR_lower_outer, nR_upper_inner, nR_upper_outer

  integer :: node_lo_start, node_lo_end
  integer :: node_li_start, node_li_end
  integer :: node_core_start, node_core_end
  integer :: node_ui_start, node_ui_end
  integer :: node_uo_start, node_uo_end

  real*8  :: accuracy
  real*8, allocatable  :: Zlines(:)

  integer :: in_section, i_z, j_r, elm_count
  logical :: debug

  ! --- Initialize local variables
  debug = .false.

  accuracy = +1.d-5
  in_section = 0
  node_uo_end = nZ

  allocate(Zlines(nZ+1))

  write(*,*)'Building grid inside wall for STEP'

  ! --- Cut the fomain with horizontal lines
  Rmin = minval(R_wall(1:n_wall)) - 1.e-3  ! slightly outside
  Rmax = maxval(R_wall(1:n_wall)) + 1.e-3  ! slightly outside

  ! --- Find lower and upper divertor split points
  ! --- Where does the core seperate into the inner/outer divertors
  call leg_split_location_step(.true., lower_split_z, lower_split_inner_r, lower_split_r, lower_split_outer_r)
  call leg_split_location_step(.false., upper_split_z, upper_split_inner_r, upper_split_r, upper_split_outer_r)

  ! --- Find lower and upper divertor z values
  call find_divertor_z_values_step(.true., lower_split_z, lower_inner_z, lower_outer_z)
  call find_divertor_z_values_step(.false., upper_split_z, upper_inner_z, upper_outer_z)

  ! --- Determine ZLines, based on z values
  call determine_zlines_step(nZ, Zlines, &
                             lower_outer_z, lower_inner_z, lower_split_z, &
                             upper_split_z, upper_inner_z, upper_outer_z)

  ! --- Determine nR inner/outer at the split points
  call determine_nr_split(nR, nR_lower_inner, nR_lower_outer, &
                          lower_split_inner_r, lower_split_r, lower_split_outer_r)
  call determine_nr_split(nR, nR_upper_inner, nR_upper_outer, &
                          upper_split_inner_r, upper_split_r, upper_split_outer_r)

  if (debug) then
    write(*, '(A, 3F7.3)')'Lower Z split, inner, outer = ',lower_split_z, lower_inner_z, lower_outer_z
    write(*, '(A, 2I5)')'Lower split nR inner, outer = ',nR_lower_inner,nR_lower_outer

    write(*, '(A, 3F7.3)')'Upper Z split, inner, outer = ',upper_split_z, upper_inner_z, upper_outer_z
    write(*, '(A, 2I5)')'Upper split nR inner, outer = ',nR_upper_inner,nR_upper_outer
  end if


  ! --- Loop over Zlines
  loop_zlines: do i_z = 1, nZ + 1

    ! --- intersect the wall polygon
    call RintersectPolygon(n_wall, R_wall(1:n_wall), Z_wall(1:n_wall), Rmin, Rmax, &
                           Zlines(i_z), accuracy, r_min1, r_max1, r_min2, r_max2)
    width1 = r_max1 - r_min1
    width2 = r_max2 - r_min2

    ! --- Which section are we in ?
    ! --- Section 1 = Lower outer divertor leg
    ! --- Section 2 = Lower inner + outer divertor legs
    ! --- Section 3 = Core
    ! --- Section 4 = Upper inner + outer divertor legs
    ! --- Section 5 = Upper outer divertor leg
    !
    ! --- Figure out start/end i_z values for each section
    ! --- Used later for creating the nodes

    single_intersection_polygon: if ((width1 > 0.d0) .and. (width2 == 0.d0)) then
      if (in_section == 0) then
        in_section = 1
        node_lo_start = i_z
      end if
      if (in_section == 2) then
        in_section = 3
        node_li_end = i_z - 2
        node_lo_end = i_z - 2
        node_core_start = i_z - 1
      end if
      if (in_section == 4) then
        in_section = 5
        node_ui_end = i_z - 2
      end if
    end if single_intersection_polygon

    double_intersection_polygon: if ((width1 > 0.d0) .and. (width2 > 0.d0)) then
      if (in_section == 1) then
        in_section = 2
        node_li_start = i_z
      end if
      if (in_section == 3) then
        in_section = 4
        node_core_end = i_z - 1
        node_ui_start = i_z
        node_uo_start = i_z
      end if
    end if double_intersection_polygon

    ! --- Set values for the R_grid and Z_grid for each section

    ! --- Section 1 - Lower outer divertor leg
    section_1: if (in_section == 1) then
      nR_grid(i_z, 1) = 0
      nR_grid(i_z, 2) = nR_lower_outer

      s1_outer: do j_r = 1, nR_lower_outer
        R_grid(j_r, i_z) = r_min1 + real(j_r - 1)/real(nR_lower_outer - 1) * (r_max1 - r_min1)
        Z_grid(j_r, i_z) = Zlines(i_z)
      end do s1_outer
    end if section_1

    ! --- Section 2 - Lower inner + outer divertor legs
    section_2: if (in_section == 2) then
      nR_grid(i_z, 1) = nR_lower_inner
      nR_grid(i_z, 2) = nR_lower_outer

      s2_inner: do j_r = 1, nR_lower_inner
        R_grid(j_r, i_z) = r_min1 + real(j_r - 1)/real(nR_lower_inner - 1) * (r_max1 - r_min1)
        Z_grid(j_r, i_z) = Zlines(i_z)
      end do s2_inner

      s2_outer: do j_r = nR_lower_inner + 1, nR_lower_inner + nR_lower_outer
        R_grid(j_r, i_z) = r_min2 + real(j_r - nR_lower_inner - 1)/real(nR_lower_outer - 1) * (r_max2 - r_min2)
        Z_grid(j_r, i_z) = Zlines(i_z)
      end do s2_outer
    end if section_2

    ! --- Section 3 - Core
    section_3: if (in_section == 3) then
      nR_grid(i_z, 1) = nR
      nR_grid(i_z, 2) = 0

      s3_core: do j_r = 1, nR
        R_grid(j_r, i_z) = r_min1 + real(j_r - 1)/real(nR - 1) * (r_max1 - r_min1)
        Z_grid(j_r, i_z) = Zlines(i_z)
      end do s3_core
    end if section_3

    ! --- Section 4 - Upper inner + outer divertor legs
    section_4: if (in_section == 4) then
      nR_grid(i_z, 1) = nR_upper_inner
      nR_grid(i_z, 2) = nR_upper_outer

      s4_inner: do j_r = 1, nR_upper_inner
        R_grid(j_r, i_z) = r_min1 + real(j_r - 1)/real(nR_upper_inner - 1) * (r_max1 - r_min1)
        Z_grid(j_r, i_z) = Zlines(i_z)
      end do s4_inner

      s4_outer: do j_r =  nR_upper_inner + 1, nR_upper_inner + nR_upper_outer
        R_grid(j_r, i_z) = r_min2 + real(j_r - nR_upper_inner - 1)/real(nR_upper_outer - 1) * (r_max2 - r_min2)
        Z_grid(j_r, i_z) = Zlines(i_z)
      end do s4_outer
    end if section_4

    ! --- Section 5 - Upper outer divertor leg
    section_5: if (in_section == 5) then
      nR_grid(i_z, 1) = 0
      nR_grid(i_z, 2) = nR_upper_outer

      s5_outer: do j_r = 1, nR_upper_outer
        R_grid(j_r, i_z) = r_min1 + real(j_r - 1)/real(nR_upper_outer - 1) * (r_max1 - r_min1)
        Z_grid(j_r, i_z) = Zlines(i_z)
      end do s5_outer
    end if section_5

  end do loop_zlines

  ! --- Create nodes for each section
  elm_count = 0

  ! --- lower inner divertor nodes
  lower_inner_divertor_nodes: do i_z = node_li_start, node_li_end
    do j_r = 1, nR_grid(i_z, 1) - 1
      elm_count = elm_count + 1
      ! -- Lower left corner
      node_index(elm_count, 1, 1) = j_r
      node_index(elm_count, 1, 2) = i_z
      ! -- Lower right corner
      node_index(elm_count, 2, 1) = j_r + 1
      node_index(elm_count, 2, 2) = i_z
      ! -- Upper right corner
      node_index(elm_count, 3, 1) = j_r + 1
      node_index(elm_count, 3, 2) = i_z + 1
      ! -- Upper left corner
      node_index(elm_count, 4, 1) = j_r
      node_index(elm_count, 4, 2) = i_z + 1

    end do
  end do lower_inner_divertor_nodes

  ! --- lower outer divertor nodes
  lower_outer_divertor_nodes: do i_z = node_lo_start, node_lo_end
    do j_r = nR_grid(i_z, 1) + 1, nR_grid(i_z, 1) + nR_grid(i_z, 2) - 1
      elm_count = elm_count + 1

      ! -- Lower left corner
      node_index(elm_count, 1, 1) = j_r
      node_index(elm_count, 1, 2) = i_z
      ! -- Lower right corner
      node_index(elm_count, 2, 1) = j_r + 1
      node_index(elm_count, 2, 2) = i_z
      ! -- Upper right corner
      if (i_z == node_li_start - 1)  then
        node_index(elm_count, 3, 1) = j_r + 1 + nR_grid(i_z+1, 1)
      else
        node_index(elm_count, 3, 1) = j_r + 1
      end if
      node_index(elm_count, 3, 2) = i_z + 1
      ! -- Upper left corner
      if (i_z == node_li_start - 1)  then
        node_index(elm_count, 4, 1) = j_r + nR_grid(i_z+1, 1)
      else
        node_index(elm_count, 4, 1) = j_r
      end if
      node_index(elm_count, 4, 2) = i_z + 1

    end do
  end do lower_outer_divertor_nodes

  ! --- Core nodes
  core_nodes: do i_z = node_core_start, node_core_end
    do j_r = 1, nR_grid(i_z, 1) + nR_grid(i_z, 2) - 1
      elm_count = elm_count + 1
      ! -- Lower left corner
      node_index(elm_count, 1, 1) = j_r
      node_index(elm_count, 1, 2) = i_z
      ! -- Lower right corner
      node_index(elm_count, 2, 1) = j_r + 1
      node_index(elm_count, 2, 2) = i_z
      ! -- Upper right corner
      node_index(elm_count, 3, 1) = j_r + 1
      node_index(elm_count, 3, 2) = i_z + 1
      ! -- Upper left corner
      node_index(elm_count, 4, 1) = j_r
      node_index(elm_count, 4, 2) = i_z + 1
    end do
  end do core_nodes

  ! --- upper inner divertor nodes
  upper_inner_nodes: do i_z = node_ui_start, node_ui_end
    do j_r = 1, nR_grid(i_z, 1) - 1
      elm_count = elm_count + 1
      ! -- Lower left corner
      node_index(elm_count, 1, 1) = j_r
      node_index(elm_count, 1, 2) = i_z
      ! -- Lower right corner
      node_index(elm_count, 2, 1) = j_r + 1
      node_index(elm_count, 2, 2) = i_z
      ! -- Upper right corner
      node_index(elm_count, 3, 1) = j_r + 1
      node_index(elm_count, 3, 2) = i_z + 1
      ! -- Upper left corner
      node_index(elm_count, 4, 1) = j_r
      node_index(elm_count, 4, 2) = i_z + 1

    end do
  end do upper_inner_nodes

  ! --- upper outer divertor nodes
  upper_outer_nodes: do i_z = node_uo_start, node_uo_end
    do j_r = nR_grid(i_z, 1) + 1, nR_grid(i_z, 1) + nR_grid(i_z, 2) - 1
      elm_count = elm_count + 1
      ! -- Lower left corner
      node_index(elm_count, 1, 1) = j_r
      node_index(elm_count, 1, 2) = i_z
      ! -- Lower right corner
      node_index(elm_count, 2, 1) = j_r + 1
      node_index(elm_count, 2, 2) = i_z
      ! -- Upper right corner
      if (i_z == node_ui_end + 1)  then
        node_index(elm_count, 3, 1) = j_r + 1 - nR_grid(i_z, 1)
      else
        node_index(elm_count, 3, 1) = j_r + 1
      end if
      node_index(elm_count, 3, 2) = i_z + 1
      ! -- Upper left corner
      if (i_z == node_ui_end + 1)  then
        node_index(elm_count, 4, 1) = j_r - nR_grid(i_z, 1)
      else
        node_index(elm_count, 4, 1) = j_r
      end if
      node_index(elm_count, 4, 2) = i_z + 1
    end do
  end do upper_outer_nodes

  n_elm = elm_count

  ! --- clean up
  deallocate(Zlines)

  return

end subroutine create_grid_inside_wall_STEP

subroutine leg_split_location_step(is_lower, common_z, inner_r, split_r, outer_r)
  ! --- Determine the leg split location
  ! --- In ASCII art below (for upper), where is (a) in r/z
  ! --- and what are the r values for the inner/outer wall at Z(a)
  !           ___
  !   ___    /  /
  !   \  \  /  /
  !    \  \/a /
  !    |      |

  use grid_xpoint_data, only: n_wall, R_wall, Z_wall

  implicit none

  ! --- Input variables
  logical, intent(in) :: is_lower
  real*8, intent(out) :: common_z, inner_r, split_r, outer_r

  ! --- Local variables
  real*8  :: Rmin, Rmax, r_min1, r_max1, r_min2, r_max2, width1, width2
  real*8  :: z_limit, z_value, accuracy
  integer :: i

  ! --- Initialize local variables
  accuracy = +1.d-5

  ! --- find the limits of the wall
  if (is_lower)        z_limit = minval(Z_wall(1:n_wall)) + 1.e-3  ! slightly inside
  if (.not. is_lower)  z_limit = maxval(Z_wall(1:n_wall)) - 1.e-3  ! slightly inside

  Rmin = minval(R_wall(1:n_wall)) - 1.e-3  ! slightly outside
  Rmax = maxval(R_wall(1:n_wall)) + 1.e-3  ! slightly outside

  ! --- Loop over Z with a very high resolution
  ! --- break loop once there are 2 intersections of the wall polygon
  find_split: do i = 300, 500
    z_value = real(i-1)/real(499) * z_limit

    ! --- intersect the wall polygon
    call RintersectPolygon(n_wall, R_wall(1:n_wall), Z_wall(1:n_wall), Rmin, Rmax, &
                           z_value, accuracy, r_min1, r_max1, r_min2, r_max2)
    width1 = r_max1 - r_min1
    width2 = r_max2 - r_min2

    ! --- Has the split been found ?
    if ((width1 > 0.d0) .and. (width2 > 0.d0)) then
      ! --- split found, determine r by taking the average of max1 and min2
      split_r = (r_max1 + r_min2) / 2.d0
      exit find_split

    else
      ! --- split not found, set the outputs required prior to split
      common_z = z_value
      inner_r = r_min1
      outer_r = r_max1
    end if

  end do find_split

  ! Apply a small buffer to the common z
  if (is_lower)        common_z = common_z + 0.01
  if (.not. is_lower)  common_z = common_z - 0.01

  return

end subroutine leg_split_location_step

subroutine find_divertor_z_values_step(is_lower, split_z, inner_z, outer_z)
  ! --- Determine the Z values of the inner and outer divertors
  ! --- In ASCII art below (for upper), find (b) and (c), given (a)
  !           ___ c
  ! b ___    /  /
  !   \  \  /  /
  !    \  \/a /
  !    |      |

  use grid_xpoint_data, only: n_wall, R_wall, Z_wall

  implicit none

  ! --- Input variables
  logical, intent(in) :: is_lower
  real*8, intent(in)  :: split_z
  real*8, intent(out) :: inner_z, outer_z

  ! --- Local variables
  real*8  :: Rmin, Rmax, r_min1, r_max1, r_min2, r_max2, width1, width2
  real*8  :: accuracy
  integer :: i

  ! --- Initialize local variables
  accuracy = +1.d-5

  ! --- find the limits of the wall
  if (is_lower)        outer_z = minval(Z_wall(1:n_wall)) + 1.e-3  ! slightly inside
  if (.not. is_lower)  outer_z = maxval(Z_wall(1:n_wall)) - 1.e-3  ! slightly inside

  Rmin = minval(R_wall(1:n_wall)) - 1.e-3  ! slightly outside
  Rmax = maxval(R_wall(1:n_wall)) + 1.e-3  ! slightly outside

  ! --- Loop over Z with a very high resolutuon
  ! --- break loop once there are 2 intersections of the wall polygon
  ! --- Start from the top/bottom and work inwards

  find_inner_z: do i= 1, 400
    inner_z = outer_z + (real(i-1)/real(399) * (split_z - outer_z))

    ! --- intersect the wall polygon
    call RintersectPolygon(n_wall, R_wall(1:n_wall), Z_wall(1:n_wall), Rmin, Rmax, &
                           inner_z, accuracy, r_min1, r_max1, r_min2, r_max2)
    width1 = r_max1 - r_min1
    width2 = r_max2 - r_min2

    if ((width1 > 0.d0) .and. (width2 > 0.d0)) exit find_inner_z

  end do find_inner_z

  ! --- Shift inner points lightly inside
  if (is_lower)   inner_z = inner_z + 1.e-3  ! slightly inside
  if (.not. is_lower)  inner_z = inner_z - 1.e-3  ! slightly inside

  return

end subroutine find_divertor_z_values_step

subroutine determine_zlines_step(nZ, Zlines, &
                                 lower_outer_z, lower_inner_z, lower_split_z, &
                                 upper_split_z, upper_inner_z, upper_outer_z)
  ! --- Determine an optimal z lines, based on multiple z values
  implicit none

  ! --- Input variables
  integer, intent(in)    :: nZ
  real*8,  intent(in)    :: lower_outer_z, lower_inner_z, lower_split_z
  real*8,  intent(in)    :: upper_split_z, upper_inner_z, upper_outer_z
  real*8,  intent(inout) :: Zlines(nZ+1)

  ! --- Local variables
  logical :: debug
  real*8  :: total_height, fraction_core
  real*8  :: fraction_lower_inner, fraction_lower_outer
  real*8  :: fraction_upper_inner, fraction_upper_outer

  integer :: nZ_core, nZ_diff
  integer :: nZ_lower_inner, nZ_lower_outer
  integer :: nZ_upper_inner, nZ_upper_outer
  integer :: i
  integer :: i_lo_start, i_lo_end
  integer :: i_li_start, i_li_end
  integer :: i_core_start, i_core_end
  integer :: i_ui_start, i_ui_end
  integer :: i_uo_start, i_uo_end

  real*8  :: delta_lower_inner, delta_lower_outer, delta_core, delta_upper_inner, delta_upper_outer

  debug = .false.

  ! --- determine fractions of Z for each section
  total_height = abs(lower_outer_z) + abs(upper_split_z)

  fraction_lower_outer = (abs(lower_outer_z) - abs(lower_inner_z)) / total_height
  fraction_lower_inner = (abs(lower_inner_z) - abs(lower_split_z)) / total_height
  fraction_core = (abs(lower_split_z) + abs(upper_split_z)) / total_height
  fraction_upper_inner = (abs(upper_inner_z) - abs(upper_split_z)) / total_height
  fraction_upper_outer = (abs(upper_outer_z) - abs(upper_inner_z)) / total_height

  ! --- determine how many lines each section should have
  ! --- truncation extreemly likely here
  nZ_lower_outer = int(fraction_lower_outer * nZ)
  nZ_lower_inner = int(fraction_lower_inner * nZ)
  nZ_core = int(fraction_core * nZ)
  nZ_upper_inner = int(fraction_upper_inner * nZ)
  nZ_upper_outer = int(fraction_upper_outer * nZ)

  ! --- closure condition: sum(nZ_section) == nZ
  ! --- add/subtract any truncation diff to the core
  nZ_diff = nZ - nZ_lower_outer - nZ_lower_inner - nZ_core - nZ_upper_inner - nZ_upper_outer
  nZ_core = nZ_core + nZ_diff
  if (debug) then
    write(*,'(A, 7I5)')'nZ lo, li, core, ui, uo, total, nZ, diff = ', &
          nZ_lower_outer, nZ_lower_inner, nZ_core, nZ_upper_inner, nZ_upper_outer, &
          nZ_lower_outer + nZ_lower_inner + nZ_core + nZ_upper_inner + nZ_upper_outer, nZ,nZ_diff
  end if

  ! --- Grid spacing deltas for each section
  delta_lower_outer = (abs(lower_outer_z) - abs(lower_inner_z)) / nZ_lower_outer
  delta_lower_inner = (abs(lower_inner_z) - abs(lower_split_z)) / nZ_lower_inner
  delta_core = (abs(lower_split_z) + abs(upper_split_z)) / nZ_core
  delta_upper_inner = (abs(upper_inner_z) - abs(upper_split_z)) / nZ_upper_inner
  delta_upper_outer = (abs(upper_outer_z) - abs(upper_inner_z)) / nZ_upper_outer

  if (debug) then
    write(*, '(A, 5F7.4)')' Deltas :: lo, li, core, ui, uo = ', &
          delta_lower_outer, delta_lower_inner, delta_core, delta_upper_inner, delta_upper_outer
  end if

  ! --- loop start/end values
  i_lo_start = 1
  i_lo_end = nZ_lower_outer

  i_li_start = i_lo_end + 1
  i_li_end = i_li_start + nZ_lower_inner - 1

  i_core_start = i_li_end + 1
  i_core_end = i_core_start + nZ_core - 1

  i_ui_start = i_core_end + 1
  i_ui_end = i_ui_start + nZ_upper_inner - 1

  i_uo_start = i_ui_end + 1
  i_uo_end = i_uo_start + nZ_upper_outer

  if (debug) then
    write(*, '(A, 2I5)')' Loop start, end :: lower outer :: ',i_lo_start, i_lo_end
    write(*, '(A, 2I5)')' Loop start, end :: lower inner :: ',i_li_start, i_li_end
    write(*, '(A, 2I5)')' Loop start, end :: core :: ',i_core_start, i_core_end
    write(*, '(A, 2I5)')' Loop start, end :: upper inner :: ',i_ui_start, i_ui_end
    write(*, '(A, 2I5)')' Loop start, end :: upper outer :: ',i_uo_start, i_uo_end
  end if

  ! --- Create Z lines
  loop_lo: do i = i_lo_start, i_lo_end
    Zlines(i) = lower_outer_z + real(i - i_lo_start) * delta_lower_outer
    if (debug) write(*,'(A, i5, F7.3)')' Lower outer i, Zlines = ',i, Zlines(i)
  end do loop_lo

  loop_li: do i = i_li_start, i_li_end
    Zlines(i) = lower_inner_z + real(i - i_li_start) * delta_lower_inner
    if (debug) write(*,'(A, i5, F7.3)')' Lower inner i, Zlines = ',i, Zlines(i)
  end do loop_li

  loop_core: do i = i_core_start, i_core_end
    Zlines(i) = lower_split_z + real(i - i_core_start) * delta_core
    if (debug) write(*,'(A, i5, F7.3)')' Core i, Zlines = ',i, Zlines(i)
  end do loop_core

  loop_ui: do i = i_ui_start, i_ui_end
    Zlines(i) = upper_split_z + real(i - i_ui_start) * delta_upper_inner
    if (debug) write(*,'(A, i5, F7.3)')' Upper inner i, Zlines = ',i, Zlines(i)
  end do loop_ui

  loop_uo: do i = i_uo_start, i_uo_end
    Zlines(i) = upper_inner_z + real(i - i_uo_start) * delta_upper_outer
    if (debug) write(*,'(A, i5, F7.3)')' Upper outer i, Zlines = ',i, Zlines(i)
  end do loop_uo

  return

end subroutine determine_zlines_step

subroutine determine_nr_split(nR, nR_inner, nR_outer, inner_r, split_r, outer_r)
  ! --- How many R grid pieces should go into the inner/outer divertors
  implicit none

  ! --- Input variables
  integer, intent(in)  :: nR
  real*8,  intent(in)  :: inner_r, split_r, outer_r
  integer, intent(out) :: nR_inner, nR_outer

  ! --- Local variables
  real*8 :: fraction_inner, fraction_outer, total_width, nR_diff

  ! --- determine fractions
  total_width = outer_r - inner_r
  fraction_inner = (split_r - inner_r) / total_width
  fraction_outer = (outer_r - split_r) / total_width

  ! --- determine how many lines each section should have
  ! --- truncation extreemly likely here
  nR_inner = int(fraction_inner * nR)
  nR_outer = int(fraction_outer * nR)

  ! --- closure condition : nR_inner + nR_outer == nR
  ! --- add/subtract any trancation diff to the inner
  nR_diff = nR - nR_inner - nR_outer
  nR_inner = nR_inner + nR_diff

  return

end subroutine determine_nr_split
