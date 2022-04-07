!> This routine finds the intersections of a flux surface with the target (ie. wall below Xpoint)
subroutine get_target_flux_surfaces(node_list, element_list, surface_list, ignore, stpts, &
                                    psi_bnd, R_axis, Z_axis, R_xpoint, Z_xpoint,          &
                                    n_int_max, n_int, R_int, Z_int, index_int,            &
                                    n_int_surf, index_int_surf, ifail)


  use data_structure
  use phys_module
  use high_resolution_wall
  use constants
  use grid_xpoint_data
  use py_plots_grids
  use mod_interp, only: interp_RZ
  implicit none
  
  ! --- Routine parameters
  type (type_node_list),        intent(in)      :: node_list
  type (type_element_list),     intent(in)      :: element_list
  type (type_surface_list),     intent(inout)   :: surface_list
  type (type_strategic_points), intent(in)      :: stpts
  logical,                      intent(in)      :: ignore(surface_list%n_psi)
  integer,                      intent(in)      :: n_int_max
  integer,                      intent(inout)   :: n_int, ifail
  real*8,                       intent(inout)   :: R_int(n_int_max), Z_int(n_int_max)
  real*8,                       intent(inout)   :: R_axis,           Z_axis
  real*8,                       intent(inout)   :: R_xpoint(2),      Z_xpoint(2)
  real*8,                       intent(inout)   :: psi_bnd
  integer,                      intent(inout)   :: index_int(n_int_max,4) ! 1->surface_index
                                                                          ! 2->surface_piece_index
                                                                          ! 3->wall_piece_index
                                                                          ! 4->target index (1:LowerLeft - 2:LowerRight - 3:UpperLeft - 4:UpperRight)
  integer,                      intent(inout)   :: n_int_surf(10*surface_list%n_psi) ! number of intersections for each surface
  integer,                      intent(inout)   :: index_int_surf(10*surface_list%n_psi,n_int_max) ! index of intersections on surface
  
  ! --- Local variables
  integer                       :: i, j, i_surf, i_int
  real*8                        :: theta_x(2), Z_line(2)
  integer                       :: n_int_tmp, index_int_tmp(n_int_max,4)
  real*8                        :: R_int_tmp(n_int_max), Z_int_tmp(n_int_max)
  integer, parameter            :: LowerLeft =1
  integer, parameter            :: LowerRight=2
  integer, parameter            :: UpperLeft =3
  integer, parameter            :: UpperRight=4
  integer                       :: index_tmp
  integer                       :: max_LowerLeft, max_LowerRight, max_UpperLeft, max_UpperRight
  real*8                        :: Z_save(4)
  integer                       :: i_save(4)
  logical                       :: debug, save_point(4), save_int, target_only
  character*256                 :: filename
  
  write(*,*) '***********************************'
  write(*,*) '*     get_target_flux_surfaces    *'
  write(*,*) '***********************************'
  
  ! --- Initialise some values
  n_int = 0
  debug = .true.
  target_only = .true.
  
  ! --- First get all intersections with wall
  do i_surf = 1,surface_list%n_psi
    if (ignore(i_surf)) cycle
    call find_wall_crossings_with_flux_surface(node_list, element_list, surface_list%flux_surfaces(i_surf), target_only, &
                                               n_int_max, n_int_tmp, R_int_tmp, Z_int_tmp, index_int_tmp, ifail)

    if (ifail .ne. 0) then
      write(*,*) 'Warning! Failed to find all wall intersections for surface',i_surf,ifail
      return
    endif
    
    ! --- Fill up target arrays
    do i=1,n_int_tmp
      R_int(n_int + i) = R_int_tmp(i)
      Z_int(n_int + i) = Z_int_tmp(i)
      index_int(n_int + i, 1) = i_surf
      index_int(n_int + i, 2) = index_int_tmp(i,2)
      index_int(n_int + i, 3) = index_int_tmp(i,3)
    enddo
    n_int = n_int + n_int_tmp
  enddo
  
  ! --- Some debug plots
  if (debug) then
    write(*,*)'number of total intersections found:',n_int 
    filename = 'plot_wall_intersections.py'
    call print_py_plot_prepare_plot(filename)
    call print_py_plot_ordered_flux_surfaces(filename, node_list, element_list, surface_list, 'r', .false.)
    call print_py_plot_points(filename,n_int,R_int,Z_int)
    call print_py_plot_wall(filename)
    call print_py_plot_finish_plot(filename)
  endif

  ! --- Copy back into temporary array
  do i=1,n_int
    R_int_tmp(i) = R_int(i)
    Z_int_tmp(i) = Z_int(i)
    index_int_tmp(i,1) = index_int(i,1)
    index_int_tmp(i,2) = index_int(i,2)
    index_int_tmp(i,3) = index_int(i,3)
  enddo
  n_int_tmp = n_int
  
  ! --- Determine angle of X-line
  if (xcase .ne. UPPER_XPOINT) theta_x(1) = atan2(Z_xpoint(1)-Z_axis,R_xpoint(1)-R_axis) + 0.5d0*PI
  if (xcase .ne. LOWER_XPOINT) theta_x(2) = atan2(Z_xpoint(2)-Z_axis,R_xpoint(2)-R_axis) + 1.5d0*PI
  if (theta_x(1) .lt. 0.d0)    theta_x(1) = theta_x(1) + 2.d0*PI
  if (theta_x(2) .lt. 0.d0)    theta_x(2) = theta_x(2) + 2.d0*PI
  if (theta_x(1) .gt. 2.d0*PI) theta_x(1) = theta_x(1) - 2.d0*PI
  if (theta_x(2) .gt. 2.d0*PI) theta_x(2) = theta_x(2) - 2.d0*PI
  
  ! --- Now extract target points. First, retain only points below Xpoint-line (or above 2nd Xpoint's)
  n_int = 0
  do i=1,n_int_tmp
    save_int = .false.
    if (xcase .ne. UPPER_XPOINT) Z_line(1) = Z_xpoint(1) + (R_int_tmp(i) - R_xpoint(1)) * tan(theta_x(1))
    if (xcase .ne. LOWER_XPOINT) Z_line(2) = Z_xpoint(2) + (R_int_tmp(i) - R_xpoint(2)) * tan(theta_x(2))
    ! --- Usual case
    if (     ((xcase .ne. UPPER_XPOINT) .and. (Z_int_tmp(i) .lt. Z_line(1)))       &
        .or. ((xcase .ne. LOWER_XPOINT) .and. (Z_int_tmp(i) .gt. Z_line(2)))        ) save_int = .true.
    ! --- If our legs go off to the side under Xpoint
    if (     ((xcase .ne. UPPER_XPOINT) .and. (R_int_tmp(i) .gt. stpts%RLimit_LowerOuterLeg)               &
                                        .and. (Z_int_tmp(i) .lt. stpts%ZLimit_LowerOuterLeg))              &
        .or. ((xcase .ne. LOWER_XPOINT) .and. (R_int_tmp(i) .gt. stpts%RLimit_UpperOuterLeg)               &
                                        .and. (Z_int_tmp(i) .gt. stpts%ZLimit_UpperOuterLeg))               ) save_int = .true.
    ! --- MAST special
    if (tokamak_device(1:4) .eq. 'MAST') then
      if (     ((xcase .ne. UPPER_XPOINT) .and. (R_int_tmp(i) .gt. stpts%RLimit_LowerMastWallBox)          &
                                          .and. (Z_int_tmp(i) .lt. stpts%ZLimit_LowerMastWallBox))         &
          .or. ((xcase .ne. LOWER_XPOINT) .and. (R_int_tmp(i) .gt. stpts%RLimit_UpperMastWallBox)          &
                                          .and. (Z_int_tmp(i) .gt. stpts%ZLimit_UpperMastWallBox))          ) save_int = .true.
    endif
    if (save_int) then
      n_int = n_int + 1
      R_int(n_int) = R_int_tmp(i)
      Z_int(n_int) = Z_int_tmp(i)
      index_int(n_int,1) = index_int_tmp(i,1)
      index_int(n_int,2) = index_int_tmp(i,2)
      index_int(n_int,3) = index_int_tmp(i,3)
      if (Z_int(n_int) .lt. Z_axis) then
        if (R_int(n_int) .lt. R_xpoint(1)) then
          index_int(n_int,4) = LowerLeft
        else
          index_int(n_int,4) = LowerRight
        endif
      else
        if (R_int(n_int) .lt. R_xpoint(2)) then
          index_int(n_int,4) = UpperLeft
        else
          index_int(n_int,4) = UpperRight
        endif
      endif
    endif
  enddo

  ! --- Copy back into temporary array
  do i=1,n_int
    R_int_tmp(i) = R_int(i)
    Z_int_tmp(i) = Z_int(i)
    index_int_tmp(i,1) = index_int(i,1)
    index_int_tmp(i,2) = index_int(i,2)
    index_int_tmp(i,3) = index_int(i,3)
    index_int_tmp(i,4) = index_int(i,4)
  enddo
  n_int_tmp = n_int
    
  ! --- Fill up surface arrays
  do i_surf = 1,surface_list%n_psi
    n_int_surf(i_surf) = 0
  enddo
  do i=1,n_int_tmp
    i_surf = index_int_tmp(i,1)
    n_int_surf(i_surf) = n_int_surf(i_surf) + 1
    index_int_surf(i_surf,n_int_surf(i_surf)) = i
  enddo
  
  ! --- Loop on all surfaces to get extrema intersections only (ie. target points)
  n_int = 0
  do i_surf=1,surface_list%n_psi
    if (ignore(i_surf)) cycle
    max_LowerLeft  = 10000
    max_LowerRight = 10000
    max_UpperLeft  = 10000
    max_UpperRight = 10000
    save_point(LowerLeft)  = .false.
    save_point(LowerRight) = .false.
    save_point(UpperLeft)  = .false.
    save_point(UpperRight) = .false.
    do i=1,n_int_surf(i_surf)
      ! --- Get the index of the intersection on the surface and its distance from the surface ends
      ! --- Note, by definition, if you want to build your grid, you need all your SOL surfaces to be continuous
      ! --- from one target to another (ie. one surface part), so we assume there is only one part to the surface...
      i_int = index_int_surf(i_surf,i)
      index_tmp = min(index_int_tmp(i_int,2), surface_list%flux_surfaces(i_surf)%n_pieces - index_int_tmp(i_int,2))
      if ( (index_int_tmp(i_int,4) .eq. LowerLeft)  .and. (index_tmp .le. max_LowerLeft) ) then
        if ( (index_tmp .eq. max_LowerLeft)  .and. (Z_int_tmp(i_int) .gt. Z_save(LowerLeft)) )  cycle ! There can be many intersections on one piece...
        max_LowerLeft = index_tmp
        save_point(LowerLeft) = .true.
        Z_save(LowerLeft) = Z_int_tmp(i_int)
        i_save(LowerLeft) = i
      endif
      if ( (index_int_tmp(i_int,4) .eq. LowerRight) .and. (index_tmp .le. max_LowerRight) ) then
        if ( (index_tmp .eq. max_LowerRight) .and. (Z_int_tmp(i_int) .gt. Z_save(LowerRight)) ) cycle ! There can be many intersections on one piece...
        max_LowerRight = index_tmp
        save_point(LowerRight) = .true.
        Z_save(LowerRight) = Z_int_tmp(i_int)
        i_save(LowerRight) = i
      endif
      if ( (index_int_tmp(i_int,4) .eq. UpperLeft)  .and. (index_tmp .le. max_UpperLeft) ) then
        if ( (index_tmp .eq. max_UpperLeft)  .and. (Z_int_tmp(i_int) .lt. Z_save(UpperLeft)) )  cycle ! There can be many intersections on one piece...
        max_UpperLeft = index_tmp
        save_point(UpperLeft) = .true.
        Z_save(UpperLeft) = Z_int_tmp(i_int)
        i_save(UpperLeft) = i
      endif
      if ( (index_int_tmp(i_int,4) .eq. UpperRight) .and. (index_tmp .le. max_UpperRight) ) then
        if ( (index_tmp .eq. max_UpperRight) .and. (Z_int_tmp(i_int) .lt. Z_save(UpperRight)) ) cycle ! There can be many intersections on one piece...
        max_UpperRight = index_tmp
        save_point(UpperRight) = .true.
        Z_save(UpperRight) = Z_int_tmp(i_int)
        i_save(UpperRight) = i
      endif
    enddo
    ! --- Now that we found extrema, save points
    do j=1,4
      if (save_point(j)) then
        n_int = n_int + 1
        i_int = index_int_surf(i_surf,i_save(j))
        R_int(n_int) = R_int_tmp(i_int)
        Z_int(n_int) = Z_int_tmp(i_int)
        index_int(n_int,1) = index_int_tmp(i_int,1)
        index_int(n_int,2) = index_int_tmp(i_int,2)
        index_int(n_int,3) = index_int_tmp(i_int,3)
        index_int(n_int,4) = index_int_tmp(i_int,4)
      endif
    enddo
  enddo
  
  ! --- Fill up surface arrays again
  do i_surf = 1,surface_list%n_psi
    n_int_surf(i_surf) = 0
  enddo
  do i=1,n_int
    i_surf = index_int(i,1)
    n_int_surf(i_surf) = n_int_surf(i_surf) + 1
    index_int_surf(i_surf,n_int_surf(i_surf)) = i
  enddo
  
  ! --- Some debug plots
  if (debug) then
    write(*,*)'number of total target intersections found:',n_int 
    filename = 'plot_target_intersections.py'
    call print_py_plot_prepare_plot(filename)
    call print_py_plot_ordered_flux_surfaces(filename, node_list, element_list, surface_list, 'r', .false.)
    call print_py_plot_points(filename,n_int,R_int,Z_int)
    call print_py_plot_wall(filename)
    call print_py_plot_finish_plot(filename)
  endif


end subroutine get_target_flux_surfaces
























!> This routine finds all intersections between a flux surface and the wall
subroutine find_wall_crossings_with_flux_surface(node_list, element_list, surface, target_only, n_int_max, n_int, R_int, Z_int, index_int, ifail)

  use data_structure
  use phys_module
  use high_resolution_wall
  use mod_interp, only: interp_RZ
  use check_point_is_inside_wall_contour
  implicit none
  
  ! --- Routine parameters
  type (type_node_list),        intent(in)      :: node_list
  type (type_element_list),     intent(in)      :: element_list
  type (type_surface),          intent(inout)   :: surface
  logical,                      intent(in)      :: target_only
  integer,                      intent(in)      :: n_int_max
  integer,                      intent(inout)   :: n_int, ifail
  real*8,                       intent(inout)   :: R_int(n_int_max), Z_int(n_int_max)
  integer,                      intent(inout)   :: index_int(n_int_max,3)
  
  ! --- Local variables
  integer                       :: i,j,i_elm,k,l
  integer                       :: istart,iend
  real*8                        :: ss,ss_int, dl_int
  real*8                        :: tt,tt_int
  real*8                        :: R_tmp,Rbeg,Rend,Rmin,Rmax
  real*8                        :: Z_tmp,Zbeg,Zend,Zmin,Zmax
  real*8                        :: R_cub1d(4)
  real*8                        :: Z_cub1d(4)
  real*8                        :: distance
  real*8                        :: accuracy
  integer                       :: n_int_tmp
  real*8                        :: R_int_tmp(n_int_max), Z_int_tmp(n_int_max)
  integer                       :: index_int_tmp(n_int_max,3)
  integer                       :: inside(4)
  logical                       :: piece_duplicated, debug
  
  ! --- Initialise some data
  n_int = 0
  accuracy = 1.d-4
  debug = .false.
  
  ! --- Step along surface
  do j=1,surface%n_pieces
    ! --- First check if piece is inside wall
    i_elm = surface%elm(j)
    
    ! --- Get approximate length of piece
    ss = surface%s(1,j)
    tt = surface%t(1,j)
    call interp_RZ(node_list,element_list,i_elm,ss,tt,Rbeg,Zbeg)

    ss = surface%s(3,j)
    tt = surface%t(3,j)
    call interp_RZ(node_list,element_list,i_elm,ss,tt,Rend,Zend)
    
    Rmin = min(Rbeg,Rend) - 0.02d0 ! extra 2cm just to make sure
    Rmax = max(Rbeg,Rend) + 0.02d0 ! extra 2cm just to make sure
    Zmin = min(Zbeg,Zend) - 0.02d0 ! extra 2cm just to make sure
    Zmax = max(Zbeg,Zend) + 0.02d0 ! extra 2cm just to make sure
    call check_point_is_inside_wall(Rmin, Zmin, inside(1))
    call check_point_is_inside_wall(Rmin, Zmax, inside(2))
    call check_point_is_inside_wall(Rmax, Zmax, inside(3))
    call check_point_is_inside_wall(Rmax, Zmin, inside(4))
    if ( (inside(1) .eq. 1) .and. (inside(2) .eq. 1) .and. (inside(3) .eq. 1) .and. (inside(4) .eq. 1) ) cycle
    
    ! --- Find any intersection with wall
    if (target_only) then
      istart = first_target_point
      iend   = last_target_point-1
      if (iend .lt. istart)    iend = n_limiter-1
      if (iend .eq. n_limiter) iend = n_limiter-1
    else
      istart = 1
      iend   = n_limiter-1
    endif
    do l=istart,iend
      R_cub1d(1) = R_limiter(l    )               ;       Z_cub1d(1) = Z_limiter(l    )
      R_cub1d(3) = R_limiter(l + 1)               ;       Z_cub1d(3) = Z_limiter(l + 1)
      R_cub1d(2) = (R_cub1d(3)-R_cub1d(1))/2.d0   ;       Z_cub1d(2) = (Z_cub1d(3)-Z_cub1d(1))/2.d0
      R_cub1d(4) = (R_cub1d(3)-R_cub1d(1))/2.d0   ;       Z_cub1d(4) = (Z_cub1d(3)-Z_cub1d(1))/2.d0
      call find_crossing_on_surface_piece(node_list,element_list,surface,j,R_cub1d,Z_cub1d, &
                                          R_tmp,Z_tmp,ss_int,tt_int,dl_int,ifail,.false.)
      if (ifail .eq. 0) then
        n_int = n_int + 1
        index_int(n_int,2) = j
        index_int(n_int,3) = l
        R_int(n_int)       = R_tmp
        Z_int(n_int)       = Z_tmp
      endif
    enddo
    if ( (target_only) .and. (last_target_point .lt. first_target_point) ) then
      istart = 1
      iend   = last_target_point-1
      do l=istart,iend
        R_cub1d(1) = R_limiter(l    )               ;       Z_cub1d(1) = Z_limiter(l    )
        R_cub1d(3) = R_limiter(l + 1)               ;       Z_cub1d(3) = Z_limiter(l + 1)
        R_cub1d(2) = (R_cub1d(3)-R_cub1d(1))/2.d0   ;       Z_cub1d(2) = (Z_cub1d(3)-Z_cub1d(1))/2.d0
        R_cub1d(4) = (R_cub1d(3)-R_cub1d(1))/2.d0   ;       Z_cub1d(4) = (Z_cub1d(3)-Z_cub1d(1))/2.d0
        call find_crossing_on_surface_piece(node_list,element_list,surface,j,R_cub1d,Z_cub1d, &
                                            R_tmp,Z_tmp,ss_int,tt_int,dl_int,ifail,.false.)
        if (ifail .eq. 0) then
          n_int = n_int + 1
          index_int(n_int,2) = j
          index_int(n_int,3) = l
          R_int(n_int)       = R_tmp
          Z_int(n_int)       = Z_tmp
        endif
      enddo
    endif
  enddo
  
  ! --- Sometimes points are duplicated (if they lie almost exactly on the edge of two neighbour pieces). Remove them...
  n_int_tmp = 0
  do i=1,n_int-1
    ! --- Loop over all other pieces to check if one is the same
    piece_duplicated = .false.
    do j=i+1,n_int
      distance   = sqrt( (R_int(i)-R_int(j))**2.d0 + (Z_int(i)-Z_int(j))**2.d0 )
      if (distance .lt. accuracy) then
        piece_duplicated = .true.
        exit
      endif
    enddo
    
    ! --- Add piece if not duplicated
    if (.not. piece_duplicated) then
      n_int_tmp = n_int_tmp + 1
      index_int_tmp(n_int_tmp,2) = index_int(i,2)
      index_int_tmp(n_int_tmp,3) = index_int(i,3)
      R_int_tmp(n_int_tmp)       = R_int(i)      
      Z_int_tmp(n_int_tmp)       = Z_int(i)      
    endif
  enddo
  if (n_int .gt. 0) then
    n_int_tmp = n_int_tmp + 1
    index_int_tmp(n_int_tmp,2) = index_int(n_int,2)
    index_int_tmp(n_int_tmp,3) = index_int(n_int,3)
    R_int_tmp(n_int_tmp)       = R_int(n_int)      
    Z_int_tmp(n_int_tmp)       = Z_int(n_int)
  endif  
  
  ! --- Copy non-duplicated pieces
  do i=1,n_int_tmp
    index_int(i,2) = index_int_tmp(i,2)
    index_int(i,3) = index_int_tmp(i,3)
    R_int(i)       = R_int_tmp(i)      
    Z_int(i)       = Z_int_tmp(i)      
    if (debug) write(*,*)'New intersection:',i,R_int(i),Z_int(i)
  enddo
  n_int = n_int_tmp
  
  
  ifail = 0
  return

end subroutine find_wall_crossings_with_flux_surface
  
  




















!> This routine finds all intersections between a flux surface and the wall
subroutine find_wall_crossings_with_flux_surface_old(node_list, element_list, surface, n_int_max, n_int, R_int, Z_int, index_int, ifail)

  use data_structure
  use phys_module
  use high_resolution_wall
  use mod_interp, only: interp_RZ
  use check_point_is_inside_wall_contour
  implicit none
  
  ! --- Routine parameters
  type (type_node_list),        intent(in)      :: node_list
  type (type_element_list),     intent(in)      :: element_list
  type (type_surface),          intent(inout)   :: surface
  integer,                      intent(in)      :: n_int_max
  integer,                      intent(inout)   :: n_int, ifail
  real*8,                       intent(inout)   :: R_int(n_int_max), Z_int(n_int_max)
  integer,                      intent(inout)   :: index_int(n_int_max,3)
  
  ! --- Local variables
  integer                       :: i,j,i_elm,n_tmp,k,l
  integer                       :: node1, node2, node3, node4
  integer                       :: inside, inside_save
  real*8                        :: dl, si, ti, dsi, dti, length, length_min
  real*8                        :: ss1, dss1, ss_int, dl_int
  real*8                        :: ss2, dss2
  real*8                        :: tt1, dtt1, tt_int
  real*8                        :: tt2, dtt2
  real*8                        :: R,R_save,R_mid,R_tmp,R2
  real*8                        :: Z,Z_save,Z_mid,Z_tmp,Z2
  real*8                        :: R_cub1d(4)
  real*8                        :: Z_cub1d(4)
  real*8                        :: distance, distance_max
  real*8                        :: surface_accuracy
  logical                       :: debug, found_point
  
  ! --- Initialise some data
  surface_accuracy = 1.d-2
  n_int = 0
  debug = .false.
  
  ! --- Step along surface with high resolution and check that points are inside wall
  do i=1,surface%n_parts
    do j=surface%parts_index(i),surface%parts_index(i+1)
      i_elm = surface%elm(j)
      
      ! --- Get approximate length of piece
      ss1 = surface%s(1,j)
      tt1 = surface%t(1,j)
      call interp_RZ(node_list,element_list,i_elm,ss1,tt1,R,Z)

      ss2 = surface%s(3,j)
      tt2 = surface%t(3,j)
      call interp_RZ(node_list,element_list,i_elm,ss2,tt2,R2,Z2)
      
      length = sqrt( (R-R2)**2.d0 + (Z-Z2)**2.d0 )
      
      ! --- Deduce number of points on piece required to get accuracy
      n_tmp = length / surface_accuracy
      
      ! --- Get variables for piece extrapolation
      node1 = element_list%element(i_elm)%vertex(1)
      node2 = element_list%element(i_elm)%vertex(2)
      node3 = element_list%element(i_elm)%vertex(3)
      node4 = element_list%element(i_elm)%vertex(4)

      ss1  = surface%s(1,j)
      dss1 = surface%s(2,j)
      ss2  = surface%s(3,j)
      dss2 = surface%s(4,j)

      tt1  = surface%t(1,j)
      dtt1 = surface%t(2,j)
      tt2  = surface%t(3,j)
      dtt2 = surface%t(4,j)

      ! --- Step on each piece point
      do k=1,n_tmp

        dl = -1.d0 + 2.d0 * real(k-1) / real(n_tmp-1)
        call CUB1D(ss1, dss1, ss2, dss2, dl, si, dsi)
        call CUB1D(tt1, dtt1, tt2, dtt2, dl, ti, dti)
        call interp_RZ(node_list,element_list,i_elm,si,ti,R,Z)
        
!write(*,*)'New int:',i,j,k,dl,R,Z
        ! --- Check if point is inside/outside wall
        call check_point_is_inside_wall(R, Z, inside)
        if (k .eq. 1) inside_save = inside

        ! --- Have we stepped accross the wall?
        if (inside .ne. inside_save) then
          ! --- Record the surface piece of the intersection
          n_int = n_int + 1
          index_int(n_int,2) = j
          ! --- Record (find) the wall piece that we've just crossed
          distance_max = 1.d10
          found_point = .false.
          do l=1,n_limiter-1
            ! --- Find the intersection itself
            R_cub1d(1) = R_limiter(l    )               ;       Z_cub1d(1) = Z_limiter(l    )
            R_cub1d(3) = R_limiter(l + 1)               ;       Z_cub1d(3) = Z_limiter(l + 1)
            R_cub1d(2) = (R_cub1d(3)-R_cub1d(1))/2.d0   ;       Z_cub1d(2) = (Z_cub1d(3)-Z_cub1d(1))/2.d0
            R_cub1d(4) = (R_cub1d(3)-R_cub1d(1))/2.d0   ;       Z_cub1d(4) = (Z_cub1d(3)-Z_cub1d(1))/2.d0
            call find_crossing_on_surface_piece(node_list,element_list,surface,j,R_cub1d,Z_cub1d, &
                                                R_tmp,Z_tmp,ss_int,tt_int,dl_int,ifail,.false.)
            if (ifail .eq. 0) then
              R_mid = (R + R_save) / 2.d0
              Z_mid = (Z + Z_save) / 2.d0
              distance = sqrt( (R_tmp - R_mid)**2.d0 + (Z_tmp - Z_mid)**2.d0 )
!write(*,*)'New before:',i,j,k,l,R_tmp,Z_tmp,distance
              if (distance .lt. distance_max) then
                distance_max       = distance
                index_int(n_int,3) = l
                R_int(n_int)       = R_tmp
                Z_int(n_int)       = Z_tmp
                found_point        = .true.
!write(*,*)'New int:   ',i,j,k,l,R_tmp,Z_tmp,distance
!write(*,*)'New int:   ',i,j,k,l,n_int,inside,inside_save,R_tmp,Z_tmp
              endif
            endif
            if ( (ifail .ne. 0) .and. (l .eq. n_limiter-1) .and. (.not. found_point) ) then
              write(*,*)'Warning! Failed to find intersection with wall at piece', index_int(n_int,2)
              write(*,*)'Increasing the size of your initial polar grid may help...'
              if (debug) then
                write(*,*)'DEBUG! Continuing with RZ : ',R,Z
                R_int(n_int) = R
                Z_int(n_int) = Z
              endif
              return
            endif
          enddo
          if (found_point) ifail = 0
        endif

        ! --- Save previous point
        inside_save = inside
        R_save = R
        Z_save = Z

      enddo
        
    enddo
  enddo
  do i=1,n_int
write(*,*)'New int:   ',i,R_int(i),Z_int(i)
  enddo
  
  return

end subroutine find_wall_crossings_with_flux_surface_old
  
  




