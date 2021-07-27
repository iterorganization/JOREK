subroutine define_flux_values(node_list, element_list, flux_list, sep_list, xcase, psi_xpoint, n_grids, sigmas)
  !-----------------------------------------------------------------------
  ! subroutine defines the flux values of the flux surfaces on which the
  ! finite element grid will be aligned
  !-----------------------------------------------------------------------

  use constants
  use tr_module 
  use data_structure
  use grid_xpoint_data
  use mod_interp
  use phys_module, only:   SDN_threshold
  use equil_info
  
  implicit none
  
  ! --- Routine parameters
  type (type_surface_list), intent(inout) :: flux_list, sep_list
  type (type_node_list),    intent(inout) :: node_list
  type (type_element_list), intent(inout) :: element_list
  integer,                  intent(in)    :: n_grids(12), xcase
  real*8,                   intent(in)    :: sigmas(16)
  real*8,                   intent(inout) :: psi_xpoint(2)
  
  ! --- local variables
  real*8, allocatable :: s_tmp(:)
  integer             :: i, j, nPieces, i_elm
  integer             :: n_flux,      n_open
  integer             :: n_outer,     n_inner
  integer             :: n_private,   n_up_priv 
  integer             :: n_leg,       n_up_leg  
  real*8              :: SIG_closed, SIG_open, SIG_outer, SIG_inner, SIG_private, SIG_up_priv
  real*8              :: SIG_leg_0, SIG_leg_1, SIG_up_leg_0, SIG_up_leg_1
  real*8              :: dPSI_open, dPSI_outer, dPSI_inner, dPSI_private, dPSI_up_priv
  real*8              :: bgf_open, bgf_closed
  real*8              :: RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss
  real*8              :: ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss
  real*8              :: rr, ss, drr, dss, tt
  real*8              :: rr1, ss1, drr1, dss1
  real*8              :: rr2, ss2, drr2, dss2
  real*8              :: psi_bnd, psi_bnd2
  logical             :: xpoint
  
  SIG_closed   = sigmas(1) 
  SIG_open     = sigmas(3) ; SIG_outer    = sigmas(4) ; SIG_inner = sigmas(5)  
  SIG_private  = sigmas(6) ; SIG_up_priv  = sigmas(7) 
  SIG_leg_0    = sigmas(8) ; SIG_leg_1    = sigmas(9) 
  SIG_up_leg_0 = sigmas(10); SIG_up_leg_1 = sigmas(11)
  dPSI_open    = sigmas(12); dPSI_outer   = sigmas(13); dPSI_inner = sigmas(14)
  dPSI_private = sigmas(15); dPSI_up_priv = sigmas(16)

  n_flux    = n_grids(1)
  n_open    = n_grids(3); n_outer   = n_grids(4); n_inner = n_grids(5)
  n_private = n_grids(6); n_up_priv = n_grids(7)
  n_leg     = n_grids(8); n_up_leg  = n_grids(9)
  
  bgf_open   = 0.6d0
  bgf_closed = 0.2d0
  
  write(*,*) '*************************************'
  write(*,*) '* X-point grid : Define Flux Values *'
  write(*,*) '*************************************'
  
  !-------------------------------- Define psi_bnd and psi_bn2
  xpoint        = .true.
  if(xcase .eq. LOWER_XPOINT) psi_bnd = psi_xpoint(1)
  if(xcase .eq. UPPER_XPOINT) psi_bnd = psi_xpoint(2)
  if(xcase .eq. DOUBLE_NULL ) then
    if ( ES%active_xpoint .eq. UPPER_XPOINT ) then
      psi_bnd  = psi_xpoint(2)
      psi_bnd2 = psi_xpoint(1)
    else
      psi_bnd  = psi_xpoint(1)
      psi_bnd2 = psi_xpoint(2)  
    endif
    ! If we have a symmetric double-null, force the single separatrix
    if ( ES%active_xpoint .eq. SYMMETRIC_XPOINT ) then
      psi_xpoint(1) = (psi_xpoint(1)+psi_xpoint(2))/2.d0
      psi_xpoint(2) = psi_xpoint(1)
      psi_bnd  = psi_xpoint(1)
      psi_bnd2 = psi_bnd
    endif
  endif
  
  !-------------------------------- Closed flux surfaces
  call tr_allocate(s_tmp,1,n_flux+1,"s_tmp",CAT_GRID)
  s_tmp = 0
  j     = 0
  call meshac2(n_flux+1,s_tmp,1.d0,9999.d0,SIG_closed,9999.d0,bgf_closed,1.0d0)
  do i=1,n_flux
    flux_list%psi_values(i+j) = ES%psi_axis + (psi_bnd - ES%psi_axis) * s_tmp(i+1)**2
  enddo
  flux_list%psi_values(n_flux) = psi_bnd
  call tr_deallocate(s_tmp,"s_tmp",CAT_GRID)
  
  !-------------------------------- Open flux surfaces (in case of single-null)
  !-------------------------------- OR Sandwich flux surfaces (in case of double-null) - in between the two separatrices
  if (ES%active_xpoint .ne. SYMMETRIC_XPOINT) then ! Ignore in case of symmetric double-null
    call tr_allocate(s_tmp,1,n_open+1,"s_tmp",CAT_GRID)
    s_tmp = 0
    j     = n_flux
    if(xcase .ne. DOUBLE_NULL) then
      call meshac2(n_open+1,s_tmp,0.d0,9999.d0,SIG_open,9999.d0,bgf_open,1.0d0)
      do i=1,n_open
        flux_list%psi_values(i+j) = ES%psi_axis + (psi_bnd - ES%psi_axis) * (1.d0 + dPSI_open*s_tmp(i+1))**2
      enddo
    else
      call meshac2(n_open+1,s_tmp,0.d0,1.d0,SIG_open,SIG_open,0.8d0,1.0d0)
      do i=1,n_open
        flux_list%psi_values(i+j) = psi_bnd + (psi_bnd2 - psi_bnd) * s_tmp(i+1)
      enddo
      flux_list%psi_values(n_flux+n_open) = psi_bnd2
    endif
    call tr_deallocate(s_tmp,"s_tmp",CAT_GRID)
  endif
  
  !-------------------------------- Outer open flux surfaces (in case of double-null)
  if(xcase .eq. DOUBLE_NULL) then
    call tr_allocate(s_tmp,1,n_outer+1,"s_tmp",CAT_GRID)
    s_tmp = 0
    j     = n_flux+n_open
    call meshac2(n_outer+1,s_tmp,0.d0,9999.d0,SIG_outer,9999.d0,0.6d0,1.0d0)
    do i=1,n_outer
      flux_list%psi_values(i+j) = ES%psi_axis + (psi_bnd2 - ES%psi_axis) * (1.d0 + dPSI_outer*s_tmp(i+1))**2
    enddo
    call tr_deallocate(s_tmp,"s_tmp",CAT_GRID)
  endif
  
  !-------------------------------- Inner open flux surfaces (in case of double-null)
  if(xcase .eq. DOUBLE_NULL) then
    call tr_allocate(s_tmp,1,n_inner+1,"s_tmp",CAT_GRID)
    s_tmp = 0
    j     = n_flux+n_open+n_outer
    call meshac2(n_inner+1,s_tmp,0.d0,9999.d0,SIG_inner,9999.d0,0.6d0,1.0d0)
    do i=1,n_inner
      flux_list%psi_values(i+j) = ES%psi_axis + (psi_bnd2 - ES%psi_axis) * (1.d0 + dPSI_inner*s_tmp(i+1))**2
    enddo
    call tr_deallocate(s_tmp,"s_tmp",CAT_GRID)
  endif
  
  !-------------------------------- Lower private flux surfaces
  if(xcase .ne. UPPER_XPOINT) then
    call tr_allocate(s_tmp,1,n_private+1,"s_tmp",CAT_GRID)
    s_tmp = 0
    j     = n_flux+n_open+n_outer+n_inner
    call meshac2(n_private+1,s_tmp,0.d0,9999.d0,SIG_private,9999.d0,0.6d0,1.0d0)
    do i=1,n_private
      flux_list%psi_values(i+j) = ES%psi_axis + (psi_xpoint(1) - ES%psi_axis) * (1.d0 - dPSI_private*s_tmp(i+1))**2
    enddo
    call tr_deallocate(s_tmp,"s_tmp",CAT_GRID)
  endif
  
  !-------------------------------- Upper private flux surfaces
  if(xcase .ne. LOWER_XPOINT) then
    call tr_allocate(s_tmp,1,n_up_priv+1,"s_tmp",CAT_GRID)
    s_tmp = 0
    j     = n_flux+n_open+n_outer+n_inner+n_private
    call meshac2(n_up_priv+1,s_tmp,0.d0,9999.d0,SIG_up_priv,9999.d0,0.6d0,1.0d0)
    do i=1,n_up_priv
      flux_list%psi_values(i+j) = ES%psi_axis + (psi_xpoint(2) - ES%psi_axis) * (1.d0 - dPSI_up_priv*s_tmp(i+1))**2
    enddo
    call tr_deallocate(s_tmp,"s_tmp",CAT_GRID)
  endif
  
  !-------------------------------- Plot separatrices with different colours
  sep_list%psi_values(1) = flux_list%psi_values(n_flux)
  sep_list%psi_values(2) = flux_list%psi_values(n_flux+n_open)
  sep_list%psi_values(3) = flux_list%psi_values(n_flux+n_open+n_private+n_up_priv)
  if(xcase .eq. DOUBLE_NULL) then
    sep_list%psi_values(1) = psi_xpoint(1)
    sep_list%psi_values(2) = psi_xpoint(2)
    sep_list%psi_values(3) = flux_list%psi_values(n_flux+n_open+n_outer)
    sep_list%psi_values(4) = flux_list%psi_values(n_flux+n_open+n_outer+n_inner)
    sep_list%psi_values(5) = flux_list%psi_values(n_flux+n_open+n_outer+n_inner+n_private)
    sep_list%psi_values(6) = flux_list%psi_values(n_flux+n_open+n_outer+n_inner+n_private+n_up_priv)
  endif
  
  call find_flux_surfaces(0,xpoint,xcase,node_list,element_list,flux_list)
  call find_flux_surfaces(0,xpoint,xcase,node_list,element_list,sep_list)  
  if(xcase .eq. DOUBLE_NULL) then
    do i=1,6
      nPieces = sep_list%flux_surfaces(i)%n_pieces
      sep_list%flux_surfaces(i)%n_pieces = 0
      do j=1,nPieces
        i_elm = sep_list%flux_surfaces(i)%elm(j)
        rr1   = sep_list%flux_surfaces(i)%s(1,j)
        drr1  = sep_list%flux_surfaces(i)%s(2,j)
        rr2   = sep_list%flux_surfaces(i)%s(3,j)
        drr2  = sep_list%flux_surfaces(i)%s(4,j)
  
        ss1   = sep_list%flux_surfaces(i)%t(1,j)
        dss1  = sep_list%flux_surfaces(i)%t(2,j)
        ss2   = sep_list%flux_surfaces(i)%t(3,j)
        dss2  = sep_list%flux_surfaces(i)%t(4,j)
        tt    = 0.d0
        call CUB1D(rr1, drr1, rr2, drr2, tt, rr, drr)
        call CUB1D(ss1, dss1, ss2, dss2, tt, ss, dss)
        call interp_RZ(node_list,element_list,i_elm,rr,ss,RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss, &
                                                          ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)
        if(i .eq. 1) then
          if( (ZZg1 .lt. ES%Z_xpoint(2)) .or. (ES%active_xpoint .eq. UPPER_XPOINT) ) then
            sep_list%flux_surfaces(i)%n_pieces = sep_list%flux_surfaces(i)%n_pieces + 1
            sep_list%flux_surfaces(i)%elm(sep_list%flux_surfaces(i)%n_pieces) = sep_list%flux_surfaces(i)%elm(j)
            sep_list%flux_surfaces(i)%s(:,sep_list%flux_surfaces(i)%n_pieces) = sep_list%flux_surfaces(i)%s(:,j)
            sep_list%flux_surfaces(i)%t(:,sep_list%flux_surfaces(i)%n_pieces) = sep_list%flux_surfaces(i)%t(:,j)
          endif
        endif
        if(i .eq. 2) then
          if( (ZZg1 .gt. ES%Z_xpoint(1)) .or. (ES%active_xpoint .eq. LOWER_XPOINT) ) then
            sep_list%flux_surfaces(i)%n_pieces = sep_list%flux_surfaces(i)%n_pieces + 1
            sep_list%flux_surfaces(i)%elm(sep_list%flux_surfaces(i)%n_pieces) = sep_list%flux_surfaces(i)%elm(j)
            sep_list%flux_surfaces(i)%s(:,sep_list%flux_surfaces(i)%n_pieces) = sep_list%flux_surfaces(i)%s(:,j)
            sep_list%flux_surfaces(i)%t(:,sep_list%flux_surfaces(i)%n_pieces) = sep_list%flux_surfaces(i)%t(:,j)
          endif
        endif
        if(i .eq. 3) then
          if( (RRg1 .gt. ES%R_xpoint(1)) .and. (RRg1 .gt. ES%R_xpoint(2)) ) then
            sep_list%flux_surfaces(i)%n_pieces = sep_list%flux_surfaces(i)%n_pieces + 1
            sep_list%flux_surfaces(i)%elm(sep_list%flux_surfaces(i)%n_pieces) = sep_list%flux_surfaces(i)%elm(j)
            sep_list%flux_surfaces(i)%s(:,sep_list%flux_surfaces(i)%n_pieces) = sep_list%flux_surfaces(i)%s(:,j)
            sep_list%flux_surfaces(i)%t(:,sep_list%flux_surfaces(i)%n_pieces) = sep_list%flux_surfaces(i)%t(:,j)
          endif
        endif      
        if(i .eq. 4) then
          if( (RRg1 .lt. ES%R_xpoint(1)) .and. (RRg1 .lt. ES%R_xpoint(2)) ) then
            sep_list%flux_surfaces(i)%n_pieces = sep_list%flux_surfaces(i)%n_pieces + 1
            sep_list%flux_surfaces(i)%elm(sep_list%flux_surfaces(i)%n_pieces) = sep_list%flux_surfaces(i)%elm(j)
            sep_list%flux_surfaces(i)%s(:,sep_list%flux_surfaces(i)%n_pieces) = sep_list%flux_surfaces(i)%s(:,j)
            sep_list%flux_surfaces(i)%t(:,sep_list%flux_surfaces(i)%n_pieces) = sep_list%flux_surfaces(i)%t(:,j)
          endif
        endif      
        if(i .eq. 5) then
          if(ZZg1 .lt. ES%Z_xpoint(1)) then
            sep_list%flux_surfaces(i)%n_pieces = sep_list%flux_surfaces(i)%n_pieces + 1
            sep_list%flux_surfaces(i)%elm(sep_list%flux_surfaces(i)%n_pieces) = sep_list%flux_surfaces(i)%elm(j)
            sep_list%flux_surfaces(i)%s(:,sep_list%flux_surfaces(i)%n_pieces) = sep_list%flux_surfaces(i)%s(:,j)
            sep_list%flux_surfaces(i)%t(:,sep_list%flux_surfaces(i)%n_pieces) = sep_list%flux_surfaces(i)%t(:,j)
          endif
        endif      
        if(i .eq. 6) then
          if(ZZg1 .gt. ES%Z_xpoint(2)) then
            sep_list%flux_surfaces(i)%n_pieces = sep_list%flux_surfaces(i)%n_pieces + 1
            sep_list%flux_surfaces(i)%elm(sep_list%flux_surfaces(i)%n_pieces) = sep_list%flux_surfaces(i)%elm(j)
            sep_list%flux_surfaces(i)%s(:,sep_list%flux_surfaces(i)%n_pieces) = sep_list%flux_surfaces(i)%s(:,j)
            sep_list%flux_surfaces(i)%t(:,sep_list%flux_surfaces(i)%n_pieces) = sep_list%flux_surfaces(i)%t(:,j)
          endif
        endif      
      enddo
    enddo
  else
    nPieces = sep_list%flux_surfaces(3)%n_pieces
    sep_list%flux_surfaces(3)%n_pieces = 0
    do j=1,nPieces
      i_elm = sep_list%flux_surfaces(3)%elm(j)
      rr1   = sep_list%flux_surfaces(3)%s(1,j)
      drr1  = sep_list%flux_surfaces(3)%s(2,j)
      rr2   = sep_list%flux_surfaces(3)%s(3,j)
      drr2  = sep_list%flux_surfaces(3)%s(4,j)
  
      ss1   = sep_list%flux_surfaces(3)%t(1,j)
      dss1  = sep_list%flux_surfaces(3)%t(2,j)
      ss2   = sep_list%flux_surfaces(3)%t(3,j)
      dss2  = sep_list%flux_surfaces(3)%t(4,j)
      tt    = 0.d0
      call CUB1D(rr1, drr1, rr2, drr2, tt, rr, drr)
      call CUB1D(ss1, dss1, ss2, dss2, tt, ss, dss)
      call interp_RZ(node_list,element_list,i_elm,rr,ss,RRg1,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss, &
                                                        ZZg1,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss)
      if( (ZZg1 .lt. ES%Z_xpoint(1)) .and. (xcase .eq. LOWER_XPOINT) ) then
        sep_list%flux_surfaces(3)%n_pieces = sep_list%flux_surfaces(3)%n_pieces + 1
        sep_list%flux_surfaces(3)%elm(sep_list%flux_surfaces(3)%n_pieces) = sep_list%flux_surfaces(3)%elm(j)
        sep_list%flux_surfaces(3)%s(:,sep_list%flux_surfaces(3)%n_pieces) = sep_list%flux_surfaces(3)%s(:,j)
        sep_list%flux_surfaces(3)%t(:,sep_list%flux_surfaces(3)%n_pieces) = sep_list%flux_surfaces(3)%t(:,j)
      endif
      if( (ZZg1 .gt. ES%Z_xpoint(2)) .and. (xcase .eq. UPPER_XPOINT) ) then
        sep_list%flux_surfaces(3)%n_pieces = sep_list%flux_surfaces(3)%n_pieces + 1
        sep_list%flux_surfaces(3)%elm(sep_list%flux_surfaces(3)%n_pieces) = sep_list%flux_surfaces(3)%elm(j)
        sep_list%flux_surfaces(3)%s(:,sep_list%flux_surfaces(3)%n_pieces) = sep_list%flux_surfaces(3)%s(:,j)
        sep_list%flux_surfaces(3)%t(:,sep_list%flux_surfaces(3)%n_pieces) = sep_list%flux_surfaces(3)%t(:,j)
      endif
    enddo
  endif
  
  return
end subroutine define_flux_values
  
  
  
  

  
  
  
  
subroutine redefine_flux_values(node_list, element_list, surface_list, xcase, n_grids, &
                                psi_xpoint, n_int_max, n_target, R_target, Z_target,   &
                                index_target, n_surf_max, n_int_surf, index_int_surf)
  !------------------------------------------------------------------------
  ! subroutine redefines the flux values of the flux surfaces to add values
  ! for the wall corners inside the grid
  !------------------------------------------------------------------------

  use constants
  use tr_module 
  use data_structure
  use phys_module, only : n_limiter, R_limiter, Z_limiter
  use py_plots_grids
  use mod_interp
  use equil_info
  
  implicit none
  
  ! --- Routine parameters
  type (type_surface_list), intent(inout)       :: surface_list
  type (type_node_list),    intent(in)          :: node_list
  type (type_element_list), intent(in)          :: element_list
  real*8,                   intent(in)          :: psi_xpoint(2)
  integer,                  intent(inout)       :: n_grids(12), xcase
  integer,                  intent(inout)       :: n_int_max, n_target, index_target(n_int_max,4)
  real*8,                   intent(inout)       :: R_target(n_int_max), Z_target(n_int_max)
  integer,                  intent(inout)       :: n_surf_max
  integer,                  intent(inout)       :: n_int_surf(n_surf_max) ! number of intersections for each surface
  integer,                  intent(inout)       :: index_int_surf(n_surf_max,n_int_max) ! index of intersections on surface
  
  ! --- local variables
  integer               :: i, j, k, i_surf
  integer               :: i_int, i_int_new
  integer               :: i_edge, index_min, index_min_new
  integer               :: istart, iend, ifail
  integer               :: n_flux,      n_open
  integer               :: n_outer,     n_inner
  integer               :: n_private,   n_up_priv 
  integer               :: n_leg,       n_up_leg  
  integer, parameter    :: LowerLeft =1
  integer, parameter    :: LowerRight=2
  integer, parameter    :: UpperLeft =3
  integer, parameter    :: UpperRight=4
  integer               :: n_pieces
  integer               :: i_beg(4), i_end(4)
  integer               :: i_elm
  real*8                :: ss
  real*8                :: tt
  real*8                :: R,Redge
  real*8                :: Z,Zedge
  real*8                :: Rmin_lower,Rmax_lower
  real*8                :: Rmin_upper,Rmax_upper
  character*256         :: filename
  logical               :: debug, target_only, more_than_one_target_point
  integer               :: n_remove_surface, i_remove_surface(n_surf_max)
  logical               :: remove_surface
  integer               :: n_int_tmp, index_int_tmp(n_int_max,4)
  real*8                :: R_int_tmp(n_int_max), Z_int_tmp(n_int_max)
  integer               :: n_int_surf_new(n_surf_max)
  integer               :: index_int_surf_new(n_surf_max,n_int_max)
  integer               :: n_target_new, index_target_new(n_int_max,4)
  real*8                :: R_target_new(n_int_max), Z_target_new(n_int_max)
  real*8                :: distance, distance_new, distance_max, accuracy
  integer               :: i_int_surf_save(4)
  integer               :: count_target_points(4)
  
  
  write(*,*) '*******************************************'
  write(*,*) '* X-point grid : Add Flux Values For Wall *'
  write(*,*) '*******************************************'
  
  ! --- Debug?
  debug = .true.
  
  n_flux    = n_grids(1)
  n_open    = n_grids(3); n_outer   = n_grids(4); n_inner = n_grids(5)
  n_private = n_grids(6); n_up_priv = n_grids(7)
  n_leg     = n_grids(8); n_up_leg  = n_grids(9)
  
  ! --- Initialise flag of normal surfaces
  do i_surf = 1,surface_list%n_psi
    surface_list%flux_surfaces(i_surf)%flag = 0 
  enddo
  
  
  ! ---------------------------------------------------------------------------------
  ! -------------- First find the begining and the end of each target ---------------
  ! ---------------------------------------------------------------------------------
  
  ! --- First find the beginings of targets (in private region)
  Rmax_lower = -1.d10
  Rmin_lower = 1.d10
  Rmax_upper = -1.d10
  Rmin_upper = 1.d10
  do i_surf = 1,surface_list%n_psi
    do i = 1,n_int_surf(i_surf)
      i_int = index_int_surf(i_surf,i)
      ! --- LowerLeft
      if (index_target(i_int,4) .eq. LowerLeft) then
        if (R_target(i_int) .gt. Rmax_lower) then
          Rmax_lower = R_target(i_int)
          i_beg(LowerLeft) = index_target(i_int,3)
        endif
      endif
      ! --- LowerRight
      if (index_target(i_int,4) .eq. LowerRight) then
        if (R_target(i_int) .lt. Rmin_lower) then
          Rmin_lower = R_target(i_int)
          i_beg(LowerRight) = index_target(i_int,3)
        endif
      endif
      ! --- UpperLeft
      if (index_target(i_int,4) .eq. UpperLeft) then
        if (R_target(i_int) .gt. Rmax_upper) then
          Rmax_upper = R_target(i_int)
          i_beg(UpperLeft) = index_target(i_int,3)
        endif
      endif
      ! --- UpperRight
      if (index_target(i_int,4) .eq. UpperRight) then
        if (R_target(i_int) .lt. Rmin_upper) then
          Rmin_upper = R_target(i_int)
          i_beg(UpperRight) = index_target(i_int,3)
        endif
      endif
    enddo
  enddo
  
  ! --- Then find end of targets in SOL regions
  if (xcase .ne. DOUBLE_NULL) then
    do i_surf = n_flux+1,n_flux+n_open
      do i = 1,n_int_surf(i_surf)
        i_int = index_int_surf(i_surf,i)
        if (xcase .eq. LOWER_XPOINT) then
          if (index_target(i_int,4) .eq. LowerLeft)  i_end(LowerLeft)  = index_target(i_int,3)
          if (index_target(i_int,4) .eq. LowerRight) i_end(LowerRight) = index_target(i_int,3)
        endif
        if (xcase .eq. UPPER_XPOINT) then
          if (index_target(i_int,4) .eq. UpperLeft)  i_end(UpperLeft)  = index_target(i_int,3)
          if (index_target(i_int,4) .eq. UpperRight) i_end(UpperRight) = index_target(i_int,3)
        endif
      enddo
    enddo
  else
    do i_surf = n_flux+n_open+1,n_flux+n_open+n_outer
      do i = 1,n_int_surf(i_surf)
        i_int = index_int_surf(i_surf,i)
        if (index_target(i_int,4) .eq. LowerRight) i_end(LowerRight) = index_target(i_int,3)
        if (index_target(i_int,4) .eq. UpperRight) i_end(UpperRight) = index_target(i_int,3)
      enddo
    enddo
    do i_surf = n_flux+n_open+n_outer+1,n_flux+n_open+n_outer+n_inner
      do i = 1,n_int_surf(i_surf)
        i_int = index_int_surf(i_surf,i)
        if (index_target(i_int,4) .eq. LowerLeft) i_end(LowerLeft) = index_target(i_int,3)
        if (index_target(i_int,4) .eq. UpperLeft) i_end(UpperLeft) = index_target(i_int,3)
      enddo
    enddo
  endif

  
  ! ---------------------------------------------------------------------------------
  ! ----------- Then add Flux surfaces for each wall corner on the targets ----------
  ! ---------------------------------------------------------------------------------
  
  ! --- Now we need to add flux surfaces at target corners
  if (xcase .ne. UPPER_XPOINT) then
    istart = i_beg(LowerLeft)
    iend   = i_end(LowerLeft)
    if (istart .gt. iend) then
      istart = i_end(LowerLeft)
      iend   = i_beg(LowerLeft)
    endif
    do i=istart+1,iend
      call add_flux_surface(node_list, element_list, surface_list, xcase, n_grids, psi_xpoint, &
                            i, n_int_max, n_target, R_target, Z_target, index_target,          &
                            n_surf_max, n_int_surf, index_int_surf)
    enddo
    istart = i_beg(LowerRight)
    iend   = i_end(LowerRight)
    if (istart .gt. iend) then
      istart = i_end(LowerRight)
      iend   = i_beg(LowerRight)
    endif
    do i=istart+1,iend
      call add_flux_surface(node_list, element_list, surface_list, xcase, n_grids, psi_xpoint, &
                            i, n_int_max, n_target, R_target, Z_target, index_target,          &
                            n_surf_max, n_int_surf, index_int_surf)
    enddo
  endif
  if (xcase .ne. LOWER_XPOINT) then
    istart = i_beg(UpperLeft)
    iend   = i_end(UpperLeft)
    if (istart .gt. iend) then
      istart = i_end(UpperLeft)
      iend   = i_beg(UpperLeft)
    endif
    do i=istart+1,iend
      call add_flux_surface(node_list, element_list, surface_list, xcase, n_grids, psi_xpoint, &
                            i, n_int_max, n_target, R_target, Z_target, index_target,          &
                            n_surf_max, n_int_surf, index_int_surf)
    enddo
    istart = i_beg(UpperRight)
    iend   = i_end(UpperRight)
    if (istart .gt. iend) then
      istart = i_end(UpperRight)
      iend   = i_beg(UpperRight)
    endif
    do i=istart+1,iend
      call add_flux_surface(node_list, element_list, surface_list, xcase, n_grids, psi_xpoint, &
                            i, n_int_max, n_target, R_target, Z_target, index_target,          &
                            n_surf_max, n_int_surf, index_int_surf)
    enddo
  endif
  
  ! --- The grid values need to be updated after we added flux surfaces
  n_flux    = n_grids(1)
  n_open    = n_grids(3); n_outer   = n_grids(4); n_inner = n_grids(5)
  n_private = n_grids(6); n_up_priv = n_grids(7)
  n_leg     = n_grids(8); n_up_leg  = n_grids(9)
  
  
  ! ------------------------------------------------------------------------------------------------------------
  ! ----------- Then make sure that we have an artificial target that goes up to the last SOL surface ----------
  ! ------------------------------------------------------------------------------------------------------------
  
  ! --- Re-find end of targets in SOL regions, but this time save it to the intersections
  if (xcase .ne. DOUBLE_NULL) then
    do i_surf = n_flux+1,n_flux+n_open
      do i = 1,n_int_surf(i_surf)
        i_int = index_int_surf(i_surf,i)
        if (xcase .eq. LOWER_XPOINT) then
          if (index_target(i_int,4) .eq. LowerLeft)  i_end(LowerLeft)  = i_int
          if (index_target(i_int,4) .eq. LowerRight) i_end(LowerRight) = i_int
        endif
        if (xcase .eq. UPPER_XPOINT) then
          if (index_target(i_int,4) .eq. UpperLeft)  i_end(UpperLeft)  = i_int
          if (index_target(i_int,4) .eq. UpperRight) i_end(UpperRight) = i_int
        endif
      enddo
    enddo
  else
    do i_surf = n_flux+n_open+1,n_flux+n_open+n_outer
      do i = 1,n_int_surf(i_surf)
        i_int = index_int_surf(i_surf,i)
        if (index_target(i_int,4) .eq. LowerRight) i_end(LowerRight) = i_int
        if (index_target(i_int,4) .eq. UpperRight) i_end(UpperRight) = i_int
      enddo
    enddo
    do i_surf = n_flux+n_open+n_outer+1,n_flux+n_open+n_outer+n_inner
      do i = 1,n_int_surf(i_surf)
        i_int = index_int_surf(i_surf,i)
        if (index_target(i_int,4) .eq. LowerLeft) i_end(LowerLeft) = i_int
        if (index_target(i_int,4) .eq. UpperLeft) i_end(UpperLeft) = i_int
      enddo
    enddo
  endif

  ! --- If the end of the target is not on the last flux surface, we need to add an artificial intersection
  if (xcase .ne. UPPER_XPOINT) then
    ! --- LowerLeft
    i_surf   = n_flux+n_open
    if (xcase .eq. DOUBLE_NULL) i_surf = n_flux+n_open+n_outer+n_inner
    if (index_target(i_end(LowerLeft),1) .ne. i_surf) then
      n_pieces = 1
      i_elm    = surface_list%flux_surfaces(i_surf)%elm(n_pieces)
      ss       = surface_list%flux_surfaces(i_surf)%s(1,n_pieces)
      tt       = surface_list%flux_surfaces(i_surf)%t(1,n_pieces)
      call interp_RZ(node_list,element_list,i_elm,ss,tt,R,Z)
      if (R .gt. ES%R_xpoint(1)) then
        n_pieces = surface_list%flux_surfaces(i_surf)%n_pieces
        i_elm    = surface_list%flux_surfaces(i_surf)%elm(n_pieces)
        ss       = surface_list%flux_surfaces(i_surf)%s(3,n_pieces)
        tt       = surface_list%flux_surfaces(i_surf)%t(3,n_pieces)
        call interp_RZ(node_list,element_list,i_elm,ss,tt,R,Z)
      endif
      n_target = n_target + 1
      R_target(n_target) = R
      Z_target(n_target) = Z
      index_target(n_target,1) = i_surf
      index_target(n_target,2) = n_pieces
      index_target(n_target,3) = -1 ! --- This is dangerous, but this is an artificial intersection, so it has no index on the wall
      index_target(n_target,4) = LowerLeft
    endif
    ! --- LowerRight
    i_surf   = n_flux+n_open
    if (xcase .eq. DOUBLE_NULL) i_surf = n_flux+n_open+n_outer
    if (index_target(i_end(LowerRight),1) .ne. i_surf) then
      n_pieces = 1
      i_elm = surface_list%flux_surfaces(i_surf)%elm(n_pieces)
      ss    = surface_list%flux_surfaces(i_surf)%s(1,n_pieces)
      tt    = surface_list%flux_surfaces(i_surf)%t(1,n_pieces)
      call interp_RZ(node_list,element_list,i_elm,ss,tt,R,Z)
      if (R .lt. ES%R_xpoint(1)) then
        n_pieces = surface_list%flux_surfaces(i_surf)%n_pieces
        i_elm    = surface_list%flux_surfaces(i_surf)%elm(n_pieces)
        ss       = surface_list%flux_surfaces(i_surf)%s(3,n_pieces)
        tt       = surface_list%flux_surfaces(i_surf)%t(3,n_pieces)
        call interp_RZ(node_list,element_list,i_elm,ss,tt,R,Z)
      endif
      n_target = n_target + 1
      R_target(n_target) = R
      Z_target(n_target) = Z
      index_target(n_target,1) = i_surf
      index_target(n_target,2) = n_pieces
      index_target(n_target,3) = -1 ! --- This is dangerous, but this is an artificial intersection, so it has no index on the wall
      index_target(n_target,4) = LowerRight
    endif
  endif
  if (xcase .ne. LOWER_XPOINT) then
    ! --- UpperLeft
    i_surf   = n_flux+n_open
    if (xcase .eq. DOUBLE_NULL) i_surf = n_flux+n_open+n_outer+n_inner
    if (index_target(i_end(UpperLeft),1) .ne. i_surf) then
      n_pieces = 1
      i_elm = surface_list%flux_surfaces(i_surf)%elm(n_pieces)
      ss    = surface_list%flux_surfaces(i_surf)%s(1,n_pieces)
      tt    = surface_list%flux_surfaces(i_surf)%t(1,n_pieces)
      call interp_RZ(node_list,element_list,i_elm,ss,tt,R,Z)
      if (R .gt. ES%R_xpoint(2)) then
        n_pieces = surface_list%flux_surfaces(i_surf)%n_pieces
        i_elm    = surface_list%flux_surfaces(i_surf)%elm(n_pieces)
        ss       = surface_list%flux_surfaces(i_surf)%s(3,n_pieces)
        tt       = surface_list%flux_surfaces(i_surf)%t(3,n_pieces)
        call interp_RZ(node_list,element_list,i_elm,ss,tt,R,Z)
      endif
      n_target = n_target + 1
      R_target(n_target) = R
      Z_target(n_target) = Z
      index_target(n_target,1) = i_surf
      index_target(n_target,2) = n_pieces
      index_target(n_target,3) = -1 ! --- This is dangerous, but this is an artificial intersection, so it has no index on the wall
      index_target(n_target,4) = UpperLeft
    endif
    ! --- UpperRight
    i_surf   = n_flux+n_open
    if (xcase .eq. DOUBLE_NULL) i_surf = n_flux+n_open+n_outer
    if (index_target(i_end(UpperRight),1) .ne. i_surf) then
      n_pieces = 1
      i_elm = surface_list%flux_surfaces(i_surf)%elm(n_pieces)
      ss    = surface_list%flux_surfaces(i_surf)%s(1,n_pieces)
      tt    = surface_list%flux_surfaces(i_surf)%t(1,n_pieces)
      call interp_RZ(node_list,element_list,i_elm,ss,tt,R,Z)
      if (R .lt. ES%R_xpoint(2)) then
        n_pieces = surface_list%flux_surfaces(i_surf)%n_pieces
        i_elm    = surface_list%flux_surfaces(i_surf)%elm(n_pieces)
        ss       = surface_list%flux_surfaces(i_surf)%s(3,n_pieces)
        tt       = surface_list%flux_surfaces(i_surf)%t(3,n_pieces)
        call interp_RZ(node_list,element_list,i_elm,ss,tt,R,Z)
      endif
      n_target = n_target + 1
      R_target(n_target) = R
      Z_target(n_target) = Z
      index_target(n_target,1) = i_surf
      index_target(n_target,2) = n_pieces
      index_target(n_target,3) = -1 ! --- This is dangerous, but this is an artificial intersection, so it has no index on the wall
      index_target(n_target,4) = UpperRight
    endif
  endif
  
  ! --- Print plot file    
  if (debug) then
    filename = 'plot_added_flux_surfaces.py'
    write(*,*)'number of target intersections found with added surfaces:',n_target
    call print_py_plot_prepare_plot(filename)
    call print_py_plot_ordered_flux_surfaces(filename, node_list, element_list, surface_list, 'r', .false.)
    call print_py_plot_wall(filename)
    call print_py_plot_points(filename,n_target,R_target,Z_target)
    call print_py_plot_finish_plot(filename)
  endif
  
  
  ! --------------------------------------------------------------------------------------------------------------------------
  ! ----------- Finally, we need to treat special cases (see grids/grid_utils/wall_grid_documentation.pdf for info) ----------
  ! --------------------------------------------------------------------------------------------------------------------------
  
  ! --- Get target intersections again
  target_only = .true.
  n_target_new = 0
  do i_surf = 1,surface_list%n_psi
    if (i_surf .lt. n_flux) cycle
    call find_wall_crossings_with_flux_surface(node_list, element_list, surface_list%flux_surfaces(i_surf), target_only, &
                                               n_int_max, n_int_tmp, R_int_tmp, Z_int_tmp, index_int_tmp, ifail)

    if (ifail .ne. 0) then
      write(*,*) 'Warning! Failed to find all wall intersections for surface',i_surf,ifail
      return
    endif
    
    ! --- Fill up target arrays
    do i=1,n_int_tmp
      R_target_new(n_target_new + i) = R_int_tmp(i)
      Z_target_new(n_target_new + i) = Z_int_tmp(i)
      index_target_new(n_target_new + i, 1) = i_surf
      index_target_new(n_target_new + i, 2) = index_int_tmp(i,2)
      index_target_new(n_target_new + i, 3) = index_int_tmp(i,3)
      if (Z_target_new(n_target_new + i) .lt. ES%Z_axis) then
        if (R_target_new(n_target_new + i) .lt. ES%R_xpoint(1)) then
          index_target_new(n_target_new + i,4) = LowerLeft
        else
          index_target_new(n_target_new + i,4) = LowerRight
        endif
      else
        if (R_target_new(n_target_new + i) .lt. ES%R_xpoint(2)) then
          index_target_new(n_target_new + i,4) = UpperLeft
        else
          index_target_new(n_target_new + i,4) = UpperRight
        endif
      endif
    enddo
    n_target_new = n_target_new + n_int_tmp
  enddo
  
  ! --- Fill up surface arrays
  do i_surf = 1,surface_list%n_psi
    n_int_surf_new(i_surf) = 0
  enddo
  do i=1,n_target_new
    i_surf = index_target_new(i,1)
    n_int_surf_new(i_surf) = n_int_surf_new(i_surf) + 1
    index_int_surf_new(i_surf,n_int_surf_new(i_surf)) = i
  enddo
  
  ! --- Make sure wall corners are not duplicated
  accuracy = 1.d-4
  n_int_tmp = 0
  do i_surf = 1,surface_list%n_psi
    if (i_surf .lt. n_flux) cycle
    ! --- If this surface has been added for a wall corner then...
    if (surface_list%flux_surfaces(i_surf)%flag .eq. 1) then
      do i=1,n_int_surf_new(i_surf)
        i_int     = index_int_surf(i_surf,1)
        i_int_new = index_int_surf_new(i_surf,i)
        distance = sqrt( (R_target_new(i_int_new) - R_target(i_int))**2.d0 + (Z_target_new(i_int_new) - Z_target(i_int))**2.d0 )
        if (distance .lt. accuracy) cycle
        ! --- If distance is not too small from wall corner, save point
        n_int_tmp = n_int_tmp + 1
        R_int_tmp(n_int_tmp) = R_target_new(i_int_new)
        Z_int_tmp(n_int_tmp) = Z_target_new(i_int_new)
        index_int_tmp(n_int_tmp,:) = index_target_new(i_int_new,:)
      enddo
      ! --- And save the target point itself
      i_int     = index_int_surf(i_surf,1)
      n_int_tmp = n_int_tmp + 1
      R_int_tmp(n_int_tmp) = R_target(i_int)
      Z_int_tmp(n_int_tmp) = Z_target(i_int)
      index_int_tmp(n_int_tmp,:) = index_target(i_int,:)
    else
      do i=1,n_int_surf_new(i_surf)
        i_int_new = index_int_surf_new(i_surf,i)
        n_int_tmp = n_int_tmp + 1
        R_int_tmp(n_int_tmp) = R_target_new(i_int_new)
        Z_int_tmp(n_int_tmp) = Z_target_new(i_int_new)
        index_int_tmp(n_int_tmp,:) = index_target_new(i_int_new,:)
      enddo
    endif
  enddo
  
  ! --- Copy back into main array
  do i=1,n_int_tmp
    R_target_new(i) = R_int_tmp(i)
    Z_target_new(i) = Z_int_tmp(i)
    index_target_new(i,:) = index_int_tmp(i,:)
  enddo
  n_target_new = n_int_tmp
  
  ! --- Fill up surface arrays again
  do i_surf = 1,surface_list%n_psi
    n_int_surf_new(i_surf) = 0
  enddo
  do i=1,n_target_new
    i_surf = index_target_new(i,1)
    n_int_surf_new(i_surf) = n_int_surf_new(i_surf) + 1
    index_int_surf_new(i_surf,n_int_surf_new(i_surf)) = i
  enddo
  
  ! --- Go over each surface and check if it is a special case
  n_remove_surface = 0
  do i_surf = 1,surface_list%n_psi
    if (i_surf .lt. n_flux) cycle
    ! --- Count how many intersections we have on each side LowerLeft/LowerRight/UpperLeft/UpperRight
    count_target_points(1:4) = 0
    do i=1,n_int_surf_new(i_surf)
      i_int_new = index_int_surf_new(i_surf,i)
      count_target_points(index_target_new(i_int_new,4)) = count_target_points(index_target_new(i_int_new,4)) + 1
    enddo
    more_than_one_target_point = .false.
    do i=1,4
      if (count_target_points(i) .gt. 1) more_than_one_target_point = .true.
    enddo
    ! --- Case 1:
    ! --- If surface is not on a wall corner and has more than one intersection with target, redefine target point
    if ( (surface_list%flux_surfaces(i_surf)%flag .eq. 0) .and. (more_than_one_target_point) ) then
      distance_max = -1.d10
      i_int_surf_save(1:4) = 0
      ! --- Loop on all intersections
      do i=1,n_int_surf_new(i_surf)
        i_int_new = index_int_surf_new(i_surf,i)
        ! --- Find the target intersection (ie. the one on the same side LowerLeft/LowerRight/UpperLeft/UpperRight)
        do j=1,n_int_surf(i_surf)
          i_int = index_int_surf(i_surf,j)
          if (index_target(i_int,4) .eq. index_target_new(i_int_new,4)) exit
        enddo
        distance  = sqrt( (R_target_new(i_int_new) - R_target(i_int))**2.d0 + (Z_target_new(i_int_new) - Z_target(i_int))**2.d0 )
        ! --- We choose the point that's the furthest away from the target point
        if (distance .gt. distance_max) then
          distance_max = distance
          i_int_surf_save(index_target_new(i_int_new,4)) = i_int_new
        endif
      enddo
      ! --- Save those points only (on surface index list)
      n_int_surf_new(i_surf) = 0
      do i=1,4
        if (i_int_surf_save(i) .ne. 0) then
          n_int_surf_new(i_surf) = n_int_surf_new(i_surf) + 1
          index_int_surf_new(i_surf,n_int_surf_new(i_surf)) = i_int_surf_save(i)
        endif
      enddo 
    endif
    ! --- Case 2:
    ! --- If surface is on a wall corner and has more than one intersection with target, but target point is last intersection, remove surface
    ! --- Case 3:
    ! --- If surface is on a wall corner and has more than one intersection with target, but target point is not last intersection, redefine target points
    if ( (surface_list%flux_surfaces(i_surf)%flag .eq. 1) .and. (more_than_one_target_point) ) then
      remove_surface = .true.
      ! --- Loop on all intersections
      do i=1,n_int_surf_new(i_surf)
        i_int_new = index_int_surf_new(i_surf,i)
        ! --- Find the target intersection (ie. the one on the same side LowerLeft/LowerRight/UpperLeft/UpperRight)
        do j=1,n_int_surf(i_surf)
          i_int = index_int_surf(i_surf,j)
          if (index_target(i_int,4) .eq. index_target_new(i_int_new,4)) exit
        enddo
        n_pieces = surface_list%flux_surfaces(i_surf)%n_pieces
        index_min     = min(index_target    (i_int    ,2),n_pieces-index_target    (i_int    ,2))
        index_min_new = min(index_target_new(i_int_new,2),n_pieces-index_target_new(i_int_new,2))
        ! --- Check if there is a point that is closer to edge of surface than target point
        if (index_min_new .lt. index_min) then
          remove_surface = .false.
          exit
        endif
        if (index_min_new .eq. index_min) then
          if (index_target(i_int,2) .lt. n_pieces-index_target(i_int,2)) then
            i_edge = 1
          else
            i_edge = 3
          endif
          i_elm = surface_list%flux_surfaces(i_surf)%elm(index_target(i_int,2))
          ss    = surface_list%flux_surfaces(i_surf)%s(i_edge,index_target(i_int,2))
          tt    = surface_list%flux_surfaces(i_surf)%t(i_edge,index_target(i_int,2))
          call interp_RZ(node_list,element_list,i_elm,ss,tt,Redge,Zedge)
          distance     = sqrt( (R_target(i_int)    -Redge)**2.d0 + (Z_target(i_int)    -Zedge)**2.d0 )
          distance_new = sqrt( (R_target(i_int_new)-Redge)**2.d0 + (Z_target(i_int_new)-Zedge)**2.d0 )
          if (distance_new .lt. distance) then
            remove_surface = .false.
            exit
          endif
        endif
      enddo
      ! --- So we need to remove that surface (Case 2)
      if (remove_surface) then
        n_remove_surface = n_remove_surface + 1
        i_remove_surface(n_remove_surface) = i_surf
      else
        ! --- If we keep the surface, this means we need to define a new polar coordinates region (Case 3). Flag surface for later...
        surface_list%flux_surfaces(i_surf)%flag = 2
      endif
    endif
  enddo
  
  ! --- Now remove the surfaces
  do i = 1,n_remove_surface
    i_surf = i_remove_surface(i)
    call remove_flux_surface(surface_list, i_surf, n_grids, &
                             n_int_max, index_target_new,   &
                             n_surf_max, n_int_surf_new, index_int_surf_new)
    ! --- Need to update the remove_surface array
    do j=i+1,n_remove_surface
      i_remove_surface(j) = i_remove_surface(j) - 1
    enddo
  enddo
  
  ! --- Fill up target arrays without the points of surfaces that were removed
  n_int_tmp = 0
  do i_surf = 1,surface_list%n_psi
    do i=1,n_int_surf_new(i_surf)
      i_int_new = index_int_surf_new(i_surf,i)
      n_int_tmp = n_int_tmp + 1
      R_int_tmp(n_int_tmp)       = R_target_new(i_int_new)
      Z_int_tmp(n_int_tmp)       = Z_target_new(i_int_new)
      index_int_tmp(n_int_tmp,:) = index_target_new(i_int_new,:)
      index_int_surf_new(i_surf,i) = n_int_tmp
    enddo
  enddo
  
  ! --- And copy into array that is sent back to main program
  do i=1,n_int_tmp
    R_target(i)       = R_int_tmp(i)      
    Z_target(i)       = Z_int_tmp(i)      
    index_target(i,:) = index_int_tmp(i,:)
  enddo
  n_target = n_int_tmp

  ! --- Consistency check
  do i_surf = 1,surface_list%n_psi
    if ( (n_int_surf_new(i_surf) .gt. 2) .and. (surface_list%flux_surfaces(i_surf)%flag .ne. 2) ) then
      write(*,*)'Warning! Found more than 2 target points on surface',i_surf
      write(*,*)'This surface is not flagged for a new polar coords region.'
      write(*,*)'Something must be wrong, please turn debug flag on and check plots...'
    endif
  enddo
  
  ! --- Print plot file    
  if (debug) then
    filename = 'plot_removed_flux_surfaces.py'
    write(*,*)'number of target intersections found after removing obsolete surfaces:',n_target_new
    call print_py_plot_prepare_plot(filename)
    call print_py_plot_ordered_flux_surfaces(filename, node_list, element_list, surface_list, 'r', .false.)
    call print_py_plot_wall(filename)
    call print_py_plot_points(filename,n_target_new,R_target_new,Z_target_new)
    call print_py_plot_finish_plot(filename)
  endif
  
  
  return
end subroutine redefine_flux_values
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

subroutine add_flux_surface(node_list, element_list, surface_list, xcase, n_grids, psi_xpoint, &
                            i_wall, n_int_max, n_target, R_target, Z_target, index_target,     &
                            n_surf_max, n_int_surf, index_int_surf)
  !-------------------------------------------------------------------------------
  ! subroutine adds a flux surface for a given wall corner with wall index i_wall
  !-------------------------------------------------------------------------------

  use constants  
  use tr_module 
  use data_structure
  use grid_xpoint_data
  use phys_module, only: n_limiter, R_limiter, Z_limiter
  use mod_interp
  use equil_info
  
  implicit none
  
  ! --- Routine parameters
  type (type_surface_list), intent(inout)       :: surface_list
  type (type_node_list),    intent(in)          :: node_list
  type (type_element_list), intent(in)          :: element_list
  real*8,                   intent(in)          :: psi_xpoint(2)
  integer,                  intent(in)          :: i_wall
  integer,                  intent(inout)       :: n_grids(12), xcase
  integer,                  intent(inout)       :: n_int_max, n_target, index_target(n_int_max,4)
  real*8,                   intent(inout)       :: R_target(n_int_max), Z_target(n_int_max)
  integer,                  intent(inout)       :: n_surf_max
  integer,                  intent(inout)       :: n_int_surf(n_surf_max) ! number of intersections for each surface
  integer,                  intent(inout)       :: index_int_surf(n_surf_max,n_int_max) ! index of intersections on surface
  
  ! --- local variables
  type (type_surface_list)      :: surface_list_tmp
  type (type_surface_list)      :: surface_list_single
  real*8                        :: psi_values_tmp(n_surf_max)
  integer                       :: n_int_surf_tmp(n_surf_max) ! number of intersections for each surface
  integer                       :: index_int_surf_tmp(n_surf_max,n_int_max) ! index of intersections on surface
  integer                       :: i, j, i_surf, i_int, i_add, istart, iend, indent, ifail, save_piece
  integer                       :: n_flux,      n_open
  integer                       :: n_outer,     n_inner
  integer                       :: n_private,   n_up_priv 
  integer                       :: n_leg,       n_up_leg  
  integer, parameter            :: LowerLeft =1
  integer, parameter            :: LowerRight=2
  integer, parameter            :: UpperLeft =3
  integer, parameter            :: UpperRight=4
  real*8                        :: R,R_out
  real*8                        :: Z,Z_out
  real*8                        :: s_out, t_out
  integer                       :: i_elm_out
  real*8                        :: psi,dpsi_ds,dpsi_dt,dpsi_dst,dpsi_dss,dpsi_dtt
  integer                       :: location
  logical, parameter            :: xpoint = .true.
  
  
  write(*,*) '*************************************************************'
  write(*,*) '* X-point grid : Adding Flux Value For Wall Point',i_wall
  write(*,*) '*************************************************************'
  
  ! --- Grid values
  n_flux    = n_grids(1)
  n_open    = n_grids(3); n_outer   = n_grids(4); n_inner = n_grids(5)
  n_private = n_grids(6); n_up_priv = n_grids(7)
  n_leg     = n_grids(8); n_up_leg  = n_grids(9)

  ! --- First find the psi-value of the wall corner
  R = R_limiter(i_wall)
  Z = Z_limiter(i_wall)
  call find_RZ(node_list,element_list,R,Z,R_out,Z_out,i_elm_out,s_out,t_out,ifail)
  if (ifail .ne. 0) then
    write(*,*)'Warning! Failed to find RZ position of wall corner ',i_wall,R,Z
    write(*,*)'         on grid when adding flux surface. Aborting...'
    return
  endif
  call interp_RZ(node_list,element_list,i_elm_out,s_out,t_out,R_out,Z_out)
  call interp(node_list,element_list,i_elm_out,1,1,s_out,t_out,psi,dpsi_ds,dpsi_dt,dpsi_dst,dpsi_dss,dpsi_dtt)

  ! --- Determine which region of the grid the wall corner belongs to
  if ( (xcase .ne. UPPER_XPOINT) .and. (abs(psi-ES%psi_axis) .lt. abs(psi_xpoint(1)-ES%psi_axis)) &
                                 .and. (Z_out .lt. ES%Z_axis)                         ) location = private
  if ( (xcase .ne. LOWER_XPOINT) .and. (abs(psi-ES%psi_axis) .lt. abs(psi_xpoint(2)-ES%psi_axis)) &
                                 .and. (Z_out .gt. ES%Z_axis)                         ) location = upper_private
  if ( (xcase .eq. LOWER_XPOINT) .and. (abs(psi-ES%psi_axis) .gt. abs(psi_xpoint(1)-ES%psi_axis)) &
                                 .and. (Z_out .lt. ES%Z_axis)                         ) location = SOL
  if ( (xcase .eq. UPPER_XPOINT) .and. (abs(psi-ES%psi_axis) .gt. abs(psi_xpoint(2)-ES%psi_axis)) &
                                 .and. (Z_out .gt. ES%Z_axis)                         ) location = SOL
  if ( (xcase .eq. DOUBLE_NULL ) .and. (abs(psi-ES%psi_axis) .gt. abs(psi_xpoint(1)-ES%psi_axis)) &
                                 .and. (abs(psi-ES%psi_axis) .lt. abs(psi_xpoint(2)-ES%psi_axis)) ) location = sandwich
  if ( (xcase .eq. DOUBLE_NULL ) .and. (abs(psi-ES%psi_axis) .lt. abs(psi_xpoint(1)-ES%psi_axis)) &
                                 .and. (abs(psi-ES%psi_axis) .gt. abs(psi_xpoint(2)-ES%psi_axis)) ) location = sandwich
  if ( (xcase .eq. DOUBLE_NULL ) .and. (abs(psi-ES%psi_axis) .gt. max(abs(psi_xpoint(1)-ES%psi_axis),abs(psi_xpoint(2)-ES%psi_axis))) &
                                 .and. (Z_out .lt. ES%Z_axis) &
                                 .and. (R_out .gt. ES%R_xpoint(1))                    ) location = outer
  if ( (xcase .eq. DOUBLE_NULL ) .and. (abs(psi-ES%psi_axis) .gt. max(abs(psi_xpoint(1)-ES%psi_axis),abs(psi_xpoint(2)-ES%psi_axis))) &
                                 .and. (Z_out .lt. ES%Z_axis) &
                                 .and. (R_out .lt. ES%R_xpoint(1))                    ) location = inner
  if ( (xcase .eq. DOUBLE_NULL ) .and. (abs(psi-ES%psi_axis) .gt. max(abs(psi_xpoint(1)-ES%psi_axis),abs(psi_xpoint(2)-ES%psi_axis))) &
                                 .and. (Z_out .gt. ES%Z_axis) &
                                 .and. (R_out .gt. ES%R_xpoint(2))                    ) location = outer
  if ( (xcase .eq. DOUBLE_NULL ) .and. (abs(psi-ES%psi_axis) .gt. max(abs(psi_xpoint(1)-ES%psi_axis),abs(psi_xpoint(2)-ES%psi_axis))) &
                                 .and. (Z_out .gt. ES%Z_axis) &
                                 .and. (R_out .lt. ES%R_xpoint(2))                    ) location = inner
  
  ! --- Determine loop parameters depending on location
  if (location .eq. private) then
    istart = n_flux+n_open+n_outer+n_inner+n_private
    iend   = n_flux+n_open+n_outer+n_inner+1
    indent = -1
  endif
  if (location .eq. upper_private) then
    istart = n_flux+n_open+n_outer+n_inner+n_private+n_up_priv
    iend   = n_flux+n_open+n_outer+n_inner+n_private+1
    indent = -1
  endif
  if ((location .eq. SOL) .or. (location .eq. sandwich)) then
    istart = n_flux+1
    iend   = n_flux+n_open
    indent = +1
  endif
  if (location .eq. outer) then
    istart = n_flux+n_open+1
    iend   = n_flux+n_open+n_outer
    indent = +1
  endif
  if (location .eq. inner) then
    istart = n_flux+n_open+n_outer+1
    iend   = n_flux+n_open+n_outer+n_inner
    indent = +1
  endif
  
  ! --- Loop over flux surfaces to determine where we need to add surface
  do i_surf=istart,iend,indent
    if (psi .lt. surface_list%psi_values(i_surf)) then
      i_add = i_surf
      exit
    endif
  enddo
  
  ! --- Find surface (we create an artificial list of surfaces)
  surface_list_single%n_psi = 1
  allocate(surface_list_single%psi_values(1))
  surface_list_single%psi_values(1) = psi
  call find_flux_surfaces(0,xpoint,xcase,node_list,element_list,surface_list_single)
  
  ! --- Record psi-values of temporary flux-surfaces individually
  surface_list_single%flux_surfaces(1)%psi = psi
  
  ! --- Clean up temporary surfaces
  call clean_single_surface(node_list,element_list,surface_list_single%flux_surfaces(1),location)
  
  ! --- Order temporary flux surfaces
  call reorder_flux_surfaces(node_list, element_list, surface_list_single, ifail)
  if (ifail .ne. 0) write(*,*)'Warning! reorder_flux_surfaces failed inside add_flux_surface:',ifail
  
  ! --- Allocate temporary surface list for copy
  surface_list_tmp%n_psi = surface_list%n_psi+1
  allocate(surface_list_tmp%psi_values   (surface_list_tmp%n_psi))
  allocate(surface_list_tmp%flux_surfaces(surface_list_tmp%n_psi))
  
  ! --- Copy all surfaces up to added surface
  do i_surf=1,i_add-1
    ! --- Copy surface
    surface_list_tmp%psi_values(i_surf)                         = surface_list%psi_values(i_surf)
    surface_list_tmp%flux_surfaces(i_surf)%flag                 = surface_list%flux_surfaces(i_surf)%flag
    surface_list_tmp%flux_surfaces(i_surf)%psi                  = surface_list%flux_surfaces(i_surf)%psi
    surface_list_tmp%flux_surfaces(i_surf)%n_pieces             = surface_list%flux_surfaces(i_surf)%n_pieces   
    surface_list_tmp%flux_surfaces(i_surf)%n_parts              = surface_list%flux_surfaces(i_surf)%n_parts    
    surface_list_tmp%flux_surfaces(i_surf)%parts_index(:)       = surface_list%flux_surfaces(i_surf)%parts_index(:)
    do j=1,surface_list%flux_surfaces(i_surf)%n_pieces
      surface_list_tmp%flux_surfaces(i_surf)%elm(j)             = surface_list%flux_surfaces(i_surf)%elm(j)
      surface_list_tmp%flux_surfaces(i_surf)%s(:,j)             = surface_list%flux_surfaces(i_surf)%s(:,j)
      surface_list_tmp%flux_surfaces(i_surf)%t(:,j)             = surface_list%flux_surfaces(i_surf)%t(:,j)
    enddo
    ! --- Copy intersections indexes
    do j=1,n_int_surf(i_surf)
      index_int_surf_tmp(i_surf,j) = index_int_surf(i_surf,j)
    enddo
    n_int_surf_tmp(i_surf) = n_int_surf(i_surf)
  enddo
  
  ! --- Copy all surfaces after added surface
  do i_surf=i_add,surface_list%n_psi
    ! --- Copy surface
    surface_list_tmp%psi_values(i_surf+1)                       = surface_list%psi_values(i_surf)
    surface_list_tmp%flux_surfaces(i_surf+1)%flag               = surface_list%flux_surfaces(i_surf)%flag
    surface_list_tmp%flux_surfaces(i_surf+1)%psi                = surface_list%flux_surfaces(i_surf)%psi
    surface_list_tmp%flux_surfaces(i_surf+1)%n_pieces           = surface_list%flux_surfaces(i_surf)%n_pieces   
    surface_list_tmp%flux_surfaces(i_surf+1)%n_parts            = surface_list%flux_surfaces(i_surf)%n_parts    
    surface_list_tmp%flux_surfaces(i_surf+1)%parts_index(:)     = surface_list%flux_surfaces(i_surf)%parts_index(:)
    do j=1,surface_list%flux_surfaces(i_surf)%n_pieces
      surface_list_tmp%flux_surfaces(i_surf+1)%elm(j)           = surface_list%flux_surfaces(i_surf)%elm(j)
      surface_list_tmp%flux_surfaces(i_surf+1)%s(:,j)           = surface_list%flux_surfaces(i_surf)%s(:,j)
      surface_list_tmp%flux_surfaces(i_surf+1)%t(:,j)           = surface_list%flux_surfaces(i_surf)%t(:,j)
    enddo
    ! --- Copy intersections indexes
    do j=1,n_int_surf(i_surf)
      index_int_surf_tmp(i_surf+1,j) = index_int_surf(i_surf,j)
    enddo
    n_int_surf_tmp(i_surf+1) = n_int_surf(i_surf)
  enddo
  
  ! --- Copy added surface
  surface_list_tmp%psi_values(i_add)                            = surface_list_single%psi_values(1)
  surface_list_tmp%flux_surfaces(i_add)%flag                    = 1
  surface_list_tmp%flux_surfaces(i_add)%psi                     = surface_list_single%flux_surfaces(1)%psi
  surface_list_tmp%flux_surfaces(i_add)%n_pieces                = surface_list_single%flux_surfaces(1)%n_pieces   
  surface_list_tmp%flux_surfaces(i_add)%n_parts                 = surface_list_single%flux_surfaces(1)%n_parts    
  surface_list_tmp%flux_surfaces(i_add)%parts_index(:)          = surface_list_single%flux_surfaces(1)%parts_index(:)
  do j=1,surface_list_single%flux_surfaces(1)%n_pieces
    surface_list_tmp%flux_surfaces(i_add)%elm(j)                = surface_list_single%flux_surfaces(1)%elm(j)
    surface_list_tmp%flux_surfaces(i_add)%s(:,j)                = surface_list_single%flux_surfaces(1)%s(:,j)
    surface_list_tmp%flux_surfaces(i_add)%t(:,j)                = surface_list_single%flux_surfaces(1)%t(:,j)
    ! --- Need to find the surface piece on which the target point lies. Find corresponding i_elm
    if (surface_list_single%flux_surfaces(1)%elm(j) .eq. i_elm_out) save_piece = j
  enddo
  n_target = n_target + 1
  R_target(n_target) = R_out
  Z_target(n_target) = Z_out
  index_target(n_target,1) = i_add
  index_target(n_target,2) = save_piece
  index_target(n_target,3) = i_wall
  if (Z_target(n_target) .lt. ES%Z_axis) then
    if (R_target(n_target) .lt. ES%R_xpoint(1)) then
      index_target(n_target,4) = LowerLeft
    else
      index_target(n_target,4) = LowerRight
    endif
  else
    if (R_target(n_target) .lt. ES%R_xpoint(2)) then
      index_target(n_target,4) = UpperLeft
    else
      index_target(n_target,4) = UpperRight
    endif
  endif
  n_int_surf_tmp(i_add) = 1
  index_int_surf_tmp(i_add,1) = n_target
  
  ! --- Reallocate surface data
  deallocate(surface_list%psi_values)
  deallocate(surface_list%flux_surfaces)
  allocate(surface_list%psi_values   (surface_list_tmp%n_psi))
  allocate(surface_list%flux_surfaces(surface_list_tmp%n_psi))
  
  ! --- Now copy back into surface_list to send back
  do i_surf=1,surface_list_tmp%n_psi
    ! --- Copy surface
    surface_list%psi_values(i_surf)                             = surface_list_tmp%psi_values(i_surf)
    surface_list%flux_surfaces(i_surf)%flag                     = surface_list_tmp%flux_surfaces(i_surf)%flag
    surface_list%flux_surfaces(i_surf)%psi                      = surface_list_tmp%flux_surfaces(i_surf)%psi
    surface_list%flux_surfaces(i_surf)%n_pieces                 = surface_list_tmp%flux_surfaces(i_surf)%n_pieces       
    surface_list%flux_surfaces(i_surf)%n_parts                  = surface_list_tmp%flux_surfaces(i_surf)%n_parts        
    surface_list%flux_surfaces(i_surf)%parts_index(:)           = surface_list_tmp%flux_surfaces(i_surf)%parts_index(:)
    do j=1,surface_list_tmp%flux_surfaces(i_surf)%n_pieces
      surface_list%flux_surfaces(i_surf)%elm(j)                 = surface_list_tmp%flux_surfaces(i_surf)%elm(j)
      surface_list%flux_surfaces(i_surf)%s(:,j)                 = surface_list_tmp%flux_surfaces(i_surf)%s(:,j)
      surface_list%flux_surfaces(i_surf)%t(:,j)                 = surface_list_tmp%flux_surfaces(i_surf)%t(:,j)
    enddo
    ! --- Copy intersections indexes
    do j=1,n_int_surf_tmp(i_surf)
      i_int = index_int_surf_tmp(i_surf,j)
      index_int_surf(i_surf,j) = i_int
      index_target(i_int,1) = i_surf
    enddo
    n_int_surf(i_surf) = n_int_surf_tmp(i_surf)
  enddo
  surface_list%n_psi = surface_list_tmp%n_psi
  
  deallocate(surface_list_tmp%psi_values)
  deallocate(surface_list_tmp%flux_surfaces)
  deallocate(surface_list_single%psi_values)
  deallocate(surface_list_single%flux_surfaces)
  
  ! --- The grid values need to be updated after we added flux surfaces
  if (location .eq. private)            n_private = n_private + 1
  if (location .eq. upper_private)      n_up_priv = n_up_priv + 1
  if (location .eq. SOL)                n_open    = n_open    + 1
  if (location .eq. sandwich)           n_open    = n_open    + 1
  if (location .eq. outer)              n_outer   = n_outer   + 1
  if (location .eq. inner)              n_inner   = n_inner   + 1
  n_grids(1) = n_flux   
  n_grids(3) = n_open   ; n_grids(4) = n_outer  ; n_grids(5) = n_inner
  n_grids(6) = n_private; n_grids(7) = n_up_priv
  n_grids(8) = n_leg    ; n_grids(9) = n_up_leg 
  
  return
end subroutine add_flux_surface
  
















subroutine remove_flux_surface(surface_list, i_remove, n_grids,  &
                               n_int_max, index_target,          &
                               n_surf_max, n_int_surf, index_int_surf)
  !--------------------------------------------------------------------------------------------------------
  ! subroutine removes a flux surface and its intersections with wall for a corresponding index i_remove
  !--------------------------------------------------------------------------------------------------------
  
  use tr_module 
  use data_structure
  use grid_xpoint_data
  use phys_module, only : n_limiter, R_limiter, Z_limiter
  
  implicit none
  
  ! --- Routine parameters
  type (type_surface_list), intent(inout)       :: surface_list
  integer,                  intent(inout)       :: i_remove, n_grids(12)
  integer,                  intent(inout)       :: n_int_max, index_target(n_int_max,4)
  integer,                  intent(inout)       :: n_surf_max
  integer,                  intent(inout)       :: n_int_surf(n_surf_max) ! number of intersections for each surface
  integer,                  intent(inout)       :: index_int_surf(n_surf_max,n_int_max) ! index of intersections on surface
  
  ! --- local variables
  integer                       :: i, j, i_surf, i_int
  integer                       :: n_min, n_max
  integer                       :: n_flux,      n_open
  integer                       :: n_outer,     n_inner
  integer                       :: n_private,   n_up_priv 
  integer                       :: n_leg,       n_up_leg  
  integer                       :: n_target_tmp, index_target_tmp(n_int_max,4)
  real*8                        :: R_target_tmp(n_int_max), Z_target_tmp(n_int_max)
  real*8                        :: psi
  integer                       :: location
  
  
  write(*,*) '*************************************************'
  write(*,*) '* X-point grid : Removing Flux Surface ',i_remove
  write(*,*) '*************************************************'
  
  ! --- Grid values
  n_flux    = n_grids(1)
  n_open    = n_grids(3); n_outer   = n_grids(4); n_inner = n_grids(5)
  n_private = n_grids(6); n_up_priv = n_grids(7)
  n_leg     = n_grids(8); n_up_leg  = n_grids(9)

  ! --- Overwrite surfaces from i_remove up to last surface
  do i_surf=i_remove+1,surface_list%n_psi
    ! --- Copy surface
    surface_list%psi_values(i_surf-1)                           = surface_list%psi_values(i_surf)
    surface_list%flux_surfaces(i_surf-1)%flag                   = surface_list%flux_surfaces(i_surf)%flag
    surface_list%flux_surfaces(i_surf-1)%psi                    = surface_list%flux_surfaces(i_surf)%psi
    surface_list%flux_surfaces(i_surf-1)%n_pieces               = surface_list%flux_surfaces(i_surf)%n_pieces   
    surface_list%flux_surfaces(i_surf-1)%n_parts                = surface_list%flux_surfaces(i_surf)%n_parts    
    surface_list%flux_surfaces(i_surf-1)%parts_index(:)         = surface_list%flux_surfaces(i_surf)%parts_index(:)
    do j=1,surface_list%flux_surfaces(i_surf)%n_pieces
      surface_list%flux_surfaces(i_surf-1)%elm(j)               = surface_list%flux_surfaces(i_surf)%elm(j)
      surface_list%flux_surfaces(i_surf-1)%s(:,j)               = surface_list%flux_surfaces(i_surf)%s(:,j)
      surface_list%flux_surfaces(i_surf-1)%t(:,j)               = surface_list%flux_surfaces(i_surf)%t(:,j)
    enddo
    ! --- Copy intersections indexes
    do j=1,n_int_surf(i_surf)
      i_int = index_int_surf(i_surf,j)
      index_int_surf(i_surf-1,j) = i_int
      index_target(i_int,1) = i_surf-1
    enddo
    n_int_surf(i_surf-1) = n_int_surf(i_surf)
  enddo
  
  ! --- The grid values need to be updated after we added flux surfaces
  
  ! --- Upper private
  n_min = n_flux+n_open+n_outer+n_inner+n_private
  n_max = n_flux+n_open+n_outer+n_inner+n_private+n_up_priv
  if ( (n_min .lt. i_remove) .and. (i_remove .lt. n_max) ) location = upper_private
  
  ! --- Lower private
  n_min = n_flux+n_open+n_outer+n_inner
  n_max = n_flux+n_open+n_outer+n_inner+n_private
  if ( (n_min .lt. i_remove) .and. (i_remove .lt. n_max) ) location = private
  
  ! --- Inner
  n_min = n_flux+n_open+n_outer
  n_max = n_flux+n_open+n_outer+n_inner
  if ( (n_min .lt. i_remove) .and. (i_remove .lt. n_max) ) location = inner
  
  ! --- Outer
  n_min = n_flux+n_open
  n_max = n_flux+n_open+n_outer
  if ( (n_min .lt. i_remove) .and. (i_remove .lt. n_max) ) location = outer
  
  ! --- Open/sandwich
  n_min = n_flux
  n_max = n_flux+n_open
  if ( (n_min .lt. i_remove) .and. (i_remove .lt. n_max) ) location = SOL
  
  ! --- Remove the appropriate one
  if (location .eq. upper_private)      n_up_priv = n_up_priv - 1
  if (location .eq. private)            n_private = n_private - 1
  if (location .eq. inner)              n_inner   = n_inner   - 1
  if (location .eq. outer)              n_outer   = n_outer   - 1
  if (location .eq. SOL)                n_open    = n_open    - 1
  
  ! --- Copy back into array
  n_grids(1) = n_flux   
  n_grids(3) = n_open   ; n_grids(4) = n_outer  ; n_grids(5) = n_inner
  n_grids(6) = n_private; n_grids(7) = n_up_priv
  n_grids(8) = n_leg    ; n_grids(9) = n_up_leg 
  
  return
end subroutine remove_flux_surface
