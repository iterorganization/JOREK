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
  real*8,                   intent(in)    :: sigmas(17)
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
  
  
