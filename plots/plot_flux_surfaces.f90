subroutine plot_flux_surfaces(node_list,element_list,surface_list,frame,every_nth,xpoint,xcase)

use constants
use tr_module 
use data_structure
use mod_interp
use phys_module, only: write_ps
use equil_info

implicit none

! --- Routine parameters
type (type_node_list),    intent(in) :: node_list
type (type_element_list), intent(in) :: element_list
type (type_surface_list), intent(in) :: surface_list
logical,                  intent(in) :: frame
integer,                  intent(in) :: every_nth     ! Plot only every_nth flux surface
integer,                  intent(in) :: xcase
logical,                  intent(in) :: xpoint

! --- internal variables
integer            :: i, j, k,ip, nplot, node1, node2, node3, node4, i_elm, found
real*8             :: t, rr1, rr2, drr1, drr2, ss1, ss2, dss1, dss2, ri, si, dri, dsi
real*8             :: R_min, R_max, Z_min, Z_max
real*8             :: psi_bnd, psi_bnd2
real*8,allocatable :: rplot(:), zplot(:)
character*13       :: LABEL


if ( .not. write_ps ) then
  write(*,*) ' Jorek2postscript deactivated. Skipping plot_flux_surfaces'
  return
endif
psi_bnd  = 0.d0
psi_bnd2 = 0.d0

if(xpoint) then
  if (xcase .eq. LOWER_XPOINT) psi_bnd = ES%psi_xpoint(1)
  if (xcase .eq. UPPER_XPOINT) psi_bnd = ES%psi_xpoint(2)
  if (xcase .eq. DOUBLE_NULL ) then
    if (ES%active_xpoint .eq. UPPER_XPOINT) then
      psi_bnd  = ES%psi_xpoint(2)
      psi_bnd2 = ES%psi_xpoint(1)
    else
      psi_bnd  = ES%psi_xpoint(1)
      psi_bnd2 = ES%psi_xpoint(2)
    endif
  endif
endif
found = 0

nplot=11
call tr_allocate(rplot,1,nplot,"rplot",CAT_GRID)
call tr_allocate(zplot,1,nplot,"zplot",CAT_GRID)

LABEL= 'Flux surfaces'

R_min = 1.d20; R_max = -1.d20; Z_min  = 1.d20; Z_max = -1.d20
do i=1,node_list%n_nodes
  R_min = min(R_min,node_list%node(i)%x(1,1,1))
  R_max = max(R_max,node_list%node(i)%x(1,1,1))
  Z_min = min(Z_min,node_list%node(i)%x(1,1,2))
  Z_max = max(Z_max,node_list%node(i)%x(1,1,2))
enddo

if (frame) CALL NFRAME(21,11,1,R_min,R_max,Z_min,Z_max,LABEL,13,'R [m]',5,'Z [m]',5)

!call plot_grid(node_list,element_list,.true.,.false.)                               ! plot the grid

do j = 1, surface_list%n_psi, every_nth

!  write(*,*) ' plot : ',j,surface_list%flux_surfaces(j)%n_pieces
  
  
  ! Make sure that we don't plot private regions twice
  if(xpoint) then
    if ((surface_list%n_psi .gt. 10) .and. (j .gt. 2)) then
      if ((surface_list%psi_values(j) .gt. psi_bnd)  .and. (found .eq. 0)) found = 1
      if ((surface_list%psi_values(j) .gt. psi_bnd2) .and. (xcase .eq. DOUBLE_NULL) .and. (found .eq. 1)) found = 2
      if ((surface_list%psi_values(j) .lt. surface_list%psi_values(j-1))  .and. (found .eq. 3)) found = 4
      if ((surface_list%psi_values(j) .lt. surface_list%psi_values(j-1))  .and. (found .eq. 2)) found = 3
      if ((surface_list%psi_values(j) .gt. surface_list%psi_values(j-1))  .and. (found .eq. 4)) found = 5
      if (xcase .gt. 1) then  ! xcase == 2(UPPER_XPOINT) or 3(DOUBLE_NULL)
        if ((surface_list%psi_values(j) .lt. ES%psi_xpoint(2)) .and. (ES%active_xpoint .eq. UPPER_XPOINT) &
            .and. (     abs(surface_list%psi_values(j)  -surface_list%psi_values(j-1)) .gt. &
             2.d0*abs(surface_list%psi_values(j-1)-surface_list%psi_values(j-2))       ) &
            .and. (found .eq. 4) ) found = 5
      endif

    endif
  endif
  
  do k=1,surface_list%flux_surfaces(j)%n_pieces

    i_elm = surface_list%flux_surfaces(j)%elm(k)

    node1 = element_list%element(i_elm)%vertex(1)
    node2 = element_list%element(i_elm)%vertex(2)
    node3 = element_list%element(i_elm)%vertex(3)
    node4 = element_list%element(i_elm)%vertex(4)

    rr1  = surface_list%flux_surfaces(j)%s(1,k)
    drr1 = surface_list%flux_surfaces(j)%s(2,k)
    rr2  = surface_list%flux_surfaces(j)%s(3,k)
    drr2 = surface_list%flux_surfaces(j)%s(4,k)

    ss1  = surface_list%flux_surfaces(j)%t(1,k)
    dss1 = surface_list%flux_surfaces(j)%t(2,k)
    ss2  = surface_list%flux_surfaces(j)%t(3,k)
    dss2 = surface_list%flux_surfaces(j)%t(4,k)

!    write(*,'(A,i5,4f10.4)') ' element : ',i_elm,rr1,ss1,rr2,ss2

    do ip = 1, nplot

      t = -1. + 2.*float(ip-1)/float(nplot-1)

      call CUB1D(rr1, drr1, rr2, drr2, t, ri, dri)
      call CUB1D(ss1, dss1, ss2, dss2, t, si, dsi)

!      call INTERP2(RR(:,node1),RR(:,node2),RR(:,node3),RR(:,node4),ri,si,rplot(ip),dummy1,dummy2)
!      call INTERP2(ZZ(:,node1),ZZ(:,node2),ZZ(:,node3),ZZ(:,node4),ri,si,zplot(ip),dummy1,dummy2)

      call interp_RZ(node_list,element_list,i_elm,ri,si,rplot(ip),zplot(ip))

    enddo

    call lincol(1)
    if (surface_list%n_psi .le. 6) call lincol(3)

    write(51,*) ' 1.5 setlinewidth '
    
    ! Make sure that we don't plot private regions twice
    if (found .eq. 0) then 
      if ((xpoint) .and. (surface_list%n_psi .gt. 6)) then
        if(     (xcase .eq. LOWER_XPOINT) &
          .and. ((surface_list%psi_values(j) .eq. ES%psi_xpoint(1)) &
          .or.  (minval(zplot) .ge. ES%Z_xpoint(1))) ) call lplot6(21,11,rplot,zplot,-nplot,' ')
        if(     (xcase .eq. UPPER_XPOINT) &
          .and. ((surface_list%psi_values(j) .eq. ES%psi_xpoint(2)) &
          .or.  (maxval(zplot) .le. ES%Z_xpoint(2))) ) call lplot6(21,11,rplot,zplot,-nplot,' ')
        if(     (xcase .eq. DOUBLE_NULL) &
          .and. (((surface_list%psi_values(j) .eq. ES%psi_xpoint(1)) .and. (maxval(zplot) .le. ES%Z_xpoint(2))) &
          .or.   ((surface_list%psi_values(j) .eq. ES%psi_xpoint(2)) .and. (minval(zplot) .ge. ES%Z_xpoint(1))) &
          .or.   ((minval(zplot) .ge. ES%Z_xpoint(1)) &
          .and.   (maxval(zplot) .le. ES%Z_xpoint(2))) ) ) call lplot6(21,11,rplot,zplot,-nplot,' ')
      else
        call lplot6(21,11,rplot,zplot,-nplot,' ')
      endif
    endif
    if (found .eq. 1) then 
      if(     (xcase .eq. LOWER_XPOINT) &
        .and. ( ((surface_list%psi_values(j) .lt. psi_bnd) .and. (maxval(zplot) .lt. ES%Z_xpoint(1))) &
        .or. (surface_list%psi_values(j) .ge. psi_bnd)) ) call lplot6(21,11,rplot,zplot,-nplot,' ')
      if(     (xcase .eq. UPPER_XPOINT) &
        .and. ( ((surface_list%psi_values(j) .lt. psi_bnd) .and. (maxval(zplot) .gt. ES%Z_xpoint(2))) &
        .or. (surface_list%psi_values(j) .ge. psi_bnd)) ) call lplot6(21,11,rplot,zplot,-nplot,' ')
      if (xcase .eq. DOUBLE_NULL) then
        if(     (ES%active_xpoint .eq. UPPER_XPOINT) &
          .and. ((surface_list%psi_values(j) .eq. ES%psi_xpoint(1)) &
          .or.   (minval(zplot) .ge. ES%Z_xpoint(1)))  ) call lplot6(21,11,rplot,zplot,-nplot,' ')
        if(     ( (ES%active_xpoint .eq. LOWER_XPOINT) .or. (ES%active_xpoint .eq. SYMMETRIC_XPOINT) ) &
          .and. ((surface_list%psi_values(j) .eq. ES%psi_xpoint(2)) &
          .or.   (maxval(zplot) .le. ES%Z_xpoint(2)))  ) call lplot6(21,11,rplot,zplot,-nplot,' ')
      endif
    endif
    if (found .eq. 2) then 
      if(     (minval(rplot) .ge. ES%R_xpoint(2)) &
        .and. (maxval(zplot) .ge. 0.d0) ) call lplot6(21,11,rplot,zplot,-nplot,' ')
      if(     (minval(rplot) .ge. ES%R_xpoint(1)) &
        .and. (minval(zplot) .le. 0.d0) ) call lplot6(21,11,rplot,zplot,-nplot,' ')
    endif
    if (found .eq. 3) then 
      if(     (maxval(rplot) .le. ES%R_xpoint(2)) &
        .and. (maxval(zplot) .ge. 0.d0) ) call lplot6(21,11,rplot,zplot,-nplot,' ')
      if(     (maxval(rplot) .le. ES%R_xpoint(1)) &
        .and. (minval(zplot) .le. 0.d0) ) call lplot6(21,11,rplot,zplot,-nplot,' ')
    endif
    if ( (found .eq. 4) .and. (maxval(zplot) .le. ES%Z_xpoint(1)) ) call lplot6(21,11,rplot,zplot,-nplot,' ')
    if ( (found .eq. 5) .and. (minval(zplot) .ge. ES%Z_xpoint(2)) ) call lplot6(21,11,rplot,zplot,-nplot,' ')

  enddo

enddo

call lincol(0)

call tr_deallocate(rplot,"rplot",CAT_GRID)
call tr_deallocate(zplot,"zplot",CAT_GRID)

return
end subroutine plot_flux_surfaces
