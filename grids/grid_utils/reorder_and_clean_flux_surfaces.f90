module reorder_and_clean_flux_surfaces

  real*8, parameter     :: psi_accuracy = 1.d-8
  integer, parameter    :: n_parts_max = 10


contains

!> This routine reorders fluxsurfaces so that pieces are one after the other
subroutine reorder_flux_surfaces(node_list, element_list, surface_list, plot_surfaces, ier)

  use data_structure
  use py_plots_grids
  use grid_xpoint_data
  use mod_interp, only: interp_RZ
  
  implicit none
  
  ! --- Routine parameters
  type (type_node_list),    intent(in)          :: node_list
  type (type_element_list), intent(in)          :: element_list
  type (type_surface_list), intent(inout)       :: surface_list
  logical,                  intent(in)          :: plot_surfaces
  integer,                  intent(inout)       :: ier
  
  ! --- Local variables
  integer       :: i_surf, i, j, nStart
  integer       :: i_piece, i_piece2, i_piece3
  integer       :: found,   found2,   found3
  integer       :: index1, index2, index_save
  integer       :: n_parts, parts_index(n_parts_max)
  integer       :: n_edge_pieces,     index_edge_pieces(2*n_parts_max)
  integer       :: n_isolated_pieces, index_isolated_pieces(2*n_parts_max)
  logical       :: invert, debug, finished, respline
  integer       :: i_elm
  integer       :: i_elm2
  real*8        :: rr,    ss
  real*8        :: rr2,   ss2
  real*8        :: R, dRR_dr, dRR_ds, dRR_drs, dRR_drr, dRR_dss
  real*8        :: R2,dRR2_dr,dRR2_ds,dRR2_drs,dRR2_drr,dRR2_dss
  real*8        :: Z, dZZ_dr, dZZ_ds, dZZ_drs, dZZ_drr, dZZ_dss
  real*8        :: Z2,dZZ2_dr,dZZ2_ds,dZZ2_drs,dZZ2_drr,dZZ2_dss
  real*8        :: R1, Z1
  real*8        :: R3, Z3
  real*8        :: distance
  real*8        :: progress
  character*256 :: filename
  
  write(*,*) '***********************************'
  write(*,*) '*     reorder_flux_surfaces       *'
  write(*,*) '***********************************'
  

  debug    = .true. ! --- Print python files for plots
  respline = .false. ! --- Respline surfaces with stan's method before plotting?
  ier      = 0
  
  ! --- Get a plot?
  if ( debug .and. plot_surfaces)  then
    filename = 'plot_unordered_flux_surfaces.py'
    call py_plot_surface(filename,node_list,element_list,respline, .false., surface_list, -1)
  endif
  
  ! --- Loop over all surfaces
  do i_surf = 1, surface_list%n_psi
    
    progress = 1.d2 * float(i_surf) / float(surface_list%n_psi)
    progress = max(0.d0,progress)
    progress = min(1.d2,progress)
    write(*,"(' Processing  ... ',I0,'%',A,$)") int(progress),CHAR(13)
    
    ! --- Make sure we don't do useless things...
    if (surface_list%flux_surfaces(i_surf)%n_pieces .eq. 0) cycle
    
    ! --- Set zero parts at first
    n_parts = 0
    
    ! --- Find all surface pieces that have no neighbour (ie. end pieces)
    call find_all_edge_pieces(node_list, element_list, surface_list%flux_surfaces(i_surf), &
                              n_edge_pieces, index_edge_pieces, n_isolated_pieces, index_isolated_pieces)
    
    ! --- If there are isolated pieces, put them at the begining
    do i = 1, n_isolated_pieces
      invert = .false.
      index1 = i
      index2 = index_isolated_pieces(i)
      call swap_surface_pieces(surface_list%flux_surfaces(i_surf), index1, index2, invert, &
                               n_edge_pieces, index_edge_pieces, n_isolated_pieces, index_isolated_pieces)
      n_parts        = i
      parts_index(i) = i
    enddo
    
    ! ---------------------------------------------------------
    ! --- In case there are no edge pieces (ie. closed surface)
    ! ---------------------------------------------------------
    if (n_edge_pieces .eq. 0) then
      
      ! --- Start at first piece (check if there are isolated pieces)
      nStart = 1
      if (n_isolated_pieces .gt. 0) nStart = parts_index(n_parts) + 1
      n_parts = n_parts + 1
      parts_index(n_parts) = nStart
      
      ! --- There might be several closed surfaces
      finished = .false.
      do while(.not. finished)
        ! --- Loop over all remaining pieces
        do i_piece = nStart, surface_list%flux_surfaces(i_surf)%n_pieces - 1
          if (i_piece .eq. surface_list%flux_surfaces(i_surf)%n_pieces - 1) finished = .true.
          call get_next_surface_piece(node_list, element_list, surface_list%flux_surfaces(i_surf), i_piece, &
                                      n_edge_pieces, index_edge_pieces, n_isolated_pieces, index_isolated_pieces, found)
          ! --- If we didn't find the next piece, this means there could be another closed surface (unlikely but try once at least)
          if (found .eq. 0) then
            rr    = surface_list%flux_surfaces(i_surf)%s(3,i_piece)
            ss    = surface_list%flux_surfaces(i_surf)%t(3,i_piece)
            i_elm = surface_list%flux_surfaces(i_surf)%elm(i_piece)
            call interp_RZ(node_list,element_list,i_elm,rr,ss,R,Z)
            
            rr2    = surface_list%flux_surfaces(i_surf)%s(1,parts_index(n_parts))
            ss2    = surface_list%flux_surfaces(i_surf)%t(1,parts_index(n_parts))
            i_elm2 = surface_list%flux_surfaces(i_surf)%elm(parts_index(n_parts))
            call interp_RZ(node_list,element_list,i_elm2,rr2,ss2,R2,Z2)
            
            distance = sqrt( (R-R2)**2.d0 + (Z-Z2)**2.d0 )
            if (distance .lt. psi_accuracy) then
              nStart = i_piece + 1
              n_parts = n_parts + 1
              parts_index(n_parts) = nStart
              exit
            else
              write(*,'(A,1i4,1e16.8)') 'Warning! Failed to find parts of the surface',i_surf,distance
              ier = 1
            endif
          endif
          
        enddo
      enddo
    
    ! -----------------------------------------------------
    ! --- In case there are edge pieces (ie. open surfaces)
    ! -----------------------------------------------------
    else
      if (mod(n_edge_pieces,2) .ne. 0) then
        write(*,'(A,1i4)') 'Warning! There are an odd number of edge pieces for surface',i_surf
        ier = 3
        return
      endif
      
      ! --- Start at first piece (check if there are isolated pieces before)
      nStart = 1
      if (n_isolated_pieces .gt. 0) nStart = parts_index(n_parts) + 1
      n_parts = n_parts + 1
      parts_index(n_parts) = nStart
      
      ! --- Loop over end pieces (two per surface)
      finished = .true.
      do i=1,n_edge_pieces/2
        
        ! --- First swap piece so that it's at the begining of our new part.
        invert = .false.
        index1 = index_edge_pieces(i)
        index2 = nStart
        call swap_surface_pieces(surface_list%flux_surfaces(i_surf), index1, index2, invert, &
                                 n_edge_pieces, index_edge_pieces, n_isolated_pieces, index_isolated_pieces)
        
        ! --- Loop over all remaining pieces
        do i_piece = nStart, surface_list%flux_surfaces(i_surf)%n_pieces - 1
          
          call get_next_surface_piece(node_list, element_list, surface_list%flux_surfaces(i_surf), i_piece, &
                                      n_edge_pieces, index_edge_pieces, n_isolated_pieces, index_isolated_pieces, found)
          
          ! --- If we didn't find the next piece, check if this should be the end
          if (found .eq. 0) then
            
            ! --- Make sure that if this piece is an edge piece, its index is > n_edge_pieces/2
            do j=1,n_edge_pieces
              if (index_edge_pieces(j) .eq. i_piece) then
                index_save = index_edge_pieces(j)
                index_edge_pieces(j) = index_edge_pieces(n_edge_pieces/2+i)
                index_edge_pieces(n_edge_pieces/2+i) = index_save
                exit
              endif
            enddo
            
            ! --- Let's get to the next surface part
            if (i .ne. n_edge_pieces/2) then
              n_parts = n_parts + 1
              parts_index(n_parts) = i_piece + 1
              exit
            else
            ! --- We ran out of edge pieces but there is more, meaning that there is probably a closed surface as well
              if (i_piece .ne. surface_list%flux_surfaces(i_surf)%n_pieces-1) then
                finished = .false.
                nStart   = i_piece + 1
                n_parts  = n_parts + 1
                parts_index(n_parts) = i_piece + 1
              endif
            endif
    
          endif
            
        enddo
        nStart = parts_index(n_parts)
      
        if (.not. finished) then
          ! --- There might be several closed surfaces
          do while(.not. finished)
            ! --- Loop over all remaining pieces
            do i_piece = nStart, surface_list%flux_surfaces(i_surf)%n_pieces - 1
              if (i_piece .eq. surface_list%flux_surfaces(i_surf)%n_pieces - 1) finished = .true.
              call get_next_surface_piece(node_list, element_list, surface_list%flux_surfaces(i_surf), i_piece, &
                                          n_edge_pieces, index_edge_pieces, n_isolated_pieces, index_isolated_pieces, found)
              ! --- If we didn't find the next piece, this means there could be another closed surface (unlikely but try once at least)
              if (found .eq. 0) then
                rr    = surface_list%flux_surfaces(i_surf)%s(3,i_piece)
                ss    = surface_list%flux_surfaces(i_surf)%t(3,i_piece)
                i_elm = surface_list%flux_surfaces(i_surf)%elm(i_piece)
                call interp_RZ(node_list,element_list,i_elm,rr,ss,R,Z)
                
                rr2    = surface_list%flux_surfaces(i_surf)%s(1,parts_index(n_parts))
                ss2    = surface_list%flux_surfaces(i_surf)%t(1,parts_index(n_parts))
                i_elm2 = surface_list%flux_surfaces(i_surf)%elm(parts_index(n_parts))
                call interp_RZ(node_list,element_list,i_elm2,rr2,ss2,R2,Z2)
                
                distance = sqrt( (R-R2)**2.d0 + (Z-Z2)**2.d0 )
                if (distance .lt. psi_accuracy) then
                  nStart = i_piece + 1
                  n_parts = n_parts + 1
                  parts_index(n_parts) = nStart
                  exit
                else
                  write(*,'(A,1i4,1e16.8)') 'Warning! Failed to find parts of the surface',i_surf,distance
                  ier = 1
                endif
              endif
              
            enddo
          enddo
          
        endif
      
      enddo
    
    endif
    
    ! --- Save parts indexes
    surface_list%flux_surfaces(i_surf)%n_parts = n_parts
    do i=1,n_parts
      surface_list%flux_surfaces(i_surf)%parts_index(i) = parts_index(i)
    enddo
    
    ! --- Print the number of pieces of the last part
    surface_list%flux_surfaces(i_surf)%parts_index(n_parts+1) = surface_list%flux_surfaces(i_surf)%n_pieces + 1
    
    ! --- Loop over ends of surfaces to make sure they are in the right order (first the beginnings)
    do i=1,surface_list%flux_surfaces(i_surf)%n_parts
      i_piece = surface_list%flux_surfaces(i_surf)%parts_index(i)
      rr    = surface_list%flux_surfaces(i_surf)%s(1,i_piece)
      ss    = surface_list%flux_surfaces(i_surf)%t(1,i_piece)
      i_elm = surface_list%flux_surfaces(i_surf)%elm(i_piece)
      call interp_RZ(node_list,element_list,i_elm,rr,ss,R1,dRR_dr,dRR_ds,dRR_drs,dRR_drr,dRR_dss, &
                                                        Z1,dZZ_dr,dZZ_ds,dZZ_drs,dZZ_drr,dZZ_dss)
      rr2    = surface_list%flux_surfaces(i_surf)%s(3,i_piece)
      ss2    = surface_list%flux_surfaces(i_surf)%t(3,i_piece)
      i_elm2 = surface_list%flux_surfaces(i_surf)%elm(i_piece)
      call interp_RZ(node_list,element_list,i_elm2,rr2,ss2,R2,dRR2_dr,dRR2_ds,dRR2_drs,dRR2_drr,dRR2_dss, &
                                                           Z2,dZZ2_dr,dZZ2_ds,dZZ2_drs,dZZ2_drr,dZZ2_dss)
      
      i_piece = surface_list%flux_surfaces(i_surf)%parts_index(i) + 1
      rr2    = surface_list%flux_surfaces(i_surf)%s(1,i_piece)
      ss2    = surface_list%flux_surfaces(i_surf)%t(1,i_piece)
      i_elm2 = surface_list%flux_surfaces(i_surf)%elm(i_piece)
      call interp_RZ(node_list,element_list,i_elm2,rr2,ss2,R3,dRR2_dr,dRR2_ds,dRR2_drs,dRR2_drr,dRR2_dss, &
                                                           Z3,dZZ2_dr,dZZ2_ds,dZZ2_drs,dZZ2_drr,dZZ2_dss)
      
      if ( sqrt( (R1-R3)**2 + (Z1-Z3)**2 ) .lt. sqrt( (R2-R3)**2 + (Z2-Z3)**2 ) ) then
        index1 = surface_list%flux_surfaces(i_surf)%parts_index(i)
        invert = .true.
        call swap_surface_pieces(surface_list%flux_surfaces(i_surf), index1, index1, invert, n_edge_pieces, index_edge_pieces, n_isolated_pieces, index_isolated_pieces)
      endif
    enddo
  
    ! --- Loop over ends of surfaces to make sure they are in the right order (then the ends)
    do i=1,surface_list%flux_surfaces(i_surf)%n_parts
      i_piece = surface_list%flux_surfaces(i_surf)%parts_index(i+1) - 1
      rr    = surface_list%flux_surfaces(i_surf)%s(3,i_piece)
      ss    = surface_list%flux_surfaces(i_surf)%t(3,i_piece)
      i_elm = surface_list%flux_surfaces(i_surf)%elm(i_piece)
      call interp_RZ(node_list,element_list,i_elm,rr,ss,R1,dRR_dr,dRR_ds,dRR_drs,dRR_drr,dRR_dss, &
                                                        Z1,dZZ_dr,dZZ_ds,dZZ_drs,dZZ_drr,dZZ_dss)
      rr2    = surface_list%flux_surfaces(i_surf)%s(1,i_piece)
      ss2    = surface_list%flux_surfaces(i_surf)%t(1,i_piece)
      i_elm2 = surface_list%flux_surfaces(i_surf)%elm(i_piece)
      call interp_RZ(node_list,element_list,i_elm2,rr2,ss2,R2,dRR2_dr,dRR2_ds,dRR2_drs,dRR2_drr,dRR2_dss, &
                                                           Z2,dZZ2_dr,dZZ2_ds,dZZ2_drs,dZZ2_drr,dZZ2_dss)
      
      i_piece = surface_list%flux_surfaces(i_surf)%parts_index(i+1) - 2
      i_piece = max(1,i_piece)
      rr2    = surface_list%flux_surfaces(i_surf)%s(3,i_piece)
      ss2    = surface_list%flux_surfaces(i_surf)%t(3,i_piece)
      i_elm2 = surface_list%flux_surfaces(i_surf)%elm(i_piece)
      call interp_RZ(node_list,element_list,i_elm2,rr2,ss2,R3,dRR2_dr,dRR2_ds,dRR2_drs,dRR2_drr,dRR2_dss, &
                                                           Z3,dZZ2_dr,dZZ2_ds,dZZ2_drs,dZZ2_drr,dZZ2_dss)
      
      if ( sqrt( (R1-R3)**2 + (Z1-Z3)**2 ) .lt. sqrt( (R2-R3)**2 + (Z2-Z3)**2 ) ) then
        index1 = surface_list%flux_surfaces(i_surf)%parts_index(i+1) - 1
        invert = .true.
        call swap_surface_pieces(surface_list%flux_surfaces(i_surf), index1, index1, invert, n_edge_pieces, index_edge_pieces, n_isolated_pieces, index_isolated_pieces)
      endif
    enddo
  
  enddo
  
  write(*,*) 'Processing  ... 100'
  write(*,*) 'finished reordering'
  if (debug .and. plot_surfaces) then
    filename = 'plot_ordered_flux_surfaces.py'
    call py_plot_surface(filename,node_list,element_list,.false., .true., surface_list, -1)
  endif
  
  return

end subroutine reorder_flux_surfaces
  
  

! ----------------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------------

!> This routine finds the next piece given a piece index of a fluxsurface
subroutine get_next_surface_piece(node_list, element_list, surface, i_piece, &
                                  n_edge_pieces, index_edge_pieces, n_isolated_pieces, index_isolated_pieces, found)
  
  use data_structure
  use mod_interp, only: interp_RZ
  implicit none
  
  ! --- Routine parameters
  type (type_node_list),    intent(in)          :: node_list
  type (type_element_list), intent(in)          :: element_list
  type (type_surface),      intent(inout)       :: surface
  integer,                  intent(in)          :: i_piece
  integer,                  intent(inout)       :: found
  integer,                  intent(inout)       :: n_edge_pieces,     index_edge_pieces(2*n_parts_max)
  integer,                  intent(inout)       :: n_isolated_pieces, index_isolated_pieces(2*n_parts_max)
  
  ! --- Internal parameters
  integer       :: i, j, k
  integer       :: index1, index2
  logical       :: invert
  integer       :: i_elm
  integer       :: i_elm2
  real*8        :: rr,    ss
  real*8        :: rr2,   ss2
  real*8        :: R, R2, Z, Z2
  real*8        :: distance
  
  found = 0
  
  ! --- Get last point of that surface piece
  rr    = surface%s(3,i_piece)
  ss    = surface%t(3,i_piece)
  i_elm = surface%elm(i_piece)
  call interp_RZ(node_list,element_list,i_elm,rr,ss,R,Z)

  ! --- Loop over all remaining pieces
  do i=i_piece+1,surface%n_pieces
  
    ! --- Try both ends of the piece
    do j=1,3,2
      rr2    = surface%s(j,i)
      ss2    = surface%t(j,i)
      i_elm2 = surface%elm(i)

      call interp_RZ(node_list,element_list,i_elm2,rr2,ss2,R2,Z2)
      
      distance = sqrt( (R-R2)**2.d0 + (Z-Z2)**2.d0 )
      if (distance .lt. psi_accuracy) then
        found = j
        exit
      endif
    enddo
  
    ! --- Have we found the next piece?
    if (found .ne. 0) then
      index1 = i_piece + 1
      index2 = i
      invert = .false.
      if (found .eq. 3) invert = .true.
      call swap_surface_pieces(surface, index1, index2, invert, n_edge_pieces, index_edge_pieces, n_isolated_pieces, index_isolated_pieces)
      ! --- We should not invert any edge pieces! If we did, reswap it with itself to invert again
      if (invert) then
        if(n_edge_pieces .gt. 0) then
          do k=1,n_edge_pieces
            if (index_edge_pieces(k) .eq. i) then
              call swap_surface_pieces(surface, index2, index2, invert, n_edge_pieces, index_edge_pieces, n_isolated_pieces, index_isolated_pieces)
              exit
            endif
          enddo
        endif
      endif
  
      exit
    endif

  enddo
  
  return

end subroutine get_next_surface_piece





! ----------------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------------


!> This routine swaps two pieces of a given fluxsurface
subroutine swap_surface_pieces(surface, index1, index2, invert, n_edge_pieces, index_edge_pieces, n_isolated_pieces, index_isolated_pieces)

  use data_structure
  implicit none
  
  ! --- Routine parameters
  type (type_surface),      intent(inout)       :: surface
  integer,                  intent(in)          :: index1, index2
  logical,                  intent(in)          :: invert
  integer,                  intent(inout)       :: n_edge_pieces,     index_edge_pieces(2*n_parts_max)
  integer,                  intent(inout)       :: n_isolated_pieces, index_isolated_pieces(2*n_parts_max)
  
  ! --- Internal parameters
  integer       :: i
  integer       :: elm_save
  real*8        :: s_save(4), t_save(4)
  
  ! --- Save firsy piece
  elm_save  = surface%elm(index1)
  s_save(:) = surface%s(:,index1)
  t_save(:) = surface%t(:,index1)
  
  ! --- And swap
  if (invert) then
    ! --- First piece
    surface%elm(index1) = surface%elm(index2)
    surface%s(1,index1) = surface%s(3,index2)
    surface%t(1,index1) = surface%t(3,index2)
    surface%s(2,index1) = -surface%s(4,index2)
    surface%t(2,index1) = -surface%t(4,index2)
    surface%s(3,index1) = surface%s(1,index2)
    surface%t(3,index1) = surface%t(1,index2)
    surface%s(4,index1) = -surface%s(2,index2)
    surface%t(4,index1) = -surface%t(2,index2)
    
    ! --- Second piece
    surface%elm(index2) = elm_save
    surface%s(1,index2) = s_save(3)
    surface%t(1,index2) = t_save(3)
    surface%s(2,index2) = -s_save(4)
    surface%t(2,index2) = -t_save(4)
    surface%s(3,index2) = s_save(1)
    surface%t(3,index2) = t_save(1)
    surface%s(4,index2) = -s_save(2)
    surface%t(4,index2) = -t_save(2)
    
  else
    ! --- First piece
    surface%elm(index1) = surface%elm(index2)
    surface%s(:,index1) = surface%s(:,index2)
    surface%t(:,index1) = surface%t(:,index2)
    
    ! --- Second piece
    surface%elm(index2) = elm_save
    surface%s(:,index2) = s_save(:)
    surface%t(:,index2) = t_save(:)
  endif
  
  ! --- Then make sure the edge_pieces indices are swapped too
  if(n_edge_pieces .gt. 0) then
    do i=1,n_edge_pieces
      if (index_edge_pieces(i) .eq. index1) then
        index_edge_pieces(i) = index2
        cycle
      endif
      if (index_edge_pieces(i) .eq. index2) then
        index_edge_pieces(i) = index1
        cycle
      endif
    enddo
  endif
  
  ! --- And make sure the isolated_pieces indices are swapped too
  if(n_isolated_pieces .gt. 0) then
    do i=1,n_isolated_pieces
      if (index_isolated_pieces(i) .eq. index1) then
        index_isolated_pieces(i) = index2
        cycle
      endif
      if (index_isolated_pieces(i) .eq. index2) then
        index_isolated_pieces(i) = index1
        cycle
      endif
    enddo
  endif
  
  return
end subroutine swap_surface_pieces
  








! ----------------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------------


!> This routine finds all the pieces of a given fluxsurface that have no neighbour
subroutine find_all_edge_pieces(node_list, element_list, surface, n_edge_pieces, index_edge_pieces, n_isolated_pieces, index_isolated_pieces)



  use data_structure
  use mod_interp, only: interp_RZ
  implicit none
  
  ! --- Routine parameters
  type (type_node_list),    intent(in)          :: node_list
  type (type_element_list), intent(in)          :: element_list
  type (type_surface),      intent(inout)       :: surface
  integer,                  intent(inout)       :: n_edge_pieces,     index_edge_pieces(2*n_parts_max)
  integer,                  intent(inout)       :: n_isolated_pieces, index_isolated_pieces(2*n_parts_max)
  
  ! --- Internal parameters
  integer       :: i1, i2
  integer       :: k1, k2
  logical       :: found(2), invert
  integer       :: i_elm
  integer       :: i_elm2
  real*8        :: rr,    ss
  real*8        :: rr2,   ss2
  real*8        :: R,R2,Z,Z2
  real*8        :: distance
  
  n_isolated_pieces = 0
  n_edge_pieces     = 0
  
  ! --- Check each piece
  do i1=1,surface%n_pieces
    found(1) = .false.
    found(2) = .false.
    ! --- Check both sides of the piece
    do k1=1,3,2
      ! --- Get edge point of that surface piece
      rr    = surface%s(k1,i1)
      ss    = surface%t(k1,i1)
      i_elm = surface%elm(i1)
      call interp_RZ(node_list,element_list,i_elm,rr,ss,R,Z)

      ! --- Loop over all other pieces
      do i2=1,surface%n_pieces
        if (i2 .ne. i1) then
          
          ! --- Check both sides of the piece
          do k2=1,3,2
            rr2    = surface%s(k2,i2)
            ss2    = surface%t(k2,i2)
            i_elm2 = surface%elm(i2)

            call interp_RZ(node_list,element_list,i_elm2,rr2,ss2,R2,Z2)
            
            distance = sqrt( (R-R2)**2.d0 + (Z-Z2)**2.d0 )
            if (distance .lt. psi_accuracy) then
              if (k1 .eq. 1) found(1) = .true.
              if (k1 .eq. 3) found(2) = .true.
              exit
            endif
          enddo
          
          if ( (k1 .eq. 1) .and. (found(1)) ) exit
          if ( (k1 .eq. 3) .and. (found(2)) ) exit
          
        endif
      enddo
      
    enddo
    
    ! --- Have we found an isolated pieces?
    if ( (.not. found(1)) .and. (.not. found(2)) ) then
      n_isolated_pieces = n_isolated_pieces + 1
      index_isolated_pieces(n_isolated_pieces) = i1
    else
      ! --- Have we found an edge pieces?
      if ( (.not. found(1)) .or. (.not. found(2)) ) then
        n_edge_pieces = n_edge_pieces + 1
        index_edge_pieces(n_edge_pieces) = i1
        
        ! --- The edge pieces need to start at the edge. Invert the piece with itself if needed
        if (.not. found(2)) then
          invert = .true.
          call swap_surface_pieces(surface, i1, i1, invert, n_edge_pieces, index_edge_pieces, n_isolated_pieces, index_isolated_pieces)
        endif
      endif
    endif
      
  enddo
  
end subroutine find_all_edge_pieces










! ----------------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------------


!> This routine removes the private region surface pieces under a given psi_value
subroutine clean_surfaces(node_list,element_list,flux_list,n_grids)


  use constants
  use data_structure
  use phys_module, only : xcase
  use grid_xpoint_data
  use py_plots_grids
  use mod_interp, only: interp_RZ
  use equil_info
  use check_point_is_inside_wall_contour
  
  implicit none
  
  ! --- Routine parameters
  type (type_node_list),        intent(in)      :: node_list
  type (type_element_list),     intent(in)      :: element_list
  type (type_surface_list),     intent(inout)   :: flux_list
  integer,                      intent(in)      :: n_grids(12) 
  
  ! --- Internal parameters
  type (type_surface_list) :: sep_list
  type (type_surface)   :: surface
  integer               :: location, ifail
  integer               :: n_flux,   n_open,   n_outer,   n_inner,   n_private,   n_up_priv  
  integer               :: i_surf, i_part, i_piece, i_part_save, i_pieces_max
  integer               :: i_elm, inside, inside2, count
  real*8                :: rr,    ss
  real*8                :: R,dR_dr, dR_ds, dR_drs, dR_drr, dR_dss
  real*8                :: Z,dZ_dr, dZ_ds, dZ_drs, dZ_drr, dZ_dss
  character*256         :: filename
  character*1           :: colour
  logical               :: dashed
  integer               :: debug
  
  debug = 1
  
  n_flux    = n_grids(1)
  n_open    = n_grids(3); n_outer   = n_grids(4); n_inner = n_grids(5)
  n_private = n_grids(6); n_up_priv = n_grids(7)

  if (debug .ge. 1) then
    filename = 'plot_cleaned_flux_surfaces.py'
    call print_py_plot_prepare_plot(filename)
    colour = 'r' ; dashed = .false.
    call print_py_plot_ordered_flux_surfaces(filename, node_list, element_list, flux_list, colour, dashed)
  endif
  
  ! --- Get separatrices
  sep_list%n_psi = 1
  if ( xcase .eq. DOUBLE_NULL ) sep_list%n_psi = 2
  if (ES%active_xpoint .eq. SYMMETRIC_XPOINT) sep_list%n_psi = 1
  if (allocated(sep_list%psi_values)) call tr_deallocate(sep_list%psi_values,"sep_list%psi_values",CAT_GRID)
  call tr_allocate(sep_list%psi_values,1,sep_list%n_psi,"sep_list%psi_values",CAT_GRID)
  sep_list%psi_values(1) = flux_list%psi_values(n_flux)
  if ((xcase .eq. DOUBLE_NULL ) .and. (ES%active_xpoint .ne. SYMMETRIC_XPOINT)) sep_list%psi_values(2) = flux_list%psi_values(n_flux+n_open)
  call find_flux_surfaces(0,.true.,xcase,node_list,element_list,sep_list)  
  call reorder_flux_surfaces(node_list, element_list, sep_list, .false., ifail)
  if (debug .eq. 2) write(*,*) 'cleaning all surfaces'
  if ((xcase .eq. DOUBLE_NULL ) .and. (ES%active_xpoint .eq. SYMMETRIC_XPOINT)) then
    call get_symmetric_separatrix_contours(node_list, element_list, sep_list)
  else
    call get_separatrix_contours(node_list, element_list, sep_list)
  endif
  
  ! --- Loop over each core surface
  if (debug .eq. 2) write(*,*) 'cleaning core surfaces'
  do i_surf=1,n_flux-1
    i_part_save = 0
    do i_part = 1, flux_list%flux_surfaces(i_surf)%n_parts
      i_piece = 0.5 * (flux_list%flux_surfaces(i_surf)%parts_index(i_part) + flux_list%flux_surfaces(i_surf)%parts_index(i_part+1)-1)
      rr    = flux_list%flux_surfaces(i_surf)%s(1,i_piece)
      ss    = flux_list%flux_surfaces(i_surf)%t(1,i_piece)
      i_elm = flux_list%flux_surfaces(i_surf)%elm(i_piece)
      call interp_RZ(node_list,element_list,i_elm,rr,ss, &
                     R,dR_dr,dR_ds,dR_drs,dR_drr,dR_dss, &
                     Z,dZ_dr,dZ_ds,dZ_drs,dZ_drr,dZ_dss)
      call check_point_is_inside_contour(R, Z, n_separatrix_contour, R_separatrix_contour, Z_separatrix_contour, inside)
      if (inside .eq. 1) then
        i_part_save = i_part
        exit
      endif
    enddo
    count = 0
    i_part = i_part_save
    do i_piece = flux_list%flux_surfaces(i_surf)%parts_index(i_part),flux_list%flux_surfaces(i_surf)%parts_index(i_part+1)-1
      count = count + 1
      flux_list%flux_surfaces(i_surf)%elm(count) = flux_list%flux_surfaces(i_surf)%elm(i_piece)
      flux_list%flux_surfaces(i_surf)%s(1,count) = flux_list%flux_surfaces(i_surf)%s(1,i_piece)
      flux_list%flux_surfaces(i_surf)%s(2,count) = flux_list%flux_surfaces(i_surf)%s(2,i_piece)
      flux_list%flux_surfaces(i_surf)%s(3,count) = flux_list%flux_surfaces(i_surf)%s(3,i_piece)
      flux_list%flux_surfaces(i_surf)%s(4,count) = flux_list%flux_surfaces(i_surf)%s(4,i_piece)
      flux_list%flux_surfaces(i_surf)%t(1,count) = flux_list%flux_surfaces(i_surf)%t(1,i_piece)
      flux_list%flux_surfaces(i_surf)%t(2,count) = flux_list%flux_surfaces(i_surf)%t(2,i_piece)
      flux_list%flux_surfaces(i_surf)%t(3,count) = flux_list%flux_surfaces(i_surf)%t(3,i_piece)
      flux_list%flux_surfaces(i_surf)%t(4,count) = flux_list%flux_surfaces(i_surf)%t(4,i_piece)
    enddo
    flux_list%flux_surfaces(i_surf)%n_parts = 1
    flux_list%flux_surfaces(i_surf)%parts_index(2) = count+1
  enddo
  
  
  ! --- The main separatrix if not symmetric
  if ( (xcase .ne. DOUBLE_NULL) .or. (ES%active_xpoint .ne. SYMMETRIC_XPOINT) ) then
    if (debug .eq. 2) write(*,*) 'cleaning separatrix'
    i_surf = n_flux
    i_part_save = 0
    do i_part = 1, flux_list%flux_surfaces(i_surf)%n_parts
      i_piece = 0.5 * (flux_list%flux_surfaces(i_surf)%parts_index(i_part) + flux_list%flux_surfaces(i_surf)%parts_index(i_part+1)-1)
      rr    = flux_list%flux_surfaces(i_surf)%s(1,i_piece)
      ss    = flux_list%flux_surfaces(i_surf)%t(1,i_piece)
      i_elm = flux_list%flux_surfaces(i_surf)%elm(i_piece)
      call interp_RZ(node_list,element_list,i_elm,rr,ss, &
                     R,dR_dr,dR_ds,dR_drs,dR_drr,dR_dss, &
                     Z,dZ_dr,dZ_ds,dZ_drs,dZ_drr,dZ_dss)
      inside2 = 0
      if (ES%active_xpoint .eq. LOWER_XPOINT) then
        call check_point_is_inside_contour(R, Z, n_private_contour, R_private_contour, Z_private_contour, inside)
        if (xcase .eq. DOUBLE_NULL) call check_point_is_inside_contour(R, Z, n_up_priv_contour, R_up_priv_contour, Z_up_priv_contour, inside2)
      else
        call check_point_is_inside_contour(R, Z, n_up_priv_contour, R_up_priv_contour, Z_up_priv_contour, inside)
        if (xcase .eq. DOUBLE_NULL) call check_point_is_inside_contour(R, Z, n_private_contour, R_private_contour, Z_private_contour, inside2)
      endif
      if ( (inside .eq. 0) .and. (inside2 .eq. 0) ) then
        i_part_save = i_part
        exit
      endif
    enddo
    count = 0
    i_part = i_part_save
    do i_piece = flux_list%flux_surfaces(i_surf)%parts_index(i_part),flux_list%flux_surfaces(i_surf)%parts_index(i_part+1)-1
      count = count + 1
      flux_list%flux_surfaces(i_surf)%elm(count) = flux_list%flux_surfaces(i_surf)%elm(i_piece)
      flux_list%flux_surfaces(i_surf)%s(1,count) = flux_list%flux_surfaces(i_surf)%s(1,i_piece)
      flux_list%flux_surfaces(i_surf)%s(2,count) = flux_list%flux_surfaces(i_surf)%s(2,i_piece)
      flux_list%flux_surfaces(i_surf)%s(3,count) = flux_list%flux_surfaces(i_surf)%s(3,i_piece)
      flux_list%flux_surfaces(i_surf)%s(4,count) = flux_list%flux_surfaces(i_surf)%s(4,i_piece)
      flux_list%flux_surfaces(i_surf)%t(1,count) = flux_list%flux_surfaces(i_surf)%t(1,i_piece)
      flux_list%flux_surfaces(i_surf)%t(2,count) = flux_list%flux_surfaces(i_surf)%t(2,i_piece)
      flux_list%flux_surfaces(i_surf)%t(3,count) = flux_list%flux_surfaces(i_surf)%t(3,i_piece)
      flux_list%flux_surfaces(i_surf)%t(4,count) = flux_list%flux_surfaces(i_surf)%t(4,i_piece)
    enddo
    flux_list%flux_surfaces(i_surf)%n_parts = 1
    flux_list%flux_surfaces(i_surf)%parts_index(2) = count+1
  endif

  ! --- Loop over each sandwich surface
  if ( (xcase .eq. DOUBLE_NULL) .and. (ES%active_xpoint .ne. SYMMETRIC_XPOINT) ) then
    if (debug .eq. 2) write(*,*) 'cleaning sandwich surfaces'
    i_part_save = 0
    do i_surf=n_flux+1,n_flux+n_open-1
      do i_part = 1, flux_list%flux_surfaces(i_surf)%n_parts
        i_piece = 0.5 * (flux_list%flux_surfaces(i_surf)%parts_index(i_part) + flux_list%flux_surfaces(i_surf)%parts_index(i_part+1)-1)
        rr    = flux_list%flux_surfaces(i_surf)%s(1,i_piece)
        ss    = flux_list%flux_surfaces(i_surf)%t(1,i_piece)
        i_elm = flux_list%flux_surfaces(i_surf)%elm(i_piece)
        call interp_RZ(node_list,element_list,i_elm,rr,ss, &
                       R,dR_dr,dR_ds,dR_drs,dR_drr,dR_dss, &
                       Z,dZ_dr,dZ_ds,dZ_drs,dZ_drr,dZ_dss)
        if (ES%active_xpoint .eq. UPPER_XPOINT) then
          call check_point_is_inside_contour(R, Z, n_private_contour, R_private_contour, Z_private_contour, inside)
        else
          call check_point_is_inside_contour(R, Z, n_up_priv_contour, R_up_priv_contour, Z_up_priv_contour, inside)
        endif
        if (inside .eq. 0) then
          i_part_save = i_part
          exit
        endif
      enddo
      count = 0
      i_part = i_part_save
      do i_piece = flux_list%flux_surfaces(i_surf)%parts_index(i_part),flux_list%flux_surfaces(i_surf)%parts_index(i_part+1)-1
        count = count + 1
        flux_list%flux_surfaces(i_surf)%elm(count) = flux_list%flux_surfaces(i_surf)%elm(i_piece)
        flux_list%flux_surfaces(i_surf)%s(1,count) = flux_list%flux_surfaces(i_surf)%s(1,i_piece)
        flux_list%flux_surfaces(i_surf)%s(2,count) = flux_list%flux_surfaces(i_surf)%s(2,i_piece)
        flux_list%flux_surfaces(i_surf)%s(3,count) = flux_list%flux_surfaces(i_surf)%s(3,i_piece)
        flux_list%flux_surfaces(i_surf)%s(4,count) = flux_list%flux_surfaces(i_surf)%s(4,i_piece)
        flux_list%flux_surfaces(i_surf)%t(1,count) = flux_list%flux_surfaces(i_surf)%t(1,i_piece)
        flux_list%flux_surfaces(i_surf)%t(2,count) = flux_list%flux_surfaces(i_surf)%t(2,i_piece)
        flux_list%flux_surfaces(i_surf)%t(3,count) = flux_list%flux_surfaces(i_surf)%t(3,i_piece)
        flux_list%flux_surfaces(i_surf)%t(4,count) = flux_list%flux_surfaces(i_surf)%t(4,i_piece)
      enddo
      flux_list%flux_surfaces(i_surf)%n_parts = 1
      flux_list%flux_surfaces(i_surf)%parts_index(2) = count+1
    enddo
  ! --- Loop over each open surface
  else
    if (debug .eq. 2) write(*,*) 'cleaning sandwich surfaces'
    i_part_save = 0
    do i_surf=n_flux+1,n_flux+n_open-1
      i_pieces_max = 0
      do i_part = 1, flux_list%flux_surfaces(i_surf)%n_parts
        ! --- Yes, this is dirty, but it should work...
        i_piece = flux_list%flux_surfaces(i_surf)%parts_index(i_part+1)-1 - flux_list%flux_surfaces(i_surf)%parts_index(i_part)
        if (i_piece .gt. i_pieces_max) then
          i_pieces_max = i_piece
          i_part_save  = i_part
        endif
      enddo
      count = 0
      i_part = i_part_save
      do i_piece = flux_list%flux_surfaces(i_surf)%parts_index(i_part),flux_list%flux_surfaces(i_surf)%parts_index(i_part+1)-1
        count = count + 1
        flux_list%flux_surfaces(i_surf)%elm(count) = flux_list%flux_surfaces(i_surf)%elm(i_piece)
        flux_list%flux_surfaces(i_surf)%s(1,count) = flux_list%flux_surfaces(i_surf)%s(1,i_piece)
        flux_list%flux_surfaces(i_surf)%s(2,count) = flux_list%flux_surfaces(i_surf)%s(2,i_piece)
        flux_list%flux_surfaces(i_surf)%s(3,count) = flux_list%flux_surfaces(i_surf)%s(3,i_piece)
        flux_list%flux_surfaces(i_surf)%s(4,count) = flux_list%flux_surfaces(i_surf)%s(4,i_piece)
        flux_list%flux_surfaces(i_surf)%t(1,count) = flux_list%flux_surfaces(i_surf)%t(1,i_piece)
        flux_list%flux_surfaces(i_surf)%t(2,count) = flux_list%flux_surfaces(i_surf)%t(2,i_piece)
        flux_list%flux_surfaces(i_surf)%t(3,count) = flux_list%flux_surfaces(i_surf)%t(3,i_piece)
        flux_list%flux_surfaces(i_surf)%t(4,count) = flux_list%flux_surfaces(i_surf)%t(4,i_piece)
      enddo
      flux_list%flux_surfaces(i_surf)%n_parts = 1
      flux_list%flux_surfaces(i_surf)%parts_index(2) = count+1
    enddo
  endif
  
  
  ! --- Loop over each outer surface
  if (xcase .eq. DOUBLE_NULL) then
    if (debug .eq. 2) write(*,*) 'cleaning outer surfaces'
    i_part_save = 0
    do i_surf=n_flux+n_open+1,n_flux+n_open+n_outer
      do i_part = 1, flux_list%flux_surfaces(i_surf)%n_parts
        i_piece = 0.5 * (flux_list%flux_surfaces(i_surf)%parts_index(i_part) + flux_list%flux_surfaces(i_surf)%parts_index(i_part+1)-1)
        rr    = flux_list%flux_surfaces(i_surf)%s(1,i_piece)
        ss    = flux_list%flux_surfaces(i_surf)%t(1,i_piece)
        i_elm = flux_list%flux_surfaces(i_surf)%elm(i_piece)
        call interp_RZ(node_list,element_list,i_elm,rr,ss, &
                       R,dR_dr,dR_ds,dR_drs,dR_drr,dR_dss, &
                       Z,dZ_dr,dZ_ds,dZ_drs,dZ_drr,dZ_dss)
        call check_point_is_inside_contour(R, Z, n_outer_contour, R_outer_contour, Z_outer_contour, inside)
        if (inside .eq. 1) then
          i_part_save = i_part
          exit
        endif
      enddo
      count = 0
      i_part = i_part_save
      do i_piece = flux_list%flux_surfaces(i_surf)%parts_index(i_part),flux_list%flux_surfaces(i_surf)%parts_index(i_part+1)-1
        count = count + 1
        flux_list%flux_surfaces(i_surf)%elm(count) = flux_list%flux_surfaces(i_surf)%elm(i_piece)
        flux_list%flux_surfaces(i_surf)%s(1,count) = flux_list%flux_surfaces(i_surf)%s(1,i_piece)
        flux_list%flux_surfaces(i_surf)%s(2,count) = flux_list%flux_surfaces(i_surf)%s(2,i_piece)
        flux_list%flux_surfaces(i_surf)%s(3,count) = flux_list%flux_surfaces(i_surf)%s(3,i_piece)
        flux_list%flux_surfaces(i_surf)%s(4,count) = flux_list%flux_surfaces(i_surf)%s(4,i_piece)
        flux_list%flux_surfaces(i_surf)%t(1,count) = flux_list%flux_surfaces(i_surf)%t(1,i_piece)
        flux_list%flux_surfaces(i_surf)%t(2,count) = flux_list%flux_surfaces(i_surf)%t(2,i_piece)
        flux_list%flux_surfaces(i_surf)%t(3,count) = flux_list%flux_surfaces(i_surf)%t(3,i_piece)
        flux_list%flux_surfaces(i_surf)%t(4,count) = flux_list%flux_surfaces(i_surf)%t(4,i_piece)
      enddo
      flux_list%flux_surfaces(i_surf)%n_parts = 1
      flux_list%flux_surfaces(i_surf)%parts_index(2) = count+1
    enddo
  endif
  
  
  
  ! --- Loop over each inner surface
  if (xcase .eq. DOUBLE_NULL) then
    if (debug .eq. 2) write(*,*) 'cleaning inner surfaces'
    i_part_save = 0
    do i_surf=n_flux+n_open+n_outer+1,n_flux+n_open+n_outer+n_inner
      do i_part = 1, flux_list%flux_surfaces(i_surf)%n_parts
        i_piece = 0.5 * (flux_list%flux_surfaces(i_surf)%parts_index(i_part) + flux_list%flux_surfaces(i_surf)%parts_index(i_part+1)-1)
        rr    = flux_list%flux_surfaces(i_surf)%s(1,i_piece)
        ss    = flux_list%flux_surfaces(i_surf)%t(1,i_piece)
        i_elm = flux_list%flux_surfaces(i_surf)%elm(i_piece)
        call interp_RZ(node_list,element_list,i_elm,rr,ss, &
                       R,dR_dr,dR_ds,dR_drs,dR_drr,dR_dss, &
                       Z,dZ_dr,dZ_ds,dZ_drs,dZ_drr,dZ_dss)
        call check_point_is_inside_contour(R, Z, n_outer_contour, R_outer_contour, Z_outer_contour, inside)
        if (inside .eq. 0) then
          i_part_save = i_part
          exit
        endif
      enddo
      count = 0
      i_part = i_part_save
      do i_piece = flux_list%flux_surfaces(i_surf)%parts_index(i_part),flux_list%flux_surfaces(i_surf)%parts_index(i_part+1)-1
        count = count + 1
        flux_list%flux_surfaces(i_surf)%elm(count) = flux_list%flux_surfaces(i_surf)%elm(i_piece)
        flux_list%flux_surfaces(i_surf)%s(1,count) = flux_list%flux_surfaces(i_surf)%s(1,i_piece)
        flux_list%flux_surfaces(i_surf)%s(2,count) = flux_list%flux_surfaces(i_surf)%s(2,i_piece)
        flux_list%flux_surfaces(i_surf)%s(3,count) = flux_list%flux_surfaces(i_surf)%s(3,i_piece)
        flux_list%flux_surfaces(i_surf)%s(4,count) = flux_list%flux_surfaces(i_surf)%s(4,i_piece)
        flux_list%flux_surfaces(i_surf)%t(1,count) = flux_list%flux_surfaces(i_surf)%t(1,i_piece)
        flux_list%flux_surfaces(i_surf)%t(2,count) = flux_list%flux_surfaces(i_surf)%t(2,i_piece)
        flux_list%flux_surfaces(i_surf)%t(3,count) = flux_list%flux_surfaces(i_surf)%t(3,i_piece)
        flux_list%flux_surfaces(i_surf)%t(4,count) = flux_list%flux_surfaces(i_surf)%t(4,i_piece)
      enddo
      flux_list%flux_surfaces(i_surf)%n_parts = 1
      flux_list%flux_surfaces(i_surf)%parts_index(2) = count+1
    enddo
  endif
  
  
  ! --- Loop over each private surface
  i_part_save = 0
  if (debug .eq. 2) write(*,*) 'cleaning private surfaces'
  do i_surf=n_flux+n_open+n_outer+n_inner+1,n_flux+n_open+n_outer+n_inner+n_private
    do i_part = 1, flux_list%flux_surfaces(i_surf)%n_parts
      i_piece = 0.5 * (flux_list%flux_surfaces(i_surf)%parts_index(i_part) + flux_list%flux_surfaces(i_surf)%parts_index(i_part+1)-1)
      rr    = flux_list%flux_surfaces(i_surf)%s(1,i_piece)
      ss    = flux_list%flux_surfaces(i_surf)%t(1,i_piece)
      i_elm = flux_list%flux_surfaces(i_surf)%elm(i_piece)
      call interp_RZ(node_list,element_list,i_elm,rr,ss, &
                     R,dR_dr,dR_ds,dR_drs,dR_drr,dR_dss, &
                     Z,dZ_dr,dZ_ds,dZ_drs,dZ_drr,dZ_dss)
      call check_point_is_inside_contour(R, Z, n_private_contour, R_private_contour, Z_private_contour, inside)
      if (inside .eq. 1) then
        i_part_save = i_part
        exit
      endif
    enddo
    count = 0
    i_part = i_part_save
    do i_piece = flux_list%flux_surfaces(i_surf)%parts_index(i_part),flux_list%flux_surfaces(i_surf)%parts_index(i_part+1)-1
      count = count + 1
      flux_list%flux_surfaces(i_surf)%elm(count) = flux_list%flux_surfaces(i_surf)%elm(i_piece)
      flux_list%flux_surfaces(i_surf)%s(1,count) = flux_list%flux_surfaces(i_surf)%s(1,i_piece)
      flux_list%flux_surfaces(i_surf)%s(2,count) = flux_list%flux_surfaces(i_surf)%s(2,i_piece)
      flux_list%flux_surfaces(i_surf)%s(3,count) = flux_list%flux_surfaces(i_surf)%s(3,i_piece)
      flux_list%flux_surfaces(i_surf)%s(4,count) = flux_list%flux_surfaces(i_surf)%s(4,i_piece)
      flux_list%flux_surfaces(i_surf)%t(1,count) = flux_list%flux_surfaces(i_surf)%t(1,i_piece)
      flux_list%flux_surfaces(i_surf)%t(2,count) = flux_list%flux_surfaces(i_surf)%t(2,i_piece)
      flux_list%flux_surfaces(i_surf)%t(3,count) = flux_list%flux_surfaces(i_surf)%t(3,i_piece)
      flux_list%flux_surfaces(i_surf)%t(4,count) = flux_list%flux_surfaces(i_surf)%t(4,i_piece)
    enddo
    flux_list%flux_surfaces(i_surf)%n_parts = 1
    flux_list%flux_surfaces(i_surf)%parts_index(2) = count+1
  enddo
  
  ! --- Loop over each private surface
  if (xcase .eq. DOUBLE_NULL) then
    if (debug .eq. 2) write(*,*) 'cleaning upper private surfaces'
    i_part_save = 0
    do i_surf=n_flux+n_open+n_outer+n_inner+n_private+1,n_flux+n_open+n_outer+n_inner+n_private+n_up_priv
      do i_part = 1, flux_list%flux_surfaces(i_surf)%n_parts
        i_piece = 0.5 * (flux_list%flux_surfaces(i_surf)%parts_index(i_part) + flux_list%flux_surfaces(i_surf)%parts_index(i_part+1)-1)
        rr    = flux_list%flux_surfaces(i_surf)%s(1,i_piece)
        ss    = flux_list%flux_surfaces(i_surf)%t(1,i_piece)
        i_elm = flux_list%flux_surfaces(i_surf)%elm(i_piece)
        call interp_RZ(node_list,element_list,i_elm,rr,ss, &
                       R,dR_dr,dR_ds,dR_drs,dR_drr,dR_dss, &
                       Z,dZ_dr,dZ_ds,dZ_drs,dZ_drr,dZ_dss)
        call check_point_is_inside_contour(R, Z, n_up_priv_contour, R_up_priv_contour, Z_up_priv_contour, inside)
        if (inside .eq. 1) then
          i_part_save = i_part
          exit
        endif
      enddo
      count = 0
      i_part = i_part_save
      do i_piece = flux_list%flux_surfaces(i_surf)%parts_index(i_part),flux_list%flux_surfaces(i_surf)%parts_index(i_part+1)-1
        count = count + 1
        flux_list%flux_surfaces(i_surf)%elm(count) = flux_list%flux_surfaces(i_surf)%elm(i_piece)
        flux_list%flux_surfaces(i_surf)%s(1,count) = flux_list%flux_surfaces(i_surf)%s(1,i_piece)
        flux_list%flux_surfaces(i_surf)%s(2,count) = flux_list%flux_surfaces(i_surf)%s(2,i_piece)
        flux_list%flux_surfaces(i_surf)%s(3,count) = flux_list%flux_surfaces(i_surf)%s(3,i_piece)
        flux_list%flux_surfaces(i_surf)%s(4,count) = flux_list%flux_surfaces(i_surf)%s(4,i_piece)
        flux_list%flux_surfaces(i_surf)%t(1,count) = flux_list%flux_surfaces(i_surf)%t(1,i_piece)
        flux_list%flux_surfaces(i_surf)%t(2,count) = flux_list%flux_surfaces(i_surf)%t(2,i_piece)
        flux_list%flux_surfaces(i_surf)%t(3,count) = flux_list%flux_surfaces(i_surf)%t(3,i_piece)
        flux_list%flux_surfaces(i_surf)%t(4,count) = flux_list%flux_surfaces(i_surf)%t(4,i_piece)
      enddo
      flux_list%flux_surfaces(i_surf)%n_parts = 1
      flux_list%flux_surfaces(i_surf)%parts_index(2) = count+1
    enddo
  endif
  
  if (debug .eq. 2) then
    write(*,'(A,6i3)')'Checking cleaned surface parts',n_flux,n_open,n_outer,n_inner,n_private,n_up_priv
    ! --- Core
    write(*,*)'Core surfaces'
    do i_surf=1,n_flux-1
      write(*,*)i_surf,flux_list%flux_surfaces(i_surf)%n_parts
    enddo
    ! --- Separatrix
    write(*,*)'Main separatrix'
    i_surf=n_flux
    write(*,*)i_surf,flux_list%flux_surfaces(i_surf)%n_parts
    ! --- Open/sandwitch
    if (xcase .eq. DOUBLE_NULL) then
      write(*,*)'Sandwitch surfaces'
    else
      write(*,*)'Open surfaces'
    endif
    do i_surf=n_flux+1,n_flux+n_open-1
      write(*,*)i_surf,flux_list%flux_surfaces(i_surf)%n_parts
    enddo
    ! --- Last Open/2nd separatrix
    if (xcase .eq. DOUBLE_NULL) then
      write(*,*)'Second separatrix surface'
    else
      write(*,*)'Last open surface'
    endif
    i_surf=n_flux+n_open
    write(*,*)i_surf,flux_list%flux_surfaces(i_surf)%n_parts
    ! --- Outer/Inner
    if (xcase .eq. DOUBLE_NULL) then
      write(*,*)'Outer surface'
      do i_surf=n_flux+n_open+1,n_flux+n_open+n_outer
        write(*,*)i_surf,flux_list%flux_surfaces(i_surf)%n_parts
      enddo
      write(*,*)'Inner surface'
      do i_surf=n_flux+n_open+n_outer+1,n_flux+n_open+n_outer+n_inner
        write(*,*)i_surf,flux_list%flux_surfaces(i_surf)%n_parts
      enddo
    endif
    ! --- Private
    write(*,*)'Private surface'
    do i_surf=n_flux+n_open+n_outer+n_inner+1,n_flux+n_open+n_outer+n_inner+n_private
      write(*,*)i_surf,flux_list%flux_surfaces(i_surf)%n_parts
    enddo
    ! --- Upper private
    if (xcase .eq. DOUBLE_NULL) then
      write(*,*)'Upper private surface'
      do i_surf=n_flux+n_open+n_outer+n_inner+n_private+1,n_flux+n_open+n_outer+n_inner+n_private+n_up_priv
        write(*,*)i_surf,flux_list%flux_surfaces(i_surf)%n_parts
      enddo
    endif
  
  endif
  
  ! --- Finish plot of cleaned surfaces
  if (debug .ge. 1) then
    filename = 'plot_cleaned_flux_surfaces.py'
    colour = 'k' ; dashed = .true.
    call print_py_plot_ordered_flux_surfaces(filename, node_list, element_list, flux_list, colour, dashed)
    call print_py_plot_wall(filename)
    call print_py_plot_finish_plot(filename)
  endif
  
  return

end subroutine clean_surfaces









! ----------------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------------------------------


!> Deprecated!!! Use the one above, much more robust...
subroutine clean_single_surface(node_list,element_list,surface,location)

  use constants
  use data_structure
  use phys_module, only: xcase
  use grid_xpoint_data
  use mod_interp, only: interp_RZ
  use equil_info

  implicit none
  
  ! --- Routine parameters
  type (type_node_list),        intent(in)      :: node_list
  type (type_element_list),     intent(in)      :: element_list
  type (type_surface),          intent(inout)   :: surface
  integer,                      intent(in)      :: location
  
  ! --- Internal parameters
  type (type_surface)   :: surface_tmp
  integer               :: i
  integer               :: i_elm
  real*8                :: rr, ss
  real*8                :: R, Z
  
  surface_tmp%n_pieces = 0
  ! --- Check each piece
  do i=1,surface%n_pieces
    ! --- Get edge point of that surface piece
    rr    = surface%s(1,i)
    ss    = surface%t(1,i)
    i_elm = surface%elm(i)
    call interp_RZ(node_list,element_list,i_elm,rr,ss,R,Z)
    
    ! --- Core region (ie. not private parts)
    if (location .eq. core) then
      if ( (xcase .eq. LOWER_XPOINT) .and. (Z .gt. ES%Z_xpoint(1)) ) then
        surface_tmp%n_pieces = surface_tmp%n_pieces + 1
        surface_tmp%elm(surface_tmp%n_pieces) = surface%elm(i)
        surface_tmp%s(:,surface_tmp%n_pieces) = surface%s(:,i)
        surface_tmp%t(:,surface_tmp%n_pieces) = surface%t(:,i)
      endif
      if ( (xcase .eq. UPPER_XPOINT) .and. (Z .lt. ES%Z_xpoint(2)) ) then
        surface_tmp%n_pieces = surface_tmp%n_pieces + 1
        surface_tmp%elm(surface_tmp%n_pieces) = surface%elm(i)
        surface_tmp%s(:,surface_tmp%n_pieces) = surface%s(:,i)
        surface_tmp%t(:,surface_tmp%n_pieces) = surface%t(:,i)
      endif
      if ( (xcase .eq. DOUBLE_NULL ) .and. (Z .lt. ES%Z_xpoint(2)) .and. (Z .gt. ES%Z_xpoint(1)) ) then
        surface_tmp%n_pieces = surface_tmp%n_pieces + 1
        surface_tmp%elm(surface_tmp%n_pieces) = surface%elm(i)
        surface_tmp%s(:,surface_tmp%n_pieces) = surface%s(:,i)
        surface_tmp%t(:,surface_tmp%n_pieces) = surface%t(:,i)
      endif
    endif
      
    ! --- SOL region (always save)
    if (location .eq. SOL) then
      surface_tmp%n_pieces = surface_tmp%n_pieces + 1
      surface_tmp%elm(surface_tmp%n_pieces) = surface%elm(i)
      surface_tmp%s(:,surface_tmp%n_pieces) = surface%s(:,i)
      surface_tmp%t(:,surface_tmp%n_pieces) = surface%t(:,i)
    endif
      
    ! --- Sandwich region (ie. not private parts)
    if (location .eq. sandwich) then ! sandwich should be DOUBLE_NULL
      if ( (ES%active_xpoint .eq. UPPER_XPOINT) .and. (Z .gt. ES%Z_xpoint(1)) ) then
        surface_tmp%n_pieces = surface_tmp%n_pieces + 1
        surface_tmp%elm(surface_tmp%n_pieces) = surface%elm(i)
        surface_tmp%s(:,surface_tmp%n_pieces) = surface%s(:,i)
        surface_tmp%t(:,surface_tmp%n_pieces) = surface%t(:,i)
      endif
      if ( (ES%active_xpoint .eq. LOWER_XPOINT) .and. (Z .lt. ES%Z_xpoint(2)) ) then
        surface_tmp%n_pieces = surface_tmp%n_pieces + 1
        surface_tmp%elm(surface_tmp%n_pieces) = surface%elm(i)
        surface_tmp%s(:,surface_tmp%n_pieces) = surface%s(:,i)
        surface_tmp%t(:,surface_tmp%n_pieces) = surface%t(:,i)
      endif
    endif
    
    ! --- Outer region (ie. not inner parts)
    if (location .eq. outer) then
      if (R .gt. min(ES%R_xpoint(1),ES%R_xpoint(2))) then
        surface_tmp%n_pieces = surface_tmp%n_pieces + 1
        surface_tmp%elm(surface_tmp%n_pieces) = surface%elm(i)
        surface_tmp%s(:,surface_tmp%n_pieces) = surface%s(:,i)
        surface_tmp%t(:,surface_tmp%n_pieces) = surface%t(:,i)
      endif
    endif
          
    ! --- Inner region (ie. not outer parts)
    if (location .eq. inner) then
      if (R .lt. max(ES%R_xpoint(1),ES%R_xpoint(2))) then
        surface_tmp%n_pieces = surface_tmp%n_pieces + 1
        surface_tmp%elm(surface_tmp%n_pieces) = surface%elm(i)
        surface_tmp%s(:,surface_tmp%n_pieces) = surface%s(:,i)
        surface_tmp%t(:,surface_tmp%n_pieces) = surface%t(:,i)
      endif
    endif
          
    ! --- Private region (ie. not core parts)
    if (location .eq. private) then
      if ( (xcase .ne. UPPER_XPOINT) .and. (Z .lt. ES%Z_xpoint(1)) ) then
        surface_tmp%n_pieces = surface_tmp%n_pieces + 1
        surface_tmp%elm(surface_tmp%n_pieces) = surface%elm(i)
        surface_tmp%s(:,surface_tmp%n_pieces) = surface%s(:,i)
        surface_tmp%t(:,surface_tmp%n_pieces) = surface%t(:,i)
      endif
      if ( (xcase .eq. UPPER_XPOINT) .and. (Z .gt. ES%Z_xpoint(2)) ) then
        surface_tmp%n_pieces = surface_tmp%n_pieces + 1
        surface_tmp%elm(surface_tmp%n_pieces) = surface%elm(i)
        surface_tmp%s(:,surface_tmp%n_pieces) = surface%s(:,i)
        surface_tmp%t(:,surface_tmp%n_pieces) = surface%t(:,i)
      endif
    endif
    
    ! --- Upper private region (ie. not core parts)
    if (location .eq. upper_private) then
      if (Z .gt. ES%Z_xpoint(2)) then
        surface_tmp%n_pieces = surface_tmp%n_pieces + 1
        surface_tmp%elm(surface_tmp%n_pieces) = surface%elm(i)
        surface_tmp%s(:,surface_tmp%n_pieces) = surface%s(:,i)
        surface_tmp%t(:,surface_tmp%n_pieces) = surface%t(:,i)
      endif
    endif
    
  enddo
  
  ! --- Copy saved pieces only
  do i=1,surface_tmp%n_pieces
    surface%elm(i) = surface_tmp%elm(i) 
    surface%s(:,i) = surface_tmp%s(:,i)
    surface%t(:,i) = surface_tmp%t(:,i)
  enddo
  
  ! --- Make sure the rest is really empty...
  do i=surface_tmp%n_pieces+1,surface%n_pieces
    surface%elm(i) = 0
    surface%s(:,i) = (/ 0.d0, 0.d0, 0.d0, 0.d0  /)
    surface%t(:,i) = (/ 0.d0, 0.d0, 0.d0, 0.d0  /)
  enddo
  surface%n_pieces = surface_tmp%n_pieces

  
  return
  
end subroutine clean_single_surface






subroutine get_symmetric_separatrix_contours(node_list, element_list, sep_list)

  use constants
  use data_structure
  use phys_module, only : xcase
  use grid_xpoint_data
  use py_plots_grids
  use mod_interp, only: interp_RZ
  use equil_info
  
  implicit none
  
  ! --- Routine parameters
  type (type_node_list),        intent(in)      :: node_list
  type (type_element_list),     intent(in)      :: element_list
  type (type_surface_list),     intent(inout)   :: sep_list
  
  ! --- Internal parameters
  integer               :: i_surf, i_part, i_piece, i_elm, ifail, n_sub, i_sub, count
  integer               :: i_part_1, i_part_2
  integer               :: i_dir_1,  i_dir_2
  integer               :: i_beg, i_end, dir
  real*8                :: rr, ss, st, dr_flux, ds_flux
  integer               :: main_xpoint, second_xpoint
  real*8                :: R,R1,R2,dR_dr,dR_ds,dR_drs,dR_drr,dR_dss
  real*8                :: Z,Z1,Z2,dZ_dr,dZ_ds,dZ_drs,dZ_drr,dZ_dss
  real*8                :: rr1 , ss1, drr1, dss1
  real*8                :: rr2 , ss2, drr2, dss2
  real*8                :: distance, distance_Xpoint
  integer               :: n_lim, index_lim(n_wall_max)
  real*8                :: R_beg, Z_beg, R_end, Z_end
  logical               :: start_stop, debug, dashed, reversed
  character*256         :: filename
  character*1           :: colour
  
  debug = .true.
  
  main_xpoint   = 1
  second_xpoint = 2
  
  n_sub = 6
  
  
  ! --- The separatrix contour
  i_surf   = 1
  i_part_1 = 0
  i_part_2 = 0
  i_dir_1  = 1
  i_dir_2  = 1
  ! --- Find x-point part
  do i_part = 1,sep_list%flux_surfaces(i_surf)%n_parts
    do i_piece = sep_list%flux_surfaces(i_surf)%parts_index(i_part),sep_list%flux_surfaces(i_surf)%parts_index(i_part+1)-1
      i_elm = sep_list%flux_surfaces(i_surf)%elm(i_piece)
      if (i_elm .eq. ES%i_elm_xpoint(main_xpoint)) then
        if (i_part_1 .eq. 0) then
          i_part_1 = i_part
          i_dir_1  = 1
        else
          i_part_2 = i_part
          i_dir_2  = -1
        endif
        exit
      endif
      if (i_elm .eq. ES%i_elm_xpoint(second_xpoint)) then
        if (i_part_1 .eq. 0) then
          i_part_1 = i_part
          i_dir_1  = -1
        else
          i_part_2 = i_part
          i_dir_2  = 1
        endif
        exit
      endif
    enddo
  enddo
  ! ---Get contour by stepping along
  start_stop = .false.
  n_separatrix_contour = 0
  do i_part = i_part_1,i_part_2,i_part_2-i_part_1
    dir = i_dir_1
    if (i_part .eq. i_part_2) dir = i_dir_2
    if (dir .eq. 1) then
      i_beg = sep_list%flux_surfaces(i_surf)%parts_index(i_part)
      i_end = sep_list%flux_surfaces(i_surf)%parts_index(i_part+1)-1
    else
      i_beg = sep_list%flux_surfaces(i_surf)%parts_index(i_part+1)-1
      i_end = sep_list%flux_surfaces(i_surf)%parts_index(i_part)
    endif
    do i_piece = i_beg, i_end, dir
      ! --- Splines
      i_elm = sep_list%flux_surfaces(i_surf)%elm(i_piece)
      rr1   = sep_list%flux_surfaces(i_surf)%s(1,i_piece);   ss1  = sep_list%flux_surfaces(i_surf)%t(1,i_piece)
      drr1  = sep_list%flux_surfaces(i_surf)%s(2,i_piece);   dss1 = sep_list%flux_surfaces(i_surf)%t(2,i_piece)
      rr2   = sep_list%flux_surfaces(i_surf)%s(3,i_piece);   ss2  = sep_list%flux_surfaces(i_surf)%t(3,i_piece)
      drr2  = sep_list%flux_surfaces(i_surf)%s(4,i_piece);   dss2 = sep_list%flux_surfaces(i_surf)%t(4,i_piece)
      ! --- Through the X-point the first time
      if (     ( (i_part .eq. i_part_1) .and. (i_elm .eq. ES%i_elm_xpoint(main_xpoint  )) ) &
          .or. ( (i_part .eq. i_part_2) .and. (i_elm .eq. ES%i_elm_xpoint(second_xpoint)) ) )then
        n_separatrix_contour = n_separatrix_contour + 1
        if (i_part .eq. i_part_1) then
          R_separatrix_contour(n_separatrix_contour) = ES%R_xpoint(main_xpoint)
          Z_separatrix_contour(n_separatrix_contour) = ES%Z_xpoint(main_xpoint)
        else
          R_separatrix_contour(n_separatrix_contour) = ES%R_xpoint(second_xpoint)
          Z_separatrix_contour(n_separatrix_contour) = ES%Z_xpoint(second_xpoint)
        endif
        if (dir .eq. 1) then
          rr  = sep_list%flux_surfaces(i_surf)%s(3,i_piece)
          ss  = sep_list%flux_surfaces(i_surf)%t(3,i_piece)
        else
          rr  = sep_list%flux_surfaces(i_surf)%s(1,i_piece)
          ss  = sep_list%flux_surfaces(i_surf)%t(1,i_piece)
        endif
        call interp_RZ(node_list,element_list,i_elm,rr,ss, &
                       R1,dR_dr,dR_ds,dR_drs,dR_drr,dR_dss, &
                       Z1,dZ_dr,dZ_ds,dZ_drs,dZ_drr,dZ_dss)
        if (i_part .eq. i_part_1) then
          distance_Xpoint = sqrt( (R1-ES%R_xpoint(main_xpoint  ))**2 + (Z1-ES%Z_xpoint(main_xpoint  ))**2 )
        else
          distance_Xpoint = sqrt( (R1-ES%R_xpoint(second_xpoint))**2 + (Z1-ES%Z_xpoint(second_xpoint))**2 )
        endif
        do i_sub = 1,n_sub
          st = -dir*1.d0 + dir*2.d0*real(i_sub-1)/real(n_sub)
          call CUB1D(rr1, drr1, rr2, drr2, st, rr, dr_flux)
          call CUB1D(ss1, dss1, ss2, dss2, st, ss, ds_flux)
          call interp_RZ(node_list,element_list,i_elm,rr,ss, &
                         R,dR_dr,dR_ds,dR_drs,dR_drr,dR_dss, &
                         Z,dZ_dr,dZ_ds,dZ_drs,dZ_drr,dZ_dss)
          distance = sqrt( (R-R1)**2 + (Z-Z1)**2 )
          if (distance .lt. distance_Xpoint) then
            n_separatrix_contour = n_separatrix_contour + 1
            R_separatrix_contour(n_separatrix_contour) = R
            Z_separatrix_contour(n_separatrix_contour) = Z
          endif
        enddo
        start_stop = .true.
      ! --- Through the X-point the second time
      elseif (     ( (i_part .eq. i_part_1) .and. (i_elm .eq. ES%i_elm_xpoint(second_xpoint)) ) &
              .or. ( (i_part .eq. i_part_2) .and. (i_elm .eq. ES%i_elm_xpoint(main_xpoint  )) ) )then
        if (dir .eq. 1) then
          rr  = sep_list%flux_surfaces(i_surf)%s(3,i_piece)
          ss  = sep_list%flux_surfaces(i_surf)%t(3,i_piece)
        else
          rr  = sep_list%flux_surfaces(i_surf)%s(1,i_piece)
          ss  = sep_list%flux_surfaces(i_surf)%t(1,i_piece)
        endif
        call interp_RZ(node_list,element_list,i_elm,rr,ss, &
                       R1,dR_dr,dR_ds,dR_drs,dR_drr,dR_dss, &
                       Z1,dZ_dr,dZ_ds,dZ_drs,dZ_drr,dZ_dss)
        if (i_part .eq. i_part_1) then
          distance_Xpoint = sqrt( (R1-ES%R_xpoint(second_xpoint))**2 + (Z1-ES%Z_xpoint(second_xpoint))**2 )
        else
          distance_Xpoint = sqrt( (R1-ES%R_xpoint(main_xpoint  ))**2 + (Z1-ES%Z_xpoint(main_xpoint  ))**2 )
        endif
        do i_sub = 1,n_sub
          st = -dir*1.d0 + dir*2.d0*real(i_sub-1)/real(n_sub)
          call CUB1D(rr1, drr1, rr2, drr2, st, rr, dr_flux)
          call CUB1D(ss1, dss1, ss2, dss2, st, ss, ds_flux)
          call interp_RZ(node_list,element_list,i_elm,rr,ss, &
                         R,dR_dr,dR_ds,dR_drs,dR_drr,dR_dss, &
                         Z,dZ_dr,dZ_ds,dZ_drs,dZ_drr,dZ_dss)
          distance = sqrt( (R-R1)**2 + (Z-Z1)**2 )
          if (distance .gt. distance_Xpoint) then
            n_separatrix_contour = n_separatrix_contour + 1
            R_separatrix_contour(n_separatrix_contour) = R
            Z_separatrix_contour(n_separatrix_contour) = Z
          endif
        enddo
        if (i_part .eq. i_part_2) then
          n_separatrix_contour = n_separatrix_contour + 1
          R_separatrix_contour(n_separatrix_contour) = ES%R_xpoint(main_xpoint)
          Z_separatrix_contour(n_separatrix_contour) = ES%Z_xpoint(main_xpoint)
        endif
        start_stop = .false.
      ! --- Through non X-point pieces
      else
        if (start_stop) then
          do i_sub = 1,n_sub
            st = -dir*1.d0 + dir*2.d0*real(i_sub-1)/real(n_sub)
            call CUB1D(rr1, drr1, rr2, drr2, st, rr, dr_flux)
            call CUB1D(ss1, dss1, ss2, dss2, st, ss, ds_flux)
            call interp_RZ(node_list,element_list,i_elm,rr,ss, &
                           R,dR_dr,dR_ds,dR_drs,dR_drr,dR_dss, &
                           Z,dZ_dr,dZ_ds,dZ_drs,dZ_drr,dZ_dss)
            n_separatrix_contour = n_separatrix_contour + 1
            R_separatrix_contour(n_separatrix_contour) = R
            Z_separatrix_contour(n_separatrix_contour) = Z
          enddo
        endif
      endif
    enddo
  enddo
  
  
  
  ! --- The lower private contour
  i_surf   = 1
  i_part_1 = 0
  i_part_2 = 0
  i_dir_1  = 1
  i_dir_2  = 1
  ! --- Find x-point part
  do i_part = 1,sep_list%flux_surfaces(i_surf)%n_parts
    do i_piece = sep_list%flux_surfaces(i_surf)%parts_index(i_part),sep_list%flux_surfaces(i_surf)%parts_index(i_part+1)-1
      i_elm = sep_list%flux_surfaces(i_surf)%elm(i_piece)
      if (i_elm .eq. ES%i_elm_xpoint(main_xpoint)) then
        if (i_part_1 .eq. 0) then
          i_part_1 = i_part
          i_dir_1  = 1
        else
          i_part_2 = i_part
          i_dir_2  = -1
        endif
        exit
      endif
      if (i_elm .eq. ES%i_elm_xpoint(second_xpoint)) then
        if (i_part_1 .eq. 0) then
          i_part_1 = i_part
          i_dir_1  = -1
        else
          i_part_2 = i_part
          i_dir_2  = 1
        endif
        exit
      endif
    enddo
  enddo
  ! ---Get contour by stepping along
  start_stop = .true.
  n_private_contour = 0
  do i_part = i_part_1,i_part_2,i_part_2-i_part_1
    dir = i_dir_1
    if (i_part .eq. i_part_2) dir = i_dir_2
    if (dir .eq. 1) then
      i_beg = sep_list%flux_surfaces(i_surf)%parts_index(i_part)
      i_end = sep_list%flux_surfaces(i_surf)%parts_index(i_part+1)-1
    else
      i_beg = sep_list%flux_surfaces(i_surf)%parts_index(i_part+1)-1
      i_end = sep_list%flux_surfaces(i_surf)%parts_index(i_part)
    endif
    do i_piece = i_beg, i_end, dir
      ! --- Splines
      i_elm = sep_list%flux_surfaces(i_surf)%elm(i_piece)
      rr1   = sep_list%flux_surfaces(i_surf)%s(1,i_piece);   ss1  = sep_list%flux_surfaces(i_surf)%t(1,i_piece)
      drr1  = sep_list%flux_surfaces(i_surf)%s(2,i_piece);   dss1 = sep_list%flux_surfaces(i_surf)%t(2,i_piece)
      rr2   = sep_list%flux_surfaces(i_surf)%s(3,i_piece);   ss2  = sep_list%flux_surfaces(i_surf)%t(3,i_piece)
      drr2  = sep_list%flux_surfaces(i_surf)%s(4,i_piece);   dss2 = sep_list%flux_surfaces(i_surf)%t(4,i_piece)
      ! --- Through the X-point the first time
      if ( (i_part .eq. i_part_1) .and. (i_elm .eq. ES%i_elm_xpoint(main_xpoint)) ) then
        if (dir .eq. 1) then
          rr  = sep_list%flux_surfaces(i_surf)%s(3,i_piece)
          ss  = sep_list%flux_surfaces(i_surf)%t(3,i_piece)
        else
          rr  = sep_list%flux_surfaces(i_surf)%s(1,i_piece)
          ss  = sep_list%flux_surfaces(i_surf)%t(1,i_piece)
        endif
        call interp_RZ(node_list,element_list,i_elm,rr,ss, &
                       R1,dR_dr,dR_ds,dR_drs,dR_drr,dR_dss, &
                       Z1,dZ_dr,dZ_ds,dZ_drs,dZ_drr,dZ_dss)
        distance_Xpoint = sqrt( (R1-ES%R_xpoint(main_xpoint  ))**2 + (Z1-ES%Z_xpoint(main_xpoint  ))**2 )
        do i_sub = 1,n_sub
          st = -dir*1.d0 + dir*2.d0*real(i_sub-1)/real(n_sub)
          call CUB1D(rr1, drr1, rr2, drr2, st, rr, dr_flux)
          call CUB1D(ss1, dss1, ss2, dss2, st, ss, ds_flux)
          call interp_RZ(node_list,element_list,i_elm,rr,ss, &
                         R,dR_dr,dR_ds,dR_drs,dR_drr,dR_dss, &
                         Z,dZ_dr,dZ_ds,dZ_drs,dZ_drr,dZ_dss)
          distance = sqrt( (R-R1)**2 + (Z-Z1)**2 )
          if (distance .gt. distance_Xpoint) then
            n_private_contour = n_private_contour + 1
            R_private_contour(n_private_contour) = R
            Z_private_contour(n_private_contour) = Z
          endif
        enddo
        n_private_contour = n_private_contour + 1
        R_private_contour(n_private_contour) = ES%R_xpoint(main_xpoint)
        Z_private_contour(n_private_contour) = ES%Z_xpoint(main_xpoint)
        start_stop = .false.
      ! --- Through the X-point the second time
      elseif ( (i_part .eq. i_part_2) .and. (i_elm .eq. ES%i_elm_xpoint(main_xpoint)) ) then
        if (dir .eq. 1) then
          rr  = sep_list%flux_surfaces(i_surf)%s(3,i_piece)
          ss  = sep_list%flux_surfaces(i_surf)%t(3,i_piece)
        else
          rr  = sep_list%flux_surfaces(i_surf)%s(1,i_piece)
          ss  = sep_list%flux_surfaces(i_surf)%t(1,i_piece)
        endif
        call interp_RZ(node_list,element_list,i_elm,rr,ss, &
                       R1,dR_dr,dR_ds,dR_drs,dR_drr,dR_dss, &
                       Z1,dZ_dr,dZ_ds,dZ_drs,dZ_drr,dZ_dss)
        distance_Xpoint = sqrt( (R1-ES%R_xpoint(main_xpoint  ))**2 + (Z1-ES%Z_xpoint(main_xpoint  ))**2 )
        do i_sub = 1,n_sub
          st = -dir*1.d0 + dir*2.d0*real(i_sub-1)/real(n_sub)
          call CUB1D(rr1, drr1, rr2, drr2, st, rr, dr_flux)
          call CUB1D(ss1, dss1, ss2, dss2, st, ss, ds_flux)
          call interp_RZ(node_list,element_list,i_elm,rr,ss, &
                         R,dR_dr,dR_ds,dR_drs,dR_drr,dR_dss, &
                         Z,dZ_dr,dZ_ds,dZ_drs,dZ_drr,dZ_dss)
          distance = sqrt( (R-R1)**2 + (Z-Z1)**2 )
          if (distance .lt. distance_Xpoint) then
            n_private_contour = n_private_contour + 1
            R_private_contour(n_private_contour) = R
            Z_private_contour(n_private_contour) = Z
          endif
        enddo
        start_stop = .true.
      ! --- Through non X-point pieces
      else
        if (start_stop) then
          do i_sub = 1,n_sub
            st = -dir*1.d0 + dir*2.d0*real(i_sub-1)/real(n_sub)
            call CUB1D(rr1, drr1, rr2, drr2, st, rr, dr_flux)
            call CUB1D(ss1, dss1, ss2, dss2, st, ss, ds_flux)
            call interp_RZ(node_list,element_list,i_elm,rr,ss, &
                           R,dR_dr,dR_ds,dR_drs,dR_drr,dR_dss, &
                           Z,dZ_dr,dZ_ds,dZ_drs,dZ_drr,dZ_dss)
            n_private_contour = n_private_contour + 1
            R_private_contour(n_private_contour) = R
            Z_private_contour(n_private_contour) = Z
          enddo
        endif
      endif
    enddo
  enddo
  
  
  
  
  ! --- Close contour with the wall
  R_beg = R_private_contour(n_private_contour)
  Z_beg = Z_private_contour(n_private_contour)
  R_end = R_private_contour(1)
  Z_end = Z_private_contour(1)
  reversed = .false.
  call close_contour_with_wall(R_beg, Z_beg, R_end, Z_end, n_lim, index_lim, reversed)
  do i_piece = 1,n_lim
    n_private_contour = n_private_contour + 1
    R_private_contour(n_private_contour) = R_wall(index_lim(i_piece))
    Z_private_contour(n_private_contour) = Z_wall(index_lim(i_piece))
  enddo
  n_private_contour = n_private_contour + 1
  R_private_contour(n_private_contour) = R_private_contour(1)
  Z_private_contour(n_private_contour) = Z_private_contour(1)
  
  
  
  
  ! --- The upper private contour
  i_surf   = 1
  i_part_1 = 0
  i_part_2 = 0
  i_dir_1  = 1
  i_dir_2  = 1
  ! --- Find x-point part
  do i_part = 1,sep_list%flux_surfaces(i_surf)%n_parts
    do i_piece = sep_list%flux_surfaces(i_surf)%parts_index(i_part),sep_list%flux_surfaces(i_surf)%parts_index(i_part+1)-1
      i_elm = sep_list%flux_surfaces(i_surf)%elm(i_piece)
      if (i_elm .eq. ES%i_elm_xpoint(main_xpoint)) then
        if (i_part_1 .eq. 0) then
          i_part_1 = i_part
          i_dir_1  = -1
        else
          i_part_2 = i_part
          i_dir_2  = 1
        endif
        exit
      endif
      if (i_elm .eq. ES%i_elm_xpoint(second_xpoint)) then
        if (i_part_1 .eq. 0) then
          i_part_1 = i_part
          i_dir_1  = 1
        else
          i_part_2 = i_part
          i_dir_2  = -1
        endif
        exit
      endif
    enddo
  enddo
  ! ---Get contour by stepping along
  start_stop = .true.
  n_up_priv_contour = 0
  do i_part = i_part_1,i_part_2,i_part_2-i_part_1
    dir = i_dir_1
    if (i_part .eq. i_part_2) dir = i_dir_2
    if (dir .eq. 1) then
      i_beg = sep_list%flux_surfaces(i_surf)%parts_index(i_part)
      i_end = sep_list%flux_surfaces(i_surf)%parts_index(i_part+1)-1
    else
      i_beg = sep_list%flux_surfaces(i_surf)%parts_index(i_part+1)-1
      i_end = sep_list%flux_surfaces(i_surf)%parts_index(i_part)
    endif
    do i_piece = i_beg, i_end, dir
      ! --- Splines
      i_elm = sep_list%flux_surfaces(i_surf)%elm(i_piece)
      rr1   = sep_list%flux_surfaces(i_surf)%s(1,i_piece);   ss1  = sep_list%flux_surfaces(i_surf)%t(1,i_piece)
      drr1  = sep_list%flux_surfaces(i_surf)%s(2,i_piece);   dss1 = sep_list%flux_surfaces(i_surf)%t(2,i_piece)
      rr2   = sep_list%flux_surfaces(i_surf)%s(3,i_piece);   ss2  = sep_list%flux_surfaces(i_surf)%t(3,i_piece)
      drr2  = sep_list%flux_surfaces(i_surf)%s(4,i_piece);   dss2 = sep_list%flux_surfaces(i_surf)%t(4,i_piece)
      ! --- Through the X-point the first time
      if ( (i_part .eq. i_part_1) .and. (i_elm .eq. ES%i_elm_xpoint(second_xpoint)) ) then
        if (dir .eq. 1) then
          rr  = sep_list%flux_surfaces(i_surf)%s(3,i_piece)
          ss  = sep_list%flux_surfaces(i_surf)%t(3,i_piece)
        else
          rr  = sep_list%flux_surfaces(i_surf)%s(1,i_piece)
          ss  = sep_list%flux_surfaces(i_surf)%t(1,i_piece)
        endif
        call interp_RZ(node_list,element_list,i_elm,rr,ss, &
                       R1,dR_dr,dR_ds,dR_drs,dR_drr,dR_dss, &
                       Z1,dZ_dr,dZ_ds,dZ_drs,dZ_drr,dZ_dss)
        distance_Xpoint = sqrt( (R1-ES%R_xpoint(second_xpoint))**2 + (Z1-ES%Z_xpoint(second_xpoint))**2 )
        do i_sub = 1,n_sub
          st = -dir*1.d0 + dir*2.d0*real(i_sub-1)/real(n_sub)
          call CUB1D(rr1, drr1, rr2, drr2, st, rr, dr_flux)
          call CUB1D(ss1, dss1, ss2, dss2, st, ss, ds_flux)
          call interp_RZ(node_list,element_list,i_elm,rr,ss, &
                         R,dR_dr,dR_ds,dR_drs,dR_drr,dR_dss, &
                         Z,dZ_dr,dZ_ds,dZ_drs,dZ_drr,dZ_dss)
          distance = sqrt( (R-R1)**2 + (Z-Z1)**2 )
          if (distance .gt. distance_Xpoint) then
            n_up_priv_contour = n_up_priv_contour + 1
            R_up_priv_contour(n_up_priv_contour) = R
            Z_up_priv_contour(n_up_priv_contour) = Z
          endif
        enddo
        n_up_priv_contour = n_up_priv_contour + 1
        R_up_priv_contour(n_up_priv_contour) = ES%R_xpoint(second_xpoint)
        Z_up_priv_contour(n_up_priv_contour) = ES%Z_xpoint(second_xpoint)
        start_stop = .false.
      ! --- Through the X-point the second time
      elseif ( (i_part .eq. i_part_2) .and. (i_elm .eq. ES%i_elm_xpoint(second_xpoint)) ) then
        if (dir .eq. 1) then
          rr  = sep_list%flux_surfaces(i_surf)%s(3,i_piece)
          ss  = sep_list%flux_surfaces(i_surf)%t(3,i_piece)
        else
          rr  = sep_list%flux_surfaces(i_surf)%s(1,i_piece)
          ss  = sep_list%flux_surfaces(i_surf)%t(1,i_piece)
        endif
        call interp_RZ(node_list,element_list,i_elm,rr,ss, &
                       R1,dR_dr,dR_ds,dR_drs,dR_drr,dR_dss, &
                       Z1,dZ_dr,dZ_ds,dZ_drs,dZ_drr,dZ_dss)
        distance_Xpoint = sqrt( (R1-ES%R_xpoint(second_xpoint))**2 + (Z1-ES%Z_xpoint(second_xpoint))**2 )
        do i_sub = 1,n_sub
          st = -dir*1.d0 + dir*2.d0*real(i_sub-1)/real(n_sub)
          call CUB1D(rr1, drr1, rr2, drr2, st, rr, dr_flux)
          call CUB1D(ss1, dss1, ss2, dss2, st, ss, ds_flux)
          call interp_RZ(node_list,element_list,i_elm,rr,ss, &
                         R,dR_dr,dR_ds,dR_drs,dR_drr,dR_dss, &
                         Z,dZ_dr,dZ_ds,dZ_drs,dZ_drr,dZ_dss)
          distance = sqrt( (R-R1)**2 + (Z-Z1)**2 )
          if (distance .lt. distance_Xpoint) then
            n_up_priv_contour = n_up_priv_contour + 1
            R_up_priv_contour(n_up_priv_contour) = R
            Z_up_priv_contour(n_up_priv_contour) = Z
          endif
        enddo
        start_stop = .true.
      ! --- Through non X-point pieces
      else
        if (start_stop) then
          do i_sub = 1,n_sub
            st = -dir*1.d0 + dir*2.d0*real(i_sub-1)/real(n_sub)
            call CUB1D(rr1, drr1, rr2, drr2, st, rr, dr_flux)
            call CUB1D(ss1, dss1, ss2, dss2, st, ss, ds_flux)
            call interp_RZ(node_list,element_list,i_elm,rr,ss, &
                           R,dR_dr,dR_ds,dR_drs,dR_drr,dR_dss, &
                           Z,dZ_dr,dZ_ds,dZ_drs,dZ_drr,dZ_dss)
            n_up_priv_contour = n_up_priv_contour + 1
            R_up_priv_contour(n_up_priv_contour) = R
            Z_up_priv_contour(n_up_priv_contour) = Z
          enddo
        endif
      endif
    enddo
  enddo
  
  
  
  
  ! --- Close contour with the wall
  R_beg = R_up_priv_contour(n_up_priv_contour)
  Z_beg = Z_up_priv_contour(n_up_priv_contour)
  R_end = R_up_priv_contour(1)
  Z_end = Z_up_priv_contour(1)
  reversed = .false.
  call close_contour_with_wall(R_beg, Z_beg, R_end, Z_end, n_lim, index_lim, reversed)
  do i_piece = 1,n_lim
    n_up_priv_contour = n_up_priv_contour + 1
    R_up_priv_contour(n_up_priv_contour) = R_wall(index_lim(i_piece))
    Z_up_priv_contour(n_up_priv_contour) = Z_wall(index_lim(i_piece))
  enddo
  n_up_priv_contour = n_up_priv_contour + 1
  R_up_priv_contour(n_up_priv_contour) = R_up_priv_contour(1)
  Z_up_priv_contour(n_up_priv_contour) = Z_up_priv_contour(1)
  
  
  
  ! --- Finally the outer contour
  i_surf   = 1
  i_part_1 = 0
  i_part_2 = 0
  i_dir_1  = 1
  i_dir_2  = 1
  ! --- Find x-point part
  do i_part = 1,sep_list%flux_surfaces(i_surf)%n_parts
    do i_piece = sep_list%flux_surfaces(i_surf)%parts_index(i_part),sep_list%flux_surfaces(i_surf)%parts_index(i_part+1)-1
      i_elm = sep_list%flux_surfaces(i_surf)%elm(i_piece)
      if (i_elm .eq. ES%i_elm_xpoint(main_xpoint)) then
        if (i_part_1 .eq. 0) then
          i_part_1 = i_part
          i_dir_1  = 1
        else
          i_part_2 = i_part
          i_dir_2  = 1
        endif
        exit
      endif
      if (i_elm .eq. ES%i_elm_xpoint(second_xpoint)) then
        if (i_part_1 .eq. 0) then
          i_part_1 = i_part
          i_dir_1  = -1
        else
          i_part_2 = i_part
          i_dir_2  = -1
        endif
        exit
      endif
    enddo
  enddo
  ! --- Make sure we are on the right (outer) side
  i_part = i_part_1
  if (i_dir_1 .eq. 1) then
    i_piece = sep_list%flux_surfaces(i_surf)%parts_index(i_part)
  else
    i_piece = sep_list%flux_surfaces(i_surf)%parts_index(i_part+1)-1
  endif
  rr    = sep_list%flux_surfaces(i_surf)%s(1,i_piece)
  ss    = sep_list%flux_surfaces(i_surf)%t(1,i_piece)
  i_elm = sep_list%flux_surfaces(i_surf)%elm(i_piece)
  call interp_RZ(node_list,element_list,i_elm,rr,ss, &
                 R1,dR_dr,dR_ds,dR_drs,dR_drr,dR_dss, &
                 Z1,dZ_dr,dZ_ds,dZ_drs,dZ_drr,dZ_dss)
  i_part = i_part_2
  if (i_dir_2 .eq. 1) then
    i_piece = sep_list%flux_surfaces(i_surf)%parts_index(i_part)
  else
    i_piece = sep_list%flux_surfaces(i_surf)%parts_index(i_part+1)-1
  endif
  i_piece = sep_list%flux_surfaces(i_surf)%parts_index(i_part)
  rr    = sep_list%flux_surfaces(i_surf)%s(1,i_piece)
  ss    = sep_list%flux_surfaces(i_surf)%t(1,i_piece)
  i_elm = sep_list%flux_surfaces(i_surf)%elm(i_piece)
  call interp_RZ(node_list,element_list,i_elm,rr,ss, &
                 R2,dR_dr,dR_ds,dR_drs,dR_drr,dR_dss, &
                 Z2,dZ_dr,dZ_ds,dZ_drs,dZ_drr,dZ_dss)
  if (R2 .gt. R1) then
    i_part   = i_part_1
    i_part_1 = i_part_2
    i_part_2 = i_part
    dir     = i_dir_1
    i_dir_1 = i_dir_2
    i_dir_2 = dir
  endif
  ! ---Get contour by stepping along
  start_stop = .true.
  n_outer_contour = 0
  do i_part = i_part_1,i_part_2,i_part_2-i_part_1
    dir = i_dir_1
    if (i_part .eq. i_part_2) dir = i_dir_2
    if (dir .eq. 1) then
      i_beg = sep_list%flux_surfaces(i_surf)%parts_index(i_part)
      i_end = sep_list%flux_surfaces(i_surf)%parts_index(i_part+1)-1
    else
      i_beg = sep_list%flux_surfaces(i_surf)%parts_index(i_part+1)-1
      i_end = sep_list%flux_surfaces(i_surf)%parts_index(i_part)
    endif
    do i_piece = i_beg, i_end, dir
      ! --- Splines
      i_elm = sep_list%flux_surfaces(i_surf)%elm(i_piece)
      rr1   = sep_list%flux_surfaces(i_surf)%s(1,i_piece);   ss1  = sep_list%flux_surfaces(i_surf)%t(1,i_piece)
      drr1  = sep_list%flux_surfaces(i_surf)%s(2,i_piece);   dss1 = sep_list%flux_surfaces(i_surf)%t(2,i_piece)
      rr2   = sep_list%flux_surfaces(i_surf)%s(3,i_piece);   ss2  = sep_list%flux_surfaces(i_surf)%t(3,i_piece)
      drr2  = sep_list%flux_surfaces(i_surf)%s(4,i_piece);   dss2 = sep_list%flux_surfaces(i_surf)%t(4,i_piece)
      if (i_elm .eq. ES%i_elm_xpoint(main_xpoint)) then
        ! --- Then end point
        if (dir .eq. 1) then
          rr  = sep_list%flux_surfaces(i_surf)%s(3,i_piece)
          ss  = sep_list%flux_surfaces(i_surf)%t(3,i_piece)
        else
          rr  = sep_list%flux_surfaces(i_surf)%s(1,i_piece)
          ss  = sep_list%flux_surfaces(i_surf)%t(1,i_piece)
        endif
        i_elm = sep_list%flux_surfaces(i_surf)%elm(i_piece)
        call interp_RZ(node_list,element_list,i_elm,rr,ss, &
                       R1,dR_dr,dR_ds,dR_drs,dR_drr,dR_dss, &
                       Z1,dZ_dr,dZ_ds,dZ_drs,dZ_drr,dZ_dss)
        distance_Xpoint = sqrt( (R1-ES%R_xpoint(main_xpoint))**2 + (Z1-ES%Z_xpoint(main_xpoint))**2 )
        ! --- Through the X-point the first time
        if (i_part .eq. i_part_1) then
          do i_sub = 1,n_sub
            st = -dir*1.d0 + dir*2.d0*real(i_sub-1)/real(n_sub)
            call CUB1D(rr1, drr1, rr2, drr2, st, rr, dr_flux)
            call CUB1D(ss1, dss1, ss2, dss2, st, ss, ds_flux)
            call interp_RZ(node_list,element_list,i_elm,rr,ss, &
                           R,dR_dr,dR_ds,dR_drs,dR_drr,dR_dss, &
                           Z,dZ_dr,dZ_ds,dZ_drs,dZ_drr,dZ_dss)
            distance = sqrt( (R-R1)**2 + (Z-Z1)**2 )
            if (distance .gt. distance_Xpoint) then
              n_outer_contour = n_outer_contour + 1
              R_outer_contour(n_outer_contour) = R
              Z_outer_contour(n_outer_contour) = Z
            endif
          enddo
          n_outer_contour = n_outer_contour + 1
          R_outer_contour(n_outer_contour) = ES%R_xpoint(main_xpoint)
          Z_outer_contour(n_outer_contour) = ES%Z_xpoint(main_xpoint)
          start_stop = .false.
          exit
        else
          do i_sub = 1,n_sub
            st = -dir*1.d0 + dir*2.d0*real(i_sub-1)/real(n_sub)
            call CUB1D(rr1, drr1, rr2, drr2, st, rr, dr_flux)
            call CUB1D(ss1, dss1, ss2, dss2, st, ss, ds_flux)
            call interp_RZ(node_list,element_list,i_elm,rr,ss, &
                           R,dR_dr,dR_ds,dR_drs,dR_drr,dR_dss, &
                           Z,dZ_dr,dZ_ds,dZ_drs,dZ_drr,dZ_dss)
            distance = sqrt( (R-R1)**2 + (Z-Z1)**2 )
            if (distance .lt. distance_Xpoint) then
              n_outer_contour = n_outer_contour + 1
              R_outer_contour(n_outer_contour) = R
              Z_outer_contour(n_outer_contour) = Z
            endif
          enddo
          start_stop = .true.
        endif
      elseif (i_elm .eq. ES%i_elm_xpoint(second_xpoint)) then
        ! --- Through the X-point the second time
        do i_sub = 1,n_sub
          st = -dir*1.d0 + dir*2.d0*real(i_sub-1)/real(n_sub)
          call CUB1D(rr1, drr1, rr2, drr2, st, rr, dr_flux)
          call CUB1D(ss1, dss1, ss2, dss2, st, ss, ds_flux)
          call interp_RZ(node_list,element_list,i_elm,rr,ss, &
                         R,dR_dr,dR_ds,dR_drs,dR_drr,dR_dss, &
                         Z,dZ_dr,dZ_ds,dZ_drs,dZ_drr,dZ_dss)
          distance = sqrt( (R-R1)**2 + (Z-Z1)**2 )
          if (distance .gt. distance_Xpoint) then
            n_outer_contour = n_outer_contour + 1
            R_outer_contour(n_outer_contour) = R
            Z_outer_contour(n_outer_contour) = Z
          endif
        enddo
        n_outer_contour = n_outer_contour + 1
        R_outer_contour(n_outer_contour) = ES%R_xpoint(second_xpoint)
        Z_outer_contour(n_outer_contour) = ES%Z_xpoint(second_xpoint)
        start_stop = .false.
        exit
      ! --- Through non X-point pieces
      else
        if (start_stop) then
          do i_sub = 1,n_sub
            st = -dir*1.d0 + dir*2.d0*real(i_sub-1)/real(n_sub)
            call CUB1D(rr1, drr1, rr2, drr2, st, rr, dr_flux)
            call CUB1D(ss1, dss1, ss2, dss2, st, ss, ds_flux)
            call interp_RZ(node_list,element_list,i_elm,rr,ss, &
                           R,dR_dr,dR_ds,dR_drs,dR_drr,dR_dss, &
                           Z,dZ_dr,dZ_ds,dZ_drs,dZ_drr,dZ_dss)
            n_outer_contour = n_outer_contour + 1
            R_outer_contour(n_outer_contour) = R
            Z_outer_contour(n_outer_contour) = Z
          enddo
        endif
      endif
    enddo
  enddo
  ! --- Do it again for the upper part!
  i_part = i_part_1
  dir = i_dir_1
  if (dir .eq. 1) then
    i_beg = sep_list%flux_surfaces(i_surf)%parts_index(i_part)
    i_end = sep_list%flux_surfaces(i_surf)%parts_index(i_part+1)-1
  else
    i_beg = sep_list%flux_surfaces(i_surf)%parts_index(i_part+1)-1
    i_end = sep_list%flux_surfaces(i_surf)%parts_index(i_part)
  endif
  do i_piece = i_beg, i_end, dir
    ! --- Splines
    i_elm = sep_list%flux_surfaces(i_surf)%elm(i_piece)
    rr1   = sep_list%flux_surfaces(i_surf)%s(1,i_piece);   ss1  = sep_list%flux_surfaces(i_surf)%t(1,i_piece)
    drr1  = sep_list%flux_surfaces(i_surf)%s(2,i_piece);   dss1 = sep_list%flux_surfaces(i_surf)%t(2,i_piece)
    rr2   = sep_list%flux_surfaces(i_surf)%s(3,i_piece);   ss2  = sep_list%flux_surfaces(i_surf)%t(3,i_piece)
    drr2  = sep_list%flux_surfaces(i_surf)%s(4,i_piece);   dss2 = sep_list%flux_surfaces(i_surf)%t(4,i_piece)
    ! --- Through the upper X-point
    if (i_elm .eq. ES%i_elm_xpoint(second_xpoint)) then
      do i_sub = 1,n_sub
        st = -dir*1.d0 + dir*2.d0*real(i_sub-1)/real(n_sub)
        call CUB1D(rr1, drr1, rr2, drr2, st, rr, dr_flux)
        call CUB1D(ss1, dss1, ss2, dss2, st, ss, ds_flux)
        call interp_RZ(node_list,element_list,i_elm,rr,ss, &
                       R,dR_dr,dR_ds,dR_drs,dR_drr,dR_dss, &
                       Z,dZ_dr,dZ_ds,dZ_drs,dZ_drr,dZ_dss)
        distance = sqrt( (R-R1)**2 + (Z-Z1)**2 )
        if (distance .lt. distance_Xpoint) then
          n_outer_contour = n_outer_contour + 1
          R_outer_contour(n_outer_contour) = R
          Z_outer_contour(n_outer_contour) = Z
        endif
      enddo
      start_stop = .true.
    ! --- Through non X-point pieces
    else
      if (start_stop) then
        do i_sub = 1,n_sub
          st = -dir*1.d0 + dir*2.d0*real(i_sub-1)/real(n_sub)
          call CUB1D(rr1, drr1, rr2, drr2, st, rr, dr_flux)
          call CUB1D(ss1, dss1, ss2, dss2, st, ss, ds_flux)
          call interp_RZ(node_list,element_list,i_elm,rr,ss, &
                         R,dR_dr,dR_ds,dR_drs,dR_drr,dR_dss, &
                         Z,dZ_dr,dZ_ds,dZ_drs,dZ_drr,dZ_dss)
          n_outer_contour = n_outer_contour + 1
          R_outer_contour(n_outer_contour) = R
          Z_outer_contour(n_outer_contour) = Z
        enddo
      endif
    endif
  enddo
  
  
  
  ! --- Close contour with the wall
  R_beg = R_outer_contour(n_outer_contour)
  Z_beg = Z_outer_contour(n_outer_contour)
  R_end = R_outer_contour(1)
  Z_end = Z_outer_contour(1)
  reversed = .true.
  call close_contour_with_wall(R_beg, Z_beg, R_end, Z_end, n_lim, index_lim, reversed)
  do i_piece = 1,n_lim
    n_outer_contour = n_outer_contour + 1
    R_outer_contour(n_outer_contour) = R_wall(index_lim(i_piece))
    Z_outer_contour(n_outer_contour) = Z_wall(index_lim(i_piece))
  enddo
  n_outer_contour = n_outer_contour + 1
  R_outer_contour(n_outer_contour) = R_outer_contour(1)
  Z_outer_contour(n_outer_contour) = Z_outer_contour(1)
  
  
  
  
  
  
  if (debug) then
    filename = 'plot_closed_contours.py'
    call print_py_plot_prepare_plot(filename)
    call print_py_plot_wall(filename)
    colour = 'k' ;  dashed = .false.
    call print_py_plot_line(filename, n_separatrix_contour, R_separatrix_contour, Z_separatrix_contour, colour, dashed)
    colour = 'r' ;  dashed = .false.
    call print_py_plot_line(filename, n_private_contour, R_private_contour, Z_private_contour, colour, dashed)
    colour = 'g' ;  dashed = .false.
    call print_py_plot_line(filename, n_up_priv_contour, R_up_priv_contour, Z_up_priv_contour, colour, dashed)
    colour = 'c' ;  dashed = .true.
    call print_py_plot_line(filename, n_outer_contour, R_outer_contour, Z_outer_contour, colour, dashed)
    call print_py_plot_finish_plot(filename)
  endif
  
  return

end subroutine get_symmetric_separatrix_contours






subroutine get_separatrix_contours(node_list, element_list, sep_list)

  use constants
  use data_structure
  use phys_module, only : xcase
  use grid_xpoint_data
  use py_plots_grids
  use mod_interp, only: interp_RZ
  use equil_info
  
  implicit none
  
  ! --- Routine parameters
  type (type_node_list),        intent(in)      :: node_list
  type (type_element_list),     intent(in)      :: element_list
  type (type_surface_list),     intent(inout)   :: sep_list
  
  ! --- Internal parameters
  integer               :: i_surf, i_part, i_piece, i_elm, ifail, n_sub, i_sub
  integer               :: i_part_1, i_part_2
  integer               :: i_dir_1,  i_dir_2
  integer               :: i_beg, i_end, dir
  real*8                :: rr, ss, st, dr_flux, ds_flux
  integer               :: main_xpoint, second_xpoint
  integer               :: n_lim, index_lim(n_wall_max)
  real*8                :: R_beg, Z_beg, R_end, Z_end
  real*8                :: R,R1,R2,dR_dr,dR_ds,dR_drs,dR_drr,dR_dss
  real*8                :: Z,Z1,Z2,dZ_dr,dZ_ds,dZ_drs,dZ_drr,dZ_dss
  real*8                :: rr1 , ss1, drr1, dss1
  real*8                :: rr2 , ss2, drr2, dss2
  real*8                :: distance, distance_Xpoint
  logical               :: start_stop, debug, dashed, reversed
  integer               :: n_contour_tmp
  real*8                :: R_contour_tmp(n_contour_max), Z_contour_tmp(n_contour_max)
  character*256         :: filename
  character*1           :: colour
  
  debug = .true.
  
  main_xpoint   = 1
  second_xpoint = 2
  if ( (xcase .eq. DOUBLE_NULL) .and. (ES%active_xpoint .eq. UPPER_XPOINT) ) then
    main_xpoint   = 2
    second_xpoint = 1
  endif
  if (xcase .eq. UPPER_XPOINT) main_xpoint   = 2
  
  n_sub = 6
  
  ! --- First the separatrix contour
  i_surf   = 1
  i_part_1 = 0
  i_part_2 = 0
  i_dir_1  = 1
  i_dir_2  = 1
  ! --- Find x-point part
  do i_part = 1,sep_list%flux_surfaces(i_surf)%n_parts
    do i_piece = sep_list%flux_surfaces(i_surf)%parts_index(i_part),sep_list%flux_surfaces(i_surf)%parts_index(i_part+1)-1
      i_elm = sep_list%flux_surfaces(i_surf)%elm(i_piece)
      if (i_elm .eq. ES%i_elm_xpoint(main_xpoint)) then
        i_part_1 = i_part
        i_part_2 = i_part
      endif
    enddo
  enddo
  ! ---Get contour by stepping along
  i_part = i_part_1
  start_stop = .false.
  n_separatrix_contour = 0
  do i_piece = sep_list%flux_surfaces(i_surf)%parts_index(i_part),sep_list%flux_surfaces(i_surf)%parts_index(i_part+1)-1
    ! --- Splines
    i_elm = sep_list%flux_surfaces(i_surf)%elm(i_piece)
    rr1   = sep_list%flux_surfaces(i_surf)%s(1,i_piece);   ss1  = sep_list%flux_surfaces(i_surf)%t(1,i_piece)
    drr1  = sep_list%flux_surfaces(i_surf)%s(2,i_piece);   dss1 = sep_list%flux_surfaces(i_surf)%t(2,i_piece)
    rr2   = sep_list%flux_surfaces(i_surf)%s(3,i_piece);   ss2  = sep_list%flux_surfaces(i_surf)%t(3,i_piece)
    drr2  = sep_list%flux_surfaces(i_surf)%s(4,i_piece);   dss2 = sep_list%flux_surfaces(i_surf)%t(4,i_piece)
    if (i_elm .eq. ES%i_elm_xpoint(main_xpoint)) then
      ! --- Then end point
      rr    = sep_list%flux_surfaces(i_surf)%s(3,i_piece)
      ss    = sep_list%flux_surfaces(i_surf)%t(3,i_piece)
      i_elm = sep_list%flux_surfaces(i_surf)%elm(i_piece)
      call interp_RZ(node_list,element_list,i_elm,rr,ss, &
                     R1,dR_dr,dR_ds,dR_drs,dR_drr,dR_dss, &
                     Z1,dZ_dr,dZ_ds,dZ_drs,dZ_drr,dZ_dss)
      distance_Xpoint = sqrt( (R1-ES%R_xpoint(main_xpoint))**2 + (Z1-ES%Z_xpoint(main_xpoint))**2 )
      ! --- Through the X-point the first time
      if (.not. start_stop) then
        n_separatrix_contour = n_separatrix_contour + 1
        R_separatrix_contour(n_separatrix_contour) = ES%R_xpoint(main_xpoint)
        Z_separatrix_contour(n_separatrix_contour) = ES%Z_xpoint(main_xpoint)
        do i_sub = 1,n_sub
          st = -1.d0 + 2.d0*real(i_sub-1)/real(n_sub)
          call CUB1D(rr1, drr1, rr2, drr2, st, rr, dr_flux)
          call CUB1D(ss1, dss1, ss2, dss2, st, ss, ds_flux)
          call interp_RZ(node_list,element_list,i_elm,rr,ss, &
                         R,dR_dr,dR_ds,dR_drs,dR_drr,dR_dss, &
                         Z,dZ_dr,dZ_ds,dZ_drs,dZ_drr,dZ_dss)
          distance = sqrt( (R-R1)**2 + (Z-Z1)**2 )
          if (distance .lt. distance_Xpoint) then
            n_separatrix_contour = n_separatrix_contour + 1
            R_separatrix_contour(n_separatrix_contour) = R
            Z_separatrix_contour(n_separatrix_contour) = Z
          endif
        enddo
        start_stop = .true.
      ! --- Through the X-point the second time
      else
        do i_sub = 1,n_sub
          st = -1.d0 + 2.d0*real(i_sub-1)/real(n_sub)
          call CUB1D(rr1, drr1, rr2, drr2, st, rr, dr_flux)
          call CUB1D(ss1, dss1, ss2, dss2, st, ss, ds_flux)
          call interp_RZ(node_list,element_list,i_elm,rr,ss, &
                         R,dR_dr,dR_ds,dR_drs,dR_drr,dR_dss, &
                         Z,dZ_dr,dZ_ds,dZ_drs,dZ_drr,dZ_dss)
          distance = sqrt( (R-R1)**2 + (Z-Z1)**2 )
          if (distance .gt. distance_Xpoint) then
            n_separatrix_contour = n_separatrix_contour + 1
            R_separatrix_contour(n_separatrix_contour) = R
            Z_separatrix_contour(n_separatrix_contour) = Z
          endif
        enddo
        n_separatrix_contour = n_separatrix_contour + 1
        R_separatrix_contour(n_separatrix_contour) = ES%R_xpoint(main_xpoint)
        Z_separatrix_contour(n_separatrix_contour) = ES%Z_xpoint(main_xpoint)
        start_stop = .false.
      endif
    ! --- Through non X-point pieces
    else
      if (start_stop) then
        do i_sub = 1,n_sub
          st = -1.d0 + 2.d0*real(i_sub-1)/real(n_sub)
          call CUB1D(rr1, drr1, rr2, drr2, st, rr, dr_flux)
          call CUB1D(ss1, dss1, ss2, dss2, st, ss, ds_flux)
          call interp_RZ(node_list,element_list,i_elm,rr,ss, &
                         R,dR_dr,dR_ds,dR_drs,dR_drr,dR_dss, &
                         Z,dZ_dr,dZ_ds,dZ_drs,dZ_drr,dZ_dss)
          n_separatrix_contour = n_separatrix_contour + 1
          R_separatrix_contour(n_separatrix_contour) = R
          Z_separatrix_contour(n_separatrix_contour) = Z
        enddo
      endif
    endif
  enddo
  
  
  ! --- Then the second separatrix contour
  if (xcase .eq. DOUBLE_NULL) then
    i_surf   = 2
    i_part_1 = 0
    i_part_2 = 0
    ! --- Find x-point part
    do i_part = 1,sep_list%flux_surfaces(i_surf)%n_parts
      do i_piece = sep_list%flux_surfaces(i_surf)%parts_index(i_part),sep_list%flux_surfaces(i_surf)%parts_index(i_part+1)-1
        i_elm = sep_list%flux_surfaces(i_surf)%elm(i_piece)
        if (i_elm .eq. ES%i_elm_xpoint(second_xpoint)) then
          if (i_part_1 .eq. 0) then
            i_part_1 = i_part
          else
            i_part_2 = i_part
          endif
          exit
        endif
      enddo
    enddo
    ! --- Get directions
    i_dir_1  = 1
    i_dir_2  = 1
    do i_part = i_part_1,i_part_2,i_part_2-i_part_1
      i_piece = sep_list%flux_surfaces(i_surf)%parts_index(i_part)
      rr    = sep_list%flux_surfaces(i_surf)%s(1,i_piece)
      ss    = sep_list%flux_surfaces(i_surf)%t(1,i_piece)
      i_elm = sep_list%flux_surfaces(i_surf)%elm(i_piece)
      call interp_RZ(node_list,element_list,i_elm,rr,ss, &
                     R1,dR_dr,dR_ds,dR_drs,dR_drr,dR_dss, &
                     Z1,dZ_dr,dZ_ds,dZ_drs,dZ_drr,dZ_dss)
      i_piece = sep_list%flux_surfaces(i_surf)%parts_index(i_part+1)-1
      rr    = sep_list%flux_surfaces(i_surf)%s(3,i_piece)
      ss    = sep_list%flux_surfaces(i_surf)%t(3,i_piece)
      i_elm = sep_list%flux_surfaces(i_surf)%elm(i_piece)
      call interp_RZ(node_list,element_list,i_elm,rr,ss, &
                     R2,dR_dr,dR_ds,dR_drs,dR_drr,dR_dss, &
                     Z2,dZ_dr,dZ_ds,dZ_drs,dZ_drr,dZ_dss)
      if (second_xpoint .eq. 1) then
        if (i_part .eq. i_part_1) then
          if (Z2 .gt. Z1) i_dir_1 = -1
        else
          if (Z2 .lt. Z1) i_dir_2 = -1
        endif
      else
        if (i_part .eq. i_part_1) then
          if (Z2 .lt. Z1) i_dir_1 = -1
        else
          if (Z2 .gt. Z1) i_dir_2 = -1
        endif
      endif
    enddo
    ! ---Get contour by stepping along
    start_stop = .true.
    n_separatrix2_contour = 0
    do i_part = i_part_1,i_part_2,i_part_2-i_part_1
      dir = i_dir_1
      if (i_part .eq. i_part_2) dir = i_dir_2
      if (dir .eq. 1) then
        i_beg = sep_list%flux_surfaces(i_surf)%parts_index(i_part)
        i_end = sep_list%flux_surfaces(i_surf)%parts_index(i_part+1)-1
      else
        i_beg = sep_list%flux_surfaces(i_surf)%parts_index(i_part+1)-1
        i_end = sep_list%flux_surfaces(i_surf)%parts_index(i_part)
      endif
      do i_piece = i_beg, i_end, dir
        ! --- Splines
        i_elm = sep_list%flux_surfaces(i_surf)%elm(i_piece)
        rr1   = sep_list%flux_surfaces(i_surf)%s(1,i_piece);   ss1  = sep_list%flux_surfaces(i_surf)%t(1,i_piece)
        drr1  = sep_list%flux_surfaces(i_surf)%s(2,i_piece);   dss1 = sep_list%flux_surfaces(i_surf)%t(2,i_piece)
        rr2   = sep_list%flux_surfaces(i_surf)%s(3,i_piece);   ss2  = sep_list%flux_surfaces(i_surf)%t(3,i_piece)
        drr2  = sep_list%flux_surfaces(i_surf)%s(4,i_piece);   dss2 = sep_list%flux_surfaces(i_surf)%t(4,i_piece)
        if (i_elm .eq. ES%i_elm_xpoint(second_xpoint)) then
          ! --- Then end point
          if (dir .eq. 1) then
            rr  = sep_list%flux_surfaces(i_surf)%s(3,i_piece)
            ss  = sep_list%flux_surfaces(i_surf)%t(3,i_piece)
          else
            rr  = sep_list%flux_surfaces(i_surf)%s(1,i_piece)
            ss  = sep_list%flux_surfaces(i_surf)%t(1,i_piece)
          endif
          i_elm = sep_list%flux_surfaces(i_surf)%elm(i_piece)
          call interp_RZ(node_list,element_list,i_elm,rr,ss, &
                         R1,dR_dr,dR_ds,dR_drs,dR_drr,dR_dss, &
                         Z1,dZ_dr,dZ_ds,dZ_drs,dZ_drr,dZ_dss)
          distance_Xpoint = sqrt( (R1-ES%R_xpoint(second_xpoint))**2 + (Z1-ES%Z_xpoint(second_xpoint))**2 )
          ! --- Through the X-point the first time
          if (i_part .eq. i_part_1) then
            do i_sub = 1,n_sub
              st = -dir*1.d0 + dir*2.d0*real(i_sub-1)/real(n_sub)
              call CUB1D(rr1, drr1, rr2, drr2, st, rr, dr_flux)
              call CUB1D(ss1, dss1, ss2, dss2, st, ss, ds_flux)
              call interp_RZ(node_list,element_list,i_elm,rr,ss, &
                             R,dR_dr,dR_ds,dR_drs,dR_drr,dR_dss, &
                             Z,dZ_dr,dZ_ds,dZ_drs,dZ_drr,dZ_dss)
              distance = sqrt( (R-R1)**2 + (Z-Z1)**2 )
              if (distance .gt. distance_Xpoint) then
                n_separatrix2_contour = n_separatrix2_contour + 1
                R_separatrix2_contour(n_separatrix2_contour) = R
                Z_separatrix2_contour(n_separatrix2_contour) = Z
              endif
            enddo
            n_separatrix2_contour = n_separatrix2_contour + 1
            R_separatrix2_contour(n_separatrix2_contour) = ES%R_xpoint(second_xpoint)
            Z_separatrix2_contour(n_separatrix2_contour) = ES%Z_xpoint(second_xpoint)
            start_stop = .false.
          ! --- Through the X-point the second time
          else
            do i_sub = 1,n_sub
              st = -dir*1.d0 + dir*2.d0*real(i_sub-1)/real(n_sub)
              call CUB1D(rr1, drr1, rr2, drr2, st, rr, dr_flux)
              call CUB1D(ss1, dss1, ss2, dss2, st, ss, ds_flux)
              call interp_RZ(node_list,element_list,i_elm,rr,ss, &
                             R,dR_dr,dR_ds,dR_drs,dR_drr,dR_dss, &
                             Z,dZ_dr,dZ_ds,dZ_drs,dZ_drr,dZ_dss)
              distance = sqrt( (R-R1)**2 + (Z-Z1)**2 )
              if (distance .lt. distance_Xpoint) then
                n_separatrix2_contour = n_separatrix2_contour + 1
                R_separatrix2_contour(n_separatrix2_contour) = R
                Z_separatrix2_contour(n_separatrix2_contour) = Z
              endif
            enddo
            start_stop = .true.
          endif
        ! --- Through non X-point pieces
        else
          if (start_stop) then
            do i_sub = 1,n_sub
              st = -dir*1.d0 + dir*2.d0*real(i_sub-1)/real(n_sub)
              call CUB1D(rr1, drr1, rr2, drr2, st, rr, dr_flux)
              call CUB1D(ss1, dss1, ss2, dss2, st, ss, ds_flux)
              call interp_RZ(node_list,element_list,i_elm,rr,ss, &
                             R,dR_dr,dR_ds,dR_drs,dR_drr,dR_dss, &
                             Z,dZ_dr,dZ_ds,dZ_drs,dZ_drr,dZ_dss)
              n_separatrix2_contour = n_separatrix2_contour + 1
              R_separatrix2_contour(n_separatrix2_contour) = R
              Z_separatrix2_contour(n_separatrix2_contour) = Z
            enddo
          endif
        endif
      enddo
    enddo
    ! --- Close contour with the wall
    R_beg = R_separatrix2_contour(n_separatrix2_contour)
    Z_beg = Z_separatrix2_contour(n_separatrix2_contour)
    R_end = R_separatrix2_contour(1)
    Z_end = Z_separatrix2_contour(1)
    reversed = .false.
    call close_contour_with_wall(R_beg, Z_beg, R_end, Z_end, n_lim, index_lim, reversed)
    do i_piece = 1,n_lim
      n_separatrix2_contour = n_separatrix2_contour + 1
      R_separatrix2_contour(n_separatrix2_contour) = R_wall(index_lim(i_piece))
      Z_separatrix2_contour(n_separatrix2_contour) = Z_wall(index_lim(i_piece))
    enddo
    n_separatrix2_contour = n_separatrix2_contour + 1
    R_separatrix2_contour(n_separatrix2_contour) = R_separatrix2_contour(1)
    Z_separatrix2_contour(n_separatrix2_contour) = Z_separatrix2_contour(1)
  endif
  
  
  ! --- Now the main private contour
  i_surf   = 1
  i_part_1 = 0
  i_dir_1  = 1
  ! --- Find x-point part
  do i_part = 1,sep_list%flux_surfaces(i_surf)%n_parts
    do i_piece = sep_list%flux_surfaces(i_surf)%parts_index(i_part),sep_list%flux_surfaces(i_surf)%parts_index(i_part+1)-1
      i_elm = sep_list%flux_surfaces(i_surf)%elm(i_piece)
      if (i_elm .eq. ES%i_elm_xpoint(main_xpoint)) then
        i_part_1 = i_part
        exit
      endif
    enddo
  enddo
  ! ---Get contour by stepping along
  start_stop = .true.
  n_private_contour = 0
  i_part = i_part_1
  do i_piece = sep_list%flux_surfaces(i_surf)%parts_index(i_part),sep_list%flux_surfaces(i_surf)%parts_index(i_part+1)-1
    ! --- Splines
    i_elm = sep_list%flux_surfaces(i_surf)%elm(i_piece)
    rr1   = sep_list%flux_surfaces(i_surf)%s(1,i_piece);   ss1  = sep_list%flux_surfaces(i_surf)%t(1,i_piece)
    drr1  = sep_list%flux_surfaces(i_surf)%s(2,i_piece);   dss1 = sep_list%flux_surfaces(i_surf)%t(2,i_piece)
    rr2   = sep_list%flux_surfaces(i_surf)%s(3,i_piece);   ss2  = sep_list%flux_surfaces(i_surf)%t(3,i_piece)
    drr2  = sep_list%flux_surfaces(i_surf)%s(4,i_piece);   dss2 = sep_list%flux_surfaces(i_surf)%t(4,i_piece)
    if (i_elm .eq. ES%i_elm_xpoint(main_xpoint)) then
      ! --- Then end point
      rr  = sep_list%flux_surfaces(i_surf)%s(3,i_piece)
      ss  = sep_list%flux_surfaces(i_surf)%t(3,i_piece)
      i_elm = sep_list%flux_surfaces(i_surf)%elm(i_piece)
      call interp_RZ(node_list,element_list,i_elm,rr,ss, &
                     R1,dR_dr,dR_ds,dR_drs,dR_drr,dR_dss, &
                     Z1,dZ_dr,dZ_ds,dZ_drs,dZ_drr,dZ_dss)
      distance_Xpoint = sqrt( (R1-ES%R_xpoint(main_xpoint))**2 + (Z1-ES%Z_xpoint(main_xpoint))**2 )
      ! --- Through the X-point the first time
      if (start_stop) then
        do i_sub = 1,n_sub
          st = -1.d0 + 2.d0*real(i_sub-1)/real(n_sub)
          call CUB1D(rr1, drr1, rr2, drr2, st, rr, dr_flux)
          call CUB1D(ss1, dss1, ss2, dss2, st, ss, ds_flux)
          call interp_RZ(node_list,element_list,i_elm,rr,ss, &
                         R,dR_dr,dR_ds,dR_drs,dR_drr,dR_dss, &
                         Z,dZ_dr,dZ_ds,dZ_drs,dZ_drr,dZ_dss)
          distance = sqrt( (R-R1)**2 + (Z-Z1)**2 )
          if (distance .gt. distance_Xpoint) then
            n_private_contour = n_private_contour + 1
            R_private_contour(n_private_contour) = R
            Z_private_contour(n_private_contour) = Z
          endif
        enddo
        n_private_contour = n_private_contour + 1
        R_private_contour(n_private_contour) = ES%R_xpoint(main_xpoint)
        Z_private_contour(n_private_contour) = ES%Z_xpoint(main_xpoint)
        start_stop = .false.
      ! --- Through the X-point the second time
      else
        do i_sub = 1,n_sub
          st = -1.d0 + 2.d0*real(i_sub-1)/real(n_sub)
          call CUB1D(rr1, drr1, rr2, drr2, st, rr, dr_flux)
          call CUB1D(ss1, dss1, ss2, dss2, st, ss, ds_flux)
          call interp_RZ(node_list,element_list,i_elm,rr,ss, &
                         R,dR_dr,dR_ds,dR_drs,dR_drr,dR_dss, &
                         Z,dZ_dr,dZ_ds,dZ_drs,dZ_drr,dZ_dss)
          distance = sqrt( (R-R1)**2 + (Z-Z1)**2 )
          if (distance .lt. distance_Xpoint) then
            n_private_contour = n_private_contour + 1
            R_private_contour(n_private_contour) = R
            Z_private_contour(n_private_contour) = Z
          endif
        enddo
        start_stop = .true.
      endif
    ! --- Through non X-point pieces
    else
      if (start_stop) then
        do i_sub = 1,n_sub
          st = -1.d0 + 2.d0*real(i_sub-1)/real(n_sub)
          call CUB1D(rr1, drr1, rr2, drr2, st, rr, dr_flux)
          call CUB1D(ss1, dss1, ss2, dss2, st, ss, ds_flux)
          call interp_RZ(node_list,element_list,i_elm,rr,ss, &
                         R,dR_dr,dR_ds,dR_drs,dR_drr,dR_dss, &
                         Z,dZ_dr,dZ_ds,dZ_drs,dZ_drr,dZ_dss)
          n_private_contour = n_private_contour + 1
          R_private_contour(n_private_contour) = R
          Z_private_contour(n_private_contour) = Z
        enddo
      endif
    endif
  enddo
  ! --- Close contour with the wall
  R_beg = R_private_contour(n_private_contour)
  Z_beg = Z_private_contour(n_private_contour)
  R_end = R_private_contour(1)
  Z_end = Z_private_contour(1)
  reversed = .false.
  call close_contour_with_wall(R_beg, Z_beg, R_end, Z_end, n_lim, index_lim, reversed)
  do i_piece = 1,n_lim
    n_private_contour = n_private_contour + 1
    R_private_contour(n_private_contour) = R_wall(index_lim(i_piece))
    Z_private_contour(n_private_contour) = Z_wall(index_lim(i_piece))
  enddo
  n_private_contour = n_private_contour + 1
  R_private_contour(n_private_contour) = R_private_contour(1)
  Z_private_contour(n_private_contour) = Z_private_contour(1)
  
  
  
  ! --- Then the second private contour
  if (xcase .eq. DOUBLE_NULL) then
    i_surf   = 2
    i_part_1 = 0
    i_part_2 = 0
    ! --- Find x-point part
    do i_part = 1,sep_list%flux_surfaces(i_surf)%n_parts
      do i_piece = sep_list%flux_surfaces(i_surf)%parts_index(i_part),sep_list%flux_surfaces(i_surf)%parts_index(i_part+1)-1
        i_elm = sep_list%flux_surfaces(i_surf)%elm(i_piece)
        if (i_elm .eq. ES%i_elm_xpoint(second_xpoint)) then
          if (i_part_1 .eq. 0) then
            i_part_1 = i_part
          else
            i_part_2 = i_part
          endif
          exit
        endif
      enddo
    enddo
    ! --- Get directions
    i_dir_1  = 1
    i_dir_2  = 1
    do i_part = i_part_1,i_part_2,i_part_2-i_part_1
      i_piece = sep_list%flux_surfaces(i_surf)%parts_index(i_part)
      rr    = sep_list%flux_surfaces(i_surf)%s(1,i_piece)
      ss    = sep_list%flux_surfaces(i_surf)%t(1,i_piece)
      i_elm = sep_list%flux_surfaces(i_surf)%elm(i_piece)
      call interp_RZ(node_list,element_list,i_elm,rr,ss, &
                     R1,dR_dr,dR_ds,dR_drs,dR_drr,dR_dss, &
                     Z1,dZ_dr,dZ_ds,dZ_drs,dZ_drr,dZ_dss)
      i_piece = sep_list%flux_surfaces(i_surf)%parts_index(i_part+1)-1
      rr    = sep_list%flux_surfaces(i_surf)%s(3,i_piece)
      ss    = sep_list%flux_surfaces(i_surf)%t(3,i_piece)
      i_elm = sep_list%flux_surfaces(i_surf)%elm(i_piece)
      call interp_RZ(node_list,element_list,i_elm,rr,ss, &
                     R2,dR_dr,dR_ds,dR_drs,dR_drr,dR_dss, &
                     Z2,dZ_dr,dZ_ds,dZ_drs,dZ_drr,dZ_dss)
      if (second_xpoint .eq. 2) then
        if (i_part .eq. i_part_1) then
          if (Z2 .gt. Z1) i_dir_1 = -1
        else
          if (Z2 .lt. Z1) i_dir_2 = -1
        endif
      else
        if (i_part .eq. i_part_1) then
          if (Z2 .lt. Z1) i_dir_1 = -1
        else
          if (Z2 .gt. Z1) i_dir_2 = -1
        endif
      endif
    enddo
    ! ---Get contour by stepping along
    start_stop = .true.
    n_up_priv_contour = 0
    do i_part = i_part_1,i_part_2,i_part_2-i_part_1
      dir = i_dir_1
      if (i_part .eq. i_part_2) dir = i_dir_2
      if (dir .eq. 1) then
        i_beg = sep_list%flux_surfaces(i_surf)%parts_index(i_part)
        i_end = sep_list%flux_surfaces(i_surf)%parts_index(i_part+1)-1
      else
        i_beg = sep_list%flux_surfaces(i_surf)%parts_index(i_part+1)-1
        i_end = sep_list%flux_surfaces(i_surf)%parts_index(i_part)
      endif
      do i_piece = i_beg, i_end, dir
        ! --- Splines
        i_elm = sep_list%flux_surfaces(i_surf)%elm(i_piece)
        rr1   = sep_list%flux_surfaces(i_surf)%s(1,i_piece);   ss1  = sep_list%flux_surfaces(i_surf)%t(1,i_piece)
        drr1  = sep_list%flux_surfaces(i_surf)%s(2,i_piece);   dss1 = sep_list%flux_surfaces(i_surf)%t(2,i_piece)
        rr2   = sep_list%flux_surfaces(i_surf)%s(3,i_piece);   ss2  = sep_list%flux_surfaces(i_surf)%t(3,i_piece)
        drr2  = sep_list%flux_surfaces(i_surf)%s(4,i_piece);   dss2 = sep_list%flux_surfaces(i_surf)%t(4,i_piece)
        if (i_elm .eq. ES%i_elm_xpoint(second_xpoint)) then
          ! --- Then end point
          if (dir .eq. 1) then
            rr  = sep_list%flux_surfaces(i_surf)%s(3,i_piece)
            ss  = sep_list%flux_surfaces(i_surf)%t(3,i_piece)
          else
            rr  = sep_list%flux_surfaces(i_surf)%s(1,i_piece)
            ss  = sep_list%flux_surfaces(i_surf)%t(1,i_piece)
          endif
          i_elm = sep_list%flux_surfaces(i_surf)%elm(i_piece)
          call interp_RZ(node_list,element_list,i_elm,rr,ss, &
                         R1,dR_dr,dR_ds,dR_drs,dR_drr,dR_dss, &
                         Z1,dZ_dr,dZ_ds,dZ_drs,dZ_drr,dZ_dss)
          distance_Xpoint = sqrt( (R1-ES%R_xpoint(second_xpoint))**2 + (Z1-ES%Z_xpoint(second_xpoint))**2 )
          ! --- Through the X-point the first time
          if (i_part .eq. i_part_1) then
            do i_sub = 1,n_sub
              st = -dir*1.d0 + dir*2.d0*real(i_sub-1)/real(n_sub)
              call CUB1D(rr1, drr1, rr2, drr2, st, rr, dr_flux)
              call CUB1D(ss1, dss1, ss2, dss2, st, ss, ds_flux)
              call interp_RZ(node_list,element_list,i_elm,rr,ss, &
                             R,dR_dr,dR_ds,dR_drs,dR_drr,dR_dss, &
                             Z,dZ_dr,dZ_ds,dZ_drs,dZ_drr,dZ_dss)
              distance = sqrt( (R-R1)**2 + (Z-Z1)**2 )
              if (distance .gt. distance_Xpoint) then
                n_up_priv_contour = n_up_priv_contour + 1
                R_up_priv_contour(n_up_priv_contour) = R
                Z_up_priv_contour(n_up_priv_contour) = Z
              endif
            enddo
            n_up_priv_contour = n_up_priv_contour + 1
            R_up_priv_contour(n_up_priv_contour) = ES%R_xpoint(second_xpoint)
            Z_up_priv_contour(n_up_priv_contour) = ES%Z_xpoint(second_xpoint)
            start_stop = .false.
          ! --- Through the X-point the second time
          else
            do i_sub = 1,n_sub
              st = -dir*1.d0 + dir*2.d0*real(i_sub-1)/real(n_sub)
              call CUB1D(rr1, drr1, rr2, drr2, st, rr, dr_flux)
              call CUB1D(ss1, dss1, ss2, dss2, st, ss, ds_flux)
              call interp_RZ(node_list,element_list,i_elm,rr,ss, &
                             R,dR_dr,dR_ds,dR_drs,dR_drr,dR_dss, &
                             Z,dZ_dr,dZ_ds,dZ_drs,dZ_drr,dZ_dss)
              distance = sqrt( (R-R1)**2 + (Z-Z1)**2 )
              if (distance .lt. distance_Xpoint) then
                n_up_priv_contour = n_up_priv_contour + 1
                R_up_priv_contour(n_up_priv_contour) = R
                Z_up_priv_contour(n_up_priv_contour) = Z
              endif
            enddo
            start_stop = .true.
          endif
        ! --- Through non X-point pieces
        else
          if (start_stop) then
            do i_sub = 1,n_sub
              st = -dir*1.d0 + dir*2.d0*real(i_sub-1)/real(n_sub)
              call CUB1D(rr1, drr1, rr2, drr2, st, rr, dr_flux)
              call CUB1D(ss1, dss1, ss2, dss2, st, ss, ds_flux)
              call interp_RZ(node_list,element_list,i_elm,rr,ss, &
                             R,dR_dr,dR_ds,dR_drs,dR_drr,dR_dss, &
                             Z,dZ_dr,dZ_ds,dZ_drs,dZ_drr,dZ_dss)
              n_up_priv_contour = n_up_priv_contour + 1
              R_up_priv_contour(n_up_priv_contour) = R
              Z_up_priv_contour(n_up_priv_contour) = Z
            enddo
          endif
        endif
      enddo
    enddo
    ! --- Close contour with the wall
    R_beg = R_up_priv_contour(n_up_priv_contour)
    Z_beg = Z_up_priv_contour(n_up_priv_contour)
    R_end = R_up_priv_contour(1)
    Z_end = Z_up_priv_contour(1)
    reversed = .false.
    call close_contour_with_wall(R_beg, Z_beg, R_end, Z_end, n_lim, index_lim, reversed)
    do i_piece = 1,n_lim
      n_up_priv_contour = n_up_priv_contour + 1
      R_up_priv_contour(n_up_priv_contour) = R_wall(index_lim(i_piece))
      Z_up_priv_contour(n_up_priv_contour) = Z_wall(index_lim(i_piece))
    enddo
    n_up_priv_contour = n_up_priv_contour + 1
    R_up_priv_contour(n_up_priv_contour) = R_up_priv_contour(1)
    Z_up_priv_contour(n_up_priv_contour) = Z_up_priv_contour(1)
  endif
  
  
  ! --- Reverse the two private contours if upper-Xpoint is the main one
  if (ES%active_xpoint .eq. UPPER_XPOINT) then
    n_contour_tmp = n_up_priv_contour
    R_contour_tmp(1:n_up_priv_contour) = R_up_priv_contour(1:n_up_priv_contour)
    Z_contour_tmp(1:n_up_priv_contour) = Z_up_priv_contour(1:n_up_priv_contour)
    n_up_priv_contour = n_private_contour 
    R_up_priv_contour(1:n_private_contour) = R_private_contour(1:n_private_contour)
    Z_up_priv_contour(1:n_private_contour) = Z_private_contour(1:n_private_contour)
    n_private_contour = n_contour_tmp
    R_private_contour(1:n_contour_tmp) = R_contour_tmp(1:n_contour_tmp)
    Z_private_contour(1:n_contour_tmp) = Z_contour_tmp(1:n_contour_tmp)
  endif

  
  
  ! --- Finally the outer contour
  if (xcase .eq. DOUBLE_NULL) then
    i_surf   = 2
    i_part_1 = 0
    i_part_2 = 0
    ! --- Find x-point part
    do i_part = 1,sep_list%flux_surfaces(i_surf)%n_parts
      do i_piece = sep_list%flux_surfaces(i_surf)%parts_index(i_part),sep_list%flux_surfaces(i_surf)%parts_index(i_part+1)-1
        i_elm = sep_list%flux_surfaces(i_surf)%elm(i_piece)
        if (i_elm .eq. ES%i_elm_xpoint(second_xpoint)) then
          if (i_part_1 .eq. 0) then
            i_part_1 = i_part
          else
            i_part_2 = i_part
          endif
          exit
        endif
      enddo
    enddo
    ! --- Get directions
    i_dir_1  = 1
    i_dir_2  = 1
    do i_part = i_part_1,i_part_2,i_part_2-i_part_1
      i_piece = sep_list%flux_surfaces(i_surf)%parts_index(i_part)
      rr    = sep_list%flux_surfaces(i_surf)%s(1,i_piece)
      ss    = sep_list%flux_surfaces(i_surf)%t(1,i_piece)
      i_elm = sep_list%flux_surfaces(i_surf)%elm(i_piece)
      call interp_RZ(node_list,element_list,i_elm,rr,ss, &
                     R1,dR_dr,dR_ds,dR_drs,dR_drr,dR_dss, &
                     Z1,dZ_dr,dZ_ds,dZ_drs,dZ_drr,dZ_dss)
      i_piece = sep_list%flux_surfaces(i_surf)%parts_index(i_part+1)-1
      rr    = sep_list%flux_surfaces(i_surf)%s(3,i_piece)
      ss    = sep_list%flux_surfaces(i_surf)%t(3,i_piece)
      i_elm = sep_list%flux_surfaces(i_surf)%elm(i_piece)
      call interp_RZ(node_list,element_list,i_elm,rr,ss, &
                     R2,dR_dr,dR_ds,dR_drs,dR_drr,dR_dss, &
                     Z2,dZ_dr,dZ_ds,dZ_drs,dZ_drr,dZ_dss)
      if (Z2 .lt. Z1) then
        if (i_part .eq. i_part_1) then
          i_dir_1 = -1
        else
          i_dir_2 = -1
        endif
      endif
    enddo
    ! --- Make sure we are on the right (outer) side
    i_part = i_part_1
    if (i_dir_1 .eq. 1) then
      i_piece = sep_list%flux_surfaces(i_surf)%parts_index(i_part)
    else
      i_piece = sep_list%flux_surfaces(i_surf)%parts_index(i_part+1)-1
    endif
    rr    = sep_list%flux_surfaces(i_surf)%s(1,i_piece)
    ss    = sep_list%flux_surfaces(i_surf)%t(1,i_piece)
    i_elm = sep_list%flux_surfaces(i_surf)%elm(i_piece)
    call interp_RZ(node_list,element_list,i_elm,rr,ss, &
                   R1,dR_dr,dR_ds,dR_drs,dR_drr,dR_dss, &
                   Z1,dZ_dr,dZ_ds,dZ_drs,dZ_drr,dZ_dss)
    i_part = i_part_2
    if (i_dir_2 .eq. 1) then
      i_piece = sep_list%flux_surfaces(i_surf)%parts_index(i_part)
    else
      i_piece = sep_list%flux_surfaces(i_surf)%parts_index(i_part+1)-1
    endif
    i_piece = sep_list%flux_surfaces(i_surf)%parts_index(i_part)
    rr    = sep_list%flux_surfaces(i_surf)%s(3,i_piece)
    ss    = sep_list%flux_surfaces(i_surf)%t(3,i_piece)
    i_elm = sep_list%flux_surfaces(i_surf)%elm(i_piece)
    call interp_RZ(node_list,element_list,i_elm,rr,ss, &
                   R2,dR_dr,dR_ds,dR_drs,dR_drr,dR_dss, &
                   Z2,dZ_dr,dZ_ds,dZ_drs,dZ_drr,dZ_dss)
    if (R2 .gt. R1) then
      i_part   = i_part_1
      i_part_1 = i_part_2
      i_part_2 = i_part
      dir     = i_dir_1
      i_dir_1 = i_dir_2
      i_dir_2 = dir
    endif
    ! ---Get contour by stepping along
    start_stop = .true.
    n_outer_contour = 0
    do i_part = i_part_1,i_part_2,i_part_2-i_part_1
      dir = i_dir_1
      if (i_part .eq. i_part_2) dir = i_dir_2
      if (dir .eq. 1) then
        i_beg = sep_list%flux_surfaces(i_surf)%parts_index(i_part)
        i_end = sep_list%flux_surfaces(i_surf)%parts_index(i_part+1)-1
      else
        i_beg = sep_list%flux_surfaces(i_surf)%parts_index(i_part+1)-1
        i_end = sep_list%flux_surfaces(i_surf)%parts_index(i_part)
      endif
      do i_piece = i_beg, i_end, dir
        ! --- Splines
        i_elm = sep_list%flux_surfaces(i_surf)%elm(i_piece)
        rr1   = sep_list%flux_surfaces(i_surf)%s(1,i_piece);   ss1  = sep_list%flux_surfaces(i_surf)%t(1,i_piece)
        drr1  = sep_list%flux_surfaces(i_surf)%s(2,i_piece);   dss1 = sep_list%flux_surfaces(i_surf)%t(2,i_piece)
        rr2   = sep_list%flux_surfaces(i_surf)%s(3,i_piece);   ss2  = sep_list%flux_surfaces(i_surf)%t(3,i_piece)
        drr2  = sep_list%flux_surfaces(i_surf)%s(4,i_piece);   dss2 = sep_list%flux_surfaces(i_surf)%t(4,i_piece)
        if (i_elm .eq. ES%i_elm_xpoint(second_xpoint)) then
          ! --- Then end point
          if (dir .eq. 1) then
            rr  = sep_list%flux_surfaces(i_surf)%s(3,i_piece)
            ss  = sep_list%flux_surfaces(i_surf)%t(3,i_piece)
          else
            rr  = sep_list%flux_surfaces(i_surf)%s(1,i_piece)
            ss  = sep_list%flux_surfaces(i_surf)%t(1,i_piece)
          endif
          i_elm = sep_list%flux_surfaces(i_surf)%elm(i_piece)
          call interp_RZ(node_list,element_list,i_elm,rr,ss, &
                         R1,dR_dr,dR_ds,dR_drs,dR_drr,dR_dss, &
                         Z1,dZ_dr,dZ_ds,dZ_drs,dZ_drr,dZ_dss)
          distance_Xpoint = sqrt( (R1-ES%R_xpoint(second_xpoint))**2 + (Z1-ES%Z_xpoint(second_xpoint))**2 )
          ! --- Through the X-point the first time
          if (i_part .eq. i_part_1) then
            do i_sub = 1,n_sub
              st = -dir*1.d0 + dir*2.d0*real(i_sub-1)/real(n_sub)
              call CUB1D(rr1, drr1, rr2, drr2, st, rr, dr_flux)
              call CUB1D(ss1, dss1, ss2, dss2, st, ss, ds_flux)
              call interp_RZ(node_list,element_list,i_elm,rr,ss, &
                             R,dR_dr,dR_ds,dR_drs,dR_drr,dR_dss, &
                             Z,dZ_dr,dZ_ds,dZ_drs,dZ_drr,dZ_dss)
              distance = sqrt( (R-R1)**2 + (Z-Z1)**2 )
              if (distance .gt. distance_Xpoint) then
                n_outer_contour = n_outer_contour + 1
                R_outer_contour(n_outer_contour) = R
                Z_outer_contour(n_outer_contour) = Z
              endif
            enddo
            n_outer_contour = n_outer_contour + 1
            R_outer_contour(n_outer_contour) = ES%R_xpoint(second_xpoint)
            Z_outer_contour(n_outer_contour) = ES%Z_xpoint(second_xpoint)
            start_stop = .false.
          ! --- Through the X-point the second time
          else
            do i_sub = 1,n_sub
              st = -dir*1.d0 + dir*2.d0*real(i_sub-1)/real(n_sub)
              call CUB1D(rr1, drr1, rr2, drr2, st, rr, dr_flux)
              call CUB1D(ss1, dss1, ss2, dss2, st, ss, ds_flux)
              call interp_RZ(node_list,element_list,i_elm,rr,ss, &
                             R,dR_dr,dR_ds,dR_drs,dR_drr,dR_dss, &
                             Z,dZ_dr,dZ_ds,dZ_drs,dZ_drr,dZ_dss)
              distance = sqrt( (R-R1)**2 + (Z-Z1)**2 )
              if (distance .lt. distance_Xpoint) then
                n_outer_contour = n_outer_contour + 1
                R_outer_contour(n_outer_contour) = R
                Z_outer_contour(n_outer_contour) = Z
              endif
            enddo
            start_stop = .true.
          endif
        ! --- Through non X-point pieces
        else
          if (start_stop) then
            do i_sub = 1,n_sub
              st = -dir*1.d0 + dir*2.d0*real(i_sub-1)/real(n_sub)
              call CUB1D(rr1, drr1, rr2, drr2, st, rr, dr_flux)
              call CUB1D(ss1, dss1, ss2, dss2, st, ss, ds_flux)
              call interp_RZ(node_list,element_list,i_elm,rr,ss, &
                             R,dR_dr,dR_ds,dR_drs,dR_drr,dR_dss, &
                             Z,dZ_dr,dZ_ds,dZ_drs,dZ_drr,dZ_dss)
              n_outer_contour = n_outer_contour + 1
              R_outer_contour(n_outer_contour) = R
              Z_outer_contour(n_outer_contour) = Z
            enddo
          endif
        endif
      enddo
    enddo
    ! --- Close contour with the wall
    R_beg = R_outer_contour(n_outer_contour)
    Z_beg = Z_outer_contour(n_outer_contour)
    R_end = R_outer_contour(1)
    Z_end = Z_outer_contour(1)
    reversed = .true.
    call close_contour_with_wall(R_beg, Z_beg, R_end, Z_end, n_lim, index_lim, reversed)
    do i_piece = 1,n_lim
      n_outer_contour = n_outer_contour + 1
      R_outer_contour(n_outer_contour) = R_wall(index_lim(i_piece))
      Z_outer_contour(n_outer_contour) = Z_wall(index_lim(i_piece))
    enddo
    n_outer_contour = n_outer_contour + 1
    R_outer_contour(n_outer_contour) = R_outer_contour(1)
    Z_outer_contour(n_outer_contour) = Z_outer_contour(1)
  endif
  
  
  
  
  
  
  
  if (debug) then
    filename = 'plot_closed_contours.py'
    call print_py_plot_prepare_plot(filename)
    call print_py_plot_wall(filename)
    colour = 'k' ; dashed = .false.
    call print_py_plot_line(filename, n_separatrix_contour, R_separatrix_contour, Z_separatrix_contour, colour, dashed)
    colour = 'r' ; dashed = .false.
    call print_py_plot_line(filename, n_separatrix2_contour, R_separatrix2_contour, Z_separatrix2_contour, colour, dashed)
    colour = 'g' ; dashed = .true.
    call print_py_plot_line(filename, n_private_contour, R_private_contour, Z_private_contour, colour, dashed)
    colour = 'c' ; dashed = .false.
    call print_py_plot_line(filename, n_up_priv_contour, R_up_priv_contour, Z_up_priv_contour, colour, dashed)
    colour = 'm' ; dashed = .true.
    call print_py_plot_line(filename, n_outer_contour, R_outer_contour, Z_outer_contour, colour, dashed)
    call print_py_plot_finish_plot(filename)
  endif
  
  return
  
end subroutine get_separatrix_contours








subroutine close_contour_with_wall(R_beg, Z_beg, R_end, Z_end, n_lim, index_lim, reversed)

  use grid_xpoint_data, only: n_wall, R_wall, Z_wall, n_wall_max
  use phys_module, only: n_limiter, R_limiter, Z_limiter
  
  implicit none
  
  ! --- Routine parameters
  real*8,  intent(in   ) :: R_beg, Z_beg, R_end, Z_end
  integer, intent(inout) :: n_lim, index_lim(n_wall_max)
  logical, intent(in   ) :: reversed
  
  ! --- Internal parameters
  integer :: count, direction, n_tmp, i
  integer :: i_wall, i_lim_next, i_lim_beg, i_lim_end
  real*8  :: length1, length2
  real*8  :: diff_min_beg, diff_min_end, diff
  real*8  :: diff_R, diff_Z
  real*8  :: R_average1, R_average2
  real*8  :: Rmin, Rmax
  real*8  :: Zmin, Zmax
  real*8  :: offset
  
  ! --- Find out which wall points are our starting/ending points
  if (n_wall .eq. 0) then
    n_wall = n_limiter
    R_wall(1:n_wall) = R_limiter(1:n_wall)
    Z_wall(1:n_wall) = Z_limiter(1:n_wall)
    if (n_wall .eq. 0) then
      write(*,*)'Error getting wall data?'
      return
    endif
  endif
  diff_min_beg = 1.d10
  diff_min_end = 1.d10
  do i_wall = 1,n_wall
    diff = sqrt( (R_wall(i_wall)-R_beg)**2 +(Z_wall(i_wall)-Z_beg)**2 )
    if (diff .lt. diff_min_beg) then
      diff_min_beg = diff
      i_lim_beg = i_wall
    endif
    diff = sqrt( (R_wall(i_wall)-R_end)**2 +(Z_wall(i_wall)-Z_end)**2 )
    if (diff .lt. diff_min_end) then
      diff_min_end = diff
      i_lim_end = i_wall
    endif
  enddo
  
  ! --- Make sure we are going in right direction
  ! --- Get length in positive direction
  count = 1
  i_wall = i_lim_beg
  direction = +1
  length1 = 0.0
  R_average1 = 0.d0
  do i=1,n_wall
    i_lim_next = i_wall + direction
    if (i_lim_next .gt. n_wall) i_lim_next = 1
    if (i_lim_next .lt. 1     ) i_lim_next = n_wall
    count = count + 1
    length1 = length1 + sqrt( (R_wall(i_wall)-R_wall(i_lim_next))**2 + (Z_wall(i_wall)-Z_wall(i_lim_next))**2 )
    R_average1 = R_average1 + 0.5 * (R_wall(i_wall) + R_wall(i_lim_next))
    if (i_lim_next .eq. i_lim_end) exit
    i_wall = i_lim_next
  enddo
  n_tmp = count
  R_average1 = R_average1 / real(count)
  
  ! --- Then length in negative direction
  count = 1
  i_wall = i_lim_beg
  direction = -1
  length2 = 0.0
  R_average2 = 0.d0
  do i=1,n_wall
    i_lim_next = i_wall + direction
    if (i_lim_next .gt. n_wall) i_lim_next = 1
    if (i_lim_next .lt. 1     ) i_lim_next = n_wall
    count = count + 1
    length2 = length2 + sqrt( (R_wall(i_wall)-R_wall(i_lim_next))**2 + (Z_wall(i_wall)-Z_wall(i_lim_next))**2 )
    R_average2 = R_average2 + 0.5 * (R_wall(i_wall) + R_wall(i_lim_next))
    if (i_lim_next .eq. i_lim_end) exit
    i_wall = i_lim_next
  enddo
  R_average2 = R_average2 / real(count)
  if (length2 .lt. length1) then
    direction = -1
  else
    direction = +1
  endif
  
  ! --- And copy indices
  count = 1
  index_lim(1) = i_lim_beg
  i_wall = i_lim_beg
  do i=1,n_wall
    i_lim_next = i_wall + direction
    if (i_lim_next .gt. n_wall) i_lim_next = 1
    if (i_lim_next .lt. 1     ) i_lim_next = n_wall
    count = count + 1
    index_lim(count) = i_lim_next
    if (i_lim_next .eq. i_lim_end) exit
    i_wall = i_lim_next
  enddo
  n_lim = count
  
  ! --- If we want reversed, this means we want the outer wall on the right
  if (reversed) then
    if (R_average1 .lt. R_average2) then
      direction = -1
    else
      direction = +1
    endif
    ! --- And copy indices
    count = 1
    index_lim(1) = i_lim_beg
    i_wall = i_lim_beg
    do i=1,n_wall
      i_lim_next = i_wall + direction
      if (i_lim_next .gt. n_wall) i_lim_next = 1
      if (i_lim_next .lt. 1     ) i_lim_next = n_wall
      count = count + 1
      index_lim(count) = i_lim_next
      if (i_lim_next .eq. i_lim_end) exit
      i_wall = i_lim_next
    enddo
    n_lim = count
  endif
  
  ! --- Check that we don't make a u-turn at the beginning or the end
  diff_R = abs(R_wall(index_lim(1))-R_wall(index_lim(2)))
  diff_Z = abs(Z_wall(index_lim(1))-Z_wall(index_lim(2)))
  if (diff_R .gt. diff_Z) then
    offset = min(0.05*diff_R, 1.d-2)
    Rmin = min(R_wall(index_lim(1)),R_wall(index_lim(2))) + offset
    Rmax = max(R_wall(index_lim(1)),R_wall(index_lim(2))) - offset
    if ( (Rmin .lt. R_beg) .and. (R_beg .lt. Rmax) ) then
      do i = 1,n_lim-1
        index_lim(i) = index_lim(i+1)
      enddo
      n_lim = n_lim-1
    endif
  else
    offset = min(0.05*diff_Z, 1.d-2)
    Zmin = min(Z_wall(index_lim(1)),Z_wall(index_lim(2))) + offset
    Zmax = max(Z_wall(index_lim(1)),Z_wall(index_lim(2))) - offset
    if ( (Zmin .lt. Z_beg) .and. (Z_beg .lt. Zmax) ) then
      do i = 1,n_lim-1
        index_lim(i) = index_lim(i+1)
      enddo
      n_lim = n_lim-1
    endif
  endif
  
  diff_R = abs(R_wall(index_lim(n_lim))-R_wall(index_lim(n_lim-1)))
  diff_Z = abs(Z_wall(index_lim(n_lim))-Z_wall(index_lim(n_lim-1)))
  if (diff_R .gt. diff_Z) then
    offset = min(0.05*diff_R, 1.d-2)
    Rmin = min(R_wall(index_lim(n_lim)),R_wall(index_lim(n_lim-1))) + offset
    Rmax = max(R_wall(index_lim(n_lim)),R_wall(index_lim(n_lim-1))) - offset
    if ( (Rmin .lt. R_end) .and. (R_end .lt. Rmax) ) then
      n_lim = n_lim-1
    endif
  else
    offset = min(0.05*diff_Z, 1.d-2)
    Zmin = min(Z_wall(index_lim(n_lim)),Z_wall(index_lim(n_lim-1))) + offset
    Zmax = max(Z_wall(index_lim(n_lim)),Z_wall(index_lim(n_lim-1))) - offset
    if ( (Zmin .lt. Z_end) .and. (Z_end .lt. Zmax) ) then
      n_lim = n_lim-1
    endif
  endif
  
  return
end subroutine close_contour_with_wall



end module reorder_and_clean_flux_surfaces
