!> Diagnostic program for nice 3D visualizations of JOREK data with the raytracing program Povray.
!!
!! * Reads jorek_restart.rst|h5 and settings from namelist file povray.nml.
!! * Produces Povray script jorek_data.pov allowing to create nice 3D visualizations.
!! * Header for Povray script needs to be added manually (camera, light_sources...).
!! * Is already working but some functionality is missing and will be added soon.
!!
program jorek2_povray
  
  use constants
  use data_structure
  use nodes_elements
  use phys_module
  use diagnostics
  use profiles
  use mod_import_restart
  use mod_interp
  use basis_at_gaussian, only: initialise_basis
  
  implicit none
  
  integer, parameter :: FH = 53 !< File handle
  
  integer :: ierr, i_elm_axis, i_elm_xpoint(2), i_elm, ifail, i_phi, i_piece, ip, jp, my_id, error,&
    i_node
  real*8  :: psi_bnd, psi_axis, R_axis, Z_axis, s_axis, t_axis, psi_xpoint(2), R_xpoint(2),        &
    Z_xpoint(2), s_xpoint(2), t_xpoint(2), phi, ss1, dss1, ss2, dss2, tt1, dtt1, tt2, dtt2, u, si, &
    dsi, ti, dti, R, R_s, R_t, R_st, R_ss, R_tt, Z, Z_s, Z_t, Z_st, Z_ss, Z_tt, P, P_s, P_t, P_st, &
    P_ss, P_tt, min_val, max_val
  real*8, allocatable :: x_cart(:,:), y_cart(:,:), z_cart(:,:), value(:,:), x_cart2(:,:,:),        &
    y_cart2(:,:,:), z_cart2(:,:,:), value2(:,:,:)
  type(type_surface_list) :: surface_list
  
  ! --- Settings
  logical :: no_zero, no_private, symmetric, volume, cover
  integer :: i_var, n_phi, nplot, i_ct
  real*8  :: psi_surf, phi_start, phi_stop
  namelist / povray / no_zero, no_private, symmetric, volume, cover, i_var, n_phi, nplot, i_ct,    &
    psi_surf, phi_start, phi_stop
  
  
  
  ! --- Preset and read settings from namelist 'povray.nml'
  volume       = .false. ! Produce volume plot? Does not work yet!
  no_private   = .true.  ! Remove private flux region?
  no_zero      = .true.  ! Exclude axisymmetric component?
  symmetric    = .true.  ! Normalize symmetric around zero?
  cover        = .true.  ! Write out cover surfaces at toroidal positions phi_start and phi_stop?
  n_phi        = 30      ! Number of toroidal surfaces to write out
  nplot        = 2       ! Number of points per Bezier element (at least 2!)
  i_var        = 7       ! Plot values of this variable
  psi_surf     = 1.0d0   ! Plot values the flux surface corresponding to Psi_N=...
  phi_start    = 0.001d0 ! Start plot at this toroidal position
  phi_stop     = 5.2d0   ! Stop plot at this toroidal position
  i_ct         = 2       ! Use this color table
  open(FH, file='povray.nml', status='old', action='read', iostat=ierr)
  if ( ierr == 0 ) read(FH, povray, iostat=ierr)
  if ( ierr /= 0 ) then
    write(*,*) 'ERROR: Could not read namelist povray from povray.nml'
    stop
  end if
  close(FH)
  
  ! --- Initialize
  my_id = 0
  call initialise_parameters(my_id, "__NO_FILENAME__")
  call det_modes()
  allocate( x_cart(nplot,n_phi), y_cart(nplot,n_phi), z_cart(nplot,n_phi), value(nplot,n_phi),     &
    x_cart2(nplot,nplot,n_phi), y_cart2(nplot,nplot,n_phi), z_cart2(nplot,nplot,n_phi),            &
    value2(nplot,nplot,n_phi) )
  open(FH, file='jorek_data.pov', status='replace', action='write', form='formatted')
  call import_restart(node_list, element_list, 'jorek_restart', rst_format, ierr, .true.)
  call initialise_basis()
  
  ! --- Locate magnetic axis and X-point.
  call find_axis(0, node_list, element_list, psi_axis, R_axis, Z_axis, i_elm_axis, s_axis, t_axis, &
    ifail)
  if ( xpoint ) then
    call find_xpoint(0, node_list, element_list, psi_xpoint, R_xpoint, Z_xpoint, i_elm_xpoint,     &
      s_xpoint, t_xpoint, xcase, ifail)
    if ( (xcase == DOUBLE_NULL) ) then
      psi_bnd = minval( psi_xpoint(:) )
    else
      psi_bnd = psi_xpoint(1)
    end if
  else
    psi_bnd = 0.d0
  end if
  
  ! --- Determine minimum and maximum variable values (in phi=0 plane) for normalization
  min_val = +1.d99
  max_val = -1.d99
  phi     = 0.d0
  do i_node = 1, node_list%n_nodes
    P = node_list%node(i_node)%values(1,1,1)
    if ( (P-psi_axis)/(psi_bnd-psi_axis) > psi_surf ) cycle
    
    R = node_list%node(i_node)%x(1,1,1)
    Z = node_list%node(i_node)%x(1,1,2)
    P = variable_value(i_var, i_node, phi, no_zero)
    min_val = min( min_val, P )
    max_val = max( max_val, P )
  end do
  if ( symmetric ) then
    max_val = max( abs(min_val), abs(max_val) )
    min_val = - max_val
  end if
  write(*,*) 'Variable value range:', min_val, max_val
  
  ! --- Write out data
  if ( volume ) then
    call plot_volume()
  else
    call plot_surface()
  end if
  
  ! --- Clean up
  close(FH)
  
  
  
  contains
  
  
  
  !> Plot data volume-like; does not work yet properly!
  subroutine plot_volume()
    
    do i_elm = 1, element_list%n_elements
      if ( mod(i_elm, element_list%n_elements/40+1) == 0 ) then
        write(*,*) 'Element', i_elm, ' of', element_list%n_elements
      end if
      
      do ip = 1, nplot
        
        si = real(ip-1)/real(nplot-1)
        
        do jp = 1, nplot
          
          ti = real(jp-1)/real(nplot-1)
          
          call interp_RZ(node_list, element_list, i_elm, si, ti, R, R_s, R_t, R_st, R_ss, R_tt,    &
            Z, Z_s, Z_t, Z_st, Z_ss, Z_tt)
          
          do i_phi = 1, n_phi
            
            phi = phi_start + (phi_stop-phi_start)*real(i_phi-1)/real(n_phi-1)
            
            value2(ip,jp,i_phi) = variable_value(i_var, i_elm, si, ti, phi, no_zero)
            
            call interp (node_list, element_list, i_elm, 1, 1, si, ti, P, P_s, P_t, P_st, P_ss,    &
              P_tt)
            
            x_cart2(ip,jp,i_phi) = + R * cos(phi)
            y_cart2(ip,jp,i_phi) = - R * sin(phi)
            z_cart2(ip,jp,i_phi) = Z
            
            if ( (ip==1) .or. (jp==1) .or. (i_phi==1)                                              &
              .or. ((P-psi_axis)/(psi_bnd-psi_axis) > psi_surf)                                    &
              .or. ( no_private .and. (Z < Z_xpoint(1)) ) ) cycle 
            
            
            P = (value2(ip-0,jp-1,i_phi-1)+value2(ip-0,jp-1,i_phi-0)+value2(ip-0,jp-0,i_phi-0))/3.d0
            if ( abs(P) > 0.75*max_val ) &
            call write_triangle( FH,                                                               &
              x_cart2(ip-0,jp-1,i_phi-1), y_cart2(ip-0,jp-1,i_phi-1), z_cart2(ip-0,jp-1,i_phi-1),  &
              x_cart2(ip-0,jp-1,i_phi-0), y_cart2(ip-0,jp-1,i_phi-0), z_cart2(ip-0,jp-1,i_phi-0),  &
              x_cart2(ip-0,jp-0,i_phi-0), y_cart2(ip-0,jp-0,i_phi-0), z_cart2(ip-0,jp-0,i_phi-0),  &
              P )
            
            P = (value2(ip-0,jp-1,i_phi-1)+value2(ip-0,jp-0,i_phi-0)+value2(ip-0,jp-0,i_phi-1))/3.d0
            if ( abs(P) > 0.75*max_val ) &
            call write_triangle( FH,                                                               &
              x_cart2(ip-0,jp-1,i_phi-1), y_cart2(ip-0,jp-1,i_phi-1), z_cart2(ip-0,jp-1,i_phi-1),  &
              x_cart2(ip-0,jp-0,i_phi-0), y_cart2(ip-0,jp-0,i_phi-0), z_cart2(ip-0,jp-0,i_phi-0),  &
              x_cart2(ip-0,jp-0,i_phi-1), y_cart2(ip-0,jp-0,i_phi-1), z_cart2(ip-0,jp-0,i_phi-1),  &
              P )
            
          end do
        end do
      end do
      
    end do
    
  end subroutine plot_volume
  
  
  
  !> Plot data on a given flux surface
  subroutine plot_surface()
    
    ! --- Find the flux surface
    surface_list%n_psi = 1
    allocate(surface_list%psi_values(surface_list%n_psi))
    surface_list%psi_values(1) = psi_axis + psi_surf*(psi_bnd-psi_axis)
    call find_flux_surfaces(0,xpoint,xcase,node_list,element_list,surface_list)
    
    ! --- Loop over all segments of the flux surface
    do i_piece = 1, surface_list%flux_surfaces(1)%n_pieces
      if ( mod(i_piece, surface_list%flux_surfaces(1)%n_pieces/40+1) == 0 ) then
        write(*,*) 'Piece', i_piece, ' of', surface_list%flux_surfaces(1)%n_pieces
      end if
      
      ! --- Bezier element, in which the current flux surface segment is located
      i_elm = surface_list%flux_surfaces(1)%elm(i_piece)
      
      ss1  = surface_list%flux_surfaces(1)%s(1,i_piece)
      dss1 = surface_list%flux_surfaces(1)%s(2,i_piece)
      ss2  = surface_list%flux_surfaces(1)%s(3,i_piece)
      dss2 = surface_list%flux_surfaces(1)%s(4,i_piece)
      
      tt1  = surface_list%flux_surfaces(1)%t(1,i_piece)
      dtt1 = surface_list%flux_surfaces(1)%t(2,i_piece)
      tt2  = surface_list%flux_surfaces(1)%t(3,i_piece)
      dtt2 = surface_list%flux_surfaces(1)%t(4,i_piece)
      
      ! --- Loop over nplot points in a flux surface segment
      do ip = 1, nplot
        
        u = -1. + 2.*float(ip-1)/float(nplot-1)
        
        ! --- Determine s and t values of the current point inside element i_elm
        call CUB1D(ss1, dss1, ss2, dss2, u, si, dsi)
        call CUB1D(tt1, dtt1, tt2, dtt2, u, ti, dti)
        
        ! --- Determine (R,Z)-coordinates of the current point on the current flux surface
        call interp_RZ(node_list, element_list, i_elm, si, ti, R, R_s, R_t, R_st, R_ss, R_tt, &
          Z, Z_s, Z_t, Z_st, Z_ss, Z_tt)
        
        if ( no_private .and. ( Z < Z_xpoint(1) ) ) cycle
        
        ! --- Loop over toroidal positions
        do i_phi = 1, n_phi
          
          phi = phi_start + (phi_stop-phi_start)*real(i_phi-1)/real(n_phi-1)
          
          x_cart(ip,i_phi) = + R * cos(phi)
          y_cart(ip,i_phi) = - R * sin(phi)
          z_cart(ip,i_phi) = Z
          value (ip,i_phi) = variable_value(i_var, R, Z, phi, no_zero, error)
          
          if ( (ip == 1) .or. (i_phi == 1) ) cycle
          
          call write_triangle( FH,                                                                 &
            x_cart(ip-1,i_phi-0), y_cart(ip-1,i_phi-0), z_cart(ip-1,i_phi-0),                      &
            x_cart(ip-0,i_phi-1), y_cart(ip-0,i_phi-1), z_cart(ip-0,i_phi-1),                      &
            x_cart(ip-1,i_phi-1), y_cart(ip-1,i_phi-1), z_cart(ip-1,i_phi-1),                      &
            ( value(ip-1,i_phi-0)+value(ip-0,i_phi-1)+value(ip-1,i_phi-1) ) / 3.d0,                &
            'finish{phong 0.4}' )
          
          call write_triangle( FH,                                                                 &
            x_cart(ip-1,i_phi-0), y_cart(ip-1,i_phi-0), z_cart(ip-1,i_phi-0),                      &
            x_cart(ip-0,i_phi-1), y_cart(ip-0,i_phi-1), z_cart(ip-0,i_phi-1),                      &
            x_cart(ip-0,i_phi-0), y_cart(ip-0,i_phi-0), z_cart(ip-0,i_phi-0),                      &
            ( value(ip-1,i_phi-0)+value(ip-0,i_phi-1)+value(ip-0,i_phi-0) ) / 3.d0,                &
            'finish{phong 0.4}' )
          
        end do
        
      end do
      
    end do
    
    ! --- Plot cover surfaces at phi_start and phi_stop
    if ( cover ) then
      
      do i_elm = 1, element_list%n_elements
        if ( mod(i_elm, element_list%n_elements/40+1) == 0 ) then
          write(*,*) 'Element', i_elm, ' of', element_list%n_elements
        end if
        
        do i_phi = 1, 2
          
          phi = phi_start + real(i_phi-1)*(phi_stop-phi_start)
          
          do ip = 1, nplot
            
            si = real(ip-1)/real(nplot-1)
            
            do jp = 1, nplot
              
              ti = real(jp-1)/real(nplot-1)
              
              call interp_RZ(node_list, element_list, i_elm, si, ti, R, R_s, R_t, R_st, R_ss, R_tt,&
                Z, Z_s, Z_t, Z_st, Z_ss, Z_tt)
            
              x_cart(ip,jp) = + R * cos(phi)
              y_cart(ip,jp) = - R * sin(phi)
              z_cart(ip,jp) = Z
              value (ip,jp) = variable_value(i_var, i_elm, si, ti, phi, no_zero)
              
              call interp (node_list, element_list, i_elm, 1, 1, si, ti, P, P_s, P_t, P_st, P_ss,  &
                P_tt)
              
              if ( (ip == 1) .or. (jp == 1) .or. ((P-psi_axis)/(psi_bnd-psi_axis) > psi_surf)      &
                .or. ( no_private .and. (Z < Z_xpoint(1)) ) ) cycle
              
              call write_triangle( FH,                                                             &
                x_cart(ip-0,jp-0), y_cart(ip-0,jp-0), z_cart(ip-0,jp-0),                           &
                x_cart(ip-1,jp-1), y_cart(ip-1,jp-1), z_cart(ip-1,jp-1),                           &
                x_cart(ip-1,jp-0), y_cart(ip-1,jp-0), z_cart(ip-1,jp-0),                           &
                ( value(ip-0,jp-0)+value(ip-1,jp-1)+value(ip-1,jp-0) ) / 3.d0,                     &
                'finish{phong 0.4}' )
  
              call write_triangle( FH,                                                             &
                x_cart(ip-0,jp-0), y_cart(ip-0,jp-0), z_cart(ip-0,jp-0),                           &
                x_cart(ip-1,jp-1), y_cart(ip-1,jp-1), z_cart(ip-1,jp-1),                           &
                x_cart(ip-0,jp-1), y_cart(ip-0,jp-1), z_cart(ip-0,jp-1),                           &
                ( value(ip-0,jp-0)+value(ip-1,jp-1)+value(ip-0,jp-1) ) / 3.d0,                     &
                'finish{phong 0.4}' )
              
            end do
            
          end do
          
        end do
        
      end do
      
    end if
    
  end subroutine plot_surface
  
  
  
  !> Write out a single triangle in povray format
  subroutine write_triangle(FH,x1,y1,z1,x2,y2,z2,x3,y3,z3,v,extra)
    integer, intent(in) :: FH
    real*8,  intent(in) :: x1,y1,z1,x2,y2,z2,x3,y3,z3,v
    character(len=*), intent(in), optional :: extra
    
    if ( ( (x1==x2) .and. (y1==y2) .and. (z1==z2) ) .or. ( (x1==x3) .and. (y1==y3) .and. (z1==z3) )&
      .or. ( (x2==x3) .and. (y2==y3) .and. (z2==z3) ) ) return ! Skip degenerate triangle
    
    150 format(a)
    151 format('  <',es23.16,',',es23.16,',',es23.16,'>')
    152 format(3a)
    write(FH,150) 'triangle{'
    write(FH,151) x1,y1,z1
    write(FH,151) x2,y2,z2
    write(FH,151) x3,y3,z3
    write(FH,152) '  pigment{color ',trim(color_for_value(v)),'}'
    if ( present(extra) ) then
      write(FH,152) '  ',trim(extra)
    end if
    write(FH,150) '}'
    
  end subroutine write_triangle
  
  
  
  !> Determine the color value for a given variable value (depending on the color table).
  character(len=80) function color_for_value(value)
    real*8, intent(in) :: value
    
    real*8 :: rgbt(4), xn
    
    xn = normalized_value(value,min_val,max_val)
    
    select case(i_ct)
    case(1) ! blue - black - red
      rgbt(:) = (/ max(0.d0, 1.d0-2.d0*xn), 0.d0, max(0.d0, -1.d0+2.d0*xn), 0.d0 /)
    case(2) ! blue - white - red
      rgbt(:) = (/ min(1.d0, 2.d0-2.d0*xn), 1.d0-2.d0*abs(xn-0.5d0), min(1.d0, 2.d0*xn), 0.d0 /)
    case(3) ! blue - purple - red
      rgbt(:) = (/ 1.d0-xn, 0.d0, xn, 0.d0 /)
    case(4) ! red - yellow
      rgbt(:) = (/ 1.d0, xn, 0.d0, 0.d0 /)
    case(5) ! yellow - red
      rgbt(:) = (/ 1.d0, 1.d0-xn, 0.d0, 0.d0 /)
    case(6) ! blue - transparent-white - red
      rgbt(:) = (/ min(1.d0, 2.d0-2.d0*xn), 1.d0-2.d0*abs(xn-0.5d0), min(1.d0, 2.d0*xn),           &
        sqrt(1.d0-2.d0*abs(xn-0.5d0)) /)
    end select
    
    160 format('rgbt<',f7.4,',',f7.4,',',f7.4,',',f7.4,'>')
    write(color_for_value,160) rgbt(:)
    
  end function color_for_value
  
  
  
  !> Normalize a variable value
  real*8 function normalized_value(value, minimum, maximum)
    real*8, intent(in) :: value, minimum, maximum
    
    normalized_value = (value-minimum)/(maximum-minimum)
    
    write(80,*) value
    write(81,*) normalized_value
    
    normalized_value = min( 1.d0, max( 0.d0, normalized_value ) )
    
  end function normalized_value
  
  
  
end program jorek2_povray
