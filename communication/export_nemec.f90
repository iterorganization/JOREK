!> Export the equilibrium for the NEMEC code.
!!
!! - iota-for-nemec contains the inverse q-profile versus the normalized toroidal flux
!! - pres-for-nemec contains the pressure profile versus the normalized toroidal flux 
!!   (flux-surface average of the instantaneous pressure)
!! - namelist-for-nemec contains some input parameters required for the NEMEC run
!! - lcms-for-descur contains a representation of the last closed magnetic surface which needs to
!!   be Fourier-expanded by the DESCUR code for NEMEC
subroutine export_nemec(node_list, element_list, xpoint, xcase)
  
  use constants
  use data_structure
  use mod_interp, only: interp, interp_RZ
  use equil_info
  use phys_module, only: with_TiTe, var_T, var_Ti, var_Te, var_rho
  
  implicit none
  
  ! --- Routine parameters
  type (type_node_list),        intent(in) :: node_list
  type (type_element_list),     intent(in) :: element_list
  logical,                      intent(in) :: xpoint
  integer,                      intent(in) :: xcase
  
  ! --- Local variables
  real*8 :: pressure, Phi_edge, sum_dl, dl, rho, temp, temp_i, dum(5), R_s, R_t, Z_s, Z_t
  integer :: ig
  real*8  :: xgs(4), wgs(4)     ! Gauss points between (-1,1), as in export_helena
  data xgs /-0.861136311594053, -0.339981043584856, 0.339981043584856,  0.861136311594053 /
  data wgs / 0.347854845137454,  0.652145154862546, 0.652145154862546,  0.347854845137454 /
  type (type_surface_list) :: surface_list
  integer :: i, i_elm_xpoint(2), i_elm_axis, ifail
  real*8  :: R_xpoint(2), Z_xpoint(2), s_xpoint(2), t_xpoint(2), psi_xpoint(2), psi_bnd
  real*8  :: psi_axis, R_axis, Z_axis, s_axis, t_axis
  real*8, allocatable :: q(:), rad(:), PhiN(:)
  integer :: nplot, j, k, i_elm, node1, node2, node3, node4, ip
  real*8  :: rr1, drr1, rr2, drr2, ss1, dss1, ss2, dss2, t, ri, dri, si, dsi
  real*8  :: R,Z
  
  write(*,*) 'BEGIN: export_nemec'
  
  ! --- Find axis and x-point
  call find_axis(0,node_list,element_list,psi_axis,R_axis,Z_axis,i_elm_axis,s_axis,t_axis,ifail)
  psi_bnd = 0.d0
  if (xpoint) then
    call find_xpoint(0,node_list,element_list,psi_xpoint,R_xpoint,Z_xpoint,i_elm_xpoint,s_xpoint,  &
      t_xpoint,xcase,ifail)
    if (ifail .ne. 1) then
      psi_bnd  = psi_xpoint(1)
      if( ES%active_xpoint .eq. UPPER_XPOINT ) then
        psi_bnd = psi_xpoint(2)
      endif
    endif
  endif
  
  ! --- Find flux surfaces
  surface_list%n_psi = 400
  call tr_allocate(surface_list%psi_values,1,surface_list%n_psi,"surface_list%psi_values",CAT_GRID)
  do i = 1, surface_list%n_psi
    surface_list%psi_values(i) = (float(i)/float(surface_list%n_psi))**2 * (psi_bnd - psi_axis)    &
      + psi_axis
  enddo
  call find_flux_surfaces(0,xpoint,xcase,node_list,element_list,surface_list)
  
  ! --- Determine the q-profile.
  call tr_allocate(q,1,surface_list%n_psi,"q",CAT_GRID)
  call tr_allocate(rad,1,surface_list%n_psi,"rad",CAT_GRID)
  call determine_q_profile(node_list,element_list,surface_list,psi_axis,psi_xpoint,Z_xpoint,q,rad)
  
  ! --- Determine the toroidal flux.
  call tr_allocate(PhiN,1,surface_list%n_psi,"PhiN",CAT_GRID)
  call determine_PhiN(surface_list, q, PhiN, Phi_edge)
  
  ! --- Write out the iota profile for NEMEC.
  write(*,*) 'Writing iota-for-nemec file.'
  open(42, file='iota-for-nemec', action='write', status='replace', form='formatted')
  write(42,*) 'iota profile written by JOREK routine export_nemec'
  write(42,*) 'diota1   diota2   iotaspli'
  write(42,*) '0.0      0.0        0'
  write(42,*) 'niotau'
  write(42,*) surface_list%n_psi
  do i = 1, surface_list%n_psi
    write(42,*) PhiN(i), 1.d0/q(i)
  end do
  close(42)
  
  ! --- Write out the pressure profile for NEMEC.
  write(*,*) 'Writing pres-for-nemec file.'
  open(42, file='pres-for-nemec', action='write', status='replace', form='formatted')
  write(42,*) 'pressure profile written by JOREK routine export_nemec'
  write(42,*) 'dpres1   dpres2   presspli'
  write(42,*) '0.0      0.0        0'
  write(42,*) 'npresu'
  write(42,*) surface_list%n_psi
  do i = 1, surface_list%n_psi
    ! --- Flux-surface average of the actual pressure rho*T (n=0 harmonic), as in export_helena
    pressure = 0.d0
    sum_dl   = 0.d0
    do k = 1, surface_list%flux_surfaces(i)%n_pieces
      i_elm = surface_list%flux_surfaces(i)%elm(k)
      do ig = 1, 4
        call CUB1D(surface_list%flux_surfaces(i)%s(1,k), surface_list%flux_surfaces(i)%s(2,k),         &
          surface_list%flux_surfaces(i)%s(3,k), surface_list%flux_surfaces(i)%s(4,k), xgs(ig), ri, dri)
        call CUB1D(surface_list%flux_surfaces(i)%t(1,k), surface_list%flux_surfaces(i)%t(2,k),         &
          surface_list%flux_surfaces(i)%t(3,k), surface_list%flux_surfaces(i)%t(4,k), xgs(ig), si, dsi)
        call interp_RZ(node_list, element_list, i_elm, ri, si, R, R_s, R_t, Z, Z_s, Z_t)
        ! --- Skip segments in the private flux region (beyond the x-point)
        if ( xpoint ) then
          if ( (xcase /= UPPER_XPOINT) .and. (Z < Z_xpoint(1)) ) cycle
          if ( (xcase /= LOWER_XPOINT) .and. (Z > Z_xpoint(2)) ) cycle
        end if
        rho = 1.d0                                         ! models without density variable
        if ( var_rho > 0 ) call interp(node_list, element_list, i_elm, var_rho, 1, ri, si, rho, dum(1),  &
          dum(2), dum(3), dum(4), dum(5))
        if ( with_TiTe ) then
          call interp(node_list, element_list, i_elm, var_Ti, 1, ri, si, temp_i, dum(1), dum(2), dum(3), dum(4), dum(5))
          call interp(node_list, element_list, i_elm, var_Te, 1, ri, si, temp, dum(1), dum(2), dum(3), dum(4), dum(5))
          temp = temp + temp_i
        else
          call interp(node_list, element_list, i_elm, var_T, 1, ri, si, temp, dum(1), dum(2), dum(3), dum(4), dum(5))
        end if
        dl       = wgs(ig) * sqrt( (R_s*dri + R_t*dsi)**2 + (Z_s*dri + Z_t*dsi)**2 )
        pressure = pressure + rho * temp * dl
        sum_dl   = sum_dl + dl
      end do
    end do
    pressure = pressure / max(sum_dl, 1.d-99) / MU_ZERO
    write(42,*) PhiN(i), pressure
  end do
  close(42)
  
  ! --- Write out some namelist fragments for NEMEC.
  write(*,*) 'Writing namelist-for-nemec file.'
  open(42, file='namelist-for-nemec', action='write', status='replace', form='formatted')
  111 format(1x,a,es14.7,a)
  write(42,*) '! *** written by JOREK routine export_nemec ***'
  write(42,*) 'pres_profile  = "user",'
  write(42,*) 'user_pressure = "./pres-for-nemec",'
  write(42,*) 'iota_profile  = "user",'
  write(42,*) 'user_iota     = "./iota-for-nemec",'
  write(42,111) 'phiedge       = ', -2.d0 * pi * Phi_edge, ','
  write(42,111) 'betascale     = ', 1.d0, ','
  write(42,111) 'raxis_co(0)   = ', R_axis, ','
  write(42,111) 'zaxis_co(0)   = ', Z_axis, ','
  close(42)
  
  ! --- Write out the lcms points for DESCUR.
  j = surface_list%n_psi
  nplot = 6
  write(*,*) 'Writing lcms-for-descur file.'
  open(42, file='lcms-for-descur', action='write', status='replace', form='formatted')
  112 format(1x,i5,a)
  write(42,*) 'ntheta  nphi   nfp  isymm  isort     (written by JOREK routine export_nemec)'
  write(42,112) surface_list%flux_surfaces(j)%n_pieces*(nplot-1), '     1     1      1      1'
  write(42,*) 'plasma boundary'
  write(42,*) ' rin            zin'
  do k = 1, surface_list%flux_surfaces(j)%n_pieces
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
    
    do ip = 1, nplot-1
      
      t = -1. + 2.*float(ip-1)/float(nplot-1)
      
      call CUB1D(rr1, drr1, rr2, drr2, t, ri, dri)
      call CUB1D(ss1, dss1, ss2, dss2, t, si, dsi)
      
      call interp_RZ(node_list,element_list,i_elm,ri,si,R,Z)
      
      write(42,'(2es15.7)') R, Z
      
    end do
    
  end do
  close(42)
  
  call tr_deallocate(q, "q",CAT_GRID)
  call tr_deallocate(rad, "rad",CAT_GRID)
  call tr_deallocate(PhiN, "PhiN",CAT_GRID)
  call tr_deallocate(surface_list%psi_values, "surface_list%psi_values",CAT_GRID)

  write(*,*) 'END: export_nemec'
  
end subroutine export_nemec
