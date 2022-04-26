module mod_bootstrap_functions

  implicit none
  integer, parameter :: n_spline = 30
  real*8             :: q_spline(n_spline), ft_spline(n_spline), B_spline(n_spline)
  real*8             :: psi_knots(n_spline), q_knots(n_spline), ft_knots(n_spline), B_knots(n_spline)
  
  integer, parameter :: n_spline_vtk = 100
  real*8             :: psi_knots_vtk(n_spline_vtk), j_knots_vtk(n_spline_vtk), j_spline_vtk(n_spline_vtk)

contains






!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
!------- THE MAIN FORMULA FROM SAUTER PoP-1999 WITH CORRECTION FROM ERRATA -------------------------
!------- http://scitation.aip.org/content/aip/journal/pop/9/12/10.1063/1.1517052 -------------------
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
subroutine bootstrap_current(R, Z,                           &
                             R_axis, Z_axis, psi_axis,       &
			     R_xpoint, Z_Xpoint, psi_bnd,    &
                             psi_norm, ps0, ps0_x, ps0_y,    &
			     r0,  r0_x,  r0_y,               &
			     Ti0, Ti0_x, Ti0_y,              &
                             Te0, Te0_x, Te0_y,            Jb)
!DEC$ ATTRIBUTES FORCEINLINE :: bootstrap_current
!---------------------------------
! calculates the bootstrap current
!---------------------------------

  use constants
  use phys_module
  use corr_neg

  implicit none
  ! --- Routine parameters
  real*8, intent(in)  :: R, Z
  real*8, intent(in)  :: R_axis, Z_axis, psi_axis
  real*8, intent(in)  :: R_xpoint(2), Z_xpoint(2), psi_bnd
  real*8, intent(in)  :: psi_norm
  real*8, intent(in)  :: ps0, ps0_x, ps0_y
  real*8, intent(in)  :: r0,  r0_x,  r0_y
  real*8              :: rho, rho_x, rho_y, drho_dpsi
  real*8, intent(in)  :: Ti0, Ti0_x, Ti0_y
  real*8              :: Ti,  Ti_x,  Ti_y,  dTi_dpsi, Ti_eV
  real*8              :: P_i, Pi_x,  Pi_y,  dPi_dpsi
  real*8, intent(in)  :: Te0, Te0_x, Te0_y
  real*8              :: Te,  Te_x,  Te_y,  dTe_dpsi, Te_eV
  real*8              :: T,   T_x,   T_y,   dT_dpsi
  real*8              :: Pe,  Pe_x,  Pe_y,  dPe_dpsi
  real*8              :: P,   P_x,   P_y,   dP_dpsi
  real*8              :: grad_psi, psi_n, X, B_tot, B_phi, B_average
  real*8              :: ZZ, Rpe, lnAi, lnAe, II
  real*8              :: F32_ee, F32_ei
  real*8              :: Nue, Nui
  real*8              :: q, ft, E_par, alpha, alpha0
  real*8              :: eps
  real*8              :: L31, L32, L34
  real*8, intent(out) :: Jb
  real*8              :: rho_norm
  real*8              :: tanh_boot, position, width
  real*8              :: distance
  real*8              :: distance_xpoint
  real*8              :: distance_xpoint_axis

  ! --- Need central_density
  if (central_density .lt. 1.d-6) then
    write(*,*)'**********************!!WARNING!!*************************'
    write(*,*)'Element_matrix asks for the bootstrap current,'
    write(*,*)'but the central_density is not defined in the input file.'
    write(*,*)'Returning zero bootstrap current...'
    write(*,*)'**********************!!WARNING!!*************************'
    Jb = 0.d0
    return
  endif
  
  ! --- Careful with convention of density (some people get scared when they see an exponent of 19-20 and prefer to just ignore it...)
  if (central_density .lt. 1.d17) then
    rho_norm = central_density*1.d20 ! (this is so clever... I hope they do it like this at NASA...)
  else
    rho_norm = central_density
  endif
  
  ! --- Renormalise temperature and density 
  ! --- Note for us density is n, not mi*n,
  ! --- but that's ok since formula is with p = n*T
  ! --- ie. rho*T for us...
  ! --- Temperature in Joules
  ! --- Density in 1/(cubic meters)
  rho   = corr_neg_dens(r0)
  rho   = rho   * rho_norm
  rho_x = r0_x  * rho_norm
  rho_y = r0_y  * rho_norm
  Ti    = corr_neg_temp(Ti0)
  Ti    = Ti    / (MU_ZERO*rho_norm)
  Ti_eV = Ti    / 1.6021765d-19
  Ti_x  = Ti0_x / (MU_ZERO*rho_norm)
  Ti_y  = Ti0_y / (MU_ZERO*rho_norm)
  Te    = corr_neg_temp(Te0)
  Te    = Te    / (MU_ZERO*rho_norm)
  Te_eV = Te    / 1.6021765d-19
  Te_x  = Te0_x / (MU_ZERO*rho_norm)
  Te_y  = Te0_y / (MU_ZERO*rho_norm)
  T     = Ti + Te
  T_x   = Ti_x + Te_x
  T_y   = Ti_y + Te_y
  P_i   = rho * Ti
  Pi_x  = rho * Ti_x + rho_x * Ti
  Pi_y  = rho * Ti_y + rho_y * Ti
  Pe    = rho * Te
  Pe_x  = rho * Te_x + rho_x * Te
  Pe_y  = rho * Te_y + rho_y * Te
  P     = rho * T
  P_x   = rho * T_x + rho_x * T
  P_y   = rho * T_y + rho_y * T
  
        
  ! --- Psi variables, including r~a*sqrt(psi) and X=sqrt(2*r/R0)
  psi_n    = psi_norm
  if (psi_n .lt. 1.d-10)  psi_n = 1.d-10 ! careful at axis since X=sqrt(psi_n)
  grad_psi = (ps0_x*ps0_x + ps0_y*ps0_y)**0.5d0
  if (grad_psi .lt. 1.d-10) grad_psi = 1.d-10

  ! --- Derivatives with respect to psi
  drho_dpsi = (rho_x*ps0_x + rho_y*ps0_y) / grad_psi**2.d0
  dTi_dpsi  = (Ti_x *ps0_x + Ti_y *ps0_y) / grad_psi**2.d0
  dTe_dpsi  = (Te_x *ps0_x + Te_y *ps0_y) / grad_psi**2.d0
  dT_dpsi   = dTi_dpsi + dTe_dpsi
  dPe_dpsi  = rho*dTe_dpsi + drho_dpsi*Te
  dPi_dpsi  = rho*dTi_dpsi + drho_dpsi*Ti
  dP_dpsi   = rho*dT_dpsi  + drho_dpsi*T
        
  ! --- Inverse aspect ratio
  eps = minRad / R_axis
        
  ! --- Ion Charge
  ZZ = 1
  
  ! --- Ratio of electron pressure Rpe
  Rpe = Pe / P
  
  ! --- Get q-profile, fraction of trapped particles ft, and averaged B_tot (all together is faster)
  call bootstrap_spline3_eval_all(psi_n,q, ft, B_average)
  
  ! --- lnAi and lnAe
  lnAi = 30   - log(ZZ**3 * sqrt(rho) / Ti_eV**1.5)
  lnAe = 31.3 - log(        sqrt(rho) / Te_eV     )
  
  ! --- Collisionality formula from Sauter paper
  Nui = 4.900d-18 * abs(q) * R_axis * rho * ZZ**4 * abs(lnAi) / (Ti_eV**2 * eps**1.5)
  Nue = 6.921d-18 * abs(q) * R_axis * rho * ZZ    * abs(lnAe) / (Te_eV**2 * eps**1.5)

  ! --- II = I(psi) = RB_phi = F0
  II = F0

  ! --- L31
  X   = ft / ( 1. + (1.-0.1*ft)*sqrt(Nue) + 0.5*(1.-ft)*Nue/ZZ )
  L31 = (1. + 1.4/(ZZ+1.))*X - 1.9/(ZZ+1.)*X**2 + 0.3/(ZZ+1.)*X**3 + 0.2/(ZZ+1.)*X**4

  ! --- L34
  X   = ft / ( 1. + (1.-0.1*ft)*sqrt(Nue) + 0.5*(1.-0.5*ft)*Nue/ZZ )
  L34 = (1. + 1.4/(ZZ+1.))*X - 1.9/(ZZ+1.)*X**2 + 0.3/(ZZ+1.)*X**3 + 0.2/(ZZ+1.)*X**4

  ! --- L32
  X      = ft / ( 1. + 0.26*(1.-ft)*sqrt(Nue) + 0.18*(1.-0.37*ft)*Nue/sqrt(ZZ) )
  F32_ee = (0.05+0.62*ZZ)/(ZZ*(1.+0.44*ZZ))*(X-X**4) + 1./(1.+0.22*ZZ)*(X**2 - X**4 - 1.2*(X**3 - X**4)) + 1.2/(1.+0.5*ZZ)*X**4
  X      = ft / ( 1. + (1.+0.6*ft)*sqrt(Nue) + 0.85*(1.-0.37*ft)*Nue*(1.+ZZ) )
  F32_ei = -(0.56+1.93*ZZ)/(ZZ*(1.+0.44*ZZ))*(X-X**4) + 4.95/(1.+2.48*ZZ)*(X**2 - X**4 - 0.55*(X**3 - X**4)) - 1.2/(1.+0.5*ZZ)*X**4
  L32    = F32_ee + F32_ei

  ! --- alpha
  alpha0 = -1.17*(1.-ft) / (1. - 0.22*ft - 0.19*ft**2)
  alpha  = ( (alpha0 + 0.25*(1.-ft**2)*sqrt(Nui)) / (1.0+0.5*sqrt(Nui)) + 0.315*Nui**2*ft**6 ) / (1.0 + 0.15*Nui**2*ft**6)
  
  ! --- Bootstrap Current (this is the formula to be compared with HELENA, which gives <j*B>)
  ! --- The Sauter estimates: L31 = 0.5 ; L34 = 0.5 ; L32 = -0.2 ; alpha = -0.5
  Jb = - II * Pe * ( L31*dP_dpsi/Pe + L32*dTe_dpsi/Te + L34*alpha*(1-Rpe)/Rpe*dTi_dpsi/Ti )

  ! --- Current with denormalisation
  ! --- (one -R because the JOREK current is defined so, mu0 of course,
  ! --- and then B_average because Sauter's formula gives <jb*B>)
  Jb = - R * Jb * MU_ZERO / B_average
  
  ! --- Project along phi (formula is for j_par, and we want j_phi)
  B_tot =  sqrt( (F0/R)**2 + (grad_psi/R)**2 )
  B_phi =  sqrt( (F0/R)**2 )
  Jb = Jb * B_phi / B_tot
  Jb = abs(Jb)
  
  ! --- There should not be any bootstrap outside plasma, the Xpoint can be noisy...
  Jb = Jb * (0.5d0 - 0.5d0 * tanh( (psi_norm - 1.01)/0.005d0 ) )
  ! --- Cut off bootstrap source around the Xpoint with a radius of 5% the distance Xpoint-axis.
  if (xpoint .and.  (xcase .ne. UPPER_XPOINT) ) then
    distance_xpoint      = sqrt( (R      - R_xpoint(1))**2 + (Z      - Z_xpoint(1))**2 )
    distance_xpoint_axis = sqrt( (R_axis - R_xpoint(1))**2 + (Z_axis - Z_xpoint(1))**2 )
    distance = 0.05 * distance_xpoint_axis
    Jb = Jb * (0.5d0 - 0.5d0 * tanh( -(distance_xpoint - distance)/0.01d0 ) )
  endif
  if (xpoint .and.  (xcase .ne. LOWER_XPOINT) ) then
    distance_xpoint      = sqrt( (R      - R_xpoint(2))**2 + (Z      - Z_xpoint(2))**2 )
    distance_xpoint_axis = sqrt( (R_axis - R_xpoint(2))**2 + (Z_axis - Z_xpoint(2))**2 )
    distance = 0.05 * distance_xpoint_axis
    Jb = Jb * (0.5d0 - 0.5d0 * tanh( -(distance_xpoint - distance)/0.01d0 ) )
  endif
  
return
end subroutine bootstrap_current















!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
!------------------------------- FIND THE MINOR RADIUS ---------------------------------------------
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
subroutine bootstrap_find_minRad(node_list, element_list, R_axis, Z_axis, psi_axis, psi_bnd)

  use data_structure
  use phys_module
  use mod_interp

  implicit none
  ! --- Routine parameters
  type (type_node_list),        intent(inout) :: node_list
  type (type_element_list),     intent(inout) :: element_list
  real*8, 			intent(in)    :: R_axis, Z_axis
  real*8, 			intent(in)    :: psi_axis, psi_bnd
  
  ! --- Internal parameters
  type (type_surface_list) 	:: surface_list, flux_list
  integer			:: n_iter, n_iter_max
  real*8			:: step
  real*8			:: R_find, Z_find
  real*8			:: R_out,  Z_out
  real*8			:: s_out,  t_out
  integer			:: i_elm_out, ifail
  real*8			:: s_find(8), t_find(8)
  integer			:: i_elm_find(8),i_find
  real*8			:: psi, psi_norm, psi_s,psi_t,psi_st,psi_ss,psi_tt
  logical			:: found

  ! --- Simplest case when we have a limiter plasma
  if (.not. xpoint) then
    flux_list%n_psi = 1
    call tr_allocate(flux_list%psi_values,1,flux_list%n_psi,"flux_list%psi_values",CAT_GRID)
    flux_list%psi_values(1) = psi_bnd
    call find_flux_surfaces(0, xpoint,xcase,node_list,element_list,flux_list)
    call find_theta_surface(node_list, element_list, flux_list, 1, 0.0, R_axis, Z_axis,i_elm_find,s_find,t_find,i_find)
    ! --- If this didn't work, it means psi=1.0 is the grid boundary, try with psi=0.99
    if (i_find .eq. 0) then
      flux_list%psi_values(1) = 0.99 * (psi_bnd - psi_axis) + psi_axis
      call find_flux_surfaces(0,xpoint,xcase,node_list,element_list,flux_list)
      call find_theta_surface(node_list, element_list, flux_list, 1, 0.0, R_axis, Z_axis,i_elm_find,s_find,t_find,i_find)
    endif
    if (i_find .ne. 0) then
      call interp_RZ(node_list,element_list,i_elm_find(1),s_find(1),t_find(1),R_find,Z_find)
      minRad = R_find - R_axis
    else
      minRad = amin
    endif
    call tr_deallocate(flux_list%psi_values,"flux_list%psi_values",CAT_GRID)
    return
  endif
  
  ! --- Step along line with 2cm resolution
  n_iter     = 0
  step       = 0.02d0
  n_iter_max = 500 ! 10m should be largely sufficient for any machine
  found      = .false.
  R_find     = R_axis
  Z_find     = Z_axis
  do while ( (n_iter .lt. n_iter_max) .and. (.not. found) )

    n_iter = n_iter + 1
    
    R_find = R_find + step
    call find_RZ(node_list,element_list, R_find,Z_find, R_out,Z_out, i_elm_out, s_out,t_out,ifail)
    if (ifail .ne. 0) then
      found = .false.
      exit
    else
      call interp(node_list,element_list,i_elm_out,1,1,s_out,t_out, psi, psi_s,psi_t,psi_st,psi_ss,psi_tt)
      psi_norm = (psi-psi_axis) / (psi_bnd-psi_axis)
      if (psi_norm .gt. 1.d0) then
        found = .true.
	exit
      endif
    endif
  
  enddo
  
  ! --- Step along line with 1mm resolution from previous location
  n_iter = 0
  step = 1.d-3
  if (found) then
    R_find = R_find - 0.02d0 ! step back 2cm
    n_iter_max = 30         ! so 3cm should be sufficient
  else
    R_find = R_axis
    n_iter_max = 5000
  endif
  found  = .false.
  do while ( (n_iter .lt. n_iter_max) .and. (.not. found) )

    n_iter = n_iter + 1
    
    R_find = R_find + step
    call find_RZ(node_list,element_list, R_find,Z_find, R_out,Z_out, i_elm_out, s_out,t_out,ifail)
    if (ifail .ne. 0) then
      found = .false.
      exit
    else
      call interp(node_list,element_list,i_elm_out,1,1,s_out,t_out, psi, psi_s,psi_t,psi_st,psi_ss,psi_tt)
      psi_norm = (psi-psi_axis) / (psi_bnd-psi_axis)
      if (psi_norm .gt. 1.d0) then
        found = .true.
	exit
      endif
    endif
  
  enddo
  
  ! --- If we still haven't found it, try with surfaces
  if (.not. found) then
    flux_list%n_psi = 1
    call tr_allocate(flux_list%psi_values,1,flux_list%n_psi,"flux_list%psi_values",CAT_GRID)
    flux_list%psi_values(1) = psi_bnd
    call find_flux_surfaces(0,xpoint,xcase,node_list,element_list,flux_list)
    call find_theta_surface(node_list, element_list, flux_list, 1, 0.0, R_axis, Z_axis,i_elm_find,s_find,t_find,i_find)
    call interp_RZ(node_list,element_list,i_elm_find(1),s_find(1),t_find(1),R_find,Z_find)
    call tr_deallocate(flux_list%psi_values,"flux_list%psi_values",CAT_GRID)
    minRad = R_find - R_axis
  else
    minRad = R_find - R_axis
  endif

end subroutine bootstrap_find_minRad










!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
!------------------------- GET SPLINES FOR averaged j-profile (or averaged-anything) ---------------
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
subroutine bootstrap_get_averaged_j_spline(node_list, element_list, psi_axis, psi_xpoint, R_xpoint, Z_xpoint)

  use data_structure
  use phys_module
  use grid_xpoint_data
  use mod_interp
  use equil_info

  implicit none
  ! --- Routine parameters
  type (type_node_list),    intent(inout) :: node_list
  type (type_element_list), intent(inout) :: element_list
  real*8,                   intent(in)    :: psi_axis, psi_xpoint(2)
  real*8,                   intent(in)    :: R_xpoint(2), Z_xpoint(2)
  
  ! --- Gaussian points between (-1.,1.) for Gauss-integration
  real*8, parameter :: xgs(4) = (/-0.861136311594053, -0.339981043584856, 0.339981043584856,  0.861136311594053 /)
  real*8, parameter :: wgs(4) = (/ 0.347854845137454,  0.652145154862546, 0.652145154862546,  0.347854845137454 /)

  ! --- Internal parameters
  type (type_surface_list) :: flux_list, sep_list
  integer                  :: i, k, ig, i_surf, i_piece, n_psi, i_elm
  real*8                   :: psi_bnd, psi_bnd2
  real*8                   :: sigmas(16)
  integer                  :: n_grids(12)
  real*8                   :: rr, s, t, ds, dt, xjac, dl, sum_dl
  real*8                   :: R, dR_ds, dR_dt, dR_dl
  real*8                   :: Z, dZ_ds, dZ_dt, dZ_dl
  real*8                   :: psi,dpsi_ds,dpsi_dt,dpsi_dst,dpsi_dss,dpsi_dtt, psi_R, psi_Z
  real*8                   :: zj, dzj_ds, dzj_dt, dzj_dst, dzj_dss, dzj_dtt
  real*8                   :: psi_xpoint_tmp(2)


  ! ------------------------------------------
  ! --- Define flux values for the spline knot
  ! ------------------------------------------
  
  ! --- Reset the grid parameters so we can use the grid_xpoint function directly
  n_psi  = n_flux + n_open + n_outer
  n_flux = n_flux * n_spline_vtk / n_psi
  if (xcase .eq. DOUBLE_NULL) then
    n_open  = n_open  * n_spline_vtk / n_psi
    n_outer = n_outer * n_spline_vtk / n_psi
  else
    n_open  = n_open  * n_spline_vtk / n_psi
    n_outer = 0
  endif
  n_psi = n_flux + n_open + n_outer
  if (n_psi .ne. n_spline_vtk) n_flux = n_flux + (n_spline_vtk - n_psi)
  n_psi = n_flux + n_open + n_outer
  n_inner   = 0
  n_up_priv = 0
  n_up_leg  = 0
  n_private = 0
  n_leg     = 0
  
  ! --- Build up some arrays to send as routine parameters (avoid long lists...)
  ! --- Take 70% of dPSI_open and dPSI_outer to make sure we stay inside domain...
  sigmas(1)  = SIG_closed   ; sigmas(2)  = SIG_theta
  sigmas(3)  = SIG_open     ; sigmas(4)  = SIG_outer     ; sigmas(5)  = SIG_inner
  sigmas(6)  = SIG_private  ; sigmas(7)  = SIG_up_priv
  sigmas(8)  = SIG_leg_0    ; sigmas(9)  = SIG_leg_1
  sigmas(10) = SIG_up_leg_0 ; sigmas(11) = SIG_up_leg_1
  sigmas(12) = dPSI_open*0.7; sigmas(13) = dPSI_outer*0.7; sigmas(14) = dPSI_inner
  sigmas(15) = dPSI_private ; sigmas(16) = dPSI_up_priv

  n_grids(1) = n_flux   ; n_grids(2) = n_tht
  n_grids(3) = n_open   ; n_grids(4) = n_outer  ; n_grids(5) = n_inner
  n_grids(6) = n_private; n_grids(7) = n_up_priv
  n_grids(8) = n_leg    ; n_grids(9) = n_up_leg
  
  ! --- Get psi_bnd
  psi_bnd  = 0.d0
  psi_bnd2 = 0.d0
  psi_xpoint_tmp = psi_xpoint
  if(xcase .eq. LOWER_XPOINT) psi_bnd = psi_xpoint(1)
  if(xcase .eq. UPPER_XPOINT) psi_bnd = psi_xpoint(2)
  if(xcase .eq. DOUBLE_NULL ) then
    if(ES%active_xpoint .eq. UPPER_XPOINT) then
      psi_bnd  = psi_xpoint(2)
      psi_bnd2 = psi_xpoint(1)
    else
      psi_bnd  = psi_xpoint(1)
      psi_bnd2 = psi_xpoint(2)  
    endif
    ! If we have a symmetric double-null, force the single separatrix
    if (abs(psi_xpoint(1)-psi_xpoint(2)) .lt. SDN_threshold) then
      psi_xpoint_tmp(1) = (psi_xpoint(1)+psi_xpoint(2))/2.d0
      psi_xpoint_tmp(2) = psi_xpoint_tmp(1)
      psi_bnd  = psi_xpoint_tmp(1)
      psi_bnd2 = psi_bnd  
      n_grids(3) = 0
    endif
  endif

  ! --- Define number of psi values and allocate flux_list structure
  flux_list%n_psi = n_psi
  call tr_allocate(flux_list%psi_values,1,flux_list%n_psi,"flux_list%psi_values",CAT_GRID)
  
  ! --- Allocate sep_list structure (that's for plotting only)
  sep_list%n_psi =3
  if(xcase .eq. DOUBLE_NULL) sep_list%n_psi =6
  call tr_allocate(sep_list%psi_values,1,sep_list%n_psi,"sep_list%psi_values",CAT_GRID)
  
  ! --- Call the routine
  call define_flux_values(node_list, element_list, flux_list, sep_list, xcase, psi_xpoint_tmp, n_grids, sigmas)
  
  
  ! ---------------------------------------------
  ! --- Compute the fraction of trapped particles
  ! ---------------------------------------------
  
  ! --- Then calculate <j>
  do i=2, flux_list%n_psi
    j_knots_vtk(i) = 0.d0
    sum_dl = 0.d0
    do k=1, flux_list%flux_surfaces(i)%n_pieces
      do ig = 1, 4
  	
	rr = xgs(ig)
  	i_elm = flux_list%flux_surfaces(i)%elm(k)
        call compute_surface_basics(flux_list, i, k, rr, s, t, ds, dt)
  	call interp_RZ(node_list,element_list,i_elm,s,t,R,dR_ds,dR_dt,Z,dZ_ds,dZ_dt)
  	call interp(node_list,element_list,i_elm,3,1,s,t,zj, dzj_ds, dzj_dt, dzj_dst, dzj_dss, dzj_dtt)
  	call interp(node_list,element_list,i_elm,1,1,s,t,psi,dpsi_ds,dpsi_dt,dpsi_dst,dpsi_dss,dpsi_dtt)

  	! --- Ignore flux surface segments in the private flux region below the x-point.
  	if ( xpoint .and. ((psi-psi_axis)/(psi_xpoint_tmp(1)-psi_axis) < 1.d0) .and. (Z < z_xpoint(1)) .and. (xcase .ne. UPPER_XPOINT)) cycle
  	if ( xpoint .and. ((psi-psi_axis)/(psi_xpoint_tmp(2)-psi_axis) < 1.d0) .and. (Z > z_xpoint(2)) .and. (xcase .ne. LOWER_XPOINT)) cycle

        dR_dl = dR_ds * ds + dR_dt * dt
        dZ_dl = dZ_ds * ds + dZ_dt * dt
        dl = sqrt(dR_dl**2 + dZ_dl**2)
  	
        sum_dl = sum_dl +  wgs(ig) * dl

        j_knots_vtk(i) = j_knots_vtk(i) +  wgs(ig) * dl * zj

      enddo
    enddo
    j_knots_vtk(i) = j_knots_vtk(i) / sum_dl
  enddo
  
  
  ! -------------------
  ! --- Fit the splines
  ! -------------------
  
  ! --- Copy flux values to array
  do i=1,n_psi
    psi_knots_vtk(i) = (flux_list%psi_values(i) - psi_axis)/(psi_bnd - psi_axis)
  enddo
  
  ! --- Fit splines
  call bootstrap_spline3_coef(n_spline_vtk-1, psi_knots_vtk, j_knots_vtk, j_spline_vtk)
  
  
    
  
  ! ---------------------
  ! --- Clean up and exit
  ! ---------------------
  
  ! --- Clear surfaces
  call tr_unregister_mem(sizeof(flux_list%flux_surfaces),"flux_list%flux_surfaces")
  deallocate(flux_list%flux_surfaces)
  call tr_unregister_mem(sizeof(flux_list%psi_values),"flux_list%psi_values")
  deallocate(flux_list%psi_values)
  call tr_unregister_mem(sizeof(sep_list%flux_surfaces),"sep_list%flux_surfaces")
  deallocate(sep_list%flux_surfaces)
  call tr_unregister_mem(sizeof(sep_list%psi_values),"sep_list%psi_values")
  deallocate(sep_list%psi_values)
  

end subroutine bootstrap_get_averaged_j_spline








!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
!------------------------- GET SPLINES FOR q-PROFILE AND FT-PROFILE --------------------------------
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
subroutine bootstrap_get_q_and_ft_splines(node_list, element_list, psi_axis, psi_xpoint, R_xpoint, Z_xpoint)

  use data_structure
  use phys_module
  use grid_xpoint_data
  use mod_interp
  use equil_info

  implicit none
  ! --- Routine parameters
  type (type_node_list),    intent(inout) :: node_list
  type (type_element_list), intent(inout) :: element_list
  real*8,                   intent(in)    :: psi_axis, psi_xpoint(2)
  real*8,                   intent(in)    :: R_xpoint(2), Z_xpoint(2)
  
  ! --- Gaussian points between (-1.,1.) for Gauss-integration
  real*8, parameter :: xgs(4) = (/-0.861136311594053, -0.339981043584856, 0.339981043584856,  0.861136311594053 /)
  real*8, parameter :: wgs(4) = (/ 0.347854845137454,  0.652145154862546, 0.652145154862546,  0.347854845137454 /)

  ! --- Internal parameters
  type (type_surface_list) :: flux_list, sep_list
  real*8                   :: rad(n_spline), Bmax(n_spline)
  integer                  :: i, k, ig, i_surf, i_piece, n_psi, i_elm, i_ft
  real*8                   :: psi_bnd, psi_bnd2
  real*8                   :: sigmas(16)
  integer                  :: n_grids(12)
  real*8                   :: rr, s, t, ds, dt, xjac, dl, sum_dl
  real*8                   :: R,  dR_ds,  dR_dt, dR_dl
  real*8                   :: Z,  dZ_ds,  dZ_dt, dZ_dl
  real*8                   :: psi,dpsi_ds,dpsi_dt,dpsi_dst,dpsi_dss,dpsi_dtt, psi_R, psi_Z
  real*8                   :: grad_psi, B_tot, B_pol
  real*8                   :: hh2(n_spline), ft_int(n_spline)
  integer                  :: n_int
  real*8                   :: lambda_ft, dlambda_ft, OneMinusLh
  real*8                   :: psi_n
  real*8                   :: psi_xpoint_tmp(2)


  ! ------------------------------------------
  ! --- Define flux values for the spline knot
  ! ------------------------------------------
  
  ! --- Reset the grid parameters so we can use the grid_xpoint function directly
  n_psi  = n_flux + n_open + n_outer ! includes the axis
  n_flux = n_flux * n_spline / n_psi
  if (xcase .eq. DOUBLE_NULL) then
    n_open  = n_open  * n_spline / n_psi
    n_outer = n_outer * n_spline / n_psi
  else
    n_open  = n_open  * n_spline / n_psi
    n_outer = 0
  endif
  n_psi = n_flux + n_open + n_outer
  if (n_psi .ne. n_spline) n_flux = n_flux + (n_spline - n_psi)
  n_psi = n_flux + n_open + n_outer
  n_inner   = 0
  n_up_priv = 0
  n_up_leg  = 0
  n_private = 0
  n_leg     = 0
  
  ! --- Build up some arrays to send as routine parameters (avoid long lists...)
  ! --- Take 70% of dPSI_open and dPSI_outer to make sure we stay inside domain...
  sigmas(1)  = SIG_closed   ; sigmas(2)  = SIG_theta
  sigmas(3)  = SIG_open     ; sigmas(4)  = SIG_outer     ; sigmas(5)  = SIG_inner
  sigmas(6)  = SIG_private  ; sigmas(7)  = SIG_up_priv
  sigmas(8)  = SIG_leg_0    ; sigmas(9)  = SIG_leg_1
  sigmas(10) = SIG_up_leg_0 ; sigmas(11) = SIG_up_leg_1
  sigmas(12) = dPSI_open*0.7; sigmas(13) = dPSI_outer*0.7; sigmas(14) = dPSI_inner
  sigmas(15) = dPSI_private ; sigmas(16) = dPSI_up_priv

  n_grids(1) = n_flux   ; n_grids(2) = n_tht
  n_grids(3) = n_open   ; n_grids(4) = n_outer  ; n_grids(5) = n_inner
  n_grids(6) = n_private; n_grids(7) = n_up_priv
  n_grids(8) = n_leg    ; n_grids(9) = n_up_leg

  ! --- Get psi_bnd
  psi_bnd  = 0.d0
  psi_bnd2 = 0.d0
  psi_xpoint_tmp = psi_xpoint
  if(xcase .eq. LOWER_XPOINT) psi_bnd = psi_xpoint(1)
  if(xcase .eq. UPPER_XPOINT) psi_bnd = psi_xpoint(2)
  if(xcase .eq. DOUBLE_NULL ) then
    if(ES%active_xpoint .eq. UPPER_XPOINT) then
      psi_bnd  = psi_xpoint(2)
      psi_bnd2 = psi_xpoint(1)
    else
      psi_bnd  = psi_xpoint(1)
      psi_bnd2 = psi_xpoint(2)  
    endif
    ! If we have a symmetric double-null, force the single separatrix
    if (abs(psi_xpoint(1)-psi_xpoint(2)) .lt. SDN_threshold) then
      psi_xpoint_tmp(1) = (psi_xpoint(1)+psi_xpoint(2))/2.d0
      psi_xpoint_tmp(2) = psi_xpoint_tmp(1)
      psi_bnd  = psi_xpoint_tmp(1)
      psi_bnd2 = psi_bnd  
      n_grids(3) = 0
    endif
  endif

  ! --- Define number of psi values and allocate flux_list structure
  flux_list%n_psi = n_psi
  call tr_allocate(flux_list%psi_values,1,flux_list%n_psi,"flux_list%psi_values",CAT_GRID)
  
  ! --- Allocate sep_list structure (that's for plotting only)
  sep_list%n_psi =3
  if(xcase .eq. DOUBLE_NULL) sep_list%n_psi =6
  call tr_allocate(sep_list%psi_values,1,sep_list%n_psi,"sep_list%psi_values",CAT_GRID)
  
  ! --- Call the routine
  call define_flux_values(node_list, element_list, flux_list, sep_list, xcase, psi_xpoint_tmp, n_grids, sigmas)
  
  ! -----------------
  ! --- Get q-profile
  ! -----------------
  
  ! --- Get q-profile on spline knots
  call determine_q_profile(node_list,element_list,flux_list,psi_axis,psi_xpoint_tmp,Z_xpoint,q_knots,rad)

  
  ! ------------------
  ! --- Get B-averaged
  ! ------------------
  
  ! --- Calculate <B>
  do i=1, flux_list%n_psi
    B_knots(i) = 0.d0
    sum_dl = 0.d0
    do k=1, flux_list%flux_surfaces(i)%n_pieces
      do ig = 1, 4
  	
	rr = xgs(ig)
  	i_elm = flux_list%flux_surfaces(i)%elm(k)
        call compute_surface_basics(flux_list, i, k, rr, s, t, ds, dt)
  	call interp_RZ(node_list,element_list,i_elm,s,t,R,dR_ds,dR_dt,Z,dZ_ds,dZ_dt)
  	call interp(node_list,element_list,i_elm,1,1,s,t,psi,dpsi_ds,dpsi_dt,dpsi_dst,dpsi_dss,dpsi_dtt)

  	! --- Ignore flux surface segments in the private flux region below the x-point.
  	if ( xpoint .and. ((psi-psi_axis)/(psi_xpoint_tmp(1)-psi_axis) < 1.d0) .and. (Z < z_xpoint(1)) .and. (xcase .ne. UPPER_XPOINT)) cycle
  	if ( xpoint .and. ((psi-psi_axis)/(psi_xpoint_tmp(2)-psi_axis) < 1.d0) .and. (Z > z_xpoint(2)) .and. (xcase .ne. LOWER_XPOINT)) cycle

  	xjac  = dR_ds * dZ_dt - dR_dt * dZ_ds

	psi_R = (   dpsi_ds * dZ_dt - dpsi_dt * dZ_ds ) / xjac
  	psi_Z = ( - dpsi_ds * dR_dt + dpsi_dt * dR_ds ) / xjac

  	grad_psi = sqrt(psi_R**2 + psi_Z**2)

  	B_tot =  sqrt( (F0/R)**2 + (grad_psi/R)**2 )
  	B_pol =  sqrt(             (grad_psi/R)**2 )
	
        dR_dl = dR_ds * ds + dR_dt * dt
        dZ_dl = dZ_ds * ds + dZ_dt * dt
        dl = sqrt(dR_dl**2 + dZ_dl**2)
  	
        sum_dl = sum_dl +  wgs(ig) * dl / B_pol

        B_knots(i) = B_knots(i) +  wgs(ig) * dl * B_tot / B_pol

      enddo
    enddo
    B_knots(i) = B_knots(i) / sum_dl
  enddo
  

  
  ! ---------------------------------------------
  ! --- Compute the fraction of trapped particles
  ! ---------------------------------------------
  
  ! --- First get Bmax for each surface
  do i=1, flux_list%n_psi
    Bmax(i) = 0.d0
    do k=1, flux_list%flux_surfaces(i)%n_pieces
      do ig = 1, 4
  	
	rr = xgs(ig)
  	i_elm = flux_list%flux_surfaces(i)%elm(k)
        call compute_surface_basics(flux_list, i, k, rr, s, t, ds, dt)
  	call interp_RZ(node_list,element_list,i_elm,s,t,R,dR_ds,dR_dt,Z,dZ_ds,dZ_dt)
  	call interp(node_list,element_list,i_elm,1,1,s,t,psi,dpsi_ds,dpsi_dt,dpsi_dst,dpsi_dss,dpsi_dtt)

  	! --- Ignore flux surface segments in the private flux region below the x-point.
  	if ( xpoint .and. ((psi-psi_axis)/(psi_xpoint_tmp(1)-psi_axis) < 1.d0) .and. (Z < z_xpoint(1)) .and. (xcase .ne. UPPER_XPOINT)) cycle
  	if ( xpoint .and. ((psi-psi_axis)/(psi_xpoint_tmp(2)-psi_axis) < 1.d0) .and. (Z > z_xpoint(2)) .and. (xcase .ne. LOWER_XPOINT)) cycle

  	xjac  = dR_ds * dZ_dt - dR_dt * dZ_ds

  	psi_R = (   dpsi_ds * dZ_dt - dpsi_dt * dZ_ds ) / xjac
  	psi_Z = ( - dpsi_ds * dR_dt + dpsi_dt * dR_ds ) / xjac

  	grad_psi = sqrt(psi_R**2 + psi_Z**2)

  	B_tot =  sqrt( (F0/R)**2 + (grad_psi/R)**2 )
	
	if (B_tot > Bmax(i)) Bmax(i) = B_tot

      enddo
    enddo
  enddo
  
  ! --- Then calculate <h**2> = <(B/Bmax)**2>
  do i=1, flux_list%n_psi
    hh2(i) = 0.d0
    sum_dl = 0.d0
    do k=1, flux_list%flux_surfaces(i)%n_pieces
      do ig = 1, 4
  	
	rr = xgs(ig)
  	i_elm = flux_list%flux_surfaces(i)%elm(k)
        call compute_surface_basics(flux_list, i, k, rr, s, t, ds, dt)
  	call interp_RZ(node_list,element_list,i_elm,s,t,R,dR_ds,dR_dt,Z,dZ_ds,dZ_dt)
  	call interp(node_list,element_list,i_elm,1,1,s,t,psi,dpsi_ds,dpsi_dt,dpsi_dst,dpsi_dss,dpsi_dtt)

  	! --- Ignore flux surface segments in the private flux region below the x-point.
  	if ( xpoint .and. ((psi-psi_axis)/(psi_xpoint_tmp(1)-psi_axis) < 1.d0) .and. (Z < z_xpoint(1)) .and. (xcase .ne. UPPER_XPOINT)) cycle
  	if ( xpoint .and. ((psi-psi_axis)/(psi_xpoint_tmp(2)-psi_axis) < 1.d0) .and. (Z > z_xpoint(2)) .and. (xcase .ne. LOWER_XPOINT)) cycle

  	xjac  = dR_ds * dZ_dt - dR_dt * dZ_ds

	psi_R = (   dpsi_ds * dZ_dt - dpsi_dt * dZ_ds ) / xjac
  	psi_Z = ( - dpsi_ds * dR_dt + dpsi_dt * dR_ds ) / xjac

  	grad_psi = sqrt(psi_R**2 + psi_Z**2)

  	B_tot =  sqrt( (F0/R)**2 + (grad_psi/R)**2 )
  	B_pol =  sqrt(             (grad_psi/R)**2 )
	
        dR_dl = dR_ds * ds + dR_dt * dt
        dZ_dl = dZ_ds * ds + dZ_dt * dt
        dl = sqrt(dR_dl**2 + dZ_dl**2)
  	
        sum_dl = sum_dl +  wgs(ig) * dl / B_pol

        hh2(i) = hh2(i) +  wgs(ig) * dl * (B_tot/Bmax(i))**2.0 / B_pol

      enddo
    enddo
    hh2(i) = hh2(i) / sum_dl
  enddo
  
  ! --- Integrate from 0 to 1 (or from 0 to Bmax, but we normalise here, as in Y. R. Lin-Liu and R. L. Miller, PoP 1995)
  n_int = 100
  dlambda_ft = real(1)/real(n_int-1)
  do i=1, flux_list%n_psi
    ft_int(i) = 0.d0
    do i_ft = 1,n_int-1 ! --- no need to go to the end, the last one is either 0 or NaN, since lambda_ft*B_tot/Bmax(i) ~ 1
      lambda_ft = real(i_ft-1)/real(n_int-1)
      OneMinusLh = 0.d0
      sum_dl = 0.d0
      do k=1, flux_list%flux_surfaces(i)%n_pieces
    	do ig = 1, 4
    	  
  	  rr = xgs(ig)
    	  i_elm = flux_list%flux_surfaces(i)%elm(k)
    	  call compute_surface_basics(flux_list, i, k, rr, s, t, ds, dt)
    	  call interp_RZ(node_list,element_list,i_elm,s,t,R,dR_ds,dR_dt,Z,dZ_ds,dZ_dt)
    	  call interp(node_list,element_list,i_elm,1,1,s,t,psi,dpsi_ds,dpsi_dt,dpsi_dst,dpsi_dss,dpsi_dtt)

    	  ! --- Ignore flux surface segments in the private flux region below the x-point.
    	  if ( xpoint .and. ((psi-psi_axis)/(psi_xpoint_tmp(1)-psi_axis) < 1.d0) .and. (Z < z_xpoint(1)) .and. (xcase .ne. UPPER_XPOINT)) cycle
    	  if ( xpoint .and. ((psi-psi_axis)/(psi_xpoint_tmp(2)-psi_axis) < 1.d0) .and. (Z > z_xpoint(2)) .and. (xcase .ne. LOWER_XPOINT)) cycle

    	  dR_dl = dR_ds * ds + dR_dt * dt
    	  dZ_dl = dZ_ds * ds + dZ_dt * dt

    	  dl = sqrt(dR_dl**2 + dZ_dl**2)

    	  xjac  = dR_ds * dZ_dt - dR_dt * dZ_ds

    	  psi_R = (   dpsi_ds * dZ_dt - dpsi_dt * dZ_ds ) / xjac
    	  psi_Z = ( - dpsi_ds * dR_dt + dpsi_dt * dR_ds ) / xjac

    	  grad_psi = sqrt(psi_R**2 + psi_Z**2)

    	  B_tot =  sqrt( (F0/R)**2 + (grad_psi/R)**2 )
  	  B_pol =  sqrt(             (grad_psi/R)**2 )

          dR_dl = dR_ds * ds + dR_dt * dt
          dZ_dl = dZ_ds * ds + dZ_dt * dt
          dl = sqrt(dR_dl**2 + dZ_dl**2)
  	
          sum_dl = sum_dl +  wgs(ig) * dl / B_pol

          OneMinusLh = OneMinusLh +  wgs(ig) * dl * sqrt(1.0 - lambda_ft*B_tot/Bmax(i)) / B_pol
    	enddo
      enddo
      OneMinusLh = OneMinusLh / sum_dl
      ft_int(i) = ft_int(i) + lambda_ft*dlambda_ft/OneMinusLh
    enddo
  enddo
  
  ! --- The final formula for the trapped fraction
  do i=1, flux_list%n_psi
    ft_knots(i) = 1 - 0.75*hh2(i)*ft_int(i)
    psi_n = (flux_list%psi_values(i) - psi_axis)/(psi_bnd - psi_axis)
    ft_knots(i) = ft_knots(i) * (0.5 - 0.5 * tanh( -(psi_n - 0.02)/0.04) )
  enddo
  
  ! --- Copy flux values to array
  do i=1,n_psi
    psi_knots(i) = (flux_list%psi_values(i) - psi_axis)/(psi_bnd - psi_axis)
  enddo
  
  ! -------------------
  ! --- Fit the splines
  ! -------------------
  
  ! --- First we reorder the flux values to have the first one at axis
  do i=n_psi,2,-1
    psi_knots(i) = psi_knots(i-1)
    q_knots(i)   = q_knots(i-1)
    B_knots(i)   = B_knots(i-1)
    ft_knots(i)  = ft_knots(i-1)
  enddo
  psi_knots(1) = 0.0
  q_knots(1)   = q_knots(2)
  B_knots(1)   = B_knots(2)
  ft_knots(1)  = 0.0
  
  
  ! --- Fit splines
  call bootstrap_spline3_coef(n_spline-1, psi_knots, q_knots,  q_spline)
  call bootstrap_spline3_coef(n_spline-1, psi_knots, B_knots,  B_spline)
  call bootstrap_spline3_coef(n_spline-1, psi_knots, ft_knots, ft_spline)
  
  
    
  
  ! ---------------------
  ! --- Clean up and exit
  ! ---------------------
  
  ! --- Clear surfaces
  call tr_unregister_mem(sizeof(flux_list%flux_surfaces),"flux_list%flux_surfaces")
  deallocate(flux_list%flux_surfaces)
  call tr_unregister_mem(sizeof(flux_list%psi_values),"flux_list%psi_values")
  deallocate(flux_list%psi_values)
  call tr_unregister_mem(sizeof(sep_list%flux_surfaces),"sep_list%flux_surfaces")
  deallocate(sep_list%flux_surfaces)
  call tr_unregister_mem(sizeof(sep_list%psi_values),"sep_list%psi_values")
  deallocate(sep_list%psi_values)
  

end subroutine bootstrap_get_q_and_ft_splines








!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
!---------------- THIS IS CALLED A LOT IN THE FLUX SURFACE INTEGRATION -----------------------------
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
subroutine compute_surface_basics(surface_list, i_surf, i_piece, rr, s, t, ds, dt)
  
  use data_structure
  implicit none
  ! --- Routine parameters
  type (type_surface_list), intent(in)    :: surface_list
  integer,                  intent(in)    :: i_surf, i_piece
  real*8,                   intent(in)    :: rr
  real*8,                   intent(inout) :: s, t, ds, dt
  
  ! --- Internal parameters
  real*8  :: rr1, rr2, drr1, drr2
  real*8  :: ss1, ss2, dss1, dss2
  
  rr1  = surface_list%flux_surfaces(i_surf)%s(1,i_piece)
  drr1 = surface_list%flux_surfaces(i_surf)%s(2,i_piece)
  rr2  = surface_list%flux_surfaces(i_surf)%s(3,i_piece)
  drr2 = surface_list%flux_surfaces(i_surf)%s(4,i_piece)

  ss1  = surface_list%flux_surfaces(i_surf)%t(1,i_piece)
  dss1 = surface_list%flux_surfaces(i_surf)%t(2,i_piece)
  ss2  = surface_list%flux_surfaces(i_surf)%t(3,i_piece)
  dss2 = surface_list%flux_surfaces(i_surf)%t(4,i_piece)

  call CUB1D(rr1, drr1, rr2, drr2, rr, s, ds)
  call CUB1D(ss1, dss1, ss2, dss2, rr, t, dt)

  return
end subroutine compute_surface_basics







!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
!-------------- CUBIC SPLINE FITTING FUNCTION FOR Q-PROFILE AND FT-PROFILE -------------------------
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
subroutine bootstrap_spline3_coef(n,t,y,z) 
  
  implicit none
  
  integer,              intent(in)   :: n
  real*8, dimension(0:n), intent(in) :: t,y
  real*8, dimension(0:n), intent(out):: z 
  real*8, dimension(0:n-1)           :: h,b
  real*8, dimension(n-1)             :: u,v
  integer                            :: i
  
  do i = 0,n-1
    h(i) = t(i+1) - t(i)
    b(i) = (y(i+1) -y(i))/h(i)    
  end do
  
  u(1) = 2.0*(h(0) + h(1))
  v(1) = 6.0*(b(1) - b(0))
  do i = 2,n-1
    u(i) = 2.0*(h(i) + h(i-1)) - h(i-1)**2/u(i-1)     
    v(i) = 6.0*(b(i) - b(i-1)) - h(i-1)*v(i-1)/u(i-1) 
  end do
  
  z(n) = 0.0  
  do i = n-1,1,-1     
    z(i) = (v(i) - h(i)*z(i+1))/u(i)
  end do
  z(0) = 0.0
  
  return
end subroutine bootstrap_spline3_coef 







!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
!-------------------------------- EVALUATE CUBIC SPLINE --------------------------------------------
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
real*8 function bootstrap_spline3_eval(n,t,y,z,x)
  integer,              intent(in)   :: n
  real*8, dimension(0:n), intent(in) :: t,y,z	 
  real*8,                 intent(in) :: x
  real*8                             :: h, temp
  integer                            :: i
  
  do i = n-1,1,-1     
    if( x - t(i) >= 0.0) exit	 
  end do
  
  h = t(i+1) - t(i)	
  temp = 0.5*z(i) + (x - t(i))*(z(i+1) - z(i))/(6.0*h) 
  temp = (y(i+1) - y(i))/h - h*(z(i+1) + 2.0*z(i))/6.0 + (x- t(i))*temp     
  bootstrap_spline3_eval = y(i) + (x - t(i))*temp  
  
end function bootstrap_spline3_eval 
                                                                   








!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
!-------------------------------- EVALUATE CUBIC SPLINE --------------------------------------------
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
subroutine bootstrap_spline3_eval_all(psi_n, q, ft, B)
  real*8,    intent(in)    :: psi_n
  real*8,    intent(inout) :: q, ft, B
  real*8                   :: h, temp
  integer                  :: n
  integer                  :: i
  
  n = n_spline-1
  do i = n-1,1,-1     
    if( psi_n - psi_knots(i) >= 0.0) exit	 
  end do
  
  h = psi_knots(i+1) - psi_knots(i)
  
  temp = 0.5*q_spline(i) + (psi_n - psi_knots(i))*(q_spline(i+1) - q_spline(i))/(6.0*h) 
  temp = (q_knots(i+1) - q_knots(i))/h - h*(q_spline(i+1) + 2.0*q_spline(i))/6.0 + (psi_n - psi_knots(i))*temp     
  q = q_knots(i) + (psi_n - psi_knots(i))*temp  
  
  temp = 0.5*ft_spline(i) + (psi_n - psi_knots(i))*(ft_spline(i+1) - ft_spline(i))/(6.0*h) 
  temp = (ft_knots(i+1) - ft_knots(i))/h - h*(ft_spline(i+1) + 2.0*ft_spline(i))/6.0 + (psi_n - psi_knots(i))*temp     
  ft = ft_knots(i) + (psi_n - psi_knots(i))*temp  
  
  temp = 0.5*B_spline(i) + (psi_n - psi_knots(i))*(B_spline(i+1) - B_spline(i))/(6.0*h) 
  temp = (B_knots(i+1) - B_knots(i))/h - h*(B_spline(i+1) + 2.0*B_spline(i))/6.0 + (psi_n - psi_knots(i))*temp     
  B = B_knots(i) + (psi_n - psi_knots(i))*temp  
  
  return
end subroutine bootstrap_spline3_eval_all 
                                                                   








!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
!-------- OLD FORMULA FROM WESSON (APPENDIX 14.12) BASED ON A H.WILSON FORMULA FROM 1994 -----------
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
subroutine bootstrap_current_wilson(R, R_axis,    &
                                    psi_axis, psi_bnd,    &
                                    psi_norm,             &
                                    ps0, ps0_x, ps0_y,    &
				    r0,  r0_x,  r0_y,     &
				    Ti0, Ti0_x, Ti0_y,    &
                                    Te0, Te0_x, Te0_y,    &
				    Jb)
!----------------------------------------------------------------------------------------------------------
! calculates the bootstrap current for the RHS of the matrix, based on the Wesson formula (Appendix 14.12)
! Note that we do not linearise all the coeffs C1-4, they are considered as local constants...
!----------------------------------------------------------------------------------------------------------

  use constants
  use phys_module

  implicit none
  ! --- Routine parameters
  real*8, intent(in)  :: R, R_axis
  real*8, intent(in)  :: psi_axis, psi_bnd
  real*8, intent(in)  :: psi_norm
  real*8, intent(in)  :: ps0, ps0_x, ps0_y
  real*8, intent(in)  :: r0,  r0_x,  r0_y
  real*8              :: rho, rho_x, rho_y, drho
  real*8, intent(in)  :: Ti0, Ti0_x, Ti0_y
  real*8              :: Ti,  Ti_x,  Ti_y,  dTi
  real*8, intent(in)  :: Te0, Te0_x, Te0_y
  real*8              :: Te,  Te_x,  Te_y,  dTe
  real*8              :: grad_psi, psi_n, X
  real*8              :: Nue, Nui
  real*8              :: DD, C1, C2, C3, C4, AC4, BC4
  real*8, intent(out) :: Jb
  real*8              :: rho_norm
  real*8              :: tanh_boot, position, width

  ! --- Need central_density
  if (central_density .lt. 1.d-6) then
    write(*,*)'**********************!!WARNING!!*************************'
    write(*,*)'Element_matrix asks for the bootstrap current,'
    write(*,*)'but the central_density is not defined in the input file.'
    write(*,*)'Returning zero bootstrap current...'
    write(*,*)'**********************!!WARNING!!*************************'
    Jb = 0.d0
    return
  endif
  
  ! --- Careful with convention of density (some people get scared when they see an exponent of 19-20 and prefer to just ignore it...)
  if (central_density .lt. 1.d17) then
    rho_norm = central_density*1.d20 ! (this is so clever... I hope they do it like this at NASA...)
  else
    rho_norm = central_density
  endif
  
  ! --- Renormalise temperature and density 
  ! --- Note for us density is n, not mi*n,
  ! --- but that's ok since formula is with p = n*T
  ! --- ie. rho*T for us...
  ! --- Temperature in Joules
  ! --- Density in 1/(cubic meters)
  rho   = r0   * rho_norm
  if (r0 .lt. rho_1*1.d-1) rho = rho_1*1.d-1 * rho_norm
  rho_x = r0_x * rho_norm
  rho_y = r0_y * rho_norm
  Ti    = Ti0   / (MU_ZERO*rho_norm)
  if (Ti0 .lt. Ti_1*1.d-1) Ti = Ti_1*1.d-1 / (MU_ZERO*rho_norm)
  Ti_x  = Ti0_x / (MU_ZERO*rho_norm)
  Ti_y  = Ti0_y / (MU_ZERO*rho_norm)
  Te    = Te0   / (MU_ZERO*rho_norm)
  if (Te0 .lt. Te_1*1.d-1) Te = Te_1*1.d-1 / (MU_ZERO*rho_norm)
  Te_x  = Te0_x / (MU_ZERO*rho_norm)
  Te_y  = Te0_y / (MU_ZERO*rho_norm)
        
  ! --- Psi variables, including r~a*sqrt(psi) and X=sqrt(2*r/R0)
  psi_n    = psi_norm
  if (psi_n .lt. 1.d-1)  psi_n = 1.d-1
  grad_psi = (ps0_x*ps0_x + ps0_y*ps0_y)**0.5d0
  if (grad_psi .lt. 1.d-1) grad_psi = 1.d-1
  X        = sqrt(2.d0*minRad*sqrt(psi_n)/R_axis)

  ! --- Derivatives with respect to psi
  drho = (rho_x*ps0_x + rho_y*ps0_y) / grad_psi**2.d0
  dTi  = (Ti_x *ps0_x + Ti_y *ps0_y) / grad_psi**2.d0
  dTe  = (Te_x *ps0_x + Te_y *ps0_y) / grad_psi**2.d0
        
  ! --- Nue* formula from Wesson : Nue* = R*q / ( eps**(3/2) * (Te/me)**(1/2) * Taue )
  ! --- where Taue is the electron collision time (formula for Nui* is very similar)
  Nui = 5.4d-56 * abs(F0) * rho / ( Ti**2.d0 * X**3.d0 * grad_psi )
  Nue = 9.3d-56 * abs(F0) * rho / ( Te**2.d0 * X**3.d0 * grad_psi )

  ! --- Coefficients
  DD  = 2.4d0 + 5.4d0*X  + 2.6d0*X*X

  C1  = (4.d0 + 2.6d0*X) / ( (1.d0 + 1.02d0*sqrt(Nue) + 1.07d0*Nue) * (1.d0 + 0.38d0*Nue*X**3.d0) )

  C2  = C1*Ti/Te
  
  C3  = (7.d0 + 6.5d0*X) / ( (1.d0 + 0.57d0*sqrt(Nue) + 0.61d0*Nue) * (1.d0 + 0.22d0*Nue*X**3.d0) ) - 2.5d0*C1

  AC4 = (  -1.17d0/(1.d0+0.46d0*X) + 0.35d0*sqrt(Nui)) / (1.d0 + 0.7d0*sqrt(Nui)) + 0.26d0*Nui*Nui*X**6.d0
  BC4 = (1.d0 - 0.125d0*Nui*Nui*X**6.d0) * (1.d0 + 0.125d0*Nue*Nue*X**6.d0)
  C4  = AC4*C2 / BC4

  ! --- Bootstrap Current
  Jb = R_axis**2.d0*X*rho*Te/DD * ( C1*(dTe/Te+drho/rho) + C2*(dTi/Ti+drho/rho) + C3*dTe/Te + C4*dTi/Ti)

  ! --- Current with denormalisation
  Jb = -Jb * MU_ZERO
  
  ! --- There should not be any bootstrap outside plasma...
  position  = max(rho_coef(5),1.d0) + 2.d0 * rho_coef(4)
  width     = rho_coef(4) / 2.d0
  tanh_boot = 0.5d0 - 0.5d0 * tanh( (psi_norm - position)/width ) 
  Jb = Jb * tanh_boot


return
end subroutine bootstrap_current_wilson



end module mod_bootstrap_functions

