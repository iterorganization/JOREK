recursive subroutine FFprime(xpoint2,xcase2,Z,Z_xpoint,psi,psi_axis,psi_bnd,FFprime_profile,dFF_dpsi,dFF_dz, &
                       dFF_dpsi2,dFF_dz2,dFF_dpsi_dz, from_F_profile_710)
!-----------------------------------------------------------------------
! Determines the F*F' value and its derivatives at the given
! position (Z, psi) from the analytical or numerical input profile.
!-----------------------------------------------------------------------
use phys_module
use vacuum, only: current_FB_fact
use mod_F_profile

implicit none

! --- Routine parameters
logical, intent(in)  :: xpoint2
integer, intent(in)  :: xcase2
real*8,  intent(in)  :: Z, Z_xpoint(2), psi, psi_axis, psi_bnd
real*8,  intent(out) :: FFprime_profile, dFF_dpsi, dFF_dz, dFF_dpsi2, dFF_dz2, dFF_dpsi_dz
logical, intent(in)  :: from_F_profile_710

! --- Internal variables.
real*8  :: prof0, prof1, dprof0_dpsi, dprof0_dpsi2, psi_barrier
real*8  :: psi_n, psi_star, delta_psi, sig_F, sigz, dprof1_dpsi, dprof1_dpsi2
real*8  :: atn,     datn,     d2atn, d3atn
real*8  :: atn_z,   datn_z,   d2atn_z
real*8  :: atn_z_u, datn_z_u, d2atn_z_u, factor
real*8  :: d_pert, d2_pert, d3_pert
real*8  :: cosh1, cosh2, cosh3, cosh3_u
real*8  :: tanh1, tanh2, tanh2_u
! for interpolating numerical profiles
integer :: left, right, mid
real*8  :: aux1, aux2, Z_star, Z_star_u
real*8  :: F_prof, dF_dpsi, dF_dz, dF_dpsi2, dF_dz2, dF_dpsi_dz
real*8  :: no_delta_psi, ffprime_out, FFprime_profile2, F_prof2


! --- the F-profile and FFprime need to be coherent. Always!!!
#ifdef fullmhd
  if (from_F_profile_710) then
    ! --- Call function
    call F_profile(xpoint2, xcase2, Z, Z_xpoint, psi,psi_axis,psi_bnd, &
                   F_prof,          dF_dpsi,  dF_dz,  dF_dpsi2,  dF_dz2,  dF_dpsi_dz , &
                   FFprime_profile, dFF_dpsi, dFF_dz, dFF_dpsi2, dFF_dz2, dFF_dpsi_dz)

    ! --- Because JOREK uses a negative FF' in the GS-equation and the current-routines
    ! --- But in Full-MHD, because we need to integrate FF', we can't do this, so we use the real F-profile and FF',
    ! --- and then reverse it for all the routines that use it.
    FFprime_profile = - FFprime_profile
    dFF_dpsi        = - dFF_dpsi
    dFF_dz          = - dFF_dz
    dFF_dpsi2       = - dFF_dpsi2
    dFF_dz2         = - dFF_dz2
    dFF_dpsi_dz     = - dFF_dpsi_dz

    return
  endif
#endif


delta_psi = psi_bnd - psi_axis
psi_n     = (psi - psi_axis) / delta_psi
no_delta_psi = 1.d0
if (FF_coef(9) .eq. 1.d0) no_delta_psi = delta_psi

psi_n = max( min(psi_n, 2.), 0. )

!factor = 1.d0
!if (xpoint2) then
!  if ((Z .lt. Z_xpoint) .and. (psi_n .lt. 1.d0) ) then
!    psi_n = 2.d0 - psi_n
!    factor = -1.d0
!  endif
!endif

! --- Profile as a function of Psi_N.
if ( .not. num_ffprime ) then ! use analytical representation
  
  d_pert  = + FF_coef(6)/cosh((psi_n - FF_coef(7))/FF_coef(8))**2 / (2.d0 * FF_coef(8)) / delta_psi * no_delta_psi
  d2_pert = - FF_coef(6)/cosh((psi_n - FF_coef(7))/FF_coef(8))**2 / (FF_coef(8)**2)  &
            * tanh((psi_n - FF_coef(7))/FF_coef(8)) / delta_psi**2 * no_delta_psi
  d3_pert = + FF_coef(6)/cosh((psi_n - FF_coef(7))/FF_coef(8))**4 / (FF_coef(8)**3)  &
            * (-2.d0 + cosh(2.d0*(psi_n-FF_coef(7))/FF_coef(8)) ) / delta_psi**3 * no_delta_psi
  
  prof0        = (FF_0 - FF_1) * ( 1.d0 + FF_coef(1) * psi_n + FF_coef(2) * psi_n**2 + FF_coef(3) * psi_n**3)
  dprof0_dpsi  = (FF_0 - FF_1) * ( FF_coef(1) + 2.d0 * FF_coef(2) * psi_n + 3.d0 * FF_coef(3) * psi_n**2)    / delta_psi
  dprof0_dpsi2 = (FF_0 - FF_1) * (2.d0 * FF_coef(2) + 6.d0 * FF_coef(3) * psi_n) / delta_psi**2
  
  prof0        = prof0        + d_pert
  dprof0_dpsi  = dprof0_dpsi  + d2_pert
  dprof0_dpsi2 = dprof0_dpsi2 + d3_pert
  
  sig_F        = FF_coef(4)
  psi_barrier  = FF_coef(5)
  
  psi_star = (psi_n - psi_barrier)/sig_F
  psi_star = min( max( psi_star, -40.d0), 40.d0) ! avoid floating-point exceptions
  
  tanh1 = tanh(psi_star)
  cosh1 = cosh(psi_star)
  cosh2 = cosh(2.d0*psi_star)
  
  atn   = (0.5d0 - 0.5d0*tanh1)
  datn  = - 1.d0/cosh1**2 / (2.d0 * sig_F) / delta_psi
  d2atn =   1.d0/cosh1**2 / sig_F**2 * tanh1 / delta_psi**2
  d3atn = - 1.d0/cosh1**4 / sig_F**3 * (-2.d0 + cosh2) / delta_psi**3
  
  prof1        = prof0        * atn
  dprof1_dpsi  = dprof0_dpsi  * atn +         prof0       * datn
  dprof1_dpsi2 = dprof0_dpsi2 * atn + 2.d0 * dprof0_dpsi  * datn + prof0       * d2atn
  
else ! use numerical representation.
  
  ! --- Interpolate profile and derivatives to position psi_n by bisections.
  left  = 1
  right = num_ffprime_len
  do
    if ( right == left + 1 ) exit
    mid = (left + right) / 2
    if ( num_ffprime_x(mid) >= psi_n ) then
      right = mid
    else
      left = mid
    end if
  end do
  aux1 = (psi_n - num_ffprime_x(left)) / (num_ffprime_x(right) - num_ffprime_x(left))
  aux2 = (1. - aux1)
  prof1        = num_ffprime_y0(left)   * aux2 + num_ffprime_y0(right) * aux1
  dprof1_dpsi  = ( num_ffprime_y1(left) * aux2 + num_ffprime_y1(right) * aux1 ) / delta_psi
  dprof1_dpsi2 = ( num_ffprime_y2(left) * aux2 + num_ffprime_y2(right) * aux1 ) / delta_psi**2
  
end if

! --- Additional explicit dependence of the profile on Z to ensure that the profile is
!     approximately zero in the private flux region below the x-point.
if ( xpoint2 ) then
  
  sigz = 0.1d0

  if (xcase2 .eq. LOWER_XPOINT) then
    atn_z_u   = 1.d0
    datn_z_u  = 0.d0
    d2atn_z_u = 0.d0
  else
    Z_star_u  = (Z-Z_xpoint(2))/sigz
    Z_star_u  = min( max( Z_star_u, -40.d0), 40.d0) ! avoid floating-point exceptions
    
    tanh2_u   = tanh(Z_star_u)
    cosh3_u   = cosh(Z_star_u)

    atn_z_u   = (0.5d0 - 0.5d0*tanh2_u)
    datn_z_u  = -0.5d0/cosh3_u**2 / sigz
    d2atn_z_u =  1.0d0/cosh3_u**2 / sigz**2 * tanh2_u
  endif
  if (xcase2 .eq. UPPER_XPOINT) then
    atn_z   = 1.d0
    datn_z  = 0.d0
    d2atn_z = 0.d0
  else
    Z_star  = (Z_xpoint(1)-Z)/sigz
    Z_star  = min( max( Z_star, -40.d0), 40.d0) ! avoid floating-point exceptions

    tanh2   = tanh(Z_star)
    cosh3   = cosh(Z_star)
      
    atn_z   = (0.5d0 - 0.5d0*tanh2)
    datn_z  =  0.5d0/cosh3**2   / sigz
    d2atn_z =  1.0d0/cosh3**2   / sigz**2 * tanh2
  endif 
  
  FFprime_profile  = prof1        *    atn_z * atn_z_u
  dFF_dpsi         = dprof1_dpsi  *    atn_z * atn_z_u
  dFF_dpsi2        = dprof1_dpsi2 *    atn_z * atn_z_u  
  dFF_dz           = prof1        * ( datn_z * atn_z_u +         atn_z * datn_z_u)
  dFF_dz2          = prof1        * (d2atn_z * atn_z_u + 2.d0 * datn_z * datn_z_u  +  atn_z * d2atn_z_u) 
  dFF_dpsi_dz      = dprof1_dpsi  * ( datn_z * atn_z_u +         atn_z * datn_z_u)

else
 
  FFprime_profile  = prof1
  dFF_dpsi         = dprof1_dpsi !   * factor
  dFF_dpsi2        = dprof1_dpsi2
  dFF_dz           = 0.d0
  dFF_dz2          = 0.d0
  dFF_dpsi_dz      = 0.d0
  
end if

if (freeboundary_equil .and. num_ffprime) then            !if the ffprime profile is given in a file and freeboundary equilibrium is on,
                                                         !the full profile is multiplied by a factor in order to iterate to a given current   
   FFprime_profile  = FFprime_profile  * current_FB_fact
   dFF_dpsi         = dFF_dpsi         * current_FB_fact
   dFF_dpsi2        = dFF_dpsi2        * current_FB_fact
   dFF_dz           = dFF_dz           * current_FB_fact
   dFF_dz2          = dFF_dz2          * current_FB_fact
   dFF_dpsi_dz      = dFF_dpsi_dz      * current_FB_fact

end if

FFprime_profile = FFprime_profile + FF_1

return
end subroutine FFprime
