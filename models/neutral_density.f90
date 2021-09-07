subroutine neutral_density(xpoint2,xcase2,Z,Z_xpoint,psi,psi_axis,psi_bnd,density_profile,dn_dpsi,dn_dz, &
                           dn_dpsi2,dn_dz2,dn_dpsi_dz,dn_dpsi3,dn_dpsi_dz2, dn_dpsi2_dz)
!-----------------------------------------------------------------------
! Determines the neutral density value and its derivatives at the given
! position (Z, psi) from the analytical or numerical input profile.
!-----------------------------------------------------------------------
use phys_module

implicit none

! --- Routine parameters
logical, intent(in)  :: xpoint2
integer, intent(in)  :: xcase2
real*8,  intent(in)  :: Z, Z_xpoint(2), psi, psi_axis, psi_bnd
real*8,  intent(out) :: density_profile, dn_dpsi, dn_dz, dn_dpsi2, dn_dz2, &
                        dn_dpsi_dz, dn_dpsi3, dn_dpsi_dz2, dn_dpsi2_dz

! --- Internal variables.
real*8  :: prof0, prof1, dprof0_dpsi, dprof0_dpsi2, dprof0_dpsi3, psi_barrier
real*8  :: psi_n, psi_star, delta_psi, sig_n, sigz, dprof1_dpsi, dprof1_dpsi2, dprof1_dpsi3
real*8  :: atn, datn, d2atn, d3atn
real*8  :: atn_z,   datn_z,   d2atn_z
real*8  :: atn_z_u, datn_z_u, d2atn_z_u, factor
real*8  :: cosh1, cosh2, cosh3, cosh3_u
real*8  :: tanh1, tanh2, tanh2_u
real*8  :: Ztan_pos
! for interpolating numerical profiles
integer :: left, right, mid
real*8  :: aux1, aux2, Z_star, Z_star_u

delta_psi = psi_bnd - psi_axis
psi_n     = (psi - psi_axis) / delta_psi

psi_n = max( min(psi_n, 2.), 0. )

!factor = 1.d0
!if (xpoint2) then
!  if ((Z .lt. Z_xpoint) .and. (psi_n .lt. 1.d0) ) then
!    psi_n = 2.d0 - psi_n
!    factor = -1.d0
!  endif
!endif

! --- Profile as a function of Psi_N.
if ( .not. num_rhon ) then ! use analytical representation
  
  prof0        = (rhon_0-rhon_1)*(1.d0 + rhon_coef(1) * psi_n + rhon_coef(2) * psi_n**2 + rhon_coef(3) * psi_n**3)
  dprof0_dpsi  = (rhon_0-rhon_1)*(rhon_coef(1) + 2.d0 * rhon_coef(2) * psi_n + 3.d0 * rhon_coef(3) * psi_n**2) / delta_psi
  dprof0_dpsi2 = (rhon_0-rhon_1)*(2.d0 * rhon_coef(2) + 6.d0 * rhon_coef(3) * psi_n)                           / delta_psi**2
  dprof0_dpsi3 = (rhon_0-rhon_1)*(6.d0 * rhon_coef(3))                                                         / delta_psi**3
  
  sig_n       = rhon_coef(4)
  psi_barrier = rhon_coef(5)
  
  psi_star = (psi_n - psi_barrier)/sig_n
  psi_star = min( max( psi_star, -40.d0), 40.d0) ! avoid floating-point exceptions
  
  tanh1 = tanh(psi_star)
  cosh1 = cosh(psi_star)
  cosh2 = cosh(2.d0*psi_star)
  
  atn   = (0.5d0 - 0.5d0*tanh1)
  datn  = - 1.d0/cosh1**2 / (2.d0 * sig_n) / delta_psi
  d2atn =   1.d0/cosh1**2 / sig_n**2 * tanh1 / delta_psi**2
  d3atn = - 1.d0/cosh1**4 / sig_n**3 * (-2.d0 + cosh2) / delta_psi**3
  
  prof1        = prof0        * atn
  dprof1_dpsi  = dprof0_dpsi  * atn +         prof0       * datn
  dprof1_dpsi2 = dprof0_dpsi2 * atn + 2.d0 * dprof0_dpsi  * datn + prof0              * d2atn
  dprof1_dpsi3 = dprof0_dpsi3 * atn + 3.d0 * dprof0_dpsi2 * datn + 3.d0 * dprof0_dpsi * d2atn + prof0 * d3atn
  
else ! use numerical representation.
  
  ! --- Interpolate profile and derivatives to position psi_n by bisections.
  left  = 1
  right = num_rhon_len
  do
    if ( right == left + 1 ) exit
    mid = (left + right) / 2
    if ( num_rhon_x(mid) >= psi_n ) then
      right = mid
    else
      left = mid
    end if
  end do
  aux1 = (psi_n - num_rhon_x(left)) / (num_rhon_x(right) - num_rhon_x(left))
  aux2 = (1. - aux1)
  prof1        = num_rhon_y0(left)   * aux2 + num_rhon_y0(right) * aux1
  dprof1_dpsi  = ( num_rhon_y1(left) * aux2 + num_rhon_y1(right) * aux1 ) / delta_psi
  dprof1_dpsi2 = ( num_rhon_y2(left) * aux2 + num_rhon_y2(right) * aux1 ) / delta_psi**2
  dprof1_dpsi3 = ( num_rhon_y3(left) * aux2 + num_rhon_y3(right) * aux1 ) / delta_psi**3
  
end if

! --- Additional explicit dependence of the profile on Z to ensure that the profile is
!     approximately zero in the private flux region below the x-point.
!     ***************
!     ***************
!     ***************
!     IMPORTANT NOTE: Everything is the same as the usual density profile, except for
!                     this Z-tanh at the X-points. Instead of setting the profile to rhon_1
!                     in the divertor regions, we set it to the coefs
!                     rhon_coef(9) in the lower divertor
!                     and rhon_coef(10) in the upper divertor
!                     For the tanh width, rhon_coef(8) is used for both divertors
!                     For the tanh Z-position, rhon_coef(6) and rhon_coef(7) are used (lower and upper respectively)
!                     If rhon_coef(6) and/or rhon_coef(7) are 0.d0, then Z_xpoint(:) are used instead.
if (xpoint2) then
  
  sigz            = rhon_coef(8)

  if (xcase2 .eq. LOWER_XPOINT) then
    atn_z_u   = 1.d0
    datn_z_u  = 0.d0
    d2atn_z_u = 0.d0
  else
    Ztan_pos  = rhon_coef(7)
    if (Ztan_pos .eq. 0.d0) Ztan_pos = Z_xpoint(2)
    Z_star_u  = (Z-Ztan_pos)/sigz
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
    Ztan_pos  = rhon_coef(6)
    if (Ztan_pos .eq. 0.d0) Ztan_pos = Z_xpoint(1)
    Z_star  = (Ztan_pos-Z)/sigz
    Z_star  = min( max( Z_star, -40.d0), 40.d0) ! avoid floating-point exceptions

    tanh2   = tanh(Z_star)
    cosh3   = cosh(Z_star)
      
    atn_z   = (0.5d0 - 0.5d0*tanh2)
    datn_z  =  0.5d0/cosh3**2 / sigz
    d2atn_z =  1.0d0/cosh3**2 / sigz**2 * tanh2
  endif
  
  density_profile = prof1        + (1.0 -   atn_z) * (rhon_coef(9)-rhon_1 -  prof1      ) + (1.0 -   atn_z_u) * (rhon_coef(10)-rhon_1 -  prof1      )
  dn_dpsi         = dprof1_dpsi  + (1.0 -   atn_z) * (                    - dprof1_dpsi ) + (1.0 -   atn_z_u) * (                     - dprof1_dpsi )
  dn_dpsi2        = dprof1_dpsi2 + (1.0 -   atn_z) * (                    - dprof1_dpsi2) + (1.0 -   atn_z_u) * (                     - dprof1_dpsi2)
  dn_dpsi3        = dprof1_dpsi3 + (1.0 -   atn_z) * (                    - dprof1_dpsi3) + (1.0 -   atn_z_u) * (                     - dprof1_dpsi3)
  dn_dz           = prof1        + (    -  datn_z) * (rhon_coef(9)-rhon_1 -  prof1      ) + (    -  datn_z_u) * (rhon_coef(10)-rhon_1 -  prof1      )
  dn_dz2          = prof1        + (    - d2atn_z) * (rhon_coef(9)-rhon_1 -  prof1      ) + (    - d2atn_z_u) * (rhon_coef(10)-rhon_1 -  prof1      )
  dn_dpsi_dz      = dprof1_dpsi  + (    -  datn_z) * (                    - dprof1_dpsi ) + (    -  datn_z_u) * (                     - dprof1_dpsi )
  dn_dpsi_dz2     = dprof1_dpsi  + (    - d2atn_z) * (                    - dprof1_dpsi ) + (    - d2atn_z_u) * (                     - dprof1_dpsi )
  dn_dpsi2_dz     = dprof1_dpsi2 + (    -  datn_z) * (                    - dprof1_dpsi2) + (    -  datn_z_u) * (                     - dprof1_dpsi2)

else
  
  density_profile = prof1
  dn_dpsi     = dprof1_dpsi!   * factor
  dn_dpsi2    = dprof1_dpsi2
  dn_dpsi3    = dprof1_dpsi3!  * factor
  dn_dz       = 0.d0
  dn_dz2      = 0.d0
  dn_dpsi_dz  = 0.d0
  dn_dpsi2_dz = 0.d0
  dn_dpsi_dz2 = 0.d0

end if

density_profile = density_profile + rhon_1

return
end subroutine neutral_density
