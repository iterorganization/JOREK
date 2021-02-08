subroutine current(xpoint2,xcase2,R,Z,Z_xpoint,psi,psi_axis,psi_bnd,zjz)
!-----------------------------------------------------------------------
! Determine the current at a given position from the density,
! temperature, and FF' input profiles
!-----------------------------------------------------------------------

use mod_parameters
use phys_module

implicit none

! --- Routine parameters
logical, intent(in)    :: xpoint2
integer, intent(in)    :: xcase2
real*8,  intent(in)    :: R, Z
real*8,  intent(in)    :: Z_xpoint(2)
real*8,  intent(in)    :: psi
real*8,  intent(in)    :: psi_axis
real*8,  intent(in)    :: psi_bnd
real*8,  intent(out)   :: zjz        ! Current at the given position.

! --- local variables
real*8  :: psi_n
real*8  :: zn,dn_dpsi,dn_dz,dn_dpsi2,dn_dz2,dn_dpsi_dz,dn_dpsi3,dn_dpsi_dz2,dn_dpsi2_dz
real*8  :: zT,dT_dpsi,dT_dz,dT_dpsi2,dT_dz2,dT_dpsi_dz,dT_dpsi3,dT_dpsi_dz2,dT_dpsi2_dz
real*8  :: zTi,zTe,dTi_dpsi,dTe_dpsi,dTi_dz,dTe_dz,dTi_dpsi2,dTe_dpsi2,dTi_dz2,dTe_dz2
real*8  :: dTi_dpsi_dz,dTe_dpsi_dz,dTi_dpsi3,dTe_dpsi3,dTi_dpsi_dz2,dTe_dpsi_dz2,dTi_dpsi2_dz,dTe_dpsi2_dz
real*8  :: zFFprime, dFFprime_dpsi, dFFprime_dz, dFFprime_dpsi_dz,dFFprime_dpsi2,dFFprime_dz2

psi_n = (psi - psi_axis) / (psi_bnd - psi_axis)

call density(    xpoint2, xcase2, Z, Z_xpoint, psi,psi_axis,psi_bnd,&
             zn,dn_dpsi,dn_dz,dn_dpsi2, dn_dz2, dn_dpsi_dz,dn_dpsi3,dn_dpsi_dz2,dn_dpsi2_dz)

if (with_TiTe) then
  
  call temperature_i(xpoint2, xcase2, Z, Z_xpoint, psi,psi_axis,psi_bnd, &
                   zTi,dTi_dpsi,dTi_dz,dTi_dpsi2,dTi_dz2,dTi_dpsi_dz,dTi_dpsi3,dTi_dpsi_dz2,dTi_dpsi2_dz)
		 
  call temperature_e(xpoint2, xcase2, Z, Z_xpoint, psi,psi_axis,psi_bnd, &
                   zTe,dTe_dpsi,dTe_dz,dTe_dpsi2,dTe_dz2,dTe_dpsi_dz,dTe_dpsi3,dTe_dpsi_dz2,dTe_dpsi2_dz)

  zT = zTi + zTe
  dT_dpsi = dTi_dpsi + dTe_dpsi

else
  
  call temperature(xpoint2, xcase2, Z, Z_xpoint, psi,psi_axis,psi_bnd, &
                   zT,dT_dpsi,dT_dz,dT_dpsi2,dT_dz2,dT_dpsi_dz,dT_dpsi3,dT_dpsi_dz2,dT_dpsi2_dz)

endif

call FFprime(    xpoint2, xcase2, Z, Z_xpoint, psi,psi_axis,psi_bnd, &
             zFFprime,dFFprime_dpsi,dFFprime_dz,dFFprime_dpsi2,dFFprime_dz2, dFFprime_dpsi_dz, .true.)

zjz   = zFFprime - R*R * (zn * dT_dpsi + dn_dpsi * zT)

!if ((bootstrap) .and. (restart)) then
!  zjz   = zjz * (0.5d0 - 0.5d0* tanh( (psi_n - (FF_coef(7)-FF_coef(8)))/FF_coef(8) ) )
!endif

return
end subroutine current
