!> Module to calculate ionization, recombination and radiation coefficients for Deuterium 
module mod_atomic_coeff_deuterium 

use mod_openadas
use constants
use phys_module, only: central_density, central_mass, gamma, deuterium_adas, deuterium_adas_1e20, old_deuterium_atomic

implicit none

type(ADF11_all) :: ad_deuterium !< ADAS structure for deuterium 

! --- Fit coefficients
real*8, parameter :: A5_ion= 0.018942148, A4_ion=-0.283846021, A3_ion= 1.694045818, A2_ion=-5.164073876, A1_ion= 7.925697878, A0_ion=-12.17397177
real*8, parameter :: A5_rec=-0.000523617, A4_rec= 0.037488784, A3_rec=-0.312405179, A2_rec= 0.883520589, A1_rec=-2.046497531, A0_rec=-11.74749202
real*8, parameter :: A5_zrb= 0.000692804, A4_zrb= 0.009293147, A3_zrb=-0.191315426, A2_zrb= 1.047258967, A1_zrb=-1.981527194, A0_zrb=-29.44512748
real*8, parameter :: A5_zlt= 0.026279151, A4_zlt=-0.367640092, A3_zlt= 1.992770100, A2_zlt=-5.378600323, A1_zlt= 7.350832521, A0_zlt=-29.45315695
real*8, parameter :: S_ion_puiss = 3.9d-1

private
public atomic_coeff_deuterium, ad_deuterium, plot_atomic_coefficients


contains

! --- This routine gives ionization, recombination and radiation coefficients for Deuterium
! ---   * The input value is the JOREK normalized electron temperature
! ---   * Outputs are the normalized coefficients
! ---   * NOTE THAT THE DERIVATIVES ARE WITH RESPECT TO THE ELECTRON TEMPERATURE (not T=Te+Ti)
! ---   * The coeffiencts are calculated for ne = 1.e20  m^-3 for the fits
subroutine atomic_coeff_deuterium(Te0, Sion_T, dSion_dT, Srec_T, dSrec_dT, LradDcont_T, dLradDcont_dT, LradDrays_T, dLradDrays_dT, ne0 ) 

  implicit none

  ! --- Routine parameters
  real*8, intent(in)    :: Te0                        ! Electron temperature in JOREK units
  real*8, intent(inout) :: Sion_T, dSion_dT           ! Normalized ionization coefficient and its temperature derivative
  real*8, intent(inout) :: Srec_T, dSrec_dT           ! Normalized recombination coefficient and its temperature derivative
  real*8, intent(inout) :: LradDcont_T, dLradDcont_dT ! Normalized Bremss and recomb radiation coefficient and its temperature derivative
  real*8, intent(inout) :: LradDrays_T, dLradDrays_dT ! Normalized line radiation coefficient and its temperature derivative
  real*8, optional, intent(in) :: ne0                 ! Electron density in JOREK units (used only for ADAS data)

  ! --- Local
  real*8 :: coef_ion_1, coef_ion_2, coef_ion_3, T0 
  real*8 :: coef_rad_1, coef_rec_1
  real*8 :: rho_norm, t_norm
  real*8 :: Te_eV, Te_evL10, dTe_eVL10_dT0, Te_eV_lim, Te_si_log10, ne_si_log10, ne_si
  real*8 :: Sion_si, dSion_si                     
  real*8 :: Srec_si, DSrec_si                     
  real*8 :: Szlt_T, dSzlt_dT, Szlt_si, dSzlt_si   
  real*8 :: Szrb_T, dSzrb_dT, Szrb_si, dSzrb_si   
  real*8 :: gamma_factor
  real*8 :: ion_log10, dion_log10
  real*8 :: rec_log10, drec_log10
  real*8 :: zlt_log10, dzlt_log10
  real*8 :: zrb_log10, dzrb_log10

  ! --- Normalization constants
  rho_norm     = central_density*1.d20 * central_mass * MASS_PROTON
  t_norm       = sqrt(MU_zero*rho_norm)
  gamma_factor = gamma-1.d0  ! Normalization factor to include terms in the pressure equation
                             ! internal_energy = pressure / (gamma - 1)

  Te_eV        = Te0 / (EL_CHG*MU_ZERO*central_density*1.d20)

  if ( (.not. old_deuterium_atomic) .and. (.not. deuterium_adas) ) then 

    ! Fit of ADAS data calculated by Guido, see comparisons in 
    ! see the public JOREK issue tracker

    Te_eV_lim = max(Te_eV,     0.2d0 )  ! ADAS fit valid between 0.2eV and 10 keV  
    Te_eV_lim = min(Te_eV_lim, 1.d4  )  
    Te_evL10  = dlog10(Te_ev_lim)

    dTe_eVL10_dT0 = 1.d0 / ( dlog(10.d0) * max(Te0, 0.2d0*EL_CHG*MU_ZERO*central_density*1.d20)  )

    ! --- Ionization
    Ion_log10   = A0_ion + A1_ion*Te_evL10 +      A2_ion*Te_evL10**2 +      A3_ion*Te_evL10**3 +      A4_ion*Te_evL10**4 +      A5_ion*Te_evL10**5
    dIon_log10  =          A1_ion          + 2.d0*A2_ion*Te_evL10    + 3.d0*A3_ion*Te_evL10**2 + 4.d0*A4_ion*Te_evL10**3 + 5.d0*A5_ion*Te_evL10**4

    Sion_si  = 10.d0**(Ion_log10 - 6.d0)                      ! ionisation rate at density 10^20 in m^3 /s
    dSion_si = 10.d0**(Ion_log10 - 6.d0) * alog(10.d0) * dion_log10
    
    Sion_T   = Sion_si  * t_norm * central_density * 1.d20
    dSion_dT = dSion_si * t_norm * central_density * 1.d20 * dTe_eVL10_dT0

    ! --- Recombination
    Rec_log10   = A0_rec + A1_rec*Te_evL10 +      A2_rec*Te_evL10**2 +      A3_rec*Te_evL10**3 +      A4_rec*Te_evL10**4 +      A5_rec*Te_evL10**5
    dRec_log10  =          A1_rec          + 2.d0*A2_rec*Te_evL10    + 3.d0*A3_rec*Te_evL10**2 + 4.d0*A4_rec*Te_evL10**3 + 5.d0*A5_rec*Te_evL10**4

    Srec_si  = 10.d0**(Rec_log10 - 6.d0)                      ! recombination rate at density 10^20 in m^3 /s
    dSrec_si = 10.d0**(Rec_log10 - 6.d0) * alog(10.d0) * drec_log10
    
    Srec_T   = Srec_si  * t_norm * central_density * 1.d20
    dSrec_dT = dSrec_si * t_norm * central_density * 1.d20 * dTe_eVL10_dT0

    ! --- Line radiation
    zlt_log10   = A0_zlt + A1_zlt*Te_evL10 +      A2_zlt*Te_evL10**2 +      A3_zlt*Te_evL10**3 +      A4_zlt*Te_evL10**4 +      A5_zlt*Te_evL10**5
    dzlt_log10  =          A1_zlt          + 2.d0*A2_zlt*Te_evL10    + 3.d0*A3_zlt*Te_evL10**2 + 4.d0*A4_zlt*Te_evL10**3 + 5.d0*A5_zlt*Te_evL10**4

    Szlt_si  = 10.d0**(zlt_log10 - 6.d0)                      !line radiation rate at density 10^20 in W m^3 /s
    dSzlt_si = 10.d0**(zlt_log10 - 6.d0) * alog(10.d0) * dzlt_log10
    
    Szlt_T   = Szlt_si  * (central_density * 1.d20)**2 * MU_zero * t_norm * gamma_factor                  
    dSzlt_dT = dSzlt_si * (central_density * 1.d20)**2 * MU_zero * t_norm * gamma_factor * dTe_eVL10_dT0

    LradDrays_T   = Szlt_T
    dLradDrays_dT = dSzlt_dT

    ! --- Bremstrahlung and recomnbination radiation
    zrb_log10   = A0_zrb + A1_zrb*Te_evL10 +      A2_zrb*Te_evL10**2 +      A3_zrb*Te_evL10**3 +      A4_zrb*Te_evL10**4 +      A5_zrb*Te_evL10**5
    dzrb_log10  =          A1_zrb          + 2.d0*A2_zrb*Te_evL10    + 3.d0*A3_zrb*Te_evL10**2 + 4.d0*A4_zrb*Te_evL10**3 + 5.d0*A5_zrb*Te_evL10**4

    Szrb_si  = 10.d0**(zrb_log10 - 6.d0)                      ! recombination and Bremstrahlung rate at density 10^20 in W m^3 /s
    dSzrb_si = 10.d0**(zrb_log10 - 6.d0) * alog(10.d0) * dzrb_log10
    
    Szrb_T   = Szrb_si  * (central_density * 1.d20)**2 * MU_zero * t_norm * gamma_factor                  
    dSzrb_dT = dSzrb_si * (central_density * 1.d20)**2 * MU_zero * t_norm * gamma_factor * dTe_eVL10_dT0

    LradDcont_T   = Szrb_T 
    dLradDcont_dT = dSzrb_dT

    if ( Te_eV < 0.2d0) then  ! --- Don't radiate or ionize below 0.2 eV, recombination allowed
      LradDcont_T   = 0.d0
      dLradDcont_dT = 0.d0
      LradDrays_T   = 0.d0
      dLradDrays_dT = 0.d0
      dSrec_dT      = 0.d0
      Sion_T        = 0.d0
      dSion_dT      = 0.d0
    endif

    if ( Te_eV > 1.d4) then   ! --- Fix values beyond 10 keV and remove derivatives
      dLradDcont_dT = 0.d0
      dLradDrays_dT = 0.d0
      dSrec_dT      = 0.d0
      dSion_dT      = 0.d0
    endif
    
  else if (.not. deuterium_adas) then  ! --- use old fitting formula

    T0 = 2.d0 * Te0   ! The total temperature was used for these old coefficients...

    ! --- Te_max=10 keV, beyond that value the fits blow up
    T0 = min( T0, 2.d0 * 1.d4 * EL_CHG*MU_ZERO*central_density*1.d20 )

    coef_ion_1  = sqrt(MU_ZERO*central_mass*MASS_PROTON) * (central_density*1.d20)**(1.5d0) * 0.2917d-13
    coef_ion_2  = 0.232d0
    coef_ion_3  = EL_CHG*MU_ZERO*central_density*1.d20 * 27.2d0

    
    coef_rad_1  = gamma_factor*MU_ZERO**1.5d0*(central_mass*MASS_PROTON)**0.5d0*(central_density*1.d20)**2.5d0
    coef_rec_1  = (MU_ZERO*central_mass*MASS_PROTON)**(0.5d0)*(central_density*1.d20)**(1.5d0)   

    ! -------------------------------------------
    ! --- Ionization rate for Deuterium
    ! --- (see Wiki for more info: http://jorek.eu/wiki/doku.php?id=model500_501_555#ionization_rate_for_deuterium)
    ! ------------------------------------------- 
    if (Te_eV .gt. 0.1) then
      Sion_T   = coef_ion_1*((coef_ion_3/T0)**S_ion_puiss)*1/(coef_ion_2+coef_ion_3/T0)*exp(-coef_ion_3/T0)
      dSion_dT = Sion_T * ( -S_ion_puiss/T0 + coef_ion_3/(T0*(coef_ion_2*T0+coef_ion_3)) + coef_ion_3*T0**(-2.d0) )
    else
      Sion_T   = 0.
      dSion_dT = 0. 
    endif

    !-------------------------------------------------
    ! --- Recombination rate for ionized Deuterium
    ! (see Wiki for more info: http://jorek.eu/wiki/doku.php?id=model500_501_555#recombination_rate_for_deuterium)
    !-------------------------------------------------
    Srec_T    =            coef_rec_1 * 0.7d-19 * (13.6d0*(2.d0*EL_CHG*MU_ZERO*central_density*1.d20))**(0.5d0) * (max(T0,1.d-6)/(2.d0))**(-0.5d0)      
    dSrec_dT  = - 0.25d0 * coef_rec_1 * 0.7d-19 * (13.6d0*(2.d0*EL_CHG*MU_ZERO*central_density*1.d20))**(0.5d0) * (max(T0,1.d-6)/(2.d0))**(-1.5d0) 

 
    if (T0 .gt. 1.d-6) then

      LradDcont_T = coef_rad_1*5.37d-37*(1.d1)**(-1.5d0)*(1.d0)**2*sqrt(Te_eV) ! Only Bremsstrahlung contribution
    
      dLradDcont_dT = coef_rad_1*5.37d-37*(1.d1)**(-1.5d0)*(1.d0)**2*(8*EL_CHG*MU_ZERO*central_density*1.d20*sqrt(Te_eV))**(-1.d0) 
    
    
      LradDrays_T = coef_rad_1*(1.d1)**(-29.44d0*exp(-(log10(Te_eV)-4.4283d0)**2.d0/(2.d0*(2.8428d0)**2.d0)) &
                                        -60.947d0*exp(-(log10(Te_eV)+2.0835d0)**2.d0/(2.d0*(0.9048d0)**2.d0)) &
                                        -24.067d0*exp(-(log10(Te_eV)+0.7363d0)**2.d0/(2.d0*(2.1700d0)**2.d0)))
    
      dLradDrays_dT = -coef_rad_1*(-29.440d0*(2.8428d0)**(-2.d0)*(log10(Te_eV)-4.4283d0)/T0*exp(-(log10(Te_eV)-4.4283d0)**2.d0/(2.d0*(2.8428d0)**2.d0)) &
                                    -60.947d0*(0.9048d0)**(-2.d0)*(log10(Te_eV)+2.0835d0)/T0*exp(-(log10(Te_eV)+2.0835d0)**2.d0/(2.d0*(0.9048d0)**2.d0)) &
                                    -24.067d0*(2.1700d0)**(-2.d0)*(log10(Te_eV)+0.7363d0)/T0*exp(-(log10(Te_eV)+0.7363d0)**2.d0/(2.d0*(2.1700d0)**2.d0)))&
                                    *(1.d1)**(-29.440d0*exp(-(log10(Te_eV)-4.4283d0)**2.d0/(2.d0*(2.8428d0)**2.d0)) &
                                    -60.947d0*exp(-(log10(Te_eV)+2.0835d0)**2.d0/(2.d0*(0.9048d0)**2.d0)) &
                                    -24.067d0*exp(-(log10(Te_eV)+0.7363d0)**2.d0/(2.d0*(2.1700d0)**2.d0))) 

    else  ! don't radiate at very low temperatures
  
      LradDcont_T   = 0.d0
      dLradDcont_dT = 0.d0
      LradDrays_T   = 0.d0
      dLradDrays_dT = 0.d0
      dSrec_dT      = 0.d0

    endif

    ! --- Transform derivatives on total T, to derivatives in Te 
    dSion_dT      = dSion_dT      * 2.d0
    dSrec_dT      = dSrec_dT      * 2.d0
    dLradDrays_dT = dLradDrays_dT * 2.d0
    dLradDcont_dT = dLradDcont_dT * 2.d0

  else  ! --- use OPEN ADAS

    Te_eV_lim  = max(Te_eV,     0.2d0)  ! ADAS data is given between 0.2eV and 10 keV  
    Te_eV_lim  = min(Te_eV_lim, 1.d4 )  

    Te_si_log10= log10( Te_eV_lim / K_BOLTZ * EL_CHG )

    ne_si      = 1.d20
    if (present(ne0) .and. (.not. deuterium_adas_1e20) ) then
      ne_si = ne0 * central_density * 1.d20
      ne_si = max(ne_si,  1.d14)    ! ADAS density is bewteen 1.d14 and 1.21 m^-3
      ne_si = min(ne_si,  1.d21) 
    endif

    ne_si_log10= log10(ne_si)

    call ad_deuterium%scd%interp( 0, ne_si_log10, Te_si_log10, Sion_T, dSion_dT)
    call ad_deuterium%acd%interp( 1, ne_si_log10, Te_si_log10, Srec_T, dSrec_dT)
    call ad_deuterium%prb%interp( 1, ne_si_log10, Te_si_log10, LradDcont_T, dLradDCont_dT)
    call ad_deuterium%plt%interp( 0, ne_si_log10, Te_si_log10, LradDrays_T, dLradDrays_dT)

    if ( Te_eV < 0.2d0) then  ! --- Don't radiate or ionize below 0.2 eV, recombination allowed
      LradDcont_T   = 0.d0
      dLradDcont_dT = 0.d0
      LradDrays_T   = 0.d0
      dLradDrays_dT = 0.d0
      dSrec_dT      = 0.d0
      Sion_T        = 0.d0
      dSion_dT      = 0.d0
    endif

    if ( Te_eV > 1.d4) then   ! --- Fix values beyond 10 keV and remove derivatives
      dLradDcont_dT = 0.d0
      dLradDrays_dT = 0.d0
      dSrec_dT      = 0.d0
      dSion_dT      = 0.d0
    endif

    ! --- Transform the coefficients to JOREK units
    Sion_T        = Sion_T   * t_norm * central_density * 1.d20
    Srec_T        = Srec_T   * t_norm * central_density * 1.d20
    LradDCont_T   = LradDCont_T   * (central_density * 1.d20)**2 * MU_zero * t_norm * gamma_factor 
    LradDrays_T   = LradDrays_T   * (central_density * 1.d20)**2 * MU_zero * t_norm * gamma_factor 

    dSion_dT      = dSion_dT * t_norm * central_density * 1.d20/ (K_BOLTZ*MU_ZERO*central_density*1.d20)
    dSrec_dT      = dSrec_dT * t_norm * central_density * 1.d20/ (K_BOLTZ*MU_ZERO*central_density*1.d20)
    dLradDCont_dT = dLradDCont_dT * (central_density * 1.d20)**2 * MU_zero * t_norm * gamma_factor &
                                  / (K_BOLTZ*MU_ZERO*central_density*1.d20)
    dLradDrays_dT = dLradDrays_dT * (central_density * 1.d20)**2 * MU_zero * t_norm * gamma_factor & 
                                  / (K_BOLTZ*MU_ZERO*central_density*1.d20) ! factor to get the T derivative in JOREK units

  endif

end subroutine


subroutine plot_atomic_coefficients()

  implicit none

  real*8 :: Te0                        
  real*8 :: Sion_T, dSion_dT           
  real*8 :: Srec_T, dSrec_dT           
  real*8 :: LradDcont_T, dLradDcont_dT 
  real*8 :: LradDrays_T, dLradDrays_dT 
  real*8 :: Tmin, Tmax, deltaT         
  integer:: i, nT

  deltaT = 2.d-6
  Tmin   = 1.d-6
  Tmax   = 1.d-1

  nT = int(Tmax/deltaT)
 
  open(unit=29, file='coefficients_atomic.txt', action='write')

  do i=1, nT

    Te0 = Tmin + (Tmax - Tmin)*float(i-1)/float(nT-1)

    call atomic_coeff_deuterium(Te0, Sion_T, dSion_dT, Srec_T, dSrec_dT,        &
                                LradDcont_T, dLradDcont_dT, LradDrays_T, dLradDrays_dT ) 

     write(29,'(9ES14.6)') Te0, Sion_T, dSion_dT, Srec_T, dSrec_dT, LradDcont_T, dLradDcont_dT, LradDrays_T, dLradDrays_dT

  enddo

  close(29)  


end subroutine






end module mod_atomic_coeff_deuterium


