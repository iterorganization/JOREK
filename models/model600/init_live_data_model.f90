!> Model-specific part of the routine init_live_data in communication/mod_live_data.f90
!!
!! Writes all input profiles to 'macroscopic_vars.dat'.
subroutine init_live_data_model(file_handle)
  
  use mod_parameters
  use phys_module,   only: xpoint, xcase
  use diffusivities, only: get_dperp, get_zkperp, get_zk_iperp, get_zk_eperp
  use mod_sources
  
  implicit none
  
  integer, intent(in) :: file_handle
  
  integer :: i
  real*8  :: psin, FFp, dFFp_dpsi, dens, dn_dpsi, T, dT_dpsi, Te, dTe_dpsi, Ti, dTi_dpsi
  real*8  :: S_rho, S_T, S_Ti, S_Te, d_perp, zk_perp, zke_perp, zki_perp
  real*8  :: d, d1, d2, d3, d4, d5, d6 ! dummies
  
  write(file_handle,'(A,I5)') '@n_input_profiles: ', 13
  write(file_handle,'(A)') '@input_profiles_xlabel: Psi_{normalized}'
  write(file_handle,'(A)') '@input_profiles_xlabel_si: Psi_{normalized}'
  write(file_handle,'(A)') '@input_profiles_ylabel: input profiles'
  write(file_handle,'(A)') '@input_profiles_ylabel_si: input profiles'
  write(file_handle,'(A,5ES17.9)') '@input_profiles_x2si: ', 1.0
  write(file_handle,'(A,5ES17.9)') '@input_profiles_y2si: ', 1.0
  write(file_handle,'(A)') '@input_profiles_logy: 0'
  if ( with_TiTe ) then
    write(file_handle,'(A)') '@input_profiles: %"psin"       "FF''"    "dFF''/dpsin"'             // &
      '    "rho"    "drho/dpsin"   "Te"     "dTe/dpsin"   "Ti"     "dTi/dpsin"   "S_rho"'         // &
      '     "S_Te"       "S_Ti"       "D_perp"    "ZKe_perp"   "ZKi_perp"'
  else
    write(file_handle,'(A)') '@input_profiles: %"psin"       "FF''"    "dFF''/dpsin"'             // &
      '    "rho"    "drho/dpsin"   "T"     "dT/dpsin"   "S_rho"'                                  // &
      '     "S_T"       "D_perp"    "ZK_perp"'
  end if
  
  do i = 0, 200
    
    psin = real(i) / real(200) * 1.2d0
    
    call density      (xpoint,xcase,0.d0,(/-99.d0,+99.d0/),psin,0.d0,1.d0,dens,dn_dpsi,d,d1,d2,d3,d4,d5,d6)
    if ( with_TiTe ) then
      call temperature_e(xpoint,xcase,0.d0,(/-99.d0,+99.d0/),psin,0.d0,1.d0,Te,dTe_dpsi,d,d1,d2,d3,d4,d5,d6)
      call temperature_i(xpoint,xcase,0.d0,(/-99.d0,+99.d0/),psin,0.d0,1.d0,Ti,dTi_dpsi,d,d1,d2,d3,d4,d5,d6)
      call sources      (xpoint,xcase,0.d0,(/-99.d0,+99.d0/),psin,0.d0,1.d0,S_rho,S_Ti,S_Te)
      zke_perp = get_zk_eperp(psin)
      zki_perp = get_zk_iperp(psin)
    else
      call temperature(xpoint,xcase,0.d0,(/-99.d0,+99.d0/),psin,0.d0,1.d0,T,dT_dpsi,d,d1,d2,d3,d4,d5,d6)
      call sources    (xpoint,xcase,0.d0,(/-99.d0,+99.d0/),psin,0.d0,1.d0,S_rho,S_T)
      zk_perp = get_zkperp(psin)
    end if
    call FFprime      (xpoint,xcase,0.d0,(/-99.d0,+99.d0/),psin,0.d0,1.d0,FFp,dFFp_dpsi,d,d1,d2,d3)
    d_perp  = get_dperp (psin)
    
    if ( with_TiTe ) then
      write(file_handle,'(a,20es13.4e3)') '@input_profiles: ', psin, FFp, dFFp_dpsi, dens, dn_dpsi,    &
        Te, dTe_dpsi, Ti, dTi_dpsi, S_rho, S_Te, S_Ti, d_perp, zke_perp, zki_perp
    else
      write(file_handle,'(a,20es13.4e3)') '@input_profiles: ', psin, FFp, dFFp_dpsi, dens, dn_dpsi,    &
        T, dT_dpsi, S_rho, S_T, d_perp, zk_perp
    end if
      
  end do
  write(file_handle,*)
  
end subroutine init_live_data_model
