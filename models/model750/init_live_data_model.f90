!> Model-specific part of the routine init_live_data in communication/mod_live_data.f90
!!
!! Writes all input profiles to 'macroscopic_vars.dat'.
subroutine init_live_data_model(file_handle)

  use mod_parameters        
  use phys_module,   only: xpoint, xcase
  use diffusivities, only: get_dperp, get_zk_iperp, get_zk_eperp
  use mod_sources

  implicit none

  integer, intent(in) :: file_handle
  
  integer :: i
  real*8  :: psin, FFp, dFFp_dpsi, dens, dn_dpsi
  real*8  :: temp, temp_i, dT_dpsi, dTi_dpsi, temp_e, dTe_dpsi, S_rho, S_T, S_Ti, S_Te, d_perp, zk_perp, zki_perp, zke_perp, rhon, drhon_dpsi
  real*8  :: d, d1, d2, d3, d4, d5, d6 ! dummies
  
  write(file_handle,'(A,I5)') '@n_input_profiles: ', 10
  write(file_handle,'(A)') '@input_profiles_xlabel: Psi_{normalized}'
  write(file_handle,'(A)') '@input_profiles_xlabel_si: Psi_{normalized}'
  write(file_handle,'(A)') '@input_profiles_ylabel: input profiles'
  write(file_handle,'(A)') '@input_profiles_ylabel_si: input profiles'
  write(file_handle,'(A,5ES17.9)') '@input_profiles_x2si: ', 1.0
  write(file_handle,'(A,5ES17.9)') '@input_profiles_y2si: ', 1.0
  write(file_handle,'(A)') '@input_profiles_logy: 0'
  if(with_TiTe)then
    write(file_handle,'(A)') '@input_profiles: %"psin"       "FF''"    "dFF''/dpsin"'             // &
      '    "rho"    "drho/dpsin"   "Ti"     "dTi/dpsin"   "Te"     "dTe/dpsin"   "S_rho"     "S_Ti"     "S_Te"    "rhon"         "D_perp"'      // &
      '    "ZKi_perp"    "ZKe_perp"'
  else
    write(file_handle,'(A)') '@input_profiles: %"psin"       "FF''" "dFF''/dpsin"'             // &
      '    "rho"    "drho/dpsin"   "T"     "dT/dpsin"  "S_rho"     "S_T"    "rhon"         "D_perp"'      // &
      '    "ZK_perp"'
  endif
  
  do i = 0, 200
    
    psin = real(i) / real(200) * 1.2d0
    
    call density        (xpoint,xcase,0.d0,(/-99.d0,-99.d0/),psin,0.d0,1.d0,dens,dn_dpsi,d,d1,d2,d3,d4,d5,d6)
    if(with_TiTe)then
      call temperature_i  (xpoint,xcase,0.d0,(/-99.d0,-99.d0/),psin,0.d0,1.d0,temp_i,dTi_dpsi,d,d1,d2,d3,d4,d5,d6)
      call temperature_e  (xpoint,xcase,0.d0,(/-99.d0,-99.d0/),psin,0.d0,1.d0,temp_e,dTe_dpsi,d,d1,d2,d3,d4,d5,d6)
    else
      call temperature    (xpoint,xcase,0.d0,(/-99.d0,-99.d0/),psin,0.d0,1.d0,temp,dT_dpsi,d,d1,d2,d3,d4,d5,d6)
    endif    
    call sources        (xpoint,xcase,0.d0,(/-99.d0,-99.d0/),psin,0.d0,1.d0,S_rho,S_Ti,S_Te)
    call FFprime        (xpoint,xcase,0.d0,(/-99.d0,-99.d0/),psin,0.d0,1.d0,FFp,dFFp_dpsi,d,d1,d2,d3, .true.)
    if(with_neutrals)then
      call neutral_density(xpoint,xcase,0.d0,(/-99.d0,+99.d0/),psin,0.d0,1.d0,rhon,drhon_dpsi,d,d1,d2,d3,d4,d5,d6)
    endif    
    d_perp  = get_dperp (psin)
    zki_perp = get_zk_iperp(psin)
    zke_perp = get_zk_eperp(psin)
    
    if(with_TiTe)then
      write(file_handle,'(a,20es13.4e3)') '@input_profiles: ', psin, FFp, dFFp_dpsi, dens, dn_dpsi,    &
        temp_i, dTi_dpsi, temp_e, dTe_dpsi, S_rho, S_Ti, S_Te, rhon, d_perp, zki_perp, zke_perp
    else
      write(file_handle,'(a,17es13.4e3)') '@input_profiles: ', psin, FFp, dFFp_dpsi, dens, dn_dpsi,    &
        temp, dT_dpsi, S_rho, S_T, rhon, d_perp, zk_perp
    endif
    
  end do
  write(file_handle,*)
  
end subroutine init_live_data_model
