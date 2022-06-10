!> Model-specific part of the routine init_live_data in communication/mod_live_data.f90
!!
!! Writes all input profiles to 'macroscopic_vars.dat'.
subroutine init_live_data_model(file_handle)
  
  use phys_module,   only: xpoint, xcase, D_perp, ZK_perp
  use mod_sources
  
  implicit none
  
  integer, intent(in) :: file_handle
  
  integer :: i
  real*8  :: psin, S_rho, S_T
  real*8  :: d, d1, d2, d3, d4, d5, d6 ! dummies
  
  write(file_handle,'(A,I5)') '@n_input_profiles: ', 10
  write(file_handle,'(A)') '@input_profiles_xlabel: Psi_{normalized}'
  write(file_handle,'(A)') '@input_profiles_xlabel_si: Psi_{normalized}'
  write(file_handle,'(A)') '@input_profiles_ylabel: input profiles'
  write(file_handle,'(A)') '@input_profiles_ylabel_si: input profiles'
  write(file_handle,'(A,5ES17.9)') '@input_profiles_x2si: ', 1.0
  write(file_handle,'(A,5ES17.9)') '@input_profiles_y2si: ', 1.0
  write(file_handle,'(A)') '@input_profiles_logy: 0'
  write(file_handle,'(A)') '@input_profiles: %"psin"       "S_rho"     "S_T"        "D_perp"'      // &
    '    "ZK_perp"'
  
  do i = 0, 200
    
    psin = real(i) / real(200) * 1.2d0
    
    write(file_handle,'(a,20es13.4e3)') '@input_profiles: ', psin, S_rho, S_T, d_perp(1), zk_perp(1)
    
  end do
  write(file_handle,*)
  
end subroutine init_live_data_model
