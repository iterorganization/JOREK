!> NUMERICAL IMPROVEMENT FOR CASES WHERE TEMPERATURES CLOSE TO OR BELOW ZERO CAN OCCUR.
!! 
!! PROBLEM:
!! 
!! T^-1.5 is undefined for negative temperatures (resistivity and other quantities).
!! Thus, abs(T) was used so far which may, however, still cause discontinuities
!! in the resistivity.
!! 
!! SOLUTION:
!! 
!! Replace abs(T) by smooth function:
!! 
!!   f(T) = T                             if T>L1+L2
!!   f(T) = L1 + L2 * exp((T-(L2+L1))/L2) otherwise
!!   
!! where L1 and L2 are derived from the input parameter T_1:
!! 
!!   L1 = T_1 * corr_neg_temp_coef(1)
!!   L2 = T_1 * corr_neg_temp_coef(2)
!!   
!! The default values corr_neg_temp_coef(:) = (/ 0.5, 0.5 /) can be
!! changed via the namelist input file. Alternatively, different values
!! can be provided via the optional routine parameter coef.
!!

real*8 function corr_neg_temp1(val)
#if _OPENMP >= 201511
  !$omp declare simd 
#endif
  use phys_module, only: T_1, T_min_neg, corr_neg_temp_coef
  
  ! --- Routine parameters
  real*8, intent(in)           :: val       !< Temperature value to be "corrected".
  
  real*8 :: L1, L2
  
  if (T_min_neg .ge. 0.d0) then
	L1 = T_min_neg * corr_neg_temp_coef(1)  
	L2 = T_min_neg * corr_neg_temp_coef(2) 
  elseif (T_min_neg .lt. 0.d0) then
	L1 = T_1 * corr_neg_temp_coef(1)  
	L2 = T_1 * corr_neg_temp_coef(2) 
  endif

  corr_neg_temp1 = val
  if ( val < L1 + L2 ) corr_neg_temp1 = L1 + L2 * exp( (val-(L1+L2)) / L2 )

end function corr_neg_temp1

real*8 function corr_neg_temp2(val, coef)
! With uniform, we declare thet coeff should be the same for all vector elements.
! Is this correct?
  use phys_module, only: T_1, T_min_neg

  ! --- Routine parameters
  real*8, intent(in)           :: val       !< Temperature value to be "corrected".
  real*8, intent(in) :: coef(2)   !< Optional coefficients, if not provided the
                                            !! input parameter corr_neg_temp2_coef is used instead.  
  real*8 :: L1, L2

#if _OPENMP >= 201511
  !$omp declare simd uniform(coef)
#endif
  
  if (T_min_neg .ge. 0.d0) then
	L1 = T_min_neg * coef(1)
	L2 = T_min_neg * coef(2)
  elseif (T_min_neg .lt. 0.d0) then
	L1 = T_1 * coef(1)
	L2 = T_1 * coef(2)
  endif
  
  corr_neg_temp2 = val
  if ( val < L1 + L2 ) corr_neg_temp2 = L1 + L2 * exp( (val-(L1+L2)) / L2 )

end function corr_neg_temp2

real*8 function corr_neg_temp3(val, coef, val_1)
! With uniform, we declare thet coeff should be the same for all vector
! elements.
! Is this correct?
  use phys_module, only: T_1,T_min_neg

  ! --- Routine parameters
  real*8, intent(in)           :: val       !< Temperature value to be "corrected".
  real*8, intent(in)           :: coef(2)   !< Optional coefficients
  real*8, intent(in)           :: val_1     !< Temperature value floor
  real*8 :: L1, L2

#if _OPENMP >= 201511
  !$omp declare simd uniform(coef)
#endif

  L1 = val_1 * coef(1)
  L2 = val_1 * coef(2)

  corr_neg_temp3 = val
  if ( val < L1 + L2 ) corr_neg_temp3 = L1 + L2 * exp( (val-(L1+L2)) / L2 )

end function corr_neg_temp3

real*8 function dcorr_neg_temp_dT1(val)
#if _OPENMP >= 201511
  !$omp declare simd 
#endif
  use phys_module, only: T_1, corr_neg_temp_coef,T_min_neg
  
  ! --- Routine parameters
  real*8, intent(in)           :: val       !< Temperature value to be "corrected".
  
  real*8 :: L1, L2
  
  if (T_min_neg .ge. 0.d0) then
	L1 = T_min_neg * corr_neg_temp_coef(1) 
	L2 = T_min_neg * corr_neg_temp_coef(2) 
  elseif (T_min_neg .lt. 0.d0) then
	L1 = T_1 * corr_neg_temp_coef(1)
	L2 = T_1 * corr_neg_temp_coef(2)
  endif

  dcorr_neg_temp_dT1 = 1.
  if ( val < L1 + L2 ) dcorr_neg_temp_dT1 = exp( (val-(L1+L2)) / L2 )

end function dcorr_neg_temp_dT1

real*8 function dcorr_neg_temp_dT2(val, coef)
! With uniform, we declare thet coeff should be the same for all vector elements.
! Is this correct?
  use phys_module, only: T_min_neg, T_1

  ! --- Routine parameters
  real*8, intent(in)           :: val       !< Temperature value to be "corrected".
  real*8, intent(in) :: coef(2)   !< Optional coefficients, if not provided the
                                            !! input parameter corr_neg_temp2_coef is used instead.  
  real*8 :: L1, L2

#if _OPENMP >= 201511
  !$omp declare simd uniform(coef)
#endif 
 
  if (T_min_neg .ge. 0.d0) then
	L1 = T_min_neg * coef(1)
	L2 = T_min_neg * coef(2)
  elseif (T_min_neg .lt. 0.d0) then
	L1 = T_1 * coef(1)
	L2 = T_1 * coef(2)
  endif

  dcorr_neg_temp_dT2 = 1.
  if ( val < L1 + L2 ) dcorr_neg_temp_dT2 = exp( (val-(L1+L2)) / L2 )

end function dcorr_neg_temp_dT2

real*8 function dcorr_neg_temp_dT3(val, coef, val_1)
! With uniform, we declare thet coeff should be the same for all vector
! elements.
! Is this correct?
  use phys_module, only: T_1

  ! --- Routine parameters
  real*8, intent(in)           :: val       !< Temperature value to be "corrected".
  real*8, intent(in)           :: coef(2)   !< Optional coefficients
  real*8, intent(in)           :: val_1     !< Temperature value floor
  real*8 :: L1, L2

#if _OPENMP >= 201511
  !$omp declare simd uniform(coef)
#endif
  
  L1 = val_1 * coef(1)
  L2 = val_1 * coef(2)

  dcorr_neg_temp_dT3 = 1.
  if ( val < L1 + L2 ) dcorr_neg_temp_dT3 = exp( (val-(L1+L2)) / L2 )

end function dcorr_neg_temp_dT3

real*8 function d2corr_neg_temp_dT21(val)
#if _OPENMP >= 201511
  !$omp declare simd 
#endif
  use phys_module, only: corr_neg_temp_coef,T_min_neg,T_1
  
  ! --- Routine parameters
  real*8, intent(in)           :: val       !< Temperature value to be "corrected".
  
  real*8 :: L1, L2
  
  if (T_min_neg .ge. 0.d0) then
	L1 = T_min_neg * corr_neg_temp_coef(1) 
	L2 = T_min_neg * corr_neg_temp_coef(2) 
  elseif (T_min_neg .lt. 0.d0) then
	L1 = T_1 * corr_neg_temp_coef(1)
	L2 = T_1 * corr_neg_temp_coef(2)
  endif


  d2corr_neg_temp_dT21 = 0.
  if ( val < L1 + L2 ) d2corr_neg_temp_dT21 = exp( (val-(L1+L2)) / L2 ) / L2

end function d2corr_neg_temp_dT21

real*8 function d2corr_neg_temp_dT22(val, coef)
! With uniform, we declare thet coeff should be the same for all vector elements.
! Is this correct?
  use phys_module, only: T_min_neg, T_1

  ! --- Routine parameters
  real*8, intent(in)           :: val       !< Temperature value to be "corrected".
  real*8, intent(in) :: coef(2)   !< Optional coefficients, if not provided the
                                            !! input parameter corr_neg_temp2_coef is used instead.
  
  real*8 :: L1, L2

#if _OPENMP >= 201511
  !$omp declare simd uniform(coef)
#endif
  
  if (T_min_neg .ge. 0.d0) then
	L1 = T_min_neg * coef(1)
	L2 = T_min_neg * coef(2)
  elseif (T_min_neg .lt. 0.d0) then
	L1 = T_1 * coef(1)
	L2 = T_1 * coef(2)
  endif
  

  d2corr_neg_temp_dT22 = 0.
  if ( val < L1 + L2 ) d2corr_neg_temp_dT22 = exp( (val-(L1+L2)) / L2 ) / L2

end function d2corr_neg_temp_dT22

real*8 function d2corr_neg_temp_dT23(val, coef, val_1)
! With uniform, we declare thet coeff should be the same for all vector
! elements.
! Is this correct?
  use phys_module, only: T_1

  ! --- Routine parameters
  real*8, intent(in)           :: val       !< Temperature value to be "corrected".
  real*8, intent(in)           :: coef(2)   !< Optional coefficients
  real*8, intent(in)           :: val_1     !< Temperature value floor
  real*8 :: L1, L2

#if _OPENMP >= 201511
  !$omp declare simd uniform(coef)
#endif
  
  L1 = val_1 * coef(1)
  L2 = val_1 * coef(2)

  d2corr_neg_temp_dT23 = 0.
  if ( val < L1 + L2 ) d2corr_neg_temp_dT23 = exp( (val-(L1+L2)) / L2 ) / L2

end function d2corr_neg_temp_dT23

!> Same for density (so far not used in element_matrix routines).
real*8 function corr_neg_dens1(val)
#if _OPENMP >= 201511
!$omp declare simd
#endif
  use phys_module, only: corr_neg_dens_coef, rho_min_neg, rho_1
  
  ! --- Routine parameters
  real*8, intent(in)           :: val       !< Density value to be "corrected".
  
  real*8 :: L1, L2

  if (rho_min_neg .ge. 0.d0) then
	L1 = rho_min_neg * corr_neg_dens_coef(1)
	L2 = rho_min_neg * corr_neg_dens_coef(2)
  elseif (rho_min_neg .lt. 0.d0) then
	L1 = rho_1 * corr_neg_dens_coef(1)
	L2 = rho_1 * corr_neg_dens_coef(2)
  endif
  

  corr_neg_dens1 = val
  if ( val < L1 + L2 ) corr_neg_dens1 = L1 + L2 * exp( (val-(L1+L2)) / L2 )

end function corr_neg_dens1

!We cannot have optional argument for a vector funct, therefore we overload it
real*8 function corr_neg_dens2(val, coef)
! With uniform, we declare thet coeff should be the same for all vector elements.
! Is this correct?
  use phys_module, only: corr_neg_dens_coef, rho_min_neg, rho_1
  
  ! --- Routine parameters
  real*8, intent(in)           :: val       !< Density value to be "corrected".
  real*8, intent(in) :: coef(2)   !< Optional coefficients, if not provided the
                                            !! input parameter corr_neg_temp_coef is used instead.
  real*8 :: L1, L2

#if _OPENMP >= 201511
  !$omp declare simd uniform(coef)
#endif
  
  if (rho_min_neg .ge. 0.d0) then
	L1 = rho_min_neg * coef(1)
	L2 = rho_min_neg * coef(2)
  elseif (rho_min_neg .lt. 0.d0) then
	L1 = rho_1 * coef(1)
	L2 = rho_1 * coef(2)
  endif  
  

  corr_neg_dens2 = val
  if ( val < L1 + L2 ) corr_neg_dens2 = L1 + L2 * exp( (val-(L1+L2)) / L2 )

end function corr_neg_dens2

real*8 function corr_neg_dens3(val, coef, val_1)
! With uniform, we declare thet coeff should be the same for all vector
! elements.
! Is this correct?
  use phys_module, only: corr_neg_dens_coef, rho_min_neg, rho_1

  ! --- Routine parameters
  real*8, intent(in)           :: val       !< Temperature value to be "corrected".
  real*8, intent(in)           :: coef(2)   !< Optional coefficients
  real*8, intent(in)           :: val_1     !< Density value floor
  real*8 :: L1, L2

#if _OPENMP >= 201511
  !$omp declare simd uniform(coef)
#endif 
  
  L1 = val_1 * coef(1)
  L2 = val_1 * coef(2)

  corr_neg_dens3 = val
  if ( val < L1 + L2 ) corr_neg_dens3 = L1 + L2 * exp( (val-(L1+L2)) / L2 )

end function corr_neg_dens3


real*8 function dcorr_neg_dens_drho1(val)
#if _OPENMP >= 201511
!$omp declare simd
#endif
  use phys_module, only: corr_neg_dens_coef, rho_min_neg, rho_1
  
  ! --- Routine parameters
  real*8, intent(in)           :: val       !< Density value to be "corrected".
  
  real*8 :: L1, L2
  
  if (rho_min_neg .ge. 0.d0) then
	L1 = rho_min_neg * corr_neg_dens_coef(1)
	L2 = rho_min_neg * corr_neg_dens_coef(2)
  elseif (rho_min_neg .lt. 0.d0) then
	L1 = rho_1 * corr_neg_dens_coef(1)
	L2 = rho_1 * corr_neg_dens_coef(2)
  endif

  dcorr_neg_dens_drho1 = 1.
  if ( val < L1 + L2 ) dcorr_neg_dens_drho1 = exp( (val-(L1+L2)) / L2 )

end function dcorr_neg_dens_drho1

!We cannot have optional argument for a vector funct, therefore we overload it
real*8 function dcorr_neg_dens_drho2(val, coef)
! With uniform, we declare thet coeff should be the same for all vector elements.
! Is this correct?
  use phys_module, only: corr_neg_dens_coef, rho_min_neg, rho_1
  
  ! --- Routine parameters
  real*8, intent(in)           :: val       !< Density value to be "corrected".
  real*8, intent(in) :: coef(2)   !< Optional coefficients, if not provided the
                                            !! input parameter corr_neg_temp_coef is used instead.
  real*8 :: L1, L2

#if _OPENMP >= 201511
  !$omp declare simd uniform(coef)
#endif
  
  if (rho_min_neg .ge. 0.d0) then
	L1 = rho_min_neg * coef(1)
	L2 = rho_min_neg * coef(2)
  elseif (rho_min_neg .lt. 0.d0) then
	L1 = rho_1 * coef(1)
	L2 = rho_1 * coef(2)
  endif  


  dcorr_neg_dens_drho2 = 1.
  if ( val < L1 + L2 ) dcorr_neg_dens_drho2 = exp( (val-(L1+L2)) / L2 )

end function dcorr_neg_dens_drho2

real*8 function dcorr_neg_dens_drho3(val, coef, val_1)
! With uniform, we declare thet coeff should be the same for all vector
! elements.
! Is this correct?
  use phys_module, only: rho_1, corr_neg_dens_coef

  ! --- Routine parameters
  real*8, intent(in)           :: val       !< Temperature value to be "corrected".
  real*8, intent(in)           :: coef(2)   !< Optional coefficients
  real*8, intent(in)           :: val_1     !< Density value floor
  real*8 :: L1, L2

#if _OPENMP >= 201511
  !$omp declare simd uniform(coef)
#endif

  L1 = val_1 * coef(1)
  L2 = val_1 * coef(2)

  dcorr_neg_dens_drho3 = 1.
  if ( val < L1 + L2 ) dcorr_neg_dens_drho3 = exp( (val-(L1+L2)) / L2 )

end function dcorr_neg_dens_drho3