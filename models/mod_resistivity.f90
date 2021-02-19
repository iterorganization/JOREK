!> Module containing functions to determine the plasma resistivity 
module mod_resistivity 
  
  use phys_module, only: eta_T_dependent, T_0, Te_0, T_min, xpoint, eta, eta_ohmic 
  use mod_parameters, only: with_TiTe
    
  implicit none
  
  private
  public resistivity, dresistivity_dT 
  
 
  contains
  
  
  !> Determine resistivity
  pure function resistivity(eta_0,T,T_max) result(eta_T)
#if _OPENMP >= 201511
    !$omp declare simd
#endif
    implicit none
    
    real*8, intent(in)           :: T
    real*8, intent(in)           :: eta_0
    real*8, intent(in)           :: T_max
    real*8                       :: eta_T
    real*8                       :: T_local

    T_local     = max(T, T_min)

    ! --- Temperature dependent resistivity
    if (with_TiTe) then
      if ( eta_T_dependent .and. (T_local <= T_max)) then
        eta_T     = eta_0 * (T_local/Te_0)**(-1.5d0)
      else if ( eta_T_dependent .and. (T_local > T_max)) then
        eta_T     = eta_0 * (T_max/Te_0)**(-1.5d0)
      else
        eta_T     = eta_0
      end if
    else
      if ( eta_T_dependent .and. (T_local <= T_max)) then
        eta_T     = eta_0 * (T_local/T_0)**(-1.5d0)
      else if ( eta_T_dependent .and. (T_local > T_max)) then
        eta_T     = eta_0 * (T_max/T_0)**(-1.5d0)
      else
        eta_T     = eta_0
      end if
    end if

  end function resistivity
    
  !> Determine eta derivative 
  pure function dresistivity_dT(eta_0,T,T_max) result(deta_dT)
#if _OPENMP >= 201511
    !$omp declare simd
#endif
    implicit none
    
    real*8, intent(in)           :: T
    real*8, intent(in)           :: eta_0
    real*8, intent(in)           :: T_max
    real*8                       :: deta_dT

    ! --- Temperature dependent resistivity
    if (with_TiTe) then
      if ( eta_T_dependent .and. (T <= T_max)) then
        deta_dT   = - eta   * (1.5d0)  * T**(-2.5d0) * Te_0**(1.5d0)
        if ( xpoint .and. (T .lt. T_min) ) then
          deta_dT   = 0.d0
        endif
      else
        deta_dT   = 0.
      end if
    else
      if ( eta_T_dependent .and. (T <= T_max)) then
        deta_dT   = - eta   * (1.5d0)  * T**(-2.5d0) * T_0**(1.5d0)
        if ( xpoint .and. (T .lt. T_min) ) then
          deta_dT   = 0.d0
        endif
      else
        deta_dT   = 0.
      end if
    endif
   
  end function dresistivity_dT 
  


end module mod_resistivity 
