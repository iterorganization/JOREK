!> Module containing plasma functions t
module mod_plasma_functions 
  
  use phys_module, only: eta_T_dependent, T_min, ZKpar_T_dependent, visco_T_dependent,  &
                         eta_num, eta_num_T_dependent, eta_num_psin_dependent, eta_num_prof, &
                         visco_num, visco_num_T_dependent, T_max_visco 
  use mod_model_settings, only:  with_impurities
    
  implicit none
  
  private
  public resistivity, conductivity_parallel, viscosity, hyper_resistivity, hyper_viscosity
  
  contains
  
 



 
  !> Determine resistivity (input/output in JOREK units)
  pure subroutine resistivity(eta_0, T_raw, T_corr, T_max, T0, Z_eff, eta_T,                       & 
                              dZ_eff_dT, dZ_eff_dr0, dZ_eff_drimp0, dr0_corr_dn, drimp0_corr_dn,   & 
                              deta_dT, d2eta_d2T, deta_dr0, deta_drimp0) 

    implicit none
    
    real*8, intent(in)             :: eta_0            ! central resistivity
    real*8, intent(in)             :: T_raw            ! temperature without correction
    real*8, intent(in)             :: T_corr           ! corrected temperature > 0
    real*8, intent(in)             :: T_max            ! max temperature to use in the function
    real*8, intent(in)             :: T0               ! central temperature at equilibrium
    real*8, intent(in)             :: Z_eff            ! effective charge (only used with_impurities at the moment)
    real*8, intent(out)            :: eta_T            ! output resistivity
    real*8, optional, intent(in)   :: dZ_eff_dT        ! Derivative of Zeff w.r.t. the temperature 
    real*8, optional, intent(in)   :: dZ_eff_dr0       ! Derivative of Zeff w.r.t. the total density 
    real*8, optional, intent(in)   :: dZ_eff_drimp0    ! Derivative of Zeff w.r.t. the impurity density
    real*8, optional, intent(in)   :: dr0_corr_dn      ! Derivative of density correction  
    real*8, optional, intent(in)   :: drimp0_corr_dn   ! Derivative of impurity density correction  
    real*8, optional, intent(out)  :: deta_dT          ! 1st derivative with respect to the temperature
    real*8, optional, intent(out)  :: d2eta_d2T        ! 2nd derivative with respect to the temperature
    real*8, optional, intent(out)  :: deta_dr0         ! 1st derivative with respect to the total density
    real*8, optional, intent(out)  :: deta_drimp0      ! 1st derivative with respect to the timpurity density

    !--- Local parameters
    real*8 :: eta_coef, deta_coef_dZeff

    if (present(deta_dT))       deta_dT     = 0.d0
    if (present(d2eta_d2T))     d2eta_d2T   = 0.d0
    if (present(deta_dr0))      deta_dr0    = 0.d0
    if (present(deta_drimp0))   deta_drimp0 = 0.d0

    if ( eta_T_dependent .and. T_corr <= T_max ) then
      eta_T     =   eta_0 * (T_corr/T0)**(-1.5d0)
      if (present(deta_dT))    deta_dT   = - eta_0 * (1.5d0)  * T_corr**(-2.5d0) * T0**(1.5d0)
      if (present(d2eta_d2T))  d2eta_d2T =   eta_0 * (3.75d0) * T_corr**(-3.5d0) * T0**(1.5d0)
    else if (eta_T_dependent .and. T_corr > T_max) then
      eta_T     = eta_0 * (T_max/T0)**(-1.5d0)
    else
      eta_T     = eta_0
    end if

    if ( eta_T_dependent .and. (T_raw .lt. T_min) ) then
      eta_T       = eta_0     * (T_min/T0)**(-1.5d0)
      if (present(deta_dT))     deta_dT    = 0.d0
      if (present(d2eta_d2T))   d2eta_d2T  = 0.d0
    endif

    if (with_impurities) then 

      eta_coef     = Z_eff*(1.+1.198*Z_eff+0.222*Z_eff**2)/(1.+2.966*Z_eff+0.753*Z_eff**2)
      eta_coef     = eta_coef / ((1.+1.198+0.222)/(1.+2.966+0.753))

      deta_coef_dZeff = (1.+1.198*Z_eff+0.222*Z_eff**2)/(1.+2.966*Z_eff+0.753*Z_eff**2)
      deta_coef_dZeff = deta_coef_dZeff + Z_eff*(1.198+2.*0.222*Z_eff)/(1.+2.966*Z_eff+0.753*Z_eff**2)
      deta_coef_dZeff = deta_coef_dZeff - Z_eff*(1.+1.198*Z_eff+0.222*Z_eff**2)*(2.966+2.*0.753*Z_eff)/((1.+2.966*Z_eff+0.753*Z_eff**2)**2)
      deta_coef_dZeff = deta_coef_dZeff / ((1.+1.198+0.222)/(1.+2.966+0.753))

      if ( eta_T_dependent ) then
        eta_T       = eta_T * eta_coef
        if (present(deta_dr0) .and. present(dZ_eff_dr0) .and. present(dr0_corr_dn)) then
          deta_dr0    = eta_T * deta_coef_dZeff * dZ_eff_dr0 * dr0_corr_dn
        endif
        if (present(deta_drimp0) .and. present(dZ_eff_drimp0) .and. present(drimp0_corr_dn)) then
            deta_drimp0 = eta_T * deta_coef_dZeff * dZ_eff_drimp0 * drimp0_corr_dn
        endif
        if (present(deta_dT) .and. present(dZ_eff_dT)) then
          deta_dT     = deta_dT * eta_coef + eta_T * deta_coef_dZeff * dZ_eff_dT
        endif
      end if

    endif

  end subroutine resistivity








  !> Determine hyper resistivity (input/output in JOREK units)
  pure subroutine hyper_resistivity(T_raw, T_corr, T0, psi_norm, eta_num_T, deta_num_dT) 

    implicit none
    
    real*8, intent(in)             :: T_raw            ! temperature without correction
    real*8, intent(in)             :: T_corr           ! corrected temperature > 0
    real*8, intent(in)             :: T0               ! central temperature at equilibrium
    real*8, intent(in)             :: psi_norm         ! normalized poloidal flux          
    real*8, intent(out)            :: eta_num_T        ! output hyper resistivity
    real*8, optional, intent(out)  :: deta_num_dT      ! Derivative w.r.t. the temperature 

    ! --- Hyper-resistivity
    if ( eta_num_psin_dependent ) then
      eta_num_T   = eta_num * 0.5d0 * ( 1.d0 - tanh( (psi_norm-eta_num_prof(1))/eta_num_prof(2)) )      
      if (present(deta_num_dT))  deta_num_dT = 0.d0      
    else if ( eta_num_T_dependent ) then
      eta_num_T     =   eta_num   * (T_corr/T0)**(-3.d0)
      if (present(deta_num_dT))  deta_num_dT   = - eta_num   * (3.d0)  * T_corr**(-4.d0) * T0**(3.d0)
      if (T_raw .lt. T_min) then
        eta_num_T     = eta_num    * (T_min/T0)**(-3.d0)
        if (present(deta_num_dT)) deta_num_dT   = 0.d0
      endif
    else
      eta_num_T     = eta_num
      if (present(deta_num_dT))   deta_num_dT   = 0.d0
    end if

  end subroutine hyper_resistivity






  ! --- Parallel conductivity (input-output in JOREK units)
  pure subroutine conductivity_parallel(ZK_par0, ZK_par_max, T_raw, T_corr, T_min_ZKpar, T0, ZK_par_T, &
                                         ZKpar_neg_thresh, ZK_par_neg, dT0_corr_dT, dZK_par_dT)

    real*8, intent(in)             :: ZK_par0          ! central parallel conduction
    real*8, intent(in)             :: ZK_par_max       ! maximum value for parallel conduction
    real*8, intent(in)             :: T_raw            ! temperature without correction
    real*8, intent(in)             :: T_corr           ! corrected temperature > 0
    real*8, intent(in)             :: T_min_ZKpar      ! min temperature to use in the function
    real*8, intent(in)             :: T0               ! central temperature at equilibrium
    real*8, intent(out)            :: ZK_par_T         ! output parallel conduction
    real*8, optional, intent(in)   :: ZKpar_neg_thresh ! threshold for negative correction             
    real*8, optional, intent(in)   :: ZK_par_neg       ! value after negative correction             
    real*8, optional, intent(in)   :: dT0_corr_dT      ! derivative of temperature correction
    real*8, optional, intent(out)  :: dZK_par_dT       ! temperature derviative of parallel conduction

    ! --- Local parameters
    real*8 :: dZK_par_dT_tmp

    if ( ZKpar_T_dependent ) then
      ZK_par_T       = ZK_par0 * (T_corr/T0)**(+2.5d0)   
      dZK_par_dT_tmp = ZK_par0 * (2.5d0)  * T_corr**(+1.5d0) * T0**(-2.5d0) 
      if (ZK_par_T .gt. ZK_par_max) then
        ZK_par_T       = Zk_par_max
        dZK_par_dT_tmp = 0.d0
      endif
      if (T_raw .lt. T_min_ZKpar) then
        ZK_par_T       = ZK_par0 * (T_min_ZKpar/T0)**(+2.5d0)
        dZK_par_dT_tmp = 0.d0
      endif
    else
      ZK_par_T      = ZK_par0 
      dZK_par_dT_tmp = 0.d0
    endif

    if (present(dZK_par_dT)) then
      dZK_par_dT = dZK_par_dT_tmp
      if (present(dT0_corr_dT)) dZK_par_dT = dZK_par_dT * dT0_corr_dT
    endif

    ! --- Increase value to avoid negative temperatures
    if (present(ZKpar_neg_thresh) .and. present(ZK_par_neg)) then
      if (T_raw .lt. ZKpar_neg_thresh) then
        ZK_par_T = ZK_par_neg
        if (present(dZK_par_dT)) dZK_par_dT = 0.d0
      endif
    endif

  end subroutine conductivity_parallel 






  !> Determine hyper viscosity (input/output in JOREK units)
  pure subroutine hyper_viscosity(T_raw, T_corr, T0, visco_num_T, dvisco_num_dT) 

    implicit none
    
    real*8, intent(in)             :: T_raw            ! temperature without correction
    real*8, intent(in)             :: T_corr           ! corrected temperature > 0
    real*8, intent(in)             :: T0               ! central temperature at equilibrium
    real*8, intent(out)            :: visco_num_T        ! output viscosity
    real*8, optional, intent(out)  :: dvisco_num_dT      ! Derivative of Zeff w.r.t. the temperature 

    if ( visco_num_T_dependent ) then
      visco_num_T     =   visco_num   * (T_corr/T0)**(-3.d0)
      if (present(dvisco_num_dT))  dvisco_num_dT   = - visco_num   * (3.d0)  * T_corr**(-4.d0) * T0**(3.d0)
      if (T_raw .lt. T_min) then
        visco_num_T     = visco_num    * (T_min/T0)**(-3.d0)
        if (present(dvisco_num_dT)) dvisco_num_dT   = 0.d0
      endif
    else
      visco_num_T     = visco_num
      if (present(dvisco_num_dT))   dvisco_num_dT   = 0.d0
    end if

  end subroutine hyper_viscosity









  ! --- Viscosity (input-output in JOREK units)
  pure subroutine viscosity(visco, T_raw, T_corr,T0, visco_T, dvisco_dT, d2visco_d2T)

    real*8, intent(in)             :: visco            ! temperature without correction
    real*8, intent(in)             :: T_raw            ! temperature without correction
    real*8, intent(in)             :: T_corr           ! corrected temperature > 0
    real*8, intent(in)             :: T0               ! central temperature at equilibrium
    real*8, intent(out)            :: visco_T          ! output viscosity          
    real*8, optional, intent(out)  :: dvisco_dT        ! 1st derivative w.r.t. temperature   
    real*8, optional, intent(out)  :: d2visco_d2T      ! 2nd derivative w.r.t. temperature           

    if ( visco_T_dependent ) then
      visco_T     =   visco * (T_corr/T0)**(-1.5d0)
      if (present(dvisco_dT))   dvisco_dT   = - visco * (1.5d0)  * T_corr**(-2.5d0) * T0**(1.5d0)
      if (present(d2visco_d2T)) d2visco_d2T =   visco * (3.75d0) * T_corr**(-3.5d0) * T0**(1.5d0)
      if (T_raw .lt. T_min) then
        visco_T     = visco  * (T_min/T0)**(-1.5d0)
        if (present(dvisco_dT  )) dvisco_dT   = 0.d0
        if (present(d2visco_d2T)) d2visco_d2T = 0.d0
      else if (T_raw .gt. T_max_visco) then
        visco_T     =   visco * (T_max_visco/T0)**(-1.5d0)
        if (present(dvisco_dT  )) dvisco_dT   = 0.d0
        if (present(d2visco_d2T)) d2visco_d2T = 0.d0
      endif
    else
      visco_T     = visco
      if (present(dvisco_dT))   dvisco_dT   = 0.d0
      if (present(d2visco_d2T)) d2visco_d2T = 0.d0
    end if

  end subroutine viscosity 







end module mod_plasma_functions
