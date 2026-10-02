!> Module containing functions to determine particle and heat diffusivities
module diffusivities

  use constants  
  use mod_parameters,  only: jorek_model
  use phys_module, only: num_d_perp,  D_perp,  D_perp_pri,  D_perp_pri_z0, D_perp_pri_width,                        &
                         num_d_perp_x, num_d_perp_y, num_d_perp_len,                                                &
                         num_zk_perp, ZK_perp, ZK_perp_pri, ZK_perp_pri_z0, ZK_perp_pri_width,                      &
                         num_zk_perp_x, num_zk_perp_y, num_zk_perp_len,                                             &
                         num_zk_i_perp, ZK_i_perp, ZK_i_perp_pri, ZK_i_perp_pri_z0, ZK_i_perp_pri_width,            &
                         num_zk_i_perp_x, num_zk_i_perp_y, num_zk_i_perp_len,                                       &
                         num_zk_e_perp, ZK_e_perp, ZK_e_perp_pri, ZK_e_perp_pri_z0, ZK_e_perp_pri_width,            &
                         num_zk_e_perp_x, num_zk_e_perp_y, num_zk_e_perp_len,                                       &
                         num_v_pinch, V_pinch_gauss, V_pinch_psin, V_pinch_sig,                                     &
                         num_v_pinch_x, num_v_pinch_y, num_v_pinch_len,                                             &
                         xpoint, xcase, rho_0, rho_coef, T_coef, Ti_coef, Te_coef
  use profiles,    only: interpolProf
    
  implicit none
  
  
  private
  public get_dperp, get_zkperp, get_zk_iperp, get_zk_eperp, get_vpinch
  
  
  interface get_dperp
    module procedure get_dperp1
    module procedure get_dperp2
    module procedure get_dperp3
    module procedure get_dperp4
    module procedure get_dperp5
  end interface get_dperp
  
  interface get_zkperp
    module procedure get_zkperp1
    module procedure get_zkperp2
    module procedure get_zkperp3
  end interface get_zkperp
  
  interface get_zk_iperp
    module procedure get_zk_iperp1
    module procedure get_zk_iperp2
    module procedure get_zk_iperp3
  end interface get_zk_iperp
    
  interface get_zk_eperp
    module procedure get_zk_eperp1
    module procedure get_zk_eperp2
    module procedure get_zk_eperp3
  end interface get_zk_eperp
  
  contains
  
  !> Determine perpendicular particle diffusivity, D_perp, as a function of Psi_N
  real*8 function get_dperp1(psin)
#if _OPENMP >= 201511
    !$omp declare simd
#endif
    implicit none
    
    real*8, intent(in) :: psin
    
    if ( num_d_perp ) then
      
      get_dperp1 = interpolProf(num_d_perp_x, num_d_perp_y, num_d_perp_len, psin)
      
    else
      
      get_dperp1 = D_perp(1) * ( (1.d0-D_perp(2)) +                                                 &
        D_perp(2)*(0.5d0 - 0.5d0*tanh((psin-D_perp(5))/D_perp(4))) )
      
      if ( jorek_model >= 300 ) then
        
        get_dperp1 = get_dperp1 + D_perp(6)*D_perp(2) *                                              &
          ((0.5d0 - 0.5d0*tanh((-psin+D_perp(5)+D_perp(3)) /D_perp(4))))
        
      end if
      
    end if
    
  end function get_dperp1  
  
  !> Determine perpendicular particle diffusivity, D_perp, as a function of Psi_N
  real*8 function get_dperp3(psin,D_perp_sp)
#if _OPENMP >= 201511
    !$omp declare simd
#endif
    implicit none
    
    real*8, intent(in)                     :: psin
    real*8, intent(in)                     :: D_perp_sp(10)

    get_dperp3 = D_perp_sp(1) * ( (1.d0-D_perp_sp(2)) +  &
      D_perp_sp(2)*(0.5d0 - 0.5d0*tanh((psin-D_perp_sp(5))/D_perp_sp(4))) )
      
    if ( jorek_model >= 300 ) then
        
      get_dperp3 = get_dperp3 + D_perp_sp(6)*D_perp_sp(2) *   &
        ((0.5d0 - 0.5d0*tanh((-psin+D_perp_sp(5)+D_perp_sp(3)) /D_perp_sp(4))))

    end if
        
  end function get_dperp3
  
  !> Determine perpendicular particle diffusivity, D_perp, as a function of Psi_N
  real*8 function get_dperp4(psin,num_d_prof_x,num_d_prof_y,num_d_prof_len)
#if _OPENMP >= 201511
    !$omp declare simd
#endif
    implicit none
    
    real*8, intent(in)                     :: psin
    real*8, intent(in), allocatable        :: num_d_prof_x(:) !<Given numerical profile
    real*8, intent(in), allocatable        :: num_d_prof_y(:) !<Given numerical profile
    integer, intent(in)                    :: num_d_prof_len  !<Length of given numerical profile

    get_dperp4 = interpolProf(num_d_prof_x, num_d_prof_y, num_d_prof_len, psin)
    
  end function get_dperp4
  
  !> Determine perpendicular particle diffusivity, D_perp, as a function of Psi_N and Z
  real*8 function get_dperp5(psin, Z)
#if _OPENMP >= 201511
    !$omp declare simd
#endif
    implicit none
    
    real*8, intent(in) :: psin, Z
    
    if ( num_d_perp ) then
      
      get_dperp5 = interpolProf(num_d_perp_x, num_d_perp_y, num_d_perp_len, psin)
      
    else
      
      get_dperp5 = D_perp(1) * ( (1.d0-D_perp(2)) +                                                 &
        D_perp(2)*(0.5d0 - 0.5d0*tanh((psin-D_perp(5))/D_perp(4))) )
      
      if ( jorek_model >= 300 ) then
        
        get_dperp5 = get_dperp5 + D_perp(6)*D_perp(2) *                                              &
          ((0.5d0 - 0.5d0*tanh((-psin+D_perp(5)+D_perp(3)) /D_perp(4))))
        
      end if
      
    end if

    if ( D_perp_pri_z0(1) .ne. 99.d0 ) then
      get_dperp5 = get_dperp5 + D_perp_pri(1) * 0.5d0*(1.d0 -                                             &
                           tanh( (Z-D_perp_pri_z0(1)) / D_perp_pri_width(1) ))
    end if
    
    if ( D_perp_pri_z0(2) .ne. 99.d0 ) then
      get_dperp5 = get_dperp5 + D_perp_pri(2) * 0.5d0*(1.d0 +                                             &
                           tanh( (Z-D_perp_pri_z0(2)) / D_perp_pri_width(2) ))
    end if
    
  end function get_dperp5  
  
  !> Determine perpendicular heat diffusivity, ZK_perp, as a function of Psi_N
  real*8 function get_zkperp1(psin)
#if _OPENMP >= 201511
    !$omp declare simd
#endif
    implicit none
    
    real*8, intent(in)  :: psin 
    
    
    if ( num_zk_perp ) then
      
      get_zkperp1 = interpolProf(num_zk_perp_x, num_zk_perp_y, num_zk_perp_len, psin)
      
    else
      
      get_zkperp1 = ZK_perp(1) * ( (1.d0-ZK_perp(2)) +                                              &
        ZK_perp(2) *(0.5d0 - 0.5d0*tanh((psin-ZK_perp(5))/ZK_perp(4))) )
      
      if ( jorek_model >= 300 ) then
        
        get_zkperp1 = get_zkperp1 + ZK_perp(6)*ZK_perp(2) *                                          &
          ((0.5d0 - 0.5d0*tanh((-psin+ZK_perp(5)+ZK_perp(3)) /ZK_perp(4))))
        
      end if
      
    end if
    
  end function get_zkperp1


  !> Determine perpendicular heat diffusivity, ZK_perp, as a function of Psi_N and Z
  real*8 function get_zkperp3(psin, Z)
#if _OPENMP >= 201511
    !$omp declare simd
#endif
    implicit none
    
    real*8, intent(in)  :: psin, Z
    
    
    if ( num_zk_perp ) then
      
      get_zkperp3 = interpolProf(num_zk_perp_x, num_zk_perp_y, num_zk_perp_len, psin)
      
    else
      
      get_zkperp3 = ZK_perp(1) * ( (1.d0-ZK_perp(2)) +                                              &
        ZK_perp(2) *(0.5d0 - 0.5d0*tanh((psin-ZK_perp(5))/ZK_perp(4))) )
      
      if ( jorek_model >= 300 ) then
        
        get_zkperp3 = get_zkperp3 + ZK_perp(6)*ZK_perp(2) *                                          &
          ((0.5d0 - 0.5d0*tanh((-psin+ZK_perp(5)+ZK_perp(3)) /ZK_perp(4))))
        
      end if
      
    end if

    if ( ZK_perp_pri_z0(1) .ne. 99.d0 ) then
      get_zkperp3 = get_zkperp3 + ZK_perp_pri(1) * 0.5d0*(1.d0-                                              &
                           tanh( (Z-ZK_perp_pri_z0(1)) / ZK_perp_pri_width(1) ))
    end if
 
    if ( ZK_perp_pri_z0(2) .ne. 99.d0 ) then
      get_zkperp3 = get_zkperp3 + ZK_perp_pri(2) * 0.5d0*(1.d0 +                                             &
                           tanh( (Z-ZK_perp_pri_z0(2)) / ZK_perp_pri_width(2) ))
    end if
  
    
  end function get_zkperp3


  !> Determine perpendicular ion heat diffusivity, ZK_perp, as a function of Psi_N and Z
  real*8 function get_zk_iperp3(psin, Z)
#if _OPENMP >= 201511
    !$omp declare simd
#endif
    implicit none
    
    real*8, intent(in)  :: psin, Z
    
    
    if ( num_zk_i_perp ) then
      
      get_zk_iperp3 = interpolProf(num_zk_i_perp_x, num_zk_i_perp_y, num_zk_i_perp_len, psin)
      
    else
      
      get_zk_iperp3 = ZK_i_perp(1) * ( (1.d0-ZK_i_perp(2)) +                                              &
        ZK_i_perp(2) *(0.5d0 - 0.5d0*tanh((psin-ZK_i_perp(5))/ZK_i_perp(4))) )
      
      if ( jorek_model >= 300 ) then
        
        get_zk_iperp3 = get_zk_iperp3 + ZK_i_perp(6)*ZK_i_perp(2) *                                          &
          ((0.5d0 - 0.5d0*tanh((-psin+ZK_i_perp(5)+ZK_i_perp(3)) /ZK_i_perp(4))))
        
      end if
      
    end if

    if ( ZK_i_perp_pri_z0(1) .ne. 99.d0 ) then
      get_zk_iperp3 = get_zk_iperp3 + ZK_i_perp_pri(1) * 0.5d0*(1.d0 -                                             &
                           tanh( (Z-ZK_i_perp_pri_z0(1)) / ZK_i_perp_pri_width(1) ))
    end if
 
    if ( ZK_i_perp_pri_z0(2) .ne. 99.d0 ) then
      get_zk_iperp3 = get_zk_iperp3 + ZK_i_perp_pri(2) * 0.5d0*(1.d0 +                                             &
                           tanh( (Z-ZK_i_perp_pri_z0(2)) / ZK_i_perp_pri_width(2) ))
    end if
  
    
  end function get_zk_iperp3



  !> Determine perpendicular electron heat diffusivity, ZK_e_perp, as a function of Psi_N and Z
  real*8 function get_zk_eperp3(psin, Z)
#if _OPENMP >= 201511
    !$omp declare simd
#endif
    implicit none
    
    real*8, intent(in)  :: psin, Z
    
    
    if ( num_zk_e_perp ) then
      
      get_zk_eperp3 = interpolProf(num_zk_e_perp_x, num_zk_e_perp_y, num_zk_e_perp_len, psin)
      
    else
      
      get_zk_eperp3 = ZK_e_perp(1) * ( (1.d0-ZK_e_perp(2)) +                                              &
        ZK_e_perp(2) *(0.5d0 - 0.5d0*tanh((psin-ZK_e_perp(5))/ZK_e_perp(4))) )
      
      if ( jorek_model >= 300 ) then
        
        get_zk_eperp3 = get_zk_eperp3 + ZK_e_perp(6)*ZK_e_perp(2) *                                          &
          ((0.5d0 - 0.5d0*tanh((-psin+ZK_e_perp(5)+ZK_e_perp(3)) /ZK_e_perp(4))))
        
      end if
      
    end if

    if ( ZK_e_perp_pri_z0(1) .ne. 99.d0 ) then
      get_zk_eperp3 = get_zk_eperp3 + ZK_e_perp_pri(1) * 0.5d0*(1.d0 -                                             &
                           tanh( (Z-ZK_e_perp_pri_z0(1)) / ZK_e_perp_pri_width(1) ))
    end if
 
    if ( ZK_e_perp_pri_z0(2) .ne. 99.d0 ) then
      get_zk_eperp3 = get_zk_eperp3 + ZK_e_perp_pri(2) * 0.5d0*(1.d0 +                                             &
                           tanh( (Z-ZK_e_perp_pri_z0(2)) / ZK_e_perp_pri_width(2) ))
    end if
  
    
  end function get_zk_eperp3


   !> Determine perpendicular heat diffusivity, ZK_perp, as a function of Psi_N, for ions
  real*8 function get_zk_iperp1(psin)
    
    implicit none
    
    real*8, intent(in)  :: psin 
    
    
    if ( num_zk_i_perp ) then
      
      get_zk_iperp1 = interpolProf(num_zk_i_perp_x, num_zk_i_perp_y, num_zk_i_perp_len, psin)
      
    else
      
      get_zk_iperp1 = ZK_i_perp(1) * ( (1.d0-ZK_i_perp(2)) +                                              &
        ZK_i_perp(2) *(0.5d0 - 0.5d0*tanh((psin-ZK_i_perp(5))/ZK_i_perp(4))) )
      
      if ( jorek_model >= 300 ) then
        
        get_zk_iperp1 = get_zk_iperp1 + ZK_i_perp(6)*ZK_i_perp(2) *                                          &
          ((0.5d0 - 0.5d0*tanh((-psin+ZK_i_perp(5)+ZK_i_perp(3)) /ZK_i_perp(4))))
        
      end if
      
    end if
    
  end function get_zk_iperp1
	
   !> Determine perpendicular heat diffusivity, ZK_perp, as a function of Psi_N, for electrons
  real*8 function get_zk_eperp1(psin)
    
    implicit none
    
    real*8, intent(in)  :: psin 
    
    
    if ( num_zk_e_perp ) then
      
      get_zk_eperp1 = interpolProf(num_zk_e_perp_x, num_zk_e_perp_y, num_zk_e_perp_len, psin)
      
    else
      
      get_zk_eperp1 = ZK_e_perp(1) * ( (1.d0-ZK_e_perp(2)) +                                              &
        ZK_e_perp(2) *(0.5d0 - 0.5d0*tanh((psin-ZK_e_perp(5))/ZK_e_perp(4))) )
      
      if ( jorek_model >= 300 ) then
        
        get_zk_eperp1 = get_zk_eperp1 + ZK_e_perp(6)*ZK_e_perp(2) *                                          &
          ((0.5d0 - 0.5d0*tanh((-psin+ZK_e_perp(5)+ZK_e_perp(3)) /ZK_e_perp(4))))
        
      end if
      
    end if
    
  end function get_zk_eperp1
	
	
	
  
  !> Determine perpendicular particle diffusivity, D_perp, as a function of Psi_N
  real*8 function get_dperp2(psi, psi_norm, psi_axis, psi_bnd, Z, Z_xpoint)
    
    implicit none
    
    real*8, intent(in) :: psi, psi_norm, psi_axis, psi_bnd, Z, Z_xpoint(2)
    real*8         :: psi_D
    real*8         :: atn_D, datn_D, atn_D_n, pol_D, dpol_D, D_min
    real*8         :: Diff(1:5)
    
    ! --- Numerical profile
    if ( num_d_perp ) then
      
      get_dperp2 = interpolProf(num_d_perp_x, num_d_perp_y, num_d_perp_len, psi_norm)
      
    ! --- Normal (analytical) profiles
    else
      
      ! --- Old profile, note, Diff(10) is our switch: =0 => use old profile, =1 => use new profile
      if (D_perp(10) .ne. 1.d0) then
        
        atn_D = (0.5d0 - 0.5d0*tanh((psi_norm-D_perp(5))/D_perp(4)))
        get_dperp2 = D_perp(1) * ( (1.d0-D_perp(2)) + D_perp(2) * atn_D )
        if (jorek_model >= 300) then
          atn_D = ((0.5d0 - 0.5d0*tanh((-psi_norm+D_perp(5)+D_perp(3)) /D_perp(4))))
          get_dperp2 = get_dperp2 + D_perp(6) * D_perp(2) * atn_D
        endif
  
      ! --- Adapted profile, going like Dperp ~ 1/grad(rho) to obtain constand perp flux in pedestal
      else  
        ! --- Take values from input file (rho_coef, T_coef...) 
        Diff(1) = rho_coef(1)
        Diff(2) = rho_coef(2)
        Diff(3) = rho_coef(3)
        Diff(4) = rho_coef(4)
        Diff(5) = rho_coef(5)
        
        ! --- Correct for hollow profiles (otherwise D_perp > infinity)
        if (Diff(1) .ge. 0.d0) Diff(1) = -0.1d0

        ! --- First compute the normalisation factor (this is the profile at a given psi)
        ! --- Note that if you choose psi_D=Diff(5) you normalise with the minimal value of D_perp
        psi_D   = 0.8       ! At the pedestal top, roughly
        !psi_D   = Diff(5)   ! At the point of steepest gradient (ie. will be the smallest D_perp)
        
        ! --- The tanh part of the initial profile and its derivative 
        atn_D   = 0.5d0 - 0.5d0 * tanh((psi_D-Diff(5))/Diff(4))
        datn_D   =   - 0.5d0 / cosh((psi_D-Diff(5))/Diff(4))**2.d0 /(Diff(4)*(psi_bnd - psi_axis))
        
        ! --- The polynomial part of the initial profile and its derivative
        pol_D   = 1 + Diff(1)*psi_D    + Diff(2)*psi_D**2.d0   + Diff(3)*psi_D**3.d0
        dpol_D   =    (Diff(1)     + 2.d0*Diff(2)*psi_D     + 3.d0*Diff(3)*psi_D**2.d0)/(psi_bnd - psi_axis)
        
        ! --- The normalisation factor
        D_min   = 1.d0 / ( dpol_D*atn_D + pol_D*datn_D )

        ! --- Then we recompute the profile derivatives for the actual local psi_norm
        ! --- Be careful at core and in SOL, where gradients are ~0 => infinite 1/grad(p)
        ! --- Take a symetric profile at from middle of pedestal outwards, because SOL is too dangerous (gradient too close to zero...)
        if (psi_norm .gt. Diff(5)) then
          psi_D = 2.d0*Diff(5) - psi_norm
        else
          psi_D = psi_norm
        endif
        
        if (psi_norm .lt. 0.5d0) psi_D = 0.5d0
        if (xcase .ne. UPPER_XPOINT)  psi_D = psi_D * (0.5d0 - 0.5d0 * tanh((Z_xpoint(1)-Z)/0.1d0))
        if (xcase .ne. LOWER_XPOINT)  psi_D = psi_D * (0.5d0 - 0.5d0 * tanh((Z-Z_xpoint(2))/0.1d0))
        
        ! --- The tanh part of the initial profile and its derivative 
        atn_D   = 0.5d0 - 0.5d0 * tanh((psi_D-Diff(5))/Diff(4))
        datn_D   =   - 0.5d0 / cosh((psi_D-Diff(5))/Diff(4))**2.d0 /(Diff(4)*(psi_bnd - psi_axis))
        
        ! --- The polynomial part of the initial profile and its derivative
        pol_D   = 1 + Diff(1)*psi_D    + Diff(2)*psi_D**2.d0   + Diff(3)*psi_D**3.d0
        dpol_D   =    (Diff(1)     + 2.d0*Diff(2)*psi_D     + 3.d0*Diff(3)*psi_D**2.d0)/(psi_bnd - psi_axis)

        ! --- Build the profile (goes like 1/grad(density_profile) in the pedestal region
        get_dperp2 = D_perp(1) / (dpol_D*atn_D + pol_D*datn_D) / D_min 
        
      end if
      
    end if
    
  end function get_dperp2
  
  
  
  !> Determine perpendicular heat diffusivity, ZK_perp, as a function of Psi_N
  real*8 function get_zkperp2(psi, psi_norm, psi_axis, psi_bnd, Z, Z_xpoint)
    
    implicit none
    
    real*8, intent(in) :: psi, psi_norm, psi_axis, psi_bnd, Z, Z_xpoint(2)
    real*8         :: psi_D
    real*8         :: atn_D, datn_D, atn_D_n, pol_D, dpol_D, D_min
    real*8         :: Diff(1:5)
    real*8             :: rho_norm
    real*8             :: zn, dn_dpsi, dn_dz, dn_dpsi2, dn_dz2, dn_dpsi_dz, dn_dpsi3, dn_dpsi_dz2, dn_dpsi2_dz 
    
    
    ! --- Numerical profile
    if ( num_zk_perp ) then
      
      get_zkperp2 = interpolProf(num_zk_perp_x, num_zk_perp_y, num_zk_perp_len, psi_norm)
      
    ! --- Normal (analytical) profiles
    else
      
      ! --- Old profile, note, Diff(10) is our switch: =0 => use old profile, =1 => use new profile
      if (ZK_perp(10) .ne. 1.d0) then
        
        atn_D = (0.5d0 - 0.5d0*tanh((psi_norm-ZK_perp(5))/ZK_perp(4)))
        get_zkperp2 = ZK_perp(1) * ( (1.d0-ZK_perp(2)) + ZK_perp(2) * atn_D )
        
        if ( jorek_model >= 300 ) then
          atn_D = ((0.5d0 - 0.5d0*tanh((-psi_norm+ZK_perp(5)+ZK_perp(3)) /ZK_perp(4))))
          get_zkperp2 = get_zkperp2 + ZK_perp(6) * ZK_perp(2) * atn_D
        endif

      ! --- Adapted profile, going like Dperp ~ 1/grad(rho) to obtain constand perp flux in pedestal
      else
        
        ! --- Take values from input file (rho_coef, T_coef...) 
        Diff(1) = T_coef(1)
        Diff(2) = T_coef(2)
        Diff(3) = T_coef(3)
        Diff(4) = T_coef(4)
        Diff(5) = T_coef(5)
        
        ! --- Correct for hollow profiles (otherwise D_perp > infinity)
        if (Diff(1) .ge. 0.d0) Diff(1) = -0.1d0

        ! --- First compute the normalisation factor (this is the profile at a given psi)
        ! --- Note that if you choose psi_D=Diff(5) you normalise with the minimal value of D_perp
        psi_D   = 0.8       ! At the pedestal top, roughly
        !psi_D   = Diff(5)   ! At the point of steepest gradient (ie. will be the smallest D_perp)
        ! --- The tanh part of the initial profile and its derivative 
        atn_D   = 0.5d0 - 0.5d0 * tanh((psi_D-Diff(5))/Diff(4))
        datn_D   =   - 0.5d0 / cosh((psi_D-Diff(5))/Diff(4))**2.d0 /(Diff(4)*(psi_bnd - psi_axis))
        ! --- The polynomial part of the initial profile and its derivative
        pol_D   = 1 + Diff(1)*psi_D    + Diff(2)*psi_D**2.d0   + Diff(3)*psi_D**3.d0
        dpol_D   =    (Diff(1)     + 2.d0*Diff(2)*psi_D     + 3.d0*Diff(3)*psi_D**2.d0)/(psi_bnd - psi_axis)
        ! --- The normalisation factor
        D_min   = 1.d0 / ( dpol_D*atn_D + pol_D*datn_D )

        ! --- Then we recompute the profile derivatives for the actual local psi_norm
        ! --- Be careful at core and in SOL, where gradients are ~0 => infinite 1/grad(p)
        ! --- Take a symetric profile at from middle of pedestal outwards, because SOL is too dangerous (gradient too close to zero...)
        if (psi_norm .gt. Diff(5)) then
          psi_D = 2.d0*Diff(5) - psi_norm
        else
          psi_D = psi_norm
        endif
        
        if (psi_norm .lt. 0.5d0) psi_D = 0.5d0
        if (xcase .ne. UPPER_XPOINT)  psi_D = psi_D * (0.5d0 - 0.5d0 * tanh((Z_xpoint(1)-Z)/0.1d0))
        if (xcase .ne. LOWER_XPOINT)  psi_D = psi_D * (0.5d0 - 0.5d0 * tanh((Z-Z_xpoint(2))/0.1d0))
        
        ! --- The tanh part of the initial profile and its derivative 
        atn_D   = 0.5d0 - 0.5d0 * tanh((psi_D-Diff(5))/Diff(4))
        datn_D   =   - 0.5d0 / cosh((psi_D-Diff(5))/Diff(4))**2.d0 /(Diff(4)*(psi_bnd - psi_axis))
        
        ! --- The polynomial part of the initial profile and its derivative
        pol_D   = 1 + Diff(1)*psi_D    + Diff(2)*psi_D**2.d0   + Diff(3)*psi_D**3.d0
        dpol_D   =    (Diff(1)     + 2.d0*Diff(2)*psi_D     + 3.d0*Diff(3)*psi_D**2.d0)/(psi_bnd - psi_axis)

        ! --- Build the profile (goes like 1/grad(density_profile) in the pedestal region
        get_zkperp2 = ZK_perp(1) / (dpol_D*atn_D + pol_D*datn_D) / D_min 

        ! --- K-perp is the pressure diffusion, not the temperature diffusion, so need to divide by rho
        psi_D = psi
        call density(xpoint, xcase, Z, Z_xpoint, psi_D,psi_axis,psi_bnd, &
               zn,dn_dpsi, dn_dz, dn_dpsi2, dn_dz2, dn_dpsi_dz, dn_dpsi3, dn_dpsi_dz2, dn_dpsi2_dz)
        rho_norm        = zn / rho_0 ! normalise to 1.0 in core
        get_zkperp2 = get_zkperp2 / rho_norm
  
      end if
      
    end if
    
  end function get_zkperp2


  !> Determine perpendicular heat diffusivity, ZK_perp, as a function of Psi_N, for ions
  real*8 function get_zk_iperp2(psi, psi_norm, psi_axis, psi_bnd, Z, Z_xpoint)
    
    implicit none
    
    real*8, intent(in) :: psi, psi_norm, psi_axis, psi_bnd, Z, Z_xpoint(2)
    real*8         :: psi_D
    real*8         :: atn_D, datn_D, atn_D_n, pol_D, dpol_D, D_min
    real*8         :: Diff(1:5)
    real*8             :: rho_norm
    real*8             :: zn, dn_dpsi, dn_dz, dn_dpsi2, dn_dz2, dn_dpsi_dz, dn_dpsi3, dn_dpsi_dz2, dn_dpsi2_dz 
    
    
    ! --- Numerical profile
    if ( num_zk_i_perp ) then
      
      get_zk_iperp2 = interpolProf(num_zk_i_perp_x, num_zk_i_perp_y, num_zk_i_perp_len, psi_norm)
      
    ! --- Normal (analytical) profiles
    else
      
      ! --- Old profile, note, Diff(10) is our switch: =0 => use old profile, =1 => use new profile
      if (ZK_i_perp(10) .ne. 1.d0) then
        
        atn_D = (0.5d0 - 0.5d0*tanh((psi_norm-ZK_i_perp(5))/ZK_i_perp(4)))
        get_zk_iperp2 = ZK_i_perp(1) * ( (1.d0-ZK_i_perp(2)) + ZK_i_perp(2) * atn_D )
        
        if ( jorek_model >= 300 ) then
          atn_D = ((0.5d0 - 0.5d0*tanh((-psi_norm+ZK_i_perp(5)+ZK_i_perp(3)) /ZK_i_perp(4))))
          get_zk_iperp2 = get_zk_iperp2 + ZK_i_perp(6) * ZK_i_perp(2) * atn_D
        endif

      ! --- Adapted profile, going like Dperp ~ 1/grad(rho) to obtain constand perp flux in pedestal
      else
        
        ! --- Take values from input file (rho_coef, T_coef...) 
        Diff(1) = Ti_coef(1)
        Diff(2) = Ti_coef(2)
        Diff(3) = Ti_coef(3)
        Diff(4) = Ti_coef(4)
        Diff(5) = Ti_coef(5)
        
        ! --- Correct for hollow profiles (otherwise D_perp > infinity)
        if (Diff(1) .ge. 0.d0) Diff(1) = -0.1d0

        ! --- First compute the normalisation factor (this is the profile at a given psi)
        ! --- Note that if you choose psi_D=Diff(5) you normalise with the minimal value of D_perp
        psi_D   = 0.8       ! At the pedestal top, roughly
        !psi_D   = Diff(5)   ! At the point of steepest gradient (ie. will be the smallest D_perp)
        ! --- The tanh part of the initial profile and its derivative 
        atn_D   = 0.5d0 - 0.5d0 * tanh((psi_D-Diff(5))/Diff(4))
        datn_D   =   - 0.5d0 / cosh((psi_D-Diff(5))/Diff(4))**2.d0 /(Diff(4)*(psi_bnd - psi_axis))
        ! --- The polynomial part of the initial profile and its derivative
        pol_D   = 1 + Diff(1)*psi_D    + Diff(2)*psi_D**2.d0   + Diff(3)*psi_D**3.d0
        dpol_D   =    (Diff(1)     + 2.d0*Diff(2)*psi_D     + 3.d0*Diff(3)*psi_D**2.d0)/(psi_bnd - psi_axis)
        ! --- The normalisation factor
        D_min   = 1.d0 / ( dpol_D*atn_D + pol_D*datn_D )

        ! --- Then we recompute the profile derivatives for the actual local psi_norm
        ! --- Be careful at core and in SOL, where gradients are ~0 => infinite 1/grad(p)
        ! --- Take a symetric profile at from middle of pedestal outwards, because SOL is too dangerous (gradient too close to zero...)
        if (psi_norm .gt. Diff(5)) then
          psi_D = 2.d0*Diff(5) - psi_norm
        else
          psi_D = psi_norm
        endif
        
        if (psi_norm .lt. 0.5d0) psi_D = 0.5d0
        if (xcase .ne. UPPER_XPOINT)  psi_D = psi_D * (0.5d0 - 0.5d0 * tanh((Z_xpoint(1)-Z)/0.1d0))
        if (xcase .ne. LOWER_XPOINT)  psi_D = psi_D * (0.5d0 - 0.5d0 * tanh((Z-Z_xpoint(2))/0.1d0))
        
        ! --- The tanh part of the initial profile and its derivative 
        atn_D   = 0.5d0 - 0.5d0 * tanh((psi_D-Diff(5))/Diff(4))
        datn_D   =   - 0.5d0 / cosh((psi_D-Diff(5))/Diff(4))**2.d0 /(Diff(4)*(psi_bnd - psi_axis))
        
        ! --- The polynomial part of the initial profile and its derivative
        pol_D   = 1 + Diff(1)*psi_D    + Diff(2)*psi_D**2.d0   + Diff(3)*psi_D**3.d0
        dpol_D   =    (Diff(1)     + 2.d0*Diff(2)*psi_D     + 3.d0*Diff(3)*psi_D**2.d0)/(psi_bnd - psi_axis)

        ! --- Build the profile (goes like 1/grad(density_profile) in the pedestal region
        get_zk_iperp2 = ZK_i_perp(1) / (dpol_D*atn_D + pol_D*datn_D) / D_min 

        ! --- K-perp is the pressure diffusion, not the temperature diffusion, so need to divide by rho
        psi_D = psi
        call density(xpoint, xcase, Z, Z_xpoint, psi_D,psi_axis,psi_bnd, &
               zn,dn_dpsi, dn_dz, dn_dpsi2, dn_dz2, dn_dpsi_dz, dn_dpsi3, dn_dpsi_dz2, dn_dpsi2_dz)
        rho_norm        = zn / rho_0 ! normalise to 1.0 in core
        get_zk_iperp2 = get_zk_iperp2 / rho_norm
  
      end if
      
    end if
    
  end function get_zk_iperp2


 !> Determine perpendicular heat diffusivity, ZK_perp, as a function of Psi_N, for electrons
  real*8 function get_zk_eperp2(psi, psi_norm, psi_axis, psi_bnd, Z, Z_xpoint)
    
    implicit none
    
    real*8, intent(in) :: psi, psi_norm, psi_axis, psi_bnd, Z, Z_xpoint(2)
    real*8         :: psi_D
    real*8         :: atn_D, datn_D, atn_D_n, pol_D, dpol_D, D_min
    real*8         :: Diff(1:5)
    real*8             :: rho_norm
    real*8             :: zn, dn_dpsi, dn_dz, dn_dpsi2, dn_dz2, dn_dpsi_dz, dn_dpsi3, dn_dpsi_dz2, dn_dpsi2_dz 
    
    
    ! --- Numerical profile
    if ( num_zk_e_perp ) then
      
      get_zk_eperp2 = interpolProf(num_zk_e_perp_x, num_zk_e_perp_y, num_zk_e_perp_len, psi_norm)
      
    ! --- Normal (analytical) profiles
    else
      
      ! --- Old profile, note, Diff(10) is our switch: =0 => use old profile, =1 => use new profile
      if (ZK_e_perp(10) .ne. 1.d0) then
        
        atn_D = (0.5d0 - 0.5d0*tanh((psi_norm-ZK_e_perp(5))/ZK_e_perp(4)))
        get_zk_eperp2 = ZK_e_perp(1) * ( (1.d0-ZK_e_perp(2)) + ZK_e_perp(2) * atn_D )
        
        if ( jorek_model >= 300 ) then
          atn_D = ((0.5d0 - 0.5d0*tanh((-psi_norm+ZK_e_perp(5)+ZK_e_perp(3)) /ZK_e_perp(4))))
          get_zk_eperp2 = get_zk_eperp2 + ZK_e_perp(6) * ZK_e_perp(2) * atn_D
        endif

      ! --- Adapted profile, going like Dperp ~ 1/grad(rho) to obtain constand perp flux in pedestal
      else
        
        ! --- Take values from input file (rho_coef, T_coef...) 
        Diff(1) = Te_coef(1)
        Diff(2) = Te_coef(2)
        Diff(3) = Te_coef(3)
        Diff(4) = Te_coef(4)
        Diff(5) = Te_coef(5)
        
        ! --- Correct for hollow profiles (otherwise D_perp > infinity)
        if (Diff(1) .ge. 0.d0) Diff(1) = -0.1d0

        ! --- First compute the normalisation factor (this is the profile at a given psi)
        ! --- Note that if you choose psi_D=Diff(5) you normalise with the minimal value of D_perp
        psi_D   = 0.8       ! At the pedestal top, roughly
        !psi_D   = Diff(5)   ! At the point of steepest gradient (ie. will be the smallest D_perp)
        ! --- The tanh part of the initial profile and its derivative 
        atn_D   = 0.5d0 - 0.5d0 * tanh((psi_D-Diff(5))/Diff(4))
        datn_D   =   - 0.5d0 / cosh((psi_D-Diff(5))/Diff(4))**2.d0 /(Diff(4)*(psi_bnd - psi_axis))
        ! --- The polynomial part of the initial profile and its derivative
        pol_D   = 1 + Diff(1)*psi_D    + Diff(2)*psi_D**2.d0   + Diff(3)*psi_D**3.d0
        dpol_D   =    (Diff(1)     + 2.d0*Diff(2)*psi_D     + 3.d0*Diff(3)*psi_D**2.d0)/(psi_bnd - psi_axis)
        ! --- The normalisation factor
        D_min   = 1.d0 / ( dpol_D*atn_D + pol_D*datn_D )

        ! --- Then we recompute the profile derivatives for the actual local psi_norm
        ! --- Be careful at core and in SOL, where gradients are ~0 => infinite 1/grad(p)
        ! --- Take a symetric profile at from middle of pedestal outwards, because SOL is too dangerous (gradient too close to zero...)
        if (psi_norm .gt. Diff(5)) then
          psi_D = 2.d0*Diff(5) - psi_norm
        else
          psi_D = psi_norm
        endif
        
        if (psi_norm .lt. 0.5d0) psi_D = 0.5d0
        if (xcase .ne. UPPER_XPOINT)  psi_D = psi_D * (0.5d0 - 0.5d0 * tanh((Z_xpoint(1)-Z)/0.1d0))
        if (xcase .ne. LOWER_XPOINT)  psi_D = psi_D * (0.5d0 - 0.5d0 * tanh((Z-Z_xpoint(2))/0.1d0))
        
        ! --- The tanh part of the initial profile and its derivative 
        atn_D   = 0.5d0 - 0.5d0 * tanh((psi_D-Diff(5))/Diff(4))
        datn_D   =   - 0.5d0 / cosh((psi_D-Diff(5))/Diff(4))**2.d0 /(Diff(4)*(psi_bnd - psi_axis))
        
        ! --- The polynomial part of the initial profile and its derivative
        pol_D   = 1 + Diff(1)*psi_D    + Diff(2)*psi_D**2.d0   + Diff(3)*psi_D**3.d0
        dpol_D   =    (Diff(1)     + 2.d0*Diff(2)*psi_D     + 3.d0*Diff(3)*psi_D**2.d0)/(psi_bnd - psi_axis)

        ! --- Build the profile (goes like 1/grad(density_profile) in the pedestal region
        get_zk_eperp2 = ZK_e_perp(1) / (dpol_D*atn_D + pol_D*datn_D) / D_min 

        ! --- K-perp is the pressure diffusion, not the temperature diffusion, so need to divide by rho
        psi_D = psi
        call density(xpoint, xcase, Z, Z_xpoint, psi_D,psi_axis,psi_bnd, &
               zn,dn_dpsi, dn_dz, dn_dpsi2, dn_dz2, dn_dpsi_dz, dn_dpsi3, dn_dpsi_dz2, dn_dpsi2_dz)
        rho_norm        = zn / rho_0 ! normalise to 1.0 in core
        get_zk_eperp2 = get_zk_eperp2 / rho_norm
  
      end if
      
    end if
    
  end function get_zk_eperp2

  !> Returns the inward pinch velocity profile for background fluid.
  !> If v_pinch_file /= 'none', the profile is interpolated from that file (two-column ASCII: psin, V_pinch).
  !> Otherwise falls back to the Gaussian: V_pinch_gauss * exp(-(psin - V_pinch_psin)^2 / V_pinch_sig^2).
  !> Positive values drive density inward (toward magnetic axis).
  real*8 function get_vpinch(psin)
    implicit none

    real*8, intent(in) :: psin

    if ( num_v_pinch ) then
      get_vpinch = interpolProf(num_v_pinch_x, num_v_pinch_y, num_v_pinch_len, psin)
    else
      get_vpinch = V_pinch_gauss * exp(-(psin - V_pinch_psin)**2 / V_pinch_sig**2)
    end if

  end function get_vpinch

end module diffusivities
