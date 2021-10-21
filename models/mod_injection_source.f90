module mod_injection_source

  use constants

  real*8 :: total_n_particles_inj     = 0.
  real*8 :: total_n_particles         = 0.
  real*8 :: total_n_particles_inj_all = 0.

  contains 


  integer function factorial(n)

    implicit none

    integer, intent(in) :: n 
    integer             :: i, Ans

    Ans = 1
 
    do i=1,n
      Ans = Ans * i
    enddo

    factorial = Ans

  end function factorial



  subroutine inj_source(ns_amplitude,ns_R,ns_Z,ns_phi,ns_radius,ns_deltaphi,ns_tor_norm,  &
                        A_Dmv,K_Dmv,V_Dmv,P_Dmv,t_ns,L_tube,R,Z,phi,rhon_source,t_now,                  &
                        JET_MGI,ASDEX_MGI,central_density,central_mass)

  !=================================================================================
  !  This subroutine computes the atom/ion number density source for a realistic Deuterium
  !  MGI in JET (if ns_timedependent is .t.).
  !  If ns_timedependent is .f., this routine computes a constant source in time
  !  where the main parameter is ns_amplitude
  !  More details in the JOREK wiki or by asking A.Fil or E.Nardon
  !=================================================================================

    use phys_module, only: imp_type

    implicit none

    real*8 :: c0_gas                   ! Sound velocity of gas in reservoir
    integer:: n_gas                    ! = 2/(gamma-1) where gamma = heat capacity ratio of gas
    real*8 :: A_gas                    ! Atomic number of gas particles
    real*8 :: mass_gas                 ! Mass of a gas particles
    real*8 :: mol_atom                 ! Number of atoms in a molecular
    real*8 :: radius
    real*8 :: ns_tor_shape
    real*8 :: ns_pol_shape
    real*8 :: dphi
    real*8 :: V_ns
    real*8 :: f_Nbar
    real*8 :: f_dNbar_dt
    real*8 :: ns_dNinj_dt
    real*8 :: ns_drhon_dt
    real*8 :: t_loc
    real*8 :: t_norm
    real*8 :: prof_temp
    real*8 :: R_Asdex
    real*8 :: mnum
    real*8 :: kst
    real*8 :: yy
    real*8 :: gam
    real*8 :: dt_open
    real*8 :: N_barlitre
    integer:: k
    real*8, intent(in)  :: R
    real*8, intent(in)  :: Z
    real*8, intent(in)  :: phi
    real*8, intent(in)  :: A_Dmv
    real*8, intent(in)  :: K_Dmv
    real*8, intent(in)  :: V_Dmv
    real*8, intent(in)  :: P_Dmv
    real*8, intent(in)  :: t_now
    real*8, intent(in)  :: t_ns
    real*8, intent(in)  :: ns_amplitude
    real*8, intent(in)  :: ns_R
    real*8, intent(in)  :: ns_Z
    real*8, intent(in)  :: ns_phi
    real*8, intent(in)  :: ns_radius
    real*8, intent(in)  :: ns_deltaphi
    real*8, intent(in)  :: L_tube
    real*8, intent(in)  :: central_density
    real*8, intent(in)  :: central_mass
    real*8              :: DMV_inj_frac
    logical, intent(in) :: JET_MGI
    logical, intent(in) :: ASDEX_MGI
    real*8, intent(out) :: rhon_source  ! This is in number desntiy
    real*8, intent(in)  :: ns_tor_norm

    select case ( trim(imp_type(1)) )
      case('D2')
        n_gas  = 5
        A_gas  = 4.
        mol_atom = 2.
        mass_gas = A_gas*MASS_PROTON
        c0_gas = sqrt(8.3145d0*293.d0/(A_gas*1.d-3)*(7.d0/5.d0))
      case('Ar')
        n_gas  = 3
        A_gas  = 40.
        mol_atom = 1.
        mass_gas = A_gas*MASS_PROTON
        c0_gas = sqrt(8.3145d0*293.d0/(A_gas*1.d-3)*(5.d0/3.d0))
      case('Ne')
        n_gas  = 3
        A_gas  = 20.
        mol_atom = 1.
        mass_gas = A_gas*MASS_PROTON
        c0_gas = sqrt(8.3145d0*293.d0/(A_gas*1.d-3)*(5.d0/3.d0))
      case default
        write(*,*) '!! Gas type "', trim(imp_type(1)), '" unknown (in mod_injection_source.f90) !!'
        write(*,*) '=> We assume the gas is D2.'
        n_gas  = 5
        A_gas  = 4.
        mol_atom = 2.
        mass_gas = A_gas*MASS_PROTON
        c0_gas = sqrt(8.3145d0*293.d0/(A_gas*1.d-3)*(7.d0/5.d0))
    end select

    ! ===================================================================
    ! Parameters related to the spatial distribution of the gas source:

    ! A gaussian shape is chosen poloidally
    radius = sqrt((R-ns_R)**2 + (Z-ns_Z)**2)
    ns_pol_shape = exp(-(radius/ns_radius)**2.d0)  

    ! A gaussian shape is chosen toroidally
    dphi = abs(phi - ns_phi)
    if (dphi .gt. PI) dphi = 2*PI - dphi  
    ns_tor_shape = exp(-(dphi/ns_deltaphi)**2.d0)

    ! Volume used for normalization, which corresponds to the integration in space 
    ! of the product of the above shape functions
    V_ns  = PI * ns_R * ns_tor_norm * ns_radius**2.d0
    ! ===================================================================

   !==================================================================================================
   ! A shifted time is used in order to start injected gas as soon as t_now = t_ns 
   ! (note: L_tube/3c0 is the time needed for the gas to propagate in the injection tube).
    t_norm = sqrt(MU_ZERO * central_mass * MASS_PROTON * central_density * 1.d20)
    t_loc = (t_now-t_ns) * t_norm + L_tube/(3.d0 * c0_gas)
   !==================================================================================================

    if (t_loc .gt. 0.) then

      if (JET_MGI) then

       !==================================================================================================
       ! We use here the formulae derived from the eq(8) in the paper of S.A. Bozhenkov - NF 51 (2011) 
       ! which gives the normalized number of particles injected at the exit of the DMV injection tube
       ! as a function of time.
       ! The parameters used are realistic:
       ! A_Dmv: cross sectional area of the injection pipe
       ! K_Dmv: Experimental correction factor to account for the gas expansion close to the tube orifice
       ! L_tube: DMV vacuum injection tube length
       ! V_Dmv: Volume of the DMV reservoir
       ! P_Dmv: Initial pressure in the DMV reservoir, directly linked to the total number of particles
       ! in the reservoir. Expressed in bar here as it is in all MGI experiments. 
       !==================================================================================================

        f_Nbar = 0.d0
        f_dNbar_dt = 0.d0

       ! Calculation of the normalized number of particles injected per unit time, following Bozhenkov
        do k = 0,n_gas+1
          f_Nbar     = f_Nbar + (-1.d0)**(k-1)*factorial(n_gas+1)/(factorial(n_gas+1-k)*factorial(k))*(1-(n_gas*c0_gas*t_loc/L_tube)**(1-k))

          f_dNbar_dt = f_dNbar_dt &
                        + (-1.d0)**(k-1)*factorial(n_gas+1)/(factorial(n_gas+1-k)*factorial(k))*(k-1)*(n_gas*c0_gas*(L_tube)**(-1.d0))**(1-k) &
                          *t_loc**(-k)
        end do

        f_Nbar     = ((1.*n_gas)**n_gas) * ((1.*(n_gas+1.))**(-n_gas-1)) * f_Nbar
        f_dNbar_dt = ((1.*n_gas)**n_gas) * ((1.*(n_gas+1.))**(-n_gas-1)) * f_dNbar_dt

        DMV_inj_frac = A_Dmv * L_tube * K_Dmv * f_Nbar/(V_Dmv)

       ! The gas injection is stopped when the initial number of particles in the reservoir is reached 
       ! if (DMV_inj_frac .gt. 1.d0) then      
       !   f_dNbar_dt = 0.d0
       ! endif

        ! Number of injected particles per unit time, normalized to reservoir content:
        ns_dNinj_dt = A_Dmv * K_Dmv * L_tube / V_Dmv * f_dNbar_dt

        ! Mass density injected per unit time (SI units):
        ns_drhon_dt = ns_dNinj_dt * (P_Dmv * 1.d5/(K_BOLTZ * 293)) * V_Dmv * mass_gas
    
        ! Distribute gas source in space
        rhon_source = ns_drhon_dt * ns_pol_shape * ns_tor_shape / V_ns

        ! Apply JOREK normalization
        rhon_source = (MU_ZERO)**(0.5d0)*(central_mass*MASS_PROTON*central_density*1.d20)**(-0.5d0) * rhon_source

        ! Converting mass density into number density
        rhon_source = rhon_source * (central_mass * MASS_PROTON / mass_gas)
      elseif (ASDEX_MGI) then

        N_barlitre = (6.02d23*1.d5*1.d-3)/(8.3144d0*293d0)

        R_Asdex = 8314.4d0
        mnum = 20.2d0
        kst = 1.666d0 
        ! kst= 5/3 for noble gas as Ne;
        ! A_Dmv = PI*0.7*0.7*1d-4
        ! V_Dmv = 80.0d-6

        yy = (2.d0/(1+kst))**(1.d0/(kst-1))*(2.d0*kst/(kst+1)*R_Asdex*293.d0/mnum)**0.5d0
    
        gam = A_Dmv/V_Dmv*yy
    
        dt_open = 1.5d-3

        if (t_loc .lt. dt_open) then

          prof_temp = - exp(-t_loc*t_loc/2.d0/dt_open*gam)

          ns_dNinj_dt = - prof_temp*t_loc*V_Dmv*1.d3*gam*P_Dmv*N_barlitre/dt_open ! Number of injected particles per unit time (not normalised)

        else

          prof_temp = - exp(-(t_loc-dt_open)*gam)*exp(-dt_open/(2*gam))

          ns_dNinj_dt = - prof_temp*gam*V_Dmv*1.d3*P_Dmv*N_barlitre ! Number of injected particles per unit time (not normalised)
    
        endif

        ns_drhon_dt =  ns_dNinj_dt * mass_gas ! Mass density injected per unit time
    
        ! Inverse of the number of particles still in the reservoir, formulae given by G. Pautasso (ASDEX-U)

        rhon_source = (MU_ZERO)**(0.5d0)*(central_mass*MASS_PROTON*central_density*1.d20)**(-0.5d0)*ns_drhon_dt * ns_pol_shape  * ns_tor_shape / V_ns

        ! Converting mass density into number density
        rhon_source = rhon_source * (central_mass * MASS_PROTON / mass_gas)

      else 

        rhon_source = ns_amplitude * ns_pol_shape * ns_tor_shape * t_norm &
                      /  (V_ns * 1.d20 * central_density)

      endif

    else

      rhon_source = 0.

    endif

  if (rhon_source < 0.) then
    rhon_source = 0.
  end if


  return
  end subroutine inj_source

  subroutine total_imp_source(R,Z,phi,source_background,source_impurity,mass_ratio) 

    use phys_module, only: using_spi, JET_MGI, ASDEX_MGI, n_spi_tot, pellets, ng_radius_ratio, ns_radius
    use phys_module, only: ng_radius_min, n_inj, n_spi, n_spi_tot, ns_deltaphi, L_tube
    use phys_module, only: ns_tor_norm, A_Dmv,K_Dmv,V_Dmv,P_Dmv,t_ns, t_now, central_density, central_mass
    use phys_module, only: ns_amplitude, ns_R, ns_Z, ns_phi

    implicit none

    real*8, intent(in)   :: R
    real*8, intent(in)   :: Z
    real*8, intent(in)   :: phi
    real*8, intent(out)  :: source_background
    real*8, intent(out)  :: source_impurity
    real*8, intent(in)   :: mass_ratio

    ! Temporary variables serving the SPI module
    integer    :: spi_i, i_inj,  n_spi_tmp
    
    real*8     :: spi_R_tmp
    real*8     :: spi_Z_tmp
    real*8     :: spi_phi_tmp
    real*8     :: spi_abl_tmp
    real*8     :: ng_radius !< Radius of neutral gas cloud as a result of the ablation
    real*8     :: source_tmp

    source_background = 0.d0
    source_impurity   = 0.d0

    if (using_spi) then

      if (JET_MGI .or. ASDEX_MGI) then
        write(*,*) "WARNING: Using SPI, disabling MGI settings"
        JET_MGI = .false.
        ASDEX_MGI = .false.
      end if

      do spi_i=1, n_spi_tot

        source_tmp = 0.d0

        if (pellets(spi_i)%spi_radius > 0.0) then
          spi_R_tmp   = pellets(spi_i)%spi_R
          spi_Z_tmp   = pellets(spi_i)%spi_Z
          spi_phi_tmp = pellets(spi_i)%spi_phi
          spi_abl_tmp = pellets(spi_i)%spi_abl

          ng_radius   = pellets(spi_i)%spi_radius * ng_radius_ratio

          if (ng_radius < ng_radius_min) then
            ng_radius = ng_radius_min
          end if
          
          n_spi_tmp = 0
          do i_inj = 1, n_inj
            n_spi_tmp = n_spi_tmp + n_spi(i_inj)
            if (spi_i <= n_spi_tmp)  exit !< Determine the injection location index of the fragment
          end do

          call inj_source(spi_abl_tmp,spi_R_tmp,spi_Z_tmp,spi_phi_tmp,ng_radius,ns_deltaphi,&
                        ns_tor_norm, A_Dmv,K_Dmv,V_Dmv,P_Dmv,t_ns(i_inj),0., R, Z,    &
                        phi,source_tmp,t_now,JET_MGI,ASDEX_MGI,central_density,central_mass)
        end if

        ! Converting number density into mass density for each species respectively
        source_background  = source_background + source_tmp * ( 1. - pellets(spi_i)%spi_species)
        source_impurity    = source_impurity + source_tmp * pellets(spi_i)%spi_species / mass_ratio

      end do

    else

      do i_inj = 1, n_inj
        source_tmp = 0.d0
        call inj_source(ns_amplitude(i_inj),ns_R(i_inj),ns_Z(i_inj),ns_phi(i_inj),   &
                        ns_radius,ns_deltaphi,ns_tor_norm, &
                        A_Dmv,K_Dmv,V_Dmv,P_Dmv,t_ns(i_inj),L_tube,R,Z,phi,source_impurity,&
                        t_now, JET_MGI,ASDEX_MGI,central_density,central_mass)

        source_impurity = source_impurity + source_tmp
      end do

      ! Converting number density into mass density for each species respectively
      source_impurity = source_impurity / mass_ratio

    end if

  end subroutine total_imp_source
end module mod_injection_source
