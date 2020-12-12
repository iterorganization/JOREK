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



  subroutine inj_source(ns_amplitude,ns_R,ns_Z,ns_phi,ns_radius,ns_sig,ns_deltaphi,ns_tor_norm,  &
                        A_Dmv,K_Dmv,V_Dmv,P_Dmv,t_ns,L_tube,R,Z,phi,rhon_source,t_now,                  &
                        JET_MGI,ASDEX_MGI,central_density,central_mass)

  !=================================================================================
  !  This subroutine computes the atom/ion number density source for a realistic Deuterium
  !  MGI in JET (if ns_timedependent is .t.).
  !  If ns_timedependent is .f., this routine computes a constant source in time
  !  where the main parameter is ns_amplitude
  !  More details in the JOREK wiki or by asking A.Fil or E.Nardon
  !=================================================================================

    use phys_module, only: gas_type

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
    real*8, intent(in)  :: ns_sig
    real*8, intent(in)  :: ns_deltaphi
    real*8, intent(in)  :: L_tube
    real*8, intent(in)  :: central_density
    real*8, intent(in)  :: central_mass
    real*8              :: DMV_inj_frac
    logical, intent(in) :: JET_MGI
    logical, intent(in) :: ASDEX_MGI
    real*8, intent(out) :: rhon_source  ! This is in number desntiy
    real*8, intent(in)  :: ns_tor_norm

    select case ( trim(gas_type) )
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
        write(*,*) '!! Gas type "', trim(gas_type), '" unknown (in mod_injection_source.f90) !!'
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

  subroutine init_imp_adas(my_id)

    use phys_module
    use mod_openadas
    use mod_coronal

    implicit none

    integer, intent(in) :: my_id
    integer             :: err_alloc, i

    character(len=512)  :: adas_suffix     !The suffix of adas data file to be read

    ! Temporary variable for charge state distribution
    integer             :: i_T, i_ion
    real*8, allocatable :: dP_imp_dT(:), P_imp(:)
    real*8              :: Z_imp

    n_adas = 1 ! For now we only trace one species, in the future probably more 

    if (allocated(imp_adas)) then
      deallocate(imp_adas)
    end if

    allocate (imp_adas(n_adas),stat=err_alloc)  !< Dynamically allocate memeries for adas data

    if (err_alloc /= 0) then
      write(*,*) "Error when trying to dynamically allocate memeries for adas data.", my_id
      stop
    else
      if (allocated(imp_cor)) then
        deallocate(imp_cor)
      end if

      allocate (imp_cor(n_adas),stat=err_alloc)  !< Dynamically allocate memeries for adas data
      if (err_alloc /= 0) then
        write(*,*) "Error when trying to dynamically allocate memeries for CE vector.", my_id
        deallocate(imp_adas)
        stop
      else
        do i=1, n_adas
          select case ( trim(gas_type) )
            case('D2')
              adas_suffix = '12_h'
            case('Ar')
              adas_suffix = '89_ar'
            case('Ne')
              adas_suffix = '96_ne'
            case default
              write(*,*) "Unrecognized species, terminating."
              adas_suffix = 'none'
              deallocate(imp_cor)
              deallocate(imp_adas)
              stop
          end select

          imp_adas(i) = read_adf11(my_id, trim(adas_suffix),trim(adas_dir))
          imp_cor(i)  = coronal(imp_adas(i))

          
          ! This is to output a coronal equilibrium charge distribution as a
          ! function of temperature assuming constant density
          if (my_id == 0) call output_coronal(imp_cor(i))
        end do
      end if

    end if

    if (allocated(xtime_radiation)) call tr_deallocate(xtime_radiation,"xtime_radiation",CAT_GRID)
    if (nstep .gt. 0) call tr_allocate(xtime_radiation,1,nstep,"xtime_radiation")
    if (allocated(xtime_rad_power)) call tr_deallocate(xtime_rad_power,"xtime_rad_power",CAT_GRID)
    if (nstep .gt. 0) call tr_allocate(xtime_rad_power,1,nstep,"xtime_rad_power")
    if (allocated(xtime_E_ion)) call tr_deallocate(xtime_E_ion,"xtime_E_ion",CAT_GRID)
    if (nstep .gt. 0) call tr_allocate(xtime_E_ion,1,nstep,"xtime_E_ion")
    if (allocated(xtime_E_ion_power)) call tr_deallocate(xtime_E_ion_power,"xtime_E_ion_power",CAT_GRID)
    if (nstep .gt. 0) call tr_allocate(xtime_E_ion_power,1,nstep,"xtime_E_ion_power")

  end subroutine init_imp_adas

  subroutine radiation_function(ad,cor, density, temperature, Lrad, dLrad_dTe, dLrad_dNe)

    use phys_module
    use mod_openadas
    use mod_coronal
    use mod_interp_splinear

    implicit none

    type(adf11_all), intent(in) :: ad
    type(coronal), intent(in)   :: cor
    real*8, intent(in)          :: density !< log10 density in m^-3
    real*8, intent(in)          :: temperature !< log10 electron temperature in K

    real*8, intent(out)         :: Lrad ! value of radiation function
    real*8, intent(out), optional :: dLrad_dTe, dLrad_dNe ! derivatives of radiation functioni

    real*8                      :: rad!Local density multiplied radiation function
    real*8                      :: radRB, dradRB_dT, dradRB_dn, radLT, dradLT_dT, dradLT_dn
    real*8, dimension(0:ad%n_Z) :: rad_p, drad_dT, drad_dn
    real*8, dimension(0:cor%n_Z):: p          !< charge state distribution
    real*8, dimension(0:cor%n_Z):: p_Te, p_Ne !< gradient of distribution of charge states (sum = 1) to Te and Ne
    integer :: iz
    
    call cor%interp(density,temperature,rad_out=rad)
    Lrad = rad / (10.0**density) ! This is to recover the radiation coefficient
    if (present(dLrad_dTe) .or. present(dLrad_dNe)) then
      call cor%interp(density,temperature,p_out=p,p_Te_out=p_Te,p_Ne_out=p_Ne)
      do iz=0,ad%n_Z
        call ad%PRB%interp(iz,density,temperature,GRC_out=radRB,dGRC_dT_out=dradRB_dT,dGRC_dn_out=dradRB_dn)
        call ad%PLT%interp(iz,density,temperature,GRC_out=radLT,dGRC_dT_out=dradLT_dT,dGRC_dn_out=dradLT_dn)
        rad_p(iz)   = radRB + radLT
        drad_dT(iz) = dradRB_dT + dradLT_dT
        drad_dn(iz) = dradRB_dn + dradLT_dn
      enddo ! radiation emitted by atoms at level iz
       if (present(dLrad_dTe)) dLrad_dTe = dot_product(p_Te,rad_p) + dot_product(p,drad_dT)
       if (present(dLrad_dNe)) dLrad_dTe = dot_product(p_Ne,rad_p) + dot_product(p,drad_dn)
    end if

  end subroutine radiation_function

  subroutine radiation_function_linear(ad,cor, density, temperature, Lrad, dLrad_dTe)

    use phys_module
    use mod_openadas
    use mod_coronal
    use mod_interp_splinear

    implicit none

    type(adf11_all), intent(in) :: ad
    type(coronal), intent(in)   :: cor
    real*8, intent(in)          :: density !< log10 density in m^-3
    real*8, intent(in)          :: temperature !< log10 electron temperature in K

    real*8, intent(out)         :: Lrad ! value of radiation function
    real*8, intent(out), optional :: dLrad_dTe ! derivatives of radiation functioni

    real*8                      :: rad!Local density multiplied radiation function
    real*8                      :: radRB, radLT
    real*8, dimension(0:ad%n_Z) :: rad_p, drad_dT, dradRB_dT, dradLT_dT
    real*8, dimension(0:cor%n_Z):: p          !< charge state distribution
    real*8, dimension(0:cor%n_Z):: p_Te       !< gradient of distribution of charge states (sum = 1) to Te and Ne
    integer :: iz
    
    call cor%interp_linear(density,temperature,rad_out=rad)
    Lrad = rad / (10.0**density) ! This is to recover the radiation coefficient
    if (present(dLrad_dTe)) then
      call cor%interp_linear(density,temperature,p_out=p,p_Te_out=p_Te)
      dradRB_dT = ad%PRB%interp_grad_T(density,temperature) !Loglog gradient still!!!
      dradLT_dT = ad%PLT%interp_grad_T(density,temperature) !Loglog gradient still!!!
      do iz=0,ad%n_Z
!        radRB     = ad%PRB%interp_linear(iz,density,temperature)
!        radLT     = ad%PLT%interp_linear(iz,density,temperature)
        call ad%PRB%interp_linear(iz,density,temperature,radRB)
        call ad%PLT%interp_linear(iz,density,temperature,radLT)
        rad_p(iz)   = radRB + radLT
        drad_dT(iz) = dradRB_dT(iz) * radRB / (10.0**temperature) &
                      + dradLT_dT(iz) * radLT / (10.0**temperature) ! Convert to normal gradient
      enddo ! radiation emitted by atoms at level iz
      if (present(dLrad_dTe)) dLrad_dTe = dot_product(p_Te,rad_p) + dot_product(p,drad_dT)
    end if

  end subroutine radiation_function_linear

end module mod_injection_source
