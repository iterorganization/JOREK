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



  subroutine inj_source(ns_amplitude,ns_R,ns_Z,ns_phi,ns_psi,ns_grad_psi,ns_radius,ns_deltaphi,&
                        ns_delta_minor_rad,ns_tor_norm,  &
                        A_Dmv,K_Dmv,V_Dmv,P_Dmv,t_ns,L_tube,R,Z,phi,psi,rhon_source,t_now,                  &
                        JET_MGI,ASDEX_MGI,central_density,central_mass,source_volume, i_main_imp)

  !=================================================================================
  !  This subroutine computes the atom/ion number density source for a realistic Deuterium
  !  MGI in JET (if ns_timedependent is .t.).
  !  If ns_timedependent is .f., this routine computes a constant source in time
  !  where the main parameter is ns_amplitude
  !  More details in the JOREK wiki or by asking A.Fil or E.Nardon
  !=================================================================================

    use phys_module, only: imp_type
    use mod_source_shape, only: source_shape

    implicit none

    real*8 :: c0_gas                   ! Sound velocity of gas in reservoir
    integer:: n_gas                    ! = 2/(gamma-1) where gamma = heat capacity ratio of gas
    real*8 :: A_gas                    ! Atomic number of gas particles
    real*8 :: mass_gas                 ! Mass of a gas particles
    real*8 :: mol_atom                 ! Number of atoms in a molecular
    real*8 :: ns_shape
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
    real*8, intent(in)  :: psi
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
    real*8, intent(in)  :: ns_psi
    real*8, intent(in)  :: ns_grad_psi
    real*8, intent(in)  :: ns_radius
    real*8, intent(in)  :: ns_deltaphi
    real*8, intent(in)  :: ns_delta_minor_rad
    real*8, intent(in)  :: L_tube
    real*8, intent(in)  :: central_density
    real*8, intent(in)  :: central_mass
    real*8              :: DMV_inj_frac
    logical, intent(in) :: JET_MGI
    logical, intent(in) :: ASDEX_MGI
    real*8, intent(out) :: rhon_source  ! This is in number density
    real*8, intent(in)  :: ns_tor_norm
    real*8, intent(in)  :: source_volume ! numerically integrated gas source volume (if larger than 0.)
    integer, intent(in) :: i_main_imp

    select case ( trim(imp_type(i_main_imp)) )
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
        write(*,*) '!! Gas type "', trim(imp_type(i_main_imp)), '" unknown (in mod_injection_source.f90) !!'
        write(*,*) '=> We assume the gas is D2.'
        n_gas  = 5
        A_gas  = 4.
        mol_atom = 2.
        mass_gas = A_gas*MASS_PROTON
        c0_gas = sqrt(8.3145d0*293.d0/(A_gas*1.d-3)*(7.d0/5.d0))
    end select

    ! ===================================================================
    ! Parameters related to the spatial distribution of the gas source:

    ! Compute the source shape
    ns_shape = source_shape(R,Z,phi,ns_R,ns_Z,ns_phi,ns_radius,ns_deltaphi,&
         psi,ns_psi,ns_grad_psi,ns_delta_minor_rad)

    ! Volume used for normalization:
    ! if finite, the input value for source_volume will be used as this will correspond to the numerically integrated gas source volume
    ! otherwise, the analytical value corresponding to the integration in space of the product of the above shape function will be used
    ! In the standard case with circular ablation cloud in the poloidal plane,
    ! the agreement between the two is very good unless the shard is just marginally inside the domain:
    ! in this case the numerical integral will be smaller than the analytical one, and the resulting total source will correctly reflect the ablation rate (although the local source will be overestimated)
    if (source_volume .gt. 0.) then ! i.e., when numerical integration of the ablation source volume is used
       V_ns = source_volume
    else ! i.e., when numerical integration of the ablation source volume is not used
       if (ns_delta_minor_rad .gt. 0.) then
    ! i.e., with poloidally elongated ablation cloud
    ! in this case the analytical formula below is approximate (usually it agrees with the numerical integral within a few percents)
          V_ns  = PI * ns_R * ns_tor_norm * ns_radius * min(ns_delta_minor_rad,ns_radius)
       else
    ! i.e., standard case with circular ablation cloud in the poloidal plane
    ! in this case the ablation source volume is given by the exact analytical formula as derived by E. Nardon
      V_ns  = PI * ns_R * ns_tor_norm * ns_radius**2.d0
       endif
    endif
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
        rhon_source = ns_drhon_dt * ns_shape / V_ns

        ! Apply JOREK normalization
        rhon_source = (MU_ZERO)**(0.5d0)*(central_mass*MASS_PROTON*central_density*1.d20)**(-0.5d0) * rhon_source

        ! Converting mass density into number density
        rhon_source = rhon_source * (central_mass * MASS_PROTON / mass_gas)
      elseif (ASDEX_MGI) then

        N_barlitre = (6.02d23*1.d5*1.d-3)/(8.3144d0*293d0)

        R_Asdex = 8314.4d0
        mnum = A_gas
        kst = 1.666d0 
        ! kst= 5/3 for noble gas as Ne;
        ! A_Dmv = PI*0.7*0.7*1d-4
        ! V_Dmv = 80.0d-6

        yy = (2.d0/(1+kst))**(1.d0/(kst-1))*(2.d0*kst/(kst+1)*R_Asdex*293.d0/mnum)**0.5d0
    
        gam = A_Dmv/V_Dmv*yy
    
        dt_open = 1.0d-3

        if (t_loc .lt. dt_open) then

          prof_temp = - exp(-t_loc*t_loc/2.d0/dt_open*gam)

          ns_dNinj_dt = - prof_temp*t_loc*V_Dmv*1.d3*gam*P_Dmv*N_barlitre/dt_open ! Number of injected particles per unit time (not normalised)

        else

          prof_temp = - exp(-(t_loc-dt_open)*gam)*exp(-dt_open/(2*gam))

          ns_dNinj_dt = - prof_temp*gam*V_Dmv*1.d3*P_Dmv*N_barlitre ! Number of injected particles per unit time (not normalised)
    
        endif

        ns_drhon_dt =  ns_dNinj_dt * mass_gas ! Mass density injected per unit time
    
        ! Inverse of the number of particles still in the reservoir, formulae given by G. Pautasso (ASDEX-U)

        rhon_source = (MU_ZERO)**(0.5d0)*(central_mass*MASS_PROTON*central_density*1.d20)**(-0.5d0)*ns_drhon_dt * ns_shape / V_ns

        ! Converting mass density into number density
        rhon_source = rhon_source * (central_mass * MASS_PROTON / mass_gas)

      else 

        rhon_source = ns_amplitude * ns_shape * t_norm &
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

  subroutine total_imp_source(R,Z,phi,psi,source_background_arr,source_impurity_arr,mass_ratio,i_main_imp, source_background_drift_arr)

    use phys_module, only: using_spi, JET_MGI, ASDEX_MGI, n_spi_tot, pellets, ns_radius_ratio, ns_radius
    use phys_module, only: ns_radius_min, n_inj_max, n_inj, n_spi, n_spi_tot, ns_deltaphi, L_tube
    use phys_module, only: ns_tor_norm, A_Dmv,K_Dmv,V_Dmv,P_Dmv,t_ns, t_now, central_density, central_mass
    use phys_module, only: ns_amplitude, ns_R, ns_Z, ns_phi
    use phys_module, only: spi_num_vol, ns_delta_minor_rad
    use phys_module, only: drift_distance

    implicit none

    real*8, intent(in)   :: R
    real*8, intent(in)   :: Z
    real*8, intent(in)   :: phi
    real*8, intent(in)   :: psi
    real*8, intent(out)  :: source_background_arr(n_inj_max)
    real*8, intent(out)  :: source_impurity_arr(n_inj_max)
    real*8, intent(out), optional :: source_background_drift_arr(n_inj_max)
    real*8, intent(in)   :: mass_ratio
    integer, intent(in)  :: i_main_imp

    ! Temporary variables serving the SPI module
    integer    :: spi_i, i_inj, n_spi_tmp, n_spi_begin, i
    real*8     :: ns_radius_loc    
    real*8     :: spi_psi_tmp
    real*8     :: spi_grad_psi_tmp
    real*8     :: source_tmp, source_tmp_drift
    real*8     :: spi_vol_tmp !< Numerically integrated gas source volume
    real*8     :: spi_vol_tmp_drift

    source_background_arr = 0.d0
    source_impurity_arr   = 0.d0
    if (present(source_background_drift_arr)) then
      source_background_drift_arr = 0.d0
    end if

    if (using_spi) then

      if (JET_MGI .or. ASDEX_MGI) then
        write(*,*) "WARNING: Using SPI, disabling MGI settings"
        JET_MGI = .false.
        ASDEX_MGI = .false.
      end if

      n_spi_begin = 1
      do i_inj = 1,n_inj

        do i = 1,n_spi(i_inj)
          spi_i = n_spi_begin + i - 1

          source_tmp = 0.d0
          source_tmp_drift = 0.d0

          if (pellets(spi_i)%spi_radius > 0.0) then

            if (spi_num_vol) then
               spi_vol_tmp = pellets(spi_i)%spi_vol
               spi_vol_tmp_drift = pellets(spi_i)%spi_vol_drift
            else
               spi_vol_tmp = 0.d0
               spi_vol_tmp_drift = 0.d0
            endif

            ns_radius_loc   = pellets(spi_i)%spi_radius * ns_radius_ratio

            if (ns_radius_loc < ns_radius_min) then
              ns_radius_loc = ns_radius_min
            end if

            call inj_source(pellets(spi_i)%spi_abl,pellets(spi_i)%spi_R,pellets(spi_i)%spi_Z,pellets(spi_i)%spi_phi, &
                          pellets(spi_i)%spi_psi,pellets(spi_i)%spi_grad_psi, &
                          ns_radius_loc,ns_deltaphi,ns_delta_minor_rad,ns_tor_norm, &
                          A_Dmv,K_Dmv,V_Dmv,P_Dmv,t_ns(i_inj),0., R, Z, phi, psi, &
                          source_tmp,t_now,JET_MGI,ASDEX_MGI,central_density,central_mass,spi_vol_tmp,i_main_imp)

            if (present(source_background_drift_arr)) then
              if (drift_distance(i_inj) /= 0.d0) then
                if (pellets(spi_i)%plasmoid_in_domain == 1) then
                  call inj_source(pellets(spi_i)%spi_abl,pellets(spi_i)%spi_R+drift_distance(i_inj),pellets(spi_i)%spi_Z,pellets(spi_i)%spi_phi, &
                              pellets(spi_i)%spi_psi_drift,pellets(spi_i)%spi_grad_psi_drift, &
                              ns_radius_loc,ns_deltaphi,ns_delta_minor_rad,ns_tor_norm, &
                              A_Dmv,K_Dmv,V_Dmv,P_Dmv,t_ns(i_inj),0., R, Z, phi, psi, &
                              source_tmp_drift,t_now,JET_MGI,ASDEX_MGI,central_density,central_mass,spi_vol_tmp_drift,i_main_imp)

                  if (pellets(spi_i)%spi_species > 0.) then
                    write(*,*) 'WARNING: Do you really want to put <plasmoid teleportation> to non-pure D SPI?'
                  end if 

                else
                  source_tmp_drift = 0.d0
                end if
              else
                source_tmp_drift = source_tmp
              end if 
              source_background_drift_arr(i_inj) = source_background_drift_arr(i_inj) + source_tmp_drift * (1. - pellets(spi_i)%spi_species)
            end if

          end if

          ! Converting number density into mass density for each species respectively
          source_background_arr(i_inj)  = source_background_arr(i_inj) + source_tmp * ( 1. - pellets(spi_i)%spi_species)
          source_impurity_arr(i_inj)    = source_impurity_arr(i_inj) + source_tmp * pellets(spi_i)%spi_species / mass_ratio

        end do

        n_spi_begin = n_spi_begin + n_spi(i_inj)

      end do

    else

      do i_inj = 1, n_inj
        source_tmp = 0.d0
        spi_vol_tmp = 0.d0
        spi_psi_tmp = 0.d0
        spi_grad_psi_tmp = 0.d0

        call inj_source(ns_amplitude(i_inj),ns_R(i_inj),ns_Z(i_inj),ns_phi(i_inj),spi_psi_tmp,spi_grad_psi_tmp, &
                        ns_radius,ns_deltaphi,ns_delta_minor_rad,ns_tor_norm, &
                        A_Dmv,K_Dmv,V_Dmv,P_Dmv,t_ns(i_inj),L_tube,R,Z,phi,psi, &
                        source_tmp, t_now, JET_MGI,ASDEX_MGI,central_density,central_mass, spi_vol_tmp,i_main_imp)

        source_impurity_arr(i_inj) = source_impurity_arr(i_inj) + source_tmp
      end do

      ! Converting number density into mass density for each species respectively
      source_impurity_arr = source_impurity_arr / mass_ratio

    end if

  end subroutine total_imp_source
end module mod_injection_source
