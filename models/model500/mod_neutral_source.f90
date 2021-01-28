!> Implements a localized neutral source, for example from MGI or a (shattered) pellet
module mod_neutral_source

  use constants

  implicit none

  real*8, save :: total_n_particles_inj     = 0.d0
  real*8, save :: total_n_particles         = 0.d0
  real*8, save :: total_n_particles_inj_all = 0.d0

  contains 



  !> Calculates the neutral source
  subroutine neutral_source(ns_amplitude,ns_R,ns_Z,ns_phi,ns_radius,ns_sig,ns_deltaphi,ns_tor_norm, &
                              A_Dmv,K_Dmv,V_Dmv,P_Dmv,t_ns,L_tube,R,Z,phi,rhon_source,t_now,               &
                              JET_MGI,ASDEX_MGI,central_density,central_mass)

    implicit none

    ! --- Routine parameters
    real*8,  intent(in)  :: R, Z, phi, A_Dmv, K_Dmv, V_Dmv, P_Dmv, t_now, t_ns, ns_amplitude
    real*8,  intent(in)  :: ns_R, ns_Z, ns_phi, ns_radius, ns_sig, ns_deltaphi, L_tube
    real*8,  intent(in)  :: central_density, central_mass, ns_tor_norm
    logical, intent(in)  :: JET_MGI, ASDEX_MGI
    real*8,  intent(out) :: rhon_source

    ! --- Local variables
    real*8  :: c0_D, radius, ns_tor_shape, ns_pol_shape, dphi, V_ns, f_Nbar, f_dNbar_dt
    real*8  :: ns_dNinj_dt, ns_drhon_dt, t_loc, t_norm, prof_temp, R_Asdex, mnum, kst, yy, gam
    real*8  :: dt_open, N_barlitre, DMV_inj_frac
    integer :: k

    radius = sqrt((R-ns_R)**2 + (Z-ns_Z)**2)

    c0_D = sqrt(8.3145d0*293.d0/4.d-3*(7.d0/5.d0))  ! Sound speed of Deuterium

    ! Poloidal gaussian shape factor
    ns_pol_shape = exp(-(radius/ns_radius)**2.d0)  

    ! Toroidal gaussian shape factor
    dphi = abs(phi - ns_phi)
    if (dphi .gt. PI) dphi = 2*PI - dphi  
    ns_tor_shape = exp(-(dphi/ns_deltaphi)**2.d0)

    V_ns  = PI * ns_R * ns_tor_norm * ns_radius**2.d0

    t_norm = sqrt(MU_ZERO * central_mass * MASS_PROTON * central_density * 1.d20) ! Time normalization factor

    t_loc = (t_now-t_ns) * t_norm

    if (t_loc .ge. 0.) then

      if (JET_MGI) then

      !! We use here the formulae derived from the eq(8) in the paper of S.A. Bozhenkov - NF 51 (2011) 
      !! which gives the normalized number of particles injected at the exit of the DMV injection tube
      !! as a function of time.
      !!
      !! The parameters used are realistic:
      !! A_Dmv: cross sectional area of the injection pipe
      !! K_Dmv: Experimental correction factor to account for the gas expansion close to the tube orifice
      !! L_tube: DMV vacuum injection tube length
      !! V_Dmv: Volume of the DMV reservoir
      !! P_Dmv: Initial pressure in the DMV reservoir, directly linked to the total number of particles
      !! in the reservoir. Expressed in bar here as it is in all MGI experiments. 

        ! Shifted t_loc so that the neutral source is turned on as soon as t_now = t_ns 
        ! (L_tube/3c0 is the time needed for the gas to propagate in the injection tube)
        t_loc = t_loc + L_tube/(3.d0 * c0_D)

        f_Nbar = 0.d0
        f_dNbar_dt = 0.d0

        ! Calculation of the normalized number of particles injected per unit time
        do k = 0,6
          f_Nbar = f_Nbar + (-1.d0)**(k-1)*factorial(6)/(factorial(5-k+1)*factorial(k))*(1-(5.d0*c0_D*t_loc/L_tube)**(1-k))

          f_dNbar_dt = f_dNbar_dt &
                       + (-1.d0)**(k-1)*factorial(6)/(factorial(5-k+1)*factorial(k))*(k-1) &
                         *(5.d0*c0_D*(L_tube)**(-1.d0))**(1-k)*t_loc**(-k)
        end do

        DMV_inj_frac = A_Dmv * L_tube * K_Dmv * f_Nbar/(V_Dmv)

        ! if (DMV_inj_frac .gt. 1.d0) then  
        !   f_dNbar_dt = 0.d0   ! The gas injection is stopped when the initial number of particles in the reservoir is reached  
        ! endif

        ns_dNinj_dt = A_Dmv * K_Dmv * L_tube / V_Dmv * (5.d0)**(5.d0) * (6.d0)**(-6.d0) * f_dNbar_dt   ! Normalised number of injected particles per unit time

        ns_drhon_dt = ns_dNinj_dt * (P_Dmv * 1.d5/(K_BOLTZ * 293)) * V_Dmv * 2.d0 * central_mass * MASS_PROTON ! Mass density per unit time (SI units)
    
        ! Apply gaussian shape (toroidally and poloidally) factor (normalized so that the number of particles injected does not depend on the shape)
        ! as well as JOREK normalization
        rhon_source = (MU_ZERO)**(0.5d0)*(central_mass*MASS_PROTON*central_density*1.d20)**(-0.5d0) * ns_drhon_dt * ns_pol_shape * ns_tor_shape / V_ns

      elseif (ASDEX_MGI) then

        N_barlitre = (6.02d23*1.d5*1.d-3)/(8.3144d0*293d0)
        R_Asdex    = 8314.4d0
        mnum       = 20.2d0
        kst        = 1.666d0 ! kst= 5/3 for noble gas like Ne

        yy = (2.d0/(1+kst))**(1.d0/(kst-1))*(2.d0*kst/(kst+1)*R_Asdex*293.d0/mnum)**0.5d0
    
        gam = A_Dmv/V_Dmv*yy
    
        dt_open = 1.5d-3

        if (t_loc .lt. dt_open) then
          prof_temp    = - exp(-t_loc*t_loc/2.d0/dt_open*gam)
          ns_dNinj_dt = - prof_temp*t_loc*V_Dmv*1.d3*gam*P_Dmv*N_barlitre/dt_open ! Number of particles injected per unit time
        else
          prof_temp    = - exp(-(t_loc-dt_open)*gam)*exp(-dt_open/(2*gam))
          ns_dNinj_dt = - prof_temp*gam*V_Dmv*1.d3*P_Dmv*N_barlitre               ! Number of particles injected per unit time
        endif

        ns_drhon_dt =  ns_dNinj_dt * central_mass * MASS_PROTON ! Mass density injected per unit time

        ! Apply JOREK normalization
        rhon_source = (MU_ZERO)**(0.5d0)*(central_mass*MASS_PROTON*central_density*1.d20)**(-0.5d0) * ns_drhon_dt * ns_pol_shape * ns_tor_shape / V_ns

      else 

        rhon_source =  ns_amplitude * ns_pol_shape * ns_tor_shape * t_norm / (V_ns * 1.d20 * central_density)  

      endif

    else ! t_loc <= 0.
      rhon_source = 0.
    endif

    if (rhon_source < 0.) then
      write(*,*) 'PROBLEM: Negative neutral source!'
    end if

    return
  end subroutine neutral_source

!============================================================!
! Important note: in order to implementing more complicated  !
!    model, we should add more arguments to inj_source       !
!============================================================!
  subroutine get_source(R,Z,phi,source_neutral) 

    use phys_module, only: using_spi, JET_MGI, ASDEX_MGI, n_spi_tot, pellets, ng_radius_ratio, ns_radius
    use phys_module, only: ng_radius_min, n_inj, n_spi, n_spi_tot, ns_sig, ns_deltaphi, L_tube
    use phys_module, only: ns_tor_norm, A_Dmv,K_Dmv,V_Dmv,P_Dmv,t_ns, t_now, central_density, central_mass
    use phys_module, only: ns_amplitude, ns_R, ns_Z, ns_phi

    implicit none

    real*8, intent(in)   :: R
    real*8, intent(in)   :: Z
    real*8, intent(in)   :: phi
    real*8, intent(out)  :: source_neutral

    ! Temporary variables serving the SPI module
    integer    :: spi_i, i_inj,  n_spi_tmp
    
    real*8     :: spi_R_tmp
    real*8     :: spi_Z_tmp
    real*8     :: spi_phi_tmp
    real*8     :: spi_abl_tmp
    real*8     :: ng_radius !< Radius of neutral gas cloud as a result of the ablation
    real*8     :: source_neutral_tmp

    if (using_spi) then

      do spi_i=1, n_spi_tot

        source_neutral_tmp = 0.d0 

        if (pellets(spi_i)%spi_radius > 0.0) then

          ng_radius   = pellets(spi_i)%spi_radius * ng_radius_ratio

          if (ng_radius < ng_radius_min) then
            ng_radius = ng_radius_min
          end if

          n_spi_tmp = 0
          do i_inj = 1, n_inj
            n_spi_tmp = n_spi_tmp + n_spi(i_inj)
            if (spi_i <= n_spi_tmp)  exit !< Determine the injection location index of the fragment
          end do

          call neutral_source(pellets(spi_i)%spi_abl,pellets(spi_i)%spi_R,pellets(spi_i)%spi_Z,pellets(spi_i)%spi_phi,&
                        ng_radius,ns_sig,ns_deltaphi,&
                        ns_tor_norm, A_Dmv,K_Dmv,V_Dmv,P_Dmv,t_ns(i_inj),0.,R,Z,     &
                        phi,source_neutral_tmp,t_now,JET_MGI,ASDEX_MGI,central_density,central_mass)

        end if

        source_neutral = source_neutral + source_neutral_tmp

      end do

    else

      do i_inj = 1, n_inj
        source_neutral_tmp = 0.d0
        call neutral_source(ns_amplitude(i_inj),ns_R(i_inj),ns_Z(i_inj),ns_phi(i_inj), &
                      ns_radius,ns_sig,ns_deltaphi,ns_tor_norm, &
                      A_Dmv,K_Dmv,V_Dmv,P_Dmv,t_ns(i_inj),L_tube,R,Z,phi,source_neutral,t_now, &
                      JET_MGI,ASDEX_MGI,central_density,central_mass)
        source_neutral = source_neutral + source_neutral_tmp
      end do
    end if

  end subroutine get_source

  !> Initialize ADAS for background impurity 
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
        if (nimp_bg .gt. 0) then
          do i=1, n_adas
            select case ( trim(imp_type) )
              case('C')
                adas_suffix = '96_c'
              case('Ar')
                adas_suffix = '89_ar'
              case('Ne')
                adas_suffix = '96_ne'
              case('W')
                adas_suffix = '50_w'
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
    end if

    if ( my_id == 0 ) then
      if (allocated(xtime_radiation)) call tr_deallocate(xtime_radiation,"xtime_radiation",CAT_GRID)
      if (nstep .gt. 0) call tr_allocate(xtime_radiation,1,nstep,"xtime_radiation")
      if (allocated(xtime_rad_power)) call tr_deallocate(xtime_rad_power,"xtime_rad_power",CAT_GRID)
      if (nstep .gt. 0) call tr_allocate(xtime_rad_power,1,nstep,"xtime_rad_power")
      if (allocated(xtime_E_ion)) call tr_deallocate(xtime_E_ion,"xtime_E_ion",CAT_GRID)
      if (nstep .gt. 0) call tr_allocate(xtime_E_ion,1,nstep,"xtime_E_ion")
      if (allocated(xtime_E_ion_power)) call tr_deallocate(xtime_E_ion_power,"xtime_E_ion_power",CAT_GRID)
      if (nstep .gt. 0) call tr_allocate(xtime_E_ion_power,1,nstep,"xtime_E_ion_power")
    end if 

  end subroutine init_imp_adas

  !> Get the radiation coefficients from adas, outputs Lrad in [wm^3] and dLrad_dTe in [Wm^3/K]
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
    real*8, intent(out), optional :: dLrad_dTe ! derivatives of radiation function

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
        radRB     = ad%PRB%interp_linear(iz,density,temperature)
        radLT     = ad%PLT%interp_linear(iz,density,temperature)
        rad_p(iz)   = radRB + radLT
        drad_dT(iz) = dradRB_dT(iz) * radRB / (10.0**temperature) &
                      + dradLT_dT(iz) * radLT / (10.0**temperature) ! Convert to normal gradient
      enddo ! radiation emitted by atoms at level iz
      if (present(dLrad_dTe)) dLrad_dTe = dot_product(p_Te,rad_p) + dot_product(p,drad_dT)
    end if

  end subroutine radiation_function_linear


  !> Calculates the total number of neutral particles injected from the start of the simulation and for each timestep.
  !> Calculates the factorial of a number (which appears in gas dynamics formulae!)
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

end module mod_neutral_source
