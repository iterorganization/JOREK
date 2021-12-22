!> Implements a localized neutral source, for example from MGI or a (shattered) pellet
module mod_neutral_source

  use constants

  implicit none

  real*8, save :: total_n_particles_inj     = 0.d0
  real*8, save :: total_n_particles         = 0.d0
  real*8, save :: total_n_particles_inj_all = 0.d0

  contains 



  !> Calculates the neutral source
  subroutine neutral_source(ns_amplitude,ns_R,ns_Z,ns_phi,ns_psi,ns_grad_psi,ns_radius,ns_deltaphi,ns_deltaminrad,ns_tor_norm, &
                              A_Dmv,K_Dmv,V_Dmv,P_Dmv,t_ns,L_tube,R,Z,phi,psi,rhon_source,t_now,               &
                              JET_MGI,ASDEX_MGI,central_density,central_mass,source_volume)

    use mod_source_shape, only: source_shape

    implicit none

    ! --- Routine parameters
    real*8,  intent(in)  :: R, Z, phi, psi, A_Dmv, K_Dmv, V_Dmv, P_Dmv, t_now, t_ns, ns_amplitude
    real*8,  intent(in)  :: ns_R, ns_Z, ns_phi, ns_psi, ns_grad_psi, ns_radius, ns_deltaphi, ns_deltaminrad, L_tube
    real*8,  intent(in)  :: central_density, central_mass, ns_tor_norm
    logical, intent(in)  :: JET_MGI, ASDEX_MGI
    real*8,  intent(out) :: rhon_source
    real*8,  intent(in)  :: source_volume ! numerically integrated gas source volume (if larger than 0.)

    ! --- Local variables
    real*8  :: c0_D, ns_shape, V_ns, f_Nbar, f_dNbar_dt
    real*8  :: ns_dNinj_dt, ns_drhon_dt, t_loc, t_norm, prof_temp, R_Asdex, mnum, kst, yy, gam
    real*8  :: dt_open, N_barlitre, DMV_inj_frac
    integer :: k

    c0_D = sqrt(8.3145d0*293.d0/4.d-3*(7.d0/5.d0))  ! Sound speed of Deuterium

    ! Compute the source shape
    ns_shape = source_shape(R,Z,phi,ns_R,ns_Z,ns_phi,ns_radius,ns_deltaphi,&
         psi,ns_psi,ns_grad_psi,ns_deltaminrad)

    ! Volume used for normalization:
    ! if finite, the input value for source_volume will be used as this will correspond to the numerically integrated gas source volume
    ! otherwise, the analytical value corresponding to the integration in space of the product of the above shape function will be used.
    ! In the standard case with circular ablation cloud in the poloidal plane,
    ! the agreement between the two is very good unless the shard is just marginally inside the domain:
    ! in this case the numerical integral will be smaller than the analytical one, and the resulting total source will correctly reflect the ablation rate (although the local source will be overestimated)
    if (source_volume .gt. 0.) then ! i.e., when numerical integration of the ablation source volume is used
       V_ns = source_volume
    else ! i.e., when numerical integration of the ablation source volume is not used
       if (ns_deltaminrad .gt. 0.) then
    ! i.e., with poloidally elongated ablation cloud
    ! in this case the analytical formula below is approximate (usually it agrees with the numerical integral within a few percents)
          V_ns  = PI * ns_R * ns_tor_norm * ns_radius * min(ns_deltaminrad,ns_radius)
       else
    ! i.e., standard case with circular ablation cloud in the poloidal plane
    ! in this case the ablation source volume is given by the exact analytical formula as derived by E. Nardon
          V_ns  = PI * ns_R * ns_tor_norm * ns_radius**2.d0
       endif
    endif
    ! ===================================================================

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
        rhon_source = (MU_ZERO)**(0.5d0)*(central_mass*MASS_PROTON*central_density*1.d20)**(-0.5d0) * ns_drhon_dt * ns_shape / V_ns

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
        rhon_source = (MU_ZERO)**(0.5d0)*(central_mass*MASS_PROTON*central_density*1.d20)**(-0.5d0) * ns_drhon_dt * ns_shape / V_ns

      else 

        rhon_source =  ns_amplitude * ns_shape * t_norm / (V_ns * 1.d20 * central_density)  

      endif

    else ! t_loc <= 0.
      rhon_source = 0.
    endif

    if (rhon_source < 0.) then
      write(*,*) 'PROBLEM: Negative neutral source!'
    end if

    return
  end subroutine neutral_source


  subroutine total_neutral_source(R,Z,phi,psi,source_neutral) 

    use phys_module, only: using_spi, JET_MGI, ASDEX_MGI, n_spi_tot, pellets, ng_radius_ratio, ns_radius
    use phys_module, only: ng_radius_min, n_inj, n_spi, n_spi_tot, ns_deltaphi, L_tube
    use phys_module, only: ns_tor_norm, A_Dmv,K_Dmv,V_Dmv,P_Dmv,t_ns, t_now, central_density, central_mass
    use phys_module, only: ns_amplitude, ns_R, ns_Z, ns_phi
    use phys_module, only: spi_num_vol, ns_deltaminrad

    implicit none

    real*8, intent(in)   :: R
    real*8, intent(in)   :: Z
    real*8, intent(in)   :: phi
    real*8, intent(in)   :: psi
    real*8, intent(out)  :: source_neutral

    ! Temporary variables serving the SPI module
    integer    :: spi_i, i_inj,  n_spi_tmp
    
    real*8     :: spi_R_tmp
    real*8     :: spi_Z_tmp
    real*8     :: spi_phi_tmp
    real*8     :: spi_abl_tmp
    real*8     :: spi_psi_tmp
    real*8     :: spi_grad_psi_tmp
    real*8     :: ng_radius !< Radius of neutral gas cloud as a result of the ablation
    real*8     :: source_neutral_tmp
    real*8     :: spi_vol_tmp !< Numerically integrated gas source volume

    if (using_spi) then

      do spi_i=1, n_spi_tot

        source_neutral_tmp = 0.d0 

        if (pellets(spi_i)%spi_radius > 0.0) then

          if (spi_num_vol) then
             spi_vol_tmp = pellets(spi_i)%spi_vol
          else
             spi_vol_tmp = 0.d0
          endif

          ng_radius   = pellets(spi_i)%spi_radius * ng_radius_ratio

          if (ng_radius < ng_radius_min) then
            ng_radius = ng_radius_min
          end if

          n_spi_tmp = 0
          do i_inj = 1, n_inj
            n_spi_tmp = n_spi_tmp + n_spi(i_inj)
            if (spi_i <= n_spi_tmp)  exit !< Determine the injection location index of the fragment
          end do

          call neutral_source(pellets(spi_i)%spi_abl,pellets(spi_i)%spi_R,pellets(spi_i)%spi_Z,pellets(spi_i)%spi_phi, &
                        pellets(spi_i)%spi_psi,pellets(spi_i)%spi_grad_psi, &
                        ng_radius,ns_deltaphi,ns_deltaminrad,ns_tor_norm, &
                        A_Dmv,K_Dmv,V_Dmv,P_Dmv,t_ns(i_inj),0.,R,Z,phi,psi, &
                        source_neutral_tmp,t_now,JET_MGI,ASDEX_MGI,central_density,central_mass,spi_vol_tmp)

        end if

        source_neutral = source_neutral + source_neutral_tmp

      end do

    else

      do i_inj = 1, n_inj
        source_neutral_tmp = 0.d0
        spi_vol_tmp = 0.d0
        spi_psi_tmp = 0.d0
        spi_grad_psi_tmp = 0.d0

        call neutral_source(ns_amplitude(i_inj),ns_R(i_inj),ns_Z(i_inj),ns_phi(i_inj),spi_psi_tmp,spi_grad_psi_tmp, &
                      ns_radius,ns_deltaphi,ns_deltaminrad,ns_tor_norm, &
                      A_Dmv,K_Dmv,V_Dmv,P_Dmv,t_ns(i_inj),L_tube,R,Z,phi,psi, &
                      source_neutral_tmp,t_now,JET_MGI,ASDEX_MGI,central_density,central_mass,spi_vol_tmp)
        source_neutral = source_neutral + source_neutral_tmp
      end do
    end if

  end subroutine total_neutral_source

  !> Calculates the total number of neutral particles injected from the start of the simulation and for each timestep.
  subroutine total_neutrals(my_id,node_list,element_list)

    use data_structure
    use phys_module
    use mpi_mod

    implicit none

    ! --- Routine parameters
    type (type_node_list),    intent(in) :: node_list
    type (type_element_list), intent(in) :: element_list
    integer,                  intent(in) :: my_id 

    ! --- Local variables
    integer :: ierr
    real*8  :: density, density_in, density_out, pressure, pressure_in,pressure_out
    real*8  :: kin_par_tot, kin_par_in, kin_par_out, mom_par_tot, mom_par_in, mom_par_out

    call Integrals_3D(my_id,node_list,element_list,density,density_in,density_out,pressure,pressure_in,pressure_out, &
                      kin_par_tot, kin_par_in, kin_par_out, mom_par_tot, mom_par_in, mom_par_out)

    total_n_particles_inj_all = total_n_particles_inj_all + total_n_particles_inj*tstep*sqrt(MU_ZERO*central_mass*MASS_PROTON*central_density*1.d20)

    if (my_id .eq. 0) then
      write(*,'(A,e14.6)') 'total neutrals particles injected per second = '                       , total_n_particles_inj
      write(*,'(A,e14.6)') 'total neutrals particles in the plasma       = '                       , total_n_particles
      write(*,'(A,e14.6)') 'total neutrals particles injected since the start of the simulation = ', total_n_particles_inj_all
    endif

  end subroutine total_neutrals

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
