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

    call Integrals_3D(my_id,node_list,element_list,density,density_in,density_out,pressure,pressure_in,pressure_out)

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
