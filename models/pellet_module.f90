module pellet_module

  use constants
  use data_structure
  use phys_module
  
  
  real*8 :: total_pellet_particles    !< the (total) pellet particles added in this timestep
  real*8 :: total_plasma_particles    !< the total plasma density (before this timestep)
  real*8 :: total_pellet_volume       !< the volume of the simulated pellet in this timestep
  
  real*8 :: phys_pellet_volume        !< the physical pellet volume (in m^3)
  real*8 :: pellet_volume             !< approximated value of simulated pellet volume
  real*8 :: pellet_atomic             !< atomic number of pellet mass
  
  real*8 :: phys_ablation             !< physical ablation rate (not normalised)
  
  real*8, allocatable  :: xtime_pellet_R(:)
  real*8, allocatable  :: xtime_pellet_Z(:)
  real*8, allocatable  :: xtime_pellet_psi(:)
  real*8, allocatable  :: xtime_pellet_particles(:)
  real*8, allocatable  :: xtime_phys_ablation(:)

  contains

  subroutine pellet_source2(pellet_amplitude,pellet_R,pellet_Z,pellet_psi,pellet_phi, &
                            pellet_radius, pellet_delta_psi, pellet_sig, pellet_length, pellet_ellipse, pellet_theta, &
                            R,Z,psi,phi, r0, T0, central_density, pellet_particles, pellet_density, pellet_volume, &
                            particle_source, volume_source)
  
  implicit none

    real*8, intent(in)  :: R, Z, psi                 ! position where the particle source is calculated
    real*8, intent(inout) :: phi                     ! toroidal position
    real*8, intent(in)  :: T0, r0                    ! local temperature and mass density (JOREK normalised)
    real*8, intent(in)  :: central_density           !< central plasma density (in units 10^20 m^-3)
    real*8, intent(inout)  :: pellet_particles       !< total number of particles in the pellet
    real*8, intent(in)  :: pellet_density            !< pellet density (units 10^20 m^-3)
    real*8, intent(in)  :: pellet_amplitude          !< amplitude of paricle source (when not using ablation model)
    real*8, intent(in)  :: pellet_R, pellet_Z        !< position of the pellet in R,Z
    real*8, intent(in)  :: pellet_phi                !< toroidal position of the pellet
    real*8, intent(in)  :: pellet_radius             !< pellet size (radius) in poloidal plane
    real*8, intent(in)  :: pellet_sig, pellet_length !< sigmas of pellet source in poloidal and toroidal direction
    real*8, intent(in)  :: pellet_ellipse            !< ellipticity of the pellet source in the poloidal plane
    real*8, intent(in)  :: pellet_theta              !< orientation of the pellet ellipse in the poloidal plane
    real*8, intent(in)  :: pellet_psi, pellet_delta_psi
    real*8, intent(in)  :: pellet_volume
    real*8, intent(out) :: particle_source           !< particle source (JOREK normalised units)
    real*8, intent(out) :: volume_source             !< volume of the pellet source (variable used to integrate total pellet volume)
    
    ! --- local variables
    real*8  :: radius, atn, atn_psi, atn_phi, atomic_mass, ablation_rate

#if _OPENMP >= 201511
  !$omp declare simd uniform(pellet_amplitude,pellet_R,pellet_Z,pellet_psi,pellet_phi,              &
  !$omp pellet_radius, pellet_delta_psi, pellet_sig, pellet_length, pellet_ellipse, pellet_theta,   &
  !$omp R, Z, r0, T0, central_density, pellet_particles, pellet_density, pellet_volume)
#endif
    
    particle_source = 0.d0
    volume_source   = 0.d0
    
    if (pellet_amplitude .gt. 0.) then             ! use the fixed source pellet model 
    
    pellet_particles = 0.0
    
      if (phi .gt. PI) phi = 2*PI - phi
    
      radius = sqrt(  (cos(pellet_theta)**2 + 1./pellet_ellipse**2 * sin(pellet_theta)**2)*(R-pellet_R)**2  &
                    + (sin(pellet_theta)**2 + 1./pellet_ellipse**2 * cos(pellet_theta)**2)*(Z-pellet_Z)**2  &
                    + 2.*(R-pellet_R)*(Z-pellet_Z)*sin(pellet_theta)*cos(pellet_theta) * (1./pellet_ellipse**2 - 1.) )
    
      atn     = (0.5d0 - 0.5d0*tanh((radius - pellet_radius)/pellet_sig))
      atn_psi = (0.5d0 - 0.5d0*tanh(abs(psi- pellet_psi)/pellet_delta_psi))
      atn_phi = (0.5d0 - 0.5d0*tanh((phi- pellet_phi)/pellet_length))
    
      particle_source = pellet_amplitude * atn * atn_phi * atn_psi
    
    !S.F. modified here for introducing moving pellet...
    else if (pellet_particles .gt. 0.) then
    
      if (phi .gt. PI) phi = 2*PI - phi
    
      radius = sqrt(  (cos(pellet_theta)**2 + 1./pellet_ellipse**2 * sin(pellet_theta)**2)*(R-pellet_R)**2  &
                    + (sin(pellet_theta)**2 + 1./pellet_ellipse**2 * cos(pellet_theta)**2)*(Z-pellet_Z)**2  &
                    + 2.*(R-pellet_R)*(Z-pellet_Z)*sin(pellet_theta)*cos(pellet_theta) * (1./pellet_ellipse**2 - 1.) )
    
      atn     = (0.5d0 - 0.5d0*tanh((radius - pellet_radius)/pellet_sig))
      atn_phi = (0.5d0 - 0.5d0*tanh((phi- pellet_phi)/pellet_length))
      atn_psi = (0.5d0 - 0.5d0*tanh(abs(psi- pellet_psi)/pellet_delta_psi))
    
    !  pellet_volume = PI * pellet_radius**2 * pellet_R * pellet_phi ! simulated pellet volume
    
      phys_pellet_volume = pellet_particles /pellet_density         ! physical pellet volume 
    
    ! the number of particles ablated from the physical pellet in units 10^20 m^-3 per unit of JOREK time
    
      pellet_atomic = 2.d0
    
    !----------------- model from Kuteev (see Polevoi, PPCF2008)
     ! ablation_rate = 1.62d5 * central_density**(-0.77) * T0**(1.72) * r0**(0.45) * phys_pellet_volume**(0.48) * pellet_atomic**0.217
    !----------------- NGS model Parks (see Gal NF2008)
      ablation_rate = 2.01d4 * central_density**(-0.81) * max(T0,0.d0)**(1.64) * max(r0,0.d0)**(0.33) * phys_pellet_volume**(0.44) * pellet_atomic**0.5
    
    ! particle source in JOREK normalisation
    
      particle_source = ablation_rate / central_density  * atn * atn_phi / pellet_volume
    
      volume_source   = atn * atn_phi
    
      if(volume_source .lt. 0.d0) then
    !    print*, 'volume_source is negative. volume_source=', volume_source
      endif
    
      if(ablation_rate .lt. 0.d0) then
    !    print*, 'ablation rate is negative. step.'
      endif
      
    
    end if
  
  return
  end subroutine pellet_source2

  !> Update the pellet position and  size of the simulated and physical pellet
  !! (from the integral of the pellet particle source)
  subroutine update_pellet(my_id,node_List,element_list)
  
    use constants
    use data_structure
    use phys_module
    use mpi_mod
    use mod_interp, only: interp
    !use mod_integrals3D
    
    implicit none
    
    integer,                      intent(in)    :: my_id
    type (type_node_list),        intent(in)    :: node_list
    type (type_element_list),     intent(in)    :: element_list
    
    real*8  :: psi_axis, psi_bnd
    integer             :: ierr, i_elm, ifail
    real*8  :: V_normalisation, density, density_in, density_out, pressure, pressure_in, pressure_out
    real*8  :: kin_par_tot, kin_par_in, kin_par_out, mom_par_tot, mom_par_in, mom_par_out

    
    real*8  :: R_out, Z_out, s_out, t_out, P0_s,P0_t,P0_st,P0_ss,P0_tt
    
    if (pellet_amplitude .gt. 0) return
    
          call Integrals_3D(my_id, node_list,element_list,density,density_in,density_out,pressure,pressure_in,pressure_out, &
                                                          kin_par_tot, kin_par_in, kin_par_out, mom_par_tot, mom_par_in, mom_par_out)
                                                          
    V_normalisation = 1.d0 / sqrt(central_density * 1d20 * mass_proton * central_mass * MU_ZERO)
    
    pellet_R = pellet_R + pellet_velocity_R * tstep / V_normalisation
    pellet_Z = pellet_Z + pellet_velocity_Z * tstep / V_normalisation
    
    
    phys_ablation = total_pellet_particles * central_density / sqrt(central_density * 1d20 * mass_proton * central_mass * MU_ZERO)
    
    total_pellet_particles = total_pellet_particles * central_density * tstep 
    total_plasma_particles = total_plasma_particles * central_density          ! undo normalisation
    
    
    call find_RZ(node_list,element_list,pellet_R,pellet_Z,R_out,Z_out,i_elm,s_out,t_out,ifail)
    call interp(node_list,element_list,i_elm,1,1,s_out,t_out,pellet_psi,P0_s,P0_t,P0_st,P0_ss,P0_tt)
    
    if (my_id .eq. 0) then
    
        pellet_particles = max(pellet_particles - total_pellet_particles, 0.d0)
    
        write(*,'(A,4e14.6)') ' pellet (R,Z) =', pellet_R, pellet_Z,pellet_velocity_R/V_normalisation,pellet_velocity_Z/V_normalisation
        write(*,'(A,6e14.6,A)') ' total particles added in this step : ', pellet_R, pellet_Z,pellet_particles,total_pellet_particles, total_plasma_particles,phys_ablation,' [10^20]'
        write(*,'(A,4e14.6)') ' remaining particles in pellet      : ', pellet_particles
        write(*,'(A,4e14.6)') ' pellet volume (sim,phys)           : ', total_pellet_volume,pellet_particles/pellet_density
    
    else 
      pellet_particles = 0.0
    end if
     
    call MPI_Bcast(pellet_particles,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    
    return 
   
  end subroutine update_pellet

  !> Update the shattered pellet position and the simulated and physical pellet sizes
  !! (from the integral of the pellet particle source)  
  subroutine update_spi(my_id,node_List,element_list,i_inj,n_spi_begin)

    use constants
    use data_structure
    use phys_module, only: pellets, imp_type, central_density, central_mass, spi_abl_model, spi_tor_rot,      &
                           ns_phi_rotate, tor_frequency, tstep, pellet_density, pellet_density_bg,            &
                           index_now, xtime_spi_ablation, xtime_spi_ablation_bg, xtime_spi_ablation_rate,&
                           xtime_spi_ablation_bg_rate, F0, R_geo, imp_cor, index_main_imp, n_adas, drift_distance
    use mpi_mod
    use corr_neg
    
    implicit none
    
    integer,                  intent(in) :: my_id
    type (type_node_list),    intent(in) :: node_list
    type (type_element_list), intent(in) :: element_list
  
    ! --- Local variables
    real*8  :: V_normalisation, density, density_in, density_out, pressure,pressure_in,pressure_out
    real*8  :: kin_par_tot, kin_par_in, kin_par_out, mom_par_tot, mom_par_in, mom_par_out
    
    real*8  :: R_out, Z_out
    integer :: i_elm, ifail, i, ierr, i_p
    
    real*8, dimension(4) :: P, P_s, P_t, P_phi
    real*8  :: R, R_s, R_t, Z, Z_s, Z_t
    real*8  :: s_out,t_out
    
    real*8  :: n_SI, T_eV, n_corr, T_corr, n_imp_SI, ne_SI
    real*8  :: t_norm, B0, nu
    real*8  :: spi_delta_phi, spi_Vel_R_tmp, spi_Vel_phi_tmp, spi_phi_inj
    real*8  :: spi_density_tmp
    real*8  :: V_ns
    real*8  :: xjac, psi_R, psi_Z
    
    !   -Mean impurity ionization state and related quantities
    real*8  :: Z_imp, beta_imp, mu_imp

    integer, intent(in) :: i_inj
    integer, intent(in) :: n_spi_begin

    ! - Extra variables when considering plasmoid drift by teleportation
    real*8  :: R_out_drift, Z_out_drift
    integer :: i_elm_drift, ifail_drift
    real*8, dimension(4) :: P_drift, P_s_drift, P_t_drift, P_phi_drift
    real*8  :: R_drift, R_s_drift, R_t_drift, Z_drift, Z_s_drift, Z_t_drift
    real*8  :: s_out_drift,t_out_drift
    real*8  :: xjac_drift, psi_R_drift, psi_Z_drift

  
    spi_delta_phi   = 0.
    spi_Vel_R_tmp   = 0.
    spi_Vel_phi_tmp = 0.

    V_normalisation = 1.d0 / sqrt(central_density * 1d20 * mass_proton * central_mass * MU_ZERO)
    t_norm          = sqrt(MU_ZERO * central_mass * MASS_PROTON * central_density * 1.d20)
  
    loop_over_shards: do i=1, n_spi(i_inj)

      i_p = i - 1 + n_spi_begin

      spi_phi_inj              = pellets(i_p)%spi_phi_init + ns_phi_rotate

      if (spi_phi_inj >= 2.*PI) then
        spi_phi_inj   = mod(spi_phi_inj,2.*PI)
      else if (spi_phi_inj < 0.) then
        spi_phi_inj   = mod(spi_phi_inj,2.*PI) + 2.*PI
      end if
      
      spi_delta_phi            = pellets(i_p)%spi_phi - spi_phi_inj
      spi_Vel_R_tmp            = pellets(i_p)%spi_Vel_R * cos(spi_delta_phi) &
                                 + pellets(i_p)%spi_Vel_RxZ * sin(spi_delta_phi)
      spi_Vel_phi_tmp          = pellets(i_p)%spi_Vel_RxZ * cos(spi_delta_phi) &
                                 - pellets(i_p)%spi_Vel_R * sin(spi_delta_phi)
      spi_Vel_phi_tmp          = spi_Vel_phi_tmp / pellets(i_p)%spi_R

      pellets(i_p)%spi_R       = pellets(i_p)%spi_R + spi_Vel_R_tmp * tstep / V_normalisation
      pellets(i_p)%spi_Z       = pellets(i_p)%spi_Z + pellets(i_p)%spi_Vel_Z * tstep / V_normalisation
      pellets(i_p)%spi_phi     = pellets(i_p)%spi_phi + spi_Vel_phi_tmp * tstep / V_normalisation
  
      if (spi_tor_rot) then
        pellets(i_p)%spi_phi   = pellets(i_p)%spi_phi + tor_frequency * 2. * PI * tstep / V_normalisation
      end if
     
      if (pellets(i_p)%spi_phi >= 2.*PI) then
        pellets(i_p)%spi_phi   = mod(pellets(i_p)%spi_phi,2.*PI)
      else if (pellets(i_p)%spi_phi < 0.) then
        pellets(i_p)%spi_phi   = mod(pellets(i_p)%spi_phi,2.*PI) + 2.*PI
      end if
  
      if (pellets(i_p)%spi_radius > 0.0) then
        if (pellets(i_p)%spi_species == 1.) then
          spi_density_tmp = pellet_density
        else if (pellets(i_p)%spi_species == 0.) then
          spi_density_tmp = pellet_density_bg
        else if (pellet_density_bg > 0. .and. pellet_density > 0.) then
          spi_density_tmp = 1./((1.-pellets(i_p)%spi_species)/pellet_density_bg + pellets(i_p)%spi_species/pellet_density)
        else
          write(*,*) "ERROR: Something is wrong when determining the pellet species, exiting"
          stop
        end if
  
        if (spi_density_tmp == 0. .or. spi_density_tmp /= spi_density_tmp) then
          write(*,*) "ERROR: Problem calculating spi_density!", spi_density_tmp
          stop
        endif
  
        pellets(i_p)%spi_radius = pellets(i_p)%spi_radius - t_norm * tstep * &
                                  (pellets(i_p)%spi_abl / (4.d0 * PI * pellets(i_p)%spi_radius**2.d0 *    &
                                  spi_density_tmp * 1.d20))
  
        if (pellets(i_p)%spi_radius < 0.d0) then
          pellets(i_p)%spi_radius = 0.d0
        end if
      end if
  
      if (my_id == 0.) then
        if (index_now > 1) then
          xtime_spi_ablation(i_p,index_now)    = xtime_spi_ablation(i_p,index_now-1) &
                                               + t_norm * tstep * pellets(i_p)%spi_abl * pellets(i_p)%spi_species
          xtime_spi_ablation_bg(i_p,index_now) = xtime_spi_ablation_bg(i_p,index_now-1) &
                                               + t_norm * tstep * pellets(i_p)%spi_abl * (1. - pellets(i_p)%spi_species)
        else
          xtime_spi_ablation(i_p,index_now)    = t_norm * tstep * pellets(i_p)%spi_abl * pellets(i_p)%spi_species
          xtime_spi_ablation_bg(i_p,index_now) = t_norm * tstep * pellets(i_p)%spi_abl * (1. - pellets(i_p)%spi_species) 
        end if
      end if
  
      if (spi_abl_model == 0) then
        pellets(i_p)%spi_abl   = ns_amplitude(i_inj)
      else if (spi_abl_model >= 1) then
  
        call find_RZ(node_list,element_list,pellets(i_p)%spi_R,pellets(i_p)%spi_Z,&
                     R_out,Z_out,i_elm,s_out,t_out,ifail)
  
        ! In case the shards are outside of the domain
        if (ifail == 99 .or. ifail == 999) then
          pellets(i_p)%spi_abl = 0.
          cycle
        else if (ifail /= 0) then
          write(*,*) "Something wrong in find_RZ!! my_id = ", my_id, i_elm, ifail
          stop
        end if

        if (drift_distance /= 0) then ! when considering plasmoid drift by shifting neutral source
          call find_RZ(node_list,element_list,pellets(i_p)%spi_R+drift_distance,pellets(i_p)%spi_Z,&
                           R_out_drift,Z_out_drift,i_elm_drift,s_out_drift,t_out_drift,ifail_drift) 
        end if

#if ((defined WITH_Impurities) || (defined WITH_Neutrals))
#ifdef WITH_TiTe
        call interp_PRZ(node_list,element_list,i_elm,[var_rho,var_Te,var_rhon,1],4,s_out,t_out,pellets(i_p)%spi_phi,&
                        P,P_s,P_t,P_phi,R,R_s,R_t,Z,Z_s,Z_t)        
        if (drift_distance /= 0) then
          call interp_PRZ(node_list,element_list,i_elm_drift,[var_rho,var_Te,var_rhon,1],4,s_out_drift,t_out_drift,pellets(i_p)%spi_phi,&
                                 P_drift,P_s_drift,P_t_drift,P_phi_drift,R_drift,R_s_drift,R_t_drift,Z_drift,Z_s_drift,Z_t_drift)
        end if

#else /* WITH_TiTe */
        call interp_PRZ(node_list,element_list,i_elm,[var_rho,var_T,var_rhon,1],4,s_out,t_out,pellets(i_p)%spi_phi,&
                        P,P_s,P_t,P_phi,R,R_s,R_t,Z,Z_s,Z_t)
        if (drift_distance /= 0) then
          call interp_PRZ(node_list,element_list,i_elm_drift,[var_rho,var_T,var_rhon,1],4,s_out_drift,t_out_drift,pellets(i_p)%spi_phi,&
                                 P_drift,P_s_drift,P_t_drift,P_phi_drift,R_drift,R_s_drift,R_t_drift,Z_drift,Z_s_drift,Z_t_drift)
        end if
#endif /* WITH_TiTe */
#endif /* ((defined WITH_Impurities) || (defined WITH_Neutrals)) */

        xjac  = R_s * Z_t - R_t * Z_s
        psi_R = (  P_s(4) * Z_t - P_t(4) * Z_s ) / xjac
        psi_Z = (- P_s(4) * R_t + P_t(4) * R_s ) / xjac
        pellets(i_p)%spi_psi = P(4)
        pellets(i_p)%spi_grad_psi = sqrt(psi_R**2 + psi_Z**2)

        if (drift_distance /= 0) then
          xjac_drift  = R_s_drift * Z_t_drift - R_t_drift * Z_s_drift
          psi_R_drift = (  P_s_drift(4) * Z_t_drift - P_t_drift(4) * Z_s_drift ) / xjac_drift
          psi_Z_drift = (- P_s_drift(4) * R_t_drift + P_t_drift(4) * R_s_drift ) / xjac_drift
          pellets(i_p)%spi_psi_drift = P_drift(4)
          pellets(i_p)%spi_grad_psi_drift = sqrt(psi_R_drift**2 + psi_Z_drift**2)
        end if

        ! Now, P(1) represents mass density and P(2) represents temperature, P(3)
        ! is the impurity density
        ! Correct any possible negative values!
  
        !n_corr         = corr_neg_dens(P(1))
        !T_corr         = corr_neg_temp(P(2))
  
        ! Reminder, temperature should be divided by 2 since T = T_e + T_i and T_e
        ! = T_i
        !n_SI           = n_corr * 1.d20 * central_density
        !T_eV           = T_corr / (2.d0* EL_CHG * MU_ZERO * central_density * 1.d20)
        
        n_SI           = P(1) * 1.d20 * central_density
        if (n_SI < 0.) n_SI = 0.
  
#ifdef WITH_TiTe
        T_eV           = P(2) / (EL_CHG * MU_ZERO * central_density * 1.d20)
#else
        T_eV           = P(2) / (2.d0* EL_CHG * MU_ZERO * central_density * 1.d20)
#endif
        if (T_eV < 0.) T_eV = 0.
  
        ! This is only used for 501 as impurity density, 500 don't use this
        ! variable
        n_imp_SI           = P(3) * 1.d20 * central_density
        if (n_imp_SI < 0.) n_imp_SI = 0.      
  
        ! NGS model
        if (spi_abl_model == 1) then
          pellets(i_p)%spi_abl    = 4.12d16 * (pellets(i_p)%spi_radius**(4.0/3.0)) * (n_SI**(1.0/3.0)) * &
                                   (T_eV**1.64)
          if (my_id == 0 .and. pellets(i_p)%spi_radius > 0.0 .and. mod(index_now,20)==0) then
            write(*,*) "Check Point, n_SI, T_eV = ", n_SI, T_eV
          end if
        else if (spi_abl_model == 2) then
          select case ( trim(imp_type(index_main_imp)) )
            case('D2')
              ne_SI   = n_SI
              ! The scaling law is in gauss unit
              pellets(i_p)%spi_abl = 3.9d14 * ((pellets(i_p)%spi_radius*1.d2)**1.455) &
                                     * ((ne_SI*1.d-6)**0.455) * (T_eV**1.679)
            case('Ar')
              if (T_eV >= 1.) then
                ! As with element_matrix, mimick density as 1.d20
                call imp_cor(index_main_imp)%interp(density=20.,temperature=log10(T_eV*EL_CHG/K_BOLTZ),z_out=Z_imp)
              else
                Z_imp = 0.
              end if
              mu_imp             = central_mass/40. ! Argon mass = 40 u and main ion mass = central_mass u
              beta_imp           = mu_imp*Z_imp - 1.
              ne_SI              = n_SI + beta_imp * n_imp_SI
  
              if (ne_SI<0.) ne_SI = 0.        
              ! The scaling law is in gauss unit
              if (pellets(i_p)%spi_species == 1.) then
                pellets(i_p)%spi_abl = 2.5d13 * ((pellets(i_p)%spi_radius*1.d2)**1.451) &
                                       * ((ne_SI*1.d-6)**0.451) * (T_eV**1.679)
              else if (pellets(i_p)%spi_species == 0.) then
                pellets(i_p)%spi_abl = 3.9d14 * ((pellets(i_p)%spi_radius*1.d2)**1.455) &
                                       * ((ne_SI*1.d-6)**0.455) * (T_eV**1.679) 
              end if
            ! Using general scaling law of Sergeev for Neon
            case('Ne')
              if (T_eV >= 1.) then
                ! As with element_matrix, mimick density as 1.d20
                call imp_cor(index_main_imp)%interp(density=20.,temperature=log10(T_eV*EL_CHG/K_BOLTZ),z_out=Z_imp)
              else
                Z_imp = 0.
              end if
              mu_imp             = central_mass/20. ! Neon mass = 20 u and main ion mass = central_mass u
              beta_imp           = mu_imp*Z_imp - 1.
              ne_SI              = n_SI + beta_imp * n_imp_SI
  
              if (ne_SI<0.) ne_SI = 0.
              ! The scaling law is in gauss unit
              ! The sublimation energy for Ne is 0.02 eV
              if (pellets(i_p)%spi_species == 1.) then
                pellets(i_p)%spi_abl = 1.94d14 * ((pellets(i_p)%spi_radius*1.d2)**1.44) &
                                       * ((ne_SI*1.d-6)**0.45) * (T_eV**1.72)         &
                                       * (0.02**(-0.16)) * (20.**(-0.28))             &
                                       * (10.**(-0.56)) * ((2./3.)**0.28)
              else if (pellets(i_p)%spi_species == 0.) then
                pellets(i_p)%spi_abl = 3.9d14 * ((pellets(i_p)%spi_radius*1.d2)**1.455) &
                                       * ((ne_SI*1.d-6)**0.455) * (T_eV**1.679)
              end if
            case default
              write(*,*) '!! Gas type "', trim(imp_type(index_main_imp)), '" unknown !!'
              write(*,*) '=> We assume the gas is D2.'
              pellets(i_p)%spi_abl = 3.9d14 * ((pellets(i_p)%spi_radius*1.d2)**1.455) &
                                     * ((n_SI*1.d-6)**0.455) * (T_eV**1.679)
          end select
          if (my_id == 0 .and. pellets(i_p)%spi_radius > 0.0 .and. mod(index_now,20)==0) then
            write(*,*) "Check Point, ne_SI, T_eV = ", ne_SI, T_eV
          end if
        else if (spi_abl_model == 3) then
          select case ( trim(imp_type(index_main_imp)) )
            case('D2') ! We temporarily wusing D2 ablation rate for H2 ablation here
              pellets(i_p)%spi_abl = 39.0023 * 2. * MOLE_NUMBER * ((pellets(i_p)%spi_radius*1.d2 / 0.2)**(4./3.)) &
                                     * ((n_SI*1.d-20)**(1./3.)) * ((T_eV/2.d3)**(5./3.)) / 4.0282
            case('Ar')  ! Argon and H2/D2 formed separately
              if (T_eV >= 1.) then
                ! As with element_matrix, mimick density as 1.d20
                call imp_cor(index_main_imp)%interp(density=20.,temperature=log10(T_eV*EL_CHG/K_BOLTZ),z_out=Z_imp)
              else
                Z_imp = 0.
              end if
              mu_imp             = central_mass/40. ! Argon mass = 40 u and main ion mass = central_mass u
              beta_imp           = mu_imp*Z_imp - 1.
              ne_SI              = n_SI + beta_imp * n_imp_SI
  
              if (ne_SI<0.) ne_SI = 0.
              ! The scaling law is in gauss unit
              if (pellets(i_p)%spi_species == 1.) then
                pellets(i_p)%spi_abl = 36.6337 * MOLE_NUMBER * ((pellets(i_p)%spi_radius*1.d2 / 0.2)**(4./3.)) &
                                       * ((ne_SI*1.d-20)**(1./3.)) * ((T_eV/2.d3)**(5./3.)) / 39.948
              else if (pellets(i_p)%spi_species == 0.) then
                pellets(i_p)%spi_abl = 39.0023 * 2. * MOLE_NUMBER * ((pellets(i_p)%spi_radius*1.d2 / 0.2)**(4./3.)) &
                                       * ((ne_SI*1.d-20)**(1./3.)) * ((T_eV/2.d3)**(5./3.)) / 4.0282
              end if
            case('Ne')  ! Neond and H2/D2 mixed together
              if (T_eV >= 1.) then
                ! As with element_matrix, mimick density as 1.d20
                call imp_cor(index_main_imp)%interp(density=20.,temperature=log10(T_eV*EL_CHG/K_BOLTZ),z_out=Z_imp)
              else
                Z_imp = 0.
              end if
              mu_imp             = central_mass/20. ! Neon mass = 20 u and main ion mass = central_mass u
              beta_imp           = mu_imp*Z_imp - 1.
              ne_SI              = n_SI + beta_imp * n_imp_SI
  
              if (ne_SI<0.) ne_SI = 0.
              ! The scaling law is in gauss unit
              ! The sublimation energy for Ne is 0.02 eV
                pellets(i_p)%spi_abl = (27.0837 + TAN(1.48709*(1.-pellets(i_p)%spi_species)/(1.+pellets(i_p)%spi_species))) &
                                       * MOLE_NUMBER * ((pellets(i_p)%spi_radius*1.d2 / 0.2)**(4./3.)) &
                                       * ((ne_SI*1.d-20)**(1./3.)) * ((T_eV/2.d3)**(5./3.)) &
                                       / (20.183*pellets(i_p)%spi_species + 2.0141*(1.-pellets(i_p)%spi_species)) 
            case default
              write(*,*) '!! Gas type "', trim(imp_type(index_main_imp)), '" unknown !!'
              write(*,*) '=> We assume the gas is D2.'
              pellets(i_p)%spi_abl = 39.0023 * 2. * MOLE_NUMBER * ((pellets(i_p)%spi_radius*1.d2 / 0.2)**(4./3.)) &
                                     * ((n_SI*1.d-20)**(1./3.)) * ((T_eV/2.d3)**(5./3.)) / 4.0282
          end select
  
          B0 = abs(F0 / R_geo)
          nu = 0.843
          if (B0 > 2.) pellets(i_p)%spi_abl = pellets(i_p)%spi_abl * (2./B0)**nu
  
        else 
          write(*,*) "Unknown ablation model, terminating now!"
          stop
        end if
        if (my_id == 0 .and. pellets(i_p)%spi_radius > 0.0 .and. mod(index_now,20)==0) then
          write(*,*) "Check Point, ne_SI, T_eV = ", ne_SI, T_eV
        end if
      else
        pellets(i_p)%spi_abl    = 0.d0
      end if
     
      if (my_id == 0) then
        xtime_spi_ablation_rate(i_p,index_now) = pellets(i_p)%spi_abl * pellets(i_p)%spi_species
        xtime_spi_ablation_bg_rate(i_p,index_now) = pellets(i_p)%spi_abl * (1. - pellets(i_p)%spi_species)
      end if
  
    end do loop_over_shards
  
    if (spi_tor_rot) then
      ns_phi_rotate  = ns_phi_rotate + tor_frequency * 2. * PI * tstep / V_normalisation
    end if
  
  return  
  end subroutine update_spi


  !> This is a blanket subroutine taking care of initialization of fragments for all injection locations
  subroutine init_spi_all()
  
    use data_structure
    use phys_module, only: pellets, n_spi, n_spi_tot, n_inj, JET_MGI, ASDEX_MGI, ns_R, ns_Z, ns_phi,&
                           ns_amplitude, spi_Vel_Rref, spi_Vel_Zref, spi_Vel_RxZref,                &
                           spi_quantity, spi_quantity_bg, spi_Vel_diff, spi_L_inj, spi_L_inj_diff,  &
                           spi_plume_file
    use mpi_mod
    
    implicit none
    
    integer             :: i, n_spi_begin
    logical             :: ferr
    
    n_spi_tot = 0
    do i = 1, n_inj
      n_spi_tot = n_spi_tot + n_spi(i)
    end do

    if (allocated(pellets)) then
      deallocate(pellets)
    end if

    allocate (pellets(n_spi_tot))  !< Dynamically allocate memeries for pellets

    if (JET_MGI .or. ASDEX_MGI) then
      write(*,*) "WARNING: Using SPI, conflicting with MGI settings"
      write(*,*) "JET_MGI:", JET_MGI
      write(*,*) "ASDEX_MGI:", ASDEX_MGI
      stop
    else      !< Do one initialization for each injection location
      n_spi_begin = 1
      do i = 1, n_inj
        if (spi_plume_file(i) /= 'none') then ! if file exists, read shards information from it
          write(*,'(A21,I2,A19,A)') " >> Initialising SPI(", i, ") from a datafile: ", trim(spi_plume_file(i))
          call init_spi_plume_file(i,n_spi(i),n_spi_begin)
          n_spi_begin = n_spi_begin + n_spi(i)
        else ! if file does not exist, initialise shards based on parameters in the JOREK input file
          write(*,'(A21,I2,A23)') " >> Initialising SPI(", i, ") from input parameters"
          call init_spi(i,ns_R(i),ns_Z(i),ns_phi(i),ns_amplitude(i),spi_Vel_Rref(i),spi_Vel_Zref(i),spi_Vel_RxZref(i),&
                        spi_quantity(i),spi_quantity_bg(i),spi_Vel_diff(i),spi_L_inj(i),spi_L_inj_diff(i),n_spi(i),n_spi_begin)
          n_spi_begin = n_spi_begin + n_spi(i)
        end if
      end do
    end if

    return
  end subroutine init_spi_all

  !> Initializes the shattered pellet position, velocity and size
  subroutine init_spi(i_inj,ns_R,ns_Z,ns_phi,ns_amplitude,spi_Vel_Rref,spi_Vel_Zref,spi_Vel_RxZref,&
                      spi_quantity,spi_quantity_bg,spi_Vel_diff,spi_L_inj,spi_L_inj_diff,n_spi,n_spi_begin)
  
    use constants
    use tr_module
    use data_structure
    use phys_module, only: pellets, imp_type, central_density, central_mass, pellet_density, pellet_density_bg,&
                           spi_rnd_seed, spi_angle, xtime_spi_ablation, xtime_spi_ablation_bg, xtime_spi_ablation_rate,&
                           xtime_spi_ablation_bg_rate, nstep, spi_shard_file, spi_abl_model, n_spi_tot, index_main_imp,&
                           n_adas
    use mpi_mod
    use corr_neg
    
    implicit none
    
    integer             :: i, i_p
    logical             :: ferr
    
    real*8  :: spi_gd_angle_01, spi_gd_angle_02        ! The dispersion angles for each shard
    real*8  :: spi_rotation_01, spi_rotation_02        ! The rotation angle from shard coordinates to (R,Z,phi) coordinates
    real*8  :: spi_Vel_totref, spi_Vel_i, spi_Vel_R_tmp, spi_Vel_Z_tmp, spi_Vel_RxZ_tmp, spi_L_inj_i
    real*8  :: spi_Vel_x, spi_Vel_y, spi_Vel_z         ! Shard velocity in injection coordinates
    real*8  :: spi_R_inj, spi_Z_inj, spi_phi_inj       ! Position of the shattering point of the pellet
                                                       ! (the apex of the spreading cone)
    real*8  :: spi_R_tmp, spi_Z_tmp, spi_phi_tmp, spi_density_tmp
    real*8  :: sign_corr
    real*8, allocatable :: rnd(:)                      ! The random number array 
    real*8, allocatable :: shard_size(:)               ! The shard size array

    real*8  :: size_beta                               ! The characteristic shard size    
    real*8  :: N_shard_norm                            ! The normalized (by size_beta) number of atoms
    real*8  :: mix_ratio                               ! Volume mixture ratio of the mixed pellet 
    real*8  :: real_spi_quantity(2)                    ! Final injection quantity

    integer,intent(in)  :: i_inj
    real*8, intent(in)  :: ns_R
    real*8, intent(in)  :: ns_Z
    real*8, intent(in)  :: ns_phi
    real*8, intent(in)  :: ns_amplitude
    real*8, intent(in)  :: spi_Vel_Rref
    real*8, intent(in)  :: spi_Vel_Zref
    real*8, intent(in)  :: spi_Vel_RxZref
    real*8, intent(in)  :: spi_quantity
    real*8, intent(in)  :: spi_quantity_bg
    real*8, intent(in)  :: spi_Vel_diff
    real*8, intent(in)  :: spi_L_inj
    real*8, intent(in)  :: spi_L_inj_diff

    integer, intent(in) :: n_spi
    integer, intent(in) :: n_spi_begin

    if (n_spi >= 1) then
      if (allocated(shard_size)) deallocate(shard_size)
      allocate (shard_size(n_spi))  !< Dynamically allocate memeries for shard sizes
      shard_size = 0.0

      size_beta    = 0.0
      N_shard_norm = 0.0
      real_spi_quantity = 0.

      ! Read normalized shard size distribution (if given in file) and calculate
      ! shard radius normalization factor size_beta
      if (spi_shard_file(i_inj) /= 'none') then 
        inquire(file=trim(spi_shard_file(i_inj)), exist=ferr) ! Check if the file exists
        if (ferr) then
          open(42,file=trim(spi_shard_file(i_inj)),status="OLD",action="READ")
          read(42,*)  shard_size(1:n_spi)
          close(42)
        else
          write(*,*) "ERROR: Shard size file does not exist, exiting now"
          stop
        end if
      else
        shard_size = 1.
      end if

#if (defined WITH_Neutrals) && (!defined WITH_Impurities)
      do i = 1, n_spi
        i_p = i - 1 + n_spi_begin
        pellets(i_p)%spi_species = 0.
        N_shard_norm = N_shard_norm + (4./3.) * PI * (shard_size(i)**3) * pellet_density_bg * 1.d20
      end do

      size_beta    = (spi_quantity_bg / N_shard_norm) ** (-1./3.)
      write(*,*) "Characteristic shard size (m):", 1./size_beta
#endif
#ifdef WITH_Impurities
      ! Determine approximately how many fragments are of the impurity, how
      ! much are of the background species. 

      if (spi_quantity_bg > 0. .and. pellet_density_bg > 0. &
          .and. spi_quantity > 0. .and. pellet_density > 0.) then
        mix_ratio = (spi_quantity/pellet_density)&
                    /((spi_quantity/pellet_density)+(spi_quantity_bg/pellet_density_bg))
      else if (spi_quantity > 0. .and. pellet_density > 0.) then
        mix_ratio = 1.
      else if (spi_quantity_bg > 0. .and. pellet_density_bg > 0.) then
        mix_ratio = 0.
      else
        write(*,*) "WARNING!!! Something is wrong in the injection quantity or pellet density, exiting."
        stop
      end if

      select case ( trim(imp_type(index_main_imp)) ) 
        case('D2')
          write(*,*) "Injection of D2 species should be done by spi_quantity_bg, please revise input file accordingly."
          stop
        case('Ne')
          ! Only Parks formula can properly treat the mixing of neon and D2/H2,
          ! otherwise we assume neon and D2/H2 formed seperately.
          if (spi_abl_model == 3 .and. mix_ratio < 1. .and. mix_ratio > 0.) then
            do i = 1, n_spi
              i_p = i - 1 + n_spi_begin
              pellets(i_p)%spi_species = spi_quantity/(spi_quantity + spi_quantity_bg)
              spi_density_tmp = 1./((1.-pellets(i_p)%spi_species)/pellet_density_bg &
                                    + pellets(i_p)%spi_species/pellet_density) 
              real_spi_quantity(2) = real_spi_quantity(2) &
                                     + (4./3.) * PI * (shard_size(i)**3) * spi_density_tmp * 1.d20 &
                                       * pellets(i_p)%spi_species
              real_spi_quantity(1) = real_spi_quantity(1) &
                                     + (4./3.) * PI * (shard_size(i)**3) * spi_density_tmp * 1.d20 &
                                       * (1. - pellets(i_p)%spi_species)
              N_shard_norm = N_shard_norm + (4./3.) * PI * (shard_size(i)**3) * spi_density_tmp *1.d20
            end do
          else
            do i = 1, n_spi
              i_p = i - 1 + n_spi_begin
              if (i <= int(real(n_spi,8)*(mix_ratio))) then
                pellets(i_p)%spi_species = 1.
                spi_density_tmp = pellet_density
                real_spi_quantity(2) = real_spi_quantity(2) + (4./3.) * PI * (shard_size(i)**3) * spi_density_tmp *1.d20
              else
                pellets(i_p)%spi_species = 0.
                spi_density_tmp = pellet_density_bg
                real_spi_quantity(1) = real_spi_quantity(1) + (4./3.) * PI * (shard_size(i)**3) * spi_density_tmp *1.d20
              end if
              N_shard_norm = N_shard_norm + (4./3.) * PI * (shard_size(i)**3) * spi_density_tmp *1.d20
            end do
          end if
        case('Ar')
          ! Argon and D2/H2 part of the pellet are always formed seperately, thus we always treat them as such.
          do i = 1, n_spi
            i_p = i - 1 + n_spi_begin
            if (i <= int(real(n_spi,8)*(mix_ratio))) then
              pellets(i_p)%spi_species = 1.
              spi_density_tmp = pellet_density
              real_spi_quantity(2) = real_spi_quantity(2) + (4./3.) * PI * (shard_size(i)**3) * spi_density_tmp *1.d20
            else
              pellets(i_p)%spi_species = 0.
              spi_density_tmp = pellet_density_bg
              real_spi_quantity(1) = real_spi_quantity(1) + (4./3.) * PI * (shard_size(i)**3) * spi_density_tmp *1.d20
            end if
            N_shard_norm = N_shard_norm + (4./3.) * PI * (shard_size(i)**3) * spi_density_tmp *1.d20
          end do
        case default
          write(*,*) '!! Gas type "', trim(imp_type(index_main_imp)), '" unknown !!'
          write(*,*) '=> We assume the gas is D2.'
          write(*,*) "Injection of D2 species should be done by spi_quantity_bg, please revise input file accordingly."
          stop
      end select
  
      size_beta    = ((spi_quantity+spi_quantity_bg) / N_shard_norm) ** (-1./3.)
      real_spi_quantity(1) = real_spi_quantity(1) / size_beta**3
      real_spi_quantity(2) = real_spi_quantity(2) / size_beta**3
      write(*,*) "Characteristic shard size (m):", 1./size_beta
      write(*,*) "Real injection quantity (atom):", real_spi_quantity(1), real_spi_quantity(2)
#endif
      ! Initialize shard radius
      do i = 1, n_spi
        i_p = i - 1 + n_spi_begin 
        pellets(i_p)%spi_radius  = shard_size(i)/size_beta !Here we are using exactly the same distribution for each injector
      end do

      !===================Determine the rotational transform of coordinate===============
      !Here, we perform the following rotational transform from the original
      !coordinate R, Z, RxZ to the so-called spi coordinate x, y ,z, with the !reference
      !direction of spi injection being the z axis, while y axis locates within the 
      !same surface as Z and z. The rotational transform from x, y, z to R, Z, RxZ is
      !as the following: first, we rotate the system around x axis clockwise, facing
      !the positive x direction, for spi_rotation_01 to get coordinate X', Y', Z'. 
      !Then we further rotate around Y' clockwise, facing the positive Y' direction !for
      !spi_rotation_02 to acquire R, Z, RxZ. Hence we have:
      !R   = cos(spi_rotation_02)*x - sin(spi_rotation_02)*(-sin(spi_rotation_01)*y + !cos(spi_rotation_01)*z)
      !Z   = cos(spi_rotation_01)*y + sin(spi_rotation_01)*z
      !RxZ = sin(spi_rotation_02)*x + cos(spi_rotation_02)*(-sin(spi_rotation_01)*y + !cos(spi_rotation_01)*z)

      spi_Vel_totref  = sqrt(spi_Vel_Rref**2+spi_Vel_Zref**2+spi_Vel_RxZref**2)

      spi_R_inj       = ns_R - spi_L_inj * (spi_Vel_Rref/spi_Vel_totref)
      spi_Z_inj       = ns_Z - spi_L_inj * (spi_Vel_Zref/spi_Vel_totref)
      spi_phi_inj     = ns_phi - spi_L_inj * (spi_Vel_RxZref/spi_Vel_totref)/ns_R

      spi_rotation_01 = asin(spi_Vel_Zref/spi_Vel_totref)
      if (cos(spi_rotation_01) == 0.) then
        spi_rotation_02 = 0.
      else
        spi_rotation_02 = acos(spi_Vel_RxZref/(spi_Vel_totref*cos(spi_rotation_01)))
      end if

      write(*,*) "Rotational transform: ", spi_rotation_01, spi_rotation_02

!==========================End of rotational angles==============================

      ! Generate a random number array rnd that contains two random angles
      ! representing the velocity direction spread, one random speed, and one random position difference.
      ! Those random numbers uniquely define a random velocity and position of the shard,
      ! which is then transformed into the R, Z, RxZ space.
      if (allocated(rnd)) deallocate(rnd)
      allocate (rnd(4*n_spi))  !< Dynamically allocate memeries for randoms

      CALL random_seed(put=spi_rnd_seed) 
      CALL random_number(rnd)

      !write(*,*) "Random number array:", rnd

      if (spi_Vel_diff < 0) then
        write(*,*) "WARNING, negative velocity spread, spi_Vel_diff = ", spi_Vel_diff
        write(*,*) "Please always use a positive spi_Vel_diff, EXITING!" 
        stop
      end if

      if (spi_L_inj_diff < 0) then
        write(*,*) "WARNING, negative position spread, spi_L_inj_diff = ", spi_L_inj_diff
        write(*,*) "Please always use a positive spi_L_inj_diff, EXITING!" 
        stop
      end if

      do i=1, n_spi

        i_p = i - 1 + n_spi_begin
        spi_gd_angle_01 =   rnd(4 * i - 3)         * spi_angle / 2.0
        spi_gd_angle_02 =   rnd(4 * i - 2)         * 2. * PI
        spi_Vel_i       = ( rnd(4 * i - 1) - 0.5 ) * spi_Vel_diff   + spi_Vel_totref
        spi_L_inj_i     = ( rnd(4 * i    ) - 0.5 ) * spi_L_inj_diff + spi_L_inj


        !write(*,*) "Random angle:", i, spi_gd_angle_01, spi_gd_angle_02

        spi_Vel_x       = spi_Vel_i * sin(spi_gd_angle_01) * cos(spi_gd_angle_02)
        spi_Vel_y       = spi_Vel_i * sin(spi_gd_angle_01) * sin(spi_gd_angle_02)
        spi_Vel_z       = spi_Vel_i * cos(spi_gd_angle_01)

        spi_Vel_R_tmp   = spi_Vel_x * cos(spi_rotation_02) &
                          - sin(spi_rotation_02) * (-sin(spi_rotation_01)*spi_Vel_y &
                          + cos(spi_rotation_01)*spi_Vel_z)

        spi_Vel_Z_tmp   = cos(spi_rotation_01) * spi_Vel_y &
                          + sin(spi_rotation_01) * spi_Vel_z
        spi_Vel_RxZ_tmp = spi_Vel_x * sin(spi_rotation_02) &
                          - cos(spi_rotation_02) * (-sin(spi_rotation_01)*spi_Vel_y &
                          + cos(spi_rotation_01)*spi_Vel_z)

        spi_R_tmp       = spi_R_inj + spi_L_inj_i * (spi_Vel_R_tmp/spi_Vel_totref)
        spi_Z_tmp       = spi_Z_inj + spi_L_inj_i * (spi_Vel_Z_tmp/spi_Vel_totref)
        spi_phi_tmp     = spi_phi_inj + spi_L_inj_i * (spi_Vel_RxZ_tmp/spi_Vel_totref)/ns_R

        pellets(i_p)%spi_R       = spi_R_tmp
        pellets(i_p)%spi_Z       = spi_Z_tmp
        pellets(i_p)%spi_phi     = spi_phi_tmp
        pellets(i_p)%spi_phi_init= spi_phi_inj
        pellets(i_p)%spi_Vel_R   = spi_Vel_R_tmp
        pellets(i_p)%spi_Vel_Z   = spi_Vel_Z_tmp
        pellets(i_p)%spi_Vel_RxZ = spi_Vel_RxZ_tmp
        pellets(i_p)%spi_abl     = 0.0
        pellets(i_p)%spi_vol     = 0.0
        pellets(i_p)%spi_psi     = 0.0
        pellets(i_p)%spi_grad_psi= 0.0
        pellets(i_p)%spi_vol_drift     = 0.0
        pellets(i_p)%spi_psi_drift     = 0.0
        pellets(i_p)%spi_grad_psi_drift= 0.0

        write(*,'(A,I5,5ES10.2)') ' *** SHATTERED PELLET PARAMETERS :',i_p, pellets(i_p)%spi_R, pellets(i_p)%spi_Z, &
                              pellets(i_p)%spi_Vel_R, pellets(i_p)%spi_Vel_Z, pellets(i_p)%spi_radius
      end do

      deallocate(rnd)

      if (allocated(xtime_spi_ablation)) call tr_deallocate(xtime_spi_ablation,"xtime_spi_ablation",CAT_UNKNOWN)
      if (nstep .gt. 0) call tr_allocate(xtime_spi_ablation,1,n_spi_tot,1,nstep,"xtime_spi_ablation")

      if (allocated(xtime_spi_ablation_rate)) &
      call tr_deallocate(xtime_spi_ablation_rate,"xtime_spi_ablation_rate",CAT_UNKNOWN)
      if (nstep .gt. 0) call tr_allocate(xtime_spi_ablation_rate,1,n_spi_tot,1,nstep,"xtime_spi_ablation_rate")

      if (allocated(xtime_spi_ablation_bg)) &
      call tr_deallocate(xtime_spi_ablation_bg,"xtime_spi_ablation_bg",CAT_UNKNOWN)
      if (nstep .gt. 0) call tr_allocate(xtime_spi_ablation_bg,1,n_spi_tot,1,nstep,"xtime_spi_ablation_bg")

      if (allocated(xtime_spi_ablation_bg_rate)) &
      call tr_deallocate(xtime_spi_ablation_bg_rate,"xtime_spi_ablation_bg_rate",CAT_UNKNOWN)
      if (nstep .gt. 0) call tr_allocate(xtime_spi_ablation_bg_rate,1,n_spi_tot,1,nstep,"xtime_spi_ablation_bg_rate")

    else
      write(*,*) "ERROR: n_spi<1"
      stop
    end if

    return
  end subroutine init_spi

  !> Initializes the shattered pellet position, velocity and size by a given data file (HDF5 flag needed to be implemented)
  !!  * non-HDF5 version of file is based on the format given by M3D-C1
  !!   - should consists of 8 columns
  !!   - 1st   2nd          3rd   4th         5th           6th         7th        8th 
  !!     R [m] phi[radians] Z [m] Vel_R [m/s] Vel_phi [m/s] Vel_Z [m/s] radius [m] mol(D2)/(mol(D2)+mol(Impurity))
  subroutine init_spi_plume_file(i_inj,n_spi,n_spi_begin)

    use iso_fortran_env

    use constants
    use tr_module
    use data_structure
    use phys_module, only: pellets, imp_type, pellet_density, pellet_density_bg,  xtime_spi_ablation,           &
                           xtime_spi_ablation_bg, xtime_spi_ablation_rate, xtime_spi_ablation_bg_rate, nstep,   &
                           spi_plume_file, spi_plume_hdf5, spi_abl_model, n_spi_tot, n_adas, index_main_imp,    &
                           spi_tor_rot, ns_phi_rotate, tor_frequency
    use mpi_mod
#ifdef USE_HDF5
    use hdf5
    use hdf5_io_module
    use mod_parameters
#endif    

    implicit none

    integer,            intent(in) :: i_inj, n_spi
    integer,            intent(in) :: n_spi_begin

    integer, parameter  :: n_col_expected = 8  ! given non-HDF5 datafile should consists of '8' columns
    real*8, allocatable :: spi_R_tmp(:),     spi_phi_tmp(:),     spi_Z_tmp(:),     &
                           spi_Vel_R_tmp(:), spi_Vel_phi_tmp(:), spi_Vel_Z_tmp(:), &
                           spi_radius_tmp(:),spi_species_molar_D2_tmp(:)

    real*8              :: spi_species_molar_D2_sum

    real*8              :: spi_species_atomic_tmp
    real*8              :: spi_density_tmp
    real*8              :: real_spi_quantity(2)

    integer             :: i, i_p
 
    logical             :: ferr

    ! variables related to check whether the shards data file is in its format
    integer             :: i_line, n_line, io, n_col
    character           :: old_char, new_char
    logical             :: beg_line
#ifdef USE_HDF5
    integer(HID_T)      :: file_id
    integer             :: error
#endif

    write(*,*) " - Following SPI-related input parameters will not play any role during the initialisation "
    write(*,'(2(A18))') "'spi_quantity',","'spi_quantity_bg'"
    write(*,'(3(A18))') "'ns_R'        ,","'ns_Z'          ,","'ns_phi'         "
    write(*,'(4(A18))') "'spi_Vel_Rref',","'spi_Vel_Zref'  ,","'spi_Vel_RxZref',","'spi_Vel_diff'   "
    write(*,'(2(A18))') "'spi_L_inj'   ,","'spi_L_inj_diff' "
    write(*,'(1(A18))') "'spi_angle'    "

    if (n_spi >= 1) then

      ! read spi shards information from Brendan's format
      if ( .not. spi_plume_hdf5) then

        ! check 1) file existence
        inquire(file=trim(spi_plume_file(i_inj)), exist=ferr)
        ! check 2) file format
        if (ferr) then
          open(32, file=trim(spi_plume_file(i_inj)), status='old', action='read', form='unformatted', access='stream')
          old_char = " "
          n_line   = 0
          beg_line = .true.
          n_col    = 0
          do
            read(32, iostat=io) new_char

            ! exit with End-Of-File
            if (io < 0) exit

            ! if empty lines exist in ASCII data file, ignore those lines and give warnings
            if (beg_line .and. new_char == new_line(new_char)) then
              write(*,'(a, i5, a)') " Warning: ", n_line, "th line is empty. It is possible that the spi shard file is defected."

            ! deal with the beginning of line for white spaces
            else if (beg_line) then
              beg_line = .false.

            ! count columns
            else if (new_char == " " .and. old_char /= " ") then
              n_col = n_col+1

            ! end of line
            else if (new_char == new_line(new_char)) then
              if (old_char /= " ") then
                n_col = n_col+1
              end if

              n_line = n_line+1

              if (n_col /= n_col_expected) then
                write(*,*) "ERROR: 'spi_plume_file' is defected."
                stop
              end if

              ! reset for next line
              n_col    = 0
              beg_line = .true.
            end if

            old_char = new_char

          end do

          close(32)

          ! check 3) Consistency of spi shard file with the parameters in the input file (especially 'n_spi')
          if (n_spi /= n_line) then
            write(*,*) "ERROR: 'spi_plume_file' contains different number of lines than the given 'n_spi' in the input file."
            stop
          end if

        else

          write(*,'(A33,I2,A)') "ERROR: 'spi_plume_file' for SPI (", i_inj, ") does not exist, exiting now"
          stop

        end if

      else ! spi_plume_hdf5 == .true.

#ifdef USE_HDF5

        call HDF5_open(trim(spi_plume_file(i_inj)),file_id,error)
        if ( error /= 0 ) then
          write(*,*) "ERROR: failed to open 'spi_plume_file (HDF5)'."
          stop
        end if

        call HDF5_integer_reading(file_id,n_line,"n_spi")

        if (n_spi /= n_line) then
            write(*,*) "ERROR: 'n_spi' value does not match between the 'spi_plume_file (HDF5)' and the input file."
            stop
        end if

#else

        write(*,*) "ERROR: trying to use 'spi_plume_file' in HDF5 format without 'USE_HDF5'"
        stop

#endif

      end if

      ! Now that we have a proper spi shard files (ASCII or HDF5), start to read it
      allocate( spi_R_tmp(n_spi),      spi_phi_tmp(n_spi),      spi_Z_tmp(n_spi),     &
                spi_Vel_R_tmp(n_spi),  spi_Vel_phi_tmp(n_spi),  spi_Vel_Z_tmp(n_spi), &
                spi_radius_tmp(n_spi), spi_species_molar_D2_tmp(n_spi) )

      if ( .not. spi_plume_hdf5) then

        open(42,file=trim(spi_plume_file(i_inj)),status="old",action="read")

        do i = 1,n_line
          read(42, *, iostat=io) spi_R_tmp(i),      spi_phi_tmp(i),      spi_Z_tmp(i),     &
                                 spi_Vel_R_tmp(i),  spi_Vel_phi_tmp(i),  spi_Vel_Z_tmp(i), &
                                 spi_radius_tmp(i), spi_species_molar_D2_tmp(i)
          if (io /= 0) then
            if (io == iostat_end) then
              exit
            else
              write(*,*) "ERROR: in reading 'spi_plume_file', iostat = ", io
              stop
            end if
          end if
        end do

        close(42)

      else ! spi_plume_hdf5 == .true.

#ifdef USE_HDF5

        call HDF5_array1D_reading(file_id, spi_R_tmp,               "spi_R")
        call HDF5_array1D_reading(file_id, spi_phi_tmp,             "spi_phi")
        call HDF5_array1D_reading(file_id, spi_Z_tmp,               "spi_Z")
        call HDF5_array1D_reading(file_id, spi_Vel_R_tmp,           "spi_Vel_R")
        call HDF5_array1D_reading(file_id, spi_Vel_phi_tmp,         "spi_Vel_phi")
        call HDF5_array1D_reading(file_id, spi_Vel_Z_tmp,           "spi_Vel_Z")
        call HDF5_array1D_reading(file_id, spi_radius_tmp,          "spi_radius")
        call HDF5_array1D_reading(file_id, spi_species_molar_D2_tmp,"spi_species_molar_D2")

        call HDF5_close(file_id)

#else

        write(*,*) "ERROR: trying to use 'spi_plume_file' in HDF5 format without 'USE_HDF5'"
        stop

#endif

      end if

      spi_species_molar_D2_sum = 0.d0
      do i = 1, n_spi
        if (spi_species_molar_D2_tmp(i) < 0. .or. spi_species_molar_D2_tmp(i) > 1.) then
          write(*,*) "ERROR: D2 molar fraction in spi data file has illegal values for the fragment ", i
        end if
        spi_species_molar_D2_sum = spi_species_molar_D2_sum + spi_species_molar_D2_tmp(i)
      end do

#if (defined WITH_Neutrals) && (!defined WITH_Impurities)
      do i = 1, n_spi
        i_p = i - 1 + n_spi_begin
        pellets(i_p)%spi_species = 0.
      end do
#endif
#ifdef WITH_Impurities

      real_spi_quantity = 0.d0

      select case ( trim(imp_type(index_main_imp)) )
        case('D2')
          write(*,*) "Injection of D2 species should be done by spi_quantity_bg, please revise input file accordingly."
          stop
        case('Ne')
          ! Only Parks formula can properly treat the mixing of neon and D2/H2.
          ! otherwise we assume neon and D2/H2 formed separately.
          if (spi_abl_model == 3) then
            do i = 1, n_spi
              i_p = i - 1 + n_spi_begin
              spi_species_atomic_tmp = 1.d0 - 2.d0 * spi_species_molar_D2_tmp(i) / (spi_species_molar_D2_tmp(i) + 1.d0)
              pellets(i_p)%spi_species = spi_species_atomic_tmp
              spi_density_tmp = 1./((1.-pellets(i_p)%spi_species)/pellet_density_bg &
                                    + pellets(i_p)%spi_species/pellet_density) 
              real_spi_quantity(2) = real_spi_quantity(2) &
                                     + (4./3.) * PI * (spi_radius_tmp(i)**3) * spi_density_tmp * 1.d20 &
                                       * pellets(i_p)%spi_species
              real_spi_quantity(1) = real_spi_quantity(1) &
                                     + (4./3.) * PI * (spi_radius_tmp(i)**3) * spi_density_tmp * 1.d20 &
                                       * (1. - pellets(i_p)%spi_species)
            end do
          else
            do i = 1, n_spi
              i_p = i - 1 + n_spi_begin
              spi_species_atomic_tmp = 1.d0 - 2.d0 * spi_species_molar_D2_tmp(i) / (spi_species_molar_D2_tmp(i) + 1.d0)
              if (spi_species_atomic_tmp == 1.) then
                pellets(i_p)%spi_species = spi_species_atomic_tmp
                spi_density_tmp = pellet_density
                real_spi_quantity(2) = real_spi_quantity(2) + (4./3.) * PI * (spi_radius_tmp(i)**3) * spi_density_tmp *1.d20
              else if (spi_species_atomic_tmp == 0.) then
                pellets(i_p)%spi_species = spi_species_atomic_tmp
                spi_density_tmp = pellet_density_bg
                real_spi_quantity(1) = real_spi_quantity(1) + (4./3.) * PI * (spi_radius_tmp(i)**3) * spi_density_tmp *1.d20
              else
                write(*,*) "ERROR: Only 'spi_abl_model = 3' can properly treat the mixing of neon and D2/H2  , exiting."
                stop
              end if
            end do
          end if
        case('Ar')
          ! Argon and D2/H2 part of the pellet are always formed separately, thus we always treat them as such.
          do i = 1, n_spi
            i_p = i - 1 + n_spi_begin
            spi_species_atomic_tmp = 1.d0 - 2.d0 * spi_species_molar_D2_tmp(i) / (spi_species_molar_D2_tmp(i) + 1.d0)
            if (spi_species_atomic_tmp == 1.) then
              pellets(i_p)%spi_species = spi_species_atomic_tmp
              spi_density_tmp = pellet_density
              real_spi_quantity(2) = real_spi_quantity(2) + (4./3.) * PI * (spi_radius_tmp(i)**3) * spi_density_tmp *1.d20
            else if (spi_species_atomic_tmp == 0.) then
              pellets(i_p)%spi_species = spi_species_atomic_tmp
              spi_density_tmp = pellet_density_bg
              real_spi_quantity(1) = real_spi_quantity(1) + (4./3.) * PI * (spi_radius_tmp(i)**3) * spi_density_tmp *1.d20
            else
              write(*,*) "ERROR: Argon and D2/H2 part of the pellet are always formed separately, exiting."
              stop
            end if
          end do
        case default
          write(*,*) '!! Gas type "', trim(imp_type(index_main_imp)), '" unknown !!'
          write(*,*) '=> We assume the gas is D2.'
          write(*,*) "Injection of D2 species should be done by spi_quantity_bg, please revise input file accordingly."
          stop
      end select

      write(*,*) "Real injection quantity (atom):", real_spi_quantity(1), real_spi_quantity(2)
#endif

      ! put values into arrays
      do i=1, n_spi

        i_p = i - 1 + n_spi_begin

        ! conversion of coordinates from M3D-C1 (R,phi,Z) to JOREK (R,Z,phi)
        pellets(i_p)%spi_R       =   spi_R_tmp(i)
        pellets(i_p)%spi_Z       =   spi_Z_tmp(i)
        pellets(i_p)%spi_phi     = - spi_phi_tmp(i)
        pellets(i_p)%spi_phi_init= - spi_phi_tmp(i)
        pellets(i_p)%spi_Vel_R   =   spi_Vel_R_tmp(i)
        pellets(i_p)%spi_Vel_Z   =   spi_Vel_Z_tmp(i)
        pellets(i_p)%spi_Vel_RxZ = - spi_Vel_phi_tmp(i)
        pellets(i_p)%spi_radius  =   spi_radius_tmp(i)
        pellets(i_p)%spi_abl     =   0.d0
        pellets(i_p)%spi_vol     =   0.d0
        pellets(i_p)%spi_psi     =   0.d0
        pellets(i_p)%spi_grad_psi=   0.d0
        pellets(i_p)%spi_vol_drift      =   0.d0
        pellets(i_p)%spi_psi_drift      =   0.d0
        pellets(i_p)%spi_grad_psi_drift =   0.d0

        write(*,'(A,I5,5ES10.2)') ' *** SHATTERED PELLET PARAMETERS :',i_p, pellets(i_p)%spi_R, pellets(i_p)%spi_Z, &
                              pellets(i_p)%spi_Vel_R, pellets(i_p)%spi_Vel_Z, pellets(i_p)%spi_radius

      end do

      if (allocated(xtime_spi_ablation)) call tr_deallocate(xtime_spi_ablation,"xtime_spi_ablation",CAT_UNKNOWN)
      if (nstep .gt. 0) call tr_allocate(xtime_spi_ablation,1,n_spi_tot,1,nstep,"xtime_spi_ablation")

      if (allocated(xtime_spi_ablation_rate)) &
      call tr_deallocate(xtime_spi_ablation_rate,"xtime_spi_ablation_rate",CAT_UNKNOWN)
      if (nstep .gt. 0) call tr_allocate(xtime_spi_ablation_rate,1,n_spi_tot,1,nstep,"xtime_spi_ablation_rate")

      if (allocated(xtime_spi_ablation_bg)) &
      call tr_deallocate(xtime_spi_ablation_bg,"xtime_spi_ablation_bg",CAT_UNKNOWN)
      if (nstep .gt. 0) call tr_allocate(xtime_spi_ablation_bg,1,n_spi_tot,1,nstep,"xtime_spi_ablation_bg")

      if (allocated(xtime_spi_ablation_bg_rate)) &
      call tr_deallocate(xtime_spi_ablation_bg_rate,"xtime_spi_ablation_bg_rate",CAT_UNKNOWN)
      if (nstep .gt. 0) call tr_allocate(xtime_spi_ablation_bg_rate,1,n_spi_tot,1,nstep,"xtime_spi_ablation_bg_rate")

    else
      write(*,*) "ERROR: n_spi<1"
      stop
    end if

    return
  end subroutine init_spi_plume_file


  !> This function creates a derived MPI type for the pellets and returns it (in honor of Daan)
  !! If it already exists the old handle is returned
  function get_pellet_derived_type() result(dtype_out)
    use mpi_mod
    use mod_parameters
  
    implicit none
  
    integer               :: ierr, dtype_out
    integer, save         :: dtype
    logical, save         :: dtype_set = .false.
  
    integer :: len(16) = (/1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1/), t(16) = (/ &
      MPI_REAL8,MPI_REAL8,MPI_REAL8,MPI_REAL8,MPI_REAL8, &
      MPI_REAL8,MPI_REAL8,MPI_REAL8,MPI_REAL8,MPI_REAL8, &
      MPI_REAL8,MPI_REAL8,MPI_REAL8,MPI_REAL8,MPI_REAL8, &
      MPI_REAL8/) ! MPI_INTEGER1 == MPI_LOGICAL1
  
    integer(kind=MPI_ADDRESS_KIND) :: base, disp(16)
    type(type_SPI) :: sample_pellet
  
    dtype_out = dtype
    if (dtype_set) return
  
    ! Get memory addresses in the type
    call MPI_Get_address(sample_pellet,             base,    ierr)
    call MPI_Get_address(sample_pellet%spi_R,       disp(1), ierr)
    call MPI_Get_address(sample_pellet%spi_Z,       disp(2), ierr)
    call MPI_Get_address(sample_pellet%spi_phi,     disp(3), ierr)
    call MPI_Get_address(sample_pellet%spi_phi_init,disp(4), ierr)
    call MPI_Get_address(sample_pellet%spi_Vel_R,   disp(5), ierr)
    call MPI_Get_address(sample_pellet%spi_Vel_Z,   disp(6), ierr)
    call MPI_Get_address(sample_pellet%spi_Vel_RxZ, disp(7), ierr)
    call MPI_Get_address(sample_pellet%spi_radius,  disp(8), ierr)
    call MPI_Get_address(sample_pellet%spi_abl,     disp(9), ierr)
    call MPI_Get_address(sample_pellet%spi_species, disp(10),ierr)
    call MPI_Get_address(sample_pellet%spi_vol,     disp(11),ierr)
    call MPI_Get_address(sample_pellet%spi_psi,     disp(12),ierr)
    call MPI_Get_address(sample_pellet%spi_grad_psi,disp(13),ierr)
    call MPI_Get_address(sample_pellet%spi_vol_drift,     disp(14),ierr)
    call MPI_Get_address(sample_pellet%spi_psi_drift,     disp(15),ierr)
    call MPI_Get_address(sample_pellet%spi_grad_psi_drift,disp(16),ierr)
  
    ! Rebase to particle memory beginning
    disp = disp - base
  
    ! Commit the structured type
    call MPI_Type_create_struct(16, len, disp, t, dtype, ierr)
    call MPI_Type_commit(dtype, ierr)
  
    ! Set the save bit
    dtype_set = .true.
    dtype_out = dtype
    return
  end function get_pellet_derived_type



  subroutine pellet_source(pellet_amplitude,pellet_R,pellet_Z,pellet_psi,pellet_phi, &
                           pellet_radius, pellet_delta_psi, pellet_sig, pellet_length, &
                           R,Z,psi,phi,particle_source)
  
    use constants
    
    implicit none
    
    real*8 :: R, Z, psi, phi, particle_source
    real*8 :: pellet_amplitude, pellet_R, pellet_Z, pellet_phi, pellet_radius, pellet_sig, pellet_length
    real*8 :: pellet_psi, pellet_delta_psi, radius, atn, atn_psi, atn_phi
    
    if (phi .gt. PI) phi = 2*PI - phi
    
    radius = sqrt((R-pellet_R)**2 + (Z-pellet_Z)**2)
    
    atn     = (0.5d0 - 0.5d0*tanh((radius - pellet_radius)/pellet_sig))
    
    atn_psi = (0.5d0 - 0.5d0*tanh(abs(psi- pellet_psi)/pellet_delta_psi))
    
    atn_phi = (0.5d0 - 0.5d0*tanh((phi- pellet_phi)/pellet_length))
    
    particle_source = pellet_amplitude * atn * atn_phi * atn_psi

    return  
  end subroutine pellet_source


end module pellet_module
