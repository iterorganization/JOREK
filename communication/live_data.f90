!> The module contains routines which write certain data toa text file while the code is running.
!!
!! - The file <b>macroscopic_vars.dat</b> is created during the code run and filled with
!!   information about certain run parameters, energy timetraces, growth rates, etc.
!! - The input parameter phys_module::produce_live_data allows to switch the functionality of
!!   this module on (default) or off.
!! - The script extract_live_data.sh in the util/ folder allows to extract certain live data from
!!   the output file. Run it with option -h for usage information.
!! - The script plot_live_data.sh allows to plot live data, e.g., the energy
!!   time traces. Run it with option -h for usage information.
!!
module live_data
  
#include "version.h"
  
  implicit none
  
  private
  public allocate_live_data, init_live_data, write_live_data, write_live_data_vacuum, finalize_live_data
  
  integer,           parameter :: LIVE_DATA_HANDLE = 43 !< File handle for live data file
  character(len=20), parameter :: LIVE_DATA_FILE   = 'macroscopic_vars.dat' !< Live data file
  
  
  
  contains
  
  
  !> Allocate arrays of live data (called in initialise_parameters)
  subroutine allocate_live_data()
    
    use mod_parameters
    use phys_module
    use pellet_module
    
    implicit none
    
    if (allocated(energies)) call tr_deallocate(energies,"energies",CAT_GRID)
    if (nstep .gt. 0) call tr_allocate(energies,1,n_tor,1,2,1,nstep,"energies",CAT_GRID)

#ifdef JECCD
    if (allocated(energies2)) call tr_deallocate(energies2,"energies2",CAT_GRID)
    if (nstep .gt. 0) call tr_allocate(energies2,1,n_tor,1,2,1,nstep,"energies2",CAT_GRID)

    if (allocated(energies3)) call tr_deallocate(energies3,"energies3",CAT_GRID)
    if (nstep .gt. 0) call tr_allocate(energies3,1,n_tor,1,2,1,nstep,"energies3",CAT_GRID)
#endif

    if (allocated(xtime)) call tr_deallocate(xtime,"xtime",CAT_GRID)
    if (nstep .gt. 0) call tr_allocate(xtime,1,nstep,"xtime",CAT_GRID)

    if (allocated(xtime_pellet_R)) call tr_deallocate(xtime_pellet_R,"xtime_pellet_R",CAT_GRID)
    if (nstep .gt. 0)              call tr_allocate(xtime_pellet_R,1,nstep,"xtime_pellet_R")
    if (allocated(xtime_pellet_Z)) call tr_deallocate(xtime_pellet_Z,"xtime_pellet_Z",CAT_GRID)
    if (nstep .gt. 0)              call tr_allocate(xtime_pellet_Z,1,nstep,"xtime_pellet_Z")
    if (allocated(xtime_pellet_psi)) call tr_deallocate(xtime_pellet_psi,"xtime_pellet_psi",CAT_GRID)
    if (nstep .gt. 0)              call tr_allocate(xtime_pellet_psi,1,nstep,"xtime_pellet_psi")
    if (allocated(xtime_pellet_particles)) call tr_deallocate(xtime_pellet_particles,"xtime_pellet_particles",CAT_GRID)
    if (nstep .gt. 0)                      call tr_allocate(xtime_pellet_particles,1,nstep,"xtime_pellet_particles")
    if (allocated(xtime_phys_ablation)) call tr_deallocate(xtime_phys_ablation,"xtime_phys_ablation",CAT_GRID)
    if (nstep .gt. 0)                   call tr_allocate(xtime_phys_ablation,1,nstep,"xtime_phys_ablation")

    if (allocated(R_axis_t)) call tr_deallocate(R_axis_t,"R_axis_t",CAT_UNKNOWN)
    if (nstep .gt. 0) call tr_allocate(R_axis_t,1,index_start+nstep,"R_axis_t",CAT_UNKNOWN)
    
    if (allocated(Z_axis_t)) call tr_deallocate(Z_axis_t,"Z_axis_t",CAT_UNKNOWN)
    if (nstep .gt. 0) call tr_allocate(Z_axis_t,1,index_start+nstep,"Z_axis_t",CAT_UNKNOWN)
    
    if (allocated(psi_axis_t)) call tr_deallocate(psi_axis_t,"psi_axis_t",CAT_UNKNOWN)
    if (nstep .gt. 0) call tr_allocate(psi_axis_t,1,index_start+nstep,"psi_axis_t",CAT_UNKNOWN)
    
    if (allocated(R_xpoint_t)) call tr_deallocate(R_xpoint_t,"R_xpoint_t",CAT_UNKNOWN)
    if (nstep .gt. 0) call tr_allocate(R_xpoint_t,1,index_start+nstep,1,2,"R_xpoint_t",CAT_UNKNOWN)
    
    if (allocated(Z_xpoint_t)) call tr_deallocate(Z_xpoint_t,"Z_xpoint_t",CAT_UNKNOWN)
    if (nstep .gt. 0) call tr_allocate(Z_xpoint_t,1,index_start+nstep,1,2,"Z_xpoint_t",CAT_UNKNOWN)

    if (allocated(psi_xpoint_t)) call tr_deallocate(psi_xpoint_t,"psi_xpoint_t",CAT_UNKNOWN)
    if (nstep .gt. 0) call tr_allocate(psi_xpoint_t,1,index_start+nstep,1,2,"psi_xpoint_t",CAT_UNKNOWN)

    if (allocated(R_bnd_t)) call tr_deallocate(R_bnd_t,"R_bnd_t",CAT_UNKNOWN)
    if (nstep .gt. 0) call tr_allocate(R_bnd_t,1,index_start+nstep,"R_bnd_t",CAT_UNKNOWN)
    
    if (allocated(Z_bnd_t)) call tr_deallocate(Z_bnd_t,"Z_bnd_t",CAT_UNKNOWN)
    if (nstep .gt. 0) call tr_allocate(Z_bnd_t,1,index_start+nstep,"Z_bnd_t",CAT_UNKNOWN)
 
    if (allocated(psi_bnd_t)) call tr_deallocate(psi_bnd_t,"psi_bnd_t",CAT_UNKNOWN)
    if (nstep .gt. 0) call tr_allocate(psi_bnd_t,1,index_start+nstep,"psi_bnd_t",CAT_UNKNOWN)
    
    if (allocated(current_t)) call tr_deallocate(current_t,"current_t",CAT_UNKNOWN)
    if (nstep .gt. 0) call tr_allocate(current_t,1,index_start+nstep,"current_t",CAT_UNKNOWN)
    
    if (allocated(beta_p_t)) call tr_deallocate(beta_p_t,"beta_p_t",CAT_UNKNOWN)
    if (nstep .gt. 0) call tr_allocate(beta_p_t,1,index_start+nstep,"beta_p_t",CAT_UNKNOWN)
    
    if (allocated(beta_t_t)) call tr_deallocate(beta_t_t,"beta_t_t",CAT_UNKNOWN)
    if (nstep .gt. 0) call tr_allocate(beta_t_t,1,index_start+nstep,"beta_t_t",CAT_UNKNOWN)
    
    if (allocated(beta_n_t)) call tr_deallocate(beta_n_t,"beta_n_t",CAT_UNKNOWN)
    if (nstep .gt. 0) call tr_allocate(beta_n_t,1,index_start+nstep,"beta_n_t",CAT_UNKNOWN)
    
    if (allocated(density_in_t)) call tr_deallocate(density_in_t,"density_in_t",CAT_UNKNOWN)
    if (nstep .gt. 0) call tr_allocate(density_in_t,1,index_start+nstep,"density_in_t",CAT_UNKNOWN)
    
    if (allocated(density_out_t)) call tr_deallocate(density_out_t,"density_out_t",CAT_UNKNOWN)
    if (nstep .gt. 0) call tr_allocate(density_out_t,1,index_start+nstep,"density_out_t",CAT_UNKNOWN)
    
    if (allocated(pressure_in_t)) call tr_deallocate(pressure_in_t,"pressure_in_t",CAT_UNKNOWN)
    if (nstep .gt. 0) call tr_allocate(pressure_in_t,1,index_start+nstep,"pressure_in_t",CAT_UNKNOWN)
    
    if (allocated(pressure_out_t)) call tr_deallocate(pressure_out_t,"pressure_out_t",CAT_UNKNOWN)
    if (nstep .gt. 0) call tr_allocate(pressure_out_t,1,index_start+nstep,"pressure_out_t",CAT_UNKNOWN)
    
    if (allocated(heat_src_in_t)) call tr_deallocate(heat_src_in_t,"heat_src_in_t",CAT_UNKNOWN)
    if (nstep .gt. 0) call tr_allocate(heat_src_in_t,1,index_start+nstep,"heat_src_in_t",CAT_UNKNOWN)
    
    if (allocated(heat_src_out_t)) call tr_deallocate(heat_src_out_t,"heat_src_out_t",CAT_UNKNOWN)
    if (nstep .gt. 0) call tr_allocate(heat_src_out_t,1,index_start+nstep,"heat_src_out_t",CAT_UNKNOWN)
    
    if (allocated(part_src_in_t)) call tr_deallocate(part_src_in_t,"part_src_in_t",CAT_UNKNOWN)
    if (nstep .gt. 0) call tr_allocate(part_src_in_t,1,index_start+nstep,"part_src_in_t",CAT_UNKNOWN)
    
    if (allocated(part_src_out_t)) call tr_deallocate(part_src_out_t,"part_src_out_t",CAT_UNKNOWN)
    if (nstep .gt. 0) call tr_allocate(part_src_out_t,1,index_start+nstep,"part_src_out_t",CAT_UNKNOWN)

    if (allocated(E_tot_t)) call tr_deallocate(E_tot_t,"E_tot_t",CAT_UNKNOWN)
    if (nstep .gt. 0) call tr_allocate(E_tot_t,1,index_start+nstep,"E_tot_t",CAT_UNKNOWN)

    if (allocated(helicity_tot_t)) call tr_deallocate(helicity_tot_t,"helicity_tot_t",CAT_UNKNOWN)
    if (nstep .gt. 0) call tr_allocate(helicity_tot_t,1,index_start+nstep,"helicity_tot_t",CAT_UNKNOWN)

    if (allocated(kin_perp_tot_t)) call tr_deallocate(kin_perp_tot_t,"kin_perp_tot_t",CAT_UNKNOWN)
    if (nstep .gt. 0) call tr_allocate(kin_perp_tot_t,1,index_start+nstep,"kin_perp_tot_t",CAT_UNKNOWN)

    if (allocated(kin_par_tot_t)) call tr_deallocate(kin_par_tot_t,"kin_par_tot_t",CAT_UNKNOWN)
    if (nstep .gt. 0) call tr_allocate(kin_par_tot_t,1,index_start+nstep,"kin_par_tot_t",CAT_UNKNOWN)

    if (allocated(thermal_tot_t)) call tr_deallocate(thermal_tot_t,"thermal_tot_t",CAT_UNKNOWN)
    if (nstep .gt. 0) call tr_allocate(thermal_tot_t,1,index_start+nstep,"thermal_tot_t",CAT_UNKNOWN)

    if (allocated(thermal_e_tot_t)) call tr_deallocate(thermal_e_tot_t,"thermal_e_tot_t",CAT_UNKNOWN)
    if (nstep .gt. 0) call tr_allocate(thermal_e_tot_t,1,index_start+nstep,"thermal_e_tot_t",CAT_UNKNOWN)

    if (allocated(thermal_i_tot_t)) call tr_deallocate(thermal_i_tot_t,"thermal_i_tot_t",CAT_UNKNOWN)
    if (nstep .gt. 0) call tr_allocate(thermal_i_tot_t,1,index_start+nstep,"thermal_i_tot_t",CAT_UNKNOWN)

    if (allocated(Wmag_tot_t)) call tr_deallocate(Wmag_tot_t,"Wmag_tot_t",CAT_UNKNOWN)
    if (nstep .gt. 0) call tr_allocate(Wmag_tot_t,1,index_start+nstep,"Wmag_tot_t",CAT_UNKNOWN)

    if (allocated(ohmic_tot_t)) call tr_deallocate(ohmic_tot_t,"ohmic_tot_t",CAT_UNKNOWN)
    if (nstep .gt. 0) call tr_allocate(ohmic_tot_t,1,index_start+nstep,"ohmic_tot_t",CAT_UNKNOWN)

    if (allocated(Magwork_tot_t)) call tr_deallocate(Magwork_tot_t,"Magwork_tot_t",CAT_UNKNOWN)
    if (nstep .gt. 0) call tr_allocate(Magwork_tot_t,1,index_start+nstep,"Magwork_tot_t",CAT_UNKNOWN)

    if (allocated(Ip_tot_t)) call tr_deallocate(Ip_tot_t,"Ip_tot_t",CAT_UNKNOWN)
    if (nstep .gt. 0) call tr_allocate(Ip_tot_t,1,index_start+nstep,"Ip_tot_t",CAT_UNKNOWN)

    if (allocated(flux_Pvn_t)) call tr_deallocate(flux_Pvn_t,"flux_Pvn_t",CAT_UNKNOWN)
    if (nstep .gt. 0) call tr_allocate(flux_Pvn_t,1,index_start+nstep,"flux_Pvn_t",CAT_UNKNOWN)

    if (allocated(flux_qpar_t)) call tr_deallocate(flux_qpar_t,"flux_qpar_t",CAT_UNKNOWN)
    if (nstep .gt. 0) call tr_allocate(flux_qpar_t,1,index_start+nstep,"flux_qpar_t",CAT_UNKNOWN)

    if (allocated(flux_qperp_t)) call tr_deallocate(flux_qperp_t,"flux_qperp_t",CAT_UNKNOWN)
    if (nstep .gt. 0) call tr_allocate(flux_qperp_t,1,index_start+nstep,"flux_qperp_t",CAT_UNKNOWN)

    if (allocated(flux_kinpar_t)) call tr_deallocate(flux_kinpar_t,"flux_kinpar_t",CAT_UNKNOWN)
    if (nstep .gt. 0) call tr_allocate(flux_kinpar_t,1,index_start+nstep,"flux_kinpar_t",CAT_UNKNOWN)

    if (allocated(flux_poynting_t)) call tr_deallocate(flux_poynting_t,"flux_poynting_t",CAT_UNKNOWN)
    if (nstep .gt. 0) call tr_allocate(flux_poynting_t,1,index_start+nstep,"flux_poynting_t",CAT_UNKNOWN)

    if (allocated(dE_tot_dt)) call tr_deallocate(dE_tot_dt,"dE_tot_dt",CAT_UNKNOWN)
    if (nstep .gt. 0) call tr_allocate(dE_tot_dt,1,index_start+nstep,"dE_tot_dt",CAT_UNKNOWN)

    if (allocated(dWmag_tot_dt)) call tr_deallocate(dWmag_tot_dt,"dWmag_tot_dt",CAT_UNKNOWN)
    if (nstep .gt. 0) call tr_allocate(dWmag_tot_dt,1,index_start+nstep,"dWmag_tot_dt",CAT_UNKNOWN)

    if (allocated(dthermal_tot_dt)) call tr_deallocate(dthermal_tot_dt,"dthermal_tot_dt",CAT_UNKNOWN)
    if (nstep .gt. 0) call tr_allocate(dthermal_tot_dt,1,index_start+nstep,"dthermal_tot_dt",CAT_UNKNOWN)

    if (allocated(dkinperp_tot_dt)) call tr_deallocate(dkinperp_tot_dt,"dkinperp_tot_dt",CAT_UNKNOWN)
    if (nstep .gt. 0) call tr_allocate(dkinperp_tot_dt,1,index_start+nstep,"dkinperp_tot_dt",CAT_UNKNOWN)

    if (allocated(dkinpar_tot_dt)) call tr_deallocate(dkinpar_tot_dt,"dkinpar_tot_dt",CAT_UNKNOWN)
    if (nstep .gt. 0) call tr_allocate(dkinpar_tot_dt,1,index_start+nstep,"dkinpar_tot_dt",CAT_UNKNOWN)

    if (allocated(thmwork_tot_t)) call tr_deallocate(thmwork_tot_t,"thmwork_tot_t",CAT_UNKNOWN)
    if (nstep .gt. 0) call tr_allocate(thmwork_tot_t,1,index_start+nstep,"thmwork_tot_t",CAT_UNKNOWN)

    if (allocated(viscopar_dissip_tot_t)) call tr_deallocate(viscopar_dissip_tot_t,"viscopar_dissip_tot_t",CAT_UNKNOWN)
    if (nstep .gt. 0) call tr_allocate(viscopar_dissip_tot_t,1,index_start+nstep,"viscopar_dissip_tot_t",CAT_UNKNOWN)

    if (allocated(friction_dissip_tot_t)) call tr_deallocate(friction_dissip_tot_t,"friction_dissip_tot_t",CAT_UNKNOWN)
    if (nstep .gt. 0) call tr_allocate(friction_dissip_tot_t,1,index_start+nstep,"friction_dissip_tot_t",CAT_UNKNOWN)

    if (allocated(viscopar_flux_t)) call tr_deallocate(viscopar_flux_t,"viscopar_flux_t",CAT_UNKNOWN)
    if (nstep .gt. 0) call tr_allocate(viscopar_flux_t,1,index_start+nstep,"viscopar_flux_t",CAT_UNKNOWN)

    if (allocated(li3_t)) call tr_deallocate(li3_t,"li3_t",CAT_UNKNOWN)
    if (nstep .gt. 0) call tr_allocate(li3_t,1,index_start+nstep,"li3_t",CAT_UNKNOWN)

    if (allocated(li3_tot_t)) call tr_deallocate(li3_tot_t,"li3_tot_t",CAT_UNKNOWN)
    if (nstep .gt. 0) call tr_allocate(li3_tot_t,1,index_start+nstep,"li3_tot_t",CAT_UNKNOWN)

    if (allocated(part_src_tot_t)) call tr_deallocate(part_src_tot_t,"part_src_tot_t",CAT_UNKNOWN)
    if (nstep .gt. 0) call tr_allocate(part_src_tot_t,1,index_start+nstep,"part_src_tot_t",CAT_UNKNOWN)

    if (allocated(heat_src_tot_t)) call tr_deallocate(heat_src_tot_t,"heat_src_tot_t",CAT_UNKNOWN)
    if (nstep .gt. 0) call tr_allocate(heat_src_tot_t,1,index_start+nstep,"heat_src_tot_t",CAT_UNKNOWN)

    if (allocated(volume_t)) call tr_deallocate(volume_t,"volume_t",CAT_UNKNOWN)
    if (nstep .gt. 0) call tr_allocate(volume_t,1,index_start+nstep,"volume_t",CAT_UNKNOWN)

    if (allocated(area_t)) call tr_deallocate(area_t,"area_t",CAT_UNKNOWN)
    if (nstep .gt. 0) call tr_allocate(area_t,1,index_start+nstep,"area_t",CAT_UNKNOWN)

    if (allocated(mag_ener_src_tot)) call tr_deallocate(mag_ener_src_tot,"mag_ener_src_tot",CAT_UNKNOWN)
    if (nstep .gt. 0) call tr_allocate(mag_ener_src_tot,1,index_start+nstep,"mag_ener_src_tot",CAT_UNKNOWN)

    if (allocated(npart_tot_t)) call tr_deallocate(npart_tot_t,"npart_tot_t",CAT_UNKNOWN)
    if (nstep .gt. 0) call tr_allocate(npart_tot_t,1,index_start+nstep,"npart_tot_t",CAT_UNKNOWN)

    if (allocated(dnpart_tot_dt)) call tr_deallocate(dnpart_tot_dt,"dnpart_tot_dt",CAT_UNKNOWN)
    if (nstep .gt. 0) call tr_allocate(dnpart_tot_dt,1,index_start+nstep,"dnpart_tot_dt",CAT_UNKNOWN)

    if (allocated(density_tot_t)) call tr_deallocate(density_tot_t,"density_tot_t",CAT_UNKNOWN)
    if (nstep .gt. 0) call tr_allocate(density_tot_t,1,index_start+nstep,"density_tot_t",CAT_UNKNOWN)

    if (allocated(part_flux_Dpar_t)) call tr_deallocate(part_flux_Dpar_t,"part_flux_Dpar_t",CAT_UNKNOWN)
    if (nstep .gt. 0) call tr_allocate(part_flux_Dpar_t,1,index_start+nstep,"part_flux_Dpar_t",CAT_UNKNOWN)

    if (allocated(part_flux_Dperp_t)) call tr_deallocate(part_flux_Dperp_t,"part_flux_Dperp_t",CAT_UNKNOWN)
    if (nstep .gt. 0) call tr_allocate(part_flux_Dperp_t,1,index_start+nstep,"part_flux_Dperp_t",CAT_UNKNOWN)
    
    if (allocated(part_flux_vpar_t)) call tr_deallocate(part_flux_vpar_t,"part_flux_vpar_t",CAT_UNKNOWN)
    if (nstep .gt. 0) call tr_allocate(part_flux_vpar_t,1,index_start+nstep,"part_flux_vpar_t",CAT_UNKNOWN)

    if (allocated(part_flux_vperp_t)) call tr_deallocate(part_flux_vperp_t,"part_flux_vperp_t",CAT_UNKNOWN)
    if (nstep .gt. 0) call tr_allocate(part_flux_vperp_t,1,index_start+nstep,"part_flux_vperp_t",CAT_UNKNOWN)

    if (allocated(npart_flux_t)) call tr_deallocate(npart_flux_t,"npart_flux_t",CAT_UNKNOWN)
    if (nstep .gt. 0) call tr_allocate(npart_flux_t,1,index_start+nstep,"npart_flux_t",CAT_UNKNOWN)

    if (allocated(dpart_tot_dt)) call tr_deallocate(dpart_tot_dt,"dpart_tot_dt",CAT_UNKNOWN)
    if (nstep .gt. 0) call tr_allocate(dpart_tot_dt,1,index_start+nstep,"dpart_tot_dt",CAT_UNKNOWN)

    return
  end subroutine allocate_live_data
  
  
  
  
  !> Open file, write out headers and some parameters.
  subroutine init_live_data()
    
    use mod_parameters,    only: n_tor, n_plane, n_period, jorek_model, variable_names, n_var
    use phys_module,   only: produce_live_data, mode, mode_type, xpoint, xcase, central_density, sqrt_mu0_rho0, sqrt_mu0_over_rho0, mu_zero
    
    implicit none
    
    logical :: opened
    integer :: n, i
    
    if ( .not. produce_live_data ) return
    
    ! --- Check, that the file handle is not already in use.
    inquire(unit=LIVE_DATA_HANDLE, opened=opened)
    if ( opened ) then
      write(*,*) 'WARNING: LIVE DATA CANNOT BE PRODUCED AS FILE HANDLE IS ALREADY IN USE!'
      produce_live_data = .false.
      return
    end if
    
    open(LIVE_DATA_HANDLE, file=LIVE_DATA_FILE, status='REPLACE', action='WRITE')
    
    ! --- Write some general information
    write(LIVE_DATA_HANDLE,*)  '@rcs_version: ', RCS_VERSION
    write(LIVE_DATA_HANDLE,'(A,I5)') '@jorek_model: ', jorek_model
    write(LIVE_DATA_HANDLE,'(A,I5)') '@n_tor: ', n_tor
    write(LIVE_DATA_HANDLE,'(A,I5)') '@n_plane: ', n_plane
    write(LIVE_DATA_HANDLE,'(A,I5)') '@n_period: ', n_period
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@sqrt_mu0_rho0: ', sqrt_mu0_rho0 
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@sqrt_mu0_over_rho0: ', sqrt_mu0_over_rho0
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@mu_zero: ', mu_zero
    write(LIVE_DATA_HANDLE,'(A)') '@plottable: energies magnetic_energies kinetic_energies growth_rates magnetic_growth_rates  &
                                    kinetic_growth_rates times input_profiles axis current betas particlecontent thermalenergy &
                                    heatingpower particlesource diag_coil_curr pf_coil_curr rmp_coil_curr integrated_energies  &
                                    bnd_fluxes dEdt helicity dissipative_terms work_terms mag_energy_balance                   &
                                    Xpoint_up Xpoint_low bnd_point                                                             &
                                    area volume li3 energy_conservation net_tor_wall_curr dparticles_dt bnd_particle_fluxes    & 
                                     vert_FB_response vert_FB_axis'
    write(LIVE_DATA_HANDLE,'(A,15(A11,1X))') '@variable_names: ', variable_names((/(i, i=1,n_var)/))
    
    ! --- Write file headers indicating what data is in the files.
    write(LIVE_DATA_HANDLE,'(A,I5)') '@n_times: ', 1
    write(LIVE_DATA_HANDLE,'(A)') '@times_xlabel: time step'
    write(LIVE_DATA_HANDLE,'(A)') '@times_xlabel_si: time step'
    write(LIVE_DATA_HANDLE,'(A)') '@times_ylabel: normalized time'
    write(LIVE_DATA_HANDLE,'(A)') '@times_ylabel_si: time [ms]'
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@times_x2si: ', 1.0
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@times_y2si: ', sqrt_mu0_rho0*1.e3
    write(LIVE_DATA_HANDLE,'(A)') '@times_logy: 0'
    write(LIVE_DATA_HANDLE,'(A)') '@times: "step"     "time"'
    
    
    write(LIVE_DATA_HANDLE,'(A,I5)') '@n_energies: ', 2*(n_tor+1)/2
    write(LIVE_DATA_HANDLE,'(A)') '@energies_xlabel: normalized time'
    write(LIVE_DATA_HANDLE,'(A)') '@energies_xlabel_si: time [ms]'
    write(LIVE_DATA_HANDLE,'(A)') '@energies_ylabel: normalized energy'
    write(LIVE_DATA_HANDLE,'(A)') '@energies_ylabel_si: normalized energy'
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@energies_x2si: ', sqrt_mu0_rho0*1.e3
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@energies_y2si: ', 1.0
    write(LIVE_DATA_HANDLE,'(A)') '@energies_logy: 1'
    write(LIVE_DATA_HANDLE,'(A)',advance='no') '@energies: %"time"           '
    do n = 1, n_tor, 2
      write(LIVE_DATA_HANDLE,'(A7,",",I2.2,A2,1x)',advance='no') '"E_{mag', mode(n), '}"'
    end do
    do n = 1, n_tor, 2
      write(LIVE_DATA_HANDLE,'(A7,",",I2.2,A2,1x)',advance='no') '"E_{kin', mode(n), '}"'
    end do
    write(LIVE_DATA_HANDLE,*)


    write(LIVE_DATA_HANDLE,'(A,I5)') '@n_magnetic_energies: ', (n_tor+1)/2
    write(LIVE_DATA_HANDLE,'(A)') '@magnetic_energies_xlabel: normalized time'
    write(LIVE_DATA_HANDLE,'(A)') '@magnetic_energies_xlabel_si: time [ms]'
    write(LIVE_DATA_HANDLE,'(A)') '@magnetic_energies_ylabel: normalized magnetic energy'
    write(LIVE_DATA_HANDLE,'(A)') '@magnetic_energies_ylabel_si: normalized magnetic energy'
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@magnetic_energies_x2si: ', sqrt_mu0_rho0*1.e3
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@magnetic_energies_y2si: ', 1.0
    write(LIVE_DATA_HANDLE,'(A)') '@magnetic_energies_logy: 1'
    write(LIVE_DATA_HANDLE,'(A)',advance='no') '@magnetic_energies: %"time"           '
    do n = 1, n_tor, 2
      write(LIVE_DATA_HANDLE,'(A7,",",I2.2,A2,1x)',advance='no') '"E_{mag', mode(n), '}"'
    end do
    write(LIVE_DATA_HANDLE,*)


    write(LIVE_DATA_HANDLE,'(A,I5)') '@n_kinetic_energies: ', (n_tor+1)/2
    write(LIVE_DATA_HANDLE,'(A)') '@kinetic_energies_xlabel: normalized time'
    write(LIVE_DATA_HANDLE,'(A)') '@kinetic_energies_xlabel_si: time [ms]'
    write(LIVE_DATA_HANDLE,'(A)') '@kinetic_energies_ylabel: normalized kinetic energy'
    write(LIVE_DATA_HANDLE,'(A)') '@kinetic_energies_ylabel_si: normalized kinetic energy'
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@kinetic_energies_x2si: ', sqrt_mu0_rho0*1.e3
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@kinetic_energies_y2si: ', 1.0
    write(LIVE_DATA_HANDLE,'(A)') '@kinetic_energies_logy: 1'
    write(LIVE_DATA_HANDLE,'(A)',advance='no') '@kinetic_energies: %"time"           '
    do n = 1, n_tor, 2
      write(LIVE_DATA_HANDLE,'(A7,",",I2.2,A2,1x)',advance='no') '"E_{kin', mode(n), '}"'
    end do
    write(LIVE_DATA_HANDLE,*)
    
    write(LIVE_DATA_HANDLE,'(A,I5)') '@n_growth_rates: ', 2*(n_tor+1)/2
    write(LIVE_DATA_HANDLE,'(A)') '@growth_rates_xlabel: normalized time'
    write(LIVE_DATA_HANDLE,'(A)') '@growth_rates_xlabel_si: time [ms]'
    write(LIVE_DATA_HANDLE,'(A)') '@growth_rates_ylabel: normalized growth rate'
    write(LIVE_DATA_HANDLE,'(A)') '@growth_rates_ylabel_si: growth rate [1/s]'
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@growth_rates_x2si: ', sqrt_mu0_rho0*1.e3
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@growth_rates_y2si: ', 1./sqrt_mu0_rho0
    write(LIVE_DATA_HANDLE,'(A)') '@growth_rates_logy: 1'
    write(LIVE_DATA_HANDLE,'(A)',advance='no') '@growth_rates: %"time"           '
    do n = 1, n_tor, 2
      write(LIVE_DATA_HANDLE,'(A7,",",I2.2,A2,1x)',advance='no') '"G_{mag', mode(n), '}"'
    end do
    do n = 1, n_tor, 2
      write(LIVE_DATA_HANDLE,'(A7,",",I2.2,A2,1x)',advance='no') '"G_{kin', mode(n), '}"'
    end do
    write(LIVE_DATA_HANDLE,*)

    write(LIVE_DATA_HANDLE,'(A,I5)') '@n_magnetic_growth_rates: ', (n_tor+1)/2
    write(LIVE_DATA_HANDLE,'(A)') '@magnetic_growth_rates_xlabel: normalized time'
    write(LIVE_DATA_HANDLE,'(A)') '@magnetic_growth_rates_xlabel_si: time [ms]'
    write(LIVE_DATA_HANDLE,'(A)') '@magnetic_growth_rates_ylabel: normalized growth rate'
    write(LIVE_DATA_HANDLE,'(A)') '@magnetic_growth_rates_ylabel_si: growth rate [1/s]'
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@magnetic_growth_rates_x2si: ', sqrt_mu0_rho0*1.e3
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@magnetic_growth_rates_y2si: ', 1./sqrt_mu0_rho0
    write(LIVE_DATA_HANDLE,'(A)') '@magnetic_growth_rates_logy: 0'
    write(LIVE_DATA_HANDLE,'(A)',advance='no') '@magnetic_growth_rates: %"time"           '
    do n = 1, n_tor, 2
      write(LIVE_DATA_HANDLE,'(A7,",",I2.2,A2,1x)',advance='no') '"G_{mag', mode(n), '}"'
    end do
    write(LIVE_DATA_HANDLE,*)
 
    write(LIVE_DATA_HANDLE,'(A,I5)') '@n_kinetic_growth_rates: ', (n_tor+1)/2
    write(LIVE_DATA_HANDLE,'(A)') '@kinetic_growth_rates_xlabel: normalized time'
    write(LIVE_DATA_HANDLE,'(A)') '@kinetic_growth_rates_xlabel_si: time [ms]'
    write(LIVE_DATA_HANDLE,'(A)') '@kinetic_growth_rates_ylabel: normalized growth rate'
    write(LIVE_DATA_HANDLE,'(A)') '@kinetic_growth_rates_ylabel_si: growth rate [1/s]'
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@kinetic_growth_rates_x2si: ', sqrt_mu0_rho0*1.e3
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@kinetic_growth_rates_y2si: ', 1./sqrt_mu0_rho0
    write(LIVE_DATA_HANDLE,'(A)') '@kinetic_growth_rates_logy: 1'
    write(LIVE_DATA_HANDLE,'(A)',advance='no') '@kinetic_growth_rates: %"time"           '
    do n = 1, n_tor, 2
      write(LIVE_DATA_HANDLE,'(A7,",",I2.2,A2,1x)',advance='no') '"G_{kin', mode(n), '}"'
    end do
    write(LIVE_DATA_HANDLE,*)

    write(LIVE_DATA_HANDLE,'(A,I5)') '@n_axis: ', 3
    write(LIVE_DATA_HANDLE,'(A)') '@axis_xlabel: normalized time'
    write(LIVE_DATA_HANDLE,'(A)') '@axis_xlabel_si: time [ms]'
    write(LIVE_DATA_HANDLE,'(A)') '@axis_ylabel: Magnetic axis properties'
    write(LIVE_DATA_HANDLE,'(A)') '@axis_ylabel_si: Magnetic axis properties'
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@axis_x2si: ', sqrt_mu0_rho0*1.e3
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@axis_y2si: ', 1.0
    write(LIVE_DATA_HANDLE,'(A)') '@axis_logy: 0'
    write(LIVE_DATA_HANDLE,'(A)') '@axis: %"time"           "R position"              "Z position"           "Psi on axis"'
    write(LIVE_DATA_HANDLE,*)

    write(LIVE_DATA_HANDLE,'(A,I5)') '@n_Xpoint_low: ', 3
    write(LIVE_DATA_HANDLE,'(A)') '@Xpoint_low_xlabel: normalized time'
    write(LIVE_DATA_HANDLE,'(A)') '@Xpoint_low_xlabel_si: time [ms]'
    write(LIVE_DATA_HANDLE,'(A)') '@Xpoint_low_ylabel: Lower X-point properties'
    write(LIVE_DATA_HANDLE,'(A)') '@Xpoint_low_ylabel_si: Lower X-point properties'
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@Xpoint_low_x2si: ', sqrt_mu0_rho0*1.e3
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@Xpoint_low_y2si: ', 1.0
    write(LIVE_DATA_HANDLE,'(A)') '@Xpoint_low_logy: 0'
    write(LIVE_DATA_HANDLE,'(A)') '@Xpoint_low: %"time"           "R position"              "Z position"           "Psi on lower X-point"'
    write(LIVE_DATA_HANDLE,*)

    write(LIVE_DATA_HANDLE,'(A,I5)') '@n_Xpoint_up: ', 3
    write(LIVE_DATA_HANDLE,'(A)') '@Xpoint_up_xlabel: normalized time'
    write(LIVE_DATA_HANDLE,'(A)') '@Xpoint_up_xlabel_si: time [ms]'
    write(LIVE_DATA_HANDLE,'(A)') '@Xpoint_up_ylabel: Upper X-point properties'
    write(LIVE_DATA_HANDLE,'(A)') '@Xpoint_up_ylabel_si: Upper X-point properties'
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@Xpoint_up_x2si: ', sqrt_mu0_rho0*1.e3
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@Xpoint_up_y2si: ', 1.0
    write(LIVE_DATA_HANDLE,'(A)') '@Xpoint_up_logy: 0'
    write(LIVE_DATA_HANDLE,'(A)') '@Xpoint_up: %"time"           "R position"              "Z position"           "Psi on upper X-point"'
    write(LIVE_DATA_HANDLE,*)

    write(LIVE_DATA_HANDLE,'(A,I5)') '@n_bnd_point: ', 3
    write(LIVE_DATA_HANDLE,'(A)') '@bnd_point_xlabel: normalized time'
    write(LIVE_DATA_HANDLE,'(A)') '@bnd_point_xlabel_si: time [ms]'
    write(LIVE_DATA_HANDLE,'(A)') '@bnd_point_ylabel: Boundary point (defining LCFS) properties'
    write(LIVE_DATA_HANDLE,'(A)') '@bnd_point_ylabel_si: Boundary point (defining LCFS) properties'
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@bnd_point_x2si: ', sqrt_mu0_rho0*1.e3
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@bnd_point_y2si: ', 1.0
    write(LIVE_DATA_HANDLE,'(A)') '@bnd_point_logy: 0'
    write(LIVE_DATA_HANDLE,'(A)') '@bnd_point: %"time"           "R position"              "Z position"           "Psi on boundary point"'
    write(LIVE_DATA_HANDLE,*)

    write(LIVE_DATA_HANDLE,'(A,I5)') '@n_integrated_energies: ', 5 
    write(LIVE_DATA_HANDLE,'(A)') '@integrated_energies_xlabel: normalized time'
    write(LIVE_DATA_HANDLE,'(A)') '@integrated_energies__xlabel_si: time [ms]'
    write(LIVE_DATA_HANDLE,'(A)') '@integrated_energies_ylabel: Total integrated energies [J]'
    write(LIVE_DATA_HANDLE,'(A)') '@integrated_energies_ylabel_si: Total integrated energies [J]'
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@integrated_energies_x2si: ', sqrt_mu0_rho0*1.e3
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@integrated_energies_y2si: ', 1.0
    write(LIVE_DATA_HANDLE,'(A)') '@integrated_energies_logy: 0'
    write(LIVE_DATA_HANDLE,'(A)') '@integrated_energies: %"time"           "Total energy"              "Magnetic"           "Kinetic parallel"    &
                                   "Kinetic perpendicular"                 "Thermal energy"     '
    write(LIVE_DATA_HANDLE,*)

    write(LIVE_DATA_HANDLE,'(A,I5)') '@n_bnd_fluxes: ', 4 
    write(LIVE_DATA_HANDLE,'(A)') '@bnd_fluxes_xlabel: normalized time'
    write(LIVE_DATA_HANDLE,'(A)') '@bnd_fluxes_xlabel_si: time [ms]'
    write(LIVE_DATA_HANDLE,'(A)') '@bnd_fluxes_ylabel: Total boundary fluxes [W]'
    write(LIVE_DATA_HANDLE,'(A)') '@bnd_fluxes_ylabel_si: Total boundary fluxes [W]'
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@bnd_fluxes_x2si: ', sqrt_mu0_rho0*1.e3
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@bnd_fluxes_y2si: ', 1.0
    write(LIVE_DATA_HANDLE,'(A)') '@bnd_fluxes_logy: 0'
    write(LIVE_DATA_HANDLE,'(A)') '@bnd_fluxes: %"time"       "p vn"  "kinpar-flux"    "qn-par"    "qn-perp"   '
    write(LIVE_DATA_HANDLE,*)

    write(LIVE_DATA_HANDLE,'(A,I5)') '@n_bnd_particle_fluxes: ', 5 
    write(LIVE_DATA_HANDLE,'(A)') '@bnd_particle_fluxes_xlabel: normalized time'
    write(LIVE_DATA_HANDLE,'(A)') '@bnd_particle_fluxes_xlabel_si: time [ms]'
    write(LIVE_DATA_HANDLE,'(A)') '@bnd_particle_fluxes_ylabel: Total particle boundary fluxes [particles/s]'
    write(LIVE_DATA_HANDLE,'(A)') '@bnd_particle_fluxes_ylabel_si: Total particle boundary fluxes [particles/s]'
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@bnd_particle_fluxes_x2si: ', sqrt_mu0_rho0*1.e3
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@bnd_particle_fluxes_y2si: ', 1.0
    write(LIVE_DATA_HANDLE,'(A)') '@bnd_particle_fluxes_logy: 0'
    write(LIVE_DATA_HANDLE,'(A)') '@bnd_particle_fluxes: %"time"    "Dpar-flux"  "Dperp-flux"    "Vpar-flux"    "Vperp-flux"  "neutral-flux"   '
    write(LIVE_DATA_HANDLE,*)
    
    write(LIVE_DATA_HANDLE,'(A,I5)') '@n_current: ', 3 
    write(LIVE_DATA_HANDLE,'(A)') '@current_xlabel: normalized time'
    write(LIVE_DATA_HANDLE,'(A)') '@current_xlabel_si: time [ms]'
    write(LIVE_DATA_HANDLE,'(A)') '@current_ylabel: plasma current [A]'
    write(LIVE_DATA_HANDLE,'(A)') '@current_ylabel_si: plasma current [MA]'
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@current_x2si: ', sqrt_mu0_rho0*1.e3
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@current_y2si: ', 1.e-6
    write(LIVE_DATA_HANDLE,'(A)') '@current_logy: 0'
    write(LIVE_DATA_HANDLE,'(A)') '@current: %"time"       "Total"    "Inside LCFS"   "Outside LCFS" '
    write(LIVE_DATA_HANDLE,*)

    write(LIVE_DATA_HANDLE,'(A,I5)') '@n_helicity: ', 1
    write(LIVE_DATA_HANDLE,'(A)') '@helicity_xlabel: normalized time'
    write(LIVE_DATA_HANDLE,'(A)') '@helicity_xlabel_si: time [ms]'
    write(LIVE_DATA_HANDLE,'(A)') '@helicity_ylabel: Total helicity [Wb^2]'
    write(LIVE_DATA_HANDLE,'(A)') '@helicity_ylabel_si: Total helicity [Wb^2]'
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@helicity_x2si: ', sqrt_mu0_rho0*1.e3
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@helicity_y2si: ', 1.0
    write(LIVE_DATA_HANDLE,'(A)') '@helicity_logy: 0'
    write(LIVE_DATA_HANDLE,'(A)') '@helicity: %"time"           "Helicity"'
    write(LIVE_DATA_HANDLE,*)

    write(LIVE_DATA_HANDLE,'(A,I5)') '@n_area: ', 1
    write(LIVE_DATA_HANDLE,'(A)') '@area_xlabel: normalized time'
    write(LIVE_DATA_HANDLE,'(A)') '@area_xlabel_si: time [ms]'
    write(LIVE_DATA_HANDLE,'(A)') '@area_ylabel: Total area [m^2]'
    write(LIVE_DATA_HANDLE,'(A)') '@area_ylabel_si: Total area [m^2]'
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@area_x2si: ', sqrt_mu0_rho0*1.e3
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@area_y2si: ', 1.0
    write(LIVE_DATA_HANDLE,'(A)') '@area_logy: 0'
    write(LIVE_DATA_HANDLE,'(A)') '@area: %"time"           "Plasma cross sectional area"'
    write(LIVE_DATA_HANDLE,*)

    write(LIVE_DATA_HANDLE,'(A,I5)') '@n_mag_energy_src: ', 1
    write(LIVE_DATA_HANDLE,'(A)') '@mag_energy_src_xlabel: normalized time'
    write(LIVE_DATA_HANDLE,'(A)') '@mag_energy_src_ylabel: Total magnetic energy source [W]'
    write(LIVE_DATA_HANDLE,'(A)') '@mag_energy_src_logy: 0'
    write(LIVE_DATA_HANDLE,'(A)') '@mag_energy_src: %"time"       "Mag source"'
    write(LIVE_DATA_HANDLE,*)

    write(LIVE_DATA_HANDLE,'(A,I5)') '@n_volume: ', 1
    write(LIVE_DATA_HANDLE,'(A)') '@volume_xlabel: normalized time'
    write(LIVE_DATA_HANDLE,'(A)') '@volume_xlabel_si: time [ms]'
    write(LIVE_DATA_HANDLE,'(A)') '@volume_ylabel: Total volume [m^3]'
    write(LIVE_DATA_HANDLE,'(A)') '@volume_ylabel_si: Total volume [m^3]'
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@volume_x2si: ', sqrt_mu0_rho0*1.e3
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@volume_y2si: ', 1.0
    write(LIVE_DATA_HANDLE,'(A)') '@volume_logy: 0'
    write(LIVE_DATA_HANDLE,'(A)') '@volume: %"time"           "Plasma volume"'
    write(LIVE_DATA_HANDLE,*)

    write(LIVE_DATA_HANDLE,'(A,I5)') '@n_li3: ', 2
    write(LIVE_DATA_HANDLE,'(A)') '@li3_xlabel: normalized time'
    write(LIVE_DATA_HANDLE,'(A)') '@li3_xlabel_si: time [ms]'
    write(LIVE_DATA_HANDLE,'(A)') '@li3_ylabel: li(3) '
    write(LIVE_DATA_HANDLE,'(A)') '@li3_ylabel_si: li(3) '
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@li3_x2si: ', sqrt_mu0_rho0*1.e3
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@li3_y2si: ', 1.0
    write(LIVE_DATA_HANDLE,'(A)') '@li3_logy: 0'
    write(LIVE_DATA_HANDLE,'(A)') '@li3: %"time"         "inside separatrix"  "All domain" '
    write(LIVE_DATA_HANDLE,*)

    write(LIVE_DATA_HANDLE,'(A,I5)') '@n_energy_conservation: ', 2
    write(LIVE_DATA_HANDLE,'(A)') '@energy_conservation_xlabel: normalized time'
    write(LIVE_DATA_HANDLE,'(A)') '@energy_conservation_xlabel_si: time [ms]'
    write(LIVE_DATA_HANDLE,'(A)') '@energy_conservation_ylabel: Total energy conservation'
    write(LIVE_DATA_HANDLE,'(A)') '@energy_conservation_ylabel_si: Total energy conservation'
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@energy_conservation_x2si: ', sqrt_mu0_rho0*1.e3
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@energy_conservation_y2si: ', 1.0
    write(LIVE_DATA_HANDLE,'(A)') '@energy_conservation_logy: 0'
    write(LIVE_DATA_HANDLE,'(A)') '@energy_conservation: %"time"       "-dEtotdt"     "Sum bnd fluxes + sources + dissipative terms" '
    write(LIVE_DATA_HANDLE,*)

    write(LIVE_DATA_HANDLE,'(A,I5)') '@n_mag_energy_balance: ', 6
    write(LIVE_DATA_HANDLE,'(A)') '@mag_energy_balance_xlabel: normalized time'
    write(LIVE_DATA_HANDLE,'(A)') '@mag_energy_balance_xlabel_si: time [ms]'
    write(LIVE_DATA_HANDLE,'(A)') '@mag_energy_balance_ylabel: Magnetic energy balance (W)'
    write(LIVE_DATA_HANDLE,'(A)') '@mag_energy_balance_ylabel_si: Magnetic energy balance (W)'
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@mag_energy_balance_x2si: ', sqrt_mu0_rho0*1.e3
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@mag_energy_balance_y2si: ', 1.0
    write(LIVE_DATA_HANDLE,'(A)') '@mag_energy_balance_logy: 0'
    write(LIVE_DATA_HANDLE,'(A)') '@mag_energy_balance: %"time"    "dWmagdt"   "Ohmic"  "Poynting"  "JxB.v"  "magSource"  "sum all losses + sources"  '
    write(LIVE_DATA_HANDLE,*)
 
#if (defined WITH_Neutrals) || (defined WITH_Impurities)
    write(LIVE_DATA_HANDLE,'(A,I5)') '@n_dissipative_terms: ', 5
#else
    write(LIVE_DATA_HANDLE,'(A,I5)') '@n_dissipative_terms: ', 3
#endif
    write(LIVE_DATA_HANDLE,'(A)') '@dissipative_terms_xlabel: normalized time'
    write(LIVE_DATA_HANDLE,'(A)') '@dissipative_terms_xlabel_si: time [ms]'
    write(LIVE_DATA_HANDLE,'(A)') '@dissipative_terms_ylabel: Dissipative powers [W]'
    write(LIVE_DATA_HANDLE,'(A)') '@dissipative_terms_ylabel_si: Dissipative powers [W]'
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@dissipative_terms_x2si: ', sqrt_mu0_rho0*1.e3
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@dissipative_terms_y2si: ', 1.0
    write(LIVE_DATA_HANDLE,'(A)') '@dissipative_terms_logy: 0'
#if (defined WITH_Neutrals) || (defined WITH_Impurities)
    write(LIVE_DATA_HANDLE,'(A)') '@dissipative_terms: %"time"         "Ohmic power"   "Frictional heating"   "Parallel viscosity power"  &
                                                        "Radiated power"  "Ionization power"'
#else
    write(LIVE_DATA_HANDLE,'(A)') '@dissipative_terms: %"time"         "Ohmic power"   "Frictional heating"   "Parallel viscosity power" '
#endif
    write(LIVE_DATA_HANDLE,*)

    write(LIVE_DATA_HANDLE,'(A,I5)') '@n_work_terms: ', 2 
    write(LIVE_DATA_HANDLE,'(A)') '@work_terms_xlabel: normalized time'
    write(LIVE_DATA_HANDLE,'(A)') '@work_terms_xlabel_si: time [ms]'
    write(LIVE_DATA_HANDLE,'(A)') '@work_terms_ylabel: Total work [W]'
    write(LIVE_DATA_HANDLE,'(A)') '@work_terms_ylabel_si: Total work [W]'
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@work_terms_x2si: ', sqrt_mu0_rho0*1.e3
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@work_terms_y2si: ', 1.0
    write(LIVE_DATA_HANDLE,'(A)') '@work_terms_logy: 0'
    write(LIVE_DATA_HANDLE,'(A)') '@work_terms: %"time"     "Magnetic = JxB~nabla p"   "Thermal = vpar*nabla p"  '
    write(LIVE_DATA_HANDLE,*)

    write(LIVE_DATA_HANDLE,'(A,I5)') '@n_dEdt: ', 5 
    write(LIVE_DATA_HANDLE,'(A)') '@dEdt_xlabel: normalized time'
    write(LIVE_DATA_HANDLE,'(A)') '@dEdt_xlabel_si: time [ms]'
    write(LIVE_DATA_HANDLE,'(A)') '@dEdt_ylabel: -dEnergydt [W]'
    write(LIVE_DATA_HANDLE,'(A)') '@dEdt_ylabel_si: -dEnergydt [W]'
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@dEdt_x2si: ', sqrt_mu0_rho0*1.e3
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@dEdt_y2si: ', 1.0
    write(LIVE_DATA_HANDLE,'(A)') '@dEdt_logy: 0'
    write(LIVE_DATA_HANDLE,'(A)') '@dEdt: %"time"    "Etot"  "Wmagtot"  "thermaltot"   "kinperptot"  "kinpartot"      '
    write(LIVE_DATA_HANDLE,*)

    write(LIVE_DATA_HANDLE,'(A,I5)') '@n_dparticles_dt: ', 3 
    write(LIVE_DATA_HANDLE,'(A)') '@dparticles_dt_xlabel: normalized time'
    write(LIVE_DATA_HANDLE,'(A)') '@dparticles_dt_xlabel_si: time [ms]'
    write(LIVE_DATA_HANDLE,'(A)') '@dparticles_dt_ylabel: dparticlesdt [1/s]'
    write(LIVE_DATA_HANDLE,'(A)') '@dparticles_dt_ylabel_si: dparticlesdt [1/s]'
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@dparticles_dt_x2si: ', sqrt_mu0_rho0*1.e3
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@dparticles_dt_y2si: ', 1.0
    write(LIVE_DATA_HANDLE,'(A)') '@dparticles_dt_logy: 0'
    write(LIVE_DATA_HANDLE,'(A)') '@dparticles_dt: %"time"  "Total"  "Ions" "Neutrals"     '
    write(LIVE_DATA_HANDLE,*)
    
    write(LIVE_DATA_HANDLE,'(A,I5)') '@n_betas: ', 3
    write(LIVE_DATA_HANDLE,'(A)') '@betas_xlabel: normalized time'
    write(LIVE_DATA_HANDLE,'(A)') '@betas_xlabel_si: time [ms]'
    write(LIVE_DATA_HANDLE,'(A)') '@betas_ylabel: plasma beta [%]'
    write(LIVE_DATA_HANDLE,'(A)') '@betas_ylabel_si: plasma beta [%]'
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@betas_x2si: ', sqrt_mu0_rho0*1.e3
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@betas_y2si: ', 1.0
    write(LIVE_DATA_HANDLE,'(A)') '@betas_logy: 0'
    write(LIVE_DATA_HANDLE,'(A)') '@betas: %"time"           "beta poloidal"       "beta toroidal"       "beta normalized"'
    write(LIVE_DATA_HANDLE,*)
    
    write(LIVE_DATA_HANDLE,'(A,I5)') '@n_particlecontent: ', 4 
    write(LIVE_DATA_HANDLE,'(A)') '@particlecontent_xlabel: normalized time'
    write(LIVE_DATA_HANDLE,'(A)') '@particlecontent_xlabel_si: time [ms]'
    write(LIVE_DATA_HANDLE,'(A)') '@particlecontent_ylabel: particle content'
    write(LIVE_DATA_HANDLE,'(A)') '@particlecontent_ylabel_si: particle content'
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@particlecontent_x2si: ', sqrt_mu0_rho0*1.e3
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@particlecontent_y2si: ', 1.0
    write(LIVE_DATA_HANDLE,'(A)') '@particlecontent_logy: 0'
#ifdef WITH_Impurities
    write(LIVE_DATA_HANDLE,'(A)') '@particlecontent: %"time"  "Total main ions" "Main ions inside LCFS"  "Main ions outside LCFS" "Total impurities"'
#else    
    write(LIVE_DATA_HANDLE,'(A)') '@particlecontent: %"time"  "Total ions" "Ions inside LCFS"  "Ions outside LCFS" "Total neutrals"'
#endif
    write(LIVE_DATA_HANDLE,*)
    
    write(LIVE_DATA_HANDLE,'(A,I5)') '@n_thermalenergy: ', 2 
    write(LIVE_DATA_HANDLE,'(A)') '@thermalenergy_xlabel: normalized time'
    write(LIVE_DATA_HANDLE,'(A)') '@thermalenergy_xlabel_si: time [ms]'
    write(LIVE_DATA_HANDLE,'(A)') '@thermalenergy_ylabel: Thermal energy [J]'
    write(LIVE_DATA_HANDLE,'(A)') '@thermalenergy_ylabel_si: thermal energy [MJ]'
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@thermalenergy_x2si: ', sqrt_mu0_rho0*1.e3
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@thermalenergy_y2si: ', 1.e-6
    write(LIVE_DATA_HANDLE,'(A)') '@thermalenergy_logy: 0'
    write(LIVE_DATA_HANDLE,'(A)') '@thermalenergy: %"time"   "inside separatrix"   "outside separatrix"'
    write(LIVE_DATA_HANDLE,*)
    
    write(LIVE_DATA_HANDLE,'(A,I5)') '@n_heatingpower: ', 3 
    write(LIVE_DATA_HANDLE,'(A)') '@heatingpower_xlabel: normalized time'
    write(LIVE_DATA_HANDLE,'(A)') '@heatingpower_xlabel_si: time [ms]'
    write(LIVE_DATA_HANDLE,'(A)') '@heatingpower_ylabel: heating power [W]'
    write(LIVE_DATA_HANDLE,'(A)') '@heatingpower_ylabel_si: heating power [MW]'
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@heatingpower_x2si: ', sqrt_mu0_rho0*1.e3
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@heatingpower_y2si: ', 1.e-6
    write(LIVE_DATA_HANDLE,'(A)') '@heatingpower_logy: 0'
    write(LIVE_DATA_HANDLE,'(A)') '@heatingpower: %"time"  "Total"     "inside separatrix"   "outside separatrix"'
    write(LIVE_DATA_HANDLE,*)
    
    write(LIVE_DATA_HANDLE,'(A,I5)') '@n_particlesource: ', 3 
    write(LIVE_DATA_HANDLE,'(A)') '@particlesource_xlabel: normalized time'
    write(LIVE_DATA_HANDLE,'(A)') '@particlesource_xlabel_si: time [ms]'
    write(LIVE_DATA_HANDLE,'(A)') '@particlesource_ylabel: particle source [10^20/m^3/s]'
    write(LIVE_DATA_HANDLE,'(A)') '@particlesource_ylabel_si: particle source [10^20/m^3/s]'
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@particlesource_x2si: ', sqrt_mu0_rho0*1.e3
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@particlesource_y2si: ', 1.0
    write(LIVE_DATA_HANDLE,'(A)') '@particlesource_logy: 0'
    write(LIVE_DATA_HANDLE,'(A)') '@particlesource: %"time"  "Total"   "inside separatrix"   "outside separatrix"   '
    write(LIVE_DATA_HANDLE,*)
    
    ! --- Call the model-specific part of the init_live_data routine
    call init_live_data_model(LIVE_DATA_HANDLE) 
    
    close(LIVE_DATA_HANDLE)
    
  end subroutine init_live_data
  
  
  
  !> Write out data to text files during the code run.
  subroutine write_live_data(index)
    
    use mod_parameters,  only: n_tor
    use phys_module, only: xtime, energies, produce_live_data, R_axis_t, Z_axis_t, Psi_axis_t,     &
      R_xpoint_t, Z_xpoint_t, psi_xpoint_t, R_bnd_t, Z_bnd_t, Psi_bnd_t,     &
      current_t, beta_p_t, beta_t_t, beta_n_t, density_in_t, density_out_t, pressure_in_t,               &
      pressure_out_t, heat_src_in_t, heat_src_out_t, part_src_in_t, part_src_out_t, &
      E_tot_t, Helicity_tot_t, Kin_perp_tot_t, thermal_tot_t, kin_par_tot_t, ohmic_tot_t,      &
      Wmag_tot_t, Ip_tot_t, flux_pvn_t, flux_qpar_t, flux_qperp_t, flux_kinpar_t, dE_tot_dt, &
      dWmag_tot_dt, dthermal_tot_dt, dkinpar_tot_dt, dkinperp_tot_dt,  Magwork_tot_t,   &
      thmwork_tot_t, viscopar_dissip_tot_t, viscopar_flux_t, li3_t, friction_dissip_tot_t,     &
      li3_tot_t, part_src_tot_t, heat_src_tot_t, volume_t, area_t, mag_ener_src_tot, eta_ohmic, eta, &
      dpart_tot_dt, part_flux_Dpar_t, part_flux_Dperp_t, part_flux_vpar_t, part_flux_vperp_t, &
      dnpart_tot_dt, npart_tot_t, npart_flux_t, density_tot_t, flux_poynting_t, xtime_rad_power, &
      xtime_E_ion_power, thermal_e_tot_t, thermal_i_tot_t


    implicit none
    
    integer, intent(in) :: index !< Timestep index to write data for
    
    integer :: i, j
    real*8  :: e1, e2, growth_rate, sum_fluxes_dissip, sum_mag_energy_terms
    
    if ( .not. produce_live_data ) return
    
    open(LIVE_DATA_HANDLE, file=LIVE_DATA_FILE, status='OLD', position='APPEND', action='WRITE')
    
    ! --- Write data to the files.
    write(LIVE_DATA_HANDLE,'(A,I6,1X,ES17.9)') '@times:', index, xtime(index)
    write(LIVE_DATA_HANDLE,'(A,ES17.9)',advance='no') '@energies:', xtime(index)
    do j = 1, 2
      do i = 1, n_tor, 2
        write(LIVE_DATA_HANDLE,'(ES17.9)',advance='no') sum(energies(max(i-1,1):i,j,index))
      end do
    end do
    write(LIVE_DATA_HANDLE,*)
    write(LIVE_DATA_HANDLE,'(A,ES17.9)',advance='no') '@magnetic_energies:', xtime(index)
    do i = 1, n_tor, 2
      write(LIVE_DATA_HANDLE,'(ES17.9)',advance='no') sum(energies(max(i-1,1):i,1,index))
    end do
    write(LIVE_DATA_HANDLE,*)
    write(LIVE_DATA_HANDLE,'(A,ES17.9)',advance='no') '@kinetic_energies:', xtime(index)
    do i = 1, n_tor, 2
      write(LIVE_DATA_HANDLE,'(ES17.9)',advance='no') sum(energies(max(i-1,1):i,2,index))
    end do
    write(LIVE_DATA_HANDLE,*)
    if ( index > 1 ) then
      write(LIVE_DATA_HANDLE,'(A,ES17.9)',advance='no') '@growth_rates:', &
        (xtime(index)+xtime(index-1))/2.d0
      do j = 1, 2
        do i = 1, n_tor, 2
          e1 = sum(energies(max(i-1,1):i,j,index))
          e2 = sum(energies(max(i-1,1):i,j,index-1))
          if ( (e1 .GT. 0.) .and. (e2 .GT. 0.) ) then
             growth_rate = 0.5d0 * ( log(e1) - log(e2) ) / (xtime(index)-xtime(index-1))
          else
             growth_rate = 0.d0
          endif
          write(LIVE_DATA_HANDLE,'(ES17.9)',advance='no') growth_rate
        end do
      end do 
      write(LIVE_DATA_HANDLE,*)
      write(LIVE_DATA_HANDLE,'(A,ES17.9)',advance='no') '@magnetic_growth_rates:', &
        (xtime(index)+xtime(index-1))/2.d0
      do i = 1, n_tor, 2
        e1 = sum(energies(max(i-1,1):i,1,index))
        e2 = sum(energies(max(i-1,1):i,1,index-1))
        if ( (e1 .GT. 0.) .and. (e2 .GT. 0.) ) then
           growth_rate = 0.5d0 * ( log(e1) - log(e2) ) / (xtime(index)-xtime(index-1))
        else
           growth_rate = 0.d0
        endif
        write(LIVE_DATA_HANDLE,'(ES17.9)',advance='no') growth_rate
      end do
      write(LIVE_DATA_HANDLE,*)
      write(LIVE_DATA_HANDLE,'(A,ES17.9)',advance='no') '@kinetic_growth_rates:', &
        (xtime(index)+xtime(index-1))/2.d0
      do i = 1, n_tor, 2
        e1 = sum(energies(max(i-1,1):i,2,index))
        e2 = sum(energies(max(i-1,1):i,2,index-1))
        if ( (e1 .GT. 0.) .and. (e2 .GT. 0.) ) then
           growth_rate = 0.5d0 * ( log(e1) - log(e2) ) / (xtime(index)-xtime(index-1))
        else
           growth_rate = 0.d0
        endif
        write(LIVE_DATA_HANDLE,'(ES17.9)',advance='no') growth_rate
      end do
    end if
    write(LIVE_DATA_HANDLE,*)
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@axis: ', xtime(index), R_axis_t(index), Z_axis_t(index), Psi_axis_t(index)
    write(LIVE_DATA_HANDLE,'(A,6ES17.9)') '@Xpoint_low: ', xtime(index), R_xpoint_t(index,1), Z_xpoint_t(index,1), Psi_xpoint_t(index,1)
    write(LIVE_DATA_HANDLE,'(A,6ES17.9)') '@Xpoint_up: ', xtime(index), R_xpoint_t(index,2), Z_xpoint_t(index,2), Psi_xpoint_t(index,2)
    write(LIVE_DATA_HANDLE,'(A,6ES17.9)') '@bnd_point: ', xtime(index), R_bnd_t(index), Z_bnd_t(index), Psi_bnd_t(index)

    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@current: ', xtime(index), Ip_tot_t(index), current_t(index), Ip_tot_t(index)-current_t(index)
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@betas: ', xtime(index), beta_p_t(index), beta_t_t(index), beta_n_t(index)
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@particlecontent: ', xtime(index),density_tot_t(index), density_in_t(index), density_out_t(index), &
                                                                npart_tot_t(index)
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@thermalenergy: ', xtime(index), pressure_in_t(index), pressure_out_t(index)
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@heatingpower: ', xtime(index), heat_src_tot_t(index), heat_src_in_t(index), heat_src_out_t(index)
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@particlesource: ', xtime(index), part_src_tot_t(index), part_src_in_t(index), part_src_out_t(index)
    write(LIVE_DATA_HANDLE,'(A,6ES17.9)') '@integrated_energies: ', xtime(index), E_tot_t(index), Wmag_tot_t(index), &
                                                     kin_par_tot_t(index),  kin_perp_tot_t(index),  thermal_tot_t(index), thermal_e_tot_t(index), thermal_i_tot_t(index) 
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@helicity: ', xtime(index), helicity_tot_t(index)
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@area: ', xtime(index), area_t(index)
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@volume: ', xtime(index), volume_t(index)
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@li3: ', xtime(index), li3_t(index), li3_tot_t(index)

#if (defined WITH_Neutrals) || (defined WITH_Impurities)
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@dissipative_terms: ', xtime(index), ohmic_tot_t(index), viscopar_dissip_tot_t(index), friction_dissip_tot_t(index),&
                                                                  xtime_rad_power(index), xtime_E_ion_power(index)
#else
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@dissipative_terms: ', xtime(index), ohmic_tot_t(index), viscopar_dissip_tot_t(index), friction_dissip_tot_t(index)
#endif

    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@mag_energy_src: ', xtime(index), mag_ener_src_tot(index)
    write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@work_terms: ', xtime(index), Magwork_tot_t(index), thmwork_tot_t(index)
    write(LIVE_DATA_HANDLE,'(A,7ES17.9)') '@bnd_fluxes: ', xtime(index), flux_Pvn_t(index), flux_kinpar_t(index), &
                                           flux_qpar_t(index), flux_qperp_t(index)
    write(LIVE_DATA_HANDLE,'(A,7ES17.9)') '@bnd_particle_fluxes: ', xtime(index), part_flux_Dpar_t(index), part_flux_Dperp_t(index), &
                                           part_flux_vpar_t(index), part_flux_vperp_t(index), -npart_flux_t(index)

   if(index>1) then
     write(LIVE_DATA_HANDLE,'(A,6ES17.9)') '@dEdt: ', xtime(index-1), -dE_tot_dt(index-1), -dWmag_tot_dt(index-1), &
                                            -dthermal_tot_dt(index-1),-dkinperp_tot_dt(index-1),-dkinpar_tot_dt(index-1)
     write(LIVE_DATA_HANDLE,'(A,6ES17.9)') '@dparticles_dt: ', xtime(index-1), dpart_tot_dt(index-1) + dnpart_tot_dt(index-1), &
                                                                      dpart_tot_dt(index-1), dnpart_tot_dt(index-1) 

     sum_fluxes_dissip = flux_Pvn_t(index-1)  + flux_kinpar_t(index-1) + flux_qpar_t(index-1) + flux_qperp_t(index-1) &
                       - heat_src_tot_t(index-1) + ohmic_tot_t(index-1)*(1.d0 - eta_ohmic/eta) - mag_ener_src_tot(index-1) &
                       - flux_poynting_t(index-1)
#if (defined WITH_Neutrals) || (defined WITH_Impurities)
     sum_fluxes_dissip = sum_fluxes_dissip + xtime_rad_power(index-1) + xtime_E_ion_power(index-1)
#endif

     sum_mag_energy_terms = -ohmic_tot_t(index-1) + flux_poynting_t(index-1) + Magwork_tot_t(index-1) + mag_ener_src_tot(index-1) 
 
     write(LIVE_DATA_HANDLE,'(A,6ES17.9)') '@energy_conservation: ', xtime(index-1), -dE_tot_dt(index-1), sum_fluxes_dissip 
     write(LIVE_DATA_HANDLE,'(A,7ES17.9)') '@mag_energy_balance: ', xtime(index-1), dWmag_tot_dt(index-1),  -ohmic_tot_t(index-1),     &
                                                            flux_poynting_t(index-1), Magwork_tot_t(index-1),mag_ener_src_tot(index-1), &
                                                            sum_mag_energy_terms 
    else

      sum_mag_energy_terms = -ohmic_tot_t(index) + flux_poynting_t(index) + Magwork_tot_t(index) + mag_ener_src_tot(index) 

      write(LIVE_DATA_HANDLE,'(A,6ES17.9)') '@dEdt: ', xtime(index), 0.d0, 0.d0, 0.d0, 0.d0, 0.d0
      write(LIVE_DATA_HANDLE,'(A,6ES17.9)') '@energy_conservation: ', xtime(index), 0.d0, 0.d0 
      write(LIVE_DATA_HANDLE,'(A,7ES17.9)') '@mag_energy_balance: ', xtime(index), 0.d0,  -ohmic_tot_t(index),     &
                                                            flux_poynting_t(index), Magwork_tot_t(index),mag_ener_src_tot(index), &
                                                            sum_mag_energy_terms 
      write(LIVE_DATA_HANDLE,'(A,6ES17.9)') '@dparticles_dt: ', xtime(index), 0.d0, 0.d0, 0.d0 
    endif
 
    close(LIVE_DATA_HANDLE)
    
  end subroutine write_live_data
  
  
  
  !> Close file.
  subroutine finalize_live_data()
    
    use phys_module, only: produce_live_data
    
    implicit none
    
    if ( .not. produce_live_data ) return
    
    ! -nothing to be done currently-
    
  end subroutine finalize_live_data
  
  
  
  subroutine write_live_data_vacuum(index)
    
    use phys_module, only: xtime, mu_zero, sqrt_mu0_rho0, Z_axis_t
    use vacuum
      
    integer,             intent(in) :: index
    logical, save :: header_written_diag = .false., header_written_pf = .false., header_written_rmp = .false., &
                     header_written_net  = .false., header_written_VFB = .false.
    integer :: n
    open(LIVE_DATA_HANDLE, file=LIVE_DATA_FILE, status='OLD', position='APPEND', action='WRITE')
    
    if ( allocated(diag_coil_curr) ) then
      if ( .not. header_written_diag ) then
        write(LIVE_DATA_HANDLE,'(A,I5)') '@n_diag_coil_curr: ', size(diag_coil_curr,2)
        write(LIVE_DATA_HANDLE,'(A)') '@diag_coil_curr_xlabel: normalized time'
        write(LIVE_DATA_HANDLE,'(A)') '@diag_coil_curr_xlabel_si: time [ms]'
        write(LIVE_DATA_HANDLE,'(A)') '@diag_coil_curr_ylabel: Diagnostic coil current'
        write(LIVE_DATA_HANDLE,'(A)') '@diag_coil_curr_ylabel_si: Diagnostic coil current [A]'
        write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@diag_coil_curr_x2si: ', sqrt_mu0_rho0*1.e3
        write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@diag_coil_curr_y2si: ', 1./mu_zero
        write(LIVE_DATA_HANDLE,'(A)') '@diag_coil_curr_logy: 0'
        write(LIVE_DATA_HANDLE,'(A)',advance='no') '@diag_coil_curr: %"time"           '
        do n = 1,size(diag_coil_curr,2)
          write(LIVE_DATA_HANDLE,'(A12,1x)',advance='no') trim(diag_coil_name(n))
        end do
        write(LIVE_DATA_HANDLE,*)
        header_written_diag = .true.
      end if
      write(LIVE_DATA_HANDLE,'(A,999ES17.9)') '@diag_coil_curr: ', xtime(index), diag_coil_curr(index,:)
    end if

    if ( allocated(pf_coil_curr) ) then
      if ( .not. header_written_pf ) then
        write(LIVE_DATA_HANDLE,'(A,I5)') '@n_pf_coil_curr: ', size(pf_coil_curr,2)
        write(LIVE_DATA_HANDLE,'(A)') '@pf_coil_curr_xlabel: normalized time'
        write(LIVE_DATA_HANDLE,'(A)') '@pf_coil_curr_xlabel_si: time [ms]'
        write(LIVE_DATA_HANDLE,'(A)') '@pf_coil_curr_ylabel: PF coil current'
        write(LIVE_DATA_HANDLE,'(A)') '@pf_coil_curr_ylabel_si: PF coil current [A]'
        write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@pf_coil_curr_x2si: ', sqrt_mu0_rho0*1.e3
        write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@pf_coil_curr_y2si: ', 1./mu_zero
        write(LIVE_DATA_HANDLE,'(A)') '@pf_coil_curr_logy: 0'
        write(LIVE_DATA_HANDLE,'(A)',advance='no') '@pf_coil_curr: %"time"           '
        do n = 1,size(pf_coil_curr,2)
          write(LIVE_DATA_HANDLE,'(A12,1x)',advance='no') trim(pf_coil_name(n))
        end do
        write(LIVE_DATA_HANDLE,*)
        header_written_pf = .true.
      end if
      write(LIVE_DATA_HANDLE,'(A,999ES17.9)') '@pf_coil_curr: ', xtime(index), pf_coil_curr(index,:)
    end if

    if ( allocated(rmp_coil_curr) ) then
      if ( .not. header_written_rmp ) then
        write(LIVE_DATA_HANDLE,'(A,I5)') '@n_rmp_coil_curr: ', size(rmp_coil_curr,2)
        write(LIVE_DATA_HANDLE,'(A)') '@rmp_coil_curr_xlabel: normalized time'
        write(LIVE_DATA_HANDLE,'(A)') '@rmp_coil_curr_xlabel_si: time [ms]'
        write(LIVE_DATA_HANDLE,'(A)') '@rmp_coil_curr_ylabel: RMP coil current'
        write(LIVE_DATA_HANDLE,'(A)') '@rmp_coil_curr_ylabel_si: RMP coil current [A]'
        write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@rmp_coil_curr_x2si: ', sqrt_mu0_rho0*1.e3
        write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@rmp_coil_curr_y2si: ', 1./mu_zero
        write(LIVE_DATA_HANDLE,'(A)') '@rmp_coil_curr_logy: 0'
        write(LIVE_DATA_HANDLE,*)
        write(LIVE_DATA_HANDLE,'(A)',advance='no') '@RMP_coil_curr: %"time"           '
        do n = 1,size(rmp_coil_curr,2)
          write(LIVE_DATA_HANDLE,'(A12,1x)',advance='no') trim(rmp_coil_name(n))
        end do
        header_written_rmp = .true.
      end if
      write(LIVE_DATA_HANDLE,'(A,999ES17.9)') '@rmp_coil_curr: ', xtime(index), rmp_coil_curr(index,:)
    end if

    if ( allocated(net_tor_wall_curr) ) then
      if ( .not. header_written_net ) then
        write(LIVE_DATA_HANDLE,'(A,I5)') '@n_net_tor_wall_curr: ', 1
        write(LIVE_DATA_HANDLE,'(A)') '@net_tor_wall_curr_xlabel: normalized time'
        write(LIVE_DATA_HANDLE,'(A)') '@net_tor_wall_curr_xlabel_si: time [ms]'
        write(LIVE_DATA_HANDLE,'(A)') '@net_tor_wall_curr_ylabel: Net toroidal wall current'
        write(LIVE_DATA_HANDLE,'(A)') '@net_tor_wall_curr_ylabel_si: Net toroidal wall current [A]'
        write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@net_tor_wall_curr_x2si: ', sqrt_mu0_rho0*1.e3
        write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@net_tor_wall_curr_y2si: ', 1./mu_zero
        write(LIVE_DATA_HANDLE,'(A)') '@net_tor_wall_curr_logy: 0'
        write(LIVE_DATA_HANDLE,'(A)') '@net_tor_wall_curr: %"time"           "I_{tor,wall}"'
        header_written_net = .true.
      end if
      write(LIVE_DATA_HANDLE,'(A,999ES17.9)') '@net_tor_wall_curr: ', xtime(index), net_tor_wall_curr(index)
    end if

    if ( allocated(vert_FB_response) ) then
      if ( .not. header_written_VFB ) then
        write(LIVE_DATA_HANDLE,'(A,I5)') '@n_vert_FB_response: ', 3
        write(LIVE_DATA_HANDLE,'(A)') '@vert_FB_response_xlabel: normalized time'
        write(LIVE_DATA_HANDLE,'(A)') '@vert_FB_response_xlabel_si: time [ms]'
        write(LIVE_DATA_HANDLE,'(A)') '@vert_FB_response_ylabel: VFB response [-]'
        write(LIVE_DATA_HANDLE,'(A)') '@vert_FB_response_ylabel_si: VFB response [-]'
        write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@vert_FB_response_x2si: ', sqrt_mu0_rho0*1.e3
        write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@vert_FB_response_y2si: ', 1.
        write(LIVE_DATA_HANDLE,'(A)') '@vert_FB_response_logy: 0'
        write(LIVE_DATA_HANDLE,'(A)') '@vert_FB_response: %"time"           "proportional_FB"              "derivative_FB"           "integral_FB"'

        write(LIVE_DATA_HANDLE,'(A,I5)') '@n_vert_FB_axis: ', 2
        write(LIVE_DATA_HANDLE,'(A)') '@vert_FB_axis_xlabel: normalized time'
        write(LIVE_DATA_HANDLE,'(A)') '@vert_FB_axis_xlabel_si: time [ms]'
        write(LIVE_DATA_HANDLE,'(A)') '@vert_FB_axis_ylabel: Z-axis [m]'
        write(LIVE_DATA_HANDLE,'(A)') '@vert_FB_axis_ylabel_si: Z-axis [m]'
        write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@vert_FB_axis_x2si: ', sqrt_mu0_rho0*1.e3
        write(LIVE_DATA_HANDLE,'(A,5ES17.9)') '@vert_FB_axis_y2si: ', 1.
        write(LIVE_DATA_HANDLE,'(A)') '@vert_FB_axis_logy: 0'
        write(LIVE_DATA_HANDLE,'(A)') '@vert_FB_axis: %"time"           "Z-axis"              "Z-axis,target"'
        header_written_VFB = .true.
      end if
      write(LIVE_DATA_HANDLE,'(A,999ES17.9)') '@vert_FB_response: ', xtime(index), vert_FB_response(index,1:3)
      write(LIVE_DATA_HANDLE,'(A,999ES17.9)') '@vert_FB_axis:     ', xtime(index), Z_axis_t(index), vert_FB_response(index,4)
    end if
   
    close(LIVE_DATA_HANDLE)
    
  end subroutine write_live_data_vacuum
  
end module live_data
