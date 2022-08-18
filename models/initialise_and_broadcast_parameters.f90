!> Initialize parameters and broadcast them to all MPI procs.
subroutine initialise_and_broadcast_parameters(my_id, filename)
  
  use constants, only: mu_zero
  use mod_parameters,  only: n_tor, n_period
  use phys_module
  
  implicit none
  
  ! --- Routine parameters
  integer,                      intent(in) :: my_id
  real*8                                   :: rho0, Te0_keV, Ti0_keV, lnA
  character(len=*),             intent(in) :: filename
  
  call initialise_parameters(my_id, filename)
  
  ! --- Broadcast input parameters from MPI thread 0 to the others.
  call broadcast_phys(my_id)
  
  ! --- Broadcast numerical input profiles from MPI thread 0 to the others.
  call broadcast_num_profiles(my_id)
  
  ! --- Initialize the time-stepping parameters.
  call update_time_evol_params()
  
  ! --- Calculate normalization factors.
  rho0               = central_density * 1.d20 * central_mass * mass_proton
  sqrt_mu0_rho0      = sqrt( mu_zero * rho0 )
  sqrt_mu0_over_rho0 = sqrt( mu_zero / rho0 )

  ! --- Calculate nominal parameters printed in the logfile for reference
  if (with_TiTe) then
    Te0_keV               = Te_0 / ( EL_CHG * mu_zero * central_density * 1.d+20 ) / 1.d+3
    Ti0_keV               = Ti_0 / ( EL_CHG * mu_zero * central_density * 1.d+20 ) / 1.d+3
    lnA                   = 14.9 - 0.5*log( central_density ) + log( Te0_keV )

    ZK_e_par_SpitzerHaerm = 4.83d+0 * central_mass*mass_proton/(mass_electron*lnA) * Te0_keV**(2.5d+0) * (gamma-1.d0) * sqrt_mu0_over_rho0
    ZK_i_par_SpitzerHaerm = 5.11d+2 * sqrt(central_mass/2.d+0)/(lnA)               * Ti0_keV**(2.5d+0) * (gamma-1.d0) * sqrt_mu0_over_rho0
  else
    Te0_keV               = T_0 / 2.d+0 / ( EL_CHG * mu_zero * central_density * 1.d+20 ) / 1.d+3
    lnA                   = 14.9 - 0.5*log( central_density ) + log( Te0_keV )

    ZK_par_SpitzerHaerm   = 4.83d+0 * central_mass*mass_proton/(mass_electron*lnA) * Te0_keV**(2.5d+0) * (gamma-1.d0) * sqrt_mu0_over_rho0
  end if
  tauIC_nominal      = central_mass * mass_proton / ( EL_CHG * F0 * sqrt_mu0_rho0 * 2.d0 )
  eta_Spitzer        = ( 1.65d-9 * lnA * Te0_keV**(-1.5d+0) ) / sqrt_mu0_over_rho0

  ! --- Assign minimum values for parallel conduction if not given
  if (T_min_ZKpar  < -1.d10) T_min_ZKpar  = T_min   
  if (Ti_min_ZKpar < -1.d10) Ti_min_ZKpar = T_min   
  if (Te_min_ZKpar < -1.d10) Te_min_ZKpar = T_min   
  
  ! --- Deprecated input parameters ---
  if ( use_murge ) then
    write(*,*) 'ERROR: use_murge=.true. is not supported any more. Remove this parameter from the namelist input file.'
    stop
  else if ( use_murge_element ) then
    write(*,*) 'ERROR: use_murge_element=.true. is not supported any more. Remove this parameter from the namelist input file.'
    stop
  end if
  ! -----------------------------------
  ! -- Set equilibrium solver if not defined by user --
  if ((.not.use_mumps_eq).and.(.not.use_pastix_eq).and.(.not.use_strumpack_eq)) then
#ifdef USE_COMPLEX_PRECOND
    use_mumps_eq = .true.
    use_pastix_eq = .false.
    use_strumpack_eq = .false.
#else
    use_mumps_eq = use_mumps
    use_pastix_eq = use_pastix
    use_strumpack_eq = use_strumpack
#endif
  endif
  ! -----------------------------------
  
  prev_FB_fact = 1.d0 ! needed to make sure current_FB_fact is applied correctly in import_restart
  
end subroutine initialise_and_broadcast_parameters
