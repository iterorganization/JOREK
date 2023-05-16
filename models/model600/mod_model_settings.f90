!> Basic model-dependend hard-coded run parameters.
module mod_model_settings

implicit none

logical, parameter :: with_vpar       = .true.
logical, parameter :: with_TiTe       = .false.
logical, parameter :: with_neutrals   = .false. 
logical, parameter :: with_impurities = .false.
logical, parameter :: with_refluid    = .false. ! not yet possible to switch


! ##################################################################################################
! ####  @USERS: Please do not change below this line ###############################################
! ##################################################################################################


! The following line is needed by ./util/config.sh:
! #SETTINGS# with_vpar with_TiTe with_neutrals with_impurities

integer, parameter :: jorek_model     = 600

logical, parameter :: hydrodynamics   = .false.
logical, parameter :: reduced_MHD     = .true.
logical, parameter :: full_MHD        = .false.

logical, parameter :: model_family      = .true.
character(len=42)  :: base_mod_descr    = 'Model family for tokamak reduced MHD'

! --- extensions to it
integer, parameter :: n_mod_ext            = 5      !< Number of model extensions
integer, parameter :: i_ext_TiTe           = 1
integer, parameter :: i_ext_vpar           = 2
integer, parameter :: i_ext_neutrals       = 3
integer, parameter :: i_ext_impurities     = 4
integer, parameter :: i_ext_refluid        = 5
logical, parameter :: with_ext(n_mod_ext) = &
  (/ with_TiTe, with_vpar, with_neutrals, with_impurities, with_refluid /)

! --- number of variables (for base model, extension, and in total)
integer, parameter :: n_var_base        = 6         !< number of variables in base model
integer, parameter :: n_var_TiTe        = sum(merge( (/1/), (/0/), with_TiTe      ))
integer, parameter :: n_var_vpar        = sum(merge( (/1/), (/0/), with_vpar      ))
integer, parameter :: n_var_neutrals    = sum(merge( (/1/), (/0/), with_neutrals  ))
integer, parameter :: n_var_impurities  = sum(merge( (/1/), (/0/), with_impurities))
integer, parameter :: n_var_refluid     = sum(merge( (/1/), (/0/), with_refluid   )) !### not yet
integer, parameter :: n_var_ext(n_mod_ext) = (/ n_var_TiTe, n_var_vpar, n_var_neutrals,         &
  n_var_impurities, n_var_refluid /)
integer, parameter :: n_var = n_var_base + sum(n_var_ext) !< total number of variables

! --- variable indices for the base model
integer, parameter :: var_psi  = 1
integer, parameter :: var_u    = 2
integer, parameter :: var_zj   = 3
integer, parameter :: var_w    = 4
integer, parameter :: var_rho  = 5
integer, parameter :: var_T    = sum(merge( (/0/), (/6/), with_TiTe ))
! --- variable indices for the model extensions
integer, parameter :: var_Ti       = sum(merge((/                               6/), (/0/), with_TiTe      ))
integer, parameter :: var_Te       = sum(merge((/n_var_base+n_var_vpar         +1/), (/0/), with_TiTe      ))
integer, parameter :: var_Vpar     = sum(merge((/                               7/), (/0/), with_vpar      ))
integer, parameter :: var_rhon     = sum(merge((/n_var_base+sum(n_var_ext(1:2))+1/), (/0/), with_neutrals  ))
integer, parameter :: var_rhoimp   = sum(merge((/n_var_base+sum(n_var_ext(1:3))+1/), (/0/), with_impurities))
integer, parameter :: var_nre      = sum(merge((/n_var_base+sum(n_var_ext(1:4))+1/), (/0/), with_refluid   ))
! --- variables not relevant to this model
integer, parameter :: var_AR   = 0
integer, parameter :: var_AZ   = 0
integer, parameter :: var_A3   = 0
integer, parameter :: var_uR   = 0
integer, parameter :: var_uZ   = 0
integer, parameter :: var_up   = 0
! --- Element matrix and element matrix fft combined?
logical, parameter :: unified_element_matrix = .true.

!> parameters for naming equation terms in the RHS diagnostic 
integer,  parameter :: max_terms    = 20
integer,  parameter :: n_terms_psi  = 5
integer,  parameter :: n_terms_u    = 11
integer,  parameter :: n_terms_zj   = 1
integer,  parameter :: n_terms_w    = 1
integer,  parameter :: n_terms_rho  = 12
integer,  parameter :: n_terms_T    = 20
integer,  parameter :: n_terms_Te   = 19
integer,  parameter :: n_terms_Ti   = 13
integer,  parameter :: n_terms_vpar = 10
integer,  parameter :: n_terms_rhon = 7
integer,  parameter :: n_terms_rhoimp = 10

character*36, dimension(n_var, max_terms) :: term_names
character*36, dimension(n_terms_psi),  parameter :: Psi_term_names=  &
                                              (/ 'psi_Eq__eta_J            ', &  ! 1: \eta j
                                                 'psi_Eq__B.grad_u         ', &  ! 2: \mathbf{B}\cdot\nabla u 
                                                 'psi_Eq__eta_num_term     ', &  ! 3: \eta_{num}\nabla^2 j 
                                                 'psi_Eq__diamag_term      ', &  ! 4: \mathbf{B}\cdot\nabla p
                                                 'psi_Eq__zeta_timevol_term'/)   ! 5: \zeta\delta\psi

character*36, dimension(n_terms_u),     parameter :: u_term_names=  &
                                              (/ 'u_Eq__rho_v.grad_v     ', &  !  1:
                                                 'u_Eq__JxB              ', &  !  2: 
                                                 'u_Eq__visco_term       ', &  !  3:
                                                 'u_Eq__grad_p           ', &  !  4:
                                                 'u_Eq__visco_num_term   ', &  !  5:
                                                 'u_Eq__tg_num_term      ', &  !  6:
                                                 'u_Eq__diamag_term      ', &  !  7:
                                                 'u_Eq__diamag_visco     ', &  !  8:
                                                 'u_Eq__zeta_timevol_term', &  !  9:
                                                 'u_Eq__ext_dens_source  ', &  ! 10:
                                                 'u_Eq__neoclassical_term'/)   ! 11:

 character*36, dimension(n_terms_zj),   parameter :: zj_term_names=  &
                                              (/ 'zj_Eq__DeltaStar_Psi   '/)  !  1:

 character*36, dimension(n_terms_w),    parameter :: w_term_names=  &
                                              (/ 'w_Eq__DeltaStar_u      '/)  !  1:

 character*36, dimension(n_terms_rho),  parameter :: rho_term_names=  &
                                              (/ 'rho_Eq__ext_dens_source  ', &  !  1:
                                                 'rho_Eq__perp_convection  ', &  !  2: 
                                                 'rho_Eq__divergence_v     ', &  !  3: 
                                                 'rho_Eq__parallel_diffus  ', &  !  4:
                                                 'rho_Eq__perp_diffusion   ', &  !  5:
                                                 'rho_Eq__parallel_convect ', &  !  6:
                                                 'rho_Eq__diamag_term      ', &  !  7:
                                                 'rho_Eq__ionization_source', &  !  8:
                                                 'rho_Eq__recombination    ', &  !  9:
                                                 'rho_Eq__zeta_time_evol   ', &  ! 10:
                                                 'rho_Eq__Dperp_num_term   ', &  ! 11:
                                                 'rho_Eq__tg_num_term      '/)   ! 12:

character*36, dimension(n_terms_T),     parameter :: T_term_names=  &
                                              (/ 'T_Eq__ext_heat_source  ', &  !  1:
                                                 'T_Eq__Vperp.grad_P     ', &  !  2: 
                                                 'T_Eq__gamma_P_div_V    ', &  !  3:
                                                 'T_Eq__Vpar.grad_P      ', &  !  4:
                                                 'T_Eq__parallel_conduct ', &  !  5:
                                                 'T_Eq__perp_conduction  ', &  !  6:
                                                 'T_Eq__ZK_perp_num_term ', &  !  7:
                                                 'T_Eq__tg_num_terms     ', &  !  8:
                                                 'T_Eq__ohmic_heating    ', &  !  9:
                                                 'T_Eq__zeta_timevol_term', &  ! 10:
                                                 'T_Eq__neutral_friction ', &  ! 11:
                                                 'T_Eq__ionization_sink  ', &  ! 12:
                                                 'T_Eq__line_radiation   ', &  ! 13:
                                                 'T_Eq__Brems_radiation  ', &  ! 14:
                                                 'T_Eq__backg_imp_radiat ', &  ! 15:
                                                 'T_Eq__main_imp_radiat  ', &  ! 16:
                                                 'T_Eq__imp_ionization   ', &  ! 17:
                                                 'T_Eq__power_teleported ', &  ! 18:
                                                 'T_Eq__viscopar_heating ', &  ! 19:
                                                 'T_Eq__impl_heating     '/)   ! 20:

character*36, dimension(n_terms_Ti),    parameter :: Ti_term_names=  &
                                              (/ 'Ti_Eq__ext_heat_source ', &  !  1:
                                                 'Ti_Eq__Vperp.grad_Pi   ', &  !  2: 
                                                 'Ti_Eq__gamma_Pi_div_V  ', &  !  3:
                                                 'Ti_Eq__Vpar.grad_Pi    ', &  !  4:
                                                 'Ti_Eq__parallel_conduct', &  !  5:
                                                 'Ti_Eq__perp_conduction ', &  !  6:
                                                 'Ti_Eq__ZK_perp_num_term', &  !  7:
                                                 'Ti_Eq__tg_num_terms    ', &  !  8:
                                                 'Ti_Eq__zeta_timevol    ', &  !  9:
                                                 'Ti_Eq__neutral_friction', &  ! 10:
                                                 'Ti_Eq__TiTe_energy_exch', &  ! 11:
                                                 'Ti_Eq__viscopar_heating', &  ! 12:
                                                 'Ti_Eq__implicit_heating'/)   ! 13:

character*36, dimension(n_terms_Te),    parameter :: Te_term_names=  &
                                              (/ 'Te_Eq__ext_heat_source ', &  !  1:
                                                 'Te_Eq__Vperp.grad_Pe   ', &  !  2: 
                                                 'Te_Eq__gamma_Pe_div_V  ', &  !  3:
                                                 'Te_Eq__Vpar.grad_Pe    ', &  !  4:
                                                 'Te_Eq__parallel_conduct', &  !  5:
                                                 'Te_Eq__perp_conduction ', &  !  6:
                                                 'Te_Eq__ZK_perp_num_term', &  !  7:
                                                 'Te_Eq__tg_num_terms    ', &  !  8:
                                                 'Te_Eq__ohmic_heating   ', &  !  9:
                                                 'Te_Eq__zeta_timevol_ter', &  ! 10:
                                                 'Ti_Eq__TiTe_energy_exch', &  ! 11:
                                                 'Te_Eq__ionization_sink ', &  ! 12:
                                                 'Te_Eq__line_radiation  ', &  ! 13:
                                                 'Te_Eq__Brems_radiation ', &  ! 14:
                                                 'Te_Eq__backg_imp_radiat', &  ! 15:
                                                 'Te_Eq__main_imp_radiat ', &  ! 16:
                                                 'Te_Eq__imp_ionization  ', &  ! 17:
                                                 'Te_Eq__power_teleported', &  ! 19:
                                                 'Te_Eq__implicit_heating'/)   ! 20:


character*36, dimension(n_terms_vpar),  parameter :: vpar_term_names=  &
                                              (/ 'vpar_Eq__B.grad_P         ', &  !  1:
                                                 'vpar_Eq__ext_part_source  ', &  !  2: 
                                                 'vpar_Eq__B._rho_v.grad_v  ', &  !  3:
                                                 'vpar_Eq__viscopar_num_term', &  !  4:
                                                 'vpar_Eq__zeta_timevol_term', &  !  5:
                                                 'vpar_Eq__tg_num_terms     ', &  !  6:
                                                 'vpar_Eq__ionization_term  ', &  !  7:
                                                 'vpar_Eq__recombin_term    ', &  !  8:
                                                 'vpar_Eq__viscopar_term    ', &  !  9:
                                                 'vpar_Eq__neoclassical_term'/)   ! 10:

 character*36, dimension(n_terms_rhon), parameter :: rhon_term_names=  &
                                              (/ 'rhon_Eq__neutral_diffusion', &  !  1:
                                                 'rhon_Eq__flow_convection  ', &  !  2: 
                                                 'rhon_Eq__ionization_sink  ', &  !  3:
                                                 'rhon_Eq__recombin_source  ', &  !  4:
                                                 'rhon_Eq__ext_neut_source  ', &  !  5:
                                                 'rhon_Eq__Dn_perp_num_term ', &  !  6:
                                                 'rhon_Eq__zeta_timevol_term'/)   !  7:

 character*36, dimension(n_terms_rhoimp), parameter :: rhoimp_term_names=  &
                                              (/ 'rhoimp_Eq__parallel_diffus ', &  !  1:
                                                 'rhoimp_Eq__perp_diffusion  ', &  !  2: 
                                                 'rhoimp_Eq__perp_convection ', &  !  3: 
                                                 'rhoimp_Eq__divergence_v    ', &  !  4: 
                                                 'rhoimp_Eq__parallel_convect', &  !  5:
                                                 'rhoimp_Eq__diamag_term     ', &  !  6:
                                                 'rhoimp_Eq__tg_num_term     ', &  !  7:
                                                 'rhoimp_Eq__ext_dens_source ', &  !  8:
                                                 'rhoimp_Eq__zeta_time_evol  ', &  !  9:
                                                 'rhoimp_Eq__Dn_perp_num_term'/)   ! 10:

contains



subroutine assign_term_names()

  implicit none

  integer :: k_var

  term_names = ''

  do k_var=1, n_var

    if (k_var == var_psi) then
      term_names(k_var, 1:n_terms_psi ) = Psi_term_names(:)
    else if (k_var == var_u   ) then
      term_names(k_var, 1:n_terms_u   ) = u_term_names(:) 
    else if (k_var == var_zj  ) then
      term_names(k_var, 1:n_terms_zj  ) = zj_term_names(:) 
    else if (k_var == var_w   ) then
      term_names(k_var, 1:n_terms_w   ) = w_term_names(:)
    else if (k_var == var_rho ) then
      term_names(k_var, 1:n_terms_rho ) = rho_term_names(:) 
    else if (k_var == var_T   ) then
      term_names(k_var, 1:n_terms_T   ) = T_term_names(:) 
    else if (k_var == var_Ti  ) then
      term_names(k_var, 1:n_terms_Ti  ) = Ti_term_names(:)
    else if (k_var == var_Te  ) then
      term_names(k_var, 1:n_terms_Te  ) = Te_term_names(:)
    else if (k_var == var_vpar) then
      term_names(k_var, 1:n_terms_vpar) = vpar_term_names(:)
    else if (k_var == var_rhon) then
      term_names(k_var, 1:n_terms_rhon) = rhon_term_names(:)
    else if (k_var == var_rhoimp) then
      term_names(k_var, 1:n_terms_rhoimp) = rhoimp_term_names(:)
    endif

  enddo

end subroutine assign_term_names 



!> Is the extension available?
elemental pure logical function ext_available(i_ext)
  
  implicit none
  
  ! --- Function parameters
  integer, intent(in) :: i_ext
  
  
  ! --- Preset to .true. unless invalid extension number
  if ( (i_ext < 1) .or. (i_ext > n_mod_ext) ) then
    ext_available = .false.
  else
    ext_available = .true.
  end if
  
  ! --- exceptions for not implemented
  if ( i_ext == i_ext_TiTe ) then
    ext_available = .true.
  else if ( i_ext == i_ext_vpar ) then
    ext_available = .true.
  else if ( i_ext == i_ext_neutrals ) then
    ext_available = .true.
  else if ( i_ext == i_ext_impurities ) then
    ext_available = .true.
  else if ( i_ext == i_ext_refluid ) then
    ext_available = .false.
  end if
  
end function ext_available



!> Are two extensions compatible?
pure logical function ext_compatible(i_ext1, i_ext2)
  
  implicit none
  
  ! --- Function parameters
  integer, intent(in) :: i_ext1, i_ext2
  
  ! --- Local variables
  integer :: iext1, iext2
  
  ! --- sort extension indices, since compatibility is symmetric
  iext1 = min(i_ext1, i_ext2)
  iext2 = max(i_ext1, i_ext2)
  
  ! --- preset to .true. unless wrong indices were provided
  if ( (iext1 < 1) .or. (iext2 > n_mod_ext) ) then
    ext_compatible = .false.
  else
    ext_compatible = .true.
  end if
  
  ! --- exceptions for compatibility
  if ( ( iext1 == i_ext_TiTe ) .and. ( iext2 == i_ext_refluid ) ) then
    ext_compatible = .false. ! ### just an example
  end if
  
end function ext_compatible


end module mod_model_settings
