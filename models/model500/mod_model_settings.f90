!> Basic model-dependend hard-coded run parameters.
module mod_model_settings

  implicit none

! ##################################################################################################
! ####  @USERS: This file should not be modified ###################################################
! ##################################################################################################

  integer, parameter :: jorek_model    = 500       !< JOREK physics model

  logical, parameter :: hydrodynamics   = .false.
  logical, parameter :: reduced_MHD     = .true.
  logical, parameter :: full_MHD        = .false.
  
  logical, parameter :: with_TiTe       = .false.
  logical, parameter :: with_neutrals   = .true.
  logical, parameter :: with_impurities = .false.
  logical, parameter :: with_Vpar       = .true.
  logical, parameter :: with_refluid    = .false.

  integer, parameter :: n_mod_ext            = 0 !< this model is not a model family => no extensions

  integer, parameter :: n_var    = 8
  
  integer, parameter :: var_A3   = 0                       ! place of variable psi/mag pot 3               (ps or A3)
  integer, parameter :: var_AR   = 0                       ! place of variable mag pot  1                  (AR)
  integer, parameter :: var_AZ   = 0                       ! place of variable mag pot  2                  (AZ)
  integer, parameter :: var_uR   = 0                       ! place of variable velocity 1                  (UR)
  integer, parameter :: var_uZ   = 0                       ! place of variable velocity 2                  (UZ)
  integer, parameter :: var_up   = 0                       ! place of variable velocity 3                  (Up)
  integer, parameter :: var_rho  = 5                       ! place of variable density                     (r or rho)
  integer, parameter :: var_T    = 6                       ! place of variable temperature                 (T)
  integer, parameter :: var_psi  = 1                       ! place of variable psi/mag pot 3               (ps or A3)  
  integer, parameter :: var_u    = 2                       ! place of variable velocity stream function    (u)
  integer, parameter :: var_zj   = 3                       ! place of variable toroidal current density    (zj)
  integer, parameter :: var_w    = 4                       ! place of variable vorticity                   (w)
  integer, parameter :: var_Vpar = 7                       ! place of variable parallel velocity           (Vpar)
  integer, parameter :: var_rhon = 8                       ! place of variable neutral or impurity density (rn)
  integer, parameter :: var_Ti   = 0                       ! place of variable ion temperature             (Ti)
  integer, parameter :: var_Te   = 0                       ! place of variable electron temperature        (Te)
  integer, parameter :: var_jec  = 0                       ! place of variable ECCD current                (jec)
  integer, parameter :: var_jec1 = 0                       ! place of variable ECCD current #1             (jec1)
  integer, parameter :: var_jec2 = 0                       ! place of variable ECCD current #2             (jec2)
  integer, parameter :: var_nre  = 0                       ! place of variable for RE number density       (nre)
  integer, parameter :: var_rhoimp = 0                     ! place of variable for impurity density        (rhoimp)

  !> element_matrix and element_matrix_fft combined into a single one?
  logical, parameter :: unified_element_matrix = .true.

  !> parameters for naming equation terms in the RHS diagnostic 
  integer,  parameter :: max_terms    = 20
  integer,  parameter :: n_terms_psi  = 5
  integer,  parameter :: n_terms_u    = 13 
  integer,  parameter :: n_terms_zj   = 1
  integer,  parameter :: n_terms_w    = 1
  integer,  parameter :: n_terms_rho  = 12
  integer,  parameter :: n_terms_T    = 17
  integer,  parameter :: n_terms_vpar = 10
  integer,  parameter :: n_terms_rhon = 7

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
                                                   'u_Eq__ext_dens_source  ', &  !  9:
                                                   'u_Eq__ionization_term  ', &  ! 10:
                                                   'u_Eq__recombin_term    ', &  ! 11:
                                                   'u_Eq__zeta_timevol_term', &  ! 12:
                                                   'u_Eq__neoclassical_term'/)   ! 13:

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
                                                   'T_Eq__Vperp.grad_p     ', &  !  2: 
                                                   'T_Eq__gamma_P_div_Vperp', &  !  3:
                                                   'T_Eq__Vpar.grad_P      ', &  !  4:
                                                   'T_Eq__gamma_P_div_Vpar ', &  !  5:
                                                   'T_Eq__parallel_conduct ', &  !  6:
                                                   'T_Eq__perp_conduction  ', &  !  7:
                                                   'T_Eq__ZK_perp_num_term ', &  !  8:
                                                   'T_Eq__ZK_par_num_term  ', &  !  9:
                                                   'T_Eq__ionization_sink  ', &  ! 10:
                                                   'T_Eq__neutral_friction ', &  ! 11:
                                                   'T_Eq__ohmic_heating    ', &  ! 12:
                                                   'T_Eq__line_radiation   ', &  ! 13:
                                                   'T_Eq__Brems_radiation  ', &  ! 14:
                                                   'T_Eq__backg_imp_radiat ', &  ! 15:
                                                   'T_Eq__tg_num_terms     ', &  ! 16:
                                                   'T_Eq__zeta_timevol_term'/)   ! 17:

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
      else if (k_var == var_vpar) then
        term_names(k_var, 1:n_terms_vpar) = vpar_term_names(:)
      else if (k_var == var_rhon) then
        term_names(k_var, 1:n_terms_rhon) = rhon_term_names(:)
      endif

    enddo

  end subroutine assign_term_names 

end module mod_model_settings
