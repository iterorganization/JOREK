!> Basic model-dependend hard-coded run parameters.
module mod_model_settings

  implicit none

! ##################################################################################################
! ####  @USERS: This file should not be modified ###################################################
! ##################################################################################################

  integer, parameter :: jorek_model    = 199       !< JOREK physics model

  logical, parameter :: hydrodynamics   = .false.
  logical, parameter :: reduced_MHD     = .true.
  logical, parameter :: full_MHD        = .false.
  
  logical, parameter :: with_TiTe       = .false.
  logical, parameter :: with_neutrals   = .false.
  logical, parameter :: with_impurities = .false.
  logical, parameter :: with_Vpar       = .false.
  logical, parameter :: with_refluid    = .false.

  integer, parameter :: n_mod_ext            = 0 !< this model is not a model family => no extensions

  integer, parameter :: n_var    = 6
  
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
  integer, parameter :: var_Vpar = 0                       ! place of variable parallel velocity           (Vpar)
  integer, parameter :: var_rhon = 0                       ! place of variable neutral or impurity density (rn)
  integer, parameter :: var_Ti   = 0                       ! place of variable ion temperature             (Ti)
  integer, parameter :: var_Te   = 0                       ! place of variable electron temperature        (Te)
  integer, parameter :: var_jec  = 0                       ! place of variable ECCD current                (jec)
  integer, parameter :: var_jec1 = 0                       ! place of variable ECCD current #1             (jec1)
  integer, parameter :: var_jec2 = 0                       ! place of variable ECCD current #2             (jec2)
  integer, parameter :: var_nre  = 0                       ! place of variable for RE number density       (nre)
  integer, parameter :: var_rhoimp = 0                     ! place of variable for impurity density        (rhoimp)

  !> element_matrix and element_matrix_fft combined into a single one?
  logical, parameter :: unified_element_matrix = .false.

end module mod_model_settings
