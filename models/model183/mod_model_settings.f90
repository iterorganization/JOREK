!> Basic model-dependend hard-coded run parameters.
module mod_model_settings

  implicit none

  logical, parameter :: with_rho        = .true.
  logical, parameter :: with_vpar       = .false.
  logical, parameter :: with_TiTe       = .false. 

! ##################################################################################################
! ####  @USERS: This file should not be modified ###################################################
! ##################################################################################################

  ! The following line is needed by ./util/config.sh:
  ! #SETTINGS# with_vpar with_TiTe

  ! JOREK stellarator reduced MHD model time evolution model: 
  !    See https://www.jorek.eu/wiki/doku.php?id=stellarator_setup
  integer, parameter :: jorek_model    = 183

  logical, parameter :: hydrodynamics   = .false.
  logical, parameter :: reduced_MHD     = .true.
  logical, parameter :: full_MHD        = .false.
  
  logical, parameter :: with_neutrals   = .false.
  logical, parameter :: with_impurities = .false.
  logical, parameter :: with_refluid    = .false.

  logical, parameter :: model_family    = .true.

! --- extensions to it
  integer, parameter :: n_mod_ext       = 2
  integer, parameter :: i_ext_TiTe      = 1
  integer, parameter :: i_ext_vpar      = 2
  logical, parameter :: with_ext(n_mod_ext) = (/ with_TiTe, with_vpar /)

! --- number of variables for base model and in case of extension
  integer, parameter :: n_var_base    = 6                         !< number of variables in model without TiTe
  integer, parameter :: n_var_TiTe    = sum(merge( (/1/), (/0/), with_TiTe) )
  integer, parameter :: n_var_vpar    = sum(merge( (/1/), (/0/), with_vpar) )
  integer, parameter :: n_var_ext(n_mod_ext) = (/ n_var_TiTe, n_var_vpar /)
  integer, parameter :: n_var = n_var_base + sum(n_var_ext)       !< total number of variables
  
! --- variable indices in general
  integer, parameter :: var_A3     = 0                       ! place of variable psi/mag pot 3               (ps or A3)
  integer, parameter :: var_AR     = 0                       ! place of variable mag pot  1                  (AR)
  integer, parameter :: var_AZ     = 0                       ! place of variable mag pot  2                  (AZ)
  integer, parameter :: var_uR     = 0                       ! place of variable velocity 1                  (UR)
  integer, parameter :: var_uZ     = 0                       ! place of variable velocity 2                  (UZ)
  integer, parameter :: var_up     = 0                       ! place of variable velocity 3                  (Up)
  integer, parameter :: var_rho    = 5                       ! place of variable density                     (r or rho)
  integer, parameter :: var_Psi    = 1                       ! place of variable psi/mag pot 3               (ps or A3)  
  integer, parameter :: var_Phi    = 2                       ! place of variable velocity stream function    (Phi)
  integer, parameter :: var_u      = var_Phi                 ! alias for velocity stream function            (u)
  integer, parameter :: var_zj     = 3                       ! place of variable toroidal current density    (zj)
  integer, parameter :: var_w      = 4                       ! place of variable vorticity                   (w)
  integer, parameter :: var_rhon   = 0                       ! place of variable neutral particle density    (rn)
  integer, parameter :: var_rhoimp = 0                       ! place of variable impurity density            (rimp)
  integer, parameter :: var_jec    = 0                       ! place of variable ECCD current                (jec)
  integer, parameter :: var_jec1   = 0                       ! place of variable ECCD current #1             (jec1)
  integer, parameter :: var_jec2   = 0                       ! place of variable ECCD current #2             (jec2)
  integer, parameter :: var_nre    = 0                       ! place of variable for RE number density       (nre)

! --- variable indices dependent on model
  integer, parameter :: var_T    = sum(merge((/                       0/), (/6/), with_TiTe ))  ! place of variable temperature          (T)
  integer, parameter :: var_Ti   = sum(merge((/                       6/), (/0/), with_TiTe ))  ! place of variable ion temperature      (Ti)
  integer, parameter :: var_Te   = sum(merge((/n_var_base+n_var_vpar  +1/), (/0/), with_TiTe ))  ! place of variable electron temperature (Te)
  integer, parameter :: var_Vpar = sum(merge((/                        7/), (/0/), with_vpar ))  ! place of variable parallel velocity    (Vpar)

  !> element_matrix and element_matrix_fft combined into a single one?
  logical, parameter :: unified_element_matrix = .true.

  contains

!> is the extension available?
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

  end function ext_available

!> Are two extensions compatible?
! --- not necessary for only one extension
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

    ! --- exceptions for compatibility (--> at the moment none)

  end function ext_compatible

end module mod_model_settings
