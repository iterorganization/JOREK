!> Basic model-dependend hard-coded run parameters.
module mod_model_settings

  implicit none

! ##################################################################################################
! ####  @USERS: This file should not be modified ###################################################
! ##################################################################################################

  integer, parameter :: jorek_model    = 502       !< JOREK physics model

  logical, parameter :: hydrodynamics   = .false.
  logical, parameter :: reduced_MHD     = .true.
  logical, parameter :: full_MHD        = .false.

  logical, parameter :: with_TiTe       = .true.
  logical, parameter :: with_neutrals   = .false.
  logical, parameter :: with_impurities = .true.
  logical, parameter :: with_Vpar       = .true.
  logical, parameter :: with_refluid    = .false.

  integer, parameter :: n_mod_ext            = 0 !< this model is not a model family => no extensions

  integer, parameter :: n_var    = 9

  integer, parameter :: var_A3   = 0                       ! place of variable psi/mag pot 3               (ps or A3)
  integer, parameter :: var_AR   = 0                       ! place of variable mag pot  1                  (AR)
  integer, parameter :: var_AZ   = 0                       ! place of variable mag pot  2                  (AZ)
  integer, parameter :: var_uR   = 0                       ! place of variable velocity 1                  (UR)
  integer, parameter :: var_uZ   = 0                       ! place of variable velocity 2                  (UZ)
  integer, parameter :: var_up   = 0                       ! place of variable velocity 3                  (Up)
  integer, parameter :: var_rho  = 5                       ! place of variable density                     (r or rho)
  integer, parameter :: var_T    = 0                       ! place of variable temperature                 (T)
  integer, parameter :: var_psi  = 1                       ! place of variable psi/mag pot 3               (ps or A3)
  integer, parameter :: var_u    = 2                       ! place of variable velocity stream function    (u)
  integer, parameter :: var_zj   = 3                       ! place of variable toroidal current density    (zj)
  integer, parameter :: var_w    = 4                       ! place of variable vorticity                   (w)
  integer, parameter :: var_Vpar = 7                       ! place of variable parallel velocity           (Vpar)
  integer, parameter :: var_rhon = 8                       ! place of variable neutral or impurity density (rn)
  integer, parameter :: var_Ti   = 6                       ! place of variable ion temperature             (Ti)
  integer, parameter :: var_Te   = 9                       ! place of variable electron temperature        (Te)
  integer, parameter :: var_jec  = 0                       ! place of variable ECCD current                (jec)
  integer, parameter :: var_jec1 = 0                       ! place of variable ECCD current #1             (jec1)
  integer, parameter :: var_jec2 = 0                       ! place of variable ECCD current #2             (jec2)
  integer, parameter :: var_rhoimp = 8                     ! place of variable for impurity density        (rhoimp)

!  integer, parameter :: n_var          = 9         !< number of variables
!  integer, parameter :: n_dim          = 2         !< number of dimensions
!  integer, parameter :: n_order        = 3         !< order of the polynomial basis
!  integer, parameter :: n_tor          = 11        !< number of toroidal harmonics in physics variables
!  integer, parameter :: n_coord_tor    = 1         !< number of toroidal harmonics in (R, Z) coordinates
!  integer, parameter :: n_period       = 1         !< periodicity in toroidal direction
!  integer, parameter :: n_plane        = 32        !< number of toroidal angles
!  integer, parameter :: n_vertex_max   = 4         !< maximum number of corners of an element
!  integer, parameter :: n_nodes_max    = 60001     !< maximum number of nodes
!  integer, parameter :: n_elements_max = 60001     !< maximum number of elements
!  integer, parameter :: n_boundary_max = 1001      !< maximum number of boundary elements
!  integer, parameter :: n_pieces_max   = 6001      !< maximum number of line pieces describing a flux surface
!  integer, parameter :: n_degrees      = n_order+1 !< degrees of freedom per variable per node
!  integer, parameter :: nref_max       = 10000     !< (refinement)
!  integer, parameter :: n_ref_list     = 10000     !< (refinement)
! 
!  !> Names of the physical variables
!  character(len=11) :: variable_names(n_var) =                       &
!    (/ 'Flux       ','Potential  ','Current    ','Vorticity  ',      &
!       'Density    ','Ions_Temp  ','V_parallel ','N_density  ','Elec_Temp  '/)

  integer, parameter :: var_nre  = 0                       ! place of variable for RE number density       (nre)

  !> element_matrix and element_matrix_fft combined into a single one?
  logical, parameter :: unified_element_matrix = .false.

end module mod_model_settings
