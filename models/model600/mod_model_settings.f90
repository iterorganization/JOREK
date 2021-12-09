!> Basic model-dependend hard-coded run parameters.
module mod_model_settings

implicit none

logical, parameter :: with_vpar       = .true.
logical, parameter :: with_TiTe       = .false.
logical, parameter :: with_neutrals   = .false. 
logical, parameter :: with_impurities = .false. ! not yet possible to switch
logical, parameter :: with_refluid    = .false. ! not yet possible to switch


! ##################################################################################################
! ####  @USERS: Please do not change below this line ###############################################
! ##################################################################################################


! The following line is needed by ./util/config.sh:
! #SETTINGS# with_vpar with_TiTe with_neutrals

integer, parameter :: jorek_model     = 600

logical, parameter :: hydrodynamics   = .false.
logical, parameter :: reduced_MHD     = .true.
logical, parameter :: full_MHD        = .false.

logical, parameter :: model_family      = .true. !< Is this 
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
integer, parameter :: n_var_neutrals    = sum(merge( (/1/), (/0/), with_neutrals  )) !### not yet
integer, parameter :: n_var_impurities  = sum(merge( (/1/), (/0/), with_impurities)) !### not yet
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
integer, parameter :: var_nre      = sum(merge((/n_var_base+sum(n_var_ext(1:3))+1/), (/0/), with_refluid   ))
integer, parameter :: var_rho_imp  = sum(merge((/n_var_base+sum(n_var_ext(1:4))+1/), (/0/), with_impurities))
! --- variables not relevant to this model
integer, parameter :: var_AR   = 0
integer, parameter :: var_AZ   = 0
integer, parameter :: var_A3   = 0
integer, parameter :: var_uR   = 0
integer, parameter :: var_uZ   = 0
integer, parameter :: var_up   = 0
! --- Element matrix and element matrix fft combined?
logical, parameter :: unified_element_matrix = .true.



contains



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
    ext_available = .false.
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
