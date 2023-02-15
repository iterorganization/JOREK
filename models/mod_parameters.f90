module mod_parameters

use mod_settings
use mod_model_settings

implicit none



contains



! --- Variable names
character(len=64) pure elemental function variable_names(i_var)

implicit none

integer,           intent(in) :: i_var !< Variable index
!logical, optional, intent(in) :: long  !< Long output with description?

character(len=20) :: name
character(len=64) :: descr
logical :: long2

long2 = .false.

if ( (i_var < 1) .or. (i_var > n_var) ) then
  name  = 'NONE'
  descr = 'Illegal variable index!'
  return
else if ( i_var == var_psi ) then
  name  = 'Psi'
  descr = 'poloidal magnetic flux'
else if ( i_var == var_u ) then
  name  = 'u'
  descr = 'velocity stream function'
else if ( i_var == var_zj ) then
  name  = 'zj'
  descr = 'toroidal plasma current'
else if ( i_var == var_w ) then
  name  = 'omega'
  descr = 'toroidal vorticity'
else if ( i_var == var_rho ) then
  name  = 'rho'
  descr = 'particle density'
else if ( i_var == var_T ) then
  name  = 'T'
  descr = 'temperature; sum of electrons and ions'
else if ( i_var == var_Ti ) then
  name  = 'T_i'
  descr = 'ion temperature'
else if ( i_var == var_Te ) then
  name  = 'T_e'
  descr = 'electron temperature'
else if ( i_var == var_Vpar ) then
  name  = 'v_par'
  descr = 'parallel velocity'
else if ( i_var == var_rhon ) then
  name  = 'rho_n'
  descr = 'neutral particle density'
else if ( i_var == var_rhoimp ) then
  name  = 'rho_imp'
  descr = 'impurity density'
else if ( i_var == var_nre ) then
  name  = 'n_RE'
  descr = 'runaway electron particle density'
else if ( i_var == var_AR ) then
  name  = "AR"
  descr = "R component of magnetic vector potential"
else if ( i_var == var_AZ ) then
  name  = "AZ"
  descr = "Z component of magnetic vector potential"
else if ( i_var == var_uR ) then
  name  = "UR"
  descr = "Velocity in R direction"
else if ( i_var == var_uZ ) then
  name  = "UZ"
  descr = "Velocity in Z direction"
else if ( i_var == var_up ) then
  name  = "Uphi"
  descr = "Velocity in phi direction"
else
  write(name, "(I4)") i_var
  descr = "Please check implementation"
end if

if ( .not. long2 ) then
  variable_names = trim(name)
else
  variable_names = trim(name)//' ('//trim(descr)//')'
end if

end function variable_names



#ifdef MODEL_FAMILY

!> Are all selected extensions available and compatible?
logical pure function ext_correct()

implicit none

! --- Local variables
integer :: iext1, iext2

! --- Preset to .true.
ext_correct = .true.

! --- Check binary compatibility of all extensions
do iext1 = 1, n_mod_ext
  if ( .not. with_ext(iext1) ) cycle
  ext_correct = ext_correct .and. ext_available(iext1)
  do iext2 = 1, n_mod_ext
    if ( .not. with_ext(iext2) ) cycle
      ext_correct = ext_correct .and. ext_compatible(iext1, iext2)
  end do
end do

end function ext_correct

#endif



end module mod_parameters
