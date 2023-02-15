module mumps_module
#ifdef USE_MUMPS  
  save
  

  include 'dmumps_struc.h'        ! MUMPS include files defining its datastructure
  
  type (DMUMPS_STRUC) :: mumps_par
  
#endif
end module mumps_module
