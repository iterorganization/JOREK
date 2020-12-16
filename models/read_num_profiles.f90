!> Read numerical input profiles
subroutine read_num_profiles(my_id)
  
  use phys_module
  use profiles, only: readProf, readProfNeo, derivProf
  use mod_F_profile
  
  implicit none
  
  integer, intent(in) :: my_id
  integer             :: i
  
  num_rho = ( rho_file /= 'none' )
  if ( num_rho .and. ( my_id == 0 ) ) then
    call readProf(num_rho_x, num_rho_y0, num_rho_len, rho_file)
    call check_num_prof(num_rho, num_rho_x, num_rho_y0, num_rho_len, 'rho', check_positive=.true.)
    rho_1 = num_rho_y0(num_rho_len)
    rho_0 = num_rho_y0(1)
    num_rho_y0 = num_rho_y0 - rho_1
  end if
  
  num_rhon = ( rhon_file /= 'none' )
  if ( num_rhon .and. ( my_id == 0 ) ) then
    call readProf(num_rhon_x, num_rhon_y0, num_rhon_len, rhon_file)
    call check_num_prof(num_rhon, num_rhon_x, num_rhon_y0, num_rhon_len, 'rhon',                   &
      check_positive=.true.)
    rhon_1 = num_rhon_y0(num_rhon_len)
    num_rhon_y0 = num_rhon_y0 - rhon_1
  end if
  
  num_T = ( T_file /= 'none' )
  if ( num_T .and. ( my_id == 0 ) ) then
    call readProf(num_T_x, num_T_y0, num_T_len, T_file)
    call check_num_prof(num_T, num_T_x, num_T_y0, num_T_len, 'T', check_positive=.true.)
    T_1 = num_T_y0(num_T_len)
    T_0 = num_T_y0(1)
    num_T_y0 = num_T_y0 - T_1
  end if
  
  num_Te = ( Te_file /= 'none' )
  if ( num_Te .and. ( my_id == 0 ) ) then
    call readProf(num_Te_x, num_Te_y0, num_Te_len, Te_file)
    call check_num_prof(num_Te, num_Te_x, num_Te_y0, num_Te_len, 'Te', check_positive=.true.)
    Te_1 = num_Te_y0(num_Te_len)
    Te_0 = num_Te_y0(1)
    num_Te_y0 = num_Te_y0 - Te_1
  end if
  
  num_Ti = ( Ti_file /= 'none' )
  if ( num_Ti .and. ( my_id == 0 ) ) then
    call readProf(num_Ti_x, num_Ti_y0, num_Ti_len, Ti_file)
    call check_num_prof(num_Ti, num_Ti_x, num_Ti_y0, num_Ti_len, 'Ti', check_positive=.true.)
    Ti_1 = num_Ti_y0(num_Ti_len)
    Ti_0 = num_Ti_y0(1)
    num_Ti_y0 = num_Ti_y0 - Ti_1
  end if
  
  if ( (jorek_model .eq. 400) .or. (jorek_model .eq. 711) ) then
    if ( my_id ==0 ) then
      T_0 = Te_0 + Ti_0
      T_1 = Te_1 + Ti_1
    end if
  end if
  
  num_ffprime = ( ffprime_file /= 'none' )
  if ( num_ffprime .and. ( my_id == 0 ) ) then
    call readProf(num_ffprime_x, num_ffprime_y0, num_ffprime_len, ffprime_file)
    call check_num_prof(num_ffprime, num_ffprime_x, num_ffprime_y0, num_ffprime_len, 'ffprime',    &
      check_positive=.false.)
    FF_0 = num_ffprime_y0(1)
    FF_1 = num_ffprime_y0(num_ffprime_len)
    num_ffprime_y0 = num_ffprime_y0 - FF_1
    ! --- Early allocation needed for full-MHD
    if ( allocated(num_ffprime_y1) ) deallocate( num_ffprime_y1 )
    if ( allocated(num_ffprime_y2) ) deallocate( num_ffprime_y2 )
    call tr_allocate(num_ffprime_y1,1,num_ffprime_len,"num_ffprime_y1",CAT_GRID)
    call tr_allocate(num_ffprime_y2,1,num_ffprime_len,"num_ffprime_y2",CAT_GRID)
  end if
  
  ! --- Special case for F-profile in model710:
  ! --- If there is no F-profile, we create one by integrating the FF' function numerically
  num_Fprofile = ( Fprofile_file /= 'none' )
  if ( num_Fprofile .and. ( my_id == 0 ) ) then
    write(*,*)'*** WARNING ***'
    write(*,*)'*** numerical Fprofile profiles are not allowed by default'
    write(*,*)'*** in order to avoid confusion with the normalisation of'
    write(*,*)'*** the profile with respect to (psi_bnd-psi_axis)'
    write(*,*)'*** If you really want to use an F-profile, you will need'
    write(*,*)'*** to normalise your F-profile as F_profile/(psi_bnd-psi_axis)'
    write(*,*)'*** Aborting...'
    stop
  endif
  if ( ( my_id == 0 ) .and. ( (jorek_model == 710) .or. (jorek_model == 711) ) ) then
    if ( .not. num_Fprofile ) then
      ! --- Numerical integration of FFprime
      call integrate_F_profile()
      ! --- Copy profile from internal profile instead of reading it from file
      num_Fprofile_len = n_Fprofile_internal
      if ( allocated(num_Fprofile_x ) ) deallocate( num_Fprofile_x  )
      if ( allocated(num_Fprofile_y0) ) deallocate( num_Fprofile_y0 )
      if ( allocated(num_Fprofile_y1) ) deallocate( num_Fprofile_y1 )
      if ( allocated(num_Fprofile_y2) ) deallocate( num_Fprofile_y2 )
      if ( allocated(num_Fprofile_y3) ) deallocate( num_Fprofile_y3 )
      call tr_allocate(num_Fprofile_x ,1,n_Fprofile_internal,"num_Fprofile_x",CAT_GRID)
      call tr_allocate(num_Fprofile_y0,1,n_Fprofile_internal,"num_Fprofile_y0",CAT_GRID)
      call tr_allocate(num_Fprofile_y1,1,n_Fprofile_internal,"num_Fprofile_y1",CAT_GRID)
      call tr_allocate(num_Fprofile_y2,1,n_Fprofile_internal,"num_Fprofile_y2",CAT_GRID)
      call tr_allocate(num_Fprofile_y3,1,n_Fprofile_internal,"num_Fprofile_y3",CAT_GRID)
      do i=1,n_Fprofile_internal
        num_Fprofile_x (i) = Fprofile_psi_max * real(i-1)/real(n_Fprofile_internal-1)
        num_Fprofile_y0(i) = Fprofile_internal(i)
      enddo
      call derivProf(num_Fprofile_x, num_Fprofile_y0, num_Fprofile_len, num_Fprofile_y1)
      call derivProf(num_Fprofile_x, num_Fprofile_y1, num_Fprofile_len, num_Fprofile_y2)
      call derivProf(num_Fprofile_x, num_Fprofile_y2, num_Fprofile_len, num_Fprofile_y3)
      ! --- Check that new profile is accurate by deriving again and comparing against input FFprime
      call check_F_profile_accuracy() ! this will abort if error is too large
    else
      call readProf(num_Fprofile_x, num_Fprofile_y0, num_Fprofile_len, Fprofile_file)
      call check_num_prof(num_Fprofile, num_Fprofile_x, num_Fprofile_y0, num_Fprofile_len, 'Fprofile',    &
        check_positive=.false.)
    end if
  end if
  
  num_d_perp = ( d_perp_file /= 'none' )
  if ( num_d_perp .and. ( my_id == 0 ) ) then
    call readProf(num_d_perp_x, num_d_perp_y, num_d_perp_len, d_perp_file)
    call check_num_prof(num_d_perp, num_d_perp_x, num_d_perp_y, num_d_perp_len, 'd_perp',          &
                        check_positive=.true.)
  end if

  num_d_perp_imp = ( d_perp_imp_file /= 'none' )
  if ( num_d_perp_imp .and. ( my_id == 0 ) ) then
    call readProf(num_d_perp_x_imp, num_d_perp_y_imp, num_d_perp_len_imp, d_perp_imp_file)
    call check_num_prof(num_d_perp, num_d_perp_x_imp, num_d_perp_y_imp, num_d_perp_len_imp, &
                        'd_perp_imp', check_positive=.true.)
  end if

  num_zk_perp = ( zk_perp_file /= 'none' )
  if ( num_zk_perp .and. ( my_id == 0 ) ) then
    call readProf(num_zk_perp_x, num_zk_perp_y, num_zk_perp_len, zk_perp_file)
    call check_num_prof(num_zk_perp, num_zk_perp_x, num_zk_perp_y, num_zk_perp_len, 'zk_perp',     &
      check_positive=.true.)
  end if

  num_zk_e_perp = ( zk_e_perp_file /= 'none' )
  if ( num_zk_e_perp .and. ( my_id == 0 ) ) then
    call readProf(num_zk_e_perp_x, num_zk_e_perp_y, num_zk_e_perp_len, zk_e_perp_file)
    call check_num_prof(num_zk_e_perp, num_zk_e_perp_x, num_zk_e_perp_y, num_zk_e_perp_len,        &
      'zk_e_perp', check_positive=.true.)
  end if
  
  num_zk_i_perp = ( zk_i_perp_file /= 'none' )
  if ( num_zk_i_perp .and. ( my_id == 0 ) ) then
    call readProf(num_zk_i_perp_x, num_zk_i_perp_y, num_zk_i_perp_len, zk_i_perp_file)
    call check_num_prof(num_zk_i_perp, num_zk_i_perp_x, num_zk_i_perp_y, num_zk_i_perp_len,        &
      'zk_i_perp', check_positive=.true.)
  end if

  num_neo_file= NEO .and. ( neo_file /= 'none')
  if ( num_neo_file .and. ( my_id == 0 ) ) then
    write(*,*) 'using ki and mui profiles from file "'//trim(neo_file)//'"'
    call readProfNeo(num_neo_psi, num_amu_value, num_aki_value, num_neo_len, neo_file)
    if ( num_neo_len <= 2 ) then 
      write(*,*) '  ERROR: Could not read the numerical profile.'
      stop
    end if
    ! (no additional checks done at present)
  end if

  num_rot = ( rot_file /= 'none' )
  if ( (num_rot) .and. ( my_id == 0 ) ) then
    call readProf(num_rot_x, num_rot_y0, num_rot_len, rot_file)
    call check_num_prof(num_rot, num_rot_x, num_rot_y0, num_rot_len, 'rot', check_positive=.false.)
  end if
  
  
  
  contains
  
  
  
  subroutine check_num_prof(num, x, y, len, name, check_positive)
    
    implicit none
    
    logical,             intent(in) :: num
    real*8, allocatable, intent(in) :: x(:)
    real*8, allocatable, intent(in) :: y(:)
    integer,             intent(in) :: len
    character(len=*),    intent(in) :: name
    logical,             intent(in) :: check_positive !< should we check for positivity?
    
    integer :: i
    
    if ( num ) then
      if ( len < 2 ) then 
        write(*,*) '  ERROR: Could not read the numerical profile for '//trim(name)
        stop
      end if
      if ( abs(x(1))>1.d-6 ) then
        write(*,*) 'WARNING: Numerical '//trim(name)//' input does not start at Psi_N=0'
      end if
      if ( x(len) < 1.d0 ) then
        write(*,*) 'WARNING: Numerical '//trim(name)//' input does not include Psi_N=1'
      end if
      if ( len < 50 ) then
        write(*,*) 'WARNING: Numerical '//trim(name)//' input has a small number of points'
        write(*,*) '  This might lead to numerical problems or inaccuracy'
        write(*,*) '  Typically recommended number of points ~200'
      else if ( len > 500 ) then
        write(*,*) 'WARNING: Numerical '//trim(name)//' input has a very large number of points'
        write(*,*) '  This might slow down your simulation'
        write(*,*) '  Typically recommended number of points ~200'
      end if
      if ( xpoint .and. ( x(len) < 1.03d0 ) ) then
        write(*,*) 'WARNING: Numerical '//trim(name)//' input in SOL possibly not defined'
        write(*,*) '  In X-point simulations, numerical profiles should be specified such that'
        write(*,*) '  they cover the whole SOL region. Otherwise, extrapolation is used.'
      end if
      do i = 1, len-1
        if ( x(i+1)<=x(i) ) then
          write(*,*) 'ERROR: Numerical '//trim(name)//' input not correct'
          write(*,*) '  Psi_N values do not increase in a strictly monotonic way'
          stop
        end if
      end do
      if ( (minval(x) /= minval(x)) .or. (minval(y) /= minval(y)) ) then
        write(*,*) 'ERROR: Numerical '//trim(name)//' input contains NaNs'
        stop
      end if
      if ( check_positive ) then
        if ( minval(y) < 0.d0 ) then
          write(*,*) 'ERROR: Numerical '//trim(name)//' input has non-positive values'
          stop
        end if
      end if
    end if
    
  end subroutine check_num_prof

end subroutine read_num_profiles
