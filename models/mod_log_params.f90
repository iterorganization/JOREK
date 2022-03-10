module mod_log_params

implicit none

contains

!> Write out all relevant parameters defined in parameters
!! and by the namelist input file.
subroutine log_parameters(my_id, short)

use phys_module
use mumps_module,  only: no_zeros_mumps, mumps_ordering
use pastix_module, only: no_zeros_pastix, pastix_smp_only, pastix_pivot, pastix_maxthrd
use vacuum
use gauss, only: n_gauss

implicit none

! --- Routine parameters
integer,           intent(in) :: my_id !< MPI proc id
logical, optional             :: short !< commandline short version or run long version

! --- Constants
character(len=512), parameter :: REAL_FMT = "(1X,A, ' = ', 99ES12.4)"
character(len=512), parameter :: REAL_FMT2 = "(1X,A, ' = ', ES12.4, A)"
character(len=512), parameter :: INTG_FMT = "(1X,A, ' = ', 100I12)"
character(len=512), parameter :: INTG_FMT2 = "(1X,A, ' = ', I12, A)"
character(len=512), parameter :: LOGI_FMT = "(1X,A, ' = ', 10L12)"
character(len=512), parameter :: REA2_FMT = "(1X,A, ' = ', 4ES12.4, '     ...    ', 4ES12.4)"
character(len=512), parameter :: REA3_FMT = "(1X,A, ' = ', 9ES12.4, '     ...')"
character(len=512), parameter :: VARI_FMT = "(3x,I3,': ',A)"
character(len=512), parameter :: MODE_FMT = "(3x,I3,': ',A,'(',A,'*phi)')"
character(len=512), parameter :: CHAR_FMT = "(1X,A, ' = ""', A, '""')"
character(len=512), parameter :: CHAR_FMT2 = "(1X,A,I2,A,' = ""',A,'""')"

! --- Local variables
integer           :: ivar, itor
integer           :: i, j, n_rows !> do loop index 
character(len=10) :: mode_num
logical           :: short2

! --- Text out format
200 format(' ',79('*'))
112 format(A, i12)
  
if (present(short)) then
  short2 = short
else
  short2 = .false.
end if

if (my_id == 0) then

  write(*,*)
  write(*,200)
  write(*,*) '* Preprocessor Options                                                        *'
  write(*,200)

  write(*,'(1x,a)',advance='no') ' USE_MUMPS           : '
#ifdef USE_MUMPS
  write(*,*) 'on'
#else
  write(*,*) 'off'
#endif

  write(*,'(1x,a)',advance='no') ' USE_FFTW            : '
#ifdef USE_FFTW
  write(*,*) 'on'
#else
  write(*,*) 'off'
#endif

  write(*,'(1x,a)',advance='no') ' USE_PASTIX          : '
#ifdef USE_PASTIX
  write(*,*) 'on'
#else
  write(*,*) 'off'
#endif

  write(*,'(1x,a)',advance='no') ' USE_WSMP            : '
#ifdef USE_WSMP
  write(*,*) 'on'
#else
  write(*,*) 'off'
#endif

#ifdef USE_MURGE
  write(*,*) 'WARNING: USE_MURGE IS NOT SUPPORTED ANY MORE'
#endif

  write(*,'(1x,a)',advance='no') ' USE_HIPS            : '
#ifdef USE_HIPS
  write(*,*) 'on'
#else
  write(*,*) 'off'
#endif

  write(*,'(1x,a)',advance='no') ' USE_BLOCK           : '
#ifdef USE_BLOCK
  write(*,*) 'on'
#else
  write(*,*) 'off'
#endif

  write(*,'(1x,a)',advance='no') ' USE_HDF5            : '
#ifdef USE_HDF5
  write(*,*) 'on'
#else
  write(*,*) 'off'
#endif

  write(*,'(1x,a)',advance='no') ' MEMTRACE            : '
#ifdef MEMTRACE
  write(*,*) 'on'
#else
  write(*,*) 'off'
#endif

  write(*,'(1x,a)',advance='no') ' JECCD               : '
#ifdef JECCD
  write(*,*) 'on'
#else
  write(*,*) 'off'
#endif

  write(*,'(1x,a)',advance='no') ' COMPARE_ELEMENT_MATRIX : '
#ifdef COMPARE_ELEMENT_MATRIX
  write(*,*) 'on'
#else
  write(*,*) 'off'
#endif

  write(*,'(1x,a)',advance='no') ' CONSTRUCT_MATRIX_OMP_ATOMIC : '
#ifdef CONSTRUCT_MATRIX_OMP_ATOMIC
  write(*,*) 'on'
#else
  write(*,*) 'off'
#endif

  write(*,'(1x,a)',advance='no') ' PRINT_ELM_RHS : '
#ifdef PRINT_ELM_RHS
  write(*,*) 'on'
#else
  write(*,*) 'off'
#endif

  write(*,'(1x,a)',advance='no') ' DIRECT_CONSTRUCTION : '
#ifdef DIRECT_CONSTRUCTION
  write(*,*) 'on'
#else
  write(*,*) 'off'
#endif

write(*,'(1x,a)',advance='no') ' USE_COMPLEX_PRECOND          : '
#ifdef USE_COMPLEX_PRECOND
  write(*,*) 'on'
#else
  write(*,*) 'off'
#endif 

  write(*,'(1x,a)',advance='no') ' GAUSS_ORDER : '
#ifdef GAUSS_ORDER
  write(*,*) 'Preprocessor flag has been set! Thus, n_gauss=', n_gauss
#else
  write(*,*) 'Preprocessor flag not set. Thus, n_gauss=', n_gauss
#endif

write(*,'(1x,a)',advance='no') ' USE_BICGSTAB : '
#ifdef USE_BICGSTAB
  write(*,*) 'on'
#else
  write(*,*) 'off'
#endif

  write(*,*)
  write(*,200)
  write(*,*) '* Hard-Coded Parameters:                                                      *'
  write(*,200)
  write(*,  112) ' jorek_model    =  ', jorek_model       
  write(*,  112) ' n_var          =  ', n_var             
  write(*,112, advance='no') ' variable_names ='
  do i = 1, n_var
    write(*,'(" ",a)', advance='no') trim(variable_names(i))
  end do
  write(*,*)
  write(*,  112) ' n_dim          =  ', n_dim             
  write(*,  112) ' n_order        =  ', n_order           
  write(*,  112) ' n_tor          =  ', n_tor             
  write(*,  112) ' n_period       =  ', n_period          
  write(*,  112) ' n_plane        =  ', n_plane           
  write(*,  112) ' n_vertex_max   =  ', n_vertex_max      
  write(*,  112) ' n_nodes_max    =  ', n_nodes_max       
  write(*,  112) ' n_elements_max =  ', n_elements_max    
  write(*,  112) ' n_boundary_max =  ', n_boundary_max    
  write(*,  112) ' n_pieces_max   =  ', n_pieces_max      
  write(*,  112) ' n_degrees      =  ', n_degrees         
  write(*,  112) ' nref_max       =  ', nref_max          
  write(*,  112) ' n_ref_list     =  ', n_ref_list        

  ! stop function when case this log function is called from command line function
  if ( short2 ) return

  write(*,*)
  write(*,200)
  write(*,*) '* NAMELIST INPUT PARAMETERS                                                   *'
  write(*,200)
  write(*,CHAR_FMT) 'time_evol_scheme      ', trim(time_evol_scheme)
  write(*,INTG_FMT) 'n_tor_fft_thresh      ', n_tor_fft_thresh
  if ( n_tor .ge. n_tor_fft_thresh ) then
    write(*,'(3x,a)') '=> fft version of element_matrix is used'
  else
    write(*,'(3x,a)') '=> no-fft version of element_matrix is used'
  end if
  write(*,REAL_FMT) 'tstep                 ', tstep
  write(*,INTG_FMT) 'nstep                 ', nstep
  write(*,REAL_FMT) 'tstep_n               ', tstep_n
  write(*,INTG_FMT) 'nstep_n               ', nstep_n
  write(*,LOGI_FMT) 'eta_T_dependent       ', eta_T_dependent
  write(*,REAL_FMT) 'eta                   ', eta
  write(*,REAL_FMT) 'eta_ohmic             ', eta_ohmic
  write(*,REAL_FMT) 'eta_Spitzer (not input parameter; printed for reference in JOREK units)', eta_Spitzer
  write(*,REAL_FMT) 'T_max_eta             ', T_max_eta
  write(*,REAL_FMT) 'T_max_eta_ohm         ', T_max_eta_ohm  
  write(*,LOGI_FMT) 'visco_T_dependent     ', visco_T_dependent
  write(*,REAL_FMT) 'visco                 ', visco
  write(*,REAL_FMT) 'visco_par             ', visco_par
  write(*,LOGI_FMT) 'restart               ', restart
  write(*,INTG_FMT) 'rst_format            ', rst_format
  write(*,INTG_FMT) 'rst_hdf5              ', rst_hdf5
  write(*,INTG_FMT) 'rst_hdf5_version      ', rst_hdf5_version
  write(*,LOGI_FMT) 'regrid                ', regrid
  write(*,LOGI_FMT) 'write_ps              ', write_ps
  write(*,INTG_FMT) 'n_R                   ', n_R
  write(*,INTG_FMT) 'n_Z                   ', n_Z
  write(*,INTG_FMT) 'n_radial              ', n_radial
  write(*,INTG_FMT) 'n_pol                 ', n_pol

  if ( n_radial > 0 ) then
    write(*,REAL_FMT) 'psi_axis_init         ', psi_axis_init
    write(*,REAL_FMT) 'xr_r                  ', xr_r(:)
    write(*,REAL_FMT) 'sig_r                 ', sig_r(:)
    write(*,REAL_FMT) 'xr_tht                ', xr_tht(:)
    write(*,REAL_FMT) 'sig_tht               ', sig_tht(:)
    write(*,REAL_FMT) 'xr_z                  ', xr_z(:)
    write(*,REAL_FMT) 'sig_z                 ', sig_z(:)
    write(*,REAL_FMT) 'bgf_r                 ', bgf_r
    write(*,REAL_FMT) 'bgf_z                 ', bgf_z
  end if

  write(*,INTG_FMT) 'n_tht                 ', n_tht
  write(*,LOGI_FMT) 'n_tht_equidistant     ', n_tht_equidistant
  write(*,INTG_FMT) 'n_flux                ', n_flux
  write(*,LOGI_FMT) 'xpoint                ', xpoint

  if ( xpoint ) then
    write(*,INTG_FMT) 'xcase                 ', xcase
    write(*,INTG_FMT) 'n_open                ', n_open
    write(*,INTG_FMT) 'n_private             ', n_private
    write(*,INTG_FMT) 'n_leg                 ', n_leg
    write(*,INTG_FMT) 'n_leg_out             ', n_leg_out
    write(*,INTG_FMT) 'n_ext                 ', n_ext
    write(*,INTG_FMT) 'n_outer               ', n_outer
    write(*,INTG_FMT) 'n_inner               ', n_inner
    write(*,LOGI_FMT) 'force_horizontal_xline', force_horizontal_xline
    write(*,INTG_FMT) 'n_up_priv             ', n_up_priv
    write(*,INTG_FMT) 'n_up_leg              ', n_up_leg
    write(*,INTG_FMT) 'n_up_leg_out          ', n_up_leg_out
    write(*,REAL_FMT) 'SIG_closed            ', SIG_closed
    write(*,REAL_FMT) 'SIG_open              ', SIG_open
    write(*,REAL_FMT) 'SIG_private           ', SIG_private
    write(*,REAL_FMT) 'SIG_theta             ', SIG_theta
    write(*,REAL_FMT) 'SIG_leg_0             ', SIG_leg_0
    write(*,REAL_FMT) 'SIG_leg_1             ', SIG_leg_1
    write(*,REAL_FMT) 'SIG_outer             ', SIG_outer
    write(*,REAL_FMT) 'SIG_inner             ', SIG_inner
    write(*,REAL_FMT) 'SIG_up_leg_0          ', SIG_up_leg_0
    write(*,REAL_FMT) 'SIG_up_leg_1          ', SIG_up_leg_1
    write(*,REAL_FMT) 'SIG_up_priv           ', SIG_up_priv
    write(*,REAL_FMT) 'dPSI_open             ', dPSI_open
    write(*,REAL_FMT) 'dPSI_private          ', dPSI_private
    write(*,REAL_FMT) 'dPSI_outer            ', dPSI_outer
    write(*,REAL_FMT) 'dPSI_inner            ', dPSI_inner
    write(*,REAL_FMT) 'dPSI_up_priv          ', dPSI_up_priv
    write(*,INTG_FMT) 'first_target_point    ', first_target_point
    write(*,INTG_FMT) 'last_target_point     ', last_target_point
    write(*,REAL_FMT) 'SDN_threshold         ', SDN_threshold
  end if

  if ( ( (grid_to_wall) .or. (extend_existing_grid) ) .and. (n_wall_blocks .gt. 0) ) then
    write(*,REAL_FMT) 'eqdsk_psi_fact        ', eqdsk_psi_fact
    write(*,LOGI_FMT) 'RZ_grid_inside_wall   ', RZ_grid_inside_wall
    write(*,REAL_FMT) 'RZ_grid_jump_thres    ', RZ_grid_jump_thres
    write(*,INTG_FMT) 'n_wall_blocks         ', n_wall_blocks
    do i=1,n_wall_blocks
      write(*,INTG_FMT) 'Wall Patch number:    ', i
      write(*,LOGI_FMT) 'n_ext_equidistant:    ', n_ext_equidistant(i)
      write(*,INTG_FMT) 'corner block:         ', corner_block(i)
      write(*,INTG_FMT) 'resolution of block:  ', n_ext_block(i)
      write(*,INTG_FMT) 'n_block_points_left   ', n_block_points_left(i)
      do j=1,n_block_points_left(i)
        write(*,INTG_FMT) 'Patch left  point:    ', j
        write(*,REAL_FMT) 'R_block_points_left   ', R_block_points_left(i,j)
        write(*,REAL_FMT) 'Z_block_points_left   ', Z_block_points_left(i,j)
      enddo
      write(*,INTG_FMT) 'n_block_points_right  ', n_block_points_right(i)
      do j=1,n_block_points_right(i)
        write(*,INTG_FMT) 'Patch right point:    ', j
        write(*,REAL_FMT) 'R_block_points_right  ', R_block_points_right(i,j)
        write(*,REAL_FMT) 'Z_block_points_right  ', Z_block_points_right(i,j)
      enddo
    enddo
    write(*,LOGI_FMT) 'use_simple_bnd_types  ', use_simple_bnd_types
  endif

  write(*,INTG_FMT) 'nout                  ', nout
  write(*,REAL_FMT) 'xr1                   ', xr1
  write(*,REAL_FMT) 'sig1                  ', sig1
  write(*,REAL_FMT) 'xr2                   ', xr2
  write(*,REAL_FMT) 'sig2                  ', sig2
  write(*,REAL_FMT) 'R_begin               ', R_begin
  write(*,REAL_FMT) 'R_end                 ', R_end
  write(*,REAL_FMT) 'Z_begin               ', Z_begin
  write(*,REAL_FMT) 'Z_end                 ', Z_end
  write(*,REAL_FMT) 'rect_grid_vac_psi     ', rect_grid_vac_psi
  write(*,REAL_FMT) 'R_geo                 ', R_geo
  write(*,REAL_FMT) 'Z_geo                 ', Z_geo
  write(*,REAL_FMT) 'amin                  ', amin
  write(*,INTG_FMT) 'mf                    ', mf

  if ( mf >= 0 ) then
    write(*,REA3_FMT) 'fbnd                  ', fbnd(1:MIN(9,mf))
    write(*,REA3_FMT) 'fpsi                  ', fpsi(1:MIN(9,mf))
  end if

  write(*,REAL_FMT) 'F0                    ', F0
  write(*,REAL_FMT) 'zjz_0                 ', zjz_0
  write(*,REAL_FMT) 'zjz_1                 ', zjz_1
  write(*,REAL_FMT) 'zj_coef               ', zj_coef

  if ( .not. num_ffprime ) then
    write(*,REAL_FMT) 'FF_0                  ', FF_0
    write(*,REAL_FMT) 'FF_1                  ', FF_1
    write(*,REAL_FMT) 'FF_coef               ', FF_coef(1:8)
  else
    write(*,CHAR_FMT) 'ffprime_file          ', trim(ffprime_file)
  end if

  if ( .not. num_rho ) then
    write(*,REAL_FMT) 'rho_0                 ', rho_0
    write(*,REAL_FMT) 'rho_1                 ', rho_1
    write(*,REAL_FMT) 'rho_coef              ', rho_coef(1:5)
  else
    write(*,CHAR_FMT) 'rho_file              ', trim(rho_file)
    write(*,REAL_FMT) 'rho_0                 ', rho_0
    write(*,REAL_FMT) 'rho_1                 ', rho_1
  end if

  if (with_neutrals) then
    if ( .not. num_rhon ) then
      write(*,REAL_FMT) 'rhon_0                ', rhon_0
      write(*,REAL_FMT) 'rhon_1                ', rhon_1
      write(*,REAL_FMT) 'rhon_coef             ', rhon_coef(1:5)
    else
      write(*,CHAR_FMT) 'rhon_file             ', trim(rhon_file)
    end if
  end if

  if ( num_rot ) then
    write(*,CHAR_FMT) 'rot_file              ', trim(rot_file)
  else
    write(*,REAL_FMT) 'V_0                   ', V_0
    write(*,REAL_FMT) 'V_1                   ', V_1
    write(*,REAL_FMT) 'V_coeff               ', V_coef(1:10)
  end if

  if ( (abs(V_0) .ge. 1.d-19) .or. (num_rot) ) then
     write(*,LOGI_FMT) 'normalized_velocity_profile', normalized_velocity_profile
  endif


  if (with_TiTe) then ! (with_TiTe), i.e. single temperature ***************************************
    if ( .not. num_Te ) then
      write(*,REAL_FMT) 'Te_0                   ', Te_0
      write(*,REAL_FMT) 'Te_1                   ', Te_1
      write(*,REAL_FMT) 'Te_coef                ', Te_coef(1:5)
    else
      write(*,CHAR_FMT) 'Te_file                ', trim(Te_file)
      write(*,REAL_FMT) 'Te_0                   ', Te_0
      write(*,REAL_FMT) 'Te_1                   ', Te_1
    end if
    if ( .not. num_Ti ) then
      write(*,REAL_FMT) 'Ti_0                   ', Ti_0
      write(*,REAL_FMT) 'Ti_1                   ', Ti_1
      write(*,REAL_FMT) 'Ti_coef                ', Ti_coef(1:5)
    else
      write(*,CHAR_FMT) 'Ti_file                ', trim(Ti_file)
      write(*,REAL_FMT) 'Ti_0                   ', Ti_0
      write(*,REAL_FMT) 'Ti_1                   ', Ti_1
    end if
    if ( .not. num_zk_e_perp ) then
      write(*,REAL_FMT) 'ZK_e_perp             ', ZK_e_perp(1:6)
    else
      write(*,CHAR_FMT) 'ZK_e_perp_file        ', trim(ZK_e_perp_file)
    end if
    if ( .not. num_zk_i_perp ) then
      write(*,REAL_FMT) 'ZK_i_perp             ', ZK_i_perp(1:6)
    else
      write(*,CHAR_FMT) 'ZK_i_perp_file        ', trim(ZK_i_perp_file)
    end if
    write(*,REAL_FMT) 'heatsource_e           ', heatsource_e
    write(*,REAL_FMT) 'heatsource_e_psin      ', heatsource_e_psin
    write(*,REAL_FMT) 'heatsource_e_sig       ', heatsource_e_sig
    write(*,REAL_FMT) 'heatsource_gauss_e     ', heatsource_gauss_e
    write(*,REAL_FMT) 'heatsource_gauss_e_psin', heatsource_gauss_e_psin
    write(*,REAL_FMT) 'heatsource_gauss_e_sig ', heatsource_gauss_e_sig
    write(*,REAL_FMT) 'heatsource_i           ', heatsource_i
    write(*,REAL_FMT) 'heatsource_i_psin      ', heatsource_i_psin
    write(*,REAL_FMT) 'heatsource_i_sig       ', heatsource_i_sig
    write(*,REAL_FMT) 'heatsource_gauss_i     ', heatsource_gauss_i
    write(*,REAL_FMT) 'heatsource_gauss_i_psin', heatsource_gauss_i_psin
    write(*,REAL_FMT) 'heatsource_gauss_i_sig ', heatsource_gauss_i_sig
    write(*,REAL_FMT) 'ZK_e_par               ', ZK_e_par
    write(*,REAL_FMT) 'ZK_i_par               ', ZK_i_par
    write(*,REAL_FMT) 'ZK_par_max            ', ZK_par_max
    write(*,REAL_FMT) 'ZK_e_par_SpitzerHaerm (not input parameter; printed for reference in JOREK units)', ZK_e_par_SpitzerHaerm
    write(*,REAL_FMT) 'ZK_i_par_SpitzerHaerm (not input parameter; printed for reference in JOREK units)', ZK_i_par_SpitzerHaerm
    write(*,LOGI_FMT) 'ZKpar_T_dependent     ', ZKpar_T_dependent
    write(*,LOGI_FMT) 'thermalization         ', thermalization

  else ! (with_TiTe), i.e. single temperature ******************************************************

    if ( .not. num_T ) then
      write(*,REAL_FMT) 'T_0                   ', T_0
      write(*,REAL_FMT) 'T_1                   ', T_1
      write(*,REAL_FMT) 'T_coef                ', T_coef(1:5)
    else
      write(*,CHAR_FMT) 'T_file                ', trim(T_file)
      write(*,REAL_FMT) 'T_0                   ', T_0
      write(*,REAL_FMT) 'T_1                   ', T_1
    end if
    if ( .not. num_zk_perp ) then
      write(*,REAL_FMT) 'ZK_perp               ', ZK_perp(1:6)
    else
      write(*,CHAR_FMT) 'ZK_perp_file          ', trim(ZK_perp_file)
    end if
    write(*,REAL_FMT) 'ZK_par                ', ZK_par
    write(*,REAL_FMT) 'ZK_par_max            ', ZK_par_max
    write(*,REAL_FMT) 'ZK_par_SpitzerHaerm (not input parameter; printed for reference in JOREK units)', ZK_par_SpitzerHaerm
    write(*,LOGI_FMT) 'ZKpar_T_dependent     ', ZKpar_T_dependent
    write(*,REAL_FMT) 'heatsource            ', heatsource
    write(*,REAL_FMT) 'heatsource_psin       ', heatsource_psin
    write(*,REAL_FMT) 'heatsource_sig        ', heatsource_sig
    write(*,REAL_FMT) 'heatsource_gauss      ', heatsource_gauss
    write(*,REAL_FMT) 'heatsource_gauss_psin ', heatsource_gauss_psin
    write(*,REAL_FMT) 'heatsource_gauss_sig  ', heatsource_gauss_sig
  end if ! (with_TiTe), i.e. single temperature ****************************************************
  write(*,REAL_FMT) 'D_par                 ', D_par
  if ( .not. num_d_perp ) then
    write(*,REAL_FMT) 'D_perp                ', D_perp(1:6)
  else
    write(*,CHAR_FMT) 'D_perp_file           ', trim(D_perp_file)
  end if
#ifdef WITH_Impurities
  write(*,REAL_FMT) 'D_par_imp               ', D_par_imp
  if ( .not. num_d_perp_imp ) then
    write(*,REAL_FMT) 'D_perp_imp            ', D_perp_imp(1:6)
  else
    write(*,CHAR_FMT) 'D_perp_imp_file       ', trim(D_perp_imp_file)
  end if
#endif
  write(*,REAL_FMT) 'particlesource        ', particlesource
  write(*,REAL_FMT) 'particlesource_psin   ', particlesource_psin
  write(*,REAL_FMT) 'particlesource_sig    ', particlesource_sig
  write(*,REAL_FMT) 'edgeparticlesource    ', edgeparticlesource
  write(*,REAL_FMT) 'edgeparticlesource_psin', edgeparticlesource_psin
  write(*,REAL_FMT) 'edgeparticlesource_sig', edgeparticlesource_sig
  write(*,REAL_FMT) 'particlesource_gauss  ', particlesource_gauss
  write(*,REAL_FMT) 'particlesource_gauss_psin', particlesource_gauss_psin
  write(*,REAL_FMT) 'particlesource_gauss_sig ', particlesource_gauss_sig
  write(*,REAL_FMT) 'gamma                 ', gamma
  write(*,REAL_FMT) 'tauIC                 ', tauIC
  write(*,REAL_FMT) 'tauIC_nominal (not input parameter; printed for reference in JOREK units)', tauIC_nominal
  write(*,LOGI_FMT) 'Wdia                  ', Wdia
  write(*,REAL_FMT) 'eta_num               ', eta_num
  write(*,LOGI_FMT) 'eta_num_T_dependent   ', eta_num_T_dependent
  write(*,LOGI_FMT) 'eta_num_psin_dependent', eta_num_psin_dependent
  write(*,REAL_FMT) 'eta_num_prof          ', eta_num_prof(1:6)
  write(*,REAL_FMT) 'visco_num             ', visco_num
  write(*,LOGI_FMT) 'visco_num_T_dependent ', visco_num_T_dependent
  write(*,REAL_FMT) 'visco_par_num         ', visco_par_num
  write(*,REAL_FMT) 'D_perp_num            ', D_perp_num
  write(*,REAL_FMT) 'Dn_perp_num           ', Dn_perp_num

  write(*,LOGI_FMT) 'use_sc                ', use_sc
  write(*,REAL_FMT) 'visco_sc_num          ', visco_sc_num
  write(*,REAL_FMT) 'D_perp_sc_num         ', D_perp_sc_num
  write(*,REAL_FMT) 'D_par_sc_num          ', D_par_sc_num
  write(*,REAL_FMT) 'ZK_perp_sc_num        ', ZK_perp_sc_num
  write(*,REAL_FMT) 'ZK_par_sc_num         ', ZK_par_sc_num
  write(*,REAL_FMT) 'ZK_i_perp_sc_num      ', ZK_i_perp_sc_num
  write(*,REAL_FMT) 'ZK_i_par_sc_num       ', ZK_i_par_sc_num
  write(*,REAL_FMT) 'ZK_e_perp_sc_num      ', ZK_e_perp_sc_num
  write(*,REAL_FMT) 'ZK_e_par_sc_num       ', ZK_e_par_sc_num
  write(*,REAL_FMT) 'visco_par_sc_num      ', visco_par_sc_num
  write(*,REAL_FMT) 'Dn_pol_sc_num         ', Dn_pol_sc_num
  write(*,REAL_FMT) 'Dn_p_sc_num           ', Dn_p_sc_num

  if(jorek_model == 004 ) then
    write(*,REAL_FMT) 'HW_coef               ', HW_coef(1:2)
  endif

  write(*,LOGI_FMT) 'add_sources_in_sc     ', add_sources_in_sc
  if (with_TiTe) then
    write(*,REAL_FMT) 'ZK_i_perp_num         ', ZK_i_perp_num
    write(*,REAL_FMT) 'ZK_e_perp_num         ', ZK_e_perp_num
  else
    write(*,REAL_FMT) 'ZK_perp_num           ', ZK_perp_num
  end if
  write(*,REAL_FMT) 'tgnum                 ', tgnum(:)
  write(*,REAL_FMT) 'tgnum_psi             ', tgnum_psi 
  write(*,REAL_FMT) 'tgnum_u               ', tgnum_u   
  write(*,REAL_FMT) 'tgnum_zj              ', tgnum_zj  
  write(*,REAL_FMT) 'tgnum_w               ', tgnum_w   
  write(*,REAL_FMT) 'tgnum_rho             ', tgnum_rho 
  if (with_TiTe) then
    write(*,REAL_FMT) 'tgnum_Ti              ', tgnum_Ti  
    write(*,REAL_FMT) 'tgnum_Te              ', tgnum_Te  
  else
    write(*,REAL_FMT) 'tgnum_T               ', tgnum_T   
  end if
  write(*,REAL_FMT) 'tgnum_vpar            ', tgnum_vpar
  write(*,REAL_FMT) 'tgnum_rhon            ', tgnum_rhon
  write(*,REAL_FMT) 'tgnum_nre             ', tgnum_nre 
  write(*,REAL_FMT) 'tgnum_AR              ', tgnum_AR  
  write(*,REAL_FMT) 'tgnum_AZ              ', tgnum_AZ  
  write(*,REAL_FMT) 'tgnum_A3              ', tgnum_A3  



  write(*,LOGI_FMT) 'keep_current_prof     ', keep_current_prof
  write(*,LOGI_FMT) 'linear_run            ', linear_run
  write(*,REAL_FMT) 'D_prof_neg            ', D_prof_neg
  write(*,REAL_FMT) 'D_prof_neg_thresh     ', D_prof_neg_thresh
  if (with_TiTe) then
    write(*,REAL_FMT) 'ZK_e_prof_neg           ', ZK_e_prof_neg
    write(*,REAL_FMT) 'ZK_e_par_neg            ', ZK_e_par_neg
    write(*,REAL_FMT) 'ZK_e_prof_neg_thresh    ', ZK_e_prof_neg_thresh
    write(*,REAL_FMT) 'ZK_e_par_neg_thresh     ', ZK_e_par_neg_thresh
    write(*,REAL_FMT) 'ZK_i_prof_neg           ', ZK_i_prof_neg
    write(*,REAL_FMT) 'ZK_i_par_neg            ', ZK_i_par_neg
    write(*,REAL_FMT) 'ZK_i_prof_neg_thresh    ', ZK_i_prof_neg_thresh
    write(*,REAL_FMT) 'ZK_i_par_neg_thresh     ', ZK_i_par_neg_thresh
  else
    write(*,REAL_FMT) 'ZK_prof_neg           ', ZK_prof_neg
    write(*,REAL_FMT) 'ZK_par_neg            ', ZK_par_neg
    write(*,REAL_FMT) 'ZK_prof_neg_thresh    ', ZK_prof_neg_thresh
    write(*,REAL_FMT) 'ZK_par_neg_thresh     ', ZK_par_neg_thresh
  endif
  write(*,REAL_FMT) 'D_imp_extra_R         ', D_imp_extra_R
  write(*,REAL_FMT) 'D_imp_extra_Z         ', D_imp_extra_Z
  write(*,REAL_FMT) 'D_imp_extra_p         ', D_imp_extra_p
  write(*,REAL_FMT) 'D_imp_extra_neg       ', D_imp_extra_neg
  write(*,REAL_FMT) 'D_imp_extra_neg_thresh', D_imp_extra_neg_thresh
  write(*,REAL_FMT) 'T_min                 ', T_min
  write(*,REAL_FMT) 'ne_SI_min             ', ne_SI_min
  write(*,REAL_FMT) 'Te_eV_min             ', Te_eV_min
  write(*,REAL_FMT) 'rn0_min               ', rn0_min
  write(*,REAL_FMT) 'rho_min               ', rho_min
  write(*,LOGI_FMT) 'use_pellet            ', use_pellet
  write(*,REAL_FMT) 'corr_neg_temp_coef    ', corr_neg_temp_coef(:)
  write(*,REAL_FMT) 'corr_neg_dens_coef    ', corr_neg_dens_coef(:)

  if ( use_pellet ) then
    write(*,REAL_FMT) 'pellet_amplitude    ', pellet_amplitude
    write(*,REAL_FMT) 'pellet_R              ', pellet_R
    write(*,REAL_FMT) 'pellet_Z              ', pellet_Z
    write(*,REAL_FMT) 'pellet_phi            ', pellet_phi
    write(*,REAL_FMT) 'pellet_radius         ', pellet_radius
    write(*,REAL_FMT) 'pellet_sig            ', pellet_sig
    write(*,REAL_FMT) 'pellet_length         ', pellet_length
    write(*,REAL_FMT) 'pellet_psi            ', pellet_psi
    write(*,REAL_FMT) 'pellet_ellipse        ', pellet_ellipse
    write(*,REAL_FMT) 'pellet_theta          ', pellet_theta
    write(*,REAL_FMT) 'pellet_delta_psi      ', pellet_delta_psi
    write(*,REAL_FMT) 'pellet_density        ', pellet_density
    write(*,REAL_FMT) 'pellet_density_bg     ', pellet_density_bg
    write(*,REAL_FMT) 'pellet_particles      ', pellet_particles
    write(*,REAL_FMT) 'pellet_velocity_R     ', pellet_velocity_R
    write(*,REAL_FMT) 'pellet_velocity_Z     ', pellet_velocity_Z
  end if

  write(*,CHAR_FMT) 'tokamak_device        ', trim(tokamak_device)
  write(*,REAL_FMT) 'ellip                 ', ellip
  write(*,REAL_FMT) 'tria_u                ', tria_u
  write(*,REAL_FMT) 'tria_l                ', tria_l
  write(*,REAL_FMT) 'quad_u                ', quad_u
  write(*,REAL_FMT) 'quad_l                ', quad_l
  write(*,REAL_FMT) 'xampl                 ', xampl
  write(*,REAL_FMT) 'xwidth                ', xwidth
  write(*,REAL_FMT) 'xsig                  ', xsig
  write(*,REAL_FMT) 'xtheta                ', xtheta
  write(*,REAL_FMT) 'xshift                ', xshift
  write(*,REAL_FMT) 'xleft                 ', xleft
  write(*,CHAR_FMT) 'R_Z_psi_bnd_file      ', trim(R_Z_psi_bnd_file)
  write(*,CHAR_FMT) 'wall_file             ', trim(wall_file)
  write(*,INTG_FMT) 'n_boundary            ', n_boundary
  if ( n_boundary > 0 ) then
    write(*,REA2_FMT) 'r_boundary            ', r_boundary(1:4), r_boundary(n_boundary-3:n_boundary)
    write(*,REA2_FMT) 'z_boundary            ', z_boundary(1:4), z_boundary(n_boundary-3:n_boundary)
    write(*,REA2_FMT) 'psi_boundary          ', psi_boundary(1:4), psi_boundary(n_boundary-3:n_boundary)
  end if
  write(*,INTG_FMT) 'n_limiter             ', n_limiter
  if ( n_limiter > 0 ) then
    write(*,REA3_FMT) 'r_limiter             ', r_limiter(1:min(9,n_limiter))
    write(*,REA3_FMT) 'z_limiter             ', z_limiter(1:min(9,n_limiter))
  end if

  write(*,LOGI_FMT) 'freeboundary_equil    ', freeboundary_equil
  write(*,LOGI_FMT) 'freeboundary          ', freeboundary
  write(*,LOGI_FMT) 'freeb_change_indices  ', freeb_change_indices
  if ( freeboundary ) then
    write(*,LOGI_FMT) 'resistive_wall        ', resistive_wall
    if ( resistive_wall ) then
      write(*,REAL_FMT2) 'wall_resistivity      ', wall_resistivity, ' (used only if STARWALL response file_version==1)'
      write(*,REAL_FMT2) 'wall_resistivity_fact ', wall_resistivity_fact, ' (used only if STARWALL response file_version>=2)'
    end if
    write(*,REAL_FMT) 'start_VFB_ts          ', start_VFB_ts
    write(*,REAL_FMT) 'vert_FB_gain          ', vert_FB_gain(:)
    write(*,REAL_FMT) 'vert_FB_amp_ts        ', vert_FB_amp_ts(1:n_pf_coils)
    write(*,REAL_FMT) 'vert_FB_tact          ', vert_FB_tact
    write(*,CHAR_FMT) 'vert_pos_file         ', trim(vert_pos_file)
    write(*,REAL_FMT) 'I_coils_max           ', I_coils_max(1:n_pf_coils)

    
  end if
  
  if ( manipulate_psi_map(1,1) /= 0.d0 ) &
    write(*,REAL_FMT) 'manipulate_psi_map(1) ', manipulate_psi_map(1,:)
  if ( manipulate_psi_map(2,1) /= 0.d0 ) &
    write(*,REAL_FMT) 'manipulate_psi_map(2) ', manipulate_psi_map(2,:)
  if ( manipulate_psi_map(3,1) /= 0.d0 ) &
    write(*,REAL_FMT) 'manipulate_psi_map(3) ', manipulate_psi_map(3,:)
  if ( manipulate_psi_map(4,1) /= 0.d0 ) &
    write(*,REAL_FMT) 'manipulate_psi_map(4) ', manipulate_psi_map(4,:)
  if ( manipulate_psi_map(5,1) /= 0.d0 ) &
    write(*,REAL_FMT) 'manipulate_psi_map(5) ', manipulate_psi_map(5,:)
  
  write(*,REAL_FMT) 'amix                  ', amix
  write(*,REAL_FMT) 'equil_accuracy        ', equil_accuracy
  write(*,REAL_FMT) 'axis_srch_radius      ', axis_srch_radius
  write(*,REAL_FMT) 'delta_psi_GS          ', delta_psi_GS
  write(*,LOGI_FMT) 'newton_GS_fixbnd      ', newton_GS_fixbnd
  write(*,LOGI_FMT) 'newton_GS_freebnd     ', newton_GS_freebnd
  
  if (freeboundary_equil) then
    write(*,LOGI_FMT) 'starwall_equil_coils  ', starwall_equil_coils
    write(*,LOGI_FMT) 'find_pf_coil_currents ', find_pf_coil_currents
    write(*,LOGI_FMT) 'freeb_equil_iterate_area    ', freeb_equil_iterate_area
    write(*,REAL_FMT) 'amix_freeb            ', amix_freeb   
    write(*,REAL_FMT) 'equil_accuracy_freeb  ', equil_accuracy_freeb
    write(*,REAL_FMT) 'current_ref           ', current_ref
    write(*,REAL_FMT) 'cte_current_FB_fact   ', cte_current_FB_fact
    write(*,REAL_FMT) 'psi_offset_freeb      ', psi_offset_freeb
    write(*,REAL_FMT) 'FB_Ip_position        ', FB_Ip_position
    write(*,REAL_FMT) 'FB_Ip_integral        ', FB_Ip_integral
    write(*,REAL_FMT) 'Z_axis_ref            ', Z_axis_ref
    write(*,REAL_FMT) 'R_axis_ref            ', R_axis_ref
    write(*,REAL_FMT) 'FB_Zaxis_position     ', FB_Zaxis_position
    write(*,REAL_FMT) 'FB_Zaxis_derivative   ', FB_Zaxis_derivative
    write(*,REAL_FMT) 'FB_Zaxis_integral     ', FB_Zaxis_integral
    write(*,INTG_FMT) 'start_VFB             ', start_VFB
    write(*,INTG_FMT) 'n_feedback_current    ', n_feedback_current
    write(*,INTG_FMT) 'n_feedback_vertical   ', n_feedback_vertical
    write(*,INTG_FMT) 'n_iter_freeb          ', n_iter_freeb
    write(*,INTG_FMT) 'n_pf_coils            ', n_pf_coils
    write(*,REAL_FMT,advance='no') 'pf_coils%current      '
    do i = 1, n_pf_coils
      write(*,'(10ES12.4)',advance='no') pf_coils(i)%current
    end do
    write(*,*)
    if ( n_pf_coils > 0 ) &
        write(*,REAL_FMT) 'vert_FB_amp           ', vert_FB_amp(1:n_pf_coils)
    if ( n_pf_coils > 0 ) &
        write(*,REAL_FMT) 'rad_FB_amp           ', rad_FB_amp(1:n_pf_coils)
    write(*,REAL_FMT,advance='no') 'pf_coils%pert         '
    do i = 1, n_pf_coils
      write(*,'(10ES12.4)',advance='no') pf_coils(i)%pert
    end do
    write(*,*)
  endif


  write(*,REAL_FMT) 'Q_bar                 ', Q_bar
  write(*,REAL_FMT) 'sigma                 ', sigma
  write(*,REAL_FMT) 'density_reflection    ', density_reflection
  write(*,REAL_FMT) 'central_density       ', central_density
  write(*,REAL_FMT) 'central_mass          ', central_mass

  if (with_TiTe) then
    write(*,REAL_FMT) 'gamma_sheath_e        ', gamma_sheath_e
    write(*,REAL_FMT) 'gamma_sheath_i        ', gamma_sheath_i
    write(*,REAL_FMT) 'gamma_e_stangeby      ', gamma_e_stangeby
    write(*,REAL_FMT) 'gamma_i_stangeby      ', gamma_i_stangeby
  else 
    write(*,REAL_FMT) 'gamma_sheath          ', gamma_sheath
    write(*,REAL_FMT) 'gamma_stangeby        ', gamma_stangeby
  end if

  write(*,LOGI_FMT) 'vpar_smoothing        ', vpar_smoothing
  if ( vpar_smoothing ) then
    write(*,REAL_FMT) 'vpar_smoothing_coef   ', vpar_smoothing_coef(:)
  end if
  write(*,REAL_FMT) 'min_sheath_angle      ', min_sheath_angle     

  write(*,LOGI_FMT) 'bc_natural_open       ', bc_natural_open
  write(*,LOGI_FMT) 'produce_live_data     ', produce_live_data
  write(*,LOGI_FMT) 'export_for_nemec      ', export_for_nemec
  write(*,LOGI_FMT) 'keep_n0_const         ', keep_n0_const
  write(*,LOGI_FMT) 'gmres                 ', gmres
  write(*,INTG_FMT) 'gmres_max_iter        ', gmres_max_iter
  write(*,REAL_FMT) 'gmres tolerance       ', gmres_tol
  write(*,INTG_FMT) 'iter_precon           ', iter_precon
  write(*,INTG_FMT) 'max_steps_noUpdate    ', max_steps_noUpdate
  write(*,INTG_FMT) 'gmres_m               ', gmres_m
  write(*,REAL_FMT) 'gmres_4               ', gmres_4
  write(*,LOGI_FMT) 'centralize_harm_mat   ', centralize_harm_mat
  write(*,LOGI_FMT) 'use_mumps             ', use_mumps
  write(*,LOGI_FMT) 'use_wsmp              ', use_wsmp
  write(*,LOGI_FMT) 'use_pastix            ', use_pastix
  write(*,LOGI_FMT) 'use_strumpack         ', use_strumpack  
  write(*,LOGI_FMT) 'use_mumps_eq          ', use_mumps_eq
  write(*,LOGI_FMT) 'use_pastix_eq         ', use_pastix_eq
  write(*,LOGI_FMT) 'use_strumpack_eq      ', use_strumpack_eq  
  write(*,LOGI_FMT) 'pastix_smp_only       ', pastix_smp_only
  write(*,REAL_FMT) 'pastix_pivot          ', pastix_pivot
  write(*,INTG_FMT) 'pastix_maxthrd        ', pastix_maxthrd
  write(*,LOGI_FMT) 'refinement            ', refinement
  write(*,LOGI_FMT) 'force_central_node    ', force_central_node
  write(*,LOGI_FMT) 'fix_axis_nodes        ', fix_axis_nodes
  write(*,LOGI_FMT) 'grid_to_wall          ', grid_to_wall
  write(*,LOGI_FMT) 'adaptive_time         ', adaptive_time
  write(*,LOGI_FMT) 'equil                 ', equil
  write(*,LOGI_FMT) 'bench_without_plot    ', bench_without_plot
  write(*,LOGI_FMT) 'no_zeros_mumps        ', no_zeros_mumps
  write(*,LOGI_FMT) 'no_zeros_pastix       ', no_zeros_pastix
  write(*,LOGI_FMT) 'mach_one_bnd_integral ', mach_one_bnd_integral
  write(*,LOGI_FMT) 'deuterium_adas        ', deuterium_adas       
  write(*,LOGI_FMT) 'deuterium_adas_1e20   ', deuterium_adas_1e20
  write(*,LOGI_FMT) 'old_deuterium_atomic  ', old_deuterium_atomic
  write(*,LOGI_FMT) 'no_mach1_bc           ', no_mach1_bc

#ifdef fullmhd
    write(*,LOGI_FMT) 'Mach1_openBC          ', Mach1_openBC
    write(*,REA3_FMT) 'eta_ARAZ_const        ', eta_ARAZ_const
    write(*,LOGI_FMT) 'eta_ARAZ_on           ', eta_ARAZ_on
    write(*,LOGI_FMT) 'eta_ARAZ_simple       ', eta_ARAZ_simple
    write(*,LOGI_FMT) 'tauIC_ARAZ_on         ', tauIC_ARAZ_on
#endif

  if (use_mumps) then
    write(*,INTG_FMT) 'mumps_ordering        ', mumps_ordering
  endif

  write(*,LOGI_FMT) 'use_BLR_compression   ', use_BLR_compression
  if (use_BLR_compression) then
    write(*,REAL_FMT) 'epsilon_BLR           ', epsilon_BLR
    write(*,LOGI_FMT) 'just_in_time_BLR      ', just_in_time_BLR
    write(*,LOGI_FMT) 'pastix_blr_abs_tol    ', pastix_blr_abs_tol
  endif

  write(*,INTG_FMT) 'n_pfc                 ', n_pfc
  if ( n_pfc > 0 ) then
    write(*,REA3_FMT) 'Rmin_pfc              ', Rmin_pfc(1:min(9,n_pfc))
    write(*,REA3_FMT) 'Rmax_pfc              ', Rmax_pfc(1:min(9,n_pfc))
    write(*,REA3_FMT) 'Zmin_pfc              ', Zmin_pfc(1:min(9,n_pfc))
    write(*,REA3_FMT) 'Zmax_pfc              ', Zmax_pfc(1:min(9,n_pfc))
    write(*,REA3_FMT) 'current_pfc           ', current_pfc(1:min(9,n_pfc))
  end if

  write(*,INTG_FMT) 'n_jropes              ', n_jropes
  if ( n_jropes > 0 ) then
    write(*,REA3_FMT) 'R_jropes              ', R_jropes(1:n_jropes)
    write(*,REA3_FMT) 'Z_jropes              ', Z_jropes(1:n_jropes)
    write(*,REA3_FMT) 'w_jropes              ', w_jropes(1:n_jropes)
    write(*,REA3_FMT) 'current_jropes        ', current_jropes(1:n_jropes)
    write(*,REA3_FMT) 'rho_jropes            ', rho_jropes(1:n_jropes)
    write(*,REA3_FMT) 'T_jropes              ', T_jropes(1:n_jropes)
  end if

  write(*,LOGI_FMT) 'RMP_on                ', RMP_on
  if (RMP_on) then
     write(*,CHAR_FMT) 'RMP_psi_cos_file      ', trim(RMP_psi_cos_file)
     write(*,CHAR_FMT) 'RMP_psi_sin_file      ', trim(RMP_psi_sin_file)
     write(*,REAL_FMT) 'RMP_growth_rate       ', RMP_growth_rate
     write(*,REAL_FMT) 'RMP_ramp_up_time      ', RMP_ramp_up_time
     write(*,INTG_FMT) 'Number_RMP_harmonics  ', Number_RMP_harmonics 
     write(*,INTG_FMT) 'RMP_har_cos_spectrum  ', RMP_har_cos_spectrum(1:Number_RMP_harmonics)
     write(*,INTG_FMT) 'RMP_har_sin_spectrum  ', RMP_har_sin_spectrum(1:Number_RMP_harmonics)
  endif
  write(*,LOGI_FMT) 'output_bnd_elements   ', output_bnd_elements
  write(*,LOGI_FMT) 'bootstrap             ', bootstrap
  write(*,LOGI_FMT) 'NEO                   ', NEO
  if (NEO) then
    write(*,LOGI_FMT) 'num_neo_file          ', num_neo_file
    if (num_neo_file) then
      write(*,CHAR_FMT) 'neo_file              ', trim(neo_file)
    else
      write(*,REAL_FMT) 'amu_neo_const         ', amu_neo_const
      write(*,REAL_FMT) 'aki_neo_const         ', aki_neo_const
    endif
  endif

  if(jorek_model == 306 ) then
     write(*,REAL_FMT) 'nu_jec1_fast        ',  nu_jec1_fast
     write(*,REAL_FMT) 'nu_jec2_fast        ',  nu_jec2_fast
     write(*,REAL_FMT) 'JJ_par              ',  JJ_par
     write(*,REAL_FMT) 'jec_pos1            ',  jec_pos1
     write(*,REAL_FMT) 'jec_width           ',  jec_width
     write(*,REAL_FMT) 'jecamp              ',  jecamp
  endif

#if (defined WITH_Neutrals) || (defined WITH_Impurities)
     write(*,REAL_FMT) 'ns_amplitude        ',  ns_amplitude
     write(*,REAL_FMT) 'ns_R                ',  ns_R
     write(*,REAL_FMT) 'ns_Z                ',  ns_Z
     write(*,REAL_FMT) 'ns_phi              ',  ns_phi
     write(*,REAL_FMT) 'ns_radius           ',  ns_radius
     write(*,REAL_FMT) 'ns_radius_min       ',  ns_radius_min
     write(*,REAL_FMT) 'ns_deltaphi         ',  ns_deltaphi
     write(*,REAL_FMT) 'ns_delta_minor_rad  ',  ns_delta_minor_rad
     write(*,REAL_FMT) 'ns_tor_norm         ',  ns_tor_norm
     write(*,REAL_FMT) 'ksi_ion             ',  ksi_ion
     write(*,LOGI_FMT) 'JET_MGI             ',  JET_MGI
     write(*,LOGI_FMT) 'ASDEX_MGI           ',  ASDEX_MGI
     write(*,REAL_FMT) 'A_Dmv               ',  A_Dmv
     write(*,REAL_FMT) 'K_Dmv               ',  K_Dmv
     write(*,REAL_FMT) 'V_Dmv               ',  V_Dmv
     write(*,REAL_FMT) 'L_tube              ',  L_tube
     write(*,REAL_FMT) 't_ns                ',  t_ns
     write(*,REAL_FMT) 'delta_n_convection  ',  delta_n_convection
     write(*,REAL_FMT) 'nimp_bg             ',  nimp_bg
     write(*,INTG_FMT) 'n_adas              ',  n_adas
     do i = 1, n_adas
       write(*,CHAR_FMT2) 'imp_type(',i,')    ', trim(imp_type(i))
     end do
     write(*,INTG_FMT) 'index_main_imp      ', index_main_imp
     write(*,REAL_FMT) 'neutral_line_source ', neutral_line_source
     write(*,REAL_FMT) 'neutral_line_R_start', neutral_line_R_start
     write(*,REAL_FMT) 'neutral_line_Z_start', neutral_line_Z_start
     write(*,REAL_FMT) 'neutral_line_R_end  ', neutral_line_R_end
     write(*,REAL_FMT) 'neutral_line_Z_end  ', neutral_line_Z_end
     write(*,REAL_FMT) 'neutral_reflection  ', neutral_reflection
     write(*,REAL_FMT) 'imp_reflection      ', imp_reflection
     write(*,LOGI_FMT) 'output_prad_phi     ', output_prad_phi
     write(*,CHAR_FMT) 'adas_dir            ', trim(adas_dir)
     write(*,LOGI_FMT) 'use_imp_adas        ', use_imp_adas
     write(*,REAL_FMT) 'drift_distance      ', drift_distance
     write(*,REAL_FMT) 'energy_teleported   ', energy_teleported

     !< Additional log for SPI model
   if(using_spi) then
     write(*,LOGI_FMT) 'using_spi           ',  using_spi
     write(*,LOGI_FMT) 'spi_tor_rot         ',  spi_tor_rot
     write(*,LOGI_FMT) 'spi_num_vol         ',  spi_num_vol
     write(*,CHAR_FMT) 'adas_dir            ',  trim(adas_dir)
     write(*,INTG_FMT) 'n_spi               ',  n_spi
     write(*,INTG_FMT) 'n_spi_tot           ',  n_spi_tot
     write(*,INTG_FMT) 'n_inj               ',  n_inj
     do i = 1,n_inj
       write(*,CHAR_FMT2) 'spi_plume_file(',i,')    ',  trim(spi_plume_file(i))
     end do
     write(*,LOGI_FMT) 'spi_plume_hdf5      ',  spi_plume_hdf5
     write(*,INTG_FMT) 'spi_rnd_seed        ',  spi_rnd_seed
     write(*,INTG_FMT) 'spi_abl_model       ',  spi_abl_model
     do i = 1,n_inj
       write(*,CHAR_FMT2) 'spi_shard_file(',i,')    ',  trim(spi_shard_file(i))
     end do
     write(*,REAL_FMT) 'spi_Vel_Rref        ',  spi_Vel_Rref
     write(*,REAL_FMT) 'spi_Vel_Zref        ',  spi_Vel_Zref
     write(*,REAL_FMT) 'spi_Vel_RxZref      ',  spi_Vel_RxZref
     write(*,REAL_FMT) 'spi_quantity        ',  spi_quantity
     write(*,REAL_FMT) 'spi_quantity_bg     ',  spi_quantity_bg
     write(*,REAL_FMT) 'pellet_density      ',  pellet_density
     write(*,REAL_FMT) 'pellet_density_bg   ',  pellet_density_bg
     write(*,REAL_FMT) 'spi_angle           ',  spi_angle
     write(*,REAL_FMT) 'spi_L_inj           ',  spi_L_inj
     write(*,REAL_FMT) 'spi_L_inj_diff      ',  spi_L_inj_diff
     write(*,REAL_FMT) 'tor_frequency       ',  tor_frequency
     write(*,REAL_FMT) 'spi_Vel_diff        ',  spi_Vel_diff
   end if
#endif

  write(*,LOGI_FMT) 'restart_particles   ',restart_particles
  write(*,REAL_FMT) 'n_particles         ',n_particles
  write(*,INTG_FMT) 'nstep_particles     ',nstep_particles
  write(*,INTG_FMT) 'nsubstep_particles  ',nsubstep_particles
  write(*,REAL_FMT) 'tstep_particles     ',tstep_particles
  write(*,REAL_FMT) 'filter_perp,        ',filter_perp
  write(*,REAL_FMT) 'filter_hyper,       ',filter_hyper
  write(*,REAL_FMT) 'filter_par,         ',filter_par
  write(*,REAL_FMT) 'filter_perp_n0,     ',filter_perp_n0
  write(*,REAL_FMT) 'filter_hyper_n0,    ',filter_hyper_n0   
  write(*,REAL_FMT) 'filter_par_n0,      ',filter_par_n0     
  write(*,LOGI_FMT) 'use_ncs,            ',use_ncs     
  write(*,LOGI_FMT) 'use_ccs,            ',use_ccs    
  write(*,LOGI_FMT) 'use_pcs,            ',use_pcs
  write(*,LOGI_FMT) 'use_ionisation,     ',use_ionisation    
  write(*,LOGI_FMT) 'use_sputtering,     ',use_sputtering    
  write(*,LOGI_FMT) 'use_cx,             ',use_cx            

  write(*,*)
  write(*,200)
  write(*,*) '* NORMALIZATION FACTORS                                                       *'
  write(*,200)
  write(*,REAL_FMT) 'sqrt(mu0*rho0)      ',  sqrt_mu0_rho0 
  write(*,REAL_FMT) 'sqrt(mu0/rho0)      ',  sqrt_mu0_over_rho0 
  write(*,*)

end if

end subroutine log_parameters

end module mod_log_params
