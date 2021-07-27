subroutine equilibrium(my_id,node_list,element_list,bnd_node_list,bnd_elm_list,xpoint2,xcase2, nice_q)
!-----------------------------------------------------------------------
! Solve the Grad-Shafranov equation to determine the plasma equilibrium
!   both freeboundary and fixed boundary solutions
!-----------------------------------------------------------------------
use tr_module 
use mod_parameters
use data_structure
use phys_module
use mod_poiss
use mod_iterate2area
use mod_plasma_response
use equil_info
use vacuum
use mpi_mod
use mod_interp, only: interp
use mod_F_profile
implicit none

          
! --- Routine parameters
integer,                      intent(in)    :: my_id
type (type_node_list),        intent(inout) :: node_list
type (type_element_list),     intent(inout) :: element_list
type (type_bnd_node_list),    intent(inout) :: bnd_node_list
type (type_bnd_element_list), intent(inout) :: bnd_elm_list
logical,                      intent(in)    :: xpoint2
integer,                      intent(in)    :: xcase2
logical,                      intent(in)    :: nice_q

! --- Local variables.
type (type_surface_list) :: surface_list, sep_list
integer    :: ierr, n_iter, iter, i, in, mm, i_elm_axis, i_elm_xpoint(2), i_elm_lim, ifail, i_elm
real*8     :: amplitude, psi, psi_bnd
real*8     :: zn,  dn_dpsi,  dn_dpsi2,  dn_dz,  dn_dz2,  dn_dpsi_dz,  dn_dpsi3,  dn_dpsi2_dz,  dn_dpsi_dz2
real*8     :: zT,  dT_dpsi,  dT_dpsi2,  dT_dz,  dT_dz2,  dT_dpsi_dz,  dT_dpsi3,  dT_dpsi2_dz,  dT_dpsi_dz2
real*8     :: zTi, dTi_dpsi, dTi_dpsi2, dTi_dz, dTi_dz2, dTi_dpsi_dz, dTi_dpsi3, dTi_dpsi2_dz, dTi_dpsi_dz2
real*8     :: zTe, dTe_dpsi, dTe_dpsi2, dTe_dz, dTe_dz2, dTe_dpsi_dz, dTe_dpsi3, dTe_dpsi2_dz, dTe_dpsi_dz2
real*8     :: Ti_prof, Te_prof
real*8     :: zFFprime,dFFprime_dpsi,dFFprime_dz, dFFprime_dpsi_dz, dFFprime_dz2, dFFprime_dpsi2
real*8     :: F_prof, dF_dpsi, dF_dz, dF_dpsi2, dF_dz2, dF_dpsi_dz
real*8     :: xx, x_s, x_t, x_st, x_ss, x_tt, yy, y_s, y_t, y_st, y_ss, y_tt
real*8     :: R_axis, Z_axis, s_axis, t_axis, psi_axis,R, Z, BigR, T0, BigR_s, T0_s
real*8     :: R_lim, Z_lim, s_lim, t_lim, psi_lim, R_out, Z_out, s_out, t_out
real*8     :: R_xpoint(2),Z_xpoint(2),s_xpoint(2),t_xpoint(2), psi_xpoint(2)
real*8     :: zjz, dj_dpsi, dj_dR, dj_dZ, dj_dR_dZ, dj_dR_DR, dj_dZ_dZ, dj_dpsi2, dj_dR_dpsi, dj_dZ_dpsi, psi_n
real*8     :: ps0_s, ps0_t, p_s, p_t, p_ss, p_st, p_tt 
real*8     :: zj0_s, zj0_t, equil_error, equil_value, ps0_x, ps0_y, Z_s, Z_t, xjac, direction, Btot
real*8     :: current_tot, current_int, diff, R_xpoint2(2), Z_xpoint2(2)
real*8     :: sigmas(16), dZ_axis, Z_axis_int, Z_axis_old, area_ref
integer    :: n_grids(12)
logical    :: freeboundary_equil2
real*8     :: T_prof, T_0_old, FF_0_old, T_1_old, FF_1_old
real*8, allocatable     :: T_profile(:)
real*8     :: density_prof
real*8, allocatable     :: density_profile(:)

if (my_id .eq. 0) then
  write(*,*) '***************************************'
  write(*,*) '*           equilibrium               *'
  write(*,*) '***************************************'
  write(*,*) '   freeboundary_equil : ',freeboundary_equil
  write(*,*) '   X-point      : ',xpoint2
  write(*,*) '   Xcase        : ',xcase2

  if ((newton_GS_fixbnd .or. newton_GS_freebnd) .and. (use_pastix_eq)) then
#ifndef USE_PASTIX6
    write(*,*) ' '
    write(*,*) ' WARNING: PASTIX 5 IS NOT EFFICIENT FOR THE GRAD-SHAFRANOV SOLVER'
    write(*,*) '           WITH THE NEWTON METHOD. PLEASE USE PASTIX 6,          '
    write(*,*) '           MUMPS OR STRUMPACK INSTEAD. For example               '
    write(*,*) '           (add use_mumps_eq=.t. to namelist and  USE_MUMPS = 1  '
    write(*,*) '           in Makefile.inc)                                      '
    write(*,*) ' '
#endif
  endif

  if ((newton_GS_fixbnd .or. newton_GS_freebnd) .and. (.not. xpoint2)) then
    write(*,*) ' '
    write(*,*) ' WARNING: THE NEWTON METHOD FOR THE GRAD-SHAFRANOV SOLVER DOES   '
    write(*,*) '           NOT TAKE EFFECT FOR LIMITER PLASMAS (XPOINT=.F.)      '
    write(*,*) '           AND PICARD ITERATIONS ARE RECOVERED                   '
    write(*,*) '           FURTHER DEVELOPMENTS ARE NEEDED FOR LIMITER PLASMAS   '
    write(*,*) ' '
  endif

endif

freeboundary_equil2 = freeboundary_equil
freeboundary_equil  = .false.

!------------------------------------ fixed boundary equilibrium
n_iter       = 200
psi_bnd      = 0.d0
ES%psi_bnd   = 0.d0
Z_xpoint(1)  = -99.d0
Z_xpoint(2)  = +99.d0
R_xpoint(:)  = R_geo
vertical_FB  = 0.d0
i_elm_xpoint = 0 
current_tot  = 0.

if (my_id == 0) then

  do iter = 1, n_iter
  
  
    call update_equil_state(my_id,node_list, element_list, bnd_elm_list, xpoint, xcase)
    call print_equil_state(.true.)

    if ((ES%ifail_axis .ne. 0) .and. (iter .le. 5)) then
      call find_RZ(node_list,element_list,R_geo,Z_geo,R_out,Z_out,i_elm,s_out,t_out,ifail)
      call interp(node_list,element_list,i_elm,1,1,s_out,t_out,psi_axis,P_s,P_t,P_st,P_ss,P_tt)
      write(*,'(A,3f10.5)')  ' changed magnetic axis to :  ', R_out,Z_out,psi_axis
      ES%R_axis     = R_out;    ES%Z_axis = Z_out;  ES%psi_axis   = psi_axis;   
      ES%s_axis     = s_out;    ES%t_axis = t_out;  ES%i_elm_axis = i_elm;
      ES%ifail_axis = ifail   
    endif
    
    if (xpoint2) then
      if (ES%ifail_xpoint == 0) then ! (otherwise, keep the values of the previous iteration as a reasonable guess)
        ES%psi_bnd  = ES%psi_xpoint(1)
        if( (xcase2 .eq. UPPER_XPOINT) .or. ((xcase2 .eq. DOUBLE_NULL) .and. (abs(ES%psi_xpoint(2)-ES%psi_axis) .lt. abs(ES%psi_xpoint(1)-ES%psi_axis))) ) then
          ES%psi_bnd = ES%psi_xpoint(2)
        endif
        psi_bnd     = ES%psi_bnd
        R_xpoint(1) = ES%R_xpoint(1)
        Z_xpoint(1) = ES%Z_xpoint(1)
        R_xpoint(2) = ES%R_xpoint(2)
        Z_xpoint(2) = ES%Z_xpoint(2)
        if(xcase2 .eq. LOWER_XPOINT) ES%Z_xpoint(2) = +99.d0
        if(xcase2 .eq. UPPER_XPOINT) ES%Z_xpoint(1) = -99.d0
      else
        ES%R_xpoint = R_xpoint
        ES%Z_xpoint = Z_xpoint
        ES%psi_bnd  = psi_bnd
        if (freeboundary_equil) then
          ES%Z_xpoint(1) = -99.d0
          ES%Z_xpoint(2) = +99.d0
        endif
      endif
    else
      ES%psi_bnd = 0.d0
    endif

    if (.not. xpoint) then
      if ( (ES%Z_lim .gt. ES%Z_xpoint(1)) .and. (ES%Z_lim .lt. ES%Z_xpoint(2)) ) then
        if (n_limiter /= 0) then   ! else n_limiter = 0 and psi_bnd is set to 0
          ES%psi_bnd = ES%psi_lim
          write(*,'(A,3f8.3)') ' LIMITER PLASMA ',psi_lim,R_lim,Z_lim
        endif
      endif
    endif
  
    if(xcase2 .eq. LOWER_XPOINT) write(*,'(A,3es14.6,i3)') ' PSI_AXIS, PSI_BND  : ',ES%psi_axis,ES%psi_bnd,ES%Z_xpoint(1),ES%ifail_xpoint
    if(xcase2 .eq. UPPER_XPOINT) write(*,'(A,3es14.6,i3)') ' PSI_AXIS, PSI_BND  : ',ES%psi_axis,ES%psi_bnd,ES%Z_xpoint(2),ES%ifail_xpoint

    write(*,'(A,1f14.8)')                       ' PSI_BND - PSI_AXIS : ', ES%psi_bnd-ES%psi_axis 
  
    call poisson(my_id,-1,node_list,element_list,bnd_node_list,bnd_elm_list,3,1,1, &
                 ES%psi_axis,ES%psi_bnd,xpoint2,xcase2,ES%Z_xpoint,freeboundary_equil,refinement,iter)   !----------- for GS use -1
  
    diff = 0.d0
    do i=1, node_list%n_nodes
      diff = diff + abs(node_list%node(i)%deltas(1,1,1))
    enddo  
    diff = diff / float(node_list%n_nodes)
    
    write(*,'(A,I4,A,ES10.3)') ' Iteration ', iter, ': diff=', diff
    
    if ( (iter > 1) .and. (diff < equil_accuracy) ) then
      write(*,'(A,I4,A)') ' Fixed boundary equilibrium converged: after', iter, ' iterations'
      exit
    else if ( iter == n_iter) then
      write(*,'(A,ES10.3)') ' WARNING: Fixed boundary equilibrium not fully converged: diff=', diff
      exit
    end if
  
  enddo

end if ! my_id == 0

!--------------------------------------- freeboundary equilibrium
freeboundary_equil = freeboundary_equil2

current_int = 0.d0; Z_axis_int = 0.d0
 
T_0_old = T_0;  FF_0_old = FF_0;  T_1_old = T_1;  FF_1_old = FF_1

if (freeboundary_equil) then

  if (my_id == 0) then

    write(*,*)
    write(*,*) '------------------------------------------------------'
    write(*,*) '--- Iterative solution of freeboundary equilibrium ---'
    write(*,*) '------------------------------------------------------'
    write(*,*)

    ! Take target delta_psi from fixed boundary if not specified
    if ((delta_psi_GS >= 10000.d0) .and. newton_GS_freebnd) then
      write(*,*) ' '
      write(*,*) ' Taking target delta_psi_GS=psi_bnd-psi_axis from fixed boundary equilibrium'
      write(*,*) ' as it has not been specified in the input file'
      write(*,*) ' '
      delta_psi_GS = ES%psi_bnd - ES%psi_axis
    endif
    
    ! Target current and axis for Picard iterations
    if (current_ref .gt. 1.d20) then    !choose fix bnd equilibrium final current in case of non specification of target current
      call integral_current(node_list,element_list,ES%psi_axis,ES%psi_bnd, xpoint2, xcase2, ES%Z_xpoint, current_ref)
    endif
   
    if (Z_axis_ref .gt. 1.d20) then     !choose fix bnd equilibrium final Zaxis in case of non specification of target Zaxis
      Z_axis_ref = ES%Z_axis
    endif
    
    ! Target poloidal cross section area for limiter plasmas
    if (freeb_equil_iterate_area .and. (.not. xpoint2)) then
      n_limiter = 0    ! Use the full domain to search psibnd enclosing given area
      call area_inside_flux_contour(node_list,element_list, xpoint2, xcase2, ES%psi_bnd, area_ref, ES%R_lim, ES%Z_lim)
      write(*,*) ' The reference area from fixed boundaray is = ', area_ref
    endif
  
  end if ! my_id == 0

  do iter=1, n_iter_freeb

    if (my_id == 0) then
      
      write(*,*)
      write(*,'(1x,a,i5,a)') '>>> ITERATION', iter, ' <<<'
 
      call update_equil_state(my_id,node_list, element_list, bnd_elm_list, xpoint, xcase)
      call print_equil_state(.true.)

      if ((ES%ifail_axis .ne. 0) .and. (iter .le. 5)) then
        call find_RZ(node_list,element_list,R_geo,Z_geo,R_out,Z_out,i_elm,s_out,t_out,ifail)
        call interp(node_list,element_list,i_elm,1,1,s_out,t_out,psi_axis,P_s,P_t,P_st,P_ss,P_tt)
        write(*,'(A,3f10.5)')  ' changed magnetic axis to :  ', R_out,Z_out,psi_axis
        ES%R_axis     = R_out;    ES%Z_axis = Z_out;  ES%psi_axis   = psi_axis;   
        ES%s_axis     = s_out;    ES%t_axis = t_out;  ES%i_elm_axis = i_elm;
        ES%ifail_axis = ifail   
      endif
      
      write(10,'(i6,9e20.12)') iter, current_tot, R_axis, Z_axis, ES%psi_bnd-ES%psi_axis
      
      ES%psi_bnd = 0.d0
   
      if (xpoint2) then
        if (ES%ifail_xpoint .ne. 1) then      
          ES%psi_bnd  = ES%psi_xpoint(1)
          if( (xcase2 .eq. 2) .or. ((xcase2 .eq. 3) .and. (abs(ES%psi_xpoint(2)-ES%psi_axis) .lt. abs(ES%psi_xpoint(1)-ES%psi_axis))) ) then
            ES%psi_bnd = ES%psi_xpoint(2)
          endif
          if(xcase2 .eq. LOWER_XPOINT) ES%Z_xpoint(2) = +99.d0
          if(xcase2 .eq. UPPER_XPOINT) ES%Z_xpoint(1) = -99.d0
        else
          ES%Z_xpoint(1) = -99.d0 
          ES%Z_xpoint(2) = +99.d0
        endif
      endif
  
      if (.not. xpoint2) then
        if ( (ES%Z_lim .gt. ES%Z_xpoint(1)) .and. (ES%Z_lim .lt. ES%Z_xpoint(2)) ) then
          call is_axis_psi_mininum(node_list, element_list, bnd_elm_list)
          if (ES%axis_is_psi_minimum) then
            ES%psi_bnd = min(ES%psi_lim,ES%psi_bnd)
          else
            ES%psi_bnd = max(ES%psi_lim,ES%psi_bnd)
          endif
          write(*,'(A,4f8.3)') ' LIMITER PLASMA ',ES%psi_lim, ES%psi_bnd, ES%R_lim,ES%Z_lim
        endif
      endif

      if (freeb_equil_iterate_area .and. (.not. xpoint2)) then
        call iterate2area(node_list,element_list, ES%psi_axis, ES%psi_lim, xpoint2, xcase2, area_ref, ES%psi_bnd)
      endif
      
      write(*,'(A,1f8.3)') ' Psi_bnd = ', ES%psi_bnd   
      
      ! Calculate current feedback
      call integral_current(node_list,element_list,ES%psi_axis, ES%psi_bnd, xpoint2, xcase2, ES%Z_xpoint, current_tot)
  
      current_int = current_int + (current_tot-current_ref)
      
      if ((mod(iter,n_feedback_current) .eq. 0) .and. (.not. newton_GS_freebnd)) then
        current_FB_fact  = current_FB_fact * (1. - FB_Ip_position * (current_tot-current_ref)/current_ref &
                                                 - FB_Ip_integral *  current_int/current_ref   )
      endif
      
      !-------------- Multiplying FF' and p' profiles by the same factor to scale total current -------------------------
      FF_0 = FF_0_old * current_FB_fact   
      FF_1 = FF_1_old * current_FB_fact      
        
      T_0  = T_0_old  * current_FB_fact    
      T_1  = T_1_old  * current_FB_fact
      !------------------------------------------------------------------------------------------------------------------
      
      write(*,'(A,1e12.4)') 'Current Feedback factor = ',  current_FB_fact
      
      !Vertical feedback - needed for vertically unstable plasmas        
      Z_axis_int = Z_axis_int + (ES%Z_axis - Z_axis_ref)
      if (iter .eq. 1) then
        dZ_axis = 0.d0
      else
        dZ_axis = ES%Z_axis - Z_axis_old
      end if
    
      if ((mod(iter,n_feedback_vertical) .eq. 0) .and. (iter .ge. start_VFB) .and. (.not. newton_GS_freebnd) ) then
        vertical_FB = FB_Zaxis_position   * (ES%Z_axis-Z_axis_ref) &   ! vertical_FB is used in vacuum_equilibrium.f90 to modify the coils current
                    + FB_Zaxis_integral   * Z_axis_int          &   
                    + FB_Zaxis_derivative * dZ_axis
      endif
        
        Z_axis_old = ES%Z_axis
       
    end if ! my_id == 0
    
    call MPI_bcast(vertical_FB, 1, MPI_DOUBLE_PRECISION,  0, MPI_COMM_WORLD,ierr)
  
    ! --- Iterate equation
    call poisson(my_id,-1,node_list,element_list,bnd_node_list,bnd_elm_list,3,1,1, &
                 ES%psi_axis,ES%psi_bnd,xpoint2,xcase2,ES%Z_xpoint,freeboundary_equil,refinement,iter)   !----------- for GS use -1
  
  !  call boundary_check
   
    if (my_id == 0) then
      diff = 0.d0
      do i=1, node_list%n_nodes
        diff = diff + abs(node_list%node(i)%deltas(1,1,1))
      enddo  
      diff = diff / float(node_list%n_nodes)
    
      write(*,'(A,i5,e14.6)') ' iteration, diff : ',iter,diff
    end if ! my_id == 0
  
    call MPI_bcast(diff, 1, MPI_DOUBLE_PRECISION,  0, MPI_COMM_WORLD,ierr)
  
    if ( (iter > 1) .and. (diff < equil_accuracy_freeb) ) then
      if (my_id == 0) write(*,'(A,I4,A)') ' Free boundary equilibrium converged: after', iter, ' iterations'
      exit
    else if (iter == n_iter_freeb) then
      if (my_id == 0) write(*,'(A,ES10.3)') ' WARNING: Free boundary equilibrium not fully converged: diff=', diff
      exit
    end if
  
  enddo
  
  if (freeb_equil_iterate_area .and. (.not. xpoint2)) then
    n_limiter = 1  ! set found limiter (defined inside iterate2area)
  endif

else
  
  psi_offset_freeb = 0.d0
  
endif

if (my_id == 0) then
  ! Update psi axis and boundary with new values from the last iteration of equilibrium solvers 
  call find_axis(my_id,node_list,element_list,psi_axis,R_axis,Z_axis,i_elm_axis,s_axis,t_axis,ifail)
        
  write(10,'(i6,9e20.12)') iter, current_tot, R_axis, Z_axis, psi_bnd-psi_axis
  
  if ((ifail .ne. 0) .and. (iter .le. 5)) then
    call find_RZ(node_list,element_list,R_geo,Z_geo,R_out,Z_out,i_elm,s_out,t_out,ifail)
    call interp(node_list,element_list,i_elm,1,1,s_out,t_out,psi_axis,P_s,P_t,P_st,P_ss,P_tt)
    write(*,*)  ' changed magnetic axis to :  ', R_out,Z_out,psi_axis
  endif
  
  psi_bnd = 0.d0
  
  if (xpoint2) then
    call find_xpoint(my_id,node_list,element_list,psi_xpoint,R_xpoint,Z_xpoint,i_elm_xpoint,s_xpoint,t_xpoint,xcase2,ifail)
    if (ifail .ne. 1) then      
      psi_bnd  = psi_xpoint(1)
      if( (xcase2 .eq. 2) .or. ((xcase2 .eq. 3) .and. (abs(psi_xpoint(2)-psi_axis) .lt. abs(psi_xpoint(1)-psi_axis))) ) then
        psi_bnd = psi_xpoint(2)
      endif
      if(xcase2 .eq. LOWER_XPOINT) Z_xpoint(2) = +99.d0
      if(xcase2 .eq. UPPER_XPOINT) Z_xpoint(1) = -99.d0
    else
      Z_xpoint(1) = -99.d0 
      Z_xpoint(2) = +99.d0
    endif
  endif
  if (.not. xpoint2) then
    call find_limiter(my_id,node_list,element_list,bnd_elm_list,psi_lim,R_lim,Z_lim)
    if ( (Z_lim .gt. Z_xpoint(1)) .and. (Z_lim .lt. Z_xpoint(2)) ) then
      call is_axis_psi_mininum(node_list, element_list, bnd_elm_list)
      if (ES%axis_is_psi_minimum) then
        psi_bnd = min(psi_lim,psi_bnd)
      else
        psi_bnd = max(psi_lim,psi_bnd)
      endif
      write(*,'(A,4f8.3)') ' LIMITER PLASMA ',psi_lim, psi_bnd, R_lim,Z_lim
    endif
  endif

  if (freeboundary_equil .and. freeb_equil_iterate_area .and. (.not. xpoint2)) then
    call iterate2area(node_list,element_list, psi_axis, psi_lim, xpoint2, xcase2, area_ref, psi_bnd)
  endif
  
  !------------------------------- end of equilibrium, start filling data
  psi_axis = psi_axis - psi_offset_freeb
  psi_bnd  = psi_bnd  - psi_offset_freeb

  do i=1,node_list%n_nodes
  
    node_list%node(i)%values(1,1,1) = node_list%node(i)%values(1,1,1) - psi_offset_freeb
    psi = node_list%node(i)%values(1,1,1)
    R   = node_list%node(i)%x(1,1,1)
    Z   = node_list%node(i)%x(1,1,2)
  
    call density(    xpoint2, xcase2, Z, Z_xpoint, psi,psi_axis,psi_bnd,zn,dn_dpsi,dn_dz,dn_dpsi2,dn_dz2,             &
                                                               dn_dpsi_dz,dn_dpsi3,dn_dpsi_dz2, dn_dpsi2_dz)
  
    if (with_TiTe) then
      call temperature_i(xpoint2, xcase2, Z, Z_xpoint, psi,psi_axis,psi_bnd, &
    		     zTi,dTi_dpsi,dTi_dz,dTi_dpsi2,dTi_dz2,dTi_dpsi_dz,dTi_dpsi3,dTi_dpsi_dz2, dTi_dpsi2_dz)
  
      call temperature_e(xpoint2, xcase2, Z, Z_xpoint, psi,psi_axis,psi_bnd, &
    		     zTe,dTe_dpsi,dTe_dz,dTe_dpsi2,dTe_dz2,dTe_dpsi_dz,dTe_dpsi3,dTe_dpsi_dz2, dTe_dpsi2_dz)
      zT  	= zTi	       + zTe
      dT_dpsi	= dTi_dpsi     + dTe_dpsi
      dT_dpsi2	= dTi_dpsi2    + dTe_dpsi2
      dT_dpsi3	= dTi_dpsi3    + dTe_dpsi3
      dT_dz	= dTi_dz       + dTe_dz
      dT_dz2	= dTi_dz2      + dTe_dz2
      dT_dpsi_dz  = dTi_dpsi_dz  + dTe_dpsi_dz
      dT_dpsi2_dz = dTi_dpsi2_dz + dTe_dpsi2_dz
      dT_dpsi_dz2 = dTi_dpsi_dz2 + dTe_dpsi_dz2 
    else
      call temperature(xpoint2, xcase2, Z, Z_xpoint, psi,psi_axis,psi_bnd, &
    		     zT,dT_dpsi,dT_dz,dT_dpsi2,dT_dz2,dT_dpsi_dz,dT_dpsi3,dT_dpsi_dz2, dT_dpsi2_dz)
    endif
  
#ifdef fullmhd
      call F_profile(xpoint2, xcase2, Z, Z_xpoint, psi,psi_axis,psi_bnd,F_prof,dF_dpsi      ,dF_dz      , &
                                                                  dF_dpsi2    ,dF_dz2       ,dF_dpsi_dz , &
                                                                  zFFprime    ,dFFprime_dpsi,dFFprime_dz, &
                                                                  dFFprime_dpsi2,dFFprime_dz2 ,dFFprime_dpsi_dz)
  
      node_list%node(i)%Fprof_eq(1) =   F_prof
      node_list%node(i)%Fprof_eq(2) =   dF_dpsi * node_list%node(i)%values(1,2,var_A3)  + dF_dz * node_list%node(i)%x(1,2,2)
      node_list%node(i)%Fprof_eq(3) =   dF_dpsi * node_list%node(i)%values(1,3,var_A3)  + dF_dz * node_list%node(i)%x(1,3,2)
      node_list%node(i)%Fprof_eq(4) =   dF_dpsi * node_list%node(i)%values(1,4,var_A3)  + dF_dz * node_list%node(i)%x(1,4,2)      &
                                      + dF_dpsi2 * node_list%node(i)%values(1,2,var_A3) * node_list%node(i)%values(1,3,var_A3)  &
                                      + dF_dz2   * node_list%node(i)%x(1,2,2) * node_list%node(i)%x(1,3,2)
#else
    call FFprime(    xpoint2, xcase2, Z, Z_xpoint, psi,psi_axis,psi_bnd,zFFprime,dFFprime_dpsi,dFFprime_dz, &
                                                               dFFprime_dpsi2,dFFprime_dz2, dFFprime_dpsi_dz, .true.)
#endif

  
    zjz     = zFFprime      - R*R *      (dn_dpsi    * zT + zn * dT_dpsi)
  
    dj_dpsi = dFFprime_dpsi - R*R *      (dn_dpsi2   * zT + zn * dT_dpsi2  + 2.d0 * dn_dpsi * dT_dpsi)
  
    dj_dR   =               - 2.d0 * R * (dn_dpsi    * zT + zn * dT_dpsi)
  
    dj_dZ   = dFFprime_dz   - R*R *      (dn_dpsi_dz * zT + dn_dpsi * dT_dz + zn * dT_dpsi_dz + dn_dz * dT_dpsi)
  
    dj_dR_dR = - 2.d0     * (dn_dpsi     * zT + zn * dT_dpsi)
  
    dj_dZ_dZ = dFFprime_dz2   - R*R * ( dn_dpsi_dz2 * zT   + dn_dpsi_dz * dT_dz  + dn_dz * dT_dpsi_dz  + dn_dz2 * dT_dpsi &
                                      +  dn_dpsi_dz  * dT_dz + dn_dpsi    * dT_dz2 + zn    * dT_dpsi_dz2 + dn_dz  * dT_dpsi_dz)
  
    dj_dpsi2 = dFFprime_dpsi2 - R*R * (dn_dpsi3 * zT + 3.d0 * dn_dpsi * dT_dpsi2 + 3.d0 * dn_dpsi2 * dT_dpsi + zn * dT_dpsi3 )
  
    dj_dR_dZ   = - 2.d0 * R * (dn_dpsi_dz * zT + dn_dpsi * dT_dz + zn * dT_dpsi_dz + dn_dz * dT_dpsi)
  
    dj_dR_dpsi = - 2.d0 * R * (dn_dpsi2   * zT + zn * dT_dpsi2   + 2.d0 * dn_dpsi * dT_dpsi)
  
    dj_dZ_dpsi = dFFprime_dpsi_dz - R*R * ( dn_dpsi2_dz * zT    + dn_dz * dT_dpsi2     + 2.d0 * dn_dpsi_dz * dT_dpsi  &
                                            + dn_dpsi2    * dT_dz + zn    * dT_dpsi2_dz  + 2.d0 * dn_dpsi    * dT_dpsi_dz)
  
  
    node_list%node(i)%values(1,1,3) = zjz
  
    node_list%node(i)%values(1,2,3) = dj_dpsi * node_list%node(i)%values(1,2,1) &
                                    + dj_dR   * node_list%node(i)%x(1,2,1)        &
                                    + dj_dZ   * node_list%node(i)%x(1,2,2)
  
    node_list%node(i)%values(1,3,3) = dj_dpsi * node_list%node(i)%values(1,3,1) &
                                    + dj_dR   * node_list%node(i)%x(1,3,1)        &
                                    + dj_dZ   * node_list%node(i)%x(1,3,2)
  
    node_list%node(i)%values(1,4,3) = dj_dpsi  * node_list%node(i)%values(1,4,1) &
                                    + dj_dR    * node_list%node(i)%x(1,4,1)        &
                                    + dj_dZ    * node_list%node(i)%x(1,4,2)        &
                                    + dj_dR_dR * node_list%node(i)%x(1,2,1) * node_list%node(i)%x(1,3,1)  &
                                    + dj_dZ_dZ * node_list%node(i)%x(1,2,2) * node_list%node(i)%x(1,3,2)  &
                                    + dj_dpsi2 * node_list%node(i)%values(1,2,1) * node_list%node(i)%values(1,3,1)  &
                                    + dj_dR_dZ * ( node_list%node(i)%x(1,2,1) * node_list%node(i)%x(1,3,2)          &
                                                 + node_list%node(i)%x(1,3,1) * node_list%node(i)%x(1,2,2) )        &
                                    + dj_dR_dpsi*( node_list%node(i)%x(1,2,1) * node_list%node(i)%values(1,3,1)   &
                                                 + node_list%node(i)%x(1,3,1) * node_list%node(i)%values(1,2,1) ) &
                                    + dj_dZ_dpsi*( node_list%node(i)%x(1,2,2) * node_list%node(i)%values(1,3,1)   &
                                                 + node_list%node(i)%x(1,3,2) * node_list%node(i)%values(1,2,1) )
  
  enddo
  
  ! --- Find flux surfaces and plot them; determine the q-profile.  
  if (xpoint2 .and. (n_flux .gt. 1)) then
    
    call find_axis(my_id,node_list,element_list,psi_axis,R_axis,Z_axis,i_elm_axis,s_axis,t_axis,ifail)
    call find_xpoint(my_id,node_list,element_list,psi_xpoint,R_xpoint2,Z_xpoint2,i_elm_xpoint,s_xpoint,t_xpoint,xcase2,ifail)
    
    n_grids = 0
    sigmas  = 0.d0
    
    ! Build up some arrays to send as routine parameters to define_flux_values
    sigmas(1)  = SIG_closed  ; sigmas(2)  = SIG_theta
    sigmas(3)  = SIG_open    ; sigmas(4)  = SIG_outer   ; sigmas(5)  = SIG_inner
    sigmas(6)  = SIG_private ; sigmas(7)  = SIG_up_priv
    sigmas(8)  = SIG_leg_0   ; sigmas(9)  = SIG_leg_1
    sigmas(10) = SIG_up_leg_0; sigmas(11) = SIG_up_leg_1
    sigmas(12) = dPSI_open   ; sigmas(13) = dPSI_outer  ; sigmas(14) = dPSI_inner
    sigmas(15) = dPSI_private; sigmas(16) = dPSI_up_priv
  
    n_grids(1) = 2*n_flux   ; n_grids(2) = n_tht
    n_grids(3) = 2*n_open   ; n_grids(4) = 2*n_outer  ; n_grids(5) = 2*n_inner
    n_grids(6) = 2*n_private; n_grids(7) = 2*n_up_priv
    n_grids(8) = n_leg      ; n_grids(9) = n_up_leg
    if (xcase .eq. LOWER_XPOINT) then
      n_grids(4) = 0
      n_grids(5) = 0
      n_grids(7) = 0
      n_grids(9) = 0
    endif
    if (xcase .eq. UPPER_XPOINT) then
      n_grids(4) = 0
      n_grids(5) = 0
      n_grids(6) = 0
      n_grids(8) = 0
    endif
  
    ! Allocate surface_list structure (that's for plotting only)
    if (xcase2 .eq. LOWER_XPOINT) surface_list%n_psi = 2*n_flux + 2*n_open + 2*n_private
    if (xcase2 .eq. UPPER_XPOINT) surface_list%n_psi = 2*n_flux + 2*n_open + 2*n_up_priv
    if (xcase2 .eq. DOUBLE_NULL ) surface_list%n_psi = 2*n_flux + 2*n_open + 2*n_outer + 2*n_inner + 2*n_private + 2*n_up_priv
    if (allocated(surface_list%psi_values)) call tr_deallocate(surface_list%psi_values,"surface_list%psi_values",CAT_GRID)
    call tr_allocate(surface_list%psi_values,1,surface_list%n_psi,"surface_list%psi_values",CAT_GRID)
    
    ! Allocate sep_list structure (that's for plotting only)  
    sep_list%n_psi =3
    if(xcase .eq. DOUBLE_NULL) sep_list%n_psi =6
    if (allocated(sep_list%psi_values)) call tr_deallocate(sep_list%psi_values,"sep_list%psi_values",CAT_GRID)
    call tr_allocate(sep_list%psi_values,1,sep_list%n_psi,"sep_list%psi_values",CAT_GRID)
    
    ! Define the flux values to be plotted...
    psi_axis = psi_axis+0.01 !Just offset a little, because finding surfaces along the side of an element (on the xpoint grid) can be hard...
    call define_flux_values(node_list, element_list, surface_list, sep_list, xcase2, psi_xpoint, n_grids, sigmas)
    psi_axis = psi_axis-0.01 !Put it back, it's not used anyway, but just for principle!
    
  else
    surface_list%n_psi = 200  
    if (allocated(surface_list%psi_values)) call tr_deallocate(surface_list%psi_values,"surface_list%psi_values",CAT_GRID)
    call tr_allocate(surface_list%psi_values,1,surface_list%n_psi,"surface_list%psi_values",CAT_GRID)
    
    do i = 1, surface_list%n_psi
      surface_list%psi_values(i) = 1.25d0*(float(i)/float(surface_list%n_psi))**2 * (psi_bnd - psi_axis) + psi_axis
    enddo
    
    call find_flux_surfaces(my_id,xpoint2,xcase2,node_list,element_list,surface_list)
  
    sep_list%n_psi =1
    if (allocated(sep_list%psi_values)) call tr_deallocate(sep_list%psi_values,"sep_list%psi_values",CAT_GRID)
    call tr_allocate(sep_list%psi_values,1,sep_list%n_psi,"sep_list%psi_values",CAT_GRID)
    sep_list%psi_values(1) = psi_bnd
  
    call find_flux_surfaces(my_id,xpoint2,xcase2,node_list,element_list,sep_list)
  endif
  
  if (freeboundary_equil) then
    !call plot_coils(.true.)
    call plot_flux_surfaces(node_list,element_list,surface_list,.false.,4,xpoint2,xcase2)
    call plot_flux_surfaces(node_list,element_list,sep_list,.false.,1,xpoint2,xcase2)
  
    call plot_flux_surfaces(node_list,element_list,surface_list,.true.,4,xpoint2,xcase2)
    call plot_flux_surfaces(node_list,element_list,sep_list,.false.,1,xpoint2,xcase2)
    !call plot_coils(.false.)
  else
    if (xpoint2 .and. (n_flux .gt. 1)) then
      call plot_flux_surfaces(node_list,element_list,surface_list,.true.,1,xpoint2,xcase2)
      call plot_flux_surfaces(node_list,element_list,sep_list,.false.,1,xpoint2,xcase2)
    else
      call plot_flux_surfaces(node_list,element_list,surface_list,.true.,1,.false.,0)
      call plot_flux_surfaces(node_list,element_list,sep_list,.false.,1,xpoint2,xcase2)
    endif
  endif
  
  if (nice_q) then
    if (xpoint2) then
      ES%xpoint = xpoint2
      ES%Z_xpoint = Z_xpoint
    endif
    ES%psi_bnd  = psi_bnd
    ES%psi_axis = psi_axis
    ES%Z_xpoint = Z_xpoint
    ES%xpoint   = xpoint
    ES%xcase    = xcase
    call q_profile(node_list,element_list,surface_list,psi_axis,psi_bnd,psi_xpoint,Z_xpoint)
  endif
  
  !================ Temperature and density profiles =f(psi_norm) similar to q(psi_norm) needed to calculate neoclassical coef===========
  if (allocated(T_profile)) call tr_deallocate(T_profile,"T_profile",CAT_GRID)
  call tr_allocate(T_profile,1,surface_list%n_psi,"T_profile",CAT_GRID)
  if (allocated(density_profile)) call tr_deallocate(density_profile,"density_profile",CAT_GRID)
  call tr_allocate(density_profile,1,surface_list%n_psi,"density_profile",CAT_GRID)
  
  do i=2,surface_list%n_psi
    psi= surface_list%psi_values(i)
    
    if (with_TiTe) then
      call temperature_i(xpoint2, xcase2, Z, Z_xpoint, psi,psi_axis,psi_bnd, &
           Ti_prof,dTi_dpsi,dTi_dz,dTi_dpsi2,dTi_dz2,dTi_dpsi_dz,dTi_dpsi3,dTi_dpsi_dz2, dTi_dpsi2_dz)
  
      call temperature_e(xpoint2, xcase2, Z, Z_xpoint, psi,psi_axis,psi_bnd, &
           Te_prof,dTe_dpsi,dTe_dz,dTe_dpsi2,dTe_dz2,dTe_dpsi_dz,dTe_dpsi3,dTe_dpsi_dz2, dTe_dpsi2_dz)
      T_prof      = Ti_prof      + Te_prof
      dT_dpsi     = dTi_dpsi     + dTe_dpsi
      dT_dpsi2    = dTi_dpsi2    + dTe_dpsi2
      dT_dpsi3    = dTi_dpsi3    + dTe_dpsi3
      dT_dz       = dTi_dz       + dTe_dz
      dT_dz2      = dTi_dz2      + dTe_dz2
      dT_dpsi_dz  = dTi_dpsi_dz  + dTe_dpsi_dz
      dT_dpsi2_dz = dTi_dpsi2_dz + dTe_dpsi2_dz
      dT_dpsi_dz2 = dTi_dpsi_dz2 + dTe_dpsi_dz2 
    else
      call temperature(.false.,xcase2,0., Z_xpoint, psi,psi_axis,psi_bnd,T_prof,dT_dpsi,dT_dz,dT_dpsi2,dT_dz2,             &
        dT_dpsi_dz,dT_dpsi3,dT_dpsi_dz2, dT_dpsi2_dz)
    endif

    call density( .false., xcase2,0., Z_xpoint, psi,psi_axis,psi_bnd,density_prof,dn_dpsi,dn_dz,dn_dpsi2,dn_dz2,             &
          dn_dpsi_dz,dn_dpsi3,dn_dpsi_dz2, dn_dpsi2_dz)
  
    T_profile(i)=T_prof
    density_profile(i)=density_prof
  end do
  
  write(*,*) '***************************************'
  write(*,*) 'output T and rho profiles (in JOREK units) for neoclassical profile calculation'
  ! --- Write out T and rho profiles to "T_rho_profiles.dat".
  open(432, file='T_rho_profiles.dat', action='write', status='replace')
  do i=2, surface_list%n_psi
     write(432,'(3ES13.5)') T_profile(i), density_profile(i)
  end do
  close(432)
  !========================= end modif ===========================================
  
  if (allocated(surface_list%psi_values))    call tr_deallocate(surface_list%psi_values,"surface_list%psi_values",CAT_GRID)
  if (allocated(surface_list%flux_surfaces)) deallocate(surface_list%flux_surfaces)
  if (allocated(sep_list%psi_values))        call tr_deallocate(sep_list%psi_values,"sep_list%psi_values",CAT_GRID)
  if (allocated(sep_list%flux_surfaces))     deallocate(sep_list%flux_surfaces)
  
  if (allocated(T_profile)) call tr_deallocate(T_profile,"T_profile",CAT_GRID)
  if (allocated(density_profile)) call tr_deallocate(density_profile,"density_profile",CAT_GRID)
  
end if ! my_id == 0

if (freeboundary_equil) then
  call broadcast_elements(my_id, element_list)
  call broadcast_nodes(my_id, node_list)  !--- This is required for boundary_check
  call broadcast_boundary(my_id, bnd_elm_list, bnd_node_list)
  call boundary_check(my_id)
endif
return
end subroutine equilibrium
