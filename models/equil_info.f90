!> Defines a data structure to collect all "equilibrium state information" of the simulation.
!!
!! This includes the location of magnetic axis, X-points, limiter point etc.
!!
!! The module also contains routines/functions for tasks like Psi_N calculation etc.
!!
!! @TODO: Better calculation of strike points (just place-holder right now) including i_elm_strike
module equil_info
  
  
  
  use constants,      only: LOWER_XPOINT, UPPER_XPOINT, DOUBLE_NULL,SYMMETRIC_XPOINT
  use data_structure, only: type_node_list, type_element_list, type_bnd_element_list
  use phys_module,    only: R_geo, Z_geo, FF_0, psi_axis_t, psi_bnd_t, Z_xpoint_t, index_now, SDN_threshold
  use mod_interp
  
  
  
  implicit none
  
  
  
  public
  
  
  
  !> This data structure contains information about the position of axis, limiter point,
  !! X-point(s), and strike points. It needs to be updated after each time step and after
  !! each iteration of the equilibrium calculation using the routine update_equil_state.
  type t_equil_state
    
    logical          :: initialized = .false.
    
    logical          :: limiter_plasma           !< Is the current state a limiter plasma?
    logical          :: axis_is_psi_minimum      !< Is psi_axis < psi_bnd or > psi_bnd?
    
    ! --- Magnetic Axis
    real*8           :: R_axis                   !< R coordinate of axis.
    real*8           :: Z_axis                   !< Z coordinate of axis.
    real*8           :: Psi_axis                 !< Poloidal flux value at axis.
    real*8           :: Psi_axis_init            !< Poloidal flux value at axis in GS-equilibrium at t=0.
    integer          :: i_elm_axis               !< Index of element containing the axis.
    real*8           :: s_axis                   !< s coordinate of axis within element.
    real*8           :: t_axis                   !< t coordinate of axis within element.
    integer          :: ifail_axis               !< Error code for axis determination.
    
    ! --- Limiter Point
    real*8           :: R_lim                    !< R coordinate of limiter point.
    real*8           :: Z_lim                    !< Z coordinate of limiter point.
    real*8           :: Psi_lim                  !< Poloidal flux value at limiter point.
    integer          :: i_elm_lim                !< Index of element containing limiter point.
    real*8           :: s_lim                    !< s coordinate of limiter point within element.
    real*8           :: t_lim                    !< t coordinate of limiter point within element.
    integer          :: ifail_lim                !< Error code for limiter determination.
    
    ! --- X-Point(s)
    logical          :: xpoint                   !< Is this an X-point case? Not necessarily active!
    integer          :: xcase                    !< Upper/lower/double X-point?
    integer          :: active_xpoint            !< Which X-point is active?
    real*8           :: R_xpoint(2)              !< R coordinate of X-point(s).
    real*8           :: Z_xpoint(2)              !< Z coordinate of X-point(s).
    real*8           :: Z_xpoint_init(2)         !< Z coordinate of X-point(s) in GS-equilibrium at t=0.
    real*8           :: Psi_xpoint(2)            !< Poloidal flux value of X-point(s).
    integer          :: i_elm_xpoint(2)          !< Index of element containing the X-point.
    real*8           :: s_xpoint(2)              !< s coordinate of X-point within element.
    real*8           :: t_xpoint(2)              !< t coordinate of X-point within element.
    integer          :: ifail_xpoint             !< Error code for X-point determination.
    
    ! --- Boundary point (point defining the plasma LCFS, either the active limiter point or X-point)
    real*8           :: R_bnd                    !< R coordinate of boundary point.
    real*8           :: Z_bnd                    !< Z coordinate of boundary point.
    real*8           :: Psi_bnd                  !< Psi of the boundary point (Psi of the LCFS)
    real*8           :: Psi_bnd_init             !< Psi of the boundary point (Psi of the LCFS) in GS-equilibrium at t=0
    integer          :: i_elm_bnd                !< Index of element containing the boundary point
    real*8           :: s_bnd                    !< s coordinate of the boundary point within element.
    real*8           :: t_bnd                    !< t coordinate of the boundary point within element.
    integer          :: ifail_bnd                !< Error code for determination of boundary point
    
    ! --- Strike Point(s) derived from axisymmetric field component.
    integer          :: num_strike               !< Number of strike points.
    real*8           :: R_strike(99)             !< R coordinate of strike point(s).
    real*8           :: Z_strike(99)             !< Z coordinate of strike point(s).
    integer          :: i_bndelm_strike(99)      !< Index of boundary element containing strike pt.
    real*8           :: s_strike(99)             !< s coordinate of strike pt within boundary elem.
    
    ! --- Inner/Outer points on the midplane close to the boundary of the computational domain.
    real*8           :: R_midpl(2)               !< R coordinate of "midplane points".

    ! --- Plasma shape parameters as defined in T. Luce, PPCF 55 (2013) 095009, equations (1-6)
    real*8           :: LCFS_Rgeo                !< Major radius
    real*8           :: LCFS_a                   !< Minor radius
    real*8           :: LCFS_epsilon             !< Inverse aspect ratio 
    real*8           :: LCFS_kappa               !< Elongation
    real*8           :: LCFS_deltaU              !< Upper triangularity
    real*8           :: LCFS_deltaL              !< Lower triangularity 
    
  end type t_equil_state
  

  type(t_equil_state), target   :: ES  
  
  
  contains
  
  
  
  !> Re-calculate the equilibrium state.
  subroutine update_equil_state(my_id, node_list, element_list, bnd_elm_list, xpoint, xcase)
    
    use phys_module,    only: freeboundary

    ! --- Routine parameters.
    integer,                     intent(in)    :: my_id
    type(type_node_list),        intent(in)    :: node_list
    type(type_element_list),     intent(in)    :: element_list
    type(type_bnd_element_list), intent(in)    :: bnd_elm_list
    logical                                    :: xpoint
    integer                                    :: xcase
    
    ! --- Local variables.
    integer :: my_id_fake, i_out, ifail, i, mv1
    real*8  :: R_out, Z_out, s_out, t_out, R1, R2, dR, R_s, R_t, Z_s, Z_t
    real*8  :: P_s, P_t, P_st, P_ss, P_tt
    
    my_id_fake  = 9999
    
    ! --- Find the magnetic axis.
    call find_axis(my_id_fake, node_list, element_list, ES%psi_axis, ES%R_axis, ES%Z_axis,              &
      ES%i_elm_axis, ES%s_axis, ES%t_axis, ES%ifail_axis)
      
    ! --- Find out if the axis is a minimum or a maximum of the poloidal flux (required for find_limiter)    
    if (.not. ES%initialized) call is_axis_psi_mininum(node_list, element_list, bnd_elm_list)
        
    ! --- Find the X-point(s).
    ES%xpoint       = xpoint
    ES%xcase        = xcase
    ES%ifail_xpoint = 0
    if ( xpoint ) call find_xpoint(my_id_fake, node_list, element_list, ES%psi_xpoint, ES%R_xpoint,     &
      ES%Z_xpoint, ES%i_elm_xpoint, ES%s_xpoint, ES%t_xpoint, ES%xcase, ES%ifail_xpoint)
    
    ! --- Find the limiter point.
    ES%ifail_lim = 0
    call find_limiter(my_id_fake, node_list, element_list, bnd_elm_list, ES%psi_lim, ES%R_lim, ES%Z_lim)
    call find_RZ(node_list, element_list, ES%R_lim, ES%Z_lim, R_out, Z_out, ES%i_elm_lim, ES%s_lim,&
      ES%t_lim, ES%ifail_lim)
    
    if ( xpoint ) then ! (X-point plasma)
      
      if ( (xcase==LOWER_XPOINT) ) then
        
        ES%psi_bnd        = ES%psi_xpoint(1)
        ES%limiter_plasma = .false.
        ES%active_xpoint  = LOWER_XPOINT
        
      else if ( (xcase==UPPER_XPOINT) ) then
        
        ES%psi_bnd        = ES%psi_xpoint(2)
        ES%limiter_plasma = .false.
        ES%active_xpoint  = UPPER_XPOINT
        
      else if ( (xcase==DOUBLE_NULL) ) then
        
        ES%limiter_plasma = .false.

        if ( abs(ES%psi_axis-ES%psi_xpoint(1)) < abs(ES%psi_axis-ES%psi_xpoint(2)) ) then
          ES%psi_bnd       = ES%psi_xpoint(1)
          ES%active_xpoint = LOWER_XPOINT

          ! --- Save boundary point inforamtion for DOUBLE_NULL cases       
          ES%R_bnd      =  ES%R_xpoint(ES%active_xpoint)
          ES%Z_bnd      =  ES%Z_xpoint(ES%active_xpoint)
          ES%i_elm_bnd  =  ES%i_elm_xpoint(ES%active_xpoint)
          ES%s_bnd      =  ES%s_xpoint(ES%active_xpoint)
          ES%t_bnd      =  ES%t_xpoint(ES%active_xpoint)
          ES%ifail_bnd  =  ES%ifail_xpoint
        else
          ES%psi_bnd       = ES%psi_xpoint(2)
          ES%active_xpoint = UPPER_XPOINT

          ! --- Save boundary point inforamtion for DOUBLE_NULL cases       
          ES%R_bnd      =  ES%R_xpoint(ES%active_xpoint)
          ES%Z_bnd      =  ES%Z_xpoint(ES%active_xpoint)
          ES%i_elm_bnd  =  ES%i_elm_xpoint(ES%active_xpoint)
          ES%s_bnd      =  ES%s_xpoint(ES%active_xpoint)
          ES%t_bnd      =  ES%t_xpoint(ES%active_xpoint)
          ES%ifail_bnd  =  ES%ifail_xpoint
        end if

        ! If one want to generate a symmetric double-null grid
        if ( abs(ES%psi_xpoint(1)-ES%psi_xpoint(2)) < SDN_threshold ) then
          ES%active_xpoint = SYMMETRIC_XPOINT
        endif
        
      else ! This should never happen.
        write(*,*) 'ERROR: ILLEGAL VALUE FOR XCASE:', xcase
        stop
      end if
      
      ! --- Has the X-plasma changed to a limiter plasma?
      if (freeboundary) then
        if ( abs(ES%psi_axis-ES%psi_lim) < abs(ES%psi_axis-ES%psi_bnd) ) then
          ES%psi_bnd        = ES%psi_lim
          ES%limiter_plasma = .true.
          ES%active_xpoint  = 0
        endif 
      endif 
      
    else ! (limiter plasma)
      
      ES%limiter_plasma = .true.
      ES%active_xpoint  = 0
      ES%R_xpoint(:)    = R_geo
      ES%Z_xpoint(1)    = -99.d0
      ES%Z_xpoint(2)    =  99.d0

      !--- If there are no X-points and find_limiter has failed, assume that the grid's boundary is a flux-surface and a limiter
      if (ES%ifail_lim /= 0) then
   
        !--- select a random boundary point and choose it as limiter
        ES%i_elm_lim = bnd_elm_list%bnd_element(1)%element
        mv1          = bnd_elm_list%bnd_element(1)%side 
        if ((mv1 .eq. 1) .or. (mv1 .eq. 4)) then
          ES%s_lim = 0.d0;  ES%t_lim = 0.d0;  
        else if (mv1 .eq. 2) then
          ES%s_lim = 1.d0;  ES%t_lim = 0.d0;  
        else
          ES%s_lim = 0.d0;  ES%t_lim = 1.d0;
        endif               
        call interp(node_list, element_list, ES%i_elm_lim, 1, 1, ES%s_lim, ES%t_lim, ES%psi_lim, P_s, P_t, P_st, P_ss, P_tt)
        call interp_RZ(node_list, element_list, ES%i_elm_lim, ES%s_lim, ES%t_lim, ES%R_lim, R_s, R_t, ES%Z_lim, Z_s, Z_t)                

      endif    
      
      ES%psi_bnd        = ES%psi_lim
            
    end if
    
    ! --- psi_axis < psi_bnd or > psi_bnd?
    ES%axis_is_psi_minimum = ( ES%psi_axis < ES%psi_bnd )
    
    ! --- Save boundary point information
    if (ES%limiter_plasma) then
      ES%R_bnd      =  ES%R_lim
      ES%Z_bnd      =  ES%Z_lim
      ES%i_elm_bnd  =  ES%i_elm_lim
      ES%s_bnd      =  ES%s_lim
      ES%t_bnd      =  ES%t_lim
      ES%ifail_bnd  =  ES%ifail_lim
    else
      if (xcase .ne. DOUBLE_NULL) then
        ES%R_bnd      =  ES%R_xpoint(ES%active_xpoint)
        ES%Z_bnd      =  ES%Z_xpoint(ES%active_xpoint)
        ES%i_elm_bnd  =  ES%i_elm_xpoint(ES%active_xpoint)
        ES%s_bnd      =  ES%s_xpoint(ES%active_xpoint)
        ES%t_bnd      =  ES%t_xpoint(ES%active_xpoint)
        ES%ifail_bnd  =  ES%ifail_xpoint
      endif
    endif  
    
    ! --- Strike points.
    ES%num_strike          = 0
    ES%R_strike(:)         = -99.d0
    ES%Z_strike(:)         = -99.d0
    ES%i_bndelm_strike(:)  = -9999
    ES%s_strike(:)         = -99.d0
    call find_strike(node_list, bnd_elm_list, ES)
    
    ! --- Inner midplane point.
    R1 = ES%R_axis
    R2 = 0.d0
    i  = 0
    do
      i = i + 1
      call find_RZ(node_list, element_list, (R1+R2)/2.d0, ES%Z_axis, R_out, Z_out, i_out, s_out,   &
        t_out, ifail)
      if ( ifail == 0 ) then
        R1 = (R1 + R2) / 2.d0
      else
        R2 = (R1 + R2) / 2.d0
      end if
      if ( abs(R2-R1) < 1.d-4 ) exit
    end do
    ES%R_midpl(1) = (R2+R1)/2.d0
    
    ! --- Outer midplane point.
    R1 = ES%R_axis
    R2 = 3.d0 * ES%R_axis
    i  = 0
    do
      i = i + 1
      call find_RZ(node_list, element_list, (R1+R2)/2.d0, ES%Z_axis, R_out, Z_out, i_out, s_out,   &
        t_out, ifail)
      if ( ifail == 0 ) then
        R1 = (R1 + R2) / 2.d0
      else
        R2 = (R1 + R2) / 2.d0
      end if
      if ( abs(R2-R1) < 1.d-4 ) exit
    end do
    ES%R_midpl(2) = (R2+R1)/2.d0
    
    ! --- Save initial psi_axis and psi_bnd (at the moment only needed by model710)
    ! --- Note: this needs to be on my_id=0, because live variables are only available on main MPI process
    ! --- If we are not on my_id=0, it means this already happened at begining of run, and these were
    ! --- already broadcasted so we are fine (since this never changes in time)
    if (my_id .eq. 0) then
      if (index_now .gt. 1) then
        ES%psi_axis_init    = psi_axis_t(1)
        ES%psi_bnd_init     = psi_bnd_t(1)
        ES%Z_xpoint_init(:) = Z_xpoint_t(1,:)
      else
        ES%psi_axis_init    = ES%psi_axis
        ES%psi_bnd_init     = ES%psi_bnd
        ES%Z_xpoint_init(:) = ES%Z_xpoint(:)
      endif
    endif

    ! --- Calculate shape parameters of the LCFS
    call LCFS_shape_parameters(node_list,element_list)
    
    ES%initialized = .true.
    
  end subroutine update_equil_state
    
  
  

  
  
  
  !> Estimate if psi_axis is a minimum or a maximum of flux
  subroutine is_axis_psi_mininum(node_list, element_list, bnd_elm_list)
    
    ! --- Routine variables
    type(type_node_list),        intent(in)    :: node_list
    type(type_element_list),     intent(in)    :: element_list
    type(type_bnd_element_list), intent(in)    :: bnd_elm_list
                                                                     
    ! --- Local variables.
    real*8  :: P, P_s, P_t, P_st, P_ss, P_tt, R_t, Z_t, R_s, Z_s
    real*8  :: R_out, Z_out, s_out, t_out, R1, Z1, R2, Z2     
    real*8  :: psi_axis, R_axis, Z_axis, s_axis, t_axis
    integer :: i_elm, i_elm_out, i_elm_axis, ifail  
    
    ! --- Get coordinates of the magnetic axis
    if (ES%initialized) then
      R_axis   = ES%R_axis
      Z_axis   = ES%Z_axis
      psi_axis = ES%psi_axis
    else
      call find_axis(99, node_list, element_list, psi_axis, R_axis, Z_axis,              &
        i_elm_axis, s_axis, t_axis, ifail)
    endif
    
    ! --- Find random point (R,Z) coordinates at computational boundary
    i_elm = bnd_elm_list%bnd_element(1)%element 
    call interp_RZ(node_list, element_list, i_elm, 0.d0, 0.d0, R1, R_s, R_t, Z1, Z_s, Z_t)

    ! --- Find point between axis and bnd point (located 25% away from axis on the connecting line)
    R2 = R_axis + 0.25d0*(R1-R_axis)
    Z2 = Z_axis + 0.25d0*(Z1-Z_axis)
    call find_RZ(node_list, element_list, R2, Z2, R_out, Z_out, i_elm_out, s_out, t_out, ifail)    
    call interp(node_list, element_list, i_elm_out, 1, 1, s_out, t_out, P, P_s, P_t, P_st, P_ss, P_tt)

    ! --- Decide whether the axis is a minimum of psi
    if ( (P - psi_axis) > 0.d0 ) then
      ES%axis_is_psi_minimum = .true.
    else 
      ES%axis_is_psi_minimum = .false.
    endif
    
    if (ifail /= 0) then   ! if the reference point was not found, use FF_0 instead
      write(*,*) ' WARNING: is_axis_psi_minimum is failing (no intermediate point found)'
      write(*,*) ' Deciding if axis is minimum using the sign of FF_0'
      ES%axis_is_psi_minimum = .false.
      if (FF_0 >= 0.d0)  ES%axis_is_psi_minimum = .true.
    endif
    
  end subroutine is_axis_psi_mininum
  
  
  
  
  
  
  
  !> Readable output of the equilibrium state for the logfile.
  subroutine print_equil_state(verbose)
    
    ! --- Routine parameters.
    logical,             intent(in) :: verbose    !< Output much additional information?
    
    ! --- Local variables
    integer :: i
    
    if ( .not. ES%initialized ) then
      write(*,*) 'ERROR in mod_equil_info|print_equil_state: ES used before initialization.'
      write(*,*) 'Call update_equil_state first.'
      stop
    end if
    
    101 format(1x,a)
    102 format(1x,5(a,f10.5))
    103 format(1x,5(a,i10))
    104 format(1x,a,i3,2f10.5)
    105 format(1x,a,i3.3,a,f10.5)
    106 format(1x,a,i3.3,a,i10)
    
    ! --- General description of the plasma.
    write(*,*)
    write(*,*) '=============================================================='
    if ( ES%limiter_plasma ) then
      write(*,101) 'Plasma_type        = Limiter Plasma'
    else
      write(*,101) 'Plasma_type        = X-Point Plasma'
      if ( ES%xcase == LOWER_XPOINT ) then
        write(*,101) 'Xpoint_type        = Lower X-Point'
      else if ( ES%xcase == UPPER_XPOINT ) then
        write(*,101) 'Xpoint_type        = Upper X-Point'
      else
        write(*,101) 'Xpoint_type        = Double X-Point'
        if ( ES%active_xpoint == LOWER_XPOINT ) then
          write(*,101) 'Active_xpoint      = Lower X-Point'
        else
          write(*,101) 'Active_xpoint      = Upper X-Point'
        end if
      end if
    end if
    if (verbose ) then
      if ( ES%axis_is_psi_minimum ) then
        write(*,101) 'Psi_axis           = Minimum'
      else
        write(*,101) 'Psi_axis           = Maximum'
      end if
    end if
    write(*,102) 'Psi_bnd            =', ES%psi_bnd
    
    ! --- Magnetic axis.
    write(*,*) '--- Magnetic Axis --------------------------------------------'
    write(*,102) 'R_axis             =', ES%R_axis
    write(*,102) 'Z_axis             =', ES%Z_axis
    write(*,102) 'Psi_axis           =', ES%Psi_axis
    if ( verbose ) then
      write(*,103) 'i_elm_axis         =', ES%i_elm_axis
      write(*,102) 's_axis             =', ES%s_axis
      write(*,102) 't_axis             =', ES%t_axis
      write(*,103) 'ifail_axis         =', ES%ifail_axis
    end if
    
    ! --- Limiter point.
    if ( ES%limiter_plasma .or. verbose ) then
      write(*,*) '--- Limiter Point --------------------------------------------'
      write(*,102) 'R_lim              =', ES%R_lim
      write(*,102) 'Z_lim              =', ES%Z_lim
      write(*,102) 'Psi_lim            =', ES%Psi_lim
      if ( verbose ) then
        write(*,103) 'i_elm_lim          =', ES%i_elm_lim
        write(*,102) 's_lim              =', ES%s_lim
        write(*,102) 't_lim              =', ES%t_lim
        write(*,103) 'ifail_lim          =', ES%ifail_lim
      end if
    end if
    
    ! --- X-point(s).
    if ( ES%xpoint ) then
      if ( ( ES%xcase == LOWER_XPOINT ) .or. ( ES%xcase == DOUBLE_NULL ) ) then
        write(*,*) '--- Lower X-Point --------------------------------------------'
        write(*,102) 'R_xpoint1          =', ES%R_xpoint(1)
        write(*,102) 'Z_xpoint1          =', ES%Z_xpoint(1)
        write(*,102) 'Psi_xpoint1        =', ES%Psi_xpoint(1)
        if ( verbose ) then
          write(*,103) 'i_elm_xpoint1      =', ES%i_elm_xpoint(1)
          write(*,102) 's_xpoint1          =', ES%s_xpoint(1)
          write(*,102) 't_xpoint1          =', ES%t_xpoint(1)
          write(*,103) 'ifail_xpoint       =', ES%ifail_xpoint
        end if
      end if
      if ( ( ES%xcase == UPPER_XPOINT ) .or. ( ES%xcase == DOUBLE_NULL ) ) then
        write(*,*) '--- Upper X-Point --------------------------------------------'
        write(*,102) 'R_xpoint2          =', ES%R_xpoint(2)
        write(*,102) 'Z_xpoint2          =', ES%Z_xpoint(2)
        write(*,102) 'Psi_xpoint2        =', ES%Psi_xpoint(2)
        if ( verbose ) then
          write(*,103) 'i_elm_xpoint2      =', ES%i_elm_xpoint(2)
          write(*,102) 's_xpoint2          =', ES%s_xpoint(2)
          write(*,102) 't_xpoint2          =', ES%t_xpoint(2)
          write(*,103) 'ifail_xpoint       =', ES%ifail_xpoint
        end if
      end if
    end if

    ! --- Boundary point (point defining LCFS)
    write(*,*) '--- Boundary point (point defining LCFS) -------------------------'
    write(*,102) 'R_bnd              =', ES%R_bnd
    write(*,102) 'Z_bnd              =', ES%Z_bnd
    write(*,102) 'Psi_bnd            =', ES%Psi_bnd
    if ( verbose ) then
      write(*,103) 'i_elm_bnd          =', ES%i_elm_bnd
      write(*,102) 's_bnd              =', ES%s_bnd
      write(*,102) 't_bnd              =', ES%t_bnd
      write(*,103) 'ifail_bnd          =', ES%ifail_bnd
    end if
    
    ! --- Strike points.
    if ( verbose ) then
      do i = 1, ES%num_strike
        write(*,'(1x,"--- Strike Point",i3,1x,42("-"))') i
        write(*,105) 'R_strike', i, '        =', ES%R_strike(i)
        write(*,105) 'Z_strike', i, '        =', ES%Z_strike(i)
        write(*,106) 'i_bndelm_strike', i, ' =', ES%i_bndelm_strike(i)
        write(*,105) 's_strike', i, '        =', ES%s_strike(i)
      end do
    end if
    
    ! --- Midplane points.
    if ( verbose ) then
      write(*,*) '--- Midplane Points ------------------------------------------'
      write(*,102) 'R_midpl1           =', ES%R_midpl(1)
      write(*,102) 'R_midpl2           =', ES%R_midpl(2)
    end if

    ! --- Shaping parameters
    if ( verbose ) then
      write(*,*) '--- LCFS shape parameters (as in PPCF 55 (2013) 095009) ------'
      write(*,102) 'R_geo              =', ES%LCFS_Rgeo    
      write(*,102) 'a_min              =', ES%LCFS_a       
      write(*,102) 'epsilon            =', ES%LCFS_epsilon 
      write(*,102) 'kappa              =', ES%LCFS_kappa   
      write(*,102) 'delta_U            =', ES%LCFS_deltaU  
      write(*,102) 'delta_L            =', ES%LCFS_deltaL  
    end if
    
    write(*,*) '=============================================================='
    write(*,*)
    
  end subroutine print_equil_state
  
  
  
  !> Nice, readable output of the equilibrium state.
  subroutine save_special_points(filename, append, ioerr)
    
    ! --- Routine parameters.
    character(len=*),    intent(in)    :: filename   !< Output to which file?
    logical,             intent(in)    :: append     !< Append to existing file?
    integer,             intent(inout) :: ioerr      !< I/O Error code.
    
    ! --- Local variables
    integer, parameter :: ifile = 33 !### better solution?
    integer :: i
    
    if ( append ) then
      open(ifile, file=filename, form='formatted', access='append', iostat=ioerr)
      if ( ioerr /= 0 ) then
        write(ifile,*)
        write(ifile,*)
      end if
    else
      open(ifile, file=filename, form='formatted', status='replace', iostat=ioerr)
    end if
    
    if ( ioerr /= 0 ) then
      write(*,*) 'ERROR in mod_equil_info|save_special_points opening file "',trim(filename),'".'
      return
    end if
    
    write(ifile,*) '# Magnetic Axis'
    write(ifile,*) ES%R_axis, ES%Z_axis
    write(ifile,*)
    write(ifile,*)
    
    write(ifile,*) '# Limiter Point'
    write(ifile,*) ES%R_lim, ES%Z_lim
    write(ifile,*)
    write(ifile,*)
    
    if ( ES%xpoint ) then
      if ( (ES%xcase==LOWER_XPOINT) .or. (ES%xcase==DOUBLE_NULL) ) then
        write(ifile,*) '# Lower X-Point'
        write(ifile,*) ES%R_xpoint(1), ES%Z_xpoint(1)
        write(ifile,*)
        write(ifile,*)
      end if
      if ( (ES%xcase==UPPER_XPOINT) .or. (ES%xcase==DOUBLE_NULL) ) then
        write(ifile,*) '# Upper X-Point'
        write(ifile,*) ES%R_xpoint(2), ES%Z_xpoint(2)
        write(ifile,*)
        write(ifile,*)
      end if
    end if

    write(ifile,*) '# Boundary point (point defining LCFS)'
    write(ifile,*) ES%R_bnd, ES%Z_bnd
    write(ifile,*)
    write(ifile,*)
    
    do i = 1, ES%num_strike
      write(ifile,*) '# Strike point', i
      write(ifile,*) ES%R_strike(i), ES%Z_strike(i)
      write(ifile,*)
      write(ifile,*)
    end do
    
    write(ifile,*) '# Inner midplane point'
    write(ifile,*) ES%R_midpl(1), ES%Z_axis
    write(ifile,*)
    write(ifile,*)
    
    write(ifile,*) '# Outer midplane point'
    write(ifile,*) ES%R_midpl(2), ES%Z_axis
    write(ifile,*)
    write(ifile,*)
    
  end subroutine save_special_points
  
  
  
  !> Calculate Psi_N for given Psi.
  pure real*8 function get_psi_n( psi, Z )
    
    ! --- Routine parameters.
    real*8,              intent(in) :: psi                      !< Poloidal flux value.
    real*8,   optional,  intent(in) :: Z                        !< vertical position coordinate Z
    
    ! --- Local
    logical  :: correct_private
    real*8   :: psi_n_xpoint_upper, psi_n_xpoint_lower
    
    ! --- If the user specifies Z, then treat private regions specially
    if (present(Z)) then
      correct_private = .true.
    else
      correct_private = .false.
    endif
    
    get_psi_n = ( psi - ES%psi_axis ) / ( ES%psi_bnd - ES%psi_axis )
    
    if (ES%xpoint .and. correct_private) then

      if ( ES%xcase .ne. UPPER_XPOINT ) then
        psi_n_xpoint_lower = ( ES%psi_xpoint(1) - ES%psi_axis ) / ( ES%psi_bnd - ES%psi_axis )       
        if ( (get_psi_n < psi_n_xpoint_lower) .and. (Z < ES%Z_xpoint(1)) ) then   ! if true is lower private region
          get_psi_n = 2.d0*psi_n_xpoint_lower - get_psi_n
        endif
      endif
      if ( ES%xcase .ne. LOWER_XPOINT ) then
        psi_n_xpoint_upper = ( ES%psi_xpoint(2) - ES%psi_axis ) / ( ES%psi_bnd - ES%psi_axis )
        if ( (get_psi_n < psi_n_xpoint_upper) .and. (Z > ES%Z_xpoint(2)) ) then   ! if true is upper private region
          get_psi_n = 2.d0*psi_n_xpoint_upper - get_psi_n
        endif
      endif

    endif
    
  end function get_psi_n
  
  
  
  
  !> Broadcast equil_state information between MPI processes
  subroutine broadcast_equil_state(my_id)
    use mpi_mod
    implicit none
    
    
    ! --- Routine parameters
    integer, intent(in) :: my_id
    
    ! --- Local variables
    integer :: ierr

    call MPI_BCAST(ES%initialized,         1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(ES%limiter_plasma,      1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(ES%axis_is_psi_minimum, 1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
        
    ! --- Magnetic Axis
    call MPI_BCAST(ES%R_axis,       1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(ES%Z_axis,       1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(ES%Psi_axis,     1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(ES%Psi_axis_init,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(ES%s_axis,       1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(ES%t_axis,       1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(ES%i_elm_axis,   1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(ES%ifail_axis,   1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

    
    ! --- Limiter Point
    call MPI_BCAST(ES%R_lim,      1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(ES%Z_lim,      1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(ES%Psi_lim,    1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(ES%s_lim,      1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(ES%t_lim,      1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(ES%i_elm_lim,  1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(ES%ifail_lim,  1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    
    ! --- X-Point(s)
    call MPI_BCAST(ES%xpoint,         1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(ES%xcase,          1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(ES%active_xpoint,  1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(ES%R_xpoint,       2,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(ES%Z_xpoint,       2,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(ES%Z_xpoint_init,  2,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(ES%Psi_xpoint,     2,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(ES%s_xpoint,       2,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(ES%t_xpoint,       2,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(ES%i_elm_xpoint,   2,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(ES%ifail_xpoint,   1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    
    ! --- Boundary Point
    call MPI_BCAST(ES%psi_bnd,     1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(ES%R_bnd,       1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(ES%Z_bnd,       1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(ES%Psi_bnd,     1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(ES%Psi_bnd_init,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(ES%s_bnd,       1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(ES%t_bnd,       1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(ES%i_elm_bnd,   1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(ES%ifail_bnd,   1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr) 

    ! --- Strike Point(s) derived from axisymmetric field component.
    call MPI_BCAST(ES%num_strike,          1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(ES%i_bndelm_strike,    99,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(ES%R_strike,           99,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(ES%Z_strike,           99,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(ES%s_strike,           99,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    
    ! --- Inner/Outer points on the midplane close to the boundary of the computational domain.
    call MPI_BCAST(ES%R_midpl,           2,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

    ! --- LCFS shape parameters
    call MPI_BCAST(ES%LCFS_Rgeo   ,      1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(ES%LCFS_a      ,      1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(ES%LCFS_epsilon,      1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(ES%LCFS_kappa  ,      1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(ES%LCFS_deltaU ,      1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(ES%LCFS_deltaL ,      1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    
  end subroutine broadcast_equil_state
  





  !> Calculates the shape parameters of the LCFS
  !> as defined in T. Luce, PPCF 55 (2013) 095009, equations (1-6)
  subroutine LCFS_shape_parameters(node_list,element_list)
  
    use data_structure
    use phys_module, only: xcase, xpoint
    use mod_interp
    
    implicit none
    
    ! --- Input parameters.
    type (type_node_list),    intent(in)    :: node_list
    type (type_element_list), intent(in)    :: element_list
    
    ! --- Local variables
    integer :: n_int
    integer :: i_elm, j, k, n1, n2, n3
    real*8  :: t,rr1, rr2, drr1, drr2, ss1, ss2, dss1, dss2, ri, si, dri, dsi, dl
    real*8  :: RRgi, dRRgi_dr, dRRgi_ds, ZZgi, dZZgi_dr, dZZgi_ds, dRRgi_dt, dZZgi_dt
    real*8  :: PSgi, dPSgi_dr, dPSgi_ds, PSI_R, PSI_Z, RZJAC, grad_psi, psi_n
    real*8  :: dRRgi_drs,dRRgi_drr,dRRgi_dss, dZZgi_drs,dZZgi_drr,dZZgi_dss, dPSgi_drs,dPSgi_drr,dPSgi_dss
    real*8  :: Rmax, Z_Rmax, Rmin, Z_Rmin, R_Zmax, Zmax, R_Zmin, Zmin
    integer :: i,m, ig, ip, npoints
   
    type (type_surface_list) :: surface_list
    
    
    surface_list%n_psi = 1 
    allocate( surface_list%psi_values(surface_list%n_psi) )
    surface_list%psi_values(1) = ES%psi_axis + (ES%psi_bnd - ES%psi_axis) * 0.9999d0 ! Not 1 to avoid legs
    
    call find_flux_surfaces(99, xpoint, xcase, node_list, element_list, surface_list)
    
    Rmax = -1.d99;   Zmax = -1.d99
    Rmin =  1.d99;   Zmin =  1.d99


    npoints = 40
    
    do i=1, surface_list%n_psi
      do k=1, surface_list%flux_surfaces(i)%n_pieces
        do ig = 1, npoints
          t = -1.0 + 2.0*float(ig-1)/float(npoints-1)
    
          rr1  = surface_list%flux_surfaces(i)%s(1,k)
          drr1 = surface_list%flux_surfaces(i)%s(2,k)
          rr2  = surface_list%flux_surfaces(i)%s(3,k)
          drr2 = surface_list%flux_surfaces(i)%s(4,k)
    
          ss1  = surface_list%flux_surfaces(i)%t(1,k)
          dss1 = surface_list%flux_surfaces(i)%t(2,k)
          ss2  = surface_list%flux_surfaces(i)%t(3,k)
          dss2 = surface_list%flux_surfaces(i)%t(4,k)
    
          call CUB1D(rr1, drr1, rr2, drr2, t, ri, dri)
          call CUB1D(ss1, dss1, ss2, dss2, t, si, dsi)
    
          i_elm = surface_list%flux_surfaces(i)%elm(k)
    
          call interp(node_list,element_list,i_elm,1,1,ri,si,PSgi,dPSgi_dr,dPSgi_ds,dPSgi_drs,dPSgi_drr,dPSgi_dss)
    
          call interp_RZ(node_list,element_list,i_elm,ri,si,RRgi,dRRgi_dr,dRRgi_ds,dRRgi_drs,dRRgi_drr,dRRgi_dss, &
                                                            ZZgi,dZZgi_dr,dZZgi_ds,dZZgi_drs,dZZgi_drr,dZZgi_dss)
                                                            
          ! --- Ignore open and private field line regions
          if ( get_psi_n(PSgi, ZZgi) > 1.d0 ) cycle
    
          if (RRgi > Rmax) then
            Rmax   = RRgi;    Z_Rmax = ZZgi;
          endif 
    
          if (ZZgi > Zmax) then
            R_Zmax = RRgi;    Zmax = ZZgi;
          endif 
    
          if (RRgi < Rmin) then
            Rmin   = RRgi;    Z_Rmin = ZZgi;
          endif 
    
          if (ZZgi < Zmin) then
            R_Zmin = RRgi;    Zmin = ZZgi;
          endif 
    
        end do
      end do
    
      ! --- As defined in T. Luce, PPCF 55 (2013) 095009, equations (1-6)
      ES%LCFS_Rgeo    = (Rmax + Rmin) / 2.0 
      ES%LCFS_a       = (Rmax - Rmin) / 2.0
      ES%LCFS_epsilon =  ES%LCFS_a / ES%LCFS_Rgeo
      ES%LCFS_kappa   = (Zmax - Zmin) / (2.0 * ES%LCFS_a ) 
      ES%LCFS_deltaU  = (ES%LCFS_Rgeo - R_Zmax) / ES%LCFS_a
      ES%LCFS_deltaL  = (ES%LCFS_Rgeo - R_Zmin) / ES%LCFS_a
    
    end do
  
  end subroutine LCFS_shape_parameters
  
  
end module equil_info 
