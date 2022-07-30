!> This module determines the poloidal straight field line angle theta*.
module mod_straight_field_line
  
  
  
  
  
  use constants
  use tr_module 
  use mod_parameters, only: n_vertex_max, n_degrees, n_plane, n_tor, n_var
  use phys_module,     only: F0, xpoint, xcase
  use equil_info
  use mod_interp
  
  
  
  
  
  implicit none
  
  
  
  
  
  save
  
  
  
  
  
  private
  public determine_theta_mag, log_mapping, t_theta_mapping, cleanup_mapping, new_determine_theta_mag!###
  
  
  
  
  
  ! --- Constants
  character(len=23), parameter, private :: THIS_MOD_NAME = 'mod_straight_field_line'
  logical,           parameter, private :: DEBUG         = .true. !< Switch on/off debugging output
  
  
  
  
  
  !> Data structure describing the mapping between theta and theta_mag.
  type :: t_theta_mapping
    integer              :: nstpts   !< Number of radial positions
    real*8,  allocatable :: psin(:)  !< Psi_normalized values at the start points
    real*8,  allocatable :: rr(:,:)  !< R coordinates of field line positions at large steps
    real*8,  allocatable :: zz(:,:)  !< Z coordinates of field line positions at large steps
    real*8,  allocatable :: tt(:,:)  !< Geometrical theta at large steps
    real*8,  allocatable :: t2(:,:)  !< Magnetic theta at large steps
    integer, allocatable :: npts(:)  !< Number of large step positions for each start point
    real*8,  allocatable :: rre(:,:) !< R coordinates at equidistant theta_mag(!) positions
    real*8,  allocatable :: zze(:,:) !< Z coordinates at equidistant theta_mag(!) positions
  end type t_theta_mapping
  
  
  
  
  
  contains
  
  
  
  
  
  !> NOTHING USEFUL YET...
  subroutine new_determine_theta_mag(node_list, element_list, eq, PsiNRange, nPsiN, nTht, &
    ierr)
    
    use data_structure
    
    ! --- Constants
    character(len=64), parameter :: THIS_ROUTINE_NAME = trim(THIS_MOD_NAME) // ':NEW_determine_theta_mag' !###
    
    ! --- Routine parameters
!    type(t_theta_mapping),  intent(inout) :: mapping
    type(type_node_list),   intent(in)    :: node_list
    type(type_element_list),intent(in)    :: element_list
    type(t_equil_state),    intent(in)    :: eq          !< Equilibrium state information
    real*8,                 intent(in)    :: PsiNRange(2)!< Radial range of Psi_N values to cover
    integer,                intent(in)    :: nPsiN       !< Number of flux surfaces
    integer,                intent(in)    :: nTht        !< Number of points in theta*
    integer,                intent(inout) :: ierr        !< Error code
    
    ! --- Local variables
    type(type_surface_list) :: surface_list
    integer :: i, i_elm, j, ip, nplot=5
    real*8  :: psi_n, R, R_s, R_t, R_st, R_ss, R_tt, Z, Z_s, Z_t, Z_st, Z_ss, Z_tt
    real*8  :: u, ss1, dss1, ss2, dss2, si, dsi, tt1, dtt1, tt2, dtt2, ti, dti
    
    surface_list%n_psi = nPsiN
    allocate( surface_list%psi_values   (nPsiN) )
    
    do i = 1, nPsiN
      psi_n = PsiNRange(1) + (PsiNRange(2)-PsiNRange(1)) * real(i-1)/real(nPsiN-1)
      surface_list%psi_values(i) = eq%psi_axis + ( eq%psi_bnd - eq%psi_axis ) * psi_n
    end do
    
    call find_flux_surfaces(0,eq%xpoint, eq%xcase, node_list, element_list, surface_list)
    
    !###
    ! --- Loop over all flux surfaces
    do i = 1, surface_list%n_psi
      
      ! --- Loop over all segments of this flux surface
      do j = 1, surface_list%flux_surfaces(i)%n_pieces
        
        ! --- Bezier element, in which the current flux surface segment is located
        i_elm = surface_list%flux_surfaces(i)%elm(j)
        
        ss1  = surface_list%flux_surfaces(i)%s(1,j)
        dss1 = surface_list%flux_surfaces(i)%s(2,j)
        ss2  = surface_list%flux_surfaces(i)%s(3,j)
        dss2 = surface_list%flux_surfaces(i)%s(4,j)
        
        tt1  = surface_list%flux_surfaces(i)%t(1,j)
        dtt1 = surface_list%flux_surfaces(i)%t(2,j)
        tt2  = surface_list%flux_surfaces(i)%t(3,j)
        dtt2 = surface_list%flux_surfaces(i)%t(4,j)
        
        ! --- Loop over nplot points in a flux surface segment
        do ip = 1, nplot
          u = -1. + 2.*float(ip-1)/float(nplot-1)
          
          ! --- Determine s and t values of the current point inside element i_elm
          call CUB1D(ss1, dss1, ss2, dss2, u, si, dsi)
          call CUB1D(tt1, dtt1, tt2, dtt2, u, ti, dti)
          
          ! --- Determine (R,Z)-coordinates of the current point on the current flux surface
          call interp_RZ(node_list, element_list, i_elm, si, ti, R, R_s, R_t, R_st, R_ss, R_tt, &
            Z, Z_s, Z_t, Z_st, Z_ss, Z_tt)
            
          ! --- Write out the (R,Z)-coordinates
          write(61,'(2ES16.7)') R, Z
        end do
        
        write(61,*)
        write(61,*)
        
      end do
      
    end do
    !###
    
  end subroutine new_determine_theta_mag
  
  
  
  
  
  !> Clean up mapping data structure.
  subroutine cleanup_mapping(mapping)
    
    ! --- Routine parameters
    type(t_theta_mapping),  intent(inout) :: mapping
    
    mapping%nstpts = 0
    if ( allocated(mapping%psin) ) deallocate(mapping%psin)
    if ( allocated(mapping%rr  ) ) deallocate(mapping%rr  )
    if ( allocated(mapping%zz  ) ) deallocate(mapping%zz  )
    if ( allocated(mapping%tt  ) ) deallocate(mapping%tt  )
    if ( allocated(mapping%t2  ) ) deallocate(mapping%t2  )
    if ( allocated(mapping%npts) ) deallocate(mapping%npts)
    if ( allocated(mapping%rre ) ) deallocate(mapping%rre )
    if ( allocated(mapping%zze ) ) deallocate(mapping%zze )
    
  end subroutine cleanup_mapping
  
  
  
  
  
  !> Allocate mapping data structure.
  subroutine alloc_mapping(mapping, nPsiN, nmaxsteps, nTht)
    
    ! --- Routine parameters
    type(t_theta_mapping),  intent(inout) :: mapping
    integer,                intent(in)    :: nPsiN       !< Number of flux surfaces
    integer,                intent(in)    :: nmaxsteps   !< Maximum steps on one flux surface
    integer,                intent(in)    :: nTht        !< Number of points in theta*
    
    call cleanup_mapping(mapping)
    
    mapping%nstpts = nPsiN
    
    allocate( mapping%rr(nPsiN,0:nmaxsteps) )
    allocate( mapping%zz(nPsiN,0:nmaxsteps) )
    allocate( mapping%tt(nPsiN,0:nmaxsteps) )
    allocate( mapping%t2(nPsiN,0:nmaxsteps) )
    allocate( mapping%npts(nPsiN)           )
    allocate( mapping%psin(nPsiN)           )
    allocate( mapping%rre(nPsiN,0:nTht-1)   )
    allocate( mapping%zze(nPsiN,0:nTht-1)   )
    
    mapping%rr   = 0.
    mapping%zz   = 0.
    mapping%tt   = 0.
    mapping%t2   = 0.
    mapping%npts = 0
    mapping%psin = 0.
    mapping%rre  = 0.
    mapping%zze  = 0.
    
  end subroutine alloc_mapping
  
  
  
  
  
  !> Determine the magnetic poloidal angle, theta_mag, by performing field line tracing in the
  !! axisymmetric part of the magnetic field.
  !!
  !! The optional parameters (nmaxsteps2, deltaphi2, nsmallsteps2) are usually not required but can
  !! help to optimize the runtime. In some cases, deltaphi2 needs to be set manually to obtain
  !! accurate results. In some cases, nmaxsteps2 might have to be increased in order to complete
  !! a full poloidal turn on a flux surface close to the separatrix.
  subroutine determine_theta_mag(mapping, node_list, element_list, equil_state, PsiNRange, nPsiN,  &
    nTht, ierr, nmaxsteps2, deltaphi2, nsmallsteps2)
    
    ! --- Constants
    character(len=64), parameter :: THIS_ROUTINE_NAME = trim(THIS_MOD_NAME) // ':determine_theta_mag'
    
    ! --- Routine parameters
    type(t_theta_mapping),  intent(inout) :: mapping
    type(type_node_list),   intent(in)    :: node_list
    type(type_element_list),intent(in)    :: element_list
    type(t_equil_state),    intent(in)    :: equil_state !< Equilibrium state information
    real*8,                 intent(in)    :: PsiNRange(2)!< Radial range of Psi_N values to cover
    integer,                intent(in)    :: nPsiN       !< Number of flux surfaces
    integer,                intent(in)    :: nTht        !< Number of points in theta*
    integer,                intent(inout) :: ierr        !< Error code
    integer, optional,      intent(in)    :: nmaxsteps2  !< Maximum number of poloidal steps
    real*8,  optional,      intent(in)    :: deltaphi2   !< Toroidal step width of the large steps
    integer, optional,      intent(in)    :: nsmallsteps2!< Number of small steps between large ones
    
    ! --- Local variables
    integer :: i, j, k, nmaxsteps, nsmallsteps
    real*8  :: deltaphi, test1, test2, test3, suggested_factor
    
    ! --- Initialisation
    suggested_factor = 1.d0
    ierr             = 0
    
    deltaphi         = min(0.3d0, 20.d0 / real(nTht)) !###
    nmaxsteps        = 400 / deltaphi !###
    nsmallsteps      = 3 !###
    if ( present(nmaxsteps2)   ) nmaxsteps   = nmaxsteps2
    if ( present(deltaphi2)    ) deltaphi    = deltaphi2
    if ( present(nsmallsteps2) ) nsmallsteps = nsmallsteps2
    if ( DEBUG ) then
      write(*,*) 'nmaxsteps   =', nmaxsteps
      write(*,*) 'deltaphi    =', deltaphi
      write(*,*) 'nsmallsteps =', nsmallsteps
      write(*,*) 'psiNrange   =', psiNrange

    end if
    
    call alloc_mapping(mapping, nPsiN, nmaxsteps, nTht)
    
    ! --- Trace field lines in the axisymmetric magnetic field component
    call trace_fieldlines(mapping, node_list, element_list, equil_state, PsiNRange, nmaxsteps,     &
      deltaphi, nsmallsteps, ierr)
    if ( ierr /= 0 ) then
      write(*,*) 'ERROR in '//trim(THIS_ROUTINE_NAME)//' calling trace_fieldlines.'
      return
    end if
    
    ! --- Interpolate theta_mag to equidistant theta_geo positions.
    do k = 1, mapping%nstpts
      if ( mapping%npts(k) < nTht ) then
        suggested_factor = max( suggested_factor, real(nTht)/real(mapping%npts(k)) )
      end if
      call interpolEquidist(mapping%t2(k,0:mapping%npts(k)), mapping%rr(k,0:mapping%npts(k)),      &
        mapping%npts(k)+1, mapping%rre(k,0:nTht-1), mapping%t2(k,0),                      &
        mapping%t2(k,0)+2.*PI*real(nTht-1)/real(nTht), nTht)
      call interpolEquidist(mapping%t2(k,0:mapping%npts(k)), mapping%zz(k,0:mapping%npts(k)),      &
        mapping%npts(k)+1, mapping%zze(k,0:nTht-1), mapping%t2(k,0),                      &
        mapping%t2(k,0)+2.*PI*real(nTht-1)/real(nTht), nTht)
    end do
    
    if ( suggested_factor > 1.d0 ) then
      write(*,*)
      write(*,*) '*** WARNING ***************************************************************'
      write(*,*) 'When calling '//trim(THIS_ROUTINE_NAME)//', deltaphi should be decreased'
      write(*,*) 'Current value for deltaphi:      deltaphi = ', deltaphi
      write(*,*) 'Recommended value for this case: deltaphi < ', deltaphi / suggested_factor
      write(*,*) '***************************************************************************'
      write(*,*)
    end if
    
    ! --- Output some data to files for debugging.
    if ( debug ) then
      open(42, FILE='determine_theta_mag.rr-zz.dat', ACTION='WRITE', STATUS='REPLACE')
      do k = 1, mapping%nstpts
        do j = 0, mapping%npts(k)
          write(42,*) mapping%rr(k,j), mapping%zz(k,j)
        end do
        write(42,*)
      end do
      close(42)
      
      open(42, FILE='determine_theta_mag.tt-t2.dat', ACTION='WRITE', STATUS='REPLACE')
      do k = 1, mapping%nstpts
        do j = 0, mapping%npts(k)
          write(42,*) mapping%tt(k,j), mapping%t2(k,j)
        end do
        write(42,*)
      end do
      close(42)
      
      open(42, FILE='determine_theta_mag.rre-zze.dat', ACTION='WRITE', STATUS='REPLACE')
      do k = 1, mapping%nstpts
        do j = 0, nTht
          write(42,*) mapping%rre(k,mod(j,nTht)), mapping%zze(k,mod(j,nTht))
        end do
        write(42,*)
      end do
      close(42)
  
      open(42, FILE='determine_theta_mag.rre-zze.2.dat', ACTION='WRITE', STATUS='REPLACE')
      do j = 0, nTht-1
        do k = 1, mapping%nstpts
          write(42,*) mapping%rre(k,j), mapping%zze(k,j)
        end do
        write(42,*)
      end do
      close(42)

      call log_mapping(mapping)
    end if
    
  end subroutine determine_theta_mag
  
  
  
  
  
  !> Trace field lines in the n=0 component of Psi.
  subroutine trace_fieldlines(mapping, node_list, element_list, equil_state, PsiNRange,     &
    nmaxsteps, deltaphi, nsmallsteps, ierr)
    use mod_interp
    
    ! --- Constants
    character(len=64), parameter :: THIS_ROUTINE_NAME = trim(THIS_MOD_NAME) // ':trace_fieldlines'
    
    ! --- Routine parameters
    type(t_theta_mapping),  intent(inout) :: mapping
    type(type_node_list),   intent(in)    :: node_list
    type(type_element_list),intent(in)    :: element_list
    type(t_equil_state),    intent(in)    :: equil_state !< Equilibrium state information
    real*8,                 intent(in)    :: PsiNRange(2)!< Radial range of Psi_N values to cover
    integer,                intent(in)    :: nmaxsteps   !< Maximum number of poloidal steps
    real*8,                 intent(in)    :: deltaphi    !< Toroidal step width of the large steps
    integer,                intent(in)    :: nsmallsteps !< Number of small steps between large ones
    integer,                intent(inout) :: ierr        !< Return code (0 if no errors)
    
    ! --- Local variables
    integer :: i, j, k, l
    real*8  :: rn, zn ! R and Z at previous small step position
    real*8  :: rp, zp ! R and Z at next small step position
    real*8  :: rh, zh ! R and Z at half small step position
    real*8  :: dpsi_dzn, dpsi_drn ! derivatives of psi at previous small step
    real*8  :: dpsi_dzh, dpsi_drh ! derivatives of psi at half small step
    real*8  :: R_out, Z_out, s_out, t_out
    integer :: i_elm_out, ifail,my_id
    real*8 :: P, P_s, P_t, P_st, P_ss, P_tt
    real*8 :: x, x_s, x_t, y, y_s, y_t
    real*8 :: xjac
    real*8 :: theta_corr
    real*8 :: smalldeltaphi
    real*8 :: phi, fact, Rleft, Rright, Rmid
    
    smalldeltaphi = deltaphi / nsmallsteps
    my_id = 0
    ierr  = 0
    
    !$omp parallel do schedule(dynamic) default(none)                                              &
    !$omp   firstprivate(theta_corr, phi, Rleft, Rright, Rmid, i_elm_out, s_out, t_out, P, P_s, &
    !$omp     P_t, rn, zn,i, R_out, Z_out, ifail, x, x_s, x_t, y, y_s, y_t, xjac, dpsi_dzn,        &
    !$omp     dpsi_drn, rh, zh, rp, zp, l, dpsi_dzh, dpsi_drh)                                     &
    !$omp   shared(mapping, node_list, element_list, equil_state, ierr, F0, smalldeltaphi,         &
    !$omp   psinrange, nmaxsteps, nsmallsteps)                                                     &
    !$omp   private(k)
    FL_STPTS: do k = mapping%nstpts, 1, -1 ! (inverse order better for OpenMP parallelization)
      if ( ierr /= 0 ) cycle
      
      theta_corr = 0.d0
      phi        = 0.d0
      
      ! --- Find field line starting position via bisections.
      if ( mapping%nstpts == 1 ) then
        mapping%psin(k) = PsiNRange(1)
      else
        mapping%psin(k) = PsiNRange(1) + real(k-1)/ real(mapping%nstpts-1) *                       &
          ( PsiNRange(2) - PsiNRange(1) )
      end if
      Rleft  = equil_state%R_axis
      Rright = equil_state%R_midpl(2)
      do
        Rmid = (Rleft + Rright) / 2.d0
        
        if ( Rright - Rleft < 1.d-7 ) then ! If converged, ...
          mapping%rr(k,0) = Rmid
          mapping%zz(k,0) = equil_state%Z_axis
          mapping%tt(k,0) = 0.d0
          exit
        end if
        
        call find_RZ(node_list, element_list, Rmid, equil_state%Z_axis, R_out, Z_out, i_elm_out,   &
          s_out, t_out, ifail)
        call interp0(node_list, element_list, i_elm_out, 1, 1, s_out, t_out, P, P_s, P_t)
        if ( get_psi_n(P, Z_out)  < mapping%psin(k) ) then
          Rleft  = Rmid
        else
          Rright = Rmid
        end if
      end do
      
      !### direction of field line tracing!
      
      ! --- Perform the field line tracing.
      FL_LARGESTEPS: do j = 1, NMAXSTEPS
        if ( ierr /= 0 ) cycle
        
        rn = mapping%rr(k,j-1)
        zn = mapping%zz(k,j-1)
        
        FL_SMALLSTEPS: do i = 1, NSMALLSTEPS
          if ( ierr /= 0 ) cycle
          
          ! --- Predictor step
          ! - Determine element number and s and t coordinates for given (R, Z) position.
          call find_RZ(node_list,element_list,rn,zn,R_out,Z_out,i_elm_out,s_out,t_out,ifail)
          if ( ifail /= 0 ) then
            write(*,*) 'ERROR in '//trim(THIS_ROUTINE_NAME)//' calling find_RZ (1).'
            ierr = 100
	    cycle
          end if
          ! - Interpolate Psi to current position.
          call interp0(node_list, element_list, i_elm_out,1,1,s_out,t_out,P,P_s,P_t)
          ! - Determine derivatives of R and Z with respect to s and t at current position.
          call interp_RZ(node_list,element_list,i_elm_out,s_out,t_out,x,x_s,x_t,y,y_s,y_t)
          ! - Determine derivatives of Psi with respect to R and Z.
          xjac     = x_s*y_t - x_t*y_s
          dpsi_dzn = ( - x_t*P_s + x_s*P_t ) / xjac
          dpsi_drn = ( + y_t*P_s - y_s*P_t ) / xjac
          ! - Advance a half step.
          rh = rn - dpsi_dzn/F0*rn * smalldeltaphi/2.
          zh = zn + dpsi_drn/F0*rn * smalldeltaphi/2.
          
          ! --- Corrector step
          ! - Determine element number and s and t coordinates for given (R, Z) position.
          call find_RZ(node_list,element_list,rh,zh,R_out,Z_out,i_elm_out,s_out,t_out,ifail)
          if ( ifail /= 0 ) then
            write(*,*) 'Error in '//trim(THIS_ROUTINE_NAME)//' calling find_RZ (2).'
            ierr = 100
	    cycle
          end if
          ! - Interpolate Psi to current position.
          call interp0(node_list, element_list, i_elm_out,1,1,s_out,t_out,P,P_s,P_t)
          ! - Determine derivatives of R and Z with respect to s and t at current position.
          call interp_RZ(node_list,element_list,i_elm_out,s_out,t_out,x,x_s,x_t,y,y_s,y_t)
          ! - Determine derivatives of Psi with respect to R and Z.
          xjac     = x_s*y_t - x_t*y_s
          dpsi_dzh = ( - x_t*P_s + x_s*P_t ) / xjac
          dpsi_drh = ( + y_t*P_s - y_s*P_t ) / xjac
          ! - Advance a full step.
          rp = rn - dpsi_dzh/F0*rh * smalldeltaphi
          zp = zn + dpsi_drh/F0*rh * smalldeltaphi
          
          phi = phi + smalldeltaphi
          
        end do FL_SMALLSTEPS
        
        ! --- Write position after series of small steps into data structure.
        mapping%rr(k,j) = rp
        mapping%zz(k,j) = zp
        mapping%tt(k,j) = atan3( zp-equil_state%Z_axis, rp-equil_state%R_axis )
        
        ! --- Compensate for theta jumps by 2*PI.
        if ( (mapping%tt(k,j-1)-theta_corr > 3.*PI/2.) .and. (mapping%tt(k,j) < PI/2.) )           &
	  theta_corr = theta_corr + 2.*PI
        if ( (mapping%tt(k,j) > 3.*PI/2.) .and. (mapping%tt(k,j-1)-theta_corr < PI/2.) )           &
	  theta_corr = theta_corr - 2.*PI
        mapping%tt(k,j) = mapping%tt(k,j) + theta_corr
        
        ! --- Stop following field line k as soon as one poloidal turn is completed.
        if ( abs(mapping%tt(k,j) - mapping%tt(k,0)) > 2*PI ) then
          mapping%npts(k) = j
          do l = 0, mapping%npts(k)
            mapping%t2(k,l) = real(l) / ( mapping%npts(k) - &
              ( ( mapping%tt(k,mapping%npts(k)) - mapping%tt(k,0) - sign(2.*PI,mapping%tt(k,j)-mapping%tt(k,0)) ) / &
              ( mapping%tt(k,mapping%npts(k)) - mapping%tt(k,mapping%npts(k)-1) ) ) ) * 2.*PI
          end do
          if ( debug ) then 
            !$omp critical
            write(*,'(" Field line",I6,":    psin=",F7.3,"    npts=",I7)') k, mapping%psin(k),     &
              mapping%npts(k)
            !$omp end critical
          end if
          exit
        end if
        
        if ( j == NMAXSTEPS ) then
          write(*,*) 'ERROR in '//trim(THIS_ROUTINE_NAME)//': nmaxsteps too small, incomplete poloidal turn.'
          ierr = 101
          cycle
        end if
        
      end do FL_LARGESTEPS
      
    end do FL_STPTS
    !$omp end parallel do
    
    if ( ierr /= 0 ) then
      write(*,*) 'Aborting '//trim(THIS_ROUTINE_NAME)//' after an error occurred.'
      stop
    end if
    
    ! --- Make theta positive.
    do k = 1, mapping%nstpts
      mapping%tt(k,0:mapping%npts(k)) = mapping%tt(k,0:mapping%npts(k))                            &
        - minval( mapping%tt(k,0:mapping%npts(k)) )
    end do
    
  end subroutine trace_fieldlines
  
  
  
  
  
  !> Interpolate y(x) data to equidistant x-positions.
  subroutine interpolEquidist(xOld, yOld, nOld, yNew, xMin, xMax, nNew)
    
    integer, intent(in)    :: nOld         !< Number of old positions
    real*8,  intent(in)    :: xOld(nOld)   !< old x values
    real*8,  intent(in)    :: yOld(nOld)   !< old y values
    integer, intent(in)    :: nNew         !< Number of new positions
    real*8,  intent(inout) :: yNew(nNew)   !< New y values
    real*8,  intent(in)    :: xMin         !< Smallest new x value
    real*8,  intent(in)    :: xMax         !< Largest new x value
    
    integer :: k, j
    real*8  :: xNew
    
    j = 1
    LOOP_OLDPTS: do k = 2, nOld
      LOOP_NEWPTS: do
        xNew = xMin + (xMax - xMin) * REAL(j-1)/REAL(nNew-1)
        if ( xNew <= xOld(k) ) then
          yNew(j) = yOld(k-1) + ( xNew - xOld(k-1) ) / ( xOld(k) - xOld(k-1) )                     &
	    * ( yOld(k) - yOld(k-1) )
          j = j + 1
        else
          cycle LOOP_OLDPTS
        end if
        if ( j > nNew ) exit LOOP_OLDPTS
      end do LOOP_NEWPTS
    end do LOOP_OLDPTS
    
  end subroutine interpolEquidist
  
  
  
  
  
  !> Same as interp, but less derivatives are calculated
  pure subroutine interp0(node_list, element_list, i_elm, i_var, i_harm, s, t, P, P_s, P_t)
    use mod_basisfunctions
  
    type(type_node_list),   intent(in)    :: node_list
    type(type_element_list),intent(in)    :: element_list
    integer, intent(in)                   :: i_elm, i_var, i_harm
    real*8, intent(in)                    :: s, t
    real*8, intent(out)                   :: P, P_s, P_t
    real*8 :: G(4,n_degrees), G_s(4,n_degrees), G_t(4,n_degrees) 
    real*8 :: G_st(4,n_degrees), G_ss(4,n_degrees), G_tt(4,n_degrees)
    integer :: kv, iv, kf
    
    call basisfunctions(s,t,G(1:4,1:n_degrees),G_s(1:4,1:n_degrees),G_t(1:4,1:n_degrees), &
      G_st(1:4,1:n_degrees),G_ss(1:4,1:n_degrees),G_tt(1:4,1:n_degrees))
  
    P = 0.d0; P_s = 0.d0; P_t = 0.d0
  
    do kv = 1,n_vertex_max  ! 4 vertices
    
      iv = element_list%element(i_elm)%vertex(kv)  ! the node number
    
      do kf = 1, n_degrees       ! basis functions
    
        P   = P   + node_list%node(iv)%values(i_harm,kf,i_var)                                     &
	  * element_list%element(i_elm)%size(kv,kf) * G(kv,kf)
        P_s = P_s + node_list%node(iv)%values(i_harm,kf,i_var)                                     &
	  * element_list%element(i_elm)%size(kv,kf) * G_s(kv,kf)
        P_t = P_t + node_list%node(iv)%values(i_harm,kf,i_var)                                     &
	  * element_list%element(i_elm)%size(kv,kf) * G_t(kv,kf)
    
      end do
    
    end do
    
  end subroutine interp0
  
  
  
  
  
  
  !> Same as atan2 (Fortran built-in), but returns only positive values by adding 2*pi if necessary
  real*8 function atan3( dy, dx )
    real*8, intent(in) :: dy, dx
    atan3 = ATAN2(dy,dx)
    if ( atan3 < 0 ) atan3 = atan3 + 2.*PI
  end function atan3
  
  
  
  
  
  !> Write information about the mapping determined by determine_theta_mag to the logfile.
  subroutine log_mapping(mapping)
    
    type(t_theta_mapping), intent(in) :: mapping
    
    write(*,'(A,I12,A)')      'nstpts     =', mapping%nstpts,	     ' (number of radial positions)'
    
    write(*,'(A)', ADVANCE='NO') 'psin       ='
    if ( allocated(mapping%psin) ) then
      write(*,'(99F7.3)', ADVANCE='NO') mapping%psin
    else
      write(*,'(A)', ADVANCE='NO') ' not allocated'
    end if
    write(*,'(A)') ' (norm. pol. flux at rad. positions)'
    
    write(*,'(A)', ADVANCE='NO') 'npts       ='
    if ( allocated(mapping%npts) ) then
      write(*,'(99I7)', ADVANCE='NO') mapping%npts
    else
      write(*,'(A)', ADVANCE='NO') ' not allocated'
    end if
    write(*,'(A)') ' (num. of points from field line tracing)'
    
  end subroutine log_mapping
  
  
  
  
  
end module mod_straight_field_line
