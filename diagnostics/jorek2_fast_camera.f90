!> Program to convert a JOREK2 restart file into binary VTK format
program jorek2_fast_camera

  use mod_parameters, only: n_var, variable_names, jorek_model
  use mod_import_restart
  use data_structure
  use phys_module
  use basis_at_gaussian
  use diffusivities, only: get_dperp, get_zkperp
  use mod_interp
  use constants, only: PI, LOWER_XPOINT, UPPER_XPOINT
  use equil_info
  use mod_element_rtree, only: populate_element_rtree

  implicit none
  !include 'mpif.h'

  type (type_node_list)   , pointer :: node_list
  type (type_element_list), pointer :: element_list
  type (type_surface_list)          :: flux_list

  integer               :: nnoel, nnos, nel, nsub, inode, ielm, n_scalars, n_vectors
  real*4,allocatable    :: xyz (:,:), scalars(:,:), vectors(:,:,:)
  integer,allocatable   :: ien (:,:)
  integer, parameter    :: ivtk = 22 ! an arbitrary unit number for the VTK output file
  integer               :: i, j, k, m, etype, irst, int, i_var, i_tor, i_tor_old, i_plane, index, index_node
  character             :: buffer*80, lf*1, str1*12, str2*12
  real*8                :: s, t
  real*8                :: P,P_s,P_t,P_st,P_ss,P_tt
  real*8                :: R,R_s,R_t,R_st,R_ss,R_tt
  real*8                :: Z,Z_s,Z_t,Z_st,Z_ss,Z_tt
  real*8                :: PPPsi,Ps_s,Ps_t,Ps_st,Ps_ss,Ps_tt
  real*8                :: ZJ,ZJ_s,ZJ_t,ZJ_st,ZJ_ss,ZJ_tt
  real*8                :: U,U_s,U_t,U_st,U_ss,U_tt
  real*8                :: W,W_s,W_t,W_st,W_ss,W_tt
  real*8                :: RRRHO,RHO_s,RHO_t,RHO_st,RHO_ss,RHO_tt
  real*8                :: TTT,TT_s,TT_t,TT_st,TT_ss,TT_tt
  real*8                :: Ti,Ti_s,Ti_t,Ti_st,Ti_ss,Ti_tt
  real*8                :: TTTe,Te_s,Te_t,Te_st,Te_ss,Te_tt
  real*8                :: V, V_s, V_t, V_st, V_ss, V_tt
  real*8                :: psi_00, rho_00, Ti_00, Te_00
  real*8                :: ps_x, ps_y
  real*8                :: u0_x, u0_y
  real*8                :: zj_x, zj_y
  real*8                :: w0_x, w0_y, w0_xx, w0_yy
  real*8                :: RHO_x, RHO_y, RHO_p
  real*8                :: TT_x, TT_y, TT_p
  real*8                :: Ti_x, Ti_y, Ti_p
  real*8                :: Te_x, Te_y, Te_p
  real*8                :: ps0, psi_norm, grad_psi
  real*8                :: xjac, xjac_x, xjac_y, v_perp, Psi_J, R_p, error, Btot, BigR
  real*8                :: particle_source, D_prof, ZK_prof, source_pellet, ZKpar_T
  integer               :: i_find, i_elm_find(8)
  real*8                :: Router,dRRg1_dr,dRRg1_ds,dRRg1_drs,dRRg1_drr,dRRg1_dss
  real*8                :: Zouter,dZZg1_dr,dZZg1_ds,dZZg1_drs,dZZg1_drr,dZZg1_dss
  real*8                :: s_find(8), t_find(8)
  real*8                :: Jb
  real*8                :: central_ne
  integer               :: k_tor
  


  ! --- MPI variables
  integer               :: my_id, n_cpu, ierr, nsend, nrecv
  integer               :: status(MPI_STATUS_SIZE)
  integer               :: pix_start, pix_end, pix_delta
  ! --- Photon Emissivity Coeff (PEC) variables
  integer               :: PEC_size, PEC_index_Ne, PEC_index_Te, PEC_index
  character*44          :: PEC_file
  character*50          :: line
  real*8,allocatable    :: PEC_dens(:), PEC_temp(:), PEC(:)
  real*8                :: PEC_tmp, Te_PEC_min
  ! --- Camera variables
  integer               :: n_pixels_hor, n_pixels_ver, n_pix, ncount
  real*8                :: pixel_dim, focus, vec_size
  real*8                :: angle_hor, angle_ver, toroidal_angle_location
  real*8                :: X_cam, Y_cam, Z_cam
  real*8,allocatable    :: Xp(:,:), Yp(:,:), Zp(:,:)
  real*8,allocatable    :: Xv(:,:), Yv(:,:), Zv(:,:)
  real*8,allocatable    :: Xref(:,:), Yref(:,:), Zref(:,:)
  real*8,allocatable    :: step_pix(:,:)
  integer,allocatable   :: n_step_pix(:,:)
  real*8,allocatable    :: Light(:)
  real*8, dimension(2)  :: Rl,Zl
  character*1           :: rgb(3)
  integer               :: itmp, icnt
  integer               :: colormap
  ! --- Integration variables
  integer               :: i_elm, ifail, inside_core, inside_wall, inside_tmp, inside_tmp2, i_int
  integer               :: n_step
  real*8                :: step, step_ref
  real*8                :: X_tmp, Y_tmp, Z_tmp
  real*8                :: X_prev,Y_prev,Z_prev
  real*8                :: RR,      ZZ,    Phi
  real*8                :: RR_prev, ZZ_prev
  real*8                :: RR_tmp, ZZ_tmp
  real*8                :: R_out, Z_out
  real*8                :: R_int, Z_int
  real*8                :: ss,    tt
  real*8                :: distance, distance_seg
  real*8,allocatable    :: HZ_tor(:)
  real*8                :: psi, rho, rho_n, Te, eV2Joules, solenoid
  logical               :: go_faster_in_core
  ! --- Reflection variables
  logical               :: include_reflections
  integer               :: n_reflections, i_ref
  real*8,allocatable    :: refl_coef(:)
  real*8                :: reflective_coefficient
  real*8, dimension(2)  :: Rvec,  Zvec
  real*8, dimension(2)  :: XXtmp, YYtmp, ZZtmp
  real*8, dimension(4)  :: XXvec, YYvec, ZZvec
  ! --- Neutrals variables
  real*8, dimension(10) :: rhon_prof
  real*8                :: tanh_psi, tanh_zmin, tanh_zpls
  ! --- Artificial diagnostic light (user defined)
  logical               :: artificial_light
  real*8                :: user_light

  real*8                :: maxlight, light_saturation
  
  ! --- User input parameters
  namelist /user_params/ n_pixels_hor, n_pixels_ver, pixel_dim, &
                         X_cam, Y_cam, Z_cam, &
                         toroidal_angle_location, angle_hor, angle_ver, &
                         focus, light_saturation, colormap, &
                         step, go_faster_in_core, &
                         include_reflections, reflective_coefficient, &
                         solenoid, rhon_prof, PEC_file, Te_PEC_min, &
                         artificial_light
  


  ! ********************************************************************** !
  !                              Start MPI                                 !
  ! ********************************************************************** !  
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)      ! id of each MPI proc
  call MPI_COMM_SIZE(MPI_COMM_WORLD, n_cpu, ierr)      ! number of MPI procs



  if (my_id .eq. 0) write(*,*) '/**************************************/'
  if (my_id .eq. 0) write(*,*) '/*******  jorek2_fast_camera  *********/'
  if (my_id .eq. 0) write(*,*) '/**************************************/'
  allocate(node_list)
  allocate(element_list)


  ! ********************************************************************** !
  !  User input parameter: Camera position, view, focal length, resolution !
  ! ********************************************************************** !
  
  ! --- MAST fast-visible camera midplane view
  ! --- Number of horizontal and vertical pixels
  n_pixels_hor = 800!400!800
  n_pixels_ver = 650!325!650
  ! --- Pixel dimension (assumed to be a square)
  pixel_dim    = 17.d-6!34.d-6!17.d-6
  
  ! --- Location of camera focus
  X_cam = 2.05
  Y_cam = 0.0 
  Z_cam = 0.0 

  ! --- Toroidal angular position of camera (easier so that Z_cam can always be zero)
  toroidal_angle_location = 0.d0 ! 0.5 * PI

  ! --- Orientation of camera (horizontal and vertical angles)
  angle_hor = 0.d0!-0.25 * PI
  angle_ver = 0.d0!-0.25 * PI
  
  ! --- Focal length of camera lens
  focus = 4.8d-3
 
  ! --- Artificial saturation of the camera: 1.0 = no saturation, 10.0 = 10x lower saturation, 0.1 = 10x higher saturation 
  light_saturation = 1.d0

  ! --- Colormap of output image
  ! --- =0  : From black to white (default, ie. natural BnW)
  ! --- =-1 : From white to black (reversed photons...)
  ! --- =1  : Heat colorbar from white to yellow to red to black
  ! --- =2  : Rainbow colorbar from blue to green to white to yellow to red to black
  colormap = 0
 
  ! --- The step size for integration and number of steps to get to the other side of plasma
  step   = 1.d-3 ! 1.d-2

  ! --- Increase integration step in the core (eg. on MAST, most light comes from the edge)
  go_faster_in_core = .false. 

  ! --- Include reflections on first-wall contour
  include_reflections = .true.
  ! --- Reflective coefficient of wall [0.0, 1.0]
  reflective_coefficient = 0.3

  ! --- Do we stop on the solenoid? (=0 if not)
  solenoid = 0.d0 !0.195

  ! --- Artificial neutrals density (if not using neutrals model500)
  rhon_prof(1) = 1.d-3  ! Core value [10^20 m^(-3)]
  rhon_prof(2) = 1.d-1  ! SOL value [10^20 m^(-3)]
  rhon_prof(3) = 3.d-1  ! Lower-divertor value [10^20 m^(-3)]
  rhon_prof(4) = 3.d-1  ! Upper-divertor value [10^20 m^(-3)]
  rhon_prof(5) = 0.10   ! Edge-psi tanh width [psi_norm]
  rhon_prof(6) = 0.95   ! Edge-psi tanh position [psi_norm]
  rhon_prof(7) = 0.20   ! Lower-divertor tanh width [m]
  rhon_prof(8) = -0.05  ! Lower-divertor tanh position relative to Z_xpoint(1) [m]
  rhon_prof(9) = 0.20   ! Upper-divertor tanh width [m]
  rhon_prof(10)= +0.05  ! Upper-divertor tanh position relative to Z_xpoint(2) [m]

  ! --- Location of emission data file
  PEC_file = "./my_pec.dat"
  ! --- Minimal temperature allowed for emission.
  ! --- This is very important for strongly non-linear cases with Te <= 0
  ! --- Because the emission spikes at very low Te.
  ! --- For neutrals simulations, where Te can be very low, maybe this can be lowered.
  Te_PEC_min = 20.0 !eV

  ! --- Artificial light instead of emission (user defined, hard-coded! default is edge gaussian)
  artificial_light = .false.





  ! --- Read parameters from namelist file 'fast_cam.nml' will overwrite all the above
  open(42, file='fast_cam.nml', action='read', status='old', iostat=ierr)
  if ( ierr == 0 ) then
    if (my_id .eq. 0) write(*,*) 'Reading parameters from fast_cam.nml namelist.'
    read(42,user_params)
    close(42)
  else
    if (my_id .eq. 0) write(*,*) 'Note: you can use an input file "fast_cam.nml" for the camera'
    if (my_id .eq. 0) write(*,*) '      parameters, see jorek/util/fast_camera/fast_cam.nml'
  end if



  ! --- Sanity checks

  ! --- only up to 1 reflection possible at present!
  n_reflections = 1
  if (.not. include_reflections) n_reflections = 0
  if ( (include_reflections) .and. (go_faster_in_core) ) then
    if (my_id .eq. 0) write(*,*)'WARNING: cannot have varying integration step with reflections'
    if (my_id .eq. 0) write(*,*)'         setting go_faster_in_core to FALSE'
    go_faster_in_core = .false. 
  endif
  ! --- Ignore isolenoid if doing wall reflections
  if (include_reflections) solenoid = 0.d0









  ! ********************************************************************** !
  !  Read data file with PEC(Ne,Te) grid (Photon Emissivity Coefficients)  !
  !  WARNING: THIS NEEDS TO BE REPLACED BY EXACT ADAS DATA, IT WAS ONLY    !
  !           MEANT FOR NICE PICTURES OF MAST, NOT QUANTITATIVE DIAGNOSTIC !
  ! ********************************************************************** !
  
  open(123, file=PEC_file, action='read', iostat=ierr)
  
  if ( ierr .ne. 0 ) then
    if (my_id .eq. 0) write(*,*) 'Failed to open PEC data file ',PEC_file,ierr
    if (my_id .eq. 0) write(*,*) 'You can find the my_pec.dat file in jorek/util/fast_camera/my_pec.dat'
    if (my_id .eq. 0) write(*,*) 'Aborting...'
    call MPI_FINALIZE(ierr)
    stop
  else
    
    ! --- File should start with a comment line
    read(123,'(A)') line
    
    ! --- Second line should be "np"
    read(123,'(A)') line
    
    ! --- Third line should be the value of "np"
    read(123,'(A)') line
    read(line,*) PEC_size
    
    ! --- Allocate vectors
    allocate(PEC_dens(PEC_size),PEC_temp(PEC_size),PEC(PEC_size*PEC_size))
    PEC_dens = 0.d0
    PEC_temp = 0.d0
    PEC      = 0.d0
    
    ! --- Then comes the density profile
    read(123,'(A)') line
    do i=1, PEC_size
      read(123,'(A)') line
      read(line,*) PEC_dens(i)
    enddo
    
    ! --- Then comes the temperature profile
    read(123,'(A)') line
    do i=1, PEC_size
      read(123,'(A)') line
      read(line,*) PEC_temp(i)
    enddo
    
    ! --- Then comes the PEC profiles
    read(123,'(A)') line
    do i=1, PEC_size*PEC_size
      read(123,'(A)') line
      read(line,*) PEC(i)
    enddo
    
    close(123)
  endif
    
    
    
    

  ! ********************************************************************** !
  !                     Build lines of sight for camera                    !
  ! ********************************************************************** !
 
  ! --- Allocate the lines of sight coords (one Point plus one Vector plus Light intensity)
  n_pix = n_pixels_hor*n_pixels_ver
  allocate(Xp  (n_reflections+1,n_pix),Yp  (n_reflections+1,n_pix),Zp  (n_reflections+1,n_pix))
  allocate(Xv  (n_reflections+1,n_pix),Yv  (n_reflections+1,n_pix),Zv  (n_reflections+1,n_pix))
  allocate(Xref(n_reflections+2,n_pix),Yref(n_reflections+2,n_pix),Zref(n_reflections+2,n_pix))
  allocate(  step_pix(n_reflections+1,n_pix))
  allocate(n_step_pix(n_reflections+1,n_pix))
  allocate(refl_coef(n_reflections+1))
  allocate(Light(n_pix))
  Light  = 0.d0
  ncount = 0
  Xp   = 0.d0 ; Yp   = 0.d0 ; Zp   = 0.d0
  Xv   = 0.d0 ; Yv   = 0.d0 ; Zv   = 0.d0
  Xref = 0.d0 ; Yref = 0.d0 ; Zref = 0.d0
  
  ! --- Calculate each Point/Vector
  ! --- (note we assume the sensor is "in front" of lens, otherwise you need to flip the image afterwards)
  do i=1, n_pixels_hor
    do j=1, n_pixels_ver
      ncount     = ncount+1
      Xp(1,ncount) = X_cam - focus
      Yp(1,ncount) = Y_cam + pixel_dim*(n_pixels_ver-1)/2 - pixel_dim*(j-1) + focus * sin(angle_ver)
      Zp(1,ncount) = Z_cam - pixel_dim*(n_pixels_hor-1)/2 + pixel_dim*(i-1) + focus * sin(angle_hor)
      Xv(1,ncount) = Xp(1,ncount) - X_cam
      Yv(1,ncount) = Yp(1,ncount) - Y_cam
      Zv(1,ncount) = Zp(1,ncount) - Z_cam
      vec_size   = ( Xv(1,ncount)**2.d0 + Yv(1,ncount)**2.d0 + Zv(1,ncount)**2.d0 )**0.5d0
      Xv(1,ncount) = Xv(1,ncount)/vec_size
      Yv(1,ncount) = Yv(1,ncount)/vec_size
      Zv(1,ncount) = Zv(1,ncount)/vec_size
    enddo
  enddo
  
  
  
  ! ********************************************************************** !
  !       Initialise, import restart, and find axis and xpoint             !
  ! ********************************************************************** !
  
  ! --- Import and initialise
  call initialise_and_broadcast_parameters(my_id, "__NO_FILENAME__")
  do k_tor=1, n_tor
    mode(k_tor) = + int(k_tor / 2) * n_period
  enddo
  !call import_restart(node_list, element_list, 'jorek_restart', rst_format, ierr, .true.)
  if (my_id .eq. 0) call import_restart(node_list, element_list, 'jorek_restart', rst_format, ierr, .true.)
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  call initialise_basis
  call broadcast_elements(my_id, element_list)                ! elements
  call broadcast_nodes(my_id, node_list)                      ! nodes
  call populate_element_rtree(node_list, element_list)        ! rtree
  call broadcast_phys(my_id)                                  ! physics parameters
  call broadcast_equil_state(my_id)                           ! equil_state
  call MPI_Barrier(MPI_COMM_WORLD,ierr)

  ! --- Sanity check for reflection and wall data
  if ( (include_reflections) .and. (n_limiter .eq. 0) ) then
    if (my_id .eq. 0) write(*,*)'WARNING: asking for wall reflections without wall contour'
    if (my_id .eq. 0) write(*,*)'         either set include_reflections to FALSE, or include'
    if (my_id .eq. 0) write(*,*)'         limiter contour (n_limiter,R_limiter,Z_limiter)'
    if (my_id .eq. 0) write(*,*)'         in your input file.'
    if (my_id .eq. 0) write(*,*)'         Aborting...'
    call MPI_FINALIZE(ierr)
    stop
  endif
  if ( (include_reflections) .and. (n_limiter .gt. 0) ) then
    if (n_limiter .lt. 4) then
      if (my_id .eq. 0) write(*,*)'WARNING: asking for wall reflections with bad wall contour'
      if (my_id .eq. 0) write(*,*)'         your wall contour has less than 4 points.'
      if (my_id .eq. 0) write(*,*)'         this doesnt make sense, include a correct limiter contour'
      if (my_id .eq. 0) write(*,*)'         in your input file'
      if (my_id .eq. 0) write(*,*)'         Aborting...'
      call MPI_FINALIZE(ierr)
      stop
    endif
    if ( (R_limiter(1) .ne. R_limiter(n_limiter)) .or. (Z_limiter(1) .ne. Z_limiter(n_limiter)) ) then
      if (my_id .eq. 0) write(*,*)'WARNING: asking for wall reflections with bad wall contour'
      if (my_id .eq. 0) write(*,*)'         your wall contour is not closed! the last point needs'
      if (my_id .eq. 0) write(*,*)'         to be the same as the first one.'
      if (my_id .eq. 0) write(*,*)'         Aborting...'
      call MPI_FINALIZE(ierr)
      stop
    endif
  endif
  
  if (xcase .eq. UPPER_XPOINT) ES%Z_xpoint(1) = -999.0
  if (xcase .eq. LOWER_XPOINT) ES%Z_xpoint(2) = +999.0
  
  ! --- Just in case
  if (central_density .gt. 1.d10) then
    central_ne = central_density
  else
    central_ne = central_density*1.d20
  endif

  ! --- Need the harmonic contributions for the toroidal location
  do i_tor=1, n_tor
    mode(i_tor) = int(i_tor / 2) * n_period
  enddo
  allocate(HZ_tor(n_tor))
  
  
  ! --- MPI loops
  pix_delta = n_pix / n_cpu
  pix_start = my_id*pix_delta + 1
  pix_end   = min(n_pix,(my_id+1)*pix_delta)
  pix_delta = pix_end-pix_start+1


  ! --- MPI barrier
  call MPI_Barrier(MPI_COMM_WORLD,ierr)

  ! --- Dummy call for useless Rtree printout...
  RR = 0.d0 ; ZZ = 0.d0
  call find_RZ(node_list,element_list,RR,ZZ,R_out,Z_out,i_elm,ss,tt,ifail)
  call sleep(3)
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
        
      
  
  ! ********************************************************************** !
  !                 Find intersections with the wall                       !
  ! ********************************************************************** !
  
  if (include_reflections) then

    if (my_id .eq. 0) write(*,*) 'Determining wall reflections...'

    ! --- The step size for integration and number of steps to get to the other side of plasma
    step_ref = 1.0 * step
    n_step = 4.d0*(2.0*X_cam/step_ref)
    
    ! --- For each line of sight...
    do i=pix_start, pix_end
      
      X_tmp  = Xp(1,i)
      Y_tmp  = Yp(1,i)
      Z_tmp  = Zp(1,i)
      inside_wall = 0
      
      ! --- For each step on the line of sight...
      do j=1, n_step
        
        ! --- No need to do anything if we have exited the wall
        if ( (include_reflections) .and. (inside_wall .eq. 10) ) exit
        
        ! --- X,Y,Z-coords
        k = max(1,inside_wall)
        X_tmp = X_tmp + step_ref*Xv(k,i)
        Y_tmp = Y_tmp + step_ref*Yv(k,i)
        Z_tmp = Z_tmp + step_ref*Zv(k,i)
        
        ! --- R,Z,Phi-coords
        RR  = ( X_tmp**2.d0 + Z_tmp**2.d0 )**0.5d0
        ZZ  = Y_tmp
        Phi = atan2(Z_tmp,X_tmp) + toroidal_angle_location
        if (Phi .lt. 0.d0) Phi = Phi + 2.d0*PI
        
        ! --- Check entry/exit of first-wall
        call check_point_is_inside_wall(RR, ZZ, inside_tmp)
        ! --- If we are inside wall
        if (inside_tmp .eq. 1) then
          ! --- This is our first entry
          if (inside_wall .eq. 0) then
            inside_wall = 1
            ! --- If camera is already inside wall
            if (j .eq. 1) then
              Xref(1,i) = X_tmp
              Yref(1,i) = Y_tmp
              Zref(1,i) = Z_tmp
            ! --- If camera is outside wall
            else
              Rl = (/RR_prev,RR/)
              Zl = (/ZZ_prev,ZZ/)
              call wall_intersection(Rl,Zl, i_int, R_int,Z_int, ifail)
              ! --- If we haven't found the intersection (which should never happen in principle)
              if (ifail .ne. 0) then
                !write(*,'(A,4e18.7)')'WARNING: wall intersection not found on entry!',Rl,Zl
                ! --- Take the last point...
                Xref(1,i) = X_tmp
                Yref(1,i) = Y_tmp
                Zref(1,i) = Z_tmp
              ! --- If we found the intersection
              else
                distance = sqrt( (R_int-RR_prev)**2 + (Z_int-ZZ_prev)**2 )
                distance_seg = sqrt( (RR-RR_prev)**2 + (ZZ-ZZ_prev)**2)
                distance = distance / distance_seg
                distance_seg = sqrt( (X_tmp-X_prev)**2 + (Y_tmp-Y_prev)**2 + (Z_tmp-Z_prev)**2 )
                distance = distance * distance_seg
                Xref(1,i) = X_prev + distance*Xv(1,i)
                Yref(1,i) = Y_prev + distance*Yv(1,i)
                Zref(1,i) = Z_prev + distance*Zv(1,i)
              endif
            endif
          endif
        ! --- If we are outside wall
        else
          ! --- If we are exiting wall for 2nd time (after reflection)
          if (inside_wall .eq. 2) then
            Rl = (/RR_prev,RR/)
            Zl = (/ZZ_prev,ZZ/)
            call wall_intersection(Rl,Zl, i_int, R_int,Z_int, ifail)
            ! --- If we haven't found the intersection (which should never happen in principle)
            if (ifail .ne. 0) then
              !write(*,'(A,4e18.7)')'WARNING: wall intersection not found on 2nd exit!',Rl,Zl
              ! --- Take the last point...
              Xref(3,i) = X_tmp
              Yref(3,i) = Y_tmp
              Zref(3,i) = Z_tmp
            ! --- If we found the intersection
            else
              distance = sqrt( (R_int-RR_prev)**2 + (Z_int-ZZ_prev)**2 )
              distance_seg = sqrt( (RR-RR_prev)**2 + (ZZ-ZZ_prev)**2)
              distance = distance / distance_seg
              distance_seg = sqrt( (X_tmp-X_prev)**2 + (Y_tmp-Y_prev)**2 + (Z_tmp-Z_prev)**2 )
              distance = distance * distance_seg
              Xref(3,i) = X_prev + distance*Xv(2,i)
              Yref(3,i) = Y_prev + distance*Yv(2,i)
              Zref(3,i) = Z_prev + distance*Zv(2,i)
            endif
            ! --- Set reflection flag for exit
            inside_wall = 10
            cycle
          endif
          ! --- If we are exiting wall
          if (inside_wall .eq. 1) then
            Rl = (/RR_prev,RR/)
            Zl = (/ZZ_prev,ZZ/)
            call wall_intersection(Rl,Zl, i_int, R_int,Z_int, ifail)
            ! --- If we haven't found the intersection (which should never happen in principle)
            if (ifail .ne. 0) then
              !write(*,'(A,4e18.7)')'WARNING: wall intersection not found on exit!',Rl,Zl
              ! --- Take the last point...
              Xref(2,i) = X_tmp
              Yref(2,i) = Y_tmp
              Zref(2,i) = Z_tmp
              ! --- Just ignore reflection, it probably means that los and wall are almost parallel...
              Xref(3,i) = Xref(2,i)
              Yref(3,i) = Yref(2,i)
              Zref(3,i) = Zref(2,i)
            ! --- If we found the intersection
            else
              distance = sqrt( (R_int-RR_prev)**2 + (Z_int-ZZ_prev)**2 )
              distance_seg = sqrt( (RR-RR_prev)**2 + (ZZ-ZZ_prev)**2)
              distance = distance / distance_seg
              distance_seg = sqrt( (X_tmp-X_prev)**2 + (Y_tmp-Y_prev)**2 + (Z_tmp-Z_prev)**2 )
              distance = distance * distance_seg
              Xref(2,i) = X_prev + distance*Xv(1,i)
              Yref(2,i) = Y_prev + distance*Yv(1,i)
              Zref(2,i) = Z_prev + distance*Zv(1,i)
              ! --- Now need to find the reflection angle
              ! --- The wall vector in (R,Z)-coords
              Phi = atan2(Zref(2,i),Xref(2,i)) + toroidal_angle_location
              if (Phi .lt. 0.d0) Phi = Phi + 2.d0*PI
              Rvec(1) = R_limiter(i_int) ; Rvec(2) = R_limiter(i_int+1)
              Zvec(1) = Z_limiter(i_int) ; Zvec(2) = Z_limiter(i_int+1)
              ! --- In (X,Y,Z)-coords
              do k=1,2
                XXtmp(k) = Rvec(k) * cos(Phi)
                YYtmp(k) = Zvec(k)
                ZZtmp(k) = Rvec(k) * sin(Phi)
              enddo
              XXvec(1) = XXtmp(2) - XXtmp(1)
              YYvec(1) = YYtmp(2) - YYtmp(1)
              ZZvec(1) = ZZtmp(2) - ZZtmp(1)
              ! --- The vector from the centre to the intersection
              XXtmp(1) = 0.d0      ; XXtmp(2) = Xref(2,i)
              YYtmp(1) = Yref(2,i) ; YYtmp(2) = Yref(2,i)
              ZZtmp(1) = 0.d0      ; ZZtmp(2) = Zref(2,i)
              XXvec(2) = XXtmp(2) - XXtmp(1)
              YYvec(2) = YYtmp(2) - YYtmp(1)
              ZZvec(2) = ZZtmp(2) - ZZtmp(1)
              ! --- The cross-product of these two vectors, ie. the vector tangent to the wall in the dPhi direction
              XXvec(3) = YYvec(1)*ZZvec(2) - YYvec(2)*ZZvec(1)
              YYvec(3) = XXvec(2)*ZZvec(1) - XXvec(1)*ZZvec(2)
              ZZvec(3) = XXvec(1)*YYvec(2) - XXvec(2)*YYvec(1)
              distance = sqrt( XXvec(3)**2 + YYvec(3)**2 + ZZvec(3)**2 )
              ! --- In case the wall is parallel to the vector from centre
              if (distance .lt. 1.d-14) then
                ! --- The vector from the centre to the intersection
                XXtmp(1) = 0.d0            ; XXtmp(2) = Xref(2,i)
                YYtmp(1) = Yref(2,i) + 1.0 ; YYtmp(2) = Yref(2,i)
                ZZtmp(1) = 0.d0            ; ZZtmp(2) = Zref(2,i)
                XXvec(2) = XXtmp(2) - XXtmp(1)
                YYvec(2) = YYtmp(2) - YYtmp(1)
                ZZvec(2) = ZZtmp(2) - ZZtmp(1)
                ! --- The cross-product of these two vectors, ie. the vector tangent to the wall in the dPhi direction
                XXvec(3) = YYvec(1)*ZZvec(2) - YYvec(2)*ZZvec(1)
                YYvec(3) = XXvec(2)*ZZvec(1) - XXvec(1)*ZZvec(2)
                ZZvec(3) = XXvec(1)*YYvec(2) - XXvec(2)*YYvec(1)
              endif
              ! --- The cross-product of the wall vector and the dPhi tangent, ie. the normal to the wall
              XXvec(4) = YYvec(1)*ZZvec(3) - YYvec(3)*ZZvec(1)
              YYvec(4) = XXvec(3)*ZZvec(1) - XXvec(1)*ZZvec(3)
              ZZvec(4) = XXvec(1)*YYvec(3) - XXvec(3)*YYvec(1)
              ! --- Normalise
              distance = sqrt( XXvec(4)**2 + YYvec(4)**2 + ZZvec(4)**2 )
              XXvec(4) = XXvec(4) / distance
              YYvec(4) = YYvec(4) / distance
              ZZvec(4) = ZZvec(4) / distance
              ! --- Check direction is inside wall (step by 1mm)
              XXtmp(1) = Xref(2,i) + XXvec(4) * 0.001
              YYtmp(1) = Yref(2,i) + YYvec(4) * 0.001
              ZZtmp(1) = Zref(2,i) + ZZvec(4) * 0.001
              RR_tmp  = ( XXtmp(1)**2.d0 + ZZtmp(1)**2.d0 )**0.5d0
              ZZ_tmp  = YYtmp(1)
              call check_point_is_inside_wall(RR_tmp, ZZ_tmp, inside_tmp2)
              if (inside_tmp2 .ne. 1) then
                XXvec(4) = -XXvec(4)
                YYvec(4) = -YYvec(4)
                ZZvec(4) = -ZZvec(4)
                ! --- Check direction is inside wall (step by 1mm)
                XXtmp(1) = Xref(2,i) + XXvec(4) * 0.001
                YYtmp(1) = Yref(2,i) + YYvec(4) * 0.001
                ZZtmp(1) = Zref(2,i) + ZZvec(4) * 0.001
                RR_tmp  = ( XXtmp(1)**2.d0 + ZZtmp(1)**2.d0 )**0.5d0
                ZZ_tmp  = YYtmp(1)
                call check_point_is_inside_wall(RR_tmp, ZZ_tmp, inside_tmp2)
                ! --- If we exit in both directions, just finish
                if (inside_tmp2 .ne. 1) then
                  Xref(3,i) = Xref(2,i)
                  Yref(3,i) = Yref(2,i)
                  Zref(3,i) = Zref(2,i)
                  ! --- Set reflection flag for exit
                  inside_wall = 10
                  cycle
                endif
              endif
              ! --- Reflection vector: r = v - 2(v.dot.n)n
              distance = Xv(1,i)*XXvec(4) + Yv(1,i)*YYvec(4) + Zv(1,i)*ZZvec(4) ! v.dot.n
              Xv(2,i) = Xv(1,i) - 2.0 * distance * XXvec(4)
              Yv(2,i) = Yv(1,i) - 2.0 * distance * YYvec(4)
              Zv(2,i) = Zv(1,i) - 2.0 * distance * ZZvec(4)
              ! --- Normalise
              distance = sqrt( Xv(2,i)**2 + Yv(2,i)**2 + Zv(2,i)**2)
              Xv(2,i) = Xv(2,i) / distance
              Yv(2,i) = Yv(2,i) / distance
              Zv(2,i) = Zv(2,i) / distance
            endif
            ! --- Set this point to be intersection
            RR    = R_int
            ZZ    = Z_int
            X_tmp = Xref(2,i)
            Y_tmp = Yref(2,i)
            Z_tmp = Zref(2,i)
            ! --- Set reflection flag
            inside_wall = 2
          endif
        endif

        ! --- Record previous point
        RR_prev = RR
        ZZ_prev = ZZ
        X_prev = X_tmp
        Y_prev = Y_tmp
        Z_prev = Z_tmp
 
      enddo ! integration steps
    
      ! --- Set up steps for each reflection
      ! --- Starting points
      Xp(1,i) = Xref(1,i) ; Xp(2,i) = Xref(2,i)
      Yp(1,i) = Yref(1,i) ; Yp(2,i) = Yref(2,i)
      Zp(1,i) = Zref(1,i) ; Zp(2,i) = Zref(2,i)
      ! --- Number/size of steps for first reflection
      distance = sqrt( (Xref(2,i)-Xref(1,i))**2 + (Yref(2,i)-Yref(1,i))**2 + (Zref(2,i)-Zref(1,i))**2 )
      n_step_pix(1,i) = int(distance / step)
      if (n_step_pix(1,i) .gt. 0) step_pix(1,i) = distance / float(n_step_pix(1,i))
      ! --- Number/size of steps for second reflection
      distance = sqrt( (Xref(3,i)-Xref(2,i))**2 + (Yref(3,i)-Yref(2,i))**2 + (Zref(3,i)-Zref(2,i))**2 )
      n_step_pix(2,i) = int(distance / step)
      if (n_step_pix(2,i) .gt. 0) step_pix(2,i) = distance / float(n_step_pix(2,i))
      
    enddo ! pixels
      
    ! --- MPI barrier
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
      
    if (my_id .eq. 0) write(*,*) 'finished reflections.'

  endif ! include_reflections



  
  ! ********************************************************************** !
  !                 Integrate radiation on each line of sight              !
  ! ********************************************************************** !
    
  if (my_id .eq. 0) write(*,*) 'Starting line-of-sight integration...'

  ! --- For each line of sight...
  do i=pix_start, pix_end
    
    ! --- The step size for integration and number of steps to get to the other side of plasma
    if (.not. include_reflections) then
      n_step_pix(1,i) = 2.d0*(2.0*X_cam/step)
      step_pix(1,i) = step
    endif

    ! --- Reflection coefficient on first travel always 1.0
    refl_coef(1) = 1.d0
    if (n_reflections .eq. 1) refl_coef(2) = reflective_coefficient

    ! --- For each reflection
    do i_ref = 1,n_reflections+1
 
      !write(*,'(A,i6,A,i6)')'n_pix : ',i,' out of ',n_pix
      X_tmp  = Xp(i_ref,i)
      Y_tmp  = Yp(i_ref,i)
      Z_tmp  = Zp(i_ref,i)
      ifail  = 1
      inside_core = 0  
      inside_wall = 0

      ! --- For each step on the line of sight...
      do j=1, n_step_pix(i_ref,i)
        
        ! --- If we are entering the plasma core, increase step size
        if ( (go_faster_in_core) .and. (ifail .eq. 0) .and. (psi .lt. 0.6) .and. (inside_core .eq. 0) ) then
          step   = 4.d0*step
          inside_core = 1      
        endif
     
        ! --- If we are exiting the plasma core, reduce step size
        if ( (go_faster_in_core) .and. (ifail .eq. 0) .and. (psi .ge. 0.6) .and. (inside_core .eq. 1) ) then
          step   = step/4.d0
          inside_core = 0      
        endif
        
        ! --- X,Y,Z-coords
        X_tmp = X_tmp + step_pix(i_ref,i)*Xv(i_ref,i)
        Y_tmp = Y_tmp + step_pix(i_ref,i)*Yv(i_ref,i)
        Z_tmp = Z_tmp + step_pix(i_ref,i)*Zv(i_ref,i)
        
        ! --- R,Z,Phi-coords
        RR  = ( X_tmp**2.d0 + Z_tmp**2.d0 )**0.5d0
        ZZ  = Y_tmp
        Phi = atan2(Z_tmp,X_tmp) + toroidal_angle_location
        if (Phi .lt. 0.d0) Phi = Phi + 2.d0*PI
        
        ! --- If we hit the solenoid, stop integrating
        if (RR .lt. solenoid) exit
        
        ! --- Look for this point in the simulation domain
        call find_RZ(node_list,element_list,RR,ZZ,R_out,Z_out,i_elm,ss,tt,ifail)
        
        ! --- Calculate emissivity at this point and add to integration
        if (ifail .eq. 0) then
     
          ! --- Need the harmonic contributions for the toroidal location
          HZ_tor(1)   = 1.d0
          do i_tor=1,(n_tor-1)/2
            HZ_tor(2*i_tor)     = cos(mode(2*i_tor)  *Phi)
            HZ_tor(2*i_tor+1)   = sin(mode(2*i_tor+1)*Phi)
          enddo
     
          ! --- Build variables
          psi   = 0.d0 
          rho   = 0.d0 
          rho_n = 0.d0 
          Te    = 0.d0 
          do i_tor=1,n_tor
            call interp(node_list,element_list,i_elm,var_psi,i_tor,ss,tt,P,P_s,P_t,P_st,P_ss,P_tt)
            psi = psi + P * HZ_tor(i_tor)
            call interp(node_list,element_list,i_elm,var_rho,i_tor,ss,tt,P,P_s,P_t,P_st,P_ss,P_tt)
            rho = rho + P * HZ_tor(i_tor)
            if ( with_TiTe ) then
              call interp(node_list,element_list,i_elm,var_Te,i_tor,ss,tt,P,P_s,P_t,P_st,P_ss,P_tt)
              Te  = Te  + 0.5 * P * HZ_tor(i_tor)
            else
              call interp(node_list,element_list,i_elm,var_T,i_tor,ss,tt,P,P_s,P_t,P_st,P_ss,P_tt)
              Te  = Te  + P * HZ_tor(i_tor)
            endif
            if ( with_neutrals ) then
              call interp(node_list,element_list,i_elm,var_rhon,i_tor,ss,tt,P,P_s,P_t,P_st,P_ss,P_tt)
              rho_n = rho_n + P * HZ_tor(i_tor)
            endif
          enddo
          
          ! --- Normalise psi and denormalise density and temperature
          eV2Joules = 1.602176487d-19
          psi   = (psi-ES%psi_axis)/(ES%psi_bnd-ES%psi_axis) ! we don't include Z-tanh on purpose!
          rho   = rho*central_ne
          rho_n = rho_n*central_ne
          Te    = Te/(central_ne*MU_ZERO*eV2Joules)
          ! --- Make sure everything is positive
          rho   = max(0.d0,rho)
          rho_n = max(1.d-12,rho_n)
          Te    = max(Te_PEC_min,Te)
          
          ! --- Use real visible emission
          if (.not. artificial_light) then
            ! --- Artificial neutral density, if not using neutrals model500
            if ( .not. with_neutrals ) then
              tanh_psi  = 0.5 - 0.5* tanh(+(psi - rhon_prof(6))/rhon_prof(5))
              tanh_zmin = 0.5 - 0.5* tanh(+(ZZ - ES%Z_xpoint(1)-rhon_prof(8 ))/rhon_prof(7))
              tanh_zpls = 0.5 - 0.5* tanh(-(ZZ - ES%Z_xpoint(2)-rhon_prof(10))/rhon_prof(9))
              rho_n = ( rhon_prof(1) - rhon_prof(2) ) * tanh_psi + rhon_prof(2)
              rho_n = rho_n * (1.0 - tanh_zmin) + rhon_prof(3) * tanh_zmin
              rho_n = rho_n * (1.0 - tanh_zpls) + rhon_prof(4) * tanh_zpls
              rho_n = rho_n * central_ne
            endif
           
            ! --- Calculate PEC(Ne,Te)
            PEC_tmp = 0.d0
            do k=2,PEC_size
              if (rho .lt. PEC_dens(k)) then
                if (abs(PEC_dens(k)-rho) .lt. abs(PEC_dens(k-1)-rho)) then
                  PEC_index_Ne = k
                else
                  PEC_index_Ne = k-1            
                endif
                exit
              endif
              if (rho .lt. PEC_dens(1       )) PEC_index_Ne = 0!1
              if (rho .gt. PEC_dens(PEC_size)) PEC_index_Ne = 0!PEC_size
            enddo
            do k=2,PEC_size
              if (Te .lt. PEC_temp(k)) then
                if (abs(PEC_temp(k)-Te) .lt. abs(PEC_temp(k-1)-Te)) then
                  PEC_index_Te = k
                else
                  PEC_index_Te = k-1            
                endif
                exit
              endif
              if (Te .lt. PEC_temp(1       )) PEC_index_Te = 0!1
              if (Te .gt. PEC_temp(PEC_size)) PEC_index_Te = 0!PEC_size
            enddo        
            if ( (PEC_index_Ne .eq. 0) .or. (PEC_index_Te .eq. 0) ) then
              PEC_tmp = 0.d0
            else
              PEC_index = (PEC_index_Ne-1)*PEC_size + PEC_index_Te
              PEC_tmp = PEC(PEC_index)
            endif
           
            ! --- Integrate Emissivity
            Light(i) = Light(i) + refl_coef(i_ref)*step_pix(i_ref,i)*rho_n*rho*PEC_tmp
          
          ! --- Use artificial "light" instead of visible emission
          else
            ! --- User defined "light"
            user_light = RR * exp( -0.5*(psi - 0.98)**2/0.1**2 ) ! gaussian at plasma edge
            ! --- Integrate Emissivity
            Light(i) = Light(i) + refl_coef(i_ref)*step_pix(i_ref,i)*user_light
          endif
                  
        endif
        
      enddo ! integration steps
         
    enddo !n_reflections+1

  enddo ! pixels
  
  ! --- MPI barrier
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
    
  if (my_id .eq. 0) write(*,*) 'Finished line-of-sight integration, gathering data...'

  ! --- MPI collect
  if (my_id .eq. 0) then
    do j=1,n_cpu-1
      call mpi_recv(pix_start,1, MPI_INTEGER, j, j, MPI_COMM_WORLD, status, ierr)
      call mpi_recv(pix_end,  1, MPI_INTEGER, j, j, MPI_COMM_WORLD, status, ierr)
      call mpi_recv(pix_delta,1, MPI_INTEGER, j, j, MPI_COMM_WORLD, status, ierr)
      nrecv = pix_delta
      call mpi_recv(Light(pix_start:pix_end),nrecv, MPI_DOUBLE_PRECISION, j, j, MPI_COMM_WORLD, status, ierr)
    enddo
  else
    call mpi_send(pix_start, 1, MPI_INTEGER, 0, my_id, MPI_COMM_WORLD, ierr)
    call mpi_send(pix_end,   1, MPI_INTEGER, 0, my_id, MPI_COMM_WORLD, ierr)
    call mpi_send(pix_delta, 1, MPI_INTEGER, 0, my_id, MPI_COMM_WORLD, ierr)
    nsend = pix_delta
    call mpi_send(Light(pix_start:pix_end), nsend, MPI_DOUBLE_PRECISION, 0, my_id, MPI_COMM_WORLD, ierr)
  endif

  
  ! --- write image file
  if (my_id .eq. 0) then
    maxlight = maxval(Light) * light_saturation
    ! --- Open image file with PPM P3 format
    open(unit=2,file='fast_camera.ppm',status='unknown')
    write(*,*) 'Now writing PPM (P3) file : ', 'fast_camera.ppm'
    ! --- header
    write(2,'(A)') 'P3'
    write(2,'(2(1x,i4),'' 255 '')')  n_pixels_hor, n_pixels_ver
    ! --- Write to image file
    icnt = 0
    do j=1, n_pixels_ver
      do i=1, n_pixels_hor
        call get_color(colormap, Light(j + (i-1)*n_pixels_ver)/maxlight, rgb)
        do k = 1, 3
          itmp = ichar(rgb(k))
          icnt = icnt + 4
          if (icnt .LT. 60) then
            write(2,fmt='(1x,i3,$)') itmp     ! "$" is not standard.
          else
            write(2,fmt='(1x,i3)') itmp
            icnt = 0
          endif
        enddo
      enddo
    enddo
    write(2,'(A)') ' '
    close(2)
  endif

  
  ! --- Data deallocation
  deallocate(Xp,Yp,Zp)
  deallocate(Xv,Yv,Zv)
  deallocate(Xref,Yref,Zref)
  deallocate(  step_pix)
  deallocate(n_step_pix)
  deallocate(refl_coef)
  deallocate(Light)
  deallocate(HZ_tor)

  ! --- MPI finilise
  call MPI_FINALIZE(ierr)

  if (my_id .eq. 0) write(*,*)'finished...'
 
end program jorek2_fast_camera



















! --- Find intersection between line and segment in Toroidal (R,Z) coordinates
subroutine segment_intersection(Rl,Zl, Rs,Zs, R_int,Z_int, ifail)
  implicit none
  ! --- Routine variables
  real*8, dimension(2), intent(in)  :: Rl,Zl, Rs,Zs
  real*8,               intent(out) :: R_int,Z_int
  integer,              intent(out) :: ifail

  ! --- Internal variables
  real*8, dimension(2)              :: R_diff, Z_diff
  real*8                            :: det, det_l, det_s, det_R, det_Z
  real*8, dimension(2)              :: det_sl
  real*8, parameter                 :: num_buff = 1.d-12

  ! --- Initialise as failure
  R_int = 0.d0
  Z_int = 0.d0
  ifail = 999

  ! --- Define difference vectors
  R_diff = (/ Rl(1) - Rl(2), Rs(1) - Rs(2) /)
  Z_diff = (/ Zl(1) - Zl(2), Zs(1) - Zs(2) /)

  ! --- Check if lines are parallel
  det = R_diff(1)*Z_diff(2) - R_diff(2)*Z_diff(1)
  if (det .eq. 0.d0) return

  ! --- Get lines intersection
  det_l  = Rl(1)*Zl(2) - Rl(2)*Zl(1)
  det_s  = Rs(1)*Zs(2) - Rs(2)*Zs(1)
  det_sl = (/ det_l, det_s /)
  det_R  = det_sl(1)*R_diff(2) - det_sl(2)*R_diff(1)
  det_Z  = det_sl(1)*Z_diff(2) - det_sl(2)*Z_diff(1)
  R_int = det_R / det
  Z_int = det_Z / det

  ! --- Copy intersection if inside segments
  if (      (R_int .le. maxval(Rs)+num_buff) .and. (R_int .ge. minval(Rs)-num_buff) &
      .and. (Z_int .le. maxval(Zs)+num_buff) .and. (Z_int .ge. minval(Zs)-num_buff) &
      .and. (R_int .le. maxval(Rl)+num_buff) .and. (R_int .ge. minval(Rl)-num_buff) &
      .and. (Z_int .le. maxval(Zl)+num_buff) .and. (Z_int .ge. minval(Zl)-num_buff) ) then
    ifail = 0
  endif

end subroutine segment_intersection





! --- Find first intersection between segment and wall (first meaning "closest to first point of segment")
subroutine wall_intersection(Rl,Zl, i_int, R_int, Z_int, ifail)
  use phys_module, only: n_limiter, R_limiter, Z_limiter
  implicit none

  ! --- Routine variables
  real*8, dimension(2), intent(in)  :: Rl,Zl
  integer,              intent(out) :: i_int
  real*8,               intent(out) :: R_int,Z_int
  integer,              intent(out) :: ifail

  ! --- Internal variables
  integer                  :: i_wall
  real*8, dimension(3)     :: int_tmp
  real*8, dimension(2)     :: Rs, Zs
  integer, parameter       :: n_max = 10 ! max 10 intersections?
  integer                  :: n_int, ifail_int
  integer,dimension(n_max) :: i_int_tmp
  real*8, dimension(n_max) :: R_int_tmp, Z_int_tmp
  real*8                   :: diff, diff_min
  integer                  :: i_min, i_find, i_tmp
  real*8                   :: R_tmp, Z_tmp

  ! --- Initialise as failure
  i_int = -1
  R_int = 0.d0
  Z_int = 0.d0
  ifail = 999
  n_int = 0

  ! --- Get all intersections with wall
  do i_wall = 1,n_limiter-1
    Rs = (/ R_limiter(i_wall), R_limiter(i_wall+1) /)
    Zs = (/ Z_limiter(i_wall), Z_limiter(i_wall+1) /)
    call segment_intersection(Rl,Zl, Rs,Zs, R_tmp, Z_tmp, ifail_int)
    if (ifail_int .eq. 0) then
      n_int = n_int + 1
      if (n_int .gt. n_max) cycle
      i_int_tmp(n_int) = i_wall
      R_int_tmp(n_int) = R_tmp
      Z_int_tmp(n_int) = Z_tmp
    endif
  enddo

  ! --- Special case
  if (n_int .eq. 0) return

  ! --- We have at least 1 intersections, take the closest one
  diff_min = 1.d10
  do i_tmp = 1,n_int
    diff = sqrt( (Rl(1)-R_int_tmp(i_tmp))**2 + (Zl(1)-Z_int_tmp(i_tmp))**2 )
    if (diff .lt. diff_min) then
      diff_min = diff
      i_min = i_tmp
    endif
  enddo
  i_int = i_int_tmp(i_min)
  R_int = R_int_tmp(i_min)
  Z_int = Z_int_tmp(i_min)
  ifail = 0
  
end subroutine wall_intersection

















! --- Given density, return red/green/blue coefs for colormap
subroutine get_color(colormap, density, rgb)

  implicit none
  
  ! --- Routine parameters
  integer, intent(in)      :: colormap
  real*8,  intent(in)      :: density
  character*1, intent(out) :: rgb(3)
  
  ! --- Internal parameters
  real*8  ::  red,  gre,  blu
  integer :: ired, igre, iblu
  

  ! --- Heat colorbar from white to yellow to red to black
  if (colormap .eq. 1) then
    if (density .lt. 1.d0/3.d0) then
      red = 1.d0     
      gre = 1.d0
      blu = 1.d0 - density * 3.d0
    elseif (density .lt. 2.d0/3.d0) then
      red = 1.d0
      gre = 1.d0 - (density-1.d0/3.d0) * 3.d0
      blu = 0.d0
    else
      red = 1.d0 - (density-2.d0/3.d0) * 3.d0
      gre = 0.d0
      blu = 0.d0
    endif
  elseif (colormap .eq. 2) then
  ! --- Rainbow colorbar from blue to green to white to yellow to red to black
    if (density .lt. 1.d0/5.d0) then
      red = 0.d0   
      gre = density * 5.d0
      blu = 1.d0 - density * 5.d0
    elseif (density .lt. 2.d0/5.d0) then
      red = (density-1.d0/5.d0) * 5.d0
      gre = 1.d0  
      blu = (density-1.d0/5.d0) * 5.d0
    elseif (density .lt. 3.d0/5.d0) then
      red = 1.d0
      gre = 1.d0  
      blu = 1.d0 - (density-2.d0/5.d0) * 5.d0
    elseif (density .lt. 4.d0/5.d0) then
      red = 1.d0
      gre = 1.d0 - (density-3.d0/5.d0) * 5.d0
      blu = 0.d0
    else
      red = 1.d0 - (density-4.d0/5.d0) * 5.d0
      gre = 0.d0
      blu = 0.d0
    endif
  elseif (colormap .eq. 0) then
  ! --- From black to white
    red = 0.d0 + density
    gre = 0.d0 + density
    blu = 0.d0 + density
  elseif (colormap .eq. -1) then
  ! --- From white to black
    red = 1.d0 - density
    gre = 1.d0 - density
    blu = 1.d0 - density
  endif
  
  ! --- Convert to 255 bit map  
  ired = int(red * 255.0D+00)
  igre = int(gre * 255.0D+00)
  iblu = int(blu * 255.0D+00)
  if (ired .GT. 255) ired = 255
  if (igre .GT. 255) igre = 255
  if (iblu .GT. 255) iblu = 255
  if (ired .LT.   0) ired =   0
  if (igre .LT.   0) igre =   0
  if (iblu .LT.   0) iblu =   0
  rgb(1) = char(ired)
  rgb(2) = char(igre)
  rgb(3) = char(iblu)
      

  return
end subroutine get_color

