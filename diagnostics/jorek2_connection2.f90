program jorek2_connection2
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
use data_structure
use phys_module
use basis_at_gaussian
use elements_nodes_neighbours
use constants
use mod_import_restart
use mod_neighbours
use mod_interp
use mpi
use equil_info

implicit none

real*8,allocatable  :: rp(:), zp(:), R_all(:), Z_all(:), C_all(:)
real*4,allocatable  :: R_strike(:),  Z_strike(:), P_strike(:)        ! position of strike points
real*8,allocatable  :: C_strike(:),  B_strike(:)                     ! connection length, boundary type at strike points
real*8,allocatable  :: T0_strike(:), T_strike(:)                     ! temperature at start and end of fieldline
real*8,allocatable  :: ZN0_strike(:), ZN_strike(:)                   ! density at start and end of fieldline
real*8, allocatable :: PS0_strike(:)                                 ! flux at starting point

real*8,allocatable  :: R_turn(:,:), Z_turn(:,:), C_turn(:,:), C_turn_tmp(:,:)
real*8,allocatable  :: T_turn(:,:), PSI_turn(:,:), ZN_turn(:,:), PSI_turn_norm(:,:), theta_turn(:,:)
integer :: i, j, iside_i, iside_j, ip, i_line, n_lines, i_tor, i_harm, i_var_psi, i_dir, k, m, ns, nt
integer :: i_elm, ifail, i_phi, n_phi, i_turn, n_turns, i_elm_out, i_elm_prev, i_elm_tmp,i_steps, n_turn_max(2)
real*8  :: R_start, Z_start, P_start, R_line, Z_line, s_line, t_line, p_line, s_mid, t_mid, p_mid, s_out, t_out
real*8  :: R, R_s, R_t, R_st, R_ss, R_tt, Z, Z_s, Z_t, Z_st, Z_ss, Z_tt, P, P_s, P_t, P_st, P_ss, P_tt
real*8  :: tol, delta_phi, Zjac, psi_s, psi_t, R_in, Z_in, R_out, Z_out
real*8  :: Rmin, Rmax, Zmin, Zmax, delta_s, delta_t, R_keep, Z_keep
real*8  :: small_delta, small_delta_s, small_delta_t, delta_phi_local, delta_phi_step, total_phi
real*8  :: Rmid,Zmid,Rmid_s,Rmid_t,Zmid_s,Zmid_t, dl2, total_length, length_max, s_ini, t_ini, zl1, zl2, partial(2)
real*8  :: value_out, psi_bnd
real*8  :: element_start_percent
integer :: elm_start, elm_end, elm_delta, local_elm_start, local_elm_end
integer :: my_id, ikeep, n_cpu, ierr, nsend, nrecv, ikeep0, inode1, inode2, i_line0
real*4,allocatable :: RZkeep(:,:),scalars(:,:)
real*4             :: ZERO
integer            :: status(MPI_STATUS_SIZE)
integer            :: nnos, n_scalars, ivtk, i_var, i_strike, i_strike0
character          :: buffer*80, lf*1, str1*12, str2*12
character*12, allocatable :: scalar_names(:)
logical :: psi_theta

integer :: n_stride

namelist /connecvtk_params/ psi_theta, n_turns, n_phi, ns, nt, n_stride, P_start, element_start_percent

call MPI_INIT(IERR)
!required=MPI_THREAD_MULTIPLE
!call MPI_Init_thread(required,provided,StatInfo)
call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)      ! id of each MPI proc
call MPI_COMM_SIZE(MPI_COMM_WORLD, n_cpu, ierr)      ! number of MPI procs
write(*,*) 'my_id = ', my_id


if (my_id .eq. 0 ) then
   write(*,*) '***************************************'
   write(*,*) '* JOREK2_poincare                     *'
   write(*,*) '***************************************'
   write(*,*) ' nperiod : ',n_period
endif

n_scalars = 3
ZERO      = 0.

call initialise_parameters(my_id, "__NO_FILENAME__")

! --- Preset parameters 
psi_theta = .true.
! step at constant delta_phi
n_turns = 100 !500             ! number of toroidal turns to follow a fieldline
n_phi   = 200 !1000            ! number of steps per toroidal turn

ns = 1                          ! number of (s) starting points within one element
nt = 1                          ! number of (t) starting points within one element
n_stride =  1                   ! interval of elements between starting points
P_start  = PI/4.!0.d0
element_start_percent = 0.25

! --- Read parameters from namelist file 'connecvtk.nml' if it exists
open(42, file='connecvtk.nml', action='read', status='old', iostat=ierr)
if ( ierr == 0 ) then
if (my_id .eq. 0 ) then
   write(*,*) 'Reading parameters from connecvtk.nml namelist.'
endif
read(42,connecvtk_params)
close(42)
end if

if (my_id .eq. 0 ) then
   write(*,*)
   write(*,*) 'Parameters:'
   write(*,*) '-----------'
   write(*,*) 'psi_theta = ', psi_theta
   write(*,*) 'n_turns = ', n_turns
   write(*,*) 'n_phi = ', n_phi
   write(*,*) 'ns = ', ns
   write(*,*) 'nt = ', nt
   write(*,*) 'n_stride = ', n_stride
   write(*,*) 'element_start = ', element_start_percent, ' percent of nb_elements'
endif


do i_tor=1, n_tor
  mode(i_tor) = + int(i_tor / 2) * n_period
  if (my_id .eq. 0 ) then
     write(*,*) ' toroidal mode numbers : ',i_tor,mode(i_tor)
  endif
enddo

call import_restart(node_list,element_list, 'jorek_restart', rst_format, ierr, .true.)

call initialise_basis                                       ! define the basis functions at the Gaussian points

if (my_id .eq. 0 ) then
   write(*,*) 'central_density = ', central_density
endif

call broadcast_elements(my_id, element_list)                ! elements
call broadcast_nodes(my_id, node_list)                      ! nodes
call broadcast_phys(my_id)                                  ! physics parameters

!-----------------------------------------------------------define element neighbours
write (*,*) 'number of elements= ', element_list%n_elements

allocate(element_neighbours(4,element_list%n_elements))

element_neighbours = 0

do i=1,element_list%n_elements

  do j=i+1,element_list%n_elements

    if (neighbours(node_list,element_list%element(i),element_list%element(j),iside_i,iside_j)) then
      element_neighbours(iside_i,i) = j
      element_neighbours(iside_j,j) = i
    endif

  enddo
enddo

delta_phi = 2.d0 * PI / float(n_period*n_phi)
tol       = 1.d-6!1.e-6

i_var_psi = 1                                  ! the index of the magnetic flux variable

n_lines = int(element_list%n_elements * ns * nt / n_stride)   ! number of starting points

allocate(R_strike(n_lines),Z_strike(n_lines),P_strike(n_lines),C_strike(n_lines),B_strike(n_lines))
allocate(T0_strike(n_lines),T_strike(n_lines),ZN0_strike(n_lines),ZN_strike(n_lines),PS0_strike(n_lines))

allocate(R_all(n_lines),Z_all(n_lines),C_all(n_lines))
allocate(R_turn(n_turns+1,2),Z_turn(n_turns+1,2),C_turn(n_turns+1,2),C_turn_tmp(n_turns+1,2))
allocate(T_turn(n_turns+1,2),PSI_turn(n_turns+1,2),ZN_turn(n_turns+1,2))
!allocate(T_turn(n_turns+1,2),PSI_turn(n_turns+1,2),ZN_turn(n_turns+1,2),PSI_turn_norm(n_turns+1,2),theta_turn(n_turns+1,2))

R_all     = 0.d0; Z_all     = 0.d0; C_all     = 0.d0
R_strike  = 0.d0; Z_strike  = 0.d0; P_strike  = 0.d0;  C_strike   = 0.d0
T0_strike = 0.d0; T_strike  = 0.d0; ZN0_strike = 0.d0; ZN_strike  = 0.d0; PS0_strike = 0.d0
R_turn    = 0.d0; Z_turn    = 0.d0; C_turn    = 0.d0;  C_turn_tmp = 0.d0

Rmin = 1.d20; Rmax = -1.d20; Zmin = 1.d20; Zmax=-1.d20
do i=1,node_list%n_nodes
  Rmin = min(Rmin,node_list%node(i)%x(1,1,1))
  Rmax = max(Rmax,node_list%node(i)%x(1,1,1))
  Zmin = min(Zmin,node_list%node(i)%x(1,1,2))
  Zmax = max(Zmax,node_list%node(i)%x(1,1,2))
enddo

!------------------------------------------------- find x-point(s)
xcase = LOWER_XPOINT  ! Why we have this line?
if (xpoint) then
  psi_bnd = ES%psi_xpoint(1)
  if( ES%active_xpoint .eq. UPPER_XPOINT ) then
    psi_bnd = ES%psi_xpoint(2)
  endif
else
  psi_bnd = 0.d0
endif

if (my_id .eq. 0 ) then
write(*,*) ' xcase,1st x-point:R,Z,psi: ',xcase, ES%R_xpoint(1),ES%Z_xpoint(1),ES%psi_xpoint(1),psi_bnd
!   write(*,*) ' PSI_XPOINT : ',ES%psi_xpoint,ES%i_elm_xpoint
   write(*,*) ' PSI_AXIS : ',ES%psi_axis,ES%i_elm_axis
   write(*,*) ' RZ_AXIS : ', ES%R_axis, ES%Z_axis
endif

!call MPI_Barrier(MPI_COMM_WORLD,ierr)
i_line   = 0
i_strike = 0

write(*,*) ' number of elements : ',element_list%n_elements

if (element_start_percent .ne. 0.) then
   elm_start = element_list%n_elements*element_start_percent
else
   elm_start = 1
end if
elm_end   = element_list%n_elements

write(*,*) ' elm start, end : ',elm_start, elm_end

elm_delta = (elm_end - elm_start) / n_cpu

local_elm_start = elm_start + my_id*elm_delta + 1
local_elm_end   = min(elm_end,elm_start+(my_id+1)*elm_delta)

write(*,*) my_id, local_elm_start, local_elm_end

!nkeep = (local_elm_end - local_elm_start) * ns * nt * n_turns

ikeep = 0

allocate(RZkeep(2,1000000),scalars(1000000,n_scalars))


do i = local_elm_start, local_elm_end, n_stride

  do k=1, ns

    s_ini = real(k)/real(ns+1)

    do m=1, nt

      t_ini = real(m)/real(nt+1)

      i_line = i_line + 1

      R_turn     = 0.d0
      Z_turn     = 0.d0
      C_turn_tmp = 0.d0

      do i_dir = -1,1,2

      s_line = s_ini
      t_line = t_ini

      delta_phi = 2.d0 * PI * float(i_dir) / float(n_period*n_phi)


      call interp_RZ(node_list,element_list,i,s_line,t_line,R_out,R_s,R_t,R_st,R_ss,R_tt,Z_out,Z_s,Z_t,Z_st,Z_ss,Z_tt)

      total_length = 0.d0
      total_phi    = 0.d0

      i_elm   = i
      R_start = R_out
      Z_start = Z_out

     ! write (*,*) 'i_line,R_start,Z_start',i_line,R_start,Z_start
      R_all(i_line) = R_start
      Z_all(i_line) = Z_start

      R_turn(1,(i_dir+1)/2+1) = R_start
      Z_turn(1,(i_dir+1)/2+1) = Z_start
      C_turn(1,(i_dir+1)/2+1) = 0.d0

      call var_value(i_elm,6,s_line,t_line,P_start,T_turn(1,(i_dir+1)/2+1))
      call var_value(i_elm,1,s_line,t_line,P_start,PSI_turn(1,(i_dir+1)/2+1))
      call var_value(i_elm,5,s_line,t_line,P_start,ZN_turn(1,(i_dir+1)/2+1))
      !PSI_turn_norm (1,(i_dir+1)/2+1)= (PSI_turn(1,(i_dir+1)/2+1)-ES%psi_axis)/(psi_bnd-ES%psi_axis)

      R_line = R_start
      Z_line = Z_start
      p_line = P_start

      i_strike = i_strike + 1
      ZN0_strike(i_strike) = ZN_turn(1,(i_dir+1)/2+1)
      T0_strike(i_strike)  = T_turn(1,(i_dir+1)/2+1)
      PS0_strike(i_strike) = PSI_turn(1,(i_dir+1)/2+1)

      do i_turn = 1, n_turns                 ! loop over toroidal turns

        n_turn_max((i_dir+1)/2+1) = i_turn
!write (*,*) 'n_turn_max = ', n_turn_max(:)
        do i_phi=1,n_phi                     ! loop over steps in toroidal angle

          delta_phi_local = 0.d0

          i_steps = 0                        ! loop inside one element

          do while ((abs(delta_phi_local) .lt. abs(delta_phi)) .and. (i_steps .lt. 100) )

            i_steps = i_steps + 1

            delta_phi_step = delta_phi - delta_phi_local


            call step(i_elm,s_line,t_line,p_line,delta_phi_step,delta_s,delta_t,R,Z,R_s,R_t,Z_s,Z_t)

            s_mid = s_line + 0.5d0 * delta_s
            t_mid = t_line + 0.5d0 * delta_t
            p_mid = p_line + 0.5d0 * delta_phi_step


            call step(i_elm,s_mid,t_mid,p_mid,delta_phi_step,delta_s,delta_t,Rmid,Zmid,Rmid_s,Rmid_t,Zmid_s,Zmid_t)

            small_delta_s = 1.d0

            if  (s_line + delta_s .gt. 1.d0) then         ! step to element boundary, not beyond

              small_delta_s = (1.d0 - s_line)/delta_s

            elseif  (s_line + delta_s .lt. 0.d0) then     ! step to element boundary, not beyond

              small_delta_s = abs(s_line/delta_s)

            endif

            small_delta_t = 1.d0

            if  (t_line + delta_t .gt. 1.d0)  then        ! step to element boundary, not beyond

              small_delta_t = (1.d0 - t_line)/delta_t

            elseif  (t_line + delta_t .lt. 0.d0)  then    ! step to element boundary, not beyond

              small_delta_t = abs(t_line/delta_t)

            endif

            small_delta = min(small_delta_s, small_delta_t)

            if (small_delta .lt. 1.d0)  then             ! this step is crossing the boundary

              s_mid = s_line + 0.5d0 * small_delta * delta_s
              t_mid = t_line + 0.5d0 * small_delta * delta_t
              p_mid = p_line + 0.5d0 * small_delta * delta_phi_step


              call step(i_elm,s_mid,t_mid,p_mid,delta_phi_step,delta_s,delta_t,Rmid,Zmid,Rmid_s,Rmid_t,Zmid_s,Zmid_t)

              if (small_delta_s .lt. small_delta_t) then

                if (s_line + delta_s .gt. 1.d0) then     ! crossing boundary 2 or 4 at s=1

                  s_line = 1.d0
                  t_line = t_line + small_delta * delta_t
                  p_line = p_line + small_delta * delta_phi_step

                  dl2 = (Rmid_s**2 + Zmid_s**2)*delta_s**2 + (Rmid_t**2 + Zmid_t**2)*delta_t**2 + Rmid**2 * delta_phi_step**2
                  dl2 = dl2 * small_delta**2

                  call interp_RZ(node_list,element_list,i_elm,s_line,t_line,R_in,R_s,R_t,R_st,R_ss,R_tt,Z_in,Z_s,Z_t,Z_st,Z_ss,Z_tt)

                  i_elm_prev = i_elm
                  i_elm      = element_neighbours(2,i_elm_prev)

                  if (i_elm .ne. 0) then

                    i_elm_tmp  = element_neighbours(4,i_elm)

                    if (i_elm_prev .ne. i_elm_tmp) write(*,*) ' WARNING : CHANGE OF ORIENTATION (1)'

                    s_line = 0.d0

                    call interp_RZ(node_list,element_list,i_elm,s_line,t_line,R_out,R_s,R_t,R_st,R_ss,R_tt,Z_out,Z_s,Z_t,Z_st,Z_ss,Z_tt)

                    if ( (abs(R_in - R_out) .gt. 1.d-8) .or. (abs(Z_in - Z_out) .gt. 1.d-8)) &
                      write(*,'(A,2i6,4f12.4)') ' error in element change (1) ',i_elm_prev,i_elm,R_in,R_out,Z_in,Z_out

                  else ! crossing an outer boundary

                    R_strike(i_strike) = R_in
                    Z_strike(i_strike) = Z_in
                    P_strike(i_strike) = p_line
                    C_strike(i_strike) = total_length + sqrt(abs(dl2))
                    
                    call var_value(i_elm_prev,6,s_line,t_line,p_line,T_strike(i_strike))
                    call var_value(i_elm_prev,5,s_line,t_line,p_line,ZN_strike(i_strike))

                    inode1 = element_list%element(i_elm_prev)%vertex(2)
                    inode2 = element_list%element(i_elm_prev)%vertex(3)

                    if ((node_list%node(inode1)%boundary .ne. 0) .and. (node_list%node(inode2)%boundary .ne. 0)) then
                      B_strike(i_strike) = min(node_list%node(inode1)%boundary,node_list%node(inode2)%boundary)
                    else
                      write(*,*) 'error : leaving domain but not at a boundary (s=1)!',inode1,inode2,R_in,Z_in
                    endif

                  endif

                elseif (s_line + delta_s .lt. 0.d0) then ! crossing boundary 2 or 4 at s=0

                  s_line = 0.d0
                  t_line = t_line + small_delta * delta_t
                  p_line = p_line + small_delta * delta_phi_step

                  dl2 = (Rmid_s**2 + Zmid_s**2)*delta_s**2 + (Rmid_t**2 + Zmid_t**2)*delta_t**2 + Rmid**2 * delta_phi_step**2
                  dl2 = dl2 * small_delta**2

                  call interp_RZ(node_list,element_list,i_elm,s_line,t_line,R_in,R_s,R_t,R_st,R_ss,R_tt,Z_in,Z_s,Z_t,Z_st,Z_ss,Z_tt)

                  i_elm_prev = i_elm
                  i_elm      = element_neighbours(4,i_elm_prev)

                  if (i_elm .ne. 0) then

                    i_elm_tmp  = element_neighbours(2,i_elm)

                    if (i_elm_prev .ne. i_elm_tmp) write(*,*) ' WARNING : CHANGE OF ORIENTATION (2)'

                    s_line = 1.d0

                    call interp_RZ(node_list,element_list,i_elm,s_line,t_line,R_out,R_s,R_t,R_st,R_ss,R_tt,Z_out,Z_s,Z_t,Z_st,Z_ss,Z_tt)

                    if ( (abs(R_in - R_out) .gt. 1.d-8) .or. (abs(Z_in - Z_out) .gt. 1.d-8))  &
                      write(*,'(A,2i6,4f12.4)') ' error in element change (2) ',i_elm_prev,i_elm,R_in,R_out,Z_in,Z_out

                  else ! crossing an outer boundary

                    R_strike(i_strike) = R_in
                    Z_strike(i_strike) = Z_in
                    P_strike(i_strike) = p_line
                    C_strike(i_strike) = total_length + sqrt(abs(dl2))
!write (*,*) 'popopop 5', my_id, i, i_elm

                    call var_value(i_elm_prev,6,s_line,t_line,p_line,T_strike(i_strike))
                    call var_value(i_elm_prev,5,s_line,t_line,p_line,ZN_strike(i_strike))

                    inode1 = element_list%element(i_elm_prev)%vertex(4)
                    inode2 = element_list%element(i_elm_prev)%vertex(1)

                    if ((node_list%node(inode1)%boundary .ne. 0) .and. (node_list%node(inode2)%boundary .ne. 0)) then
                      B_strike(i_strike) = min(node_list%node(inode1)%boundary,node_list%node(inode2)%boundary)
                    else
                      write(*,*) 'error : leaving domain but not at a boundary (s=0)!',inode1,inode2,R_in,Z_in
                    endif

                  endif

                endif

              else

                if (t_line + delta_t .gt. 1.d0) then  ! crossing boundary 1 or 3 at t=1

                  s_line = s_line + small_delta * delta_s
                  t_line = 1.d0
                  p_line = p_line + small_delta * delta_phi_step

                  dl2 = (Rmid_s**2 + Zmid_s**2)*delta_s**2 + (Rmid_t**2 + Zmid_t**2)*delta_t**2 + Rmid**2 * delta_phi_step**2
                  dl2 = dl2 * small_delta**2

                  call interp_RZ(node_list,element_list,i_elm,s_line,t_line,R_in,R_s,R_t,R_st,R_ss,R_tt,Z_in,Z_s,Z_t,Z_st,Z_ss,Z_tt)

                  i_elm_prev = i_elm
                  i_elm      = element_neighbours(3,i_elm_prev)

                  if (i_elm .ne. 0) then

                    i_elm_tmp  = element_neighbours(1,i_elm)

                    if (i_elm_prev .ne. i_elm_tmp) write(*,*) ' WARNING : CHANGE OF ORIENTATION (3)'

                    t_line = 0.d0

                    call interp_RZ(node_list,element_list,i_elm,s_line,t_line,R_out,R_s,R_t,R_st,R_ss,R_tt,Z_out,Z_s,Z_t,Z_st,Z_ss,Z_tt)

                    if ( (abs(R_in - R_out) .gt. 1.d-8) .or. (abs(Z_in - Z_out) .gt. 1.d-8))  &
                      write(*,'(A,2i6,4f12.4)') ' error in element change (3) ',i_elm_prev,i_elm,R_in,R_out,Z_in,Z_out

                  else ! crossing an outer boundary

                    R_strike(i_strike) = R_in
                    Z_strike(i_strike) = Z_in
                    P_strike(i_strike) = p_line
                    C_strike(i_strike) = total_length + sqrt(abs(dl2))

                    call var_value(i_elm_prev,6,s_line,t_line,p_line,T_strike(i_strike))
                    call var_value(i_elm_prev,5,s_line,t_line,p_line,ZN_strike(i_strike))
!write (*,*) 'popopop 6', my_id, i, i_elm

                    inode1 = element_list%element(i_elm_prev)%vertex(3)
                    inode2 = element_list%element(i_elm_prev)%vertex(4)

                    if ((node_list%node(inode1)%boundary .ne. 0) .and. (node_list%node(inode2)%boundary .ne. 0)) then
                      B_strike(i_strike) = min(node_list%node(inode1)%boundary,node_list%node(inode2)%boundary)
                    else
                      write(*,*) 'error : leaving domain but not at a boundary (t=1)!',inode1,inode2,R_in,Z_in
                    endif

                  endif

                elseif (t_line + delta_t .lt. 0.d0) then  ! crossing boundary 1 or 3 at t=0

                  s_line = s_line + small_delta * delta_s	
                  t_line = 0.d0
                  p_line = p_line + small_delta * delta_phi_step

                  dl2 = (Rmid_s**2 + Zmid_s**2)*delta_s**2 + (Rmid_t**2 + Zmid_t**2)*delta_t**2 + Rmid**2 * delta_phi_step**2
                  dl2 = dl2 * small_delta**2

                  call interp_RZ(node_list,element_list,i_elm,s_line,t_line,R_in,R_s,R_t,R_st,R_ss,R_tt,Z_in,Z_s,Z_t,Z_st,Z_ss,Z_tt)

                  i_elm_prev = i_elm
                  i_elm      = element_neighbours(1,i_elm_prev)

                  if (i_elm .ne. 0) then

                    i_elm_tmp  = element_neighbours(3,i_elm)

                    if (i_elm_prev .ne. i_elm_tmp) write(*,*) ' WARNING : CHANGE OF ORIENTATION (4)'

                    t_line = 1.d0

                    call interp_RZ(node_list,element_list,i_elm,s_line,t_line,R_out,R_s,R_t,R_st,R_ss,R_tt,Z_out,Z_s,Z_t,Z_st,Z_ss,Z_tt)

                    if ( (abs(R_in - R_out) .gt. 1.d-8) .or. (abs(Z_in - Z_out) .gt. 1.d-8))  &
                      write(*,'(A,2i6,4f12.4)') ' error in element change (4) ',i_elm_prev,i_elm,R_in,R_out,Z_in,Z_out

                  else ! crossing an outer boundary

                    R_strike(i_strike) = R_in
                    Z_strike(i_strike) = Z_in
                    P_strike(i_strike) = p_line
                    C_strike(i_strike) = total_length + sqrt(abs(dl2))
!write (*,*) 'popopop 7', my_id, i, i_elm

                    call var_value(i_elm_prev,6,s_line,t_line,p_line,T_strike(i_strike))
                    call var_value(i_elm_prev,5,s_line,t_line,p_line,ZN_strike(i_strike))
!write (*,*) 'popopop 8', my_id, i, i_elm

                    inode1 = element_list%element(i_elm_prev)%vertex(1)
                    inode2 = element_list%element(i_elm_prev)%vertex(2)

                    if ((node_list%node(inode1)%boundary .ne. 0) .and. (node_list%node(inode2)%boundary .ne. 0)) then
                      B_strike(i_strike) = min(node_list%node(inode1)%boundary,node_list%node(inode2)%boundary)
                    else
                      write(*,*) 'error : leaving domain but not at a boundary (t=0)!',inode1,inode2,R_in,Z_in
                    endif

                  endif

                endif

              endif

            else  ! this step remains within the element

              s_line = s_line + delta_s
              t_line = t_line + delta_t
              p_line = p_line + delta_phi_step

              dl2 = (Rmid_s**2 + Zmid_s**2)*delta_s**2 + (Rmid_t**2 + Zmid_t**2)*delta_t**2 + Rmid**2 * delta_phi_step**2

              small_delta = 1.d0

              call interp_RZ(node_list,element_list,i_elm,s_line,t_line,R_in,R_s,R_t,R_st,R_ss,R_tt,Z_in,Z_s,Z_t,Z_st,Z_ss,Z_tt)

            endif

            delta_phi_local = delta_phi_local + small_delta * delta_phi_step

            total_length = total_length + sqrt(abs(dl2))
            total_phi    = total_phi    + small_delta * delta_phi_step

            if (i_elm .eq. 0) exit

          enddo  ! end of loop over steps within one element

          if (i_elm .eq. 0) exit

        enddo    ! end of a 2Pi turn (or before if end of open field line)

        if (i_elm .eq. 0) exit

        call interp_RZ(node_list,element_list,i_elm,s_line,t_line,R_in,R_s,R_t,R_st,R_ss,R_tt,Z_in,Z_s,Z_t,Z_st,Z_ss,Z_tt)

        R_turn(i_turn+1,(i_dir+1)/2+1) = R_in
        Z_turn(i_turn+1,(i_dir+1)/2+1) = Z_in
        C_turn_tmp(i_turn+1,(i_dir+1)/2+1) = total_length

        call var_value(i_elm,6,s_line,t_line,p_line,T_turn(i_turn+1,(i_dir+1)/2+1))
        call var_value(i_elm,1,s_line,t_line,p_line,PSI_turn(i_turn+1,(i_dir+1)/2+1))
        call var_value(i_elm,5,s_line,t_line,p_line,ZN_turn(i_turn+1,(i_dir+1)/2+1))
        !PSI_turn_norm (1,(i_dir+1)/2+1)= (PSI_turn(1,(i_dir+1)/2+1)-ES%psi_axis)/(psi_bnd-ES%psi_axis)

      enddo  ! end of loop over toroidal turns

      if (i_elm .ne. 0) then  ! field line still in domain, after n_turn turns
        R_strike(i_strike) = R_in
        Z_strike(i_strike) = Z_in
        P_strike(i_strike) = p_line
        C_strike(i_strike) = 0.d0          ! to be done (total_length needs correction)
        B_strike(i_strike) = 0
        call var_value(i_elm,6,s_line,t_line,p_line,T_strike(i_strike))
        call var_value(i_elm,5,s_line,t_line,p_line,ZN_strike(i_strike))
      endif


      if (i_dir .eq. -1) then  
        C_all(i_line) = total_length
        partial(1)    = total_length
      else
        C_all(i_line) = min(C_all(i_line),total_length)
        partial(2)    = total_length
      endif

      enddo  ! end of two directions

!------------------------------- correct the connection lengths
      do i_turn = 1, n_turn_max(1)
        C_turn(i_turn,1) = partial(1) - c_turn_tmp(i_turn,1)
      enddo
      do i_turn = 1, n_turn_max(2)
        C_turn(i_turn,2) = partial(2) - c_turn_tmp(i_turn,2)
      enddo

      do i_turn=1,n_turn_max(1)+1                  ! keep only field lines starting inside the plasma
!write (*,*) 'test_boucle1 = ', (R_turn(i_turn,1) .gt. 0.d0)
        if (R_turn(i_turn,1) .gt. 0.d0) then

          zl1 = C_turn(i_turn,1)
          zl2 = C_turn(1,1) - C_turn(i_turn,1) + C_turn(1,2) 
!!$write (*,*) 'test_boucle2.1 = ', (PSI_turn(1,1)    .lt. ES%psi_xpoint(1))
!!$write (*,*) 'test_boucle2.2 = ', (Z_turn(1,1)      .gt. ES%Z_xpoint(1))
!!$write (*,*) 'test_boucle2.3 = ', (Z_turn(i_turn,1) .lt. -2.d0)
          if ( (  (PSI_turn(1,1).le. psi_bnd)  &
               .and. (Z_turn(1,1) .ge. ES%Z_xpoint(1)) ) &    !.and. (Z_turn(1,1).le.ES%Z_xpoint(2))) then    
               .and. (Z_turn(i_turn,1) .lt. 2.d0)  )  then

            if (n_turn_max(1) .lt. n_turns) then
              ikeep = ikeep + 1
!write(*,*) 'ikeep1 = ', ikeep
              if(psi_theta) then
                 RZkeep(1,ikeep) = ( PSI_turn(i_turn,1) - ES%psi_axis ) / (psi_bnd - ES%psi_axis )
                 RZkeep(2,ikeep) = atan2( (Z_turn(i_turn,1) - ES%Z_axis) , (R_turn(i_turn,1) - ES%R_axis) ) / (2.d0*PI)
!write(*,*) '1my_id, psi, theta = ', my_id, RZkeep(1,ikeep), RZkeep(2,ikeep)
              else
                 RZkeep(1,ikeep)            = R_turn(i_turn,1)
                 RZkeep(2,ikeep)            = Z_turn(i_turn,1)
              endif
              scalars(ikeep,1:n_scalars) = (/ min(zl1,zl2),T_turn(1,1), PSI_turn(i_turn,1) /)
            else
              ikeep = ikeep + 1
!write(*,*) 'ikeep1 = ', ikeep
              if(psi_theta) then
                 RZkeep(1,ikeep) = ( PSI_turn(i_turn,1) - ES%psi_axis ) / (psi_bnd - ES%psi_axis )
                 RZkeep(2,ikeep) = atan2( (Z_turn(i_turn,1) - ES%Z_axis) , (R_turn(i_turn,1) - ES%R_axis) ) / (2.d0*PI)
!write(*,*) '2my_id, psi, theta = ', my_id, RZkeep(1,ikeep), RZkeep(2,ikeep)
              else
                 RZkeep(1,ikeep)            = R_turn(i_turn,1)
                 RZkeep(2,ikeep)            = Z_turn(i_turn,1)
              endif
              scalars(ikeep,1:n_scalars) = (/ maxval(partial),T_turn(1,1), PSI_turn(i_turn,1) /)
            endif 

         endif ! attention: ici c'est commente dans l'ancienne version

        Endif
      enddo

      do i_turn=1,n_turn_max(2)+1          ! keep only field lines starting inside the plasma

        if (R_turn(i_turn,2) .gt. 0.d0) then

          zl1 = C_turn(i_turn,2)
          zl2 = C_turn(1,2) - C_turn(i_turn,2) + C_turn(1,1) 

          if (     ( (PSI_turn(1,2).le. psi_bnd)  &
               .and. (Z_turn(1,2).ge. ES%Z_xpoint(1)) ) & !.and.(Z_turn(1,2).le. ES%Z_xpoint(2))) then 
               .and. (Z_turn(i_turn,2) .lt. 2.d0) )  then
!           if (Z_turn(i_turn,1) .lt. -2.d0) then
            if (n_turn_max(2) .lt. n_turns) then
              ikeep = ikeep + 1
!write(*,*) 'ikeep2 = ', ikeep
              if(psi_theta) then
                 RZkeep(1,ikeep) = ( PSI_turn(i_turn,2)  - ES%psi_axis ) / (psi_bnd - ES%psi_axis )
                 RZkeep(2,ikeep) = atan2( (Z_turn(i_turn,2) - ES%Z_axis) , (R_turn(i_turn,2) - ES%R_axis) ) / (2.d0*PI)
!write(*,*) '3my_id, psi, theta = ', my_id, RZkeep(1,ikeep), RZkeep(2,ikeep)
              else
                 RZkeep(1,ikeep)            = R_turn(i_turn,2)
                 RZkeep(2,ikeep)            = Z_turn(i_turn,2)
              endif
              scalars(ikeep,1:n_scalars) = (/ min(zl1,zl2),T_turn(1,2), PSI_turn(i_turn,2) /)
            else
              ikeep = ikeep + 1
!write(*,*) 'ikeep2 = ', ikeep
              if(psi_theta) then
                 RZkeep(1,ikeep) = ( PSI_turn(i_turn,2) - ES%psi_axis ) / (psi_bnd - ES%psi_axis )
                 RZkeep(2,ikeep) = atan2( (Z_turn(i_turn,2) - ES%Z_axis) , (R_turn(i_turn,2) - ES%R_axis) ) / (2.d0*PI)
!write(*,*) '4my_id, psi, theta = ', my_id, RZkeep(1,ikeep), RZkeep(2,ikeep)
              else
                 RZkeep(1,ikeep)            = R_turn(i_turn,2)
                 RZkeep(2,ikeep)            = Z_turn(i_turn,2)
              endif
              scalars(ikeep,1:n_scalars) = (/ maxval(partial),T_turn(1,2), PSI_turn(i_turn,2) /)
            endif

          endif

        endif
      enddo
if ((i == local_elm_start) .or.(i >= local_elm_end)) then
write (*,*) 'popopop 9', my_id, i, scalars(ikeep,1:n_scalars)
endif

    enddo  ! end over loop over starting points within one element ( ns)
  enddo    ! end over loop over starting points within one element ( nt)


enddo ! end of loop over elements

!----------------------------------------------- write to VTK file (one after the other)
ikeep0  = ikeep
write (*,*) 'ikeep = ', ikeep, 'my_id = ', my_id
write (*,*) 'ikeep0 = ', ikeep0, 'my_id = ', my_id

call MPI_Barrier(MPI_COMM_WORLD,ierr)

call MPI_Reduce(ikeep,nnos,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)

if (my_id .eq. 0) write(*,*) ' number of points_ikeep : ',nnos

ivtk = 22                 ! an arbitrary unit number for the VTK output file
n_scalars = 3             ! number of scalars to write to the VTK output file

allocate(scalar_names(n_scalars))

scalar_names  = (/ 'length_m    ','T_start_keV ','psi_norm    ' /)

lf = char(10) ! line feed character

if (my_id .eq. 0) then
#ifdef IBM_MACHINE
  open(unit=ivtk,file='connection_new.vtk',form='unformatted',access='stream',status='replace')
#else
  open(unit=ivtk,file='connection_new.vtk',form='unformatted',access='stream',convert='BIG_ENDIAN',status='replace')
#endif

  buffer = '# vtk DataFile Version 3.0'//lf                                             ; write(ivtk) trim(buffer)
  buffer = 'vtk output'//lf                                                             ; write(ivtk) trim(buffer)
  buffer = 'BINARY'//lf                                                                 ; write(ivtk) trim(buffer)
  buffer = 'DATASET UNSTRUCTURED_GRID'//lf//lf                                          ; write(ivtk) trim(buffer)
 
  ! POINTS SECTION
  write(str1(1:12),'(i12)') nnos
  buffer = 'POINTS '//str1//'  float'//lf                                               ; write(ivtk) trim(buffer)
endif

if (my_id .eq. 0) then
   write(ivtk) ( (/RZkeep(1,i), RZkeep(2,i), ZERO /),i=1,ikeep0)
endif

if (my_id .eq. 0) then
  do j=1,n_cpu-1
    call mpi_recv(ikeep,1, MPI_INTEGER, j, j, MPI_COMM_WORLD, status, ierr)
    if (ikeep .gt. 0) then
      nrecv = 2*ikeep
      call mpi_recv(RZkeep,nrecv, MPI_DOUBLE_PRECISION, j, j, MPI_COMM_WORLD, status, ierr)
      write(ivtk) ( (/RZkeep(1,i), RZkeep(2,i), ZERO /),i=1,ikeep)
    endif
  enddo
else
  call mpi_send(ikeep, 1, MPI_INTEGER, 0, my_id, MPI_COMM_WORLD, ierr)
  if (ikeep .gt. 0) then
    nsend = 2*ikeep
    call mpi_send(RZkeep, nsend, MPI_DOUBLE_PRECISION, 0, my_id, MPI_COMM_WORLD, ierr)
  endif
endif

if (my_id .eq. 0) then
  ! POINT_DATA SECTION
  write(str1(1:12),'(i12)') nnos
  buffer = lf//lf//'POINT_DATA '//str1//lf                                              ; write(ivtk) trim(buffer)
endif
!===========================================Temperature in keV
scalars(:,2) = scalars(:,2) / MU_zero / (central_density * 1d20) / 1.602d-19 /2.*1.e-3 !(assumes Te=Ti=T/2)
! ------- normalisation of psi
scalars(:,3) = (scalars(:,3) - ES%psi_axis ) / (psi_bnd - ES%psi_axis )
!=============================================
do i_var =1, n_scalars

  if (my_id .eq. 0) then
    buffer = 'SCALARS '//scalar_names(i_var)//' float'//lf                              ; write(ivtk) trim(buffer)
    buffer = 'LOOKUP_TABLE default'//lf                                                 ; write(ivtk) trim(buffer)
    write(ivtk) (scalars(i,i_var),i=1,ikeep0)
  endif

  if (my_id .eq. 0) then
    do j=1,n_cpu-1
      call mpi_recv(ikeep,1, MPI_INTEGER, j, j, MPI_COMM_WORLD, status, ierr)
      if (ikeep .gt. 0) then
        nrecv = ikeep
        call mpi_recv(scalars(1:ikeep,i_var),nrecv, MPI_DOUBLE_PRECISION, j, j, MPI_COMM_WORLD, status, ierr)
        write(ivtk) (scalars(i,i_var),i=1,ikeep)
      endif
    enddo
  else
    call mpi_send(ikeep, 1, MPI_INTEGER, 0, my_id, MPI_COMM_WORLD, ierr)
    if (ikeep .gt. 0) then
      nsend = ikeep
      call mpi_send(scalars(1:ikeep,i_var), nsend, MPI_DOUBLE_PRECISION, 0, my_id, MPI_COMM_WORLD, ierr)
    endif
  endif

enddo

close(ivtk)
if (my_id .eq. 0) then
   write(*,*) 'file connection.vtk written'
endif

deallocate(RZkeep,scalars,scalar_names)

!======================================== write strike point data to file
open(23,file='strikes_coordinates.txt')
open(24,file='strikes_values.txt')

do i=1,i_strike
  if (abs(R_strike(i)) .gt. 10.d0) R_strike(i) = 0.d0
  if (abs(Z_strike(i)) .gt. 10.d0) Z_strike(i) = 0.d0
!  if (PS0_strike(i) .gt. ES%psi_xpoint(1)) then       ! to exclude points started outside the plasma
!    R_strike(i) = 0.d0
!    Z_strike(i) = 0.d0
!  endif
enddo

i_strike0 = i_strike

call MPI_Reduce(i_strike,nnos,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)

if (my_id .eq. 0) write(*,*) ' number of points_i_strike : ',nnos

n_scalars = 3             ! number of scalars to write to the VTK output file

allocate(scalar_names(n_scalars),scalars(100000,n_scalars))

scalar_names  = (/ 'length_m    ','psi_start   ','T_start_keV '/)

do i=1,i_strike
!  write(*,*) i
!  write(*,*) i,R_strike(i), Z_strike(i),PS0_strike(i), T0_strike(i)
  scalars(i,1) = C_strike(i)                 ! needs correction !!! see above
  scalars(i,2) = PS0_strike(i)
  scalars(i,3) = T0_strike(i)
enddo

if (my_id .eq. 0) then
#ifdef IBM_MACHINE
  open(unit=ivtk,file='strikes.vtk',form='unformatted',access='stream')
#else
  open(unit=ivtk,file='strikes.vtk',form='unformatted',access='stream',convert='BIG_ENDIAN')
#endif

  buffer = '# vtk DataFile Version 3.0'//lf                                             ; write(ivtk) trim(buffer)
  buffer = 'vtk output'//lf                                                             ; write(ivtk) trim(buffer)
  buffer = 'BINARY'//lf                                                                 ; write(ivtk) trim(buffer)
  buffer = 'DATASET UNSTRUCTURED_GRID'//lf//lf                                          ; write(ivtk) trim(buffer)

  ! POINTS SECTION
  write(str1(1:12),'(i12)') nnos
  buffer = 'POINTS '//str1//'  float'//lf                                               ; write(ivtk) trim(buffer)
endif

if (my_id .eq. 0) then
  write(ivtk) ( (/R_strike(i)*cos(P_strike(i)), Z_strike(i), R_strike(i)*sin(P_strike(i)) /),i=1,i_strike0)
!  write(23,'(3e16.8)') ( (/ R_strike(i)*cos(P_strike(i)), Z_strike(i), R_strike(i)*sin(P_strike(i)) /),i=1,i_strike0)
   write(23,'(3e16.8)') ( (/R_strike(i), Z_strike(i), P_strike(i)/),i=1,i_strike0)
   write(24,'(3e16.8)') ( (/T0_strike(i),C_strike(i),PS0_strike(i)/),i=1,i_strike0)
endif

if (my_id .eq. 0) then
  do j=1,n_cpu-1
    call mpi_recv(i_strike,1, MPI_INTEGER, j, j, MPI_COMM_WORLD, status, ierr)
    if (i_strike .gt. 0) then
      nrecv = i_strike
      call mpi_recv(R_strike,nrecv, MPI_DOUBLE_PRECISION, j, j, MPI_COMM_WORLD, status, ierr)
      call mpi_recv(Z_strike,nrecv, MPI_DOUBLE_PRECISION, j, j, MPI_COMM_WORLD, status, ierr)
      call mpi_recv(P_strike,nrecv, MPI_DOUBLE_PRECISION, j, j, MPI_COMM_WORLD, status, ierr)
      write(ivtk) ( (/R_strike(i)*cos(P_strike(i)), Z_strike(i), R_strike(i)*sin(P_strike(i)) /),i=1,i_strike)
!      write(23,'(3e16.8)') ( (/R_strike(i)*cos(P_strike(i)), Z_strike(i), R_strike(i)*sin(P_strike(i)) /),i=1,i_strike)
    write(23,'(3e16.8)') ( (/R_strike(i), Z_strike(i), P_strike(i)/),i=1,i_strike0)
    write(24,'(3e16.8)') ( (/T0_strike(i),C_strike(i),PS0_strike(i)/),i=1,i_strike0)
    endif
  enddo
else
  call mpi_send(i_strike, 1, MPI_INTEGER, 0, my_id, MPI_COMM_WORLD, ierr)
  if (i_strike .gt. 0) then
    nsend = i_strike
    call mpi_send(R_strike, nsend, MPI_DOUBLE_PRECISION, 0, my_id, MPI_COMM_WORLD, ierr)
    call mpi_send(Z_strike, nsend, MPI_DOUBLE_PRECISION, 0, my_id, MPI_COMM_WORLD, ierr)
    call mpi_send(P_strike, nsend, MPI_DOUBLE_PRECISION, 0, my_id, MPI_COMM_WORLD, ierr)
  endif
endif

if (my_id .eq. 0) then
  ! POINT_DATA SECTION
  write(str1(1:12),'(i12)') nnos
  buffer = lf//lf//'POINT_DATA '//str1//lf                                              ; write(ivtk) trim(buffer)
endif
!===========================================Temperature in keV
scalars(:,3) = scalars(:,3) / MU_zero / (central_density * 1d20) / 1.602d-19 /2.*1.e-3 !(assumes Te=Ti=T/2)

do i_var =1, n_scalars

  if (my_id .eq. 0) then
    buffer = 'SCALARS '//scalar_names(i_var)//' float'//lf                              ; write(ivtk) trim(buffer)
    buffer = 'LOOKUP_TABLE default'//lf                                                 ; write(ivtk) trim(buffer)
    write(ivtk) (scalars(i,i_var),i=1,i_strike0)
  endif

  if (my_id .eq. 0) then
    do j=1,n_cpu-1
      call mpi_recv(i_strike,1, MPI_INTEGER, j, j, MPI_COMM_WORLD, status, ierr)
      if (i_strike .gt. 0) then
        nrecv = i_strike
        call mpi_recv(scalars(1:i_strike,i_var),nrecv, MPI_DOUBLE_PRECISION, j, j, MPI_COMM_WORLD, status, ierr)
        write(ivtk) (scalars(i,i_var),i=1,i_strike)
      endif
    enddo
  else
    call mpi_send(i_strike, 1, MPI_INTEGER, 0, my_id, MPI_COMM_WORLD, ierr)
    if (i_strike .gt. 0) then
      nsend = i_strike
      call mpi_send(scalars(1:i_strike,i_var), nsend, MPI_DOUBLE_PRECISION, 0, my_id, MPI_COMM_WORLD, ierr)
    endif
  endif

enddo

close(ivtk)
close(23)
if (my_id .eq. 0) then
   write(*,*) 'file strikes.vtk written'
endif

call MPI_FINALIZE(IERR)                                ! clean up MPI

end program jorek2_connection2

subroutine step(i_elm,s_in,t_in,p_in,delta_p,delta_s,delta_t,R,Z,R_s,R_t,Z_s,Z_t)
use mod_parameters
use elements_nodes_neighbours
use phys_module
use mod_interp

implicit none

integer :: i_var_psi, i_elm, i_tor, i_harm

real*8 :: s_in, t_in, p_in, delta_p, delta_s, delta_t
real*8 :: R_out, Z_out, Rs_out, Rt_out, Zs_out, Zt_out
real*8 :: R,R_s,R_t,R_st,R_ss,R_tt,Z,Z_s,Z_t,Z_st,Z_ss,Z_tt
real*8 :: Pcos,Pcos_s,Pcos_t,Pcos_st,Pcos_ss,Pcos_tt, Psin,Psin_s,Psin_t,Psin_st,Psin_ss,Psin_tt
real*8 :: P0,P0_s,P0_t,P0_st,P0_ss,P0_tt, psi_s, psi_t, Zjac

i_var_psi = 1

call interp_RZ(node_list,element_list,i_elm,s_in,t_in,R,R_s,R_t,Z,Z_s,Z_t)

Zjac = (R_s * Z_t - R_t * Z_s)

call interp(node_list,element_list,i_elm,i_var_psi,1,s_in,t_in,P0,P0_s,P0_t,P0_st,P0_ss,P0_tt)

psi_s = P0_s 
psi_t = P0_t 

do i_tor = 1, (n_tor-1)/2

  i_harm = 2*i_tor

  call interp(node_list,element_list,i_elm,i_var_psi,i_harm,s_in,t_in,Pcos,Pcos_s,Pcos_t,Pcos_st,Pcos_ss,Pcos_tt)

  psi_s = psi_s + Pcos_s * cos(mode(i_harm)*p_in)
  psi_t = psi_t + Pcos_t * cos(mode(i_harm)*p_in)

  call interp(node_list,element_list,i_elm,i_var_psi,i_harm+1,s_in,t_in,Psin,Psin_s,Psin_t,Psin_st,Psin_ss,Psin_tt)

  psi_s = psi_s + Psin_s * sin(mode(i_harm+1)*p_in)
  psi_t = psi_t + Psin_t * sin(mode(i_harm+1)*p_in)

enddo

delta_s =   psi_t * R / (Zjac * F0) * delta_p
delta_t = - psi_s * R / (Zjac * F0) * delta_p

return
end subroutine step

subroutine var_value(i_elm,i_var,s_in,t_in,p_in,value_out)
use mod_parameters
use elements_nodes_neighbours
use phys_module
use mod_interp, only: interp

implicit none

integer :: i_var, i_elm, i_tor, i_harm

real*8 :: s_in, t_in, p_in
real*8 :: Pcos,Pcos_s,Pcos_t,Pcos_st,Pcos_ss,Pcos_tt, Psin,Psin_s,Psin_t,Psin_st,Psin_ss,Psin_tt
real*8 :: P0,P0_s,P0_t,P0_st,P0_ss,P0_tt
real*8 :: value_out


!call interp_RZ(node_list,element_list,i_elm,s_in,t_in,R,R_s,R_t,R_st,R_ss,R_tt,Z,Z_s,Z_t,Z_st,Z_ss,Z_tt)
!Zjac = (R_s * Z_t - R_t * Z_s)

call interp(node_list,element_list,i_elm,i_var,1,s_in,t_in,P0,P0_s,P0_t,P0_st,P0_ss,P0_tt)

value_out = P0

do i_tor = 1, (n_tor-1)/2

  i_harm = 2*i_tor

  call interp(node_list,element_list,i_elm,i_var,i_harm,s_in,t_in,Pcos,Pcos_s,Pcos_t,Pcos_st,Pcos_ss,Pcos_tt)

  value_out = value_out + Pcos * cos(mode(i_harm)*p_in)

  call interp(node_list,element_list,i_elm,i_var,i_harm+1,s_in,t_in,Psin,Psin_s,Psin_t,Psin_st,Psin_ss,Psin_tt)

  value_out = value_out + Psin * sin(mode(i_harm+1)*p_in)

enddo

return
end subroutine var_value
