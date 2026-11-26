!> Compute q-profile for a stellarator JOREK restart file
program q_stel

use constants
use mod_parameters
use data_structure
use basis_at_gaussian
use nodes_elements, only:bnd_node_list,bnd_elm_list
use equil_info
use phys_module
use mod_import_restart
use mod_log_params
use mod_interp
use mod_chi
use elements_nodes_neighbours
use mod_neighbours
use mod_find_rz_nearby
use mpi_mod
use mod_element_rtree, only: populate_element_rtree

implicit none

!--- Input parameters --------------------!
!-----------------------------------------!
integer, parameter :: points_per_turn = 2500
real*8, parameter  :: delta_phi = 2*PI/float(n_period*points_per_turn)
integer, parameter :: n_lines = 20
integer, parameter :: num_pol_turns = 200
integer, parameter :: assumed_max_q = 8
!-----------------------------------------!
!-----------------------------------------!

real*8    :: Rstart(n_lines), Zstart(n_lines)                      
integer   :: i, j, k, i_elm, ifail, my_id, ierr, inode, i_mpi            
integer   :: iside_i, iside_j
real*8    :: Rout, Zout, polturns, torturns
logical   :: stop_tracing
real*8    :: R_axis, Z_axis, R_max

real*8    :: s, t, t_global, phi, phinew, phiold, snew, tnew, sold, told
real*8    :: RR, R_s, R_t, R_p, ZZ, Z_s, Z_t, Z_p
real*8    :: BR, BZ, Bp, xjac, Fprof
real*8    :: dummy, dum01, dum02, dum03, dum04, dum05
real*8    :: Rold, Zold

integer   :: nsend, nrecv
integer   :: required,provided,StatInfo, n_mpi
integer*4 :: rank, comm_size 
logical   :: responsible(n_lines)
integer   :: status(MPI_STATUS_SIZE)

real*8    :: rphin_arr(n_lines) = 0.d0, polturns_arr(n_lines) = 0.d0, torturns_arr(n_lines) = 0.d0, phi_arr(n_lines) = 0.d0, R_arr(n_lines) = 0.d0
real*8    :: rphin_arr_tot(n_lines) = 0.d0, polturns_arr_tot(n_lines) = 0.0, torturns_arr_tot(n_lines) = 0.0, phi_arr_tot(n_lines) = 0.d0, R_arr_tot(n_lines) = 0.d0
real*8    :: R_poinc_tot(n_lines*num_pol_turns*n_period*assumed_max_q) = 0.d0, Z_poinc_tot(n_lines*num_pol_turns*n_period*assumed_max_q), phi_poinc_tot(n_lines*num_pol_turns*n_period*assumed_max_q)

! --- Initialise constants
integer   :: v_s0_t0   = 1    ! the vertex and edge indices follow and anti-clockwise convention
integer   :: v_s1_t0   = 2    ! around the element: https://www.jorek.eu/wiki/doku.php?id=grids
integer   :: v_s1_t1   = 3
integer   :: v_s0_t1   = 4
integer   :: e_splus   = 2
integer   :: e_sminus  = 4
integer   :: e_tplus   = 3
integer   :: e_tminus  = 1

integer   :: i_var_psi = 1

required = MPI_THREAD_FUNNELED
call MPI_Init_thread(required, provided, StatInfo)
call init_threads()  ! on some systems init_threads needs to come after mpi_init_thread
call MPI_COMM_SIZE(MPI_COMM_WORLD, comm_size, ierr)
n_mpi = comm_size

! --- Determine ID of each MPI proc
call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
my_id = rank

if ( my_id == 0 ) then
  write(*,*) '***************************************'
  write(*,*) '* Calculate q-profile...              *'
  write(*,*) '***************************************'
end if

! Read and broadcast restart file
call det_modes()
call initialise_basis
call init_chi_basis
call initialise_parameters(my_id,  "__NO_FILENAME__")
call log_parameters(my_id)
if ( my_id == 0 ) call import_restart(node_list,element_list, 'jorek_restart', rst_format, ierr, .true.)
call broadcast_phys(my_id)  
call broadcast_elements(my_id, element_list)                ! elements
call broadcast_nodes(my_id, node_list)                      ! nodes
call broadcast_boundary(my_id,bnd_elm_list,bnd_node_list)
call populate_element_rtree(node_list, element_list)
call broadcast_equil_state(my_id)

! Determine element neighbours
allocate(element_neighbours(4,element_list%n_elements))
element_neighbours = 0
! --- Neighbour search is O(N) with a bounding-box grid (see build_element_neighbours in mod_neighbours)
call build_element_neighbours(node_list, element_list, element_neighbours)

! Get R_axis, Z_axis and R_max - WARNING: this initialisation assumes a flux aligned grid ordering
phi = 2.d0*pi*float(i_plane_rtree - 1)/float(n_period*n_plane)
call interp_RZP(node_list,element_list,1,0.0,0.0,phi,R_axis,Z_axis)
call interp_RZP(node_list,element_list,element_list%n_elements,1.0,1.0,phi,R_max,dummy)
R_max = (R_max-R_axis)/bloating_factor + R_axis   ! This will calculate approximately up to the LCFS.
if ( my_id == 0 ) then
  write(*,*) '*** ...start tracing... ***'
  write(*,*) 'delta_phi        = ', delta_phi
  write(*,*) 'num_pts          = ', n_lines
  write(*,*) 'num_pol_turns    = ', num_pol_turns
  write(*,*) '(R_axis, Z_axis) = ', R_axis, Z_axis
  write(*,*) '(R_max, Z_max)   = ', R_max, dummy
#ifdef POINC_GVEC
  write(*,*) 'Using GVEC field for tracing'
#endif
endif

! --- Define starting points for field lines along midplane
responsible = .false.
do i = 1, n_lines
  
  if ( ( real(my_id)/real(n_mpi)*n_lines < i ) .and. ( real(my_id+1)/real(n_mpi)*n_lines >= i ) ) then
    responsible(i) = .true.
    Rstart(i) = R_axis + float(i) * 0.99 * (R_max - R_axis)/n_lines
    Zstart(i) = Z_axis
  end if
  
end do

! --- Write out distribution of field lines among tasks
do i = 0, n_mpi-1
  if ( my_id == i ) then
    write(*,*) 'Task ', my_id, 'responsible for (field line number & radius):'
    do j = 1, n_lines
      if (responsible(j)) write(*,*) j, Rstart(j), Zstart(j)
    end do
  end if
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
end do

! Loop through starting points, identifying which field lines MPI process is responsible for
R_poinc_tot = 0.0; Z_poinc_tot = 0.0; phi_poinc_tot = 0.0
do i = 1, n_lines
  if ( .not. responsible(i) ) then
    rphin_arr(i) = 0.d0
    phi_arr(i)  = 0.d0
    R_arr(i)  = 0.d0
    cycle
  end if
  
  ! Get starting point R, Z, phi, s and t
  phi = 2.d0*pi*float(i_plane_rtree - 1)/float(n_period*n_plane)
  t_global = 0.0
  call find_RZ(node_list,element_list,Rstart(i),Zstart(i),RR,ZZ,i_elm,s,t,ifail)
  if (ifail .ne. 0) then
    write(*,*) "Can not find RZ,", ifail 
    exit
  endif
  R_arr(i)        = RR
  
  ! Get radial coordinate from current element and local s coordinate
  call interp_gvec(node_list,element_list,i_elm,4,1,1,s,t,rphin_arr(i),dummy,dummy,dummy,dummy,dummy)

  if ( i_elm < 1 ) then
    write(*,*) 'Illegal starting point for field line ', i, '. Skipping.'
    cycle
  end if
  
  ! --- Trace field line
  stop_tracing = .false.
  polturns     = 0.d0; torturns     = 0.d0
  j            = 0
  ifail = 0
  do while( .not. stop_tracing )
    call do_step()

    ! Check if tracing failed and exit 
    if (ifail .ne. 0) then
       write(*,*) 'Field line tracing failed for line: ', i
       exit
    endif
    j = j + 1
    
    ! Check the assumed maximum number of toroidal turns is not exceeded
    if (abs(phi) / (2.d0*PI) .gt. assumed_max_q*num_pol_turns) then
      write(*, *) "ERROR: Assumed maximum q has been exceeded for line", i
      stop_tracing = .true.
    endif
    
    ! Include points when single field period is crossed in Poincare
    if ( mod(j, points_per_turn) .eq. 0 ) then
      R_poinc_tot((i-1)*num_pol_turns*n_period*assumed_max_q+j/points_per_turn) = RR
      Z_poinc_tot((i-1)*num_pol_turns*n_period*assumed_max_q+j/points_per_turn) = ZZ
      phi_poinc_tot((i-1)*num_pol_turns*n_period*assumed_max_q+j/points_per_turn) = phi
    endif

    ! Determine if poloidal turn is made
    if (abs(t_global) > n_tht) then
      t_global = t_global - sign(float(n_tht), t_global)
      polturns = polturns + 1.d0
      if (polturns .eq. num_pol_turns) stop_tracing = .true.
    endif

  end do

  ! Get total toroidal distance traced
  phi_arr(i)      = phi - 2.d0*pi*float(i_plane_rtree - 1)/float(n_period*n_plane)
  torturns_arr(i) = phi_arr(i)/float(n_period*n_plane)/(2.d0*PI)
  polturns_arr(i) = polturns + t_global / float(n_tht) 
  write(*,*) 'Finished tracing line ', i

end do
call MPI_Barrier(MPI_COMM_WORLD,ierr)

! Fill output arrays for q profile data
call MPI_Reduce(rphin_arr,  rphin_arr_tot,  n_lines,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
call MPI_Reduce(phi_arr,   phi_arr_tot,   n_lines,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
call MPI_Reduce(polturns_arr,   polturns_arr_tot,   n_lines,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
call MPI_Reduce(torturns_arr,   torturns_arr_tot,   n_lines,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
call MPI_Reduce(R_arr,     R_arr_tot,     n_lines,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

! --- Open the output files to which the Poincare data will be written in ascii format
if (my_id .eq. 0) then
  open(21,file='q_stel_poinc_R-Z.dat')
  write(21,*) '#  R       Z       % Completion        Safety Factor       Field Line No'
  
  ! --- Write points for local MPI (id=0)
  do i = 1, n_lines
    do j = 1, num_pol_turns*n_period*assumed_max_q 
      if (R_poinc_tot((i-1)*num_pol_turns*n_period*assumed_max_q+j) .ne. 0.0) then
        write(21,'(4e18.8,i6)') R_poinc_tot((i-1)*num_pol_turns*n_period*assumed_max_q+j), Z_poinc_tot((i-1)*num_pol_turns*n_period*assumed_max_q+j), 100.0 * polturns_arr_tot(i) / float(num_pol_turns), phi_arr_tot(i)/(2.d0*PI)/polturns_arr_tot(i), int(i)
      endif
    enddo
    write(21,*)
    write(21,*)
  enddo
  ! --- Write points for all other MPIs
  ! --- If this is mpi_0, we receive data from the other MPIs and print it
  do i_mpi=1,n_mpi-1
    nrecv = n_lines*num_pol_turns*n_period*assumed_max_q
    call mpi_recv(R_poinc_tot,nrecv, MPI_DOUBLE_PRECISION, i_mpi, i_mpi, MPI_COMM_WORLD, status, ierr)
    call mpi_recv(Z_poinc_tot,nrecv, MPI_DOUBLE_PRECISION, i_mpi, i_mpi, MPI_COMM_WORLD, status, ierr)
    call mpi_recv(phi_poinc_tot,nrecv, MPI_DOUBLE_PRECISION, i_mpi, i_mpi, MPI_COMM_WORLD, status, ierr)
    do i = 1, n_lines
      do j = 1, num_pol_turns*n_period*assumed_max_q 
        if (R_poinc_tot((i-1)*num_pol_turns*n_period*assumed_max_q+j) .ne. 0.0) then
          write(21,'(4e18.8,i6)') R_poinc_tot((i-1)*num_pol_turns*n_period*assumed_max_q+j), Z_poinc_tot((i-1)*num_pol_turns*n_period*assumed_max_q+j), 100.0 * polturns_arr_tot(i) / float(num_pol_turns), phi_arr_tot(i)/(2.d0*PI)/polturns_arr_tot(i), int(i)
        endif
      enddo
      write(21,*)
      write(21,*)
    enddo
  enddo
  close(21)
else
  ! --- If this is not mpi_0, we send data to the main MPI 0
  nsend = n_lines*num_pol_turns*n_period*assumed_max_q
  call mpi_send(R_poinc_tot, nsend,MPI_DOUBLE_PRECISION, 0, my_id, MPI_COMM_WORLD, ierr)
  call mpi_send(Z_poinc_tot, nsend,MPI_DOUBLE_PRECISION, 0, my_id, MPI_COMM_WORLD, ierr)
  call mpi_send(phi_poinc_tot, nsend,MPI_DOUBLE_PRECISION, 0, my_id, MPI_COMM_WORLD, ierr)
endif

! --- Write out q profile result
if ( my_id == 0 ) then
  open(99, file='q_profile.dat', action='write', status='replace')
  write(99,*) '# SQRT{Phi_N} | Safety factor | Radius at LFS midplane [m] | % Completion '
  do i = 1, n_lines
    write(99,'(4es23.5)') rphin_arr_tot(i), phi_arr_tot(i)/(2.d0*PI)/polturns_arr_tot(i), R_arr_tot(i), 100 * polturns_arr_tot(i) / float(num_pol_turns)
  end do
  close(99)
end if

call MPI_FINALIZE(IERR)                                ! clean up MPI

contains

  !> Field line tracing using simple half-step scheme
  subroutine do_step()
  
  integer    :: i_steps   ! Number of iterations carried out in performing step
  real*8     :: delta_s, delta_t, delta_phi_local, delta_phi_step, small_delta, small_delta_s, small_delta_t
  real*8     :: s_mid, t_mid, p_mid, R_mid, Z_mid
  logical    :: debug

  ! Store old location
  Rold   = RR
  Zold   = ZZ
  phiold = phi

  ! --- Loop within one element
  delta_phi_local = 0.d0
  i_steps = 0
  do while ((abs(delta_phi_local) .lt. abs(delta_phi)) .and. (i_steps .le. 10) )
    ! --- Count element steps (lt 10)
    i_steps = i_steps + 1

    ! --- Find the new position after element step
    delta_phi_step = delta_phi - delta_phi_local
    call determine_step_in_element(delta_phi_step, s, t, phi, delta_s, delta_t)
  
    ! --- Take another step from middle point (ie. 1.5*step)
    s_mid = s + 0.5d0 * delta_s
    t_mid = t + 0.5d0 * delta_t
    p_mid = phi + 0.5d0 * delta_phi_step
    call determine_step_in_element(delta_phi_step, s_mid, t_mid, p_mid, delta_s, delta_t)
  
    ! --- Step to element boundary, not beyond (we do smaller steps to approach boundary)
    small_delta_s = 1.d0
    if (s + delta_s .gt. 1.d0) small_delta_s = (1.d0 - s)/delta_s
    if (s + delta_s .lt. 0.d0) small_delta_s = abs(s/delta_s)
  
    ! --- Step to element boundary, not beyond (we do smaller steps to approach boundary)
    small_delta_t = 1.d0
    if (t + delta_t .gt. 1.d0) small_delta_t = (1.d0 - t)/delta_t
    if (t + delta_t .lt. 0.d0) small_delta_t = abs(t/delta_t)
  
    ! --- Do we require a smaller step (< 1.0) 
    small_delta = min(small_delta_s, small_delta_t)
    debug = .false.
    if (small_delta .le. 0) then
      write(*,*) 'ERROR: Something went wrong when calculating small_delta: ', i, i_steps, small_delta, small_delta_s, small_delta_t, s, t, delta_s, delta_t
      debug = .true.
    endif  

    if (small_delta .lt. 1.d0)  then       ! this step is crossing the boundary
      s_mid = s + 0.5d0 * small_delta * delta_s
      t_mid = t + 0.5d0 * small_delta * delta_t
      p_mid = phi + 0.5d0 * small_delta * delta_phi_step
  
      call determine_step_in_element(delta_phi_step, s_mid, t_mid, p_mid, delta_s, delta_t)
  
      if (small_delta_s .lt. small_delta_t) then
  
        if (s + delta_s .gt. 1.d0) then     ! crossing boundary 2 or 4 at s=1
          s = 1.0
          t_global = t_global + small_delta * delta_t
          call find_new_element(.true., i_elm, &
                                s, t, phi, small_delta, delta_s, delta_t, delta_phi_step, &
                                e_splus, e_sminus, v_s1_t0, v_s1_t1)
        elseif (s + delta_s .lt. 0.d0) then ! crossing boundary 2 or 4 at s=0
          s = 0.d0
          t_global = t_global + small_delta * delta_t
          call find_new_element(.true., i_elm, &
                                s, t, phi, small_delta, delta_s, delta_t, delta_phi_step, &
                                e_sminus, e_splus, v_s0_t1, v_s0_t0)
        else !if (debug) then
          write(*,*) 'ERROR: It should not be possible to get here!', s, delta_s, s+delta_s
        endif
      else
        if (t + delta_t .gt. 1.d0) then  ! crossing boundary 1 or 3 at t=1
          t = 1.0
          t_global = t_global + small_delta * delta_t
          call find_new_element(.false., i_elm, &
                                s, t, phi, small_delta, delta_s, delta_t, delta_phi_step, &
                                e_tplus, e_tminus, v_s1_t1, v_s0_t1)
        elseif (t + delta_t .lt. 0.d0) then  ! crossing boundary 1 or 3 at t=0
          t = 0.d0
          t_global = t_global + small_delta * delta_t
          call find_new_element(.false., i_elm, &
                                s, t, phi, small_delta, delta_s, delta_t, delta_phi_step, &
                                e_tminus, e_tplus, v_s0_t0, v_s0_t1)
        else !if (debug) then
          write(*,*) 'ERROR: It should not be possible to get here!', t, delta_t, t+delta_t
        endif
      endif
    else  ! this step remains within the element
      s = s + delta_s
      t = t + delta_t
      t_global = t_global + delta_t
      phi = phi + delta_phi_step
  
      small_delta = 1.d0
  
      call interp_RZP(node_list,element_list,i_elm,s,t,phi,RR,ZZ)
    endif

    ! --- Reset the toroial step
    delta_phi_local = delta_phi_local + small_delta * delta_phi_step

    ! --- Exit if we stepped out of domain
    if (i_elm .eq. 0) then
      write(*,*) "ERROR: New element found in do_step is outside the domain"
      stop
    endif

  enddo ! end loop for single step
  
  if (i_steps .gt. 10) then
    write(*,*) "WARNING: Maximum number of iterations exceeded in do_step!", i, j
    ifail = 1
  endif

  end subroutine do_step
 

 !> Determine change in s and t for a given step size
  subroutine determine_step_in_element(delta_p, s_in, t_in, p_in, delta_s, delta_t)
  use mod_parameters
  use elements_nodes_neighbours
  use phys_module
  use mod_interp
  use mod_chi
  implicit none
  
  real*8, intent(in)     :: delta_p, s_in, t_in, p_in
  real*8, intent(inout)  :: delta_s, delta_t
  
  integer :: i_tor, i_harm
  real*8 :: Pcos,Pcos_s,Pcos_t,Pcos_st,Pcos_ss,Pcos_tt, Psin,Psin_s,Psin_t,Psin_st,Psin_ss,Psin_tt
  real*8 :: P0,P0_s,P0_t,P0_st,P0_ss,P0_tt
  real*8 :: psi_s, psi_t, psi_R, psi_z, psi_p, st_psi_p
  real*8 :: BR, BZ, BP, BR0cos,BR0sin,BZ0cos,BZ0sin,Bp0cos,Bp0sin
  real*8 :: Zjac
  real*8, dimension(0:n_order-1,0:n_order-1,0:n_order-1) :: chi
  
  ! Get current location of field line and necessary derivatives
  call interp_RZP(node_list,element_list,i_elm,s_in,t_in,p_in,RR,R_s,R_t,R_p,dummy,dummy,dummy,dummy,dummy,dummy, &
                                                              ZZ,Z_s,Z_t,Z_p,dummy,dummy,dummy,dummy,dummy,dummy)
  chi  = get_chi(RR,ZZ,p_in,node_list,element_list,i_elm,s_in,t_in)
  Zjac = (R_s * Z_t - R_t * Z_s)

#ifdef POINC_GVEC
  call interp_gvec(node_list,element_list,i_elm,1,1,1,s_in,t_in,BR,dummy,dummy,dummy,dummy,dummy)
  call interp_gvec(node_list,element_list,i_elm,1,2,1,s_in,t_in,BZ,dummy,dummy,dummy,dummy,dummy)
  call interp_gvec(node_list,element_list,i_elm,1,3,1,s_in,t_in,Bp,dummy,dummy,dummy,dummy,dummy)
  do i_tor=1,(n_coord_tor-1)/2
    i_harm = 2*i_tor
    
    call interp_gvec(node_list,element_list,i_elm,1,1,i_harm,s_in,t_in,BR0cos,dummy,dummy,dummy,dummy,dummy)
    call interp_gvec(node_list,element_list,i_elm,1,2,i_harm,s_in,t_in,BZ0cos,dummy,dummy,dummy,dummy,dummy)
    call interp_gvec(node_list,element_list,i_elm,1,3,i_harm,s_in,t_in,Bp0cos,dummy,dummy,dummy,dummy,dummy)
    
    BR = BR + BR0cos*cos(mode_coord(i_harm)*p_in)
    BZ = BZ + BZ0cos*cos(mode_coord(i_harm)*p_in)
    Bp = Bp + Bp0cos*cos(mode_coord(i_harm)*p_in)
    
    call interp_gvec(node_list,element_list,i_elm,1,1,i_harm+1,s_in,t_in,BR0sin,dummy,dummy,dummy,dummy,dummy)
    call interp_gvec(node_list,element_list,i_elm,1,2,i_harm+1,s_in,t_in,BZ0sin,dummy,dummy,dummy,dummy,dummy)
    call interp_gvec(node_list,element_list,i_elm,1,3,i_harm+1,s_in,t_in,Bp0sin,dummy,dummy,dummy,dummy,dummy)
    
    BR = BR - BR0sin*sin(mode_coord(i_harm+1)*p_in)
    BZ = BZ - BZ0sin*sin(mode_coord(i_harm+1)*p_in)
    Bp = Bp - Bp0sin*sin(mode_coord(i_harm+1)*p_in)
  end do
#else
  ! Get n=0 component of Psi and derivatives
  call interp(node_list,element_list,i_elm,i_var_psi,1,s_in,t_in,P0,P0_s,P0_t,P0_st,P0_ss,P0_tt)
  psi_s = P0_s 
  psi_t = P0_t 
  st_psi_p = 0.d0

  ! Get higher toroidal harmonics of Psi and derivatives
  do i_tor = 1, (n_tor-1)/2
  
    i_harm = 2*i_tor
    call interp(node_list,element_list,i_elm,i_var_psi,i_harm,s_in,t_in,Pcos,Pcos_s,Pcos_t,Pcos_st,Pcos_ss,Pcos_tt)
    psi_s = psi_s + Pcos_s * cos(mode(i_harm)*p_in)
    psi_t = psi_t + Pcos_t * cos(mode(i_harm)*p_in)
    st_psi_p = st_psi_p - Pcos*mode(i_harm)*sin(mode(i_harm)*p_in)
  
    call interp(node_list,element_list,i_elm,i_var_psi,i_harm+1,s_in,t_in,Psin,Psin_s,Psin_t,Psin_st,Psin_ss,Psin_tt)
    psi_s = psi_s + Psin_s * sin(mode(i_harm+1)*p_in)
    psi_t = psi_t + Psin_t * sin(mode(i_harm+1)*p_in)
    st_psi_p = st_psi_p + Psin*mode(i_harm+1)*cos(mode(i_harm+1)*p_in)
  enddo
  
  ! Get magnetic field components
  psi_R = ( Z_t*psi_s - Z_s*psi_t)/Zjac
  psi_z = (-R_t*psi_s + R_s*psi_t)/Zjac
  psi_p = st_psi_p - R_p*psi_R - Z_p*psi_z
  BR = chi(1,0,0)    + (psi_z*chi(0,0,1) - psi_p*chi(0,1,0))/(F0*RR)
  BZ = chi(0,1,0)    - (psi_R*chi(0,0,1) - psi_p*chi(1,0,0))/(F0*RR) 
  Bp = chi(0,0,1)/RR + (psi_R*chi(0,1,0) - psi_z*chi(1,0,0))/F0     
#endif

  ! Determine distance moved in s-t space
  delta_s = (-Z_t*R_p + R_t*Z_p + RR*(Z_t*BR - R_t*BZ)/Bp)*delta_p/Zjac
  delta_t = ( Z_s*R_p - R_s*Z_p - RR*(Z_s*BR - R_s*BZ)/Bp)*delta_p/Zjac

  end subroutine determine_step_in_element


  subroutine find_new_element(s_bnd, i_elm,  &
                            s_line, t_line, p_line, small_delta, delta_s, delta_t, delta_phi_step, &
                            element_to_neighbour_idx, neighbour_to_element_idx, vertex_1, vertex_2)
  use mod_parameters
  use elements_nodes_neighbours
  use phys_module
  use mod_interp

  implicit none
  ! Determines the new element after performing a step along a field line that has crossed a boundary
  ! 
  ! 1. Find boundary point in original element 
  ! 2. Get neighbour element
  ! 3. Confirm boundary point is the same for both elements
  ! 4. If new element is outside domain, check the element edge being crossed is on the boundary
  
  ! Input parameters
  logical, intent(in)     :: s_bnd                                                                    ! Determine whether crossing is from the s coordinate  
  integer, intent(inout)  :: i_elm                                                                    ! Element index of current element
  integer, intent(in)     :: element_to_neighbour_idx, neighbour_to_element_idx                       ! Indices of current element to neighbour and expected neighbour to current element without orientation change
  integer, intent(in)     :: vertex_1, vertex_2                                                       ! Vertex indices for boundary nodes
  real*8, intent(inout)   :: s_line, t_line, p_line, small_delta, delta_s, delta_t, delta_phi_step    ! Current values of s, t, and phi and deltas from current step

  ! Local variables
  integer                 :: i_elm_prev                                                               ! Element index of previous element
  real*8                  :: R_in, R_out                                                              ! Variables for interpolation of R
  real*8                  :: Z_in, Z_out                                                              ! Variables for interpolation of Z
  integer                 :: i_elm_tmp                                                                ! Temporary element to check neighbour element orientation is consistent
  integer                 :: inode1, inode2                                                           ! Nodes corresponding to boundary edge

  ! Get new location on element boundary
  if (s_bnd) then
    t_line = t_line + small_delta * delta_t
  else
    s_line = s_line + small_delta * delta_s
  endif
  p_line = p_line + small_delta * delta_phi_step
  call interp_RZP(node_list,element_list,i_elm,s_line,t_line,p_line,R_in,Z_in)

  ! Find element neighbour
  i_elm_prev = i_elm
  i_elm = element_neighbours(element_to_neighbour_idx,i_elm_prev)
  
  if (i_elm .ne. 0) then ! check new element neighbour is consistent and update s and t
    i_elm_tmp  = element_neighbours(neighbour_to_element_idx,i_elm)
    if (i_elm_prev .ne. i_elm_tmp) write(*,*) ' WARNING : CHANGE OF ORIENTATION'

    ! Switch boundary coordinate from 1.0 to 0.0 or vice versa
    if (s_bnd) then
      s_line = abs(1.d0 - s_line)
    else
      t_line = abs(1.d0 - t_line)
    endif

    ! Check if R and Z from both elements are close
    call interp_RZP(node_list,element_list,i_elm,s_line,t_line,p_line,R_out,Z_out)
    if ( (abs(R_in - R_out) .gt. 1.d-8) .or. (abs(Z_in - Z_out) .gt. 1.d-8)) then
      write(*,'(A,2i6,4f12.4)') ' error in element change ',i_elm_prev,i_elm,R_in,R_out,Z_in,Z_out
    endif

  else ! crossing an outer boundary
    inode1 = element_list%element(i_elm_prev)%vertex(vertex_1)
    inode2 = element_list%element(i_elm_prev)%vertex(vertex_2)

    if ((node_list%node(inode1)%boundary .ne. 0) .and. (node_list%node(inode2)%boundary .ne. 0)) then
      write(*,*) 'Crossing an outer boundary'
    else
      write(*,*) 'ERROR : leaving domain but not at correct boundary!',inode1, inode2, R_in, Z_in
      stop
    endif
  endif
end


end program q_stel
