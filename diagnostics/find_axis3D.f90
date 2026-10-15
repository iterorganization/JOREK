program find_axis3D
!---------------------------------------------------------------------
! This routine attempts to find the location of the magnetic axis in a
! given poloidal plane of a 3D configuration (not restricted to stellarators), 
! by field line tracing. The steps of the algorithm are:
!  1) Guess (provided as input) a point for the magnetic axis in a poloidal plane
!  2) Iterate over the following procedure n_iter times:
!    1) Trace field line from guessed point for n_turn toroidal turns
!    2) Check if the error (1/10 of geometric mean max(Rp) - min(Rp) and max(zp) - min(zp),
!       where {(Rp(i),zp(i))} are the points where field line intersects the poloidal plane)
!       is less than the tolerance
!    3) If so, write out axis location, and exit; else iterate again starting from
!       ((max(Rp)+min(Rp))/2,(max(zp)+min(zp))/2)
!  3) If n_iter is exceeded, terminate and write out best approximation
!---------------------------------------------------------------------
  use data_structure
  use phys_module
  use basis_at_gaussian
  use elements_nodes_neighbours
  use mod_neighbours
  use mod_import_restart
  use mod_interp
  use mod_chi
  implicit none
  
  integer :: n_turn, n_iter, ierr, i, j, iside_i, iside_j, ielm, ifail, i_turn, n_phi, i_phi, i_steps, ielm_prev, ielm_tmp, ip, checked_elms
  real*8  :: R_init, z_init, P_init, R_out, z_out, s_out, t_out, R_line, z_line, p_line, s_line, t_line, delta_phi, delta_phi_local, tol
  real*8  :: delta_phi_step, delta_s, delta_t, s_mid, t_mid, p_mid, small_delta_s, small_delta_t, small_delta, R_in, z_in, R, z
  real*8  :: Rc, zc, error, R_start, z_start, P_start, mix
  logical :: pos_out, success
  
  real*8, dimension(:), allocatable :: Rp, zp
  
  character(:), allocatable :: filename, outfile
  
  namelist /find_axis_params/ filename, outfile, R_init, z_init, P_init, n_turn, n_iter, mix, tol, pos_out
  call det_modes
  call initialise_basis
  call init_chi_basis
  call initialise_parameters(0,  "__NO_FILENAME__")
  
  filename = "jorek_restart"
  outfile = "axis_position.dat"
  R_init = R_geo
  z_init = z_geo
  P_init = 2.d0*pi*float(i_plane_rtree - 1)/float(n_period*n_plane)
  n_turn = 10
  n_iter = 10
  mix    = 0.5
  tol    = 1.d-6
  pos_out = .true.
  
  open(11, file="find_axis.nml", action='read', status='old', iostat=ierr)
  if (ierr .eq. 0) then
    read(11, find_axis_params)
    close(11)
  end if
  
  if (mix .gt. 1.d0 .or. mix .le. 0.d0) mix = 0.5
  
  write(*,*) "***************************************"
  write(*,*) "*             find_axis3D             *"
  write(*,*) "***************************************"
  write(*,*) " Searching for axis in file ", filename
  write(*,*) " Output file: ", outfile
  write(*,'(A,2E14.6)') "  Initial guess: ", R_init, z_init
  write(*,*) " On poloidal plane ", P_init
  write(*,*) " Number of toroidal turns to be used: ", n_turn
  write(*,*) " Number of iterations: ", n_iter
  write(*,'(A,E14.6)') "  Mixing factor: ", mix
  write(*,'(A,E14.6)') "  Tolerance: ", tol
  write(*,*) " Output field line trace results for each iteration: ", pos_out
  write(*,*)
  
  call import_restart(node_list,element_list, filename, rst_format, ierr, .true.)
  
  allocate(element_neighbours(4,element_list%n_elements))
  element_neighbours = 0
  ! --- Neighbour search is O(N) with a bounding-box grid (see build_element_neighbours in mod_neighbours)
  call build_element_neighbours(node_list, element_list, element_neighbours)
  
  n_phi     = 1500
  delta_phi = 2.d0*pi/float(n_period*n_phi)
  allocate(Rp(n_turn),zp(n_turn))
  R_start = R_init
  z_start = z_init
  P_start = P_init
  success = .false.
  
  if (pos_out) open(13,file="position_R-z.dat")
  ! Begin iterations
  L_IT: do i=1,n_iter
    ip = 0
    call find_RZP(node_list,element_list,R_start,z_start,P_start,R_out,z_out,ielm,s_out,t_out,ifail,checked_elms)
    if (ifail .ne. 0) then
      write(*,*) "Can not find RZ,", ifail, R_start, z_start, P_start
      stop 1
    end if
    
    write(*,'(A,I4,A,2E14.6)') "Starting iteration ", i, "; position: ", R_start, z_start
  
    R_line = R_start
    z_line = z_start
    p_line = P_start
    s_line = s_out
    t_line = t_out
  
    ! Complete n_turn toroidal turns around the device
    do i_turn=1,n_turn
      do i_phi=1,n_phi
        delta_phi_local = 0.d0
        i_steps = 0
      
        do while ((delta_phi_local .lt. delta_phi) .and. (i_steps .lt.10))
          i_steps = i_steps + 1
          delta_phi_step = delta_phi - delta_phi_local
        
          call step(ielm,s_line,t_line,p_line,delta_phi_step,delta_s,delta_t)
          s_mid = s_line + 0.5d0*delta_s
          t_mid = t_line + 0.5d0*delta_t
          p_mid = p_line + 0.5d0*delta_phi_step
          call step(ielm,s_mid,t_mid,p_mid,delta_phi_step,delta_s,delta_t)
        
          small_delta_s = 1.d0
          if (s_line + delta_s .gt. 1.d0) then
            small_delta_s = (1.d0 - s_line)/delta_s
          else if (s_line + delta_s .lt. 0.d0) then
            small_delta_s = abs(s_line/delta_s)
          end if
        
          small_delta_t = 1.d0
          if (t_line + delta_t .gt. 1.d0) then
            small_delta_t = (1.d0 - t_line)/delta_t
          else if (t_line + delta_t .lt. 0.d0) then
            small_delta_t = abs(t_line/delta_t)
          end if
        
          small_delta = min(small_delta_s,small_delta_t)
          if (small_delta .lt. 1.d0) then 
            s_mid = s_line + 0.5d0*small_delta*delta_s
            t_mid = t_line + 0.5d0*small_delta*delta_t
            p_mid = p_line + 0.5d0*small_delta*delta_phi_step
            call step(ielm,s_mid,t_mid,p_mid,delta_phi_step,delta_s,delta_t)
          
            if (small_delta_s .lt. small_delta_t) then
              if (s_line + delta_s .gt. 1.d0) then
                s_line = 1.d0
                t_line = t_line + small_delta*delta_t
                p_line = p_line + small_delta*delta_phi_step
                call interp_RZP(node_list,element_list,ielm,s_line,t_line,p_line,R_in,z_in)

                ielm_prev = ielm      
                ielm      = element_neighbours(2,ielm_prev)
                ielm_tmp  = element_neighbours(4,ielm)
	      
                if (ielm_prev .ne. ielm_tmp) write(*,*) ' WARNING : CHANGE OF ORIENTATION (1)'
	
                s_line = 0.d0
                call interp_RZP(node_list,element_list,ielm,s_line,t_line,p_line,R_out,z_out)
              
                if ((abs(R_in - R_out) .gt. 1.d-8) .or. (abs(z_in - z_out) .gt. 1.d-8)) &
                  write(*,'(A,2i6,4f8.4)') ' error in element change (1) ',ielm_prev,ielm,R_in,R_out,z_in,z_out
              else if (s_line + delta_s .lt. 0.d0) then
                s_line = 0.d0
                t_line = t_line + small_delta*delta_t
                p_line = p_line + small_delta*delta_phi_step
                call interp_RZP(node_list,element_list,ielm,s_line,t_line,p_line,R_in,z_in)
              
                ielm_prev = ielm      
                ielm      = element_neighbours(4,ielm_prev)
                if (ielm .eq. 0) then
                  if (i .gt. 1) then
                    exit L_IT
                  else
                    write(*,*) "Failed to find axis"
                    stop 1
                  end if
                end if
                ielm_tmp  = element_neighbours(2,ielm)
              
                if (ielm_prev .ne. ielm_tmp) write(*,*) ' WARNING : CHANGE OF ORIENTATION (2)'
                
                s_line = 1.d0
                call interp_RZP(node_list,element_list,ielm,s_line,t_line,p_line,R_out,z_out)

                if ((abs(R_in - R_out) .gt. 1.d-8) .or. (abs(z_in - z_out) .gt. 1.d-8)) &
                  write(*,'(A,2i6,4f8.4)') ' error in element change (2) ',ielm_prev,ielm,R_in,R_out,z_in,z_out
              end if
            else
              if (t_line + delta_t .gt. 1.d0) then
                s_line = s_line + small_delta*delta_s
                t_line = 1.d0
                p_line = p_line + small_delta*delta_phi_step
                call interp_RZP(node_list,element_list,ielm,s_line,t_line,p_line,R_in,z_in)

                ielm_prev = ielm      
                ielm      = element_neighbours(3,ielm_prev)
                if (ielm .eq. 0) then
                  if (i .gt. 1) then
                      exit L_IT
                    else
                      write(*,*) "Failed to find axis"
                      stop 1
                    end if
                end if
                ielm_tmp  = element_neighbours(1,ielm)

                if (ielm_prev .ne. ielm_tmp) write(*,*) ' WARNING : CHANGE OF ORIENTATION (3)'
              
                t_line = 0.d0
                call interp_RZP(node_list,element_list,ielm,s_line,t_line,p_line,R_out,z_out)

                if ((abs(R_in - R_out) .gt. 1.d-8) .or. (abs(z_in - z_out) .gt. 1.d-8)) &
                  write(*,'(A,2i6,4f8.4)') ' error in element change (3) ',ielm_prev,ielm,R_in,R_out,z_in,z_out
              else if (t_line + delta_t .lt. 0.d0) then
                s_line = s_line + small_delta*delta_s	
                t_line = 0.d0
                p_line = p_line + small_delta*delta_phi_step
                call interp_RZP(node_list,element_list,ielm,s_line,t_line,p_line,R_in,z_in)

                ielm_prev = ielm      
                ielm      = element_neighbours(1,ielm_prev)
                if (ielm .eq. 0) then
                  if (i .gt. 1) then
                      exit L_IT
                    else
                      write(*,*) "Failed to find axis"
                      stop 1
                    end if
                end if
                ielm_tmp  = element_neighbours(3,ielm)

                if (ielm_prev .ne. ielm_tmp) write(*,*) ' WARNING : CHANGE OF ORIENTATION (4)',ielm_prev,ielm
              
                t_line = 1.d0
                call interp_RZP(node_list,element_list,ielm,s_line,t_line,p_line,R_out,z_out)

                if ((abs(R_in - R_out) .gt. 1.d-8) .or. (abs(z_in - z_out) .gt. 1.d-8)) &
                  write(*,'(A,2i6,4f8.4)') ' error in element change (4) ',ielm_prev,ielm,R_in,R_out,z_in,z_out
              end if
            end if
          else
            s_line = s_line + delta_s
            t_line = t_line + delta_t
            p_line = p_line + delta_phi_step
            small_delta = 1.d0
          end if
         
          delta_phi_local = delta_phi_local + small_delta*delta_phi_step
         
          if (ielm .eq. 0) then
            if (i .gt. 1) then
              exit L_IT
            else
              write(*,*) "Failed to find axis"
              stop 1
            end if
          end if
        end do
       
        if (ielm .eq. 0) then
          if (i .gt. 1) then
            exit L_IT
          else
            write(*,*) "Failed to find axis"
            stop 1
          end if
        end if
      end do ! End of one toroidal turn
    
      call interp_RZP(node_list,element_list,ielm,s_line,t_line,p_line,R,z)
      R_line = R
      z_line = z
    
      ip = ip + 1
      Rp(ip) = R_line
      zp(ip) = z_line
    
      if (ielm .eq. 0) then
        if (i .gt. 1) then
          exit L_IT
        else
          write(*,*) "Failed to find axis"
          stop 1
        end if
      end if
    end do
    
    Rc = (maxval(Rp) + minval(Rp))/2.d0
    zc = (maxval(zp) + minval(zp))/2.d0
    error = sqrt((maxval(Rp) - minval(Rp))*(maxval(zp) - minval(zp)))/10.0
    write(*,'(A,I4,A,2E14.6,A,E14.6)') "Iteration ", i, "; center: ", Rc, zc, "; error: ", error
    
    ! If the error is less than the tolerance, the axis has been found.
    if (error .lt. tol) then
      success = .true.
      exit ! End the iterations, write results to a file and finish the program
    end if
    
    if (pos_out) then
      do j=1,n_turn
        write(13,'(2E18.8)') Rp(j), Zp(j)
      end do
      write(13,*)
    end if
    
    R_start = (1.d0 - mix)*R_start + mix*Rc
    z_start = (1.d0 - mix)*z_start + mix*zc
  end do L_IT ! End of one iteration
  
  ! If the iterations have exceeded n_iter and the axis has still not been found, print an error message and the best guess
  if (.not. success) then
    write(*,*) "Failed to find axis to desired tolerance"
    write(*,*) "Best guess: ", Rc, zc
    write(*,*) "Error: ", error
    
    if (pos_out) then
      do j=1,ip
        write(13,'(2E18.8)') Rp(j), Zp(j)
      end do
      write(13,*)
    end if
  end if
  if (pos_out) close(13)
  
  open(12,file=outfile,action="write",status="unknown",position="append")
  write(12,'(4E18.8)') t_start, Rc, zc, error
  close(12)
  
  contains
  subroutine step(i_elm,s_in,t_in,p_in,delta_p,delta_s,delta_t)
    use mod_parameters
    implicit none

    integer, intent(in)  :: i_elm
    real*8,  intent(in)  :: s_in, t_in, p_in, delta_p
    real*8,  intent(out) :: delta_s, delta_t

    integer :: i_tor, i_harm
    real*8  :: R,R_s,R_t,R_p,z,z_s,z_t,z_p,dummy,BR0cos,BR0sin,Bz0cos,Bz0sin,Bp0cos,Bp0sin
    real*8  :: Pcos,Pcos_s,Pcos_t,Pcos_st,Pcos_ss,Pcos_tt, Psin,Psin_s,Psin_t,Psin_st,Psin_ss,Psin_tt
    real*8  :: P0,P0_s,P0_t,P0_st,P0_ss,P0_tt, psi_s, psi_t, psi_R, psi_z, psi_p, st_psi_p, zjac
    real*8  :: AR0_z, AR0_p, AR0_s, AR0_t, Az0_R, Az0_p, Az0_s, Az0_t, A30_R, A30_z, BR0, Bz0, Bp0, Fprof

    real*8, dimension(0:n_order-1,0:n_order-1,0:n_order-1) :: chi

    call interp_RZP(node_list,element_list,i_elm,s_in,t_in,p_in,R,R_s,R_t,R_p,dummy,dummy,dummy,dummy,dummy,dummy, &
                                                            z,z_s,z_t,z_p,dummy,dummy,dummy,dummy,dummy,dummy)

    chi  = get_chi(R,z,p_in,node_list,element_list,i_elm,s_in,t_in)
    zjac = (R_s*z_t - R_t*z_s)

    call interp(node_list,element_list,i_elm,var_Psi,1,s_in,t_in,P0,P0_s,P0_t,P0_st,P0_ss,P0_tt)

    psi_s = P0_s 
    psi_t = P0_t 
    st_psi_p = 0.d0

#ifdef POINC_GVEC
    call interp_gvec(node_list,element_list,i_elm,1,1,1,s_in,t_in,BR0,dummy,dummy,dummy,dummy,dummy)
    call interp_gvec(node_list,element_list,i_elm,1,2,1,s_in,t_in,Bz0,dummy,dummy,dummy,dummy,dummy)
    call interp_gvec(node_list,element_list,i_elm,1,3,1,s_in,t_in,Bp0,dummy,dummy,dummy,dummy,dummy)
#endif

#ifdef fullmhd
    call interp(node_list,element_list,i_elm,var_AR,1,s_in,t_in,P0,P0_s,P0_t,P0_st,P0_ss,P0_tt)
    AR0_s = P0_s 
    AR0_t = P0_t 

    call interp(node_list,element_list,i_elm,var_Az,1,s_in,t_in,P0,P0_s,P0_t,P0_st,P0_ss,P0_tt)
    Az0_s = P0_s 
    Az0_t = P0_t 

    AR0_p = 0.d0
    Az0_p = 0.d0

    call interp(node_list,element_list,i_elm,710,1,s_in,t_in,P0,P0_s,P0_t,P0_st,P0_ss,P0_tt)
    Fprof = P0
#endif

    do i_tor=1,(n_tor-1)/2
      i_harm = 2*i_tor

      call interp(node_list,element_list,i_elm,var_Psi,i_harm,s_in,t_in,Pcos,Pcos_s,Pcos_t,Pcos_st,Pcos_ss,Pcos_tt)
      psi_s = psi_s + Pcos_s*cos(mode(i_harm)*p_in)
      psi_t = psi_t + Pcos_t*cos(mode(i_harm)*p_in)
      st_psi_p = st_psi_p - Pcos*mode(i_harm)*sin(mode(i_harm)*p_in)

      call interp(node_list,element_list,i_elm,var_Psi,i_harm+1,s_in,t_in,Psin,Psin_s,Psin_t,Psin_st,Psin_ss,Psin_tt)
      psi_s = psi_s + Psin_s*sin(mode(i_harm+1)*p_in)
      psi_t = psi_t + Psin_t*sin(mode(i_harm+1)*p_in)
      st_psi_p = st_psi_p + Psin*mode(i_harm+1)*cos(mode(i_harm+1)*p_in)

#ifdef fullmhd
      call interp(node_list,element_list,i_elm,var_AR,i_harm,s_in,t_in,Pcos,Pcos_s,Pcos_t,Pcos_st,Pcos_ss,Pcos_tt)
      AR0_s = AR0_s + Pcos_s*cos(mode(i_harm)*p_in)
      AR0_t = AR0_t + Pcos_t*cos(mode(i_harm)*p_in)
      AR0_p = AR0_p - Pcos  *sin(mode(i_harm)*p_in)*mode(i_harm)
      call interp(node_list,element_list,i_elm,var_AR,i_harm+1,s_in,t_in,Psin,Psin_s,Psin_t,Psin_st,Psin_ss,Psin_tt)
      AR0_s = AR0_s + Psin_s*sin(mode(i_harm+1)*p_in)
      AR0_t = AR0_t + Psin_t*sin(mode(i_harm+1)*p_in)
      AR0_p = AR0_p + Psin  *cos(mode(i_harm+1)*p_in)*mode(i_harm+1)

      call interp(node_list,element_list,i_elm,var_Az,i_harm,s_in,t_in,Pcos,Pcos_s,Pcos_t,Pcos_st,Pcos_ss,Pcos_tt)
      Az0_s = Az0_s + Pcos_s*cos(mode(i_harm)*p_in)
      Az0_t = Az0_t + Pcos_t*cos(mode(i_harm)*p_in)
      Az0_p = Az0_p - Pcos  *sin(mode(i_harm)*p_in)*mode(i_harm)
      call interp(node_list,element_list,i_elm,var_Az,i_harm+1,s_in,t_in,Psin,Psin_s,Psin_t,Psin_st,Psin_ss,Psin_tt)
      Az0_s = Az0_s + Psin_s*sin(mode(i_harm+1)*p_in)
      Az0_t = Az0_t + Psin_t*sin(mode(i_harm+1)*p_in)
      Az0_p = Az0_p + Psin  *cos(mode(i_harm+1)*p_in)*mode(i_harm+1)
#endif
    end do

#ifdef POINC_GVEC
    do i_tor=1,(n_coord_tor-1)/2
      i_harm = 2*i_tor

      call interp_gvec(node_list,element_list,i_elm,1,1,i_harm,s_in,t_in,BR0cos,dummy,dummy,dummy,dummy,dummy)
      call interp_gvec(node_list,element_list,i_elm,1,2,i_harm,s_in,t_in,Bz0cos,dummy,dummy,dummy,dummy,dummy)
      call interp_gvec(node_list,element_list,i_elm,1,3,i_harm,s_in,t_in,Bp0cos,dummy,dummy,dummy,dummy,dummy)
      BR0 = BR0 + BR0cos*cos(mode_coord(i_harm)*p_in)
      Bz0 = Bz0 + Bz0cos*cos(mode_coord(i_harm)*p_in)
      Bp0 = Bp0 + Bp0cos*cos(mode_coord(i_harm)*p_in)
  
      call interp_gvec(node_list,element_list,i_elm,1,1,i_harm+1,s_in,t_in,BR0sin,dummy,dummy,dummy,dummy,dummy)
      call interp_gvec(node_list,element_list,i_elm,1,2,i_harm+1,s_in,t_in,Bz0sin,dummy,dummy,dummy,dummy,dummy)
      call interp_gvec(node_list,element_list,i_elm,1,3,i_harm+1,s_in,t_in,Bp0sin,dummy,dummy,dummy,dummy,dummy)
      BR0 = BR0 - BR0sin*sin(mode_coord(i_harm+1)*p_in)
      Bz0 = Bz0 - Bz0sin*sin(mode_coord(i_harm+1)*p_in)
      Bp0 = Bp0 - Bp0sin*sin(mode_coord(i_harm+1)*p_in)
    end do
#endif

#ifdef fullmhd
    AR0_z = (-R_t*AR0_s + R_s*AR0_t)/zjac
    Az0_R = ( z_t*Az0_s - z_s*Az0_t)/zjac
    A30_R = ( z_t*psi_s - z_s*psi_t)/zjac
    A30_z = (-R_t*psi_s + R_s*psi_t)/zjac

    BR0 = (A30_z - Az0_p)/R
    Bz0 = (AR0_p - A30_R)/R
    Bp0 = (Az0_R - AR0_z) + Fprof/R
#else
    psi_R = ( z_t*psi_s - z_s*psi_t)/zjac
    psi_z = (-R_t*psi_s + R_s*psi_t)/zjac
    psi_p = st_psi_p - R_p*psi_R - z_p*psi_z

#ifndef POINC_GVEC
    BR0 = chi(1,0,0)   + (psi_z*chi(0,0,1) - psi_p*chi(0,1,0))/(F0*R) ! comment out these lines to use the
    Bz0 = chi(0,1,0)   - (psi_R*chi(0,0,1) - psi_p*chi(1,0,0))/(F0*R) !   GVEC magnetic field instead of
    Bp0 = chi(0,0,1)/R + (psi_R*chi(0,1,0) - psi_z*chi(1,0,0))/F0     !   the reduced MHD magnetic field
#endif
#endif

! dR/Rdphi = B_R / B_phi ; dz/Rdphi = B_z / B_phi
! ds/dphi = s_phi + s_R dR/dphi + s_z dz/dphi = (-z_t R_p + R_t z_p + z_t dR/dphi - R_t dz/dphi)/zjac
! dt/dphi = t_phi + t_R dR/dphi + t_z dz/dphi = ( z_s R_p - R_s z_p - z_s dR/dphi + R_s dz/dphi)/zjac
    delta_s = (-z_t*R_p + R_t*z_p + R*(z_t*BR0 - R_t*Bz0)/Bp0)*delta_p/zjac
    delta_t = ( z_s*R_p - R_s*z_p - R*(z_s*BR0 - R_s*Bz0)/Bp0)*delta_p/zjac
  end subroutine step
end program find_axis3D
