!> Create a Poincare plot for a JOREK restart file
program jorek2_poincare

use data_structure
use phys_module
use basis_at_gaussian
use elements_nodes_neighbours
use mod_neighbours
use mod_import_restart
use mod_log_params
use equil_info, only : get_psi_n, ES
use mod_interp

implicit none

character(len=512) :: s
real*8,allocatable  :: rp(:), zp(:), tp(:), pp(:)
integer, allocatable :: n_turn(:)
integer :: my_id, nr, ntour, curr
real*8  :: rr, zz, psi
integer :: i, j, iside_i, iside_j, ip, i_lines, n_lines, i_tor, i_harm, i_var_psi, iplot_type
integer :: i_elm, ifail, i_phi, n_phi, i_turn, i_elm_out, i_elm_prev, i_elm_tmp,i_steps
real*8  :: R_line, Z_line, s_line, t_line, p_line, R_mid, Z_mid, s_mid, t_mid, p_mid, s_out, t_out
real*8, allocatable :: R_start(:), Z_start(:), P_start(:)
real*8  :: R, Z, P, P_s, P_t, P_st, P_ss, P_tt
real*8  :: tol, delta_phi, Zjac, psi_s, psi_t, R_in, Z_in, R_out, Z_out, Rmin, Rmax, Zmin, Zmax, delta_s, delta_t, R_keep, Z_keep
real*8  :: small_delta, small_delta_s, small_delta_t, delta_phi_local, delta_phi_step
real*8  :: atmp, cur_pert
real*8  :: psi_out
integer :: ierr


write(*,*) '***************************************'
write(*,*) '* JOREK2_poincare                     *'
write(*,*) '***************************************'

my_id=0

! --- Initialize mode and mode_type arrays
call det_modes()

call initialise_parameters(my_id,  "__NO_FILENAME__")
call log_parameters(my_id)

iplot_type = 2 ! 1: Poincare plot in (R,Z) coordinates, 2: in (R,theta) coordinates

call import_restart(node_list,element_list, 'jorek_restart', rst_format, ierr, .true.)

call initialise_basis

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

!-------- possibilities
! - total length of traced field line (connection length)
! - keep all RZ positions of field lines
! - keep only crossing with fixed plane phi=constant
! - (local) q-profile

! step at constant delta_phi

! --- Read start points from file 'stpts'.
!
! Example for a stpts file:
!   Eleven field lines will be started, the first ten between (1.7,0.0) and (1.8,0.0), the eleventh at (1.85,0.2)
! 
!   +-------------------------------------------------------
!   |# n_lines
!   |  11
!   |# nr   R_start   Z_start    psi_start   n_turns
!   |   1    1.700      0.000     0.000      100
!   |  10    1.800      0.000     0.000      200
!   |  11    1.850      0.200     0.000      800
!   +-------------------------------------------------------
!
open(21, file='stpts', status='old', action='read', iostat=ierr)

if ( ierr == 0 ) then ! stpts file exists, use it.

  read(21, '(a)') s ! read comment line (ignored)
  read(21,*) n_lines
  if ( n_lines < 1 ) then
    write(*,*) 'ERROR in stpts file: n_lines must be >= 1.'
    stop
  end if
  read(21, '(A)') s ! read comment line (ignored)
  
  allocate( R_start(n_lines), Z_start(n_lines), P_start(n_lines), n_turn(n_lines) )
  
  curr = 0
  do
    if ( curr >= n_lines ) exit
    
    read(21, *) nr, rr, zz, psi, ntour
    
    if ( ( nr == 1 ) .and. ( curr == 0 ) ) then
      R_start(1) = rr
      Z_start(1) = zz
      P_start(1) = psi
      n_turn(1)  = ntour
    else if ( curr == 0 ) then
      write(*,*) 'ERROR in stpts file: first start point must be nr=1.'
      stop
    else if ( ( nr < 1 ) .or. ( nr > n_lines ) ) then
      write(*,*) 'ERROR in stpts file: nr must be > 0 and < n_lines.'
      stop
    else if ( nr <= curr ) then
      write(*,*) 'ERROR in stpts file: start points must be sorted in ascending nr-order.'
      stop
    else
      
      do i_lines = curr + 1, nr
        R_start(i_lines) = R_start(curr) + ( rr - R_start(curr) ) * ( real(i_lines-curr) / real(nr-curr) )
        Z_start(i_lines) = Z_start(curr) + ( zz - Z_start(curr) ) * ( real(i_lines-curr) / real(nr-curr) )
        P_start(i_lines) = P_start(curr) + ( psi- P_start(curr) ) * ( real(i_lines-curr) / real(nr-curr) )
        n_turn(i_lines)  = nint( n_turn(curr) + real( ntour - n_turn(curr) ) * ( real(i_lines-curr) / real(nr-curr) ) )
      end do
      
    end if
    
    curr = nr
    
  end do
  
  close(21)

else ! if no stpts file exists, use the following hard-coded default startpoints

  n_lines = 50

  allocate( R_start(n_lines), Z_start(n_lines), P_start(n_lines), n_turn(n_lines) )

  do i_lines = 1, n_lines
    R_start(i_lines) = 1.7156 + (2.18-1.7156) * float(i_lines-1)/float(n_lines-1)
    Z_start(i_lines) = 0.12237
    P_start(i_lines) = 0.d0
  end do

  n_turn  = 500
 
end if

n_phi   = 1500
delta_phi = 2.d0 * PI / float(n_period*n_phi)
tol       = 1.d-6

i_var_psi = 1

allocate(Rp(maxval(n_turn)),Zp(maxval(n_turn)), Tp(maxval(n_turn)),Pp(maxval(n_turn)))

Rmin = 1.d20; Rmax = -1.d20; Zmin = 1.d20; Zmax=-1.d20
do i=1,node_list%n_nodes
  Rmin = min(Rmin,node_list%node(i)%x(1,1,1))
  Rmax = max(Rmax,node_list%node(i)%x(1,1,1))
  Zmin = min(Zmin,node_list%node(i)%x(1,1,2))
  Zmax = max(Zmax,node_list%node(i)%x(1,1,2))
enddo

mode(1) = 0
do i=1,(n_tor-1)/2
  mode(2*i)   = i * n_period
  mode(2*i+1) = i * n_period
enddo
write(*,*) ' modes   : ',mode
write(*,*) ' nperiod : ',n_period

  
call begplt('poincare.ps')

! --- Open the output files to which the Poincare data will be written in ascii format
open(21,file='poinc_R-Z.dat')
write(21,*) '#  R                 Z'
open(22,file='poinc_rho-theta.dat')
write(22,*) '# rho=sqrt(psi_n)'
write(22,*) '# psi_n=(psi - psi_axis)/(psi_bnd - psi_axis)'
write(22,*) '#'
write(22,*) '#  rho               theta'

if (iplot_type .eq. 1) then
  call nframe(21,11,1,Rmin,Rmax,Zmin,Zmax,'Poincare',8,'R [m]',4,'Z [m]',4)
else
  call nframe(1,11,1,0.d0,1.2d0,-PI,PI,'Poincare',8,'r [m]',4,'theta',5)
endif

! --- Trace the fieldlines
L_IL: do i_lines=1,n_lines
  ip = 0

  write(*,*)
  write(*,'(1x,2(a,i6),a,2f8.3)') 'Line',i_lines,' of',n_lines,' started at',R_start(i_lines),Z_start(i_lines)

  call find_RZ(node_list,element_list,R_start(i_lines),Z_start(i_lines),R_out,Z_out,i_elm,s_out,t_out,ifail)
 
  if (ifail .ne. 0) write(*,*) "Can not find RZ,", ifail 
  if (ifail .ne. 0) exit

  R_line = R_start(i_lines)
  Z_line = Z_start(i_lines)
  p_line = P_start(i_lines)
  s_line = s_out
  t_line = t_out
  
  L_IT: do i_turn = 1, n_turn(i_lines)
    if ( mod(i_turn-1,max(n_turn(i_lines)/6+1,5)) == 0 ) then
      write(*,'(3x,2(a,i6))') 'Turn',i_turn,' of',n_turn(i_lines)
    end if

    do i_phi=1,n_phi
    
      delta_phi_local = 0.d0
      
      i_steps = 0
    
      do while ((delta_phi_local .lt. delta_phi) .and. (i_steps .lt.10) )
      
        i_steps = i_steps + 1
      
        delta_phi_step = delta_phi - delta_phi_local
	
!	write(*,'(5i6,3e16.8)') i_lines,i_turn,i_phi,i_steps,i_elm,delta_phi_step, delta_phi, delta_phi_local
      
        call step(i_elm,s_line,t_line,p_line,delta_phi_step,delta_s,delta_t)

        s_mid = s_line + 0.5d0 * delta_s
        t_mid = t_line + 0.5d0 * delta_t
        p_mid = p_line + 0.5d0 * delta_phi_step
            
        call step(i_elm,s_mid,t_mid,p_mid,delta_phi_step,delta_s,delta_t)
        
        small_delta_s = 1.d0
      
        if  (s_line + delta_s .gt. 1.d0) then
      
          small_delta_s = (1.d0 - s_line)/delta_s
 
        elseif  (s_line + delta_s .lt. 0.d0) then

          small_delta_s = abs(s_line/delta_s)
	
        endif
      
        small_delta_t = 1.d0

        if  (t_line + delta_t .gt. 1.d0)  then
      
          small_delta_t = (1.d0 - t_line)/delta_t
	
        elseif  (t_line + delta_t .lt. 0.d0)  then

          small_delta_t = abs(t_line/delta_t)
      
        endif      

        small_delta = min(small_delta_s, small_delta_t)
	
!	write(*,'(A,5e16.8)') ' small delta : ',small_delta,delta_s,delta_t,s_line,t_line

        if (small_delta .lt. 1.d0)  then 

          s_mid = s_line + 0.5d0 * small_delta * delta_s
          t_mid = t_line + 0.5d0 * small_delta * delta_t
          p_mid = p_line + 0.5d0 * small_delta * delta_phi_step
            
          call step(i_elm,s_mid,t_mid,p_mid,delta_phi_step,delta_s,delta_t)
      
          if (small_delta_s .lt. small_delta_t) then

            if (s_line + delta_s .gt. 1.d0) then
	    
	      s_line = 1.d0
              t_line = t_line + small_delta * delta_t
              p_line = p_line + small_delta * delta_phi_step
	      
              call interp_RZ(node_list,element_list,i_elm,s_line,t_line,R_in,Z_in)

	      i_elm_prev = i_elm      
              i_elm      = element_neighbours(2,i_elm_prev)
              if ( i_elm == 0 ) exit L_IT
	      i_elm_tmp  = element_neighbours(4,i_elm)
	      
	      if (i_elm_prev .ne. i_elm_tmp) write(*,*) ' WARNING : CHANGE OF ORIENTATION (1)'
	
              s_line = 0.d0

              call interp_RZ(node_list,element_list,i_elm,s_line,t_line,R_out,Z_out)
	      
	      if ( (abs(R_in - R_out) .gt. 1.d-8) .or. (abs(Z_in - Z_out) .gt. 1.d-8)) &
	        write(*,'(A,2i6,4f8.4)') ' error in element change (1) ',i_elm_prev,i_elm,R_in,R_out,Z_in,Z_out
	
            elseif (s_line + delta_s .lt. 0.d0) then
	
              s_line = 0.d0
	      t_line = t_line + small_delta * delta_t
	      p_line = p_line + small_delta * delta_phi_step

              call interp_RZ(node_list,element_list,i_elm,s_line,t_line,R_in,Z_in)
      
	      i_elm_prev = i_elm      
              i_elm      = element_neighbours(4,i_elm_prev)
              if ( i_elm == 0 ) exit L_IT
	      i_elm_tmp  = element_neighbours(2,i_elm)
	      
	      if (i_elm_prev .ne. i_elm_tmp) write(*,*) ' WARNING : CHANGE OF ORIENTATION (2)'

              s_line = 1.d0

              call interp_RZ(node_list,element_list,i_elm,s_line,t_line,R_out,Z_out)
	      
	      if ( (abs(R_in - R_out) .gt. 1.d-8) .or. (abs(Z_in - Z_out) .gt. 1.d-8)) &
	        write(*,'(A,2i6,4f8.4)') ' error in element change (2) ',i_elm_prev,i_elm,R_in,R_out,Z_in,Z_out
	  
	    endif
	
	  else
	
            if (t_line + delta_t .gt. 1.d0) then
      
              s_line = s_line + small_delta * delta_s
              t_line = 1.d0
              p_line = p_line + small_delta * delta_phi_step

              call interp_RZ(node_list,element_list,i_elm,s_line,t_line,R_in,Z_in)

	      i_elm_prev = i_elm      
              i_elm      = element_neighbours(3,i_elm_prev)
              if ( i_elm == 0 ) exit L_IT
	      i_elm_tmp  = element_neighbours(1,i_elm)
	      
	      if (i_elm_prev .ne. i_elm_tmp) write(*,*) ' WARNING : CHANGE OF ORIENTATION (3)'
	
              t_line = 0.d0

              call interp_RZ(node_list,element_list,i_elm,s_line,t_line,R_out,Z_out)
	      
	      if ( (abs(R_in - R_out) .gt. 1.d-8) .or. (abs(Z_in - Z_out) .gt. 1.d-8)) &
	        write(*,'(A,2i6,4f8.4)') ' error in element change (3) ',i_elm_prev,i_elm,R_in,R_out,Z_in,Z_out
	
            elseif (t_line + delta_t .lt. 0.d0) then

              s_line = s_line + small_delta * delta_s	
	      t_line = 0.d0
	      p_line = p_line + small_delta * delta_phi_step
 
              call interp_RZ(node_list,element_list,i_elm,s_line,t_line,R_in,Z_in)

	      i_elm_prev = i_elm      
              i_elm      = element_neighbours(1,i_elm_prev)
              if ( i_elm == 0 ) exit L_IT
	      i_elm_tmp  = element_neighbours(3,i_elm)
	      
	      if (i_elm_prev .ne. i_elm_tmp) write(*,*) ' WARNING : CHANGE OF ORIENTATION (4)',i_elm_prev,i_elm
 
              t_line = 1.d0

              call interp_RZ(node_list,element_list,i_elm,s_line,t_line,R_out,Z_out)
	      
	      if ( (abs(R_in - R_out) .gt. 1.d-8) .or. (abs(Z_in - Z_out) .gt. 1.d-8)) &
	        write(*,'(A,2i6,4f8.4)') ' error in element change (4) ',i_elm_prev,i_elm,R_in,R_out,Z_in,Z_out
	
            endif
        
	  endif
         
        else      

          s_line = s_line + delta_s
          t_line = t_line + delta_t
          p_line = p_line + delta_phi_step
	
	  small_delta = 1.d0
      
        endif
      
        delta_phi_local = delta_phi_local + small_delta * delta_phi_step
     
        if (i_elm .eq. 0) exit
      
!        write(*,'(A,5e16.8)') ' s,t : ',s_line,t_line
      
      enddo

      if (i_elm .eq. 0) exit
      
!      if (i_steps .gt. 8) write(*,'(A,5i6)') ' WARNING : isteps ',i_lines,i_turn,i_phi,i_steps,i_elm
            
    enddo ! end of a 2Pi turn

    call interp_RZ(node_list,element_list,i_elm,s_line,t_line,R,Z)

    R_line = R
    Z_line = Z
            
    ip = ip+1

    call var_value(i_elm,1,s_line,t_line,p_line,psi_out)

    Rp(ip) = R_line
    Zp(ip) = Z_line
    Tp(ip)  = atan2( Z_line - ES%Z_axis, R_line - ES%R_axis)
    Pp(ip)  = get_psi_n(psi_out, Z_line)

    if (i_elm .eq. 0) exit
     
  enddo L_IT
  
  write(*,'(3x,a,i6,a)') '=>',ip,' points'

  do i=1,ip
    write(21,'(4e18.8)') Rp(i),Zp(i)
    write(22,'(4e18.8)') SQRT( MAX(Pp(i), 0.) ),Tp(i)
  enddo

  write(21,*)
  write(21,*)
  write(22,*)
  write(22,*)
  
  call lincol(mod(i_lines,8))
   
  if (iplot_type .eq. 1) then
    call pplot(1,1,Rp,Zp,ip,1)
  else
    call pplot(1,1,Pp,Tp,ip,1)
  endif
  
end do L_IL

close(21)
close(22)
call finplt

end program jorek2_poincare



subroutine step(i_elm,s_in,t_in,p_in,delta_p,delta_s,delta_t)
use mod_parameters
use elements_nodes_neighbours
use phys_module
use mod_interp

implicit none

integer :: i_var_psi, i_elm, i_tor, i_harm

real*8 :: s_in, t_in, p_in, delta_p, delta_s, delta_t
real*8 :: R,R_s,R_t,Z,Z_s,Z_t
real*8 :: Pcos,Pcos_s,Pcos_t,Pcos_st,Pcos_ss,Pcos_tt, Psin,Psin_s,Psin_t,Psin_st,Psin_ss,Psin_tt
real*8 :: P0,P0_s,P0_t,P0_st,P0_ss,P0_tt, psi_s, psi_t, Zjac
real*8 :: AR0_Z, AR0_p, AR0_s, AR0_t, AZ0_R, AZ0_p, AZ0_s, AZ0_t, A30_R, A30_Z, BR0, BZ0, Bp0, Fprof

i_var_psi = 1

call interp_RZ(node_list,element_list,i_elm,s_in,t_in,R,R_s,R_t,Z,Z_s,Z_t)

Zjac = (R_s * Z_t - R_t * Z_s)

call interp(node_list,element_list,i_elm,i_var_psi,1,s_in,t_in,P0,P0_s,P0_t,P0_st,P0_ss,P0_tt)

psi_s = P0_s 
psi_t = P0_t 

#ifdef fullmhd
  call interp(node_list,element_list,i_elm,var_AR,1,s_in,t_in,P0,P0_s,P0_t,P0_st,P0_ss,P0_tt)
  AR0_s = P0_s 
  AR0_t = P0_t 

  call interp(node_list,element_list,i_elm,var_AZ,1,s_in,t_in,P0,P0_s,P0_t,P0_st,P0_ss,P0_tt)
  AZ0_s = P0_s 
  AZ0_t = P0_t 

  AR0_p = 0.d0
  AZ0_p = 0.d0

  call interp(node_list,element_list,i_elm,710,1,s_in,t_in,P0,P0_s,P0_t,P0_st,P0_ss,P0_tt)
  Fprof = P0
#endif

do i_tor = 1, (n_tor-1)/2

  i_harm = 2*i_tor

  call interp(node_list,element_list,i_elm,i_var_psi,i_harm,s_in,t_in,Pcos,Pcos_s,Pcos_t,Pcos_st,Pcos_ss,Pcos_tt)

  psi_s = psi_s + Pcos_s * cos(mode(i_harm)*p_in)
  psi_t = psi_t + Pcos_t * cos(mode(i_harm)*p_in)

  call interp(node_list,element_list,i_elm,i_var_psi,i_harm+1,s_in,t_in,Psin,Psin_s,Psin_t,Psin_st,Psin_ss,Psin_tt)

  psi_s = psi_s + Psin_s * sin(mode(i_harm+1)*p_in)
  psi_t = psi_t + Psin_t * sin(mode(i_harm+1)*p_in)

#ifdef fullmhd
  call interp(node_list,element_list,i_elm,var_AR,i_harm,s_in,t_in,Pcos,Pcos_s,Pcos_t,Pcos_st,Pcos_ss,Pcos_tt)
  AR0_s = AR0_s + Pcos_s * cos(mode(i_harm)*p_in)
  AR0_t = AR0_t + Pcos_t * cos(mode(i_harm)*p_in)
  AR0_p = AR0_p - Pcos   * sin(mode(i_harm)*p_in) * mode(i_harm)
  call interp(node_list,element_list,i_elm,var_AR,i_harm+1,s_in,t_in,Psin,Psin_s,Psin_t,Psin_st,Psin_ss,Psin_tt)
  AR0_s = AR0_s + Psin_s * sin(mode(i_harm+1)*p_in)
  AR0_t = AR0_t + Psin_t * sin(mode(i_harm+1)*p_in)
  AR0_p = AR0_p + Psin   * cos(mode(i_harm+1)*p_in) * mode(i_harm+1)

  call interp(node_list,element_list,i_elm,var_AZ,i_harm,s_in,t_in,Pcos,Pcos_s,Pcos_t,Pcos_st,Pcos_ss,Pcos_tt)
  AZ0_s = AZ0_s + Pcos_s * cos(mode(i_harm)*p_in)
  AZ0_t = AZ0_t + Pcos_t * cos(mode(i_harm)*p_in)
  AZ0_p = AZ0_p - Pcos   * sin(mode(i_harm)*p_in) * mode(i_harm)
  call interp(node_list,element_list,i_elm,var_AZ,i_harm+1,s_in,t_in,Psin,Psin_s,Psin_t,Psin_st,Psin_ss,Psin_tt)
  AZ0_s = AZ0_s + Psin_s * sin(mode(i_harm+1)*p_in)
  AZ0_t = AZ0_t + Psin_t * sin(mode(i_harm+1)*p_in)
  AZ0_p = AZ0_p + Psin   * cos(mode(i_harm+1)*p_in) * mode(i_harm+1)
#endif

enddo

#ifdef fullmhd
AR0_Z = ( - R_t * AR0_s  + R_s * AR0_t ) / Zjac
AZ0_R = (   Z_t * AZ0_s  - Z_s * AZ0_t ) / Zjac
A30_R = (   Z_t * psi_s  - Z_s * psi_t ) / Zjac
A30_Z = ( - R_t * psi_s  + R_s * psi_t ) / Zjac

BR0 = ( A30_Z - AZ0_p )/ R
BZ0 = ( AR0_p - A30_R )/ R
Bp0 = ( AZ0_R - AR0_Z )       +   Fprof / R

! dR/Rdphi = B_R / B_phi ; dZ/Rdphi = B_Z / B_phi
! ds = (Z_t dR - R_t dZ) / Zjac ; dt = ( -Z_s dR + R_s dZ) / Zjac
delta_s =  ( Z_t*BR0 - R_t*BZ0) / ( Bp0 * Zjac ) * R * delta_p
delta_t =  (-Z_s*BR0 + R_s*BZ0) / ( Bp0 * Zjac ) * R * delta_p     

#else
delta_s =   psi_t * R / (Zjac * F0) * delta_p
delta_t = - psi_s * R / (Zjac * F0) * delta_p
#endif

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
